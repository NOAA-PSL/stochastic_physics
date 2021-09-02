module cellular_automata_global_mod

use update_ca, only : rstate_global
implicit none

contains

subroutine cellular_automata_global(kstep,restart,first_time_step,ca1_cpl,ca2_cpl,ca3_cpl,ca1_diag,ca2_diag,ca3_diag, &
            domain_for_coupler,nblks,isc,iec,jsc,jec,npx,npy,nlev, &
            nca,ncells,nlives,nfracseed,nseed,iseed_ca, mytile, &
            ca_smooth,nspinup,blocksize,nsmooth,ca_amplitude,mpiroot,mpicomm,noise_option)

use kinddef,           only: kind_phys
use update_ca,         only: update_cells_global,define_ca_domain
use halo_exchange,     only: atmosphere_scalar_field_halo
use mersenne_twister,  only: random_setseed,random_gauss,random_stat,random_number
use mpp_domains_mod,   only: domain2D,mpp_get_global_domain,CENTER, mpp_get_data_domain, mpp_get_compute_domain
use block_control_mod, only: block_control_type, define_blocks_packed
use mpi_wrapper,       only: mp_reduce_sum,mp_bcst,mp_reduce_max,mp_reduce_min, &
                             mpi_wrapper_initialize,mype
use mpp_mod

implicit none

!L.Bengtsson, 2017-06

!This program evolves a cellular automaton uniform over the globe 

integer,              intent(in)    :: kstep,ncells,nca,nlives,nseed,nspinup,nsmooth,mpiroot,mpicomm
integer*8,            intent(in)    :: iseed_ca
integer,              intent(in)    :: mytile
real(kind=kind_phys), intent(in)    :: nfracseed,ca_amplitude
logical,              intent(in)    :: ca_smooth,first_time_step, restart
integer,              intent(in)    :: nblks,isc,iec,jsc,jec,npx,npy,nlev,blocksize
real(kind=kind_phys), intent(out)   :: ca1_cpl(:,:),ca2_cpl(:,:),ca3_cpl(:,:)
real(kind=kind_phys), intent(out)   :: ca1_diag(:,:),ca2_diag(:,:),ca3_diag(:,:)
type(domain2D),       intent(inout) :: domain_for_coupler
type(domain2D),save                 :: domain_ncellx    
integer,              intent(in)    :: noise_option
type(block_control_type) :: Atm_block
integer :: nlon, nlat, isize,jsize,nf,nn
integer :: inci, incj, nxc, nyc, nxch, nych
integer :: halo, k_in, i, j, k
integer :: seed, ierr7,blk, ix, iix, count4,ih,jh
integer :: blocksz,levs
integer,save :: isdnx,iednx,jsdnx,jednx
integer,save :: iscnx,iecnx,jscnx,jecnx
integer :: nxncells, nyncells
integer(8) :: count, count_rate, count_max, count_trunc,nx_full,ny_full
integer(8) :: iscale = 10000000000
integer, allocatable :: iini_g(:,:,:),ilives_g(:,:)
real(kind=kind_phys), allocatable :: field_out(:,:,:), field_smooth(:,:)
real(kind=kind_phys), allocatable :: CA(:,:),CA1(:,:),CA2(:,:),CA3(:,:)
real*8              , allocatable :: noise1D(:),noise2D(:,:),noise(:,:,:)
real*8               :: psum,CAmean,sq_diff,CAstdv
real*8               :: Detmax,Detmin
integer*8            :: csum
logical,save         :: block_message=.true.
logical              :: global_rescale=.true.
integer*8            :: i1,j1


!nca         :: switch for number of cellular automata to be used.
!nfracseed   :: switch for number of random cells initially seeded
!nlives      :: switch for maximum number of lives a cell can have
!nspinup     :: switch for number of itterations to spin up the ca
!ncells      :: switch for higher resolution grid e.g ncells=4
!               gives 4x4 times the FV3 model grid resolution.
!ca_smooth   :: switch to smooth the cellular automata
if (nca .LT. 1) return
! Initialize MPI and OpenMP
if (first_time_step) then
   call mpi_wrapper_initialize(mpiroot,mpicomm)
end if

halo=1
k_in=1

!----------------------------------------------------------------------------
! Get information about the compute domain, allocate fields on this
! domain

! WRITE(*,*)'Entering cellular automata calculations'

! Some security checks for namelist combinations:
 if(nca > 3)then
 write(0,*)'Namelist option nca cannot be larger than 3 - exiting'
 stop
 endif


 nlon=iec-isc+1
 nlat=jec-jsc+1
 isize=nlon+2*halo
 jsize=nlat+2*halo

 inci=ncells
 incj=ncells

 
 !--- get params from domain_ncellx for building board and board_halo
 !Get CA domain
                                                                                                                                                     
 if(first_time_step)then 
    call define_ca_domain(domain_for_coupler,domain_ncellx,ncells,nxncells,nyncells)
    call mpp_get_data_domain    (domain_ncellx,isdnx,iednx,jsdnx,jednx)
    call mpp_get_compute_domain (domain_ncellx,iscnx,iecnx,jscnx,jecnx)
 endif

 nxc = iecnx-iscnx+1
 nyc = jecnx-jscnx+1
 nxch = iednx-isdnx+1
 nych = jednx-jsdnx+1


 !Allocate fields:
 allocate(field_out(isize,jsize,1))
 allocate(field_smooth(nlon,nlat))
 allocate(iini_g(nxc,nyc,nca))
 allocate(ilives_g(nxc,nyc))
 allocate(CA(nlon,nlat))
 allocate(CA1(nlon,nlat))
 allocate(CA2(nlon,nlat))
 allocate(CA3(nlon,nlat))
 allocate(noise(nxc,nyc,nca))
 nx_full=int(npx-1,kind=8)
 ny_full=int(npy-1,kind=8)
if (noise_option.eq.0) then
    allocate(noise1D(nxc*nyc))
else
    allocate(noise1D(nx_full*ny_full))
    allocate(noise2D(nx_full,ny_full))
endif

 !Initialize:

 noise(:,:,:) = 0.0
 noise1D(:) = 0.0
if (noise_option.eq.1) then
 noise2D(:,:) = 0.0
endif
 iini_g(:,:,:) = 0
 ilives_g(:,:) = 0
 CA1(:,:) = 0.0
 CA2(:,:) = 0.0
 CA3(:,:) = 0.0

 !Put the blocks of model fields into a 2d array - can't use nlev and blocksize directly,
 !because the arguments to define_blocks_packed are intent(inout) and not intent(in).
 levs=nlev
 blocksz=blocksize

 call define_blocks_packed('cellular_automata', Atm_block, isc, iec, jsc, jec, levs, &
                              blocksz, block_message)

if(first_time_step.and. .not.restart)then
!Generate random number, following stochastic physics code:
      print*,'initialize sgs noise_option',noise_option
  if (noise_option .EQ. 0) then

     if (iseed_ca == 0) then
       ! generate a random seed from system clock and ens member number
       call system_clock(count, count_rate, count_max)
       ! iseed is elapsed time since unix epoch began (secs)
       ! truncate to 4 byte integer
       count_trunc = iscale*(count/iscale)
       count4 = count - count_trunc
     else if (iseed_ca > 0) then
       ! don't rely on compiler to truncate integer(8) to integer(4) on
       ! overflow, do wrap around explicitly.
       count4 = mod(mype + iseed_ca +1+ 2147483648, 4294967296) - 2147483648
     endif
   
     call random_setseed(count4,rstate_global)
     do nf=1,nca
     !Set seed (to be different) on all tasks. Save random state.
       call random_number(noise1D,rstate_global)
     !Put on 2D:
       do j=1,nyc
         do i=1,nxc
           noise(i,j,nf)=noise1D(i+(j-1)*nxc)
         enddo
       enddo
      enddo
  else
     if (iseed_ca == 0) then
       ! generate a random seed from system clock and ens member number
       call system_clock(count, count_rate, count_max)
       ! iseed is elapsed time since unix epoch began (secs)
       ! truncate to 4 byte integer
       count_trunc = iscale*(count/iscale)
       count4 = count - count_trunc
     else if (iseed_ca > 0) then
       ! don't rely on compiler to truncate integer(8) to integer(4) on
       ! overflow, do wrap around explicitly.
       count4 = mod(mytile + iseed_ca+1 + 2147483648, 4294967296) - 2147483648
     endif
   
     call random_setseed(count4,rstate_global)
     do nf=1,nca
     !Set seed (to be different) on all tasks. Save random state.
       call random_number(noise1D,rstate_global)
     if (noise_option .EQ. 1) then
       noise2D(:,:)=RESHAPE(noise1D,(/nx_full,ny_full/))
       !pick out points on this task
       do j=1,nyc
         do i=1,nxc
           noise(i,j,nf)=noise2D( ncells*isc-ncells+i,ncells*jsc-ncells+j)
         enddo
       enddo
      enddo
      else
         !pick out points on my task
         do j=1,nyc
            j1=j-1+isc*ncells
            do i=1,nxc
                i1=i-1+isc*ncells
                noise(i,j)=NOISE1D( i1+nx_full*j1)
            enddo
         enddo
         deallocate(noise1D)
       endif
  endif

!Initiate the cellular automaton with random numbers larger than nfracseed

   do nf=1,nca
    do j = 1,nyc
      do i = 1,nxc
        if (noise(i,j,nf) > nfracseed ) then
          iini_g(i,j,nf)=1
        else
          iini_g(i,j,nf)=0
        endif
      enddo
    enddo
  enddo !nf
endif ! 

!In case we want to condition the cellular automaton on a large scale field
!we here set the "condition" variable to a different model field depending
!on nf. (this is not used if ca_global = .true.)


  do nf=1,nca !update each ca

   do j = 1,nyc
    do i = 1,nxc
     ilives_g(i,j)=int(real(nlives)*1.5*noise(i,j,nf))
    enddo
   enddo



!Calculate neighbours and update the automata
!If ca-global is used, then nca independent CAs are called and weighted together to create one field; CA


  CA(:,:)=0.

  call update_cells_global(kstep,first_time_step,restart,nca,nxc,nyc,nxch,nych,nlon,nlat,nxncells,nyncells,isc,iec,jsc,jec, &
                           npx,npy,iscnx,iecnx,jscnx,jecnx,domain_ncellx,CA,iini_g,ilives_g,         &
                           nlives,ncells,nfracseed,nseed,nspinup,nf,noise_option)


if (ca_smooth) then

   do nn=1,nsmooth !number of iterations for the smoothing.

      field_out=0.
      field_out(1+halo:nlon+halo,1+halo:nlat+halo,1) = real(CA(1:nlon,1:nlat),kind=8)

      call atmosphere_scalar_field_halo(field_out,halo,isize,jsize,k_in,isc,iec,jsc,jec,npx,npy,domain_for_coupler)

      do j=1,nlat
         do i=1,nlon
            ih=i+halo
            jh=j+halo
            field_smooth(i,j)=(2.0*field_out(ih,jh,1)+2.0*field_out(ih-1,jh,1)+ &
                               2.0*field_out(ih,jh-1,1)+2.0*field_out(ih+1,jh,1)+&
                               2.0*field_out(ih,jh+1,1)+2.0*field_out(ih-1,jh-1,1)+&
                               2.0*field_out(ih-1,jh+1,1)+2.0*field_out(ih+1,jh+1,1)+&
                               2.0*field_out(ih+1,jh-1,1))/18.
            enddo
         enddo
      do j=1,nlat
         do i=1,nlon
            CA(i,j)=field_smooth(i,j)
         enddo
      enddo

   enddo !nn
endif !smooth

!!!!Post processing, should be made into a separate subroutine
if (global_rescale) then
   
   !mean:
   psum=SUM(CA)
   csum=int(nlon,kind=8)*int(nlat,kind=8)
   
   call mp_reduce_sum(psum)
   call mp_reduce_sum(csum)
   
   CAmean=psum/csum
   
   !std:
   CAstdv=0.
   sq_diff = 0.
   do j=1,nlat
    do i=1,nlon
     sq_diff = sq_diff + (CA(i,j)-CAmean)**2.0
    enddo
   enddo
   
   call mp_reduce_sum(sq_diff)
   
   CAstdv = sqrt(sq_diff/csum)
   
   !Transform to mean of 1 and ca_amplitude standard deviation
      
   do j=1,nlat
    do i=1,nlon
     CA(i,j)=1.0 + (CA(i,j)-CAmean)*(ca_amplitude/CAstdv)
    enddo
   enddo
   
   do j=1,nlat
    do i=1,nlon
       CA(i,j)=min(max(CA(i,j),0.),2.0)
    enddo
   enddo
endif
!Put back into blocks 1D array to be passed to physics
!or diagnostics output

    if(nf==1)then
       CA1(:,:)=CA(:,:)
    elseif(nf==2)then
       CA2(:,:)=CA(:,:)
    else
       CA3(:,:)=CA(:,:)
    endif

enddo !nf

  do blk = 1, Atm_block%nblks
  do ix = 1,Atm_block%blksz(blk)
     i = Atm_block%index(blk)%ii(ix) - isc + 1
     j = Atm_block%index(blk)%jj(ix) - jsc + 1
     ca1_diag(blk,ix)=CA1(i,j)
     ca2_diag(blk,ix)=CA2(i,j)
     ca3_diag(blk,ix)=CA3(i,j)
     ca1_cpl(blk,ix)=CA1(i,j)
     ca2_cpl(blk,ix)=CA2(i,j)
     ca3_cpl(blk,ix)=CA3(i,j)
  enddo
  enddo


 deallocate(field_out)
 deallocate(field_smooth)
 deallocate(iini_g)
 deallocate(ilives_g)
 deallocate(CA)
 deallocate(CA1)
 deallocate(CA2)
 deallocate(CA3)
 deallocate(noise)
 deallocate(noise1D)
if (noise_option.eq.1) then
 deallocate(noise2D)
endif

end subroutine cellular_automata_global

end module cellular_automata_global_mod
