subroutine cellular_automata_global(kstep,Statein,Coupling,Diag,nblks,nlev, &
            nca,ncells,nlives,nfracseed,nseed,nthresh,ca_global,ca_sgs,iseed_ca, &
            ca_smooth,nspinup,blocksize,nsmooth,ca_amplitude)

use machine
use update_ca,         only: update_cells_sgs, update_cells_global
#ifdef STOCHY_UNIT_TEST
use standalone_stochy_module,      only: GFS_Coupling_type, GFS_diag_type, GFS_statein_type
use atmosphere_stub_mod,    only: atmosphere_resolution, atmosphere_domain, &
                                   atmosphere_scalar_field_halo, atmosphere_control_data
#else
use GFS_typedefs,      only: GFS_Coupling_type, GFS_diag_type, GFS_statein_type
use atmosphere_mod,    only: atmosphere_resolution, atmosphere_domain, &
                             atmosphere_scalar_field_halo, atmosphere_control_data
#endif
use mersenne_twister,  only: random_setseed,random_gauss,random_stat,random_number
use mpp_domains_mod,   only: domain2D
use block_control_mod, only: block_control_type, define_blocks_packed
use fv_mp_mod,         only : mp_reduce_sum,mp_bcst,mp_reduce_max,mp_reduce_min,is_master


implicit none

!L.Bengtsson, 2017-06

!This program evolves a cellular automaton uniform over the globe given
!the flag ca_global

integer,intent(in) :: kstep,ncells,nca,nlives,nseed,iseed_ca,nspinup,nsmooth
real,intent(in) :: nfracseed,nthresh,ca_amplitude
logical,intent(in) :: ca_global, ca_sgs, ca_smooth
integer, intent(in) :: nblks,nlev,blocksize
type(GFS_coupling_type),intent(inout) :: Coupling(nblks)
type(GFS_diag_type),intent(inout) :: Diag(nblks)
type(GFS_statein_type),intent(in) :: Statein(nblks)
type(block_control_type)          :: Atm_block
type(random_stat) :: rstate
integer :: nlon, nlat, isize,jsize,nf,nn
integer :: inci, incj, nxc, nyc, nxch, nych
integer :: halo, k_in, i, j, k, iec, jec, isc, jsc
integer :: seed, ierr7,blk, ix, iix, count4,ih,jh
integer :: blocksz,levs
integer(8) :: count, count_rate, count_max, count_trunc
integer(8) :: iscale = 10000000000
integer, allocatable :: iini_g(:,:,:),ilives_g(:,:)
real(kind=kind_phys), allocatable :: field_out(:,:,:), field_in(:,:),field_smooth(:,:)
real(kind=kind_phys), allocatable :: CA(:,:),CA1(:,:),CA2(:,:),CA3(:,:)
real(kind=kind_phys), allocatable :: noise1D(:),vertvelhigh(:,:),noise(:,:,:)
real(kind=kind_phys) :: psum,csum,CAmean,sq_diff,CAstdv
real(kind=kind_phys) :: Detmax(nca),Detmin(nca),Detmean(nca)
logical,save         :: block_message=.true.


!nca         :: switch for number of cellular automata to be used.
!ca_global   :: switch for global cellular automata
!ca_sgs      :: switch for cellular automata conditioned on SGS perturbed vertvel. 
!nfracseed   :: switch for number of random cells initially seeded
!nlives      :: switch for maximum number of lives a cell can have
!nspinup     :: switch for number of itterations to spin up the ca
!ncells      :: switch for higher resolution grid e.g ncells=4 
!               gives 4x4 times the FV3 model grid resolution.                
!ca_smooth   :: switch to smooth the cellular automata
!nthresh     :: threshold of perturbed vertical velocity used in case of sgs


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

! if(ca_global == .true. .and. ca_sgs == .true.)then
! write(0,*)'Namelist options ca_global and ca_sgs cannot both be true - exiting'
! stop
! endif

! if(ca_sgs == .true. .and. ca_smooth == .true.)then
! write(0,*)'Currently ca_smooth does not work with ca_sgs - exiting'
! stop
! endif


 call atmosphere_resolution (nlon, nlat, global=.false.)
 isize=nlon+2*halo
 jsize=nlat+2*halo
 !nlon,nlat is the compute domain - without haloes       
 !mlon,mlat is the cubed-sphere tile size. 

 inci=ncells
 incj=ncells
 
 nxc=nlon*ncells
 nyc=nlat*ncells
 
 nxch=nxc+2*halo
 nych=nyc+2*halo


 !Allocate fields:

 allocate(field_in(nlon*nlat,1))
 allocate(field_out(isize,jsize,1))
 allocate(field_smooth(nlon,nlat))
 allocate(iini_g(nxc,nyc,nca))
 allocate(ilives_g(nxc,nyc))
 allocate(vertvelhigh(nxc,nyc))
 allocate(CA(nlon,nlat))
 allocate(CA1(nlon,nlat))
 allocate(CA2(nlon,nlat))
 allocate(CA3(nlon,nlat))
 allocate(noise(nxc,nyc,nca))
 allocate(noise1D(nxc*nyc))
  
 !Initialize:
 
 noise(:,:,:) = 0.0 
 noise1D(:) = 0.0
 iini_g(:,:,:) = 0
 ilives_g(:,:) = 0
 CA1(:,:) = 0.0
 CA2(:,:) = 0.0 
 CA3(:,:) = 0.0
 
!Put the blocks of model fields into a 2d array
 levs=nlev
 blocksz=blocksize

 call atmosphere_control_data(isc, iec, jsc, jec, levs)
 call define_blocks_packed('cellular_automata', Atm_block, isc, iec, jsc, jec, levs, &
                              blocksz, block_message)

  isc = Atm_block%isc
  iec = Atm_block%iec
  jsc = Atm_block%jsc
  jec = Atm_block%jec 

!Generate random number, following stochastic physics code:
do nf=1,nca

  if (iseed_ca == 0) then
    ! generate a random seed from system clock and ens member number
    call system_clock(count, count_rate, count_max)
    ! iseed is elapsed time since unix epoch began (secs)
    ! truncate to 4 byte integer
    count_trunc = iscale*(count/iscale)
    count4 = count - count_trunc + nf*201 
  else
    ! don't rely on compiler to truncate integer(8) to integer(4) on
    ! overflow, do wrap around explicitly.
    count4 = mod(iseed_ca + 2147483648, 4294967296) - 2147483648 + nf*201
  endif

 !Set seed (to be the same) on all tasks. Save random state.
 call random_setseed(count4,rstate)
 call random_number(noise1D,rstate)
 !Put on 2D:
 do j=1,nyc
  do i=1,nxc
  noise(i,j,nf)=noise1D(i+(j-1)*nxc)
  enddo
 enddo

!Initiate the cellular automaton with random numbers larger than nfracseed
 
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

  call update_cells_global(kstep,nca,nxc,nyc,nxch,nych,nlon,nlat,CA,iini_g,ilives_g, &
                   nlives, ncells, nfracseed, nseed,nthresh, nspinup,nf) 

   
if (ca_smooth) then

do nn=1,nsmooth !number of itterations for the smoothing. 

field_in=0.

!get halo
do j=1,nlat
 do i=1,nlon
 field_in(i+(j-1)*nlon,1)=CA(i,j)
 enddo
enddo

field_out=0.

call atmosphere_scalar_field_halo(field_out,halo,isize,jsize,k_in,field_in)

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

Detmax(1)=maxval(CA)
call mp_reduce_max(Detmax(1))
Detmin(1)=minval(CA)
call mp_reduce_min(Detmin(1))

do j=1,nlat
 do i=1,nlon
    CA(i,j) = ((CA(i,j) - Detmin(1))/(Detmax(1)-Detmin(1)))
 enddo
enddo

!mean:
CAmean=0.
psum=0.
csum=0.
do j=1,nlat
 do i=1,nlon
  psum=psum+(CA(i,j))
  csum=csum+1
 enddo
enddo

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

!Put back into blocks 1D array to be passed to physics
!or diagnostics output
if(kstep < 1)then
CA(:,:)=1.
endif

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
     Diag(blk)%ca1(ix)=CA1(i,j)
     Diag(blk)%ca2(ix)=CA2(i,j)
     Diag(blk)%ca3(ix)=CA3(i,j)
     Coupling(blk)%ca1(ix)=CA1(i,j)
     Coupling(blk)%ca2(ix)=CA2(i,j)
     Coupling(blk)%ca3(ix)=CA3(i,j)
  enddo
  enddo


 deallocate(field_in)
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

end subroutine cellular_automata_global
