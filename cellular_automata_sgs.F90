module cellular_automata_sgs_mod

implicit none

contains

subroutine cellular_automata_sgs(kstep,dtf,restart,first_time_step,sst,lsmsk,lake,condition_cpl, &
            ca_deep_cpl,ca_turb_cpl,ca_shal_cpl,ca_deep_diag,ca_turb_diag,ca_shal_diag,domain, &
            nblks,isc,iec,jsc,jec,npx,npy,nlev,nthresh,rcell, &
            nca,scells,tlives,nfracseed,nseed,ca_global,ca_sgs,iseed_ca, &
            ca_smooth,nspinup,ca_trigger,blocksize,mpiroot,mpicomm)

use kinddef,           only: kind_phys
use update_ca,         only: update_cells_sgs, update_cells_global, define_ca_domain
use mersenne_twister,  only: random_setseed,random_gauss,random_stat,random_number
use mpp_domains_mod,   only: domain2D
use block_control_mod, only: block_control_type, define_blocks_packed
use time_manager_mod, only: time_type
use mpi_wrapper,       only: mype,mp_reduce_sum,mp_bcst,mp_reduce_max,mp_reduce_min, &
                             mpi_wrapper_initialize
use mpp_domains_mod
use mpp_mod



implicit none
!L.Bengtsson, 2017-06

!L.Bengtsson, 2021-05
!Significant cleaning of old ideas
!Inclusion of restart capability
!Setting control variables as a function of dx and dt
!for scale adaptation. 

!This routine produces an output field CA_DEEP for coupling to convection (saSAS).
!CA_DEEP can be either number of plumes in a cluster (nca_plumes=true) or updraft 
!area fraction (nca_plumes=false)

integer,intent(in) :: kstep,scells,nca,tlives,nseed,iseed_ca,nspinup,mpiroot,mpicomm
real(kind=kind_phys), intent(in)    :: nfracseed,dtf,rcell
logical,intent(in) :: ca_global, ca_sgs, ca_smooth, restart,ca_trigger,first_time_step
integer, intent(in) :: nblks,isc,iec,jsc,jec,npx,npy,nlev,blocksize
real  , intent(out) :: nthresh
real(kind=kind_phys), intent(in)    :: sst(:,:),lsmsk(:,:),lake(:,:)
real(kind=kind_phys), intent(inout) :: condition_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_deep_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_turb_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_shal_cpl(:,:)
real(kind=kind_phys), intent(out)   :: ca_deep_diag(:,:)
real(kind=kind_phys), intent(out)   :: ca_turb_diag(:,:)
real(kind=kind_phys), intent(out)   :: ca_shal_diag(:,:)
type(domain2D),       intent(inout) :: domain

type(block_control_type)          :: Atm_block
type(random_stat) :: rstate
integer :: nlon, nlat, isize,jsize,nf,nn
integer :: inci, incj, nxc, nyc, nxch, nych, nx, ny
integer :: nxncells, nyncells
integer :: halo, k_in, i, j, k
integer :: seed, ierr7,blk, ix, iix, count4,ih,jh
integer :: blocksz,levs
integer :: isdnx,iednx,jsdnx,jednx
integer :: iscnx,iecnx,jscnx,jecnx
integer :: ncells,nlives
integer, save :: initialize_ca
integer(8) :: count, count_rate, count_max, count_trunc
integer(8) :: iscale = 10000000000
integer, allocatable :: iini(:,:,:),ilives_in(:,:,:),ca_plumes(:,:)
real(kind=kind_phys), allocatable :: ssti(:,:),lsmski(:,:),lakei(:,:)
real(kind=kind_phys), allocatable :: CA(:,:),condition(:,:),conditiongrid(:,:)
real(kind=kind_phys), allocatable :: CA_DEEP(:,:)
real(kind=kind_phys), allocatable :: noise1D(:),noise(:,:,:)
real(kind=kind_phys) :: condmax,livesmax,factor,dx,pi,re
type(domain2D)       :: domain_ncellx
logical,save         :: block_message=.true.
logical              :: nca_plumes
logical,save         :: first_flag

!nca         :: switch for number of cellular automata to be used.
!            :: for the moment only 1 CA can be used if ca_sgs = true
!ca_global   :: switch for global cellular automata
!ca_sgs      :: switch for cellular automata for deep convection
!nfracseed   :: switch for number of random cells initially seeded
!tlives      :: switch for time scale (s)
!nspinup     :: switch for number of itterations to spin up the ca
!scells      :: switch for CA cell size (m)
!ca_smooth   :: switch to smooth the cellular automata
!nca_plumes   :: compute number of CA-cells ("plumes") within a NWP gridbox.

! Initialize MPI and OpenMP
if (first_time_step) then
   call mpi_wrapper_initialize(mpiroot,mpicomm)
end if

halo=1
k_in=1

nca_plumes = .true.

if(first_time_step)then
   first_flag = .false.
   initialize_ca = 100000
endif

!----------------------------------------------------------------------------
! Get information about the compute domain, allocate fields on this
! domain

! Some security checks for namelist combinations:
 if(nca > 1)then
 write(0,*)'When ca_sgs=.True., nca has to be 1 - exiting'
 stop
 endif

 nlon=iec-isc+1
 nlat=jec-jsc+1
 isize=nlon+2*halo
 jsize=nlat+2*halo

 !Set time and length scales:
 call mpp_get_global_domain(domain,xsize=nx,ysize=ny,position=CENTER)
 pi=3.14159
 re=6371000.
 dx=0.5*pi*re/real(nx)
 ncells=int(dx/real(scells))
 nlives=int(real(tlives)/dtf)
 ncells = MIN(ncells,10)
 nlives = MAX(nlives,5)

 nthresh=rcell*real(ncells)*real(ncells)

 if(mype == 1)then
 write(*,*)'ncells=',ncells
 write(*,*)'nlives=',nlives
 write(*,*)'nthresh=',nthresh
 endif

 inci=ncells
 incj=ncells

!--- get params from domain_ncellx for building board and board_halo                                                                                

  !Get CA domain                                                                                                                                       
  call define_ca_domain(domain,domain_ncellx,ncells,nxncells,nyncells)
  call mpp_get_data_domain    (domain_ncellx,isdnx,iednx,jsdnx,jednx)
  call mpp_get_compute_domain (domain_ncellx,iscnx,iecnx,jscnx,jecnx)
  !write(1000+mpp_pe(),*) "nxncells,nyncells: ",nxncells,nyncells
  !write(1000+mpp_pe(),*) "iscnx,iecnx,jscnx,jecnx: ",iscnx,iecnx,jscnx,jecnx
  !write(1000+mpp_pe(),*) "isdnx,iednx,jsdnx,jednx: ",isdnx,iednx,jsdnx,jednx

  nxc = iecnx-iscnx+1
  nyc = jecnx-jscnx+1
  nxch = iednx-isdnx+1
  nych = jednx-jsdnx+1



 !Allocate fields:

 allocate(ssti(nlon,nlat))
 allocate(lsmski(nlon,nlat))
 allocate(lakei(nlon,nlat))
 allocate(iini(nxc,nyc,nca))
 allocate(ilives_in(nxc,nyc,nca))
 allocate(condition(nxc,nyc))
 allocate(conditiongrid(nlon,nlat))
 allocate(CA(nlon,nlat))
 allocate(ca_plumes(nlon,nlat))
 allocate(CA_DEEP(nlon,nlat))
 allocate(noise(nxc,nyc,nca))
 allocate(noise1D(nxc*nyc))

 !Initialize:
 condition(:,:)=0.
 conditiongrid(:,:)=0.
 ca_plumes(:,:) = 0
 noise(:,:,:) = 0.0
 noise1D(:) = 0.0
 iini(:,:,:) = 0
 ilives_in(:,:,:) = 0
 CA_DEEP(:,:) = 0.

 !Put the blocks of model fields into a 2d array - can't use nlev and blocksize directly,
 !because the arguments to define_blocks_packed are intent(inout) and not intent(in).
 levs=nlev
 blocksz=blocksize

 call define_blocks_packed('cellular_automata', Atm_block, isc, iec, jsc, jec, levs, &
                              blocksz, block_message)

 do blk = 1,Atm_block%nblks
  do ix = 1, Atm_block%blksz(blk)
      i = Atm_block%index(blk)%ii(ix) - isc + 1
      j = Atm_block%index(blk)%jj(ix) - jsc + 1
      conditiongrid(i,j) = condition_cpl(blk,ix)
      ssti(i,j)          = sst(blk,ix)
      lsmski(i,j)        = lsmsk(blk,ix)
      lakei(i,j)         = lake(blk,ix)
  enddo
 enddo

!Initialize the CA when the condition field is populated

  do j=1,nyc
   do i=1,nxc
     condition(i,j)=conditiongrid(inci/ncells,incj/ncells)
     if(i.eq.inci)then
     inci=inci+ncells
     endif
   enddo
   inci=ncells
   if(j.eq.incj)then
   incj=incj+ncells
   endif
  enddo

  condmax=maxval(condition)
  call mp_reduce_max(condmax)
  if(condmax > 0.)then
     if(.not. first_flag)then
        first_flag = .true.
        initialize_ca = kstep
     endif
  endif

if(kstep >=initialize_ca)then
  do nf=1,nca
     do j = 1,nyc
        do i = 1,nxc
           ilives_in(i,j,nf)=int(real(nlives)*(condition(i,j)/condmax))
        enddo
     enddo
  enddo

else

   do nf=1,nca
      do j = 1,nyc
         do i = 1,nxc
            ilives_in(i,j,nf)=0
         enddo
      enddo
   enddo

endif
                                                                                                                                        
!Generate random number, following stochastic physics code:
if(kstep == initialize_ca) then
   if (iseed_ca == 0) then
    ! generate a random seed from system clock and ens member number
    call system_clock(count, count_rate, count_max)
    ! iseed is elapsed time since unix epoch began (secs)
    ! truncate to 4 byte integer
    count_trunc = iscale*(count/iscale)
    count4 = count - count_trunc
  else
    ! don't rely on compiler to truncate integer(8) to integer(4) on
    ! overflow, do wrap around explicitly.
    count4 = mod(mype + iseed_ca + 2147483648, 4294967296) - 2147483648
  endif

  call random_setseed(count4)

  do nf=1,nca
    call random_number(noise1D)
    !Put on 2D:
    do j=1,nyc
      do i=1,nxc
        noise(i,j,nf)=noise1D(i+(j-1)*nxc)
      enddo
    enddo
   enddo

!Initiate the cellular automaton with random numbers larger than nfracseed
   do nf=1,nca
    do j = 1,nyc
      do i = 1,nxc
        if (noise(i,j,nf) > nfracseed ) then
          iini(i,j,nf)=1
        else
          iini(i,j,nf)=0
        endif
      enddo
    enddo
  enddo !nf

endif ! 

!Calculate neighbours and update the automata
 do nf=1,nca
  call update_cells_sgs(kstep,initialize_ca,first_flag,restart,first_time_step,iseed_ca,nca,nxc,nyc, &
                        nxch,nych,nlon,nlat,nxncells,nyncells,isc,iec,jsc,jec, &
                        npx,npy,isdnx,iednx,jsdnx,jednx,iscnx,iecnx,jscnx,jecnx,domain_ncellx,CA,ca_plumes,iini,ilives_in,        &
                        nlives,nfracseed,nseed,nspinup,nf,nca_plumes,ncells)

    if(nca_plumes)then
    do j=1,nlat
       do i=1,nlon
          CA_DEEP(i,j)=ca_plumes(i,j)
       enddo
    enddo
    else
    livesmax=maxval(ilives_in)
    call mp_reduce_max(livesmax)
    do j=1,nlat
       do i=1,nlon
          CA_DEEP(i,j)=CA(i,j)/livesmax
       enddo
    enddo
    endif

 enddo !nf (nca)

!Limit CA activity to the Tropical Ocean

do j=1,nlat
   do i=1,nlon
      if(ssti(i,j) < 300. .or. lsmski(i,j) /= 0. .or. lakei(i,j) > 0.0)then
      CA_DEEP(i,j)=0.
      endif
   enddo
enddo

!Put back into blocks 1D array to be passed to physics
!or diagnostics output

  do blk = 1, Atm_block%nblks
  do ix = 1,Atm_block%blksz(blk)
     i = Atm_block%index(blk)%ii(ix) - isc + 1
     j = Atm_block%index(blk)%jj(ix) - jsc + 1
     ca_deep_diag(blk,ix)=CA_DEEP(i,j)
     ca_deep_cpl(blk,ix)=CA_DEEP(i,j) 
  enddo
  enddo

 deallocate(conditiongrid)
 deallocate(ssti)
 deallocate(lsmski)
 deallocate(lakei)
 deallocate(iini)
 deallocate(ilives_in)
 deallocate(condition)
 deallocate(CA)
 deallocate(ca_plumes)
 deallocate(CA_DEEP)
 deallocate(noise)
 deallocate(noise1D)

end subroutine cellular_automata_sgs

end module cellular_automata_sgs_mod
