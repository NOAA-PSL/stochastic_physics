module cellular_automata_sgs_mod

implicit none

contains

subroutine cellular_automata_sgs(kstep,ugrs,qgrs,pgr,vvl,prsl,condition_cpl, &
            ca_deep_cpl,ca_turb_cpl,ca_shal_cpl,ca_deep_diag,ca_turb_diag,ca_shal_diag,domain_for_coupler, &
            nblks,isc,iec,jsc,jec,npx,npy,nlev, &
            nca,ncells,nlives,nfracseed,nseed,nthresh,ca_global,ca_sgs,iseed_ca, &
            ca_smooth,nspinup,blocksize,mpiroot, mpicomm)

use kinddef,           only: kind_phys
use update_ca,         only: update_cells_sgs, update_cells_global
use mersenne_twister,  only: random_setseed,random_gauss,random_stat,random_number
use mpp_domains_mod,   only: domain2D
use block_control_mod, only: block_control_type, define_blocks_packed
use mpi_wrapper,       only: mype,mp_reduce_sum,mp_bcst,mp_reduce_max,mp_reduce_min, &
                             mpi_wrapper_initialize


implicit none

!L.Bengtsson, 2017-06

!This program evolves a cellular automaton uniform over the globe given
!the flag ca_global, if instead ca_sgs is .true. it evolves a cellular automata conditioned on
!perturbed grid-box mean field. The perturbations to the mean field are given by a
!stochastic gaussian skewed (SGS) distribution.

!If ca_global is .true. it weighs the number of ca (nca) together to produce 1 output pattern
!If instead ca_sgs is given, it produces nca ca:
! 1 CA_DEEP = deep convection
! 2 CA_SHAL = shallow convection
! 3 CA_TURB = turbulence

!PLEASE NOTE: This is considered to be version 0 of the cellular automata code for FV3GFS, some functionally
!is missing/limited.

integer,intent(in) :: kstep,ncells,nca,nlives,nseed,iseed_ca,nspinup,mpiroot,mpicomm
real(kind=kind_phys), intent(in)    :: nfracseed,nthresh
logical,intent(in) :: ca_global, ca_sgs, ca_smooth
integer, intent(in) :: nblks,isc,iec,jsc,jec,npx,npy,nlev,blocksize
real(kind=kind_phys), intent(in)    :: ugrs(:,:,:)
real(kind=kind_phys), intent(in)    :: qgrs(:,:,:)
real(kind=kind_phys), intent(in)    :: pgr(:,:)
real(kind=kind_phys), intent(in)    :: vvl(:,:,:)
real(kind=kind_phys), intent(in)    :: prsl(:,:,:)
real(kind=kind_phys), intent(inout) :: condition_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_deep_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_turb_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_shal_cpl(:,:)
real(kind=kind_phys), intent(out)   :: ca_deep_diag(:,:)
real(kind=kind_phys), intent(out)   :: ca_turb_diag(:,:)
real(kind=kind_phys), intent(out)   :: ca_shal_diag(:,:)
type(domain2D), intent(inout)       :: domain_for_coupler

type(block_control_type)          :: Atm_block
type(random_stat) :: rstate
integer :: nlon, nlat, isize,jsize,nf,nn
integer :: inci, incj, nxc, nyc, nxch, nych
integer :: halo, k_in, i, j, k
integer :: seed, ierr7,blk, ix, iix, count4,ih,jh
integer :: blocksz,levs,k350,k850
integer(8) :: count, count_rate, count_max, count_trunc
integer(8) :: iscale = 10000000000
integer, allocatable :: iini(:,:,:),ilives(:,:,:),iini_g(:,:,:),ilives_g(:,:),ca_plumes(:,:)
real(kind=kind_phys), allocatable :: field_out(:,:,:), field_in(:,:),field_smooth(:,:),Detfield(:,:,:)
real(kind=kind_phys), allocatable :: omega(:,:,:),pressure(:,:,:),cloud(:,:),humidity(:,:),uwind(:,:)
real(kind=kind_phys), allocatable :: vertvelsum(:,:),vertvelmean(:,:),dp(:,:,:),surfp(:,:),shalp(:,:),gamt(:,:)
real(kind=kind_phys), allocatable :: CA(:,:),condition(:,:),rho(:,:),conditiongrid(:,:)
real(kind=kind_phys), allocatable :: CA_DEEP(:,:),CA_TURB(:,:),CA_SHAL(:,:)
real(kind=kind_phys), allocatable :: noise1D(:),vertvelhigh(:,:),noise(:,:,:)
real(kind=kind_phys) :: psum,csum,CAmean,sq_diff,CAstdv,count1,lambda
real(kind=kind_phys) :: Detmax(nca),Detmin(nca),Detmean(nca),phi,stdev,delt,condmax
logical,save         :: block_message=.true.
logical              :: nca_plumes

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
!nca_plumes   :: compute number of CA-cells ("plumes") within a NWP gridbox.

! Initialize MPI and OpenMP
if (kstep==0) then
   call mpi_wrapper_initialize(mpiroot,mpicomm)
end if

halo=1
k_in=1

if (nlev .EQ. 64) then
   k350=29
   k850=13
elseif (nlev .EQ. 127) then
   k350=61
   k850=28
else ! make a guess
   k350=int(nlev/2)
   k850=int(nlev/5)
   print*,'this level selection is not supported, making an approximation for k350 and k850'
endif

nca_plumes = .true.
!----------------------------------------------------------------------------
! Get information about the compute domain, allocate fields on this
! domain

! Some security checks for namelist combinations:
 if(nca > 5)then
 write(0,*)'Namelist option nca cannot be larger than 5 - exiting'
 stop
 endif

 nlon=iec-isc+1
 nlat=jec-jsc+1
 isize=nlon+2*halo
 jsize=nlat+2*halo

 inci=ncells
 incj=ncells

 nxc=nlon*ncells
 nyc=nlat*ncells

 nxch=nxc+2*halo
 nych=nyc+2*halo

 !Allocate fields:

 allocate(cloud(nlon,nlat))
 allocate(omega(nlon,nlat,nlev))
 allocate(pressure(nlon,nlat,nlev))
 allocate(humidity(nlon,nlat))
 allocate(uwind(nlon,nlat))
 allocate(dp(nlon,nlat,nlev))
 allocate(rho(nlon,nlat))
 allocate(surfp(nlon,nlat))
 allocate(vertvelmean(nlon,nlat))
 allocate(vertvelsum(nlon,nlat))
 allocate(field_in(nlon*nlat,1))
 allocate(field_out(isize,jsize,1))
 allocate(field_smooth(nlon,nlat))
 allocate(iini(nxc,nyc,nca))
 allocate(ilives(nxc,nyc,nca))
 allocate(iini_g(nxc,nyc,nca))
 allocate(ilives_g(nxc,nyc))
 allocate(vertvelhigh(nxc,nyc))
 allocate(condition(nxc,nyc))
 allocate(conditiongrid(nlon,nlat))
 allocate(shalp(nlon,nlat))
 allocate(gamt(nlon,nlat))
 allocate(Detfield(nlon,nlat,nca))
 allocate(CA(nlon,nlat))
 allocate(ca_plumes(nlon,nlat))
 allocate(CA_TURB(nlon,nlat))
 allocate(CA_DEEP(nlon,nlat))
 allocate(CA_SHAL(nlon,nlat))
 allocate(noise(nxc,nyc,nca))
 allocate(noise1D(nxc*nyc))

 !Initialize:
 Detfield(:,:,:)=0.
 vertvelmean(:,:) =0.
 vertvelsum(:,:)=0.
 cloud(:,:)=0.
 humidity(:,:)=0.
 uwind(:,:) = 0.
 condition(:,:)=0.
 conditiongrid(:,:)=0.
 vertvelhigh(:,:)=0.
 ca_plumes(:,:) = 0
 noise(:,:,:) = 0.0
 noise1D(:) = 0.0
 iini(:,:,:) = 0
 ilives(:,:,:) = 0
 iini_g(:,:,:) = 0
 ilives_g(:,:) = 0
 Detmax(:)=0.
 Detmin(:)=0.

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
      uwind(i,j)         = ugrs(blk,ix,k350)
      conditiongrid(i,j) = condition_cpl(blk,ix)
      shalp(i,j)         = ca_shal_cpl(blk,ix)
      gamt(i,j)          = ca_turb_cpl(blk,ix)
      surfp(i,j)         = pgr(blk,ix)
      humidity(i,j)      = qgrs(blk,ix,k850) !about 850 hpa
      do k = 1,k350 !Lower troposphere
      omega(i,j,k)       = vvl(blk,ix,k) ! layer mean vertical velocity in pa/sec
      pressure(i,j,k)    = prsl(blk,ix,k) ! layer mean pressure in Pa
      enddo
  enddo
 enddo

 CA_TURB(:,:) = 0.0
 CA_DEEP(:,:) = 0.0
 CA_SHAL(:,:) = 0.0

!Compute layer averaged vertical velocity (Pa/s)
 vertvelsum=0.
 vertvelmean=0.
 do j=1,nlat
  do i =1,nlon
    dp(i,j,1)=(surfp(i,j)-pressure(i,j,1))
    do k=2,k350
     dp(i,j,k)=(pressure(i,j,k-1)-pressure(i,j,k))
    enddo
    count1=0.
    do k=1,k350
     count1=count1+1.
     vertvelsum(i,j)=vertvelsum(i,j)+(omega(i,j,k)*dp(i,j,k))
   enddo
  enddo
 enddo

 do j=1,nlat
  do i=1,nlon
   vertvelmean(i,j)=vertvelsum(i,j)/(surfp(i,j)-pressure(i,j,k350))
  enddo
 enddo

!Generate random number, following stochastic physics code:
!if(kstep==2) then
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


!Initiate the cellular automaton with random numbers larger than nfracseed

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
!endif ! kstep=0

!In case we want to condition the cellular automaton on a large scale field
!we here set the "condition" variable to a different model field depending
!on nf. (this is not used if ca_global = .true.)


do nf=1,nca !update each ca


  if(nf==1)then
  inci=ncells
  incj=ncells
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

   do j = 1,nyc
    do i = 1,nxc
      ilives(i,j,nf)=real(nlives)*(condition(i,j)/condmax)
    enddo
   enddo

 elseif(nf==2)then
  inci=ncells
  incj=ncells
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

   do j = 1,nyc
    do i = 1,nxc
      ilives(i,j,nf)=real(nlives)*(condition(i,j)/condmax)
    enddo
   enddo

 else

  inci=ncells
  incj=ncells
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

  do j = 1,nyc
    do i = 1,nxc
      ilives(i,j,nf)=real(nlives)*(condition(i,j)/condmax)
    enddo
   enddo

 endif !nf


!Vertical velocity has its own variable in order to condition on combination
!of "condition" and vertical velocity.

  inci=ncells
  incj=ncells
  do j=1,nyc
   do i=1,nxc
     vertvelhigh(i,j)=vertvelmean(inci/ncells,incj/ncells)
     if(i.eq.inci)then
     inci=inci+ncells
     endif
   enddo
   inci=ncells
   if(j.eq.incj)then
   incj=incj+ncells
   endif
  enddo

!Calculate neighbours and update the automata
!If ca-global is used, then nca independent CAs are called and weighted together to create one field; CA

  call update_cells_sgs(kstep,nca,nxc,nyc,nxch,nych,nlon,nlat,isc,iec,jsc,jec, &
                   npx,npy,domain_for_coupler,CA,ca_plumes,iini,ilives,        &
                   nlives,ncells,nfracseed,nseed,nthresh,nspinup,nf,nca_plumes)

   if(nf==1)then
    CA_DEEP(:,:)=CA(:,:)
   elseif(nf==2)then
    CA_SHAL(:,:)=CA(:,:)
   else
    CA_TURB(:,:)=CA(:,:)
   endif


 enddo !nf (nca)

!!Post-processesing - could be made into a separate sub-routine

!Deep convection ====================================

if(kstep > 1)then

!Use min-max method to normalize range
Detmax(1)=maxval(CA_DEEP,CA_DEEP.NE.0.)
call mp_reduce_max(Detmax(1))
Detmin(1)=minval(CA_DEEP,CA_DEEP.NE.0.)
call mp_reduce_min(Detmin(1))


do j=1,nlat
 do i=1,nlon
 if(CA_DEEP(i,j).NE.0.)then
    CA_DEEP(i,j) =(CA_DEEP(i,j) - Detmin(1))/(Detmax(1)-Detmin(1))
 endif
 enddo
enddo

!Compute the mean of the new range and subtract
CAmean=0.
psum=0.
csum=0.
do j=1,nlat
 do i=1,nlon
  if(CA_DEEP(i,j).NE.0.)then
  psum=psum+(CA_DEEP(i,j))
  csum=csum+1
  endif
 enddo
enddo

call mp_reduce_sum(psum)
call mp_reduce_sum(csum)

CAmean=psum/csum

do j=1,nlat
 do i=1,nlon
 if(CA_DEEP(i,j).NE.0.)then
  CA_DEEP(i,j)=(CA_DEEP(i,j)-CAmean)
 endif
 enddo
enddo

Detmin(1) = minval(CA_DEEP,CA_DEEP.NE.0)
call mp_reduce_min(Detmin(1))

!Shallow convection ============================================================

!Use min-max method to normalize range
Detmax(2)=maxval(CA_SHAL,CA_SHAL.NE.0)
call mp_reduce_max(Detmax(2))
Detmin(2)=minval(CA_SHAL,CA_SHAL.NE.0)
call mp_reduce_min(Detmin(2))

do j=1,nlat
 do i=1,nlon
 if(CA_SHAL(i,j).NE.0.)then
    CA_SHAL(i,j)=(CA_SHAL(i,j) - Detmin(2))/(Detmax(2)-Detmin(2))
 endif
 enddo
enddo

!Compute the mean of the new range and subtract
CAmean=0.
psum=0.
csum=0.
do j=1,nlat
 do i=1,nlon
  if(CA_SHAL(i,j).NE.0.)then
  psum=psum+(CA_SHAL(i,j))
  csum=csum+1
  endif
 enddo
enddo

call mp_reduce_sum(psum)
call mp_reduce_sum(csum)

CAmean=psum/csum

do j=1,nlat
 do i=1,nlon
 if(CA_SHAL(i,j).NE.0.)then
 CA_SHAL(i,j)=(CA_SHAL(i,j)-CAmean)
 endif
 enddo
enddo

!Turbulence =============================================================================

!Use min-max method to normalize range
Detmax(3)=maxval(CA_TURB,CA_TURB.NE.0)
call mp_reduce_max(Detmax(3))
Detmin(3)=minval(CA_TURB,CA_TURB.NE.0)
call mp_reduce_min(Detmin(3))

do j=1,nlat
 do i=1,nlon
 if(CA_TURB(i,j).NE.0.)then
    CA_TURB(i,j)=(CA_TURB(i,j) - Detmin(3))/(Detmax(3)-Detmin(3))
 endif
 enddo
enddo

!Compute the mean of the new range and subtract
CAmean=0.
psum=0.
csum=0.
do j=1,nlat
 do i=1,nlon
  if(CA_TURB(i,j).NE.0.)then
  psum=psum+(CA_TURB(i,j))
  csum=csum+1
  endif
 enddo
enddo

call mp_reduce_sum(psum)
call mp_reduce_sum(csum)

CAmean=psum/csum

do j=1,nlat
 do i=1,nlon
 if(CA_TURB(i,j).NE.0.)then
 CA_TURB(i,j)=(CA_TURB(i,j)-CAmean)
 endif
 enddo
enddo

endif !kstep >1

do j=1,nlat
 do i=1,nlon
    if(conditiongrid(i,j) == 0)then
     CA_DEEP(i,j)=0.
     ca_plumes(i,j)=0.
   endif
 enddo
enddo

if(kstep == 1)then
do j=1,nlat
 do i=1,nlon
   ca_plumes(i,j)=0.
 enddo
enddo
endif

!Put back into blocks 1D array to be passed to physics
!or diagnostics output

  do blk = 1, Atm_block%nblks
  do ix = 1,Atm_block%blksz(blk)
     i = Atm_block%index(blk)%ii(ix) - isc + 1
     j = Atm_block%index(blk)%jj(ix) - jsc + 1
     ca_deep_diag(blk,ix)=ca_plumes(i,j)
     ca_turb_diag(blk,ix)=conditiongrid(i,j)
     ca_shal_diag(blk,ix)=CA_SHAL(i,j)
     ca_deep_cpl(blk,ix)=ca_plumes(i,j)
     ca_turb_cpl(blk,ix)=CA_TURB(i,j)
     ca_shal_cpl(blk,ix)=CA_SHAL(i,j)
  enddo
  enddo


 deallocate(omega)
 deallocate(pressure)
 deallocate(humidity)
 deallocate(dp)
 deallocate(conditiongrid)
 deallocate(shalp)
 deallocate(gamt)
 deallocate(rho)
 deallocate(surfp)
 deallocate(vertvelmean)
 deallocate(vertvelsum)
 deallocate(field_in)
 deallocate(field_out)
 deallocate(field_smooth)
 deallocate(iini)
 deallocate(ilives)
 deallocate(condition)
 deallocate(Detfield)
 deallocate(CA)
 deallocate(ca_plumes)
 deallocate(CA_TURB)
 deallocate(CA_DEEP)
 deallocate(CA_SHAL)
 deallocate(noise)
 deallocate(noise1D)

end subroutine cellular_automata_sgs

end module cellular_automata_sgs_mod
