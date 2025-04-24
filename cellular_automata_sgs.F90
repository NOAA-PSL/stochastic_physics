module cellular_automata_sgs_mod

use update_ca, only : domain_global,domain_sgs,iscnx,iecnx,jscnx,jecnx,isdnx,iednx,jsdnx,jednx,nxncells,nyncells,cold_start_ca_sgs
implicit none


contains

subroutine cellular_automata_sgs(kstep,dtf,restart,first_time_step,sst,lsmsk,lake,uwind,vwind,height,dx,condition_cpl, &
            ca_deep_cpl,ca_turb_cpl,ca_shal_cpl,domain_in, &
            nblks,isc,iec,jsc,jec,npx,npy,nlev,nthresh, mytile, &
            nca,ncells,nlives,nfracseed,nseed,iseed_ca,ca_advect, &
            nspinup,ca_trigger,blocksize,mpiroot,mpicomm)

use mpi_f08
use kinddef,           only: kind_phys,kind_dbl_prec
use update_ca,         only: update_cells_sgs, define_ca_domain
use random_numbers,    only: random_01_CB
use mpp_domains_mod,   only: domain2D,mpp_get_global_domain,CENTER, mpp_get_data_domain, mpp_get_compute_domain,&
                             mpp_define_io_domain,mpp_get_io_domain_layout,mpp_get_ntile_count
use block_control_mod, only: block_control_type, define_blocks_packed
use time_manager_mod, only: time_type
use mpi_wrapper,       only: mype,mp_reduce_max, &
                             mpi_wrapper_initialize



implicit none
!L.Bengtsson, 2017-06

!L.Bengtsson, 2021-05
!Significant cleaning of old ideas
!Inclusion of restart capability
!Setting control variables as a function of dx and dt
!for scale adaptation. 

!P.Pegion, 2021-09
! swtich to new random number generator and improve computational efficiency
! and remove unsued code

!L.Bengtsson, 2023-05
!Add horizontal advection of CA cells

!This routine produces an output field CA_DEEP for coupling to convection (saSAS).
!CA_DEEP can be either number of plumes in a cluster (nca_plumes=true) or updraft 
!area fraction (nca_plumes=false)

integer,intent(in) :: kstep,ncells,nca,nlives,nseed,nspinup,mpiroot,mytile
type(MPI_Comm),intent(in) :: mpicomm
integer(kind=kind_dbl_prec),           intent(in)    :: iseed_ca
real(kind=kind_phys), intent(in)    :: nfracseed,dtf,nthresh
logical,intent(in) :: restart,ca_trigger,first_time_step,ca_advect
integer, intent(in) :: nblks,isc,iec,jsc,jec,npx,npy,nlev,blocksize
real(kind=kind_phys), intent(in)    :: sst(:,:),lsmsk(:,:),lake(:,:),uwind(:,:,:),dx(:,:),vwind(:,:,:),height(:,:,:)
real(kind=kind_phys), intent(inout) :: condition_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_deep_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_turb_cpl(:,:)
real(kind=kind_phys), intent(inout) :: ca_shal_cpl(:,:)
type(domain2D),       intent(inout) :: domain_in

type(block_control_type)          :: Atm_block
integer :: nlon, nlat, isize,jsize,nf,nn
integer :: inci, incj, nxc, nyc, nxch, nych, nx, ny
integer :: halo, k_in, i, j, k
integer :: seed, ierr7,blk, ix, iix, count4,ih,jh
integer :: blocksz,levs,u200,u850
integer, save :: initialize_ca
integer(8) :: count, count_rate, count_max, count_trunc,nx_full
integer(8) :: iscale = 10000000000_8
integer, allocatable :: iini(:,:,:),ilives_in(:,:,:),ca_plumes(:,:),io_layout(:)
real(kind=kind_phys), allocatable :: ssti(:,:),lsmski(:,:),lakei(:,:),uwindi(:,:),vwindi(:,:),dxi(:,:),heighti(:,:,:)
real(kind=kind_phys), allocatable :: CA(:,:),condition(:,:),uhigh(:,:),vhigh(:,:),dxhigh(:,:),conditiongrid(:,:)
real(kind=kind_phys), allocatable :: CA_DEEP(:,:),zi(:,:,:),sumx(:,:)
real*8              , allocatable :: noise(:,:,:)
real(kind=kind_phys) :: condmax,condmaxinv,livesmax,livesmaxinv,factor,pi,re
logical,save         :: block_message=.true.
logical              :: nca_plumes
logical,save         :: first_flag
integer*8            :: i1,j1
integer              :: ct,ntiles
real                 :: dz,invgrav

!nca         :: switch for number of cellular automata to be used.
!            :: for the moment only 1 CA can be used 
!nfracseed   :: switch for number of random cells initially seeded
!nlives      :: switch for time scale (s)
!nspinup     :: switch for number of itterations to spin up the ca
!ncells      :: switch for CA cell size (m)
!nca_plumes   :: compute number of CA-cells ("plumes") within a NWP gridbox.

if (nca .LT. 1) return
! Initialize MPI and OpenMP
if (first_time_step) then
   call mpi_wrapper_initialize(mpiroot,mpicomm)
end if

halo=3
k_in=1
!Right now these values are experimental:                                                                                                                                                                                     
u200=56
u850=13
!Gravitational acceleration: 
invgrav=1./9.81

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
 call mpp_get_global_domain(domain_in,xsize=nx,ysize=ny,position=CENTER)
 ntiles = mpp_get_ntile_count(domain_in)

 if(mype == 1)then
 write(*,*)'ncells=',ncells
 write(*,*)'nlives=',nlives
 write(*,*)'nthresh=',nthresh
 endif

 inci=ncells
 incj=ncells

!--- get params from domain_sgs for building board and board_halo                                                                                

 if(first_time_step)then 
  !Get CA domain                                                                                                                                       
  if (.not.restart) then
      allocate(io_layout(2))
      io_layout=mpp_get_io_domain_layout(domain_in)
      call define_ca_domain(domain_in,domain_sgs,halo,ncells,nxncells,nyncells)
      call mpp_define_io_domain(domain_sgs, io_layout)
  endif
  call mpp_get_data_domain    (domain_sgs,isdnx,iednx,jsdnx,jednx)
  call mpp_get_compute_domain (domain_sgs,iscnx,iecnx,jscnx,jecnx)
 endif

  nxc = iecnx-iscnx+1
  nyc = jecnx-jscnx+1
  nxch = iednx-isdnx+1
  nych = jednx-jsdnx+1
 
 !Allocate fields:

 allocate(ssti(nlon,nlat))
 allocate(lsmski(nlon,nlat))
 allocate(lakei(nlon,nlat))
 allocate(uwindi(nlon,nlat))
 allocate(vwindi(nlon,nlat))
 allocate(heighti(nlon,nlat,nlev))
 allocate(zi(nlon,nlat,nlev))
 allocate(sumx(nlon,nlat))
 allocate(dxi(nlon,nlat))
 allocate(iini(nxc,nyc,nca))
 allocate(ilives_in(nxc,nyc,nca))
 allocate(condition(nxc,nyc))
 allocate(uhigh(nxc,nyc))
 allocate(dxhigh(nxc,nyc))
 allocate(vhigh(nxc,nyc))
 allocate(conditiongrid(nlon,nlat))
 allocate(CA(nlon,nlat))
 allocate(ca_plumes(nlon,nlat))
 allocate(CA_DEEP(nlon,nlat))

 !Initialize:
 ilives_in(:,:,:) = 0
 iini(:,:,:) = 0
 sumx(:,:) = 0.
 uwindi(:,:) = 0.
 vwindi(:,:) =0.

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
      dxi(i,j)           = dx(blk,ix)
      do k = 1,levs
         heighti(i,j,k) = height(blk,ix,k)*invgrav
      enddo
      do k = 2,levs
         zi(i,j,k)=0.5*(heighti(i,j,k)+heighti(i,j,k-1))
      enddo
      do k = u850,u200
         dz=zi(i,j,k)-zi(i,j,k-1)
         uwindi(i,j)   = uwindi(i,j) + uwind(blk,ix,k)*dz
         vwindi(i,j)   = vwindi(i,j) + vwind(blk,ix,k)*dz
         sumx(i,j) = sumx(i,j) + dz
      enddo
  enddo
 enddo

 do blk = 1,Atm_block%nblks
  do ix = 1, Atm_block%blksz(blk)
      i = Atm_block%index(blk)%ii(ix) - isc + 1
      j = Atm_block%index(blk)%jj(ix) - jsc + 1
      uwindi(i,j)=uwindi(i,j)/sumx(i,j)
      vwindi(i,j)=vwindi(i,j)/sumx(i,j)
   enddo
 enddo
         

!Initialize the CA when the condition field is populated
  do j=1,nyc
   do i=1,nxc
     condition(i,j)=conditiongrid(inci/ncells,incj/ncells)
     uhigh(i,j)=uwindi(inci/ncells,incj/ncells)
     vhigh(i,j)=vwindi(inci/ncells,incj/ncells)
     dxhigh(i,j)=dxi(inci/ncells,incj/ncells)/ncells !dx on the finer grid
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
  condmaxinv=0.0
  call mp_reduce_max(condmax)
  if(condmax > 0.)then
     if(.not. first_flag)then
        first_flag = .true.
        initialize_ca = kstep
     endif
     condmaxinv=1.0/condmax
  endif

if(kstep >=initialize_ca)then
  do nf=1,nca
     do j = 1,nyc
        do i = 1,nxc
           ilives_in(i,j,nf)=int(real(nlives)*(condition(i,j)*condmaxinv))
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
if (cold_start_ca_sgs) then
   if(kstep == initialize_ca) then
      nx_full=int(ncells,kind=8)*int(npx-1,kind=8)
      allocate(noise(nxc,nyc,nca))
      do j=1,nyc
         j1=j+(jsc-1)*ncells
         do i=1,nxc
            i1=i+(isc-1)*ncells
            if (iseed_ca <= 0) then
               ! generate a random seed from system clock and ens member number
               call system_clock(count, count_rate, count_max)
               ! iseed is elapsed time since unix epoch began (secs)
               ! truncate to 4 byte integer
               count_trunc = iscale*(count/iscale)
               count4 = count - count_trunc + mytile *( i1+nx_full*(j1-1)) ! no need to multply by 7 since time will be different in sgs
            else
               ! don't rely on compiler to truncate integer(8) to integer(4) on
               ! overflow, do wrap around explicitly.
               count4 = int(mod(int((iseed_ca+mytile)*(i1+nx_full*(j1-1)), 8) + 2147483648_8, 4294967296_8) - 2147483648_8)
            endif
            ct=1
            do nf=1,nca
               noise(i,j,nf)=real(random_01_CB(ct,count4),kind=8)
               ct=ct+1
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
   
   deallocate(noise)
   endif ! 
endif !  cold_start_ca_sgs

!Calculate neighbours and update the automata
 do nf=1,nca
  call update_cells_sgs(kstep,halo,dtf,initialize_ca,iseed_ca,first_flag,restart,first_time_step,nca,nxc,nyc, &
                        nxch,nych,nlon,nlat,isc,iec,jsc,jec,ca_advect, &
                        npx,npy,CA,ca_plumes,iini,ilives_in,uhigh,vhigh,dxhigh, &
                        nlives,nfracseed,nseed,nspinup,nf,nca_plumes,ncells,mytile)

    if(nca_plumes)then
    do j=1,nlat
       do i=1,nlon
          CA_DEEP(i,j)=ca_plumes(i,j)
       enddo
    enddo
    else
    livesmax=maxval(ilives_in)
    call mp_reduce_max(livesmax)
    livesmaxinv=1.0/livesmax
    do j=1,nlat
       do i=1,nlon
          CA_DEEP(i,j)=CA(i,j)*livesmaxinv
       enddo
    enddo
    endif

 enddo !nf (nca)

!Limit CA activity to the Tropical Ocean

if(ntiles==6)then 
do j=1,nlat
   do i=1,nlon
      if(ssti(i,j) < 300. .or. lsmski(i,j) /= 0. .or. lakei(i,j) > 0.0)then
      CA_DEEP(i,j)=0.
      endif
   enddo
enddo
endif

!Put back into blocks 1D array to be passed to physics
!or diagnostics output

 do blk = 1, Atm_block%nblks
    do ix = 1,Atm_block%blksz(blk)
       i = Atm_block%index(blk)%ii(ix) - isc + 1
       j = Atm_block%index(blk)%jj(ix) - jsc + 1
       ca_deep_cpl(blk,ix)=CA_DEEP(i,j) 
    enddo
 enddo

 deallocate(conditiongrid)
 deallocate(ssti)
 deallocate(lsmski)
 deallocate(lakei)
 deallocate(iini)
 deallocate(ilives_in)
 deallocate(uwindi)
 deallocate(vwindi)
 deallocate(uhigh)
 deallocate(vhigh)
 deallocate(dxhigh)
 deallocate(condition)
 deallocate(CA)
 deallocate(ca_plumes)
 deallocate(CA_DEEP)

end subroutine cellular_automata_sgs

end module cellular_automata_sgs_mod
