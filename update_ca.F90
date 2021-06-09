module update_ca
!Main module for evolving CA in time, also includes
!read and write restart routines, to restart fields 
!on the ncellsxncells CA grid

use halo_exchange,    only: atmosphere_scalar_field_halo
use mersenne_twister, only: random_setseed,random_gauss,random_stat,random_number
use mpi_wrapper,      only: mype,mp_reduce_sum,mp_bcst,mp_reduce_min,mp_reduce_max
use mpp_domains_mod
use mpp_mod

implicit none

public  write_ca_restart
public  read_ca_restart
public  update_cells_sgs
public  update_cells_global

integer,allocatable,save :: board(:,:,:), lives(:,:,:)

contains

!Compute CA domain:--------------------------------------------------------------------------
subroutine define_ca_domain(fv_domain,domain_ncellx,ncells,nxncells,nyncells)
implicit none

type(domain2D),intent(inout) :: fv_domain
type(domain2D),intent(inout) :: domain_ncellx
integer,intent(in)           :: ncells
integer,intent(out) :: nxncells, nyncells
integer :: halo1 = 1
integer :: layout(2)
integer :: ntiles
integer, allocatable :: pe_start(:), pe_end(:)

integer :: i, j, k, n
integer :: nx, ny
integer :: isc,iec,jsc,jec
integer :: isd,ied,jsd,jed
integer :: iscnx,iecnx,jscnx,jecnx
integer :: isdnx,iednx,jsdnx,jednx
integer :: nxc,nyc,nxch,nych

!--- get params from fv domain mosaic for building domain_ncellx
  call mpp_get_global_domain(fv_domain,xsize=nx,ysize=ny,position=CENTER)
  call mpp_get_layout(fv_domain,layout)
  ntiles = mpp_get_ntile_count(fv_domain)
  !write(1000+mpp_pe(),*) "nx,ny: ",nx,ny
  !write(1000+mpp_pe(),*) "layout: ",layout

!--- define mosaic for domain_ncellx refined by 'ncells' from fv_domain
  nxncells=nx*ncells+1
  nyncells=ny*ncells+1

  allocate(pe_start(ntiles))
  allocate(pe_end(ntiles))
  do n = 1, ntiles
    pe_start(n) = mpp_root_pe() + (n-1)*layout(1)*layout(2)
    pe_end(n)   = mpp_root_pe() +     n*layout(1)*layout(2)-1
  enddo
  call define_cubic_mosaic(domain_ncellx, nxncells-1, nyncells-1, layout, pe_start, pe_end, halo1 )
  deallocate(pe_start)
  deallocate(pe_end)

end subroutine define_ca_domain
!---------------------------------------------------------------------------------------------

subroutine write_ca_restart(fv_domain,scells,timestamp)
!Write restart files 

use fms_io_mod,          only: restart_file_type, free_restart_type, &
                               register_restart_field,               &
                               restore_state, save_restart

implicit none
type(domain2d),intent(inout) :: fv_domain
type(domain2d) :: domain_ncellx
integer,intent(in) :: scells
character(len=32), optional, intent(in)    :: timestamp
character(len=32)  :: fn_phy = 'ca_data.nc'

type(restart_file_type) :: CA_restart
real    :: pi,re,dx
integer :: id_restart,ncells,nx,ny
integer :: nxncells, nyncells

!Return if not allocated:
if(.not. allocated(board) .or. .not. allocated(lives))return

call mpp_get_global_domain(fv_domain,xsize=nx,ysize=ny,position=CENTER)
!Set time and length scales:                                                                                                                          
 pi=3.14159
 re=6371000.
 dx=0.5*pi*re/real(nx)
 ncells=int(dx/real(scells))
 ncells= MIN(ncells,10)

!Get CA domain
call define_ca_domain(fv_domain,domain_ncellx,ncells,nxncells, nyncells)

!Register restart field
id_restart = register_restart_field (CA_restart, fn_phy, "board", &
     board(:,:,1), domain = domain_ncellx, mandatory=.false.)

id_restart = register_restart_field (CA_restart, fn_phy, "lives", &
     lives(:,:,1), domain = domain_ncellx,  mandatory=.false.)

call save_restart(CA_restart, timestamp)

end subroutine write_ca_restart

subroutine read_ca_restart(fv_domain,scells)
!Read restart files
use fms_io_mod,          only: restart_file_type, free_restart_type, &
                               register_restart_field,               &
                               restore_state, save_restart
use mpp_mod,             only: mpp_error,  mpp_pe, mpp_root_pe, &
                               mpp_chksum, NOTE,   FATAL
use fms_mod,             only: file_exist, stdout

implicit none
type(domain2d) :: domain_ncellx
type(restart_file_type) :: CA_restart
type(domain2D), intent(inout) :: fv_domain
integer,intent(in) :: scells
character(len=32)  :: fn_phy = 'ca_data.nc'

character(len=64) :: fname
integer :: id_restart
integer :: nxncells, nyncells
integer :: isdnx,iednx,jsdnx,jednx
integer :: iscnx,iecnx,jscnx,jecnx
integer :: nxc,nyc,nca
real    :: pi,re,dx
integer :: ncells,nx,ny

nca=1

call mpp_get_global_domain(fv_domain,xsize=nx,ysize=ny,position=CENTER)
!Set time and length scales:                                                                                                                          
 pi=3.14159
 re=6371000.
 dx=0.5*pi*re/real(nx)
 ncells=int(dx/real(scells))
 ncells= MIN(ncells,10)

!Get CA domain                                                                                                                                                                    
call define_ca_domain(fv_domain,domain_ncellx,ncells,nxncells,nyncells)
call mpp_get_data_domain    (domain_ncellx,isdnx,iednx,jsdnx,jednx)
call mpp_get_compute_domain (domain_ncellx,iscnx,iecnx,jscnx,jecnx)

nxc = iecnx-iscnx+1
nyc = jecnx-jscnx+1

if (.not. allocated(board))then
   allocate(board(nxc,nyc,nca))
endif
if (.not. allocated(lives))then
   allocate(lives(nxc,nyc,nca))
endif

!Read restart
id_restart = register_restart_field (CA_restart, fn_phy, "board", &
     board(:,:,1), domain = domain_ncellx, mandatory=.false.)

id_restart = register_restart_field (CA_restart, fn_phy, "lives", &
     lives(:,:,1), domain = domain_ncellx,  mandatory=.false.)

fname = 'INPUT/'//trim(fn_phy)
if (file_exist(fname)) then
   !--- read the CA restart data
   call mpp_error(NOTE,'reading CA restart data from INPUT/ca_data.tile*.nc')
   call restore_state(CA_restart)
else
   call mpp_error(NOTE,'No CA restarts - cold starting CA')
   return
endif

end subroutine read_ca_restart

subroutine update_cells_sgs(kstep,initialize_ca,first_flag,restart,first_time_step,iseed_ca,nca,nxc,nyc,nxch,nych,nlon,&
                            nlat,nxncells,nyncells,isc,iec,jsc,jec, npx,npy,isdnx,iednx,jsdnx,jednx,   &
                            iscnx,iecnx,jscnx,jecnx,domain_ncellx,CA,ca_plumes,iini,ilives_in,nlives,     &
                            nfracseed,nseed,nspinup,nf,nca_plumes,ncells)

implicit none

integer, intent(in)  :: kstep,nxc,nyc,nlon,nlat,nxch,nych,nca,isc,iec,jsc,jec,npx,npy
integer, intent(in)  :: iini(nxc,nyc,nca),iseed_ca,initialize_ca,ilives_in(nxc,nyc,nca)
integer, intent(in)  :: nxncells,nyncells,isdnx,iednx,jsdnx,jednx,iscnx,iecnx,jscnx,jecnx
real,    intent(out) :: CA(nlon,nlat)
integer, intent(out) :: ca_plumes(nlon,nlat)
integer, intent(in)  :: nlives,nseed, nspinup, nf,ncells
real,    intent(in)  :: nfracseed
logical, intent(in)  :: nca_plumes,restart,first_flag,first_time_step
type(domain2D), intent(inout) :: domain_ncellx
real,    dimension(nlon,nlat) :: frac
integer, allocatable  :: V(:),L(:),B(:)
integer, allocatable  :: AG(:,:)
integer              :: inci, incj, i, j, k,sub,spinup,it,halo,k_in,isize,jsize
integer              :: ih, jh,kend, count4,boardmax,livemax,ilivemax
real,    allocatable :: board_halo(:,:,:)
integer, dimension(nxc,nyc) :: neighbours, birth, newlives,thresh
integer, dimension(nxc,nyc) :: neg, newcell, oldlives, newval,temp,newseed
integer, dimension(ncells,ncells) :: onegrid
integer(8)           :: count, count_rate, count_max, count_trunc
integer(8)           :: iscale = 10000000000
logical, save        :: start_from_restart

real, dimension(nxc,nyc) :: NOISE_B
real, dimension(nxc*nyc) :: noise1D2


!------------------------------------------------------------------------------------------------

if(first_time_step)then
start_from_restart = .False.
endif

!-------------------------------------------------------------------------------------------------
halo=1
isize=nlon+2*halo
jsize=nlat+2*halo
k_in=1
 
  if (.not. allocated(board))then
     allocate(board(nxc,nyc,nca))
  endif
  if (.not. allocated(lives))then
     allocate(lives(nxc,nyc,nca))
  endif
  if(.not. allocated(board_halo))then
     allocate(board_halo(nxch,nych,1))
  endif
 
 !Step 1: Generate a new random number each time-step for "random seeding"
 !each nseed timestep where mod(kstep,nseed) = 0

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
    count4 = mod(kstep + mype + iseed_ca + 2147483648, 4294967296) - 2147483648
 endif
 call random_setseed(count4)
 
 noise1D2 = 0.0
 call random_number(noise1D2)

 !Put on 2D:
 do j=1,nyc
  do i=1,nxc
     NOISE_B(i,j)=noise1D2(i+(j-1)*nxc)    
  enddo
 enddo

 !Step 2: Initialize CA, if restart data exist (board,lives > 0) initialize from restart file, otherwise initialize at time-
 !step initialize_ca.
 boardmax=maxval(board)
 call mp_reduce_max(boardmax)
 livemax=maxval(lives)
 call mp_reduce_max(livemax)
 ilivemax=maxval(ilives_in)
 call mp_reduce_max(ilivemax)


 if(restart .and. first_time_step .and. boardmax > 0 .and. livemax > 0)then
    !restart
    start_from_restart = .true.
    spinup = 1
 else

   if(kstep < initialize_ca .and. .not. start_from_restart)then
    do j=1,nyc
     do i=1,nxc
      board(i,j,nf) = 0
      lives(i,j,nf) = 0
     enddo
    enddo
   endif

  if(kstep == initialize_ca .and. .not. start_from_restart)then 
   do j=1,nyc
    do i=1,nxc
    board(i,j,nf) = iini(i,j,nf)
    lives(i,j,nf) = ilives_in(i,j,nf)*iini(i,j,nf)
   enddo
   enddo
   spinup=nspinup
  else
   spinup=1
  endif

 endif

  newseed = 0

 !seed with new active cells each nseed time-step regardless of restart/cold start

  if(mod(kstep,nseed)==0. .and. (kstep >= initialize_ca .or. start_from_restart))then
   do j=1,nyc
    do i=1,nxc
     if(board(i,j,nf) == 0 .and. NOISE_B(i,j)>0.90 )then
       newseed(i,j) = 1
     endif
     board(i,j,nf) = board(i,j,nf) + newseed(i,j)
    enddo
   enddo
  endif

 
 !Step 3: Evolve CA
 do it = 1,spinup
 
 CA=0
 neighbours=0
 birth=0
 newlives=0
 neg=0
 newcell=0
 oldlives=0
 newval=0
 frac=0
 board_halo=0

 
 !--- copy board into the halo-augmented board_halo                                                         
 board_halo(1+halo:nxc+halo,1+halo:nyc+halo,1) = real(board(1:nxc,1:nyc,1),kind=8)
 !write(1000+mpp_pe(),*) "board_halo pre: ",board_halo(:,:,1)

 !--- perform halo update
 call atmosphere_scalar_field_halo (board_halo, halo, nxch, nych, 1, &
                                     iscnx, iecnx, jscnx, jecnx, &
                                     nxncells, nyncells, domain_ncellx)

 !--- output data to ensure proper update                                                                            
 !write(1000+mpp_pe(),*) "board_halo post: ",board_halo(:,:,1)

 !--- Count the neighbours
  do j=1,nyc
     do i=1,nxc
        ih=i+halo
        jh=j+halo
        neighbours(i,j)=board_halo(ih-1,jh-1,1)+board_halo(ih-1,jh,1)+ &
                        board_halo(ih-1,jh+1,1)+board_halo(ih,jh+1,1)+board_halo(ih+1,jh+1,1)+&
                        board_halo(ih+1,jh,1)+board_halo(ih+1,jh-1,1)+board_halo(ih,jh-1,1)
     enddo
  enddo

 !--- Check rules; 

 !birth
  do j=1,nyc
   do i=1,nxc
     if((neighbours(i,j) == 3 .or. neighbours(i,j) == 2))then 
     birth(i,j)=1
     endif
   enddo
  enddo
 
 !death                                                                                                                        
  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j) < 2 .or. neighbours(i,j) > 3)then  
     lives(i,j,nf)=lives(i,j,nf) - 1
     endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
   if(lives(i,j,nf) < 0)then
     lives(i,j,nf)=0
   endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
    if(birth(i,j)==1 .and. lives(i,j,nf)==0)then
    newcell(i,j)=1
    endif
   enddo
  enddo


  do j=1,nyc
     do i=1,nxc
        lives(i,j,nf)=lives(i,j,nf)+newcell(i,j)*ilives_in(i,j,nf)
     enddo
  enddo

  do j=1,nyc
   do i=1,nxc
    if(neighbours(i,j)==3 .or. (board(i,j,nf)==1 .and. neighbours(i,j)==2))then
    board(i,j,nf)=1
    else
    board(i,j,nf)=0
    endif
   enddo
  enddo


 enddo !spinup


!COARSE-GRAIN BACK TO NWP GRID
 
  inci=ncells
  incj=ncells
  sub=ncells-1
  DO j=1,nlat
     DO i=1,nlon
        CA(i,j)=(SUM(lives(inci-sub:inci,incj-sub:incj,nf)))/real(ncells*ncells)
        inci=inci+ncells
     ENDDO
     inci=ncells
     incj=incj+ncells
  ENDDO

if(nca_plumes) then
!COMPUTE NUMBER OF CLUSTERS (CONVECTIVE PLUMES) IN EACH CA-CELL
!Note, at the moment we only use the count of the plumes found in a grid-cell
!In the future the routine "plumes" can also be used to give the size of 
!each individual plume for better coupling to the convection scheme.

  temp=0
  do j=1,nyc
   do i=1,nxc
     if(lives(i,j,1) > 0)then
      temp(i,j)=1  
     endif
   enddo
  enddo

  kend=ceiling((ncells*ncells)/2.)
  if (.not. allocated(V))then
  allocate(V(kend))
  endif
  if (.not. allocated(L))then
  allocate(L(kend))
  endif
  if (.not. allocated(B))then
  allocate(B(kend))
  endif
  if (.not. allocated(AG))then
  allocate(AG(ncells,ncells))
  endif
  
  ca_plumes(:,:)=0
  inci=ncells
  incj=ncells
  sub=ncells-1
  DO j=1,nlat
     DO i=1,nlon
        B(:)=0
        L(:)=0
        V(:)=0
        onegrid(1:ncells,1:ncells)=temp(inci-sub:inci,incj-sub:incj)
        call plumes(V,L,AG,onegrid,ncells,ncells,kend)
        do k=1,kend
           if(V(k)==1)then
              B(k)=L(k) !to avoid considering clusters of 0
           endif
        enddo
        ca_plumes(i,j)=MAXVAL(B(1:kend))
        inci=inci+ncells
     ENDDO
     inci=ncells
     incj=incj+ncells
  ENDDO

else

ca_plumes(:,:)=0.

endif ! nca_plumes

end subroutine update_cells_sgs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_cells_global(kstep,first_time_step,iseed_ca,nca,nxc,nyc,nxch,nych,nlon,nlat,nxncells,nyncells,isc,iec,jsc,jec, &
                        npx,npy,iscnx,iecnx,jscnx,jecnx,domain_ncellx,CA,iini_g,ilives_g,                &
                        nlives,ncells,nfracseed,nseed,nspinup,nf)

implicit none

integer, intent(in) :: kstep,nxc,nyc,nlon,nlat,nxch,nych,nca,isc,iec,jsc,jec,npx,npy
integer, intent(in) :: iini_g(nxc,nyc,nca), ilives_g(nxc,nyc), iseed_ca,nxncells,nyncells
real, intent(out) :: CA(nlon,nlat)
logical, intent(in) :: first_time_step
integer, intent(in) :: nlives, ncells, nseed, nspinup, nf,iscnx,iecnx,jscnx,jecnx
real, intent(in) :: nfracseed
type(domain2D), intent(inout) :: domain_ncellx
real, dimension(nlon,nlat) :: frac
integer,allocatable,save :: board_g(:,:,:), lives_g(:,:,:)
integer,allocatable :: V(:),L(:)
integer :: inci, incj, i, j, k ,sub,spinup,it,halo,k_in,isize,jsize
integer :: ih, jh, count4
real, allocatable :: board_halo(:,:,:)
integer, dimension(nxc,nyc) :: neighbours, birth, newlives, thresh
integer, dimension(nxc,nyc) :: neg, newcell, oldlives, newval,temp,newseed
real, dimension(nxc,nyc) :: NOISE_B
real, dimension(nxc*nyc) :: noise1D2
integer(8) :: count, count_rate, count_max, count_trunc
integer(8) :: iscale = 10000000000

!-------------------------------------------------------------------------------------------------
halo=1
isize=nlon+2*halo
jsize=nlat+2*halo
k_in=1

 if (.not. allocated(board_g))then
 allocate(board_g(nxc,nyc,nca))
 endif
 if (.not. allocated(lives_g))then
 allocate(lives_g(nxc,nyc,nca))
 endif
 if(.not. allocated(board_halo))then                                                                      
 allocate(board_halo(nxch,nych,1))   
 endif

 !Generate a new random number each time-step for "random seeding"
 !each nseed timestep
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

 noise1D2 = 0.0
 call random_number(noise1D2)
 
 !random numbers:
 noise1D2 = 0.0

 call random_number(noise1D2)

 !Put on 2D:
 do j=1,nyc
  do i=1,nxc
  NOISE_B(i,j)=noise1D2(i+(j-1)*nxc)
  enddo
 enddo                        


  if(first_time_step)then
   do j=1,nyc
    do i=1,nxc
     board_g(i,j,nf) = iini_g(i,j,nf)
     lives_g(i,j,nf) = ilives_g(i,j)*iini_g(i,j,nf)
    enddo
   enddo

  endif

 !Seed with new CA cells at each nseed step
  newseed=0

  if(mod(kstep,nseed) == 0)then
   do j=1,nyc
    do i=1,nxc
     if(board_g(i,j,nf) == 0 .and. NOISE_B(i,j)>0.75 )then
        newseed(i,j)=1
     endif
     board_g(i,j,nf) = board_g(i,j,nf) + newseed(i,j)
    enddo
   enddo
  endif

  if(first_time_step)then
  spinup=nspinup
  else
  spinup = 1
  endif
 

do it=1,spinup
!Step 2 - Initialize variables to 0 and extract the halo
 
 neighbours=0
 birth=0
 newlives=0
 neg=0
 newcell=0
 oldlives=0
 newval=0
 frac=0
 board_halo=0

!The input to scalar_field_halo needs to be 1D.
!take the updated board_g fields and extract the halo
! in order to have updated values in the halo region. 

 !--- copy board into the halo-augmented board_halo                                                                                                    
 board_halo(1+halo:nxc+halo,1+halo:nyc+halo,1) = real(board(1:nxc,1:nyc,1),kind=8)
 !write(1000+mpp_pe(),*) "board_halo pre: ",board_halo(:,:,1)                                                                                          

 !--- perform halo update                                                                                                                              
 call atmosphere_scalar_field_halo (board_halo, halo, nxch, nych, 1, &
                                     iscnx, iecnx, jscnx, jecnx, &
                                     nxncells, nyncells, domain_ncellx)

  do j=1,nyc
     do i=1,nxc
        ih=i+halo
        jh=j+halo
        neighbours(i,j)=board_halo(ih-1,jh-1,1)+board_halo(ih-1,jh,1)+ &
                        board_halo(ih-1,jh+1,1)+board_halo(ih,jh+1,1)+board_halo(ih+1,jh+1,1)+&
                        board_halo(ih+1,jh,1)+board_halo(ih+1,jh-1,1)+board_halo(ih,jh-1,1)
     enddo
  enddo



  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j)==2 .or. neighbours(i,j)==3)then
     birth(i,j)=1
     endif
   enddo
  enddo


  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j)<2 .or. neighbours(i,j)>3)then
     lives_g(i,j,nf)=lives_g(i,j,nf) - 1
     endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
   if(lives_g(i,j,nf)<0)then
     lives_g(i,j,nf)=0
   endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
    if(birth(i,j)==1 .and. lives_g(i,j,nf)==0)then
    newcell(i,j)=1
    endif
   enddo
  enddo


  do j=1,nyc
   do i=1,nxc
    lives_g(i,j,nf)=lives_g(i,j,nf)+newcell(i,j)*ilives_g(i,j)
   enddo
  enddo


   do j=1,nyc
   do i=1,nxc
    if( (board_g(i,j,nf) ==1 .and. (neighbours(i,j)==3 .or. neighbours(i,j)==2) ).or. (board_g(i,j,nf)==0 .and. neighbours(i,j)==3) )then
    board_g(i,j,nf)=1
    else
    board_g(i,j,nf)=0
    endif
   enddo
  enddo

enddo !spinup

!COARSE-GRAIN BACK TO NWP GRID
 
  inci=ncells
  incj=ncells
  sub=ncells-1
  DO j=1,nlat
     DO i=1,nlon
        CA(i,j)=(SUM(lives_g(inci-sub:inci,incj-sub:incj,nf)))/(ncells*ncells)
        inci=inci+ncells
     ENDDO
     inci=ncells
     incj=incj+ncells
  ENDDO

end subroutine update_cells_global

!================================
 ! This subroutine is copied from FMS/test_fms/test_mpp_domains.F90
  ! and modified to make it simpler to use.
  ! domain_decomp in fv_mp_mod.F90 does something similar, but it does a
  ! few other unnecessary things (and requires more arguments).
  subroutine define_cubic_mosaic(domain, ni, nj, layout, pe_start, pe_end, halo)
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: ni, nj
    integer,        intent(in)    :: layout(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    integer,        intent(in)    :: halo
    !--- local variables
    integer                       :: global_indices(4,6), layout2D(2,6)
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact
    integer                       :: i

    ntiles = 6
    num_contact = 12
    if(size(pe_start(:)) .NE. 6 .OR. size(pe_end(:)) .NE. 6 ) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of pe_start and pe_end should be 6")
    if(size(layout) .NE. 2) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of layout should be 2")
    do i = 1, 6
       layout2D(:,i) = layout(:)
       global_indices(1,i) = 1
       global_indices(2,i) = ni
       global_indices(3,i) = 1
       global_indices(4,i) = nj
    enddo
!--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1;     tile2(1) = 2
    istart1(1) = ni;  iend1(1) = ni;  jstart1(1) = 1;      jend1(1) = nj
    istart2(1) = 1;   iend2(1) = 1;   jstart2(1) = 1;      jend2(1) = nj

    !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;      iend1(2) = ni;  jstart1(2) = nj;  jend1(2) = nj
    istart2(2) = 1;      iend2(2) = 1;   jstart2(2) = nj;  jend2(2) = 1

    !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1;     tile2(3) = 5
    istart1(3) = 1;   iend1(3) = 1;      jstart1(3) = 1;   jend1(3) = nj
    istart2(3) = ni;  iend2(3) = 1;      jstart2(3) = nj;  jend2(3) = nj

    !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1; tile2(4) = 6
    istart1(4) = 1;      iend1(4) = ni;  jstart1(4) = 1;   jend1(4) = 1
    istart2(4) = 1;      iend2(4) = ni;  jstart2(4) = nj;  jend2(4) = nj

    !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2;        tile2(5) = 3
    istart1(5) = 1;      iend1(5) = ni;  jstart1(5) = nj;  jend1(5) = nj
    istart2(5) = 1;      iend2(5) = ni;  jstart2(5) = 1;   jend2(5) = 1

    !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = ni;  iend1(6) = ni;  jstart1(6) = 1;      jend1(6) = nj
    istart2(6) = ni;  iend2(6) = 1;   jstart2(6) = 1;      jend2(6) = 1

    !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2; tile2(7) = 6
    istart1(7) = 1;   iend1(7) = ni;  jstart1(7) = 1;   jend1(7) = 1
    istart2(7) = ni;  iend2(7) = ni;  jstart2(7) = nj;  jend2(7) = 1

    !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = ni;  iend1(8) = ni;  jstart1(8) = 1;      jend1(8) = nj
    istart2(8) = 1;   iend2(8) = 1;   jstart2(8) = 1;      jend2(8) = nj

    !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3; tile2(9) = 5
    istart1(9) = 1;      iend1(9) = ni;  jstart1(9) = nj;  jend1(9) = nj
    istart2(9) = 1;      iend2(9) = 1;   jstart2(9) = nj;  jend2(9) = 1

    !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4; tile2(10) = 5
    istart1(10) = 1;     iend1(10) = ni; jstart1(10) = nj; jend1(10) = nj
    istart2(10) = 1;     iend2(10) = ni; jstart2(10) = 1;  jend2(10) = 1

    !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4; tile2(11) = 6
    istart1(11) = ni; iend1(11) = ni; jstart1(11) = 1;     jend1(11) = nj
    istart2(11) = ni; iend2(11) = 1;  jstart2(11) = 1;     jend2(11) = 1

    !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5; tile2(12) = 6
    istart1(12) = ni; iend1(12) = ni; jstart1(12) = 1;     jend1(12) = nj
    istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;     jend2(12) = nj

    call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, &
         num_contact, tile1, tile2, istart1, iend1, jstart1, jend1, &
         istart2, iend2, jstart2, jend2, pe_start, pe_end, symmetry=.true., &
         whalo=halo, ehalo=halo, shalo=halo, nhalo=halo, &
         name='CA cubic mosaic')

  end subroutine define_cubic_mosaic

end module update_ca
