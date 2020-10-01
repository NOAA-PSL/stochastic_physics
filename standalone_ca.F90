program  standalone_stochy_new

use standalone_stochy_module

use atmosphere_stub_mod, only: Atm,atmosphere_init_stub
!use mpp_domains
use mpp_mod,             only: mpp_set_current_pelist,mpp_get_current_pelist,mpp_init,mpp_pe,mpp_npes ,mpp_declare_pelist
use mpp_domains_mod,     only: mpp_broadcast_domain,MPP_DOMAIN_TIME,mpp_domains_init ,mpp_domains_set_stack_size
use fms_mod,             only:  fms_init
!use time_manager_mod,    only: time_type
use xgrid_mod,           only: grid_box_type
use netcdf


implicit none
type(GFS_control_type)  :: Model
integer, parameter      :: nlevs=64
integer                 :: ntasks,fid,ct
integer                 :: nthreads,omp_get_num_threads
integer                 :: ncid,xt_dim_id,yt_dim_id,time_dim_id,xt_var_id,yt_var_id,time_var_id,ca_out_id
integer                 :: ca1_id,ca2_id,ca3_id
character*1             :: strid
type(GFS_grid_type),allocatable     :: Grid(:)
type(GFS_diag_type),allocatable :: Diag(:)
type(GFS_statein_type),allocatable :: Statein(:)
type(GFS_coupling_type),allocatable :: Coupling(:)
include 'mpif.h'
include 'netcdf.inc'
real(kind=4) :: ts,undef

integer     :: cres,blksz,nblks,ierr,my_id,i,j,nx2,ny2,nx,ny,id
integer,target :: npx,npy
integer     :: ng,layout(2),io_layout(2),commID,grid_type,ntiles
integer :: halo_update_type = 1
logical,target :: nested
integer  :: pe,npes,stackmax=4000000

real(kind=4),allocatable,dimension(:,:) :: workg
real(kind=4),allocatable,dimension(:) :: grid_xt,grid_yt
real(kind=8),pointer    ,dimension(:,:) :: area
type(grid_box_type)           :: grid_box
!type(time_type)               :: Time               ! current time
!type(time_type)               :: Time_step          ! atmospheric time step.
!type(time_type)               :: Time_init          ! reference time.
!---cellular automata control parameters
integer              :: nca             !< number of independent cellular automata
integer              :: nlives          !< cellular automata lifetime
integer              :: ncells          !< cellular automata finer grid
real                 :: nfracseed       !< cellular automata seed probability
integer              :: nseed           !< cellular automata seed frequency
logical              :: do_ca           !< cellular automata main switch
logical              :: ca_sgs          !< switch for sgs ca
logical              :: ca_global       !< switch for global ca
logical              :: ca_smooth       !< switch for gaussian spatial filter
logical              :: isppt_deep      !< switch for combination with isppt_deep. OBS! Switches off SPPT on other tendencies!
logical              :: isppt_pbl
logical              :: isppt_shal
logical              :: pert_flux   
logical              :: pert_trigger
integer              :: iseed_ca        !< seed for random number generation in ca scheme
integer              :: nspinup         !< number of iterations to spin up the ca
integer              :: ca_amplitude
real                 :: nthresh         !< threshold used for perturbed vertical velocity

NAMELIST /gfs_physics_nml/ nca, ncells, nlives, nfracseed,nseed, nthresh, &
         do_ca,ca_sgs, ca_global,iseed_ca,ca_smooth,isppt_pbl,isppt_shal,isppt_deep,nspinup,&
         pert_trigger,pert_flux,ca_amplitude

! default values
nca            = 1
ncells         = 5
nlives         = 10
nfracseed      = 0.5
nseed          = 100000
iseed_ca       = 0
nspinup        = 1
do_ca          = .false.
ca_sgs         = .false.
ca_global      = .false.
ca_smooth      = .false.
isppt_deep     = .false.
isppt_shal     = .false.
isppt_pbl      = .false.
pert_trigger   = .false.
pert_flux      = .false.
ca_amplitude   = 500.
nthresh        = 0.0

! open namelist file
open (unit=565, file='input.nml', READONLY, status='OLD', iostat=ierr)
read(565,gfs_physics_nml)
close(565)
! define stuff
ng=3  ! ghost region
undef=9.99e+20

! initialize fms
call fms_init()
call mpp_init()
call fms_init
my_id=mpp_pe()

call atmosphere_init_stub (grid_box, area)
!define domain
isd=Atm(1)%bd%isd
ied=Atm(1)%bd%ied
jsd=Atm(1)%bd%jsd
jed=Atm(1)%bd%jed
isc=Atm(1)%bd%isc
iec=Atm(1)%bd%iec
jsc=Atm(1)%bd%jsc
jec=Atm(1)%bd%jec
nx=Atm(1)%npx-1
ny=Atm(1)%npy-1
allocate(workg(nx,ny))

! for this simple test, nblocks = ny, blksz=ny
nblks=ny
blksz=nx
nthreads = omp_get_num_threads()
Model%me=my_id

Model%nca            = nca
Model%ncells         = ncells
Model%nlives         = nlives
Model%nfracseed      = nfracseed
Model%nseed          = nseed  
Model%iseed_ca       = iseed_ca
Model%nspinup        = nspinup
Model%do_ca          = do_ca
Model%ca_sgs         = ca_sgs  
Model%ca_global      = ca_global
Model%ca_smooth      = ca_smooth
Model%isppt_deep     = isppt_deep
Model%isppt_pbl      = isppt_pbl 
Model%isppt_shal     = isppt_shal
Model%pert_flux      = pert_flux 
Model%pert_trigger   = pert_trigger
Model%nthresh        = nthresh

! setup GFS_init parameters

!define model grid

allocate(grid_xt(nx),grid_yt(ny))
do i=1,nx
  grid_xt(i)=i
enddo
do i=1,ny
  grid_yt(i)=i
enddo

!setup GFS_coupling
allocate(Diag(nblks))
allocate(Coupling(nblks))
allocate(Statein(nblks))
write(strid,'(I1.1)') my_id+1
fid=30+my_id
ierr=nf90_create('ca_out.tile'//strid//'.nc',cmode=NF90_CLOBBER,ncid=ncid)
ierr=NF90_DEF_DIM(ncid,"grid_xt",nx,xt_dim_id)
ierr=NF90_DEF_DIM(ncid,"grid_yt",ny,yt_dim_id)
ierr=NF90_DEF_DIM(ncid,"time",NF90_UNLIMITED,time_dim_id)
  !> - Define the dimension variables.
ierr=NF90_DEF_VAR(ncid,"grid_xt",NF90_FLOAT,(/ xt_dim_id /), xt_var_id)
ierr=NF90_PUT_ATT(ncid,xt_var_id,"long_name","T-cell longitude")
ierr=NF90_PUT_ATT(ncid,xt_var_id,"cartesian_axis","X")
ierr=NF90_PUT_ATT(ncid,xt_var_id,"units","degrees_E")
ierr=NF90_DEF_VAR(ncid,"grid_yt",NF90_FLOAT,(/ yt_dim_id /), yt_var_id)
ierr=NF90_PUT_ATT(ncid,yt_var_id,"long_name","T-cell latitude")
ierr=NF90_PUT_ATT(ncid,yt_var_id,"cartesian_axis","Y")
ierr=NF90_PUT_ATT(ncid,yt_var_id,"units","degrees_N")
ierr=NF90_DEF_VAR(ncid,"time",NF90_FLOAT,(/ time_dim_id /), time_var_id)
ierr=NF90_PUT_ATT(ncid,time_var_id,"long_name","time")
ierr=NF90_PUT_ATT(ncid,time_var_id,"units","hours since 2014-08-01 00:00:00")
ierr=NF90_PUT_ATT(ncid,time_var_id,"cartesian_axis","T")
ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar_type","JULIAN")
ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar","JULIAN")
!ierr=NF90_DEF_VAR(ncid,"ca_out",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), ca_out_id)
!ierr=NF90_PUT_ATT(ncid,ca_out_id,"long_name","random pattern")
!ierr=NF90_PUT_ATT(ncid,ca_out_id,"units","None")
!ierr=NF90_PUT_ATT(ncid,ca_out_id,"missing_value",undef)
!ierr=NF90_PUT_ATT(ncid,ca_out_id,"_FillValue",undef)
!ierr=NF90_PUT_ATT(ncid,ca_out_id,"cell_methods","time: point")
ierr=NF90_DEF_VAR(ncid,"ca1",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), ca1_id)
ierr=NF90_PUT_ATT(ncid,ca1_id,"long_name","random pattern")
ierr=NF90_PUT_ATT(ncid,ca1_id,"units","None")
ierr=NF90_PUT_ATT(ncid,ca1_id,"missing_value",undef)
ierr=NF90_PUT_ATT(ncid,ca1_id,"_FillValue",undef)
ierr=NF90_PUT_ATT(ncid,ca1_id,"cell_methods","time: point")
ierr=NF90_DEF_VAR(ncid,"ca2",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), ca2_id)
ierr=NF90_PUT_ATT(ncid,ca2_id,"long_name","random pattern")
ierr=NF90_PUT_ATT(ncid,ca2_id,"units","None")
ierr=NF90_PUT_ATT(ncid,ca2_id,"missing_value",undef)
ierr=NF90_PUT_ATT(ncid,ca2_id,"_FillValue",undef)
ierr=NF90_PUT_ATT(ncid,ca2_id,"cell_methods","time: point")
ierr=NF90_DEF_VAR(ncid,"ca3",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), ca3_id)
ierr=NF90_PUT_ATT(ncid,ca3_id,"long_name","random pattern")
ierr=NF90_PUT_ATT(ncid,ca3_id,"units","None")
ierr=NF90_PUT_ATT(ncid,ca3_id,"missing_value",undef)
ierr=NF90_PUT_ATT(ncid,ca3_id,"_FillValue",undef)
ierr=NF90_PUT_ATT(ncid,ca3_id,"cell_methods","time: point")
ierr=NF90_ENDDEF(ncid)
ierr=NF90_PUT_VAR(ncid,xt_var_id,grid_xt)
ierr=NF90_PUT_VAR(ncid,yt_var_id,grid_yt)
! allocate diagnostics
DO i =1,nblks
   allocate(Diag(i)%ca_out(blksz))
   allocate(Diag(i)%ca_deep(blksz))
   allocate(Diag(i)%ca_turb(blksz))
   allocate(Diag(i)%ca_shal(blksz))
   allocate(Diag(i)%ca_rad(blksz))
   allocate(Diag(i)%ca_micro(blksz))
   allocate(Diag(i)%ca1(blksz))
   allocate(Diag(i)%ca2(blksz))
   allocate(Diag(i)%ca3(blksz))
! allocate coupling
   allocate(Coupling(i)%cape(blksz))
   allocate(Coupling(i)%ca_out(blksz))
   allocate(Coupling(i)%ca_deep(blksz))
   allocate(Coupling(i)%ca_turb(blksz))
   allocate(Coupling(i)%ca_shal(blksz))
   allocate(Coupling(i)%ca_rad(blksz))
   allocate(Coupling(i)%ca_micro(blksz))
   allocate(Coupling(i)%ca1(blksz))
   allocate(Coupling(i)%ca2(blksz))
   allocate(Coupling(i)%ca3(blksz))
! allocate coupling
   allocate(Statein(i)%pgr(blksz))
   allocate(Statein(i)%qgrs(blksz,nlevs,1))
   allocate(Statein(i)%vvl(blksz,nlevs))
   allocate(Statein(i)%prsl(blksz,nlevs))
ENDDO
ct=1
do i=1,600
   ts=i/8.0  ! hard coded to write out hourly based on a 450 second time-step
   call cellular_automata_global(i-1, Statein, Coupling, Diag, &
                          nblks, Model%levs, Model%nca, Model%ncells,          &
                          Model%nlives, Model%nfracseed, Model%nseed,                    &
                          Model%nthresh, Model%ca_global, Model%ca_sgs,                  &
                          Model%iseed_ca, Model%ca_smooth, Model%nspinup,                &
                          blksz)
   if (mod(i,8).EQ.0) then
      do j=1,ny
         workg(:,j)=Diag(j)%ca1(:)   
      enddo
      ierr=NF90_PUT_VAR(ncid,ca1_id,workg,(/1,1,ct/))
      do j=1,ny
         workg(:,j)=Diag(j)%ca2(:)   
      enddo
      ierr=NF90_PUT_VAR(ncid,ca2_id,workg,(/1,1,ct/))
      do j=1,ny
         workg(:,j)=Diag(j)%ca3(:)   
      enddo
      ierr=NF90_PUT_VAR(ncid,ca3_id,workg,(/1,1,ct/))
      ierr=NF90_PUT_VAR(ncid,time_var_id,ts,(/ct/))
      ct=ct+1
   if (my_id.EQ.0) write(6,fmt='(a,i5,4f6.3)') 'ca=',i,Diag(1)%ca1(1:4)
   endif
enddo
!close(fid)
ierr=NF90_CLOSE(ncid)
end
