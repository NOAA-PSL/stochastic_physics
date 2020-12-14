program  standalone_stochy

use standalone_stochy_module
use stochastic_physics,  only : init_stochastic_physics,run_stochastic_physics

use atmosphere_stub_mod, only: Atm,atmosphere_init_stub
!use mpp_domains
use mpp_mod,             only: mpp_set_current_pelist,mpp_get_current_pelist,mpp_init,mpp_pe,mpp_npes ,mpp_declare_pelist
use mpp_domains_mod,     only: mpp_broadcast_domain,MPP_DOMAIN_TIME,mpp_domains_init ,mpp_domains_set_stack_size
use fms_mod,             only:  fms_init
use xgrid_mod,           only: grid_box_type
use netcdf

implicit none
type(GFS_control_type)  :: Model
type(GFS_init_type)     :: Init_parm
integer, parameter      :: nlevs=64
integer                 :: ntasks,fid
integer                 :: nthreads,omp_get_num_threads
integer                 :: ncid,xt_dim_id,yt_dim_id,time_dim_id,xt_var_id,yt_var_id,time_var_id,var_id_lat,var_id_lon,var_id_tile
integer                 :: varid1,varid2,varid3,varid4,varid_lon,varid_lat,varid_tile
integer                 :: zt_dim_id,zt_var_id
character*1             :: strid
type(GFS_grid_type),allocatable     :: Grid(:)
type(GFS_coupling_type),allocatable :: Coupling(:)
! stochastic namelist fields
integer nssppt,nsshum,nsskeb,lon_s,lat_s,ntrunc
integer skeb_varspect_opt,skeb_npass
logical sppt_sfclimit

real(kind=kind_dbl_prec) :: skeb_sigtop1,skeb_sigtop2,          &
                   sppt_sigtop1,sppt_sigtop2,shum_sigefold, &
                   skeb_vdof
real(kind=kind_dbl_prec) skeb_diss_smooth,spptint,shumint,skebint,skebnorm
real(kind=kind_dbl_prec), dimension(5) :: skeb,skeb_lscale,skeb_tau
real(kind=kind_dbl_prec), dimension(5) :: sppt,sppt_lscale,sppt_tau
real(kind=kind_dbl_prec), dimension(5) :: shum,shum_lscale,shum_tau
integer,dimension(5) ::skeb_vfilt
integer(8),dimension(5) ::iseed_sppt,iseed_shum,iseed_skeb
logical stochini,sppt_logit,new_lscale
logical use_zmtnblck
include 'mpif.h'
include 'netcdf.inc'
real :: ak(nlevs+1),bk(nlevs+1)
real(kind=4) :: ts,undef

data ak(:) /0.000, 0.000, 0.575, 5.741, 21.516, 55.712, 116.899, 214.015, 356.223, 552.720, 812.489, &
   1143.988, 1554.789, 2051.150, 2637.553, 3316.217, 4086.614, 4945.029, 5884.206, 6893.117,    &
   7956.908, 9057.051, 10171.712, 11276.348, 12344.490, 13348.671, 14261.435, 15056.342,        &
   15708.893, 16197.315, 16503.145, 16611.604, 16511.736, 16197.967, 15683.489, 14993.074,      &
   14154.316, 13197.065, 12152.937, 11054.853, 9936.614, 8832.537, 7777.150, 6804.874, 5937.050,&
   5167.146, 4485.493, 3883.052, 3351.460, 2883.038, 2470.788, 2108.366, 1790.051, 1510.711,    &
   1265.752, 1051.080, 863.058, 698.457, 554.424, 428.434, 318.266, 221.958, 137.790, 64.247,0.0 /
data bk(:) /1.00000000, 0.99467117, 0.98862660, 0.98174226, 0.97386760, 0.96482760, 0.95443410, 0.94249105, &
  0.92879730, 0.91315103, 0.89535499, 0.87522358, 0.85259068, 0.82731885, 0.79930973, 0.76851469, &
  0.73494524, 0.69868290, 0.65988702, 0.61879963, 0.57574666, 0.53113484, 0.48544332, 0.43921080, &
  0.39301825, 0.34746850, 0.30316412, 0.26068544, 0.22057019, 0.18329623, 0.14926878, 0.11881219, &
  0.09216691, 0.06947458, 0.05064684, 0.03544162, 0.02355588, 0.01463712, 0.00829402, 0.00410671, &
  0.00163591, 0.00043106, 0.00003697, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, &
  0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, &
  0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, &
  0.00000000 /
integer     :: cres,blksz,nblks,ierr,my_id,i,j,k,nx2,ny2,nx,ny,id
integer,target :: npx,npy
integer     :: ng,layout(2),io_layout(2),commID,grid_type,ntiles
integer :: halo_update_type = 1
real        :: dx,dy,pi,rd,cp
logical,target :: nested
logical   :: write_this_tile
integer  :: nargs,ntile_out,nlunit,pe,npes,stackmax=4000000
character*80 :: fname
character*1  :: ntile_out_str
integer :: iret

real(kind=4),allocatable,dimension(:,:) :: workg,tile_number
real(kind=4),allocatable,dimension(:,:,:) :: workg3d
real(kind=4),allocatable,dimension(:) :: grid_xt,grid_yt
real(kind=8),pointer    ,dimension(:,:) :: area
real(kind=8)                    :: ex3d(nlevs+1),pressi(nlevs+1),pressl(nlevs),p1000,exn

type(grid_box_type)           :: grid_box

      namelist /nam_stochy/ntrunc,lon_s,lat_s,sppt,sppt_tau,sppt_lscale,sppt_logit, &
      iseed_shum,iseed_sppt,shum,shum_tau,&
      shum_lscale,stochini,skeb_varspect_opt,sppt_sfclimit, &
      skeb,skeb_tau,skeb_vdof,skeb_lscale,iseed_skeb,skeb_vfilt,skeb_diss_smooth, &
      skeb_sigtop1,skeb_sigtop2,skebnorm,sppt_sigtop1,sppt_sigtop2,&
      shum_sigefold,spptint,shumint,skebint,skeb_npass,use_zmtnblck,new_lscale
write_this_tile=.false.
ntile_out_str='0'
nargs=iargc()
if (nargs.EQ.1) then
   call getarg(1,ntile_out_str)
endif
read(ntile_out_str,'(I1.1)') ntile_out
open (unit=nlunit, file='input.nml', READONLY, status='OLD')
read(nlunit,nam_stochy)
close(nlunit)
Model%do_sppt=.false.
Model%do_shum=.false.
Model%do_skeb=.false.
if (sppt(1).GT.0) Model%do_sppt=.true.
if (shum(1).GT.0) Model%do_shum=.true.
if (skeb(1).GT.0) Model%do_skeb=.true.
! define stuff
ng=3
pi=3.14159265359
undef=9.99e+20
p1000=100000.0
!define mid-layer pressure
rd=287.0
cp=1004.0
DO k=1,nlevs
   pressi(k)=ak(k)+p1000*bk(k)
ENDDO
ex3d=cp*(pressi/p1000)**(rd/cp)
DO k=1,nlevs
   exn = (ex3d(k)*pressi(k)-ex3d(k+1)*pressi(k+1))/((cp+rd)*(pressi(k)-pressi(k+1)))
   pressl(k)=p1000*exn**(cp/rd)
ENDDO

call fms_init()
call mpp_init()
call fms_init
my_id=mpp_pe()
ntasks=mpp_npes()

call atmosphere_init_stub (grid_box, area)
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
allocate(tile_number(nx,ny))
allocate(workg3d(nx,ny,nlevs))
nblks=ny
blksz=nx
Allocate(Model%blksz(nblks))
Model%blksz(:)=blksz
nthreads = omp_get_num_threads()
Model%me=my_id
Model%phour=0
Model%kdt=1
Model%dtp=900
Model%fn_nml='input.nml'
Model%levs=nlevs
allocate(Init_parm%blksz(nblks))
Init_parm%blksz(:)=blksz
! setup GFS_init parameters
allocate(Init_parm%ak(nlevs+1))
allocate(Init_parm%bk(nlevs+1))
Init_parm%ak=ak
Init_parm%bk=bk
Init_parm%nlunit=21

!define model grid
Model%nx=nx
Model%ny=ny
dx=360.0/Model%nx
dy=180.0/Model%ny
allocate(Init_parm%xlon(Model%nx,Model%ny))
allocate(Init_parm%xlat(Model%nx,Model%ny))
Init_parm%xlon(:,:)=Atm(1)%gridstruct%agrid(:,:,1)
Init_parm%xlat(:,:)=Atm(1)%gridstruct%agrid(:,:,2)

allocate(Grid(nblks))
do i=1,nblks
   allocate(Grid(i)%xlat(blksz))
   allocate(Grid(i)%xlon(blksz))
enddo
do j=1,nblks
     Grid(j)%xlat(:)=Init_parm%xlat(:,j)
     Grid(j)%xlon(:)=Init_parm%xlon(:,j)
enddo
allocate(grid_xt(nx),grid_yt(ny))
do i=1,nx
  grid_xt(i)=i
enddo
do i=1,ny
  grid_yt(i)=i
enddo
!setup GFS_coupling
allocate(Coupling(nblks))
call init_stochastic_physics(Model, Init_parm, ntasks, nthreads, iret)
if (iret .ne. 0) print *, 'ERROR init_stochastic_physics call' ! Draper - need proper error trapping here
call get_outfile(fname)
write(strid,'(I1.1)') my_id+1
if (ntile_out.EQ.0) write_this_tile=.true.
if ((my_id+1).EQ.ntile_out) write_this_tile=.true.
print*,trim(fname)//'.tile'//strid//'.nc',write_this_tile
if (write_this_tile) then
fid=30+my_id
!ierr=nf90_create(trim(fname)//'.tile'//strid//'.nc',cmode=NF90_CLOBBER,ncid=ncid)
ierr=nf90_create(trim(fname)//'.tile'//strid//'.nc',cmode=NF90_CLOBBER,ncid=ncid)
ierr=NF90_DEF_DIM(ncid,"grid_xt",nx,xt_dim_id)
ierr=NF90_DEF_DIM(ncid,"grid_yt",ny,yt_dim_id)
if (Model%do_skeb)ierr=NF90_DEF_DIM(ncid,"p_ref",nlevs,zt_dim_id)
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
ierr=NF90_DEF_VAR(ncid,"grid_lat",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_lat)
ierr=NF90_PUT_ATT(ncid,var_id_lat,"long_name","T-cell latitudes")
ierr=NF90_PUT_ATT(ncid,var_id_lat,"units","degrees_N")
ierr=NF90_PUT_ATT(ncid,var_id_lat,"missing_value",undef)
ierr=NF90_PUT_ATT(ncid,var_id_lat,"_FillValue",undef)
ierr=NF90_DEF_VAR(ncid,"grid_lon",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_lon)
ierr=NF90_PUT_ATT(ncid,var_id_lon,"long_name","T-cell longitudes")
ierr=NF90_PUT_ATT(ncid,var_id_lon,"units","degrees_N")
ierr=NF90_PUT_ATT(ncid,var_id_lon,"missing_value",undef)
ierr=NF90_PUT_ATT(ncid,var_id_lon,"_FillValue",undef)
ierr=NF90_DEF_VAR(ncid,"tile_num",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_tile)
ierr=NF90_PUT_ATT(ncid,var_id_tile,"long_name","tile number")
ierr=NF90_PUT_ATT(ncid,var_id_tile,"missing_value",undef)
ierr=NF90_PUT_ATT(ncid,var_id_tile,"_FillValue",undef)
if (Model%do_skeb)then
   ierr=NF90_DEF_VAR(ncid,"p_ref",NF90_FLOAT,(/ zt_dim_id /), zt_var_id)
   ierr=NF90_PUT_ATT(ncid,zt_var_id,"long_name","reference pressure")
   ierr=NF90_PUT_ATT(ncid,zt_var_id,"cartesian_axis","Z")
   ierr=NF90_PUT_ATT(ncid,zt_var_id,"units","Pa")
endif
ierr=NF90_DEF_VAR(ncid,"time",NF90_FLOAT,(/ time_dim_id /), time_var_id)
ierr=NF90_DEF_VAR(ncid,"time",NF90_FLOAT,(/ time_dim_id /), time_var_id)
ierr=NF90_PUT_ATT(ncid,time_var_id,"long_name","time")
ierr=NF90_PUT_ATT(ncid,time_var_id,"units","hours since 2014-08-01 00:00:00")
ierr=NF90_PUT_ATT(ncid,time_var_id,"cartesian_axis","T")
ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar_type","JULIAN")
ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar","JULIAN")
if (Model%do_sppt)then
   ierr=NF90_DEF_VAR(ncid,"sppt_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), varid1)
   ierr=NF90_PUT_ATT(ncid,varid1,"long_name","sppt pattern")
   ierr=NF90_PUT_ATT(ncid,varid1,"units","None")
   ierr=NF90_PUT_ATT(ncid,varid1,"missing_value",undef)
   ierr=NF90_PUT_ATT(ncid,varid1,"_FillValue",undef)
   ierr=NF90_PUT_ATT(ncid,varid1,"cell_methods","time: point")
endif
if (Model%do_shum)then
   ierr=NF90_DEF_VAR(ncid,"shum_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), varid2)
   ierr=NF90_PUT_ATT(ncid,varid2,"long_name","shum pattern")
   ierr=NF90_PUT_ATT(ncid,varid2,"units","None")
   ierr=NF90_PUT_ATT(ncid,varid2,"missing_value",undef)
   ierr=NF90_PUT_ATT(ncid,varid2,"_FillValue",undef)
   ierr=NF90_PUT_ATT(ncid,varid2,"cell_methods","time: point")
endif
if (Model%do_skeb)then
   ierr=NF90_DEF_VAR(ncid,"skebu_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,zt_dim_id,time_dim_id/), varid3)
   ierr=NF90_DEF_VAR(ncid,"skebv_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,zt_dim_id,time_dim_id/), varid4)
   ierr=NF90_PUT_ATT(ncid,varid3,"long_name","skeb u pattern")
   ierr=NF90_PUT_ATT(ncid,varid3,"units","None")
   ierr=NF90_PUT_ATT(ncid,varid3,"missing_value",undef)
   ierr=NF90_PUT_ATT(ncid,varid3,"_FillValue",undef)
   ierr=NF90_PUT_ATT(ncid,varid3,"cell_methods","time: point")
   ierr=NF90_PUT_ATT(ncid,varid4,"long_name","skeb v pattern")
   ierr=NF90_PUT_ATT(ncid,varid4,"units","None")
   ierr=NF90_PUT_ATT(ncid,varid4,"missing_value",undef)
   ierr=NF90_PUT_ATT(ncid,varid4,"_FillValue",undef)
   ierr=NF90_PUT_ATT(ncid,varid4,"cell_methods","time: point")
endif
ierr=NF90_ENDDEF(ncid)
ierr=NF90_PUT_VAR(ncid,xt_var_id,grid_xt)
ierr=NF90_PUT_VAR(ncid,yt_var_id,grid_yt)
if (Model%do_skeb)then
   ierr=NF90_PUT_VAR(ncid,zt_var_id,pressl)
endif
endif
! put lat lon and tile number
ierr=NF90_PUT_VAR(ncid,var_id_lon,Init_parm%xlon,(/1,1,1/))
ierr=NF90_PUT_VAR(ncid,var_id_lat,Init_parm%xlat,(/1,1,1/))
tile_number=my_id+1
ierr=NF90_PUT_VAR(ncid,var_id_tile,tile_number,(/1,1,1/))
do i=1,nblks
   if (Model%do_sppt)allocate(Coupling(i)%sppt_wts(blksz,nlevs))
   if (Model%do_shum)allocate(Coupling(i)%shum_wts(blksz,nlevs))
   if (Model%do_skeb)allocate(Coupling(i)%skebu_wts(blksz,nlevs))
   if (Model%do_skeb)allocate(Coupling(i)%skebv_wts(blksz,nlevs))
enddo
do i=1,200
   Model%kdt=i
   ts=i/4.0
   call run_stochastic_physics(Model, Grid, Coupling, nthreads)
   if (Model%me.EQ.0) print*,'sppt_wts=',i,Coupling(1)%sppt_wts(1,20)
   if (write_this_tile) then
   if (Model%do_sppt)then
      do j=1,ny
         workg(:,j)=Coupling(j)%sppt_wts(:,20)   
      enddo
      ierr=NF90_PUT_VAR(ncid,varid1,workg,(/1,1,i/))
   endif
   if (Model%do_shum)then
      do j=1,ny
         workg(:,j)=Coupling(j)%shum_wts(:,1)
      enddo
      ierr=NF90_PUT_VAR(ncid,varid2,workg,(/1,1,i/))
   endif
   if (Model%do_skeb)then
      do k=1,nlevs
         do j=1,ny
            workg3d(:,j,k)=Coupling(j)%skebu_wts(:,k)
         enddo
      enddo
      ierr=NF90_PUT_VAR(ncid,varid3,workg3d,(/1,1,1,i/))
      do k=1,nlevs
         do j=1,ny
            workg3d(:,j,k)=Coupling(j)%skebv_wts(:,k)
         enddo
      enddo
      ierr=NF90_PUT_VAR(ncid,varid4,workg3d,(/1,1,1,i/))
   endif
   ierr=NF90_PUT_VAR(ncid,time_var_id,ts,(/i/))
   endif
enddo
if (write_this_tile) ierr=NF90_CLOSE(ncid)
end
subroutine get_outfile(fname)
use stochy_namelist_def
character*80,intent(out) :: fname
character*4   :: s_ntrunc,s_lat,s_lon
   write(s_ntrunc,'(I4)') ntrunc
   write(s_lat,'(I4)') lat_s 
   write(s_lon,'(I4)') lon_s  
   fname=trim('workg_T'//trim(adjustl(s_ntrunc))//'_'//trim(adjustl(s_lon))//'x'//trim(adjustl(s_lat)))
   return
end
