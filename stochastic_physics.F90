!>@brief The module 'stochastic_physics' is for initialization and running of
!! the stochastic physics random pattern generators
module stochastic_physics

implicit none

private

public :: init_stochastic_physics
public :: run_stochastic_physics

contains

!>@brief The subroutine 'init_stochastic_physics' initializes the stochastic
!!pattern genertors
!>@details It reads the stochastic physics namelist (nam_stoch and nam_sfcperts)
!allocates and polulates the necessary arrays

subroutine init_stochastic_physics(Model, Init_parm, ntasks, nthreads)
!\callgraph
use fv_mp_mod, only : is_master
use stochy_internal_state_mod
use stochy_data_mod, only : nshum,rpattern_shum,init_stochdata,rpattern_sppt,nsppt,rpattern_skeb,nskeb,gg_lats,gg_lons,&
                            rad2deg,INTTYP,wlon,rnlat,gis_stochy,vfact_skeb,vfact_sppt,vfact_shum,skeb_vpts,skeb_vwts,sl
use stochy_resol_def , only : latg,lonf,skeblevs
use stochy_gg_def,only : colrad_a
use stochy_namelist_def
use physcons, only: con_pi
use spectral_layout_mod,only:me,ompthreads,nodes
use mpp_mod
#ifdef STOCHY_UNIT_TEST
 use standalone_stochy_module,   only: GFS_control_type, GFS_init_type
# else
 use GFS_typedefs,       only: GFS_control_type, GFS_init_type
#endif

implicit none
type(GFS_control_type),   intent(inout) :: Model
type(GFS_init_type),      intent(in)    :: Init_parm
integer,                  intent(in)    :: ntasks
integer,                  intent(in)    :: nthreads

integer :: nblks
integer :: iret
real*8 :: PRSI(Model%levs),PRSL(Model%levs),dx
real, allocatable :: skeb_vloc(:)
integer :: k,kflip,latghf,blk,k2
character*2::proc

! Set/update shared variables in spectral_layout_mod
ompthreads  = nthreads

! ------------------------------------------

nblks = size(Model%blksz)

! replace
rad2deg=180.0/con_pi
INTTYP=0 ! bilinear interpolation
me=Model%me
nodes=ntasks
gis_stochy%me=me
gis_stochy%nodes=nodes
call init_stochdata(Model%levs,Model%dtp,Model%input_nml_file,Model%fn_nml,Init_parm%nlunit,iret)
! check to see decomposition
!if(Model%isppt_deep == .true.)then
!do_sppt = .true.
!endif
! check namelist entries for consistency
if (Model%do_sppt.neqv.do_sppt) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings do_sppt and sppt'
   stop
else if (Model%do_shum.neqv.do_shum) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings do_shum and shum'
   stop
else if (Model%do_skeb.neqv.do_skeb) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings do_skeb and skeb'
   stop
else if (Model%lndp_type.neqv.lndp_type) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings lndp_type in physics and nam_sfcperts'
   stop
else if (Model%n_var_lndp .ne. n_var_lndp) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings n_var_lndp in physics nml, and lndp_* in nam_sfcperts'
   stop !  return is not being error-trapped. Need to either kill the execution here, 
        ! or set iret non-zero, then trap it from the calling routine (also, in many other places)
end if
! update remaining model configuration parameters from namelist
Model%use_zmtnblck=use_zmtnblck
Model%skeb_npass=skeb_npass
Model%n_var_lndp=n_var_lndp         ! mg, sfc-perts
Model%lndp_var_list=lndp_var_list
Model%lndp_prt_list=lndp_prt_list
if ( (.NOT. do_sppt) .AND. (.NOT. do_shum) .AND. (.NOT. do_skeb)  .AND. (lndp_type==0) ) return
allocate(sl(Model%levs))
do k=1,Model%levs
   sl(k)= 0.5*(Init_parm%ak(k)/101300.+Init_parm%bk(k)+Init_parm%ak(k+1)/101300.0+Init_parm%bk(k+1)) ! si are now sigmas
!   if(is_master())print*,'sl(k)',k,sl(k),Init_parm%ak(k),Init_parm%bk(k)
enddo
if (Model%do_sppt) then
   allocate(vfact_sppt(Model%levs))
   do k=1,Model%levs
      if (sl(k) .lt. sppt_sigtop1 .and. sl(k) .gt. sppt_sigtop2) then
         vfact_sppt(k) = (sl(k)-sppt_sigtop2)/(sppt_sigtop1-sppt_sigtop2)
      else if (sl(k) .lt. sppt_sigtop2) then
          vfact_sppt(k) = 0.0
      else
          vfact_sppt(k) = 1.0
      endif
   enddo
   if (sppt_sfclimit) then
       vfact_sppt(2)=vfact_sppt(3)*0.5
       vfact_sppt(1)=0.0
   endif
   if (is_master()) then
      do k=1,MOdel%levs
         print *,'sppt vert profile',k,sl(k),vfact_sppt(k)
      enddo
   endif
endif
if (Model%do_skeb) then
   !print*,'allocating skeb stuff',skeblevs
   allocate(vfact_skeb(Model%levs))
   allocate(skeb_vloc(skeblevs)) ! local
   allocate(skeb_vwts(Model%levs,2)) ! save for later
   allocate(skeb_vpts(Model%levs,2)) ! save for later
   do k=1,Model%levs
      if (sl(k) .lt. skeb_sigtop1 .and. sl(k) .gt. skeb_sigtop2) then
         vfact_skeb(k) = (sl(k)-skeb_sigtop2)/(skeb_sigtop1-skeb_sigtop2)
      else if (sl(k) .lt. skeb_sigtop2) then
          vfact_skeb(k) = 0.0
      else
          vfact_skeb(k) = 1.0
      endif
      if (is_master())  print *,'skeb vert profile',k,sl(k),vfact_skeb(k)
   enddo
! calculate vertical interpolation weights
   do k=1,skeblevs
      skeb_vloc(k)=sl(1)-real(k-1)/real(skeblevs-1.0)*(sl(1)-sl(Model%levs))
   enddo
! surface
skeb_vwts(1,2)=0
skeb_vpts(1,1)=1
! top
skeb_vwts(Model%levs,2)=1
skeb_vpts(Model%levs,1)=skeblevs-2
! internal
DO k=2,Model%levs-1
   DO k2=1,skeblevs-1
      IF (sl(k) .LE. skeb_vloc(k2) .AND. sl(k) .GT. skeb_vloc(k2+1)) THEN
        skeb_vpts(k,1)=k2
        skeb_vwts(k,2)=(skeb_vloc(k2)-sl(k))/(skeb_vloc(k2)-skeb_vloc(k2+1))
      ENDIF
   ENDDO
ENDDO
deallocate(skeb_vloc)
if (is_master()) then
DO k=1,Model%levs
   print*,'skeb vpts ',skeb_vpts(k,1),skeb_vwts(k,2)
ENDDO
endif
skeb_vwts(:,1)=1.0-skeb_vwts(:,2)
skeb_vpts(:,2)=skeb_vpts(:,1)+1.0
endif

if (Model%do_shum) then
   allocate(vfact_shum(Model%levs))
   do k=1,Model%levs
      vfact_shum(k) = exp((sl(k)-1.)/shum_sigefold)
      if (sl(k).LT. 2*shum_sigefold) then
         vfact_shum(k)=0.0
      endif
      if (is_master())  print *,'shum vert profile',k,sl(k),vfact_shum(k)
   enddo
endif
! get interpolation weights
! define gaussian grid lats and lons
latghf=latg/2
!print *,'define interp weights',latghf,lonf
!print *,allocated(gg_lats),allocated(gg_lons)
allocate(gg_lats(latg))
!print *,'aloocated lats'
allocate(gg_lons(lonf))
!print *,'aloocated lons'
do k=1,latghf
   gg_lats(k)=-1.0*colrad_a(latghf-k+1)*rad2deg
   gg_lats(latg-k+1)=-1*gg_lats(k)
enddo
dx=360.0/lonf
!print*,'dx=',dx
do k=1,lonf
  gg_lons(k)=dx*(k-1)
enddo
WLON=gg_lons(1)-(gg_lons(2)-gg_lons(1))
RNLAT=gg_lats(1)*2-gg_lats(2)


end subroutine init_stochastic_physics

!>@brief The subroutine 'run_stochastic_physics' updates the random patterns if
!!necessary
!>@details It updates the AR(1) in spectral space
!allocates and polulates the necessary arrays

subroutine run_stochastic_physics(Model, Grid, Coupling, nthreads)
!\callgraph
use fv_mp_mod, only : is_master
use stochy_internal_state_mod
use stochy_data_mod, only : nshum,rpattern_shum,rpattern_sppt,nsppt,rpattern_skeb,nskeb,&
                            rad2deg,INTTYP,wlon,rnlat,gis_stochy,vfact_sppt,vfact_shum,vfact_skeb, & 
                            rpattern_sfc, nlndp
use get_stochy_pattern_mod,only : get_random_pattern_fv3,get_random_pattern_fv3_vect,dump_patterns, & 
                                  get_random_pattern_fv3_sfc
use stochy_resol_def , only : latg,lonf
use stochy_namelist_def
use spectral_layout_mod,only:me,ompthreads
use mpp_mod
#ifdef STOCHY_UNIT_TEST
use standalone_stochy_module,   only: GFS_control_type, GFS_grid_type, GFS_Coupling_type
#else
use GFS_typedefs,       only: GFS_control_type, GFS_grid_type, GFS_Coupling_type
#endif
implicit none
type(GFS_control_type),   intent(in) :: Model
type(GFS_grid_type),      intent(in) :: Grid(:)
type(GFS_coupling_type),  intent(inout) :: Coupling(:)
integer,                  intent(in)    :: nthreads

real,allocatable :: tmp_wts(:,:),tmpu_wts(:,:,:),tmpv_wts(:,:,:), tmps_wts(:,:,:)
!D-grid
integer :: k
integer j,ierr,i
integer :: nblks, blk, len, maxlen
character*120 :: sfile
character*6   :: STRFH
logical :: do_advance_pattern

if ( (.NOT. Model%do_sppt) .AND. (.NOT. Model%do_shum) .AND. (.NOT. Model%do_skeb) .AND. (Model%lndp_type==0 ) ) return

! Update number of threads in shared variables in spectral_layout_mod and set block-related variables
ompthreads = nthreads
nblks = size(Model%blksz)
maxlen = maxval(Model%blksz(:))

! check to see if it is time to write out random patterns
if (fhstoch.GE. 0 .AND. MOD(Model%phour,fhstoch) .EQ. 0) then
   write(STRFH,FMT='(I6.6)') nint(Model%phour)
   sfile='stoch_out.F'//trim(STRFH)
   call dump_patterns(sfile)
endif
allocate(tmp_wts(nblks,maxlen))
allocate(tmpu_wts(nblks,maxlen,Model%levs))
allocate(tmpv_wts(nblks,maxlen,Model%levs))
allocate(tmps_wts(nblks,maxlen,Model%n_var_lndp))
if (Model%do_sppt) then
   if (mod(Model%kdt,nssppt) == 1 .or. nssppt == 1) then
      call get_random_pattern_fv3(rpattern_sppt,nsppt,gis_stochy,Model,Grid,nblks,maxlen,tmp_wts)
      DO blk=1,nblks
         len=size(Grid(blk)%xlat,1)
         DO k=1,Model%levs
            Coupling(blk)%sppt_wts(:,k)=tmp_wts(blk,1:len)*vfact_sppt(k)
         ENDDO
         if (sppt_logit) Coupling(blk)%sppt_wts(:,:) = (2./(1.+exp(Coupling(blk)%sppt_wts(:,:))))-1.
          Coupling(blk)%sppt_wts(:,:)= Coupling(blk)%sppt_wts(:,:)+1.0
      ENDDO
   endif
endif
if (Model%do_shum) then
   if (mod(Model%kdt,nsshum) == 1 .or. nsshum == 1) then
      call get_random_pattern_fv3(rpattern_shum,nshum,gis_stochy,Model,Grid,nblks,maxlen,tmp_wts)
      DO blk=1,nblks
         len=size(Grid(blk)%xlat,1)
         DO k=1,Model%levs
            Coupling(blk)%shum_wts(:,k)=tmp_wts(blk,1:len)*vfact_shum(k)
         ENDDO
      ENDDO
   endif
endif
if (Model%do_skeb) then
   if (mod(Model%kdt,nsskeb) == 1 .or. nsskeb == 1) then
      call get_random_pattern_fv3_vect(rpattern_skeb,nskeb,gis_stochy,Model,Grid,nblks,maxlen,tmpu_wts,tmpv_wts)
      DO blk=1,nblks
         len=size(Grid(blk)%xlat,1)
         DO k=1,Model%levs
            Coupling(blk)%skebu_wts(:,k)=tmpu_wts(blk,1:len,k)*vfact_skeb(k)
            Coupling(blk)%skebv_wts(:,k)=tmpv_wts(blk,1:len,k)*vfact_skeb(k)
         ENDDO
      ENDDO
   endif
endif
if (lndp_type .GT. 0) then
    if (lndp_type==1) then  ! this scheme applies same pert every time step
        do_advance_pattern=.false. 
    else  ! this scheme evolves the patterns through time
        do_advance_pattern=.true. 
    endif
    call get_random_pattern_fv3_sfc(rpattern_sfc,nlndp,gis_stochy,Model,Grid,nblks,maxlen,do_advance_pattern,tmps_wts)
    DO blk=1,nblks
       len=size(Grid(blk)%xlat,1)
       ! for perturbing vars or states, saved value is N(0,1)  and apply scaling later.
       DO k=1,Model%n_var_lndp
           Coupling(blk)%sfc_wts(:,k) = tmps_wts(blk,1:len,k)
       ENDDO
    ENDDO
endif
deallocate(tmp_wts)
deallocate(tmpu_wts)
deallocate(tmpv_wts)
deallocate(tmps_wts)

end subroutine run_stochastic_physics

end module stochastic_physics
