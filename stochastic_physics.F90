!>@brief The module 'stochastic_physics' is for initialization and running of
!! the stochastic physics random pattern generators
module stochastic_physics

use kinddef, only : kind_dbl_prec

implicit none

private

public :: init_stochastic_physics
public :: run_stochastic_physics

contains

!>@brief The subroutine 'init_stochastic_physics' initializes the stochastic
!!pattern genertors
!>@details It reads the stochastic physics namelist (nam_stoch and nam_sfcperts)
!allocates and polulates the necessary arrays

subroutine init_stochastic_physics(levs, blksz, dtp, input_nml_file_in, fn_nml, nlunit, &
    do_sppt_in, do_shum_in, do_skeb_in, lndp_type_in, n_var_lndp_in, use_zmtnblck_out, skeb_npass_out,    &
    lndp_var_list_out, lndp_prt_list_out, ak, bk, nthreads, mpiroot, mpicomm, iret) 
!\callgraph
!use stochy_internal_state_mod
use stochy_data_mod, only : nshum,rpattern_shum,init_stochdata,rpattern_sppt,nsppt,rpattern_skeb,nskeb,gg_lats,gg_lons,&
                            rad2deg,INTTYP,wlon,rnlat,gis_stochy,vfact_skeb,vfact_sppt,vfact_shum,skeb_vpts,skeb_vwts,sl
use stochy_resol_def , only : latg,lonf,skeblevs
use stochy_gg_def,only : colrad_a
use stochy_namelist_def
use spectral_layout_mod,only:me,master,nodes,ompthreads
use mpi_wrapper, only : mpi_wrapper_initialize,mype,npes,is_master

implicit none
integer, intent(out)                    :: iret

! Interface variables

integer,                  intent(in)    :: levs, nlunit, nthreads, mpiroot, mpicomm
integer,                  intent(in)    :: blksz(:)
real(kind=kind_dbl_prec), intent(in)    :: dtp
character(len=*),         intent(in)    :: input_nml_file_in(:)
character(len=*),         intent(in)    :: fn_nml
logical,                  intent(in)    :: do_sppt_in, do_shum_in, do_skeb_in 
integer,                  intent(in)    :: lndp_type_in, n_var_lndp_in
real(kind=kind_dbl_prec), intent(in)    :: ak(:), bk(:) 
logical,                  intent(out)   :: use_zmtnblck_out
integer,                  intent(out)   :: skeb_npass_out
character(len=3), dimension(max_n_var_lndp),         intent(out) :: lndp_var_list_out
real(kind=kind_dbl_prec), dimension(max_n_var_lndp), intent(out) :: lndp_prt_list_out


! Local variables
real(kind=kind_dbl_prec), parameter     :: con_pi =4.0d0*atan(1.0d0)
integer :: nblks
real*8 :: PRSI(levs),PRSL(levs),dx
real, allocatable :: skeb_vloc(:)
integer :: k,kflip,latghf,blk,k2
character*2::proc

! Initialize MPI and OpenMP
call mpi_wrapper_initialize(mpiroot,mpicomm)
me         = mype
nodes      = npes
master     = mpiroot
ompthreads = nthreads

! ------------------------------------------

nblks = size(blksz)

! replace
rad2deg=180.0/con_pi
INTTYP=0 ! bilinear interpolation
gis_stochy%me=me
gis_stochy%nodes=nodes
call init_stochdata(levs,dtp,input_nml_file_in,fn_nml,nlunit,iret)
if (iret .ne. 0) return
! check namelist entries for consistency
if (do_sppt_in.neqv.do_sppt) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings do_sppt and sppt'
   iret = 20 
   return
else if (do_shum_in.neqv.do_shum) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings do_shum and shum'
   iret = 20 
   return
else if (do_skeb_in.neqv.do_skeb) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings do_skeb and skeb'
   iret = 20 
   return
else if (lndp_type_in /= lndp_type) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings lndp_type in physics and nam_sfcperts'
   iret = 20 
   return
else if (n_var_lndp_in /=  n_var_lndp) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings n_var_lndp in physics nml, and lndp_* in nam_sfcperts'
   iret = 20 
   return
end if
! update remaining model configuration parameters from namelist
use_zmtnblck_out=use_zmtnblck
skeb_npass_out=skeb_npass
lndp_var_list_out=lndp_var_list
lndp_prt_list_out=lndp_prt_list
if ( (.NOT. do_sppt) .AND. (.NOT. do_shum) .AND. (.NOT. do_skeb)  .AND. (lndp_type==0) ) return
allocate(sl(levs))
do k=1,levs
   sl(k)= 0.5*(ak(k)/101300.+bk(k)+ak(k+1)/101300.0+bk(k+1)) ! si are now sigmas
enddo
if (do_sppt) then
   allocate(vfact_sppt(levs))
   do k=1,levs
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
      do k=1,levs
         print *,'sppt vert profile',k,sl(k),vfact_sppt(k)
      enddo
   endif
endif
if (do_skeb) then
   !print*,'allocating skeb stuff',skeblevs
   allocate(vfact_skeb(levs))
   allocate(skeb_vloc(skeblevs)) ! local
   allocate(skeb_vwts(levs,2)) ! save for later
   allocate(skeb_vpts(levs,2)) ! save for later
   do k=1,levs
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
      skeb_vloc(k)=sl(1)-real(k-1)/real(skeblevs-1.0)*(sl(1)-sl(levs))
   enddo
! surface
skeb_vwts(1,2)=0
skeb_vpts(1,1)=1
! top
skeb_vwts(levs,2)=1
skeb_vpts(levs,1)=skeblevs-2
! internal
DO k=2,levs-1
   DO k2=1,skeblevs-1
      IF (sl(k) .LE. skeb_vloc(k2) .AND. sl(k) .GT. skeb_vloc(k2+1)) THEN
        skeb_vpts(k,1)=k2
        skeb_vwts(k,2)=(skeb_vloc(k2)-sl(k))/(skeb_vloc(k2)-skeb_vloc(k2+1))
      ENDIF
   ENDDO
ENDDO
deallocate(skeb_vloc)
if (is_master()) then
DO k=1,levs
   print*,'skeb vpts ',skeb_vpts(k,1),skeb_vwts(k,2)
ENDDO
endif
skeb_vwts(:,1)=1.0-skeb_vwts(:,2)
skeb_vpts(:,2)=skeb_vpts(:,1)+1.0
endif

if (do_shum) then
   allocate(vfact_shum(levs))
   do k=1,levs
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

!print *,'done with init_stochastic_physics'

end subroutine init_stochastic_physics

!>@brief The subroutine 'run_stochastic_physics' updates the random patterns if
!!necessary
!>@details It updates the AR(1) in spectral space
!allocates and polulates the necessary arrays

subroutine run_stochastic_physics(levs, kdt, phour, blksz, xlat, xlon, sppt_wts, shum_wts, skebu_wts,  & 
                                  skebv_wts, sfc_wts,nthreads)

!\callgraph
!use stochy_internal_state_mod
use stochy_data_mod, only : nshum,rpattern_shum,rpattern_sppt,nsppt,rpattern_skeb,nskeb,&
                            rad2deg,INTTYP,wlon,rnlat,gis_stochy,vfact_sppt,vfact_shum,vfact_skeb, & 
                            rpattern_sfc, nlndp
use get_stochy_pattern_mod,only : get_random_pattern_fv3,get_random_pattern_fv3_vect,dump_patterns, & 
                                  get_random_pattern_fv3_sfc
use stochy_resol_def , only : latg,lonf
use stochy_namelist_def, only : do_shum,do_sppt,do_skeb,fhstoch,nssppt,nsshum,nsskeb,sppt_logit,    & 
                                lndp_type, n_var_lndp
use mpi_wrapper, only: is_master
use spectral_layout_mod,only:ompthreads
implicit none

! Interface variables
integer,                  intent(in) :: levs, kdt
real(kind=kind_dbl_prec), intent(in) :: phour
integer,                  intent(in) :: blksz(:)
real(kind=kind_dbl_prec), intent(in) :: xlat(:,:)
real(kind=kind_dbl_prec), intent(in) :: xlon(:,:)
real(kind=kind_dbl_prec), intent(inout) :: sppt_wts(:,:,:)
real(kind=kind_dbl_prec), intent(inout) :: shum_wts(:,:,:)
real(kind=kind_dbl_prec), intent(inout) :: skebu_wts(:,:,:)
real(kind=kind_dbl_prec), intent(inout) :: skebv_wts(:,:,:)
real(kind=kind_dbl_prec), intent(inout) :: sfc_wts(:,:,:)
integer,                  intent(in)    :: nthreads

real,allocatable :: tmp_wts(:,:),tmpu_wts(:,:,:),tmpv_wts(:,:,:), tmpl_wts(:,:,:)
!D-grid
integer :: k
integer j,ierr,i
integer :: nblks, blk, len, maxlen
character*120 :: sfile
character*6   :: STRFH
logical :: do_advance_pattern

if ( (.NOT. do_sppt) .AND. (.NOT. do_shum) .AND. (.NOT. do_skeb) .AND. (lndp_type==0 ) ) return

! Update number of threads in shared variables in spectral_layout_mod and set block-related variables
ompthreads = nthreads
nblks = size(blksz)
maxlen = maxval(blksz(:))


if ( (lndp_type==1) .and. (kdt==0) ) then ! old land pert scheme called once at start
        write(0,*) 'calling get_random_pattern_fv3_sfc'  
        allocate(tmpl_wts(nblks,maxlen,n_var_lndp))
        call get_random_pattern_fv3_sfc(rpattern_sfc,nlndp,gis_stochy,xlat,xlon,blksz,nblks,maxlen,tmpl_wts)
        DO blk=1,nblks
           len=blksz(blk)
           ! for perturbing vars or states, saved value is N(0,1)  and apply scaling later.
           DO k=1,n_var_lndp
               sfc_wts(blk,1:len,k) = tmpl_wts(blk,1:len,k)
           ENDDO
        ENDDO
        deallocate(tmpl_wts)
else 
    ! check to see if it is time to write out random patterns
    if (fhstoch.GE. 0 .AND. MOD(phour,fhstoch) .EQ. 0) then
       write(STRFH,FMT='(I6.6)') nint(phour)
       sfile='stoch_out.F'//trim(STRFH)
       call dump_patterns(sfile)
    endif
    allocate(tmp_wts(nblks,maxlen))
    allocate(tmpu_wts(nblks,maxlen,levs))
    allocate(tmpv_wts(nblks,maxlen,levs))
    if (do_sppt) then
       if (mod(kdt,nssppt) == 1 .or. nssppt == 1) then
          call get_random_pattern_fv3(rpattern_sppt,nsppt,gis_stochy,xlat,xlon,blksz,nblks,maxlen,tmp_wts)
          DO blk=1,nblks
             len=blksz(blk)
             DO k=1,levs
                sppt_wts(blk,1:len,k)=tmp_wts(blk,1:len)*vfact_sppt(k)
             ENDDO
             if (sppt_logit) sppt_wts(blk,:,:) = (2./(1.+exp(sppt_wts(blk,:,:))))-1.
             sppt_wts(blk,:,:) = sppt_wts(blk,:,:)+1.0
          ENDDO
       endif
    endif
    if (do_shum) then
       if (mod(kdt,nsshum) == 1 .or. nsshum == 1) then
          call get_random_pattern_fv3(rpattern_shum,nshum,gis_stochy,xlat,xlon,blksz,nblks,maxlen,tmp_wts)
          DO blk=1,nblks
             len=blksz(blk)
             DO k=1,levs
                shum_wts(blk,1:len,k)=tmp_wts(blk,1:len)*vfact_shum(k)
             ENDDO
          ENDDO
       endif
    endif
    if (do_skeb) then
       if (mod(kdt,nsskeb) == 1 .or. nsskeb == 1) then
          call get_random_pattern_fv3_vect(rpattern_skeb,nskeb,gis_stochy,levs,xlat,xlon,blksz,nblks,maxlen,tmpu_wts,tmpv_wts)
          DO blk=1,nblks
             len=blksz(blk)
             DO k=1,levs
                skebu_wts(blk,1:len,k)=tmpu_wts(blk,1:len,k)*vfact_skeb(k)
                skebv_wts(blk,1:len,k)=tmpv_wts(blk,1:len,k)*vfact_skeb(k)
             ENDDO
          ENDDO
       endif
    endif
    if ( lndp_type .EQ. 2  ) then 
        ! add time check?
        allocate(tmpl_wts(nblks,maxlen,n_var_lndp))
        call get_random_pattern_fv3_sfc(rpattern_sfc,nlndp,gis_stochy,xlat,xlon,blksz,nblks,maxlen,tmpl_wts)
        DO blk=1,nblks
           len=blksz(blk)
           ! for perturbing vars or states, saved value is N(0,1)  and apply scaling later.
           DO k=1,n_var_lndp
               sfc_wts(blk,1:len,k) = tmpl_wts(blk,1:len,k)
           ENDDO
        ENDDO
        deallocate(tmpl_wts)
    endif
    deallocate(tmp_wts)
    deallocate(tmpu_wts)
    deallocate(tmpv_wts)

end if 

end subroutine run_stochastic_physics

end module stochastic_physics
