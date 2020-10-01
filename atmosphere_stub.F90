module atmosphere_stub_mod

#include <fms_platform.h>

!-----------------
! FMS modules:
!-----------------
use time_manager_mod,       only: time_type, get_time, set_time, operator(+), &
                                  operator(-), operator(/), time_type_to_real
use fms_mod,                only: file_exist, open_namelist_file,    &
                                  close_file, error_mesg, FATAL,     &
                                  check_nml_error, stdlog,           &
                                  write_version_number,              &
                                  set_domain,                        &
                                  read_data,                         &
                                  mpp_clock_id, mpp_clock_begin,     &
                                  mpp_clock_end, CLOCK_SUBCOMPONENT, &
                                  clock_flag_default, nullify_domain
use mpp_mod,                only: mpp_error, stdout, FATAL, NOTE,    &
                                  input_nml_file, mpp_root_pe,       &
                                  mpp_npes, mpp_pe, mpp_chksum,      &
                                  mpp_get_current_pelist,            &
                                  mpp_set_current_pelist
use mpp_parameter_mod,      only: EUPDATE, WUPDATE, SUPDATE, NUPDATE
use mpp_domains_mod,        only: domain2d, mpp_update_domains
use xgrid_mod,              only: grid_box_type

!-----------------
! FV core modules:
!-----------------
use fv_arrays_mod,      only: fv_atmos_type,fv_grid_bounds_type,fv_grid_type
use fv_control_stub_mod,     only: fv_init, ngrids
use fv_timing_mod,      only: timing_on, timing_off
use fv_sg_mod,          only: fv_subgrid_z
use fv_update_phys_mod, only: fv_update_phys
use mpp_domains_mod,    only:  mpp_get_data_domain, mpp_get_compute_domain
use tp_core_mod,        only: copy_corners
use a2b_edge_mod,       only: a2b_ord4

implicit none
private

!--- driver routines
public :: atmosphere_init_stub

!--- utility routines
!public ::  atmosphere_return_winds, atmosphere_smooth_noise
public ::  atmosphere_resolution,atmosphere_domain,&
           atmosphere_scalar_field_halo,atmosphere_control_data

!--- physics/radiation data exchange routines

!-----------------------------------------------------------------------
! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>
character(len=20)   :: mod_name = 'fvGFS/atmosphere_mod'

!---- private data ----
  public Atm, mytile

  !These are convenience variables for local use only, and are set to values in Atm%
  real    :: dt_atmos
  integer :: npx, npy, npz, ncnst, pnats
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: sec, seconds, days
  integer :: id_dynam, id_fv_diag, id_subgridz


  integer :: mytile  = 1
  integer :: p_split = 1
  integer, allocatable :: pelist(:)
  logical, allocatable :: grids_on_this_pe(:)
  type(fv_atmos_type), allocatable, target :: Atm(:)

  integer :: id_udt_dyn, id_vdt_dyn


!---dynamics tendencies for use in fv_subgrid_z and during fv_update_phys
  real, allocatable, dimension(:,:,:)   :: u_dt, v_dt, t_dt
  real, allocatable                     :: pref(:,:), dum1d(:)

  logical :: first_diag = .true.

contains


!>@brief The subroutine 'atmosphere_init' is an API to initialize the FV3 dynamical core,
!! including the grid structures, memory, initial state (self-initialization or restart),
!! and diagnostics.
 subroutine atmosphere_init_stub (Grid_box, area)
   type(grid_box_type), intent(inout) :: Grid_box
   real*8, pointer, dimension(:,:), intent(inout) :: area
!--- local variables ---
   integer :: i, n

                    call timing_on('ATMOS_INIT')
   allocate(pelist(mpp_npes()))
   call mpp_get_current_pelist(pelist)


!---- compute physics/atmos time step in seconds ----

   dt_atmos = real(sec)

   call fv_init( Atm, dt_atmos, grids_on_this_pe, p_split )  ! allocates Atm components

   do n=1,ngrids
      if (grids_on_this_pe(n)) mytile = n
   enddo

   npx   = Atm(mytile)%npx
   npy   = Atm(mytile)%npy
   npz   = Atm(mytile)%npz
   ncnst = Atm(mytile)%ncnst
   pnats = Atm(mytile)%flagstruct%pnats

   isc = Atm(mytile)%bd%isc
   iec = Atm(mytile)%bd%iec
   jsc = Atm(mytile)%bd%jsc
   jec = Atm(mytile)%bd%jec

   isd = isc - Atm(mytile)%bd%ng
   ied = iec + Atm(mytile)%bd%ng
   jsd = jsc - Atm(mytile)%bd%ng
   jed = jec + Atm(mytile)%bd%ng



   ! Allocate grid variables to be used to calculate gradient in 2nd order flux exchange
   ! This data is only needed for the COARSEST grid.

   allocate(Grid_box%dx    (   isc:iec  , jsc:jec+1))
   allocate(Grid_box%dy    (   isc:iec+1, jsc:jec  ))
   allocate(Grid_box%area  (   isc:iec  , jsc:jec  ))
   allocate(Grid_box%edge_w(              jsc:jec+1))
   allocate(Grid_box%edge_e(              jsc:jec+1))
   allocate(Grid_box%edge_s(   isc:iec+1           ))
   allocate(Grid_box%edge_n(   isc:iec+1           ))
   allocate(Grid_box%en1   (3, isc:iec  , jsc:jec+1))
   allocate(Grid_box%en2   (3, isc:iec+1, jsc:jec  ))
   allocate(Grid_box%vlon  (3, isc:iec  , jsc:jec  ))
   allocate(Grid_box%vlat  (3, isc:iec  , jsc:jec  ))
   Grid_box%dx    (   isc:iec  , jsc:jec+1) = Atm(mytile)%gridstruct%dx    (   isc:iec,   jsc:jec+1)
   Grid_box%dy    (   isc:iec+1, jsc:jec  ) = Atm(mytile)%gridstruct%dy    (   isc:iec+1, jsc:jec  )
   Grid_box%area  (   isc:iec  , jsc:jec  ) = Atm(mytile)%gridstruct%area  (   isc:iec  , jsc:jec  )
   Grid_box%edge_w(              jsc:jec+1) = Atm(mytile)%gridstruct%edge_w(              jsc:jec+1)
   Grid_box%edge_e(              jsc:jec+1) = Atm(mytile)%gridstruct%edge_e(              jsc:jec+1)
   Grid_box%edge_s(   isc:iec+1           ) = Atm(mytile)%gridstruct%edge_s(   isc:iec+1)
   Grid_box%edge_n(   isc:iec+1           ) = Atm(mytile)%gridstruct%edge_n(   isc:iec+1)
   Grid_box%en1   (:, isc:iec  , jsc:jec+1) = Atm(mytile)%gridstruct%en1   (:, isc:iec  , jsc:jec+1)
   Grid_box%en2   (:, isc:iec+1, jsc:jec  ) = Atm(mytile)%gridstruct%en2   (:, isc:iec+1, jsc:jec  )
   do i = 1,3
     Grid_box%vlon(i, isc:iec  , jsc:jec  ) = Atm(mytile)%gridstruct%vlon  (isc:iec ,  jsc:jec, i )
     Grid_box%vlat(i, isc:iec  , jsc:jec  ) = Atm(mytile)%gridstruct%vlat  (isc:iec ,  jsc:jec, i )
   enddo
   allocate (area(isc:iec  , jsc:jec  ))
   area(isc:iec,jsc:jec) = Atm(mytile)%gridstruct%area_64(isc:iec,jsc:jec)


   call set_domain ( Atm(mytile)%domain )

!----- initialize atmos_axes and fv_dynamics diagnostics
       !I've had trouble getting this to work with multiple grids at a time; worth revisiting?
!  --- initialize clocks for dynamics, physics_down and physics_up
   id_dynam     = mpp_clock_id ('FV dy-core',  flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )
   id_fv_diag   = mpp_clock_id ('FV Diag',     flags = clock_flag_default, grain=CLOCK_SUBCOMPONENT )

                    call timing_off('ATMOS_INIT')


 end subroutine atmosphere_init_stub

! subroutine atmosphere_smooth_noise (wnoise,npass,ns_type,renorm_type)
!
!   !--- interface variables ---
!   real,intent(inout)     :: wnoise(isd:ied,jsd:jed,1)
!   integer, intent(in) :: npass,ns_type,renorm_type
!   !--- local variables
!   integer:: i,j,nloops,nlast
!   real ::inflation(isc:iec,jsc:jec),inflation2
!   ! scale factor for restoring inflation
!   ! logic:
!   ! if box mean: scalar get basic scaling, vector gets 1/grid dependent scaling  0-0 ; 0 - 1
!   ! if box mean2: no scaling
!   ! if del2   : scalar gets grid dependent scaling,vector get basic scaling  1  0; 1 1
!   if(npass.GT.0) then
!      if (ns_type.NE.2) then
!         if (ns_type.EQ. 0) then
!             !inflation2=1.0/sqrt(1.0/(4.0*npass))
!             inflation2=1.0/sqrt(1.0/(9.0*npass))
!         else
!             inflation2=1.0/sqrt(1.0/(11.0/3.0*npass))
!         endif
!        if ( ns_type.EQ.1) then ! del2 smoothing needs to be scaled by grid-size
!           do j=jsc,jec
!              do i=isc,iec
!                 inflation(i,j)=inflation2*Atm(mytile)%gridstruct%dxAV/(0.5*(Atm(mytile)%gridstruct%dx(i,j)+Atm(mytile)%gridstruct%dy(i,j)))
!              enddo
!           enddo
!        else
!           if ( renorm_type.EQ.1) then  ! box smooth does not need scaling for scalar
!               do j=jsc,jec
!                  do i=isc,iec
!                   inflation(i,j)=inflation2
!                  enddo
!               enddo
!           else
!              ! box mean needs inversize grid-size scaling for vector
!              do j=jsc,jec
!                 do i=isc,iec
!                    inflation(i,j)=inflation2*(0.5*(Atm(mytile)%gridstruct%dx(i,j)+Atm(mytile)%gridstruct%dy(i,j)))/Atm(mytile)%gridstruct%dxAV
!                 enddo
!              enddo
!           endif
!        endif
!     endif
!     nloops=npass/3
!     nlast=mod(npass,3)
!     do j=1,nloops
!        if (ns_type.EQ.1) then
!           !call del2_cubed(wnoise , 0.25*Atm(mytile)%gridstruct%da_min, Atm(mytile)%gridstruct, &
!           call del2_cubed(wnoise , 0.20*Atm(mytile)%gridstruct%da_min, Atm(mytile)%gridstruct, &
!                           Atm(mytile)%domain, npx, npy, 1, 3, Atm(mytile)%bd)
!        else if (ns_type .EQ. 0) then
!           call box_mean(wnoise , Atm(mytile)%gridstruct, Atm(mytile)%domain, Atm(mytile)%npx, Atm(mytile)%npy, 1, 3, Atm(mytile)%bd)
!        else if (ns_type .EQ. 2) then
!           call box_mean2(wnoise , Atm(mytile)%gridstruct, Atm(mytile)%domain, Atm(mytile)%npx, Atm(mytile)%npy, 1, 3, Atm(mytile)%bd)
!        endif
!     enddo
!     if(nlast>0) then
!        if (ns_type.EQ.1) then
!           !call del2_cubed(wnoise , 0.25*Atm(mytile)%gridstruct%da_min, Atm(mytile)%gridstruct, &
!           call del2_cubed(wnoise , 0.20*Atm(mytile)%gridstruct%da_min, Atm(mytile)%gridstruct, &
!                           Atm(mytile)%domain, npx, npy, 1, nlast, Atm(mytile)%bd)
!        else if (ns_type .EQ. 0) then
!           call box_mean(wnoise , Atm(mytile)%gridstruct, Atm(mytile)%domain, Atm(mytile)%npx, Atm(mytile)%npy, 1, nlast, Atm(mytile)%bd)
!        else if (ns_type .EQ. 2) then
!           call box_mean2(wnoise , Atm(mytile)%gridstruct, Atm(mytile)%domain, Atm(mytile)%npx, Atm(mytile)%npy, 1, nlast, Atm(mytile)%bd)
!        endif
!     endif
!  ! restore amplitude
!    if (ns_type.NE.2) then
!       do j=jsc,jec
!          do i=isc,iec
!             wnoise(i,j,1)=wnoise(i,j,1)*inflation(i,j)
!          enddo
!       enddo
!      endif
!   endif
! end subroutine atmosphere_smooth_noise

! subroutine atmosphere_return_winds (psi,ua,va,edge,km,vwts)
! integer,intent(in) :: edge
! real,intent(inout) :: psi(isd:ied,jsd:jed)
! real,intent(inout) :: ua(isc:iec+edge,jsc:jec)
! real,intent(inout) :: va(isc:iec,jsc:jec+edge)
! integer, optional,intent(in):: km
! real, optional,intent(in):: vwts(:)
! integer :: k
! call timing_on('COMM_TOTAL')
! call mpp_update_domains(psi, Atm(mytile)%domain, complete=.true.)
! call timing_off('COMM_TOTAL')
! if (edge.EQ.0) then
!    call make_a_winds(ua, va, psi,Atm(mytile)%ng,Atm(mytile)%gridstruct,Atm(mytile)%bd,Atm(mytile)%npx,Atm(mytile)%npy)
! endif
! if (edge.EQ.1) then
!    call make_c_winds(ua, va, psi,Atm(mytile)%ng,Atm(mytile)%gridstruct,Atm(mytile)%bd,Atm(mytile)%npx,Atm(mytile)%npy)
!! populate wind perturbations right here
!    do k=1,km
!       Atm(mytile)%urandom_c(isc:iec+edge,jsc:jec     ,k)=ua*vwts(k)
!       Atm(mytile)%vrandom_c(isc:iec     ,jsc:jec+edge,k)=va*vwts(k)
!    enddo
!    !call mpp_update_domains(Atm(mytile)%urandom_c, Atm(mytile)%domain, complete=.true.)
!    !call mpp_update_domains(Atm(mytile)%vrandom_c, Atm(mytile)%domain, complete=.true.)
! endif
! end subroutine atmosphere_return_winds
!
  subroutine del2_cubed(q, cd, gridstruct, domain, npx, npy, km, nmax, bd)
      !---------------------------------------------------------------
      ! This routine is for filtering the omega field for the physics
      !---------------------------------------------------------------
      integer, intent(in):: npx, npy, km, nmax
      real(kind=8),   intent(in):: cd            !< cd = K * da_min;   0 < K < 0.25
      type(fv_grid_bounds_type), intent(IN) :: bd
      real, intent(inout):: q(bd%isd:bd%ied,bd%jsd:bd%jed,km)
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain
      real, parameter:: r3  = 1./3.
      real :: fx(bd%isd:bd%ied+1,bd%jsd:bd%jed), fy(bd%isd:bd%ied,bd%jsd:bd%jed+1)
      real :: q2(bd%isd:bd%ied,bd%jsd:bd%jed)
      integer i,j,k, n, nt, ntimes
      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      !Local routine pointers
!     real, pointer, dimension(:,:) :: rarea
!     real, pointer, dimension(:,:) :: del6_u, del6_v
!     logical, pointer :: sw_corner, se_corner, ne_corner, nw_corner

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

      ntimes = min(3, nmax)

      call timing_on('COMM_TOTAL')
      call mpp_update_domains(q, domain, complete=.true.)
      call timing_off('COMM_TOTAL')


      do n=1,ntimes
         nt = ntimes - n

!$OMP parallel do default(none) shared(km,q,is,ie,js,je,npx,npy, &
!$OMP                                  nt,isd,jsd,gridstruct,bd, &
!$OMP                                  cd) &
!$OMP                          private(fx, fy)
         do k=1,km

            if ( gridstruct%sw_corner ) then
               q(1,1,k) = (q(1,1,k)+q(0,1,k)+q(1,0,k)) * r3
               q(0,1,k) =  q(1,1,k)
               q(1,0,k) =  q(1,1,k)
            endif
            if ( gridstruct%se_corner ) then
               q(ie, 1,k) = (q(ie,1,k)+q(npx,1,k)+q(ie,0,k)) * r3
               q(npx,1,k) =  q(ie,1,k)
               q(ie, 0,k) =  q(ie,1,k)
            endif
            if ( gridstruct%ne_corner ) then
               q(ie, je,k) = (q(ie,je,k)+q(npx,je,k)+q(ie,npy,k)) * r3
               q(npx,je,k) =  q(ie,je,k)
               q(ie,npy,k) =  q(ie,je,k)
            endif
            if ( gridstruct%nw_corner ) then
               q(1, je,k) = (q(1,je,k)+q(0,je,k)+q(1,npy,k)) * r3
               q(0, je,k) =  q(1,je,k)
               q(1,npy,k) =  q(1,je,k)
            endif

            if(nt>0 .and. (.not. gridstruct%regional)) call copy_corners(q(isd,jsd,k), npx, npy, 1, gridstruct%nested, bd, &
                 gridstruct%sw_corner, gridstruct%se_corner, gridstruct%nw_corner, gridstruct%ne_corner )
            do j=js-nt,je+nt
               do i=is-nt,ie+1+nt
#ifdef USE_SG
                  fx(i,j) = gridstruct%dy(i,j)*gridstruct%sina_u(i,j)*(q(i-1,j,k)-q(i,j,k))*gridstruct%rdxc(i,j)
#else
                  fx(i,j) = gridstruct%del6_v(i,j)*(q(i-1,j,k)-q(i,j,k))
#endif
                enddo
            enddo

            if(nt>0 .and. (.not. gridstruct%regional)) call copy_corners(q(isd,jsd,k), npx, npy, 2, gridstruct%nested, bd, &
                 gridstruct%sw_corner, gridstruct%se_corner, gridstruct%nw_corner, gridstruct%ne_corner)
            do j=js-nt,je+1+nt
               do i=is-nt,ie+nt
#ifdef USE_SG
                  fy(i,j) = gridstruct%dx(i,j)*gridstruct%sina_v(i,j)*(q(i,j-1,k)-q(i,j,k))*gridstruct%rdyc(i,j)
#else
                  fy(i,j) = gridstruct%del6_u(i,j)*(q(i,j-1,k)-q(i,j,k))
#endif
               enddo
            enddo

            do j=js-nt,je+nt
               do i=is-nt,ie+nt
                  q(i,j,k) = q(i,j,k) + cd*gridstruct%rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))
               enddo
            enddo
         enddo
      enddo

 end subroutine del2_cubed

!>@brief The subroutine 'box_mean' filters a field with a 9-point mean stencil

 subroutine box_mean(q, gridstruct, domain, npx, npy, km, nmax, bd)
      !---------------------------------------------------------------
      ! This routine is for filtering the omega field for the physics
      !---------------------------------------------------------------
      integer, intent(in):: npx, npy, km, nmax
      type(fv_grid_bounds_type), intent(IN) :: bd
      real, intent(inout):: q(bd%isd:bd%ied,bd%jsd:bd%jed,km)
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain
      real, parameter:: r3  = 1./3.,r9=1./9.
      real :: q2(bd%isd:bd%ied,bd%jsd:bd%jed)
      integer i,j,k, n, nt, ntimes
      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      !Local routine pointers
!     real, pointer, dimension(:,:) :: rarea
!     real, pointer, dimension(:,:) :: del6_u, del6_v
!     logical, pointer :: sw_corner, se_corner, ne_corner, nw_corner

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

      ntimes = min(3, nmax)

      call timing_on('COMM_TOTAL')
      call mpp_update_domains(q, domain, complete=.true.)
      call timing_off('COMM_TOTAL')


      do n=1,ntimes
         nt = ntimes !- n

!$OMP parallel do default(none) shared(km,is,ie,js,je,npx,npy, &
!$OMP                                  q,nt,isd,jsd,gridstruct,bd) &
!$OMP                          private(q2)
         do k=1,km

            if ( gridstruct%sw_corner ) then
               q(1,1,k) = (q(1,1,k)+q(0,1,k)+q(1,0,k)) * r3
               q(0,1,k) =  q(1,1,k)
               q(1,0,k) =  q(1,1,k)
            endif
            if ( gridstruct%se_corner ) then
               q(ie, 1,k) = (q(ie,1,k)+q(npx,1,k)+q(ie,0,k)) * r3
               q(npx,1,k) =  q(ie,1,k)
               q(ie, 0,k) =  q(ie,1,k)
            endif
            if ( gridstruct%ne_corner ) then
               q(ie, je,k) = (q(ie,je,k)+q(npx,je,k)+q(ie,npy,k)) * r3
               q(npx,je,k) =  q(ie,je,k)
               q(ie,npy,k) =  q(ie,je,k)
            endif
            if ( gridstruct%nw_corner ) then
               q(1, je,k) = (q(1,je,k)+q(0,je,k)+q(1,npy,k)) * r3
               q(0, je,k) =  q(1,je,k)
               q(1,npy,k) =  q(1,je,k)
            endif

            if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 1, gridstruct%nested, bd, &
                 gridstruct%sw_corner, gridstruct%se_corner, gridstruct%nw_corner, gridstruct%ne_corner )

            if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 2, gridstruct%nested, bd, &
                 gridstruct%sw_corner, gridstruct%se_corner, gridstruct%nw_corner, gridstruct%ne_corner)

            !do j=js-nt,je+nt
            !   do i=is-nt,ie+nt
            do j=jsd+1,jed-1
               do i=isd+1,ied-1
                  !q2(i,j) = (gridstruct%area(i-1,j+1)*q(i-1,j+1,k) + gridstruct%area(i,j+1)*q(i,j+1,k) + gridstruct%area(i+1,j+1)*q(i+1,j+1,k) +&
                  !           gridstruct%area(i-1,j  )*q(i-1,j,k)   + gridstruct%area(i,j  )*q(i,j  ,k) + gridstruct%area(i+1,j  )*q(i+1,j  ,k) +&
                  !           gridstruct%area(i-1,j-1)*q(i-1,j-1,k) + gridstruct%area(i,j-1)*q(i,j-1,k) + gridstruct%area(i+1,j-1)*q(i+1,j-1,k))/SUM(gridstruct%area(i-1:i+1,j-1:j+1))
                  q2(i,j) = r9*(q(i-1,j+1,k)+q(i,j+1,k)+q(i+1,j+1,k)+q(i-1,j,k)+q(i,j,k)+q(i+1,j,k)+q(i-1,j-1,k)+q(i,j-1,k)+q(i+1,j-1,k))
                  !if (j.GE. je .AND. i.GE. ie) print*,'area +1=',gridstruct%area(i-1:i+1,j+1)
                  !if (j.GE. je .AND. i.GE. ie) print*,'area   =',gridstruct%area(i-1:i+1,j)
                  !if (j.GE. je .AND. i.GE. ie) print*,'area -1=',gridstruct%area(i-1:i+1,j-1)
                  !if (j.GE. je .AND. i.GE. ie) print*,'q    +1=',q(i-1:i+1,j+1,k)
                  !if (j.GE. je .AND. i.GE. ie) print*,'q      =',q(i-1:i+1,j,k)
                  !if (j.GE. je .AND. i.GE. ie) print*,'q    -1=',q(i-1:i+1,j-1,k)
               enddo
            enddo
            do j=js-nt,je+nt
               do i=is-nt,ie+nt
                  q(i,j,k)=q2(i,j)
               enddo
            enddo
         enddo
      enddo
 end subroutine box_mean

 subroutine box_mean2(q, gridstruct, domain, npx, npy, km, nmax, bd)
      !---------------------------------------------------------------
      ! This routine is for filtering the omega field for the physics
      !---------------------------------------------------------------
      integer, intent(in):: npx, npy, km, nmax
      type(fv_grid_bounds_type), intent(IN) :: bd
      real, intent(inout):: q(bd%isd:bd%ied,bd%jsd:bd%jed,km)
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain
      real, parameter:: r3  = 1./3.,r10=0.1
      real :: q2(bd%isd:bd%ied,bd%jsd:bd%jed)
      integer i,j,k, n, nt, ntimes
      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      !Local routine pointers
!     real, pointer, dimension(:,:) :: rarea
!     real, pointer, dimension(:,:) :: del6_u, del6_v
!     logical, pointer :: sw_corner, se_corner, ne_corner, nw_corner

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

      ntimes = min(3, nmax)

      call timing_on('COMM_TOTAL')
      call mpp_update_domains(q, domain, complete=.true.)
      call timing_off('COMM_TOTAL')


      do n=1,ntimes
         nt = ntimes !- n

!$OMP parallel do default(none) shared(km,is,ie,js,je,npx,npy, &
!$OMP                                  q,nt,isd,jsd,gridstruct,bd) &
!$OMP                          private(q2)
         do k=1,km

            if ( gridstruct%sw_corner ) then
               q(1,1,k) = (q(1,1,k)+q(0,1,k)+q(1,0,k)) * r3
               q(0,1,k) =  q(1,1,k)
               q(1,0,k) =  q(1,1,k)
            endif
            if ( gridstruct%se_corner ) then
               q(ie, 1,k) = (q(ie,1,k)+q(npx,1,k)+q(ie,0,k)) * r3
               q(npx,1,k) =  q(ie,1,k)
               q(ie, 0,k) =  q(ie,1,k)
            endif
            if ( gridstruct%ne_corner ) then
               q(ie, je,k) = (q(ie,je,k)+q(npx,je,k)+q(ie,npy,k)) * r3
               q(npx,je,k) =  q(ie,je,k)
               q(ie,npy,k) =  q(ie,je,k)
            endif
            if ( gridstruct%nw_corner ) then
               q(1, je,k) = (q(1,je,k)+q(0,je,k)+q(1,npy,k)) * r3
               q(0, je,k) =  q(1,je,k)
               q(1,npy,k) =  q(1,je,k)
            endif

            if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 1, gridstruct%nested, bd, &
                 gridstruct%sw_corner, gridstruct%se_corner, gridstruct%nw_corner, gridstruct%ne_corner )

            if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 2, gridstruct%nested, bd, &
                 gridstruct%sw_corner, gridstruct%se_corner, gridstruct%nw_corner, gridstruct%ne_corner)

            do j=jsd+1,jed-1
               do i=isd+1,ied-1
                  q2(i,j) = r10*(q(i-1,j+1,k)+q(i,j+1,k)+q(i+1,j+1,k)+q(i-1,j,k)+2*q(i,j,k)+q(i+1,j,k)+q(i-1,j-1,k)+q(i,j-1,k)+q(i+1,j-1,k))
               enddo
            enddo
            do j=js-nt,je+nt
               do i=is-nt,ie+nt
                  q(i,j,k)=q2(i,j)
               enddo
            enddo
         enddo
      enddo

 end subroutine box_mean2
subroutine make_a_winds(ua, va, psi, ng, gridstruct, bd, npx, npy)

integer, intent(IN) :: ng, npx, npy
type(fv_grid_bounds_type), intent(IN) :: bd
real,    intent(inout) :: psi(bd%isd:bd%ied,bd%jsd:bd%jed)
real, intent(inout) ::      ua(bd%isc:bd%iec  ,bd%jsc:bd%jec )
real, intent(inout) ::      va(bd%isc:bd%iec  ,bd%jsc:bd%jec )
type(fv_grid_type), intent(IN), target :: gridstruct
! Local:
real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed) :: wk
real, dimension(bd%isc:bd%iec,bd%jsc:bd%jec) :: u,v
integer i,j

integer :: is,  ie,  js,  je
is  = bd%is
ie  = bd%ie
js  = bd%js
je  = bd%je

call a2b_ord4( psi, wk, gridstruct, npx, npy, is, ie, js, je, ng, .false.)
do j=js,je
   do i=is,ie
      u(i,j) = gridstruct%rdy(i,j)*0.5*(wk(i,j+1)+wk(i+1,j+1)-(wk(i,j)+wk(i+1,j)))
   enddo
enddo
do j=js,je
   do i=is,ie
      v(i,j) = gridstruct%rdx(i,j)*0.5*(wk(i,j)+wk(i,j+1)-(wk(i+1,j)+wk(i+1,j+1)))
   enddo
enddo
do j=js,je
   do i=is,ie
      ua(i,j) = 0.5*(gridstruct%a11(i,j)+gridstruct%a11(i,j+1))*u(i,j) + 0.5*(gridstruct%a12(i,j)+gridstruct%a12(i,j+1))*v(i,j)
      va(i,j) = 0.5*(gridstruct%a21(i,j)+gridstruct%a21(i+1,j))*u(i,j) + 0.5*(gridstruct%a22(i,j)+gridstruct%a22(i+1,j))*v(i,j)
   enddo
enddo

end subroutine make_a_winds

subroutine make_c_winds(uc, vc, psi, ng, gridstruct, bd, npx, npy)

integer, intent(IN) :: ng, npx, npy
type(fv_grid_bounds_type), intent(IN) :: bd
real,    intent(inout) :: psi(bd%isd:bd%ied,bd%jsd:bd%jed)
real, intent(inout) ::      uc(bd%isc:bd%iec+1 ,bd%jsc:bd%jec )
real, intent(inout) ::      vc(bd%isc:bd%iec   ,bd%jsc:bd%jec+1)
type(fv_grid_type), intent(IN), target :: gridstruct
! Local:
real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed) :: wk
real, dimension(bd%isc:bd%iec,bd%jsc:bd%jec) :: u,v
integer i,j

integer :: is,  ie,  js,  je
is  = bd%is
ie  = bd%ie
js  = bd%js
je  = bd%je

call a2b_ord4( psi, wk, gridstruct, npx, npy, is, ie, js, je, ng, .false.)
do j=js,je
   do i=is,ie+1
      uc(i,j) = gridstruct%rdy(i,j)*(wk(i,j+1)-wk(i,j))
   enddo
enddo
do j=js,je+1
   do i=is,ie
      vc(i,j) = gridstruct%rdx(i,j)*(wk(i,j)-wk(i+1,j))
   enddo
enddo

end subroutine make_c_winds

!>@brief The subroutine 'atmospehre_resolution' is an API to return the local
!! extents of the current MPI-rank or the global extents of the current
!! cubed-sphere tile.
 subroutine atmosphere_resolution (i_size, j_size, global)
   integer, intent(out)          :: i_size, j_size
   logical, intent(in), optional :: global
   logical :: local

   local = .true.
   if( PRESENT(global) ) local = .NOT.global

   if( local ) then
       i_size = iec - isc + 1
       j_size = jec - jsc + 1
   else
       i_size = npx - 1
       j_size = npy - 1
   end if
 end subroutine atmosphere_resolution
!>@brief The subroutine 'atmosphere_domain' is an API to return
!! the "domain2d" variable associated with the coupling grid and the
!! decomposition for the current cubed-sphere tile.
!>@detail Coupling is done using the mass/temperature grid with no halos.
 subroutine atmosphere_domain ( fv_domain, layout, regional, nested, pelist )
   type(domain2d), intent(out) :: fv_domain
   integer, intent(out) :: layout(2)
   logical, intent(out) :: regional
   logical, intent(out) :: nested
   integer, pointer, intent(out) :: pelist(:)
!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

   fv_domain = Atm(mytile)%domain_for_coupler
   layout(1:2) =  Atm(mytile)%layout(1:2)
   regional = Atm(mytile)%flagstruct%regional
   nested = ngrids > 1
   call set_atmosphere_pelist()
   pelist => Atm(mytile)%pelist

 end subroutine atmosphere_domain

 subroutine set_atmosphere_pelist ()
   call mpp_set_current_pelist(Atm(mytile)%pelist, no_sync=.TRUE.)
 end subroutine set_atmosphere_pelist


!>@brief The subroutine 'atmosphere_scalar_field_halo' is an API to return halo information
!! of the current MPI_rank for an input scalar field.
!>@detail Up to three point haloes can be returned by this API which includes special handling for
!! the cubed-sphere tile corners. Output will be in (i,j,k) while input can be in (i,j,k) or
!! horizontally-packed form (ix,k).
 subroutine atmosphere_scalar_field_halo (data, halo, isize, jsize, ksize, data_p)
   !--------------------------------------------------------------------
   ! data   - output array to return the field with halo (i,j,k)
   !          optionally input for field already in (i,j,k) form
   !          sized to include the halo of the field (+ 2*halo)
   ! halo   - size of the halo (must be less than 3)
   ! ied    - horizontal resolution in i-dir with haloes
   ! jed    - horizontal resolution in j-dir with haloes
   ! ksize  - vertical resolution
   ! data_p - optional input field in packed format (ix,k)
   !--------------------------------------------------------------------
   !--- interface variables ---
   real*8, dimension(1:isize,1:jsize,ksize), intent(inout) :: data !< output array to return the field with halo (i,j,k)
                                                                                 !< optionally input for field already in (i,j,k) form
                                                                                 !< sized to include the halo of the field (+ 2*halo)
   integer, intent(in) :: halo  !< size of the halo (must be less than 3)
   integer, intent(in) :: isize !< horizontal resolution in i-dir with haloes
   integer, intent(in) :: jsize !< horizontal resolution in j-dir with haloes
   integer, intent(in) :: ksize !< vertical resolution
   real*8, dimension(:,:), optional, intent(in) :: data_p !< optional input field in packed format (ix,k)
   !--- local variables ---
   integer :: i, j, k
   integer :: ic, jc
   character(len=44) :: modname = 'atmosphere_mod::atmosphere_scalar_field_halo'
   integer :: mpp_flags

   !--- perform error checking
   if (halo .gt. 3) call mpp_error(FATAL, modname//' - halo.gt.3 requires extending the MPP domain to support')
   ic = isize - 2 * halo
   jc = jsize - 2 * halo

   !--- if packed data is present, unpack it into the two-dimensional data array
   if (present(data_p)) then
     if (ic*jc .ne. size(data_p,1)) call mpp_error(FATAL, modname//' - incorrect sizes for incoming &
                                                  &variables data and data_p')
     data = 0.
!$OMP parallel do default (none) &
!$OMP              shared (data, data_p, halo, ic, jc, ksize) &
!$OMP             private (i, j, k)
     do k = 1, ksize
       do j = 1, jc
         do i = 1, ic
           data(i+halo, j+halo, k) = data_p(i + (j-1)*ic, k)
         enddo
       enddo
     enddo
   endif

   mpp_flags = EUPDATE + WUPDATE + SUPDATE + NUPDATE
   if (halo == 1) then
     call mpp_update_domains(data, Atm(mytile)%domain_for_coupler, flags=mpp_flags, complete=.true.)
   elseif (halo == 3) then
     call mpp_update_domains(data, Atm(mytile)%domain, flags=mpp_flags, complete=.true.)
   else
     call mpp_error(FATAL, modname//' - unsupported halo size')
   endif

   !--- fill the halo points when at a corner of the cubed-sphere tile
   !--- interior domain corners are handled correctly
   if ( (isc==1) .or. (jsc==1) .or. (iec==npx-1) .or. (jec==npy-1) ) then
     do k = 1, ksize
       do j=1,halo
         do i=1,halo
           if ((isc==    1) .and. (jsc==    1)) data(halo+1-j ,halo+1-i ,k) = data(halo+i     ,halo+1-j ,k)  !SW Corner
           if ((isc==    1) .and. (jec==npy-1)) data(halo+1-j ,halo+jc+i,k) = data(halo+i     ,halo+jc+j,k)  !NW Corner
           if ((iec==npx-1) .and. (jsc==    1)) data(halo+ic+j,halo+1-i ,k) = data(halo+ic-i+1,halo+1-j ,k)  !SE Corner
           if ((iec==npx-1) .and. (jec==npy-1)) data(halo+ic+j,halo+jc+i,k) = data(halo+ic-i+1,halo+jc+j,k)  !NE Corner
         enddo
       enddo
     enddo
   endif

   return
 end subroutine atmosphere_scalar_field_halo


 subroutine atmosphere_control_data (i1, i2, j1, j2, kt, p_hydro, hydro, tile_num)
   integer, intent(out)           :: i1, i2, j1, j2, kt
   logical, intent(out), optional :: p_hydro, hydro
   integer, intent(out), optional :: tile_num
   i1 = Atm(mytile)%bd%isc
   i2 = Atm(mytile)%bd%iec
   j1 = Atm(mytile)%bd%jsc
   j2 = Atm(mytile)%bd%jec
   kt = Atm(mytile)%npz

   if (present(tile_num)) tile_num = Atm(mytile)%tile
   if (present(p_hydro)) p_hydro   = Atm(mytile)%flagstruct%phys_hydrostatic
   if (present(  hydro))   hydro   = Atm(mytile)%flagstruct%hydrostatic

 end subroutine atmosphere_control_data

end module atmosphere_stub_mod
