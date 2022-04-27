module lndp_apply_perts_mod

    use kinddef, only : kind_dbl_prec
    use stochy_namelist_def

    implicit none

    private

    public :: lndp_apply_perts

    contains

!====================================================================
! lndp_apply_perts
!====================================================================
! Driver for applying perturbations to specified land states or parameters, 
! following the LNDP_TYPE=2 scheme.
! Draper, July 2020.

! Can select perturbations by specifying lndp_var_list and lndp_pert_list

! Notes on how the perturbations are added.
! 1. Model prognostic variables.
! If running a long forecast or cycling DA system (as in global UFS), 
! perturbing the prognostic variables only at the initial conditions will 
! have very limited impact, and they should instead be perturbed at every time step.  
! In this case, the pertrubations should be specified as a rate (pert/hr)
! to avoid the ensemble spread being dependent on the model time step. 
!
! For a short forecast (~days, as in regional HRRR), can see impact from 
! perturbing only the initial conditions. In this case, the perturbation 
! is specified as an absolute value (not a rate). 
! 
! 2. Model parameters: 
! The timing of how to perturb the parameters depends on how / whether 
! the parameters are updated over time. For the UFS global system, global_cycle
! is periodically called to update the parameters (controlled by FHCYC). 
! Each time it's called global_cycle overwrites most of the 
! prior parameters (overwriting any perturbations applied to those
! parameters). Hence, the perturbations are applied only immediately after global_cycle  
! has been called, and the parameters are not applied as a rate (since they 
! don't accumulate).  
! For the regional models, FHCYC is 0, and the global_cycle is not called, so 
! can perturb parameters every time step. Hence, need to specify the perturbations 
! as a rate. 
!
! The above cases are controlled by the lndp_model_type variable. 
!
! If adding perturbations to new parameters, need to check how/whether 
! the parameters are updated by the model.

    subroutine lndp_apply_perts(blksz, lsm, lsm_noah, lsm_ruc, lsm_noahmp, iopt_dveg, & 
                lsoil, dtf, kdt, n_var_lndp, lndp_var_list, lndp_prt_list,            &
                sfc_wts, xlon, xlat, stype, smcmax, smcmin, param_update_flag,  &
                smc, slc, stc, vfrac, alnsf, alnwf, snoalb, semis, zorll, ierr)

        implicit none

        ! intent(in)
        integer,                      intent(in) :: blksz(:)
        integer,                      intent(in) :: n_var_lndp, lsoil, kdt, iopt_dveg
        integer,                      intent(in) :: lsm, lsm_noah, lsm_ruc, lsm_noahmp
        character(len=3),             intent(in) :: lndp_var_list(:)
        real(kind=kind_dbl_prec),     intent(in) :: lndp_prt_list(:)
        real(kind=kind_dbl_prec),     intent(in) :: dtf
        real(kind=kind_dbl_prec),     intent(in) :: sfc_wts(:,:,:)
        real(kind=kind_dbl_prec),     intent(in) :: xlon(:,:)
        real(kind=kind_dbl_prec),     intent(in) :: xlat(:,:)
        logical,                      intent(in) :: param_update_flag
                                        ! true =  parameters have just been updated by global_cycle
        integer,     intent(in) :: stype(:,:)
        real(kind=kind_dbl_prec),     intent(in) :: smcmax(:)
        real(kind=kind_dbl_prec),     intent(in) :: smcmin(:)

        ! intent(inout)
        real(kind=kind_dbl_prec),     intent(inout) :: smc(:,:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: slc(:,:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: stc(:,:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: vfrac(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: snoalb(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: alnsf(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: alnwf(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: semis(:,:)
        real(kind=kind_dbl_prec),     intent(inout) :: zorll(:,:)

        ! intent(out)
        integer,                        intent(out) :: ierr

        ! local
        integer         :: nblks, print_i, print_nb, i, nb
        !integer         :: this_im, v, soiltyp, k
        integer         :: this_im, v, k
        logical         :: print_flag, do_pert_state, do_pert_param

        real(kind=kind_dbl_prec) :: p, min_bound, max_bound, pert
        real(kind=kind_dbl_prec) :: tmp_smc
        real(kind=kind_dbl_prec) :: conv_hr2tstep, tfactor_state, tfactor_param
        real(kind=kind_dbl_prec), dimension(lsoil) :: zslayer, smc_vertscale, stc_vertscale

        ! decrease in applied pert with depth
        !-- Noah lsm
        real(kind=kind_dbl_prec), dimension(4), parameter  :: smc_vertscale_noah = (/1.0,0.5,0.25,0.125/)
        real(kind=kind_dbl_prec), dimension(4), parameter  :: stc_vertscale_noah = (/1.0,0.5,0.25,0.125/)
        real(kind=kind_dbl_prec), dimension(4), parameter  :: zs_noah = (/0.1, 0.3, 0.6, 1.0/)
        !-- RUC lsm
        real(kind=kind_dbl_prec), dimension(9), parameter :: smc_vertscale_ruc = (/1.0,0.9,0.8,0.6,0.4,0.2,0.1,0.05,0./)
        real(kind=kind_dbl_prec), dimension(9), parameter :: stc_vertscale_ruc = (/1.0,0.9,0.8,0.6,0.4,0.2,0.1,0.05,0./)
        real(kind=kind_dbl_prec), dimension(9), parameter :: zs_ruc = (/0.05, 0.15, 0.20, 0.20, 0.40, 0.60, 0.60, 0.80, 1.00/)

        ierr = 0

        if (lsm/=lsm_noah .and. lsm/=lsm_ruc .and. lsm/=lsm_noahmp) then
          write(6,*) 'ERROR: lndp_apply_pert assumes LSM is Noah, Noah-MP, or RUC,', &
                     ' may need to adapt variable names for a different LSM'
          ierr=10
          return
        endif

        if (lsm==lsm_noahmp) then 
            do v = 1,n_var_lndp
                select case (trim(lndp_var_list(v)))
                case ('alb','sal','emi','zol') 
                    print*, &
                     'ERROR:  lndp_prt_list option in lndp_apply_pert', trim(lndp_var_list(v)) , & 
                     ' has not been checked for Noah-MP. Please check how the parameter is set/updated ', & 
                     ' before applying. Note: in Noah-MP many variables that have traditionally been', & 
                     ' externally specified parameters are now prognostic. Also, many parameters are', & 
                     ' set at each Noah-MP model call from data tables' 
                    ierr = 10
                    return
                case ('veg') 
                    if (iopt_dveg .NE. 4 ) then 
                    print*, &
                     'ERROR:  veg perturbations have not yet been coded for dveg options other than 4' 
                      ierr = 10 
                      return 
                    endif
                end select
            enddo 
        endif

        ! for perturbations applied as a rate, lndp_prt_list input is per hour. Converts to per timestep
        conv_hr2tstep = dtf/3600. ! conversion factor from per hour to per tstep.

        ! determine whether updating state variables and/or parameters

        do_pert_state=.false. 
        do_pert_param=.false. 

        select case (lndp_model_type)
        case(1)  ! global, perturb states every time step (pert applied as a rate) 
                 !         perturb parameters only when they've been update (pert is not a rate) 
            do_pert_state=.true.
            tfactor_state=conv_hr2tstep
            if (param_update_flag) then 
                do_pert_param=.true.
                tfactor_param=1.
            endif
        case(2)  ! regional, perturb states only at first time step (pert is not a rate) 
                 !           perurb parameters at every time step (pert is a rate) 
            if ( kdt == 2 ) then 
                do_pert_state=.true.
                tfactor_state=1. 
            endif
            do_pert_param = .true.
            tfactor_param = conv_hr2tstep
        case(3)  ! special case to apply perturbations at initial time step only (pert is not a rate) 
            if ( kdt == 2 ) then 
                do_pert_state=.true.
                tfactor_state=1.
                do_pert_param=.true.
                tfactor_param=1.
            endif
        case default
            print*, &
             'ERROR: unrecognised lndp_model_type option in lndp_apply_pert, exiting', trim(lndp_var_list(v))
            ierr = 10
            return
        end select

        if (lsm == lsm_noah .or. lsm == lsm_noahmp) then
          do k = 1, lsoil
            zslayer(k) = zs_noah(k)
            smc_vertscale(k) = smc_vertscale_noah(k)
            stc_vertscale(k) = stc_vertscale_noah(k)
          enddo
        elseif (lsm == lsm_ruc) then
          do k = 1, lsoil
            zslayer(k) = zs_ruc(k)
            smc_vertscale(k) = smc_vertscale_ruc(k)
            stc_vertscale(k) = stc_vertscale_ruc(k)
          enddo
        endif

        nblks = size(blksz)

        call  set_printing_nb_i(blksz,xlon,xlat,print_i,print_nb)

        do nb =1,nblks
           do i = 1, blksz(nb)

             if ( smc(nb,i,1) .EQ. 1.) cycle ! skip  non-soil (land-ice and non-land)
             ! set printing
             if ( (i==print_i)  .and. (nb==print_nb) ) then
                print_flag = .true.
             else
                print_flag=.false.
             endif

             do v = 1,n_var_lndp
                select case (trim(lndp_var_list(v)))
                !=================================================================
                ! State updates - performed every cycle
                !=================================================================
                case('smc')
                    if (do_pert_state) then
                        p=5.
                        min_bound = smcmin(stype(nb,i))
                        max_bound = smcmax(stype(nb,i))

                      ! with RUC LSM perturb smc only at time step = 2, as in HRRR
                        do k=1,lsoil
                             ! apply perts to a copy of smc, retain original smc
                             ! for later update to liquid soil moisture.
                             ! note: previously we were saving the ice water content 
                             ! (smc-slc) and subtracting this from the perturbed smc to 
                             ! get the perturbed slc. This was introducing small errors in the slc 
                             ! when passing back to the calling program, I think due to precision issues, 
                             ! as the ice content is typically zero. Clara Draper, March, 2022. 
                             tmp_smc = smc(nb,i,k)

                             ! perturb total soil moisture
                             ! factor of sldepth*1000 converts from mm to m3/m3
                             pert = sfc_wts(nb,i,v)*smc_vertscale(k)*lndp_prt_list(v)/(zslayer(k)*1000.)
                             pert = pert*tfactor_state

                             call apply_pert('smc',pert,print_flag, tmp_smc,ierr,p,min_bound, max_bound)

                             ! assign all of applied pert to the liquid soil moisture
                             slc(nb,i,k) = slc(nb,i,k) + tmp_smc - smc(nb,i,k)
                             smc(nb,i,k) = tmp_smc
                             
                        enddo
                    endif

                case('stc')
                    if (do_pert_state) then
                        do k=1,lsoil
                             pert = sfc_wts(nb,i,v)*stc_vertscale(k)*lndp_prt_list(v)
                             pert = tfactor_state
                             call apply_pert('stc',pert,print_flag, stc(nb,i,k),ierr)
                        enddo
                    endif 

                !=================================================================
                ! Parameter updates 
                !=================================================================

                ! are all of the params below included in noah?  

                case('vgf')  ! vegetation fraction
                    if (do_pert_param) then
                         p =5.
                         min_bound=0.
                         max_bound=1.

                         pert = sfc_wts(nb,i,v)*lndp_prt_list(v)
                         pert = pert*tfactor_param
                         call apply_pert ('vfrac',pert,print_flag, vfrac(nb,i), ierr,p,min_bound, max_bound)
                     endif
                case('alb')  ! albedo
                    if (do_pert_param) then
                         p =5.
                         min_bound=0.0
                         max_bound=0.4

                         pert = sfc_wts(nb,i,v)*lndp_prt_list(v)
                         pert = pert*tfactor_param
                         call apply_pert ('alnsf',pert,print_flag, alnsf(nb,i), ierr,p,min_bound, max_bound)
                         call apply_pert ('alnwf',pert,print_flag, alnwf(nb,i), ierr,p,min_bound, max_bound)
                     endif
                case('sal')  ! snow albedo
                    if (do_pert_param) then
                         p =5.
                         min_bound=0.3
                         max_bound=0.85

                         pert = sfc_wts(nb,i,v)*lndp_prt_list(v)
                         pert = pert*tfactor_param
                         call apply_pert ('snoalb',pert,print_flag, snoalb(nb,i), ierr,p,min_bound, max_bound)
                     endif
                case('emi')  ! emissivity
                    if (do_pert_param) then
                         p =5.
                         min_bound=0.8
                         max_bound=1.

                         pert = sfc_wts(nb,i,v)*lndp_prt_list(v)
                         pert = pert*tfactor_param
                         call apply_pert ('semis',pert,print_flag, semis(nb,i), ierr,p,min_bound, max_bound)
                     endif
                case('zol')  ! land roughness length
                    if (do_pert_param) then
                         p =5.
                         min_bound=0.
                         max_bound=300.

                         pert = sfc_wts(nb,i,v)*lndp_prt_list(v)
                         pert = pert*tfactor_param
                         call apply_pert ('zol',pert,print_flag, zorll(nb,i), ierr,p,min_bound, max_bound)
                     endif
                case default
                    print*, &
                     'ERROR: unrecognised lndp_prt_list option in lndp_apply_pert, exiting', trim(lndp_var_list(v))
                    ierr = 10
                    return
                end select
             enddo
           enddo
        enddo
    end subroutine  lndp_apply_perts

!====================================================================
! apply_perts
!====================================================================
! Apply perturbations to selected land states or parameters

  subroutine apply_pert(vname,pert,print_flag, state,ierr,p,vmin, vmax)

   ! intent in
    logical, intent(in)                 :: print_flag
    real(kind=kind_dbl_prec), intent(in)    :: pert
    character(len=*), intent(in)        :: vname ! name of variable being perturbed

    real(kind=kind_dbl_prec), optional, intent(in)    :: p ! flat-top paramater, 0 = no flat-top
                                                       ! flat-top function is used for bounded variables
                                                       ! to reduce the magnitude of perturbations  near boundaries.
    real(kind=kind_dbl_prec), optional, intent(in) :: vmin, vmax ! min,max bounds of variable being perturbed

    ! intent (inout)
    real(kind=kind_dbl_prec), intent(inout) :: state

    ! intent out
    integer                             :: ierr

    !local
    real(kind=kind_dbl_prec) :: z

       if ( print_flag ) then
              write(*,*) 'LNDP - applying lndp to ',vname, ', initial value', state
       endif

       ! apply perturbation
       if (present(p) ) then
           if ( .not. (present(vmin) .and. present(vmax) )) then
              ierr=20
              print*, 'error, flat-top function requires min & max to be specified'
           endif

           z = -1. + 2*(state  - vmin)/(vmax - vmin) ! flat-top function
           state =  state  + pert*(1-abs(z**p))
       else
           state =  state  + pert
       endif

       if (present(vmax)) state =  min( state , vmax )
       if (present(vmin)) state =  max( state , vmin )

       if ( print_flag ) then
              write(*,*) 'LNDP - applying lndp to ',vname, ', final value', state
       endif

  end subroutine apply_pert


!====================================================================
! set_printing_nb_i
!====================================================================
! routine to turn on print flag for selected location
!
    subroutine set_printing_nb_i(blksz,xlon,xlat,print_i,print_nb)

        implicit none

        ! intent (in)
        integer,                  intent(in) :: blksz(:)
        real(kind=kind_dbl_prec),     intent(in) :: xlon(:,:)
        real(kind=kind_dbl_prec),     intent(in) :: xlat(:,:)


        ! intent (out)
        integer, intent(out) :: print_i, print_nb

        ! local
        integer :: nblks,nb,i
        real, parameter :: plon_trunc =  114.9
        real, parameter :: plat_trunc =  -26.6
        real, parameter  :: delta  = 1.

        nblks = size(blksz)

        print_i = -9
        print_nb = -9
        do nb = 1,nblks
         do i = 1,blksz(nb)
        if ( (xlon(nb,i)*57.29578 > plon_trunc) .and.  (xlon(nb,i)*57.29578 < plon_trunc+delta ) .and. &
           (xlat(nb,i)*57.29578 >  plat_trunc ) .and.  (xlat(nb,i)*57.29578 < plat_trunc+delta ) ) then
                      print_i=i
                      print_nb=nb
                      write(*,*) 'LNDP -print flag is on', xlon(nb,i)*57.29578, xlat(nb,i)*57.29578, nb, i
                      return
         endif
         enddo
        enddo

    end subroutine set_printing_nb_i

end module lndp_apply_perts_mod


