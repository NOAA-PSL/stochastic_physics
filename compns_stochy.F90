module compns_stochy_mod

   implicit none

   contains

!-----------------------------------------------------------------------
      subroutine compns_stochy (me,sz_nml,input_nml_file,fn_nml,nlunit,deltim,iret)
!$$$  Subprogram Documentation Block
!
! Subprogram:  compns     Check and compute namelist frequencies
!   Prgmmr: Iredell       Org: NP23          Date: 1999-01-26
!
! Abstract: This subprogram checks global spectral model namelist
!           frequencies in hour units for validity.  If they are valid,
!           then the frequencies are computed in timestep units.
!           The following rules are applied:
!             1. the timestep must be positive;
!
! Program History Log:
!   2016-10-11  Phil Pegion  make the stochastic physics stand alone
!
! Usage:    call compns_stochy (me,deltim,nlunit, stochy_namelist,iret)
!   Input Arguments:
!     deltim   - real timestep in seconds
!   Output Arguments:
!     iret     - integer return code (0 if successful or
!                between 1 and 8 for which rule above was broken)
!     stochy_namelist
!
! Attributes:
!   Language: Fortran 90
!
!$$$


      use stochy_namelist_def

      implicit none


      integer,              intent(out) :: iret
      integer,              intent(in)  :: nlunit,me,sz_nml
      character(len=*),     intent(in)  :: input_nml_file(sz_nml)
      character(len=64),    intent(in)  :: fn_nml
      real,                 intent(in)  :: deltim
      real tol,l_min
      real :: rerth,circ,tmp_lat
      integer k,ios
      integer,parameter :: four=4

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      namelist /nam_stochy/ntrunc,lon_s,lat_s,sppt,sppt_tau,sppt_lscale,sppt_logit, &
      iseed_shum,iseed_sppt,shum,shum_tau,&
      shum_lscale,fhstoch,stochini,skeb_varspect_opt,sppt_sfclimit, &
      skeb,skeb_tau,skeb_vdof,skeb_lscale,iseed_skeb,skeb_vfilt,skeb_diss_smooth, &
      skeb_sigtop1,skeb_sigtop2,skebnorm,sppt_sigtop1,sppt_sigtop2,&
      shum_sigefold,spptint,shumint,skebint,skeb_npass,use_zmtnblck,new_lscale
      namelist /nam_sfcperts/lndp_type,lndp_z0,lndp_hc,lndp_zt,lndp_la, & ! mg, sfcperts
      lndp_vf,lndp_al,iseed_lndp,lndp_tau,lndp_lscale 

      rerth  =6.3712e+6      ! radius of earth (m)
      tol=0.01  ! tolerance for calculations
!     spectral resolution defintion
      ntrunc=-999
      lon_s=-999
      lat_s=-999
      ! can specify up to 5 values for the stochastic physics parameters
      ! (each is an array of length 5)
      sppt             = -999.  ! stochastic physics tendency amplitude
      shum             = -999.  ! stochastic boundary layer spf hum amp
      skeb             = -999.  ! stochastic KE backscatter amplitude
      ! mg, sfcperts
      lndp_z0          = -999.  ! momentum roughness length amplitude
      lndp_hc          = -999.  ! soil hydraulic conductivity amp
      lndp_zt          = -999.  ! mom/heat roughness length amplitude
      lndp_la          = -999.  ! leaf area index amplitude
      lndp_vf          = -999.  ! vegetation fraction amplitude
      lndp_al          = -999.  ! albedo perturbations amplitude
! logicals
      do_sppt = .false.
      use_zmtnblck = .false.
      new_lscale = .false.
      do_shum = .false.
      do_skeb = .false.
! input land pert variables: 
! LNDP_TYPE = 0
! no explicit land perturbations
! LNDP_Type =  1 
! this is the initial land sfc pert scheme, introduced and tested for impact on GEFS forecasts.
! see https://journals.ametsoc.org/doi/full/10.1175/MWR-D-18-0057.1
! perturbations are assigned once at the start of the forecast
      lndp_type = 0 ! 0 -none, 1 - lndp_fcst, 2  - lndp_edas
      lndp_lscale  = -999.       ! length scales
      lndp_tau     = -999.       ! time scales
      iseed_lndp   = 0           ! random seeds (if 0 use system clock)
! derived land pert variables
      lndp_ind_z0 = 0
      lndp_ind_hc = 0
      lndp_ind_zt = 0
      lndp_ind_la = 0
      lndp_ind_vf = 0
      lndp_ind_al = 0
! for SKEB random patterns.
      skeb_vfilt       = 0
      skebint          = 0
      spptint          = 0
      shumint          = 0
      skeb_npass       = 11  ! number of passes of smoother for dissipation estiamte
      sppt_tau         = -999.  ! time scales
      shum_tau         = -999.
      skeb_tau         = -999.
      skeb_vdof        = 5 ! proxy for vertical correlation, 5 is close to 40 passes of the 1-2-1 filter in the GFS
      skebnorm         = 0  ! 0 - random pattern is stream function, 1- pattern is kenorm, 2- pattern is vorticity
      sppt_lscale      = -999.  ! length scales
      shum_lscale      = -999.
      skeb_lscale      = -999.
      iseed_sppt       = 0      ! random seeds (if 0 use system clock)
      iseed_shum       = 0
      iseed_skeb       = 0
! parameters to control vertical tapering of stochastic physics with
! height
      sppt_sigtop1 = 0.1
      sppt_sigtop2 = 0.025
      skeb_sigtop1 = 0.1
      skeb_sigtop2 = 0.025
      shum_sigefold = 0.2
! reduce amplitude of sppt near surface (lowest 2 levels)
      sppt_sfclimit = .false.
! gaussian or power law variance spectrum for skeb (0: gaussian, 1:
! power law). If power law, skeb_lscale interpreted as a power not a
! length scale.
      skeb_varspect_opt = 0
      sppt_logit        = .false. ! logit transform for sppt to bounded interval [-1,+1]
      fhstoch           = -999.0  ! forecast interval (in hours) to dump random patterns
      stochini          = .false. ! true= read in pattern, false=initialize from seed

#ifdef INTERNAL_FILE_NML
      read(input_nml_file, nml=nam_stochy)
#else
      rewind (nlunit)
      open (unit=nlunit, file=fn_nml, READONLY, status='OLD', iostat=ios)
      read(nlunit,nam_stochy)
#endif
#ifdef INTERNAL_FILE_NML
      read(input_nml_file, nml=nam_sfcperts)
#else
      rewind (nlunit)
      open (unit=nlunit, file=fn_nml, READONLY, status='OLD', iostat=ios)
      read(nlunit,nam_sfcperts)
#endif

      if (me == 0) then
      print *,' in compns_stochy'
      print*,'skeb=',skeb
      endif

! PJP stochastic physics additions
      IF (sppt(1) > 0 ) THEN
        do_sppt=.true.
      ENDIF
      IF (shum(1) > 0 ) THEN
        do_shum=.true.
!     shum parameter has units of 1/hour, to remove time step
!     dependence.
!     change shum parameter units from per hour to per timestep
         DO k=1,5
            IF (shum(k) .gt. 0.0) shum(k)=shum(k)*deltim/3600.0
         ENDDO
      ENDIF
      IF (skeb(1) > 0 ) THEN
         do_skeb=.true.
         if (skebnorm==0) then ! stream function norm
            skeb=skeb*1.111e3*sqrt(deltim)
            !skeb=skeb*5.0e5/sqrt(deltim)
         endif
         if (skebnorm==1) then ! stream function norm
            skeb=skeb*0.00222*sqrt(deltim)
            !skeb=skeb*1/sqrt(deltim)
         endif
         if (skebnorm==2) then ! vorticty function norm
            skeb=skeb*1.111e-9*sqrt(deltim)
            !skeb=skeb*5.0e-7/sqrt(deltim)
         endif
!      adjust skeb values for resolution.
!      scaling is such that a value of 1.0 at T574 with a 900 second
!      timestep produces well-calibrated values of forecast spread.
!         DO k=1,5
!            IF  (skeb(k) .gt. 0.0) THEN
!               skeb(k)=skeb(k)*deltim/(ntrunc*(ntrunc+1))*365765.0  ! 365765 is a scale factor so the base SKEB value in the namelist is 1.0
!               skeb(k)=skeb(k)*deltim/(ntrunc*(ntrunc+1))*2000.0  ! 2000 is new scale factor so the base SKEB value in the namelist is 1.0
!            ENDIF
!         ENDDO
      ENDIF
!    compute frequencty to estimate dissipation timescale
      IF (skebint == 0.) skebint=deltim
      nsskeb=nint(skebint/deltim)                              ! skebint in seconds
      IF(nsskeb<=0 .or. abs(nsskeb-skebint/deltim)>tol) THEN
         WRITE(0,*) "SKEB interval is invalid",skebint
        iret=9
        return
      ENDIF
      IF (spptint == 0.) spptint=deltim
      nssppt=nint(spptint/deltim)                              ! spptint in seconds
      IF(nssppt<=0 .or. abs(nssppt-spptint/deltim)>tol) THEN
         WRITE(0,*) "SPPT interval is invalid",spptint
        iret=9
        return
      ENDIF
      IF (shumint == 0.) shumint=deltim
      nsshum=nint(shumint/deltim)                              ! shumint in seconds
      IF(nsshum<=0 .or. abs(nsshum-shumint/deltim)>tol) THEN
         WRITE(0,*) "SHUM interval is invalid",shumint
        iret=9
        return
      ENDIF
!calculate ntrunc if not supplied
     if (ntrunc .LT. 1) then  
        if (me==0) print*,'ntrunc not supplied, calculating'
        circ=2*3.1415928*rerth ! start with lengthscale that is circumference of the earth
        l_min=circ
        do k=1,5
           if (sppt(k).GT.0) l_min=min(sppt_lscale(k),l_min)
           if (shum(k).GT.0) l_min=min(shum_lscale(k),l_min)
           if (skeb(k).GT.0) l_min=min(skeb_lscale(k),l_min)
       enddo
       !ntrunc=1.5*circ/l_min
       ntrunc=circ/l_min
       if (me==0) print*,'ntrunc calculated from l_min',l_min,ntrunc
     endif
     ! ensure lat_s is a mutiple of 4 with a reminader of two
     ntrunc=INT((ntrunc+1)/four)*four+2
     if (me==0) print*,'NOTE ntrunc adjusted for even nlats',ntrunc

! set up gaussian grid for ntrunc if not already defined. 
     if (lon_s.LT.1 .OR. lat_s.LT.1) then
        lat_s=ntrunc*1.5+1
        lon_s=lat_s*2+4
! Grid needs to be larger since interpolation is bi-linear
        lat_s=lat_s*2
        lon_s=lon_s*2
        if (me==0) print*,'gaussian grid not set, defining here',lon_s,lat_s
     endif

! 
! land perts  - parse nml input
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     select case (lndp_type)
     case (0) 
        if (me==0) print*, & 
           'no land perturbations selected'
     case (1) 
        if (me==0) print*, & 
          'land perturbations will be applied to selected paramaters, using older scheme designed for S2S fcst spread' 
          ! set indexes for perturbed params  
          if (any(lndp_z0 > 0)  ) then
            n_var_lndp = n_var_lndp+1
            lndp_ind_z0 =  n_var_lndp
          endif
          if (any(lndp_zt  > 0) ) then
            n_var_lndp = n_var_lndp+1
            lndp_ind_zt =  n_var_lndp
          endif
          if (any(lndp_hc > 0) ) then
            n_var_lndp = n_var_lndp+1
            lndp_ind_hc =  n_var_lndp
          endif 
          if ( any(lndp_la > 0) ) then
            n_var_lndp = n_var_lndp+1
            lndp_ind_la =  n_var_lndp
          endif
          if ( any(lndp_al > 0) ) then
            n_var_lndp = n_var_lndp+1
            lndp_ind_al =  n_var_lndp
          endif
          if ( any(lndp_vf > 0) ) then
            n_var_lndp = n_var_lndp+1
            lndp_ind_vf =  n_var_lndp
          endif
     case (2) 
        if (me==0) print*, & 
         'land perturbations will be applied to selected paramaters, using newer scheme designed for DA ens spread'
     case default 
        if (me==0) print*, & 
         'lndp_type out of range, set to 0 (none), 1 (for fcst spread), 2 (for cycling DA spread)'
         iret = 10 
         return 
     end select 
!
!  All checks are successful.
!
      if (me == 0) then
         print *, 'stochastic physics'
         print *, ' do_sppt : ', do_sppt
         print *, ' do_shum : ', do_shum
         print *, ' do_skeb : ', do_skeb
         print *, ' lndp_type : ', lndp_type
         if (lndp_type .NE. 0) print *, ' n_var_lndp : ', n_var_lndp
      endif
      iret = 0
!
      return
      end subroutine compns_stochy

end module compns_stochy_mod
