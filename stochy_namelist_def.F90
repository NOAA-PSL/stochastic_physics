      module stochy_namelist_def
!
! program log
! 11 Oct 2016:    Philip Pegion create standalone stochastic physics
!
      use machine
      implicit none

      public
      integer nssppt,nsshum,nsskeb,lon_s,lat_s,ntrunc

! pjp stochastic phyics
      integer skeb_varspect_opt,skeb_npass
      logical sppt_sfclimit

      real(kind=kind_dbl_prec) :: skeb_sigtop1,skeb_sigtop2,          &
                         sppt_sigtop1,sppt_sigtop2,shum_sigefold, &
                         skeb_vdof
      real(kind=kind_dbl_prec) fhstoch,skeb_diss_smooth,spptint,shumint,skebint,skebnorm
      real(kind=kind_dbl_prec), dimension(5) :: skeb,skeb_lscale,skeb_tau
      real(kind=kind_dbl_prec), dimension(5) :: sppt,sppt_lscale,sppt_tau
      real(kind=kind_dbl_prec), dimension(5) :: shum,shum_lscale,shum_tau
      integer,dimension(5) ::skeb_vfilt
      integer(8),dimension(5) ::iseed_sppt,iseed_shum,iseed_skeb
      logical stochini,sppt_logit,new_lscale
      logical use_zmtnblck
      logical do_shum,do_sppt,do_skeb

! mg surface perturbations
      real(kind=kind_dbl_prec), dimension(5) :: lndp_lscale,lndp_tau
      real(kind=kind_dbl_prec), dimension(5) :: lndp_z0,lndp_hc,lndp_zt
      real(kind=kind_dbl_prec), dimension(5) :: lndp_la,lndp_vf,lndp_al
      integer n_var_lndp
      integer(8),dimension(5) ::iseed_lndp
      integer lndp_type
      integer lndp_ind_z0,lndp_ind_hc,lndp_ind_zt
      integer lndp_ind_la,lndp_ind_vf,lndp_ind_al

      end module stochy_namelist_def
