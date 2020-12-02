!>@brief The module 'four_to_grid_mod' contains the subroute four_to_grid
      module four_to_grid_mod


      implicit none

      contains

!>@brief The subroutine 'four_to_grd' calculate real values form fourrier coefficients
!>@details This code is taken from the legacy spectral GFS
      subroutine four_to_grid(syn_gr_a_1,syn_gr_a_2,
     &                        lon_dim_coef,lon_dim_grid,lons_lat,lot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use machine
      implicit none
!!
      real(kind=kind_dbl_prec)     syn_gr_a_1(lon_dim_coef,lot)
      real(kind=kind_dbl_prec)     syn_gr_a_2(lon_dim_grid,lot)
      integer                  lon_dim_coef
      integer                  lon_dim_grid
      integer                  lons_lat
      integer                  lot
!________________________________________________________
#ifdef MKL
      integer*8                plan
#else
      real(kind=kind_dbl_prec) aux1crs(42002)
      real(kind=kind_dbl_prec)     scale_ibm
      integer                  ibmsign
      integer                  init
#endif
#ifdef MKL
      include                  "fftw/fftw3.f"
      integer                  NULL
#else
      external                 dcrft
      external                 scrft
#endif
!________________________________________________________


#ifdef MKL
             call dfftw_plan_many_dft_c2r(                              &
     &                    plan, 1, lons_lat, lot,                       &
     &                    syn_gr_a_1, NULL, 1, size(syn_gr_a_1,dim=1),  &
     &                    syn_gr_a_2, NULL, 1, size(syn_gr_a_2,dim=1),  &
     &                    FFTW_ESTIMATE)
             call dfftw_execute(plan)
             call dfftw_destroy_plan(plan)
#else
             init      = 1
             ibmsign   = -1
             scale_ibm = 1.0d0

             call dcrft(init,
     &               syn_gr_a_1(1,1)   ,lon_dim_coef/2,
     &               syn_gr_a_2(1,1)   ,lon_dim_grid,
     &               lons_lat,lot,ibmsign,scale_ibm,
     &               aux1crs,22000,
     &               aux1crs(22001),20000)

             init = 0
             call dcrft(init,
     &               syn_gr_a_1(1,1)   ,lon_dim_coef/2,
     &               syn_gr_a_2(1,1)   ,lon_dim_grid,
     &               lons_lat,lot,ibmsign,scale_ibm,
     &               aux1crs,22000,
     &               aux1crs(22001),20000)
#endif

      return
      end
      end module four_to_grid_mod
