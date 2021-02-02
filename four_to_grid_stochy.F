!>@brief The module 'four_to_grid_mod' contains the subroute four_to_grid
      module four_to_grid_mod


      implicit none

      contains

!>@brief The subroutine 'four_to_grd' calculate real values form fourrier coefficients
!>@details This code is taken from the legacy spectral GFS
      subroutine four_to_grid(syn_gr_a_1,syn_gr_a_2,
     &                        lon_dim_coef,lon_dim_grid,lons_lat,lot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use kinddef
      implicit none
!!
      real(kind=kind_dbl_prec)     syn_gr_a_1(lon_dim_coef,lot)
      real(kind=kind_dbl_prec)     syn_gr_a_2(lon_dim_grid,lot)
      integer                  lon_dim_coef
      integer                  lon_dim_grid
      integer                  lons_lat
      integer                  lot
!________________________________________________________
      real(kind=kind_dbl_prec) aux1crs(44002)
      real(kind=kind_dbl_prec)     scale_ibm
      integer                  ibmsign
      integer                  init
!________________________________________________________


             init      = 1
             ibmsign   = -1
             scale_ibm = 1.0d0

             call dcrft_stochy(init,
     &               syn_gr_a_1(1,1)   ,lon_dim_coef/2,
     &               syn_gr_a_2(1,1)   ,lon_dim_grid,
     &               lons_lat,lot,ibmsign,scale_ibm,
     &               aux1crs,22000,
     &               aux1crs(22001),20000)

             init = 0
             call dcrft_stochy(init,
     &               syn_gr_a_1(1,1)   ,lon_dim_coef/2,
     &               syn_gr_a_2(1,1)   ,lon_dim_grid,
     &               lons_lat,lot,ibmsign,scale_ibm,
     &               aux1crs,22000,
     &               aux1crs(22001),20000)

      return
      end


      end module four_to_grid_mod
