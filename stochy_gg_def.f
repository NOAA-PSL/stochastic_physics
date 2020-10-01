!>@brief The module 'stochy_gg_def' declares array defining the gaussian grid attributes
      module stochy_gg_def
      use kinddef
      implicit none

      real(kind=kind_dbl_prec), allocatable, dimension(:) ::  colrad_a,
     &                      wgt_a, wgtcs_a, rcs2_a, sinlat_a, coslat_a
!
      integer ,allocatable, dimension(:) :: lats_nodes_h,global_lats_h
      end module stochy_gg_def
