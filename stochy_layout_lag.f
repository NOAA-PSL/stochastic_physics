!>@brief The module 'stochy_layout_lag' contains the decomposition attributes of the gaussian grid
      module stochy_layout_lag
      use kinddef
      implicit none
      save
cc
      integer lats_dim_h,
     x        lats_node_h,
     x        lats_node_h_max,
     x        ipt_lats_node_h,
     x        lon_dim_h
cc
      INTEGER ,ALLOCATABLE :: lat1s_h(:)
      end module stochy_layout_lag
