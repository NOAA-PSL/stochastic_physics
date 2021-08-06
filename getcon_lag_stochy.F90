!>@brief The module 'getcon_lag_stochy_mod' contains the subroute getcon_lag_stochy
      module getcon_lag_stochy_mod

      implicit none

      contains

!>@brief The subroutine 'getcon_lag' calculates grid properties and domain decompostion for the gaussian grid
!>@details This code is taken from the legacy spectral GFS
      subroutine getcon_lag_stochy(gis_stochy, lats_nodes_h,global_lats_h_sn)

      use spectral_layout_mod, only : jcap,latg,latg2,lonf, &
                                colrad_a,sinlat_a,                   &
                                ipt_lats_node_h,lats_dim_h,          &
                                lats_node_h,lats_node_h_max
      use setlats_lag_stochy_mod, only: setlats_lag_stochy
      use stochy_internal_state_mod, only : stochy_internal_state
      implicit none

      type(stochy_internal_state), intent(inout) :: gis_stochy
      integer, intent(inout), dimension(gis_stochy%nodes) :: lats_nodes_h
      integer, intent(inout), dimension(latg+2*gis_stochy%yhalo*gis_stochy%nodes) :: global_lats_h_sn
!
      integer  i,j,l,n,lat,i1,i2,node,nodesio
      integer, dimension(latg+2*gis_stochy%yhalo*gis_stochy%nodes) :: global_lats_h_ns
!
      do lat = 1, latg2
         gis_stochy%lonsperlat(latg+1-lat) = gis_stochy%lonsperlat(lat)
      end do
      nodesio = gis_stochy%nodes

      call setlats_lag_stochy(gis_stochy,lats_nodes_h,global_lats_h_ns)

!  reverse order for use in set_halos

      i1 = 1
      i2 = 0
      do n=1,gis_stochy%nodes
         j  = 0
         i2 = i2 + lats_nodes_h(n)
         do i=i1,i2
            j = j + 1
            global_lats_h_sn(i) = global_lats_h_ns(i2+1-j)
         enddo
         i1 = i2 + 1
      enddo

 830   format(10(i4,1x))
      lats_dim_h = 0
      do node=1,gis_stochy%nodes
         lats_dim_h = max(lats_dim_h, lats_nodes_h(node))
      enddo
      lats_node_h     = lats_nodes_h(gis_stochy%mype+1)
      lats_node_h_max = 0
      do i=1,gis_stochy%nodes
        lats_node_h_max  = max(lats_node_h_max, lats_nodes_h(i))
      enddo
      ipt_lats_node_h = 1
      if ( gis_stochy%mype > 0 ) then
         do node=1,gis_stochy%mype
            ipt_lats_node_h = ipt_lats_node_h + lats_nodes_h(node)
         enddo
      endif
      do j=1,latg2
        sinlat_a(j) = cos(colrad_a(j))
      enddo
      return
      end

      end module getcon_lag_stochy_mod
