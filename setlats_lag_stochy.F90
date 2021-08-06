!>@brief The module 'setlats_lag_stochy_mod' contains the subroutine setlats_lag_stochy
      module setlats_lag_stochy_mod

      implicit none

      contains

!>@brief The subroutine 'setlats_a_stochy' selects the latitude points on this task
! and halos
!>@details This code is taken from the legacy spectral GFS
      subroutine setlats_lag_stochy(gis_stochy, lats_nodes_h, global_lats_h)
!
      use spectral_layout_mod,   only : latg
      use stochy_internal_state_mod, only : stochy_internal_state

      implicit none
      type(stochy_internal_state), intent(inout) :: gis_stochy
!
!
      integer            lats_nodes_a(gis_stochy%nodes), lats_nodes_h(gis_stochy%nodes), &
                         global_lats_a(latg),                      &
                         global_lats_h(latg+2*gis_stochy%yhalo*gis_stochy%nodes)
!
      integer              jj,jpt_a,jpt_h,lat_val,nn,nodes_lats,   &
                           j1, j2, iprint
!
      lats_nodes_h = 0
!
      nodes_lats = 0
      do nn=1,gis_stochy%nodes
         if (gis_stochy%lats_nodes_a(nn) > 0) then
             lats_nodes_h(nn) = gis_stochy%lats_nodes_a(nn) + gis_stochy%yhalo + gis_stochy%yhalo
             nodes_lats       = nodes_lats + 1
         endif
      enddo
!
      global_lats_h = 0
!
!    set non-yhalo latitudes
!
      jpt_a = 0
      jpt_h = gis_stochy%yhalo
      do nn=1,gis_stochy%nodes
         if (gis_stochy%lats_nodes_a(nn) > 0) then
            do jj=1,gis_stochy%lats_nodes_a(nn)
               jpt_a = jpt_a + 1
               jpt_h = jpt_h + 1
               global_lats_h(jpt_h) = gis_stochy%global_lats_a(jpt_a)
            enddo
            jpt_h = jpt_h + gis_stochy%yhalo + gis_stochy%yhalo
         endif
      enddo
!
      j1 = latg + (gis_stochy%yhalo+gis_stochy%yhalo) * nodes_lats
      do jj=1,gis_stochy%yhalo
        j2 = gis_stochy%yhalo - jj
         global_lats_h(jj)    = gis_stochy%global_lats_a(1)    + j2     ! set north pole yhalo
         global_lats_h(j1-j2) = gis_stochy%global_lats_a(latg) + 1 - jj ! set south pole yhalo
      enddo
!
      if (gis_stochy%lats_nodes_a(1) /= latg) then
!
!       set non-polar south yhalos
         jpt_h = 0
         do nn=1,gis_stochy%nodes-1
            if (lats_nodes_h(nn).GT.0) then
               jpt_h   = jpt_h + lats_nodes_h(nn)
               lat_val = global_lats_h(jpt_h-gis_stochy%yhalo)
               do jj=1,gis_stochy%yhalo
                  global_lats_h(jpt_h-gis_stochy%yhalo+jj) = min(lat_val+jj,latg)
               enddo
            endif
         enddo
!
!       set non-polar north yhalos
         jpt_h = 0
         do nn=1,gis_stochy%nodes-1
            if (lats_nodes_h(nn).GT.0) then
               jpt_h   = jpt_h + lats_nodes_h(nn)
               lat_val = global_lats_h(jpt_h+gis_stochy%yhalo+1)
               do jj=1,gis_stochy%yhalo
                  global_lats_h(jpt_h+gis_stochy%yhalo-(jj-1)) = max(lat_val-jj,1)
               enddo
            endif
         enddo
!
      endif
!

      return
      end

      end module setlats_lag_stochy_mod
