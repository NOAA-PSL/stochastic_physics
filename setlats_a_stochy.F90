!>@brief The module 'setlats_a_stochy_mod' contains the subroutine setlats_a_stochy
      module setlats_a_stochy_mod

      implicit none

      contains
!>@brief The subroutine 'setlats_a_stochy' selects the latitude points on this task
!>@details This code is taken from the legacy spectral GFS
      subroutine setlats_a_stochy(gis_stochy)
!
      use spectral_layout_mod   , only : latg,lonf
      use stochy_internal_state_mod, only : stochy_internal_state
!
      implicit none
!
      type(stochy_internal_state), intent(inout) :: gis_stochy

      integer :: ifin,nodesio,                       &
                 jcount,jpt,lat,lats_sum,node,i,ii,  &
                 ngrptg,ngrptl,ipe,irest,idp,        &
                 ngrptgh,nodesioh           
!
      integer,allocatable :: lats_hold(:,:)
!
      allocate ( lats_hold(latg,gis_stochy%nodes) )
!
      gis_stochy%lats_nodes_a = 0
      nodesio = gis_stochy%nodes
!
      ngrptg = 0
      do lat=1,latg
         do i=1,gis_stochy%lonsperlat(lat)
           ngrptg = ngrptg + 1
         enddo
      enddo

!
!   ngrptg contains total number of grid points.
!
!     distribution of the grid
      nodesioh = nodesio / 2

      if (nodesioh*2 /= nodesio) then
        ngrptl = 0
        ipe    = 0
        irest  = 0
        idp    = 1

        do lat=1,latg
          ifin   = gis_stochy%lonsperlat(lat)
          ngrptl = ngrptl + ifin

          if (ngrptl*nodesio <= ngrptg+irest) then
            gis_stochy%lats_nodes_a(ipe+1)  = gis_stochy%lats_nodes_a(ipe+1) + 1
            lats_hold(idp,ipe+1) = lat
            idp = idp + 1
          else
            ipe = ipe + 1
            if (ipe <= nodesio) lats_hold(1,ipe+1) = lat
            idp    = 2
            irest  = irest + ngrptg - (ngrptl-ifin)*nodesio
            ngrptl = ifin
            gis_stochy%lats_nodes_a(ipe+1) = gis_stochy%lats_nodes_a(ipe+1) + 1
          endif
        enddo
      else
        nodesioh = nodesio/2
        ngrptgh  = ngrptg/2
        ngrptl = 0
        ipe    = 0
        irest  = 0
        idp    = 1

        do lat=1,latg/2
          ifin   = gis_stochy%lonsperlat(lat)
          ngrptl = ngrptl + ifin

          if (ngrptl*nodesioh <= ngrptgh+irest .or. lat == latg/2) then
            gis_stochy%lats_nodes_a(ipe+1)  = gis_stochy%lats_nodes_a(ipe+1) + 1
            lats_hold(idp,ipe+1) = lat
            idp = idp + 1
          else
            ipe = ipe + 1
            if (ipe <= nodesioh) then
              lats_hold(1,ipe+1) = lat
            endif
            idp    = 2
            irest  = irest + ngrptgh - (ngrptl-ifin)*nodesioh
            ngrptl = ifin
            gis_stochy%lats_nodes_a(ipe+1) = gis_stochy%lats_nodes_a(ipe+1) + 1
          endif
        enddo
        do node=1, nodesioh
          ii = nodesio-node+1
          jpt = gis_stochy%lats_nodes_a(node)
          gis_stochy%lats_nodes_a(ii) = jpt
          do i=1,jpt
            lats_hold(jpt+1-i,ii) = latg+1-lats_hold(i,node)
          enddo
        enddo


      endif
!!
!!........................................................
!!
      jpt = 0
      do node=1,nodesio
        if ( gis_stochy%lats_nodes_a(node) > 0 ) then
          do jcount=1,gis_stochy%lats_nodes_a(node)
            gis_stochy%global_lats_a(jpt+jcount) = lats_hold(jcount,node)
          enddo
        endif
        jpt = jpt + gis_stochy%lats_nodes_a(node)
      enddo

      deallocate (lats_hold)

      return
      end

      end module setlats_a_stochy_mod
