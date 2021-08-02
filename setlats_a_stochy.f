!>@brief The module 'setlats_a_stochy_mod' contains the subroutine setlats_a_stochy
      module setlats_a_stochy_mod

      implicit none

      contains
!>@brief The subroutine 'setlats_a_stochy' selects the latitude points on this task
!>@details This code is taken from the legacy spectral GFS
      subroutine setlats_a_stochy(lats_nodes_a,global_lats_a,
     &                         iprint,lonsperlat)
!
      use spectral_layout_mod   , only : nodes,me,latg,lonf
!
      implicit none
!
      integer, dimension(latg) :: global_lats_a, lonsperlat
      integer                     lats_nodes_a(nodes)

      integer              iprint,ifin,nodesio
     &,                    jcount,jpt,lat,lats_sum,node,i,ii
     &,                    ngrptg,ngrptl,ipe,irest,idp
     &,                    ngrptgh,nodesioh
!    &,                    ilatpe,ngrptg,ngrptl,ipe,irest,idp
!
      integer,allocatable :: lats_hold(:,:)
!
      allocate ( lats_hold(latg,nodes) )
!
      iprint = 0
      opt    = 1                       ! reduced grid
      lats_nodes_a = 0
        nodesio = nodes
!
      ngrptg = 0
      do lat=1,latg
         do i=1,lonsperlat(lat)
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
          ifin   = lonsperlat(lat)
          ngrptl = ngrptl + ifin

          if (ngrptl*nodesio <= ngrptg+irest) then
            lats_nodes_a(ipe+1)  = lats_nodes_a(ipe+1) + 1
            lats_hold(idp,ipe+1) = lat
            idp = idp + 1
          else
            ipe = ipe + 1
            if (ipe <= nodesio) lats_hold(1,ipe+1) = lat
            idp    = 2
            irest  = irest + ngrptg - (ngrptl-ifin)*nodesio
            ngrptl = ifin
            lats_nodes_a(ipe+1) = lats_nodes_a(ipe+1) + 1
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
          ifin   = lonsperlat(lat)
          ngrptl = ngrptl + ifin

          if (ngrptl*nodesioh <= ngrptgh+irest .or. lat == latg/2) then
            lats_nodes_a(ipe+1)  = lats_nodes_a(ipe+1) + 1
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
            lats_nodes_a(ipe+1) = lats_nodes_a(ipe+1) + 1
          endif
        enddo
        do node=1, nodesioh
          ii = nodesio-node+1
          jpt = lats_nodes_a(node)
          lats_nodes_a(ii) = jpt
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
        if ( lats_nodes_a(node) > 0 ) then
          do jcount=1,lats_nodes_a(node)
            global_lats_a(jpt+jcount) = lats_hold(jcount,node)
          enddo
        endif
        jpt = jpt + lats_nodes_a(node)
      enddo
!!
      deallocate (lats_hold)
      if ( iprint /= 1 ) return
!!
      if (me == 0) then
      jpt=0
      do node=1,nodesio
         if ( lats_nodes_a(node) > 0 ) then
            print 600
            lats_sum=0
            do jcount=1,lats_nodes_a(node)
               lats_sum=lats_sum + lonsperlat(global_lats_a(jpt+jcount))
               print 700, node-1,
     x                    node,    lats_nodes_a(node),
     x                    jpt+jcount, global_lats_a(jpt+jcount),
     x                     lonsperlat(global_lats_a(jpt+jcount)),
     x                    lats_sum
            enddo
         endif
         jpt=jpt+lats_nodes_a(node)
      enddo
!
      print 600
!
  600 format ( ' ' )
!
  700 format (  'setlats  me=', i4,
     x          '  lats_nodes_a(',  i4, ' )=', i4,
     x          '  global_lats_a(', i4, ' )=', i4,
     x          '  lonsperlat=', i5,
     x          '  lats_sum=',   i6 )
!
      endif

      return
      end

      end module setlats_a_stochy_mod
