!>@brief The module 'get_ls_node_stochy_mod' contains the subroutine get_ls_node_stochy
      module get_ls_node_stochy_mod

      implicit none

      contains

!>@brief The subroutine 'get_ls_node_stochy' calculates the decomposition of the spherical harmonics based on the processor layout
      subroutine get_ls_node_stochy(me_fake,ls_node,ls_max_node_fake, nodes)
!>@details This code is taken from the legacy spectral GFS
!
      use spectral_layout_mod, only : jcap,ls_dim
      implicit none
!
      integer   me_fake, ls_max_node_fake, nodes
      integer   ls_node(ls_dim)

      integer   ijk, jptls, l, node, nodesio, jcap1
!
      nodesio = nodes

      ls_node = -1
      jcap1=jcap+1
      print*,'in ls_node',jcap1,nodes,me_fake
!
      jptls =  0
      l = 0
!.............................................
      do ijk=1,jcap1
!
         do node=1,nodesio
            if (node == me_fake+1) then
               jptls = jptls + 1
               ls_node(jptls) = l
            endif
            l = l + 1
            if (l > jcap) go to 200
         enddo
!
         do node=nodesio,1,-1
            if (node == me_fake+1) then
               jptls = jptls + 1
               ls_node(jptls) = l
            endif
            l = l + 1
            if (l > jcap) go to 200
         enddo
!
      enddo
!.............................................
!
  200 continue
      print*,'in ls_node',minval(ls_node),maxval(ls_node)
!
!.............................................
!
      ls_max_node_fake = 0
      do ijk=1,ls_dim
         if(ls_node(ijk) >= 0) then
            ls_max_node_fake = ijk
          endif
      enddo
!
      return
      end

      end module get_ls_node_stochy_mod
