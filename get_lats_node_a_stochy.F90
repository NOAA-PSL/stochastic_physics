!>@brief The module 'get_lats_node_a_stochy_mod' contains the subroutine get_lats_node_a_stochy
      module get_lats_node_a_stochy_mod

      implicit none

      contains
!>@brief The subroutine 'get_lats_node_a_stochy' calculates the decomposition of the gaussian grid based on the processor layout
!>@details This code is taken from the legacy spectral GFS
      subroutine get_lats_node_a_stochy(me_fake,global_lats_a, &
                lats_nodes_a_fake,gl_lats_index,               &
                global_time_sort_index,nodes)
!
      use spectral_layout_mod
      implicit none

      integer,intent(in)  :: me_fake
      integer,intent(in)  :: nodes
      integer,intent(in)  :: lats_nodes_a_fake
      integer,intent(out) :: gl_lats_index
      integer,intent(out) :: global_lats_a(latg)
      integer, intent(in) :: global_time_sort_index(latg)

      integer :: ijk
      integer :: jptlats
      integer :: lat
      integer :: node

      lat = 1

!.............................................
      do ijk=1,latg
         do node=1,nodes
            if (node.eq.me_fake+1) then
               gl_lats_index=gl_lats_index+1
               global_lats_a(gl_lats_index) = global_time_sort_index(lat)
            endif
            lat = lat + 1
            if (lat .gt. latg) go to 200
         enddo

         do node=nodes,1,-1
            if (node.eq.me_fake+1) then
               gl_lats_index=gl_lats_index+1
               global_lats_a(gl_lats_index) = global_time_sort_index(lat)
            endif
            lat = lat + 1
            if (lat .gt. latg) go to 200
         enddo
      enddo
  200 continue
      return
      end

      end module get_lats_node_a_stochy_mod
