!>@brief The module 'getcon_spectral_mod' contains the subroutine getcon_spectral
module getcon_spectral_mod

   implicit none

   contains

!>@brief The subroutine 'getcon_spectral' gets various constants for the spectral and related gaussian grid
!! and caluated the assoicate legendre polynomials
!>@details This code is taken from the legacy spectral GFS
      subroutine getcon_spectral ( ls_node,ls_nodes,max_ls_nodes,  &
                                  lats_nodes_a,global_lats_a,      &
                                  lonsperlat,latsmax,              &
                                  epse,epso,epsedn,epsodn,         &
                                  snnp1ev,snnp1od,                 &
                                  plnev_a,plnod_a,                 &
                                  plnew_a,plnow_a)

! program log:
! 20110220    henry juang update code to fit mass_dp and ndslfv
! 20201002    philip pegion  clean up of code
!
      use epslon_stochy_mod, only: epslon_stochy
      use get_lats_node_a_stochy_mod, only: get_lats_node_a_stochy
      use get_ls_node_stochy_mod, only: get_ls_node_stochy
      use glats_stochy_mod, only: glats_stochy
      use gozrineo_a_stochy_mod, only: gozrineo_a_stochy
      use pln2eo_a_stochy_mod, only: pln2eo_a_stochy
      use setlats_a_stochy_mod, only: setlats_a_stochy
      use spectral_layout_mod
      use stochy_internal_state_mod
      use kinddef 

      implicit none
!
      integer              i,j,l,lat,n
      integer              ls_node(ls_dim,3)
!
!     ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!     ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!     ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!
      integer                      ls_nodes(ls_dim,nodes)
      integer, dimension(nodes) :: max_ls_nodes,  lats_nodes_a
      integer, dimension(latg)  :: global_lats_a, lonsperlat
!
!
      real(kind=kind_dbl_prec), dimension(len_trie_ls) :: epse, epsedn, snnp1ev
      real(kind=kind_dbl_prec), dimension(len_trio_ls) :: epso, epsodn, snnp1od
!
      real(kind=kind_dbl_prec), dimension(len_trie_ls,latg2) :: plnev_a, plnew_a
      real(kind=kind_dbl_prec), dimension(len_trio_ls,latg2) :: plnod_a, plnow_a
!
      integer       iprint,locl,node,&
                    len_trie_ls_nod, len_trio_ls_nod,&
                    indev, indod, indlsev,jbasev,indlsod,jbasod
!
      integer gl_lats_index, latsmax
      integer global_time_sort_index_a(latg)
!
      include 'function2'
!
      real(kind=kind_dbl_prec) global_time_a(latg)
!
      real(kind=kind_dbl_prec), parameter :: cons0 = 0.d0, cons0p5  = 0.5d0,&
                                         cons1 = 1.d0, cons0p92 = 0.92d0
!
      gl_lats_index = 0
      global_lats_a = -1
      do lat = 1,latg                  !my intialize global_time_a to lonsperlat
          global_time_a(lat) = lonsperlat(lat)
      enddo

      do lat = 1, latg2
         lonsperlat(latg+1-lat) = lonsperlat(lat)
      end do
      do node=1,nodes
          call get_lats_node_a_stochy( node-1, global_lats_a,lats_nodes_a(node),&
                               gl_lats_index,global_time_sort_index_a, iprint)
      enddo
      call setlats_a_stochy(lats_nodes_a,global_lats_a,iprint, lonsperlat)

      iprint = 0
      do node=1,nodes
         call get_ls_node_stochy( node-1, ls_nodes(1,node),max_ls_nodes(node), iprint )
      enddo
!
      len_trie_ls_max = 0
      len_trio_ls_max = 0
      do node=1,nodes
!
         len_trie_ls_nod = 0
         len_trio_ls_nod = 0
         do locl=1,max_ls_nodes(node)
            l=ls_nodes(locl,node)
            len_trie_ls_nod = len_trie_ls_nod+(jcap+3-l)/2
            len_trio_ls_nod = len_trio_ls_nod+(jcap+2-l)/2
         enddo
         len_trie_ls_max = max(len_trie_ls_max,len_trie_ls_nod)
         len_trio_ls_max = max(len_trio_ls_max,len_trio_ls_nod)
!
      enddo
!
      iprint = 0
!
      lats_dim_a = 0
      do node=1,nodes
         lats_dim_a = max(lats_dim_a,lats_nodes_a(node))
      enddo
      lats_node_a = lats_nodes_a(me+1)

      lats_node_a_max = 0
      do i=1,nodes
        lats_node_a_max = max(lats_node_a_max, lats_nodes_a(i))
      enddo
      latsmax = lats_node_a_max

!
      ipt_lats_node_a   = 1
      if ( me > 0 ) then
        do node=1,me
          ipt_lats_node_a = ipt_lats_node_a + lats_nodes_a(node)
        enddo
      endif

!
      iprint = 0
!
           call glats_stochy(latg2,colrad_a,wgt_a,rcs2_a)
           call epslon_stochy(epse,epso,epsedn,epsodn,ls_node)
           call pln2eo_a_stochy(plnev_a,plnod_a,epse,epso,ls_node,latg2)
           call gozrineo_a_stochy(plnev_a,plnod_a, &
                plnew_a,plnow_a,epse,epso,ls_node,latg2)
!
!
      do locl=1,ls_max_node
              l = ls_node(locl,1)
         jbasev = ls_node(locl,2)
         indev  = indlsev(l,l)
         do n = l, jcap, 2
            snnp1ev(indev) = n*(n+1)
              indev        = indev+1
         end do
      end do
!
!
      do locl=1,ls_max_node
              l = ls_node(locl,1)
         jbasod = ls_node(locl,3)
         if ( l <= jcap-1 ) then
            indod = indlsod(l+1,l)
            do n = l+1, jcap, 2
               snnp1od(indod) = n*(n+1)
                 indod        = indod+1
            end do
         end if
      end do
!
!
      do locl=1,ls_max_node
              l = ls_node(locl,1)
         jbasev = ls_node(locl,2)
         jbasod = ls_node(locl,3)
         if (mod(L,2) == mod(jcap+1,2)) then ! set even (n-l) terms of top row to zero
            snnp1ev(indlsev(jcap+1,l)) = cons0
         else                                ! set odd (n-l) terms of top row to zero
            snnp1od(indlsod(jcap+1,l)) = cons0
         endif
      enddo
!
      do j=1,latg
        if( j <= latg2 ) then
          sinlat_a(j) =  cos(colrad_a(j))
        else
          sinlat_a(j) = -cos(colrad_a(latg+1-j))
        endif
        coslat_a(j) = sqrt(1.-sinlat_a(j)*sinlat_a(j))
      enddo
!
      do L=0,jcap
         do lat = 1, latg2
            if ( L <= min(jcap,lonsperlat(lat)/2) ) then
               lat1s_a(L) = lat
               go to 200
            endif
         end do
  200    continue
      end do
!

      do j=1,lats_node_a
         lat = global_lats_a(ipt_lats_node_a-1+j)
         if ( lonsperlat(lat) == lonf ) then
            lon_dims_a(j) = lonfx
         else
            lon_dims_a(j) = lonsperlat(lat) + 2
         endif
      enddo
!
      return
      end

end module getcon_spectral_mod
