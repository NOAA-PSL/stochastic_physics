!>@brief The module 'gozrineo_a_stochy_mod' contains the subroutine 'gozrineo_a_stochy'
      module gozrineo_a_stochy_mod

      implicit none

      contains
!>@brief The subroutine 'gozrineo_a_stochy' calculates the deriviates of assoicate legendre polynomials
!>@details This code is taken from the legacy spectral GFS
      subroutine gozrineo_a_stochy(plnev_a,plnod_a,
     &                      plnew_a,plnow_a,
     &                      epse,epso,ls_node,num_lat)
cc
      use spectral_layout_mod
      use kinddef
      implicit none
cc
      real(kind=kind_dbl_prec) plnev_a(len_trie_ls,latg2)
      real(kind=kind_dbl_prec) plnod_a(len_trio_ls,latg2)
      real(kind=kind_dbl_prec) plnew_a(len_trie_ls,latg2)
      real(kind=kind_dbl_prec) plnow_a(len_trio_ls,latg2)
cc
      real(kind=kind_dbl_prec)    epse(len_trie_ls)
      real(kind=kind_dbl_prec)    epso(len_trio_ls)
cc
cc
      integer                  ls_node(ls_dim,3)
cc
      integer                  num_lat
cc
      integer                  l,lat,locl,n
      integer                  indev,indev1,indev2
      integer                  indod,indod1,indod2
      integer                  inddif
cc
      real(kind=kind_dbl_prec) rn,rnp1,wcsa
cc
      real(kind=kind_dbl_prec) cons0     !constant
      real(kind=kind_dbl_prec) cons2     !constant
      real  rerth
cc
      integer                  indlsev,jbasev
      integer                  indlsod,jbasod
cc
      include 'function2'
cc
cc
      cons0 = 0.d0     !constant
      cons2 = 2.d0     !constant
      rerth  =6.3712e+6      ! radius of earth (m)
cc
cc
      do lat=1,num_lat
cc
         wcsa=rcs2_a(lat)/rerth
cc
cc
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
            jbasev=ls_node(locl,2)
            jbasod=ls_node(locl,3)
            indev1 = indlsev(L,L)
            indod1 = indlsod(L+1,L)
            if (mod(L,2).eq.mod(jcap+1,2)) then
               indev2 = indlsev(jcap+1,L)
               indod2 = indlsod(jcap  ,L)
            else
               indev2 = indlsev(jcap  ,L)
               indod2 = indlsod(jcap+1,L)
            endif
            do indev = indev1 , indev2
               plnew_a(indev,lat) = plnev_a(indev,lat) * wgt_a(lat)
            enddo
cc
            do indod = indod1 , indod2
               plnow_a(indod,lat) = plnod_a(indod,lat) * wgt_a(lat)
            enddo
cc
         enddo
cc
      enddo
cc
      return
      end

      end module gozrineo_a_stochy_mod
