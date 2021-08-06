!>@brief The module 'gozrineo_a_stochy_mod' contains the subroutine 'gozrineo_a_stochy'
      module gozrineo_a_stochy_mod

      implicit none

      contains
!>@brief The subroutine 'gozrineo_a_stochy' calculates the deriviates of assoicate legendre polynomials
!>@details This code is taken from the legacy spectral GFS
      subroutine gozrineo_a_stochy(gis_stochy, num_lat)

      use spectral_layout_mod
      use kinddef
      use stochy_internal_state_mod, only : stochy_internal_state
      implicit none

      type(stochy_internal_state), intent(inout) :: gis_stochy
      integer,  intent(in)                       :: num_lat

      integer                  l,lat,locl,n
      integer                  indev,indev1,indev2
      integer                  indod,indod1,indod2
      integer                  inddif

      real(kind=kind_dbl_prec) rn,rnp1,wcsa

      real(kind=kind_dbl_prec) cons0     !constant
      real(kind=kind_dbl_prec) cons2     !constant
      real  rerth

      integer                  indlsev,jbasev
      integer                  indlsod,jbasod

      include 'function2'


      cons0 = 0.d0     !constant
      cons2 = 2.d0     !constant
      rerth  =6.3712e+6      ! radius of earth (m)


      do lat=1,num_lat

         wcsa=rcs2_a(lat)/rerth

         do locl=1,ls_max_node
                 l=gis_stochy%ls_node(locl,1)
            jbasev=gis_stochy%ls_node(locl,2)
            jbasod=gis_stochy%ls_node(locl,3)
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
               gis_stochy%plnew_a(indev,lat) = gis_stochy%plnev_a(indev,lat) * wgt_a(lat)
            enddo

            do indod = indod1 , indod2
               gis_stochy%plnow_a(indod,lat) = gis_stochy%plnod_a(indod,lat) * wgt_a(lat)
            enddo
         enddo
      enddo
      return
      end

      end module gozrineo_a_stochy_mod
