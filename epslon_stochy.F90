!>@brief The module 'epslon_stochy_mod' contains the subroute epslon_stochy
      module epslon_stochy_mod

      implicit none

      contains

!>@brief The subroutine 'epslon_stochy' calculate coeffients for use in spectral space
!>@details This code is taken from the legacy spectral GFS
      subroutine epslon_stochy(gis_stochy)

      use spectral_layout_mod
      use stochy_internal_state_mod, only : stochy_internal_state
      use kinddef
      implicit none

      type(stochy_internal_state), intent(inout) :: gis_stochy

      integer                  ls_node(ls_dim,3)

!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod

      integer                  l,locl,n

      integer                  indev
      integer                  indod

      real(kind_dbl_prec) f1,f2,rn,val

      real(kind_dbl_prec) cons0     !constant

      integer                  indlsev,jbasev
      integer                  indlsod,jbasod

      include 'function2'

      cons0=0.0d0     !constant
!c
!c......................................................................
!c
      do locl=1,ls_max_node
              l=gis_stochy%ls_node(locl,1)
         jbasev=gis_stochy%ls_node(locl,2)
         indev=indlsev(l,l)
         gis_stochy%epse  (indev)=cons0     !constant
         gis_stochy%epsedn(indev)=cons0     !constant
          indev=indev+1

         do n=l+2,jcap+1,2
            rn=n
            f1=n*n-l*l
            f2=4*n*n-1
            val=sqrt(f1/f2)
            gis_stochy%epse  (indev)=val
            gis_stochy%epsedn(indev)=val/rn
             indev=indev+1
         enddo
      enddo
      do locl=1,ls_max_node
              l=gis_stochy%ls_node(locl,1)
         jbasod=gis_stochy%ls_node(locl,3)
         indod=indlsod(l+1,l)

         do n=l+1,jcap+1,2
            rn=n
            f1=n*n-l*l
            f2=4*n*n-1
            val=sqrt(f1/f2)
            gis_stochy%epso  (indod)=val
            gis_stochy%epsodn(indod)=val/rn
             indod=indod+1
         enddo
      enddo

      return
      end

      end module epslon_stochy_mod
