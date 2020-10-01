!>@brief The module 'spectral_layout_mod' contains the gaussian grid domain decompostion
! and the subroutine to interpolate from the gaussian grid to cubed-sphere (or any lat-lon pair)
module spectral_layout_mod

      implicit none

!
! program log:
! 20161011   philip pegion :  make stochastic pattern generator standalone
!
! 20190503   dom heinzeller : add ompthreads and stochy_la2ga; todo: cleanup nodes, me, ... (defined multiple times in several files)
!
   integer :: nodes,             &
              me,                &
              master,            &
              lon_dim_a,         &
              ls_dim,            &
              ls_max_node,       &
              lats_dim_a,        &
              lats_node_a,       &
              lats_node_a_max,   &
              ipt_lats_node_a,   &
              len_trie_ls,       &
              len_trio_ls,       &
              len_trie_ls_max,   &
              len_trio_ls_max,   &
              me_l_0,            &
              lats_dim_ext,      &
              lats_node_ext,     &
              ipt_lats_node_ext

   integer :: ompthreads = 1
!
   integer, allocatable :: lat1s_a(:), lon_dims_a(:),lon_dims_ext(:)

contains

   logical function is_master()
      if (me == master) then
         is_master = .true.
      else
         is_master = .false.
      end if
   end function is_master

   !  interpolation from lat/lon or gaussian grid to other lat/lon grid
   !
!>@brief The subroutine 'stochy_la2ga' intepolates from the global gaussian grid
!! to the cubed sphere points
!>@details This code is taken from the legacy spectral GFS
   subroutine stochy_la2ga(regin,imxin,jmxin,rinlon,rinlat,rlon,rlat, &
                           gauout,len,rslmsk, outlat, outlon)
      use kinddef , only : kind_io8, kind_io4
      implicit none
      ! interface variables
      real (kind=kind_io8), intent(in)  :: regin(imxin,jmxin)
      integer,              intent(in)  :: imxin
      integer,              intent(in)  :: jmxin
      real (kind=kind_io8), intent(in)  :: rinlon(imxin)
      real (kind=kind_io8), intent(in)  :: rinlat(jmxin)
      real (kind=kind_io8), intent(in)  :: rlon
      real (kind=kind_io8), intent(in)  :: rlat
      real (kind=kind_io8), intent(out) :: gauout(len)
      integer,              intent(in)  :: len
      real (kind=kind_io8), intent(in)  :: rslmsk(imxin,jmxin)
      real (kind=kind_io8), intent(in)  :: outlat(len)
      real (kind=kind_io8), intent(in)  :: outlon(len)
      ! local variables
      real (kind=kind_io8) :: wei4,wei3,wei2,sum2,sum1,sum3,wei1,sum4
      real (kind=kind_io8) :: wsum,wsumiv,sums,sumn,wi2j2,x,y,wi1j1
      real (kind=kind_io8) :: wi1j2,wi2j1,aphi,rnume,alamd,denom
      integer              :: jy,ix,i,j,jq,jx
      integer              :: j1,j2,ii,i1,i2,kmami,it
      integer              :: nx,kxs,kxt
      integer              :: iindx1(len)
      integer              :: iindx2(len)
      integer              :: jindx1(len)
      integer              :: jindx2(len)
      real(kind=kind_io8)  :: ddx(len)
      real(kind=kind_io8)  :: ddy(len)
      real(kind=kind_io8)  :: wrk(len)
      integer              :: len_thread_m
      integer              :: len_thread
      integer              :: i1_t
      integer              :: i2_t
!
      len_thread_m  = (len+ompthreads-1) / ompthreads
!
      !$omp parallel do num_threads(ompthreads) default(none)               &
      !$omp private(i1_t,i2_t,len_thread,it,i,ii,i1,i2)                     &
      !$omp private(j,j1,j2,jq,ix,jy,nx,kxs,kxt,kmami)                      &
      !$omp private(alamd,denom,rnume,aphi,x,y,wsum,wsumiv,sum1,sum2)       &
      !$omp private(sum3,sum4,wi1j1,wi2j1,wi1j2,wi2j2,wei1,wei2,wei3,wei4)  &
      !$omp private(sumn,sums)                                              &
      !$omp shared(imxin,jmxin)                                             &
      !$omp shared(outlon,outlat,wrk,iindx1,rinlon,jindx1,rinlat,ddx,ddy)   &
      !$omp shared(rlon,rlat,regin,gauout)                                  &
      !$omp shared(ompthreads,len_thread_m,len,iindx2,jindx2,rslmsk)
      do it=1,ompthreads ! start of threaded loop
        i1_t       = (it-1)*len_thread_m+1
        i2_t       = min(i1_t+len_thread_m-1,len)
        len_thread = i2_t-i1_t+1
!
!       find i-index for interpolation
!
        do i=i1_t, i2_t
          alamd = outlon(i)
          if (alamd .lt. rlon)   alamd = alamd + 360.0
          if (alamd .gt. 360.0+rlon) alamd = alamd - 360.0
          wrk(i)    = alamd
          iindx1(i) = imxin
        enddo
        do i=i1_t,i2_t
          do ii=1,imxin
            if(wrk(i) .ge. rinlon(ii)) iindx1(i) = ii
          enddo
        enddo
        do i=i1_t,i2_t
          i1 = iindx1(i)
          if (i1 .lt. 1) i1 = imxin
          i2 = i1 + 1
          if (i2 .gt. imxin) i2 = 1
          iindx1(i) = i1
          iindx2(i) = i2
          denom     = rinlon(i2) - rinlon(i1)
          if(denom.lt.0.) denom = denom + 360.
          rnume = wrk(i) - rinlon(i1)
          if(rnume.lt.0.) rnume = rnume + 360.
          ddx(i) = rnume / denom
        enddo
!
!  find j-index for interplation
!
        if(rlat.gt.0.) then
          do j=i1_t,i2_t
            jindx1(j)=0
          enddo
          do jx=1,jmxin
            do j=i1_t,i2_t
              if(outlat(j).le.rinlat(jx)) jindx1(j) = jx
            enddo
          enddo
          do j=i1_t,i2_t
            jq = jindx1(j)
            aphi=outlat(j)
            if(jq.ge.1 .and. jq .lt. jmxin) then
              j2=jq+1
              j1=jq
             ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
            elseif (jq .eq. 0) then
              j2=1
              j1=1
              if(abs(90.-rinlat(j1)).gt.0.001) then
                ddy(j)=(aphi-rinlat(j1))/(90.-rinlat(j1))
              else
                ddy(j)=0.0
              endif
            else
              j2=jmxin
              j1=jmxin
              if(abs(-90.-rinlat(j1)).gt.0.001) then
                ddy(j)=(aphi-rinlat(j1))/(-90.-rinlat(j1))
              else
                ddy(j)=0.0
              endif
            endif
            jindx1(j)=j1
            jindx2(j)=j2
          enddo
        else
          do j=i1_t,i2_t
            jindx1(j) = jmxin+1
          enddo
          do jx=jmxin,1,-1
            do j=i1_t,i2_t
              if(outlat(j).le.rinlat(jx)) jindx1(j) = jx
            enddo
          enddo
          do j=i1_t,i2_t
            jq = jindx1(j)
            aphi=outlat(j)
            if(jq.gt.1 .and. jq .le. jmxin) then
              j2=jq
              j1=jq-1
              ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
            elseif (jq .eq. 1) then
              j2=1
              j1=1
              if(abs(-90.-rinlat(j1)).gt.0.001) then
                 ddy(j)=(aphi-rinlat(j1))/(-90.-rinlat(j1))
              else
                 ddy(j)=0.0
              endif
            else
              j2=jmxin
              j1=jmxin
              if(abs(90.-rinlat(j1)).gt.0.001) then
                 ddy(j)=(aphi-rinlat(j1))/(90.-rinlat(j1))
              else
                 ddy(j)=0.0
              endif
            endif
            jindx1(j)=j1
            jindx2(j)=j2
          enddo
        endif
!
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        do i=1,imxin
          sum1 = sum1 + regin(i,1)
          sum2 = sum2 + regin(i,jmxin)
        enddo
        sum1 = sum1 / imxin
        sum2 = sum2 / imxin
        sum3 = sum1
        sum4 = sum2
!
!  quasi-bilinear interpolation
!
        do i=i1_t,i2_t
          y  = ddy(i)
          j1 = jindx1(i)
          j2 = jindx2(i)
          x  = ddx(i)
          i1 = iindx1(i)
          i2 = iindx2(i)
!
          wi1j1 = (1.-x) * (1.-y)
          wi2j1 =     x  *( 1.-y)
          wi1j2 = (1.-x) *      y
          wi2j2 =     x  *      y
!
          wsum   = wi1j1 + wi2j1 + wi1j2 + wi2j2
          wrk(i) = wsum
          if(wsum.ne.0.) then
            wsumiv = 1./wsum
            if(j1.ne.j2) then
              gauout(i) = (wi1j1*regin(i1,j1) + wi2j1*regin(i2,j1) + &
                           wi1j2*regin(i1,j2) + wi2j2*regin(i2,j2))  &
                        *wsumiv
            else
              if (rlat .gt. 0.0) then
                sumn = sum3
                sums = sum4
                if( j1 .eq. 1) then
                  gauout(i) = (wi1j1*sumn        +wi2j1*sumn        + &
                               wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2)) &
                            * wsumiv
                elseif (j1 .eq. jmxin) then
                  gauout(i) = (wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+ &
                               wi1j2*sums        +wi2j2*sums        ) &
                            * wsumiv
                endif
              else
                sums = sum3
                sumn = sum4
                if( j1 .eq. 1) then
                  gauout(i) = (wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+ &
                               wi1j2*sums        +wi2j2*sums        ) &
                            * wsumiv
                elseif (j1 .eq. jmxin) then
                  gauout(i) = (wi1j1*sumn        +wi2j1*sumn        + &
                               wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2)) &
                            * wsumiv
                endif
              endif
            endif            ! if j1 .ne. j2
          endif
        enddo
        do i=i1_t,i2_t
          j1 = jindx1(i)
          j2 = jindx2(i)
          i1 = iindx1(i)
          i2 = iindx2(i)
          if(wrk(i) .eq. 0.0) then
            write(6,*) ' la2ga: bad rslmsk given'
            call sleep(2)
            stop
          endif
        enddo
      enddo ! end of threaded loop
!$omp end parallel do
!
      return
!
   end subroutine stochy_la2ga

end module spectral_layout_mod
