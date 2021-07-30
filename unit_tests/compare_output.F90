program compare_output


use netcdf

implicit none

integer, parameter :: CRES=96
integer, parameter :: ntiles=6
integer, parameter :: nexpt=3
integer, allocatable :: layout_x(:),layout_y(:)

real, allocatable :: truth(:,:,:),expt(:,:,:),diff(:,:,:),accum_error(:)

integer :: i1,i2,j1,j2,nx,ny,ierr1,ierr2,ierr3,i,j,k,ncid,varid,ncid2,t,t2
character*2,tile1,tile2
character*1,lx,ly

allocate(layout_x(nexpt))
allocate(layout_y(nexpt))
allocate(accum_error(nexpt))
allocate(truth(CRES,CRES,2))
layout_x(1)=1
layout_x(2)=2
layout_x(3)=4
layout_y(1)=4
layout_y(2)=2
layout_y(3)=1

accum_error(:)=0.0

do k=1,nexpt
   nx=CRES/layout_x(k)
   ny=CRES/layout_y(k)
   allocate(expt(nx,ny,2))
   allocate(diff(nx,ny,2))
   t2=1
   do t=1,ntiles
      write(tile1,fmt='(I2.2)') t
      ierr1=nf90_open('layout_1x1/workg_T162_984x488.tile'//tile1//'.nc',mode=nf90_nowrite,ncid=ncid)
      ierr2=nf90_inq_varid(ncid,'sppt_wts',varid)
      ierr3=nf90_get_var(ncid,varid,truth,count=(/CRES,CRES,2/))
      if (ierr1+ierr2+ierr3.NE.0) then
         print*,'error reading in truth files'
         call exit(1)
      endif
      i1=1
      i2=i1+nx-1
      j1=1
      j2=j1+ny-1
      do j=1,layout_y(k)
         do i=1,layout_x(k)
            write(tile2,fmt='(I2.2)') t2
            write(lx,fmt='(I1)') layout_x(k)
            write(ly,fmt='(I1)') layout_y(k)
            ierr1=nf90_open('layout_'//lx//'x'//ly//'/workg_T162_984x488.tile'//tile2//'.nc',mode=nf90_nowrite,ncid=ncid2)
            ierr2=nf90_inq_varid(ncid2,'sppt_wts',varid)
            ierr3=nf90_get_var(ncid2,varid,expt,count=(/nx,ny,2/))
            if (ierr1+ierr2+ierr3.NE.0) then
               print*,'error reading in expt files'
               call exit(2)
            endif
!            print*,'i1,i2,j1,j2,tile=',i1,i2,j1,j2,t,t2
!            print*,lbound(expt,1),ubound(expt,1),lbound(expt,2),ubound(expt,2),lbound(truth,1),ubound(truth,1),lbound(truth,2),ubound(truth,2)
            diff(:,:,:)=expt(:,:,:)-truth(i1:i2,j1:j2,:)
!            print*,sum(abs(expt)),sum(abs(truth(i1:i2,j1:j2)))
            accum_error(k)=accum_error(k)+sum(abs(diff))
            i1=i1+nx
            i2=i2+nx
            if (i2.GT.CRES) then
               i1=1
               i2=i1+nx-1
               j1=j1+ny
               j2=j2+ny
            endif
            ierr1=nf90_close(ncid2)
            t2=t2+1
         enddo
      enddo
      ierr1=nf90_close(ncid)
   enddo
   deallocate(expt)
   deallocate(diff)
enddo

if (sum(accum_error(:)) .EQ. 0.0 )then
   print*,'all tests pass'
else
   print*,'decomposition test fail',accum_error
   call exit(3)
endif

end
