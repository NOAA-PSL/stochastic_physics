module update_ca

use halo_exchange,    only: atmosphere_scalar_field_halo
use mersenne_twister, only: random_gauss,random_stat,random_number
use mpi_wrapper,      only: mp_reduce_sum,mp_bcst,mp_reduce_min,mp_reduce_max
use mpp_domains_mod,  only: domain2D

implicit none

!L. Bengtsson 2017-06
!Evolve the cellular automata in time


contains

subroutine update_cells_sgs(kstep,nca,nxc,nyc,nxch,nych,nlon,nlat,isc,iec,jsc,jec, &
                        npx,npy,domain_for_coupler,CA,ca_plumes,iini,ilives,       &
                        nlives,ncells,nfracseed,nseed,nthresh,nspinup,nf,nca_plumes)

implicit none

integer, intent(in) :: kstep,nxc,nyc,nlon,nlat,nxch,nych,nca,isc,iec,jsc,jec,npx,npy
integer, intent(in) :: iini(nxc,nyc,nca)
integer, intent(inout) :: ilives(nxc,nyc,nca)
real, intent(out) :: CA(nlon,nlat)
integer, intent(out) :: ca_plumes(nlon,nlat)
integer, intent(in) :: nlives, ncells, nseed, nspinup, nf
real, intent(in) :: nfracseed, nthresh
logical,intent(in) :: nca_plumes
type(domain2D), intent(inout) :: domain_for_coupler
real, dimension(nlon,nlat) :: frac
integer, dimension(nlon,nlat) :: maxlives
integer,allocatable,save :: board(:,:,:), lives(:,:,:)
integer,allocatable :: V(:),L(:),B(:)
integer,allocatable :: AG(:,:)
integer :: inci, incj, i, j, k,sub,spinup,it,halo,k_in,isize,jsize
integer :: ih, jh,kend
real, allocatable :: field_in(:,:),board_halo(:,:,:)
integer, dimension(nxc,nyc) :: neighbours, birth, newlives,thresh,maxliveshigh
integer, dimension(nxc,nyc) :: neg, newcell, oldlives, newval,temp,newseed
integer, dimension(ncells,ncells) :: onegrid

real, dimension(nxc,nyc) :: NOISE_B
real, dimension(nxc*nyc) :: noise1D2


!-------------------------------------------------------------------------------------------------
halo=1
isize=nlon+2*halo
jsize=nlat+2*halo
k_in=1

 if (.not. allocated(board))then
 allocate(board(nxc,nyc,nca))
 endif
 if (.not. allocated(lives))then
 allocate(lives(nxc,nyc,nca))
 endif
 if (.not. allocated(field_in))then
 allocate(field_in(nxc*nyc,1))
 endif
 if(.not. allocated(board_halo))then
 allocate(board_halo(nxch,nych,1))   
 endif
  
 noise1D2 = 0.0
!
 call random_number(noise1D2)

 !Put on 2D:
 do j=1,nyc
  do i=1,nxc
  NOISE_B(i,j)=noise1D2(i+(j-1)*nxc)    
  enddo
 enddo

  if(kstep <= 1)then
   do j=1,nyc
    do i=1,nxc
     board(i,j,nf) = 0
     lives(i,j,nf) = 0
    enddo
   enddo
  endif

  if(kstep == 2)then !Initiate CA at kstep 2 as physics field is empty at 0 and 1.
   do j=1,nyc
    do i=1,nxc
    board(i,j,nf) = iini(i,j,nf)
    lives(i,j,nf) = ilives(i,j,nf)*iini(i,j,nf)
   enddo
   enddo
  endif

 !Seed with new CA cells at each nseed step

  newseed = 0
  if(mod(kstep,nseed) == 0 .and. kstep >= 2)then
   do j=1,nyc
    do i=1,nxc
     if(board(i,j,nf) == 0 .and. NOISE_B(i,j)>0.95 )then
       newseed(i,j) = 1
     endif
     board(i,j,nf) = board(i,j,nf) + newseed(i,j)
    enddo
   enddo

  endif

  if(kstep == 2)then
  spinup=nspinup
  else
  spinup = 1
  endif


 do it = 1,spinup
 
!Step 1 - Initialize variables to 0 and extract the halo

 CA=0
 neighbours=0
 birth=0
 newlives=0
 neg=0
 newcell=0
 oldlives=0
 newval=0
 frac=0
 board_halo=0
 field_in=0
 maxlives = 0
 maxliveshigh =0
 

!Step 4 - Compute the neighbourhood
 
  do j=1,nyc
   do i=1,nxc
   field_in(i+(j-1)*nxc,1)=board(i,j,nf)
   enddo
  enddo

  call atmosphere_scalar_field_halo(board_halo,halo,nxch,nych,k_in,field_in, &
                                    isc,iec,jsc,jec,npx,npy,domain_for_coupler)


  do j=1,nyc
     do i=1,nxc
        ih=i+halo
        jh=j+halo
        neighbours(i,j)=board_halo(ih-1,jh-1,1)+board_halo(ih-1,jh,1)+ &
                        board_halo(ih-1,jh+1,1)+board_halo(ih,jh+1,1)+board_halo(ih+1,jh+1,1)+&
                        board_halo(ih+1,jh,1)+board_halo(ih+1,jh-1,1)+board_halo(ih,jh-1,1)
     enddo
  enddo

! Step 5 - Check rules; 

 !birth
 
  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j) == 3 .or. neighbours(i,j) ==2)then 
     birth(i,j)=1
     endif
   enddo
  enddo
 
 !death                                                                                                                        
  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j) < 2 .or. neighbours(i,j) > 3)then  
     lives(i,j,nf)=lives(i,j,nf) - 1
     endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
   if(lives(i,j,nf) < 0)then
     lives(i,j,nf)=0
   endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
    if(birth(i,j)==1 .and. lives(i,j,nf)==0)then
    newcell(i,j)=1
    endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
    lives(i,j,nf)=lives(i,j,nf)+newcell(i,j)*ilives(i,j,nf)
   enddo
  enddo

   do j=1,nyc
   do i=1,nxc
    if(neighbours(i,j)==3 .or. (board(i,j,nf)==1 .and. neighbours(i,j)==2))then
    board(i,j,nf)=1
    else
    board(i,j,nf)=0
    endif
   enddo
  enddo


enddo !spinup

!COARSE-GRAIN BACK TO NWP GRID
 
  inci=ncells
  incj=ncells
  sub=ncells-1
  DO j=1,nlat
     DO i=1,nlon
        CA(i,j)=(SUM(lives(inci-sub:inci,incj-sub:incj,nf)))/real(ncells*ncells)
        inci=inci+ncells
     ENDDO
     inci=ncells
     incj=incj+ncells
  ENDDO

if(nca_plumes) then
!COMPUTE NUMBER OF CLUSTERS (CONVECTIVE PLUMES) IN EACH CA-CELL
!Note, at the moment we only use the count of the plumes found in a grid-cell
!In the future the routine "plumes" can also be used to give the size of 
!each individual plume for better coupling to the convection scheme.

  temp=0
  do j=1,nyc
   do i=1,nxc
     if(lives(i,j,1) > 0)then
      temp(i,j)=1  
     endif
   enddo
  enddo

  kend=ceiling((ncells*ncells)/2.)
  if (.not. allocated(V))then
  allocate(V(kend))
  endif
  if (.not. allocated(L))then
  allocate(L(kend))
  endif
  if (.not. allocated(B))then
  allocate(B(kend))
  endif
  if (.not. allocated(AG))then
  allocate(AG(ncells,ncells))
  endif
  
  ca_plumes(:,:)=0
  inci=ncells
  incj=ncells
  sub=ncells-1
  DO j=1,nlat
     DO i=1,nlon
        B(:)=0
        L(:)=0
        V(:)=0
        onegrid(1:ncells,1:ncells)=temp(inci-sub:inci,incj-sub:incj)
        call plumes(V,L,AG,onegrid,ncells,ncells,kend)
        do k=1,kend
           if(V(k)==1)then
              B(k)=L(k) !to avoid considering clusters of 0
           endif
        enddo
        ca_plumes(i,j)=MAXVAL(B(1:kend))
        inci=inci+ncells
     ENDDO
     inci=ncells
     incj=incj+ncells
  ENDDO

else

ca_plumes(:,:)=0.

endif ! nca_plumes

end subroutine update_cells_sgs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_cells_global(kstep,nca,nxc,nyc,nxch,nych,nlon,nlat,isc,iec,jsc,jec, &
                        npx,npy,domain_for_coupler,CA,iini_g,ilives_g,                &
                        nlives,ncells,nfracseed,nseed,nthresh,nspinup,nf)

implicit none

integer, intent(in) :: kstep,nxc,nyc,nlon,nlat,nxch,nych,nca,isc,iec,jsc,jec,npx,npy
integer, intent(in) :: iini_g(nxc,nyc,nca), ilives_g(nxc,nyc)
real, intent(out) :: CA(nlon,nlat)
integer, intent(in) :: nlives, ncells, nseed, nspinup, nf
real, intent(in) :: nfracseed, nthresh
type(domain2D), intent(inout) :: domain_for_coupler
real, dimension(nlon,nlat) :: frac
integer,allocatable,save :: board_g(:,:,:), lives_g(:,:,:)
integer,allocatable :: V(:),L(:)
integer :: inci, incj, i, j, k ,sub,spinup,it,halo,k_in,isize,jsize
integer :: ih, jh
real, allocatable :: field_in(:,:),board_halo(:,:,:)
integer, dimension(nxc,nyc) :: neighbours, birth, newlives, thresh
integer, dimension(nxc,nyc) :: neg, newcell, oldlives, newval,temp,newseed
real, dimension(nxc,nyc) :: NOISE_B
real, dimension(nxc*nyc) :: noise1D2


!-------------------------------------------------------------------------------------------------
halo=1
isize=nlon+2*halo
jsize=nlat+2*halo
k_in=1

 if (.not. allocated(board_g))then
 allocate(board_g(nxc,nyc,nca))
 endif
 if (.not. allocated(lives_g))then
 allocate(lives_g(nxc,nyc,nca))
 endif
 if (.not. allocated(field_in))then
 allocate(field_in(nxc*nyc,1))
 endif
 if(.not. allocated(board_halo))then                                                                      
 allocate(board_halo(nxch,nych,1))   
 endif

 
 !random numbers:
 noise1D2 = 0.0

 call random_number(noise1D2)

 !Put on 2D:                                                                                                                                                                          
 do j=1,nyc
  do i=1,nxc
  NOISE_B(i,j)=noise1D2(i+(j-1)*nxc)
  enddo
 enddo                        


  if(kstep == 0)then

   do j=1,nyc
    do i=1,nxc
     board_g(i,j,nf) = iini_g(i,j,nf)
     lives_g(i,j,nf) = ilives_g(i,j)*iini_g(i,j,nf)
    enddo
   enddo

  endif

 !Seed with new CA cells at each nseed step
  newseed=0

  if(mod(kstep,nseed) == 0 .and. kstep > 0.)then
   do j=1,nyc
    do i=1,nxc
     if(board_g(i,j,nf) == 0 .and. NOISE_B(i,j)>0.95 )then
        newseed(i,j)=1
     endif
     board_g(i,j,nf) = board_g(i,j,nf) + newseed(i,j)
    enddo
   enddo
  endif

  if(kstep == 0)then
  spinup=nspinup
  else
  spinup = 1
  endif
 

do it=1,spinup

 
!Step 2 - Initialize variables to 0 and extract the halo
 
 neighbours=0
 birth=0
 newlives=0
 neg=0
 newcell=0
 oldlives=0
 newval=0
 frac=0
 board_halo=0
 field_in=0

!The input to scalar_field_halo needs to be 1D.
!take the updated board_g fields and extract the halo
! in order to have updated values in the halo region. 


  do j=1,nyc
   do i=1,nxc
   field_in(i+(j-1)*nxc,1)=board_g(i,j,nf)
   enddo
  enddo

  call atmosphere_scalar_field_halo(board_halo,halo,nxch,nych,k_in,field_in, &
                                    isc,iec,jsc,jec,npx,npy,domain_for_coupler)


  do j=1,nyc
     do i=1,nxc
        ih=i+halo
        jh=j+halo
        neighbours(i,j)=board_halo(ih-1,jh-1,1)+board_halo(ih-1,jh,1)+ &
                        board_halo(ih-1,jh+1,1)+board_halo(ih,jh+1,1)+board_halo(ih+1,jh+1,1)+&
                        board_halo(ih+1,jh,1)+board_halo(ih+1,jh-1,1)+board_halo(ih,jh-1,1)
     enddo
  enddo



  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j)==2 .or. neighbours(i,j)==3)then
     birth(i,j)=1
     endif
   enddo
  enddo


  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j)<2 .or. neighbours(i,j)>3)then
     lives_g(i,j,nf)=lives_g(i,j,nf) - 1
     endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
   if(lives_g(i,j,nf)<0)then
     lives_g(i,j,nf)=0
   endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
    if(birth(i,j)==1 .and. lives_g(i,j,nf)==0)then
    newcell(i,j)=1
    endif
   enddo
  enddo


  do j=1,nyc
   do i=1,nxc
    lives_g(i,j,nf)=lives_g(i,j,nf)+newcell(i,j)*ilives_g(i,j)
   enddo
  enddo


   do j=1,nyc
   do i=1,nxc
    if( (board_g(i,j,nf) ==1 .and. (neighbours(i,j)==3 .or. neighbours(i,j)==2) ).or. (board_g(i,j,nf)==0 .and. neighbours(i,j)==3) )then
    board_g(i,j,nf)=1
    else
    board_g(i,j,nf)=0
    endif
   enddo
  enddo

enddo !spinup

!COARSE-GRAIN BACK TO NWP GRID
 
  inci=ncells
  incj=ncells
  sub=ncells-1
  DO j=1,nlat
     DO i=1,nlon
        CA(i,j)=(SUM(lives_g(inci-sub:inci,incj-sub:incj,nf)))/(ncells*ncells)
        inci=inci+ncells
     ENDDO
     inci=ncells
     incj=incj+ncells
  ENDDO

end subroutine update_cells_global

end module update_ca
