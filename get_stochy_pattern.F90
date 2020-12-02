!>@brief The module 'get_stochy_pattern_mod' contains the subroutines to retrieve the random pattern in the cubed-sphere grid
module get_stochy_pattern_mod
 use kinddef, only : kind_dbl_prec, kind_evod
 use spectral_layout_mod, only : ipt_lats_node_a, lat1s_a, lats_dim_a,      &
                                 lats_node_a, lon_dim_a, len_trie_ls,       &
                                 len_trio_ls, ls_dim, nodes, stochy_la2ga,  &
                                 coslat_a, latg, latg2, levs, lonf, skeblevs
 use stochy_namelist_def, only : n_var_lndp, ntrunc, stochini
 use stochy_data_mod, only : gg_lats, gg_lons, inttyp, nskeb, nshum, nsppt, &
                             rnlat, rpattern_sfc, rpattern_skeb,   &
                             rpattern_shum, rpattern_sppt, skebu_save,      &
                             skebv_save, skeb_vwts, skeb_vpts, wlon
 use stochy_patterngenerator_mod, only: random_pattern, ndimspec,           &
                                        patterngenerator_advance
 use stochy_internal_state_mod, only: stochy_internal_state
 use mpp_mod,only: mpp_npes,mpp_pe,mpp_root_pe,mpp_sum
#ifdef STOCHY_UNIT_TEST
 use standalone_stochy_module,   only: GFS_control_type, GFS_grid_type,random_seed
#else
 use GFS_typedefs,       only: GFS_control_type, GFS_grid_type
 use mersenne_twister, only: random_seed
#endif
 use dezouv_stochy_mod, only: dezouv_stochy
 use dozeuv_stochy_mod, only: dozeuv_stochy
 use four_to_grid_mod, only: four_to_grid
 use sumfln_stochy_mod, only: sumfln_stochy
 implicit none
 private

 public  get_random_pattern_vector   
 public  get_random_pattern_sfc,get_random_pattern_scalar
 public  dump_patterns
 logical :: first_call=.true.
 contains

!>@brief The subroutine 'get_random_pattern_sfc' converts spherical harmonics to the gaussian grid then interpolates to the target grid
!>@details This subroutine is for a 2-D (lat-lon) scalar field
subroutine get_random_pattern_sfc(rpattern,npatterns,&
           gis_stochy,pattern_3d)
!\callgraph

! generate a random pattern for stochastic physics
 implicit none
 type(random_pattern), intent(inout) :: rpattern(npatterns)
 type(stochy_internal_state), intent(in) :: gis_stochy
 integer,intent(in)::   npatterns

 integer i,j,lat,n,k
 real(kind=kind_dbl_prec), dimension(lonf,gis_stochy%lats_node_a,1):: wrk2d

! logical lprint

 real(kind=kind_dbl_prec), allocatable, dimension(:,:) :: workg
 real (kind=kind_dbl_prec)   glolal(lonf,gis_stochy%lats_node_a)
 integer kmsk0(lonf,gis_stochy%lats_node_a)
 real(kind=kind_dbl_prec),intent(out) :: pattern_3d(gis_stochy%nx,gis_stochy%ny,n_var_lndp)
 real(kind=kind_dbl_prec) :: pattern_1d(gis_stochy%nx)

 do k=1,n_var_lndp
   kmsk0 = 0
   glolal = 0.
   do n=1,npatterns
     if (mpp_pe()==mpp_root_pe()) print *, 'Random pattern for SFC-PERTS in get_random_pattern_sfc: k, min, max ',k,minval(rpattern_sfc(n)%spec_o(:,:,k)), maxval(rpattern_sfc(n)%spec_o(:,:,k))
     call scalarspect_to_gaugrid(                       &
         rpattern(n)%spec_e(:,:,k),rpattern(n)%spec_o(:,:,k),wrk2d,&
         gis_stochy%ls_node,gis_stochy%ls_nodes,gis_stochy%max_ls_nodes,&
         gis_stochy%lats_nodes_a,gis_stochy%global_lats_a,gis_stochy%lonsperlat,&
         gis_stochy%plnev_a,gis_stochy%plnod_a,1)
     glolal = glolal + wrk2d(:,:,1)
   enddo

   allocate(workg(lonf,latg))
   workg = 0.
   do j=1,gis_stochy%lats_node_a
     lat=gis_stochy%global_lats_a(ipt_lats_node_a-1+j)
     do i=1,lonf
        workg(i,lat) = glolal(i,j)
     enddo
   enddo

   call mpp_sum(workg,lonf*latg)
   if (mpp_pe()==mpp_root_pe()) print *, 'workg after mpp_sum for SFC-PERTS in get_random_pattern_sfc: k, min, max ',k,minval(workg), maxval(workg)

! interpolate to cube grid

   do j=1,gis_stochy%ny
      pattern_1d = 0
      associate( tlats=>gis_stochy%parent_lats(:,j),&
                 tlons=>gis_stochy%parent_lons(:,j))
      call stochy_la2ga(workg,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                        pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
      pattern_3d(:,j,k)=pattern_1d(:)
      end associate
   enddo
   if (mpp_pe()==mpp_root_pe()) print *, '3D pattern for SFC-PERTS in get_random_pattern_sfc: k, min, max ',k,minval(pattern_3d(:,:,k)), maxval(pattern_3d(:,:,k))
   deallocate(workg)

 enddo  ! loop over k, n_var_lndp

end subroutine get_random_pattern_sfc


!>@brief The subroutine 'get_random_pattern_fv3_vect' converts spherical harmonics to a vector on gaussian grid then interpolates to the target grid 
!>@details This subroutine is for a 2-D (lat-lon) vector field
subroutine get_random_pattern_vector(rpattern,npatterns,&
           gis_stochy,upattern_3d,vpattern_3d)
!\callgraph

! generate a random pattern for stochastic physics
 implicit none
 type(stochy_internal_state), intent(in) :: gis_stochy
 type(random_pattern), intent(inout) :: rpattern(npatterns)

 real(kind=kind_evod), dimension(len_trie_ls,2,1) ::  vrtspec_e,divspec_e
 real(kind=kind_evod), dimension(len_trio_ls,2,1) ::  vrtspec_o,divspec_o
 integer::   npatterns

 real(kind=kind_dbl_prec) :: upattern_3d(gis_stochy%nx,gis_stochy%ny,levs)
 real(kind=kind_dbl_prec) :: vpattern_3d(gis_stochy%nx,gis_stochy%ny,levs)
 real(kind=kind_dbl_prec) :: pattern_1d(gis_stochy%nx)
 integer i,j,lat,n,nn,k
 real(kind_dbl_prec), dimension(lonf,gis_stochy%lats_node_a,1):: wrk2du,wrk2dv

! logical lprint

 real, allocatable, dimension(:,:) :: workgu,workgv
 integer kmsk0(lonf,gis_stochy%lats_node_a)
 kmsk0 = 0
 allocate(workgu(lonf,latg))
 allocate(workgv(lonf,latg))
 print*,'before first_call',skeblevs
 if (first_call) then
    allocate(skebu_save(gis_stochy%nx,gis_stochy%ny,skeblevs))
    allocate(skebv_save(gis_stochy%nx,gis_stochy%ny,skeblevs))
    do k=2,skeblevs
       workgu = 0.
       workgv = 0.
       do n=1,npatterns
          if (.not. stochini) call patterngenerator_advance(rpattern(n),k,first_call)
      !   ke norm (convert streamfunction forcing to vorticity forcing)
          divspec_e = 0; divspec_o = 0.
          do nn=1,2
             vrtspec_e(:,nn,1) = gis_stochy%kenorm_e*rpattern(n)%spec_e(:,nn,k)
             vrtspec_o(:,nn,1) = gis_stochy%kenorm_o*rpattern(n)%spec_o(:,nn,k)
          enddo
        ! convert to winds
          call vrtdivspect_to_uvgrid(&
                 divspec_e,divspec_o,vrtspec_e,vrtspec_o,&
                 wrk2du,wrk2dv,&
                 gis_stochy%ls_node,gis_stochy%ls_nodes,gis_stochy%max_ls_nodes,&
                 gis_stochy%lats_nodes_a,gis_stochy%global_lats_a,gis_stochy%lonsperlat,&
                 gis_stochy%epsedn,gis_stochy%epsodn,gis_stochy%snnp1ev,gis_stochy%snnp1od,&
                 gis_stochy%plnev_a,gis_stochy%plnod_a,1)
       do i=1,lonf
          do j=1,gis_stochy%lats_node_a
             lat=gis_stochy%global_lats_a(ipt_lats_node_a-1+j)
             workgu(i,lat) = workgu(i,lat) + wrk2du(i,j,1)
             workgv(i,lat) = workgv(i,lat) + wrk2dv(i,j,1)
          enddo
       enddo
    enddo
    call mpp_sum(workgu,lonf*latg)
    call mpp_sum(workgv,lonf*latg)
! interpolate to cube grid
    do j=1,gis_stochy%ny
       pattern_1d = 0
       associate( tlats=>gis_stochy%parent_lats(1:gis_stochy%len(j),j),&
                  tlons=>gis_stochy%parent_lons(1:gis_stochy%len(j),j))
       call stochy_la2ga(workgu,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                        pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
       skebu_save(:,j,k)=pattern_1d(:)
       call stochy_la2ga(workgv,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                        pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
       skebv_save(:,j,k)=-1*pattern_1d(:)
       end associate
    enddo
  enddo
 endif
 do k=1,skeblevs-1
    skebu_save(:,:,k)=skebu_save(:,:,k+1)
    skebv_save(:,:,k)=skebv_save(:,:,k+1)
    do n=1,npatterns
       rpattern(n)%spec_e(:,:,k)=rpattern(n)%spec_e(:,:,k+1)
       rpattern(n)%spec_o(:,:,k)=rpattern(n)%spec_o(:,:,k+1)
    enddo
 enddo

! get pattern for last level
 workgu = 0.
 workgv = 0.
 do n=1,npatterns
!    if (stochini.AND. first_call) then
!       print*,'skipping advance'
!     else
       print*,'advancing'
       call patterngenerator_advance(rpattern(n),skeblevs,first_call)
!    endif
!   ke norm (convert streamfunction forcing to vorticity forcing)
    divspec_e = 0; divspec_o = 0.
    do nn=1,2
       vrtspec_e(:,nn,1) = gis_stochy%kenorm_e*rpattern(n)%spec_e(:,nn,skeblevs)
       vrtspec_o(:,nn,1) = gis_stochy%kenorm_o*rpattern(n)%spec_o(:,nn,skeblevs)
    enddo
  ! convert to winds
    call vrtdivspect_to_uvgrid(&
           divspec_e,divspec_o,vrtspec_e,vrtspec_o,&
           wrk2du,wrk2dv,&
           gis_stochy%ls_node,gis_stochy%ls_nodes,gis_stochy%max_ls_nodes,&
           gis_stochy%lats_nodes_a,gis_stochy%global_lats_a,gis_stochy%lonsperlat,&
           gis_stochy%epsedn,gis_stochy%epsodn,gis_stochy%snnp1ev,gis_stochy%snnp1od,&
           gis_stochy%plnev_a,gis_stochy%plnod_a,1)
    do i=1,lonf
       do j=1,gis_stochy%lats_node_a
          lat=gis_stochy%global_lats_a(ipt_lats_node_a-1+j)
          workgu(i,lat) = workgu(i,lat) + wrk2du(i,j,1)
          workgv(i,lat) = workgv(i,lat) + wrk2dv(i,j,1)
       enddo
    enddo
 enddo
 call mpp_sum(workgu,lonf*latg)
 call mpp_sum(workgv,lonf*latg)
! interpolate to cube grid
 do j=1,gis_stochy%ny
    pattern_1d = 0
    associate( tlats=>gis_stochy%parent_lats(:,j),&
               tlons=>gis_stochy%parent_lons(:,j))
    call stochy_la2ga(workgu,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                      pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
    skebu_save(:,j,skeblevs)=pattern_1d(:)
    call stochy_la2ga(workgv,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                      pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
    skebv_save(:,j,skeblevs)=-1*pattern_1d(:)
    end associate
  enddo
  deallocate(workgu)
  deallocate(workgv)
! interpolate in the vertical  ! consider moving to cubed sphere side,  more memory, but less interpolations
 print*,'looping through levs'
 do k=1,levs
    do j=1,gis_stochy%ny
       upattern_3d(:,j,k) = skeb_vwts(k,1)*skebu_save(:,j,skeb_vpts(k,1))+skeb_vwts(k,2)*skebu_save(:,j,skeb_vpts(k,2))
       vpattern_3d(:,j,k) = skeb_vwts(k,1)*skebv_save(:,j,skeb_vpts(k,1))+skeb_vwts(k,2)*skebv_save(:,j,skeb_vpts(k,2))
    enddo
 enddo
 first_call=.false.

end subroutine get_random_pattern_vector   

!>@brief The subroutine 'get_random_pattern_scalar' converts spherical harmonics to the gaussian grid then interpolates to the target grid
!>@details This subroutine is for a 2-D (lat-lon) scalar field
subroutine get_random_pattern_scalar(rpattern,npatterns,&
           gis_stochy,pattern_2d)

! generate a random pattern for stochastic physics
 implicit none
 type(random_pattern), intent(inout)  :: rpattern(npatterns)
 type(stochy_internal_state)          :: gis_stochy
 integer,intent(in)::   npatterns

 integer i,j,lat,n,pe_print
 real(kind=kind_dbl_prec), dimension(lonf,gis_stochy%lats_node_a,1):: wrk2d

! logical lprint

 real(kind=kind_dbl_prec), allocatable, dimension(:,:) :: workg
 real (kind=kind_dbl_prec)   glolal(lonf,gis_stochy%lats_node_a)
 integer kmsk0(lonf,gis_stochy%lats_node_a)
 real(kind=kind_dbl_prec) :: pattern_2d(gis_stochy%nx,gis_stochy%ny)
 real(kind=kind_dbl_prec) :: pattern_1d(gis_stochy%nx)
 pe_print=111
 !if (mpp_pe()==pe_print) then
 !   print*,'in get_random_pattern_scalar',npatterns
 !   print*,'nx,ny',gis_stochy%nx,gis_stochy%ny
 !endif

 kmsk0 = 0
 glolal = 0.
 do n=1,npatterns
    call patterngenerator_advance(rpattern(n),1,.false.)
    call scalarspect_to_gaugrid(                       &
         rpattern(n)%spec_e,rpattern(n)%spec_o,wrk2d,&
         gis_stochy%ls_node,gis_stochy%ls_nodes,gis_stochy%max_ls_nodes,&
         gis_stochy%lats_nodes_a,gis_stochy%global_lats_a,gis_stochy%lonsperlat,&
         gis_stochy%plnev_a,gis_stochy%plnod_a,1)
    glolal = glolal + wrk2d(:,:,1)
 enddo

 allocate(workg(lonf,latg))
 workg = 0.
  do j=1,gis_stochy%lats_node_a
     lat=gis_stochy%global_lats_a(ipt_lats_node_a-1+j)
     do i=1,lonf
        workg(i,lat) = glolal(i,j)
     enddo
  enddo

   call mpp_sum(workg,lonf*latg)
   !if (mpp_pe()==pe_print) then
   !   print*,'before interpolate',maxval(workg),minval(workg)
   !endif

! interpolate to cube grid
 !     if (mpp_pe()==pe_print) then
 !           write(34) real(workg,kind=4)
 !     endif
   do j=1,gis_stochy%ny
      pattern_1d = 0
      associate( tlats=>gis_stochy%parent_lats(1:gis_stochy%len(j),j),&
                 tlons=>gis_stochy%parent_lons(1:gis_stochy%len(j),j))
   !   if (mpp_pe()==pe_print) then
   !      if (j.EQ.1.OR.j.EQ.gis_stochy%ny) then
   !      !if (j.GT.100) then
   !      do i=1,gis_stochy%len(j) 
   !         print*,'lat,lon',j,gis_stochy%parent_lats(i,j),gis_stochy%parent_lons(i,j)
   !         print*,'lats',j,minval(gis_stochy%parent_lats(1:gis_stochy%len(j),j)),maxval(gis_stochy%parent_lats(1:gis_stochy%len(j),j))
   !         print*,'lons',j,minval(gis_stochy%parent_lons(1:gis_stochy%len(j),j)),maxval(gis_stochy%parent_lons(1:gis_stochy%len(j),j))
   !         !print*,'lonf,latg',lonf,latg,gis_stochy%nx,gis_stochy%ny,gis_stochy%len(j)
   !      enddo
   !      endif
   !   endif
      call stochy_la2ga(workg,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                        pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
      pattern_2d(:,j)=pattern_1d(:)
   !   if (mpp_pe()==pe_print) then
   !      print*,'gg_lats=',gg_lats(1),gg_lats(latg)
   !      print*,'after interpolate',j,maxval(pattern_1d),minval(pattern_1d)
   !      if (j.GT.110) then
   !      do i=1,gis_stochy%len(j) 
   !         print*,'lat,lon',j,gis_stochy%parent_lats(i,j),gis_stochy%parent_lons(i,j),pattern_1d(i)
   !      enddo
   !      endif
   !   endif
      end associate
   enddo
   deallocate(workg)
   !if (mpp_pe()==pe_print) then
   !   print*,'outside loop',mpp_pe(),maxval(pattern_2d),minval(pattern_2d)
   !endif

end subroutine get_random_pattern_scalar
!

!>@brief The subroutine 'scalarspect_to_gaugrid' converts scalar spherical harmonics to a scalar on a gaussian grid
!>@details This subroutine is for a 2-D (lat-lon) scalar field
subroutine scalarspect_to_gaugrid(&
           trie_ls,trio_ls,datag,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlat,&
           plnev_a,plnod_a,nlevs)
!\callgraph


      implicit none
      real(kind=kind_dbl_prec), intent(in) :: trie_ls(len_trie_ls,2,nlevs)
      real(kind=kind_dbl_prec), intent(in) :: trio_ls(len_trio_ls,2,nlevs)
      real(kind=kind_dbl_prec),  intent(out) :: datag(lonf,lats_node_a,nlevs)
      integer, intent(in) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),&
        nlevs,max_ls_nodes(nodes),lats_nodes_a(nodes),global_lats_a(latg),lonsperlat(latg)
      real(kind=kind_dbl_prec),intent(in) :: plnev_a(len_trie_ls,latg2),plnod_a(len_trio_ls,latg2)
! local vars
      real(kind=kind_dbl_prec) for_gr_a_1(lon_dim_a,nlevs,lats_dim_a)
      real(kind=kind_dbl_prec) for_gr_a_2(lonf,nlevs,lats_dim_a)
      integer              i,k
      integer              lan,lat
      integer              lons_lat

      call sumfln_stochy(trie_ls,&
                  trio_ls,&
                  lat1s_a,&
                  plnev_a,plnod_a,&
                  nlevs,ls_node,latg2,&
                  lats_dim_a,nlevs,for_gr_a_1,&
                  ls_nodes,max_ls_nodes,&
                  lats_nodes_a,global_lats_a,&
                  lats_node_a,ipt_lats_node_a,&
                  lonsperlat,lon_dim_a,latg,0)

      do lan=1,lats_node_a
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)
         CALL FOUR_TO_GRID(for_gr_a_1(1,1,lan),for_gr_a_2(1,1,lan),&
                           lon_dim_a,lonf,lons_lat,nlevs)
      enddo

      datag = 0.
      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        do k=1,nlevs
          do i=1,lons_lat
            datag(i,lan,k) = for_gr_a_2(i,k,lan)
          enddo
        enddo
      enddo

      return
      end subroutine scalarspect_to_gaugrid

!>@brief The subroutine 'dump_patterns' writes out the speherical harmonics to a file, controlled by FHSTOCH
!>@details Only the active patterns are written out
subroutine dump_patterns(sfile)
!\callgraph
    implicit none
    character*120 :: sfile
    integer :: stochlun,k,n
    stochlun=99
    if (mpp_pe()==mpp_root_pe()) then
       if (nsppt > 0 .OR. nshum > 0 .OR. nskeb > 0) then
          OPEN(stochlun,file=sfile,form='unformatted')
          print*,'open ',sfile,' for output'
       endif
    endif
    if (nsppt > 0) then
       do n=1,nsppt
       call write_pattern(rpattern_sppt(n),1,stochlun)
       enddo
    endif
    if (nshum > 0) then
       do n=1,nshum
       call write_pattern(rpattern_shum(n),1,stochlun)
       enddo
    endif
    if (nskeb > 0) then
       do n=1,nskeb
       do k=1,skeblevs
          call write_pattern(rpattern_skeb(n),k,stochlun)
       enddo
       enddo
    endif
    close(stochlun)
 end subroutine dump_patterns
!>@brief The subroutine 'write_patterns' writes out a single pattern and the seed associated with the random number sequence
!>@details Spherical harminoncs are stored with trianglular truncation
 subroutine write_pattern(rpattern,lev,lunptn)
!\callgraph
   implicit none
   type(random_pattern), intent(inout) :: rpattern
   integer, intent(in) :: lunptn,lev
   real(kind_dbl_prec), allocatable  :: pattern2d(:)
   integer nm,nn,arrlen,isize
   integer,allocatable :: isave(:)
   arrlen=2*ndimspec

   allocate(pattern2d(arrlen))
   pattern2d=0.0
   ! fill in apprpriate pieces of array
   !print*,'before collection...',me,maxval(rpattern%spec_e),maxval(rpattern%spec_o) &
   ! ,minval(rpattern%spec_e),minval(rpattern%spec_o)
   do nn=1,len_trie_ls
      nm = rpattern%idx_e(nn)
      if (nm == 0) cycle
      pattern2d(nm)          = rpattern%spec_e(nn,1,lev)
      pattern2d(ndimspec+nm) = rpattern%spec_e(nn,2,lev)
   enddo
   do nn=1,len_trio_ls
      nm = rpattern%idx_o(nn)
      if (nm == 0) cycle
      pattern2d(nm)          = rpattern%spec_o(nn,1,lev)
      pattern2d(ndimspec+nm) = rpattern%spec_o(nn,2,lev)
   enddo
   call mpp_sum(pattern2d,arrlen)
  !  write only on root process
   if (mpp_pe()==mpp_root_pe()) then
      print*,'writing out random pattern (min/max/size)',&
      minval(pattern2d),maxval(pattern2d),size(pattern2d)
      !print*,'max/min pattern=',maxval(pattern2d),minval(pattern2d)
      write(lunptn) ntrunc
      call random_seed(size=isize) ! get seed size
      allocate(isave(isize)) ! get seed
      call random_seed(get=isave,stat=rpattern%rstate) ! write seed
      write(lunptn) isave
      write(lunptn) pattern2d
   endif
   deallocate(pattern2d)
 end subroutine write_pattern
!>@brief The subroutine 'vrtdivspect_to_uvgrid' converts vorticty and divergence spherical harmonics to 
! zonal and meridional winds on the gaussian grid
!>@details This subroutine is for a 2-D (lat-lon) vector field
 subroutine vrtdivspect_to_uvgrid(&
           trie_di,trio_di,trie_ze,trio_ze,&
           uug,vvg,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           epsedn,epsodn,snnp1ev,snnp1od,plnev_a,plnod_a,nlevs)
!\callgraph

      implicit none
      real(kind=kind_dbl_prec), intent(in) :: trie_di(len_trie_ls,2,nlevs)
      real(kind=kind_dbl_prec), intent(in) :: trio_di(len_trio_ls,2,nlevs)
      real(kind=kind_dbl_prec), intent(in) :: trie_ze(len_trie_ls,2,nlevs)
      real(kind=kind_dbl_prec), intent(in) :: trio_ze(len_trio_ls,2,nlevs)
      real(kind=kind_dbl_prec),  intent(out) :: uug(lonf,lats_node_a,nlevs)
      real(kind=kind_dbl_prec),  intent(out) :: vvg(lonf,lats_node_a,nlevs)
      integer, intent(in) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),&
        nlevs,max_ls_nodes(nodes),lats_nodes_a(nodes),global_lats_a(latg),lonsperlar(latg)
      real(kind=kind_dbl_prec),intent(in) ::  epsedn(len_trie_ls),&
       epsodn(len_trio_ls),snnp1ev(len_trie_ls),snnp1od(len_trio_ls),&
       plnev_a(len_trie_ls,latg2),plnod_a(len_trio_ls,latg2)
! local vars
      real(kind=kind_dbl_prec) trie_ls(len_trie_ls,2,2*nlevs)
      real(kind=kind_dbl_prec) trio_ls(len_trio_ls,2,2*nlevs)
      real(kind=kind_dbl_prec) for_gr_a_1(lon_dim_a,2*nlevs,lats_dim_a)
      real(kind=kind_dbl_prec) for_gr_a_2(lonf,2*nlevs,lats_dim_a)
      integer              i,k
      integer              lan,lat
      integer              lons_lat
      real (kind=kind_dbl_prec) tx1

      do k=1,nlevs
        call dezouv_stochy(trie_di(1,1,k),       trio_ze(1,1,k),&
                    trie_ls(1,1,k), trio_ls(1,1,nlevs+k),&
                    epsedn,epsodn,snnp1ev,snnp1od,ls_node)
        call dozeuv_stochy(trio_di(1,1,k),       trie_ze(1,1,k),&
                    trio_ls(1,1,k), trie_ls(1,1,nlevs+k),&
                    epsedn,epsodn,snnp1ev,snnp1od,ls_node)
      enddo

      call sumfln_stochy(trie_ls,&
                  trio_ls,&
                  lat1s_a,&
                  plnev_a,plnod_a,&
                  2*nlevs,ls_node,latg2,&
                  lats_dim_a,2*nlevs,for_gr_a_1,&
                  ls_nodes,max_ls_nodes,&
                  lats_nodes_a,global_lats_a,&
                  lats_node_a,ipt_lats_node_a,&
                  lonsperlar,lon_dim_a,latg,0)

      do lan=1,lats_node_a
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlar(lat)
         CALL FOUR_TO_GRID(for_gr_a_1(1,1,lan),for_gr_a_2(1,1,lan),&
                           lon_dim_a,lonf,lons_lat,2*nlevs)
      enddo

      uug = 0.; vvg = 0.
      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlar(lat)
        tx1      = 1. / coslat_a(lat)
        do k=1,nlevs
          do i=1,lons_lat
            uug(i,lan,k) = for_gr_a_2(i,k,lan) * tx1
            vvg(i,lan,k) = for_gr_a_2(i,nlevs+k,lan) * tx1
          enddo
        enddo
      enddo

      return
 end subroutine vrtdivspect_to_uvgrid
end module get_stochy_pattern_mod
