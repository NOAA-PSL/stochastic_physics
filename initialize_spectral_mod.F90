!>@brief The module 'initialize_spectral_mod' cotains the subroutine initialize_spectral
! !module: stochy_initialize_spectral
!          --- initialize module of the
!              gridded component of the stochastic physics patteern
!              generator, which is in spectral space
!
! !description: gfs dynamics gridded component initialize module.
!
! !revision history:
!
!  oct 11  2016 P.Pegion   copy of gsm/dynamics to create stand alone version
!
! !interface:
!
      module initialize_spectral_mod
!
!!uses:
!
      use kinddef
      use spectral_layout_mod,      only : len_trie_ls,len_trio_ls &
                                      ,ls_max_node,ls_dim,lat1s_a
      use stochy_internal_state_mod, only : stochy_internal_state
      use spectral_layout_mod,only:jcap,lon_dims_a,wgt_a,sinlat_a,coslat_a,colrad_a,rcs2_a,lats_nodes_h,global_lats_h,&
                                   latg,latg2,lonf
      use stochy_namelist_def
      use mpp_mod, only : mpp_pe,mpp_root_pe
      use getcon_spectral_mod, only: getcon_spectral
      use get_ls_node_stochy_mod, only: get_ls_node_stochy
      use getcon_lag_stochy_mod, only: getcon_lag_stochy
      !use mpp_mod
#ifndef IBM
      USE omp_lib
#endif

      implicit none

      contains

!>@brief The subroutine 'initialize_spectral' initializes the
!gridded component of the stochastic physics pattern
!>@details This code is taken from the legacy spectral GFS
      subroutine initialize_spectral(gis_stochy, rc)

! this subroutine set up the internal state variables,
! allocate internal state arrays for initializing the gfs system.
!----------------------------------------------------------------
!
      implicit none
!
!      type(stochy_internal_state), pointer, intent(inout) :: gis_stochy
      type(stochy_internal_state), intent(inout) :: gis_stochy
      integer,                                    intent(out)   :: rc
      integer           :: npe_single_member
      integer           :: i, l, locl

!-------------------------------------------------------------------

! set up gfs internal state dimension and values for dynamics etc
!-------------------------------------------------------------------

      npe_single_member = gis_stochy%npe_single_member

      gis_stochy%lon_dim_a = lon_s + 2
      jcap=ntrunc
      latg   = lat_s
      latg2  = latg/2
      lonf   = lon_s

      allocate(lat1s_a(0:jcap))
      allocate(lon_dims_a(latg))

      allocate(wgt_a(latg2))
      allocate(rcs2_a(latg2))

      ls_dim = (jcap)/gis_stochy%nodes+1
      allocate(gis_stochy%lonsperlat(latg))

      gis_stochy%lonsperlat(:)=lonf
!!
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      allocate (      gis_stochy%ls_node (ls_dim,3) )
      allocate (      gis_stochy%ls_nodes(ls_dim,gis_stochy%nodes) )
      allocate (  gis_stochy%max_ls_nodes(gis_stochy%nodes) )
!
      allocate (  gis_stochy%lats_nodes_a(gis_stochy%nodes) )
      allocate ( gis_stochy%global_lats_a(latg) )
!

! internal parallel structure.   Weiyu.
!---------------------------------------------------
      ALLOCATE(gis_stochy%TRIE_LS_SIZE      (npe_single_member))
      ALLOCATE(gis_stochy%TRIO_LS_SIZE      (npe_single_member))
      ALLOCATE(gis_stochy%TRIEO_LS_SIZE     (npe_single_member))
      ALLOCATE(gis_stochy%LS_MAX_NODE_GLOBAL(npe_single_member))
      ALLOCATE(gis_stochy%LS_NODE_GLOBAL    (LS_DIM*3, npe_single_member))

      gis_stochy%LS_NODE_GLOBAL     = 0
      gis_stochy%LS_MAX_NODE_GLOBAL = 0
      gis_stochy%TRIEO_TOTAL_SIZE   = 0

      DO i = 1, npe_single_member
          CALL GET_LS_NODE_STOCHY(i-1, gis_stochy%LS_NODE_GLOBAL(1, i),               &
                            gis_stochy%LS_MAX_NODE_GLOBAL(i), gis_stochy%IPRINT)
          gis_stochy%TRIE_LS_SIZE(i) = 0
          gis_stochy%TRIO_LS_SIZE(i) = 0
          DO LOCL = 1, gis_stochy%LS_MAX_NODE_GLOBAL(i)
              gis_stochy%LS_NODE_GLOBAL(LOCL+  LS_DIM, i)   = gis_stochy%TRIE_LS_SIZE(i)
              gis_stochy%LS_NODE_GLOBAL(LOCL+  2*LS_DIM, i) = gis_stochy%TRIO_LS_SIZE(i)

              L = gis_stochy%LS_NODE_GLOBAL(LOCL, i)

              gis_stochy%TRIE_LS_SIZE(i) = gis_stochy%TRIE_LS_SIZE(i) + (JCAP+3-L)/2
              gis_stochy%TRIO_LS_SIZE(i) = gis_stochy%TRIO_LS_SIZE(i) + (JCAP+2-L)/2
          END DO
          gis_stochy%TRIEO_LS_SIZE(i) = gis_stochy%TRIE_LS_SIZE(i)  + gis_stochy%TRIO_LS_SIZE(i) + 3
          gis_stochy%TRIEO_TOTAL_SIZE = gis_stochy%TRIEO_TOTAL_SIZE + gis_stochy%TRIEO_LS_SIZE(i)
      END DO


!---------------------------------------------------
!
      gis_stochy%iprint = 0
      call get_ls_node_stochy( gis_stochy%mype, gis_stochy%ls_node(:,1), ls_max_node, gis_stochy%nodes)
!
!
      len_trie_ls = 0
      len_trio_ls = 0
      do locl=1,ls_max_node
         gis_stochy%ls_node(locl,2) = len_trie_ls
         gis_stochy%ls_node(locl,3) = len_trio_ls
         l = gis_stochy%ls_node(locl,1)
         len_trie_ls = len_trie_ls+(jcap+3-l)/2
         len_trio_ls = len_trio_ls+(jcap+2-l)/2
      enddo
!
      allocate ( gis_stochy%epse  (len_trie_ls) )
      allocate ( gis_stochy%epso  (len_trio_ls) )
      allocate ( gis_stochy%epsedn(len_trie_ls) )
      allocate ( gis_stochy%epsodn(len_trio_ls) )
      allocate ( gis_stochy%kenorm_e(len_trie_ls) )
      allocate ( gis_stochy%kenorm_o(len_trio_ls) )
!
      allocate ( gis_stochy%snnp1ev(len_trie_ls) )
      allocate ( gis_stochy%snnp1od(len_trio_ls) )
!
      allocate ( gis_stochy%plnev_a(len_trie_ls,latg2) )
      allocate ( gis_stochy%plnod_a(len_trio_ls,latg2) )
      allocate ( gis_stochy%plnew_a(len_trie_ls,latg2) )
      allocate ( gis_stochy%plnow_a(len_trio_ls,latg2) )

      allocate(colrad_a(latg2))
      allocate(sinlat_a(latg))
      allocate(coslat_a(latg))
!!
      call getcon_spectral(gis_stochy)
!
      gis_stochy%lats_node_a     = gis_stochy%lats_nodes_a(gis_stochy%mype+1)

        if (.not. allocated(lats_nodes_h))  allocate (lats_nodes_h(gis_stochy%nodes))
        if (.not. allocated(global_lats_h)) allocate (global_lats_h(latg+2*gis_stochy%yhalo*gis_stochy%nodes))
        call getcon_lag_stochy(gis_stochy, lats_nodes_h, global_lats_h)

      allocate ( gis_stochy%trie_ls (len_trie_ls,2,gis_stochy%lotls) )
      allocate ( gis_stochy%trio_ls (len_trio_ls,2,gis_stochy%lotls) )

      rc=0

      end subroutine initialize_spectral

end module initialize_spectral_mod
