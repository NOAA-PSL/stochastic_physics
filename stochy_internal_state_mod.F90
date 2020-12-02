!>@brief The module 'stochy_internal_state_mod' contains the spherical harmonic definitions and arrays
!! describing the target gaussian grid
!
! !module: stochy_internal_state_mod
!                         --- internal state definition of the
!                             gridded component of the spectral random patterns
!
! !description:  define the spectral internal state used to
!                                             create the internal state.
!---------------------------------------------------------------------------
! !revision history:
!
!  Oct 11 2016     P Pegion port of gfs_dynamics_interal_state
!
! !interface:
!

      module stochy_internal_state_mod

!!uses:
!------
      use spectral_layout_mod


      implicit none
      private

! -----------------------------------------------
      type,public::stochy_internal_state		! start type define
! -----------------------------------------------

      integer                   :: nodes

!
      integer lonf,latg,lats_node_a_max

      integer npe_single_member

      character(16)                     ::  cfhour1
!jws
      integer                           ::  num_file
      character(32)        ,allocatable ::  filename_base(:)
      integer                           ::  ipt_lats_node_a
      integer                           ::  lats_node_a
      integer                           ::  me
!jwe

      integer                           ::  nblck,kdt
!      real                              ::  deltim

      integer              ,allocatable ::      lonsperlat (:)
      integer              ,allocatable ::      ls_node    (:)
      integer              ,allocatable ::      ls_nodes   (:, :)
      integer              ,allocatable ::  max_ls_nodes   (:)

      integer              ,allocatable ::  lats_nodes_a   (:)
      integer              ,allocatable ::  global_lats_a  (:)
      integer              ,allocatable ::  global_lats_h  (:)
      integer                           :: xhalo,yhalo

      integer              ,allocatable ::  lats_nodes_a_fix (:)

      real,allocatable ::        epse  (:)
      real,allocatable ::        epso  (:)
      real,allocatable ::        epsedn(:)
      real,allocatable ::        epsodn(:)
      real,allocatable ::        kenorm_e(:)
      real,allocatable ::        kenorm_o(:)

      real,allocatable ::       snnp1ev(:)
      real,allocatable ::       snnp1od(:)

      real,allocatable ::       plnev_a(:,:)
      real,allocatable ::       plnod_a(:,:)
      real,allocatable ::       pddev_a(:,:)
      real,allocatable ::       pddod_a(:,:)
      real,allocatable ::       plnew_a(:,:)
      real,allocatable ::       plnow_a(:,:)


      real,allocatable ::       trie_ls(:,:,:)
      real,allocatable ::       trio_ls(:,:,:)

      INTEGER                               :: TRIEO_TOTAL_SIZE
      INTEGER, ALLOCATABLE, DIMENSION(:)    :: TRIE_LS_SIZE
      INTEGER, ALLOCATABLE, DIMENSION(:)    :: TRIO_LS_SIZE
      INTEGER, ALLOCATABLE, DIMENSION(:)    :: TRIEO_LS_SIZE
      INTEGER, ALLOCATABLE, DIMENSION(:)    :: LS_MAX_NODE_GLOBAL
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: LS_NODE_GLOBAL


!

!!
      integer              init,jpt,node,ibmsign,lon_dim

      integer              lotls

!      integer              jdt,ksout,maxstp
!      integer              mdt,idt
!      integer              mods,n1,n2,ndgf,ndgi,nfiles,nflps
      integer              nlunit

      integer              iret,ierr,iprint,k,l,locl,n
      integer              lan,lat
      integer              nx,ny,nz
      integer, allocatable :: len(:)
      real, allocatable :: parent_lons(:,:),parent_lats(:,:)


!
! -----------------------------------------------------
      end type stochy_internal_state		! end type define
! -----------------------------------------------------

! this state is supported by c pointer not f90 pointer, thus
! need this wrap.
!-----------------------------------------------------------
      type stochy_wrap		! begin type define
          type (stochy_internal_state), pointer :: int_state
      end type stochy_wrap	! end type define

      end module stochy_internal_state_mod
