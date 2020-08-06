module  standalone_stochy_module

use kinddef
implicit none
public
integer :: isc,jsc,iec,jec,isd,ied,jsd,jed

type GFS_diag_type
    real (kind=kind_phys), allocatable :: ca_out  (:)      !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca_deep  (:)     !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca_turb  (:)     !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca_shal  (:)     !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca_rad   (:)     !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca_micro (:)     !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca1   (:)  
    real (kind=kind_phys), allocatable :: ca2   (:)  
    real (kind=kind_phys), allocatable :: ca3   (:)  
end type GFS_diag_type
type GFS_control_type
   integer  :: levs,me,nx,ny
   integer,allocatable  :: blksz(:)        !< for explicit data blocking: block sizes of all blocks
   real(kind=kind_phys) :: dtp             !< physics timestep in seconds
   real(kind=kind_phys) :: phour           !< previous forecast hour
   real(kind=kind_phys) :: sppt_amp        !< amplitude of sppt (to go to cld scheme)
   integer              :: kdt             !< current forecast iteration
   logical  :: do_sppt,do_shum,do_skeb,do_sfcperts,use_zmtnblck,do_rnda
   integer  ::  skeb_npass,nsfcpert
   character(len=65) :: fn_nml                   !< namelist filename
   character(len=256),allocatable :: input_nml_file(:) !< character string containing full namelist
   real(kind=kind_phys) ::  pertz0(5),pertzt(5),pertshc(5),pertlai(5),pertvegf(5),pertalb(5)
 !---cellular automata control parameters
    integer              :: nca             !< number of independent cellular automata
    integer              :: nlives          !< cellular automata lifetime
    integer              :: ncells          !< cellular automata finer grid
    real(kind=kind_phys) :: nfracseed       !< cellular automata seed probability
    integer              :: nseed           !< cellular automata seed frequency
    logical              :: do_ca           !< cellular automata main switch
    logical              :: ca_sgs          !< switch for sgs ca
    logical              :: ca_global       !< switch for global ca
    logical              :: ca_smooth       !< switch for gaussian spatial filter
    logical              :: isppt_deep      !< switch for combination with isppt_deep. OBS! Switches off SPPT on other tendencies!
    logical              :: isppt_pbl
    logical              :: isppt_shal      !
    logical              :: pert_flux       !
    logical              :: pert_trigger    !
    integer              :: iseed_ca        !< seed for random number generation in ca scheme
    integer              :: ca_amplitude    !< seed for random number generation in ca scheme
    integer              :: nspinup         !< number of iterations to spin up the ca
    real(kind=kind_phys) :: nthresh         !< threshold used for perturbed vertical velocity
end type GFS_control_type

 type GFS_statein_type
   real (kind=kind_phys), allocatable :: pgr  (:)      !< surface pressure (Pa) real
   real (kind=kind_phys), allocatable :: qgrs (:,:,:)  !< layer mean tracer concentration
   real (kind=kind_phys), allocatable :: vvl  (:,:)    !< layer mean vertical velocity in pa/sec
   real (kind=kind_phys), allocatable :: prsl  (:,:)   !< model layer mean pressure Pa
end type GFS_statein_type


type GFS_init_type
   integer :: nlunit
   real(kind=kind_phys),allocatable :: ak(:),bk(:),xlon(:,:),xlat(:,:)
   integer,allocatable  :: blksz(:)        !< for explicit data blocking: block sizes of all blocks
end type GFS_init_type

type GFS_grid_type
    real (kind=kind_phys),allocatable :: xlat   (:)    !< grid latitude in radians, default to pi/2 ->
    real (kind=kind_phys),allocatable :: xlon   (:)    !< grid longitude in radians, default to pi/2 ->
end type GFS_grid_type

type GFS_coupling_type
    real (kind=kind_phys),allocatable :: shum_wts  (:,:)
    real (kind=kind_phys),allocatable :: sppt_wts  (:,:)
    real (kind=kind_phys),allocatable :: sppt_pattern(:)
    real (kind=kind_phys),allocatable :: skebu_wts (:,:)
    real (kind=kind_phys),allocatable :: skebv_wts (:,:)
    real (kind=kind_phys),allocatable :: sfc_wts   (:,:)
    real (kind=kind_phys), allocatable :: cape    (:)      !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca_out  (:)      !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca_deep  (:)     !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca_turb  (:)     !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca_shal  (:)     !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca_rad   (:)     !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca_micro (:)     !< cellular automata fraction
    real (kind=kind_phys), allocatable :: ca1   (:)  
    real (kind=kind_phys), allocatable :: ca2   (:)  
    real (kind=kind_phys), allocatable :: ca3   (:)  
    integer              :: nsfcpert=6                             !< number of sfc perturbations
end type GFS_coupling_type
end module standalone_stochy_module

