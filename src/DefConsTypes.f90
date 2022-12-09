!===================================================================================================
MODULE DefConsTypes

!*****************************************************
!*   Definition of all global constants              *
!*   and derived variable types                      *
!*                                                   *
!*   R. Petrie,    2.0:  10-06-2011                  *
!*   R. Petrie,    3.0:  30-07-2013                  *
!*   R. Bannister, 3.1:  30-07-2013                  *
!*   R. Bannister, 1.4da 20-10-2017                  *
!*   R. Bannister, 1.6da 03-11/2020                  *
!*   J. Zhu,       2.0   2022                        *
!*                                                   *
!*****************************************************

IMPLICIT NONE

!Definition of some new data types
!Define an integer that describes a unified double precision
!-----------------------------------------------------------
INTEGER, PARAMETER        :: ZREAL8         = SELECTED_REAL_KIND(15,307)
INTEGER, PARAMETER        :: wp             = KIND(1.0D0)
INTEGER, PARAMETER        :: Nensmax        = 50
INTEGER, PARAMETER        :: Nlatsmax       = 512
INTEGER, PARAMETER        :: Npointsmax     = 24
INTEGER, PARAMETER        :: maxbatches     = 100
REAL(ZREAL8)              :: ConditionFudge = 0.000000001
REAL(ZREAL8)              :: small          = 0.00000000001
REAL(ZREAL8)              :: zero           = 0.0
REAL(ZREAL8)              :: unity          = 1.0
INTEGER                   :: random_seed    = 0

! Model parameters
!-----------------
INTEGER                   :: nlongs = 360            ! total number of longitude points
INTEGER                   :: nlevs  = 60             ! total number of vertical levels
INTEGER                   :: ntimesteps              ! number of timesteps to integrate
REAL(ZREAL8)              :: dt     = 4.             ! full timestep
REAL(ZREAL8)              :: dx     = 1500.          ! dx grid resolution
REAL(ZREAL8)              :: H      = 14862.01       ! UM model height
REAL(ZREAL8)              :: deltat
REAL(ZREAL8)              :: Lx
REAL(ZREAL8)              :: dz

! Mathematical and physical constants
!------------------------------------
REAL(ZREAL8), PARAMETER   :: Rd    = 287.058         ! Gas constant for dry air
REAL(ZREAL8), PARAMETER   :: Re    = 6.371E6         ! Mean radius of the earth
REAL(ZREAL8), PARAMETER   :: Cp    = 1005.7          ! Specific heat at constant pressure
REAL(ZREAL8), PARAMETER   :: Cv    = 719.0           ! Specific heat at constant volume
REAL(ZREAL8), PARAMETER   :: g     = 9.81            ! Acceleration due to gravity
REAL(ZREAL8), PARAMETER   :: p00   = 1.0E5           ! Reference surface pressure 1000 hPa
REAL(ZREAL8), PARAMETER   :: rho0  = 1.225           ! Reference density (value may need to change)
REAL(ZREAL8), PARAMETER   :: rho00 = 1.225           ! Density of air (at surface in real atmos)
REAL(ZREAL8), PARAMETER   :: kappa = 0.286           ! Rd/Cp
REAL(ZREAL8), PARAMETER   :: pi    = 3.141592654     ! pi

! Tuneable model parameters
!--------------------------
REAL(ZREAL8)              :: f     = 1.0E-4          ! Coriolis Parameter
REAL(ZREAL8)              :: A     = 0.02            ! A is the buoyancy frequency
REAL(ZREAL8)              :: B     = 0.01            ! B accoustic wave speed modulator
REAL(ZREAL8)              :: C     = 1.0E5           ! Constant relating pressure and density
REAL(ZREAL8)              :: BoundSpread = 50.       ! No of grid points to spread boundary discontinuity info
REAL(ZREAL8)              :: theta_r     = 273.

! Useful constants
!-----------------
REAL(ZREAL8)              :: third, half
REAL(ZREAL8)              :: recippi
REAL(ZREAL8)              :: recipsqrt2
REAL(ZREAL8)              :: recipdx
REAL(ZREAL8)              :: recipdx2
REAL(ZREAL8)              :: recip2dx2
REAL(ZREAL8)              :: recip2dx
REAL(ZREAL8)              :: alpha_f
REAL(ZREAL8)              :: alpha_N
REAL(ZREAL8)              :: beta_f
REAL(ZREAL8)              :: beta_N
REAL(ZREAL8)              :: recip_alpha_f
REAL(ZREAL8)              :: recip_alpha_N
REAL(ZREAL8)              :: bdiva_f
REAL(ZREAL8)              :: bdiva_N
REAL(ZREAL8)              :: fourpi2
REAL(ZREAL8)              :: sr_nlongs
REAL(ZREAL8)              :: half_sr_nlongs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Moist parameters by Hydro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(ZREAL8)              :: tau = 1E4                 ! evaporation timescale
REAL(ZREAL8)              :: Lv = 2500.                ! latent heat constant
REAL(ZREAL8)              :: theta00 = 300             ! surface potential temperature
REAL(ZREAL8)              :: SH = 9000                 ! Scale height
REAL(ZREAL8), ALLOCATABLE :: p0_FL(:)                  ! hydrostatic pressure at full level
REAL(ZREAL8), ALLOCATABLE :: rho0_FL(:)                ! hydrostatic density at full level
REAL(ZREAL8)              :: CO_thresh = 0.            ! threshold for condensation
REAL(ZREAL8)              :: q_amp     = 0.            ! amplitude for initiation of q
INTEGER                   :: Ev_option = 1             ! options for evaporation (sign of b')
INTEGER                   :: IQ_option = 1             ! options for q initiation
LOGICAL                   :: Moist_On  = .TRUE.        ! Turn on Moisture process by Hydro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Moist parameters by Hydro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!								

!**************************************************************************************************
! Declare Compound Types
!**************************************************************************************************

!-----------------------------------------------------------------------
! To store raw UM data slice
TYPE UM_type
  REAL(ZREAL8), ALLOCATABLE :: longs_u(:)
  REAL(ZREAL8), ALLOCATABLE :: longs_v(:)
  REAL(ZREAL8), ALLOCATABLE :: half_levs(:)
  REAL(ZREAL8), ALLOCATABLE :: full_levs(:)
  REAL(ZREAL8), ALLOCATABLE :: u(:,:)
  REAL(ZREAL8), ALLOCATABLE :: v(:,:)
  REAL(ZREAL8), ALLOCATABLE :: w(:,:)
  REAL(ZREAL8), ALLOCATABLE :: density(:,:)
  REAL(ZREAL8), ALLOCATABLE :: theta(:,:)
  REAL(ZREAL8), ALLOCATABLE :: exner_pressure(:,:)
  REAL(ZREAL8), ALLOCATABLE :: orog_height(:)
END TYPE UM_type


!-----------------------------------------------------------------------
! To store information about the dimensions (axes)
TYPE dims_type
  REAL(ZREAL8), ALLOCATABLE :: longs_u(:)
  REAL(ZREAL8), ALLOCATABLE :: longs_v(:)
  REAL(ZREAL8), ALLOCATABLE :: half_levs(:)
  REAL(ZREAL8), ALLOCATABLE :: full_levs(:)
  ! Variables required for vertical interpolation
  REAL(ZREAL8), ALLOCATABLE :: a1(:)
  REAL(ZREAL8), ALLOCATABLE :: b1(:)
  REAL(ZREAL8), ALLOCATABLE :: a2(:)
  REAL(ZREAL8), ALLOCATABLE :: b2(:)
  REAL(ZREAL8), ALLOCATABLE :: recip_half_kp1_k(:)
  REAL(ZREAL8), ALLOCATABLE :: recip_half_k_km1(:)
  REAL(ZREAL8), ALLOCATABLE :: recip_full_kp1_k(:)
  REAL(ZREAL8), ALLOCATABLE :: recip_full_k_km1(:)
END TYPE dims_type


!-----------------------------------------------------------------------
! To store ABC model fields and some diagnostics (single time)
TYPE ABC_type
  ! Horizontal grid is Awakara C grid
  ! Vertical grid is Charney-Phillips
  REAL(ZREAL8), ALLOCATABLE :: u(:,:)               ! Zonal wind perturbation
  REAL(ZREAL8), ALLOCATABLE :: v(:,:)               ! Meridional wind perturbation
  REAL(ZREAL8), ALLOCATABLE :: w(:,:)               ! Vertical wind perturbation
  REAL(ZREAL8), ALLOCATABLE :: r(:,:)               ! rho density perturbation
  REAL(ZREAL8), ALLOCATABLE :: b(:,:)               ! buoyancy perturbation
  REAL(ZREAL8), ALLOCATABLE :: tracer(:,:)          ! Tracer
  REAL(ZREAL8), ALLOCATABLE :: rho(:,:)             ! rho full field
  REAL(ZREAL8), ALLOCATABLE :: b_ef(:,:)            ! Effective buoyancy
  REAL(ZREAL8), ALLOCATABLE :: hydro_imbal(:,:)     ! Hydrostatic imbalance
  REAL(ZREAL8), ALLOCATABLE :: geost_imbal(:,:)     ! Geostrophic imbalance
  REAL(ZREAL8), ALLOCATABLE :: vert_mom_source(:,:) ! Vertical momentum source
  REAL(ZREAL8), ALLOCATABLE :: horiz_div(:,:)       ! Horizontal divergence
  REAL(ZREAL8), ALLOCATABLE :: horiz_vort(:,:)      ! Horizontal divergence
  REAL(ZREAL8)              :: Kinetic_Energy
  REAL(ZREAL8)              :: Buoyant_Energy
  REAL(ZREAL8)              :: Elastic_Energy
  REAL(ZREAL8)              :: Total_Energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Moist parameters by Hydro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(ZREAL8), ALLOCATABLE :: q(:,:)               ! Scaled Absolute humidity, rho * q
  REAL(ZREAL8), ALLOCATABLE :: qc(:,:)              ! rho * qc
  REAL(ZREAL8), ALLOCATABLE :: RH(:,:)              ! relative humidity
  REAL(ZREAL8)              :: Latent_Energy        ! total latent heat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Moist parameters by Hydro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END TYPE ABC_type


!-----------------------------------------------------------------------
! Used with model integration scheme
TYPE Averages_type
  REAL(ZREAL8), ALLOCATABLE :: u_1(:,:)
  REAL(ZREAL8), ALLOCATABLE :: u_2(:,:)
  REAL(ZREAL8), ALLOCATABLE :: u_m(:,:)
  REAL(ZREAL8), ALLOCATABLE :: w_1(:,:)
  REAL(ZREAL8), ALLOCATABLE :: w_2(:,:)
  REAL(ZREAL8), ALLOCATABLE :: w_m(:,:)
END TYPE Averages_type


!-----------------------------------------------------------------------
! Scheme to store control variable fields (classic use commented)
TYPE CV_type
                                        ! CLASSIC CVT INTERPRETATON  NORMAL MODE INTERPRETATION
  REAL(ZREAL8), ALLOCATABLE :: v1(:,:)  ! streamfunction             velocity potential
  REAL(ZREAL8), ALLOCATABLE :: v2(:,:)  ! velocity potential         streamfunction
  REAL(ZREAL8), ALLOCATABLE :: v3(:,:)  ! (unbalanced) r             w
  REAL(ZREAL8), ALLOCATABLE :: v4(:,:)  ! (unbalanced) b             r
  REAL(ZREAL8), ALLOCATABLE :: v5(:,:)  ! (unbalanced) w             b
  REAL(ZREAL8), ALLOCATABLE :: v6(:,:)  ! tracer
END TYPE CV_type



!-----------------------------------------------------------------------
! Scheme to store the control variable transform data (CVT)
TYPE CVT_type

  ! Options for the transforms
  INTEGER      :: CVT_order          ! 1 = as original MetO
                                     ! 2 = reversed horiz/vert
                                     ! 3 = normal mode-based CVT
  INTEGER      :: CVT_param_opt_gb   ! 1 = analytical balance (geostrophic)
                                     ! 2 = statistical balance
                                     ! 3 = no geostrophic balance
  INTEGER      :: CVT_param_opt_hb   ! 1 = analytical balance (hydrostatic)
                                     ! 2 = statistical balance
                                     ! 3 = no hydrostatic balance
  INTEGER      :: CVT_param_opt_ab   ! 1 = analytical balance (anelastic for w)
                                     ! 2 = no anelastic balance
  INTEGER      :: CVT_param_opt_reg  ! 1 = use vertical regression of the gb r
                                     ! 2 = no vertical regression
  INTEGER      :: CVT_vert_opt_sym   ! 1 = non-symmetric transform
                                     ! 2 = symmetric transform

  INTEGER      :: CVT_stddev_opt     ! 1 = Stddev constant for each control variable
                                     ! 2 = Level dependent only
                                     ! 3 = Longitude and level dependent

  ! Data structures to hold the transforms
  ! Standard deviations of the 6 control parameters
  REAL(ZREAL8), ALLOCATABLE :: sigma1(:,:)
  REAL(ZREAL8), ALLOCATABLE :: sigma2(:,:)
  REAL(ZREAL8), ALLOCATABLE :: sigma3(:,:)
  REAL(ZREAL8), ALLOCATABLE :: sigma4(:,:)
  REAL(ZREAL8), ALLOCATABLE :: sigma5(:,:)
  REAL(ZREAL8), ALLOCATABLE :: sigma6(:,:)

  !! The following are used for CVT_order=1 or 2 only (conventional CVT)
  !! -------------------------------------------------------------------
  ! Vertical modes of the 6 control parameters
  REAL(ZREAL8), ALLOCATABLE :: VertMode1(:,:,:)
  REAL(ZREAL8), ALLOCATABLE :: VertMode2(:,:,:)
  REAL(ZREAL8), ALLOCATABLE :: VertMode3(:,:,:)
  REAL(ZREAL8), ALLOCATABLE :: VertMode4(:,:,:)
  REAL(ZREAL8), ALLOCATABLE :: VertMode5(:,:,:)
  REAL(ZREAL8), ALLOCATABLE :: VertMode6(:,:,:)
  ! Vertical eigenvalues of the 6 control parameters (these are actually the square-roots)
  REAL(ZREAL8), ALLOCATABLE :: VertEV1(:,:)
  REAL(ZREAL8), ALLOCATABLE :: VertEV2(:,:)
  REAL(ZREAL8), ALLOCATABLE :: VertEV3(:,:)
  REAL(ZREAL8), ALLOCATABLE :: VertEV4(:,:)
  REAL(ZREAL8), ALLOCATABLE :: VertEV5(:,:)
  REAL(ZREAL8), ALLOCATABLE :: VertEV6(:,:)
  ! Horizontal eigenvalues of the 6 control parameters (these are actually the square-roots)
  REAL(ZREAL8), ALLOCATABLE :: HorizEV1(:,:)
  REAL(ZREAL8), ALLOCATABLE :: HorizEV2(:,:)
  REAL(ZREAL8), ALLOCATABLE :: HorizEV3(:,:)
  REAL(ZREAL8), ALLOCATABLE :: HorizEV4(:,:)
  REAL(ZREAL8), ALLOCATABLE :: HorizEV5(:,:)
  REAL(ZREAL8), ALLOCATABLE :: HorizEV6(:,:)
  ! Regression data for balanced density
  REAL(ZREAL8), ALLOCATABLE :: Cov_rbalrbal(:,:)
  REAL(ZREAL8), ALLOCATABLE :: Cov_rtotrbal(:,:)
  REAL(ZREAL8), ALLOCATABLE :: Regression(:,:)
  REAL(ZREAL8)              :: LB_optimal_fac
  REAL(ZREAL8)              :: HB_optimal_fac

  !! The following are used or CVT_order=3 only (normal mode-based CVT)
  !! ------------------------------------------------------------------
  ! Storage of the vertical mode structures
  REAL(ZREAL8),    ALLOCATABLE :: sin_mh_lh(:,:)  ! For chi, psi
  REAL(ZREAL8),    ALLOCATABLE :: sin_m_l(:,:)    ! For w, bp
  REAL(ZREAL8),    ALLOCATABLE :: cos_m_lh(:,:)   ! For rp
  REAL(ZREAL8),    ALLOCATABLE :: sin_m_lh(:,:)
  REAL(ZREAL8),    ALLOCATABLE :: cos_m_l(:,:)
  ! Storage of the vertical overlap matrices
  REAL(ZREAL8),    ALLOCATABLE :: X(:,:)
  REAL(ZREAL8),    ALLOCATABLE :: X_rho(:,:)
  REAL(ZREAL8),    ALLOCATABLE :: W(:)            ! Diag matrix, so diagonal elements only
  REAL(ZREAL8),    ALLOCATABLE :: W_rho(:,:)
  REAL(ZREAL8),    ALLOCATABLE :: P(:,:)
  REAL(ZREAL8),    ALLOCATABLE :: P_chi(:,:)
  REAL(ZREAL8),    ALLOCATABLE :: P_w(:,:)
  ! Storage of eigenvalues and eigenvectors of system matrix (m>0)
  REAL(ZREAL8),    ALLOCATABLE :: NM_evecs(:,:,:)  ! (k, component, vec no)
  REAL(ZREAL8),    ALLOCATABLE :: NM_evals(:,:)    ! (k, vec no)
  ! Switch set if the standard deviations have been set (or not)
  LOGICAL                      :: NM_stddev_set

END TYPE CVT_type


!-----------------------------------------------------------------------
! Storing the observations
TYPE Obs_type
  ! Batch index (for later development, e.g. correlated obs errs)
  INTEGER      :: batch
  INTEGER      :: obnumber_thisfile

  ! Observation time
  !INTEGER      :: year
  !INTEGER      :: month
  !INTEGER      :: day
  !INTEGER      :: hour
  !INTEGER      :: min
  !INTEGER      :: sec
  ! Absoute (seconds)
  INTEGER      :: t

  ! Location
  REAL(ZREAL8) :: longitude_m
  REAL(ZREAL8) :: level_ht
  INTEGER      :: xbox_lower    ! Closest lower grid index (x)
  INTEGER      :: xbox_lower_ws ! Auxillary info needed for wind speed obs
  INTEGER      :: zbox_lower    ! Closest lower grid index (z)
  INTEGER      :: zbox_lower_ws ! Auxillary info needed for wind speed obs
  INTEGER      :: tstep_lower   ! Closest lower time index of corresponding model state sequence

  ! Observations
  LOGICAL      :: ob_ok
  INTEGER      :: ob_of_what    ! 1 (u), 2 (v), 3 (w), 4 (r), 5 (b), 6 (tracer)
                                ! 7 (horizontal wind speed)
                                ! 8 (total wind speed)
  LOGICAL      :: y_true_known
  REAL(ZREAL8) :: y_true        ! True value of ob
  REAL(ZREAL8) :: y             ! Observation value
  REAL(ZREAL8) :: stddev        ! Error standard deviation
  REAL(ZREAL8) :: variance      ! Error variance (stddev squared)
  REAL(ZREAL8) :: y_ref         ! Model observation (reference)
  REAL(ZREAL8) :: d             ! y - y_ref
  REAL(ZREAL8) :: deltay_m      ! Model observation increment
  REAL(ZREAL8) :: hxmy          ! deltay_m - d = y_ref + deltay_m - y
  REAL(ZREAL8) :: deltay_m_hat  ! d(JO)/d(deltay_m) = (R^-1) hxmy

  ! Pointer to the next observation record
  TYPE(Obs_type), POINTER :: next

END TYPE Obs_type


!-----------------------------------------------------------------------
! Specification of the observations for generation
TYPE ObsSpec_type
  ! Reference time for all observations
  INTEGER      :: year0
  INTEGER      :: month0
  INTEGER      :: day0
  INTEGER      :: hour0
  INTEGER      :: min0
  INTEGER      :: sec0

  INTEGER      :: NumBatches                  ! Number of obervation batches
  INTEGER      :: batch(1:maxbatches)         ! The batch numbers
  INTEGER      :: seconds(1:maxbatches)       ! The absolute time of this batch (seconds)
  INTEGER      :: ob_of_what(1:maxbatches)    ! What is observed (see same variable name in Obs_type)
  INTEGER      :: NumObs_long(1:maxbatches)   ! Number of observations in longitude direction
  INTEGER      :: NumObs_height(1:maxbatches) ! Number of observations in the height direction
  REAL(ZREAL8) :: long_min(1:maxbatches)      ! min longitude of obs patch
  REAL(ZREAL8) :: long_max(1:maxbatches)      ! max longitude of obs patch
  REAL(ZREAL8) :: height_min(1:maxbatches)    ! min height of obs patch
  REAL(ZREAL8) :: height_max(1:maxbatches)    ! max height of obs patch
  REAL(ZREAL8) :: stddev(1:maxbatches)        ! error standard dev
END TYPE ObsSpec_type





! Variables to do with preparing the ABC init state
! -------------------------------------------------
INTEGER                   :: Init_ABC_opt                 ! 1=Take UM data
                                                          ! 2=Pressure blob only
                                                          ! 3=Sum of 1 and 2 (UM + press blob)
                                                          ! 4=Buoyancy blob only
                                                          ! 5=Sum of 1 and 4 (UM + buoy blob)
                                                          ! 6=Sum of 2 and 4 (press + buoy blob)
                                                          ! 7=Sum of 1, 2, 4 (UM + press + buoy blob)
CHARACTER(LEN=256)        :: datadirUM=''                 ! Directory to do with UM data
CHARACTER(LEN=256)        :: init_um_file=''              ! Input UM filename
INTEGER                   :: latitude = 144               ! The latitude to be extracted
LOGICAL                   :: Regular_vert_grid = .TRUE.   ! Set to use regular vertical levels
LOGICAL                   :: gravity_wave_switch = .FALSE.! To set u=0 (simulate gws)
CHARACTER(LEN=256)        :: init_ABC_file=''             ! Initial ABC filename
INTEGER                   :: source_x = 180               ! Long pos of centre of press/buoy blob
INTEGER                   :: source_z = 30                ! Level pos of centre of press/buoy blob
INTEGER                   :: x_scale = 80                 ! Long size of press/buoy blob
INTEGER                   :: z_scale = 3                  ! Level size of press/buoy blob

! Variables to do with running the forward model
! ----------------------------------------------
CHARACTER(LEN=256)        :: datadirABC_in=''             ! Directory to do with input of simplified model data
CHARACTER(LEN=256)        :: datadirABC_out=''            ! Directory to do with output of simplified model data
CHARACTER(LEN=256)        :: output_ABC_file=''           ! Dump file
CHARACTER(LEN=256)        :: diagnostics_file=''          ! For diagnostics
REAL(ZREAL8)              :: runlength = 60.0             ! Runlength in seconds
INTEGER                   :: ndumps = 10                  ! number of dump times
LOGICAL                   :: convection_switch = .FALSE.  ! Set to ?
REAL(ZREAL8)              :: press_amp = 0.01             ! Amplitude of pressure blob
REAL(ZREAL8)              :: buoy_amp = 0.1               ! Amplitude of buoyancy blob
LOGICAL                   :: Adv_tracer = .FALSE.         ! Set to advect tracer in calculations
LOGICAL                   :: Lengthscale_diagnostics = .FALSE.
                                                          ! Set to do lengthscale diagnostics at final time
LOGICAL                   :: ExtraDiagFields = .FALSE.    ! Set to output extra disgnostic fields with model output

! Variables to do with the calibration of the CVT, etc.
! ----------------------------------------
INTEGER                   :: CalibRunStage = 1            ! Stage 1 is to convert UM to ABC forecasts
                                                          ! Stage 2 is to compute perturbations
                                                          ! Stage 3 is to determine the regression
                                                          ! Stage 4 is parameter transform
                                                          ! Stage 5 is calibration of spatial stats
LOGICAL                   :: CreateCVTFileOnly = .FALSE.  ! Set to use stage 1 calib only to create skeleton CVT file
INTEGER                   :: NEns = 50                    ! Number of ensembles (0=do not use ensembles)
CHARACTER(LEN=256)        :: EnsDirs(1:Nensmax)=''        ! Directories containing the ensembles
INTEGER                   :: NEnsMems = 24                ! Number of ensemble members

INTEGER                   :: NNMC = 0                     ! Number of NMC pairs (0=do not use NMC)
CHARACTER(LEN=256)        :: NMCDir = ''                  ! Directories containing the NMC pairs

INTEGER                   :: NCanadianQuick = 0           ! Number of Canadian quick covs ens members produced

INTEGER                   :: Nlats = 1                    ! Number of different lats to extract from
                                                          ! UM files
INTEGER                   :: latindex(1:Nlatsmax)         ! The latitude indices

CHARACTER(LEN=256)        :: datadirABCfcs=''             ! Directory of ABC model forecasts (stage 1 calibration)
CHARACTER(LEN=256)        :: datadirABCperts=''           ! Directory of ABC model forecasts (stage 2 calibration)
CHARACTER(LEN=256)        :: datadirRegression=''         ! Directory of regression diagnostics (stage 3 calibration)
CHARACTER(LEN=256)        :: datadirConParams=''          ! Directory of control parameters (stage 4 calibration)
CHARACTER(LEN=256)        :: datadirConNMs=''             ! Directory of control normal modes (stage 8 calibration)
CHARACTER(LEN=256)        :: datadirCVT=''                ! Directory of the CVT data (stage 5 calibration)
CHARACTER(LEN=256)        :: CVT_file=''                  ! Name of file containing CVT data

INTEGER                   :: VertSmoothPoints = 0         ! Number of points in vertical to average for standard dev.
INTEGER                   :: HorizSmoothPoints = 0        ! Number of points in horizontal to average for standard dev.
LOGICAL                   :: ForceCor = .FALSE.           ! Set to adjust transforms so that horiz. and vert. transforms imply correlations
LOGICAL                   :: LevMeanBalr = .TRUE.         ! Set to make level mean r balanced
INTEGER                   :: CVT_order = 1                ! 1=original MetO trans order, 2=reversed horiz/vert
INTEGER                   :: CVT_param_opt_gb = 1         ! Geo bal opt, 1=anal, 2=stat, 3=off
INTEGER                   :: CVT_param_opt_hb = 1         ! Hyd bal opt, 1=anal, 2=stat, 3=off
INTEGER                   :: CVT_param_opt_ab = 1         ! Anelbal opt, 1=anal, 2=off
INTEGER                   :: CVT_param_opt_reg = 1        ! Vert regress opt, 1=on, 2=off
INTEGER                   :: CVT_vert_opt_sym = 1         ! 1=non symm, 2=sym
INTEGER                   :: CVT_stddev_opt = 2           ! 1=stddev cons, 2=level dep, 3=long/lev dep
REAL(ZREAL8)              :: LB_optimal_fac = 1.0         ! Multiplies the analytic LB operator
REAL(ZREAL8)              :: HB_optimal_fac = 1.0         ! Multiplies the analytic HB operator


! Variables to do with testing the DA
! ----------------------------------------------
CHARACTER(LEN=256)        :: datadirTestDA=''             ! Directory to do with testing the DA components
CHARACTER(LEN=256)        :: LS_file=''                   ! Name of file containing LS
CHARACTER(LEN=256)        :: Pert_file=''                 ! Name of file containing pert data


! Variables to do with differences between LS files
! -------------------------------------------------
CHARACTER(LEN=256)        :: LS_file1=''                  ! Name of file containing LS 1
CHARACTER(LEN=256)        :: LS_file2=''                  ! Name of file containing LS 2


LOGICAL                   :: RunAdjTests_CVT = .FALSE.
LOGICAL                   :: RunAdjTests_obs = .FALSE.
LOGICAL                   :: RunInvTests     = .FALSE.


! Variables to do with generating a background state, and the observations
! ----------------------------------------------
INTEGER                   :: Generate_mode = 1            ! 1 = Generate file that specifies obs positions/times/etc
                                                          ! 2 = Generate synthetic observations
                                                          ! 3 = Generate synthetic background state
TYPE(ObsSpec_type)        :: ObsSpec                      ! Initialised in SetOptions
CHARACTER(LEN=256)        :: datadir_ObsSpec=''           ! Directory to place obs spec file
CHARACTER(LEN=256)        :: ObsSpec_file=''              ! Filename of observation specifications
CHARACTER(LEN=256)        :: datadir_Bg=''                ! Directory containing background
CHARACTER(LEN=256)        :: Bg_file=''                   ! Background state filename
REAL(ZREAL8)              :: Bg_inflation=1.0             ! Inflation factor for generation of bg state
REAL(ZREAL8)              :: Bg_fac_v1=1.0                ! As above, but just for v1
REAL(ZREAL8)              :: Bg_fac_v2=1.0                ! As above, but just for v2
REAL(ZREAL8)              :: Bg_fac_v3=1.0                ! As above, but just for v3
REAL(ZREAL8)              :: Bg_fac_v4=1.0                ! As above, but just for v4
REAL(ZREAL8)              :: Bg_fac_v5=1.0                ! As above, but just for v5
REAL(ZREAL8)              :: Bg_fac_v6=1.0                ! As above, but just for v6
INTEGER                   :: v_modes=0                    ! For NM scheme only.  Use only first v_modes in U_v_NM(_adj).
CHARACTER(LEN=256)        :: datadir_Obs=''               ! Directory containing observations
CHARACTER(LEN=256)        :: Obs_file=''                  ! Observations filename
REAL(ZREAL8)              :: dt_da = 60.0                 ! Time-step of the DA (seconds)


! Variables to do with implied/raw covariances
! ----------------------------------------------
INTEGER                   :: ImplCov_npoints = 0          ! Number of different source points to consider
INTEGER                   :: longindex(1:Npointsmax)      ! The longitude indices
INTEGER                   :: levindex(1:Npointsmax)       ! The level indices
CHARACTER(LEN=256)        :: datadirImpliedCov=''         ! Directory to output implied covariances
CHARACTER(LEN=256)        :: datadirRawCov=''             ! Directory to output raw covariances


! Variables to do with the DA
! ----------------------------------------------
CHARACTER(LEN=256)        :: datadirAnal=''               ! Directory for the analysis files
CHARACTER(LEN=256)        :: anal_file=''                 ! Analysis file
CHARACTER(LEN=256)        :: analinc_file=''              ! Analysis increment file
INTEGER                   :: t0 = 0                       ! Time of start of this DA cycle (seconds)
INTEGER                   :: Hybrid_opt = 1               ! 1 = standard B
                                                          ! 2 = pure EnVar
                                                          ! 3 = hybrid EnVar
                                                          ! 4 = reduced rank KF-type hybrid
INTEGER                   :: Vartype = 3                  ! 3 = 3DVar
                                                          ! 35 = 3DFGAT
                                                          ! 4 = 4DVar
INTEGER                   :: N_outerloops = 1             ! Number of outer loops
INTEGER                   :: N_innerloops_max = 10        ! Maximum number of inner loops
REAL(ZREAL8)              :: mu = 0.001                   ! Small number for perturbing search direction
REAL(ZREAL8)              :: minus_mu                     ! Negative of above
REAL(ZREAL8)              :: crit_inner = 0.01            ! Stopping criterion for inner loop

! Variables to do with FFTs
! -------------------------
INTEGER                   :: fft_worklen_x
LOGICAL                   :: fft_init_x = .FALSE.
REAL(ZREAL8), ALLOCATABLE :: fft_wsave_x(:)
REAL(ZREAL8), ALLOCATABLE :: fft_work_x(:)
INTEGER                   :: fft_worklen_z
LOGICAL                   :: fft_init_z = .FALSE.
REAL(ZREAL8), ALLOCATABLE :: fft_wsave_z(:)
REAL(ZREAL8), ALLOCATABLE :: fft_work_z(:)


! Variables to do with linear analysis of the model
! -------------------------------------------------
CHARACTER(LEN=256)        :: datadirLinearAnal=''         ! Separate directory for the linear analysis data




! Define namelist
! ---------------
NAMELIST / UserOptions /                                                           &
! --- Set initial state for ABC model ---
  Init_ABC_opt, datadirUM, init_um_file, latitude, Regular_vert_grid,              &
  gravity_wave_switch, f, A, B, C, BoundSpread, init_ABC_file, nlongs, nlevs,      &
  tau, theta00, SH, Moist_On, IQ_option, Lv, Ev_option, CO_thresh, q_amp,          &  ! Moist by Hydro
  
! --- Run the forward model ---
  datadirABC_in, datadirABC_out, output_ABC_file, diagnostics_file, dt, dx, H,     &
  runlength, ndumps,                                                               &
  convection_switch, source_x, source_z, x_scale, z_scale, press_amp, buoy_amp,    &
  Adv_tracer, Lengthscale_diagnostics, ExtraDiagFields,                            &
! --- Testing the DA ---
  datadirTestDA, RunAdjTests_CVT, RunAdjTests_obs, RunInvTests, LS_file, Pert_file,&
  LS_file1, LS_file2,                                                              &
! --- DA ---
  datadirCVT, CVT_file, Hybrid_opt, Vartype, datadirAnal, anal_file, analinc_file, &
  N_outerloops, N_innerloops_max, crit_inner,                                      &
! --- Linear analysis ---
  datadirLinearAnal,                                                               &
! --- Calibration, CVT, etc.
  CalibRunStage, CreateCVTFileOnly, NEns, EnsDirs, NEnsMems, NNMC, NMCDir,         &
  datadirConParams,                                                                &
  datadirABCfcs, datadirABCperts, datadirRegression, datadirConNMs, NCanadianQuick,&
  Nlats, latindex, VertSmoothPoints, HorizSmoothPoints, ForceCor, LevMeanBalr,     &
  CVT_order, CVT_param_opt_gb, CVT_param_opt_hb, CVT_param_opt_ab,                 &
  CVT_param_opt_reg, CVT_vert_opt_sym, CVT_stddev_opt,                             &
  LB_optimal_fac, HB_optimal_fac,                                                  &
! --- Generate background, obs, etc.
  Generate_mode, ObsSpec, datadir_ObsSpec, ObsSpec_file, datadir_Bg,               &
  datadir_Obs, Bg_file, Bg_inflation, Obs_file, dt_da, t0, random_seed,            &
  Bg_fac_v1, Bg_fac_v2, Bg_fac_v3, Bg_fac_v4, Bg_fac_v5, Bg_fac_v6, v_modes,       &
! Implied covs
  datadirImpliedCov, datadirRawCov, ImplCov_npoints, longindex, levindex


END MODULE DefConsTypes
