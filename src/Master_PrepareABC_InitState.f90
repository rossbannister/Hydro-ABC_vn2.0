PROGRAM Master_PrepareABC_InitState

!*****************************************************
!*   Code to prepare the initial ABC model slice     *
!*                                                   *
!*   R. Petrie,    2.0:  10-06-2011                  *
!*   R. Petrie,    3.0:  30-07-2013                  *
!*   R. Bannister, 3.1:  30-07-2013                  *
!*   R. Bannister, 1.4da 20-10-2017                  *
!*   J. Zhu,       2.0   2022                        *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    UM_type,                     &
    dims_type,                   &
    ABC_type,                    &
    Init_ABC_opt,                &
    latitude,                    &
    Regular_vert_grid,           &
    datadirUM,                   &
    datadirABC_out,              &
    init_um_file,                &
    init_ABC_file,               &
    Adv_tracer,                  &
    gravity_wave_switch,         &
    nlongs,                      &
    nlevs,                       &
    x_scale,                     &
    z_scale,                     &
    source_x,                    &
    source_z,                    &
    press_amp,                   &
    buoy_amp,                    &
    q_amp,                       &   ! Hydro
    Moist_On,                    &   ! Hydro
    IQ_option,                   &   ! Hydro
    ExtraDiagFields


IMPLICIT NONE

INCLUDE "Boundaries.interface"


! Declare variables
!==========================
TYPE(UM_type)            :: um_data
TYPE(dims_type)          :: dims
TYPE(ABC_type)           :: ABC_data

CHARACTER(LEN=320)       :: UMfile, ABCfile

INTEGER                  :: i, k
REAL(ZREAL8)             :: x_scale2, z_scale2, dist_x, dist_x2, dist_z, dist_z2, bneg
REAL(ZREAL8), ALLOCATABLE :: qs(:)  !!! Hydro


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_PrepareABC_InitState'
PRINT*, '*************************************************************************'

! Read namelist
CALL SetOptions

ABCfile = TRIM(datadirABC_out) // '/' // TRIM(init_ABC_file)

CALL Initialise_dims (dims)
CALL Initialise_model_vars (ABC_data, .FALSE.)


!Create initial conditions from a slice of UM data
! ---------------------------------------------------------------
! This still needs to be done even if initialising with a prescribed blob
! (to get dimension data)

! Set derived parameters, and do initialisation

UMfile  = TRIM(datadirUM)  // '/' // TRIM(init_um_file)

CALL Initialise_um_data (um_data)

! Read in raw UM data, store in init_um_data
PRINT*, 'Reading UM data ...'
CALL Read_um_data_2d (um_data, UMfile, latitude)
PRINT*, '--done'

! Process um data
! Set grid
PRINT*, 'Setting Grid ...'
CALL Set_grid (um_data, dims, Regular_vert_grid)
CALL Set_ht_dep_cons (dims)
PRINT*, '--done'

IF ((Init_ABC_opt == 1) .OR. (Init_ABC_opt == 3) .OR. (Init_ABC_opt == 5) .OR. (Init_ABC_opt == 7)) THEN
  ! Define variables for simplified model from UM data, store in init_state
  PRINT*, 'Processing UM data to make it compatible with simplified model ...'
  CALL Process_um_data (um_data, ABC_data, dims)
ELSE
  ! Only setting pressure or buoyancy blobs, so reset fields to zero
  CALL Initialise_model_vars (ABC_data, .FALSE.)
END IF

IF ((Init_ABC_opt == 2) .OR. (Init_ABC_opt == 3) .OR. (Init_ABC_opt == 6) .OR. (Init_ABC_opt == 7)) THEN
  !Create initial conditions as a pressure perturbation
  ! ---------------------------------------------------------------
  x_scale2 = REAL(x_scale * x_scale)
  z_scale2 = REAL(z_scale * z_scale)
  DO i = 1, nlongs
    dist_x  = source_x - i
    dist_x2 = REAL(dist_x * dist_x)
    DO k = 1, nlevs
      dist_z  = source_z - k
      dist_z2 = REAL(dist_z * dist_z)
      ABC_data % r(i,k)   = ABC_data % r(i,k) + press_amp * EXP ( -1. * (dist_x2/x_scale2 + dist_z2/z_scale2) )
      ABC_data % rho(i,k) = 1.0 + ABC_data % r(i,k)
    END DO
  END DO
  CALL Boundaries (ABC_data, set_r=.TRUE., set_rho=.TRUE.)
END IF


IF ((Init_ABC_opt == 4) .OR. (Init_ABC_opt == 5) .OR. (Init_ABC_opt == 6) .OR. (Init_ABC_opt == 7)) THEN
  !Create initial conditions as a buoyancy perturbation
  ! ---------------------------------------------------------------
  x_scale2 = REAL(x_scale * x_scale)
  z_scale2 = REAL(z_scale * z_scale)
  DO i = 1, nlongs
    dist_x  = source_x - i
    dist_x2 = REAL(dist_x * dist_x)
    DO k = 1, nlevs
      dist_z  = source_z - k
      dist_z2 = REAL(dist_z * dist_z)
      ABC_data % b(i,k)   = ABC_data % b(i,k) + buoy_amp * EXP ( -1. * (dist_x2/x_scale2 + dist_z2/z_scale2) )
    END DO
  END DO
  CALL Boundaries (ABC_data, set_b=.TRUE.)
END IF


IF (Init_ABC_opt == 8) THEN
  ! Find the most negative value of b' (call this value bneg), and add 3*|bneg| to all b' values
  bneg = MINVAL(ABC_data % b)
  PRINT *, "The minimum value of the b' field is", bneg
  IF (bneg < 0.) THEN
    ! Add 3*bneg to the entire b' field
    PRINT *, "Making entire b' field positive"
    bneg = 3. * DABS(bneg)
    ABC_data % b(0:nlongs+1, 0:nlevs+1) = ABC_data % b(0:nlongs+1, 0:nlevs+1) + bneg
  ELSE
    PRINT *, "The b' field is non-negative so nothing done"
  END IF
END IF

! Set q by Hydro
IF (Moist_On .eqv. .TRUE.) THEN
 
  ! Initiate p0, rho0, saturated q profiles
  CALL  Init_p0rho0(dims) 

  ALLOCATE(qs(1:nlevs))
  CALL  Init_qs(dims, qs)  !! assuming all perturbations are zero

  x_scale2 = REAL(x_scale * x_scale)
  z_scale2 = REAL(z_scale * z_scale)
  
  SELECT CASE ( IQ_option )
      
    CASE(1)  !!!! add vapor as a bubble, q_amp is used as relative huimidity
      DO i = 1, nlongs

        dist_x  = abs(source_x - i)

        if( dist_x .gt. 0.5 * nlongs ) then
          dist_x = nlongs - dist_x   !!! make it round
        end if

        dist_x2 = REAL(dist_x * dist_x)

        DO k = 1, nlevs
          dist_z  = source_z - k
          dist_z2 = REAL(dist_z * dist_z)
          ABC_data % q(i,k) = qs(k) * q_amp * EXP ( -1. * (dist_x2/x_scale2 + dist_z2/z_scale2) )
        END DO
      END DO
  
    CASE(2)  !!!! add moisture as even layers, q_amp is used as mixing ratio
      DO k = 1, nlevs
        dist_z  = k - 1
        dist_z2 = REAL(dist_z * dist_z)
        ABC_data % q(:,k) = q_amp * EXP ( -1. * dist_z2 / z_scale2 )
      END DO
  
    CASE(3)  !!!! add moisture as even layers, q_amp is used as relative huimidity
  
      DO k = 1, nlevs
        dist_z  = source_z - k
        dist_z2 = REAL(dist_z * dist_z)
        ABC_data % q(:,k) = qs(k) * q_amp * EXP ( -1. * dist_z2 / z_scale2 )
      END DO
  
  END SELECT

  CALL Boundaries (ABC_data, set_q=.TRUE., set_qc=.TRUE.)

  DEALLOCATE(qs)

END IF  
! End Set q by Hydro

! Set the tracer
! --------------
IF (Adv_tracer) THEN
  PRINT*, 'Setting the tracer ...'
  DO i = 1, 4
    DO k = 1, 5
      ABC_data % tracer(INT(REAL(nlongs*i)/5.0)-2:INT(REAL(nlongs*i)/5.0)+2, k*10) = 1.0
      WRITE (*,*) 'Tracer at position',  i, INT(REAL(nlongs*i)/5.0)
    END DO
  END DO
  CALL Boundaries (ABC_data, set_tracer=.TRUE.)
END IF


! Add some gravity wave noise by setting u=0
! ------------------------------------------
IF (gravity_wave_switch) ABC_data % u(0:nlongs+1, 0:nlevs+1) = 0.


! Calculate diagnostics
!------------------------
CALL Calc_hydro(ABC_data, dims)
CALL Calc_geost(ABC_data)
CALL Energy(ABC_data, dims)
CALL Effective_buoyancy(ABC_data, dims)
CALL Calc_vert_mom_source(ABC_data, dims)
CALL Calc_horiz_div(ABC_data, dims)
CALL Calc_horiz_vort(ABC_data, dims)

! Output the result
CALL Write_state_2d (ABCfile, ABC_data, dims, 1, 0, 0, .not.ExtraDiagFields)  !! changed core_only, Hydro



! Clean up
CALL Deallocate_dims (dims)
CALL Deallocate_model_vars (ABC_data)
CALL Deallocate_um_data (um_data)

PRINT*,'-------------------------------------------------------------------------------'
PRINT*,'    Netcdf inital conds file: '
PRINT*,'    xconv -i ', ABCfile, ' &'
PRINT*,'-------------------------------------------------------------------------------'!


END PROGRAM Master_PrepareABC_InitState
