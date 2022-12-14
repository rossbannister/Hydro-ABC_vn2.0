SUBROUTINE Boundaries (state, set_u, set_v, set_w, set_r, set_b, set_rho, set_beff, set_tracer, set_q, set_qc)

!************************************
!* Subroutine to apply boundary     *
!* conditions to model variables    *
!*                                  *
!* Input flags to choose which      *
!* variables to apply boundary      *
!* conditions to                    *
!*                                  *
!* Set no flags to do them all      *
!************************************

USE DefConsTypes, ONLY :     &
    ABC_type,                &
    nlongs, nlevs,           &
    zero, Moist_On

IMPLICIT NONE

TYPE(ABC_type),        INTENT(INOUT) :: state
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_u
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_v
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_w
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_r
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_b
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_rho
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_beff
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_tracer
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_q     ! moist by Hydro
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_qc	  ! moist by Hydro										 
!Local variables
!---------------
LOGICAL                              :: flag_u, flag_v, flag_w, flag_r
LOGICAL                              :: flag_b, flag_rho, flag_beff, flag_tracer
LOGICAL                              :: flag_q, flag_qc		! moist by Hydro											   
LOGICAL                              :: some_flags

IF (PRESENT(set_u)) THEN
  flag_u = set_u
ELSE
  flag_u = .FALSE.
END IF
IF (PRESENT(set_v)) THEN
  flag_v = set_v
ELSE
  flag_v = .FALSE.
END IF
IF (PRESENT(set_w)) THEN
  flag_w = set_w
ELSE
  flag_w = .FALSE.
END IF
IF (PRESENT(set_r)) THEN
  flag_r = set_r
ELSE
  flag_r = .FALSE.
END IF
IF (PRESENT(set_b)) THEN
  flag_b = set_b
ELSE
  flag_b = .FALSE.
END IF
IF (PRESENT(set_rho)) THEN
  flag_rho = set_rho
ELSE
  flag_rho = .FALSE.
END IF
IF (PRESENT(set_beff)) THEN
  flag_beff = set_beff
ELSE
  flag_beff = .FALSE.
END IF
IF (PRESENT(set_tracer)) THEN
  flag_tracer = set_tracer
ELSE
  flag_tracer = .FALSE.
END IF
!! moist by Hydro
IF (PRESENT(set_q)) THEN
  flag_q = set_q
ELSE
  flag_q = .FALSE.
END IF
IF (PRESENT(set_qc)) THEN
  flag_qc = set_qc
ELSE
  flag_qc = .FALSE.
END IF	  
!! end moist by Hydro

some_flags = PRESENT(set_u) .OR. PRESENT(set_v) .OR. PRESENT(set_w) .OR. PRESENT(set_r) .OR. &
             PRESENT(set_b) .OR. PRESENT(set_rho) .OR. PRESENT(set_beff) .OR. PRESENT(set_tracer) .OR. &
             PRESENT(set_Q) .OR. PRESENT(set_QC)   ! Hydro												

IF (.NOT.some_flags) THEN
  ! No flags set.  This is shorthand for all flags
  flag_u      = .TRUE.
  flag_v      = .TRUE.
  flag_w      = .TRUE.
  flag_r      = .TRUE.
  flag_b      = .TRUE.
  flag_rho    = .TRUE.
  flag_beff   = .TRUE.
  flag_tracer = .TRUE.
  IF(Moist_On) THEN
    flag_q      = .TRUE.  ! Hydro
    flag_qc     = .TRUE.  ! Hydro
  END IF  
END IF

IF (flag_u) THEN
  ! Horizontal boundaries
  state % u(0, 1:nlevs)          = state % u(nlongs, 1:nlevs)
  state % u(nlongs+1, 1:nlevs)   = state % u(1, 1:nlevs)
  ! Vertical boundaries
  state % u(0:nlongs+1,0)        = - 1.0 * state % u(0:nlongs+1,1)
  state % u(0:nlongs+1,nlevs+1)  = state % u(0:nlongs+1, nlevs)
END IF

IF (flag_v) THEN
  ! Horizontal boundaries
  state % v(0, 1:nlevs)          = state % v(nlongs, 1:nlevs)
  state % v(nlongs+1, 1:nlevs)   = state % v(1, 1:nlevs)
  ! Vertical boundaries
  state % v(0:nlongs+1,0)        = - 1.0 * state % v(0:nlongs+1,1)
  state % v(0:nlongs+1,nlevs+1)  = state % v(0:nlongs+1, nlevs)
END IF

IF (flag_w) THEN
  ! Horizontal boundaries
  state % w(0, 1:nlevs)          = state % w(nlongs, 1:nlevs)
  state % w(nlongs+1, 1:nlevs)   = state % w(1, 1:nlevs)
  ! Vertical boundaries
  state % w(0:nlongs+1,0)        = zero
  state % w(0:nlongs+1,nlevs)    = zero
  state % w(0:nlongs+1,nlevs+1)  = zero
END IF

IF (flag_r) THEN
  ! Horizontal boundaries
  state % r(0, 1:nlevs)          = state % r(nlongs, 1:nlevs)
  state % r(nlongs+1, 1:nlevs)   = state % r(1, 1:nlevs)
  ! Vertical boundaries
  state % r(0:nlongs+1,0)        = state % r(0:nlongs+1,1)
  state % r(0:nlongs+1,nlevs+1)  = state % r(0:nlongs+1,nlevs)
END IF

IF (flag_b) THEN
  ! Horizontal boundaries
  state % b(0, 1:nlevs)          = state % b(nlongs, 1:nlevs)
  state % b(nlongs+1, 1:nlevs)   = state % b(1, 1:nlevs)
  ! Vertical boundaries
  state % b(0:nlongs+1,0)        = zero
  state % b(0:nlongs+1,nlevs)    = zero
  state % b(0:nlongs+1,nlevs+1)  = zero
END IF

IF (flag_rho) THEN
  ! Horizontal boundaries
  state % rho(0, 1:nlevs)          = state % rho(nlongs, 1:nlevs)
  state % rho(nlongs+1, 1:nlevs)   = state % rho(1, 1:nlevs)
  ! Vertical boundaries
  state % rho(0:nlongs+1,0)        = state % rho(0:nlongs+1,1)
  state % rho(0:nlongs+1,nlevs+1)  = state % rho(0:nlongs+1,nlevs)
END IF

IF (flag_beff) THEN
  ! Horizontal boundaries
  state % b_ef(0, 1:nlevs)          = state % b_ef(nlongs, 1:nlevs)
  state % b_ef(nlongs+1, 1:nlevs)   = state % b_ef(1, 1:nlevs)
  ! Vertical boundaries
  state % b_ef(0:nlongs+1,0)        = zero
  state % b_ef(0:nlongs+1,nlevs)    = zero
  state % b_ef(0:nlongs+1,nlevs+1)  = zero
END IF

IF (flag_tracer) THEN
  ! Horizontal boundaries
  state % tracer(0, 1:nlevs)          = state % tracer(nlongs, 1:nlevs)
  state % tracer(nlongs+1, 1:nlevs)   = state % tracer(1, 1:nlevs)
  ! Vertical boundaries
  state % tracer(0:nlongs+1, 0)       = state % tracer(0:nlongs+1, 1)
  state % tracer(0:nlongs+1, nlevs+1) = state % tracer(0:nlongs+1, nlevs)
END IF


!! moist by Hydro
IF (flag_q) THEN
  ! Horizontal boundaries
  state % q(0, 1:nlevs)          = state % q(nlongs, 1:nlevs)
  state % q(nlongs+1, 1:nlevs)   = state % q(1, 1:nlevs)
  ! Vertical boundaries
  state % q(0:nlongs+1, 0)       = state % q(0:nlongs+1, 1)
  state % q(0:nlongs+1, nlevs+1) = state % q(0:nlongs+1, nlevs)
END IF

IF (flag_qc) THEN
  ! Horizontal boundaries
  state % qc(0, 1:nlevs)          = state % qc(nlongs, 1:nlevs)
  state % qc(nlongs+1, 1:nlevs)   = state % qc(1, 1:nlevs)
  ! Vertical boundaries
  state % qc(0:nlongs+1, 0)       = state % qc(0:nlongs+1, 1)
  state % qc(0:nlongs+1, nlevs+1) = state % qc(0:nlongs+1, nlevs)
END IF
!! end of moist by Hydro
END SUBROUTINE Boundaries
