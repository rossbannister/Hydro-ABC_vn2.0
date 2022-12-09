SUBROUTINE Energy (state, dims)

! To calculate components of energy

USE DefConsTypes, ONLY :  &
    ABC_type,             &
    dims_type,            &
    ZREAL8,               &
    nlongs, nlevs,        &
    half,                 &
    A, B, C, dx, rho00, Lv, Moist_On

IMPLICIT NONE

! Declare subroutine parameters
TYPE(ABC_type),  INTENT(INOUT)   :: state
TYPE(dims_type), INTENT(IN)      :: dims

! Declare local variables
INTEGER                          :: x, z
REAL(ZREAL8)                     :: u, v, w, rho, bp, r, vol, mass0
REAL(ZREAL8)                     :: E_k, E_b, E_e
REAL(ZREAL8)                     :: q, qc, E_L

! Declare functions
REAL(ZREAL8)                  :: INT_FH

! Initialise
E_k = 0.0
E_e = 0.0
E_b = 0.0
E_L  = 0.0   !! Moist by Hydro


! Calculate energy
DO z = 1, nlevs
  vol  = dx * (dims % full_levs(z) - dims % full_levs(z-1))
  DO x = 1, nlongs
    ! Interpolate to rho points on the grid
    r     = state % r(x,z)
    u     = (state % u(x-1,z) + state % u(x,z)) * half
    v     = state % v(x,z)
    w     = INT_FH (state % w(x,z-1), state % w(x,z), z, dims)
    bp    = INT_FH (state % b(x,z-1), state % b(x,z), z, dims)
    rho   = state % rho(x,z)
    mass0 = vol * rho00

								 
    E_k  = E_k + mass0 * rho * (u*u + v*v + w*w) * half
    E_e  = E_e + mass0 * C * r*r * half / B
    E_b  = E_b + mass0 * rho * bp*bp * half / (A*A)


    IF ( Moist_On .eqv. .TRUE. )  E_L  = E_L + mass0 * rho * state % q(x,z) * Lv   !! Hydro

  END DO
END DO

state % Kinetic_Energy = E_k
state % Buoyant_Energy = E_b
state % Elastic_Energy = E_e

IF ( Moist_On .eqv. .FALSE. ) THEN
  state % Total_Energy   = E_k + E_e + E_b   !Hydro
ELSE
  state % Latent_Energy = E_L   !Hydro
  state % Total_Energy   = E_k + E_e + E_b + E_L   !Hydro
END IF

END SUBROUTINE Energy
