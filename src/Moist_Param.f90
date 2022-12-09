!********************************************************************************
!*                                                                              *
!*  Moist parametrisation routines for hydro_ABC                                *
!*                                                                              *
!*                                                                              *
!*   J. Zhu, vn2.0, 2022                                                        *
!*                                                                              *
!********************************************************************************



!!!!!!!!!  Initialize parameterized p0 and rho0 at full levels
SUBROUTINE Init_p0rho0(dims)

USE DefConsTypes, ONLY :   &
    dims_type, nlongs, nlevs, &
    rho00, p00, p0_FL, rho0_FL, SH, g, ZREAL8

IMPLICIT NONE

TYPE(dims_type), INTENT(IN)    :: dims

INTEGER                        :: z

IF (.NOT.ALLOCATED(p0_FL))                   ALLOCATE(p0_FL(1:nlevs))
IF (.NOT.ALLOCATED(rho0_FL))                 ALLOCATE(rho0_FL(1:nlevs))


DO z = 1, nlevs
  rho0_FL(z) = rho00 * EXP( - dims % full_levs(z) / SH )
  p0_FL(z) = p00 - SH * g * ( rho00 - rho0_FL(z) )
ENDDO

END SUBROUTINE Init_p0rho0

!!!!!!!!!  Initialize qs for Master_PrepareABC_InitState, assumuing all perturbation 0.
SUBROUTINE Init_qs(dims, qs)

USE DefConsTypes, ONLY :   &
    dims_type, nlongs, nlevs, &
    theta_r, g, kappa, p00, theta00, p0_FL, ZREAL8 !, rho0_FL

IMPLICIT NONE

TYPE(dims_type), INTENT(IN)    :: dims
REAL(ZREAL8),    INTENT(OUT)   :: qs(1:nlevs)

INTEGER                        :: x, z

REAL(ZREAL8)                   :: p, temp       ! pressure, temperature, saturation specific humidity

!OPEN(UNIT = 100, FILE = "profile.txt")

DO z = 1, nlevs

    p = p0_FL(z)

    temp = ( theta00 + 1E-4 * theta_r / g * dims % full_levs(z)  ) *  &
           ( p / p00 ) ** kappa

    qs(z) = 380000 / p * EXP( 17.3*(temp-273.2)/(temp-35.9) )

    !write(100, '(I2 F12.3 F12.3 F12.3 F12.3 F12.3)')  z, dims % full_levs(z), p/100, temp, rho0_FL(z), qs(z)

ENDDO

!CLOSE(100)

END SUBROUTINE Init_qs


SUBROUTINE q_qs(q, f, df_dq, z, q0, b0, p0)

  USE DefConsTypes, ONLY :   &
    A, theta_r, g, kappa, p00, theta00,  ZREAL8, Lv

  REAL(ZREAL8), INTENT(IN)  ::  q, z, q0, b0, p0
  REAL(ZREAL8), INTENT(OUT) ::  f, df_dq
  
  b = sqrt( b0 * b0 - 2 * A*A  * Lv * (q - q0) )

  temp = ( theta00 + 1E-4 * theta_r / g * z + theta_r / g * b ) * ( p0 / p00 ) ** kappa

  qs = 380000. / p0 * exp( 17.3*(temp-273.2)/(temp-35.9) )

  f = q - qs
  
  dqs_dtemp = 156001020000. * exp(17.3*(temp-273.2)/(temp-35.9)) / (p0 * (temp-35.9)**2)

  dtemp_db = theta_r / g * ( p0 / p00 ) ** kappa

  db_dq = - A*A*Lv / b

  dqs_dq = dqs_dtemp * dtemp_db * db_dq

  df_dq = 1 - dqs_dq

END SUBROUTINE q_qs


FUNCTION find_q(qs, q0, tol, z, b0, p0)

  USE DefConsTypes, ONLY : ZREAL8

  IMPLICIT NONE

  REAL(ZREAL8), INTENT(IN) :: qs, z, q0, b0, p0
  REAL(ZREAL8), INTENT(IN) :: tol
  REAL(ZREAL8) :: find_q, f, df_dq, dx
  INTEGER, PARAMETER :: MAXIT = 200

  INTEGER :: j

  !find_q = 0.5 * (qs + q0) !Initial guess.
  find_q = q0 !Initial guess.

  do j=1, MAXIT

    call q_qs(find_q, f, df_dq, z, q0, b0, p0)

    dx = f / df_dq
    find_q = find_q - dx

    if ((qs - find_q) * (find_q - q0) < 0.0) then
      write(*,*) 'Warning, unsolvable find_q, b = ', b0, qs, find_q, q0
      find_q = -999.   ! unsolvable, set to negative
      write(*,*) 'Set find_q to ', find_q
      RETURN
    end if

    if (abs(dx) < tol) RETURN  !Convergence.

  end do

  stop 'find_q exceeded maximum iterations'

END FUNCTION find_q


!!! Moist Parameterization, first calculate CO. EV. rate, 
!!! then net condensation and latent heat, and update model state b, q, qc, RH
SUBROUTINE Moist(state, dims)

USE DefConsTypes, ONLY :   &
    ABC_type, dims_type, nlongs, nlevs, dt, tau, &
    A, C, theta_r, g, kappa, p00, theta00, p0_FL, rho0_FL, ZREAL8,  &
    Moist_On, Ev_option, CO_thresh, Lv

IMPLICIT NONE

TYPE(ABC_type),  INTENT(INOUT) :: state
TYPE(dims_type), INTENT(IN)    :: dims

INTEGER                        :: i, k   ! coordinates of grid point

REAL(ZREAL8)                   :: r_at_b
REAL(ZREAL8)                   :: p, temp, qs, q, qc, z, b, find_q, trial_q   ! pressure, temperature, saturation specific humidity
REAL(ZREAL8)                   :: INT_HF
REAL(ZREAL8)                   :: dq, delta_LE, b2, b2_new

DO k = 1, nlevs
  DO i = 1,nlongs

    z = dims % full_levs(k)

    b = state % b(i,k) 

    r_at_b = INT_HF( state % r(i,k), state % r(i,k+1), k, Dims)
    
    p = p0_FL(k) + C * rho0_FL(k) * r_at_b
    
    temp = ( theta00 + 1E-4 * theta_r / g * z + theta_r / g * b ) *  &
           ( p / p00 ) ** kappa
    
    qs = 380000. / p * EXP( 17.3*(temp-273.2)/(temp-35.9) )
    
    q  = state % q(i,k)
    qc = state % qc(i,k)
    
    state % RH(i,k) = q / qs   !! update relative humidity for model output

    !! calculate trial dq
    IF ( q .gt. qs ) THEN

      !! only condense when b is postive, or negative but exceeding threshold
      IF ( b .ge. 0. ) THEN          !! if exceeding threshold let it condense

        trial_q = find_q( qs, q, 1d-1, z, b, p )

        if(trial_q .ge. 0.) then 

          dq = trial_q - q

        else

          dq = qs - q

        end if
 
        delta_LE = 2. * A * A * Lv * dq 
        b2 = b*b
        b2_new = b2 - delta_LE
  
        state % b (i,k)  = sqrt(b2_new)
        state % q (i,k)  = state % q (i,k)  + dq
        state % qc (i,k) = state % qc (i,k) - dq

      ELSE  ! b is negative

        dq = qs - q

        delta_LE = 2. * A * A * Lv * dq 
        b2 = b*b
        b2_new = b2 - delta_LE

        IF ( -delta_LE/b2 .ge. CO_thresh ) THEN

          state % b (i,k) = sqrt(b2_new)
          state % q (i,k)  = state % q (i,k)  + dq
          state % qc (i,k) = state % qc (i,k) - dq

        END IF

      END IF

    ELSE IF ( ( qc .gt. 0. ) .and. (b .ge. 0.) ) THEN   !! evaporation, only when b greater than 0

      dq = min((qs - q) * (1-exp(-dt/tau)), qc)

      delta_LE = 2. * A * A * Lv * dq 
      b2 = b*b
      b2_new = b2 - delta_LE

      if ( delta_LE .le. b2 ) then  !! b energy is enough for all evaporation
        state % b (i,k) = sqrt(b2_new)
      else         !! b energy is not enough, use all b energy, recalculate dq
        state % b (i,k) = 0.
        dq = b2 / (2. * A * A * Lv)
      end if

      state % q (i,k)  = state % q (i,k)  + dq
      state % qc (i,k) = state % qc (i,k) - dq

    END IF

  END DO
END DO

END SUBROUTINE Moist
