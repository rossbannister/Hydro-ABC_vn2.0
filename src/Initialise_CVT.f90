SUBROUTINE Initialise_CVT (CVT)

! Initialise control variable transform structure
! Does not alter the ABC parameters and options though

USE DefConsTypes, ONLY :     &
    CVT_type,                &
    nlongs, nlevs,           &
    ZREAL8,                  &
    pi

IMPLICIT NONE

! Subroutine parameters
TYPE(CVT_type), INTENT(INOUT)   :: CVT

! Local variables
INTEGER                         :: l, m, mp
REAL(ZREAL8)                    :: ld, ldh, md, mdh, mpd, mpdh, nlevd




IF (CVT % CVT_order < 3) THEN

  ! ===== CONVENTIONAL CVT WITH PARAMETER, HORIZ, AND VERTICAL TRANSFORMS =====

  ! Standard deviations of the 6 control parameters
  ! -----------------------------------------------
  IF (.NOT.ALLOCATED(CVT % sigma1)) ALLOCATE (CVT % sigma1(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma2)) ALLOCATE (CVT % sigma2(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma3)) ALLOCATE (CVT % sigma3(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma4)) ALLOCATE (CVT % sigma4(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma5)) ALLOCATE (CVT % sigma5(1:nlongs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma6)) ALLOCATE (CVT % sigma6(1:nlongs, 1:nlevs))
  CVT % sigma1(1:nlongs, 1:nlevs) = 0.0
  CVT % sigma2(1:nlongs, 1:nlevs) = 0.0
  CVT % sigma3(1:nlongs, 1:nlevs) = 0.0
  CVT % sigma4(1:nlongs, 1:nlevs) = 0.0
  CVT % sigma5(1:nlongs, 1:nlevs) = 0.0
  CVT % sigma6(1:nlongs, 1:nlevs) = 0.0

  ! Vertical modes of the 6 control parameters
  ! ------------------------------------------
  IF (.NOT.ALLOCATED(CVT % VertMode1)) ALLOCATE (CVT % VertMode1(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertMode2)) ALLOCATE (CVT % VertMode2(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertMode3)) ALLOCATE (CVT % VertMode3(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertMode4)) ALLOCATE (CVT % VertMode4(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertMode5)) ALLOCATE (CVT % VertMode5(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertMode6)) ALLOCATE (CVT % VertMode6(1:nlevs, 1:nlevs, 1:nlongs/2+1))
  CVT % VertMode1(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode2(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode3(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode4(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode5(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode6(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0

  ! Vertical eigenvalues of the 6 control parameters (these are actually the square-roots)
  ! --------------------------------------------------------------------------------------
  IF (.NOT.ALLOCATED(CVT % VertEV1)) ALLOCATE (CVT % VertEV1(1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertEV2)) ALLOCATE (CVT % VertEV2(1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertEV3)) ALLOCATE (CVT % VertEV3(1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertEV4)) ALLOCATE (CVT % VertEV4(1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertEV5)) ALLOCATE (CVT % VertEV5(1:nlevs, 1:nlongs/2+1))
  IF (.NOT.ALLOCATED(CVT % VertEV6)) ALLOCATE (CVT % VertEV6(1:nlevs, 1:nlongs/2+1))
  CVT % VertEV1(1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertEV2(1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertEV3(1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertEV4(1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertEV5(1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertEV6(1:nlevs, 1:nlongs/2+1) = 0.0

  ! Horizontal eigenvalues of the 6 control parameters (these are actually the square-roots)
  ! ----------------------------------------------------------------------------------------
  IF (.NOT.ALLOCATED(CVT % HorizEV1)) ALLOCATE (CVT % HorizEV1(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % HorizEV2)) ALLOCATE (CVT % HorizEV2(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % HorizEV3)) ALLOCATE (CVT % HorizEV3(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % HorizEV4)) ALLOCATE (CVT % HorizEV4(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % HorizEV5)) ALLOCATE (CVT % HorizEV5(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % HorizEV6)) ALLOCATE (CVT % HorizEV6(1:nlongs/2+1, 1:nlevs))
  CVT % HorizEV1(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % HorizEV2(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % HorizEV3(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % HorizEV4(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % HorizEV5(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % HorizEV6(1:nlongs/2+1, 1:nlevs) = 0.0

  ! Regression data for balanced density
  ! ------------------------------------
  IF (.NOT.ALLOCATED(CVT % Cov_rbalrbal)) ALLOCATE (CVT % Cov_rbalrbal(1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % Cov_rtotrbal)) ALLOCATE (CVT % Cov_rtotrbal(1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % Regression))   ALLOCATE (CVT % Regression(1:nlevs, 1:nlevs))
  CVT % Cov_rbalrbal(1:nlevs, 1:nlevs) = 0.0
  CVT % Cov_rtotrbal(1:nlevs, 1:nlevs) = 0.0
  CVT % Regression(1:nlevs, 1:nlevs)   = 0.0


ELSE

  ! ===== NORMAL MODE-BASED CVT =====

  ! Storage of the vertical mode structures
  IF (.NOT.ALLOCATED(CVT % sin_mh_lh)) ALLOCATE (CVT % sin_mh_lh(1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sin_m_l))   ALLOCATE (CVT % sin_m_l(1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % cos_m_lh))  ALLOCATE (CVT % cos_m_lh(1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sin_m_lh))  ALLOCATE (CVT % sin_m_lh(1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % cos_m_l))   ALLOCATE (CVT % cos_m_l(1:nlevs, 1:nlevs))

  nlevd = REAL(nlevs)
  DO l = 1, nlevs
    ld  = REAL(l)
    ldh = ld - 0.5
    DO m = 1, nlevs
      md  = REAL(m-1)
      mdh = md + 0.5
      CVT % sin_mh_lh(l,m) = SIN(mdh*pi*ldh/nlevd)  ! On half levels
      CVT % sin_m_l(l,m)   = SIN(md*pi*ld/nlevd)    ! On full levels
      CVT % cos_m_lh(l,m)  = COS(md*pi*ldh/nlevd)   ! On half levels
      CVT % sin_m_lh(l,m)  = SIN(md*pi*ldh/nlevd)   ! On half levels
      CVT % cos_m_l(l,m)   = COS(md*pi*ld/nlevd)    ! On full levels
    END DO
  END DO

  ! Storage of the vertical overlap matrices
  IF (.NOT.ALLOCATED(CVT % X))         ALLOCATE (CVT % X    (1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % X_rho))     ALLOCATE (CVT % X_rho(1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % W))         ALLOCATE (CVT % W    (1:nlevs))
  IF (.NOT.ALLOCATED(CVT % W_rho))     ALLOCATE (CVT % W_rho(1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % P))         ALLOCATE (CVT % P    (1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % P_chi))     ALLOCATE (CVT % P_chi(1:nlevs, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % P_w))       ALLOCATE (CVT % P_w  (1:nlevs, 1:nlevs))

  DO mp = 1, nlevs
    IF (mp == 0) THEN
      CVT % W(mp) = 0.
    ELSE
      CVT % W(mp) = nlevd/2.
    END IF
    DO m = 1, nlevs
      CVT % X(mp,m)     = SUM(CVT % sin_mh_lh(1:nlevs,mp) * CVT % sin_mh_lh(1:nlevs,m))
      CVT % X_rho(mp,m) = SUM(CVT % sin_mh_lh(1:nlevs,mp) * CVT % cos_m_lh(1:nlevs,m))
      CVT % W_rho(mp,m) = SUM(CVT % sin_m_l(1:nlevs,mp) * CVT % sin_m_lh(1:nlevs,m))
      CVT % P(mp,m)     = SUM(CVT % cos_m_lh(1:nlevs,mp) * CVT % cos_m_lh(1:nlevs,m))
      CVT % P_chi(mp,m) = SUM(CVT % cos_m_lh(1:nlevs,mp) * CVT % sin_mh_lh(1:nlevs,m))
      CVT % P_w(mp,m)   = SUM(CVT % cos_m_lh(1:nlevs,mp) * CVT % cos_m_l(1:nlevs,m))
    END DO
  END DO

  ! Standard deviations of the 5 control parameters.  These are the 5 types of normal mode
  IF (.NOT.ALLOCATED(CVT % sigma1)) ALLOCATE (CVT % sigma1(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma2)) ALLOCATE (CVT % sigma2(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma3)) ALLOCATE (CVT % sigma3(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma4)) ALLOCATE (CVT % sigma4(1:nlongs/2+1, 1:nlevs))
  IF (.NOT.ALLOCATED(CVT % sigma5)) ALLOCATE (CVT % sigma5(1:nlongs/2+1, 1:nlevs))
  CVT % sigma1(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % sigma2(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % sigma3(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % sigma4(1:nlongs/2+1, 1:nlevs) = 0.0
  CVT % sigma5(1:nlongs/2+1, 1:nlevs) = 0.0


  ! Storage of eigenvalues and eigenvectors of system matrix (m>0)
  IF (.NOT.ALLOCATED(CVT % NM_evecs))  ALLOCATE (CVT % NM_evecs(1:nlongs/2+1,   &  !k
                                                                1:5*(nlevs-1),  &  !component
                                                                1:5*(nlevs-1)))    !vector no
  IF (.NOT.ALLOCATED(CVT % NM_evals))  ALLOCATE (CVT % NM_evals(1:nlongs/2+1,   &  !k
                                                                1:5*(nlevs-1)))    !vector no

  ! The following will be set to .TRUE. when the standard deviations of the normal modes are set
  CVT % NM_stddev_set = .FALSE.

END IF


END SUBROUTINE Initialise_CVT
