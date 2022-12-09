! ================================================================================================
SUBROUTINE RMS_diff_cvcv (ControlVar1, ControlVar2, unit)

! Compute the RMS of the difference between the two inputs
! Ross Bannister, March 2021

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  CV_type,                &
  nlongs,                 &
  nlevs


IMPLICIT NONE

TYPE(CV_type),  INTENT(IN) :: ControlVar1, ControlVar2
INTEGER,        INTENT(IN) :: unit

REAL(ZREAL8)               :: total, total1, total2, total3, total4, total5


total1 = SUM( (ControlVar1 % v1(1:nlongs, 1:nlevs) - ControlVar2 % v1(1:nlongs, 1:nlevs)) *   &
              (ControlVar1 % v1(1:nlongs, 1:nlevs) - ControlVar2 % v1(1:nlongs, 1:nlevs)) ) / &
              (REAL(nlongs+2) * REAL(nlevs))
total2 = SUM( (ControlVar1 % v2(1:nlongs, 1:nlevs) - ControlVar2 % v2(1:nlongs, 1:nlevs)) *   &
              (ControlVar1 % v2(1:nlongs, 1:nlevs) - ControlVar2 % v2(1:nlongs, 1:nlevs)) ) / &
              (REAL(nlongs+2) * REAL(nlevs))
total3 = SUM( (ControlVar1 % v3(1:nlongs, 1:nlevs) - ControlVar2 % v3(1:nlongs, 1:nlevs)) *   &
              (ControlVar1 % v3(1:nlongs, 1:nlevs) - ControlVar2 % v3(1:nlongs, 1:nlevs)) ) / &
              (REAL(nlongs+2) * REAL(nlevs))
total4 = SUM( (ControlVar1 % v4(1:nlongs, 1:nlevs) - ControlVar2 % v4(1:nlongs, 1:nlevs)) *   &
              (ControlVar1 % v4(1:nlongs, 1:nlevs) - ControlVar2 % v4(1:nlongs, 1:nlevs)) ) / &
              (REAL(nlongs+2) * REAL(nlevs))
total5 = SUM( (ControlVar1 % v5(1:nlongs, 1:nlevs) - ControlVar2 % v5(1:nlongs, 1:nlevs)) *   &
              (ControlVar1 % v5(1:nlongs, 1:nlevs) - ControlVar2 % v5(1:nlongs, 1:nlevs)) ) / &
              (REAL(nlongs+2) * REAL(nlevs))

WRITE (unit,*) '  MS diff field 1', total1
WRITE (unit,*) '  MS diff field 2', total2
WRITE (unit,*) '  MS diff field 3', total3
WRITE (unit,*) '  MS diff field 4', total4
WRITE (unit,*) '  MS diff field 5', total5
total = total1 + total2 + total3 + total4 + total5
WRITE (unit,*) '  RMS diff is    ', SQRT(total)

END SUBROUTINE RMS_diff_cvcv




! ================================================================================================
SUBROUTINE RMS_diff_mvmv (ModelVar1, ModelVar2, unit)

! Compute the RMS of the difference between the two inputs
! Ross Bannister, March 2021

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  nlongs,                 &
  nlevs


IMPLICIT NONE

TYPE(ABC_type), INTENT(IN) :: ModelVar1, ModelVar2
INTEGER,        INTENT(IN) :: unit

REAL(ZREAL8)               :: total, total1, total2, total3, total4, total5

total1 = SUM( (ModelVar1 % u(0:nlongs+1, 1:nlevs) - ModelVar2 % u(0:nlongs+1, 1:nlevs)) *   &
              (ModelVar1 % u(0:nlongs+1, 1:nlevs) - ModelVar2 % u(0:nlongs+1, 1:nlevs)) ) / &
              (REAL(nlongs+2) * REAL(nlevs))
total2 = SUM( (ModelVar1 % v(0:nlongs+1, 1:nlevs) - ModelVar2 % v(0:nlongs+1, 1:nlevs)) *   &
              (ModelVar1 % v(0:nlongs+1, 1:nlevs) - ModelVar2 % v(0:nlongs+1, 1:nlevs)) ) / &
              (REAL(nlongs+2) * REAL(nlevs))
total3 = SUM( (ModelVar1 % w(0:nlongs+1, 1:nlevs) - ModelVar2 % w(0:nlongs+1, 1:nlevs)) *   &
              (ModelVar1 % w(0:nlongs+1, 1:nlevs) - ModelVar2 % w(0:nlongs+1, 1:nlevs)) ) / &
              (REAL(nlongs+2) * REAL(nlevs))
total4 = SUM( (ModelVar1 % r(0:nlongs+1, 1:nlevs) - ModelVar2 % r(0:nlongs+1, 1:nlevs)) *   &
              (ModelVar1 % r(0:nlongs+1, 1:nlevs) - ModelVar2 % r(0:nlongs+1, 1:nlevs)) ) / &
              (REAL(nlongs+2) * REAL(nlevs))
total5 = SUM( (ModelVar1 % b(0:nlongs+1, 1:nlevs) - ModelVar2 % b(0:nlongs+1, 1:nlevs)) *   &
              (ModelVar1 % b(0:nlongs+1, 1:nlevs) - ModelVar2 % b(0:nlongs+1, 1:nlevs)) ) / &
              (REAL(nlongs+2) * REAL(nlevs))

WRITE (unit,*) '  MS diff field 1', total1
WRITE (unit,*) '  MS diff field 2', total2
WRITE (unit,*) '  MS diff field 3', total3
WRITE (unit,*) '  MS diff field 4', total4
WRITE (unit,*) '  MS diff field 5', total5
total = total1 + total2 + total3 + total4 + total5
WRITE (unit,*) '  RMS diff is    ', SQRT(total)

END SUBROUTINE RMS_diff_mvmv
