! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
! Purpose: Real function to do linear interpolation.
!
!          Called from UKCA_INTERPJ.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_flupj_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_FLUPJ_MOD'

CONTAINS

FUNCTION ukca_flupj(x,xs,ys,ln)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ukca_error_mod,  ONLY: maxlen_message, i_error_method_abort, error_report

IMPLICIT NONE


INTEGER, INTENT(IN) :: ln         ! Length of input arrays xs an

REAL, INTENT(IN) :: x             ! Input value of point
REAL, INTENT(IN) :: xs(ln)        ! Input array of points
REAL, INTENT(IN) :: ys(ln)        ! Input array of values at poi

!       Local variables

INTEGER :: errcode                ! Variable passed to ereport
INTEGER ::  i                     ! Loop variable

REAL :: delta_xs                  ! Variable used in interpolati
REAL :: temp_xs(ln)               ! Temporary array of sx
REAL :: temp_ys(ln)               ! Temporary array of ys
REAL :: ukca_flupj                ! Interpolated value at point

LOGICAL :: L_in_range             ! Flag to indicate if x is wit

INTEGER, POINTER :: err_code_ptr
CHARACTER(LEN=maxlen_message) :: cmessage     ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_FLUPJ'


!       Check if xs array is increasing or decreasing

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (xs(1) < xs(2)) THEN
  temp_ys(1:ln) = ys(1:ln)
  temp_xs(1:ln) = xs(1:ln)
ELSE
  DO i = ln,1,-1
    temp_xs(ln-i+1) = xs(i)
    temp_ys(ln-i+1) = ys(i)
  END DO
END IF

!       Do linear interpolation

L_in_range = .FALSE.
DO i = 2,ln
  IF (temp_xs(i-1) <= x .AND. x <= temp_xs(i)) THEN
    delta_xs = temp_xs(i) - temp_xs(i-1)
    ukca_flupj = (x - temp_xs(i-1))                                            &
               * temp_ys(i)/delta_xs                                           &
               + (temp_xs(i) - x  )                                            &
               * temp_ys(i-1)/delta_xs
    L_in_range = .TRUE.
    EXIT
  END IF
END DO

IF (.NOT. (L_in_range)) THEN
  err_code_ptr = 123
  WRITE(cmessage, '(A,3(1x,E12.5))') 'Value of x not within range of xs', x,   &
    MINVAL(xs), MAXVAL(xs)
  CALL error_report(i_error_method_abort, err_code_ptr, cmessage, RoutineName)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION ukca_flupj
END MODULE ukca_flupj_mod

