! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE ei1_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE umPrintMgr, ONLY: umMessage, umPrint
IMPLICIT NONE


! Description:
!     This function calculates the first exponential integral using a
!     series expansion.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EI1_MOD'

CONTAINS
REAL FUNCTION ei1(x)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Function interface
REAL, INTENT(IN) :: x  ! Argument of the first exponential integral

! Local variables
REAL, PARAMETER :: thresh=1.0e-35
REAL, PARAMETER :: xmax=50.0

REAL, PARAMETER :: a011= -0.57721566
REAL, PARAMETER :: a012=  0.99999193
REAL, PARAMETER :: a013= -0.24991055
REAL, PARAMETER :: a014=  0.05519968
REAL, PARAMETER :: a015= -0.00976004
REAL, PARAMETER :: a016=  0.00107857

REAL, PARAMETER :: a1101= 8.5733285301
REAL, PARAMETER :: a1102=18.0590169520
REAL, PARAMETER :: a1103= 8.6347608925
REAL, PARAMETER :: a1104= 0.2677737343

REAL, PARAMETER :: b1101= 9.5733223454
REAL, PARAMETER :: b1102=25.6329561486
REAL, PARAMETER :: b1103=21.0996530827
REAL, PARAMETER :: b1104= 3.9584969228

REAL, PARAMETER :: a10inf1=4.03640
REAL, PARAMETER :: a10inf2=1.15198
REAL, PARAMETER :: b10inf1=5.03637
REAL, PARAMETER :: b10inf2=4.19160

REAL :: top
REAL :: bot
REAL :: x2
REAL :: x3
REAL :: x4

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EI1'



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF ( x < 0.0 ) THEN
  WRITE(umMessage,'(1X,A,1X,E12.5)')                                           &
  '** EI1 WARNING ** Negative argument >',x
  CALL umPrint(umMessage,src='ei1_mod')
END IF

!     Check for minimum value
IF ( x <= thresh ) THEN
  ei1 = 86.9210129
ELSE IF ( (x > thresh) .AND. (x <= 1.0) ) THEN
  ei1 = a011                                                                   &
    + a012* x                                                                  &
    + a013*(x**(2))                                                            &
    + a014*(x**(3))                                                            &
    + a015*(x**(4))                                                            &
    + a016*(x**(5))                                                            &
    - LOG(x)
ELSE IF ( (x > 1.0) .AND. (x <= 10.0) ) THEN
  x2 = x*x
  x3 = x2*x
  x4 = x3*x
  top = x4 + a1101*x3 + a1102*x2 + a1103*x + a1104
  bot = x4 + b1101*x3 + b1102*x2 + b1103*x + b1104
  ei1 = top/(bot*x*EXP(x))
ELSE IF ( (x > 10.0) .AND. (x < xmax) ) THEN
  x2 = x*x
  top = x2 + a10inf1*x + a10inf2
  bot = x2 + b10inf1*x + b10inf2
  ei1 = top/(bot*x*EXP(x))
ELSE IF ( x >= xmax ) THEN
  ei1 = 0.0
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END FUNCTION ei1
END MODULE ei1_mod
