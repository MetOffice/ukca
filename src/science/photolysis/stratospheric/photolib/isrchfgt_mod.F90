! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE isrchfgt_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Given an array XX of length N, and given a value X, this returns a
!     value J such that X lies between XX(J) and XX(J+1). XX must be
!     monotonically  increasing. J=0 or J=N is returned
!     to indicate that X is out of range.

!     The table entry J is found by bisection.


!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ISRCHFGT_MOD'

CONTAINS
INTEGER FUNCTION isrchfgt(n,xx,x)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) :: n       ! Length of array XX
REAL, INTENT(IN)    :: xx(n)   ! Monotonic array of length N
REAL, INTENT(IN)    :: x       ! Value whose position in XX is required

! Local variables
INTEGER :: jl
INTEGER :: ju
INTEGER :: jm

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ISRCHFGT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Initialize upper & lower limits.
jl = 0
ju = n + 1

DO WHILE ( ju-jl > 1 )
  jm = (ju + jl) / 2
  IF ( (xx(n) > xx(1)) .EQV. (x > xx(jm)) ) THEN
    jl = jm
  ELSE
    ju = jm
  END IF
END DO
isrchfgt = jl + 1
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END FUNCTION isrchfgt
END MODULE isrchfgt_mod
