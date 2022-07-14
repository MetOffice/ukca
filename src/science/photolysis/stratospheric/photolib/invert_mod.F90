! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE invert_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Invert photolysis matrix

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INVERT_MOD'

CONTAINS
SUBROUTINE invert(a,ainv,acopy,indx,np)
USE ludcmp_mod, ONLY: ludcmp
USE lubksb_mod, ONLY: lubksb
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN)    :: np
INTEGER, INTENT(IN OUT) :: indx(np)
REAL, INTENT(IN)       :: a(np,np)
REAL, INTENT(OUT)      :: ainv(np,np)
REAL, INTENT(OUT)      :: acopy(np,np)

! Local variables
REAL :: d
INTEGER :: i
INTEGER :: j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INVERT'


!     Set up the identity matrix.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO i = 1 , np
  DO j = 1 , np
    ainv(i,j) = 0.0
  END DO
  ainv(i,i) = 1.0
END DO

!     Take copy
DO i = 1 , np
  DO j = 1 , np
    acopy(i,j) = a(i,j)
  END DO
END DO

!     Decompose the matrix once.
CALL ludcmp(acopy,np,np,indx,d)

!     Find the inverse by columns.
DO j = 1 , np
  CALL lubksb(acopy,np,np,indx,ainv(1,j))
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE invert
END MODULE invert_mod
