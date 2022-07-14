! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE ludcmp_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     For matrix inversion

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LUDCMP_MOD'

CONTAINS
SUBROUTINE ludcmp(a,n,np,indx,d)
USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: np
INTEGER, INTENT(IN) :: n
REAL, INTENT(IN OUT) :: d
INTEGER, INTENT(OUT):: indx(n)
REAL, INTENT(IN OUT) :: a(np,np)

! Local variables
INTEGER, PARAMETER :: nmax=100
REAL, PARAMETER :: var_tiny=1.0e-20

INTEGER :: i
INTEGER :: imax
INTEGER :: j
INTEGER :: k
INTEGER :: icode

REAL :: aamax
REAL :: dum
REAL :: var_sum
REAL :: vv(nmax)

CHARACTER(LEN=errormessagelength) :: cmessage           ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LUDCMP'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
d=1.0
DO i=1,n
  aamax=0.0
  DO j=1,n
    IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
  END DO
  IF (aamax == 0.0) THEN
    ! Fatal error
    cmessage='Singular matrix in subroutine ludcmp'
    icode=1
    CALL ereport('LUDCMP',icode,cmessage)
  END IF
  vv(i)=1.0/aamax
END DO
DO j=1,n
  IF (j > 1) THEN
    DO i=1,j-1
      var_sum=a(i,j)
      IF (i > 1) THEN
        DO k=1,i-1
          var_sum=var_sum-a(i,k)*a(k,j)
        END DO
        a(i,j)=var_sum
      END IF
    END DO
  END IF
  aamax=0.0
  DO i=j,n
    var_sum=a(i,j)
    IF (j > 1) THEN
      DO k=1,j-1
        var_sum=var_sum-a(i,k)*a(k,j)
      END DO
      a(i,j)=var_sum
    END IF
    dum=vv(i)*ABS(var_sum)
    IF (dum >= aamax) THEN
      imax=i
      aamax=dum
    END IF
  END DO
  IF (j /= imax) THEN
    DO k=1,n
      dum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum
    END DO
    d=-d
    vv(imax)=vv(j)
  END IF
  indx(j)=imax
  IF (j /= n) THEN
    IF (a(j,j) == 0.0) a(j,j)=var_tiny
    dum=1.0/a(j,j)
    DO i=j+1,n
      a(i,j)=a(i,j)*dum
    END DO
  END IF
END DO
IF (a(n,n) == 0.0) a(n,n)=var_tiny
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ludcmp
END MODULE ludcmp_mod
