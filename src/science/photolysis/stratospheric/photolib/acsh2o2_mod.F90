! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsh2o2_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate H2O2 temperature dependent cross sections,

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSH2O2_MOD'

CONTAINS
SUBROUTINE acsh2o2(temp,jpwav,wavenm,ah2o2)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: temp
REAL, INTENT(IN)    :: wavenm(jpwav) ! Contains the wavelength intervals
! Contains the cross sections at model wavelengths.Contains the result on exit
REAL, INTENT(IN OUT) :: ah2o2(jpwav)

! Local variables
REAL, PARAMETER :: a0= 6.4761e4
REAL, PARAMETER :: a1=-9.2170972e2
REAL, PARAMETER :: a2= 4.535649
REAL, PARAMETER :: a3=-4.4589016e-3
REAL, PARAMETER :: a4=-4.035101e-5
REAL, PARAMETER :: a5= 1.6878206e-7
REAL, PARAMETER :: a6=-2.652014e-10
REAL, PARAMETER :: a7= 1.5534675e-13
REAL, PARAMETER :: b0= 6.8123e3
REAL, PARAMETER :: b1=-5.1351e1
REAL, PARAMETER :: b2= 1.1522e-1
REAL, PARAMETER :: b3=-3.0493e-5
REAL, PARAMETER :: b4=-1.0924e-7

REAL :: chi
REAL :: lam
REAL :: lam2
REAL :: lam3
REAL :: lam4
REAL :: lam5
REAL :: lam6
REAL :: lam7
REAL :: tc

INTEGER :: jw  ! Loop counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSH2O2'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = temp
IF ( tc > 400.0 ) tc = 400.0
IF ( tc < 200.0 ) tc = 200.0

chi=1.0/(1.0 + EXP(-1265.0/tc))

!     Wavelength loop 260nm - 350nm
DO jw = 83,103
  lam =wavenm(jw)
  lam2=lam*lam
  lam3=lam*lam2
  lam4=lam*lam3
  lam5=lam*lam4
  lam6=lam*lam5
  lam7=lam*lam6
  ah2o2(jw)= 1.0e-21*(chi*(a0 + a1*lam  + a2*lam2 + a3*lam3                    &
    +                 a4*lam4 + a5*lam5 + a6*lam6 + a7*lam7)                   &
    + (1.0-chi)*(b0 + b1*lam  + b2*lam2 + b3*lam3 + b4*lam4))
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsh2o2
END MODULE acsh2o2_mod
