! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsn2o5_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!      Calculate the temperature dependent N2O5 cross section

! Method:
!      using the parameterisation of Yao et al. (1982), given in JPL
!      EVAL 8 (1987) pp. 110-111. .

!      This applies for wavelengths between 285 and 380 nm, and
!      temperatures 225 K to 300 K.

!      If the temperature is below 225 K, the value for 225 K is used.
!      If the temperature is above 300 K, the value for 300 K is used.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSN2O5_MOD'

CONTAINS
SUBROUTINE acsn2o5(t,jpwav,wavenm,an2o5)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: t              ! Temperature in kelvin
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section of N2O5 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(IN OUT) :: an2o5(jpwav)

! Local variables
INTEGER :: jw
REAL :: lam

!     Temperature in kelvin.
REAL :: tc
REAL :: arg

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSN2O5'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check is not out of parameterisation range.
tc = t
IF ( tc < 225.0 ) tc = 225.0
IF ( tc > 300.0 ) tc = 300.0

!     Wavelengths between 285 nm & 380 nm.
DO jw = 89 , 109

  !        Wavelength in nm.
  lam = wavenm(jw)

  !        Evaluate the parameterisation expression.
  arg = 2.735 + ((4728.5-17.127*lam)/tc)

  !        Calculate the cross section.
  an2o5(jw) = 1.0e-20*EXP(arg)

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsn2o5
END MODULE acsn2o5_mod
