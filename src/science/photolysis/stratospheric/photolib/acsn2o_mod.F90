! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsn2o_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!      Calculate the temperature dependent N2O absorption cross
!      section between 173 nm and 240 nm.

! Method:
!      This parameterisation
!      applies for temperatures between 194 and 320 K.  It is taken
!      from the work of Selwyn et al. (1977) presented in JPL EVAL 8
!      (1987) pp. 108. .

!      When the temperature is out of the parameterisation range;
!      < 194 K, the cross section for 194 K is used,
!      > 320 K, the cross section for 320 K is used.


!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSN2O_MOD'

CONTAINS
SUBROUTINE acsn2o(t,jpwav,wavenm,an2o)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: t              ! Temperature in kelvin
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section of N2O in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(IN OUT) :: an2o(jpwav)

! Local variables
!     Polynomial coefficients used.
REAL, PARAMETER :: a1n2o=68.210230
REAL, PARAMETER :: a2n2o=-4.071805
REAL, PARAMETER :: a3n2o= 4.301146e-2
REAL, PARAMETER :: a4n2o=-1.777846e-4
REAL, PARAMETER :: a5n2o= 2.520672e-7

REAL, PARAMETER :: b1n2o=123.401400
REAL, PARAMETER :: b2n2o=-2.116255
REAL, PARAMETER :: b3n2o= 1.111572e-2
REAL, PARAMETER :: b4n2o=-1.881058e-5

INTEGER :: jw

!     Temperature in kelvin.
REAL :: tc
REAL :: tm300
REAL :: lam1
REAL :: lam2
REAL :: lam3
REAL :: lam4
REAL :: arg
REAL :: arga
REAL :: argb

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSN2O'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check is not out of parameterisation range.
tc = t
IF ( tc < 194.0 ) tc = 194.0
IF ( tc > 320.0 ) tc = 320.0

!     Temperature - 300.
tm300 = tc - 300.0

!     Wavelengths between 173 nm & 240 nm.
DO jw = 44 , 76

  !        Various powers of wavelength in nm.
  lam1 = wavenm(jw)
  lam2 = lam1*lam1
  lam3 = lam2*lam1
  lam4 = lam3*lam1

  !        Evaluate the two polynomial expressions.
  arga = a1n2o + a2n2o*lam1 + a3n2o*lam2 + a4n2o*lam3 +                        &
    a5n2o*lam4
  argb = b1n2o + b2n2o*lam1 + b3n2o*lam2 + b4n2o*lam3

  arg = arga + tm300*EXP(argb)

  an2o(jw) = EXP(arg)

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsn2o
END MODULE acsn2o_mod
