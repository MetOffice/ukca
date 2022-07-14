! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsmc_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Subroutine which calculates the methyl chloride (CH3Cl) absorption
!     cross section based on P. C. Simon et al. (1988).

! Method:
!     Journal of atmospheric chemistry, Vol. 7, pp. 107-135, 1988.

!     This is done via a polynomial expression of the form;
!     log10(sigma)=A(lamda) + T B(lamda)

!     The expression is valid for the wavelength range; 174-226 nm
!                            and the temperature range; 210-300 K.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSMC_MOD'

CONTAINS
SUBROUTINE acsmc(t,jpwav,wavenm,amc)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t                ! Temperature in kelvin.
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)    ! Wavelength of each interval in nm
! Absorption cross-section of CH3Cl in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(IN OUT) :: amc(jpwav)

! Local variables

!     Polynomial coefficients.
REAL, PARAMETER :: a0=-299.80
REAL, PARAMETER :: a1= 5.1047
REAL, PARAMETER :: a2=-3.363e-2
REAL, PARAMETER :: a3= 9.5805e-5
REAL, PARAMETER :: a4=-1.0135e-7

REAL, PARAMETER :: b0=-7.1727
REAL, PARAMETER :: b1= 1.4837e-1
REAL, PARAMETER :: b2=-1.1463e-3
REAL, PARAMETER :: b3= 3.9188e-6
REAL, PARAMETER :: b4=-4.9994e-9

INTEGER :: jw

!     Temperature in kelvin.
REAL :: tc
REAL :: lam
REAL :: lam2
REAL :: lam3
REAL :: lam4
REAL :: arg

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSMC'


!     Check that temperature is in range.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
tc = t
IF ( tc > 300.0 ) tc = 300.0
IF ( tc < 210.0 ) tc = 210.0

!     Wavelength.
DO jw = 45 , 67

  !        Wavelength in nm.
  lam = wavenm(jw)

  lam2 = lam*lam
  lam3 = lam*lam2
  lam4 = lam*lam3
  arg = a0 + a1*lam + a2*lam2 + a3*lam3 + a4*lam4 +                            &
    tc*(b0+b1*lam+b2*lam2+b3*lam3+b4*lam4)
  amc(jw) = 10.0**arg

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsmc
END MODULE acsmc_mod
