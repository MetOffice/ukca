! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsccl4_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     A subroutine which calculates the carbon tetra chloride (CCl4)
!     absorption cross section

! Method:
!     based on P. C. Simon et al. (1988).

!     Journal of atmospheric chemistry, Vol. 7, pp. 107-135, 1988.

!     This is done via a polynomial expression of the form;
!     log10(sigma)=A(lamda) + T B(lamda)

!     The expression is valid for the wavelength range; 194-250 nm
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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSCCL4_MOD'

CONTAINS

SUBROUTINE acsccl4(t,jpwav,wavenm,accl4)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: t              ! Temperature in kelvin.
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section of CCl4 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(IN OUT) :: accl4(jpwav)

! Local variables.
!     Polynomial coefficients.
REAL, PARAMETER :: a0=-37.104
REAL, PARAMETER :: a1=-5.8218e-1
REAL, PARAMETER :: a2= 9.9974e-3
REAL, PARAMETER :: a3=-4.6765e-5
REAL, PARAMETER :: a4= 6.8501e-8

REAL, PARAMETER :: b0= 1.0739
REAL, PARAMETER :: b1=-1.6275e-2
REAL, PARAMETER :: b2= 8.8141e-5
REAL, PARAMETER :: b3=-1.9811e-7
REAL, PARAMETER :: b4= 1.5022e-10

INTEGER :: jw
REAL :: lam
REAL :: lam2
REAL :: lam3
REAL :: lam4
REAL :: arg
REAL :: tc

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSCCL4'


!     Check that temperature is in range.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
tc = t
IF ( tc > 300.0 ) tc = 300.0
IF ( tc < 210.0 ) tc = 210.0

!     Wavelength loop.
DO jw = 57 , 79

  !        Wavelength in nm.
  lam = wavenm(jw)

  lam2 = lam*lam
  lam3 = lam*lam2
  lam4 = lam*lam3

  arg = a0 + a1*lam + a2*lam2 + a3*lam3 + a4*lam4 +                            &
    tc*(b0+b1*lam+b2*lam2+b3*lam3+b4*lam4)

  accl4(jw) = 10.0**arg

  !     End of the wavelength loop.
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsccl4
END MODULE acsccl4_mod
