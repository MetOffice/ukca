! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsf12_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSF12_MOD'

CONTAINS
SUBROUTINE acsf12(t,wavenm,af12)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN) :: t                ! Temperature in kelvin
REAL, INTENT(IN) :: wavenm(jpwav)    ! Wavelength of each interval in nm
! Absorption cross-section of F12 in cm2 for temperature T for each
! interval wavelength.
REAL, INTENT(IN OUT) :: af12(jpwav)

! Local variables

INTEGER :: i
INTEGER :: jw
REAL :: tc
REAL :: arg
! Absorption cross-section of F12 at 298K
! CF2Cl2: JPL 1992 T=298K values.
REAL,PARAMETER :: af12t(jpwav) = [                                             &
  (0.0,i=1,36),                                                                &
  0.000e+00,0.000e+00,0.000e+00,0.000e+00,9.716e-20,8.639e-19,                 &
  1.370e-18,1.630e-18,1.752e-18,1.846e-18,1.894e-18,1.778e-18,                 &
  1.655e-18,1.515e-18,1.312e-18,1.030e-18,8.615e-19,6.682e-19,                 &
  4.953e-19,3.567e-19,2.494e-19,1.659e-19,1.088e-19,7.081e-20,                 &
  4.327e-20,2.667e-20,1.753e-20,9.740e-21,5.336e-21,2.976e-21,                 &
  2.572e-21,3.840e-21,5.644e-22,3.270e-22,1.769e-22,8.850e-23,                 &
  4.328e-23,2.236e-23,1.040e-23,3.751e-24,1.146e-24,0.000e+00,                 &
  (0.0,i=1,125)]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSF12'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = MIN (t, 300.0)
tc = MAX (t, 210.0)

!     Wavelength loop.
DO jw = 45 , 71
  arg = 4.1e-4*(wavenm(jw)-184.9)*(tc-298.0)
  af12(jw) = af12t(jw)*EXP(arg)
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsf12
END MODULE acsf12_mod
