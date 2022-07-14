! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsf11_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate T-dependent F11 cross sections
!     The expression is valid for the wavelength range; 174-230 nm
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


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSF11_MOD'

CONTAINS
SUBROUTINE acsf11(t,wavenm,af11)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN) :: t                ! Temperature in kelvin
REAL, INTENT(IN) :: wavenm(jpwav)    ! Wavelength of each interval in nm
! Absorption cross-section of F11 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(IN OUT) :: af11(jpwav)

! Local variables

INTEGER :: i
INTEGER :: jw

REAL :: tc
REAL :: arg
! Absorption cross-section of F11 at 298K
! CFCl3: JPL 1992 T=298K values.
REAL,PARAMETER :: af11t(jpwav) = [                                             &
  (0.0,i=1,36),                                                                &
  0.000e+00,0.000e+00,0.000e+00,0.000e+00,3.139e-19,2.478e-18,                 &
  3.182e-18,3.168e-18,3.141e-18,3.098e-18,3.042e-18,3.096e-18,                 &
  2.968e-18,2.765e-18,2.558e-18,2.319e-18,2.107e-18,1.839e-18,                 &
  1.574e-18,1.332e-18,1.092e-18,8.911e-19,7.221e-19,5.751e-19,                 &
  4.389e-19,3.340e-19,2.377e-19,1.700e-19,1.171e-19,7.662e-20,                 &
  5.082e-20,3.184e-20,1.970e-20,1.206e-20,8.000e-21,4.834e-21,                 &
  2.831e-21,1.629e-21,9.327e-22,5.209e-22,3.013e-22,1.617e-22,                 &
  9.035e-23,5.427e-23,3.474e-23,2.141e-23,9.102e-24,1.499e-25,                 &
  (0.0,i=1,119)]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSF11'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!     Check that temperature is in range.
tc = MIN (t, 300.0)
tc = MAX (tc,210.0)

!     Wavelength loop.
DO jw = 45 , 72
  arg = 1.0e-4*(wavenm(jw)-184.9)*(tc-298.0)
  af11(jw) = af11t(jw)*EXP(arg)
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsf11
END MODULE acsf11_mod
