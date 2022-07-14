! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsso3_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate T-dependent SO3 cross sections, following JPL (2002)

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSSO3_MOD'

CONTAINS
SUBROUTINE acsso3(aso3)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface

INTEGER :: i

! Absorption cross-section of SO3 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(OUT) :: aso3(jpwav)

! Local variables
!     SO3: Burkholder & McKeen. GRL. 24. 3201-3204. 1997
!     196-330 nm @ 298K
REAL, PARAMETER :: aso3d(jpwav) = [                                            &
  (0.0,i=1,57),                                                                &
  78.10, 75.00, 71.15, 66.57, 61.36, 54.77, 48.26,                             &
  41.28, 34.21, 27.80, 22.04, 17.35, 13.50, 10.60, 8.46,                       &
  6.90,  5.66,  4.65,  3.83,  3.13,  2.52,  2.02,  1.61,                       &
  1.27,  0.99,  0.79,  0.61,  0.47,  0.36,  0.27,  0.20,                       &
  0.15,  0.11,  0.08,  0.06,  0.04,  0.03,  0.02,  0.01,                       &
  0.01,  0.01,                                                                 &
  (0.0,i=1,105) ]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSSO3'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Linear interpolation between the cross sections given.
aso3 = 1.0e-20 * aso3d

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE acsso3
END MODULE acsso3_mod
