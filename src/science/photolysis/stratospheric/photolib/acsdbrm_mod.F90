! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsdbrm_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Description:
!     Calculate T-dependent CH2Br2 cross sections

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSDBRM_MOD'

CONTAINS

SUBROUTINE acsdbrm(t,adbrm)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

REAL, INTENT(IN) :: t               ! Temperature in kelvin.
! Absorption cross-section of CH2BR2 in cm^2 for temperature T for each
! interval wavelength. On first entry, ADBRM needs to have the cross section
! at 298 K.
REAL, INTENT(IN OUT) :: adbrm(jpwav)
! Local variables

LOGICAL, SAVE :: first = .TRUE.
REAL :: tc
REAL,SAVE :: adbrm298(29)

REAL, SAVE :: b(29) = [                                                        &
  -0.0006464,-0.0015756,-0.0019648,-0.0018544,-0.0017320,                      &
  -0.0015870,-0.0013920,-0.0011112,-0.0007392,-0.0002424,                      &
  0.0002150, 0.0006350, 0.0012332, 0.0018720, 0.0024300,                       &
  0.0030872, 0.0038616, 0.0046850, 0.0054876, 0.0063066,                       &
  0.0073240, 0.0081658, 0.0092548, 0.0114284, 0.0132280,                       &
  0.0146400, 0.0157800, 0.0179040, 0.0227640 ]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSDBRM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (first) THEN
  adbrm298 = adbrm(65:93)
  first = .FALSE.
END IF

tc = MIN(MAX(t,250.0),348.0) - 298.0
adbrm(65:93) = adbrm298 * EXP(b*tc)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE acsdbrm
END MODULE acsdbrm_mod
