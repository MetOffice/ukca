! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsh2so4_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate T-dependent H2SO4 cross sections, following JPL (2002)

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSH2SO4_MOD'

CONTAINS
SUBROUTINE acsh2so4(ah2so4)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Absorption cross-section of H2SO4 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(OUT) :: ah2so4(jpwav)

! Local variables

INTEGER :: i

! Taken from Slimane Bekki's model. Based on HCl data from JPL
REAL,PARAMETER :: ah2so4d(jpwav) = [                                           &
  (0.0,i=1,32),                                                                &
  211.00, 281.00, 281.00, 345.00, 345.00, 382.00, 382.00, 332.00,              &
  332.00, 248.00, 248.00, 163.00, 163.00, 109.00, 109.00, 58.80,               &
  58.80,  58.80,  31.30,  31.30,  31.30,  14.50,  14.50,  14.50,               &
  6.18,   6.18,   2.56,   2.46,   2.56,   0.98,   0.98,   0.98,                &
  0.39,   0.39,   0.14,   0.14,   0.05,   0.05,                                &
  (0.0,i=1,133) ]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSH2SO4'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Linear interpolation between the cross sections given.
ah2so4 = 1.0e-20 * ah2so4d *0.016  ! S. Bekki data

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE acsh2so4
END MODULE acsh2so4_mod
