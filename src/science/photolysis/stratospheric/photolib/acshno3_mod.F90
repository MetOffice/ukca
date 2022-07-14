! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acshno3_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate HNO3 temperature dependent cross sections

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSHNO3_MOD'

CONTAINS
SUBROUTINE acshno3(temp,ahno3)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: temp
! Contains the cross sections at model wavelengths.Contains the result on exit.
REAL, INTENT(IN OUT) :: ahno3(jpwav)

! Local variables

INTEGER :: i

REAL :: tc            ! Contains the wavelength intervals.
! Contains the cross sections at 298K
REAL,PARAMETER :: ahno3t(jpwav)  = [(0.0,i=1,48),                              &
  0.000e+00,8.305e-18,1.336e-17,1.575e-17,1.491e-17,1.385e-17,                 &
  1.265e-17,1.150e-17,1.012e-17,8.565e-18,6.739e-18,5.147e-18,                 &
  3.788e-18,2.719e-18,1.796e-18,1.180e-18,7.377e-19,4.487e-19,                 &
  2.810e-19,1.826e-19,1.324e-19,1.010e-19,8.020e-20,6.479e-20,                 &
  5.204e-20,4.178e-20,3.200e-20,2.657e-20,2.298e-20,2.086e-20,                 &
  1.991e-20,1.962e-20,1.952e-20,1.929e-20,1.882e-20,1.804e-20,                 &
  1.681e-20,1.526e-20,1.335e-20,1.136e-20,9.242e-21,7.186e-21,                 &
  5.320e-21,3.705e-21,2.393e-21,1.442e-21,8.140e-22,4.131e-22,                 &
  1.970e-22,9.434e-23,4.310e-23,2.204e-23,1.030e-23,5.841e-24,                 &
  4.170e-24,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,                 &
  (0.0,i=1,95) ]
! Contains B coefficient. Intercept B from Burkholder
REAL,PARAMETER :: b     (jpwav) = [ (0.0,i=1,54),                              &
  0.000e+00,0.000e+00,1.152e-03,1.668e-03,1.653e-03,1.673e-03,                 &
  1.720e-03,1.750e-03,1.817e-03,1.935e-03,2.060e-03,2.168e-03,                 &
  2.178e-03,2.195e-03,2.106e-03,1.987e-03,1.840e-03,1.782e-03,                 &
  1.838e-03,1.897e-03,1.970e-03,1.978e-03,1.855e-03,1.655e-03,                 &
  1.416e-03,1.247e-03,1.162e-03,1.121e-03,1.136e-03,1.199e-03,                 &
  1.315e-03,1.493e-03,1.637e-03,1.767e-03,1.928e-03,2.139e-03,                 &
  2.380e-03,2.736e-03,3.139e-03,3.695e-03,4.230e-03,5.151e-03,                 &
  6.450e-03,7.327e-03,9.750e-03,1.013e-02,1.180e-02,1.108e-02,                 &
  9.300e-03,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,                 &
  (0.0,i=1,95) ]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSHNO3'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
tc=temp
tc=MAX(tc, 200.0)
tc=MIN(tc, 360.0)

ahno3(50:103) = ahno3t(50:103)*EXP(b(50:103)*(tc-298.0))
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acshno3
END MODULE acshno3_mod
