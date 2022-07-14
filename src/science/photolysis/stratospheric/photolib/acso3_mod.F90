! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acso3_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     A subroutine which calculates the ozone absorption cross section
!     based on John E. Frederick (1985).

! Method:
!     The temperature dependent cross
!     section data set is that of A. M. Bass of the National Bureau of
!     Standards provided by R. D. McPeters.

!     The temperature dependence is in the 3rd significant figure between
!     263.158 - 266.167 nm. (Intervals 84-102).

!     The temperature range covered is 203 to 298 K.

!     The fit used is a quadratic fit of the form;
!     sigma(O3,t)={C0(i)+C1(i)(T-230)+C2(i)(T-230)^2}10^-n

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSO3_MOD'

CONTAINS
SUBROUTINE acso3(t,jpwav,ao3)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t           ! Temperature in kelvin
INTEGER, INTENT(IN) :: jpwav       ! Wavelength of each interval in nm
! Absorption cross-section of O3 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(IN OUT) :: ao3(jpwav)

! Local varaiables

INTEGER :: i

REAL :: tc                         !     Temperature in kelvin.

!     Local variables & Polynomial coefficients.
REAL :: tm230
REAL :: tm2302

REAL,PARAMETER :: c(3,19) = RESHAPE([                                          &
  9.6312e+0 , 1.1875e-3 , -1.7386e-5 , 8.3211e+0 ,                             &
  3.6495e-4 , 2.4691e-6 , 6.8810e+0 , 2.4598e-4 , 1.1692e-5 ,                  &
  5.3744e+0 , 1.0325e-3 , 1.2573e-6 , 3.9575e+0 , 1.6851e-3 ,                  &
 -6.8648e-6 , 2.7095e+0 , 1.4502e-3 ,-2.8925e-6 , 1.7464e+0 ,                  &
  8.9350e-4 , 3.5914e-6 , 1.0574e+0 , 7.8270e-4 , 2.0024e-6 ,                  &
  5.9574e+0 , 4.9448e-3 , 3.6589e-5 , 3.2348e+0 , 3.5392e-3 ,                  &
  2.4769e-5 , 1.7164e+0 , 2.4542e-3 , 1.6913e-5 , 8.9612e+0 ,                  &
  1.4121e-2 , 1.2498e-4 , 4.5004e+0 , 8.4327e-3 , 7.8903e-5 ,                  &
  2.1866e+0 , 4.8343e-3 , 5.1970e-5 , 1.0071e+1 , 3.3409e-2 ,                  &
  2.6621e-4 , 5.0848e+0 , 1.8178e-2 , 1.6301e-4 , 2.1233e+0 ,                  &
  8.8453e-3 , 1.2633e-4 , 8.2861e+0 , 4.2692e-2 , 8.7057e-4 ,                  &
  2.9415e+0 , 5.3051e-2 , 3.4964e-4],[3,19])
REAL,PARAMETER :: n(19) = [ (18.0,i=1,8), (19.0,i=1,3), (20.0,i=1,3),          &
 (21.0,i=1,3), (22.0,i=1,2) ]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSO3'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = t
IF ( tc > 298.0 ) tc = 298.0
IF ( tc < 203.0 ) tc = 203.0

tm230 = tc - 230.0
tm2302 = tm230*tm230

ao3(84:102) = (c(1,1:19)+c(2,1:19)*tm230+c(3,1:19)*tm2302)                     &
  *(10.0**(-n(1:19)))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acso3
END MODULE acso3_mod
