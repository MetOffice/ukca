! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE quanto1d_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     A subroutine which calculates the quantum yield of O(1D) from O3
!     photolysis at wavelengths less than 310 nm.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='QUANTO1D_MOD'

CONTAINS
SUBROUTINE quanto1d(t,jpwav,wavenm,qeo1d)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: t                ! Temperature in Kelvin
REAL, INTENT(IN)    :: wavenm(jpwav)    ! Wavelength of each interval in nm
! The quantum yield of O(1D) from O3 photolysis at wavelengths less than
! 310 nm.
REAL, INTENT(IN OUT) :: qeo1d(jpwav)

! Local variables
INTEGER :: jw
REAL :: tc
REAL :: a
REAL :: arg
REAL :: b
REAL :: c
REAL :: tau
REAL :: tau2
REAL :: tau3

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='QUANTO1D'



!     Set values of PHI at values not covered by the parameterisation.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
qeo1d(1:91) = 0.9

qeo1d(98:jpwav) = 0.0

!     Intervals 92 to 97.
!     Set up parameter TAU used in the parameterisation.
!     [A Temperature parameter.]
tc = t
tau = tc - 230.0
tau2 = tau*tau
tau3 = tau2*tau

!     Set up the quantum yield of O(1D) from O3 photolysis at
!     wavelengths less than 310 nm using the JPL (evaluation 8)
!     parameterisation due to Moorgat & Kudzus (1978).

DO jw = 92 , 97

  a =  0.332 + 2.5650e-4*tau + 1.152e-5*tau2 + 2.3130e-8*tau3
  b = -0.575 + 5.5900e-3*tau - 1.439e-5*tau2 - 3.2700e-8*tau3
  c =  0.466 + 8.8830e-4*tau - 3.546e-5*tau2 + 3.5190e-7*tau3
  arg= 308.2 + 4.4871e-2*tau + 6.938e-5*tau2 - 2.5452e-6*tau3

  qeo1d(jw) = a*ATAN(b*(wavenm(jw)-arg)) + c
  qeo1d(jw) = MIN(qeo1d(jw),0.9)
  qeo1d(jw) = MAX(qeo1d(jw),0.0)

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE quanto1d
END MODULE quanto1d_mod
