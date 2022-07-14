! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE quanto12_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     A subroutine which calculates the quantum yield of O(1D) from O3
!     photolysis based on parameteristn of Michelson et al GRL Oct 1994.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='QUANTO12_MOD'

CONTAINS
SUBROUTINE quanto12(t,jpwav,wavenm,qeo1d)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t                ! Temperature in Kelvin
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)    ! Wavelength of each interval in nm
REAL, INTENT(IN OUT) :: qeo1d(jpwav)     ! Quantum yield of O(1D) from O3 photol

! Local variables
REAL :: phifac
!     Data from Michelsen et al. 1994
REAL,PARAMETER :: phi1da(5) = [1.01,  2.45, 16.55,  13.82,   11.8]
REAL,PARAMETER :: phi1db(5) = [3.93, 294.2, 846.3, 1061.9, 1435.0]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='QUANTO12'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Factor = hc (planck (Js) x speed of light (cms-1))
phifac = 6.626e-34*2.99e10

qeo1d(1:85) = 0.9

qeo1d(86:93) = 1.98 - (301.0/wavenm(86:93))

qeo1d(94:98) = phi1da(1:5)*EXP((-phi1db(1:5)*phifac)                           &
  /(1.38e-23*t))

qeo1d(99:jpwav) = 0.0
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE quanto12
END MODULE quanto12_mod
