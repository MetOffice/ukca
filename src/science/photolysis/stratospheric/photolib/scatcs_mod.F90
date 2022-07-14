! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE scatcs_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Subroutine which calculates the Rayleigh scattering cross section
!     based on Nicolet (1984).

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SCATCS_MOD'

CONTAINS
SUBROUTINE scatcs(jpwav,scs,wavenm)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! The Rayleigh scattering cross section for air molecules in
! molecules per cm^2.
REAL, INTENT(IN OUT) :: scs(jpwav)

! Local variables
INTEGER :: jw
REAL :: xlamda
REAL :: chi

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SCATCS'



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
scs(1:45) = 0.0
DO jw = 46 , 143
  xlamda  = wavenm(jw)*1.0e-3
  chi     = 0.389*xlamda + (0.09426/xlamda) - 0.3228
  scs(jw) = 4.02e-28*(xlamda**(-4.0 - chi))
END DO
DO jw = 144 , jpwav
  xlamda  = wavenm(jw)*1.0e-3
  scs(jw) = 4.02e-28*(xlamda**(-4.04))
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE scatcs
END MODULE scatcs_mod
