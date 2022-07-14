! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acscos_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate T-dependent COS cross sections, following JPL (2002)
!     The expression is valid for the wavelength range; 186.1-296.3 nm
!                            and the temperature range; 225-295 K.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSCOS_MOD'

CONTAINS
SUBROUTINE acscos(t, aocs)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)  :: t               ! Temperature in kelvin
! Absorption cross-section of COS in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(OUT) :: aocs(jpwav)

! Local variables

INTEGER :: i

REAL :: arg
REAL :: tc
!     COS: JPL 2002 T=295K values.
REAL, PARAMETER :: acos295(jpwav) = [                                          &
  (0.0,i=1,51),                                                                &
  18.9, 8.33, 3.75, 2.21, 1.79, 1.94, 2.48, 3.30, 4.48, 6.12,                  &
  8.19, 10.8, 14.1, 17.6, 21.8, 25.5, 28.2, 30.5, 31.9, 30.2,                  &
  26.8, 22.1, 17.1, 12.5, 8.54, 5.61, 3.51, 2.11, 1.21, .674,                  &
  .361, .193,.0941,.0486,.0248,.0119,.0584,.0264,.0012,.0005,                  &
  .0002,                                                                       &
  (0.0,i=1,111)]
!     COS: JPL 2002 T=225K values.
REAL, PARAMETER :: acos225(jpwav) = [                                          &
  (0.0,i=1,51),                                                                &
  13.0, 5.63, 2.50, 1.61, 1.53, 1.84, 2.44, 3.30, 4.50, 6.17,                  &
  8.27, 10.9, 14.2, 17.6, 21.8, 25.3, 27.7, 29.4, 29.5, 27.4,                  &
  23.7, 18.8, 14.0, 9.72, 6.24, 3.89, 2.29, 1.29, .679, .353,                  &
  .178,.0900,.0419,.0199,.0101,.0048,.0021,.0009,.0005,.0002,                  &
  (0.0,i=1,112)]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSCOS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = MAX(MIN (t, 295.0), 225.0)

arg = (tc - 225.0)/70.0
! Linear interpolation between the cross sections given.
aocs = 1.0e-20 *                                                               &
  (acos295 * arg + acos225 * (1.0 - arg))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE acscos
END MODULE acscos_mod
