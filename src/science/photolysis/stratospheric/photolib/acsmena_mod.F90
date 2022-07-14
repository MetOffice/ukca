! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsmena_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate T-dependent MeONO2 cross sections
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


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSMENA_MOD'

CONTAINS
SUBROUTINE acsmena(t, amena)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

REAL, INTENT(IN)    :: t              ! Temperature in kelvin
! Absorption cross-section of COS in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(IN OUT) :: amena(jpwav)
! Local variables

INTEGER :: i

REAL :: tc
!     JPL 2002 T=295K values.
REAL, PARAMETER :: amena298(jpwav) = [                                         &
  (0.0,i=1,75),                                                                &
  6.049e-20  , 5.0688e-20, 4.142e-20  , 3.77e-20  ,                            &
  3.4972e-20 , 3.3116e-20, 3.1512e-20 , 2.9788e-20,                            &
  2.7758e-20 , 2.504e-20 , 2.2262e-20 , 1.9244e-20,                            &
  1.6052e-20 , 1.2914e-20, 9.99601e-21, 7.372e-21 ,                            &
  5.13921e-21, 3.3664e-21, 2.092e-21  , 1.34e-21  ,                            &
  6.33e-22   , 3.16e-22  , 1.44e-22   , 6.61e-23  ,                            &
  2.74e-23   , 1.22e-23  , (0.0,i=1,102) ]
!     JPL 2002 T=225K values.
REAL, PARAMETER :: b(jpwav) = [                                                &
  (0.0,i=1,75),                                                                &
  0.003499 , 0.0033888, 0.0032636, 0.003059 ,                                  &
  0.0029152, 0.0028256, 0.0028262, 0.0028552,                                  &
  0.0029182, 0.003032 , 0.003164 , 0.0033214,                                  &
  0.0034962, 0.0037098, 0.0039256, 0.004212 ,                                  &
  0.0045922, 0.0050392, 0.0056062, 0.00633  ,                                  &
  0.00734  , 0.00874  , 0.00997  , 0.0136   ,                                  &
  0.0136   , 0.0136   , (0.0,i=1,102) ]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSMENA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!     Check that temperature is in range.
tc = MAX(MIN (t, 330.0), 240.0)

! ln sigma = ln sigma(298K) + b * (T - 298K)
amena(76:101) = EXP(LOG(amena298(76:101)) +                                    &
  b(76:101) * (tc - 298.0))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE acsmena
END MODULE acsmena_mod
