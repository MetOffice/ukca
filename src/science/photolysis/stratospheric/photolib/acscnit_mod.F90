! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acscnit_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate ClONO2 temperature dependent cross sections,

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSCNIT_MOD'

CONTAINS
SUBROUTINE acscnit(temp,wavenm,acnita, acnitb)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: temp
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Contains the wavelength intervals
! ACNITA, ACNITB Contain the cross sections times quantum yields
! at model wavelengths.Contains the result on exit.
REAL, INTENT(IN OUT) :: acnita(jpwav)
REAL, INTENT(OUT) :: acnitb(jpwav)

! Local variables

INTEGER :: i

REAL :: tc
LOGICAL, SAVE :: first = .TRUE.
!     ACNIT at 296K
REAL,PARAMETER :: acnitt(jpwav) = [(0.0,i=1,54),                               &
  4.320e-19,1.980e-18,2.911e-18,3.015e-18,2.871e-18,2.785e-18,                 &
  2.780e-18,2.839e-18,2.956e-18,3.097e-18,3.264e-18,3.386e-18,                 &
  3.448e-18,3.392e-18,3.236e-18,2.971e-18,2.640e-18,2.268e-18,                 &
  1.922e-18,1.591e-18,1.314e-18,1.086e-18,8.967e-19,7.444e-19,                 &
  6.074e-19,5.129e-19,4.352e-19,3.703e-19,3.152e-19,2.662e-19,                 &
  2.213e-19,1.840e-19,1.498e-19,1.211e-19,9.519e-20,7.333e-20,                 &
  5.500e-20,4.007e-20,2.969e-20,2.190e-20,1.600e-20,1.142e-20,                 &
  8.310e-21,6.114e-21,4.660e-21,3.657e-21,3.020e-21,2.576e-21,                 &
  2.290e-21,2.079e-21,2.000e-21,1.795e-21,1.590e-21,1.414e-21,                 &
  1.210e-21,1.056e-21,9.090e-22,7.588e-22,6.380e-22,5.376e-22,                 &
  4.440e-22,3.672e-22,3.160e-22,2.314e-22,1.890e-22,5.264e-23,                 &
  (0.0,i=1,83)]
!     Coeffs A1, A2 from Burkholder et al GRL 1994
REAL,PARAMETER :: a1(jpwav) = [(0.0,i=1,54),                                   &
  1.73e-05, 7.50e-05, 1.00e-04, 8.82e-05, 3.61e-05,-5.88e-05,                  &
  -1.95e-04,-3.44e-04,-5.11e-04,-6.59e-04,-7.85e-04,-8.71e-04,                 &
  -9.03e-04,-8.73e-04,-7.83e-04,-6.38e-04,-4.53e-04,-2.35e-04,                 &
  -6.97e-06, 2.19e-04, 4.16e-04, 5.64e-04, 6.78e-04, 7.81e-04,                 &
  9.08e-04, 1.08e-03, 1.26e-03, 1.44e-03, 1.59e-03, 1.74e-03,                  &
  1.88e-03, 2.03e-03, 2.21e-03, 2.37e-03, 2.55e-03, 2.74e-03,                  &
  2.95e-03, 3.28e-03, 3.72e-03, 4.12e-03, 4.53e-03, 4.98e-03,                  &
  5.40e-03, 5.75e-03, 5.92e-03, 5.85e-03, 5.51e-03, 4.92e-03,                  &
  4.02e-03, 3.27e-03, 2.70e-03, 2.08e-03, 1.33e-03, 7.65e-04,                  &
  3.53e-04, 2.39e-04, 4.10e-04, 7.77e-04, 1.38e-03, 2.15e-03,                  &
  3.38e-03, 4.88e-03, 6.70e-03, 8.10e-03, 9.72e-03, 9.96e-03,                  &
  (0.0,i=1,83)]
REAL,PARAMETER :: a2(jpwav) = [(0.0,i=1,54),                                   &
  -1.16e-06,-5.32e-06,-7.85e-06,-8.21e-06,-7.81e-06,-7.52e-06,                 &
  -7.46e-06,-7.62e-06,-7.93e-06,-8.30e-06,-8.69e-06,-9.03e-06,                 &
  -9.25e-06,-9.37e-06,-9.37e-06,-9.27e-06,-9.06e-06,-8.65e-06,                 &
  -7.98e-06,-7.13e-06,-6.36e-06,-6.00e-06,-6.09e-06,-6.45e-06,                 &
  -6.57e-06,-6.03e-06,-5.08e-06,-4.20e-06,-3.50e-06,-2.74e-06,                 &
  -1.87e-06,-9.85e-07, 6.15e-08, 1.13e-06, 2.14e-06, 3.05e-06,                 &
  3.74e-06, 5.29e-06, 7.78e-06, 9.88e-06, 1.20e-05, 1.49e-05,                  &
  1.84e-05, 2.27e-05, 2.70e-05, 3.01e-05, 3.11e-05, 2.86e-05,                  &
  2.07e-05, 1.38e-05, 8.59e-06, 2.01e-06,-7.40e-06,-1.44e-05,                  &
  -1.91e-05,-2.11e-05,-2.05e-05,-1.87e-05,-1.42e-05,-7.14e-06,                 &
  4.47e-06, 1.93e-05, 3.87e-05, 5.57e-05, 7.52e-05, 7.81e-05,                  &
  (0.0,i=1,83)]
REAL, SAVE :: quanty(jpwav) ! Quantum yield for Cl + NO3 channel

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSCNIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (first) THEN
  quanty = MIN(MAX(7.143e-3 * wavenm - 1.6, 0.6), 1.0)
  first = .FALSE.
END IF

tc=temp
tc=MAX(tc, 220.0)
tc=MIN(tc, 298.0)

! Initialise acnitb here as acnita already initialised in fill_spectra
acnitb(:) = 0.0

! ClONO2 cross section, following JPL (2002)
acnita(55:120) = acnitt(55:120)*(1.0 + a1(55:120)*(tc-296.0)                   &
  + a2(55:120)*((tc-296.0)**2))
! ClONO2 cross section, times quantum yield for ClO + NO2 channel
acnitb(55:120) = (1.0 - quanty(55:120)) * acnita(55:120)
! ClONO2 cross section, times quantum yield for Cl  + NO3 channel
acnita(55:120) = quanty(55:120) * acnita(55:120)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acscnit
END MODULE acscnit_mod
