! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acssr_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     O2 cross section in Schumann-Runge region.

! Method:
!     Origin of data: Table
!     7.7 in WMO (1985) gives transmission in the Schumann-Runge window
!     (175-205 nm). d ln T/ d (column oxygen) is displayed here, which
!     is the O2 cross section.


!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSSR_MOD'

CONTAINS
SUBROUTINE acssr(tc,jpwav,ao2sr)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
! Total slant path O2 column above the point in molecules per cm^2.
REAL, INTENT(IN)    :: tc
INTEGER, INTENT(IN) :: jpwav !     Number of wavelength intervals
! Absorption cross-section of O2 in cm^2 for each wavelength interval
! in the Schumann-Runge bands.
REAL, INTENT(IN OUT) :: ao2sr(jpwav)

! Local parameters and variables
!     Number of Schumann-Runge wavelength intervals.
INTEGER, PARAMETER :: jpwavesr=17
! Number of O2 column intervals
INTEGER, PARAMETER :: no2col = 20

!     Oxygen column intervol
INTEGER :: jo2
REAL :: frac
LOGICAL, SAVE :: first = .TRUE.

! O2 cross section in Schumann-Runge bands
REAL, SAVE :: sr(jpwavesr, no2col)
REAL, SAVE :: logsr(jpwavesr, no2col)
REAL, SAVE :: o2col(no2col)
REAL, SAVE :: logo2col(no2col)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSSR'


! Main block: Initialize data upon first entry
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (first) THEN
  sr(:,1) = [                                                                  &
    2.0596e-19,1.3273e-19,6.6113e-20,5.1048e-20,4.5025e-20,1.4983e-20,         &
    8.9858e-21,2.9953e-21,2.9932e-21,1.0000e-40,1.0000e-40,1.0000e-40,         &
    1.0000e-40,1.0000e-40,1.0000e-40,1.0000e-40,1.0000e-40 ]
  sr(:,2) = [                                                                  &
    1.9882e-19,1.1754e-19,5.8335e-20,5.0977e-20,4.0015e-20,1.4504e-20,         &
    6.0382e-21,3.6212e-21,1.2071e-21,1.0000e-40,1.0000e-40,1.0000e-40,         &
    1.0000e-40,1.0000e-40,1.0000e-40,1.0000e-40,1.0000e-40 ]
  sr(:,3) = [                                                                  &
    1.8704e-19,9.5822e-20,5.1962e-20,4.8400e-20,3.6629e-20,1.2836e-20,         &
    6.4058e-21,2.9866e-21,1.2793e-21,8.5257e-22,4.2614e-22,1.0000e-40,         &
    1.0000e-40,1.0000e-40,1.0000e-40,1.0000e-40,1.0000e-40 ]
  sr(:,4) = [                                                                  &
    1.6133e-19,6.4581e-20,4.2558e-20,4.2452e-20,2.9852e-20,1.1903e-20,         &
    6.1489e-21,2.8427e-21,1.1954e-21,4.4801e-22,5.9735e-22,1.4923e-22,         &
    1.0000e-40,1.0000e-40,1.0000e-40,1.0000e-40,1.0000e-40 ]
  sr(:,5) = [                                                                  &
    1.2039e-19,3.7382e-20,2.9637e-20,3.1606e-20,2.0773e-20,1.0706e-20,         &
    5.7048e-21,2.7256e-21,1.1954e-21,5.4263e-22,5.4263e-22,1.0845e-22,         &
    1.0000e-40,1.0000e-40,1.0000e-40,1.0000e-40,1.0000e-40 ]
  sr(:,6) = [                                                                  &
    8.1089e-20,2.2320e-20,1.7193e-20,2.0043e-20,1.2370e-20,8.3474e-21,         &
    4.9413e-21,2.4020e-21,1.1291e-21,4.5832e-22,5.2091e-22,1.4555e-22,         &
    4.1557e-23,1.0000e-40,1.0000e-40,1.0000e-40,1.0000e-40 ]
  sr(:,7) = [                                                                  &
    5.2569e-20,1.4185e-20,9.5606e-21,1.1560e-20,7.1082e-21,5.4135e-21,         &
    3.7658e-21,1.8750e-21,1.0425e-21,4.5744e-22,5.1708e-22,1.2648e-22,         &
    3.3682e-23,8.4137e-24,1.0000e-40,1.0000e-40,1.0000e-40 ]
  sr(:,8) = [                                                                  &
    3.3180e-20,9.4479e-21,5.6342e-21,6.7482e-21,4.2844e-21,3.1560e-21,         &
    2.4997e-21,1.3218e-21,8.5520e-22,4.1620e-22,4.9829e-22,1.3120e-22,         &
    3.6328e-23,1.0000e-40,3.6300e-24,3.6300e-24,3.6274e-24  ]
  sr(:,9) = [                                                                  &
    2.1783e-20,6.5782e-21,3.6155e-21,4.2548e-21,2.7253e-21,1.9454e-21,         &
    1.5785e-21,9.1690e-22,6.2044e-22,3.6428e-22,4.5727e-22,1.2660e-22,         &
    3.4726e-23,3.2997e-24,1.6499e-24,1.0000e-40,1.0000e-40  ]
  sr(:,10) = [                                                                 &
    1.6102e-20,4.7863e-21,2.5785e-21,2.8755e-21,1.8324e-21,1.2969e-21,         &
    1.0247e-21,6.3919e-22,4.0578e-22,2.7130e-22,3.8879e-22,1.1894e-22,         &
    3.2220e-23,4.6906e-24,1.5631e-24,1.0000e-40,1.0000e-40  ]
  sr(:,11) = [                                                                 &
    1.4089e-20,3.6881e-21,2.0761e-21,3.7937e-21,1.3054e-21,9.1906e-22,         &
    7.2141e-22,4.5887e-22,2.7029e-22,2.1300e-22,3.0038e-22,1.0554e-22,         &
    2.9401e-23,4.2685e-24,3.8789e-24,3.8765e-25,7.7525e-25  ]
  sr(:,12) = [                                                                 &
    1.4608e-20,3.1322e-21,1.8465e-21,8.0188e-22,8.9807e-22,6.8966e-22,         &
    5.5373e-22,3.4768e-22,1.9326e-22,1.6128e-22,2.1394e-22,8.7938e-23,         &
    2.4211e-23,5.0370e-24,1.4084e-24,6.0272e-25,2.0093e-25  ]
  sr(:,13) = [                                                                 &
    1.4608e-20,2.9988e-21,1.7814e-21,1.3742e-21,7.9193e-22,5.3877e-22,         &
    4.6246e-22,2.8948e-22,1.4895e-22,1.2973e-22,1.5169e-22,6.8160e-23,         &
    2.0371e-23,5.7585e-24,2.4894e-24,6.4781e-25,5.3968e-25  ]
  sr(:,14) = [                                                                 &
    1.4608e-20,3.1296e-21,1.7853e-21,1.1892e-21,6.3662e-22,4.2117e-22,         &
    3.8823e-22,1.8296e-22,1.1938e-22,1.0509e-22,1.0961e-22,5.1871e-23,         &
    1.7125e-23,6.7396e-24,2.6245e-24,7.7167e-25,4.1525e-25  ]
  sr(:,15) = [                                                                 &
    1.4608e-20,2.9777e-21,1.5883e-21,9.6037e-22,5.0665e-22,3.0382e-22,         &
    2.7925e-22,1.3249e-22,8.6811e-23,7.6305e-23,7.3698e-23,3.7378e-23,         &
    1.3751e-23,5.5964e-24,2.6709e-24,7.9697e-25,4.1381e-25  ]
  sr(:,16) = [                                                                 &
    1.4608e-20,2.9777e-21,1.2210e-21,7.0899e-22,3.6453e-22,2.0328e-22,         &
    1.7512e-22,8.4507e-23,5.6116e-23,4.8982e-23,4.3889e-23,2.4642e-23,         &
    9.6018e-24,4.3275e-24,2.4343e-24,6.6242e-25,4.1886e-25  ]
  sr(:,17) = [                                                                 &
    1.4608e-20,2.9777e-21,1.2210e-21,5.5897e-22,2.7411e-22,1.4980e-22,         &
    1.1815e-22,5.8790e-23,3.6523e-23,3.2044e-23,2.5235e-23,1.5842e-23,         &
    6.5220e-24,3.2660e-24,2.2322e-24,6.0985e-25,4.1236e-25  ]
  sr(:,18) = [                                                                 &
    1.4608e-20,2.9777e-21,1.2210e-21,5.5897e-22,2.3540e-22,1.2884e-22,         &
    1.7942e-22,4.6834e-23,2.6335e-23,2.2983e-23,1.5504e-23,1.0513e-23,         &
    4.5950e-24,2.6892e-24,2.1383e-24,5.7770e-25,3.9394e-25  ]
  sr(:,19) = [                                                                 &
    1.4608e-20,2.9777e-21,1.2210e-21,5.5897e-22,2.3540e-22,1.1905e-22,         &
    5.0663e-23,4.0295e-23,2.1386e-23,1.7647e-23,1.0693e-23,7.3588e-24,         &
    3.3738e-24,2.3488e-24,2.0769e-24,5.6322e-25,3.9207e-25  ]
  sr(:,20) = [                                                                 &
    1.4608e-20,2.9777e-21,1.2210e-21,5.5897e-22,2.3540e-22,1.1905e-22,         &
    5.0663e-23,3.5399e-23,1.8521e-23,1.4160e-23,8.0985e-24,5.4163e-24,         &
    2.5977e-24,2.1190e-24,2.0413e-24,5.4415e-25,3.8368e-25  ]

  logsr = LOG(sr)

  o2col = [                                                                    &
    5.3368e+16,1.0627e+17,2.4629e+17,6.4304e+17,1.7548e+18,4.7351e+18,         &
    1.2299e+19,3.0403e+19,7.1301e+19,1.5943e+20,3.4126e+20,6.9993e+20,         &
    1.3797e+21,2.6309e+21,4.9365e+21,9.3681e+21,1.8360e+22,3.7372e+22,         &
    7.8501e+22,1.6851e+23 ]

  logo2col = LOG(o2col)
  first = .FALSE.
END IF

! Perform linear interpolation in log-log space of cross section in
! Schumann-Runge window. In case the O2 column is outside the window
! covered, take

IF (tc < o2col(1)) THEN
  jo2 = 1
  frac = 0.0
ELSE IF (tc  >=  o2col(no2col)) THEN
  jo2 = no2col-1
  frac = 1.0
ELSE
  jo2 = no2col-1
  DO WHILE (o2col(jo2) > tc)
    jo2 = jo2 - 1
  END DO
  frac = (LOG(tc)          - logo2col(jo2))/                                   &
    (logo2col(jo2 + 1) - logo2col(jo2))
END IF

ao2sr(46:45 + jpwavesr) =                                                      &
  EXP((1.0 - frac) * logsr(:,jo2) + frac * logsr(:,jo2 + 1))
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acssr
END MODULE acssr_mod
