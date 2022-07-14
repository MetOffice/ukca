! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsno2_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate NO2 temperature dependent cross sections,

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSNO2_MOD'

CONTAINS
SUBROUTINE acsno2(temp, ano2)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface

INTEGER :: i

REAL, INTENT(IN)    :: temp
! Contains the cross sections at model
! wavelengths.Contains the result on exit
REAL, INTENT(IN OUT) :: ano2(jpwav)

! Local variables
! Contains the cross sections at 298K
REAL,PARAMETER ::  ano2t (jpwav) = [                                           &
  (0.0,i=1,48),                                                                &
  0.000e+00,0.000e+00,0.000e+00,2.670e-19,2.780e-19,2.900e-19,                 &
  2.790e-19,2.600e-19,2.420e-19,2.450e-19,2.480e-19,2.750e-19,                 &
  4.145e-19,4.478e-19,4.454e-19,4.641e-19,4.866e-19,4.818e-19,                 &
  5.022e-19,4.441e-19,4.713e-19,3.772e-19,3.929e-19,2.740e-19,                 &
  2.778e-19,1.689e-19,1.618e-19,8.812e-20,7.472e-20,3.909e-20,                 &
  2.753e-20,2.007e-20,1.973e-20,2.111e-20,2.357e-20,2.698e-20,                 &
  3.247e-20,3.785e-20,5.030e-20,5.880e-20,7.000e-20,8.150e-20,                 &
  9.720e-20,1.154e-19,1.344e-19,1.589e-19,1.867e-19,2.153e-19,                 &
  2.477e-19,2.807e-19,3.133e-19,3.425e-19,3.798e-19,4.065e-19,                 &
  4.313e-19,4.717e-19,4.833e-19,5.166e-19,5.315e-19,5.508e-19,                 &
  5.644e-19,5.757e-19,5.927e-19,5.845e-19,6.021e-19,5.781e-19,                 &
  5.999e-19,5.651e-19,5.812e-19,0.000e+00,0.000e+00,0.000e+00,                 &
  (0.0,i=1,83) ]
!     slope A
REAL,PARAMETER ::  a     (jpwav) = [                                           &
  (0.0,i=1,86),0.075, 0.082,-0.053,-0.043,-0.031,-0.162,-0.284,                &
  -0.357,-0.536,-0.686,-0.786,-1.105,-1.355,-1.277,-1.612,                     &
  -1.890,-1.219,-1.921,-1.095,-1.322,-1.102,-0.806,-0.867,                     &
  -0.945,-0.923,-0.738,-0.599,-0.545,-1.129, 0.001,-1.208,                     &
  (0.0,i=1,86) ]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSNO2'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ano2(87:117) = ano2t(87:117) + 1.0e-22*a(87:117)*(temp-273.0)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsno2
END MODULE acsno2_mod
