! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsno_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate the NO absorption cross section accounting for
!     the NO absorption in the DEL(0-0) & the DEL(0-1) bands.

!  Method:
!     Calculated from the Mark Allen & John E Frederick
!     parameterisation, taken from;

!     The Journal of atmospheric sciences, Vol. 39, September 82.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSNO_MOD'

CONTAINS
SUBROUTINE acsno(angle,p,vc,ano)
USE photol_constants_mod, ONLY: pi_over_180 => const_pi_over_180
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: angle      ! Solar zenith angle in radians
REAL, INTENT(IN)    :: p          ! Pressure in mb
REAL, INTENT(IN)    :: vc         ! Vertical O2 column above the point
!                                   in molecules per cm^2.
! Absorption cross-section of NO in cm^2 for pressure P and solar zenith
! angle ANGLE.
REAL, INTENT(IN OUT) :: ano(jpwav)

! Local variables
INTEGER, PARAMETER :: jpnob=9
INTEGER, PARAMETER :: jpnoz=5

REAL :: cz0
REAL :: cz1
REAL :: sec
REAL :: sigmae0
REAL :: sigmae1
REAL :: zlogn
REAL :: zlogn2
REAL :: zlogn3
REAL :: zlogn4
REAL :: zlogp
REAL :: zlogp2
REAL :: zlogp3
REAL :: zlogp4
REAL :: zlogp5
REAL :: zlogp6
REAL :: zlogp7
REAL :: zlogp8
REAL :: c00
REAL :: c01


!     NO parameterisation Allen & Frederick 1983.
!     Polynomial coefficients for the NO effective absorption cross
!     section Zenith angle dependence.
!     Band del(0-0)
REAL,PARAMETER :: band00(jpnob) = [ -1.790868e1, -1.924701e-1, -7.217717e-2,   &
  5.648282e-2 , 4.569175e-2 , 8.353572e-3 , 0.0 , 0.0 , 0.0 ]
REAL,PARAMETER  :: zen00(jpnoz) = [7836.832, -1549.88, 114.8342, -3.777754,    &
  4.655696e-2 ]
!     Band del(0-1)
REAL,PARAMETER  :: band01(jpnob) = [ -1.654245e1, 5.836899e-1, 3.449436e-1,    &
  1.700653e-1 , -3.324717e-2, -4.952424e-2, 1.579306e-2,                       &
  1.835462e-2 , 3.368125e-3]
REAL,PARAMETER  :: zen01(jpnoz) =[12975.81 , -2582.981 , 192.7709, -6.393008,  &
  7.949835e-2 ]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSNO'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Initialise cross sections
ano(49:50)=0.0
ano(54:55)=0.0

!     Set up the two logical conditions.
!     i.e. Only continue if the parameterisation applies.

!          > 20 KM           < 84 degrees.
IF ( (p <= 50.0) .AND. (angle <= (84.0*pi_over_180)) ) THEN

  !       Set values of the NO crossection for the TWO wavelength
  !       intervals where the parameterisation applies over the
  !       altitudes which it applies, i.e. ABOVE 20 kilometres ONLY.
  !       This parameterisation was performed by Allen & Frederick
  !       using high spectral resolution calculations of the delta
  !       band predissociation of nitric oxide.

  !       N.B. It is very important that the clause (ANGLE <= (84.0 deg))
  !           is included. This limits the calculation to angles less than
  !           84 degrees. At angles greater than this the
  !           parameterisation is NOT valid.

  !        Set Log(P), and Log(N), P being Pressure, N being the VERTICAL
  !        O2 Column above the given point.
  zlogp = LOG10(p)
  zlogn = LOG10(vc)

  !        Initialise Zenith angle dependence variables for both bands.
  c00 = 0.0
  c01 = 0.0

  sec = (1/COS(angle))

  zlogp2 = zlogp*zlogp
  zlogp3 = zlogp2*zlogp
  zlogp4 = zlogp3*zlogp
  zlogp5 = zlogp4*zlogp
  zlogp6 = zlogp5*zlogp
  zlogp7 = zlogp6*zlogp
  zlogp8 = zlogp7*zlogp

  zlogn2 = zlogn*zlogn
  zlogn3 = zlogn2*zlogn
  zlogn4 = zlogn3*zlogn


  !        For BAND DEL(1-0) calculate the NO effective cross section for
  !        an OVERHEAD sun using a simple polynomial expression.
  sigmae1 = band01(1) + band01(2)*zlogp + band01(3)                            &
    *zlogp2 + band01(4)*zlogp3 + band01(5)                                     &
    *zlogp4 + band01(6)*zlogp5 + band01(7)                                     &
    *zlogp6 + band01(8)*zlogp7 + band01(9)*zlogp8

  sigmae1 = 10.0**sigmae1

  !        For BAND DEL(1-0) calculate the Zenith angle dependence of the
  !        NO effective cross section using a simple polynomial expression.
  cz1 = zen01(1) + zen01(2)*zlogn + zen01(3)*zlogn2 + zen01(4)                 &
    *zlogn3 + zen01(5)*zlogn4

  c01 = (sec)**cz1

  !        Effective NO absorption cross sections for BAND DEL(1-0)
  ano(49) = sigmae1*c01
  ano(50) = ano(49)


  !        For BAND DEL(0-0) calculate the NO effective cross section for
  !        an OVERHEAD sun using a simple polynomial expression.
  sigmae0 = band00(1) + band00(2)*zlogp + band00(3)                            &
    *zlogp2 + band00(4)*zlogp3 + band00(5)                                     &
    *zlogp4 + band00(6)*zlogp5

  sigmae0 = 10.0**sigmae0

  !        For BAND DEL(0-0) calculate the Zenith angle dependence of the
  !        NO effective cross section using a simple polynomial expression.
  cz0 = zen00(1) + zen00(2)*zlogn + zen00(3)*zlogn2 + zen00(4)                 &
    *zlogn3 + zen00(5)*zlogn4

  c00 = (sec)**cz0

  !        Effective NO absorption cross sections for BAND DEL(0-0)
  ano(54) = sigmae0*c00
  ano(55) = ano(54)

  !     End of in range if statement.
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsno
END MODULE acsno_mod
