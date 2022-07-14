! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Initialise the reaction rate co-efficients for use in the
!  Backward Euler solver with RAQ chemistry (based on STOCHEM ).
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  Fortran 95
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_chemco_raq_init_mod

IMPLICIT NONE

!
! Set up rate coefficient data type
! Rate constant k = A.(T/300)**n.exp(E/T)
! E is -Ea/R, where Ea is the activation energy and R the ideal gas constant
! Hence Ea has units of K.
! r is the reaction type:
!  0 - no reaction
!  1 - k = A  (temperature-independent)
!  2 - k = A.(T/300)**n
!  3 - k = A.exp(E/T)
!  4 - k = A.(T/300)**n.exp(E/T)
!
TYPE, PRIVATE :: rr
  REAL    :: a    ! Pre-exponential factor
  REAL    :: n    ! Power for temperature
  REAL    :: e    ! Activation Energy (actually -Ea/R)
  INTEGER :: r    ! Reaction type (see above)
END TYPE rr

! Set up type to hold information needed to calculate pressure-dependent
! rate constants
TYPE, PRIVATE :: pd
  INTEGER :: klo  ! Reaction nr holding low-pressure rate constant
  INTEGER :: khi  ! Reaction nr holding high-pressure rate constant
  INTEGER :: res  ! Rate constant to store final result in
  REAL    :: fc   ! Fc factor
END TYPE pd


! Holds A, n, Ea for each rate constant
TYPE(rr), ALLOCATABLE, SAVE, PUBLIC :: rk(:)

! Pressure-dependent reactions
TYPE(pd), ALLOCATABLE, SAVE, PUBLIC :: pdep(:)


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CHEMCO_RAQ_INIT_MOD'

CONTAINS

SUBROUTINE ukca_chemco_raq_init (nr, npdep)

! Method: Set derived type arrays holding rate constants for thermal reactions
!         (rk) as well as information relevant for the calculation of pressure-
!         dependent reaction rates (pdep).

USE ukca_config_specification_mod, ONLY: ukca_config
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umprintmgr, ONLY: umPrint, umMessage
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

INTEGER, INTENT(IN) :: nr     ! number of thermal reactions
INTEGER, INTENT(IN) :: npdep  ! number of pressure-dependent reactions

! Local variables

INTEGER :: icode           ! Error code
INTEGER :: i               ! Loop counters
INTEGER :: r               ! Reaction type
REAL :: ar                 ! Arrhenius pre-exponential factor
REAL :: n                  ! Temperature power (K)
REAL :: ea                 ! Activation energy (as -Ea/R, units: K)

CHARACTER(LEN=errormessagelength) :: cmessage  ! Error return message
CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMCO_RAQ_INIT'

TYPE(rr) :: zero_rate = rr(0.0, 0.0, 0.0, 0)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Safety checks
IF (ukca_config%l_ukca_raqaero) THEN
  IF (nr /= 197) THEN
    cmessage = 'Invalid number of thermal reactions for RAQ-AERO chemistry'
    icode = INT(nr)
    CALL ereport(routinename, icode, cmessage)
  END IF
ELSE
  IF (nr /= 192) THEN
    cmessage = 'Invalid number of thermal reactions for RAQ chemistry'
    icode = INT(nr)
    CALL ereport(routinename, icode, cmessage)
  END IF
END IF

IF (npdep /= 10) THEN
  cmessage = 'Invalid number of pressure-dependent reactions'
  icode = INT(npdep)
  CALL ereport(routinename, icode, cmessage)
END IF


ALLOCATE(rk(nr))  ! Holds A, n, Ea for each rate constant

!              A       n     -Ea/R  rtype   Num   Reaction Details

!   1 : O + O2 + M = O3 + M
rk(1)   = rr( 6.0e-34, -2.6,     0.0,  2)
rk(2)   = zero_rate
rk(3)   = zero_rate
!   4 : O + NO + M = NO2 + M  klo
rk(4)   = rr( 1.0e-31, -1.6,     0.0,  2)
!   5 : O + NO + M = NO2 + M  khi
rk(5)   = rr( 3.0e-11,  0.3,     0.0,  2)
!   6 : O(1D) + M = O(3P) + M
rk(6)   = rr( 3.2e-11,  0.0,    70.0,  3)
!   7 : O(1D) + M 2nd expression
rk(7)   = rr( 1.8e-11,  0.0,   110.0,  3)
!   8 : O(1D) + H2O = 2 OH
rk(8)   = rr( 2.2e-10,  0.0,     0.0,  1)
rk(9)   = zero_rate
rk(10)  = zero_rate

!  11 : NO + O3 = NO2 + O2
rk(11)  = rr( 1.4e-12,  0.0, -1310.0,  3)
!  12 : NO2 + O3 = NO3 + O2
rk(12)  = rr( 1.4e-13,  0.0, -2470.0,  3)
!  13 : OH + O3 = HO2 + O2
rk(13)  = rr( 1.7e-12,  0.0,  -940.0,  3)
!  14 : HO2 + O3 = OH + O2 + O2
rk(14)  = rr( 2.03e-16, 4.57, 693.0,   4)
!  15 : NO + NO3 = NO2 + NO2
rk(15)  = rr( 1.8e-11,  0.0,   110.0,  3)
!  16 : NO2 + O = NO + O2
rk(16)  = rr( 5.5e-12,  0.0,   188.0,  3)
!  17 : NO + HO2 = OH + NO2
rk(17)  = rr( 3.6e-12,  0.0,   270.0,  3)
rk(18)  = zero_rate
!  19 : NO2 + NO3 = NO + NO2 + O2
rk(19)  = rr( 4.5e-14,  0.0, -1260.0,  3)
!  20 : NO2 + NO3 + M = N2O5 + M  klo see R46
rk(20)  = rr( 3.6e-30, -4.1,     0.0,  2)

!  21 : NO2 + OH + M = HNO3 + M  klo see R47
rk(21)  = rr( 3.3e-30, -3.0,     0.0,  2)
!  22 : NO2 + HO2 + M = HO2NO2 + M  klo see R48
rk(22)  = rr( 1.8e-31, -3.2,     0.0,  2)
!  23 : HO2NO2 + M = HO2 + NO2 + M   klo see R25
rk(23)  = rr( 4.1e-5,  0.0, -10650.0,  3)
!  24 : OH + HO2NO2 = H2O + NO2 + O2
rk(24)  = rr( 3.2e-13, 0.0,    690.0,  3)
!  25 : HO2NO2 + M = HO2 + NO2 + M   khi see R23
rk(25)  = rr( 4.8e15,  0.0, -11170.0,  3)
rk(26)  = zero_rate
!  27 : NO3 + NO3 = NO2 + NO2 + O2
rk(27)  = rr( 8.5e-13,  0.0, -2450.0,  3)
!  28 : OH + NH3  = NULL (Note that products NH2 + H2O are neglected)
rk(28)  = rr( 3.5e-12,  0.0,  -925.0,  3)
!  29 : N2O5 + M = NO2 + NO3 + M  klo
rk(29)  = rr( 1.3e-03,-3.5,-11000.0,  4)

!  30 : HO2 + OH = H2O + O2
rk(30)  = rr( 4.8e-11,  0.0,   250.0,  3)
!  31 : OH + H2O2 = H2O + HO2
rk(31)  = rr( 2.9e-12,  0.0,  -160.0,  3)
!  32 : NO3 + HO2 = HNO3 + O2
rk(32)  = rr( 4.2e-12,  0.0,     0.0,  1)
!  33 : OH + H2 (+O2) = HO2 + H2O
rk(33)  = rr( 5.5e-12,  0.0, -2000.0,  3)
!  34 : NO3 + HO2 = OH + NO2 + O2
rk(34)  = rr( 3.5e-12,  0.0,     0.0,  1)
!  35 : OH + HNO3 = NO3 + H2O  See R50,R51
rk(35)  = rr( 2.4e-14,  0.0,   460.0,  3)
!  36 : HO2 + HO2 = H2O2 + O2
rk(36)  = rr( 2.2e-13,  0.0,   600.0,  3)
!  37 : HO2 + HO2 (+M) = H2O2 + O2
rk(37)  = rr( 1.9e-33,  0.0,   980.0,  3)
!  38 : HO2 + HO2 (+H2O) = H2O2 + O2
rk(38)  = rr( 1.4e-21,  0.0,  2200.0,  3)
rk(39)  = zero_rate

rk(40)  = zero_rate
rk(41)  = zero_rate
rk(42)  = zero_rate
rk(43)  = zero_rate
!  44 : OH + CH3OOH = CH3O2 + H2O
rk(44)  = rr( 2.66e-12, 0.0,   200.0,  3)
!  45 : OH + CH3OOH = OH + HCHO + H2O
rk(45)  = rr( 1.14e-12, 0.0,   200.0,  3)
!  46 : NO2 + NO3 + M = N2O5 + M  khi for R20
rk(46)  = rr( 1.9e-12,  0.2,     0.0,  2)
!  47 : NO2 + OH + M = HNO3 + M  khi see R21
rk(47)  = rr( 4.1e-11,  0.0,     0.0,  1)
!  48 : NO2 + HO2 + M = HO2NO2 + M  khi see R22
rk(48)  = rr( 4.7e-12,  0.0,     0.0,  1)
!  49 : N2O5 + M = NO2 + NO3 + M  khi see R29
rk(49)  = rr( 9.7e+14,  0.1,-11080.0,  4)

!  50 : OH + HNO3 = NO3 + H2O 2nd term for R35
rk(50)  = rr( 2.7e-17,  0.0,  2199.0,  3)
!  51 : OH + HNO3 = NO3 + H2O 3rd term for R35
rk(51)  = rr( 6.5e-34,  0.0,  1335.0,  3)
rk(52)  = zero_rate
rk(53)  = zero_rate
rk(54)  = zero_rate
rk(55)  = zero_rate
rk(56)  = zero_rate
rk(57)  = zero_rate
rk(58)  = zero_rate
!  59 : OH + CH4 = CH3O2 + H2O
rk(59)  = rr(1.85e-12,  0.0, -1690.0,  3)

!  60 : NO + CH3O2 +O2 = HCHO + HO2 + NO2
rk(60)  = rr(2.3e-12,   0.0,   360.0,  3)
!  61 : CH3O2 + CH3O2 = 2 HCHO + 2 HO2
rk(61)  = rr(7.4e-13,   0.0,  -520.0,  3)
!  62 : Total rate const: CH3O2 + CH3O2 = Products
rk(62)  = rr(1.03e-13,  0.0,   365.0,  3)
!  63 : CH3OH + OH = HO2 + HCHO + H2O
rk(63)  = rr(7.3e-12,   0.0,  -620.0,  3)
rk(64)  = zero_rate
!  65 : CH3O2 + HO2 = CH3OOH + O2
rk(65)  = rr( 4.1e-13,  0.0,   750.0,  3)
!  66 : OH + HCHO = HO2 + CO + H2O
rk(66)  = rr( 5.4e-12,  0.0,   135.0,  3)
!  67 : NO3 + HCHO = HO2 + CO + HNO3
rk(67)  = rr( 5.8e-16,  0.0,     0.0,  1)
rk(68)  = zero_rate
!  69 : OH + CO Press dependent term
rk(69)  = rr(3.54e-33,  0.0,     0.0,  1)

!  70 : OH + CO = HO2 + CO2
rk(70)  = rr( 1.5e-13,  0.0,     0.0,  1)
!  71 : OH + C2H6 = C2H5O2
rk(71)  = rr(6.9e-12,   0.0, -1000.0,  3)
!  72 : C2H5O2 + NO = CH3CHO + HO2 + NO2
rk(72)  = rr( 2.6e-12,  0.0,     0.0,  1)
!  73 : C2H5O2 + CH3O2 = CH3CHO + 2 HO2 + HCHO + O2
rk(73)  = rr( 2.0e-13,  0.0,     0.0,  1)
!  74 : branching ratio for R80
!       CH3O2 + CH3COO2 = 2 HCHO
rk(74)  = rr( 4.4e5,    0.0, -3910.0,  3)
!  75 : OH + CH3CHO = CH3COO2 + H2O
rk(75)  = rr( 4.4e-12,  0.0,   365.0,  3)
!  76 : CH3COO2 + NO2 + M = CH3COO2NO2 + M  klo
rk(76)  = rr( 2.7e-28, -7.1,     0.0,  2)
!  77 : CH3COO2 + NO2 + M = CH3COO2NO2 + M  khi
rk(77)  = rr( 1.2e-11, -0.9,     0.0,  2)
!  78 : CH3COO2NO2 + M = CH3COO2 + NO2 + M  klo
rk(78)  = rr( 4.9e-03,  0.0,-12100.0,  3)
!  79 : CH3COO2 + NO = CH3O2 + CO2 + NO2
rk(79)  = rr( 2.0e-11,  0.0,     0.0,  1)

!  80 : CH3O2 + CH3COO2 = HCHO + HO2 + CH3O2 + CO2 + O2
rk(80)  = rr( 1.1e-11,  0.0,     0.0,  1)
!  81 : OH + n-C4H10 (+O2) = s-C4H9O2 + H2O
rk(81)  = rr(7.9e-13,   2.0,   300.0,  4)
!  82 : CH3COO2NO2 + M = CH3COO2 + NO2 + M  khi
rk(82)  = rr( 5.4e+16,  0.0,-13830.0,  3)
!  83 : s-C4H9O2 + NO = CH3COC2H5 + HO2 + NO2
rk(83)  = rr( 2.54e-12, 0.0,   360.0,  3)
!  84 : s-C4H9O2 + CH3O2 = CH3COC2H5+2 HO2 +HCHO+O2
rk(84)  = rr( 2.5e-13,  0.0,     0.0,  1)
rk(85) = zero_rate
!  86 : OH + CH3COC2H5 (+O2) = MEKO2 + H2O
rk(86) = rr(1.3e-12,   0.0,   -25.0,  3)
rk(87) = zero_rate
rk(88) = zero_rate
!  89 : CH3COCH3 + OH = CH3COCH2O2 + H2O  See R94
rk(89)  = rr( 8.8e-12,  0.0, -1320.0,  3)

!  90 : C2H5O2 + C2H5O2 = C2H5O + C2H5O + O2
rk(90)  = rr( 6.4e-14,  0.0,     0.0,  1)
!  91 : CH3COO2 + CH3COO2 = 2 CH3O2 + 2 CO2 + O2
rk(91)  = rr( 2.9e-12,  0.0,   500.0,  3)
!  92 : OH + C3H8 = i-C3H7O2 + H2O
rk(92)  = rr( 7.6e-12,  0.0,  -585.0,  3)
!  93 : i-C3H7O2 + NO = NO2 + HO2 + CH3COCH3
rk(93)  = rr( 2.7e-12 , 0.0,   360.0,  3)
!  94 : CH3COCH3 + OH = CH3COCH2O2 + H2O  2nd term for R89
rk(94)  = rr( 1.7e-14,  0.0,   423.0,  3)
!  95 : CH3COCH2O2 + NO = NO2 + CH3COO2 + HCHO
rk(95)  = rr(2.45e-12,  0.0,   360.0,  3)
!  96 : CH3COCH2O2 + CH3O2 = HO2 + 2 HCHO + CH3COO2
rk(96)  = rr( 3.8e-12,  0.0,     0.0,  1)
!  97 : i-C3H7O2 + CH3O2 = HCHO + 2 HO2 + CH3COCH3
rk(97)  = rr( 4.0e-14,  0.0,     0.0,  1)
!  98 : OH + PAN = NO3 + HCHO
rk(98)  = rr( 9.5e-13,  0.0,  -650.0,  3)
!  99 : HO2 + C2H5O2 = C2H5OOH + O2
rk(99)  = rr( 3.8e-13,  0.0,   900.0,  3)

! 100 : OH + C2H5OOH = CH3CHO + OH
rk(100) = rr( 8.0e-12,  0.0,     0.0,  1)
! 101 : HO2 + C3H7O2 = i-C3H7OOH + O2
rk(101) = rr(1.51e-13,  0.0,  1300.0,  3)
! 102 : OH + i-C3H7OOH = CH3COCH3 + OH
rk(102) = rr(1.66e-11,  0.0,     0.0,  1)
! 103 : HO2 + s-C4H9O2 = s-C4H9OOH + O2
rk(103) = rr(1.82e-13,  0.0,  1300.0,  3)
! 104 : OH + s-C4H9OOH = CH3COC2H5 + OH
rk(104) = rr(2.15e-11,  0.0,     0.0,  1)
! 105 : NO + CH3COCH(O2)CH3 = CH3COO2 + CH3CHO + NO2
rk(105) = rr(2.54e-12,  0.0,   360.0,  3)
! 106 : CH3O2 + CH3COCH(O2)CH3 =
!                    HCHO + HO2 + CH3COO2 + CH3CHO + O2
rk(106) = rr( 8.8e-13,  0.0,     0.0,  1)
rk(107) = zero_rate
! 108 : OH + C2H4 (+O2) = HOC2H4O2  klo
rk(108)  = rr( 8.6e-29, -3.1,     0.0,  2)
! 109 : OH + C2H4 (+O2) = HOC2H4O2  khi
rk(109)  = rr( 9.0e-12,-0.85,     0.0,  2)

! 110 : HOC2H4O2 + NO = 2 HCHO + HO2 + NO2
rk(110)  = rr( 9.0e-12,  0.0,     0.0,  1)
! 111 : CH3O2 + HOC2H4O2 = 3 HCHO + 2 HO2 + O2
rk(111)  = rr( 2.0e-12,  0.0,     0.0,  1)
! 112 : O3 + C2H4 = HCHO + 0.47 CH2O2 + 0.31 CO
!                 + 0.22 CO2 + 0.31 H2O + 0.13 H2 + 0.20 HO2
rk(112)  = rr( 1.2e-14,  0.0, -2630.0,  3)
rk(113) = zero_rate
rk(114) = zero_rate
rk(115) = zero_rate
rk(116) = zero_rate
rk(117) = zero_rate
rk(118) = zero_rate
rk(119) = zero_rate

rk(120) = zero_rate
rk(121) = zero_rate
rk(122) = zero_rate
! 123 : O3 + C3H6 = HCHO + 0.30 CH4 + 0.40 CO
!    + 0.60 CO2 + 0.28 OH + 0.12 CH3OH + 0.30 HO2 + 0.58 CH3O2
rk(123) = rr(2.75e-15,  0.0, -1878.0,  3)
! 124 : O3 + C3H6 = CH3CHO + 0.24 H2 + 0.58 CO
!          + 0.42 CO2 + 0.58 H2O +0.18 HO2
rk(124) = rr(2.75e-15,  0.0, -1878.0,  3)
! 125 : OH + C3H6 (+O2) = HOC3H6O2  klo see R130
rk(125)= rr( 8.0e-27, -3.5,     0.0,  2)
! 126 : HOC3H6O2 + NO = HCHO + HO2 + CH3CHO + NO2
rk(126) = rr(2.54e-12,  0.0,   360.0,  3)
! 127 : CH3O2 + HOC3H6O2 = 2 HCHO + 2 HO2 + CH3CHO + O2
rk(127) = rr( 6.0e-13,  0.0,     0.0,  1)
! 128 : O3 + C5H8 = MVK + 0.78 CO + 0.22 CH2OO + 0.27 HO2
!                       + 0.27 OH
rk(128) = rr(7.86e-15,  0.0, -1913.0,  3)
! 129 : O3 + MVK = MGLYOX + 0.76 CO + 0.24 CH2OO + 0.36 HO2
!                         + 0.36 OH
rk(129) = rr(7.56e-16,  0.0, -1521.0,  3)

! 130 : OH + C3H6 (+O2) = HOC3H6O2  khi see R125
rk(130) = rr( 3.0e-11,  -1.0,    0.0,  2)
! 131 : CH3CO3 + HO2 = 0.3 O3 + 0.8CH3O2 + 0.2CH3COO2
rk(131) = rr( 5.2e-13,  0.0,   980.0,  3)
! 132 : HOIPO2 + CH3O2 = 2HO2 + HCHO + MVK
rk(132) = rr( 5.0e-13,  0.0,     0.0,  1)
! 133 : HOMVKO2+ CH3O2 = 2HO2 + HCHO + MGLYOX
rk(133) = rr( 2.0e-12,  0.0,     0.0,  1)
! 134 : HOIPO2 + HO2 = ISOOH + O2
rk(134) = rr(2.45e-13,  0.0,  1250.0,  3)
! 135 : ISOOH  + OH = MVK + HCHO + OH
rk(135)  = rr( 4.2e-11,  0.0,     0.0,  1)
! 136 : HOMVKO2 + HO2 = MVKOOH + O2
rk(136) = rr(2.23e-13,  0.0,  1250.0,  3)
! 137 : MVKOOH + OH = MGLYOX + HCHO + OH
rk(137) = rr(5.77e-11,  0.0,     0.0,  1)
! 138 : MGLYOX + OH = CH3COO2 + CO
rk(138) = rr(1.72e-11,  0.0,     0.0,  1)
! 139 : GLYOX + OH = HO2 + 2 CO
rk(139) = rr(1.14e-11,  0.0,     0.0,  1)

rk(140) = zero_rate
! 141 : OH + C5H8 = HOIPO2 + H2O
rk(141) = rr(2.54e-11,  0.0,   410.0,  3)
! 142 : HOC5H8O2 + NO = MVK + HO2 + HCHO + NO2
rk(142) = rr(2.08e-12,  0.0,   180.0,  3)
! 143 : OH + MVK = HOMVKO2 + H2O
rk(143) = rr(4.13e-12,  0.0,   452.0,  3)
! 144 : HOMVKO2 + NO = CH3COCHO + CH2O + HO2 + NO2
rk(144) = rr(2.5e-12,   0.0,   360.0,  3)
rk(145) = zero_rate
rk(146) = zero_rate
rk(147) = zero_rate
rk(148) = zero_rate
rk(149) = zero_rate

rk(150) = zero_rate
! 151 : NO3 + C2H6 = C2H5O2 + HNO3
rk(151)  = rr( 5.7e-12,  0.0, -4426.0,  3)
! 152 : NO3 + n-C4H10 = s-C4H9O2 + HNO3
rk(152) = rr(2.8e-12,   0.0, -3280.0,  3)
! 153 : NO3 + C2H4 = C2H4NO3 = CH2(NO3)CHO + HO2
rk(153) = rr(3.3e-12,   2.0, -2880.0,  4)
! 154 : CH2(NO3)CHO + OH = HCHO + NO2 + CO2
rk(154) = rr(4.95e-12,  0.0,     0.0,  1)
! 155 : NO3 + C3H6 = C3H6NO3 = CH3CH(NO3)CHO + HO2
rk(155) = rr(4.59e-13,  0.0, -1156.0,  3)
! 156 : CH3CH(NO3)CHO + OH = CH3CHO + NO2 + CO2
rk(156) = rr(5.25e-12,  0.0,     0.0,  1)
rk(157) = zero_rate
! 158 : NO3 + CH3CHO = CH3COO2 + HNO3
rk(158) = rr( 1.4e-12,  0.0, -1860.0,  3)
! 159 : NO3 + C5H8 = (NO3)C4H6CHO + HO2
rk(159) = rr(3.03e-12,  0.0,  -446.0,  3)

! 160 : (NO3)C4H6CHO + OH = MVK + NO2
rk(160) = rr(4.16e-11,  0.0,     0.0,  1)
rk(161) = zero_rate
rk(162) = zero_rate
rk(163) = zero_rate
rk(164) = zero_rate
rk(165) = zero_rate
rk(166) = zero_rate
rk(167) = zero_rate
rk(168) = zero_rate
rk(169) = zero_rate

! 170 : OH + oXYL = HO2 + 0.8 MEMALD + 0.8 MGLYOX
rk(170) = rr(1.36e-11,  0.0,     0.0,  1)
! 171 : OXYL1 + NO2 = ORGNIT
rk(171) = rr( 1.0e-11,  0.0,     0.0,  1)
! 172 : OH + MEMALD = MEMALD1
rk(172) = rr( 5.6e-11,  0.0,     0.0,  1)
! 173 : MEMALD1 + NO = HO2 + CH3COCHO + CHOCHO + NO2
rk(173) = rr( 2.54e-12, 0.0,   360.0,  3)
! 174 : OH + TOLUENE = MEMALD + GLYOX + HO2
rk(174) = rr(1.18e-12,  0.0,   338.0,  3)
! 175 : OH + TOLUENE = TOLP1
rk(175)  = rr( 3.6e-13,  0.0,     0.0,  1)
! 176 : OH + oXYL    = OXYL1
rk(176) = rr(1.36e-11,  0.0,     0.0,  1)
! 177 : OH + ORGNIT  = MEMALD + GLYOX + NO2
rk(177) = rr( 2.7e-12,  0.0,     0.0,  1)
! 178 : NO3 + ORGNIT = MEMALD + GLYOX + 2NO2
rk(178)  = rr( 7.0e-14,  0.0,     0.0,  1)
rk(179) = zero_rate

! 180 : OXYL1 + HO2 = MGLYOX + MEMALD
rk(180)  = rr( 2.5e-13,  0.0,  1300.0,  3)
! 181 : MEMALD1 + CH3O2 = 2HO2 + HCHO + MGLYOX + GLYOX
rk(181)  = rr( 1.0e-13,  0.0,     0.0,  1)
! 182 : HO2 + TOLP1 = MEMALD + GLYOX + OH
rk(182)  = rr( 1.0e-11,  0.0,     0.0,  1)
! 183 : NO2 + TOLP1 = ORGNIT
rk(183)  = rr( 1.0e-11,  0.0,     0.0,  1)
rk(184) = zero_rate ! reserved for RAQ-AERO: DMS + OH
rk(185) = zero_rate ! reserved for RAQ-AERO: DMS + OH
rk(186) = zero_rate ! reserved for RAQ-AERO: DMS + NO3
rk(187) = zero_rate ! reserved for RAQ-AERO: DMSO + OH
rk(188) = zero_rate ! reserved for RAQ-AERO: SO2 + OH
rk(189) = zero_rate

rk(190) = zero_rate ! reserved for RAQ-AERO: SO2 + O3
rk(191) = zero_rate ! reserved for RAQ-AERO: SO2 + H2O2
rk(192) = zero_rate ! reserved for RAQ-AERO: Monoterp + OH

! The additional RAQ-AERO rates are defined in ukca_chemco_raq; we initialise
! to zero here to enable the error checking on rk(i) to pass.
IF (ukca_config%l_ukca_raqaero) THEN
  rk(193:197) = zero_rate
END IF

! Check that reaction types are correct
DO i = 1, nr
  r = rk(i)%r
  ar = rk(i)%a
  n = rk(i)%n
  ea = rk(i)%e
  IF (    ((r == 1) .AND. (n /= 0.0 .OR. ea /= 0.0))      .OR.                 &
          ((r == 2) .AND. (n == 0.0 .OR. ea /= 0.0))      .OR.                 &
          ((r == 3) .AND. (n /= 0.0 .OR. ea == 0.0))      .OR.                 &
          ((r == 4) .AND. (n == 0.0 .OR. ea == 0.0))   ) THEN
    cmessage = 'Faulty rate constant'
    WRITE(umMessage,'(A,A,I8,A,I8,A,F16.8,A,F16.8,A,F16.8)')                   &
                                      'Faulty rate constant',                  &
                                      ' i=', i,                                &
                                      ' r=', r,                                &
                                      ' ar=', ar,                              &
                                      ' n=', n,                                &
                                      ' ea=', ea
    CALL umPrint(umMessage, src=routinename)
    icode = i
    CALL ereport(routinename, icode, cmessage)
  END IF
END DO

! Pressure-dependent reactions
ALLOCATE(pdep(npdep))

!              klo khi res  Fc
! O + NO + M = NO2 + M
pdep( 1) =  pd(  4,  5,  5, 0.85)

! NO2 + NO3 + M = N2O5 + M
pdep( 2) =  pd( 20, 46, 20, 0.35)

! NO2 + OH + M = HNO3 + M
pdep( 3) =  pd( 21, 47, 21, 0.4 )

! NO2 + HO2 + M = HO2NO2 + M
pdep( 4) =  pd( 22, 48, 22, 0.6 )

! HO2NO2 + M = HO2 + NO2 + M
pdep( 5) =  pd( 23, 25, 23, 0.6 )

! N2O5 + M = NO2 + NO3 + M
pdep( 6) =  pd( 29, 49, 29, 0.35)

! CH3COO2 + NO2 + M = CH3COO2NO2 (PAN) + M
pdep( 7) =  pd( 76, 77, 77, 0.3 )

! CH3COO2NO2 + M = CH3COO2 + NO2 + M
pdep( 8) =  pd( 78, 82, 78, 0.3 )

! OH + C2H4 + M (+O2) = HOC2H4O2 + M
pdep( 9) =  pd(108,109,109, 0.48)

! OH + C3H6 + M (+O2) = HOC3H6O2 + M
pdep(10) =  pd(125,130,125, 0.5)



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chemco_raq_init

END MODULE ukca_chemco_raq_init_mod
