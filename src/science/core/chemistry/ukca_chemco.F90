! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  To calculate the reaction rate co-efficients for use in the
!  Backward Euler solver
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method: Calculates rate coefficients
!                 Units : cm, molecule, s
!
!   Inputs  : tc,m,h2o,o2
!   Outputs : rc
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_chemco_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CHEMCO_MOD'

CONTAINS

SUBROUTINE ukca_chemco(nr,n_pnts,tc,m,h2o,o2,                                  &
                       clw,fcloud,fdiss,k_dms,rc,H_plus_1d_arr)

USE water_constants_mod, ONLY: rho_water
USE chemistry_constants_mod, ONLY: avogadro
USE ukca_constants,   ONLY: H_plus, m_air
USE ukca_config_specification_mod, ONLY: ukca_config
USE asad_mod,         ONLY: peps, jpeq, jpspec
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE


INTEGER, INTENT(IN) :: nr                ! No of reactions
INTEGER, INTENT(IN) :: n_pnts            ! No of points

REAL, INTENT(IN)    :: tc(n_pnts)        ! Temperature
REAL, INTENT(IN)    :: m(n_pnts)         ! Air density
REAL, INTENT(IN)    :: h2o(n_pnts)       ! Water vapour
REAL, INTENT(IN)    :: o2(n_pnts)        ! Oxygen
REAL, INTENT(IN)    :: clw(n_pnts)       ! Cloud Liquid water
REAL, INTENT(IN)    :: fcloud(n_pnts)    ! Cloud fraction
! Dissolved fraction - allows for jpeq+1 ions:
REAL, INTENT(IN)    :: fdiss(n_pnts,jpspec,jpeq+1)

! 1-D pH array to be used in calculating rate constants
REAL, INTENT(IN) :: H_plus_1d_arr(n_pnts)

! Rate coeffs to allow parameterisation of DMS products:
REAL, INTENT(OUT)   :: k_dms(n_pnts,5)

!     Thermal rate coefficients:
REAL, INTENT(OUT)   :: rc(n_pnts,nr)

!     Local variables

INTEGER :: i, j

REAL :: vr(n_pnts)                ! Volume ratio of cloud water
REAL :: z1(n_pnts)
REAL :: z2(n_pnts)
REAL :: z3(n_pnts)
REAL :: z4(n_pnts)
REAL :: ratiob2total(n_pnts)
REAL :: ratioa2b(n_pnts)
REAL :: at(7)
REAL :: t300(n_pnts)
REAL :: zo(n_pnts), zi(n_pnts), arg2(n_pnts)
REAL :: zfc(n_pnts), zr(n_pnts)
REAL :: faq(n_pnts,jpspec)           ! total dissolved fraction

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMCO'



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ukca_config%l_ukca_aerchem) THEN
  faq(:,:) = fdiss(:,:,1) + fdiss(:,:,2)
  ! convert clw in kg/kg to volume ratio
  vr(:) = clw(:)*m(:)*m_air*(1e6/avogadro)/rho_water
END IF

!     Initialise reaction rate coefficients

DO j = 1, nr
  DO i = 1, n_pnts
    rc(i,j) = 0.0
  END DO
END DO

!     Calculate bimolecular reaction rate coefficients

!     R1: HO2 + NO = OH + NO2              3.60E-12  0.00   -270.0
!     IUPAC no change

rc(:,1) = 3.60e-12*EXP(270.0/tc(:))

!     R2: HO2 + NO3 = OH + NO2             4.00E-12  0.00      0.0
!     IUPAC no change

rc(:,2) = 4.00e-12

!     R3: HO2 + O3 = OH + O2               2.03E-16  4.57   -693.0
!     IUPAC no change

rc(:,3) = 2.03e-16*(tc(:)/300.0)**4.57*EXP(693.0/tc(:))

!     R4: HO2 + HO2 = H2O2                 2.20E-13  0.00   -600.0
!     IUPAC no change - Note 4

rc(:,4) = 2.20e-13*EXP(600.0/tc(:))
rc(:,4) = rc(:,4)*(1.0 +                                                       &
          1.4e-21*h2o(:)*EXP(2200.0/tc(:)) )

!     R5: HO2 + MeOO = MeOOH               3.80E-13  0.00   -780.0
!     IUPAC 2005 update - Note 5 - Brnch A

ratiob2total(:) = 1.0 / (1.0 + 498.0*EXP(-1160.0/tc(:)))
rc(:,5)         = 3.80e-13*EXP(780.0/tc(:))
rc(:,5)         = rc(:,5)*(1.0-ratiob2total(:))

!     R6: HO2 + MeOO = HCHO                3.80E-13  0.00   -780.0
!     IUPAC 2005 update - Note 6 - Brnch B

rc(:,6) = 3.80e-13*EXP(780.0/tc(:))
rc(:,6) = rc(:,6)*ratiob2total(:)

!     R7: HO2 + EtOO = EtOOH               3.80E-13  0.00   -900.0
!     IUPAC no change

rc(:,7) = 3.80e-13*EXP(900.0/tc(:))

!     R8: HO2 + MeCO3 = MeCO3H             2.08E-13  0.00   -980.0
!     IUPAC 2005 update - Note 8

rc(:,8) = 2.08e-13*EXP(980.0/tc(:))

!     R9: HO2 + MeCO3 = MeCO2H + O3        1.04E-13  0.00   -980.0
!     IUPAC 2005 update - note 9

rc(:,9) = 1.04e-13*EXP(980.0/tc(:))

!     R10 HO2 + MeCO3 = OH + MeOO          2.08E-13  0.00   -980.0
!     IUPAC 2005 update - Note 10

rc(:,10) = 2.08e-13*EXP(980.0/tc(:))

!     R11 HO2 + n-PrOO = n-PrOOH           1.51E-13  0.00  -1300.0
!     MCM no change

rc(:,11) = 1.51e-13*EXP(1300.0/tc(:))

!     R12 HO2 + i-PrOO = i-PrOOH           1.51E-13  0.00  -1300.0
!     MCM no change

rc(:,12) = 1.51e-13*EXP(1300.0/tc(:))

!     R13 HO2 + EtCO3 = O2 + EtCO3H        3.05E-13  0.00  -1040.0
!     MCM no change

rc(:,13) = 3.05e-13*EXP(1040.0/tc(:))

!     R14 HO2 + EtCO3 = O3 + EtCO2H        1.25E-13  0.00  -1040.0
!     MCM no change

rc(:,14) = 1.25e-13*EXP(1040.0/tc(:))

!     R15 HO2 + MeCOCH2OO = MeCOCH2OOH     1.36E-13  0.00  -1250.0
!     MCM no change - Note 15

rc(:,15) = 1.36e-13*EXP(1250.0/tc(:))

!     R16 MeOO + NO = HO2 + HCHO + NO2     2.95E-12  0.00   -285.0
!     IUPAC no change

rc(:,16) = 2.95e-12*EXP(285.0/tc(:))

!     R17 MeOO + NO = MeONO2               2.95E-15  0.00   -285.0
!     IUPAC no change - note 17

rc(:,17) = 2.95e-15*EXP(285.0/tc(:))

!     R18 MeOO + NO3 = HO2 + HCHO + NO2    1.30E-12  0.00      0.0
!     IUPAC no change

rc(:,18) = 1.30e-12

!     R19 MeOO + MeOO = MeOH + HCHO        1.03E-13  0.00   -365.0
!     IUPAC no change - note 19 - Brnch A

ratiob2total(:) = 1.0/(1.0+(EXP(1300.0/tc(:)))/33.0)
rc(:,19)        = 1.03e-13*EXP(365.0/tc(:))
rc(:,19)        = rc(:,19)*(1.0-ratiob2total(:))

!     R20 MeOO + MeOO = 2HO2 + 2HCHO       1.03E-13  0.00   -365.0
!     IUPAC no change - note 20 - Brnch B

rc(:,20) = 1.03e-13*EXP(365.0/tc(:))
rc(:,20) = rc(:,20)*ratiob2total(:)

!     R21 MeOO + MeCO3 = HO2 + HCHO + MeOO 1.80E-12  0.00   -500.0
!     IUPAC no change - note 21

rc(:,21) = 1.80e-12*EXP(500.0/tc(:))

!     R22 MeOO + MeCO3 = MeCO2H + HCHO     2.00E-13  0.00   -500.0
!     IUPAC no change - note 22

rc(:,22) = 2.00e-13*EXP(500.0/tc(:))

!     R23 EtOO + NO = MeCHO + HO2 + NO2    2.60E-12  0.00   -380.0
!     IUPAC no change - note 23

rc(:,23) = 2.60e-12*EXP(380.0/tc(:))

!     R24 EtOO + NO3 = MeCHO + HO2 + NO2   2.30E-12  0.00      0.0
!     IUPAC no change

rc(:,24) = 2.30e-12

!     R25 EtOO + MeCO3 = MeCHO + HO2 + MeOO 4.40E-13  0.00  -1070.0
!     IUPAC no change - note 25

rc(:,25) = 4.40e-13*EXP(1070.0/tc(:))

!     R26 MeCO3 + NO = MeOO + CO2 + NO2    7.50E-12  0.00   -290.0
!     IUPAC no change

rc(:,26) = 7.50e-12*EXP(290.0/tc(:))

!     R27 MeCO3 + NO3 = MeOO + CO2 + NO2   4.00E-12  0.00      0.0
!     MCM no change

rc(:,27) = 4.00e-12

!     R28 n-PrOO + NO = EtCHO + HO2 + NO2  2.90E-12  0.00   -350.0
!     IUPAC no change - note 28

rc(:,28) = 2.90e-12*EXP(350.0/tc(:))

!     R29 n-PrOO + NO3 = EtCHO + HO2 + NO2 2.50E-12  0.00      0.0
!     MCM no change

rc(:,29) = 2.50e-12

!     R30 i-PrOO + NO = Me2CO + HO2 + NO2  2.70E-12  0.00   -360.0
!     IUPAC no change - note 30

rc(:,30) = 2.70e-12*EXP(360.0/tc(:))

!     R31 i-PrOO + NO3 = Me2CO + HO2 + NO2 2.50E-12  0.00      0.0
!     MCM no change

rc(:,31) = 2.50e-12

!     R32 EtCO3 + NO = EtOO + CO2 + NO2    6.70E-12  0.00   -340.0
!     IUPAC no change

rc(:,32) = 6.70e-12*EXP(340.0/tc(:))

!     R33 EtCO3 + NO3 = EtOO + CO2 + NO2   4.00E-12  0.00      0.0
!     MCM no change

rc(:,33) = 4.00e-12

!     R34 MeCOCH2OO + NO = MeCO3+HCHO+NO2  2.80E-12  0.00   -300.0
!     Tyndall et al. - note 34

rc(:,34) = 2.80e-12*EXP(300.0/tc(:))

!     R35 MeCOCH2OO + NO3 = MeCO3+HCHO+NO2 2.50E-12  0.00      0.0
!     MCM no change

rc(:,35) = 2.50e-12

!     R36 NO + NO3 = NO2 + NO2             1.80E-11  0.00   -110.0
!     IUPAC no change

rc(:,36) = 1.80e-11*EXP(110.0/tc(:))

!     R37 NO + O3 = NO2                    1.40E-12  0.00   1310.0
!     IUPAC no change

rc(:,37) = 1.40e-12*EXP(-1310.0/tc(:))

!     R38 NO2 + O3 = NO3                   1.40E-13  0.00   2470.0
!     IUPAC no change

rc(:,38) = 1.40e-13*EXP(-2470.0/tc(:))

!     R39 NO3 + HCHO = HONO2 + HO2 + CO    2.00E-12  0.00   2440.0
!     IUPAC no change - note 39

rc(:,39) = 2.00e-12*EXP(-2440.0/tc(:))

!     R40 NO3 + MeCHO = HONO2 + MeCO3      1.40E-12  0.00   1860.0
!     IUPAC no change

rc(:,40) = 1.40e-12*EXP(-1860.0/tc(:))

!     R41 NO3 + EtCHO = HONO2 + EtCO3      3.46E-12  0.00   1862.0
!     MCM no change - note 41

rc(:,41) = 3.46e-12*EXP(-1862.0/tc(:))

!     R42 NO3 + Me2CO = HONO2 + MeCOCH2OO  3.00E-17  0.00      0.0
!     IUPAC no change - note 42

rc(:,42) = 3.00e-17

!     R43 N2O5 + H2O = HONO2 + HONO2       2.50E-22  0.00      0.0
!     IUPAC no change - note 43

rc(:,43) = 2.50e-22

!     R44 O(3P) + O3 = O2 + O2             8.00E-12  0.00   2060.0
!     IUPAC no change

rc(:,44) = 8.00e-12*EXP(-2060.0/tc(:))

!     R45 O(1D) + CH4 = OH + MeOO          1.05E-10  0.00      0.0
!     IUPAC no change - note 45

rc(:,45) = 1.05e-10

!     R46 O(1D) + CH4 = HCHO + H2          7.50E-12  0.00      0.0
!     IUPAC no change - note 46

rc(:,46) = 7.50e-12

!     R47 O(1D) + CH4 = HCHO + HO2 + HO2   3.45E-11  0.00      0.0
!     IUPAC no change - note 47

rc(:,47) = 3.45e-11

!     R48 O(1D) + H2O = OH + OH            2.20E-10  0.00      0.0
!     IUPAC no change

rc(:,48) = 2.20e-10

!     R49 O(1D) + N2 = O(3P) + N2          2.10E-11  0.00   -115.0
!     Ravishankara et al. - note 49

rc(:,49) = 2.10e-11*EXP(115.0/tc(:))

!     R50 O(1D) + O2 = O(3P) + O2          3.20E-11  0.00    -67.0
!     IUPAC no change

rc(:,50) = 3.20e-11*EXP(67.0/tc(:))

!     R51 OH + CH4 = H2O + MeOO            1.85E-12  0.00   1690.0
!     IUPAC no change

rc(:,51) = 1.85e-12*EXP(-1690.0/tc(:))

!     R52 OH + C2H6 = H2O + EtOO           6.90E-12  0.00   1000.0
!     IUPAC no change

rc(:,52) = 6.90e-12*EXP(-1000.0/tc(:))

!     R53 OH + C3H8 = n-PrOO + H2O         7.60E-12  0.00    585.0
!     IUPAC no change - note 53 - Brnch A

ratioa2b(:) = 226.0 * tc(:)**(-0.64)*EXP(-816.0/tc(:))
rc(:,53) = 7.60e-12*EXP(-585.0/tc(:))
rc(:,53) = rc(:,53)*(ratioa2b(:)/(ratioa2b(:)+1.0))

!     R54 OH + C3H8 = i-PrOO + H2O         7.60E-12  0.00    585.0
!     IUPAC no change - note 54 - Brnch B

rc(:,54) = 7.60e-12*EXP(-585.0/tc(:))
rc(:,54) = rc(:,54)/(ratioa2b(:)+1.0)

!     R55 OH + CO = HO2                    1.44E-13  0.00      0.0
!     IUPAC 2005 update - note 55

rc(:,55) = 1.44e-13*(1.0+m(:)/4.2e19)

!     R56 OH + EtCHO = H2O + EtCO3         5.10E-12  0.00   -405.0
!     IUPAC no change

rc(:,56) = 5.10e-12*EXP(405.0/tc(:))

!     R57 OH + EtOOH = H2O + MeCHO + OH    8.01E-12  0.00      0.0
!     MCM no change

rc(:,57) = 8.01e-12

!     R58 OH + EtOOH = H2O + EtOO          1.90E-12  0.00   -190.0
!     MCM no change

rc(:,58) = 1.90e-12*EXP(190.0/tc(:))

!     R59 OH + H2 = H2O + HO2              7.70E-12  0.00   2100.0
!     IUPAC no change

rc(:,59) = 7.70e-12*EXP(-2100.0/tc(:))

!     R60 OH + H2O2 = H2O + HO2            2.90E-12  0.00    160.0
!     IUPAC no change

rc(:,60) = 2.90e-12*EXP(-160.0/tc(:))

IF (ukca_config%l_ukca_aerchem) THEN

  !       Reduce rate to account for dissolved H2O2
  rc(:,60) = rc(:,60)*(1.0-(faq(:,13)*fcloud(:)))

END IF

!     R61 OH + HCHO = H2O + HO2 + CO       5.40E-12  0.00   -135.0
!     IUPAC 2004 update - note 61

rc(:,61) = 5.40e-12*EXP(135.0/tc(:))

!     R62 OH + HO2 = H2O                   4.80E-11  0.00   -250.0
!     IUPAC no change

rc(:,62) = 4.80e-11*EXP(250.0/tc(:))

!     R63 OH + HO2NO2 = H2O + NO2          1.90E-12  0.00   -270.0
!     IUPAC no change

rc(:,63) = 1.90e-12*EXP(270.0/tc(:))

!     R64 OH + HONO2 = H2O + NO3           1.50E-13  0.00      0.0
!     IUPAC no change - note 64

z1(:)    = 2.4e-14 * EXP( 460.0/tc(:))
z3(:)    = 6.5e-34 * EXP(1335.0/tc(:))
z4(:)    = 2.7e-17 * EXP(2199.0/tc(:))
z2(:)    = z3(:)*m(:) / (1.0+z3(:)*m(:)/z4(:))
rc(:,64) = z1(:) + z2(:)

!     R65 OH + HONO = H2O + NO2            2.50E-12  0.00   -260.0
!     IUPAC no change

rc(:,65) = 2.50e-12*EXP(260.0/tc(:))

!     R66 OH + MeOOH = H2O + HCHO + OH     1.02E-12  0.00   -190.0
!     IUPAC no change - note 66

rc(:,66) = 1.02e-12*EXP(190.0/tc(:))

!     R67 OH + MeOOH = H2O + MeOO          1.89E-12  0.00   -190.0
!     IUPAC no change - note 67

rc(:,67) = 1.89e-12*EXP(190.0/tc(:))

!     R68 OH + MeONO2 = HCHO + NO2 + H2O   4.00E-13  0.00    845.0
!     IUPAC no change

rc(:,68) = 4.00e-13*EXP(-845.0/tc(:))

!     R69 OH + Me2CO = H2O + MeCOCH2OO     8.80E-12  0.00   1320.0
!     IUPAC no change - note 69

rc(:,69) = 8.80e-12*EXP(-1320.0/tc(:))

!     R70 OH + Me2CO = H2O + MeCOCH2OO     1.70E-14  0.00   -420.0
!     IUPAC no change - note 70

rc(:,70) = 1.70e-14*EXP(420.0/tc(:))

!     R71 OH + MeCOCH2OOH = H2O+MeCOCH2OO  1.90E-12  0.00   -190.0
!     MCM no change - note 71

rc(:,71) = 1.90e-12*EXP(190.0/tc(:))

!     R72 OH + MeCOCH2OOH = OH + MGLY      8.39E-12  0.00      0.0
!     MCM no change - note 72

rc(:,72) = 8.39e-12

!     R73 OH + MeCHO = H2O + MeCO3         4.40E-12  0.00   -365.0
!     IUPAC no change

rc(:,73) = 4.40e-12*EXP(365.0/tc(:))

!     R74 OH + NO3 = HO2 + NO2             2.00E-11  0.00      0.0
!     IUPAC no change

rc(:,74) = 2.00e-11

!     R75 OH + O3 = HO2 + O2               1.70E-12  0.00    940.0
!     IUPAC no change

rc(:,75) = 1.70e-12*EXP(-940.0/tc(:))

!     R76 OH + OH = H2O + O(3P)            6.31E-14  2.60   -945.0
!     IUPAC no change - note 76

rc(:,76) = 6.31e-14*((tc(:)/300.0)**2.6)*EXP(945.0/tc(:))

!     R77 OH + PAN = HCHO + NO2 + H2O      3.00E-14  0.00      0.0
!     IUPAC no change - note 77

rc(:,77) = 3.00e-14

!     R78 OH + PPAN = MeCHO + NO2 + H2O    1.27E-12  0.00      0.0
!     MCM no change

rc(:,78) = 1.27e-12

!     R79 OH + n-PrOOH = n-PrOO + H2O      1.90E-12  0.00   -190.0
!     MCM no change

rc(:,79) = 1.90e-12*EXP(190.0/tc(:))

!     R80 OH + n-PrOOH = EtCHO + H2O + OH  1.10E-11  0.00      0.0
!     MCM no change

rc(:,80) = 1.10e-11

!     R81 OH + i-PrOOH = i-PrOO + H2O      1.90E-12  0.00   -190.0
!     MCM no change

rc(:,81) = 1.90e-12*EXP(190.0/tc(:))

!     R82 OH + i-PrOOH = Me2CO + OH        1.66E-11  0.00      0.0
!     MCM no change

rc(:,82) = 1.66e-11

!     R83 O(3P) + NO2 = NO + O2            5.50E-12  0.00   -188.0
!     IUPAC no change

rc(:,83) = 5.50e-12*EXP(188.0/tc(:))

!     R84 HO2S + O3S = HO2S + O2           2.03E-16  4.57   -693.0
!     IUPAC no change - note 84

rc(:,84) = 2.03e-16*((tc(:)/300.0)**4.57)*EXP(693.0/tc(:))

!     R85 OHS + O3S = OHS + O2             1.70E-12  0.00    940.0
!     IUPAC no change

rc(:,85) = 1.70e-12*EXP(-940.0/tc(:))

!     R86 O(1D)S + H2O = H2O               2.20E-10  0.00      0.0
!     IUPAC no change

rc(:,86) = 2.20e-10

!     R87 O(1D)S + N2 = O(3P)S + N2        2.10E-11  0.00   -115.0
!     Ravishankara et al. - note 49

rc(:,87) = 2.10e-11*EXP(115.0/tc(:))

!     R88 O(1D)S + O2 = O(3P)S + O2        3.20E-11  0.00    -67.0
!     IUPAC no change

rc(:,88) = 3.20e-11*EXP(67.0/tc(:))

! R89 HO2 + HO2 = H2O2 + O2     0.00 1.90E-33  0.00   -980.0 0.00E+00  0.00  0.0
!     IUPAC no change - note 1

at(1) = 0.00
at(2) = 1.90e-33
at(3) = 0.00
at(4) = -980.0
at(5) = 0.00e+00
at(6) = 0.00
at(7) = 0.0

t300(:) = tc(:)/300.0
zo(:)   = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:)   = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,89) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,89) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)    = zo(i) / zi(i)
    arg2(i)  = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,89) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO
rc(:,89) = rc(:,89)*(1.0 + 1.4e-21*h2o(:)*EXP(2200.0/tc(:)) )


!     R90 HO2 + NO2 = HO2NO2 + m    0.60 1.80E-31 -3.20  0.0 4.70E-12  0.00  0.0
!     IUPAC no change

at(1) = 0.60
at(2) = 1.80e-31
at(3) = -3.20
at(4) = 0.0
at(5) = 4.70e-12
at(6) = 0.00
at(7) = 0.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,90) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,90) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)    = zo(i) / zi(i)
    arg2(i)  = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,90) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

! R91 HO2NO2 + m = HO2 + NO2  0.60 4.10E-05  0.00  10650.0 4.80E+15 0.00 11170.0
!     IUPAC no change

at(1) = 0.60
at(2) = 4.10e-5
at(3) = 0.0
at(4) = 10650.0
at(5) = 4.80e+15
at(6) = 0.00
at(7) = 11170.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,91) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,91) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)    = zo(i) / zi(i)
    arg2(i)  = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,91) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

!     R92 MeCO3 + NO2 = PAN + m  0.30 2.70E-28 -7.10   0.0 1.20E-11 -0.90    0.0
!     IUPAC no change

at(1) = 0.30
at(2) = 2.70e-28
at(3) = -7.10
at(4) = 0.0
at(5) = 1.21e-11
at(6) = -0.90
at(7) = 0.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,92) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,92) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)    = zo(i) / zi(i)
    arg2(i)  = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,92) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

!   R93 PAN + m = MeCO3 + NO2  0.30 4.90E-03  0.00 12100.0 5.40E+16 0.00 13830.0
!     IUPAC no change

at(1) = 0.30
at(2) = 4.90e-03
at(3) = 0.00
at(4) = 12100.0
at(5) = 5.40e+16
at(6) = 0.00
at(7) = 13830.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,93) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,93) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)    = zo(i) / zi(i)
    arg2(i)  = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,93) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

!   R94 N2O5 + m = NO2 + NO3  0.35 1.30E-03 -3.50 11000.0 9.70E+14  0.10 11080.0
!     IUPAC no change - note 6

at(1) = 0.35
at(2) = 1.30e-03
at(3) = -3.50
at(4) = 11000.0
at(5) = 9.70e+14
at(6) = 0.10
at(7) = 11080.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,94) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,94) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)    = zo(i) / zi(i)
    arg2(i)  = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,94) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

!     R95 NO2 + NO3 = N2O5 + m   0.35 3.60E-30 -4.10   0.0 1.90E-12  0.20   0.0
!     IUPAC no change - note 7

at(1) = 0.35
at(2) = 3.60e-30
at(3) = -4.10
at(4) = 0.0
at(5) = 1.90e-12
at(6) = 0.20
at(7) = 0.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,95) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,95) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)    = zo(i) / zi(i)
    arg2(i)  = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,95) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

!     R96 O(3P) + O2 = O3 + m      0.00 5.70E-34 -2.60  0.0 0.00E+00  0.00   0.0
!     IUPAC no change - note 8

at(1) = 0.00
at(2) = 5.70e-34
at(3) = -2.60
at(4) = 0.0
at(5) = 0.00e+00
at(6) = 0.00
at(7) = 0.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,96) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,96) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)    = zo(i) / zi(i)
    arg2(i)  = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,96) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

! R97 OH + NO = HONO + m       1420.00 7.40E-31 -2.40   0.0 3.30E-11 -0.30  0.0
! IUPAC no change - note 9

at(1) = 1420.0
at(2) = 7.40e-31
at(3) = -2.40
at(4) = 0.0
at(5) = 3.30e-11
at(6) = -0.30
at(7) = 0.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,97) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,97) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)    = zo(i) / zi(i)
    arg2(i)  = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,97) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

!     R98 OH + NO2 = HONO2 + m   0.40 3.30E-30 -3.00   0.0 4.10E-11  0.00    0.0
!     IUPAC no change - note 10

at(1) = 0.40
at(2) = 3.30e-30
at(3) = -3.00
at(4) = 0.0
at(5) = 4.10e-11
at(6) = 0.00
at(7) = 0.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,98) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,98) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)    = zo(i) / zi(i)
    arg2(i)  = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,98) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

!     R99 OH + OH = H2O2 + m      0.50 6.90E-31 -0.80   0.0 2.60E-11  0.00   0.0
!     IUPAC no change

at(1) = 0.50
at(2) = 6.90e-31
at(3) = -0.80
at(4) = 0.0
at(5) = 2.60e-11
at(6) = 0.00
at(7) = 0.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,99) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,99) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)    = zo(i) / zi(i)
    arg2(i)  = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,99) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

!     R100 EtCO3 + NO2 = PPAN + m   0.30 2.70E-28 -7.10  0.0 1.20E-11 -0.90  0.0
!     MCM v3.1 no change - note 12

at(1) = 0.30
at(2) = 2.70e-28
at(3) = -7.10
at(4) = 0.0
at(5) = 1.20e-11
at(6) = -0.90
at(7) = 0.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,100) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,100) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)     = zo(i) / zi(i)
    arg2(i)   = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,100) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

! R101 PPAN + m = EtCO3 + NO2 0.36 1.70E-03 0.00  11280.0 8.30E+16  0.00 13940.0
!     IUPAC no change

at(1) = 0.36
at(2) = 1.70e-03
at(3) = 0.00
at(4) = 11280.0
at(5) = 8.30e+16
at(6) = 0.00
at(7) = 13940.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,101) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,101) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)     = zo(i) / zi(i)
    arg2(i)   = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,101) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

!     R102 O(3P)S + O2 = O3S + m   0.00 5.70E-34 -2.60   0.0 0.00E+00  0.00 0.0
!     IUPAC no change - note 8

at(1) = 0.00
at(2) = 5.70e-34
at(3) = -2.60
at(4) = 0.0
at(5) = 0.00e+00
at(6) = 0.00
at(7) = 0.0

zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

DO i=1,n_pnts
  IF ( zo(i) < peps ) THEN
    rc(i,102) = zi(i)
  ELSE IF ( zi(i) < peps ) THEN
    rc(i,102) = zo(i)
  ELSE
    IF ( at(1) <= 1.0 ) THEN
      zfc(i) = at(1)
    ELSE
      zfc(i) = EXP( -tc(i)/at(1) )
    END IF
    zr(i)     = zo(i) / zi(i)
    arg2(i)   = 1.0/(1.0+((LOG10(zr(i)))**2))
    rc(i,102) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
  END IF
END DO

IF (ukca_config%l_ukca_aerchem) THEN

  !       R103  SO2      OH     m    H2SO4   Rate data from NASA/JPL.

  at(1) = 0.60
  at(2) = 3.00e-31
  at(3) = -3.30
  at(4) = 0.0
  at(5) = 1.50e-12
  at(6) = 0.00
  at(7) = 0.0

  zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
  zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

  DO i=1,n_pnts
    IF (zo(i) < peps ) THEN
      rc(i,103) = zi(i)
    ELSE IF ( zi(i) < peps ) THEN
      rc(i,103) = zo(i)
    ELSE
      IF ( at(1) <= 1.0 ) THEN
        zfc(i) = at(1)
      ELSE
        zfc(i) = EXP( -tc(i)/at(1) )
      END IF
      zr(i) = zo(i) / zi(i)
      arg2(i) = 1.0/(1.0+((LOG10(zr(i)))**2))
      rc(i,103) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*LOG(zfc(i)))
      rc(i,103) = rc(i,103) * (1.0 - (faq(i,48) * fcloud(i)))
    END IF
  END DO

  !       R104:  HSO3(-)(aq) + H2O2(aq) =>  H2SO4
  !       Rate data from Eqn. 7 of Bower et al (1991, Atm.Env. 25A, 2401-2418)
  !       Adjustments to allow for use of aqueous chemistry in model with
  !       gas-phase reactions from Bergel et al (2004, JGR 109).

  WHERE (clw(:) < 1e-10)
    rc(:,104) = 0.0
  ELSE WHERE
    ! Replace Global pH with 1D values from H_plus array
    rc(:,104) = 1.10e+12 * (EXP(-3655.3 / tc(:))) *                            &
                (H_plus_1d_arr(:)/(0.1 + H_plus_1d_arr(:)))*fcloud(:)          &
                *fdiss(:,48,2)*fdiss(:,13,1) * 1000.0 / (avogadro * vr(:))
  END WHERE

  !       DMS Oxidation from IUPAC March 2007.
  !       R105 : DMS + OH = DMSO2

  rc(:,105) = 1.12e-11*EXP(-250.0/tc(:))

  !       R106 : DMS + OH = DMSO2

  rc(:,106) = 9.3e-39*o2(:)*EXP(5270.0/tc(:))/                                 &
            (1.0 + 7.4e-29*o2(:)*EXP(5610.0/tc(:)))

  !       R107 : DMS + NO3 = DMSO2 + HNO3

  rc(:,107) = 1.9e-13*EXP(520.0/tc(:))

  !       R108, R109, R110 : parameterised DMS

  !       R111 : NH3 + OH = NH2 + H2O  (NH2 + NO = N2 + H2O)
  rc(:,111) = 3.5e-12*EXP(-925.0/tc(:))

  !       R112 : MonoTerp + OH = Sec_Org
  rc(:,112) = 1.2e-11*EXP(444.0/tc(:))            ! MM (2008)

  !       R113 : MonoTerp + O3 = Sec_Org
  rc(:,113) = 1.01e-15*EXP(-732.0/tc(:))          ! MM (2008)

  !       R114 : MonoTerp + NO3 = Sec_Org
  rc(:,114) = 1.19e-12*EXP(490.0/tc(:))           ! MM (2008)

  !       R115 : SO3 + O3 (Aq) =  SO4                     !
  WHERE (clw(:) < 1e-10)
    rc(:,115) = 0.0
  ELSE WHERE
    !          rc(:,115) = ((4.28E+13*(EXP(-5532.8/tc(:)))*fdiss(:,48,2))+   &
    ! Only so3 + o3 for now.
    rc(:,115) = (7.43e+16*(EXP(-5280.0/tc(:)))*fdiss(:,48,3))*                 &
               fcloud(:)*fdiss(:,3,1)*1000.0/(avogadro*vr(:))
  END WHERE

  !       These are not used directly, but to calculate product ratios

  !       CH3SO2 => CH3 + SO2
  k_dms(:,1) = 100.0                       ! Karl et al 2007

  !       CH3SO2 + O3 => CH3SO3
  k_dms(:,2) = 6.3e-13                     ! Karl et al 2007

  !       CH3SO2 + NO2 => CH3SO3
  k_dms(:,3) = 2.2e-11                     ! Karl et al 2007

  !       CH3SO3 + HO2 => MSA
  k_dms(:,4) = 5.0e-11                     ! Karl et al 2007

  !       CH3SO3 => CH3 + SO3
  k_dms(:,5) = 1.2e-3                      ! Karl et al 2007

END IF    ! l_ukca_aerchem


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chemco
END MODULE ukca_chemco_mod
