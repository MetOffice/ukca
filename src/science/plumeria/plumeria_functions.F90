! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module containing the functions used in the Plumeria model.

!  Methods:
!    This module contains the below functions:
!    plumeria_psat, plumeria_Airtemp, plumeria_AirPres,
!    plumeria_AirHumid, plumeria_zeta, plumeria_wind, plumeria_vx,
!    plumeria_vy, plumeria_h_a, plumeria_cp_a, plumeria_cp_l,
!    plumeria_h_i, plumeria_h_l, plumeria_h_m, plumeria_h_v,
!    plumeria_Tsat
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE plumeria_functions_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PLUMERIA_FUNCTIONS_MOD'

CONTAINS

!=======================================================================!

FUNCTION plumeria_psat(t)
! ----------------------------------------------------------------------
! Function that returns the partial pressure of water (in Pascals) at
! saturation, given a temperature (t) in Kelvin. Taken from  Ghiorso.
! ----------------------------------------------------------------------

USE plumeria_param_mod, ONLY: T_ice, T_ColdWater

IMPLICIT NONE

! Input variable
REAL, INTENT(IN) :: t   ! Temperature

! Function return value
REAL :: plumeria_psat   ! Partial pressure of water at saturation

! Local variables
REAL :: a(8)
REAL :: v, ff, e0, MassFrac, w, Zhaar
INTEGER  :: i

DATA a/-7.8889166,2.5514255,-6.716169,33.239495,-105.38479,174.35319,          &
       -148.39348,48.631602/

e0 = 611.0 !partial pressure of water vapor at 273.15 K

IF (t < T_ice) THEN
  plumeria_psat = EXP(LOG(e0) + 22.49 - 6142 / t)
  RETURN

ELSE IF (t < T_ColdWater) THEN
  !in the temperature range (T_ice<t<=T_coldwater) where liquid & ice coexist,
  !the partial pressure is weighted between the partial pressure of
  !ice and of liquid.
  MassFrac = (t - T_ice) / (T_ColdWater - T_ice)
  plumeria_psat = (1.0 - MassFrac) * EXP( LOG(e0) + 22.49 - 6142 / t ) +       &
                   MassFrac * ( 100000.0 * EXP( 6.3573118 - (8858.843 / t) +   &
                   ( 607.56335 / (t ** 0.6) ) ) )
  RETURN

ELSE IF (t < 314.0) THEN
  plumeria_psat = 100000.0 * EXP(6.3573118 - (8858.843 / t) +                  &
                  ( 607.56335 / (t ** 0.6) ) )
  RETURN

ELSE
  v = t / 647.25
  w = ABS(1.0 - v * 0.9997)
  ff = 0.0

  DO i = 1,8
    Zhaar = i
    ff = ff + a(i) * w ** ((Zhaar + 1.0) / 2.0)
  END DO
  plumeria_psat = 100000.0 * 1.001 * 220.93 * EXP(ff / v)
END IF

END FUNCTION plumeria_psat

!=======================================================================!

FUNCTION plumeria_AirTemp(znow)
! ----------------------------------------------------------------------
! Function that calculates air temperature in a layered atmosphere,
! given the elevation znow above the vent
! ----------------------------------------------------------------------

USE plumeria_param_mod, ONLY: vent_elevation, ZairLayer, iatmlayers,           &
                              TairLayer

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: znow    ! Elevation above vent

! Function return value
REAL :: plumeria_AirTemp    ! Air temperature

! Local variable
INTEGER  :: ilayer

IF (((znow + vent_elevation) > ZairLayer(1)) .AND.                             &
     ((znow + vent_elevation) < ZairLayer(iatmlayers))) THEN
  check_airtemp: DO ilayer = 2, iatmlayers
    IF ((znow + vent_elevation) < ZairLayer(ilayer)) THEN
      plumeria_AirTemp = TairLayer(ilayer - 1) +                               &
                         ( (znow + vent_elevation) - ZairLayer(ilayer - 1) ) * &
                         ( TairLayer(ilayer) - TairLayer(ilayer - 1) ) /       &
                         ( ZairLayer(ilayer) - ZairLayer(ilayer - 1) )
      EXIT check_airtemp
    END IF
  END DO check_airtemp
  !If we!re below the lowest sounding
ELSE IF ((znow + vent_elevation) <= ZairLayer(1)) THEN
  plumeria_AirTemp = TairLayer(1) +                                            &
                     ( (znow + vent_elevation) - ZairLayer(1)) *               &
                     ( TairLayer(2) - TairLayer(1) ) /                         &
                     ( ZairLayer(2) - ZairLayer(1) )
  !If we!re above the highest sounding
ELSE
  plumeria_AirTemp = TairLayer(iatmlayers) +                                   &
                     ( (znow + vent_elevation) - ZairLayer(iatmlayers - 1) ) * &
                     ( TairLayer(iatmlayers) - TairLayer(iatmlayers - 1) ) /   &
                     ( ZairLayer(iatmlayers) - ZairLayer(iatmlayers - 1) )
END IF
RETURN

END FUNCTION plumeria_AirTemp

!=======================================================================!

FUNCTION plumeria_AirPres(VentElevation)
! ----------------------------------------------------------------------
! Function is called only to find air pressure at the vent
! ----------------------------------------------------------------------

USE plumeria_param_mod,    ONLY: ZairLayer, TairLayer, iatmlayers, pAirLayer
USE comorph_constants_mod, ONLY: R_dry, gravity

IMPLICIT NONE

! Argument
REAL, INTENT(IN):: VentElevation   ! Vent altitude

! Function return value
REAL :: plumeria_AirPres           ! Air pressure at the vent

! Local variables
INTEGER   :: ilayer
REAL      :: Tnow, p0, LapseRate, t0, z0

Tnow = plumeria_AirTemp(0.0) !gives air temperature at the vent elevation

IF ((VentElevation >= ZairLayer(1)) .AND.                                      &
    (VentElevation < ZairLayer(iatmlayers))) THEN
  check_airpres: DO ilayer = 2, iatmlayers
    IF (VentElevation < ZairLayer(ilayer)) THEN
      LapseRate = (TairLayer(ilayer) - TairLayer(ilayer - 1)) /                &
                  (ZairLayer(ilayer) - ZairLayer(ilayer - 1))
      p0 = pAirLayer(ilayer - 1)
      z0 = ZairLayer(ilayer - 1)
      t0 = plumeria_AirTemp(ZairLayer(ilayer - 1) - VentElevation)
      EXIT check_airpres
    END IF
  END DO check_airpres
  !If we!re below the lowest sounding
ELSE IF (VentElevation < ZairLayer(1)) THEN
  LapseRate = ( TairLayer(2) - TairLayer(1) ) /                                &
              ( ZairLayer(2) - ZairLayer(1) )
  p0 = pAirLayer(1)
  z0 = ZairLayer(1)
  t0 = plumeria_AirTemp(ZairLayer(1) - VentElevation)
  !If we!re above the highest sounding
ELSE
  LapseRate = ( TairLayer(iatmlayers) - TairLayer(iatmlayers - 1) ) /          &
              ( ZairLayer(iatmlayers) - ZairLayer(iatmlayers - 1) )
  p0 = pAirLayer(iatmlayers)
  z0 = ZairLayer(iatmlayers)
  t0 = plumeria_AirTemp(ZairLayer(iatmlayers) - VentElevation)
END IF

IF (ABS(LapseRate) > 1.0e-08) THEN
  plumeria_AirPres = p0 * (t0 / Tnow) ** (gravity / (R_dry * LapseRate))
ELSE
  plumeria_AirPres = p0 * EXP( -gravity *                                      &
                      (VentElevation - z0) / (R_dry * Tnow) )
END IF

END FUNCTION plumeria_AirPres

!=======================================================================!

FUNCTION plumeria_AirHumid(znow)
! ----------------------------------------------------------------------
! Function that calculates humidity in a layered atmosphere,
! given the elevation z above the vent
! ----------------------------------------------------------------------

USE plumeria_param_mod, ONLY: vent_elevation, ZairLayer, iatmlayers,           &
                              HumidAirLayer

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: znow      ! Elevation above the vent

! Function return value
REAL :: plumeria_AirHumid     ! Humidity in the atmosphere

! Local variable
INTEGER       :: ilayer

IF (((znow + vent_elevation) > ZairLayer(1)) .AND.                             &
  ((znow + vent_elevation) < ZairLayer(iatmlayers))) THEN
  check_airhumid: DO ilayer = 2, iatmlayers
    IF ((znow + vent_elevation) < ZairLayer(ilayer)) THEN
      plumeria_AirHumid = HumidAirLayer(ilayer - 1) +                          &
                         ( (znow + vent_elevation) - ZairLayer(ilayer - 1) ) * &
                         ( HumidAirLayer(ilayer) - HumidAirLayer(ilayer - 1)) /&
                         ( ZairLayer(ilayer) - ZairLayer(ilayer - 1) )
      EXIT check_airhumid
    END IF
  END DO check_airhumid
  !If we!re below the lowest sounding
ELSE IF ((znow + vent_elevation) <= ZairLayer(1)) THEN
  plumeria_AirHumid = HumidAirLayer(1) +                                       &
                      ( (znow + vent_elevation) - ZairLayer(1) ) *             &
                      ( HumidAirLayer(2) - HumidAirLayer(1) ) /                &
                      ( ZairLayer(2) - ZairLayer(1) )
  !If we!re above the highest sounding
ELSE
  plumeria_AirHumid = HumidAirLayer(iatmlayers) +                              &
             ((znow + vent_elevation) - ZairLayer(iatmlayers - 1)) *           &
             (HumidAirLayer(iatmlayers) - HumidAirLayer(iatmlayers - 1)) /     &
             (ZairLayer(iatmlayers) - ZairLayer(iatmlayers - 1))
END IF

END FUNCTION plumeria_AirHumid

!=======================================================================!

FUNCTION plumeria_zeta(znow)
! ----------------------------------------------------------------------
! Function that calculates the wind direction parallel to the plume
! ----------------------------------------------------------------------

USE plumeria_param_mod, ONLY: vent_elevation, ZairLayer, iatmlayers,           &
                              WindspeedLayer, WinddirLayer
USE conversions_mod, ONLY: pi_over_180

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: znow      ! Elevation above the vent

! Function return value
REAL :: plumeria_zeta         ! Wind directoin parallel to the plume

! Local variables
INTEGER :: ilayer
REAL    :: xwind1, xwind2, xwindnow, ywind1, ywind2, ywindnow

ilayer = 0

IF (((znow + vent_elevation) > ZairLayer(1)) .AND.                             &
   ((znow + vent_elevation) < ZairLayer(iatmlayers))) THEN
  check_zeta: DO ilayer = 2, iatmlayers
    IF ((znow + vent_elevation) < ZairLayer(ilayer)) THEN
      xwind1 = WindspeedLayer(ilayer - 1) *                                    &
               SIN(pi_over_180 * (WinddirLayer(ilayer - 1) + 180.0))
      xwind2 = WindspeedLayer(ilayer) *                                        &
               SIN(pi_over_180 * (WinddirLayer(ilayer) + 180.0))
      ywind1 = WindspeedLayer(ilayer - 1) *                                    &
               COS(pi_over_180 * (WinddirLayer(ilayer - 1) + 180.0))
      ywind2 = WindspeedLayer(ilayer) *                                        &
               COS(pi_over_180 * (WinddirLayer(ilayer) + 180.0))
      xwindnow = xwind1 +                                                      &
                ( (znow + vent_elevation) - ZairLayer(ilayer - 1) ) *          &
                ( xwind2 - xwind1 ) /                                          &
                ( ZairLayer(ilayer) - ZairLayer(ilayer - 1) )
      ywindnow = ywind1 +                                                      &
                ( (znow + vent_elevation) - ZairLayer(ilayer-1) ) *            &
                ( ywind2 - ywind1 ) /                                          &
                ( ZairLayer(ilayer) - ZairLayer(ilayer-1) )
      plumeria_zeta = ATAN2(xwindnow, ywindnow)
      EXIT check_zeta
    END IF
  END DO check_zeta
  !If we!re below the lowest sounding
ELSE IF ((znow + vent_elevation) <= ZairLayer(1)) THEN
  plumeria_zeta = pi_over_180 * (WinddirLayer(1) + 180.0)
  !If we!re above the highest sounding
ELSE
  plumeria_zeta = pi_over_180 * (WinddirLayer(iatmlayers) + 180.0)
END IF

END FUNCTION plumeria_zeta

!=======================================================================!

FUNCTION plumeria_wind(znow)
! ----------------------------------------------------------------------
! Function that calculates the wind speed
! ----------------------------------------------------------------------

USE plumeria_param_mod, ONLY: vent_elevation, ZairLayer, iatmlayers,           &
                              WindspeedLayer, WinddirLayer
USE conversions_mod, ONLY: pi_over_180

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: znow     ! Elevation above vent

! Function return value
REAL :: plumeria_wind        ! Wind speed

! Local variables
INTEGER        :: ilayer
REAL           :: xwind1, xwind2, xwindnow, ywind1, ywind2, ywindnow

ilayer = 0

IF (((znow + vent_elevation) > ZairLayer(1)) .AND.                             &
  ((znow + vent_elevation) < ZairLayer(iatmlayers))) THEN
  check_wind: DO ilayer = 2, iatmlayers
    IF ((znow + vent_elevation) < ZairLayer(ilayer)) THEN
      xwind1 = WindspeedLayer(ilayer - 1) *                                    &
               SIN(pi_over_180 * (WinddirLayer(ilayer - 1) + 180.0))
      xwind2 = WindspeedLayer(ilayer) *                                        &
               SIN(pi_over_180 * (WinddirLayer(ilayer) + 180.0))
      ywind1 = WindspeedLayer(ilayer - 1) *                                    &
               COS(pi_over_180 * (WinddirLayer(ilayer - 1) + 180.0))
      ywind2 = WindspeedLayer(ilayer) *                                        &
               COS(pi_over_180 * (WinddirLayer(ilayer) + 180.0))
      xwindnow = xwind1 + ((znow + vent_elevation) - ZairLayer(ilayer - 1)) *  &
                ( xwind2 - xwind1 ) /                                          &
                ( ZairLayer(ilayer) - ZairLayer(ilayer - 1))
      ywindnow = ywind1 + ((znow + vent_elevation) - ZairLayer(ilayer - 1)) *  &
                (ywind2 - ywind1) /                                            &
                (ZairLayer(ilayer) - ZairLayer(ilayer - 1))
      plumeria_wind = SQRT(xwindnow ** 2 + ywindnow ** 2)
      EXIT check_wind
    END IF
  END DO check_wind
  !If we!re below the lowest sounding
ELSE IF ((znow + vent_elevation) <= ZairLayer(1)) THEN
  plumeria_wind = WindspeedLayer(1)
  !If we!re above the highest sounding
ELSE
  plumeria_wind = WindspeedLayer(iatmlayers)
END IF

END FUNCTION plumeria_wind

!=======================================================================!

FUNCTION plumeria_vx(znow)
! ----------------------------------------------------------------------
! Function that calculates the x component of the wind vector
! ----------------------------------------------------------------------

USE conversions_mod, ONLY: pi

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: znow     ! Elevation above vent

! Function return value
REAL :: plumeria_vx          ! winter vector in x direction

plumeria_vx = plumeria_wind(znow) * COS(pi / 2.0 - plumeria_zeta(znow))

END FUNCTION plumeria_vx

!=======================================================================!

FUNCTION plumeria_vy(znow)
! ----------------------------------------------------------------------
! Function that calculates the y component of the wind vector
! ----------------------------------------------------------------------

USE conversions_mod, ONLY: pi

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: znow      ! Elevation above vent

! Function return value
REAL :: plumeria_vy           ! wind vector in y direction

plumeria_vy = plumeria_wind(znow) * SIN(pi / 2.0 - plumeria_zeta(znow))

END FUNCTION plumeria_vy

!=======================================================================!

FUNCTION plumeria_h_a(t_k)
! ----------------------------------------------------------------------
! Function that calculates the enthalpy of air using polynomial
! coefficients for cp_air on p. 718 of Moran & Shapiro (1992).
! ----------------------------------------------------------------------

USE comorph_constants_mod, ONLY: R_dry

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: t_k    ! Temperature

! Function return value
REAL :: plumeria_h_a       ! Enthalpy of air

! Local variables
REAL :: h0, t0, alpha, beta, gamma_coef, delta, epsilon_coef

h0 = 270110.0     !enthalpy at T=270 K (from Table A-15 of Moran & Shapiro, 1992)
t0 = 270.0
alpha = 3.653
beta = -0.001337
gamma_coef = 0.000003294
delta = -0.000000001913
epsilon_coef = 2.763e-13

plumeria_h_a = h0 + R_dry * ( ( alpha * (t_k - t0) ) +                         &
                              ( beta / 2) * (t_k ** 2 - t0 ** 2) +             &
                              ( gamma_coef / 3) * (t_k ** 3 - t0 ** 3) +       &
                              ( delta / 4) * (t_k ** 4 - t0 ** 4) +            &
                              ( epsilon_coef / 5) * (t_k ** 5 - t0 ** 5) )

END FUNCTION plumeria_h_a

!=======================================================================!

FUNCTION plumeria_cp_a(t_k)
! ----------------------------------------------------------------------
! Function that calculates the cp of air using polynomial coefficients
! on p. 718 of Moran & Shapiro (1992).
! Based on NASA SP-273, U.S. Government Printing Office, Washington, D.C. 1971.
! Moran, M.J., and Shapiro, H.N., Fundamentals of Engineering Thermodynamnics,
! 2nd Ed., John Wiley & Sons, New York, 804pp.
! ----------------------------------------------------------------------

USE comorph_constants_mod, ONLY: R_dry

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: t_k   ! Temperature

! Function return value
REAL :: plumeria_cp_a     ! Specific heat of air

! Local variables
REAL :: alpha, beta, gamma_coef, delta, epsilon_coef

!polynomial fitting coefficients
alpha = 3.653
beta = -0.001337
gamma_coef = 0.000003294
delta = -0.000000001913
epsilon_coef = 2.763e-13

plumeria_cp_a = R_dry * ( alpha +                                              &
                          beta * t_k +                                         &
                          gamma_coef * t_k ** 2 +                              &
                          delta * t_k ** 3 +                                   &
                          epsilon_coef * t_k ** 4 )

END FUNCTION plumeria_cp_a

!=======================================================================!

FUNCTION plumeria_cp_l(t_k)
! ----------------------------------------------------------------------
! Function that calculates the cp of liquid water.
! For a given temperature, the enthalpy is calculated by taking the enthalpy
! value at the next lower temperature in the table (at 10-degree increments,
! from 0.1 C to 99 C at saturation pressure), and then use the formula:
!                     h = h(last T) + cp*delta(T)
! This gives values that are 0.3 to 0.0005% from those in Haar et al.
! with the lowest accuracies at lower temperatures.
!
! Reference:
! Haar, L., Gallagher, J.S., and Kell, G.S., 1984, NBS/NRC Steam Tables,
! Hemisphere Publishing Corporation, New York, 320 pp.
! ----------------------------------------------------------------------

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: t_k     ! Temperature

! Function return value
REAL :: plumeria_cp_l       ! Specific heat of liquid water

! Local variables
REAL :: cp(11), tk(11)
INTEGER       :: i

INTEGER :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage   ! Error message

CHARACTER(LEN=*), PARAMETER :: RoutineName='PLUMERIA_CP_L'

DATA cp/4227.9,4188.0,4183.3,4183.3,4182.4,4181.7,4182.9,4187.0,4194.3,        &
     4204.5,4217.2/
DATA tk/273.25,283.25,293.25,303.25,313.25,323.25,333.25,343.25,353.25,        &
     363.25,373.25/

IF (t_k < tk(1)) THEN
  plumeria_cp_l = cp(1)                   !needed for some cases
  cmessage = "Warning: temperature below 0.1 C"
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)

ELSE IF (t_k > tk(11)) THEN
  plumeria_cp_l = cp(11)
  cmessage = "Warning: liquid water temperature above 100 C"
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)

ELSE
  check_cp_l: DO i=1,11
    IF (tk(i) > t_k) THEN
      plumeria_cp_l = cp(i - 1) +                                              &
                      ( t_k - tk(i - 1) ) * ( cp(i - 1) + cp(i) ) /            &
                      ( tk(i) - tk(i - 1) )
      EXIT check_cp_l
    END IF
  END DO check_cp_l
END IF

END FUNCTION plumeria_cp_l

!=======================================================================!

FUNCTION plumeria_h_i(t_k)
! ----------------------------------------------------------------------
! Function that gives enthalpy of ice as a function of temperature.
! It assumes that the enthalpy is independent of pressure and that the
! enthalpy is simply the enthalpy of ice at T=0 C, p=1 atm (-333430), plus
! cp of ice (~1850 J/kg K) times the temperature in Celsius
! ----------------------------------------------------------------------

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: t_k       ! Temperature
! Function return value
REAL :: plumeria_h_i          ! Enthalpy of ice

plumeria_h_i = -333430 + 1850.0 * (t_k - 273.15)

END FUNCTION plumeria_h_i

!=======================================================================!

FUNCTION plumeria_h_l(t_k)
! ----------------------------------------------------------------------
! Function that calculates the enthalpy of liquid water.
! For a given temperature, the enthalpy is calculated by taking the enthalpy
! value at the next lower temperature in the table (at 10-degree increments,
! from 0.1 C to 99 C at saturation pressure), and then use the formula:
!                     h = h(last T) + cp*delta(T)
! This gives values that are 0.3 to 0.0005% from those in Haar et al.
! with the lowest accuracies at lower temperatures.
!
! Reference:
! Haar, L., Gallagher, J.S., and Kell, G.S., 1984, NBS/NRC Steam Tables,
! Hemisphere Publishing Corporation, New York, 320 pp.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: t_k       ! Temperature

! Function return value
REAL :: plumeria_h_l          ! Enthalpy of liquid water

! Local variables
REAL :: cp(11), tk(11), hnow(11)
INTEGER  :: i

DATA cp/4227.9,4188.0,4183.3,4183.3,4182.4,4181.7,4182.9,4187.0,4194.3,        &
     4204.5,4217.2/
DATA tk/273.25,283.25,293.25,303.25,313.25,323.25,333.25,343.25,353.25,        &
     363.25,373.25/
DATA hnow/381.14,42406.0,84254.0,126090.0,167920.0,209750.0,251570.0,293430.0, &
     335350.0,377350.0,419490.0/

IF (tk(1) >= t_k) THEN
  plumeria_h_l = hnow(1) + 4180.0 * (t_k - tk(1))
ELSE IF (tk(11) <= t_k) THEN
  plumeria_h_l = hnow(11) + 4200.0 * (t_k - tk(11))
ELSE
  check_h_l: DO i = 1,11
    IF (tk(i) > t_k) THEN
      plumeria_h_l = hnow(i - 1) + (t_k - tk(i - 1)) * (cp(i - 1) + cp(i)) / 2
      EXIT check_h_l
    END IF
  END DO check_h_l
END IF

END FUNCTION plumeria_h_l

!=======================================================================!

FUNCTION plumeria_h_m(t_k, p_Pa)
! ----------------------------------------------------------------------
! Function that calculates magma enthalpy, assuming constant specific heat
! and density
! ----------------------------------------------------------------------

USE plumeria_param_mod, ONLY: Cp_m, rho_m

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: t_k          ! Temperature
REAL, INTENT(IN) :: p_Pa         ! Partial pressure

! Function return value
REAL :: plumeria_h_m             ! Enthalpy of magma

! Local variables
REAL :: t0, p0

t0 = 273.15  !reference temperature
p0 = 101300.0 !reference pressure

plumeria_h_m = Cp_m * (t_k - t0) + ((p_Pa - p0) / rho_m)

END FUNCTION plumeria_h_m

!=======================================================================!

FUNCTION plumeria_h_v(t_k)
! ----------------------------------------------------------------------
! Function that calculates the enthalpy of water vapor
! ----------------------------------------------------------------------

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: t_k       ! Temperature

! Function return value
REAL :: plumeria_h_v          ! Enthalpy of water vapor

!Polynomial coefficients for specific heat of water vapor
REAL :: alpha, beta, gamma_coef, delta, epsilon_coef

! Local variables
REAL :: t0, h0, cp

!Values taken from Moran & Shapiro, "Engineering Thermodynamics", 2nd. Ed.,
!John Wiley & Sons, 1984, p. 718.
alpha = 4.07
beta = -0.001108
gamma_coef = 0.000004152
delta = -0.000000002964
epsilon_coef = 0.000000000000807

h0 = 2500700.0                !enthalpy at T=0.1 C, p=0.1013 MPa
t0 = 273.15                  !reference temperature (K)

plumeria_h_v = h0 + (8.314 / 0.0180152) *                                      &
               ( (alpha * (t_k - t0) ) +                                       &
               ( (beta / 2.0) * (t_k ** 2 - t0 ** 2) ) +                       &
               ( (gamma_coef / 3.0) * (t_k ** 3 - t0 ** 3) ) +                 &
               ( (delta / 4.0) * (t_k ** 4 - t0 ** 4) ) +                      &
               ( (epsilon_coef / 5.0) * (t_k ** 5 - t0 ** 5)) )

END FUNCTION plumeria_h_v

!=======================================================================!

FUNCTION plumeria_Tsat(p_Pa)
! ----------------------------------------------------------------------
! Function that calculates saturation temperature of water at the
! given pressure
! ----------------------------------------------------------------------

USE comorph_constants_mod, ONLY: R_vap

IMPLICIT NONE

! Argument
REAL, INTENT(IN) :: p_Pa        ! partial pressure

! Function return value
REAL :: plumeria_Tsat           ! Saturation temperature of water

! Local variables
REAL :: T_boiling, pnow, t0

pnow = 100000.0    !default value
t0 = 273.15        !reference temperature (K)

!first guess, using best-fit curve
T_boiling = 182.82 * p_Pa ** 0.0611
pnow = plumeria_psat(T_boiling)

!Adjust until we get within 10 Pascals
check_pnow: DO
  IF (ABS(pnow - p_Pa) >= 10.0) THEN

    IF (T_boiling < t0) THEN
      !This statement keeps the function h_l from blowing up when it!s given a
      ! temperature value less than freezing.
      T_boiling = T_boiling +                                                  &
                 ( p_Pa - plumeria_psat(T_boiling) ) * R_vap * T_boiling ** 2 /&
                 ( ( plumeria_h_v(t0) - plumeria_h_l(t0) ) * p_Pa )
    ELSE
      T_boiling = T_boiling +                                                  &
                ( p_Pa - plumeria_psat(T_boiling) ) * R_vap * T_boiling ** 2 / &
                ( ( plumeria_h_v(T_boiling) - plumeria_h_l(T_boiling)) * p_Pa )
    END IF
    pnow = plumeria_psat(T_boiling)

  ELSE
    EXIT check_pnow
  END IF
END DO check_pnow

plumeria_Tsat = T_boiling

END FUNCTION plumeria_Tsat

END MODULE plumeria_functions_mod
