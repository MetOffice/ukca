! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module containing the main subroutine for running the Plumeria model
!    during explosive volcanic eruption (called from UKCA_VOLCANIC_SO2) to
!    calculate the neutral buoyancy level of volcanic plume.
!
!  Methods:
!    1. Prescribed eruption properties (i.e., vent altitude, mass eruption
!       rate) and the vertical column of atmospheric variables at the
!       location of the volcano are passed as arguments to 'plume_main'
!       with INTENT IN OUT.
!
!    2. Plumeria subroutines are called for the calculation of the
!       neutral buoyancy level of the volcanic plume
!       - plumeria_cashkarp
!       - plumeria_functions
!       - plumeria_derivs
!       - plumeria_FindT
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

MODULE plumeria_main_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PLUMERIA_MAIN_MOD'

CONTAINS

SUBROUTINE plumeria_main(model_levels, wind_speed, wind_dir, t_lev, p_lev,     &
                         geopH_lev, rh_lev, vent_alt, mer_erup, z_hb)

USE plumeria_param_mod,         ONLY: TairLayer, ZairLayer, HumidAirLayer,     &
                                  pAirLayer, WinddirLayer, WindspeedLayer,     &
                                  mdot, vent_elevation, iatmlayers,            &
                                  stopatthetop, m_gas, mw, n_0, n_0air,        &
                                  rho_0, T_mag, T_water, hmix, c_mix,          &
                                  zstep, hnext, hdid, dydx, yarr,              &
                                  m_a, m_i, m_l, m_m, m_w, m_v, p, phi,        &
                                  r, rho_air, rho_mix, s, T_air, T_mix,        &
                                  TIME, u, ux, uy, uz, x, y, z, rho_diff,      &
                                  alpha, gamma_w, n_exp, kgmole_air, kgmole_w, &
                                  Cp_m, rho_m, rho_w, rho_ice,                 &
                                  T_ColdWater, T_ice, gamma_gas, idensflag
USE plumeria_derivs_mod,        ONLY: plumeria_derivs
USE plumeria_cashkarp_mod,      ONLY: plumeria_cashkarpqs
USE plumeria_FindT_mod,         ONLY: plumeria_FindT_init, plumeria_FindT
USE plumeria_functions_mod,     ONLY: plumeria_AirTemp, plumeria_AirPres,      &
                                      plumeria_h_a, plumeria_h_l,              &
                                      plumeria_h_m, plumeria_h_v

USE comorph_constants_mod,  ONLY: R_dry, R_vap, gravity
USE conversions_mod,        ONLY: pi
USE umPrintMgr,             ONLY: umPrint, umMessage
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE, INTRINSIC :: ieee_arithmetic

USE parkind1,             ONLY: jpim, jprb      ! DrHook
USE yomhook,              ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! ----------------------------------------------------------------------
! Declarations
! ----------------------------------------------------------------------
! Arguments
INTEGER, INTENT(IN) :: model_levels
REAL, INTENT(IN) :: wind_speed(model_levels)   ! wind speed, in m/s
REAL, INTENT(IN) :: wind_dir(model_levels)     ! wind direction
REAL, INTENT(IN) :: t_lev(model_levels)        ! temperature, in K
REAL, INTENT(IN) :: p_lev(model_levels)        ! pressure level, in Pa
REAL, INTENT(IN) :: rh_lev(model_levels)       ! relative humidity,in fraction
REAL, INTENT(IN) :: geopH_lev(model_levels)    ! geopotential height, in m
REAL, INTENT(IN) :: vent_alt                   ! vent altitude, in m
REAL, INTENT(IN) :: mer_erup                   ! mass eruption rate, in kg/S
REAL, INTENT(OUT) :: z_hb                      ! Neutral Buoyancy Height,in km

! Error messages
INTEGER :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage

! Local variables
REAL :: ratio_ur        ! richardson number
REAL :: gperm3_l(10000) ! density g/m3 of liquid water
REAL :: gperm3_i(10000) ! density g/m3 of ice
REAL :: pathmax         ! Value of s/s_zmax at which calculations stop
REAL :: zmax            ! Maximum plume height
REAL :: s_zmax          ! Value of s at zmax
REAL :: uzchgmin        ! Minimum change in uz between successive steps
REAL :: hdid_min        ! Execution stops when hdid<hdid_min*r(1)
REAL :: umin            ! Velocity below which the simulation stops
! Other local variables
REAL :: snow, snow_out, du, sstep, uzlast, yscale(13)

! Logical to check whether plume is still rising
LOGICAL       :: StillRising

! Integer variables:
INTEGER ::  maxsteps         ! maximum number of integration steps
INTEGER ::  imax             ! Index at top of plume
INTEGER ::  istep            ! model timestep
INTEGER ::  hb_istep         ! timestep of inversion level
INTEGER ::  count_inversion  ! number of inversion levels

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PLUMERIA_MAIN'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! ----------------------------------------------------------------------
! Set default values
! ----------------------------------------------------------------------
kgmole_air = 8.314/ R_dry       ! Specific gas constant for air
kgmole_w   = 0.0180152          ! Molar weight of water
iatmlayers = 1                  ! Number of layers in atmosphere
imax       = 0                  ! Index at top of plume
pathmax    = 2.0                ! Value of s/s_zmax at which calculations stop
uzchgmin   = 1.0e-07            ! Minimum change in uz between successive steps
hdid_min   = 0.000001           ! Execution stops when hdid<hdid_min*r(1)
umin       = 0.01               ! Velocity below which the simulation stops
s_zmax     = 99999.0            ! Value of s at zmax
zmax       = 0.0                ! Maximum plume height
stopatthetop = .TRUE.           ! =.true if we want to stop at max elevation
StillRising  = .TRUE.           ! =.true if the plume has not yet reached
                                ! its peak height

! ----------------------------------------------------------------------
! Input arguements to local variable names
! ----------------------------------------------------------------------
pAirLayer = p_lev
ZairLayer = geopH_lev
TairLayer = t_lev
HumidAirLayer = rh_lev
WinddirLayer = wind_dir
WindspeedLayer = wind_speed
vent_elevation = vent_alt
mdot = mer_erup
iatmlayers = model_levels

! ----------------------------------------------------------------------
! Vent properties
! ----------------------------------------------------------------------
T_air(1) = plumeria_AirTemp(0.0)
mw = 0                                  ! Mass fraction added water
z(1) = 0.0                               ! height above the vent
s(1) = 0.0                               ! Distance along plume axis
x(1) = 0.0
y(1) = 0.0
TIME(1) = 0.0                            ! Time (in s) after leaving the vent

! ----------------------------------------------------------------------
! Magma properties
! ----------------------------------------------------------------------
T_mag = 1100                            ! Magma temperature (Celsius)
T_mag = T_mag + 273.15                  ! Convert from C to K
n_0 = 0.05                              ! Mass fraction gas in magma
n_0air = 0.001                          ! Mass fraction air
Cp_m = 1250                             ! Magma specific heat, J/kg K
rho_m = 2350                            ! Magma density (DRE), kg/m3
gamma_gas = 1.25                        ! cp/cv for magmatic gas

! ----------------------------------------------------------------------
! Troposphere properties
! ----------------------------------------------------------------------
p(1) = plumeria_AirPres(vent_elevation)    !air pressure at the vent
rho_air(1) = p(1) / (R_dry * T_air(1))         !air density at the vent

! ----------------------------------------------------------------------
! Water properties
! ----------------------------------------------------------------------
rho_w = 1000.0            ! Density of water
T_water = 273.15          ! Temperature of external water mixed
                          ! with magma at beginning
rho_ice = 900.0           ! Ice density
T_ColdWater = 266.65      ! Top of temperature range at which
                          ! liquid water & ice coexist
T_ice = 258.15            ! Bottom of temperature range which
                          ! liquid water & ice coexist

! ----------------------------------------------------------------------
! Set integration parameters
! ----------------------------------------------------------------------
maxsteps = 10000          ! Program will stop if maxsteps is exceeded
alpha = 0.1               ! Entrainment coefficient
gamma_w = 0.25            ! Crossflow entrainmnent coefficient
n_exp   = 1.0             ! Devenish exponent to entrainment equation

! ----------------------------------------------------------------------
! Mass fractions
! ----------------------------------------------------------------------
m_gas = (1.0 - mw) * n_0                ! mass fraction magmatic gas
m_w(1) = (1.0 - n_0air) * (mw + m_gas)  ! mass fraction total water in mixture
m_m(1) = (1.0 - n_0air) - m_w(1)        ! mass fraction of magma
m_a(1) = n_0air                         ! mass fraction air

! Calculate new mixture temperature
hmix = ( m_m(1) * plumeria_h_m(T_mag, p(1)) ) +                                &
       ( m_gas * plumeria_h_v(T_mag) )        +                                &
       ( m_a(1) * plumeria_h_a(T_mag) )       +                                &
       ( mw * plumeria_h_l(T_water) )
CALL plumeria_FindT_init(m_m(1), m_a(1), m_w(1), hmix, p(1), T_mix(1), m_v(1), &
                         m_l(1), m_i(1))
rho_mix(1) = 1 / ( ( m_l(1) / rho_w )         +                                &
                   ( m_m(1) / rho_m )         +                                &
                   ( m_i(1) / rho_ice )       +                                &
                   ( m_a(1) * R_dry + m_v(1) * R_vap ) * T_mix(1) / p(1) )

ratio_ur = 0.002
r(1) = ( mdot / rho_mix(1) / pi *                                              &
         SQRT(ratio_ur) ) ** (2.0 / 5.0)  ! initial vent radius
u(1) = SQRT(r(1) / ratio_ur)              ! initial exit velocity
gperm3_l(1) = 1000.0 * m_l(1) * rho_mix(1)
gperm3_i(1) = 1000.0 * m_i(1) * rho_mix(1)

! ----------------------------------------------------------------------
! Set initial conditions and trap errors
! ----------------------------------------------------------------------
istep = 1
sstep = 1.0
IF (r(1) < 10.0) THEN
  sstep = r(1) / 10.0
END IF
hnext = sstep
ux(1) = 0.0
uy(1) = 0.0
uz(1) = u(1)
TIME(1) = 0.0
phi(1) = pi / 2.0                 !initial plume angle from horizontal
idensflag = 0                     !set this flag to zero at the beginning
uzlast = 0.0
count_inversion = 0.0

! Calculate mixture sound speed
c_mix = SQRT(gamma_gas * p(1) / rho_mix(1))

! Trap possible errors
IF (u(1) < 5.0) THEN
  cmessage = 'This model can!t handle such low exit velocities.'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Make sure there's enough momentum to get the plume to the first dz step
du = -( ( rho_mix(1) - ( rho_air(1) / rho_mix(1) ) ) * gravity / u(1) +        &
        ( 2.0 * alpha * u(1) / r(1)) * (rho_air(1) / rho_mix(1)) ) * sstep
IF ((u(1) + du) < 0.0) THEN
  cmessage = 'Insufficient momentum to lift the plume. Program stopped'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! ----------------------------------------------------------------------
! Begin calculations
! ----------------------------------------------------------------------
begin_calc: DO
  ! End loop if the timestep > maximum timestep
  IF (istep < maxsteps) THEN

    yarr(1) = pi * r(istep) ** 2 *                                             &
               u(istep) * rho_mix(istep)             ! column mass flux
    yarr(2) = pi * r(istep) ** 2 *                                             &
               u(istep) * ux(istep) * rho_mix(istep) ! column momentum,x
    yarr(3) = pi * r(istep) ** 2 *                                             &
               u(istep) * uy(istep) * rho_mix(istep) ! column momentum,y
    yarr(4) = pi * r(istep) ** 2 *                                             &
               u(istep) * uz(istep) * rho_mix(istep) ! column momentum,z
    yarr(5) = yarr(1) * ( u(istep) ** 2 / 2 +                                  &
               ( gravity * z(istep) ) + hmix )       ! column tot energy
    yarr(6) = p(istep)                               ! column pressure
    yarr(7) = m_m(istep)                             ! mass frac magma
    yarr(8) = m_a(istep)                             ! mass frac air
    yarr(9) = m_w(istep)                             ! mass frac w.vapor
    yarr(10) = TIME(istep)                           ! time
    yarr(11) = x(istep)                              ! x, meters north
    yarr(12) = y(istep)                              ! y, meters east
    yarr(13) = z(istep)                              ! z, meters height

    !Values by which rkqs normalizes errors to check for accuracy
    yscale(1) = 100.0 * yarr(1)
    yscale(2) = 100.0 * (yarr(2) + yarr(3) + yarr(4))
    yscale(3) = 100.0 * (yarr(2) + yarr(3) + yarr(4))
    yscale(4) = 100.0 * (yarr(2) + yarr(3) + yarr(4))
    yscale(5) = 100.0 * yarr(5)
    yscale(6) = yarr(6)
    yscale(7) = 1.0
    yscale(8) = 1.0
    yscale(9) = 1.0
    yscale(10) = 100.0
    yscale(11)= 1000.0
    yscale(12)= 1000.0
    yscale(13)= 1000.0

    snow = s(istep)

    sstep = hnext

    CALL plumeria_derivs(snow,yarr,dydx)
    CALL plumeria_cashkarpqs(yarr, dydx, snow, snow_out,                       &
                             sstep, yscale, hdid, hnext)

    istep = istep + 1
    s(istep) = snow_out
    uzlast = uz(istep - 1)

    ! Assign plume properties
    p(istep)        = yarr(6)
    m_m(istep)      = yarr(7)
    m_a(istep)      = yarr(8)
    m_w(istep)      = yarr(9)
    TIME(istep)     = yarr(10)
    x(istep)        = yarr(11)
    y(istep)        = yarr(12)
    z(istep)        = yarr(13)
    m_m(istep)      = m_m(1) * mdot / yarr(1)

    ! Find temperature, mass fractions of aqueous phases at this step
    ! based on new enthalpy
    IF (z(istep) > zmax)  zmax = z(istep)
    ux(istep)       = yarr(2) / yarr(1)
    uy(istep)       = yarr(3) / yarr(1)
    uz(istep)       = yarr(4) / yarr(1)

    IF (IEEE_IS_NAN(uz(istep))) THEN
      cmessage = 'uz(istep) is nan. Program stopped.'
      errcode = 1
      CALL ereport(RoutineName,errcode,cmessage)
    END IF

    u(istep) = SQRT(ux(istep) ** 2 + uy(istep) ** 2 + uz(istep) ** 2)
    hmix = ( yarr(5) / yarr(1) ) -                                             &
            ( u(istep) ** 2 / 2 ) - ( gravity * z(istep) )
    CALL plumeria_FindT(m_m(istep), m_a(istep), m_w(istep), hmix, p(istep),    &
                        T_mix(istep), m_v(istep), m_l(istep), m_i(istep))
    rho_mix(istep) = 1.0 /                                                     &
                     ( ( m_m(istep) / rho_m )    +                             &
                       ( m_l(istep) / rho_w )    +                             &
                       ( m_i(istep) / rho_ice )  +                             &
                       ( m_a(istep) * R_dry      +                             &
                         m_v(istep) * R_vap) * T_mix(istep) / p(istep))
    r(istep)  = SQRT( yarr(1) / (pi * u(istep) * rho_mix(istep) ) )

    ! Atmospheric Properties
    T_air(istep)    = plumeria_AirTemp(z(istep))
    rho_air(istep)  = p(istep) / (R_dry * T_air(istep))
    gperm3_l(istep) = 1000.0 * m_l(istep) * rho_mix(istep)
    gperm3_i(istep) = 1000.0 * m_i(istep) * rho_mix(istep)

    ! Check if the plume is still rising
    IF (StillRising .AND. (uz(istep) < 0.0)) THEN
      StillRising = .FALSE.
      imax        = istep
      s_zmax      = s(istep)
    END IF

    IF ((stopatthetop) .AND. (uz(istep) < 0.0)) THEN
      WRITE(umMessage, '(A)') 'Stopped because uz < 0'
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    IF (ABS(x(istep)) > 600.0) THEN
      WRITE(umMessage,'(A)') 'Stopped because x > 600m'
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    IF ((s(istep)/s_zmax) > pathmax) THEN
      WRITE(umMessage,'(A)') 'Stopped because s/z > pathmax'
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    ! Calculate difference in density between plume and air
    rho_diff(1) = rho_mix(1) - rho_air(1)
    rho_diff(istep) = rho_mix(istep) - rho_air(istep)

    ! Calculate the number of inversion level and save timestep of inversion
    IF ((rho_diff(istep) > 0) .AND. (rho_diff(istep - 1) < 0)) THEN
      count_inversion = count_inversion + 1
      hb_istep = istep
    ELSE IF ((rho_diff(istep) < 0) .AND. (rho_diff(istep - 1) > 0)) THEN
      count_inversion = count_inversion + 1
      hb_istep = istep
    END IF

    ! Exit the calculation if istep reached maxsteps
  ELSE IF (istep >= maxsteps) THEN
    EXIT begin_calc
  END IF
END DO begin_calc

! ----------------------------------------------------------------------
! Output Neutral Buoyancy Level based on the inversion levels
! ----------------------------------------------------------------------
IF ((rho_diff(1) > 0) .AND. (count_inversion == 1)) THEN
  z_hb = zmax / 1000.0
ELSE IF ((rho_diff(1) > 0) .AND. (count_inversion == 2)) THEN
  z_hb = z(hb_istep)/1000.0
ELSE IF ((rho_diff(1) > 0) .AND. (count_inversion == 0)) THEN
  z_hb = zmax / 1000.0
ELSE IF ((rho_diff(1) < 0) .AND. (count_inversion == 1)) THEN
  z_hb = z(hb_istep) / 1000.0
ELSE IF ((rho_diff(1) > 0) .AND. (count_inversion > 2)) THEN
  z_hb = z(hb_istep) / 1000.0
END IF

IF (istep == maxsteps) THEN
  WRITE(umMessage,'(A,I0)') 'Maximum number of integration steps reached:',    &
                             maxsteps
  CALL umPrint(umMessage,src=RoutineName)
END IF

IF (u(istep) < umin) THEN
  WRITE(umMessage,'(A)') 'u(istep) < umin'
  CALL umPrint(umMessage,src=RoutineName)
END IF

WRITE(umMessage,'(A,F0.5)') 'Neutral buoyancy level (in km) =', z_hb
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)') 'Successful completion of Plumeria.'
CALL umPrint(umMessage,src=RoutineName)

END SUBROUTINE plumeria_main
END MODULE plumeria_main_mod
