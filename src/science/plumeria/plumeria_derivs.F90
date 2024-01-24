! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module containing the subroutine PLUMERIA_DERIVS.
!
!  Methods:
!    The subroutine PLUMERIA_DERIVS calculates the mass flux, momentum
!    flux and total energy flux at a given point along the plume's path.
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

MODULE plumeria_derivs_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PLUMERIA_DERIVS_MOD'

CONTAINS

SUBROUTINE plumeria_derivs(x, yarr, dydx)

USE plumeria_param_mod,     ONLY: rho_m, rho_w, rho_ice, n_exp, gamma_w,       &
                                  alpha, idensflag, kgmole_w, kgmole_air
USE plumeria_functions_mod, ONLY: plumeria_psat, plumeria_AirTemp,             &
                                  plumeria_AirHumid, plumeria_zeta,            &
                                  plumeria_wind, plumeria_vx, plumeria_vy,     &
                                  plumeria_h_a, plumeria_h_v
USE plumeria_FindT_mod,     ONLY: plumeria_FindT
USE conversions_mod,        ONLY: pi
USE comorph_constants_mod,  ONLY: R_dry, R_vap, gravity
USE parkind1,               ONLY: jpim, jprb      ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Input variables
REAL, INTENT(IN) :: x, yarr(13)

! Return value
REAL, INTENT(OUT) :: dydx(13)


! Local variables
! Real variables in input array, yarr(13)
REAL :: Massflux            ! Column mass, yarr(1)
REAL :: Mx, My, Mz          ! Column momentum in x, y, z directions, yarr(2:4)
REAL :: pnow                ! Column pressure, yarr(5)
REAL :: m_mnow              ! Mass fraction of magma in mixture column, yarr(6)
REAL :: m_anow              ! Mass fraction of dry air in mixture column, yarr(7)
REAL :: m_wnow              ! Mass fraction of water in mixture column, yarr(8)
REAL :: TimeNow             ! Time, yarr(9)
REAL :: xnow, ynow, znow    ! Distance in x, y, z directions (in m), yarr(10:13)

REAL :: dQds                ! Mass flux, dydx(1)
REAL :: dMxds, dMyds, dMzds ! Moment components in x, y, z directions, dydx(2:4)
REAL :: h_airnow            ! Enthalpy in air
REAL :: h_mixnow            ! Enthalpy of mixture column
REAL :: HumidityNow         ! Humidity in the atmosphere
REAL :: m_inow              ! Mass fraction of ice in mixture column
REAL :: m_lnow              ! Mass fraction of liquid in mixture column
REAL :: ma_airnow           ! Mass fraction in ambient air
REAL :: mv_airnow           ! Mass fraction water vapor in ambient air
REAL :: m_vnow              ! Mass fraction water vapor in column
REAL :: w_airnow            ! Specific humidity in air
REAL :: ux, uy, uz          ! Velocity components in x, y, z directions
REAL :: unow                ! Magnitude of Velocity
REAL :: sigma               ! Azimuth angle of velocity vector
REAL :: phinow              ! Angle between vector and xy plane
REAL :: MomentumFlux_now    ! Momentum flux in the mixture column
REAL :: rnow                ! Radius of mixture column
REAL :: rho_airnow          ! Air density, kg/m3
REAL :: rho_mixnow          ! Density of mixture column
REAL :: T_airnow            ! Temperature of air
REAL :: T_mixnow            ! Temperature of mixture column
REAL :: v_perpnow           ! Maximum wind component perpendicular to plume
REAL :: v_s1, v_s2, v_s3    ! Components of wind parallel and perpendicular
                            ! to the plume


INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PLUMERIA_DERIVS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Get variables from yarr
Massflux = yarr(1)
Mx       = yarr(2)
My       = yarr(3)
Mz       = yarr(4)
pnow     = yarr(6)
m_mnow   = yarr(7)
m_anow   = yarr(8)
m_wnow   = yarr(9)
TimeNow  = yarr(10)
xnow     = yarr(11)
ynow     = yarr(12)
znow     = yarr(13)

! Calculate additional variables
ux       = yarr(2) / yarr(1)
uy       = yarr(3) / yarr(1)
uz       = yarr(4) / yarr(1)
phinow   = ATAN2(uz, SQRT(ux ** 2 + uy ** 2))

IF ((ux /= 0.0) .OR. (uy /= 0.0)) THEN
  sigma = ATAN2(uy, ux)
ELSE
  sigma = pi / 2.0 - plumeria_zeta(znow)
END IF

unow     = SQRT(ux ** 2 + uy ** 2 + uz ** 2)
h_mixnow = (yarr(5) / Massflux) - ( unow ** 2 / 2.0 ) - ( gravity * znow )

!when u-->0, MassFlux-->0, and h_mixnow and unow can get very unstable.
!the lines below stabilize these values until u increases enough to stabilize
!the calculations
IF (h_mixnow < 0.0) THEN
  uz = 0.001
  unow = SQRT(ux ** 2 + uy ** 2 + uz ** 2)
  h_mixnow = (yarr(5) / Massflux) - ( unow ** 2 / 2.0 ) - ( gravity * znow )
END IF

MomentumFlux_now = MassFlux * unow

! Calculate Component Properties
T_airnow = plumeria_AirTemp(znow)
HumidityNow = plumeria_AirHumid(znow)
rho_airnow = pnow / (R_dry * T_airnow)
w_airnow = ( kgmole_w / kgmole_air ) *                                         &
           ( plumeria_psat(T_airnow) * HumidityNow / pnow )
ma_airnow = 1.0 / (1.0 + w_airnow)
mv_airnow = w_airnow / (1.0 + w_airnow)
h_airnow = ( ma_airnow * plumeria_h_a(T_airnow) ) +                            &
           ( mv_airnow * plumeria_h_v(T_airnow) )

! Calculate temperature, radius, density at current elevation
! Find temperature
CALL plumeria_FindT(m_mnow, m_anow, m_wnow, h_mixnow, pnow, T_mixnow,          &
                    m_vnow, m_lnow, m_inow)
! Find density
rho_mixnow = 1 / ( ( m_lnow / rho_w)     +                                     &
                   ( m_mnow / rho_m)     +                                     &
                   ( m_inow / rho_ice )  +                                     &
                   ( m_anow * R_dry + m_vnow * R_vap ) * T_mixnow / pnow )
! Find radius
rnow = SQRT(yarr(1) / (pi * unow * rho_mixnow))

! Calculate components of wind parallel and perpendicular to the plume
v_s1 = ( plumeria_vx(znow) * COS(sigma) +                                      &
         plumeria_vy(znow) * SIN(sigma) ) * COS(phinow)
v_s2 = ( -plumeria_vx(znow) * SIN(sigma) +                                     &
         plumeria_vy(znow) * COS(sigma) )
v_s3 = ( plumeria_vx(znow) * COS(sigma) +                                      &
         plumeria_vy(znow) * SIN(sigma) ) * SIN(phinow)
v_perpnow = SQRT(v_s2 ** 2 + v_s3 ** 2)

! Calculate new values of dydx
IF ((idensflag == 0) .AND. (rho_mixnow > rho_airnow)) THEN
  !Used in Woods 1993 & later
  dydx(1) = 2.0 * pi * rnow *                                                  &
            ( ( alpha * ABS(unow - v_s1)) ** n_exp +                           &
              ( gamma_w * v_perpnow ) ** n_exp ) ** (1.0 / n_exp) *            &
              SQRT(rho_airnow * rho_mixnow)
ELSE
  idensflag = 1
  dydx(1) = 2.0 * rho_airnow * pi * rnow *                                     &
            ( ( alpha * ABS(unow - v_s1)) ** n_exp +                           &
              ( gamma_w * v_perpnow) ** n_exp ) ** (1.0 / n_exp)
END IF

dQds = dydx(1)

! Calculate moment components in x, y, z direction
dydx(2) = dQds * plumeria_vx(znow)                              !dMx/ds
dydx(3) = dQds * plumeria_vy(znow)                              !dMy/ds
dydx(4) = pi * rnow ** 2 * gravity * (rho_airnow - rho_mixnow)  !dMz/ds
dMxds   = dydx(2)
dMyds   = dydx(3)
dMzds   = dydx(4)
dydx(5) = dQds * ( h_airnow + plumeria_wind(znow) ** 2 / 2.0 +                 &
           (gravity * znow) )                                   !energy grad
dydx(6) = -rho_airnow * gravity * SIN(phinow)                   !pressure grad
dydx(7) = -(m_mnow / Massflux) * dQds                           !dm_m/ds
dydx(8) = (ma_airnow - m_anow) * dQds / Massflux                !dm_a/ds
dydx(9) = (mv_airnow - m_wnow) * dQds / Massflux                !dm_w/ds
dydx(10) = 1.0 / unow                                           !dt/ds
dydx(11) = ux / unow                                            !dx/ds
dydx(12) = uy / unow                                            !dy/ds
dydx(13) = uz / unow                                            !dz/ds

END SUBROUTINE plumeria_derivs
END MODULE plumeria_derivs_mod
