! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module declaring parameters used in Plumeria
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

MODULE plumeria_param_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PLUMERIA_PARAM_MOD'

! Logical
LOGICAL       :: stopatthetop !.true. to stop model at the maximum elevation

! Flag
INTEGER       :: idensflag    !Set to zero at beginning of each run
                              !then set to 1 once we cross into the
                              !convective thrust region.

! Variable names used in Plumeria for input arguments from UKCA
REAL          :: TairLayer(85)       ! Temperature in atm layers
REAL          :: Zairlayer(85)       ! Geopotential height in atm layers
REAL          :: HumidAirLayer(85)   ! Relative humidity in atm layers
REAL          :: pAirLayer(85)       ! Pressure in atm layers
REAL          :: WinddirLayer(85)    ! Wind direction in atm layers
REAL          :: WindspeedLayer(85)  ! Wind speed in atm layers
REAL          :: mdot                ! Mass eruption rate
REAL          :: vent_elevation      ! Vent altitude
INTEGER       :: iatmlayers          ! Number of atmospheric layers

! Velocity
REAL          :: u_0           ! initial velocity at vent
REAL          :: c_mix         ! mixture sound speed

! Mass fraction
REAL          :: m_gas         ! mass fraction magmatic gas
REAL          :: mw            ! mass fraction water added to erupting mixture
REAL          :: n_0           ! gas mass fraction at the vent
REAL          :: n_0air        ! gas mass fraction at the vent,
                               ! assuming all gas is air.

! Integration variables
REAL          :: dydx(13)      ! variable used for integration
REAL          :: yarr(13)      ! variable used for integration
REAL          :: hnext         ! variable used for integration
REAL          :: hdid          ! variable used for integration
REAL          :: zstep         ! vertical step of integration

! Coefficients
REAL          :: alpha         ! Entrainment coefficient
REAL          :: gamma_w       ! Crossflow entrainmnent coefficient
REAL          :: n_exp         ! Devenish exponent to entrainment equation

! Molar weight
REAL          :: kgmole_air    ! molar weight of air (kg/mole)
REAL          :: kgmole_w      ! molar weight of water (kg/mole)

! Enthalpy
REAL          :: hmix          ! mixture enthalpy

! Magma properties
REAL          :: Cp_m          ! Magma specific heat, J/kg K
REAL          :: gamma_gas     ! cp/cv for magmatic gas

! Density
REAL          :: rho_0         ! bulk density of mixture at the vent
REAL          :: rho_m         ! magma density (kg/m3)
REAL          :: rho_w         ! density of water (1,000 kg/m3)
REAL          :: rho_ice       ! density of ice

! Temperature
REAL          :: T_mag          ! magma temperature (K)
REAL          :: T_water        ! ambient water temperature
REAL          :: T_ColdWater    ! Top temperature range which
                                ! liquid water & ice coexist
REAL          :: T_ice          ! Bottom temperature range which
                                ! liquid water & ice coexist

!Real 1-D variable arrays that vary with elevation:
!Dimension of variables set as 10000, = maxsteps in plume_main

! Mass fraction
REAL         :: m_a(10000)           ! mass fraction dry air in column
REAL         :: m_i(10000)           ! mass fraction ice in water column
REAL         :: m_l(10000)           ! mass fraction liquid water in column
REAL         :: m_m(10000)           ! mass fraction magma in column
REAL         :: m_w(10000)           ! total water in column
REAL         :: m_v(10000)           ! mass fraction water vapor in column

! Density
REAL         :: rho_air(10000)       ! air density (kg/m3)
REAL         :: rho_mix(10000)       ! bulk density (kg/m3)
REAL         :: rho_diff(10000)      ! difference in density between rho_air
                                     ! and rho_mix (kg/m3)

! Temperature
REAL         :: T_air(10000)         ! air temperature (K)
REAL         :: T_mix(10000)         ! temperature of column (K)

! Air pressure
REAL         :: p(10000)             ! air pressure (Pascals)

! Plume-related variables
REAL         :: TIME(10000)          ! Time (s) after leaving the vent
REAL         :: phi(10000)           ! angle of plume from horizontal
REAL         :: r(10000)             ! plume radius (m)
REAL         :: s(10000)             ! Distance along plume axis
REAL         :: u(10000)             ! velocity (m/s)
REAL         :: ux(10000)            ! velocity in x direction
REAL         :: uy(10000)            ! velocity in y direction
REAL         :: uz(10000)            ! velocity in z direction
REAL         :: x(10000)             ! Distance in x direction
REAL         :: y(10000)             ! Distance in y direction
REAL         :: z(10000)             ! Distance in z direction,
                                     ! height above the vent

END MODULE plumeria_param_mod
