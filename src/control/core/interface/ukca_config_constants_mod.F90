! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module providing various constants for UKCA that are potentially
!   configurable via the ukca_setup API call.
!
!   The following public procedure is also provided for use within UKCA
!     init_config_constants  - initialises/resets all UKCA configurable
!                              constants ready for a new UKCA configuration to
!                              be set up
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_config_constants_mod

USE ukca_missing_data_mod, ONLY: imdi, rmdi

IMPLICIT NONE
PUBLIC

! ---------------------------------------------------------------------------
! -- Configurable set of constants used in UKCA --
! ---------------------------------------------------------------------------
! These are constants that can potentially be overridden by a parent
! application (e.g. if required for consistency with values used elsewhere)
! or are derived from those constants.
! Values are set in the routine 'init_config_constants' defined below.

!!!! Override functionality via the 'ukca_setup' argument list is not yet
!!!! available for all of these values.

!!!! For some:
! UM values are currently provided via 'ukca_um_legacy_mod'.
! These replace the default values from this module in UM builds.

! Molar universal gas constant (J/K/mol)
REAL, SAVE :: rmol = rmdi

! Planet Constants (Earth values by default)
REAL, SAVE :: planet_radius = rmdi  ! Planet radius (m)
REAL, SAVE :: g = rmdi              ! Mean acceleration due to gravity at
                                    ! surface (m/s^2)
REAL, SAVE :: r = rmdi              ! Gas constant for dry air (J/kg/K)
REAL, SAVE :: cp = rmdi             ! Specific heat of dry air at constant
                                    ! pressure (J/kg/K)
REAL, SAVE :: pref = rmdi           ! Reference surface pressure (Pa)
REAL, SAVE :: vkman = rmdi          ! Von Karman's constant
REAL, SAVE :: repsilon = rmdi       ! Ratio of molecular weights of water
                                    ! and dry air
REAL, SAVE :: rv = rmdi             ! Gas constant for water vapour (J/kg/K)
REAL, SAVE :: kappa = rmdi          ! r/cp
REAL, SAVE :: c_virtual = rmdi      ! (1/repsilon) - 1

! Water constants
REAL, SAVE :: tfs = rmdi            ! Temperature at which sea water freezes (K)
REAL, SAVE :: rho_water = rmdi      ! Density of pure water (kg/m3)
REAL, SAVE :: rhosea = rmdi         ! Density of sea water
REAL, SAVE :: lc = rmdi             ! Latent heat of condensation of water at
                                    ! 0 deg C

! Chemistry constants
REAL, SAVE :: avogadro = rmdi       ! No.of molecules in 1 mole
REAL, SAVE :: boltzmann = rmdi      ! Boltzmann's constant (J K-1)
REAL, SAVE :: rho_so4 = rmdi        ! Density of SO4 particle (kg/m3)

! CLASSIC aerosol parameters for heterogeneous PSC chemistry
REAL, SAVE :: rad_ait = rmdi  ! Mean radius of Aitken mode particles (m)
REAL, SAVE :: rad_acc = rmdi  ! Mean radius of accumulation mode particles (m)
REAL, SAVE :: chi = rmdi      ! Mole fraction of S in particle
REAL, SAVE :: sigma = rmdi    ! Std. dev. of accumulation mode particle size
                              ! distribution

! ---------------------------------------------------------------------------
! -- Flag to indicate whether configurable UKCA constants are set up --
! ---------------------------------------------------------------------------

LOGICAL, SAVE :: l_ukca_constants_available = .FALSE.

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_CONFIG_CONSTANTS_MOD'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE init_config_constants()
! ----------------------------------------------------------------------
! Description:
!   Initialises/resets the UKCA configurable constants to their default
!   values ready for a new configuration to be set up.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Planet Constants (Earth values by default)
planet_radius = 6371229.0
g = 9.80665
r = 287.05
cp = 1005.0
pref = 100000.0
vkman = 0.4
repsilon = 0.62198
rv = r / repsilon
kappa = r / cp
c_virtual = 1.0 / repsilon - 1.0
! Molar universal gas constant (J/K/mol)
rmol = 8.314
! Water constants
tfs = 271.35
rho_water = 1000.0
rhosea = 1026.0
lc = 2.501e6
! Chemistry constants
avogadro = 6.022e23
boltzmann = 1.3804e-23
rho_so4 = 1769.0
! CLASSIC aerosol parameters for heterogeneous PSC chemistry
rad_ait = 6.5e-9
rad_acc = 95.0e-9
chi = 32.0 / 132.0
sigma = 1.4

END SUBROUTINE init_config_constants

END MODULE ukca_config_constants_mod
