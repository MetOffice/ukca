! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module providing data structures to specify details of the active
!   photolysis configuration:
!     photol_config  - contains the configuration variables i.e. control/
!                      scheme choices, static scalar values
!
!   The module also provides the following procedure for the Photolysis_API
!     photol_get_config - returns values for photolysis configuration variables.
!
!   The following additional public procedures are provided for use
!     init_photol_config   - initialises/resets all configuration
!                            data ready for a new photolysis configuration to
!                            be set up
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_Photolysis
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE photol_config_specification_mod

USE missing_data_mod, ONLY: imdi, rmdi

IMPLICIT NONE
PUBLIC

! ---------------------------------------------------------------------------
! -- Type for holding Photolysis configuration data --
! ---------------------------------------------------------------------------

! Components in each section other than internal variables are set by the
! parent application via optional keyword arguments with matching names
! in 'photol_setup'.
! In the event of any change to keyword argument names, old components may
! be retained in the structure as duplicates of the new input variables
! (to support retention of old keyword argument names for backward
! compatibility and/or to avoid the need to change variable names throughout
! the photolysis code), but should be re-categorised as internal variables.

TYPE :: photol_config_spec_type

  ! -- Context information --
  ! Includes array sizes for input fields, indices for accessing specific
  ! items in the input arrays, model timestep etc

  ! ----- -- Integer items ------------------------------
  INTEGER :: chem_timestep           ! Chemical timestep in seconds
                                     ! (for N-R and Offline oxidant schemes)

  INTEGER :: fastjx_mode             ! Method to use above cutoff pressure level
                                     ! 1 = use only 2D/ lookup table approach
                                     ! 2 = merge 2D and FastJX,
                                     ! 3 = Use only FastJX
  INTEGER :: fastjx_numwl            ! No. of wavelengths to use (8, 12, 18)

  INTEGER :: global_row_length       ! Points per global row

  INTEGER :: i_photol_scheme         ! Photolysis scheme to use. Choices are
                                     ! 0 = Photolysis Off;
                                     ! 1 = Stratospheric photolysis only
                                     ! 2 = Offline/ 2-D photolysis scheme
                                     ! 3 = FastJX scheme

  INTEGER :: i_solcylc_type          ! Solar cycle option for Fast-JX
                                     ! 0 = no solar cycle
                                     ! 1 = use observed cycle for a range of
                                     !     years, give initial year as input
                                     !     and number of months in input file
                                     ! 2 = use an average cycle for all times

  INTEGER :: model_levels            ! Z dimension of domain (levels)

  INTEGER :: n_cca_lev               ! Number of Conv Cloud Amount levels


  INTEGER :: solcylc_start_year      ! 1st year for the solar cycle obs data

  INTEGER :: i_error_method          ! Governs process to follow on error, i.e.
                                     ! - Write error message and abort, or
                                     ! - Return to parent
                                     ! - Write error as warning and return
  ! ----- -- Logical items ------------------------------
  LOGICAL :: l_cal360                ! True if using a 360-day calendar

  LOGICAL :: l_cloud_pc2             ! True if the atmosphere model is
                                     ! using the PC2 cloud scheme
  LOGICAL :: l_3d_cca                ! True if Convective Cloud Amount is
                                     ! calculated as a 3-D field, i.e.
                                     ! loop over n_cca_lev, else use only
                                     ! first level from cca array.

  ! Use of external photolysis rates (e.g. from radiation scheme)
  LOGICAL :: l_environ_jo2           ! True if using ext O2 -> O(3P) rate
  LOGICAL :: l_environ_jo2b          ! True if using ext O2 -> O(1D) rate

  LOGICAL :: l_environ_ztop          ! True if z_top_of_model is being
                                     ! supplied by parent model
  LOGICAL :: l_enable_diag_um        ! True to enable diagnostic output via
                                     ! the Unified Model STASH system
  LOGICAL :: l_strat_chem            ! True if UKCA is using a scheme that
                                     ! involves stratospheric chemistry i.e.
                                     ! Std Strat, Strattrop or CRI-strat

  ! ----- -- Real items ------------------------------
  REAL    :: fastjx_prescutoff       ! Pressure level (Pa) above which a
                                     ! simplified/ 2D/ look-up table approach
                                     ! can be used for photolysis rates

  REAL    :: timestep                ! Model (dynamical) timestep in seconds
                                     ! (an integer no of timesteps is
                                     !  expected in 1 hour)

END TYPE photol_config_spec_type

TYPE(photol_config_spec_type),   SAVE :: photol_config

LOGICAL :: l_photol_config_available   ! Flag to determine if photolysis
                                       ! configuration has been set.

! ---------------------------------------------------------------------------
! -- Option values for configuration variables --
! ---------------------------------------------------------------------------

! Option codes for 'i_photol_scheme'
INTEGER, PARAMETER :: i_scheme_nophot = 0       ! photolysis off
INTEGER, PARAMETER :: i_scheme_photol_strat = 1 ! stratospheric photolysis only
INTEGER, PARAMETER :: i_scheme_phot2d = 2       ! offline 2D photolysis scheme
INTEGER, PARAMETER :: i_scheme_fastjx = 3       ! Fast-JX

! Oprion codes for 'fastjx_mode' method to use above prescutoff level
INTEGER, PARAMETER :: fjx_mode_2Donly = 1  ! use only 2D/ lookup table approach
INTEGER, PARAMETER :: fjx_mode_merged = 2  ! merge 2D and FastJX,
INTEGER, PARAMETER :: fjx_mode_fastjx = 3  ! Use only FastJX

! Option codes for 'i_solcylc_type'
INTEGER, PARAMETER :: i_no_solcylc = 0          ! No solar cycle applied
INTEGER, PARAMETER :: i_obs_solcylc = 1         ! Use observed cycle for a
                                                ! given range of years
INTEGER, PARAMETER :: i_avg_solcylc = 2         ! Use an average cycle for
                                                ! the whole run.

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE init_photol_configuration()
! ----------------------------------------------------------------------
! Description:
!   Initialises/resets the Photol configuration data ready for a new
!   configuration to be set up.
!
!   In addition, all configurable constants are set to their UM
!   default values.
! ----------------------------------------------------------------------

USE photol_constants_mod,  ONLY: const_pi, const_pi_over_180,                  &
                              const_recip_pi_over_180, const_rhour_per_day,    &
                              const_o3_mmr_vmr, const_molemass_sulp,           &
                              const_molemass_nh42so4, const_molemass_air,      &
                              const_planet_radius, const_g, const_r,           &
                              const_avogadro, const_tm, const_s2r

USE ukca_error_mod,        ONLY: i_error_method_abort

IMPLICIT NONE

l_photol_config_available = .FALSE.

! -- Set integer items --
photol_config%chem_timestep = imdi

photol_config%fastjx_mode = imdi
photol_config%fastjx_numwl = imdi

photol_config%global_row_length = imdi

photol_config%i_photol_scheme  = imdi
photol_config%i_solcylc_type = imdi

photol_config%model_levels = imdi
photol_config%n_cca_lev = imdi

photol_config%solcylc_start_year = imdi

photol_config%i_error_method = i_error_method_abort

! -- Set Logicals ---
photol_config%l_cal360 = .FALSE.

photol_config%l_cloud_pc2 = .FALSE.
photol_config%l_3d_cca = .FALSE.

photol_config%l_enable_diag_um = .FALSE.

photol_config%l_environ_jo2 = .FALSE.
photol_config%l_environ_jo2b = .FALSE.

photol_config%l_environ_ztop = .FALSE.
photol_config%l_strat_chem = .FALSE.

! -- Set Real Items ---
photol_config%fastjx_prescutoff = rmdi

photol_config%timestep = rmdi

! -- Constants used by Photolysis --
! Set to UM values by default.

! Configurable by parent
const_pi = 3.14159265358979323846
const_pi_over_180      = const_pi / 180.0
const_recip_pi_over_180 = 180.0 / const_pi

const_o3_mmr_vmr       = 1.657
const_molemass_sulp    = 32.07
const_molemass_nh42so4 = 132.16
const_molemass_air     = 0.02897

const_planet_radius    = 6371229.0
! Seconds-to-radians calculated as :
! earth_dha (Increment to Earth's hour angle per day number from epoch) /
!  rsec_per_day
const_s2r             = (2.0*const_pi)/ 86400.0

! Not currently configurable
const_rhour_per_day    = 24.0
const_g                = 9.80665
const_r                = 287.05
const_avogadro         = 6.022e23
const_tm               = 273.15

RETURN
END SUBROUTINE init_photol_configuration

! ----------------------------------------------------------------------
SUBROUTINE photol_get_config(                                                  &
   fastjx_mode, fastjx_numwavel, i_photol_scheme, i_solcylc_type,              &
   solcylc_start_year, l_config_available, l_environ_jo2, l_environ_jo2b,      &
   l_strat_chem, fastjx_prescutoff                                             &
   )
! ----------------------------------------------------------------------
! Description:
!   Returns values of configuration variables that are supplied as
!   optional arguments.
!   These are values that the parent or a coupled model might need to know
!    for internal scheme or processing choices.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
! (follow order as defined in 'photol_config_spec_type' group)

INTEGER, OPTIONAL, INTENT(OUT) :: fastjx_numwavel
INTEGER, OPTIONAL, INTENT(OUT) :: fastjx_mode

INTEGER, OPTIONAL, INTENT(OUT) :: i_photol_scheme

INTEGER, OPTIONAL, INTENT(OUT) :: i_solcylc_type
INTEGER, OPTIONAL, INTENT(OUT) :: solcylc_start_year

LOGICAL, OPTIONAL, INTENT(OUT) :: l_config_available

LOGICAL, OPTIONAL, INTENT(OUT) :: l_environ_jo2
LOGICAL, OPTIONAL, INTENT(OUT) :: l_environ_jo2b

LOGICAL, OPTIONAL, INTENT(OUT) :: l_strat_chem

REAL, OPTIONAL, INTENT(OUT)    :: fastjx_prescutoff


! -- Context information --
IF (PRESENT(fastjx_numwavel)) fastjx_numwavel = photol_config%fastjx_numwl
IF (PRESENT(fastjx_mode)) fastjx_mode = photol_config%fastjx_mode
IF (PRESENT(i_photol_scheme)) i_photol_scheme = photol_config%i_photol_scheme

IF (PRESENT(i_solcylc_type)) i_solcylc_type = photol_config%i_solcylc_type
IF (PRESENT(solcylc_start_year)) solcylc_start_year =                          &
                                          photol_config%solcylc_start_year

IF (PRESENT(l_environ_jo2)) l_environ_jo2 = photol_config%l_environ_jo2
IF (PRESENT(l_environ_jo2b)) l_environ_jo2b = photol_config%l_environ_jo2b
IF (PRESENT(l_strat_chem)) l_strat_chem = photol_config%l_strat_chem
IF (PRESENT(fastjx_prescutoff))                                                &
                       fastjx_prescutoff = photol_config%fastjx_prescutoff

! -- Availability of a valid configuration --
IF (PRESENT(l_config_available)) l_config_available = l_photol_config_available

RETURN
END SUBROUTINE photol_get_config

END MODULE photol_config_specification_mod
