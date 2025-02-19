!*****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module containing subroutine photol_setup for receiving configuration
!   options for photolysis and setting up internal data to
!   define the configuration.
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

MODULE photol_setup_mod

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='PHOTOL_SETUP_MOD'

! Public procedures
PUBLIC :: photol_setup

CONTAINS

! ----------------------------------------------------------------------

SUBROUTINE photol_setup(i_photol_scheme,                                       &
                      error_code,                                              &
                      chem_timestep,                                           &
                      fastjx_mode,                                             &
                      fastjx_numwl,                                            &
                      global_row_length,                                       &
                      i_solcylc_type,                                          &
                      ip_aerosol_param_moist,                                  &
                      ip_accum_sulphate,                                       &
                      ip_aitken_sulphate,                                      &
                      model_levels,                                            &
                      n_cca_lev,                                               &
                      solcylc_start_year,                                      &
                      i_error_method,                                          &
                      l_cal360,                                                &
                      l_cloud_pc2,                                             &
                      l_3d_cca,                                                &
                      l_enable_diag_um,                                        &
                      l_environ_jo2,                                           &
                      l_environ_jo2b,                                          &
                      l_environ_ztop,                                          &
                      l_strat_chem,                                            &
                      fastjx_prescutoff,                                       &
                      timestep,                                                &
                      pi,                                                      &
                      o3_mmr_vmr,                                              &
                      molemass_sulp,                                           &
                      molemass_nh42so4,                                        &
                      molemass_air,                                            &
                      planet_radius,                                           &
                      error_message, error_routine)

! ----------------------------------------------------------------------
! Description:
!
!  Given the input configuration control data, check its validity and
!  set up the Photolysis internal configuration data accordingly.
!  This includes some basic initialisation and everything required to
!  establish details of the selected configuration that will determine
!  the required environmental drivers.
!
! Method:
!
!  1. Copy the configuration variables provided as keyword arguments into
!     component variables with matching names in the photol_config structures.
!     For certain variables, default values are set here for use if input
!     values are not provided.
!     ------------------------------------------------------------------
!     Note: Values of variables that are inactive in the current
!     configuration will normally be ignored and not set.
!     ------------------------------------------------------------------
!  2. Check that configuration values are consistent.
!  3. FUTURE: Initialise photolysis rate names.
!  4. FUTURE: Set up lists of environmental driver fields required.
!  5. FUTURE: Initialise master diagnostics list and determine availability of
!     each diagnostic given the configuration.
!
! ----------------------------------------------------------------------

USE photol_config_specification_mod, ONLY: photol_config,                      &
       l_photol_config_available, i_scheme_nophot, i_scheme_photol_strat,      &
       i_scheme_phot2d, i_scheme_fastjx, i_obs_solcylc, i_avg_solcylc,         &
       init_photol_configuration

USE photol_constants_mod,  ONLY: const_pi, const_pi_over_180,                  &
                              const_recip_pi_over_180,                         &
                              const_o3_mmr_vmr, const_molemass_sulp,           &
                              const_molemass_nh42so4, const_molemass_air,      &
                              const_planet_radius, const_s2r

USE photol_environment_mod, ONLY: photol_init_environ_req

USE ukca_error_mod,         ONLY: maxlen_message, maxlen_procname,             &
                                  error_report, errcode_value_unknown,         &
                                  errcode_value_invalid

USE parkind1,               ONLY: jpim, jprb      ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Subroutine arguments.

! Except for the top-level photolysis scheme choice, each input configuration
! variable is an optional keyword argument with a matching component in the
! Photolysis configuration structure.
! Based on the scheme choices, certain variables may be expected to be
!  always defined by the parent as they are accessed by default in the
!  current photolysis workflow and not setting these can lead to unexpected
!  behaviour.

INTEGER, INTENT(IN) :: i_photol_scheme
INTEGER, TARGET, INTENT(OUT) :: error_code

! Optional arguments
INTEGER, OPTIONAL, INTENT(IN) :: chem_timestep

INTEGER, OPTIONAL, INTENT(IN) :: fastjx_mode
INTEGER, OPTIONAL, INTENT(IN) :: fastjx_numwl
INTEGER, OPTIONAL, INTENT(IN) :: global_row_length

INTEGER, OPTIONAL, INTENT(IN) :: i_solcylc_type

INTEGER, OPTIONAL, INTENT(IN) :: ip_aerosol_param_moist
INTEGER, OPTIONAL, INTENT(IN) :: ip_accum_sulphate
INTEGER, OPTIONAL, INTENT(IN) :: ip_aitken_sulphate

INTEGER, OPTIONAL, INTENT(IN) :: model_levels
INTEGER, OPTIONAL, INTENT(IN) :: n_cca_lev
INTEGER, OPTIONAL, INTENT(IN) :: solcylc_start_year
INTEGER, OPTIONAL, INTENT(IN) :: i_error_method

LOGICAL, OPTIONAL, INTENT(IN) :: l_cal360

LOGICAL, OPTIONAL, INTENT(IN) :: l_cloud_pc2
LOGICAL, OPTIONAL, INTENT(IN) :: l_3d_cca
LOGICAL, OPTIONAL, INTENT(IN) :: l_enable_diag_um

LOGICAL, OPTIONAL, INTENT(IN) :: l_environ_jo2
LOGICAL, OPTIONAL, INTENT(IN) :: l_environ_jo2b
LOGICAL, OPTIONAL, INTENT(IN) :: l_environ_ztop
LOGICAL, OPTIONAL, INTENT(IN) :: l_strat_chem

REAL, OPTIONAL, INTENT(IN)    :: fastjx_prescutoff
REAL, OPTIONAL, INTENT(IN)    :: timestep
! Configurable constants
REAL, OPTIONAL, INTENT(IN)    :: pi
REAL, OPTIONAL, INTENT(IN)    :: o3_mmr_vmr
REAL, OPTIONAL, INTENT(IN)    :: molemass_sulp
REAL, OPTIONAL, INTENT(IN)    :: molemass_nh42so4
REAL, OPTIONAL, INTENT(IN)    :: molemass_air
REAL, OPTIONAL, INTENT(IN)    :: planet_radius

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables
INTEGER, POINTER :: error_code_ptr
CHARACTER(LEN=maxlen_message) :: err_message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHOTOL_SETUP'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Set defaults for output arguments
error_code_ptr => error_code
error_code_ptr = 0
err_message = ''
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Set all configuration data to default values
CALL init_photol_configuration()

! First, check if parent has specified a method for error handling, in case of
! any errors further on.
IF (PRESENT(i_error_method)) photol_config%i_error_method = i_error_method

! Check that a known photolysis scheme is specified.
IF ( i_photol_scheme /= i_scheme_nophot        .AND.                           &
     i_photol_scheme /= i_scheme_photol_strat  .AND.                           &
     i_photol_scheme /= i_scheme_phot2d        .AND.                           &
     i_photol_scheme /= i_scheme_fastjx ) THEN
  error_code_ptr = errcode_value_unknown
  WRITE(err_message, '(A,I0)') 'Unknown Photolysis scheme specified ',         &
    i_photol_scheme
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_message, &
         RoutineName, msg_out=error_message, locn_out=error_routine)

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
ELSE
  photol_config%i_photol_scheme = i_photol_scheme
END IF

! Collate input data specifying the UKCA configuration
! (loosely following 'photol_config_spec_type'

IF (PRESENT(chem_timestep))  photol_config%chem_timestep = chem_timestep
IF (PRESENT(fastjx_mode))    photol_config%fastjx_mode  = fastjx_mode
IF (PRESENT(fastjx_numwl))   photol_config%fastjx_numwl = fastjx_numwl

IF (PRESENT(global_row_length)) photol_config%global_row_length                &
                                                        = global_row_length
IF (PRESENT(i_solcylc_type)) photol_config%i_solcylc_type = i_solcylc_type

IF (PRESENT(model_levels)) photol_config%model_levels = model_levels
IF (PRESENT(n_cca_lev))    photol_config%n_cca_lev    = n_cca_lev

IF (PRESENT(solcylc_start_year)) photol_config%solcylc_start_year              &
                                                      = solcylc_start_year

IF (PRESENT(l_cal360))     photol_config%l_cal360     = l_cal360

IF (PRESENT(l_cloud_pc2))    photol_config%l_cloud_pc2 = l_cloud_pc2
IF (PRESENT(l_3d_cca))       photol_config%l_3d_cca    = l_3d_cca
IF (PRESENT(l_environ_jo2))  photol_config%l_environ_jo2  = l_environ_jo2
IF (PRESENT(l_environ_jo2b)) photol_config%l_environ_jo2b = l_environ_jo2b

IF (PRESENT(l_environ_ztop)) photol_config%l_environ_ztop = l_environ_ztop

IF (PRESENT(l_enable_diag_um))  photol_config%l_enable_diag_um                 &
                                                      = l_enable_diag_um
IF (PRESENT(l_strat_chem))   photol_config%l_strat_chem   = l_strat_chem

IF (PRESENT(fastjx_prescutoff)) photol_config%fastjx_prescutoff                &
                                                      = fastjx_prescutoff
IF (PRESENT(timestep))       photol_config%timestep   = timestep

IF (PRESENT(pi)) THEN
  const_pi = pi
  const_pi_over_180      = const_pi / 180.0
  const_recip_pi_over_180 = 180.0 / const_pi
  const_s2r             = (2.0*const_pi)/ 86400.0
END IF

IF (PRESENT(o3_mmr_vmr)) const_o3_mmr_vmr = o3_mmr_vmr
IF (PRESENT(molemass_sulp)) const_molemass_sulp = molemass_sulp
IF (PRESENT(molemass_nh42so4)) const_molemass_nh42so4 = molemass_nh42so4
IF (PRESENT(molemass_air)) const_molemass_air = molemass_air
IF (PRESENT(planet_radius)) const_planet_radius = planet_radius

! Set flag to show that a valid photolysis configuration is set up
l_photol_config_available = .TRUE.

! Routine to set up list of driving fields required based on user choices
CALL photol_init_environ_req(error_code_ptr, error_message=error_message,      &
                             error_routine=error_routine)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE photol_setup

END MODULE photol_setup_mod
