! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module for handling the Photolysis environmental driver fields.
!   Environmental drivers are input fields defined on the current
!   model grid that may be varied by the parent application during the run.
!
!   The module provides the following procedures for the Photolysis API.
!
!     photol_get_environ_varlist - Returns list of names of required fields.
!
!   The following additional public procedures are provided for use within
!   photolysis
!
!     photol_init_environ_req        - Determines environment data requirement
!
!     photol_check_environ_availability - Checks availability of required
!                                      environment fields.
!     photol_print_environ_list      - Writes a summary of the current
!                                      environment data to the log file.
!     photol_clear_environ_req       - Resets all data relating to environment
!                                      data requirement to its initial state
!                                      for a new configuration.
!
!   The module also provides a public logical 'l_environ_req_available'
!   that indicates the availability status of the environment data requirement
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds,
! University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA/Photolysis
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE photol_environment_mod

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

USE photol_fieldname_mod,  ONLY: fieldname_len,                                &
                                 fldname_aod_sulph_aitk,                       &
                                 fldname_aod_sulph_accum,                      &
                                 fldname_area_cloud_fraction,                  &
                                 fldname_conv_cloud_amount,                    &
                                 fldname_conv_cloud_base,                      &
                                 fldname_conv_cloud_lwp,                       &
                                 fldname_conv_cloud_top,                       &
                                 fldname_cos_latitude,                         &
                                 fldname_equation_of_time,                     &
                                 fldname_land_fraction,                        &
                                 fldname_longitude,                            &
                                 fldname_ozone_mmr,                            &
                                 fldname_p_layer_boundaries,                   &
                                 fldname_p_theta_levels,                       &
                                 fldname_photol_rates_2d,                      &
                                 fldname_qcf,                                  &
                                 fldname_qcl,                                  &
                                 fldname_rad_ctl_jo2,                          &
                                 fldname_rad_ctl_jo2b,                         &
                                 fldname_r_rho_levels,                         &
                                 fldname_r_theta_levels,                       &
                                 fldname_sec_since_midnight,                   &
                                 fldname_sin_declination,                      &
                                 fldname_sin_latitude,                         &
                                 fldname_so4_aitken,                           &
                                 fldname_so4_accum,                            &
                                 fldname_surf_albedo,                          &
                                 fldname_t_theta_levels,                       &
                                 fldname_tan_latitude,                         &
                                 fldname_z_top_of_model

USE ukca_error_mod,          ONLY: maxlen_message, maxlen_procname,            &
                                   errcode_env_req_uninit, error_report,       &
                                   errcode_value_invalid

IMPLICIT NONE

PRIVATE

! Public procedures
PUBLIC photol_init_environ_req, photol_get_environ_varlist,                    &
       photol_print_environ_list, photol_clear_environ_req

! Public flag for use within Photolysis to indicate whether the environment
! fields requirement has been initialised
LOGICAL, SAVE, PUBLIC :: l_photol_environ_req_available = .FALSE.

! Maximum number of envrionment fields in current photolysis schemes
INTEGER, PARAMETER :: n_max_fields = 30

! Actual number of env fields required for this configuration
INTEGER, SAVE, PUBLIC :: n_fields

! Field group codes for collating driver fields in arrays, as required for
! the API routine 'photol_step_control' :
! 'scalar' group comprises fields defined by a scalar value internally
! 'flat' groups comprise fields defined on a flat spatial grid with 2D
! representation internally (input can have reduced dimension)
! 'fullht' group comprises fields defined on a full height spatial grid with 3D
! representation internally (input can have reduced dimension)
! 'fullht0' group comprises fields defined on an extended full height spatial
! grid with an additional level 0 below and 3D representation internally
! (input can have reduced dimension)

INTEGER, PARAMETER :: group_undefined    = 0      ! Not assigned to a group
INTEGER, PARAMETER :: group_scalar_real  = 1      ! Scalar real
INTEGER, PARAMETER :: group_flat_integer = 2      ! 2D spatial integer
INTEGER, PARAMETER :: group_flat_real    = 3      ! 2D spatial real
INTEGER, PARAMETER :: group_fullht_real  = 4      ! 3D spatial real on all
                                                  ! model levels (theta levels)
INTEGER, PARAMETER :: group_fullht0_real = 5      ! 3D spatial real on all model
                                                  ! levels + zero level below
INTEGER, PARAMETER :: group_fullhtphot_real = 6   ! 3D spatial real on all model
                                                  ! levels + photolytic species
! Any environmental driver not assigned to a group will be ignored by
! API routine 'photol_step_control' but will appear in the full list returned by
! 'photol_get_environ_varlist'.

! List of environment fields required, by subgroup
CHARACTER(LEN=fieldname_len), ALLOCATABLE, TARGET, SAVE, PUBLIC ::             &
      environ_fldnames_scalar_real(:)      ! Field names of real scalars
CHARACTER(LEN=fieldname_len), ALLOCATABLE, TARGET, SAVE, PUBLIC ::             &
      environ_fldnames_flat_integer(:)     ! Field names of 2D spatial
                                           ! integers
CHARACTER(LEN=fieldname_len), ALLOCATABLE, TARGET, SAVE, PUBLIC ::             &
      environ_fldnames_flat_real(:)        ! Field names of 2D spatial reals
CHARACTER(LEN=fieldname_len), ALLOCATABLE, TARGET, SAVE, PUBLIC ::             &
      environ_fldnames_fullht_real(:)      ! Field names of full height 3D
                                           ! spatial reals
CHARACTER(LEN=fieldname_len), ALLOCATABLE, TARGET, SAVE, PUBLIC ::             &
      environ_fldnames_fullht0_real(:)     ! Field names of 3D spatial reals
                                           ! on all model levels plus zero
                                           ! level below
CHARACTER(LEN=fieldname_len), ALLOCATABLE, TARGET, SAVE, PUBLIC ::             &
      environ_fldnames_fullhtphot_real(:)  ! Field names of 3D spatial reals
                                           ! on all model levels plus number   &
                                           ! of photolytic species
! Data structure holding information on photolysis environment fields.
TYPE, PUBLIC :: photol_envfld_type
  CHARACTER(LEN=fieldname_len) :: env_fieldname    ! Env field name
  INTEGER                      :: fld_group_type   ! Field group type
END TYPE photol_envfld_type

TYPE(photol_envfld_type), PUBLIC, SAVE :: photol_env_fields(n_max_fields)

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName='PHOTOL_ENVIRONMENT_MOD'

CONTAINS

  ! ----------------------------------------------------------------------
  ! Routine which sets up list of photol environment driver fields that
  ! are needed for this configuration
  !
  ! ----------------------------------------------------------------------
SUBROUTINE photol_init_environ_req(error_code_ptr, error_message,              &
                                   error_routine)

USE photol_config_specification_mod,  ONLY: photol_config,                     &
                                            i_scheme_nophot,                   &
                                            i_scheme_photol_strat,             &
                                            i_scheme_phot2d,                   &
                                            i_scheme_fastjx

IMPLICIT NONE

! Arguments and local variables for error handling
INTEGER, POINTER, INTENT(IN) :: error_code_ptr
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

CHARACTER(LEN=maxlen_message) :: cmessage

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'PHOTOL_INIT_ENVIRON_REQ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in,                &
                        zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Clear environment requirements list if set previously
IF (l_photol_environ_req_available)  CALL photol_clear_environ_req()

n_fields = 0
photol_env_fields(:)%env_fieldname = 'empty'
photol_env_fields(:)%fld_group_type = group_undefined
! Determine requirement of driving fields based on scheme choices

! If no photolysis scheme chosen - no environment fields required
IF ( photol_config%i_photol_scheme == i_scheme_nophot ) THEN
  l_photol_environ_req_available = .TRUE.

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,             &
                          zhook_handle)
  RETURN
END IF

! Check that n_fields does not exceed n_max_fields before adding to list
!
!-- Fields that are required for all Photolysis schemes/ combinations

! Sin declination angle
IF ( n_fields <= n_max_fields ) THEN
  CALL register_env_field(n_fields, fldname_sin_declination,                   &
                          group_scalar_real)
END IF
! SIN of latitude
IF ( n_fields <= n_max_fields ) THEN
  CALL register_env_field(n_fields, fldname_sin_latitude,                      &
                          group_flat_real)
END IF
! TAN of latitude
IF ( n_fields <= n_max_fields ) THEN
  CALL register_env_field(n_fields,  fldname_tan_latitude,                     &
                          group_flat_real)
END IF
!-- Fields required for a single scheme or option
! JO2 rates from Radiation
IF ( n_fields <= n_max_fields .AND. photol_config%l_environ_jo2 ) THEN
  CALL register_env_field(n_fields,  fldname_rad_ctl_jo2,                      &
                          group_fullht_real)
END IF
! JO2B rates from Radiation
IF ( n_fields <= n_max_fields .AND. photol_config%l_environ_jo2b ) THEN
  CALL register_env_field(n_fields,  fldname_rad_ctl_jo2b,                     &
                          group_fullht_real)
END IF

!-- Fields required only for 2-D Photolysis
IF ( photol_config%i_photol_scheme == i_scheme_phot2d ) THEN
  ! Photolysis rates (from file)
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_photol_rates_2d,                 &
                            group_fullhtphot_real)
  END IF
END IF

!-- Fields required only for Stratospheric schemes - by photol_solang
IF ( photol_config%l_strat_chem ) THEN
  ! Equation of time
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_equation_of_time,                &
                            group_scalar_real)
  END IF
  ! Seconds since midnight
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_sec_since_midnight,              &
                            group_scalar_real)
  END IF
  ! COS latitude
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_cos_latitude,                    &
                            group_flat_real)
  END IF
END IF   ! l_strat_chem

! -- Fields required only for FastJX -- by type
IF ( photol_config%i_photol_scheme == i_scheme_fastjx ) THEN
  ! Conv Cloud Base
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_conv_cloud_base,                 &
                            group_flat_integer)
  END IF
  ! Conv Cloud Top
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_conv_cloud_top,                  &
                            group_flat_integer)
  END IF
  ! Land Fraction
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_land_fraction,                   &
                            group_flat_real)
  END IF
  ! Surface albedo
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_surf_albedo,                     &
                            group_flat_real)
  END IF
  ! Aerosol Optical Depth - Sulphate, Accumulation mode
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_aod_sulph_accum,                 &
                            group_fullht_real)
  END IF
  ! Aerosol Optical Depth - Sulphate, Aitken mode
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_aod_sulph_aitk,                  &
                            group_fullht_real)
  END IF
  ! Area Cloud fraction
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_area_cloud_fraction,             &
                            group_fullht_real)
  END IF
  ! Frozen Cloud fraction
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_qcf, group_fullht_real)
  END IF
  ! Liquid Cloud fraction
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_qcl, group_fullht_real)
  END IF
  ! Height of Rho levels (from earth's centre)
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_r_rho_levels,                    &
                            group_fullht_real)
  END IF
  ! Sulphate mmr - accumulation mode
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_so4_accum,                       &
                            group_fullht_real)
  END IF
  ! Sulphate mmr - aitken mode
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_so4_aitken,                      &
                            group_fullht_real)
  END IF

  !-- Fields only required if FastJX and not using PC2 cloud scheme
  IF ( .NOT. photol_config%l_cloud_pc2 ) THEN
    ! Conv Cloud Amount
    IF ( n_fields <= n_max_fields ) THEN
      CALL register_env_field(n_fields, fldname_conv_cloud_amount,             &
                              group_fullht_real)
    END IF
    ! Conv Cloud liquid water path
    IF ( n_fields <= n_max_fields ) THEN
      CALL register_env_field(n_fields, fldname_conv_cloud_lwp,                &
                              group_flat_real)
    END IF
  END IF

END IF   ! scheme= fastjx

!-- Fields required if either FastJX or a Stratospheric scheme is used.
IF ( photol_config%i_photol_scheme == i_scheme_fastjx  .OR.                    &
      photol_config%l_strat_chem ) THEN
  ! z_top_of_model, if being supplied by parent
  IF ( n_fields <= n_max_fields .AND. photol_config%l_environ_ztop ) THEN
    CALL register_env_field(n_fields, fldname_z_top_of_model, group_scalar_real)
  END IF
  ! Longitude
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_longitude, group_flat_real)
  END IF
  ! Ozone Mass mixing ratio
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_ozone_mmr,                       &
                            group_fullht_real)
  END IF
  ! Temperature on theta levels
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_t_theta_levels,                  &
                            group_fullht_real)
  END IF
  ! Pressure at layer boundaries
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_p_layer_boundaries,              &
                            group_fullht0_real)
  END IF
END IF   ! Fastjx or l_strat_chem

! Height of theta levels will be required if using Fast-JX, or to calculate
! model top for stratospheric schemes if latter is not supplied by parent
IF ( photol_config%i_photol_scheme == i_scheme_fastjx  .OR.                    &
  ( photol_config%l_strat_chem .AND. .NOT. photol_config%l_environ_ztop) ) THEN
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_r_theta_levels,                  &
                            group_fullht0_real)
  END IF
END IF

!-- Fields required for 2-D, FastJX or Stratospheric scheme
IF ( photol_config%i_photol_scheme == i_scheme_fastjx  .OR.                    &
     photol_config%i_photol_scheme == i_scheme_phot2d  .OR.                    &
     photol_config%l_strat_chem ) THEN
  ! Pressure on theta levels
  IF ( n_fields <= n_max_fields ) THEN
    CALL register_env_field(n_fields, fldname_p_theta_levels,                  &
                            group_fullht_real)
  END IF
END IF

IF ( n_fields > n_max_fields ) THEN
  error_code_ptr = errcode_value_invalid
  cmessage = 'Number of required environment fields exceeds maximum'
  CALL error_report(photol_config%i_error_method, error_code_ptr, cmessage,    &
         RoutineName, msg_out=error_message, locn_out=error_routine)

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,             &
                          zhook_handle)
  RETURN
END IF

! Flag to denote environment requirements list has been set
l_photol_environ_req_available = .TRUE.

! Added for now to debug
CALL photol_print_environ_list()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,               &
                        zhook_handle)
RETURN

END SUBROUTINE photol_init_environ_req

! ---------------------------------------------------------------------------
! Routine that returns the list of environment fields required for this
! configuration - grouped by field type
! ---------------------------------------------------------------------------
SUBROUTINE photol_get_environ_varlist( errcode, varnames_scalar_real_ptr,      &
                      varnames_flat_integer_ptr, varnames_flat_real_ptr,       &
                      varnames_fullht_real_ptr, varnames_fullht0_real_ptr,     &
                      varnames_fullhtphot_real_ptr,                            &
                      error_message, error_routine )

USE photol_config_specification_mod,  ONLY: photol_config

IMPLICIT NONE

! Optional character arrays for names of each type of environment field
! - defined as pointers since size is not known to calling routine
CHARACTER(LEN=fieldname_len), POINTER, OPTIONAL, INTENT(OUT) ::                &
                              varnames_scalar_real_ptr(:)
CHARACTER(LEN=fieldname_len), POINTER, OPTIONAL, INTENT(OUT) ::                &
                              varnames_flat_integer_ptr(:)
CHARACTER(LEN=fieldname_len), POINTER, OPTIONAL, INTENT(OUT) ::                &
                              varnames_flat_real_ptr(:)
CHARACTER(LEN=fieldname_len), POINTER, OPTIONAL, INTENT(OUT) ::                &
                              varnames_fullht_real_ptr(:)
CHARACTER(LEN=fieldname_len), POINTER, OPTIONAL, INTENT(OUT) ::                &
                              varnames_fullht0_real_ptr(:)
CHARACTER(LEN=fieldname_len), POINTER, OPTIONAL, INTENT(OUT) ::                &
                              varnames_fullhtphot_real_ptr(:)

! Error handling
INTEGER, TARGET, INTENT(IN OUT) :: errcode

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

INTEGER, POINTER :: error_code_ptr
CHARACTER(LEN=maxlen_message) :: err_msg

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'PHOTOL_GET_ENVIRON_VARLIST'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in,                &
                        zhook_handle)

error_code_ptr => errcode

error_code_ptr = 0
err_msg = ''
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check that environment list has been set up, else return error
IF ( .NOT. l_photol_environ_req_available ) THEN

  error_code_ptr = errcode_env_req_uninit
  err_msg = 'Photolysis Environment Field requirements not set'
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
         RoutineName, msg_out=error_message, locn_out=error_routine)

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,             &
                          zhook_handle)
  RETURN
END IF

IF (PRESENT(varnames_scalar_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_fldnames_scalar_real)) THEN
    CALL set_fldgroup_varlist(group_scalar_real,                               &
                              environ_fldnames_scalar_real)
  END IF
  varnames_scalar_real_ptr => environ_fldnames_scalar_real
END IF
IF (PRESENT(varnames_flat_integer_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_fldnames_flat_integer)) THEN
    CALL set_fldgroup_varlist(group_flat_integer,                              &
                              environ_fldnames_flat_integer)
  END IF
  varnames_flat_integer_ptr => environ_fldnames_flat_integer
END IF
IF (PRESENT(varnames_flat_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_fldnames_flat_real)) THEN
    CALL set_fldgroup_varlist(group_flat_real,                                 &
                              environ_fldnames_flat_real)
  END IF
  varnames_flat_real_ptr => environ_fldnames_flat_real
END IF
IF (PRESENT(varnames_fullht_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_fldnames_fullht_real)) THEN
    CALL set_fldgroup_varlist(group_fullht_real,                               &
                              environ_fldnames_fullht_real)
  END IF
  varnames_fullht_real_ptr => environ_fldnames_fullht_real
END IF
IF (PRESENT(varnames_fullht0_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_fldnames_fullht0_real)) THEN
    CALL set_fldgroup_varlist(group_fullht0_real,                              &
                              environ_fldnames_fullht0_real)
  END IF
  varnames_fullht0_real_ptr => environ_fldnames_fullht0_real
END IF

IF (PRESENT(varnames_fullhtphot_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_fldnames_fullhtphot_real)) THEN
    CALL set_fldgroup_varlist(group_fullhtphot_real,                           &
                              environ_fldnames_fullhtphot_real)
  END IF
  varnames_fullhtphot_real_ptr => environ_fldnames_fullhtphot_real
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,               &
                        zhook_handle)
RETURN
END SUBROUTINE photol_get_environ_varlist

! ---------------------------------------------------------------------------
! Routine to print out full list of environment fields required for
! this configuration.
! ---------------------------------------------------------------------------
SUBROUTINE photol_print_environ_list()

USE umPrintMgr, ONLY: umMessage, UmPrint

IMPLICIT NONE

INTEGER :: i

! Array of charcter field types - indexed as per fld_group codes
CHARACTER(LEN=*), PARAMETER :: cgroup_type(6) = ['scalar_real    ',            &
   'flat_integer   ','flat_real      ','fullht_real    ','fullht0_real   ',    &
   'fullhtphot_real']

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHOTOL_PRINT_ENVIRON_LIST'

! Check that environment list has been set up
IF ( .NOT. l_photol_environ_req_available )  RETURN

WRITE(umMessage,'(A)')' %%%%%%% PHOTOL ENVIRONMENT FIELD LIST %%%%%%%%%%% '
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)')' ================================================= '
CALL umPrint(umMessage,src=RoutineName)

DO i = 1, n_fields
  WRITE(umMessage,'(I3,2(1x,A,A))') i, photol_env_fields(i)%env_fieldname,     &
        cgroup_type(photol_env_fields(i)%fld_group_type)
  CALL umPrint(umMessage,src=RoutineName)
END DO
WRITE(umMessage,'(A)')' ================================================= '
CALL umPrint(umMessage,src=RoutineName)

RETURN
END SUBROUTINE photol_print_environ_list

! ---------------------------------------------------------------------------
! Routine to clear the Environment requirements list
! ---------------------------------------------------------------------------
SUBROUTINE photol_clear_environ_req()

IMPLICIT NONE

IF ( l_photol_environ_req_available ) THEN
  n_fields = 0
  photol_env_fields(:)%env_fieldname = 'empty'
  photol_env_fields(:)%fld_group_type = group_undefined

  IF (ALLOCATED(environ_fldnames_scalar_real))                                 &
                DEALLOCATE(environ_fldnames_scalar_real)
  IF (ALLOCATED(environ_fldnames_flat_integer))                                &
                DEALLOCATE(environ_fldnames_flat_integer)
  IF (ALLOCATED(environ_fldnames_flat_real))                                   &
                DEALLOCATE(environ_fldnames_flat_real)
  IF (ALLOCATED(environ_fldnames_fullht_real))                                 &
                DEALLOCATE(environ_fldnames_fullht_real)
  IF (ALLOCATED(environ_fldnames_fullht0_real))                                &
                DEALLOCATE(environ_fldnames_fullht0_real)
  IF (ALLOCATED(environ_fldnames_fullhtphot_real))                             &
                DEALLOCATE(environ_fldnames_fullhtphot_real)

  l_photol_environ_req_available = .FALSE.
END IF

RETURN
END SUBROUTINE photol_clear_environ_req

! Routine to add given fieldname to photol environ fields list and
! increment the index, as long as it is within max_fields
! ---------------------------------------------------------------------------
SUBROUTINE register_env_field(n, fieldname, field_group)

IMPLICIT NONE

INTEGER, INTENT(IN OUT)      :: n           ! index
CHARACTER(LEN=*), INTENT(IN) :: fieldname   ! field to register
INTEGER, INTENT(IN)          :: field_group ! field type

n = n + 1
IF ( n <= n_max_fields ) THEN
  photol_env_fields(n)%env_fieldname = fieldname
  photol_env_fields(n)%fld_group_type = field_group
END IF

RETURN
END SUBROUTINE register_env_field

!----------------------------------------------------------------------------
! Routine to populate the variable names list for given field group
!----------------------------------------------------------------------------
SUBROUTINE set_fldgroup_varlist(fld_group, fldnames)

IMPLICIT NONE

INTEGER, INTENT(IN) :: fld_group
CHARACTER(LEN=fieldname_len), ALLOCATABLE, INTENT(OUT) :: fldnames(:)
INTEGER :: n, i
LOGICAL, ALLOCATABLE :: mask(:)  ! To count fields of input type

n = 0
ALLOCATE(mask(n_fields))
mask = ( photol_env_fields(1:n_fields)%fld_group_type == fld_group )
n = COUNT(mask)
ALLOCATE(fldnames(n))
DEALLOCATE(mask)

i = 0  ! index for fldnames array
DO n = 1, n_fields
  IF ( photol_env_fields(n)%fld_group_type == fld_group ) THEN
    i = i + 1
    fldnames(i) = photol_env_fields(n)%env_fieldname
  END IF
END DO

RETURN
END SUBROUTINE set_fldgroup_varlist

END MODULE photol_environment_mod
