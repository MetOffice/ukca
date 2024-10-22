! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module for handling the UKCA time step
!
!   The module provides the following procedure for the UKCA API.
!
!     ukca_step - Perform one UKCA time step (overloaded for domains
!                 with different domain dimensions)
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds,
! University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_step_mod

USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_environment_fields_mod, ONLY: environ_field_ptrs
USE ukca_diagnostics_type_mod, ONLY: diagnostics_type,                         &
                                     dgroup_flat_real, dgroup_fullht_real
USE ukca_diagnostics_set_ptrs_mod, ONLY: set_diag_ptrs_1d_domain,              &
                                         set_diag_ptrs_3d_domain
USE ukca_diagnostics_output_mod, ONLY: output_diag_status,                     &
                                       init_diag_status, diag_status_dealloc
USE ukca_tracers_mod, ONLY: tracer_copy_in_1d, tracer_copy_in_3d,              &
                            tracer_copy_out_1d, tracer_copy_out_3d,            &
                            tracer_dealloc, all_tracers
USE ukca_ntp_mod, ONLY: ntp_copy_in_1d, ntp_copy_in_3d, ntp_copy_out_1d,       &
                        ntp_copy_out_3d, ntp_dealloc, all_ntp
USE ukca_main1_mod, ONLY: ukca_main1
USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname

! Dr Hook modules
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC ukca_step

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_STEP_MOD'

! Generic interface for UKCA time step subroutine - overloaded according to
! the dimension of the tracer and NTP data
INTERFACE ukca_step
  MODULE PROCEDURE ukca_step_1d_domain
  MODULE PROCEDURE ukca_step_3d_domain
END INTERFACE ukca_step


CONTAINS


! ----------------------------------------------------------------------
SUBROUTINE ukca_step_1d_domain(timestep_number, current_time,                  &
                               tracer_data_parent, ntp_data_parent,            &
                               r_theta_levels, r_rho_levels,                   &
                               error_code, previous_time, eta_theta_levels,    &
                               diag_status_flat_real, diag_status_fullht_real, &
                               diag_data_flat_real, diag_data_fullht_real,     &
                               error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Variant of the UKCA API generic procedure ukca_step.
!   Performs one UKCA time step.
!   Each input tracer and NTP field is defined for a single column.
!
! Method:
!   1) Set up pointers to access the current diagnostic request
!      info and their associated parent-supplied output arrays and
!      initialise diagnostic status flags from the request status flags.
!   2) Copy tracer and NTP data from the parent-supplied 4D arrays to
!      UKCA's internal data structures, ignoring any data outside the
!      required domain.
!   3) Perform the UKCA time step, copying any diagnostic data to the
!      output arrays and updating diagnostic status flags.
!   4) Copy the updated tracers and NTPs back to the parent arrays.
!   5) Update parent-supplied diagnostic status flag arrays if present.
!   6) Deallocate internal workspace.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

! Model timestep number (counted from basis time at start of run)
INTEGER, INTENT(IN) :: timestep_number

! Current model time (year, month, day, hour, minute, second, day of year)
INTEGER, INTENT(IN) :: current_time(7)

! Height of theta and rho levels from Earth centre
REAL, INTENT(IN) :: r_theta_levels(:,:,0:), r_rho_levels(:,:,:)

! UKCA tracers from the parent model. Dimensions: Z,N
! where Z is no. of levels in tracer fields
!       N is number of tracers
REAL, ALLOCATABLE, INTENT(IN OUT) :: tracer_data_parent(:,:)

! Non-transported prognostics from the parent model. Dimensions: Z,N
REAL, ALLOCATABLE, INTENT(IN OUT) :: ntp_data_parent(:,:)

! Error code for status reporting
INTEGER, TARGET, INTENT(OUT) :: error_code

! Model time at previous timestep (required for chemistry)
INTEGER, OPTIONAL, INTENT(IN) :: previous_time(7)

! Non-dimensional coordinate vector for theta levels (0.0 at planet radius,
! 1.0 at top of model), used to define level height without orography effect.
! Allocatable to preserve bounds (may or may not include Level 0).
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: eta_theta_levels(:)

! Diagnostic status flags
INTEGER, OPTIONAL, INTENT(OUT) :: diag_status_flat_real(:)
INTEGER, OPTIONAL, INTENT(OUT) :: diag_status_fullht_real(:)

! Diagnostic data
REAL, TARGET, OPTIONAL, INTENT(OUT) :: diag_data_flat_real(:)
REAL, TARGET, OPTIONAL, INTENT(OUT) :: diag_data_fullht_real(:,:)

! Further arguments for status reporting
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER, POINTER :: error_code_ptr

TYPE(diagnostics_type) :: diagnostics

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_STEP_1D_DOMAIN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Use parent supplied argument for return code.
! Note that this argument is redundant if UKCA is configured to abort on error
! and may be made optional in future, being replaced by an internal error code
! variable when absent.
error_code_ptr => error_code

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Set up pointers in the 'diagnostics' structure to access the diagnostic
! request informations and the output arrays passed by the parent.
CALL set_diag_ptrs_1d_domain(error_code_ptr, ukca_config, diagnostics,         &
                             data_flat_real=diag_data_flat_real,               &
                             data_fullht_real=diag_data_fullht_real,           &
                             error_message=error_message,                      &
                             error_routine=error_routine)

IF (error_code_ptr <= 0) THEN

  ! Allocate diagnostic status flags and initialise from the request flags.
  CALL init_diag_status(diagnostics)

  ! Populate the all_tracers array using the corresponding tracer data from
  ! the parent model
  CALL tracer_copy_in_1d(error_code_ptr, tracer_data_parent,                   &
                         ukca_config%model_levels,                             &
                         error_message=error_message,                          &
                         error_routine=error_routine)

END IF

IF (error_code_ptr <= 0) THEN

  ! Populate the all_ntp array using the corresponding non-transported
  ! prognostic data from the parent model
  CALL ntp_copy_in_1d(error_code_ptr, ntp_data_parent,                         &
                      ukca_config%model_levels, error_message=error_message,   &
                      error_routine=error_routine)

END IF

! Do the time step
IF (error_code_ptr <= 0) THEN

  CALL ukca_main1(error_code_ptr, timestep_number, current_time,               &
                  environ_field_ptrs, r_theta_levels, r_rho_levels,            &
                  diagnostics, all_tracers, all_ntp,                           &
                  previous_time=previous_time,                                 &
                  eta_theta_levels=eta_theta_levels,                           &
                  error_message=error_message, error_routine=error_routine)

END IF

IF (error_code_ptr <= 0) THEN

  ! Update the tracer_data_parent and ntp_data_parent arrays from the UKCA state
  ! at the end of the time step. These are then passed back to the parent.
  CALL tracer_copy_out_1d(ukca_config%model_levels, tracer_data_parent)
  CALL ntp_copy_out_1d(ukca_config%model_levels, ntp_data_parent)

  ! Output parent copy of diagnostics status flags if required (flat real group)
  IF (PRESENT(diag_status_flat_real)) THEN
    CALL output_diag_status(error_code_ptr, dgroup_flat_real, diagnostics,     &
                            status_flags=diag_status_flat_real,                &
                            error_message=error_message,                       &
                            error_routine=error_routine)
  END IF

END IF

IF (error_code_ptr <= 0) THEN

  ! Output parent copy of diagnostic status flags if required
  ! (full height real group)
  IF (PRESENT(diag_status_fullht_real)) THEN
    CALL output_diag_status(error_code_ptr, dgroup_fullht_real, diagnostics,   &
                            status_flags=diag_status_fullht_real,              &
                            error_message=error_message,                       &
                            error_routine=error_routine)
  END IF

END IF

! Clear internal UKCA state data and diagnostic flags ready for next time step
! (Do this irrespective of error code)
CALL ntp_dealloc()
CALL tracer_dealloc()
CALL diag_status_dealloc(diagnostics)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_step_1d_domain

! ----------------------------------------------------------------------
SUBROUTINE ukca_step_3d_domain(timestep_number, current_time,                  &
                               tracer_data_parent, ntp_data_parent,            &
                               r_theta_levels, r_rho_levels,                   &
                               error_code, previous_time, eta_theta_levels,    &
                               diag_status_flat_real, diag_status_fullht_real, &
                               diag_data_flat_real, diag_data_fullht_real,     &
                               error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Variant of the UKCA API generic procedure ukca_step.
!   Performs one UKCA time step.
!   Each input tracer and NTP field is defined on a 3D grid.
!
! Method:
!   See ukca_step_1d_domain above
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

! Model timestep number (counted from basis time at start of run)
INTEGER, INTENT(IN) :: timestep_number

! Current model time (year, month, day, hour, minute, second, day of year)
INTEGER, INTENT(IN) :: current_time(7)

REAL, INTENT(IN) :: r_theta_levels(:,:,0:), r_rho_levels(:,:,:)

! UKCA tracers from the parent model. Dimensions: X,Y,Z,N
! where X is row length of tracer field (= no. of columns)
!       Y is no. of rows in tracer field
!       Z is no. of levels in tracer fields
!       N is number of tracers
REAL, ALLOCATABLE, INTENT(IN OUT) :: tracer_data_parent(:,:,:,:)

! Non-transported prognostics from the parent model. Dimensions: X,Y,Z,N
REAL, ALLOCATABLE, INTENT(IN OUT) :: ntp_data_parent(:,:,:,:)

! Error code for status reporting
INTEGER, TARGET, INTENT(OUT) :: error_code

! Model time at previous timestep (required for chemistry)
INTEGER, OPTIONAL, INTENT(IN) :: previous_time(7)

! Non-dimensional coordinate vector for theta levels (0.0 at planet radius,
! 1.0 at top of model), used to define level height without orography effect.
! Allocatable to preserve bounds (may or may not include Level 0).
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: eta_theta_levels(:)

! Diagnostic status flags
INTEGER, OPTIONAL, INTENT(OUT) :: diag_status_flat_real(:)
INTEGER, OPTIONAL, INTENT(OUT) :: diag_status_fullht_real(:)

! Diagnostic data
REAL, TARGET, OPTIONAL, INTENT(OUT) :: diag_data_flat_real(:,:,:)
REAL, TARGET, OPTIONAL, INTENT(OUT) :: diag_data_fullht_real(:,:,:,:)

! Further arguments for status reporting
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER, POINTER :: error_code_ptr

TYPE(diagnostics_type) :: diagnostics

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_STEP_3D_DOMAIN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Use parent supplied argument for return code.
! Note that this argument is redundant if UKCA is configured to abort on error
! and may be made optional in future, being replaced by an internal error code
! variable when absent
error_code_ptr => error_code

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Set up pointers in the 'diagnostics' structure to access the diagnostic
! request informations and the output arrays passed by the parent.
CALL set_diag_ptrs_3d_domain(error_code_ptr, ukca_config, diagnostics,         &
                             data_flat_real=diag_data_flat_real,               &
                             data_fullht_real=diag_data_fullht_real,           &
                             error_message=error_message,                      &
                             error_routine=error_routine)

IF (error_code_ptr <= 0) THEN

  ! Allocate diagnostic status flags and initialise from the request flags.
  CALL init_diag_status(diagnostics)

  ! Populate the all_tracers array using the corresponding tracer data from
  ! the parent model
  CALL tracer_copy_in_3d(error_code_ptr, tracer_data_parent,                   &
                         ukca_config%row_length, ukca_config%rows,             &
                         ukca_config%model_levels,                             &
                         error_message=error_message,                          &
                         error_routine=error_routine)

END IF

IF (error_code_ptr <= 0) THEN

  ! Populate the all_ntp array using the corresponding non-transported
  ! prognostic data from the parent model
  CALL ntp_copy_in_3d(error_code_ptr, ntp_data_parent, ukca_config%row_length, &
                      ukca_config%rows, ukca_config%model_levels,              &
                      error_message=error_message, error_routine=error_routine)

END IF

! Do the time step
IF (error_code_ptr <= 0) THEN

  CALL ukca_main1(error_code_ptr, timestep_number, current_time,               &
                  environ_field_ptrs, r_theta_levels, r_rho_levels,            &
                  diagnostics, all_tracers, all_ntp,                           &
                  previous_time=previous_time,                                 &
                  eta_theta_levels=eta_theta_levels,                           &
                  error_message=error_message, error_routine=error_routine)

END IF

IF (error_code_ptr <= 0) THEN

  ! Update the tracer_data_parent and ntp_data_parent arrays from the UKCA state
  ! at the end of the time step. These are then passed back to the parent.
  CALL tracer_copy_out_3d(ukca_config%row_length, ukca_config%rows,            &
                          ukca_config%model_levels, tracer_data_parent)
  CALL ntp_copy_out_3d(ukca_config%row_length, ukca_config%rows,               &
                       ukca_config%model_levels, ntp_data_parent)

  ! Output parent copy of diagnostics status flags if required (flat real group)
  IF (PRESENT(diag_status_flat_real)) THEN
    CALL output_diag_status(error_code_ptr, dgroup_flat_real, diagnostics,     &
                            status_flags=diag_status_flat_real,                &
                            error_message=error_message,                       &
                            error_routine=error_routine)
  END IF

END IF

IF (error_code_ptr <= 0) THEN

  ! Output parent copy of diagnostic status flags if required
  ! (full height real group)
  IF (PRESENT(diag_status_fullht_real)) THEN
    CALL output_diag_status(error_code_ptr, dgroup_fullht_real, diagnostics,   &
                            status_flags=diag_status_fullht_real,              &
                            error_message=error_message,                       &
                            error_routine=error_routine)
  END IF

END IF

! Clear internal UKCA state data and diagnostic flags ready for next time step
! (Do this irrespective of error code)
CALL ntp_dealloc()
CALL tracer_dealloc()
CALL diag_status_dealloc(diagnostics)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_step_3d_domain

END MODULE ukca_step_mod
