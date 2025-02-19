! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module containing top-level subroutine for output of non-ASAD
!   diagnostics
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

MODULE ukca_diags_output_ctl_mod

IMPLICIT NONE

PRIVATE

PUBLIC ukca_diags_output_ctl

CHARACTER(LEN=*), PARAMETER :: ModuleName='ASAD_DIAGS_OUTPUT_CTL_MOD'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE ukca_diags_output_ctl(error_code_ptr, row_length, rows,             &
                                 model_levels, n_use_tracers,                  &
                                 z_top_of_model, do_chemistry,                 &
                                 p_tropopause, p_layer_boundaries,             &
                                 p_theta_levels, plumeria_height, all_tracers, &
                                 photol_rates, diagnostics,                    &
                                 error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Output diagnostics supported by UKCA's diagnostic handling system
!   that are available on the current timestep with the exception of
!   ASAD framework chemical flux diagnostics.
!   Note: While it is generally recommended that output of any supported
!   non-ASAD diagnostic is handled here, the possibility of handling such
!   diagnostics by calls to 'update_diagnostics_2d_real' or
!   'update_diagnostics_3d_real' elsewhere (e.g. at the point of
!   generation) is not precluded.
! ----------------------------------------------------------------------

USE ukca_diagnostics_type_mod, ONLY: diagnostics_type, n_diag_group,           &
                                     dgroup_flat_real, dgroup_fullht_real
USE ukca_diagnostics_init_mod, ONLY: available_diag_varnames,                  &
                                     l_available_diag_chem_timestep,           &
                                     available_diag_asad_ids
USE ukca_diagnostics_output_mod, ONLY: seek_active_requests,                   &
                                       update_diagnostics_2d_real,             &
                                       update_diagnostics_3d_real

USE ukca_fieldname_mod, ONLY: maxlen_diagname,                                 &
                              diagname_jrate_no2, diagname_jrate_o3a,          &
                              diagname_jrate_o3b, diagname_jrate_o2b,          &
                              diagname_p_tropopause, diagname_o3_column_du,    &
                              diagname_plumeria_height

USE ukca_environment_diags_mod, ONLY: photol_rate_diag
USE ukca_config_specification_mod, ONLY:  calc_ozonecol
USE ukca_cspecies, ONLY: n_o3
USE ukca_constants, ONLY: c_o3, dobson
USE ukca_missing_data_mod, ONLY: imdi

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine arguments

! Return code for status reporting
INTEGER, POINTER, INTENT(IN) :: error_code_ptr

! Array dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: n_use_tracers

! height at top of model (m)
REAL, INTENT(IN) :: z_top_of_model

! Chemistry time step indicator
LOGICAL, INTENT(IN) :: do_chemistry

! PV-theta tropopause surface (Pa)
REAL, INTENT(IN) :: p_tropopause(row_length, rows)
! Interface pressures
REAL, INTENT(IN) :: p_layer_boundaries(row_length, rows, 0:model_levels)
! Pressure on theta levels
REAL, INTENT(IN) :: p_theta_levels(row_length, rows, model_levels)
! Plume height from exp eruptions (m)
REAL, INTENT(IN) :: plumeria_height(row_length, rows)
! Tracers
REAL, INTENT(IN) :: all_tracers(row_length, rows, model_levels, n_use_tracers)
! Photolysis rates (may only be allocated on chemistry time steps)
REAL, ALLOCATABLE, INTENT(IN) :: photol_rates(:,:,:,:)

! Diagnostic request info and pointers to parent arrays for diagnostic output
TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics

! Optional arguments for status reporting
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: i
LOGICAL :: l_process_here
LOGICAL :: l_active_requests
CHARACTER(LEN=maxlen_diagname) :: diagname

INTEGER :: i_diag_req(n_diag_group)      ! Indices of diagnostic requests
LOGICAL :: l_check_group(n_diag_group)   ! True to check group for active
                                         ! requests

! Ozone column (from model level to top of atmosphere) in Dobson Unit
REAL :: o3col_du(row_length, rows, model_levels)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DIAGS_OUTPUT_CTL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Process all available diagnostics that are applicable
DO i = 1, SIZE(available_diag_varnames)

  ! Determine whether we need to process the diagnostic:
  ! Ignore diagnostics only available on chemistry time steps if not in a
  ! chemistry time step and always ignore ASAD diagnostics
  IF (l_available_diag_chem_timestep(i)) THEN
    l_process_here = (do_chemistry .AND. available_diag_asad_ids(i) == imdi)
  ELSE
    l_process_here = .TRUE.
  END IF

  IF (l_process_here) THEN

    ! Do diagnostic-specific output processing

    ! In active request checks, avoid checking any request group implicitly
    l_check_group(:) = .FALSE.

    diagname = available_diag_varnames(i)
    SELECT CASE(diagname)

    CASE (diagname_jrate_no2, diagname_jrate_o3a, diagname_jrate_o3b,          &
          diagname_jrate_o2b)

      ! ----------------------------------------------------------------
      ! Photolysis rate from environmental driver data supplied to UKCA
      ! ----------------------------------------------------------------

      ! Check for active requests first to avoid unnecessary processing
      l_check_group(dgroup_flat_real) = .TRUE.
      l_check_group(dgroup_fullht_real) = .TRUE.
      CALL seek_active_requests(diagname, diagnostics, l_check_group,          &
                                i_diag_req, l_active_requests)

      IF (l_active_requests) THEN
        CALL photol_rate_diag(error_code_ptr, diagname, i_diag_req,            &
                              photol_rates, diagnostics,                       &
                              error_message=error_message,                     &
                              error_routine=error_routine)
      END IF

    CASE (diagname_p_tropopause)

      ! ----------------------------------------------------------------
      ! Tropopause pressure
      ! ----------------------------------------------------------------
      CALL update_diagnostics_2d_real(error_code_ptr, diagname,                &
                                      p_tropopause, diagnostics,               &
                                      error_message=error_message,             &
                                      error_routine=error_routine)

    CASE (diagname_o3_column_du)

      ! ----------------------------------------------------------------
      ! Ozone column in Dobson units
      ! ----------------------------------------------------------------

      ! Check for active requests first to avoid unnecessary calculations
      l_check_group(dgroup_flat_real) = .TRUE.
      l_check_group(dgroup_fullht_real) = .TRUE.
      CALL seek_active_requests(diagname, diagnostics, l_check_group,          &
                                i_diag_req, l_active_requests)

      IF (l_active_requests) THEN

        ! Calculate ozone column from actual post-chemistry
        ! ozone, for diagnostic purposes. As ozone is tracer,
        ! this makes the ozone column available on all timesteps
        CALL calc_ozonecol(error_code_ptr, row_length, rows, model_levels,     &
                           z_top_of_model, p_layer_boundaries, p_theta_levels, &
                           all_tracers(:,:,:,n_o3) / c_o3, o3col_du(:,:,:),    &
                           error_message=error_message,                        &
                           error_routine=error_routine)

        IF ( error_code_ptr <= 0 ) THEN
          ! Convert to Dobson Units from molecules/cm2
          o3col_du(:,:,:) = o3col_du(:,:,:)/(dobson * 1.0e-4)

          CALL update_diagnostics_3d_real(error_code_ptr, diagname,            &
                                        o3col_du, diagnostics,                 &
                                        i_diag_req=i_diag_req,                 &
                                        error_message=error_message,           &
                                        error_routine=error_routine)
        END IF

      END IF

    CASE (diagname_plumeria_height)

      CALL update_diagnostics_2d_real(error_code_ptr, diagname,                &
                                      plumeria_height, diagnostics,            &
                                      error_message=error_message,             &
                                      error_routine=error_routine)

    END SELECT

  END IF  ! l_process_here

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_diags_output_ctl

END MODULE ukca_diags_output_ctl_mod
