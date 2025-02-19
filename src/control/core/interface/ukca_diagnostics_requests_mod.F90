! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module providing procedures to create and update requests for UKCA
!   diagnostics.
!
!   UKCA API procedures provided:
!
!     ukca_set_diagnostic_requests
!       - Create a new set of diagnostic requests
!     ukca_update_diagnostic_requests
!       - Activate or deactivate existing requests via their status flags
!     ukca_get_diagnostic_request_info
!       - Return details of current diagnostic requests
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

MODULE ukca_diagnostics_requests_mod

USE ukca_diagnostics_type_mod, ONLY: n_diag_group, dgroup_flat_real,           &
                                     dgroup_fullht_real,                       &
                                     diag_requests_type, diagnostics_type,     &
                                     diag_status_inactive,                     &
                                     diag_status_requested,                    &
                                     diag_status_unavailable

USE ukca_diagnostics_master_mod, ONLY: master_diag_list

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_ukca_uninit, errcode_diag_req_unknown,       &
                          errcode_diag_req_duplicate,                          &
                          errcode_diag_req_unsupported_use,                    &
                          errcode_diag_mismatch, errcode_value_unknown,        &
                          errcode_value_invalid, error_report,                 &
                          i_error_method_abort

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

! Public procedures
PUBLIC ukca_set_diagnostic_requests, ukca_update_diagnostic_requests,          &
       ukca_get_diagnostic_request_info

! Current diagnostic requests
TYPE(diag_requests_type), TARGET, PUBLIC :: diag_requests(n_diag_group)

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_DIAGNOSTIC_REQUESTS_MOD'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE ukca_set_diagnostic_requests(error_code,                            &
                                        names_flat_real, names_fullht_real,    &
                                        dreq_status_flat_real,                 &
                                        dreq_status_fullht_real,               &
                                        error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!
!   UKCA API procedure for creating a set of diagnostic requests
!   comprising a list of requests for one or both diagnostic groups.
!   A diagnostic request consists of the diagnostic name and a request
!   status flag which is set to 'diag_status_requested' or
!   'diag_status_inactive' to indicate whether the request is to be
!   active currently. The input status is overridden for any diagnostic
!   not available in the current UKCA configuration, being set to
!   'diag_status_unavailable' on return.
!
!   Names and status flags for each group are supplied as separate arrays.
!   Each array is optional but the pair of arrays for a group must be
!   consistent.
!
!   This routine may be called more than once, with the new diagnostic
!   request set replacing the old, but only if ASAD diagnostics have not
!   been requested previously; repetition of setup processing is not
!   currently supported for these diagnostics.
!
!   Note that a request set can be empty as indicated by omission of
!   arguments and/or use of zero length arrays. The latter are allowed to
!   support use cases where the number of diagnostics in each group is
!   not fixed and may sometimes be zero.
!
!   To use time-varying diagnostic requests that include ASAD diagnostics
!   this routine should be called once to set up requests for all
!   diagnostics to be used during the run and status flags updated
!   subsequently by calls to 'ukca_update_diagnostic_requests' to
!   activate/deactivate individual diagnostic requests as needed.
!
!   N.B. When UKCA is coupled with the Unified Model, this subroutine will
!   handle UM legacy requests for ASAD diagnostics made via the STASH system
!   as well as requests made explicitly via the argument list.
!   It must be called to support the use of ASAD diagnostics in the UM even
!   if there are no explicit requests to set up.
!
! Developer's note:
!   The 'error_code' argument can be made optional since it is not required
!   if UKCA is configured to abort on error. It remains mandatory at
!   present for consistency with other API procedures that do not yet
!   support this functionality.
! ----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: ukca_config, l_ukca_config_available
USE asad_chem_flux_diags, ONLY: asad_setstash_chemdiag,                        &
                                asad_init_chemdiag, n_asad_diags

IMPLICIT NONE

! Subroutine arguments

INTEGER, TARGET, INTENT(OUT) :: error_code

CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: names_flat_real(:)
                                            ! Names of requested diagnostics in
                                            ! flat real group
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: names_fullht_real(:)
                                            ! Names of requested diagnostics in
                                            ! full height real group

INTEGER, OPTIONAL, INTENT(IN OUT) :: dreq_status_flat_real(:)
                                            ! Status flags for requested
                                            ! diagnostics in flat real group
INTEGER, OPTIONAL, INTENT(IN OUT) :: dreq_status_fullht_real(:)
                                            ! Status flags for requested
                                            ! diagnostics in full height real
                                            ! group

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER, POINTER :: error_code_ptr
INTEGER :: group
INTEGER :: n_req_flat_real
INTEGER :: n_req_fullht_real

CHARACTER(LEN=maxlen_message) :: message_txt  ! Buffer for output message

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SET_DIAGNOSTIC_REQUESTS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Use parent supplied argument for return code
error_code_ptr => error_code

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

message_txt = ''

! Check for availability of UKCA configuration data
IF (.NOT. l_ukca_config_available) THEN
  error_code_ptr = errcode_ukca_uninit
  CALL error_report(i_error_method_abort, error_code_ptr,                      &
         'No UKCA configuration has been set up', RoutineName,                 &
         msg_out=error_message, locn_out=error_routine)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Cancel any existing requests for all diagnostic groups
DO group = 1, n_diag_group
  IF (ALLOCATED(diag_requests(group)%varnames))                                &
    DEALLOCATE(diag_requests(group)%varnames)
  IF (ALLOCATED(diag_requests(group)%status_flags))                            &
    DEALLOCATE(diag_requests(group)%status_flags)
  IF (ALLOCATED(diag_requests(group)%i_master))                                &
    DEALLOCATE(diag_requests(group)%i_master)
END DO

n_req_flat_real = 0
n_req_fullht_real = 0

! Check availability and create new request data for each diagnostic group.
! Status flags will be updated for any unavailable diagnostics if necessary.

IF (PRESENT(names_flat_real)) THEN
  IF (PRESENT(dreq_status_flat_real)) THEN
    CALL set_diag_requests(error_code_ptr, dgroup_flat_real,                   &
                           names_flat_real, dreq_status_flat_real,             &
                           diag_requests, n_req_flat_real,                     &
                           error_message=error_message,                        &
                           error_routine=error_routine)
  ELSE
    error_code_ptr = errcode_diag_mismatch
    message_txt =                                                              &
      'Status flags missing for diagnostic requests (flat real group)'
  END IF
ELSE
  IF (PRESENT(dreq_status_flat_real)) THEN
    error_code_ptr = errcode_diag_mismatch
    message_txt =                                                              &
      'Status flag array present but no diagnostic name array ' //             &
      '(flat real group)'
  END IF
END IF

IF (error_code_ptr <= 0) THEN

  IF (PRESENT(names_fullht_real)) THEN
    IF (PRESENT(dreq_status_fullht_real)) THEN
      CALL set_diag_requests(error_code_ptr, dgroup_fullht_real,               &
                             names_fullht_real, dreq_status_fullht_real,       &
                             diag_requests, n_req_fullht_real,                 &
                             error_message=error_message,                      &
                             error_routine=error_routine)
    ELSE
      error_code_ptr = errcode_diag_mismatch
      message_txt =                                                            &
        'Status flags missing for diagnostic requests (full height real group)'
    END IF
  ELSE
    IF (PRESENT(dreq_status_fullht_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt =                                                            &
        'Status flag array present but no diagnostic name array ' //           &
        '(full height real group)'
    END IF
  END IF

END IF

! Initialise ASAD chemical diagnostics if supported for the current UKCA
! configuration (or skip if we know there are no requests to process)

IF (error_code_ptr <=0 .AND. ukca_config%l_asad_chem_diags_support) THEN
  ! Trap the case where ASAD diagnostic requests have already been set up in
  ! a previous call since repeating the ASAD diagnostic setup during the run
  ! is not currently supported. It is intended that this restriction will be
  ! removed in future.
  IF (n_asad_diags /= 0) THEN
    error_code_ptr = errcode_diag_req_unsupported_use
    message_txt =                                                              &
      'ASAD diagnostic requests are already set up; ' //                       &
      'subsequent changes are not currently supported'
  END IF
  ! Check requests (including UM legacy requests if parent is the UM) to see
  ! which ASAD chemical diagnostics are required (if any) and set up STASH
  ! codes and other info for processing them.
  IF (error_code_ptr <= 0 .AND. (ukca_config%l_enable_diag_um .OR.             &
      n_req_flat_real > 0 .OR. n_req_fullht_real > 0)) THEN
    CALL asad_setstash_chemdiag(ukca_config%row_length, ukca_config%rows,      &
                                ukca_config%model_levels, diag_requests,       &
                                master_diag_list)
    ! Initialise info for chemical fluxes contributing to any required
    ! diagnostics
    IF (n_asad_diags > 0) CALL asad_init_chemdiag()
  END IF
END IF

! Call error report here if an error has been trapped locally (i.e. not in a
! subroutine call) as indicated by presence of message text.
IF (error_code_ptr > 0 .AND. message_txt /= '') THEN
  CALL error_report(ukca_config%i_error_method, error_code_ptr, message_txt,   &
                    RoutineName, msg_out=error_message, locn_out=error_routine)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE ukca_set_diagnostic_requests

! ----------------------------------------------------------------------
SUBROUTINE set_diag_requests(error_code_ptr, group, varnames,                  &
                             status_flags, diag_requests, n_req,               &
                             error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Check availability and create new request data for the given group.
!   Status is updated for any requested diagnostics that are unavailable
!   in the current UKCA configuration.
! ----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_fieldname_mod, ONLY: maxlen_diagname

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr  ! Pointer to return code
INTEGER, INTENT(IN) :: group                    ! Diagnostic group

CHARACTER(LEN=*), INTENT(IN) :: varnames(:)     ! List of diagnostic names

INTEGER, INTENT(IN OUT) :: status_flags(:)      ! Status flag for each
                                                ! diagnostic in list

TYPE(diag_requests_type), INTENT(IN OUT) :: diag_requests(:)
                                                ! Diagnostic request info
                                                ! by group

INTEGER, INTENT(OUT) :: n_req                   ! Number of diagnostic requests
                                                ! in group

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: i
INTEGER :: j

INTEGER, ALLOCATABLE :: i_master(:)  ! Index in master diagnostics list

LOGICAL :: l_found  ! True if diagnostic name found in master list

CHARACTER(LEN=maxlen_diagname), ALLOCATABLE :: varnames_checked(:)
                    ! List of accepted names against which to check for
                    ! duplication.
CHARACTER(LEN=maxlen_message) :: message_txt  ! Buffer for output message

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_DIAG_REQUESTS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

message_txt = ''

n_req = SIZE(varnames)

! Check number of status flags matches number of requests
IF (SIZE(status_flags) /= n_req) THEN
  error_code_ptr = errcode_diag_mismatch
  WRITE(message_txt, '(A,I0,A)')                                               &
    'Wrong number of status flags for group ', group,                          &
    ' diagnostic requests'
  CALL error_report(ukca_config%i_error_method, error_code_ptr, message_txt,   &
                    RoutineName, msg_out=error_message, locn_out=error_routine)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Validate each request in turn and find index in master list

ALLOCATE(i_master(n_req))
IF (n_req > 0) THEN
  i_master(:) = 0
  ALLOCATE(varnames_checked(n_req))
  varnames_checked = ''
END IF

DO i = 1, n_req

  ! Diagnostic name must not be a duplicate
  IF (ANY(varnames_checked(:) == varnames(i))) THEN
    error_code_ptr = errcode_diag_req_duplicate
    WRITE(message_txt, '(A,I0,A,A,A)')                                         &
      'Duplicate group ', group, ' request found for diagnostic ''',           &
      TRIM(varnames(i)), ''''
    CALL error_report(ukca_config%i_error_method, error_code_ptr, message_txt, &
                    RoutineName, msg_out=error_message, locn_out=error_routine)
    IF (lhook) THEN
      CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
    END IF
    RETURN
  ELSE
    varnames_checked(i) = varnames(i)
  END IF

  ! Check that status flag has an expected value

  IF (status_flags(i) /= diag_status_inactive .AND.                            &
      status_flags(i) /= diag_status_requested) THEN
    error_code_ptr = errcode_value_invalid
    WRITE(message_txt, '(A,I0,A,I0,A)')                                        &
      'Unexpected value of status flag (', status_flags(i), ') in group ',     &
      group, ' diagnostic request'
    CALL error_report(ukca_config%i_error_method, error_code_ptr, message_txt, &
                    RoutineName, msg_out=error_message, locn_out=error_routine)
    IF (lhook) THEN
      CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
    END IF
    RETURN
  END IF

  ! Check for diagnostic name in master list

  l_found = .FALSE.
  j = 0
  DO WHILE (j < SIZE(master_diag_list) .AND. .NOT. l_found)
    j = j + 1
    IF (master_diag_list(j)%varname == varnames(i)) THEN
      ! Diagnostic name is recognised so check compatibility
      IF (master_diag_list(j)%group /= group .AND.                             &
          master_diag_list(j)%group_alt /= group) THEN
        error_code_ptr = errcode_diag_mismatch
        WRITE(message_txt, '(A,A,A,I0)')                                       &
          'Diagnostic ''', TRIM(varnames(i)),                                  &
          ''' is not compatible with request group ', group
        CALL error_report(ukca_config%i_error_method, error_code_ptr,          &
                          message_txt, RoutineName, msg_out=error_message,     &
                          locn_out=error_routine)
        IF (lhook) THEN
          CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
        END IF
        RETURN
      END IF
      ! Compatible - Add master list index to request for later referencing
      l_found = .TRUE.
      i_master(i) = j
      ! Modify status flag if diagnostic is unavailable
      IF (.NOT. master_diag_list(j)%l_available)                               &
        status_flags(i) = diag_status_unavailable
    END IF
  END DO

  IF (.NOT. l_found) THEN
    error_code_ptr = errcode_diag_req_unknown
    CALL error_report(ukca_config%i_error_method, error_code_ptr,              &
           'Diagnostic name ''' // TRIM(varnames(i)) //                        &
           ''' is not recognised',                                             &
           RoutineName, msg_out=error_message, locn_out=error_routine)
    IF (lhook) THEN
      CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
    END IF
    RETURN
  END IF

END DO

! Copy diagnostic requests to UKCA data structure
ALLOCATE(diag_requests(group)%varnames(n_req))
ALLOCATE(diag_requests(group)%status_flags(n_req))
ALLOCATE(diag_requests(group)%i_master(n_req))
diag_requests(group)%varnames = varnames
diag_requests(group)%status_flags = status_flags
diag_requests(group)%i_master = i_master

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE set_diag_requests

! ----------------------------------------------------------------------
SUBROUTINE ukca_update_diagnostic_requests(error_code,                         &
                                           dreq_status_flat_real,              &
                                           dreq_status_fullht_real,            &
                                           error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!
!   UKCA API procedure to activate or deactivate individual diagnostic
!   requests via their status flags. UKCA's internal request status flags
!   for the diagnostic groups are updated according to the array(s) of
!   status flags specified.
!
!   Status flag values 'diag_status_requested' and 'diag_status_inactive'
!   are used to activate or deactivate requests. They are ignored if given
!   for any unavailable diagnostics.
!
!   The value 'diag_status_unavailable' is also allowed on input but only
!   where it corresponds to diagnostic requests with that status internally;
!   it cannot be used to update the availability because this is
!   determined by UKCA itself.
!
!   In addition, input status flag values can be negative to avoid updating
!   the status of specific requests (the actual value being arbitrary).
!
!   The new status flag values for each request are returned in the status
!   flag arrays for reference.
!
! Developer's note:
!   The 'error_code' argument can be made optional since it is not required
!   if UKCA is configured to abort on error. It remains mandatory at
!   present for consistency with other API procedures that do not yet
!   support this functionality.
! ----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: ukca_config, l_ukca_config_available

IMPLICIT NONE

! Subroutine arguments

INTEGER, TARGET, INTENT(OUT) :: error_code

INTEGER, OPTIONAL, INTENT(IN OUT) :: dreq_status_flat_real(:)
                                       ! Status flags for requested
                                       ! diagnostics in flat real group
INTEGER, OPTIONAL, INTENT(IN OUT) :: dreq_status_fullht_real(:)
                                       ! Status flags for requested
                                       ! diagnostics in full height real group

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER, POINTER :: error_code_ptr

LOGICAL :: l_status_flags_present

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_UPDATE_DIAGNOSTIC_REQUESTS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Use parent supplied argument for return code
error_code_ptr => error_code

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check for availability of UKCA configuration data
IF (.NOT. l_ukca_config_available) THEN
  error_code_ptr = errcode_ukca_uninit
  CALL error_report(i_error_method_abort, error_code_ptr,                      &
         'No UKCA configuration has been set up', RoutineName,                 &
         msg_out=error_message, locn_out=error_routine)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Check compatibility of each status flag array with the existing request
! data for each diagnostic group and update flags as indicated ignoring any
! requests for unavailable diagnostics.
! At least one status flag array must be present for the call to make sense

l_status_flags_present = .FALSE.

IF (PRESENT(dreq_status_flat_real)) THEN
  l_status_flags_present = .TRUE.
  CALL update_diag_requests(error_code_ptr, dgroup_flat_real,                  &
                            dreq_status_flat_real,                             &
                            diag_requests,                                     &
                            error_message=error_message,                       &
                            error_routine=error_routine)
END IF

IF (error_code_ptr <= 0 .AND. PRESENT(dreq_status_fullht_real)) THEN
  l_status_flags_present = .TRUE.
  CALL update_diag_requests(error_code_ptr, dgroup_fullht_real,                &
                            dreq_status_fullht_real,                           &
                            diag_requests,                                     &
                            error_message=error_message,                       &
                            error_routine=error_routine)
END IF

IF (.NOT. l_status_flags_present) THEN
  error_code_ptr = errcode_diag_mismatch
  CALL error_report(ukca_config%i_error_method, error_code_ptr,                &
                    'No status flag arrays provided', RoutineName,             &
                     msg_out=error_message, locn_out=error_routine)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE ukca_update_diagnostic_requests

! ----------------------------------------------------------------------
SUBROUTINE update_diag_requests(error_code_ptr, group,                         &
                                status_flags, diag_requests,                   &
                                error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Check compatibility of status flag array with existing request data
!   for the given group and update flags as indicated. Negative status
!   flags in the input array imply no update. Any updates for
!   unavailable diagnostics are ignored.
! ----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: ukca_config

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr  ! Pointer to return code
INTEGER, INTENT(IN) :: group                    ! Diagnostic group

INTEGER, INTENT(IN OUT) :: status_flags(:)      ! Status flags for requested
                                                ! diagnostics in group

TYPE(diag_requests_type), INTENT(IN OUT) :: diag_requests(:)
                                                ! Diagnostic request info
                                                ! by group

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: n_req  ! Number of diagnostics requested in group
INTEGER :: i

CHARACTER(LEN=maxlen_message) :: message_txt  ! Buffer for output message

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UPDATE_DIAG_REQUESTS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

message_txt = ''

! Check whether requests have been set up
IF (ALLOCATED(diag_requests(group)%varnames)) THEN
  n_req = SIZE(diag_requests(group)%varnames)
ELSE
  n_req = 0
  error_code_ptr = errcode_diag_mismatch
  WRITE(message_txt,'(A,I0,A)')                                                &
    'Status flag array present but diagnostic requests are not set ' //        &
    '(group ', group, ')'
END IF

! Check number of status flags matches number of requests
IF (error_code_ptr <=0 .AND. SIZE(status_flags) /= n_req) THEN
  error_code_ptr = errcode_diag_mismatch
  WRITE(message_txt, '(A,I0)')                                                 &
    'Wrong number of status flags in diagnostic request update for group ',    &
    group
END IF

IF (error_code_ptr > 0) THEN
  CALL error_report(ukca_config%i_error_method, error_code_ptr, message_txt,   &
                    RoutineName, msg_out=error_message, locn_out=error_routine)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Validate each input status flag value in turn and update the corresponding
! internal status flag if required

DO i = 1, n_req

  ! Ignore negative status flags
  IF (status_flags(i) >= 0) THEN

    ! Check that status flag has an expected value
    IF (status_flags(i) /= diag_status_inactive .AND.                          &
        status_flags(i) /= diag_status_requested .AND.                         &
        status_flags(i) /= diag_status_unavailable) THEN
      error_code_ptr = errcode_value_invalid
      WRITE(message_txt, '(A,I0,A,I0,A)')                                      &
        'Unexpected value of status flag (', status_flags(i), ') in group ',   &
         group, ' diagnostic request'
      CALL error_report(ukca_config%i_error_method, error_code_ptr,            &
                        message_txt, RoutineName,                              &
                        msg_out=error_message, locn_out=error_routine)
      IF (lhook) THEN
        CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      END IF
      RETURN
    END IF

    ! Update status flag if appropriate.
    ! The parent is not allowed to update the availability status because
    ! this is determined by UKCA itself. An attempt to mark an available
    ! diagnostic as unavailable is trapped as an error but an attempt to
    ! request a diagnostic already marked as unavailable is silently ignored,
    ! consistent with the behaviour of 'ukca_set_diagnostic_requests'.
    IF (status_flags(i) == diag_status_unavailable) THEN
      IF (diag_requests(group)%status_flags(i) /=                              &
          diag_status_unavailable) THEN
        error_code_ptr = errcode_value_invalid
        WRITE(message_txt, '(A,A,A,I0)')                                       &
          'Status flag set to ''unavailable'' for available diagnostic ''',    &
          TRIM(diag_requests(group)%varnames(i)), ''' in group ', group
        CALL error_report(ukca_config%i_error_method, error_code_ptr,          &
                          message_txt, RoutineName,                            &
                          msg_out=error_message, locn_out=error_routine)
        IF (lhook) THEN
          CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
        END IF
        RETURN
      END IF
    ELSE IF (diag_requests(group)%status_flags(i) /=                           &
             diag_status_unavailable) THEN
      ! Update internal status flag
      diag_requests(group)%status_flags(i) = status_flags(i)
    END IF

  END IF  ! status_flags(i) >= 0

END DO

! Return updated status flag array
status_flags = diag_requests(group)%status_flags

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE update_diag_requests

! ----------------------------------------------------------------------
SUBROUTINE ukca_get_diagnostic_request_info(error_code,                        &
                                            names_flat_real_ptr,               &
                                            names_fullht_real_ptr,             &
                                            dreq_status_flat_real_ptr,         &
                                            dreq_status_fullht_real_ptr,       &
                                            error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   UKCA API procedure that returns information about the current set of
!   diagnostic requests. Names and/or request status flags can be
!   retrieved for one or both diagnostic groups.
!   This routine may be called at any time after the call to 'ukca_setup'.
!
! Developer's note:
!   The 'error_code' can be made optional since it is not required if
!   UKCA is configured to abort on error. It remains mandatory at
!   present for consistency with other API procedures that do not yet
!   support this functionality.
! ----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: l_ukca_config_available

IMPLICIT NONE

! Subroutine arguments

INTEGER, TARGET, INTENT(OUT) :: error_code

CHARACTER(LEN=*), POINTER, OPTIONAL, INTENT(OUT) :: names_flat_real_ptr(:)
                                    ! Names of requested diagnostics in
                                    ! flat real group
CHARACTER(LEN=*), POINTER, OPTIONAL, INTENT(OUT) :: names_fullht_real_ptr(:)
                                    ! Names of requested diagnostics in
                                    ! full height real group

INTEGER, POINTER, OPTIONAL, INTENT(OUT) :: dreq_status_flat_real_ptr(:)
                                    ! Status flags for requested
                                    ! diagnostics in flat real group
INTEGER, POINTER, OPTIONAL, INTENT(OUT) :: dreq_status_fullht_real_ptr(:)
                                    ! Status flags for requested
                                    ! diagnostics in full height real group

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER, POINTER :: error_code_ptr

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_GET_DIAGNOSTIC_REQUEST_INFO'

! End of header

! Use parent error code argument
error_code_ptr => error_code

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check for availability of UKCA configuration data
IF (.NOT. l_ukca_config_available) THEN
  error_code_ptr = errcode_ukca_uninit
  CALL error_report(i_error_method_abort, error_code_ptr,                      &
         'No UKCA configuration has been set up', RoutineName,                 &
         msg_out=error_message, locn_out=error_routine)
  IF (PRESENT(names_flat_real_ptr)) NULLIFY(names_flat_real_ptr)
  IF (PRESENT(names_fullht_real_ptr)) NULLIFY(names_fullht_real_ptr)
  IF (PRESENT(dreq_status_flat_real_ptr)) NULLIFY(dreq_status_flat_real_ptr)
  IF (PRESENT(dreq_status_fullht_real_ptr)) NULLIFY(dreq_status_fullht_real_ptr)
  RETURN
END IF

! Assign pointers to the reference lists. Nullify if no list is set up.
IF (PRESENT(names_flat_real_ptr)) THEN
  IF (ALLOCATED(diag_requests(dgroup_flat_real)%varnames)) THEN
    names_flat_real_ptr => diag_requests(dgroup_flat_real)%varnames
  ELSE
    NULLIFY(names_flat_real_ptr)
  END IF
END IF
IF (PRESENT(names_fullht_real_ptr)) THEN
  IF (ALLOCATED(diag_requests(dgroup_fullht_real)%varnames)) THEN
    names_fullht_real_ptr => diag_requests(dgroup_fullht_real)%varnames
  ELSE
    NULLIFY(names_fullht_real_ptr)
  END IF
END IF
IF (PRESENT(dreq_status_flat_real_ptr)) THEN
  IF (ALLOCATED(diag_requests(dgroup_flat_real)%status_flags)) THEN
    dreq_status_flat_real_ptr => diag_requests(dgroup_flat_real)%status_flags
  ELSE
    NULLIFY(dreq_status_flat_real_ptr)
  END IF
END IF
IF (PRESENT(dreq_status_fullht_real_ptr)) THEN
  IF (ALLOCATED(diag_requests(dgroup_fullht_real)%status_flags)) THEN
    dreq_status_fullht_real_ptr =>                                             &
      diag_requests(dgroup_fullht_real)%status_flags
  ELSE
    NULLIFY(dreq_status_fullht_real_ptr)
  END IF
END IF

RETURN
END SUBROUTINE ukca_get_diagnostic_request_info

END MODULE ukca_diagnostics_requests_mod
