! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Diagnostic handling module containing routines to copy diagnostic
!   field data and the corresponding diagnostic status flags to parent
!   arrays during the UKCA time step.
!
!   Public procedures provided:
!
!     update_skipped_diag_flags
!      - set status to 'skipped' where applicable for a non-chemistry time step
!     update_diagnostics_2d_real
!      - service any diagnostic requests satisfied by a given 2d field
!     update_diagnostics_3d_real
!      - service any diagnostic requests satisfied by a given 3d field
!     blank_out_missing_diags
!      - overwrite any non-valid diagnostic fields in the output array(s)
!     init_diag_status
!      - initialise internal diagnostic status flags from request status flags
!     output_diag_status
!      - copy the diagnostic status flag values to parent arrays
!     diag_status_dealloc
!      - deallocate the internal diagnostic status flag arrays
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

MODULE ukca_diagnostics_output_mod

USE ukca_diagnostics_type_mod, ONLY: dgroup_flat_real, dgroup_fullht_real,     &
                                     n_diag_group, diagnostics_type,           &
                                     diag_status_inactive,                     &
                                     diag_status_requested,                    &
                                     diag_status_skipped, diag_status_valid

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_diag_mismatch, errcode_ukca_internal_fault,  &
                          errcode_value_unknown, error_report

USE yomhook,        ONLY: lhook, dr_hook
USE parkind1,       ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_DIAGNOSTICS_OUTPUT_MOD'

! Public procedures
PUBLIC update_skipped_diag_flags, seek_active_requests,                        &
       update_diagnostics_2d_real, update_diagnostics_3d_real,                 &
       blank_out_missing_diags, init_diag_status, output_diag_status,          &
       diag_status_dealloc

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE update_skipped_diag_flags(diagnostics)
! ----------------------------------------------------------------------
! Description:
!   Set status of all requests for diagnostics that are only output on
!   chemistry time steps to 'diag_status_skipped'. These diagnostics are
!   identified by the 'l_chem_timesteps' switch in the master diagnostics
!   list (see ukca_init_diagnostics_mod).
! ----------------------------------------------------------------------

USE ukca_diagnostics_master_mod, ONLY: master_diag_list

IMPLICIT NONE

! Subroutine arguments
TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics  ! Diagnostics request
                                                       ! info and pointers
                                                       ! to output arrays

! Local variables

INTEGER :: group
INTEGER :: i
INTEGER :: i_master

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UPDATE_SKIPPED_DIAG_FLAGS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO group = 1, n_diag_group
  DO i = 1, diagnostics%n_request(group)
    IF (diagnostics%outvalue_status(group)%status_flags(i) ==                  &
        diag_status_requested) THEN
      i_master = diagnostics%requests_ptr(group)%i_master(i)
      IF (master_diag_list(i_master)%l_chem_timestep)                          &
        diagnostics%outvalue_status(group)%status_flags(i) =                   &
          diag_status_skipped
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE update_skipped_diag_flags

! ----------------------------------------------------------------------
SUBROUTINE seek_active_requests(diag_name, diagnostics, l_check_group,         &
                                i_diag_req, l_active_requests)
! ----------------------------------------------------------------------
! Description:
!   Check each applicable request group to determine whether a request
!   exists for the named diagnostic field (or a sub-field thereof) and
!   whether any requests found are active. Request indices are returned
!   in 'i_diag_req' for later use. A zero index indicates that the
!   diagnostic is not present in the request list.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

CHARACTER(LEN=*), INTENT(IN) :: diag_name           ! Diagnostic name

TYPE(diagnostics_type), INTENT(IN) :: diagnostics
                                                    ! Diagnostics request info

LOGICAL, INTENT(IN) :: l_check_group(n_diag_group)  ! True to check the request
                                                    ! group (should be true for
                                                    ! all apllicable groups)

INTEGER, INTENT(OUT) :: i_diag_req(n_diag_group)    ! Indices of requests for
                                                    ! the named diagnostic

LOGICAL, INTENT(OUT) :: l_active_requests           ! True if active requests
                                                    ! found

! Local variables

INTEGER :: group

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SEEK_ACTIVE_REQUESTS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

i_diag_req = 0
l_active_requests = .FALSE.


! Check each diagnostic request group
DO group = 1, n_diag_group

  ! Search the list of requests if group is applicable and has a non-empty
  ! request list
  IF (l_check_group(group) .AND. diagnostics%n_request(group) > 0) THEN

    ! Check for a request for the named diagnostic
    i_diag_req(group) =                                                        &
      diag_index(diag_name, diagnostics%requests_ptr(group)%varnames)

    ! If found, check status flag to
    ! update status flag to show output is valid (or ignore if request is
    ! inactive)
    IF (i_diag_req(group) > 0) l_active_requests = .TRUE.

  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE seek_active_requests

! ----------------------------------------------------------------------
SUBROUTINE update_diagnostics_2d_real(error_code_ptr, diag_name,               &
                                      field, diagnostics, i_diag_req,          &
                                      error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Update the diagnostics output by copying data from the given field
!   to the output array if the named diagnostic is requested. Presence
!   of a request is determined from request indices 'i_diag_req' if
!   supplied or from checking the appropriate request list if not.
!   (Only check the flat real group requests since other requests are
!   incompatible with 2D fields.)
!   Update internal status flags accordingly.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr  ! Pointer to return code

CHARACTER(LEN=*), INTENT(IN) :: diag_name       ! Diagnostic name


REAL, INTENT(IN) :: field(:,:)                  ! Field data for output

TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics
                                                ! Diagnostics request info
                                                ! and pointers to output arrays

INTEGER, OPTIONAL, INTENT(IN) :: i_diag_req(n_diag_group)
                                                ! Indices of requests for the
                                                ! named diagnostic

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: i

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UPDATE_DIAGNOSTICS_2D_REAL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Search the list of flat real group requests if present and not empty
IF (diagnostics%n_request(dgroup_flat_real) > 0) THEN

  ! Check for a request for the named diagnostic or use the results of a
  ! previous check if available
  IF (PRESENT(i_diag_req)) THEN
    i = i_diag_req(dgroup_flat_real)
  ELSE
    i = diag_index(diag_name,                                                  &
                   diagnostics%requests_ptr(dgroup_flat_real)%varnames)
  END IF

  ! If found, copy 2d field to output scalar or array provided by parent and
  ! update status flag to show output is valid (or ignore if request is
  ! inactive)
  IF (i > 0) THEN
    CALL update_diag_field(error_code_ptr, dgroup_flat_real, i, diagnostics,   &
                           field_2d=field, error_message=error_message,        &
                           error_routine=error_routine)
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE update_diagnostics_2d_real

! ----------------------------------------------------------------------
SUBROUTINE update_diagnostics_3d_real(error_code_ptr, diag_name,               &
                                      field, diagnostics, i_diag_req,          &
                                      error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Update the diagnostics output by copying data from the given 3D field
!   to the output array (or arrays) if the named diagnostic is requested.
!   Presence of a request is determined from request indices 'i_diag_req'
!   if supplied or from checking the request lists for each group if not.
!   Allow for a request for the full height field and/or a request for
!   the surface level subfield. Update internal status flags accordingly.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr  ! Pointer to return code

CHARACTER(LEN=*), INTENT(IN) :: diag_name       ! Diagnostic name

REAL, INTENT(IN) :: field(:,:,:)                ! Field data for output

TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics
                                                ! Diagnostic request info
                                                ! and pointers to output arrays

INTEGER, OPTIONAL, INTENT(IN) :: i_diag_req(n_diag_group)
                                                ! Indices of requests for the
                                                ! named diagnostic

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: group
INTEGER :: i

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UPDATE_DIAGNOSTICS_3D_REAL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Process flat real and full height real group requests. For a flat real group
! request the surface-level 2d subfield is copied. For a full height real group
! request the 3d field is copied.
DO group = 1, n_diag_group

  ! Search the list of requests if present and not empty
  IF (diagnostics%n_request(group) > 0) THEN

    ! Check for a request for the named diagnostic in the current group or use
    ! the results of a previous check if available
    IF (PRESENT(i_diag_req)) THEN
      i = i_diag_req(group)
    ELSE
      i = diag_index(diag_name, diagnostics%requests_ptr(group)%varnames)
    END IF

    ! If found, copy 3d field or 2d subfield to output scalar or array provided
    ! by parent and update status flag to show output is valid (or ignore if
    ! request is inactive)
    IF (i > 0) THEN
      CALL update_diag_field(error_code_ptr, group, i, diagnostics,            &
                             field_3d=field, error_message=error_message,      &
                             error_routine=error_routine)
    END IF

  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE update_diagnostics_3d_real

! ----------------------------------------------------------------------
SUBROUTINE update_diag_field(error_code_ptr, group, i_diag, diagnostics,       &
                             field_2d, field_3d, error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Copy the 2D or 3D field provided to the appropriate parent array if
!   the request is active and update the request's status flag to show
!   the array entry contains valid data.
! ----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: diag2d_copy_out, diag3d_copy_out
USE ukca_fieldname_mod, ONLY: maxlen_diagname

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr  ! Pointer to return code
INTEGER, INTENT(IN) :: group                    ! Diagnostic group
INTEGER, INTENT(IN) :: i_diag                   ! Index of diagnostic in
                                                ! request list

TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics
                                                ! Diagnostic request info
                                                ! and pointers to output
                                                ! arrays

REAL, OPTIONAL, INTENT(IN) :: field_2d(:,:)     ! Field data for output if 2D
REAL, OPTIONAL, INTENT(IN) :: field_3d(:,:,:)   ! Field data for output if 3D

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

CHARACTER(LEN=maxlen_message) :: message_txt  ! Buffer for output message
CHARACTER(LEN=maxlen_diagname) :: diag_name

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UPDATE_DIAG_FIELD'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

diag_name = diagnostics%requests_ptr(group)%varnames(i_diag)
message_txt = ''

SELECT CASE (diagnostics%outvalue_status(group)%status_flags(i_diag))

CASE (diag_status_requested)

  ! Copy the diagnostic field to the output array
  SELECT CASE (group)
  CASE (dgroup_flat_real)

    IF (PRESENT(field_2d)) THEN
      IF (ASSOCIATED(diag2d_copy_out)) THEN
        ! Use parent callback routine
        CALL diag2d_copy_out(diag_name, field_2d)
      ELSE
        ! Copy to array supplied by parent via the time step API call
        SELECT CASE (diagnostics%dimension_out(group))
        CASE (0)  ! Output diagnostic is scalar
          diagnostics%value_0d_real_ptr(i_diag) = field_2d(1,1)
        CASE (2)  ! Output diagnostic is 2D
          diagnostics%value_2d_real_ptr(:,:,i_diag) = field_2d(:,:)
        CASE DEFAULT
          error_code_ptr = errcode_ukca_internal_fault
          message_txt =                                                        &
            'Invalid no. of dimensions for diagnostic output (flat real group)'
        END SELECT
      END IF
    ELSE IF (PRESENT(field_3d)) THEN
      ! Select surface sub-field
      IF (ASSOCIATED(diag2d_copy_out)) THEN
        ! Use parent callback routine
        CALL diag2d_copy_out(diag_name, field_3d(:,:,1))
      ELSE
        ! Copy to array supplied by parent via the time step API call
        SELECT CASE (diagnostics%dimension_out(group))
        CASE (0)  ! Output diagnostic is scalar
          diagnostics%value_0d_real_ptr(i_diag) = field_3d(1,1,1)
        CASE (2)  ! Output diagnostic is 2D
          diagnostics%value_2d_real_ptr(:,:,i_diag) = field_3d(:,:,1)
        CASE DEFAULT
          error_code_ptr = errcode_ukca_internal_fault
          message_txt =                                                        &
            'Invalid no. of dimensions for diagnostic output (flat real group)'
        END SELECT
      END IF
    ELSE
      error_code_ptr = errcode_ukca_internal_fault
      message_txt = 'No field for diagnostic output (flat real group)'
    END IF

  CASE (dgroup_fullht_real)

    IF (PRESENT(field_3d)) THEN
      IF (ASSOCIATED(diag3d_copy_out)) THEN
        ! Use parent callback routine
        CALL diag3d_copy_out(diag_name, field_3d)
      ELSE
        ! Copy to array supplied by parent via the time step API call
        SELECT CASE (diagnostics%dimension_out(group))
        CASE (1)  ! Output diagnostic is 1D column
          diagnostics%value_1d_real_ptr(:,i_diag) = field_3d(1,1,:)
        CASE (3)  ! Output diagnostic is 3D
          diagnostics%value_3d_real_ptr(:,:,:,i_diag) = field_3d(:,:,:)
        CASE DEFAULT
          error_code_ptr = errcode_ukca_internal_fault
          message_txt =                                                        &
            'Invalid no. of dimensions for diagnostic output ' //              &
            '(full height real group)'
        END SELECT
      END IF
    ELSE
      error_code_ptr = errcode_ukca_internal_fault
      message_txt = 'No field for diagnostic output (full height real group)'
    END IF

  END SELECT

  ! Update the diagnostic status to indicate that the output field is valid
  IF (error_code_ptr <= 0) THEN
    diagnostics%outvalue_status(group)%status_flags(i_diag) = diag_status_valid
  END IF

CASE (diag_status_inactive, diag_status_skipped)
  ! Do nothing since request is either deactivated by the parent or being
  ! skipped on a non-chemistry time step

CASE DEFAULT
  error_code_ptr = errcode_ukca_internal_fault
  WRITE(message_txt,'(A,I0,A,A,A,I0,A)')                                       &
    'Unexpected diagnostic status code (',                                     &
    diagnostics%outvalue_status(group)%status_flags(i_diag),                   &
    ') for diagnostic request (''', TRIM(diag_name), ''' in group ', group, ')'

END SELECT

IF (error_code_ptr > 0) THEN
  CALL error_report(error_code_ptr, message_txt, RoutineName,                  &
                    msg_out=error_message, locn_out=error_routine)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE update_diag_field

! ----------------------------------------------------------------------
INTEGER FUNCTION diag_index(diag_name, varnames)
! Description:
!   Return index of the given name in the given list of names.
!   Return 0 if not present.
! ----------------------------------------------------------------------

USE ukca_fieldname_mod, ONLY: maxlen_diagname

IMPLICIT NONE

! Function arguments
CHARACTER(LEN=*), INTENT(IN) :: diag_name                   ! Name to find
CHARACTER(LEN=maxlen_diagname), INTENT(IN) :: varnames(:)   ! List of names

! Local variables
INTEGER :: i
LOGICAL :: l_found

! End of header

l_found = .FALSE.
i = 0
DO WHILE (i < SIZE(varnames) .AND. .NOT. l_found)
  i = i + 1
  IF (varnames(i) == diag_name) l_found = .TRUE.
END DO

IF (l_found) THEN
  diag_index = i
ELSE
  diag_index = 0
END IF

RETURN
END FUNCTION diag_index

! ----------------------------------------------------------------------
SUBROUTINE blank_out_missing_diags(error_code_ptr, diagnostics,                &
                                   error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Check status of all diagnostic requests and blank out the values in
!   the output array(s) for any diagnostics that do not contain valid
!   data.
! ----------------------------------------------------------------------

USE missing_data_mod, ONLY: rmdi

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr     ! Pointer to return code

TYPE(diagnostics_type), INTENT(IN) :: diagnostics  ! Diagnostic request info
                                                   ! and pointers to output
                                                   ! arrays

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: group
INTEGER :: i
INTEGER :: i_master

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'BLANK_OUT_MISSING_DIAGS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

DO group = 1, n_diag_group
  DO i = 1, diagnostics%n_request(group)
    IF (diagnostics%outvalue_status(group)%status_flags(i) /=                  &
        diag_status_valid) THEN
      SELECT CASE (diagnostics%dimension_out(group))
      CASE (0)
        diagnostics%value_0d_real_ptr(i) = rmdi
      CASE (1)
        diagnostics%value_1d_real_ptr(:,i) = rmdi
      CASE (2)
        diagnostics%value_2d_real_ptr(:,:,i) = rmdi
      CASE (3)
        diagnostics%value_3d_real_ptr(:,:,:,i) = rmdi
      CASE DEFAULT
        error_code_ptr = errcode_ukca_internal_fault
        CALL error_report(error_code_ptr,                                      &
          'Invalid no. of dimensions for diagnostic output',                   &
          RoutineName, msg_out=error_message, locn_out=error_routine)
        IF (lhook) THEN
          CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
        END IF
        RETURN
      END SELECT
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE blank_out_missing_diags

! ----------------------------------------------------------------------
SUBROUTINE init_diag_status(diagnostics)
! ----------------------------------------------------------------------
! Description:
!   Allocate diagnostic status flag arrays and initialise them from the
!   request status flags.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics ! Diagnostic request info
                                                      ! and pointers to output
                                                      ! arrays

! Local variables
INTEGER :: group

DO group = 1, n_diag_group
  IF (diagnostics%n_request(group) > 0) THEN
    ALLOCATE(diagnostics%outvalue_status(group)%status_flags(                  &
             SIZE(diagnostics%requests_ptr(group)%status_flags)))
    diagnostics%outvalue_status(group)%status_flags =                          &
      diagnostics%requests_ptr(group)%status_flags
  END IF
END DO

RETURN
END SUBROUTINE init_diag_status

! ----------------------------------------------------------------------
SUBROUTINE output_diag_status(error_code_ptr, group, diagnostics,              &
                              status_flags, error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Copy diagnostic status flags for the selected diagnostic group
!   to the given output status flag array.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr     ! Pointer to return code
INTEGER :: group                                   ! Diagnostic group

TYPE(diagnostics_type), INTENT(IN) :: diagnostics  ! Diagnostic request info
                                                   ! and pointers to output
                                                   ! arrays

INTEGER, INTENT(OUT) :: status_flags(:)            ! Returned status flags

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: n_req

CHARACTER(LEN=maxlen_message) :: message_txt  ! Buffer for output message

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'OUTPUT_DIAG_STATUS'

! End of header

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

IF (ALLOCATED(diagnostics%requests_ptr(dgroup_flat_real)%varnames)) THEN
  n_req = diagnostics%n_request(group)
  IF (SIZE(status_flags) /= n_req) THEN
    error_code_ptr = errcode_diag_mismatch
    WRITE(message_txt,'(A,I0,A)')                                              &
      'Status flag array size inconsistent with no. of requested ' //          &
      'diagnostics (group ', group, ')'
  END IF
ELSE
  error_code_ptr = errcode_diag_mismatch
  WRITE(message_txt,'(A,I0,A)')                                                &
    'Status flag array present but diagnostic requests are not set ' //        &
    '(group ', group, ')'
END IF

IF (error_code_ptr > 0) THEN
  CALL error_report(error_code_ptr, message_txt, RoutineName,                  &
                    msg_out=error_message, locn_out=error_routine)
ELSE IF (n_req > 0) THEN
  status_flags = diagnostics%outvalue_status(group)%status_flags
END IF

RETURN
END SUBROUTINE output_diag_status

! ----------------------------------------------------------------------
SUBROUTINE diag_status_dealloc(diagnostics)
! ----------------------------------------------------------------------
! Description:
!   Deallocate the diagnostic status flag arrays.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics ! Diagnostic request info
                                                      ! and pointers to output
                                                      ! arrays

! Local variables
INTEGER :: group

DO group = n_diag_group, 1, -1
  IF (ALLOCATED(diagnostics%outvalue_status(group)%status_flags)) THEN
    DEALLOCATE(diagnostics%outvalue_status(group)%status_flags)
  END IF
END DO

RETURN
END SUBROUTINE diag_status_dealloc

END MODULE ukca_diagnostics_output_mod
