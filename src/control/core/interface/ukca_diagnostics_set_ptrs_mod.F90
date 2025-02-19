! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module containing subroutines for assigning pointers used to access
!   diagnostic output arrays provided by the parent application.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds,
! University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_diagnostics_set_ptrs_mod

USE ukca_diagnostics_type_mod, ONLY: dgroup_flat_real, dgroup_fullht_real,     &
                                     diagnostics_type

USE ukca_diagnostics_requests_mod, ONLY: diag_requests

USE ukca_config_specification_mod, ONLY: ukca_config_spec_type,                &
                                         diag2d_copy_out, diag3d_copy_out
USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_diag_mismatch, error_report

USE ukca_missing_data_mod, ONLY: imdi

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_DIAGNOSTICS_SET_PTRS_MOD'

! Public procedures
PUBLIC set_diag_ptrs_1d_domain, set_diag_ptrs_3d_domain

! Common error message text used in multiple places
CHARACTER(LEN=*), PARAMETER :: message_txt_data_array_redundant =              &
  'Supplied diagnostic data array is redundant when using callback ' //        &
  'routine for output'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE set_diag_ptrs_1d_domain(error_code_ptr, ukca_config, diagnostics,   &
                                   data_flat_real, data_fullht_real,           &
                                   error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Set up pointers in the diagnostics structure to access the request
!   information and the output arrays (if any) for each diagnostic group.
!   This variant handles the case where the output arrays are for
!   a column domain.
!
! N.B. The code allows for the case of zero requests for a particular group
! being indicated by the request array being absent, i.e. unallocated
! or having zero size. If the latter then consistent diagnostic data
! arrays are expected to be present (with the exception noted below).
! If the request array is absent, the corresponding data array must also
! be absent.
! The diagnostic data array for a group is redundant when UKCA is
! configured to copy its diagnostics to an alternative workspace using
! a callback routine. In that case the data array must be absent.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr  ! Pointer to return code

TYPE(ukca_config_spec_type), INTENT(IN) :: ukca_config
                                                ! UKCA configuration data

TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics
                                                ! Diagnostic request info and
                                                ! pointers to output arrays

REAL, TARGET, OPTIONAL, INTENT(IN) :: data_flat_real(:)
                                                ! Output array for flat real
                                                ! group diagnostics
REAL, TARGET, OPTIONAL, INTENT(IN) :: data_fullht_real(:,:)
                                                ! Output array for full height
                                                ! real group diagnostics

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: group

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=maxlen_message) :: message_txt  ! Buffer for output message

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_DIAG_PTRS_1D_DOMAIN'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Use the global diagnostic request information
diagnostics%requests_ptr => diag_requests

! Check that presence of required diagnostic data arrays and their outer
! dimensions are consistent with current diagnostic requests.
! For full height diagnostics, also check extent of array's spatial dimension.
! (Record current number of requests for reference.)

! Flat real group:

IF (ALLOCATED(diagnostics%requests_ptr(dgroup_flat_real)%varnames)) THEN

  ! A request array (size 0 or greater) has been set up for this group

  diagnostics%n_request(dgroup_flat_real) =                                    &
    SIZE(diagnostics%requests_ptr(dgroup_flat_real)%varnames)

  ! Check data array

  IF (ASSOCIATED(diag2d_copy_out)) THEN
    ! Parent provides a callback routine to copy diagnostic output to
    ! alternative workspace so there shouldn't be a specified data array
    ! for output
    IF (PRESENT(data_flat_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt = message_txt_data_array_redundant // ' (flat real group)'
    END IF
  ELSE
    ! Check data array is present and consistent with request
    IF (.NOT. PRESENT(data_flat_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt = 'Missing diagnostic output array (flat real group)'
    ELSE IF (SIZE(data_flat_real) /=                                           &
             diagnostics%n_request(dgroup_flat_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      WRITE(message_txt,'(A,I0,A,I0,A)')                                       &
        'Output array size (', SIZE(data_flat_real),                           &
        ') does not match no. of requested diagnostics in flat real group (',  &
        diagnostics%n_request(dgroup_flat_real), ')'
    END IF
  END IF

ELSE

  ! No request array set up for this group

  diagnostics%n_request(dgroup_flat_real) = 0

  ! Check data array is absent
  IF (PRESENT(data_flat_real)) THEN
    error_code_ptr = errcode_diag_mismatch
    IF (ASSOCIATED(diag2d_copy_out)) THEN
      message_txt = message_txt_data_array_redundant // ' (flat real group)'
    ELSE
      message_txt =                                                            &
        'Array present for unrequested diagnostic output (flat real group)'
    END IF
  END IF

END IF

IF (error_code_ptr > 0) THEN
  CALL error_report(ukca_config%i_error_method, error_code_ptr, message_txt,   &
                    RoutineName, msg_out=error_message, locn_out=error_routine)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Full height real group:

IF (ALLOCATED(diagnostics%requests_ptr(dgroup_fullht_real)%varnames)) THEN

  ! A request array (size 0 or greater) has been set up for this group

  diagnostics%n_request(dgroup_fullht_real) =                                  &
    SIZE(diagnostics%requests_ptr(dgroup_fullht_real)%varnames)

  ! Check data array

  IF (ASSOCIATED(diag3d_copy_out)) THEN
    ! Parent provides a callback routine to copy diagnostic output to
    ! alternative workspace so there shouldn't be a specified data array
    ! for output
    IF (PRESENT(data_fullht_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt = message_txt_data_array_redundant //                        &
                    ' (full height real group)'
    END IF
  ELSE
    ! Check data array is present and consistent with request
    IF (.NOT. PRESENT(data_fullht_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt = 'Missing diagnostic output array (full height real group)'
    ELSE IF (SIZE(data_fullht_real, DIM=2) /=                                  &
             diagnostics%n_request(dgroup_fullht_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      WRITE(message_txt,'(A,I0,A,I0,A)')                                       &
        'Output array DIM 2 size (', SIZE(data_fullht_real, DIM=2),            &
        ') does not match no. of requested diagnostics in full height ' //     &
        'real group (', diagnostics%n_request(dgroup_fullht_real), ')'
    ELSE IF (SIZE(data_fullht_real, DIM=1) /= ukca_config%model_levels) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt =                                                            &
        'Output array is inconsistent with no. of vertical levels ' //         &
        '(full height real group)'
    END IF
  END IF

ELSE

  ! No request array set up for this group

  diagnostics%n_request(dgroup_fullht_real) = 0

  ! Check data array is absent
  IF (PRESENT(data_fullht_real)) THEN
    error_code_ptr = errcode_diag_mismatch
    IF (ASSOCIATED(diag3d_copy_out)) THEN
      message_txt = message_txt_data_array_redundant //                        &
                    ' (full height real group)'
    ELSE
      message_txt =                                                            &
        'Array present for unrequested diagnostic output ' //                  &
        '(full height real group)'
    END IF
  END IF

END IF

IF (error_code_ptr > 0) THEN
  CALL error_report(ukca_config%i_error_method, error_code_ptr, message_txt,   &
                    RoutineName, msg_out=error_message, locn_out=error_routine)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Record number of spatial dimensions in use for each diagnostic group
! and set pointers

IF (ASSOCIATED(diag2d_copy_out)) THEN
  diagnostics%dimension_out(dgroup_flat_real) = imdi
ELSE
  diagnostics%dimension_out(dgroup_flat_real) = 0
END IF

IF (ASSOCIATED(diag3d_copy_out)) THEN
  diagnostics%dimension_out(dgroup_fullht_real) = imdi
ELSE
  diagnostics%dimension_out(dgroup_fullht_real) = 1
END IF

IF (PRESENT(data_flat_real)) THEN
  diagnostics%value_0d_real_ptr => data_flat_real
ELSE
  NULLIFY(diagnostics%value_0d_real_ptr)
END IF

IF (PRESENT(data_fullht_real)) THEN
  diagnostics%value_1d_real_ptr => data_fullht_real
ELSE
  NULLIFY(diagnostics%value_1d_real_ptr)
END IF

NULLIFY(diagnostics%value_2d_real_ptr)
NULLIFY(diagnostics%value_3d_real_ptr)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE set_diag_ptrs_1d_domain

! ----------------------------------------------------------------------
SUBROUTINE set_diag_ptrs_3d_domain(error_code_ptr, ukca_config, diagnostics,   &
                                   data_flat_real, data_fullht_real,           &
                                   error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Set up pointers in the diagnostics structure to access the request
!   information and the output arrays (if any) for each diagnostic group.
!   This variant handles the case where the output arrays are for
!   a 3D domain.
!
! N.B. The code allows for the case of zero requests for a particular group
! being indicated by the request array being absent, i.e. unallocated
! or having zero size. If the latter then consistent diagnostic data
! arrays are expected to be present (with the exception noted below).
! If the request array is absent, the corresponding data array must also
! be absent.
! The diagnostic data array for a group is redundant when UKCA is
! configured to copy its diagnostics to an alternative workspace using
! a callback routine. In that case the data array must be absent.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr  ! Pointer to return code

TYPE(ukca_config_spec_type), INTENT(IN) :: ukca_config
                                                ! UKCA configuration data

TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics
                                                ! Diagnostic request info and
                                                ! pointers to output arrays

REAL, TARGET, OPTIONAL, INTENT(IN) :: data_flat_real(:,:,:)
                                                ! Output array for flat real
                                                ! group diagnostics
REAL, TARGET, OPTIONAL, INTENT(IN) :: data_fullht_real(:,:,:,:)
                                                ! Output array for full height
                                                ! real group diagnostics

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: group

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=maxlen_message) :: message_txt  ! Buffer for output message

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_DIAG_PTRS_3D_DOMAIN'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Use the global diagnostic request information
diagnostics%requests_ptr => diag_requests

! Check that presence of required diagnostic data arrays and their outer
! dimensions are consistent with current diagnostic requests and check
! extent of the arrays' spatial dimensions.
! (Record current number of requests for reference.)

! Flat real group:

IF (ALLOCATED(diagnostics%requests_ptr(dgroup_flat_real)%varnames)) THEN

  ! A request array (size 0 or greater) has been set up for this group

  diagnostics%n_request(dgroup_flat_real) =                                    &
    SIZE(diagnostics%requests_ptr(dgroup_flat_real)%varnames)

  ! Check data array

  IF (ASSOCIATED(diag2d_copy_out)) THEN
    ! Parent provides a callback routine to copy diagnostic output to
    ! alternative workspace so there shouldn't be a specified data array
    ! for output
    IF (PRESENT(data_flat_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt = message_txt_data_array_redundant // ' (flat real group)'
    END IF
  ELSE
    ! Check data array is present and consistent with request
    IF (.NOT. PRESENT(data_flat_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt = 'Missing diagnostic output array (flat real group)'
    ELSE IF (SIZE(data_flat_real, DIM=3) /=                                    &
             diagnostics%n_request(dgroup_flat_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      WRITE(message_txt,'(A,I0,A,I0,A)')                                       &
        'Output array DIM 3 size (', SIZE(data_flat_real, DIM=3),              &
        ') does not match no. of requested diagnostics in flat real group (',  &
        diagnostics%n_request(dgroup_flat_real), ')'
    ELSE IF (SIZE(data_flat_real, DIM=1) /= ukca_config%row_length .OR.        &
             SIZE(data_flat_real, DIM=2) /= ukca_config%rows) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt =                                                            &
        'Output array is inconsistent with domain extent (flat real group)'
    END IF
  END IF

ELSE

  ! No request array set up for this group

  diagnostics%n_request(dgroup_flat_real) = 0

  ! Check data array is absent
  IF (PRESENT(data_flat_real)) THEN
    error_code_ptr = errcode_diag_mismatch
    IF (ASSOCIATED(diag2d_copy_out)) THEN
      message_txt = message_txt_data_array_redundant // ' (flat real group)'
    ELSE
      message_txt =                                                            &
        'Array present for unrequested diagnostic output (flat real group)'
    END IF
  END IF

END IF

IF (error_code_ptr > 0) THEN
  CALL error_report(ukca_config%i_error_method, error_code_ptr, message_txt,   &
                    RoutineName, msg_out=error_message, locn_out=error_routine)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Full height real group:

IF (ALLOCATED(diagnostics%requests_ptr(dgroup_fullht_real)%varnames)) THEN

  ! A request array (size 0 or greater) has been set up for this group

  diagnostics%n_request(dgroup_fullht_real) =                                  &
    SIZE(diagnostics%requests_ptr(dgroup_fullht_real)%varnames)

  ! Check data array

  IF (ASSOCIATED(diag3d_copy_out)) THEN
    ! Parent provides a callback routine to copy diagnostic output to
    ! alternative workspace so there shouldn't be a specified data array
    ! for output
    IF (PRESENT(data_fullht_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt = message_txt_data_array_redundant //                        &
                    ' (full height real group)'
    END IF
  ELSE
    ! Check data array is present and consistent with request
    IF (.NOT. PRESENT(data_fullht_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt = 'Missing diagnostic output (full height real group)'
    ELSE IF (SIZE(data_fullht_real, DIM=4) /=                                  &
             diagnostics%n_request(dgroup_fullht_real)) THEN
      error_code_ptr = errcode_diag_mismatch
      WRITE(message_txt,'(A,I0,A,I0,A)')                                       &
        'Output array DIM 4 size (', SIZE(data_fullht_real, DIM=4),            &
        ') does not match no. of requested diagnostics in full height ' //     &
        'real group (', diagnostics%n_request(dgroup_fullht_real), ')'
    ELSE IF (SIZE(data_fullht_real, DIM=1) /= ukca_config%row_length .OR.      &
             SIZE(data_fullht_real, DIM=2) /= ukca_config%rows .OR.            &
             SIZE(data_fullht_real, DIM=3) /= ukca_config%model_levels) THEN
      error_code_ptr = errcode_diag_mismatch
      message_txt =                                                            &
        'Output array inconsistent with domain extent ' //                     &
        '(full height real group)'
    END IF
  END IF

ELSE

  ! No request array set up for this group

  diagnostics%n_request(dgroup_fullht_real) = 0

  ! Check data array is absent
  IF (PRESENT(data_fullht_real)) THEN
    error_code_ptr = errcode_diag_mismatch
    IF (ASSOCIATED(diag3d_copy_out)) THEN
      message_txt = message_txt_data_array_redundant //                        &
                    ' (full height real group)'
    ELSE
      message_txt =                                                            &
        'Array present for unrequested diagnostic output ' //                  &
        '(full height real group)'
    END IF
  END IF

END IF

IF (error_code_ptr > 0) THEN
  CALL error_report(ukca_config%i_error_method, error_code_ptr, message_txt,   &
                    RoutineName, msg_out=error_message, locn_out=error_routine)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Record number of spatial dimensions in use for each diagnostic group
! and set pointers

IF (ASSOCIATED(diag2d_copy_out)) THEN
  diagnostics%dimension_out(dgroup_flat_real) = imdi
ELSE
  diagnostics%dimension_out(dgroup_flat_real) = 2
END IF

IF (ASSOCIATED(diag3d_copy_out)) THEN
  diagnostics%dimension_out(dgroup_fullht_real) = imdi
ELSE
  diagnostics%dimension_out(dgroup_fullht_real) = 3
END IF

IF (PRESENT(data_flat_real)) THEN
  diagnostics%value_2d_real_ptr => data_flat_real
ELSE
  NULLIFY(diagnostics%value_2d_real_ptr)
END IF

IF (PRESENT(data_fullht_real)) THEN
  diagnostics%value_3d_real_ptr => data_fullht_real
ELSE
  NULLIFY(diagnostics%value_3d_real_ptr)
END IF

NULLIFY(diagnostics%value_0d_real_ptr)
NULLIFY(diagnostics%value_1d_real_ptr)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE set_diag_ptrs_3d_domain

END MODULE ukca_diagnostics_set_ptrs_mod

