! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module providing parameters for UKCA error handling.
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

MODULE ukca_error_mod

IMPLICIT NONE
PUBLIC

INTEGER, PARAMETER :: maxlen_message = 256
INTEGER, PARAMETER :: maxlen_procname = 80
INTEGER, PARAMETER :: errcode_ukca_uninit = 1
INTEGER, PARAMETER :: errcode_tracer_req_uninit = 2
INTEGER, PARAMETER :: errcode_tracer_mismatch = 3
INTEGER, PARAMETER :: errcode_ntp_uninit = 4
INTEGER, PARAMETER :: errcode_ntp_mismatch = 5
INTEGER, PARAMETER :: errcode_env_req_uninit = 6
INTEGER, PARAMETER :: errcode_env_field_unknown = 7
INTEGER, PARAMETER :: errcode_env_field_mismatch = 8
INTEGER, PARAMETER :: errcode_env_field_missing = 9
INTEGER, PARAMETER :: errcode_diag_req_unknown = 10
INTEGER, PARAMETER :: errcode_diag_req_duplicate = 11
INTEGER, PARAMETER :: errcode_diag_req_unsupported_use = 12
INTEGER, PARAMETER :: errcode_diag_mismatch = 13
INTEGER, PARAMETER :: errcode_ukca_internal_fault = 14
INTEGER, PARAMETER :: errcode_value_unknown = 15
INTEGER, PARAMETER :: errcode_value_invalid = 16
INTEGER, PARAMETER :: errcode_value_missing = 17

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE error_report(error_code_ptr, msg_in, locn_in,                       &
                        msg_out, locn_out)
! ----------------------------------------------------------------------
! Description:
!   Handle UKCA error trap.
!
! Method:
!
!   For fatal errors (positive error code):
!   If UKCA is configured to abort on trapping an error then write error
!   information and abort (via call to 'ereport'). Otherwise, return
!   control to the calling routine. The error message can be printed as
!   a warning before returning control if this action is specified in the
!   configuration.
!   A copy of the error information is returned as required according to
!   the presence or absence of the arguments 'msg_out' and/or 'locn_out'.
!
!   For warnings (negative error code):
!   Write error information (via call to 'ereport') and return control
!   to the calling routine.
!
!   If a zero error code is passed, treat as a warning with handling of
!   any additional information defined by `ereport`.
! ----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: ukca_config,                          &
                                         i_error_method_abort,                 &
                                         i_error_method_return,                &
                                         i_error_method_warn_and_return

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr

CHARACTER(LEN=*), INTENT(IN) :: msg_in
CHARACTER(LEN=*), INTENT(IN) :: locn_in

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: msg_out
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: locn_out

! Local variables

INTEGER :: ereport_code  ! Non-pointer error code for 'ereport' call

CHARACTER(LEN=maxlen_message) :: ereport_msg  ! Modified error message

CHARACTER(LEN=*), PARAMETER :: RoutineName='ERROR_REPORT'

! End of header

IF (error_code_ptr <= 0) THEN

  ! Print warning message
  ereport_code = error_code_ptr
  CALL ereport(locn_in, ereport_code, msg_in)
  ! Control is returned from 'ereport' with error code
  ! reset to zero. No error info needed.
  error_code_ptr = ereport_code
  IF (PRESENT(msg_out)) msg_out = ''
  IF (PRESENT(locn_out)) locn_out = ''

ELSE

  SELECT CASE (ukca_config%i_error_method)

  CASE (i_error_method_abort)
    ! ----------------------------------------------------
    ! Error handling method: write error message and abort
    ! ----------------------------------------------------
    ereport_code = error_code_ptr
    CALL ereport(locn_in, ereport_code, msg_in)

  CASE (i_error_method_return)
    ! -----------------------------------------------
    ! Error handling method: return control to parent
    ! -----------------------------------------------
    ! Nothing to do here

  CASE (i_error_method_warn_and_return)
    ! -----------------------------------------------------
    ! Error handling method: Return control to parent after
    ! printing error message as a warning
    ! -----------------------------------------------------
    ereport_code = - error_code_ptr
    ereport_msg = TRIM(msg_in) // ' - UKCA TERMINATING'
    CALL ereport(locn_in, ereport_code, ereport_msg)

  END SELECT

  ! Copy error information for returning to parent via calling chain
  IF (PRESENT(msg_out)) msg_out = msg_in
  IF (PRESENT(locn_out)) locn_out = locn_in

END IF

RETURN
END SUBROUTINE error_report

END MODULE ukca_error_mod
