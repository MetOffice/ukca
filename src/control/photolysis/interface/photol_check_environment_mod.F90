! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
! Module containing procedures to check that environment fields provided by
! the parent model are complete and compatible for this configuration
!
! This provides the following interface procedures for the Photolysis API.
!
!   check_env_fields_list: Compares the full list of required environment
!     fields against that received from parent and prints out those that are
!     missing.
!
!   check_environ_group - Processes the env fields of a given type
!     (e.g scalar_real, flat_real) at one time and checks for all required
!     fields as well as the spatial dimensions
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds,
! University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_Photolysis
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE photol_check_environment_mod

USE photol_fieldname_mod,   ONLY:  fieldname_len
USE photol_environment_mod, ONLY: n_fields, photol_env_fields,                 &
                                  environ_fldnames_scalar_real,                &
                                  environ_fldnames_flat_integer,               &
                                  environ_fldnames_flat_real,                  &
                                  environ_fldnames_fullht_real,                &
                                  environ_fldnames_fullht0_real,               &
                                  environ_fldnames_fullhtphot_real
USE photol_config_specification_mod, ONLY: photol_config

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname, error_report,       &
                          errcode_value_unknown, errcode_env_req_uninit,       &
                          errcode_env_field_mismatch, errcode_env_field_missing
! Dr Hook modules
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC check_env_fields_list, check_environ_group

! Generic interface for procedures that can handle grouped fields of
! different types and dimensions.
! Note: fullht_real and fullht0_real groups contain fields with 4 dimensions
! (3 spatial + 1 num_fields - only difference is 3rd dimension having one
! extra level) so can be processed via single procedure
! Returns l_fields_this_group=False if group is empty (no fields of this type)
INTERFACE check_environ_group
  MODULE PROCEDURE check_envgroup_scalar_real
  MODULE PROCEDURE check_envgroup_flat_integer
  MODULE PROCEDURE check_envgroup_flat_real
  MODULE PROCEDURE check_envgroup_fullht_real
  MODULE PROCEDURE check_envgroup_fullhtphot_real
END INTERFACE check_environ_group

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'PHOTOL_CHECK_ENVIRONMENT_MOD'

CONTAINS

! Compares input fields list against fields required for this run
!------------------------------------------------------------------------------
SUBROUTINE check_env_fields_list(envfield_names_in, error_code_ptr,            &
                                 error_message, error_routine)

USE umPrintMgr, ONLY: umMessage, UmPrint

IMPLICIT NONE

CHARACTER(LEN=fieldname_len), INTENT(IN)   :: envfield_names_in(:)
INTEGER, POINTER, INTENT(IN) :: error_code_ptr
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

CHARACTER(LEN=maxlen_message) :: err_msg
INTEGER :: n, n_missing

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CHK_ENV_FIELDS_LIST'
! ------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
err_msg = ''
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check number of fields passed
n_missing = 0
DO n = 1, n_fields
  IF ( ANY(envfield_names_in == photol_env_fields(n)%env_fieldname) ) THEN
    ! Do nothing - found field name
  ELSE
    error_code_ptr = errcode_env_field_missing
    WRITE(umMessage,'(A,A)')' MISSING: ',photol_env_fields(n)%env_fieldname
    CALL umPrint(umMessage)
    n_missing = n_missing + 1
  END IF
END DO

IF ( error_code_ptr /= 0 .OR. n_missing > 0 ) THEN
  WRITE(err_msg,'(I0,A)') n_missing,                                           &
  ' required Photolysis environ fields MISSING from provided list.'
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
          RoutineName, msg_out=error_message, locn_out=error_routine)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_env_fields_list

! Checks that provided fields of type 'scalar real' match the requirements
!------------------------------------------------------------------------------
SUBROUTINE check_envgroup_scalar_real(flds_in, group_name, error_code_ptr,     &
             l_fields_this_group, error_message, error_routine)

IMPLICIT NONE

REAL, INTENT(IN)    :: flds_in(:)
CHARACTER(LEN=*), INTENT(IN) :: group_name
LOGICAL, INTENT(OUT) :: l_fields_this_group

INTEGER, POINTER, INTENT(IN) :: error_code_ptr
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

CHARACTER(LEN=maxlen_message) :: err_msg
INTEGER :: n_flds_in, n_flds_req

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CHK_ENVGROUP_SCALAR_REAL'
! ------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr= 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check that any fields are required/ have been set for this group
IF ( .NOT. ALLOCATED(environ_fldnames_scalar_real) .OR.                        &
     SIZE(environ_fldnames_scalar_real) == 0 ) THEN
  l_fields_this_group = .FALSE.   ! Assume no fields of this type required

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check number of fields passed
n_flds_in = SIZE(flds_in)
n_flds_req = SIZE(environ_fldnames_scalar_real)
IF ( n_flds_in /= n_flds_req ) THEN
  error_code_ptr = errcode_env_field_mismatch
  WRITE(err_msg, '(A,I0,A,I0,A)') 'Wrong number of environment fields (',      &
  n_flds_in, 'for group: '//group_name//'. Expecting ', n_flds_req
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
          RoutineName, msg_out=error_message, locn_out=error_routine)

ELSE
  l_fields_this_group = .TRUE.   ! Some fields of this type required
END IF

! No other dimensions to check
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_envgroup_scalar_real

! Checks that provided fields of type 'flat integer' match the requirements
!------------------------------------------------------------------------------
SUBROUTINE check_envgroup_flat_integer(flds_in, group_name, dim_x, dim_y,      &
             error_code_ptr, l_fields_this_group, error_message, error_routine)

IMPLICIT NONE

INTEGER, INTENT(IN)    :: flds_in(:,:,:)
INTEGER, INTENT(IN) :: dim_x, dim_y
CHARACTER(LEN=*), INTENT(IN) :: group_name

INTEGER, POINTER, INTENT(IN) :: error_code_ptr
LOGICAL, INTENT(OUT) :: l_fields_this_group

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

CHARACTER(LEN=maxlen_message) :: err_msg
INTEGER :: n_flds_in, n_flds_req

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CHK_ENVGROUP_FLAT_INTEGER'
! ------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr= 0
err_msg = ''
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check that any fields are required/ have been set for this group
IF ( .NOT. ALLOCATED(environ_fldnames_flat_integer) .OR.                       &
     SIZE(environ_fldnames_flat_integer) == 0 ) THEN
  l_fields_this_group = .FALSE.   ! Assume no fields of this type required

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check number of fields passed
n_flds_in = SIZE(flds_in, DIM=3)
n_flds_req = SIZE(environ_fldnames_flat_integer)
IF ( n_flds_in /= n_flds_req ) THEN
  error_code_ptr = errcode_env_field_mismatch
  WRITE(err_msg, '(A,I0,A,I0,A)') 'Wrong number of environment fields (',      &
    n_flds_in, 'for group: '//group_name//'. Expecting ', n_flds_req
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
          RoutineName, msg_out=error_message, locn_out=error_routine)

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
ELSE
  l_fields_this_group = .TRUE.   ! Some fields of this type required
END IF

! Check first two dimensions
IF ( SIZE(flds_in, DIM=1) /= dim_x .OR. SIZE(flds_in, DIM=2) /= dim_y) THEN
  error_code_ptr = errcode_env_field_mismatch
  WRITE(err_msg, '(A,2(A,2I4))')'Incorrect dimensions for group '//group_name, &
    '. Expected ',dim_x, dim_y, '. Found ', SIZE(flds_in, DIM=1),              &
    SIZE(flds_in, DIM=2)
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
         RoutineName, msg_out=error_message, locn_out=error_routine)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_envgroup_flat_integer

! Checks that provided fields of type 'flat real' match the requirements
!------------------------------------------------------------------------------
SUBROUTINE check_envgroup_flat_real(flds_in, group_name, dim_x, dim_y,         &
             error_code_ptr, l_fields_this_group, error_message, error_routine)

IMPLICIT NONE

REAL, INTENT(IN)    :: flds_in(:,:,:)
INTEGER, INTENT(IN) :: dim_x, dim_y
CHARACTER(LEN=*), INTENT(IN) :: group_name

INTEGER, POINTER, INTENT(IN) :: error_code_ptr
LOGICAL, INTENT(OUT) :: l_fields_this_group
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

CHARACTER(LEN=maxlen_message) :: err_msg
INTEGER :: n_flds_in, n_flds_req

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CHK_ENVGROUP_FLAT_REAL'
! ------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr= 0
err_msg = ''
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check that any fields are required/ have been set for this group
IF ( .NOT. ALLOCATED(environ_fldnames_flat_real) .OR.                          &
     SIZE(environ_fldnames_flat_real) == 0 ) THEN
  l_fields_this_group = .FALSE.   ! Assume no fields of this type required

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check number of fields passed
n_flds_in = SIZE(flds_in, DIM=3)
n_flds_req = SIZE(environ_fldnames_flat_real)
IF ( n_flds_in /= n_flds_req ) THEN
  error_code_ptr = errcode_env_field_mismatch
  WRITE(err_msg, '(A,I0,A,I0,A)')'Wrong number of environment fields (',       &
    n_flds_in, 'for group: '//group_name//'. Expecting ', n_flds_req
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
         RoutineName, msg_out=error_message, locn_out=error_routine)

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
ELSE
  l_fields_this_group = .TRUE.   ! Some fields of this type required
END IF

! Check first two dimensions
IF ( SIZE(flds_in, DIM=1) /= dim_x .OR. SIZE(flds_in, DIM=2) /= dim_y) THEN
  error_code_ptr = errcode_env_field_mismatch
  WRITE(err_msg, '(A,2(A,2I4))')'Incorrect dimensions for group '//group_name, &
    '. Expected ',dim_x, dim_y, '. Found ', SIZE(flds_in, DIM=1),              &
    SIZE(flds_in, DIM=2)
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
         RoutineName, msg_out=error_message, locn_out=error_routine)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_envgroup_flat_real

! Checks that provided fields of 'fullht or fullht0 real' match requirements
!------------------------------------------------------------------------------
SUBROUTINE check_envgroup_fullht_real(flds_in, group_name, dim_x, dim_y,       &
             dim_z, error_code_ptr, l_fields_this_group, error_message,        &
             error_routine)

IMPLICIT NONE

REAL, INTENT(IN)    :: flds_in(:,:,:,:)
INTEGER, INTENT(IN) :: dim_x, dim_y, dim_z
CHARACTER(LEN=*), INTENT(IN) :: group_name

INTEGER, POINTER, INTENT(IN) :: error_code_ptr
LOGICAL, INTENT(OUT) :: l_fields_this_group

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

CHARACTER(LEN=maxlen_message) :: err_msg
INTEGER :: n_flds_in, n_flds_req, zdim

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CHK_ENVGROUP_FULLHT_REAL'
! ------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr= 0
err_msg = ''
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check that any fields are required and have been set for this group
IF ( (group_name == 'fullht_real' .AND.                                        &
       (.NOT. ALLOCATED(environ_fldnames_fullht_real) .OR.                     &
       SIZE(environ_fldnames_fullht_real) == 0) ) .OR.                         &
     (group_name == 'fullht0_real' .AND.                                       &
       (.NOT. ALLOCATED(environ_fldnames_fullht0_real) .OR.                    &
       SIZE(environ_fldnames_fullht0_real) == 0) ) ) THEN
  l_fields_this_group = .FALSE.   ! Assume no fields of this type required

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check number of fields passed
n_flds_in = SIZE(flds_in, DIM=4)

! Different expected values for 'fullht_real' and 'fullht0_real'
IF ( group_name == 'fullht_real' ) THEN
  n_flds_req = SIZE(environ_fldnames_fullht_real)
  zdim = dim_z
ELSE IF (group_name == 'fullht0_real' ) THEN
  n_flds_req = SIZE(environ_fldnames_fullht0_real)
  zdim = dim_z + 1
ELSE
  ! Unrecognised group name
  error_code_ptr = errcode_value_unknown
  WRITE(err_msg,'(A,A)')' Unknown group name provided '//group_name,           &
    'Expecting fullht_real or fullht0_real'
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
         RoutineName, msg_out=error_message, locn_out=error_routine)

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
IF ( n_flds_in /= n_flds_req ) THEN
  error_code_ptr = errcode_env_field_mismatch
  WRITE(err_msg, '(A,I0,A,I0,A)')'Wrong number of environment fields (',       &
    n_flds_in, 'for group: '//group_name//'. Expecting ', n_flds_req
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
         RoutineName, msg_out=error_message, locn_out=error_routine)

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
ELSE
  l_fields_this_group = .TRUE.   ! Some fields of this type required
END IF

! Check first three dimensions - compare levels against local zdim
IF ( SIZE(flds_in, DIM=1) /= dim_x .OR. SIZE(flds_in, DIM=2) /= dim_y .OR.     &
     SIZE(flds_in, DIM=3) /= zdim ) THEN
  error_code_ptr = errcode_env_field_mismatch
  WRITE(err_msg, '(A,2(A,3I4))')'Incorrect dimensions for group '//group_name, &
   '. Expected ',dim_x, dim_y, zdim, '. Found ', SIZE(flds_in, DIM=1),         &
    SIZE(flds_in, DIM=2), SIZE(flds_in, DIM=3)
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
         RoutineName, msg_out=error_message, locn_out=error_routine)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_envgroup_fullht_real

! Checks that provided fields of 'fullhtphot_ real' match requirements
!------------------------------------------------------------------------------
SUBROUTINE check_envgroup_fullhtphot_real(flds_in, group_name, dim_x, dim_y,   &
             dim_z, jppj, error_code_ptr, l_fields_this_group, error_message,  &
             error_routine)

IMPLICIT NONE

REAL, INTENT(IN)    :: flds_in(:,:,:,:,:)
INTEGER, INTENT(IN) :: dim_x, dim_y, dim_z, jppj
CHARACTER(LEN=*), INTENT(IN) :: group_name

INTEGER, POINTER, INTENT(IN) :: error_code_ptr
LOGICAL, INTENT(OUT) :: l_fields_this_group

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

CHARACTER(LEN=maxlen_message) :: err_msg
INTEGER :: n_flds_in, n_flds_req

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CHK_ENVGROUP_FULLHTPHOT_REAL'
! ------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr= 0
err_msg = ''
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check that any fields are required and have been set for this group
IF ( group_name == 'fullhtphot_real' .AND.                                     &
       (.NOT. ALLOCATED(environ_fldnames_fullhtphot_real) .OR.                 &
       SIZE(environ_fldnames_fullhtphot_real) == 0) ) THEN
  l_fields_this_group = .FALSE.   ! Assume no fields of this type required

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check number of fields passed
n_flds_in = SIZE(flds_in, DIM=5)

n_flds_req = SIZE(environ_fldnames_fullhtphot_real)

IF ( n_flds_in /= n_flds_req ) THEN
  error_code_ptr = errcode_env_field_mismatch
  WRITE(err_msg, '(A,I0,A,I0,A)')'Wrong number of environment fields (',       &
    n_flds_in, 'for group: '//group_name//'. Expecting ', n_flds_req
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
         RoutineName, msg_out=error_message, locn_out=error_routine)

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
ELSE
  l_fields_this_group = .TRUE.   ! Some fields of this type required
END IF

! Check first four dimensions
IF ( SIZE(flds_in, DIM=1) /= dim_x .OR. SIZE(flds_in, DIM=2) /= dim_y .OR.     &
     SIZE(flds_in, DIM=3) /= dim_z .OR. SIZE(flds_in, DIM=4) /= jppj) THEN
  error_code_ptr = errcode_env_field_mismatch
  WRITE(err_msg, '(A,2(A,4I4))')'Incorrect dimensions for group '//group_name, &
   '. Expected ',dim_x, dim_y, dim_z, jppj, '. Found ', SIZE(flds_in, DIM=1),  &
    SIZE(flds_in, DIM=2), SIZE(flds_in, DIM=3), SIZE(flds_in, DIM=4)
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
         RoutineName, msg_out=error_message, locn_out=error_routine)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_envgroup_fullhtphot_real

END MODULE photol_check_environment_mod

