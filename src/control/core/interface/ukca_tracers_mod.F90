! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module for handling UKCA's tracer fields.
!
!   The module provides the following procedure for the UKCA API.
!
!     ukca_get_tracer_varlist - Returns list of names of required UKCA tracers
!
!   The following additional public procedures are provided for use within UKCA.
!
!     init_tracer_req    - Determine tracer requirement
!     tracer_copy_in_1d  - Copy tracers for a single column domain from
!                          2D parent array to internal array
!     tracer_copy_in_3d  - Copy tracers for a 3D domain from 4D parent array to
!                          internal array
!     tracer_copy_out_1d - Copy tracers from internal array to 2D parent array
!     tracer_copy_out_3d - Copy tracers from internal array to 4D parent array
!     tracer_dealloc     - Deallocate internal tracer array
!     clear_tracer_req   - Reset all tracer-related data to its initial state
!                          for a new UKCA configuration.
!
!   The module also provides a public array 'all_tracers' to hold
!   the UKCA tracer data (the internal array referred to above) and
!   a corresponding public array 'all_tracers_names' holding the names of
!   the tracers therein. The number of tracers is available as a public
!   integer 'n_tracers'.
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

MODULE ukca_tracers_mod

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

USE ukca_fieldname_mod,  ONLY: maxlen_fieldname,                               &
  fldname_nuc_sol_n,                                                           &
  fldname_nuc_sol_su,                                                          &
  fldname_ait_sol_n,                                                           &
  fldname_ait_sol_su,                                                          &
  fldname_ait_sol_bc,                                                          &
  fldname_ait_sol_om,                                                          &
  fldname_acc_sol_n,                                                           &
  fldname_acc_sol_su,                                                          &
  fldname_acc_sol_bc,                                                          &
  fldname_acc_sol_om,                                                          &
  fldname_acc_sol_ss,                                                          &
  fldname_acc_sol_du,                                                          &
  fldname_cor_sol_n,                                                           &
  fldname_cor_sol_su,                                                          &
  fldname_cor_sol_bc,                                                          &
  fldname_cor_sol_om,                                                          &
  fldname_cor_sol_ss,                                                          &
  fldname_cor_sol_du,                                                          &
  fldname_ait_ins_n,                                                           &
  fldname_ait_ins_bc,                                                          &
  fldname_ait_ins_om,                                                          &
  fldname_acc_ins_n,                                                           &
  fldname_acc_ins_du,                                                          &
  fldname_cor_ins_n,                                                           &
  fldname_cor_ins_du,                                                          &
  fldname_sup_ins_n,                                                           &
  fldname_sup_ins_du,                                                          &
  fldname_nuc_sol_om,                                                          &
  fldname_ait_sol_ss,                                                          &
  fldname_nuc_sol_so,                                                          &
  fldname_ait_sol_so,                                                          &
  fldname_acc_sol_so,                                                          &
  fldname_cor_sol_so,                                                          &
  fldname_nuc_sol_nh,                                                          &
  fldname_ait_sol_nh,                                                          &
  fldname_acc_sol_nh,                                                          &
  fldname_cor_sol_nh,                                                          &
  fldname_nuc_sol_nt,                                                          &
  fldname_ait_sol_nt,                                                          &
  fldname_acc_sol_nt,                                                          &
  fldname_cor_sol_nt,                                                          &
  fldname_acc_sol_nn,                                                          &
  fldname_cor_sol_nn,                                                          &
  fldname_passive_o3,                                                          &
  fldname_age_of_air

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_tracer_req_uninit, errcode_tracer_mismatch

IMPLICIT NONE

PRIVATE

! Public procedures
PUBLIC init_tracer_req, ukca_get_tracer_varlist, tracer_copy_in_1d,            &
       tracer_copy_in_3d, tracer_copy_out_1d, tracer_copy_out_3d,              &
       tracer_dealloc, clear_tracer_req

! Public UKCA tracer and tracer_names arrays for the current UKCA configuration
REAL, ALLOCATABLE, SAVE, PUBLIC :: all_tracers(:,:,:,:)
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE, PUBLIC ::          &
  all_tracers_names(:)

! Number of tracers required for the current UKCA configuration
INTEGER, SAVE, PUBLIC :: n_tracers = 0

! Flag to indicate whether the tracer requirement has been initialised
LOGICAL, SAVE :: l_tracer_req_available = .FALSE.

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_TRACERS_MOD'


CONTAINS


! ----------------------------------------------------------------------
SUBROUTINE init_tracer_req(config, advt,                                       &
                           error_code, error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Determines the tracers required for the current UKCA configuration.
!
! Method:
!   Create and save a reference list containing the names of the
!   required tracers. The chemistry tracers are listed first in the same
!   order as the ASAD tracer array 'advt' and the MODE tracers are
!   listed as a block immediately after the chemistry tracers. These
!   restrictions are required to ensure compatibility with other UKCA
!   modules.
! ----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: ukca_config_spec_type,                &
                                         glomap_variables

USE asad_mod, ONLY: jpctr

USE ukca_mode_setup, ONLY: nmodes,                                             &
                           mode_nuc_sol, mode_ait_sol, mode_acc_sol,           &
                           mode_cor_sol, mode_ait_insol, mode_acc_insol,       &
                           mode_cor_insol, mode_sup_insol,                     &
                           cp_su, cp_bc, cp_oc, cp_cl,                         &
                           cp_du, cp_so, cp_no3, cp_nh4, cp_nn

IMPLICIT NONE

! Subroutine arguments
TYPE(ukca_config_spec_type), INTENT(IN) :: config   ! UKCA configuration info
CHARACTER(LEN=10), INTENT(IN) :: advt(jpctr)        ! UKCA chemistry tracers
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

! Field counts
INTEGER :: n_max        ! Maximum number of tracers
INTEGER :: n            ! Count of environment fields selected

! Temporary field name array for use in collating data that will subsequently
! be copied to the field list array having the correct size allocation
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE :: fld_names(:)

! Table of mode and component associations for all recognised MODE tracers

TYPE :: mode_tr_entry
  CHARACTER(LEN=maxlen_fieldname) :: fieldname  ! Tracer name
  INTEGER :: imode                              ! Index for tracer's mode
  INTEGER :: icp                                ! Index for MMR component
                                                ! associated with tracer
                                                ! (0 for the mode's NMR tracer)
END TYPE mode_tr_entry

INTEGER, PARAMETER :: n_mode_tracers = 43       ! No. of recognised MODE tracers
TYPE(mode_tr_entry):: mode_info(n_mode_tracers) ! Table of associations

! Local variables
INTEGER :: i                                    ! Loop counter
INTEGER :: imode                                ! mode index
INTEGER :: icp                                  ! component index
LOGICAL :: l_required                           ! True if tracer is required

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'INIT_TRACER_REQ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Set up the association table for mode tracers {UM:34101-34139} to indicate
! the mode and component associated with each.

mode_info=[                                                                    &
  mode_tr_entry(fldname_nuc_sol_n, mode_nuc_sol, 0),                           &
  mode_tr_entry(fldname_nuc_sol_su, mode_nuc_sol, cp_su),                      &
  mode_tr_entry(fldname_ait_sol_n, mode_ait_sol, 0),                           &
  mode_tr_entry(fldname_ait_sol_su, mode_ait_sol, cp_su),                      &
  mode_tr_entry(fldname_ait_sol_bc, mode_ait_sol, cp_bc),                      &
  mode_tr_entry(fldname_ait_sol_om, mode_ait_sol, cp_oc),                      &
  mode_tr_entry(fldname_acc_sol_n, mode_acc_sol, 0),                           &
  mode_tr_entry(fldname_acc_sol_su, mode_acc_sol, cp_su),                      &
  mode_tr_entry(fldname_acc_sol_bc, mode_acc_sol, cp_bc),                      &
  mode_tr_entry(fldname_acc_sol_om, mode_acc_sol, cp_oc),                      &
  mode_tr_entry(fldname_acc_sol_ss, mode_acc_sol, cp_cl),                      &
  mode_tr_entry(fldname_acc_sol_du, mode_acc_sol, cp_du),                      &
  mode_tr_entry(fldname_cor_sol_n, mode_cor_sol, 0),                           &
  mode_tr_entry(fldname_cor_sol_su, mode_cor_sol, cp_su),                      &
  mode_tr_entry(fldname_cor_sol_bc, mode_cor_sol, cp_bc),                      &
  mode_tr_entry(fldname_cor_sol_om, mode_cor_sol, cp_oc),                      &
  mode_tr_entry(fldname_cor_sol_ss, mode_cor_sol, cp_cl),                      &
  mode_tr_entry(fldname_cor_sol_du, mode_cor_sol, cp_du),                      &
  mode_tr_entry(fldname_ait_ins_n, mode_ait_insol, 0),                         &
  mode_tr_entry(fldname_ait_ins_bc, mode_ait_insol, cp_bc),                    &
  mode_tr_entry(fldname_ait_ins_om, mode_ait_insol, cp_oc),                    &
  mode_tr_entry(fldname_acc_ins_n, mode_acc_insol, 0),                         &
  mode_tr_entry(fldname_acc_ins_du, mode_acc_insol, cp_du),                    &
  mode_tr_entry(fldname_cor_ins_n, mode_cor_insol, 0),                         &
  mode_tr_entry(fldname_cor_ins_du, mode_cor_insol, cp_du),                    &
  mode_tr_entry(fldname_nuc_sol_om, mode_nuc_sol, cp_oc),                      &
  mode_tr_entry(fldname_ait_sol_ss, mode_ait_sol, cp_cl),                      &
  mode_tr_entry(fldname_nuc_sol_so, mode_nuc_sol, cp_so),                      &
  mode_tr_entry(fldname_ait_sol_so, mode_ait_sol, cp_so),                      &
  mode_tr_entry(fldname_acc_sol_so, mode_acc_sol, cp_so),                      &
  mode_tr_entry(fldname_cor_sol_so, mode_cor_sol, cp_so),                      &
  mode_tr_entry(fldname_nuc_sol_nh, mode_nuc_sol, cp_nh4),                     &
  mode_tr_entry(fldname_ait_sol_nh, mode_ait_sol, cp_nh4),                     &
  mode_tr_entry(fldname_acc_sol_nh, mode_acc_sol, cp_nh4),                     &
  mode_tr_entry(fldname_cor_sol_nh, mode_cor_sol, cp_nh4),                     &
  mode_tr_entry(fldname_nuc_sol_nt, mode_nuc_sol, cp_no3),                     &
  mode_tr_entry(fldname_ait_sol_nt, mode_ait_sol, cp_no3),                     &
  mode_tr_entry(fldname_acc_sol_nt, mode_acc_sol, cp_no3),                     &
  mode_tr_entry(fldname_cor_sol_nt, mode_cor_sol, cp_no3),                     &
  mode_tr_entry(fldname_acc_sol_nn, mode_acc_sol, cp_nn),                      &
  mode_tr_entry(fldname_cor_sol_nn, mode_cor_sol, cp_nn),                      &
  mode_tr_entry(fldname_sup_ins_n, mode_sup_insol, 0),                         &
  mode_tr_entry(fldname_sup_ins_du, mode_sup_insol, cp_du)                     &
]

! Ensure all tracer-related data are in uninitialised state
IF (l_tracer_req_available) CALL clear_tracer_req()

! Allocate temporary name array (allow for chemistry, mode, 'Passive O3' and
! 'Age of Air' tracers)
n_max = jpctr + n_mode_tracers + 2
ALLOCATE(fld_names(n_max))

! Add all tracers from advt (may include H2O)
n = jpctr
IF (jpctr > 0) fld_names(1:n) = advt(:)

! For each MODE tracer, the tracer is required if it's associated mode is
! active and either it is the number mixing ratio tracer for that mode
! (indicated by icp=0) or it's associated component (expressed as a mass mixing
! ratio) is active for that mode.

DO i = 1, n_mode_tracers

  l_required = .FALSE.
  imode = mode_info(i)%imode
  icp = mode_info(i)%icp
  ! Check whether mode is active (assume inactive if status is unspecified)
  IF (imode <= nmodes) THEN
    IF (glomap_variables%mode(imode)) THEN
      IF (icp == 0) THEN
        l_required = .TRUE. ! NMR for an active mode
      ELSE
        ! Check whether component is active for this mode
        ! (assume inactive if status is unspecified)
        IF (icp <= glomap_variables%ncp) THEN
          IF (glomap_variables%component(imode,icp)) l_required = .TRUE.
                            ! MMR for an active mode/component combination
        END IF
      END IF
    END IF
  END IF

  IF (l_required) THEN
    n = n + 1
    IF (n <= n_max) fld_names(n) = mode_info(i)%fieldname
  END IF

END DO

! Add 'Passive O3' tracer if required
IF (config%l_ukca_strattrop .OR. config%l_ukca_strat) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_passive_o3
END IF

! Add 'Age of Air' tracer if required
IF (config%l_ukca_ageair) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_age_of_air
END IF

! Check number of fields required against maximum
IF (n > n_max) THEN
  error_code = errcode_tracer_req_uninit
  IF (PRESENT(error_message)) WRITE(error_message,'(A,I0,A,I0)')               &
    'Number of required tracers (', n,                                         &
    ') exceeds maximum: n_max = ', n_max
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Create reference list of required tracers
ALLOCATE(all_tracers_names(n))
all_tracers_names = fld_names(1:n)
n_tracers = n

! Set public flag to show availability of tracer requirement
l_tracer_req_available = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_tracer_req


! ----------------------------------------------------------------------
SUBROUTINE ukca_get_tracer_varlist(varnames_ptr, error_code,                   &
                                   error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   UKCA API procedure that returns a list of field names identifying
!   the tracers required for the current UKCA configuration.
!
! Method:
!   Return pointer to the reference list giving the names of required
!   tracers.
!   A non-zero error code is returned if the requirement for the current
!   UKCA configuration has not been initialised and the pointer will be
!   disassociated.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=maxlen_fieldname), POINTER, INTENT(OUT) :: varnames_ptr(:)
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_GET_TRACER_VARLIST'

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of tracer requirement
IF (.NOT. l_tracer_req_available) THEN
  error_code = errcode_tracer_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Tracer requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  NULLIFY(varnames_ptr)
  RETURN
END IF

! Assign pointer to the reference list
varnames_ptr => all_tracers_names

RETURN
END SUBROUTINE ukca_get_tracer_varlist


! ----------------------------------------------------------------------
SUBROUTINE tracer_copy_in_1d(error_code_ptr, tracer_data, model_levels,        &
                             error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
! Copies tracer data for a single column domain from the given
! tracer_data array to the all_tracers array.
! Bounds for the vertical array dimension are checked for validity based
! on the UKCA extent i.e. 1:model_levels.
! The input fields must span the UKCA extent but may optionally extend
! beyond it in either direction. The vertical extent of the data copied
! is restricted to the UKCA bounds.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr ! Return status code

! Input tracer field array. Dimensions: Z,N
! where Z is no. of levels in tracer fields
!       N is number of tracers
REAL, ALLOCATABLE, INTENT(IN) :: tracer_data(:,:)

! Dimension of UKCA domain
INTEGER, INTENT(IN) :: model_levels         ! Size of UKCA z dimension

! Error information
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
                                            ! Return error message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
                                            ! Routine name where error trapped

! Local variables
! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'TRACER_COPY_IN_1D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of tracer requirement
IF (.NOT. l_tracer_req_available) THEN
  error_code_ptr = errcode_tracer_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Tracer requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check the data field array contains the expected number of tracers
IF (SIZE(tracer_data,DIM=2) /= n_tracers) THEN
  error_code_ptr = errcode_tracer_mismatch
  IF (PRESENT(error_message)) WRITE(error_message,'(A,I0,A,I0)')               &
    'Number of tracer fields (', SIZE(tracer_data,DIM=2),                      &
    ') does not match requirement: n_tracers = ', n_tracers
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (n_tracers > 0) THEN

  ! Check array bounds for data fields are compatible with UKCA.
  ! The field data supplied must fill the required domain but may extend
  ! beyond it to avoid the need for pre-trimming by the parent model.
  IF (LBOUND(tracer_data,DIM=1) > 1 .OR.                                       &
      UBOUND(tracer_data,DIM=1) < model_levels)                                &
  THEN
    error_code_ptr = errcode_tracer_mismatch
    IF (PRESENT(error_message))                                                &
      error_message =                                                          &
        'The tracer fields have one or more invalid array bounds'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  ! Copy required data.
  ! Any data outside the required bounds (e.g. halos) are discarded.
  ALLOCATE(all_tracers(1,1,1:model_levels,1:n_tracers))
  all_tracers(1,1,:,:) = tracer_data(1:model_levels,:)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE tracer_copy_in_1d


! ----------------------------------------------------------------------
SUBROUTINE tracer_copy_in_3d(error_code_ptr, tracer_data, row_length,          &
                             rows, model_levels, error_message,                &
                             error_routine)
! ----------------------------------------------------------------------
! Description:
! Copies tracer data for a 3D domain from the given tracer_data array
! to the all_tracers array.
! Bounds are checked for validity based on the UKCA horizontal and
! vertical extents i.e 1:row_length, 1:rows and 1:model_levels.
! The input fields must span the UKCA extents but may optionally extend
! beyond them. The horizontal and vertical extents of the data copied
! are restricted to the UKCA bounds.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr ! Return status code

! Input tracer field array. Dimensions: X,Y,Z,N
! where X is row length of tracer field
!       Y is no. of rows in tracer field
!       Z is no. of levels in tracer fields
!       N is number of tracers
REAL, ALLOCATABLE, INTENT(IN) :: tracer_data(:,:,:,:)

! Dimensions of UKCA domain
INTEGER, INTENT(IN) :: row_length           ! Size of UKCA x dimension
INTEGER, INTENT(IN) :: rows                 ! Size of UKCA y dimension
INTEGER, INTENT(IN) :: model_levels         ! Size of UKCA z dimension

! Error information
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
                                            ! Return error message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
                                            ! Routine name where error trapped

! Local variables
! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'TRACER_COPY_IN_3D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of tracer requirement
IF (.NOT. l_tracer_req_available) THEN
  error_code_ptr = errcode_tracer_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Tracer requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check the data field array contains the expected number of tracers
IF (SIZE(tracer_data,DIM=4) /= n_tracers) THEN
  error_code_ptr = errcode_tracer_mismatch
  IF (PRESENT(error_message)) WRITE(error_message,'(A,I0,A,I0)')               &
    'Number of tracer fields (', SIZE(tracer_data,DIM=4),                      &
    ') does not match requirement: n_tracers = ', n_tracers
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (n_tracers > 0) THEN

  ! Check array bounds for data fields are compatible with UKCA.
  ! The field data supplied must fill the required domain but may extend
  ! beyond it to avoid the need for pre-trimming (e.g. halo removal) by
  ! the parent model.
  IF (LBOUND(tracer_data,DIM=1) > 1 .OR.                                       &
      UBOUND(tracer_data,DIM=1) < row_length .OR.                              &
      LBOUND(tracer_data,DIM=2) > 1 .OR.                                       &
      UBOUND(tracer_data,DIM=2) < rows .OR.                                    &
      LBOUND(tracer_data,DIM=3) > 1 .OR.                                       &
      UBOUND(tracer_data,DIM=3) < model_levels)                                &
  THEN
    error_code_ptr = errcode_tracer_mismatch
    IF (PRESENT(error_message))                                                &
      error_message =                                                          &
        'The tracer fields have one or more invalid array bounds'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  ! Copy required data.
  ! Any data outside the required bounds (e.g. halos) are discarded.
  ALLOCATE(all_tracers(1:row_length,1:rows,1:model_levels,1:n_tracers))
  all_tracers = tracer_data(1:row_length,1:rows,1:model_levels,:)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE tracer_copy_in_3d


! ----------------------------------------------------------------------
SUBROUTINE tracer_copy_out_1d(model_levels, tracer_data)
! ----------------------------------------------------------------------
! Description:
! Copies the tracer data from the all_tracers array to the given
! tracer_data array for a single column domain.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN) :: model_levels        ! Size of UKCA z dimension

! Tracer field array to update. Dimensions: Z,N
! where Z is no. of levels in tracer fields
!       N is number of tracers
REAL, ALLOCATABLE, INTENT(IN OUT) :: tracer_data(:,:)

! Local variables
! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'TRACER_COPY_OUT_1D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

tracer_data(1:model_levels,:) = all_tracers(1,1,:,:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE tracer_copy_out_1d


! ----------------------------------------------------------------------
SUBROUTINE tracer_copy_out_3d(row_length, rows, model_levels, tracer_data)
! ----------------------------------------------------------------------
! Description:
! Copies the tracer data from the all_tracers array to the given
! tracer_data array for a 3D domain.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: row_length          ! Size of UKCA x dimension
INTEGER, INTENT(IN) :: rows                ! Size of UKCA y dimension
INTEGER, INTENT(IN) :: model_levels        ! Size of UKCA z dimension

! Tracer field array to update. Dimensions: X,Y,Z,N
! where X is row length of tracer field
!       Y is no. of rows in tracer field
!       Z is no. of levels in tracer fields
!       N is number of tracers
REAL, ALLOCATABLE, INTENT(IN OUT) :: tracer_data(:,:,:,:)

! Local variables
! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'TRACER_COPY_OUT_3D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

tracer_data(1:row_length,1:rows,1:model_levels,:) = all_tracers(:,:,:,:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE tracer_copy_out_3d


! ----------------------------------------------------------------------
SUBROUTINE tracer_dealloc()
! ----------------------------------------------------------------------
! Description:
! Deallocates the all_tracers array.
! ----------------------------------------------------------------------

IMPLICIT NONE

! CHARACTER(LEN=*), PARAMETER :: RoutineName = 'TRACER_DEALLOC'

IF (ALLOCATED(all_tracers)) DEALLOCATE(all_tracers)

END SUBROUTINE tracer_dealloc


! ----------------------------------------------------------------------
SUBROUTINE clear_tracer_req()
! ----------------------------------------------------------------------
! Description:
!   Resets all tracer-related data to its initial state for a new
!   UKCA configuration.
!
! Method:
!   Ensure the UKCA tracer data is cleared.
!   Deallocate required tracers array and reset flag showing
!   availability status of tracer requirement.
! ----------------------------------------------------------------------

IMPLICIT NONE

! CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CLEAR_TRACER_REQ'

! Clear any UKCA tracer data
IF (ALLOCATED(all_tracers)) DEALLOCATE(all_tracers)

! Clear tracer requirement data
n_tracers = 0
IF (ALLOCATED(all_tracers_names)) DEALLOCATE(all_tracers_names)
l_tracer_req_available = .FALSE.

END SUBROUTINE clear_tracer_req

END MODULE ukca_tracers_mod
