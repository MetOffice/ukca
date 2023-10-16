! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module containing code for setting up a list of diagnostics
!   recognised by UKCA's diagnostic handling system and identifying
!   those that are available in the current UKCA configuration.
!
!   UKCA API procedure provided:
!
!     ukca_get_available_diag_varlist
!       - return the list of available diagnostics by name
!
!   Other public procedures provided:
!
!     init_diagnostics
!       - initialise master diagnostics list and determine availability
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

MODULE ukca_diagnostics_init_mod

USE ukca_config_specification_mod, ONLY: ukca_config_spec_type,                &
                                         l_ukca_config_available

USE ukca_diagnostics_type_mod, ONLY: dgroup_flat_real, dgroup_fullht_real

USE ukca_fieldname_mod, ONLY: maxlen_diagname

USE missing_data_mod, ONLY: imdi

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_ukca_uninit, error_report

IMPLICIT NONE

PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_DIAGNOSTICS_INIT_MOD'

! Public procedures
PUBLIC init_diagnostics, ukca_get_available_diag_varlist

! --- Module variables ---

! Info relating to available diagnostics
CHARACTER(LEN=maxlen_diagname), ALLOCATABLE, TARGET, SAVE, PUBLIC ::           &
  available_diag_varnames(:)          ! Names of available diagnostics
INTEGER, ALLOCATABLE, TARGET, SAVE ::                                          &
  available_diag_groups(:)            ! Groups nos. of available diagnostics
LOGICAL, ALLOCATABLE, TARGET, SAVE, PUBLIC ::                                  &
  l_available_diag_chem_timestep(:)   ! True if available on chem timestep only
INTEGER, ALLOCATABLE, SAVE, PUBLIC ::                                          &
  available_diag_asad_ids(:)          ! ASAD diagnostic identifiers

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE init_diagnostics(error_code_ptr, ukca_config, advt,                 &
                            error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Initialise master diagnostics list and determine availability of
!   each diagnostic given the configuration.
! ----------------------------------------------------------------------

USE ukca_diagnostics_master_mod,  ONLY: master_diag_list, bound_info,          &
                                        create_master_diagnostics_list,        &
                                        set_diagnostic_availabilities
USE asad_flux_dat,                ONLY: asad_load_default_fluxes

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr  ! Pointer to return code

TYPE(ukca_config_spec_type), INTENT(IN OUT) :: ukca_config
                                                ! UKCA configuration data

CHARACTER(LEN=*), INTENT(IN) :: advt(:)         ! Advected chemical species

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: i
INTEGER :: j
INTEGER :: n
INTEGER :: n_available

! Dr Hook

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'INIT_DIAGNOSTICS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Set up array bounds information for diagnostic groups

bound_info(dgroup_flat_real)%dim1_lower = 1
bound_info(dgroup_flat_real)%dim1_upper = ukca_config%row_length
bound_info(dgroup_flat_real)%dim2_lower = 1
bound_info(dgroup_flat_real)%dim2_upper = ukca_config%rows
bound_info(dgroup_flat_real)%dim3_lower = imdi
bound_info(dgroup_flat_real)%dim3_upper = imdi

bound_info(dgroup_fullht_real)%dim1_lower = 1
bound_info(dgroup_fullht_real)%dim1_upper = ukca_config%row_length
bound_info(dgroup_fullht_real)%dim2_lower = 1
bound_info(dgroup_fullht_real)%dim2_upper = ukca_config%rows
bound_info(dgroup_fullht_real)%dim3_lower = 1
bound_info(dgroup_fullht_real)%dim3_upper = ukca_config%model_levels

! Create the master diagnostics list
CALL create_master_diagnostics_list(error_code_ptr, error_message,             &
                                    error_routine)

IF (error_code_ptr <= 0) THEN

  ! Determine whether the ASAD diagnostic framework is applicable to the
  ! selected chemistry scheme and load the list of available ASAD chemical
  ! flux fields if it is. The list created ('asad_chemical_fluxes') will be
  ! specific to the chosen configuration settings.
  ! Note that there can be multiple entries in the ASAD chemical fluxes
  ! list for each requested diagnostic. The diagnostic output will be
  ! the sum of all those flux fields matching a given request.

  ukca_config%l_asad_chem_diags_support =                                      &
    ukca_config%l_use_gridbox_volume .AND.                                     &
    (ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat .OR.        &
     ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be)

  IF (ukca_config%l_asad_chem_diags_support) CALL asad_load_default_fluxes()

  ! Set availability of each diagnostic in the master diagnostics list given
  ! the UKCA configuration settings and count number of diagnostics available
  CALL set_diagnostic_availabilities(error_code_ptr, ukca_config, advt,        &
                                     n_available, error_message, error_routine)

END IF

IF (error_code_ptr <= 0) THEN

  ! Create list of available diagnostics and supporting information
  ALLOCATE(available_diag_varnames(n_available))
  ALLOCATE(available_diag_groups(n_available))
  ALLOCATE(l_available_diag_chem_timestep(n_available))
  ALLOCATE(available_diag_asad_ids(n_available))
  j = 0
  DO i = 1, SIZE(master_diag_list)
    IF (master_diag_list(i)%l_available) THEN
      j = j + 1
      available_diag_varnames(j) = master_diag_list(i)%varname
      available_diag_groups(j) = master_diag_list(i)%group
      l_available_diag_chem_timestep(j) = master_diag_list(i)%l_chem_timestep
      available_diag_asad_ids(j) = master_diag_list(i)%asad_id
    END IF
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE init_diagnostics

! ----------------------------------------------------------------------
SUBROUTINE ukca_get_available_diag_varlist(error_code, varnames_ptr,           &
                                           groups_ptr, l_chem_timestep_ptr,    &
                                           error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   UKCA API procedure that returns an array of field names identifying
!   diagnostics available in the current UKCA configuration.
!   Optionally, corresponding arrays are returned giving the group that
!   each diagnostic belongs to and whether or not it is only available
!   on chemistry time steps. This routine may be called at any time
!   after the call to 'ukca_setup'.
!
! Developer's note:
!   The 'error_code' can be made optional since it is not required if
!   UKCA is configured to abort on error. It remains mandatory at
!   present for consistency with other API procedures that do not yet
!   support this functionality.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

INTEGER, TARGET, INTENT(OUT) :: error_code

CHARACTER(LEN=maxlen_diagname), POINTER, INTENT(OUT) :: varnames_ptr(:)
                                   ! Pointer to list of available diagnostics

INTEGER, POINTER, OPTIONAL, INTENT(OUT) :: groups_ptr(:)
                                   ! Group of each diagnostic in the list

LOGICAL, POINTER, OPTIONAL, INTENT(OUT) :: l_chem_timestep_ptr(:)
                                   ! Indicator for each diagnostic in the list
                                   ! set to T if the diagnostic is only
                                   ! available on chemistry time steps

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER, POINTER :: error_code_ptr

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_GET_AVAILABLE_DIAG_VARLIST'

! End of header

! Use parent error code argument
error_code_ptr => error_code

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check for availability of UKCA configuration data
IF (.NOT. l_ukca_config_available) THEN
  error_code_ptr = errcode_ukca_uninit
  CALL error_report(error_code_ptr,                                            &
         'No UKCA configuration has been set up', RoutineName,                 &
         msg_out=error_message, locn_out=error_routine)
  NULLIFY(varnames_ptr)
  IF (PRESENT(groups_ptr)) NULLIFY(groups_ptr)
  IF (PRESENT(l_chem_timestep_ptr)) NULLIFY(l_chem_timestep_ptr)
  RETURN
END IF

! Assign pointers to the reference lists
varnames_ptr => available_diag_varnames
IF (PRESENT(groups_ptr)) groups_ptr => available_diag_groups
IF (PRESENT(l_chem_timestep_ptr))                                              &
  l_chem_timestep_ptr => l_available_diag_chem_timestep

RETURN
END SUBROUTINE ukca_get_available_diag_varlist

END MODULE ukca_diagnostics_init_mod
