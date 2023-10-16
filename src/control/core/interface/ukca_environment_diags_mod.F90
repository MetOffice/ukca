! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module containing routines to handle the servicing of diagnostic
!   requests for environmental driver inputs.
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

MODULE ukca_environment_diags_mod

USE yomhook,        ONLY: lhook, dr_hook
USE parkind1,       ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_ENVIRONMENT_DIAGS_MOD'

! Public procedures
PUBLIC photol_rate_diag

CONTAINS
! ----------------------------------------------------------------------
SUBROUTINE photol_rate_diag(error_code_ptr, diag_name, i_diag_req,             &
                            photol_rates, diagnostics,                         &
                            error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Check if the given diagnostic name is a photolysis rate present in
!   the environmental driver data. If it is then service any active
!   diagnostic requests for this photolysis rate.
! ----------------------------------------------------------------------

USE ukca_diagnostics_type_mod, ONLY: diagnostics_type, n_diag_group
USE ukca_diagnostics_output_mod, ONLY: update_diagnostics_3d_real
USE ukca_fieldname_mod, ONLY: diagname_jrate_no2, diagname_jrate_o3a,          &
                              diagname_jrate_o3b, diagname_jrate_o2b
USE ukca_chem_defs_mod, ONLY: ratj_defs
USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_ukca_internal_fault, error_report

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr  ! Pointer to return code

CHARACTER(LEN=*), INTENT(IN) :: diag_name       ! Diagnostic name

INTEGER, INTENT(IN) :: i_diag_req(n_diag_group) ! Indices of requests for the
                                                ! named diagnostic

REAL, INTENT(IN) :: photol_rates(:,:,:,:)       ! Photolysis rate fields

TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics
                                                ! Diagnostic request info
                                                ! and pointers to output arrays

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: i

LOGICAL :: l_found

CHARACTER(LEN=10) :: jrate_reactant
CHARACTER(LEN=10) :: jrate_prod1
CHARACTER(LEN=10) :: jrate_prod2

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'PHOTOL_RATE_DIAG'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Get the reactant and products for identifying the reaction in the chemistry
! definitions associated with the specified photolysis rate diagnostic.
! The item no. indicated below identifies the photolysis reaction in the
! master chemistry definitions (see module ukca_chem_master_mod)
SELECT CASE (diag_name)
CASE (diagname_jrate_no2)
  ! Item no. 11
  jrate_reactant = 'NO2'
  jrate_prod1 = 'NO'
  jrate_prod2 = 'O(3P)'
CASE (diagname_jrate_o3a)
  ! Item no. 15
  jrate_reactant = 'O3'
  jrate_prod1 = 'O2'
  jrate_prod2 = 'O(1D)'
CASE (diagname_jrate_o3b)
  ! Item no. 16
  jrate_reactant = 'O3'
  jrate_prod1 = 'O2'
  jrate_prod2 = 'O(3P)'
CASE (diagname_jrate_o2b)
  ! Item no. 40
  jrate_reactant = 'O2'
  jrate_prod1 = 'O(3P)'
  jrate_prod2 = 'O(1D)'
CASE DEFAULT
  error_code_ptr = errcode_ukca_internal_fault
  CALL error_report(error_code_ptr,                                            &
    'The name given (' // TRIM(diag_name) //                                   &
    ') does not match a recognised photolysis rate diagnostic',                &
    RoutineName, msg_out=error_message, locn_out=error_routine)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END SELECT

! If the photolysis reaction is present in the chemistry definitions for
! the current UKCA configuration, the diagnostic update routine is called
! passing the relevant photolysis rate field from the environmental driver
! array. The diagnostic update routine will then service any requests for
! that photolysis rate.
i = 0
l_found = .FALSE.
DO WHILE (i < SIZE(ratj_defs) .AND. .NOT. l_found)
  i = i + 1
  IF (ratj_defs(i)%react1 == jrate_reactant .AND.                              &
      ratj_defs(i)%prod1 == jrate_prod1 .AND.                                  &
      ratj_defs(i)%prod2 == jrate_prod2) THEN
    l_found = .TRUE.
    CALL update_diagnostics_3d_real(error_code_ptr, diag_name,                 &
                                    photol_rates(:,:,:,i), diagnostics,        &
                                    i_diag_req=i_diag_req,                     &
                                    error_message=error_message,               &
                                    error_routine=error_routine)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE photol_rate_diag

END MODULE ukca_environment_diags_mod
