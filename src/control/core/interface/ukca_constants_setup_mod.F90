! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module providing UKCA API routine to set up UKCA's configurable constants
!   i.e. those that can be overriden by parent values.
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

MODULE ukca_constants_setup_mod

IMPLICIT NONE
PRIVATE

PUBLIC ukca_constants_setup

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_CONSTANTS_SETUP_MOD'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE ukca_constants_setup(error_code,                                    &
                                const_rmol,                                    &
                                const_tfs,                                     &
                                const_rho_water,                               &
                                const_rhosea,                                  &
                                const_lc,                                      &
                                const_avogadro,                                &
                                const_boltzmann,                               &
                                const_rho_so4,                                 &
                                error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!
!   Initialise UKCA's configurable constants to their default values and
!   override any of these with alternatives values if present in the
!   argument list.
!
! Developer's note:
!   The 'error_code' argument can be made optional since it is not required
!   if UKCA is configured to abort on error. It remains mandatory at
!   present for consistency with other API procedures that do not yet
!   support this functionality.
! ----------------------------------------------------------------------

USE ukca_config_constants_mod, ONLY: init_config_constants,                    &
                                     rmol,                                     &
                                     tfs,                                      &
                                     rho_water,                                &
                                     rhosea,                                   &
                                     lc,                                       &
                                     avogadro,                                 &
                                     boltzmann,                                &
                                     rho_so4,                                  &
                                     l_ukca_constants_available

USE ukca_config_specification_mod, ONLY: ukca_config

USE ukca_error_mod, ONLY: error_report, maxlen_message, maxlen_procname,       &
                          errcode_unexpected_api_call

USE parkind1, ONLY: jpim, jprb      ! DrHook
USE yomhook,  ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Subroutine arguments.

! Optional arguments are used to override values of configurable constants
! defined in 'ukca_config_constants_mod'.

INTEGER, TARGET, INTENT(OUT) :: error_code

REAL, OPTIONAL, INTENT(IN) :: const_rmol
REAL, OPTIONAL, INTENT(IN) :: const_tfs
REAL, OPTIONAL, INTENT(IN) :: const_rho_water
REAL, OPTIONAL, INTENT(IN) :: const_rhosea
REAL, OPTIONAL, INTENT(IN) :: const_lc
REAL, OPTIONAL, INTENT(IN) :: const_avogadro
REAL, OPTIONAL, INTENT(IN) :: const_boltzmann
REAL, OPTIONAL, INTENT(IN) :: const_rho_so4

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER, POINTER :: error_code_ptr

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CONSTANTS_SETUP'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Use parent-supplied argument for return code
error_code_ptr => error_code

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check for potential conflict with a previous call (possibly via 'ukca_setup')
IF (l_ukca_constants_available) THEN
  error_code_ptr = errcode_unexpected_api_call
  CALL error_report(ukca_config%i_error_method, error_code_ptr,                &
         'UKCA configurable constants have already been set', RoutineName,     &
         msg_out=error_message, locn_out=error_routine)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Set constants to default values
CALL init_config_constants()

! Override default values with any values provided in argument list

IF (PRESENT(const_rmol))                                                       &
  rmol = const_rmol
IF (PRESENT(const_tfs))                                                        &
  tfs = const_tfs
IF (PRESENT(const_rho_water))                                                  &
  rho_water = const_rho_water
IF (PRESENT(const_rhosea))                                                     &
  rhosea = const_rhosea
IF (PRESENT(const_lc))                                                         &
  lc = const_lc
IF (PRESENT(const_avogadro))                                                   &
  avogadro = const_avogadro
IF (PRESENT(const_boltzmann))                                                  &
  boltzmann = const_boltzmann
IF (PRESENT(const_rho_so4))                                                    &
  rho_so4 = const_rho_so4

l_ukca_constants_available = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE ukca_constants_setup

END MODULE ukca_constants_setup_mod
