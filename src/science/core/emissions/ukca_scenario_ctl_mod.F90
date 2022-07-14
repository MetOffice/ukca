! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
! Description:
!  Specify surface boundary conditions
!  - This routine will call one of three others, depending on the value
!    of i_strat_lbc_source in the UKCA configuration structure:
!      i_strat_lbc_env   - UKCA_SCENARIO_PRESCRIBED which uses external
!                          environment values
!      i_strat_lbc_wmoa1 - UKCA_SCENARIO_WMOA1 which uses the WMOA1(b) scenario
!      i_strat_lbc_rcp   - UKCA_SCENARIO_RCP which reads-in the scenario from
!         an input datafile, in a format corresponding to those which describe
!         the representative concentration pathways defined for CMIP5.
!  - It is also possible to test the UKCA_SCENARIO_RCP code by setting the
!    logical L_UKCA_TEST_SCENARIO_RCP to .TRUE. in UKCA_SCENARIO_CTL_MOD
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
!  Called from UKCA_MAIN1.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_scenario_ctl_mod
IMPLICIT NONE

PRIVATE

PUBLIC :: ukca_scenario_ctl

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_SCENARIO_CTL_MOD'

CONTAINS

SUBROUTINE ukca_scenario_ctl(n_boundary_vals, lbc_spec,                        &
                             i_year, i_day_number,                             &
                             lbc_mmr, perform_lumping)

USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim
USE ereport_mod,           ONLY: ereport
USE umPrintMgr, ONLY: umPrint, umMessage, PrStatus_Oper, PrintStatus
USE ukca_config_specification_mod, ONLY: ukca_config,                          &
                                         i_strat_lbc_wmoa1, i_strat_lbc_env
USE errormessagelength_mod, ONLY: errormessagelength

USE ukca_scenario_prescribed_mod, ONLY: ukca_scenario_prescribed
USE ukca_scenario_wmoa1_mod, ONLY: ukca_scenario_wmoa1
IMPLICIT NONE

! Number of lower BC species
INTEGER, INTENT(IN) :: n_boundary_vals
! LBC species
CHARACTER(LEN=10), INTENT(IN) :: lbc_spec(n_boundary_vals)
INTEGER, INTENT(IN) :: i_year                    ! model year
INTEGER, INTENT(IN) :: i_day_number              ! model day
! Lower BC mass mixing ratios
REAL, INTENT(IN OUT) :: lbc_mmr(n_boundary_vals)
! T for lumping of lower boundary conditions
LOGICAL, INTENT(IN) :: perform_lumping

LOGICAL, SAVE :: L_first=.TRUE.

! error handling
INTEGER :: ierr
CHARACTER(LEN=errormessagelength) :: cmessage

! loop variable for diagnostic output
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SCENARIO_CTL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE (ukca_config%i_strat_lbc_source)
CASE (i_strat_lbc_env)
  IF (l_first .AND. PrintStatus >= Prstatus_Oper) THEN
    WRITE(umMessage,'(A)')                                                     &
      'Taking UKCA Lower BC values from external environment'
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  ! take the LBC values from the driving/ parent model, but assume these numbers
  ! are only valid for the surface and that the parent model deals with the
  ! time evolution of these values
  CALL ukca_scenario_prescribed(n_boundary_vals, lbc_spec,                     &
       lbc_mmr, perform_lumping)
CASE (i_strat_lbc_wmoa1)
  IF (l_first .AND. PrintStatus >= Prstatus_Oper) THEN
    WRITE(umMessage,'(A)')'Taking UKCA Lower BC values from WMOA1 routine'
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  ! take LBC values from the WMOA1b specification (2006)
  CALL ukca_scenario_wmoa1(n_boundary_vals, lbc_spec,                          &
       i_year, i_day_number,                                                   &
       lbc_mmr, perform_lumping)
CASE DEFAULT
  ! should not have come here, so exit with an error
  WRITE(cmessage,'(A,I0)')                                                     &
    'No valid source selected for lower BC values: i_strat_lbc_source = ',     &
    ukca_config%i_strat_lbc_source
  ierr = 1
  CALL ereport(ModuleName//':'//RoutineName,ierr,cmessage)
END SELECT

! diagnostic output of the LBC values at this timestep (if requested)
IF (PrintStatus >= Prstatus_Oper) THEN
  DO i=1,n_boundary_vals
    WRITE(umMessage,'(3A,E18.8)')' UKCA Lower BC: ', TRIM(lbc_spec(i)), ' = ', &
                                 lbc_mmr(i)
    CALL umPrint(umMessage,src=RoutineName)
  END DO
END IF

L_first = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_scenario_ctl

END MODULE ukca_scenario_ctl_mod
