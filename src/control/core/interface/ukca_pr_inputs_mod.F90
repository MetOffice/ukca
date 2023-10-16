! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!  Print inputs to UKCA to help with debugging
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  Fortran 2003
!
! ----------------------------------------------------------------------
MODULE ukca_pr_inputs_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: ukca_pr_inputs

CONTAINS

SUBROUTINE ukca_pr_inputs(error_code_ptr, timestep_number, environ_ptrs,       &
                          land_fraction, thick_bl_levels, t_theta_levels,      &
                          rel_humid_frac, z_half, error_message, error_routine)

USE ukca_um_legacy_mod,   ONLY: mype
USE umPrintMgr,           ONLY: umMessage, UmPrint, PrintStatus,               &
                                PrStatus_Diag, PrStatus_Oper
USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_environment_req_mod,      ONLY: print_environment_summary
USE ukca_environment_fields_mod,   ONLY: env_field_ptrs_type
USE ukca_ntp_mod,         ONLY: print_all_ntp
USE ukca_error_mod,       ONLY: maxlen_message, maxlen_procname

IMPLICIT NONE

INTEGER, POINTER, INTENT(IN) :: error_code_ptr
INTEGER, INTENT(IN) :: timestep_number

TYPE(env_field_ptrs_type), INTENT(IN) :: environ_ptrs(:)

REAL, INTENT(IN) :: land_fraction(ukca_config%row_length, ukca_config%rows)
REAL, ALLOCATABLE, INTENT(IN) :: thick_bl_levels(:,:,:)
REAL, ALLOCATABLE, INTENT(IN) :: t_theta_levels(:,:,:)
REAL, ALLOCATABLE, INTENT(IN) :: rel_humid_frac(:,:,:)
REAL, ALLOCATABLE, INTENT(IN) :: z_half(:,:,:)

! Arguments for error reporting
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

INTEGER :: k ! loop counter

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PR_INPUTS'

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Debug Non-transported prognostics every TS, depending on PrintStatus

IF (PrintStatus >= PrStatus_Diag) THEN

  ! print out info on NTP
  CALL print_all_ntp

END IF

IF (PrintStatus >= PrStatus_Oper .AND.                                         &
    timestep_number == ukca_config%env_log_step) THEN

  WRITE(umMessage,'(A)') '================================================='
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I0)') 'MAX and MIN of UKCA INPUTS at Timestep ',         &
                            timestep_number
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A)') '================================================='
  CALL umPrint(umMessage,src=RoutineName)

  ! Print summary stats for environmental drivers

  CALL print_environment_summary(error_code_ptr, mype, environ_ptrs,           &
                                 error_message, error_routine)
  IF (error_code_ptr > 0) RETURN

  ! Print some additional stats for derived values

  IF (ALLOCATED(thick_bl_levels)) THEN
    WRITE(umMessage,'(A,I8,2E12.4)') 'Thick_bl_levels (level 1): ', mype,      &
      MAXVAL(thick_bl_levels(:,:,1)), MINVAL(thick_bl_levels(:,:,1))
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  IF (ALLOCATED(t_theta_levels)) THEN
    DO k = 1, ukca_config%model_levels
      WRITE(umMessage,'(A,I6)') 'LEVEL: ',k
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I8,I6,2E12.4)') 't_theta_levels: ', mype, k,         &
        MAXVAL(t_theta_levels(:,:,k)), MINVAL(t_theta_levels(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF

  IF (ALLOCATED(rel_humid_frac)) THEN
    DO k = 1, ukca_config%model_levels
      WRITE(umMessage,'(A,I6)') 'WET LEVEL: ',k
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I8,I6,2E12.4)') 'rel_humid_frac: ', mype, k,         &
        MAXVAL(rel_humid_frac(:,:,k)), MINVAL(rel_humid_frac(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF

  IF (ALLOCATED(z_half)) THEN
    DO k = 1, ukca_config%bl_levels
      WRITE(umMessage,'(A,I6)') 'BL LEVEL: ',k
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I8,I6,2E12.4)') 'z_half:     ', mype, k,             &
        MAXVAL(z_half(:,:,k)), MINVAL(z_half(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF

  WRITE(umMessage,'(A16,I5,2E12.4)') 'land_fraction: ', mype,                  &
    MAXVAL(land_fraction), MINVAL(land_fraction)

END IF  ! PrintStatus >= PrStatus_Oper .AND. timestep_number == ...'

END SUBROUTINE ukca_pr_inputs

END MODULE ukca_pr_inputs_mod
