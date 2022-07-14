! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module for handling the UKCA time step
!
!   The module provides the following procedure for the UKCA API.
!
!     ukca_step - Perform one UKCA time step (overloaded for domains
!                 with different domain dimensions)
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

MODULE ukca_step_mod

USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_tracers_mod, ONLY: tracer_copy_in, tracer_copy_out, tracer_dealloc,   &
                            all_tracers
USE ukca_ntp_mod, ONLY: ntp_copy_in, ntp_copy_out, ntp_dealloc, all_ntp
USE ukca_main1_mod, ONLY: ukca_main1
USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname

! Dr Hook modules
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC ukca_step

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_STEP_MOD'

! Generic interface for UKCA time step subroutine - overloaded according to
! the dimension of the tracer and NTP data
INTERFACE ukca_step
  MODULE PROCEDURE ukca_step_1d_domain
  MODULE PROCEDURE ukca_step_3d_domain
END INTERFACE ukca_step


CONTAINS


! ----------------------------------------------------------------------
SUBROUTINE ukca_step_1d_domain(timestep_number, current_time,                  &
                               tracer_data_parent, ntp_data_parent,            &
                               error_code, previous_time,                      &
                               error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Variant of the UKCA API generic procedure ukca_step.
!   Performs one UKCA time step.
!   Each input tracer and NTP field is defined for a single column.
!
! Method:
!   1) Copy tracer and NTP data from the parent-supplied 2D arrays to
!   UKCA's internal data structures, ignoring any data outside the
!   required domain.
!   2) Perform the UKCA time step
!   3) Copy the updated tracers and NTPs back to the parent arrays.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

! Model timestep number (counted from basis time at start of run)
INTEGER, INTENT(IN) :: timestep_number

! Current model time (year, month, day, hour, minute, second, day of year)
INTEGER, INTENT(IN) :: current_time(7)

! UKCA tracers from the parent model. Dimensions: Z,N
! where Z is no. of levels in tracer fields
!       N is number of tracers
REAL, ALLOCATABLE, INTENT(IN OUT) :: tracer_data_parent(:, :)

! Non-transported prognostics from the parent model. Dimensions: Z,N
REAL, ALLOCATABLE, INTENT(IN OUT) :: ntp_data_parent(:, :)

! Error code for status reporting
INTEGER, INTENT(OUT) :: error_code

! Model time at previous timestep (required for chemistry)
INTEGER, OPTIONAL, INTENT(IN) :: previous_time(7)

! Further arguments for status reporting
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_STEP_1D_DOMAIN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Populate the all_tracers and all_ntp arrays defining the UKCA state using
! the corresponding data from the parent model
CALL tracer_copy_in(tracer_data_parent, ukca_config%model_levels,              &
                    error_code, error_message=error_message,                   &
                    error_routine=error_routine)
IF (error_code > 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
CALL ntp_copy_in(ntp_data_parent,ukca_config%model_levels,                     &
                 error_code, error_message=error_message,                      &
                 error_routine=error_routine)
IF (error_code > 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Do the time step
CALL ukca_main1(timestep_number, current_time, all_tracers, all_ntp,           &
                error_code, previous_time=previous_time,                       &
                error_message=error_message, error_routine=error_routine)

! Update the tracer_data_parent and ntp_data_parent arrays from the UKCA state
! at the end of the time step. These are then passed back to the parent.
CALL tracer_copy_out(ukca_config%model_levels, tracer_data_parent)
CALL ntp_copy_out(ukca_config%model_levels, ntp_data_parent)

! Clear UKCA state data ready for next time step
CALL tracer_dealloc()
CALL ntp_dealloc()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_step_1d_domain

! ----------------------------------------------------------------------
SUBROUTINE ukca_step_3d_domain(timestep_number, current_time,                  &
                               tracer_data_parent, ntp_data_parent,            &
                               error_code, previous_time,                      &
                               error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Variant of the UKCA API generic procedure ukca_step.
!   Performs one UKCA time step.
!   Each input tracer and NTP field is defined on a 3D grid.
!
! Method:
!   1) Copy tracer and NTP data from the parent-supplied 4D arrays to
!   UKCA's internal data structures, ignoring any data outside the
!   required domain.
!   2) Perform the UKCA time step
!   3) Copy the updated tracers and NTPs back to the parent arrays.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

! Model timestep number (counted from basis time at start of run)
INTEGER, INTENT(IN) :: timestep_number

! Current model time (year, month, day, hour, minute, second, day of year)
INTEGER, INTENT(IN) :: current_time(7)

! UKCA tracers from the parent model. Dimensions: X,Y,Z,N
! where X is row length of tracer field (= no. of columns)
!       Y is no. of rows in tracer field
!       Z is no. of levels in tracer fields
!       N is number of tracers
REAL, ALLOCATABLE, INTENT(IN OUT) :: tracer_data_parent(:, :, :, :)

! Non-transported prognostics from the parent model. Dimensions: X,Y,Z,N
REAL, ALLOCATABLE, INTENT(IN OUT) :: ntp_data_parent(:, :, :, :)

! Error code for status reporting
INTEGER, INTENT(OUT) :: error_code

! Model time at previous timestep (required for chemistry)
INTEGER, OPTIONAL, INTENT(IN) :: previous_time(7)

! Further arguments for status reporting
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_STEP_3D_DOMAIN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Populate the all_tracers and all_ntp arrays defining the UKCA state using
! the corresponding data from the parent model
CALL tracer_copy_in(tracer_data_parent, ukca_config%row_length,                &
                    ukca_config%rows, ukca_config%model_levels,                &
                    error_code, error_message=error_message,                   &
                    error_routine=error_routine)
IF (error_code > 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
CALL ntp_copy_in(ntp_data_parent, ukca_config%row_length,                      &
                 ukca_config%rows, ukca_config%model_levels,                   &
                 error_code, error_message=error_message,                      &
                 error_routine=error_routine)
IF (error_code > 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Do the time step
CALL ukca_main1(timestep_number, current_time, all_tracers, all_ntp,           &
                error_code, previous_time, error_message, error_routine)

! Update the tracer_data_parent and ntp_data_parent arrays from the UKCA state
! at the end of the time step. These are then passed back to the parent.
CALL tracer_copy_out(ukca_config%row_length, ukca_config%rows,                 &
                     ukca_config%model_levels, tracer_data_parent)
CALL ntp_copy_out(ukca_config%row_length, ukca_config%rows,                    &
                  ukca_config%model_levels, ntp_data_parent)

! Clear UKCA state data ready for next time step
CALL tracer_dealloc()
CALL ntp_dealloc()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_step_3d_domain

END MODULE ukca_step_mod
