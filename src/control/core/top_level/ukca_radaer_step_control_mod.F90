! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module for controlling the UKCA radaer time step in a parent application
!   where all data transfer between the parent and UKCA during a run is
!   restricted to a single API call for each time step.
!   This is a necessary condition for running multiple UKCA radaer time step
!   calls in parallel on the same node with shared memory (one time step call
!   for each column of a parent model domain).
!
!   The module provides the following procedure for the UKCA API.
!
!     ukca_radaer_step_control - Obtain all required environmental driver data
!                                and perform one UKCA radaer time step for a
!                                single column.
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

MODULE ukca_radaer_step_control_mod

! Dr Hook modules
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC ukca_radaer_step_control

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_RADAER_STEP_CONTROL_MOD'

! Generic interface for UKCA time step subroutine - overloaded according to
! the dimension of the tracer and NTP data
INTERFACE ukca_radaer_step_control
  MODULE PROCEDURE ukca_radaer_step_control_1d_domain ! Column mode
  ! Include a full domain mode later
END INTERFACE ukca_radaer_step_control

CONTAINS

! ----------------------------------------------------------------------

SUBROUTINE ukca_radaer_step_control_1d_domain(                                 &
                                          timestep_number,                     &
                                          ! Scalar environment field groups
                                          envgroup_flat_integer,               &
                                          envgroup_scalar_real,                &
                                          envgroup_flat_real,                  &
                                          emissions_flat,                      &
                                          envgroup_flat_logical,               &
                                          ! 1D environment field groups
                                          envgroup_fullht_real,                &
                                          envgroup_fullht0_real,               &
                                          ! Diagnostics output
                                          diag_status_fullht_real,             &
                                          diag_data_fullht_real,               &
                                          ! Error return info
                                          error_message,                       &
                                          error_routine )
! ----------------------------------------------------------------------
! Description:
!   Obtains environmental driver and emission fields and performs
!   one UKCA radaer time step for a single column.
!
! Method:
!   Environmental driver fields are grouped according to their
!   dimensionality, size and type.
!   The procedure for handling the UKCA radaer time step is:
!   1) Set up environmental drivers by calling 'ukca_radaer_set_environment'
!   for each field in each group.
!   2) Execute the time step.
!   Processing of environmental drivers by groups allows the same arguments
!   to be used for different sets of input fields so that the same
!   'ukca_radaer_step_control' call can be used for different configurations.
!   It also reduces the number of arguments that need to be processed
!   (fields do not need to be processed individually).
!
!   N.B. Handling of environmental drivers does not yet
!   support multi-thread calls. Local copies of the driver
!   data will be required for thread-safe operation.
! ----------------------------------------------------------------------
















! Do the time step
IF (error_code <= 0) THEN

  CALL ukca_radaer_step(timestep_number, current_time,                         &
                 tracer_data, ntp_data, r_theta_levels, r_rho_levels,          &
                 error_code, previous_time=previous_time,                      &
                 eta_theta_levels=eta_theta_levels,                            &
                 diag_status_flat_real=diag_status_flat_real,                  &
                 diag_status_fullht_real=diag_status_fullht_real,              &
                 diag_data_flat_real=diag_data_flat_real,                      &
                 diag_data_fullht_real=diag_data_fullht_real,                  &
                 error_message=error_message, error_routine=error_routine)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_radaer_step_control_1d_domain

END MODULE ukca_radaer_step_control_mod
