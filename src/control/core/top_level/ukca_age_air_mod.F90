! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Increment the age-of-air tracer in the UKCA tracers array
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: Fortran
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
MODULE ukca_age_air_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_AGE_AIR_MOD'

CONTAINS

SUBROUTINE ukca_age_air(row_length, rows, model_levels, timestep,              &
                        z_top_of_model, all_tracers_names, all_tracers,        &
                        eta_theta_levels)

USE ukca_fieldname_mod,     ONLY: maxlen_fieldname, fldname_age_of_air

USE ukca_config_specification_mod, ONLY: ukca_config, i_age_reset_by_level

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,               ONLY: jprb, jpim
USE ukca_tracers_mod,       ONLY: n_tracers
USE UmPrintMgr,             ONLY: umPrint

USE yomhook,                ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN)    :: row_length
INTEGER, INTENT(IN)    :: rows
INTEGER, INTENT(IN)    :: model_levels
REAL, INTENT(IN)    :: timestep
REAL, INTENT(IN)    :: z_top_of_model ! Model top height

! List of UKCA tracer names
CHARACTER(LEN=maxlen_fieldname), INTENT(IN) :: all_tracers_names(n_tracers)

! UKCA tracer fields corresponding to the given list of tracer names
REAL, INTENT(IN OUT) :: all_tracers(row_length, rows, model_levels, n_tracers)

! Non-dimensional coordinate vector for theta levels (0.0 at planet radius,
! 1.0 at top of model), used to define level height without orography effect
! Allocatable to preserve bounds (may or may not include Level 0).
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: eta_theta_levels(:)

! Local variables

! Index of Age-of-Air tracer in all_tracers array
INTEGER, SAVE :: n_age

! Model level upto which to reset the tracer values
INTEGER, SAVE :: max_zero_age

LOGICAL, SAVE :: l_first_call=.TRUE.
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message

! loop counters
INTEGER :: i,k
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_AGE_AIR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set n_age and maximum reset level on first call
IF (l_first_call) THEN

  ! Set n_age to the index of the Age-of-Air tracer in the UKCA tracer array
  n_age = 0
  i_loop: DO i = 1, n_tracers
    IF (all_tracers_names(i) == fldname_age_of_air) THEN
      n_age = i
      EXIT i_loop
    END IF
  END DO i_loop
  IF ( n_age == 0 ) THEN
    cmessage = ' Error setting index for Age-of-Air tracer.'
    errcode = 1
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF

  ! Set maximum reset level
  max_zero_age = -1
  IF ( ukca_config%i_ageair_reset_method == i_age_reset_by_level ) THEN
    ! Check input namelist value
    IF ( ukca_config%max_ageair_reset_level > model_levels ) THEN
      WRITE(cmessage,'(A,2I8)')'Age-of-air tracer reset level > model_levels ',&
            ukca_config%max_ageair_reset_level, model_levels
      errcode = 2
      CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
    END IF
    ! Reset based on user-specified level
    max_zero_age = ukca_config%max_ageair_reset_level
  ELSE
    ! Reset based on user-specified height
    IF (.NOT. PRESENT(eta_theta_levels)) THEN
      cmessage =                                                               &
        'Missing eta_theta_levels for finding max. reset level by height.'
      errcode = 3
      CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
    END IF
    DO k = 1, model_levels-1
      IF ( (eta_theta_levels(k) * z_top_of_model)                              &
            <= ukca_config%max_ageair_reset_height ) max_zero_age = k
    END DO
  END IF   ! reset method
  IF ( max_zero_age < 0 ) THEN
    cmessage = ' Error setting maximum reset level for Age-of-Air tracer.'
    errcode = 4
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  ELSE
    WRITE(cmessage,'(2(A,I5))') 'UKCA AGE-OF-AIR: Reset method= ',             &
       ukca_config%i_ageair_reset_method,'. Tracer will be reset upto level ', &
       max_zero_age
    CALL umPrint(cmessage,src=ModuleName//':'//RoutineName)
  END IF

  l_first_call =.FALSE.

END IF   ! l_first_call

! Increment the 'age-of-air' tracer by the length of the current
! timestep. This makes this tracer equal to the time since
! the air was last in the lowest model levels
all_tracers(:,:,1:model_levels,n_age)=                                         &
    all_tracers(:,:,1:model_levels,n_age) + timestep

! set age in lower levels to zero on all timesteps
DO k = 1, max_zero_age
  all_tracers(:,:,k,n_age) = 0.0
END DO

! enforce upper boundary condition
all_tracers(:,:,model_levels,n_age) =                                          &
    all_tracers(:,:,model_levels-1,n_age)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_age_air

END MODULE ukca_age_air_mod
