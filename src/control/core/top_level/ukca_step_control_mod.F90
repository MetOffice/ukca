! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module for controlling the UKCA time step in a parent application
!   where all data transfer between the parent and UKCA during a run is
!   restricted to a single API call for each time step.
!   This is a necessary condition for running multiple UKCA time step calls in
!   parallel on the same node with shared memory (one time step call for
!   each column of a parent model domain).
!   A further condition, not currently met, is that all fields having a
!   horizontal dimension in the parent must have a separate instance in
!   each time step call, so must be passed by argument within UKCA instead
!   of being shared module variables. Further work is planned to implement
!   this.
!
!   The module provides the following procedure for the UKCA API.
!
!     ukca_step_control - Obtain all required environmental driver data
!                         and perform one UKCA time step for a single column.
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

MODULE ukca_step_control_mod

! Dr Hook modules
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC ukca_step_control

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_STEP_CONTROL_MOD'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE ukca_step_control(timestep_number, current_time,                    &
                             tracer_data, ntp_data,                            &
                             error_code, previous_time,                        &
                             ! Scalar environment field groups
                             envgroup_flat_integer,                            &
                             envgroup_scalar_real,                             &
                             envgroup_flat_real,                               &
                             emissions_flat,                                   &
                             envgroup_flat_logical,                            &
                             ! 1D environment field groups
                             envgroup_flatpft_real,                            &
                             envgroup_fullht_real,                             &
                             envgroup_fullht0_real,                            &
                             envgroup_fullhtp1_real,                           &
                             envgroup_bllev_real,                              &
                             envgroup_entlev_real,                             &
                             envgroup_land_real,                               &
                             emissions_fullht,                                 &
                             ! 2D environment field groups
                             envgroup_landtile_real,                           &
                             envgroup_landpft_real,                            &
                             envgroup_landtile_logical,                        &
                             ! Error return info
                             error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Obtains environmental driver and emission fields and performs
!   one UKCA time step for a single column.
!
! Method:
!   Environmental driver fields are grouped according to their
!   dimensionality, size and type. Emissions fields are treated as a
!   special case of environmental drivers and are grouped separately
!   as they require different treatment.
!   The procedure for handling the UKCA time step is:
!   1) Set up environmental drivers are setup by calling 'ukca_set_emission'
!   for each field in each group.
!   2) Set up emissions by calling 'ukca_set_emission' for each field
!   in flat and full height emissions groups.
!   3) Execute the time step.
!   Processing of environmental drivers by groups allows the same arguments
!   to be used for different sets of input fields so that the same
!   'ukca_step_control' call can be used for different UKCA configurations.
!   It also reduces the number of arguments that need to be processed
!   (fields do not need to be processed individually).
! ----------------------------------------------------------------------

USE ukca_environment_req_mod, ONLY: ukca_get_envgroup_varlists
USE ukca_environment_mod, ONLY: ukca_set_environment
USE ukca_emiss_api_mod, ONLY: get_registered_ems_info, ukca_set_emission
USE ukca_step_mod, ONLY: ukca_step
USE ukca_fieldname_mod, ONLY: maxlen_fieldname
USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_env_field_mismatch

IMPLICIT NONE

! Subroutine arguments

! Model timestep number (counted from basis time at start of run)
INTEGER, INTENT(IN) :: timestep_number

! Current model time (year, month, day, hour, minute, second, day of year)
INTEGER, INTENT(IN) :: current_time(7)

! UKCA tracers. Dimensions: X,Y,Z,N
! where X is row length of tracer field (= no. of columns)
!       Y is no. of rows in tracer field
!       Z is no. of levels in tracer fields
!       N is number of tracers
REAL, ALLOCATABLE, INTENT(IN OUT) :: tracer_data(:, :)

! Non-transported prognostics. Dimensions: X,Y,Z,N
REAL, ALLOCATABLE, INTENT(IN OUT) :: ntp_data(:, :)

! Error code for status reporting
INTEGER, INTENT(OUT) :: error_code

! Model time at previous timestep (required for chemistry)
INTEGER, OPTIONAL, INTENT(IN) :: previous_time(7)

! Environmental driver field groups (ordered by dimension and type)
! Outer dimension is field index determined with reference to list of
! required variables from the group

! Scalar fields (scalar & flat grid groups)

INTEGER, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_flat_integer(:)

REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_scalar_real(:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_flat_real(:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: emissions_flat(:)

LOGICAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_flat_logical(:)

! 1D fields (flat grid plant function type tile group, vertical grid groups &
! land point group)

REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_flatpft_real(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_fullht_real(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_fullht0_real(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_fullhtp1_real(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_bllev_real(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_entlev_real(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_land_real(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: emissions_fullht(:,:)

! 2D fields (land-point tile groups)

REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_landtile_real(:,:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_landpft_real(:,:,:)

LOGICAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: envgroup_landtile_logical(:,:,:)

! Further arguments for status reporting
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: i
INTEGER :: n
INTEGER :: n_flat_ems
INTEGER :: n_fullht_ems
INTEGER :: emiss_id

REAL, ALLOCATABLE :: tmp_1d_real(:)
REAL, ALLOCATABLE :: tmp_2d_real(:,:)
REAL, ALLOCATABLE :: tmp_3d_real(:,:,:)

LOGICAL :: l_ndim_order
LOGICAL, ALLOCATABLE :: tmp_2d_logical(:,:)

CHARACTER(LEN=maxlen_fieldname), POINTER :: varnames(:)

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_STEP_CONTROL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Set environmental drivers in scalar group

IF (PRESENT(envgroup_scalar_real)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_scalar_real_ptr=varnames,           &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_scalar_real)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of environment fields (', n,                           &
          ') for scalar real group. Expecting ', SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE
      DO i = 1, n
        CALL ukca_set_environment(varnames(i), envgroup_scalar_real(i),        &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
    END IF
  END IF
END IF

! Set environmental drivers in flat grid groups
! (Input values are scalar for the column model)

IF (error_code == 0 .AND. PRESENT(envgroup_flat_integer)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_flat_integer_ptr=varnames,          &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_flat_integer)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of scalar environment fields (', n,                    &
          ') for integer group on flat grid. Expecting ', SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE
      DO i = 1, n
        CALL ukca_set_environment(varnames(i), envgroup_flat_integer(i),       &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
    END IF
  END IF
END IF

IF (error_code == 0 .AND. PRESENT(envgroup_flat_real)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_flat_real_ptr=varnames,             &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_flat_real)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of scalar environment fields (', n,                    &
          ') for real group on flat grid. Expecting ', SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE
      DO i = 1, n
        CALL ukca_set_environment(varnames(i), envgroup_flat_real(i),          &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
    END IF
  END IF
END IF

IF (error_code == 0 .AND. PRESENT(envgroup_flat_logical)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_flat_logical_ptr=varnames,          &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_flat_logical)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of scalar environment fields (', n,                    &
          ') for logical group on flat grid. Expecting ', SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE
      DO i = 1, n
        CALL ukca_set_environment(varnames(i), envgroup_flat_logical(i),       &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
    END IF
  END IF
END IF

! Set environmental drivers in flat grid plant functional type tile group
! (Input values are 1D for the column model)

IF (error_code == 0 .AND. PRESENT(envgroup_flatpft_real)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_flatpft_real_ptr=varnames,          &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_flatpft_real,DIM=2)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of environment fields (', n,                           &
          ') for real group on flat grid PFT tiles. Expecting ', SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE IF (n > 0) THEN
      ALLOCATE(tmp_1d_real(SIZE(envgroup_flatpft_real,DIM=1)))
      DO i = 1, n
        tmp_1d_real = envgroup_flatpft_real(:,i)
        CALL ukca_set_environment(varnames(i), tmp_1d_real,                    &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
      DEALLOCATE(tmp_1d_real)
    END IF
  END IF
END IF

! Set environmental drivers in full-height grid group
! (Input values are 1D for the column model)

IF (error_code == 0 .AND. PRESENT(envgroup_fullht_real)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_fullht_real_ptr=varnames,           &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_fullht_real,DIM=2)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of environment fields (', n,                           &
          ') for real group on full height grid. Expecting ', SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE IF (n > 0) THEN
      ALLOCATE(tmp_1d_real(SIZE(envgroup_fullht_real,DIM=1)))
      DO i = 1, n
        tmp_1d_real = envgroup_fullht_real(:,i)
        CALL ukca_set_environment(varnames(i), tmp_1d_real,                    &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
      DEALLOCATE(tmp_1d_real)
    END IF
  END IF
END IF

! Set environmental drivers in full-height plus zeroth level grid group
! (Input values are 1D for the column model)

IF (error_code == 0 .AND. PRESENT(envgroup_fullht0_real)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_fullht0_real_ptr=varnames,          &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_fullht0_real,DIM=2)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of environment fields (', n,                           &
          ') for real group on full height plus zero level grid. Expecting ',  &
          SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE IF (n > 0) THEN
      ALLOCATE(tmp_1d_real(0:SIZE(envgroup_fullht0_real,DIM=1)-1))
      DO i = 1, n
        tmp_1d_real = envgroup_fullht0_real(:,i)
        CALL ukca_set_environment(varnames(i), tmp_1d_real,                    &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
      DEALLOCATE(tmp_1d_real)
    END IF
  END IF
END IF

! Set environmental drivers in full-height plus one grid group
! (Input values are 1D for the column model)

IF (error_code == 0 .AND. PRESENT(envgroup_fullhtp1_real)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_fullhtp1_real_ptr=varnames,         &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_fullhtp1_real,DIM=2)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of environment fields (', n,                           &
          ') for real group on full height plus one grid. Expecting ',         &
          SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE IF (n > 0) THEN
      ALLOCATE(tmp_1d_real(SIZE(envgroup_fullhtp1_real,DIM=1)))
      DO i = 1, n
        tmp_1d_real = envgroup_fullhtp1_real(:,i)
        CALL ukca_set_environment(varnames(i), tmp_1d_real,                    &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
      DEALLOCATE(tmp_1d_real)
    END IF
  END IF
END IF

! Set environmental drivers in boundary layer levels group
! (Input values are 1D for the column model)

IF (error_code == 0 .AND. PRESENT(envgroup_bllev_real)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_bllev_real_ptr=varnames,            &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_bllev_real,DIM=2)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of environment fields (', n,                           &
          ') for real group on boundary layer levels. Expecting ',             &
          SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE IF (n > 0) THEN
      ALLOCATE(tmp_1d_real(SIZE(envgroup_bllev_real,DIM=1)))
      DO i = 1, n
        tmp_1d_real = envgroup_bllev_real(:,i)
        CALL ukca_set_environment(varnames(i), tmp_1d_real,                    &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
      DEALLOCATE(tmp_1d_real)
    END IF
  END IF
END IF

! Set environmental drivers in entrainment levels group
! (Input values are 1D for the column model)

IF (error_code == 0 .AND. PRESENT(envgroup_entlev_real)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_entlev_real_ptr=varnames,           &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_entlev_real,DIM=2)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of environment fields (', n,                           &
          ') for real group on entrainment levels. Expecting ', SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE IF (n > 0) THEN
      ALLOCATE(tmp_1d_real(SIZE(envgroup_entlev_real,DIM=1)))
      DO i = 1, n
        tmp_1d_real = envgroup_entlev_real(:,i)
        CALL ukca_set_environment(varnames(i), tmp_1d_real,                    &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
      DEALLOCATE(tmp_1d_real)
    END IF
  END IF
END IF

! Set environmental drivers in land point group
! (Input values are 1D for the column model, length 0 for sea or 1 for land)

IF (error_code == 0 .AND. PRESENT(envgroup_land_real)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_land_real_ptr=varnames,             &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_land_real,DIM=2)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of environment fields (', n,                           &
          ') for real group on land points. Expecting ', SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE IF (n > 0) THEN
      ALLOCATE(tmp_1d_real(SIZE(envgroup_land_real,DIM=1)))
      DO i = 1, n
        tmp_1d_real = envgroup_land_real(:,i)
        CALL ukca_set_environment(varnames(i), tmp_1d_real,                    &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
      DEALLOCATE(tmp_1d_real)
    END IF
  END IF
END IF

! Set environmental drivers in land-point tile groups
! (Input values are 2D for the column model, dim1 as for land point group)

IF (error_code == 0 .AND. PRESENT(envgroup_landtile_real)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_landtile_real_ptr=varnames,         &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_landtile_real,DIM=3)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of environment fields (', n,                           &
          ') for real group on land-point tiles. Expecting ', SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE IF (n > 0) THEN
      ALLOCATE(tmp_2d_real(SIZE(envgroup_landtile_real,DIM=1),                 &
                           SIZE(envgroup_landtile_real,DIM=2)))
      DO i = 1, n
        tmp_2d_real = envgroup_landtile_real(:,:,i)
        CALL ukca_set_environment(varnames(i), tmp_2d_real,                    &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
      DEALLOCATE(tmp_2d_real)
    END IF
  END IF
END IF

IF (error_code == 0 .AND. PRESENT(envgroup_landtile_logical)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_landtile_logical_ptr=varnames,      &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_landtile_logical,DIM=3)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of environment fields (', n,                           &
          ') for logical group on land-point tiles. Expecting ', SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE IF (n > 0) THEN
      ALLOCATE(tmp_2d_logical(SIZE(envgroup_landtile_logical,DIM=1),           &
                              SIZE(envgroup_landtile_logical,DIM=2)))
      DO i = 1, n
        tmp_2d_logical = envgroup_landtile_logical(:,:,i)
        CALL ukca_set_environment(varnames(i), tmp_2d_logical,                 &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
      DEALLOCATE(tmp_2d_logical)
    END IF
  END IF
END IF

! Set environmental drivers in land-point plant functional type tile group
! (Input values are 2D for the column model, dim1 as for land point group)

IF (error_code == 0 .AND. PRESENT(envgroup_landpft_real)) THEN
  CALL ukca_get_envgroup_varlists(error_code,                                  &
                                  varnames_landpft_real_ptr=varnames,          &
                                  error_message=error_message,                 &
                                  error_routine=error_routine)
  IF (error_code == 0) THEN
    n = SIZE(envgroup_landpft_real,DIM=3)
    IF (n /= SIZE(varnames)) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,I0,A,I0,A)')                                  &
          'Wrong number of environment fields (', n,                           &
          ') for real group on land-point PFT tiles. Expecting ', SIZE(varnames)
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE IF (n > 0) THEN
      ALLOCATE(tmp_2d_real(SIZE(envgroup_landpft_real,DIM=1),                  &
                           SIZE(envgroup_landpft_real,DIM=2)))
      DO i = 1, n
        tmp_2d_real = envgroup_landpft_real(:,:,i)
        CALL ukca_set_environment(varnames(i), tmp_2d_real,                    &
                                  error_code, error_message=error_message,     &
                                  error_routine=error_routine)
      END DO
      DEALLOCATE(tmp_2d_real)
    END IF
  END IF
END IF

! Get information about registered emissions

CALL get_registered_ems_info(n_flat_ems, n_fullht_ems, l_ndim_order)

IF (.NOT. l_ndim_order) THEN
  error_code = errcode_env_field_mismatch
  IF (PRESENT(error_message))                                                  &
    error_message =                                                            &
      'Emissions have not been registered in order of their dimensionality'
  IF (PRESENT(error_routine)) error_routine = RoutineName
END IF

! Initialise sequential emissions id
emiss_id = 0

! Set flat emissions fields
! (Input values are scalar for the column model but must be converted to 3D
! since 'ukca_set_emission' does not currently allow a scalar argument)

IF (error_code == 0 .AND. n_flat_ems > 0) THEN
  IF (PRESENT(emissions_flat)) THEN
    n = SIZE(emissions_flat)
  ELSE
    n = 0
  END IF
  IF (n /= n_flat_ems) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message))                                                &
      WRITE(error_message, '(A,I0,A,I0,A)')                                    &
        'Wrong no. of flat emissions provided (', n, '). Expecting ',          &
        n_flat_ems
    IF (PRESENT(error_routine)) error_routine = RoutineName
  ELSE
    ALLOCATE(tmp_3d_real(1,1,1))
    DO i = 1, n
      emiss_id = emiss_id + 1
      tmp_3d_real(1,1,1) = emissions_flat(i)
      CALL ukca_set_emission(emiss_id, tmp_3d_real)
    END DO
    DEALLOCATE(tmp_3d_real)
  END IF
END IF

! Set full height emissions fields
! (Input values are 1D for the column model nut must be converted to 3D since
! 'ukca_set_emission' does not currently allow a 1D argument)

IF (error_code == 0 .AND. n_fullht_ems > 0) THEN
  IF (PRESENT(emissions_fullht)) THEN
    n = SIZE(emissions_fullht,DIM=2)
  ELSE
    n = 0
  END IF
  IF (n /= n_fullht_ems) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message))                                                &
      WRITE(error_message, '(A,I0,A,I0,A)')                                    &
        'Wrong no. of full height emissions provided (', n, '). Expecting ',   &
        n_fullht_ems
    IF (PRESENT(error_routine)) error_routine = RoutineName
  ELSE
    ALLOCATE(tmp_3d_real(1,1,SIZE(emissions_fullht,DIM=1)))
    DO i = 1, n
      emiss_id = emiss_id + 1
      tmp_3d_real(1,1,:) = emissions_fullht(:,i)
      CALL ukca_set_emission(emiss_id, tmp_3d_real)
    END DO
    DEALLOCATE(tmp_3d_real)
  END IF
END IF

! If any errors encountered above then bail out
IF (error_code /= 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Do the time step
CALL ukca_step(timestep_number, current_time,                                  &
               tracer_data, ntp_data, error_code,                              &
               previous_time=previous_time,                                    &
               error_message=error_message, error_routine=error_routine)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_step_control

END MODULE ukca_step_control_mod

