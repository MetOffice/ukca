! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module for controlling each Photolysis call in the parent application
!   where all data transfer between the parent and Photolyis during a run is
!   restricted to a single API call for each time step.
!   This is a necessary condition for running multiple Photolysis step calls in
!   parallel on the same node with shared memory (one time step call for
!   each parent model domain or sub-domain).
!
!   The module provides the following procedure for the Photolysis API.
!
!     photol_step_control - Obtain all required environmental driver data
!     and perform one Photolysis step for given sub-domain.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds,
! University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_Photolysis
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE photol_step_control_mod

USE photol_fieldname_mod,   ONLY:  fieldname_len,                              &
                                   photol_varname_len,                         &
                                   fldname_aod_sulph_aitk,                     &
                                   fldname_aod_sulph_accum,                    &
                                   fldname_area_cloud_fraction,                &
                                   fldname_conv_cloud_amount,                  &
                                   fldname_conv_cloud_base,                    &
                                   fldname_conv_cloud_lwp,                     &
                                   fldname_conv_cloud_top,                     &
                                   fldname_cos_latitude,                       &
                                   fldname_equation_of_time,                   &
                                   fldname_land_fraction,                      &
                                   fldname_longitude,                          &
                                   fldname_ozone_mmr,                          &
                                   fldname_p_layer_boundaries,                 &
                                   fldname_p_theta_levels,                     &
                                   fldname_photol_rates_2d,                    &
                                   fldname_qcf,                                &
                                   fldname_qcl,                                &
                                   fldname_rad_ctl_jo2,                        &
                                   fldname_rad_ctl_jo2b,                       &
                                   fldname_r_rho_levels,                       &
                                   fldname_r_theta_levels,                     &
                                   fldname_sec_since_midnight,                 &
                                   fldname_sin_declination,                    &
                                   fldname_sin_latitude,                       &
                                   fldname_so4_aitken,                         &
                                   fldname_so4_accum,                          &
                                   fldname_surf_albedo,                        &
                                   fldname_t_theta_levels,                     &
                                   fldname_tan_latitude,                       &
                                   fldname_z_top_of_model

USE photol_environment_mod, ONLY: l_photol_environ_req_available,              &
                                  environ_fldnames_scalar_real,                &
                                  environ_fldnames_flat_integer,               &
                                  environ_fldnames_flat_real,                  &
                                  environ_fldnames_fullht_real,                &
                                  environ_fldnames_fullht0_real,               &
                                  environ_fldnames_fullhtphot_real

USE photol_check_environment_mod, ONLY: check_env_fields_list,                 &
                                        check_environ_group

USE photol_ctl_mod, ONLY: photol_ctl
USE photol_config_specification_mod, ONLY: photol_config

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname, error_report,       &
                          errcode_env_field_mismatch, errcode_env_req_uninit,  &
                          errcode_env_field_missing, errcode_value_invalid
! Dr Hook modules
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC photol_step_control

! Driving fields for Photolysis call -arranged by type
! Scalar -real
REAL  :: seconds_since_midnight
REAL  :: sin_declination
REAL  :: equation_of_time
REAL  :: z_top_of_model
! Flat Int and Real
INTEGER, ALLOCATABLE   :: conv_cloud_base(:,:)
INTEGER, ALLOCATABLE   :: conv_cloud_top(:,:)
REAL, ALLOCATABLE   :: conv_cloud_lwp(:,:)
REAL, ALLOCATABLE   :: surf_albedo(:,:)
REAL, ALLOCATABLE   :: land_fraction(:,:)
REAL, ALLOCATABLE   :: longitude(:,:)
REAL, ALLOCATABLE   :: sin_latitude(:,:)
REAL, ALLOCATABLE   :: cos_latitude(:,:)
REAL, ALLOCATABLE   :: tan_latitude(:,:)
! Fullht Real
REAL, ALLOCATABLE   :: p_theta_levels(:,:,:)
REAL, ALLOCATABLE   :: r_rho_levels(:,:,:)
REAL, ALLOCATABLE   :: qcl(:,:,:)
REAL, ALLOCATABLE   :: qcf(:,:,:)
REAL, ALLOCATABLE   :: area_cloud_fraction(:,:,:)
REAL, ALLOCATABLE   :: conv_cloud_amount(:,:,:)
REAL, ALLOCATABLE   :: ozone_mmr(:,:,:)
REAL, ALLOCATABLE   :: so4_aitken(:,:,:)
REAL, ALLOCATABLE   :: so4_accum(:,:,:)
REAL, ALLOCATABLE   :: aod_sulph_aitk(:,:,:)
REAL, ALLOCATABLE   :: aod_sulph_accum(:,:,:)
REAL, ALLOCATABLE   :: rad_ctl_jo2(:,:,:)
REAL, ALLOCATABLE   :: rad_ctl_jo2b(:,:,:)
REAL, ALLOCATABLE   :: t_theta_levels(:,:,:)

! Fullht0 Real
REAL, ALLOCATABLE   :: p_layer_boundaries(:,:,:)
REAL, ALLOCATABLE   :: r_theta_levels(:,:,:)

! Fullhtphot Real
REAL, ALLOCATABLE   :: photol_rates_2d(:,:,:,:)

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'PHOTOL_STEP_CONTROL_MOD'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE photol_step_control(current_time, x_dim, y_dim, z_dim, jppj,        &
                               ratj_data, ratj_varnames, error_code,           &
                               photol_rates,                                   &
                               ! Names of environ fields supplied, optional
                               envfield_names_in,                              &
                               ! Environment Field groups, optional
                               envgroup_scalar_real,                           &
                               envgroup_flat_integer,                          &
                               envgroup_flat_real,                             &
                               envgroup_fullht_real,                           &
                               envgroup_fullht0_real,                          &
                               envgroup_fullhtphot_real,                       &
                               ! Error return info
                               error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Unpacks the grouped environmental driver fields and performs
!   one Photolysis time step for the given model sub-domain.
!
! Method:
!   Environmental driver fields are grouped according to their
!   dimensionality, size and type.
!   The procedure for handling the Photolysis call is:
!   1) Verify that all required fields are provided by comparing the full
!      environment fieldnames, as well as individual field group items
!   2) Unpack all the field groups into individual driving fields
!   3) For environment fields not required from parent for this configuration
!      allocate minimal sizes with default values
!   4) Call the top-level Photol_ctl routine by passing the dimensions as well
!      as all the driving fields as arguments
!
!  Passing the (horizontal) spatial dimensions as arguments for each
!  call and not pre-setting anywhere in the photolysis routines will be an
!  useful step towards eveuentually calling multiple instances of photolysis
!  in a thread-safe manner.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments

! Current model time (year, month, day, hour, minute, second, day of year)
INTEGER, INTENT(IN) :: current_time(7)

INTEGER, INTENT(IN) :: x_dim, y_dim, z_dim  ! Dimensions of the data supplied
INTEGER, INTENT(IN) :: jppj                 ! Number of photolytic species

! Photolysis species names
CHARACTER(LEN=photol_varname_len), POINTER, INTENT(IN) :: ratj_data(:,:)
CHARACTER(LEN=photol_varname_len), POINTER, INTENT(IN) :: ratj_varnames(:)

! Names of environment fields provided by the parent driving routine
CHARACTER(LEN=fieldname_len), OPTIONAL, INTENT(IN) :: envfield_names_in(:)

! Environmental driver field groups (ordered by dimension and type)
! Outer dimension is field index determined with reference to list of
! required variables from the group
INTEGER, OPTIONAL, INTENT(IN) :: envgroup_flat_integer(:,:,:)
REAL, OPTIONAL, INTENT(IN) :: envgroup_scalar_real(:)
REAL, OPTIONAL, INTENT(IN) :: envgroup_flat_real(:,:,:)
REAL, OPTIONAL, INTENT(IN) :: envgroup_fullht_real(:,:,:,:)
REAL, OPTIONAL, INTENT(IN) :: envgroup_fullht0_real(:,:,:,:)
REAL, OPTIONAL, INTENT(IN) :: envgroup_fullhtphot_real(:,:,:,:,:)

! Error code for status reporting
INTEGER, TARGET, INTENT(OUT) :: error_code

! Calculated Photolysis rates, OUT
REAL, ALLOCATABLE, INTENT(OUT) :: photol_rates(:,:,:,:)

! Further arguments for status reporting
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: n
! Flag that denotes if any fields of particular type have been provided
LOGICAL :: l_fields_this_group

! Local copy for receiving error output from called routines
INTEGER, POINTER :: error_code_ptr
CHARACTER(LEN=maxlen_message) :: err_msg

! Dr Hook data
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PHOTOL_STEP_CONTROL'
! ------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr => error_code
error_code_ptr = 0
err_msg = ''
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check that environment fields requirements have been set up
IF ( .NOT. l_photol_environ_req_available ) THEN
  error_code_ptr = errcode_env_req_uninit
  err_msg = 'Photolysis environment request not set'
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
         RoutineName, msg_out=error_message, locn_out=error_routine)

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (.NOT. ALLOCATED(photol_rates))                                             &
  ALLOCATE(photol_rates(x_dim,y_dim,z_dim,jppj))
photol_rates(:,:,:,:) = 0.0

! If envfield_names_in is supplied by parent, check all required fields are
! being provided- irrespective of group
IF ( PRESENT(envfield_names_in) ) THEN
  CALL check_env_fields_list(envfield_names_in, error_code_ptr,                &
         error_message=error_message, error_routine=error_routine)
  IF ( error_code_ptr > 0 ) THEN
    IF (ALLOCATED(photol_rates)) DEALLOCATE(photol_rates)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                            zhook_handle)
    RETURN
  END IF
END IF

! For each field group, check number and dimensions via chk_environ_group
! Proceed to next group only if no errors in previous step (errcode == 0)
! ** Assumption while unpacking = Fields are grouped in same order as names **

! Unpack environmental drivers in scalar real group
IF (PRESENT(envgroup_scalar_real)) THEN
  CALL check_environ_group(envgroup_scalar_real, 'scalar_real',                &
          error_code_ptr, l_fields_this_group, error_message=error_message,    &
          error_routine=error_routine)
  IF ( l_fields_this_group .AND. error_code <= 0 ) THEN
    DO n = 1, SIZE(environ_fldnames_scalar_real)
      SELECT CASE(environ_fldnames_scalar_real(n))
      CASE (fldname_equation_of_time)
        equation_of_time = envgroup_scalar_real(n)
      CASE (fldname_sec_since_midnight)
        seconds_since_midnight = envgroup_scalar_real(n)
      CASE (fldname_sin_declination)
        sin_declination = envgroup_scalar_real(n)
      CASE (fldname_z_top_of_model)
        z_top_of_model = envgroup_scalar_real(n)
      CASE DEFAULT
        error_code_ptr = errcode_env_field_mismatch
        WRITE(err_msg,'(A,A)') 'Unknown SCALAR_REAL env field request ',       &
          environ_fldnames_scalar_real(n)
      END SELECT
    END DO
  END IF
END IF

! Unpack environmental drivers in flat integer group
IF ( error_code_ptr <= 0 .AND. PRESENT(envgroup_flat_integer) ) THEN
  CALL check_environ_group(envgroup_flat_integer, 'flat_integer',              &
             x_dim, y_dim, error_code_ptr, l_fields_this_group,                &
             error_message=error_message, error_routine=error_routine)
  IF ( l_fields_this_group .AND. error_code_ptr <= 0 ) THEN
    DO n = 1, SIZE(environ_fldnames_flat_integer)
      SELECT CASE(environ_fldnames_flat_integer(n))
      CASE (fldname_conv_cloud_base)
        ALLOCATE(conv_cloud_base(x_dim,y_dim))
        conv_cloud_base(:,:) = envgroup_flat_integer(:,:,n)
      CASE (fldname_conv_cloud_top)
        ALLOCATE(conv_cloud_top(x_dim,y_dim))
        conv_cloud_top(:,:) = envgroup_flat_integer(:,:,n)
      CASE DEFAULT
        error_code_ptr = errcode_env_field_mismatch
        WRITE(err_msg,'(A,A)')'Unknown FLAT INTEGER env field request ',       &
          environ_fldnames_flat_integer(n)
      END SELECT
    END DO
  END IF
END IF

! Unpack environmental drivers in flat real group
IF ( error_code_ptr <= 0 .AND. PRESENT(envgroup_flat_real) ) THEN
  CALL check_environ_group(envgroup_flat_real, 'flat_real',                    &
             x_dim, y_dim, error_code_ptr, l_fields_this_group,                &
             error_message=error_message, error_routine=error_routine)
  IF ( l_fields_this_group .AND. error_code_ptr <= 0 ) THEN
    DO n = 1, SIZE(environ_fldnames_flat_real)
      SELECT CASE(environ_fldnames_flat_real(n))
      CASE (fldname_conv_cloud_lwp)
        ALLOCATE(conv_cloud_lwp(x_dim,y_dim))
        conv_cloud_lwp(:,:) = envgroup_flat_real(:,:,n)
      CASE (fldname_cos_latitude)
        ALLOCATE(cos_latitude(x_dim,y_dim))
        cos_latitude(:,:) = envgroup_flat_real(:,:,n)
      CASE (fldname_land_fraction)
        ALLOCATE(land_fraction(x_dim,y_dim))
        land_fraction(:,:) = envgroup_flat_real(:,:,n)
      CASE (fldname_longitude)
        ALLOCATE(longitude(x_dim,y_dim))
        longitude(:,:) = envgroup_flat_real(:,:,n)
      CASE (fldname_sin_latitude)
        ALLOCATE(sin_latitude(x_dim,y_dim))
        sin_latitude(:,:) = envgroup_flat_real(:,:,n)
      CASE (fldname_surf_albedo)
        ALLOCATE(surf_albedo(x_dim,y_dim))
        surf_albedo(:,:) = envgroup_flat_real(:,:,n)
      CASE (fldname_tan_latitude)
        ALLOCATE(tan_latitude(x_dim,y_dim))
        tan_latitude(:,:) = envgroup_flat_real(:,:,n)
      CASE DEFAULT
        error_code_ptr = errcode_env_field_mismatch
        WRITE(err_msg,'(A,A)') 'Unknown FLAT REAL env field request ',         &
          environ_fldnames_flat_real(n)
      END SELECT
    END DO
  END IF
END IF

! Unpack environmental drivers in fullht real group
IF ( error_code_ptr <= 0 .AND. PRESENT(envgroup_fullht_real) ) THEN
  CALL check_environ_group(envgroup_fullht_real, 'fullht_real',                &
       x_dim, y_dim, z_dim, error_code_ptr, l_fields_this_group,               &
       error_message=error_message, error_routine=error_routine)
  IF ( l_fields_this_group .AND. error_code_ptr <= 0 ) THEN
    DO n = 1, SIZE(environ_fldnames_fullht_real)
      SELECT CASE(environ_fldnames_fullht_real(n))
      CASE (fldname_aod_sulph_aitk)
        ALLOCATE(aod_sulph_aitk(x_dim,y_dim,z_dim))
        aod_sulph_aitk(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_aod_sulph_accum)
        ALLOCATE(aod_sulph_accum(x_dim,y_dim,z_dim))
        aod_sulph_accum(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_area_cloud_fraction)
        ALLOCATE(area_cloud_fraction(x_dim,y_dim,z_dim))
        area_cloud_fraction(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_conv_cloud_amount)
        ALLOCATE(conv_cloud_amount(x_dim,y_dim,z_dim))
        conv_cloud_amount(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_ozone_mmr)
        ALLOCATE(ozone_mmr(x_dim,y_dim,z_dim))
        ozone_mmr(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_p_theta_levels)
        ALLOCATE(p_theta_levels(x_dim,y_dim,z_dim))
        p_theta_levels(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_qcf)
        ALLOCATE(qcf(x_dim,y_dim,z_dim))
        qcf(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_qcl)
        ALLOCATE(qcl(x_dim,y_dim,z_dim))
        qcl(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_rad_ctl_jo2)
        ALLOCATE(rad_ctl_jo2(x_dim,y_dim,z_dim))
        rad_ctl_jo2(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_rad_ctl_jo2b)
        ALLOCATE(rad_ctl_jo2b(x_dim,y_dim,z_dim))
        rad_ctl_jo2b(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_r_rho_levels)
        ALLOCATE(r_rho_levels(x_dim,y_dim,z_dim))
        r_rho_levels(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_so4_aitken)
        ALLOCATE(so4_aitken(x_dim,y_dim,z_dim))
        so4_aitken(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_so4_accum)
        ALLOCATE(so4_accum(x_dim,y_dim,z_dim))
        so4_accum(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE (fldname_t_theta_levels)
        ALLOCATE(t_theta_levels(x_dim,y_dim,z_dim))
        t_theta_levels(:,:,:) = envgroup_fullht_real(:,:,:,n)
      CASE DEFAULT
        error_code_ptr = errcode_env_field_mismatch
        WRITE(err_msg,'(A,A)')'Unknown FULLHT REAL env field request ',        &
          environ_fldnames_fullht_real(n)
      END SELECT
    END DO
  END IF
END IF

! Unpack environmental drivers in fullht0 real group
IF ( error_code_ptr <= 0 .AND. PRESENT(envgroup_fullht0_real) ) THEN
  CALL check_environ_group(envgroup_fullht0_real, 'fullht0_real',              &
      x_dim, y_dim, z_dim, error_code_ptr, l_fields_this_group,                &
      error_message=error_message, error_routine=error_routine)
  IF ( l_fields_this_group .AND. error_code_ptr <= 0 ) THEN
    DO n = 1, SIZE(environ_fldnames_fullht0_real)
      SELECT CASE(environ_fldnames_fullht0_real(n))
      CASE (fldname_p_layer_boundaries)
        ALLOCATE(p_layer_boundaries(x_dim,y_dim,z_dim+1))
        p_layer_boundaries(:,:,:) = envgroup_fullht0_real(:,:,:,n)
      CASE (fldname_r_theta_levels)
        ALLOCATE(r_theta_levels(x_dim,y_dim,z_dim+1))
        r_theta_levels(:,:,:) = envgroup_fullht0_real(:,:,:,n)
      CASE DEFAULT
        error_code_ptr = errcode_env_field_mismatch
        WRITE(err_msg,'(A,A)')'Unknown FULLHT0 REAL env field request ',       &
          environ_fldnames_fullht0_real(n)
      END SELECT
    END DO
  END IF
END IF

! Unpack environmental drivers in fullhtphot real group
IF ( error_code_ptr <= 0 .AND. PRESENT(envgroup_fullhtphot_real) ) THEN
  CALL check_environ_group(envgroup_fullhtphot_real, 'fullhtphot_real',        &
      x_dim, y_dim, z_dim, jppj, error_code_ptr, l_fields_this_group,          &
      error_message=error_message, error_routine=error_routine)
  IF ( l_fields_this_group .AND. error_code_ptr <= 0 ) THEN
    DO n = 1, SIZE(environ_fldnames_fullhtphot_real)
      SELECT CASE(environ_fldnames_fullhtphot_real(n))
      CASE (fldname_photol_rates_2d)
        ALLOCATE(photol_rates_2d(x_dim,y_dim,z_dim,jppj))
        photol_rates_2d(:,:,:,:) = envgroup_fullhtphot_real(:,:,:,:,n)
      CASE DEFAULT
        error_code_ptr = errcode_env_field_mismatch
        WRITE(err_msg,'(A,A)')'Unknown FULLHTPHOT REAL env field request ',    &
          environ_fldnames_fullhtphot_real(n)
      END SELECT
    END DO
  END IF
END IF

! Check if unpacking of the groups above has raised an error
IF (error_code_ptr > 0 .AND. err_msg /= '') THEN
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_msg,     &
         RoutineName, msg_out=error_message, locn_out=error_routine)

  CALL deallocate_drv_fields()
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (error_code_ptr <= 0) THEN
  ! Allocate any remaining arrays to minimum size to enable passing as
  ! arguments to photol_ctl.
  ! The unrequired arrays should not be accessed in photolysis routines.
  IF (.NOT. ALLOCATED(conv_cloud_base)) THEN
    ALLOCATE(conv_cloud_base(1,1))
    conv_cloud_base(:,:) = 0
  END IF
  IF (.NOT. ALLOCATED(conv_cloud_top)) THEN
    ALLOCATE(conv_cloud_top(1,1))
    conv_cloud_top(:,:) = 0
  END IF
  IF (.NOT. ALLOCATED(cos_latitude)) THEN
    ALLOCATE(cos_latitude(1,1))
    cos_latitude(:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(land_fraction)) THEN
    ALLOCATE(land_fraction(1,1))
    land_fraction(:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(longitude)) THEN
    ALLOCATE(longitude(1,1))
    longitude(:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(sin_latitude)) THEN
    ALLOCATE(sin_latitude(1,1))
    sin_latitude(:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(surf_albedo)) THEN
    ALLOCATE(surf_albedo(1,1))
    surf_albedo(:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(tan_latitude)) THEN
    ALLOCATE(tan_latitude(1,1))
    tan_latitude(:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(aod_sulph_aitk)) THEN
    ALLOCATE(aod_sulph_aitk(1,1,1))
    aod_sulph_aitk(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(aod_sulph_accum)) THEN
    ALLOCATE(aod_sulph_accum(1,1,1))
    aod_sulph_accum(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(area_cloud_fraction)) THEN
    ALLOCATE(area_cloud_fraction(1,1,1))
    area_cloud_fraction(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(conv_cloud_amount)) THEN
    ALLOCATE(conv_cloud_amount(1,1,1))
    conv_cloud_amount(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(conv_cloud_lwp)) THEN
    ALLOCATE(conv_cloud_lwp(1,1))
    conv_cloud_lwp(:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(ozone_mmr)) THEN
    ALLOCATE(ozone_mmr(1,1,1))
    ozone_mmr(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(p_theta_levels)) THEN
    ALLOCATE(p_theta_levels(1,1,1))
    p_theta_levels(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(qcf)) THEN
    ALLOCATE(qcf(1,1,1))
    qcf(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(qcl)) THEN
    ALLOCATE(qcl(1,1,1))
    qcl(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(rad_ctl_jo2)) THEN
    ALLOCATE(rad_ctl_jo2(1,1,1))
    rad_ctl_jo2(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(rad_ctl_jo2b)) THEN
    ALLOCATE(rad_ctl_jo2b(1,1,1))
    rad_ctl_jo2b(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(r_rho_levels)) THEN
    ALLOCATE(r_rho_levels(1,1,1))
    r_rho_levels(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(so4_aitken)) THEN
    ALLOCATE(so4_aitken(1,1,1))
    so4_aitken(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(so4_accum)) THEN
    ALLOCATE(so4_accum(1,1,1))
    so4_accum(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(t_theta_levels)) THEN
    ALLOCATE(t_theta_levels(1,1,1))
    t_theta_levels(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(p_layer_boundaries)) THEN
    ALLOCATE(p_layer_boundaries(1,1,1))
    p_layer_boundaries(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(r_theta_levels)) THEN
    ALLOCATE(r_theta_levels(1,1,1))
    r_theta_levels(:,:,:) = 0.0
  END IF
  IF (.NOT. ALLOCATED(photol_rates_2d)) THEN
    ALLOCATE(photol_rates_2d(1,1,1,1))
    photol_rates_2d(:,:,:,:) = 0.0
  END IF
  ! Calculate Photolysis rates
  CALL photol_ctl(error_code_ptr, x_dim, y_dim, z_dim, jppj, z_top_of_model,   &
                  seconds_since_midnight, sin_declination, equation_of_time,   &
                  current_time, ratj_data, ratj_varnames,                      &
                  conv_cloud_base, conv_cloud_top, land_fraction, surf_albedo, &
                  longitude, sin_latitude, cos_latitude, tan_latitude,         &
                  p_theta_levels, p_layer_boundaries, r_theta_levels,          &
                  r_rho_levels, qcl, qcf, area_cloud_fraction,                 &
                  conv_cloud_amount, conv_cloud_lwp, ozone_mmr, so4_aitken,    &
                  so4_accum, aod_sulph_aitk, aod_sulph_accum, rad_ctl_jo2,     &
                  rad_ctl_jo2b, t_theta_levels, photol_rates_2d, photol_rates, &
                  error_message=error_message, error_routine=error_routine)

END IF   ! Error code <= 0

CALL deallocate_drv_fields()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE photol_step_control

! ----------------------------------------------------------------------------
!
! Routine to Deallocate driving field arrays
! ----------------------------------------------------------------------------

SUBROUTINE deallocate_drv_fields()

IMPLICIT NONE

IF (ALLOCATED(photol_rates_2d)) DEALLOCATE(photol_rates_2d)
IF (ALLOCATED(r_theta_levels)) DEALLOCATE(r_theta_levels)
IF (ALLOCATED(p_layer_boundaries)) DEALLOCATE(p_layer_boundaries)

IF (ALLOCATED(t_theta_levels)) DEALLOCATE(t_theta_levels)
IF (ALLOCATED(so4_accum)) DEALLOCATE(so4_accum)
IF (ALLOCATED(so4_aitken)) DEALLOCATE(so4_aitken)
IF (ALLOCATED(r_rho_levels)) DEALLOCATE(r_rho_levels)
IF (ALLOCATED(rad_ctl_jo2b)) DEALLOCATE(rad_ctl_jo2b)
IF (ALLOCATED(rad_ctl_jo2)) DEALLOCATE(rad_ctl_jo2)
IF (ALLOCATED(qcl)) DEALLOCATE(qcl)
IF (ALLOCATED(qcf)) DEALLOCATE(qcf)
IF (ALLOCATED(p_theta_levels)) DEALLOCATE(p_theta_levels)
IF (ALLOCATED(ozone_mmr)) DEALLOCATE(ozone_mmr)
IF (ALLOCATED(conv_cloud_amount)) DEALLOCATE(conv_cloud_amount)
IF (ALLOCATED(area_cloud_fraction)) DEALLOCATE(area_cloud_fraction)
IF (ALLOCATED(aod_sulph_accum)) DEALLOCATE(aod_sulph_accum)
IF (ALLOCATED(aod_sulph_aitk)) DEALLOCATE(aod_sulph_aitk)

IF (ALLOCATED(tan_latitude)) DEALLOCATE(tan_latitude)
IF (ALLOCATED(surf_albedo)) DEALLOCATE(surf_albedo)
IF (ALLOCATED(sin_latitude)) DEALLOCATE(sin_latitude)
IF (ALLOCATED(longitude)) DEALLOCATE(longitude)
IF (ALLOCATED(land_fraction)) DEALLOCATE(land_fraction)
IF (ALLOCATED(cos_latitude)) DEALLOCATE(cos_latitude)
IF (ALLOCATED(conv_cloud_lwp)) DEALLOCATE(conv_cloud_lwp)

IF (ALLOCATED(conv_cloud_top)) DEALLOCATE(conv_cloud_top)
IF (ALLOCATED(conv_cloud_base)) DEALLOCATE(conv_cloud_base)

RETURN
END SUBROUTINE deallocate_drv_fields

END MODULE photol_step_control_mod
