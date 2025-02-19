! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module handling the reading in and processing of photolysis data,
!  before this data is passed to UKCA.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_Photolysis
!
!  Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE photol_ctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PHOTOL_CTL_MOD'

CONTAINS

SUBROUTINE photol_ctl(error_code_ptr, row_length, rows, model_levels, jppj,    &
                      z_top_of_model, seconds_since_midnight, sin_declination, &
                      equation_of_time, current_time,                          &
                      ratj_data, ratj_varnames, conv_cloud_base,               &
                      conv_cloud_top, land_fraction, surf_albedo,              &
                      longitude, sin_latitude, cos_latitude, tan_latitude,     &
                      p_theta_levels, p_layer_boundaries, r_theta_levels,      &
                      r_rho_levels, qcl, qcf, area_cloud_fraction,             &
                      conv_cloud_amount, conv_cloud_lwp, ozone_mmr,            &
                      so4_aitken, so4_accum, aod_sulph_aitk, aod_sulph_accum,  &
                      rad_ctl_jo2, rad_ctl_jo2b, t_theta_levels,               &
                      photol_rates_2d, photol_rates,                           &
                      error_message, error_routine)
!
! Purpose: Subroutine to handle reading in and processing of photolysis
! data before this data is passed to UKCA.
!
!          Called from PHOTOL_STEP_CONTROL
!
! ---------------------------------------------------------------------
!

USE ukca_error_mod,       ONLY: maxlen_message, maxlen_procname,               &
                                errcode_value_invalid, error_report

USE ukca_fastjx_mod,    ONLY: ukca_fastjx
USE photol_calc_ozonecol_mod,  ONLY: photol_calc_ozonecol
USE photol_fieldname_mod, ONLY: photol_varname_len
USE photol_config_specification_mod,                                           &
                          ONLY: photol_config, fjx_mode_fastjx,                &
                              fjx_mode_merged, fjx_mode_2Donly,                &
                              photolysis_off => i_scheme_nophot,               &
                              photolysis_2d => i_scheme_phot2d,                &
                              photolysis_fastjx => i_scheme_fastjx
USE photol_constants_mod, ONLY: rhour_per_day => const_rhour_per_day,          &
                              pi => const_pi, c_o3 => const_o3_mmr_vmr,        &
                              planet_radius => const_planet_radius
USE photol_solang_mod,  ONLY: photol_solang

USE ukca_um_strat_photol_mod,  ONLY: strat_photol
USE ukca_um_dissoc_mod,        ONLY: strat_photol_init, strat_photol_dealloc

! UM profiling and error handling
USE yomhook,                  ONLY: lhook, dr_hook
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

! -- Input arguments -- !
INTEGER, POINTER, INTENT(IN) :: error_code_ptr
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: jppj
CHARACTER(LEN=photol_varname_len), POINTER, INTENT(IN) :: ratj_data(:,:)
CHARACTER(LEN=photol_varname_len), POINTER, INTENT(IN) :: ratj_varnames(:)
REAL, INTENT(IN)    :: land_fraction(row_length,rows)
REAL, INTENT(IN)    :: longitude(row_length, rows)
REAL, INTENT(IN)    :: sin_latitude(row_length, rows)
REAL, INTENT(IN)    :: cos_latitude(row_length, rows)
REAL, INTENT(IN)    :: tan_latitude(row_length, rows)
REAL, INTENT(IN)    :: p_theta_levels(row_length, rows, model_levels)
REAL, INTENT(IN)    :: p_layer_boundaries(row_length, rows, 0:model_levels)
REAL, INTENT(IN)    :: r_theta_levels(row_length,rows,0:model_levels)
REAL, INTENT(IN)    :: r_rho_levels(row_length,rows, model_levels)
REAL, INTENT(IN)    :: qcl(row_length, rows, model_levels)
REAL, INTENT(IN)    :: qcf(row_length, rows, model_levels)
REAL, INTENT(IN)    :: area_cloud_fraction(row_length, rows, model_levels)
REAL, INTENT(IN)    :: conv_cloud_amount(row_length, rows, model_levels)
REAL, INTENT(IN)    :: conv_cloud_lwp(row_length, rows)
INTEGER, INTENT(IN) :: conv_cloud_base(row_length, rows)
INTEGER, INTENT(IN) :: conv_cloud_top(row_length, rows)
REAL, INTENT(IN)    :: surf_albedo(row_length, rows)

REAL, INTENT(IN)    :: z_top_of_model
REAL, INTENT(IN)    :: seconds_since_midnight
REAL, INTENT(IN)    :: sin_declination
REAL, INTENT(IN)    :: equation_of_time
INTEGER, INTENT(IN) :: current_time(7)
REAL, INTENT(IN)    :: ozone_mmr(row_length, rows, model_levels)
! Aerosol MMRs and Optical Depths required for Fast-JX
REAL, INTENT(IN)    :: so4_aitken(row_length, rows, model_levels)
REAL, INTENT(IN)    :: so4_accum(row_length, rows, model_levels)
REAL, INTENT(IN)    :: aod_sulph_aitk(row_length, rows, model_levels)
REAL, INTENT(IN)    :: aod_sulph_accum(row_length, rows, model_levels)
! Photol rates drectly from Radiation scheme
REAL, INTENT(IN)    :: rad_ctl_jo2(row_length, rows, model_levels)
REAL, INTENT(IN)    :: rad_ctl_jo2b(row_length, rows, model_levels)
! t_theta_levels also required for strat_photol and Fast-JX
REAL, INTENT(IN)    :: t_theta_levels(row_length, rows, model_levels)
! 2-D (tabulated/ Offline) Photolysis rates
REAL, INTENT(IN)    :: photol_rates_2d(row_length,rows,model_levels,jppj)
! Out: Photolysis rates
REAL, INTENT(OUT)   :: photol_rates(row_length, rows, model_levels, jppj)

! error handling arguments
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
                                                       ! Error return message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
                                         ! Routine in which error was trapped

! -- Local variables -- !
REAL, ALLOCATABLE, SAVE :: photol_rates_fastjx(:,:,:,:)
REAL                    :: photol2d(row_length, rows, jppj)
REAL                    :: photol_strat(row_length, rows, jppj)
! Temporary SO4 aerosol arrays
REAL                    :: photol_so4_aitken(row_length, rows, model_levels)
REAL                    :: photol_so4_accum(row_length, rows, model_levels)
! cos of solar zenith angle, used for strat_photol
REAL                    :: cos_zenith_angle(row_length, rows)

! ozonecol for stratospheric chemistry
REAL                    :: ozonecol(row_length, rows, model_levels)

! local copy of model top height - if needed to calculate here
REAL                    :: loc_z_top_model

! loop variables
INTEGER                 :: i
INTEGER                 :: j
INTEGER                 :: k
INTEGER                 :: l
! size of theta levels
INTEGER                 :: theta_field_size

! Local time variables
REAL                    :: r_secs_per_step

! if ukca_photin called, set to true (to prevent 2D photolysis data being
! read twice)
LOGICAL, SAVE           :: l_ukca_photin_called = .FALSE.

! DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHOTOL_CTL'

! Error handling for routines called here
CHARACTER(LEN=maxlen_message) :: err_message   ! Error return message

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! -- body of subroutine begins here -- !
error_code_ptr = 0
err_message = ''
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

photol_rates(:,:,:,:) = 0.0
! if photolysis is off exit photol_ctl
IF (photol_config%i_photol_scheme == photolysis_off) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! get number of seconds in parent model timestep as real
r_secs_per_step = REAL(photol_config%chem_timestep)

theta_field_size = row_length * rows

! IF model top height is not provided by parent, derive by subtracting
! planet radius from topmost theta level height
IF ( .NOT. photol_config%l_environ_ztop ) THEN
  ! Crude check for inputs
  IF ( planet_radius <= 0.0 .OR. MAXVAL(r_theta_levels) <= planet_radius ) THEN
    error_code_ptr = errcode_value_invalid
    WRITE(err_message, '(2(A,E10.3))') 'Invalid values for planet_radius ',    &
        planet_radius,' or r_theta_levels ',MAXVAL(r_theta_levels)

    CALL error_report(photol_config%i_error_method, error_code_ptr,            &
                      err_message, RoutineName, msg_out=error_message,         &
                      locn_out=error_routine)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  ! Derive model top. r_theta_levels at top should be horizontally uniform so
  ! any grid cell can be used on that level.
  loc_z_top_model = r_theta_levels(1,1,model_levels) - planet_radius
ELSE
  ! Use parent supplied value
  loc_z_top_model = z_top_of_model
END IF

IF (photol_config%l_strat_chem) THEN
  ! Set up arrays for later strat_photol calculations
  CALL strat_photol_init(theta_field_size)
END IF

! Calculate photolysis rates from Fast-JX
IF (photol_config%i_photol_scheme == photolysis_fastjx) THEN

  ! allocate array to store Fast-JX photolysis rates
  IF (.NOT. ALLOCATED(photol_rates_fastjx)) THEN
    ALLOCATE(photol_rates_fastjx(row_length, rows, model_levels, jppj))
  END IF
  photol_rates_fastjx = 0.0

  ! now call ukca_fastjx to compute photol_rates_fastjx
  CALL ukca_fastjx(error_code_ptr, row_length, rows, model_levels, jppj,       &
                   p_layer_boundaries, t_theta_levels,                         &
                   r_theta_levels, r_rho_levels,                               &
                   longitude, sin_latitude, loc_z_top_model,                   &
                   so4_aitken, so4_accum,                                      &
                   qcl, qcf, area_cloud_fraction,                              &
                   conv_cloud_lwp, conv_cloud_top, conv_cloud_base,            &
                   conv_cloud_amount, aod_sulph_aitk, aod_sulph_accum,         &
                   surf_albedo, ozone_mmr, land_fraction, current_time,        &
                   photol_rates_fastjx, error_message=error_message,           &
                   error_routine=error_routine)

  IF (error_code_ptr > 0) THEN
    IF (ALLOCATED(photol_rates_fastjx)) DEALLOCATE(photol_rates_fastjx)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                            zhook_handle)
    RETURN
  END IF

END IF

! Do prerequisite calculations for call to strat_photol
IF (photol_config%l_strat_chem) THEN
  ! Calculate cos_zenith_angle
  CALL photol_solang(sin_declination, seconds_since_midnight,                  &
                 r_secs_per_step, equation_of_time,                            &
                 sin_latitude, cos_latitude, longitude, theta_field_size,      &
                 cos_zenith_angle)

  ! Calculate ozone column
  CALL photol_calc_ozonecol(error_code_ptr, row_length, rows, model_levels,    &
                   loc_z_top_model, p_layer_boundaries, p_theta_levels,        &
                   ozone_mmr/c_o3, ozonecol, error_message=error_message,      &
                   error_routine=error_routine)
  IF (error_code_ptr > 0 ) THEN
    IF (ALLOCATED(photol_rates_fastjx)) DEALLOCATE(photol_rates_fastjx)
    CALL strat_photol_dealloc()

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                            zhook_handle)
    RETURN
  END IF
END IF

! Fill photol_rates array to pass to UKCA
DO k = 1, model_levels

  ! reset photol2d to zero here
  photol2d = 0.0

  ! If using 2-D (tabulated) scheme pass in values from environment field
  IF (photol_config%i_photol_scheme == photolysis_2d) THEN
    photol2d(:,:,:) = photol_rates_2d(:,:,k,:)
  END IF

  IF (photol_config%l_strat_chem) THEN

    IF (photol_config%i_photol_scheme == photolysis_fastjx) THEN
      ! if any pressure point on domain is below cutoff value call
      ! strat_photol
      photol_strat = 0.0

      IF (MINVAL(p_theta_levels(:,:,k)) < photol_config%fastjx_prescutoff      &
          .AND. photol_config%fastjx_mode /= fjx_mode_fastjx) THEN

        CALL strat_photol(row_length, rows, ratj_data,                         &
              p_theta_levels(:,:,k), t_theta_levels(:,:,k), ozonecol(:,:,k),   &
              cos_zenith_angle, current_time, photol_strat)
      END IF

      ! fill in photol2d if pressure is below minimum
      ! and depending on fastjx_mode
      DO i = 1, rows
        DO j = 1, row_length
          IF (p_theta_levels(j,i,k) < photol_config%fastjx_prescutoff) THEN
            IF (photol_config%fastjx_mode == fjx_mode_2Donly) THEN
              ! fastjx_mode == 1: from strat_photol
              photol2d(j,i,:) = photol_strat(j,i,:)
            ELSE IF (photol_config%fastjx_mode == fjx_mode_merged) THEN
              ! fastjx_mode == 2: from strat_photol & Fast-JX
              photol2d(j,i,:) =                                                &
                photol_strat(j,i,:) + photol_rates_fastjx(j,i,k,:)
            ELSE IF (photol_config%fastjx_mode == fjx_mode_fastjx) THEN
              ! fastjx_mode == 3: from Fast-JX
              photol2d(j,i,:) = photol_rates_fastjx(j,i,k,:)
            ELSE
              ! no other options should be available
              error_code_ptr = errcode_value_invalid
              WRITE(err_message,'(A,I0)')'Incorrect option for fastjx_mode: ', &
                photol_config%fastjx_mode
              IF ( ALLOCATED(photol_rates_fastjx) )                            &
                DEALLOCATE( photol_rates_fastjx )
              CALL strat_photol_dealloc()

              CALL error_report(photol_config%i_error_method, error_code_ptr,  &
                     err_message, RoutineName, msg_out=error_message,          &
                     locn_out=error_routine)

              IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, &
                                      zhook_handle)
              RETURN
            END IF
          ELSE
            ! only take Fast-JX here
            photol2d(j,i,:) = photol_rates_fastjx(j,i,k,:)
          END IF
        END DO
      END DO
    ELSE
      ! apply stratospheric photolysis to previously calculated photol2d
      ! (either 2D photolysis rates or zero if
      ! i_ukca_photol == ukca_photolysis_strat_only)
      CALL strat_photol(row_length, rows, ratj_data,                           &
             p_theta_levels(:,:,k), t_theta_levels(:,:,k), ozonecol(:,:,k),    &
             cos_zenith_angle, current_time, photol2d)
    END IF
  ELSE
    ! tropospheric chemistry selected. Use Fast-JX rates or 2D photolysis
    ! rates. If using Fast-JX, overwrite with Fast-JX values
    IF (photol_config%i_photol_scheme == photolysis_fastjx) THEN
      photol2d = photol_rates_fastjx(:,:,k,:)
    END IF
  END IF ! l_strat_chem

  ! construct photol_rates array to pass to UKCA
  photol_rates(:,:,k,:) = photol2d(:,:,:)
END DO

! Deallocate remaining local arrays
IF (ALLOCATED(photol_rates_fastjx)) DEALLOCATE(photol_rates_fastjx)

IF (photol_config%l_strat_chem) THEN
  ! Deallocate arrays used for strat_photol calculations
  CALL strat_photol_dealloc()
END IF

! If required, use photolysis rates from external sources (rad_ctl)
! to overwrite photol_rates of species jo2 and jo2b
IF (photol_config%l_environ_jo2) THEN
  DO i = 1, jppj
      ! O2 + hv --> O3P + O3P
    IF (ratj_varnames(i) == 'jo2       ') THEN
      photol_rates(:,:,:,i) = rad_ctl_jo2
    END IF
  END DO
END IF
IF (photol_config%l_environ_jo2b) THEN
  DO i = 1, jppj
      ! O2 + hv --> O3P + O3P
    IF (ratj_varnames(i) == 'jo2b       ') THEN
      photol_rates(:,:,:,i) = rad_ctl_jo2b
    END IF
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE photol_ctl

END MODULE photol_ctl_mod
