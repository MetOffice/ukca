! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
! Module to extend environment field handling to a reduced dimension
! domain. Supports native vector and scalar drivers for a single column
! model.
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

MODULE ukca_environment_rdim_mod

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

USE ukca_fieldname_mod,  ONLY:                                                 &
  fldname_latitude,                                                            &
  fldname_longitude,                                                           &
  fldname_sin_latitude,                                                        &
  fldname_cos_latitude,                                                        &
  fldname_tan_latitude,                                                        &
  fldname_conv_cloud_lwp,                                                      &
  fldname_tstar,                                                               &
  fldname_zbl,                                                                 &
  fldname_rough_length,                                                        &
  fldname_seaice_frac,                                                         &
  fldname_pstar,                                                               &
  fldname_surf_albedo,                                                         &
  fldname_zhsc,                                                                &
  fldname_u_scalar_10m,                                                        &
  fldname_surf_hf,                                                             &
  fldname_u_s,                                                                 &
  fldname_ch4_wetl_emiss,                                                      &
  fldname_dms_sea_conc,                                                        &
  fldname_chloro_sea,                                                          &
  fldname_dust_flux_div1,                                                      &
  fldname_dust_flux_div2,                                                      &
  fldname_dust_flux_div3,                                                      &
  fldname_dust_flux_div4,                                                      &
  fldname_dust_flux_div5,                                                      &
  fldname_dust_flux_div6,                                                      &
  fldname_surf_wetness,                                                        &
  fldname_kent,                                                                &
  fldname_kent_dsc,                                                            &
  fldname_conv_cloud_base,                                                     &
  fldname_conv_cloud_top,                                                      &
  fldname_ext_cg_flash,                                                        &
  fldname_ext_ic_flash,                                                        &
  fldname_land_sea_mask,                                                       &
  fldname_theta,                                                               &
  fldname_q,                                                                   &
  fldname_qcf,                                                                 &
  fldname_conv_cloud_amount,                                                   &
  fldname_rho_r2,                                                              &
  fldname_qcl,                                                                 &
  fldname_exner_rho_levels,                                                    &
  fldname_area_cloud_fraction,                                                 &
  fldname_cloud_frac,                                                          &
  fldname_cloud_liq_frac,                                                      &
  fldname_exner_theta_levels,                                                  &
  fldname_p_rho_levels,                                                        &
  fldname_p_theta_levels,                                                      &
  fldname_rhokh_rdz,                                                           &
  fldname_dtrdz,                                                               &
  fldname_we_lim,                                                              &
  fldname_t_frac,                                                              &
  fldname_zrzi,                                                                &
  fldname_we_lim_dsc,                                                          &
  fldname_t_frac_dsc,                                                          &
  fldname_zrzi_dsc,                                                            &
  fldname_stcon,                                                               &
  fldname_ls_rain3d,                                                           &
  fldname_ls_snow3d,                                                           &
  fldname_autoconv,                                                            &
  fldname_accretion,                                                           &
  fldname_pv_on_theta_mlevs,                                                   &
  fldname_conv_rain3d,                                                         &
  fldname_conv_snow3d,                                                         &
  fldname_so4_sa_clim,                                                         &
  fldname_so4_aitken,                                                          &
  fldname_so4_accum,                                                           &
  fldname_soot_fresh,                                                          &
  fldname_soot_aged,                                                           &
  fldname_ocff_fresh,                                                          &
  fldname_ocff_aged,                                                           &
  fldname_biogenic,                                                            &
  fldname_sea_salt_film,                                                       &
  fldname_sea_salt_jet,                                                        &
  fldname_co2_interactive,                                                     &
  fldname_rim_cry,                                                             &
  fldname_rim_agg,                                                             &
  fldname_vertvel,                                                             &
  fldname_bl_tke,                                                              &
  fldname_h2o2_offline,                                                        &
  fldname_ho2_offline,                                                         &
  fldname_no3_offline,                                                         &
  fldname_o3_offline,                                                          &
  fldname_oh_offline,                                                          &
  fldname_dust_div1,                                                           &
  fldname_dust_div2,                                                           &
  fldname_dust_div3,                                                           &
  fldname_dust_div4,                                                           &
  fldname_dust_div5,                                                           &
  fldname_dust_div6,                                                           &
  fldname_interf_z,                                                            &
  fldname_grid_surf_area,                                                      &
  fldname_lscat_zhang,                                                         &
  fldname_photol_rates,                                                        &
  fldname_grid_area_fullht

USE ukca_environment_fields_mod, ONLY:                                         &
  environ_field_info,                                                          &
  k1_dust_flux, k2_dust_flux,                                                  &
  locate_land_points,                                                          &
  clear_land_only_fields,                                                      &
  latitude,                                                                    &
  longitude,                                                                   &
  sin_latitude,                                                                &
  cos_latitude,                                                                &
  tan_latitude,                                                                &
  conv_cloud_lwp,                                                              &
  tstar,                                                                       &
  zbl,                                                                         &
  rough_length,                                                                &
  seaice_frac,                                                                 &
  pstar,                                                                       &
  surf_albedo,                                                                 &
  zhsc,                                                                        &
  u_scalar_10m,                                                                &
  surf_hf,                                                                     &
  u_s,                                                                         &
  ch4_wetl_emiss,                                                              &
  dms_sea_conc,                                                                &
  chloro_sea,                                                                  &
  dust_flux,                                                                   &
  surf_wetness,                                                                &
  kent,                                                                        &
  kent_dsc,                                                                    &
  conv_cloud_base,                                                             &
  conv_cloud_top,                                                              &
  ext_cg_flash,                                                                &
  ext_ic_flash,                                                                &
  land_sea_mask,                                                               &
  theta,                                                                       &
  q,                                                                           &
  qcf,                                                                         &
  conv_cloud_amount,                                                           &
  rho_r2,                                                                      &
  qcl,                                                                         &
  exner_rho_levels,                                                            &
  area_cloud_fraction,                                                         &
  cloud_frac,                                                                  &
  cloud_liq_frac,                                                              &
  exner_theta_levels,                                                          &
  p_rho_levels,                                                                &
  p_theta_levels,                                                              &
  rhokh_rdz,                                                                   &
  dtrdz,                                                                       &
  we_lim,                                                                      &
  t_frac,                                                                      &
  zrzi,                                                                        &
  we_lim_dsc,                                                                  &
  t_frac_dsc,                                                                  &
  zrzi_dsc,                                                                    &
  stcon,                                                                       &
  ls_rain3d,                                                                   &
  ls_snow3d,                                                                   &
  autoconv,                                                                    &
  accretion,                                                                   &
  pv_on_theta_mlevs,                                                           &
  conv_rain3d,                                                                 &
  conv_snow3d,                                                                 &
  so4_sa_clim,                                                                 &
  so4_aitken,                                                                  &
  so4_accum,                                                                   &
  soot_fresh,                                                                  &
  soot_aged,                                                                   &
  ocff_fresh,                                                                  &
  ocff_aged,                                                                   &
  biogenic,                                                                    &
  sea_salt_film,                                                               &
  sea_salt_jet,                                                                &
  co2_interactive,                                                             &
  rim_cry,                                                                     &
  rim_agg,                                                                     &
  vertvel,                                                                     &
  bl_tke,                                                                      &
  h2o2_offline,                                                                &
  ho2_offline,                                                                 &
  no3_offline,                                                                 &
  o3_offline,                                                                  &
  oh_offline,                                                                  &
  dust_div1,                                                                   &
  dust_div2,                                                                   &
  dust_div3,                                                                   &
  dust_div4,                                                                   &
  dust_div5,                                                                   &
  dust_div6,                                                                   &
  interf_z,                                                                    &
  grid_surf_area,                                                              &
  lscat_zhang,                                                                 &
  photol_rates,                                                                &
  grid_area_fullht

USE ukca_error_mod,  ONLY: errcode_env_field_unknown,                          &
                           errcode_env_field_mismatch

IMPLICIT NONE

PRIVATE

! Public procedures
PUBLIC set_env_2d_from_0d_real, set_env_2d_from_0d_integer,                    &
       set_env_2d_from_0d_logical, set_env_3d_from_1d_real,                    &
       set_env_4d_from_2d_real

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_ENVIRONMENT_RDIM_MOD'


CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE set_env_2d_from_0d_real(varname, i_field, field_data, error_code)
! ----------------------------------------------------------------------
! Description:
!   Set a 2D environment field from a scalar input value of type real.
!
! Method:
!   The field to be set is identified by 'varname'. The argument 'i_field'
!   is expected to give its position in the list of required fields or be
!   zero if it is not a required field.
!   Irrespective of whether the field is required, a non-zero error code
!   is returned if 'varname' does not refer to a valid 2D field that can be
!   set in this this way. Valid fields are all defined internally on a 2D
!   horizontal spatial grid.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: i_field
REAL, INTENT(IN) :: field_data
INTEGER, INTENT(OUT) :: error_code

! Local variables

INTEGER :: i1
INTEGER :: i2
INTEGER :: j1
INTEGER :: j2

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_ENV_2D_FROM_0D_REAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0

! Get bounds for the 2D internal field
IF (i_field /= 0) THEN
  i1 = environ_field_info(i_field)%lbound_dim1
  i2 = environ_field_info(i_field)%ubound_dim1
  j1 = environ_field_info(i_field)%lbound_dim2
  j2 = environ_field_info(i_field)%ubound_dim2
END IF

! Copy real field value to the appropriate UKCA internal array if required
SELECT CASE (varname)
CASE (fldname_latitude)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(latitude)) ALLOCATE(latitude(i1:i2,j1:j2))
    latitude = field_data
  END IF
CASE (fldname_longitude)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(longitude)) ALLOCATE(longitude(i1:i2,j1:j2))
    longitude = field_data
  END IF
CASE (fldname_sin_latitude)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(sin_latitude)) ALLOCATE(sin_latitude(i1:i2,j1:j2))
    sin_latitude = field_data
  END IF
CASE (fldname_cos_latitude)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(cos_latitude)) ALLOCATE(cos_latitude(i1:i2,j1:j2))
    cos_latitude = field_data
  END IF
CASE (fldname_tan_latitude)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(tan_latitude)) ALLOCATE(tan_latitude(i1:i2,j1:j2))
    tan_latitude = field_data
  END IF
CASE (fldname_conv_cloud_lwp)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(conv_cloud_lwp)) ALLOCATE(conv_cloud_lwp(i1:i2,j1:j2))
    conv_cloud_lwp = field_data
  END IF
CASE (fldname_tstar)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(tstar)) ALLOCATE(tstar(i1:i2,j1:j2))
    tstar = field_data
  END IF
CASE (fldname_zbl)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(zbl)) ALLOCATE(zbl(i1:i2,j1:j2))
    zbl = field_data
  END IF
CASE (fldname_rough_length)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(rough_length)) ALLOCATE(rough_length(i1:i2,j1:j2))
    rough_length = field_data
  END IF
CASE (fldname_seaice_frac)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(seaice_frac)) ALLOCATE(seaice_frac(i1:i2,j1:j2))
    seaice_frac = field_data
  END IF
CASE (fldname_pstar)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(pstar)) ALLOCATE(pstar(i1:i2,j1:j2))
    pstar = field_data
  END IF
CASE (fldname_surf_albedo)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(surf_albedo)) ALLOCATE(surf_albedo(i1:i2,j1:j2))
    surf_albedo = field_data
  END IF
CASE (fldname_zhsc)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(zhsc)) ALLOCATE(zhsc(i1:i2,j1:j2))
    zhsc = field_data
  END IF
CASE (fldname_u_scalar_10m)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(u_scalar_10m)) ALLOCATE(u_scalar_10m(i1:i2,j1:j2))
    u_scalar_10m = field_data
  END IF
CASE (fldname_surf_hf)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(surf_hf)) ALLOCATE(surf_hf(i1:i2,j1:j2))
    surf_hf = field_data
  END IF
CASE (fldname_u_s)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(u_s)) ALLOCATE(u_s(i1:i2,j1:j2))
    u_s = field_data
  END IF
CASE (fldname_ch4_wetl_emiss)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ch4_wetl_emiss)) ALLOCATE(ch4_wetl_emiss(i1:i2,j1:j2))
    ch4_wetl_emiss = field_data
  END IF
CASE (fldname_dms_sea_conc)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dms_sea_conc)) ALLOCATE(dms_sea_conc(i1:i2,j1:j2))
    dms_sea_conc = field_data
  END IF
CASE (fldname_chloro_sea)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(chloro_sea)) ALLOCATE(chloro_sea(i1:i2,j1:j2))
    chloro_sea = field_data
  END IF
CASE (fldname_dust_flux_div1)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,1) = field_data
  END IF
CASE (fldname_dust_flux_div2)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,2) = field_data
  END IF
CASE (fldname_dust_flux_div3)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,3) = field_data
  END IF
CASE (fldname_dust_flux_div4)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,4) = field_data
  END IF
CASE (fldname_dust_flux_div5)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,5) = field_data
  END IF
CASE (fldname_dust_flux_div6)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,6) = field_data
  END IF
CASE (fldname_surf_wetness)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(surf_wetness)) ALLOCATE(surf_wetness(i1:i2,j1:j2))
    surf_wetness = field_data
  END IF
CASE (fldname_grid_surf_area)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(grid_surf_area)) ALLOCATE(grid_surf_area(i1:i2,j1:j2))
    grid_surf_area = field_data
  END IF
CASE (fldname_ext_cg_flash)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ext_cg_flash)) ALLOCATE(ext_cg_flash(i1:i2,j1:j2))
    ext_cg_flash = field_data
  END IF
CASE (fldname_ext_ic_flash)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ext_ic_flash)) ALLOCATE(ext_ic_flash(i1:i2,j1:j2))
    ext_ic_flash = field_data
  END IF
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_env_2d_from_0d_real


! ----------------------------------------------------------------------
SUBROUTINE set_env_2d_from_0d_integer(varname, i_field, field_data, error_code)
! ----------------------------------------------------------------------
! Description:
!   Set a 2D environment field from a scalar input value of type integer.
!
! Method:
!   See set_env_2d_from_0d_real.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: i_field
INTEGER, INTENT(IN) :: field_data
INTEGER, INTENT(OUT) :: error_code

! Local variables

INTEGER :: i1
INTEGER :: i2
INTEGER :: j1
INTEGER :: j2

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_ENV_2D_FROM_0D_INTEGER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0

! Get bounds for the 2D internal field
IF (i_field /= 0) THEN
  i1 = environ_field_info(i_field)%lbound_dim1
  i2 = environ_field_info(i_field)%ubound_dim1
  j1 = environ_field_info(i_field)%lbound_dim2
  j2 = environ_field_info(i_field)%ubound_dim2
END IF

! Copy integer field value to the appropriate UKCA internal array if required
SELECT CASE (varname)
CASE (fldname_kent)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(kent)) ALLOCATE(kent(i1:i2,j1:j2))
    kent = field_data
  END IF
CASE (fldname_kent_dsc)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(kent_dsc)) ALLOCATE(kent_dsc(i1:i2,j1:j2))
    kent_dsc = field_data
  END IF
CASE (fldname_conv_cloud_base)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(conv_cloud_base)) ALLOCATE(conv_cloud_base(i1:i2,j1:j2))
    conv_cloud_base = field_data
  END IF
CASE (fldname_conv_cloud_top)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(conv_cloud_top)) ALLOCATE(conv_cloud_top(i1:i2,j1:j2))
    conv_cloud_top = field_data
  END IF
CASE (fldname_lscat_zhang)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(lscat_zhang)) ALLOCATE(lscat_zhang(i1:i2,j1:j2))
    lscat_zhang = field_data
  END IF
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_env_2d_from_0d_integer


! ----------------------------------------------------------------------
SUBROUTINE set_env_2d_from_0d_logical(varname, i_field, field_data, error_code)
! ----------------------------------------------------------------------
! Description:
! Set a 2D environment field from a scalar input value of type logical.
!
! Method:
!   See set_env_2d_from_0d_real.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: i_field
LOGICAL, INTENT(IN) :: field_data
INTEGER, INTENT(OUT) :: error_code

! Local variables

INTEGER :: i1
INTEGER :: i2
INTEGER :: j1
INTEGER :: j2

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_ENV_2D_FROM_0D_LOGICAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0

! Get bounds for the 2D internal field
IF (i_field /= 0) THEN
  i1 = environ_field_info(i_field)%lbound_dim1
  i2 = environ_field_info(i_field)%ubound_dim1
  j1 = environ_field_info(i_field)%lbound_dim2
  j2 = environ_field_info(i_field)%ubound_dim2
END IF

! Copy logical field value to the appropriate UKCA internal array if required
SELECT CASE (varname)
CASE (fldname_land_sea_mask)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(land_sea_mask)) ALLOCATE(land_sea_mask(i1:i2,j1:j2))
    land_sea_mask = field_data
  END IF
  ! Special processing to set up index of land points for locating
  ! land-only environment fields on a 2-D grid
  CALL locate_land_points()
  ! Clear any existing land-only fields as these fields may be inconsistent
  ! with the new land sea mask
  CALL clear_land_only_fields()
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_env_2d_from_0d_logical


! ----------------------------------------------------------------------
SUBROUTINE set_env_3d_from_1d_real(varname, i_field, field_data, error_code)
! ----------------------------------------------------------------------
! Description:
! Set a 3D environment field from a 1D input field of type real.
!
! Method:
!   The field to be set is identified by 'varname'. The argument 'i_field'
!   is expected to give its position in the list of required fields.
!   If i_field is zero, control is returned without any action.
!   A non-zero error code is returned if 'varname' does not refer to a valid 3D
!   field that can be set in this this way.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: i_field
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:)
INTEGER, INTENT(OUT) :: error_code

! Local variables

INTEGER :: i1
INTEGER :: i2
INTEGER :: j1
INTEGER :: j2
INTEGER :: k1
INTEGER :: k2
INTEGER :: i
INTEGER :: j

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_ENV_3D_FROM_1D_REAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0

! The routine should not be called if i_field is 0, but check anyway.
IF (i_field == 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Get bounds for the 3D internal field
i1 = environ_field_info(i_field)%lbound_dim1
i2 = environ_field_info(i_field)%ubound_dim1
j1 = environ_field_info(i_field)%lbound_dim2
j2 = environ_field_info(i_field)%ubound_dim2
k1 = environ_field_info(i_field)%lbound_dim3
k2 = environ_field_info(i_field)%ubound_dim3

! Copy 1D real field to the appropriate UKCA internal array.
! Any data outside the required bounds are discarded.
SELECT CASE (varname)
CASE (fldname_stcon)
  IF (.NOT. ALLOCATED(stcon)) ALLOCATE(stcon(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      stcon(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_theta)
  IF (.NOT. ALLOCATED(theta)) ALLOCATE(theta(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      theta(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_q)
  IF (.NOT. ALLOCATED(q)) ALLOCATE(q(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      q(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_qcf)
  IF (.NOT. ALLOCATED(qcf)) ALLOCATE(qcf(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      qcf(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_conv_cloud_amount)
  IF (.NOT. ALLOCATED(conv_cloud_amount))                                      &
    ALLOCATE(conv_cloud_amount(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      conv_cloud_amount(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_rho_r2)
  IF (.NOT. ALLOCATED(rho_r2)) ALLOCATE(rho_r2(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      rho_r2(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_qcl)
  IF (.NOT. ALLOCATED(qcl)) ALLOCATE(qcl(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      qcl(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_exner_rho_levels)
  IF (.NOT. ALLOCATED(exner_rho_levels))                                       &
    ALLOCATE(exner_rho_levels(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      exner_rho_levels(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_area_cloud_fraction)
  IF (.NOT. ALLOCATED(area_cloud_fraction))                                    &
    ALLOCATE(area_cloud_fraction(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      area_cloud_fraction(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_cloud_frac)
  IF (.NOT. ALLOCATED(cloud_frac)) ALLOCATE(cloud_frac(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      cloud_frac(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_cloud_liq_frac)
  IF (.NOT. ALLOCATED(cloud_liq_frac))                                         &
    ALLOCATE(cloud_liq_frac(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      cloud_liq_frac(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_exner_theta_levels)
  IF (.NOT. ALLOCATED(exner_theta_levels))                                     &
    ALLOCATE(exner_theta_levels(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      exner_theta_levels(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_p_rho_levels)
  IF (.NOT. ALLOCATED(p_rho_levels)) ALLOCATE(p_rho_levels(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      p_rho_levels(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_p_theta_levels)
  IF (.NOT. ALLOCATED(p_theta_levels))                                         &
    ALLOCATE(p_theta_levels(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      p_theta_levels(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_rhokh_rdz)
  IF (.NOT. ALLOCATED(rhokh_rdz)) ALLOCATE(rhokh_rdz(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      rhokh_rdz(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_dtrdz)
  IF (.NOT. ALLOCATED(dtrdz)) ALLOCATE(dtrdz(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      dtrdz(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_we_lim)
  IF (.NOT. ALLOCATED(we_lim)) ALLOCATE(we_lim(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      we_lim(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_t_frac)
  IF (.NOT. ALLOCATED(t_frac)) ALLOCATE(t_frac(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      t_frac(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_zrzi)
  IF (.NOT. ALLOCATED(zrzi)) ALLOCATE(zrzi(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      zrzi(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_we_lim_dsc)
  IF (.NOT. ALLOCATED(we_lim_dsc)) ALLOCATE(we_lim_dsc(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      we_lim_dsc(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_t_frac_dsc)
  IF (.NOT. ALLOCATED(t_frac_dsc)) ALLOCATE(t_frac_dsc(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      t_frac_dsc(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_zrzi_dsc)
  IF (.NOT. ALLOCATED(zrzi_dsc)) ALLOCATE(zrzi_dsc(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      zrzi_dsc(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_ls_rain3d)
  IF (.NOT. ALLOCATED(ls_rain3d)) ALLOCATE(ls_rain3d(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      ls_rain3d(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_ls_snow3d)
  IF (.NOT. ALLOCATED(ls_snow3d)) ALLOCATE(ls_snow3d(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      ls_snow3d(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_autoconv)
  IF (.NOT. ALLOCATED(autoconv)) ALLOCATE(autoconv(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      autoconv(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_accretion)
  IF (.NOT. ALLOCATED(accretion)) ALLOCATE(accretion(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      accretion(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_pv_on_theta_mlevs)
  IF (.NOT. ALLOCATED(pv_on_theta_mlevs))                                      &
    ALLOCATE(pv_on_theta_mlevs(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      pv_on_theta_mlevs(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_conv_rain3d)
  IF (.NOT. ALLOCATED(conv_rain3d)) ALLOCATE(conv_rain3d(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      conv_rain3d(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_conv_snow3d)
  IF (.NOT. ALLOCATED(conv_snow3d)) ALLOCATE(conv_snow3d(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      conv_snow3d(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_so4_sa_clim)
  IF (.NOT. ALLOCATED(so4_sa_clim)) ALLOCATE(so4_sa_clim(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      so4_sa_clim(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_so4_aitken)
  IF (.NOT. ALLOCATED(so4_aitken)) ALLOCATE(so4_aitken(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      so4_aitken(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_so4_accum)
  IF (.NOT. ALLOCATED(so4_accum)) ALLOCATE(so4_accum(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      so4_accum(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_soot_fresh)
  IF (.NOT. ALLOCATED(soot_fresh)) ALLOCATE(soot_fresh(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      soot_fresh(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_soot_aged)
  IF (.NOT. ALLOCATED(soot_aged)) ALLOCATE(soot_aged(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      soot_aged(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_ocff_fresh)
  IF (.NOT. ALLOCATED(ocff_fresh)) ALLOCATE(ocff_fresh(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      ocff_fresh(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_ocff_aged)
  IF (.NOT. ALLOCATED(ocff_aged)) ALLOCATE(ocff_aged(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      ocff_aged(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_biogenic)
  IF (.NOT. ALLOCATED(biogenic)) ALLOCATE(biogenic(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      biogenic(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_dust_div1)
  IF (.NOT. ALLOCATED(dust_div1)) ALLOCATE(dust_div1(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      dust_div1(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_dust_div2)
  IF (.NOT. ALLOCATED(dust_div2)) ALLOCATE(dust_div2(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      dust_div2(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_dust_div3)
  IF (.NOT. ALLOCATED(dust_div3)) ALLOCATE(dust_div3(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      dust_div3(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_dust_div4)
  IF (.NOT. ALLOCATED(dust_div4)) ALLOCATE(dust_div4(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      dust_div4(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_dust_div5)
  IF (.NOT. ALLOCATED(dust_div5)) ALLOCATE(dust_div5(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      dust_div5(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_dust_div6)
  IF (.NOT. ALLOCATED(dust_div6)) ALLOCATE(dust_div6(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      dust_div6(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_sea_salt_film)
  IF (.NOT. ALLOCATED(sea_salt_film)) ALLOCATE(sea_salt_film(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      sea_salt_film(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_sea_salt_jet)
  IF (.NOT. ALLOCATED(sea_salt_jet)) ALLOCATE(sea_salt_jet(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      sea_salt_jet(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_co2_interactive)
  IF (.NOT. ALLOCATED(co2_interactive))                                        &
    ALLOCATE(co2_interactive(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      co2_interactive(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_rim_cry)
  IF (.NOT. ALLOCATED(rim_cry)) ALLOCATE(rim_cry(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      rim_cry(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_rim_agg)
  IF (.NOT. ALLOCATED(rim_agg)) ALLOCATE(rim_agg(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      rim_agg(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_vertvel)
  IF (.NOT. ALLOCATED(vertvel)) ALLOCATE(vertvel(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      vertvel(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_bl_tke)
  IF (.NOT. ALLOCATED(bl_tke)) ALLOCATE(bl_tke(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      bl_tke(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_interf_z)
  IF (.NOT. ALLOCATED(interf_z)) ALLOCATE(interf_z(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      interf_z(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_h2o2_offline)
  IF (.NOT. ALLOCATED(h2o2_offline)) ALLOCATE(h2o2_offline(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      h2o2_offline(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_ho2_offline)
  IF (.NOT. ALLOCATED(ho2_offline)) ALLOCATE(ho2_offline(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      ho2_offline(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_no3_offline)
  IF (.NOT. ALLOCATED(no3_offline)) ALLOCATE(no3_offline(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      no3_offline(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_o3_offline)
  IF (.NOT. ALLOCATED(o3_offline)) ALLOCATE(o3_offline(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      o3_offline(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_oh_offline)
  IF (.NOT. ALLOCATED(oh_offline)) ALLOCATE(oh_offline(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      oh_offline(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE (fldname_grid_area_fullht)
  IF (.NOT. ALLOCATED(grid_area_fullht))                                     &
    ALLOCATE(grid_area_fullht(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      grid_area_fullht(i,j,:) = field_data(k1:k2)
    END DO
  END DO
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_env_3d_from_1d_real

! ----------------------------------------------------------------------
SUBROUTINE set_env_4d_from_2d_real(varname, i_field, field_data, error_code)
! ----------------------------------------------------------------------
! Description:
! Set a 4D environment field from a 2D input field of type real.
!
! Method:
!   The field to be set is identified by 'varname'. The argument 'i_field'
!   is expected to give its position in the list of required fields.
!   If i_field is zero control is returned wthout any action.
!   A non-zero error code is returned if 'varname' does not refer to a valid 4D
!   field that can be set in this this way.
!   Currently the only 4-D environment variable active is 'photol_rates'
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: i_field
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:,:)
INTEGER, INTENT(OUT) :: error_code

! Local variables

INTEGER :: i1
INTEGER :: i2
INTEGER :: j1
INTEGER :: j2
INTEGER :: k1
INTEGER :: k2
INTEGER :: n1
INTEGER :: n2
INTEGER :: i
INTEGER :: j
INTEGER :: k

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_ENV_4D_FROM_2D_REAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0

! The routine should not called if i_field is 0, but check anyway.
IF (i_field == 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Get bounds for the 4D internal field
i1 = environ_field_info(i_field)%lbound_dim1
i2 = environ_field_info(i_field)%ubound_dim1
j1 = environ_field_info(i_field)%lbound_dim2
j2 = environ_field_info(i_field)%ubound_dim2
k1 = environ_field_info(i_field)%lbound_dim3
k2 = environ_field_info(i_field)%ubound_dim3
n1 = environ_field_info(i_field)%lbound_dim4
n2 = environ_field_info(i_field)%ubound_dim4

! Copy 2D real field to the appropriate UKCA internal array (currently only
!  photol_rates). Any data outside the required bounds are discarded.
SELECT CASE (varname)
CASE (fldname_photol_rates)
  IF (.NOT. ALLOCATED(photol_rates))                                           &
    ALLOCATE(photol_rates(i1:i2,j1:j2,k1:k2,n1:n2))
  DO i = i1,i2
    DO j = j1,j2
      DO k = k1,k2
        photol_rates(i,j,k,:) = field_data(k,n1:n2)
      END DO
    END DO
  END DO
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_env_4d_from_2d_real

END MODULE ukca_environment_rdim_mod

