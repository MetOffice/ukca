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
  fldname_grid_area_fullht,                                                    &
  fldname_grid_volume

USE ukca_environment_fields_mod, ONLY:                                         &
  environ_field_info,                                                          &
  environ_field_ptrs,                                                          &
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
  grid_area_fullht,                                                            &
  grid_volume

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
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 latitude)
CASE (fldname_longitude)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 longitude)
CASE (fldname_sin_latitude)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 sin_latitude)
CASE (fldname_cos_latitude)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 cos_latitude)
CASE (fldname_tan_latitude)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 tan_latitude)
CASE (fldname_conv_cloud_lwp)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 conv_cloud_lwp)
CASE (fldname_tstar)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 tstar)
CASE (fldname_zbl)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 zbl)
CASE (fldname_rough_length)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 rough_length)
CASE (fldname_seaice_frac)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 seaice_frac)
CASE (fldname_pstar)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 pstar)
CASE (fldname_surf_albedo)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 surf_albedo)
CASE (fldname_zhsc)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 zhsc)
CASE (fldname_u_scalar_10m)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 u_scalar_10m)
CASE (fldname_surf_hf)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 surf_hf)
CASE (fldname_u_s)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 u_s)
CASE (fldname_ch4_wetl_emiss)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 ch4_wetl_emiss)
CASE (fldname_dms_sea_conc)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 dms_sea_conc)
CASE (fldname_chloro_sea)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 chloro_sea)
  ! Dust fluxes are stored in a 3D array so need special treatment
CASE (fldname_dust_flux_div1)
  CALL set_field_2d_from_0d_real_k(i_field, i1, i2, j1, j2,                    &
                                   k1_dust_flux, k2_dust_flux, 1,              &
                                   field_data, dust_flux)
CASE (fldname_dust_flux_div2)
  CALL set_field_2d_from_0d_real_k(i_field, i1, i2, j1, j2,                    &
                                   k1_dust_flux, k2_dust_flux, 2,              &
                                   field_data, dust_flux)
CASE (fldname_dust_flux_div3)
  CALL set_field_2d_from_0d_real_k(i_field, i1, i2, j1, j2,                    &
                                   k1_dust_flux, k2_dust_flux, 3,              &
                                   field_data, dust_flux)
CASE (fldname_dust_flux_div4)
  CALL set_field_2d_from_0d_real_k(i_field, i1, i2, j1, j2,                    &
                                   k1_dust_flux, k2_dust_flux, 4,              &
                                   field_data, dust_flux)
CASE (fldname_dust_flux_div5)
  CALL set_field_2d_from_0d_real_k(i_field, i1, i2, j1, j2,                    &
                                   k1_dust_flux, k2_dust_flux, 5,              &
                                   field_data, dust_flux)
CASE (fldname_dust_flux_div6)
  CALL set_field_2d_from_0d_real_k(i_field, i1, i2, j1, j2,                    &
                                   k1_dust_flux, k2_dust_flux, 6,              &
                                   field_data, dust_flux)
  ! End of dust flux cases
CASE (fldname_surf_wetness)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 surf_wetness)
CASE (fldname_grid_surf_area)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 grid_surf_area)
CASE (fldname_ext_cg_flash)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 ext_cg_flash)
CASE (fldname_ext_ic_flash)
  CALL set_field_2d_from_0d_real(i_field, i1, i2, j1, j2, field_data,          &
                                 ext_ic_flash)
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
  CALL set_field_2d_from_0d_integer(i_field, i1, i2, j1, j2, field_data,       &
                                    kent)
CASE (fldname_kent_dsc)
  CALL set_field_2d_from_0d_integer(i_field, i1, i2, j1, j2, field_data,       &
                                    kent_dsc)
CASE (fldname_conv_cloud_base)
  CALL set_field_2d_from_0d_integer(i_field, i1, i2, j1, j2, field_data,       &
                                    conv_cloud_base)
CASE (fldname_conv_cloud_top)
  CALL set_field_2d_from_0d_integer(i_field, i1, i2, j1, j2, field_data,       &
                                    conv_cloud_top)
CASE (fldname_lscat_zhang)
  CALL set_field_2d_from_0d_integer(i_field, i1, i2, j1, j2, field_data,       &
                                    lscat_zhang)
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
    CALL set_field_2d_from_0d_logical(i_field, i1, i2, j1, j2, field_data,     &
                                      land_sea_mask)
    ! Special processing to set up index of land points for locating
    ! land-only environment fields on a 2-D grid
    CALL locate_land_points()
    ! Clear any existing land-only fields as these fields may be inconsistent
    ! with the new land sea mask
    CALL clear_land_only_fields()
  END IF
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
!   is expected to give its position in the list of required fields or may
!   be zero if it is not a required field. (Note that the latter case is not
!   expected but is allowed for to be consistent with other subroutines.)
!   Irrespective of whether the field is required, a non-zero error code
!   is returned if 'varname' does not refer to a valid 3D field that can be
!   set in this this way.
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

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_ENV_3D_FROM_1D_REAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0

IF (i_field /= 0) THEN
  ! Get bounds for the 3D internal field
  i1 = environ_field_info(i_field)%lbound_dim1
  i2 = environ_field_info(i_field)%ubound_dim1
  j1 = environ_field_info(i_field)%lbound_dim2
  j2 = environ_field_info(i_field)%ubound_dim2
  k1 = environ_field_info(i_field)%lbound_dim3
  k2 = environ_field_info(i_field)%ubound_dim3
END IF

! Copy 1D real field to the appropriate UKCA internal array.
! Any data outside the required bounds are discarded.
SELECT CASE (varname)
CASE (fldname_stcon)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 stcon)
CASE (fldname_theta)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 theta)
CASE (fldname_q)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 q)
CASE (fldname_qcf)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 qcf)
CASE (fldname_conv_cloud_amount)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 conv_cloud_amount)
CASE (fldname_rho_r2)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 rho_r2)
CASE (fldname_qcl)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 qcl)
CASE (fldname_exner_rho_levels)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 exner_rho_levels)
CASE (fldname_area_cloud_fraction)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 area_cloud_fraction)
CASE (fldname_cloud_frac)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 cloud_frac)
CASE (fldname_cloud_liq_frac)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 cloud_liq_frac)
CASE (fldname_exner_theta_levels)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 exner_theta_levels)
CASE (fldname_p_rho_levels)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 p_rho_levels)
CASE (fldname_p_theta_levels)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 p_theta_levels)
CASE (fldname_rhokh_rdz)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 rhokh_rdz)
CASE (fldname_dtrdz)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 dtrdz)
CASE (fldname_we_lim)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 we_lim)
CASE (fldname_t_frac)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 t_frac)
CASE (fldname_zrzi)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 zrzi)
CASE (fldname_we_lim_dsc)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 we_lim_dsc)
CASE (fldname_t_frac_dsc)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 t_frac_dsc)
CASE (fldname_zrzi_dsc)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 zrzi_dsc)
CASE (fldname_ls_rain3d)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 ls_rain3d)
CASE (fldname_ls_snow3d)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 ls_snow3d)
CASE (fldname_autoconv)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 autoconv)
CASE (fldname_accretion)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 accretion)
CASE (fldname_pv_on_theta_mlevs)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 pv_on_theta_mlevs)
CASE (fldname_conv_rain3d)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 conv_rain3d)
CASE (fldname_conv_snow3d)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 conv_snow3d)
CASE (fldname_so4_sa_clim)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 so4_sa_clim)
CASE (fldname_so4_aitken)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 so4_aitken)
CASE (fldname_so4_accum)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 so4_accum)
CASE (fldname_soot_fresh)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 soot_fresh)
CASE (fldname_soot_aged)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 soot_aged)
CASE (fldname_ocff_fresh)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 ocff_fresh)
CASE (fldname_ocff_aged)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 ocff_aged)
CASE (fldname_biogenic)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 biogenic)
CASE (fldname_dust_div1)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 dust_div1)
CASE (fldname_dust_div2)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 dust_div2)
CASE (fldname_dust_div3)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 dust_div3)
CASE (fldname_dust_div4)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 dust_div4)
CASE (fldname_dust_div5)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 dust_div5)
CASE (fldname_dust_div6)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 dust_div6)
CASE (fldname_sea_salt_film)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 sea_salt_film)
CASE (fldname_sea_salt_jet)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 sea_salt_jet)
CASE (fldname_co2_interactive)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 co2_interactive)
CASE (fldname_rim_cry)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 rim_cry)
CASE (fldname_rim_agg)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 rim_agg)
CASE (fldname_vertvel)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 vertvel)
CASE (fldname_bl_tke)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 bl_tke)
CASE (fldname_interf_z)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 interf_z)
CASE (fldname_h2o2_offline)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 h2o2_offline)
CASE (fldname_ho2_offline)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 ho2_offline)
CASE (fldname_no3_offline)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 no3_offline)
CASE (fldname_o3_offline)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 o3_offline)
CASE (fldname_oh_offline)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 oh_offline)
CASE (fldname_grid_area_fullht)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 grid_area_fullht)
CASE (fldname_grid_volume)
  CALL set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,  &
                                 grid_volume)
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
!   is expected to give its position in the list of required fields or may
!   be zero if it is not a required field. (Note that the latter case is not
!   expected but is allowed for to be consistent with other subroutines.)
!   Irrespective of whether the field is required, a non-zero error code
!   is returned if 'varname' does not refer to a valid 3D field that can be
!   set in this this way.
!   Currently the only 4-D environment variable active is 'photol_rates'.
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

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_ENV_4D_FROM_2D_REAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0

IF (i_field /= 0) THEN
  ! Get bounds for the 4D internal field
  i1 = environ_field_info(i_field)%lbound_dim1
  i2 = environ_field_info(i_field)%ubound_dim1
  j1 = environ_field_info(i_field)%lbound_dim2
  j2 = environ_field_info(i_field)%ubound_dim2
  k1 = environ_field_info(i_field)%lbound_dim3
  k2 = environ_field_info(i_field)%ubound_dim3
  n1 = environ_field_info(i_field)%lbound_dim4
  n2 = environ_field_info(i_field)%ubound_dim4
END IF

! Copy 2D real field to the appropriate UKCA internal array (currently only
!  photol_rates). Any data outside the required bounds are discarded.
SELECT CASE (varname)
CASE (fldname_photol_rates)
  CALL set_field_4d_from_2d_real(i_field, i1, i2, j1, j2, k1, k2, n1, n2,      &
                                 field_data, photol_rates)
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_env_4d_from_2d_real


! ----------------------------------------------------------------------
SUBROUTINE set_field_2d_from_0d_real(i_field, i1, i2, j1, j2,                  &
                                     field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 2d real field from a 0D field and updates the
!   relevant pointer in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1       ! Lower bound of field dimension 1
INTEGER, INTENT(IN) :: i2       ! Upper bound of field dimension 1
INTEGER, INTENT(IN) :: j1       ! Lower bound of field dimension 2
INTEGER, INTENT(IN) :: j2       ! Upper bound of field dimension 2
REAL, INTENT(IN) :: field_data  ! Field data supplied
REAL, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:) ! Environment field

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2))
  env_field = field_data
  environ_field_ptrs(i_field)%value_2d_real => env_field
END IF

RETURN
END SUBROUTINE set_field_2d_from_0d_real


! ----------------------------------------------------------------------
SUBROUTINE set_field_2d_from_0d_real_k(i_field, i1, i2, j1, j2, k1, k2,        &
                                       k, field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 2d real field from a 0D field and updates the
!   relevant pointer in the environment fields pointers array.
!   This routine handles a special case when the 2D environment field
!   to be set is stored internally with related 2D fields in a 3D array.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1       ! Lower bound of field dimension 1
INTEGER, INTENT(IN) :: i2       ! Upper bound of field dimension 1
INTEGER, INTENT(IN) :: j1       ! Lower bound of field dimension 2
INTEGER, INTENT(IN) :: j2       ! Upper bound of field dimension 2
INTEGER, INTENT(IN) :: k1       ! Lower bound in dimension 3
INTEGER, INTENT(IN) :: k2       ! Upper bound in dimension 3
INTEGER, INTENT(IN) :: k        ! Index of field in dimension 3
REAL, INTENT(IN) :: field_data  ! Field data supplied
REAL, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:,:)
                                ! Environment field

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2,k1:k2))
  env_field(:,:,k) = field_data
  environ_field_ptrs(i_field)%value_2d_real => env_field(:,:,k)
END IF

RETURN
END SUBROUTINE set_field_2d_from_0d_real_k


! ----------------------------------------------------------------------
SUBROUTINE set_field_2d_from_0d_integer(i_field, i1, i2, j1, j2,               &
                                        field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 2d integer field from a 0D field and updates the
!   relevant pointer in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field     ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1          ! Lower bound of field dimension 1
INTEGER, INTENT(IN) :: i2          ! Upper bound of field dimension 1
INTEGER, INTENT(IN) :: j1          ! Lower bound of field dimension 2
INTEGER, INTENT(IN) :: j2          ! Upper bound of field dimension 2
INTEGER, INTENT(IN) :: field_data  ! Field data supplied
INTEGER, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:)
                                   ! Environment field

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2))
  env_field = field_data
  environ_field_ptrs(i_field)%value_2d_integer => env_field
END IF

RETURN
END SUBROUTINE set_field_2d_from_0d_integer


! ----------------------------------------------------------------------
SUBROUTINE set_field_2d_from_0d_logical(i_field, i1, i2, j1, j2,               &
                                        field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 2d logical field from a 0D field and updates the
!   relevant pointer in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field     ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1          ! Lower bound of field dimension 1
INTEGER, INTENT(IN) :: i2          ! Upper bound of field dimension 1
INTEGER, INTENT(IN) :: j1          ! Lower bound of field dimension 2
INTEGER, INTENT(IN) :: j2          ! Upper bound of field dimension 2
LOGICAL, INTENT(IN) :: field_data  ! Field data supplied
LOGICAL, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:)
                                   ! Environment field

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2))
  env_field = field_data
  environ_field_ptrs(i_field)%value_2d_logical => env_field
END IF

RETURN
END SUBROUTINE set_field_2d_from_0d_logical


! ----------------------------------------------------------------------
SUBROUTINE set_field_3d_from_1d_real(i_field, i1, i2, j1, j2, k1, k2,          &
                                     field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 3d real field from a 1D field and updates the
!   relevant pointer in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1       ! Lower bound of field dimension 1
INTEGER, INTENT(IN) :: i2       ! Upper bound of field dimension 1
INTEGER, INTENT(IN) :: j1       ! Lower bound of field dimension 2
INTEGER, INTENT(IN) :: j2       ! Upper bound of field dimension 2
INTEGER, INTENT(IN) :: k1       ! Lower bound of field dimension 3
INTEGER, INTENT(IN) :: k2       ! Upper bound of field dimension 3
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:) ! Field data supplied
REAL, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:,:)
                                               ! Environment field

! Local variables
INTEGER :: i
INTEGER :: j

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2,k1:k2))
  DO i = i1,i2
    DO j = j1,j2
      env_field(i,j,:) = field_data(k1:k2)
    END DO
  END DO
  environ_field_ptrs(i_field)%value_3d_real => env_field
END IF

RETURN
END SUBROUTINE set_field_3d_from_1d_real


! ----------------------------------------------------------------------
SUBROUTINE set_field_4d_from_2d_real(i_field, i1, i2, j1, j2, k1, k2,          &
                                     n1, n2, field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 4d real field from a 2D field and updates the
!   relevant pointer in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1       ! Lower bound of field dimension 1
INTEGER, INTENT(IN) :: i2       ! Upper bound of field dimension 1
INTEGER, INTENT(IN) :: j1       ! Lower bound of field dimension 2
INTEGER, INTENT(IN) :: j2       ! Upper bound of field dimension 2
INTEGER, INTENT(IN) :: k1       ! Lower bound of field dimension 3
INTEGER, INTENT(IN) :: k2       ! Upper bound of field dimension 3
INTEGER, INTENT(IN) :: n1       ! Lower bound of field dimension 4
INTEGER, INTENT(IN) :: n2       ! Upper bound of field dimension 4
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:,:) ! Field data supplied
REAL, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:,:,:)
                                                 ! Environment field

! Local variables
INTEGER :: i
INTEGER :: j
INTEGER :: k

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2,k1:k2,n1:n2))
  DO i = i1,i2
    DO j = j1,j2
      DO k = k1,k2
        env_field(i,j,k,:) = field_data(k,n1:n2)
      END DO
    END DO
  END DO
  environ_field_ptrs(i_field)%value_4d_real => env_field
END IF

RETURN
END SUBROUTINE set_field_4d_from_2d_real

END MODULE ukca_environment_rdim_mod

