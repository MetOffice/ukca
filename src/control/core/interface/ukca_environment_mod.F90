! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module for handling UKCA's environmental driver fields. These are input
!   fields on the current UKCA model grid that may be varied by the parent
!   application during the run.
!
!   The module provides the following procedure for the UKCA API.
!
!     ukca_set_environment      - Sets or updates a named environment
!                                 field (overloaded for different field
!                                 dimensions and types).
!
!   Usage note: The land sea mask must be set before setting any environment
!   fields that are defined on land points only.
!
!   The following additional public procedure is provided for use within UKCA.
!
!     clear_environment_fields  - Clears data for all environment fields
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

MODULE ukca_environment_mod

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

USE ukca_fieldname_mod,  ONLY:                                                 &
  maxlen_fieldname,                                                            &
  fldname_sin_declination,                                                     &
  fldname_equation_of_time,                                                    &
  fldname_atmospheric_ch4,                                                     &
  fldname_atmospheric_co2,                                                     &
  fldname_atmospheric_h2,                                                      &
  fldname_atmospheric_n2,                                                      &
  fldname_atmospheric_o2,                                                      &
  fldname_atmospheric_n2o,                                                     &
  fldname_atmospheric_cfc11,                                                   &
  fldname_atmospheric_cfc12,                                                   &
  fldname_atmospheric_cfc113,                                                  &
  fldname_atmospheric_hcfc22,                                                  &
  fldname_atmospheric_hfc125,                                                  &
  fldname_atmospheric_hfc134a,                                                 &
  fldname_atmospheric_mebr,                                                    &
  fldname_atmospheric_mecl,                                                    &
  fldname_atmospheric_ch2br2,                                                  &
  fldname_atmospheric_chbr3,                                                   &
  fldname_atmospheric_cfc114,                                                  &
  fldname_atmospheric_cfc115,                                                  &
  fldname_atmospheric_ccl4,                                                    &
  fldname_atmospheric_meccl3,                                                  &
  fldname_atmospheric_hcfc141b,                                                &
  fldname_atmospheric_hcfc142b,                                                &
  fldname_atmospheric_h1211,                                                   &
  fldname_atmospheric_h1202,                                                   &
  fldname_atmospheric_h1301,                                                   &
  fldname_atmospheric_h2402,                                                   &
  fldname_atmospheric_cos,                                                     &
  fldname_soil_moisture_layer1,                                                &
  fldname_fland,                                                               &
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
  fldname_frac_types,                                                          &
  fldname_laift_lp,                                                            &
  fldname_canhtft_lp,                                                          &
  fldname_tstar_tile,                                                          &
  fldname_z0tile_lp,                                                           &
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
  fldname_l_tile_active,                                                       &
  fldname_u_rho_levels,                                                        &
  fldname_v_rho_levels,                                                        &
  fldname_geopH_on_theta_mlevs,                                                &
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
  fldname_photol_rates,                                                        &
  fldname_ibvoc_isoprene,                                                      &
  fldname_ibvoc_terpene,                                                       &
  fldname_ibvoc_methanol,                                                      &
  fldname_ibvoc_acetone,                                                       &
  fldname_inferno_bc,                                                          &
  fldname_inferno_ch4,                                                         &
  fldname_inferno_co,                                                          &
  fldname_inferno_nox,                                                         &
  fldname_inferno_oc,                                                          &
  fldname_inferno_so2,                                                         &
  fldname_lscat_zhang,                                                         &
  fldname_grid_area_fullht,                                                    &
  fldname_grid_volume

USE ukca_environment_fields_mod, ONLY:                                         &
  environ_field_info,                                                          &
  environ_field_ptrs,                                                          &
  l_environ_field_available,                                                   &
  no_bound_value,                                                              &
  k1_dust_flux, k2_dust_flux,                                                  &
  locate_land_points,                                                          &
  clear_land_only_fields,                                                      &
  land_points, land_index,                                                     &
  no_data_value,                                                               &
  sin_declination,                                                             &
  equation_of_time,                                                            &
  atmospheric_ch4,                                                             &
  atmospheric_co2,                                                             &
  atmospheric_h2,                                                              &
  atmospheric_n2,                                                              &
  atmospheric_o2,                                                              &
  atmospheric_n2o,                                                             &
  atmospheric_cfc11,                                                           &
  atmospheric_cfc12,                                                           &
  atmospheric_cfc113,                                                          &
  atmospheric_hcfc22,                                                          &
  atmospheric_hfc125,                                                          &
  atmospheric_hfc134a,                                                         &
  atmospheric_mebr,                                                            &
  atmospheric_mecl,                                                            &
  atmospheric_ch2br2,                                                          &
  atmospheric_chbr3,                                                           &
  atmospheric_cfc114,                                                          &
  atmospheric_cfc115,                                                          &
  atmospheric_ccl4,                                                            &
  atmospheric_meccl3,                                                          &
  atmospheric_hcfc141b,                                                        &
  atmospheric_hcfc142b,                                                        &
  atmospheric_h1211,                                                           &
  atmospheric_h1202,                                                           &
  atmospheric_h1301,                                                           &
  atmospheric_h2402,                                                           &
  atmospheric_cos,                                                             &
  soil_moisture_layer1,                                                        &
  fland,                                                                       &
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
  frac_types,                                                                  &
  laift_lp,                                                                    &
  canhtft_lp,                                                                  &
  tstar_tile,                                                                  &
  z0tile_lp,                                                                   &
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
  l_tile_active,                                                               &
  u_rho_levels,                                                                &
  v_rho_levels,                                                                &
  geopH_on_theta_mlevs,                                                        &
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
  photol_rates,                                                                &
  ibvoc_isoprene,                                                              &
  ibvoc_terpene,                                                               &
  ibvoc_methanol,                                                              &
  ibvoc_acetone,                                                               &
  inferno_bc,                                                                  &
  inferno_ch4,                                                                 &
  inferno_co,                                                                  &
  inferno_nox,                                                                 &
  inferno_oc,                                                                  &
  inferno_so2,                                                                 &
  lscat_zhang,                                                                 &
  grid_area_fullht,                                                            &
  grid_volume

USE ukca_environment_req_mod, ONLY: environ_field_index,                       &
                                    l_environ_req_available

USE ukca_environment_rdim_mod, ONLY: set_env_2d_from_0d_real,                  &
                                     set_env_2d_from_0d_integer,               &
                                     set_env_2d_from_0d_logical,               &
                                     set_env_3d_from_1d_real,                  &
                                     set_env_4d_from_2d_real

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_env_req_uninit, errcode_env_field_unknown,   &
                          errcode_env_field_mismatch

IMPLICIT NONE

PRIVATE

! Public procedures
PUBLIC ukca_set_environment, clear_environment_fields

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_ENVIRONMENT_MOD'

! Generic interface for subroutines to set or update a named environmental
! input field - overloaded according to the dimension and type of the field
! data.
INTERFACE ukca_set_environment
  MODULE PROCEDURE ukca_set_environment_0d_real
  MODULE PROCEDURE ukca_set_environment_0d_integer
  MODULE PROCEDURE ukca_set_environment_0d_logical
  MODULE PROCEDURE ukca_set_environment_1d_real
  MODULE PROCEDURE ukca_set_environment_2d_real
  MODULE PROCEDURE ukca_set_environment_2d_integer
  MODULE PROCEDURE ukca_set_environment_2d_logical
  MODULE PROCEDURE ukca_set_environment_3d_real
  MODULE PROCEDURE ukca_set_environment_4d_real
END INTERFACE ukca_set_environment


CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE ukca_set_environment_0d_real(varname, field_data, error_code,       &
                                        error_message, error_routine,          &
                                        field_index)
! ----------------------------------------------------------------------
! Description:
!   Variant of UKCA API procedure ukca_set_environment.
!   Sets or updates a named scalar environment field of type real.
!
! Method (for all variants of ukca_set_environment):
!
!   Look up the input variable name in the list of required fields and
!   use the associated field info to allocate the appropriate UKCA
!   internal array and copy in the field data.
!   Update the field availability array accordingly.
!   An optional field index argument will be set to the index in the
!   required field array or zero if the field is not included.
!
!   The 'field_data' argument must be a single value of the expected type
!   or an allocatable array having an appropriate dimension and the expected
!   type. Note that it is possible for fields with spatial dimensions
!   to be set from scalar input fields or fields with fewer dimensions
!   e.g. 3D fields can be set from 1D fields where there is no horizontal
!   variation.
!
!   Use of allocatable arrays preserves the array bounds when the
!   lower bound is not 1. The array bounds of the field data supplied are
!   allowed to extend beyond the UKCA grid and the fields may include halos
!   or other regions not relevant to UKCA. Any data outside the required
!   bounds are discarded.
!
!   A non-zero error code is returned if the requirement for the current
!   UKCA configuration has not been initialised, if the field name is not
!   recognised as a potential UKCA environment field of the dimension and
!   type supplied or if the field data supplied does not span the required
!   grid points.
!   In addition, a non-zero code is returned if an attempt is made to set a
!   field for which the required grid is not fully defined (occurring for
!   land-only fields in the absence of a land sea mask).
!   If an attempt is made to set a recognised field that is not required,
!   no action is taken. This occurence is indicated by both the error code
!   and field index being zero on return.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
REAL, INTENT(IN) :: field_data
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
INTEGER, OPTIONAL, INTENT(OUT) :: field_index

! Local variables

INTEGER :: i_field            ! Index of field in required fields array

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SET_ENVIRONMENT_0D_REAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of environment fields requirement
IF (.NOT. l_environ_req_available) THEN
  error_code = errcode_env_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Environment fields requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Find field index in array of required fields
i_field = environ_field_index(varname)
IF (PRESENT(field_index)) field_index = i_field

! Copy real field value to the appropriate UKCA internal variable if required
SELECT CASE (varname)
CASE (fldname_sin_declination)
  CALL set_field_0d_real(i_field, field_data, sin_declination)
CASE (fldname_equation_of_time)
  CALL set_field_0d_real(i_field, field_data, equation_of_time)
CASE (fldname_atmospheric_ch4)
  CALL set_field_0d_real(i_field, field_data, atmospheric_ch4)
CASE (fldname_atmospheric_co2)
  CALL set_field_0d_real(i_field, field_data, atmospheric_co2)
CASE (fldname_atmospheric_h2)
  CALL set_field_0d_real(i_field, field_data, atmospheric_h2)
CASE (fldname_atmospheric_n2)
  CALL set_field_0d_real(i_field, field_data, atmospheric_n2)
CASE (fldname_atmospheric_o2)
  CALL set_field_0d_real(i_field, field_data, atmospheric_o2)
CASE (fldname_atmospheric_n2o)
  CALL set_field_0d_real(i_field, field_data, atmospheric_n2o)
CASE (fldname_atmospheric_cfc11)
  CALL set_field_0d_real(i_field, field_data, atmospheric_cfc11)
CASE (fldname_atmospheric_cfc12)
  CALL set_field_0d_real(i_field, field_data, atmospheric_cfc12)
CASE (fldname_atmospheric_cfc113)
  CALL set_field_0d_real(i_field, field_data, atmospheric_cfc113)
CASE (fldname_atmospheric_hcfc22)
  CALL set_field_0d_real(i_field, field_data, atmospheric_hcfc22)
CASE (fldname_atmospheric_hfc125)
  CALL set_field_0d_real(i_field, field_data, atmospheric_hfc125)
CASE (fldname_atmospheric_hfc134a)
  CALL set_field_0d_real(i_field, field_data, atmospheric_hfc134a)
CASE (fldname_atmospheric_mebr)
  CALL set_field_0d_real(i_field, field_data, atmospheric_mebr)
CASE (fldname_atmospheric_mecl)
  CALL set_field_0d_real(i_field, field_data, atmospheric_mecl)
CASE (fldname_atmospheric_ch2br2)
  CALL set_field_0d_real(i_field, field_data, atmospheric_ch2br2)
CASE (fldname_atmospheric_chbr3)
  CALL set_field_0d_real(i_field, field_data, atmospheric_chbr3)
CASE (fldname_atmospheric_cfc114)
  CALL set_field_0d_real(i_field, field_data, atmospheric_cfc114)
CASE (fldname_atmospheric_cfc115)
  CALL set_field_0d_real(i_field, field_data, atmospheric_cfc115)
CASE (fldname_atmospheric_ccl4)
  CALL set_field_0d_real(i_field, field_data, atmospheric_ccl4)
CASE (fldname_atmospheric_meccl3)
  CALL set_field_0d_real(i_field, field_data, atmospheric_meccl3)
CASE (fldname_atmospheric_hcfc141b)
  CALL set_field_0d_real(i_field, field_data, atmospheric_hcfc141b)
CASE (fldname_atmospheric_hcfc142b)
  CALL set_field_0d_real(i_field, field_data, atmospheric_hcfc142b)
CASE (fldname_atmospheric_h1211)
  CALL set_field_0d_real(i_field, field_data, atmospheric_h1211)
CASE (fldname_atmospheric_h1202)
  CALL set_field_0d_real(i_field, field_data, atmospheric_h1202)
CASE (fldname_atmospheric_h1301)
  CALL set_field_0d_real(i_field, field_data, atmospheric_h1301)
CASE (fldname_atmospheric_h2402)
  CALL set_field_0d_real(i_field, field_data, atmospheric_h2402)
CASE (fldname_atmospheric_cos)
  CALL set_field_0d_real(i_field, field_data, atmospheric_cos)
CASE DEFAULT
  ! The named field may be a 2D field internally
  CALL set_env_2d_from_0d_real(varname, i_field, field_data, error_code)
END SELECT

IF (error_code == errcode_env_field_unknown) THEN
  IF (PRESENT(error_message)) error_message =                                  &
    'Unknown name for 0D real environmental input field: ''' //                &
    TRIM(varname) // ''''
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Update status to show that the field is available
IF (i_field /= 0) l_environ_field_available(i_field) = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_environment_0d_real


! ----------------------------------------------------------------------
SUBROUTINE ukca_set_environment_0d_integer(varname, field_data, error_code,    &
                                           error_message, error_routine,       &
                                           field_index)
! ----------------------------------------------------------------------
! Description:
!   Variant of UKCA API procedure ukca_set_environment.
!   Sets or updates a named scalar environment field of type integer.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: field_data
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
INTEGER, OPTIONAL, INTENT(OUT) :: field_index

! Local variables

INTEGER :: i_field            ! Index of field in required fields array

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SET_ENVIRONMENT_0D_INTEGER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of environment fields requirement
IF (.NOT. l_environ_req_available) THEN
  error_code = errcode_env_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Environment fields requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Find field index in array of required fields
i_field = environ_field_index(varname)
IF (PRESENT(field_index)) field_index = i_field

! Copy integer field value to the appropriate UKCA internal variable
! if required.
! The named field may be a 2D field internally.
CALL set_env_2d_from_0d_integer(varname, i_field, field_data, error_code)

IF (error_code == errcode_env_field_unknown) THEN
  IF (PRESENT(error_message)) error_message =                                  &
    'Unknown name for 0D integer environmental input field: ''' //             &
    TRIM(varname) // ''''
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Update status to show that the field is available
IF (i_field /= 0) l_environ_field_available(i_field) = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_environment_0d_integer


! ----------------------------------------------------------------------
SUBROUTINE ukca_set_environment_0d_logical(varname, field_data, error_code,    &
                                           error_message, error_routine,       &
                                           field_index)
! ----------------------------------------------------------------------
! Description:
!   Variant of UKCA API procedure ukca_set_environment.
!   Sets or updates a named scalar environment field of type logical.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, INTENT(IN) :: field_data
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
INTEGER, OPTIONAL, INTENT(OUT) :: field_index

! Local variables

INTEGER :: i_field            ! Index of field in required fields array

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SET_ENVIRONMENT_0D_LOGICAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of environment fields requirement
IF (.NOT. l_environ_req_available) THEN
  error_code = errcode_env_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Environment fields requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Find field index in array of required fields
i_field = environ_field_index(varname)
IF (PRESENT(field_index)) field_index = i_field

! Copy logical field value to the appropriate UKCA internal variable
! if required.
! The named field may be a 2D field internally
CALL set_env_2d_from_0d_logical(varname, i_field, field_data, error_code)

IF (error_code == errcode_env_field_unknown) THEN
  IF (PRESENT(error_message)) error_message =                                  &
    'Unknown name for 0D logical environmental input field: ''' //             &
    TRIM(varname) // ''''
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Update status to show that the field is available
IF (i_field /= 0) l_environ_field_available(i_field) = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_environment_0d_logical


! ----------------------------------------------------------------------
SUBROUTINE ukca_set_environment_1d_real(varname, field_data, error_code,       &
                                        error_message, error_routine,          &
                                        field_index)
! ----------------------------------------------------------------------
! Description:
!   Variant of UKCA API procedure ukca_set_environment.
!   Sets or updates a named 1D environment field of type real.
!
! Method:
!   See ukca_set_environment_0d_real.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:)
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
INTEGER, OPTIONAL, INTENT(OUT) :: field_index

! Local variables

! List of 3-D fields that can be populated from 1D real data arrays
! Have to use DATA statement as names are not of uniform length
CHARACTER(LEN=maxlen_fieldname), DIMENSION(58) :: fldnames_3d_real

DATA fldnames_3d_real / fldname_stcon, fldname_theta, fldname_q, fldname_qcf,  &
  fldname_conv_cloud_amount, fldname_rho_r2, fldname_qcl,                      &
  fldname_exner_rho_levels, fldname_area_cloud_fraction, fldname_cloud_frac,   &
  fldname_cloud_liq_frac, fldname_exner_theta_levels, fldname_p_rho_levels,    &
  fldname_p_theta_levels, fldname_rhokh_rdz, fldname_dtrdz, fldname_we_lim,    &
  fldname_t_frac, fldname_zrzi, fldname_we_lim_dsc, fldname_t_frac_dsc,        &
  fldname_zrzi_dsc, fldname_ls_rain3d, fldname_ls_snow3d, fldname_autoconv,    &
  fldname_accretion, fldname_pv_on_theta_mlevs, fldname_conv_rain3d,           &
  fldname_conv_snow3d, fldname_so4_sa_clim, fldname_so4_aitken,                &
  fldname_so4_accum, fldname_soot_fresh, fldname_soot_aged,                    &
  fldname_ocff_fresh, fldname_ocff_aged, fldname_biogenic, fldname_dust_div1,  &
  fldname_dust_div2, fldname_dust_div3, fldname_dust_div4, fldname_dust_div5,  &
  fldname_dust_div6, fldname_sea_salt_film, fldname_sea_salt_jet,              &
  fldname_co2_interactive, fldname_rim_cry, fldname_rim_agg, fldname_vertvel,  &
  fldname_bl_tke, fldname_interf_z, fldname_h2o2_offline,                      &
  fldname_ho2_offline, fldname_no3_offline, fldname_o3_offline,                &
  fldname_oh_offline, fldname_grid_area_fullht, fldname_grid_volume /

INTEGER :: i_field  ! Index of field in required fields array

! Required bounds of environment field data
INTEGER :: i1
INTEGER :: i2
INTEGER :: k1
INTEGER :: k2

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SET_ENVIRONMENT_1D_REAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of environment fields requirement
IF (.NOT. l_environ_req_available) THEN
  error_code = errcode_env_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Environment fields requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Find field index in array of required fields
i_field = environ_field_index(varname)
IF (PRESENT(field_index)) field_index = i_field

! If field is required, check the supplied field data array is allocated
IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(field_data)) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '1D real environment field for ''' // TRIM(varname) // ''' is unallocated'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Check whether this might be a 1-D field to populate a 3-D field
IF ( ANY(fldnames_3d_real == varname) ) THEN

  ! Field name is recognised, process it if it is a required field, else ignore
  IF (i_field /= 0) THEN
    ! Verify that the input array is not smaller than the last dimension
    k1 = environ_field_info(i_field)%lbound_dim3
    k2 = environ_field_info(i_field)%ubound_dim3
    IF (LBOUND(field_data,DIM=1) > k1 .OR. UBOUND(field_data,DIM=1) < k2) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message)) error_message =                              &
        'input 1D real environment field for 3-D field ''' // TRIM(varname) // &
        ''' has one or more invalid array bounds'
      IF (PRESENT(error_routine)) error_routine = RoutineName
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                    &
                              zhook_out, zhook_handle)
      RETURN
    END IF
    ! Set the field
    CALL set_env_3d_from_1d_real(varname, i_field, field_data, error_code)
    IF ( error_code /= 0 ) THEN
      IF (PRESENT(error_message)) error_message =                              &
        'Error populating 3-D field from 1-D real : ''' // TRIM(varname) // ''''
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE
      l_environ_field_available(i_field) = .TRUE.
    END IF
  END IF

  ! Return in any case with/ without errors
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN

END IF

! If field is required, check that its bounds are compatible with the UKCA
! configuration. The field data supplied must fill the required domain but may
! extend beyond it to avoid the need for pre-trimming by the parent model.
! Allow for the possibility that the expected bounds for land only fields may
! undefined if the UKCA environment field 'land_sea_mask' has not been set.
IF (i_field /= 0) THEN
  i1 = environ_field_info(i_field)%lbound_dim1
  i2 = environ_field_info(i_field)%ubound_dim1
  IF (i2 == no_bound_value) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      'The required dimension 1 upper bound is undefined for ' //              &
      '1D real environment field ''' // TRIM(varname) // ''''
  ELSE IF (LBOUND(field_data,DIM=1) > i1 .OR.                                  &
           UBOUND(field_data,DIM=1) < i2 .OR.                                  &
           (environ_field_info(i_field)%l_land_only .AND.                      &
            (LBOUND(field_data,DIM=1) /= i1 .OR.                               &
             UBOUND(field_data,DIM=1) /= i2))) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '1D real environment field for ''' // TRIM(varname) //                   &
      ''' has one or more invalid array bounds'
  END IF
  IF (error_code /= 0) THEN
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Copy 1D real field to the appropriate UKCA internal array if required.
! Any data outside the required bounds are discarded.
SELECT CASE (varname)
CASE (fldname_soil_moisture_layer1)
  CALL set_field_1d_real(i_field, i1, i2, field_data, soil_moisture_layer1)
CASE (fldname_fland)
  CALL set_field_1d_real(i_field, i1, i2, field_data, fland)
CASE (fldname_ibvoc_isoprene)
  CALL set_field_1d_real(i_field, i1, i2, field_data, ibvoc_isoprene)
CASE (fldname_ibvoc_terpene)
  CALL set_field_1d_real(i_field, i1, i2, field_data, ibvoc_terpene)
CASE (fldname_ibvoc_methanol)
  CALL set_field_1d_real(i_field, i1, i2, field_data, ibvoc_methanol)
CASE (fldname_ibvoc_acetone)
  CALL set_field_1d_real(i_field, i1, i2, field_data, ibvoc_acetone)
CASE (fldname_inferno_bc)
  CALL set_field_1d_real(i_field, i1, i2, field_data, inferno_bc)
CASE (fldname_inferno_ch4)
  CALL set_field_1d_real(i_field, i1, i2, field_data, inferno_ch4)
CASE (fldname_inferno_co)
  CALL set_field_1d_real(i_field, i1, i2, field_data, inferno_co)
CASE (fldname_inferno_nox)
  CALL set_field_1d_real(i_field, i1, i2, field_data, inferno_nox)
CASE (fldname_inferno_oc)
  CALL set_field_1d_real(i_field, i1, i2, field_data, inferno_oc)
CASE (fldname_inferno_so2)
  CALL set_field_1d_real(i_field, i1, i2, field_data, inferno_so2)
CASE DEFAULT
  ! Field not a recognised 1-D (or 3-D) field
  IF ( .NOT. (ANY(fldnames_3d_real == varname)) ) THEN
    IF (PRESENT(error_message)) error_message =                                &
      'Unknown name for 1D real environmental input field: ''' //              &
      TRIM(varname) // ''''
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END SELECT

! Update status to show that the field is available
IF (i_field /= 0) l_environ_field_available(i_field) = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_environment_1d_real


! ----------------------------------------------------------------------
SUBROUTINE ukca_set_environment_2d_real(varname, field_data, error_code,       &
                                        error_message, error_routine,          &
                                        field_index)
! ----------------------------------------------------------------------
! Description:
!   Variant of UKCA API procedure ukca_set_environment.
!   Sets or updates a named 2D environment field of type real.
!
! Method:
!   See ukca_set_environment_0d_real.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:,:)
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
INTEGER, OPTIONAL, INTENT(OUT) :: field_index

! Local variables

! List of 4-D fields that can be populated from 2D real arrays
CHARACTER(LEN=maxlen_fieldname), DIMENSION(1) :: fldnames_4d_real
DATA fldnames_4d_real / fldname_photol_rates /

INTEGER :: i_field  ! Index of field in required fields array

! Required bounds of environment field data
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
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SET_ENVIRONMENT_2D_REAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of environment fields requirement
IF (.NOT. l_environ_req_available) THEN
  error_code = errcode_env_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Environment fields requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Find field index in array of required fields
i_field = environ_field_index(varname)
IF (PRESENT(field_index)) field_index = i_field

! If field is required, check the supplied field data array is allocated
IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(field_data)) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '2D real environment field for ''' // TRIM(varname) // ''' is unallocated'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Check whether this might be a 2-D field to populate a 4-D field
IF ( ANY(fldnames_4d_real == varname) ) THEN

  ! Field name is recognised, process it if it is a required field, else ignore
  IF (i_field /= 0) THEN
    ! Verify that the spatial dimensions of the input array are not smaller than
    ! the corresponding dimensions of the target variable. The 4th dimension
    ! (usu. non-spatial, e.g. photol species) should match the bounds exactly.
    k1 = environ_field_info(i_field)%lbound_dim3
    k2 = environ_field_info(i_field)%ubound_dim3
    n1 = environ_field_info(i_field)%lbound_dim4
    n2 = environ_field_info(i_field)%ubound_dim4
    IF (LBOUND(field_data,DIM=1) > k1 .OR. UBOUND(field_data,DIM=1) < k2 .OR.  &
       LBOUND(field_data,DIM=2) /= n1 .OR. UBOUND(field_data,DIM=2) /= n2) THEN
      error_code = errcode_env_field_mismatch
      IF (PRESENT(error_message)) error_message =                              &
        'input 2D real environment field for 4-D field ''' // TRIM(varname) // &
        ''' has one or more invalid array bounds'
      IF (PRESENT(error_routine)) error_routine = RoutineName
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                    &
                              zhook_out, zhook_handle)
      RETURN
    END IF
    ! Set the 4-D field
    CALL set_env_4d_from_2d_real(varname, i_field, field_data, error_code)
    IF ( error_code /= 0 ) THEN
      IF (PRESENT(error_message)) error_message =                              &
       'Error populating 4-D field from 2-D real : ''' // TRIM(varname) // ''' '
      IF (PRESENT(error_routine)) error_routine = RoutineName
    ELSE
      l_environ_field_available(i_field) = .TRUE.
    END IF
  END IF

  ! Return in any case with/ without errors
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN

END IF

! If field is required, check that the input data bounds are compatible with
! the UKCA configuration. The field data supplied must fill the required domain
! but may extend beyond it to avoid the need for pre-trimming (e.g. halo
! removal) by the parent model.
! Allow for the possibility that the expected bounds for land only fields may
! undefined if the UKCA environment field 'land_sea_mask' has not been set.
IF (i_field /= 0) THEN
  i1 = environ_field_info(i_field)%lbound_dim1
  i2 = environ_field_info(i_field)%ubound_dim1
  j1 = environ_field_info(i_field)%lbound_dim2
  j2 = environ_field_info(i_field)%ubound_dim2
  IF (i2 == no_bound_value) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      'The required dimension 1 upper bound is undefined for ' //              &
      '2D real environment field ''' // TRIM(varname) // ''''
  ELSE IF (LBOUND(field_data,DIM=1) > i1 .OR.                                  &
           UBOUND(field_data,DIM=1) < i2 .OR.                                  &
           LBOUND(field_data,DIM=2) > j1 .OR.                                  &
           UBOUND(field_data,DIM=2) < j2 .OR.                                  &
           (environ_field_info(i_field)%l_land_only .AND.                      &
            (LBOUND(field_data,DIM=1) /= i1 .OR.                               &
             UBOUND(field_data,DIM=1) /= i2))) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '2D real environment field for ''' // TRIM(varname) //                   &
      ''' has one or more invalid array bounds'
  END IF
  IF (error_code /= 0) THEN
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Copy 2D real field to the appropriate UKCA internal array if required
! Any data outside the required bounds (e.g. halos) are discarded.
SELECT CASE (varname)
CASE (fldname_latitude)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, latitude)
CASE (fldname_longitude)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, longitude)
CASE (fldname_sin_latitude)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, sin_latitude)
CASE (fldname_cos_latitude)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, cos_latitude)
CASE (fldname_tan_latitude)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, tan_latitude)
CASE (fldname_conv_cloud_lwp)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, conv_cloud_lwp)
CASE (fldname_tstar)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, tstar)
CASE (fldname_zbl)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, zbl)
CASE (fldname_rough_length)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, rough_length)
CASE (fldname_seaice_frac)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, seaice_frac)
CASE (fldname_frac_types)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, frac_types)
CASE (fldname_laift_lp)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, laift_lp)
CASE (fldname_canhtft_lp)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, canhtft_lp)
CASE (fldname_tstar_tile)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, tstar_tile)
CASE (fldname_z0tile_lp)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, z0tile_lp)
CASE (fldname_pstar)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, pstar)
CASE (fldname_surf_albedo)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, surf_albedo)
CASE (fldname_zhsc)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, zhsc)
CASE (fldname_u_scalar_10m)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, u_scalar_10m)
CASE (fldname_surf_hf)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, surf_hf)
CASE (fldname_u_s)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, u_s)
CASE (fldname_ch4_wetl_emiss)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, ch4_wetl_emiss)
CASE (fldname_dms_sea_conc)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, dms_sea_conc)
CASE (fldname_chloro_sea)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, chloro_sea)
  ! Dust fluxes are stored in a 3D array so need special treatment
CASE (fldname_dust_flux_div1)
  CALL set_field_2d_real_k(i_field, i1, i2, j1, j2,                            &
                           k1_dust_flux, k2_dust_flux, 1, field_data, dust_flux)
CASE (fldname_dust_flux_div2)
  CALL set_field_2d_real_k(i_field, i1, i2, j1, j2,                            &
                           k1_dust_flux, k2_dust_flux, 2, field_data, dust_flux)
CASE (fldname_dust_flux_div3)
  CALL set_field_2d_real_k(i_field, i1, i2, j1, j2,                            &
                           k1_dust_flux, k2_dust_flux, 3, field_data, dust_flux)
CASE (fldname_dust_flux_div4)
  CALL set_field_2d_real_k(i_field, i1, i2, j1, j2,                            &
                           k1_dust_flux, k2_dust_flux, 4, field_data, dust_flux)
CASE (fldname_dust_flux_div5)
  CALL set_field_2d_real_k(i_field, i1, i2, j1, j2,                            &
                           k1_dust_flux, k2_dust_flux, 5, field_data, dust_flux)
CASE (fldname_dust_flux_div6)
  CALL set_field_2d_real_k(i_field, i1, i2, j1, j2,                            &
                           k1_dust_flux, k2_dust_flux, 6, field_data, dust_flux)
  ! End of dust flux cases
CASE (fldname_surf_wetness)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, surf_wetness)
CASE (fldname_grid_surf_area)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, grid_surf_area)
CASE (fldname_ext_cg_flash)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, ext_cg_flash)
CASE (fldname_ext_ic_flash)
  CALL set_field_2d_real(i_field, i1, i2, j1, j2, field_data, ext_ic_flash)
CASE DEFAULT
  ! Error: Not a recognised 2-D (or 4-D) field
  IF ( .NOT. (ANY(fldnames_4d_real == varname)) ) THEN
    error_code = errcode_env_field_unknown
    IF (PRESENT(error_message)) error_message =                                &
      'Unknown name for 2D real environmental input field: ''' //              &
      TRIM(varname) // ''''
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END SELECT

! Update status to show that the field is available
IF (i_field /= 0) l_environ_field_available(i_field) = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_environment_2d_real


! ----------------------------------------------------------------------
SUBROUTINE ukca_set_environment_2d_integer(varname, field_data, error_code,    &
                                           error_message, error_routine,       &
                                           field_index)
! ----------------------------------------------------------------------
! Description:
!   Variant of UKCA API procedure ukca_set_environment.
!   Sets or updates a named 2D environment field of type integer.
!
! Method:
!   See ukca_set_environment_0d_real.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
INTEGER, ALLOCATABLE, INTENT(IN) :: field_data(:,:)
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
INTEGER, OPTIONAL, INTENT(OUT) :: field_index

! Local variables

INTEGER :: i_field  ! Index of field in required fields array

! Required bounds of environment field data
INTEGER :: i1
INTEGER :: i2
INTEGER :: j1
INTEGER :: j2

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SET_ENVIRONMENT_2D_INTEGER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of environment fields requirement
IF (.NOT. l_environ_req_available) THEN
  error_code = errcode_env_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Environment fields requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Find field index in array of required fields
i_field = environ_field_index(varname)
IF (PRESENT(field_index)) field_index = i_field

! If field is required, check the supplied field data array is allocated and
! that its bounds are compatible with the UKCA configuration.
! The field data supplied must fill the required domain but may extend beyond
! it to avoid the need for pre-trimming (e.g. halo removal) by the parent model.
IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(field_data)) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '2D integer environment field for ''' // TRIM(varname) //                &
      ''' is unallocated'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
  i1 = environ_field_info(i_field)%lbound_dim1
  i2 = environ_field_info(i_field)%ubound_dim1
  j1 = environ_field_info(i_field)%lbound_dim2
  j2 = environ_field_info(i_field)%ubound_dim2
  IF (LBOUND(field_data,DIM=1) > i1 .OR. UBOUND(field_data,DIM=1) < i2 .OR.    &
      LBOUND(field_data,DIM=2) > j1 .OR. UBOUND(field_data,DIM=2) < j2) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '2D integer environment field for ''' // TRIM(varname) //                &
      ''' has one or more invalid array bounds'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Copy 2D integer field to the appropriate UKCA internal array if required.
! Any data outside the required bounds (e.g. halos) are discarded.
SELECT CASE (varname)
CASE (fldname_kent)
  CALL set_field_2d_integer(i_field, i1, i2, j1, j2, field_data, kent)
CASE (fldname_kent_dsc)
  CALL set_field_2d_integer(i_field, i1, i2, j1, j2, field_data, kent_dsc)
CASE (fldname_conv_cloud_base)
  CALL set_field_2d_integer(i_field, i1, i2, j1, j2, field_data,               &
                            conv_cloud_base)
CASE (fldname_conv_cloud_top)
  CALL set_field_2d_integer(i_field, i1, i2, j1, j2, field_data, conv_cloud_top)
CASE (fldname_lscat_zhang)
  CALL set_field_2d_integer(i_field, i1, i2, j1, j2, field_data, lscat_zhang)
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
  IF (PRESENT(error_message)) error_message =                                  &
    'Unknown name for 2D integer environmental input field: ''' //             &
    TRIM(varname) // ''''
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END SELECT

! Update status to show that the field is available
IF (i_field /= 0) l_environ_field_available(i_field) = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_environment_2d_integer


! ----------------------------------------------------------------------
SUBROUTINE ukca_set_environment_2d_logical(varname, field_data, error_code,    &
                                           error_message, error_routine,       &
                                           field_index)
! ----------------------------------------------------------------------
! Description:
!   Variant of UKCA API procedure ukca_set_environment.
!   Sets or updates a named 2D environment field of type logical
!
! Method:
!   See ukca_set_environment_0d_real.
!   If the field being set is the land sea mask, special processing is
!   called to set up an index of land points that is needed within UKCA
!   for locating land-only data on a 2D grid.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, ALLOCATABLE, INTENT(IN) :: field_data(:,:)
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
INTEGER, OPTIONAL, INTENT(OUT) :: field_index

! Local variables

INTEGER :: i_field  ! Index of field in required fields array

! Required bounds of environment field data
INTEGER :: i1
INTEGER :: i2
INTEGER :: j1
INTEGER :: j2

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SET_ENVIRONMENT_2D_LOGICAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of environment fields requirement
IF (.NOT. l_environ_req_available) THEN
  error_code = errcode_env_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Environment fields requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Find field index in array of required fields
i_field = environ_field_index(varname)
IF (PRESENT(field_index)) field_index = i_field

! If field is required, check the supplied field data array is allocated and
! that its bounds are compatible with the UKCA configuration.
! The field data supplied must fill the required domain but may extend beyond
! it to avoid the need for pre-trimming (e.g. halo removal) by the parent model.
IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(field_data)) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '2D logical environment field for ''' // TRIM(varname) //                &
      ''' is unallocated'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
  i1 = environ_field_info(i_field)%lbound_dim1
  i2 = environ_field_info(i_field)%ubound_dim1
  j1 = environ_field_info(i_field)%lbound_dim2
  j2 = environ_field_info(i_field)%ubound_dim2
  IF (LBOUND(field_data,DIM=1) > i1 .OR. UBOUND(field_data,DIM=1) < i2 .OR.    &
      LBOUND(field_data,DIM=2) > j1 .OR. UBOUND(field_data,DIM=2) < j2) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '2D logical environment field for ''' // TRIM(varname) //                &
      ''' has one or more invalid array bounds'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Copy 2D logical field to the appropriate UKCA internal array if required.
! Any data outside the required bounds (e.g. halos) are discarded.
SELECT CASE (varname)
CASE (fldname_land_sea_mask)
  IF (i_field /= 0) THEN
    CALL set_field_2d_logical(i_field, i1, i2, j1, j2, field_data,             &
                              land_sea_mask)
    ! Special processing to set up index of land points for locating
    ! land-only environment fields on a 2-D grid
    CALL locate_land_points()
    ! Clear any existing land-only fields as these fields may be inconsistent
    ! with the new land sea mask
    CALL clear_land_only_fields()
  END IF
CASE (fldname_l_tile_active)
  CALL set_field_2d_logical(i_field, i1, i2, j1, j2, field_data, l_tile_active)
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
  IF (PRESENT(error_message)) error_message =                                  &
    'Unknown name for 2D logical environmental input field: ''' //             &
    TRIM(varname) // ''''
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END SELECT

! Update status to show that the field is available
IF (i_field /= 0) l_environ_field_available(i_field) = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_environment_2d_logical


! ----------------------------------------------------------------------
SUBROUTINE ukca_set_environment_3d_real(varname, field_data, error_code,       &
                                        error_message, error_routine,          &
                                        field_index)
! ----------------------------------------------------------------------
! Description:
!   Variant of UKCA API procedure ukca_set_environment.
!   Sets or updates a named 3D environment field of type real.
!
! Method:
!   See ukca_set_environment_0d_real.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:,:,:)
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
INTEGER, OPTIONAL, INTENT(OUT) :: field_index

! Local variables

INTEGER :: i_field  ! Index of field in required fields array

! Required bounds of environment field data
INTEGER :: i1
INTEGER :: i2
INTEGER :: j1
INTEGER :: j2
INTEGER :: k1
INTEGER :: k2

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SET_ENVIRONMENT_3D_REAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of environment fields requirement
IF (.NOT. l_environ_req_available) THEN
  error_code = errcode_env_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Environment fields requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Find field index in array of required fields
i_field = environ_field_index(varname)
IF (PRESENT(field_index)) field_index = i_field

! If field is required, check the supplied field data array is allocated and
! that its bounds are compatible with the UKCA configuration.
! The field data supplied must fill the required domain but may extend beyond
! it to avoid the need for pre-trimming (e.g. halo removal) by the parent model.
IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(field_data)) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '3D real environment field for ''' // TRIM(varname) // ''' is unallocated'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
  i1 = environ_field_info(i_field)%lbound_dim1
  i2 = environ_field_info(i_field)%ubound_dim1
  j1 = environ_field_info(i_field)%lbound_dim2
  j2 = environ_field_info(i_field)%ubound_dim2
  k1 = environ_field_info(i_field)%lbound_dim3
  k2 = environ_field_info(i_field)%ubound_dim3
  IF (LBOUND(field_data,DIM=1) > i1 .OR. UBOUND(field_data,DIM=1) < i2 .OR.    &
      LBOUND(field_data,DIM=2) > j1 .OR. UBOUND(field_data,DIM=2) < j2 .OR.    &
      LBOUND(field_data,DIM=3) > k1 .OR. UBOUND(field_data,DIM=3) < k2) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '3D real environment field for ''' // TRIM(varname) //                   &
      ''' has one or more invalid array bounds'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Copy 3D real field to the appropriate UKCA internal array if required.
! Any data outside the required bounds (e.g. halos) are discarded.
SELECT CASE (varname)
CASE (fldname_u_rho_levels)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         u_rho_levels)
CASE (fldname_v_rho_levels)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         v_rho_levels)
CASE (fldname_geopH_on_theta_mlevs)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         geopH_on_theta_mlevs)
CASE (fldname_theta)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, theta)
CASE (fldname_q)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, q)
CASE (fldname_qcf)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, qcf)
CASE (fldname_conv_cloud_amount)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         conv_cloud_amount)
CASE (fldname_rho_r2)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, rho_r2)
CASE (fldname_qcl)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, qcl)
CASE (fldname_exner_rho_levels)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         exner_rho_levels)
CASE (fldname_area_cloud_fraction)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         area_cloud_fraction)
CASE (fldname_cloud_frac)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         cloud_frac)
CASE (fldname_cloud_liq_frac)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         cloud_liq_frac)
CASE (fldname_exner_theta_levels)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         exner_theta_levels)
CASE (fldname_p_rho_levels)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         p_rho_levels)
CASE (fldname_p_theta_levels)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         p_theta_levels)
CASE (fldname_rhokh_rdz)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, rhokh_rdz)
CASE (fldname_dtrdz)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, dtrdz)
CASE (fldname_we_lim)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, we_lim)
CASE (fldname_t_frac)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, t_frac)
CASE (fldname_zrzi)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, zrzi)
CASE (fldname_we_lim_dsc)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         we_lim_dsc)
CASE (fldname_t_frac_dsc)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         t_frac_dsc)
CASE (fldname_zrzi_dsc)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, zrzi_dsc)
CASE (fldname_stcon)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, stcon)
CASE (fldname_ls_rain3d)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, ls_rain3d)
CASE (fldname_ls_snow3d)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, ls_snow3d)
CASE (fldname_autoconv)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, autoconv)
CASE (fldname_accretion)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, accretion)
CASE (fldname_pv_on_theta_mlevs)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         pv_on_theta_mlevs)
CASE (fldname_conv_rain3d)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         conv_rain3d)
CASE (fldname_conv_snow3d)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         conv_snow3d)
CASE (fldname_so4_sa_clim)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         so4_sa_clim)
CASE (fldname_so4_aitken)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         so4_aitken)
CASE (fldname_so4_accum)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, so4_accum)
CASE (fldname_soot_fresh)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         soot_fresh)
CASE (fldname_soot_aged)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, soot_aged)
CASE (fldname_ocff_fresh)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         ocff_fresh)
CASE (fldname_ocff_aged)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, ocff_aged)
CASE (fldname_biogenic)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, biogenic)
CASE (fldname_dust_div1)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, dust_div1)
CASE (fldname_dust_div2)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, dust_div2)
CASE (fldname_dust_div3)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, dust_div3)
CASE (fldname_dust_div4)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, dust_div4)
CASE (fldname_dust_div5)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, dust_div5)
CASE (fldname_dust_div6)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, dust_div6)
CASE (fldname_sea_salt_film)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         sea_salt_film)
CASE (fldname_sea_salt_jet)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         sea_salt_jet)
CASE (fldname_co2_interactive)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         co2_interactive)
CASE (fldname_rim_cry)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, rim_cry)
CASE (fldname_rim_agg)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, rim_agg)
CASE (fldname_vertvel)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, vertvel)
CASE (fldname_bl_tke)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, bl_tke)
CASE (fldname_interf_z)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data, interf_z)
CASE (fldname_h2o2_offline)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         h2o2_offline)
CASE (fldname_ho2_offline)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         ho2_offline)
CASE (fldname_no3_offline)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         no3_offline)
CASE (fldname_o3_offline)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         o3_offline)
CASE (fldname_oh_offline)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         oh_offline)
CASE (fldname_grid_area_fullht)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         grid_area_fullht)
CASE (fldname_grid_volume)
  CALL set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2, field_data,          &
                         grid_volume)
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
  IF (PRESENT(error_message)) error_message =                                  &
    'Unknown name for 3D real environmental input field: ''' //                &
    TRIM(varname) // ''''
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END SELECT

! Update status to show that the field is available
IF (i_field /= 0) l_environ_field_available(i_field) = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_environment_3d_real


! ----------------------------------------------------------------------
SUBROUTINE ukca_set_environment_4d_real(varname, field_data, error_code,       &
                                        error_message, error_routine,          &
                                        field_index)
! ----------------------------------------------------------------------
! Description:
!   Variant of UKCA API procedure ukca_set_environment.
!   Sets or updates a named 4D environment field of type real.
!
! Method:
!   See ukca_set_environment_0d.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:,:,:,:)
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
INTEGER, OPTIONAL, INTENT(OUT) :: field_index

! Local variables

INTEGER :: i_field  ! Index of field in required fields array

! Required bounds of environment field data
INTEGER :: i1
INTEGER :: i2
INTEGER :: j1
INTEGER :: j2
INTEGER :: k1
INTEGER :: k2
INTEGER :: n1
INTEGER :: n2

! Dr hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_SET_ENVIRONMENT_4D_REAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of environment fields requirement
IF (.NOT. l_environ_req_available) THEN
  error_code = errcode_env_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Environment fields requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Find field index in array of required fields
i_field = environ_field_index(varname)
IF (PRESENT(field_index)) field_index = i_field

! If field is required, check the supplied field data array is allocated and
! that its bounds are compatible with the UKCA configuration.
! The field data supplied must fill the required domain but may extend beyond
! it to avoid the need for pre-trimming (e.g. halo removal) by the parent model.
! The bounds are expected to match in the last dimension.
IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(field_data)) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '4D real environment field for ''' // TRIM(varname) // ''' is unallocated'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
  i1 = environ_field_info(i_field)%lbound_dim1
  i2 = environ_field_info(i_field)%ubound_dim1
  j1 = environ_field_info(i_field)%lbound_dim2
  j2 = environ_field_info(i_field)%ubound_dim2
  k1 = environ_field_info(i_field)%lbound_dim3
  k2 = environ_field_info(i_field)%ubound_dim3
  n1 = environ_field_info(i_field)%lbound_dim4
  n2 = environ_field_info(i_field)%ubound_dim4
  IF (LBOUND(field_data,DIM=1) > i1 .OR. UBOUND(field_data,DIM=1) < i2 .OR.    &
      LBOUND(field_data,DIM=2) > j1 .OR. UBOUND(field_data,DIM=2) < j2 .OR.    &
      LBOUND(field_data,DIM=3) > k1 .OR. UBOUND(field_data,DIM=3) < k2 .OR.    &
      LBOUND(field_data,DIM=4) /= n1 .OR. UBOUND(field_data,DIM=4) /= n2) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '4D real environment field for ''' // TRIM(varname) //                   &
      ''' has one or more invalid array bounds'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Copy 4D real field to the appropriate UKCA internal array if required
! Any data outside the required bounds (e.g. halos) are discarded.
SELECT CASE (varname)
CASE (fldname_photol_rates)
  CALL set_field_4d_real(i_field, i1, i2, j1, j2, k1, k2, n1, n2, field_data,  &
                         photol_rates)
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
  IF (PRESENT(error_message)) error_message =                                  &
    'Unknown name for 4D real environmental input field: ''' //                &
    TRIM(varname) // ''''
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END SELECT

! Update status to show that the field is available
IF (i_field /= 0) l_environ_field_available(i_field) = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_environment_4d_real


! ----------------------------------------------------------------------
SUBROUTINE set_field_0d_real(i_field, field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 0D field and updates the relevant pointer
!   in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
REAL, INTENT(IN) :: field_data  ! Field data supplied
REAL, TARGET, INTENT(IN OUT) :: env_field  ! Environment field

IF (i_field /= 0) THEN
  env_field = field_data
  environ_field_ptrs(i_field)%value_0d_real => env_field
END IF

RETURN
END SUBROUTINE set_field_0d_real


! ----------------------------------------------------------------------
SUBROUTINE set_field_1d_real(i_field, i1, i2, field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 1D real field and updates the relevant pointer
!   in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1       ! Lower bound of field
INTEGER, INTENT(IN) :: i2       ! Upper bound of field
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:)  ! Field data supplied
REAL, TARGET, ALLOCATABLE, INTENT(IN OUT) :: env_field(:)  ! Environment field

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2))
  env_field = field_data(i1:i2)
  environ_field_ptrs(i_field)%value_1d_real => env_field
END IF

RETURN
END SUBROUTINE set_field_1d_real


! ----------------------------------------------------------------------
SUBROUTINE set_field_2d_real(i_field, i1, i2, j1, j2, field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 2d real field and updates the relevant pointer
!   in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1       ! Lower bound of field in dimension 1
INTEGER, INTENT(IN) :: i2       ! Upper bound of field in dimension 1
INTEGER, INTENT(IN) :: j1       ! Lower bound of field in dimension 2
INTEGER, INTENT(IN) :: j2       ! Upper bound of field in dimension 2
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:,:)  ! Field data supplied
REAL, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:)  ! Environment field

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2))
  env_field = field_data(i1:i2,j1:j2)
  environ_field_ptrs(i_field)%value_2d_real => env_field
END IF

RETURN
END SUBROUTINE set_field_2d_real


! ----------------------------------------------------------------------
SUBROUTINE set_field_2d_real_k(i_field, i1, i2, j1, j2, k1, k2, k,             &
                               field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 2d real field and updates the relevant pointer
!   in the environment fields pointers array.
!   This routine handles a special case when the 2D environment field
!   to be set is stored internally with related 2D fields in a 3D array.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1       ! Lower bound of field in dimension 1
INTEGER, INTENT(IN) :: i2       ! Upper bound of field in dimension 1
INTEGER, INTENT(IN) :: j1       ! Lower bound of field in dimension 2
INTEGER, INTENT(IN) :: j2       ! Upper bound of field in dimension 2
INTEGER, INTENT(IN) :: k1       ! Lower bound in dimension 3
INTEGER, INTENT(IN) :: k2       ! Upper bound in dimension 3
INTEGER, INTENT(IN) :: k        ! Index of field in dimension 3
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:,:)    ! Field data supplied
REAL, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:,:)
                                                    ! Environment field

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2,k1:k2))
  env_field(:,:,k) = field_data(i1:i2,j1:j2)
  environ_field_ptrs(i_field)%value_2d_real => env_field(:,:,k)
END IF

RETURN
END SUBROUTINE set_field_2d_real_k


! ----------------------------------------------------------------------
SUBROUTINE set_field_2d_integer(i_field, i1, i2, j1, j2, field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 2d integer field and updates the relevant pointer
!   in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1       ! Lower bound of field in dimension 1
INTEGER, INTENT(IN) :: i2       ! Upper bound of field in dimension 1
INTEGER, INTENT(IN) :: j1       ! Lower bound of field in dimension 2
INTEGER, INTENT(IN) :: j2       ! Upper bound of field in dimension 2
INTEGER, ALLOCATABLE, INTENT(IN) :: field_data(:,:)  ! Field data supplied
INTEGER, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:)
                                                     ! Environment field

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2))
  env_field = field_data(i1:i2,j1:j2)
  environ_field_ptrs(i_field)%value_2d_integer => env_field
END IF

RETURN
END SUBROUTINE set_field_2d_integer


! ----------------------------------------------------------------------
SUBROUTINE set_field_2d_logical(i_field, i1, i2, j1, j2, field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 2D logical field and updates the relevant pointer
!   in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1       ! Lower bound of field in dimension 1
INTEGER, INTENT(IN) :: i2       ! Upper bound of field in dimension 1
INTEGER, INTENT(IN) :: j1       ! Lower bound of field in dimension 2
INTEGER, INTENT(IN) :: j2       ! Upper bound of field in dimension 2
LOGICAL, ALLOCATABLE, INTENT(IN) :: field_data(:,:)  ! Field data supplied
LOGICAL, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:)
                                                     ! Environment field

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2))
  env_field = field_data(i1:i2,j1:j2)
  environ_field_ptrs(i_field)%value_2d_logical => env_field
END IF

RETURN
END SUBROUTINE set_field_2d_logical


! ----------------------------------------------------------------------
SUBROUTINE set_field_3d_real(i_field, i1, i2, j1, j2, k1, k2,                  &
                             field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 3d real field and updates the relevant pointer
!   in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1       ! Lower bound of field of dimension 1
INTEGER, INTENT(IN) :: i2       ! Upper bound of field of diemnsion 1
INTEGER, INTENT(IN) :: j1       ! Lower bound of field of dimension 2
INTEGER, INTENT(IN) :: j2       ! Upper bound of field of diemnsion 2
INTEGER, INTENT(IN) :: k1       ! Lower bound of field of dimension 3
INTEGER, INTENT(IN) :: k2       ! Upper bound of field of diemnsion 3
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:,:,:)  ! Field data supplied
REAL, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:,:)
                                                    ! Environment field

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2,k1:k2))
  env_field = field_data(i1:i2,j1:j2,k1:k2)
  environ_field_ptrs(i_field)%value_3d_real => env_field
END IF

RETURN
END SUBROUTINE set_field_3d_real


! ----------------------------------------------------------------------
SUBROUTINE set_field_4d_real(i_field, i1, i2, j1, j2, k1, k2, n1, n2,          &
                             field_data, env_field)
! ----------------------------------------------------------------------
! Description:
!   Sets the value of a 4d real field and updates the relevant pointer
!   in the environment fields pointers array
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_field  ! Index of field in required fields array
INTEGER, INTENT(IN) :: i1       ! Lower bound of field of dimension 1
INTEGER, INTENT(IN) :: i2       ! Upper bound of field of diemnsion 1
INTEGER, INTENT(IN) :: j1       ! Lower bound of field of dimension 2
INTEGER, INTENT(IN) :: j2       ! Upper bound of field of diemnsion 2
INTEGER, INTENT(IN) :: k1       ! Lower bound of field of dimension 3
INTEGER, INTENT(IN) :: k2       ! Upper bound of field of diemnsion 3
INTEGER, INTENT(IN) :: n1       ! Lower bound of field of dimension 4
INTEGER, INTENT(IN) :: n2       ! Upper bound of field of diemnsion 4
REAL, ALLOCATABLE, INTENT(IN) :: field_data(:,:,:,:) ! Field data supplied
REAL, ALLOCATABLE, TARGET, INTENT(IN OUT) :: env_field(:,:,:,:)
                                                     ! Environment field

IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(env_field)) ALLOCATE(env_field(i1:i2,j1:j2,k1:k2,n1:n2))
  env_field = field_data(i1:i2,j1:j2,k1:k2,n1:n2)
  environ_field_ptrs(i_field)%value_4d_real => env_field
END IF

RETURN
END SUBROUTINE set_field_4d_real


! ----------------------------------------------------------------------
SUBROUTINE clear_environment_fields()
! ----------------------------------------------------------------------
! Description:
!   Resets scalar fields values, deallocates environment field arrays
!   and land point index array and updates availability flags
! ----------------------------------------------------------------------

IMPLICIT NONE

! Local variables

INTEGER :: i

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CLEAR_ENVIRONMENT_FIELDS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Clear environment field pointers before re-setting/de-allocating the fields
DO i = 1, SIZE(environ_field_ptrs)
  NULLIFY(environ_field_ptrs(i)%value_0d_real)
  NULLIFY(environ_field_ptrs(i)%value_1d_real)
  NULLIFY(environ_field_ptrs(i)%value_2d_integer)
  NULLIFY(environ_field_ptrs(i)%value_2d_real)
  NULLIFY(environ_field_ptrs(i)%value_2d_logical)
  NULLIFY(environ_field_ptrs(i)%value_3d_real)
  NULLIFY(environ_field_ptrs(i)%value_4d_real)
END DO

! Clear the environment fields defined on land points only
CALL clear_land_only_fields()

! Reset/dellocate the remaining environment fields
sin_declination = no_data_value
equation_of_time = no_data_value
atmospheric_ch4 = no_data_value
atmospheric_co2 = no_data_value
atmospheric_h2 = no_data_value
atmospheric_n2 = no_data_value
atmospheric_o2 = no_data_value
atmospheric_n2o = no_data_value
atmospheric_cfc11 = no_data_value
atmospheric_cfc12 = no_data_value
atmospheric_cfc113 = no_data_value
atmospheric_hcfc22 = no_data_value
atmospheric_hfc125 = no_data_value
atmospheric_hfc134a = no_data_value
atmospheric_mebr = no_data_value
atmospheric_mecl = no_data_value
atmospheric_ch2br2 = no_data_value
atmospheric_chbr3 = no_data_value
atmospheric_cfc114 = no_data_value
atmospheric_cfc115 = no_data_value
atmospheric_ccl4 = no_data_value
atmospheric_meccl3 = no_data_value
atmospheric_hcfc141b = no_data_value
atmospheric_hcfc142b = no_data_value
atmospheric_h1211 = no_data_value
atmospheric_h1202 = no_data_value
atmospheric_h1301 = no_data_value
atmospheric_h2402 = no_data_value
atmospheric_cos = no_data_value
IF (ALLOCATED(latitude)) DEALLOCATE(latitude)
IF (ALLOCATED(longitude)) DEALLOCATE(longitude)
IF (ALLOCATED(sin_latitude)) DEALLOCATE(sin_latitude)
IF (ALLOCATED(cos_latitude)) DEALLOCATE(cos_latitude)
IF (ALLOCATED(tan_latitude)) DEALLOCATE(tan_latitude)
IF (ALLOCATED(conv_cloud_lwp)) DEALLOCATE(conv_cloud_lwp)
IF (ALLOCATED(tstar)) DEALLOCATE(tstar)
IF (ALLOCATED(zbl)) DEALLOCATE(zbl)
IF (ALLOCATED(rough_length)) DEALLOCATE(rough_length)
IF (ALLOCATED(seaice_frac)) DEALLOCATE(seaice_frac)
IF (ALLOCATED(pstar)) DEALLOCATE(pstar)
IF (ALLOCATED(surf_albedo)) DEALLOCATE(surf_albedo)
IF (ALLOCATED(zhsc)) DEALLOCATE(zhsc)
IF (ALLOCATED(u_scalar_10m)) DEALLOCATE(u_scalar_10m)
IF (ALLOCATED(surf_hf)) DEALLOCATE(surf_hf)
IF (ALLOCATED(u_s)) DEALLOCATE(u_s)
IF (ALLOCATED(ch4_wetl_emiss)) DEALLOCATE(ch4_wetl_emiss)
IF (ALLOCATED(dms_sea_conc)) DEALLOCATE(dms_sea_conc)
IF (ALLOCATED(chloro_sea)) DEALLOCATE(chloro_sea)
IF (ALLOCATED(dust_flux)) DEALLOCATE(dust_flux)
IF (ALLOCATED(surf_wetness)) DEALLOCATE(surf_wetness)
IF (ALLOCATED(kent)) DEALLOCATE(kent)
IF (ALLOCATED(kent_dsc)) DEALLOCATE(kent_dsc)
IF (ALLOCATED(conv_cloud_base)) DEALLOCATE(conv_cloud_base)
IF (ALLOCATED(conv_cloud_top)) DEALLOCATE(conv_cloud_top)
IF (ALLOCATED(land_sea_mask)) DEALLOCATE(land_sea_mask)
IF (ALLOCATED(u_rho_levels)) DEALLOCATE(u_rho_levels)
IF (ALLOCATED(v_rho_levels)) DEALLOCATE(v_rho_levels)
IF (ALLOCATED(geopH_on_theta_mlevs)) DEALLOCATE(geopH_on_theta_mlevs)
IF (ALLOCATED(theta)) DEALLOCATE(theta)
IF (ALLOCATED(q)) DEALLOCATE(q)
IF (ALLOCATED(qcf)) DEALLOCATE(qcf)
IF (ALLOCATED(conv_cloud_amount)) DEALLOCATE(conv_cloud_amount)
IF (ALLOCATED(rho_r2)) DEALLOCATE(rho_r2)
IF (ALLOCATED(qcl)) DEALLOCATE(qcl)
IF (ALLOCATED(exner_rho_levels)) DEALLOCATE(exner_rho_levels)
IF (ALLOCATED(area_cloud_fraction)) DEALLOCATE(area_cloud_fraction)
IF (ALLOCATED(cloud_frac)) DEALLOCATE(cloud_frac)
IF (ALLOCATED(cloud_liq_frac)) DEALLOCATE(cloud_liq_frac)
IF (ALLOCATED(exner_theta_levels)) DEALLOCATE(exner_theta_levels)
IF (ALLOCATED(p_rho_levels)) DEALLOCATE(p_rho_levels)
IF (ALLOCATED(p_theta_levels)) DEALLOCATE(p_theta_levels)
IF (ALLOCATED(rhokh_rdz)) DEALLOCATE(rhokh_rdz)
IF (ALLOCATED(dtrdz)) DEALLOCATE(dtrdz)
IF (ALLOCATED(we_lim)) DEALLOCATE(we_lim)
IF (ALLOCATED(t_frac)) DEALLOCATE(t_frac)
IF (ALLOCATED(zrzi)) DEALLOCATE(zrzi)
IF (ALLOCATED(we_lim_dsc)) DEALLOCATE(we_lim_dsc)
IF (ALLOCATED(t_frac_dsc)) DEALLOCATE(t_frac_dsc)
IF (ALLOCATED(zrzi_dsc)) DEALLOCATE(zrzi_dsc)
IF (ALLOCATED(stcon)) DEALLOCATE(stcon)
IF (ALLOCATED(ls_rain3d)) DEALLOCATE(ls_rain3d)
IF (ALLOCATED(ls_snow3d)) DEALLOCATE(ls_snow3d)
IF (ALLOCATED(autoconv)) DEALLOCATE(autoconv)
IF (ALLOCATED(accretion)) DEALLOCATE(accretion)
IF (ALLOCATED(pv_on_theta_mlevs)) DEALLOCATE(pv_on_theta_mlevs)
IF (ALLOCATED(conv_rain3d)) DEALLOCATE(conv_rain3d)
IF (ALLOCATED(conv_snow3d)) DEALLOCATE(conv_snow3d)
IF (ALLOCATED(so4_sa_clim)) DEALLOCATE(so4_sa_clim)
IF (ALLOCATED(so4_aitken)) DEALLOCATE(so4_aitken)
IF (ALLOCATED(so4_accum)) DEALLOCATE(so4_accum)
IF (ALLOCATED(soot_fresh)) DEALLOCATE(soot_fresh)
IF (ALLOCATED(soot_aged)) DEALLOCATE(soot_aged)
IF (ALLOCATED(ocff_fresh)) DEALLOCATE(ocff_fresh)
IF (ALLOCATED(ocff_aged)) DEALLOCATE(ocff_aged)
IF (ALLOCATED(biogenic)) DEALLOCATE(biogenic)
IF (ALLOCATED(dust_div1)) DEALLOCATE(dust_div1)
IF (ALLOCATED(dust_div2)) DEALLOCATE(dust_div2)
IF (ALLOCATED(dust_div3)) DEALLOCATE(dust_div3)
IF (ALLOCATED(dust_div4)) DEALLOCATE(dust_div4)
IF (ALLOCATED(dust_div5)) DEALLOCATE(dust_div5)
IF (ALLOCATED(dust_div6)) DEALLOCATE(dust_div6)
IF (ALLOCATED(sea_salt_film)) DEALLOCATE(sea_salt_film)
IF (ALLOCATED(sea_salt_jet)) DEALLOCATE(sea_salt_jet)
IF (ALLOCATED(co2_interactive)) DEALLOCATE(co2_interactive)
IF (ALLOCATED(rim_cry)) DEALLOCATE(rim_cry)
IF (ALLOCATED(rim_agg)) DEALLOCATE(rim_agg)
IF (ALLOCATED(vertvel)) DEALLOCATE(vertvel)
IF (ALLOCATED(bl_tke)) DEALLOCATE(bl_tke)
IF (ALLOCATED(h2o2_offline)) DEALLOCATE(h2o2_offline)
IF (ALLOCATED(ho2_offline)) DEALLOCATE(ho2_offline)
IF (ALLOCATED(no3_offline)) DEALLOCATE(no3_offline)
IF (ALLOCATED(o3_offline)) DEALLOCATE(o3_offline)
IF (ALLOCATED(oh_offline)) DEALLOCATE(oh_offline)
IF (ALLOCATED(interf_z)) DEALLOCATE(interf_z)
IF (ALLOCATED(grid_surf_area)) DEALLOCATE(grid_surf_area)
IF (ALLOCATED(ext_cg_flash)) DEALLOCATE(ext_cg_flash)
IF (ALLOCATED(ext_ic_flash)) DEALLOCATE(ext_ic_flash)
IF (ALLOCATED(photol_rates)) DEALLOCATE(photol_rates)
IF (ALLOCATED(lscat_zhang)) DEALLOCATE(lscat_zhang)
IF (ALLOCATED(grid_area_fullht)) DEALLOCATE(grid_area_fullht)
IF (ALLOCATED(grid_volume)) DEALLOCATE(grid_volume)

! Update field availability
l_environ_field_available(:) = .FALSE.

! The land sea mask is no longer valid so clear land point indices
IF (ALLOCATED(land_index)) DEALLOCATE(land_index)
land_points = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE clear_environment_fields

END MODULE ukca_environment_mod
