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
  fldname_lscat_zhang

USE ukca_environment_fields_mod, ONLY:                                         &
  environ_field_info,                                                          &
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
  lscat_zhang

USE ukca_environment_req_mod, ONLY: environ_field_index,                       &
                                    l_environ_req_available

USE ukca_environment_rdim_mod, ONLY: set_env_2d_from_0d_real,                  &
                                     set_env_2d_from_0d_integer,               &
                                     set_env_2d_from_0d_logical,               &
                                     set_env_3d_from_1d_real

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
  IF (i_field /= 0) sin_declination = field_data
CASE (fldname_equation_of_time)
  IF (i_field /= 0) equation_of_time = field_data
CASE (fldname_atmospheric_ch4)
  IF (i_field /= 0) atmospheric_ch4 = field_data
CASE (fldname_atmospheric_co2)
  IF (i_field /= 0) atmospheric_co2 = field_data
CASE (fldname_atmospheric_h2)
  IF (i_field /= 0) atmospheric_h2 = field_data
CASE (fldname_atmospheric_n2)
  IF (i_field /= 0) atmospheric_n2 = field_data
CASE (fldname_atmospheric_o2)
  IF (i_field /= 0) atmospheric_o2 = field_data
CASE (fldname_atmospheric_n2o)
  IF (i_field /= 0) atmospheric_n2o = field_data
CASE (fldname_atmospheric_cfc11)
  IF (i_field /= 0) atmospheric_cfc11 = field_data
CASE (fldname_atmospheric_cfc12)
  IF (i_field /= 0) atmospheric_cfc12 = field_data
CASE (fldname_atmospheric_cfc113)
  IF (i_field /= 0) atmospheric_cfc113 = field_data
CASE (fldname_atmospheric_hcfc22)
  IF (i_field /= 0) atmospheric_hcfc22 = field_data
CASE (fldname_atmospheric_hfc125)
  IF (i_field /= 0) atmospheric_hfc125 = field_data
CASE (fldname_atmospheric_hfc134a)
  IF (i_field /= 0) atmospheric_hfc134a = field_data
CASE (fldname_atmospheric_mebr)
  IF (i_field /= 0) atmospheric_mebr = field_data
CASE (fldname_atmospheric_mecl)
  IF (i_field /= 0) atmospheric_mecl = field_data
CASE (fldname_atmospheric_ch2br2)
  IF (i_field /= 0) atmospheric_ch2br2 = field_data
CASE (fldname_atmospheric_chbr3)
  IF (i_field /= 0) atmospheric_chbr3 = field_data
CASE (fldname_atmospheric_cfc114)
  IF (i_field /= 0) atmospheric_cfc114 = field_data
CASE (fldname_atmospheric_cfc115)
  IF (i_field /= 0) atmospheric_cfc115 = field_data
CASE (fldname_atmospheric_ccl4)
  IF (i_field /= 0) atmospheric_ccl4 = field_data
CASE (fldname_atmospheric_meccl3)
  IF (i_field /= 0) atmospheric_meccl3 = field_data
CASE (fldname_atmospheric_hcfc141b)
  IF (i_field /= 0) atmospheric_hcfc141b = field_data
CASE (fldname_atmospheric_hcfc142b)
  IF (i_field /= 0) atmospheric_hcfc142b = field_data
CASE (fldname_atmospheric_h1211)
  IF (i_field /= 0) atmospheric_h1211 = field_data
CASE (fldname_atmospheric_h1202)
  IF (i_field /= 0) atmospheric_h1202 = field_data
CASE (fldname_atmospheric_h1301)
  IF (i_field /= 0) atmospheric_h1301 = field_data
CASE (fldname_atmospheric_h2402)
  IF (i_field /= 0) atmospheric_h2402 = field_data
CASE (fldname_atmospheric_cos)
  IF (i_field /= 0) atmospheric_cos = field_data
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

INTEGER :: i_field  ! Index of field in required fields array

! Required bounds of environment field data
INTEGER :: i1
INTEGER :: i2

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

! If field is required, check the supplied field data array is allocated and
! that its bounds are compatible with the UKCA configuration.
! The field data supplied must fill the required domain but may extend beyond
! it to avoid the need for pre-trimming by the parent model.
! Allow for the possibility that the expected bounds for land only fields may
! undefined if the UKCA environment field 'land_sea_mask' has not been set.
IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(field_data)) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '1D real environment field for ''' // TRIM(varname) // ''' is unallocated'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
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
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(soil_moisture_layer1))                                 &
      ALLOCATE(soil_moisture_layer1(i1:i2))
    soil_moisture_layer1 = field_data(i1:i2)
  END IF
CASE (fldname_fland)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(fland)) ALLOCATE(fland(i1:i2))
    fland = field_data(i1:i2)
  END IF
CASE (fldname_ibvoc_isoprene)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ibvoc_isoprene)) ALLOCATE(ibvoc_isoprene(i1:i2))
    ibvoc_isoprene = field_data(i1:i2)
  END IF
CASE (fldname_ibvoc_terpene)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ibvoc_terpene)) ALLOCATE(ibvoc_terpene(i1:i2))
    ibvoc_terpene = field_data(i1:i2)
  END IF
CASE (fldname_ibvoc_methanol)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ibvoc_methanol)) ALLOCATE(ibvoc_methanol(i1:i2))
    ibvoc_methanol = field_data(i1:i2)
  END IF
CASE (fldname_ibvoc_acetone)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ibvoc_acetone)) ALLOCATE(ibvoc_acetone(i1:i2))
    ibvoc_acetone = field_data(i1:i2)
  END IF
CASE (fldname_inferno_bc)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(inferno_bc)) ALLOCATE(inferno_bc(i1:i2))
    inferno_bc = field_data(i1:i2)
  END IF
CASE (fldname_inferno_ch4)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(inferno_ch4)) ALLOCATE(inferno_ch4(i1:i2))
    inferno_ch4 = field_data(i1:i2)
  END IF
CASE (fldname_inferno_co)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(inferno_co)) ALLOCATE(inferno_co(i1:i2))
    inferno_co = field_data(i1:i2)
  END IF
CASE (fldname_inferno_nox)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(inferno_nox)) ALLOCATE(inferno_nox(i1:i2))
    inferno_nox = field_data(i1:i2)
  END IF
CASE (fldname_inferno_oc)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(inferno_oc)) ALLOCATE(inferno_oc(i1:i2))
    inferno_oc = field_data(i1:i2)
  END IF
CASE (fldname_inferno_so2)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(inferno_so2)) ALLOCATE(inferno_so2(i1:i2))
    inferno_so2 = field_data(i1:i2)
  END IF
CASE DEFAULT
  ! The named field may be a 3D field internally
  CALL set_env_3d_from_1d_real(varname, i_field, field_data, error_code)
END SELECT

IF (error_code == errcode_env_field_unknown) THEN
  IF (PRESENT(error_message)) error_message =                                  &
    'Unknown name for 1D environmental input field: ''' // TRIM(varname) // ''''
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

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

INTEGER :: i_field  ! Index of field in required fields array

! Required bounds of environment field data
INTEGER :: i1
INTEGER :: i2
INTEGER :: j1
INTEGER :: j2

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

! If field is required, check the supplied field data array is allocated and
! that its bounds are compatible with the UKCA configuration.
! The field data supplied must fill the required domain but may extend beyond
! it to avoid the need for pre-trimming (e.g. halo removal) by the parent model.
! Allow for the possibility that the expected bounds for land only fields may
! undefined if the UKCA environment field 'land_sea_mask' has not been set.
IF (i_field /= 0) THEN
  IF (.NOT. ALLOCATED(field_data)) THEN
    error_code = errcode_env_field_mismatch
    IF (PRESENT(error_message)) error_message =                                &
      '2D real environment field for ''' // TRIM(varname) // ''' is unallocated'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
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
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(latitude)) ALLOCATE(latitude(i1:i2,j1:j2))
    latitude = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_longitude)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(longitude)) ALLOCATE(longitude(i1:i2,j1:j2))
    longitude = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_sin_latitude)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(sin_latitude)) ALLOCATE(sin_latitude(i1:i2,j1:j2))
    sin_latitude = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_cos_latitude)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(cos_latitude)) ALLOCATE(cos_latitude(i1:i2,j1:j2))
    cos_latitude = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_tan_latitude)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(tan_latitude)) ALLOCATE(tan_latitude(i1:i2,j1:j2))
    tan_latitude = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_conv_cloud_lwp)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(conv_cloud_lwp)) ALLOCATE(conv_cloud_lwp(i1:i2,j1:j2))
    conv_cloud_lwp = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_tstar)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(tstar)) ALLOCATE(tstar(i1:i2,j1:j2))
    tstar = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_zbl)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(zbl)) ALLOCATE(zbl(i1:i2,j1:j2))
    zbl = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_rough_length)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(rough_length)) ALLOCATE(rough_length(i1:i2,j1:j2))
    rough_length = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_seaice_frac)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(seaice_frac)) ALLOCATE(seaice_frac(i1:i2,j1:j2))
    seaice_frac = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_frac_types)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(frac_types)) ALLOCATE(frac_types(i1:i2,j1:j2))
    frac_types = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_laift_lp)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(laift_lp)) ALLOCATE(laift_lp(i1:i2,j1:j2))
    laift_lp = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_canhtft_lp)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(canhtft_lp)) ALLOCATE(canhtft_lp(i1:i2,j1:j2))
    canhtft_lp = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_tstar_tile)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(tstar_tile)) ALLOCATE(tstar_tile(i1:i2,j1:j2))
    tstar_tile = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_z0tile_lp)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(z0tile_lp)) ALLOCATE(z0tile_lp(i1:i2,j1:j2))
    z0tile_lp = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_pstar)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(pstar)) ALLOCATE(pstar(i1:i2,j1:j2))
    pstar = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_surf_albedo)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(surf_albedo)) ALLOCATE(surf_albedo(i1:i2,j1:j2))
    surf_albedo = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_zhsc)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(zhsc)) ALLOCATE(zhsc(i1:i2,j1:j2))
    zhsc = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_u_scalar_10m)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(u_scalar_10m)) ALLOCATE(u_scalar_10m(i1:i2,j1:j2))
    u_scalar_10m = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_surf_hf)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(surf_hf)) ALLOCATE(surf_hf(i1:i2,j1:j2))
    surf_hf = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_u_s)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(u_s)) ALLOCATE(u_s(i1:i2,j1:j2))
    u_s = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_ch4_wetl_emiss)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ch4_wetl_emiss)) ALLOCATE(ch4_wetl_emiss(i1:i2,j1:j2))
    ch4_wetl_emiss = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_dms_sea_conc)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dms_sea_conc)) ALLOCATE(dms_sea_conc(i1:i2,j1:j2))
    dms_sea_conc = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_chloro_sea)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(chloro_sea)) ALLOCATE(chloro_sea(i1:i2,j1:j2))
    chloro_sea = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_dust_flux_div1)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,1) = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_dust_flux_div2)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,2) = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_dust_flux_div3)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,3) = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_dust_flux_div4)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,4) = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_dust_flux_div5)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,5) = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_dust_flux_div6)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_flux))                                            &
      ALLOCATE(dust_flux(i1:i2,j1:j2,k1_dust_flux:k2_dust_flux))
    dust_flux(:,:,6) = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_surf_wetness)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(surf_wetness))                                         &
      ALLOCATE(surf_wetness(i1:i2,j1:j2))
    surf_wetness = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_grid_surf_area)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(grid_surf_area))                                       &
      ALLOCATE(grid_surf_area(i1:i2,j1:j2))
    grid_surf_area = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_ext_cg_flash)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ext_cg_flash))                                         &
      ALLOCATE(ext_cg_flash(i1:i2,j1:j2))
    ext_cg_flash = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_ext_ic_flash)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ext_ic_flash))                                         &
      ALLOCATE(ext_ic_flash(i1:i2,j1:j2))
    ext_ic_flash = field_data(i1:i2,j1:j2)
  END IF
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
  IF (PRESENT(error_message)) error_message =                                  &
    'Unknown name for 2D real environmental input field: ''' //                &
    TRIM(varname) // ''''
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
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
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(kent)) ALLOCATE(kent(i1:i2,j1:j2))
    kent = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_kent_dsc)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(kent_dsc)) ALLOCATE(kent_dsc(i1:i2,j1:j2))
    kent_dsc = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_conv_cloud_base)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(conv_cloud_base)) ALLOCATE(conv_cloud_base(i1:i2,j1:j2))
    conv_cloud_base = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_conv_cloud_top)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(conv_cloud_top)) ALLOCATE(conv_cloud_top(i1:i2,j1:j2))
    conv_cloud_top = field_data(i1:i2,j1:j2)
  END IF
CASE (fldname_lscat_zhang)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(lscat_zhang)) ALLOCATE(lscat_zhang(i1:i2,j1:j2))
    lscat_zhang = field_data(i1:i2,j1:j2)
  END IF
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
    IF (.NOT. ALLOCATED(land_sea_mask)) ALLOCATE(land_sea_mask(i1:i2,j1:j2))
    land_sea_mask = field_data(i1:i2,j1:j2)
  END IF
  ! Special processing to set up index of land points for locating
  ! land-only environment fields on a 2-D grid
  CALL locate_land_points()
  ! Clear any existing land-only fields as these fields may be inconsistent
  ! with the new land sea mask
  CALL clear_land_only_fields()
CASE (fldname_l_tile_active)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(l_tile_active)) ALLOCATE(l_tile_active(i1:i2,j1:j2))
    l_tile_active = field_data(i1:i2,j1:j2)
  END IF
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
CASE (fldname_theta)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(theta)) ALLOCATE(theta(i1:i2,j1:j2,k1:k2))
    theta = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_q)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(q)) ALLOCATE(q(i1:i2,j1:j2,k1:k2))
    q = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_qcf)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(qcf)) ALLOCATE(qcf(i1:i2,j1:j2,k1:k2))
    qcf = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_conv_cloud_amount)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(conv_cloud_amount))                                    &
      ALLOCATE(conv_cloud_amount(i1:i2,j1:j2,k1:k2))
    conv_cloud_amount = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_rho_r2)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(rho_r2)) ALLOCATE(rho_r2(i1:i2,j1:j2,k1:k2))
    rho_r2 = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_qcl)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(qcl)) ALLOCATE(qcl(i1:i2,j1:j2,k1:k2))
    qcl = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_exner_rho_levels)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(exner_rho_levels))                                     &
      ALLOCATE(exner_rho_levels(i1:i2,j1:j2,k1:k2))
    exner_rho_levels = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_area_cloud_fraction)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(area_cloud_fraction))                                  &
      ALLOCATE(area_cloud_fraction(i1:i2,j1:j2,k1:k2))
    area_cloud_fraction = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_cloud_frac)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(cloud_frac)) ALLOCATE(cloud_frac(i1:i2,j1:j2,k1:k2))
    cloud_frac = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_cloud_liq_frac)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(cloud_liq_frac))                                       &
      ALLOCATE(cloud_liq_frac(i1:i2,j1:j2,k1:k2))
    cloud_liq_frac = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_exner_theta_levels)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(exner_theta_levels))                                   &
      ALLOCATE(exner_theta_levels(i1:i2,j1:j2,k1:k2))
    exner_theta_levels = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_p_rho_levels)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(p_rho_levels)) ALLOCATE(p_rho_levels(i1:i2,j1:j2,k1:k2))
    p_rho_levels = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_p_theta_levels)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(p_theta_levels))                                       &
      ALLOCATE(p_theta_levels(i1:i2,j1:j2,k1:k2))
    p_theta_levels = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_rhokh_rdz)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(rhokh_rdz)) ALLOCATE(rhokh_rdz(i1:i2,j1:j2,k1:k2))
    rhokh_rdz = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_dtrdz)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dtrdz))                                                &
      ALLOCATE(dtrdz(i1:i2,j1:j2,k1:k2))
    dtrdz = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_we_lim)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(we_lim)) ALLOCATE(we_lim(i1:i2,j1:j2,k1:k2))
    we_lim = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_t_frac)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(t_frac)) ALLOCATE(t_frac(i1:i2,j1:j2,k1:k2))
    t_frac = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_zrzi)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(zrzi)) ALLOCATE(zrzi(i1:i2,j1:j2,k1:k2))
    zrzi = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_we_lim_dsc)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(we_lim_dsc)) ALLOCATE(we_lim_dsc(i1:i2,j1:j2,k1:k2))
    we_lim_dsc = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_t_frac_dsc)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(t_frac_dsc)) ALLOCATE(t_frac_dsc(i1:i2,j1:j2,k1:k2))
    t_frac_dsc = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_zrzi_dsc)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(zrzi_dsc)) ALLOCATE(zrzi_dsc(i1:i2,j1:j2,k1:k2))
    zrzi_dsc = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_stcon)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(stcon)) ALLOCATE(stcon(i1:i2,j1:j2,k1:k2))
    stcon = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_ls_rain3d)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ls_rain3d)) ALLOCATE(ls_rain3d(i1:i2,j1:j2,k1:k2))
    ls_rain3d = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_ls_snow3d)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ls_snow3d)) ALLOCATE(ls_snow3d(i1:i2,j1:j2,k1:k2))
    ls_snow3d = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_autoconv)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(autoconv)) ALLOCATE(autoconv(i1:i2,j1:j2,k1:k2))
    autoconv = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_accretion)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(accretion)) ALLOCATE(accretion(i1:i2,j1:j2,k1:k2))
    accretion = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_pv_on_theta_mlevs)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(pv_on_theta_mlevs))                                    &
      ALLOCATE(pv_on_theta_mlevs(i1:i2,j1:j2,k1:k2))
    pv_on_theta_mlevs = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_conv_rain3d)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(conv_rain3d)) ALLOCATE(conv_rain3d(i1:i2,j1:j2,k1:k2))
    conv_rain3d = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_conv_snow3d)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(conv_snow3d)) ALLOCATE(conv_snow3d(i1:i2,j1:j2,k1:k2))
    conv_snow3d = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_so4_sa_clim)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(so4_sa_clim)) ALLOCATE(so4_sa_clim(i1:i2,j1:j2,k1:k2))
    so4_sa_clim = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_so4_aitken)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(so4_aitken)) ALLOCATE(so4_aitken(i1:i2,j1:j2,k1:k2))
    so4_aitken = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_so4_accum)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(so4_accum)) ALLOCATE(so4_accum(i1:i2,j1:j2,k1:k2))
    so4_accum = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_soot_fresh)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(soot_fresh)) ALLOCATE(soot_fresh(i1:i2,j1:j2,k1:k2))
    soot_fresh = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_soot_aged)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(soot_aged)) ALLOCATE(soot_aged(i1:i2,j1:j2,k1:k2))
    soot_aged = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_ocff_fresh)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ocff_fresh)) ALLOCATE(ocff_fresh(i1:i2,j1:j2,k1:k2))
    ocff_fresh = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_ocff_aged)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ocff_aged)) ALLOCATE(ocff_aged(i1:i2,j1:j2,k1:k2))
    ocff_aged = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_biogenic)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(biogenic)) ALLOCATE(biogenic(i1:i2,j1:j2,k1:k2))
    biogenic = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_dust_div1)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_div1)) ALLOCATE(dust_div1(i1:i2,j1:j2,k1:k2))
    dust_div1 = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_dust_div2)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_div2)) ALLOCATE(dust_div2(i1:i2,j1:j2,k1:k2))
    dust_div2 = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_dust_div3)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_div3)) ALLOCATE(dust_div3(i1:i2,j1:j2,k1:k2))
    dust_div3 = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_dust_div4)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_div4)) ALLOCATE(dust_div4(i1:i2,j1:j2,k1:k2))
    dust_div4 = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_dust_div5)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_div5)) ALLOCATE(dust_div5(i1:i2,j1:j2,k1:k2))
    dust_div5 = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_dust_div6)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(dust_div6)) ALLOCATE(dust_div6(i1:i2,j1:j2,k1:k2))
    dust_div6 = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_sea_salt_film)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(sea_salt_film))                                        &
      ALLOCATE(sea_salt_film(i1:i2,j1:j2,k1:k2))
    sea_salt_film = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_sea_salt_jet)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(sea_salt_jet)) ALLOCATE(sea_salt_jet(i1:i2,j1:j2,k1:k2))
    sea_salt_jet = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_co2_interactive)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(co2_interactive))                                      &
      ALLOCATE(co2_interactive(i1:i2,j1:j2,k1:k2))
    co2_interactive = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_rim_cry)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(rim_cry)) ALLOCATE(rim_cry(i1:i2,j1:j2,k1:k2))
    rim_cry = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_rim_agg)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(rim_agg)) ALLOCATE(rim_agg(i1:i2,j1:j2,k1:k2))
    rim_agg = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_vertvel)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(vertvel)) ALLOCATE(vertvel(i1:i2,j1:j2,k1:k2))
    vertvel = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_bl_tke)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(bl_tke)) ALLOCATE(bl_tke(i1:i2,j1:j2,k1:k2))
    bl_tke = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_interf_z)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(interf_z)) ALLOCATE(interf_z(i1:i2,j1:j2,k1:k2))
    interf_z = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_h2o2_offline)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(h2o2_offline))                                         &
      ALLOCATE(h2o2_offline(i1:i2,j1:j2,k1:k2))
    h2o2_offline = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_ho2_offline)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(ho2_offline)) ALLOCATE(ho2_offline(i1:i2,j1:j2,k1:k2))
    ho2_offline = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_no3_offline)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(no3_offline)) ALLOCATE(no3_offline(i1:i2,j1:j2,k1:k2))
    no3_offline = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_o3_offline)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(o3_offline)) ALLOCATE(o3_offline(i1:i2,j1:j2,k1:k2))
    o3_offline = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE (fldname_oh_offline)
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(oh_offline)) ALLOCATE(oh_offline(i1:i2,j1:j2,k1:k2))
    oh_offline = field_data(i1:i2,j1:j2,k1:k2)
  END IF
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
  IF (PRESENT(error_message)) error_message =                                  &
    'Unknown name for 3D environmental input field: ''' // TRIM(varname) //    &
    ''''
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
INTEGER :: l1
INTEGER :: l2

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
  l1 = environ_field_info(i_field)%lbound_dim4
  l2 = environ_field_info(i_field)%ubound_dim4
  IF (LBOUND(field_data,DIM=1) > i1 .OR. UBOUND(field_data,DIM=1) < i2 .OR.    &
      LBOUND(field_data,DIM=2) > j1 .OR. UBOUND(field_data,DIM=2) < j2 .OR.    &
      LBOUND(field_data,DIM=3) > k1 .OR. UBOUND(field_data,DIM=3) < k2 .OR.    &
      LBOUND(field_data,DIM=4) > l1 .OR. UBOUND(field_data,DIM=4) < l2) THEN
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
  IF (i_field /= 0) THEN
    IF (.NOT. ALLOCATED(photol_rates)) THEN
      ALLOCATE(photol_rates(i1:i2,j1:j2,k1:k2,l1:l2))
    END IF
    photol_rates = field_data(i1:i2,j1:j2,k1:k2,l1:l2)
  END IF
CASE DEFAULT
  ! Error: Not a recognised field
  error_code = errcode_env_field_unknown
  IF (PRESENT(error_message)) error_message =                                  &
    'Unknown name for 4D environmental input field: ''' // TRIM(varname) //    &
    ''''
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
SUBROUTINE clear_environment_fields()
! ----------------------------------------------------------------------
! Description:
!   Resets scalar fields values, deallocates environment field arrays
!   and land point index array and updates availability flags
! ----------------------------------------------------------------------

IMPLICIT NONE

! Local variables

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CLEAR_ENVIRONMENT_FIELDS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

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

! Update field availability
l_environ_field_available(:) = .FALSE.

! The land sea mask is no longer valid so clear land point indices
IF (ALLOCATED(land_index)) DEALLOCATE(land_index)
land_points = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE clear_environment_fields

END MODULE ukca_environment_mod
