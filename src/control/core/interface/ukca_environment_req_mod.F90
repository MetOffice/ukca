! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module for handling UKCA's environmental driver requirement.
!   Environmental drivers are input fields defined on the current UKCA
!   model grid that may be varied by the parent application during the run.
!
!   The module provides the following procedures for the UKCA API.
!
!     ukca_get_environment_varlist   - Returns list of names of required fields.
!     ukca_get_envgroup_varlists     - Returns lists of names of required fields
!                                      by group.
!
!   The following additional public procedures are provided for use within UKCA.
!
!     init_environment_req           - Determines environment data requirement
!     environ_field_index            - Return index of a specific environment
!                                      field in the list of names of required
!                                      fields.
!     check_environment_availability - Checks availability of required
!                                      environment fields.
!     environ_field_available        - Return T or F depending on whether a
!                                      specific environment field is available.
!     print_environment_summary      - Writes a summary of the current
!                                      environment data to the log file.
!     clear_environment_req          - Resets all data relating to environment
!                                      data requirement to its initial state
!                                      for a new UKCA configuration.
!
!   The module also provides a public logical 'l_environ_req_available'
!   that indicates the availability status of the environment data requirement
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

MODULE ukca_environment_req_mod

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

USE ukca_fieldname_mod,  ONLY: maxlen_fieldname,                               &
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
  fldname_lscat_zhang,                                                         &
  fldname_grid_area_fullht,                                                    &
  fldname_grid_volume

USE ukca_environment_fields_mod, ONLY: environ_field_info, environ_field_ptrs, &
                                       l_environ_field_available,              &
                                       no_bound_value, env_field_info_type,    &
                                       env_field_ptrs_type

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_env_req_uninit,                              &
                          errcode_env_field_mismatch,                          &
                          errcode_ukca_internal_fault

IMPLICIT NONE

PRIVATE

! Public procedures
PUBLIC init_environment_req, ukca_get_environment_varlist,                     &
       ukca_get_envgroup_varlists, environ_field_index,                        &
       check_environment_availability, environ_field_available,                &
       print_environment_summary, clear_environment_req

! Public flag for use within UKCA to indicate whether the environment fields
! requirement has been initialised
LOGICAL, SAVE, PUBLIC :: l_environ_req_available = .FALSE.

! Field group codes for collating driver fields in arrays, as required for
! the API routine 'ukca_step_control' (an alternative to using separate
! 'ukca_set_environment' and 'ukca_step' calls):
! 'scalar' group comprises fields defined by a scalar value internally
! 'flat' groups comprise fields defined on a flat spatial grid with 2D
! representation internally (input can have reduced dimension)
! 'flatpft' group has an additional dimension for plant functional type tiles
! 'fullht' group comprises fields defined on a full height spatial grid with 3D
! representation internally (input can have reduced dimension)
! 'fullht0' group comprises fields defined on an extended full height spatial
! grid with an additional level 0 below and 3D representation internally
! (input can have reduced dimension)
! 'fullhtp1' group comprises fields defined on an extended full height spatial
! grid with an additional level above and 3D representation internally
! (input can have reduced dimension)
! 'bllev' group comprises fields defined on a spatial grid for boundary
! layer levels with 3D representation internally (input can have reduced
! dimension)
! 'entlev' group comprises fields defined on a spatial grid for entrainment
! levels used in tr_mix with 3D representation internally (input can have
! reduced dimension)
! 'land' group comprises fields defined on land points only
! 'landtile' and 'landpft' groups comprise fields defined on land-point tiles
INTEGER, PARAMETER :: group_undefined = 0         ! Not assigned to a group
INTEGER, PARAMETER :: group_scalar_real = 1       ! Scalar
INTEGER, PARAMETER :: group_flat_integer = 2      ! 2D spatial integer
INTEGER, PARAMETER :: group_flat_real = 3         ! 2D spatial real
INTEGER, PARAMETER :: group_flat_logical = 4      ! 2D spatial logical
INTEGER, PARAMETER :: group_flatpft_real = 5      ! 3D real on PFT tiles on 2D
                                                  ! spatial grid
INTEGER, PARAMETER :: group_fullht_real = 6       ! 3D spatial real on all
                                                  ! model levels (theta levels)
INTEGER, PARAMETER :: group_fullht0_real = 7      ! 3D spatial real on all model
                                                  ! levels + zero level below
INTEGER, PARAMETER :: group_fullhtp1_real = 8     ! 3D spatial real on all model
                                                  ! levels + one level above
INTEGER, PARAMETER :: group_bllev_real = 9        ! 3D spatial real on boundary
                                                  ! layer levels
INTEGER, PARAMETER :: group_entlev_real = 10      ! 3D spatial real at
                                                  ! entrainment_levels
INTEGER, PARAMETER :: group_land_real = 11        ! 1D real on land points
INTEGER, PARAMETER :: group_landtile_real = 12    ! 2D real on land pt. tiles
INTEGER, PARAMETER :: group_landtile_logical = 13 ! 2D logical on land pt. tiles
INTEGER, PARAMETER :: group_landpft_real = 14     ! 2D real on PFT tiles on
                                                  ! land points
INTEGER, PARAMETER :: group_fullhtphot_real = 15  ! 4D: spatial on all model
                                                  ! levels and photolysis reacn
! Any environmental driver not assigned to a group will be ignored by
! API routines 'ukca_get_envgroup_varlists' and 'ukca_step_control' but will
! appear in the full list returned by 'ukca_get_environment_varlist' and must be
! set via a separate 'ukca_set_environment' call if required.

! List of environment fields required for the current UKCA configuration
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames(:)

! Lists of environment fields required for the current UKCA configuration
! by subgroup
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_scalar_real(:)      ! Field names of scalars
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_flat_integer(:)     ! Field names of 2D spatial
                                             ! integers
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_flat_real(:)        ! Field names of 2D spatial reals
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_flat_logical(:)     ! Field names of 2D spatial
                                             ! logicals
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_flatpft_real(:)     ! Field names of 3D reals on plant
                                             ! functional type tiles on 2D
                                             ! spatial grid
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_fullht_real(:)      ! Field names of full height 3D
                                             ! spatial reals
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_fullht0_real(:)     ! Field names of 3D spatial reals
                                             ! on all model levels plus zero
                                             ! level below
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_fullhtp1_real(:)    ! Field names of 3D spatial reals
                                             ! on all model levels plus one
                                             ! above
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_bllev_real(:)       ! Field names of 3D spatial reals
                                             ! on boundary layer levels
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_entlev_real(:)      ! Field names of 3D spatial reals
                                             ! on entrainment layer levels
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_land_real(:)        ! Field names of 1D reals on land
                                             ! points
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_landtile_real(:)    ! Field names of 2D reals on land-
                                             ! point tiles
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_landtile_logical(:) ! Field names of 2D logicals on
                                             ! land-point tiles
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_landpft_real(:)     ! Field names of 2D reals on plant
                                             ! functional type tiles on land
                                             ! points
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, TARGET, SAVE ::                  &
  environ_field_varnames_fullhtphot_real(:)  ! Field names of full height 3D
                                             ! spatial + photol species reals

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_ENVIRONMENT_REQ_MOD'


CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE init_environment_req(ukca_config, glomap_config,                    &
                                speci, advt, lbc_spec, cfc_lumped,             &
                                em_chem_spec, ctype, error_code,               &
                                error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Determines the environment data required for the current UKCA
!   configuration.
!
! Method:
!   Create and save a reference list containing the names of the
!   required environment fields and the corresponding field info and
!   field availability arrays.
!   The list entries are selected by determining, for each
!   UKCA-recognised environmental driver field in turn, whether there
!   is a requirement for that field. The requirement is established
!   with reference to the given UKCA configuration specification (passed
!   as input arguments 'ukca_config' and 'glomap_config'), species lists
!   defining the species in the selected chemistry scheme ('speci') and
!   the subsets of species that are treated as tracers ('advt') or have
!   prescribed lower boundary conditions in stratospheric schemes
!   (listed in 'lbc_spec' or 'cfc_lumped'), or defined as a prescribed
!   field (CF in 'ctype'), or if land surface emissions are required by
!   this configuration (listed in 'em_chem_spec').
! ----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: ukca_config_spec_type,                &
                                         glomap_config_spec_type,              &
                                         i_ukca_fastjx,                        &
                                         i_ukca_photol_strat,                  &
                                         i_ukca_nophot,                        &
                                         i_light_param_pr,                     &
                                         i_light_param_luhar,                  &
                                         i_strat_lbc_env,                      &
                                         i_strat_lbc_off,                      &
                                         i_ukca_activation_arg,                &
                                         i_light_param_ext,                    &
                                         i_light_param_off,                    &
                                         i_top_bc

USE ukca_chem_defs_mod, ONLY: ratj_defs
! Need ratj_defs to calculate the dimensions of photol_rates array

USE ukca_um_legacy_mod, ONLY: tdims
!!!! tdims required to replicate results impacted by cloud_frac level offset
!!!! bug pending removal of temporary logical l_fix_ukca_cloud_frac

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

! Subroutine arguments
TYPE(ukca_config_spec_type), INTENT(IN) :: ukca_config
TYPE(glomap_config_spec_type), INTENT(IN) :: glomap_config
CHARACTER(LEN=*), INTENT(IN) :: speci(:)      ! Names of all active species
CHARACTER(LEN=*), INTENT(IN) :: advt(:)       ! Advected species
CHARACTER(LEN=*), INTENT(IN) :: lbc_spec(:)   ! Species requiring lower boundary
CHARACTER(LEN=*), INTENT(IN) :: cfc_lumped(:) ! CFCs lumped into major LBC ones
CHARACTER(LEN=*), INTENT(IN) :: em_chem_spec(:)  ! Species requiring emissions
CHARACTER(LEN=*), INTENT(IN) :: ctype(:)      ! Type of chemical species

INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

! Field counts
INTEGER, PARAMETER :: n_max = 144  ! Maximum number of environment fields
INTEGER :: n                       ! Count of environment fields selected
INTEGER :: i                       ! Counter for loops

! Temporary field name array and corresponding field info for use in
! collating data that will subsequently be copied to the field list arrays
! having the correct size allocation
CHARACTER(LEN=maxlen_fieldname) :: fld_names(n_max) = ''
TYPE(env_field_info_type) :: fld_info(n_max)

! Default values for field info
TYPE(env_field_info_type) :: fld_info_default

! Requirement for solar angle calculations
LOGICAL :: l_solang

! Requirement for external lower BC values for a stratospheric chemistry scheme
LOGICAL :: l_req_strat_lbc

! Requirement of JULES emission field - temp variables to simplify the logic
LOGICAL :: l_req_emiss
CHARACTER(LEN=maxlen_fieldname) :: use_fldname = ''

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'INIT_ENVIRONMENT_REQ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Ensure all environment-related data are in uninitialised state
IF (l_environ_req_available) CALL clear_environment_req()

! If run does not include chemistry then no fields are required so set array
! sizes to zero, set public flag to show that the environment fields
! requirement is available and return
IF (.NOT. ukca_config%l_ukca_chem) THEN
  ALLOCATE(environ_field_varnames(0))
  ALLOCATE(environ_field_info(0))
  ALLOCATE(environ_field_ptrs(0))
  ALLOCATE(l_environ_field_available(0))
  l_environ_req_available = .TRUE.
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Setup a temporary field info array with default bounds matching
! the UKCA model grid
fld_info_default%group = group_fullht_real
fld_info_default%lbound_dim1 = 1
fld_info_default%ubound_dim1 = ukca_config%row_length
fld_info_default%lbound_dim2 = 1
fld_info_default%ubound_dim2 = ukca_config%rows
fld_info_default%lbound_dim3 = 1
fld_info_default%ubound_dim3 = ukca_config%model_levels
fld_info_default%lbound_dim4 = 1
fld_info_default%ubound_dim4 = 1
fld_info_default%l_land_only = .FALSE.
fld_info(:) = fld_info_default

! Add each required field to temporary field name array below.
! The order in which fields are added determines the order in which they
! will appear in the required field list that is accessible to a parent
! application via the 'ukca_get_environment_varlist' API call.
! By convention:
! - Fields should be in group order.
! - Fields that are required in all configurations should appear first in
!   each group.
! Note that the land sea mask preceeds any land-only fields because it
! determines the expected size of these fields and must be available when they
! are set. If the parent calls 'ukca_set_environment' for fields in the order
! listed, the order of the groups ensures that this dependency is satisfied.

! Set logical to show whether solar angle calculations are required.
! These are required for application of a diurnal cycle to offline oxidants
! if offline oxidants chemistry is selected or to isoprene emissions.
! They are also required for photolysis if Fast-JX and/or stratospheric
! photolysis scheme is used. This is temporary pending completion of the
! separation of photolysis from UKCA. (In that case, they are not used
! internally.)
l_solang = ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be .OR.  &
           (ukca_config%l_diurnal_isopems .AND.                                &
            .NOT. ukca_config%l_ukca_emissions_off) .OR.                       &
           ukca_config%i_ukca_photol == i_ukca_fastjx .OR.                     &
           ukca_config%i_ukca_photol == i_ukca_photol_strat

n = 0

! -- Environmental drivers in scalar group --

! Sin(declination) and equation of time are required for solar angle
! calculations.
IF (l_solang) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_sin_declination
    fld_info(n)%group = group_scalar_real
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_equation_of_time
    fld_info(n)%group = group_scalar_real
  END IF
END IF

! The following scalar gas mixing ratios are required if the option to use an
! external value for driving the chemistry is set or we are using external
! lower boundary conditions for a stratospheric chemistry scheme.
! Additionally, CH4 may be required as a prescribed value in place of CH4
! emissions.

l_req_strat_lbc = (ukca_config%i_strat_lbc_source == i_strat_lbc_env)

IF ((l_req_strat_lbc .AND. ANY(lbc_spec(:) == 'CH4       ')) .OR.              &
     ukca_config%l_chem_environ_ch4_scalar                   .OR.              &
    (ukca_config%l_ukca_prescribech4 .AND. ANY(advt(:) == 'CH4       '))) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_atmospheric_ch4
    fld_info(n)%group = group_scalar_real
  END IF
END IF

IF ((l_req_strat_lbc .AND. ANY(lbc_spec(:) == 'CO2       ')) .OR.              &
     ukca_config%l_chem_environ_co2_scalar) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_atmospheric_co2
    fld_info(n)%group = group_scalar_real
  END IF
END IF

IF ((l_req_strat_lbc .AND. ANY(lbc_spec(:) == 'H2        ')) .OR.              &
     ukca_config%l_chem_environ_h2_scalar) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_atmospheric_h2
    fld_info(n)%group = group_scalar_real
  END IF
END IF

IF ((l_req_strat_lbc .AND. ANY(lbc_spec(:) == 'N2        ')) .OR.              &
     ukca_config%l_chem_environ_n2_scalar) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_atmospheric_n2
    fld_info(n)%group = group_scalar_real
  END IF
END IF

IF ((l_req_strat_lbc .AND. ANY(lbc_spec(:) == 'O2        ')) .OR.              &
     ukca_config%l_chem_environ_o2_scalar) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_atmospheric_o2
    fld_info(n)%group = group_scalar_real
  END IF
END IF

! These scalar gas mixing ratios are required if using external lower
! boundary conditions for a stratospheric chemistry scheme.
IF (l_req_strat_lbc) THEN
  IF (ANY(lbc_spec == 'N2O       ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_n2o
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CFCl3     ') .OR. ANY(cfc_lumped == 'CFCl3     ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_cfc11
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CF2Cl2    ') .OR. ANY(cfc_lumped == 'CF2Cl2    ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_cfc12
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CF2ClCFCl2') .OR. ANY(cfc_lumped == 'CF2ClCFCl2')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_cfc113
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CHF2Cl    ') .OR. ANY(cfc_lumped == 'CHF2Cl    ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_hcfc22
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CF3CHF2   ') .OR. ANY(cfc_lumped == 'CF3CHF2   ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_hfc125
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CHF2FCF3  ') .OR. ANY(cfc_lumped == 'CHF2FCF3  ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_hfc134a
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'MeBr      ') .OR. ANY(cfc_lumped == 'MeBr      ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_mebr
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CF2Cl2    ') .OR. ANY(cfc_lumped == 'CF2Cl2    ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_mecl
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CH2Br2    ') .OR. ANY(cfc_lumped == 'CH2Br2    ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_ch2br2
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CHBr3     ') .OR. ANY(cfc_lumped == 'CHBr3     ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_chbr3
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CF2ClCF2Cl') .OR. ANY(cfc_lumped == 'CF2ClCF2Cl')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_cfc114
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CF2ClCF3  ') .OR. ANY(cfc_lumped == 'CF2ClCF3  ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_cfc115
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CCl4      ') .OR. ANY(cfc_lumped == 'CCl4      ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_ccl4
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'MeCCl3    ') .OR. ANY(cfc_lumped == 'MeCCl3    ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_meccl3
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'MeCFCl2   ') .OR. ANY(cfc_lumped == 'MeCFCl2   ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_hcfc141b
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'MeCF2Cl   ') .OR. ANY(cfc_lumped == 'MeCF2Cl   ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_hcfc142b
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CF2ClBr   ') .OR. ANY(cfc_lumped == 'CF2ClBr   ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_h1211
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CF2Br2    ') .OR. ANY(cfc_lumped == 'CF2Br2    ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_h1202
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CF3Br     ') .OR. ANY(cfc_lumped == 'CF3Br     ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_h1301
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'CF2BrCF2Br') .OR. ANY(cfc_lumped == 'CF2BrCF2Br')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_h2402
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
  IF (ANY(lbc_spec == 'COS       ')) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_atmospheric_cos
      fld_info(n)%group = group_scalar_real
    END IF
  END IF
END IF  ! l_strat_req_lbc = True

! -- Environmental drivers in flat grid integer group --

! Entrainment level fields are always required unless emissions are off or
! tracer updates from emissions and boundary layer mixing are suppressed
IF (.NOT. (ukca_config%l_ukca_emissions_off .OR.                               &
          ukca_config%l_suppress_ems)) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_kent
    fld_info(n)%group = group_flat_integer
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_kent_dsc
    fld_info(n)%group = group_flat_integer
  END IF
END IF

! Calculation of lightning NOx emissions requires convection diagnostics
! e.g. from a parameterized convection scheme in the parent model.
! These are also required for photolysis if Fast-JX is used. This is
! temporary pending completion of the separation of photolysis from UKCA.
! (In that case, they are not used internally.)
IF (ukca_config%i_ukca_light_param == i_light_param_pr .OR.                    &
    ukca_config%i_ukca_light_param == i_light_param_luhar .OR.                 &
    ukca_config%i_ukca_photol == i_ukca_fastjx) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_conv_cloud_base
    fld_info(n)%group = group_flat_integer
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_conv_cloud_top
    fld_info(n)%group = group_flat_integer
  END IF
END IF

! Index of dominant surface category in grid-box as per Zhang - required for
! GLOMAP if logical for 'new' method is On.
IF (glomap_config%l_improve_aero_drydep) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_lscat_zhang
    fld_info(n)%group = group_flat_integer
  END IF
END IF


! -- Environmental drivers in flat grid real group --

n = n + 1
IF (n <= n_max) THEN
  fld_names(n) = fldname_pstar
  fld_info(n)%group = group_flat_real
END IF

! Latitude is required if a top boundary scheme denoted by code 'i_top_bc' or
! higher is in use. In addition, it is required for calculation of tropopause
! height (when not using a fixed tropopause level), for both dry and wet
! deposition (if not switched off) and for any internally calculated lightning
! emissions.
IF (ukca_config%i_ukca_topboundary >= i_top_bc .OR.                            &
    .NOT. (ukca_config%l_fix_tropopause_level .AND.                            &
          ukca_config%l_ukca_drydep_off .AND.                                  &
          ukca_config%l_ukca_wetdep_off .AND.                                  &
          (ukca_config%l_ukca_emissions_off .OR.                               &
           ukca_config%i_ukca_light_param == i_light_param_off .OR.            &
           ukca_config%i_ukca_light_param == i_light_param_ext))) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_latitude
    fld_info(n)%group = group_flat_real
  END IF
END IF

! Longitude is always required unless emissions are off since it may be
! required for applying hourly scaling to emissions in general or applying a
! diurnal cycle to isoprene emissions.
! If emissions are off, it is still required if solar angle calculations are
! needed (i.e. for application of a diurnal cycle to offline oxidants when
! offline oxidants chemistry is selected) or for dry deposition when using
! the non-interactive scheme.
IF ((.NOT. ukca_config%l_ukca_emissions_off) .OR. l_solang .OR.                &
    .NOT. (ukca_config%l_ukca_intdd .OR. ukca_config%l_ukca_drydep_off)) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_longitude
    fld_info(n)%group = group_flat_real
  END IF
END IF

! Cos(latitude) is required for solar angle calculations
IF (l_solang) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_cos_latitude
    fld_info(n)%group = group_flat_real
  END IF
END IF

! Sin(latitude) is required for solar angle calculations.
! It is also required if the interactive dry deposition scheme is selected.
IF (l_solang .OR. ukca_config%l_ukca_intdd) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_sin_latitude
    fld_info(n)%group = group_flat_real
  END IF
END IF

! Tan(latitude) is required for solar angle calculations.
! It is also required for dry deposition when using the non-interactive scheme.
IF (l_solang .OR.                                                              &
    .NOT. (ukca_config%l_ukca_intdd .OR. ukca_config%l_ukca_drydep_off)) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_tan_latitude
    fld_info(n)%group = group_flat_real
  END IF
END IF

! Boundary layer height is required for dry deposition if using the interactive
! scheme, for the GLOMAP-mode boundary layer nucleation scheme (if active) or
! for emissions unless tracer updates from emissions and boundary layer mixing
! are suppressed
IF (ukca_config%l_ukca_intdd .OR. glomap_config%l_mode_bln_on .OR.             &
    .NOT. (ukca_config%l_ukca_emissions_off .OR.                               &
           ukca_config%l_suppress_ems)) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_zbl
    fld_info(n)%group = group_flat_real
  END IF
END IF

! Surface temperature is required for dry deposition or marine DMS emissions
IF ((.NOT. ukca_config%l_ukca_drydep_off) .OR. ukca_config%l_seawater_dms) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_tstar
    fld_info(n)%group = group_flat_real
  END IF
END IF

! Sea ice fraction is required for dry deposition if using the interactive
! scheme or if GLOMAP-mode dry deposition is active without the surface type
! being supplied by the parent. It is also required if marine DMS emissions or
! sea-salt emissions are on.
IF (ukca_config%l_ukca_intdd .OR.                                              &
    (glomap_config%l_ddepaer .AND.                                             &
    .NOT. glomap_config%l_improve_aero_drydep) .OR.                            &
    ukca_config%l_seawater_dms .OR. glomap_config%l_ukca_primss) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_seaice_frac
    fld_info(n)%group = group_flat_real
  END IF
END IF

! Fields in flat real group that are required for dry deposition
IF (.NOT. ukca_config%l_ukca_drydep_off) THEN

  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_rough_length
    fld_info(n)%group = group_flat_real
  END IF

  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_u_s
    fld_info(n)%group = group_flat_real
  END IF

  ! Surface sensible heat flux is required if using the interactive
  ! dry deposition scheme and surface wetness is also required if its
  ! impact on dry deposition is to be considered
  IF (ukca_config%l_ukca_intdd) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_surf_hf
      fld_info(n)%group = group_flat_real
    END IF
    IF (ukca_config%l_ukca_dry_dep_so2wet) THEN
      n = n + 1
      IF (n <= n_max) THEN
        fld_names(n) = fldname_surf_wetness
        fld_info(n)%group = group_flat_real
      END IF
    END IF
  END IF

END IF  ! .NOT. ukca_config%l_ukca_drydep_off

! Convective cloud liquid water path and surface albedo are required for
! for photolysis if Fast-JX is used. This is temporary pending completion of
! the separation of photolysis from UKCA. (They are not used internally.)
IF (ukca_config%i_ukca_photol == i_ukca_fastjx) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_conv_cloud_lwp
    fld_info(n)%group = group_flat_real
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_surf_albedo
    fld_info(n)%group = group_flat_real
  END IF
END IF

! Fields in flat grid real group that are required for emissions
IF (.NOT. ukca_config%l_ukca_emissions_off) THEN

  ! Height at top of decoupled stratocumulus layer is always required for
  ! emissions unless tracer updates from emissions and boundary layer mixing
  ! are suppressed
  IF (.NOT. ukca_config%l_suppress_ems) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_zhsc
      fld_info(n)%group = group_flat_real
    END IF
  END IF

  ! 10m wind speed is required if seawater-DMS or sea-salt emissions are enabled
  IF (ukca_config%l_seawater_dms .OR. glomap_config%l_ukca_primss) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_u_scalar_10m
      fld_info(n)%group = group_flat_real
    END IF
  END IF

  ! Gridbox surface area is required for lightning emissions, for application
  ! of lower boundary conditions in stratospheric chemistry schemes, for
  ! emissions diagnostics and for converting offline emissions provided in
  ! gridbox units.
  ! Note: No emissions diags. are currently supported independently of the UM
  IF (ukca_config%i_ukca_light_param /= i_light_param_off .OR.                 &
      (ukca_config%i_strat_lbc_source /= imdi .AND.                            &
       ukca_config%i_strat_lbc_source /= i_strat_lbc_off) .OR.                 &
      ukca_config%l_enable_diag_um .OR.                                        &
      ukca_config%l_support_ems_gridbox_units) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_grid_surf_area
      fld_info(n)%group = group_flat_real
    END IF
  END IF

  ! CH4 wetland flux is required for online wetland CH4 emissions option
  IF (ukca_config%l_ukca_qch4inter) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_ch4_wetl_emiss
      fld_info(n)%group = group_flat_real
    END IF
  END IF

  ! DMS conc. in seawater is required if marine DMS emissions are to be modelled
  IF (ukca_config%l_seawater_dms) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_dms_sea_conc
      fld_info(n)%group = group_flat_real
    END IF
  END IF

  ! Ocean near-surface chlorophyll required if online marine organic carbon
  ! emissions are on in GLOMAP-mode
  IF (glomap_config%l_ukca_prim_moc) THEN
    n = n + 1
    IF (n <= n_max) THEN
      fld_names(n) = fldname_chloro_sea
      fld_info(n)%group = group_flat_real
    END IF
  END IF

  ! Dust emissions by size bin are required if dust is modelled in GLOMAP-mode
  IF (glomap_config%l_ukca_primdu) THEN
    n = n + glomap_config%n_dust_emissions
    IF (n <= n_max) THEN
      IF (glomap_config%n_dust_emissions == 6) THEN
        fld_names(n-5) = fldname_dust_flux_div1
        fld_names(n-4) = fldname_dust_flux_div2
        fld_names(n-3) = fldname_dust_flux_div3
        fld_names(n-2) = fldname_dust_flux_div4
        fld_names(n-1) = fldname_dust_flux_div5
        fld_names(n) = fldname_dust_flux_div6
        fld_info(n-5:n)%group = group_flat_real
      ELSE
        error_code = errcode_env_req_uninit
        IF (PRESENT(error_message)) error_message =                            &
          'Unexpected number of dust emissions required'
        IF (PRESENT(error_routine)) error_routine = RoutineName
        IF (lhook)                                                             &
          CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
        RETURN
      END IF
    END IF
  END IF

  ! External lightning flash sources
  IF (ukca_config%i_ukca_light_param == i_light_param_ext) THEN
    n = n + 2
    IF (n <= n_max) THEN
      fld_names(n-1) = fldname_ext_cg_flash
      fld_names(n)   = fldname_ext_ic_flash
      fld_info(n-1:n)%group = group_flat_real
    END IF
  END IF

END IF  !.NOT. ukca_config%l_ukca_emissions_off

! -- Environmental drivers in flat grid logical group --

! Land-sea mask is required for emissions, interactive dry deposition acheme
! and aerosol dry deposition.
! It is also required for photolysis if Fast-JX is used. This is
! temporary pending completion of the separation of photolysis from UKCA.
IF ((.NOT. ukca_config%l_ukca_emissions_off) .OR.                              &
    ukca_config%l_ukca_intdd .OR. glomap_config%l_ddepaer .OR.                 &
    ukca_config%i_ukca_photol == i_ukca_fastjx) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_land_sea_mask
    fld_info(n)%group = group_flat_logical
  END IF
END IF

! -- Environmental drivers in flat grid plant functional type tile group --

! Stomatal conductance is required if using the interactive
! dry deposition scheme
IF (ukca_config%l_ukca_intdd) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_stcon
    fld_info(n)%group = group_flatpft_real
    fld_info(n)%ubound_dim3 = ukca_config%npft
  END IF
END IF

! -- Environmental drivers in full-height grid group --

n = n + 1
IF (n <= n_max) fld_names(n) = fldname_theta
n = n + 1
IF (n <= n_max) fld_names(n) = fldname_q
n = n + 1
IF (n <= n_max) fld_names(n) = fldname_qcf
n = n + 1
IF (n <= n_max) fld_names(n) = fldname_qcl
n = n + 1
IF (n <= n_max) fld_names(n) = fldname_exner_theta_levels
n = n + 1
IF (n <= n_max) fld_names(n) = fldname_p_rho_levels
n = n + 1
IF (n <= n_max) fld_names(n) = fldname_p_theta_levels
n = n + 1
IF (n <= n_max) THEN
  fld_names(n) = fldname_cloud_frac
  IF (.NOT. ukca_config%l_fix_ukca_cloud_frac) THEN
    fld_info(n)%lbound_dim3 = tdims%k_start
    !!!! Required to replicate results impacted by cloud_frac level offset bug
    !!!! pending removal of temporary logical l_fix_ukca_cloud_frac
  END IF
END IF

! This field is always required unless emissions are off or tracer updates
! from emissions and boundary layer mixing are suppressed.
! It is required anyway if nitrate emissions are produced (irrespective of
! whether emissions updates are suppressed) or if PSC heterogeneous chemistry
! is used with climatological surface area and CLASSIC SO4.
IF ((.NOT. (ukca_config%l_ukca_emissions_off .OR.                              &
           ukca_config%l_suppress_ems)) .OR.                                   &
    glomap_config%l_ukca_fine_no3_prod .OR.                                    &
    glomap_config%l_ukca_coarse_no3_prod .OR.                                  &
    (ukca_config%l_ukca_het_psc .AND. ukca_config%l_ukca_sa_clim .AND.         &
     ukca_config%l_use_classic_so4)) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_rho_r2
END IF

! Potential vorticity is always required unless the option to use an arbitrary
! fixed tropopause level is selected
IF (.NOT. ukca_config%l_fix_tropopause_level) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_pv_on_theta_mlevs
END IF

! Convective cloud amount and cloud fraction are required for photolysis
! if Fast-JX is used. This is temporary pending completion of the separation
! of photolysis from UKCA. (They are not used internally.)
IF (ukca_config%i_ukca_photol == i_ukca_fastjx) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_conv_cloud_amount
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_area_cloud_fraction
       !!!! WARNING This field may be updated internally (before Fast-JX call)
END IF

! Precipitation fields are required for wet deposition or for conservation of
! certain species if running a stratospheric scheme
IF ((.NOT. ukca_config%l_ukca_wetdep_off) .OR.                                 &
    ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_stratcfc .OR.             &
    ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat) THEN

  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_ls_rain3d
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_ls_snow3d

  ! Precipitation diagnosed by convection scheme are required if
  ! parameterized convection is used in the parent model
  IF (ukca_config%l_param_conv) THEN
    n = n + 1
    IF (n <= n_max) fld_names(n) = fldname_conv_rain3d
    n = n + 1
    IF (n <= n_max) fld_names(n) = fldname_conv_snow3d
  END IF

END IF

! SO4 surface aerosol climatology required if option to use
! climatological aerosol for surface area is selected
IF (ukca_config%l_ukca_sa_clim) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_so4_sa_clim
END IF

! SO4 in Aitken and accumulation modes are required if using CLASSIC sulphate
! aerosols from an external data source for heterogeneous chemistry
IF (ukca_config%l_use_classic_so4) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_so4_aitken
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_so4_accum
END IF

! CLASSIC aerosol fields are required if using option for heterogeneous
! chemistry on CLASSIC aerosols and the relevant species type is selected
IF (ukca_config%l_ukca_classic_hetchem .AND.                                   &
    ukca_config%l_use_classic_soot) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_soot_fresh
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_soot_aged
END IF
IF (ukca_config%l_ukca_classic_hetchem .AND.                                   &
    ukca_config%l_use_classic_ocff) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_ocff_fresh
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_ocff_aged
END IF
IF (ukca_config%l_ukca_classic_hetchem .AND.                                   &
    ukca_config%l_use_classic_biogenic) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_biogenic
END IF
IF (ukca_config%l_ukca_classic_hetchem .AND.                                   &
    ukca_config%l_use_classic_seasalt) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_sea_salt_film
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_sea_salt_jet
END IF

! CO2 field is required if option to use an external CO2 field is set
IF (ukca_config%l_chem_environ_co2_fld .AND. ANY(speci(:) == 'CO2       ')) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_co2_interactive
END IF

! Oxidant species are required for Offline chemistry schemes if they are
!  classed as 'CF' (constant field) in the species definition
IF (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) THEN
  DO i = 1, SIZE(ctype)
    IF (ctype(i) == 'CF') THEN
      SELECT CASE(speci(i))
      CASE ('HO2      ')
        n = n + 1
        IF (n <= n_max) fld_names(n) = fldname_ho2_offline
      CASE ('NO3      ')
        n = n + 1
        IF (n <= n_max) fld_names(n) = fldname_no3_offline
      CASE ('O3       ')
        n = n + 1
        IF (n <= n_max) fld_names(n) = fldname_o3_offline
      CASE ('OH       ')
        n = n + 1
        IF (n <= n_max) fld_names(n) = fldname_oh_offline
      CASE DEFAULT
        ! Do nothing
      END SELECT
    END IF
  END DO   ! Loop over ctype
  ! Special case for H2O2: this is a tracer but limited by external values
  IF (ANY(speci == 'H2O2     ')) THEN
    n = n + 1
    IF (n <= n_max) fld_names(n) = fldname_h2o2_offline
  END IF
  ! O3 field is also required for use in stratosphere if using a
  ! troposphere-only chemistry scheme unless overwrite of O3 & NO3 in the
  ! stratosphere is effectively switched off by setting the no. of levels
  ! above the tropopause at which overwrite occurs to be above the model domain
ELSE IF ((ukca_config%l_ukca_trop .OR. ukca_config%l_ukca_raq .OR.             &
          ukca_config%l_ukca_raqaero .OR. ukca_config%l_ukca_tropisop) .AND.   &
         ukca_config%nlev_above_trop_o3_env < ukca_config%model_levels) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_o3_offline
    IF (ukca_config%l_zon_av_ozone) fld_info(n)%ubound_dim1=1
  END IF
END IF

! Cloud liquid fraction is always required if GLOMAP-mode is on
IF (ukca_config%l_ukca_mode) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_cloud_liq_frac
END IF

! Autoconversion, accretion and riming rates are required if
! nucleation scavenging is on in GLOMAP-mode
IF (glomap_config%l_rainout) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_autoconv
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_accretion
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_rim_cry
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_rim_agg
END IF

! DUST 6 BINS NEEDED FOR NITRATE SCHEME
IF (glomap_config%l_6bin_dust_no3) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_dust_div1
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_dust_div2
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_dust_div3
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_dust_div4
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_dust_div5
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_dust_div6
END IF

! DUST 2 BINS NEEDED FOR NITRATE SCHEME
IF (glomap_config%l_2bin_dust_no3) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_dust_div1
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_dust_div2
END IF

! Grid box area on model_levels is required for certain MODE diagnostics.
! Note: no MODE diagnostics are currently supported independently of the UM
IF (ukca_config%l_enable_diag_um .AND. ukca_config%l_ukca_mode) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_grid_area_fullht
  END IF
END IF

! Grid box volume is required for various chemistry scheme diagnostics
! Note: a limited number of these diagnostics are currently supported
! independently of the UM. Of this subset, grid box volume is required for
! ASAD framework diagnostics.
IF (ukca_config%l_asad_chem_diags_support .OR.                                 &
    ukca_config%l_enable_diag_um) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_grid_volume
  END IF
END IF

! -- Environmental drivers in full-height plus zeroth level grid group --

! Altitudes of grid-cell interfaces are required if full support for vertical
! scaling of emissions is enabled or if MODE diagnostics are to be produced.
! Note: no MODE diagnostics are currently supported independently of the UM
IF (ukca_config%l_support_ems_vertprof .OR.                                    &
    (ukca_config%l_enable_diag_um .AND. ukca_config%l_ukca_mode)) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_interf_z
    fld_info(n)%group = group_fullht0_real
    fld_info(n)%lbound_dim3 = 0
  END IF
END IF

! -- Environmental drivers in full-height plus one grid group --

! Exner pressure on rho levels is always required unless emssions are off or
! tracer updates from emissions and boundary layer mixing are suppressed.
! It is required anyway if nitrate emissions are produced.
IF ((.NOT. (ukca_config%l_ukca_emissions_off .OR.                              &
           ukca_config%l_suppress_ems)) .OR.                                   &
    glomap_config%l_ukca_fine_no3_prod .OR.                                    &
    glomap_config%l_ukca_coarse_no3_prod) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_exner_rho_levels
    fld_info(n)%group = group_fullhtp1_real
    fld_info(n)%ubound_dim3 = ukca_config%model_levels + 1
  END IF
END IF

! -- Environmental drivers in boundary layer levels group --

! These fields are always required unless emissions are off or tracer updates
! from emissions and boundary layer mixing are suppressed
IF (.NOT. (ukca_config%l_ukca_emissions_off .OR.                               &
          ukca_config%l_suppress_ems)) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_rhokh_rdz
    fld_info(n)%group = group_bllev_real
    fld_info(n)%lbound_dim3 = 2
    fld_info(n)%ubound_dim3 = ukca_config%bl_levels
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_dtrdz
    fld_info(n)%group = group_bllev_real
    fld_info(n)%ubound_dim3 = ukca_config%bl_levels
  END IF
END IF

! Vertical component of wind speed and TKE are required for the Activate scheme
! Note: TKE is assumed to be unavailable at the top boundary layer level
! (as is the case in the UM parent model)
IF (glomap_config%i_ukca_activation_scheme == i_ukca_activation_arg) THEN
  n = n + 1
  IF (n <= n_max) fld_names(n) = fldname_vertvel
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_bl_tke
    fld_info(n)%group = group_bllev_real
    fld_info(n)%ubound_dim3 = ukca_config%bl_levels-1
  END IF
END IF

! -- Environmental drivers in entrainment levels group --

! These fields are always required unless emissions are off or tracer updates
! from emissions and boundary layer mixing are suppressed
IF (.NOT. (ukca_config%l_ukca_emissions_off .OR.                               &
          ukca_config%l_suppress_ems)) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_we_lim
    fld_info(n)%group = group_entlev_real
    fld_info(n)%ubound_dim3 = ukca_config%nlev_ent_tr_mix
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_t_frac
    fld_info(n)%group = group_entlev_real
    fld_info(n)%ubound_dim3 = ukca_config%nlev_ent_tr_mix
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_zrzi
    fld_info(n)%group = group_entlev_real
    fld_info(n)%ubound_dim3 = ukca_config%nlev_ent_tr_mix
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_we_lim_dsc
    fld_info(n)%group = group_entlev_real
    fld_info(n)%ubound_dim3 = ukca_config%nlev_ent_tr_mix
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_t_frac_dsc
    fld_info(n)%group = group_entlev_real
    fld_info(n)%ubound_dim3 = ukca_config%nlev_ent_tr_mix
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_zrzi_dsc
    fld_info(n)%group = group_entlev_real
    fld_info(n)%ubound_dim3 = ukca_config%nlev_ent_tr_mix
  END IF
END IF

! -- Environmental drivers in land point group --

! Soil moisture is required if using the interactive dry deposition scheme
IF (ukca_config%l_ukca_intdd) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_soil_moisture_layer1
    fld_info(n)%group = group_land_real
    fld_info(n)%ubound_dim1 = no_bound_value
    fld_info(n)%l_land_only = .TRUE.
  END IF
END IF

! Land fraction is required if using partial land points at coasts
IF (ukca_config%l_ctile) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_fland
    fld_info(n)%group = group_land_real
    fld_info(n)%ubound_dim1 = no_bound_value
    fld_info(n)%l_land_only = .TRUE.
  END IF
END IF

! Emissions from the land surface scheme ---- landpoints only -------
IF (.NOT. ukca_config%l_ukca_emissions_off) THEN

  DO i = 1, SIZE(em_chem_spec)
    l_req_emiss = .FALSE.
    use_fldname = ''

    SELECT CASE(TRIM(ADJUSTL(em_chem_spec(i))))
      ! Interactive biogenic emissions
    CASE ('C5H8')
      IF (ukca_config%l_ukca_ibvoc) THEN
        l_req_emiss = .TRUE.
        use_fldname = fldname_ibvoc_isoprene
      END IF

    CASE ('Monoterp')
      IF (ukca_config%l_ukca_ibvoc) THEN
        l_req_emiss = .TRUE.
        use_fldname = fldname_ibvoc_terpene
      END IF

    CASE ('MeOH', 'CH3OH')
      IF (ukca_config%l_ukca_ibvoc) THEN
        l_req_emiss = .TRUE.
        use_fldname = fldname_ibvoc_methanol
      END IF

    CASE ('Me2CO')
      IF (ukca_config%l_ukca_ibvoc) THEN
        l_req_emiss = .TRUE.
        use_fldname = fldname_ibvoc_acetone
      END IF

      !  Interactive Fire emissions
    CASE ('BC_biomass')
      IF (ukca_config%l_ukca_inferno) THEN
        l_req_emiss = .TRUE.
        use_fldname = fldname_inferno_bc
      END IF

    CASE ('CH4')
      IF ( ukca_config%l_ukca_inferno .AND. ukca_config%l_ukca_inferno_ch4     &
           .AND. .NOT. ukca_config%l_ukca_prescribech4 ) THEN
        l_req_emiss = .TRUE.
        use_fldname = fldname_inferno_ch4
      END IF

    CASE ('CO')
      IF (ukca_config%l_ukca_inferno) THEN
        l_req_emiss = .TRUE.
        use_fldname = fldname_inferno_co
      END IF

    CASE ('NO')
      IF (ukca_config%l_ukca_inferno) THEN
        l_req_emiss = .TRUE.
        use_fldname = fldname_inferno_nox
      END IF

    CASE ('OM_biomass')
      IF (ukca_config%l_ukca_inferno) THEN
        l_req_emiss = .TRUE.
        use_fldname = fldname_inferno_oc
      END IF

    CASE ('SO2_nat')
      IF (ukca_config%l_ukca_inferno) THEN
        l_req_emiss = .TRUE.
        use_fldname = fldname_inferno_so2
      END IF
    END SELECT

    ! Populate fld_name and fld_info if emission is required.
    IF ( l_req_emiss ) THEN
      n = n + 1
      IF (n <= n_max) THEN
        fld_names(n) = use_fldname
        fld_info(n)%group = group_land_real
        fld_info(n)%ubound_dim1 = no_bound_value
        fld_info(n)%l_land_only = .TRUE.
      END IF
    END IF  ! l_req_emiss

  END DO   ! loop over em_chem_spec

END IF   ! .NOT. ukca_config%l_ukca_emissions_off

! -- Environmental drivers in land-point tile real group --

! The following fields are required if using the interactive
! dry deposition scheme
IF (ukca_config%l_ukca_intdd) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_frac_types
    fld_info(n)%group = group_landtile_real
    fld_info(n)%ubound_dim1=no_bound_value
    fld_info(n)%ubound_dim2=ukca_config%ntype
    fld_info(n)%l_land_only = .TRUE.
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_tstar_tile
    fld_info(n)%group = group_landtile_real
    fld_info(n)%ubound_dim1=no_bound_value
    fld_info(n)%ubound_dim2=ukca_config%ntype
    fld_info(n)%l_land_only = .TRUE.
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_z0tile_lp
    fld_info(n)%group = group_landtile_real
    fld_info(n)%ubound_dim1=no_bound_value
    fld_info(n)%ubound_dim2=ukca_config%ntype
    fld_info(n)%l_land_only = .TRUE.
  END IF
END IF

! -- Environmental drivers in land-point tile logical group --

! The active tile indicator is required if using the interactive
! dry deposition scheme
IF (ukca_config%l_ukca_intdd) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_l_tile_active
    fld_info(n)%group = group_landtile_logical
    fld_info(n)%ubound_dim1=no_bound_value
    fld_info(n)%ubound_dim2=ukca_config%ntype
    fld_info(n)%l_land_only = .TRUE.
  END IF
END IF

! -- Environmental drivers in land-point plant functional type tile group --

! The following fields are required if using the interactive
! dry deposition scheme
IF (ukca_config%l_ukca_intdd) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_laift_lp
    fld_info(n)%group = group_landpft_real
    fld_info(n)%ubound_dim1=no_bound_value
    fld_info(n)%ubound_dim2=ukca_config%npft
    fld_info(n)%l_land_only = .TRUE.
  END IF
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_canhtft_lp
    fld_info(n)%group = group_landpft_real
    fld_info(n)%ubound_dim1=no_bound_value
    fld_info(n)%ubound_dim2=ukca_config%npft
    fld_info(n)%l_land_only = .TRUE.
  END IF
END IF

! Environment field to pass photolysis rates to UKCA
IF (ukca_config%i_ukca_photol /= i_ukca_nophot) THEN
  n = n + 1
  IF (n <= n_max) THEN
    fld_names(n) = fldname_photol_rates
    fld_info(n)%group = group_fullhtphot_real
    fld_info(n)%ubound_dim1 = ukca_config%row_length
    fld_info(n)%ubound_dim2 = ukca_config%rows
    fld_info(n)%ubound_dim3 = ukca_config%model_levels
    fld_info(n)%ubound_dim4 = SIZE(ratj_defs)
  END IF
END IF

! Check number of fields required against maximum
IF (n > n_max) THEN
  error_code = errcode_env_req_uninit
  IF (PRESENT(error_message)) WRITE(error_message,'(A,I0,A,I0)')               &
    'Number of required environment fields (', n,                              &
    ') exceeds maximum: n_max = ', n_max
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Create reference list of required fields with corresponding field info,
! field pointers and availability status flags.
ALLOCATE(environ_field_varnames(n))
environ_field_varnames = fld_names(1:n)
ALLOCATE(environ_field_info(n))
environ_field_info = fld_info(1:n)
ALLOCATE(environ_field_ptrs(n))
DO i = 1, n
  NULLIFY(environ_field_ptrs(i)%value_0d_real)
  NULLIFY(environ_field_ptrs(i)%value_1d_real)
  NULLIFY(environ_field_ptrs(i)%value_2d_integer)
  NULLIFY(environ_field_ptrs(i)%value_2d_real)
  NULLIFY(environ_field_ptrs(i)%value_2d_logical)
  NULLIFY(environ_field_ptrs(i)%value_3d_real)
  NULLIFY(environ_field_ptrs(i)%value_4d_real)
END DO
ALLOCATE(l_environ_field_available(n))
l_environ_field_available(:) = .FALSE.

! Set public flag to show availability of environment fields requirement
l_environ_req_available = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_environment_req


! ----------------------------------------------------------------------
SUBROUTINE ukca_get_environment_varlist(varnames_ptr, error_code,              &
                                        error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   UKCA API procedure that returns a list of field names identifying
!   the environment data required for the current UKCA configuration.
!
! Method:
!   Return pointer to the reference list giving the names of required
!   environment fields.
!   A non-zero error code is returned if the requirement for the current
!   UKCA configuration has not been initialised and the pointer will be
!   disassociated.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=maxlen_fieldname), POINTER, INTENT(OUT) :: varnames_ptr(:)
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_GET_ENVIRONMENT_VARLIST'

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check availability of environment fields requirement
IF (.NOT. l_environ_req_available) THEN
  error_code = errcode_env_req_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Environment fields requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  NULLIFY(varnames_ptr)
  RETURN
END IF

! Assign pointer to the reference list
varnames_ptr => environ_field_varnames

RETURN
END SUBROUTINE ukca_get_environment_varlist


! ----------------------------------------------------------------------
SUBROUTINE ukca_get_envgroup_varlists(error_code,                              &
                                      varnames_scalar_real_ptr,                &
                                      varnames_flat_integer_ptr,               &
                                      varnames_flat_real_ptr,                  &
                                      varnames_flat_logical_ptr,               &
                                      varnames_flatpft_real_ptr,               &
                                      varnames_fullht_real_ptr,                &
                                      varnames_fullht0_real_ptr,               &
                                      varnames_fullhtp1_real_ptr,              &
                                      varnames_bllev_real_ptr,                 &
                                      varnames_entlev_real_ptr,                &
                                      varnames_land_real_ptr,                  &
                                      varnames_landtile_real_ptr,              &
                                      varnames_landtile_logical_ptr,           &
                                      varnames_landpft_real_ptr,               &
                                      varnames_fullhtphot_real_ptr,            &
                                      error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   UKCA API procedure that returns lists of field names identifying
!   the environment data required for the current UKCA configuration by
!   group. Each group list is derived from the master list when first
!   requested.
!
! Method:
!   Return the list of names for each subgroup of the required
!   environment fields for which a pointer argument is present.
!   A non-zero error code is returned if the requirement for the current
!   UKCA configuration has not been initialised and any pointers passed
!   will be disassociated.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_scalar_real_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_flat_integer_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_flat_real_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_flat_logical_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_flatpft_real_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_fullht_real_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_fullht0_real_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_fullhtp1_real_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_bllev_real_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_entlev_real_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_land_real_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_landtile_real_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_landtile_logical_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_landpft_real_ptr(:)
CHARACTER(LEN=maxlen_fieldname), POINTER, OPTIONAL, INTENT(OUT) ::             &
  varnames_fullhtphot_real_ptr(:)
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables
! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_GET_ENVGROUP_VARLISTS'

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
  IF (PRESENT(varnames_scalar_real_ptr)) NULLIFY(varnames_scalar_real_ptr)
  IF (PRESENT(varnames_flat_integer_ptr)) NULLIFY(varnames_flat_integer_ptr)
  IF (PRESENT(varnames_flat_real_ptr)) NULLIFY(varnames_flat_real_ptr)
  IF (PRESENT(varnames_flat_logical_ptr)) NULLIFY(varnames_flat_logical_ptr)
  IF (PRESENT(varnames_flatpft_real_ptr)) NULLIFY(varnames_flatpft_real_ptr)
  IF (PRESENT(varnames_fullht_real_ptr)) NULLIFY(varnames_fullht_real_ptr)
  IF (PRESENT(varnames_fullht0_real_ptr)) NULLIFY(varnames_fullht0_real_ptr)
  IF (PRESENT(varnames_fullhtp1_real_ptr)) NULLIFY(varnames_fullhtp1_real_ptr)
  IF (PRESENT(varnames_bllev_real_ptr)) NULLIFY(varnames_bllev_real_ptr)
  IF (PRESENT(varnames_entlev_real_ptr)) NULLIFY(varnames_entlev_real_ptr)
  IF (PRESENT(varnames_land_real_ptr)) NULLIFY(varnames_land_real_ptr)
  IF (PRESENT(varnames_landtile_real_ptr)) NULLIFY(varnames_landtile_real_ptr)
  IF (PRESENT(varnames_landtile_logical_ptr))                                  &
    NULLIFY(varnames_landtile_logical_ptr)
  IF (PRESENT(varnames_landpft_real_ptr)) NULLIFY(varnames_landpft_real_ptr)
  IF (PRESENT(varnames_fullhtphot_real_ptr))                                   &
    NULLIFY(varnames_fullhtphot_real_ptr)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Assign pointers to subgroup lists for any subgroup pointer arguments present.
! Subgroup lists are set up first if not already available.

IF (PRESENT(varnames_scalar_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_scalar_real)) THEN
    CALL setup_envgroup_varlist(group_scalar_real,                             &
                                environ_field_varnames_scalar_real)
  END IF
  varnames_scalar_real_ptr => environ_field_varnames_scalar_real
END IF

IF (PRESENT(varnames_flat_integer_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_flat_integer)) THEN
    CALL setup_envgroup_varlist(group_flat_integer,                            &
                                environ_field_varnames_flat_integer)
  END IF
  varnames_flat_integer_ptr => environ_field_varnames_flat_integer
END IF

IF (PRESENT(varnames_flat_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_flat_real)) THEN
    CALL setup_envgroup_varlist(group_flat_real,                               &
                                environ_field_varnames_flat_real)
  END IF
  varnames_flat_real_ptr => environ_field_varnames_flat_real
END IF

IF (PRESENT(varnames_flat_logical_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_flat_logical)) THEN
    CALL setup_envgroup_varlist(group_flat_logical,                            &
                                environ_field_varnames_flat_logical)
  END IF
  varnames_flat_logical_ptr => environ_field_varnames_flat_logical
END IF

IF (PRESENT(varnames_flatpft_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_flatpft_real)) THEN
    CALL setup_envgroup_varlist(group_flatpft_real,                            &
                                environ_field_varnames_flatpft_real)
  END IF
  varnames_flatpft_real_ptr => environ_field_varnames_flatpft_real
END IF

IF (PRESENT(varnames_fullht_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_fullht_real)) THEN
    CALL setup_envgroup_varlist(group_fullht_real,                             &
                                environ_field_varnames_fullht_real)
  END IF
  varnames_fullht_real_ptr => environ_field_varnames_fullht_real
END IF

IF (PRESENT(varnames_fullht0_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_fullht0_real)) THEN
    CALL setup_envgroup_varlist(group_fullht0_real,                            &
                                environ_field_varnames_fullht0_real)
  END IF
  varnames_fullht0_real_ptr => environ_field_varnames_fullht0_real
END IF

IF (PRESENT(varnames_fullhtp1_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_fullhtp1_real)) THEN
    CALL setup_envgroup_varlist(group_fullhtp1_real,                           &
                                environ_field_varnames_fullhtp1_real)
  END IF
  varnames_fullhtp1_real_ptr => environ_field_varnames_fullhtp1_real
END IF

IF (PRESENT(varnames_bllev_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_bllev_real)) THEN
    CALL setup_envgroup_varlist(group_bllev_real,                              &
                                environ_field_varnames_bllev_real)
  END IF
  varnames_bllev_real_ptr => environ_field_varnames_bllev_real
END IF

IF (PRESENT(varnames_entlev_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_entlev_real)) THEN
    CALL setup_envgroup_varlist(group_entlev_real,                             &
                                environ_field_varnames_entlev_real)
  END IF
  varnames_entlev_real_ptr => environ_field_varnames_entlev_real
END IF


IF (PRESENT(varnames_land_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_land_real)) THEN
    CALL setup_envgroup_varlist(group_land_real,                               &
                                environ_field_varnames_land_real)
  END IF
  varnames_land_real_ptr => environ_field_varnames_land_real
END IF

IF (PRESENT(varnames_landtile_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_landtile_real)) THEN
    CALL setup_envgroup_varlist(group_landtile_real,                           &
                                environ_field_varnames_landtile_real)
  END IF
  varnames_landtile_real_ptr => environ_field_varnames_landtile_real
END IF

IF (PRESENT(varnames_landtile_logical_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_landtile_logical)) THEN
    CALL setup_envgroup_varlist(group_landtile_logical,                        &
                                environ_field_varnames_landtile_logical)
  END IF
  varnames_landtile_logical_ptr => environ_field_varnames_landtile_logical
END IF

IF (PRESENT(varnames_landpft_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_landpft_real)) THEN
    CALL setup_envgroup_varlist(group_landpft_real,                            &
                                environ_field_varnames_landpft_real)
  END IF
  varnames_landpft_real_ptr => environ_field_varnames_landpft_real
END IF

IF (PRESENT(varnames_fullhtphot_real_ptr)) THEN
  IF (.NOT. ALLOCATED(environ_field_varnames_fullhtphot_real)) THEN
    CALL setup_envgroup_varlist(group_fullhtphot_real,                         &
                                environ_field_varnames_fullhtphot_real)
  END IF
  varnames_fullhtphot_real_ptr => environ_field_varnames_fullhtphot_real
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_get_envgroup_varlists


! ----------------------------------------------------------------------
SUBROUTINE setup_envgroup_varlist(group,varnames)
! ----------------------------------------------------------------------
! Description:
!   Returns a list of the field names for the specified subgroup of
!   environment fields
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: group
CHARACTER(LEN=maxlen_fieldname), ALLOCATABLE, INTENT(OUT) :: varnames(:)

! Local variables
INTEGER :: n_req
INTEGER :: n
INTEGER :: i

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SETUP_ENVGROUP_VARLIST'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

n_req = SIZE(environ_field_varnames)

! Allocate space for the list of names
n = 0
DO i = 1, n_req
  IF (environ_field_info(i)%group == group) n = n + 1
END DO
ALLOCATE(varnames(n))

! Populate the list from the master list of required fields
n = 0
DO i = 1, n_req
  IF (environ_field_info(i)%group == group) THEN
    n = n + 1
    varnames(n) = environ_field_varnames(i)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE setup_envgroup_varlist


! ----------------------------------------------------------------------
INTEGER FUNCTION environ_field_index(varname)
! ----------------------------------------------------------------------
! Description:
!   Returns the index of a variable name in the required environment
!   field array or 0 if it is not found (or if the environment fields
!   requirement has not been initialised).
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname

! Local variables

INTEGER :: i
LOGICAL :: found

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ENVIRON_FIELD_INDEX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

environ_field_index = 0
IF (l_environ_req_available) THEN
  found = .FALSE.
  i = 0
  DO WHILE (i < SIZE(environ_field_varnames) .AND. (.NOT. found))
    i = i + 1
    found = (varname == environ_field_varnames(i))
  END DO
  IF (found) environ_field_index = i
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION environ_field_index


! ----------------------------------------------------------------------
SUBROUTINE check_environment_availability(n_fld_present, n_fld_missing,        &
                                          availability_ptr)
! ----------------------------------------------------------------------
! Description:
!   Checks availability of required UKCA environment fields.
!
! Method:
!   Count number of required fields present and number of required fields
!   missing as indicated by the field availability flags.
!   If an availability pointer is provided as an optional argument,
!   assign this to the availability flag array to give access to the
!   availability status of each required field.
!   If the requirement for the current UKCA configuration has not been
!   initialised, the field count will be zero and the availability
!   pointer (if present) will be disassociated.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(OUT) :: n_fld_present
INTEGER, INTENT(OUT) :: n_fld_missing
LOGICAL, POINTER, OPTIONAL, INTENT(OUT) :: availability_ptr(:)

! Local variables

INTEGER :: i

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CHECK_ENVIRONMENT_AVAILABILITY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_environ_req_available) THEN

  ! Check availability of required fields
  n_fld_present = 0
  n_fld_missing = 0
  DO i = 1, SIZE(l_environ_field_available)
    IF (l_environ_field_available(i)) THEN
      n_fld_present = n_fld_present + 1
    ELSE
      n_fld_missing = n_fld_missing + 1
    END IF
  END DO

  ! Set pointer argument if present to return supporting info
  IF (PRESENT(availability_ptr))                                               &
    availability_ptr => l_environ_field_available

ELSE
  n_fld_present = 0
  n_fld_missing = 0
  IF (PRESENT(availability_ptr)) NULLIFY(availability_ptr)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_environment_availability


! ----------------------------------------------------------------------
LOGICAL FUNCTION environ_field_available(varname)
! ----------------------------------------------------------------------
! Description:
!   Returns the availability of a variable name in the required environment
!   field array. T if available, F if not available (or if not required or
!   if the environment fields requirement has not been initialised).
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: varname

! Local variables
INTEGER :: i

! Dr hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ENVIRON_FIELD_AVAILABLE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

environ_field_available = .FALSE.

! Proceed only if environment fields requirement has been set
IF ( l_environ_req_available ) THEN
  fields_loop: DO i = 1, SIZE(environ_field_varnames)
    IF (varname == environ_field_varnames(i)) THEN
      environ_field_available = l_environ_field_available(i)
      EXIT fields_loop
    END IF
  END DO fields_loop
  ! If we are here the supplied varname does not match any required env field
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION environ_field_available


! ----------------------------------------------------------------------
SUBROUTINE print_environment_summary(error_code_ptr, id_code, environ_ptrs,    &
                                     error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Writes summary info on UKCA environment field values to the log file.
!
! Method:
!
!   Summary data are printed for each environmental driver field in the
!   form of four tables:
!   Table 1 gives max, min and mean for each field in 0D, 1D, or 2D groups
!   and for 3D fields without a vertical dimension.
!   Table 2 gives max, min and mean by level for each 3D field with a
!   vertical dimension.
!   Table 3 gives max, min and mean by reaction for the photolysis rate
!   4D field.
!   Table 4 gives the count of true and false values for logical fields
!
!   An id value passed as an input argument is included as the first
!   entry in each record tabluated. This can be set by the calling
!   routine to indicate, for example, which processor the data are from
!   when running in a distributed memory system.
!
!   A character is printed at the beginning of each table entry to
!   prevent any automatic left justification that might be applied
!   on output from affecting the column alignment.
! ----------------------------------------------------------------------

USE ukca_chem_defs_mod, ONLY: ratj_defs
USE umPrintMgr, ONLY: umMessage, UmPrint

IMPLICIT NONE

! Subroutine arguments

INTEGER, POINTER, INTENT(IN) :: error_code_ptr
INTEGER, INTENT(IN) :: id_code

TYPE(env_field_ptrs_type), INTENT(IN) :: environ_ptrs(:)

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER :: i_jrat                                  ! Photolysis reaction index
INTEGER :: n_true                                  ! No. of .TRUE. values in
                                                   ! a field array
INTEGER :: level_min                               ! Min level of all 3D fields
INTEGER :: level_max                               ! Max level of all 3D fields
INTEGER :: k1                                      ! 3D field - lowest level
INTEGER :: k2                                      ! 3D field - highest level
INTEGER :: i, k                                    ! Loop counters

CHARACTER(LEN=1), PARAMETER :: l_border = '*'      ! Character used for left
                                                   ! border of table
CHARACTER(LEN=43) :: products_txt                  ! Buffer for names of up to 4
                                                   ! photolysis products plus
                                                   ! separators

LOGICAL :: l_print_logicals_table                  ! T if table of logical
                                                   ! fields summary table
                                                   ! needs printing
LOGICAL :: l_print_photol_rates_table              ! T if photol_rates summary
                                                   ! table needs printing

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'PRINT_ENVIRONMENT_SUMMARY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check size of field pointers array is as expected
IF (SIZE(environ_ptrs) /= SIZE(environ_field_varnames)) THEN
  error_code_ptr = errcode_env_field_mismatch
  IF (PRESENT(error_message))                                                  &
    WRITE(error_message, '(A)')                                                &
      'Size of environ_ptrs does not match no. of required environment fields'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook)                                                                   &
    CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF


! Write column headers for main table
WRITE(umMessage,'(A,A5,1X,A3,1X,A20,3(1X,A11))')                               &
      l_border, 'ID   ', 'GRP', 'FIELD_NAME          ',                        &
      'MAXIMUM    ', 'MINIMUM    ', 'MEAN       '
CALL umPrint(umMessage,src=RoutineName)

! Write data summary records for each 0D, 1D and 2D driver field except
! logicals. Also include 3D fields without vertical dimension.
! For logicals and 4D fields, flag presence for later tabulation.
! For 3D fields with a vertical dimension, determine minimum lower and maximum
! upper levels for later tabulation.

level_min = HUGE(0)
level_max = -HUGE(0)
l_print_logicals_table = .FALSE.
l_print_photol_rates_table = .FALSE.

DO i = 1, SIZE(environ_field_varnames)

  SELECT CASE (environ_field_info(i)%group)
  CASE (group_scalar_real)
    ! 0D field
    IF (ASSOCIATED(environ_ptrs(i)%value_0d_real)) THEN
      WRITE(umMessage,'(A,I5,1X,I3,1X,A20,3(1X,E11.4))')                       &
            l_border, id_code,                                                 &
            environ_field_info(i)%group,                                       &
            environ_field_varnames(i),                                         &
            environ_ptrs(i)%value_0d_real,                                     &
            environ_ptrs(i)%value_0d_real,                                     &
            environ_ptrs(i)%value_0d_real
      CALL umPrint(umMessage,src=RoutineName)
    ELSE
      error_code_ptr = errcode_ukca_internal_fault
    END IF
  CASE (group_land_real)
    ! 1D field
    IF (ASSOCIATED(environ_ptrs(i)%value_1d_real)) THEN
      WRITE(umMessage,'(A,I5,1X,I3,1X,A20,3(1X,E11.4))')                       &
            l_border, id_code,                                                 &
            environ_field_info(i)%group,                                       &
            environ_field_varnames(i),                                         &
            MAXVAL(environ_ptrs(i)%value_1d_real),                             &
            MINVAL(environ_ptrs(i)%value_1d_real),                             &
            SUM(environ_ptrs(i)%value_1d_real) /                               &
            SIZE(environ_ptrs(i)%value_1d_real)
      CALL umPrint(umMessage,src=RoutineName)
    ELSE
      error_code_ptr = errcode_ukca_internal_fault
    END IF
  CASE (group_flat_integer)
    ! 2D integer field
    IF (ASSOCIATED(environ_ptrs(i)%value_2d_integer)) THEN
      WRITE(umMessage,'(A,I5,1X,I3,1X,A20,2(1X,I11),1X,E11.4)')                &
            l_border, id_code,                                                 &
            environ_field_info(i)%group,                                       &
            environ_field_varnames(i),                                         &
            MAXVAL(environ_ptrs(i)%value_2d_integer),                          &
            MINVAL(environ_ptrs(i)%value_2d_integer),                          &
            REAL(SUM(environ_ptrs(i)%value_2d_integer)) /                      &
            SIZE(environ_ptrs(i)%value_2d_integer)
      CALL umPrint(umMessage,src=RoutineName)
    ELSE
      error_code_ptr = errcode_ukca_internal_fault
    END IF
  CASE (group_flat_real, group_landtile_real, group_landpft_real)
    ! 2D real field
    IF (ASSOCIATED(environ_ptrs(i)%value_2d_real)) THEN
      WRITE(umMessage,'(A,I5,1X,I3,1X,A20,3(1X,E11.4))')                       &
            l_border, id_code,                                                 &
            environ_field_info(i)%group,                                       &
            environ_field_varnames(i),                                         &
            MAXVAL(environ_ptrs(i)%value_2d_real),                             &
            MINVAL(environ_ptrs(i)%value_2d_real),                             &
            SUM(environ_ptrs(i)%value_2d_real) /                               &
            SIZE(environ_ptrs(i)%value_2d_real)
      CALL umPrint(umMessage,src=RoutineName)
    ELSE
      error_code_ptr = errcode_ukca_internal_fault
    END IF
  CASE (group_flat_logical)
    ! 2D logical field (print separately later)
    l_print_logicals_table = .TRUE.
  CASE (group_flatpft_real)
    ! Flat 3D real field, defined on PFTs
    IF (ASSOCIATED(environ_ptrs(i)%value_3d_real)) THEN
      WRITE(umMessage,'(A,I5,1X,I3,1X,A20,3(1X,E11.4))')                       &
            l_border, id_code,                                                 &
            environ_field_info(i)%group,                                       &
            environ_field_varnames(i),                                         &
            MAXVAL(environ_ptrs(i)%value_3d_real),                             &
            MINVAL(environ_ptrs(i)%value_3d_real),                             &
            SUM(environ_ptrs(i)%value_3d_real) /                               &
            SIZE(environ_ptrs(i)%value_3d_real)
      CALL umPrint(umMessage,src=RoutineName)
    ELSE
      error_code_ptr = errcode_ukca_internal_fault
    END IF
  CASE (group_fullht_real, group_fullht0_real, group_fullhtp1_real,            &
        group_bllev_real, group_entlev_real)
    ! 3D real field (print summary by level)
    IF (ASSOCIATED(environ_ptrs(i)%value_3d_real)) THEN
      k1 = LBOUND(environ_ptrs(i)%value_3d_real, DIM=3)
      k2 = UBOUND(environ_ptrs(i)%value_3d_real, DIM=3)
      level_min = MIN(level_min, k1)
      level_max = MAX(level_max, k2)
    ELSE
      error_code_ptr = errcode_ukca_internal_fault
    END IF
  CASE (group_fullhtphot_real)
    ! 4D real field - photol_rates (print separately later)
    l_print_photol_rates_table = .TRUE.
  END SELECT

  IF (error_code_ptr > 0) THEN
    IF (PRESENT(error_message))                                                &
      WRITE(error_message, '(A,A)')                                            &
        'Missing pointer for environmental driver field ',                     &
        environ_field_varnames(i)
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook)                                                                 &
      CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

END DO  ! i

! Write separate table for 3D fields, ordered by level
! (Note: pointer association has already been checked above, no need to repeat)

! Write column headers
WRITE(umMessage,'(A,A5,1X,A3,1X,A5,1X,A20,3(1X,A11))')                         &
      l_border, 'ID   ', 'GRP', 'LEVEL', 'FIELD_NAME          ',               &
      'MAXIMUM    ', 'MINIMUM    ', 'MEAN       '
CALL umPrint(umMessage,src=RoutineName)

DO k = level_min, level_max
  DO i = 1, SIZE(environ_field_varnames)
    SELECT CASE(environ_field_info(i)%group)
    CASE (group_fullht_real, group_fullht0_real, group_fullhtp1_real,          &
          group_bllev_real, group_entlev_real)
      IF (k >= LBOUND(environ_ptrs(i)%value_3d_real, DIM=3) .AND.              &
          k <= UBOUND(environ_ptrs(i)%value_3d_real, DIM=3)) THEN
        WRITE(umMessage,'(A,I5,1X,I3,1X,I5,1X,A20,3(1X,E11.4))')               &
              l_border, id_code,                                               &
              environ_field_info(i)%group, k,                                  &
              environ_field_varnames(i),                                       &
              MAXVAL(environ_ptrs(i)%value_3d_real(:,:,k)),                    &
              MINVAL(environ_ptrs(i)%value_3d_real(:,:,k)),                    &
              SUM(environ_ptrs(i)%value_3d_real(:,:,k)) /                      &
              SIZE(environ_ptrs(i)%value_3d_real(:,:,k))
        CALL umPrint(umMessage,src=RoutineName)
      END IF
    END SELECT
  END DO  ! i
END DO  ! k

! Write separate table for photolysis rates
! 4D real field - photol_rates (print summary by reaction, combined levels)

IF (l_print_photol_rates_table) THEN

  ! Write column headers
  WRITE(umMessage,'(A,A5,1X,A3,1X,A20,1X,A5,3(1X,A11),2(1X,A10),1X,A8)')       &
        l_border, 'ID   ', 'GRP', 'FIELD_NAME          ', 'INDEX',             &
        'MAXIMUM    ', 'MINIMUM    ', 'MEAN       ', 'RATE_LABEL',             &
        'REACTANT  ', 'PRODUCTS'
  CALL umPrint(umMessage,src=RoutineName)

  DO i = 1, SIZE(environ_field_varnames)

    IF (environ_field_info(i)%group == group_fullhtphot_real) THEN

      IF (ASSOCIATED(environ_ptrs(i)%value_4d_real)) THEN

        ! Write data summary records
        DO i_jrat = 1, SIZE(environ_ptrs(i)%value_4d_real, DIM=4)

          ! Collate products of reaction
          IF (ratj_defs(i_jrat)%prod4 /= '') THEN
            products_txt = TRIM(ratj_defs(i_jrat)%prod1) // ',' //             &
                           TRIM(ratj_defs(i_jrat)%prod2) // ',' //             &
                           TRIM(ratj_defs(i_jrat)%prod3) // ',' //             &
                           TRIM(ratj_defs(i_jrat)%prod4)
          ELSE IF (ratj_defs(i_jrat)%prod3 /= '') THEN
            products_txt = TRIM(ratj_defs(i_jrat)%prod1) // ',' //             &
                           TRIM(ratj_defs(i_jrat)%prod2) // ',' //             &
                           TRIM(ratj_defs(i_jrat)%prod3)
          ELSE IF (ratj_defs(i_jrat)%prod2 /= '') THEN
            products_txt = TRIM(ratj_defs(i_jrat)%prod1) // ',' //             &
                           TRIM(ratj_defs(i_jrat)%prod2)
          ELSE
            products_txt = TRIM(ratj_defs(i_jrat)%prod1)
          END IF

          WRITE(umMessage,                                                     &
                '(A,I5,1X,I3,1X,A20,1X,I5,3(1X,E11.4),2(1X,A10),1X,A43)')      &
                l_border, id_code,                                             &
                environ_field_info(i)%group,                                   &
                environ_field_varnames(i), i_jrat,                             &
                MAXVAL(environ_ptrs(i)%value_4d_real(:,:,:,i_jrat)),           &
                MINVAL(environ_ptrs(i)%value_4d_real(:,:,:,i_jrat)),           &
                SUM(environ_ptrs(i)%value_4d_real(:,:,:,i_jrat)) /             &
                SIZE(environ_ptrs(i)%value_4d_real(:,:,:,i_jrat)),             &
                ratj_defs(i_jrat)%fname,                                       &
                ratj_defs(i_jrat)%react1,                                      &
                products_txt
          CALL umPrint(umMessage,src=RoutineName)

        END DO  ! i_jrat

      ELSE
        error_code_ptr = errcode_ukca_internal_fault
      END IF

    END IF

    IF (error_code_ptr /= 0) THEN
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,A)')                                          &
          'Missing pointer for environmental driver field ',                   &
          environ_field_varnames(i)
      IF (PRESENT(error_routine)) error_routine = RoutineName
      IF (lhook)                                                               &
        CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF

  END DO  ! i

END IF

! Write separate table for logical fields
! 2D logical fields (print count of .TRUE. and .FALSE. values)

IF (l_print_logicals_table) THEN

  ! Write column headers
  WRITE(umMessage,'(A,A5,1X,A3,1X,A20,1X,2(1X,A7))')                           &
        l_border, 'ID   ', 'GRP', 'FIELD_NAME          ', 'N_TRUE ', 'N_FALSE'
  CALL umPrint(umMessage,src=RoutineName)

  DO i = 1, SIZE(environ_field_varnames)

    IF (environ_field_info(i)%group == group_flat_logical) THEN
      IF (ASSOCIATED(environ_ptrs(i)%value_2d_logical)) THEN
        n_true = COUNT(environ_ptrs(i)%value_2d_logical)
        WRITE(umMessage,'(A,I5,1X,I3,1X,A20,1X,2(1X,I7))')                     &
              l_border, id_code,                                               &
              environ_field_info(i)%group,                                     &
              environ_field_varnames(i), n_true,                               &
              SIZE(environ_ptrs(i)%value_2d_logical) - n_true
        CALL umPrint(umMessage,src=RoutineName)
      ELSE
        error_code_ptr = errcode_ukca_internal_fault
      END IF
    END IF

    IF (error_code_ptr /= 0) THEN
      IF (PRESENT(error_message))                                              &
        WRITE(error_message, '(A,A)')                                          &
          'Missing pointer for environmental driver field ',                   &
          environ_field_varnames(i)
      IF (PRESENT(error_routine)) error_routine = RoutineName
      IF (lhook)                                                               &
        CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF

  END DO  ! i

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE print_environment_summary


! ----------------------------------------------------------------------
SUBROUTINE clear_environment_req()
! ----------------------------------------------------------------------
! Description:
!   Resets all data relating to environment data requirement to its
!   initial state for a new UKCA configuration.
!
! Method:
!   Deallocate required field list arrays and reset the flag showing
!   availability status of the environment data requirement.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Local variables
! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CLEAR_ENVIRONMENT_REQ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(environ_field_varnames)) DEALLOCATE(environ_field_varnames)

IF (ALLOCATED(environ_field_varnames_scalar_real))                             &
  DEALLOCATE(environ_field_varnames_scalar_real)
IF (ALLOCATED(environ_field_varnames_flat_integer))                            &
  DEALLOCATE(environ_field_varnames_flat_integer)
IF (ALLOCATED(environ_field_varnames_flat_real))                               &
  DEALLOCATE(environ_field_varnames_flat_real)
IF (ALLOCATED(environ_field_varnames_flat_logical))                            &
  DEALLOCATE(environ_field_varnames_flat_logical)
IF (ALLOCATED(environ_field_varnames_flatpft_real))                            &
  DEALLOCATE(environ_field_varnames_flatpft_real)
IF (ALLOCATED(environ_field_varnames_fullht_real))                             &
  DEALLOCATE(environ_field_varnames_fullht_real)
IF (ALLOCATED(environ_field_varnames_fullht0_real))                            &
  DEALLOCATE(environ_field_varnames_fullht0_real)
IF (ALLOCATED(environ_field_varnames_fullhtp1_real))                           &
  DEALLOCATE(environ_field_varnames_fullhtp1_real)
IF (ALLOCATED(environ_field_varnames_bllev_real))                              &
  DEALLOCATE(environ_field_varnames_bllev_real)
IF (ALLOCATED(environ_field_varnames_entlev_real))                             &
  DEALLOCATE(environ_field_varnames_entlev_real)
IF (ALLOCATED(environ_field_varnames_land_real))                               &
  DEALLOCATE(environ_field_varnames_land_real)
IF (ALLOCATED(environ_field_varnames_landtile_real))                           &
  DEALLOCATE(environ_field_varnames_landtile_real)
IF (ALLOCATED(environ_field_varnames_landtile_logical))                        &
  DEALLOCATE(environ_field_varnames_landtile_logical)
IF (ALLOCATED(environ_field_varnames_landpft_real))                            &
  DEALLOCATE(environ_field_varnames_landpft_real)
IF (ALLOCATED(environ_field_varnames_fullhtphot_real))                         &
  DEALLOCATE(environ_field_varnames_fullhtphot_real)

IF (ALLOCATED(environ_field_info)) DEALLOCATE(environ_field_info)
IF (ALLOCATED(environ_field_ptrs)) DEALLOCATE(environ_field_ptrs)
IF (ALLOCATED(l_environ_field_available)) DEALLOCATE(l_environ_field_available)

l_environ_req_available = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE clear_environment_req

! ----------------------------------------------------------------------

END MODULE ukca_environment_req_mod
