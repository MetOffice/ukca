!*****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module containing subroutine ukca_setup for checking input
!   configuration options and setting up internal UKCA data to
!   define the configuration.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_setup_mod

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_SETUP_MOD'

! Public procedures
PUBLIC :: ukca_setup

CONTAINS

! ----------------------------------------------------------------------

SUBROUTINE ukca_setup(error_code,                                              &
                      row_length,                                              &
                      rows,                                                    &
                      model_levels,                                            &
                      bl_levels,                                               &
                      nlev_ent_tr_mix,                                         &
                      ntype,                                                   &
                      npft,                                                    &
                      i_brd_leaf,                                              &
                      i_brd_leaf_dec,                                          &
                      i_brd_leaf_eg_trop,                                      &
                      i_brd_leaf_eg_temp,                                      &
                      i_ndl_leaf,                                              &
                      i_ndl_leaf_dec,                                          &
                      i_ndl_leaf_eg,                                           &
                      i_c3_grass,                                              &
                      i_c3_crop,                                               &
                      i_c3_pasture,                                            &
                      i_c4_grass,                                              &
                      i_c4_crop,                                               &
                      i_c4_pasture,                                            &
                      i_shrub,                                                 &
                      i_shrub_dec,                                             &
                      i_shrub_eg,                                              &
                      i_urban,                                                 &
                      i_lake,                                                  &
                      i_soil,                                                  &
                      i_ice,                                                   &
                      i_elev_ice,                                              &
                      i_ukca_chem,                                             &
                      fixed_tropopause_level,                                  &
                      i_ageair_reset_method,                                   &
                      max_ageair_reset_level,                                  &
                      i_error_method,                                          &
                      i_ukca_chem_version,                                     &
                      nrsteps,                                                 &
                      chem_timestep,                                           &
                      dts0,                                                    &
                      nit,                                                     &
                      i_ukca_quasinewton_start,                                &
                      i_ukca_quasinewton_end,                                  &
                      ukca_chem_seg_size,                                      &
                      nlev_above_trop_o3_env,                                  &
                      nlev_ch4_stratloss,                                      &
                      i_ukca_topboundary,                                      &
                      i_ukca_hetconfig,                                        &
                      i_inferno_emi,                                           &
                      i_ukca_dms_flux,                                         &
                      i_ukca_light_param,                                      &
                      i_strat_lbc_source,                                      &
                      i_primss_method,                                         &
                      env_log_step,                                            &
                      i_mode_nzts,                                             &
                      ukca_mode_seg_size,                                      &
                      i_mode_setup,                                            &
                      i_mode_bln_param_method,                                 &
                      i_mode_nucscav,                                          &
                      i_ukca_activation_scheme,                                &
                      i_ukca_nwbins,                                           &
                      i_ukca_tune_bc,                                          &
                      i_photol_scheme,                                         &
                      i_photol_scheme_off,                                     &
                      i_photol_scheme_strat_only,                              &
                      i_photol_scheme_2d,                                      &
                      i_photol_scheme_fastjx,                                  &
                      dzsoil_layer1,                                           &
                      timestep,                                                &
                      max_ageair_reset_height,                                 &
                      max_z_for_offline_chem,                                  &
                      soa_yield_scaling_mt,                                    &
                      soa_yield_scaling_isop,                                  &
                      dry_depvel_so2_scaling,                                  &
                      fastjx_prescutoff,                                       &
                      mode_parfrac,                                            &
                      seadms_ems_scaling,                                      &
                      sea_salt_ems_scaling,                                    &
                      marine_pom_ems_scaling,                                  &
                      lightnox_scale_fac,                                      &
                      anth_so2_ems_scaling,                                    &
                      mode_activation_dryr,                                    &
                      dry_depvel_acc_scaling,                                  &
                      acc_cor_scav_scaling,                                    &
                      mode_incld_so2_rfrac,                                    &
                      biom_aer_ems_scaling,                                    &
                      hno3_uptake_coeff,                                       &
                      ph_fit_coeff_a,                                          &
                      ph_fit_coeff_b,                                          &
                      ph_fit_intercept,                                        &
                      sigwmin,                                                 &
                      sigma_updraught_scaling,                                 &
                      const_rmol,                                              &
                      const_tfs,                                               &
                      const_rho_water,                                         &
                      const_rhosea,                                            &
                      const_lc,                                                &
                      const_avogadro,                                          &
                      const_boltzmann,                                         &
                      const_rho_so4,                                           &
                      l_cal360,                                                &
                      l_ukca_chem_aero,                                        &
                      l_ukca_mode,                                             &
                      l_fix_tropopause_level,                                  &
                      l_ukca_ageair,                                           &
                      l_blankout_invalid_diags,                                &
                      l_enable_diag_um,                                        &
                      l_ukca_persist_off,                                      &
                      l_timer,                                                 &
                      l_ukca_emissions_off,                                    &
                      l_ukca_drydep_off,                                       &
                      l_ukca_wetdep_off,                                       &
                      l_ukca_scale_ppe,                                        &
                      l_ukca_asad_columns,                                     &
                      l_ukca_asad_full,                                        &
                      l_ukca_debug_asad,                                       &
                      l_ukca_intdd,                                            &
                      l_ukca_ddepo3_ocean,                                     &
                      l_ukca_ddep_lev1,                                        &
                      l_ukca_dry_dep_so2wet,                                   &
                      l_deposition_jules,                                      &
                      l_ukca_quasinewton,                                      &
                      l_tracer_lumping,                                        &
                      l_ukca_ro2_ntp,                                          &
                      l_ukca_ro2_perm,                                         &
                      l_ukca_intph,                                            &
                      l_ukca_scale_soa_yield_mt,                               &
                      l_ukca_scale_soa_yield_isop,                             &
                      l_ukca_het_psc,                                          &
                      l_ukca_limit_nat,                                        &
                      l_ukca_sa_clim,                                          &
                      l_ukca_trophet,                                          &
                      l_ukca_classic_hetchem,                                  &
                      l_use_photolysis,                                        &
                      l_ukca_ibvoc,                                            &
                      l_ukca_inferno,                                          &
                      l_ukca_inferno_ch4,                                      &
                      l_ukca_so2ems_expvolc,                                   &
                      l_ukca_so2ems_plumeria,                                  &
                      l_ukca_qch4inter,                                        &
                      l_ukca_emsdrvn_ch4,                                      &
                      l_ukca_enable_seadms_ems,                                &
                      l_ukca_scale_seadms_ems,                                 &
                      l_ukca_linox_scaling,                                    &
                      l_support_ems_vertprof,                                  &
                      l_support_ems_gridbox_units,                             &
                      l_suppress_ems,                                          &
                      l_ukca_h2o_feedback,                                     &
                      l_ukca_conserve_h,                                       &
                      l_param_conv,                                            &
                      l_ctile,                                                 &
                      l_zon_av_ozone,                                          &
                      l_chem_environ_gas_scalars,                              &
                      l_chem_environ_co2_fld,                                  &
                      l_ukca_prescribech4,                                     &
                      l_use_classic_so4,                                       &
                      l_use_classic_soot,                                      &
                      l_use_classic_ocff,                                      &
                      l_use_classic_biogenic,                                  &
                      l_use_classic_seasalt,                                   &
                      l_use_gridbox_volume,                                    &
                      l_use_gridbox_mass,                                      &
                      l_environ_rel_humid,                                     &
                      l_environ_z_top,                                         &
                      l_fix_ukca_cloud_frac,                                   &
                      l_fix_improve_drydep,                                    &
                      l_fix_drydep_so2_water,                                  &
                      l_fix_ukca_h2dd_x,                                       &
                      l_fix_ukca_offox_h2o_fac,                                &
                      l_fix_ukca_h2so4_ystore,                                 &
                      l_mode_bhn_on,                                           &
                      l_mode_bln_on,                                           &
                      l_ddepaer,                                               &
                      l_aero_rainout,                                          &
                      l_cv_rainout,                                            &
                      l_impc_scav,                                             &
                      l_dust_slinn_impc_scav,                                  &
                      l_ukca_primss,                                           &
                      l_ukca_primsu,                                           &
                      l_ukca_primdu,                                           &
                      l_ukca_primbcoc,                                         &
                      l_ukca_prim_moc,                                         &
                      l_bcoc_bf,                                               &
                      l_bcoc_bm,                                               &
                      l_bcoc_ff,                                               &
                      l_ukca_scale_biom_aer_ems,                               &
                      l_ukca_fine_no3_prod,                                    &
                      l_ukca_coarse_no3_prod,                                  &
                      l_no3_prod_in_aero_step,                                 &
                      l_ukca_scale_sea_salt_ems,                               &
                      l_ukca_scale_marine_pom_ems,                             &
                      l_ukca_radaer,                                           &
                      l_ntpreq_n_activ_sum,                                    &
                      l_ntpreq_dryd_nuc_sol,                                   &
                      l_ukca_sfix,                                             &
                      l_fix_neg_pvol_wat,                                      &
                      l_fix_ukca_impscav,                                      &
                      l_fix_nacl_density,                                      &
                      l_improve_aero_drydep,                                   &
                      l_fix_ukca_water_content,                                &
                      l_fix_ukca_activate_pdf,                                 &
                      l_fix_ukca_activate_vert_rep,                            &
                      l_bug_repro_tke_index,                                   &
                      l_dust_ageing_on,                                        &
                      l_fix_ukca_hygroscopicities,                             &
                      l_skip_const_setup,                                      &
                      proc_bl_tracer_mix,                                      &
                      proc_calc_ozonecol,                                      &
                      proc_diag2d_copy_out,                                    &
                      proc_diag3d_copy_out,                                    &
                      error_message, error_routine)

! ----------------------------------------------------------------------
! Description:
!
!  Given the input configuration control data, check its validity and
!  set up UKCA's internal configuration data accordingly.
!  This includes some basic initialisation and everything required to
!  establish details of the selected configuration that will determine
!  the required tracers, non-transported prognostics (NTPs) and
!  environmental drivers.
!
! Method:
!
!  1. Set up UKCA's configurable constants to their default values or
!     values provided by keyword arguments with the prefix 'const_'.
!     This step should be skipped if constants have been set up by a previous
!     API call to 'ukca_constants_setup'. Keyword argument 'l_skip_const_setup'
!     sholuld be set to .TRUE. for this.
!  2. Copy the configuration variables provided as keyword arguments into
!     component variables with matching names in the UKCA (non-GLOMAP)
!     and/or GLOMAP configuration specification structures.
!     For keyword argument with the prefix 'proc_', set the appropriate
!     UKCA procedure pointer to access the given callback procedure.
!     For certain variables, default values are set here for use if input
!     values are not provided.
!     ------------------------------------------------------------------
!     Note: Values of variables that are inactive in the current
!     configuration will normally be ignored and not set. This is
!     intended to avoid the proliferation of spurious values that may be
!     misleading if the UKCA configuration data are interrogated by the
!     parent.
!     Exceptions are made for inactivated option variables that might
!     reasonably be expected to be relevant. Silently ignoring supplied
!     values for such option variables could be misleading so the values
!     are set and allowed through to subsequent validation checks
!     designed to trap un-supported option combinations.
!     Exceptions are also made for some 'general purpose' variables
!     with uses that are not specifically associated with particular
!     schemes or other options that determine whether they are in use.
!     ------------------------------------------------------------------
!  3. Check UKCA logicals are consistent and set internal UKCA values
!     based on the configuration data provided.
!  4. Initialise chemistry definition arrays.
!  5. Set up lists of tracer and NTP species required for the selected
!     chemical scheme.
!  6. Set up structure holding details of all NTPs required for the
!     selected UKCA configuration.
!  7. Set up lists of emitted species.
!  8. Set up lists of environmental driver fields required.
!  9. Clear environmental field data to ensure it is properly
!     initialised for subsequent set up via ukca_set_environment calls.
! 10. Initialise master diagnostics list and determine availability of
!     each diagnostic given the configuration.
! ----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: init_ukca_configuration,              &
    copy_config_vector, ukca_config, glomap_config,                            &
    i_ukca_chem_off, i_ukca_chem_trop, i_ukca_chem_raq,                        &
    i_ukca_chem_offline_be, i_ukca_chem_tropisop, i_ukca_chem_strattrop,       &
    i_ukca_chem_strat, i_ukca_chem_offline, i_ukca_chem_cristrat,              &
    i_age_reset_by_level, i_age_reset_by_height, i_strat_lbc_off,              &
    i_light_param_off, i_light_param_pr,                                       &
    i_top_none, i_top_bc, i_ukca_activation_off, i_ukca_activation_arg,        &
    i_dms_flux_off,                                                            &
    l_ukca_config_available,                                                   &
    template_proc_bl_tracer_mix, bl_tracer_mix,                                &
    template_proc_calc_ozonecol, calc_ozonecol,                                &
    template_proc_diag2d_copy_out, diag2d_copy_out,                            &
    template_proc_diag3d_copy_out, diag3d_copy_out

USE ukca_config_constants_mod, ONLY: l_ukca_constants_available
USE ukca_constants_setup_mod, ONLY: ukca_constants_setup

USE ukca_um_legacy_mod, ONLY: l_dust, l_twobin_dust, l_um_infrastructure

USE ukca_constants, ONLY: l_ukca_diurnal_isopems
USE ukca_init_mod, ONLY: ukca_init
USE ukca_chem1_dat, ONLY: ukca_chem1_init
USE ukca_ntp_mod, ONLY: ntp_init
USE asad_mod, ONLY: asad_mod_pre_setup_init, advt, speci, ctype
USE asad_inrats_mod, ONLY: asad_inrats_set_sp_lists
USE ukca_config_defs_mod, ONLY: ukca_set_config_defs, em_chem_spec, lbc_spec,  &
                                cfc_lumped
USE ukca_tracers_mod, ONLY: init_tracer_req
USE ukca_environment_req_mod, ONLY: init_environment_req
USE ukca_environment_mod, ONLY: clear_environment_fields
USE ukca_diagnostics_init_mod, ONLY: init_diagnostics
USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_ukca_uninit, errcode_value_missing,          &
                          i_error_method_abort, i_error_method_return,         &
                          i_error_method_warn_and_return

USE parkind1,               ONLY: jpim, jprb      ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Subroutine arguments.

! Each input configuration variable is an optional keyword argument with a
! matching component in UKCA or GLOMAP configuration structures.
! Refer to 'ukca_config_spec_type' & 'glomap_config_spec_type' definitions for
! descriptions (see 'ukca_config_specification_mod').
! Order of arguments within each type group should match 'ukca_config_spec_type'
! & 'glomap_config_spec_type' order.

! Additional optional arguments with the prefix 'const_' are used to override
! values of configurable constants defined in 'ukca_config_constants_mod'.

! Additional optional arguments with the prefix 'proc_' are used to supply
! parent call back routines complying with templates defined in
! 'ukca_config_specifcation_mod'.

! WARNING: keyword argument names are part of the UKCA API and any changes to
! existing names may affect backwards compatibility. To avoid this, new names
! can instead be added. Values set using the new names should then override
! corresponding values set using the old names.

INTEGER, TARGET, INTENT(OUT) :: error_code

INTEGER, OPTIONAL, INTENT(IN) :: row_length
INTEGER, OPTIONAL, INTENT(IN) :: rows
INTEGER, OPTIONAL, INTENT(IN) :: model_levels
INTEGER, OPTIONAL, INTENT(IN) :: bl_levels
INTEGER, OPTIONAL, INTENT(IN) :: nlev_ent_tr_mix
INTEGER, OPTIONAL, INTENT(IN) :: ntype
INTEGER, OPTIONAL, INTENT(IN) :: npft
INTEGER, OPTIONAL, INTENT(IN) :: i_brd_leaf
INTEGER, OPTIONAL, INTENT(IN) :: i_brd_leaf_dec
INTEGER, OPTIONAL, INTENT(IN) :: i_brd_leaf_eg_trop
INTEGER, OPTIONAL, INTENT(IN) :: i_brd_leaf_eg_temp
INTEGER, OPTIONAL, INTENT(IN) :: i_ndl_leaf
INTEGER, OPTIONAL, INTENT(IN) :: i_ndl_leaf_dec
INTEGER, OPTIONAL, INTENT(IN) :: i_ndl_leaf_eg
INTEGER, OPTIONAL, INTENT(IN) :: i_c3_grass
INTEGER, OPTIONAL, INTENT(IN) :: i_c3_crop
INTEGER, OPTIONAL, INTENT(IN) :: i_c3_pasture
INTEGER, OPTIONAL, INTENT(IN) :: i_c4_grass
INTEGER, OPTIONAL, INTENT(IN) :: i_c4_crop
INTEGER, OPTIONAL, INTENT(IN) :: i_c4_pasture
INTEGER, OPTIONAL, INTENT(IN) :: i_shrub
INTEGER, OPTIONAL, INTENT(IN) :: i_shrub_dec
INTEGER, OPTIONAL, INTENT(IN) :: i_shrub_eg
INTEGER, OPTIONAL, INTENT(IN) :: i_urban
INTEGER, OPTIONAL, INTENT(IN) :: i_lake
INTEGER, OPTIONAL, INTENT(IN) :: i_soil
INTEGER, OPTIONAL, INTENT(IN) :: i_ice
INTEGER, ALLOCATABLE, OPTIONAL, INTENT(IN) :: i_elev_ice(:)
INTEGER, OPTIONAL, INTENT(IN) :: i_ukca_chem
INTEGER, OPTIONAL, INTENT(IN) :: fixed_tropopause_level
INTEGER, OPTIONAL, INTENT(IN) :: i_ageair_reset_method
INTEGER, OPTIONAL, INTENT(IN) :: max_ageair_reset_level
INTEGER, OPTIONAL, INTENT(IN) :: i_error_method
INTEGER, OPTIONAL, INTENT(IN) :: i_ukca_chem_version
INTEGER, OPTIONAL, INTENT(IN) :: nrsteps
INTEGER, OPTIONAL, INTENT(IN) :: chem_timestep
INTEGER, OPTIONAL, INTENT(IN) :: dts0
INTEGER, OPTIONAL, INTENT(IN) :: nit
INTEGER, OPTIONAL, INTENT(IN) :: i_ukca_quasinewton_start
INTEGER, OPTIONAL, INTENT(IN) :: i_ukca_quasinewton_end
INTEGER, OPTIONAL, INTENT(IN) :: ukca_chem_seg_size
INTEGER, OPTIONAL, INTENT(IN) :: nlev_above_trop_o3_env
INTEGER, OPTIONAL, INTENT(IN) :: nlev_ch4_stratloss
INTEGER, OPTIONAL, INTENT(IN) :: i_ukca_topboundary
INTEGER, OPTIONAL, INTENT(IN) :: i_ukca_hetconfig
INTEGER, OPTIONAL, INTENT(IN) :: i_inferno_emi
INTEGER, OPTIONAL, INTENT(IN) :: i_ukca_dms_flux
INTEGER, OPTIONAL, INTENT(IN) :: i_ukca_light_param
INTEGER, OPTIONAL, INTENT(IN) :: i_strat_lbc_source
INTEGER, OPTIONAL, INTENT(IN) :: env_log_step
INTEGER, OPTIONAL, INTENT(IN) :: i_mode_nzts
INTEGER, OPTIONAL, INTENT(IN) :: ukca_mode_seg_size
INTEGER, OPTIONAL, INTENT(IN) :: i_mode_setup
INTEGER, OPTIONAL, INTENT(IN) :: i_mode_bln_param_method
INTEGER, OPTIONAL, INTENT(IN) :: i_mode_nucscav
INTEGER, OPTIONAL, INTENT(IN) :: i_ukca_activation_scheme
INTEGER, OPTIONAL, INTENT(IN) :: i_ukca_nwbins
INTEGER, OPTIONAL, INTENT(IN) :: i_ukca_tune_bc
INTEGER, OPTIONAL, INTENT(IN) :: i_primss_method
INTEGER, OPTIONAL, INTENT(IN) :: i_photol_scheme
INTEGER, OPTIONAL, INTENT(IN) :: i_photol_scheme_off
INTEGER, OPTIONAL, INTENT(IN) :: i_photol_scheme_strat_only
INTEGER, OPTIONAL, INTENT(IN) :: i_photol_scheme_2d
INTEGER, OPTIONAL, INTENT(IN) :: i_photol_scheme_fastjx

REAL, OPTIONAL, INTENT(IN) :: dzsoil_layer1
REAL, OPTIONAL, INTENT(IN) :: timestep
REAL, OPTIONAL, INTENT(IN) :: max_ageair_reset_height
REAL, OPTIONAL, INTENT(IN) :: max_z_for_offline_chem
REAL, OPTIONAL, INTENT(IN) :: fastjx_prescutoff
REAL, OPTIONAL, INTENT(IN) :: mode_parfrac
REAL, OPTIONAL, INTENT(IN) :: seadms_ems_scaling
REAL, OPTIONAL, INTENT(IN) :: sea_salt_ems_scaling
REAL, OPTIONAL, INTENT(IN) :: marine_pom_ems_scaling
REAL, OPTIONAL, INTENT(IN) :: lightnox_scale_fac
REAL, OPTIONAL, INTENT(IN) :: anth_so2_ems_scaling
REAL, OPTIONAL, INTENT(IN) :: soa_yield_scaling_mt
REAL, OPTIONAL, INTENT(IN) :: soa_yield_scaling_isop
REAL, OPTIONAL, INTENT(IN) :: dry_depvel_so2_scaling
REAL, OPTIONAL, INTENT(IN) :: mode_activation_dryr
REAL, OPTIONAL, INTENT(IN) :: dry_depvel_acc_scaling
REAL, OPTIONAL, INTENT(IN) :: acc_cor_scav_scaling
REAL, OPTIONAL, INTENT(IN) :: mode_incld_so2_rfrac
REAL, OPTIONAL, INTENT(IN) :: biom_aer_ems_scaling
REAL, OPTIONAL, INTENT(IN) :: ph_fit_coeff_a
REAL, OPTIONAL, INTENT(IN) :: ph_fit_coeff_b
REAL, OPTIONAL, INTENT(IN) :: ph_fit_intercept
REAL, OPTIONAL, INTENT(IN) :: sigwmin
REAL, OPTIONAL, INTENT(IN) :: hno3_uptake_coeff
REAL, OPTIONAL, INTENT(IN) :: sigma_updraught_scaling
REAL, OPTIONAL, INTENT(IN) :: const_rmol
REAL, OPTIONAL, INTENT(IN) :: const_tfs
REAL, OPTIONAL, INTENT(IN) :: const_rho_water
REAL, OPTIONAL, INTENT(IN) :: const_rhosea
REAL, OPTIONAL, INTENT(IN) :: const_lc
REAL, OPTIONAL, INTENT(IN) :: const_avogadro
REAL, OPTIONAL, INTENT(IN) :: const_boltzmann
REAL, OPTIONAL, INTENT(IN) :: const_rho_so4

LOGICAL, OPTIONAL, INTENT(IN) :: l_cal360
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_chem_aero
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_mode
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_tropopause_level
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_ageair
LOGICAL, OPTIONAL, INTENT(IN) :: l_blankout_invalid_diags
LOGICAL, OPTIONAL, INTENT(IN) :: l_enable_diag_um
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_persist_off
LOGICAL, OPTIONAL, INTENT(IN) :: l_timer
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_emissions_off
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_drydep_off
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_wetdep_off
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_scale_ppe
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_asad_columns
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_asad_full
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_debug_asad
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_intdd
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_ddepo3_ocean
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_ddep_lev1
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_dry_dep_so2wet
LOGICAL, OPTIONAL, INTENT(IN) :: l_deposition_jules
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_quasinewton
LOGICAL, OPTIONAL, INTENT(IN) :: l_tracer_lumping
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_ro2_ntp
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_ro2_perm
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_intph
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_het_psc
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_limit_nat
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_sa_clim
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_trophet
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_classic_hetchem
LOGICAL, OPTIONAL, INTENT(IN) :: l_use_photolysis
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_ibvoc
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_inferno
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_inferno_ch4
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_so2ems_expvolc
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_so2ems_plumeria
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_qch4inter
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_emsdrvn_ch4
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_enable_seadms_ems
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_scale_seadms_ems
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_linox_scaling
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_scale_soa_yield_mt
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_scale_soa_yield_isop
LOGICAL, OPTIONAL, INTENT(IN) :: l_support_ems_vertprof
LOGICAL, OPTIONAL, INTENT(IN) :: l_support_ems_gridbox_units
LOGICAL, OPTIONAL, INTENT(IN) :: l_suppress_ems
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_h2o_feedback
LOGICAL, OPTIONAL, INTENT(IN) :: l_param_conv
LOGICAL, OPTIONAL, INTENT(IN) :: l_ctile
LOGICAL, OPTIONAL, INTENT(IN) :: l_zon_av_ozone
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_conserve_h
LOGICAL, OPTIONAL, INTENT(IN) :: l_chem_environ_gas_scalars
LOGICAL, OPTIONAL, INTENT(IN) :: l_chem_environ_co2_fld
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_prescribech4
LOGICAL, OPTIONAL, INTENT(IN) :: l_use_classic_so4
LOGICAL, OPTIONAL, INTENT(IN) :: l_use_classic_soot
LOGICAL, OPTIONAL, INTENT(IN) :: l_use_classic_ocff
LOGICAL, OPTIONAL, INTENT(IN) :: l_use_classic_biogenic
LOGICAL, OPTIONAL, INTENT(IN) :: l_use_classic_seasalt
LOGICAL, OPTIONAL, INTENT(IN) :: l_use_gridbox_volume
LOGICAL, OPTIONAL, INTENT(IN) :: l_use_gridbox_mass
LOGICAL, OPTIONAL, INTENT(IN) :: l_environ_rel_humid
LOGICAL, OPTIONAL, INTENT(IN) :: l_environ_z_top
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_ukca_cloud_frac
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_improve_drydep
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_drydep_so2_water
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_ukca_h2dd_x
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_ukca_offox_h2o_fac
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_ukca_h2so4_ystore
LOGICAL, OPTIONAL, INTENT(IN) :: l_mode_bhn_on
LOGICAL, OPTIONAL, INTENT(IN) :: l_mode_bln_on
LOGICAL, OPTIONAL, INTENT(IN) :: l_ddepaer
LOGICAL, OPTIONAL, INTENT(IN) :: l_aero_rainout
LOGICAL, OPTIONAL, INTENT(IN) :: l_cv_rainout
LOGICAL, OPTIONAL, INTENT(IN) :: l_impc_scav
LOGICAL, OPTIONAL, INTENT(IN) :: l_dust_slinn_impc_scav
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_primss
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_primsu
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_primdu
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_primbcoc
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_prim_moc
LOGICAL, OPTIONAL, INTENT(IN) :: l_bcoc_bf
LOGICAL, OPTIONAL, INTENT(IN) :: l_bcoc_bm
LOGICAL, OPTIONAL, INTENT(IN) :: l_bcoc_ff
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_scale_biom_aer_ems
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_fine_no3_prod
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_coarse_no3_prod
LOGICAL, OPTIONAL, INTENT(IN) :: l_no3_prod_in_aero_step
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_scale_sea_salt_ems
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_scale_marine_pom_ems
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_radaer
LOGICAL, OPTIONAL, INTENT(IN) :: l_ntpreq_n_activ_sum
LOGICAL, OPTIONAL, INTENT(IN) :: l_ntpreq_dryd_nuc_sol
LOGICAL, OPTIONAL, INTENT(IN) :: l_ukca_sfix
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_neg_pvol_wat
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_ukca_impscav
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_nacl_density
LOGICAL, OPTIONAL, INTENT(IN) :: l_improve_aero_drydep
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_ukca_water_content
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_ukca_activate_pdf
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_ukca_activate_vert_rep
LOGICAL, OPTIONAL, INTENT(IN) :: l_bug_repro_tke_index
LOGICAL, OPTIONAL, INTENT(IN) :: l_fix_ukca_hygroscopicities
LOGICAL, OPTIONAL, INTENT(IN) :: l_dust_ageing_on

! Control argument for skipping set up of constants if already done
LOGICAL, OPTIONAL, INTENT(IN) :: l_skip_const_setup

PROCEDURE(template_proc_bl_tracer_mix), OPTIONAL :: proc_bl_tracer_mix
PROCEDURE(template_proc_calc_ozonecol), OPTIONAL :: proc_calc_ozonecol
PROCEDURE(template_proc_diag2d_copy_out), OPTIONAL :: proc_diag2d_copy_out
PROCEDURE(template_proc_diag3d_copy_out), OPTIONAL :: proc_diag3d_copy_out

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

INTEGER, POINTER :: error_code_ptr ! Pointer for selecting error code variable
INTEGER :: i                       ! loop counter

LOGICAL :: l_setup_constants       ! True if constants need setting up here
LOGICAL :: l_be_scheme_selected    ! True if B-E solver required for chemistry
LOGICAL :: l_nr_scheme_selected    ! True if N-R solver required for chemistry
LOGICAL :: l_strat_scheme_selected ! True if a Stratospheric scheme is
                                   ! selected for chemistry

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SETUP'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Use parent supplied argument for return code. Note that this argument is
! redundant if UKCA is configured to abort on error and may be made optional
! in future, being replaced by an internal error code variable when absent.
error_code_ptr => error_code

! Set defaults for output arguments
error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Set all configurable constants to defaults or values provided
! if not asked to skip this step

IF (PRESENT(l_skip_const_setup)) THEN
  l_setup_constants = .NOT. l_skip_const_setup
ELSE
  l_setup_constants = .TRUE.
END IF

IF (l_setup_constants) THEN

  CALL ukca_constants_setup(error_code_ptr,                                    &
                            const_rmol=const_rmol,                             &
                            const_tfs=const_tfs,                               &
                            const_rho_water=const_rho_water,                   &
                            const_rhosea=const_rhosea,                         &
                            const_lc=const_lc,                                 &
                            const_avogadro=const_avogadro,                     &
                            const_boltzmann=const_boltzmann,                   &
                            const_rho_so4=const_rho_so4,                       &
                            error_message=error_message,                       &
                            error_routine=error_routine)

ELSE

  ! Asked to skip, so check that constants have been set up previously
  IF (.NOT. l_ukca_constants_available) THEN
    error_code_ptr = errcode_ukca_uninit
    IF (PRESENT(error_message))                                                &
      error_message = 'Configurable UKCA constants have not been set up'
    IF (PRESENT(error_routine)) error_routine = RoutineName
  END IF

END IF

IF (error_code_ptr > 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Set all configuration data to default values
CALL init_ukca_configuration()

! Collate input data specifying the UKCA configuration
! (loosely following 'ukca_config_spec_type' & 'glomap_config_spec_type' order
! where practical given the required logic).
! Default values are set here if applicable.
! Input values for configuration variables that are inactive in the current
! configuration are ignored if appropriate (see subroutine method Step 1 for
! details).

! -- Context information required for all configurations ---------------

ukca_config%row_length = 1
ukca_config%rows = 1
ukca_config%model_levels = 1

IF (PRESENT(row_length)) ukca_config%row_length = row_length
IF (PRESENT(rows)) ukca_config%rows = rows
IF (PRESENT(model_levels)) ukca_config%model_levels = model_levels
IF (PRESENT(timestep)) ukca_config%timestep = timestep

! -- General UKCA configuration options --------------------------------

ukca_config%i_ukca_chem = i_ukca_chem_off
ukca_config%i_error_method = i_error_method_abort

IF (PRESENT(i_ukca_chem)) ukca_config%i_ukca_chem = i_ukca_chem
IF (PRESENT(l_ukca_ageair)) ukca_config%l_ukca_ageair = l_ukca_ageair
IF (PRESENT(l_blankout_invalid_diags))                                         &
 ukca_config%l_blankout_invalid_diags = l_blankout_invalid_diags
IF (PRESENT(l_enable_diag_um)) ukca_config%l_enable_diag_um = l_enable_diag_um
IF (PRESENT(l_ukca_persist_off))                                               &
  ukca_config%l_ukca_persist_off = l_ukca_persist_off
IF (PRESENT(l_ukca_mode)) ukca_config%l_ukca_mode = l_ukca_mode
IF (PRESENT(l_timer)) ukca_config%l_timer = l_timer
IF (PRESENT(i_error_method))                                                   &
  ukca_config%i_error_method = i_error_method

IF (ukca_config%i_ukca_chem /= i_ukca_chem_off) THEN

  ukca_config%bl_levels = 1

  IF (PRESENT(bl_levels)) ukca_config%bl_levels = bl_levels

  IF (ukca_config%i_ukca_chem == i_ukca_chem_offline .OR.                      &
      ukca_config%i_ukca_chem == i_ukca_chem_offline_be) THEN
    ukca_config%l_ukca_chem_aero = .TRUE.
  ELSE
    IF (PRESENT(l_ukca_chem_aero))                                             &
      ukca_config%l_ukca_chem_aero = l_ukca_chem_aero
  END IF

  IF (PRESENT(l_fix_tropopause_level))                                         &
    ukca_config%l_fix_tropopause_level = l_fix_tropopause_level
  IF (ukca_config%l_fix_tropopause_level) THEN
    ukca_config%fixed_tropopause_level = 1
    IF (PRESENT(fixed_tropopause_level))                                       &
      ukca_config%fixed_tropopause_level = fixed_tropopause_level
  END IF

  IF (PRESENT(l_ukca_emissions_off))                                           &
    ukca_config%l_ukca_emissions_off = l_ukca_emissions_off
  IF (PRESENT(l_ukca_drydep_off))                                              &
    ukca_config%l_ukca_drydep_off = l_ukca_drydep_off
  IF (PRESENT(l_ukca_wetdep_off))                                              &
    ukca_config%l_ukca_wetdep_off = l_ukca_wetdep_off

END IF

! PPE scalings
IF (PRESENT(l_ukca_scale_ppe)) ukca_config%l_ukca_scale_ppe = l_ukca_scale_ppe

! Age-of-air configuration

IF (ukca_config%l_ukca_ageair) THEN

  IF (PRESENT(i_ageair_reset_method))                                          &
    ukca_config%i_ageair_reset_method = i_ageair_reset_method

  SELECT CASE(ukca_config%i_ageair_reset_method)
  CASE (i_age_reset_by_level)
    IF (PRESENT(max_ageair_reset_level))                                       &
      ukca_config%max_ageair_reset_level = max_ageair_reset_level
  CASE (i_age_reset_by_height)
    IF (PRESENT(max_ageair_reset_height))                                      &
      ukca_config%max_ageair_reset_height = max_ageair_reset_height
  END SELECT

END IF

! Determine which solver is to be used for chemistry; set switches for later use
IF (ukca_config%i_ukca_chem == i_ukca_chem_off) THEN
  l_be_scheme_selected = .FALSE.
  l_nr_scheme_selected = .FALSE.
ELSE IF (ukca_config%i_ukca_chem == i_ukca_chem_trop .OR.                      &
    ukca_config%i_ukca_chem == i_ukca_chem_raq .OR.                            &
    ukca_config%i_ukca_chem == i_ukca_chem_offline_be) THEN
  l_be_scheme_selected = .TRUE.
  l_nr_scheme_selected = .FALSE.
ELSE
  l_be_scheme_selected = .FALSE.
  l_nr_scheme_selected = .TRUE.
END IF

! -- Chemistry configuration options -----------------------------------

l_strat_scheme_selected = (ukca_config%i_ukca_chem == i_ukca_chem_strat        &
                      .OR. ukca_config%i_ukca_chem == i_ukca_chem_strattrop    &
                      .OR. ukca_config%i_ukca_chem == i_ukca_chem_cristrat )

IF (ukca_config%i_ukca_chem /= i_ukca_chem_off) THEN

  ukca_config%ntype = 0
  ukca_config%npft = 0

  IF (PRESENT(chem_timestep)) ukca_config%chem_timestep = chem_timestep
  IF (PRESENT(l_ukca_asad_columns))                                            &
    ukca_config%l_ukca_asad_columns = l_ukca_asad_columns
  IF (PRESENT(l_ukca_asad_full))                                               &
    ukca_config%l_ukca_asad_full = l_ukca_asad_full
  IF (PRESENT(l_ukca_debug_asad))                                              &
    ukca_config%l_ukca_debug_asad = l_ukca_debug_asad

  IF (.NOT. (ukca_config%l_ukca_drydep_off)) THEN
    IF (PRESENT(l_ukca_intdd)) ukca_config%l_ukca_intdd = l_ukca_intdd
  END IF

  IF (PRESENT(l_ukca_intph)) ukca_config%l_ukca_intph = l_ukca_intph

  ! Solver-specific configuration

  IF (l_be_scheme_selected) THEN

    ukca_config%dts0 = 300
    ukca_config%nit = 8

    IF (PRESENT(dts0)) ukca_config%dts0 = dts0
    IF (PRESENT(nit)) ukca_config%nit = nit

  ELSE IF (l_nr_scheme_selected) THEN

    IF (PRESENT(i_ukca_chem_version))                                          &
      ukca_config%i_ukca_chem_version = i_ukca_chem_version
    IF (PRESENT(nrsteps)) ukca_config%nrsteps = nrsteps
    IF (PRESENT(l_ukca_quasinewton))                                           &
      ukca_config%l_ukca_quasinewton = l_ukca_quasinewton

    IF (ukca_config%l_ukca_quasinewton) THEN
      IF (PRESENT(i_ukca_quasinewton_start))                                   &
        ukca_config%i_ukca_quasinewton_start = i_ukca_quasinewton_start
      IF (PRESENT(i_ukca_quasinewton_end))                                     &
        ukca_config%i_ukca_quasinewton_end = i_ukca_quasinewton_end
    END IF

  END IF

  ! Column-based run configuration

  IF (ukca_config%l_ukca_asad_columns) THEN

    ukca_config%ukca_chem_seg_size = ukca_config%model_levels

    IF (PRESENT(ukca_chem_seg_size))                                           &
      ukca_config%ukca_chem_seg_size = ukca_chem_seg_size

  END IF

  ! Configuration specific to explicit B-E Offline Oxidants scheme
  IF (ukca_config%i_ukca_chem == i_ukca_chem_offline_be) THEN
    IF (PRESENT(max_z_for_offline_chem))                                       &
      ukca_config%max_z_for_offline_chem = max_z_for_offline_chem
  END IF

  ! Configuration specific to troposphere only schemes

  IF (ukca_config%i_ukca_chem == i_ukca_chem_trop .OR.                         &
      ukca_config%i_ukca_chem == i_ukca_chem_raq) THEN

    ! By default, set
    ! - no. of levels above tropopause for overwrite of O3 & HNO3 = model_levels
    ! - no. of top levels at which stratopsheric loss of CH4 is applied = 0
    ! such that no levels are overwritten.
    ! Exception: no defaults are set if using the UM infrastructure. Default
    ! number of levels above tropopause for overwrite is then set in subroutine
    ! ukca_stratf (based on the total number of model levels) and default number
    ! of top levels for CH4 loss is set in ukca_chemistry_ctl.
    ! In both cases, provision of that functionality within UKCA is deprecated
    ! and the UM should be modified to set both values explicitly.
    IF (.NOT. l_um_infrastructure) THEN
      ukca_config%nlev_above_trop_o3_env = ukca_config%model_levels
      ukca_config%nlev_ch4_stratloss = 0
    END IF

    IF (PRESENT(nlev_above_trop_o3_env))                                       &
      ukca_config%nlev_above_trop_o3_env = nlev_above_trop_o3_env
    IF (PRESENT(nlev_ch4_stratloss))                                           &
      ukca_config%nlev_ch4_stratloss = nlev_ch4_stratloss

  END IF

  ! Configuration specific to stratospheric schemes

  IF ( l_strat_scheme_selected ) THEN

    ukca_config%l_tracer_lumping = .TRUE.
    ukca_config%i_ukca_topboundary = i_top_none

    IF (PRESENT(l_tracer_lumping))                                             &
      ukca_config%l_tracer_lumping = l_tracer_lumping
    IF (PRESENT(i_ukca_topboundary))                                           &
      ukca_config%i_ukca_topboundary = i_ukca_topboundary
    IF (PRESENT(l_ukca_ro2_ntp)) ukca_config%l_ukca_ro2_ntp = l_ukca_ro2_ntp
    IF (PRESENT(l_ukca_ro2_perm)) ukca_config%l_ukca_ro2_perm = l_ukca_ro2_perm

  END IF

  ! Configuration specific to interactive dry-deposition scheme
  IF (ukca_config%l_ukca_intdd) THEN
    ! Context information
    IF (PRESENT(ntype)) ukca_config%ntype = ntype
    IF (PRESENT(npft)) ukca_config%npft = npft
    IF (PRESENT(i_brd_leaf)) ukca_config%i_brd_leaf = i_brd_leaf
    IF (PRESENT(i_brd_leaf_dec)) ukca_config%i_brd_leaf_dec = i_brd_leaf_dec
    IF (PRESENT(i_brd_leaf_eg_trop))                                           &
      ukca_config%i_brd_leaf_eg_trop = i_brd_leaf_eg_trop
    IF (PRESENT(i_brd_leaf_eg_temp))                                           &
      ukca_config%i_brd_leaf_eg_temp = i_brd_leaf_eg_temp
    IF (PRESENT(i_ndl_leaf)) ukca_config%i_ndl_leaf = i_ndl_leaf
    IF (PRESENT(i_ndl_leaf_dec)) ukca_config%i_ndl_leaf_dec = i_ndl_leaf_dec
    IF (PRESENT(i_ndl_leaf_eg)) ukca_config%i_ndl_leaf_eg = i_ndl_leaf_eg
    IF (PRESENT(i_c3_grass)) ukca_config%i_c3_grass = i_c3_grass
    IF (PRESENT(i_c3_crop)) ukca_config%i_c3_crop = i_c3_crop
    IF (PRESENT(i_c3_pasture)) ukca_config%i_c3_pasture = i_c3_pasture
    IF (PRESENT(i_c4_grass)) ukca_config%i_c4_grass = i_c4_grass
    IF (PRESENT(i_c4_crop)) ukca_config%i_c4_crop = i_c4_crop
    IF (PRESENT(i_c4_pasture)) ukca_config%i_c4_pasture = i_c4_pasture
    IF (PRESENT(i_shrub)) ukca_config%i_shrub = i_shrub
    IF (PRESENT(i_shrub_dec)) ukca_config%i_shrub_dec = i_shrub_dec
    IF (PRESENT(i_shrub_eg)) ukca_config%i_shrub_eg = i_shrub_eg
    IF (PRESENT(i_urban)) ukca_config%i_urban = i_urban
    IF (PRESENT(i_lake)) ukca_config%i_lake = i_lake
    IF (PRESENT(i_soil)) ukca_config%i_soil = i_soil
    IF (PRESENT(i_ice)) ukca_config%i_ice = i_ice
    IF (PRESENT(i_elev_ice))                                                   &
      CALL copy_config_vector(i_elev_ice,ukca_config%i_elev_ice)
    IF (PRESENT(dzsoil_layer1)) ukca_config%dzsoil_layer1 = dzsoil_layer1
    ! Chemistry options
    IF (PRESENT(l_ukca_ddepo3_ocean))                                          &
      ukca_config%l_ukca_ddepo3_ocean = l_ukca_ddepo3_ocean
    IF (PRESENT(l_ukca_ddep_lev1))                                             &
      ukca_config%l_ukca_ddep_lev1 = l_ukca_ddep_lev1
    IF (PRESENT(l_ukca_dry_dep_so2wet))                                        &
      ukca_config%l_ukca_dry_dep_so2wet = l_ukca_dry_dep_so2wet
    IF (PRESENT(l_deposition_jules))                                           &
      ukca_config%l_deposition_jules = l_deposition_jules
  END IF

  ! Configuration specific to interactive cloud pH scheme
  IF (ukca_config%l_ukca_intph) THEN
    IF (PRESENT(ph_fit_coeff_a)) ukca_config%ph_fit_coeff_a = ph_fit_coeff_a
    IF (PRESENT(ph_fit_coeff_b)) ukca_config%ph_fit_coeff_b = ph_fit_coeff_b
    IF (PRESENT(ph_fit_intercept))                                             &
      ukca_config%ph_fit_intercept = ph_fit_intercept
  END IF

  ! Configuration specific to aerosol chemistry
  IF (ukca_config%l_ukca_chem_aero) THEN
    IF (PRESENT(l_ukca_scale_soa_yield_mt))                                    &
      ukca_config%l_ukca_scale_soa_yield_mt = l_ukca_scale_soa_yield_mt
    IF (ukca_config%l_ukca_scale_soa_yield_mt) THEN
      IF (PRESENT(soa_yield_scaling_mt))                                       &
        ukca_config%soa_yield_scaling_mt = soa_yield_scaling_mt
    END IF
    IF (PRESENT(l_ukca_scale_soa_yield_isop))                                  &
      ukca_config%l_ukca_scale_soa_yield_isop = l_ukca_scale_soa_yield_isop
    IF (ukca_config%l_ukca_scale_soa_yield_isop) THEN
      IF (PRESENT(soa_yield_scaling_isop))                                     &
        ukca_config%soa_yield_scaling_isop = soa_yield_scaling_isop
    END IF

  END IF

  IF (ukca_config%l_ukca_scale_ppe) THEN
    IF (PRESENT(dry_depvel_so2_scaling))                                       &
      ukca_config%dry_depvel_so2_scaling = dry_depvel_so2_scaling
  END IF

  ! -- Chemistry - Heterogeneous chemistry --

  IF (PRESENT(l_ukca_het_psc)) ukca_config%l_ukca_het_psc = l_ukca_het_psc
  IF (PRESENT(l_ukca_trophet)) ukca_config%l_ukca_trophet = l_ukca_trophet
  IF (PRESENT(l_ukca_classic_hetchem))                                         &
    ukca_config%l_ukca_classic_hetchem = l_ukca_classic_hetchem

  ! Stratospheric heterogeneous chemistry options

  IF (ukca_config%l_ukca_het_psc) THEN

    ukca_config%i_ukca_hetconfig = 0  ! i.e. 5 reactions

    IF (PRESENT(i_ukca_hetconfig))                                             &
      ukca_config%i_ukca_hetconfig = i_ukca_hetconfig
    IF (PRESENT(l_ukca_limit_nat))                                             &
      ukca_config%l_ukca_limit_nat = l_ukca_limit_nat
    IF (PRESENT(l_ukca_sa_clim)) ukca_config%l_ukca_sa_clim = l_ukca_sa_clim

  END IF

  ! -- Photolysis options --
  IF (PRESENT(l_use_photolysis)) ukca_config%l_use_photolysis = l_use_photolysis

END IF  ! i_ukca_chem /= i_ukca_chem_off

! -- Context information with dependencies on chem. setup and/or emissions --

IF ((ukca_config%i_ukca_chem /= i_ukca_chem_off .AND.                          &
     .NOT. ukca_config%l_ukca_emissions_off) .OR.                              &
    ukca_config%i_ukca_topboundary >= i_top_bc .OR.                            &
    ukca_config%l_use_photolysis) THEN
  IF (PRESENT(l_cal360)) ukca_config%l_cal360 = l_cal360
END IF

! -- UKCA emissions configuration options ------------------------------

IF (ukca_config%i_ukca_chem /= i_ukca_chem_off .AND.                           &
    .NOT. ukca_config%l_ukca_emissions_off) THEN

  ukca_config%nlev_ent_tr_mix = 0
  ukca_config%i_ukca_dms_flux = i_dms_flux_off

  IF (PRESENT(l_ukca_ibvoc)) ukca_config%l_ukca_ibvoc = l_ukca_ibvoc
  IF (PRESENT(l_ukca_inferno)) ukca_config%l_ukca_inferno = l_ukca_inferno
  IF (PRESENT(l_ukca_so2ems_expvolc))                                          &
    ukca_config%l_ukca_so2ems_expvolc = l_ukca_so2ems_expvolc
  IF (PRESENT(l_ukca_qch4inter)) ukca_config%l_ukca_qch4inter = l_ukca_qch4inter
  IF (PRESENT(l_ukca_emsdrvn_ch4))                                             &
    ukca_config%l_ukca_emsdrvn_ch4 = l_ukca_emsdrvn_ch4
  IF (PRESENT(l_support_ems_vertprof))                                         &
    ukca_config%l_support_ems_vertprof = l_support_ems_vertprof
  IF (PRESENT(l_support_ems_gridbox_units))                                    &
    ukca_config%l_support_ems_gridbox_units = l_support_ems_gridbox_units
  IF (PRESENT(l_suppress_ems))                                                 &
    ukca_config%l_suppress_ems = l_suppress_ems

  ! Explosive volcanic emissions options
  IF (ukca_config%l_ukca_so2ems_expvolc) THEN
    IF (PRESENT(l_ukca_so2ems_plumeria))                                       &
      ukca_config%l_ukca_so2ems_plumeria = l_ukca_so2ems_plumeria
  END IF

  ! Lightning NOx emissions options
  ! Always off for Offline Oxidants chemistry; default to Price & Rind for full
  ! chemistry
  IF (ukca_config%i_ukca_chem == i_ukca_chem_offline .OR.                      &
      ukca_config%i_ukca_chem == i_ukca_chem_offline_be) THEN
    ! Set NOx emissions off for Offline Oxidants chemistry
    ukca_config%i_ukca_light_param = i_light_param_off
  ELSE
    ! Default to Prince & Rind scheme
    ukca_config%i_ukca_light_param = i_light_param_pr
    IF (PRESENT(i_ukca_light_param))                                           &
      ukca_config%i_ukca_light_param = i_ukca_light_param
  END IF
  IF (ukca_config%i_ukca_light_param /= i_light_param_off) THEN
    IF (PRESENT(l_ukca_linox_scaling))                                         &
      ukca_config%l_ukca_linox_scaling = l_ukca_linox_scaling
  END IF

  ! Inferno scheme options
  IF ( ukca_config%l_ukca_inferno) THEN
    IF (PRESENT(l_ukca_inferno_ch4))                                           &
      ukca_config%l_ukca_inferno_ch4 = l_ukca_inferno_ch4
    IF (PRESENT(i_inferno_emi)) ukca_config%i_inferno_emi = i_inferno_emi
  END IF

  ! Marine DMS emissions options.
  ! Note that Marine DMS emissions should normally be controlled by setting
  ! 'i_ukca_dms_flux'.
  ! Use of 'l_ukca_enable_seadms_ems' is deprecated but overrides any potential
  ! setting of 'i_ukca_dms_flux' if it is present and set to .FALSE.
  IF (PRESENT(l_ukca_enable_seadms_ems)) THEN
    IF (l_ukca_enable_seadms_ems) THEN
      ukca_config%l_ukca_enable_seadms_ems = .TRUE.
      IF (PRESENT(i_ukca_dms_flux))                                            &
        ukca_config%i_ukca_dms_flux = i_ukca_dms_flux
    END IF
  ELSE
    IF (PRESENT(i_ukca_dms_flux)) ukca_config%i_ukca_dms_flux = i_ukca_dms_flux
  END IF
  IF (ukca_config%i_ukca_dms_flux /= i_dms_flux_off) THEN
    IF (PRESENT(l_ukca_scale_seadms_ems))                                      &
      ukca_config%l_ukca_scale_seadms_ems = l_ukca_scale_seadms_ems
  END IF

  ! Partitioning & scaling factors

  IF (ukca_config%l_ukca_chem_aero) THEN
    IF (PRESENT(mode_parfrac)) ukca_config%mode_parfrac = mode_parfrac
  END IF

  IF (ukca_config%l_ukca_scale_seadms_ems) THEN
    IF (PRESENT(seadms_ems_scaling))                                           &
      ukca_config%seadms_ems_scaling = seadms_ems_scaling
  END IF

  IF (ukca_config%i_ukca_light_param /= i_light_param_off) THEN
    IF (PRESENT(lightnox_scale_fac))                                           &
      ukca_config%lightnox_scale_fac = lightnox_scale_fac
  END IF

  ! Number of grid levels for tr_mix fields is required if emissions updates
  ! are not suppressed
  IF (.NOT. ukca_config%l_suppress_ems .AND. PRESENT(nlev_ent_tr_mix))         &
    ukca_config%nlev_ent_tr_mix = nlev_ent_tr_mix

  IF (ukca_config%l_ukca_scale_ppe) THEN
    IF (PRESENT(anth_so2_ems_scaling))                                         &
      ukca_config%anth_so2_ems_scaling = anth_so2_ems_scaling
  END IF

END IF

! -- UKCA feedback configuration options -------------------------------

IF (ukca_config%i_ukca_chem /= i_ukca_chem_off) THEN

  IF (PRESENT(l_ukca_h2o_feedback))                                            &
    ukca_config%l_ukca_h2o_feedback = l_ukca_h2o_feedback

  IF (ukca_config%l_ukca_h2o_feedback) THEN
    IF (PRESENT(l_ukca_conserve_h))                                            &
      ukca_config%l_ukca_conserve_h = l_ukca_conserve_h
  END IF

END IF

! -- UKCA environmental driver configuration options -------------------
!    (not including options with GLOMAP dependencies)

ukca_config%env_log_step = 1

IF (ukca_config%i_ukca_chem /= i_ukca_chem_off) THEN

  IF (.NOT. ukca_config%l_ukca_wetdep_off .OR. l_strat_scheme_selected ) THEN
    IF (PRESENT(l_param_conv)) ukca_config%l_param_conv = l_param_conv
  END IF

  IF (PRESENT(l_chem_environ_gas_scalars))                                     &
    ukca_config%l_chem_environ_gas_scalars = l_chem_environ_gas_scalars
  IF (PRESENT(l_chem_environ_co2_fld))                                         &
    ukca_config%l_chem_environ_co2_fld = l_chem_environ_co2_fld

  IF (.NOT. ukca_config%l_ukca_emissions_off) THEN
    IF (PRESENT(l_ukca_prescribech4))                                          &
      ukca_config%l_ukca_prescribech4 = l_ukca_prescribech4
  END IF

  IF (ukca_config%l_ukca_classic_hetchem .OR. ukca_config%l_ukca_sa_clim) THEN
    IF (PRESENT(l_use_classic_so4))                                            &
      ukca_config%l_use_classic_so4 = l_use_classic_so4
  END IF

  IF (ukca_config%l_ukca_classic_hetchem) THEN
    IF (PRESENT(l_use_classic_soot))                                           &
      ukca_config%l_use_classic_soot = l_use_classic_soot
    IF (PRESENT(l_use_classic_ocff))                                           &
      ukca_config%l_use_classic_ocff = l_use_classic_ocff
    IF (PRESENT(l_use_classic_biogenic))                                       &
      ukca_config%l_use_classic_biogenic = l_use_classic_biogenic
    IF (PRESENT(l_use_classic_seasalt))                                        &
      ukca_config%l_use_classic_seasalt = l_use_classic_seasalt
  END IF

  IF (ukca_config%i_ukca_chem == i_ukca_chem_strattrop .OR.                    &
      ukca_config%i_ukca_chem == i_ukca_chem_cristrat .OR.                     &
      ukca_config%i_ukca_chem == i_ukca_chem_offline .OR.                      &
      ukca_config%i_ukca_chem == i_ukca_chem_offline_be .OR.                   &
      ukca_config%l_enable_diag_um) THEN
    IF (PRESENT(l_use_gridbox_volume))                                         &
      ukca_config%l_use_gridbox_volume = l_use_gridbox_volume
  END IF

  ! Gridbox airmass is required for some diagnostics, the GLOMAP aerosol scheme
  ! and for applying lower boundary conditions during processing of emissions
  ! for Stratospheric chemistry schemes
  IF ( ukca_config%l_enable_diag_um .OR. ukca_config%l_ukca_mode .OR.          &
   (l_strat_scheme_selected .AND. .NOT. ukca_config%l_ukca_emissions_off)) THEN
    IF (PRESENT(l_use_gridbox_mass))                                           &
       ukca_config%l_use_gridbox_mass = l_use_gridbox_mass
  END IF

  IF (PRESENT(l_environ_rel_humid))                                            &
    ukca_config%l_environ_rel_humid = l_environ_rel_humid

END IF

IF (PRESENT(l_environ_z_top))                                                  &
  ukca_config%l_environ_z_top = l_environ_z_top
IF (PRESENT(env_log_step))                                                     &
  ukca_config%env_log_step = env_log_step

! Configuration specific to tropospheric schemes
IF (ukca_config%i_ukca_chem == i_ukca_chem_trop .OR.                           &
    ukca_config%i_ukca_chem == i_ukca_chem_raq .OR.                            &
    ukca_config%i_ukca_chem == i_ukca_chem_tropisop) THEN
  IF (PRESENT(l_zon_av_ozone)) ukca_config%l_zon_av_ozone = l_zon_av_ozone
END IF

! Configuration specific to stratospheric schemes

IF (l_strat_scheme_selected .AND. .NOT. ukca_config%l_ukca_emissions_off) THEN

  ukca_config%i_strat_lbc_source = i_strat_lbc_off

  IF (PRESENT(i_strat_lbc_source))                                             &
    ukca_config%i_strat_lbc_source = i_strat_lbc_source

END IF

! -- UKCA temporary logicals -------------------------------------------

IF (ukca_config%i_ukca_chem /= i_ukca_chem_off) THEN

  ukca_config%l_fix_ukca_cloud_frac = .TRUE.

  ! Removing this fix is not applicable when building without UM infrastructure
  IF (PRESENT(l_fix_ukca_cloud_frac) .AND. l_um_infrastructure)                &
    ukca_config%l_fix_ukca_cloud_frac = l_fix_ukca_cloud_frac

END IF

! Temporary logicals for interactive dry deposition scheme

IF (ukca_config%l_ukca_intdd) THEN

  ukca_config%l_fix_improve_drydep = .TRUE.
  ukca_config%l_fix_ukca_h2dd_x = .TRUE.

  IF (PRESENT(l_fix_improve_drydep))                                           &
    ukca_config%l_fix_improve_drydep = l_fix_improve_drydep
  IF (PRESENT(l_fix_ukca_h2dd_x))                                              &
    ukca_config%l_fix_ukca_h2dd_x = l_fix_ukca_h2dd_x

  IF ((ukca_config%ntype == 13 .OR. ukca_config%ntype == 17 .OR.               &
       ukca_config%ntype == 27) .AND.                                          &
      (.NOT. ukca_config%l_ukca_dry_dep_so2wet)) THEN
    ukca_config%l_fix_drydep_so2_water = .TRUE.
    IF (PRESENT(l_fix_drydep_so2_water))                                       &
      ukca_config%l_fix_drydep_so2_water = l_fix_drydep_so2_water
  END IF

END IF

! Temporary logicals for specific chemistry schemes

IF (ukca_config%i_ukca_chem == i_ukca_chem_offline_be) THEN
  ukca_config%l_fix_ukca_offox_h2o_fac = .TRUE.
  IF (PRESENT(l_fix_ukca_offox_h2o_fac))                                       &
    ukca_config%l_fix_ukca_offox_h2o_fac = l_fix_ukca_offox_h2o_fac
END IF

IF (ukca_config%l_ukca_mode .AND. l_nr_scheme_selected) THEN
  ukca_config%l_fix_ukca_h2so4_ystore = .TRUE.
  IF (PRESENT(l_fix_ukca_h2so4_ystore))                                        &
    ukca_config%l_fix_ukca_h2so4_ystore = l_fix_ukca_h2so4_ystore
END IF

! Settings for managing photolysis environmental driver
! requirements on behalf of external UKCA Photolysis code

IF (ukca_config%i_ukca_chem /= i_ukca_chem_off) THEN

  ukca_config%i_photol_scheme = 0
  ukca_config%i_photol_scheme_off = 0

  IF (PRESENT(i_photol_scheme)) THEN

    ukca_config%i_photol_scheme = i_photol_scheme

    IF (PRESENT(i_photol_scheme_off)) THEN
      ukca_config%i_photol_scheme_off = i_photol_scheme_off
    ELSE
      error_code_ptr = errcode_value_missing
    END IF

    IF (PRESENT(i_photol_scheme_strat_only)) THEN
      ukca_config%i_photol_scheme_strat_only = i_photol_scheme_strat_only
    ELSE
      error_code_ptr = errcode_value_missing
    END IF

    IF (PRESENT(i_photol_scheme_2d)) THEN
      ukca_config%i_photol_scheme_2d = i_photol_scheme_2d
    ELSE
      error_code_ptr = errcode_value_missing
    END IF

    IF (PRESENT(i_photol_scheme_fastjx)) THEN
      ukca_config%i_photol_scheme_fastjx = i_photol_scheme_fastjx
    ELSE
      error_code_ptr = errcode_value_missing
    END IF

    IF (error_code_ptr > 0) THEN
      IF (PRESENT(error_message))                                              &
        error_message =                                                        &
          'Missing one or more option codes for interpreting i_photol_scheme'
      IF (PRESENT(error_routine)) error_routine = RoutineName
      IF (lhook)                                                               &
        CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF

  END IF

END IF

! -- GLOMAP-mode configuration ----------------------------------------

IF (ukca_config%l_ukca_mode) THEN
  IF (ukca_config%i_ukca_chem /= i_ukca_chem_off) THEN

    ! -- General GLOMAP configuration options ---------------------------

    glomap_config%ukca_mode_seg_size = ukca_config%row_length *                &
                                       ukca_config%rows
    glomap_config%l_mode_bhn_on = .TRUE.

    IF (PRESENT(i_mode_nzts)) glomap_config%i_mode_nzts = i_mode_nzts
    IF (PRESENT(ukca_mode_seg_size))                                           &
      glomap_config%ukca_mode_seg_size = ukca_mode_seg_size
    IF (PRESENT(i_mode_setup)) glomap_config%i_mode_setup = i_mode_setup
    IF (PRESENT(l_mode_bhn_on)) glomap_config%l_mode_bhn_on = l_mode_bhn_on
    IF (PRESENT(mode_activation_dryr))                                         &
      glomap_config%mode_activation_dryr = mode_activation_dryr
    IF (PRESENT(l_dust_ageing_on))                                             &
      glomap_config%l_dust_ageing_on = l_dust_ageing_on

    IF (ukca_config%l_ukca_scale_ppe) THEN

      IF (PRESENT(dry_depvel_acc_scaling))                                     &
        glomap_config%dry_depvel_acc_scaling = dry_depvel_acc_scaling

      IF (PRESENT(acc_cor_scav_scaling))                                       &
        glomap_config%acc_cor_scav_scaling = acc_cor_scav_scaling

    END IF

    ! -- GLOMAP deposition configuration options --

    ! Dry deposition

    IF (.NOT. (ukca_config%l_ukca_drydep_off)) THEN

      glomap_config%l_ddepaer = .TRUE.

      IF (PRESENT(l_ddepaer)) glomap_config%l_ddepaer = l_ddepaer

    END IF

    ! Wet deposition

    ! Fraction of in-cloud oxidised SO2 removed by precip. is zero by default
    ! and is always zero if wet deposition is off. Note that a value of 0.25
    ! is recommended.
    glomap_config%mode_incld_so2_rfrac = 0.0

    IF (.NOT. (ukca_config%l_ukca_wetdep_off)) THEN

      glomap_config%l_aero_rainout = .TRUE.
      glomap_config%l_impc_scav = .TRUE.

      IF (PRESENT(mode_incld_so2_rfrac))                                       &
        glomap_config%mode_incld_so2_rfrac = mode_incld_so2_rfrac
      IF (PRESENT(l_aero_rainout)) glomap_config%l_aero_rainout = l_aero_rainout
      IF (PRESENT(l_impc_scav)) glomap_config%l_impc_scav = l_impc_scav

      ! Nucleation scavenging options
      glomap_config%i_mode_nucscav = 3
      glomap_config%l_cv_rainout = .TRUE.

      IF (PRESENT(l_cv_rainout)) glomap_config%l_cv_rainout = l_cv_rainout
      IF (PRESENT(i_mode_nucscav)) glomap_config%i_mode_nucscav = i_mode_nucscav

      ! Impaction scavenging options
      IF (glomap_config%l_impc_scav) THEN

        glomap_config%l_dust_slinn_impc_scav = .TRUE.

        IF (PRESENT(l_dust_slinn_impc_scav))                                   &
          glomap_config%l_dust_slinn_impc_scav = l_dust_slinn_impc_scav

      END IF

    END IF

    ! Boundary layer nucleation options

    IF (glomap_config%l_mode_bhn_on) THEN

      glomap_config%l_mode_bln_on = .TRUE.
      glomap_config%i_mode_bln_param_method = 1

      IF (PRESENT(l_mode_bln_on)) glomap_config%l_mode_bln_on = l_mode_bln_on
      IF (PRESENT(i_mode_bln_param_method))                                    &
        glomap_config%i_mode_bln_param_method = i_mode_bln_param_method

    END IF

    ! -- GLOMAP emissions configuration options --------------------------

    IF (.NOT. ukca_config%l_ukca_emissions_off) THEN

      IF (PRESENT(l_ukca_primss)) glomap_config%l_ukca_primss = l_ukca_primss
      IF (PRESENT(l_ukca_primsu)) glomap_config%l_ukca_primsu = l_ukca_primsu
      IF (PRESENT(l_ukca_primdu)) glomap_config%l_ukca_primdu = l_ukca_primdu
      IF (PRESENT(l_ukca_primbcoc))                                            &
        glomap_config%l_ukca_primbcoc = l_ukca_primbcoc

      ! Carbonaceous emissions configuration
      IF (glomap_config%l_ukca_primbcoc) THEN
        IF (PRESENT(l_ukca_prim_moc))                                          &
          glomap_config%l_ukca_prim_moc = l_ukca_prim_moc
        IF (PRESENT(l_bcoc_bf)) glomap_config%l_bcoc_bf = l_bcoc_bf
        IF (PRESENT(l_bcoc_bm)) glomap_config%l_bcoc_bm = l_bcoc_bm
        IF (PRESENT(l_bcoc_ff)) glomap_config%l_bcoc_ff = l_bcoc_ff
        IF (PRESENT(l_ukca_scale_biom_aer_ems))                                &
          glomap_config%l_ukca_scale_biom_aer_ems = l_ukca_scale_biom_aer_ems
      END IF

      ! Sea salt scaling
      IF (glomap_config%l_ukca_primss) THEN
        IF (PRESENT(l_ukca_scale_sea_salt_ems))                                &
          glomap_config%l_ukca_scale_sea_salt_ems = l_ukca_scale_sea_salt_ems
        IF (PRESENT(i_primss_method))                                          &
          glomap_config%i_primss_method = i_primss_method
      END IF
      ! Marine POM scaling
      IF (glomap_config%l_ukca_prim_moc .AND.                                  &
          glomap_config%l_ukca_primbcoc .AND.                                  &
          glomap_config%l_ukca_primss) THEN
        IF (PRESENT(l_ukca_scale_marine_pom_ems))                              &
          glomap_config%l_ukca_scale_marine_pom_ems =                          &
          l_ukca_scale_marine_pom_ems
      END IF

      ! Scaling factor for biomass burning emissions
      IF (glomap_config%l_ukca_scale_biom_aer_ems) THEN
        IF (PRESENT(biom_aer_ems_scaling))                                     &
          glomap_config%biom_aer_ems_scaling = biom_aer_ems_scaling
      END IF
      ! Scaling factor for sea-salt emissions
      IF (glomap_config%l_ukca_scale_sea_salt_ems) THEN
        IF (PRESENT(sea_salt_ems_scaling))                                     &
          glomap_config%sea_salt_ems_scaling = sea_salt_ems_scaling
      END IF
      ! Scaling factor for POM emissions
      IF (glomap_config%l_ukca_scale_marine_pom_ems) THEN
        IF (PRESENT(marine_pom_ems_scaling))                                   &
          glomap_config%marine_pom_ems_scaling = marine_pom_ems_scaling
      END IF

    END IF

    ! Nitrate emissions configuration
    IF (PRESENT(l_ukca_fine_no3_prod))                                         &
      glomap_config%l_ukca_fine_no3_prod = l_ukca_fine_no3_prod
    IF (PRESENT(l_ukca_coarse_no3_prod))                                       &
      glomap_config%l_ukca_coarse_no3_prod = l_ukca_coarse_no3_prod
    IF (PRESENT(l_no3_prod_in_aero_step))                                      &
      glomap_config%l_no3_prod_in_aero_step = l_no3_prod_in_aero_step
    IF (PRESENT(hno3_uptake_coeff))                                            &
      glomap_config%hno3_uptake_coeff = hno3_uptake_coeff

    ! -- GLOMAP feedback configuration options ---------------------------

    IF (PRESENT(l_ukca_radaer)) glomap_config%l_ukca_radaer = l_ukca_radaer
    IF (glomap_config%l_ukca_radaer) THEN
      IF (PRESENT(i_ukca_tune_bc)) glomap_config%i_ukca_tune_bc =              &
                                   i_ukca_tune_bc
    END IF

    glomap_config%i_ukca_activation_scheme = i_ukca_activation_off
    IF (PRESENT(i_ukca_activation_scheme))                                     &
    glomap_config%i_ukca_activation_scheme = i_ukca_activation_scheme

    ! Activation scheme configuration

    IF (glomap_config%i_ukca_activation_scheme /= i_ukca_activation_off) THEN

      IF (PRESENT(l_ntpreq_n_activ_sum))                                       &
        glomap_config%l_ntpreq_n_activ_sum = l_ntpreq_n_activ_sum
      IF (PRESENT(l_ntpreq_dryd_nuc_sol))                                      &
        glomap_config%l_ntpreq_dryd_nuc_sol = l_ntpreq_dryd_nuc_sol

      ! Configuration specific to UKCA Activate

      IF (glomap_config%i_ukca_activation_scheme == i_ukca_activation_arg) THEN

        glomap_config%sigwmin = 0.01

        IF (PRESENT(i_ukca_nwbins)) glomap_config%i_ukca_nwbins = i_ukca_nwbins
        IF (PRESENT(sigwmin)) glomap_config%sigwmin = sigwmin
        IF (PRESENT(l_ukca_sfix)) glomap_config%l_ukca_sfix = l_ukca_sfix

      END IF

    END IF

    IF (ukca_config%l_ukca_scale_ppe) THEN
      IF (PRESENT(sigma_updraught_scaling))                                    &
        glomap_config%sigma_updraught_scaling = sigma_updraught_scaling
    END IF

    ! -- GLOMAP temporary logicals ---------------------------------------

    glomap_config%l_fix_neg_pvol_wat = .TRUE.
    glomap_config%l_fix_nacl_density = .TRUE.

    IF (PRESENT(l_fix_neg_pvol_wat))                                           &
      glomap_config%l_fix_neg_pvol_wat = l_fix_neg_pvol_wat

    IF (glomap_config%l_impc_scav) THEN

      glomap_config%l_fix_ukca_impscav = .TRUE.

      IF (PRESENT(l_fix_ukca_impscav))                                         &
        glomap_config%l_fix_ukca_impscav = l_fix_ukca_impscav

    END IF

    IF (PRESENT(l_fix_nacl_density))                                           &
      glomap_config%l_fix_nacl_density = l_fix_nacl_density

    IF (glomap_config%l_ddepaer) THEN
      IF (PRESENT(l_improve_aero_drydep))                                      &
        glomap_config%l_improve_aero_drydep = l_improve_aero_drydep
    END IF

    IF (PRESENT(l_fix_ukca_water_content))                                     &
      glomap_config%l_fix_ukca_water_content = l_fix_ukca_water_content

    ! Temporary logicals for UKCA Activate scheme

    IF (glomap_config%i_ukca_activation_scheme == i_ukca_activation_arg) THEN

      glomap_config%l_fix_ukca_activate_pdf = .TRUE.
      glomap_config%l_fix_ukca_activate_vert_rep = .TRUE.

      IF (PRESENT(l_fix_ukca_activate_pdf))                                    &
        glomap_config%l_fix_ukca_activate_pdf = l_fix_ukca_activate_pdf
      IF (PRESENT(l_fix_ukca_activate_vert_rep))                               &
        glomap_config%l_fix_ukca_activate_vert_rep =                           &
        l_fix_ukca_activate_vert_rep
      IF (PRESENT(l_bug_repro_tke_index))                                      &
        glomap_config%l_bug_repro_tke_index = l_bug_repro_tke_index

    END IF

    IF (PRESENT(l_fix_ukca_hygroscopicities))                                  &
        glomap_config%l_fix_ukca_hygroscopicities = l_fix_ukca_hygroscopicities

  ELSE
    !
    !  Below is the standalone GLOMAP-mode setup without chemistry
    !  including the wider UKCA and environment options needed to run
    !  dust-only
    !
    ! UKCA settings - emissions, environment, etc
    !
    ukca_config%bl_levels = 1
    IF (PRESENT(bl_levels)) ukca_config%bl_levels = bl_levels
    IF (PRESENT(chem_timestep)) ukca_config%chem_timestep = chem_timestep
    ukca_config%l_ukca_chem_aero = .FALSE.
    IF (PRESENT(l_fix_tropopause_level))                                       &
      ukca_config%l_fix_tropopause_level = l_fix_tropopause_level
    IF (ukca_config%l_fix_tropopause_level) THEN
      ukca_config%fixed_tropopause_level = 1
      IF (PRESENT(fixed_tropopause_level))                                     &
        ukca_config%fixed_tropopause_level = fixed_tropopause_level
    END IF
    IF (PRESENT(l_ukca_emissions_off))                                         &
      ukca_config%l_ukca_emissions_off = l_ukca_emissions_off
    IF (PRESENT(l_ukca_drydep_off))                                            &
      ukca_config%l_ukca_drydep_off = l_ukca_drydep_off
    IF (PRESENT(l_ukca_wetdep_off))                                            &
      ukca_config%l_ukca_wetdep_off = l_ukca_wetdep_off
    IF (PRESENT(l_support_ems_gridbox_units))                                  &
      ukca_config%l_support_ems_gridbox_units = l_support_ems_gridbox_units
    IF (PRESENT(l_suppress_ems))                                               &
      ukca_config%l_suppress_ems = l_suppress_ems
    ukca_config%i_ukca_light_param = i_light_param_off
    IF (PRESENT(l_use_gridbox_volume))                                         &
      ukca_config%l_use_gridbox_volume = l_use_gridbox_volume
    ukca_config%nlev_ent_tr_mix = 0
    ukca_config%i_ukca_dms_flux = i_dms_flux_off
    IF (.NOT. ukca_config%l_suppress_ems .AND. PRESENT(nlev_ent_tr_mix))       &
      ukca_config%nlev_ent_tr_mix = nlev_ent_tr_mix
    IF (.NOT. ukca_config%l_ukca_wetdep_off) THEN
      IF (PRESENT(l_param_conv)) ukca_config%l_param_conv = l_param_conv
    END IF
    IF (PRESENT(l_use_gridbox_mass))                                           &
       ukca_config%l_use_gridbox_mass = l_use_gridbox_mass
    ukca_config%l_fix_ukca_cloud_frac = .TRUE.
    IF (PRESENT(l_fix_ukca_cloud_frac))                                        &
      ukca_config%l_fix_ukca_cloud_frac = l_fix_ukca_cloud_frac

    IF (PRESENT(l_environ_rel_humid))                                          &
      ukca_config%l_environ_rel_humid = l_environ_rel_humid

    IF (ukca_config%l_ukca_scale_ppe) THEN

      IF (PRESENT(dry_depvel_acc_scaling))                                     &
        glomap_config%dry_depvel_acc_scaling = dry_depvel_acc_scaling

      IF (PRESENT(acc_cor_scav_scaling))                                       &
        glomap_config%acc_cor_scav_scaling = acc_cor_scav_scaling

    END IF


    !
    ! GLOMAP settings
    !
    glomap_config%ukca_mode_seg_size = ukca_config%row_length * ukca_config%rows
    glomap_config%l_mode_bhn_on = .FALSE.
    glomap_config%l_mode_bln_on = .FALSE.
    IF (PRESENT(i_mode_nzts)) glomap_config%i_mode_nzts = i_mode_nzts



    IF (PRESENT(ukca_mode_seg_size))                                           &
      glomap_config%ukca_mode_seg_size = ukca_mode_seg_size
    IF (PRESENT(i_mode_setup)) glomap_config%i_mode_setup = i_mode_setup
    IF (PRESENT(l_dust_ageing_on))                                             &
      glomap_config%l_dust_ageing_on = l_dust_ageing_on
    IF (.NOT. (ukca_config%l_ukca_drydep_off)) THEN
      glomap_config%l_ddepaer = .TRUE.
      IF (PRESENT(l_ddepaer)) glomap_config%l_ddepaer = l_ddepaer
    END IF
    IF (.NOT. (ukca_config%l_ukca_wetdep_off)) THEN
      glomap_config%l_aero_rainout = .TRUE.
      glomap_config%l_impc_scav = .TRUE.
      IF (PRESENT(mode_incld_so2_rfrac))                                       &
        glomap_config%mode_incld_so2_rfrac = mode_incld_so2_rfrac
      IF (PRESENT(l_aero_rainout)) glomap_config%l_aero_rainout = l_aero_rainout
      IF (PRESENT(l_impc_scav)) glomap_config%l_impc_scav = l_impc_scav
      ! Nucleation scavenging options
      glomap_config%i_mode_nucscav = 3
      glomap_config%l_cv_rainout = .TRUE.
      IF (PRESENT(l_cv_rainout)) glomap_config%l_cv_rainout = l_cv_rainout
      IF (PRESENT(i_mode_nucscav)) glomap_config%i_mode_nucscav = i_mode_nucscav
      ! Impaction scavenging options
      IF (glomap_config%l_impc_scav) THEN
        glomap_config%l_dust_slinn_impc_scav = .TRUE.
        IF (PRESENT(l_dust_slinn_impc_scav))                                   &
          glomap_config%l_dust_slinn_impc_scav = l_dust_slinn_impc_scav
      END IF
    END IF
    IF (.NOT. ukca_config%l_ukca_emissions_off) THEN
      IF (PRESENT(l_ukca_primdu)) glomap_config%l_ukca_primdu = l_ukca_primdu
    END IF
    IF (PRESENT(l_ukca_radaer)) glomap_config%l_ukca_radaer = l_ukca_radaer
    IF (glomap_config%l_ukca_radaer) THEN
      IF (PRESENT(i_ukca_tune_bc)) glomap_config%i_ukca_tune_bc =              &
                                   i_ukca_tune_bc
    END IF
    glomap_config%i_ukca_activation_scheme = i_ukca_activation_off

    IF (ukca_config%l_ukca_scale_ppe) THEN
      IF (PRESENT(sigma_updraught_scaling))                                    &
        glomap_config%sigma_updraught_scaling = sigma_updraught_scaling
    END IF

    glomap_config%l_fix_neg_pvol_wat = .TRUE.
    glomap_config%l_fix_nacl_density = .TRUE.
    IF (PRESENT(l_fix_neg_pvol_wat))                                           &
      glomap_config%l_fix_neg_pvol_wat = l_fix_neg_pvol_wat
    IF (glomap_config%l_impc_scav) THEN
      glomap_config%l_fix_ukca_impscav = .TRUE.
      IF (PRESENT(l_fix_ukca_impscav))                                         &
        glomap_config%l_fix_ukca_impscav = l_fix_ukca_impscav
    END IF
    IF (glomap_config%l_ddepaer) THEN
      IF (PRESENT(l_improve_aero_drydep))                                      &
        glomap_config%l_improve_aero_drydep = l_improve_aero_drydep
    END IF
    IF (PRESENT(l_fix_ukca_water_content))                                     &
        glomap_config%l_fix_ukca_water_content = l_fix_ukca_water_content
    IF (.NOT. (ukca_config%l_ukca_emissions_off .OR.                           &
               ukca_config%l_suppress_ems)) THEN
      IF (PRESENT(proc_bl_tracer_mix)) bl_tracer_mix => proc_bl_tracer_mix
    END IF
  END IF  ! ukca_config%i_ukca_chem /= i_ukca_chem_off
END IF  ! l_ukca_mode

! -- UKCA environmental driver options with GLOMAP dependencies --------

IF (ukca_config%l_ukca_intdd .OR. glomap_config%l_ddepaer .OR.                 &
    ukca_config%l_ukca_qch4inter .OR. ukca_config%l_ukca_ibvoc .OR.            &
    ukca_config%l_ukca_inferno .OR. glomap_config%l_ukca_primss) THEN
  IF (PRESENT(l_ctile)) ukca_config%l_ctile = l_ctile
END IF

! -- Parent callback procedures ----------------------------------------

IF (ukca_config%i_ukca_chem /= i_ukca_chem_off .AND.                           &
    .NOT. (ukca_config%l_ukca_emissions_off .OR.                               &
          ukca_config%l_suppress_ems)) THEN
  IF (PRESENT(proc_bl_tracer_mix)) bl_tracer_mix => proc_bl_tracer_mix
END IF

IF (PRESENT(proc_calc_ozonecol)) calc_ozonecol => proc_calc_ozonecol
IF (PRESENT(proc_diag2d_copy_out)) diag2d_copy_out => proc_diag2d_copy_out
IF (PRESENT(proc_diag3d_copy_out)) diag3d_copy_out => proc_diag3d_copy_out

! ----------------------------------------------------------------------

! Check UKCA logicals are consistent and set internal UKCA configuration values
! based on configuration data provided
CALL ukca_init()

! -- Dust scheme needed for nitrate emissions --

IF (ukca_config%l_ukca_mode .AND. .NOT. ukca_config%l_ukca_emissions_off) THEN

  glomap_config%l_6bin_dust_no3 =                                              &
    ( ( glomap_config%l_ukca_coarse_no3_prod .OR.                              &
        glomap_config%l_ukca_fine_no3_prod ) .AND.                             &
      l_dust .AND. ( .NOT. l_twobin_dust ) .AND.                               &
      ( .NOT. glomap_config%l_ukca_primdu ) .AND.                              &
      ( .NOT. glomap_config%l_no3_prod_in_aero_step) )

  glomap_config%l_2bin_dust_no3 =                                              &
    ( ( glomap_config%l_ukca_coarse_no3_prod .OR.                              &
        glomap_config%l_ukca_fine_no3_prod ) .AND.                             &
      l_dust .AND. l_twobin_dust .AND.                                         &
      ( .NOT. glomap_config%l_ukca_primdu ) .AND.                              &
      ( .NOT. glomap_config%l_no3_prod_in_aero_step) )

  IF ( glomap_config%l_6bin_dust_no3 ) THEN
    glomap_config%i_dust_scheme=1                      ! CLASSIC 6 bin dust
  ELSE IF ( glomap_config%l_2bin_dust_no3 ) THEN
    glomap_config%i_dust_scheme=2                      ! CLASSIC 2 bin dust
  ELSE IF ( glomap_config%l_ukca_primdu ) THEN
    glomap_config%i_dust_scheme=3                      ! GLOMAP-mode dust
  ELSE
    glomap_config%i_dust_scheme=-1                     ! No dust
  END IF

END IF

! Initialise chemical definition arrays
CALL ukca_chem1_init()

! Set up lists of tracer and non-transported prognostic species for chemistry
CALL asad_mod_pre_setup_init()
CALL asad_inrats_set_sp_lists()

! Determine which species are specified as 'constants' in the scheme and need
! values to be passed from the parent
IF ( ukca_config%l_chem_environ_gas_scalars ) THEN
  DO i = 1, SIZE(speci)
    SELECT CASE(speci(i))
    CASE ('CH4')
      ukca_config%l_chem_environ_ch4_scalar = (ctype(i) == 'CT')
    CASE ('CO2')
      ukca_config%l_chem_environ_co2_scalar = (ctype(i) == 'CT')
    CASE ('H2')
      ukca_config%l_chem_environ_h2_scalar  = (ctype(i) == 'CT')
    CASE ('N2')
      ukca_config%l_chem_environ_n2_scalar  = (ctype(i) == 'CT')
    CASE ('O2')
      ukca_config%l_chem_environ_o2_scalar  = (ctype(i) == 'CT')
    END SELECT
  END DO
END IF

! Specify tracer requirement based on details of the configuration
CALL init_tracer_req( ukca_config,                                             &
                      advt, error_code_ptr,                                    &
                      error_message=error_message,                             &
                      error_routine=error_routine)
IF (error_code_ptr > 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Set up structure holding all NTPs and specify NTP requirement
! based on details of the configuration
CALL ntp_init()

! Set emitted species lists and various array sizes according to which
! chemistry scheme is selected
CALL ukca_set_config_defs()

! Set configuration details that depend on emissions settings
IF (ukca_config%i_ukca_chem /= i_ukca_chem_off .AND.                           &
    .NOT. ukca_config%l_ukca_emissions_off) THEN
  ukca_config%l_diurnal_isopems =                                              &
    l_ukca_diurnal_isopems .AND. ANY(em_chem_spec == 'C5H8')
  ukca_config%l_seawater_dms =                                                 &
    (ukca_config%i_ukca_dms_flux /= i_dms_flux_off) .AND.                      &
    ANY(em_chem_spec == 'DMS')
END IF

! Initialise master diagnostics list and determine availability of
! each diagnostic given the configuration.
CALL init_diagnostics(error_code_ptr, ukca_config, advt,                       &
                      error_message=error_message,                             &
                      error_routine=error_routine)
IF (error_code_ptr > 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Specify environment field requirement based on details of the configuration
CALL init_environment_req(ukca_config, glomap_config,                          &
                          speci, advt, lbc_spec, cfc_lumped, em_chem_spec,     &
                          ctype, error_code_ptr, error_message=error_message,  &
                          error_routine=error_routine)
IF (error_code_ptr > 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Ensure environment fields and associated data are in an uninitialised state
CALL clear_environment_fields()

! Set flag to show that a valid UKCA configuration is set up
l_ukca_config_available = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE ukca_setup

END MODULE ukca_setup_mod
