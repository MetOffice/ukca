! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Application program interface (API) module for UKCA.
!
! Method:
!
!  This module provides access to subroutines and parameters
!  required by a parent model for running UKCA.
!  It acts as a collation point for components of the API defined
!  in other UKCA modules rather than including any definitions itself.
!
!  Note that 'ukca_api_mod' should be the only UKCA module used by a parent
!  application. Where exceptions exist in UM legacy code, these are considered
!  deprecated and may not be robust to future UKCA developments.
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

MODULE ukca_api_mod

! Procedures and parameters made available below constitute the formal
! UKCA API. All names made available begin with 'ukca_' for clarity.
! This is important to avoid pollution of a parent application's namespace.

USE ukca_setup_mod, ONLY: ukca_setup
USE ukca_step_control_mod, ONLY: ukca_step_control
USE ukca_step_mod, ONLY: ukca_step
USE ukca_config_specification_mod, ONLY:                                       &
  ukca_get_config,                                                             &
  ukca_error_method_abort => i_error_method_abort,                             &
  ukca_error_method_return => i_error_method_return,                           &
  ukca_error_method_warn_and_return => i_error_method_warn_and_return,         &
  ukca_chem_off => i_ukca_chem_off,                                            &
  ukca_chem_trop => i_ukca_chem_trop,                                          &
  ukca_chem_raq => i_ukca_chem_raq,                                            &
  ukca_chem_offline_be => i_ukca_chem_offline_be,                              &
  ukca_chem_tropisop => i_ukca_chem_tropisop,                                  &
  ukca_chem_strattrop => i_ukca_chem_strattrop,                                &
  ukca_chem_strat => i_ukca_chem_strat,                                        &
  ukca_chem_offline => i_ukca_chem_offline,                                    &
  ukca_chem_cristrat => i_ukca_chem_cristrat,                                  &
  ukca_age_reset_by_level => i_age_reset_by_level,                             &
  ukca_age_reset_by_height => i_age_reset_by_height,                           &
  ukca_strat_lbc_off => i_strat_lbc_off,                                       &
  ukca_strat_lbc_wmoa1 => i_strat_lbc_wmoa1,                                   &
  ukca_strat_lbc_env => i_strat_lbc_env,                                       &
  ukca_strat_lbc_rcp => i_strat_lbc_rcp,                                       &
  ukca_int_method_impact => int_method_impact,                                 &
  ukca_int_method_nr => int_method_nr,                                         &
  ukca_int_method_be => int_method_be,                                         &
  ukca_int_method_be_explicit => int_method_be_explicit,                       &
  ukca_activation_off => i_ukca_activation_off,                                &
  ukca_activation_arg => i_ukca_activation_arg,                                &
  ukca_activation_jones => i_ukca_activation_jones,                            &
  ukca_light_param_off => i_light_param_off,                                   &
  ukca_light_param_pr => i_light_param_pr,                                     &
  ukca_light_param_luhar => i_light_param_luhar,                               &
  ukca_light_param_ext => i_light_param_ext,                                   &
  ukca_top_none => i_top_none,                                                 &
  ukca_top_2levh20 => i_top_2levh2o,                                           &
  ukca_top_1lev => i_top_1lev,                                                 &
  ukca_top_bc => i_top_bc,                                                     &
  ukca_top_bc_h2o => i_top_bc_h2o,                                             &
  ukca_dms_flux_off => i_dms_flux_off,                                         &
  ukca_liss_merlivat => i_liss_merlivat,                                       &
  ukca_wanninkhof => i_wanninkhof,                                             &
  ukca_nightingale => i_nightingale,                                           &
  ukca_i_primss_method_smith => i_primss_method_smith,                         &
  ukca_i_primss_method_monahan => i_primss_method_monahan,                     &
  ukca_i_primss_method_combined => i_primss_method_combined,                   &
  ukca_i_primss_method_jaegle => i_primss_method_jaegle
USE ukca_constants_setup_mod, ONLY: ukca_constants_setup
USE ukca_tracers_mod, ONLY: ukca_get_tracer_varlist
USE ukca_ntp_mod, ONLY: ukca_get_ntp_varlist
USE ukca_chem_defs_mod, ONLY:                                                  &
  ukca_get_photol_reaction_data,                                               &
  ukca_photol_varname_len => photol_varname_len
USE ukca_environment_req_mod, ONLY:                                            &
  ukca_get_environment_varlist,                                                &
  ukca_get_envgroup_varlists
USE ukca_environment_mod, ONLY: ukca_set_environment
USE ukca_fieldname_mod, ONLY:                                                  &
  ukca_maxlen_fieldname => maxlen_fieldname,                                   &
  ukca_maxlen_diagname => maxlen_diagname,                                     &
  ukca_diagname_jrate_no2 => diagname_jrate_no2,                               &
  ukca_diagname_jrate_o3a => diagname_jrate_o3a,                               &
  ukca_diagname_jrate_o3b => diagname_jrate_o3b,                               &
  ukca_diagname_jrate_o2b => diagname_jrate_o2b,                               &
  ukca_diagname_rxnflux_oh_ch4_trop => diagname_rxnflux_oh_ch4_trop,           &
  ukca_diagname_p_tropopause => diagname_p_tropopause,                         &
  ukca_diagname_o3_column_du => diagname_o3_column_du,                         &
  ukca_diagname_plumeria_height => diagname_plumeria_height
USE ukca_emiss_api_mod, ONLY:                                                  &
  ukca_get_emission_varlist,                                                   &
  ukca_register_emission,                                                      &
  ukca_set_emission
USE ukca_emiss_struct_mod, ONLY:                                               &
  ukca_maxlen_emiss_var_name => maxlen_emiss_var_name,                         &
  ukca_maxlen_emiss_tracer_name => maxlen_emiss_tracer_name,                   &
  ukca_maxlen_emiss_std_name => maxlen_emiss_std_name,                         &
  ukca_maxlen_emiss_long_name => maxlen_emiss_long_name,                       &
  ukca_maxlen_emiss_units => maxlen_emiss_units,                               &
  ukca_maxlen_emiss_hourly_fact => maxlen_emiss_hourly_fact,                   &
  ukca_maxlen_emiss_daily_fact => maxlen_emiss_daily_fact,                     &
  ukca_maxlen_emiss_vert_fact => maxlen_emiss_vert_fact
USE ukca_diagnostics_type_mod, ONLY:                                           &
  ukca_diag_status_inactive => diag_status_inactive,                           &
  ukca_diag_status_requested => diag_status_requested,                         &
  ukca_diag_status_valid => diag_status_valid,                                 &
  ukca_diag_status_skipped => diag_status_skipped,                             &
  ukca_diag_status_unavailable => diag_status_unavailable
USE ukca_diagnostics_init_mod, ONLY: ukca_get_available_diag_varlist
USE ukca_diagnostics_requests_mod, ONLY:                                       &
  ukca_set_diagnostic_requests,                                                &
  ukca_update_diagnostic_requests,                                             &
  ukca_get_diagnostic_request_info
USE ukca_error_mod, ONLY:                                                      &
  ukca_maxlen_message => maxlen_message,                                       &
  ukca_maxlen_procname => maxlen_procname,                                     &
  ukca_errcode_ukca_uninit => errcode_ukca_uninit,                             &
  ukca_errcode_tracer_req_uninit => errcode_tracer_req_uninit,                 &
  ukca_errcode_tracer_mismatch => errcode_tracer_mismatch,                     &
  ukca_errcode_ntp_uninit => errcode_ntp_uninit,                               &
  ukca_errcode_ntp_mismatch => errcode_ntp_mismatch,                           &
  ukca_errcode_env_req_uninit => errcode_env_req_uninit,                       &
  ukca_errcode_env_field_unknown => errcode_env_field_unknown,                 &
  ukca_errcode_env_field_mismatch => errcode_env_field_mismatch,               &
  ukca_errcode_env_field_missing => errcode_env_field_missing,                 &
  ukca_errcode_diag_req_unknown => errcode_diag_req_unknown,                   &
  ukca_errcode_diag_req_duplicate => errcode_diag_req_duplicate,               &
  ukca_errcode_diag_req_mismatch => errcode_diag_mismatch,                     &
  ukca_errcode_ukca_internal_fault => errcode_ukca_internal_fault,             &
  ukca_errcode_value_unknown => errcode_value_unknown,                         &
  ukca_errcode_value_invalid => errcode_value_invalid,                         &
  ukca_errcode_value_missing => errcode_value_missing,                         &
  ukca_errcode_unexpected_api_call => errcode_unexpected_api_call
USE ukca_ddepaer_coeff_mod, ONLY:                                                 ukca_zhg_eg_nedleaf => zhg_eg_nedleaf, &
  ukca_zhg_eg_brdleaf => zhg_eg_brdleaf,                                       &
  ukca_zhg_dec_nedleaf => zhg_dec_nedleaf,                                     &
  ukca_zhg_dec_brdleaf => zhg_dec_brdleaf,                                     &
  ukca_zhg_mix_brdned_leaf => zhg_mix_brdned_leaf,                             &
  ukca_zhg_grass => zhg_grass,                                                 &
  ukca_zhg_crop => zhg_crop,                                                   &
  ukca_zhg_desert => zhg_desert,                                               &
  ukca_zhg_tundra => zhg_tundra,                                               &
  ukca_zhg_shrub => zhg_shrub,                                                 &
  ukca_zhg_wetl_veg => zhg_wetl_veg,                                           &
  ukca_zhg_ice => zhg_ice,                                                     &
  ukca_zhg_inl_water => zhg_inl_water,                                         &
  ukca_zhg_ocean => zhg_ocean,                                                 &
  ukca_zhg_urban => zhg_urban,                                                 &
  ukca_zhg_ned_leaf => zhg_ned_leaf,                                           &
  ukca_zhg_brd_leaf => zhg_brd_leaf

IMPLICIT NONE

END MODULE ukca_api_mod
