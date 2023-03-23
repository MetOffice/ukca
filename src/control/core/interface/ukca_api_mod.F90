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
!  The UKCA API is currently under development. Once completed, this
!  should be the only UKCA module to be used by a parent application.
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
  ukca_photolysis_off => i_ukca_nophot,                                        &
  ukca_photolysis_strat_only => i_ukca_photol_strat,                           &
  ukca_photolysis_2d => i_ukca_phot2d,                                         &
  ukca_photolysis_fastjx => i_ukca_fastjx,                                     &
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
  ukca_nightingale => i_nightingale

USE ukca_tracers_mod, ONLY: ukca_get_tracer_varlist
USE ukca_ntp_mod, ONLY: ukca_get_ntp_varlist
USE ukca_chem_defs_mod, ONLY: ukca_get_photol_reaction_data,                   &
  ukca_photol_varname_len => photol_varname_len
USE ukca_environment_req_mod, ONLY: ukca_get_environment_varlist,              &
                                    ukca_get_envgroup_varlists
USE ukca_environment_mod, ONLY: ukca_set_environment
USE ukca_fieldname_mod, ONLY: ukca_maxlen_fieldname => maxlen_fieldname
USE ukca_emiss_api_mod, ONLY: ukca_get_emission_varlist,                       &
                              ukca_register_emission, ukca_set_emission
USE ukca_emiss_struct_mod, ONLY:                                               &
  ukca_maxlen_emiss_var_name => maxlen_emiss_var_name,                         &
  ukca_maxlen_emiss_tracer_name => maxlen_emiss_tracer_name,                   &
  ukca_maxlen_emiss_std_name => maxlen_emiss_std_name,                         &
  ukca_maxlen_emiss_long_name => maxlen_emiss_long_name,                       &
  ukca_maxlen_emiss_units => maxlen_emiss_units,                               &
  ukca_maxlen_emiss_hourly_fact => maxlen_emiss_hourly_fact,                   &
  ukca_maxlen_emiss_daily_fact => maxlen_emiss_daily_fact,                     &
  ukca_maxlen_emiss_vert_fact => maxlen_emiss_vert_fact
USE ukca_error_mod, ONLY: ukca_maxlen_message => maxlen_message,               &
                          ukca_maxlen_procname => maxlen_procname
USE ukca_ddepaer_coeff_mod, ONLY: ukca_zhg_eg_nedleaf => zhg_eg_nedleaf,       &
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
