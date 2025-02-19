! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Application program interface (API) module for UKCA Photolysis.
!
! Method:
!
!  This module provides access to subroutines and parameters
!  required by a parent or component model for running Photolysis.
!  It acts as a collation point for components of the API defined
!  in other Photolysis modules rather than including any definitions
!  itself.
!
!  Note that 'photol_api_mod' should be the only Photolysis module used
!  by
!  a parent application.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_Photolysis
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE photol_api_mod

! Procedures and parameters made available below constitute the formal
! Photolysis API. All names made available begin with 'photol_' for
! clarity.
! This is important to avoid pollution of a parent application's
! namespace.

USE photol_calc_ozonecol_mod, ONLY: photol_calc_ozonecol
USE photol_setup_mod, ONLY: photol_setup

USE photol_config_specification_mod, ONLY:                                     &
  photol_get_config, photol_off => i_scheme_nophot,                            &
  photol_strat_only => i_scheme_photol_strat,                                  &
  photol_2d => i_scheme_phot2d, photol_fastjx => i_scheme_fastjx

USE photol_fieldname_mod,   ONLY:                                              &
  photol_fieldname_len =>  fieldname_len,                                      &
  photol_fldname_aod_sulph_aitk => fldname_aod_sulph_aitk,                     &
  photol_fldname_aod_sulph_accum => fldname_aod_sulph_accum,                   &
  photol_fldname_area_cloud_fraction => fldname_area_cloud_fraction,           &
  photol_fldname_conv_cloud_amount => fldname_conv_cloud_amount,               &
  photol_fldname_conv_cloud_base => fldname_conv_cloud_base,                   &
  photol_fldname_conv_cloud_lwp => fldname_conv_cloud_lwp,                     &
  photol_fldname_conv_cloud_top => fldname_conv_cloud_top,                     &
  photol_fldname_cos_latitude => fldname_cos_latitude,                         &
  photol_fldname_equation_of_time => fldname_equation_of_time,                 &
  photol_fldname_land_fraction => fldname_land_fraction,                       &
  photol_fldname_longitude => fldname_longitude,                               &
  photol_fldname_ozone_mmr => fldname_ozone_mmr,                               &
  photol_fldname_p_layer_boundaries => fldname_p_layer_boundaries,             &
  photol_fldname_p_theta_levels => fldname_p_theta_levels,                     &
  photol_fldname_photol_rates_2d => fldname_photol_rates_2d,                   &
  photol_fldname_qcf => fldname_qcf, photol_fldname_qcl => fldname_qcl,        &
  photol_fldname_rad_ctl_jo2 => fldname_rad_ctl_jo2,                           &
  photol_fldname_rad_ctl_jo2b => fldname_rad_ctl_jo2b,                         &
  photol_fldname_r_rho_levels => fldname_r_rho_levels,                         &
  photol_fldname_r_theta_levels => fldname_r_theta_levels,                     &
  photol_fldname_sec_since_midnight => fldname_sec_since_midnight,             &
  photol_fldname_sin_declination => fldname_sin_declination,                   &
  photol_fldname_sin_latitude => fldname_sin_latitude,                         &
  photol_fldname_so4_aitken => fldname_so4_aitken,                             &
  photol_fldname_so4_accum => fldname_so4_accum,                               &
  photol_fldname_surf_albedo => fldname_surf_albedo,                           &
  photol_fldname_t_theta_levels => fldname_t_theta_levels,                     &
  photol_fldname_tan_latitude => fldname_tan_latitude,                         &
  photol_fldname_z_top_of_model => fldname_z_top_of_model

USE ukca_error_mod, ONLY: photol_error_method_abort => i_error_method_abort,   &
  photol_error_method_return => i_error_method_return,                         &
  photol_error_method_warn_and_return => i_error_method_warn_and_return

! Helper routines to process 2-D photolysis data
USE ukca_phot2d,        ONLY: ukca_photin
USE photol_curve_mod,   ONLY: photol_curve

USE photol_environment_mod, ONLY: photol_get_environ_varlist
USE photol_step_control_mod, ONLY: photol_step_control

IMPLICIT NONE

END MODULE photol_api_mod
