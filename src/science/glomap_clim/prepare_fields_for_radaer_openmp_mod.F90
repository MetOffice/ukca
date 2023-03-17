! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution..
! *****************************COPYRIGHT*******************************
!
! Purpose:
!   This module exists so that local variables within the OpenMP region
!   do not need to be part of the segmentation code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE prepare_fields_for_radaer_openmp_mod

USE um_types,                        ONLY:                                     &
    real_umphys

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='PREPARE_FIELDS_FOR_RADAER_OPENMP_MOD'

PUBLIC :: prepare_fields_for_radaer_openmp

CONTAINS

SUBROUTINE prepare_fields_for_radaer_openmp( n_points, ncp,                    &
  p_theta_levels_1d, t_theta_levels_1d,                                        &
  gc_nd_nuc_sol_1d, gc_nuc_sol_su_1d, gc_nuc_sol_oc_1d,                        &
  gc_nd_ait_sol_1d, gc_ait_sol_su_1d, gc_ait_sol_bc_1d, gc_ait_sol_oc_1d,      &
  gc_nd_acc_sol_1d, gc_acc_sol_su_1d, gc_acc_sol_bc_1d, gc_acc_sol_oc_1d,      &
  gc_acc_sol_ss_1d,                                                            &
  gc_nd_cor_sol_1d, gc_cor_sol_su_1d, gc_cor_sol_bc_1d, gc_cor_sol_oc_1d,      &
  gc_cor_sol_ss_1d,                                                            &
  gc_nd_ait_ins_1d, gc_ait_ins_bc_1d, gc_ait_ins_oc_1d,                        &
  rhcrit_1d, q_1d, qcf_1d, cloud_liq_frac_1d, cloud_blk_frac_1d,               &
  drydp, wetdp, rhopar, pvol_wat, pvol )

USE glomap_clim_calc_md_mdt_nd_mod,     ONLY:                                  &
    glomap_clim_calc_md_mdt_nd

USE glomap_clim_calc_rh_frac_clear_mod, ONLY:                                  &
    glomap_clim_calc_rh_frac_clear

USE glomap_clim_option_mod,             ONLY:                                  &
    i_glomap_clim_setup

USE parkind1,                           ONLY:                                  &
    jprb,                                                                      &
    jpim

USE ukca_calc_drydiam_mod,              ONLY:                                  &
    ukca_calc_drydiam

USE ukca_config_specification_mod,      ONLY:                                  &
    glomap_variables_climatology

USE ukca_mode_setup,                    ONLY:                                  &
    nmodes

USE ukca_volume_mode_mod,               ONLY:                                  &
    ukca_volume_mode

USE yomhook,                            ONLY:                                  &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: n_points
INTEGER, INTENT(IN) :: ncp
REAL(KIND=real_umphys), INTENT(IN) :: p_theta_levels_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) :: t_theta_levels_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_nuc_sol_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nuc_sol_su_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nuc_sol_oc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_ait_sol_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_ait_sol_su_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_ait_sol_bc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_ait_sol_oc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_acc_sol_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_acc_sol_su_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_acc_sol_bc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_acc_sol_oc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_acc_sol_ss_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_cor_sol_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_cor_sol_su_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_cor_sol_bc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_cor_sol_oc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_cor_sol_ss_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_ait_ins_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_ait_ins_bc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_ait_ins_oc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::         rhcrit_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::              q_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::            qcf_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) :: cloud_liq_frac_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) :: cloud_blk_frac_1d( n_points )

! Median particle dry diameter for each mode (m)
REAL(KIND=real_umphys), INTENT(OUT) ::    drydp( n_points, nmodes )

! Avg wet diameter of size mode (m)
REAL(KIND=real_umphys), INTENT(OUT) ::    wetdp( n_points, nmodes )

! Particle density (kg/m^3) [includes H2O & insoluble cpts]
REAL(KIND=real_umphys), INTENT(OUT) ::   rhopar( n_points, nmodes )

! Partial volume of water in each mode (m3)
REAL(KIND=real_umphys), INTENT(OUT) :: pvol_wat( n_points, nmodes )

! Partial volumes of each component in each mode (m3)
REAL(KIND=real_umphys), INTENT(OUT) ::     pvol( n_points, nmodes, ncp )

! Local variables

! Aerosol ptcl no. concentration (ptcls per cc)
REAL(KIND=real_umphys) ::    nd( n_points, nmodes )

! Total median aerosol mass (molecules per ptcl)
REAL(KIND=real_umphys) ::   mdt( n_points, nmodes )

! Component median aerosol mass (molecules per ptcl)
REAL(KIND=real_umphys) ::    md( n_points, nmodes, ncp )

! Avg dry volume of size mode (cubic metres)
REAL(KIND=real_umphys) ::  dvol( n_points, nmodes )

! Clear sky relative humidity as a fraction 1-D
REAL(KIND=real_umphys) :: rh_clr_1d ( n_points )

! These fields are only required by UKCA (not climatology).
! These fields arguments passed OUT of ukca_volume_mode.
REAL(KIND=real_umphys) ::  wvol( n_points, nmodes )
REAL(KIND=real_umphys) :: mdwat( n_points, nmodes )

! Dummy fields required so that seven modes functionality can be called
REAL(KIND=real_umphys) :: dummy_acc_sol_du_1d( n_points )
REAL(KIND=real_umphys) :: dummy_cor_sol_du_1d( n_points )
REAL(KIND=real_umphys) :: dummy_nd_acc_ins_1d( n_points )
REAL(KIND=real_umphys) :: dummy_acc_ins_du_1d( n_points )
REAL(KIND=real_umphys) :: dummy_nd_cor_ins_1d( n_points )
REAL(KIND=real_umphys) :: dummy_cor_ins_du_1d( n_points )

CHARACTER(LEN=*), PARAMETER :: RoutineName='PREPARE_FIELDS_FOR_RADAER_OPENMP'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise dummy fields to zero
dummy_acc_sol_du_1d(:) = 0.0
dummy_cor_sol_du_1d(:) = 0.0
dummy_nd_acc_ins_1d(:) = 0.0
dummy_acc_ins_du_1d(:) = 0.0
dummy_nd_cor_ins_1d(:) = 0.0
dummy_cor_ins_du_1d(:) = 0.0

! Calculate clear sky relative humidity as a fraction 1-D
CALL glomap_clim_calc_rh_frac_clear( n_points,                                 &
                                     q_1d,                                     &
                                     qcf_1d,                                   &
                                     cloud_liq_frac_1d,                        &
                                     cloud_blk_frac_1d,                        &
                                     t_theta_levels_1d,                        &
                                     p_theta_levels_1d,                        &
                                     rhcrit_1d,                                &
                                     rh_clr_1d )

! Calculate fields nd, md, mdt
CALL glomap_clim_calc_md_mdt_nd( n_points, ncp, i_glomap_clim_setup,           &
  p_theta_levels_1d, t_theta_levels_1d,                                        &
  gc_nd_nuc_sol_1d, gc_nuc_sol_su_1d, gc_nuc_sol_oc_1d,                        &
  gc_nd_ait_sol_1d, gc_ait_sol_su_1d, gc_ait_sol_bc_1d, gc_ait_sol_oc_1d,      &
  gc_nd_acc_sol_1d, gc_acc_sol_su_1d, gc_acc_sol_bc_1d, gc_acc_sol_oc_1d,      &
  gc_acc_sol_ss_1d, dummy_acc_sol_du_1d,                                       &
  gc_nd_cor_sol_1d, gc_cor_sol_su_1d, gc_cor_sol_bc_1d, gc_cor_sol_oc_1d,      &
  gc_cor_sol_ss_1d, dummy_cor_sol_du_1d,                                       &
  gc_nd_ait_ins_1d, gc_ait_ins_bc_1d, gc_ait_ins_oc_1d,                        &
  dummy_nd_acc_ins_1d, dummy_acc_ins_du_1d,                                    &
  dummy_nd_cor_ins_1d, dummy_cor_ins_du_1d,                                    &
  nd, md, mdt )

! Calculate the dry diameters and volumes
CALL ukca_calc_drydiam( n_points, glomap_variables_climatology,                &
                        nd, md, mdt, drydp, dvol )

! Calculate wet diameters, densities, partial volumes...
CALL ukca_volume_mode( glomap_variables_climatology, n_points, nd, md, mdt,    &
                       rh_clr_1d, dvol, drydp,                                 &
                       t_theta_levels_1d, p_theta_levels_1d, q_1d,             &
                       mdwat, wvol, wetdp, rhopar, pvol, pvol_wat )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE prepare_fields_for_radaer_openmp
END MODULE prepare_fields_for_radaer_openmp_mod
