! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Interface between LFRic and GLOMAP_CLIM subroutines
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! Procedure:
!   1) CALL glomap_clim_calc_rh_frac_clear
!   2) CALL glomap_clim_calc_md_mdt_nd
!   3) CALL ukca_calc_drydiam
!   4) CALL ukca_volume_mode
!   5) CALL ukca_cdnc_jones
!
! ---------------------------------------------------------------------

MODULE glomap_clim_interface_mod

USE um_types,                           ONLY:                                  &
    real_umphys

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE ::                                          &
ModuleName='GLOMAP_CLIM_INTERFACE_MOD'

CONTAINS

SUBROUTINE glomap_clim_interface ( n_points, i_glomap_clim_setup,              &
  rad_this_tstep, l_radaer, act_radius,                                        &
  p_theta_levels_1d, t_theta_levels_1d,                                        &
  gc_nd_nuc_sol_1d, gc_nuc_sol_su_1d, gc_nuc_sol_om_1d,                        &
  gc_nd_ait_sol_1d, gc_ait_sol_su_1d, gc_ait_sol_bc_1d, gc_ait_sol_om_1d,      &
  gc_nd_acc_sol_1d, gc_acc_sol_su_1d, gc_acc_sol_bc_1d, gc_acc_sol_om_1d,      &
  gc_acc_sol_ss_1d, gc_acc_sol_du_1d,                                          &
  gc_nd_cor_sol_1d, gc_cor_sol_su_1d, gc_cor_sol_bc_1d, gc_cor_sol_om_1d,      &
  gc_cor_sol_ss_1d, gc_cor_sol_du_1d,                                          &
  gc_nd_ait_ins_1d, gc_ait_ins_bc_1d, gc_ait_ins_om_1d,                        &
  gc_nd_acc_ins_1d, gc_acc_ins_du_1d,                                          &
  gc_nd_cor_ins_1d, gc_cor_ins_du_1d,                                          &
  rhcrit_1d, q_1d, qcf_1d, cloud_liq_frac_1d, cloud_blk_frac_1d,               &
  cdnc_1d, drydp, wetdp, rhopar, pvol_wat,                                     &
  pvol )

USE glomap_clim_calc_md_mdt_nd_mod,     ONLY:                                  &
    glomap_clim_calc_md_mdt_nd

USE glomap_clim_calc_rh_frac_clear_mod, ONLY:                                  &
    glomap_clim_calc_rh_frac_clear

USE parkind1,                           ONLY:                                  &
    jpim,                                                                      &
    jprb

USE ukca_calc_drydiam_mod,              ONLY:                                  &
    ukca_calc_drydiam

USE ukca_cdnc_jones_mod,                ONLY:                                  &
    ukca_cdnc_jones

USE ukca_config_specification_mod,      ONLY:                                  &
    glomap_variables_climatology

USE ukca_mode_setup,                    ONLY:                                  &
    cp_su, cp_bc, cp_oc, cp_cl, cp_du,                                         &
    nmodes

USE ukca_volume_mode_mod,               ONLY:                                  &
    ukca_volume_mode

USE yomhook,                            ONLY:                                  &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: n_points
INTEGER, INTENT(IN) :: i_glomap_clim_setup

LOGICAL, INTENT(IN) :: rad_this_tstep
LOGICAL, INTENT(IN) :: l_radaer

REAL(KIND=real_umphys), INTENT(IN) :: act_radius

! pressure on theta levels
REAL(KIND=real_umphys), INTENT(IN) :: p_theta_levels_1d( n_points )

! temperature on theta levels
REAL(KIND=real_umphys), INTENT(IN) :: t_theta_levels_1d( n_points )

! aerosol fields
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_nuc_sol_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nuc_sol_su_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nuc_sol_om_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_ait_sol_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_ait_sol_su_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_ait_sol_bc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_ait_sol_om_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_acc_sol_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_acc_sol_su_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_acc_sol_bc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_acc_sol_om_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_acc_sol_ss_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_acc_sol_du_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_cor_sol_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_cor_sol_su_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_cor_sol_bc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_cor_sol_om_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_cor_sol_ss_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_cor_sol_du_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_ait_ins_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_ait_ins_bc_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_ait_ins_om_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_acc_ins_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_acc_ins_du_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_nd_cor_ins_1d( n_points )
REAL(KIND=real_umphys), INTENT(IN) ::  gc_cor_ins_du_1d( n_points )

! Critical relative humidity
REAL(KIND=real_umphys), INTENT(IN) ::         rhcrit_1d( n_points )

! relative humidity on theta levels
REAL(KIND=real_umphys), INTENT(IN) ::              q_1d( n_points )

! qcf on theta levels
REAL(KIND=real_umphys), INTENT(IN) ::            qcf_1d( n_points )

! Cloud liquid fraction
REAL(KIND=real_umphys), INTENT(IN) :: cloud_liq_frac_1d( n_points )

! Cloud bulk fraction
REAL(KIND=real_umphys), INTENT(IN) :: cloud_blk_frac_1d( n_points )

! Cloud droplet number concentration from Jones method (m^-3^)
REAL(KIND=real_umphys), INTENT(OUT) ::          cdnc_1d( n_points )

! Median particle dry diameter for each mode (m)
REAL(KIND=real_umphys), INTENT(OUT) ::    drydp( n_points, nmodes )

! Avg wet diameter of size mode (m)
REAL(KIND=real_umphys), INTENT(OUT) ::    wetdp( n_points, nmodes )

! Particle density (kg/m^3) [includes H2O & insoluble cpts]
REAL(KIND=real_umphys), INTENT(OUT) ::   rhopar( n_points, nmodes )

! Partial volume of water in each mode (m3)
REAL(KIND=real_umphys), INTENT(OUT) :: pvol_wat( n_points, nmodes )

! Partial volumes of each component in each mode (m^3^)
REAL(KIND=real_umphys), INTENT(OUT) ::  pvol( n_points, nmodes,                &
                                              glomap_variables_climatology%ncp )

! Local variables

! Aerosol ptcl no. concentration (ptcls per cc)
REAL(KIND=real_umphys) ::    nd( n_points, nmodes )

! Total median aerosol mass (molecules per ptcl)
REAL(KIND=real_umphys) ::   mdt( n_points, nmodes )

! Component median aerosol mass (molecules per ptcl)
REAL(KIND=real_umphys) :: md(n_points,nmodes,glomap_variables_climatology%ncp)

! Avg dry volume of size mode (cubic metres)
REAL(KIND=real_umphys) ::  dvol( n_points, nmodes )

! Clear sky relative humidity as a fraction 1-D
REAL(KIND=real_umphys) :: rh_clr_1d ( n_points )

! These fields are only required by UKCA (not climatology).
! These fields arguments passed OUT of ukca_volume_mode.
REAL(KIND=real_umphys) ::  wvol( n_points, nmodes )
REAL(KIND=real_umphys) :: mdwat( n_points, nmodes )

! Field ccn_1d not required here but is a required argument elsewhere
REAL(KIND=real_umphys) ::  ccn_1d( n_points )

CHARACTER(LEN=*), PARAMETER   :: RoutineName='GLOMAP_CLIM_INTERFACE'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialise fields
cdnc_1d(:)    = 0.0
drydp(:,:)    = 0.0
wetdp(:,:)    = 0.0
pvol_wat(:,:) = 0.0
nd(:,:)       = 0.0
mdt(:,:)      = 0.0
md(:,:,:)     = 0.0
dvol(:,:)     = 0.0
rh_clr_1d(:)  = 0.0
wvol(:,:)     = 0.0
mdwat(:,:)    = 0.0
ccn_1d(:)     = 0.0
pvol(:,:,:)   = 0.0

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
CALL glomap_clim_calc_md_mdt_nd( n_points, glomap_variables_climatology%ncp,   &
  i_glomap_clim_setup,                                                         &
  p_theta_levels_1d, t_theta_levels_1d,                                        &
  gc_nd_nuc_sol_1d, gc_nuc_sol_su_1d, gc_nuc_sol_om_1d,                        &
  gc_nd_ait_sol_1d, gc_ait_sol_su_1d, gc_ait_sol_bc_1d, gc_ait_sol_om_1d,      &
  gc_nd_acc_sol_1d, gc_acc_sol_su_1d, gc_acc_sol_bc_1d, gc_acc_sol_om_1d,      &
  gc_acc_sol_ss_1d, gc_acc_sol_du_1d,                                          &
  gc_nd_cor_sol_1d, gc_cor_sol_su_1d, gc_cor_sol_bc_1d, gc_cor_sol_om_1d,      &
  gc_cor_sol_ss_1d, gc_cor_sol_du_1d,                                          &
  gc_nd_ait_ins_1d, gc_ait_ins_bc_1d, gc_ait_ins_om_1d,                        &
  gc_nd_acc_ins_1d, gc_acc_ins_du_1d,                                          &
  gc_nd_cor_ins_1d, gc_cor_ins_du_1d,                                          &
  nd, md, mdt )

! Calculate the dry diameters and volumes
CALL ukca_calc_drydiam( n_points, glomap_variables_climatology,                &
                        nd, md, mdt, drydp, dvol )

! The following aerosol properties are required by radaer, and thus only need
! calculating if both radaer is being used, and it is a radiation timestep:
!     wetdp , rhopar , pvol , pvol_wat
IF ( rad_this_tstep .AND. l_radaer ) THEN
  ! Calculate wet diameters, densities, partial volumes...
  CALL ukca_volume_mode( glomap_variables_climatology, n_points, nd, md, mdt,  &
                         rh_clr_1d, dvol, drydp,                               &
                         t_theta_levels_1d, p_theta_levels_1d, q_1d,           &
                         mdwat, wvol, wetdp, rhopar, pvol, pvol_wat )
END IF

! Obtain CDNC using Jones method doi:10.1038/370450a0
CALL ukca_cdnc_jones( n_points, act_radius, drydp, nd,                         &
                      glomap_variables_climatology, ccn_1d, cdnc_1d )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE glomap_clim_interface

END MODULE glomap_clim_interface_mod
