! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!   Calculate md, mdt, nd fields
!   mmr1d, nmr1d, aird are local fields to save memory
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE glomap_clim_calc_md_mdt_nd_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE
PRIVATE

PUBLIC :: glomap_clim_calc_md_mdt_nd

CHARACTER(LEN=*),PARAMETER,PRIVATE:: ModuleName='GLOMAP_CLIM_CALC_MD_MDT_ND_MOD'

CONTAINS

SUBROUTINE glomap_clim_calc_md_mdt_nd ( n_points, i_glomap_clim_setup_in,      &
  p_theta_levels_1d, t_theta_levels_1d,                                        &
  gc_nd_nuc_sol_1d, gc_nuc_sol_su_1d, gc_nuc_sol_oc_1d,                        &
  gc_nd_ait_sol_1d, gc_ait_sol_su_1d, gc_ait_sol_bc_1d, gc_ait_sol_oc_1d,      &
  gc_nd_acc_sol_1d, gc_acc_sol_su_1d, gc_acc_sol_bc_1d, gc_acc_sol_oc_1d,      &
  gc_acc_sol_ss_1d, gc_acc_sol_du_1d,                                          &
  gc_nd_cor_sol_1d, gc_cor_sol_su_1d, gc_cor_sol_bc_1d, gc_cor_sol_oc_1d,      &
  gc_cor_sol_ss_1d, gc_cor_sol_du_1d,                                          &
  gc_nd_ait_ins_1d, gc_ait_ins_bc_1d, gc_ait_ins_oc_1d,                        &
  gc_nd_acc_ins_1d, gc_acc_ins_du_1d,                                          &
  gc_nd_cor_ins_1d, gc_cor_ins_du_1d,                                          &
  nd, md, mdt )

USE get_gc_aerosol_fields_1d_mod,    ONLY:                                     &
    get_gc_aerosol_fields_1d

USE glomap_clim_calc_aird_mod,       ONLY:                                     &
    glomap_clim_calc_aird

USE glomap_clim_pop_md_mdt_nd_mod,   ONLY:                                     &
    glomap_clim_pop_md_mdt_nd

USE parkind1,                        ONLY:                                     &
    jprb,                                                                      &
    jpim

USE ukca_mode_setup,                 ONLY:                                     &
    ncp,                                                                       &
    nmodes

USE yomhook,                         ONLY:                                     &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments


INTEGER, INTENT(IN) :: n_points

INTEGER, INTENT(IN) :: i_glomap_clim_setup_in

REAL(KIND=real_umphys),INTENT(IN) :: p_theta_levels_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: t_theta_levels_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_nuc_sol_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nuc_sol_su_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nuc_sol_oc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_ait_sol_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_sol_su_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_sol_bc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_sol_oc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_acc_sol_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_su_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_bc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_oc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_ss_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_du_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_cor_sol_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_su_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_bc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_oc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_ss_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_du_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_ait_ins_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_ins_bc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_ins_oc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_acc_ins_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_ins_du_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_cor_ins_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_ins_du_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN OUT) :: nd( n_points, nmodes )

REAL(KIND=real_umphys),INTENT(IN OUT) :: md( n_points, nmodes, ncp )

REAL(KIND=real_umphys),INTENT(IN OUT) :: mdt( n_points, nmodes )

! Local variables

! Dry air density
REAL(KIND=real_umphys) ::  aird( n_points )

! Aerosol number mixing ratio ( particles / molecule of air )
REAL(KIND=real_umphys) :: nmr1d( n_points, nmodes )

! Aerosol mass mixing ratio ( kg / kg of air )
REAL(KIND=real_umphys) :: mmr1d( n_points, nmodes, ncp)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName='GLOMAP_CLIM_CALC_MD_MDT_ND'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialise arrays
nd(:,:)   = 0.0
md(:,:,:) = 0.0
mdt(:,:) = 0.0

! Calculate aird
CALL glomap_clim_calc_aird( n_points, p_theta_levels_1d, t_theta_levels_1d,    &
                            aird )

! Copy required fields from climatology aerosol array pointer into nmr1d & mmr1d
CALL get_gc_aerosol_fields_1d( n_points, i_glomap_clim_setup_in,               &
  gc_nd_nuc_sol_1d, gc_nuc_sol_su_1d, gc_nuc_sol_oc_1d,                        &
  gc_nd_ait_sol_1d, gc_ait_sol_su_1d, gc_ait_sol_bc_1d, gc_ait_sol_oc_1d,      &
  gc_nd_acc_sol_1d, gc_acc_sol_su_1d, gc_acc_sol_bc_1d, gc_acc_sol_oc_1d,      &
  gc_acc_sol_ss_1d, gc_acc_sol_du_1d,                                          &
  gc_nd_cor_sol_1d, gc_cor_sol_su_1d, gc_cor_sol_bc_1d, gc_cor_sol_oc_1d,      &
  gc_cor_sol_ss_1d, gc_cor_sol_du_1d,                                          &
  gc_nd_ait_ins_1d, gc_ait_ins_bc_1d, gc_ait_ins_oc_1d,                        &
  gc_nd_acc_ins_1d, gc_acc_ins_du_1d,                                          &
  gc_nd_cor_ins_1d, gc_cor_ins_du_1d,                                          &
  nmr1d, mmr1d )

! Populate the md, mdt & nd arrays
CALL glomap_clim_pop_md_mdt_nd ( i_glomap_clim_setup_in, n_points,             &
                                 aird, mmr1d, nmr1d, md, mdt, nd )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE glomap_clim_calc_md_mdt_nd
END MODULE glomap_clim_calc_md_mdt_nd_mod
