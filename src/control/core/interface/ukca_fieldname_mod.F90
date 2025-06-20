! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Field name definitions for UKCA.
!
! Method:
!
!  This module provides field name parameters for referring to fields
!  within UKCA. These parameters should be used in place of character literals
!  so that any changes made to names in this module will propagate throughout
!  the code base. The literals defined here will be used in communications with
!  parent applications via the UKCA API. They are thus considered part of the
!  API definition and cannot be changed without affecting the API.
!
!  N.B. IMPLEMENTATION OF THIS METHOD IS NOT YET COMPLETE: NAMES SHOULD NOT
!  BE MODIFIED (EVEN AS PART OF AN API CHANGE) UNTIL THE UKCA CODE BASE IS
!  FULLY COMPLIANT
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

MODULE ukca_fieldname_mod

IMPLICIT NONE

PUBLIC

INTEGER, PARAMETER :: maxlen_fieldname = 20
INTEGER, PARAMETER :: maxlen_diagname = 160

! -----------------------------------------------------------------
! Prognostic Fields (i.e. Tracers and Non-transported Prognostics)
! -----------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: fldname_nuc_sol_n = 'Nuc_SOL_N'
CHARACTER(LEN=*), PARAMETER :: fldname_nuc_sol_su = 'Nuc_SOL_SU'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_sol_n = 'Ait_SOL_N'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_sol_su = 'Ait_SOL_SU'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_sol_bc = 'Ait_SOL_BC'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_sol_om = 'Ait_SOL_OM'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_sol_n = 'Acc_SOL_N'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_sol_su = 'Acc_SOL_SU'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_sol_bc = 'Acc_SOL_BC'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_sol_om = 'Acc_SOL_OM'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_sol_ss = 'Acc_SOL_SS'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_sol_du = 'Acc_SOL_DU'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_sol_n = 'Cor_SOL_N'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_sol_su = 'Cor_SOL_SU'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_sol_bc = 'Cor_SOL_BC'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_sol_om = 'Cor_SOL_OM'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_sol_ss = 'Cor_SOL_SS'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_sol_du = 'Cor_SOL_DU'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_ins_n = 'Ait_INS_N'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_ins_bc = 'Ait_INS_BC'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_ins_om = 'Ait_INS_OM'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_ins_n = 'Acc_INS_N'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_ins_du = 'Acc_INS_DU'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_ins_n = 'Cor_INS_N'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_ins_du = 'Cor_INS_DU'
CHARACTER(LEN=*), PARAMETER :: fldname_nuc_sol_om = 'Nuc_SOL_OM'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_sol_ss = 'Ait_SOL_SS'
CHARACTER(LEN=*), PARAMETER :: fldname_nuc_sol_so = 'Nuc_SOL_SO'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_sol_so = 'Ait_SOL_SO'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_sol_so = 'Acc_SOL_SO'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_sol_so = 'Cor_SOL_SO'
CHARACTER(LEN=*), PARAMETER :: fldname_nuc_sol_nh = 'Nuc_SOL_NH'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_sol_nh = 'Ait_SOL_NH'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_sol_nh = 'Acc_SOL_NH'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_sol_nh = 'Cor_SOL_NH'
CHARACTER(LEN=*), PARAMETER :: fldname_nuc_sol_nt = 'Nuc_SOL_NT'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_sol_nt = 'Ait_SOL_NT'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_sol_nt = 'Acc_SOL_NT'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_sol_nt = 'Cor_SOL_NT'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_sol_nn = 'Acc_SOL_NN'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_sol_nn = 'Cor_SOL_NN'
CHARACTER(LEN=*), PARAMETER :: fldname_sup_ins_n = 'Sup_INS_N'
CHARACTER(LEN=*), PARAMETER :: fldname_sup_ins_du = 'Sup_INS_DU'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_sol_mp = 'Ait_SOL_MP'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_sol_mp = 'Acc_SOL_MP'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_sol_mp = 'Cor_SOL_MP'
CHARACTER(LEN=*), PARAMETER :: fldname_ait_ins_mp = 'Ait_INS_MP'
CHARACTER(LEN=*), PARAMETER :: fldname_acc_ins_mp = 'Acc_INS_MP'
CHARACTER(LEN=*), PARAMETER :: fldname_cor_ins_mp = 'Cor_INS_MP'
CHARACTER(LEN=*), PARAMETER :: fldname_sup_ins_mp = 'Sup_INS_MP'
CHARACTER(LEN=*), PARAMETER :: fldname_passive_o3 = 'PASSIVE O3'
CHARACTER(LEN=*), PARAMETER :: fldname_age_of_air = 'AGE OF AIR'
CHARACTER(LEN=*), PARAMETER :: fldname_RTX22O2='RTX22O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RTX24O2='RTX24O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_NRTX28O2='NRTX28O2  '
CHARACTER(LEN=*), PARAMETER :: fldname_RTX28O2='RTX28O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RTN10O2='RTN10O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RTN14O2='RTN14O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RTN23O2='RTN23O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RTN24O2='RTN24O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RTN25O2='RTN25O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RTN26O2='RTN26O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_NRTN28O2='NRTN28O2  '
CHARACTER(LEN=*), PARAMETER :: fldname_RTN28O2='RTN28O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_NRU12O2='NRU12O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_NRU14O2='NRU14O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_NRN12O2='NRN12O2   '
CHARACTER(LEN=*), PARAMETER :: fldname_NRN9O2='NRN9O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_NRN6O2='NRN6O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RU10O2='RU10O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RU12O2='RU12O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RU14O2='RU14O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RN17O2='RN17O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RN14O2='RN14O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RN11O2='RN11O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RN8O2='RN8O2     '
CHARACTER(LEN=*), PARAMETER :: fldname_HOCH2CO3='HOCH2CO3  '
CHARACTER(LEN=*), PARAMETER :: fldname_RN18AO2='RN18AO2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RN15AO2='RN15AO2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RN18O2='RN18O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RN15O2='RN15O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RN12O2='RN12O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RN9O2='RN9O2     '
CHARACTER(LEN=*), PARAMETER :: fldname_HOCH2CH2O2='HOCH2CH2O2'
CHARACTER(LEN=*), PARAMETER :: fldname_RA19CO2='RA19CO2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RA19AO2='RA19AO2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RA16O2='RA16O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RA13O2='RA13O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RN16AO2='RN16AO2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RN13AO2='RN13AO2   '
CHARACTER(LEN=*), PARAMETER :: fldname_RN19O2='RN19O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RN16O2='RN16O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RN13O2='RN13O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_RN10O2='RN10O2    '
CHARACTER(LEN=*), PARAMETER :: fldname_DHPR12O2='DHPR12O2  '
CHARACTER(LEN=*), PARAMETER :: fldname_RU10AO2='RU10AO2   '
CHARACTER(LEN=*), PARAMETER :: fldname_MACO3='MACO3     '
CHARACTER(LEN=*), PARAMETER :: fldname_n_activ_sum = 'n_activ_sum'
CHARACTER(LEN=*), PARAMETER :: fldname_MACRO2 = 'MACRO2'
CHARACTER(LEN=*), PARAMETER :: fldname_ISO2 = 'ISO2'
CHARACTER(LEN=*), PARAMETER :: fldname_surfarea = 'surfarea'
CHARACTER(LEN=*), PARAMETER :: fldname_cdnc = 'cdnc'
CHARACTER(LEN=*), PARAMETER :: fldname_HO2S = 'HO2S'
CHARACTER(LEN=*), PARAMETER :: fldname_OHS = 'OHS'
CHARACTER(LEN=*), PARAMETER :: fldname_O1DS = 'O(1D)S'
CHARACTER(LEN=*), PARAMETER :: fldname_O3PS = 'O(3P)S'
CHARACTER(LEN=*), PARAMETER :: fldname_het_ho2 = 'het_ho2'
CHARACTER(LEN=*), PARAMETER :: fldname_het_n2o5 = 'het_n2o5'
CHARACTER(LEN=*), PARAMETER :: fldname_TOLP1 = 'TOLP1'
CHARACTER(LEN=*), PARAMETER :: fldname_HOIPO2 = 'HOIPO2'
CHARACTER(LEN=*), PARAMETER :: fldname_HOMVKO2 = 'HOMVKO2'
CHARACTER(LEN=*), PARAMETER :: fldname_MEMALD1 = 'MEMALD1'
CHARACTER(LEN=*), PARAMETER :: fldname_OXYL1 = 'OXYL1'
CHARACTER(LEN=*), PARAMETER :: fldname_HOC3H6O2 = 'HOC3H6O2'
CHARACTER(LEN=*), PARAMETER :: fldname_HOC2H4O2 = 'HOC2H4O2'
CHARACTER(LEN=*), PARAMETER :: fldname_MEKO2 = 'MEKO2'
CHARACTER(LEN=*), PARAMETER :: fldname_MeCOCH2OO = 'MeCOCH2OO'
CHARACTER(LEN=*), PARAMETER :: fldname_MeCOC2OO = 'MeCOC2OO'
CHARACTER(LEN=*), PARAMETER :: fldname_EtCO3 = 'EtCO3'
CHARACTER(LEN=*), PARAMETER :: fldname_i_PrOO = 'i-PrOO'
CHARACTER(LEN=*), PARAMETER :: fldname_s_BuOO = 's-BuOO'
CHARACTER(LEN=*), PARAMETER :: fldname_n_PrOO = 'n-PrOO'
CHARACTER(LEN=*), PARAMETER :: fldname_MeCO3 = 'MeCO3'
CHARACTER(LEN=*), PARAMETER :: fldname_EtOO = 'EtOO'
CHARACTER(LEN=*), PARAMETER :: fldname_MeOO = 'MeOO'
CHARACTER(LEN=*), PARAMETER :: fldname_HCl = 'HCl'
CHARACTER(LEN=*), PARAMETER :: fldname_HO2 = 'HO2'
CHARACTER(LEN=*), PARAMETER :: fldname_BrO = 'BrO'
CHARACTER(LEN=*), PARAMETER :: fldname_OH = 'OH'
CHARACTER(LEN=*), PARAMETER :: fldname_NO2 = 'NO2'
CHARACTER(LEN=*), PARAMETER :: fldname_O1D = 'O(1D)'
CHARACTER(LEN=*), PARAMETER :: fldname_O3P = 'O(3P)'
CHARACTER(LEN=*), PARAMETER :: fldname_drydiam_nuc_sol = 'drydiam_nuc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_drydiam_ait_sol = 'drydiam_ait_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_drydiam_acc_sol = 'drydiam_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_drydiam_cor_sol = 'drydiam_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_drydiam_ait_insol = 'drydiam_ait_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_drydiam_acc_insol = 'drydiam_acc_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_drydiam_cor_insol = 'drydiam_cor_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_drydiam_sup_insol = 'drydiam_sup_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_wetdiam_ait_sol = 'wetdiam_ait_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_wetdiam_acc_sol = 'wetdiam_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_wetdiam_cor_sol = 'wetdiam_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_aerdens_ait_sol = 'aerdens_ait_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_aerdens_acc_sol = 'aerdens_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_aerdens_cor_sol = 'aerdens_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_aerdens_ait_insol = 'aerdens_ait_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_aerdens_acc_insol = 'aerdens_acc_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_aerdens_cor_insol = 'aerdens_cor_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_aerdens_sup_insol = 'aerdens_sup_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_su_ait_sol = 'pvol_su_ait_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_bc_ait_sol = 'pvol_bc_ait_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_oc_ait_sol = 'pvol_oc_ait_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_so_ait_sol = 'pvol_so_ait_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_no3_ait_sol = 'pvol_no3_ait_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_nh4_ait_sol = 'pvol_nh4_ait_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_mp_ait_sol = 'pvol_mp_ait_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_h2o_ait_sol = 'pvol_h2o_ait_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_su_acc_sol = 'pvol_su_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_bc_acc_sol = 'pvol_bc_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_oc_acc_sol = 'pvol_oc_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_ss_acc_sol = 'pvol_ss_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_no3_acc_sol = 'pvol_no3_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_nh4_acc_sol = 'pvol_nh4_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_nn_acc_sol = 'pvol_nn_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_du_acc_sol = 'pvol_du_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_so_acc_sol = 'pvol_so_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_mp_acc_sol = 'pvol_mp_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_h2o_acc_sol = 'pvol_h2o_acc_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_su_cor_sol = 'pvol_su_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_bc_cor_sol = 'pvol_bc_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_oc_cor_sol = 'pvol_oc_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_ss_cor_sol = 'pvol_ss_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_no3_cor_sol = 'pvol_no3_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_nh4_cor_sol = 'pvol_nh4_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_nn_cor_sol = 'pvol_nn_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_du_cor_sol = 'pvol_du_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_so_cor_sol = 'pvol_so_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_mp_cor_sol = 'pvol_mp_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_h2o_cor_sol = 'pvol_h2o_cor_sol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_bc_ait_insol = 'pvol_bc_ait_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_oc_ait_insol = 'pvol_oc_ait_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_mp_ait_insol = 'pvol_mp_ait_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_du_acc_insol = 'pvol_du_acc_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_mp_acc_insol = 'pvol_mp_acc_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_du_cor_insol = 'pvol_du_cor_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_mp_cor_insol = 'pvol_mp_cor_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_du_sup_insol = 'pvol_du_sup_insol'
CHARACTER(LEN=*), PARAMETER :: fldname_pvol_mp_sup_insol = 'pvol_mp_sup_insol'

!--------------------
! Environment fields
!--------------------

! --- Scalar values of type real ---
CHARACTER(LEN=*), PARAMETER :: fldname_sin_declination =                       &
                                      'sin_declination'
CHARACTER(LEN=*), PARAMETER :: fldname_equation_of_time =                      &
                                      'equation_of_time'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_ch4 =                       &
                                      'atmospheric_ch4'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_co2 =                       &
                                      'atmospheric_co2'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_h2 =                        &
                                      'atmospheric_h2'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_n2 =                        &
                                      'atmospheric_n2'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_o2 =                        &
                                      'atmospheric_o2'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_n2o =                       &
                                      'atmospheric_n2o'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_cfc11 =                     &
                                      'atmospheric_cfc11'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_cfc12 =                     &
                                      'atmospheric_cfc12'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_cfc113 =                    &
                                      'atmospheric_cfc113'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_hcfc22 =                    &
                                      'atmospheric_hcfc22'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_hfc125 =                    &
                                      'atmospheric_hfc125'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_hfc134a =                   &
                                      'atmospheric_hfc134a'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_mebr =                      &
                                      'atmospheric_mebr'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_mecl =                      &
                                      'atmospheric_mecl'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_ch2br2 =                    &
                                      'atmospheric_ch2br2'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_cfc114 =                    &
                                      'atmospheric_cfc114'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_cfc115 =                    &
                                      'atmospheric_cfc115'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_ccl4 =                      &
                                      'atmospheric_ccl4'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_meccl3 =                    &
                                      'atmospheric_meccl3'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_hcfc141b =                  &
                                      'atmospheric_hcfc141b'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_hcfc142b =                  &
                                      'atmospheric_hcfc142b'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_h1211 =                     &
                                      'atmospheric_h1211'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_h1202 =                     &
                                      'atmospheric_h1202'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_h1301 =                     &
                                      'atmospheric_h1301'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_h2402 =                     &
                                      'atmospheric_h2402'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_cos =                       &
                                      'atmospheric_cos'
CHARACTER(LEN=*), PARAMETER :: fldname_atmospheric_chbr3 =                     &
                                      'atmospheric_chbr3'
! --- 1D fields of type real ---
CHARACTER(LEN=*), PARAMETER :: fldname_soil_moisture_layer1 =                  &
                                      'soil_moisture_layer1'
CHARACTER(LEN=*), PARAMETER :: fldname_fland =                                 &
                                      'fland'

! --- 2D fields of type real ---
CHARACTER(LEN=*), PARAMETER :: fldname_latitude =                              &
                                      'latitude'
CHARACTER(LEN=*), PARAMETER :: fldname_longitude =                             &
                                      'longitude'
CHARACTER(LEN=*), PARAMETER :: fldname_sin_latitude =                          &
                                      'sin_latitude'
CHARACTER(LEN=*), PARAMETER :: fldname_cos_latitude =                          &
                                      'cos_latitude'
CHARACTER(LEN=*), PARAMETER :: fldname_tan_latitude =                          &
                                      'tan_latitude'
CHARACTER(LEN=*), PARAMETER :: fldname_conv_cloud_lwp =                        &
                                      'conv_cloud_lwp'
CHARACTER(LEN=*), PARAMETER :: fldname_tstar =                                 &
                                      'tstar'
CHARACTER(LEN=*), PARAMETER :: fldname_zbl =                                   &
                                      'zbl'
CHARACTER(LEN=*), PARAMETER :: fldname_rough_length =                          &
                                      'rough_length'
CHARACTER(LEN=*), PARAMETER :: fldname_seaice_frac =                           &
                                      'seaice_frac'
CHARACTER(LEN=*), PARAMETER :: fldname_frac_types =                            &
                                      'frac_types'
CHARACTER(LEN=*), PARAMETER :: fldname_laift_lp =                              &
                                      'laift_lp'
CHARACTER(LEN=*), PARAMETER :: fldname_canhtft_lp =                            &
                                      'canhtft_lp'
CHARACTER(LEN=*), PARAMETER :: fldname_tstar_tile =                            &
                                      'tstar_tile'
CHARACTER(LEN=*), PARAMETER :: fldname_z0tile_lp =                             &
                                      'z0tile_lp'
CHARACTER(LEN=*), PARAMETER :: fldname_pstar =                                 &
                                      'pstar'
CHARACTER(LEN=*), PARAMETER :: fldname_surf_albedo =                           &
                                      'surf_albedo'
CHARACTER(LEN=*), PARAMETER :: fldname_zhsc =                                  &
                                      'zhsc'
CHARACTER(LEN=*), PARAMETER :: fldname_u_scalar_10m =                          &
                                      'u_scalar_10m'
CHARACTER(LEN=*), PARAMETER :: fldname_surf_hf =                               &
                                      'surf_hf'
CHARACTER(LEN=*), PARAMETER :: fldname_u_s =                                   &
                                      'u_s'
CHARACTER(LEN=*), PARAMETER :: fldname_ch4_wetl_emiss =                        &
                                      'ch4_wetl_emiss'
CHARACTER(LEN=*), PARAMETER :: fldname_dms_sea_conc =                          &
                                      'dms_sea_conc'
CHARACTER(LEN=*), PARAMETER :: fldname_chloro_sea =                            &
                                      'chloro_sea'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_flux_div1 =                        &
                                      'dust_flux_div1'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_flux_div2 =                        &
                                      'dust_flux_div2'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_flux_div3 =                        &
                                      'dust_flux_div3'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_flux_div4 =                        &
                                      'dust_flux_div4'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_flux_div5 =                        &
                                      'dust_flux_div5'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_flux_div6 =                        &
                                      'dust_flux_div6'
CHARACTER(LEN=*), PARAMETER :: fldname_surf_wetness =                          &
                                      'surf_wetness'
CHARACTER(LEN=*), PARAMETER :: fldname_grid_surf_area =                        &
                                      'grid_surf_area'
CHARACTER(LEN=*), PARAMETER :: fldname_ext_cg_flash  =                         &
                                      'ext_cg_flash'
CHARACTER(LEN=*), PARAMETER :: fldname_ext_ic_flash  =                         &
                                      'ext_ic_flash'

! Interactive biogenic emissions
CHARACTER(LEN=*), PARAMETER :: fldname_ibvoc_isoprene = 'ibvoc_isoprene'
CHARACTER(LEN=*), PARAMETER :: fldname_ibvoc_terpene  = 'ibvoc_terpene'
CHARACTER(LEN=*), PARAMETER :: fldname_ibvoc_methanol = 'ibvoc_methanol'
CHARACTER(LEN=*), PARAMETER :: fldname_ibvoc_acetone  = 'ibvoc_acetone'

! Interactive fire emissions
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_bc    = 'inferno_bc'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_ch4   = 'inferno_ch4'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_co    = 'inferno_co'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_nox   = 'inferno_nox'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_oc    = 'inferno_oc'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_so2   = 'inferno_so2'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_c2h4  = 'inferno_c2h4'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_c2h6  = 'inferno_c2h6'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_c3h8  = 'inferno_c3h8'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_hcho  = 'inferno_hcho'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_mecho = 'inferno_mecho'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_nh3   = 'inferno_nh3'
CHARACTER(LEN=*), PARAMETER :: fldname_inferno_dms   = 'inferno_dms'


! --- 2D fields of type integer ---
CHARACTER(LEN=*), PARAMETER :: fldname_kent =                                  &
                                      'kent'
CHARACTER(LEN=*), PARAMETER :: fldname_kent_dsc =                              &
                                      'kent_dsc'
CHARACTER(LEN=*), PARAMETER :: fldname_conv_cloud_base =                       &
                                      'conv_cloud_base'
CHARACTER(LEN=*), PARAMETER :: fldname_conv_cloud_top =                        &
                                      'conv_cloud_top'
CHARACTER(LEN=*), PARAMETER :: fldname_lscat_zhang    =                        &
                                      'lscat_zhang'

! --- 2D fields of type logical ---
CHARACTER(LEN=*), PARAMETER :: fldname_land_sea_mask =                         &
                                      'land_sea_mask'
CHARACTER(LEN=*), PARAMETER :: fldname_l_tile_active =                         &
                                      'l_tile_active'

! --- 3D fields of type real ---
! Note: field names for species used in chemistry schemes are capitalised for
! consistency with the corresponding tracer names. They are currently
! defined this way in module 'ukca_chem_master_mod'. This is seen as a
! temporary constraint pending replacement of the field names throughout UKCA
! by parameters defined in this module, at which point the names in this module
! may be considered for review (taking into account that any changes are likely
! to break backwards compatibility with parent models and require
! corresponding changes thereof).

CHARACTER(LEN=*), PARAMETER :: fldname_u_rho_levels =                          &
                                      'u_rho_levels'
CHARACTER(LEN=*), PARAMETER :: fldname_v_rho_levels =                          &
                                      'v_rho_levels'
CHARACTER(LEN=*), PARAMETER :: fldname_geopH_on_theta_mlevs =                  &
                                      'geopH_on_theta_mlevs'
CHARACTER(LEN=*), PARAMETER :: fldname_theta =                                 &
                                      'theta'
CHARACTER(LEN=*), PARAMETER :: fldname_q =                                     &
                                      'q'
CHARACTER(LEN=*), PARAMETER :: fldname_qcf =                                   &
                                      'qcf'
CHARACTER(LEN=*), PARAMETER :: fldname_conv_cloud_amount =                     &
                                      'conv_cloud_amount'
CHARACTER(LEN=*), PARAMETER :: fldname_rho_r2 =                                &
                                      'rho_r2'
CHARACTER(LEN=*), PARAMETER :: fldname_qcl =                                   &
                                      'qcl'
CHARACTER(LEN=*), PARAMETER :: fldname_exner_rho_levels =                      &
                                      'exner_rho_levels'
CHARACTER(LEN=*), PARAMETER :: fldname_area_cloud_fraction =                   &
                                      'area_cloud_fraction'
CHARACTER(LEN=*), PARAMETER :: fldname_cloud_frac =                            &
                                      'cloud_frac'
CHARACTER(LEN=*), PARAMETER :: fldname_cloud_liq_frac =                        &
                                      'cloud_liq_frac'
CHARACTER(LEN=*), PARAMETER :: fldname_exner_theta_levels =                    &
                                      'exner_theta_levels'
CHARACTER(LEN=*), PARAMETER :: fldname_p_rho_levels =                          &
                                      'p_rho_levels'
CHARACTER(LEN=*), PARAMETER :: fldname_p_theta_levels =                        &
                                      'p_theta_levels'
CHARACTER(LEN=*), PARAMETER :: fldname_rhokh_rdz =                             &
                                      'rhokh_rdz'
CHARACTER(LEN=*), PARAMETER :: fldname_dtrdz =                                 &
                                      'dtrdz'
CHARACTER(LEN=*), PARAMETER :: fldname_we_lim =                                &
                                      'we_lim'
CHARACTER(LEN=*), PARAMETER :: fldname_t_frac =                                &
                                      't_frac'
CHARACTER(LEN=*), PARAMETER :: fldname_zrzi =                                  &
                                      'zrzi'
CHARACTER(LEN=*), PARAMETER :: fldname_we_lim_dsc =                            &
                                      'we_lim_dsc'
CHARACTER(LEN=*), PARAMETER :: fldname_t_frac_dsc =                            &
                                      't_frac_dsc'
CHARACTER(LEN=*), PARAMETER :: fldname_zrzi_dsc =                              &
                                      'zrzi_dsc'
CHARACTER(LEN=*), PARAMETER :: fldname_stcon =                                 &
                                      'stcon'
CHARACTER(LEN=*), PARAMETER :: fldname_ls_rain3d =                             &
                                      'ls_rain3d'
CHARACTER(LEN=*), PARAMETER :: fldname_ls_snow3d =                             &
                                      'ls_snow3d'
CHARACTER(LEN=*), PARAMETER :: fldname_autoconv =                              &
                                      'autoconv'
CHARACTER(LEN=*), PARAMETER :: fldname_accretion =                             &
                                      'accretion'
CHARACTER(LEN=*), PARAMETER :: fldname_pv_on_theta_mlevs =                     &
                                      'pv_on_theta_mlevs'
CHARACTER(LEN=*), PARAMETER :: fldname_conv_rain3d =                           &
                                      'conv_rain3d'
CHARACTER(LEN=*), PARAMETER :: fldname_conv_snow3d =                           &
                                      'conv_snow3d'
CHARACTER(LEN=*), PARAMETER :: fldname_so4_sa_clim =                           &
                                      'so4_sa_clim'
CHARACTER(LEN=*), PARAMETER :: fldname_so4_aitken =                            &
                                      'so4_aitken'
CHARACTER(LEN=*), PARAMETER :: fldname_so4_accum =                             &
                                      'so4_accum'
CHARACTER(LEN=*), PARAMETER :: fldname_soot_fresh =                            &
                                      'soot_fresh'
CHARACTER(LEN=*), PARAMETER :: fldname_soot_aged =                             &
                                      'soot_aged'
CHARACTER(LEN=*), PARAMETER :: fldname_ocff_fresh =                            &
                                      'ocff_fresh'
CHARACTER(LEN=*), PARAMETER :: fldname_ocff_aged =                             &
                                      'ocff_aged'
CHARACTER(LEN=*), PARAMETER :: fldname_biogenic =                              &
                                      'biogenic'
CHARACTER(LEN=*), PARAMETER :: fldname_sea_salt_film =                         &
                                      'sea_salt_film'
CHARACTER(LEN=*), PARAMETER :: fldname_sea_salt_jet =                          &
                                      'sea_salt_jet'
CHARACTER(LEN=*), PARAMETER :: fldname_co2_interactive =                       &
                                      'CO2'
CHARACTER(LEN=*), PARAMETER :: fldname_rim_cry =                               &
                                      'rim_cry'
CHARACTER(LEN=*), PARAMETER :: fldname_rim_agg =                               &
                                      'rim_agg'
CHARACTER(LEN=*), PARAMETER :: fldname_vertvel =                               &
                                      'vertvel'
CHARACTER(LEN=*), PARAMETER :: fldname_bl_tke =                                &
                                      'bl_tke'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_div1 =                             &
                                      'dust_div1'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_div2 =                             &
                                      'dust_div2'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_div3 =                             &
                                      'dust_div3'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_div4 =                             &
                                      'dust_div4'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_div5 =                             &
                                      'dust_div5'
CHARACTER(LEN=*), PARAMETER :: fldname_dust_div6 =                             &
                                      'dust_div6'
CHARACTER(LEN=*), PARAMETER :: fldname_interf_z =                              &
                                      'interf_z'
CHARACTER(LEN=*), PARAMETER :: fldname_grid_area_fullht =                      &
                                      'grid_area_fullht'
CHARACTER(LEN=*), PARAMETER :: fldname_grid_volume =                           &
                                      'grid_volume'
CHARACTER(LEN=*), PARAMETER :: fldname_grid_airmass =                          &
                                      'grid_airmass'
CHARACTER(LEN=*), PARAMETER :: fldname_rel_humid_frac =                        &
                                      'rel_humid_frac'
CHARACTER(LEN=*), PARAMETER :: fldname_rel_humid_frac_clr =                    &
                                      'rel_humid_frac_clr'
CHARACTER(LEN=*), PARAMETER :: fldname_qsvp =                                  &
                                      'qsvp'

! Oxidants for Offline Oxidants chemistry - 3-D of type Real
! (can be used for other chemistry schemes if needed)
CHARACTER(LEN=*), PARAMETER :: fldname_h2o2_offline = 'H2O2'
CHARACTER(LEN=*), PARAMETER :: fldname_ho2_offline  = 'HO2'
CHARACTER(LEN=*), PARAMETER :: fldname_no3_offline  = 'NO3'
CHARACTER(LEN=*), PARAMETER :: fldname_o3_offline   = 'O3'
CHARACTER(LEN=*), PARAMETER :: fldname_oh_offline   = 'OH'

! --- 4D fields of type real ---
CHARACTER(LEN=*), PARAMETER :: fldname_photol_rates = 'photol_rates'

!------------------------------------------------------------------
! Diagnostic fields
! Note 1: CF names or names constructed using CF guidelines are used.
!         The notation [CF] in the comments after a name defined below
!         indicates that the name is listed in the CF Standard Names
!         Table at cfconventions.org.
! Note 2: Parameter names for the diagnostic fields are included in
!         the API definition so neither the parameter names or the
!         corresponding literals can be changed without affecting
!         parent applications. A parent may use either in
!         communications with UKCA via the API.
!------------------------------------------------------------------

! Photolysis rates (s-1) 'jrate_'<phot rate label>. Item numbers below refer to
! entries in 'ratj_defs_master'; see module ukca_chem_master_mod.
CHARACTER(LEN=*), PARAMETER :: diagname_jrate_no2 =                            &
  'photolysis_rate_of_nitrogen_dioxide'
                                 ! [CF] Item 11: NO2 -> NO + O(3P)
CHARACTER(LEN=*), PARAMETER :: diagname_jrate_o3a =                            &
  'photolysis_rate_of_ozone_to_1d_oxygen_atom'
                                 ! [CF] Item 15: O3 -> O2 + O(1D)
CHARACTER(LEN=*), PARAMETER :: diagname_jrate_o3b =                            &
  'photolysis_rate_of_ozone_to_3p_oxygen_atom'
                                 ! Item 16: O3 -> O2 + O(3P)
CHARACTER(LEN=*), PARAMETER :: diagname_jrate_o2b =                            &
  'photolysis_rate_of_molecular_oxygen_to_1d_oxygen_atom'
                                 ! Item 40: O2 -> O(3P) + O(1D)

! ASAD framework diagnostics by STASH code (only 1 currently)
! A UM STASH code is used as an identifier for each ASAD diagnostic. Those
! listed below can be requested by name via the UKCA API by non-UM parent
! models and by the UM
CHARACTER(LEN=*), PARAMETER :: diagname_rxnflux_oh_ch4_trop =                  &
  'minus_tendency_of_troposphere_moles_of_methane_due_to_' //                  &
  'reaction_with_hydroxyl_radical'
                                 ! 50041: Reaction rate OH + CH4 -> H2O + MeOO
                                 ! in the troposphere (mol s-1)

! Pressure of PV-theta tropopause surface (Pa)
CHARACTER(LEN=*), PARAMETER :: diagname_p_tropopause =                         &
  'tropopause_air_pressure'      ! [CF]

! [CF] Ozone column in Dobson units
CHARACTER(LEN=*), PARAMETER :: diagname_o3_column_du =                         &
  'equivalent_thickness_at_stp_of_atmosphere_ozone_content'  ! [CF]

! Plume height from exp eruptions (m)
CHARACTER(LEN=*), PARAMETER :: diagname_plumeria_height =                      &
  'plumeria_height_of_explosive_eruptions'

CONTAINS

! ----------------------------------------------------------------------
LOGICAL FUNCTION is_mode_ntp(varname)
! ----------------------------------------------------------------------
! Description:
! Determine whether field name is one for a GLOMAP-mode non-transported
! prognostic field.
! ----------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: varname

IF (LEN(varname) < 8) THEN
  is_mode_ntp = .FALSE.
ELSE
  IF (varname == 'surfarea' .OR.                                               &
      varname (1:7) == 'drydiam' .OR.                                          &
      varname (1:7) == 'wetdiam' .OR.                                          &
      varname (1:7) == 'aerdens' .OR.                                          &
      varname (1:4) == 'pvol' ) THEN
    is_mode_ntp = .TRUE.
  ELSE
    is_mode_ntp = .FALSE.
  END IF
END IF

RETURN

END FUNCTION is_mode_ntp

END MODULE ukca_fieldname_mod
