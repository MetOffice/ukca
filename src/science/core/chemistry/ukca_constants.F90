! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to contain constants used in UKCA
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_constants

USE ukca_um_legacy_mod, ONLY: pi, recip_pi_over_180, rhour_per_day
USE chemistry_constants_mod, ONLY: avogadro

IMPLICIT NONE

! UM definitions of MODE constants
REAL, PARAMETER :: recip_avogad=1.0/Avogadro

! MODE/aerosol chemistry constants not already in UM
REAL, PARAMETER :: nmol=1.0e2          !
REAL, PARAMETER :: tdays=0.0           !
REAL, PARAMETER :: ems_eps=1.0e-8      !
REAL, PARAMETER :: conc_eps=1.0e-8     !
REAL, PARAMETER :: dn_eps=1.0e-8       !
REAL, PARAMETER :: bconst=3.44e13      ! Volume_mode
REAL, PARAMETER :: nu_h2so4=3.0        ! Volume_mode
REAL, PARAMETER :: H_plus=1.0e-5       ! cloud/rain H+ concentration

! Dobson Unit (molecules per m2)
REAL, PARAMETER :: dobson = 2.685e20

! molecular masses in kg/mol
LOGICAL, PARAMETER  ::  l_ukca_diurnal_isopems   = .TRUE.
REAL, PARAMETER :: m_air        = 0.02897       ! Air
REAL, PARAMETER :: mmsul        = 0.09808       ! H2SO4
REAL, PARAMETER :: mmw          = 0.0180154     ! H2O


!  -------------------------------------------------------------------
!  Conversion factor from vmr to mmr for each species.
!             vmr*c_species = mmr
!             c_species = m_species/m_air  (m_air = 28.97)
!
!  Followed by molecular masses for each species (e.g. m_ch4) in g/mol
!  -------------------------------------------------------------------

REAL, PARAMETER :: c_o3p        = 0.5523
REAL, PARAMETER :: c_o1d        = 0.5523
REAL, PARAMETER :: c_o3         = 1.657
REAL, PARAMETER :: c_no         = 1.036
REAL, PARAMETER :: c_no3        = 2.140
REAL, PARAMETER :: c_no2        = 1.588
REAL, PARAMETER :: c_n2o5       = 3.728
REAL, PARAMETER :: c_ho2no2     = 2.727
REAL, PARAMETER :: c_hono2      = 2.175
REAL, PARAMETER :: c_hno3       = c_hono2   ! used in nitrate.F90
REAL, PARAMETER :: c_oh         = 0.5868
REAL, PARAMETER :: c_ho2        = 1.139
REAL, PARAMETER :: c_h2         = 0.06904
REAL, PARAMETER :: c_h2o2       = 1.174
REAL, PARAMETER :: c_ch4        = 0.5523
REAL, PARAMETER :: c_c          = 0.4142
REAL, PARAMETER :: c_co         = 0.9665
REAL, PARAMETER :: c_co2        = 1.5188
REAL, PARAMETER :: c_hcho       = 1.036
REAL, PARAMETER :: C_MeOO       = 1.622
REAL, PARAMETER :: c_h2o        = 0.6213
REAL, PARAMETER :: c_h2os       = 0.6213
REAL, PARAMETER :: C_MeOOH      = 1.657
REAL, PARAMETER :: c_hono       = 1.622
REAL, PARAMETER :: c_o2         = 1.105
REAL, PARAMETER :: c_n2         = 0.9665
REAL, PARAMETER :: c_c2h6       = 1.036
REAL, PARAMETER :: C_EtOO       = 2.106
REAL, PARAMETER :: C_EtOOH      = 2.140
REAL, PARAMETER :: C_MeCHO      = 1.519
REAL, PARAMETER :: c_toth       = 1.000
REAL, PARAMETER :: C_MeCO3      = 2.589
REAL, PARAMETER :: c_pan        = 4.177
REAL, PARAMETER :: c_c3h8       = 1.519
REAL, PARAMETER :: C_PrOO       = 2.589
REAL, PARAMETER :: C_PrOOH      = 2.623
REAL, PARAMETER :: C_EtCHO      = 2.002
REAL, PARAMETER :: C_EtCO3      = 3.072
REAL, PARAMETER :: C_Me2CO      = 2.002
REAL, PARAMETER :: C_MeCOCH2OO  = 3.072
REAL, PARAMETER :: C_MeCOCH2OOH = 3.107
REAL, PARAMETER :: c_ppan       = 4.660
REAL, PARAMETER :: C_MeONO2     = 2.658
REAL, PARAMETER :: c_n          = 0.48325
REAL, PARAMETER :: c_h          = 0.03452
REAL, PARAMETER :: c_n2o        = 1.5188
REAL, PARAMETER :: C_CFCl3      = 4.7480
REAL, PARAMETER :: C_CF2Cl2     = 4.1783
REAL, PARAMETER :: C_ClO        = 1.7784
REAL, PARAMETER :: C_HCl        = 1.2604
REAL, PARAMETER :: C_ClONO2     = 3.3668
REAL, PARAMETER :: C_HOCl       = 1.8129
REAL, PARAMETER :: C_OClO       = 2.3309
REAL, PARAMETER :: C_BrO        = 3.315
REAL, PARAMETER :: C_BrONO2     = 4.9034
REAL, PARAMETER :: C_HBr        = 2.7970
REAL, PARAMETER :: C_HOBr       = 3.3495
REAL, PARAMETER :: C_BrCl       = 3.9884
REAL, PARAMETER :: C_MeBr       = 3.2805
REAL, PARAMETER :: c_so2        = 2.2112
REAL, PARAMETER :: c_so3        = 2.7615
REAL, PARAMETER :: C_Me2S       = 2.145
REAL, PARAMETER :: c_dms        = 2.145
REAL, PARAMETER :: c_dmso       = 2.6965
REAL, PARAMETER :: c_ocs        = 2.0711
REAL, PARAMETER :: c_cos        = 2.0711
REAL, PARAMETER :: c_h2s        = 1.1766
REAL, PARAMETER :: c_cs2        = 2.6282
REAL, PARAMETER :: c_sad        = 4.1255
REAL, PARAMETER :: c_msa        = 3.317
REAL, PARAMETER :: c_s          = 1.1046
REAL, PARAMETER :: c_h2so4      = 3.385
REAL, PARAMETER :: C_CF2ClCFCl2 = 6.4722
REAL, PARAMETER :: C_CHF2Cl     = 2.9858
REAL, PARAMETER :: C_MeCCl3     = 4.6082
REAL, PARAMETER :: C_CCl4       = 5.3158
REAL, PARAMETER :: C_MeCl       = 1.7432
REAL, PARAMETER :: C_CF2ClBr    = 5.7128
REAL, PARAMETER :: C_CF3Br      = 5.1432
REAL, PARAMETER :: C_Cl         = 1.2261
REAL, PARAMETER :: C_Cl2O2      = 3.5568
REAL, PARAMETER :: C_Br         = 2.7627
REAL, PARAMETER :: C_CH2Br2     = 6.0013
REAL, PARAMETER :: c_mecf2cl    = 3.4673
REAL, PARAMETER :: c_cf2br2     = 7.2489
REAL, PARAMETER :: c_cf2brcf2br = 8.9748
REAL, PARAMETER :: c_cf2clcf3   = 5.3314
REAL, PARAMETER :: c_cf2clcf2cl = 5.8992
REAL, PARAMETER :: c_mecfcl2    = 4.0352
REAL, PARAMETER :: c_c5h8       = 2.3473
REAL, PARAMETER :: c_iso2       = 4.0387
REAL, PARAMETER :: c_isooh      = 4.0732
REAL, PARAMETER :: c_ison       = 5.3504  ! Might be revised to 5.0052
!                       for RAQ chem where ISON = (NO3)C4H6CHO: C5H7NO4: 145
REAL, PARAMETER :: c_macr       = 2.4163
REAL, PARAMETER :: c_macro2     = 4.1077
REAL, PARAMETER :: c_macrooh    = 4.1422
REAL, PARAMETER :: c_mpan       = 5.0742
REAL, PARAMETER :: c_hacet      = 2.5544
REAL, PARAMETER :: c_mgly       = 2.4853
REAL, PARAMETER :: c_nald       = 3.6244
REAL, PARAMETER :: c_hcooh      = 1.5878
REAL, PARAMETER :: c_meco3h     = 2.6234
REAL, PARAMETER :: c_meco2h     = 2.0711
REAL, PARAMETER :: c_nh3        = 0.5879
REAL, PARAMETER :: c_monoterp   = 4.7034
REAL, PARAMETER :: c_sec_org    = 5.1782    ! Molecular weight=150.

!     species required for ukca_scenario_rcp
REAL, PARAMETER :: C_cf3chf2    = 4.1429 ! HFC125
REAL, PARAMETER :: C_ch2fcf3    = 3.5219 ! HFC134a
REAL, PARAMETER :: C_chbr3      = 8.7238


!     Extra species for RAQ chemistry
REAL, PARAMETER :: c_c4h10      = 2.0021
REAL, PARAMETER :: c_c2h4       = 0.9665
REAL, PARAMETER :: c_c3h6       = 1.4498
REAL, PARAMETER :: c_rnc2h4     = 3.6244
REAL, PARAMETER :: c_rnc3h6     = 4.1077
!     rnc2h4 & rnc3h6 are CH2(NO3)CHO & CH3CH(NO3)CHO
REAL, PARAMETER :: c_ch3oh      = 1.1046
REAL, PARAMETER :: c_meoh       = 1.1046
REAL, PARAMETER :: c_toluene    = 3.1757
REAL, PARAMETER :: c_oxylene    = 3.6590
REAL, PARAMETER :: c_memald     = 3.3828
!     MEMALDIAL is CH3-CO-CH=CH-CHO, ring fragmentation
!     product from degradation of toluene and o-xylene
REAL, PARAMETER :: c_buooh      = 3.1067
REAL, PARAMETER :: c_mek        = 2.4853
REAL, PARAMETER :: c_mvk        = 2.4163
REAL, PARAMETER :: c_mvkooh     = 4.1422
REAL, PARAMETER :: c_gly        = 2.0021
REAL, PARAMETER :: c_rnc5h8     = 5.0052
REAL, PARAMETER :: c_orgnit     = 5.5230
REAL, PARAMETER :: c_buoo       = 3.0721
REAL, PARAMETER :: c_meko2      = 3.5554
REAL, PARAMETER :: c_hoc2h4o2   = 2.6579
REAL, PARAMETER :: c_hoc3h6o2   = 3.1412
REAL, PARAMETER :: c_hoipo2     = 4.0387
REAL, PARAMETER :: c_tolp1      = 3.7280
REAL, PARAMETER :: c_oxyl1      = 4.2113
REAL, PARAMETER :: c_memald1    = 3.9351
REAL, PARAMETER :: c_homvko2    = 4.1077

!     Extra species for EXTTC chemistry
REAL, PARAMETER :: c_isoo       = 4.0387
REAL, PARAMETER :: c_macroo     = 4.1077
REAL, PARAMETER :: c_apin       = 4.9645
REAL, PARAMETER :: C_PropeOO    = 2.5198
REAL, PARAMETER :: C_PropeOOH   = 2.5544
REAL, PARAMETER :: C_OnitU      = 3.5554
REAL, PARAMETER :: c_mekoo      = 3.5554
REAL, PARAMETER :: c_mekooh     = 3.5889
REAL, PARAMETER :: C_EteOO      = 2.0366
REAL, PARAMETER :: c_alka       = 2.0021
REAL, PARAMETER :: c_alkaoo     = 3.0721
REAL, PARAMETER :: c_alkaooh    = 3.1067
REAL, PARAMETER :: c_arom       = 3.4173
REAL, PARAMETER :: c_aromoo     = 4.4874
REAL, PARAMETER :: c_aromooh    = 4.5219
REAL, PARAMETER :: c_bsvoc1     = 4.9645   ! as APIN
REAL, PARAMETER :: c_bsvoc2     = 4.9645   ! as APIN
REAL, PARAMETER :: C_B2ndry     = 4.9645   ! as APIN
REAL, PARAMETER :: c_asvoc1     = 3.4173   ! as AROM
REAL, PARAMETER :: c_asvoc2     = 3.4173   ! as AROM
REAL, PARAMETER :: C_A2ndry     = 3.4173   ! as AROM
REAL, PARAMETER :: c_bsoa       = 5.1778   ! 150.0
REAL, PARAMETER :: c_asoa       = 5.1778   ! 150.0
REAL, PARAMETER :: c_isosvoc1   = 2.3473   ! as C5H8
REAL, PARAMETER :: c_isosvoc2   = 2.3473   ! as C5H8
REAL, PARAMETER :: c_isosoa     = 4.4874   ! 130.0

!     Extra species for CRI chemistry
REAL, PARAMETER :: c_noa        = 4.110
REAL, PARAMETER :: c_meo2no2    = 3.212
REAL, PARAMETER :: c_etono2     = 3.143
REAL, PARAMETER :: c_prono2     = 3.628
REAL, PARAMETER :: c_c2h2       = 0.899
REAL, PARAMETER :: c_benzene    = 2.696
REAL, PARAMETER :: c_tbut2ene   = 1.937
REAL, PARAMETER :: c_proh       = 2.074
REAL, PARAMETER :: c_etoh       = 1.590
REAL, PARAMETER :: c_etco3h     = 3.109
REAL, PARAMETER :: c_hoch2co3   = 3.143
REAL, PARAMETER :: c_hoch2co3h  = 3.177
REAL, PARAMETER :: c_hoch2ch2o2 = 2.660
REAL, PARAMETER :: c_hoch2cho   = 2.073
REAL, PARAMETER :: c_hoc2h4ooh  = 2.695
REAL, PARAMETER :: c_hoc2h4no3  = 3.696
REAL, PARAMETER :: c_phan       = 4.731
REAL, PARAMETER :: c_ch3sch2oo  = 3.215
REAL, PARAMETER :: c_ch3s       = 1.626
REAL, PARAMETER :: c_ch3so      = 2.178
REAL, PARAMETER :: c_ch3so2     = 2.730
REAL, PARAMETER :: c_ch3so3     = 3.283
REAL, PARAMETER :: c_msia       = 2.765

! CRI intermediate species with well defined masses that are wet
! and/or dry deposited
REAL, PARAMETER :: c_rn13no3    = 4.112
REAL, PARAMETER :: c_rn16no3    = 4.625
REAL, PARAMETER :: c_rn19no3    = 5.074
REAL, PARAMETER :: c_ra13no3    = 6.528
REAL, PARAMETER :: c_ra16no3    = 7.007
REAL, PARAMETER :: c_ra19no3    = 7.491
REAL, PARAMETER :: c_rtx24no3   = 6.800
REAL, PARAMETER :: c_rn9no3     = 4.180
REAL, PARAMETER :: c_rn12no3    = 4.664
REAL, PARAMETER :: c_rn15no3    = 5.147
REAL, PARAMETER :: c_rn18no3    = 5.631
REAL, PARAMETER :: c_ru14no3    = 5.079
REAL, PARAMETER :: c_rtn28no3   = 7.430
REAL, PARAMETER :: c_rtn25no3   = 6.945
REAL, PARAMETER :: c_rtn23no3   = 8.050
REAL, PARAMETER :: c_rtx28no3   = 7.428
REAL, PARAMETER :: c_rtx22no3   = 7.359
REAL, PARAMETER :: c_rn16ooh    = 3.594
REAL, PARAMETER :: c_rn19ooh    = 4.077
REAL, PARAMETER :: c_rn8ooh     = 3.109
REAL, PARAMETER :: c_rn11ooh    = 3.593
REAL, PARAMETER :: c_rn14ooh    = 4.077
REAL, PARAMETER :: c_rn17ooh    = 4.56
REAL, PARAMETER :: c_nru14ooh   = 5.631
REAL, PARAMETER :: c_nru12ooh   = 6.735
REAL, PARAMETER :: c_rn9ooh     = 3.179
REAL, PARAMETER :: c_rn12ooh    = 3.663
REAL, PARAMETER :: c_rn15ooh    = 4.146
REAL, PARAMETER :: c_rn18ooh    = 4.63
REAL, PARAMETER :: c_nrn6ooh    = 4.248
REAL, PARAMETER :: c_nrn9ooh    = 4.732
REAL, PARAMETER :: c_nrn12ooh   = 5.216
REAL, PARAMETER :: c_ra13ooh    = 5.527
REAL, PARAMETER :: c_ra16ooh    = 6.01
REAL, PARAMETER :: c_ra19ooh    = 6.493
REAL, PARAMETER :: c_rtn28ooh   = 6.429
REAL, PARAMETER :: c_nrtn28ooh  = 7.982
REAL, PARAMETER :: c_rtn26ooh   = 6.912
REAL, PARAMETER :: c_rtn25ooh   = 5.944
REAL, PARAMETER :: c_rtn24ooh   = 6.496
REAL, PARAMETER :: c_rtn23ooh   = 7.049
REAL, PARAMETER :: c_rtn14ooh   = 5.561
REAL, PARAMETER :: c_rtn10ooh   = 5.043
REAL, PARAMETER :: c_rtx28ooh   = 6.429
REAL, PARAMETER :: c_rtx24ooh   = 5.875
REAL, PARAMETER :: c_rtx22ooh   = 6.358
REAL, PARAMETER :: c_nrtx28ooh  = 7.982
REAL, PARAMETER :: c_carb14     = 2.973
REAL, PARAMETER :: c_carb17     = 3.457
REAL, PARAMETER :: c_carb11a    = 2.489
REAL, PARAMETER :: c_carb10     = 3.041
REAL, PARAMETER :: c_carb13     = 3.525
REAL, PARAMETER :: c_carb16     = 4.008
REAL, PARAMETER :: c_carb9      = 2.972
REAL, PARAMETER :: c_carb12     = 3.456
REAL, PARAMETER :: c_carb15     = 3.940
REAL, PARAMETER :: c_ucarb12    = 3.456
REAL, PARAMETER :: c_nucarb12   = 5.009
REAL, PARAMETER :: c_udcarb8    = 2.902
REAL, PARAMETER :: c_udcarb11   = 3.385
REAL, PARAMETER :: c_udcarb14   = 3.868
REAL, PARAMETER :: c_tncarb26   = 5.807
REAL, PARAMETER :: c_tncarb10   = 3.938
REAL, PARAMETER :: c_tncarb12   = 4.974
REAL, PARAMETER :: c_tncarb11   = 4.905
REAL, PARAMETER :: c_ccarb12    = 4.349
REAL, PARAMETER :: c_tncarb15   = 3.383
REAL, PARAMETER :: c_rcooh25    = 6.358
REAL, PARAMETER :: c_txcarb24   = 4.771
REAL, PARAMETER :: c_txcarb22   = 5.253
REAL, PARAMETER :: c_ru12pan    = 6.114
REAL, PARAMETER :: c_rtn26pan   = 8.465
REAL, PARAMETER :: c_aroh14     = 3.249
REAL, PARAMETER :: c_aroh17     = 3.733
REAL, PARAMETER :: c_arnoh14    = 4.802
REAL, PARAMETER :: c_arnoh17    = 5.286
REAL, PARAMETER :: c_anhy       = 3.385
REAL, PARAMETER :: c_dhpcarb9   = 4.698
REAL, PARAMETER :: c_hpucarb12  = 4.008
REAL, PARAMETER :: c_hucarb9    = 2.972
REAL, PARAMETER :: c_iepox      = 4.078
REAL, PARAMETER :: c_hmml       = 3.524
REAL, PARAMETER :: c_dhpr12ooh  = 6.287
REAL, PARAMETER :: c_dhcarb9    = 3.594
REAL, PARAMETER :: c_ru12no3    = 6.183
REAL, PARAMETER :: c_ru10no3    = 5.147
REAL, PARAMETER :: c_dhpr12o2   = 6.252
REAL, PARAMETER :: c_ru10ao2    = 4.111
REAL, PARAMETER :: c_maco3      = 3.489


! Define conversion factor for all other common representative
! intermediate (CRI) species, as these do not have defined
! molecular masses.
! - Using m_cri = 150g/mol ~ Rough average of intermediate species
!                            Same as Sec_org
REAL, PARAMETER :: c_cri        = 5.1782

!     molecular masses in g/mol of emitted species,
!     for budget calculations
REAL, PARAMETER :: m_ho2     =  33.007
REAL, PARAMETER :: m_ch4     =  16.0
REAL, PARAMETER :: m_co      =  28.0
REAL, PARAMETER :: m_hcho    =  30.0
REAL, PARAMETER :: m_c2h6    =  30.0
REAL, PARAMETER :: m_c3h8    =  44.0
REAL, PARAMETER :: m_mecho   =  44.0
REAL, PARAMETER :: m_no2     =  46.0
REAL, PARAMETER :: m_n2o5    = 108.01
REAL, PARAMETER :: m_me2co   =  58.0
REAL, PARAMETER :: m_isop    =  68.0
REAL, PARAMETER :: m_no      =  30.0
REAL, PARAMETER :: m_n       =  14.0
REAL, PARAMETER :: m_c       =  12.0
REAL, PARAMETER :: m_monoterp =  136.24

!     molecular masses of stratospheric species, for which surface
!     mmrs are prescribed
REAL, PARAMETER :: m_hcl        =  36.5
REAL, PARAMETER :: m_n2o        =  44.0
REAL, PARAMETER :: m_clo        =  51.5
REAL, PARAMETER :: m_hocl       =  52.5
REAL, PARAMETER :: m_oclo       =  67.5
REAL, PARAMETER :: m_clono2     =  97.5
REAL, PARAMETER :: m_cf2cl2     = 121.0
REAL, PARAMETER :: m_cfcl3      = 137.5
REAL, PARAMETER :: m_hbr        =  81.0
REAL, PARAMETER :: m_mebr       =  95.0
REAL, PARAMETER :: m_bro        =  96.0
REAL, PARAMETER :: m_hobr       =  97.0
REAL, PARAMETER :: m_brcl       = 115.5
REAL, PARAMETER :: m_brono2     = 142.0
REAL, PARAMETER :: m_cf2clcfcl2 = 187.5
REAL, PARAMETER :: m_chf2cl     =  86.5
REAL, PARAMETER :: m_meccl3     = 133.5
REAL, PARAMETER :: m_ccl4       = 154.0
REAL, PARAMETER :: m_mecl       =  50.5
REAL, PARAMETER :: m_cf2clbr    = 165.5
REAL, PARAMETER :: m_cf3br      = 149.0
REAL, PARAMETER :: m_ch2br2     = 173.835

! sulphur containing, etc.
REAL, PARAMETER :: m_ocs        =  60.0
REAL, PARAMETER :: m_cos        =  60.0
REAL, PARAMETER :: m_h2s        =  34.086
REAL, PARAMETER :: m_cs2        =  76.14
REAL, PARAMETER :: m_dms        =  62.1
REAL, PARAMETER :: m_dmso       =  78.13
REAL, PARAMETER :: m_me2s       =  62.1
REAL, PARAMETER :: m_msa        =  96.1
REAL, PARAMETER :: m_sec_org    =  150.0
REAL, PARAMETER :: m_s          =  32.07
REAL, PARAMETER :: m_so2        =  64.06
REAL, PARAMETER :: m_so3        =  80.06
REAL, PARAMETER :: m_so4        =  96.06
REAL, PARAMETER :: m_h2so4      =  98.07
REAL, PARAMETER :: m_nh3        =  17.03
REAL, PARAMETER :: m_nh42so4    =  132.16

REAL, PARAMETER :: m_cl         = 35.5
REAL, PARAMETER :: m_cl2o2      = 103.0
REAL, PARAMETER :: m_br         = 80.0
REAL, PARAMETER :: m_h2         = 2.016
REAL, PARAMETER :: m_h2o        = 18.0
REAL, PARAMETER :: m_mecoch2ooh = 90.0
REAL, PARAMETER :: m_isooh      = 118.0
REAL, PARAMETER :: m_mpan       = 147.0
REAL, PARAMETER :: m_ppan       = 135.0
REAL, PARAMETER :: m_pan        = 121.0
REAL, PARAMETER :: m_hno3       = 63.0
REAL, PARAMETER :: m_hono2      = 63.0

!     Extra masses for RAQ or other chemistries
REAL, PARAMETER :: m_c5h8    =  68.0
REAL, PARAMETER :: m_c4h10   =  58.0
REAL, PARAMETER :: m_c2h4    =  28.0
REAL, PARAMETER :: m_c3h6    =  42.0
REAL, PARAMETER :: m_toluene =  92.0
REAL, PARAMETER :: m_oxylene = 106.0
REAL, PARAMETER :: m_ch3oh   =  32.0
REAL, PARAMETER :: m_meoh    =  32.0
REAL, PARAMETER :: m_buooh   =  90.0
REAL, PARAMETER :: m_mvkooh  = 120.0
REAL, PARAMETER :: m_orgnit  = 160.0
REAL, PARAMETER :: m_macrooh = 120.0
REAL, PARAMETER :: m_gly     =  58.0

REAL, PARAMETER :: m_hono = 47.0
!     Extra masses for Wesely scheme
REAL, PARAMETER :: m_macr    = 70.0
REAL, PARAMETER :: m_etcho   = 58.0
REAL, PARAMETER :: m_nald    = 105.0
REAL, PARAMETER :: m_mgly    = 72.0
REAL, PARAMETER :: m_hacet   = 74.0
REAL, PARAMETER :: m_hcooh   = 46.0
REAL, PARAMETER :: m_meco2h  = 60.0

!     Extra masses for EXTTC chemistry
REAL, PARAMETER :: m_apin     =  136.0
REAL, PARAMETER :: m_mvk      =  70.0
REAL, PARAMETER :: m_mek      =  72.0
REAL, PARAMETER :: m_alka     =  58.0        ! as butane
REAL, PARAMETER :: m_arom     =  99.0        ! (toluene + xylene)/2
REAL, PARAMETER :: m_bsvoc1   = 144.0
REAL, PARAMETER :: m_bsvoc2   = 144.0
REAL, PARAMETER :: m_asvoc1   = 99.0
REAL, PARAMETER :: m_asvoc2   = 99.0
REAL, PARAMETER :: m_isosvoc1 = 68.0
REAL, PARAMETER :: m_isosvoc2 = 68.0
REAL, PARAMETER :: m_onitu    = 102.0
REAL, PARAMETER :: m_bsoa     = 150.0
REAL, PARAMETER :: m_asoa     = 150.0
REAL, PARAMETER :: m_isosoa   = 130.0
REAL, PARAMETER :: m_alkaooh  = 90.0
REAL, PARAMETER :: m_aromooh  = 130.0
REAL, PARAMETER :: m_mekooh   = 104.0

!     The mass of organic nitrate is an approximation,
!     calculated as the average of ORGNIT formed by two
!     reacs. in UKCA_CHEMCO_RAQ:
!      NO2 + TOLP1 --> ORGNIT (A)
!      NO2 + OXYL1 --> ORGNIT (B)
!      * TOL  = methylbenzene       = C6H5(CH3)
!        OXYL = 1,2-dimethylbenzene = C6H4(CH3)2
!      * TOL  + OH --> TOLP1: C6H4(OH)(CH3)  = methyl phenol
!        OXYL + OH --> OXYL1: C6H3(OH)(CH3)2 = dimethyl phenol
!      * ORGNIT A: TOLP1 + NO2 ~ C6H3(CH3)(OH)NO2  ~
!                  C7H7NO3: methyl nitrophenol   -> 153
!        ORGNIT B: OXYL1 + NO2 ~ C6H2(CH3)2(OH)NO2 ~
!                  C8H9NO3: dimethyl nitrophenol -> 167
!  -------------------------------------------------------------------

!     Extra species for CRI chemistry
REAL, PARAMETER :: m_noa        = 119.08
REAL, PARAMETER :: m_meo2no2    = 93.039
REAL, PARAMETER :: m_etono2     = 91.066
REAL, PARAMETER :: m_prono2     = 105.09
REAL, PARAMETER :: m_c2h2       = 26.037
REAL, PARAMETER :: m_benzene    = 78.112
REAL, PARAMETER :: m_tbut2ene   = 56.106
REAL, PARAMETER :: m_proh       = 60.095
REAL, PARAMETER :: m_etoh       = 46.068
REAL, PARAMETER :: m_etco3h     = 90.078
REAL, PARAMETER :: m_hoch2co3   = 91.043
REAL, PARAMETER :: m_hoch2co3h  = 92.051
REAL, PARAMETER :: m_hoch2ch2o2 = 77.059
REAL, PARAMETER :: m_hoch2cho   = 60.052
REAL, PARAMETER :: m_hoc2h4ooh  = 78.067
REAL, PARAMETER :: m_hoc2h4no3  = 107.07
REAL, PARAMETER :: m_phan       = 137.05
REAL, PARAMETER :: m_ch3sch2oo  = 93.125
REAL, PARAMETER :: m_ch3s       = 47.01
REAL, PARAMETER :: m_ch3so      = 63.099
REAL, PARAMETER :: m_ch3so2     = 79.098
REAL, PARAMETER :: m_ch3so3     = 95.098
REAL, PARAMETER :: m_msia       = 80.106
! CRI intermediate species with well defined masses:
REAL, PARAMETER :: m_rn13no3    = 119.1
REAL, PARAMETER :: m_rn16no3    = 134.0
REAL, PARAMETER :: m_rn19no3    = 147.0
REAL, PARAMETER :: m_ra13no3    = 189.1
REAL, PARAMETER :: m_ra16no3    = 203.1
REAL, PARAMETER :: m_ra19no3    = 217.1
REAL, PARAMETER :: m_rtx24no3   = 197.1
REAL, PARAMETER :: m_hoc2h4no2  = 107.0
REAL, PARAMETER :: m_rn9no3     = 121.1
REAL, PARAMETER :: m_rn12no3    = 135.1
REAL, PARAMETER :: m_rn15no3    = 149.1
REAL, PARAMETER :: m_rn18no3    = 163.1
REAL, PARAMETER :: m_ru14no3    = 147.1
REAL, PARAMETER :: m_rtn28no3   = 215.3
REAL, PARAMETER :: m_rtn25no3   = 201.2
REAL, PARAMETER :: m_rtn23no3   = 233.2
REAL, PARAMETER :: m_rtx28no3   = 215.3
REAL, PARAMETER :: m_rtx22no3   = 213.2
REAL, PARAMETER :: m_rn16ooh    = 104.1
REAL, PARAMETER :: m_rn19ooh    = 118.1
REAL, PARAMETER :: m_rn14ooh    = 118.1 ! ~ m_isooh
REAL, PARAMETER :: m_rn17ooh    = 132.1
REAL, PARAMETER :: m_nru14ooh   = 163.1
REAL, PARAMETER :: m_nru12ooh   = 195.1
REAL, PARAMETER :: m_rn9ooh     = 92.09
REAL, PARAMETER :: m_rn12ooh    = 106.1
REAL, PARAMETER :: m_rn15ooh    = 120.1
REAL, PARAMETER :: m_rn18ooh    = 134.1
REAL, PARAMETER :: m_nrn6ooh    = 123.1
REAL, PARAMETER :: m_nrn9ooh    = 137.1
REAL, PARAMETER :: m_nrn12ooh   = 151.1
REAL, PARAMETER :: m_ra13ooh    = 160.1
REAL, PARAMETER :: m_ra16ooh    = 174.1
REAL, PARAMETER :: m_ra19ooh    = 188.1
REAL, PARAMETER :: m_rtn28ooh   = 186.3
REAL, PARAMETER :: m_nrtn28ooh  = 231.3
REAL, PARAMETER :: m_rtn26ooh   = 200.2
REAL, PARAMETER :: m_rtn25ooh   = 172.2
REAL, PARAMETER :: m_rtn24ooh   = 188.2
REAL, PARAMETER :: m_rtn23ooh   = 204.2
REAL, PARAMETER :: m_rtn14ooh   = 161.1
REAL, PARAMETER :: m_rtn10ooh   = 146.1
REAL, PARAMETER :: m_rtx28ooh   = 186.2
REAL, PARAMETER :: m_rtx24ooh   = 170.2
REAL, PARAMETER :: m_rtx22ooh   = 184.2
REAL, PARAMETER :: m_nrtx28ooh  = 231.2
REAL, PARAMETER :: m_carb14     = 86.13
REAL, PARAMETER :: m_carb17     = 100.2
REAL, PARAMETER :: m_carb11a    = 72.11
REAL, PARAMETER :: m_carb10     = 88.11
REAL, PARAMETER :: m_carb13     = 102.1
REAL, PARAMETER :: m_carb16     = 116.1
REAL, PARAMETER :: m_carb9      = 86.09
REAL, PARAMETER :: m_carb12     = 100.1
REAL, PARAMETER :: m_carb15     = 114.1
REAL, PARAMETER :: m_ucarb12    = 100.1
REAL, PARAMETER :: m_nucarb12   = 145.1
REAL, PARAMETER :: m_udcarb8    = 84.07
REAL, PARAMETER :: m_udcarb11   = 98.07
REAL, PARAMETER :: m_udcarb14   = 112.1
REAL, PARAMETER :: m_tncarb26   = 168.2
REAL, PARAMETER :: m_tncarb10   = 114.1
REAL, PARAMETER :: m_tncarb12   = 144.1
REAL, PARAMETER :: m_tncarb11   = 142.1
REAL, PARAMETER :: m_ccarb12    = 126.0
REAL, PARAMETER :: m_tncarb15   = 98.0
REAL, PARAMETER :: m_rcooh25    = 184.2
REAL, PARAMETER :: m_txcarb24   = 138.2
REAL, PARAMETER :: m_txcarb22   = 152.2
REAL, PARAMETER :: m_ru12pan    = 177.1
REAL, PARAMETER :: m_rtn26pan   = 245.2
REAL, PARAMETER :: m_anhy       = 98.06 ! = MALANHY
REAL, PARAMETER :: m_aroh14     = 94.11 ! = PHENOL
REAL, PARAMETER :: m_aroh17     = 108.1 ! = CRESOL
REAL, PARAMETER :: m_arnoh14    = 139.1
REAL, PARAMETER :: m_arnoh17    = 153.1
REAL, PARAMETER :: m_dhpcarb9   = 136.102
REAL, PARAMETER :: m_hpucarb12  = 116.115
REAL, PARAMETER :: m_hucarb9    = 86.089
REAL, PARAMETER :: m_iepox      = 118.131
REAL, PARAMETER :: m_hmml       = 102.089
REAL, PARAMETER :: m_dhpr12ooh  = 182.127
REAL, PARAMETER :: m_dhcarb9    = 104.105
REAL, PARAMETER :: m_ru12no3    = 179.128
REAL, PARAMETER :: m_ru10no3    = 149.102
REAL, PARAMETER :: m_dhpr12o2   = 181.119
REAL, PARAMETER :: m_ru10ao2    = 119.096
REAL, PARAMETER :: m_maco3      = 101.081


!   Define dummy molar mass for all common representative intermediate species
!   in CRI mechanism, as these do not have defined molecular masses.
!   - Using m_cri = 150g/mol (same as Sec_org)
REAL, PARAMETER :: m_cri = 150.0

! For solar calculations
REAL :: fxb  ! Latitude of tropic of Capricorn (radians)
REAL :: fxc  ! Hours per day / pi

! For ACTIVATE, used in hygroscopy and curvature parameter

REAL, PARAMETER :: zsten = 75.0e-3 ! surface tension of H2O [J m-2]
!                                  !   neglecting salts and temperature
                                   ! from Vargaftik et al, JPCRD 1983
REAL, PARAMETER :: zosm = 1.0      ! Osmotic coefficient, currently fixed


CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE set_derived_constants()
! ----------------------------------------------------------------------
! Description:
!   Sets values of UKCA constants derived from other constants that are
!   potentially user configurable (so cannot be determined at compile time)
! ----------------------------------------------------------------------

USE ukca_um_legacy_mod, ONLY: rhour_per_day, pi, recip_pi_over_180

IMPLICIT NONE

fxb = 23.45 / recip_pi_over_180
fxc = rhour_per_day / pi

RETURN
END SUBROUTINE set_derived_constants

END MODULE ukca_constants
