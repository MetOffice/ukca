! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module to contain tracer and species numbers
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: Fortran
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_cspecies

USE yomhook,         ONLY: lhook, dr_hook
USE parkind1,        ONLY: jprb, jpim
IMPLICIT NONE
PRIVATE

! Index location in tracer array
INTEGER, SAVE, PUBLIC :: n_ox    ! tracer numbers
INTEGER, SAVE, PUBLIC :: n_o3
INTEGER, SAVE, PUBLIC :: n_o3s
INTEGER, SAVE, PUBLIC :: n_nox
INTEGER, SAVE, PUBLIC :: n_no
INTEGER, SAVE, PUBLIC :: n_no2
INTEGER, SAVE, PUBLIC :: n_no3
INTEGER, SAVE, PUBLIC :: n_n2o5
INTEGER, SAVE, PUBLIC :: n_ho2no2
INTEGER, SAVE, PUBLIC :: n_hono2
INTEGER, SAVE, PUBLIC :: n_pan
INTEGER, SAVE, PUBLIC :: n_ison
INTEGER, SAVE, PUBLIC :: n_orgnit
INTEGER, SAVE, PUBLIC :: n_rnc2h4
INTEGER, SAVE, PUBLIC :: n_rnc3h6
INTEGER, SAVE, PUBLIC :: n_h2o2
INTEGER, SAVE, PUBLIC :: n_ch4
INTEGER, SAVE, PUBLIC :: n_sx
INTEGER, SAVE, PUBLIC :: n_h2
INTEGER, SAVE, PUBLIC :: n_h2o
INTEGER, SAVE, PUBLIC :: n_cl
INTEGER, SAVE, PUBLIC :: n_ox1
INTEGER, SAVE, PUBLIC :: n_sx1
INTEGER, SAVE, PUBLIC :: n_dms
INTEGER, SAVE, PUBLIC :: n_so2
INTEGER, SAVE, PUBLIC :: n_so3
INTEGER, SAVE, PUBLIC :: n_h2so4
INTEGER, SAVE, PUBLIC :: n_sec_org
INTEGER, SAVE, PUBLIC :: n_sec_org_i  ! Secondary organic from isoprene
INTEGER, SAVE, PUBLIC :: n_cfcl3    ! CFC-11
INTEGER, SAVE, PUBLIC :: n_cf2cl2   ! CFC-12
INTEGER, SAVE, PUBLIC :: n_bro
INTEGER, SAVE, PUBLIC :: n_hcl
INTEGER, SAVE, PUBLIC :: n_o1d
INTEGER, SAVE, PUBLIC :: n_passive ! passive o3
INTEGER, SAVE, PUBLIC :: n_hcho
INTEGER, SAVE, PUBLIC :: n_c2h6
INTEGER, SAVE, PUBLIC :: n_mecho
INTEGER, SAVE, PUBLIC :: n_c3h8
INTEGER, SAVE, PUBLIC :: n_me2co
INTEGER, SAVE, PUBLIC :: n_c5h8
INTEGER, SAVE, PUBLIC :: n_mgly
INTEGER, SAVE, PUBLIC :: n_c4h10
INTEGER, SAVE, PUBLIC :: n_mek
INTEGER, SAVE, PUBLIC :: n_c2h4
INTEGER, SAVE, PUBLIC :: n_c3h6
INTEGER, SAVE, PUBLIC :: n_oxylene
INTEGER, SAVE, PUBLIC :: n_ch3oh
INTEGER, SAVE, PUBLIC :: n_gly
INTEGER, SAVE, PUBLIC :: n_mvk
INTEGER, SAVE, PUBLIC :: n_toluene
INTEGER, SAVE, PUBLIC :: n_co       ! CO
INTEGER, SAVE, PUBLIC :: n_n2o      ! N2O
INTEGER, SAVE, PUBLIC :: n_mebr     ! CH3Br
INTEGER, SAVE, PUBLIC :: n_nh3      ! Ammonia

! Index location in species array
INTEGER, SAVE, PUBLIC :: nn_o3    ! Species numbers
INTEGER, SAVE, PUBLIC :: nn_o3s
INTEGER, SAVE, PUBLIC :: nn_oh
INTEGER, SAVE, PUBLIC :: nn_ho2
INTEGER, SAVE, PUBLIC :: nn_h2o2
INTEGER, SAVE, PUBLIC :: nn_no
INTEGER, SAVE, PUBLIC :: nn_no2
INTEGER, SAVE, PUBLIC :: nn_o1d
INTEGER, SAVE, PUBLIC :: nn_o3p
INTEGER, SAVE, PUBLIC :: nn_meoo
INTEGER, SAVE, PUBLIC :: nn_meco3
INTEGER, SAVE, PUBLIC :: nn_etoo
INTEGER, SAVE, PUBLIC :: nn_etco3
INTEGER, SAVE, PUBLIC :: nn_nproo
INTEGER, SAVE, PUBLIC :: nn_nprooh
INTEGER, SAVE, PUBLIC :: nn_iProo
INTEGER, SAVE, PUBLIC :: nn_iProoh
INTEGER, SAVE, PUBLIC :: nn_mecoch2oo
INTEGER, SAVE, PUBLIC :: nn_no3
INTEGER, SAVE, PUBLIC :: nn_n2o5
INTEGER, SAVE, PUBLIC :: nn_ho2no2
INTEGER, SAVE, PUBLIC :: nn_hono2
INTEGER, SAVE, PUBLIC :: nn_ch4
INTEGER, SAVE, PUBLIC :: nn_so2
INTEGER, SAVE, PUBLIC :: nn_so3
INTEGER, SAVE, PUBLIC :: nn_h2so4
INTEGER, SAVE, PUBLIC :: nn_ohs
INTEGER, SAVE, PUBLIC :: nn_ho2s
INTEGER, SAVE, PUBLIC :: nn_o1ds
INTEGER, SAVE, PUBLIC :: nn_o3ps
INTEGER, SAVE, PUBLIC :: nn_n
INTEGER, SAVE, PUBLIC :: nn_cl
INTEGER, SAVE, PUBLIC :: nn_clo
INTEGER, SAVE, PUBLIC :: nn_hcl
INTEGER, SAVE, PUBLIC :: nn_cl2o2
INTEGER, SAVE, PUBLIC :: nn_meo
INTEGER, SAVE, PUBLIC :: nn_bro
INTEGER, SAVE, PUBLIC :: nn_br
INTEGER, SAVE, PUBLIC :: nn_buoo        !for RAQ chemistry
INTEGER, SAVE, PUBLIC :: nn_meko2
INTEGER, SAVE, PUBLIC :: nn_hoc2h4o2
INTEGER, SAVE, PUBLIC :: nn_hoc3h6o2
INTEGER, SAVE, PUBLIC :: nn_oxyl1
INTEGER, SAVE, PUBLIC :: nn_memald1
INTEGER, SAVE, PUBLIC :: nn_hoipo2
INTEGER, SAVE, PUBLIC :: nn_homvko2
INTEGER, SAVE, PUBLIC :: nn_tolp1
INTEGER, SAVE, PUBLIC :: nn_cfcl3
INTEGER, SAVE, PUBLIC :: nn_cf2cl2
INTEGER, SAVE, PUBLIC :: nn_nh3

! Advected tracers
REAL, SAVE, ALLOCATABLE, PUBLIC :: c_species(:)
! Non-advected tracers
REAL, SAVE, ALLOCATABLE, PUBLIC :: c_na_species(:)

PUBLIC ukca_calc_cspecies

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CSPECIES'

CONTAINS

SUBROUTINE ukca_calc_cspecies

USE ukca_constants,  ONLY: c_br, c_brcl, c_bro, c_brono2, c_buoo,              &
    c_buooh, c_c2h4, c_c2h6, c_c3h6, c_c3h8, c_c4h10, c_c5h8, c_ccl4,          &
    c_cf2cl2, c_cf2clbr, c_cf2clcfcl2, c_cf3br, c_cfcl3, c_ch2br2,             &
    c_ch3oh, c_ch4, c_chf2cl, c_cl, c_cl2o2, c_clo, c_clono2, c_co,            &
    c_cos, c_cs2, c_dms, c_dmso, c_etcho, c_etco3, c_etoo, c_etooh,            &
    c_gly, c_h, c_h2, c_h2o, c_h2o2, c_h2s, c_h2so4, c_hacet, c_hbr,           &
    c_hcho, c_hcl, c_hcooh, c_ho2, c_ho2no2, c_hobr, c_hoc2h4o2,               &
    c_hoc3h6o2, c_hocl, c_hoipo2, c_homvko2, c_hono, c_hono2, c_iso2,          &
    c_ison, c_isooh, c_macr, c_macro2, c_macrooh, c_me2co, c_me2s,             &
    c_mebr, c_meccl3, c_mecho, c_mecl, c_meco2h, c_meco3, c_meco3h,            &
    c_mecoch2oo, c_mecoch2ooh, c_mek, c_meko2, c_memald, c_memald1,            &
    c_meoh, c_meono2, c_meoo, c_meooh, c_mgly, c_monoterp, c_mpan,             &
    c_msa, c_mvk, c_mvkooh, c_n, c_n2o, c_n2o5, c_nald, c_nh3, c_no,           &
    c_no2, c_no3, c_o1d, c_o3, c_o3p, c_oclo, c_oh, c_orgnit,                  &
    c_oxyl1, c_oxylene, c_pan, c_ppan, c_proo, c_prooh, c_rnc2h4,              &
    c_rnc3h6, c_sec_org, c_so2, c_so3, c_tolp1, c_toluene, c_noa,              &
    c_rnc3h6, c_sec_org, c_sec_org_i, c_so2, c_so3, c_tolp1, c_toluene, c_noa, &
    c_meo2no2, c_etono2, c_prono2, c_c2h2, c_benzene, c_tbut2ene,              &
    c_proh, c_etoh, c_etco3h, c_hoch2co3, c_hoch2co3h, c_hoch2ch2o2,           &
    c_hoch2cho, c_hoc2h4ooh, c_hoc2h4no3, c_phan, c_ch3sch2oo,                 &
    c_ch3s, c_ch3so, c_ch3so2, c_ch3so3, c_msia, c_alkaooh, c_mekooh,          &
    c_rn13no3, c_rn16no3, c_rn19no3, c_ra13no3, c_ra16no3, c_ra19no3,          &
    c_rtx24no3, c_rn9no3, c_rn12no3, c_rn15no3, c_rn18no3, c_ru14no3,          &
    c_rtn28no3, c_rtn25no3, c_rtn23no3, c_rtx28no3, c_rtx22no3, c_rn16ooh,     &
    c_rn19ooh, c_rn14ooh, c_rn17ooh, c_nru14ooh, c_nru12ooh, c_rn9ooh,         &
    c_rn12ooh, c_rn15ooh, c_rn18ooh, c_nrn6ooh, c_nrn9ooh, c_nrn12ooh,         &
    c_ra13ooh, c_ra16ooh, c_ra19ooh, c_rtn28ooh, c_nrtn28ooh, c_rtn26ooh,      &
    c_rtn25ooh, c_rtn24ooh, c_rtn23ooh, c_rtn14ooh, c_rtn10ooh, c_rtx28ooh,    &
    c_rtx24ooh, c_rtx22ooh, c_nrtx28ooh, c_carb14, c_carb17, c_carb11a,        &
    c_carb10, c_carb13, c_carb16, c_carb9, c_carb12, c_carb15, c_ucarb12,      &
    c_nucarb12, c_udcarb8, c_udcarb11, c_udcarb14, c_tncarb26, c_tncarb10,     &
    c_tncarb12, c_tncarb11, c_ccarb12, c_tncarb15, c_rcooh25, c_txcarb24,      &
    c_txcarb22, c_ru12pan, c_rtn26pan, c_aroh14, c_aroh17, c_arnoh14,          &
    c_arnoh17, c_anhy, c_cri,c_dhpcarb9, c_hpucarb12, c_hucarb9, c_iepox,      &
    c_hmml, c_dhpr12ooh, c_dhcarb9, c_ru12no3, c_ru10no3, c_dhpr12o2,          &
    c_ru10ao2, c_maco3


USE asad_mod,        ONLY: advt, nadvt, speci, jpctr, jpspec
USE ukca_config_specification_mod, ONLY: ukca_config

USE ereport_mod,     ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER :: errcode                ! Variable passed to ereport
INTEGER :: i, m
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALC_CSPECIES'


!     Compute list of conversion factors from vmr to mmr.
!     Add new entry if new tracers are introduced, and add value for c_xx
!      in UKCA_CONSTANTS, where xx is the species.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! advected tracers
c_species = 0.0
WHERE (advt == 'Ox        ') c_species = c_o3
WHERE (advt == 'O3        ') c_species = c_o3
WHERE (advt == 'NO        ') c_species = c_no
WHERE (advt == 'NO2       ') c_species = c_no2
WHERE (advt == 'NO3       ') c_species = c_no3
WHERE (advt == 'NOx       ') c_species = c_no2
WHERE (advt == 'N2O5      ') c_species = c_n2o5
WHERE (advt == 'HO2NO2    ') c_species = c_ho2no2
WHERE (advt == 'HONO2     ') c_species = c_hono2
WHERE (advt == 'H2O2      ') c_species = c_h2o2
WHERE (advt == 'CH4       ') c_species = c_ch4
WHERE (advt == 'CO        ') c_species = c_co
WHERE (advt == 'HCHO      ') c_species = c_hcho
WHERE (advt == 'MeOOH     ') c_species = c_meooh
WHERE (advt == 'HONO      ') c_species = c_hono
WHERE (advt == 'C2H6      ') c_species = c_c2h6
WHERE (advt == 'EtOOH     ') c_species = c_etooh
WHERE (advt == 'MeCHO     ') c_species = c_mecho
WHERE (advt == 'PAN       ') c_species = c_pan
WHERE (advt == 'C3H8      ') c_species = c_c3h8
WHERE (advt == 'i-PrOOH   ') c_species = c_prooh
WHERE (advt == 'n-PrOOH   ') c_species = c_prooh
WHERE (advt == 'EtCHO     ') c_species = c_etcho
WHERE (advt == 'Me2CO     ') c_species = c_me2co
WHERE (advt == 'MeCOCH2OOH') c_species = c_mecoch2ooh
WHERE (advt == 'PPAN      ') c_species = c_ppan
WHERE (advt == 'MeONO2    ') c_species = c_meono2
WHERE (advt == 'Sx        ') c_species = c_o3
WHERE (advt == 'O3S       ') c_species = c_o3
WHERE (advt == 'HOx       ') c_species = c_ho2
WHERE (advt == 'N2O       ') c_species = c_n2o
WHERE (advt == 'CFCl3     ') c_species = c_cfcl3
WHERE (advt == 'CF2Cl2    ') c_species = c_cf2cl2
WHERE (advt == 'H2O       ') c_species = c_h2o
WHERE (advt == 'ClONO2    ') c_species = c_clono2
WHERE (advt == 'Clx       ') c_species = c_clo
WHERE (advt == 'ClO       ') c_species = c_clo
WHERE (advt == 'Cl2O2     ') c_species = c_cl2o2
WHERE (advt == 'Cl        ') c_species = c_cl
WHERE (advt == 'HCl       ') c_species = c_hcl
WHERE (advt == 'HOCl      ') c_species = c_hocl
WHERE (advt == 'OClO      ') c_species = c_oclo
WHERE (advt == 'Br        ') c_species = c_br
WHERE (advt == 'BrO       ') c_species = c_bro
WHERE (advt == 'Brx       ') c_species = c_bro
WHERE (advt == 'HOBr      ') c_species = c_hobr
WHERE (advt == 'BrONO2    ') c_species = c_brono2
WHERE (advt == 'BrCl      ') c_species = c_brcl
WHERE (advt == 'MeBr      ') c_species = c_mebr
WHERE (advt == 'HBr       ') c_species = c_hbr
WHERE (advt == 'CF2ClCFCl2') c_species = c_cf2clcfcl2
WHERE (advt == 'CHF2Cl    ') c_species = c_chf2cl
WHERE (advt == 'MeCCl3    ') c_species = c_meccl3
WHERE (advt == 'CCl4      ') c_species = c_ccl4
WHERE (advt == 'MeCl      ') c_species = c_mecl
WHERE (advt == 'CF2ClBr   ') c_species = c_cf2clbr
WHERE (advt == 'CF3Br     ') c_species = c_cf3br
WHERE (advt == 'CH2Br2    ') c_species = c_ch2br2
WHERE (advt == 'C5H8      ') c_species = c_c5h8
WHERE (advt == 'ISO2      ') c_species = c_iso2
WHERE (advt == 'ISOOH     ') c_species = c_isooh
WHERE (advt == 'ISON      ') c_species = c_ison
WHERE (advt == 'MACR      ') c_species = c_macr
WHERE (advt == 'MACRO2    ') c_species = c_macro2
WHERE (advt == 'MACROOH   ') c_species = c_macrooh
WHERE (advt == 'MPAN      ') c_species = c_mpan
WHERE (advt == 'HACET     ') c_species = c_hacet
WHERE (advt == 'MGLY      ') c_species = c_mgly
WHERE (advt == 'NALD      ') c_species = c_nald
WHERE (advt == 'HCOOH     ') c_species = c_hcooh
WHERE (advt == 'MeCO3H    ') c_species = c_meco3h
WHERE (advt == 'MeCO2H    ') c_species = c_meco2h
WHERE (advt == 'SO2       ') c_species = c_so2
WHERE (advt == 'SO3       ') c_species = c_so3
WHERE (advt == 'DMS       ') c_species = c_dms
WHERE (advt == 'DMSO      ') c_species = c_dmso
WHERE (advt == 'Me2S      ') c_species = c_me2s
WHERE (advt == 'COS       ') c_species = c_cos
WHERE (advt == 'H2S       ') c_species = c_h2s
WHERE (advt == 'CS2       ') c_species = c_cs2
WHERE (advt == 'MSA       ') c_species = c_msa
WHERE (advt == 'H2SO4     ') c_species = c_h2so4
WHERE (advt == 'NH3       ') c_species = c_nh3
WHERE (advt == 'MeOH      ') c_species = c_meoh
WHERE (advt == 'Monoterp  ') c_species = c_monoterp
WHERE (advt == 'Sec_Org   ') c_species = c_sec_org
WHERE (advt == 'SEC_ORG_I ') c_species = c_sec_org_i
WHERE (advt == 'O(3P)     ') c_species = c_o3p
WHERE (advt == 'OH        ') c_species = c_oh
WHERE (advt == 'HO2       ') c_species = c_ho2
WHERE (advt == 'MeOO      ') c_species = c_meoo
WHERE (advt == 'EtOO      ') c_species = c_etoo
WHERE (advt == 'MeCO3     ') c_species = c_meco3
WHERE (advt == 'n-PrOO    ') c_species = c_proo
WHERE (advt == 'i-PrOO    ') c_species = c_proo
WHERE (advt == 'EtCO3     ') c_species = c_etco3
WHERE (advt == 'MeCOCH2OO ') c_species = c_mecoch2oo
WHERE (advt == 'O(3P)S    ') c_species = c_o3p
WHERE (advt == 'O(1D)S    ') c_species = c_o1d
WHERE (advt == 'OHS       ') c_species = c_oh
WHERE (advt == 'HO2S      ') c_species = c_ho2
WHERE (advt == 'N         ') c_species = c_n
WHERE (advt == 'H         ') c_species = c_h
WHERE (advt == 'H2        ') c_species = c_h2       !for RAQ chemistry
WHERE (advt == 'C4H10     ') c_species = c_c4h10
WHERE (advt == 'MEK       ') c_species = c_mek
WHERE (advt == 'C2H4      ') c_species = c_c2h4
WHERE (advt == 'C3H6      ') c_species = c_c3h6
WHERE (advt == 'oXYLENE   ') c_species = c_oxylene
WHERE (advt == 's-BuOOH   ') c_species = c_buooh
WHERE (advt == 'CH3OH     ') c_species = c_ch3oh
WHERE (advt == 'MeOH      ') c_species = c_ch3oh
WHERE (advt == 'GLY       ') c_species = c_gly
WHERE (advt == 'MEMALD    ') c_species = c_memald
WHERE (advt == 'MVK       ') c_species = c_mvk
WHERE (advt == 'MVKOOH    ') c_species = c_mvkooh
WHERE (advt == 'TOLUENE   ') c_species = c_toluene
WHERE (advt == 'RNC2H4    ') c_species = c_rnc2h4
WHERE (advt == 'RNC3H6    ') c_species = c_rnc3h6
WHERE (advt == 'ORGNIT    ') c_species = c_orgnit
WHERE (advt == 'NOA       ') c_species = c_noa      ! for CRI chemistry
WHERE (advt == 'MeO2NO2   ') c_species = c_meo2no2
WHERE (advt == 'EtONO2    ') c_species = c_etono2
WHERE (advt == 'i-PrONO2  ') c_species = c_prono2
WHERE (advt == 'C2H2      ') c_species = c_c2h2
WHERE (advt == 'BENZENE   ') c_species = c_benzene
WHERE (advt == 'APINENE   ') c_species = c_monoterp !A/B-pin both monoterpenes
WHERE (advt == 'BPINENE   ') c_species = c_monoterp
WHERE (advt == 'TBUT2ENE  ') c_species = c_tbut2ene
WHERE (advt == 'n-PrOH    ') c_species = c_proh
WHERE (advt == 'i-PrOH    ') c_species = c_proh
WHERE (advt == 'EtOH      ') c_species = c_etoh
WHERE (advt == 'EtCO3H    ') c_species = c_etco3h
WHERE (advt == 'HOCH2CO3  ') c_species = c_hoch2co3
WHERE (advt == 'HOCH2CO3H ') c_species = c_hoch2co3h
WHERE (advt == 'HOCH2CH2O2') c_species = c_hoch2ch2o2
WHERE (advt == 'HOCH2CHO  ') c_species = c_hoch2cho
WHERE (advt == 'HOC2H4OOH ') c_species = c_hoc2h4ooh
WHERE (advt == 'HOC2H4NO3 ') c_species = c_hoc2h4no3
WHERE (advt == 'PHAN      ') c_species = c_phan
WHERE (advt == 'MeSCH2OO ') c_species = c_ch3sch2oo
WHERE (advt == 'MeS      ') c_species = c_ch3s
WHERE (advt == 'MeSO     ') c_species = c_ch3so
WHERE (advt == 'MeSO2    ') c_species = c_ch3so2
WHERE (advt == 'MeSO3    ') c_species = c_ch3so3
WHERE (advt == 'MSIA      ') c_species = c_msia
! CRI intermediates with direct analogues to other species have
! exact mass defined where relevant for dry deposition rates
WHERE (advt == 'CARB7     ') c_species = c_hacet
WHERE (advt == 'CARB3     ') c_species = c_gly
WHERE (advt == 'CARB6     ') c_species = c_mgly
WHERE (advt == 'UCARB10   ') c_species = c_macr
WHERE (advt == 'RN10NO3   ') c_species = c_prono2
WHERE (advt == 'RN13NO3   ') c_species = c_rn13no3
WHERE (advt == 'RN16NO3   ') c_species = c_rn16no3
WHERE (advt == 'RN19NO3   ') c_species = c_rn19no3
WHERE (advt == 'RA13NO3   ') c_species = c_ra13no3
WHERE (advt == 'RA16NO3   ') c_species = c_ra16no3
WHERE (advt == 'RA19NO3   ') c_species = c_ra19no3
WHERE (advt == 'RTX24NO3  ') c_species = c_rtx24no3
WHERE (advt == 'RN9NO3    ') c_species = c_rn9no3
WHERE (advt == 'RN12NO3   ') c_species = c_rn12no3
WHERE (advt == 'RN15NO3   ') c_species = c_rn15no3
WHERE (advt == 'RN18NO3   ') c_species = c_rn18no3
WHERE (advt == 'RU14NO3   ') c_species = c_ru14no3
WHERE (advt == 'RTN28NO3  ') c_species = c_rtn28no3
WHERE (advt == 'RTN25NO3  ') c_species = c_rtn25no3
WHERE (advt == 'RTN23NO3  ') c_species = c_rtn23no3
WHERE (advt == 'RTX28NO3  ') c_species = c_rtx28no3
WHERE (advt == 'RTX22NO3  ') c_species = c_rtx22no3
WHERE (advt == 'RN10OOH   ') c_species = c_prooh
WHERE (advt == 'RN13OOH   ') c_species = c_buooh
WHERE (advt == 'RN16OOH   ') c_species = c_rn16ooh
WHERE (advt == 'RN19OOH   ') c_species = c_rn19ooh
WHERE (advt == 'RN8OOH    ') c_species = c_alkaooh
WHERE (advt == 'RN11OOH   ') c_species = c_mekooh
WHERE (advt == 'RN14OOH   ') c_species = c_rn14ooh
WHERE (advt == 'RN17OOH   ') c_species = c_rn17ooh
WHERE (advt == 'RU14OOH   ') c_species = c_isooh
WHERE (advt == 'RU12OOH   ') c_species = c_cri ! actually has mass = 150.
WHERE (advt == 'RU10OOH   ') c_species = c_macrooh
WHERE (advt == 'NRU14OOH  ') c_species = c_nru14ooh
WHERE (advt == 'NRU12OOH  ') c_species = c_nru12ooh
WHERE (advt == 'RN9OOH    ') c_species = c_rn9ooh
WHERE (advt == 'RN12OOH   ') c_species = c_rn12ooh
WHERE (advt == 'RN15OOH   ') c_species = c_rn15ooh
WHERE (advt == 'RN18OOH   ') c_species = c_rn18ooh
WHERE (advt == 'NRN6OOH   ') c_species = c_nrn6ooh
WHERE (advt == 'NRN9OOH   ') c_species = c_nrn9ooh
WHERE (advt == 'NRN12OOH  ') c_species = c_nrn12ooh
WHERE (advt == 'RA13OOH   ') c_species = c_ra13ooh
WHERE (advt == 'RA16OOH   ') c_species = c_ra16ooh
WHERE (advt == 'RA19OOH   ') c_species = c_ra19ooh
WHERE (advt == 'RTN28OOH  ') c_species = c_rtn28ooh
WHERE (advt == 'NRTN28OOH ') c_species = c_nrtn28ooh
WHERE (advt == 'RTN26OOH  ') c_species = c_rtn26ooh
WHERE (advt == 'RTN25OOH  ') c_species = c_rtn25ooh
WHERE (advt == 'RTN24OOH  ') c_species = c_rtn24ooh
WHERE (advt == 'RTN23OOH  ') c_species = c_rtn23ooh
WHERE (advt == 'RTN14OOH  ') c_species = c_rtn14ooh
WHERE (advt == 'RTN10OOH  ') c_species = c_rtn10ooh
WHERE (advt == 'RTX28OOH  ') c_species = c_rtx28ooh
WHERE (advt == 'RTX24OOH  ') c_species = c_rtx24ooh
WHERE (advt == 'RTX22OOH  ') c_species = c_rtx22ooh
WHERE (advt == 'NRTX28OOH ') c_species = c_nrtx28ooh
WHERE (advt == 'CARB14    ') c_species = c_carb14
WHERE (advt == 'CARB17    ') c_species = c_carb17
WHERE (advt == 'CARB11A   ') c_species = c_carb11a
WHERE (advt == 'CARB10    ') c_species = c_carb10
WHERE (advt == 'CARB13    ') c_species = c_carb13
WHERE (advt == 'CARB16    ') c_species = c_carb16
WHERE (advt == 'CARB9     ') c_species = c_carb9
WHERE (advt == 'CARB12    ') c_species = c_carb12
WHERE (advt == 'CARB15    ') c_species = c_carb15
WHERE (advt == 'UCARB12   ') c_species = c_ucarb12
WHERE (advt == 'NUCARB12  ') c_species = c_nucarb12
WHERE (advt == 'UDCARB8   ') c_species = c_udcarb8
WHERE (advt == 'UDCARB11  ') c_species = c_udcarb11
WHERE (advt == 'UDCARB14  ') c_species = c_udcarb14
WHERE (advt == 'TNCARB26  ') c_species = c_tncarb26
WHERE (advt == 'TNCARB10  ') c_species = c_tncarb10
WHERE (advt == 'TNCARB12  ') c_species = c_tncarb12
WHERE (advt == 'TNCARB11  ') c_species = c_tncarb11
WHERE (advt == 'CCARB12   ') c_species = c_ccarb12
WHERE (advt == 'TNCARB15  ') c_species = c_tncarb15
WHERE (advt == 'RCOOH25   ') c_species = c_rcooh25
WHERE (advt == 'TXCARB24  ') c_species = c_txcarb24
WHERE (advt == 'TXCARB22  ') c_species = c_txcarb22
WHERE (advt == 'RU12PAN   ') c_species = c_ru12pan
WHERE (advt == 'RTN26PAN  ') c_species = c_rtn26pan
WHERE (advt == 'AROH14    ') c_species = c_aroh14
WHERE (advt == 'AROH17    ') c_species = c_aroh17
WHERE (advt == 'ARNOH14   ') c_species = c_arnoh14
WHERE (advt == 'ARNOH17   ') c_species = c_arnoh17
WHERE (advt == 'ANHY      ') c_species = c_anhy
WHERE (advt == 'DHPCARB9  ') c_species = c_dhpcarb9
WHERE (advt == 'HPUCARB12 ') c_species = c_hpucarb12
WHERE (advt == 'HUCARB9   ') c_species = c_hucarb9
WHERE (advt == 'IEPOX     ') c_species = c_iepox
WHERE (advt == 'HMML      ') c_species = c_hmml
WHERE (advt == 'DHPR12OOH ') c_species = c_dhpr12ooh
WHERE (advt == 'DHCARB9   ') c_species = c_dhcarb9
WHERE (advt == 'RU12NO3   ') c_species = c_ru12no3
WHERE (advt == 'RU10NO3   ') c_species = c_ru10no3

! Mass for these species does not matter as not emitted or deposited
! Use default mass 150g/mole
WHERE (advt == 'RAROH14   ') c_species = c_cri
WHERE (advt == 'RAROH17   ') c_species = c_cri
WHERE (advt == 'PASSIVE O3') c_species = 1.0
WHERE (advt == 'AGE OF AIR') c_species = 1.0

! non-advected tracers
c_na_species=0.0
WHERE (nadvt == 'OH        ') c_na_species = c_oh
WHERE (nadvt == 'HO2       ') c_na_species = c_ho2
WHERE (nadvt == 'O(3P)     ') c_na_species = c_o3p
WHERE (nadvt == 'O3P       ') c_na_species = c_o3p
WHERE (nadvt == 'O(1D)     ') c_na_species = c_o1d
WHERE (nadvt == 'O1D       ') c_na_species = c_o1d
WHERE (nadvt == 'MeOO      ') c_na_species = c_meoo
WHERE (nadvt == 'EtOO      ') c_na_species = c_etoo
WHERE (nadvt == 'MeCO3     ') c_na_species = c_meco3
WHERE (nadvt == 'n-PrOO    ') c_na_species = c_proo
WHERE (nadvt == 'i-PrOO    ') c_na_species = c_proo
WHERE (nadvt == 'EtCO3     ') c_na_species = c_etco3
WHERE (nadvt == 'MeCOCH2OO ') c_na_species = c_mecoch2oo
WHERE (nadvt == 'O(3P)S    ') c_na_species = c_o3p
WHERE (nadvt == 'O(1D)S    ') c_na_species = c_o1d
WHERE (nadvt == 'OHS       ') c_na_species = c_oh
WHERE (nadvt == 'HO2S      ') c_na_species = c_ho2
WHERE (nadvt == 's-BuOO    ') c_na_species = c_buoo     !for RAQ chemistry
WHERE (nadvt == 'MEKO2     ') c_na_species = c_meko2
WHERE (nadvt == 'HOC2H4O2  ') c_na_species = c_hoc2h4o2
WHERE (nadvt == 'HOC3H6O2  ') c_na_species = c_hoc3h6o2
WHERE (nadvt == 'OXYL1     ') c_na_species = c_oxyl1
WHERE (nadvt == 'MEMALD1   ') c_na_species = c_memald1
WHERE (nadvt == 'HOIPO2    ') c_na_species = c_hoipo2
WHERE (nadvt == 'HOMVKO2   ') c_na_species = c_homvko2
WHERE (nadvt == 'TOLP1     ') c_na_species = c_tolp1
WHERE (nadvt == 'ISO2      ') c_na_species = c_iso2
WHERE (nadvt == 'MACRO2    ') c_na_species = c_macro2
! Non-advected peroxy radicals for CRI chemistry:
WHERE (nadvt == 'HOCH2CO3  ') c_na_species = c_hoch2co3
WHERE (nadvt == 'HOCH2CH2O2') c_na_species = c_hoch2ch2o2
! Mass for these species does not matter as not emitted, transported
! or deposited
WHERE (nadvt == 'RN10O2    ') c_na_species = c_cri
WHERE (nadvt == 'RN13O2    ') c_na_species = c_cri
WHERE (nadvt == 'RN16O2    ') c_na_species = c_cri
WHERE (nadvt == 'RN19O2    ') c_na_species = c_cri
WHERE (nadvt == 'RN13AO2   ') c_na_species = c_cri
WHERE (nadvt == 'RN16AO2   ') c_na_species = c_cri
WHERE (nadvt == 'RA13O2    ') c_na_species = c_cri
WHERE (nadvt == 'RA16O2    ') c_na_species = c_cri
WHERE (nadvt == 'RA19AO2   ') c_na_species = c_cri
WHERE (nadvt == 'RA19CO2   ') c_na_species = c_cri
WHERE (nadvt == 'RN9O2     ') c_na_species = c_cri
WHERE (nadvt == 'RN12O2    ') c_na_species = c_cri
WHERE (nadvt == 'RN15O2    ') c_na_species = c_cri
WHERE (nadvt == 'RN18O2    ') c_na_species = c_cri
WHERE (nadvt == 'RN15AO2   ') c_na_species = c_cri
WHERE (nadvt == 'RN18AO2   ') c_na_species = c_cri
WHERE (nadvt == 'RN8O2     ') c_na_species = c_cri
WHERE (nadvt == 'RN11O2    ') c_na_species = c_cri
WHERE (nadvt == 'RN14O2    ') c_na_species = c_cri
WHERE (nadvt == 'RN17O2    ') c_na_species = c_cri
WHERE (nadvt == 'RU14O2    ') c_na_species = c_cri
WHERE (nadvt == 'RU12O2    ') c_na_species = c_cri
WHERE (nadvt == 'RU10O2    ') c_na_species = c_cri
WHERE (nadvt == 'NRN6O2    ') c_na_species = c_cri
WHERE (nadvt == 'NRN9O2    ') c_na_species = c_cri
WHERE (nadvt == 'NRN12O2   ') c_na_species = c_cri
WHERE (nadvt == 'NRU14O2   ') c_na_species = c_cri
WHERE (nadvt == 'NRU12O2   ') c_na_species = c_cri
WHERE (nadvt == 'RTN28O2   ') c_na_species = c_cri
WHERE (nadvt == 'NRTN28O2  ') c_na_species = c_cri
WHERE (nadvt == 'RTN26O2   ') c_na_species = c_cri
WHERE (nadvt == 'RTN25O2   ') c_na_species = c_cri
WHERE (nadvt == 'RTN24O2   ') c_na_species = c_cri
WHERE (nadvt == 'RTN23O2   ') c_na_species = c_cri
WHERE (nadvt == 'RTN14O2   ') c_na_species = c_cri
WHERE (nadvt == 'RTN10O2   ') c_na_species = c_cri
WHERE (nadvt == 'RTX28O2   ') c_na_species = c_cri
WHERE (nadvt == 'NRTX28O2  ') c_na_species = c_cri
WHERE (nadvt == 'RTX24O2   ') c_na_species = c_cri
WHERE (nadvt == 'RTX22O2   ') c_na_species = c_cri
WHERE (nadvt == 'RO2       ') c_na_species = c_cri
WHERE (nadvt == 'DHPR12O2  ') c_na_species = c_dhpr12o2
WHERE (nadvt == 'RU10AO2   ') c_na_species = c_ru10ao2
WHERE (nadvt == 'MACO3     ') c_na_species = c_maco3


!     Initialise tracer numbers

n_h2o    = 0
n_cl     = 0
n_sx     = 0
n_n2o5   = 0
n_ox     = 0
n_o3     = 0
n_o3s    = 0
n_nox    = 0
n_no     = 0
n_no2    = 0
n_no3    = 0
n_ho2no2 = 0
n_hono2  = 0
n_pan    = 0
n_ison   = 0
n_orgnit = 0
n_rnc2h4 = 0
n_rnc3h6 = 0
n_h2o2   = 0
n_ch4    = 0
n_h2     = 0
n_dms    = 0
n_so2    = 0
n_so3    = 0
n_h2so4  = 0
n_sec_org= 0
n_sec_org_i = 0
n_passive = 0
n_hcho   = 0
n_c2h6   = 0
n_mecho  = 0
n_c3h8   = 0
n_me2co  = 0
n_c5h8   = 0
n_mgly   = 0
n_c4h10  = 0
n_mek    = 0
n_c2h4   = 0
n_c3h6   = 0
n_oxylene = 0
n_ch3oh  = 0
n_gly    = 0
n_mvk    = 0
n_toluene = 0
n_co = 0
n_n2o = 0
n_cf2cl2 = 0  ! CFC-12
n_cfcl3 = 0   ! CFC-11
n_mebr = 0    ! CH3Br
n_nh3 = 0

!     Find tracer numbers

DO m=1,jpctr
  SELECT CASE (advt(m))
  CASE ('Ox        ')
    n_ox     = m
    n_ox1    = m
  CASE ('O3        ')
    n_o3     = m
  CASE ('O3S       ')
    n_o3s    = m
  CASE ('Sx        ')   ! stratospheric ozone tracer
    n_sx     = m
    n_sx1    = m
  CASE ('NOx       ')
    n_nox    = m
  CASE ('NO        ')
    n_no     = m
  CASE ('NO2       ')
    n_no2    = m
  CASE ('NO3       ')
    n_no3    = m
  CASE ('N2O5      ')
    n_n2o5   = m
  CASE ('HO2NO2    ')
    n_ho2no2 = m
  CASE ('HONO2     ')
    n_hono2  = m
  CASE ('PAN       ')
    n_pan    = m
  CASE ('ISON      ')
    n_ison   = m
  CASE ('ORGNIT    ')
    n_orgnit = m
  CASE ('RNC2H4    ')
    n_rnc2h4 = m
  CASE ('RNC3H6    ')
    n_rnc3h6 = m
  CASE ('H2O2      ')
    n_h2o2   = m
  CASE ('CH4       ')
    n_ch4    = m
  CASE ('H2O       ')
    n_h2o    = m
  CASE ('Cl        ')
    n_cl     = m
  CASE ('H2        ')
    n_h2     = m
  CASE ('DMS       ')
    n_dms    = m
  CASE ('SO2       ')
    n_so2    = m
  CASE ('SO3       ')
    n_so3    = m
  CASE ('H2SO4     ')
    n_h2so4  = m
  CASE ('Sec_Org   ')
    n_sec_org= m
  CASE ('SEC_ORG_I ')
    n_sec_org_i = m
  CASE ('BrO       ')
    n_bro    = m
  CASE ('HCl       ')
    n_hcl    = m
  CASE ('CFCl3     ')   ! CFC-11
    n_cfcl3  = m
  CASE ('CF2Cl2    ')   ! CFC-12
    n_cf2cl2 = m
  CASE ('O(1D)     ')
    n_o1d    = m
  CASE ('HCHO      ')
    n_hcho   = m
  CASE ('C2H6      ')
    n_c2h6   = m
  CASE ('MeCHO     ')
    n_mecho  = m
  CASE ('C3H8      ')
    n_c3h8   = m
  CASE ('Me2CO     ')
    n_me2co  = m
  CASE ('C5H8      ')
    n_c5h8   = m
  CASE ('MGLY      ')
    n_mgly   = m
  CASE ('C4H10     ')
    n_c4h10  = m
  CASE ('MEK       ')
    n_mek    = m
  CASE ('C2H4      ')
    n_c2h4   = m
  CASE ('C3H6      ')
    n_c3h6   = m
  CASE ('oXYLENE   ')
    n_oxylene = m
  CASE ('CH3OH     ')
    n_ch3oh  = m
  CASE ('GLY       ')
    n_gly    = m
  CASE ('MVK       ')
    n_mvk    = m
  CASE ('TOLUENE   ')
    n_toluene = m
  CASE ('CO        ')
    n_co = m
  CASE ('MeBr      ')   ! CH3Br
    n_mebr = m
  CASE ('N2O       ')
    n_n2o = m
  CASE ('NH3       ')
    n_nh3 = m
  END SELECT
END DO
DO i=1,jpctr
  IF ((advt(i) == 'Ox        ') .OR. (advt(i) == 'O3        '))                &
    n_o3 = i
END DO

IF (n_o3 == 0 .AND. .NOT. (ukca_config%l_ukca_offline .OR.                     &
                           ukca_config%l_ukca_offline_be)) THEN
  cmessage='Ozone not found among chemical tracers.'
  errcode=1
  CALL ereport('UKCA_CSPECIES',errcode,cmessage)
END IF

!     Initialise species numbers

nn_o3        = 0
nn_o3s       = 0
nn_oh        = 0
nn_ho2       = 0
nn_h2o2      = 0
nn_no        = 0
nn_no2       = 0
nn_o1d       = 0
nn_o3p       = 0
nn_meoo      = 0
nn_meco3     = 0
nn_etoo      = 0
nn_etco3     = 0
nn_nproo     = 0
nn_nprooh    = 0
nn_iproo     = 0
nn_iprooh    = 0
nn_mecoch2oo = 0
nn_so2       = 0
nn_so3       = 0
nn_h2so4     = 0
nn_ohs       = 0
nn_ho2s      = 0
nn_o1ds      = 0
nn_o3ps      = 0
nn_n         = 0
nn_meo       = 0
nn_cl        = 0
nn_clo       = 0
nn_cl2o2     = 0
nn_hcl       = 0
nn_br        = 0
nn_bro       = 0
nn_buoo      = 0   ! species for RAQ chemistry
nn_meko2     = 0
nn_hoc2h4o2  = 0
nn_hoc3h6o2  = 0
nn_oxyl1     = 0
nn_memald1   = 0
nn_hoipo2    = 0
nn_homvko2   = 0
nn_tolp1     = 0
nn_nh3 = 0

!     Find species numbers

DO m=1,jpspec
  SELECT CASE(speci(m))
  CASE ('O3        ')
    nn_o3 = m
  CASE ('O3S       ')
    nn_o3s = m
  CASE ('OH        ')
    nn_oh = m
  CASE ('HO2       ')
    nn_ho2 = m
  CASE ('H2O2      ')
    nn_h2o2 = m
  CASE ('NO        ')
    nn_no = m
  CASE ('NO2       ')
    nn_no2 = m
  CASE ('O(1D)     ')
    nn_o1d = m
  CASE ('O(3P)     ')
    nn_o3p = m
  CASE ('MeOO      ')
    nn_meoo = m
  CASE ('MeCO3     ')
    nn_meco3 = m
  CASE ('EtOO      ')
    nn_etoo = m
  CASE ('EtCO3     ')
    nn_etco3 = m
  CASE ('n-PrOO    ')
    nn_nProo = m
  CASE ('n-PrOOH   ')
    nn_nprooh = m
  CASE ('i-PrOO    ')
    nn_iProo = m
  CASE ('i-PrOOH   ')
    nn_iprooh = m
  CASE ('MeCOCH2OO ')
    nn_mecoch2oo = m
  CASE ('NO3       ')
    nn_no3 = m
  CASE ('N2O5      ')
    nn_n2o5 = m
  CASE ('HO2NO2    ')
    nn_ho2no2 = m
  CASE ('HONO2     ')
    nn_hono2 = m
  CASE ('CH4       ')
    nn_ch4 = m
  CASE ('SO2       ')
    nn_so2 = m
  CASE ('SO3       ')
    nn_so3 = m
  CASE ('H2SO4     ')
    nn_h2so4 = m
  CASE ('O(1D)S    ')
    nn_o1ds = m
  CASE ('O(3P)S    ')
    nn_o3ps = m
  CASE ('OHS       ')
    nn_ohs = m
  CASE ('HO2S      ')
    nn_ho2s = m
  CASE ('N         ')
    nn_n = m
  CASE ('MeO       ')
    nn_meo = m
  CASE ('Cl        ')
    nn_cl = m
  CASE ('ClO       ')
    nn_clo = m
  CASE ('Cl2O2     ')
    nn_cl2o2 = m
  CASE ('HCl       ')
    nn_hcl = m
  CASE ('Br        ')
    nn_br = m
  CASE ('BrO       ')
    nn_bro = m
  CASE ('s-BuOO    ')
    nn_buoo = m
  CASE ('MEKO2     ')
    nn_meko2 = m
  CASE ('HOC2H4O2  ')
    nn_hoc2h4o2 = m
  CASE ('HOC3H6O2  ')
    nn_hoc3h6o2 = m
  CASE ('OXYL1     ')
    nn_oxyl1 = m
  CASE ('MEMALD1   ')
    nn_memald1 = m
  CASE ('HOIPO2    ')
    nn_hoipo2 = m
  CASE ('HOMVKO2   ')
    nn_homvko2 = m
  CASE ('TOLP1     ')
    nn_tolp1 = m
  CASE ('CFCl3     ')
    nn_cfcl3 = m
  CASE ('CF2Cl2    ')
    nn_cf2cl2 = m
  CASE ('NH3       ')
    nn_nh3 = m
  END SELECT
END DO
IF (printstatus > prstatus_oper) THEN
  CALL umPrint( 'Species indices:',src='ukca_cspecies')
  WRITE(umMessage,'(A,1X,I6,A,1X,I6)') 'O3 = ',nn_o3,' OH = ',nn_oh
  CALL umPrint(umMessage,src='ukca_cspecies')
  WRITE(umMessage,'(A,1X,I6,A,1X,I6)') 'HO2 = ',nn_ho2,' NO = ',nn_no
  CALL umPrint(umMessage,src='ukca_cspecies')
  WRITE(umMessage,'(A,1X,I6,A,1X,I6)') 'NO2 = ',nn_no2,' O(1D) = ',nn_o1d
  CALL umPrint(umMessage,src='ukca_cspecies')
  WRITE(umMessage,'(A,1X,I6,A,1X,I6)') 'MeOO = ',nn_meoo,' MeCO3 = ',nn_meco3
  CALL umPrint(umMessage,src='ukca_cspecies')
  WRITE(umMessage,'(A,1X,I6,A,1X,I6)') 'EtOO = ',nn_etoo,' EtCO3 = ',nn_etco3
  CALL umPrint(umMessage,src='ukca_cspecies')
  WRITE(umMessage,'(A,1X,I6,A,1X,I6)') 'n-Proo = ',                            &
                                    nn_nproo,'n-Prooh = ',nn_nprooh
  CALL umPrint(umMessage,src='ukca_cspecies')
  WRITE(umMessage,'(A,1X,I6,A,1X,I6)') 'i-Proo = ',                            &
                                    nn_iproo,'i-Prooh = ',nn_iprooh
  CALL umPrint(umMessage,src='ukca_cspecies')
  WRITE(umMessage,'(A,1X,I6,A,1X,I6)') 'MeCOCH2OO = ',nn_mecoch2oo
  CALL umPrint(umMessage,src='ukca_cspecies')
  WRITE(umMessage,'(A,1X,I6,A,1X,I6)') 'SO2: ',nn_so2,' H2O2: ',nn_h2o2
  CALL umPrint(umMessage,src='ukca_cspecies')
  WRITE(umMessage,'(A,1X,I6,A,1X,I6)') 'SO3: ',nn_so3,' H2SO4: ',nn_h2so4
  CALL umPrint(umMessage,src='ukca_cspecies')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calc_cspecies

END MODULE ukca_cspecies
