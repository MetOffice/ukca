! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module to define setup and indices for gas and aerosol tracers.
!    Many of these are not used in UM but kept for consistency with
!    set-up in TOMCAT.
!    Contains public subroutines:
!      UKCA_INDICES_ORGV1_SOto3
!      UKCA_INDICES_ORGV1_SOto3_SOLINSOL
!      UKCA_INDICES_ORGV1_SOto3_COUPL
!      UKCA_INDICES_ORGV1_SOto6
!      UKCA_INDICES_ORGV1_SOto6_COUPL
!      UKCA_INDICES_SV1
!      UKCA_INDICES_SV1_COUPLED
!      UKCA_INDICES_NOCHEM
!      UKCA_INDICES_TRAQU38
!      UKCA_INDICES_TRAQU9
!      UKCA_INDICES_SUSS_4MODE
!      UKCA_INDICES_SUSSBCOC_5MODE
!      UKCA_INDICES_SUSSBCOC_4MODE
!      UKCA_INDICES_SUSSBCOCSO_5MODE
!      UKCA_INDICES_SUSSBCOCSO_4MODE
!      UKCA_INDICES_DUonly_2MODE
!      UKCA_INDICES_DUonly_3MODE (needs to be added at some point)
!      UKCA_INDICES_SUSSBCOCDU_7MODE
!      UKCA_INDICES_SUSSBCOCDU_4MODE
!      UKCA_INDICES_SUSSBCOCNTNH_5MODE
!      UKCA_INDICES_SOLINSOL_6MODE
!    which define particular setups.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
! Subroutine Interface:
MODULE ukca_setup_indices
! ---------------------------------------------------------------------|
!  Module to define setup and indices for gas and aerosol tracers.
!  Many of these are not used in UM but kept for consistency with
!  set-up in TOMCAT.



USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER :: ntraer     ! # of aerosol advected tracers
INTEGER :: nchemg     ! # of gas phase chemistry tracers
INTEGER :: ichem      ! 1/0 = gas phase chemistry tracers on/off
INTEGER :: nadvg      ! # of gas phase advected tracers
INTEGER :: noffox     ! # of offline oxidant species
INTEGER :: ntrag      ! total # of gas phase species
INTEGER :: budget     ! 1/0 = budget fields on/off
INTEGER :: nbudchem   ! # of gas chem budget fields
INTEGER :: idustdep   ! 1/0 = size-resolved dust dep. fields on/off
INTEGER :: ndustdep   ! # of size-resolved dust deposition fields
INTEGER :: nbudaer    ! # of aerosol budget fields
INTEGER :: nbudaertot ! total # of aersol budget fields
INTEGER :: nbudget    ! total # of budget terms
INTEGER :: traqu      ! 1/0 = extra diagnostic terms on/off
INTEGER :: ntraqu     ! # of extra diagnostic terms
INTEGER :: gasbudget  ! 1/0 = gas phase chemistry fluxes on/off
INTEGER :: ngasbudget ! # of gas phase chemistry fluxes
!
! .. (indices for 40 gas phase tropospheric chemistry species in S0)
INTEGER :: mox
INTEGER :: mnox
INTEGER :: mn2o5
INTEGER :: mhno4
INTEGER :: mhno3
INTEGER :: mh2o2
INTEGER :: mch4
INTEGER :: mco
INTEGER :: mch2o
INTEGER :: mmhp
INTEGER :: mhono
INTEGER :: mc2h6
INTEGER :: metooh
INTEGER :: mmecho
INTEGER :: mpan
INTEGER :: mc3h8
INTEGER :: mpnooh
INTEGER :: mpiooh
INTEGER :: metcho
INTEGER :: mme2co
INTEGER :: mmecoh
INTEGER :: mppan
INTEGER :: mmeno3
INTEGER :: moxs
INTEGER :: mnoys
INTEGER :: misop
INTEGER :: mc2h4
INTEGER :: mc2h2
INTEGER :: misooh
INTEGER :: mison
INTEGER :: mmacr
INTEGER :: mmacrooh
INTEGER :: mmpan
INTEGER :: mhacet
INTEGER :: mmgly
INTEGER :: mnald
INTEGER :: mhcooh
INTEGER :: mmeco3h
INTEGER :: MMeCO2H
INTEGER :: mmeoh
!
! .. (indices for advected gas phase tracers in array S0)
! .. (8 sulfur species, H2O2, MONOTER, Sec_Org, Q3D, PT)
INTEGER :: msotwo     ! for SO2
INTEGER :: MMeSMe     ! for DMS
INTEGER :: mh2so4     ! for H2SO4
INTEGER :: mdmso      ! for DMSO
INTEGER :: mmsa       ! for MSA
INTEGER :: mcs2       ! for CS2
INTEGER :: mh2s       ! for H2S
INTEGER :: mcos       ! for COS
INTEGER :: mmonoter   ! for lumped monoterpene species
INTEGER :: msec_org   ! for involatile organic species
INTEGER :: mh2o2f     ! for H2O2 semi-prognostic
INTEGER :: mq3d       ! for water vapour
INTEGER :: mpt        ! for potential temperature
!
! .. (indices for tropospheric chemistry species in array ST)
INTEGER :: no
INTEGER :: no1d
INTEGER :: no3
INTEGER :: nno
INTEGER :: nno3
INTEGER :: nno2
INTEGER :: nn2o5
INTEGER :: nhno4
INTEGER :: nhno3
INTEGER :: noh
INTEGER :: nho2
INTEGER :: nh2o2
INTEGER :: nch4
INTEGER :: nco
INTEGER :: nch2o
INTEGER :: nmeoo
INTEGER :: nh2o
INTEGER :: nmhp
INTEGER :: nhono
INTEGER :: nc2h6
INTEGER :: netoo
INTEGER :: netooh
INTEGER :: nmecho
INTEGER :: nmeco3
INTEGER :: npan
INTEGER :: nc3h8
INTEGER :: npnoo
INTEGER :: npioo
INTEGER :: npnooh
INTEGER :: npiooh
INTEGER :: netcho
INTEGER :: netco3
INTEGER :: nme2co
INTEGER :: nmecoo
INTEGER :: nmecoh
INTEGER :: nppan
INTEGER :: nmeno3
INTEGER :: nos
INTEGER :: no1ds
INTEGER :: no3s
INTEGER :: nnoxs
INTEGER :: nhno3s
INTEGER :: nnoys
INTEGER :: nisop
INTEGER :: nc2h4
INTEGER :: nc2h2
INTEGER :: niso2
INTEGER :: nisooh
INTEGER :: nison
INTEGER :: nmacr
INTEGER :: nmacro2
INTEGER :: nmacrooh
INTEGER :: nmpan
INTEGER :: nhacet
INTEGER :: nmgly
INTEGER :: nnald
INTEGER :: nhcooh
INTEGER :: nmeco3h
INTEGER :: nmeco2h
INTEGER :: nmeoh
!
! .. (indices for all gas phase tracers in array ST)
INTEGER :: nsotwo     ! for SO2
INTEGER :: NMeSMe     ! for DMS
INTEGER :: nh2so4     ! for H2SO4
INTEGER :: ndmso      ! for DMSO
INTEGER :: nmsa       ! for MSA
INTEGER :: ncs2       ! for CS2
INTEGER :: nh2s       ! for H2S
INTEGER :: ncos       ! for COS
INTEGER :: nmonoter   ! for lumped monoterpene species
INTEGER :: nsec_org   ! for involatile organic species
INTEGER :: nh2o2f     ! for H2O2 (semi-prognostic)
INTEGER :: no3f       ! for O3 (offline)
INTEGER :: nohf       ! for OH (offline)
INTEGER :: nno3f      ! for NO3 (offline)
INTEGER :: nq3d       ! for water vapour
INTEGER :: npt        ! for potential temperature
!
! .. (indices for gas phase budget mass fluxes)
INTEGER :: ndmsemoc   ! for DMS emissions flux (oceanic sources)
INTEGER :: ndmstend   ! for DMS ASAD tendency (combined chem,ddep,wdep)
INTEGER :: nso2eman   ! for SO2 emissions flux (anthropogenic sources)
INTEGER :: nso2embm   ! for SO2 emissions flux (biomass burning sources)
INTEGER :: nso2emvl   ! for SO2 emissions flux (volcanic sources)
INTEGER :: nso2tend   ! for SO2 ASAD tendency (combined chem,ddep,wdep)
INTEGER :: nso2ddep   ! for SO2 dry deposition flux
INTEGER :: nso2wdep   ! for SO2 wet deposition flux
INTEGER :: nh2so4tend ! for H2SO4 ASAD tendency (combined chem,ddep,wdep)
INTEGER :: nh2so4ddep ! for H2SO4 dry deposition flux
INTEGER :: ncoseman   ! for COS emissions flux (anthropogenic sources)
INTEGER :: ncosemoc   ! for COS emissions flux (oceanic sources)
INTEGER :: ncostend   ! for COS ASAD tendency (combined chem,ddep,wdep)
INTEGER :: ncs2eman   ! for CS2 emissions flux (anthropogenic sources)
INTEGER :: ncs2emoc   ! for CS2 emissions flux (oceanic sources)
INTEGER :: ncs2tend   ! for CS2 ASAD tendency (combined chem,ddep,wdep)
INTEGER :: ndmsotend  ! for DMSO ASAD tendency (combined chem,ddep,wdep)
INTEGER :: ndmsoddep  ! for DMSO dry deposition flux
INTEGER :: nmsatend   ! for MSA ASAD tendency (combined chem,ddep,wdep)
INTEGER :: nmsaddep   ! for MSA dry deposition flux
INTEGER :: nterp_em   ! for terpene emissions flux (biogenic sources)
INTEGER :: nterp_tend ! for terpene ASAD tendency (combined chem,ddep,wdep)
INTEGER :: nterp_ddep ! for terpene dry deposition flux
INTEGER :: nsorg_tend ! for Sec_Org ASAD tendency (combined chem,ddep,wdep)
INTEGER :: nsorg_ddep ! for Sec_Org dry deposition flux
INTEGER :: nsorg_wdep ! for Sec_Org wet deposition flux
!
! .. (indices for gas phase chemistry fluxes)
INTEGER :: iohdms1  ! DMS +OH -->   SO2(+   MeOO +HCHO)
INTEGER :: iohdms2  ! DMS +OH -->0.6SO2 +0.4DMSO(+MeOO)
INTEGER :: ino3dms  ! DMS +NO3-->   SO2(+  HONO2 +MeOO+HCHO)
INTEGER :: idmsooh1 ! DMSO+OH -->   SO2
INTEGER :: idmsooh2 ! DMSO+OH -->0.6SO2 +0.4MSA
INTEGER :: ics2oh   ! DMS +OH -->   SO2 +COS
INTEGER :: ih2soh   ! DMS +OH -->   SO2
INTEGER :: icosoh   ! DMS +OH -->   SO2
!
! .. (indices for aerosol budget mass fluxes)
INTEGER :: nmasprimsuaitsol ! for primary SU ems (all srces) to Aitsol
INTEGER :: nmasprimsuaccsol ! for primary SU ems (all srces) to accsol
INTEGER :: nmasprimsucorsol ! for primary SU ems (all srces) to corsol
INTEGER :: nmasprimssaccsol ! for primary SS ems (ocean) to accsol
INTEGER :: nmasprimsscorsol ! for primary SS ems (ocean) to corsol
INTEGER :: nmasprimbcaitsol ! for primary BC ems (all srces) to Aitsol
INTEGER :: nmasprimbcaitins ! for primary BC ems (all srces) to Aitins
INTEGER :: nmasprimocaitsol ! for primary OC ems (all srces) to Aitsol
INTEGER :: nmasprimocaitins ! for primary OC ems (all srces) to Aitins
INTEGER :: nmasddepsunucsol ! for drydep of SU from nucsol
INTEGER :: nmasddepsuaitsol ! for drydep of SU from Aitsol
INTEGER :: nmasddepsuaccsol ! for drydep of SU from accsol
INTEGER :: nmasddepsucorsol ! for drydep of SU from corsol
INTEGER :: nmasddepssaccsol ! for drydep of SS from accsol
INTEGER :: nmasddepsscorsol ! for drydep of SS from corsol
INTEGER :: nmasddepbcaitsol ! for drydep of BC from Aitsol
INTEGER :: nmasddepbcaccsol ! for drydep of BC from accsol
INTEGER :: nmasddepbccorsol ! for drydep of BC from corsol
INTEGER :: nmasddepbcaitins ! for drydep of BC from Aitins
INTEGER :: nmasddepocnucsol ! for drydep of OC from nucsol
INTEGER :: nmasddepocaitsol ! for drydep of OC from Aitsol
INTEGER :: nmasddepocaccsol ! for drydep of OC from accsol
INTEGER :: nmasddepoccorsol ! for drydep of OC from corsol
INTEGER :: nmasddepocaitins ! for drydep of OC from Aitins
INTEGER :: nmasddepsonucsol ! for drydep of SO from nucsol
INTEGER :: nmasddepsoaitsol ! for drydep of SO from Aitsol
INTEGER :: nmasddepsoaccsol ! for drydep of SO from accsol
INTEGER :: nmasddepsocorsol ! for drydep of SO from corsol
INTEGER :: nmasnuscsunucsol ! for nuscav of SU from nucsol
INTEGER :: nmasnuscsuaitsol ! for nuscav of SU from Aitsol
INTEGER :: nmasnuscsuaccsol ! for nuscav of SU from accsol
INTEGER :: nmasnuscsucorsol ! for nuscav of SU from corsol
INTEGER :: nmasnuscssaccsol ! for nuscav of SS from accsol
INTEGER :: nmasnuscsscorsol ! for nuscav of SS from corsol
INTEGER :: nmasnuscbcaitsol ! for nuscav of BC from Aitsol
INTEGER :: nmasnuscbcaccsol ! for nuscav of BC from accsol
INTEGER :: nmasnuscbccorsol ! for nuscav of BC from corsol
INTEGER :: nmasnuscbcaitins ! for nuscav of BC from Aitins
INTEGER :: nmasnuscocnucsol ! for nuscav of OC from nucsol
INTEGER :: nmasnuscocaitsol ! for nuscav of OC from Aitsol
INTEGER :: nmasnuscocaccsol ! for nuscav of OC from accsol
INTEGER :: nmasnuscoccorsol ! for nuscav of OC from corsol
INTEGER :: nmasnuscocaitins ! for nuscav of OC from Aitins
INTEGER :: nmasnuscsonucsol ! for nuscav of SO from nucsol
INTEGER :: nmasnuscsoaitsol ! for nuscav of SO from Aitsol
INTEGER :: nmasnuscsoaccsol ! for nuscav of SO from accsol
INTEGER :: nmasnuscsocorsol ! for nuscav of SO from corsol
INTEGER :: nmasimscsunucsol ! for imscav of SU from nucsol
INTEGER :: nmasimscsuaitsol ! for imscav of SU from Aitsol
INTEGER :: nmasimscsuaccsol ! for imscav of SU from accsol
INTEGER :: nmasimscsucorsol ! for imscav of SU from corsol
INTEGER :: nmasimscssaccsol ! for imscav of SS from accsol
INTEGER :: nmasimscsscorsol ! for imscav of SS from corsol
INTEGER :: nmasimscbcaitsol ! for imscav of BC from Aitsol
INTEGER :: nmasimscbcaccsol ! for imscav of BC from accsol
INTEGER :: nmasimscbccorsol ! for imscav of BC from corsol
INTEGER :: nmasimscbcaitins ! for imscav of BC from Aitins
INTEGER :: nmasimscocnucsol ! for imscav of OC from nucsol
INTEGER :: nmasimscocaitsol ! for imscav of OC from Aitsol
INTEGER :: nmasimscocaccsol ! for imscav of OC from accsol
INTEGER :: nmasimscoccorsol ! for imscav of OC from corsol
INTEGER :: nmasimscocaitins ! for imscav of OC from Aitins
INTEGER :: nmasimscsonucsol ! for imscav of SO from nucsol
INTEGER :: nmasimscsoaitsol ! for imscav of SO from Aitsol
INTEGER :: nmasimscsoaccsol ! for imscav of SO from accsol
INTEGER :: nmasimscsocorsol ! for imscav of SO from corsol
INTEGER :: nmasclprsuaitsol1 ! for in-cl. ox. of SO2 by H2O2 to SU in Aitsol
INTEGER :: nmasclprsuaccsol1 ! for in-cl. ox. of SO2 by H2O2 to SU in accsol
INTEGER :: nmasclprsucorsol1 ! for in-cl. ox. of SO2 by H2O2 to SU in corsol
INTEGER :: nmasclprsuaitsol2 ! for in-cl. ox. of SO2 by O3   to SU in Aitsol
INTEGER :: nmasclprsuaccsol2 ! for in-cl. ox. of SO2 by O3   to SU in accsol
INTEGER :: nmasclprsucorsol2 ! for in-cl. ox. of SO2 by O3   to SU in corsol
INTEGER :: nmascondsunucsol ! for conden. of SU to nucsol
INTEGER :: nmascondsuaitsol ! for conden. of SU to Aitsol
INTEGER :: nmascondsuaccsol ! for conden. of SU to accsol
INTEGER :: nmascondsucorsol ! for conden. of SU to corsol
INTEGER :: nmascondsuaitins ! for conden. of SU to Aitins
INTEGER :: nmasnuclsunucsol ! for nucln. of SU to nucsol
INTEGER :: nmascondocnucsol ! for conden. of OC to nucsol
INTEGER :: nmascondocaitsol ! for conden. of OC to Aitsol
INTEGER :: nmascondocaccsol ! for conden. of OC to accsol
INTEGER :: nmascondoccorsol ! for conden. of OC to corsol
INTEGER :: nmascondocaitins ! for conden. of OC to Aitins
INTEGER :: nmascondsonucsol ! for conden. of SO to nucsol
INTEGER :: nmascondsoaitsol ! for conden. of SO to Aitsol
INTEGER :: nmascondsoaccsol ! for conden. of SO to accsol
INTEGER :: nmascondsocorsol ! for conden. of SO to corsol
INTEGER :: nmascondsoaitins ! for conden. of SO to Aitins
INTEGER :: nmascoagsuintr12 ! for inter-modal coag SU nucsol->Aitsol
INTEGER :: nmascoagsuintr13 ! for inter-modal coag SU nucsol->accsol
INTEGER :: nmascoagsuintr14 ! for inter-modal coag SU nucsol->corsol
INTEGER :: nmascoagsuintr15 ! for inter-modal coag SU nucsol->Aitins
INTEGER :: nmascoagocintr12 ! for inter-modal coag OC nucsol->Aitsol
INTEGER :: nmascoagocintr13 ! for inter-modal coag OC nucsol->accsol
INTEGER :: nmascoagocintr14 ! for inter-modal coag OC nucsol->corsol
INTEGER :: nmascoagocintr15 ! for inter-modal coag OC nucsol->Aitins
INTEGER :: nmascoagsointr12 ! for inter-modal coag SO nucsol->Aitsol
INTEGER :: nmascoagsointr13 ! for inter-modal coag SO nucsol->accsol
INTEGER :: nmascoagsointr14 ! for inter-modal coag SO nucsol->corsol
INTEGER :: nmascoagsointr15 ! for inter-modal coag SO nucsol->Aitins
INTEGER :: nmascoagsuintr23 ! for inter-modal coag SU Aitsol->accsol
INTEGER :: nmascoagbcintr23 ! for inter-modal coag BC Aitsol->accsol
INTEGER :: nmascoagocintr23 ! for inter-modal coag OC Aitsol->accsol
INTEGER :: nmascoagsointr23 ! for inter-modal coag SO Aitsol->accsol
INTEGER :: nmascoagsuintr24 ! for inter-modal coag SU Aitsol->corsol
INTEGER :: nmascoagbcintr24 ! for inter-modal coag BC Aitsol->corsol
INTEGER :: nmascoagocintr24 ! for inter-modal coag OC Aitsol->corsol
INTEGER :: nmascoagsointr24 ! for inter-modal coag SO Aitsol->corsol
INTEGER :: nmascoagsuintr34 ! for inter-modal coag SU accsol->corsol
INTEGER :: nmascoagbcintr34 ! for inter-modal coag BC accsol->corsol
INTEGER :: nmascoagocintr34 ! for inter-modal coag OC accsol->corsol
INTEGER :: nmascoagssintr34 ! for inter-modal coag SS accsol->corsol
INTEGER :: nmascoagsointr34 ! for inter-modal coag SO accsol->corsol
INTEGER :: nmascoagbcintr53 ! for inter-modal coag BC Aitins->accsol
INTEGER :: nmascoagocintr53 ! for inter-modal coag OC Aitins->accsol
INTEGER :: nmascoagbcintr54 ! for inter-modal coag BC Aitins->corsol
INTEGER :: nmascoagocintr54 ! for inter-modal coag OC Aitins->corsol
INTEGER :: nmasagedsuintr52 ! for SU ageing flux Aitins->Aitsol
INTEGER :: nmasagedbcintr52 ! for BC ageing flux Aitins->Aitsol
INTEGER :: nmasagedocintr52 ! for OC ageing flux Aitins->Aitsol
INTEGER :: nmasagedsointr52 ! for SO ageing flux Aitins->Aitsol
INTEGER :: nmasmergsuintr12 ! for SU mode-merging flux nucsol->Aitsol
INTEGER :: nmasmergocintr12 ! for OC mode-merging flux nucsol->Aitsol
INTEGER :: nmasmergsointr12 ! for SO mode-merging flux nucsol->Aitsol
INTEGER :: nmasmergsuintr23 ! for SU mode-merging flux Aitsol->accsol
INTEGER :: nmasmergbcintr23 ! for BC mode-merging flux Aitsol->accsol
INTEGER :: nmasmergocintr23 ! for OC mode-merging flux Aitsol->accsol
INTEGER :: nmasmergsointr23 ! for SO mode-merging flux Aitsol->accsol
INTEGER :: nmasmergsuintr34 ! for SU mode-merging flux accsol->corsol
INTEGER :: nmasmergssintr34 ! for SS mode-merging flux accsol->corsol
INTEGER :: nmasmergbcintr34 ! for BC mode-merging flux accsol->corsol
INTEGER :: nmasmergocintr34 ! for OC mode-merging flux accsol->corsol
INTEGER :: nmasmergsointr34 ! for SO mode-merging flux accsol->corsol
INTEGER :: nmasprocsuintr23 ! for SU cloud-processing Aitsol->accsol
INTEGER :: nmasprocbcintr23 ! for BC cloud-processing Aitsol->accsol
INTEGER :: nmasprococintr23 ! for OC cloud-processing Aitsol->accsol
INTEGER :: nmasprocsointr23 ! for SO cloud-processing Aitsol->accsol
!
! .. below are new ones for dust & modes 6/7 to be integrated
INTEGER :: nmasprimduaccsol ! for dust emissions to accsol
INTEGER :: nmasprimducorsol ! for dust emissions to corsol
INTEGER :: nmasprimduaccins ! for dust emissions to accins
INTEGER :: nmasprimducorins ! for dust emissions to corins
INTEGER :: nmasddepduaccsol ! for dust dry dep from accsol
INTEGER :: nmasddepducorsol ! for dust dry dep from corsol
INTEGER :: nmasddepduaccins ! for dust dry dep from accins
INTEGER :: nmasddepducorins ! for dust dry dep from corins
INTEGER :: nmasnuscduaccsol ! for dust nucscav from accsol
INTEGER :: nmasnuscducorsol ! for dust nucscav from corsol
INTEGER :: nmasnuscduaccins ! for dust nucscav from accins
INTEGER :: nmasnuscducorins ! for dust nucscav from corins
INTEGER :: nmasimscduaccsol ! for dust impscav from accsol
INTEGER :: nmasimscducorsol ! for dust impscav from corsol
INTEGER :: nmasimscduaccins ! for dust impscav from accins
INTEGER :: nmasimscducorins ! for dust impscav from corins
INTEGER :: nmascondsuaccins ! for conden. of SU to accins
INTEGER :: nmascondsucorins ! for conden. of SU to corins
INTEGER :: nmascondocaccins ! for conden. of OC to accins
INTEGER :: nmascondoccorins ! for conden. of OC to corins
INTEGER :: nmascondsoaccins ! for conden. of SO to accins
INTEGER :: nmascondsocorins ! for conden. of SO to corins
INTEGER :: nmascoagsuintr16 ! for inter-modal coag SU nucsol->accins
INTEGER :: nmascoagsuintr17 ! for inter-modal coag SU nucsol->corins
INTEGER :: nmascoagocintr16 ! for inter-modal coag OC nucsol->accins
INTEGER :: nmascoagocintr17 ! for inter-modal coag OC nucsol->corins
INTEGER :: nmascoagsointr16 ! for inter-modal coag SO nucsol->accins
INTEGER :: nmascoagsointr17 ! for inter-modal coag SO nucsol->corins
INTEGER :: nmascoagduintr34 ! for inter-modal coag DU accsol->corsol
INTEGER :: nmascoagduintr64 ! for inter-modal coag DU accins->corsol
INTEGER :: nmasagedsuintr63 ! for SU ageing flux accins->accsol
INTEGER :: nmasagedduintr63 ! for DU ageing flux accins->accsol
INTEGER :: nmasagedocintr63 ! for OC ageing flux accins->accsol
INTEGER :: nmasagedsointr63 ! for SO ageing flux accins->accsol
INTEGER :: nmasagedsuintr74 ! for SU ageing flux corins->corsol
INTEGER :: nmasagedduintr74 ! for DU ageing flux corins->corsol
INTEGER :: nmasagedocintr74 ! for OC ageing flux corins->corsol
INTEGER :: nmasagedsointr74 ! for SO ageing flux corins->corsol
INTEGER :: nmasmergduintr34 ! for DU mode-merging flux accsol->corsol

! Nitrate scheme budgets
! .. (indices for aerosol budget mass fluxes)
INTEGER :: nmasprimntnucsol !* for primary NO3 ems (all srces) to Nucsol
INTEGER :: nmasprimntaitsol !* for primary NO3 ems (all srces) to Aitsol
INTEGER :: nmasprimntaccsol !* for primary NO3 ems (all srces) to Accsol
INTEGER :: nmasprimntcorsol !* for primary NO3 ems (all srces) to corsol
INTEGER :: nmasprimnhnucsol !* for primary NH4 ems (all srces) to Nucsol
INTEGER :: nmasprimnhaitsol !* for primary NH4 ems (all srces) to Aitsol
INTEGER :: nmasprimnhaccsol !* for primary NH4 ems (all srces) to Accsol
INTEGER :: nmasprimnhcorsol !* for primary NH4 ems (all srces) to corsol
INTEGER :: nmascondnnaccsol !* for conden NaNO3 ems (all srces) to Accsol
INTEGER :: nmascondnncorsol !* for conden NaNO3 ems (all srces) to corsol
INTEGER :: nmasddepntaitsol !* for drydep of NO3 from aitsol
INTEGER :: nmasddepntaccsol !* for drydep of NO3 from accsol
INTEGER :: nmasddepntcorsol !* for drydep of NO3 from corsol
INTEGER :: nmasddepnhaitsol !* for drydep of NH4 from aitsol
INTEGER :: nmasddepnhaccsol !* for drydep of NH4 from accsol
INTEGER :: nmasddepnhcorsol !* for drydep of NH4 from corsol
INTEGER :: nmasddepnnaccsol !* for drydep of NaNO3 from accsol
INTEGER :: nmasddepnncorsol !* for drydep of NaNO3 from corsol
INTEGER :: nmasnuscntaitsol !* for nuscav of NO3 from aitsol
INTEGER :: nmasnuscntaccsol !* for nuscav of NO3 from accsol
INTEGER :: nmasnuscntcorsol !* for nuscav of NO3 from corsol
INTEGER :: nmasnuscnhaitsol !* for nuscav of NH4 from aitsol
INTEGER :: nmasnuscnhaccsol !* for nuscav of NH4 from accsol
INTEGER :: nmasnuscnhcorsol !* for nuscav of NH4 from corsol
INTEGER :: nmasnuscnnaccsol !* for nuscav of NaNO3 from accsol
INTEGER :: nmasnuscnncorsol !* for nuscav of NaNO3 from corsol
INTEGER :: nmasimscntaitsol !* for imscav of NO3 from aitsol
INTEGER :: nmasimscntaccsol !* for imscav of NO3 from accsol
INTEGER :: nmasimscntcorsol !* for imscav of NO3 from corsol
INTEGER :: nmasimscnhaitsol !* for imscav of NH4 from aitsol
INTEGER :: nmasimscnhaccsol !* for imscav of NH4 from accsol
INTEGER :: nmasimscnhcorsol !* for imscav of NH4 from corsol
INTEGER :: nmasimscnnaccsol !* for imscav of NaNO3 from accsol
INTEGER :: nmasimscnncorsol !* for imscav of NaNO3 from corsol
INTEGER :: nmascoagntintr23 !* for inter-modal coag NO3 Aitsol->accsol
INTEGER :: nmascoagnhintr23 !* for inter-modal coag NH4 Aitsol->accsol
INTEGER :: nmascoagntintr24 !* for inter-modal coag NO3 Aitsol->corsol
INTEGER :: nmascoagnhintr24 !* for inter-modal coag NH4 Aitsol->corsol
INTEGER :: nmascoagntintr34 !* for inter-modal coag NO3 accsol->corsol
INTEGER :: nmascoagnhintr34 !* for inter-modal coag NH4 accsol->corsol
INTEGER :: nmascoagnnintr34 !* for inter-modal coag NaNO3 accsol->corsol
INTEGER :: nmasmergntintr23 !* for NO3 mode-merging flux Aitsol->accsol
INTEGER :: nmasmergnhintr23 !* for NH4 mode-merging flux Aitsol->accsol
INTEGER :: nmasmergntintr34 !* for NO3 mode-merging flux accsol->corsol
INTEGER :: nmasmergnhintr34 !* for NH4 mode-merging flux accsol->corsol
INTEGER :: nmasmergnnintr34 !* for NaNO3 mode-merging flux accsol->corsol
INTEGER :: nmasprocntintr23 !* for NO3 cloud-processing Aitsol->accsol
INTEGER :: nmasprocnhintr23 !* for NH4 cloud-processing Aitsol->accsol

!
! .. (indices for additional diagnostics)
INTEGER :: nwtcnt1
INTEGER :: nwtcnt2
INTEGER :: nwtcnt3
INTEGER :: nwtcnt4
INTEGER :: nlcfrac
INTEGER :: nlwc
INTEGER :: nrainls
INTEGER :: nraincs
INTEGER :: nwindsp
INTEGER :: ndrydp1
INTEGER :: ndrydp2
INTEGER :: ndrydp3
INTEGER :: ndrydp4
INTEGER :: ndrydp5
INTEGER :: ndrydp6
INTEGER :: ndrydp7
! .. note below are partial volumes for SUSSBCOCSO_5MODE
INTEGER :: npvol11
INTEGER :: npvol16
INTEGER :: npvol1w
INTEGER :: npvol21
INTEGER :: npvol22
INTEGER :: npvol23
INTEGER :: npvol26
INTEGER :: npvol2w
INTEGER :: npvol31
INTEGER :: npvol32
INTEGER :: npvol33
INTEGER :: npvol34
INTEGER :: npvol36
INTEGER :: npvol3w
INTEGER :: npvol41
INTEGER :: npvol42
INTEGER :: npvol43
INTEGER :: npvol44
INTEGER :: npvol46
INTEGER :: npvol4w
INTEGER :: npvol52
INTEGER :: npvol53
!
INTEGER, PARAMETER :: nchemgmax=50 ! max value for NCHEMG
!
! Indices of aerosol components into which condensable gases go to
INTEGER :: condensable_choice(nchemgmax)

! Gas phase species which are condensable (T/F)
LOGICAL :: condensable(nchemgmax)

! Molecular masses of gas phase species (kg/mol)
REAL :: mm_gas(nchemgmax)

! Molecular diameter of condensable gas phase species (others = 0)
REAL :: dimen(nchemgmax)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_SETUP_INDICES'

CONTAINS

! ######################################################################
SUBROUTINE ukca_indices_nochem

IMPLICIT NONE

!-----------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_NOCHEM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nchemg=0           ! # of gas phase chemistry tracers
ichem=0            ! 1/0 = gas phase chemistry tracers on/off
noffox=0           ! # of offline oxidant species
nbudchem=0         ! # of gas chem budget fields
gasbudget=0        ! 1/0 = gas phase chemistry fluxes on/off
ngasbudget=0       ! # of gas phase chemistry fluxes
!
nadvg=2+nchemg     ! # of gas phase advected tracers
ntrag=nadvg+noffox ! total # of gas phase species

! For nochem, NTRAG=2, NADVG= 2
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
!
! .. below are the 40 tropospheric chemistry species
mox     = 0 ! not included
mnox    = 0 ! not included
mn2o5   = 0 ! not included
mhno4   = 0 ! not included
mhno3   = 0 ! not included
mh2o2   = 0 ! not included
mch4    = 0 ! not included
mco     = 0 ! not included
mch2o   = 0 ! not included
mmhp    = 0 ! not included
mhono   = 0 ! not included
mc2h6   = 0 ! not included
metooh  = 0 ! not included
mmecho  = 0 ! not included
mpan    = 0 ! not included
mc3h8   = 0 ! not included
mpnooh  = 0 ! not included
mpiooh  = 0 ! not included
metcho  = 0 ! not included
mme2co  = 0 ! not included
mmecoh  = 0 ! not included
mppan   = 0 ! not included
mmeno3  = 0 ! not included
moxs    = 0 ! not included
mnoys   = 0 ! not included
misop   = 0 ! not included
mc2h4   = 0 ! not included
mc2h2   = 0 ! not included
misooh  = 0 ! not included
mison   = 0 ! not included
mmacr   = 0 ! not included
mmacrooh= 0 ! not included
mmpan   = 0 ! not included
mhacet  = 0 ! not included
mmgly   = 0 ! not included
mnald   = 0 ! not included
mhcooh  = 0 ! not included
mmeco3h = 0 ! not included
mmeco2h = 0 ! not included
mmeoh   = 0 ! not included
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
msotwo  = 0 ! not included
MMeSMe  = 0 ! not included
mh2so4  = 0 ! not included
mdmso   = 0 ! not included
mmsa    = 0 ! not included
mcs2    = 0 ! not included
mh2s    = 0 ! not included
mcos    = 0 ! not included
mmonoter= 0 ! not included
msec_org= 0 ! not included
mh2o2f  = 0 ! not included
mq3d    = 1
mpt     = 2
!
! .. molar masses (kg/mol) for gases for nochem
! .. (48 dummy values)
mm_gas=[0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                       &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000]
!
condensable_choice=[0,0,0,0,0,0,0,0,0,0,0,0,                                   &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0]

condensable=(condensable_choice > 0)
!
dimen=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                        &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0]
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY
!
! .. below are the 60 indices for tropospheric chemistry species (ST)
no      = 0 ! not included
no1d    = 0 ! not included
no3     = 0 ! not included
nno     = 0 ! not included
nno3    = 0 ! not included
nno2    = 0 ! not included
nn2o5   = 0 ! not included
nhno4   = 0 ! not included
nhno3   = 0 ! not included
noh     = 0 ! not included
nho2    = 0 ! not included
nh2o2   = 0 ! not included
nch4    = 0 ! not included
nco     = 0 ! not included
nch2o   = 0 ! not included
nmeoo   = 0 ! not included
nh2o    = 0 ! not included
nmhp    = 0 ! not included
nhono   = 0 ! not included
nc2h6   = 0 ! not included
netoo   = 0 ! not included
netooh  = 0 ! not included
nmecho  = 0 ! not included
nmeco3  = 0 ! not included
npan    = 0 ! not included
nc3h8   = 0 ! not included
npnoo   = 0 ! not included
npioo   = 0 ! not included
npnooh  = 0 ! not included
npiooh  = 0 ! not included
netcho  = 0 ! not included
netco3  = 0 ! not included
nme2co  = 0 ! not included
nmecoo  = 0 ! not included
nmecoh  = 0 ! not included
nppan   = 0 ! not included
nmeno3  = 0 ! not included
nos     = 0 ! not included
no1ds   = 0 ! not included
no3s    = 0 ! not included
nnoxs   = 0 ! not included
nhno3s  = 0 ! not included
nnoys   = 0 ! not included
nisop   = 0 ! not included
nc2h4   = 0 ! not included
nc2h2   = 0 ! not included
niso2   = 0 ! not included
nisooh  = 0 ! not included
nison   = 0 ! not included
nmacr   = 0 ! not included
nmacro2 = 0 ! not included
nmacrooh= 0 ! not included
nmpan   = 0 ! not included
nhacet  = 0 ! not included
nmgly   = 0 ! not included
nnald   = 0 ! not included
nhcooh  = 0 ! not included
nmeco3h = 0 ! not included
nmeco2h = 0 ! not included
nmeoh   = 0 ! not included
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
nsotwo  = 0 ! not included
NMeSMe  = 0 ! not included
nh2so4  = 0 ! not included
ndmso   = 0 ! not included
nmsa    = 0 ! not included
ncs2    = 0 ! not included
nh2s    = 0 ! not included
ncos    = 0 ! not included
nmonoter= 0 ! not included
nsec_org= 0 ! not included
nh2o2f  = 0 ! not included
no3f    = 0 ! not included
nohf    = 0 ! not included
nno3f   = 0 ! not included
nq3d    = 1
npt     = 2
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for nochem
ndmsemoc  = 0 ! not included
ndmstend  = 0 ! not included
nso2eman  = 0 ! not included
nso2embm  = 0 ! not included
nso2emvl  = 0 ! not included
nso2tend  = 0 ! not included
nso2ddep  = 0 ! not included
nso2wdep  = 0 ! not included
nh2so4tend= 0 ! not included
nh2so4ddep= 0 ! not included
ncoseman  = 0 ! not included
ncosemoc  = 0 ! not included
ncostend  = 0 ! not included
ncs2eman  = 0 ! not included
ncs2emoc  = 0 ! not included
ncs2tend  = 0 ! not included
ndmsotend = 0 ! not included
ndmsoddep = 0 ! not included
nmsatend  = 0 ! not included
nmsaddep  = 0 ! not included
nterp_em  = 0 ! not included
nterp_tend= 0 ! not included
nterp_ddep= 0 ! not included
nsorg_tend= 0 ! not included
nsorg_ddep= 0 ! not included
nsorg_wdep= 0 ! not included
!
! REACTION INDICES for gas phase chemistry fluxes
iohdms1 =  0
iohdms2 =  0
ino3dms =  0
idmsooh1=  0
idmsooh2=  0
ics2oh  =  0
ih2soh  =  0
icosoh  =  0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_nochem

! ######################################################################
SUBROUTINE UKCA_INDICES_ORGV1_SOto3

IMPLICIT NONE

!-----------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_ORGV1_SOTO3'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nchemg=11          ! # of gas phase chemistry tracers
ichem=1            ! 1/0 = gas phase chemistry tracers on/off
noffox=3           ! # of offline oxidant species
nbudchem=26        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
gasbudget=0        ! 1/0 = gas phase chemistry fluxes on/off
ngasbudget=8       ! # of gas phase chemistry fluxes
!
nadvg=2+nchemg     ! # of gas phase advected tracers
ntrag=nadvg+noffox ! total # of gas phase species

! For orgv1, NTRAG=16, NADVG=13
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
!
! .. below are the 40 tropospheric chemistry species
mox     = 0 ! not included
mnox    = 0 ! not included
mn2o5   = 0 ! not included
mhno4   = 0 ! not included
mhno3   = 0 ! not included
mh2o2   = 0 ! not included
mch4    = 0 ! not included
mco     = 0 ! not included
mch2o   = 0 ! not included
mmhp    = 0 ! not included
mhono   = 0 ! not included
mc2h6   = 0 ! not included
metooh  = 0 ! not included
mmecho  = 0 ! not included
mpan    = 0 ! not included
mc3h8   = 0 ! not included
mpnooh  = 0 ! not included
mpiooh  = 0 ! not included
metcho  = 0 ! not included
mme2co  = 0 ! not included
mmecoh  = 0 ! not included
mppan   = 0 ! not included
mmeno3  = 0 ! not included
moxs    = 0 ! not included
mnoys   = 0 ! not included
misop   = 0 ! not included
mc2h4   = 0 ! not included
mc2h2   = 0 ! not included
misooh  = 0 ! not included
mison   = 0 ! not included
mmacr   = 0 ! not included
mmacrooh= 0 ! not included
mmpan   = 0 ! not included
mhacet  = 0 ! not included
mmgly   = 0 ! not included
mnald   = 0 ! not included
mhcooh  = 0 ! not included
mmeco3h = 0 ! not included
mmeco2h = 0 ! not included
mmeoh   = 0 ! not included
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
msotwo  = 1
MMeSMe  = 2
mh2so4  = 3
mdmso   = 4
mmsa    = 5
mcs2    = 6
mh2s    = 7
mcos    = 8
mmonoter= 9
msec_org=10
mh2o2f  =11
mq3d    =12
mpt     =13
!
! .. molar masses (kg/mol) for gases for orgv1
! .. (8 S species, terp, Sec_Org, H2O2F then 37 dummy values)
mm_gas=[0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,                       &
         0.136,0.150,0.034,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000]
!
condensable_choice=[0,0,1,0,0,0,0,0,0,3,0,0,                                   &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0]
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
! .. Sec_Org to condense into 3rd aerosol component (CP_OC)

condensable=(condensable_choice > 0)
!
dimen=[0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0,0.0,                &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0]
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY
!
! .. below are the 60 indices for tropospheric chemistry species (ST)
no      = 0 ! not included
no1d    = 0 ! not included
no3     = 0 ! not included
nno     = 0 ! not included
nno3    = 0 ! not included
nno2    = 0 ! not included
nn2o5   = 0 ! not included
nhno4   = 0 ! not included
nhno3   = 0 ! not included
noh     = 0 ! not included
nho2    = 0 ! not included
nh2o2   = 0 ! not included
nch4    = 0 ! not included
nco     = 0 ! not included
nch2o   = 0 ! not included
nmeoo   = 0 ! not included
nh2o    = 0 ! not included
nmhp    = 0 ! not included
nhono   = 0 ! not included
nc2h6   = 0 ! not included
netoo   = 0 ! not included
netooh  = 0 ! not included
nmecho  = 0 ! not included
nmeco3  = 0 ! not included
npan    = 0 ! not included
nc3h8   = 0 ! not included
npnoo   = 0 ! not included
npioo   = 0 ! not included
npnooh  = 0 ! not included
npiooh  = 0 ! not included
netcho  = 0 ! not included
netco3  = 0 ! not included
nme2co  = 0 ! not included
nmecoo  = 0 ! not included
nmecoh  = 0 ! not included
nppan   = 0 ! not included
nmeno3  = 0 ! not included
nos     = 0 ! not included
no1ds   = 0 ! not included
no3s    = 0 ! not included
nnoxs   = 0 ! not included
nhno3s  = 0 ! not included
nnoys   = 0 ! not included
nisop   = 0 ! not included
nc2h4   = 0 ! not included
nc2h2   = 0 ! not included
niso2   = 0 ! not included
nisooh  = 0 ! not included
nison   = 0 ! not included
nmacr   = 0 ! not included
nmacro2 = 0 ! not included
nmacrooh= 0 ! not included
nmpan   = 0 ! not included
nhacet  = 0 ! not included
nmgly   = 0 ! not included
nnald   = 0 ! not included
nhcooh  = 0 ! not included
nmeco3h = 0 ! not included
nmeco2h = 0 ! not included
nmeoh   = 0 ! not included
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
nsotwo  = 1
NMeSMe  = 2
nh2so4  = 3
ndmso   = 4
nmsa    = 5
ncs2    = 6
nh2s    = 7
ncos    = 8
nmonoter= 9
nsec_org=10
nh2o2f  =11
no3f    =12
nohf    =13
nno3f   =14
nq3d    =15
npt     =16
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for orgv1
ndmsemoc  = 1
ndmstend  = 2
nso2eman  = 3
nso2embm  = 4
nso2emvl  = 5
nso2tend  = 6
nso2ddep  = 7
nso2wdep  = 8
nh2so4tend= 9
nh2so4ddep=10
ncoseman  =11
ncosemoc  =12
ncostend  =13
ncs2eman  =14
ncs2emoc  =15
ncs2tend  =16
ndmsotend =17
ndmsoddep =18
nmsatend  =19
nmsaddep  =20
nterp_em  =21
nterp_tend=22
nterp_ddep=23
nsorg_tend=24
nsorg_ddep=25
nsorg_wdep=26
!
! REACTION INDICES for gas phase chemistry fluxes
iohdms1 =  1
iohdms2 =  2
ino3dms =  4
idmsooh1=  0
idmsooh2=  3
ics2oh  =  5
ih2soh  =  6
icosoh  =  7
!
!----------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE UKCA_INDICES_ORGV1_SOto3

! ######################################################################
SUBROUTINE UKCA_INDICES_ORGV1_SOto3_SOLINSOL

IMPLICIT NONE

!-----------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_ORGV1_SOTO3'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nchemg=11          ! # of gas phase chemistry tracers
ichem=1            ! 1/0 = gas phase chemistry tracers on/off
noffox=3           ! # of offline oxidant species
nbudchem=26        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
gasbudget=0        ! 1/0 = gas phase chemistry fluxes on/off
ngasbudget=8       ! # of gas phase chemistry fluxes
!
nadvg=2+nchemg     ! # of gas phase advected tracers
ntrag=nadvg+noffox ! total # of gas phase species

! For orgv1, NTRAG=16, NADVG=13
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
!
! .. below are the 40 tropospheric chemistry species
mox     = 0 ! not included
mnox    = 0 ! not included
mn2o5   = 0 ! not included
mhno4   = 0 ! not included
mhno3   = 0 ! not included
mh2o2   = 0 ! not included
mch4    = 0 ! not included
mco     = 0 ! not included
mch2o   = 0 ! not included
mmhp    = 0 ! not included
mhono   = 0 ! not included
mc2h6   = 0 ! not included
metooh  = 0 ! not included
mmecho  = 0 ! not included
mpan    = 0 ! not included
mc3h8   = 0 ! not included
mpnooh  = 0 ! not included
mpiooh  = 0 ! not included
metcho  = 0 ! not included
mme2co  = 0 ! not included
mmecoh  = 0 ! not included
mppan   = 0 ! not included
mmeno3  = 0 ! not included
moxs    = 0 ! not included
mnoys   = 0 ! not included
misop   = 0 ! not included
mc2h4   = 0 ! not included
mc2h2   = 0 ! not included
misooh  = 0 ! not included
mison   = 0 ! not included
mmacr   = 0 ! not included
mmacrooh= 0 ! not included
mmpan   = 0 ! not included
mhacet  = 0 ! not included
mmgly   = 0 ! not included
mnald   = 0 ! not included
mhcooh  = 0 ! not included
mmeco3h = 0 ! not included
mmeco2h = 0 ! not included
mmeoh   = 0 ! not included
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
msotwo  = 1
MMeSMe  = 2
mh2so4  = 3
mdmso   = 4
mmsa    = 5
mcs2    = 6
mh2s    = 7
mcos    = 8
mmonoter= 9
msec_org=10
mh2o2f  =11
mq3d    =12
mpt     =13
!
! .. molar masses (kg/mol) for gases for orgv1
! .. (8 S species, terp, Sec_Org, H2O2F then 37 dummy values)
mm_gas=[0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,                       &
         0.136,0.150,0.034,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000]
!
condensable_choice=[0,0,1,0,0,0,0,0,0,1,0,0,                                   &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0]
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
! .. Sec_Org to condense into 1st aerosol component (CP_SU)

condensable=(condensable_choice > 0)
!
dimen=[0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0,0.0,                &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0]
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY
!
! .. below are the 60 indices for tropospheric chemistry species (ST)
no      = 0 ! not included
no1d    = 0 ! not included
no3     = 0 ! not included
nno     = 0 ! not included
nno3    = 0 ! not included
nno2    = 0 ! not included
nn2o5   = 0 ! not included
nhno4   = 0 ! not included
nhno3   = 0 ! not included
noh     = 0 ! not included
nho2    = 0 ! not included
nh2o2   = 0 ! not included
nch4    = 0 ! not included
nco     = 0 ! not included
nch2o   = 0 ! not included
nmeoo   = 0 ! not included
nh2o    = 0 ! not included
nmhp    = 0 ! not included
nhono   = 0 ! not included
nc2h6   = 0 ! not included
netoo   = 0 ! not included
netooh  = 0 ! not included
nmecho  = 0 ! not included
nmeco3  = 0 ! not included
npan    = 0 ! not included
nc3h8   = 0 ! not included
npnoo   = 0 ! not included
npioo   = 0 ! not included
npnooh  = 0 ! not included
npiooh  = 0 ! not included
netcho  = 0 ! not included
netco3  = 0 ! not included
nme2co  = 0 ! not included
nmecoo  = 0 ! not included
nmecoh  = 0 ! not included
nppan   = 0 ! not included
nmeno3  = 0 ! not included
nos     = 0 ! not included
no1ds   = 0 ! not included
no3s    = 0 ! not included
nnoxs   = 0 ! not included
nhno3s  = 0 ! not included
nnoys   = 0 ! not included
nisop   = 0 ! not included
nc2h4   = 0 ! not included
nc2h2   = 0 ! not included
niso2   = 0 ! not included
nisooh  = 0 ! not included
nison   = 0 ! not included
nmacr   = 0 ! not included
nmacro2 = 0 ! not included
nmacrooh= 0 ! not included
nmpan   = 0 ! not included
nhacet  = 0 ! not included
nmgly   = 0 ! not included
nnald   = 0 ! not included
nhcooh  = 0 ! not included
nmeco3h = 0 ! not included
nmeco2h = 0 ! not included
nmeoh   = 0 ! not included
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
nsotwo  = 1
NMeSMe  = 2
nh2so4  = 3
ndmso   = 4
nmsa    = 5
ncs2    = 6
nh2s    = 7
ncos    = 8
nmonoter= 9
nsec_org=10
nh2o2f  =11
no3f    =12
nohf    =13
nno3f   =14
nq3d    =15
npt     =16
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for orgv1
ndmsemoc  = 1
ndmstend  = 2
nso2eman  = 3
nso2embm  = 4
nso2emvl  = 5
nso2tend  = 6
nso2ddep  = 7
nso2wdep  = 8
nh2so4tend= 9
nh2so4ddep=10
ncoseman  =11
ncosemoc  =12
ncostend  =13
ncs2eman  =14
ncs2emoc  =15
ncs2tend  =16
ndmsotend =17
ndmsoddep =18
nmsatend  =19
nmsaddep  =20
nterp_em  =21
nterp_tend=22
nterp_ddep=23
nsorg_tend=24
nsorg_ddep=25
nsorg_wdep=26
!
! REACTION INDICES for gas phase chemistry fluxes
iohdms1 =  1
iohdms2 =  2
ino3dms =  4
idmsooh1=  0
idmsooh2=  3
ics2oh  =  5
ih2soh  =  6
icosoh  =  7
!
!----------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE UKCA_INDICES_ORGV1_SOto3_SOLINSOL

! ######################################################################
SUBROUTINE UKCA_INDICES_ORGV1_SOto3_COUPL

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_ORGV1_SOTO3_COUPL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Main array lengths and switches

nchemg=50          ! # of gas phase chemistry tracers
ichem=1            ! 1/0 = gas phase chemistry tracers on/off
noffox=24          ! # of offline oxidant species
nbudchem=26        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
gasbudget=0        ! 1/0 = gas phase chemistry fluxes on/off
ngasbudget=8       ! # of gas phase chemistry fluxes

nadvg=2+nchemg     ! # of gas phase advected tracers
ntrag=nadvg+noffox ! total # of gas phase species


!   orgv1_soto3_coupled: NTRAG=76, NADVG=52
!
!-----------------------------------------------------------

! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
! .. below are the 40 tropospheric chemistry species
mox     = 1
mnox    = 2
mn2o5   = 3
mhno4   = 4
mhno3   = 5
mh2o2   = 6
mch4    = 7
mco     = 8
mch2o   = 9
mmhp    =10
mhono   =11
mc2h6   =12
metooh  =13
mmecho  =14
mpan    =15
mc3h8   =16
mpnooh  =17
mpiooh  =18
metcho  =19
mme2co  =20
mmecoh  =21
mppan   =22
mmeno3  =23
moxs    =24
mnoys   =25
misop   =26
mc2h4   =27
mc2h2   =28
misooh  =29
mison   =30
mmacr   =31
mmacrooh=32
mmpan   =33
mhacet  =34
mmgly   =35
mnald   =36
mhcooh  =37
mmeco3h =38
mmeco2h =39
mmeoh   =40
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
msotwo  =41
MMeSMe  =42
mh2so4  =43
mdmso   =44
mmsa    =45
mcs2    =46
mh2s    =47
mcos    =48
mmonoter=49
msec_org=50
mh2o2f  = 0 ! not included
mq3d    =51
mpt     =52
!
! .. molar masses (kg/mol) for gases for sv1_coupled
! .. (40 tropospheric gases then 8 sulphur species)
!     Dummy values used for OX =1, NOX  =2, MHP  =10,ETOOH=13, MECHO=14
!     species & no          PAN=15 PNOOH=17,PIOOH=18,ETCHO=19, ME2CO=20
!                           MECOH=21, PPAN=22, MENO3=23, OXS=24
!                           NOYS=25, ISOP=26
mm_gas=[0.048,0.030,0.108,0.079,0.053,0.034,0.016,0.013,                       &
         0.030,0.033,0.047,0.030,0.056,0.045,0.088,0.044,                      &
         0.045,0.067,0.056,0.087,0.067,0.051,0.034,0.013,                      &
         0.018,0.045,0.028,0.026,0.015,0.015,0.015,0.015,                      &
         0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,                      &
         0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,                      &
         0.136,0.150]

condensable_choice=[0,0,0,0,0,0,0,0,0,0,0,0,                                   &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,1,0,0,0,0,0,0,3]
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
! .. Sec_Org to condense into 3rd aerosol component (CP_OC)
!
condensable=(condensable_choice > 0)
!
dimen=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                        &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,                   &
        0.0,4.5e-10]
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY

! .. below are the 60 indices for tropospheric chemistry species (ST)
no      = 1
no1d    = 2
no3     = 3
nno     = 4
nno3    = 5
nno2    = 6
nn2o5   = 7
nhno4   = 8
nhno3   = 9
noh     =10
nho2    =11
nh2o2   =12
nch4    =13
nco     =14
nch2o   =15
nmeoo   =16
nh2o    =17
nmhp    =18
nhono   =19
nc2h6   =20
netoo   =21
netooh  =22
nmecho  =23
nmeco3  =24
npan    =25
nc3h8   =26
npnoo   =27
npioo   =28
npnooh  =29
npiooh  =30
netcho  =31
netco3  =32
nme2co  =33
nmecoo  =34
nmecoh  =35
nppan   =36
nmeno3  =37
nos     =38
no1ds   =39
no3s    =40
nnoxs   =41
nhno3s  =42
nnoys   =43
nisop   =44
nc2h4   =45
nc2h2   =46
niso2   =47
nisooh  =48
nison   =49
nmacr   =50
nmacro2 =51
nmacrooh=52
nmpan   =53
nhacet  =54
nmgly   =55
nnald   =56
nhcooh  =57
nmeco3h =58
nmeco2h =59
nmeoh   =60
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
nsotwo  =61
NMeSMe  =62
nh2so4  =63
ndmso   =64
nmsa    =65
ncs2    =66
nh2s    =67
ncos    =68
nmonoter=69
nsec_org=70
nh2o2f  =76
no3f    =73
nohf    =74
nno3f   =75
nq3d    =71
npt     =72
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for orgv1
ndmsemoc  = 1
ndmstend  = 2
nso2eman  = 3
nso2embm  = 4
nso2emvl  = 5
nso2tend  = 6
nso2ddep  = 7
nso2wdep  = 8
nh2so4tend= 9
nh2so4ddep=10
ncoseman  =11
ncosemoc  =12
ncostend  =13
ncs2eman  =14
ncs2emoc  =15
ncs2tend  =16
ndmsotend =17
ndmsoddep =18
nmsatend  =19
nmsaddep  =20
nterp_em  =21
nterp_tend=22
nterp_ddep=23
nsorg_tend=24
nsorg_ddep=25
nsorg_wdep=26
!
! REACTION INDICES for gas phase chemistry fluxes
iohdms1 =123
iohdms2 =124
ino3dms =127
idmsooh1=125
idmsooh2=126
ics2oh  =128
ih2soh  =129
icosoh  =130
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE UKCA_INDICES_ORGV1_SOto3_COUPL

! ######################################################################
SUBROUTINE UKCA_INDICES_ORGV1_SOto6

IMPLICIT NONE

!-----------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_ORGV1_SOTO6'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nchemg=11          ! # of gas phase chemistry tracers
ichem=1            ! 1/0 = gas phase chemistry tracers on/off
noffox=3           ! # of offline oxidant species
nbudchem=26        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
gasbudget=0        ! 1/0 = gas phase chemistry fluxes on/off
ngasbudget=8       ! # of gas phase chemistry fluxes
!
nadvg=2+nchemg     ! # of gas phase advected tracers
ntrag=nadvg+noffox ! total # of gas phase species

! For orgv1, NTRAG=16, NADVG=13
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
! .. below are the 40 tropospheric chemistry species
mox     = 0 ! not included
mnox    = 0 ! not included
mn2o5   = 0 ! not included
mhno4   = 0 ! not included
mhno3   = 0 ! not included
mh2o2   = 0 ! not included
mch4    = 0 ! not included
mco     = 0 ! not included
mch2o   = 0 ! not included
mmhp    = 0 ! not included
mhono   = 0 ! not included
mc2h6   = 0 ! not included
metooh  = 0 ! not included
mmecho  = 0 ! not included
mpan    = 0 ! not included
mc3h8   = 0 ! not included
mpnooh  = 0 ! not included
mpiooh  = 0 ! not included
metcho  = 0 ! not included
mme2co  = 0 ! not included
mmecoh  = 0 ! not included
mppan   = 0 ! not included
mmeno3  = 0 ! not included
moxs    = 0 ! not included
mnoys   = 0 ! not included
misop   = 0 ! not included
mc2h4   = 0 ! not included
mc2h2   = 0 ! not included
misooh  = 0 ! not included
mison   = 0 ! not included
mmacr   = 0 ! not included
mmacrooh= 0 ! not included
mmpan   = 0 ! not included
mhacet  = 0 ! not included
mmgly   = 0 ! not included
mnald   = 0 ! not included
mhcooh  = 0 ! not included
mmeco3h = 0 ! not included
mmeco2h = 0 ! not included
mmeoh   = 0 ! not included
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
msotwo  = 1
MMeSMe  = 2
mh2so4  = 3
mdmso   = 4
mmsa    = 5
mcs2    = 6
mh2s    = 7
mcos    = 8
mmonoter= 9
msec_org=10
mh2o2f  =11
mq3d    =12
mpt     =13
!
! .. molar masses (kg/mol) for gases for orgv1
! .. (8 S species, terp, Sec_Org, H2O2F then 37 dummy values)
mm_gas=[0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,                       &
         0.136,0.150,0.034,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000]
!
condensable_choice=[0,0,1,0,0,0,0,0,0,6,0,0,                                   &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0]
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
! .. Sec_Org to condense into 6th aerosol component (CP_SO)

condensable=(condensable_choice > 0)
!
dimen=[0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0,0.0,                &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0]
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY
!
! .. below are the 60 indices for tropospheric chemistry species (ST)
no      = 0 ! not included
no1d    = 0 ! not included
no3     = 0 ! not included
nno     = 0 ! not included
nno3    = 0 ! not included
nno2    = 0 ! not included
nn2o5   = 0 ! not included
nhno4   = 0 ! not included
nhno3   = 0 ! not included
noh     = 0 ! not included
nho2    = 0 ! not included
nh2o2   = 0 ! not included
nch4    = 0 ! not included
nco     = 0 ! not included
nch2o   = 0 ! not included
nmeoo   = 0 ! not included
nh2o    = 0 ! not included
nmhp    = 0 ! not included
nhono   = 0 ! not included
nc2h6   = 0 ! not included
netoo   = 0 ! not included
netooh  = 0 ! not included
nmecho  = 0 ! not included
nmeco3  = 0 ! not included
npan    = 0 ! not included
nc3h8   = 0 ! not included
npnoo   = 0 ! not included
npioo   = 0 ! not included
npnooh  = 0 ! not included
npiooh  = 0 ! not included
netcho  = 0 ! not included
netco3  = 0 ! not included
nme2co  = 0 ! not included
nmecoo  = 0 ! not included
nmecoh  = 0 ! not included
nppan   = 0 ! not included
nmeno3  = 0 ! not included
nos     = 0 ! not included
no1ds   = 0 ! not included
no3s    = 0 ! not included
nnoxs   = 0 ! not included
nhno3s  = 0 ! not included
nnoys   = 0 ! not included
nisop   = 0 ! not included
nc2h4   = 0 ! not included
nc2h2   = 0 ! not included
niso2   = 0 ! not included
nisooh  = 0 ! not included
nison   = 0 ! not included
nmacr   = 0 ! not included
nmacro2 = 0 ! not included
nmacrooh= 0 ! not included
nmpan   = 0 ! not included
nhacet  = 0 ! not included
nmgly   = 0 ! not included
nnald   = 0 ! not included
nhcooh  = 0 ! not included
nmeco3h = 0 ! not included
nmeco2h = 0 ! not included
nmeoh   = 0 ! not included
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
nsotwo  = 1
NMeSMe  = 2
nh2so4  = 3
ndmso   = 4
nmsa    = 5
ncs2    = 6
nh2s    = 7
ncos    = 8
nmonoter= 9
nsec_org=10
nh2o2f  =11
no3f    =12
nohf    =13
nno3f   =14
nq3d    =15
npt     =16
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for orgv1
ndmsemoc  = 1
ndmstend  = 2
nso2eman  = 3
nso2embm  = 4
nso2emvl  = 5
nso2tend  = 6
nso2ddep  = 7
nso2wdep  = 8
nh2so4tend= 9
nh2so4ddep=10
ncoseman  =11
ncosemoc  =12
ncostend  =13
ncs2eman  =14
ncs2emoc  =15
ncs2tend  =16
ndmsotend =17
ndmsoddep =18
nmsatend  =19
nmsaddep  =20
nterp_em  =21
nterp_tend=22
nterp_ddep=23
nsorg_tend=24
nsorg_ddep=25
nsorg_wdep=26

! REACTION INDICES for gas phase chemistry fluxes
iohdms1 =  1
iohdms2 =  2
ino3dms =  4
idmsooh1=  0
idmsooh2=  3
ics2oh  =  5
ih2soh  =  6
icosoh  =  7
!
!----------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE UKCA_INDICES_ORGV1_SOto6

! ######################################################################
SUBROUTINE UKCA_INDICES_ORGV1_SOto6_COUPL

IMPLICIT NONE

!-----------------------------------------------------------
! Main array lengths and switches

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_ORGV1_SOTO6_COUPL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nchemg=50          ! # of gas phase chemistry tracers
ichem=1            ! 1/0 = gas phase chemistry tracers on/off
noffox=24          ! # of offline oxidant species
nbudchem=26        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
gasbudget=0        ! 1/0 = gas phase chemistry fluxes on/off
ngasbudget=8       ! # of gas phase chemistry fluxes
!
nadvg=2+nchemg     ! # of gas phase advected tracers
ntrag=nadvg+noffox ! total # of gas phase species
!
!   orgv1_soto3_coupled: NTRAG=76, NADVG=52
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
! .. below are the 40 tropospheric chemistry species
mox     = 1
mnox    = 2
mn2o5   = 3
mhno4   = 4
mhno3   = 5
mh2o2   = 6
mch4    = 7
mco     = 8
mch2o   = 9
mmhp    =10
mhono   =11
mc2h6   =12
metooh  =13
mmecho  =14
mpan    =15
mc3h8   =16
mpnooh  =17
mpiooh  =18
metcho  =19
mme2co  =20
mmecoh  =21
mppan   =22
mmeno3  =23
moxs    =24
mnoys   =25
misop   =26
mc2h4   =27
mc2h2   =28
misooh  =29
mison   =30
mmacr   =31
mmacrooh=32
mmpan   =33
mhacet  =34
mmgly   =35
mnald   =36
mhcooh  =37
mmeco3h =38
mmeco2h =39
mmeoh   =40
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
msotwo  =41
MMeSMe  =42
mh2so4  =43
mdmso   =44
mmsa    =45
mcs2    =46
mh2s    =47
mcos    =48
mmonoter=49
msec_org=50
mh2o2f  = 0 ! not included
mq3d    =51
mpt     =52
!
! .. molar masses (kg/mol) for gases for sv1_coupled
! .. (40 tropospheric gases then 8 sulphur species)
!     Dummy values used for OX =1, NOX  =2, MHP  =10,ETOOH=13, MECHO=14
!     species & no          PAN=15 PNOOH=17,PIOOH=18,ETCHO=19, ME2CO=20
!                           MECOH=21, PPAN=22, MENO3=23, OXS=24
!                           NOYS=25, ISOP=26
mm_gas=[0.048,0.030,0.108,0.079,0.053,0.034,0.016,0.013,                       &
         0.030,0.033,0.047,0.030,0.056,0.045,0.088,0.044,                      &
         0.045,0.067,0.056,0.087,0.067,0.051,0.034,0.013,                      &
         0.018,0.045,0.028,0.026,0.015,0.015,0.015,0.015,                      &
         0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,                      &
         0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,                      &
         0.136,0.150]

condensable_choice=[0,0,0,0,0,0,0,0,0,0,0,0,                                   &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,1,0,0,0,0,0,0,6]
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
! .. Sec_Org to condense into 6th aerosol component (CP_SO)
!
condensable=(condensable_choice > 0)
!
dimen=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                        &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,                   &
        0.0,4.5e-10]
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY

! .. below are the 60 indices for tropospheric chemistry species (ST)
no      = 1
no1d    = 2
no3     = 3
nno     = 4
nno3    = 5
nno2    = 6
nn2o5   = 7
nhno4   = 8
nhno3   = 9
noh     =10
nho2    =11
nh2o2   =12
nch4    =13
nco     =14
nch2o   =15
nmeoo   =16
nh2o    =17
nmhp    =18
nhono   =19
nc2h6   =20
netoo   =21
netooh  =22
nmecho  =23
nmeco3  =24
npan    =25
nc3h8   =26
npnoo   =27
npioo   =28
npnooh  =29
npiooh  =30
netcho  =31
netco3  =32
nme2co  =33
nmecoo  =34
nmecoh  =35
nppan   =36
nmeno3  =37
nos     =38
no1ds   =39
no3s    =40
nnoxs   =41
nhno3s  =42
nnoys   =43
nisop   =44
nc2h4   =45
nc2h2   =46
niso2   =47
nisooh  =48
nison   =49
nmacr   =50
nmacro2 =51
nmacrooh=52
nmpan   =53
nhacet  =54
nmgly   =55
nnald   =56
nhcooh  =57
nmeco3h =58
nmeco2h =59
nmeoh   =60
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
nsotwo  =61
NMeSMe  =62
nh2so4  =63
ndmso   =64
nmsa    =65
ncs2    =66
nh2s    =67
ncos    =68
nmonoter=69
nsec_org=70
nh2o2f  =76
no3f    =73
nohf    =74
nno3f   =75
nq3d    =71
npt     =72
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for orgv1
ndmsemoc  = 1
ndmstend  = 2
nso2eman  = 3
nso2embm  = 4
nso2emvl  = 5
nso2tend  = 6
nso2ddep  = 7
nso2wdep  = 8
nh2so4tend= 9
nh2so4ddep=10
ncoseman  =11
ncosemoc  =12
ncostend  =13
ncs2eman  =14
ncs2emoc  =15
ncs2tend  =16
ndmsotend =17
ndmsoddep =18
nmsatend  =19
nmsaddep  =20
nterp_em  =21
nterp_tend=22
nterp_ddep=23
nsorg_tend=24
nsorg_ddep=25
nsorg_wdep=26

! REACTION INDICES for gas phase chemistry fluxes
iohdms1 =123
iohdms2 =124
ino3dms =127
idmsooh1=125
idmsooh2=126
ics2oh  =128
ih2soh  =129
icosoh  =130

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE UKCA_INDICES_ORGV1_SOto6_COUPL

! ######################################################################
SUBROUTINE ukca_indices_sv1
!
IMPLICIT NONE
!
!-----------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_SV1'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nchemg=9           ! # of gas phase chemistry tracers
ichem=1            ! 1/0 = gas phase chemistry tracers on/off
noffox=2           ! # of offline oxidant species
nbudchem=20        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
gasbudget=0        ! 1/0 = gas phase chemistry fluxes on/off
ngasbudget=8       ! # of gas phase chemistry fluxes
!
nadvg=2+nchemg     ! # of gas phase advected tracers
ntrag=nadvg+noffox ! total # of gas phase species
!
!   sv1: NTRAG=13, NADVG=11
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
! .. below are the 40 tropospheric chemistry species
mox     = 0 ! not included
mnox    = 0 ! not included
mn2o5   = 0 ! not included
mhno4   = 0 ! not included
mhno3   = 0 ! not included
mh2o2   = 0 ! not included
mch4    = 0 ! not included
mco     = 0 ! not included
mch2o   = 0 ! not included
mmhp    = 0 ! not included
mhono   = 0 ! not included
mc2h6   = 0 ! not included
metooh  = 0 ! not included
mmecho  = 0 ! not included
mpan    = 0 ! not included
mc3h8   = 0 ! not included
mpnooh  = 0 ! not included
mpiooh  = 0 ! not included
metcho  = 0 ! not included
mme2co  = 0 ! not included
mmecoh  = 0 ! not included
mppan   = 0 ! not included
mmeno3  = 0 ! not included
moxs    = 0 ! not included
mnoys   = 0 ! not included
misop   = 0 ! not included
mc2h4   = 0 ! not included
mc2h2   = 0 ! not included
misooh  = 0 ! not included
mison   = 0 ! not included
mmacr   = 0 ! not included
mmacrooh= 0 ! not included
mmpan   = 0 ! not included
mhacet  = 0 ! not included
mmgly   = 0 ! not included
mnald   = 0 ! not included
mhcooh  = 0 ! not included
mmeco3h = 0 ! not included
mmeco2h = 0 ! not included
mmeoh   = 0 ! not included
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
msotwo  = 1
MMeSMe  = 2
mh2so4  = 3
mdmso   = 4
mmsa    = 5
mcs2    = 6
mh2s    = 7
mcos    = 8
mmonoter= 0 ! not included
msec_org= 0 ! not included
mh2o2f  = 9
mq3d    =10
mpt     =11
!
! .. molar masses (kg/mol) for gases for sv1
! .. (8 S species, H2O2F then 39 dummy values)
mm_gas=[0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,                       &
         0.034,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,                      &
         0.000,0.000]
!
condensable_choice=[0,0,1,0,0,0,0,0,0,0,0,0,                                   &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0]
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
!
condensable=(condensable_choice > 0)
!
dimen=[0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                    &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0]
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY

! .. below are the 60 indices for tropospheric chemistry species (ST)
no      = 0 ! not included
no1d    = 0 ! not included
no3     = 0 ! not included
nno     = 0 ! not included
nno3    = 0 ! not included
nno2    = 0 ! not included
nn2o5   = 0 ! not included
nhno4   = 0 ! not included
nhno3   = 0 ! not included
noh     = 0 ! not included
nho2    = 0 ! not included
nh2o2   = 0 ! not included
nch4    = 0 ! not included
nco     = 0 ! not included
nch2o   = 0 ! not included
nmeoo   = 0 ! not included
nh2o    = 0 ! not included
nmhp    = 0 ! not included
nhono   = 0 ! not included
nc2h6   = 0 ! not included
netoo   = 0 ! not included
netooh  = 0 ! not included
nmecho  = 0 ! not included
nmeco3  = 0 ! not included
npan    = 0 ! not included
nc3h8   = 0 ! not included
npnoo   = 0 ! not included
npioo   = 0 ! not included
npnooh  = 0 ! not included
npiooh  = 0 ! not included
netcho  = 0 ! not included
netco3  = 0 ! not included
nme2co  = 0 ! not included
nmecoo  = 0 ! not included
nmecoh  = 0 ! not included
nppan   = 0 ! not included
nmeno3  = 0 ! not included
nos     = 0 ! not included
no1ds   = 0 ! not included
no3s    = 0 ! not included
nnoxs   = 0 ! not included
nhno3s  = 0 ! not included
nnoys   = 0 ! not included
nisop   = 0 ! not included
nc2h4   = 0 ! not included
nc2h2   = 0 ! not included
niso2   = 0 ! not included
nisooh  = 0 ! not included
nison   = 0 ! not included
nmacr   = 0 ! not included
nmacro2 = 0 ! not included
nmacrooh= 0 ! not included
nmpan   = 0 ! not included
nhacet  = 0 ! not included
nmgly   = 0 ! not included
nnald   = 0 ! not included
nhcooh  = 0 ! not included
nmeco3h = 0 ! not included
nmeco2h = 0 ! not included
nmeoh   = 0 ! not included
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
nsotwo  = 1
NMeSMe  = 2
nh2so4  = 3
ndmso   = 4
nmsa    = 5
ncs2    = 6
nh2s    = 7
ncos    = 8
nmonoter= 0 ! not included
nsec_org= 0 ! not included
nh2o2f  = 9
no3f    = 0
nohf    =10
nno3f   =11
nq3d    =12
npt     =13
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 20 gas phase budget quantity indices for sv1
ndmsemoc  = 1
ndmstend  = 2
nso2eman  = 3
nso2embm  = 4
nso2emvl  = 5
nso2tend  = 6
nso2ddep  = 7
nso2wdep  = 8
nh2so4tend= 9
nh2so4ddep=10
ncoseman  =11
ncosemoc  =12
ncostend  =13
ncs2eman  =14
ncs2emoc  =15
ncs2tend  =16
ndmsotend =17
ndmsoddep =18
nmsatend  =19
nmsaddep  =20
nterp_em  = 0 ! not included
nterp_tend= 0 ! not included
nterp_ddep= 0 ! not included
nsorg_tend= 0 ! not included
nsorg_ddep= 0 ! not included
nsorg_wdep= 0 ! not included
!
! REACTION INDICES for gas phase chemistry fluxes
iohdms1 =  1
iohdms2 =  2
ino3dms =  4
idmsooh1=  0
idmsooh2=  3
ics2oh  =  5
ih2soh  =  6
icosoh  =  7

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_sv1

! ######################################################################
SUBROUTINE ukca_indices_sv1_coupled

IMPLICIT NONE
!
!-----------------------------------------------------------

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_SV1_COUPLED'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Main array lengths and switches
!
nchemg=48          ! # of gas phase chemistry tracers
ichem=1            ! 1/0 = gas phase chemistry tracers on/off
noffox=24          ! # of offline oxidant species
nbudchem=20        ! # of gas chem budget fields
!      GASBUDGET=1        ! 1/0 = gas phase chemistry fluxes on/off
!      NGASBUDGET=8       ! # of gas phase chemistry fluxes
gasbudget=0        ! 1/0 = gas phase chemistry fluxes on/off
ngasbudget=8       ! # of gas phase chemistry fluxes
!
nadvg=2+nchemg     ! # of gas phase advected tracers
ntrag=nadvg+noffox ! total # of gas phase species
!
!   sv1_coupled: NTRAG=74, NADVG=50
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
! .. below are the 40 tropospheric chemistry species
mox     = 1
mnox    = 2
mn2o5   = 3
mhno4   = 4
mhno3   = 5
mh2o2   = 6
mch4    = 7
mco     = 8
mch2o   = 9
mmhp    =10
mhono   =11
mc2h6   =12
metooh  =13
mmecho  =14
mpan    =15
mc3h8   =16
mpnooh  =17
mpiooh  =18
metcho  =19
mme2co  =20
mmecoh  =21
mppan   =22
mmeno3  =23
moxs    =24
mnoys   =25
misop   =26
mc2h4   =27
mc2h2   =28
misooh  =29
mison   =30
mmacr   =31
mmacrooh=32
mmpan   =33
mhacet  =34
mmgly   =35
mnald   =36
mhcooh  =37
mmeco3h =38
mmeco2h =39
mmeoh   =40
!
! .. below are the  8 sulphur species, 2 organics + Q3D,PT
msotwo  =41
MMeSMe  =42
mh2so4  =43
mdmso   =44
mmsa    =45
mcs2    =46
mh2s    =47
mcos    =48
mmonoter= 0 ! not included
msec_org= 0 ! not included
mh2o2f  = 0 ! not included
mq3d    =49
mpt     =50
!
! .. molar masses (kg/mol) for gases for sv1_coupled
! .. (40 tropospheric gases then 8 sulphur species)
!     Dummy values used for OX =1, NOX  =2, MHP  =10,ETOOH=13, MECHO=14
!     species & no          PAN=15 PNOOH=17,PIOOH=18,ETCHO=19, ME2CO=20
!                           MECOH=21, PPAN=22, MENO3=23, OXS=24
!                           NOYS=25, ISOP=26
mm_gas=[0.048,0.030,0.108,0.079,0.053,0.034,0.016,0.013,                       &
         0.030,0.033,0.047,0.030,0.056,0.045,0.088,0.044,                      &
         0.045,0.067,0.056,0.087,0.067,0.051,0.034,0.013,                      &
         0.018,0.045,0.028,0.026,0.015,0.015,0.015,0.015,                      &
         0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,                      &
         0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,                      &
         0.000,0.000]

condensable_choice=[0,0,0,0,0,0,0,0,0,0,0,0,                                   &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,0,0,0,0,0,0,                                  &
                     0,0,0,0,0,0,1,0,0,0,0,0,0,0]
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
!
condensable=(condensable_choice > 0)
!
dimen=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                        &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,                       &
        0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,                   &
        0.0,0.0]
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY

! .. below are the 60 indices for tropospheric chemistry species (ST)
no      = 1
no1d    = 2
no3     = 3
nno     = 4
nno3    = 5
nno2    = 6
nn2o5   = 7
nhno4   = 8
nhno3   = 9
noh     =10
nho2    =11
nh2o2   =12
nch4    =13
nco     =14
nch2o   =15
nmeoo   =16
nh2o    =17
nmhp    =18
nhono   =19
nc2h6   =20
netoo   =21
netooh  =22
nmecho  =23
nmeco3  =24
npan    =25
nc3h8   =26
npnoo   =27
npioo   =28
npnooh  =29
npiooh  =30
netcho  =31
netco3  =32
nme2co  =33
nmecoo  =34
nmecoh  =35
nppan   =36
nmeno3  =37
nos     =38
no1ds   =39
no3s    =40
nnoxs   =41
nhno3s  =42
nnoys   =43
nisop   =44
nc2h4   =45
nc2h2   =46
niso2   =47
nisooh  =48
nison   =49
nmacr   =50
nmacro2 =51
nmacrooh=52
nmpan   =53
nhacet  =54
nmgly   =55
nnald   =56
nhcooh  =57
nmeco3h =58
nmeco2h =59
nmeoh   =60
!
! .. below are the 8 sulphur species, 2 organics
! .. H2O2F, 3 offline oxidants and Q3D,PT indices for ST
nsotwo  =61
NMeSMe  =62
nh2so4  =63
ndmso   =64
nmsa    =65
ncs2    =66
nh2s    =67
ncos    =68
nmonoter= 0 ! not included
nsec_org= 0 ! not included
nh2o2f  =74
no3f    =71
nohf    =72
nno3f   =73
nq3d    =69
npt     =70
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 20 gas phase budget quantity indices for sv1
ndmsemoc  = 1
ndmstend  = 2
nso2eman  = 3
nso2embm  = 4
nso2emvl  = 5
nso2tend  = 6
nso2ddep  = 7
nso2wdep  = 8
nh2so4tend= 9
nh2so4ddep=10
ncoseman  =11
ncosemoc  =12
ncostend  =13
ncs2eman  =14
ncs2emoc  =15
ncs2tend  =16
ndmsotend =17
ndmsoddep =18
nmsatend  =19
nmsaddep  =20
nterp_em  = 0 ! not included
nterp_tend= 0 ! not included
nterp_ddep= 0 ! not included
nsorg_tend= 0 ! not included
nsorg_ddep= 0 ! not included
nsorg_wdep= 0 ! not included

! REACTION INDICES for gas phase chemistry fluxes
iohdms1 =123
iohdms2 =124
ino3dms =127
idmsooh1=125
idmsooh2=126
ics2oh  =128
ih2soh  =129
icosoh  =130

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_sv1_coupled

! ######################################################################
SUBROUTINE ukca_indices_traqu38

IMPLICIT NONE

!-----------------------------------------------------------

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_TRAQU38'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Main array lengths and switches

traqu=1        ! 1/0 = to include extra diagnostics in ST
ntraqu=traqu*38! # of extra diagnostic fields
budget=1       ! 1/0 = to include budget terms in ST
nbudget=nbudchem+nbudaer+ngasbudget ! total # of budget terms
!
! EXTRA DIAGNOSTIC INDICES
!
! .. below are 4 water content and 5 cloud, precip & windspd fields (TRAQU)
nwtcnt1    = 1
nwtcnt2    = 2
nwtcnt3    = 3
nwtcnt4    = 4
nlcfrac    = 5
nlwc       = 6
nrainls    = 7
nraincs    = 8
nwindsp    = 9
ndrydp1    =10
ndrydp2    =11
ndrydp3    =12
ndrydp4    =13
ndrydp5    =14
ndrydp6    =15
ndrydp7    =16
! .. note below are partial volumes for SUSSBCOCSO_5MODE
npvol11    =17
npvol16    =18
npvol1w    =19
npvol21    =20
npvol22    =21
npvol23    =22
npvol26    =23
npvol2w    =24
npvol31    =25
npvol32    =26
npvol33    =27
npvol34    =28
npvol36    =29
npvol3w    =30
npvol41    =31
npvol42    =32
npvol43    =33
npvol44    =34
npvol46    =35
npvol4w    =36
npvol52    =37
npvol53    =38
!
!----------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_traqu38

! ######################################################################
SUBROUTINE ukca_indices_traqu9

IMPLICIT NONE

!-----------------------------------------------------------

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_TRAQU9'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Main array lengths and switches

traqu=1        ! 1/0 = to include extra diagnostics in ST
ntraqu=traqu*9 ! # of extra diagnostic fields
budget=1       ! 1/0 = to include budget terms in ST
nbudget=nbudchem+nbudaer+ngasbudget ! total # of budget terms
!
! EXTRA DIAGNOSTIC INDICES
!
! .. below are 4 water content and 5 cloud, precip & windspd fields (TRAQU)
nwtcnt1    = 1
nwtcnt2    = 2
nwtcnt3    = 3
nwtcnt4    = 4
nlcfrac    = 5
nlwc       = 6
nrainls    = 7
nraincs    = 8
nwindsp    = 9
ndrydp1    = 0 ! not included
ndrydp2    = 0 ! not included
ndrydp3    = 0 ! not included
ndrydp4    = 0 ! not included
ndrydp5    = 0 ! not included
ndrydp6    = 0 ! not included
ndrydp7    = 0 ! not included
! .. note below are partial volumes for SUSSBCOCSO_5MODE
npvol11    = 0 ! not included
npvol16    = 0 ! not included
npvol1w    = 0 ! not included
npvol21    = 0 ! not included
npvol22    = 0 ! not included
npvol23    = 0 ! not included
npvol26    = 0 ! not included
npvol2w    = 0 ! not included
npvol31    = 0 ! not included
npvol32    = 0 ! not included
npvol33    = 0 ! not included
npvol34    = 0 ! not included
npvol36    = 0 ! not included
npvol3w    = 0 ! not included
npvol41    = 0 ! not included
npvol42    = 0 ! not included
npvol43    = 0 ! not included
npvol44    = 0 ! not included
npvol46    = 0 ! not included
npvol4w    = 0 ! not included
npvol52    = 0 ! not included
npvol53    = 0 ! not included

!
!----------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_traqu9

! ######################################################################
SUBROUTINE ukca_indices_sussbcocdu_7mode

IMPLICIT NONE
!---------------------------------------------------------------

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_SUSSBCOCDU_7MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Main array lengths and switches

ntraer=26          ! # of aerosol advected tracers
nbudaer=138        ! # of aerosol budget fields
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  26    16    164     9 = 42+164+ 9=  215 (orgv1 ,traqu9 )
! ..  SUSSBCOC  26    16    164    38 = 42+164+38=  244 (orgv1 ,traqu38)
! ..  SUSSBCOC  26    76    164     9 =102+164+ 9=  275 (orgv1c,traqu9 )
! ..  SUSSBCOC  26    76    164    38 =102+164+38=  304 (orgv1c,traqu38)
!                        (138+26)
!
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 26    13  =  39 (orgv1 )
! ..  SUSSBCOC 26    52  =  78 (orgv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. 138 aerosol budget indices for SUSSBCOCDU_7MODE [SO in OC]
!
! .. redo these to follow Dominicks approach as categorized under
!
! .. BUD_PRIM, BUD_DDEP, BUD_NUSC, BUD_IMSC which have NBOX,NMODES,NCP
! ..
! .. then also have:
! ..
! .. BUD_CLPR --- NBOX,NMODES,2 (production of sulfate by H2O2, O3)
! .. BUD_COND --- NBOX,NMODES,2 (conden to each modes by H2SO4,Sec_Org)
! .. BUD_NUCL --- NBOX,2 (BHN and BLN)
! .. BUD_COAG --- NBOX,NMODES,NMODES,NCP
! .. BUD_AGED --- NBOX,3,NCP,2 (ageing of 3 ins modes by H2SO4,Sec_Org)
! .. BUD_MERG --- NBOX,NMODES,NCP
! .. BUD_PROC --- NBOX,NCP (processing of Aitsol mode to accsol mode)
! ..
!
nmasprimsuaitsol= 1
nmasprimsuaccsol= 2
nmasprimsucorsol= 3
nmasprimssaccsol= 4
nmasprimsscorsol= 5
nmasprimbcaitsol= 0 ! BC only emitted to insoluble
nmasprimbcaitins= 6
nmasprimocaitsol= 7
nmasprimocaitins= 8
!
nmasddepsunucsol= 9
nmasddepsuaitsol=10
nmasddepsuaccsol=11
nmasddepsucorsol=12
nmasddepssaccsol=13
nmasddepsscorsol=14
nmasddepbcaitsol=15
nmasddepbcaccsol=16
nmasddepbccorsol=17
nmasddepbcaitins=18
nmasddepocnucsol=19
nmasddepocaitsol=20
nmasddepocaccsol=21
nmasddepoccorsol=22
nmasddepocaitins=23
nmasddepsonucsol= 0 ! SO stored in OC cpt
nmasddepsoaitsol= 0 ! SO stored in OC cpt
nmasddepsoaccsol= 0 ! SO stored in OC cpt
nmasddepsocorsol= 0 ! SO stored in OC cpt
!
nmasnuscsunucsol=24
nmasnuscsuaitsol=25
nmasnuscsuaccsol=26
nmasnuscsucorsol=27
nmasnuscssaccsol=28
nmasnuscsscorsol=29
nmasnuscbcaitsol=30
nmasnuscbcaccsol=31
nmasnuscbccorsol=32
nmasnuscbcaitins=33
nmasnuscocnucsol=34
nmasnuscocaitsol=35
nmasnuscocaccsol=36
nmasnuscoccorsol=37
nmasnuscocaitins=38
nmasnuscsonucsol= 0 ! SO stored in OC cpt
nmasnuscsoaitsol= 0 ! SO stored in OC cpt
nmasnuscsoaccsol= 0 ! SO stored in OC cpt
nmasnuscsocorsol= 0 ! SO stored in OC cpt
!
nmasimscsunucsol=39
nmasimscsuaitsol=40
nmasimscsuaccsol=41
nmasimscsucorsol=42
nmasimscssaccsol=43
nmasimscsscorsol=44
nmasimscbcaitsol=45
nmasimscbcaccsol=46
nmasimscbccorsol=47
nmasimscbcaitins=48
nmasimscocnucsol=49
nmasimscocaitsol=50
nmasimscocaccsol=51
nmasimscoccorsol=52
nmasimscocaitins=53
nmasimscsonucsol= 0 ! SO stored in OC cpt
nmasimscsoaitsol= 0 ! SO stored in OC cpt
nmasimscsoaccsol= 0 ! SO stored in OC cpt
nmasimscsocorsol= 0 ! SO stored in OC cpt
!
nmasclprsuaitsol1=54
nmasclprsuaccsol1=55
nmasclprsucorsol1=56
nmasclprsuaitsol2=57
nmasclprsuaccsol2=58
nmasclprsucorsol2=59
!
nmascondsunucsol=60
nmascondsuaitsol=61
nmascondsuaccsol=62
nmascondsucorsol=63
nmascondsuaitins=64
nmasnuclsunucsol=65
nmascondocnucsol=66
nmascondocaitsol=67
nmascondocaccsol=68
nmascondoccorsol=69
nmascondocaitins=70
nmascondsonucsol= 0 ! SO stored in OC cpt
nmascondsoaitsol= 0 ! SO stored in OC cpt
nmascondsoaccsol= 0 ! SO stored in OC cpt
nmascondsocorsol= 0 ! SO stored in OC cpt
nmascondsoaitins= 0 ! SO stored in OC cpt
!
nmascoagsuintr12=71
nmascoagsuintr13=72
nmascoagsuintr14=73
nmascoagsuintr15=74
nmascoagocintr12=75
nmascoagocintr13=76
nmascoagocintr14=77
nmascoagocintr15=78
nmascoagsointr12= 0 ! stored in NMASCOAGOCINTR12
nmascoagsointr13= 0 ! stored in NMASCOAGOCINTR13
nmascoagsointr14= 0 ! stored in NMASCOAGOCINTR14
nmascoagsointr15= 0 ! stored in NMASCOAGOCINTR15
nmascoagsuintr23=79
nmascoagbcintr23=80
nmascoagocintr23=81
nmascoagsointr23= 0 ! stored in NMASCOAGOCINTR23
nmascoagsuintr24=82
nmascoagbcintr24=83
nmascoagocintr24=84
nmascoagsointr24= 0 ! stored in NMASCOAGOCINTR24
nmascoagsuintr34=85
nmascoagbcintr34=86
nmascoagocintr34=87
nmascoagssintr34=88
nmascoagsointr34= 0 ! stored in NMASCOAGOCINTR34
!
nmascoagbcintr53=89
nmascoagocintr53=90
nmascoagbcintr54=91
nmascoagocintr54=92
!
nmasagedsuintr52=93
nmasagedbcintr52=94
nmasagedocintr52=95
nmasagedsointr52= 0 ! stored in NMASAGEDOCINTR52
!
nmasmergsuintr12=96
nmasmergocintr12=97
nmasmergsointr12= 0 ! stored in NMASMERGOCINTR12
nmasmergsuintr23=98
nmasmergbcintr23=99
nmasmergocintr23=100
nmasmergsointr23=  0 ! stored in NMASMERGOCINTR12
nmasmergsuintr34=101
nmasmergssintr34=102
nmasmergbcintr34=103
nmasmergocintr34=104
nmasmergsointr34=  0 ! stored in NMASMERGOCINTR34
nmasprocsuintr23=105
nmasprocbcintr23=106
nmasprococintr23=107
nmasprocsointr23=  0 ! stored in NMASPROCOCINTR23
!
! .. below are new ones for dust & modes 6/7 to be integrated
nmasprimduaccsol=  0 ! DU only emitted to modes 6/7 here
nmasprimducorsol=  0 ! DU only emitted to modes 6/7 here
nmasprimduaccins=108
nmasprimducorins=109
nmasddepduaccsol=110
nmasddepducorsol=111
nmasddepduaccins=112
nmasddepducorins=113
nmasnuscduaccsol=114
nmasnuscducorsol=115
nmasnuscduaccins=116
nmasnuscducorins=117
nmasimscduaccsol=118
nmasimscducorsol=119
nmasimscduaccins=120
nmasimscducorins=121
nmascondsuaccins=122
nmascondsucorins=123
nmascondocaccins=124
nmascondoccorins=125
nmascondsoaccins=  0 ! secondary organic in OC cpt in this setup
nmascondsocorins=  0 ! secondary organic in OC cpt in this setup
nmascoagsuintr16=126
nmascoagsuintr17=127
nmascoagocintr16=128
nmascoagocintr17=129
nmascoagsointr16=  0 ! secondary organic in OC cpt in this setup
nmascoagsointr17=  0 ! secondary organic in OC cpt in this setup
nmascoagduintr34=130
nmascoagduintr64=131
nmasagedsuintr63=132
nmasagedduintr63=133
nmasagedocintr63=134
nmasagedsointr63=  0 ! secondary organic in OC cpt in this setup
nmasagedsuintr74=135
nmasagedduintr74=136
nmasagedocintr74=137
nmasagedsointr74=  0 ! secondary organic in OC cpt in this setup
nmasmergduintr34=138

! .. below are new ones for no3, nh4, and nano3 to be integrated
nmasprimntnucsol= 0
nmasprimntaitsol= 0
nmasprimntaccsol= 0
nmasprimntcorsol= 0
nmasprimnhnucsol= 0
nmasprimnhaitsol= 0
nmasprimnhaccsol= 0
nmasprimnhcorsol= 0
nmascondnnaccsol= 0
nmascondnncorsol= 0
nmasddepntaitsol= 0
nmasddepntaccsol= 0
nmasddepntcorsol= 0
nmasddepnhaitsol= 0
nmasddepnhaccsol= 0
nmasddepnhcorsol= 0
nmasddepnnaccsol= 0
nmasddepnncorsol= 0
nmasnuscntaitsol= 0
nmasnuscntaccsol= 0
nmasnuscntcorsol= 0
nmasnuscnhaitsol= 0
nmasnuscnhaccsol= 0
nmasnuscnhcorsol= 0
nmasnuscnnaccsol= 0
nmasnuscnncorsol= 0
nmasimscntaitsol= 0
nmasimscntaccsol= 0
nmasimscntcorsol= 0
nmasimscnhaitsol= 0
nmasimscnhaccsol= 0
nmasimscnhcorsol= 0
nmasimscnnaccsol= 0
nmasimscnncorsol= 0
nmascoagntintr23= 0
nmascoagnhintr23= 0
nmascoagntintr24= 0
nmascoagnhintr24= 0
nmascoagntintr34= 0
nmascoagnhintr34= 0
nmascoagnnintr34= 0
nmasmergntintr23= 0
nmasmergnhintr23= 0
nmasmergntintr34= 0
nmasmergnhintr34= 0
nmasmergnnintr34= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0
!
!----------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_sussbcocdu_7mode

! ######################################################################
SUBROUTINE ukca_indices_sussbcocdu_4mode

IMPLICIT NONE
!---------------------------------------------------------------

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_SUSSBCOCDU_4MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! MAIN ARRAY LENGTHS AND SWITCHES
!
ntraer=19          ! # of aerosol advected tracers
nbudaer=99         ! # of aerosol budget fields
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  19    16    125     9 = 35+125+ 9=  169 (orgv1 ,traqu9 )
! ..  SUSSBCOC  19    16    125    38 = 35+125+38=  198 (orgv1 ,traqu38)
! ..  SUSSBCOC  19    76    125     9 = 95+125+ 9=  229 (orgv1c,traqu9 )
! ..  SUSSBCOC  19    76    125    38 = 95+125+38=  258 (orgv1c,traqu38)
!                        (99+26)
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 19    13  =  32 (orgv1 )
! ..  SUSSBCOC 19    52  =  71 (orgv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! ..  99 aerosol budget indices for SUSSBCOCDU_4mode [SO in OC]
!
! .. redo these to follow Dominicks approach as categorized under
!
! .. BUD_PRIM, BUD_DDEP, BUD_NUSC, BUD_IMSC which have NBOX,NMODES,NCP
! ..
! .. then also have:
! ..
! .. BUD_CLPR --- NBOX,NMODES,2 (production of sulfate by H2O2, O3)
! .. BUD_COND --- NBOX,NMODES,2 (conden to each modes by H2SO4,Sec_Org)
! .. BUD_NUCL --- NBOX,2 (BHN and BLN)
! .. BUD_COAG --- NBOX,NMODES,NMODES,NCP
! .. BUD_AGED --- NBOX,3,NCP,2 (ageing of 3 ins modes by H2SO4,Sec_Org)
! .. BUD_MERG --- NBOX,NMODES,NCP
! .. BUD_PROC --- NBOX,NCP (processing of Aitsol mode to accsol mode)
! ..
!
nmasprimsuaitsol= 1
nmasprimsuaccsol= 2
nmasprimsucorsol= 3
nmasprimssaccsol= 4
nmasprimsscorsol= 5
nmasprimbcaitsol= 6
nmasprimbcaitins= 0 ! BC only emitted to soluble
nmasprimocaitsol= 7
nmasprimocaitins= 0 ! OC only emitted to soluble
!
nmasddepsunucsol= 8
nmasddepsuaitsol= 9
nmasddepsuaccsol=10
nmasddepsucorsol=11
nmasddepssaccsol=12
nmasddepsscorsol=13
nmasddepbcaitsol=14
nmasddepbcaccsol=15
nmasddepbccorsol=16
nmasddepbcaitins= 0 ! BC only present in soluble
nmasddepocnucsol=17
nmasddepocaitsol=18
nmasddepocaccsol=19
nmasddepoccorsol=20
nmasddepocaitins= 0 ! OC only present in soluble
nmasddepsonucsol= 0 ! SO stored in OC cpt
nmasddepsoaitsol= 0 ! SO stored in OC cpt
nmasddepsoaccsol= 0 ! SO stored in OC cpt
nmasddepsocorsol= 0 ! SO stored in OC cpt
!
nmasnuscsunucsol=21
nmasnuscsuaitsol=22
nmasnuscsuaccsol=23
nmasnuscsucorsol=24
nmasnuscssaccsol=25
nmasnuscsscorsol=26
nmasnuscbcaitsol=27
nmasnuscbcaccsol=28
nmasnuscbccorsol=29
nmasnuscbcaitins= 0 ! BC only present in soluble
nmasnuscocnucsol=30
nmasnuscocaitsol=31
nmasnuscocaccsol=32
nmasnuscoccorsol=33
nmasnuscocaitins= 0 ! OC only present in soluble
nmasnuscsonucsol= 0 ! SO stored in OC cpt
nmasnuscsoaitsol= 0 ! SO stored in OC cpt
nmasnuscsoaccsol= 0 ! SO stored in OC cpt
nmasnuscsocorsol= 0 ! SO stored in OC cpt
!
nmasimscsunucsol=34
nmasimscsuaitsol=35
nmasimscsuaccsol=36
nmasimscsucorsol=37
nmasimscssaccsol=38
nmasimscsscorsol=39
nmasimscbcaitsol=40
nmasimscbcaccsol=41
nmasimscbccorsol=42
nmasimscbcaitins= 0 ! BC only present in soluble
nmasimscocnucsol=43
nmasimscocaitsol=44
nmasimscocaccsol=45
nmasimscoccorsol=46
nmasimscocaitins= 0 ! OC only present in soluble
nmasimscsonucsol= 0 ! SO stored in OC cpt
nmasimscsoaitsol= 0 ! SO stored in OC cpt
nmasimscsoaccsol= 0 ! SO stored in OC cpt
nmasimscsocorsol= 0 ! SO stored in OC cpt
!
nmasclprsuaitsol1=47
nmasclprsuaccsol1=48
nmasclprsucorsol1=49
nmasclprsuaitsol2=50
nmasclprsuaccsol2=51
nmasclprsucorsol2=52
!
nmascondsunucsol=53
nmascondsuaitsol=54
nmascondsuaccsol=55
nmascondsucorsol=56
nmascondsuaitins= 0 ! only soluble modes
nmasnuclsunucsol=57
nmascondocnucsol=58
nmascondocaitsol=59
nmascondocaccsol=60
nmascondoccorsol=61
nmascondocaitins= 0 ! only soluble modes
nmascondsonucsol= 0 ! SO stored in OC cpt
nmascondsoaitsol= 0 ! SO stored in OC cpt
nmascondsoaccsol= 0 ! SO stored in OC cpt
nmascondsocorsol= 0 ! SO stored in OC cpt
nmascondsoaitins= 0 ! SO stored in OC cpt
!
nmascoagsuintr12=62
nmascoagsuintr13=63
nmascoagsuintr14=64
nmascoagsuintr15= 0 ! only soluble modes
nmascoagocintr12=65
nmascoagocintr13=66
nmascoagocintr14=67
nmascoagocintr15= 0 ! only soluble modes
nmascoagsointr12= 0 ! stored in NMASCOAGOCINTR12
nmascoagsointr13= 0 ! stored in NMASCOAGOCINTR13
nmascoagsointr14= 0 ! stored in NMASCOAGOCINTR14
nmascoagsointr15= 0 ! stored in NMASCOAGOCINTR15
nmascoagsuintr23=68
nmascoagbcintr23=69
nmascoagocintr23=70
nmascoagsointr23= 0 ! stored in NMASCOAGOCINTR23
nmascoagsuintr24=71
nmascoagbcintr24=72
nmascoagocintr24=73
nmascoagsointr24= 0 ! stored in NMASCOAGOCINTR24
nmascoagsuintr34=74
nmascoagbcintr34=75
nmascoagocintr34=76
nmascoagssintr34=77
nmascoagsointr34= 0 ! stored in NMASCOAGOCINTR34
!
nmascoagbcintr53= 0 ! only soluble modes
nmascoagocintr53= 0 ! only soluble modes
nmascoagbcintr54= 0 ! only soluble modes
nmascoagocintr54= 0 ! only soluble modes
!
nmasagedsuintr52= 0 ! only soluble modes
nmasagedbcintr52= 0 ! only soluble modes
nmasagedocintr52= 0 ! only soluble modes
nmasagedsointr52= 0 ! only soluble modes
!
nmasmergsuintr12=78
nmasmergocintr12=79
nmasmergsointr12= 0 ! stored in NMASMERGOCINTR12
nmasmergsuintr23=80
nmasmergbcintr23=81
nmasmergocintr23=82
nmasmergsointr23= 0 ! stored in NMASMERGOCINTR12
nmasmergsuintr34=83
nmasmergssintr34=84
nmasmergbcintr34=85
nmasmergocintr34=86
nmasmergsointr34= 0 ! stored in NMASMERGOCINTR34
nmasprocsuintr23=87
nmasprocbcintr23=88
nmasprococintr23=89
nmasprocsointr23= 0 ! stored in NMASPROCOCINTR23
!
! .. below are new ones for dust & modes 6/7 to be integrated
nmasprimduaccsol= 90
nmasprimducorsol= 91
nmasprimduaccins=  0 ! no acc-ins cor-ins in this setup
nmasprimducorins=  0 ! no acc-ins cor-ins in this setup
nmasddepduaccsol= 92
nmasddepducorsol= 93
nmasddepduaccins=  0 ! no acc-ins cor-ins in this setup
nmasddepducorins=  0 ! no acc-ins cor-ins in this setup
nmasnuscduaccsol= 94
nmasnuscducorsol= 95
nmasnuscduaccins=  0 ! no acc-ins cor-ins in this setup
nmasnuscducorins=  0 ! no acc-ins cor-ins in this setup
nmasimscduaccsol= 96
nmasimscducorsol= 97
nmasimscduaccins=  0 ! no acc-ins cor-ins in this setup
nmasimscducorins=  0 ! no acc-ins cor-ins in this setup
nmascondsuaccins=  0 ! no acc-ins cor-ins in this setup
nmascondsucorins=  0 ! no acc-ins cor-ins in this setup
nmascondocaccins=  0 ! no acc-ins cor-ins in this setup
nmascondoccorins=  0 ! no acc-ins cor-ins in this setup
nmascondsoaccins=  0 ! no acc-ins cor-ins in this setup
nmascondsocorins=  0 ! no acc-ins cor-ins in this setup
nmascoagsuintr16=  0 ! no acc-ins cor-ins in this setup
nmascoagsuintr17=  0 ! no acc-ins cor-ins in this setup
nmascoagocintr16=  0 ! no acc-ins cor-ins in this setup
nmascoagocintr17=  0 ! no acc-ins cor-ins in this setup
nmascoagsointr16=  0 ! no acc-ins cor-ins in this setup
nmascoagsointr17=  0 ! no acc-ins cor-ins in this setup
nmascoagduintr34= 98
nmascoagduintr64=  0 ! no acc-ins cor-ins in this setup
nmasagedsuintr63=  0 ! no acc-ins cor-ins in this setup
nmasagedduintr63=  0 ! no acc-ins cor-ins in this setup
nmasagedocintr63=  0 ! no acc-ins cor-ins in this setup
nmasagedsointr63=  0 ! no acc-ins cor-ins in this setup
nmasagedsuintr74=  0 ! no acc-ins cor-ins in this setup
nmasagedduintr74=  0 ! no acc-ins cor-ins in this setup
nmasagedocintr74=  0 ! no acc-ins cor-ins in this setup
nmasagedsointr74=  0 ! no acc-ins cor-ins in this setup
nmasmergduintr34= 99

! .. below are new ones for no3, nh4, and nano3 to be integrated
nmasprimntnucsol= 0
nmasprimntaitsol= 0
nmasprimntaccsol= 0
nmasprimntcorsol= 0
nmasprimnhnucsol= 0
nmasprimnhaitsol= 0
nmasprimnhaccsol= 0
nmasprimnhcorsol= 0
nmascondnnaccsol= 0
nmascondnncorsol= 0
nmasddepntaitsol= 0
nmasddepntaccsol= 0
nmasddepntcorsol= 0
nmasddepnhaitsol= 0
nmasddepnhaccsol= 0
nmasddepnhcorsol= 0
nmasddepnnaccsol= 0
nmasddepnncorsol= 0
nmasnuscntaitsol= 0
nmasnuscntaccsol= 0
nmasnuscntcorsol= 0
nmasnuscnhaitsol= 0
nmasnuscnhaccsol= 0
nmasnuscnhcorsol= 0
nmasnuscnnaccsol= 0
nmasnuscnncorsol= 0
nmasimscntaitsol= 0
nmasimscntaccsol= 0
nmasimscntcorsol= 0
nmasimscnhaitsol= 0
nmasimscnhaccsol= 0
nmasimscnhcorsol= 0
nmasimscnnaccsol= 0
nmasimscnncorsol= 0
nmascoagntintr23= 0
nmascoagnhintr23= 0
nmascoagntintr24= 0
nmascoagnhintr24= 0
nmascoagntintr34= 0
nmascoagnhintr34= 0
nmascoagnnintr34= 0
nmasmergntintr23= 0
nmasmergnhintr23= 0
nmasmergntintr34= 0
nmasmergnhintr34= 0
nmasmergnnintr34= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0
!
!----------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_sussbcocdu_4mode

! ######################################################################
SUBROUTINE ukca_indices_sussbcoc_5mode

IMPLICIT NONE
!---------------------------------------------------------------
!
! Main array lengths and switches
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_SUSSBCOC_5MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ntraer=20          ! # of aerosol advected tracers
nbudaer=107        ! # of aerosol budget fields
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  20    16    133     9 = 36+133+ 9=  178 (orgv1 ,traqu9 )
! ..  SUSSBCOC  20    16    133    38 = 36+133+38=  207 (orgv1 ,traqu38)
! ..  SUSSBCOC  20    76    133     9 = 96+133+ 9=  238 (orgv1c,traqu9 )
! ..  SUSSBCOC  20    76    133    38 = 96+133+38=  267 (orgv1c,traqu38)
!                         (107+26)
!
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 20    13  =  33 (orgv1 )
! ..  SUSSBCOC 20    52  =  72 (orgv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. below are 107 aerosol budget indices for SUSSBCOC [SO in OC]
!
! .. redo these to follow Dominicks approach as categorized under
!
! .. BUD_PRIM, BUD_DDEP, BUD_NUSC, BUD_IMSC which have NBOX,NMODES,NCP
! ..
! .. then also have:
! ..
! .. BUD_CLPR --- NBOX,NMODES,2 (production of sulfate by H2O2, O3)
! .. BUD_COND --- NBOX,NMODES,2 (conden to each modes by H2SO4,Sec_Org)
! .. BUD_NUCL --- NBOX,2 (BHN and BLN)
! .. BUD_COAG --- NBOX,NMODES,NMODES,NCP
! .. BUD_AGED --- NBOX,3,NCP,2 (ageing of 3 ins modes by H2SO4,Sec_Org)
! .. BUD_MERG --- NBOX,NMODES,NCP
! .. BUD_PROC --- NBOX,NCP (processing of Aitsol mode to accsol mode)
! ..
!
nmasprimsuaitsol= 1
nmasprimsuaccsol= 2
nmasprimsucorsol= 3
nmasprimssaccsol= 4
nmasprimsscorsol= 5
nmasprimbcaitsol= 0 ! BC only emitted to insoluble
nmasprimbcaitins= 6
nmasprimocaitsol= 7
nmasprimocaitins= 8
!
nmasddepsunucsol= 9
nmasddepsuaitsol=10
nmasddepsuaccsol=11
nmasddepsucorsol=12
nmasddepssaccsol=13
nmasddepsscorsol=14
nmasddepbcaitsol=15
nmasddepbcaccsol=16
nmasddepbccorsol=17
nmasddepbcaitins=18
nmasddepocnucsol=19
nmasddepocaitsol=20
nmasddepocaccsol=21
nmasddepoccorsol=22
nmasddepocaitins=23
nmasddepsonucsol= 0 ! SO stored in OC cpt
nmasddepsoaitsol= 0 ! SO stored in OC cpt
nmasddepsoaccsol= 0 ! SO stored in OC cpt
nmasddepsocorsol= 0 ! SO stored in OC cpt
!
nmasnuscsunucsol=24
nmasnuscsuaitsol=25
nmasnuscsuaccsol=26
nmasnuscsucorsol=27
nmasnuscssaccsol=28
nmasnuscsscorsol=29
nmasnuscbcaitsol=30
nmasnuscbcaccsol=31
nmasnuscbccorsol=32
nmasnuscbcaitins=33
nmasnuscocnucsol=34
nmasnuscocaitsol=35
nmasnuscocaccsol=36
nmasnuscoccorsol=37
nmasnuscocaitins=38
nmasnuscsonucsol= 0 ! SO stored in OC cpt
nmasnuscsoaitsol= 0 ! SO stored in OC cpt
nmasnuscsoaccsol= 0 ! SO stored in OC cpt
nmasnuscsocorsol= 0 ! SO stored in OC cpt
!
nmasimscsunucsol=39
nmasimscsuaitsol=40
nmasimscsuaccsol=41
nmasimscsucorsol=42
nmasimscssaccsol=43
nmasimscsscorsol=44
nmasimscbcaitsol=45
nmasimscbcaccsol=46
nmasimscbccorsol=47
nmasimscbcaitins=48
nmasimscocnucsol=49
nmasimscocaitsol=50
nmasimscocaccsol=51
nmasimscoccorsol=52
nmasimscocaitins=53
nmasimscsonucsol= 0 ! SO stored in OC cpt
nmasimscsoaitsol= 0 ! SO stored in OC cpt
nmasimscsoaccsol= 0 ! SO stored in OC cpt
nmasimscsocorsol= 0 ! SO stored in OC cpt
!
nmasclprsuaitsol1=54
nmasclprsuaccsol1=55
nmasclprsucorsol1=56
nmasclprsuaitsol2=57
nmasclprsuaccsol2=58
nmasclprsucorsol2=59
!
nmascondsunucsol=60
nmascondsuaitsol=61
nmascondsuaccsol=62
nmascondsucorsol=63
nmascondsuaitins=64
nmasnuclsunucsol=65
nmascondocnucsol=66
nmascondocaitsol=67
nmascondocaccsol=68
nmascondoccorsol=69
nmascondocaitins=70
nmascondsonucsol= 0 ! SO stored in OC cpt
nmascondsoaitsol= 0 ! SO stored in OC cpt
nmascondsoaccsol= 0 ! SO stored in OC cpt
nmascondsocorsol= 0 ! SO stored in OC cpt
nmascondsoaitins= 0 ! SO stored in OC cpt
!
nmascoagsuintr12=71
nmascoagsuintr13=72
nmascoagsuintr14=73
nmascoagsuintr15=74
nmascoagocintr12=75
nmascoagocintr13=76
nmascoagocintr14=77
nmascoagocintr15=78
nmascoagsointr12= 0 ! stored in NMASCOAGOCINTR12
nmascoagsointr13= 0 ! stored in NMASCOAGOCINTR13
nmascoagsointr14= 0 ! stored in NMASCOAGOCINTR14
nmascoagsointr15= 0 ! stored in NMASCOAGOCINTR15
nmascoagsuintr23=79
nmascoagbcintr23=80
nmascoagocintr23=81
nmascoagsointr23= 0 ! stored in NMASCOAGOCINTR23
nmascoagsuintr24=82
nmascoagbcintr24=83
nmascoagocintr24=84
nmascoagsointr24= 0 ! stored in NMASCOAGOCINTR24
nmascoagsuintr34=85
nmascoagbcintr34=86
nmascoagocintr34=87
nmascoagssintr34=88
nmascoagsointr34= 0 ! stored in NMASCOAGOCINTR34
!
nmascoagbcintr53=89
nmascoagocintr53=90
nmascoagbcintr54=91
nmascoagocintr54=92
!
nmasagedsuintr52=93
nmasagedbcintr52=94
nmasagedocintr52=95
nmasagedsointr52= 0 ! stored in NMASAGEDOCINTR52
!
nmasmergsuintr12=96
nmasmergocintr12=97
nmasmergsointr12= 0 ! stored in NMASMERGOCINTR12
nmasmergsuintr23=98
nmasmergbcintr23=99
nmasmergocintr23=100
nmasmergsointr23=  0 ! stored in NMASMERGOCINTR12
nmasmergsuintr34=101
nmasmergssintr34=102
nmasmergbcintr34=103
nmasmergocintr34=104
nmasmergsointr34=  0 ! stored in NMASMERGOCINTR34
nmasprocsuintr23=105
nmasprocbcintr23=106
nmasprococintr23=107
nmasprocsointr23=  0 ! stored in NMASPROCOCINTR23

! .. below are new ones for dust & modes 6/7 to be integrated
nmasprimduaccsol= 0 ! no DU in this setup
nmasprimducorsol= 0 ! no DU in this setup
nmasprimduaccins= 0 ! no DU in this setup
nmasprimducorins= 0 ! no DU in this setup
nmasddepduaccsol= 0 ! no DU in this setup
nmasddepducorsol= 0 ! no DU in this setup
nmasddepduaccins= 0 ! no DU in this setup
nmasddepducorins= 0 ! no DU in this setup
nmasnuscduaccsol= 0 ! no DU in this setup
nmasnuscducorsol= 0 ! no DU in this setup
nmasnuscduaccins= 0 ! no DU in this setup
nmasnuscducorins= 0 ! no DU in this setup
nmasimscduaccsol= 0 ! no DU in this setup
nmasimscducorsol= 0 ! no DU in this setup
nmasimscduaccins= 0 ! no DU in this setup
nmasimscducorins= 0 ! no DU in this setup
nmascondsuaccins= 0 ! no acc-ins cor-ins in this setup
nmascondsucorins= 0 ! no acc-ins cor-ins in this setup
nmascondocaccins= 0 ! no acc-ins cor-ins in this setup
nmascondoccorins= 0 ! no acc-ins cor-ins in this setup
nmascondsoaccins= 0 ! no acc-ins cor-ins in this setup
nmascondsocorins= 0 ! no acc-ins cor-ins in this setup
nmascoagsuintr16= 0 ! no acc-ins cor-ins in this setup
nmascoagsuintr17= 0 ! no acc-ins cor-ins in this setup
nmascoagocintr16= 0 ! no acc-ins cor-ins in this setup
nmascoagocintr17= 0 ! no acc-ins cor-ins in this setup
nmascoagsointr16= 0 ! no acc-ins cor-ins in this setup
nmascoagsointr17= 0 ! no acc-ins cor-ins in this setup
nmascoagduintr34= 0 ! no DU in this setup
nmascoagduintr64= 0 ! no DU in this setup
nmasagedsuintr63= 0 ! no acc-ins cor-ins in this setup
nmasagedduintr63= 0 ! no DU in this setup
nmasagedocintr63= 0 ! no BC/OC/SO in this setup
nmasagedsointr63= 0 ! no BC/OC/SO in this setup
nmasagedsuintr74= 0 ! no acc-ins cor-ins in this setup
nmasagedduintr74= 0 ! no DU in this setup
nmasagedocintr74= 0 ! no BC/OC/SO in this setup
nmasagedsointr74= 0 ! no BC/OC/SO in this setup
nmasmergduintr34= 0 ! no DU in this setup

! .. below are new ones for no3, nh4, and nano3 to be integrated
nmasprimntnucsol= 0
nmasprimntaitsol= 0
nmasprimntaccsol= 0
nmasprimntcorsol= 0
nmasprimnhnucsol= 0
nmasprimnhaitsol= 0
nmasprimnhaccsol= 0
nmasprimnhcorsol= 0
nmascondnnaccsol= 0
nmascondnncorsol= 0
nmasddepntaitsol= 0
nmasddepntaccsol= 0
nmasddepntcorsol= 0
nmasddepnhaitsol= 0
nmasddepnhaccsol= 0
nmasddepnhcorsol= 0
nmasddepnnaccsol= 0
nmasddepnncorsol= 0
nmasnuscntaitsol= 0
nmasnuscntaccsol= 0
nmasnuscntcorsol= 0
nmasnuscnhaitsol= 0
nmasnuscnhaccsol= 0
nmasnuscnhcorsol= 0
nmasnuscnnaccsol= 0
nmasnuscnncorsol= 0
nmasimscntaitsol= 0
nmasimscntaccsol= 0
nmasimscntcorsol= 0
nmasimscnhaitsol= 0
nmasimscnhaccsol= 0
nmasimscnhcorsol= 0
nmasimscnnaccsol= 0
nmasimscnncorsol= 0
nmascoagntintr23= 0
nmascoagnhintr23= 0
nmascoagntintr24= 0
nmascoagnhintr24= 0
nmascoagntintr34= 0
nmascoagnhintr34= 0
nmascoagnnintr34= 0
nmasmergntintr23= 0
nmasmergnhintr23= 0
nmasmergntintr34= 0
nmasmergnhintr34= 0
nmasmergnnintr34= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_sussbcoc_5mode

! ######################################################################
SUBROUTINE ukca_indices_sussbcoc_4mode

IMPLICIT NONE
!---------------------------------------------------------------
!
! Main array lengths and switches

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_SUSSBCOC_4MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ntraer=17          ! # of aerosol advected tracers
nbudaer=89         ! # of aerosol budget fields
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  17    16    115     9 = 33+115+ 9=  157 (orgv1 ,traqu9 )
! ..  SUSSBCOC  17    16    115    38 = 33+115+38=  186 (orgv1 ,traqu38)
! ..  SUSSBCOC  17    76    115     9 = 63+115+ 9=  217 (orgv1c,traqu9 )
! ..  SUSSBCOC  17    76    115    38 = 63+115+38=  246 (orgv1c,traqu38)
!                        (89+26)
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 17    13  =  30
! ..  SUSSBCOC 17    52  =  69
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! ..  89 aerosol budget indices for SUSSBCOC_4MODE [SO in OC]
!
! .. redo these to follow Dominicks approach as categorized under
!
! .. BUD_PRIM, BUD_DDEP, BUD_NUSC, BUD_IMSC which have NBOX,NMODES,NCP
! ..
! .. then also have:
! ..
! .. BUD_CLPR --- NBOX,NMODES,2 (production of sulfate by H2O2, O3)
! .. BUD_COND --- NBOX,NMODES,2 (conden to each modes by H2SO4,Sec_Org)
! .. BUD_NUCL --- NBOX,2 (BHN and BLN)
! .. BUD_COAG --- NBOX,NMODES,NMODES,NCP
! .. BUD_AGED --- NBOX,3,NCP,2 (ageing of 3 ins modes by H2SO4,Sec_Org)
! .. BUD_MERG --- NBOX,NMODES,NCP
! .. BUD_PROC --- NBOX,NCP (processing of Aitsol mode to accsol mode)
! ..
!
nmasprimsuaitsol= 1
nmasprimsuaccsol= 2
nmasprimsucorsol= 3
nmasprimssaccsol= 4
nmasprimsscorsol= 5
nmasprimbcaitsol= 6
nmasprimbcaitins= 0 ! BC only emitted to soluble
nmasprimocaitsol= 7
nmasprimocaitins= 0 ! OC only emitted to soluble
!
nmasddepsunucsol= 8
nmasddepsuaitsol= 9
nmasddepsuaccsol=10
nmasddepsucorsol=11
nmasddepssaccsol=12
nmasddepsscorsol=13
nmasddepbcaitsol=14
nmasddepbcaccsol=15
nmasddepbccorsol=16
nmasddepbcaitins= 0 ! BC only present in soluble
nmasddepocnucsol=17
nmasddepocaitsol=18
nmasddepocaccsol=19
nmasddepoccorsol=20
nmasddepocaitins= 0 ! OC only present in soluble
nmasddepsonucsol= 0 ! SO stored in OC cpt
nmasddepsoaitsol= 0 ! SO stored in OC cpt
nmasddepsoaccsol= 0 ! SO stored in OC cpt
nmasddepsocorsol= 0 ! SO stored in OC cpt
!
nmasnuscsunucsol=21
nmasnuscsuaitsol=22
nmasnuscsuaccsol=23
nmasnuscsucorsol=24
nmasnuscssaccsol=25
nmasnuscsscorsol=26
nmasnuscbcaitsol=27
nmasnuscbcaccsol=28
nmasnuscbccorsol=29
nmasnuscbcaitins= 0 ! BC only present in soluble
nmasnuscocnucsol=30
nmasnuscocaitsol=31
nmasnuscocaccsol=32
nmasnuscoccorsol=33
nmasnuscocaitins= 0 ! OC only present in soluble
nmasnuscsonucsol= 0 ! SO stored in OC cpt
nmasnuscsoaitsol= 0 ! SO stored in OC cpt
nmasnuscsoaccsol= 0 ! SO stored in OC cpt
nmasnuscsocorsol= 0 ! SO stored in OC cpt
!
nmasimscsunucsol=34
nmasimscsuaitsol=35
nmasimscsuaccsol=36
nmasimscsucorsol=37
nmasimscssaccsol=38
nmasimscsscorsol=39
nmasimscbcaitsol=40
nmasimscbcaccsol=41
nmasimscbccorsol=42
nmasimscbcaitins= 0 ! BC only present in soluble
nmasimscocnucsol=43
nmasimscocaitsol=44
nmasimscocaccsol=45
nmasimscoccorsol=46
nmasimscocaitins= 0 ! OC only present in soluble
nmasimscsonucsol= 0 ! SO stored in OC cpt
nmasimscsoaitsol= 0 ! SO stored in OC cpt
nmasimscsoaccsol= 0 ! SO stored in OC cpt
nmasimscsocorsol= 0 ! SO stored in OC cpt
!
nmasclprsuaitsol1=47
nmasclprsuaccsol1=48
nmasclprsucorsol1=49
nmasclprsuaitsol2=50
nmasclprsuaccsol2=51
nmasclprsucorsol2=52
!
nmascondsunucsol=53
nmascondsuaitsol=54
nmascondsuaccsol=55
nmascondsucorsol=56
nmascondsuaitins= 0 ! only soluble modes
nmasnuclsunucsol=57
nmascondocnucsol=58
nmascondocaitsol=59
nmascondocaccsol=60
nmascondoccorsol=61
nmascondocaitins= 0 ! only soluble modes
nmascondsonucsol= 0 ! SO stored in OC cpt
nmascondsoaitsol= 0 ! SO stored in OC cpt
nmascondsoaccsol= 0 ! SO stored in OC cpt
nmascondsocorsol= 0 ! SO stored in OC cpt
nmascondsoaitins= 0 ! SO stored in OC cpt
!
nmascoagsuintr12=62
nmascoagsuintr13=63
nmascoagsuintr14=64
nmascoagsuintr15= 0 ! only soluble modes
nmascoagocintr12=65
nmascoagocintr13=66
nmascoagocintr14=67
nmascoagocintr15= 0 ! only soluble modes
nmascoagsointr12= 0 ! stored in NMASCOAGOCINTR12
nmascoagsointr13= 0 ! stored in NMASCOAGOCINTR13
nmascoagsointr14= 0 ! stored in NMASCOAGOCINTR14
nmascoagsointr15= 0 ! stored in NMASCOAGOCINTR15
nmascoagsuintr23=68
nmascoagbcintr23=69
nmascoagocintr23=70
nmascoagsointr23= 0 ! stored in NMASCOAGOCINTR23
nmascoagsuintr24=71
nmascoagbcintr24=72
nmascoagocintr24=73
nmascoagsointr24= 0 ! stored in NMASCOAGOCINTR24
nmascoagsuintr34=74
nmascoagbcintr34=75
nmascoagocintr34=76
nmascoagssintr34=77
nmascoagsointr34= 0 ! stored in NMASCOAGOCINTR34
!
nmascoagbcintr53= 0 ! only soluble modes
nmascoagocintr53= 0 ! only soluble modes
nmascoagbcintr54= 0 ! only soluble modes
nmascoagocintr54= 0 ! only soluble modes
!
nmasagedsuintr52= 0 ! only soluble modes
nmasagedbcintr52= 0 ! only soluble modes
nmasagedocintr52= 0 ! only soluble modes
nmasagedsointr52= 0 ! only soluble modes
!
nmasmergsuintr12=78
nmasmergocintr12=79
nmasmergsointr12= 0 ! stored in NMASMERGOCINTR12
nmasmergsuintr23=80
nmasmergbcintr23=81
nmasmergocintr23=82
nmasmergsointr23= 0 ! stored in NMASMERGOCINTR12
nmasmergsuintr34=83
nmasmergssintr34=84
nmasmergbcintr34=85
nmasmergocintr34=86
nmasmergsointr34= 0 ! stored in NMASMERGOCINTR34
nmasprocsuintr23=87
nmasprocbcintr23=88
nmasprococintr23=89
nmasprocsointr23= 0 ! stored in NMASPROCOCINTR23
!
! .. below are new ones for dust & modes 6/7 to be integrated
nmasprimduaccsol= 0 ! no DU in this setup
nmasprimducorsol= 0 ! no DU in this setup
nmasprimduaccins= 0 ! no DU in this setup
nmasprimducorins= 0 ! no DU in this setup
nmasddepduaccsol= 0 ! no DU in this setup
nmasddepducorsol= 0 ! no DU in this setup
nmasddepduaccins= 0 ! no DU in this setup
nmasddepducorins= 0 ! no DU in this setup
nmasnuscduaccsol= 0 ! no DU in this setup
nmasnuscducorsol= 0 ! no DU in this setup
nmasnuscduaccins= 0 ! no DU in this setup
nmasnuscducorins= 0 ! no DU in this setup
nmasimscduaccsol= 0 ! no DU in this setup
nmasimscducorsol= 0 ! no DU in this setup
nmasimscduaccins= 0 ! no DU in this setup
nmasimscducorins= 0 ! no DU in this setup
nmascondsuaccins= 0 ! no acc-ins cor-ins in this setup
nmascondsucorins= 0 ! no acc-ins cor-ins in this setup
nmascondocaccins= 0 ! no acc-ins cor-ins in this setup
nmascondoccorins= 0 ! no acc-ins cor-ins in this setup
nmascondsoaccins= 0 ! no acc-ins cor-ins in this setup
nmascondsocorins= 0 ! no acc-ins cor-ins in this setup
nmascoagsuintr16= 0 ! no acc-ins cor-ins in this setup
nmascoagsuintr17= 0 ! no acc-ins cor-ins in this setup
nmascoagocintr16= 0 ! no acc-ins cor-ins in this setup
nmascoagocintr17= 0 ! no acc-ins cor-ins in this setup
nmascoagsointr16= 0 ! no acc-ins cor-ins in this setup
nmascoagsointr17= 0 ! no acc-ins cor-ins in this setup
nmascoagduintr34= 0 ! no DU in this setup
nmascoagduintr64= 0 ! no DU in this setup
nmasagedsuintr63= 0 ! no acc-ins cor-ins in this setup
nmasagedduintr63= 0 ! no DU in this setup
nmasagedocintr63= 0 ! no BC/OC/SO in this setup
nmasagedsointr63= 0 ! no BC/OC/SO in this setup
nmasagedsuintr74= 0 ! no acc-ins cor-ins in this setup
nmasagedduintr74= 0 ! no DU in this setup
nmasagedocintr74= 0 ! no BC/OC/SO in this setup
nmasagedsointr74= 0 ! no BC/OC/SO in this setup
nmasmergduintr34= 0 ! no DU in this setup

! .. below are new ones for no3, nh4, and nano3 to be integrated
nmasprimntnucsol= 0
nmasprimntaitsol= 0
nmasprimntaccsol= 0
nmasprimntcorsol= 0
nmasprimnhnucsol= 0
nmasprimnhaitsol= 0
nmasprimnhaccsol= 0
nmasprimnhcorsol= 0
nmascondnnaccsol= 0
nmascondnncorsol= 0
nmasddepntaitsol= 0
nmasddepntaccsol= 0
nmasddepntcorsol= 0
nmasddepnhaitsol= 0
nmasddepnhaccsol= 0
nmasddepnhcorsol= 0
nmasddepnnaccsol= 0
nmasddepnncorsol= 0
nmasnuscntaitsol= 0
nmasnuscntaccsol= 0
nmasnuscntcorsol= 0
nmasnuscnhaitsol= 0
nmasnuscnhaccsol= 0
nmasnuscnhcorsol= 0
nmasnuscnnaccsol= 0
nmasnuscnncorsol= 0
nmasimscntaitsol= 0
nmasimscntaccsol= 0
nmasimscntcorsol= 0
nmasimscnhaitsol= 0
nmasimscnhaccsol= 0
nmasimscnhcorsol= 0
nmasimscnnaccsol= 0
nmasimscnncorsol= 0
nmascoagntintr23= 0
nmascoagnhintr23= 0
nmascoagntintr24= 0
nmascoagnhintr24= 0
nmascoagntintr34= 0
nmascoagnhintr34= 0
nmascoagnnintr34= 0
nmasmergntintr23= 0
nmasmergnhintr23= 0
nmasmergntintr34= 0
nmasmergnhintr34= 0
nmasmergnnintr34= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_sussbcoc_4mode

! ######################################################################
SUBROUTINE ukca_indices_sussbcocso_5mode

IMPLICIT NONE

!---------------------------------------------------------------
!
! Main array lengths and switches

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_SUSSBCOCSO_5MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ntraer=23          ! # of aerosol advected tracers
nbudaer=123        ! # of aerosol budget fields

! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  23    16    149     9 = 39+149+ 9=  197 (orgv1 ,traqu9 )
! ..  SUSSBCOC  23    16    149    38 = 39+149+38=  226 (orgv1 ,traqu38)
! ..  SUSSBCOC  23    76    149     9 = 99+149+ 9=  257 (orgv1c,traqu9 )
! ..  SUSSBCOC  23    76    149    38 = 99+149+38=  286 (orgv1c,traqu38)
!                        (123+26)
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 23    13  =  36 (orgv1 )
! ..  SUSSBCOC 23    52  =  75 (orgv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. below are 123 aerosol budget indices for SUSSBCOCSO [SO in SO]
!
nmasprimsuaitsol= 1
nmasprimsuaccsol= 2
nmasprimsucorsol= 3
nmasprimssaccsol= 4
nmasprimsscorsol= 5
nmasprimbcaitsol= 0 ! BC only emitted to insoluble
nmasprimbcaitins= 6
nmasprimocaitsol= 7
nmasprimocaitins= 8
!
nmasddepsunucsol= 9
nmasddepsuaitsol=10
nmasddepsuaccsol=11
nmasddepsucorsol=12
nmasddepssaccsol=13
nmasddepsscorsol=14
nmasddepbcaitsol=15
nmasddepbcaccsol=16
nmasddepbccorsol=17
nmasddepbcaitins=18
nmasddepocnucsol= 0  ! stored in NMASDDEPSONUCSOL
nmasddepocaitsol=19
nmasddepocaccsol=20
nmasddepoccorsol=21
nmasddepocaitins=22
nmasddepsonucsol=23
nmasddepsoaitsol=24
nmasddepsoaccsol=25
nmasddepsocorsol=26
!
nmasnuscsunucsol=27
nmasnuscsuaitsol=28
nmasnuscsuaccsol=29
nmasnuscsucorsol=30
nmasnuscssaccsol=31
nmasnuscsscorsol=32
nmasnuscbcaitsol=33
nmasnuscbcaccsol=34
nmasnuscbccorsol=35
nmasnuscbcaitins=36
nmasnuscocnucsol= 0  ! stored in NMASNUSCSONUCSOL
nmasnuscocaitsol=37
nmasnuscocaccsol=38
nmasnuscoccorsol=39
nmasnuscocaitins=40
nmasnuscsonucsol=41
nmasnuscsoaitsol=42
nmasnuscsoaccsol=43
nmasnuscsocorsol=44
!
nmasimscsunucsol=45
nmasimscsuaitsol=46
nmasimscsuaccsol=47
nmasimscsucorsol=48
nmasimscssaccsol=49
nmasimscsscorsol=50
nmasimscbcaitsol=51
nmasimscbcaccsol=52
nmasimscbccorsol=53
nmasimscbcaitins=54
nmasimscocnucsol= 0 ! stored in NMASIMSCSONUCSOL
nmasimscocaitsol=55
nmasimscocaccsol=56
nmasimscoccorsol=57
nmasimscocaitins=58
nmasimscsonucsol=59
nmasimscsoaitsol=60
nmasimscsoaccsol=61
nmasimscsocorsol=62
!
nmasclprsuaitsol1=63
nmasclprsuaccsol1=64
nmasclprsucorsol1=65
nmasclprsuaitsol2=66
nmasclprsuaccsol2=67
nmasclprsucorsol2=68
!
nmascondsunucsol=69
nmascondsuaitsol=70
nmascondsuaccsol=71
nmascondsucorsol=72
nmascondsuaitins=73
nmasnuclsunucsol=74
nmascondocnucsol= 0 ! stored in NMASCONDSONUCSOL
nmascondocaitsol= 0 ! stored in NMASCONDSOAITSOL
nmascondocaccsol= 0 ! stored in NMASCONDSOACCSOL
nmascondoccorsol= 0 ! stored in NMASCONDSOCORSOL
nmascondocaitins= 0 ! stored in NMASCONDSOAITINS
nmascondsonucsol=75
nmascondsoaitsol=76
nmascondsoaccsol=77
nmascondsocorsol=78
nmascondsoaitins=79
!
nmascoagsuintr12=80
nmascoagsuintr13=81
nmascoagsuintr14=82
nmascoagsuintr15=83
nmascoagocintr12= 0 ! stored in NMASCOAGSOINTR12
nmascoagocintr13= 0 ! stored in NMASCOAGSOINTR13
nmascoagocintr14= 0 ! stored in NMASCOAGSOINTR14
nmascoagocintr15= 0 ! stored in NMASCOAGSOINTR15
nmascoagsointr12=84
nmascoagsointr13=85
nmascoagsointr14=86
nmascoagsointr15=87
nmascoagsuintr23=88
nmascoagbcintr23=89
nmascoagocintr23=90
nmascoagsointr23=91
nmascoagsuintr24=92
nmascoagbcintr24=93
nmascoagocintr24=94
nmascoagsointr24=95
nmascoagsuintr34=96
nmascoagbcintr34=97
nmascoagocintr34=98
nmascoagssintr34=99
nmascoagsointr34=100
!
nmascoagbcintr53=101
nmascoagocintr53=102
nmascoagbcintr54=103
nmascoagocintr54=104
!
nmasagedsuintr52=105
nmasagedbcintr52=106
nmasagedocintr52=107
nmasagedsointr52=108
!
nmasmergsuintr12=109
nmasmergocintr12=  0 ! separate SO component
nmasmergsointr12=110
nmasmergsuintr23=111
nmasmergbcintr23=112
nmasmergocintr23=113
nmasmergsointr23=114
nmasmergsuintr34=115
nmasmergssintr34=116
nmasmergbcintr34=117
nmasmergocintr34=118
nmasmergsointr34=119
nmasprocsuintr23=120
nmasprocbcintr23=121
nmasprococintr23=122
nmasprocsointr23=123
!
! .. below are new ones for dust & modes 6/7 to be integrated
nmasprimduaccsol= 0 ! no DU in this setup
nmasprimducorsol= 0 ! no DU in this setup
nmasprimduaccins= 0 ! no DU in this setup
nmasprimducorins= 0 ! no DU in this setup
nmasddepduaccsol= 0 ! no DU in this setup
nmasddepducorsol= 0 ! no DU in this setup
nmasddepduaccins= 0 ! no DU in this setup
nmasddepducorins= 0 ! no DU in this setup
nmasnuscduaccsol= 0 ! no DU in this setup
nmasnuscducorsol= 0 ! no DU in this setup
nmasnuscduaccins= 0 ! no DU in this setup
nmasnuscducorins= 0 ! no DU in this setup
nmasimscduaccsol= 0 ! no DU in this setup
nmasimscducorsol= 0 ! no DU in this setup
nmasimscduaccins= 0 ! no DU in this setup
nmasimscducorins= 0 ! no DU in this setup
nmascondsuaccins= 0 ! no acc-ins cor-ins in this setup
nmascondsucorins= 0 ! no acc-ins cor-ins in this setup
nmascondocaccins= 0 ! no acc-ins cor-ins in this setup
nmascondoccorins= 0 ! no acc-ins cor-ins in this setup
nmascondsoaccins= 0 ! no acc-ins cor-ins in this setup
nmascondsocorins= 0 ! no acc-ins cor-ins in this setup
nmascoagsuintr16= 0 ! no acc-ins cor-ins in this setup
nmascoagsuintr17= 0 ! no acc-ins cor-ins in this setup
nmascoagocintr16= 0 ! no acc-ins cor-ins in this setup
nmascoagocintr17= 0 ! no acc-ins cor-ins in this setup
nmascoagsointr16= 0 ! no acc-ins cor-ins in this setup
nmascoagsointr17= 0 ! no acc-ins cor-ins in this setup
nmascoagduintr34= 0 ! no DU in this setup
nmascoagduintr64= 0 ! no DU in this setup
nmasagedsuintr63= 0 ! no acc-ins cor-ins in this setup
nmasagedduintr63= 0 ! no DU in this setup
nmasagedocintr63= 0 ! no BC/OC/SO in this setup
nmasagedsointr63= 0 ! no BC/OC/SO in this setup
nmasagedsuintr74= 0 ! no acc-ins cor-ins in this setup
nmasagedduintr74= 0 ! no DU in this setup
nmasagedocintr74= 0 ! no BC/OC/SO in this setup
nmasagedsointr74= 0 ! no BC/OC/SO in this setup
nmasmergduintr34= 0 ! no DU in this setup

! .. below are new ones for no3, nh4, and nano3 to be integrated
nmasprimntnucsol= 0
nmasprimntaitsol= 0
nmasprimntaccsol= 0
nmasprimntcorsol= 0
nmasprimnhnucsol= 0
nmasprimnhaitsol= 0
nmasprimnhaccsol= 0
nmasprimnhcorsol= 0
nmascondnnaccsol= 0
nmascondnncorsol= 0
nmasddepntaitsol= 0
nmasddepntaccsol= 0
nmasddepntcorsol= 0
nmasddepnhaitsol= 0
nmasddepnhaccsol= 0
nmasddepnhcorsol= 0
nmasddepnnaccsol= 0
nmasddepnncorsol= 0
nmasnuscntaitsol= 0
nmasnuscntaccsol= 0
nmasnuscntcorsol= 0
nmasnuscnhaitsol= 0
nmasnuscnhaccsol= 0
nmasnuscnhcorsol= 0
nmasnuscnnaccsol= 0
nmasnuscnncorsol= 0
nmasimscntaitsol= 0
nmasimscntaccsol= 0
nmasimscntcorsol= 0
nmasimscnhaitsol= 0
nmasimscnhaccsol= 0
nmasimscnhcorsol= 0
nmasimscnnaccsol= 0
nmasimscnncorsol= 0
nmascoagntintr23= 0
nmascoagnhintr23= 0
nmascoagntintr24= 0
nmascoagnhintr24= 0
nmascoagntintr34= 0
nmascoagnhintr34= 0
nmascoagnnintr34= 0
nmasmergntintr23= 0
nmasmergnhintr23= 0
nmasmergntintr34= 0
nmasmergnhintr34= 0
nmasmergnnintr34= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_sussbcocso_5mode

! ######################################################################
SUBROUTINE ukca_indices_sussbcocso_4mode

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_SUSSBCOCSO_4MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Main array lengths and switches

ntraer=20          ! # of aerosol advected tracers
nbudaer=104        ! # of aerosol budget fields

! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
! For orgv1 , NTRAG=16, NADVG=13
! For orgv1c, NTRAG=76, NADVG=52
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU         = NVTOT
! ..  SUSSBCOC  20    16    130     9 = 36+130+ 9=  175 (orgv1 ,traqu9 )
! ..  SUSSBCOC  20    16    130    38 = 36+130+38=  204 (orgv1 ,traqu38)
! ..  SUSSBCOC  20    76    130     9 = 36+190+ 9=  235 (orgv1c,traqu9 )
! ..  SUSSBCOC  20    76    130    38 = 36+190+38=  264 (orgv1c,traqu38)
!                         (104+26)
!            NTRAER+NADVG  NTRA
! ..  SUSSBCOC 20    13  =  33 (orgv1 )
! ..  SUSSBCOC 20    52  =  72 (orgv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. 104 aerosol budget indices for SUSSBCOCSO_4MODE [SO in SO]
!
nmasprimsuaitsol= 1
nmasprimsuaccsol= 2
nmasprimsucorsol= 3
nmasprimssaccsol= 4
nmasprimsscorsol= 5
nmasprimbcaitsol= 6
nmasprimbcaitins= 0 ! BC only emitted to soluble
nmasprimocaitsol= 7
nmasprimocaitins= 0 ! OC only emitted to soluble
!
nmasddepsunucsol= 8
nmasddepsuaitsol= 9
nmasddepsuaccsol=10
nmasddepsucorsol=11
nmasddepssaccsol=12
nmasddepsscorsol=13
nmasddepbcaitsol=14
nmasddepbcaccsol=15
nmasddepbccorsol=16
nmasddepbcaitins= 0  ! only soluble modes
nmasddepocnucsol= 0  ! stored in NMASDDEPSONUCSOL
nmasddepocaitsol=17
nmasddepocaccsol=18
nmasddepoccorsol=19
nmasddepocaitins= 0  ! only soluble modes
nmasddepsonucsol=20
nmasddepsoaitsol=21
nmasddepsoaccsol=22
nmasddepsocorsol=23
!
nmasnuscsunucsol=24
nmasnuscsuaitsol=25
nmasnuscsuaccsol=26
nmasnuscsucorsol=27
nmasnuscssaccsol=28
nmasnuscsscorsol=29
nmasnuscbcaitsol=30
nmasnuscbcaccsol=31
nmasnuscbccorsol=32
nmasnuscbcaitins= 0 ! only soluble modes
nmasnuscocnucsol= 0 ! stored in NMASNUSCSONUCSOL
nmasnuscocaitsol=33
nmasnuscocaccsol=34
nmasnuscoccorsol=35
nmasnuscocaitins= 0 ! only soluble modes
nmasnuscsonucsol=36
nmasnuscsoaitsol=37
nmasnuscsoaccsol=38
nmasnuscsocorsol=39
!
nmasimscsunucsol=40
nmasimscsuaitsol=41
nmasimscsuaccsol=42
nmasimscsucorsol=43
nmasimscssaccsol=44
nmasimscsscorsol=45
nmasimscbcaitsol=46
nmasimscbcaccsol=47
nmasimscbccorsol=48
nmasimscbcaitins= 0 ! only soluble modes
nmasimscocnucsol= 0 ! stored in NMASIMSCSONUCSOL
nmasimscocaitsol=49
nmasimscocaccsol=50
nmasimscoccorsol=51
nmasimscocaitins= 0 ! only soluble modes
nmasimscsonucsol=52
nmasimscsoaitsol=53
nmasimscsoaccsol=54
nmasimscsocorsol=55
!
nmasclprsuaitsol1=56
nmasclprsuaccsol1=57
nmasclprsucorsol1=58
nmasclprsuaitsol2=59
nmasclprsuaccsol2=60
nmasclprsucorsol2=61
!
nmascondsunucsol=62
nmascondsuaitsol=63
nmascondsuaccsol=64
nmascondsucorsol=65
nmascondsuaitins= 0 ! only soluble modes
nmasnuclsunucsol=66
nmascondocnucsol= 0 ! stored in NMASCONDSONUCSOL
nmascondocaitsol= 0 ! stored in NMASCONDSOAITSOL
nmascondocaccsol= 0 ! stored in NMASCONDSOACCSOL
nmascondoccorsol= 0 ! stored in NMASCONDSOCORSOL
nmascondocaitins= 0 ! stored in NMASCONDSOAITINS
nmascondsonucsol=67
nmascondsoaitsol=68
nmascondsoaccsol=69
nmascondsocorsol=70
nmascondsoaitins= 0 ! only soluble modes
!
nmascoagsuintr12=71
nmascoagsuintr13=72
nmascoagsuintr14=73
nmascoagsuintr15= 0 ! only soluble modes
nmascoagocintr12= 0 ! stored in NMASCOAGSOINTR12
nmascoagocintr13= 0 ! stored in NMASCOAGSOINTR13
nmascoagocintr14= 0 ! stored in NMASCOAGSOINTR14
nmascoagocintr15= 0 ! stored in NMASCOAGSOINTR15
nmascoagsointr12=74
nmascoagsointr13=75
nmascoagsointr14=76
nmascoagsointr15= 0 ! only soluble modes
nmascoagsuintr23=77
nmascoagbcintr23=78
nmascoagocintr23=79
nmascoagsointr23=80
nmascoagsuintr24=81
nmascoagbcintr24=82
nmascoagocintr24=83
nmascoagsointr24=84
nmascoagsuintr34=85
nmascoagbcintr34=86
nmascoagocintr34=87
nmascoagssintr34=88
nmascoagsointr34=89
!
nmascoagbcintr53= 0 ! only soluble modes
nmascoagocintr53= 0 ! only soluble modes
nmascoagbcintr54= 0 ! only soluble modes
nmascoagocintr54= 0 ! only soluble modes
!
nmasagedsuintr52= 0 ! only soluble modes
nmasagedbcintr52= 0 ! only soluble modes
nmasagedocintr52= 0 ! only soluble modes
nmasagedsointr52= 0 ! only soluble modes
!
nmasmergsuintr12=90
nmasmergocintr12= 0 ! separate SO component
nmasmergsointr12=91
nmasmergsuintr23=92
nmasmergbcintr23=93
nmasmergocintr23=94
nmasmergsointr23=95
nmasmergsuintr34=96
nmasmergssintr34=97
nmasmergbcintr34=98
nmasmergocintr34=99
nmasmergsointr34=100
nmasprocsuintr23=101
nmasprocbcintr23=102
nmasprococintr23=103
nmasprocsointr23=104
!
! .. below are new ones for dust & modes 6/7 to be integrated
nmasprimduaccsol= 0 ! no DU in this setup
nmasprimducorsol= 0 ! no DU in this setup
nmasprimduaccins= 0 ! no DU in this setup
nmasprimducorins= 0 ! no DU in this setup
nmasddepduaccsol= 0 ! no DU in this setup
nmasddepducorsol= 0 ! no DU in this setup
nmasddepduaccins= 0 ! no DU in this setup
nmasddepducorins= 0 ! no DU in this setup
nmasnuscduaccsol= 0 ! no DU in this setup
nmasnuscducorsol= 0 ! no DU in this setup
nmasnuscduaccins= 0 ! no DU in this setup
nmasnuscducorins= 0 ! no DU in this setup
nmasimscduaccsol= 0 ! no DU in this setup
nmasimscducorsol= 0 ! no DU in this setup
nmasimscduaccins= 0 ! no DU in this setup
nmasimscducorins= 0 ! no DU in this setup
nmascondsuaccins= 0 ! no acc-ins cor-ins in this setup
nmascondsucorins= 0 ! no acc-ins cor-ins in this setup
nmascondocaccins= 0 ! no acc-ins cor-ins in this setup
nmascondoccorins= 0 ! no acc-ins cor-ins in this setup
nmascondsoaccins= 0 ! no acc-ins cor-ins in this setup
nmascondsocorins= 0 ! no acc-ins cor-ins in this setup
nmascoagsuintr16= 0 ! no acc-ins cor-ins in this setup
nmascoagsuintr17= 0 ! no acc-ins cor-ins in this setup
nmascoagocintr16= 0 ! no acc-ins cor-ins in this setup
nmascoagocintr17= 0 ! no acc-ins cor-ins in this setup
nmascoagsointr16= 0 ! no acc-ins cor-ins in this setup
nmascoagsointr17= 0 ! no acc-ins cor-ins in this setup
nmascoagduintr34= 0 ! no DU in this setup
nmascoagduintr64= 0 ! no DU in this setup
nmasagedsuintr63= 0 ! no acc-ins cor-ins in this setup
nmasagedduintr63= 0 ! no DU in this setup
nmasagedocintr63= 0 ! no BC/OC/SO in this setup
nmasagedsointr63= 0 ! no BC/OC/SO in this setup
nmasagedsuintr74= 0 ! no acc-ins cor-ins in this setup
nmasagedduintr74= 0 ! no DU in this setup
nmasagedocintr74= 0 ! no BC/OC/SO in this setup
nmasagedsointr74= 0 ! no BC/OC/SO in this setup
nmasmergduintr34= 0 ! no DU in this setup

! .. below are new ones for no3, nh4, and nano3 to be integrated
nmasprimntnucsol= 0
nmasprimntaitsol= 0
nmasprimntaccsol= 0
nmasprimntcorsol= 0
nmasprimnhnucsol= 0
nmasprimnhaitsol= 0
nmasprimnhaccsol= 0
nmasprimnhcorsol= 0
nmascondnnaccsol= 0
nmascondnncorsol= 0
nmasddepntaitsol= 0
nmasddepntaccsol= 0
nmasddepntcorsol= 0
nmasddepnhaitsol= 0
nmasddepnhaccsol= 0
nmasddepnhcorsol= 0
nmasddepnnaccsol= 0
nmasddepnncorsol= 0
nmasnuscntaitsol= 0
nmasnuscntaccsol= 0
nmasnuscntcorsol= 0
nmasnuscnhaitsol= 0
nmasnuscnhaccsol= 0
nmasnuscnhcorsol= 0
nmasnuscnnaccsol= 0
nmasnuscnncorsol= 0
nmasimscntaitsol= 0
nmasimscntaccsol= 0
nmasimscntcorsol= 0
nmasimscnhaitsol= 0
nmasimscnhaccsol= 0
nmasimscnhcorsol= 0
nmasimscnnaccsol= 0
nmasimscnncorsol= 0
nmascoagntintr23= 0
nmascoagnhintr23= 0
nmascoagntintr24= 0
nmascoagnhintr24= 0
nmascoagntintr34= 0
nmascoagnhintr34= 0
nmascoagnnintr34= 0
nmasmergntintr23= 0
nmasmergnhintr23= 0
nmasmergntintr34= 0
nmasmergnhintr34= 0
nmasmergnnintr34= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_sussbcocso_4mode

! ######################################################################
SUBROUTINE ukca_indices_suss_4mode

IMPLICIT NONE


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_SUSS_4MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Main array lengths and switches
ntraer=10          ! # of aerosol advected tracers
nbudaer=46         ! # of aerosol budget fields

!
! For sv1         : NTRAG=13, NADVG=11
! For sv1_coupled : NTRAG=74, NADVG=50
!
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU             = NVTOT
! ..  SUSS      10    13     74     9 = 23 +  74 + 9 = 106(sv1 ,traqu9 )
! ..  SUSS      10    13     74     38= 23 +  74 + 38= 135(sv1 ,traqu38)
! ..  SUSS      10    74     74     9 = 84 +  74 + 9 = 167(sv1c,traqu9 )
! ..  SUSS      10    74     74     38= 84 +  74 + 38= 196(sv1c,traqu38)
!
!            NTRAER+NADVG  NTRA
! ..  SUSS      10   11  =  21 (sv1 )
! ..  SUSS      10   50  =  60 (sv1c)
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. below are 46 aerosol budget variables for SUSS aerosol system
!
nmasprimsuaitsol= 1
nmasprimsuaccsol= 2
nmasprimsucorsol= 3
nmasprimssaccsol= 4
nmasprimsscorsol= 5
nmasprimbcaitsol= 0 ! no BC/OC/SO in this setup
nmasprimbcaitins= 0 ! no BC/OC/SO in this setup
nmasprimocaitsol= 0 ! no BC/OC/SO in this setup
nmasprimocaitins= 0 ! no BC/OC/SO in this setup
!
nmasddepsunucsol= 6
nmasddepsuaitsol= 7
nmasddepsuaccsol= 8
nmasddepsucorsol= 9
nmasddepssaccsol=10
nmasddepsscorsol=11
nmasddepbcaitsol= 0 ! no BC/OC/SO in this setup
nmasddepbcaccsol= 0 ! no BC/OC/SO in this setup
nmasddepbccorsol= 0 ! no BC/OC/SO in this setup
nmasddepbcaitins= 0 ! no BC/OC/SO in this setup
nmasddepocnucsol= 0 ! no BC/OC/SO in this setup
nmasddepocaitsol= 0 ! no BC/OC/SO in this setup
nmasddepocaccsol= 0 ! no BC/OC/SO in this setup
nmasddepoccorsol= 0 ! no BC/OC/SO in this setup
nmasddepocaitins= 0 ! no BC/OC/SO in this setup
nmasddepsonucsol= 0 ! no BC/OC/SO in this setup
nmasddepsoaitsol= 0 ! no BC/OC/SO in this setup
nmasddepsoaccsol= 0 ! no BC/OC/SO in this setup
nmasddepsocorsol= 0 ! no BC/OC/SO in this setup
!
nmasnuscsunucsol=12
nmasnuscsuaitsol=13
nmasnuscsuaccsol=14
nmasnuscsucorsol=15
nmasnuscssaccsol=16
nmasnuscsscorsol=17
nmasnuscbcaitsol= 0
nmasnuscbcaccsol= 0 ! no BC/OC/SO in this setup
nmasnuscbccorsol= 0 ! no BC/OC/SO in this setup
nmasnuscbcaitins= 0 ! no BC/OC/SO in this setup
nmasnuscocnucsol= 0 ! no BC/OC/SO in this setup
nmasnuscocaitsol= 0 ! no BC/OC/SO in this setup
nmasnuscocaccsol= 0 ! no BC/OC/SO in this setup
nmasnuscoccorsol= 0 ! no BC/OC/SO in this setup
nmasnuscocaitins= 0 ! no BC/OC/SO in this setup
nmasnuscsonucsol= 0 ! no BC/OC/SO in this setup
nmasnuscsoaitsol= 0 ! no BC/OC/SO in this setup
nmasnuscsoaccsol= 0 ! no BC/OC/SO in this setup
nmasnuscsocorsol= 0 ! no BC/OC/SO in this setup
!
nmasimscsunucsol=18
nmasimscsuaitsol=19
nmasimscsuaccsol=20
nmasimscsucorsol=21
nmasimscssaccsol=22
nmasimscsscorsol=23
nmasimscbcaitsol= 0 ! no BC/OC/SO in this setup
nmasimscbcaccsol= 0 ! no BC/OC/SO in this setup
nmasimscbccorsol= 0 ! no BC/OC/SO in this setup
nmasimscbcaitins= 0 ! no BC/OC/SO in this setup
nmasimscocnucsol= 0 ! no BC/OC/SO in this setup
nmasimscocaitsol= 0 ! no BC/OC/SO in this setup
nmasimscocaccsol= 0 ! no BC/OC/SO in this setup
nmasimscoccorsol= 0 ! no BC/OC/SO in this setup
nmasimscocaitins= 0 ! no BC/OC/SO in this setup
nmasimscsonucsol= 0 ! no BC/OC/SO in this setup
nmasimscsoaitsol= 0 ! no BC/OC/SO in this setup
nmasimscsoaccsol= 0 ! no BC/OC/SO in this setup
nmasimscsocorsol= 0 ! no BC/OC/SO in this setup
!
nmasclprsuaitsol1=24
nmasclprsuaccsol1=25
nmasclprsucorsol1=26
nmasclprsuaitsol2=27
nmasclprsuaccsol2=28
nmasclprsucorsol2=29
!
nmascondsunucsol=30
nmascondsuaitsol=31
nmascondsuaccsol=32
nmascondsucorsol=33
nmascondsuaitins= 0 ! only soluble modes in this setup
nmasnuclsunucsol=34
nmascondocnucsol= 0 ! no BC/OC/SO in this setup
nmascondocaitsol= 0 ! no BC/OC/SO in this setup
nmascondocaccsol= 0 ! no BC/OC/SO in this setup
nmascondoccorsol= 0 ! no BC/OC/SO in this setup
nmascondocaitins= 0 ! no BC/OC/SO in this setup
nmascondsonucsol= 0 ! no BC/OC/SO in this setup
nmascondsoaitsol= 0 ! no BC/OC/SO in this setup
nmascondsoaccsol= 0 ! no BC/OC/SO in this setup
nmascondsocorsol= 0 ! no BC/OC/SO in this setup
nmascondsoaitins= 0 ! no BC/OC/SO in this setup
!
nmascoagsuintr12=35
nmascoagsuintr13=36
nmascoagsuintr14=37
nmascoagsuintr15= 0 ! only soluble modes in this setup
nmascoagocintr12= 0 ! no BC/OC/SO in this setup
nmascoagocintr13= 0 ! no BC/OC/SO in this setup
nmascoagocintr14= 0 ! no BC/OC/SO in this setup
nmascoagocintr15= 0 ! no BC/OC/SO in this setup
nmascoagsointr12= 0 ! no BC/OC/SO in this setup
nmascoagsointr13= 0 ! no BC/OC/SO in this setup
nmascoagsointr14= 0 ! no BC/OC/SO in this setup
nmascoagsointr15= 0 ! no BC/OC/SO in this setup
nmascoagsuintr23=38
nmascoagbcintr23= 0 ! no BC/OC/SO in this setup
nmascoagocintr23= 0 ! no BC/OC/SO in this setup
nmascoagsointr23= 0 ! no BC/OC/SO in this setup
nmascoagsuintr24=39
nmascoagbcintr24= 0 ! no BC/OC/SO in this setup
nmascoagocintr24= 0 ! no BC/OC/SO in this setup
nmascoagsointr24= 0 ! no BC/OC/SO in this setup
nmascoagsuintr34=40
nmascoagbcintr34= 0 ! no BC/OC/SO in this setup
nmascoagocintr34= 0 ! no BC/OC/SO in this setup
nmascoagssintr34=41
nmascoagsointr34= 0 ! no BC/OC/SO in this setup
!
nmascoagbcintr53= 0 ! no BC/OC/SO in this setup
nmascoagocintr53= 0 ! no BC/OC/SO in this setup
nmascoagbcintr54= 0 ! no BC/OC/SO in this setup
nmascoagocintr54= 0 ! no BC/OC/SO in this setup
!
nmasagedsuintr52= 0 ! only soluble modes in this setup
nmasagedbcintr52= 0 ! no BC/OC/SO in this setup
nmasagedocintr52= 0 ! no BC/OC/SO in this setup
nmasagedsointr52= 0 ! no BC/OC/SO in this setup
!
nmasmergsuintr12=42
nmasmergocintr12= 0 ! no BC/OC/SO in this setup
nmasmergsointr12= 0 ! no BC/OC/SO in this setup
nmasmergsuintr23=43
nmasmergbcintr23= 0 ! no BC/OC/SO in this setup
nmasmergocintr23= 0 ! no BC/OC/SO in this setup
nmasmergsointr23= 0 ! no BC/OC/SO in this setup
nmasmergsuintr34=44
nmasmergssintr34=45
nmasmergbcintr34= 0 ! no BC/OC/SO in this setup
nmasmergocintr34= 0 ! no BC/OC/SO in this setup
nmasmergsointr34= 0 ! no BC/OC/SO in this setup
nmasprocsuintr23=46
nmasprocbcintr23= 0 ! no BC/OC/SO in this setup
nmasprococintr23= 0 ! no BC/OC/SO in this setup
nmasprocsointr23= 0 ! no BC/OC/SO in this setup
!
! .. below are new ones for dust & modes 6/7 to be integrated
nmasprimduaccsol= 0 ! no DU in this setup
nmasprimducorsol= 0 ! no DU in this setup
nmasprimduaccins= 0 ! no DU in this setup
nmasprimducorins= 0 ! no DU in this setup
nmasddepduaccsol= 0 ! no DU in this setup
nmasddepducorsol= 0 ! no DU in this setup
nmasddepduaccins= 0 ! no DU in this setup
nmasddepducorins= 0 ! no DU in this setup
nmasnuscduaccsol= 0 ! no DU in this setup
nmasnuscducorsol= 0 ! no DU in this setup
nmasnuscduaccins= 0 ! no DU in this setup
nmasnuscducorins= 0 ! no DU in this setup
nmasimscduaccsol= 0 ! no DU in this setup
nmasimscducorsol= 0 ! no DU in this setup
nmasimscduaccins= 0 ! no DU in this setup
nmasimscducorins= 0 ! no DU in this setup
nmascondsuaccins= 0 ! only soluble modes in this setup
nmascondsucorins= 0 ! only soluble modes in this setup
nmascondocaccins= 0 ! only soluble modes in this setup
nmascondoccorins= 0 ! only soluble modes in this setup
nmascondsoaccins= 0 ! only soluble modes in this setup
nmascondsocorins= 0 ! only soluble modes in this setup
nmascoagsuintr16= 0 ! only soluble modes in this setup
nmascoagsuintr17= 0 ! only soluble modes in this setup
nmascoagocintr16= 0 ! no BC/OC/SO in this setup
nmascoagocintr17= 0 ! no BC/OC/SO in this setup
nmascoagsointr16= 0 ! no BC/OC/SO in this setup
nmascoagsointr17= 0 ! no BC/OC/SO in this setup
nmascoagduintr34= 0 ! no DU in this setup
nmascoagduintr64= 0 ! no DU in this setup
nmasagedsuintr63= 0 ! only soluble modes in this setup
nmasagedduintr63= 0 ! no DU in this setup
nmasagedocintr63= 0 ! no BC/OC/SO in this setup
nmasagedsointr63= 0 ! no BC/OC/SO in this setup
nmasagedsuintr74= 0 ! only soluble modes in this setup
nmasagedduintr74= 0 ! no DU in this setup
nmasagedocintr74= 0 ! no BC/OC/SO in this setup
nmasagedsointr74= 0 ! no BC/OC/SO in this setup
nmasmergduintr34= 0 ! no DU in this setup

! .. below are new ones for no3, nh4, and nano3 to be integrated
nmasprimntnucsol= 0
nmasprimntaitsol= 0
nmasprimntaccsol= 0
nmasprimntcorsol= 0
nmasprimnhnucsol= 0
nmasprimnhaitsol= 0
nmasprimnhaccsol= 0
nmasprimnhcorsol= 0
nmascondnnaccsol= 0
nmascondnncorsol= 0
nmasddepntaitsol= 0
nmasddepntaccsol= 0
nmasddepntcorsol= 0
nmasddepnhaitsol= 0
nmasddepnhaccsol= 0
nmasddepnhcorsol= 0
nmasddepnnaccsol= 0
nmasddepnncorsol= 0
nmasnuscntaitsol= 0
nmasnuscntaccsol= 0
nmasnuscntcorsol= 0
nmasnuscnhaitsol= 0
nmasnuscnhaccsol= 0
nmasnuscnhcorsol= 0
nmasnuscnnaccsol= 0
nmasnuscnncorsol= 0
nmasimscntaitsol= 0
nmasimscntaccsol= 0
nmasimscntcorsol= 0
nmasimscnhaitsol= 0
nmasimscnhaccsol= 0
nmasimscnhcorsol= 0
nmasimscnnaccsol= 0
nmasimscnncorsol= 0
nmascoagntintr23= 0
nmascoagnhintr23= 0
nmascoagntintr24= 0
nmascoagnhintr24= 0
nmascoagntintr34= 0
nmascoagnhintr34= 0
nmascoagnnintr34= 0
nmasmergntintr23= 0
nmasmergnhintr23= 0
nmasmergntintr34= 0
nmasmergnhintr34= 0
nmasmergnnintr34= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_suss_4mode

! ######################################################################
SUBROUTINE UKCA_INDICES_DUonly_2MODE

IMPLICIT NONE


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_DUONLY_2MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Main array lengths and switches
ntraer=4           ! # of aerosol advected tracers
nbudaer=8          ! # of aerosol budget fields
!
! For nochem: NTRAG= 2, NADVG= 2
!
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
!            NTRAER+NTRAG+NBUDGET+NTRAQU             = NVTOT
! ..  DU-only    4     2      8     9 =  6 +   8 + 9 =  23 (nochem)
! ..  DU-only    4     2      8     38=  6 +   8 + 38=  52 (nochem)
!
!            NTRAER+NADVG  NTRA
! ..  DU-only    4    2  =   6
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. below are  8 aerosol budget variables for DU-only aerosol system
!
nmasprimsuaitsol= 0 ! no SO4 or SS in this setup
nmasprimsuaccsol= 0 ! no SO4 or SS in this setup
nmasprimsucorsol= 0 ! no SO4 or SS in this setup
nmasprimssaccsol= 0 ! no SO4 or SS in this setup
nmasprimsscorsol= 0 ! no SO4 or SS in this setup
nmasprimbcaitsol= 0 ! no BC/OC/SO in this setup
nmasprimbcaitins= 0 ! no BC/OC/SO in this setup
nmasprimocaitsol= 0 ! no BC/OC/SO in this setup
nmasprimocaitins= 0 ! no BC/OC/SO in this setup
!
nmasddepsunucsol= 0 ! no SO4 or SS in this setup
nmasddepsuaitsol= 0 ! no SO4 or SS in this setup
nmasddepsuaccsol= 0 ! no SO4 or SS in this setup
nmasddepsucorsol= 0 ! no SO4 or SS in this setup
nmasddepssaccsol= 0 ! no SO4 or SS in this setup
nmasddepsscorsol= 0 ! no SO4 or SS in this setup
nmasddepbcaitsol= 0 ! no BC/OC/SO in this setup
nmasddepbcaccsol= 0 ! no BC/OC/SO in this setup
nmasddepbccorsol= 0 ! no BC/OC/SO in this setup
nmasddepbcaitins= 0 ! no BC/OC/SO in this setup
nmasddepocnucsol= 0 ! no BC/OC/SO in this setup
nmasddepocaitsol= 0 ! no BC/OC/SO in this setup
nmasddepocaccsol= 0 ! no BC/OC/SO in this setup
nmasddepoccorsol= 0 ! no BC/OC/SO in this setup
nmasddepocaitins= 0 ! no BC/OC/SO in this setup
nmasddepsonucsol= 0 ! no BC/OC/SO in this setup
nmasddepsoaitsol= 0 ! no BC/OC/SO in this setup
nmasddepsoaccsol= 0 ! no BC/OC/SO in this setup
nmasddepsocorsol= 0 ! no BC/OC/SO in this setup
!
nmasnuscsunucsol= 0 ! no SO4 or SS in this setup
nmasnuscsuaitsol= 0 ! no SO4 or SS in this setup
nmasnuscsuaccsol= 0 ! no SO4 or SS in this setup
nmasnuscsucorsol= 0 ! no SO4 or SS in this setup
nmasnuscssaccsol= 0 ! no SO4 or SS in this setup
nmasnuscsscorsol= 0 ! no SO4 or SS in this setup
nmasnuscbcaitsol= 0 ! no BC/OC/SO in this setup
nmasnuscbcaccsol= 0 ! no BC/OC/SO in this setup
nmasnuscbccorsol= 0 ! no BC/OC/SO in this setup
nmasnuscbcaitins= 0 ! no BC/OC/SO in this setup
nmasnuscocnucsol= 0 ! no BC/OC/SO in this setup
nmasnuscocaitsol= 0 ! no BC/OC/SO in this setup
nmasnuscocaccsol= 0 ! no BC/OC/SO in this setup
nmasnuscoccorsol= 0 ! no BC/OC/SO in this setup
nmasnuscocaitins= 0 ! no BC/OC/SO in this setup
nmasnuscsonucsol= 0 ! no BC/OC/SO in this setup
nmasnuscsoaitsol= 0 ! no BC/OC/SO in this setup
nmasnuscsoaccsol= 0 ! no BC/OC/SO in this setup
nmasnuscsocorsol= 0 ! no BC/OC/SO in this setup
!
nmasimscsunucsol= 0 ! no SO4 or SS in this setup
nmasimscsuaitsol= 0 ! no SO4 or SS in this setup
nmasimscsuaccsol= 0 ! no SO4 or SS in this setup
nmasimscsucorsol= 0 ! no SO4 or SS in this setup
nmasimscssaccsol= 0 ! no SO4 or SS in this setup
nmasimscsscorsol= 0 ! no SO4 or SS in this setup
nmasimscbcaitsol= 0 ! no BC/OC/SO in this setup
nmasimscbcaccsol= 0 ! no BC/OC/SO in this setup
nmasimscbccorsol= 0 ! no BC/OC/SO in this setup
nmasimscbcaitins= 0 ! no BC/OC/SO in this setup
nmasimscocnucsol= 0 ! no BC/OC/SO in this setup
nmasimscocaitsol= 0 ! no BC/OC/SO in this setup
nmasimscocaccsol= 0 ! no BC/OC/SO in this setup
nmasimscoccorsol= 0 ! no BC/OC/SO in this setup
nmasimscocaitins= 0 ! no BC/OC/SO in this setup
nmasimscsonucsol= 0 ! no BC/OC/SO in this setup
nmasimscsoaitsol= 0 ! no BC/OC/SO in this setup
nmasimscsoaccsol= 0 ! no BC/OC/SO in this setup
nmasimscsocorsol= 0 ! no BC/OC/SO in this setup
!
nmasclprsuaitsol1=0 ! no SO4 or SS in this setup
nmasclprsuaccsol1=0 ! no SO4 or SS in this setup
nmasclprsucorsol1=0 ! no SO4 or SS in this setup
nmasclprsuaitsol2=0 ! no SO4 or SS in this setup
nmasclprsuaccsol2=0 ! no SO4 or SS in this setup
nmasclprsucorsol2=0 ! no SO4 or SS in this setup
!
nmascondsunucsol= 0 ! no SO4 or SS in this setup
nmascondsuaitsol= 0 ! no SO4 or SS in this setup
nmascondsuaccsol= 0 ! no SO4 or SS in this setup
nmascondsucorsol= 0 ! no SO4 or SS in this setup
nmascondsuaitins= 0 ! only soluble modes in this setup
nmasnuclsunucsol= 0 ! no SO4 or SS in this setup
nmascondocnucsol= 0 ! no BC/OC/SO in this setup
nmascondocaitsol= 0 ! no BC/OC/SO in this setup
nmascondocaccsol= 0 ! no BC/OC/SO in this setup
nmascondoccorsol= 0 ! no BC/OC/SO in this setup
nmascondocaitins= 0 ! no BC/OC/SO in this setup
nmascondsonucsol= 0 ! no BC/OC/SO in this setup
nmascondsoaitsol= 0 ! no BC/OC/SO in this setup
nmascondsoaccsol= 0 ! no BC/OC/SO in this setup
nmascondsocorsol= 0 ! no BC/OC/SO in this setup
nmascondsoaitins= 0 ! no BC/OC/SO in this setup
!
nmascoagsuintr12= 0 ! no SO4 or SS in this setup
nmascoagsuintr13= 0 ! no SO4 or SS in this setup
nmascoagsuintr14= 0 ! no SO4 or SS in this setup
nmascoagsuintr15= 0 ! only soluble modes in this setup
nmascoagocintr12= 0 ! no BC/OC/SO in this setup
nmascoagocintr13= 0 ! no BC/OC/SO in this setup
nmascoagocintr14= 0 ! no BC/OC/SO in this setup
nmascoagocintr15= 0 ! no BC/OC/SO in this setup
nmascoagsointr12= 0 ! no BC/OC/SO in this setup
nmascoagsointr13= 0 ! no BC/OC/SO in this setup
nmascoagsointr14= 0 ! no BC/OC/SO in this setup
nmascoagsointr15= 0 ! no BC/OC/SO in this setup
nmascoagsuintr23= 0 ! no SO4 or SS in this setup
nmascoagbcintr23= 0 ! no BC/OC/SO in this setup
nmascoagocintr23= 0 ! no BC/OC/SO in this setup
nmascoagsointr23= 0 ! no BC/OC/SO in this setup
nmascoagsuintr24= 0 ! no SO4 or SS in this setup
nmascoagbcintr24= 0 ! no BC/OC/SO in this setup
nmascoagocintr24= 0 ! no BC/OC/SO in this setup
nmascoagsointr24= 0 ! no BC/OC/SO in this setup
nmascoagsuintr34= 0 ! no SO4 or SS in this setup
nmascoagbcintr34= 0 ! no BC/OC/SO in this setup
nmascoagocintr34= 0 ! no BC/OC/SO in this setup
nmascoagssintr34= 0 ! no SO4 or SS in this setup
nmascoagsointr34= 0 ! no BC/OC/SO in this setup
!
nmascoagbcintr53= 0 ! no BC/OC/SO in this setup
nmascoagocintr53= 0 ! no BC/OC/SO in this setup
nmascoagbcintr54= 0 ! no BC/OC/SO in this setup
nmascoagocintr54= 0 ! no BC/OC/SO in this setup
!
nmasagedsuintr52= 0 ! only soluble modes in this setup
nmasagedbcintr52= 0 ! no BC/OC/SO in this setup
nmasagedocintr52= 0 ! no BC/OC/SO in this setup
nmasagedsointr52= 0 ! no BC/OC/SO in this setup
!
nmasmergsuintr12= 0 ! no SO4 or SS in this setup
nmasmergocintr12= 0 ! no BC/OC/SO in this setup
nmasmergsointr12= 0 ! no BC/OC/SO in this setup
nmasmergsuintr23= 0 ! no SO4 or SS in this setup
nmasmergbcintr23= 0 ! no BC/OC/SO in this setup
nmasmergocintr23= 0 ! no BC/OC/SO in this setup
nmasmergsointr23= 0 ! no BC/OC/SO in this setup
nmasmergsuintr34= 0 ! no SO4 or SS in this setup
nmasmergssintr34= 0 ! no SO4 or SS in this setup
nmasmergbcintr34= 0 ! no BC/OC/SO in this setup
nmasmergocintr34= 0 ! no BC/OC/SO in this setup
nmasmergsointr34= 0 ! no BC/OC/SO in this setup
nmasprocsuintr23= 0 ! no SO4 or SS in this setup
nmasprocbcintr23= 0 ! no BC/OC/SO in this setup
nmasprococintr23= 0 ! no BC/OC/SO in this setup
nmasprocsointr23= 0 ! no BC/OC/SO in this setup
!
! .. below are new ones for dust & modes 6/7 to be integrated
nmasprimduaccsol= 0 ! DU emitted into insoluble modes
nmasprimducorsol= 0 ! DU emitted into insoluble modes
nmasprimduaccins= 1
nmasprimducorins= 2
nmasddepduaccsol= 0 ! no aged DU in this setup
nmasddepducorsol= 0 ! no aged DU in this setup
nmasddepduaccins= 3
nmasddepducorins= 4
nmasnuscduaccsol= 0 ! no aged DU in this setup
nmasnuscducorsol= 0 ! no aged DU in this setup
nmasnuscduaccins= 5
nmasnuscducorins= 6
nmasimscduaccsol= 0 ! no aged DU in this setup
nmasimscducorsol= 0 ! no aged DU in this setup
nmasimscduaccins= 7
nmasimscducorins= 8
nmascondsuaccins= 0 ! no SO4 or SS in this setup
nmascondsucorins= 0 ! no SO4 or SS in this setup
nmascondocaccins= 0 ! no BC/OC/SO in this setup
nmascondoccorins= 0 ! no BC/OC/SO in this setup
nmascondsoaccins= 0 ! no BC/OC/SO in this setup
nmascondsocorins= 0 ! no BC/OC/SO in this setup
nmascoagsuintr16= 0 ! no SO4 or SS in this setup
nmascoagsuintr17= 0 ! no SO4 or SS in this setup
nmascoagocintr16= 0 ! no BC/OC/SO in this setup
nmascoagocintr17= 0 ! no BC/OC/SO in this setup
nmascoagsointr16= 0 ! no BC/OC/SO in this setup
nmascoagsointr17= 0 ! no BC/OC/SO in this setup
nmascoagduintr34= 0 ! no aged DU in this setup
nmascoagduintr64= 0 ! no aged DU in this setup
nmasagedsuintr63= 0 ! no aged DU in this setup
nmasagedduintr63= 0 ! no aged DU in this setup
nmasagedocintr63= 0 ! no aged DU in this setup
nmasagedsointr63= 0 ! no aged DU in this setup
nmasagedsuintr74= 0 ! no aged DU in this setup
nmasagedduintr74= 0 ! no aged DU in this setup
nmasagedocintr74= 0 ! no aged DU in this setup
nmasagedsointr74= 0 ! no aged DU in this setup
nmasmergduintr34= 0 ! no aged DU in this setup

! .. below are new ones for no3, nh4, and nano3 to be integrated
nmasprimntnucsol= 0
nmasprimntaitsol= 0
nmasprimntaccsol= 0
nmasprimntcorsol= 0
nmasprimnhnucsol= 0
nmasprimnhaitsol= 0
nmasprimnhaccsol= 0
nmasprimnhcorsol= 0
nmascondnnaccsol= 0
nmascondnncorsol= 0
nmasddepntaitsol= 0
nmasddepntaccsol= 0
nmasddepntcorsol= 0
nmasddepnhaitsol= 0
nmasddepnhaccsol= 0
nmasddepnhcorsol= 0
nmasddepnnaccsol= 0
nmasddepnncorsol= 0
nmasnuscntaitsol= 0
nmasnuscntaccsol= 0
nmasnuscntcorsol= 0
nmasnuscnhaitsol= 0
nmasnuscnhaccsol= 0
nmasnuscnhcorsol= 0
nmasnuscnnaccsol= 0
nmasnuscnncorsol= 0
nmasimscntaitsol= 0
nmasimscntaccsol= 0
nmasimscntcorsol= 0
nmasimscnhaitsol= 0
nmasimscnhaccsol= 0
nmasimscnhcorsol= 0
nmasimscnnaccsol= 0
nmasimscnncorsol= 0
nmascoagntintr23= 0
nmascoagnhintr23= 0
nmascoagntintr24= 0
nmascoagnhintr24= 0
nmascoagntintr34= 0
nmascoagnhintr34= 0
nmascoagnnintr34= 0
nmasmergntintr23= 0
nmasmergnhintr23= 0
nmasmergntintr34= 0
nmasmergnhintr34= 0
nmasmergnnintr34= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE UKCA_INDICES_DUonly_2MODE

! ######################################################################
SUBROUTINE ukca_indices_sussbcocntnh_5mode

IMPLICIT NONE
!---------------------------------------------------------------
!
! Main array lengths and switches
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_SUSSBCOCNTNH_5MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Main array lengths and switches
ntraer=28          ! # of aerosol advected tracers
nbudaer=153        ! # of aerosol budget fields

nmasprimsuaitsol= 1
nmasprimsuaccsol= 2
nmasprimsucorsol= 3
nmasprimssaccsol= 4
nmasprimsscorsol= 5
nmasprimbcaitsol= 0 ! BC only emitted to insoluble
nmasprimbcaitins= 6
nmasprimocaitsol= 7
nmasprimocaitins= 8
nmasprimntnucsol= 0
nmasprimntaitsol= 9 ! NO3 emissions into ait - cor sol
nmasprimntaccsol=10
nmasprimntcorsol=11
nmasprimnhnucsol= 0
nmasprimnhaitsol=12 ! NH4 emissions into ait - cor sol
nmasprimnhaccsol=13
nmasprimnhcorsol=14
nmascondnnaccsol=15 ! NaNO3 emissions scheme into acc/coa sol
nmascondnncorsol=16 ! NaNO3 emissions scheme into acc/coa sol
!
nmasddepsunucsol=17
nmasddepsuaitsol=18
nmasddepsuaccsol=19
nmasddepsucorsol=20
nmasddepssaccsol=21
nmasddepsscorsol=22
nmasddepbcaitsol=23
nmasddepbcaccsol=24
nmasddepbccorsol=25
nmasddepbcaitins=26
nmasddepocnucsol=27
nmasddepocaitsol=28
nmasddepocaccsol=29
nmasddepoccorsol=30
nmasddepocaitins=31
nmasddepsonucsol= 0 ! SO stored in OC cpt
nmasddepsoaitsol= 0 ! SO stored in OC cpt
nmasddepsoaccsol= 0 ! SO stored in OC cpt
nmasddepsocorsol= 0 ! SO stored in OC cpt
nmasddepntaitsol=32
nmasddepntaccsol=33
nmasddepntcorsol=34
nmasddepnhaitsol=35
nmasddepnhaccsol=36
nmasddepnhcorsol=37
nmasddepnnaccsol=38
nmasddepnncorsol=39
!
nmasnuscsunucsol=40
nmasnuscsuaitsol=41
nmasnuscsuaccsol=42
nmasnuscsucorsol=43
nmasnuscssaccsol=44
nmasnuscsscorsol=45
nmasnuscbcaitsol=46
nmasnuscbcaccsol=47
nmasnuscbccorsol=48
nmasnuscbcaitins=49
nmasnuscocnucsol=50
nmasnuscocaitsol=51
nmasnuscocaccsol=52
nmasnuscoccorsol=53
nmasnuscocaitins=54
nmasnuscsonucsol= 0 ! SO stored in OC cpt
nmasnuscsoaitsol= 0 ! SO stored in OC cpt
nmasnuscsoaccsol= 0 ! SO stored in OC cpt
nmasnuscsocorsol= 0 ! SO stored in OC cpt
nmasnuscntaitsol=55
nmasnuscntaccsol=56
nmasnuscntcorsol=57
nmasnuscnhaitsol=58
nmasnuscnhaccsol=59
nmasnuscnhcorsol=60
nmasnuscnnaccsol=61
nmasnuscnncorsol=62
!
nmasimscsunucsol=63
nmasimscsuaitsol=64
nmasimscsuaccsol=65
nmasimscsucorsol=66
nmasimscssaccsol=67
nmasimscsscorsol=68
nmasimscbcaitsol=69
nmasimscbcaccsol=70
nmasimscbccorsol=71
nmasimscbcaitins=72
nmasimscocnucsol=73
nmasimscocaitsol=74
nmasimscocaccsol=75
nmasimscoccorsol=76
nmasimscocaitins=77
nmasimscsonucsol= 0 ! SO stored in OC cpt
nmasimscsoaitsol= 0 ! SO stored in OC cpt
nmasimscsoaccsol= 0 ! SO stored in OC cpt
nmasimscsocorsol= 0 ! SO stored in OC cpt
nmasimscntaitsol=78
nmasimscntaccsol=79
nmasimscntcorsol=80
nmasimscnhaitsol=81
nmasimscnhaccsol=82
nmasimscnhcorsol=83
nmasimscnnaccsol=84
nmasimscnncorsol=85
!
nmasclprsuaitsol1=86
nmasclprsuaccsol1=87
nmasclprsucorsol1=88
nmasclprsuaitsol2=89
nmasclprsuaccsol2=90
nmasclprsucorsol2=91
!
nmascondsunucsol=92
nmascondsuaitsol=93
nmascondsuaccsol=94
nmascondsucorsol=95
nmascondsuaitins=96
nmasnuclsunucsol=97
nmascondocnucsol=98
nmascondocaitsol=99
nmascondocaccsol=100
nmascondoccorsol=101
nmascondocaitins=102
nmascondsonucsol= 0 ! SO stored in OC cpt
nmascondsoaitsol= 0 ! SO stored in OC cpt
nmascondsoaccsol= 0 ! SO stored in OC cpt
nmascondsocorsol= 0 ! SO stored in OC cpt
nmascondsoaitins= 0 ! SO stored in OC cpt
!
nmascoagsuintr12=103
nmascoagsuintr13=104
nmascoagsuintr14=105
nmascoagsuintr15=106
nmascoagocintr12=107
nmascoagocintr13=108
nmascoagocintr14=109
nmascoagocintr15=110
nmascoagsointr12= 0 ! stored in NMASCOAGOCINTR12
nmascoagsointr13= 0 ! stored in NMASCOAGOCINTR13
nmascoagsointr14= 0 ! stored in NMASCOAGOCINTR14
nmascoagsointr15= 0 ! stored in NMASCOAGOCINTR15
nmascoagsuintr23=111
nmascoagbcintr23=112
nmascoagocintr23=113
nmascoagsointr23= 0 ! stored in NMASCOAGOCINTR23
nmascoagntintr23=114
nmascoagnhintr23=115
nmascoagsuintr24=116
nmascoagbcintr24=117
nmascoagocintr24=118
nmascoagsointr24= 0 ! stored in NMASCOAGOCINTR24
nmascoagntintr24=119
nmascoagnhintr24=120
nmascoagsuintr34=121
nmascoagbcintr34=122
nmascoagocintr34=123
nmascoagssintr34=124
nmascoagsointr34= 0 ! stored in NMASCOAGOCINTR34
nmascoagntintr34=125
nmascoagnhintr34=126
nmascoagnnintr34=127
!
nmascoagbcintr53=128
nmascoagocintr53=129
nmascoagbcintr54=130
nmascoagocintr54=131
!
nmasagedsuintr52=132
nmasagedbcintr52=133
nmasagedocintr52=134
nmasagedsointr52= 0 ! stored in NMASAGEDOCINTR52
!
nmasmergsuintr12=135
nmasmergocintr12=136
nmasmergsointr12= 0 ! stored in NMASMERGOCINTR12
nmasmergsuintr23=137
nmasmergbcintr23=138
nmasmergocintr23=139
nmasmergsointr23=  0 ! stored in NMASMERGOCINTR12
nmasmergntintr23=140
nmasmergnhintr23=141
nmasmergsuintr34=142
nmasmergssintr34=143
nmasmergbcintr34=144
nmasmergocintr34=145
nmasmergsointr34=  0 ! stored in NMASMERGOCINTR34
nmasmergntintr34=146
nmasmergnhintr34=147
nmasmergnnintr34=148
nmasprocsuintr23=149
nmasprocbcintr23=150
nmasprococintr23=151
nmasprocsointr23=  0 ! stored in NMASPROCOCINTR23
nmasprocntintr23=152
nmasprocnhintr23=153

! .. below are new ones for dust & modes 6/7 to be integrated
nmasprimduaccsol= 0 ! no DU in this setup
nmasprimducorsol= 0 ! no DU in this setup
nmasprimduaccins= 0 ! no DU in this setup
nmasprimducorins= 0 ! no DU in this setup
nmasddepduaccsol= 0 ! no DU in this setup
nmasddepducorsol= 0 ! no DU in this setup
nmasddepduaccins= 0 ! no DU in this setup
nmasddepducorins= 0 ! no DU in this setup
nmasnuscduaccsol= 0 ! no DU in this setup
nmasnuscducorsol= 0 ! no DU in this setup
nmasnuscduaccins= 0 ! no DU in this setup
nmasnuscducorins= 0 ! no DU in this setup
nmasimscduaccsol= 0 ! no DU in this setup
nmasimscducorsol= 0 ! no DU in this setup
nmasimscduaccins= 0 ! no DU in this setup
nmasimscducorins= 0 ! no DU in this setup
nmascondsuaccins= 0 ! no acc-ins cor-ins in this setup
nmascondsucorins= 0 ! no acc-ins cor-ins in this setup
nmascondocaccins= 0 ! no acc-ins cor-ins in this setup
nmascondoccorins= 0 ! no acc-ins cor-ins in this setup
nmascondsoaccins= 0 ! no acc-ins cor-ins in this setup
nmascondsocorins= 0 ! no acc-ins cor-ins in this setup
nmascoagsuintr16= 0 ! no acc-ins cor-ins in this setup
nmascoagsuintr17= 0 ! no acc-ins cor-ins in this setup
nmascoagocintr16= 0 ! no acc-ins cor-ins in this setup
nmascoagocintr17= 0 ! no acc-ins cor-ins in this setup
nmascoagsointr16= 0 ! no acc-ins cor-ins in this setup
nmascoagsointr17= 0 ! no acc-ins cor-ins in this setup
nmascoagduintr34= 0 ! no DU in this setup
nmascoagduintr64= 0 ! no DU in this setup
nmasagedsuintr63= 0 ! no acc-ins cor-ins in this setup
nmasagedduintr63= 0 ! no DU in this setup
nmasagedocintr63= 0 ! no BC/OC/SO in this setup
nmasagedsointr63= 0 ! no BC/OC/SO in this setup
nmasagedsuintr74= 0 ! no acc-ins cor-ins in this setup
nmasagedduintr74= 0 ! no DU in this setup
nmasagedocintr74= 0 ! no BC/OC/SO in this setup
nmasagedsointr74= 0 ! no BC/OC/SO in this setup
nmasmergduintr34= 0 ! no DU in this setup

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_sussbcocntnh_5mode

! ######################################################################
SUBROUTINE ukca_indices_solinsol_6mode

IMPLICIT NONE
!---------------------------------------------------------------

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDICES_SOLINSOL_6MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Main array lengths and switches

ntraer=12         ! # of aerosol advected tracers
nbudaer=44        ! # of aerosol budget fields

!
! This configuration only considers two components:
! -- soluble, which uses sulfate tracers and diagnostics
! -- insoluble, which uses dust tracers and diagnostics
!
nmasprimsuaitsol= 1
nmasprimsuaccsol= 2
nmasprimsucorsol= 3
nmasprimssaccsol= 0
nmasprimsscorsol= 0
nmasprimbcaitsol= 0
nmasprimbcaitins= 0
nmasprimocaitsol= 0
nmasprimocaitins= 0
!
nmasddepsunucsol= 4
nmasddepsuaitsol= 5
nmasddepsuaccsol= 6
nmasddepsucorsol= 7
nmasddepssaccsol= 0
nmasddepsscorsol= 0
nmasddepbcaitsol= 0
nmasddepbcaccsol= 0
nmasddepbccorsol= 0
nmasddepbcaitins= 0
nmasddepocnucsol= 0
nmasddepocaitsol= 0
nmasddepocaccsol= 0
nmasddepoccorsol= 0
nmasddepocaitins= 0
nmasddepsonucsol= 0
nmasddepsoaitsol= 0
nmasddepsoaccsol= 0
nmasddepsocorsol= 0
!
nmasnuscsunucsol= 8
nmasnuscsuaitsol= 9
nmasnuscsuaccsol=10
nmasnuscsucorsol=11
nmasnuscssaccsol= 0
nmasnuscsscorsol= 0
nmasnuscbcaitsol= 0
nmasnuscbcaccsol= 0
nmasnuscbccorsol= 0
nmasnuscbcaitins= 0
nmasnuscocnucsol= 0
nmasnuscocaitsol= 0
nmasnuscocaccsol= 0
nmasnuscoccorsol= 0
nmasnuscocaitins= 0
nmasnuscsonucsol= 0
nmasnuscsoaitsol= 0
nmasnuscsoaccsol= 0
nmasnuscsocorsol= 0
!
nmasimscsunucsol=12
nmasimscsuaitsol=13
nmasimscsuaccsol=14
nmasimscsucorsol=15
nmasimscssaccsol= 0
nmasimscsscorsol= 0
nmasimscbcaitsol= 0
nmasimscbcaccsol= 0
nmasimscbccorsol= 0
nmasimscbcaitins= 0
nmasimscocnucsol= 0
nmasimscocaitsol= 0
nmasimscocaccsol= 0
nmasimscoccorsol= 0
nmasimscocaitins= 0
nmasimscsonucsol= 0
nmasimscsoaitsol= 0
nmasimscsoaccsol= 0
nmasimscsocorsol= 0
!
nmasclprsuaitsol1=16
nmasclprsuaccsol1=17
nmasclprsucorsol1=18
nmasclprsuaitsol2=19
nmasclprsuaccsol2=20
nmasclprsucorsol2=21
!
nmascondsunucsol=22
nmascondsuaitsol=23
nmascondsuaccsol=24
nmascondsucorsol=25
nmascondsuaitins= 0
nmasnuclsunucsol=26
nmascondocnucsol= 0
nmascondocaitsol= 0
nmascondocaccsol= 0
nmascondoccorsol= 0
nmascondocaitins= 0
nmascondsonucsol= 0
nmascondsoaitsol= 0
nmascondsoaccsol= 0
nmascondsocorsol= 0
nmascondsoaitins= 0
!
nmascoagsuintr12=27
nmascoagsuintr13=28
nmascoagsuintr14=29
nmascoagsuintr15= 0
nmascoagocintr12= 0
nmascoagocintr13= 0
nmascoagocintr14= 0
nmascoagocintr15= 0
nmascoagsointr12= 0
nmascoagsointr13= 0
nmascoagsointr14= 0
nmascoagsointr15= 0
nmascoagsuintr23=30
nmascoagbcintr23= 0
nmascoagocintr23= 0
nmascoagsointr23= 0
nmascoagsuintr24=31
nmascoagbcintr24= 0
nmascoagocintr24= 0
nmascoagsointr24= 0
nmascoagsuintr34=32
nmascoagbcintr34= 0
nmascoagocintr34= 0
nmascoagssintr34= 0
nmascoagsointr34= 0
!
nmascoagbcintr53= 0
nmascoagocintr53= 0
nmascoagbcintr54= 0
nmascoagocintr54= 0
!
nmasagedsuintr52= 0
nmasagedbcintr52= 0
nmasagedocintr52= 0
nmasagedsointr52= 0
!
nmasmergsuintr12=33
nmasmergocintr12= 0
nmasmergsointr12= 0
nmasmergsuintr23=34
nmasmergbcintr23= 0
nmasmergocintr23= 0
nmasmergsointr23= 0
nmasmergsuintr34=35
nmasmergssintr34= 0
nmasmergbcintr34= 0
nmasmergocintr34= 0
nmasmergsointr34= 0
nmasprocsuintr23=36
nmasprocbcintr23= 0
nmasprococintr23= 0
nmasprocsointr23= 0
!
! .. below are new ones for dust & modes 6/7 to be integrated
nmasprimduaccsol=  0 ! DU only emitted to modes 6/7 here
nmasprimducorsol=  0 ! DU only emitted to modes 6/7 here
nmasprimduaccins=37
nmasprimducorins=38
nmasddepduaccsol= 0
nmasddepducorsol= 0
nmasddepduaccins=39
nmasddepducorins=40
nmasnuscduaccsol= 0
nmasnuscducorsol= 0
nmasnuscduaccins=41
nmasnuscducorins=42
nmasimscduaccsol= 0
nmasimscducorsol= 0
nmasimscduaccins=43
nmasimscducorins=44
nmascondsuaccins= 0
nmascondsucorins= 0
nmascondocaccins= 0
nmascondoccorins= 0
nmascondsoaccins= 0
nmascondsocorins= 0
nmascoagsuintr16= 0
nmascoagsuintr17= 0
nmascoagocintr16= 0
nmascoagocintr17= 0
nmascoagsointr16= 0
nmascoagsointr17= 0
nmascoagduintr34= 0
nmascoagduintr64= 0
nmasagedsuintr63= 0
nmasagedduintr63= 0
nmasagedocintr63= 0
nmasagedsointr63= 0
nmasagedsuintr74= 0
nmasagedduintr74= 0
nmasagedocintr74= 0
nmasagedsointr74= 0
nmasmergduintr34= 0

! .. below are new ones for no3, nh4, and nano3 to be integrated
nmasprimntnucsol= 0
nmasprimntaitsol= 0
nmasprimntaccsol= 0
nmasprimntcorsol= 0
nmasprimnhnucsol= 0
nmasprimnhaitsol= 0
nmasprimnhaccsol= 0
nmasprimnhcorsol= 0
nmascondnnaccsol= 0
nmascondnncorsol= 0
nmasddepntaitsol= 0
nmasddepntaccsol= 0
nmasddepntcorsol= 0
nmasddepnhaitsol= 0
nmasddepnhaccsol= 0
nmasddepnhcorsol= 0
nmasddepnnaccsol= 0
nmasddepnncorsol= 0
nmasnuscntaitsol= 0
nmasnuscntaccsol= 0
nmasnuscntcorsol= 0
nmasnuscnhaitsol= 0
nmasnuscnhaccsol= 0
nmasnuscnhcorsol= 0
nmasnuscnnaccsol= 0
nmasnuscnncorsol= 0
nmasimscntaitsol= 0
nmasimscntaccsol= 0
nmasimscntcorsol= 0
nmasimscnhaitsol= 0
nmasimscnhaccsol= 0
nmasimscnhcorsol= 0
nmasimscnnaccsol= 0
nmasimscnncorsol= 0
nmascoagntintr23= 0
nmascoagnhintr23= 0
nmascoagntintr24= 0
nmascoagnhintr24= 0
nmascoagntintr34= 0
nmascoagnhintr34= 0
nmascoagnnintr34= 0
nmasmergntintr23= 0
nmasmergnhintr23= 0
nmasmergntintr34= 0
nmasmergnhintr34= 0
nmasmergnnintr34= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0
nmasprocntintr23= 0
nmasprocnhintr23= 0

!
!----------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_indices_solinsol_6mode

END MODULE ukca_setup_indices
