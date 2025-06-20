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
!  Does one chemistry time step of the UKCA-MODE aerosol model.
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
MODULE ukca_aero_step_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_AERO_STEP_MOD'

CONTAINS

SUBROUTINE ukca_aero_step(nbox,nchemg,nadvg,nbudaer,                           &
 nd,mdt,md,mdwat,s0g,drydp,wetdp,rhopar,dvol,wvol,sm,                          &
 aird,airdm3,rhoa,mfpa,dvisc,t,tsqrt,rh,RH_clr,                                &
 s,pmid,pupper,plower,                                                         &
 zo3,zho2,zh2o2,ustr,znot,surtp,                                               &
 crain,drain,crain_up,csnow,dsnow,                                             &
 fconv_conv,lowcloud,vfac,clf,                                                 &
 autoconv1d,                                                                   &
 dtc,dtz,nmts,nzts,lday,act,bud_aer_mas,                                       &
 rainout_on,iextra_checks,                                                     &
 imscav_on,wetox_on,ddepaer_on,sedi_on,iso2wetoxbyo3,                          &
 dryox_in_aer,wetox_in_aer,delso2,delso2_2,                                    &
 cond_on,nucl_on,coag_on,bln_on,icoag,                                         &
 imerge,fine_no3_prod_on,coarse_no3_prod_on,hno3_uptake_coeff,                 &
 ifuchs,idcmfp,icondiam,ibln,i_nuc_method,                                     &
 iactmethod,iddepaer,inucscav,ichem,lcvrainout,l_dust_mp_slinn_impc_scav,      &
 verbose,checkmd_nd,intraoff,                                                  &
 interoff,s0g_dot_condensable,lwc,clwc,pvol,pvol_wat,                          &
 jlabove,ilscat,n_merge_1d,height,htpblg)

!-----------------------------------------------------------------------
!  Inputs
!  ------
!  NBOX        : Number of grid boxes
!  nchemg      : Number of gas phase chemistry tracers
!  nadvg       : Number of gas phase advected tracers
!  nbudaer     : Number of aerosol budget fields
!  ND          : Aerosol ptcl number density for mode (cm^-3)
!  MDT         : Avg tot mass of aerosol ptcl in mode (particle^-1)
!  MD          : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!  MDWAT       : Molecular concentration of water (molecules per ptcl)
!  S0G         : Partial masses of gas phase species (kg per gridbox)
!  RHOPAR      : Particle density [incl. H2O & insoluble cpts] (kgm^-3)
!  DVOL        : Geometric mean dry volume for each mode (m^3)
!  WVOL        : Geometric mean wet volume for each mode (m^3)
!  SM          : Grid box mass of air (kg)
!  AIRD        : Number density of air (per cm3)
!  AIRDM3      : Number density of air (per m3)
!  RHOA        : Air density (kg/m3)
!  MFPA        : Mean free path of air (m)
!  DVISC       : Dynamic viscosity of air (kg m^-1 s^-1)
!  T           : Centre level temperature (K)
!  TSQRT       : Square-root of centre level temperature (K)
!  RH          : Relative humidity (dimensionless 0-1)
!  RH_clr      : Relative humidity of clear-sky portion (dimensionless 0-1)
!  S           : Specific humidity (kg/kg)
!  PMID        : Centre level pressure (Pa)
!  PUPPER      : Upper interface pressure (Pa)
!  PLOWER      : Lower interface pressure (Pa)
!  ZO3         : Backgrnd vmr of O3
!  ZHO2        : Backgrnd conc. of HO2 (molecules per cc)
!  ZH2O2       : Backgrnd conc. of H2O2 (molecules per cc)
!  USTR        : Surface friction velocity (m/s)
!  ZNOT        : Roughness length (m)
!  SURTP       : Surface type: 0=seasurf,1=landsurf,2=oversea,3=overland
!  CRAIN       : Rain rate for conv precip. in box (kgm^-2s^-1)
!  DRAIN       : Rain rate for dyn. precip. in box (kgm^-2s^-1)
!  csnow       : Snowfall rate for conv precip. in box (kgm^-2s^-1)
!  dsnow       : Snowfall rate for dyn. precip. in box (kgm^-2s^-1)
!  CRAIN_UP    : Rain rate for conv precip. in box above (kgm^-2s^-1)
!  FCONV_CONV  : Fraction of box condensate --> rain in 6 hours (conv)
!  LOWCLOUD    : Horizontal low cloud fraction
!  VFAC        : Vertical low cloud fraction
!  DTC         : Chemistry time step (s)
!  DTZ         : Competition (cond/nucl) time step (s)
!  NMTS        : Number of microphysics timesteps per DTC
!  NZTS        : Number of competition timesteps per microphysics time step
!  LDAY        : 1-Day, 0-Night
!  MODESOL     : Specifies whether mode is soluble/insoluble (=1/0)
!  ACT         : Particle dry radius above which activation is assumed
!  RAINOUT_ON  : Switch : rainout (nucl. scav.) is on/off (1/0)
!  IMSCAV_ON   : Switch : impaction scavenging is on/off (1/0)
!  WETOX_ON    : Switch : wet oxidation [+ cloudproc] is on/off (1/0)
!  DDEPAER_ON  : Switch : aerosol dry deposition is on/off (1/0)
!  SEDI_ON     : Switch : aerosol sedimentation is on/off (1/0)
!  ISO2WETOXBYO3:Switch : so2in-cloud ox by o3 (assumed pH) on/off (1/0)
!  DRYOX_IN_AER: Switch : update gas phase condensible concentrations
!                         on competition tstep in aerosol code? (1/0)
!  WETOX_IN_AER: Switch : calc. wet ox. of SO2 in aerosol code? (1/0)
!  COND_ON     : Switch : vapour condensation is  on/off (1/0)
!  NUCL_ON     : Switch : binary nucleation is on/off (1/0)
!  COAG_ON     : Switch : coagulation is on/off (1/0)
!  BLN_ON      : Switch : Boundary layer nucleation is on/off (1/0)
!  ICOAG       : KIJ method (1:GLOMAP, 2: M7, 3: UMorig, 4:UMorig MFPP)
!  IMERGE      : Switch to use mid-pts (=1), edges (2) or dynamic (=3)
!  FINE_NO3_PROD_ON   : Switch: controls if fine NO3 production is on
!  COARSE_NO3_PROD_ON : Switch: controls if coarse NO3 production is on
!  HNO3_UPTAKE_COEFF  : Uptake coefficient for HNO3 in fine NO3 production
!  IFUCHS      : Switch : Fuchs (1964) or Fuchs-Sutugin (1971) for CC
!  IACTMETHOD  : Switch : activation method (0=off,1=fixed ract,2=NSO3)
!  IDDEPAER    : Switch : dry dep method (1=as in Spr05, 2=incl. sedi)
!  IDCMFP      : Switch : diffusion/mfp  (1=as Gbin v1, 2=as Gbin v1_1)
!  ICONDIAM    : Switch : wet diam in UKCA_CONDEN (1=g.mean,2=condiam.)
!  IBLN        : Switch : BLN method (1=activation,2=kinetic,3=PNAS)
!  I_NUC_METHOD: Switch for nucleation (how to combine BHN and BLN)
! (1=initial Pandis94 approach (no BLN even if switched on) -- Do not use!!
! (2=binary homogeneous nucleation applying BLN to BL only if switched on)
!   note there is an additional switch i_bhn_method (local to CALCNUCRATE)
!   to switch between using Kulmala98 or Vehkamaki02 for BHN rate
! (3=use same rate at all levels either activation(IBLN=1), kinetic(IBLN=2),
!  PNAS(IBLN=3)
!  note that if I_NUC_METHOD=3 and IBLN=3 then also add on BHN rate as in PNAS.
!  I_NUC_METHOD: Switch : 1=use Pandis94 (no BLN), 2=use BHN with BLN if on
!                       : 3=apply BLN parameterization throughout column
!  INUCSCAV    : Switch : scheme for removal by nucl scav
!                (1=as GLOMAP [accsol,corsol], 2=scav params as Stier05]
!                 3=as (1) but no nucl scav of modes 6 & 7 )
!  ICHEM              : Switch for gas phase chemistry scheme (0=off,1=on)
!  LCVRAINOUT  : Switch : .FALSE. to disable convective rainout if done
!                          alongside convective transport
!  l_dust_mp_slinn_impc_scav : Switch : .TRUE. for new dust and microplastics
!                                              impaction scavenging scheme
!  VERBOSE     : If =1 prints min/max (ND,MDT etc) after each process
!                for 1st grid box (for box model tests)
!  INTRAOFF    : Switch to turn off intra-modal coagulation
!  INTEROFF    : Switch to turn off intra-modal coagulation
!  CHECKMD_ND  : Switch : check for values of MD, ND out of range? (1/0)
!  S0G_DOT_CONDENSABLE : Gas phase chemistry tendencies for change in
!                        condensable gas phase species (vmr per s)
!                        due to chem, dry and wet deposition.
!  LWC         : Cloud liquid water content [kg/m3]
!  CLWC        : Cloud liquid water content [kg/kg]
!  autoconv1d  : Autoconversion rate including accretion [kg.kg^-1.s^-1]
!  JLABOVE     : Index of box directly above this grid box
!  ILSCAT      : Land surface category (based on 9 UM landsurf types)
!  N_MERGE_1D  : Count number of mode-merges in each box, each mode
!  DELSO2      : S(IV) --> S(VI) by H2O2 (molecules per cc)
!                [input from gas phase chem. scheme if WETOX_IN_AER=0]
!  DELSO2_2    : S(IV) --> S(VI) by O3u  (molecules per cc)
!                [input from gas phase chem. scheme if WETOX_IN_AER=0]
!  HEIGHT      : Mid-level height of gridbox
!  HTPBLG      : Height of boundary-layer in gridbox vertical-column
!
!  Outputs
!  -------
!  ND          : Aerosol ptcl number density for mode (cm^-3)
!  MDT         : Avg tot mass of aerosol ptcl in mode (particle^-1)
!  MD          : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!  DRYDP       : Geometric mean dry diameter for each mode (m)
!  WETDP       : Geometric mean wet diameter for each mode (m)
!  MDWAT       : molecular concentration of water (molecules per ptcl)
!  S0G         : Partial masses of gas phase species (kg per gridbox)
!  BUD_AER_MAS: Aerosol budget terms for mass
!  DELSO2      : S(IV) oxidised to S(VI) by H2O2 (molecules per cc)
!                [input as zero if WETOX_IN_AER=1]
!  DELSO2_2    : S(IV) oxidised to S(VI) by O3 (molecules per cc)
!                [input as zero if WETOX_IN_AER=1]
!  PVOL        : Partial volumes of each component in each mode (m3)
!  PVOL_WAT    : Partial volume of water in each mode (m3)
!
!  Local Variables
!  ---------------
!  SO2         : Sulfur dioxide conc (molecules per cc)
!  SO2_VMR     : Sulfur dioxode vmr
!  H2SO4       : Sulfuric acid vapour conc (moleculesH2SO4/cc)
!  Sec_Org     : 1st stage oxidation product from terpene (molecules/cc)
!  H2O2        : Hygrogen peroxide conc (per cc)
!  GC          : Condensable vapour conc (molecules cpt per cc)
!  GCOLD       : Condensable vapour conc b4 process (molecules cpt/cm3)
!  DELH2O2     : Change in hygrogen peroxide conc due to process (/cm3)
!  FRAC_AQ_ACC : Fraction of wet oxidised SO2 -> soluble accum. mode
!  FRAC_AQ_COR : Fraction of wet oxidised SO2 -> soluble coarse mode
!  TOTWETOX    : Total wet oxidation rate of SO2 (molecules/cm3)
!  DELGC_COND: Change in vapour conc due to cond (molecules cpt/cm3)
!  DELGC_NUCL: Change in vapour conc due to nucl (molecules cpt/cm3)
!  DELH2SO4_NUCL: Change in H2SO4 conc due to nucl (molecules/cm3)
!  DELTAS0G    : Overall change in partial masses of condensable gases
!                after all NZTS competition steps (kg/box).
!  S0G_TO_GC   : Molar mass ratio : condensing gas to aerosol phase cpt
!                e.g. MM_GAS(msec_org)=0.150 kg/mol but it may condense
!                into organic carbon component with MM(CP_OC)=0.0168.
!                Need to apply conversion when going from gas-->aerosol
!  S0G_TO_GC_SEC_ORG  Molar mass ratio specifically for Sec_Org
!                Required as Sec_Org conc needs rescaling for ukca_calnucrate
!  AGETERM1    : Depletion rate of each component (molecules cpt/cc/DTZ)
!                from condensation onto the 4 insoluble modes.
!                Used for calculation of ageing rate in UKCA_AGEING.
!  AGETERM2    : Rate of accomodation of material to each insoluble mode
!                as a result of coagulation with smaller soluble modes
!                (in molecules cpt /cm3/DTZ)
!                Used for calculation of ageing rate in UKCA_AGEING.
!  PROCESS     : Character string to store process currently being done
!                for error checking in box model
!  KII_ARR     : Coag coeff for intra-modal coag (IMODE-IMODE) (cm^3/s)
!  KIJ_ARR     : Coag coeff for inter-modal coag (IMODE-JMODE) (cm^3/s)
!
!  Inputted by module UKCA_MODE_SETUP
!  ----------------------------------
!  NMODES      : Number of possible aerosol modes
!  NMODES_SOL  : Number of possible soluble aerosol modes
!  NMODES_INS  : Number of possible insoluble aerosol modes
!  NCP         : Number of possible aerosol components
!  MM          : Molar masses of components (kg/mole)
!  NUM_EPS     : Value of NEWN below which do not carry out process
!  CP_SU       : Component where sulfate is stored
!
!  Inputted by module UKCA_SETUP_INDICES
!  -------------------------------------
!  CONDENSABLE        : Logical variable defining which cpts are condensable
!  CONDENSABLE_CHOICE : Indices of aerosol components into which condensable
!                         gases go to
!  MH2O2              : Index of MM_GAS, WTRATC and S0G for H2O2
!  MH2O2F             : Index of MM_GAS, WTRATC and S0G for H2O2F
!  MH2SO4             : Index of MM_GAS, WTRATC and S0G for H2SO4
!  MM_GAS             : Molar masses for gas phase species (kg/mol)
!  msec_org           : For involatile organic species
!  MSOTWO             : Index of MM_GAS, WTRATC and S0G for SO2
!  nmasclprsuaccsol1  : For in-cl. ox. of SO2 by H2O2 to SU in accsol
!  nmasclprsuaccsol2  : For in-cl. ox. of SO2 by O3   to SU in accsol
!  nmasclprsucorsol1  : For in-cl. ox. of SO2 by H2O2 to SU in corsol
!  nmasclprsucorsol2  : For in-cl. ox. of SO2 by O3   to SU in corsol
!  nmasnuclsunucsol   : For nucln. of SU to nucsol
!
!-----------------------------------------------------------------------

USE ereport_mod,                             ONLY:                             &
    ereport

USE errormessagelength_mod,                  ONLY:                             &
    errormessagelength

USE parkind1,                                ONLY:                             &
    jprb,                                                                      &
    jpim

USE ukca_ageing_mod,                         ONLY:                             &
    ukca_ageing

USE ukca_calc_coag_kernel_mod,               ONLY:                             &
    ukca_calc_coag_kernel

USE ukca_calc_drydiam_mod,                   ONLY:                             &
    ukca_calc_drydiam

USE ukca_calcminmaxgc_mod,                   ONLY:                             &
    ukca_calcminmaxgc

USE ukca_calcminmaxndmdt_mod,                ONLY:                             &
    ukca_calcminmaxndmdt

USE ukca_calcnucrate_mod,                    ONLY:                             &
    ukca_calcnucrate

USE ukca_check_md_nd_mod,                    ONLY:                             &
    ukca_check_md_nd

USE ukca_cloudproc_mod,                      ONLY:                             &
    ukca_cloudproc

USE ukca_coagwithnucl_mod,                   ONLY:                             &
    ukca_coagwithnucl

USE ukca_coarse_no3_mod,                     ONLY:                             &
    ukca_coarse_no3

USE ukca_conden_mod,                         ONLY:                             &
    ukca_conden

USE ukca_config_specification_mod,           ONLY:                             &
    glomap_variables

USE ukca_ddepaer_mod,                        ONLY:                             &
    ukca_ddepaer

USE ukca_ddepaer_incl_sedi_mod,              ONLY:                             &
    ukca_ddepaer_incl_sedi

USE ukca_fine_no3_mod,                     ONLY:                               &
    ukca_fine_no3

USE ukca_impc_scav_mod,                      ONLY:                             &
    ukca_impc_scav

USE ukca_impc_scav_dust_mod,                 ONLY:                             &
    ukca_impc_scav_dust

USE ukca_mode_setup,                         ONLY:                             &
    nmodes,                                                                    &
    nmodes_sol,                                                                &
    nmodes_ins,                                                                &
    cp_su

USE ukca_rainout_mod,                        ONLY:                             &
    ukca_rainout

USE ukca_remode_mod,                         ONLY:                             &
    ukca_remode

USE ukca_setup_indices,                      ONLY:                             &
    condensable,                                                               &
    condensable_choice,                                                        &
    mh2o2,                                                                     &
    mh2o2f,                                                                    &
    mh2so4,                                                                    &
    mm_gas,                                                                    &
    msec_org,                                                                  &
    msotwo,                                                                    &
    nmasclprsuaccsol1,                                                         &
    nmasclprsuaccsol2,                                                         &
    nmasclprsucorsol1,                                                         &
    nmasclprsucorsol2,                                                         &
    nmasnuclsunucsol

USE ukca_volume_mode_mod,                    ONLY:                             &
    ukca_volume_mode

USE ukca_wetox_mod,                          ONLY:                             &
    ukca_wetox

USE ukca_types_mod,                          ONLY:                             &
    log_small,                                                                 &
    integer_32

USE umPrintMgr,                              ONLY:                             &
    umPrint,                                                                   &
    umMessage

USE yomhook,                                 ONLY:                             &
    lhook,                                                                     &
    dr_hook


IMPLICIT NONE
!
! .. Subroutine Interface
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: nchemg
INTEGER, INTENT(IN) :: nadvg
INTEGER, INTENT(IN) :: nbudaer
INTEGER, INTENT(IN) :: nmts
INTEGER, INTENT(IN) :: nzts
INTEGER, INTENT(IN) :: lday(nbox)
INTEGER, INTENT(IN) :: rainout_on
INTEGER, INTENT(IN) :: iextra_checks
INTEGER, INTENT(IN) :: imscav_on
INTEGER, INTENT(IN) :: wetox_on
INTEGER, INTENT(IN) :: ddepaer_on
INTEGER, INTENT(IN) :: sedi_on
INTEGER, INTENT(IN) :: iso2wetoxbyo3
INTEGER, INTENT(IN) :: dryox_in_aer
INTEGER, INTENT(IN) :: wetox_in_aer
INTEGER, INTENT(IN) :: cond_on
INTEGER, INTENT(IN) :: nucl_on
INTEGER, INTENT(IN) :: coag_on
INTEGER, INTENT(IN) :: bln_on
INTEGER, INTENT(IN) :: icoag
INTEGER, INTENT(IN) :: imerge
INTEGER, INTENT(IN) :: ifuchs
INTEGER, INTENT(IN) :: iactmethod
INTEGER, INTENT(IN) :: iddepaer
INTEGER, INTENT(IN) :: idcmfp
INTEGER, INTENT(IN) :: icondiam
INTEGER, INTENT(IN) :: ibln
INTEGER, INTENT(IN) :: i_nuc_method
INTEGER, INTENT(IN) :: inucscav
INTEGER, INTENT(IN) :: ichem
INTEGER, INTENT(IN) :: verbose
INTEGER, INTENT(IN) :: checkmd_nd
INTEGER, INTENT(IN) :: intraoff
INTEGER, INTENT(IN) :: interoff
INTEGER, INTENT(IN) :: fine_no3_prod_on
INTEGER, INTENT(IN) :: coarse_no3_prod_on
INTEGER, INTENT(IN) :: jlabove(nbox)
INTEGER, INTENT(IN) :: ilscat(nbox)
INTEGER (KIND=integer_32), INTENT(IN OUT) :: n_merge_1d(nbox,nmodes)
LOGICAL, INTENT(IN) :: lcvrainout
LOGICAL, INTENT(IN) :: l_dust_mp_slinn_impc_scav

REAL, INTENT(IN OUT) :: rhopar(nbox,nmodes)
REAL, INTENT(IN OUT) :: dvol(nbox,nmodes)
REAL, INTENT(IN OUT) :: wvol(nbox,nmodes)
REAL, INTENT(IN) :: sm(nbox)
REAL, INTENT(IN) :: aird(nbox)
REAL, INTENT(IN) :: airdm3(nbox)
REAL, INTENT(IN) :: rhoa(nbox)
REAL, INTENT(IN) :: mfpa(nbox)
REAL, INTENT(IN) :: dvisc(nbox)
REAL, INTENT(IN) :: t(nbox)
REAL, INTENT(IN) :: tsqrt(nbox)
REAL, INTENT(IN) :: rh(nbox)
REAL, INTENT(IN) :: RH_clr(nbox)
REAL, INTENT(IN) :: s(nbox)
REAL, INTENT(IN) :: pmid(nbox)
REAL, INTENT(IN) :: pupper(nbox)
REAL, INTENT(IN) :: plower(nbox)
REAL, INTENT(IN) :: zo3(nbox)
REAL, INTENT(IN) :: zho2(nbox)
REAL, INTENT(IN) :: zh2o2(nbox)
REAL, INTENT(IN) :: ustr(nbox)
REAL, INTENT(IN) :: znot(nbox)
REAL, INTENT(IN) :: surtp(nbox)
REAL, INTENT(IN) :: crain(nbox)
REAL, INTENT(IN) :: drain(nbox)
REAL, INTENT(IN) :: crain_up(nbox)
REAL, INTENT(IN) :: csnow(nbox)
REAL, INTENT(IN) :: dsnow(nbox)
REAL, INTENT(IN) :: fconv_conv(nbox)
REAL, INTENT(IN) :: lowcloud(nbox)
REAL, INTENT(IN) :: vfac(nbox)
REAL, INTENT(IN) :: clf(nbox)
REAL, INTENT(IN) :: autoconv1d(nbox)
REAL, INTENT(IN) :: dtc
REAL, INTENT(IN) :: dtz
REAL, INTENT(IN) :: act
REAL, INTENT(IN) :: s0g_dot_condensable(nbox,nchemg)
REAL, INTENT(IN) :: lwc(nbox)
REAL, INTENT(IN) :: clwc(nbox)
REAL, INTENT(IN) :: htpblg(nbox)
REAL, INTENT(IN) :: height(nbox)
REAL, INTENT(IN) :: hno3_uptake_coeff
!
! .. Outputs
REAL, INTENT(IN OUT) :: nd(nbox,nmodes)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables%ncp)
REAL, INTENT(IN OUT) :: drydp(nbox,nmodes)
REAL, INTENT(IN OUT) :: wetdp(nbox,nmodes)
REAL, INTENT(IN OUT) :: mdwat(nbox,nmodes)
REAL, INTENT(IN OUT) :: s0g(nbox,nadvg)
REAL, INTENT(IN OUT) :: delso2(nbox)
REAL, INTENT(IN OUT) :: delso2_2(nbox)
REAL, INTENT(IN OUT) :: pvol(nbox,nmodes,glomap_variables%ncp)
REAL, INTENT(IN OUT) :: pvol_wat(nbox,nmodes)
REAL, INTENT(IN OUT) :: bud_aer_mas(nbox,0:nbudaer)
!
!     Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
REAL,    POINTER :: mm(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: num_eps(:)

INTEGER :: errcode                ! Variable passed to ereport

INTEGER :: jv
INTEGER :: icp
INTEGER :: cp_so4                 ! index of aerosol cpt to condense
INTEGER :: cp_sec_org             ! index of cpt to condense (for sec_org)
INTEGER :: imts
INTEGER :: izts
REAL :: so2(nbox)
REAL :: h2o2(nbox)
REAL :: h2so4(nbox)
REAL :: sec_org(nbox)
REAL :: delh2o2(nbox)
REAL :: so2_vmr(nbox)
REAL :: gc(nbox,nchemg)
REAL :: gcold(nbox,nchemg)
REAL :: delgc_cond(nbox,nchemg)
REAL :: delgc_nucl(nbox,nchemg)
REAL :: deltas0g(nbox)
REAL :: delh2so4_nucl(nbox)
REAL :: ageterm1(nbox,nmodes_ins,nchemg)
REAL :: ageterm2(nbox,nmodes_sol,nmodes_ins,glomap_variables%ncp)
REAL :: s0g_to_gc
REAL :: s0g_to_gc_sec_org
REAL :: frac_aq_acc(nbox)
REAL :: frac_aq_cor(nbox)
REAL :: totwetox(nbox)
REAL :: kii_arr(nbox,nmodes)
REAL :: kij_arr(nbox,nmodes,nmodes)
CHARACTER(LEN=30) :: process
CHARACTER(LEN=2)  :: strizts
LOGICAL (KIND=log_small) :: mask(nbox)
LOGICAL (KIND=log_small) :: mask1(nbox)
LOGICAL (KIND=log_small) :: mask2(nbox)
LOGICAL (KIND=log_small) :: mask3(nbox)
REAL :: s_cond_s(nbox)

CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_AERO_STEP'

!
!     CHECKMD_ND set to 1/0 for whether/not to check for
!                values of ND, MD and MDT out of range.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
mm          => glomap_variables%mm
ncp         => glomap_variables%ncp
num_eps     => glomap_variables%num_eps

IF (checkmd_nd == 1) THEN
  process='At start of UKCA_AERO_STEP    '
  CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
END IF
!
!     VERBOSE included so that max,min,mean of ND/MD can be
!             written out and monitored after each process
!             when error checking in box model.
IF (verbose >= 2) THEN
  WRITE(umMessage,'(A40)') 'At start of UKCA_AERO_STEP :   ND,MD,MDT'
  CALL umPrint(umMessage,src='ukca_aero_step')
  CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
END IF
!
!     Calculate dry diameter & volume for initial remode check
CALL ukca_calc_drydiam( nbox, glomap_variables,                                &
                        nd, md, mdt, drydp, dvol )
!
!     Initial remode call (check if advection taken DRYDP out of bounds)
CALL ukca_remode(nbox,nbudaer,nd,md,mdt,drydp,                                 &
                 imerge,bud_aer_mas,n_merge_1d,pmid)

!
IF (checkmd_nd == 1) THEN
  process='Done UKCA_REMODE (1st call)   '
  CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
END IF
!
IF (verbose >= 2) THEN
  WRITE(umMessage,'(A40)') 'After UKCA_REMODE1 has updated ND,MD,MDT'
  CALL umPrint(umMessage,src='ukca_aero_step')
  CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
END IF
!
!     Recalculate dry and wet diameter and volume
CALL ukca_calc_drydiam( nbox, glomap_variables,                                &
                        nd, md, mdt, drydp, dvol )
!
CALL ukca_volume_mode(glomap_variables,nbox, nd,md,mdt,                        &
                      RH_clr,dvol,drydp,t,pmid,s,                              &
                      mdwat,wvol,wetdp,rhopar,pvol,pvol_wat)
!
!     Calculate impaction scavenging of aerosol (washout)
IF (imscav_on == 1) THEN
  !
  CALL ukca_impc_scav(nbox,nbudaer,nd,md,                                      &
 crain,drain,csnow,dsnow,wetdp,dtc,l_dust_mp_slinn_impc_scav,bud_aer_mas)

  ! Run new impaction scavenging scheme for UKCA dust
  IF (l_dust_mp_slinn_impc_scav) THEN
    CALL ukca_impc_scav_dust(nbox,nbudaer,nd,md,mdt,crain,drain,               &
                             wetdp,dtc,bud_aer_mas,iextra_checks)
  END IF
  !
  IF (checkmd_nd == 1) THEN
    process='Done impaction scavenging     '
    CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
  END IF
  !
  IF (verbose >= 2) THEN
    WRITE(umMessage,'(A35)') 'After UKCA_IMPC_SCAV hasupdatedND'
    CALL umPrint(umMessage,src='ukca_aero_step')
    CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
  END IF
  !
END IF
!
!     Calculate in cloud, Nucleation Scavenging
IF (rainout_on == 1) THEN
  !
  CALL ukca_rainout(nbox,nbudaer,nd,md,mdt,fconv_conv,                         &
                    crain,crain_up,clwc,clf,autoconv1d,t,dtc,                  &
                    bud_aer_mas,inucscav,lcvrainout,wetdp)
  !
  IF (checkmd_nd == 1) THEN
    process='Done nucleation scavenging    '
    CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
  END IF
  !
  IF (verbose >= 2) THEN
    WRITE(umMessage,'(A35)') 'After UKCA_RAINOUT has updated ND'
    CALL umPrint(umMessage,src='ukca_aero_step')
    CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
  END IF
  !
END IF
!
IF (ichem == 1) THEN
  !
  ! Calculate aqueous chemistry
  ! Initialise variables
  frac_aq_acc(:) = 0.0
  frac_aq_cor(:) = 0.0

  IF (wetox_on == 1) THEN
    !
    !       Initialise variables for aqueous chemistry
    so2     (:)=0.0
    h2o2    (:)=0.0
    delh2o2 (:)=0.0
    !
    IF (wetox_in_aer == 1) THEN
      !
      IF (msotwo > 0) THEN
        mask(:)=(s0g(:,msotwo) > 0.0)
        WHERE (mask(:))
          so2(:)=s0g(:,msotwo)*aird(:)/sm(:)
        END WHERE
        so2_vmr(:)=so2(:)/aird(:)
      END IF
      !
      IF (mh2o2 > 0) THEN

        ! .. if using coupled H2O2
        mask(:)=(s0g(:,mh2o2) > 0.0)
        WHERE (mask(:))
          h2o2(:)=s0g(:,mh2o2)*aird(:)/sm(:)
        END WHERE
        !
      ELSE
        !
        IF (mh2o2f > 0) THEN
          !
          ! .. if using semi-prognostic H2O2
          mask(:)=(s0g(:,mh2o2f) > 0.0)
          WHERE (mask(:))
            h2o2(:)=s0g(:,mh2o2f)*aird(:)/sm(:)
          END WHERE
          !
        END IF
        !
      END IF
      !
      ! .. below calculates DELSO2, DELSO2_2, DELH2O2 in the
      ! .. aerosol module (WETOX_IN_AER=1)
      !
      CALL ukca_wetox(nbox,nd,drydp,delso2,delso2_2,                           &
         delh2o2,lowcloud,vfac,so2_vmr,h2o2,zh2o2,zo3,                         &
         zho2,pmid,t,aird,s,dtc,lday,lwc,iso2wetoxbyo3)
      !
      IF (msotwo > 0) THEN
        !
        mask(:)=(s0g(:,msotwo) > 0.0)
        WHERE (mask(:))
          !
          !          Update SO2 partial mass after reaction with H2O2.
          s0g(:,msotwo)=s0g(:,msotwo)-delso2  (:)*sm(:)/aird(:)
          !          Update SO2 partial mass after reaction with O3.
          s0g(:,msotwo)=s0g(:,msotwo)-delso2_2(:)*sm(:)/aird(:)
          !
        END WHERE
        !
        mask(:)=(s0g(:,msotwo) < 0.0)
        WHERE (mask(:)) s0g(:,msotwo)=0.0
        !
      END IF
      !
      !        Update H2O2 partial mass after SO2 wetox & replenishment
      !                                      (both included in DELH2O2)
      IF (mh2o2 > 0) THEN
        ! .. if using fully coupled H2O2
        mask(:)=(s0g(:,mh2o2) > 0.0)
        WHERE (mask(:))
          s0g(:,mh2o2)=s0g(:,mh2o2)+delh2o2(:)*sm(:)/aird(:)
        END WHERE
      END IF
      !
      IF (mh2o2f > 0) THEN
        ! .. if using semi-prognostic H2O2
        mask(:)=(s0g(:,mh2o2f) > 0.0)
        WHERE (mask(:))
          s0g(:,mh2o2f)=s0g(:,mh2o2f)+delh2o2(:)*sm(:)/aird(:)
        END WHERE
      END IF
      !
      IF ( (mh2o2 <= 0) .AND. (mh2o2f <= 0) ) THEN
        cmessage = ' AQUEOUS on but MH2O2<=0 and MH2O2F <=0'
        errcode = 1
        CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
      END IF
      !
    END IF ! IF WETOX_IN_AER=1
    ! n.b. if WETOX_IN_AER.NE.1, then DELSO2,DELSO2_2 passed in to AERO_STEP
    !
  ELSE
    delso2  (:)=0.0
    delso2_2(:)=0.0
  END IF ! if WETOX_ON=1
  !
  totwetox(:)=delso2(:)+delso2_2(:) ! this is zero if WETOX_ON=0
  !
  IF (iactmethod > 0) THEN ! if cloud processing on
    !
    ! .. below cloud-processes those aerosol in the Aitken soluble
    ! .. mode which are larger than ACT to accumulation mode
    ! .. representing the fact that they have activated so
    ! .. making minimum between Aitsol and accsol respond to activation
    !
    CALL ukca_cloudproc(nbox,nbudaer,nd,md,mdt,drydp,                          &
                        lowcloud,vfac,act,iactmethod,bud_aer_mas)
    !
  END IF ! IF IACTMETHOD > 0
  !
  mask1(:)=((nd(:,3)+nd(:,4)) > num_eps(4))
  mask2(:)=(mask1(:) .AND. (nd(:,3) > num_eps(3)))
  mask3(:)=(mask1(:) .AND. (nd(:,4) > num_eps(4)))
  !
  WHERE (mask1(:))
    !
    !       Calculate number fractions to partition between acc and cor
    frac_aq_acc(:)=nd(:,3)/(nd(:,3)+nd(:,4))
    frac_aq_cor(:)=nd(:,4)/(nd(:,3)+nd(:,4))
    !
  END WHERE
  !
  WHERE (mask2(:))
    !       Update soluble accum mode H2SO4 mass due to aqueous chem.
    md(:,3,cp_su)=(md(:,3,cp_su)*nd(:,3)+                                      &
                   frac_aq_acc(:)*totwetox(:))/nd(:,3)
    mdt(:,3)=(mdt(:,3)*nd(:,3)+frac_aq_acc(:)*totwetox(:))/nd(:,3)
  END WHERE ! if some particles in soluble accum mode
  !
  IF (nmasclprsuaccsol1 > 0) THEN
    WHERE (mask2(:))
      bud_aer_mas(:,nmasclprsuaccsol1)=                                        &
        frac_aq_acc(:)*delso2  (:)*sm(:)/aird(:)
    END WHERE
  END IF
  !
  IF (nmasclprsuaccsol2 > 0) THEN
    WHERE (mask2(:))
      bud_aer_mas(:,nmasclprsuaccsol2)=                                        &
        frac_aq_acc(:)*delso2_2(:)*sm(:)/aird(:)
    END WHERE
  END IF
  !
  WHERE (mask3(:))
    !       Update soluble coarse mode H2SO4 mass due to aqueous chem.
    md(:,4,cp_su)=(md(:,4,cp_su)*nd(:,4)+                                      &
                   frac_aq_cor(:)*totwetox(:))/nd(:,4)
    mdt(:,4)=(mdt(:,4)*nd(:,4)+frac_aq_cor(:)*totwetox(:))/nd(:,4)
  END WHERE ! if some particles in soluble coarse mode
  !
  IF (nmasclprsucorsol1 > 0) THEN
    WHERE (mask3(:))
      bud_aer_mas(:,nmasclprsucorsol1)=                                        &
        frac_aq_cor(:)*delso2  (:)*sm(:)/aird(:)
    END WHERE
  END IF
  !
  IF (nmasclprsucorsol2 > 0) THEN
    WHERE (mask3(:))
      bud_aer_mas(:,nmasclprsucorsol2)=                                        &
        frac_aq_cor(:)*delso2_2(:)*sm(:)/aird(:)
    END WHERE
  END IF
  !
  IF (checkmd_nd == 1) THEN
    process='Done aqueous phase chemistry  '
    CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
  END IF
  !
  IF (verbose >= 2) THEN
    WRITE(umMessage,'(A35)') 'After UKCA_WETOX has updated MD,MDT'
    CALL umPrint(umMessage,src='ukca_aero_step')
    CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
  END IF
  !
END IF ! ICHEM=1
!
!     Calculate aerosol dry deposition
IF (ddepaer_on == 1) THEN
  IF (iddepaer == 1) THEN
    CALL ukca_ddepaer(nbox,nbudaer,nd,md,mdt,rhopar,znot,                      &
                      dtc,wetdp,ustr,pmid,pupper,plower,t,                     &
                      surtp,rhoa,mfpa,dvisc,bud_aer_mas,ilscat)

  END IF
  IF (iddepaer == 2) THEN
    CALL ukca_ddepaer_incl_sedi(nbox,nbudaer,nd,md,mdt,rhopar,znot,            &
                                dtc,wetdp,ustr,pmid,pupper,                    &
                                plower,t,surtp,rhoa,                           &
                                mfpa,dvisc,bud_aer_mas,                        &
                                jlabove,ilscat,sedi_on,sm)
  END IF
  !
  IF (checkmd_nd == 1) THEN
    process='Done aerosol dry deposition   '
    CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
  END IF
  !
  IF (verbose >= 2) THEN
    WRITE(umMessage,'(A30)') 'After DDEPAER has updated ND'
    CALL umPrint(umMessage,src='ukca_aero_step')
    CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
  END IF
  !
END IF ! DDEPAER_ON=1
!
!     Loop over number of microphysics timesteps NMTS
DO imts=1,nmts
  !
  !     Calculate rate of condensation of gases onto existing particles
  !     and nucleation rate
  !
  IF (ichem == 1) THEN
    DO jv=1,nchemg ! loop over gas phase components
      IF (condensable(jv)) THEN
        !
        !        Set index of cpt into which condensable gas to be stored
        icp=condensable_choice(jv)
        !
        !        Calculate ratio of gas phase to aerosol component molar masses
        s0g_to_gc=mm_gas(jv)/mm(icp)
        !
        mask(:)=(s0g(:,jv) > 0.0)
        !
        WHERE (mask(:))
          !
          !         GC is gas phase molecular concentration of condensables
          !         but molecules refer to molar mass of aerosol cpt
          gc(:,jv)=s0g_to_gc*s0g(:,jv)*aird(:)/sm(:)
          !
        END WHERE
        !
        WHERE (.NOT. mask(:))
          !
          gc(:,jv)=0.0
          !
        END WHERE
        !
        gcold(:,jv)=gc(:,jv)
        !
      END IF
    END DO
  END IF ! ichem == 1
  !
  !      Split microphysics substep into NZTS subsubsteps to allow for
  !      competition between nucleation and condensation (and gas phase
  !      production if DRYOX_IN_AER=1) to compete on short timesteps.
  !      Given that nitrate formation is via condensation, include new
  !      routines for coarse and fine mode nitrate production here
  !
  IF (coag_on == 1) THEN
    CALL ukca_calc_coag_kernel(nbox,kii_arr,kij_arr,                           &
     drydp,dvol,wetdp,wvol,rhopar,mfpa,dvisc,t,                                &
     coag_on,icoag)
  ELSE
    kii_arr(:,:) = 0.0
    kij_arr(:,:,:) = 0.0
  END IF
  !
  DO izts=1,nzts
    !
    IF (verbose >= 2) THEN
      WRITE(strizts,'(1i2.2)') izts
      process='Start of loop, IZTS='//strizts//' (GC)   '
      CALL ukca_calcminmaxgc(process,nbox,nchemg,gc)
      process='Start of loop, IZTS='//strizts//' (GCOLD)'
      CALL ukca_calcminmaxgc(process,nbox,nchemg,gcold)
    END IF
    !
    IF (ichem == 1) THEN
      !
      IF (dryox_in_aer > 0) THEN
        !
        DO jv=1,nchemg
          IF (condensable(jv)) THEN
            !
            ! Set index of cpt into which condensable gas will be stored
            icp=condensable_choice(jv)
            !
            ! Calculate ratio of gas phase to aerosol cpt molar masses
            s0g_to_gc=mm_gas(jv)/mm(icp)
            !
            gc(:,jv)=gc(:,jv)+dtz*s0g_dot_condensable(:,jv)                    &
                             *aird(:)*s0g_to_gc
            !
            ! .. update condensable gas phase species by chemical tendencies on
            ! .. competition tstep
            !
            ! .. n.b. S0G_DOT_CONDENSABLE is in vmr of S0G(:,JV) per s,
            ! .. so need to * by S0G_TO_GC
            !
          END IF ! if CONDENSABLE(JV)
        END DO ! loop JV=1,NCHEMG
        !
        IF (verbose >= 2) THEN
          WRITE(strizts,'(1i2.2)') izts
          process='Done chem upd, IZTS='//strizts//' (GC)   '
          CALL ukca_calcminmaxgc(process,nbox,nchemg,gc)
          process='Done chem upd, IZTS='//strizts//' (GCOLD)'
          CALL ukca_calcminmaxgc(process,nbox,nchemg,gcold)
        END IF
        !
      END IF ! DRYOX_IN_AER
      !
      !        Carry out uptake of condensable gas phase species onto aerosol
      IF (cond_on == 1) THEN
        CALL ukca_conden(nbox,nchemg,nbudaer,ifuchs,idcmfp,icondiam,           &
                         nd,tsqrt,rhoa,airdm3,dtz,wetdp,pmid,t,                &
                         md,mdt,gc,bud_aer_mas,                                &
                         delgc_cond,ageterm1,s_cond_s)
        !
        IF (verbose >= 2) THEN
          WRITE(strizts,'(1i2.2)') izts
          process='Done cond upd, IZTS='//strizts//' (GC)   '
          CALL ukca_calcminmaxgc(process,nbox,nchemg,gc)
          process='Done cond upd, IZTS='//strizts//' (GCOLD)'
          CALL ukca_calcminmaxgc(process,nbox,nchemg,gcold)
        END IF
        !
        IF (checkmd_nd == 1) THEN
          process='Done conden of H2SO4 & Sec_Org'
          CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
        END IF
        !
        IF (verbose >= 2) THEN
          WRITE(umMessage,'(A40)') 'After UKCA_CONDEN has updated MD,MDT'
          CALL umPrint(umMessage,src='ukca_aero_step')
          CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
        END IF
        !
      ELSE
        ageterm1(:,:,:)=0.0
      END IF ! if COND_ON=1
      !
      !        Do nitrate chemistry - split into fine and coarse production
      IF (fine_no3_prod_on == 1) THEN
        CALL ukca_fine_no3(nbox,nadvg,nbudaer,hno3_uptake_coeff,dtz,rhoa,aird, &
                           wetdp,RH_clr,t,sm,mfpa,nd,md,mdt,s0g,bud_aer_mas)
        !
        IF (checkmd_nd == 1) THEN
          process='Done fine nitrate production'
          CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
        END IF
        !
        IF (verbose >= 2) THEN
          WRITE(umMessage,'(A40)') 'After UKCA_FINE_NO3 has updated MD,MDT'
          CALL umPrint(umMessage,src='ukca_aero_step')
          CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
        END IF
        !
      END IF ! if FINE_NO3_PROD_ON=1
      !
      IF (coarse_no3_prod_on == 1) THEN
        CALL ukca_coarse_no3(nbox,nadvg,nbudaer,dtz,rhoa,aird,wetdp,RH_clr,t,  &
                             sm,mfpa,nd,md,mdt,s0g,bud_aer_mas)
        !
        IF (checkmd_nd == 1) THEN
          process='Done coarse nitrate production'
          CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
        END IF
        !
        IF (verbose >= 2) THEN
          WRITE(umMessage,'(A40)') 'After UKCA_COARSE_NO3 has updated MD,MDT'
          CALL umPrint(umMessage,src='ukca_aero_step')
          CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
        END IF
        !
      END IF ! if COARSE_NO3_PROD_ON=1
      !
      !        Calculate rate of binary H2SO4-H2O nucleation
      IF (nucl_on == 1) THEN
        !
        IF (mh2so4 > 0) THEN
          !
          ! Set index of aerosol cpt which condensed H2SO4 to be stored
          cp_so4=condensable_choice(mh2so4)
          !
          IF (cp_so4 > 0) THEN
            !
            ! Calculate ratio of gas phase to aerosol cpt molar mass (H2SO4)
            s0g_to_gc=mm_gas(mh2so4)/mm(cp_so4)
            !
            mask(:)=(gc(:,mh2so4) > 0.0)
            !
            WHERE (mask(:))
              h2so4(:)=gc(:,mh2so4)/s0g_to_gc
            END WHERE
            ! .. Set (competition-step-updated) H2SO4 molecular concentration
            ! .. (using MM_GAS) for calculation of BHN rate in UKCA_CALCNUCRATE
            !
            WHERE (.NOT. mask(:))
              h2so4(:)=0.0
            END WHERE
            !
            ! Note that gc has been calculated above and has taken into account
            !  the molecular weight ratio of sec_org compared to OC

            ! While above statement is true, ukca_calcnucrate requires
            ! concentration of gas phase Sec_Org (mr 150 g/mol) not organic
            ! carbon component (mr 16.8 g/mol)
            ! Hence rescaling is applied solely for ukca_calcnurate using
            ! factor s0g_to_gc_sec_org = 150/16.8

            cp_sec_org=condensable_choice(msec_org)
            s0g_to_gc_sec_org=mm_gas(msec_org)/mm(cp_sec_org) ! scaling factor


            IF (msec_org > 0) THEN
              sec_org(:) = gc(:,msec_org)/s0g_to_gc_sec_org
            ELSE
              sec_org(:) = 0.0
            END IF

            CALL ukca_calcnucrate(nbox, dtz, t, s, rh, aird, h2so4,            &
                                  delh2so4_nucl, sec_org, bln_on, ibln,        &
                                  i_nuc_method, height, htpblg, s_cond_s)
            !
            mask(:)=(h2so4(:) > 0.0)
            !
            WHERE (mask(:))
              gc(:,mh2so4)=h2so4(:)*s0g_to_gc
            END WHERE
            ! .. Update GC due to updated value of H2SO4 (bug fixed in gm5c)
            !
            WHERE (.NOT. mask(:))
              gc(:,mh2so4)=0.0
            END WHERE
            !
            delgc_nucl(:,mh2so4)=delh2so4_nucl(:)*s0g_to_gc
            !
            IF (nmasnuclsunucsol > 0)                                          &
             bud_aer_mas(:,nmasnuclsunucsol)=                                  &
             bud_aer_mas(:,nmasnuclsunucsol)+delgc_nucl(:,mh2so4)
            !
            IF (verbose >= 2) THEN
              WRITE(strizts,'(1i2.2)') izts
              process='Done nucl upd, IZTS='//strizts//' (GC)   '
              !
              CALL ukca_calcminmaxgc(process,nbox,nchemg,gc)
              process='Done nucl upd, IZTS='//strizts//' (GCOLD)'
              CALL ukca_calcminmaxgc(process,nbox,nchemg,gcold)
            END IF
          ELSE
            errcode=1
            cmessage='CP_SO4 <= 0'
            WRITE(umMessage,'(A25,I6)') cmessage,'CP_SO4 = ',cp_so4
            CALL umPrint(umMessage,src='ukca_aero_step')

            CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
          END IF
          !
        ELSE
          errcode=1
          cmessage='MH2SO4 <= 0'
          WRITE(umMessage,'(A25,E12.3)') cmessage,'MH2SO4 = ',mh2so4
          CALL umPrint(umMessage,src='ukca_aero_step')

          CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
        END IF ! if MH2SO4 <=0
        !
      ELSE
        IF (mh2so4 > 0) delgc_nucl(:,mh2so4)=0.0
      END IF
    ELSE
      IF (mh2so4 > 0) delgc_nucl(:,mh2so4)=0.0
    END IF
    !
    IF (verbose >= 2) THEN
      WRITE(umMessage,'(A35)') 'After DRYDIAM & VOL has updated ND'
      CALL umPrint(umMessage,src='ukca_aero_step')
      CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
    END IF
    !

    IF ((coag_on == 1) .OR. (nucl_on == 1)) THEN
      IF (verbose >= 2) THEN
        WRITE(umMessage,'(A35)') 'About to call UKCA_COAGWITHNUCL'
        CALL umPrint(umMessage,src='ukca_aero_step')
        CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
      END IF

      !        Update ND & MD due to combined coagulation-nucleation
      CALL ukca_coagwithnucl(nbox,nchemg,nbudaer,nd,md,mdt,delgc_nucl,dtz,     &
                             ageterm2,intraoff,interoff,                       &
                             bud_aer_mas,kii_arr,kij_arr,                      &
                             iextra_checks)

      IF (checkmd_nd == 1) THEN
        process='Done combined coag./nucleation'
        CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
      END IF

      IF (verbose >= 2) THEN
        IF ((coag_on == 1) .AND. (nucl_on == 1)) THEN
          WRITE(umMessage,'(A40)') 'After COAG & NUCL have updatedND,MD,MDT'
          CALL umPrint(umMessage,src='ukca_aero_step')
        END IF
        IF ((coag_on == 0) .AND. (nucl_on == 1)) THEN
          WRITE(umMessage,'(A35)') 'After NUCL has updated ND,MD,MDT'
          CALL umPrint(umMessage,src='ukca_aero_step')
        END IF
        IF ((coag_on == 1) .AND. (nucl_on == 0)) THEN
          WRITE(umMessage,'(A35)') 'After COAG has updated ND,MD,MDT'
          CALL umPrint(umMessage,src='ukca_aero_step')
        END IF
        CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
      END IF
    ELSE
      ageterm2(:,:,:,:)=0.0
    END IF ! if COAG_ON = 1 or NUCL_ON = 1
    !
    !       Apply ageing -- transfer ND,MD from insol. to sol. modes
    CALL ukca_ageing(nbox, nchemg, nbudaer, nd, md, mdt, ageterm1, ageterm2,   &
                     wetdp, bud_aer_mas)
    !
    IF (checkmd_nd == 1) THEN
      process='Done ageing of insoluble modes'
      CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
    END IF
    !
    IF (verbose >= 2) THEN
      WRITE(umMessage,'(A40)') 'After UKCA_AGEING has updated ND,MD,MDT'
      CALL umPrint(umMessage,src='ukca_aero_step')
      CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
    END IF
    !
  END DO ! end loop over competition subsubtimesteps
  !
  !      Recalculate dry and wet diameter and volume
  CALL ukca_calc_drydiam( nbox, glomap_variables,                              &
                          nd, md, mdt, drydp, dvol )
  !
  CALL ukca_volume_mode(glomap_variables,nbox, nd,md,mdt,                      &
                        RH_clr,dvol,drydp,t,pmid,s,                            &
                        mdwat,wvol,wetdp,rhopar,pvol,pvol_wat)


  !
  !      Apply mode-merging where necessary
  CALL ukca_remode(nbox,nbudaer,nd,md,mdt,drydp,                               &
                   imerge,bud_aer_mas,n_merge_1d,pmid)
  !
  IF (checkmd_nd == 1) THEN
    process='Done UKCA_REMODE              '
    CALL ukca_check_md_nd(nbox,process,nd,md,mdt)
  END IF
  !
  IF (verbose >= 2) THEN
    WRITE(umMessage,'(A40)') 'After UKCA_REMODE has updated ND,MD,MDT'
    CALL umPrint(umMessage,src='ukca_aero_step')
    CALL ukca_calcminmaxndmdt(nbox,nd,mdt,verbose)
  END IF
  !
  ! Recalculate dry and wet diameter and volume after re-moding
  CALL ukca_calc_drydiam( nbox, glomap_variables,                              &
                          nd, md, mdt, drydp, dvol )

  !
  CALL ukca_volume_mode(glomap_variables,nbox, nd,md,mdt,                      &
                        RH_clr,dvol,drydp,t,pmid,s,                            &
                        mdwat,wvol,wetdp,rhopar,pvol,pvol_wat)

  IF (ichem == 1) THEN
    !
    DO jv=1,nchemg
      IF (condensable(jv)) THEN
        !
        icp=condensable_choice(jv)
        !
        ! Calculate ratio of gas phase to aerosol cpt molar masses
        s0g_to_gc=mm_gas(jv)/mm(icp)
        !
        deltas0g(:)=(gc(:,jv)-gcold(:,jv))*(sm(:)/aird(:))/s0g_to_gc
        !
        mask(:)=(deltas0g(:) < -s0g(:,jv))
        ! above limits deltaS0G to be -S0G (stop -ves)
        WHERE (mask(:))
          deltas0g(:)=-s0g(:,jv)
        END WHERE
        !
        s0g(:,jv)=s0g(:,jv)+deltas0g(:)
        !
      END IF ! IF CONDENSABLE(JV)
    END DO ! loop JV=1,NCHEMG
    !
  END IF ! if ICHEM = 1
  !
END DO ! end loop over NMTS

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_aero_step
END MODULE ukca_aero_step_mod
