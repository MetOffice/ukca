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
!    Calculate (wet) volume corresponding to mid-pt particles in each
!    aerosol mode.
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
MODULE ukca_volume_mode_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_VOLUME_MODE_MOD'

CONTAINS

SUBROUTINE ukca_volume_mode(nbox,nd,md,mdt,                                    &
 rh,wvol,wetdp,rhopar,                                                         &
 dvol,drydp,mdwat,pvol,pvol_wat,                                               &
 t,pmid,s)
!----------------------------------------------------------------------
!
! Calculate (wet) volume corresponding to mid-pt particles in each mode.
!
!                     calculate the water-content and density of the
!                     aerosol using ZSR and water activity coeffs for
!                     H+,SO42-,Na+,Cl- from Jacobsen (pg 610).
!                     Complete dissociation of solute ions is assumed
!                     to give the electrolyte concentrations and
!                     associated water content associated with each.
!                     The dry volume of insoluble molecules is then
!                     added to give overall (wet) ptcl volume.
!                     Note that OC is assumed water-insoluble in the
!                     insoluble mode but is assumed to have aged
!                     chemically in the aerosol to become hygroscopic
!                     once transferred to the soluble distribution.
!                     To respresent this, in the ZSR calculation,
!                     the concentration of SO4 ions is incremented
!                     by FHYG_AOM*MD(:,:,CP_OC)/F_AO -- i.e. the aged
!                     OC is assumed to take up water at a fraction
!                     FHYG_AOM (set at 0.65) of SO4.
!
! Sphericity is assumed to give the ptcl wet radius.
!
! Purpose
! -------
! Calculate avg wet volume WVOL & wet diameter WETDP for each aerosol
! particle size mode from the relative humidity and avg. number of
! molecules per particle (MD) of each cpt in each mode.
!
! (N.b. Local rel. hum. values are corrected to lie between 0.1 & 0.9)
!
! Inputs
! ------
! NBOX     : Number of grid boxes
! ND       : Aerosol ptcl no. concentration (ptcls per cc)
! MD       : Component median aerosol mass (molecules per ptcl)
! MDT      : Total median aerosol mass (molecules per ptcl)
! RH       : Relative humidity (corrected to lie within range 0.1-0.9)
! DVOL     : Dry volume of particle (m^3)
! DRYDP    : Dry diameter of particle (m)
! T        : Temperature
! PMID     : Pressure at mid level (Pa)
! S        : Specific humidity
!
! Outputs
! ------
! WVOL     : Avg wet volume of size mode (m3)
! WETDP    : Avg wet diameter of size mode (m)
! MDWAT    : Molecular concentration of water (molecules per particle)
! RHOPAR   : Particle density [incl. H2O & insoluble cpts] (kgm^-3)
! PVOL     : Partial volumes of each component in each mode (m3)
! PVOL_WAT : Partial volume of water in each mode (m3)
!
! Local Variables
! ---------------
! CORRH    : Locally corrected RH (to lie between 0.1 and 0.9)
! FHYG_AOM : Hygroscopicity of aged organic (fraction relative to SO4)
! F_AO     : No. of "moles" of POM in 1 mole of aged organic species
! MM_AGE_ORG: Molar mass of aged organic species (kg/mol)
! MM_POM   : Molar mass of particulate organic matter [for CP_OC,CP_SO]
! IONS     : Logical indicating presence of ions
! CL       : Ion concentrations (moles/cm3 of air)
! MASK     : Mask where in domain to calculate values
! MDSOL    : Mass per particle (total over soluble cpts) (mlcls/ptcl)
! RHOSOL   : Density of particle solution [excl. insoluble cpts] (kg/m3)
! DENOM    : Temporary variable calculating denominator of expression
! DENOM2   : Temporary variable calculating denominator of expression
! RHOTMP   : Temporary variable in calculation of particle density
! RHOTMP2  : Temporary variable in calculation of particle density
! WC       : Water content for aerosol (moles/cm3 of air)
! SIXOVRPIX: (6.0/pi)*{ 1.0/EXP((9/2)*LOG^2(SIGMA_G)) }
!
! Inputted by module CHEMISTRY_CONSTANTS, UKCA_CONSTANTS, WATER_CONST
! ---------------------------------
! AVOGADRO : Avogadro's constant (molecules per mole)
! RHO_WATER  : Density of water (=1000.0 kgm^-3)
! MMW      : Molecular mass of water (=0.018 kg/mol)
! RHO_SO4  : Mass density of a pure H2SO4 aerosol (kg per m^3)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES   : Number of possible aerosol modes
! NCP      : Number of possible aerosol components
! MODE     : Logical variable denoting where mode is defined
! COMPONENT: Logical variable denoting where cpt is defined
! MFRAC_0  : Initial mass fraction to set when no particles.
! SOLUBLE  : Logical variable defining which cpts are soluble
! MM       : Molar masses of components (kg per mole)
! MODESOL  : Defines which modes are soluble (integer)
! RHOCOMP  : Densities (dry) of each component (kg/m^3)
! X        : EXP((9/2)*LOG^2(SIGMA_G))
! NUM_EPS  : Value of NEWN below which don't recalculate MD (per cc)
!                                            or carry out process
! NANION   : Number of possible anion species
! NCATION  : Number of possible cation species
! CP_SU    : Index of component containing SO4    component
! CP_OC    : Index of component containing 1st OC component
! CP_CL    : Index of component containing NaCl   component
! CP_SO    : Index of component containing 2nd OC component
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! NCHEMG   : Number of gas phase tracers for gas phase chemistry scheme
! MM_GAS   : Array of molar masses for gas phase species (kg/mol)
!
!--------------------------------------------------------------------

USE ukca_um_legacy_mod,   ONLY: pi, l_glomap_clim_radaer, cubrt_v
USE ukca_types_mod,       ONLY: log_small
USE water_constants_mod,  ONLY: rho_water
USE chemistry_constants_mod, ONLY: avogadro, rho_so4
USE ukca_constants,       ONLY: mmw
USE ukca_mode_setup,      ONLY: nmodes, ncp, mode, component,                  &
                                mfrac_0, soluble, mm, modesol,                 &
                                rhocomp, x, num_eps, nanion,                   &
                                ncation, cp_su, cp_oc, cp_cl,                  &
                                cp_so, cp_nh4, cp_no3, cp_nn
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
USE ereport_mod,          ONLY: ereport
USE umPrintMgr,           ONLY: umPrint, umMessage, newline

USE errormessagelength_mod, ONLY: errormessagelength

USE ukca_vapour_mod, ONLY: ukca_vapour
USE ukca_water_content_v_mod, ONLY: ukca_water_content_v

USE ukca_config_specification_mod, ONLY: glomap_config

IMPLICIT NONE


! Subroutine interface
INTEGER, INTENT(IN) :: nbox
REAL, INTENT(IN)    :: nd(nbox,nmodes)
REAL, INTENT(IN)    :: md(nbox,nmodes,ncp)
REAL, INTENT(IN)    :: mdt(nbox,nmodes)
REAL, INTENT(IN)    :: rh(nbox)
REAL, INTENT(IN)    :: dvol(nbox,nmodes)
REAL, INTENT(IN)    :: drydp(nbox,nmodes)
REAL, INTENT(IN)    :: t(nbox)
REAL, INTENT(IN)    :: pmid(nbox)
REAL, INTENT(IN)    :: s(nbox)
REAL, INTENT(OUT)   :: mdwat(nbox,nmodes)
REAL, INTENT(OUT)   :: wvol(nbox,nmodes)
REAL, INTENT(OUT)   :: wetdp(nbox,nmodes)
REAL, INTENT(OUT)   :: rhopar(nbox,nmodes)
REAL, INTENT(OUT)   :: pvol(nbox,nmodes,ncp)
REAL, INTENT(OUT)   :: pvol_wat(nbox,nmodes)

! Local variables
INTEGER :: errcode                ! Variable passed to ereport
INTEGER :: i
INTEGER :: imode
INTEGER :: icp
INTEGER :: jjl(1)
INTEGER :: iimode
INTEGER :: iicp
INTEGER :: ierr(nbox)
LOGICAL (KIND=log_small) :: mask(nbox)
LOGICAL (KIND=log_small) :: mask_sol(nbox)
LOGICAL (KIND=log_small) :: mask_nosol(nbox)
LOGICAL (KIND=log_small) :: mask_error(nbox)
LOGICAL (KIND=log_small) :: ions(nbox,-nanion:ncation) !ION PRESENCE SWITCHES
REAL    :: corrh(nbox)
REAL    :: mdsol(nbox)
REAL    :: rhosol(nbox)
REAL    :: wc(nbox)
REAL    :: rhotmp(nbox)
REAL    :: rhotmp2(nbox)
REAL    :: denom(nbox)
REAL    :: denom2(nbox)
REAL    :: piovrsix
REAL    :: sixovrpix(nmodes)
REAL    :: mm_ovravc(ncp)
REAL    :: mmwovravc
REAL    :: mm_ovravcrhocp(ncp)
REAL    :: mm_rhocp(ncp)
REAL    :: mmwrhow
REAL    :: f_ao
REAL    :: cl(nbox,-nanion:ncation) !ION CONCS (MOL/CC OF AIR)

REAL    :: tmp1(nbox,nmodes)
!
! extra local vars for stratospheric calculation
REAL    :: rp(nbox)
REAL    :: wts(nbox)
REAL    :: rhosol_strat(nbox)
REAL    :: massh2so4kg(nbox)
REAL    :: masswaterkg(nbox)

CHARACTER(LEN=errormessagelength) :: cmessage
REAL, PARAMETER :: fhyg_aom=0.65
REAL, PARAMETER :: mm_age_org=0.150
REAL, PARAMETER :: mm_pom=0.0168
REAL, PARAMETER :: putls=1.5e4 ! set max pressure for UTLS region to 150hPa
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_VOLUME_MODE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!at this point in the code, the value of RP does not matter
rp(:)=100.0e-9 ! dummy value
CALL ukca_vapour(nbox,t,pmid,s,rp,wts,rhosol_strat)
! Note that here we only want to get WTS and RHOSOL_STRAT which
!  are independent of particle size and composition -- so we only
!  need to call this once here -0- .

piovrsix=pi/6.0
sixovrpix(:)=1.0/(x(:)*piovrsix)
mm_ovravc(:)=mm(:)/avogadro
mmwovravc   =mmw/avogadro
mm_ovravcrhocp(:)=mm_ovravc(:)/rhocomp(:)
mm_rhocp(:)=mm(:)*rhocomp(:)
mmwrhow    =mmw*rho_water

denom(:)=0.0         ! define on all points to avoid error
denom2(:)=0.0

! Correct relative humidities to lie within the range of 10-90%
corrh(:)=rh(:)
WHERE (corrh(:) > 0.9) corrh(:)=0.9
WHERE (corrh(:) < 0.1) corrh(:)=0.1

DO imode=1,nmodes
  IF (mode(imode)) THEN

    mask(:) =(nd(:,imode) > num_eps(imode))

    IF (modesol(imode) == 1) THEN

      ! .. first calculate the total mass per particle over all soluble cpts
      mdsol(:)=0.0
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          IF (soluble(icp)) THEN
            ! .. calculate total mass (in molecules per particle) of all solutes
            WHERE (mask(:)) mdsol(:)=mdsol(:)+md(:,imode,icp)
          END IF
        END IF
      END DO

      ! .. make mask "MASK_SOL"   for where ND>NUM_EPS and some soluble mass
      ! .. make mask "MASK_NOSOL" for where ND>NUM_EPS and no   soluble mass
      mask_sol  (:) = (mask(:) .AND. (mdsol(:) > 0.0))
      mask_nosol(:) = (mask(:) .AND. (mdsol(:) == 0.0))

      f_ao=mm_age_org/mm_pom

        ! .. Use composition information to calculate water uptake by
        ! .. each component according to ZSR using water activity data from
        ! .. Jacobsen page 610 (Table B.10) for binary electrolyte molalities

        ! **********************************************************************
        ! *  Liquid Phase Species:
        ! *
        ! *  **Cations**           **Anions**             **Neutrals**
        ! *  1: H                  -1: HSO4               0: H2O
        ! *  2: NH4                -2: SO4
        ! *  3: Na                 -3: NO3
        ! *                        -4: Cl
        ! **********************************************************************

      mask(:) =(nd(:,imode) > num_eps(imode))

      DO i=-nanion,ncation
        cl(:,i)=0.0 ! set all concentrations to zero initially
      END DO

      IF (component(imode,cp_su)) THEN ! assume all H2SO4 --> SO4
        cl(:,-2)=md(:,imode,cp_su)/avogadro  ! [SO4] in moles/cc (air)
      END IF

      IF (component(imode,cp_so)) THEN
        cl(:,-2)=cl(:,-2)+(fhyg_aom/avogadro)*(md(:,imode,cp_so)/f_ao)
        ! .. Increment concentration of SO4 ions to represent the
        ! .. presence of hygroscopic aged organic aerosol mass in CP_SO.
        ! .. Assume it has uptake behaviour at fraction FHYG_AOM of SO4.
        ! .. Need to divide by F_AO because MD of CP_SO is in "moles" of POM
        ! .. whereas CL needs to be in moles of aged organic (MM=0.15kg/mol).
      END IF

      IF (component(imode,cp_oc)) THEN
        cl(:,-2)=cl(:,-2)+(fhyg_aom/avogadro)*(md(:,imode,cp_oc)/f_ao)
        ! .. Increment concentration of SO4 ions to represent the
        ! .. presence of hygroscopic aged organic aerosol mass in CP_OC.
        ! .. Assume it has uptake behaviour at fraction FHYG_AOM of SO4.
        ! .. Need to divide by F_AO because MD of CP_OC is in "moles" of POM
        ! .. whereas CL needs to be in moles of aged organic (MM=0.15kg/mol).
        !
        ! .. This effectively says that by the time the primary carbonaceous
        ! .. aerosol has been microphysically aged to the soluble mode, the
        ! .. organic component has been chemically aged to become hygroscopic.
        ! .. (it is assumed to be water-insoluble in the insoluble mode)
        !
      END IF

      IF (component(imode,cp_cl)) THEN ! assume complete dissociation
        cl(:,3)=md(:,imode,cp_cl)/avogadro ! [Na] in moles per cc (air)
        cl(:,-4)=md(:,imode,cp_cl)/avogadro ! [Cl] in moles per cc (air)
      END IF

      ! Nitrate scheme
      IF (UBOUND(component,DIM=2) >= cp_no3) THEN

        IF (component(imode,cp_nh4)) THEN ! assume complete dissociation
          cl(:,2)=md(:,imode,cp_nh4)/avogadro ! [NH4] in moles per cc (air)
        END IF

        IF (component(imode,cp_no3)) THEN ! assume complete dissociation
          cl(:,-3)=md(:,imode,cp_no3)/avogadro ! [NO3] in moles per cc (air)
        END IF

        IF (component(imode,cp_nn)) THEN ! assume complete dissociation
          cl(:,3)=cl(:,3)+md(:,imode,cp_nn)/avogadro
                    ! [Na] in moles per cc (air)
          cl(:,-3)=cl(:,-3)+md(:,imode,cp_nn)/avogadro
                    ! [NO3] in moles per cc (air)
        END IF

      END IF

      ! SET H+ FOR CHARGE BALANCE  -- CL(1) is [H] in moles per cc(air)
      cl(:,1)=MAX((2.0*cl(:,-2)+cl(:,-1)+cl(:,-3)+cl(:,-4)                     &
                   -cl(:,2)-cl(:,3)),0.0)

      DO i=-nanion,ncation
        ions(:,i)=(cl(:,i) > 0.0)
      END DO

      CALL ukca_water_content_v(nbox,mask,cl,corrh,ions,wc)

      WHERE (mask(:)) mdwat(:,imode)=wc(:)*avogadro
      !
      ! over-write MDWAT in strat with value from WTS retured from UKCA_VAPOUR
      WHERE (mask(:) .AND. (pmid(:) <  putls)) ! P<PUTLS
        massh2so4kg(:)=md(:,imode,cp_su)*mm(cp_su)/avogadro ! in kg/particle
        masswaterkg(:)=(100.0/wts(:)-1.0)*massh2so4kg(:) ! kg/ptcl
        mdwat(:,imode)=masswaterkg(:)/mmwovravc ! in molecules per particle
      END WHERE
      !
      WHERE (mask(:))
        ! .. calculate solution density (avg over each cpt mass contribution)
        rhotmp (:)=mdwat(:,imode)*mmwrhow
        rhotmp2(:)=rhotmp(:)
        denom (:)=mdwat(:,imode)*mmw
        denom2(:)=denom(:)
        ! ..note only need to define DENOM,DENOM2,RHOTMP,RHOTMP2 where MASK(:)
      ELSE WHERE
        mdwat(:,imode)=0.0
        ! .. set MDWAT to be zero where ND<NUM_EPS
        rhotmp(:)=0.0
        denom (:)=0.0
        rhotmp2(:)=0.0
        denom2 (:)=0.0
        ! .. also set RHOTMP,DENOM,RHOTMP2,DENOM2=0 (they will not be used)
      END WHERE

      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          IF (soluble(icp)) THEN
            WHERE (mask(:))
              rhotmp(:)=rhotmp(:)+md(:,imode,icp)*mm_rhocp(icp)
              denom (:)=denom (:)+md(:,imode,icp)*mm(icp)
            END WHERE
          END IF
          WHERE (mask(:))
            rhotmp2(:)=rhotmp2(:)+md(:,imode,icp)*mm_rhocp(icp)
            denom2 (:)=denom2 (:)+md(:,imode,icp)*mm(icp)
          END WHERE
        END IF
      END DO

      ! .. RHOSOL is density of ptcl solution [exclud. insoluble cpts] (kgm/3)
      !
      !  .. section below added to check for DENOM <= 0.0 or DENOM2 <= 0.0

      mask_error(:)= (mask_sol(:) .AND. (denom(:) <= 0.0))
      ierr(:)=0
      WHERE (mask_error(:))
        ierr(:)=1
      END WHERE
      IF (SUM(ierr(:)) > 0) THEN ! error (DENOM<=0 when some soluble cpt)
        DO i=1,nbox
          IF ((mask(i)) .AND. (denom(i) <= 0.0)) THEN

            WRITE(umMessage,'(A)') 'DENOM(i) <= 0.0'
            CALL umPrint(umMessage,src=RoutineName)

            WRITE(umMessage,'(A,I0,2E15.6,4X,L1)') 'i, rhotmp(i), ' //         &
                                                   'denom(i), mask(i) : ',     &
                                                    i, rhotmp(i),              &
                                                    denom(i), mask(i)
            CALL umPrint(umMessage,src=RoutineName)

            WRITE(umMessage,'(A,I0,4X,I0,2E15.6)') 'imode, ' //                &
                                                   'modesol(imode), ' //       &
                                                   'nd(i,imode), ' //          &
                                                   'num_eps(imode) : ',        &
                                                    imode,                     &
                                                    modesol(imode),            &
                                                    nd(i,imode),               &
                                                    num_eps(imode)
            CALL umPrint(umMessage,src=RoutineName)

            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                WRITE(umMessage,'(A,I0,3E15.6)') 'icp, mm(icp), ' //           &
                                                 'mm_rhocp(icp), ' //          &
                                                 'md(i,imode,icp) : ',         &
                                                  icp, mm(icp),                &
                                                  mm_rhocp(icp),               &
                                                  md(i,imode,icp)
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            END DO

            WRITE(umMessage,'(A,2E15.6)') 'mdwat(i,imode), mmw : ',            &
                                           mdwat(i,imode), mmw
            CALL umPrint(umMessage,src=RoutineName)

            errcode=ABS(i)
            WRITE(cmessage,'(3(A,I0),4(A,E15.6))')                             &
                 'Demoninator <= 0 for i=', i, ' , imode=',imode,              &
                  newline//                                                    &
                 'modesol(imode)=', modesol(imode),' , denom(i) = ',denom(i),  &
                  newline//                                                    &
                 'rhotmp(i) = ',rhotmp(i),' , nd(i,imode) = ', nd(i,imode),    &
                  newline//                                                    &
                 'num_eps(imode) = ',num_eps(imode)
            CALL ereport(RoutineName,errcode,cmessage)

          END IF ! if DENOM(I) <= 0.0

          IF ((mask(i)) .AND. (denom2(i) <= 0.0)) THEN

            WRITE(umMessage,'(A)') 'DENOM2(I) <= 0.0'
            CALL umPrint(umMessage,src=RoutineName)

            WRITE(umMessage,'(A,I0,2E15.6,4X,L1)') 'i, rhotmp2(i), ' //        &
                                                   'denom2(i), mask(i) : ',    &
                                                    i, rhotmp2(i),             &
                                                    denom2(i), mask(i)
            CALL umPrint(umMessage,src=RoutineName)

            WRITE(umMessage,'(A,I0,4X,I0)') 'imode, modesol(imode) : ',        &
                                             imode, modesol(imode)
            CALL umPrint(umMessage,src=RoutineName)

            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                WRITE(umMessage,'(A,I0,3E15.6)') 'icp, mm(icp), ' //           &
                                                 'mm_rhocp(icp), ' //          &
                                                 'md(i,imode,icp) : ',         &
                                                  icp, mm(icp),                &
                                                  mm_rhocp(icp),               &
                                                  md(i,imode,icp)
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            END DO

            WRITE(umMessage,'(A,2E15.6)') 'mdwat(i,imode), mmw : ',            &
                                           mdwat(i,imode), mmw
            CALL umPrint(umMessage,src=RoutineName)

            errcode=ABS(i)
            WRITE(cmessage,'(3(A,I0),2(A,E15.6))')                             &
                 'Demoninator2 <= 0 for i=', i, ' , imode=', imode,            &
                  newline//                                                    &
                 'modesol(imode)=',modesol(imode),', denom2(i) = ',denom2(i),  &
                  newline//                                                    &
                  'rhotmp2(i) = ',rhotmp2(i)
            CALL ereport(RoutineName,errcode,cmessage)

          END IF ! if DENOM2(I) <= 0.0
        END DO ! loop over NBOX (I)
      END IF

      ! .. where ND>NUM_EPS and there is some soluble material
      WHERE (mask_sol(:))
        rhosol(:)=rhotmp(:)/denom(:)
        ! ..only need to define RHOSOL where ND>NUM_EPS and MDSOL>0 (MASK_SOL)
      END WHERE
      !
      !  over-write RHOSOL in strat with RHOSOL_STRAT from UKCA_VAPOUR
      WHERE (mask_sol(:) .AND. (pmid(:) <  putls)) ! P<PUTLS
        rhosol(:)=rhosol_strat(:)
      END WHERE
      !
      wvol(:,imode)=0.0        ! initialise wet volume to zero
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          IF (soluble(icp)) THEN
            WHERE (mask_sol(:))
              pvol(:,imode,icp)=md(:,imode,icp)*mm_ovravc(icp)/rhosol(:)
              ! .. for soluble cpts set PVOL according to cpt mass and
              !     solution density
              wvol(:,imode)=wvol(:,imode)+pvol(:,imode,icp)
            ELSE WHERE
              pvol(:,imode,icp)=dvol(:,imode)*mfrac_0(imode,icp)
              ! where ND<NUM_EPS, over-write PVOL with default value
            END WHERE
            WHERE (mask_nosol(:))
              pvol(:,imode,icp)=0.0
              ! .. In case where soluble mode has some insoluble cpt mass
              !     but no soluble component mass, DENOM=0.0
              !     so just set PVOL for soluble cpt =0.0 here
            END WHERE
          ELSE
            WHERE (mask(:))
              pvol(:,imode,icp)=md(:,imode,icp)*mm_ovravcrhocp(icp)
              ! .. for insoluble cpts set PVOL according to cpt mass and cpt
              !     density
              wvol(:,imode)=wvol(:,imode)+pvol(:,imode,icp)
            ELSE WHERE
              pvol(:,imode,icp)=dvol(:,imode)*mfrac_0(imode,icp)
              ! where ND<NUM_EPS, over-write PVOL with default value
            END WHERE
          END IF
          ! .. initially set WVOL to sum of partial volumes of each component
        END IF
      END DO

      WHERE (mask_sol(:))
        pvol_wat(:,imode)=mdwat(:,imode)*mmwovravc/rhosol(:)
        ! .. for water, set PVOL according to water mass and solution density
        wvol(:,imode)=wvol(:,imode)+pvol_wat(:,imode)
        ! .. add on partial volume of water to give total wet volume
        rhopar(:,imode)=rhotmp2(:)/denom2(:)
        ! .. RHOPAR is total particle density [incl H2O & insoluble cpts]
        !     (kgm^-3)
        ! .. calculated as mass-weighted mean of component & water densities
      ELSE WHERE
        pvol_wat(:,imode)=0.0
        wvol(:,imode)=dvol(:,imode)
        rhopar(:,imode)=rho_so4
        ! where ND<NUM_EPS, set PVOL_WAT=0.0, set WVOL=DVOL and RHOPAR=RHOSO4
      END WHERE

    ELSE  ! if mode not soluble

      ! where mode not soluble, set PVOL_WAT=0.0, set WVOL=DVOL and MDWAT=0.0
      pvol_wat(:,imode)=0.0
      wvol(:,imode)=dvol(:,imode)
      mdwat(:,imode)=0.0

      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          pvol(:,imode,icp)=md(:,imode,icp)*mm_ovravcrhocp(icp)
          ! .. all cpts insoluble in insoluble modes,
          ! .. so set PVOL according to cpt mass & cpt density
        END IF
      END DO

      rhotmp2(:)=0.0
      denom2 (:)=0.0
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WHERE (mask(:))
            rhotmp2(:)=rhotmp2(:)+md(:,imode,icp)*mm_rhocp(icp)
            denom2 (:)=denom2 (:)+md(:,imode,icp)*mm(icp)
          END WHERE
        END IF
      END DO

      WHERE (mask(:))
        rhopar(:,imode)=rhotmp2(:)/denom2(:)
        ! .. RHOPAR is total particle density [here only insoluble cpts](kgm^-3)
        ! .. calculated as mass-weighted mean of component densities
      ELSE WHERE
        rhopar(:,imode)=rho_so4
        ! .. if insignificant # of particles, set density to default value
      END WHERE

    END IF ! if mode is soluble

  ELSE

    ! .. set PVOL_WAT=MDWAT=0.0 and WVOL=DVOL if mode not present
    pvol_wat(:,imode)=0.0
    wvol(:,imode)=dvol(:,imode)
    mdwat(:,imode)=0.0

    ! .. set density to default value if mode not present
    rhopar(:,imode)=rho_so4
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        pvol(:,imode,icp)=dvol(:,imode)*mfrac_0(imode,icp)
        ! .. set partial volume to default value if mode not present
      END IF
    END DO
  END IF ! if(mode(imode))
END DO

DO imode=1,nmodes
  tmp1(:,imode)=sixovrpix(imode)*wvol(:,imode)
END DO ! loop over modes
CALL cubrt_v(nmodes*nbox, tmp1, wetdp)
DO imode=1,nmodes
  IF (.NOT. mode(imode)) THEN
    wetdp(:,imode)=drydp(:,imode)
  END IF  ! if(.not. mode(imode))
END DO ! loop over modes

DO imode=1,nmodes
  IF ( ( MINVAL(wetdp (:,imode)) <= 0.0 )   .OR.                               &
       ( MINVAL(drydp (:,imode)) <= 0.0 )   .OR.                               &
       ( MINVAL(wvol  (:,imode)) <= 0.0 )   .OR.                               &
       ( MINVAL(dvol  (:,imode)) <= 0.0 )   .OR.                               &
       ( MINVAL(rhopar(:,imode)) <= 0.0 ) ) THEN

    WRITE(umMessage,'(A)') 'Error in VOLUME_MODE'
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A)') 'Minval <= 0'
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A)') 'wetdp/drydp/wvol/dvol/rhopar <= 0'
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'imode : ', imode
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MINVAL(mdwat(:,imode)) : ',                  &
                                  MINVAL(mdwat(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MAXVAL(mdwat(:,imode)) : ',                  &
                                  MAXVAL(mdwat(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'SUM(mdwat(:,imode)) : ', SUM(mdwat(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MINLOC(mdwat(:,imode)) : ',                     &
                               MINLOC(mdwat(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MAXLOC(mdwat(:,imode)) : ',                     &
                               MAXLOC(mdwat(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MINVAL(wetdp(:,imode)) : ',                  &
                                  MINVAL(wetdp(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MAXVAL(wetdp(:,imode)) : ',                  &
                                  MAXVAL(wetdp(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'SUM(wetdp(:,imode)) : ', SUM(wetdp(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MINLOC(wetdp(:,imode)) : ',                     &
                               MINLOC(wetdp(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MAXLOC(wetdp(:,imode)) : ',                     &
                               MAXLOC(wetdp(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MINVAL(drydp(:,imode)) : ',                  &
                                  MINVAL(drydp(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MAXVAL(drydp(:,imode)) : ',                  &
                                  MAXVAL(drydp(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'SUM(drydp(:,imode)) : ', SUM(drydp(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MINLOC(drydp(:,imode)) : ',                     &
                               MINLOC(drydp(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MAXLOC(drydp(:,imode)) : ',                     &
                               MAXLOC(drydp(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MINVAL(rhopar(:,imode)) : ',                 &
                                  MINVAL(rhopar(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MAXVAL(rhopar(:,imode)) : ',                 &
                                  MAXVAL(rhopar(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'SUM(rhopar(:,imode)) : ', SUM(rhopar(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MINLOC(rhopar(:,imode)) : ',                    &
                               MINLOC(rhopar(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MAXLOC(rhopar(:,imode)) : ',                    &
                               MAXLOC(rhopar(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MINVAL(dvol(:,imode)) : ',                   &
                                  MINVAL(dvol(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MAXVAL(dvol(:,imode)) : ',                   &
                                  MAXVAL(dvol(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'SUM(dvol(:,imode)) : ', SUM(dvol(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MINLOC(dvol(:,imode)) : ',                      &
        MINLOC(dvol(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MAXLOC(dvol(:,imode)) : ',                      &
        MAXLOC(dvol(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MINVAL(wvol(:,imode)) : ',                   &
                                  MINVAL(wvol(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'MAXVAL(wvol(:,imode)) : ',                   &
                                  MAXVAL(wvol(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E15.6)') 'SUM(wvol(:,imode)) : ', SUM(wvol(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MINLOC(wvol (:,imode)) : ',                     &
        MINLOC(wvol (:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'MAXLOC(wvol (:,imode)) : ',                     &
        MAXLOC(wvol (:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    jjl=MINLOC(wvol (:,imode), DIM=1)
    DO iimode=1,nmodes
      IF (mode(iimode)) THEN
        WRITE(umMessage,'(A,I0,2E15.6)') 'iimode, nd(jjl(1),iimode), ' //      &
                                         'mdt(jjl(1),iimode) : ',              &
                                          iimode, nd(jjl(1),iimode),           &
                                          mdt(jjl(1),iimode)
        CALL umPrint(umMessage,src=RoutineName)
        DO iicp=1,ncp
          IF (component(iimode,iicp)) THEN
            WRITE(umMessage,'(A,I0,4X,I0,E15.6)') 'iimode, iicp, MD : ',       &
                                                   iimode, iicp,               &
                                                   md(jjl(1),iimode,iicp)
            CALL umPrint(umMessage,src=RoutineName)
          END IF ! if component present in mode (iicp)
        END DO ! loop over cpts (iicp)
      END IF ! if mode present (iimode)
    END DO ! loop over modes (iimode)

    errcode=1
    WRITE(cmessage,'(5(A,E15.6,A,I0),A)')                                      &
         'Minimum value of "wetdp" or "drydp" or "wvol"'                       &
         //newline//                                                           &
         ' or "wvol" or "dvol" or "rhopar"  < 0.0 '                            &
         //newline//                                                           &
         'MINVAL(wetdp(:,imode)) : ',MINVAL(wetdp(:,imode)),                   &
         ' at location : ', MINLOC(wetdp(:,imode), DIM=1),                     &
         newline//                                                             &
         'MINVAL(drydp(:,imode)) : ',MINVAL(drydp(:,imode)),                   &
         ' at location : ', MINLOC(drydp(:,imode), DIM=1),                     &
         newline//                                                             &
         'MINVAL(wvol(:,imode)) : ',MINVAL(wvol(:,imode)),                     &
         ' at location : ', MINLOC(wvol(:,imode), DIM=1),                      &
         newline//                                                             &
         'MINVAL(dvol(:,imode)) : ',MINVAL(dvol(:,imode)),                     &
         ' at location : ', MINLOC(dvol(:,imode), DIM=1),                      &
         newline//                                                             &
         'MINVAL(rhopar(:,imode)) : ',MINVAL(rhopar(:,imode)),                 &
         ' at location : ', MINLOC(rhopar(:,imode), DIM=1),                    &
         newline//                                                             &
         'Further debugging output in stdout.'
    CALL ereport(RoutineName,errcode,cmessage)

  END IF ! if wetdp/drydp/wvol/dvol/rhopar <= 0
END DO ! loop over modes

IF (glomap_config%l_fix_neg_pvol_wat .OR. l_glomap_clim_radaer) THEN
  DO imode=1,nmodes
    IF ( ( MINVAL(pvol_wat(:,imode)) < 0.0 )   .OR.                            &
         ( MINVAL(mdwat(:,imode))    < 0.0 ) ) THEN
      errcode=1
      WRITE(cmessage,'(2(A,E15.6,A,I0))')                                      &
           'Minimum value of "pvol_wat" or "mdwat" < 0.0 '                     &
           //newline//                                                         &
           'MINVAL(mdwat(:,imode)) : ', MINVAL(mdwat(:,imode)),                &
           ' at location : ', MINLOC(mdwat(:,imode), DIM=1),                   &
           newline//                                                           &
           'MINVAL(pvol_wat(:,imode)) : ', MINVAL(pvol_wat(:,imode)),          &
           ' at location : ', MINLOC(pvol_wat(:,imode), DIM=1)
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_volume_mode
END MODULE ukca_volume_mode_mod
