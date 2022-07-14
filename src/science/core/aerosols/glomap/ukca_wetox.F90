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
!    Calculates conversion of dissolved SO2 to sulfate by reaction with
!    H2O2 and O3. Uses Henry's law constants & oxidation rate constants
!    from Seinfeld & Pandis "Atmospheric Chemistry & Physics"
!    Note that this code is for use in the CTM and box model framework
!    or where the SO2 wet oxidation is to be done in the aerosol
!    code (WETOX_IN_AER=1).
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
MODULE ukca_wetox_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_WETOX_MOD'

CONTAINS

SUBROUTINE ukca_wetox(nbox,nd,drydp,delso2,delso2_2,                           &
 delh2o2,lowcloud,vfac,so2vmr,h2o2,zh2o2,zo3,zho2,pmid,                        &
 t,aird,s,dtc,lday,lwc,iso2wetoxbyo3)
!----------------------------------------------------------------------
!
!     Purpose
!     -------
!     Calculates conversion of dissolved SO2 to sulfate by reaction with
!     H2O2 and O3. Uses Henry's law constants & oxidation rate constants
!     from Seinfeld & Pandis "Atmospheric Chemistry & Physics"
!     (see pgs 337-367).
!     Assumes cloud LWC=0.2 gm^-3 and pH cloud = 4/5 SO2 dependant.
!
!     Also replenishes H2O2 concentration by HO2 self-reaction according
!     to expression in Jones et al (2001), JGR Atmos. vol. 106, no.D17,
!     pp 20,293--20,310.
!
!     Henry's law constant values at 298K from Seinfeld & Pandis pg 341.
!     Heats of dissolution values for Henry's law constants at other
!     temperatures from Seinfeld & Pandis page 342.
!
!     Note that this code is for use in the CTM and box model framework
!     or where the SO2 wet oxidation is to be done in the aerosol
!     code (WETOX_IN_AER=1).
!
!     For usual UKCA runs in the UM, this code is not called.
!
!     Note also, that this code currently has an assumption that the
!     clouds indicated by LOWCLOUD and VFAC are low clouds --- in fact
!     in the current version a liquid water content of 0.2 g/m3 is set
!     for all such clouds. It is intended to be used when climatological
!     monthly-mean cloud fields from ISCCP satellite database are used.
!
!     Inputs
!     ------
!     NBOX       : Number of grid boxes
!     ND         : Initial no. concentration of aerosol mode (ptcls/cc)
!     DRYDP      : Geometric mean dry diameter for each mode (m)
!     LOWCLOUD   : Horizontal low cloud fraction
!     VFAC       : Vertical low cloud fraction
!     SO2VMR     : Volume mixing ratio of SO2(g)
!     H2O2       : Number density of H2O2(g) [molecules per cc]
!     ZO3        : Background vmr of O3
!     ZH2O2      : Background concentration of H2O2 [molecules per cc]
!     ZHO2       : Background concentration of HO2  [molecules per cc]
!     PMID       : Air pressure (Pa)
!     T          : Air temperature (K)
!     AIRD       : Number density of air molecules
!     S          : Air specific humidity (kg/kg)
!     DTC        : Time step (s)
!     LDAY       : 1-Day, 0-Night
!     LWC        : Cloud liquid water content [kg/m3]
!     ISO2WETOXBYO3: Switch (1/0) for whether to include in-cloud oxdtn
!                    by O3 (pH value prescribed according to SO2 concn)
!
!     Outputs
!     -------
!     DELSO2     : Mass of S(IV) --> S(VI) by H2O2 (molecules/cc/DTC)
!     DELSO2_2   : Mass of S(IV) --> S(VI) by O3 (molecules/cc/DTC)
!     DELH2O2    : Total change in H2O2 mass due to reaction with S(IV)
!                : and replenishment via HO2 (molecules/cc/DTC)
!
!     Local Variables
!     ---------------
!     F          : cloud fraction
!     H2O2VMR    : volume mixing ratio of H2O2(g)
!     PHCLOUD    : pH of cloud (assumed to be 4/5 if SO2>or<0.5ppbv)
!     HSO2_298   : Henry's law constant for SO2 dissol. at 298K (M/atm)
!     DELHSO2_KCAL: Heat of dissolution for SO2 (kCal per mole)
!     DELHSO2_R  : Heat of dissolution for HSO2 / R (K)
!     HH2O2_298  : Henry's law constant for H2O2 dissol. at 298K (M/atm)
!     DELHH2O2_KCAL: Heat of dissolution for H2O2 (kCal per mole)
!     DELHH2O2_R : Heat of dissolution for HH2O2 / R (K)
!     K_S1       : constant in dissociation of SO2.H2O to H+ + HSO3-
!     K_S2       : constant in dissociation of HSO3 to H+ + SO3 2-
!     K1         : constant in equation for S(IV) oxidation
!     K2         : constant in equation for S(IV) oxidation
!     JPERKCAL   : number of Joules per kilocalorie (J kcal^-1)
!     HSO2       : Henry's law constant for SO2 dissol. at T (M/atm)
!     HH2O2      : Henry's law constant for H2O2 dissol. at T (M/atm)
!     HPLUSM     : Concentration of H+ ions in solution (M)
!     HSO3M      : Concentration of HSO3- ions in solution (M)
!     PSO2ATM    : Partial pressure of SO2(g) in atmospheres
!     PH2O2ATM   : Partial pressure of H2O2(g) in atmospheres
!     TERM       : Temporary variable in Henry's law constant calcn
!     SO2        : Number concentration of SO2 (molecules per cc)
!     MAXDSL     : Maximum change in sulfate ( =min{SO2,H2O2} )
!     NWP        : No. concentration of water vapour (molecules per cc)
!     DH2O2TMP   : Rate of re-formation of H2O2 via HO2
!     OLDH2O2    : Temporary variable to store H2O2 concentration
!
!     Inputs from Shared Modules
!     ---------------------------------
!     AVOGADRO   : Avogadro's constant
!     RMOL       : Molar gas constant (J/mol/K)
!     NMODES     : Number of possible aerosol modes
!     MODE       : Defines which modes are set
!     NUM_EPS    : Value of NEWN below which don't apply process
!     MH2O2F     : Index for semi-prognostic H2O2
!--------------------------------------------------------------------
USE chemistry_constants_mod, ONLY: avogadro
USE ukca_um_legacy_mod,   ONLY: rmol, pi
USE ukca_mode_setup,      ONLY: nmodes, mode, num_eps
USE ukca_setup_indices,   ONLY: mh2o2f
USE ukca_config_specification_mod, ONLY: ukca_config
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: lday(nbox)
INTEGER, INTENT(IN) :: iso2wetoxbyo3
REAL, INTENT(IN)    :: nd(nbox,nmodes)
REAL, INTENT(IN)    :: drydp(nbox,nmodes)
REAL, INTENT(IN)    :: lowcloud(nbox)
REAL, INTENT(IN)    :: vfac(nbox)
REAL, INTENT(IN)    :: so2vmr(nbox)
REAL, INTENT(IN)    :: zo3(nbox)
REAL, INTENT(IN)    :: zh2o2(nbox)
REAL, INTENT(IN)    :: zho2(nbox)
REAL, INTENT(IN)    :: pmid(nbox)
REAL, INTENT(IN)    :: t(nbox)
REAL, INTENT(IN)    :: aird(nbox)
REAL, INTENT(IN)    :: s(nbox)
REAL, INTENT(IN)    :: dtc
REAL, INTENT(IN)    :: lwc(nbox)
REAL, INTENT(IN OUT) :: h2o2(nbox)
REAL, INTENT(OUT)   :: delso2(nbox)
REAL, INTENT(OUT)   :: delso2_2(nbox)
REAL, INTENT(OUT)   :: delh2o2(nbox)
!
!     Local variables
INTEGER :: jl
REAL    :: f
REAL    :: phcloud
REAL    :: delhso2_r
REAL    :: delho3_r
REAL    :: delhh2o2_r
REAL    :: hso2
REAL    :: ho3
REAL    :: hh2o2
REAL    :: hplusm
REAL    :: hso3m
REAL    :: o3m
REAL    :: dsivmdt2
REAL    :: h2o2vmr
REAL    :: o3vmr
REAL    :: pso2atm
REAL    :: po3atm
REAL    :: ph2o2atm
REAL    :: term
REAL    :: so2
REAL    :: maxdsl
REAL    :: nwp
REAL    :: dh2o2tmp
REAL    :: oldh2o2
REAL    :: so2vmr2
REAL    :: pso2atm2
REAL    :: so2h2o
REAL    :: so32m
REAL    :: lwc_gm3(nbox)
! .. below are new local variables added
REAL    :: dc
REAL    :: fac_c
REAL    :: dropr_3
REAL    :: c_3
REAL    :: r_3
REAL    :: dropr_4
REAL    :: c_4
REAL    :: r_4
REAL, PARAMETER :: hso2_298     =  1.23
REAL, PARAMETER :: delhso2_kcal = -6.23
REAL, PARAMETER :: ho3_298      =  1.13e-2
REAL, PARAMETER :: delho3_kcal  = -5.04
REAL, PARAMETER :: hh2o2_298    =  7.45e4
REAL, PARAMETER :: delhh2o2_kcal=-14.5
REAL, PARAMETER :: k_s1         =  1.3e-2
REAL, PARAMETER :: k_s2         =  6.6e-8
REAL, PARAMETER :: k0_o3        =  2.4e4
REAL, PARAMETER :: k1_o3        =  3.7e5
REAL, PARAMETER :: k2_o3        =  1.5e9
REAL, PARAMETER :: jperkcal     =  4.187e3

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_WETOX'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!      LWC(:)=0.2e-3 ! set to constant at 0.2 g/m3 [of air] at present
!                      NOW INTENT(IN)
lwc_gm3(:)=lwc(:)*1.0e3
! .. LWC is in kg/m3 [of air], LWC_GM3 is in g/m3 [of air]

DO jl=1,nbox
  !
  oldh2o2=h2o2(jl)
  delso2(jl)=0.0
  delso2_2(jl)=0.0
  !      Calculate in-cloud (wet) oxidation of SO2
  f=lowcloud(jl)*vfac(jl)
  ! .. F is the cloud fraction in the gridbox
  IF (f > 0.0) THEN
    IF (mode(3) .OR. mode(4)) THEN ! if accum-sol or coarse-sol mode
      ! .. below is for if some particles in either mode in this gridbox.
      IF ((nd(jl,3) > num_eps(3)) .OR. (nd(jl,4) > num_eps(4))) THEN
        IF (so2vmr(jl) > 5e-10) phcloud=4.0
        IF (so2vmr(jl) < 5e-10) phcloud=5.0
        ! .. assumes pH of cloud is either 4 or 5 depending on if SO2>/<0.5ppb
        h2o2vmr=h2o2(jl)/aird(jl)
        o3vmr=zo3(jl) ! ZO3 is already vmr

        ! .. below calculates partial pressures of SO2,O3,H2O2 in atmospheres
        pso2atm=so2vmr(jl)*pmid(jl)/1.0e5
        po3atm=o3vmr*pmid(jl)/1.0e5
        ph2o2atm=h2o2vmr*pmid(jl)/1.0e5

        ! .. below calculates heats of dissolution for Henry's law coeffs
        ! .. (values taken from page 342 Seinfeld & Pandis).
        delhso2_r=delhso2_kcal*jperkcal/rmol
        delho3_r=delho3_kcal*jperkcal/rmol
        delhh2o2_r=delhh2o2_kcal*jperkcal/rmol

        ! .. below stores temperature correction for Henry's law coeffs
        term=(1.0/t(jl))-(1.0/298.0)

        ! .. below calculates Henry's law coefficients for SO2,O3,H2O2
        hso2=hso2_298*EXP(-delhso2_r*term)
        ho3=ho3_298*EXP(-delho3_r*term)
        hh2o2=hh2o2_298*EXP(-delhh2o2_r*term)

        ! .. below stores dissolved ion concentrations for H+,HSO3-
        hplusm=10.0**(-phcloud)
        hso3m=k_s1*hso2*pso2atm/hplusm       ! [6.37 S&P pg 349]

        ! .. below calculates dissolved O3 concentration
        o3m=ho3*po3atm

        !---------------------------------------------------------------------
        ! This section below calculates in-cloud oxidation of SO2 by H2O2
        !

        ! set SO2 molecular concentration (per cc)
        so2=aird(jl)*so2vmr(jl)

        dc=1.0e-5
        fac_c=4.0*pi*dc*1.0e6
        ! 1.0E6 converts ND cm^-3 ---> m^-3 ---> ensures R_3 is dimensionless
        !        n.b. DROPR,DC,ND,DTC are in m,m^2s^-1,cm^-3,s

        ! below is for acc-sol mode
        dropr_3=40.0*(0.5*drydp(jl,3)) ! assume cloud radius = 40 * g.m.radius
        c_3=fac_c*dropr_3
        r_3=-nd(jl,3)*c_3*dtc

        ! below is for cor-sol mode
        dropr_4=40.0*(0.5*drydp(jl,4)) ! assume cloud radius = 40 * g.m.radius
        c_4=fac_c*dropr_4
        r_4=-nd(jl,4)*c_4*dtc

        ! below calculates as in GLOMAP-bin -- here include F as cloud fraction
        delso2(jl)=f*(1.0-EXP(r_3+r_4))*so2

        IF (so2 <= h2o2(jl)) maxdsl=so2
        IF (so2 > h2o2(jl)) maxdsl=h2o2(jl)
        IF (delso2(jl) > maxdsl) delso2(jl)=maxdsl
        ! .. above limits DELSO2 to be max(SO2,H2O2)

        ! .. below reduces SO2, H2O2 gas phase according to above rate
        h2o2(jl)=h2o2(jl)-delso2(jl)
        so2vmr2=so2vmr(jl)-delso2(jl)/aird(jl)

        IF (iso2wetoxbyo3 == 1) THEN ! switch for SO2 wetox by O3
          !---------------------------------------------------------------------
          ! This section below then calculates in-cloud oxidation of SO2 by O3

          ! .. below calculates new partial pressure of SO2 in atmospheres
          pso2atm2=so2vmr2*pmid(jl)/1.0e5

          ! .. below calculates new dissolved concentration of SO2
          so2h2o=hso2*pso2atm2

          ! .. below then calculates new dissolved ions conc of HSO3^-
          !                                                     & SO3^{2-}
          hso3m=k_s1*so2h2o/hplusm
          so32m=k_s2*hso3m/hplusm

          ! .. below calculates oxidation rate of S(IV) by O3 (M/s)
          ! .. (see Seinfeld & Pandis, page 363, equation 6.80)
          dsivmdt2=(k0_o3*so2h2o+k1_o3*hso3m+k2_o3*so32m)*o3m

          ! .. below then converts this rate to a change in concentration
          ! .. in SO2 in the gas phase over the chemistry time step DTC.
          delso2_2(jl)=dsivmdt2*f*(lwc_gm3(jl)*1.0e-6)*avogadro*dtc*1.0e-3
          ! .. 0.001 per liter = 1 per cc
          ! .. 1.0e-6 converts from g/m3 [of air] to g/cm3 [of air]
          ! .. 1.0e-3 is 1/rho_water where rho_water = 1000 g/dm3 (=1000kg/m3)

          IF (so2vmr2 <= zo3(jl)) maxdsl=so2vmr2*aird(jl)
          IF (so2vmr2 > zo3(jl)) maxdsl=zo3(jl)*aird(jl)
          IF (delso2_2(jl) > maxdsl) delso2_2(jl)=maxdsl
          ! .. above limits DELSO2_2 to available SO2,O3

          ! .. Note that SO2 does not need to be updated again here as
          ! .. we are not updating the main SO2 concentrations in this
          ! .. subroutine.
          ! .. This is done using the returned arrays DELSO2 and DELSO2_2 in the
          ! .. UKCA_AERO_STEP routine.

        END IF ! if switch for SO2 in-cloud oxidation by O3 = 1

      END IF ! if some particles in either soluble accum or coarse
    END IF ! if either soluble accum or coarse are defined
  END IF ! if grid box is in cloud

  !---------------------------------------------------------------------
  ! This section calculates replenishment of H2O2 by HO2 self-reaction.
  ! (only does this for off-line version with semi-prognostic H2O2)

  IF (mh2o2f > 0 .AND. (.NOT. ukca_config%l_ukca_chem)) THEN
                                                ! if using semi-prognostic H2O2
    IF (lday(jl) == 1) THEN ! if daytime then H2O2 is reformed

      !        Below calculates water vapour concentration (cm-3) NWP
      !            (1.609*S(JL)) is vmr of H2O (g)
      nwp = aird(jl)*(1.609*s(jl))

      !        Calculate rate of H2O2 replenishment according to background
      !        HO2 concentration ZHO2 (follows expression in Jones et al 2001)
      dh2o2tmp = (2.3e-13*EXP(600.0/t(jl)) +                                   &
          1.9e-33*aird(jl)*EXP(890.0/t(jl))) *                                 &
          (1.0 + 1.4e-21*nwp*EXP(2200.0/t(jl)))*                               &
          zho2(jl)*zho2(jl)

      h2o2(jl)=h2o2(jl)+dh2o2tmp*dtc
      ! .. above replenishes H2O2 due to production by HO2 self-reaction
      IF (h2o2(jl) > zh2o2(jl)) h2o2(jl)=zh2o2(jl)
      ! .. above limits H2O2 to be at most ZH2O2 (background concentration)
    END IF ! if day-time
  END IF ! if MH2O2F>0 (using semi-prognostic H2O2)

  !! n.b. the line below has been moved here cf v1_gm4_coupled (bug-fix).
  delh2o2(jl)=h2o2(jl)-oldh2o2 ! calculate net change in H2O2

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_wetox
END MODULE ukca_wetox_mod
