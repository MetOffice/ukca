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
!    Calculate the H2SO4-H2O binary nucleation rate of sulfate particles.
!    The parameterizations are based on a variety of sources with boundary
!     layer nucleation method selected from the UKCA namelist.
!     Method: 0) Pandis et al (1994), JGR, vol.99, no. D8, pp 16,945-16,957.
!               [Unsupported now]
!             1) Kulmala et al (1998), JGR, vol.103, no. D7, pp 8,301-8,307.
!                [Recommended]
!             2) Vehkamaki et at (2002), JGR, vol.107, no.D22, 4622,
!                doi:10.1029/2002JD002184, 2002
!                [Available for testing]
!     Kulmala/Vehkamaki methodology can be changed by setting i_bhn_method
!      in the code below.
!    The BLN parameterisation is selected from ibln which is derived from
!     the UKCA namelist (i_mode_bln_param_method):
!             1) Activation - Spracklen et al (2008)
!             2) Kinetic    - Spracklen et al (2010)
!             3) Organically mediated - Metzger (2010)
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
MODULE ukca_calcnucrate_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CALCNUCRATE_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE ukca_calcnucrate(nbox, dtz, t, s, rh, aird, h2so4,                  &
                            delh2so4_nucl, sec_org, bln_on,                    &
                            ibln ,i_nuc_method, height, htpblg, s_cond_s)
!---------------------------------------------------------------
!
! Purpose
! -------
! Calculate the H2SO4-H2O binary nucleation rate of sulfate particles.
! Parameterization based on
! Method 1) Pandis et al (1994), JGR, vol.99, no. D8, pp 16,945-16,957.
!        2) Kulmala et al (1998), JGR, vol.103, no. D7, pp 8,301-8,307.
!        3) Vehkamaki et at (2002), JGR, vol.107, no. D22, 4622, doi:
!           10.1029/2002JD002184
!
! Parameters
! ----------
! None
!
! Inputs:
! ------
! NBOX    : Number of grid boxes
! DTZ     : Timestep for nucl/cond competition (s)
! T       : Mid-level temperature (K)
! S       : Specific humidity (kg/kg)
! RH      : Relative humidity (dimensionless 0-1)
! AIRD    : Number density of air (per cc)
! H2SO4   : H2SO4   number density (cm-3)
! Sec_Org : Sec_Org number density (cm-3)
! BLN_ON  : Switch for whether BLN is on (1) or off (0)
! IBLN    : Switch for whether BLN is Activation (1), Kinetic (2), or
!           Organically mediated (Metzger) (3)
! I_NUC_METHOD: Switch for nucleation (how to combine BHN and BLN)
! (1=initial Pandis94 approach (no BLN even if switched on) -- Do not use!!
! (2=binary homogeneous nucleation applying BLN to BL only if switched on)
!   note there is an additional switch i_bhn_method (local to CALCNUCRATE)
!   to switch between using Kulmala98 or Vehkamaki02 for BHN rate
! (3=use same rate at all levels either activation(IBLN=1), kinetic(IBLN=2),
!  Metzger(IBLN=3)
!  Note that if I_NUC_METHOD=3 and IBLN=3 then also add on BHN rate
!   as in Metzger.
! HEIGHT  : Mid-level height of gridbox
! HTPBLG  : Height of boundary-layer in gridbox vertical-column
! S_COND_S :  Condensation sink
!
! Outputs:
! -------
! H2SO4      : H2SO4 number density (cm-3)  (INOUT)
! DELH2SO4_NUCL : Change in H2SO4 due to nucleation (molecules/DTZ)
!
!
! Local Variables
! ---------------
! JRATE      : H2SO4 depletion rate by nucleation (/cc/s)
! CORT      : Air temperature corrected to lie in range (K)
! CORRH     : Relative humidity corr. to lie in range (0-1)
! JKUL      : Kulmala H2SO4/H2O nucleation rate (cm^-3 s^-1)
! JBLN      : BL nucleation rate (cm^-3 s^-1)
! H2SO4OLD2 : Old number density of H2SO4(g) (cm^-3)
! ACONS     : Constant in calculation of JKUL
! BCONS     : Constant in calculation of JKUL
! DELTA     : Parameter in calculation of JKUL
! A2        : Parameter in calculation of JKUL
! B2        : Parameter in calculation of JKUL
! C2        : Parameter in calculation of JKUL
! D2        : Parameter in calculation of JKUL
! E2        : Parameter in calculation of JKUL
! NCRIT     : H2SO4 no. dens. giving JKUL=1 cm^-3s^-1 (cm^-3)
! NSULF     : ln(H2SO4/NCRIT)
! PP        : partial pressure of H2SO4 (Pa)
! SVP       : saturation vapour pressure of H2SO4 (Pa)
! RAC       : relative acidity = PP/SVP
! NWP       : water vapour concentration (cm-3)
! XAL       : H2SO4 mole fraction
! THETA     : =chi in Kulmala et al (98) eq 20 = ln(JKUL)
! L1,L2     : Logicals for BLN methods and switches
! i_bhn_method : method for binary nucleation (Kulmala or Vehkamaki)
!
! Inputted by module CHEMISTRY_CONSTANTS,UKCA_CONSTANTS
! ---------------------------------
! BOLTZMANN  : Boltzman Constant (kg m2 s-2 K-1 molec-1)
! NMOL    : Number of molecules per particle at nucleation
! CONC_EPS: Threshold for concentration (molecules per cc)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! CP_SU   : Index of component in which sulfate is stored
!
! References
! ----------
! * Pandis et al (1994) "The relationship between
!   DMS flux and CCN concentration in remote marine
!   regions", JGR, vol. 99, no. D8, pp. 16,945-16,957.
!
! * Kulmala et al (1998) "Parameterizations for
!   sulfuric acid/water nucleation rates",
!   JGR, vol. 103, no. D7, pp. 8301-8307.
!
! * Vehkamaki et al (2002) "An improved parameterization for
!   sulfuric acid-water nucleation rates for tropospheric
!   and stratospheric conditions!, JGR, vol 107, No D22
!   4622, doi:10.1029/2002JD002184
!
! * Metzger et al (2010) "Evidence for the role of organics in aerosol
!     particle formation under atmospheric conditions", Proc. Nat.
!     Acad. Sci., 107, 6646-6651.
!
! * Spacklen et al (2008) "Contribution of particle formation to cloud
!     condensation nuclei concentrations", GRL, 35, MAR 29 2008.
!
! * Spacklen et al (2010) "Explaining global surface aerosol number
!     concentrations in terms of primary emissions and particle
!     formation, ACP, 10:4775-4793.
!
!----------------------------------------------------------------------
USE ukca_config_constants_mod, ONLY: boltzmann
USE ukca_constants,     ONLY: nmol, conc_eps
USE ukca_binapara_mod,  ONLY: ukca_binapara
USE yomhook,            ONLY: lhook, dr_hook
USE parkind1,           ONLY: jprb, jpim
USE ereport_mod,        ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!     Arguments
INTEGER, INTENT(IN) :: nbox             ! No. of boxes
INTEGER, INTENT(IN) :: bln_on           ! switch: BLN is on (1) or off (0)
INTEGER, INTENT(IN) :: ibln             ! Switch for whether BLN is
                                        ! Activation(1), Kinetic(2),
                                        ! or Metzger(3)
                                        ! IBLN is derived from namelist
                                        ! variable I_MODE_BLN_PARAM_METHOD
INTEGER, INTENT(IN) :: i_nuc_method     ! Switch for nucleation
                                        ! (how to combine BHN and BLN)

REAL, INTENT(IN)    :: sec_org(nbox)    ! Sec_Org number density (cm-3)
REAL, INTENT(IN)    :: dtz              ! timestep for nucl/condens (s)
REAL, INTENT(IN)    :: t(nbox)          ! Temperature (K)
REAL, INTENT(IN)    :: s(nbox)          ! Specific humidity (kg/kg)
REAL, INTENT(IN)    :: rh(nbox)         ! Relative humidity fraction
REAL, INTENT(IN)    :: aird(nbox)       ! Number density of air (per cc)
REAL, INTENT(IN)    :: height(nbox)     ! Mid-level height of gridbox (m)
REAL, INTENT(IN)    :: htpblg(nbox)     ! Boundary layer height (m)
REAL, INTENT(IN)    :: s_cond_s(nbox)   ! Condensation sink

REAL, INTENT(IN OUT) :: h2so4(nbox)      ! H2SO4 number density (cm-3)
REAL, INTENT(OUT)   :: delh2so4_nucl(nbox) ! Change in H2SO4 due to
                                           ! nucleation (molecules/DTZ)

!     Local variables
CHARACTER(LEN=errormessagelength) :: cmessage     ! error message
INTEGER, PARAMETER :: i_bhn_method_kulmala   = 1
INTEGER, PARAMETER :: i_bhn_method_vekhamaki = 2
INTEGER :: i_bhn_method                 ! Set to one of the above
INTEGER :: errcode                      ! Variable passed to ereport
INTEGER :: jl                           ! loop counter

REAL, PARAMETER :: zmaxbln = 6000.0     ! Max height (m) for BL nucleation

REAL, PARAMETER :: afac_pna = 5.0e-13   ! Metzger et al (2010)
REAL, PARAMETER :: afac_act = 5.0e-7    ! Activation coefficient
REAL, PARAMETER :: afac_kin = 4.0e-13   ! Kinetic coefficient

REAL :: jrate(nbox)     ! H2SO4 depletion rate by nucleation (/cc/s)
REAL :: cort(nbox)      ! Corrected temperature
REAL :: corrh(nbox)     ! Corrected RH
REAL :: h2so4old2       ! Initial h2so4 number concentration
REAL :: acons           ! Constant for jkul calculation
REAL :: bcons           ! Constant for jkul calculation
REAL :: delta           ! Intermediate results used in Kulmala
REAL :: a2              !        "
REAL :: b2              !        "
REAL :: c2              !        "
REAL :: d2              !        "
REAL :: e2              !        "
REAL :: ncrit           ! H2SO4 no. dens. giving JKUL=1 cm^-3s^-1 (cm^-3)
REAL :: nsulf           ! ln(H2SO4/NCRIT)
REAL :: pp              ! partial pressure of H2SO4 (Pa)
REAL :: svp             ! saturation vapour pressure of H2SO4 (Pa)
REAL :: nwp             ! water vapour concentration (cm-3)
REAL :: xal             ! H2SO4 mole fraction
REAL :: rac             ! relative acidity = PP/SVP
REAL :: theta           ! =chi in Kulmala et al (98) eq 20 = ln(JKUL)
REAL :: jkul            ! Kulmala H2SO4/H2O nucleation rate (cm^-3 s^-1)
REAL :: jbln            ! BL nucleation rate (cm^-3 s^-1)
REAL :: japp_bln        ! Binary nucleation rate
REAL :: dpbln           ! Diameter (nm) for calculating rate of BLN
LOGICAL :: l1           ! (i_nuc_method == 2 & height > zBL) or bln_on == 0
LOGICAL :: l2           ! i_nuc_method == 3 & ibln == 3

REAL :: jveh(nbox)      ! Vehkamaki nucleation rate (cm^-3 s^-1)
REAL :: rc(nbox)        ! Radius of the critical cluster in nm
REAL :: japp(nbox)      ! Binary nucleation rate
REAL :: zbl(nbox)       ! Boundary layer height from htpblg, but filtered
                        ! by MIN(htpblg,zmaxbln)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALCNUCRATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     i_bhn_method = i_bhn_method_kulmala   ! set for Kulmala
i_bhn_method = i_bhn_method_vekhamaki ! set for Vehkamakki

! Filter the BL heights to exclude those > zmaxbln
DO jl=1,nbox
  zbl(jl) = MIN(htpblg(jl),zmaxbln)
END DO

IF ((ibln >= 1) .AND. (ibln <= 3)) dpbln=1.5 ! Metzger  (2010,pnas)

! .. local copies of rh and t with bounds
corrh(:) = rh(:)
cort(:)  = t(:)
WHERE (corrh < 0.1) corrh = 0.1
WHERE (corrh > 0.9) corrh = 0.9
WHERE (cort < 233.0) cort = 233.0
WHERE (cort > 298.0) cort = 298.0

! Initialise JRATE and DELH2SO4_NUCL to zero
jrate(:) = 0.0
delh2so4_nucl(:) = 0.0

! I_NUC_METHOD must be either 2 or 3
! IBLN must be in range [1-3]
IF ( ( i_nuc_method < 2 .OR. i_nuc_method > 3) .OR.                            &
     ( ibln < 1 .OR. ibln > 3 ) ) THEN
  cmessage='I_NUC_METHOD (2 - 3) or IBLN (1 - 3) '//'out of valid range'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF ((i_nuc_method == 2) .OR. (i_nuc_method == 3)) THEN

  IF (i_bhn_method == i_bhn_method_vekhamaki) THEN ! use Vehkamaki for BHN
    CALL ukca_binapara(nbox, t, rh, h2so4, jveh, rc)
    !      Returns nucleation rate JVEH at critical-cluster-radius RC
    !        using Vehkamaki et al (2002)
  END IF

  l2 = (i_nuc_method == 3 .AND. ibln == 3)
  ! .. L2 is logical for having chosen BLN throughout the column but
  ! ..   also have chosen the Metzger approach which also includes
  ! ..   Kulmala nucleation rate at all levels

  DO jl = 1,nbox

    l1 = (i_nuc_method == 2 .AND. (height(jl) > zbl(jl) .OR. bln_on == 0))
    ! .. L1 is logical for using BHN when BLN is off or BLN is on
    ! ..   and the box is above the height of the BL

    IF (l1 .OR. l2) THEN ! include Kulmala in JRATE and DELH2SO4_NUCL

      IF (i_bhn_method == i_bhn_method_kulmala) THEN
        !       Calculate H2SO4-H2O binary nucleation rate using Kulmala(98)
        IF (h2so4(jl) > conc_eps) THEN
          !          Calculate parametrised variables
          delta = 1.0+(cort(jl)-273.15)/273.15
          a2 = 25.1289-(4890.8/cort(jl))-2.2479*delta*corrh(jl)
          b2 = (7643.4/cort(jl))-1.9712*delta/corrh(jl)
          c2 = -1743.3/cort(jl)
          !          Calculate H2SO4 needed to get J=1cm-3 s-1
          ncrit = EXP(-14.5125 + 0.1335*cort(jl) -                             &
                      10.5462*corrh(jl) + 1958.4*corrh(jl)/cort(jl))
          nsulf = LOG(h2so4(jl)/ncrit)
          !          Calculate saturation vapour pressure of H2SO4 (Pa)
          svp = EXP(27.78492066 - 10156.0/t(jl))
          !          Calculate water vapour concentration (cm-3)
          nwp = aird(jl)*1.609*s(jl)
          !          Calculate partial pressure of H2SO4 (Pa)
          pp = h2so4(jl)*1.0e6*boltzmann*t(jl)
          !          Calculate relative acidity
          rac = pp/svp
          IF (nwp > 0.0) THEN
            d2 = 1.2233 - (0.0154*rac/(rac+corrh(jl)))-                        &
                   0.0415*LOG(nwp) + 0.0016*cort(jl)
          ELSE
            d2 = 1.2233 - (0.0154*rac/(rac+corrh(jl)))                         &
                    + 0.0016*cort(jl)
          END IF
          e2 = 0.0102
          xal = d2 + e2*LOG(h2so4(jl))
          ! .. NSULF is ln(H2SO4/NCRIT)
          theta = a2*nsulf+b2*xal+c2
          ! .. JKUL is nucleation rate in particles/cc/s
          jkul = EXP(theta)

          !          Calculate new H2SO4
          IF (jkul > 1.0e-3) THEN
            acons = (EXP(b2*d2+c2))*(1.0/ncrit)**a2
            bcons = a2 + b2*e2
            !           Calculate new H2S04 concentrations
            h2so4old2 = h2so4(jl)
            h2so4(jl) = ( h2so4(jl)**(1.0-bcons) +                             &
                      (bcons - 1.0)*nmol*acons*dtz)                            &
                      **((1.0 / (1.0 - bcons)))
            !           Calculate change in sulfuric acid
            IF (h2so4(jl) > h2so4old2) h2so4(jl) = h2so4old2
            IF (h2so4(jl) < 0.0) h2so4(jl) = 0.0
            delh2so4_nucl(jl) = delh2so4_nucl(jl) + h2so4old2 - h2so4(jl)
            !           Calculate depletion rate of H2SO4 (molecules/cc/s)
            jrate(jl) = jrate(jl) + nmol*acons*h2so4old2**bcons
          END IF ! JKUL > 1.0E-3
        END IF ! H2SO4 > CONC_EPS

      ELSE IF (i_bhn_method == i_bhn_method_vekhamaki) THEN
        !         Calculate H2SO4-H2O binary nucleation rate using Vehkamaki

        IF ( (h2so4(jl) > conc_eps) .AND. (jveh(jl) > 0.0)                     &
           .AND. (s_cond_s(jl) > 0.0) .AND. (rc(jl) > 0.0) ) THEN
          !
          japp(jl) = jveh(jl)*EXP(0.23*(1.0/3.0 - 1.0/(2*rc(jl)))*             &
                              s_cond_s(jl)/(1.0e6*h2so4(jl)*1.0e-13))
          !
          !          Calculate new H2SO4 and particle number density
          IF (japp(jl) > 1.0e-3) THEN
            jrate(jl) = jrate(jl)+japp(jl)
            !
            !           Calculate new H2S04 concentrations
            h2so4old2 = h2so4(jl)
            h2so4(jl) = h2so4(jl)-jrate(jl)*nmol*dtz
            !
            !  Calculate change in sulfuric acid
            IF (h2so4(jl) > h2so4old2) h2so4(jl) = h2so4old2
            IF (h2so4(jl) < 0.0) h2so4(jl) = 0.0
            !
            delh2so4_nucl(jl) = delh2so4_nucl(jl) + h2so4old2 - h2so4(jl)
            !
          END IF ! if JAPP > 1.0E-3
        END IF ! if H2SO4 > 0 & JVEH > 0 & S_COND_S > 0 & RC > 0
        !
      ELSE
        ! Invalid i_bhn_method
        cmessage = ' I_BHN_METHOD is incorrect'
        errcode = 1
        CALL ereport(RoutineName,errcode,cmessage)
      END IF ! if i_bhn_method=....

    END IF ! if L1 or L2
    !
    IF (.NOT. l1) THEN    ! include BLN in JRATE and DELH2SO4_NUCL
      ! .. .NOT.L1 means either I_NUC_METHOD=3 or I_NUC_METHOD=2 and BLN_ON=1
      !     and in BL
      IF (h2so4(jl) > conc_eps) THEN
        ! .. JBLN is the nucleation rate at diameter DPBLN (nm) in cm^-3 s^-1
        IF (ibln == 1) THEN
          jbln = afac_act*h2so4(jl)             ! activation
        ELSE IF (ibln == 2) THEN
          jbln = afac_kin*h2so4(jl)*h2so4(jl)   ! kinetic
        ELSE IF (ibln == 3) THEN
          jbln = afac_pna*h2so4(jl)*sec_org(jl) ! Metzger
        ELSE
          cmessage = ' Value of IBLN outside range'
          errcode = 1
          CALL ereport(RoutineName,errcode,cmessage)
        END IF

        japp_bln = jbln*EXP(0.23*(1.0/3.0 - 1.0/dpbln)*                        &
                   s_cond_s(jl)/(1.0e6*h2so4(jl)*1.0e-13))

        IF (japp_bln > 1.0e-3) THEN
          jrate(jl) = jrate(jl) + japp_bln
          h2so4old2 = h2so4(jl)
          h2so4(jl) = h2so4(jl) - japp_bln*nmol*dtz
          ! Calculate change in sulfuric acid
          IF (h2so4(jl) > h2so4old2) h2so4(jl) = h2so4old2
          IF (h2so4(jl) < 0.0) h2so4(jl) = 0.0
          delh2so4_nucl(jl) = delh2so4_nucl(jl) + h2so4old2 - h2so4(jl)
        END IF    ! japp_bln > 1.0E-3
      END IF    ! h2so4 > conc_eps
    END IF    ! method=2 and in BL with bln_on=1 or method=3

  END DO    ! Loop over NBOX
END IF     ! METHOD 2 OR 3

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calcnucrate
END MODULE ukca_calcnucrate_mod
