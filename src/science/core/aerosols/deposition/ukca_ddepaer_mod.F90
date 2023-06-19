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
!    Calculates aerosol dry deposition (and sedimentation).
!    Based on the parameterisation of Zhang et al (2001) which
!    uses the method in the model of Slinn (1982).
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
MODULE ukca_ddepaer_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_DDEPAER_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE ukca_ddepaer (nbox,nbudaer,nd,md,mdt,                               &
                         rhopar,znot,                                          &
                         dtc,wetdp,ustr,pmid,pupper,plower,t,surtp,            &
                         rhoa,mfpa,dvisc,bud_aer_mas,ilscat)
!
! Purpose
! -------
! Calculates aerosol dry deposition (and sedimentation).
! Based on the parameterisation of Zhang et al (2001) which
! uses the method in the model of Slinn (1982).
!
! Evaluate deposition velocity in lowest level as:
!
! V_dep = V_g + 1/(A_r + S_r)
!
! where V_dep is the deposition velocity
!       V_g   is the gravitational velocity = rho_p*Dp^2*g*CF/(18*DVISC)
!       CF    is the Cunningham slip correction
!       DVISC is the dynamic viscosity
!       A_r   is the aerodynamic resitance
!       S_r   is the surface resistance
!       Dp    is the particle diameter
!       rho_p is the particle density
!       g     is the gravitational acceleration
!
! Evaluate S_r=1/{ 3 * ustar * (EB + EIM + EIN) }
!
! following parameterization by Zhang et al (2001) where
!
! EB,EIM,EIN are collection efficiencies for Brownian diffusion,
! impaction and interception respectively.
!
! EB = Sc^-YR where Sc is the particle Schmidt number = nu/D
!                                where nu = kinematic viscosity of air
!                                      D =  particle diffusion coeff.
!
!  and YR is surface-dependent constant, values as in Table 3 (Zhang01)
!         0.50 over water          (Land use category 13-14)
!         0.56 over forest         (Land use category  1- 5)
!         0.54 over grass and ice  (Land use category  6,12)
!
! EIM = { St/(ALPHA+St) }^2
!
!    where St is the Stokes number = V_g * ustar^2 / DVISC  (z0<1mm)
!                                  = V_g * ustar   / (g*CR) (z0>1mm)
!
!                                    [smooth & rough flow regimes]
!
!      and ALPHA,CR are surface-dependent constant, values as Table 3:
!         ALPHA=100.0, CR=0.0 over water [only divide by CR for veg]
!         ALPHA= 50.0, CR=0.0 over ice   [only divide by CR for veg]
!         ALPHA=  1.0, CR=0.005 over grass
!         ALPHA=  1.2, CR=0.002 over forest
!
! EIN = 0.5*Dp/CR
!
! Evaluates drydep & sedimentation for number & mass using 0th & 3rd
! order moment specific coefficients for modal aerosol as in Appendix 4
! Binkowski & Shankar (1995) JGR, vol 100, no D12, pp. 26,191--26,209.
!
! Note --- only evaluates sedimentation at lowest gridbox if this
!          routine is used --- sedimentation at higher levels neglected.
!
! Inputs :
! ------
! NBOX      : Number of grid boxes
! nbudaer   : Number of aerosol budget fields
! ND        : Initial no. concentration of aerosol mode (ptcls/cc)
! MD        : Avg cpt mass of aerosol ptcl in size mode (particle^-1)
! MDT       : Avg tot mass of aerosol ptcl in size mode (particle^-1)
! RHOPAR    : Particle density [incl. H2O & insoluble cpts] (kgm^-3)
! ZNOT      : Roughness length (m)
! DTC       : Chemical timestep (s)
! WETDP     : Wet diameter for ptcl with dry diameter DRYDP (m)
! USTR      : Friction velocity(ms-1)
! PMID      : Centre level pressure (Pa)
! PUPPER    : Pressure at box upper interface (Pa)
! PLOWER    : Pressure at box lower interface (Pa)
! T         : Centre level temperature (K)
! SURTP     : Surface type [0=sea-surf,1=land-surf,2=above-surf]
! RHOA      : Air density (kg/m3)
! MFPA      : Mean free path of air (m)
! DVISC     : Dynamic viscosity of air (kg m-1 s-1)
! ILSCAT    : Land surface category (based on 9 landsurf types)
!
! Outputs
! -------
! Updated particle number density ND (/cm3)
! Updated particle avg mass MD (molecules/particle)
! Updated Avg tot mass of aerosol ptcl in size mode MDT (particle^-1)
! BUD_AER_MAS : Aerosol mass budgets (mlcls/cc/tstep)
!
! Local Variables
! ---------------
! PS_AV_0    : 0th moment avg particle Schmidt Number
! PS_AV_3    : 3rd moment avg particle Schmidt Number
! KVISC      : Kinematic viscosity of air (m2 s-1)
! VGRAV_AV_0 : 0th moment avg grav. settling vel. (m/s)
! VGRAV_AV_3 : 3rd moment avg grav. settling vel. (m/s)
! VDEP_AV_0  : 0th moment avg deposition velocity (m/s)
! VDEP_AV_3  : 3rd moment avg deposition velocity (m/s)
! DCOEF_AV_0 : 0th moment avg particle diffusion coefficient(m2/s)
! DCOEF_AV_3 : 3rd moment avg particle diffusion coefficient(m2/s)
! SN_AV_0    : 0th moment avg Stokes number
! SN_AV_3    : 3rd moment avg Stokes number
! SR_AV_0    : 0th moment avg surface resistance
! SR_AV_3    : 3rd moment avg surface resistance
! EB_AV_0    : 0th moment avg collection eff. for Brownian diffusion
! EB_AV_3    : 3rd moment avg collection eff. for Brownian diffusion
! EIM_AV_0   : 0th moment avg collection eff. for impaction
! EIM_AV_3   : 3rd moment avg collection eff. for impaction
! EIN        : Collection eff. for interception
! AR         : Aerodynamic resistance
! MTOT       : Total aerosol mass conc [all cpts] (molecules/cm3)
! MCPTOT     : Total aersool mass conc [1 cpt] (molecules/cm3)
! NEWN       : Updated number concentration (/cm3)
! DZ         : Ht difference between box vertical interfaces (m)
! DZMID      : Ht difference between box lower interface & mid-level (m)
! SIGMA      : Geometric standard deviation of mode
! MASK      : Logical to define regions of domain to work on.
! MASK_SMOO :Logical to define regsions over "smooth"  surface categories
! MASK_VEGE :Logical to define regsions over vegetated surface categories
!
! Inputted by module UKCA_UM_LEGACY_MOD
! -------------------------------------
! GG        : Gravitational acceleration (ms^-2)
! VKMAN     : Von Karman's constant
! RGAS      : Dry air gas constant (Jkg^-1 K^-1)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! NCP       : Number of possible aerosol components
! MODE      : Defines which modes are set
! COMPONENT : Defines which cpts are allowed in each mode
! SIGMAG    : Geometric standard deviation of mode
! MM        : Molar masses of components (kg/mole)
! NUM_EPS   : Value of NEWN below which do not recalculate MD (per cc)
!                                              or carry out process
! CP_SU     : Index of component in which H2SO4 is stored
! CP_BC     : Index of component in which BC is stored
! CP_OC     : Index of component in which 1st OC cpt is stored
! CP_CL     : Index of component in which NaCl is stored
! CP_DU     : Index of component in which dust   cpt is stored
! CP_SO     : Index of component in which 2nd OC cpt is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! Various indices for budget terms in BUD_AER_MAS
!
! Inputted by module UKCA_DDEPAER_COEFF_MOD
!-------------------------------------------------
! CR,Y,ALPHA: aerosol deposition coefficients
!     [vary with land category & input via DATA statements]
! CR        : Characteristic radius of collectors (m)
! Y         : Parameter for calculating Brownian diffusion
! ALPHA     : Parameter for calculating EIM

! References
! ----------
! Slinn, Atmos. En., 1982, 16, 1785-1794
! Zhang et al, Atmos. En., 2001, 35, 549-560
!
!----------------------------------------------------------------------

USE parkind1,                ONLY: jpim, jprb
USE ukca_dcoff_par_av_k_mod, ONLY: ukca_dcoff_par_av_k

USE ukca_config_specification_mod, ONLY: glomap_variables

USE ukca_mode_setup,         ONLY: nmodes,                                     &
                                   cp_su, cp_bc, cp_oc, cp_cl,                 &
                                   cp_du, cp_so, moment_number, moment_mass,   &
                                   mode_nuc_sol, mode_ait_sol, mode_acc_sol,   &
                                   mode_cor_sol, mode_ait_insol,               &
                                   mode_acc_insol, mode_cor_insol,             &
                                   cp_no3, cp_nh4, cp_nn,                      &
                                   mode_sup_insol

USE ukca_setup_indices,      ONLY: nmasddepbcaccsol,                           &
                                   nmasddepbcaitins, nmasddepbcaitsol,         &
                                   nmasddepbccorsol, nmasddepduaccins,         &
                                   nmasddepduaccsol, nmasddepducorins,         &
                                   nmasddepducorsol, nmasddepocaccsol,         &
                                   nmasddepocaitins, nmasddepocaitsol,         &
                                   nmasddepoccorsol, nmasddepocnucsol,         &
                                   nmasddepsoaccsol, nmasddepsoaitsol,         &
                                   nmasddepsocorsol, nmasddepsonucsol,         &
                                   nmasddepssaccsol, nmasddepsscorsol,         &
                                   nmasddepsuaccsol, nmasddepsuaitsol,         &
                                   nmasddepsucorsol, nmasddepsunucsol,         &
                                   nmasddepntaitsol, nmasddepntaccsol,         &
                                   nmasddepntcorsol, nmasddepnhaitsol,         &
                                   nmasddepnhaccsol, nmasddepnhcorsol,         &
                                   nmasddepnnaccsol, nmasddepnncorsol,         &
                                   nmasddepdusupins
USE ukca_ddepaer_coeff_mod,  ONLY: alpha, cr, yr, ls_ice, ls_ocean, ls_soil,   &
                                   ls_water

USE ukca_vgrav_av_k_mod,     ONLY: ukca_vgrav_av_k
USE ukca_um_legacy_mod,      ONLY: vkman, gg => g, rgas => r
USE ukca_types_mod,          ONLY: log_small
USE yomhook,                 ONLY: lhook, dr_hook

IMPLICIT NONE

! .. Subroutine interface
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: nbudaer
INTEGER, INTENT(IN) :: ilscat(nbox)
REAL, INTENT(IN)    :: rhopar(nbox,nmodes)
REAL, INTENT(IN)    :: znot(nbox)
REAL, INTENT(IN)    :: dtc
REAL, INTENT(IN)    :: wetdp(nbox,nmodes)
REAL, INTENT(IN)    :: ustr(nbox)
REAL, INTENT(IN)    :: pmid(nbox)
REAL, INTENT(IN)    :: pupper(nbox)
REAL, INTENT(IN)    :: plower(nbox)
REAL, INTENT(IN)    :: t(nbox)
REAL, INTENT(IN)    :: surtp(nbox)
REAL, INTENT(IN)    :: rhoa(nbox)
REAL, INTENT(IN)    :: mfpa(nbox)
REAL, INTENT(IN)    :: dvisc(nbox)
REAL, INTENT(IN OUT) :: nd(nbox,nmodes)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables%ncp)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
REAL, INTENT(IN OUT) :: bud_aer_mas(nbox,0:nbudaer)

!    Local Variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: num_eps(:)
REAL,    POINTER :: sigmag(:)

INTEGER :: imode
INTEGER :: icp
LOGICAL (KIND=log_small) :: mask3(nbox)
LOGICAL (KIND=log_small) :: mask4(nbox)
LOGICAL (KIND=log_small) :: masksurf(nbox)
LOGICAL (KIND=log_small) :: mask_smoo(nbox)
LOGICAL (KIND=log_small) :: mask_vege(nbox)
REAL    :: delnddep(nbox)
REAL    :: delmddep(nbox)
REAL    :: ps_av_0(nbox)
REAL    :: ps_av_3(nbox)
REAL    :: kvisc(nbox)
REAL    :: vgrav_av_0(nbox)
REAL    :: vgrav_av_3(nbox)
REAL    :: dcoef_av_0(nbox)
REAL    :: dcoef_av_3(nbox)
REAL    :: eb_av_0(nbox)
REAL    :: eb_av_3(nbox)
REAL    :: eim_av_0(nbox)
REAL    :: eim_av_3(nbox)
REAL    :: ein(nbox)
REAL    :: sn_av_0(nbox)
REAL    :: sn_av_3(nbox)
REAL    :: ar(nbox)
REAL    :: sr_av_0(nbox)
REAL    :: sr_av_3(nbox)
REAL    :: vdep_av_0(nbox)
REAL    :: vdep_av_3(nbox)
REAL    :: mtot(nbox)
REAL    :: mcptot(nbox)
REAL    :: term1(nbox)
REAL    :: newn(nbox)
REAL    :: dzmid(nbox)
REAL    :: dz(nbox)

!--------------------------------------------------------------------

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DDEPAER'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
num_eps     => glomap_variables%num_eps
sigmag      => glomap_variables%sigmag

dzmid(:)=(rgas*t(:)/gg)*LOG(plower(:)/pmid(:))
dz(:)=(rgas*t(:)/gg)*LOG(plower(:)/pupper(:))

! .. Calculate aerodynamic resistance
ar(:)=LOG(dzmid(:)/znot(:))/(vkman*ustr(:))

!    Loop over modes
DO imode=1,nmodes
  IF (mode(imode)) THEN

    !       Calculate 0th moment avg. grav. settling velocities
    CALL ukca_vgrav_av_k(nbox,moment_number,wetdp(:,imode),sigmag(imode),      &
        dvisc(:),mfpa(:),rhopar(:,imode),vgrav_av_0(:))

    !       Calculate 3rd moment avg. grav. settling velocities
    CALL ukca_vgrav_av_k(nbox,moment_mass,wetdp(:,imode),sigmag(imode),        &
        dvisc(:),mfpa(:),rhopar(:,imode),vgrav_av_3(:))


    !       Calculate 0th moment avg particle diffusion coeffs
    CALL ukca_dcoff_par_av_k(nbox,moment_number,wetdp(:,imode),sigmag(imode),  &
             t(:),dvisc(:),mfpa(:),dcoef_av_0(:))

    !       Calculate 3rd moment avg particle diffusion coeffs
    CALL ukca_dcoff_par_av_k(nbox,moment_mass,wetdp(:,imode),sigmag(imode),    &
             t(:),dvisc(:),mfpa(:),dcoef_av_3(:))

    !      Calculate kinematic viscosity of air
    kvisc(:)=dvisc(:)/rhoa(:)

    !      Calculate 0th and 3rd moment avg. particle Schmidt number
    ps_av_0(:)=kvisc(:)/dcoef_av_0(:)
    ps_av_3(:)=kvisc(:)/dcoef_av_3(:)
    !      Calculate particle collection efficiencies
    !       -- For Brownian Diffusion
    eb_av_0(:)=ps_av_0(:)**(-yr(ilscat(:)))
    eb_av_3(:)=ps_av_3(:)**(-yr(ilscat(:)))

    ! Set smooth surfaces to be water, soil or ice
    ! All other surfaces are vegetated (have CR>0)
    mask_smoo=( (ilscat(:) == ls_ice) .OR. (ilscat(:) == ls_ocean) .OR.        &
                (ilscat(:) == ls_soil) .OR. (ilscat(:) == ls_water) )
    mask_vege= .NOT. mask_smoo(:)

    !       -- For Impaction
    WHERE (mask_smoo(:))
      !        Calculate stokes number for smooth surfaces
      sn_av_0(:)=vgrav_av_0(:)*ustr(:)*ustr(:)/dvisc(:)
      sn_av_3(:)=vgrav_av_3(:)*ustr(:)*ustr(:)/dvisc(:)
    END WHERE
    WHERE (mask_vege(:))
      !        Calculate stokes number for vegetated surfcaes
      sn_av_0(:)=vgrav_av_0(:)*ustr(:)/(gg*cr(ilscat(:)))
      sn_av_3(:)=vgrav_av_3(:)*ustr(:)/(gg*cr(ilscat(:)))
    END WHERE

    eim_av_0(:)=(sn_av_0(:)/(alpha(ilscat(:))+sn_av_0(:)))**2
    eim_av_3(:)=(sn_av_3(:)/(alpha(ilscat(:))+sn_av_3(:)))**2

    !       -- For Interception
    WHERE (mask_smoo(:))
      ein(:)=0.0
    END WHERE
    WHERE (mask_vege(:))
      ein(:)=0.5*(wetdp(:,imode)*wetdp(:,imode)                                &
                 /cr(ilscat(:))/cr(ilscat(:)))
    END WHERE

    !       Calculate surface resistance
    sr_av_0(:)=1.0/(3.0*ustr(:)*(eb_av_0(:)+eim_av_0(:)+ein(:)))
    sr_av_3(:)=1.0/(3.0*ustr(:)*(eb_av_3(:)+eim_av_3(:)+ein(:)))
    !       Calculate deposition velocity
    vdep_av_0(:)=vgrav_av_0(:)+1.0/(ar(:)+sr_av_0(:))
    vdep_av_3(:)=vgrav_av_3(:)+1.0/(ar(:)+sr_av_3(:))

    masksurf(:)=(surtp(:) < 2.0) ! boxes at surface.
    mask3(:)=(nd(:,imode) > num_eps(imode))
    mask4(:)=mask3(:) .AND. masksurf(:) ! also at surface

    WHERE (mask4(:)) ! only do at surface & where some particles

      delnddep(:)=nd(:,imode)*(1.0-EXP(-vdep_av_0(:)*dtc/dz(:)))

      !        Set updated particle concentration to NEWN
      newn(:)=nd(:,imode)-delnddep(:)
      !
      !        Update total mass per particle MDT
      mtot(:)=nd(:,imode)*mdt(:,imode)
      mdt(:,imode)=mtot*EXP(-vdep_av_3(:)*dtc/dz(:))/newn(:)

    END WHERE

    DO icp=1,ncp
      IF (component(imode,icp)) THEN

        WHERE (mask4(:)) ! only do at surface & where some particles

          mcptot(:)=nd(:,imode)*md(:,imode,icp)
          term1(:)=EXP(-vdep_av_3(:)*dtc/dz(:))
          md(:,imode,icp)=mcptot(:)*term1(:)/newn(:)
          delmddep(:)=mcptot(:)*(1.0-term1(:))

        END WHERE

        ! .. only store budgets at surface & where some particles
        IF (icp == cp_su) THEN
          IF ((imode == mode_nuc_sol) .AND. (nmasddepsunucsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepsunucsol)=                                  &
             bud_aer_mas(:,nmasddepsunucsol)+delmddep(:)
          END IF
          IF ((imode == mode_ait_sol) .AND. (nmasddepsuaitsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepsuaitsol)=                                  &
             bud_aer_mas(:,nmasddepsuaitsol)+delmddep(:)
          END IF
          IF ((imode == mode_acc_sol) .AND. (nmasddepsuaccsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepsuaccsol)=                                  &
             bud_aer_mas(:,nmasddepsuaccsol)+delmddep(:)
          END IF
          IF ((imode == mode_cor_sol) .AND. (nmasddepsucorsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepsucorsol)=                                  &
             bud_aer_mas(:,nmasddepsucorsol)+delmddep(:)
          END IF
        END IF
        IF (icp == cp_bc) THEN
          IF ((imode == mode_ait_sol) .AND. (nmasddepbcaitsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepbcaitsol)=                                  &
             bud_aer_mas(:,nmasddepbcaitsol)+delmddep(:)
          END IF
          IF ((imode == mode_acc_sol) .AND. (nmasddepbcaccsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepbcaccsol)=                                  &
             bud_aer_mas(:,nmasddepbcaccsol)+delmddep(:)
          END IF
          IF ((imode == mode_cor_sol) .AND. (nmasddepbccorsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepbccorsol)=                                  &
             bud_aer_mas(:,nmasddepbccorsol)+delmddep(:)
          END IF
          IF ((imode == mode_ait_insol) .AND. (nmasddepbcaitins > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepbcaitins)=                                  &
             bud_aer_mas(:,nmasddepbcaitins)+delmddep(:)
          END IF
        END IF
        IF (icp == cp_oc) THEN
          IF ((imode == mode_nuc_sol) .AND. (nmasddepocnucsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepocnucsol)=                                  &
             bud_aer_mas(:,nmasddepocnucsol)+delmddep(:)
          END IF
          IF ((imode == mode_ait_sol) .AND. (nmasddepocaitsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepocaitsol)=                                  &
             bud_aer_mas(:,nmasddepocaitsol)+delmddep(:)
          END IF
          IF ((imode == mode_acc_sol) .AND. (nmasddepocaccsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepocaccsol)=                                  &
             bud_aer_mas(:,nmasddepocaccsol)+delmddep(:)
          END IF
          IF ((imode == mode_cor_sol) .AND. (nmasddepoccorsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepoccorsol)=                                  &
             bud_aer_mas(:,nmasddepoccorsol)+delmddep(:)
          END IF
          IF ((imode == mode_ait_insol) .AND. (nmasddepocaitins > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepocaitins)=                                  &
             bud_aer_mas(:,nmasddepocaitins)+delmddep(:)
          END IF
        END IF
        IF (icp == cp_cl) THEN
          IF ((imode == mode_acc_sol) .AND. (nmasddepssaccsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepssaccsol)=                                  &
             bud_aer_mas(:,nmasddepssaccsol)+delmddep(:)
          END IF
          IF ((imode == mode_cor_sol) .AND. (nmasddepsscorsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepsscorsol)=                                  &
             bud_aer_mas(:,nmasddepsscorsol)+delmddep(:)
          END IF
        END IF
        IF (icp == cp_so) THEN
          IF ((imode == mode_nuc_sol) .AND. (nmasddepsonucsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepsonucsol)=                                  &
             bud_aer_mas(:,nmasddepsonucsol)+delmddep(:)
          END IF
          IF ((imode == mode_ait_sol) .AND. (nmasddepsoaitsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepsoaitsol)=                                  &
             bud_aer_mas(:,nmasddepsoaitsol)+delmddep(:)
          END IF
          IF ((imode == mode_acc_sol) .AND. (nmasddepsoaccsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepsoaccsol)=                                  &
             bud_aer_mas(:,nmasddepsoaccsol)+delmddep(:)
          END IF
          IF ((imode == mode_cor_sol) .AND. (nmasddepsocorsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepsocorsol)=                                  &
             bud_aer_mas(:,nmasddepsocorsol)+delmddep(:)
          END IF
        END IF
        IF (icp == cp_du) THEN
          IF ((imode == mode_acc_sol) .AND. (nmasddepduaccsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepduaccsol)=                                  &
             bud_aer_mas(:,nmasddepduaccsol)+delmddep(:)
          END IF
          IF ((imode == mode_cor_sol) .AND. (nmasddepducorsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepducorsol)=                                  &
             bud_aer_mas(:,nmasddepducorsol)+delmddep(:)
          END IF
          IF ((imode == mode_acc_insol) .AND. (nmasddepduaccins > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepduaccins)=                                  &
             bud_aer_mas(:,nmasddepduaccins)+delmddep(:)
          END IF
          IF ((imode == mode_cor_insol) .AND. (nmasddepducorins > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepducorins)=                                  &
             bud_aer_mas(:,nmasddepducorins)+delmddep(:)
          END IF
          IF ((imode == mode_sup_insol) .AND. (nmasddepdusupins > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepdusupins)=                                  &
             bud_aer_mas(:,nmasddepdusupins)+delmddep(:)
          END IF
        END IF
        ! Nitrate scheme (MS10)
        IF (icp == cp_no3) THEN
          !IF ((imode == mode_nuc_sol) .AND. (nmasddepntnucsol > 0)) THEN
          !  WHERE (mask4(:))                                                  &
          !   bud_aer_mas(:,nmasddepntnucsol)=                                 &
          !   bud_aer_mas(:,nmasddepntnucsol)+delmddep(:)
          !END IF
          IF ((imode == mode_ait_sol) .AND. (nmasddepntaitsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepntaitsol)=                                  &
             bud_aer_mas(:,nmasddepntaitsol)+delmddep(:)
          END IF
          IF ((imode == mode_acc_sol) .AND. (nmasddepntaccsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepntaccsol)=                                  &
             bud_aer_mas(:,nmasddepntaccsol)+delmddep(:)
          END IF
          IF ((imode == mode_cor_sol) .AND. (nmasddepntcorsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepntcorsol)=                                  &
             bud_aer_mas(:,nmasddepntcorsol)+delmddep(:)
          END IF
        END IF
        IF (icp == cp_nh4) THEN
          !IF ((imode == mode_nuc_sol) .AND. (nmasddepnhnucsol > 0)) THEN
          !  WHERE (mask4(:))                                                  &
          !   bud_aer_mas(:,nmasddepnhnucsol)=                                 &
          !   bud_aer_mas(:,nmasddepnhnucsol)+delmddep(:)
          !END IF
          IF ((imode == mode_ait_sol) .AND. (nmasddepnhaitsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepnhaitsol)=                                  &
             bud_aer_mas(:,nmasddepnhaitsol)+delmddep(:)
          END IF
          IF ((imode == mode_acc_sol) .AND. (nmasddepnhaccsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepnhaccsol)=                                  &
             bud_aer_mas(:,nmasddepnhaccsol)+delmddep(:)
          END IF
          IF ((imode == mode_cor_sol) .AND. (nmasddepnhcorsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepnhcorsol)=                                  &
             bud_aer_mas(:,nmasddepnhcorsol)+delmddep(:)
          END IF
        END IF
        IF (icp == cp_nn) THEN
          IF ((imode == mode_acc_sol) .AND. (nmasddepnnaccsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepnnaccsol)=                                  &
             bud_aer_mas(:,nmasddepnnaccsol)+delmddep(:)
          END IF
          IF ((imode == mode_cor_sol) .AND. (nmasddepnncorsol > 0)) THEN
            WHERE (mask4(:))                                                   &
             bud_aer_mas(:,nmasddepnncorsol)=                                  &
             bud_aer_mas(:,nmasddepnncorsol)+delmddep(:)
          END IF
        END IF
      END IF ! if component present in mode
    END DO ! loop over components

    !       Update number concentration to NEW
    !        (only do at surface & where some particles)
    WHERE (mask4(:)) nd(:,imode)=newn(:)

  END IF ! if mode present
END DO ! loop over modes

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_ddepaer

END MODULE ukca_ddepaer_mod
