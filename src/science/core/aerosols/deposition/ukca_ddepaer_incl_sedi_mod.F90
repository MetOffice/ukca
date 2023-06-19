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
!    Calculates aerosol dry deposition and sedimentation.
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
MODULE ukca_ddepaer_incl_sedi_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
                  ModuleName='UKCA_DDEPAER_INCL_SEDI_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE ukca_ddepaer_incl_sedi(nbox,nbudaer,nd,md,mdt,rhopar,znot,          &
 dtc,wetdp,ustr,pmid,pupper,plower,t,surtp,                                    &
 rhoa,mfpa,dvisc,bud_aer_mas,jlabove,ilscat,sedi_on,sm)

!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculates aerosol dry deposition and sedimentation.
! Based on the parameterisation of Zhang et al (2001) which
! uses the method in the model of Slinn (1982).
!
! Sedimentation is done using a simple explicit discretization
! which should be adequate for this process-split method.
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
! Note --- in this routine, sedimentation is included at all levels.
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
! JLABOVE   : Index of box directly above this grid box
! ILSCAT    : Land-use category (1-9 based on UM landsurf types)
! SEDI_ON   : Switch for whether aerosol sedimentation is on/off
! SM        : Grid box mass of air (kg)
!
! Calls subroutine GETROUGH to read roughness length
!
! Outputs
! -------
! Updated particle number density ND (/cm3)
! Updated particle avg mass MD (molecules/particle)
! Updated Avg tot mass of aerosol ptcl in size mode MDT (particle^-1)
! BUD_AER_MAS: Aerosol mass budgets (mlcls/cc/tstep)
!
! Local Variables
! ---------------
! PS_AV_0    : 0th moment avg particle Schmidt Number
! PS_AV_3    : 3rd moment avg particle Schmidt Number
! KVISC      : Kinematic viscosity of air (m2 s-1)
! VGRAV_AV_0 : 0th moment avg grav. settling vel. (m s^-1)
! VGRAV_AV_3 : 3rd moment avg grav. settling vel. (m s^-1)
! VDEP_AV_0  : 0th moment avg deposition velocity (m s^-1)
! VDEP_AV_3  : 3rd moment avg deposition velocity (m s^-1)
! DCOEF_AV_0 : 0th moment avg particle diffusion coefficient(m2 s-1)
! DCOEF_AV_3 : 3rd moment avg particle diffusion coefficient(m2 s-1)
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
! NEWN       : Updated number concentration (/cm3)
! DZ         : Ht difference between box vertical interfaces (m)
! DZMID      : Ht difference between box lower interface & mid-level (m)
! MASK1     : Logical to define regions of domain for roughness categories
! MASK2     : Logical to define regions of domain for JLABOVE
! MASK3     : Logical to define regions with some aerosol in box or 1above
! MASK4     : Logical to define surface regions with aerosol in box/1above
! MASK_SMOO :Logical to define regsions over "smooth"  surface categories
! MASK_VEGE :Logical to define regsions over vegetated surface categories
! MASK_ABOVE_LIM: Logical to define boxes with VGRAV > allowed max value
! cfl_fraction: Safety margin by which to stay within CFL limits
! dtsedi    : Desired sedimentation timestep for each mode
! dtmode    : Actual sedimentation timestep for IMODE (integer
!             fraction of DTC to ensure nsedi is an integer)
! nsedi     : Number of sedimentation substeps for IMODE
! isedi     : Sedimentation substep interation index
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
! NUM_EPS   : Value of ND_0 below which do not recalculate MD (per cc)
!                                              or carry out process
! CP_SU     : Index of component in which H2SO4 is stored
! CP_BC     : Index of component in which BC is stored
! CP_OC     : Index of component in which 1st OC cpt is stored
! CP_CL     : Index of component in which NaCl is stored
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
!
! References
! ----------
! Slinn, Atmos. En., 1982, 16, 1785-1794
! Zhang et al, Atmos. En., 2001, 35, 549-560
!
!----------------------------------------------------------------------

USE parkind1,                ONLY:                                             &
    jpim,                                                                      &
    jprb

USE ukca_dcoff_par_av_k_mod, ONLY:                                             &
    ukca_dcoff_par_av_k

USE ukca_config_specification_mod, ONLY:                                       &
    glomap_variables

USE ukca_mode_setup,         ONLY:                                             &
    nmodes,                                                                    &
    cp_su, cp_bc, cp_oc, cp_cl, cp_du, cp_so, cp_no3, cp_nh4, cp_nn,           &
    moment_number, moment_mass,                                                &
    mode_nuc_sol, mode_ait_sol,                                                &
    mode_acc_sol, mode_cor_sol,                                                &
    mode_ait_insol, mode_acc_insol,                                            &
    mode_cor_insol, mode_sup_insol

USE ukca_setup_indices,      ONLY:                                             &
    nmasddepsunucsol,                                                          &
    nmasddepsuaitsol, nmasddepsuaccsol,                                        &
    nmasddepsucorsol, nmasddepbcaitsol,                                        &
    nmasddepbcaccsol, nmasddepbccorsol,                                        &
    nmasddepbcaitins, nmasddepocnucsol,                                        &
    nmasddepocaitsol, nmasddepocaccsol,                                        &
    nmasddepoccorsol, nmasddepocaitins,                                        &
    nmasddepssaccsol, nmasddepsscorsol,                                        &
    nmasddepsonucsol, nmasddepsoaitsol,                                        &
    nmasddepsoaccsol, nmasddepsocorsol,                                        &
    nmasddepduaccsol, nmasddepducorsol,                                        &
    nmasddepduaccins, nmasddepducorins,                                        &
    nmasddepntaitsol, nmasddepntaccsol,                                        &
    nmasddepntcorsol, nmasddepnhaitsol,                                        &
    nmasddepnhaccsol, nmasddepnhcorsol,                                        &
    nmasddepnnaccsol, nmasddepnncorsol,                                        &
    nmasddepdusupins

USE ukca_ddepaer_coeff_mod,  ONLY: alpha, cr, yr, ls_ice, ls_ocean, ls_soil,   &
                                   ls_water

USE ukca_vgrav_av_k_mod,     ONLY:                                             &
    ukca_vgrav_av_k

USE ukca_um_legacy_mod,      ONLY:                                             &
    log_v,                                                                     &
    vkman,                                                                     &
    gg => g,                                                                   &
    rgas => r

USE ukca_types_mod,          ONLY:                                             &
    log_small


USE yomhook,                 ONLY:                                             &
    lhook,                                                                     &
    dr_hook


IMPLICIT NONE

! .. Subroutine interface
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: nbudaer
INTEGER, INTENT(IN) :: sedi_on
INTEGER, INTENT(IN) :: jlabove(nbox)
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
REAL, INTENT(IN)    :: sm(nbox)
REAL, INTENT(IN OUT) :: bud_aer_mas(nbox,0:nbudaer)
REAL, INTENT(IN OUT) :: nd(nbox,nmodes)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables%ncp)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)

! .. Local Variables

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
INTEGER :: isedi
INTEGER :: nsedi
LOGICAL (KIND=log_small) :: mask2(nbox)
LOGICAL (KIND=log_small) :: mask3(nbox)
LOGICAL (KIND=log_small) :: mask4(nbox)
LOGICAL (KIND=log_small) :: masksurf(nbox)
LOGICAL (KIND=log_small) :: mask_smoo(nbox)
LOGICAL (KIND=log_small) :: mask_vege(nbox)
LOGICAL (KIND=log_small) :: mask_above_lim(nbox)
REAL    :: dtmode
REAL    :: ps_av_0(nbox)
REAL    :: ps_av_3(nbox)
REAL    :: kvisc(nbox)
REAL    :: vgrav_av_0(nbox)
REAL    :: vgrav_av_3(nbox)
REAL    :: vgrav_av_0_up(nbox)
REAL    :: vgrav_av_3_up(nbox)
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
REAL    :: dzmid(nbox)
REAL    :: dz(nbox)
REAL    :: dz_up(nbox)
REAL    :: t_up(nbox)
REAL    :: plo_up(nbox)
REAL    :: pup_up(nbox)
REAL    :: sm_up(nbox)
REAL    :: rhoa_up(nbox)
REAL    :: nd0(nbox,nmodes)
REAL    :: md0(nbox,nmodes,glomap_variables%ncp)
REAL    :: nd0_up(nbox,nmodes)
REAL    :: md0_up(nbox,nmodes,glomap_variables%ncp)
REAL    :: ndnew(nbox)
REAL    :: termin_1(nbox)
REAL    :: termin_2(nbox)
REAL    :: termin_n(nbox)
REAL    :: termout_n(nbox)
REAL    :: termin_m(nbox,glomap_variables%ncp)
REAL    :: termout_m(nbox,glomap_variables%ncp)
REAL    :: delnsedi(nbox)
REAL    :: delmsedi(nbox,glomap_variables%ncp)
REAL    :: delmddep(nbox)
REAL    :: vgrav_lim(nbox)
REAL, PARAMETER :: cfl_fraction = 0.9
REAL, PARAMETER :: dtsedi(nmodes) =                                            &
                [ 3600.0, 3600.0, 1800.0, 900.0, 3600.0, 1800.0, 900.0, 300.0 ]

REAL :: p_ratio(nbox)
REAL :: p_ratio_out(nbox)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName='UKCA_DDEPAER_INCL_SEDI'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
num_eps     => glomap_variables%num_eps
sigmag      => glomap_variables%sigmag

mask2(:)=(jlabove(:) > 0) ! where using JLABOVE

WHERE (mask2(:))
  t_up(:)=t(jlabove(:))
  plo_up(:)=plower(jlabove(:))
  pup_up(:)=pupper(jlabove(:))
  sm_up(:)=sm(jlabove(:))
  rhoa_up(:)=rhoa(jlabove(:))
ELSE WHERE
  t_up(:)=t(:)
  plo_up(:)=plower(:)
  pup_up(:)=pupper(:)
  sm_up(:)=sm(:)
  rhoa_up(:)=rhoa(:)
END WHERE

p_ratio(:)=plower(:)/pupper(:)
CALL log_v(nbox,p_ratio,p_ratio_out)
dz(:)=rgas*t(:)* p_ratio_out(:)/gg

p_ratio(:)=plower(:)/pmid(:)
CALL log_v(nbox,p_ratio,p_ratio_out)
dzmid(:)=rgas*t(:)* p_ratio_out(:)/gg

p_ratio(:)=plo_up(:)/pup_up(:)
CALL log_v(nbox,p_ratio,p_ratio_out)
dz_up(:)=rgas*t_up(:)* p_ratio_out(:)/gg

kvisc(:)=dvisc(:)/rhoa(:) ! Calculate kinematic viscosity of air

masksurf(:)=(surtp(:) < 2.0) ! create mask for boxes at surface.

! .. Calculate aerodynamic resistance (only in gridboxes at surface)
WHERE (masksurf(:)) ar(:)=LOG(dzmid(:)/znot(:))/(vkman*ustr(:))

DO imode=1,nmodes
  IF (mode(imode)) THEN

     ! .. Compute num iters from timestep dtmode, then recompute dtmode to
     ! .. force it to be an exact divisor of DTC.  Ensure dtmode is not larger
     ! .. than DTC.  Round nsedi up (instead of down or nearest) to ensure
     ! .. timestep is no longer than dtsedi(IMODE)
    dtmode = MIN(dtsedi(imode), dtc)
    nsedi = CEILING(dtc/dtmode)
    dtmode = dtc / REAL(nsedi)

    DO isedi = 1, nsedi

      ! .. First copy original values of ND,MD to ND0,MD0
      nd0(:,imode)=nd(:,imode)
      md0(:,imode,:)=md(:,imode,:)

      WHERE (mask2(:))
        nd0_up(:,imode)=nd0(jlabove(:),imode)
      ELSE WHERE
        nd0_up(:,imode)=0.0
      END WHERE
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WHERE (mask2(:))
            md0_up(:,imode,icp)=md0(jlabove(:),imode,icp)
          ELSE WHERE
            md0_up(:,imode,icp)=md0(:,imode,icp)
          END WHERE
        END IF
      END DO

      ! .. Calculate 0th moment avg. grav. settling velocities
      CALL ukca_vgrav_av_k(nbox,moment_number,wetdp(:,imode),sigmag(imode),    &
             dvisc(:),mfpa(:),rhopar(:,imode),vgrav_av_0(:))
      ! .. Calculate 3rd moment avg. grav. settling velocities
      CALL ukca_vgrav_av_k(nbox,moment_mass,wetdp(:,imode),sigmag(imode),      &
             dvisc(:),mfpa(:),rhopar(:,imode),vgrav_av_3(:))
      !
      ! .. Store values from box above so that can get flux into box
      WHERE (mask2(:))
        vgrav_av_0_up(:)=vgrav_av_0(jlabove(:))
        vgrav_av_3_up(:)=vgrav_av_3(jlabove(:))
      ELSE WHERE
        vgrav_av_0_up(:)=0.0
        vgrav_av_3_up(:)=0.0
      END WHERE

      ! .. If sedimentation switched off, set all calculated VGRAV to zero
      IF (sedi_on == 0) THEN
        vgrav_av_0(:)=0.0
        vgrav_av_3(:)=0.0
        vgrav_av_0_up(:)=0.0
        vgrav_av_3_up(:)=0.0
      END IF

      !       Calculate 0th moment avg particle diffusion coeffs
      CALL ukca_dcoff_par_av_k(nbox,moment_number,wetdp(:,imode),              &
                sigmag(imode),t(:),dvisc(:),mfpa(:),dcoef_av_0(:))

      !       Calculate 3rd moment avg particle diffusion coeffs
      CALL ukca_dcoff_par_av_k(nbox,moment_mass,wetdp(:,imode),                &
                sigmag(imode),t(:),dvisc(:),mfpa(:),dcoef_av_3(:))

      ! .. only calculate Schmidt number and collection eff at surface
      WHERE (masksurf(:))

        !       Calculate 0th and 3rd moment avg. particle Schmidt number
        ps_av_0(:)=kvisc(:)/dcoef_av_0(:)
        ps_av_3(:)=kvisc(:)/dcoef_av_3(:)
        !       Calculate particle collection efficiencies
        !       -- For Brownian Diffusion
        eb_av_0(:)=ps_av_0(:)**(-yr(ilscat(:)))
        eb_av_3(:)=ps_av_3(:)**(-yr(ilscat(:)))
        !
      END WHERE

      ! Set smooth surfaces to be water, soil or ice
      ! All other surfaces are vegetated (have CR>0)
      mask_smoo=( (ilscat(:) == ls_ice) .OR. (ilscat(:) == ls_ocean) .OR.      &
                  (ilscat(:) == ls_soil) .OR. (ilscat(:) == ls_water) )
      mask_vege= .NOT. mask_smoo(:)

      ! .. only calculate Stokes number at surface
      WHERE (mask_smoo(:) .AND. masksurf(:))
        !        Calculate stokes number for smooth surfaces
        sn_av_0(:)=vgrav_av_0(:)*ustr(:)*ustr(:)/dvisc(:)
        sn_av_3(:)=vgrav_av_3(:)*ustr(:)*ustr(:)/dvisc(:)
      END WHERE
      WHERE (mask_vege(:) .AND. masksurf(:))
        !        Calculate stokes number for vegetated surfcaes
        sn_av_0(:)=vgrav_av_0(:)*ustr(:)/(gg*cr(ilscat(:)))
        sn_av_3(:)=vgrav_av_3(:)*ustr(:)/(gg*cr(ilscat(:)))
      END WHERE

      ! .. only calculate impaction collection efficiency at surface
      WHERE (masksurf(:))

        !       -- For Impaction
        eim_av_0(:)=(sn_av_0(:)/(alpha(ilscat(:))+sn_av_0(:)))**2
        eim_av_3(:)=(sn_av_3(:)/(alpha(ilscat(:))+sn_av_3(:)))**2

      END WHERE

      ! .. only calculate interception collection eff (smooth) at surface
      WHERE (mask_smoo(:) .AND. masksurf(:))

        !       -- For Interception (smooth surfaces)
        ein(:)=0.0

      END WHERE

      ! .. only calculate interception collection eff (vegetd) at surface
      WHERE (mask_vege(:) .AND. masksurf(:))

        !       -- For Interception (vegetd surfaces)
        ein(:)=0.5*(wetdp(:,imode)*wetdp(:,imode)                              &
                             /cr(ilscat(:))/cr(ilscat(:)))

      END WHERE

      ! .. only calculate surface resistance and deposition vel. at surface
      ! .. this section also increases VGRAV due to dry dep (at surface)
      WHERE (masksurf(:))

        !  Calculate surface resistance
        sr_av_0(:)=1.0/(3.0*ustr(:)*(eb_av_0(:)+eim_av_0(:)+ein(:)))
        sr_av_3(:)=1.0/(3.0*ustr(:)*(eb_av_3(:)+eim_av_3(:)+ein(:)))

        !  Calculate deposition velocity
        vdep_av_0(:)=vgrav_av_0(:)+1.0/(ar(:)+sr_av_0(:))
        vdep_av_3(:)=vgrav_av_3(:)+1.0/(ar(:)+sr_av_3(:))

        !  Set gravitational velocity to deposition velocity if in lowest box
        vgrav_av_0(:)=vdep_av_0(:)
        vgrav_av_3(:)=vdep_av_3(:)

        !  VGRAV_AV_UP never at surface so no need to set to dep. vel.

      END WHERE

      ! .. limit 0&3 grav. settling vel so only falls
      ! .. cfl_fraction*box max [numerical]
      vgrav_lim(:)=cfl_fraction*dz(:)/dtmode

      mask_above_lim(:)=(vgrav_av_0(:) > vgrav_lim(:))
      WHERE (mask_above_lim(:)) vgrav_av_0(:)=vgrav_lim(:)

      mask_above_lim(:)=(vgrav_av_3(:) > vgrav_lim(:))
      WHERE (mask_above_lim(:)) vgrav_av_3(:)=vgrav_lim(:)

      ! .. limit 0&3 grav. settling vel so only falls
      ! .. cfl_fraction*box max [box above]
      vgrav_lim(:)=cfl_fraction*dz_up(:)/dtmode

      mask_above_lim(:)=(vgrav_av_0_up(:) > vgrav_lim(:))
      WHERE (mask_above_lim(:)) vgrav_av_0_up(:)=vgrav_lim(:)

      mask_above_lim(:)=(vgrav_av_3_up(:) > vgrav_lim(:))
      WHERE (mask_above_lim(:)) vgrav_av_3_up(:)=vgrav_lim(:)

      ! .. Calculate sedimenting in term (number) using 0th moment VGRAV_AV_0
      WHERE (mask2(:)) ! where not top model level
        termin_1(:)=nd0_up(:,imode)/dz_up(:)
        termin_2(:)=sm_up(:)*rhoa(:)/sm(:)/rhoa_up(:)
        termin_n(:)=termin_1(:)*termin_2(:)*dtmode*vgrav_av_0_up(:)
      ELSE WHERE
        !    If top model level TERMIN_N(:)=0.0
        termin_n(:)=0.0
      END WHERE

      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WHERE (mask2(:))
            termin_1(:)=nd0_up(:,imode)*md0_up(:,imode,icp)/dz_up(:)
            termin_m(:,icp)=termin_1(:)*termin_2(:)*dtmode*vgrav_av_3_up(:)
          ELSE WHERE
            !    If top model level TERMIN_M(:,ICP)=0.0
            termin_m(:,icp)=0.0
          END WHERE
        END IF
      END DO

      ! .. Calculate sedimenting out term (number) using 0th moment VGRAV_AV_0
      termout_n(:)=(nd0(:,imode)/dz(:))*dtmode*vgrav_av_0(:)

      ! .. Calculate sedimenting in term (mass  ) using 3rd moment VGRAV_AV_3
      ! .. Calculate sedimenting out term (mass  ) using 3rd moment VGRAV_AV_3
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          termout_m(:,icp)=(nd0(:,imode)*md0(:,imode,icp)/dz(:))               &
                        *dtmode*vgrav_av_3(:)
        END IF
      END DO

      ! .. below calculates net change in number and mass concentration
      delnsedi(:)=-(termin_n(:)-termout_n(:))
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          delmsedi(:,icp)=-(termin_m(:,icp)-termout_m(:,icp))
        END IF
      END DO

      ! .. below masks boxes with some updates to apply
      mask3(:)=( (nd0   (:,imode) > num_eps(imode)) .OR.                       &
                 (nd0_up(:,imode) > num_eps(imode)) )
      mask4(:)=mask3(:) .AND. masksurf(:) ! some ptcls & also at surface

      ! .. below sets NDNEW to new value and re-sets MDT in boxes to update
      WHERE (mask3(:)) ! only do where some particles in/out
        ndnew(:)=nd0(:,imode)-delnsedi(:)
        mdt(:,imode)=0.0
      END WHERE

      ! .. below updates component masses and MDT for each mode
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WHERE (mask3(:)) ! only do where some particles in/out
            md(:,imode,icp)=                                                   &
              (md0(:,imode,icp)*nd0(:,imode)-delmsedi(:,icp))/ndnew(:)
            mdt(:,imode)=mdt(:,imode)+md(:,imode,icp)
          END WHERE
        END IF ! IF COMPONENT(ICP)
      END DO ! loop over cpts

      ! .. below updates number concentration to NDNEW
      ! ..  (only do where some particles in/out)
      WHERE (mask3(:)) nd(:,imode)=ndnew(:)

      ! .. below stores ddep/sedi fluxes to BUD_AER_MAS
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WHERE (mask4(:)) delmddep(:)=termout_m(:,icp)
          IF (icp == cp_su) THEN
            IF ((imode == mode_nuc_sol) .AND. (nmasddepsunucsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsunucsol)=                 &
                           bud_aer_mas(:,nmasddepsunucsol)+delmddep(:)
            IF ((imode == mode_ait_sol) .AND. (nmasddepsuaitsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsuaitsol)=                 &
                           bud_aer_mas(:,nmasddepsuaitsol)+delmddep(:)
            IF ((imode == mode_acc_sol) .AND. (nmasddepsuaccsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsuaccsol)=                 &
                           bud_aer_mas(:,nmasddepsuaccsol)+delmddep(:)
            IF ((imode == mode_cor_sol) .AND. (nmasddepsucorsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsucorsol)=                 &
                           bud_aer_mas(:,nmasddepsucorsol)+delmddep(:)
          END IF
          IF (icp == cp_bc) THEN
            IF ((imode == mode_ait_sol) .AND. (nmasddepbcaitsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepbcaitsol)=                 &
                        bud_aer_mas(:,nmasddepbcaitsol)+delmddep(:)
            IF ((imode == mode_acc_sol) .AND. (nmasddepbcaccsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepbcaccsol)=                 &
                        bud_aer_mas(:,nmasddepbcaccsol)+delmddep(:)
            IF ((imode == mode_cor_sol) .AND. (nmasddepbccorsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepbccorsol)=                 &
                        bud_aer_mas(:,nmasddepbccorsol)+delmddep(:)
            IF ((imode == mode_ait_insol) .AND. (nmasddepbcaitins > 0))        &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepbcaitins)=                 &
                        bud_aer_mas(:,nmasddepbcaitins)+delmddep(:)
          END IF
          IF (icp == cp_oc) THEN
            IF ((imode == mode_nuc_sol) .AND. (nmasddepocnucsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepocnucsol)=                 &
                        bud_aer_mas(:,nmasddepocnucsol)+delmddep(:)
            IF ((imode == mode_ait_sol) .AND. (nmasddepocaitsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepocaitsol)=                 &
                        bud_aer_mas(:,nmasddepocaitsol)+delmddep(:)
            IF ((imode == mode_acc_sol) .AND. (nmasddepocaccsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepocaccsol)=                 &
                        bud_aer_mas(:,nmasddepocaccsol)+delmddep(:)
            IF ((imode == mode_cor_sol) .AND. (nmasddepoccorsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepoccorsol)=                 &
                        bud_aer_mas(:,nmasddepoccorsol)+delmddep(:)
            IF ((imode == mode_ait_insol) .AND. (nmasddepocaitins > 0))        &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepocaitins)=                 &
                        bud_aer_mas(:,nmasddepocaitins)+delmddep(:)
          END IF
          IF (icp == cp_cl) THEN
            IF ((imode == mode_acc_sol) .AND. (nmasddepssaccsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepssaccsol)=                 &
                        bud_aer_mas(:,nmasddepssaccsol)+delmddep(:)
            IF ((imode == mode_cor_sol) .AND. (nmasddepsscorsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsscorsol)=                 &
                        bud_aer_mas(:,nmasddepsscorsol)+delmddep(:)
          END IF
          IF (icp == cp_so) THEN
            IF ((imode == mode_nuc_sol) .AND. (nmasddepsonucsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsonucsol)=                 &
                        bud_aer_mas(:,nmasddepsonucsol)+delmddep(:)
            IF ((imode == mode_ait_sol) .AND. (nmasddepsoaitsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsoaitsol)=                 &
                        bud_aer_mas(:,nmasddepsoaitsol)+delmddep(:)
            IF ((imode == mode_acc_sol) .AND. (nmasddepsoaccsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsoaccsol)=                 &
                        bud_aer_mas(:,nmasddepsoaccsol)+delmddep(:)
            IF ((imode == mode_cor_sol) .AND. (nmasddepsocorsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsocorsol)=                 &
                        bud_aer_mas(:,nmasddepsocorsol)+delmddep(:)
          END IF
          IF (icp == cp_du) THEN
            IF ((imode == mode_acc_sol) .AND. (nmasddepduaccsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepduaccsol)=                 &
                        bud_aer_mas(:,nmasddepduaccsol)+delmddep(:)
            IF ((imode == mode_cor_sol) .AND. (nmasddepducorsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepducorsol)=                 &
                        bud_aer_mas(:,nmasddepducorsol)+delmddep(:)
            IF ((imode == mode_acc_insol) .AND. (nmasddepduaccins > 0))        &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepduaccins)=                 &
                        bud_aer_mas(:,nmasddepduaccins)+delmddep(:)
            IF ((imode == mode_cor_insol) .AND. (nmasddepducorins > 0))        &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepducorins)=                 &
                        bud_aer_mas(:,nmasddepducorins)+delmddep(:)
            IF ((imode == mode_sup_insol) .AND. (nmasddepdusupins > 0))        &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepdusupins)=                 &
                        bud_aer_mas(:,nmasddepdusupins)+delmddep(:)
          END IF
          IF (icp == cp_no3) THEN
            !IF ((imode == mode_nuc_sol) .AND. (nmasddepntnucsol > 0))         &
            ! WHERE (mask4(:)) bud_aer_mas(:,nmasddepntnucsol)=                &
            !               bud_aer_mas(:,nmasddepntnucsol)+delmddep(:)
            IF ((imode == mode_ait_sol) .AND. (nmasddepntaitsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepntaitsol)=                 &
                           bud_aer_mas(:,nmasddepntaitsol)+delmddep(:)
            IF ((imode == mode_acc_sol) .AND. (nmasddepntaccsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepntaccsol)=                 &
                           bud_aer_mas(:,nmasddepntaccsol)+delmddep(:)
            IF ((imode == mode_cor_sol) .AND. (nmasddepntcorsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepntcorsol)=                 &
                           bud_aer_mas(:,nmasddepntcorsol)+delmddep(:)
          END IF
          IF (icp == cp_nh4) THEN
            !IF ((imode == mode_nuc_sol) .AND. (nmasddepnhnucsol > 0))         &
            ! WHERE (mask4(:)) bud_aer_mas(:,nmasddepnhnucsol)=                &
            !               bud_aer_mas(:,nmasddepnhnucsol)+delmddep(:)
            IF ((imode == mode_ait_sol) .AND. (nmasddepnhaitsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepnhaitsol)=                 &
                           bud_aer_mas(:,nmasddepnhaitsol)+delmddep(:)
            IF ((imode == mode_acc_sol) .AND. (nmasddepnhaccsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepnhaccsol)=                 &
                           bud_aer_mas(:,nmasddepnhaccsol)+delmddep(:)
            IF ((imode == mode_cor_sol) .AND. (nmasddepnhcorsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepnhcorsol)=                 &
                           bud_aer_mas(:,nmasddepnhcorsol)+delmddep(:)
          END IF
          IF (icp == cp_nn) THEN
            IF ((imode == mode_acc_sol) .AND. (nmasddepnnaccsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepnnaccsol)=                 &
                           bud_aer_mas(:,nmasddepnnaccsol)+delmddep(:)
            IF ((imode == mode_cor_sol) .AND. (nmasddepnncorsol > 0))          &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepnncorsol)=                 &
                           bud_aer_mas(:,nmasddepnncorsol)+delmddep(:)
          END IF
        END IF ! if component present in mode
      END DO ! loop over components
    END DO ! loop over sedimentation timesteps
  END IF ! IF MODE(IMODE)
END DO ! loop over modes

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_ddepaer_incl_sedi

END MODULE ukca_ddepaer_incl_sedi_mod
