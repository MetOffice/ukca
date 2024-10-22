! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Calculate the reaction rate coefficients for use in the
!  Backward Euler solver with RAQ chemistry (based on STOCHEM ).
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method: Calculate reaction rates for RAQ chemistry by using
!   data taken from ukca_chemco_raq_init.
!   This routine was adapted from the modset stochem_chemistry.mf77
!   and some reaction rate coefficients were corrected following
!   Atkinson (2004, 2006), Atm. Chem. Phys.
!   If RAQ-AERO is used then add the relevant reaction rates
!   (e.g. sulphur chemistry and formation of secondary organic
!    aerosol, SOA).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  Fortran 95
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_chemco_raq_mod

IMPLICIT NONE

INTEGER, SAVE :: in2o5 = 0    ! index for heterog hydrolysis of N2O5
INTEGER, SAVE :: in2o5_ss = 0 ! index for heterog hydrolysis of N2O5 on seasalt

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CHEMCO_RAQ_MOD'

CONTAINS

SUBROUTINE ukca_chemco_raq (nr, n_pnts, tc, m, h2o, o2,                        &
                            clw,           fcloud    , fdiss,                  &
                            rho_air,       rh_frac,                            &
                            so4_aitken,    so4_accum ,                         &
                            soot_fresh,    soot_aged ,                         &
                            ocff_fresh,    ocff_aged , biogenic,               &
                            sea_salt_film, sea_salt_jet,                       &
                            stratflag, rc, hrc, H_plus_1d_arr)
!
USE ukca_chemco_raq_init_mod, ONLY: ukca_chemco_raq_init, rk, pdep
USE asad_mod, ONLY: peps, jpeq, jpspec, jphk
USE ukca_constants, ONLY: m_air, m_n2o5

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ukca_config_constants_mod, ONLY: rho_water, avogadro
USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_mode_setup,        ONLY: cp_su, cp_bc, cp_oc, cp_cl, cp_so
USE ukca_um_legacy_mod,     ONLY: calc_surf_area, l_um_calc_surf_area
USE ukca_trop_hetchem_mod,  ONLY: cal_n2o5, cal_loss1k

IMPLICIT NONE

INTEGER, INTENT(IN) :: nr
INTEGER, INTENT(IN) :: n_pnts

REAL, INTENT(IN)    :: tc      (n_pnts)  ! Temperature (K)
REAL, INTENT(IN)    :: m       (n_pnts)  ! Air density (molecules cm-3)
REAL, INTENT(IN)    :: h2o     (n_pnts)  ! Water concentration (molec cm-3)
REAL, INTENT(IN)    :: o2      (n_pnts)  ! Oxygen concn (molec cm-3)
REAL, INTENT(IN)    :: clw(n_pnts)       ! Cloud liquid water
REAL, INTENT(IN)    :: fcloud(n_pnts)    ! Cloud fraction
! Dissolved fraction, allowing for jpeq+1 ions:
REAL, INTENT(IN)    :: fdiss(n_pnts, jpspec, jpeq+1)
REAL, INTENT(IN)    :: rho_air (n_pnts)  ! Air density (kg m-3)
REAL, INTENT(IN)    :: rh_frac (n_pnts)  ! RH (fraction: 0.0 - 0.999)

! Mass mixing ratios (kg kg-1) or aerosol number (m-3) of CLASSIC aerosols.
REAL, INTENT(IN)    :: so4_aitken    (n_pnts)  ! kg kg-1
REAL, INTENT(IN)    :: so4_accum     (n_pnts)
REAL, INTENT(IN)    :: soot_fresh    (n_pnts)
REAL, INTENT(IN)    :: soot_aged     (n_pnts)
REAL, INTENT(IN)    :: ocff_fresh    (n_pnts)
REAL, INTENT(IN)    :: ocff_aged     (n_pnts)
REAL, INTENT(IN)    :: biogenic      (n_pnts)
REAL, INTENT(IN)    :: sea_salt_film (n_pnts)  ! m-3
REAL, INTENT(IN)    :: sea_salt_jet  (n_pnts)

LOGICAL, INTENT(IN) :: stratflag     (n_pnts)  ! True for cells in stratosphere

! 1-D pH array to be used in calculating rate constants
REAL, INTENT(IN) :: H_plus_1d_arr(n_pnts)

REAL, INTENT(OUT)   :: rc  (n_pnts,nr)   ! Gas phase rate constants (s-1)
REAL, INTENT(OUT)   :: hrc (n_pnts,jphk) ! Heterogeneous rate constants (s-1)

! Local variables

INTEGER, PARAMETER :: npdep = 10        ! No. press-dep reactions
REAL,    PARAMETER :: o2frac= 0.20948   ! o2 molar fraction in air
REAL,    PARAMETER :: n2frac= 0.78084   ! n2 molar fraction in air
REAL,    PARAMETER :: qcl_min = 1.0e-12 ! do calculations when qcl > qcl_min

! Indices for soluble species
INTEGER, PARAMETER :: ih_o3   =  4
INTEGER, PARAMETER :: ih_hno3 = 10
INTEGER, PARAMETER :: ih_h2o2 = 11
INTEGER, PARAMETER :: ih_so2  = 45

!
INTEGER :: i, j, k                     ! Loop counts
INTEGER :: r                           ! Reaction type
INTEGER :: nklo, nkhi, nres ! Reaction nos of klo, khi and location
!                                   of result for pressure-dependent reacs
!
REAL :: tfc = 1.0 / 300.0  ! Used in rates of the type (T/300)**n
REAL :: ar                 ! Arrhenius pre-exponential factor
REAL :: n                  ! Temperature power (K)
REAL :: ea                 ! Activation energy (as -Ea/R, units: K)
REAL :: fac1               ! Temporary store
REAL :: fac2               ! Temporary store
REAL :: fac3               ! Temporary store
REAL :: brn                ! Temporary store
REAL :: f                  ! Appropriate value of Fc factor
REAL :: fc                 ! Fc factor for R29
REAL :: fac1v(n_pnts)      ! Temporary store vector
REAL :: faq(n_pnts,jpspec) ! Total dissolved fraction
REAL :: vr                 ! Volume ratio of cloud water
!
! Surface areas from CLASSIC aerosols (cm2 cm-3)
REAL :: sa_so4_ait    (n_pnts) ! ammonium sulphate Aitken
REAL :: sa_so4_acc    (n_pnts) ! ammonium sulphate accumulation
REAL :: sa_bc_fresh   (n_pnts) ! black carbon (soot) fresh
REAL :: sa_bc_aged    (n_pnts) ! black carbon (soot) aged
REAL :: sa_ocff_fresh (n_pnts) ! OCFF fresh
REAL :: sa_ocff_aged  (n_pnts) ! OCFF aged
REAL :: sa_soa        (n_pnts) ! biogenic aerosol (climatology)
REAL :: sa_ss_film    (n_pnts) ! sea-salt film
REAL :: sa_ss_jet     (n_pnts) ! sea-salt jet

! Wet radii for CLASSSIC aerosols (cm)
REAL :: wr_so4_ait    (n_pnts)
REAL :: wr_so4_acc    (n_pnts)
REAL :: wr_bc_fresh   (n_pnts)
REAL :: wr_bc_aged    (n_pnts)
REAL :: wr_ocff_fresh (n_pnts)
REAL :: wr_ocff_aged  (n_pnts)
REAL :: wr_soa        (n_pnts)
REAL :: wr_ss_film    (n_pnts)
REAL :: wr_ss_jet     (n_pnts)

! Uptake coefficients (dimensionless) of N2O5 on different aerosol types,
! parameterised depending on T & RH.
REAL :: gamma_n2o5_so4  (n_pnts)   ! ammonium sulphate (Aitken & accumulation)
REAL :: gamma_n2o5_bc   (n_pnts)   ! black carbon (fresh & aged)
REAL :: gamma_n2o5_ocff (n_pnts)   ! organic carbon fossil fuel (fresh & aged)
REAL :: gamma_n2o5_soa  (n_pnts)   ! biogenic secondary organic aerosol
REAL :: gamma_n2o5_ss   (n_pnts)   ! sea-salt aerosol (film & jet)

! Temporary variables for SO2 oxidation
REAL :: arg2
REAL :: zfc
REAL :: zr
REAL :: zi(n_pnts)
REAL :: zo(n_pnts)
REAL :: alpha(n_pnts)      ! factor for DMS oxidation rates
REAL :: t300(n_pnts)       ! T/300 factor
REAL :: at(7)              ! Factors in SO2 oxidation rates

LOGICAL, SAVE :: first = .TRUE.
LOGICAL, SAVE :: first_pass = .TRUE.

INTEGER :: errcode                             ! Error code for ereport
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMCO_RAQ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! OMP CRITICAL will only allow one thread through this code at a time,
! while the other threads are held until completion.
!$OMP CRITICAL (init_ukca_chemco_raq)
IF (first_pass) THEN
  IF (first) THEN
    ! Set up derived type arrays rk (rate coefficients of thermal reactions)
    ! and pdep (data to calculate pressure-dependent reactions) on first call.
    CALL ukca_chemco_raq_init(nr, npdep)

    ! Set indices for heterogeneous reactions if used. For the moment
    ! only the heterogeneous hydrolysis of N2O5 on CLASSIC aerosols
    ! is considered, but other reactions can be added in the future.
    IF (ukca_config%l_ukca_classic_hetchem) THEN
      in2o5    = 1
      in2o5_ss = 2
    END IF

    first = .FALSE.
  END IF  ! If first = .TRUE.
  first_pass = .FALSE.
END IF    ! IF first_pass = .TRUE.
!$OMP END CRITICAL (init_ukca_chemco_raq)

! Initialise reaction rates
rc  (:,:) = 0.0
hrc (:,:) = 0.0

! Set up rate coefficient data type
! Rate constant k = A.(T/300)**n.exp(E/T)
! E is -Ea/R, where Ea is the activation energy and R the ideal gas constant
! Hence Ea has units of K.
! r is the reaction type:
!  0 - no reaction
!  1 - k = A  (temperature-independent)
!  2 - k = A * (T/300)**n
!  3 - k = A * exp(E/T)
!  4 - k = A * (T/300)**n * exp(E/T)
DO i = 1, nr
  r  = rk(i)%r
  ar = rk(i)%a
  n  = rk(i)%n
  ea = rk(i)%e
  IF (r == 1) THEN
    rc(:,i) = ar
  ELSE IF (r == 2) THEN
    rc(:,i) = ar * ((tc(:)*tfc) ** n)
  ELSE IF (r == 3) THEN
    rc(:,i) = ar * EXP(ea / tc(:))
  ELSE IF (r == 4) THEN
    rc(:,i) = ar * ((tc(:)*tfc) ** n) * EXP(ea / tc(:))
  END IF

END DO


! Pressure-dependent reactions

DO k = 1, npdep
  nklo = pdep(k)%klo                  ! Reaction no. for klo
  nkhi = pdep(k)%khi                  ! Reaction no. for khi
  nres = pdep(k)%res                  ! Reaction no. for result
  f    = pdep(k)%fc
  IF (f > 0.01) THEN
    fc = f
    DO j = 1, n_pnts
      rc(j,nklo) = rc(j,nklo) * m(j)          ! klo * M
      brn        = 0.75 - 1.27 * LOG10(fc)
      fac1       = rc(j,nklo) / rc(j,nkhi)    ! khi
      fac2       = rc(j,nklo) / (1.0 + fac1)
      fac3       = 1.0 + ((LOG10(fac1) / brn) ** 2)
      rc(j,nres) = fac2 * (fc ** (1.0 / fac3))
    END DO
  ELSE ! In stochem_chemistry.mf77 f was set to zero for R29
       ! to flag this reaction and recalculate f later.
       ! Following Atkinson (2004), f changed to 0.35, but
       ! we leave the code below in case of any new change.
    DO j = 1, n_pnts
      fc         = EXP(-tc(j)/250.0) + EXP(-1050.0/tc(j)) ! For R29
      rc(j,nklo) = rc(j,nklo) * m(j)          ! klo * M
      brn        = 0.75 - 1.27 * LOG10(fc)
      fac1       = rc(j,nklo) / rc(j,nkhi)    ! khi
      fac2       = rc(j,nklo) / (1.0 + fac1)
      fac3       = 1.0 + ((LOG10(fac1) / brn) ** 2)
      rc(j,nres) = fac2 * (fc ** (1.0 / fac3))
    END DO
  END IF
END DO


! Do remaining calculations for certain reactions. This section
! needs to be edited by hand

! Reaction   7 : O(1D) + M
rc(:, 7) = (o2frac*rc(:,6) + n2frac*rc(:,7)) * m(:)

! Reaction  35 : OH + HNO3
fac1v(:) = rc(:,51) * m(:)
rc(:,35) = rc(:,35) + fac1v(:) / (1.0 + fac1v(:) / rc(:,50))

! Reaction  36 : HO2 + HO2 (+ M, H2O)
rc(:,36) = (rc(:,36) + rc(:,37) * m(:)) *                                      &
     (1.0 + rc(:,38) * h2o(:))

! Reactions  61 and 62 : CH3O2 + CH3O2
! (use both branches in UKCA_DERIV)
rc(:,62) = rc(:,62) - rc(:,61)

! Reaction  70 : OH + CO
rc(:,70) = rc(:,70) + (rc(:,69) * m(:))

! Reactions  74 and 80 : CH3O2 + CH3COO2
! (consider both of them in UKCA_DERIV)
fac1v(:) = rc(:,74) / (1.0 + rc(:,74))
rc(:,74) = rc(:,80) * (1.0-fac1v(:))  ! -> 2 HCHO + O2
rc(:,80) = rc(:,80) * fac1v(:)        ! -> HCHO + HO2 + CH3O2 + CO2 + O2

! Reaction 89 and 94: CH3COCH3 + OH = CH3COCH2O2 + H2O
rc(:,94) = rc(:,89) + rc(:,94)

! Reaction   1 : O + O2 (+ M) = O3 + M
rc(:, 1) = rc(:, 1) * m(:) * o2(:)

! -------------------------------------------------------------
! Heterogeneous reactions on the surface of CLASSIC aerosols
!
IF (ukca_config%l_ukca_classic_hetchem) THEN

  IF (.NOT. (l_um_calc_surf_area)) THEN
    errcode = 1
    cmessage =                                                                 &
      'UM support unavailable for heterogenous chemistry on CLASSIC particles'
    CALL ereport(RoutineName, errcode, cmessage)
  END IF

  ! Call routine to calculate aerosol surface areas and wet radii
  CALL calc_surf_area (n_pnts,                                                 &
                    rho_air       (1:n_pnts), rh_frac      (1:n_pnts),         &
                    so4_aitken    (1:n_pnts), so4_accum    (1:n_pnts),         &
                    soot_fresh    (1:n_pnts), soot_aged    (1:n_pnts),         &
                    ocff_fresh    (1:n_pnts), ocff_aged    (1:n_pnts),         &
                    biogenic      (1:n_pnts),                                  &
                    sea_salt_film (1:n_pnts), sea_salt_jet (1:n_pnts),         &
                    ! Output arguments: surface areas (cm2/cm3)
                    sa_so4_ait    (1:n_pnts), sa_so4_acc   (1:n_pnts),         &
                    sa_bc_fresh   (1:n_pnts), sa_bc_aged   (1:n_pnts),         &
                    sa_ocff_fresh (1:n_pnts), sa_ocff_aged (1:n_pnts),         &
                    sa_soa        (1:n_pnts),                                  &
                    sa_ss_film    (1:n_pnts), sa_ss_jet    (1:n_pnts),         &
                    ! Output arguments: wet radii (cm)
                    wr_so4_ait    (1:n_pnts), wr_so4_acc   (1:n_pnts),         &
                    wr_bc_fresh   (1:n_pnts), wr_bc_aged   (1:n_pnts),         &
                    wr_ocff_fresh (1:n_pnts), wr_ocff_aged (1:n_pnts),         &
                    wr_soa        (1:n_pnts),                                  &
                    wr_ss_film    (1:n_pnts), wr_ss_jet  (1:n_pnts))

  ! Calculate the rate of the following reaction if it is used:
  !    N2O5 (CLASSIC aerosols) + H2O --> 2 HNO3
  !
  ! -------------------------------------------------------------------
  ! Step 1: get uptake coefficients for N2O5 (gamma_n2o5) by calling the
  ! function cal_n2o5 (aerotype, nbox, temp, rh) separately for each
  ! aerosol type. This function forms part of a module which was
  ! originally introduced for GLOMAP-mode: ukca_trop_hetchem_mod.
  ! The following analogy is established here for some of the aerosol
  ! types present in the two aerosol schemes:
  !
  !   CLASSIC aerosols                           GLOMAP-mode (aerotype)
  ! -------------------------------------       ------------------------
  !  ammonium sulphate          (NH4)2-SO4 <-->    SO4      (cp_su)
  !  black carbon               (BC)       <-->    BC       (cp_bc)
  !  organic carbon fossil fuel (OCFF)     <-->    OC       (cp_oc)
  !  biogenic secondary aerosol (BSOA)     <-->    SO       (cp_so)
  !  sea-salt                   (SS)       <-->    Cl       (cp_cl)
  !
  ! This analogy is imperfect. As an example, there is no separate sulphate
  ! component in CLASSIC so ammonium sulphate has to be used instead.
  ! Note also that the following aerosol types have not been used for the
  ! calculation of the heterogeneous rate, but the code could be extended
  ! in the future to consider some of them:
  !
  ! * Dust: Not too relevant over the regional domain (UK and part of Western
  !   Europe) for which the RAQ chemistry scheme is used at the moment. It
  !   could be used for larger domains in the future. In order to do that it
  !   would be necessary to work out a formulation for the surface area of
  !   dust, since unlike the other aerosol types it does not follow a
  !   log-normal distribution.
  !
  ! * Biomass burning (BB) aerosol: Again not too relevant over the current
  !   domain for which RAQ is used. Note that the function cal_n2o5 cannot
  !   use BB aerosol directly. Before extending the code for BB, it would be
  !   necessary to check how this aerosol type can be apportioned to both
  !   BC and OC, and then include separate calls to cal_n2o5 for them.
  !
  ! * Ammonium nitrate: This species is not included because it does not
  !   favour the heterogeneous hydrolysis of N2O5. Note that it is not
  !   considered by the function cal_n2o5 either.
  !
  ! Step 2: calculate the final rate by calling the function
  ! cal_loss1k (nbox, sarea, wetdp, aird, mode_gamma, temp, molec_wt_s),
  ! which is also part of the module ukca_trop_hetchem_mod, for
  ! each aerosol type and mode.
  ! The units of the arguments passed here are: sarea --> cm2/cm3,
  ! wetdp --> cm, aird --> cm-3, temp --> K, molec_wt_s --> g/mol.
  ! The other two arguments (nbox and mode_gamma) are dimensionless.
  !
  ! Note that it is necessary to pass different surface areas and
  ! wet radii for each aerosol type and mode, while the gamma_n2o5
  ! value is the same for all modes belonging to the same aerosol type.
  ! If a given aerosol type is not present then the corresponding
  ! surface area will be filled with zeros and cal_loss1k will
  ! also return zeros.
  !
  ! Note that reaction on sea salt is treated independently. This is
  ! because HNO3 formed in the heterogeneous reaction of N2O5 on sea
  ! salt is unlikely to contribute to PM2.5 and is treated separately.

  IF (in2o5 > 0) THEN

    ! The rh_frac values are limited to < 1.0. The function
    ! cal_n2o5 converts rh_frac to RH expressed as percentage
    ! before doing the calculations. Note that gamma_n2o5_ss
    ! is determined in "IF (in2o5_ss > 0)"
    gamma_n2o5_so4  (:) = cal_n2o5 ( cp_su, n_pnts, tc(:), rh_frac(:) )
    gamma_n2o5_bc   (:) = cal_n2o5 ( cp_bc, n_pnts, tc(:), rh_frac(:) )
    gamma_n2o5_ocff (:) = cal_n2o5 ( cp_oc, n_pnts, tc(:), rh_frac(:) )
    gamma_n2o5_soa  (:) = cal_n2o5 ( cp_so, n_pnts, tc(:), rh_frac(:) )

    ! Add contribution from SO4 Aitken
    hrc (:, in2o5) = cal_loss1k (n_pnts,                                       &
                     sa_so4_ait (:), wr_so4_ait (:), m (:),                    &
                     gamma_n2o5_so4 (:), tc(:), m_n2o5 )

    ! Add contribution from SO4 accumulation
    hrc (:, in2o5) = hrc (:, in2o5) + cal_loss1k ( n_pnts,                     &
                     sa_so4_acc (:), wr_so4_acc (:), m (:),                    &
                     gamma_n2o5_so4 (:), tc(:), m_n2o5 )

    ! Add contribution from BC fresh
    hrc (:, in2o5) = hrc (:, in2o5) + cal_loss1k ( n_pnts,                     &
                     sa_bc_fresh (:), wr_bc_fresh (:), m (:),                  &
                     gamma_n2o5_bc (:), tc(:), m_n2o5 )

    ! Add contribution from BC aged
    hrc (:, in2o5) = hrc (:, in2o5) + cal_loss1k ( n_pnts,                     &
                     sa_bc_aged (:), wr_bc_aged (:), m (:),                    &
                     gamma_n2o5_bc (:), tc(:), m_n2o5 )

    ! Add contribution from OCFF fresh
    hrc (:, in2o5) = hrc (:, in2o5) + cal_loss1k ( n_pnts,                     &
                     sa_ocff_fresh (:), wr_ocff_fresh (:), m (:),              &
                     gamma_n2o5_ocff (:), tc(:), m_n2o5 )

    ! Add contribution from OCFF aged
    hrc (:, in2o5) = hrc (:, in2o5) + cal_loss1k ( n_pnts,                     &
                     sa_ocff_aged (:), wr_ocff_aged (:), m (:),                &
                     gamma_n2o5_ocff (:), tc(:), m_n2o5 )

    ! Add contribution from BSOA
    hrc (:, in2o5) = hrc (:, in2o5) + cal_loss1k ( n_pnts,                     &
                     sa_soa (:), wr_soa (:), m (:),                            &
                     gamma_n2o5_soa (:), tc(:), m_n2o5 )

    ! Set the heterogeneous rates to zero in the stratosphere to
    ! avoid model instability.
    WHERE (stratflag) hrc (:, in2o5)  = 0.0

  END IF  ! in2o5 > 0

  IF (in2o5_ss > 0) THEN

    gamma_n2o5_ss   (:) = cal_n2o5 ( cp_cl, n_pnts, tc(:), rh_frac(:) )

    ! Add contribution from sea-salt film
    hrc (:, in2o5_ss) = cal_loss1k ( n_pnts,                                   &
                        sa_ss_film (:), wr_ss_film (:), m (:),                 &
                        gamma_n2o5_ss (:), tc(:), m_n2o5 )

    ! Add contribution from sea-salt jet
    hrc (:, in2o5_ss) = hrc (:, in2o5_ss) + cal_loss1k ( n_pnts,               &
                        sa_ss_jet (:), wr_ss_jet (:), m (:),                   &
                        gamma_n2o5_ss (:), tc(:), m_n2o5 )

    ! Set the heterogeneous rates to zero in the stratosphere to
    ! avoid model instability.
    WHERE (stratflag) hrc (:, in2o5_ss)  = 0.0

  END IF  ! in2o5_ss > 0
END IF    ! l_ukca_classic_hetchem = .TRUE.


! Finally do the additional reactions for aerosol chemistry if it is on
IF (ukca_config%l_ukca_raqaero) THEN
  ! Calculate the total dissolved fraction by adding the fraction of the
  ! species in each dissolved state
  faq(:,:) = fdiss(:,:,1) + fdiss(:,:,2)

  ! DMS oxidation based on Pham et al (1995).
  ! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 100, NO. D12, PP. 26,061-26,092, 1995.
  ! Reaction rates and products adopted from ukca_chem_tropisop.F90 at vn10.2.
  !
  ! R184: OH-initiated chemistry following hydrogen abstraction produces SO2
  !       as the sole sulphur compound.
  ! R185: OH-initiated chemistry following OH addition produces both SO2
  !       and DMSO. The latter is oxidised to MSA (see R187).
  ! R186: NO3-initiated chemistry produces SO2.
  !
  ! RC(184) %DMS     +OH      =SO2     +MeOO    +HCHO
  rc(:,184) =  9.60e-12*EXP(-240.0/tc(:))
  ! RC(185) %DMS     +OH      =0.6 SO2     +0.4 DMSO    +MeOO
  alpha(:)=1.106e-31*EXP(7460.0/tc(:))*m(:)
  rc(:,185) = 3.04e-12*EXP(350.00/tc(:))*alpha(:)/(1 + alpha(:))
  ! RC(186) %DMS     +NO3     =SO2     +HONO2   +MeOO    +HCHO
  rc(:,186) = 1.90e-13*EXP(500.00/tc(:))
  ! RC(187) %DMSO    +OH      =0.6 SO2     +0.4 MSA
  rc(:,187) = 5.80e-11

  ! Dry oxidation of SO2 by OH adapted from ukca_chemco.F90 at vn10.2
  ! RC(188) %SO2     +OH      =H2SO4   +HO2 Rate data from NASA/JPL.
  at(1) = 0.60
  at(2) = 3.00e-31
  at(3) = -3.30
  at(4) = 0.0
  at(5) = 1.50e-12
  at(6) = 0.00
  at(7) = 0.0

  t300(:) = tc(:)/300.0
  zo(:) = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
  zi(:) = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

  DO i=1,n_pnts
    IF (zo(i) < peps ) THEN
      rc(i,188) = zi(i)
    ELSE IF ( zi(i) < peps ) THEN
      rc(i,188) = zo(i)
    ELSE
      IF ( at(1) <= 1.0 ) THEN
        zfc = at(1)
      ELSE
        zfc = EXP( -tc(i)/at(1) )
      END IF
      zr = zo(i) / zi(i)
      arg2 = 1.0/(1.0+((LOG10(zr))**2))
      rc(i,188) = (zo(i)/(1.0+zr))*EXP(arg2*LOG(zfc))
      ! reduce reaction rate to account for dissolved fraction of SO2
      rc(i,188) = rc(i,188) * (1.0 - (faq(i,ih_so2) * fcloud(i)))
    END IF
  END DO

  !
  !          Calculate aqueous rates and reduce rates due to aqueous fraction
  !          as done in asad_hetero at vn10.2
  !          ----------------------------------------------------------------

  DO i=1,n_pnts
    ! Only do these calculations for grid cells that have a minimum amount of
    ! water
    IF (clw(i) > qcl_min) THEN

      ! convert clw in kg/kg to volume ratio
      vr = clw(i)*m(i)*m_air*(1e6/avogadro)/rho_water

      ! Reduce rates of reactions H2O2+OH (R31) and HNO3+OH (R35) to take
      ! account of the dissolved fraction of H2O2 and HNO3, respectively.
      ! R31: H2O2 + OH
      rc(i,31) = rc(i,31)*                                                     &
        (1.0 - (fdiss(i,ih_h2o2,1)+fdiss(i,ih_h2o2,2))*fcloud(i))
      ! R35: HNO3 + OH
      rc(i,35) = rc(i,35)*                                                     &
        (1.0 - (fdiss(i,ih_hno3,1)+fdiss(i,ih_hno3,2))*fcloud(i))


      ! Aqueous phase chemistry of SO2 is given by reactions R190, 191 & 197.
      ! Their rates have been optimised from Table 5 of Kreidenweis et al.,
      ! J. Geophys. Res. VOL. 108, NO. D7, 4213, doi:10.1029/2002JD002697, 2003:
      !   R190: SO3-- + O3(aq)   => S(VI)
      !   R191: HSO3- + H2O2(aq) => S(VI)
      !   R197: HSO3- + O3(aq)   => S(VI)
      ! where S(VI) = SO4-- + HSO4-.
      !
      ! Note that the RAQ-AERO mechanism uses the reactant SO2 (instead of
      ! SO3-- and HSO3-), while the products are disregarded because it is
      ! assumed they are rained out.

      ! R190 : SO3-- + O3(aq) => SO4-- [Kreidenweis et al. (2003), optimised]
      rc(i,190) = 7.43e+16*EXP(-5280.0/tc(i))*fcloud(i)*                       &
                     fdiss(i,ih_so2,3)*fdiss(i,ih_o3,1)*                       &
                     1000.0/(avogadro*vr)

      ! R191:  HSO3- + H2O2(aq) => SO4--  [Kreidenweis et al. (2003), optimised]
      ! Replace Global pH with 1D values from H_plus array
      rc(i,191) = 2.1295e+14*EXP(-4430.0/tc(i))*                               &
        (H_plus_1d_arr(i)/(1.0 + 13.0*H_plus_1d_arr(i)))*fcloud(i)             &
         *fdiss(i,ih_so2,2)*fdiss(i,ih_h2o2,1)*1000.0/(avogadro*vr)

      ! R197: HSO3- + O3(aq) => SO4--  [Kreidenweis et al. (2003), optimised]
      rc(i,197) = 4.0113e+13*EXP(-5530.0/tc(i))*                               &
                   fcloud(i)*fdiss(i,ih_so2,2)*fdiss(i,ih_o3,1)*               &
                   1000.0/(avogadro*vr)

    ELSE

      ! If no clouds then aqueous rates are zero
      rc(i,190) = 0.0
      rc(i,191) = 0.0
      rc(i,197) = 0.0

    END IF ! clw > qcl_min
  END DO
  !
  ! Monoterpene reaction rates adapted from ukca_chemco at vn10.2.
  ! A fixed yield of 0.14 is considered for secondary organic aerosol (SOA),
  ! but this could be controlled in the future via namelist.
  ! RC(192) %Monoterp+OH      =0.14 Sec_Org
  rc(:,192) = 1.2e-11*EXP(444.0/tc(:))
  ! RC(193) %Monoterp+O3      =0.14 Sec_Org
  rc(:,193) = 1.01e-15*EXP(-732.0/tc(:))
  ! RC(194) %Monoterp+NO3     =0.14 Sec_Org
  rc(:,194) = 1.19e-12*EXP(490.0/tc(:))
  ! Heterogeneous reaction rates (zero as not currently in RAQ-AERO)
  ! RC(195)   %N2O5    +M       =2 HONO2
  rc(:, 195) = 0.0
  ! RC(196)   %HO2     +HO2     =H2O2
  rc(:, 196) = 0.0

END IF ! l_ukca_raqaero


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chemco_raq

END MODULE ukca_chemco_raq_mod
