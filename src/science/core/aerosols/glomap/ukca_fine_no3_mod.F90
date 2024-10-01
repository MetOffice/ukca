! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Calculates fine NO3/NH4 emissions - the same as ukca_prod_no3_fine
!    in emiss_ctl but instead called by UKCA_AERO_STEP within the
!    aerosol scheme
!
!  References:
!    Hauglustaine et al ACP 2014
!    Code from Samuel Remy (samuel.remy at lmd.jussieu.fr)
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and The Met Office.
!  See:  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ######################################################################
!
! Subroutine Interface:
MODULE ukca_fine_no3_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_FINE_NO3_MOD'

CONTAINS

SUBROUTINE ukca_fine_no3(nbox,nadvg,nbudaer,zghno3,dtz,rhoa,aird,wetdp,RH_clr,t,sm,   &
                         mfpa,nd,md,mdt,s0g,bud_aer_mas)

!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculates condensation HNO3 and NH3 gases onto existing particles
! using thermodynamic equilibrium theory
!-----------------------------------------------------------------------
!  Inputs
!  NBOX        : Number of grid boxes
!  nadvg       : Number of gas phase advected tracers
!  nbudaer     : Number of aerosol budget fields
!  ZGHNO3      : HNO3 uptake coefficient
!  RHOA        : Air density (kg/m3)
!  AIRD        : Number density of air (per cm3)
!  WETDP       : Geometric mean wet diameter for each mode (m)
!  RH_clr      : Relative humidity of clear-sky portion (dimensionless 0-1)
!  T           : Centre level temperature (K)
!  SM          : Grid box mass of air (kg)
!  MFPA        : Mean free path of air (m)
!  DTZ         : Competition (cond/nucl) time step (s)
!
!  Outputs
!  ND          : Aerosol ptcl number density for mode (cm^-3)
!  MDT         : Avg tot mass of aerosol ptcl in mode (particle^-1)
!  MD          : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!  S0G         : Partial masses of gas phase species (kg per gridbox)
!-----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: glomap_variables
USE ukca_mode_setup,         ONLY: nmodes, cp_no3, cp_nh4, cp_su,              &
                                   mode_ait_sol, mode_acc_sol, mode_cor_sol,   &
                                   mode_acc_insol
USE ukca_um_legacy_mod,      ONLY: pi, rgas => r
USE chemistry_constants_mod, ONLY: avc => avogadro, zboltz => boltzmann
USE ukca_setup_indices,      ONLY: mm_gas, mhno3, mnh3,                        &
                                   nmasprimntaitsol, nmasprimntaccsol,         &
                                   nmasprimntcorsol, nmasprimnhaitsol,         &
                                   nmasprimnhaccsol, nmasprimnhcorsol
USE parkind1,   ONLY: jpim, jprb      ! DrHook
USE yomhook,    ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! IN/OUT arguments
INTEGER, INTENT(IN)  :: nbox
INTEGER, INTENT(IN)  :: nadvg
INTEGER, INTENT(IN)  :: nbudaer
REAL,    INTENT(IN)  :: zghno3
REAL,    INTENT(IN)  :: dtz
REAL,    INTENT(IN)  :: wetdp(nbox,nmodes)
REAL,    INTENT(IN)  :: t(nbox)
REAL,    INTENT(IN)  :: rhoa(nbox)
REAL,    INTENT(IN)  :: aird(nbox)
REAL,    INTENT(IN)  :: sm(nbox)
REAL,    INTENT(IN)  :: mfpa(nbox)
REAL,    INTENT(IN)  :: RH_clr(nbox)
REAL, INTENT(IN OUT) :: nd(nbox,nmodes)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables%ncp)
REAL, INTENT(IN OUT) :: s0g(nbox,nadvg)
REAL, INTENT(IN OUT) :: bud_aer_mas(nbox,0:nbudaer)

! Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: ddplim1(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: rhocomp(:)
REAL,    POINTER :: sigmag(:)
REAL,    POINTER :: x(:)
REAL,    POINTER :: mm(:)
REAL,    POINTER :: mlo(:)
REAL,    POINTER :: mmid(:)
REAL,    POINTER :: mfrac_0(:,:)
REAL,    POINTER :: num_eps(:)

! Molar masses - these should be consistent with other UKCA definitions
REAL :: zmwnh3    !17.0e-3 !Kg/mol
REAL :: zmwhno3   !63.0e-3
REAL :: zmwnh4    !18.0e-3
REAL :: zmwno3    !62.0e-3
REAL :: zmwair

! Local arrays
REAL :: dmd_no3(nbox,nmodes)    ! NO3 md tendency
REAL :: dmd_nh4(nbox,nmodes)    ! NH4 md tendency
REAL :: ptno3(nbox,3)           ! NO3 molar tendency
REAL :: ptnh4(nbox,3)           ! NH4 molar tendency
REAL :: pthno3(nbox)            ! HNO3 molar tendency
REAL :: ptnh3(nbox)             ! NH3 molar tendency
REAL :: zdghno3(nbox)           ! Diffusion coefficient
REAL :: zghno3_dd(nbox)         ! Diffusion coefficient for dust
REAL :: no3_mol_mode(nbox,3)    ! NO3 molar concentration
REAL :: nh4_mol_mode(nbox,3)    ! NH4 molar concentration
REAL :: so4_mol_mode(nbox,3)    ! SO4 molar concentration
REAL :: nh4inso4(nbox,3)        ! NH4 associated with SO4
REAL :: nh4inno3(nbox,3)        ! NH4 associated with NO3
REAL :: relrate(nbox,3)         ! Relative uptake rate
REAL :: mode_khno3(nbox,2)      ! Uptake coefficient
REAL :: ndnew(nbox,nmodes)      ! Updated nd
REAL :: sixovrpix(nmodes)

! Variables for main primary production calculation
REAL :: zdg0, zkn, nk1, zrh1580
REAL :: t_inv, t_log, rh_local, rh_inv, aird_local
REAL :: drh, kps, kpl1, kpl2, kpl3, kpl
REAL :: no3_mol_tot, nh4_mol_tot, so4_mol_tot, hno3_mol_tot, nh3_mol_tot
REAL :: ztn, zts, zta, zso4, ztadisp, ztnta
REAL :: nh4plus, nh4inso4tot, nh3inso4
REAL :: zwrk1, zwrk2, nh3_mol_up, hno3_mol_eq, nh3_mol_eq
REAL :: hno3_mol_tdep, nh3_mol_tdep, no3_mol_tdep, nh4_mol_tdep
REAL :: ptno3_tot, pt_no3_mol
REAL :: ptnh4_tot,zhno3_test, zthno3_tmp, znh3_test, ztnh3_tmp
REAL :: dtno3, dtno3lim, dtnh4, dtnh4lim
REAL :: totrate, tau_tot, tau_fac

! Variables for tendency conversion from kg/kg/s to model units
REAL :: test_tot_mdt, test_nd_min, vol2num, mdtmin, meanfac

! dummy variables
INTEGER :: i, icp, imode, updmode, idx

! Parameters (consistent with fine NO3 production in emissions)
REAL, PARAMETER :: zrgas = 8.314         ! Gas constant [J/mol/K]
REAL, PARAMETER :: zmgas = 29.0e-3        ! Molecular weight dry air [Kg mol-1]
REAL, PARAMETER :: zdmolec = 4.5e-10     ! Molec diameter air [m]
REAL, PARAMETER :: zratiohno3 = (63.0e-3+29e-3) / 63.0e-3
REAL, PARAMETER :: zghno315_dd = 1.0e-5
REAL, PARAMETER :: zghno330_dd = 1.0e-4
REAL, PARAMETER :: zghno370_dd = 6.0e-4
REAL, PARAMETER :: zghno380_dd = 1.05e-3

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_FINE_NO3'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------------
! 0.0 Setup constants

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
ddplim1     => glomap_variables%ddplim1
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
rhocomp     => glomap_variables%rhocomp
sigmag      => glomap_variables%sigmag
x           => glomap_variables%x
mm          => glomap_variables%mm
mlo         => glomap_variables%mlo
mmid        => glomap_variables%mmid
mfrac_0     => glomap_variables%mfrac_0
num_eps     => glomap_variables%num_eps

! Molar masses
zmwnh3  = mm_gas(mnh3)
zmwhno3 = mm_gas(mhno3)
zmwnh4  = mm(cp_nh4)
zmwno3  = mm(cp_no3)
zmwair  = avc * zboltz / rgas ! Molar mass of dry air (kg/mol)

!---------------------------------------
! 1.0 Calculate the diffusion constant and mean free path
DO i = 1,nbox
  zdg0 = 3.0/(8.0*avc*rhoa(i)*(zdmolec**2.0))
  zdghno3(i) = zdg0 * SQRT(zrgas*zmgas*t(i)*zratiohno3/(2.0*pi))
END DO

!---------------------------------------
! 2.0 Calculate the modal uptake coefficient
mode_khno3(:,:) = 0.0
DO imode = mode_ait_sol, mode_acc_sol
  sixovrpix(imode) = 6.0/(pi*x(imode))
END DO

DO imode = mode_ait_sol, mode_acc_sol
  IF (mode(imode)) THEN
    idx=imode-mode_ait_sol+1
    meanfac=EXP(0.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
    DO i = 1,nbox
      IF (nd(i,imode) > num_eps(imode)) THEN
        ! Calculate the uptake coefficient
        zkn  = 2.0*mfpa(i) / (wetdp(i,imode)*meanfac)
        nk1  = (2.0*pi*(wetdp(i,imode)*meanfac)*zdghno3(i)) / (1.0+(4.0*zkn/   &
                    zghno3/3.0)*(1.0-0.47*zghno3/(1.0+zkn)))
        mode_khno3(i,idx) = mode_khno3(i,idx) + nd(i,imode)*1e6*nk1
      END IF
    END DO
  END IF ! mode
END DO !imode

! 2.1 Also include contribution from dust although note that at present
!     ageing of dust is not permitted and mass is simply transferred to
!     the equivalent soluble mode
!     Note that modal dust may not be on, in which case this is not used
imode = mode_acc_insol
zghno3_dd(:) = 0.0
IF (mode(imode)) THEN
  idx=mode_acc_sol-1
  meanfac=EXP(0.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
  !RH dependent HNO3 uptake coefficient based on Fairlie et al. 2010
  DO i = 1,nbox
    zrh1580 = RH_clr(i)
    IF (zrh1580 < 0.15) zrh1580 = 0.15
    IF (zrh1580 > 0.80) zrh1580 = 0.80
    IF (zrh1580 < 0.30 ) THEN
      zghno3_dd(i) = zghno315_dd + (zghno330_dd-zghno315_dd)/0.15*             &
                     (zrh1580-0.15)
    ELSE IF ((zrh1580 > 0.30) .AND. (zrh1580 < 0.70)) THEN
      zghno3_dd(i) = zghno330_dd + (zghno370_dd-zghno330_dd)/0.40*             &
                   (zrh1580-0.30)
    ELSE IF (zrh1580 > 0.70 ) THEN
      zghno3_dd(i) = zghno370_dd + (zghno380_dd-zghno370_dd)/0.10*             &
                   (zrh1580-0.70)
    END IF

    ! Calculate the uptake coefficient and add to accumulation soluble mode
    IF (nd(i,imode) > num_eps(imode)) THEN
      ! Calculate the uptake coefficient
      zkn  = 2.0*mfpa(i) / (wetdp(i,imode)*meanfac)
      nk1  = (2.0*pi*(wetdp(i,imode)*meanfac)*zdghno3(i)) / (1.0+(4.0*zkn/     &
              zghno3_dd(i)/3.0)*(1.0-0.47*zghno3_dd(i)/(1.0+zkn)))
      mode_khno3(i,idx) = mode_khno3(i,idx) + nd(i,imode)*1e6*nk1
    END IF
  END DO ! i = 1,nbox
END IF ! mode

! Set a minimum value for KHNO3 to avoid dividing by zero
DO imode = mode_ait_sol, mode_acc_sol
  idx=imode-mode_ait_sol+1
  DO i = 1,nbox
    IF (mode_khno3(i,idx) < 1.0e-40) mode_khno3(i,idx) = 1.0e-40
  END DO
END DO

!---------------------------------------
! 3. Calculate equilibrium concentrations of HNO3 and NH3 and apportion to
!    ammonium, nitrate, and sulfate
!    Scheme follows Hauglustaine (2014) Atmos. Chem. Phys., 14, 11031-11063

! Initialise arrays to hold concentration tendencies
no3_mol_mode(:,:)=0.0
nh4_mol_mode(:,:)=0.0
so4_mol_mode(:,:)=0.0
ptno3(:,:) = 0.0
ptnh4(:,:) = 0.0
pthno3(:) = 0.0
ptnh3(:) = 0.0
nh4inso4(:,:) = 0.0
nh4inno3(:,:) = 0.0
relrate(:,:) = 0.0

updmode = mode_ait_sol - 1 ! any remaining NH4 is emitted into Aitken mode

! Begin loop over gridcells - hopefully not striding too much
DO i = 1,nbox

  ! 3.1 Parameters related to temperature / humidity, loop over nbox
  t_inv = 1.0 / t(i)
  t_log = LOG(t(i))
  rh_local = RH_clr(i)
  aird_local = rhoa(i)
  IF (rh_local < 0.0) rh_local = 0
  IF (rh_local > 0.98) rh_local = 0.98
  rh_inv = 1.0 - rh_local

  ! 3.2 Equilibrium constant for [NH4.NO3] based on Mozurkewich, 1993
  drh  = EXP(723.7 * t_inv + 1.6954) * 1.0e-2
  kps  = EXP( 118.87 - 24084.0 * t_inv - 6.025 * t_log)
  kpl1 = EXP(-135.94 + 8763.0  * t_inv + 19.12 * t_log)
  kpl2 = EXP(-122.65 + 9969.0  * t_inv + 16.22 * t_log)
  kpl3 = EXP(-182.61 + 13875.0 * t_inv + 24.46 * t_log)
  kpl  = kps
  IF (rh_local >= drh ) THEN
    kpl = (kpl1-kpl2*rh_inv+kpl3*rh_inv**2.0) *rh_inv**1.75 * kpl
  END IF

  ! 3.3 Convert aerosol and gas fields to molar concentrations
  !     Aerosol md is in molecules #-1 (molecules of x per air molecule)
  !     Aerosol nd is in # cm-3 (molecules of air per cm3)
  !     avc is molecules of substance per mole
  !     1e15 = 1e6 (mol->nmol) * 1e9 (cm-3 to m-3)
  !     Upshot is nmol(substance) per m3
  no3_mol_tot = 0.0
  nh4_mol_tot = 0.0
  so4_mol_tot = 0.0
  DO imode = mode_ait_sol, mode_cor_sol
    idx = imode - mode_ait_sol + 1
    no3_mol_mode(i,idx) = nd(i,imode)*md(i,imode,cp_no3)*1e15/avc
    nh4_mol_mode(i,idx) = nd(i,imode)*md(i,imode,cp_nh4)*1e15/avc
    so4_mol_mode(i,idx) = nd(i,imode)*md(i,imode,cp_su)*1e15/avc
    no3_mol_tot = no3_mol_tot + no3_mol_mode(i,idx)
    nh4_mol_tot = nh4_mol_tot + nh4_mol_mode(i,idx)
    so4_mol_tot = so4_mol_tot + so4_mol_mode(i,idx)
  END DO
  hno3_mol_tot = s0g(i,mhno3)*aird_local*1.0e9/(zmwair*sm(i))
  nh3_mol_tot = s0g(i,mnh3)*aird_local*1.0e9/(zmwair*sm(i))
  ztn = hno3_mol_tot + no3_mol_tot
  zta = nh4_mol_tot + nh3_mol_tot
  zts = so4_mol_tot

  ! 3.4 SO4 state. Metzger et al. 2002
  zso4 = 2.0
  IF (zts > 0.5*zta) zso4 = 1.5
  IF (zts > zta) zso4 = 1.0

  ! 3.5 Neutralize SO4 using NH4 and if necessary NH3
  nh4plus = zta
  nh4inso4tot = 0.0
  DO imode = mode_ait_sol, mode_cor_sol
    idx = imode - mode_ait_sol + 1
    nh4inso4(i,idx) = MIN(nh4plus,zso4*so4_mol_mode(i,idx))
    nh4plus=MAX(0.0,nh4plus-nh4inso4(i,idx))
    nh4inso4tot = nh4inso4tot + nh4inso4(i,idx)
  END DO
  ztadisp = MAX(0.0,zta-nh4inso4tot)
  nh3inso4 = MAX(0.0,nh4inso4tot-nh4_mol_tot)

  ! 3.6 Determine whether NO3 is formed or disassociates
  ztnta = ztn * ztadisp
  ptno3_tot = 0.0
  IF (ztnta > kpl ) THEN
    ! 3.6.1 Equilibrium NH4.NO3 concentrations
    zwrk1  = (ztadisp+ztn)**2.0-4.0*(ztnta-kpl)
    zwrk1  = MAX(zwrk1,0.0)
    zwrk2  = 0.5*(ztadisp+ztn-SQRT(zwrk1))
    zwrk2  = MAX(zwrk2,0.0)
    zwrk2  = MIN(zwrk2,ztn)
    zwrk2  = MIN(zwrk2,ztadisp)

    ! 3.6.2 Calculate time dependence of uptake in Aitken
    !       and Accumulation modes = khno3 * number conc
    totrate = 0.0
    DO imode = mode_ait_sol, mode_acc_sol
      idx = imode - mode_ait_sol + 1
      totrate = totrate + mode_khno3(i,idx)
    END DO
    DO imode = mode_ait_sol, mode_acc_sol
      idx = imode - mode_ait_sol + 1
      relrate(i,idx) = mode_khno3(i,idx)/totrate
    END DO
    tau_tot = 1.0/totrate
    tau_fac = 1.0-EXP(-1.0*dtz/tau_tot)

    ! 3.6.3 Time dependence gas phase concentrations
    !       Firstly update NH3 to account for SO4 neutralisation
    !       Then calculate equilibrium gas concentrations
    !       Then add time-dependence
    nh3_mol_up   = MAX(nh3_mol_tot-nh3inso4, 0.0)
    hno3_mol_eq  = MAX(ztn-zwrk2, 0.0 )
    nh3_mol_eq   = MAX(ztadisp-zwrk2, 0.0)
    hno3_mol_tdep= MAX(hno3_mol_tot-tau_fac*(hno3_mol_tot-hno3_mol_eq),0.0)
    nh3_mol_tdep = MAX(nh3_mol_up-tau_fac*(nh3_mol_up-nh3_mol_eq),0.0)

    ! 3.6.4 Update aerosol phase - firstly in moles, then convert to MMRs
    no3_mol_tdep = MAX(ztn-hno3_mol_tdep, 0.0)
    nh4_mol_tdep = MAX(ztadisp-nh3_mol_tdep, 0.0)
    nh4plus = nh4_mol_tdep
    DO imode = mode_ait_sol, mode_cor_sol
      idx = imode - mode_ait_sol + 1
      pt_no3_mol = MAX(0.0,(no3_mol_tdep-no3_mol_tot))*relrate(i,idx)/dtz
      nh4inno3(i,idx) = MIN(nh4plus,no3_mol_mode(i,idx)+pt_no3_mol*dtz)
      nh4plus = MAX(0.0,nh4plus-nh4inno3(i,idx))
      ptno3(i,idx) = pt_no3_mol*zmwno3/(aird_local*1.0e9)
      ptno3_tot = ptno3_tot+ptno3(i,idx)
    END DO
    nh4inno3(i,updmode)=nh4inno3(i,updmode)+nh4plus

  ELSE
    ! 3.6.5 Remove all nitrate and set NH4(NO3) concentration to zero
    DO imode = mode_ait_sol, mode_cor_sol
      idx = imode - mode_ait_sol + 1
      ptno3(i,idx) = -1.0*no3_mol_mode(i,idx)*zmwno3/(aird_local*1.0e9*dtz)
      ptno3_tot = ptno3_tot + ptno3(i,idx)
    END DO
  END IF ! ztnta > kpl

  ! 3.7 Absorb remaining NH4 concentration differences int0
  !    total NH4 tendency
  ptnh4_tot = 0.0
  DO imode = mode_ait_sol, mode_cor_sol
    idx = imode - mode_ait_sol + 1
    ptnh4(i,idx) = (nh4inso4(i,idx)+nh4inno3(i,idx)-nh4_mol_mode(i,idx)) *     &
                    zmwnh4/(aird_local*1.0e9*dtz)
    ptnh4_tot = ptnh4_tot + ptnh4(i,idx)
  END DO

  ! 3.8 Update gas phase tendencies and check for negative masses
  pthno3(i) = -1.0*ptno3_tot*zmwhno3/zmwno3
  ptnh3(i)  = -1.0*ptnh4_tot*zmwnh3/zmwnh4

  zhno3_test = hno3_mol_tot*zmwhno3/aird_local/1.0e9 + pthno3(i)*dtz
  IF (zhno3_test < 0.0) THEN
    zthno3_tmp = pthno3(i)
    pthno3(i)  = -1.0*s0g(i,mhno3)*zmwhno3/(sm(i)*zmwair)/dtz
    zthno3_tmp = pthno3(i) -zthno3_tmp
    DO imode = mode_ait_sol, mode_acc_sol
      idx = imode - mode_ait_sol + 1
      ptno3(i,idx) = ptno3(i,idx)-zthno3_tmp*(zmwno3/zmwhno3)*relrate(i,idx)
    END DO
  END IF ! HNO3 test

  znh3_test = nh3_mol_tot*zmwnh3/aird_local/1.0e9 + ptnh3(i) * dtz
  IF (znh3_test < 0.0) THEN
    ztnh3_tmp = ptnh3(i)
    ptnh3(i) = -1.0*s0g(i,mnh3)*zmwnh3/(sm(i)*zmwair)/dtz
    ztnh3_tmp = ptnh3(i)-ztnh3_tmp
    DO imode = mode_ait_sol, mode_acc_sol
      idx = imode - mode_ait_sol + 1
      ptnh4(i,idx) = ptnh4(i,idx)-ztnh3_tmp*(zmwnh4/zmwnh3)*relrate(i,idx)
    END DO
  END IF ! NH3 test

END DO ! nbox

!---------------------------------------------------------------------------
! 4. Convert NO3 and NH4 mass tendencies to Md units from kg/kg/s
!    Then fill output mass tendency arrays
!    Currently removing all mass from modes if negative emissions
!    and removing number as a ratio of increment to total mass
!    to conserve diameter
ndnew(:,:) = 0.0
dmd_no3(:,:) = 0.0
dmd_nh4(:,:) = 0.0

DO imode = mode_ait_sol, mode_cor_sol
  idx = imode - mode_ait_sol + 1
  vol2num = (zboltz/rgas)*sixovrpix(imode)/(ddplim1(imode)**3)

  DO i = 1,nbox
    IF (nd(i,imode) > num_eps(imode)) THEN
      dtno3 = ptno3(i,idx)*(zmwair/zmwno3)*aird(i)/nd(i,imode)
      dtno3lim = -1.0*md(i,imode,cp_no3)/dtz
      dtnh4 = ptnh4(i,idx)*(zmwair/zmwnh4)*aird(i)/nd(i,imode)
      dtnh4lim = -1.0*md(i,imode,cp_nh4)/dtz

      dmd_no3(i,imode) = MAX(dtno3lim, dtno3)
      dmd_nh4(i,imode) = MAX(dtnh4lim, dtnh4)

      ! CHECK NUMBER 1
      ! If all aerosol mass is removed from mode (test_tot_mdt)
      ! then we must also remove number concentration and reset the mode
      test_tot_mdt = mdt(i,imode) + (dmd_nh4(i,imode)+dmd_no3(i,imode))*dtz
      mdtmin=mlo(imode)*0.001
      IF (test_tot_mdt < mdtmin) THEN
        DO icp = 1,ncp
          IF (component(imode,icp)) THEN
            md(i,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
          END IF
        END DO
        mdt(i,imode) = mmid(imode)
        nd(i,imode) = 0.0
      END IF

      ! CHECK NUMBER 2 - Emissions mode
      ! If adding NH4NO3 mass where current number concentrations
      ! are negligible we also add number (assuming the growth of
      ! nucleation size particles to form Aitken sizes)
      test_nd_min = (test_tot_mdt*mm(cp_su)/(avc*rhocomp(cp_su))) * vol2num
      IF ((nd(i,imode) < test_nd_min) .AND. (imode < mode_cor_sol)) THEN
        ndnew(i,imode) = test_nd_min
      ELSE
        ndnew(i,imode) = nd(i,imode)
      END IF

      ! Update budget aer mass (diagnostic output)
      IF ((imode == mode_ait_sol) .AND. (nmasprimntaitsol > 0)) THEN
        bud_aer_mas(i,nmasprimntaitsol) = bud_aer_mas(i,nmasprimntaitsol) +    &
                                          dmd_no3(i,imode)*ndnew(i,imode)*dtz
      END IF
      IF ((imode == mode_ait_sol) .AND. (nmasprimnhaitsol > 0)) THEN
        bud_aer_mas(i,nmasprimnhaitsol) = bud_aer_mas(i,nmasprimnhaitsol) +    &
                                          dmd_nh4(i,imode)*ndnew(i,imode)*dtz
      END IF
      IF ((imode == mode_acc_sol) .AND. (nmasprimntaccsol > 0)) THEN
        bud_aer_mas(i,nmasprimntaccsol) = bud_aer_mas(i,nmasprimntaccsol) +    &
                                          dmd_no3(i,imode)*ndnew(i,imode)*dtz
      END IF
      IF ((imode == mode_acc_sol) .AND. (nmasprimnhaccsol > 0)) THEN
        bud_aer_mas(i,nmasprimnhaccsol) = bud_aer_mas(i,nmasprimnhaccsol) +    &
                                          dmd_nh4(i,imode)*ndnew(i,imode)*dtz
      END IF
      IF ((imode == mode_cor_sol) .AND. (nmasprimntcorsol > 0)) THEN
        bud_aer_mas(i,nmasprimntcorsol) = bud_aer_mas(i,nmasprimntcorsol) +    &
                                          dmd_no3(i,imode)*ndnew(i,imode)*dtz
      END IF
      IF ((imode == mode_cor_sol) .AND. (nmasprimnhcorsol > 0)) THEN
        bud_aer_mas(i,nmasprimnhcorsol) = bud_aer_mas(i,nmasprimnhcorsol) +    &
                                          dmd_nh4(i,imode)*ndnew(i,imode)*dtz
      END IF

      ! Update md, nd and mdt
      md(i,imode,cp_no3) = MAX(0.0,md(i,imode,cp_no3)+dmd_no3(i,imode)*dtz)
      md(i,imode,cp_nh4) = MAX(0.0,md(i,imode,cp_nh4)+dmd_nh4(i,imode)*dtz)
      mdt(i,imode) = 0.0
      DO icp = 1,ncp
        IF (component(imode,icp)) THEN
          mdt(i,imode) = mdt(i,imode) + md(i,imode,icp)
        END IF
      END DO
      nd(i,imode) = ndnew(i,imode)

    ELSE

      ! Reset mode if nd < num_eps
      DO icp = 1,ncp
        IF (component(imode,icp)) THEN
          md(i,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
        END IF
      END DO
      mdt(i,imode) = mmid(imode)
      nd(i,imode) = 0.0

    END IF ! nd > num_eps

  END DO ! i=1,nbox
END DO ! imode

! Update gas concentrations
DO i = 1,nbox
  s0g(i,mhno3) = s0g(i,mhno3) + pthno3(i)*sm(i)*(zmwair/zmwhno3)*dtz
  s0g(i,mnh3) = s0g(i,mnh3) + ptnh3(i)*sm(i)*(zmwair/zmwnh3)*dtz
  IF (s0g(i,mhno3) < 0.0) s0g(i,mhno3) = 0.0
  IF (s0g(i,mnh3) < 0.0) s0g(i,mnh3) = 0.0
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_fine_no3
END MODULE ukca_fine_no3_mod
