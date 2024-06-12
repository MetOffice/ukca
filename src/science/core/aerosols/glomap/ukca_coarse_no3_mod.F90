! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Calculates coarse-NO3 emissions - the same as ukca_prod_no3_coarse
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
MODULE ukca_coarse_no3_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_COARSE_NO3_MOD'

CONTAINS

SUBROUTINE ukca_coarse_no3(nbox,nadvg,nbudaer,dtz,rhoa,aird,wetdp,RH_clr,t,    &
                           sm,mfpa,nd,md,mdt,s0g,bud_aer_mas)

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
USE ukca_mode_setup,         ONLY: nmodes, cp_nn, cp_du, cp_cl,                &
                                   mode_acc_sol, mode_cor_sol,                 &
                                   mode_acc_insol, mode_cor_insol
USE ukca_um_legacy_mod,      ONLY: pi, rgas => r
USE chemistry_constants_mod, ONLY: avc => avogadro, zboltz => boltzmann
USE ukca_setup_indices,      ONLY: mm_gas, mhno3,                              &
                                   nmascondnnaccsol, nmascondnncorsol

USE parkind1,   ONLY: jpim, jprb      ! DrHook
USE yomhook,    ONLY: lhook, dr_hook  ! DrHook


IMPLICIT NONE

! IN/OUT arguments
INTEGER, INTENT(IN)  :: nbox
INTEGER, INTENT(IN)  :: nadvg
INTEGER, INTENT(IN)  :: nbudaer
REAL,    INTENT(IN)  :: dtz
REAL,    INTENT(IN)  :: wetdp(nbox,nmodes)
REAL,    INTENT(IN)  :: t(nbox)
REAL,    INTENT(IN)  :: rhoa(nbox)
REAL,    INTENT(IN)  :: aird(nbox)
REAL,    INTENT(IN)  :: sm(nbox)
REAL,    INTENT(IN)  :: RH_clr(nbox)
REAL,    INTENT(IN)  :: mfpa(nbox)
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
REAL,    POINTER :: ddplim0(:)
REAL,    POINTER :: ddplim1(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: rhocomp(:)
REAL,    POINTER :: sigmag(:)
REAL,    POINTER :: x(:)
REAL,    POINTER :: mm(:)
REAL,    POINTER :: mmid(:)
REAL,    POINTER :: mfrac_0(:,:)
REAL,    POINTER :: num_eps(:)

! Local Molar weights
REAL :: zmwnacl   != 58.44e-3
REAL :: zmwhno3   != 63.0e-3
REAL :: zmwnano3  != 84.0e-3
REAL :: zmwcaco3  != 100.0e-3
REAL :: zmwair
REAL, PARAMETER :: zmwcano3_2  = 164.0e-3
REAL, PARAMETER :: zmwca = 40.0e-3
REAL, PARAMETER :: zratiohno3  = (63.0e-3+29e-3) / 63.0e-3
INTEGER, PARAMETER :: nssmodes = 2 ! # of sea-salt modes
INTEGER, PARAMETER :: ndumodes = 2 ! # of dust modes
INTEGER, PARAMETER :: nno3modes = 2 ! # of coarse-NO3 modes

! dummy variables
INTEGER :: i, icp, imode, idx

! Local multi-dimension arrays
REAL    :: seasalt_diam(nbox,nssmodes)
REAL    :: seasalt_mconc(nbox,nssmodes)
REAL    :: seasalt_nconc(nbox,nssmodes)
REAL    :: dust_diam(nbox,ndumodes)
REAL    :: dust_mconc(nbox,ndumodes)
REAL    :: dust_nconc(nbox,ndumodes)
REAL    :: hno3_mconc(nbox)
REAL    :: ptdu(nbox,ndumodes)    ! Dust  MMR tendency (kg/kg/s)
REAL    :: ptss(nbox,nssmodes)    ! Sea-salt MMR tendency (kg/kg/s)
REAL    :: pthno3(nbox)           ! HNO3 MMR tendency (kg/kg/s)
REAL    :: ptnn(nbox,nno3modes)   ! coarse NO3 MMR tendency (kg/kg/s)
REAL    :: dmd_nn(nbox,nmodes)    ! coarse NO3 md tendency
REAL    :: dmd_cl(nbox,nmodes)    ! SS md tendency

! Local variables used in aer_no3_ukcadu
INTEGER :: iirh, jss, jdu, jtab
REAL    :: rssgrowth_rhtab(12) = [ 1.0, 1.0, 1.0, 1.0, 1.442, 1.555,           &
                                 1.666, 1.799, 1.988, 2.131, 2.361, 2.876]
REAL    :: rrhtab(12) = [ 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, &
                        85.0, 90.0, 95.0 ]
REAL    :: zdg0, zdghno3, zrh1580
REAL    :: zghno3_dd, zghno3_ss, sumdu
REAL    :: zkn_dd, zhno3_test, zthno3_tmp
REAL    :: zkn_ss, zrhloc
REAL    :: zkhno3_ss(nssmodes)
REAL    :: zkhno3_dd(ndumodes)
REAL    :: zaccno3, zcorno3, dtno3, zfracm1, zfracm2
REAL    :: zdd1, zdd2, zss1, zss2, sstmp, ddtmp, dtno3old, mcf
REAL    :: sixovrpix(nmodes)
REAL    :: meanfac, tmpdiam
REAL, PARAMETER :: zghno315_dd = 1.0e-5
REAL, PARAMETER :: zghno330_dd = 1.0e-4
REAL, PARAMETER :: zghno370_dd = 6.0e-4
REAL, PARAMETER :: zghno380_dd = 1.05e-3
REAL, PARAMETER :: zghno315_ss = 1.0e-3
REAL, PARAMETER :: zghno330_ss = 1.0e-2
REAL, PARAMETER :: zghno370_ss = 6.0e-2
REAL, PARAMETER :: zghno380_ss = 1.05e-1
REAL, PARAMETER :: zrgas = 8.314            !J/mol/K
REAL, PARAMETER :: zmgas = 29.0e-3          !Molecular weight in Kg
REAL, PARAMETER :: zdmolec = 4.5e-10        !molec diameter in m
REAL, PARAMETER :: zavg = 6.02217e+23       ! Avogadro number [mol-1]
REAL, PARAMETER :: pfrac_ca=0.05            ! Fraction of dust that is Ca2+
LOGICAL :: do_dust

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_COARSE_NO3'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------------
! 0.0 Setup constants

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
ddplim0     => glomap_variables%ddplim0
ddplim1     => glomap_variables%ddplim1
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
rhocomp     => glomap_variables%rhocomp
sigmag      => glomap_variables%sigmag
x           => glomap_variables%x
mm          => glomap_variables%mm
mmid        => glomap_variables%mmid
mfrac_0     => glomap_variables%mfrac_0
num_eps     => glomap_variables%num_eps

! Molar masses
zmwhno3  = mm_gas(mhno3)
zmwnacl  = mm(cp_cl)
zmwnano3 = mm(cp_nn)
zmwcaco3 = mm(cp_du)
zmwair   = avc * zboltz / rgas ! Molar mass of dry air (kg/mol)
DO i = 1,nmodes
  sixovrpix(:) = 6.0/(pi*x(:))
END DO

! HNO3 concentration to MMR
hno3_mconc(:) = 0.0
DO i = 1,nbox
  hno3_mconc(i) = s0g(i,mhno3)*zmwhno3/(sm(i)*zmwair)
END DO

! Determine if modal dust is used in this setup
do_dust = ANY(component(mode_acc_insol:mode_cor_insol,cp_du))

!---------------------------------------------------------------------------
! 1.0 Sea-salt into local arrays
! Initialise arrays to zero
seasalt_mconc(:,:)=0.0
seasalt_nconc(:,:)=0.0
seasalt_diam(:,:)=0.0

DO imode = mode_acc_sol, mode_cor_sol
  IF (mode(imode)) THEN
    idx=imode-mode_acc_sol+1
    meanfac=EXP(0.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
    ! Determine geometric mean diameter for mode and sea-salt number
    ! concentration assuming sea-salt is externally mixed
    DO i = 1,nbox
      IF (nd(i,imode) > num_eps(imode)) THEN
        seasalt_mconc(i,idx) = md(i,imode,cp_cl)*nd(i,imode)*                  &
                               (zmwnacl/zmwair)/aird(i)
        tmpdiam = MAX(ddplim0(imode),wetdp(i,imode))
        tmpdiam = MIN(tmpdiam,ddplim1(imode))
        seasalt_diam(i,idx) = tmpdiam*meanfac
        seasalt_nconc(i,idx) = (seasalt_mconc(i,idx) / rhocomp(cp_cl)) *       &
                               sixovrpix(imode) * rhoa(i) / (tmpdiam**3.0)
      ELSE
        seasalt_mconc(i,idx) = 0.0
        seasalt_diam(i,idx) = ddplim0(imode)
        seasalt_nconc(i,idx) = 0.0
      END IF
    END DO
  END IF
END DO !imode

!---------------------------------------------------------------------------
! 2.0 Dust into local arrays (if applicable)
! Initialise arrays to zero
dust_mconc(:,:)=0.0
dust_nconc(:,:)=0.0
dust_diam(:,:)=0.0

IF (do_dust) THEN
  DO imode = mode_acc_insol, mode_cor_insol
    IF (mode(imode)) THEN
      idx=imode-mode_acc_insol+1
      meanfac=EXP(0.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
      ! Determine geometric mean diameter for mode and sea-salt number
      ! concentration assuming sea-salt is externally mixed
      DO i = 1,nbox
        IF (nd(i,imode) > num_eps(imode)) THEN
          dust_mconc(i,idx) = md(i,imode,cp_du)*nd(i,imode)*                     &
                              (zmwcaco3/zmwair)/aird(i)
          tmpdiam = MAX(ddplim0(imode),wetdp(i,imode))
          tmpdiam = MIN(tmpdiam,ddplim1(imode))
          dust_diam(i,idx) = tmpdiam*meanfac
          dust_nconc(i,idx) = (dust_mconc(i,idx) / rhocomp(cp_du)) *             &
                              sixovrpix(imode) * rhoa(i) / (tmpdiam**3.0)
        ELSE
          dust_mconc(i,idx) = 0.0
          dust_diam(i,idx) = ddplim0(imode)
          dust_nconc(i,idx) = 0.0
        END IF
      END DO
    END IF
  END DO !imode
ELSE
  DO imode = mode_acc_insol, mode_cor_insol
    idx=imode-mode_acc_insol+1
    DO i = 1,nbox
      dust_diam(i,idx) = ddplim0(imode)
    END DO
  END DO
END IF ! do_dust

!---------------------------------------------------------------------------
! 3.0 Run internal copy of aer_no3_ukcadu
ptss(:,:) = 0.0
pthno3(:) = 0.0
ptnn(:,:) = 0.0
ptdu(:,:) = 0.0

! Loop over boxes - in ukca_aer_no3_mod, a tropopause limit is applied
! here. This could be replaced with an air density or altitude constraint
DO i =1,nbox

  ! Relative humidity considerations
  iirh=1
  DO jtab=1,12
    IF (RH_clr(i)*100.0 > rrhtab(jtab)) THEN
      iirh=jtab
    END IF
  END DO
  zrhloc = RH_clr(i)
  IF (zrhloc < 0.0 ) zrhloc = 0.0
  IF (zrhloc > 0.98 ) zrhloc = 0.98

  ! needed for HNO3 update on dust and ss
  zdg0= 3.0/(8.0*zavg*rhoa(i)*(zdmolec**2.0))
  zdghno3 = zdg0 * (zrgas*zmgas*t(i)*zratiohno3/(2.0*pi))**0.5

  !RH dependent HNO3 uptake coefficient based on Fairlie et al. 2010
  !and scale for SS
  zrh1580 = MAX(MIN(zrhloc,0.80),0.15)
  IF ( zrh1580 < 0.30 ) THEN
    zghno3_dd = zghno315_dd + (zghno330_dd-zghno315_dd)/0.15*(zrh1580-0.15)
    zghno3_ss = zghno315_ss + (zghno330_ss-zghno315_ss)/0.15*(zrh1580-0.15)
  ELSE IF (zrh1580 >= 0.30 .AND. zrh1580 <= 0.70) THEN
    zghno3_dd = zghno330_dd + (zghno370_dd-zghno330_dd)/0.40*(zrh1580-0.30)
    zghno3_ss = zghno330_ss + (zghno370_ss-zghno330_ss)/0.40*(zrh1580-0.30)
  ELSE IF (zrh1580 > 0.70 ) THEN
    zghno3_dd = zghno370_dd + (zghno380_dd-zghno370_dd)/0.10*(zrh1580-0.70)
    zghno3_ss = zghno370_ss + (zghno380_ss-zghno370_ss)/0.10*(zrh1580-0.70)
  END IF

  ! Uptake of HNO3 on DUST
  IF (do_dust) THEN
    DO jdu = 1,ndumodes
      zkn_dd = 2.0 * mfpa(i) / dust_diam(i,jdu)
      zkhno3_dd(jdu) = (2.0*pi*dust_diam(i,jdu)*zdghno3) /                     &
                       (1.0+(4.0*zkn_dd/zghno3_dd/3.0)*(1.0-0.47*zghno3_dd/    &
                       (1.0+zkn_dd)))
      ddtmp = MAX(0.0,(dust_nconc(i,jdu)*zmwca/zmwhno3*hno3_mconc(i)*dtz))
      IF ((dust_mconc(i,jdu)*pfrac_ca-ddtmp*zkhno3_dd(jdu)<0.0) .AND.          &
          (ddtmp>0.0))                                                         &
        zkhno3_dd(jdu) = dust_mconc(i,jdu)*pfrac_ca/ddtmp
    END DO
  ELSE
    DO jdu = 1,ndumodes
      zkhno3_dd(jdu) = 0.0
    END DO
  END IF

  ! Uptake of HNO3 on SEA-SALT
  DO jss = 1,nssmodes
    zkn_ss = 2.0 * mfpa(i) / (seasalt_diam(i,jss) * rssgrowth_rhtab(iirh))
    zkhno3_ss(jss) = (2.0*pi*seasalt_diam(i,jss)*rssgrowth_rhtab(iirh)*        &
                     zdghno3) /(1.0+(4.0*zkn_ss/zghno3_ss/3.0)*                &
                     (1.0-0.47*zghno3_ss/(1.0+zkn_ss)))
    sstmp = MAX(0.0,(seasalt_nconc(i,jss)*(zmwnacl/zmwhno3)*hno3_mconc(i)*dtz))
    IF ((seasalt_mconc(i,jss)-sstmp*zkhno3_ss(jss)<0.0) .AND. (sstmp>0.0))     &
      zkhno3_ss(jss) = seasalt_mconc(i,jss)/sstmp
  END DO

  ! Divide NO3 tendencies into accum and coarse modes
  zdd1 = MAX(0.0,zkhno3_dd(1)*dust_nconc(i,1))*(zmwcano3_2/zmwhno3)
  zdd2 = MAX(0.0,zkhno3_dd(2)*dust_nconc(i,2))*(zmwcano3_2/zmwhno3)
  zss1 = MAX(0.0,zkhno3_ss(1)*seasalt_nconc(i,1))*(zmwnano3/zmwhno3)
  zss2 = MAX(0.0,zkhno3_ss(2)*seasalt_nconc(i,2))*(zmwnano3/zmwhno3)
  zaccno3 = zdd1 + zss1
  zcorno3 = zdd2 + zss2

  ! If the NO3 tendency is greater than 0 then store in output array
  IF ((zaccno3+zcorno3) > 0.0) THEN
    zfracm1 = zaccno3 / (zaccno3+zcorno3)
    zfracm2 = zcorno3 / (zaccno3+zcorno3)

    ! Update NO3 and HNO3 tendencies
    dtno3=(zaccno3+zcorno3)*MAX(0.0,hno3_mconc(i))
    pthno3(i) = -1.0*MAX(0.0,hno3_mconc(i))*                                   &
                       ((zdd1+zdd2)*zmwhno3*2.0/zmwcano3_2 +                   &
                        (zss1+zss2)*zmwhno3/zmwnano3)

    ! Check for negative HNO3
    zhno3_test = hno3_mconc(i) + pthno3(i)*dtz
    mcf=1.0
    IF (zhno3_test < 0.0) THEN
      dtno3old=dtno3
      zthno3_tmp = pthno3(i)
      pthno3(i) = -1.0*hno3_mconc(i)/dtz
      zthno3_tmp = pthno3(i)-zthno3_tmp
      dtno3 = dtno3-zthno3_tmp*zmwnano3/zmwhno3
      mcf=dtno3/dtno3old
    END IF
    ptnn(i,1) = dtno3 * zfracm1
    ptnn(i,2) = dtno3 * zfracm2

    ! Update dust tendencies
    DO jdu = 1,ndumodes
      ptdu(i,jdu)= -1.0 * zkhno3_dd(jdu) * dust_nconc(i,jdu) *                 &
                          hno3_mconc(i) * zmwcaco3 / zmwhno3 * mcf
    END DO

    ! Update sea-salt tendencies
    DO jss = 1,nssmodes
      ptss(i,jss)= -1.0 * zkhno3_ss(jss) * seasalt_nconc(i,jss) *              &
                          hno3_mconc(i) * zmwnacl / zmwhno3 * mcf
    END DO
  END IF  ! Rates are greater than zero
END DO !jk

!---------------------------------------------------------------------------
! 4. Convert coarse-NO3 and sea-salt MMR tendencies to md tendencies
!    and then update md, and mdt
dmd_nn(:,:) = 0.0
dmd_cl(:,:) = 0.0

DO imode = mode_acc_sol, mode_cor_sol
  idx = imode - mode_acc_sol + 1
  DO i = 1,nbox
    IF (nd(i,imode) > num_eps(imode)) THEN
      ! Convert mass tendencies to md tendencies
      dmd_nn(i,imode) = ptnn(i,idx)*(zmwair/zmwnano3)*aird(i)/nd(i,imode)
      dmd_cl(i,imode) = ptss(i,idx)*(zmwair/zmwnacl)*aird(i)/nd(i,imode)

      ! At present, just update sea-salt and not dust which could probably
      ! be handled by the ageing. Just budget for coarse-NO3 produced
      IF ((imode == mode_acc_sol) .AND. (nmascondnnaccsol > 0)) THEN
        bud_aer_mas(i,nmascondnnaccsol) = bud_aer_mas(i,nmascondnnaccsol) +    &
                                          dmd_nn(i,imode)*nd(i,imode)*dtz
      END IF
      IF ((imode == mode_cor_sol) .AND. (nmascondnncorsol > 0)) THEN
        bud_aer_mas(i,nmascondnncorsol) = bud_aer_mas(i,nmascondnncorsol) +    &
                                          dmd_nn(i,imode)*nd(i,imode)*dtz
      END IF

      ! Update md and mdt - note that cp_cl is modified
      md(i,imode,cp_nn) = md(i,imode,cp_nn)+dmd_nn(i,imode)*dtz
      md(i,imode,cp_cl) = md(i,imode,cp_cl)+dmd_cl(i,imode)*dtz
      mdt(i,imode) = 0.0
      DO icp = 1,ncp
        IF (component(imode,icp)) THEN
          mdt(i,imode) = mdt(i,imode) + md(i,imode,icp)
        END IF
      END DO
    ELSE
      DO icp = 1,ncp
        IF (component(imode,icp)) THEN
          md(i,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
        END IF
      END DO
      mdt(i,imode) = mmid(imode)
      nd(i,imode) = 0.0
    END IF ! nd > num_eps

  END DO !  i = 1,nbox
END DO ! imode

! Update gas concentrations
DO i = 1,nbox
  s0g(i,mhno3) = s0g(i,mhno3) + pthno3(i)*sm(i)*(zmwair/zmwhno3)*dtz
  IF (s0g(i,mhno3) < 0.0) s0g(i,mhno3) = 0.0
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_coarse_no3
END MODULE ukca_coarse_no3_mod
