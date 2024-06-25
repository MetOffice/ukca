! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE inijtab_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE umPrintMgr, ONLY:                                                          &
    umPrint,                                                                   &
    umMessage

IMPLICIT NONE



! Description:
!     Subroutine to initialise photolysis tables
!      AJ2A:     O2       +  hv             ->   O     +   O
!      AJ2B:     O2       +  hv             ->   O     +   O(1D)
!      AJ3:      O3       +  hv             ->   O2    +   O(3P)
!      AJ3A:     O3       +  hv(<310nm)     ->   O2    +   O(1D)
!      AJNO:     NO       +  hv             ->   N     +   O
!      AJNO2:    NO2      +  hv             ->   NO    +   O(3P)
!      AJNO31:   NO3      +  hv             ->   NO    +   O2
!      AJNO32:   NO3      +  hv             ->   NO2   +   O
!      AJHNO3:   HNO3     +  hv             ->   NO2   +   OH
!      AJPNA:    HO2NO2   +  hv             ->   NO2   +   HO2 / NO3 + OH
!      AJN2O:    N2O      +  hv             ->   N2    +   O(1D)
!      AJN2O5:   N2O5     +  hv             ->   NO2   +   NO3
!      AJH2O:    H2O      +  hv             ->   OH    +   H
!      AJH2O2:   H2O2     +  hv             ->   OH    +   OH
!      AJCNITA:  ClONO2   +  hv             ->   Cl    +   NO3
!      AJCNITB:  ClONO2   +  hv             ->   ClO   +   NO2
!      AJF11:    CFCl3    +  hv             ->   3Cl
!      AJF12:    CF2Cl2   +  hv             ->   2Cl
!      AJF22:    CHF2Cl   +  hv             ->   Cl
!      AJF113:   CF2CFCl3 +  hv             ->   3Cl
!      AJCH3CL:  CH3Cl    +  hv             ->   CH3   +   Cl
!      AJCCL4:   CCl4     +  hv             ->   4Cl
!      AJHCL:    HCl      +  hv             ->   H     +   Cl
!      AJHOCL:   HOCl     +  hv             ->   OH    +   Cl
!      AJOCLO:   OClO     +  hv             ->   O     +   ClO
!      AJCL2O2:  Cl2O2    +  hv             ->   O2    +   2Cl
!      AJBRO:    BrO      +  hv             ->   Br    +   O
!      AJHOBR:   HOBr     +  hv             ->   Br    +   OH
!      AJBRNO3:  BrNO3    +  hv             ->   Br    +   NO3 / Bro + NO2
!      AJBRCL:   BrCl     +  hv             ->   Br    +   Cl
!      AJC2OA:   CH2O     +  hv             ->   H     +   CHO
!      AJC2OB:   CH2O     +  hv             ->   H2    +   CO
!      AJMHP :   CH3OOH   +  hv             ->   CH3O  +   OH
!      AJCH3BR:  CH3Br    +  hv             ->   CH3   +   Br
!      AJMCFM:   CH3CCl3  +  hv             ->   CHx   +   3Cl
!      AJCH4:    CH4      +  hv             ->   CH3   +   H
!      AJF12B1:  CBrClF2  +  hv             ->   Cxx   +   Cl   + Br
!      AJF13B1:  CBrF3    +  hv             ->   Cxx   +   Br
!      AJCOF2:   COF2     +  hv             ->   Cx    +   2HF
!      AJCOFCL:  COFCl    +  hv             ->   Cx    +   Cl + HF
!      AJCO2:    CO2      +  hv             ->   CO    +   O(3P)
!      AJCOS:    COS      +  hv             ->   CO    +   S
!      AJHONO:   HONO     +  hv             ->   OH    +   NO
!      AJMENA:   MeONO2   +  hv             ->   HO2   +   HCHO + NO2
!      AJCHBr3:  CHBr3    +  hv             ->   HBr   +   2Br
!      AJDBRM :  CH2Br2   +  hv             ->   2Br   +   H2O
!      AJCS2:    CS2      +  hv             ->   COS   +   SO2
!      AJH2SO4:  H2SO4    +  hv             ->   SO3   +   OH
!      AJSO3:    SO3      +  hv             ->   SO2   +   O(3P)


!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIJTAB_MOD'

CONTAINS
SUBROUTINE inijtab(mype, current_time, l_ukca_fastjx)
USE ukca_stdto3_mod, ONLY: stdto3, nin
USE ukca_tbjs_mod, ONLY:     tabj2a   , tabj2b   , tabj3    ,                  &
       tabj3a   , tabjno   , tabjno31 , tabjno32 , tabjno2  ,                  &
       tabjn2o5 , tabjhno3 , tabjh2o  , tabjh2o2 , tabjf11  ,                  &
       tabjf12  , tabjf22  , tabjf113 , tabjch3cl, tabjccl4 ,                  &
       tabjcnita, tabjcnitb, tabjhcl  , tabjhocl , tabjpna  ,                  &
       tabjcl2o2, tabjn2o  , tabjbrcl , tabjbrno3, tabjhobr ,                  &
       tabjbro  , tabjoclo , tabjc2oa , tabjc2ob , tabjmhp  ,                  &
       tabjmcfm , tabjch3br, tabjf12b1, tabjf13b1, tabjcof2 ,                  &
       tabjcofcl, tabjch4  , tabjcos  , tabjco2  ,                             &
       tabjhono , tabjmena , tabjchbr3, tabjdbrm , tabjcs2  ,                  &
       tabjh2so4, tabjso3  , tabpres  , tabang,    tabt     ,                  &
       tabo3,     sao3c
USE ukca_crossec_mod, ONLY: accl4, acl2o2, amcfm, acnita,                      &
    acnitb, aco2 , af11 , af113, af12 , af22 , ah2o , ah2o2,                   &
    ahcl , ahno3, ahocl, ach3cl, an2o , an2o5, ano  ,                          &
    ano2 , ano31, ano32, ao2  , ao2a , ao2b , abrcl, abrno3,                   &
    abro , ahobr, aoclo, ac2oa, ac2ob, amhp , ao2sr, ao3  ,                    &
    apna , acof2, acofcl, ach3br, af12b1, af13b1, scs  , qeno2,                &
    qeo1d, quanta, ach4 , wavenm, ahono, aocs ,                                &
    amena, achbr3, adbrm, acs2 , ah2so4, aso3
USE ukca_parpho_mod, ONLY: jplevp1, jplev, jpchi,                              &
                           jpwav, jplo, jphi, jptem, jpo3p

! Module procedures used by this subroutine
USE fill_spectra_mod, ONLY:fill_spectra
USE scatcs_mod, ONLY: scatcs
USE setzen_mod, ONLY: setzen
USE settab_mod, ONLY: settab
USE lymana_mod, ONLY: lymana
USE acsn2o5_mod, ONLY: acsn2o5
USE acsn2o_mod, ONLY: acsn2o
USE acsccl4_mod, ONLY: acsccl4
USE acsf11_mod, ONLY: acsf11
USE acsf12_mod, ONLY: acsf12
USE acsf22_mod, ONLY: acsf22
USE acsmc_mod, ONLY: acsmc
USE quanto12_mod, ONLY: quanto12
USE acso3_mod, ONLY: acso3
USE acshno3_mod, ONLY: acshno3
USE acscnit_mod, ONLY: acscnit
USE acsh2o2_mod, ONLY: acsh2o2
USE acsbrcl_mod, ONLY: acsbrcl
USE acsno2_mod, ONLY: acsno2
USE acscos_mod, ONLY: acscos
USE acsmena_mod, ONLY: acsmena
USE acsdbrm_mod, ONLY: acsdbrm
USE acscs2_mod, ONLY: acscs2
USE acsh2so4_mod, ONLY: acsh2so4
USE acsso3_mod, ONLY: acsso3
USE acsno_mod, ONLY: acsno
USE acssr_mod, ONLY: acssr

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE fastjx_data, ONLY: solcyc_quanta
USE photol_solflux_mod, ONLY: photol_solflux
USE photol_config_specification_mod, ONLY: photol_config
USE photol_constants_mod,  ONLY:  gg => const_g, rgas => const_r

USE missing_data_mod,      ONLY: rmdi

#if !defined(LFRIC)
USE mpl, ONLY: mpl_real
USE ukca_option_mod, ONLY: ukcasto3_file => ukcasto3
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport
#endif

IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: mype  ! processor number
INTEGER, INTENT(IN) :: current_time(7) ! model time for photol_solflux
LOGICAL, INTENT(IN) :: l_ukca_fastjx

! Local variables
! do not read STDTO3 input file. Instead, use values from the module.
! READ_STDO3 controls whether to read temp & o3 data from file. However,
! this mode is not currently supported (filename no longer in RUN_UKCA
! namelist) and should be removed in the near future.
LOGICAL, PARAMETER :: read_stdto3 = .FALSE.

! used for solar cycle to inform calulation method
LOGICAL, PARAMETER :: l_lookup = .TRUE.

REAL, PARAMETER :: colfac=2.132e20
!     O2 mixing ratio
REAL, PARAMETER :: r2o2=0.209
!     Ground albedo
REAL, PARAMETER :: albedo=0.30

!     Whether scattering is on or not.
LOGICAL, PARAMETER :: lscat =.TRUE.

INTEGER :: jc
INTEGER :: jl
INTEGER :: jlu
INTEGER :: jll
INTEGER :: jo
INTEGER :: jt
INTEGER :: jw
INTEGER :: l

REAL :: o3above

REAL :: cinf
REAL :: csup
REAL :: dpdif
REAL :: pdif
REAL :: qu

REAL :: factor(jpwav)
REAL :: do2c  (jplev)    ! Average number density of O2 within
!                          each level in molecules per cm^3.
REAL :: do3c  (jplev)    ! Average number density of O3 within
!                          each level in molecules per cm^3.
REAL :: drsc  (jplev)    ! [M] within level
REAL :: do3cu (jplev)    ! Scaled O3 profiles
REAL :: pres  (jplevp1)  ! Pressure (mb) at the edge of each level
REAL :: presc (jplev)    ! Pressure (mb) at the centre of each level
REAL :: o3vmrc(jplev)
REAL :: sao2  (jplevp1)  ! Vertical O2 column above the edge of
!                          each level in molecules per cm^2

REAL :: tabs  (jplev,jpchi,jpwav)
REAL :: tempc (jplev)            ! Temperature (K) at the centre of each level
REAL :: tspo2 (jplev,jpchi)
REAL :: zenmax(jplev)
REAL :: alt   (jplevp1)          ! Altitude at the edge of each level in km
REAL :: altc  (jplev)            ! Altitude at the centre of each level in km
REAL :: dalt  (jplev)            ! Thickness of level (km)

REAL :: rdtemp(nin)
REAL :: rdpres(nin)
REAL :: rdo3v (nin)

INTEGER :: last_wn

REAL :: tau                ! dummy variable for solar cycle routine
REAL :: quanta_cyc (jpwav) ! solar cyclical component of quanta
REAL :: quanta_in (jpwav)

#if !defined(LFRIC)
INTEGER :: my_comm, icode
INTEGER :: ukcasto3_unit
CHARACTER(LEN=errormessagelength) :: cmessage
#endif

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIJTAB'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if !defined(LFRIC)

CALL gc_get_communicator(my_comm, icode)

IF (read_stdto3) THEN

  ! This mode is not currently supported and the code below needs
  ! to be checked and removed.
  icode = 1
  cmessage='Std temp+O3 file read not supported. Needs a namelist change'
  CALL ereport(RoutineName, icode, cmessage)

  IF ( mype == 0 ) THEN
    ! Read in standard atmosphere temperature and O3
    CALL assign_file_unit(ukcasto3_file, ukcasto3_unit, handler="fortran")
    OPEN(ukcasto3_unit,FILE=ukcasto3_file,ACTION='READ')
    REWIND(ukcasto3_unit)
    DO l=1,nin
      READ(UNIT=ukcasto3_unit,FMT=*) rdpres(l),rdtemp(l),rdo3v(l)
    END DO
    CLOSE(ukcasto3_unit)
    CALL release_file_unit(ukcasto3_unit, handler="fortran")
  END IF !mype==0
  CALL mpl_bcast(rdpres, nin, mpl_real, 0, my_comm, icode)
  CALL mpl_bcast(rdtemp, nin, mpl_real, 0, my_comm, icode)
  CALL mpl_bcast(rdo3v , nin, mpl_real, 0, my_comm, icode)

END IF

#endif

IF (.NOT. read_stdto3) THEN
  ! Read standard atmosphere data from string variable
  DO l=1,nin
    READ(UNIT=stdto3(l),FMT=*) rdpres(l),rdtemp(l),rdo3v(l)
  END DO
END IF

!     Set up j table inter levels equally spaced in log(p)
pres(      1)=rdpres(  1)
pres(jplevp1)=rdpres(nin)
pdif=LOG(rdpres(1)) - LOG(rdpres(nin))
dpdif=pdif/(REAL(jplevp1-1))
DO l=2,jplevp1-1
  pres (l)=EXP(LOG(pres(l-1)) - dpdif)
END DO
!     Pressure at jtable mid levels
DO l=1,jplev
  presc  (l)=EXP(0.5*(LOG(pres(l)) + LOG(pres(l+1))))
  tabpres(l)=presc(l)
END DO

!     Interpolate temperature and O3 to level centres
DO l=1,jplev
  inner_loop: DO jl=1,nin
    IF (rdpres(jl) < presc(l)) EXIT inner_loop
  END DO inner_loop
  jlu=jl
  jll=jl-1
  csup=(LOG(rdpres(jll)) - LOG(presc (l  )))/                                  &
    (LOG(rdpres(jll)) - LOG(rdpres(jlu)))
  cinf=1.0 - csup
  tempc (l)=csup*rdtemp(jlu) + cinf*rdtemp(jll)
  o3vmrc(l)=csup*rdo3v (jlu) + cinf*rdo3v (jll)
END DO

!     Work out geopotential height (km). Bottom level is ground.
pres(1)=1013.0
alt (1)=0.0
DO l=2,jplevp1
  alt(l)=alt(l-1) +0.001*rgas*tempc(l-1)*LOG(pres(l-1)/pres(l))/gg
END DO
!     Height at centre of level and thickness of level (km).
DO l=1,jplev
  altc(l)=0.5*(alt(l+1) + alt(l))
  dalt(l)=     alt(l+1) - alt(l)
END DO

!     Work out O2 column above interfaces
DO l=1,jplevp1
  sao2(l)=r2o2*100.0*pres(l)*colfac
END DO

!     Average [O2] and [M] within a level
DO l=1,jplev
  do2c(l)=(sao2(l)-sao2(l+1))/(1.0e5*dalt(l))
  drsc(l)= do2c(l)/r2o2
END DO
!     Add O2 above top level to top level
do2c(jplev)=do2c(jplev) + sao2(jplevp1)/(1.0e5*dalt(jplev))

!     Average [O3] within a level
DO l=1,jplev
  do3c(l)=o3vmrc(l)*do2c(l)/r2o2
END DO

!     O3 above top level with scale height of 3.5 km
o3above=0.5*o3vmrc(jplev)*sao2(jplevp1)/r2o2
do3c (jplev)=do3c(jplev) + o3above/(1.0e5*dalt(jplev))

!     Column O3 above centre:
sao3c(jplev)=0.5*(o3above + 1.0e5*do3c(jplev)*dalt(jplev))
DO l=jplev-1,1,-1
  sao3c(l)=sao3c(l+1) + 0.5e5*(dalt(l+1)*do3c(l+1)                             &
    +dalt(l  )*do3c(l  ))
END DO
! Call fill_spectra routine here to define spectral information etc.
CALL fill_spectra


! IF FAST-JX is selected, only do first 45 wavelength (up to 177 nm)
IF (l_ukca_fastjx) THEN
  DO jw = 1, jphi
    IF (wavenm(jw) < 177.0) last_wn = jw
  END DO
ELSE
  last_wn = jphi
END IF

WRITE(umMessage,'(1X,A,1X,I3,1X,I3,1X,E12.5)') 'last_wn = ',last_wn, jphi,     &
       wavenm(last_wn)
CALL umPrint(umMessage,src='inijtab_mod')

!     Calculate the Rayleigh scattering cross section.
CALL scatcs(jpwav,scs,wavenm)

!     Set zenith angle, max za, temperature and O3 grid.
CALL setzen(tabang,zenmax,altc,tabt,tabo3)

! If using solar cycle calculate cyclical compontent of quanta
IF (photol_config%i_solcylc_type > 0) THEN
  ! tau should not be used so set to missing data
  tau = rmdi
  CALL photol_solflux(current_time, tau, l_lookup, solcyc_quanta, quanta_cyc)
  quanta_in = quanta*(1.0+quanta_cyc)
ELSE
  quanta_in = quanta
END IF

!     Loop over O3 profiles
DO jo=1,jpo3p

  !     Scale O3 profile
  DO jl=1,jplev
    do3cu(jl)=do3c(jl)*tabo3(jo)
  END DO

  CALL settab(alt,altc,do2c,do3cu,tempc,                                       &
    drsc,tabs,tabang,albedo,lscat,                                             &
    tspo2,scs,ao2,ao2sr,ao3,dalt)

  !  Set up the photolysis rate tabulations. 39 rates.

  tabj2a   (:,:,:,jo)=0.0
  tabj2b   (:,:,:,jo)=0.0
  tabj3    (:,:,:,jo)=0.0
  tabj3a   (:,:,:,jo)=0.0
  tabjno2  (:,:,:,jo)=0.0
  tabjno   (:,:,:,jo)=0.0
  tabjno31 (:,:,:,jo)=0.0
  tabjno32 (:,:,:,jo)=0.0
  tabjn2o  (:,:,:,jo)=0.0
  tabjn2o5 (:,:,:,jo)=0.0
  tabjhno3 (:,:,:,jo)=0.0
  tabjpna  (:,:,:,jo)=0.0
  tabjcnita(:,:,:,jo)=0.0
  tabjcnitb(:,:,:,jo)=0.0
  tabjh2o2 (:,:,:,jo)=0.0
  tabjh2o  (:,:,:,jo)=0.0
  tabjhocl (:,:,:,jo)=0.0
  tabjhcl  (:,:,:,jo)=0.0
  tabjcl2o2(:,:,:,jo)=0.0
  tabjch3cl(:,:,:,jo)=0.0
  tabjccl4 (:,:,:,jo)=0.0
  tabjf11  (:,:,:,jo)=0.0
  tabjf12  (:,:,:,jo)=0.0
  tabjf22  (:,:,:,jo)=0.0
  tabjf113 (:,:,:,jo)=0.0
  tabjbro  (:,:,:,jo)=0.0
  tabjhobr (:,:,:,jo)=0.0
  tabjbrno3(:,:,:,jo)=0.0
  tabjbrcl (:,:,:,jo)=0.0
  tabjoclo (:,:,:,jo)=0.0
  tabjc2oa (:,:,:,jo)=0.0
  tabjc2ob (:,:,:,jo)=0.0
  tabjmhp  (:,:,:,jo)=0.0
  tabjmcfm (:,:,:,jo)=0.0
  tabjch3br(:,:,:,jo)=0.0
  tabjf12b1(:,:,:,jo)=0.0
  tabjf13b1(:,:,:,jo)=0.0
  tabjcof2 (:,:,:,jo)=0.0
  tabjcofcl(:,:,:,jo)=0.0
  tabjch4  (:,:,:,jo)=0.0
  tabjco2  (:,:,:,jo)=0.0
  tabjcos  (:,:,:,jo)=0.0
  tabjcos  (:,:,:,jo)=0.0
  tabjhono (:,:,:,jo)=0.0
  tabjmena (:,:,:,jo)=0.0
  tabjchbr3(:,:,:,jo)=0.0
  tabjdbrm (:,:,:,jo)=0.0
  tabjcs2  (:,:,:,jo)=0.0
  tabjh2so4(:,:,:,jo)=0.0
  tabjso3  (:,:,:,jo)=0.0

  !     Lyman alpha photolysis
  CALL lymana(tspo2,tabj2b,tabjh2o,tabjch4,quanta_in,jo)

  !     Temperature loop.
  DO jt = 1,jptem

    !        Set up the cross-sections for this pressure & temperature.
    !        N2O5
    CALL acsn2o5(tabt(jt),jpwav,wavenm,an2o5)

    !        N2O
    CALL acsn2o (tabt(jt),jpwav,wavenm,an2o)

    !        CCl4
    CALL acsccl4(tabt(jt),jpwav,wavenm,accl4)

    !        F11 (CFCl3)
    CALL acsf11 (tabt(jt),wavenm,af11)

    !        F12 (CF2Cl2)
    CALL acsf12 (tabt(jt),wavenm,af12)

    !        F22 (CHF2Cl)
    CALL acsf22 (tabt(jt),jpwav,wavenm,af22)

    !        CH3Cl
    CALL acsmc  (tabt(jt),jpwav,wavenm,ach3cl)

    !        O[1D] Quantum yield.
    CALL quanto12(tabt(jt),jpwav,wavenm,qeo1d)

    !        O3
    CALL acso3  (tabt(jt),jpwav,ao3)

    !        HNO3
    CALL acshno3(tabt(jt),ahno3)

    !        ClONO2
    CALL acscnit(tabt(jt),wavenm,acnita, acnitb)

    !        H2O2
    CALL acsh2o2(tabt(jt),jpwav,wavenm,ah2o2)

    !        BrCl
    CALL acsbrcl(tabt(jt),jpwav,wavenm,abrcl)

    !        NO2
    CALL acsno2 (tabt(jt),ano2 )

    !        COS
    CALL acscos (tabt(jt),aocs )

    !        MeONO2
    CALL acsmena(tabt(jt), amena )

    !        CHBr3
    ! Do not calculate CHBr3. Not needed for CCMVal
    !         CALL ACSCHBR3(TABT(JT),WAVENM,ACHBR3 )

    !        CH2Br2
    CALL acsdbrm(tabt(jt),adbrm )

    !        CS2
    CALL acscs2 (acs2 )
    !        H2SO4
    CALL acsh2so4 (ah2so4 )
    !        SO3
    CALL acsso3 (aso3 )

    !        Level loop.
    DO jl = 1,jplev

      !        Zenith angle loop.
      DO jc = 1 , jpchi

        !          Calculate the NO pressure dependent cross section.
        CALL acsno(tabang(jc),presc(jl),tspo2(jl,1),ano)

        !          Schumann-Runge bands.
        IF ( tabang(jc) <= zenmax(jl) ) THEN
          CALL acssr(tspo2(jl,jc),jpwav,ao2sr)
        ELSE
          DO jw = jplo, jphi
            ao2sr(jw) = 0.0
          END DO
        END IF

        !          Set up enhancement factor array.
        DO jw = jplo , jphi
          factor(jw) = MAX(0.0,tabs(jl,jc,jw))
        END DO

        !          Wavelength loop.
        DO jw = jplo , last_wn

          !           Number of photons into the volume element.
          !           (need to attenuate quanta above top level?)
          qu = quanta_in(jw) * factor(jw)

          tabj2a (jl,jc,jt,jo)=tabj2a(jl,jc,jt,jo)+qu*(ao2a(jw)+ao2sr(jw))
          tabj2b (jl,jc,jt,jo)=tabj2b(jl,jc,jt,jo)+qu*ao2b(jw)
          tabj3  (jl,jc,jt,jo)=tabj3(jl,jc,jt,jo)+qu*ao3(jw)*(1.0-qeo1d(jw))
          tabj3a (jl,jc,jt,jo)=tabj3a (jl,jc,jt,jo)+qu*ao3 (jw)*qeo1d(jw)
          tabjno2(jl,jc,jt,jo)=tabjno2(jl,jc,jt,jo)+qu*ano2(jw)*qeno2(jw)
          tabjno   (jl,jc,jt,jo)=tabjno   (jl,jc,jt,jo) + qu*ano   (jw)
          tabjno31 (jl,jc,jt,jo)=tabjno31 (jl,jc,jt,jo) + qu*ano31 (jw)
          tabjno32 (jl,jc,jt,jo)=tabjno32 (jl,jc,jt,jo) + qu*ano32 (jw)
          tabjn2o  (jl,jc,jt,jo)=tabjn2o  (jl,jc,jt,jo) + qu*an2o  (jw)
          tabjn2o5 (jl,jc,jt,jo)=tabjn2o5 (jl,jc,jt,jo) + qu*an2o5 (jw)
          tabjhno3 (jl,jc,jt,jo)=tabjhno3 (jl,jc,jt,jo) + qu*ahno3 (jw)
          tabjpna  (jl,jc,jt,jo)=tabjpna  (jl,jc,jt,jo) + qu*apna  (jw)
          tabjcnita(jl,jc,jt,jo)=tabjcnita(jl,jc,jt,jo) + qu*acnita(jw)
          tabjcnitb(jl,jc,jt,jo)=tabjcnitb(jl,jc,jt,jo) + qu*acnitb(jw)
          tabjh2o2 (jl,jc,jt,jo)=tabjh2o2 (jl,jc,jt,jo) + qu*ah2o2 (jw)
          tabjh2o  (jl,jc,jt,jo)=tabjh2o  (jl,jc,jt,jo) + qu*ah2o  (jw)
          tabjhocl (jl,jc,jt,jo)=tabjhocl (jl,jc,jt,jo) + qu*ahocl (jw)
          tabjhcl  (jl,jc,jt,jo)=tabjhcl  (jl,jc,jt,jo) + qu*ahcl  (jw)
          tabjcl2o2(jl,jc,jt,jo)=tabjcl2o2(jl,jc,jt,jo) + qu*acl2o2(jw)
          tabjch3cl(jl,jc,jt,jo)=tabjch3cl(jl,jc,jt,jo) + qu*ach3cl(jw)
          tabjccl4 (jl,jc,jt,jo)=tabjccl4 (jl,jc,jt,jo) + qu*accl4 (jw)
          tabjf11  (jl,jc,jt,jo)=tabjf11  (jl,jc,jt,jo) + qu*af11  (jw)
          tabjf12  (jl,jc,jt,jo)=tabjf12  (jl,jc,jt,jo) + qu*af12  (jw)
          tabjf22  (jl,jc,jt,jo)=tabjf22  (jl,jc,jt,jo) + qu*af22  (jw)
          tabjf113 (jl,jc,jt,jo)=tabjf113 (jl,jc,jt,jo) + qu*af113 (jw)
          tabjbro  (jl,jc,jt,jo)=tabjbro  (jl,jc,jt,jo) + qu*abro  (jw)
          tabjhobr (jl,jc,jt,jo)=tabjhobr (jl,jc,jt,jo) + qu*ahobr (jw)
          tabjbrno3(jl,jc,jt,jo)=tabjbrno3(jl,jc,jt,jo) + qu*abrno3(jw)
          tabjbrcl (jl,jc,jt,jo)=tabjbrcl (jl,jc,jt,jo) + qu*abrcl (jw)
          tabjoclo (jl,jc,jt,jo)=tabjoclo (jl,jc,jt,jo) + qu*aoclo (jw)
          tabjc2oa (jl,jc,jt,jo)=tabjc2oa (jl,jc,jt,jo) + qu*ac2oa (jw)
          tabjc2ob (jl,jc,jt,jo)=tabjc2ob (jl,jc,jt,jo) + qu*ac2ob (jw)
          tabjmhp  (jl,jc,jt,jo)=tabjmhp  (jl,jc,jt,jo) + qu*amhp  (jw)
          tabjch3br(jl,jc,jt,jo)=tabjch3br(jl,jc,jt,jo) + qu*ach3br(jw)
          tabjf12b1(jl,jc,jt,jo)=tabjf12b1(jl,jc,jt,jo) + qu*af12b1(jw)
          tabjf13b1(jl,jc,jt,jo)=tabjf13b1(jl,jc,jt,jo) + qu*af13b1(jw)
          tabjcof2 (jl,jc,jt,jo)=tabjcof2 (jl,jc,jt,jo) + qu*acof2 (jw)
          tabjcofcl(jl,jc,jt,jo)=tabjcofcl(jl,jc,jt,jo) + qu*acofcl(jw)
          tabjmcfm (jl,jc,jt,jo)=tabjmcfm (jl,jc,jt,jo) + qu*amcfm (jw)
          tabjch4  (jl,jc,jt,jo)=tabjch4  (jl,jc,jt,jo) + qu*ach4  (jw)
          tabjco2  (jl,jc,jt,jo)=tabjco2  (jl,jc,jt,jo) + qu*aco2  (jw)
          tabjcos  (jl,jc,jt,jo)=tabjcos  (jl,jc,jt,jo) + qu*aocs  (jw)
          tabjhono (jl,jc,jt,jo)=tabjhono (jl,jc,jt,jo) + qu*ahono (jw)
          tabjmena (jl,jc,jt,jo)=tabjmena (jl,jc,jt,jo) + qu*amena (jw)
          tabjchbr3(jl,jc,jt,jo)=tabjchbr3(jl,jc,jt,jo) + qu*achbr3(jw)
          tabjdbrm (jl,jc,jt,jo)=tabjdbrm (jl,jc,jt,jo) + qu*adbrm (jw)
          tabjcs2  (jl,jc,jt,jo)=tabjcs2  (jl,jc,jt,jo) + qu*acs2  (jw)
          tabjh2so4(jl,jc,jt,jo)=tabjh2so4(jl,jc,jt,jo) + qu*ah2so4(jw)
          tabjso3  (jl,jc,jt,jo)=tabjso3  (jl,jc,jt,jo) + qu*aso3  (jw)

        END DO


        !        End of the zenith angle loop.
      END DO

      !     End of level loop.
    END DO

    !     End of temperature loop.
  END DO

  !     End of O3 profile loop.
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN



END SUBROUTINE inijtab
END MODULE inijtab_mod
