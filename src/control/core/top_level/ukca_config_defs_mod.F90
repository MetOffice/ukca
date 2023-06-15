! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Subroutine to define emitted species lists and various
!   configuration-related array sizes used within UKCA.
!
! Method:
!   Set emissions and array sizes according to which chemistry scheme is
!   selected.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_config_defs_mod

IMPLICIT NONE

PUBLIC

INTEGER, SAVE :: n_use_tracers        ! No. of tracers used
INTEGER, SAVE :: n_use_emissions      ! No. of emissions used
INTEGER, SAVE :: n_3d_emissions       ! No. of 3-D emissions used
INTEGER, PARAMETER :: n_boundary_vals  = 7    ! No. lower boundary vals for
                                              ! stratospheric chemistry
INTEGER, PARAMETER :: n_cfc_lumped  = 13 ! No. of CFC species that are not part
                                         ! of boundary values but are lumped
                                         ! into some of the boundary val CFCs
INTEGER, SAVE   :: n_chem_tracers        ! No. tracers for chemistry
INTEGER, SAVE   :: n_aero_tracers        ! No. tracers for aerosol chemistry
INTEGER, SAVE   :: n_mode_tracers        ! No. tracers for MODE
INTEGER, SAVE   :: n_nonchem_tracers     ! No. tracers non-chemistry
INTEGER, SAVE   :: n_chem_emissions      ! No. emissions for chemistry
INTEGER, SAVE   :: nmax_strat_fluxdiags  ! Max no strat flux diags
INTEGER, SAVE   :: nmax_mode_diags       ! Max no MODE diags
INTEGER, SAVE   :: nr_therm              ! No. thermal reactions
INTEGER, SAVE   :: nr_phot               ! No. photolytic reactions

! Names for tracers which have surface emissions
CHARACTER(LEN=10), ALLOCATABLE, SAVE :: em_chem_spec(:)
! Lower BCs for stratospheric species
CHARACTER(LEN=10), SAVE :: lbc_spec(n_boundary_vals)  ! species names
CHARACTER(LEN=10), SAVE :: cfc_lumped(n_cfc_lumped)  ! CFC species lumped
REAL, SAVE :: lbc_mmr(n_boundary_vals)               ! mixing ratios

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CONFIG_DEFS_MOD'

CONTAINS

! ----------------------------------------------------------------------

SUBROUTINE ukca_set_config_defs()

USE ukca_config_specification_mod, ONLY: ukca_config, glomap_config
USE asad_mod,                ONLY: jpro2
USE parkind1,                ONLY: jprb, jpim
USE yomhook,                 ONLY: lhook, dr_hook
USE ereport_mod,             ONLY: ereport
USE missing_data_mod,        ONLY: rmdi
USE errormessagelength_mod,  ONLY: errormessagelength

IMPLICIT NONE

! Local variables

INTEGER, PARAMETER :: ichem_ver132 = 132
INTEGER :: n_mode_emissions       ! No. of emissions for MODE
INTEGER :: errcode                ! Variable passed to ereport
CHARACTER (LEN=errormessagelength) :: cmessage        ! Error message

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_CONFIG_DEFS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Set max fluxdiags to zero initially
nmax_strat_fluxdiags = 0
nmax_mode_diags      = 0
n_nonchem_tracers    = 0

! Initialise LBC arrays used in emiss_ctl
lbc_mmr(:) = rmdi
!  Species with potential lower boundary conditions (same for all flavours)
lbc_spec =                                                                     &
     ['N2O       ','CF2Cl2    ',                                               &
       'CFCl3     ','MeBr      ',                                              &
       'H2        ','CH4       ',                                              &
       'COS       ']
! Species that are lumped into the major CFCs used as boundary values
cfc_lumped = ['CF2ClCFCl2','CF2ClCF2Cl','CF2ClCF3  ','CCl4      ',             &
               'MeCCl3    ','CHF2Cl    ','MeCFCl2   ','MeCF2Cl   ',            &
               'CF2ClBr   ','CF2Br2    ','CF3Br     ','CF2BrCF2Br',            &
               'MeCl      ']

IF (ukca_config%l_ukca_trop) THEN

  ! Standard tropospheric chemistry for B-E solver
  ! ==============================================
  n_chem_emissions = 8
  n_3d_emissions = 1       ! aircraft NOX
  n_aero_tracers = 0
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  n_chem_tracers =  26
  nr_therm       = 102        ! thermal reactions
  nr_phot        = 27         ! photolytic ---"---
  nmax_strat_fluxdiags = n_chem_tracers
  em_chem_spec =                                                               &
  ['NO        ','CH4       ','CO        ','HCHO      ',                        &
    'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',                       &
    'NO_aircrft']
ELSE IF (ukca_config%l_ukca_tropisop .AND. ukca_config%l_ukca_achem) THEN

  ! Std tropospheric chemistry + MIM with aerosol scheme (N-R)
  ! ==========================================================
  n_chem_emissions = 19    ! 2D emission fields
  n_3d_emissions = 4       ! SO2_nat, BC & OC biomass, aircraft NOX
  n_aero_tracers = 9       ! DMS, SO2... aerosol precursor species
  n_chem_tracers = 51
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  nr_therm       = 113
  nr_phot        = 37
  nmax_strat_fluxdiags = n_chem_tracers
  ! Table refers to emissions, more species are emitted using surrogates
  em_chem_spec =                                                               &
  ['NO        ','CH4       ','CO        ','HCHO      ',                        &
    'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',                       &
    'C5H8      ','BC_fossil ','BC_biofuel','OM_fossil ',                       &
    'OM_biofuel','Monoterp  ','MeOH      ','SO2_low   ',                       &
    'SO2_high  ','NH3       ','DMS       ','SO2_nat   ',                       &
    'BC_biomass','OM_biomass','NO_aircrft']
ELSE IF (ukca_config%l_ukca_tropisop .AND. .NOT. ukca_config%l_ukca_achem) THEN

  ! Std tropospheric chemistry with MIM isoprene scheme (N-R)
  ! =========================================================
  n_chem_emissions = 9
  n_3d_emissions = 1       ! aircraft NOX
  n_aero_tracers = 0
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  n_chem_tracers = 49
  nr_therm       = 132
  nr_phot        = 35
  nmax_strat_fluxdiags = n_chem_tracers
  em_chem_spec =                                                               &
  ['NO        ','CH4       ','CO        ','HCHO      ',                        &
    'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',                       &
    'C5H8      ','NO_aircrft']
ELSE IF (ukca_config%l_ukca_aerchem) THEN

  ! Std trop chem with SO2, DMS, NH3, and monoterpene (BE)
  ! ======================================================
  n_chem_emissions = 18       ! Surface/ high-level emissions
  n_3d_emissions = 4          ! SO2_nat, aircraft NOX, OC & BC Biomass
  n_chem_tracers = 26         ! advected chemical tracers
  n_aero_tracers =  7         ! advected aerochem ---"---
  nr_therm       = 137        ! thermal reactions
  nr_phot        = 27         ! photolytic ---"---
  nmax_strat_fluxdiags = n_chem_tracers
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  em_chem_spec =                                                               &
  ['NO        ','CH4       ','CO        ','HCHO      ',                        &
    'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',                       &
    'C5H8      ','BC_fossil ','BC_biofuel','OM_fossil ',                       &
    'Monoterp  ','MeOH      ','SO2_low   ','SO2_high  ',                       &
    'NH3       ','DMS       ','SO2_nat   ','BC_biomass',                       &
    'OM_biomass','NO_aircrft']
ELSE IF (ukca_config%l_ukca_raq) THEN

  ! Regional air quality chemistry (RAQ), based on STOCHEM
  ! ========================================================
  n_chem_emissions  = 16
  n_3d_emissions    = 1       ! aircraft NOx
  n_chem_tracers    = 40      ! advected chemical tracers
  n_aero_tracers    = 0       ! advected aerochem tracers
  nr_therm          = 192     ! thermal reactions
  nr_phot           = 23      ! photolytic reacs
  nmax_strat_fluxdiags = n_chem_tracers
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  em_chem_spec =                                                               &
    ['NO        ','CH4       ','CO        ','HCHO      ',                      &
      'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',                     &
      'C5H8      ','C4H10     ','C2H4      ','C3H6      ',                     &
      'TOLUENE   ','oXYLENE   ','CH3OH     ','H2        ',                     &
      'NO_aircrft' ]
ELSE IF (ukca_config%l_ukca_raqaero) THEN

  ! Regional air quality chemistry plus aerosols RAQ-AERO
  ! ========================================================
  n_chem_emissions  = 20
  n_3d_emissions    = 1       ! aircraft NOx
  n_chem_tracers    = 40      ! advected chemical tracers
  n_aero_tracers    = 8       ! advected aerochem tracers
  nr_therm          = 197     ! thermal reactions
  nr_phot           = 23      ! photolytic reacs
  nmax_strat_fluxdiags = n_chem_tracers

  IF (ukca_config%l_ukca_mode) THEN
    ! BC & OC fossil, biofuel, & biomass emissions
    n_chem_emissions = n_chem_emissions + 6
  END IF

  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  em_chem_spec(1:21) =                                                         &
    ['NO        ','CH4       ','CO        ','HCHO      ',        & ! 4
      'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',        & ! 8
      'C5H8      ','C4H10     ','C2H4      ','C3H6      ',        & ! 12
      'TOLUENE   ','oXYLENE   ','CH3OH     ','Monoterp  ',        & ! 16
      'SO2_low   ','SO2_high  ','DMS       ','NH3       ',        & ! 20
      'NO_aircrft'                                                & ! 21
      ]

  IF (ukca_config%l_ukca_mode) THEN
    ! The fossil/biofuel pairs will be 2D;
    ! the biomass emissions will also be 2D
    em_chem_spec(22:23) = [ 'BC_fossil ', 'BC_biofuel' ]
    em_chem_spec(24:25) = [ 'OM_fossil ', 'OM_biofuel' ]
    em_chem_spec(26:27) = [ 'OM_biomass', 'BC_biomass' ]
  END IF


ELSE IF (ukca_config%l_ukca_offline_be) THEN

  ! Offline oxidants scheme with aerosol chemistry
  ! ==============================================
  n_chem_emissions = 8     ! 2D emission fields
  n_3d_emissions = 3       ! SO2_nat, BC & OC biomass
  n_aero_tracers = 7       ! DMS, SO2... aerosol precursor species
  n_chem_tracers = 0
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  nr_therm       = 11      ! ratb + ratt
  nr_phot        = 0
  nmax_strat_fluxdiags = n_chem_tracers

  ! Table refers to emissions, more species may be emitted using surrogates
  em_chem_spec =                                                               &
  ['BC_fossil ','BC_biofuel','OM_fossil ','OM_biofuel',                        &
    'Monoterp  ','SO2_low   ','SO2_high  ','DMS       ',                       &
    'SO2_nat   ','BC_biomass','OM_biomass']

ELSE IF (ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_strattrop .OR.       &
         ukca_config%l_ukca_stratcfc) THEN

  ! Stratospheric chemistry
  ! =======================
  n_nonchem_tracers  = 1            ! Passive O3 included by default

  IF (ukca_config%l_ukca_strat) THEN
    IF (.NOT. ukca_config%l_ukca_achem) THEN ! NOT using aerosol chemistry
       ! emissions:
      n_chem_emissions = 4
      n_3d_emissions = 1       ! aircraft NOX
      n_aero_tracers =  0
      ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
      em_chem_spec =                                                           &
           ['NO        ','CH4       ',                                         &
             'CO        ','HCHO      ',                                        &
             'NO_aircrft']
      ! tracers and reactions:
      n_chem_tracers  = 37   ! CCMVal !!No H2OS, but does have H2O
      nr_therm       = 135
      nr_phot        = 34
    ELSE ! USING AEROSOL CHEMISTRY
       ! emissions:
      n_chem_emissions = 7       ! em_chem_spec below
      n_3d_emissions   = 2       ! volc SO2 & aircraft NOX
      ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
      em_chem_spec =                                                           &
           ['NO        ','CH4       ',                                         &
             'CO        ','HCHO      ',                                        &
             'SO2_low   ','SO2_high  ',                                        &
             'DMS       ','SO2_nat   ',                                        &
             'NO_aircrft']
      ! tracers and reactions:
      n_chem_tracers  = 45
      n_aero_tracers  = 8
      nr_therm        = 149
      nr_phot         = 38
    END IF

    ! Strat-trop chemistry
    ! =======================
  ELSE IF (ukca_config%l_ukca_strattrop) THEN
    IF (.NOT. ukca_config%l_ukca_achem) THEN  ! If NOT using aerosol chemistry
      n_chem_emissions = 10
      n_3d_emissions = 1       ! aircraft NOX
      n_aero_tracers =  0
      ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
      em_chem_spec =                                                           &
          ['NO        ','CH4       ','CO        ','HCHO      ',                &
            'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',               &
            'C5H8      ','MeOH      ','NO_aircrft']
      n_chem_tracers = 71         ! No chem tracers
      nr_therm       = 220        ! thermal reactions
      nr_phot        = 55         ! photolytic (ATA)

    ELSE  ! If using aerosol chemistry
      n_chem_emissions = 19      ! em_chem_spec below
      n_3d_emissions   = 4       ! BC, OC, volc SO2 & aircraft NOX
      ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
      em_chem_spec =                                                           &
          ['NO        ','CH4       ','CO        ','HCHO      ',                &
            'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',               &
            'C5H8      ','BC_fossil ','BC_biofuel','OM_fossil ',               &
            'OM_biofuel','Monoterp  ','MeOH      ','SO2_low   ',               &
            'SO2_high  ','NH3       ','DMS       ','SO2_nat   ',               &
            'BC_biomass','OM_biomass','NO_aircrft']
      IF (ukca_config%i_ukca_chem_version >= ichem_ver132) THEN
        ! Include secondary organic species from isoprene oxidation (SEC_ORG_I)
        n_aero_tracers = 13
      ELSE
        n_aero_tracers = 12
      END IF
      n_chem_tracers = 71         ! No chem tracers
      IF (ukca_config%l_ukca_trophet) THEN
        nr_therm     = 241        ! thermal reactions
      ELSE
        nr_therm     = 239        ! thermal reactions
      END IF
      nr_phot        = 59         ! photolytic (ATA)

    END IF

    ! Make correction to number of tracers if turning off RO2 transport:
    ! reduce No. of transported tracers by No. of peroxy-radical species
    IF (ukca_config%l_ukca_ro2_ntp) THEN
      n_chem_tracers = n_chem_tracers - jpro2
    END IF

  ELSE IF (ukca_config%l_ukca_stratcfc) THEN
    n_chem_tracers = 43
  END IF
  nr_therm       = 102        ! thermal reactions
  nr_phot        = 27         ! photolytic (ATA)
  nmax_strat_fluxdiags = n_chem_tracers
ELSE IF (ukca_config%l_ukca_offline) THEN

  ! Offline oxidants scheme with aerosol chemistry
  ! ==============================================
  n_chem_emissions = 8     ! 2D emission fields
  n_3d_emissions = 3       ! SO2_nat, BC & OC biomass
  n_aero_tracers = 7       ! DMS, SO2... aerosol precursor species
  n_chem_tracers = 0
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  nr_therm       = 11      ! ratb + ratt, unused ?
  nr_phot        = 0
  nmax_strat_fluxdiags = n_chem_tracers

  ! Table refers to emissions, more species may be emitted using surrogates
  em_chem_spec =                                                               &
  ['BC_fossil ','BC_biofuel','OM_fossil ','OM_biofuel',                        &
    'Monoterp  ','SO2_low   ','SO2_high  ','DMS       ',                       &
    'SO2_nat   ','BC_biomass','OM_biomass']



  ! CRI-Strat and CRI-Strat 2 chemistry
  ! =======================
ELSE IF (ukca_config%l_ukca_cristrat) THEN
  IF (.NOT. ukca_config%l_ukca_achem) THEN ! Without aerosol
    IF (ukca_config%i_ukca_chem_version >= 119) THEN  ! CRI-Strat 2 w/o aerosol
      n_chem_tracers  = 171   ! advected chemical tracers (+9, -4)
      n_aero_tracers  = 0    ! DMS, SO2... aerosol precursor species
      nr_therm        = 625  ! thermal reactions (+45 bimol, +8 termol)
      nr_phot         = 135  ! photolytic reacs
    ELSE  ! CRI-Strat without aerosol
      n_chem_tracers  = 166    ! advected chemical tracers
      n_aero_tracers  = 0    ! DMS, SO2... aerosol precursor species
      nr_therm        = 572  ! thermal reactions
      nr_phot         = 121  ! photolytic reacs
    END IF
    n_chem_emissions  = 26
    n_3d_emissions    = 1       ! aircraft NOX
    ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
    em_chem_spec =                                                             &
      ['NO        ','CH4       ','CO        ','C2H6      ',                    &
        'C3H8      ','C4H10     ','C2H4      ','C3H6      ',                   &
        'TBUT2ENE  ','C2H2      ','C5H8      ','APINENE   ',                   &
        'BPINENE   ','BENZENE   ','TOLUENE   ','oXYLENE   ',                   &
        'MeOH      ','EtOH      ','HCHO      ','MeCHO     ',                   &
        'EtCHO     ','Me2CO     ','MEK       ','HCOOH     ',                   &
        'MeCO2H    ','HOCH2CHO  ','NO_aircrft']
  ELSE ! With aerosol
    IF (ukca_config%i_ukca_chem_version >= 119) THEN  ! CRI-Strat 2 with aerosol
      n_chem_tracers  = 171   ! advected chemical tracers (+9, -4)
      n_aero_tracers  = 17    ! DMS, SO2... aerosol precursor species
      nr_phot         = 137     ! photolytic reacs
      IF (ukca_config%l_ukca_trophet) THEN
        nr_therm      = 644        ! thermal reactions  (+45 bimol, +8 termol)
      ELSE
        nr_therm      = 646        ! thermal reactions  (+45 bimol, +8 termol)
      END IF
    ELSE ! CRI-Strat with aerosol
      n_chem_tracers  = 166    ! advected chemical tracers
      n_aero_tracers  = 17    ! DMS, SO2... aerosol precursor species
      nr_phot         = 123     ! photolytic reacs
      IF (ukca_config%l_ukca_trophet) THEN
        nr_therm      = 591        ! thermal reactions
      ELSE
        nr_therm      = 593        ! thermal reactions
      END IF
    END IF
    n_chem_emissions = 34      ! em_chem_spec below
    n_3d_emissions   = 4       ! BC, OC, volc SO2 & aircraft NOX
    ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
    em_chem_spec =                                                             &
        ['NO        ','CH4       ','CO        ','C2H6      ',                  &
          'C3H8      ','C4H10     ','C2H4      ','C3H6      ',                 &
          'TBUT2ENE  ','C2H2      ','C5H8      ','APINENE   ',                 &
          'BPINENE   ','BENZENE   ','TOLUENE   ','oXYLENE   ',                 &
          'MeOH      ','EtOH      ','HCHO      ','MeCHO     ',                 &
          'EtCHO     ','Me2CO     ','MEK       ','HCOOH     ',                 &
          'MeCO2H    ','HOCH2CHO  ','BC_fossil ','BC_biofuel',                 &
          'OM_fossil ','OM_biofuel','SO2_low   ','SO2_high  ',                 &
          'NH3       ','DMS       ','SO2_nat   ','BC_biomass',                 &
          'OM_biomass','NO_aircrft']
  END IF
ELSE

  ! No chemistry
  ! ============
  n_chem_emissions = 0
  n_chem_tracers   = 0
  ALLOCATE(em_chem_spec(1))
  em_chem_spec = '          '
END IF

IF (ukca_config%l_ukca_ageair) THEN
  ! Include Age of air tracer
  n_nonchem_tracers = n_nonchem_tracers + 1
END IF

IF (ukca_config%l_ukca_mode) THEN
  IF (.NOT. (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_nr_aqchem .OR. &
             ukca_config%l_ukca_offline_be .OR.                                &
             ukca_config%l_ukca_raqaero)) THEN
    cmessage=' the driving aerosol chemistry is needed for GLOMAP-MODE'
    errcode=1
    CALL ereport('UKCA_SET_CONFIG_DEFS',errcode,cmessage)
  END IF
  n_mode_emissions = 0          ! See aerosol chemistry section
  ! n_mode_tracers is now calculated in UKCA_INIT
  nmax_mode_diags  = 359
ELSE
  n_mode_emissions = 0
END IF

IF (glomap_config%l_ukca_primdu) THEN
  glomap_config%n_dust_emissions = 6  ! Use 6 bin dust emissions (Woodward,2001)
ELSE
  glomap_config%n_dust_emissions = 0
END IF

n_use_tracers   = n_chem_tracers + n_mode_tracers +                            &
                  n_aero_tracers + n_nonchem_tracers
n_use_emissions = n_chem_emissions + n_mode_emissions +                        &
                  n_3d_emissions

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_set_config_defs

! ----------------------------------------------------------------------

END MODULE ukca_config_defs_mod
