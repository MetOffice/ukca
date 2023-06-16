! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  A routine to set up values such as the number of reactions used
!  by each chemistry scheme
!
! Method:
!  This routine uses UKCA configuration information to set the number of
!  reactants and the number of each type of reaction, plus the number of
!  dry and wet deposited species
!  Also includes error checking for schemes which have not yet
!  been set up
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_setup_chem_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_SETUP_CHEM_MOD'

CONTAINS

SUBROUTINE ukca_setup_chem

USE ukca_config_specification_mod, ONLY:                                       &
   ukca_config, int_method_be_explicit, int_method_nr,                         &
   i_ukca_chem_off, i_ukca_chem_trop, i_ukca_chem_raq,                         &
   i_ukca_chem_tropisop, i_ukca_chem_strattrop, i_ukca_chem_strat,             &
   i_ukca_chem_offline, i_ukca_chem_offline_be, i_ukca_chem_cristrat

USE asad_mod, ONLY: jpctr, jpspec, jpbk, jptk, jppj, jphk, jpdd, jpdw,         &
                    jpcspf, jpnr

USE ukca_config_specification_mod, ONLY: i_liss_merlivat, i_wanninkhof,        &
                              i_nightingale, i_blomquist
USE missing_data_mod, ONLY: imdi
USE ereport_mod,     ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Diag
USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook

USE ukca_chem_offline, ONLY: ndry_offline, nwet_offline
USE ukca_chem_aer,     ONLY: ndry_aer_aer, nwet_aer_aer
USE ukca_chem1_dat,    ONLY: ndry_trop, nwet_trop
USE ukca_chem_raq,     ONLY: ndry_raq, nwet_raq
USE ukca_chem_raqaero_mod, ONLY: ndry_raqaero, nwet_raqaero

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CHARACTER (LEN=errormessagelength) :: cmessage        ! Error message
INTEGER            :: errcode         ! Variable passed to ereport

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SETUP_CHEM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


SELECT CASE (ukca_config%i_ukca_chem)

CASE (i_ukca_chem_off)
  ! Chemistry off completely
  ukca_config%l_ukca_chem     = .FALSE.
  ukca_config%ukca_int_method = 0
  jpctr           = 0
  jpspec          = 0
  jpbk            = 0
  jptk            = 0
  jppj            = 0
  jphk            = 0
  jpdd            = 0
  jpdw            = 0
  jpcspf          = 0

CASE (i_ukca_chem_trop)

  ukca_config%l_ukca_chem     = .TRUE.
  ukca_config%l_ukca_trop     = .TRUE.
  ukca_config%ukca_int_method = int_method_BE_explicit
  IF (ukca_config%l_ukca_chem_aero) THEN
    ! Tropospheric chemistry with Aerosols (BE)
    ukca_config%l_ukca_aerchem = .TRUE.
    jpctr           = 33
    jpspec          = 53
    jpbk            = 88
    jptk            = 14
    jppj            = 20
    jphk            = 0
    jpdd            = ndry_aer_aer
    jpdw            = nwet_aer_aer

  ELSE
    ! Tropopspheric chemistry (BE) without aerosols
    jpctr            = 26
    jpspec           = 46
    jpbk             = 88
    jptk             = 14
    jppj             = 20
    jphk             = 0
    jpdd             = ndry_trop
    jpdw             = nwet_trop
  END IF

  ! Set No. of chemically active species equal to no. of tracers
  jpcspf          = jpctr

CASE (i_ukca_chem_raq)
  ! Regional Air Quality (BE)
  ukca_config%l_ukca_chem     = .TRUE.
  ukca_config%l_ukca_raq      = .TRUE.
  ukca_config%l_ukca_raqaero  = .FALSE.
  ukca_config%ukca_int_method = int_method_BE_explicit
  jpctr           = 40
  jpspec          = 58
  jpbk            = 113
  jptk            = 12
  jppj            = 23
  jphk            = 0   ! By default zero heterogeneous reactions
  jpdd            = ndry_raq
  jpdw            = nwet_raq

  ! Increase number of heterogeneous reactions if needed. For the moment
  ! only the heterog. hydrolysis of N2O5 is used, but more reactions
  ! can be added in the future.
  IF (ukca_config%l_ukca_classic_hetchem) THEN
    jphk = jphk + 2
  END IF

  IF (ukca_config%l_ukca_chem_aero) THEN
    IF (ukca_config%l_ukca_classic_hetchem) THEN
      cmessage='RAQ-AERO is not compatibile with the CLASSIC ' //              &
               'heterogeneous chemistry'
      errcode=113
      CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF
    ukca_config%l_ukca_raq     = .FALSE.
    ukca_config%l_ukca_raqaero = .TRUE.
    jpctr          = jpctr + 8
    jpspec         = jpspec + 8
    jpbk           = jpbk + 10
    jptk           = jptk + 2
    jppj           = jppj + 0
    jphk           = jphk + 2
    jpdd           = ndry_raqaero
    jpdw           = nwet_raqaero

  END IF

  ! Set No. of chemically active species equal to no. of tracers
  jpcspf          = jpctr

CASE (i_ukca_chem_offline_be)
  ! Offline oxidants aerosol precursor chemistry with backward-Euler solver
  ukca_config%l_ukca_chem       = .TRUE.
  ukca_config%l_ukca_offline_be = .TRUE.
  ukca_config%ukca_int_method   = int_method_BE_explicit
  jpctr           = 7
  jpspec          = 11
  jpbk            = 9
  jptk            = 1
  jppj            = 0
  jphk            = 3
  jpdd            = ndry_offline
  jpdw            = nwet_offline

  ! Set No. of chemically active species equal to no. of tracers
  jpcspf          = jpctr

CASE (i_ukca_chem_tropisop)
  ukca_config%l_ukca_chem     = .TRUE.
  ukca_config%l_ukca_tropisop = .TRUE.
  ukca_config%ukca_int_method = int_method_nr

  ! If aerosol chemistry on, add additional reactions
  IF (ukca_config%l_ukca_chem_aero) THEN
    ukca_config%l_ukca_achem     = .TRUE.
    ukca_config%l_ukca_nr_aqchem = .TRUE.

  ELSE
     ! Can't have trophet reactions without aerosol chemistry
    IF (ukca_config%l_ukca_trophet) THEN
      CALL umPrint( 'trop het chem requires aerosol chemistry on',             &
          src='ukca_setup_chem_mod')
      cmessage='Unsupported option choice'
      errcode=ABS(ukca_config%i_ukca_chem)
      CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF
  END IF

CASE (i_ukca_chem_strattrop)
  ukca_config%l_ukca_chem      = .TRUE.
  ukca_config%l_ukca_strattrop = .TRUE.
  ukca_config%l_ukca_advh2o    = .TRUE.
  ukca_config%ukca_int_method  = int_method_nr

  ! If aerosol chemistry on, add additional reactions
  IF (ukca_config%l_ukca_chem_aero) THEN
    ukca_config%l_ukca_achem     = .TRUE.
    ukca_config%l_ukca_nr_aqchem = .TRUE.

  ELSE
     ! Can't have trophet reactions without aerosol chemistry
    IF (ukca_config%l_ukca_trophet) THEN
      CALL umPrint( 'trop het chem requires aerosol chemistry on',             &
          src='ukca_setup_chem_mod')
      cmessage='Unsupported option choice'
      errcode=ABS(ukca_config%i_ukca_chem)
      CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF
  END IF

CASE (i_ukca_chem_strat)
  ukca_config%l_ukca_chem     = .TRUE.
  ukca_config%l_ukca_strat    = .TRUE.
  ukca_config%l_ukca_advh2o   = .TRUE.
  ukca_config%ukca_int_method = int_method_nr
  !  ! If aerosol chemistry on, add additional reactions
  IF (ukca_config%l_ukca_chem_aero) THEN
    ukca_config%l_ukca_achem     = .TRUE.
    ukca_config%l_ukca_nr_aqchem = .TRUE.

  END IF

  IF (ukca_config%l_ukca_trophet) THEN
    CALL umPrint( 'Strat chem does not support trop het chem',                 &
        src='ukca_setup_chem_mod')
    cmessage='Unsupported option choice'
    errcode=ABS(ukca_config%i_ukca_chem)
    CALL ereport('ukca_setup_chem',errcode,cmessage)
  END IF

CASE (i_ukca_chem_offline)
  ukca_config%l_ukca_chem      = .TRUE.
  ukca_config%l_ukca_offline   = .TRUE.
  ukca_config%l_ukca_nr_aqchem = .TRUE.
  ukca_config%l_ukca_advh2o    = .FALSE.
  ukca_config%ukca_int_method  = int_method_nr

  ! Case select for CRI+Stratosphere chemistry scheme
CASE (i_ukca_chem_cristrat)
  ukca_config%l_ukca_chem     = .TRUE.
  ukca_config%l_ukca_cristrat = .TRUE.
  ukca_config%l_ukca_advh2o   = .TRUE.
  ukca_config%ukca_int_method = int_method_nr

  ! If aerosol chemistry on, add additional reactions
  IF (ukca_config%l_ukca_chem_aero) THEN
    ukca_config%l_ukca_achem     = .TRUE.
    ukca_config%l_ukca_nr_aqchem = .TRUE.

  ELSE
     ! Can't have trophet reactions without aerosol chemistry
    IF (ukca_config%l_ukca_trophet) THEN
      CALL umPrint( 'trop het chem requires aerosol chemistry on',             &
          src='ukca_setup_chem_mod')
      cmessage='Unsupported option choice'
      errcode=ABS(ukca_config%i_ukca_chem)
      CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF
  END IF


CASE DEFAULT
  ! If we get here, we don't know how to deal with this chemistry scheme,
  ! stop
  WRITE(umMessage,'(A,1X,I6)') 'Unknown chemistry scheme: ',                   &
    ukca_config%i_ukca_chem
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  cmessage='Unknown chemistry scheme'
  errcode=ABS(ukca_config%i_ukca_chem)
  CALL ereport('ukca_setup_chem',errcode,cmessage)

END SELECT

! calculate jpnr from the sum of the reactions
jpnr = jpbk + jptk + jppj + jphk


! Error checking
IF (ukca_config%ukca_int_method == int_method_nr .OR.                          &
    ukca_config%i_ukca_chem ==  i_ukca_chem_offline_be) THEN
  ! The chemical timestep must be set for N-R and offline oxidants (BE) schemes
  IF (ukca_config%chem_timestep == imdi) THEN
    cmessage = ' The chemical timestep has not been set'
    errcode = 1
    CALL ereport('ukca_setup_chem',errcode,cmessage)
  END IF
END IF

IF (ukca_config%ukca_int_method  == int_method_BE_explicit .AND.               &
      ukca_config%i_ukca_chem /= i_ukca_chem_offline_be .AND.                  &
      ukca_config%i_ukca_chem /= i_ukca_chem_raq) THEN
    ! If aerosol chemistry or tropospheric het chem on abort - not supported
  IF (ukca_config%l_ukca_chem_aero .OR. ukca_config%l_ukca_trophet) THEN
    CALL umPrint( 'This BE scheme does not support add-on aerosol '            &
       // 'chemistry or tropospheric heterogeneous chemistry',                 &
       src='ukca_setup_chem_mod')
    cmessage='Unsupported option combination'
    errcode=ABS(ukca_config%i_ukca_chem)
    CALL ereport('ukca_setup_chem',errcode,cmessage)
  END IF
END IF

! Tropospheric heterogeneous chemistry requires GLOMAP-mode
IF (ukca_config%l_ukca_trophet .AND. .NOT. ukca_config%l_ukca_mode) THEN
  CALL umPrint( 'Need Glomap MODE on for '                                     &
     // 'tropospheric heterogeneous chemistry',src='ukca_setup_chem_mod')
  cmessage='Unsupported option combination'
  errcode=1
  CALL ereport('ukca_setup_chem',errcode,cmessage)
END IF

! Check chemistry supports GLOMAP-mode
IF (ukca_config%l_ukca_mode) THEN

   ! Newton Raphson
   ! Require that l_ukca_achem true or offline chemistry
  IF (ukca_config%ukca_int_method == int_method_nr) THEN
    IF (.NOT. ukca_config%l_ukca_nr_aqchem) THEN
      CALL umPrint( 'Need aerosol chemistry on for '                           &
     // 'Glomap MODE',src='ukca_setup_chem_mod')
      cmessage='Unsupported option combination'
      errcode=2
      CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF

    ! Backward Euler - only Trop + Aerosols, offline and RAQ-Aero support
    ! GLOMAP-mode
  ELSE
    IF (.NOT. (ukca_config%l_ukca_aerchem .OR.                                 &
               ukca_config%l_ukca_offline_be .OR.                              &
               ukca_config%l_ukca_raqaero)) THEN
      CALL umPrint( 'Need aerosol chemistry on for '                           &
     // 'Glomap MODE',src='ukca_setup_chem_mod')
      cmessage='Unsupported option combination'
      errcode=3
      CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF
  END IF

END IF

! Check specified DMS scheme logicals, only if marine DMS emissions are enabled
IF (ukca_config%l_ukca_chem_aero .AND. ukca_config%l_seawater_dms) THEN
  SELECT CASE (ukca_config%i_ukca_dms_flux )
  CASE (i_liss_merlivat)
    !    Do nothing
  CASE (i_wanninkhof)
    !    Do nothing
  CASE (i_nightingale)
    !    Do nothing
  CASE (i_blomquist)
    !    Do nothing

  CASE DEFAULT
    ! If not set or set to unknown value, throw up an error
    CALL umPrint(                                                              &
      'Marine DMS emissions expected but no UKCA DMS scheme is selected '      &
      // 'i_ukca_dms_flux should be 1,2,3 or 4')
    errcode = 4
    cmessage = 'RUN_UKCA: DMS flux scheme not specified'
    CALL ereport('UKCA_SETUP_CHEM_MOD',errcode,cmessage)
  END SELECT
END IF              ! chem_aero and seawater DMS

! Print the values we have set
IF (PrintStatus >= PrStatus_Diag) THEN
  CALL umPrint(       'Logicals and parameters set by ukca_setup_chem:',       &
      src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_chem      = ', ukca_config%l_ukca_chem
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_trop      = ', ukca_config%l_ukca_trop
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_achem   = ',   ukca_config%l_ukca_achem
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_nr_aqchem = ',                          &
    ukca_config%l_ukca_nr_aqchem
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_aerchem   = ',                          &
    ukca_config%l_ukca_aerchem
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_tropisop  = ',                          &
    ukca_config%l_ukca_tropisop
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_strattrop = ',                          &
    ukca_config%l_ukca_strattrop
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_strat     = ', ukca_config%l_ukca_strat
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_offline   = ', ukca_config%l_ukca_offline
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_cristrat  = ',                          &
    ukca_config%l_ukca_cristrat
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'ukca_int_method  = ',                          &
    ukca_config%ukca_int_method
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpctr            = ', jpctr
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpcspf           = ', jpcspf
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpspec           = ', jpspec
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpbk             = ', jpbk
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jptk             = ', jptk
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jppj             = ', jppj
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jphk             = ', jphk
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpnr             = ', jpnr
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpdd             = ', jpdd
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpdw             = ', jpdw
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I2)') 'i_ukca_dms_flux  = ',                          &
    ukca_config%i_ukca_dms_flux
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_setup_chem

END MODULE ukca_setup_chem_mod
