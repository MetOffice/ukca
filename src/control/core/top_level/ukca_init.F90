! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: To initialize internal values, addresses and
! other information needed for UKCA
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: Fortran
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_init_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_INIT_MOD'

CONTAINS

SUBROUTINE ukca_init

USE ukca_um_legacy_mod,    ONLY: isec_per_day, isec_per_hour
USE asad_mod,              ONLY: cdt, interval, ncsteps, tslimit

USE ukca_config_specification_mod, ONLY:                                       &
                                 ukca_config, glomap_config,                   &
                                 glomap_variables,                             &
                                 i_ukca_chem_off, int_method_nr,               &
                                 int_method_be_explicit,                       &
                                 int_method_impact,                            &
                                 i_suss_4mode,                                 &
                                 i_sussbcoc_5mode,                             &
                                 i_sussbcoc_4mode,                             &
                                 i_sussbcocso_5mode,                           &
                                 i_sussbcocso_4mode,                           &
                                 i_du_2mode,                                   &
                                 i_sussbcocdu_7mode,                           &
                                 i_sussbcocntnh_5mode_7cpt,                    &
                                 i_solinsol_6mode,                             &
                                 i_sussbcocduntnh_8mode_8cpt

USE ukca_mode_setup_interface_mod, ONLY: ukca_mode_setup_interface

USE ukca_constants,        ONLY: set_derived_constants
USE ukca_setup_chem_mod,   ONLY: ukca_setup_chem
USE ukca_config_defs_mod,  ONLY: n_mode_tracers

USE ukca_mode_setup,       ONLY: nmodes

USE ukca_setup_indices,    ONLY: ukca_indices_sv1,                             &
                                 ukca_indices_suss_4mode,                      &
                                 ukca_indices_orgv1_soto3,                     &
                                 ukca_indices_orgv1_soto3_isop,                &
                                 ukca_indices_sussbcoc_5mode,                  &
                                 ukca_indices_sussbcoc_5mode_isop,             &
                                 ukca_indices_sussbcoc_4mode,                  &
                                 ukca_indices_orgv1_soto6,                     &
                                 ukca_indices_sussbcocso_5mode,                &
                                 ukca_indices_sussbcocso_4mode,                &
                                 ukca_indices_duonly_2mode,                    &
                                 ukca_indices_sussbcocdu_7mode,                &
                                 ukca_indices_sussbcocntnh_5mode,              &
                                 ukca_indices_orgv1_soto3_solinsol,            &
                                 ukca_indices_solinsol_6mode,                  &
                                 ukca_indices_sussbcocduntnh_8mode_8cpt

USE umPrintMgr,            ONLY: umPrint, umMessage,                           &
                                 PrintStatus, PrStatus_Oper
USE ereport_mod,           ONLY: ereport
USE parkind1,              ONLY: jprb, jpim
USE yomhook,               ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Local variables

INTEGER                       :: imode     ! loop counter for modes
INTEGER                       :: icp       ! loop counter for components
INTEGER                       :: n_reqd_tracers ! no. of required tracers
INTEGER                       :: errcode=0     ! Error code: ereport
INTEGER                       :: timestep      ! Dynamical timestep
INTEGER, PARAMETER            :: ichem_ver132 = 132  ! To identify chemical vn
CHARACTER(LEN=errormessagelength)     :: cmessage=' '  ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check the logical switches
CALL check_settings(ukca_config, glomap_config)

! Set derived constants that depend on potentially user configurable values
CALL set_derived_constants()

! Set internal UKCA values for chemistry scheme
CALL ukca_setup_chem()

! Set up timestep counting. Interval depends on solver: Backward-Euler
! and other solvers are run every dynamical timestep. Newton-Raphson and
! Offline oxidants (backward-Euler) may be run every 2 / 3 timesteps for
! a 30/20 minutes dynamical timestep.

! dynamical timestep
timestep = INT(ukca_config%timestep)
IF (timestep <= 0 .OR. MOD(isec_per_hour, timestep) /= 0) THEN
  cmessage='No. of timesteps in an hour is not a positive integer'
  errcode=1
  CALL ereport(RoutineName,errcode,cmessage)
END IF
ukca_config%timesteps_per_day = isec_per_day / timestep
ukca_config%timesteps_per_hour = isec_per_hour / timestep

! Do not check for solver type and timestep if none of the chemistry
! schemes is selected e.g. for an Age-of-air-only configuration.
! In that case, set values to default in case they are used elsewhere
!$OMP PARALLEL
IF ( ukca_config%i_ukca_chem == i_ukca_chem_off ) THEN
  interval = 1
  cdt = REAL(timestep)
  ncsteps = 1
ELSE
  ! Calculate interval depending on chemistry scheme
  ukca_method_type:IF (ukca_config%ukca_int_method == int_method_nr) THEN
    ! Newton-Raphson solver
    interval = ukca_config%chem_timestep/timestep
    cdt = REAL(ukca_config%chem_timestep)
    ncsteps = 1
  ELSE IF (ukca_config%ukca_int_method == int_method_impact) THEN
    ! IMPACT solver use about 15 or 10 minutes, depending on dynamical timestep
    IF (timestep < tslimit) THEN
      cdt = REAL(timestep)
      ncsteps = 1
      interval = 1
    ELSE
      cdt = 0.5*REAL(timestep)
      ncsteps = 2
      interval = 1
    END IF
  ELSE IF (ukca_config%ukca_int_method == int_method_be_explicit) THEN
    ! explicit Backward-Euler solver
    ! solver interval derived from namelist value of chemical timestep
    interval = ukca_config%chem_timestep/timestep
    cdt = REAL(ukca_config%chem_timestep)
    ncsteps = 1
  ELSE
    ! Unknown solver type
    WRITE(cmessage, '(A,I0,A)')                                                &
      'Type of solver (ukca_int_method = ',ukca_config%ukca_int_method,        &
      ') not recognised.'
    errcode = 2
    CALL ereport('UKCA_INIT',errcode,cmessage)
  END IF ukca_method_type

  IF (printstatus >= prstatus_oper) THEN
    WRITE(umMessage,'(A40,I6)') 'Interval for chemical solver set to: ',       &
                                                               interval
    CALL umPrint(umMessage,src='ukca_init')
    WRITE(umMessage,'(A40,E12.4)') 'Timestep for chemical solver set to: ',    &
                                                               cdt
    CALL umPrint(umMessage,src='ukca_init')
    WRITE(umMessage,'(A40,I6)') 'No. steps for chemical solver set to: ',      &
                                                               ncsteps
    CALL umPrint(umMessage,src='ukca_init')
  END IF

  ! Verify that the interval and timestep values have been set correctly
  IF (ABS(cdt*ncsteps - REAL(timestep*interval)) > 1e-4) THEN
    cmessage=' chemical timestep does not fit dynamical timestep'
    WRITE(umMessage,'(A)') cmessage
    CALL umPrint(umMessage,src='ukca_init')
    WRITE(umMessage,'(A,I6,A,I6)') ' timestep: ',timestep,                     &
                                   ' interval: ',interval
    CALL umPrint(umMessage,src='ukca_init')
    errcode = ukca_config%chem_timestep
    CALL ereport('UKCA_INIT',errcode,cmessage)
  END IF

END IF
!$OMP END PARALLEL

IF (ukca_config%l_ukca_mode) THEN
  ! Call appropriate MODE setup routine
  IF ( ukca_config%i_ukca_chem_version >= ichem_ver132) THEN
    IF (glomap_config%i_mode_setup == i_sussbcoc_5mode) THEN
      CALL ukca_indices_orgv1_soto3_isop
      CALL ukca_indices_sussbcoc_5mode_isop
    ELSE
      errcode     = 4
      cmessage    = 'Isoprene SOA (i_chem_version > 132)   '  //               &
                   'only works with i_mode_setup=2 '
      CALL ereport('UKCA_INIT', errcode, cmessage)
    END IF
  ELSE
    IF ( glomap_config%i_mode_setup == i_suss_4mode ) THEN ! 1
      CALL ukca_indices_sv1
      CALL ukca_indices_suss_4mode
    ELSE IF ( glomap_config%i_mode_setup == i_sussbcoc_5mode ) THEN ! 2
      CALL ukca_indices_orgv1_soto3
      CALL ukca_indices_sussbcoc_5mode
    ELSE IF ( glomap_config%i_mode_setup == i_sussbcoc_4mode ) THEN ! 3
      CALL ukca_indices_orgv1_soto3
      CALL ukca_indices_sussbcoc_4mode
    ELSE IF ( glomap_config%i_mode_setup == i_sussbcocso_5mode ) THEN ! 4
      CALL ukca_indices_orgv1_soto6
      CALL ukca_indices_sussbcocso_5mode
    ELSE IF ( glomap_config%i_mode_setup == i_sussbcocso_4mode ) THEN ! 5
      CALL ukca_indices_orgv1_soto6
      CALL ukca_indices_sussbcocso_4mode
    ELSE IF ( glomap_config%i_mode_setup == i_du_2mode ) THEN ! 6
      !!    CALL ukca_indices_nochem
      !! temporarily run 2-mode dust only with chemistry, though it's not needed
      CALL ukca_indices_orgv1_soto3
      CALL ukca_indices_duonly_2mode
      !!  ELSE IF ( glomap_config%i_mode_setup == 7 ) THEN
        !!    CALL ukca_indices_nochem
        !!    CALL ukca_indices_duonly_3mode
    ELSE IF ( glomap_config%i_mode_setup == i_sussbcocdu_7mode ) THEN ! 8
      CALL ukca_indices_orgv1_soto3
      CALL ukca_indices_sussbcocdu_7mode
      !!  ELSE IF ( glomap_config%i_mode_setup == 9 ) THEN
        !!    CALL ukca_indices_orgv1_soto3
        !!    CALL ukca_indices_sussbcocdu_4mode
    ELSE IF ( glomap_config%i_mode_setup == i_sussbcocntnh_5mode_7cpt ) THEN ! 10
      CALL ukca_indices_orgv1_soto3
      CALL ukca_indices_sussbcocntnh_5mode
    ELSE IF ( glomap_config%i_mode_setup == i_solinsol_6mode ) THEN ! 11
      CALL ukca_indices_orgv1_soto3_solinsol
      CALL ukca_indices_solinsol_6mode
    ELSE IF (glomap_config%i_mode_setup == i_sussbcocduntnh_8mode_8cpt) THEN ! 12
      CALL ukca_indices_orgv1_soto3
      CALL ukca_indices_sussbcocduntnh_8mode_8cpt
    ELSE
      cmessage=' i_mode_setup has unrecognised value'
      WRITE(umMessage,'(A,I4)') cmessage,glomap_config%i_mode_setup
      CALL umPrint(umMessage,src='ukca_init')
      errcode = 4
      CALL ereport('UKCA_INIT',errcode,cmessage)
    END IF
  END IF       ! i_mode_setup

  CALL ukca_mode_setup_interface ( glomap_config%i_mode_setup,                 &
                                   glomap_config%l_ukca_radaer,                &
                                   glomap_config%i_ukca_tune_bc,               &
                                   glomap_config%l_fix_nacl_density,           &
                                   glomap_config%l_fix_ukca_hygroscopicities )

  ! Calculate number of aerosol tracers required for components and number
  n_reqd_tracers = 0
  DO imode=1,nmodes
    IF (glomap_variables%mode(imode)) THEN
      DO icp=1,glomap_variables%ncp
        IF (glomap_variables%component(imode,icp)) THEN
          n_reqd_tracers = n_reqd_tracers + 1
        END IF
      END DO
    END IF
  END DO
  n_mode_tracers = n_reqd_tracers + SUM(glomap_variables%mode_choice)

ELSE
  ! Allocate arrays that are referred to outside GLOMAP
  ! to avoid compiler errors
  IF (.NOT. ALLOCATED(glomap_variables%component)) THEN
    ALLOCATE(glomap_variables%component(1,1))
  END IF

END IF    ! ukca_config%_ukca_mode

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_init


SUBROUTINE check_settings(ukca_config, glomap_config)

! Description:
!   Subroutine to apply logic checks based on the UKCA
!   options selected and perform other validation of configuration variables.
!   May apply automatic adjustments of some variables if out of range.
!   Note that some chemistry scheme specific checks
!   are done in 'ukca_setup_chem'

USE ukca_um_legacy_mod,     ONLY: l_um_infrastructure, l_um_emissions_updates

USE ukca_config_specification_mod, ONLY:                                       &
                                 ukca_config_spec_type,                        &
                                 glomap_config_spec_type,                      &
                                 i_ukca_chem_off,                              &
                                 i_ukca_chem_raq,                              &
                                 i_ukca_chem_offline_be,                       &
                                 i_ukca_chem_tropisop,                         &
                                 i_ukca_chem_strattrop,                        &
                                 i_ukca_chem_strat,                            &
                                 i_ukca_chem_offline,                          &
                                 i_ukca_chem_cristrat,                         &
                                 i_age_reset_by_level,                         &
                                 i_age_reset_by_height,                        &
                                 i_ukca_nophot,                                &
                                 i_ukca_phot2d,                                &
                                 i_ukca_photol_strat,                          &
                                 i_ukca_fastjx,                                &
                                 i_light_param_off,                            &
                                 i_light_param_pr,                             &
                                 i_light_param_luhar,                          &
                                 i_ukca_activation_arg,                        &
                                 i_top_BC,                                     &
                                 i_top_BC_H2O,                                 &
                                 bl_tracer_mix

USE asad_mod,              ONLY: nrsteps_max

USE umPrintMgr,            ONLY: PrStatus_Normal,PrintStatus,newline,          &
                                 umPrint
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod,           ONLY: ereport
USE parkind1,              ONLY: jprb, jpim
USE yomhook,               ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine arguments
TYPE(ukca_config_spec_type), INTENT(IN OUT) :: ukca_config
TYPE(glomap_config_spec_type), INTENT(IN OUT) :: glomap_config

! Local variables

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_SETTINGS'
CHARACTER (LEN=errormessagelength) :: cmessage   ! Error message
INTEGER                            :: errcode    ! Variable passed to ereport

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

errcode = 0                   ! Initialise

! Check photolysis switches
IF ( ukca_config%i_ukca_photol /= i_ukca_nophot .AND.                          &
     ukca_config%i_ukca_photol /= i_ukca_photol_strat .AND.                    &
     ukca_config%i_ukca_photol /= i_ukca_phot2d .AND.                          &
     ukca_config%i_ukca_photol /= i_ukca_fastjx .AND.                          &
     ukca_config%i_ukca_chem /= i_ukca_chem_off ) THEN
  cmessage='Unknown photolysis scheme.'
  errcode=ABS(ukca_config%i_ukca_photol)
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! biom_aer_ems_scaling - valid range 0. - 30.
IF ( (glomap_config%biom_aer_ems_scaling < 0.0 .OR.                            &
      glomap_config%biom_aer_ems_scaling > 30.0) .AND.                         &
     glomap_config%l_ukca_scale_biom_aer_ems ) THEN
  cmessage='biom_aer_ems_scaling should be between 0.0 - 30.0'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! seadms_ems_scaling - valid range 0. - 10.
IF ( (ukca_config%seadms_ems_scaling < 0.0 .OR.                                &
      ukca_config%seadms_ems_scaling > 10.0) .AND.                             &
     ukca_config%l_ukca_scale_seadms_ems ) THEN
  cmessage='seadms_ems_scaling should be between 0.0 - 10.0'
  errcode = 2
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! sea_salt_ems_scaling - valid range 0. - 10.
IF ( (glomap_config%sea_salt_ems_scaling < 0.0 .OR.                            &
      glomap_config%sea_salt_ems_scaling > 10.0) .AND.                         &
     glomap_config%l_ukca_scale_sea_salt_ems ) THEN
  cmessage='sea_salt_ems_scaling should be between 0.0 - 10.0'
  errcode = 2
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! marine_pom_ems_scaling - valid range 0. - 10.
IF ( (glomap_config%marine_pom_ems_scaling < 0.0 .OR.                          &
      glomap_config%marine_pom_ems_scaling > 10.0) .AND.                       &
     glomap_config%l_ukca_scale_marine_pom_ems ) THEN
  cmessage='marine_pom_ems_scaling should be between 0.0 - 10.0'
  errcode = 2
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! soa_yield_scaling - valid range 0. - 5.
IF ( (ukca_config%soa_yield_scaling_mt < 0.0 .OR.                              &
      ukca_config%soa_yield_scaling_mt > 5.0) .AND.                            &
     ukca_config%l_ukca_scale_soa_yield_mt) THEN
  cmessage='soa_yield_scaling should be between 0.0 - 5.0'
  errcode = 3
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! soa_yield_scaling_isop - valid range 0. - 10.
IF ( (ukca_config%soa_yield_scaling_isop < 0.0 .OR.                            &
      ukca_config%soa_yield_scaling_isop > 10.0) .AND.                         &
     ukca_config%l_ukca_scale_soa_yield_isop ) THEN
  cmessage='soa_yield_scaling_isop should be between 0.0 -10.0'
  errcode = 37
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! mode_activation_dryr - valid range 20. - 100.
IF ( (glomap_config%mode_activation_dryr < 20.0 .OR.                           &
      glomap_config%mode_activation_dryr > 100.0) .AND.                        &
     ukca_config%l_ukca_mode) THEN
  cmessage=' mode_activation_dryr should be between 20.0 and 100.0'
  errcode = 4
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! mode_incld_so2_rfrac - valid range 0. - 1.
IF ( (glomap_config%mode_incld_so2_rfrac < 0.0 .OR.                            &
      glomap_config%mode_incld_so2_rfrac > 1.0) .AND.                          &
     ukca_config%l_ukca_mode) THEN
  cmessage=' mode_incld_so2_rfrac should be between 0.0 and 1.0'
  errcode = 5
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! If lightning emissions of NOx are on with chemistry check that a full
! chemistry scheme is selected and validate lightning scheme configuration
IF (ukca_config%i_ukca_chem /= i_ukca_chem_off .AND.                           &
    (.NOT. ukca_config%l_ukca_emissions_off) .AND.                             &
    ukca_config%i_ukca_light_param /= i_light_param_off) THEN
  ! Check the chemistry scheme is not an Offline scheme
  IF (ukca_config%i_ukca_chem == i_ukca_chem_offline .OR.                      &
      ukca_config%i_ukca_chem == i_ukca_chem_offline_be) THEN
    cmessage = 'Lightning NOx emissions cannot be used with an Offline scheme'
    errcode = 6
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
  ! Check values of Lightning NOX scaling factor and parameterisation.
  IF (ukca_config%lightnox_scale_fac < 0.0 .OR.                                &
      ukca_config%lightnox_scale_fac > 10.0) THEN
    cmessage=' Light-NOx scale factor should be between 0.0 and 10.0'
    errcode = 7
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
  IF (ukca_config%i_ukca_light_param < 1 .OR.                                  &
      ukca_config%i_ukca_light_param > 3) THEN
    cmessage=' Light-NOx parameterision should be 1, 2 or 3'
    errcode = 8
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF

! Warning that convective scavenging won't work with no convection scheme
IF ( (.NOT. ukca_config%l_param_conv) .AND.                                    &
     (ukca_config%i_ukca_chem /= i_ukca_chem_off)) THEN
  WRITE(cmessage,'(A)')' Convection parametrisation not used, no '             &
      //newline//                                                              &
      'convective scavenging of chemical species will be considered.'
  errcode=-1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Check the values of parameters that control resetting of age-of-air tracer
IF ( ukca_config%l_ukca_ageair ) THEN
  SELECT CASE(ukca_config%i_ageair_reset_method)
  CASE (i_age_reset_by_level)
    IF ( ukca_config%max_ageair_reset_level < 1 ) THEN
      errcode = 9
      WRITE(cmessage,'(A,I0)')'Inconsistent Age-of-air reset level: ',         &
          ukca_config%max_ageair_reset_level
    END IF
  CASE (i_age_reset_by_height)
    IF ( ukca_config%max_ageair_reset_height < 0.0 .OR.                        &
          ukca_config%max_ageair_reset_height > 20000.0 ) THEN
      errcode = 10
      WRITE(cmessage,'(A,F16.4)')'Inconsistent Age-of-air reset height: ',     &
          ukca_config%max_ageair_reset_height
    END IF
  CASE DEFAULT
    errcode = 11
    WRITE(cmessage,'(A,I0)') ' Inconsistent Age-of-air reset method: ',        &
           ukca_config%i_ageair_reset_method
  END SELECT
  IF ( errcode > 0 ) CALL ereport(RoutineName,errcode,cmessage)

END IF     ! l_ukca_ageair

! Check settings for Newton-Raphson schemes
IF ( ukca_config%i_ukca_chem == i_ukca_chem_strat .OR.                         &
     ukca_config%i_ukca_chem == i_ukca_chem_strattrop .OR.                     &
     ukca_config%i_ukca_chem == i_ukca_chem_tropisop .OR.                      &
     ukca_config%i_ukca_chem == i_ukca_chem_offline .OR.                       &
     ukca_config%i_ukca_chem == i_ukca_chem_cristrat) THEN
  ! check number of N-R iterations
  IF (ukca_config%nrsteps < 0 .OR. ukca_config%nrsteps > nrsteps_max) THEN
    WRITE(cmessage,'(2(A,1x,I0))') 'NRSTEPS is not in range - 0:',nrsteps_max, &
          '. NRSTEPS=',ukca_config%nrsteps
    errcode = 12
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
  ! check settings in quasi-Newton step to ensure they are sensible
  IF (ukca_config%l_ukca_quasinewton) THEN
    IF (ukca_config%i_ukca_quasinewton_end <                                   &
        ukca_config%i_ukca_quasinewton_start) THEN
      WRITE(cmessage,'(A,A,I4,I4)')                                            &
           ' i_ukca_quasinewton_start must be less than or equal ',            &
           ' to i_ukca_quasinewton_end', ukca_config%i_ukca_quasinewton_start, &
           ukca_config%i_ukca_quasinewton_end
      errcode = 13
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
    IF ((ukca_config%i_ukca_quasinewton_start < 2) .OR.                        &
        (ukca_config%i_ukca_quasinewton_start > 50)) THEN
      WRITE(cmessage,'(A,I4)')                                                 &
           ' i_ukca_quasinewton_start must be between 2 & 50',                 &
           ukca_config%i_ukca_quasinewton_start
      errcode = 14
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
    IF ((ukca_config%i_ukca_quasinewton_end < 2) .OR.                          &
        (ukca_config%i_ukca_quasinewton_end > 50)) THEN
      WRITE(cmessage,'(A,I4)')                                                 &
           ' i_ukca_quasinewton_end must be between 2 & 50',                   &
           ukca_config%i_ukca_quasinewton_end
      errcode = 15
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  END IF
END IF

IF ( .NOT. (ukca_config%i_ukca_chem == i_ukca_chem_strat .OR.                  &
      ukca_config%i_ukca_chem == i_ukca_chem_strattrop .OR.                    &
      ukca_config%i_ukca_chem == i_ukca_chem_tropisop .OR.                     &
      ukca_config%i_ukca_chem == i_ukca_chem_offline .OR.                      &
      ukca_config%i_ukca_chem == i_ukca_chem_cristrat)) THEN
  ! check settings in quasi-Newton step to ensure they are sensible
  IF (ukca_config%l_ukca_asad_columns) THEN
    errcode = 16
    WRITE(cmessage,'(A,L1)')                                                   &
         ' Column-call can only be for Newton-Raphson schemes: ',              &
         ukca_config%l_ukca_asad_columns
  END IF
END IF

IF ( .NOT. (ukca_config%i_ukca_chem == i_ukca_chem_strat .OR.                  &
      ukca_config%i_ukca_chem == i_ukca_chem_strattrop .OR.                    &
      ukca_config%i_ukca_chem == i_ukca_chem_tropisop .OR.                     &
      ukca_config%i_ukca_chem == i_ukca_chem_offline .OR.                      &
      ukca_config%i_ukca_chem == i_ukca_chem_cristrat)) THEN
  IF (ukca_config%l_ukca_debug_asad) THEN
    errcode = 17
    WRITE(cmessage,'(A,L1)')                                                   &
         ' ASAD debugging can only be for Newton-Raphson schemes: ',           &
         ukca_config%l_ukca_debug_asad
  END IF
END IF

! Check validity of top boundary option for a stratospheric scheme
IF ((ukca_config%i_ukca_chem == i_ukca_chem_strattrop .OR.                     &
     ukca_config%i_ukca_chem == i_ukca_chem_strat .OR.                         &
     ukca_config%i_ukca_chem == i_ukca_chem_cristrat) .AND.                    &
    (ukca_config%i_ukca_topboundary < 0 .OR.                                   &
     ukca_config%i_ukca_topboundary > 4)) THEN
  WRITE(cmessage,'(A,I0,A)')                                                   &
       ' Incorrect value for i_ukca_topboundary ( ',                           &
       ukca_config%i_ukca_topboundary,                                         &
       ' ) - should be between 0 and 4'
  errcode = 18
  CALL ereport(routinename,errcode,cmessage)
END IF

IF ((.NOT. ukca_config%l_ukca_h2o_feedback) .AND.                              &
    (ukca_config%i_ukca_topboundary == i_top_BC_H2O)) THEN
  ! Cannot impose top boundary condition if H2O is not interactive.
  cmessage = 'Cannot impose top boundary for H2O if it is not interactive.'
  errcode = 19
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF ((ukca_config%i_ukca_topboundary == i_top_BC_H2O) .AND.                     &
    (ukca_config%i_ukca_chem /= i_ukca_chem_strat) .AND.                       &
    (ukca_config%i_ukca_chem /= i_ukca_chem_strattrop)) THEN
  ! Cannot conserve hydrogen because H, OH, or H2 don't exist
  cmessage='Cannot impose top boundary for H2O if H, OH, or H2 do not exist.'
  errcode = 20
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF ((ukca_config%i_ukca_topboundary >= i_top_BC) .AND.                         &
    (ukca_config%i_ukca_chem /= i_ukca_chem_strattrop) .AND.                   &
    (ukca_config%i_ukca_chem /= i_ukca_chem_strat) .AND.                       &
    (ukca_config%i_ukca_chem /= i_ukca_chem_cristrat)) THEN
  ! Cannot impose top boundary condition for species
  cmessage='Can only impose top boundary for O3, NO and CO in' // newline //   &
           'schemes with stratospheric chemistry'
  errcode = 21
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF (ukca_config%l_ukca_inferno_ch4 .AND. ukca_config%l_ukca_prescribech4) THEN
  WRITE(cmessage,'(A,A,L1)')                                                   &
       ' l_ukca_inferno_ch4 is .true and l_ukca_prescribech4 is .true.',       &
       ' CH4 interactive fire emissions are not compatible with prescribed'
  errcode = 22
  CALL ereport(routinename,errcode,cmessage)
END IF

! Make sure RO2 transport option is only used with StratTrop mechanism
IF (ukca_config%l_ukca_ro2_ntp) THEN
  IF ((ukca_config%i_ukca_chem == i_ukca_chem_strattrop) .OR.                  &
      (ukca_config%i_ukca_chem == i_ukca_chem_cristrat)) THEN
    IF ( PrintStatus > PrStatus_Normal )                                       &
      CALL umPrint('Transport of peroxy-radicals turned off.',                 &
        src=RoutineName)
  ELSE
    WRITE(cmessage,'(A)')                                                      &
      ' l_ukca_ro2_ntp can only be T with StratTrop or CRI-Strat chemistry'
    errcode = 23
    CALL ereport(routinename,errcode,cmessage)
  END IF
ELSE
  ! Model will crash if NOT using RO2 NTP with CRI
  IF (ukca_config%i_ukca_chem == i_ukca_chem_cristrat) THEN
    WRITE(cmessage,'(A)')                                                      &
        ' l_ukca_ro2_ntp must be TRUE if running with CRI chemistry'
    errcode = 24
    CALL ereport(routinename,errcode,cmessage)
  END IF
END IF

! Make sure RO2-permutation reactions are only used with StratTrop
IF (ukca_config%l_ukca_ro2_perm) THEN
  IF ((ukca_config%i_ukca_chem == i_ukca_chem_strattrop) .OR.                  &
       (ukca_config%i_ukca_chem == i_ukca_chem_cristrat)) THEN
    IF ( PrintStatus > PrStatus_Normal )                                       &
      CALL umPrint('RO2-permutation reactions turned on.',                     &
        src=RoutineName)
  ELSE
    WRITE(cmessage,'(A)')                                                      &
      ' l_ukca_ro2_perm can only be T with StratTrop or CRI-Strat chemistry'
    errcode = 25
    CALL ereport(routinename,errcode,cmessage)
  END IF
ELSE
  ! Model will crash if NOT using RO2-permutation chemistry with CRI
  IF (ukca_config%i_ukca_chem == i_ukca_chem_cristrat) THEN
    WRITE(cmessage,'(A)')                                                      &
       ' l_ukca_ro2_perm must be TRUE if running with CRI chemistry'
    errcode = 26
    CALL ereport(routinename,errcode,cmessage)
  END IF
END IF

IF (ukca_config%l_ukca_ddepo3_ocean .AND. .NOT. ukca_config%l_ukca_intdd) THEN
  WRITE(cmessage,'(A)')                                                        &
       ' l_ukca_ddepo3_ocean is .true. but l_ukca_intdd is .false.'
  errcode = 27
  CALL ereport(routinename,errcode,cmessage)
END IF

! JULES-based deposition only available for interactive dry deposition
IF (ukca_config%l_deposition_jules .AND. .NOT. ukca_config%l_ukca_intdd) THEN
  WRITE(cmessage,'(A)')                                                        &
       ' l_deposition_jules is .true. but l_ukca_intdd is .false.'
  errcode = 28
  CALL ereport(routinename,errcode,cmessage)
END IF

IF ((glomap_config%i_ukca_activation_scheme == i_ukca_activation_arg) .AND.    &
    ((glomap_config%i_ukca_nwbins < 1) .OR.                                    &
     (glomap_config%i_ukca_nwbins > 20))) THEN
  cmessage = 'Cannot set i_ukca_nwbins less than one or greater than 20.'
  errcode  = 29
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Check settings for i_ukca_chem_version
! - only do for Newton-Raphson schemes
IF ( ukca_config%i_ukca_chem == i_ukca_chem_strat      .OR.                    &
     ukca_config%i_ukca_chem == i_ukca_chem_strattrop  .OR.                    &
     ukca_config%i_ukca_chem == i_ukca_chem_tropisop   .OR.                    &
     ukca_config%i_ukca_chem == i_ukca_chem_offline    .OR.                    &
     ukca_config%i_ukca_chem == i_ukca_chem_cristrat ) THEN
  IF ( ukca_config%i_ukca_chem_version < 107 ) THEN
    cmessage = 'Value of i_ukca_chem_version cannot be less than 107.'
    errcode  = 30
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF

! Check configuration supports heterogeneous chemistry on CLASSIC aerosols if
! selected
IF (ukca_config%l_ukca_classic_hetchem .AND.                                   &
    ukca_config%i_ukca_chem /= i_ukca_chem_raq) THEN
  cmessage = 'Heterogeneous chemistry on CLASSIC aerosols requires RAQ'
  errcode  = 31
  CALL ereport(RoutineName,errcode,cmessage)
END IF
IF (ukca_config%l_ukca_classic_hetchem .AND. ukca_config%l_ukca_chem_aero) THEN
  cmessage = 'Heterogeneous chemistry is not supported in RAQ-Aero'
  errcode  = 32
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Check configuration supports PSC heterogeneous chemistry if selected
IF (ukca_config%l_ukca_het_psc .AND.                                           &
    ukca_config%i_ukca_chem /= i_ukca_chem_strattrop .AND.                     &
    ukca_config%i_ukca_chem /= i_ukca_chem_strat .AND.                         &
    ukca_config%i_ukca_chem /= i_ukca_chem_cristrat) THEN
  cmessage = 'PSC heterogeneous chemistry requires a scheme with'              &
             // newline // 'stratospheric chemistry'
  errcode  = 33
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Check configuration supports tropospheric heterogeneous chemistry if selected
IF (ukca_config%l_ukca_trophet .AND.                                           &
    ukca_config%i_ukca_chem /= i_ukca_chem_tropisop .AND.                      &
    ukca_config%i_ukca_chem /= i_ukca_chem_strattrop .AND.                     &
    ukca_config%i_ukca_chem /= i_ukca_chem_cristrat) THEN
  cmessage = 'Tropospheric heterogeneous chemistry requires an N-R scheme with'&
             // newline // 'tropospheric chemistry'
  errcode  = 34
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Check that any CLASSIC aerosol species selected are supported by an
! appropriate heterogeneous chemistry scheme
IF (ukca_config%l_use_classic_so4 .AND.                                        &
   .NOT. (ukca_config%l_ukca_het_psc .OR.                                      &
          ukca_config%l_ukca_classic_hetchem)) THEN
  cmessage = 'CLASSIC SO4 can only be used with PSC heterogeneous chemistry'   &
             // newline // 'or heterogeneous chemistry on CLASSIC aerosols'
  errcode  = 35
  CALL ereport(RoutineName,errcode,cmessage)
END IF
IF ((ukca_config%l_use_classic_soot .OR. ukca_config%l_use_classic_ocff .OR.   &
     ukca_config%l_use_classic_biogenic .OR.                                   &
     ukca_config%l_use_classic_seasalt) .AND.                                  &
     .NOT. ukca_config%l_ukca_classic_hetchem) THEN
  cmessage = 'CLASSIC aerosols other than SO4 can only be used with'           &
             // newline // 'heterogeneous chemistry on CLASSIC aerosols'
  errcode  = 36
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Check that use of actual gridbox mass is enabled if UM diagnostics are
! enabled with a chemistry scheme being active or if a stratospheric chemistry
! scheme is used without turning off emissions (since it is then needed in
! stratospheric schemes for lower boundary conditions).
IF ((.NOT. ukca_config%l_use_gridbox_mass) .AND.                               &
    ((ukca_config%i_ukca_chem /= i_ukca_chem_off .AND.                         &
      ukca_config%l_enable_diag_um) .OR.                                       &
     ((ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_stratcfc .OR.          &
      ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat) .AND.     &
      .NOT. ukca_config%l_ukca_emissions_off))) THEN
  cmessage = 'Cannot run without mass of air in grid box as it is needed'      &
             // newline // 'for UM diagnostics and/or a stratospheric'         &
             // newline // 'chemistry scheme'
  errcode  = 37
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! check that cloud pH fitting parameters are correct
! a and b should all be between -10 to 10
! y-intercept should be between 0 and 14
IF (ukca_config%l_ukca_intph) THEN
  ! Cloud pH Fitting parameter a
  IF (ukca_config%ph_fit_coeff_a < -10.0 .OR.                                  &
      ukca_config%ph_fit_coeff_a > 10.0) THEN
    cmessage='Check ph_fit_coeff_a value as should be between -10.0 - 10.0'
    errcode = 38
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
  ! Cloud pH Fitting parameter b
  IF (ukca_config%ph_fit_coeff_b < -10.0 .OR.                                  &
      ukca_config%ph_fit_coeff_b > 10.0) THEN
    cmessage='Check ph_fit_coeff_b value as should be between -10.0 - 10.0'
    errcode = 39
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
  ! Cloud pH intercept
  IF (ukca_config%ph_fit_intercept < 0.0 .OR.                                  &
      ukca_config%ph_fit_intercept > 14.0) THEN
    cmessage='Check ph_fit_intercept as should be between 0.0 - 14.0'
    errcode = 40
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF

! If UM infrastructure code is not available, check that it is not needed.
IF (.NOT. l_um_infrastructure) THEN
  IF (ukca_config%i_ukca_light_param == i_light_param_pr .OR.                  &
      ukca_config%i_ukca_light_param == i_light_param_luhar) THEN
    cmessage = 'Lightning NOx parameterization is only supported for UM grids'
    errcode  = 41
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
  IF (ukca_config%l_enable_diag_um) THEN
    cmessage = 'No support for UM diagnostics outside the UM'
    errcode  = 42
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
  IF (ukca_config%l_use_gridbox_mass) THEN
    cmessage = 'Calculation of gridbox mass is only supported for UM grids'
    errcode  = 43
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
  IF (ukca_config%l_environ_z_top) THEN
    cmessage = 'Override of top-of-model height is not allowed outside the UM'
    errcode  = 44
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF

! If code to support emissions updating of tracer is not available,
! check that it is not needed.
IF ((.NOT. l_um_emissions_updates) .OR. (.NOT. ASSOCIATED(bl_tracer_mix))) THEN
  IF (ukca_config%i_ukca_chem /= i_ukca_chem_off .AND.                         &
      (.NOT. (ukca_config%l_ukca_emissions_off .OR.                            &
             ukca_config%l_suppress_ems))) THEN
    cmessage = 'No support code available for tracer updates from emissions'
    errcode  = 45
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF

IF ( (glomap_config%i_ukca_tune_bc < 0 .OR.                                    &
      glomap_config%i_ukca_tune_bc > 2)  .AND.                                 &
     glomap_config%l_ukca_radaer) THEN
  cmessage='i_ukca_tune_bc should be 0, 1 or 2'
  errcode = 46
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Fixed tropopause level - valid range 1 - model_levels
IF (ukca_config%l_fix_tropopause_level) THEN
  IF (ukca_config%fixed_tropopause_level < 1 .OR.                              &
      ukca_config%fixed_tropopause_level > ukca_config%model_levels) THEN
    cmessage='fixed_tropopause_level is out of range 1 - model_levels'
    errcode = 47
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF

! Check that use of gridbox volume is enabled if UM diagnostics are
! enabled with a chemistry scheme being active.
IF ((.NOT. ukca_config%l_use_gridbox_volume) .AND.                             &
    ukca_config%i_ukca_chem /= i_ukca_chem_off .AND.                           &
    ukca_config%l_enable_diag_um) THEN
  cmessage = 'Cannot run without grid box volume as it is needed'              &
             // newline // 'for UM diagnostics'
  errcode  = 48
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF (ukca_config%l_ukca_so2ems_plumeria .AND.                                   &
    .NOT. ukca_config%l_ukca_so2ems_expvolc) THEN
  cmessage='l_ukca_so2ems_plumeria is true but l_ukca_so2ems_expvolc is false'
  errcode = 49
  CALL ereport(routinename,errcode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_settings

END MODULE ukca_init_mod
