!*****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module containing subroutine photol_setup for receiving configuration
!   options for photolysis and setting up internal data to
!   define the configuration.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_Photolysis
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE photol_setup_mod

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='PHOTOL_SETUP_MOD'

! Public procedures
PUBLIC :: photol_setup

CONTAINS

! ----------------------------------------------------------------------

SUBROUTINE photol_setup(i_photol_scheme,                                       &
                      error_code,                                              &
                      chem_timestep,                                           &
                      fastjx_mode,                                             &
                      fastjx_numwl,                                            &
                      global_row_length,                                       &
                      i_solcylc_type,                                          &
                      ip_aerosol_param_moist,                                  &
                      ip_accum_sulphate,                                       &
                      ip_aitken_sulphate,                                      &
                      model_levels,                                            &
                      n_cca_lev,                                               &
                      solcylc_start_year,                                      &
                      i_error_method,                                          &
                      n_phot_spc,                                              &
                      njval,                                                   &
                      nw1,                                                     &
                      nw2,                                                     &
                      jtaumx,                                                  &
                      naa,                                                     &
                      n_solcyc_ts,                                             &
                      jind,                                                    &
                      l_cal360,                                                &
                      l_cloud_pc2,                                             &
                      l_3d_cca,                                                &
                      l_enable_diag_um,                                        &
                      l_environ_jo2,                                           &
                      l_environ_jo2b,                                          &
                      l_environ_ztop,                                          &
                      l_strat_chem,                                            &
                      fastjx_prescutoff,                                       &
                      timestep,                                                &
                      atau,                                                    &
                      atau0,                                                   &
                      fl,                                                      &
                      q1d,                                                     &
                      qo2,                                                     &
                      qo3,                                                     &
                      qqq,                                                     &
                      qrayl,                                                   &
                      tqq,                                                     &
                      wl,                                                      &
                      jfacta,                                                  &
                      daa,                                                     &
                      paa,                                                     &
                      qaa,                                                     &
                      raa,                                                     &
                      saa,                                                     &
                      waa,                                                     &
                      solcyc_av,                                               &
                      solcyc_quanta,                                           &
                      solcyc_ts,                                               &
                      solcyc_spec,                                             &
                      jlabel,                                                  &
                      titlej,                                                  &
                      pi,                                                      &
                      o3_mmr_vmr,                                              &
                      molemass_sulp,                                           &
                      molemass_nh42so4,                                        &
                      molemass_air,                                            &
                      planet_radius,                                           &
                      error_message, error_routine)

! ----------------------------------------------------------------------
! Description:
!
!  Given the input configuration control data, check its validity and
!  set up the Photolysis internal configuration data accordingly.
!  This includes some basic initialisation and everything required to
!  establish details of the selected configuration that will determine
!  the required environmental drivers.
!
! Method:
!
!  1. Copy the configuration variables provided as keyword arguments into
!     component variables with matching names in the photol_config structures.
!     For certain variables, default values are set here for use if input
!     values are not provided.
!     ------------------------------------------------------------------
!     Note: Values of variables that are inactive in the current
!     configuration will normally be ignored and not set.
!     ------------------------------------------------------------------
!  2. Check that configuration values are consistent.
!  3. FUTURE: Initialise photolysis rate names.
!  4. FUTURE: Set up lists of environmental driver fields required.
!  5. FUTURE: Initialise master diagnostics list and determine availability of
!     each diagnostic given the configuration.
!
! ----------------------------------------------------------------------

USE photol_config_specification_mod, ONLY: photol_config,                      &
       l_photol_config_available, i_scheme_nophot, i_scheme_photol_strat,      &
       i_scheme_phot2d, i_scheme_fastjx, i_obs_solcylc, i_avg_solcylc,         &
       init_photol_configuration, copy_config_value

USE photol_constants_mod,  ONLY: const_pi, const_pi_over_180,                  &
                              const_recip_pi_over_180,                         &
                              const_o3_mmr_vmr, const_molemass_sulp,           &
                              const_molemass_nh42so4, const_molemass_air,      &
                              const_planet_radius, const_s2r

USE photol_environment_mod, ONLY: photol_init_environ_req

USE photol_fieldname_mod,  ONLY: photol_jlabel_len
USE fastjx_data,           ONLY: fastjx_set_data_from_config

USE ukca_error_mod,        ONLY: maxlen_message, maxlen_procname,              &
                                 errcode_value_unknown, errcode_value_invalid, &
                                 errcode_value_missing, error_report
USE umPrintMgr,            ONLY: umMessage, umPrint

USE parkind1,               ONLY: jpim, jprb      ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Subroutine arguments.

! Except for the top-level photolysis scheme choice, each input configuration
! variable is an optional keyword argument with a matching component in the
! Photolysis configuration structure.
! Based on the scheme choices, certain variables may be expected to be
!  always defined by the parent as they are accessed by default in the
!  current photolysis workflow and not setting these can lead to unexpected
!  behaviour.

INTEGER, INTENT(IN) :: i_photol_scheme
INTEGER, TARGET, INTENT(OUT) :: error_code

! Optional arguments
INTEGER, OPTIONAL, INTENT(IN) :: chem_timestep

INTEGER, OPTIONAL, INTENT(IN) :: fastjx_mode
INTEGER, OPTIONAL, INTENT(IN) :: fastjx_numwl
INTEGER, OPTIONAL, INTENT(IN) :: global_row_length

INTEGER, OPTIONAL, INTENT(IN) :: i_solcylc_type

INTEGER, OPTIONAL, INTENT(IN) :: ip_aerosol_param_moist
INTEGER, OPTIONAL, INTENT(IN) :: ip_accum_sulphate
INTEGER, OPTIONAL, INTENT(IN) :: ip_aitken_sulphate

INTEGER, OPTIONAL, INTENT(IN) :: model_levels
INTEGER, OPTIONAL, INTENT(IN) :: n_cca_lev
INTEGER, OPTIONAL, INTENT(IN) :: solcylc_start_year
INTEGER, OPTIONAL, INTENT(IN) :: i_error_method

INTEGER, OPTIONAL, INTENT(IN) :: n_phot_spc
INTEGER, OPTIONAL, INTENT(IN) :: njval
INTEGER, OPTIONAL, INTENT(IN) :: nw1, nw2
INTEGER, OPTIONAL, INTENT(IN) :: jtaumx
INTEGER, OPTIONAL, INTENT(IN) :: naa
INTEGER, OPTIONAL, INTENT(IN) :: n_solcyc_ts
INTEGER, ALLOCATABLE, OPTIONAL, INTENT(IN) :: jind(:)

LOGICAL, OPTIONAL, INTENT(IN) :: l_cal360

LOGICAL, OPTIONAL, INTENT(IN) :: l_cloud_pc2
LOGICAL, OPTIONAL, INTENT(IN) :: l_3d_cca
LOGICAL, OPTIONAL, INTENT(IN) :: l_enable_diag_um

LOGICAL, OPTIONAL, INTENT(IN) :: l_environ_jo2
LOGICAL, OPTIONAL, INTENT(IN) :: l_environ_jo2b
LOGICAL, OPTIONAL, INTENT(IN) :: l_environ_ztop
LOGICAL, OPTIONAL, INTENT(IN) :: l_strat_chem

REAL, OPTIONAL, INTENT(IN)    :: fastjx_prescutoff
REAL, OPTIONAL, INTENT(IN)    :: timestep

! Fast-JX spectral file data
REAL, OPTIONAL, INTENT(IN) :: atau
REAL, OPTIONAL, INTENT(IN) :: atau0
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: fl(:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: q1d(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: qo2(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: qo3(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: qqq(:,:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: qrayl(:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: tqq(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: wl(:)

REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: jfacta(:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: daa(:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: paa(:,:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: qaa(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: raa(:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: saa(:,:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: waa(:,:)

REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: solcyc_av(:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: solcyc_quanta(:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: solcyc_ts(:)
REAL, ALLOCATABLE, OPTIONAL, INTENT(IN) :: solcyc_spec(:)

CHARACTER(LEN=photol_jlabel_len), ALLOCATABLE, OPTIONAL, INTENT(IN) :: jlabel(:)
CHARACTER(LEN=photol_jlabel_len), ALLOCATABLE, OPTIONAL, INTENT(IN) :: titlej(:)

! Configurable constants
REAL, OPTIONAL, INTENT(IN)    :: pi
REAL, OPTIONAL, INTENT(IN)    :: o3_mmr_vmr
REAL, OPTIONAL, INTENT(IN)    :: molemass_sulp
REAL, OPTIONAL, INTENT(IN)    :: molemass_nh42so4
REAL, OPTIONAL, INTENT(IN)    :: molemass_air
REAL, OPTIONAL, INTENT(IN)    :: planet_radius

CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables
INTEGER, POINTER :: error_code_ptr
CHARACTER(LEN=maxlen_message) :: err_message

CHARACTER(LEN=maxlen_message) :: var_missing
LOGICAL :: l_missing
INTEGER :: nvar_missing

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHOTOL_SETUP'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Set defaults for output arguments
error_code_ptr => error_code
error_code_ptr = 0
err_message = ''
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Set all configuration data to default values
CALL init_photol_configuration()

! First, check if parent has specified a method for error handling, in case of
! any errors further on.
IF (PRESENT(i_error_method)) photol_config%i_error_method = i_error_method

! Check that a known photolysis scheme is specified.
IF ( i_photol_scheme /= i_scheme_nophot        .AND.                           &
     i_photol_scheme /= i_scheme_photol_strat  .AND.                           &
     i_photol_scheme /= i_scheme_phot2d        .AND.                           &
     i_photol_scheme /= i_scheme_fastjx ) THEN
  error_code_ptr = errcode_value_unknown
  WRITE(err_message, '(A,I0)') 'Unknown Photolysis scheme specified ',         &
    i_photol_scheme
  CALL error_report(photol_config%i_error_method, error_code_ptr, err_message, &
         RoutineName, msg_out=error_message, locn_out=error_routine)

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
ELSE
  photol_config%i_photol_scheme = i_photol_scheme
END IF

! Collate input data specifying the UKCA configuration
! (loosely following 'photol_config_spec_type'

IF (PRESENT(chem_timestep))  photol_config%chem_timestep = chem_timestep
IF (PRESENT(fastjx_mode))    photol_config%fastjx_mode  = fastjx_mode
IF (PRESENT(fastjx_numwl))   photol_config%fastjx_numwl = fastjx_numwl

IF (PRESENT(global_row_length)) photol_config%global_row_length                &
                                                        = global_row_length
IF (PRESENT(i_solcylc_type)) photol_config%i_solcylc_type = i_solcylc_type

IF (PRESENT(model_levels)) photol_config%model_levels = model_levels
IF (PRESENT(n_cca_lev))    photol_config%n_cca_lev    = n_cca_lev

IF (PRESENT(solcylc_start_year)) photol_config%solcylc_start_year              &
                                                      = solcylc_start_year

IF (PRESENT(l_cal360))     photol_config%l_cal360     = l_cal360

IF (PRESENT(l_cloud_pc2))    photol_config%l_cloud_pc2 = l_cloud_pc2
IF (PRESENT(l_3d_cca))       photol_config%l_3d_cca    = l_3d_cca
IF (PRESENT(l_environ_jo2))  photol_config%l_environ_jo2  = l_environ_jo2
IF (PRESENT(l_environ_jo2b)) photol_config%l_environ_jo2b = l_environ_jo2b

IF (PRESENT(l_environ_ztop)) photol_config%l_environ_ztop = l_environ_ztop

IF (PRESENT(l_enable_diag_um))  photol_config%l_enable_diag_um                 &
                                                      = l_enable_diag_um
IF (PRESENT(l_strat_chem))   photol_config%l_strat_chem   = l_strat_chem

IF (PRESENT(fastjx_prescutoff)) photol_config%fastjx_prescutoff                &
                                                      = fastjx_prescutoff
IF (PRESENT(timestep))       photol_config%timestep   = timestep

IF (PRESENT(n_phot_spc)) photol_config%n_phot_spc = n_phot_spc

! Transfer spectral and solar cycle data if using Fast-JX scheme.
! Variables that are vectors or arrays are copied using a bespoke function.
IF ( i_photol_scheme == i_scheme_fastjx ) THEN

  ! Check that all the spectral / solar cycle variables have been provided.
  var_missing = 'FastJX variable missing:'
  nvar_missing = 0
  IF (.NOT. PRESENT(njval)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" njval"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(nw1)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" nw1"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(nw2)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" nw2"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(jtaumx)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" jtaumx"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(naa)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" naa"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(jind)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" jind"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(atau)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" atau"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(atau0)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" atau0"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(fl)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" fl"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(q1d)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" q1d"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(qo2)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" qo2"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(qo3)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" qo3"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(qqq)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" qqq"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(qrayl)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" qrayl"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(tqq)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" tqq"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(wl)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" wl"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(jfacta)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" jfacta"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(daa)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" daa"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(paa)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" paa"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(qaa)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" qaa"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(raa)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" raa"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(saa)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" saa"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(waa)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" waa"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(jlabel)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" jlabel"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (.NOT. PRESENT(titlej)) THEN
    nvar_missing = nvar_missing + 1
    WRITE(umMessage,'(A35)') TRIM(var_missing)//" titlej"
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  IF (nvar_missing>0) THEN
    error_code = errcode_value_missing
    WRITE(err_message, '(I3,A)') nvar_missing, ' required variable/s ' //      &
                                     'missing from spectral file data'
    CALL error_report(photol_config%i_error_method, error_code_ptr,            &
                      err_message, RoutineName, msg_out=error_message,         &
                      locn_out=error_routine)
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                            zhook_handle)
    RETURN
  END IF
  ! Solar cylc data - if scheme chosen
  l_missing = .FALSE.
  var_missing = ''
  IF ( photol_config%i_solcylc_type > 0 ) THEN
    IF (.NOT. PRESENT(n_solcyc_ts)) THEN
      l_missing = .TRUE.
      var_missing = TRIM(var_missing)//",n_solcyc_ts"
    END IF
    IF (.NOT. PRESENT(solcyc_av)) THEN
      l_missing = .TRUE.
      var_missing = TRIM(var_missing)//",solcyc_av"
    END IF
    IF (.NOT. PRESENT(solcyc_quanta)) THEN
      l_missing = .TRUE.
      var_missing = TRIM(var_missing)//",solcyc_quanta"
    END IF
    IF (.NOT. PRESENT(solcyc_ts)) THEN
      l_missing = .TRUE.
      var_missing = TRIM(var_missing)//",solcyc_ts"
    END IF
    IF (.NOT. PRESENT(solcyc_spec)) THEN
      l_missing = .TRUE.
      var_missing = TRIM(var_missing)//",solcyc_spec"
    END IF
  END IF  ! solcyc variables
  IF (l_missing) THEN
    error_code = errcode_value_missing
    WRITE(err_message, '(A,A)') 'Required solar cycle data missing: ',         &
      TRIM(var_missing)
    CALL error_report(photol_config%i_error_method, error_code_ptr,            &
                      err_message, RoutineName, msg_out=error_message,         &
                      locn_out=error_routine)
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                            zhook_handle)
    RETURN
  END IF
  photol_config%njval = njval
  photol_config%nw1 = nw1
  photol_config%nw2 = nw2
  photol_config%naa = naa
  photol_config%jtaumx = jtaumx
  photol_config%atau = atau
  photol_config%atau0 = atau0

  CALL copy_config_value(fl, photol_config%fl)
  CALL copy_config_value(q1d, photol_config%q1d)
  CALL copy_config_value(qo2, photol_config%qo2)
  CALL copy_config_value(qo3, photol_config%qo3)
  CALL copy_config_value(qqq, photol_config%qqq)
  CALL copy_config_value(qrayl, photol_config%qrayl)
  CALL copy_config_value(tqq, photol_config%tqq)
  CALL copy_config_value(wl, photol_config%wl)
  CALL copy_config_value(daa, photol_config%daa)
  CALL copy_config_value(paa, photol_config%paa)
  CALL copy_config_value(qaa, photol_config%qaa)
  CALL copy_config_value(raa, photol_config%raa)
  CALL copy_config_value(saa, photol_config%saa)
  CALL copy_config_value(waa, photol_config%waa)

  CALL copy_config_value(jind, photol_config%jind)
  CALL copy_config_value(jfacta, photol_config%jfacta)
  CALL copy_config_value(jlabel, photol_config%jlabel)
  CALL copy_config_value(titlej, photol_config%titlej)

  ! Transfer Solar Cycle data if to be read from file
  IF ( photol_config%i_solcylc_type > 0 ) THEN
    photol_config%n_solcyc_ts = n_solcyc_ts
    CALL copy_config_value(solcyc_av, photol_config%solcyc_av)
    CALL copy_config_value(solcyc_quanta, photol_config%solcyc_quanta)
    CALL copy_config_value(solcyc_ts, photol_config%solcyc_ts)
    CALL copy_config_value(solcyc_spec, photol_config%solcyc_spec)
  END IF

  ! Call routine to transfer spectral and solar data to internal variables
  CALL fastjx_set_data_from_config()

END IF   ! Fast-JX

! Copy constant variables if passed in
IF (PRESENT(pi)) THEN
  const_pi = pi
  const_pi_over_180      = const_pi / 180.0
  const_recip_pi_over_180 = 180.0 / const_pi
  const_s2r             = (2.0*const_pi)/ 86400.0
END IF

IF (PRESENT(o3_mmr_vmr)) const_o3_mmr_vmr = o3_mmr_vmr
IF (PRESENT(molemass_sulp)) const_molemass_sulp = molemass_sulp
IF (PRESENT(molemass_nh42so4)) const_molemass_nh42so4 = molemass_nh42so4
IF (PRESENT(molemass_air)) const_molemass_air = molemass_air
IF (PRESENT(planet_radius)) const_planet_radius = planet_radius

! Set flag to show that a valid photolysis configuration is set up
l_photol_config_available = .TRUE.

! Routine to set up list of driving fields required based on user choices
CALL photol_init_environ_req(error_code_ptr, error_message=error_message,      &
                             error_routine=error_routine)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE photol_setup

END MODULE photol_setup_mod
