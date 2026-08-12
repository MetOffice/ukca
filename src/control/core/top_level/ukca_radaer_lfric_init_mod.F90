! -----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! -----------------------------------------------------------------------------
!
! Purpose: An interface routine to initialise radaer setting with lfric parent
!
! Procedure:
!   1) CALL the relevant mode setup subroutine from ukca_mode_setup
!
!   2) CALL 
!
! This file belongs in section: UKCA
!
MODULE ukca_radaer_lfric_init_mod

IMPLICIT NONE

! Default private
PRIVATE

PUBLIC :: ukca_radaer_lfric_init

PUBLIC :: n_ukca_mode , n_ukca_cpnt

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName = 'UKCA_RADAER_LFRIC_INIT_MOD'

CONTAINS

SUBROUTINE ukca_radaer_lfric_init( i_mode_setup_in,                            &
                                   i_tune_bc_in,                               &
                                   l_dust_mp_ageing,                           &
                                   l_ukca_radaer_sustrat )

USE ereport_mod,                       ONLY:                                   &
    ereport

USE errormessagelength_mod,            ONLY:                                   &
    errormessagelength

USE parkind1,                          ONLY:                                   &
    jpim,                                                                      &
    jprb

USE ukca_config_specification_mod,     ONLY:                                   &
    i_suss_4mode,                                                              &
    i_sussbcoc_5mode,                                                          &
    i_sussbcoc_4mode,                                                          &
    i_sussbcocso_5mode,                                                        &
    i_sussbcocso_4mode,                                                        &
    i_du_2mode,                                                                &
    i_sussbcocdu_7mode,                                                        &
    i_sussbcocntnh_5mode_7cpt,                                                 &
    i_solinsol_6mode,                                                          &
    i_sussbcocduntnh_8mode_8cpt,                                               &
    i_sussbcocdump_8mode
    glomap_variables_radaer

USE ukca_mode_setup,                   ONLY:                                   &
    ukca_mode_suss_4mode,                                                      &
    ukca_mode_sussbcoc_5mode,                                                  &
    ukca_mode_sussbcoc_4mode,                                                  &
    ukca_mode_sussbcocso_5mode,                                                &
    ukca_mode_sussbcocso_4mode,                                                &
    ukca_mode_duonly_2mode,                                                    &
    ukca_mode_sussbcocdu_7mode,                                                &
    ukca_mode_sussbcocntnh_5mode_7cpt,                                         &
    ukca_mode_solinsol_6mode,                                                  &
    ukca_mode_sussbcocduntnh_8mode_8cpt,                                       &
    ukca_mode_sussbcocdump_8mode,                                              &
    glomap_variables_type

USE ukca_radaer_lfric_settings_mod,    ONLY:                                   &
    ukca_radaer_lfric_settings

USE ukca_radaer_lfric_saved_mod,       ONLY:                                   &
    ukca_radaer_lfric

USE umPrintMgr,                        ONLY:                                   &
    umPrint,                                                                   &
    umMessage

USE yomhook,                           ONLY:                                   &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: i_mode_setup_in
INTEGER, INTENT(IN) :: i_tune_bc_in
LOGICAL, INTENT(IN) :: l_dust_mp_ageing

! Local variables

! These are used in radaer_kernel_mod
INTEGER, SAVE :: n_ukca_mode = 0
INTEGER, SAVE :: n_ukca_cpnt = 0

! We are calling this module because radaer is on. Forced true
LOGICAL, PARAMETER :: l_radaer_in = .true.

! Force existing temporary logicals to always be true from LFRic
LOGICAL, PARAMETER :: l_fix_nacl_density_in = .true.
LOGICAL, PARAMETER :: l_fix_ukca_hygroscopicities_in = .true.

INTEGER                           :: errcode  ! error code
CHARACTER(LEN=errormessagelength) :: cmessage ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*),   PARAMETER :: RoutineName='UKCA_RADAER_LFRIC_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

SELECT CASE(i_mode_setup_in)
CASE (i_suss_4mode) ! 1
  CALL ukca_mode_suss_4mode( glomap_variables_radaer,                          &
                             l_radaer_in,                                      &
                             i_tune_bc_in,                                     &
                             l_fix_nacl_density_in,                            &
                             l_fix_ukca_hygroscopicities_in,                   &
                             l_dust_mp_ageing)

CASE (i_sussbcoc_5mode) ! 2
  CALL ukca_mode_sussbcoc_5mode( glomap_variables_radaer,                      &
                                 l_radaer_in,                                  &
                                 i_tune_bc_in,                                 &
                                 l_fix_nacl_density_in,                        &
                                 l_fix_ukca_hygroscopicities_in,               &
                                 l_dust_mp_ageing)

CASE (i_sussbcoc_4mode) ! 3
  CALL ukca_mode_sussbcoc_4mode( glomap_variables_radaer,                      &
                                 l_radaer_in,                                  &
                                 i_tune_bc_in,                                 &
                                 l_fix_nacl_density_in,                        &
                                 l_fix_ukca_hygroscopicities_in,               &
                                 l_dust_mp_ageing )

CASE (i_sussbcocso_5mode) ! 4
  CALL ukca_mode_sussbcocso_5mode( glomap_variables_radaer,                    &
                                   l_radaer_in,                                &
                                   i_tune_bc_in,                               &
                                   l_fix_nacl_density_in,                      &
                                   l_fix_ukca_hygroscopicities_in,             &
                                   l_dust_mp_ageing )

CASE (i_sussbcocso_4mode) ! 5
  CALL ukca_mode_sussbcocso_4mode( glomap_variables_radaer,                    &
                                   l_radaer_in,                                &
                                   i_tune_bc_in,                               &
                                   l_fix_nacl_density_in,                      &
                                   l_fix_ukca_hygroscopicities_in,             &
                                   l_dust_mp_ageing )

CASE (i_du_2mode) ! 6
  CALL ukca_mode_duonly_2mode( glomap_variables_radaer,                        &
                               l_radaer_in,                                    &
                               i_tune_bc_in,                                   &
                               l_fix_nacl_density_in,                          &
                               l_fix_ukca_hygroscopicities_in,                 &
                               l_dust_mp_ageing )

CASE (i_sussbcocdu_7mode) ! 8
  CALL ukca_mode_sussbcocdu_7mode( glomap_variables_radaer,                    &
                                   l_radaer_in,                                &
                                   i_tune_bc_in,                               &
                                   l_fix_nacl_density_in,                      &
                                   l_fix_ukca_hygroscopicities_in,             &
                                   l_dust_mp_ageing )

CASE (i_sussbcocntnh_5mode_7cpt) ! 10
  CALL ukca_mode_sussbcocntnh_5mode_7cpt( glomap_variables_radaer,             &
                                          l_radaer_in,                         &
                                          i_tune_bc_in,                        &
                                          l_fix_nacl_density_in,               &
                                          l_fix_ukca_hygroscopicities_in,      &
                                          l_dust_mp_ageing )

CASE (i_solinsol_6mode) ! 11
  CALL ukca_mode_solinsol_6mode( glomap_variables_radaer,                      &
                                 l_radaer_in,                                  &
                                 i_tune_bc_in,                                 &
                                 l_fix_nacl_density_in,                        &
                                 l_fix_ukca_hygroscopicities_in,               &
                                 l_dust_mp_ageing )

CASE (i_sussbcocduntnh_8mode_8cpt) ! 12
  CALL ukca_mode_sussbcocduntnh_8mode_8cpt( glomap_variables_radaer,           &
                                            l_radaer_in,                       &
                                            i_tune_bc_in,                      &
                                            l_fix_nacl_density_in,             &
                                            l_fix_ukca_hygroscopicities_in,    &
                                            l_dust_mp_ageing )

CASE (i_sussbcocdump_8mode) ! 13
  CALL ukca_mode_sussbcocdump_8mode( glomap_variables_radaer,                  &
                                     l_radaer_in,                              &
                                     i_tune_bc_in,                             &
                                     l_fix_nacl_density_in,                    &
                                     l_fix_ukca_hygroscopicities_in,           &
                                     l_dust_mp_ageing )

CASE DEFAULT
  cmessage='i_mode_setup_in has unrecognised value'
  WRITE(umMessage,'(A,I0)') cmessage, i_mode_setup_in
  CALL umPrint(umMessage,src=RoutineName)
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END SELECT

! Now call ukca_radaer_lfric_settings
CALL ukca_radaer_lfric_settings ( l_ukca_radaer_sustrat,                       &
                                  glomap_variables_radaer,                     &
                                  ukca_radaer_lfric,                           &
                                  n_ukca_mode,                                 &
                                  n_ukca_cpnt )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE ukca_radaer_lfric_init

END MODULE ukca_radaer_lfric_init_mod
