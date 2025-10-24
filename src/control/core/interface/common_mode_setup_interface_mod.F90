! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: A wrapper routine to call the appropriate mode variables
!          configuration routine from ukca_mode_setup
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!
! This file belongs in section: INTERFACE
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! Procedure:
!   1) CALL the relevant mode setup subroutine from ukca_mode_setup
!
! ---------------------------------------------------------------------


MODULE common_mode_setup_interface_mod

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE ::                                          &
                                   ModuleName='COMMON_MODE_SETUP_INTERFACE_MOD'

CONTAINS

SUBROUTINE common_mode_setup_interface ( glomap_variables_local,               &
                                         i_mode_setup_in,                      &
                                         l_radaer_in,                          &
                                         i_tune_bc_in,                         &
                                         l_fix_nacl_density_in,                &
                                         l_fix_ukca_hygroscopicities_in,       &
                                         l_dust_mp_ageing)

USE ereport_mod,                     ONLY:                                     &
    ereport

USE errormessagelength_mod,          ONLY:                                     &
    errormessagelength

USE parkind1,                        ONLY:                                     &
    jpim,                                                                      &
    jprb

USE ukca_mode_setup,                 ONLY:                                     &
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

USE ukca_config_specification_mod,   ONLY:                                     &
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

USE umPrintMgr,                      ONLY:                                     &
    umPrint,                                                                   &
    umMessage

USE yomhook,                         ONLY:                                     &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
INTEGER, INTENT(IN) :: i_mode_setup_in
LOGICAL, INTENT(IN) :: l_radaer_in
INTEGER, INTENT(IN) :: i_tune_bc_in
LOGICAL, INTENT(IN) :: l_fix_nacl_density_in
LOGICAL, INTENT(IN) :: l_fix_ukca_hygroscopicities_in
LOGICAL, INTENT(IN) :: l_dust_mp_ageing

! Local variables

INTEGER                           :: errcode  ! error code
CHARACTER(LEN=errormessagelength) :: cmessage ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='COMMON_MODE_SETUP_INTERFACE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

SELECT CASE(i_mode_setup_in)
CASE (i_suss_4mode) ! 1
  CALL ukca_mode_suss_4mode( glomap_variables_local,                           &
                             l_radaer_in,                                      &
                             i_tune_bc_in,                                     &
                             l_fix_nacl_density_in,                            &
                             l_fix_ukca_hygroscopicities_in,                   &
                             l_dust_mp_ageing)

CASE (i_sussbcoc_5mode) ! 2
  CALL ukca_mode_sussbcoc_5mode( glomap_variables_local,                       &
                                 l_radaer_in,                                  &
                                 i_tune_bc_in,                                 &
                                 l_fix_nacl_density_in,                        &
                                 l_fix_ukca_hygroscopicities_in,               &
                                 l_dust_mp_ageing)

CASE (i_sussbcoc_4mode) ! 3
  CALL ukca_mode_sussbcoc_4mode( glomap_variables_local,                       &
                                 l_radaer_in,                                  &
                                 i_tune_bc_in,                                 &
                                 l_fix_nacl_density_in,                        &
                                 l_fix_ukca_hygroscopicities_in,               &
                                 l_dust_mp_ageing )

CASE (i_sussbcocso_5mode) ! 4
  CALL ukca_mode_sussbcocso_5mode( glomap_variables_local,                     &
                                   l_radaer_in,                                &
                                   i_tune_bc_in,                               &
                                   l_fix_nacl_density_in,                      &
                                   l_fix_ukca_hygroscopicities_in,             &
                                   l_dust_mp_ageing )

CASE (i_sussbcocso_4mode) ! 5
  CALL ukca_mode_sussbcocso_4mode( glomap_variables_local,                     &
                                   l_radaer_in,                                &
                                   i_tune_bc_in,                               &
                                   l_fix_nacl_density_in,                      &
                                   l_fix_ukca_hygroscopicities_in,             &
                                   l_dust_mp_ageing )

CASE (i_du_2mode) ! 6
  CALL ukca_mode_duonly_2mode( glomap_variables_local,                         &
                               l_radaer_in,                                    &
                               i_tune_bc_in,                                   &
                               l_fix_nacl_density_in,                          &
                               l_fix_ukca_hygroscopicities_in,                 &
                               l_dust_mp_ageing )

CASE (i_sussbcocdu_7mode) ! 8
  CALL ukca_mode_sussbcocdu_7mode( glomap_variables_local,                     &
                                   l_radaer_in,                                &
                                   i_tune_bc_in,                               &
                                   l_fix_nacl_density_in,                      &
                                   l_fix_ukca_hygroscopicities_in,             &
                                   l_dust_mp_ageing )

CASE (i_sussbcocntnh_5mode_7cpt) ! 10
  CALL ukca_mode_sussbcocntnh_5mode_7cpt( glomap_variables_local,              &
                                          l_radaer_in,                         &
                                          i_tune_bc_in,                        &
                                          l_fix_nacl_density_in,               &
                                          l_fix_ukca_hygroscopicities_in,      &
                                          l_dust_mp_ageing )

CASE (i_solinsol_6mode) ! 11
  CALL ukca_mode_solinsol_6mode( glomap_variables_local,                       &
                                 l_radaer_in,                                  &
                                 i_tune_bc_in,                                 &
                                 l_fix_nacl_density_in,                        &
                                 l_fix_ukca_hygroscopicities_in,               &
                                 l_dust_mp_ageing )

CASE (i_sussbcocduntnh_8mode_8cpt) ! 12
  CALL ukca_mode_sussbcocduntnh_8mode_8cpt( glomap_variables_local,            &
                                            l_radaer_in,                       &
                                            i_tune_bc_in,                      &
                                            l_fix_nacl_density_in,             &
                                            l_fix_ukca_hygroscopicities_in,    &
                                            l_dust_mp_ageing )

CASE (i_sussbcocdump_8mode) ! 13
  CALL ukca_mode_sussbcocdump_8mode( glomap_variables_local,                   &
                                            l_radaer_in,                       &
                                            i_tune_bc_in,                      &
                                            l_fix_nacl_density_in,             &
                                            l_fix_ukca_hygroscopicities_in,    &
                                            l_dust_mp_ageing )

CASE DEFAULT
  cmessage='i_mode_setup_in has unrecognised value'
  WRITE(umMessage,'(A,I0)') cmessage, i_mode_setup_in
  CALL umPrint(umMessage,src=RoutineName)
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE common_mode_setup_interface

END MODULE common_mode_setup_interface_mod
