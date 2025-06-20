! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Interface between ukca_mode_setup.F90 and external repository
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
!   1) CALL common_mode_setup_interface
!
! ---------------------------------------------------------------------


MODULE glomap_clim_mode_setup_interface_mod

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE ::                                          &
ModuleName='GLOMAP_CLIM_MODE_SETUP_INTERFACE_MOD'

CONTAINS

SUBROUTINE glomap_clim_mode_setup_interface ( i_mode_setup_in,                 &
                                              l_radaer_in,                     &
                                              i_tune_bc_in,                    &
                                              l_fix_nacl_density_in,           &
                                              l_fix_ukca_hygroscopicities_in,  &
                                              l_dust_mp_ageing )

USE common_mode_setup_interface_mod, ONLY:                                     &
    common_mode_setup_interface

USE parkind1,                        ONLY:                                     &
    jpim,                                                                      &
    jprb

USE ukca_config_specification_mod,   ONLY:                                     &
    glomap_variables_climatology

USE ukca_config_constants_mod,        ONLY:                                    &
    l_ukca_constants_available

USE ereport_mod,                     ONLY:                                     &
    ereport

USE errormessagelength_mod,          ONLY:                                     &
    errormessagelength

USE yomhook,                         ONLY:                                     &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: i_mode_setup_in
LOGICAL, INTENT(IN) :: l_radaer_in
INTEGER, INTENT(IN) :: i_tune_bc_in
LOGICAL, INTENT(IN) :: l_fix_nacl_density_in
LOGICAL, INTENT(IN) :: l_fix_ukca_hygroscopicities_in
LOGICAL, INTENT(IN) :: l_dust_mp_ageing

! Local variables

CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER                           :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='GLOMAP_CLIM_MODE_SETUP_INTERFACE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! This routine is called from outside UKCA, by-passing the UKCA API.
! This is only allowed if UKCA constants are properly set up.
IF (.NOT. l_ukca_constants_available) THEN
  cmessage = 'Configurable UKCA constants have not been set up'
  errcode = 1
  CALL ereport(RoutineName, errcode, cmessage)
END IF

CALL common_mode_setup_interface ( glomap_variables_climatology,               &
                                   i_mode_setup_in,                            &
                                   l_radaer_in,                                &
                                   i_tune_bc_in,                               &
                                   l_fix_nacl_density_in,                      &
                                   l_fix_ukca_hygroscopicities_in,             &
                                   l_dust_mp_ageing )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE glomap_clim_mode_setup_interface

END MODULE glomap_clim_mode_setup_interface_mod
