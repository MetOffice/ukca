! -----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! -----------------------------------------------------------------------------
!
! Purpose:
!   Allocatable lists used in radaer kernel are set here
!
! ---------------------------------------------------------------------

MODULE ukca_radaer_lfric_list_mod

IMPLICIT NONE
PRIVATE

! Maximum number of characters allowed for names in list
INTEGER, PARAMETER, PUBLIC :: max_fldname_len = 40

CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: pvol_comp_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: comp_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: mode_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: rhopar_mode_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: dry_diam_mode_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: modal_volume_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: ait_sol_volume_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: acc_sol_volume_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: cor_sol_volume_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: ait_ins_volume_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: acc_ins_volume_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: cor_ins_volume_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: modal_wtv_names(:)
CHARACTER(LEN=max_fldname_len),PUBLIC,ALLOCATABLE,SAVE:: wet_diam_mode_names(:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_RADAER_LFRIC_LIST_MOD'

PUBLIC :: ukca_radaer_lfric_list

CONTAINS

SUBROUTINE ukca_radaer_lfric_list( i_mode_setup_in )

USE ereport_mod,                       ONLY:                                   &
    ereport

USE errormessagelength_mod,            ONLY:                                   &
     errormessagelength

USE parkind1,                          ONLY:                                   &
    jpim,                                                                      &
    jprb

USE ukca_config_specification_mod,     ONLY:                                   &
    i_du_2mode,                                                                &
    i_sussbcocdu_7mode

USE umPrintMgr,                        ONLY:                                   &
    umPrint,                                                                   &
    umMessage

USE yomhook,                           ONLY:                                   &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: i_mode_setup_in

! Local variables

INTEGER                           :: errcode  ! error code
CHARACTER(LEN=errormessagelength) :: cmessage ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*),   PARAMETER :: RoutineName='UKCA_RADAER_LFRIC_LIST'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

SELECT CASE(i_mode_setup_in)
CASE (i_du_2mode)

  IF (.NOT. ALLOCATED(pvol_comp_names))      ALLOCATE(pvol_comp_names(2))
  IF (.NOT. ALLOCATED(comp_names))           ALLOCATE(comp_names(2))
  IF (.NOT. ALLOCATED(mode_names))           ALLOCATE(mode_names(6))
  IF (.NOT. ALLOCATED(rhopar_mode_names))    ALLOCATE(rhopar_mode_names(6))
  IF (.NOT. ALLOCATED(dry_diam_mode_names))  ALLOCATE(dry_diam_mode_names(6))
  IF (.NOT. ALLOCATED(modal_volume_names))   ALLOCATE(modal_volume_names(5))
  IF (.NOT. ALLOCATED(ait_sol_volume_names)) ALLOCATE(ait_sol_volume_names(1))
  IF (.NOT. ALLOCATED(acc_sol_volume_names)) ALLOCATE(acc_sol_volume_names(1))
  IF (.NOT. ALLOCATED(cor_sol_volume_names)) ALLOCATE(cor_sol_volume_names(1))
  IF (.NOT. ALLOCATED(ait_ins_volume_names)) ALLOCATE(ait_ins_volume_names(1))
  IF (.NOT. ALLOCATED(acc_ins_volume_names)) ALLOCATE(acc_ins_volume_names(1))
  IF (.NOT. ALLOCATED(cor_ins_volume_names)) ALLOCATE(cor_ins_volume_names(1))
  IF (.NOT. ALLOCATED(modal_wtv_names))      ALLOCATE(modal_wtv_names(6))
  IF (.NOT. ALLOCATED(wet_diam_mode_names))  ALLOCATE(wet_diam_mode_names(6))

  pvol_comp_names = [ 'fldname_pvol_du_acc_ins' ,                              &
                      'fldname_pvol_du_cor_ins' ]

  comp_names = [ 'fldname_acc_ins_du' ,                                        &
                 'fldname_cor_ins_du' ]

  mode_names = [ 'fldname_n_ait_sol' ,                                         &
                 'fldname_n_acc_sol' ,                                         &
                 'fldname_n_cor_sol' ,                                         &
                 'fldname_n_ait_ins' ,                                         &
                 'fldname_n_acc_ins' ,                                         &
                 'fldname_n_cor_ins' ]

  rhopar_mode_names = [ 'fldname_rhopar_ait_sol' ,                             &
                        'fldname_rhopar_acc_sol' ,                             &
                        'fldname_rhopar_cor_sol' ,                             &
                        'fldname_rhopar_ait_ins' ,                             &
                        'fldname_rhopar_acc_ins' ,                             &
                        'fldname_rhopar_cor_ins' ]

  dry_diam_mode_names = [ 'fldname_drydp_ait_sol' ,                            &
                          'fldname_drydp_acc_sol' ,                            &
                          'fldname_drydp_cor_sol' ,                            &
                          'fldname_drydp_ait_ins' ,                            &
                          'fldname_drydp_acc_ins' ,                            &
                          'fldname_drydp_cor_ins' ]

  modal_volume_names = [ 'fldname_mod_vol_ait_sol' ,                           &
                         'fldname_mod_vol_acc_sol' ,                           &
                         'fldname_mod_vol_cor_sol' ,                           &
                         'fldname_mod_vol_acc_ins' ,                           &
                         'fldname_mod_vol_cor_ins' ]

  ait_sol_volume_names = [ 'pvol_wat_ait_sol' ]

  acc_sol_volume_names = [ 'pvol_wat_acc_sol' ]

  cor_sol_volume_names = [ 'pvol_wat_cor_sol' ]

  ait_ins_volume_names = [ 'null' ]

  acc_ins_volume_names = [ 'pvol_du_acc_ins' ]

  cor_ins_volume_names = [ 'pvol_du_cor_ins' ]

  modal_wtv_names = [ 'fldname_pvol_wat_ait_sol' ,                             &
                      'fldname_pvol_wat_acc_sol' ,                             &
                      'fldname_pvol_wat_cor_sol' ,                             &
                      'fldname_pvol_wat_ait_ins' ,                             &
                      'fldname_pvol_wat_acc_ins' ,                             &
                      'fldname_pvol_wat_cor_ins' ]

  wet_diam_mode_names = [ 'fldname_wetdp_ait_sol' ,                            &
                          'fldname_wetdp_acc_sol' ,                            &
                          'fldname_wetdp_cor_sol' ,                            &
                          'fldname_wetdp_ait_ins' ,                            &
                          'fldname_wetdp_acc_ins' ,                            &
                          'fldname_wetdp_cor_ins' ]

CASE (i_sussbcocdu_7mode)

  IF (.NOT. ALLOCATED(pvol_comp_names))       ALLOCATE(pvol_comp_names(17))
  IF (.NOT. ALLOCATED(comp_names))            ALLOCATE(comp_names(17))
  IF (.NOT. ALLOCATED(mode_names))            ALLOCATE(mode_names(6))
  IF (.NOT. ALLOCATED(rhopar_mode_names))     ALLOCATE(rhopar_mode_names(6))
  IF (.NOT. ALLOCATED(dry_diam_mode_names))   ALLOCATE(dry_diam_mode_names(6))
  IF (.NOT. ALLOCATED(modal_volume_names))    ALLOCATE(modal_volume_names(6))
  IF (.NOT. ALLOCATED(ait_sol_volume_names))  ALLOCATE(ait_sol_volume_names(4))
  IF (.NOT. ALLOCATED(acc_sol_volume_names))  ALLOCATE(acc_sol_volume_names(5))
  IF (.NOT. ALLOCATED(cor_sol_volume_names))  ALLOCATE(cor_sol_volume_names(5))
  IF (.NOT. ALLOCATED(ait_ins_volume_names))  ALLOCATE(ait_ins_volume_names(2))
  IF (.NOT. ALLOCATED(acc_ins_volume_names))  ALLOCATE(acc_ins_volume_names(1))
  IF (.NOT. ALLOCATED(cor_ins_volume_names))  ALLOCATE(cor_ins_volume_names(1))
  IF (.NOT. ALLOCATED(modal_wtv_names))       ALLOCATE(modal_wtv_names(6))
  IF (.NOT. ALLOCATED(wet_diam_mode_names))   ALLOCATE(wet_diam_mode_names(6))

  pvol_comp_names = [ 'fldname_pvol_su_ait_sol' ,                              &
                      'fldname_pvol_bc_ait_sol' ,                              &
                      'fldname_pvol_om_ait_sol' ,                              &
                      'fldname_pvol_su_acc_sol' ,                              &
                      'fldname_pvol_bc_acc_sol' ,                              &
                      'fldname_pvol_om_acc_sol' ,                              &
                      'fldname_pvol_ss_acc_sol' ,                              &
                      'fldname_pvol_du_acc_sol' ,                              &
                      'fldname_pvol_su_cor_sol' ,                              &
                      'fldname_pvol_bc_cor_sol' ,                              &
                      'fldname_pvol_om_cor_sol' ,                              &
                      'fldname_pvol_ss_cor_sol' ,                              &
                      'fldname_pvol_du_cor_sol' ,                              &
                      'fldname_pvol_bc_ait_ins' ,                              &
                      'fldname_pvol_om_ait_ins' ,                              &
                      'fldname_pvol_du_acc_ins' ,                              &
                      'fldname_pvol_du_cor_ins' ]

  comp_names = [ 'fldname_ait_sol_su' ,                                        &
                 'fldname_ait_sol_bc' ,                                        &
                 'fldname_ait_sol_om' ,                                        &
                 'fldname_acc_sol_su' ,                                        &
                 'fldname_acc_sol_bc' ,                                        &
                 'fldname_acc_sol_om' ,                                        &
                 'fldname_acc_sol_ss' ,                                        &
                 'fldname_acc_sol_du' ,                                        &
                 'fldname_cor_sol_su' ,                                        &
                 'fldname_cor_sol_bc' ,                                        &
                 'fldname_cor_sol_om' ,                                        &
                 'fldname_cor_sol_ss' ,                                        &
                 'fldname_cor_sol_du' ,                                        &
                 'fldname_ait_ins_bc' ,                                        &
                 'fldname_ait_ins_om' ,                                        &
                 'fldname_acc_ins_du' ,                                        &
                 'fldname_cor_ins_du' ]

  mode_names = [ 'fldname_n_ait_sol' ,                                         &
                 'fldname_n_acc_sol' ,                                         &
                 'fldname_n_cor_sol' ,                                         &
                 'fldname_n_ait_ins' ,                                         &
                 'fldname_n_acc_ins' ,                                         &
                 'fldname_n_cor_ins' ]

  rhopar_mode_names = [ 'fldname_rhopar_ait_sol' ,                             &
                        'fldname_rhopar_acc_sol' ,                             &
                        'fldname_rhopar_cor_sol' ,                             &
                        'fldname_rhopar_ait_ins' ,                             &
                        'fldname_rhopar_acc_ins' ,                             &
                        'fldname_rhopar_cor_ins' ]

  dry_diam_mode_names = [ 'fldname_drydp_ait_sol' ,                            &
                          'fldname_drydp_acc_sol' ,                            &
                          'fldname_drydp_cor_sol' ,                            &
                          'fldname_drydp_ait_ins' ,                            &
                          'fldname_drydp_acc_ins' ,                            &
                          'fldname_drydp_cor_ins' ]

  modal_volume_names = [ 'fldname_mod_vol_ait_sol' ,                           &
                         'fldname_mod_vol_acc_sol' ,                           &
                         'fldname_mod_vol_cor_sol' ,                           &
                         'fldname_mod_vol_ait_ins' ,                           &
                         'fldname_mod_vol_acc_ins' ,                           &
                         'fldname_mod_vol_cor_ins' ]

  ait_sol_volume_names = [ 'pvol_wat_ait_sol' ,                                &
                           'pvol_su_ait_sol ' ,                                &
                           'pvol_bc_ait_sol ' ,                                &
                           'pvol_om_ait_sol ' ]

  acc_sol_volume_names = [ 'pvol_wat_acc_sol' ,                                &
                           'pvol_su_acc_sol ' ,                                &
                           'pvol_bc_acc_sol ' ,                                &
                           'pvol_om_acc_sol ' ,                                &
                           'pvol_ss_acc_sol ' ]

  cor_sol_volume_names = [ 'pvol_wat_cor_sol' ,                                &
                           'pvol_su_cor_sol ' ,                                &
                           'pvol_bc_cor_sol ' ,                                &
                           'pvol_om_cor_sol ' ,                                &
                           'pvol_ss_cor_sol ' ]

  ait_ins_volume_names = [ 'pvol_bc_ait_ins' ,                                 &
                           'pvol_om_ait_ins' ]

  acc_ins_volume_names = [ 'pvol_du_acc_ins' ]

  cor_ins_volume_names = [ 'pvol_du_cor_ins' ]

  modal_wtv_names = [ 'fldname_pvol_wat_ait_sol' ,                             &
                      'fldname_pvol_wat_acc_sol' ,                             &
                      'fldname_pvol_wat_cor_sol' ,                             &
                      'fldname_pvol_wat_ait_ins' ,                             &
                      'fldname_pvol_wat_acc_ins' ,                             &
                      'fldname_pvol_wat_cor_ins' ]

  wet_diam_mode_names = [ 'fldname_wetdp_ait_sol' ,                            &
                          'fldname_wetdp_acc_sol' ,                            &
                          'fldname_wetdp_cor_sol' ,                            &
                          'fldname_wetdp_ait_ins' ,                            &
                          'fldname_wetdp_acc_ins' ,                            &
                          'fldname_wetdp_cor_ins' ]

CASE DEFAULT
  cmessage='i_mode_setup_in has unrecognised value'
  WRITE(umMessage,'(A,I0)') cmessage, i_mode_setup_in
  CALL umPrint(umMessage,src=RoutineName)
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE ukca_radaer_lfric_list

END MODULE ukca_radaer_lfric_list_mod
