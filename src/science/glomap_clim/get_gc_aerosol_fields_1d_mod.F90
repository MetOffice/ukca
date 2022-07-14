! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution..
! *****************************COPYRIGHT*******************************
!
! Purpose:
!   To put aerosol fields into nmr1d and mmr1d
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE get_gc_aerosol_fields_1d_mod

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='GET_GC_AEROSOL_FIELDS_1D_MOD'

PUBLIC :: get_gc_aerosol_fields_1d

CONTAINS

SUBROUTINE get_gc_aerosol_fields_1d( n_points, i_glomap_clim_setup_in,         &
  gc_nd_nuc_sol_1d, gc_nuc_sol_su_1d, gc_nuc_sol_oc_1d,                        &
  gc_nd_ait_sol_1d, gc_ait_sol_su_1d, gc_ait_sol_bc_1d, gc_ait_sol_oc_1d,      &
  gc_nd_acc_sol_1d, gc_acc_sol_su_1d, gc_acc_sol_bc_1d, gc_acc_sol_oc_1d,      &
  gc_acc_sol_ss_1d, gc_acc_sol_du_1d,                                          &
  gc_nd_cor_sol_1d, gc_cor_sol_su_1d, gc_cor_sol_bc_1d, gc_cor_sol_oc_1d,      &
  gc_cor_sol_ss_1d, gc_cor_sol_du_1d,                                          &
  gc_nd_ait_ins_1d, gc_ait_ins_bc_1d, gc_ait_ins_oc_1d,                        &
  gc_nd_acc_ins_1d, gc_acc_ins_du_1d,                                          &
  gc_nd_cor_ins_1d, gc_cor_ins_du_1d,                                          &
  nmr1d, mmr1d )

USE ereport_mod,            ONLY:                                              &
    ereport

USE errormessagelength_mod, ONLY:                                              &
    errormessagelength

USE glomap_clim_option_mod, ONLY:                                              &
    i_gc_sussocbc_5mode,                                                       &
    i_gc_sussocbcdu_7mode

USE parkind1,               ONLY:                                              &
    jprb,                                                                      &
    jpim

USE ukca_mode_setup,        ONLY:                                              &
    cp_su,                                                                     &
    cp_bc,                                                                     &
    cp_oc,                                                                     &
    cp_cl,                                                                     &
    cp_du,                                                                     &
    mode_nuc_sol,                                                              &
    mode_ait_sol,                                                              &
    mode_acc_sol,                                                              &
    mode_cor_sol,                                                              &
    mode_ait_insol,                                                            &
    mode_acc_insol,                                                            &
    mode_cor_insol,                                                            &
    ncp,                                                                       &
    nmodes

USE umPrintMgr,             ONLY:                                              &
    newline

USE um_types,               ONLY:                                              &
    real_umphys

USE yomhook,                ONLY:                                              &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: n_points

INTEGER, INTENT(IN) :: i_glomap_clim_setup_in

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_nuc_sol_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nuc_sol_su_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nuc_sol_oc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_ait_sol_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_sol_su_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_sol_bc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_sol_oc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_acc_sol_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_su_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_bc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_oc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_ss_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_du_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_cor_sol_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_su_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_bc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_oc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_ss_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_du_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_ait_ins_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_ins_bc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_ins_oc_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_acc_ins_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_ins_du_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_cor_ins_1d( n_points )

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_ins_du_1d( n_points )

REAL(KIND=real_umphys),INTENT(OUT) :: nmr1d( n_points, nmodes )

REAL(KIND=real_umphys),INTENT(OUT) :: mmr1d( n_points, nmodes, ncp )

! Local Variables

INTEGER                           :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='GET_GC_AEROSOL_FIELDS_1D'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nmr1d(:,:)   = 0.0
mmr1d(:,:,:) = 0.0

SELECT CASE(i_glomap_clim_setup_in)
CASE (i_gc_sussocbc_5mode, i_gc_sussocbcdu_7mode)

  nmr1d(:,mode_nuc_sol) = gc_nd_nuc_sol_1d(:)
  mmr1d(:,mode_nuc_sol,cp_su)   = gc_nuc_sol_su_1d(:)
  mmr1d(:,mode_nuc_sol,cp_oc)   = gc_nuc_sol_oc_1d(:)
  nmr1d(:,mode_ait_sol) = gc_nd_ait_sol_1d(:)
  mmr1d(:,mode_ait_sol,cp_su)   = gc_ait_sol_su_1d(:)
  mmr1d(:,mode_ait_sol,cp_bc)   = gc_ait_sol_bc_1d(:)
  mmr1d(:,mode_ait_sol,cp_oc)   = gc_ait_sol_oc_1d(:)
  nmr1d(:,mode_acc_sol) = gc_nd_acc_sol_1d(:)
  mmr1d(:,mode_acc_sol,cp_su)   = gc_acc_sol_su_1d(:)
  mmr1d(:,mode_acc_sol,cp_bc)   = gc_acc_sol_bc_1d(:)
  mmr1d(:,mode_acc_sol,cp_oc)   = gc_acc_sol_oc_1d(:)
  mmr1d(:,mode_acc_sol,cp_cl)   = gc_acc_sol_ss_1d(:)
  nmr1d(:,mode_cor_sol) = gc_nd_cor_sol_1d(:)
  mmr1d(:,mode_cor_sol,cp_su)   = gc_cor_sol_su_1d(:)
  mmr1d(:,mode_cor_sol,cp_bc)   = gc_cor_sol_bc_1d(:)
  mmr1d(:,mode_cor_sol,cp_oc)   = gc_cor_sol_oc_1d(:)
  mmr1d(:,mode_cor_sol,cp_cl)   = gc_cor_sol_ss_1d(:)
  nmr1d(:,mode_ait_insol) = gc_nd_ait_ins_1d(:)
  mmr1d(:,mode_ait_insol,cp_bc) = gc_ait_ins_bc_1d(:)
  mmr1d(:,mode_ait_insol,cp_oc) = gc_ait_ins_oc_1d(:)

CASE DEFAULT
  errcode = 1
  WRITE(cmessage,'(A,I0,A)') 'i_glomap_clim_setup_in = ',                      &
                              i_glomap_clim_setup_in,                          &
                              newline // 'This option not available.'
  CALL ereport(RoutineName,errcode,cmessage)
END SELECT


SELECT CASE(i_glomap_clim_setup_in)
CASE (i_gc_sussocbcdu_7mode)
  mmr1d(:,mode_acc_sol,cp_du)   = gc_acc_sol_du_1d(:)
  mmr1d(:,mode_cor_sol,cp_du)   = gc_cor_sol_du_1d(:)
  nmr1d(:,mode_acc_insol)       = gc_nd_acc_ins_1d(:)
  mmr1d(:,mode_acc_insol,cp_du) = gc_acc_ins_du_1d(:)
  nmr1d(:,mode_cor_insol)       = gc_nd_cor_ins_1d(:)
  mmr1d(:,mode_cor_insol,cp_du) = gc_cor_ins_du_1d(:)
CASE DEFAULT
  ! Do nothing
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_gc_aerosol_fields_1d
END MODULE get_gc_aerosol_fields_1d_mod
