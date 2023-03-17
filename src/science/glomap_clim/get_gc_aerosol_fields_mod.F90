! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution..
! *****************************COPYRIGHT*******************************
!
! Purpose:
!   To put fields required by RADAER into stash.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE get_gc_aerosol_fields_mod

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='GET_GC_AEROSOL_FIELDS_MOD'

PUBLIC :: get_gc_aerosol_fields

CONTAINS

SUBROUTINE get_gc_aerosol_fields( gc_nd_nuc_sol, gc_nuc_sol_su, gc_nuc_sol_oc, &
                                  gc_nd_ait_sol, gc_ait_sol_su, gc_ait_sol_bc, &
                                  gc_ait_sol_oc,                               &
                                  gc_nd_acc_sol, gc_acc_sol_su, gc_acc_sol_bc, &
                                  gc_acc_sol_oc, gc_acc_sol_ss,                &
                                  gc_nd_cor_sol, gc_cor_sol_su, gc_cor_sol_bc, &
                                  gc_cor_sol_oc, gc_cor_sol_ss,                &
                                  gc_nd_ait_ins, gc_ait_ins_bc, gc_ait_ins_oc, &
                                  n_points, mmr1d, nmr1d)

! Copy required fields from climatology aerosol array pointer into nmr1d & mmr1d

USE atm_fields_bounds_mod,  ONLY:                                              &
    tdims

USE ereport_mod,            ONLY:                                              &
    ereport

USE errormessagelength_mod, ONLY:                                              &
    errormessagelength

USE glomap_clim_option_mod, ONLY:                                              &
    i_glomap_clim_setup

USE parkind1,               ONLY:                                              &
    jprb,                                                                      &
    jpim

USE ukca_mode_setup,        ONLY:                                              &
    cp_su,                                                                     &
    cp_bc,                                                                     &
    cp_oc,                                                                     &
    cp_cl,                                                                     &
    mode_nuc_sol,                                                              &
    mode_ait_sol,                                                              &
    mode_acc_sol,                                                              &
    mode_cor_sol,                                                              &
    mode_ait_insol,                                                            &
    nmodes

USE ukca_config_specification_mod, ONLY:                                       &
    glomap_variables_climatology,                                              &
    i_sussbcoc_5mode

USE umPrintMgr,             ONLY:                                              &
    newline

USE um_types,               ONLY:                                              &
    real_umphys

USE yomhook,                ONLY:                                              &
    lhook,                                                                     &
    dr_hook


IMPLICIT NONE

! Arguments

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_nuc_sol(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_nuc_sol_su(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_nuc_sol_oc(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_ait_sol(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_sol_su(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_sol_bc(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_sol_oc(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_acc_sol(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_su(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_bc(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_oc(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_acc_sol_ss(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_cor_sol(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_su(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_bc(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_oc(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_cor_sol_ss(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_nd_ait_ins(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_ins_bc(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

REAL(KIND=real_umphys),INTENT(IN) :: gc_ait_ins_oc(tdims%i_start:tdims%i_end,  &
                                                   tdims%j_start:tdims%j_end,  &
                                                   tdims%k_start:tdims%k_end)

INTEGER, INTENT(IN)                 :: n_points
REAL(KIND=real_umphys), INTENT(OUT) :: mmr1d(n_points, nmodes,                 &
                                             glomap_variables_climatology%ncp)
REAL(KIND=real_umphys), INTENT(OUT) :: nmr1d(n_points, nmodes)

! Local Variables

INTEGER                           :: i,j,k,loop
INTEGER                           :: jj, kk
INTEGER                           :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='GET_GC_AEROSOL_FIELDS'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nmr1d(:,:)   = 0.0
mmr1d(:,:,:) = 0.0

SELECT CASE(i_glomap_clim_setup)
CASE (i_sussbcoc_5mode)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,loop,jj,kk)     &
!$OMP SHARED(tdims, nmr1d, mmr1d, gc_nd_nuc_sol, gc_nuc_sol_su, gc_nuc_sol_oc, &
!$OMP        gc_nd_ait_sol, gc_ait_sol_su, gc_ait_sol_bc, gc_ait_sol_oc,       &
!$OMP        gc_nd_acc_sol, gc_acc_sol_su, gc_acc_sol_bc, gc_acc_sol_oc,       &
!$OMP        gc_acc_sol_ss, gc_nd_cor_sol, gc_cor_sol_su, gc_cor_sol_bc,       &
!$OMP        gc_cor_sol_oc, gc_cor_sol_ss, gc_nd_ait_ins, gc_ait_ins_bc,       &
!$OMP        gc_ait_ins_oc)
  DO k = 1, tdims%k_end
    kk= (k-1) * tdims%j_len * tdims%i_len
    DO j = 1, tdims%j_end
      jj = (j-1) * tdims%i_len
      DO i = 1, tdims%i_end

        ! Calculate vector position
        loop = i + jj + kk

        nmr1d(loop,mode_nuc_sol)   = gc_nd_nuc_sol(i,j,k)
        nmr1d(loop,mode_ait_sol)   = gc_nd_ait_sol(i,j,k)
        nmr1d(loop,mode_acc_sol)   = gc_nd_acc_sol(i,j,k)
        nmr1d(loop,mode_cor_sol)   = gc_nd_cor_sol(i,j,k)
        nmr1d(loop,mode_ait_insol) = gc_nd_ait_ins(i,j,k)

        mmr1d(loop,mode_nuc_sol,cp_su)   = gc_nuc_sol_su(i,j,k)
        mmr1d(loop,mode_nuc_sol,cp_oc)   = gc_nuc_sol_oc(i,j,k)
        mmr1d(loop,mode_ait_sol,cp_su)   = gc_ait_sol_su(i,j,k)
        mmr1d(loop,mode_ait_sol,cp_bc)   = gc_ait_sol_bc(i,j,k)
        mmr1d(loop,mode_ait_sol,cp_oc)   = gc_ait_sol_oc(i,j,k)
        mmr1d(loop,mode_acc_sol,cp_su)   = gc_acc_sol_su(i,j,k)
        mmr1d(loop,mode_acc_sol,cp_bc)   = gc_acc_sol_bc(i,j,k)
        mmr1d(loop,mode_acc_sol,cp_oc)   = gc_acc_sol_oc(i,j,k)
        mmr1d(loop,mode_acc_sol,cp_cl)   = gc_acc_sol_ss(i,j,k)
        mmr1d(loop,mode_cor_sol,cp_su)   = gc_cor_sol_su(i,j,k)
        mmr1d(loop,mode_cor_sol,cp_bc)   = gc_cor_sol_bc(i,j,k)
        mmr1d(loop,mode_cor_sol,cp_oc)   = gc_cor_sol_oc(i,j,k)
        mmr1d(loop,mode_cor_sol,cp_cl)   = gc_cor_sol_ss(i,j,k)
        mmr1d(loop,mode_ait_insol,cp_bc) = gc_ait_ins_bc(i,j,k)
        mmr1d(loop,mode_ait_insol,cp_oc) = gc_ait_ins_oc(i,j,k)

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO


CASE DEFAULT
  errcode = 1
  WRITE(cmessage,'(A,I0,A)') 'i_glomap_clim_setup = ', i_glomap_clim_setup,    &
                              newline// 'This option not available.'
  CALL ereport(RoutineName,errcode,cmessage)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_gc_aerosol_fields

END MODULE get_gc_aerosol_fields_mod
