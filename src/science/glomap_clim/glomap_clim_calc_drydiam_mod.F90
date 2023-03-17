! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Populate drydiam field
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! Procedure:
!   1) CALL glomap_clim_identify_fields to identify item numbers required
!
!   2) CALL get_gc_aerosol_fields to populate the nmr1d and mmr1d fields
!
!   3) CALL ukca_calc_drydiam routine to calculate the drydiam field required
!       by ACTIVATE
!
! ---------------------------------------------------------------------
MODULE glomap_clim_calc_drydiam_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='GLOMAP_CLIM_CALC_DRYDIAM_MOD'

PUBLIC :: glomap_clim_calc_drydiam

CONTAINS

SUBROUTINE glomap_clim_calc_drydiam ( n_points, t_theta_levels_1d,             &
                                      p_theta_levels_1d )

USE atm_fields_mod,             ONLY:                                          &
    gc_nd_nuc_sol,                                                             &
    gc_nuc_sol_su,                                                             &
    gc_nuc_sol_oc,                                                             &
    gc_nd_ait_sol,                                                             &
    gc_ait_sol_su,                                                             &
    gc_ait_sol_bc,                                                             &
    gc_ait_sol_oc,                                                             &
    gc_nd_acc_sol,                                                             &
    gc_acc_sol_su,                                                             &
    gc_acc_sol_bc,                                                             &
    gc_acc_sol_oc,                                                             &
    gc_acc_sol_ss,                                                             &
    gc_nd_cor_sol,                                                             &
    gc_cor_sol_su,                                                             &
    gc_cor_sol_bc,                                                             &
    gc_cor_sol_oc,                                                             &
    gc_cor_sol_ss,                                                             &
    gc_nd_ait_ins,                                                             &
    gc_ait_ins_bc,                                                             &
    gc_ait_ins_oc

USE get_gc_aerosol_fields_mod,       ONLY:                                     &
    get_gc_aerosol_fields

USE glomap_clim_calc_aird_mod,       ONLY:                                     &
    glomap_clim_calc_aird

USE glomap_clim_fields_mod,          ONLY:                                     &
    aird,                                                                      &
    drydp,                                                                     &
    dvol,                                                                      &
    md,                                                                        &
    mdt,                                                                       &
    mmr1d,                                                                     &
    nd,                                                                        &
    nmr1d

USE glomap_clim_identify_fields_mod, ONLY:                                     &
    glomap_clim_identify_fields

USE glomap_clim_option_mod,          ONLY:                                     &
    i_glomap_clim_setup

USE glomap_clim_pop_md_mdt_nd_mod,   ONLY:                                     &
    glomap_clim_pop_md_mdt_nd

USE nlsizes_namelist_mod,            ONLY:                                     &
    row_length,                                                                &
    rows,                                                                      &
    model_levels

USE parkind1,                        ONLY:                                     &
    jpim,                                                                      &
    jprb

USE ukca_drydiam_field_mod,          ONLY:                                     &
    drydiam

USE ukca_calc_drydiam_mod,           ONLY:                                     &
    ukca_calc_drydiam

USE ukca_config_specification_mod,   ONLY:                                     &
    glomap_variables_climatology

USE ukca_mode_setup,                 ONLY:                                     &
    mode_list_sussbcoc_5mode,                                                  &
    nmodes,                                                                    &
    nmodes_list_sussbcoc_5mode

USE yomhook,                         ONLY:                                     &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: n_points

!                                   Pressure on theta levels (1D)
REAL(KIND=real_umphys),INTENT(IN) :: p_theta_levels_1d( n_points )

!                                   Temperature on theta levels (1D)
REAL(KIND=real_umphys),INTENT(IN) :: t_theta_levels_1d( n_points )

! Local variables

LOGICAL, SAVE :: firstcall = .TRUE.  ! first call to this routine t/f

INTEGER :: imode                      ! counter for modes
INTEGER :: i,j,k,loop,m               ! counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName='GLOMAP_CLIM_CALC_DRYDIAM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! This is not deallocated until ukca_activate
IF (.NOT. ALLOCATED(drydiam))                                                  &
           ALLOCATE(drydiam(row_length,rows,model_levels,nmodes))

! Find the information for items required from D1
IF ( firstcall ) THEN
  CALL glomap_clim_identify_fields()
  firstcall = .FALSE.
END IF

IF (.NOT. ALLOCATED(nd))            ALLOCATE(nd(n_points,nmodes))
IF (.NOT. ALLOCATED(md))                                                       &
          ALLOCATE( md(n_points,nmodes,glomap_variables_climatology%ncp))
IF (.NOT. ALLOCATED(mdt))           ALLOCATE(mdt(n_points,nmodes))

! Initialise arrays
nd(:,:)   = 0.0
md(:,:,:) = 0.0
mdt(:,:) = 0.0

IF (.NOT. ALLOCATED(aird))          ALLOCATE(aird(n_points))
IF (.NOT. ALLOCATED(nmr1d))         ALLOCATE(nmr1d(n_points,nmodes))
IF (.NOT. ALLOCATED(mmr1d))                                                    &
          ALLOCATE( mmr1d(n_points,nmodes,glomap_variables_climatology%ncp))

! Calculate aird
CALL glomap_clim_calc_aird( n_points, p_theta_levels_1d,                       &
                            t_theta_levels_1d, aird )

! Copy required fields from climatology aerosol array pointer into nmr1d & mmr1d
CALL get_gc_aerosol_fields( gc_nd_nuc_sol, gc_nuc_sol_su, gc_nuc_sol_oc,       &
                            gc_nd_ait_sol, gc_ait_sol_su, gc_ait_sol_bc,       &
                            gc_ait_sol_oc,                                     &
                            gc_nd_acc_sol, gc_acc_sol_su, gc_acc_sol_bc,       &
                            gc_acc_sol_oc, gc_acc_sol_ss,                      &
                            gc_nd_cor_sol, gc_cor_sol_su, gc_cor_sol_bc,       &
                            gc_cor_sol_oc, gc_cor_sol_ss,                      &
                            gc_nd_ait_ins, gc_ait_ins_bc, gc_ait_ins_oc,       &
                            n_points, mmr1d, nmr1d)

! Populate the md, mdt & nd arrays
CALL glomap_clim_pop_md_mdt_nd ( i_glomap_clim_setup, n_points, aird,          &
                                 mmr1d, nmr1d, md, mdt, nd )

IF (ALLOCATED(mmr1d))               DEALLOCATE(mmr1d)
IF (ALLOCATED(nmr1d))               DEALLOCATE(nmr1d)
IF (ALLOCATED(aird))                DEALLOCATE(aird)

IF (.NOT. ALLOCATED(drydp))         ALLOCATE(drydp(n_points,nmodes))
IF (.NOT. ALLOCATED(dvol))          ALLOCATE(dvol(n_points,nmodes))

! Calculate the dry diameters and volumes
CALL ukca_calc_drydiam( n_points, glomap_variables_climatology,                &
                        nd, md, mdt, drydp, dvol )


drydiam(:,:,:,:) = 0.0
DO m = 1, nmodes_list_sussbcoc_5mode
  imode = mode_list_sussbcoc_5mode(m)

  loop = 0
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        !Calculate vector position
        loop = loop + 1

        drydiam(i,j,k,imode) = drydp(loop,imode)

      END DO
    END DO
  END DO

END DO

IF (ALLOCATED(dvol))                DEALLOCATE(dvol)
IF (ALLOCATED(drydp))               DEALLOCATE(drydp)

IF (ALLOCATED(mdt))                 DEALLOCATE(mdt)
IF (ALLOCATED(md))                  DEALLOCATE(md)
IF (ALLOCATED(nd))                  DEALLOCATE(nd)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE glomap_clim_calc_drydiam

END MODULE glomap_clim_calc_drydiam_mod
