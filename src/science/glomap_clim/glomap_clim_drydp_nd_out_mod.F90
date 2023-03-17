! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Output drydp and nd, populate drydiam field
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

MODULE glomap_clim_drydp_nd_out_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='GLOMAP_CLIM_DRYDP_ND_OUT_MOD'

CONTAINS

SUBROUTINE glomap_clim_drydp_nd_out (                                          &
                                  gc_nd_nuc_sol, gc_nuc_sol_su, gc_nuc_sol_oc, &
                                  gc_nd_ait_sol, gc_ait_sol_su, gc_ait_sol_bc, &
                                  gc_ait_sol_oc,                               &
                                  gc_nd_acc_sol, gc_acc_sol_su, gc_acc_sol_bc, &
                                  gc_acc_sol_oc, gc_acc_sol_ss,                &
                                  gc_nd_cor_sol, gc_cor_sol_su, gc_cor_sol_bc, &
                                  gc_cor_sol_oc, gc_cor_sol_ss,                &
                                  gc_nd_ait_ins, gc_ait_ins_bc, gc_ait_ins_oc, &
                                  n_points,                                    &
                                  p_theta_levels_1d, t_theta_levels_1d,        &
                                  drydp, nd )

USE atm_fields_bounds_mod,           ONLY:                                     &
    tdims

USE get_gc_aerosol_fields_mod,       ONLY:                                     &
    get_gc_aerosol_fields

USE glomap_clim_calc_aird_mod,       ONLY:                                     &
    glomap_clim_calc_aird

#if !defined(LFRIC)
USE glomap_clim_identify_fields_mod, ONLY:                                     &
    glomap_clim_identify_fields
#endif

USE glomap_clim_option_mod,          ONLY:                                     &
    i_glomap_clim_setup

USE glomap_clim_pop_md_mdt_nd_mod,   ONLY:                                     &
    glomap_clim_pop_md_mdt_nd

USE parkind1,                        ONLY:                                     &
    jpim,                                                                      &
    jprb

#if !defined(LFRIC)
USE ukca_drydiam_field_mod,          ONLY:                                     &
    drydiam
#endif

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


INTEGER, INTENT(IN) :: n_points ! Number of points in 3D field

!                                   Pressure on theta levels (1D)
REAL(KIND=real_umphys),INTENT(IN) :: p_theta_levels_1d( n_points )

!                                   Temperature on theta levels (1D)
REAL(KIND=real_umphys),INTENT(IN) :: t_theta_levels_1d( n_points )

REAL, INTENT(OUT)   :: drydp(n_points,nmodes)
REAL, INTENT(OUT)   :: nd(n_points,nmodes)

! Local variables

#if !defined(LFRIC)
LOGICAL, SAVE        :: firstcall = .TRUE.  ! first call to this routine t/f

INTEGER :: imode                      ! counter for modes
INTEGER :: i,j,k,loop,m               ! counter
INTEGER :: jj, kk
#endif

REAL(KIND=real_umphys), ALLOCATABLE :: aird(:)
REAL(KIND=real_umphys), ALLOCATABLE :: dvol(:,:)
REAL(KIND=real_umphys), ALLOCATABLE :: md(:,:,:)
REAL(KIND=real_umphys), ALLOCATABLE :: mdt(:,:)
REAL(KIND=real_umphys), ALLOCATABLE :: mmr1d(:,:,:)
REAL(KIND=real_umphys), ALLOCATABLE :: nmr1d(:,:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName='GLOMAP_CLIM_DRYDP_ND_OUT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

#if !defined(LFRIC)
! Find the information for items required from D1
IF ( firstcall ) THEN
  CALL glomap_clim_identify_fields()
  firstcall = .FALSE.
END IF
#endif

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
          ALLOCATE( mmr1d(n_points,nmodes, glomap_variables_climatology%ncp))

! Calculate aird
CALL glomap_clim_calc_aird( n_points, p_theta_levels_1d, t_theta_levels_1d,    &
                            aird )

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

IF (.NOT. ALLOCATED(dvol))          ALLOCATE(dvol(n_points,nmodes))

! Calculate the dry diameters and volumes
CALL ukca_calc_drydiam( n_points, glomap_variables_climatology,                &
                        nd, md, mdt, drydp, dvol )

#if !defined(LFRIC)

DO m = 1, nmodes
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k)                    &
!$OMP SHARED(tdims, drydiam, m)
  DO k = 1, tdims%k_end
    drydiam(:,:,k,m) = 0.0
  END DO
!$OMP END PARALLEL DO
END DO

DO m = 1, nmodes_list_sussbcoc_5mode
  imode = mode_list_sussbcoc_5mode(m)


!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,loop,jj,kk)     &
!$OMP SHARED(tdims, drydiam, drydp, imode)
  DO k = 1, tdims%k_end
    kk = (k-1) * tdims%j_len * tdims%i_len
    DO j = 1, tdims%j_end
      jj = (j-1) * tdims%i_len
      DO i = 1, tdims%i_end

        ! Calculate vector position
        loop = i + jj +kk

        drydiam(i,j,k,imode) = drydp(loop,imode)

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END DO
#endif

IF (ALLOCATED(dvol))                DEALLOCATE(dvol)

IF (ALLOCATED(mdt))                 DEALLOCATE(mdt)
IF (ALLOCATED(md))                  DEALLOCATE(md)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE glomap_clim_drydp_nd_out

END MODULE glomap_clim_drydp_nd_out_mod
