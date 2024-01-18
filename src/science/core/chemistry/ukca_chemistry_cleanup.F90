! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Description:
!  Cleanup routine for chemistry
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_CHEMISTRY_CTL, UKCA_CHEMISTRY_CTL_COL and
!   UKCA_CHEMISTRY_CTL_TROPRAQ.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_chemistry_cleanup_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CHEMISTRY_CLEANUP_MOD'

CONTAINS

SUBROUTINE ukca_chemistry_cleanup(row_length, rows, model_levels, ntracers,    &
                                  secs_per_step, pres, drain, crain,           &
                                  latitude, env_ozone3d, qcf, r_theta_levels,  &
                                  mass, z_top_of_model, shno3_3d, nat_psc,     &
                                  tracer)

USE asad_mod,                      ONLY: jpctr

USE ukca_config_specification_mod, ONLY: ukca_config, i_top_2levH2O,           &
                                         i_top_1lev, i_top_BC
USE ukca_conserve_mod,             ONLY: ukca_conserve
USE ukca_cspecies,                 ONLY: n_h2o, n_hono2, nn_cl
USE ukca_sediment_mod,             ONLY: ukca_sediment
USE ukca_stratf_mod,               ONLY: ukca_stratf
USE ukca_topboundary_mod,          ONLY: ukca_top_boundary
USE ukca_tropopause,               ONLY: l_stratosphere

USE ereport_mod,                   ONLY: ereport
USE errormessagelength_mod,        ONLY: errormessagelength
USE parkind1,                      ONLY: jprb, jpim
USE yomhook,                       ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: ntracers

REAL, INTENT(IN)     :: secs_per_step                       ! time step
REAL, INTENT(IN)     :: pres(row_length,rows,model_levels)  ! pressure
REAL, INTENT(IN)     :: drain(row_length,rows,model_levels) ! 3-D LS rain
REAL, INTENT(IN)     :: crain(row_length,rows,model_levels) ! 3-D convec
REAL, INTENT(IN)     :: latitude(row_length,rows)           ! latitude (degrees)
REAL, INTENT(IN)     :: env_ozone3d(row_length,rows,model_levels) ! O3
REAL, INTENT(IN)     :: qcf(row_length,rows,model_levels)
REAL, INTENT(IN)     :: r_theta_levels(row_length,rows,0:model_levels)
REAL, INTENT(IN)     :: mass(row_length,rows,model_levels)  ! cell mass
REAL, INTENT(IN)     :: z_top_of_model                      ! top of model (m)
REAL, INTENT(IN OUT) :: shno3_3d(row_length,rows,model_levels)
REAL, INTENT(OUT)    :: nat_psc(row_length,                                    &
                                rows,model_levels)    ! Nitric acid trihydrate
                                                      ! (kg(nat)/kg(air))
REAL, INTENT(IN OUT) :: tracer(row_length,rows,                                &
                               model_levels,ntracers) ! tracer MMR
! Local variables

INTEGER :: i
INTEGER :: j
INTEGER :: k

INTEGER :: errcode ! Error code: ereport

! The call to ukca_conserve requires a logical to be set.
! ukca_conserve calculates and conserves total chlorine, bromine, and
! hydrogen. For these elements closed chemistry should be prescribed.
! Called after chemistry, after_chem, it rescales the chlorine, bromine
! and hydrogen containing compounds so that total chlorine, bromine
! and hydrogen are conserved under chemistry.
LOGICAL, PARAMETER :: after_chem = .FALSE.

CHARACTER(LEN=errormessagelength) :: cmessage         ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMISTRY_CLEANUP'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ukca_config%l_ukca_offline_be) RETURN

! Rescale bromine and chlorine tracers to guarantee conservation of total
! chlorine, bromine, and hydrogen over timestep. Only makes sense if at least
! chlorine chemistry is present.
IF (nn_cl > 0) THEN
  CALL ukca_conserve(row_length, rows, model_levels, ntracers,                 &
       tracer, pres, drain, crain, after_chem)
END IF

IF (ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_stratcfc .OR.             &
    ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat) THEN

  IF (ukca_config%l_ukca_het_psc) THEN
    ! Do NAT PSC sedimentation

    ! take NAT out of gasphase again
    tracer(:,:,:,n_hono2) = tracer(:,:,:,n_hono2) - shno3_3d

    CALL ukca_sediment(rows, row_length, model_levels, shno3_3d, qcf,          &
             r_theta_levels, mass, secs_per_step, l_stratosphere(:,:,:))

    ! add solid-phase HNO3 back to gasphase HNO3
    tracer(:,:,:,n_hono2) = tracer(:,:,:,n_hono2) + shno3_3d
  END IF

  ! i_ukca_topboundary==0 (i_top_none) corresponds to no overwriting of
  ! top level(s) or any top boundary condition
  IF (ukca_config%i_ukca_topboundary == i_top_2levH2O) THEN
    ! Tracer overwrites required to stop accumulation of tracer mass
    ! in the uppermost layers.  Exclude water vapour.
    DO i=1,rows
      DO j=1,row_length
        DO k=1,ntracers
          IF (k /= n_h2o) THEN
            tracer(j,i,model_levels  ,k) = tracer(j,i,model_levels-2,k)
            tracer(j,i,model_levels-1,k) = tracer(j,i,model_levels-2,k)
          END IF
        END DO
      END DO
    END DO
  ELSE IF (ukca_config%i_ukca_topboundary == i_top_1lev) THEN
    ! over-write top level for all tracers with 2nd-highest level
    DO i=1,rows
      DO j=1,row_length
        DO k=1,ntracers
          tracer(j,i,model_levels,k) = tracer(j,i,model_levels-1,k)
        END DO
      END DO
    END DO
  ELSE IF (ukca_config%i_ukca_topboundary >= i_top_BC) THEN
    ! Apply top boundary condition for NO, CO, O3, optionally H2O, using
    ! ACE-FTS climatologies (assumes constant latitude on each row)
    ! no tracer over-writing at top levels
    IF ((z_top_of_model > 85500.0) .OR. (z_top_of_model < 79000.0)) THEN
      cmessage='Can only impose a top boundary condition at 85 km.'
      errcode=25
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
    IF (.NOT. ALL(ABS(latitude(row_length,:) - latitude(1,:))                  &
                                                          < EPSILON(0.0))) THEN
      cmessage=                                                                &
        'Can only impose top boundary condition if latitude is constant on rows'
      errcode=1
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
    CALL ukca_top_boundary(row_length, rows, model_levels, ntracers,           &
                           latitude(1,:), tracer)
  END IF

  ! Copy NAT MMR into user_diagostics
  IF (ukca_config%l_ukca_het_psc) THEN
    nat_psc(:,:,:) = shno3_3d(:,:,:)
  END IF

ELSE IF (.NOT. ukca_config%l_ukca_offline) THEN   ! tropospheric chemistry
  ! Call routine to overwrite O3 and HNO3 species once per day
  ! above tropopause. Only for tropospheric chemistry

  CALL ukca_stratf(row_length,rows, model_levels,                              &
                   jpctr, env_ozone3d,                                         &
                   tracer(1:row_length,1:rows,1:model_levels,1:jpctr))

END IF     ! l_ukca_strat etc

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_chemistry_cleanup

END MODULE ukca_chemistry_cleanup_mod
