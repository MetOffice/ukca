! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Local copy of trsrce-trsrce2a.F90 to calculate air density and burden
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------

MODULE ukca_calc_rho_mod

USE ukca_um_legacy_mod, ONLY: kappa, c_virtual, pref, cp


IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CALC_RHO_MOD'

CONTAINS

SUBROUTINE ukca_calc_rho(                                                      &
 row_length, rows, model_levels, theta, q, qcl, qcf, exner, rho_r2,            &
 r_theta_levels, r_rho_levels, air_density, air_burden)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN)    :: row_length
INTEGER, INTENT(IN)    :: rows
INTEGER, INTENT(IN)    :: model_levels

REAL, INTENT(IN)       ::                                                      &
      theta (1:row_length,1:rows,model_levels)                                 &
                                                 ! pot temp
,     q     (1:row_length,1:rows,model_levels)                                 &
                                                 ! Q on theta levs
,     qcl   (1:row_length,1:rows,model_levels)                                 &
                                                 ! Qcl on theta levs
,     qcf   (1:row_length,1:rows,model_levels)                                 &
                                                 ! Qcf on theta levs
,     exner (1:row_length,1:rows,model_levels+1)                               &
                                                 ! exner on rho levs
,     rho_r2(1:row_length,1:rows,model_levels)                                 &
                                                 ! density * r * r
                                                 ! on rho levs
,     r_theta_levels(row_length, rows, 0:model_levels)                         &
                                                 ! height of theta levels
,     r_rho_levels(row_length, rows, model_levels)
                                                 ! height of rho levels

REAL, INTENT(IN OUT)       ::                                                  &
      air_density (1:row_length,1:rows,model_levels)                           &
,     air_burden  (1:row_length,1:rows,model_levels)

! Local variables

REAL           :: dm, am
REAL           :: thetav       ! virtual potential temperature
REAL           :: exner_ave    ! an averaged exner term
REAL           :: rho_theta    ! rho on theta level
REAL           :: rho1         ! } values of rho after the
REAL           :: rho2         ! } r-squared factor removed

INTEGER :: i, j, k   ! Loop countera

! Error Reporting
CHARACTER (LEN=*), PARAMETER :: RoutineName='UKCA_CALC_RHO'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

air_density(:,:,:) = 0.0
air_burden (:,:,:) = 0.0

DO k = 1, model_levels

  IF (k < model_levels) THEN

    DO j = 1, rows
      DO i = 1, row_length

        ! Remove the r squared factor from rho before interpolation
        rho1 =  rho_r2(i,j,k) / (r_rho_levels(i,j,k) *                         &
                                 r_rho_levels(i,j,k) )
        rho2 =  rho_r2(i,j,k+1) / (r_rho_levels(i,j,k+1) *                     &
                                   r_rho_levels(i,j,k+1))

        ! DM = density (interpolated on to theta levels) * delta r
        !    = burden (kg m-2)
        dm = rho2 * (r_theta_levels(i,j,k) - r_rho_levels(i,j,k)) +            &
             rho1 * (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))

        ! AM = density interpolated on to theta levels (kg m-3)
        am = dm / (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))

        ! Special case for lowest layer to get correct mass
        IF (k == 1) THEN
          dm = dm * (r_rho_levels(i,j,2) - r_theta_levels(i,j,0)) /            &
                    (r_rho_levels(i,j,2) - r_rho_levels(i,j,1))
        END IF
        !
        ! Convert DM to DRY density and store in output arrays
        air_burden(i,j,k)  = dm * (1.0 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k))
        air_density(i,j,k) = am * (1.0 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k))

      END DO ! i
    END DO ! j

  ELSE ! k = model_level
       !--------------------------------------------------------------------
       ! Cannot average here to get rho_theta. Hence calculate by using
       ! the equation of state vis
       !         ___r                     P0
       !         rho   = -----------------------------------------
       !                 kappa * Cp * ____________________r * theta_v
       !                              [          kappa-1 ]
       !                              [ exner ** ------- ]
       !                               [          kappa   ]
       ! -------------------------------------------------------------------

    DO j = 1, rows
      DO i = 1, row_length

        thetav = theta(i,j,k) *                                                &
                 (1.0 + (q(i,j,k)*c_Virtual) - qcl(i,j,k) - qcf(i,j,k))

        exner_ave = (exner(i,j,k) ** ((kappa - 1.0)/kappa) +                   &
                     exner(i,j,k+1) ** ((kappa - 1.0)/kappa))/2.0
        rho_theta = pref / ( kappa * cp * exner_ave * thetav )

        ! rho_theta is at the top theta level. We also need the value
        ! at the top rho level. This will be rho1
        rho1 = rho_r2(i,j,model_levels) / (r_rho_levels(i,j,model_levels) *    &
                                           r_rho_levels(i,j,model_levels))

        ! rho2 will be the average of rho1 and rho_theta
        rho2 = ( rho1 + rho_theta ) * 0.5

        air_burden(i,j,k)  = rho2 * (r_theta_levels(i,j,k) -                   &
                                     r_rho_levels(i,j,k))
        air_density(i,j,k) = rho2

      END DO ! i
    END DO ! j

  END IF  ! k < model_levels

END DO ! k

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_calc_rho
END MODULE ukca_calc_rho_mod
