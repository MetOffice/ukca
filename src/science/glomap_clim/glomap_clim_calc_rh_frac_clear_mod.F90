! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution..
! *****************************COPYRIGHT*******************************
!
! Purpose:
!   To calculate clear sky relative humidity
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE glomap_clim_calc_rh_frac_clear_mod

USE um_types,                        ONLY:                                     &
    real_umphys

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='GLOMAP_CLIM_CALC_RH_FRAC_CLEAR_MOD'

PUBLIC :: glomap_clim_calc_rh_frac_clear

CONTAINS

SUBROUTINE glomap_clim_calc_rh_frac_clear( n_points,                           &
                                           q_1d,                               &
                                           qcf_1d,                             &
                                           cloud_liq_frac_1d,                  &
                                           cloud_blk_frac_1d,                  &
                                           t_theta_levels_1d,                  &
                                           p_theta_levels_1d,                  &
                                           rhcrit_1d,                          &
                                           rh_clr_1d )

USE lsp_subgrid_mod,                 ONLY:                                     &
    lsp_qclear

USE parkind1,                        ONLY:                                     &
    jprb,                                                                      &
    jpim

USE qsat_mod,                        ONLY:                                     &
    qsat,                                                                      &
    qsat_wat_mix

USE yomhook,                         ONLY:                                     &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

!                                   Number of grid boxes
INTEGER, INTENT(IN) :: n_points

!                                   Relative Humidity
REAL(KIND=real_umphys), INTENT(IN) :: q_1d( n_points )

!                                   Gridbox mean specific ice content
REAL(KIND=real_umphys), INTENT(IN) :: qcf_1d( n_points )

!                                   Cloud liquid fraction
REAL(KIND=real_umphys), INTENT(IN) :: cloud_liq_frac_1d( n_points )

!                                   Cloud bulk fraction
REAL(KIND=real_umphys), INTENT(IN) :: cloud_blk_frac_1d( n_points )

!                                   Temperature on theta levels
REAL(KIND=real_umphys), INTENT(IN) :: t_theta_levels_1d( n_points )

!                                   Pressure on theta levels
REAL(KIND=real_umphys), INTENT(IN) :: p_theta_levels_1d( n_points )

!                                   rhcrit
REAL(KIND=real_umphys), INTENT(IN) :: rhcrit_1d( n_points )

!                                   Clear sky relative humidity as a fraction
REAL(KIND=real_umphys), INTENT(IN OUT) :: rh_clr_1d( n_points )

! Local variables

! Saturated specific humidity 1-D
REAL(KIND=real_umphys) :: qsatmr_1d( n_points )

! Clear sky specific humidity 1-D
REAL(KIND=real_umphys) :: q_clr_1d( n_points )

! Sat. mixing ratio with respect to liquid water irrespective of temperature
REAL(KIND=real_umphys) :: qsatmr_wat_1d( n_points )

INTEGER :: loop

! Arbitrary value just less than one to prevent arithmatic issues later
REAL, PARAMETER :: less_than_one = 0.999

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*),   PARAMETER :: RoutineName='GLOMAP_CLIM_CALC_RH_FRAC_CLEAR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL qsat_wat_mix ( qsatmr_wat_1d, t_theta_levels_1d, p_theta_levels_1d,       &
                    n_points )

CALL qsat ( qsatmr_1d, t_theta_levels_1d, p_theta_levels_1d, n_points )

! Calculate clear sky relative humidity
CALL lsp_qclear( q_1d, qsatmr_1d, qsatmr_wat_1d,                               &
                 qcf_1d, cloud_liq_frac_1d, cloud_blk_frac_1d,                 &
                 rhcrit_1d, q_clr_1d, n_points )

! Calculate clear sky relative humidity as a fraction
DO loop = 1, n_points
  rh_clr_1d(loop) = q_clr_1d(loop) / qsatmr_wat_1d(loop)
END DO

DO loop = 1, n_points
  IF ( rh_clr_1d(loop) < 0.0 ) THEN
    rh_clr_1d(loop) = 0.0    ! remove negatives
  END IF
END DO

DO loop = 1, n_points
  IF ( rh_clr_1d(loop) > less_than_one ) THEN
    rh_clr_1d(loop) = less_than_one  ! remove values greater or equal to one
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE glomap_clim_calc_rh_frac_clear
END MODULE glomap_clim_calc_rh_frac_clear_mod
