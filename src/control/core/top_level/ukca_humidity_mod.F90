! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Module for humidity-related calculations.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_humidity_mod

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

! Public procedures
PUBLIC ukca_vmrsat_liq, ukca_vmr_clear_sky

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_HUMIDITY_MOD'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE ukca_vmrsat_liq(row_length, rows, model_levels,                     &
                           temp, pres, vmr_sat, svp)
! ----------------------------------------------------------------------
! Description:
!
!   Returns the water vapour saturation mixing ratio and/or saturation vapour
!   pressure respect to liquid water.
!
! Method:
!
!   Calculations are based on Tetens formula giving saturation vapour
!   pressure in kPa as
!     svp = 0.611*exp(17.27*T/(T+237.3))
!   where T is temperature in degrees C. The equivalent expression for
!   saturation vapour pressure in Pa in terms of temperature in K is
!     svp = 611*exp(17.27*(T-273.15)/(T+35.85))                    Eq(1)
!   The equivalent mass mixing ratio is given by
!     vmr_sat = 0.62197*svp/(P-svp)                                Eq(2)
!   where P is atmospheric pressure.
!
!   (Eq(2) from re-arranging Eq A4.3 in Gill, A., Atmosphere-Ocean Dynamics
!    1982).
!
!   Combining Eq(1) and Eq(2) above gives the following expression for the
!   saturation mixing ratio in terms of pressure and temperature
!     vmr_sat = 380/(P/exp(17.27*(T-273.15)/(T-35.85))-611)       Eq (3)
! ----------------------------------------------------------------------

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE ukca_config_constants_mod, ONLY: repsilon
USE ukca_constants, ONLY: zerodegc

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
REAL, INTENT(IN) :: temp(row_length, rows, model_levels)
  ! Temperature (K)
REAL, OPTIONAL, INTENT(IN) :: pres(row_length, rows, model_levels)
  ! Pressure (Pa)
REAL, OPTIONAL, INTENT(OUT) :: vmr_sat(row_length, rows, model_levels)
  ! Water vapour saturation mixing ratio w.r.t. water (kg/kg dry air)
REAL, OPTIONAL, INTENT(OUT) :: svp(row_length, rows, model_levels)
  ! Saturation vapour pressure (Pa)

! Local variables

INTEGER :: i
INTEGER :: j
INTEGER :: k

REAL, PARAMETER :: t_factor = 17.27        ! Temperature function scaling
                                           ! factor in Tetens formula
REAL, PARAMETER :: t_offset = 35.85        ! Temperature offset in Tetens
                                           ! formula (K)
REAL, PARAMETER :: svp_fp = 611.0          ! Saturation vapour pressure at
                                           ! freezing point (Pa)
REAL, PARAMETER :: vmr_sat_factor = 380.0  ! Factor for vapour saturation m.r.
                                           ! calculation: svp_fp * 0.62197

REAL :: vmr_sat_limit
REAL :: temp_fn
REAL :: p_over_temp_fn

! ErrorStatus
INTEGER :: errcode                               ! Error flag
CHARACTER(LEN=errormessagelength) :: cmessage    ! Error return message

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_VMRSAT_LIQ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(vmr_sat) .AND. (.NOT. PRESENT(pres))) THEN
  errcode = 1
  cmessage = 'Pressure unavailable for calculating saturation mixing ratio'
  CALL ereport(RoutineName, errcode, cmessage)
END IF

! An upper limit for `vmr_sat` is set to match the maximum returned by the
! equivalent subroutine in the Unified Model ('qsat_wat_mix')
vmr_sat_limit = repsilon * 10.0

DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      IF (temp(i,j,k) > t_offset) THEN
        temp_fn = EXP(t_factor * (temp(i,j,k) - zerodegc) /                    &
                                 (temp(i,j,k) - t_offset))
        IF (PRESENT(vmr_sat)) THEN
          ! Calculate saturation mixing ratio using Eq(3)
          p_over_temp_fn = pres(i,j,k) / temp_fn
          IF (p_over_temp_fn > svp_fp) THEN
            vmr_sat(i,j,k) = MIN(vmr_sat_factor / (p_over_temp_fn - svp_fp),   &
                                 vmr_sat_limit)
          ELSE
            vmr_sat(i,j,k) = vmr_sat_limit
          END IF
        END IF
        IF (PRESENT(svp)) THEN
          ! Calculate saturation vapour pressure using Eq(1)
          svp(i,j,k) = svp_fp * temp_fn
        END IF
      ELSE
        errcode = 1
        WRITE(cmessage, '(A,E12.3,A,I0,A,I0,A,I0,A,E12.3)')                    &
          'Temperature (K) = ', temp(i,j,k),                                   &
          ' - too low for saturation vapour pressure calculation at (',        &
          i, ',', j, ',', k, '); T = ', temp(i,j,k)
        CALL ereport(RoutineName, errcode, cmessage)
      END IF
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_vmrsat_liq

! ----------------------------------------------------------------------
FUNCTION ukca_vmr_clear_sky(row_length, rows, model_levels, vmr, vmr_sat,      &
                            cloud_liq_frac) RESULT (vmr_clear)
! ----------------------------------------------------------------------
! Description:
!   Return the average water vapour mixing ratio (kg/kg dry air) for the
!   clear-sky fraction of each grid box.
!
! Method:
!
!   A simplified sub-grid partitioning of water vapour is assumed such
!   such that specific humidity is uniform everywhere outside of cloud
!   liquid water. i.e. it is the same for the clear-sky and ice-only
!   fractions of the grid box. An assumption of instantaneous condensation
!   fixes the vapour content within the grid box fraction containing
!   liquid cloud to the saturation content with respect to liquid water.
!   The clear sky mixing ratio vmr_clear is then given by re-arranging
!   the equation for the grid-box average:
!
!   vmr = (1 - cloud_liq_frac) * vmr_clear + cloud_liq_frac * vmr_sat
! ----------------------------------------------------------------------

IMPLICIT NONE

! Function arguments
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
REAL, INTENT(IN) :: vmr(row_length, rows, model_levels)
  ! Grid box average vapour mixing ratio (kg/kg)
REAL, INTENT(IN) :: vmr_sat(row_length, rows, model_levels)
  ! Saturation mixing ratio with respect to liquid water (kg/kg)
REAL, INTENT(IN) :: cloud_liq_frac(row_length, rows, model_levels)
  ! Cloud liquid fraction by volume

! Function result
REAL :: vmr_clear(row_length, rows, model_levels)

! Local variables

INTEGER :: i
INTEGER :: j
INTEGER :: k

! Dr Hook
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_VMR_CLEAR_SKY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      IF (cloud_liq_frac(i,j,k) < 1.0) THEN
        vmr_clear(i,j,k) =                                                     &
          (vmr(i,j,k) - cloud_liq_frac(i,j,k) * vmr_sat(i,j,k)) /              &
          (1.0 - cloud_liq_frac(i,j,k))
      ELSE
        vmr_clear(i,j,k) = vmr(i,j,k)
      END IF
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION ukca_vmr_clear_sky

END MODULE ukca_humidity_mod

