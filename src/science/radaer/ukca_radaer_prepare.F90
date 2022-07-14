! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Re-arrange UKCA-MODE input to match the expectations of
!  routine ukca_radaer_band_average().
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
MODULE ukca_radaer_prepare_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_RADAER_PREPARE_MOD'

CONTAINS

SUBROUTINE ukca_radaer_prepare(                                                &
      ! Input Actual array dimensions
      n_profile, n_layer, n_ukca_mode, n_ukca_cpnt,                            &
      ! Input Fixed array dimensions
      npd_profile, npd_layer, npd_aerosol_mode,                                &
      ! Input from the UKCA_RADAER structure
      nmodes, ncp_max,                                                         &
      i_cpnt_index, n_cpnt_in_mode,                                            &
      ! Input Component mass-mixing ratios
      ukca_mix_ratio,                                                          &
      ! Input modal number concentrations
      ukca_modal_nbr,                                                          &
      ! Input Pressure and temperature
      pressure, temperature,                                                   &
      ! Output modal mass-mixing ratios
      ukca_modal_mixr,                                                         &
      ! Output modal number concentrations
      ukca_modal_number )

USE chemistry_constants_mod, ONLY: boltzmann
USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

! Arguments

!
! Actual array dimensions
!
INTEGER, INTENT(IN) :: n_profile
INTEGER, INTENT(IN) :: n_layer
INTEGER, INTENT(IN) :: n_ukca_mode
INTEGER, INTENT(IN) :: n_ukca_cpnt
!
! Fixed array dimensions
!
INTEGER, INTENT(IN) :: npd_profile
INTEGER, INTENT(IN) :: npd_layer
INTEGER, INTENT(IN) :: npd_aerosol_mode
!
! Structure for UKCA/radiation interaction
!
INTEGER, INTENT(IN) :: nmodes
INTEGER, INTENT(IN) :: ncp_max
INTEGER, INTENT(IN) :: i_cpnt_index( ncp_max, nmodes )
INTEGER, INTENT(IN) :: n_cpnt_in_mode( nmodes )
!
! Component mass-mixing ratios
!
REAL, INTENT(IN) :: ukca_mix_ratio(n_ukca_cpnt, npd_profile, npd_layer)
!
! Modal number concentrations divided by molecular concentration of air
!
REAL, INTENT(IN) :: ukca_modal_nbr(npd_profile, npd_layer, n_ukca_mode)
!
! Pressure and temperature fields.
REAL, INTENT(IN) :: pressure(npd_profile, npd_layer)
REAL, INTENT(IN) :: temperature(npd_profile, npd_layer)
!
! Modal mass-mixing ratios
!
REAL, INTENT(OUT) :: ukca_modal_mixr(npd_profile, npd_layer, npd_aerosol_mode)
!
! Modal number concentrations (in m-3)
!
REAL, INTENT(OUT) :: ukca_modal_number(npd_profile, npd_layer, n_ukca_mode)

!
! Local variables
!

INTEGER :: i,                                                                  &
        j,                                                                     &
        k,                                                                     &
        profile
INTEGER :: this_cpnt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_PREPARE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)




!
! Modal mass-mixing ratios.
!
! Simply sum the mixing ratios of all components included in a
! given mode to get the mixing ratio for that mode.
!
DO j = 1, n_ukca_mode

  DO k = 1, n_layer

    DO profile = 1, n_profile

      ukca_modal_mixr(profile, k, j) = 0.0

    END DO ! profile

  END DO ! k

END DO ! j

DO j = 1, n_ukca_mode

  DO i = 1, n_cpnt_in_mode(j)

    this_cpnt = i_cpnt_index(i, j)

    DO k = 1, n_layer

      DO profile = 1, n_profile

        ukca_modal_mixr(profile, k, j) =                                       &
          ukca_modal_mixr(profile, k, j) +                                     &
          ukca_mix_ratio(this_cpnt, profile, k)

      END DO ! profile

    END DO ! k

  END DO ! i

END DO ! j

!
! Modal number concentrations
!
! Multiply by the molecular concentration of air (p/kT) to obtain
! the acutal aerosol number concentrations.
!
DO j = 1, n_ukca_mode

  DO k = 1, n_layer

    DO profile = 1, n_profile

      ukca_modal_number(profile, k, j) = ukca_modal_nbr(profile, k, j) *       &
                                         pressure(profile, k) /                &
                                         (boltzmann * temperature(profile, k))

    END DO ! profile

  END DO ! k

END DO ! j

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_prepare
END MODULE ukca_radaer_prepare_mod
