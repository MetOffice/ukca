! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  An interface routine to pass chunked columns from lfric to radaer
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
MODULE ukca_radaer_lfric_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_RADAER_LFRIC_MOD'

CONTAINS
       
SUBROUTINE ukca_radaer_lfric_interface(                                        &
    ! Fixed array dimensions (input)
    npd_profile,                                                               &
    npd_layer,                                                                 &
    npd_exclude_lw,                                                            &
    npd_exclude_sw,                                                            &
    npd_ukca_aod_wavel,                                                        &
    ! Spectral information (input)
    ip_infra_red,                                                              &
    ip_solar,                                                                  &
    ! Actual array dimensions (input)
    n_ukca_mode,                                                               &
    n_ukca_cpnt,                                                               &
    ! Prescribed SSA dimensions
    nd_prof_ssa,                                                               &
    nd_layr_ssa,                                                               &
    nd_band_ssa,                                                               &
    ! Variables related to waveband exclusion
    l_exclude_lw,                                                              &
    l_exclude_sw,                                                              &
    ! UKCA_RADAER structure (input)
    nmodes,                                                                    &
    ncp_max,                                                                   &
    ncp_max_x_nmodes,                                                          &
    i_cpnt_index,                                                              &
    i_cpnt_type,                                                               &
    i_mode_type,                                                               &
    l_nitrate,                                                                 &
    l_soluble,                                                                 &
    l_sustrat,                                                                 &
    l_cornarrow_ins,                                                           &
    n_cpnt_in_mode,                                                            &
    ! Modal diameters from UKCA module (input)
    ukca_dry_diam_um,                                                          &
    ukca_wet_diam_um,                                                          &
    ! Other inputs from UKCA module (input)
    ukca_comp_vol_um,                                                          &
    ukca_modal_vol_um,                                                         &
    ukca_modal_rho_um,                                                         &
    ukca_modal_wtv_um,                                                         &
    ! Logical to describe orientation
    l_inverted,                                                                &
    ! Control option for prescribed single scattering albedo array
    i_ukca_radaer_prescribe_ssa,                                               &
    ! Model level of the tropopause (input)
    trindxrad_um,                                                              &
    ! Whether we need to run shortwave band_average because lit or not
    l_any_lit_points_um,                                                       &
    ! Prescription of single-scattering albedo
    ukca_radaer_presc_ssa,                                                     &
    ! Input Component mass-mixing ratios
    ukca_mix_ratio_um,                                                         &
    ! Input modal number concentrations
    ukca_modal_nbr_um,                                                         &
    ! Input Pressure and temperature
    p_theta_levels, t_theta_levels,                                            &
    ! Maxwell-Garnett mixing approach logical control switches
    i_ukca_tune_bc, i_glomap_clim_tune_bc,                                     &
    ! Type selection
    soluble_wanted,                                                            &
    soluble_unwanted,                                                          &
    ! Which aerosol optical depth diagnostics to calculate
    l_aod_ukca_ait_sol, l_aaod_ukca_ait_sol,                                   &
    l_aod_ukca_acc_sol, l_aaod_ukca_acc_sol,                                   &
    l_aod_ukca_cor_sol, l_aaod_ukca_cor_sol,                                   &
    l_aod_ukca_ait_ins, l_aaod_ukca_ait_ins,                                   &
    l_aod_ukca_acc_ins, l_aaod_ukca_acc_ins,                                   &
    l_aod_ukca_cor_ins, l_aaod_ukca_cor_ins,                                   &
    ! Mass thickness of layers
    d_mass_theta_levels_um,                                                    &
    ! Modal mass-mixing ratios (input output)
    ukca_mode_mix_ratio_um,                                                    &
    ! Band-averaged optical properties (output)
    aer_lw_absorption_um,                                                      &
    aer_sw_absorption_um,                                                      &
    aer_lw_scattering_um,                                                      &
    aer_sw_scattering_um,                                                      &
    aer_lw_asymmetry_um,                                                       &
    aer_sw_asymmetry_um,                                                       &
    aod_ukca_all_modes_um,                                                     &
    aaod_ukca_all_modes_um )

USE socrates_init_mod,                 ONLY: n_sw_band,                        &
                                             sw_n_band_exclude,                &
                                             sw_index_exclude,                 &
                                             n_lw_band,                        &
                                             lw_n_band_exclude,                &
                                             lw_index_exclude

USE ukca_mode_setup,                   ONLY: mode_ait_sol, mode_acc_sol,       &
                                             mode_cor_sol, mode_ait_insol,     &
                                             mode_acc_insol, mode_cor_insol,   &
                                             ip_ukca_mode_aitken,              &
                                             ip_ukca_mode_accum,               &
                                             ip_ukca_mode_coarse

USE ukca_radaer_band_average_mod,      ONLY: ukca_radaer_band_average

USE ukca_radaer_compute_aod_mod,       ONLY: ukca_radaer_compute_aod

USE ukca_radaer_prepare_mod,           ONLY: ukca_radaer_prepare

USE um_physics_init_mod,               ONLY: n_radaer_mode

USE parkind1,                          ONLY: jpim, jprb
USE yomhook,                           ONLY: lhook, dr_hook

IMPLICIT NONE

! Arguments

! Fixed array dimensions
INTEGER, INTENT(IN) :: npd_profile, npd_layer, npd_exclude_lw,                 &
                       npd_exclude_sw, npd_ukca_aod_wavel

! Spectral information (input)
INTEGER, INTENT(IN) :: ip_infra_red, ip_solar

! RADAER array dimensions (note that nucleation mode is excluded)
INTEGER, INTENT(IN) :: n_ukca_mode, n_ukca_cpnt

! Fixed array dimensions for prescribed SSA
INTEGER, INTENT(IN) :: nd_prof_ssa, nd_layr_ssa, nd_band_ssa

! Variables related to waveband exclusion
LOGICAL, INTENT(IN) :: l_exclude_lw, l_exclude_sw

! From ukca_radaer Structure for UKCA/radiation interaction
INTEGER, INTENT(IN) :: nmodes
INTEGER, INTENT(IN) :: ncp_max
INTEGER, INTENT(IN) :: ncp_max_x_nmodes
INTEGER, INTENT(IN) :: i_cpnt_index( ncp_max, nmodes )
INTEGER, INTENT(IN) :: i_cpnt_type( ncp_max_x_nmodes )
INTEGER, INTENT(IN) :: i_mode_type( nmodes )
LOGICAL, INTENT(IN) :: l_nitrate
LOGICAL, INTENT(IN) :: l_soluble( nmodes )
LOGICAL, INTENT(IN) :: l_sustrat
LOGICAL, INTENT(IN) :: l_cornarrow_ins
INTEGER, INTENT(IN) :: n_cpnt_in_mode( nmodes )

! Modal diameters from UKCA module (input)
REAL, INTENT(IN) ::  ukca_dry_diam_um(npd_profile, npd_layer, n_ukca_mode)
REAL, INTENT(IN) ::  ukca_wet_diam_um(npd_profile, npd_layer, n_ukca_mode)

! Component volume
REAL, INTENT(IN) ::  ukca_comp_vol_um( n_ukca_cpnt, npd_profile, npd_layer )

! Modal volumes and densities
REAL, INTENT(IN) :: ukca_modal_vol_um( npd_profile, npd_layer, n_ukca_mode )
REAL, INTENT(IN) :: ukca_modal_rho_um( npd_profile, npd_layer, n_ukca_mode )

! Volume of water in modes
REAL, INTENT(IN) :: ukca_modal_wtv_um( npd_profile, npd_layer, n_ukca_mode )

! Logical to describe orientation
LOGICAL, INTENT(IN) :: l_inverted

! When > 0, use a prescribed single scattering albedo field
INTEGER, INTENT(IN) :: i_ukca_radaer_prescribe_ssa

! Model level of tropopause
INTEGER, INTENT(IN) :: trindxrad_um(npd_profile)

! Whether we need to run shortwave band_average because lit or not
LOGICAL, INTENT(IN) :: l_any_lit_points_um

! Prescription of single-scattering albedo
REAL, INTENT(IN) :: ukca_radaer_presc_ssa( nd_prof_ssa, nd_layr_ssa,           &
                                           nd_band_ssa )

! Component mass-mixing ratios
REAL, INTENT(IN) :: ukca_mix_ratio_um( n_ukca_cpnt, npd_profile, npd_layer )

! Modal number concentrations divided by molecular concentration of air
REAL, INTENT(IN) :: ukca_modal_nbr_um( n_ukca_cpnt, npd_profile, npd_layer )

! pressure on theta levels
REAL, INTENT(IN) :: p_theta_levels( npd_profile, npd_layer )

! temperature on theta levels
REAL, INTENT(IN) :: t_theta_levels( npd_profile, npd_layer )

! Maxwell-Garnett mixing approach logical control switches
INTEGER, INTENT(IN) :: i_ukca_tune_bc
INTEGER, INTENT(IN) :: i_glomap_clim_tune_bc

! Type selection
LOGICAL, INTENT(IN) :: soluble_wanted
LOGICAL, INTENT(IN) :: soluble_unwanted

! Which aerosol optical depth diagnostics to calculate
LOGICAL, INTENT(IN) :: l_aod_ukca_ait_sol, l_aaod_ukca_ait_sol,                &
                       l_aod_ukca_acc_sol, l_aaod_ukca_acc_sol,                &
                       l_aod_ukca_cor_sol, l_aaod_ukca_cor_sol,                &
                       l_aod_ukca_ait_ins, l_aaod_ukca_ait_ins,                &
                       l_aod_ukca_acc_ins, l_aaod_ukca_acc_ins,                &
                       l_aod_ukca_cor_ins, l_aaod_ukca_cor_ins

! Mass thickness of layers
REAL, INTENT(IN) :: d_mass_theta_levels_um( npd_profile, npd_layer )

! Modal mass-mixing ratios
REAL, INTENT(IN OUT) ::  ukca_mode_mix_ratio_um( npd_profile, npd_layer,       &
                                             n_radaer_mode )

! Band-averaged modal optical properties
REAL, INTENT(IN OUT) :: aer_lw_absorption_um( npd_profile, npd_layer,          &
                                              n_radaer_mode, n_lw_band )

REAL, INTENT(IN OUT) :: aer_sw_absorption_um( npd_profile, npd_layer,          &
                                              n_radaer_mode, n_sw_band )

REAL, INTENT(IN OUT) :: aer_lw_scattering_um( npd_profile, npd_layer,          &
                                              n_radaer_mode, n_lw_band )

REAL, INTENT(IN OUT) :: aer_sw_scattering_um( npd_profile, npd_layer,          &
                                              n_radaer_mode, n_sw_band )

REAL, INTENT(IN OUT) :: aer_lw_asymmetry_um( npd_profile, npd_layer,           &
                                             n_radaer_mode, n_lw_band )

REAL, INTENT(IN OUT) :: aer_sw_asymmetry_um( npd_profile, npd_layer,           &
                                             n_radaer_mode, n_sw_band )

! Aerosol Optical Depth diagnostics
REAL, INTENT(IN OUT) :: aod_ukca_all_modes_um( npd_profile, npd_ukca_aod_wavel,&
                                               n_ukca_mode )

REAL, INTENT(IN OUT) :: aaod_ukca_all_modes_um(npd_profile, npd_ukca_aod_wavel,&
                                               n_ukca_mode )

! Local variables

! Loop variables
INTEGER :: i, k

! Modal number concentrations (m-3)
REAL ::  ukca_modal_number_um( npd_profile, npd_layer, n_ukca_mode)

! Local AOD diagnostics
REAL ::  aod_ukca_this_mode_um( npd_profile, npd_ukca_aod_wavel )
REAL :: aaod_ukca_this_mode_um( npd_profile, npd_ukca_aod_wavel )
REAL ::  sod_ukca_this_mode_um( npd_profile, npd_ukca_aod_wavel )

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*),   PARAMETER :: RoutineName='UKCA_RADAER_LFRIC_INTERFACE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

  CALL ukca_radaer_prepare(                                                    &
    ! Input Actual array dimensions
    npd_profile, npd_layer, n_ukca_mode, n_ukca_cpnt,                          &
    ! Input Fixed array dimensions
    npd_profile, npd_layer, n_radaer_mode,                                     &
    ! Input from the UKCA_RADAER structure
    nmodes, ncp_max, i_cpnt_index, n_cpnt_in_mode,                             &
    ! Input Component mass-mixing ratios
    ukca_mix_ratio_um,                                                         &
    ! Input modal number concentrations
    ukca_modal_nbr_um,                                                         &
    ! Input Pressure and temperature
    p_theta_levels, t_theta_levels,                                            &
    ! Output Modal mass-mixing ratios
    ukca_mode_mix_ratio_um,                                                    &
    ! Output modal number concentrations
    ukca_modal_number_um                                                       &
    )

  ! Long wave ( e.g. ip_infra_red )
  CALL ukca_radaer_band_average(                                               &
    ! Fixed array dimensions (input)
    npd_profile,                                                               &
    npd_layer,                                                                 &
    n_radaer_mode,                                                             &
    n_lw_band,                                                                 &
    npd_exclude_lw,                                                            &
    ! Spectral information (input)
    n_lw_band,                                                                 &
    ip_infra_red,                                                              &
    l_exclude_lw,                                                              &
    lw_n_band_exclude,                                                         &
    lw_index_exclude,                                                          &
    ! Actual array dimensions (input)
    npd_profile,                                                               &
    npd_layer,                                                                 &
    n_ukca_mode,                                                               &
    n_ukca_cpnt,                                                               &
    ! Prescribed SSA dimensions
    nd_prof_ssa,                                                               &
    nd_layr_ssa,                                                               &
    nd_band_ssa,                                                               &
    ! UKCA_RADAER structure (input)
    nmodes,                                                                    &
    ncp_max,                                                                   &
    ncp_max_x_nmodes,                                                          &
    i_cpnt_index,                                                              &
    i_cpnt_type,                                                               &
    i_mode_type,                                                               &
    l_nitrate,                                                                 &
    l_soluble,                                                                 &
    l_sustrat,                                                                 &
    l_cornarrow_ins,                                                           &
    n_cpnt_in_mode,                                                            &
    ! Modal mass-mixing ratios (input)
    ukca_mode_mix_ratio_um,                                                    &
    ! Modal number concentrations (input)
    ukca_modal_number_um,                                                      &
    ! Modal diameters from UKCA module (input)
    ukca_dry_diam_um,                                                          &
    ukca_wet_diam_um,                                                          &
    ! Other inputs from UKCA module (input)
    ukca_comp_vol_um,                                                          &
    ukca_modal_vol_um,                                                         &
    ukca_modal_rho_um,                                                         &
    ukca_modal_wtv_um,                                                         &
    ! Logical to describe orientation
    l_inverted,                                                                &
    ! Control option for prescribed single scattering albedo array
    i_ukca_radaer_prescribe_ssa,                                               &
    ! Model level of the tropopause (input)
    trindxrad_um,                                                              &
    ! Prescription of single-scattering albedo
    ukca_radaer_presc_ssa,                                                     &
    ! Maxwell-Garnett mixing approach logical control switches
    i_ukca_tune_bc, i_glomap_clim_tune_bc,                                     &
    ! Band-averaged optical properties (output)
    aer_lw_absorption_um,                                                      &
    aer_lw_scattering_um,                                                      &
    aer_lw_asymmetry_um                                                        &
    )

  ! Short wave (e.g. ip_solar ) - only calculate on lit points
  IF ( l_any_lit_points_um ) THEN

    CALL ukca_radaer_band_average(                                             &
        ! Fixed array dimensions (input)
        npd_profile,                                                           &
        npd_layer,                                                             &
        n_radaer_mode,                                                         &
        n_sw_band,                                                             &
        npd_exclude_sw,                                                        &
        ! Spectral information (input)
        n_sw_band,                                                             &
        ip_solar,                                                              &
        l_exclude_sw,                                                          &
        sw_n_band_exclude,                                                     &
        sw_index_exclude,                                                      &
        ! Actual array dimensions (input)
        npd_profile,                                                           &
        npd_layer,                                                             &
        n_ukca_mode,                                                           &
        n_ukca_cpnt,                                                           &
        ! Prescribed SSA dimensions
        nd_prof_ssa,                                                           &
        nd_layr_ssa,                                                           &
        nd_band_ssa,                                                           &
        ! UKCA_RADAER structure (input)
        nmodes,                                                                &
        ncp_max,                                                               &
        ncp_max_x_nmodes,                                                      &
        i_cpnt_index,                                                          &
        i_cpnt_type,                                                           &
        i_mode_type,                                                           &
        l_nitrate,                                                             &
        l_soluble,                                                             &
        l_sustrat,                                                             &
        l_cornarrow_ins,                                                       &
        n_cpnt_in_mode,                                                        &
        ! Modal mass-mixing ratios (input)
        ukca_mode_mix_ratio_um,                                                &
        ! Modal number concentrations (input)
        ukca_modal_number_um,                                                  &
        ! Modal diameters from UKCA module (input)
        ukca_dry_diam_um,                                                      &
        ukca_wet_diam_um,                                                      &
        ! Other inputs from UKCA module (input)
        ukca_comp_vol_um,                                                      &
        ukca_modal_vol_um,                                                     &
        ukca_modal_rho_um,                                                     &
        ukca_modal_wtv_um,                                                     &
        ! Logical to describe orientation
        l_inverted,                                                            &
        ! Switch for prescribed single scattering albedo array
        i_ukca_radaer_prescribe_ssa,                                           &
        ! Model level of the tropopause (input)
        trindxrad_um,                                                          &
        ! Prescription of single-scattering albedo
        ukca_radaer_presc_ssa,                                                 &
        ! Maxwell-Garnett mixing approach logical control switches
        i_ukca_tune_bc, i_glomap_clim_tune_bc,                                 &
        ! Band-averaged optical properties (output)
        aer_sw_absorption_um,                                                  &
        aer_sw_scattering_um,                                                  &
        aer_sw_asymmetry_um )

  END IF

  !------------------------------------------------
  ! Now calculate aod and aaod for Aitken Soluble mode

  IF ( l_aod_ukca_ait_sol .or. l_aaod_ukca_ait_sol ) THEN

    CALL ukca_radaer_compute_aod(                                              &
         ! Fixed array dimensions (input)
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_ukca_mode,                                                          &
         n_ukca_cpnt,                                                          &
         npd_ukca_aod_wavel,                                                   &
         ! Fixed array Prescribed ssa dimensions (input)
         nd_prof_ssa,                                                          &
         nd_layr_ssa,                                                          &
         nd_band_ssa,                                                          &
         ! UKCA_RADAER structure (input)
         nmodes,                                                               &
         ncp_max,                                                              &
         ncp_max_x_nmodes,                                                     &
         i_cpnt_index,                                                         &
         i_cpnt_type,                                                          &
         n_cpnt_in_mode,                                                       &
         l_nitrate,                                                            &
         l_soluble,                                                            &
         l_sustrat,                                                            &
         i_mode_type,                                                          &
         l_cornarrow_ins,                                                      &
         ! Modal diameters from UKCA module
         ukca_dry_diam_um,                                                     &
         ukca_wet_diam_um,                                                     &
         ! Mass thickness of layers
         d_mass_theta_levels_um,                                               &
         ! Component volumes
         ukca_comp_vol_um,                                                     &
         ! Modal volumes, densities, and water content
         ukca_modal_vol_um,                                                    &
         ukca_modal_rho_um,                                                    &
         ukca_modal_wtv_um,                                                    &
         ! Modal mass-mixing ratios
         ukca_mode_mix_ratio_um,                                               &
         ! Modal number concentrations
         ukca_modal_number_um,                                                 &
         ! Type selection
         ip_ukca_mode_aitken,                                                  &
         soluble_wanted,                                                       &
         ! Switch for if prescribed SSA is on
         i_ukca_radaer_prescribe_ssa,                                          &
         ! Model level of the tropopause
         trindxrad_um,                                                         &
         ! Prescription of single-scattering albedo
         ukca_radaer_presc_ssa,                                                &
         ! Modal extinction aerosol opt depth - column (output)
         aod_ukca_this_mode_um,                                                &
         ! Modal extinction aerosol opt depth - stratosphere (output)
         sod_ukca_this_mode_um,                                                &
         ! Modal absorption aerosol opt depth (output)
         aaod_ukca_this_mode_um,                                               &
         ! Fixed array dimensions
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_radaer_mode,                                                        &
         npd_ukca_aod_wavel )

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aod_ukca_all_modes_um(i,k,mode_ait_sol-1) = aod_ukca_this_mode_um(i,k)
      END DO
    END DO

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aaod_ukca_all_modes_um(i,k,mode_ait_sol-1) = aaod_ukca_this_mode_um(i,k)
      END DO
    END DO

  END IF ! Calculate AOD Aitken Soluble mode

  !------------------------------------------------
  ! Now calculate aod and aaod for Accumulation Soluble mode

  IF ( l_aod_ukca_acc_sol .OR. l_aaod_ukca_acc_sol ) THEN

    CALL ukca_radaer_compute_aod(                                              &
         ! Fixed array dimensions (input)
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_ukca_mode,                                                          &
         n_ukca_cpnt,                                                          &
         npd_ukca_aod_wavel,                                                   &
         ! Fixed array Prescribed ssa dimensions (input)
         nd_prof_ssa,                                                          &
         nd_layr_ssa,                                                          &
         nd_band_ssa,                                                          &
         ! UKCA_RADAER structure (input)
         nmodes,                                                               &
         ncp_max,                                                              &
         ncp_max_x_nmodes,                                                     &
         i_cpnt_index,                                                         &
         i_cpnt_type,                                                          &
         n_cpnt_in_mode,                                                       &
         l_nitrate,                                                            &
         l_soluble,                                                            &
         l_sustrat,                                                            &
         i_mode_type,                                                          &
         l_cornarrow_ins,                                                      &
         ! Modal diameters from UKCA module
         ukca_dry_diam_um,                                                     &
         ukca_wet_diam_um,                                                     &
         ! Mass thickness of layers
         d_mass_theta_levels_um,                                               &
         ! Component volumes
         ukca_comp_vol_um,                                                     &
         ! Modal volumes, densities, and water content
         ukca_modal_vol_um,                                                    &
         ukca_modal_rho_um,                                                    &
         ukca_modal_wtv_um,                                                    &
         ! Modal mass-mixing ratios
         ukca_mode_mix_ratio_um,                                               &
         ! Modal number concentrations
         ukca_modal_number_um,                                                 &
         ! Type selection
         ip_ukca_mode_accum,                                                   &
         soluble_wanted,                                                       &
         ! Switch for if prescribed SSA is on
         i_ukca_radaer_prescribe_ssa,                                          &
         ! Model level of the tropopause
         trindxrad_um,                                                         &
         ! Prescription of single-scattering albedo
         ukca_radaer_presc_ssa,                                                &
         ! Modal extinction aerosol opt depth - column (output)
         aod_ukca_this_mode_um,                                                &
         ! Modal extinction aerosol opt depth - stratosphere (output)
         sod_ukca_this_mode_um,                                                &
         ! Modal absorption aerosol opt depth (output)
         aaod_ukca_this_mode_um,                                               &
         ! Fixed array dimensions
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_radaer_mode,                                                        &
         npd_ukca_aod_wavel )

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aod_ukca_all_modes_um(i,k,mode_acc_sol-1) = aod_ukca_this_mode_um(i,k)
      END DO
    END DO

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aaod_ukca_all_modes_um(i,k,mode_acc_sol-1) = aaod_ukca_this_mode_um(i,k)
      END DO
    END DO

  END IF ! Calculate AOD Accumulation Soluble mode

  !------------------------------------------------
  ! Now calculate aod and aaod for Coarse Soluble mode

  IF ( l_aod_ukca_cor_sol .OR. l_aaod_ukca_cor_sol ) THEN

    CALL ukca_radaer_compute_aod(                                              &
         ! Fixed array dimensions (input)
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_ukca_mode,                                                          &
         n_ukca_cpnt,                                                          &
         npd_ukca_aod_wavel,                                                   &
         ! Fixed array Prescribed ssa dimensions (input)
         nd_prof_ssa,                                                          &
         nd_layr_ssa,                                                          &
         nd_band_ssa,                                                          &
         ! UKCA_RADAER structure (input)
         nmodes,                                                               &
         ncp_max,                                                              &
         ncp_max_x_nmodes,                                                     &
         i_cpnt_index,                                                         &
         i_cpnt_type,                                                          &
         n_cpnt_in_mode,                                                       &
         l_nitrate,                                                            &
         l_soluble,                                                            &
         l_sustrat,                                                            &
         i_mode_type,                                                          &
         l_cornarrow_ins,                                                      &
         ! Modal diameters from UKCA module
         ukca_dry_diam_um,                                                     &
         ukca_wet_diam_um,                                                     &
         ! Mass thickness of layers
         d_mass_theta_levels_um,                                               &
         ! Component volumes
         ukca_comp_vol_um,                                                     &
         ! Modal volumes, densities, and water content
         ukca_modal_vol_um,                                                    &
         ukca_modal_rho_um,                                                    &
         ukca_modal_wtv_um,                                                    &
         ! Modal mass-mixing ratios
         ukca_mode_mix_ratio_um,                                               &
         ! Modal number concentrations
         ukca_modal_number_um,                                                 &
         ! Type selection
         ip_ukca_mode_coarse,                                                  &
         soluble_wanted,                                                       &
         ! Switch for if prescribed SSA is on
         i_ukca_radaer_prescribe_ssa,                                          &
         ! Model level of the tropopause
         trindxrad_um,                                                         &
         ! Prescription of single-scattering albedo
         ukca_radaer_presc_ssa,                                                &
         ! Modal extinction aerosol opt depth - column (output)
         aod_ukca_this_mode_um,                                                &
         ! Modal extinction aerosol opt depth - stratosphere (output)
         sod_ukca_this_mode_um,                                                &
         ! Modal absorption aerosol opt depth (output)
         aaod_ukca_this_mode_um,                                               &
         ! Fixed array dimensions
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_radaer_mode,                                                        &
         npd_ukca_aod_wavel )

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aod_ukca_all_modes_um(i,k,mode_cor_sol-1) = aod_ukca_this_mode_um(i,k)
      END DO
    END DO

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aaod_ukca_all_modes_um(i,k,mode_cor_sol-1) = aaod_ukca_this_mode_um(i,k)
      END DO
    END DO

  END IF ! Calculate AOD Coarse Soluble mode

  !------------------------------------------------
  ! Now calculate aod and aaod for Aitken Insoluble mode
 
  IF ( l_aod_ukca_ait_ins .OR. l_aaod_ukca_ait_ins ) THEN

    CALL ukca_radaer_compute_aod(                                              &
         ! Fixed array dimensions (input)
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_ukca_mode,                                                          &
         n_ukca_cpnt,                                                          &
         npd_ukca_aod_wavel,                                                   &
         ! Fixed array Prescribed ssa dimensions (input)
         nd_prof_ssa,                                                          &
         nd_layr_ssa,                                                          &
         nd_band_ssa,                                                          &
         ! UKCA_RADAER structure (input)
         nmodes,                                                               &
         ncp_max,                                                              &
         ncp_max_x_nmodes,                                                     &
         i_cpnt_index,                                                         &
         i_cpnt_type,                                                          &
         n_cpnt_in_mode,                                                       &
         l_nitrate,                                                            &
         l_soluble,                                                            &
         l_sustrat,                                                            &
         i_mode_type,                                                          &
         l_cornarrow_ins,                                                      &
         ! Modal diameters from UKCA module
         ukca_dry_diam_um,                                                     &
         ukca_wet_diam_um,                                                     &
         ! Mass thickness of layers
         d_mass_theta_levels_um,                                               &
         ! Component volumes
         ukca_comp_vol_um,                                                     &
         ! Modal volumes, densities, and water content
         ukca_modal_vol_um,                                                    &
         ukca_modal_rho_um,                                                    &
         ukca_modal_wtv_um,                                                    &
         ! Modal mass-mixing ratios
         ukca_mode_mix_ratio_um,                                               &
         ! Modal number concentrations
         ukca_modal_number_um,                                                 &
         ! Type selection
         ip_ukca_mode_aitken,                                                  &
         soluble_unwanted,                                                     &
         ! Switch for if prescribed SSA is on
         i_ukca_radaer_prescribe_ssa,                                          &
         ! Model level of the tropopause
         trindxrad_um,                                                         &
         ! Prescription of single-scattering albedo
         ukca_radaer_presc_ssa,                                                &
         ! Modal extinction aerosol opt depth - column (output)
         aod_ukca_this_mode_um,                                                &
         ! Modal extinction aerosol opt depth - stratosphere (output)
         sod_ukca_this_mode_um,                                                &
         ! Modal absorption aerosol opt depth (output)
         aaod_ukca_this_mode_um,                                               &
         ! Fixed array dimensions
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_radaer_mode,                                                        &
         npd_ukca_aod_wavel )

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aod_ukca_all_modes_um(i,k,mode_ait_insol-1) = aod_ukca_this_mode_um(i,k)
      END DO
    END DO

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aaod_ukca_all_modes_um(i,k,mode_ait_insol-1)=aaod_ukca_this_mode_um(i,k)
      END DO
    END DO

  END IF ! Calculate AOD Aitkin Insoluble mode

  !------------------------------------------------
  ! Now calculate aod and aaod for Accumulation Insoluble mode

  IF ( l_aod_ukca_acc_ins .OR. l_aaod_ukca_acc_ins ) THEN

    call ukca_radaer_compute_aod(                                              &
         ! Fixed array dimensions (input)
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_ukca_mode,                                                          &
         n_ukca_cpnt,                                                          &
         npd_ukca_aod_wavel,                                                   &
         ! Fixed array Prescribed ssa dimensions (input)
         nd_prof_ssa,                                                          &
         nd_layr_ssa,                                                          &
         nd_band_ssa,                                                          &
         ! UKCA_RADAER structure (input)
         nmodes,                                                               &
         ncp_max,                                                              &
         ncp_max_x_nmodes,                                                     &
         i_cpnt_index,                                                         &
         i_cpnt_type,                                                          &
         n_cpnt_in_mode,                                                       &
         l_nitrate,                                                            &
         l_soluble,                                                            &
         l_sustrat,                                                            &
         i_mode_type,                                                          &
         l_cornarrow_ins,                                                      &
         ! Modal diameters from UKCA module
         ukca_dry_diam_um,                                                     &
         ukca_wet_diam_um,                                                     &
         ! Mass thickness of layers
         d_mass_theta_levels_um,                                               &
         ! Component volumes
         ukca_comp_vol_um,                                                     &
         ! Modal volumes, densities, and water content
         ukca_modal_vol_um,                                                    &
         ukca_modal_rho_um,                                                    &
         ukca_modal_wtv_um,                                                    &
         ! Modal mass-mixing ratios
         ukca_mode_mix_ratio_um,                                               &
         ! Modal number concentrations
         ukca_modal_number_um,                                                 &
         ! Type selection
         ip_ukca_mode_accum,                                                   &
         soluble_unwanted,                                                     &
         ! Switch for if prescribed SSA is on
         i_ukca_radaer_prescribe_ssa,                                          &
         ! Model level of the tropopause
         trindxrad_um,                                                         &
         ! Prescription of single-scattering albedo
         ukca_radaer_presc_ssa,                                                &
         ! Modal extinction aerosol opt depth - column (output)
         aod_ukca_this_mode_um,                                                &
         ! Modal extinction aerosol opt depth - stratosphere (output)
         sod_ukca_this_mode_um,                                                &
         ! Modal absorption aerosol opt depth (output)
         aaod_ukca_this_mode_um,                                               &
         ! Fixed array dimensions
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_radaer_mode,                                                        &
         npd_ukca_aod_wavel )

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aod_ukca_all_modes_um(i,k,mode_acc_insol-1) = aod_ukca_this_mode_um(i,k)
      END DO
    END DO

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aaod_ukca_all_modes_um(i,k,mode_acc_insol-1)=aaod_ukca_this_mode_um(i,k)
      END DO
    END DO

  END IF ! Calculate AOD Accumulation Insoluble mode

  !------------------------------------------------
  ! Now calculate aod and aaod for Coarse Insoluble mode
  IF ( l_aod_ukca_cor_ins .OR. l_aaod_ukca_cor_ins ) THEN

    call ukca_radaer_compute_aod(                                              &
         ! Fixed array dimensions (input)
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_ukca_mode,                                                          &
         n_ukca_cpnt,                                                          &
         npd_ukca_aod_wavel,                                                   &
         ! Fixed array Prescribed ssa dimensions (input)
         nd_prof_ssa,                                                          &
         nd_layr_ssa,                                                          &
         nd_band_ssa,                                                          &
         ! UKCA_RADAER structure (input)
         nmodes,                                                               &
         ncp_max,                                                              &
         ncp_max_x_nmodes,                                                     &
         i_cpnt_index,                                                         &
         i_cpnt_type,                                                          &
         n_cpnt_in_mode,                                                       &
         l_nitrate,                                                            &
         l_soluble,                                                            &
         l_sustrat,                                                            &
         i_mode_type,                                                          &
         l_cornarrow_ins,                                                      &
         ! Modal diameters from UKCA module
         ukca_dry_diam_um,                                                     &
         ukca_wet_diam_um,                                                     &
         ! Mass thickness of layers
         d_mass_theta_levels_um,                                               &
         ! Component volumes
         ukca_comp_vol_um,                                                     &
         ! Modal volumes, densities, and water content
         ukca_modal_vol_um,                                                    &
         ukca_modal_rho_um,                                                    &
         ukca_modal_wtv_um,                                                    &
         ! Modal mass-mixing ratios
         ukca_mode_mix_ratio_um,                                               &
         ! Modal number concentrations
         ukca_modal_number_um,                                                 &
         ! Type selection
         ip_ukca_mode_coarse,                                                  &
         soluble_unwanted,                                                     &
         ! Switch for if prescribed SSA is on
         i_ukca_radaer_prescribe_ssa,                                          &
         ! Model level of the tropopause
         trindxrad_um,                                                         &
         ! Prescription of single-scattering albedo
         ukca_radaer_presc_ssa,                                                &
         ! Modal extinction aerosol opt depth - column (output)
         aod_ukca_this_mode_um,                                                &
         ! Modal extinction aerosol opt depth - stratosphere (output)
         sod_ukca_this_mode_um,                                                &
         ! Modal absorption aerosol opt depth (output)
         aaod_ukca_this_mode_um,                                               &
         ! Fixed array dimensions
         npd_profile,                                                          &
         npd_layer,                                                            &
         n_radaer_mode,                                                        &
         npd_ukca_aod_wavel )

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aod_ukca_all_modes_um(i,k,mode_cor_insol-1) = aod_ukca_this_mode_um(i,k)
      END DO
    END DO

    DO k = 1, npd_ukca_aod_wavel
      DO i = 1, npd_profile
        aaod_ukca_all_modes_um(i,k,mode_cor_insol-1)=aaod_ukca_this_mode_um(i,k)
      END DO
    END DO

  END IF ! Calculate AOD Coarse Insoluble mode

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE ukca_radaer_lfric_interface

END MODULE ukca_radaer_lfric_mod
