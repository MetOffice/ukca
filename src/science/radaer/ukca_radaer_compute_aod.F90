! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Obtain UKCA-MODE aerosol optical properties at optical depth
!  wavelengths using UKCA_RADAER look-up tables, and compute the
!  UKCA-MODE modal optical depth for the mode requested.
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
MODULE ukca_radaer_compute_aod_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName = 'UKCA_RADAER_COMPUTE_AOD_MOD'

CONTAINS

SUBROUTINE ukca_radaer_compute_aod(                                            &
      ! Actual array dimensions
      n_profile, n_layer, n_ukca_mode, n_ukca_cpnt, n_aod_wavel,               &
      ! Prescribed ssa dimensions (Fixed array)
      npd_prof_ssa, npd_layr_ssa, npd_naod_ssa,                                &
      ! UKCA_RADAER structure
      ! From the structure ukca_radaer for UKCA/radiation interaction
      nmodes, ncp_max, ncp_max_x_nmodes,                                       &
      i_cpnt_index, i_cpnt_type, n_cpnt_in_mode,                               &
      l_nitrate, l_soluble, l_sustrat, i_mode_type, l_cornarrow_ins,           &
      ! Modal diameters from UKCA module
      ukca_dry, ukca_wet,                                                      &
      ! Mass thickness of layers
      d_mass,                                                                  &
      ! Component volumes
      ukca_cpnt_volume,                                                        &
      ! Modal volumes, densities, and water content
      ukca_modal_volume, ukca_modal_density, ukca_water_volume,                &
      ! Modal mass-mixing ratios
      ukca_modal_mmr,                                                          &
      ! Modal number concentrations
      ukca_modal_number,                                                       &
      ! Type selection
      type_wanted, soluble_wanted,                                             &
      ! Logical for if prescribed SSA is on
      l_ukca_radaer_prescribe_ssa,                                             &
      ! Model level of the tropopause
      trindxrad,                                                               &
      ! Prescription of single-scattering albedo
      ukca_radaer_presc_ssa_aod,                                               &
      ! Modal extinction aerosol optical depths:
      !   full column and stratosphere
      ukca_modal_aod, ukca_modal_strat_aod,                                    &
      ! Modal absorption aerosol optical depth:  full column
      ukca_modal_aaod,                                                         &
      ! Fixed array dimensions
      npd_profile,  npd_layer, npd_aerosol_mode, npd_aod_wavel                 &
   )

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE conversions_mod, ONLY: pi
USE ukca_radaer_lut_read_in, ONLY: ukca_radaer_get_lut_index

! UKCA look-up tables
USE ukca_radaer_lut,        ONLY:                                              &
    ip_ukca_lut_accum,                                                         &
    ip_ukca_lut_coarse,                                                        &
    ip_ukca_lut_accnarrow,                                                     &
    ip_ukca_lut_cornarrow,                                                     &
    ip_ukca_lut_supercoarse,                                                   &
    ip_ukca_lut_sw,                                                            &
    ukca_lut

! UKCA pre-computed values
USE ukca_radaer_precalc,    ONLY:                                              &
    precalc

USE ukca_mode_setup,        ONLY:                                              &
    ip_ukca_mode_aitken,                                                       &
    ip_ukca_mode_accum,                                                        &
    ip_ukca_mode_coarse,                                                       &
    ip_ukca_mode_supercoarse

USE ukca_radaer_struct_mod, ONLY:                                              &
    ukca_radaer_struct,                                                        &
    threshold_mmr, threshold_vol,                                              &
    threshold_nbr

USE ukca_radaer_ri_calc_mod, ONLY:                                             &
    ukca_radaer_ri_calc

USE ukca_option_mod,         ONLY:                                             &
    i_ukca_tune_bc

USE glomap_clim_option_mod,  ONLY:                                             &
    i_glomap_clim_tune_bc

IMPLICIT NONE

!
! Arguments with intent in
!

! Fixed array dimensions
INTEGER, INTENT(IN) :: npd_profile,                                            &
        npd_layer,                                                             &
        npd_aerosol_mode,                                                      &
        npd_aod_wavel

! Actual array dimensions
INTEGER, INTENT(IN) :: n_profile,                                              &
        n_layer,                                                               &
        n_ukca_mode,                                                           &
        n_ukca_cpnt,                                                           &
        n_aod_wavel

!
! Fixed array dimensions for prescribed SSA
!
INTEGER, INTENT(IN) :: npd_prof_ssa,                                           &
        npd_layr_ssa,                                                          &
        npd_naod_ssa

! From radaer structure
INTEGER, INTENT(IN) :: nmodes
INTEGER, INTENT(IN) :: ncp_max
INTEGER, INTENT(IN) :: ncp_max_x_nmodes
INTEGER, INTENT(IN) :: i_cpnt_index( ncp_max, nmodes )
INTEGER, INTENT(IN) :: i_cpnt_type( ncp_max_x_nmodes )
INTEGER, INTENT(IN) :: n_cpnt_in_mode( nmodes )
LOGICAL, INTENT(IN) :: l_nitrate
LOGICAL, INTENT(IN) :: l_soluble( nmodes )
LOGICAL, INTENT(IN) :: l_sustrat
INTEGER, INTENT(IN) :: i_mode_type( nmodes )
LOGICAL, INTENT(IN) :: l_cornarrow_ins

! Modal dry and wet diameters
REAL, INTENT(IN) :: ukca_dry(npd_profile, npd_layer, n_ukca_mode)
REAL, INTENT(IN) :: ukca_wet(npd_profile, npd_layer, n_ukca_mode)

! Mass-thickness of vertical layers
REAL, INTENT(IN) :: d_mass(npd_profile, npd_layer)

! Component volumes
REAL, INTENT(IN) :: ukca_cpnt_volume(n_ukca_cpnt, npd_profile, npd_layer)

! Modal volumes
REAL, INTENT(IN) :: ukca_modal_volume(npd_profile, npd_layer, n_ukca_mode)

! Modal density
REAL, INTENT(IN) :: ukca_modal_density(npd_profile, npd_layer, n_ukca_mode)

! Volume of water in each mode
REAL, INTENT(IN) :: ukca_water_volume(npd_profile, npd_layer, n_ukca_mode)

! Modal mass-mixing ratios
REAL, INTENT(IN) :: ukca_modal_mmr(npd_profile, npd_layer, npd_aerosol_mode)

! Modal number concentrations
REAL, INTENT(IN) :: ukca_modal_number(npd_profile, npd_layer, n_ukca_mode)

! Type selection
INTEGER, INTENT(IN) :: type_wanted
LOGICAL, INTENT(IN) :: soluble_wanted

!
! When true, use a prescribed single scattering albedo field
!
LOGICAL, INTENT(IN) :: l_ukca_radaer_prescribe_ssa

! Model level of the tropopause
INTEGER, INTENT(IN) :: trindxrad(npd_profile)

! Prescription of single-scattering albedo
REAL, INTENT(IN) :: ukca_radaer_presc_ssa_aod( npd_prof_ssa,                   &
                                               npd_layr_ssa,                   &
                                               npd_naod_ssa )

!
! Arguments with intent out
!

! Modal aerosol optical depths: full column and stratosphere
REAL , INTENT(OUT):: ukca_modal_aod       (npd_profile, npd_aod_wavel),        &
                     ukca_modal_aaod      (npd_profile, npd_aod_wavel),        &
                     ukca_modal_strat_aod (npd_profile, npd_aod_wavel)
!
! Local variables
!
INTEGER :: i_mode,                                                             &
        i_layr,                                                                &
        i_prof,                                                                &
        n

!
! Values at the AOD wavelengths:
!      Mie parameter for the wet and dry diameters and
!      the indices of their nearest neighbour
!      Complex refractive index and the index of its nearest neighbour
!
REAL :: x
INTEGER :: n_x
REAL :: x_dry
INTEGER :: n_x_dry
REAL :: re_m
INTEGER :: n_nr
REAL :: im_m
INTEGER :: n_ni
INTEGER, PARAMETER :: n_ni_fix = 1

!
! Values at given wavelength
!
REAL :: loc_abs
REAL :: loc_sca
REAL :: loc_vol
REAL :: factor

!
! Single-scattering albedo to prescribe
!
REAL :: this_ssa

!
! Local copies of typedef members
!
INTEGER :: nx
REAL :: logxmin         ! log(xmin)
REAL :: logxmaxmlogxmin ! log(xmax) - log(xmin)
INTEGER :: nnr
REAL :: nrmin
REAL :: incr_nr
INTEGER :: nni
REAL :: ni_min
REAL :: ni_max
REAL :: ni_c

!
! Local copies of mode type, component index and component type
!
INTEGER :: this_mode_type

!
! Thresholds on the modal mass-mixing ratio, volume, and modal number
! concentrations above which aerosol optical properties are to be
! computed - threshold_mmr, threshold_vol, threshold_nbr - specified
! in ukca_radaer_struct_mod
!
! Indicates whether current level is above the tropopause.
!
LOGICAL :: l_in_stratosphere

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_COMPUTE_AOD'

!
!
!

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

DO n = 1, npd_aod_wavel
  DO i_prof = 1, npd_profile
    ukca_modal_aod(       i_prof, n ) = 0.0e+00
    ukca_modal_aaod(      i_prof, n ) = 0.0e+00
    ukca_modal_strat_aod( i_prof, n ) = 0.0e+00
  END DO
END DO

DO i_mode = 1, n_ukca_mode

  IF ( (i_mode_type(i_mode) == type_wanted) .AND.                              &
       (l_soluble(i_mode) .EQV. soluble_wanted) ) THEN

    !
    ! Mode type. From a look-up table point of view, Aitken and
    ! accumulation types are treated in the same way.
    ! Accumulation soluble mode may use a narrower width (i.e. another
    ! look-up table) than other Aitken and accumulation modes.
    ! The coarse insoluble mode in the 3 dust mode setup may have a
    ! narrower width than the default 2.0 which is the case if the
    ! super-coarse insoluble mode is selected
    ! Once we know which look-up table to select, make local copies
    ! of info needed for nearest-neighbour calculations.
    ! For aerosol optical depth calculations, the shortwave spectrum
    ! look-up tables are used.
    !

    SELECT CASE (i_mode_type(i_mode))

    CASE (ip_ukca_mode_aitken)
      this_mode_type = ip_ukca_lut_accum

    CASE (ip_ukca_mode_accum)
      IF (l_soluble(i_mode)) THEN
        this_mode_type = ip_ukca_lut_accnarrow
      ELSE
        this_mode_type = ip_ukca_lut_accum
      END IF

    CASE (ip_ukca_mode_coarse)
      IF ( (.NOT. l_soluble(i_mode)) .AND. l_cornarrow_ins) THEN
        this_mode_type = ip_ukca_lut_cornarrow
      ELSE
        this_mode_type = ip_ukca_lut_coarse
      END IF

    CASE (ip_ukca_mode_supercoarse)
      this_mode_type = ip_ukca_lut_supercoarse

    END SELECT

    nx      = ukca_lut(this_mode_type, ip_ukca_lut_sw)%n_x
    logxmin = LOG(ukca_lut(this_mode_type, ip_ukca_lut_sw)%x_min)
    logxmaxmlogxmin =                                                          &
            LOG(ukca_lut(this_mode_type, ip_ukca_lut_sw)%x_max) - logxmin

    nnr     = ukca_lut(this_mode_type, ip_ukca_lut_sw)%n_nr
    nrmin   = ukca_lut(this_mode_type, ip_ukca_lut_sw)%nr_min
    incr_nr = ukca_lut(this_mode_type, ip_ukca_lut_sw)%incr_nr

    nni     = ukca_lut(this_mode_type, ip_ukca_lut_sw)%n_ni
    ni_min  = ukca_lut(this_mode_type, ip_ukca_lut_sw)%ni_min
    ni_max  = ukca_lut(this_mode_type, ip_ukca_lut_sw)%ni_max
    ni_c    = ukca_lut(this_mode_type, ip_ukca_lut_sw)%ni_c

    DO i_layr = 1, n_layer

      DO i_prof = 1, n_profile

        l_in_stratosphere = i_layr <= trindxrad(i_prof)

        ! Only make calculations if there is significant aerosol mass, and
        ! if the number concentration is large enough.
        ! This test is especially important for the first timestep,
        ! as UKCA has not run yet and its output is therefore
        ! not guaranteed to be valid. Mass mixing ratios and numbers
        ! are initialised to zero as prognostics.
        ! Also, at low number concentrations, the size informations
        ! given by UKCA are unreliable and might produce erroneous
        ! optical properties.
        !
        ! The threshold on ukca_modal_volume is a way of ensuring
        ! that UKCA-mode has actually been called
        ! (ukca_modal_volume will be zero by default first time step)
        !

        IF (ukca_modal_mmr   (i_prof, i_layr, i_mode) > threshold_mmr  .AND.   &
            ukca_modal_number(i_prof, i_layr, i_mode) > threshold_nbr  .AND.   &
            ukca_modal_volume(i_prof, i_layr, i_mode) > threshold_vol) THEN

          DO n = 1, n_aod_wavel

            x = pi * ukca_wet(i_prof, i_layr, i_mode) / precalc%aod_wavel(n)
            n_x = NINT( (LOG(x)    - logxmin) /                                &
                         logxmaxmlogxmin * (nx-1) ) + 1
            n_x = MIN(nx, MAX(1, n_x))

            x_dry = pi * ukca_dry(i_prof, i_layr, i_mode)                      &
                             / precalc%aod_wavel(n)
            n_x_dry = NINT( (LOG(x_dry) - logxmin) /                           &
                             logxmaxmlogxmin * (nx-1) ) + 1
            n_x_dry = MIN(nx, MAX(1, n_x_dry))

            ! Compute the modal complex refractive index via
            ! volume-weighting for non-BC components. If i_glomap_clim_tune_bc
            ! or i_ukca_tune_bc are set to i_ukca_bc_mg_mix the BC component
            ! will be incorporated via the Maxwell-Garnett mixing approach
            ! otherwise use volume-weighting for BC also.
            !
            CALL ukca_radaer_ri_calc(                                          &
              ! From the structure ukca_radaer for UKCA/radiation interaction
              nmodes,                                                          &
              ncp_max,                                                         &
              ncp_max_x_nmodes,                                                &
              i_cpnt_index,                                                    &
              i_cpnt_type,                                                     &
              n_cpnt_in_mode,                                                  &
              l_nitrate,                                                       &
              l_soluble,                                                       &
              l_sustrat,                                                       &
              ! Refractive index
              precalc%aod_realrefr( :, n ),                                    &
              precalc%aod_imagrefr( :, n ),                                    &
              ! Modal properties
              ukca_cpnt_volume( :, i_prof, i_layr ),                           &
              ukca_modal_volume( i_prof, i_layr, i_mode ),                     &
              ukca_water_volume( i_prof, i_layr, i_mode ),                     &
              ! Indicies for arrays
              i_mode, n_ukca_cpnt,                                             &
              ! Stratospheric aerosol treated as sulphuric acid?
              l_in_stratosphere,                                               &
              ! Logical control switches
              i_ukca_tune_bc, i_glomap_clim_tune_bc,                           &
              l_ukca_radaer_prescribe_ssa,                                     &
              ! Output refractive index real and imag parts
              re_m, im_m )

            n_nr = NINT( (re_m - nrmin) / incr_nr ) + 1
            n_nr = MIN(nnr, MAX(1, n_nr))

            ! The AOD calculations depend on whether SSA is prescribed
            IF (l_ukca_radaer_prescribe_ssa) THEN

              ! Fix the imaginary index to 1 (no absorption) as absorptivity
              ! is prescribed
              n_ni = n_ni_fix

              loc_sca = ukca_lut(this_mode_type, ip_ukca_lut_sw)%              &
                        ukca_scattering(n_x, n_ni, n_nr)

              loc_vol = ukca_lut(this_mode_type, ip_ukca_lut_sw)%              &
                        volume_fraction(n_x_dry)

              factor = (ukca_modal_density(i_prof, i_layr, i_mode) * loc_vol * &
                        precalc%aod_wavel(n))

              loc_sca = loc_sca / factor

              !
              ! The single-scattering albedo is prescribed by distributing
              ! extinction (which is equal to scattering in the non-absorption
              ! case) to absorption and scattering coefficients in the
              ! proportion indicated by the prescription.
              !
              this_ssa = ukca_radaer_presc_ssa_aod(i_prof, i_layr, n)

              !
              ! aerosol optical depth for extinction with prescribed SSA
              !
              ukca_modal_aod(i_prof, n) =                                      &
                ukca_modal_aod(i_prof, n) +                                    &
                ukca_modal_mmr(i_prof, i_layr, i_mode) *                       &
                d_mass(i_prof, i_layr) * loc_sca

              !
              ! absorption aerosol optical depth with prescribed SSA
              !
              ukca_modal_aaod(i_prof, n) =                                     &
                ukca_modal_aaod(i_prof, n) +                                   &
                ukca_modal_mmr(i_prof, i_layr, i_mode) *                       &
                d_mass(i_prof, i_layr) * loc_sca * (1.0 - this_ssa)

              !
              ! aerosol optical depth for extinction in the stratosphere
              ! with prescribed SSA
              !
              IF (l_in_stratosphere) THEN

                ukca_modal_strat_aod(i_prof, n) =                              &
                  ukca_modal_strat_aod(i_prof, n) +                            &
                  ukca_modal_mmr(i_prof, i_layr, i_mode) *                     &
                  d_mass(i_prof, i_layr) * loc_sca

              END IF

            ELSE

              CALL ukca_radaer_get_lut_index(                                  &
                   nni, im_m, ni_min, ni_max, ni_c, n_ni)

              loc_abs = ukca_lut(this_mode_type, ip_ukca_lut_sw)%              &
                        ukca_absorption(n_x, n_ni, n_nr)

              loc_sca = ukca_lut(this_mode_type, ip_ukca_lut_sw)%              &
                        ukca_scattering(n_x, n_ni, n_nr)

              loc_vol = ukca_lut(this_mode_type, ip_ukca_lut_sw)%              &
                        volume_fraction(n_x_dry)

              factor = (ukca_modal_density(i_prof, i_layr, i_mode) * loc_vol * &
                        precalc%aod_wavel(n))

              loc_abs = loc_abs / factor
              loc_sca = loc_sca / factor

              !
              ! aerosol optical depth for extinction
              !
              ukca_modal_aod(i_prof, n) =                                      &
              ukca_modal_aod(i_prof, n) +                                      &
              ukca_modal_mmr(i_prof, i_layr, i_mode) *                         &
              d_mass(i_prof, i_layr) * (loc_abs + loc_sca)

              !
              ! absorption aerosol optical depth
              !
              ukca_modal_aaod(i_prof, n) =                                     &
              ukca_modal_aaod(i_prof, n) +                                     &
              ukca_modal_mmr(i_prof, i_layr, i_mode) *                         &
              d_mass(i_prof, i_layr) * loc_abs

              !
              ! aerosol optical depth for extinction in the stratosphere
              !
              IF (l_in_stratosphere) THEN

                ukca_modal_strat_aod(i_prof, n) =                              &
                ukca_modal_strat_aod(i_prof, n) +                              &
                ukca_modal_mmr(i_prof, i_layr, i_mode) *                       &
                d_mass(i_prof, i_layr) * (loc_abs + loc_sca)

              END IF

            END IF ! IF (l_ukca_radaer_prescribe_ssa)

          END DO ! n

        END IF

      END DO ! i_prof

    END DO ! i_layr

  END IF

END DO ! i_mode

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_compute_aod
END MODULE ukca_radaer_compute_aod_mod
