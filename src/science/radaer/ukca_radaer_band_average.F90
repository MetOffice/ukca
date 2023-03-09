! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Average optical properties of UKCA-MODE aerosols, as obtained from
!  look-up tables, over spectral wavebands.
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
MODULE ukca_radaer_band_average_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName = 'UKCA_RADAER_BAND_AVERAGE_MOD'

CONTAINS

SUBROUTINE ukca_radaer_band_average(                                           &
      ! Fixed array dimensions
      npd_profile, npd_layer, npd_aerosol_mode, npd_band, npd_exclude          &
      ! Spectral information
   ,  n_band, isolir, l_exclude, n_band_exclude, index_exclude                 &
      ! Actual array dimensions
   ,  n_profile, n_layer, n_ukca_mode, n_ukca_cpnt                             &
      ! Prescribed ssa dimensions (Fixed array)
   ,  npd_prof_ssa, npd_layr_ssa, npd_band_ssa                                 &
      ! From the structure ukca_radaer for UKCA/radiation interaction
   ,  nmodes                                                                   &
   ,  ncp_max                                                                  &
   ,  ncp_max_x_nmodes                                                         &
   ,  i_cpnt_index                                                             &
   ,  i_cpnt_type                                                              &
   ,  i_mode_type                                                              &
   ,  l_nitrate                                                                &
   ,  l_soluble                                                                &
   ,  l_sustrat                                                                &
   ,  n_cpnt_in_mode                                                           &
      ! Modal mass-mixing ratios from UKCA module
   ,  ukca_modal_mmr                                                           &
      ! Modal number concentrations from UKCA module
   ,  ukca_modal_number                                                        &
      ! Modal diameters from UKCA module
   ,  ukca_dry_diam, ukca_wet_diam                                             &
      ! Other inputs from UKCA module
   ,  ukca_cpnt_volume, ukca_modal_volume, ukca_modal_density                  &
   ,  ukca_water_volume                                                        &
      ! Logical to describe orientation
   ,  l_inverted                                                               &
      ! Logical for prescribed single scattering albedo array
   ,  l_ukca_radaer_prescribe_ssa                                              &
      ! Model level of tropopause
   ,  trindxrad                                                                &
      ! Prescription of single-scattering albedo
   ,  ukca_radaer_presc_ssa                                                    &
      ! Logical control switches
   ,  i_ukca_tune_bc, i_glomap_clim_tune_bc                                    &
      ! Band-averaged optical properties (outputs)
   ,  ukca_absorption, ukca_scattering, ukca_asymmetry                         &
   )

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE conversions_mod, ONLY: pi
USE ukca_radaer_lut_read_in, ONLY: ukca_radaer_get_lut_index

! UKCA look-up tables
USE ukca_radaer_lut,        ONLY:                                              &
    ip_ukca_lut_accum,                                                         &
    ip_ukca_lut_coarse,                                                        &
    ip_ukca_lut_accnarrow,                                                     &
    ukca_lut

USE ukca_radaer_precalc,    ONLY:                                              &
    precalc

USE ukca_mode_setup,        ONLY:                                              &
    ip_ukca_mode_aitken,                                                       &
    ip_ukca_mode_accum,                                                        &
    ip_ukca_mode_coarse

USE ukca_radaer_struct_mod, ONLY:                                              &
    threshold_mmr,                                                             &
    threshold_vol,                                                             &
    threshold_nbr

USE ukca_radaer_ri_calc_mod, ONLY:                                             &
    ukca_radaer_ri_calc

IMPLICIT NONE

!
! Arguments with intent(in)
!


!
! Current spectrum
!
INTEGER, INTENT(IN) :: isolir
!
! Fixed array dimensions
!
INTEGER, INTENT(IN) :: npd_profile,                                            &
                       npd_layer,                                              &
                       npd_aerosol_mode,                                       &
                       npd_band,                                               &
                       npd_exclude

!
! Actual array dimensions
!
INTEGER, INTENT(IN) :: n_profile,                                              &
                       n_layer,                                                &
                       n_band,                                                 &
                       n_ukca_mode,                                            &
                       n_ukca_cpnt

!
! Fixed array dimensions for prescribed SSA
!
INTEGER, INTENT(IN) :: npd_prof_ssa,                                           &
                       npd_layr_ssa,                                           &
                       npd_band_ssa

!
! Variables related to waveband exclusion
!
LOGICAL, INTENT(IN) :: l_exclude
INTEGER, INTENT(IN) :: n_band_exclude(npd_band)
INTEGER, INTENT(IN) :: index_exclude(npd_exclude, npd_band)

!
! From ukca_radaer Structure for UKCA/radiation interaction
!
INTEGER, INTENT(IN) :: nmodes
INTEGER, INTENT(IN) :: ncp_max
INTEGER, INTENT(IN) :: ncp_max_x_nmodes
INTEGER, INTENT(IN) :: i_cpnt_index( ncp_max, nmodes )
INTEGER, INTENT(IN) :: i_cpnt_type( ncp_max_x_nmodes )
INTEGER, INTENT(IN) :: i_mode_type( nmodes )
LOGICAL, INTENT(IN) :: l_nitrate
LOGICAL, INTENT(IN) :: l_soluble( nmodes )
LOGICAL, INTENT(IN) :: l_sustrat
INTEGER, INTENT(IN) :: n_cpnt_in_mode( nmodes )

!
! Modal mass-mixing ratios
!
REAL, INTENT(IN) :: ukca_modal_mmr (npd_profile, npd_layer, npd_aerosol_mode)

!
! Modal number concentrations (m-3)
!
REAL, INTENT(IN) :: ukca_modal_number (npd_profile, npd_layer, n_ukca_mode)

!
! Dry and wet modal diameters
!
REAL, INTENT(IN) :: ukca_dry_diam (npd_profile, npd_layer, n_ukca_mode)
REAL, INTENT(IN) :: ukca_wet_diam (npd_profile, npd_layer, n_ukca_mode)

!
! Component volumes
!
REAL, INTENT(IN) :: ukca_cpnt_volume (n_ukca_cpnt, npd_profile, npd_layer)

!
! Modal volumes and densities
!
REAL, INTENT(IN) :: ukca_modal_volume  (npd_profile, npd_layer, n_ukca_mode)
REAL, INTENT(IN) :: ukca_modal_density (npd_profile, npd_layer, n_ukca_mode)

!
! Volume of water in modes
!
REAL, INTENT(IN) :: ukca_water_volume (npd_profile, npd_layer, n_ukca_mode)

!
! When true, arrays have been inverted
!
LOGICAL, INTENT(IN) :: l_inverted

!
! When true, use a prescribed single scattering albedo field
!
LOGICAL, INTENT(IN) :: l_ukca_radaer_prescribe_ssa

!
! Model level of tropopause
! Note levels are inverted in LFRic so we have to do something different here
!
INTEGER, INTENT(IN) :: trindxrad (npd_profile)

INTEGER, INTENT(IN) :: i_glomap_clim_tune_bc
INTEGER, INTENT(IN) :: i_ukca_tune_bc

!
! Prescription of single-scattering albedo
!
REAL, INTENT(IN) :: ukca_radaer_presc_ssa(npd_prof_ssa, npd_layr_ssa,          &
                                          npd_band_ssa)

!
! Arguments with intent(out)
!
! Band-averaged modal optical properties
!
REAL, INTENT(IN OUT) :: ukca_absorption (npd_profile, npd_layer,               &
                                        npd_aerosol_mode, npd_band)
REAL, INTENT(IN OUT) :: ukca_scattering (npd_profile, npd_layer,               &
                                        npd_aerosol_mode, npd_band)
REAL, INTENT(IN OUT) :: ukca_asymmetry  (npd_profile, npd_layer,               &
                                        npd_aerosol_mode, npd_band)

!
!
! Local variables
!
!
!
! Spectrum definitions
!

!
! Values at the point of integration:
!      Mie parameter for the wet and dry diameters and the indices of
!      their nearest neighbour
!      Complex refractive index and the index of its nearest neighbour
!
REAL :: x
INTEGER :: n_x
REAL :: x_dry
INTEGER :: n_x_dry
INTEGER :: n_nr

REAL    :: re_m(precalc%n_integ_pts)
REAL    :: im_m(precalc%n_integ_pts)
INTEGER :: n_ni(precalc%n_integ_pts)

!
! Integrals
!
REAL :: integrated_abs(npd_profile, npd_layer, n_ukca_mode, npd_band)
REAL :: integrated_sca(npd_profile, npd_layer, n_ukca_mode, npd_band)
REAL :: integrated_asy(npd_profile, npd_layer, n_ukca_mode, npd_band)
REAL :: loc_abs(precalc%n_integ_pts)
REAL :: loc_sca(precalc%n_integ_pts)
REAL :: loc_asy(precalc%n_integ_pts)
REAL :: loc_vol
REAL :: factor

!
! Waveband-integrated flux corrected for exclusions
!
REAL :: exclflux

!
! Local copy of single-scattering albedo to prescribe.
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
REAL :: ni_c_power
INTEGER, PARAMETER :: n_ni_fix = 1

!
! Local copies of mode type, component index and component type
!
INTEGER :: this_mode_type

!
! Loop variables
!
INTEGER :: i_band, & ! loop on wavebands
        i_mode, & ! loop on aerosol modes
        i_layr, & ! loop on vertical dimension
        i_prof, & ! loop on horizontal dimension
        i_intg    ! loop on integration points and excluded bands

!
! Thresholds on the modal mass-mixing ratio, volume, and modal number
! concentrations above which aerosol optical properties are to be
! computed - threshold_mmr, threshold_vol, threshold_nbr - specified
! in ukca_radaer_struct_mod
!
! Limits for the asymmetry parameter, since values of
! exactly -1.0 or +1.0 can cause div-by-zero errors
! further on in the Radiation code.
REAL, PARAMETER :: minus1_plus_epsi1 = -1.0 + EPSILON(1.0)
REAL, PARAMETER :: one_minus_epsi1 = 1.0 - EPSILON(1.0)

!
! Indicates whether current level is above the tropopause.
!
LOGICAL :: l_in_stratosphere

! error message
CHARACTER (LEN=errormessagelength) :: cmessage
! error indicator
INTEGER:: icode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_BAND_AVERAGE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!
! To band-average modal optical properties, we need first to compute
! adequate indices in the look-up tables. For that, we need:
! *** the modal dry radius (we've got the diameter as input)
! *** the modal wet radius (we've got the diameter as input)
! *** the modal refractive index
! In addition, in order to output specific coefficients for absorption
! and scattering (in m2/kg from m-1), we need the modal density.
!

! If the single scattering albedo is prescribed then n_ni == 1
! and initialise loc_abs to zero
IF (l_ukca_radaer_prescribe_ssa) THEN
  DO i_intg = 1, precalc%n_integ_pts
    n_ni(i_intg) = n_ni_fix
    loc_abs(i_intg) = 0.0
  END DO
END IF

DO i_band = 1, n_band

  DO i_mode = 1, n_ukca_mode

    !
    ! Mode type. From a look-up table point of view, Aitken and
    ! accumulation types are treated in the same way.
    ! Accumulation soluble mode may use a narrower width (i.e. another
    ! look-up table) than other Aitken and accumulation modes.
    ! Once we know which look-up table to select, make local copies
    ! of info needed for nearest-neighbour calculations.
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
      this_mode_type = ip_ukca_lut_coarse

    CASE DEFAULT
      ! Likely developer is trying to pass nucleation mode to radaer
      icode = 1
      cmessage = 'Mode is not one of aitken , accumulation , coarse'
      CALL ereport(RoutineName,icode,cmessage)
    END SELECT

    nx      = ukca_lut(this_mode_type, isolir)%n_x
    logxmin = LOG(ukca_lut(this_mode_type, isolir)%x_min)
    logxmaxmlogxmin =                                                          &
              LOG(ukca_lut(this_mode_type, isolir)%x_max) - logxmin

    nnr     = ukca_lut(this_mode_type, isolir)%n_nr
    nrmin   = ukca_lut(this_mode_type, isolir)%nr_min
    incr_nr = ukca_lut(this_mode_type, isolir)%incr_nr

    nni     = ukca_lut(this_mode_type, isolir)%n_ni
    ni_min  = ukca_lut(this_mode_type, isolir)%ni_min
    ni_max  = ukca_lut(this_mode_type, isolir)%ni_max
    ni_c    = ukca_lut(this_mode_type, isolir)%ni_c
    ni_c_power = 10.0**ni_c

    !
    ! Wavelength-dependent calculations.
    ! Waveband-integration is done at the same time, so most computed
    ! items are not stored in arrays.
    !
    DO i_layr = 1, n_layer

      DO i_prof = 1, n_profile

        IF (l_inverted) THEN
          l_in_stratosphere = i_layr <= trindxrad(i_prof)
        ELSE
          l_in_stratosphere = i_layr >= trindxrad(i_prof)
        END IF

        !
        ! Only make calculations if there are some aerosols, and
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

        IF (ukca_modal_mmr   (i_prof, i_layr, i_mode) > threshold_mmr .AND.    &
            ukca_modal_number(i_prof, i_layr, i_mode) > threshold_nbr .AND.    &
            ukca_modal_volume(i_prof, i_layr, i_mode) > threshold_vol) THEN

          DO i_intg = 1, precalc%n_integ_pts

             ! Compute the modal complex refractive index via
             ! volume-weighting for non-BC components. If
             ! i_ukca_tune_bc or i_glomap_clim_tune_bc are set to
             ! i_ukca_bc_mg_mix then the BC component will be added
             ! via the Maxwell-Garnet mixing approach, otherwise use
             ! volume-weighting for BC also.
             !
            CALL ukca_radaer_ri_calc(                                          &
                 ! From the structure ukca_radaer for UKCA/radiation interaction
                 nmodes, ncp_max, ncp_max_x_nmodes,                            &
                 i_cpnt_index, i_cpnt_type, n_cpnt_in_mode,                    &
                 l_nitrate, l_soluble, l_sustrat,                              &
                 ! Refractive index
                 precalc%realrefr( :, i_intg, i_band, isolir ),                &
                 precalc%imagrefr( :, i_intg, i_band, isolir ),                &
                 ! Modal properties
                 ukca_cpnt_volume( :, i_prof, i_layr ),                        &
                 ukca_modal_volume( i_prof, i_layr, i_mode ),                  &
                 ukca_water_volume( i_prof, i_layr, i_mode ),                  &
                 ! Indicies for arrays
                 i_mode, n_ukca_cpnt,                                          &
                 ! Stratospheric aerosol treated as sulphuric acid?
                 l_in_stratosphere,                                            &
                 ! Logical control switches
                 i_ukca_tune_bc, i_glomap_clim_tune_bc,                        &
                 l_ukca_radaer_prescribe_ssa,                                  &
                 ! Output refractive index real and imag parts
                 re_m(i_intg), im_m(i_intg) )

          END DO  ! i_intg

          ! Do not calculate the index of the imaginary component if
          ! SSA is prescribed
          IF (.NOT. l_ukca_radaer_prescribe_ssa) THEN
            CALL ukca_radaer_get_lut_index(                                    &
                 nni, im_m, ni_min, ni_max, ni_c, n_ni,                        &
                 precalc%n_integ_pts, ni_c_power=ni_c_power)
          END IF

          DO i_intg = 1, precalc%n_integ_pts

            !
            ! Compute the Mie parameter from the wet diameter
            ! and get the LUT-array index of its nearest neighbour.
            !
            x = pi * ukca_wet_diam(i_prof, i_layr, i_mode) /                   &
                precalc%wavelength(i_intg, i_band, isolir)
            n_x = NINT( (LOG(x)    - logxmin) /                                &
                         logxmaxmlogxmin * (nx-1) ) + 1
            n_x = MIN(nx, MAX(1, n_x))

            !
            ! Same for the dry diameter (needed to access the volume
            ! fraction)
            !
            x_dry = pi * ukca_dry_diam(i_prof, i_layr, i_mode) /               &
                    precalc%wavelength(i_intg, i_band, isolir)
            n_x_dry = NINT( (LOG(x_dry) - logxmin) /                           &
                             logxmaxmlogxmin * (nx-1) ) + 1
            n_x_dry = MIN(nx, MAX(1, n_x_dry))

            !
            ! Compute the modal complex refractive index as
            ! volume-weighted component refractive indices.
            ! Get the LUT-array index of their nearest neighbours.
            !

            n_nr = NINT( (re_m(i_intg) - nrmin) / incr_nr ) + 1
            n_nr = MIN(nnr, MAX(1, n_nr))

            IF (l_ukca_radaer_prescribe_ssa) THEN

              !
              ! Get local copies of the relevant look-up table entries.
              !
              loc_sca(i_intg) = ukca_lut(this_mode_type, isolir)%              &
                           ukca_scattering(n_x, n_ni(i_intg), n_nr)

              loc_asy(i_intg) = ukca_lut(this_mode_type, isolir)%              &
                           ukca_asymmetry(n_x, n_ni(i_intg), n_nr)

              loc_vol = ukca_lut(this_mode_type, isolir)%                      &
                        volume_fraction(n_x_dry)

              !
              ! Offline Mie calculations were integrated using the Mie
              ! parameter. Compared to an integration using the particle
              ! radius, extra factors are introduced. Absorption and
              ! scattering efficiencies must be multiplied by the squared
              ! wavelength, and the volume fraction by the cubed wavelength.
              ! Consequently, ratios abs/volfrac and sca/volfrac have then
              ! to be divided by the wavelength.
              ! We also weight by the solar irradiance or Planckian
              ! irradiance.
              !
              factor = precalc%irrad(i_intg, i_band, isolir) /                 &
                 (ukca_modal_density(i_prof, i_layr, i_mode) * loc_vol *       &
                  precalc%wavelength(i_intg, i_band, isolir))
              !
              ! The single-scattering albedo is prescribed by distributing
              ! extinction (which is equal to scattering in the non-absorption
              ! case) to absorption and scattering coefficients in the
              ! proportion indicated by the prescription.
              !
              this_ssa = ukca_radaer_presc_ssa(i_prof, i_layr, i_band)
              loc_abs(i_intg) = loc_sca(i_intg) * factor * (1.0 - this_ssa)
              loc_sca(i_intg) = loc_sca(i_intg) * factor * this_ssa
              loc_asy(i_intg) = loc_asy(i_intg) * loc_sca(i_intg)

            ELSE

              !
              ! Get local copies of the relevant look-up table entries.
              !
              loc_abs(i_intg) = ukca_lut(this_mode_type, isolir)%              &
                           ukca_absorption(n_x, n_ni(i_intg), n_nr)

              loc_sca(i_intg) = ukca_lut(this_mode_type, isolir)%              &
                           ukca_scattering(n_x, n_ni(i_intg), n_nr)

              loc_asy(i_intg) = ukca_lut(this_mode_type, isolir)%              &
                           ukca_asymmetry(n_x, n_ni(i_intg), n_nr)

              loc_vol = ukca_lut(this_mode_type, isolir)%                      &
                        volume_fraction(n_x_dry)

              !
              ! Offline Mie calculations were integrated using the Mie
              ! parameter. Compared to an integration using the particle
              ! radius, extra factors are introduced. Absorption and
              ! scattering efficiencies must be multiplied by the squared
              ! wavelength, and the volume fraction by the cubed wavelength.
              ! Consequently, ratios abs/volfrac and sca/volfrac have then
              ! to be divided by the wavelength.
              ! We also weight by the solar irradiance or Planckian
              ! irradiance.
              !
              factor = precalc%irrad(i_intg, i_band, isolir) /                 &
                 (ukca_modal_density(i_prof, i_layr, i_mode) * loc_vol *       &
                  precalc%wavelength(i_intg, i_band, isolir))
              loc_abs(i_intg) = loc_abs(i_intg) * factor
              loc_sca(i_intg) = loc_sca(i_intg) * factor
              loc_asy(i_intg) = loc_asy(i_intg) * loc_sca(i_intg)

            END IF ! IF (l_ukca_radaer_prescribe_ssa)

          END DO ! i_intg

          !
          ! Trapezoidal integration
          !
          integrated_abs(i_prof, i_layr, i_mode, i_band) = 0.0e+00
          integrated_sca(i_prof, i_layr, i_mode, i_band) = 0.0e+00
          integrated_asy(i_prof, i_layr, i_mode, i_band) = 0.0e+00
          DO i_intg = 1, precalc%n_integ_pts - 1
            integrated_abs(i_prof, i_layr, i_mode, i_band) =                   &
               integrated_abs(i_prof, i_layr, i_mode, i_band) +                &
               (precalc%wavelength(i_intg+1, i_band, isolir) -                 &
                precalc%wavelength(i_intg, i_band, isolir)) *                  &
               (loc_abs(i_intg+1) + loc_abs(i_intg))
            integrated_sca(i_prof, i_layr, i_mode, i_band) =                   &
               integrated_sca(i_prof, i_layr, i_mode, i_band) +                &
               (precalc%wavelength(i_intg+1, i_band, isolir) -                 &
                precalc%wavelength(i_intg, i_band, isolir)) *                  &
               (loc_sca(i_intg+1) + loc_sca(i_intg))
            integrated_asy(i_prof, i_layr, i_mode, i_band) =                   &
               integrated_asy(i_prof, i_layr, i_mode, i_band) +                &
               (precalc%wavelength(i_intg+1, i_band, isolir) -                 &
                precalc%wavelength(i_intg, i_band, isolir)) *                  &
               (loc_asy(i_intg+1) + loc_asy(i_intg))
          END DO ! i_intg
          integrated_abs(i_prof, i_layr, i_mode, i_band) =                     &
             integrated_abs(i_prof, i_layr, i_mode, i_band) * 0.5
          integrated_sca(i_prof, i_layr, i_mode, i_band) =                     &
             integrated_sca(i_prof, i_layr, i_mode, i_band) * 0.5
          integrated_asy(i_prof, i_layr, i_mode, i_band) =                     &
             integrated_asy(i_prof, i_layr, i_mode, i_band) * 0.5

        ELSE

          integrated_abs(i_prof, i_layr, i_mode, i_band) = 0.0e+00
          integrated_sca(i_prof, i_layr, i_mode, i_band) = 0.0e+00
          integrated_asy(i_prof, i_layr, i_mode, i_band) = 0.0e+00

        END IF

      END DO ! i_prof

    END DO ! i_layr

  END DO ! i_mode

END DO ! i_band

!
! Final integrals. Depend on excluded bands.
!

DO i_band = 1, n_band

  IF (l_exclude) THEN

    IF (n_band_exclude(i_band) > 0) THEN

      !
      ! Remove contribution from excluded bands.
      !
      DO i_intg = 1, n_band_exclude(i_band)

        DO i_mode = 1, n_ukca_mode

          DO i_layr = 1, n_layer

            DO i_prof = 1, n_profile

              integrated_abs(i_prof, i_layr, i_mode, i_band) =                 &
                        integrated_abs(i_prof, i_layr, i_mode, i_band) -       &
                        integrated_abs(i_prof, i_layr, i_mode,                 &
                                       index_exclude(i_intg, i_band))
              integrated_sca(i_prof, i_layr, i_mode, i_band) =                 &
                        integrated_sca(i_prof, i_layr, i_mode, i_band) -       &
                        integrated_sca(i_prof, i_layr, i_mode,                 &
                                       index_exclude(i_intg, i_band))
              integrated_asy(i_prof, i_layr, i_mode, i_band) =                 &
                        integrated_asy(i_prof, i_layr, i_mode, i_band) -       &
                        integrated_asy(i_prof, i_layr, i_mode,                 &
                                       index_exclude(i_intg, i_band))

            END DO ! i_prof

          END DO ! i_layr

        END DO ! i_mode

        exclflux = precalc%flux(i_band, isolir) -                              &
                   precalc%flux(index_exclude(i_intg, i_band), isolir)

      END DO ! i_intg

    ELSE

      exclflux = precalc%flux(i_band, isolir)

    END IF

  ELSE

    exclflux = precalc%flux(i_band, isolir)

  END IF

  DO i_mode = 1, n_ukca_mode

    DO i_layr = 1, n_layer

      DO i_prof = 1, n_profile

        !
        ! Pathological combinations of Mie parameters and refractive index
        ! may cause unphysical values, especially for accumulation-mode
        ! aerosols in the longwave spectrum. Also, band exclusion can yield
        ! negative (albeit small) scattering or absorption coefficients.
        !
        ! Here, we make sure that optical properties remain within sensible
        ! bounds: specific scattering and absorption coefficients must be
        ! positive, and asymmetry parameter must be within [-1,+1].
        !

        ! First check absorption and scattering
        !

        IF (integrated_abs(i_prof, i_layr, i_mode, i_band) < 0.0e+00) THEN

          integrated_abs(i_prof, i_layr, i_mode, i_band) = 0.0e+00
        END IF

        IF (integrated_sca(i_prof, i_layr, i_mode, i_band) < 0.0e+00) THEN

          integrated_sca(i_prof, i_layr, i_mode, i_band) = 0.0e+00
        END IF

        ! Calculate ukca_absorption, ukca_scatterig using exclflux
        ! If exclflux  <= 0 then skip the calculation and set
        ! ukca_absorption, ukca_scattering to zero.

        IF (exclflux > 0.0e+00) THEN
          !
          ukca_absorption(i_prof, i_layr, i_mode, i_band) =                    &
            integrated_abs(i_prof, i_layr, i_mode, i_band) / exclflux

          ukca_scattering(i_prof, i_layr, i_mode, i_band) =                    &
            integrated_sca(i_prof, i_layr, i_mode, i_band) / exclflux
          !
        ELSE
          !
          ukca_absorption(i_prof, i_layr, i_mode, i_band) = 0.0e+00
          ukca_scattering(i_prof, i_layr, i_mode, i_band) = 0.0e+00

        END IF

        ! Calculate asymmetry parameter

        IF (integrated_sca(i_prof, i_layr, i_mode, i_band) >  0.0e+00) THEN

          ukca_asymmetry(i_prof, i_layr, i_mode, i_band)  =                    &
            integrated_asy(i_prof, i_layr, i_mode, i_band) /                   &
            integrated_sca(i_prof, i_layr, i_mode, i_band)

        ELSE

          ukca_asymmetry(i_prof, i_layr, i_mode, i_band) = 0.0e+00

        END IF

        ! Check that asymmetry parameter has physical values [-1, 1]
        ! but do not allow exactly 1 or -1 as this can cause
        ! divide by zero elsewhere in the radiation code. Uses a
        ! deviation of EPSILON(1.0) from +/- 1.0

        IF (ukca_asymmetry(i_prof, i_layr, i_mode, i_band) <                   &
            minus1_plus_epsi1) THEN

          ukca_asymmetry(i_prof, i_layr, i_mode, i_band) = minus1_plus_epsi1

        ELSE IF (ukca_asymmetry(i_prof, i_layr, i_mode, i_band) >              &
            one_minus_epsi1) THEN

          ukca_asymmetry(i_prof, i_layr, i_mode, i_band) = one_minus_epsi1

        END IF

      END DO ! i_prof

    END DO ! i_layr

  END DO  ! i_mode

END DO ! i_band


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_band_average

END MODULE ukca_radaer_band_average_mod
