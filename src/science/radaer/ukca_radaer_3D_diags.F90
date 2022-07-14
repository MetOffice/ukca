! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Calculates the 3D field of aerosol extinction coefficient, absorption,
!  scattering, and gsca for UKCA-MODE aerosols.
!
MODULE ukca_radaer_3d_diags_mod

IMPLICIT NONE

! Description:
!   Outputs extinction, absorption, scattering and asymmetry * scattering
!   as 3D fields in units (m-1) forUKCA 3D aerosol-radiation diagnostics
!
! Method:
!   Loops over all radiatively active modes.
!   Specific scattering and absorption are multiplied by mass mixing ratio and
!   air density to calculate scattering, absorption and extinction
!   in units (m-1). Asymmetry is multiplied by scattering to give gsca in
!   units (m-1).
!
! Output:  ukca_aerosol_ext, ukca_aerosol_abs
!          ukca_aerosol_abs, ukca_aerosol_gsca
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
! Code description:
!   Language: Fortran 95.
!   This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_RADAER_3D_DIAGS_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE ukca_radaer_3d_diags(                                               &
      ! Fixed array dimensions
      nd_profile,  nd_layer, nd_aerosol_mode,                                  &
      ! Actual array dimensions
      n_profile, n_layer, n_ukca_mode, n_ukca_cpnt,                            &
      ! UKCA_RADAER structure
      ukca_radaer,                                                             &
      ! Modal diameters from UKCA module
      ukca_dry, ukca_wet,                                                      &
      ! Air density
      air_density,                                                             &
      ! Component volumes
      ukca_cpnt_volume,                                                        &
      ! Modal volumes, densities, and water content
      ukca_modal_volume, ukca_modal_density, ukca_water_volume,                &
      ! Modal mass-mixing ratios
      ukca_modal_mmr,                                                          &
      ! Modal number concentrations
      ukca_modal_number,                                                       &
      ! Model level of the tropopause
      trindxrad,                                                               &
      ! Index of wavelength to consider
      i_wavel,                                                                 &
      ! 3D aerosol extinction and absorption
      ukca_aerosol_ext, ukca_aerosol_abs,                                      &
      ! 3D aerosol scattering and scattering*asymmetry
      ukca_aerosol_sca, ukca_aerosol_gsca)

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE conversions_mod, ONLY: pi

USE ukca_radaer_lut,         ONLY:                                             &
    ip_ukca_lut_accum,                                                         &
    ip_ukca_lut_coarse,                                                        &
    ip_ukca_lut_accnarrow,                                                     &
    ip_ukca_lut_sw,                                                            &
    ukca_lut

USE ukca_radaer_precalc,     ONLY:                                             &
    precalc

USE ukca_mode_setup,         ONLY:                                             &
    ip_ukca_mode_aitken,                                                       &
    ip_ukca_mode_accum,                                                        &
    ip_ukca_mode_coarse

USE ukca_radaer_struct_mod,  ONLY:                                             &
    ukca_radaer_struct,                                                        &
    threshold_mmr, threshold_vol,                                              &
    threshold_nbr

USE ukca_radaer_lut_read_in, ONLY:                                             &
    ukca_radaer_get_lut_index

USE ukca_radaer_ri_calc_mod, ONLY:                                             &
    ukca_radaer_ri_calc

USE ukca_option_mod,         ONLY:                                             &
    l_ukca_tune_bc

USE glomap_clim_option_mod,  ONLY:                                             &
    l_glomap_clim_tune_bc

IMPLICIT NONE

! Subroutine arguments

! Arguments with intent(in)

! Fixed array dimensions
INTEGER, INTENT(IN) ::                                                         &
     nd_profile,                                                               &
     nd_layer,                                                                 &
     nd_aerosol_mode
!
! Actual array dimensions
INTEGER, INTENT(IN) ::                                                         &
     n_profile,                                                                &
     n_layer,                                                                  &
     n_ukca_mode,                                                              &
     n_ukca_cpnt
!
! Structure for UKCA/radiation interaction
TYPE (ukca_radaer_struct), INTENT(IN) :: ukca_radaer

!
! UKCA-MODE aerosol properties
!
REAL, INTENT(IN) ::                                                            &
! Dry and wet modal diameters
     ukca_dry (nd_profile, nd_layer, n_ukca_mode),                             &
     ukca_wet (nd_profile, nd_layer, n_ukca_mode),                             &
! Air density
     air_density(nd_profile, nd_layer),                                        &
! Component volumes
     ukca_cpnt_volume (n_ukca_cpnt, nd_profile, nd_layer),                     &
! Modal volumes
     ukca_modal_volume (nd_profile, nd_layer, n_ukca_mode),                    &
! Modal densities
     ukca_modal_density (nd_profile, nd_layer, n_ukca_mode),                   &
! Volume of water in modes
     ukca_water_volume (nd_profile, nd_layer, n_ukca_mode),                    &
! Modal mass-mixing ratios
     ukca_modal_mmr (nd_profile, nd_layer, nd_aerosol_mode),                   &
! Modal number concentrations (m-3)
     ukca_modal_number (nd_profile, nd_layer, n_ukca_mode)

! Model level of tropopause
INTEGER, INTENT(IN) :: trindxrad (nd_profile)

! Index of wavelength to consider
INTEGER, INTENT(IN) :: i_wavel

!
! Arguments with intent out
!

! Aerosol extinction scattering, absorption and scattering*asymmetry profiles
REAL, INTENT(OUT) ::                                                           &
     ukca_aerosol_ext (nd_profile, nd_layer),                                  &
     ukca_aerosol_abs (nd_profile, nd_layer),                                  &
     ukca_aerosol_sca (nd_profile, nd_layer),                                  &
     ukca_aerosol_gsca (nd_profile, nd_layer)

!
! Loop variables
!
INTEGER :: i_mode, & ! loop on aerosol modes
           i_layr, & ! loop on vertical dimension
           i_prof    ! loop on horizontal dimension

! Mie parameter for the wet and dry diameters
REAL    :: x, x_dry

! Index for accessing x, x_dry in look-up-table
INTEGER :: n_x, n_x_dry

!      Complex refractive index and its nearest neighbour
REAL    :: re_m, im_m

! Indices for accessing refractive index in look-up-table
INTEGER :: n_nr, n_ni

! Values at given wavelength
REAL :: loc_abs
REAL :: loc_sca
REAL :: loc_gsca

! Volume fraction from look-up-table
REAL :: loc_vol

! Local factor used in optical calculation
REAL :: factor

! Local copies of typedef members
INTEGER :: nx
REAL    :: logxmin         ! log(xmin)
REAL    :: logxmaxmlogxmin ! log(xmax) - log(xmin)
INTEGER :: nnr
REAL    :: nrmin
REAL    :: incr_nr
INTEGER :: nni
REAL    :: ni_min
REAL    :: ni_max
REAL    :: ni_c

! Local copies of mode type, component index and component type
INTEGER :: this_mode_type

! Thresholds on the modal mass-mixing ratio and modal number
! concentrations above which aerosol optical properties are to be
! computed - threshold_mmr, threshold_vol, threshold_nbr - specified
! in ukca_radaer_struct_mod

! Indicates whether the model level is above the tropopause.
!
LOGICAL :: l_in_stratosphere

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_3D_DIAGS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialize output arrays
ukca_aerosol_ext(:,:) = 0.0
ukca_aerosol_abs(:,:) = 0.0
ukca_aerosol_sca(:,:) = 0.0
ukca_aerosol_gsca(:,:) = 0.0

! Begin large loop over UKCA modes
DO i_mode = 1, n_ukca_mode
  ! Mode type. From a look-up table point of view, Aitken and
  ! accumulation types are treated in the same way.
  ! Accumulation soluble mode may use a narrower width (i.e. another
  ! look-up table) than other Aitken and accumulation modes.
  ! Once we know which look-up table to select, make local copies
  ! of info needed for nearest-neighbour calculations.
  SELECT CASE (ukca_radaer%i_mode_type(i_mode))
  CASE (ip_ukca_mode_aitken)
    this_mode_type = ip_ukca_lut_accum
  CASE (ip_ukca_mode_accum)
    IF (ukca_radaer%l_soluble(i_mode)) THEN
      this_mode_type = ip_ukca_lut_accnarrow
    ELSE
      this_mode_type = ip_ukca_lut_accum
    END IF
  CASE (ip_ukca_mode_coarse)
    this_mode_type = ip_ukca_lut_coarse
  END SELECT

  nx      = ukca_lut(this_mode_type, ip_ukca_lut_sw)%n_x
  logxmin = LOG(ukca_lut(this_mode_type, ip_ukca_lut_sw)%x_min)
  logxmaxmlogxmin = LOG(ukca_lut(this_mode_type, ip_ukca_lut_sw)%x_max) -      &
                    logxmin

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
      IF ((ukca_modal_mmr   (i_prof, i_layr, i_mode) > threshold_mmr) .AND.    &
          (ukca_modal_number(i_prof, i_layr, i_mode) > threshold_nbr) .AND.    &
          (ukca_modal_volume(i_prof, i_layr, i_mode) > threshold_vol)) THEN
        ! Compute the Mie parameter from the wet diameter
        ! and get the LUT-array index of its nearest neighbour.
        x = pi * ukca_wet(i_prof, i_layr, i_mode) /                            &
                 precalc%aod_wavel(i_wavel)
        n_x = NINT( (LOG(x)    - logxmin) / logxmaxmlogxmin * (nx-1) ) + 1
        n_x = MIN(nx, MAX(1, n_x))

        ! Same for the dry diameter (needed to access the volume
        ! fraction)
        x_dry = pi * ukca_dry(i_prof, i_layr, i_mode)                          &
             / precalc%aod_wavel(i_wavel)
        n_x_dry = NINT( (LOG(x_dry) - logxmin) /                               &
             logxmaxmlogxmin * (nx-1) ) + 1
        n_x_dry = MIN(nx, MAX(1, n_x_dry))

        ! Compute the modal complex refractive index via
        ! volume-weighting for non-BC components. If
        ! l_ukca_tune_bc is set add the BC component
        ! via the Maxwell-Garnett mixing approach, otherwise use
        ! volume-weighting for BC also.
        !
        CALL ukca_radaer_ri_calc(                                              &
             ! From the structure ukca_radaer for UKCA/radiation interaction
             ukca_radaer%nmodes,                                               &
             ukca_radaer%ncp_max,                                              &
             ukca_radaer%ncp_max_x_nmodes,                                     &
             ukca_radaer%i_cpnt_index,                                         &
             ukca_radaer%i_cpnt_type,                                          &
             ukca_radaer%n_cpnt_in_mode,                                       &
             ukca_radaer%l_nitrate,                                            &
             ukca_radaer%l_soluble,                                            &
             ukca_radaer%l_sustrat,                                            &
             ! Refractive index
             precalc%aod_realrefr( :, i_wavel ),                               &
             precalc%aod_imagrefr( :, i_wavel ),                               &
             ! Modal properties
             ukca_cpnt_volume( :, i_prof, i_layr ),                            &
             ukca_modal_volume( i_prof, i_layr, i_mode ),                      &
             ukca_water_volume( i_prof, i_layr, i_mode ),                      &
             ! Indicies for arrays
             i_mode, n_ukca_cpnt,                                              &
             ! Stratospheric aerosol treated as sulphuric acid?
             l_in_stratosphere,                                                &
             ! Logical control switches
             l_ukca_tune_bc, l_glomap_clim_tune_bc,                            &
             ! Output refractive index real and imag parts
             re_m, im_m )

        n_nr = NINT( (re_m - nrmin) / incr_nr ) + 1
        n_nr = MIN(nnr, MAX(1, n_nr))

        CALL ukca_radaer_get_lut_index( nni, im_m, ni_min, ni_max, ni_c, n_ni)
        !
        ! Get local copies of the relevant look-up table entries.
        !
        loc_abs = ukca_lut(this_mode_type, ip_ukca_lut_sw)%                    &
             ukca_absorption(n_x, n_ni, n_nr)

        loc_sca = ukca_lut(this_mode_type, ip_ukca_lut_sw)%                    &
             ukca_scattering(n_x, n_ni, n_nr)

        loc_gsca = ukca_lut(this_mode_type, ip_ukca_lut_sw)%                   &
             ukca_asymmetry (n_x, n_ni, n_nr)                                  &
             * ukca_lut(this_mode_type, ip_ukca_lut_sw)%                       &
             ukca_scattering(n_x, n_ni, n_nr)

        loc_vol = ukca_lut(this_mode_type, ip_ukca_lut_sw)%                    &
             volume_fraction(n_x_dry)

        ! Offline Mie calculations were integrated using the Mie
        ! parameter. Compared to an integration using the particle
        ! radius, extra factors are introduced. Absorption and
        ! scattering efficiencies must be multiplied by the squared
        ! wavelength, and the volume fraction by the cubed wavelength.
        ! Consequently, ratios abs/volfrac and sca/volfrac have then
        ! to be divided by the wavelength.
        !
        factor = (ukca_modal_density(i_prof, i_layr, i_mode) * loc_vol *       &
                  precalc%aod_wavel(i_wavel))

        loc_abs = loc_abs / factor
        loc_sca = loc_sca / factor
        loc_gsca = loc_gsca / factor
        !
        ! aerosol extinction profile
        !
        ukca_aerosol_ext(i_prof, i_layr) =                                     &
             ukca_aerosol_ext(i_prof, i_layr) +                                &
             ukca_modal_mmr(i_prof, i_layr, i_mode)                            &
             * air_density(i_prof, i_layr)                                     &
             * (loc_abs + loc_sca)

        !
        ! aerosol absorption profile
        !
        ukca_aerosol_abs(i_prof, i_layr) =                                     &
             ukca_aerosol_abs(i_prof, i_layr) +                                &
             ukca_modal_mmr(i_prof, i_layr, i_mode)                            &
             * air_density(i_prof, i_layr)                                     &
             * loc_abs
        !
        ! aerosol scattering profile
        !
        ukca_aerosol_sca(i_prof, i_layr) =                                     &
             ukca_aerosol_ext(i_prof, i_layr) -                                &
             ukca_aerosol_abs(i_prof, i_layr)

        !
        ! aerosol scattering * asymmetry profile
        !
        ukca_aerosol_gsca(i_prof, i_layr) =                                    &
             ukca_aerosol_gsca(i_prof, i_layr) +                               &
             ukca_modal_mmr(i_prof, i_layr, i_mode)                            &
             * air_density(i_prof, i_layr)                                     &
             * loc_gsca
      END IF ! mmr, number and volume thresholds
    END DO ! i_prof
  END DO ! i_layr
END DO ! i_mode

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE ukca_radaer_3d_diags

END MODULE ukca_radaer_3d_diags_mod
