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
MODULE ukca_radaer_band_average_faster_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName = 'UKCA_RADAER_BAND_AVERAGE_FASTER_MOD'

CONTAINS

SUBROUTINE ukca_radaer_band_average_faster()

USE conversions_mod,        ONLY: pi
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,               ONLY: jpim, jprb
USE vectlib_mod,            ONLY: log_v
USE yomhook,                ONLY: lhook, dr_hook

IMPLICIT NONE

! Arguments

! Current spectrum
INTEGER, INTENT(IN) :: isolir

! Fixed array dimensions
INTEGER, INTENT(IN) :: npd_profile,                                            &
                       npd_layer,                                              &
                       npd_aerosol_mode,                                       &
                       npd_band,                                               &
                       npd_exclude

! Actual array dimensions
INTEGER, INTENT(IN) :: n_profile,                                              &
                       n_layer,                                                &
                       n_band,                                                 &
                       n_ukca_mode,                                            &
                       n_ukca_cpnt

! Fixed array dimensions for prescribed SSA
INTEGER, INTENT(IN) :: npd_prof_ssa,                                           &
                       npd_layr_ssa,                                           &
                       npd_band_ssa

! Variables related to waveband exclusion
LOGICAL, INTENT(IN) :: l_exclude
INTEGER, INTENT(IN) :: n_band_exclude(npd_band)
INTEGER, INTENT(IN) :: index_exclude(npd_exclude, npd_band)

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

! Modal mass-mixing ratios
REAL, INTENT(IN) :: ukca_modal_mmr (npd_profile, npd_layer, npd_aerosol_mode)

! Modal number concentrations (m-3)
REAL, INTENT(IN) :: ukca_modal_number (npd_profile, npd_layer, n_ukca_mode)

! Dry and wet modal diameters
REAL, INTENT(IN) :: ukca_dry_diam (npd_profile, npd_layer, n_ukca_mode)
REAL, INTENT(IN) :: ukca_wet_diam (npd_profile, npd_layer, n_ukca_mode)

! Component volumes
REAL, INTENT(IN) :: ukca_cpnt_volume (n_ukca_cpnt, npd_profile, npd_layer)

! Modal volumes and densities
REAL, INTENT(IN) :: ukca_modal_volume  (npd_profile, npd_layer, n_ukca_mode)
REAL, INTENT(IN) :: ukca_modal_density (npd_profile, npd_layer, n_ukca_mode)

! Volume of water in modes
REAL, INTENT(IN) :: ukca_water_volume (npd_profile, npd_layer, n_ukca_mode)

! When true, arrays have been inverted
LOGICAL, INTENT(IN) :: l_inverted

! When > 0, use a prescribed single scattering albedo field
INTEGER, INTENT(IN) :: i_ukca_radaer_prescribe_ssa

! Model level of tropopause
! Note levels are inverted in LFRic so we have to do something different here
INTEGER, INTENT(IN) :: trindxrad (npd_profile)

! Get rid of these arguments
!INTEGER, INTENT(IN) :: i_glomap_clim_tune_bc
!INTEGER, INTENT(IN) :: i_ukca_tune_bc

! Prescription of single-scattering albedo
REAL, INTENT(IN) :: ukca_radaer_presc_ssa( npd_prof_ssa, npd_layr_ssa,         &
                                           npd_band_ssa)

! Band-averaged modal optical properties
REAL, INTENT(IN OUT) :: ukca_absorption ( npd_profile, npd_layer,              &
                                          npd_aerosol_mode, npd_band)

REAL, INTENT(IN OUT) :: ukca_scattering ( npd_profile, npd_layer,              &
                                          npd_aerosol_mode, npd_band)

REAL, INTENT(IN OUT) :: ukca_asymmetry  ( npd_profile, npd_layer,              &
                                          npd_aerosol_mode, npd_band)

!
! Local variables
!

INTEGER, PARAMETER :: one = 1

! Values at the point of integration:
!      Mie parameter for the wet and dry diameters and the indices of
!      their nearest neighbour
!      Complex refractive index and the index of its nearest neighbour
REAL :: x
INTEGER :: n_x
REAL :: x_dry
INTEGER :: n_x_dry
INTEGER :: n_nr

! Real part of refractive index
REAL    :: re_m( npd_profile, npd_layer, npd_band, npd_aerosol_mode )
! Imaginary part of refractive index
REAL    :: im_m( npd_profile, npd_layer, npd_band, npd_aerosol_mode )
! Index
INTEGER :: ni_ind( npd_profile, npd_layer, npd_band, npd_aerosol_mode )

! Integrals
REAL :: loc_abs
REAL :: loc_sca( npd_profile, npd_layer, npd_band, npd_aerosol_mode )
REAL :: loc_asy( npd_profile, npd_layer, npd_band, npd_aerosol_mode )
REAL :: loc_vol
REAL :: factor( npd_profile, npd_layer, npd_band, npd_aerosol_mode )

! Local copy of single-scattering albedo to prescribe.
REAL :: this_ssa

! Local copies of typedef members
INTEGER :: nx( npd_aerosol_mode )
REAL :: logxmin( npd_aerosol_mode )         ! log(xmin)
REAL :: logxmaxmlogxmin( npd_aerosol_mode ) ! log(xmax) - log(xmin)
INTEGER :: nnr( npd_aerosol_mode )
REAL :: nrmin( npd_aerosol_mode )
REAL :: incr_nr( npd_aerosol_mode )
INTEGER :: nni( npd_aerosol_mode )
REAL :: ni_min( npd_aerosol_mode )
REAL :: ni_max( npd_aerosol_mode )
REAL :: ni_c( npd_aerosol_mode )
REAL :: ni_c_power( npd_aerosol_mode )
INTEGER, PARAMETER :: n_ni_fix = 1

! Local copies of mode type, component index and component type
INTEGER :: this_mode_type( npd_aerosol_mode )

! Loop variables
INTEGER :: i_mode ! loop on aerosol modes
INTEGER :: i_band ! loop on wavebands
INTEGER :: i_layr ! loop on vertical dimension
INTEGER :: i_prof ! loop on horizontal dimension

! Index for SSA array
INTEGER :: i_band_ssa

REAL :: logs_array_in(one)
REAL :: logs_array_out(one)
REAL :: incr_ni(npd_aerosol_mode)

REAL :: re_m( npd_profile, npd_layer, npd_band, npd_aerosol_mode )
REAL :: im_m( npd_profile, npd_layer, npd_band, npd_aerosol_mode )

REAL, PARAMETER :: min_ni_c = 0.001 ! Lowest value of ni_c to accept
REAL, PARAMETER :: max_ni_c = 5.0   ! Highest value of ni_c to accept
REAL, PARAMETER :: inv_ln_10 = 1.0 / LOG(10.0)

! Limits for the asymmetry parameter, since values of
! exactly -1.0 or +1.0 can cause div-by-zero errors
! further on in the Radiation code.
REAL, PARAMETER :: minus1_plus_epsi1 = -1.0 + EPSILON(1.0)
REAL, PARAMETER :: one_minus_epsi1 = 1.0 - EPSILON(1.0)

! Indicates whether current level is above the tropopause.
LOGICAL :: l_in_stratosphere( npd_profile, npd_layer )

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_BAND_AVERAGE_FASTER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)


DO i_mode = 1, n_ukca_mode

  ! Mode type. From a look-up table point of view, Aitken and
  ! accumulation types are treated in the same way.
  ! Accumulation soluble mode may use a narrower width (i.e. another
  ! look-up table) than other Aitken and accumulation modes.
  ! The coarse insoluble mode in the 3 dust mode setup may have a
  ! narrower width than the default 2.0 which is the case if the
  ! super-coarse insoluble mode is selected
  ! Once we know which look-up table to select, make local copies
  ! of info needed for nearest-neighbour calculations.
  !
  SELECT CASE (i_mode_type(i_mode))

  CASE (ip_ukca_mode_aitken)
    this_mode_type(i_mode) = ip_ukca_lut_accum

  CASE (ip_ukca_mode_accum)
    IF (l_soluble(i_mode)) THEN
      this_mode_type(i_mode) = ip_ukca_lut_accnarrow
    ELSE
      this_mode_type(i_mode) = ip_ukca_lut_accum
    END IF

  CASE (ip_ukca_mode_coarse)
    IF ((.NOT. l_soluble(i_mode)) .AND. l_cornarrow_ins) THEN
      this_mode_type(i_mode) = ip_ukca_lut_cornarrow
    ELSE
      this_mode_type(i_mode) = ip_ukca_lut_coarse
    END IF

  CASE (ip_ukca_mode_supercoarse)
    this_mode_type(i_mode) = ip_ukca_lut_supercoarse

  CASE DEFAULT
    ! Likely developer is trying to pass nucleation mode to radaer
    icode = 1
    cmessage = 'Mode is not one of aitken , accumulation , coarse'
    CALL ereport(RoutineName,icode,cmessage)
  END SELECT

END DO ! i_mode = 1, n_ukca_mode

DO i_mode = 1, n_ukca_mode
  
  nx(i_mode)      = ukca_lut(this_mode_type(i_mode), isolir)%n_x
  logxmin(i_mode) = LOG(ukca_lut(this_mode_type(i_mode), isolir)%x_min)
  logxmaxmlogxmin(i_mode) =                                                    &
                    LOG(ukca_lut(this_mode_type(i_mode), isolir)%x_max) -      &
                    LOG(ukca_lut(this_mode_type(i_mode), isolir)%x_min)

  nnr(i_mode)     = ukca_lut(this_mode_type(i_mode), isolir)%n_nr
  nrmin(i_mode)   = ukca_lut(this_mode_type(i_mode), isolir)%nr_min
  incr_nr(i_mode) = ukca_lut(this_mode_type(i_mode), isolir)%incr_nr

  nni(i_mode)     = ukca_lut(this_mode_type(i_mode), isolir)%n_ni
  ni_min(i_mode)  = ukca_lut(this_mode_type(i_mode), isolir)%ni_min
  ni_max(i_mode)  = ukca_lut(this_mode_type(i_mode), isolir)%ni_max
  ni_c(i_mode)    = ukca_lut(this_mode_type(i_mode), isolir)%ni_c
  ni_c_power(i_mode) = 10.0**( ukca_lut(this_mode_type(i_mode), isolir)%ni_c )

END DO ! i_mode = 1, n_ukca_mode

DO i_mode = 1, n_ukca_mode
  IF (ni_c > max_ni_c) THEN

    icode = 1
    cmessage='UKCA RADAER Look-up table'//newline//'NI_C exceeds upper limit'
    CALL ereport(RoutineName,icode,cmessage)

  END IF   
END DO

DO i_layr = 1, n_layer
  DO i_prof = 1, n_profile
    IF (l_inverted) THEN
      l_in_stratosphere(i_prof,i_layr) = i_layr <= trindxrad(i_prof)
    ELSE
      l_in_stratosphere(i_prof,i_layr) = i_layr >= trindxrad(i_prof)
    END IF
  END DO
END DO

DO i_mode = 1, n_ukca_mode
  DO i_layr = 1, n_layer
    DO i_prof = 1, n_profile
      DO i_cmpt = 1, n_cpnt_in_mode(i_mode)

        this_cpnt_type( i_cmpt,i_prof,i_layr,i_mode ) =                        &
                                                    i_cpnt_index(i_cmpt,i_mode)

      END DO ! i_cmpt
    END DO ! i_prof
  END DO ! i_layr
END DO ! i_mode
    
IF ( l_sustrat ) THEN

  DO i_mode = 1, n_ukca_mode
    DO i_layr = 1, n_layer
      DO i_prof = 1, n_profile
        DO i_cmpt = 1, n_cpnt_in_mode(i_mode)
          !
          ! If requested, switch the refractive index of the
          ! sulphate component to that for sulphuric acid
          ! for levels above the tropopause.
          !
          IF ( ( i_cpnt_type(this_cpnt) == cp_su ) .AND.                       &
               l_in_stratosphere(i_prof,i_layr) .AND.                          &
               ( .NOT. l_nitrate ) ) THEN

            this_cpnt_type( i_cmpt,i_prof,i_layr,i_mode ) = ip_ukca_h2so4

          END IF
        END DO ! i_cmpt
      END DO ! i_prof
    END DO ! i_layr
  END DO ! i_mode

END IF

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Part two - calculate the refractive index for real and optionally imaginary

re_m( i_prof, i_layr, i_band, i_mode ) = 0.0
im_m( i_prof, i_layr, i_band, i_mode ) = 0.0

! If single scattering albedo is prescribed, calculate both the real and
! imagniary components of refractive index
IF ( i_ukca_radaer_prescribe_ssa == do_not_prescribe ) THEN

  DO i_mode = 1, n_ukca_mode
    DO i_band = 1, n_band
      DO i_layr = 1, n_layer
        DO i_prof = 1, n_profile
            
          IF (ukca_modal_mmr   (i_prof, i_layr, i_mode) > threshold_mmr .AND.  &
              ukca_modal_number(i_prof, i_layr, i_mode) > threshold_nbr .AND.  &
              ukca_modal_volume(i_prof, i_layr, i_mode) > threshold_vol) THEN

            DO i_cmpt = 1, n_cpnt_in_mode(i_mode)

              ! Sum up refractive index, weighting by component volume
              re_m( i_prof, i_layr, i_band, i_mode ) =                         &
                   re_m( i_prof, i_layr, i_band, i_mode ) +                    &
                   ( ukca_cpnt_volume( i_cmpt, i_prof, i_layr ) *              &
                     precalc%realrefr( i_cmpt, one, i_band, isolir ) )

              ! Sum up refractive index, weighting by component volume
              im_m(i_prof,i_layr,i_band,i_mode) =                              &
                   im_m(i_prof,i_layr,i_band,i_mode) +                         &
                   ( ukca_cpnt_volume( i_cmpt, i_prof, i_layr ) *              &
                   precalc%imagrefr( i_cmpt, one, i_band, isolir ) )

            END DO ! i_cmpt

          END IF

        END DO ! i_prof
      END DO ! i_layr
    END DO ! i_band
  END DO ! i_mode

  DO i_mode = 1, n_ukca_mode
    IF ( l_soluble(i_mode) ) THEN
      DO i_band = 1, n_band
        DO i_layr = 1, n_layer
          DO i_prof = 1, n_profile

            IF (ukca_modal_mmr   (i_prof,i_layr,i_mode) > threshold_mmr .AND.  &
                ukca_modal_number(i_prof,i_layr,i_mode) > threshold_nbr .AND.  &
                ukca_modal_volume(i_prof,i_layr,i_mode) > threshold_vol) THEN

              ! Account for refractive index of water
              re_m(i_prof,i_layr,i_band,i_mode) =                              &
                      re_m(i_prof,i_layr,i_band,i_mode) +                      &
                      ( ukca_water_volume( i_prof, i_layr, i_mode ) *          &
                        precalc%realrefr(ip_ukca_water, one, i_band, isolir ) )

            END IF
            
          END DO ! i_prof
        END DO ! i_layr
      END DO ! i_band
    END IF ! l_soluble(i_mode)
  END DO ! i_mode

ELSE
! If single scattering albedo is not prescribed, calculate only the real
! component of refractive index

  DO i_mode = 1, n_ukca_mode
    DO i_band = 1, n_band
      DO i_layr = 1, n_layer
        DO i_prof = 1, n_profile

          IF (ukca_modal_mmr   (i_prof, i_layr, i_mode) > threshold_mmr .AND.  &
              ukca_modal_number(i_prof, i_layr, i_mode) > threshold_nbr .AND.  &
              ukca_modal_volume(i_prof, i_layr, i_mode) > threshold_vol) THEN

            DO i_cmpt = 1, n_cpnt_in_mode(i_mode)

              ! Sum up refractive index, weighting by component volume
              re_m( i_prof, i_layr, i_band, i_mode ) =                         &
                   re_m( i_prof, i_layr, i_band, i_mode ) +                    &
                   ( ukca_cpnt_volume( i_cmpt, i_prof, i_layr ) *              &
                     precalc%realrefr( i_cmpt, one, i_band, isolir ) )

            END DO ! i_cmpt

          END IF

        END DO ! i_prof
      END DO ! i_layr
    END DO ! i_band
  END DO ! i_mode

  DO i_mode = 1, n_ukca_mode
    IF ( l_soluble(i_mode) ) THEN
      DO i_band = 1, n_band
        DO i_layr = 1, n_layer
          DO i_prof = 1, n_profile

            IF (ukca_modal_mmr   (i_prof,i_layr,i_mode) > threshold_mmr .AND.  &
                ukca_modal_number(i_prof,i_layr,i_mode) > threshold_nbr .AND.  &
                ukca_modal_volume(i_prof,i_layr,i_mode) > threshold_vol) THEN

              ! Account for refractive index of water
              re_m(i_prof,i_layr,i_band,i_mode) =                              &
                      re_m(i_prof,i_layr,i_band,i_mode) +                      &
                      ( ukca_water_volume( i_prof, i_layr, i_mode ) *          &
                        precalc%realrefr(ip_ukca_water, one, i_band, isolir ) )

              im_m(i_prof,i_layr,i_band,i_mode) =                              &
                      im_m(i_prof,i_layr,i_band,i_mode) +                      &
                      ( ukca_water_volume( i_prof, i_layr, i_mode ) *          &
                        precalc%imagrefr(ip_ukca_water, one, i_band, isolir ) )

            END IF

          END DO ! i_prof
        END DO ! i_layr
      END DO ! i_band
    END IF ! l_soluble(i_mode)
  END DO ! i_mode
 
END IF       

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Part three - optionally obtain the imaginary component of SSA
!              and the nearest neighbour index of ni

DO i_mode = 1, n_ukca_mode

  a(i_mode) = ni_max(i_mode) / ( ni_c_power(i_mode) ) - 1.0 )
  b(i_mode) = REAL( nni(i_mode) ) / ni_c(i_mode)

  incr_ni(i_mode) = ( ni_max(i_mode) - ni_min(i_mode) ) / REAL(nni(i_mode)-1)

END DO

IF (i_ukca_radaer_prescribe_ssa == do_not_prescribe) THEN

  DO i_mode = 1, n_ukca_mode

    IF (ni_c(i_mode) > min_ni_c) THEN

      DO i_band = 1, n_band
        DO i_layr = 1, n_layer
          DO i_prof = 1, n_profile

            logs_array_in(one) = ( im_m( i_prof, i_layr, i_band, i_mode ) /    &
                                   a( i_mode ) ) + 1.0

            CALL log_v( one, logs_array_in(one), logs_array_out(one) )

            ni_ind( i_prof, i_layr, i_band, i_mode ) =                         &
                        NINT( b(i_mode) * logs_array_out(one) * inv_ln_10 ) + 1

          END DO ! i_prof
        END DO ! i_layr
      END DO ! i_band

   ELSE ! (ni_c(i_mode) > min_ni_c) THEN

      DO i_band = 1, n_band
        DO i_layr = 1, n_layer
          DO i_prof = 1, n_profile

            ni_ind( i_prof, i_layr, i_band, i_mode ) =                         &
                 NINT( ( im_m( i_prof, i_layr, i_band, i_mode ) -              &
                         ni_min( i_mode ) ) /                                  &
                       incr_ni( i_mode ) ) + 1

          END DO ! i_prof
        END DO ! i_layr
      END DO ! i_band

    END IF ! (ni_c(i_mode) > min_ni_c) THEN

  END DO ! i_mode

  DO i_mode = 1, n_ukca_mode
    DO i_band = 1, n_band
      DO i_layr = 1, n_layer
        DO i_prof = 1, n_profile

          ni_ind( i_prof, i_layr, i_band, i_mode ) =                           &
             MIN( nni(i_mode), MAX( 1, ni_indi_prof, i_layr, i_band, i_mode ) )

        END DO ! i_prof
      END DO ! i_layr
    END DO ! i_band
  END DO ! i_mode

END IF

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Part four - Compute the Mie parameter

IF (i_ukca_radaer_prescribe_ssa /= do_not_prescribe) THEN
  DO i_mode = 1, n_ukca_mode
    DO i_band = 1, n_band
      DO i_layr = 1, n_layer
        DO i_prof = 1, n_profile

          ! Fix index to one for prescribed SSA
          ni_ind( i_prof, i_layr, i_band, i_mode ) = n_ni_fix

        END DO ! i_prof
      END DO ! i_layr
    END DO ! i_band
  END DO ! i_mode
END IF

          
  DO i_mode = 1, n_ukca_mode
    DO i_band = 1, n_band
      DO i_layr = 1, n_layer
        DO i_prof = 1, n_profile

          ! Compute the Mie parameter from the wet diameter
          ! and get the LUT-array index of its nearest neighbour.
          x = pi * ukca_wet_diam( i_prof, i_layr, i_mode ) /                   &
                   precalc%wavelength( one, i_band, isolir )

          n_x = NINT( ( LOG(x) - logxmin ) / logxmaxmlogxmin * (nx-1) ) + 1

          n_x = MIN( nx, MAX( 1, n_x ) )

          ! Same for the dry diameter (needed to access the volume fraction)
          x_dry = pi * ukca_dry_diam( i_prof, i_layr, i_mode ) /               &
                       precalc%wavelength( one, i_band, isolir )

          n_x_dry = NINT( ( LOG(x_dry) - logxmin ) /                           &
                            logxmaxmlogxmin * (nx-1) ) + 1

          n_x_dry = MIN( nx, MAX( 1, n_x_dry ) )

          ! Compute the modal complex refractive index as
          ! volume-weighted component refractive indices.
          ! Get the LUT-array index of their nearest neighbours.
          n_nr( i_prof, i_layr, i_band, i_mode ) =                             &
                                 NINT( ( re_m(i_intg) - nrmin ) / incr_nr ) + 1

          n_nr( i_prof, i_layr, i_band, i_mode ) = &
               MIN( nnr, MAX( 1, n_nr( i_prof, i_layr, i_band, i_mode ) ) )

          ! Get local copies of the relevant look-up table entries.
          loc_sca( i_prof, i_layr, i_band, i_mode ) =                          &
                              ukca_lut(this_mode_type, isolir)%                &
               ukca_scattering( n_x, ni_ind( i_prof, i_layr, i_band, i_mode ), &
                                n_nr )

          loc_asy( i_prof, i_layr, i_mode, i_band ) =                          &
                              ukca_lut(this_mode_type, isolir)%                &
               ukca_asymmetry(  n_x, ni_ind( i_prof, i_layr, i_band, i_mode ), &
               n_nr )

          loc_vol = ukca_lut(this_mode_type, isolir)%                        &
               volume_fraction( n_x_dry )

          factor( i_prof, i_layr, i_mode, i_band ) = 1.0 /                     &
                   ( ukca_modal_density( i_prof, i_layr, i_mode) *         &
                     loc_vol *                                             &
                     precalc%wavelength( 1 , i_band, isolir) )


        END DO ! i_prof
      END DO ! i_layr
    END DO ! i_band
  END DO ! i_mode


IF ( i_ukca_radaer_prescribe_ssa == do_not_prescribe) THEN

  DO i_mode = 1, n_ukca_mode
    DO i_band = 1, n_band
      DO i_layr = 1, n_layer
        DO i_prof = 1, n_profile

          IF (ukca_modal_mmr   (i_prof,i_layr,i_mode) > threshold_mmr .AND.    &
              ukca_modal_number(i_prof,i_layr,i_mode) > threshold_nbr .AND.    &
              ukca_modal_volume(i_prof,i_layr,i_mode) > threshold_vol) THEN
           
             loc_abs = ukca_lut( this_mode_type(i_mode), isolir )%             &
                ukca_absorption(    n_x( i_prof, i_layr, i_band, i_mode ),     &
                                 ni_ind( i_prof, i_layr, i_band, i_mode ),     &
                                 n_nr( i_mode ) )

             ukca_absorption( i_prof, i_layr, i_mode, i_band ) = MAX( 0.0,     &
                           loc_abs * factor( i_prof, i_layr, i_mode, i_band ) )

             ukca_scattering( i_prof, i_layr, i_mode, i_band ) = MAX( 0.0,     &
                                   loc_sca( i_prof, i_layr, i_mode, i_band ) * &
                                   factor( i_prof, i_layr, i_mode, i_band ) )

          ELSE ! Below threshold

            ukca_absorption( i_prof, i_layr, i_mode, i_band ) = 0.0

            ukca_scattering( i_prof, i_layr, i_mode, i_band ) = 0.0
             
          END IF

        END DO ! i_prof
      END DO ! i_layr
    END DO ! i_band
  END DO ! i_mode

ELSE ! ( i_ukca_radaer_prescribe_ssa == do_not_prescribe)

  DO i_mode = 1, n_ukca_mode
    DO i_band = 1, n_band
      DO i_layr = 1, n_layer
        DO i_prof = 1, n_profile

          IF (ukca_modal_mmr   (i_prof,i_layr,i_mode) > threshold_mmr .AND.    &
              ukca_modal_number(i_prof,i_layr,i_mode) > threshold_nbr .AND.    &
              ukca_modal_volume(i_prof,i_layr,i_mode) > threshold_vol) THEN
           
            i_band_ssa = MIN( npd_band_ssa, i_band )

            this_ssa = ukca_radaer_presc_ssa( i_prof, i_layr, i_band_ssa )

            ukca_absorption( i_prof, i_layr, i_mode, i_band ) = MAX( 0.0,      &
               loc_sca( i_prof, i_layr, i_mode, i_band ) *                     &
               factor( i_prof, i_layr, i_mode, i_band ) *                      &
               ( 1.0 - this_ssa ) )

            ukca_scattering( i_prof, i_layr, i_mode, i_band ) = MAX( 0.0,      &
               loc_sca( i_prof, i_layr, i_mode, i_band ) *                     &
               factor( i_prof, i_layr, i_mode, i_band ) *                      &
               this_ssa )

          ELSE ! Below threshold

            ukca_absorption( i_prof, i_layr, i_mode, i_band ) = 0.0

            ukca_scattering( i_prof, i_layr, i_mode, i_band ) = 0.0

          END IF
            
        END DO ! i_prof
      END DO ! i_layr
    END DO ! i_band
  END DO ! i_mode           

END IF ! ( i_ukca_radaer_prescribe_ssa == do_not_prescribe)


DO i_mode = 1, n_ukca_mode
  DO i_band = 1, n_band
    DO i_layr = 1, n_layer
      DO i_prof = 1, n_profile

        IF (ukca_modal_mmr   (i_prof,i_layr,i_mode) > threshold_mmr .AND.    &
            ukca_modal_number(i_prof,i_layr,i_mode) > threshold_nbr .AND.    &
            ukca_modal_volume(i_prof,i_layr,i_mode) > threshold_vol) THEN
           
           ukca_asymmetry( i_prof, i_layr, i_mode, i_band ) =                  &
                MAX( minus1_plus_epsi1, MIN( one_minus_epsi1,                  &
                loc_asy( i_prof, i_layr, i_mode, i_band ) ) )

         ELSE ! Below threshold

           ukca_asymmetry( i_prof, i_layr, i_mode, i_band ) = 0.0

         END IF
                
      END DO ! i_prof
    END DO ! i_layr
  END DO ! i_band
END DO ! i_mode  

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_band_average_faster

END MODULE ukca_radaer_band_average_faster_mod
