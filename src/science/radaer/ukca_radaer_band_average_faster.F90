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

! USE

IMPLICIT NONE

!
! Arguments with intent(in)
!

! INTEGER, INTENT(IN) :: blah

!
! Local variables
!
INTEGER, PARAMETER :: one = 1


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
