! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Populates the ukca_lut structure from information previously obtained
!  from a namelist read
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM

MODULE ukca_radaer_populate_lut_mod

USE errormessagelength_mod,  ONLY:                                             &
    errormessagelength

USE parkind1,                ONLY:                                             &
    jpim, jprb

USE ukca_radaer_lut,         ONLY:                                             &
    ukca_lut

USE ukca_option_mod,         ONLY:                                             &
    l_ukca_radaer_prescribe_ssa

USE yomhook,                 ONLY:                                             &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
                                      ModuleName='UKCA_RADAER_POPULATE_LUT_MOD'

CONTAINS

SUBROUTINE ukca_radaer_populate_lut( aerosol_band, wavelength_band,            &
                                     n_x, n_nr, n_ni,                          &
                                     npd_x, npd_nr, npd_ni,                    &
                                     stdev, x_min, x_max, nr_min,              &
                                     nr_max, ni_min, ni_max, ni_c,             &
                                     ukca_absorption,                          &
                                     ukca_scattering,                          &
                                     ukca_asymmetry,                           &
                                     volume_fraction,                          &
                                     icode, cmessage )



IMPLICIT NONE

!
! Arguments
!

! This is one of ip_ukca_lut_accum , ip_ukca_lut_coarse , ip_ukca_lut_accnarrow,
!                ip_ukca_lut_cornarrow, ip_ukca_lut_supercoarse
INTEGER, INTENT(IN) :: aerosol_band

! This is one of ip_ukca_lut_sw , ip_ukca_lut_lw
INTEGER, INTENT(IN) :: wavelength_band

INTEGER, INTENT(IN OUT) :: n_x
INTEGER, INTENT(IN OUT) :: n_nr
INTEGER, INTENT(IN OUT) :: n_ni
INTEGER, INTENT(IN) :: npd_x
INTEGER, INTENT(IN) :: npd_nr
INTEGER, INTENT(IN) :: npd_ni
REAL, INTENT(IN)    :: stdev
REAL, INTENT(IN)    :: x_min
REAL, INTENT(IN)    :: x_max
REAL, INTENT(IN)    :: nr_min
REAL, INTENT(IN)    :: nr_max
REAL, INTENT(IN)    :: ni_min
REAL, INTENT(IN)    :: ni_max
REAL, INTENT(IN)    :: ni_c
REAL, INTENT(IN)    :: ukca_absorption(0:npd_x, 0:npd_ni, 0:npd_nr)
REAL, INTENT(IN)    :: ukca_scattering(0:npd_x, 0:npd_ni, 0:npd_nr)
REAL, INTENT(IN)    :: ukca_asymmetry(0:npd_x, 0:npd_ni, 0:npd_nr)
REAL, INTENT(IN)    :: volume_fraction(0:npd_x)

! Error reporting arguments
INTEGER, INTENT(IN OUT) :: icode
CHARACTER(LEN=errormessagelength), INTENT(OUT) :: cmessage


! Local variables

INTEGER :: i, j, k


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*),   PARAMETER :: RoutineName='UKCA_RADAER_POPULATE_LUT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialise error code
icode = 0

!
! Check that actual array dimensions do not exceed
! those that were fixed at compile time.
!
IF (n_x > npd_x) THEN
  icode = 1
  cmessage='Look-up table dimension exceeds built-in limit (X).'
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

IF (n_nr > npd_nr) THEN
  icode = 1
  cmessage='Look-up table dimension exceeds built-in limit (NR).'
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

IF (n_ni > npd_ni) THEN
  icode = 1
  cmessage='Look-up table dimension exceeds built-in limit (NI).'
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Copy the look-up table into the structure passed as an argument
ukca_lut(aerosol_band,wavelength_band)%stdev = stdev

!
! Namelist arrays indices start at 0. For the sake of ease of use,
! online arrays start at 1 -> simply shift the content by 1 and add
! 1 to each dimension.
! All calculations within UKCA_RADAER expect array indexing starting
! at 1.
!
! If l_ukca_radaer_prescribe_ssa is .true., then the n_ni dimension
! becomes 1 and only the first LUT element in that dimension, which
! corresponds to a purely scattering case, is read in. Absorption
! coefficients are not stored, since they are all zero.
!
n_x = n_x + 1
n_nr = n_nr + 1
IF (l_ukca_radaer_prescribe_ssa) THEN
  n_ni = 1
ELSE
  n_ni = n_ni + 1
END IF

ukca_lut(aerosol_band,wavelength_band)%n_x  = n_x
ukca_lut(aerosol_band,wavelength_band)%n_nr = n_nr
ukca_lut(aerosol_band,wavelength_band)%n_ni = n_ni

ukca_lut(aerosol_band,wavelength_band)%x_min  = x_min
ukca_lut(aerosol_band,wavelength_band)%x_max  = x_max
ukca_lut(aerosol_band,wavelength_band)%nr_min = nr_min
ukca_lut(aerosol_band,wavelength_band)%nr_max = nr_max

! Check that we are not dividing by zero
IF (n_nr > 1) THEN
  ukca_lut(aerosol_band,wavelength_band)%incr_nr =                             &
    (nr_max - nr_min) / REAL(n_nr-1)
ELSE
  ukca_lut(aerosol_band,wavelength_band)%incr_nr = 0.0
END IF

IF (l_ukca_radaer_prescribe_ssa) THEN
  ukca_lut(aerosol_band,wavelength_band)%ni_min = 0.0
  ukca_lut(aerosol_band,wavelength_band)%ni_max = 0.0
  ukca_lut(aerosol_band,wavelength_band)%ni_c   = 0.0
ELSE
  ukca_lut(aerosol_band,wavelength_band)%ni_min = ni_min
  ukca_lut(aerosol_band,wavelength_band)%ni_max = ni_max
  ukca_lut(aerosol_band,wavelength_band)%ni_c   = ni_c
END IF

! Allocate the dynamic arrays
IF (l_ukca_radaer_prescribe_ssa) THEN

  ALLOCATE(ukca_lut(aerosol_band,wavelength_band)%ukca_absorption(1,1,1))
                                               ! Unused, so smallest allocation
  ALLOCATE(ukca_lut(aerosol_band,                                              &
                    wavelength_band)%ukca_scattering(n_x,n_ni,n_nr))
  ALLOCATE(ukca_lut(aerosol_band,                                              &
                    wavelength_band)%ukca_asymmetry( n_x,n_ni,n_nr))
  ALLOCATE(ukca_lut(aerosol_band,wavelength_band)%volume_fraction(n_x))

ELSE

  ALLOCATE(ukca_lut(aerosol_band,                                              &
                    wavelength_band)%ukca_absorption(n_x,n_ni,n_nr))
  ALLOCATE(ukca_lut(aerosol_band,                                              &
                    wavelength_band)%ukca_scattering(n_x,n_ni,n_nr))
  ALLOCATE(ukca_lut(aerosol_band,                                              &
                    wavelength_band)%ukca_asymmetry( n_x,n_ni,n_nr))
  ALLOCATE(ukca_lut(aerosol_band,wavelength_band)%volume_fraction(n_x))

END IF


! Populate the dynamic arrays
IF (l_ukca_radaer_prescribe_ssa) THEN

  DO k = 1, n_nr
    DO j = 1, n_ni
      DO i = 1, n_x
        ukca_lut(aerosol_band,wavelength_band)%ukca_scattering(i,j,k) =        &
            ukca_scattering( i-1, j-1, k-1 )
        ukca_lut(aerosol_band,wavelength_band)%ukca_asymmetry(i,j,k)  =        &
            ukca_asymmetry(  i-1, j-1, k-1 )
      END DO
    END DO
  END DO

ELSE

  DO k = 1, n_nr
    DO j = 1, n_ni
      DO i = 1, n_x
        ukca_lut(aerosol_band,wavelength_band)%ukca_absorption(i,j,k) =        &
            ukca_absorption( i-1, j-1, k-1 )
        ukca_lut(aerosol_band,wavelength_band)%ukca_scattering(i,j,k) =        &
            ukca_scattering( i-1, j-1, k-1 )
        ukca_lut(aerosol_band,wavelength_band)%ukca_asymmetry(i,j,k)  =        &
            ukca_asymmetry(  i-1, j-1, k-1 )
      END DO
    END DO
  END DO

END IF

DO i = 1, n_x
  ukca_lut(aerosol_band,wavelength_band)%volume_fraction(i) =                  &
                                         volume_fraction(i-1)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_populate_lut
END MODULE ukca_radaer_populate_lut_mod
