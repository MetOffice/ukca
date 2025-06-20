! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! UKCA_RADAER pre-computed values.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
! Those values depend only on the definitions of SW and LW
! wavebands, and on the UKCA look-up tables.
!
MODULE ukca_radaer_precalc

IMPLICIT NONE
SAVE

!
! Maximum number of integration points.
!
INTEGER, PARAMETER :: npd_integ_pts = 6

!
! Maximum number of wavebands in a spectrum.
! npd_ukca_band must be as large as the number of SW or LW bands
! (whichever is the largest), as given in the spectral files.
!
INTEGER, PARAMETER :: npd_ukca_band = 9

!
! Number of spectra (SW and LW, typically).
! Must be equal to npd_ukca_lut_spectrum in module ukca_radaer_lut_mod.
!
INTEGER, PARAMETER :: npd_ukca_spectrum = 2

!
! Maximum number of component types. Must be equal to the largest
! component number (typically, ip_ukca_water) defined in
! ukca_radaer_struct_mod.
!
INTEGER, PARAMETER :: npd_ukca_maxcomptype = 11

!
! Maximum number of wavelengths for the aerosol optical depth
!
INTEGER, PARAMETER :: npd_ukca_aod_wavel = 6

TYPE :: ukca_radaer_tprecalc

  INTEGER :: n_integ_pts

  REAL :: wavelength(npd_integ_pts, npd_ukca_band, npd_ukca_spectrum)

  REAL :: irrad(npd_integ_pts, npd_ukca_band, npd_ukca_spectrum)

  REAL :: weight(npd_integ_pts)

  REAL :: flux(npd_ukca_band, npd_ukca_spectrum)

  !
  ! Component refractive indices.
  !
  REAL :: realrefr(npd_ukca_maxcomptype, npd_integ_pts,                        &
                  npd_ukca_band, npd_ukca_spectrum)
  REAL :: imagrefr(npd_ukca_maxcomptype, npd_integ_pts,                        &
                  npd_ukca_band, npd_ukca_spectrum)

  !
  ! Wavelengths for the aerosol optical depth (metre).
  !
  INTEGER :: n_ukca_aod_wavel

  REAL :: aod_wavel(npd_ukca_aod_wavel)

  !
  ! Component refractive index for the aerosol optical depth
  !
  REAL ::  aod_realrefr(npd_ukca_maxcomptype, npd_ukca_aod_wavel)
  REAL ::  aod_imagrefr(npd_ukca_maxcomptype, npd_ukca_aod_wavel)

END TYPE ukca_radaer_tprecalc

TYPE (ukca_radaer_tprecalc) :: precalc

END MODULE ukca_radaer_precalc
