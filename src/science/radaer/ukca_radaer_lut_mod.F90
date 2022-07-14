! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! UKCA_RADAER look-up tables.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
! There are six look-up tables, for wide accumulation mode, narrow
! accumulation mode, and coarse mode aerosols (first dimension) and
! for SW and LW spectra (second dimension).
!
MODULE ukca_radaer_lut

USE ukca_radaer_tlut_mod, ONLY: ukca_radaer_tlut

IMPLICIT NONE
SAVE

INTEGER, PARAMETER :: npd_ukca_lut_mode     = 3
INTEGER, PARAMETER :: ip_ukca_lut_accum     = 1
INTEGER, PARAMETER :: ip_ukca_lut_coarse    = 2
INTEGER, PARAMETER :: ip_ukca_lut_accnarrow = 3

INTEGER, PARAMETER :: npd_ukca_lut_spectrum = 2
!                 those two values must be consistent with
!                 the parameters in spcrg3a_mod
INTEGER, PARAMETER :: ip_ukca_lut_sw        = 1
INTEGER, PARAMETER :: ip_ukca_lut_lw        = 2

TYPE (ukca_radaer_tlut) :: ukca_lut(npd_ukca_lut_mode, npd_ukca_lut_spectrum)

END MODULE ukca_radaer_lut
