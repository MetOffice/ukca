! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module defines the elements of the structure
! gathering the contents of a UKCA look-up table
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
MODULE ukca_radaer_tlut_mod

IMPLICIT NONE

TYPE :: ukca_radaer_tlut

  !
  ! Standard deviation of the lognormal size
  ! distribution (differs between accumulation-
  ! and coarse-mode aerosols).
  !
  REAL :: stdev

  !
  ! Array dimensions
  ! (resp. Mie parameter, real part of the
  !        refractive index, imaginary part
  !        of the refractive index)
  !
  INTEGER :: n_x
  INTEGER :: n_nr
  INTEGER :: n_ni

  !
  ! Min and max values for each dimension
  !
  REAL :: x_min
  REAL :: x_max
  REAL :: nr_min
  REAL :: nr_max
  REAL :: ni_min
  REAL :: ni_max
  REAL :: ni_c

  !
  ! Increments for linearly-varying dimensions
  !
  REAL :: incr_nr

  !
  ! Look-up tables themselves
  ! (resp. absorption and scattering coefficients,
  !        asymmetry parameter and volume fraction)
  ! Those will be dynamically allocated.
  !
  REAL, POINTER :: ukca_absorption(:, :, :)
  REAL, POINTER :: ukca_scattering(:, :, :)
  REAL, POINTER :: ukca_asymmetry(:, :, :)
  REAL, POINTER :: volume_fraction(:)

END TYPE ukca_radaer_tlut

END MODULE ukca_radaer_tlut_mod
