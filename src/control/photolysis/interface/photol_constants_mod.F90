! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module providing constant values needed by the photolysis component.
!
!   These are currently initialised to UM values in photol_init_configuration
!   but some of the constants are configurable by parent models to have
!   consistent values across the components.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_Photolysis
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE photol_constants_mod

USE missing_data_mod, ONLY: rmdi

IMPLICIT NONE

PUBLIC

! --- Constants used by Photolysis ---
! These will be initialised to current UM values but available to be
! overwritten by any parent model if necessary.
! These have a 'const_' prefix to avoid conflict with original variables in
! the codebase, but aliased in photol routines where these are actually used.

REAL :: const_pi = rmdi               ! Pi
REAL :: const_pi_over_180 = rmdi      ! pi / 180. deg
REAL :: const_recip_pi_over_180 = rmdi ! reciprocal of pi / 180.
REAL :: const_rhour_per_day = rmdi   ! Hours in a day

REAL :: const_o3_mmr_vmr = rmdi      ! Factor to convert ozone between
                                     ! mass mix ratio to volume mix ratio
REAL :: const_molemass_sulp = rmdi   ! Molecular mass of sulphur (g/mol)
REAL :: const_molemass_nh42so4 = rmdi  ! Molecular mass of ammonium sulphate
REAL :: const_molemass_air = rmdi    ! Molecular mass of dry air (kg/mol)
REAL :: const_planet_radius = rmdi   ! Planet radius (m)
REAL :: const_g = rmdi               ! Planet accel due to gravity (m s-2)
REAL :: const_r = rmdi               ! Gas constant for dry air
REAL :: const_avogadro = rmdi        ! Number of molecules per mole
REAL :: const_tm = rmdi              ! Freezing point - water (melting pt -ice)

REAL :: const_s2r = rmdi             ! Seconds-to-radians converter,
                                     ! calculated as planet_dha/rsec_per_day

END MODULE photol_constants_mod
