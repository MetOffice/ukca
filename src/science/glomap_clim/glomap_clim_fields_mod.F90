! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!   To put fields required by RADAER into stash.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
MODULE glomap_clim_fields_mod

USE um_types,           ONLY: real_umphys

IMPLICIT NONE

REAL(KIND=real_umphys), ALLOCATABLE, PUBLIC :: aird(:)
! Air density number concentration (cm^-3) (n_points)

REAL(KIND=real_umphys), ALLOCATABLE, PUBLIC :: drydp(:,:)
! Median particle dry diameter for each mode (m)      (n_points,nmodes)

REAL(KIND=real_umphys), ALLOCATABLE, PUBLIC :: dvol(:,:)
! Median particle dry volume for each mode (m^3)      (n_points,nmodes)

REAL(KIND=real_umphys), ALLOCATABLE, PUBLIC :: md(:,:,:)
! Avg cpt mass of aerosol particle in mode (particle^-1) (n_points,nmodes,ncp)

REAL(KIND=real_umphys), ALLOCATABLE, PUBLIC :: mdt(:,:)
! Avg tot mass of aerosol ptcl in mode (particle^-1)  (n_points,nmodes)

REAL(KIND=real_umphys), ALLOCATABLE, PUBLIC :: mmr1d(:,:,:)
! Avg cpt mass mixing ratio of aerosol particle in mode (particle^-1)
!                                                       (n_points,nmodes,ncp)

REAL(KIND=real_umphys), ALLOCATABLE, PUBLIC :: nd(:,:)
! Aerosol ptcl number density for mode (cm^-3) (n_points,nmodes)

REAL(KIND=real_umphys), ALLOCATABLE, PUBLIC :: nmr1d(:,:)
! Aerosol ptcl (number density/ air density) for mode (cm^-3)
!                                                     (n_points,nmodes)

REAL(KIND=real_umphys), ALLOCATABLE, PUBLIC :: pmid_1d(:)
! Pressure at mid levels (Pa)

REAL(KIND=real_umphys), ALLOCATABLE, PUBLIC :: temp_1d(:)
! Temperature 1-D

END MODULE glomap_clim_fields_mod
