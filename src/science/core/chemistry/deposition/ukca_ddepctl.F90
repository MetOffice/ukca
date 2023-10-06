! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Purpose: To assign dry deposition rates in s-1 to array dpd in
!           ASAD module asad_mod. This is passed for use in UKCA chemistry.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!           Called from UKCA_CHEMISTRY_CTRL.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_ddepctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_DDEPCTL_MOD'

CONTAINS

SUBROUTINE ukca_ddepctl(row_length, rows, bl_levels, ntype,                    &
  npft, land_points, land_index, tile_pts, tile_index, timestep,               &
  sinlat, frac_types, t_surf, p_surf, dzl, zbl, surf_hf, u_s,                  &
  rh, stcon, soilmc_lp, fland, seaice_frac, laift_lp,                          &
  canhtft_lp, z0tile_lp, t0tile_lp,                                            &
  nlev_with_ddep, zdryrt, len_stashwork, stashwork)

USE ukca_config_specification_mod, ONLY: ukca_config
USE asad_mod, ONLY: jpdd
USE water_constants_mod, ONLY: tfs, rho_water
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

USE ukca_aerod_mod, ONLY: ukca_aerod
USE ukca_ddcalc_mod, ONLY: ukca_ddcalc
USE ukca_surfddr_mod, ONLY: ukca_surfddr

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: bl_levels
INTEGER, INTENT(IN) :: land_points
INTEGER, INTENT(IN) :: ntype
INTEGER, INTENT(IN) :: npft
INTEGER, INTENT(IN) :: land_index(land_points)

INTEGER, INTENT(IN) :: tile_pts(ntype)
INTEGER, INTENT(IN) :: tile_index(land_points,ntype)

REAL, INTENT(IN) :: timestep

REAL, INTENT(IN) :: sinlat(row_length,rows)
REAL, INTENT(IN) :: frac_types(:,:)
REAL, INTENT(IN) :: t_surf(row_length,rows)
REAL, INTENT(IN) :: p_surf(row_length,rows)
REAL, INTENT(IN) :: dzl(row_length,rows,bl_levels)
REAL, INTENT(IN) :: zbl(row_length,rows)
REAL, INTENT(IN) :: surf_hf(row_length,rows)
REAL, INTENT(IN) :: u_s(row_length,rows)
REAL, INTENT(IN) :: rh(row_length,rows)
REAL, INTENT(IN) :: seaice_frac(row_length,rows)
REAL, INTENT(IN) :: stcon(row_length,rows,npft)
REAL, INTENT(IN) :: soilmc_lp(land_points)
REAL, INTENT(IN) :: fland(land_points)
REAL, INTENT(IN) :: laift_lp(land_points,npft)
REAL, INTENT(IN) :: canhtft_lp(land_points,npft)
REAL, INTENT(IN) :: z0tile_lp(land_points,ntype)
REAL, INTENT(IN) :: t0tile_lp(land_points,ntype)

! Diagnostics array
INTEGER, INTENT(IN)  :: len_stashwork
REAL, INTENT(IN OUT) :: stashwork (len_stashwork)

INTEGER, INTENT(OUT) :: nlev_with_ddep(row_length,rows)

REAL, INTENT(OUT) :: zdryrt(row_length,rows,jpdd)

!     Local variables

INTEGER :: i, j, k, l, m, n

INTEGER :: lake ! Surface type index for type 'lake'

REAL :: seafrac ! Fraction of sea in grid square
REAL :: lftotal ! Sum of land tile fractions

REAL :: smr(row_length,rows)           ! soil moisture as volume fraction
REAL :: land_fraction(row_length,rows)
REAL :: lai_ft(row_length,rows,npft)
REAL :: canht_ft(row_length,rows,npft)
REAL :: z0tile(row_length,rows,ntype)
REAL :: t0tile(row_length,rows,ntype)
REAL :: o3_stom_frac(row_length,rows)

! Aerodynamic resistance (s m-1)
REAL :: resa(row_length,rows,ntype)

! Quasi-laminar resistance (s m-1)
REAL :: rb(row_length,rows,jpdd)

! Surface resistance (s m-1)
REAL :: rc(row_length,rows,ntype,jpdd)

! Global surface fraction array
REAL :: gsf(row_length,rows,ntype)

! Sulphate aerosol dep velocity (m s-1)
REAL:: so4_vd(row_length,rows)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DDEPCTL'

!
! Set up global surface fraction array, gsf
! Tile fractions (in array frac_types) are only defined over land.
! Set water fractions over the sea. Add sea ice and sea fractions to gsf.
! Adjust tile fractions so they add up to 1. If gsf has the water tile
! fraction = 1.0, set the land fraction to 0.0
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
land_fraction = 0.0
smr           = 0.0
lai_ft        = 0.0
canht_ft      = 0.0
z0tile        = 0.0
t0tile        = 0.0
gsf           = 0.0
!
! Expand arrays from land points to lon-lat grid
!
DO l = 1, land_points
  j = (land_index(l)-1)/row_length + 1
  i = land_index(l) - (j-1)*row_length
  smr(i,j) = soilmc_lp(l)/(ukca_config%dzsoil_layer1*rho_water)
  land_fraction(i,j) = fland(l)
END DO
!
DO m = 1, npft
  DO n = 1, tile_pts(m)
    l = tile_index(n,m)
    j = (land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    lai_ft(i,j,m) = laift_lp(l,m)
    canht_ft(i,j,m) = canhtft_lp(l,m)
  END DO
END DO
!
DO m = 1, ntype
  DO n = 1, tile_pts(m)
    l = tile_index(n,m)
    j = (land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    z0tile(i,j,m) = z0tile_lp(l,m)
    t0tile(i,j,m) = t0tile_lp(l,m)
    gsf(i,j,m) = frac_types(l,m)
  END DO
END DO
!
! Add in sea ice data and coastal land fractions to gsf array.
!
lake = ukca_config%i_lake
DO k = 1, rows
  DO i = 1, row_length
    IF (gsf(i,k,lake) == 1.0) land_fraction(i,k) = 0.0
    seafrac = 1.0 - land_fraction(i,k)
    IF (land_fraction(i,k) < 1.0 .AND. gsf(i,k,lake) < 1.0) THEN
      gsf(i,k,lake) = seafrac
      IF (land_fraction(i,k) > 0.0) THEN
        lftotal = 0.0
        n = 0
        DO j = 1, ntype
          n = n + 1
          IF (j == lake) CYCLE      ! ignore water
          lftotal = lftotal + gsf(i,k,n)
        END DO
        n = 0
        DO j = 1, ntype
          n = n + 1
          IF (j == lake) CYCLE      ! ignore water
          gsf(i,k,n) = gsf(i,k,n) * land_fraction(i,k) / lftotal
        END DO
      END IF
    END IF
    IF (seaice_frac(i,k) > 0.0) THEN
      gsf(i,k,lake) = (1.0 - seaice_frac(i,k)) * seafrac
      gsf(i,k,ntype) = gsf(i,k,ntype) + seaice_frac(i,k) * seafrac

      ! The equation above, involving sea ice, implicitly assumes that
      ! the last surface type is ice. This is valid for 3 out of the 4 allowed
      ! UKCA surface tile configurations, 5 pft/9 tile, 9 pft/13 tile and
      ! 13 pft/17 tile, as ice is usually the last surface type. However, it is
      ! not valid for the fourth UKCA surface tile configuration
      ! (13 pft/27 tile), where the last surface type is element 10 of
      ! elevated ice.

      ! To maintain outputs and comparisons with known good output for UKESM1
      ! and UKESM1.1, the UKCA code is left unchanged and this comment added.
      ! The equation has been corrected in the equivalent JULES-based
      ! deposition routines, with ntype replaced by ice. This will ensure
      ! that the allocation of sea ice to ice is correct in all surface tile
      ! configurations, especially those where ice is not the last surface
      ! type, e.g., in UKESM2, where element 10 of elevated rock will be
      ! the last surface type.

    END IF
  END DO
END DO
!
! Set all tile temperatures to t0 where undefined.
!
DO n = 1, ntype
  DO k = 1, rows
    DO i = 1, row_length
      IF (t0tile(i,k,n) < 100.0) t0tile(i,k,n) = t_surf(i,k)
    END DO
  END DO
END DO
!
! Set up tile temperatures where a mixture of sea and sea ice is
! present. Set sea to freezing temperature (tfs) and ice to sea ice
! temperature.
!
DO k = 1, rows
  DO i = 1, row_length
    IF (seaice_frac(i,k) > 0.0) THEN
      !           t0tile(i,k,ntype) = sicetemp(i,k)
      t0tile(i,k,lake) = tfs
    END IF
  END DO
END DO
!
! Calculate aerodynamic and quasi-laminar resistances (Resa, Rb).
!
CALL ukca_aerod(row_length, rows, ntype, npft, lake, t_surf, p_surf,           &
  surf_hf, u_s, canht_ft, gsf, zbl, z0tile, resa, rb, so4_vd)
!
! Calculate surface resistance term Rc
!
CALL ukca_surfddr(row_length, rows, ntype, npft,                               &
        sinlat, t_surf, p_surf, rh, smr, u_s, gsf, stcon, t0tile,              &
        lai_ft, so4_vd, rc, o3_stom_frac)
!
! Combine resistance terms to calculate overall dry deposition
! velocity, and hence first-order rate constant describing
! dry deposition rate
!
CALL ukca_ddcalc(row_length, rows, bl_levels, ntype,                           &
  npft, timestep, dzl, zbl, gsf, resa, rb, rc,                                 &
  nlev_with_ddep, zdryrt, len_stashwork, stashwork)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
!
END SUBROUTINE ukca_ddepctl
END MODULE ukca_ddepctl_mod
