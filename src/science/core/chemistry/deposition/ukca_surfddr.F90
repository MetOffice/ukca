! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_surfddr_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_SURFDDR_MOD'

CONTAINS

SUBROUTINE ukca_surfddr( row_length, rows, ntype, npft, sinlat, t0, p0, rh,    &
                         smr, u_s, gsf, stcon, t0tile, lai_ft, so4_vd,         &
                         rc, o3_stom_frac )

USE asad_mod,                ONLY:                                             &
    ndepd,                                                                     &
    nldepd,                                                                    &
    speci,                                                                     &
    jpdd

USE ereport_mod,             ONLY:                                             &
    ereport

USE errormessagelength_mod,  ONLY:                                             &
    errormessagelength

USE ukca_ddepo3_ocean_mod,   ONLY:                                             &
    ukca_ddepo3_ocean

USE ukca_environment_fields_mod, ONLY:                                         &
    surf_wetness

USE parkind1,                ONLY:                                             &
    jprb,                                                                      &
    jpim

USE ukca_um_legacy_mod,      ONLY:                                             &
    rmol

USE ukca_config_specification_mod, ONLY:                                       &
    ukca_config

USE umPrintMgr,              ONLY:                                             &
    PrintStatus,                                                               &
    PrStatus_Diag,                                                             &
    umPrint,                                                                   &
    umMessage

USE yomhook,                 ONLY:                                             &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

!     Null resistance for deposition (1/r_null ~ 0)

REAL, PARAMETER :: r_null = 1.0e50

INTEGER, INTENT(IN) :: row_length         ! number columns
INTEGER, INTENT(IN) :: rows               ! number of rows
INTEGER, INTENT(IN) :: ntype              ! number of surface types
INTEGER, INTENT(IN) :: npft               ! number of plant functional types
REAL, INTENT(IN) :: sinlat(row_length,rows)
                ! Sine(latitude)
REAL, INTENT(IN) :: t0(row_length,rows)
                ! Surface temperature (K)
REAL, INTENT(IN) :: p0(row_length,rows)
                ! Surface pressure (Pa)
REAL, INTENT(IN) :: rh(row_length,rows)
                ! Relative humidity (fraction)
REAL, INTENT(IN) :: smr(row_length,rows)
                ! Soil moisture content (Fraction by volume)
REAL, INTENT(IN) :: u_s(row_length,rows)
                ! Surface friction velocity (m s-1)
REAL, INTENT(IN) :: gsf(row_length,rows,ntype)
                ! Global surface fractions
REAL, INTENT(IN) :: t0tile(row_length,rows,ntype)
                ! Surface temperature on tiles (K)

REAL, INTENT(IN) ::  stcon(row_length,rows,npft)
                ! Stomatal conductance (m s-1)
REAL, INTENT(IN) :: lai_ft(row_length,rows,npft)
                ! Leaf area index (m2 leaf m-2)

REAL, INTENT(IN) :: so4_vd(row_length,rows)
                ! Aerosol deposition velocity (m s-1)
                !(assumed to be the same for SO4 and other aerosols)

!     Surface resistance on tiles (s m-1).

REAL, INTENT(OUT) :: rc(row_length,rows,ntype,jpdd)
REAL, INTENT(OUT) :: o3_stom_frac(row_length,rows)

!     Local variables
INTEGER :: errcode                   ! error code
LOGICAL, SAVE :: first = .TRUE.
CHARACTER(LEN=errormessagelength) :: cmessage

!     Surface type indices (local copies)
INTEGER, PARAMETER :: n_elev_ice = 10
INTEGER, SAVE :: brd_leaf
INTEGER, SAVE :: brd_leaf_dec
INTEGER, SAVE :: brd_leaf_eg_trop
INTEGER, SAVE :: brd_leaf_eg_temp
INTEGER, SAVE :: ndl_leaf
INTEGER, SAVE :: ndl_leaf_dec
INTEGER, SAVE :: ndl_leaf_eg
INTEGER, SAVE :: c3_grass
INTEGER, SAVE :: c3_crop
INTEGER, SAVE :: c3_pasture
INTEGER, SAVE :: c4_grass
INTEGER, SAVE :: c4_crop
INTEGER, SAVE :: c4_pasture
INTEGER, SAVE :: shrub
INTEGER, SAVE :: shrub_dec
INTEGER, SAVE :: shrub_eg
INTEGER, SAVE :: urban
INTEGER, SAVE :: lake
INTEGER, SAVE :: soil
INTEGER, SAVE :: ice
INTEGER, SAVE :: elev_ice(n_elev_ice)

INTEGER :: i         ! Loop count over longitudes
INTEGER :: k         ! Loop count over latitudes
INTEGER :: j         ! Loop count over species that deposit
INTEGER :: n         ! Loop count over tiles

INTEGER, PARAMETER :: first_term  = 1
INTEGER, PARAMETER :: second_term = 2
INTEGER, PARAMETER :: third_term  = 3
INTEGER, PARAMETER :: fourth_term = 4

INTEGER, PARAMETER :: r_stom_spec = 5 ! NO2 , O3 , PAN , SO2, NH3
INTEGER, PARAMETER :: r_stom_no2  = 1
INTEGER, PARAMETER :: r_stom_o3   = 2
INTEGER, PARAMETER :: r_stom_pan  = 3
INTEGER, PARAMETER :: r_stom_so2  = 4
INTEGER, PARAMETER :: r_stom_nh3  = 5

! Tolerance (below which the number is assumed to be zero)
REAL, PARAMETER :: tol               = 1.0e-10
! Thresholds for determining if soil is wet or not (0.3 is a fairly
! arbitrary choice)
REAL, PARAMETER :: soil_moist_thresh = 0.3
! Relative humidity threshold
REAL, PARAMETER :: rh_thresh         = 0.813
! -1 degrees C in Kelvin
REAL, PARAMETER :: minus1degc        = 272.15
! -5 degrees C in Kelvin
REAL, PARAMETER :: minus5degc        = 268.15
! Wet soil and cuticular surface resistance for SO2 (s m-1)
REAL, PARAMETER :: r_wet_so2         = 1.0
! Dry cuticular surface resistance for SO2 (s m-1)
REAL, PARAMETER :: r_dry_so2         = 2000.0
! Wet cuticular surface resistance for SO2 when < -5 degrees C
REAL, PARAMETER :: r_cut_so2_5degc   = 500.0
! Wet cuticular surface resistance for SO2 when > -5 and < -1 degrees C
REAL, PARAMETER :: r_cut_so2_5to_1   = 200.0
! The constants used in Equation 9 of Erisman, Pul and Wyers (1994).
REAL, PARAMETER :: so2con1           = 2.5e4
REAL, PARAMETER :: so2con2           = 0.58e12
REAL, PARAMETER :: so2con3           = -6.93
REAL, PARAMETER :: so2con4           = -27.8
! The constants used in Reay et al. (2001) for CH4 soil uptake efficiency.
REAL, PARAMETER :: smfrac1           = 0.16
REAL, PARAMETER :: smfrac2           = 0.30
REAL, PARAMETER :: smfrac3           = 0.50
REAL, PARAMETER :: smfrac4           = 0.20


REAL :: sm            ! Soil moisture content of gridbox
REAL :: rr            ! General temporary store for resistances.
REAL :: ts            ! Temperature of a particular tile.
REAL :: f_ch4_uptk    ! Factor to modify CH4 uptake fluxes

REAL :: t_surf1       ! surface temperature
REAL :: p_surf1       ! surface pressure
REAL :: ustar         ! Friction velocity
REAL :: sst           ! Sea surface temperature (K)
REAL :: rc_ocean      ! ocean surface resistance term (s/m)

REAL :: mml = 1.008e5 ! Factor to convert methane flux to dry dep vel.
REAL :: r_wet_o3 = 500.0  ! Wet soil surface resistance for ozone (s m-1)
REAL :: cuticle_o3 = 5000.0 ! Constant for caln of O3 cuticular resistance
REAL :: tundra_s_limit = 0.866 ! Southern limit for tundra (SIN(60))

!     MML - Used to convert methane flux in ug m-2 h-1 to dry dep vel in
!     m s-1; MML=3600*0.016*1.0E9*1.75E-6, where 0.016=RMM methane (kg),
!     1.0E9 converts ug -> kg, 1.75E-6 = assumed CH4 vmr

REAL :: r_cut_o3(row_length,npft) ! Cuticular Resistance for O3

! The cuticular resistance for SO2
REAL :: r_cut_so2

REAL :: r_stom(row_length,rows,npft,r_stom_spec) ! Stomatal resistance
REAL, ALLOCATABLE, SAVE :: rsurf(:,:)        ! Standard surface resistances

!     Scaling of CH4 soil uptake to match present-day TAR value of 30 Tg/year

REAL, PARAMETER :: TAR_scaling = 15.0
REAL, PARAMETER :: glmin       = 1.0e-6  ! Minimum leaf conductance

!     Following arrays used to set up surface resistance array rsurf.
!      Must take the same dimension as ntype
REAL, ALLOCATABLE, SAVE :: r_null_tile(:)
REAL, ALLOCATABLE, SAVE :: rooh(:)
REAL, ALLOCATABLE, SAVE :: aerosol(:)
REAL, ALLOCATABLE, SAVE :: tenpointzero(:)

!     CH4 uptake fluxes, in ug m-2 hr-1
REAL, ALLOCATABLE, SAVE :: ch4_up_flux(:)


!     Hydrogen - linear dependence on soil moisture (except savannah)
!      Must take the same dimension as npft
REAL, ALLOCATABLE, SAVE :: h2dd_c(:)
REAL, ALLOCATABLE, SAVE :: h2dd_m(:)
REAL, ALLOCATABLE, SAVE :: h2dd_q(:) ! Quadratic term ( savannah only )

!     Resistances for Tundra if different to standard value in rsurf().

INTEGER, PARAMETER :: r_tundra_spec = 5 ! NO2 , CO , O3 , H2, PAN
INTEGER, PARAMETER :: r_tundra_no2  = 1
INTEGER, PARAMETER :: r_tundra_co   = 2
INTEGER, PARAMETER :: r_tundra_o3   = 3
INTEGER, PARAMETER :: r_tundra_pan  = 5

REAL :: r_tundra(r_tundra_spec) = [ 1200.0, 25000.0, 800.0, 3850.0, 1100.0]
!                                    NO2     CO       O3     H2      PAN

!     CH4 loss to tundra - Cubic polynomial fit to data.
!     N.B. Loss flux is in units of ug(CH4) m-2 s-1
REAL :: ch4dd_tun(4) = [ -4.757e-6, 4.0288e-3, -1.13592, 106.636 ]


!     HNO3 dry dep to ice; quadratic dependence
REAL :: hno3dd_ice(3) = [ -13.57, 6841.9, -857410.6 ]

!   Diffusion correction for stomatal conductance Wesley (1989)
!                                     Atmos. Env. 23, 1293.
REAL :: dif(5) = [ 1.6, 1.6, 2.6, 1.9, 0.97 ]
!                   NO2   O3  PAN  SO2  NH3

REAL, PARAMETER :: min_tile_frac = 0.0 ! minimum of tile fraction to
!                  include tile in CH4 soil uptake calculations
!                  (filter on todo-list)
LOGICAL :: todo(row_length,rows,ntype) ! True if tile fraction > min_tile_frac
LOGICAL :: microb(row_length,rows)     ! True if T > 5 C and RH > 40%
!                  (i.e. microbes in soil are active).

! Number of species for which surface resistance values are not set
INTEGER  :: n_nosurf

! Temporary logicals (local copies)
LOGICAL, SAVE :: l_fix_improve_drydep   ! True to fix dry deposition velocities
LOGICAL, SAVE :: l_fix_ukca_h2dd_x      ! True to fix H2 depos. to shrub/soil
LOGICAL, SAVE :: l_fix_drydep_so2_water ! True to use correct surface resistance
                                        ! of water when calculating dry
                                        ! deposition of SO2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SURFDDR'


! Set up standard resistance array rsurf on first call only
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first) THEN

  l_fix_improve_drydep = ukca_config%l_fix_improve_drydep
  l_fix_ukca_h2dd_x = ukca_config%l_fix_ukca_h2dd_x
  l_fix_drydep_so2_water = ukca_config%l_fix_drydep_so2_water

  brd_leaf         = ukca_config%i_brd_leaf
  brd_leaf_dec     = ukca_config%i_brd_leaf_dec
  brd_leaf_eg_trop = ukca_config%i_brd_leaf_eg_trop
  brd_leaf_eg_temp = ukca_config%i_brd_leaf_eg_temp
  ndl_leaf         = ukca_config%i_ndl_leaf
  ndl_leaf_dec     = ukca_config%i_ndl_leaf_dec
  ndl_leaf_eg      = ukca_config%i_ndl_leaf_eg
  c3_grass         = ukca_config%i_c3_grass
  c3_crop          = ukca_config%i_c3_crop
  c3_pasture       = ukca_config%i_c3_pasture
  c4_grass         = ukca_config%i_c4_grass
  c4_crop          = ukca_config%i_c4_crop
  c4_pasture       = ukca_config%i_c4_pasture
  shrub            = ukca_config%i_shrub
  shrub_dec        = ukca_config%i_shrub_dec
  shrub_eg         = ukca_config%i_shrub_eg
  urban            = ukca_config%i_urban
  lake             = ukca_config%i_lake
  soil             = ukca_config%i_soil
  ice              = ukca_config%i_ice

  IF (ntype == 27) THEN
    IF (SIZE(ukca_config%i_elev_ice) /= n_elev_ice) THEN
      WRITE(cmessage,'(A,I0,A,I0,A)') 'No of elevated ice surface types (',    &
        SIZE(ukca_config%i_elev_ice), ') does not match expected (',           &
        n_elev_ice, ')'
      errcode=1
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
    elev_ice = ukca_config%i_elev_ice
  END IF

  ALLOCATE(rsurf(ntype,jpdd))
  ALLOCATE(r_null_tile(ntype))
  ALLOCATE(rooh(ntype))
  ALLOCATE(aerosol(ntype))
  ALLOCATE(tenpointzero(ntype))
  ALLOCATE(ch4_up_flux(ntype))
  ALLOCATE(h2dd_c(ntype))
  ALLOCATE(h2dd_m(ntype))
  ALLOCATE(h2dd_q(ntype))
  rsurf(:,:) = r_null
  r_null_tile(:)    = r_null
  tenpointzero(:) = 10.0

  ! Check if we have standard JULES tile configuration and fail if not as
  ! Dry deposition in UKCA has not been set up for this.
  ! Use standard tile set up or extend UKCA dry depostion
  WRITE(cmessage,'(A)')                                                        &
       'UKCA does not handle flexible tiles yet: '//                           &
       'Dry deposition needs extending. '         //                           &
       'Please use standard tile configuration'
  SELECT CASE (ntype)
  CASE (9)
    IF ( brd_leaf /= 1  .OR.                                                   &
         ndl_leaf /= 2  .OR.                                                   &
         c3_grass /= 3  .OR.                                                   &
         c4_grass /= 4  .OR.                                                   &
         shrub    /= 5  .OR.                                                   &
         urban    /= 6  .OR.                                                   &
         lake     /= 7  .OR.                                                   &
         soil     /= 8  .OR.                                                   &
         ice      /= 9 ) THEN

      ! Tile order changed from standard setup i.e. MOSES-II.
      WRITE(umMessage,'(A)') 'Nine tile order changed from standard setup'
      CALL umPrint(umMessage,src=RoutineName)
      errcode=1009
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  CASE (13,17,27)
    IF ( brd_leaf_dec     /= 1   .OR.                                          &
         brd_leaf_eg_trop /= 2   .OR.                                          &
         brd_leaf_eg_temp /= 3   .OR.                                          &
         ndl_leaf_dec     /= 4   .OR.                                          &
         ndl_leaf_eg      /= 5   .OR.                                          &
         c3_grass         /= 6 ) THEN

      ! Tile order does not match that expected
      WRITE(umMessage,'(A)') 'Tile order changed from standard setup'
      CALL umPrint(umMessage,src=RoutineName)
      errcode=1001
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  CASE DEFAULT
    ! ntype must equal 9, 13, 17, 27
    WRITE(umMessage,'(A,I0)') 'ntype must be set to 9, 13, 17, 27.' //         &
                              'ntype = ', ntype
    CALL umPrint(umMessage,src=RoutineName)
    errcode=1002
    CALL ereport(RoutineName,errcode,cmessage)
  END SELECT

  SELECT CASE (ntype)
  CASE (13)
    IF ( c4_grass         /= 7    .OR.                                         &
         shrub_dec        /= 8    .OR.                                         &
         shrub_eg         /= 9    .OR.                                         &
         urban            /= 10   .OR.                                         &
         lake             /= 11   .OR.                                         &
         soil             /= 12   .OR.                                         &
         ice              /= 13 ) THEN

      ! Tile order does not match that expected
      WRITE(umMessage,'(A)') 'Thirteen tile order changed from standard setup'
      CALL umPrint(umMessage,src=RoutineName)
      errcode=1013
      CALL ereport(RoutineName,errcode,cmessage)
    END IF

  CASE (17, 27)
    IF ( c3_crop          /= 7    .OR.                                         &
         c3_pasture       /= 8    .OR.                                         &
         c4_grass         /= 9    .OR.                                         &
         c4_crop          /= 10   .OR.                                         &
         c4_pasture       /= 11   .OR.                                         &
         shrub_dec        /= 12   .OR.                                         &
         shrub_eg         /= 13   .OR.                                         &
         urban            /= 14   .OR.                                         &
         lake             /= 15   .OR.                                         &
         soil             /= 16   .OR.                                         &
         ice              /= 17 ) THEN

      ! Tile order does not match that expected
      WRITE(umMessage,'(A)') 'Tile order changed from standard setup'
      CALL umPrint(umMessage,src=RoutineName)
      errcode=1727
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  END SELECT

  IF ( ntype == 27) THEN
    IF ( elev_ice(1)      /= 18   .OR.                                         &
         elev_ice(2)      /= 19   .OR.                                         &
         elev_ice(3)      /= 20   .OR.                                         &
         elev_ice(4)      /= 21   .OR.                                         &
         elev_ice(5)      /= 22   .OR.                                         &
         elev_ice(6)      /= 23   .OR.                                         &
         elev_ice(7)      /= 24   .OR.                                         &
         elev_ice(8)      /= 25   .OR.                                         &
         elev_ice(9)      /= 26   .OR.                                         &
         elev_ice(10)     /= 27 ) THEN

      WRITE(umMessage,'(A)') 'Tile order changed from standard setup'
      CALL umPrint(umMessage,src=RoutineName)
      ! Tile order does not match that expected
      errcode=1027
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  END IF

  SELECT CASE (ntype)
  CASE (9)
    IF (l_fix_improve_drydep) THEN
      rooh =        [279.2,238.2,366.3,322.9,362.5,424.9,933.1,585.4,1156.1]
    ELSE
      rooh =        [30.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0]
    END IF

    aerosol =     [r_null,r_null,r_null,r_null,r_null,r_null,                  &
                    1000.0,r_null,20000.0]

    IF (l_fix_improve_drydep) THEN
      ! soil should behave as c3_grass not as shrub
      ch4_up_flux = [39.5,50.0,30.0,37.0,27.5,0.0,0.0,30.0,0.0]
    ELSE
      ch4_up_flux = [39.5,50.0,30.0,37.0,27.5,0.0,0.0,27.5,0.0]
    END IF

    IF ( l_fix_ukca_h2dd_x ) THEN
      ! from Table 1 Sanderson, J. Atmos. Chem., v46, 15-28, 2003
      ! note : conversion factor 1.0e-4 (from paper) will only be applied
      !        in equation. Prior to the temporary logical switch, this was
      !        being applied twice in a few places.

      ! urban, lake and ice not used
      ! shrub/soil same as c3_grass below tundra limit else same as shrub
      h2dd_c =      [  19.70,  19.70,  17.70,  1.235,  1.000,                  &
                        0.000,  0.000,  0.000,  0.000 ]

      ! urban, lake and ice not used
      ! shrub/soil same as c3_grass below tundra limit else same as shrub
      h2dd_m =      [ -41.90, -41.90, -41.39, -0.472,  0.000,                  &
                        0.000,  0.000, -0.000,  0.000 ]

      ! Quadratic term for H2 loss to savannah only
      h2dd_q =      [   0.00,   0.00,   0.00,   0.27,   0.00,                  &
                         0.00,   0.00,   0.00,   0.00 ]
    ELSE
      ! from Table 1 Sanderson, J. Atmos. Chem., v46, 15-28, 2003
      !  with some errors

      ! urban, lake and ice not used
      ! shrub/soil same as c3_grass below tundra limit else same as shrub
      ! incorrect value for c4_grass
      h2dd_c =      [  0.00197,  0.00197,  0.00177,  1.23460, 0.00010,         &
                        0.00000,  0.00000,  0.00000,  0.00000 ]

      ! urban, lake and ice not used
      ! shrub/soil same as c3_grass below tundra limit else same as shrub
      h2dd_m =      [ -0.00419, -0.00419, -0.00414, -0.47200, 0.00000,         &
                        0.00000,  0.00000,  0.00000,  0.00000 ]

      ! Quadratic term for H2 loss (incorrect as should be savannah only)
      h2dd_q =      [     0.27,     0.27,     0.27,     0.27,    0.27,         &
                           0.27,     0.27,     0.27,     0.27 ]
    END IF

  CASE (13,17,27)
    rooh(1:6)       = [  300.3,  270.3,  266.9,  238.0,  238.5,  366.3 ]
    aerosol(1:6)    = [ r_null, r_null, r_null, r_null, r_null, r_null ]
    ch4_up_flux(1:6)= [   39.5,   39.5,   39.5,   50.0,   50.0,   30.0 ]
    h2dd_c(1:6)     = [  19.70,  19.70,  19.70,  19.70,  19.70,  17.70 ]
    h2dd_m(1:6)     = [ -41.90, -41.90, -41.90, -41.90, -41.90, -41.39 ]
    IF ( l_fix_ukca_h2dd_x ) THEN
      h2dd_q(1:6)   = [   0.00,   0.00,   0.00,   0.00,   0.00,   0.00 ]
    ELSE
      h2dd_q(1:6)   = [   0.27,   0.27,   0.27,   0.27,   0.27,   0.27 ]
    END IF
  CASE DEFAULT
    WRITE(umMessage,'(A,I0)') 'ntype = ', ntype
    CALL umPrint(umMessage,src=RoutineName)
    errcode=1
    WRITE(cmessage,'(A)')                                                      &
         'UKCA does not handle flexible tiles yet: ' //                        &
         'Dry deposition needs extending. '          //                        &
         'Please use standard tile configuration'
    CALL ereport(RoutineName,errcode,cmessage)
  END SELECT

  SELECT CASE (ntype)
  CASE (13)
    rooh(7:13)       = [   322.9, 332.8, 392.2, 424.9, 933.1, 585.4, 1156.1 ]
    aerosol(7:13)    = [  r_null,r_null,r_null,r_null,1000.0,r_null,20000.0 ]

    IF ( l_fix_ukca_h2dd_x ) THEN
      !                   c4_grass,shrub_dec,shrub_eg,urban, lake, soil, ice
      ch4_up_flux(7:13)= [   37.0,     27.5,    27.5,  0.0,  0.0, 30.0, 0.0 ]
      h2dd_c(7:13)     = [  1.235,      1.0,     1.0,  0.0,  0.0,  0.0, 0.0 ]
      ! shrub/soil same as c3_grass below tundra limit else same as shrub
      h2dd_m(7:13)     = [ -0.472,      0.0,     0.0,  0.0,  0.0,  0.0, 0.0 ]
      h2dd_q(7:13)     = [   0.27,      0.0,     0.0,  0.0,  0.0,  0.0, 0.0 ]
    ELSE
      !                   c4_grass,shrub_dec,shrub_eg,urban, lake, soil, ice
      ch4_up_flux(7:13)= [   37.0,     27.5,    27.5,  0.0,  0.0, 27.5, 0.0 ]
      h2dd_c(7:13)     = [ 1.2346,   0.0001,  0.0001,  0.0,  0.0,  0.0, 0.0 ]
      ! shrub/soil same as c3_grass below tundra limit else same as shrub
      h2dd_m(7:13)     = [ -0.472,      0.0,     0.0,  0.0,  0.0,  0.0, 0.0 ]
      h2dd_q(7:13)     = [   0.27,     0.27,    0.27, 0.27, 0.27, 0.27,0.27 ]
    END IF
  CASE (17, 27)
    rooh(7:17)       = [  366.3,  366.3,  322.9,  322.9,  322.9,               &
                           332.8,  392.2,  424.9,  933.1,  585.4,  1156.1 ]
    aerosol(7:17)    = [ r_null, r_null, r_null, r_null, r_null,               &
                          r_null, r_null, r_null, 1000.0, r_null, 20000.0]

    IF ( l_fix_ukca_h2dd_x ) THEN
      ch4_up_flux(7:17)= [30.0,     30.0,     37.0,     37.0,     37.0,        &
                           27.5,     27.5,      0.0,      0.0,     30.0, 0.0 ]
      ! from Table 1 Sanderson, J. Atmos. Chem., v46, 15-28, 2003
      ! urban, lake and ice not used & soil same as c3_grass
      h2dd_c(7:17) = [   17.70,    17.70,    1.235,    1.235,    1.235,        &
                          1.000,    1.000,    0.000,    0.000,    0.000, 0.0 ]

      ! urban, lake and ice not used
      ! shrub/soil same as c3_grass below tundra limit else same as shrub
      h2dd_m(7:17) = [  -41.39,   -41.39,   -0.472,   -0.472,   -0.472,        &
                          0.000,    0.000,    0.000,    0.000,   -41.39, 0.0 ]

      ! c4_grass, c4_crop, c4_pasture contain quadratic terms for H2 loss
      h2dd_q(7:17) = [    0.00,     0.00,     0.27,     0.27,     0.27,        &
                           0.00,     0.00,     0.00,     0.00,     0.00, 0.0 ]
    ELSE
      ch4_up_flux(7:17)= [30.0,     30.0,     37.0,     37.0,     37.0,        &
                           27.5,     27.5,      0.0,      0.0,     27.5, 0.0 ]

      h2dd_c(7:17) = [ 0.00177,  0.00177,  1.23460,  1.23460,  1.23460,        &
                        0.00010,  0.00010,  0.00000,  0.00000,  0.00000, 0.0 ]

      h2dd_m(7:17) = [-0.00414, -0.00414, -0.47200, -0.47200, -0.47200,        &
                        0.00000,  0.00000,  0.00000,  0.00000,  0.00000, 0.0 ]

      h2dd_q(7:17) = [    0.27,     0.27,     0.27,     0.27,     0.27,        &
                           0.27,     0.27,     0.27,     0.27,     0.27, 0.27 ]
    END IF
  END SELECT

  IF (ntype == 27) THEN
    rooh(18:27)        = [  1156.1, 1156.1, 1156.1, 1156.1, 1156.1,            &
                             1156.1, 1156.1, 1156.1, 1156.1, 1156.1 ]
    aerosol(18:27)     = [ 20000.0,20000.0,20000.0,20000.0,20000.0,            &
                            20000.0,20000.0,20000.0,20000.0,20000.0 ]
    ch4_up_flux(18:27) = [     0.0,    0.0,    0.0,    0.0,    0.0,            &
                                0.0,    0.0,    0.0,    0.0,    0.0 ]
    h2dd_c(18:27)      = [     0.0,    0.0,    0.0,    0.0,    0.0,            &
                                0.0,    0.0,    0.0,    0.0,    0.0 ]
    h2dd_m(18:27)      = [     0.0,    0.0,    0.0,    0.0,    0.0,            &
                                0.0,    0.0,    0.0,    0.0,    0.0 ]
    IF ( l_fix_ukca_h2dd_x ) THEN
      h2dd_q(18:27)      = [   0.0,    0.0,    0.0,    0.0,    0.0,            &
                                0.0,    0.0,    0.0,    0.0,    0.0 ]
    ELSE
      h2dd_q(18:27)      = [  0.27,   0.27,   0.27,   0.27,   0.27,            &
                               0.27,   0.27,   0.27,   0.27,   0.27 ]
    END IF
  END IF

  SELECT CASE (ntype)
  CASE (9)
    ! Standard surface resistances (s m-1). Values are for 9 tiles in
    ! order: Broadleaved trees, Needleleaf trees, C3 Grass, C4 Grass,
    ! Shrub, Urban, Water, Bare Soil, Ice.
    !
    ! Nine tile dry deposition surface resistance have been updated April 2018
    ! (behind l_fix_improve_drydep logical) to reflect 13/17/27 tiles
    n_nosurf = 0
    DO n = 1, ndepd
      SELECT CASE (speci(nldepd(n)))
      CASE ('O3        ','O3S       ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 219.3, 233.0, 355.0, 309.3, 358.2,                      &
                        444.4,2000.0, 645.2,2000.0 ]
        ELSE
          rsurf(:,n)=[ 200.0, 200.0, 200.0, 200.0, 400.0,                      &
                        800.0,2200.0, 800.0,2500.0 ]
        END IF

      CASE ('NO2       ','NO3       ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 364.1, 291.3, 443.8, 386.6, 447.8,                      &
                        555.6,2500.0, 806.5,2500.0 ]
        ELSE
          rsurf(:,n)=[ 225.0, 225.0, 400.0, 400.0, 600.0,                      &
                       1200.0,2600.0,1200.0,3500.0 ]
        END IF
      CASE ('NO        ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 2184.5, 1747.6, 3662.7, 2319.6, 2686.8,                 &
                        3333.3,15000.0, 4838.7,15000.0 ]
        ELSE
          rsurf(:,n)=[ 1350.0, 1350.0, 2400.0, 2400.0, 3600.0,                 &
                       72000.0, r_null,72000.0,21000.0 ]
        END IF
      CASE ('HNO3      ','HONO2     ','B2ndry    ','A2ndry    ','N2O5      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 8.5, 8.4, 13.2, 12.0, 62.3, 12.8, 13.9, 16.0, 19.4 ]
        ELSE
          rsurf(:,n)=tenpointzero
        END IF
      CASE ('HNO4      ','HO2NO2    ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 17.0, 16.8, 26.4, 24.0, 24.9, 25.7, 27.7, 32.1, 38.8 ]
        ELSE
          rsurf(:,n)=tenpointzero
        END IF
        ! CRI organic nitrates deposit as ISON
      CASE ('ISON      ',                                                      &
            'HOC2H4NO3 ','RTX24NO3  ','RN9NO3    ',                            &
            'RN12NO3   ','RN15NO3   ','RN18NO3   ','RU14NO3   ',               &
            'RTN28NO3  ','RTN25NO3  ','RTX28NO3  ','RTX22NO3  ',               &
            'RA22NO3   ','RA25NO3   ','RTN23NO3  ','RU12NO3   ',               &
            'RU10NO3   ','RA13NO3   ','RA16NO3   ','RA19NO3   ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 582.5, 466.0, 710.1, 618.6, 716.5,                      &
                        888.9,4000.0,1290.3,4000.0 ]
        ELSE
          rsurf(:,n)=tenpointzero
        END IF
      CASE ('HCl       ','HOCl      ','HBr       ','HOBr      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 8.5, 8.4, 13.2, 12.0, 62.3, 12.8, 13.9, 16.0, 19.4 ]
        ELSE
          rsurf(:,n)=r_null_tile(:)
        END IF
      CASE ('H2SO4     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)= [  84.9,  83.8, 131.9, 120.0, 124.4,                     &
                         128.5, 138.6, 160.4, 194.2 ]
        ELSE
          rsurf(:,n)=r_null_tile(:)
        END IF
      CASE ('H2O2      ','HOOH      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 84.9,83.8,131.9,120.0,124.4,128.5,138.6,160.4,194.2 ]
        ELSE
          rsurf(:,n)=tenpointzero
        END IF
        ! CRI R-OOH species deposit as MeOOH
        ! Also including hydroxy-ketones (CARB7-16) which deposit as HACET
      CASE ('CH3OOH    ','MeOOH     ','C2H5OOH   ','EtOOH     ',               &
            'n_C3H7OOH ','i_C3H7OOH ','n-PrOOH   ','i-PrOOH   ',               &
            'MeCOCH2OOH','ISOOH     ','MACROOH   ','MeCO3H    ',               &
            'MeCO2H    ','HCOOH     ','PropeOOH  ','MEKOOH    ',               &
            'ALKAOOH   ','AROMOOH   ','BSVOC1    ','BSVOC2    ',               &
            'ASVOC1    ','ASVOC2    ','ISOSVOC1  ','ISOSVOC2  ',               &
            's-BuOOH   ','MVKOOH    ','HACET     ',                            &
            'EtCO3H    ','HOCH2CO3H ','RCOOH25   ',                            &
            'HOC2H4OOH ','RN10OOH   ','RN13OOH   ','RN16OOH   ',               &
            'RN19OOH   ','RN8OOH    ','RN11OOH   ','RN14OOH   ',               &
            'RN17OOH   ','RU14OOH   ','RU12OOH   ','RU10OOH   ',               &
            'NRU14OOH  ','NRU12OOH  ','RN9OOH    ','RN12OOH   ',               &
            'RN15OOH   ','RN18OOH   ','NRN6OOH   ','NRN9OOH   ',               &
            'NRN12OOH  ','RA13OOH   ','RA16OOH   ','RA19OOH   ',               &
            'RTN28OOH  ','NRTN28OOH ','RTN26OOH  ','RTN25OOH  ',               &
            'RTN24OOH  ','RTN23OOH  ','RTN14OOH  ','RTN10OOH  ',               &
            'RTX28OOH  ','RTX24OOH  ','RTX22OOH  ','NRTX28OOH ',               &
            'CARB7     ','CARB10    ','CARB13    ','CARB16    ',               &
            'RA22OOH   ','RA25OOH   ','HPUCARB12 ','DHPCARB9  ',               &
            'DHPR12OOH')
        rsurf(:,n)=rooh
        ! CRI PAN-type species
      CASE ('PAN       ','PPAN      ',                                         &
            'PHAN      ','RU12PAN   ','RTN26PAN  ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 485.4, 388.4, 591.7, 515.5, 597.1,                      &
                        740.7,3333.3,1075.3,3333.3 ]
        ELSE
          rsurf(:,n)=[ 500.0,  500.0, 500.0,  500.0, 500.0,                    &
                       r_null,12500.0, 500.0,12500.0 ]
        END IF
      CASE ('MPAN      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[  970.9,  776.7, 1183.4, 1030.9, 1194.2,                 &
                        1481.5, 6666.7, 2150.5, 6666.7 ]
        ELSE
          rsurf(:,n)=[ 500.0,  500.0, 500.0,  500.0, 500.0,                    &
                       r_null,12500.0, 500.0,12500.0 ]
        END IF
      CASE ('OnitU     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 582.5,  466.0, 710.1,  618.6, 716.5,                    &
                        888.9, 4000.0,1290.3, 4000.0 ]
        ELSE
          rsurf(:,n)=[ 500.0,  500.0, 500.0,  500.0, 500.0,                    &
                       r_null,12500.0, 500.0,12500.0 ]
        END IF
      CASE ('NH3       ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 120.0, 130.9, 209.8, 196.1, 191.0,                      &
                        180.7, 148.9, 213.5, 215.1 ]
        ELSE
          rsurf(:,n)=tenpointzero
        END IF
      CASE ('CO        ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 3700.0, 7300.0, 4550.0, 1960.0, 4550.0,                 &
                        r_null, r_null, 4550.0, r_null ]
        ELSE
          rsurf(:,n)=[ 3700.0, 7300.0, 4550.0, 1960.0, 4550.0,                 &
                        r_null, r_null, 4550.0, r_null ]
        END IF
        ! Shrub+bare soil set to C3 grass (guess)
      CASE ('CH4       ')
        IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
          rsurf(:,n)=1.0/r_null_tile(:)
        ELSE
          rsurf(:,n)=r_null_tile(:)
        END IF
      CASE ('HONO      ')
        rsurf(:,n)=r_null_tile(:)
      CASE ('H2        ')
        rsurf(:,n)=r_null_tile(:)
      CASE ('SO2       ')
        IF (ukca_config%l_ukca_dry_dep_so2wet) THEN
          ! With the implementation of the Erisman, Pul and Wyers (1994)
          ! parameterization we additionally reduce the surface resistance
          ! value of water to 1.0 s m-1.
          rsurf(:,n)=[ 120.0, 130.9, 209.8, 196.1, 191.0,                      &
                        180.7,   1.0, 213.5, 215.1 ]
        ELSE
          IF (l_fix_improve_drydep) THEN
            rsurf(:,n)=[ 120.0, 130.9, 209.8, 196.1, 191.0,                    &
                          180.7,  10.0, 213.5, 215.1 ]
          ELSE
            rsurf(:,n)=[ 100.0, 100.0, 150.0, 350.0, 400.0,                    &
                          400.0,  10.0, 700.0, r_null ]
          END IF
        END IF
      CASE ('BSOA      ','ASOA      ','ISOSOA    ')
        rsurf(:,n)=aerosol
      CASE ('ORGNIT    ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 582.5, 466.0, 710.1, 618.6, 716.5,                      &
                        888.9,4000.0,1290.3,4000.0 ]
        ELSE
          rsurf(:,n)=aerosol
        END IF
      CASE ('Sec_Org   ','SEC_ORG_I ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=aerosol
        ELSE
          rsurf(:,n)=r_null_tile(:)
        END IF
      CASE ('HCHO      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 136.0, 143.5, 228.5, 211.6, 210.5,                      &
                        205.1, 182.7, 246.5, 261.8 ]
        ELSE
          rsurf(:,n)=[ 100.0, 100.0, 150.0, 350.0, 600.0,                      &
                        400.0, 200.0, 700.0, 700.0 ]
        END IF
        ! CRI alcohols deposit as MeOH
      CASE ('MeOH      ','EtOH      ','i-PrOH    ','n-PrOH    ',               &
            'AROH14    ','ARNOH14   ','AROH17    ','ARNOH17   ',               &
            'IEPOX     ','HMML      ','HUCARB9   ','DHCARB9   ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 187.1, 199.4, 318.3, 295.6, 292.2,                      &
                        282.1, 245.1, 337.2, 352.2 ]
        ELSE
          rsurf(:,n)=r_null_tile(:)
        END IF
        ! CRI Carbonyls, copying MeCHO rates
        ! Second generation nitrates deposit as NALD
      CASE ('MeCHO     ','EtCHO     ','MACR      ','NALD      ',               &
            'HOCH2CHO  ',                                                      &
            'CARB14    ','CARB17    ','CARB11A   ','UCARB10   ',               &
            'CARB15    ','UCARB12   ','UDCARB8   ','UDCARB11  ',               &
            'UDCARB14  ','TNCARB26  ','TNCARB10  ','TNCARB12  ',               &
            'TNCARB11  ','CCARB12   ','TNCARB15  ','TXCARB24  ',               &
            'TXCARB22  ','UDCARB17  ','NOA       ','NUCARB12  ',               &
            'ANHY      '                                                       &
            )
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 5825.2, 4660.3, 7100.6, 6185.6, 7164.8,                 &
                        8888.9,40000.0,12903.2,40000.0 ]
        ELSE
          rsurf(:,n)=[ 1200.0, 1200.0, 1200.0, 1200.0, 1000.0,                 &
                        2400.0, r_null, r_null, r_null ]
        END IF
        ! CRI: CARB3 = GLY, CARB6 = MGLY etc.
      CASE ('MGLY      ','CARB3     ','CARB6     ','CARB9     ','CARB12    ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=[ 12001.2, 13086.3, 20979.0, 19607.8, 19091.9,            &
                        18072.3, 14888.3, 21352.3, 21505.4 ]
        ELSE
          rsurf(:,n)=[ 1200.0, 1200.0, 1200.0, 1200.0, 1000.0,                 &
                        2400.0, r_null, r_null, r_null ]
        END IF
      CASE ('DMSO      ','Monoterp  ','APINENE   ','BPINENE   ')
        rsurf(:,n)=r_null_tile(:)
      CASE DEFAULT
        n_nosurf = n_nosurf + 1
        IF (first .AND. PrintStatus == PrStatus_Diag ) THEN
          WRITE(umMessage,'(A)') 'Surface resistance values not set for '//    &
                                  speci(nldepd(n))
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      END SELECT
    END DO
    ! Fail if surface resistance values unavailable for some species
    IF (first .AND. n_nosurf > 0) THEN
      WRITE(cmessage,'(A,I0,A)')                                               &
        ' Surface resistance values not found for ',n_nosurf,' species.'
      errcode = ABS(n_nosurf)
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  CASE (13, 17, 27)
    ! Standard surface resistances (s m-1). Values are for 13/17/27 tiles in
    ! order: Broadleaf deciduous trees, Broadleaf evergreen tropical trees,
    !        Broadleaf evergreem temperate trees, Needleleaf deciduous trees,
    !        Needleleaf evergreen trees, C3 Grass,
    !
    n_nosurf = 0
    DO n = 1, ndepd
      SELECT CASE (speci(nldepd(n)))
      CASE ('O3        ','O3S       ')
        rsurf(1:6,n)=[ 307.7,285.7,280.4,232.6,233.5,355.0 ]
      CASE ('SO2       ')
        rsurf(1:6,n)=[ 137.0,111.1,111.9,131.3,130.4,209.8 ]
      CASE ('NO2       ','NO3       ')
        rsurf(1:6,n)=[ 384.6,357.1,350.5,290.7,291.8,443.8 ]
      CASE ('NO        ')
        rsurf(1:6,n)=[ 2307.7,2142.9,2102.8,1744.2,1751.0,3662.7 ]
      CASE ('HNO3      ','HONO2     ','B2ndry    ','A2ndry    ','N2O5      ')
        rsurf(1:6,n)=[ 9.5,8.0,8.0,8.4,8.4,13.2 ]
      CASE ('HCl       ','HOCl      ','HBr       ','HOBr      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(1:6,n)=[ 9.5,8.0,8.0,8.4,8.4,13.2 ]
        ELSE
          rsurf(:,n)=r_null_tile(:)
        END IF
      CASE ('H2SO4     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(1:6,n)=[ 94.8, 80.0, 80.0, 83.9, 83.7, 131.9 ]
        ELSE
          rsurf(:,n)=r_null_tile(:)
        END IF
      CASE ('HNO4      ','HO2NO2    ')
        rsurf(1:6,n)=[ 19.0,16.0,16.0,16.8,16.7,26.4 ]
      CASE ('HONO      ')
        rsurf(1:6,n)=[ 47.4,40.0,40.0,42.0,41.8,65.9 ]
      CASE ('H2O2      ','HOOH      ')
        rsurf(1:6,n)=[ 94.8,80.0,80.0,83.9,83.7,131.9 ]
        ! Extra CRI PAN species
      CASE ('PAN       ','PPAN      ',                                         &
            'PHAN      ','RU12PAN   ','RTN26PAN  ')
        rsurf(1:6,n)=[ 512.8,476.2,467.3,387.6,389.1,591.7 ]
      CASE ('MPAN      ')
        rsurf(1:6,n)=[ 1025.6,952.4,934.6,775.2,778.2,1183.4 ]
      CASE ('OnitU     ','ISON      ','ORGNIT    ',                            &
             'HOC2H4NO3 ','RN9NO3    ','RN12NO3   ','RTN25NO3  ',              &
             'RN15NO3   ','RN18NO3   ','RU14NO3   ','RTN28NO3  ',              &
             'RTX28NO3  ','RA22NO3   ','RA25NO3   ','RU12NO3   ',              &
             'RU10NO3   ','RTX22NO3  ','RTN23NO3  ','RTX24NO3  ',              &
             'RA13NO3   ','RA16NO3   ','RA19NO3   ')
        rsurf(1:6,n)=[ 615.4,571.4,560.7,465.1,466.9,710.1 ]
        ! CRI R-OOH species.
        ! Also including hydroxy-ketones (CARB7-16) which deposit as HACET
      CASE ('CH3OOH    ','MeOOH     ','C2H5OOH   ','EtOOH     ',               &
            'n_C3H7OOH ','i_C3H7OOH ','n-PrOOH   ','i-PrOOH   ',               &
            'MeCOCH2OOH','ISOOH     ','MACROOH   ','MeCO3H    ',               &
            'MeCO2H    ','HCOOH     ','PropeOOH  ','MEKOOH    ',               &
            'ALKAOOH   ','AROMOOH   ','BSVOC1    ','BSVOC2    ',               &
            'ASVOC1    ','ASVOC2    ','ISOSVOC1  ','ISOSVOC2  ',               &
            's-BuOOH   ','MVKOOH    ','HACET     ',                            &
            'EtCO3H    ','HOCH2CO3H ','RCOOH25   ',                            &
            'HOC2H4OOH ','RN10OOH   ','RN13OOH   ','RN16OOH   ',               &
            'RN19OOH   ','RN8OOH    ','RN11OOH   ','RN14OOH   ',               &
            'RN17OOH   ','RU14OOH   ','RU12OOH   ','RU10OOH   ',               &
            'NRU14OOH  ','NRU12OOH  ','RN9OOH    ','RN12OOH   ',               &
            'RN15OOH   ','RN18OOH   ','NRN6OOH   ','NRN9OOH   ',               &
            'NRN12OOH  ','RA13OOH   ','RA16OOH   ','RA19OOH   ',               &
            'RTN28OOH  ','NRTN28OOH ','RTN26OOH  ','RTN25OOH  ',               &
            'RTN24OOH  ','RTN23OOH  ','RTN14OOH  ','RTN10OOH  ',               &
            'RTX28OOH  ','RTX24OOH  ','RTX22OOH  ','NRTX28OOH ',               &
            'CARB7     ','CARB10    ','CARB13    ','CARB16    ',               &
            'RA22OOH   ','RA25OOH   ','HPUCARB12 ','DHPCARB9  ',               &
            'DHPR12OOH ')
        rsurf(:,n)=rooh
      CASE ('NH3       ')
        rsurf(1:6,n)=[ 137.0,111.1,111.9,131.3,130.4,209.8 ]
      CASE ('CO        ')
        rsurf(1:6,n)=[ 3700.0,3700.0,3700.0,7300.0,7300.0,4550.0 ]
        ! Shrub+bare soil set to C3 grass (guess)
      CASE ('CH4       ')
        IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
          rsurf(:,n)=1.0/r_null_tile(:) ! removal rate, NOT a resistance!
        ELSE
          rsurf(:,n)=r_null_tile(:)
        END IF
      CASE ('H2        ')
        rsurf(:,n)=r_null_tile(:)
      CASE ('HCHO      ')
        rsurf(1:6,n)=[ 154.1,126.6,127.2,143.8,143.1,228.5 ]
      CASE ('MeOH      ','EtOH      ','i-PrOH    ','n-PrOH    ',               &
            'AROH14    ','ARNOH14   ','AROH17    ','ARNOH17   ',               &
            'IEPOX     ','HMML      ','HUCARB9   ','DHCARB9   ')
        IF (l_fix_improve_drydep) THEN
          rsurf(1:6,n)=[ 212.6, 173.9, 174.9, 200.0, 198.0, 318.3 ]
        ELSE
          rsurf(:,n)=r_null_tile(:)
        END IF
        ! CRI Carbonyls, copying MeCHO rates
        ! Second generation nitrates deposit as NALD
      CASE ('MeCHO     ','EtCHO     ','MACR      ','NALD      ',               &
            'HOCH2CHO  ',                                                      &
            'CARB14    ','CARB17    ','CARB11A   ','UCARB10   ',               &
            'UCARB12   ','UDCARB8   ','UDCARB11  ',                            &
            'UDCARB14  ','TNCARB26  ','TNCARB10  ','TNCARB12  ',               &
            'TNCARB11  ','CCARB12   ','TNCARB15  ','TXCARB24  ',               &
            'TXCARB22  ','UDCARB17  ','NOA       ','NUCARB12  ',               &
            'ANHY      '                                                       &
            )
        rsurf(1:6,n)=[ 6153.8,5714.3,5607.5,4651.2,4669.3,7100.6 ]
        ! Glyoxal-type CRI species, deposit as MGLY
        !   n.b. CARB6 ~ MGLOX == MGLY
      CASE ('MGLY      ',                                                      &
            'CARB3     ','CARB6     ','CARB9     ','CARB12    ',               &
            'CARB15    '                                                       &
            )

        rsurf(1:6,n)=[ 13698.6,11111.1,11194.0,13129.1,13043.5,20979.0 ]
      CASE ('BSOA      ','ASOA      ','ISOSOA    ')
        rsurf(:,n)=aerosol
      CASE ('Sec_Org   ','SEC_ORG_I ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=aerosol
        ELSE
          rsurf(:,n)=r_null_tile(:)
        END IF
      CASE ('DMSO      ','Monoterp  ','APINENE   ','BPINENE   ')
        rsurf(:,n)=r_null_tile(:)
      CASE DEFAULT
        n_nosurf = n_nosurf + 1
        IF (first .AND. PrintStatus == PrStatus_Diag ) THEN
          umMessage = ' Surface resistance values not set for '//              &
                        speci(nldepd(n))
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      END SELECT
    END DO
    ! Fail if surface resistance values unavailable for some species
    IF (first .AND. n_nosurf > 0) THEN
      WRITE(cmessage,'(A,I0,A)')                                               &
      ' Surface resistance values not found for ',n_nosurf,' species.'
      errcode = ABS(n_nosurf)
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  CASE DEFAULT
    errcode=319
    WRITE(cmessage,'(A)')                                                      &
         'UKCA does not handle flexible tiles yet: '//                         &
         'Dry deposition needs extending. '         //                         &
         'Please use standard tile configuration'
    CALL ereport(RoutineName,errcode,cmessage)
  END SELECT

  SELECT CASE (ntype)
  CASE (13)
    ! Standard surface resistances (s m-1). Values are for 13 tiles in
    ! order: C4 Grass, Shrub deciduous, Shrub evergreen,
    !        Urban, Water, Bare Soil, Ice.
    DO n = 1, ndepd
      SELECT CASE (speci(nldepd(n)))
      CASE ('O3        ','O3S       ')
        rsurf(7:13,n)=[ 309.3,324.3,392.2,444.4,2000.0,645.2,2000.0 ]
      CASE ('SO2       ')
        IF (ukca_config%l_ukca_dry_dep_so2wet) THEN
          ! With the implementation of the Erisman, Pul and Wyers (1994)
          ! parameterization we additionally reduce the surface resistance
          ! value of water to 1.0 s m-1.
          rsurf(7:13,n)=[ 196.1,185.8,196.1,180.7,  1.0,213.5,215.1 ]
        ELSE
          IF (l_fix_drydep_so2_water) THEN
            ! Setting resistance of water to 148.9 s m-1 was an error, so
            ! this change reduces it to its previous value (10.0 s m-1).
            rsurf(7:13,n)=[ 196.1,185.8,196.1,180.7, 10.0,213.5,215.1 ]
          ELSE
            ! Retains previous suspect value of 148.9 s m-1 for water
            ! resistance.
            rsurf(7:13,n)=[ 196.1,185.8,196.1,180.7,148.9,213.5,215.1 ]
          END IF
        END IF
      CASE ('NO2       ','NO3       ')
        rsurf(7:13,n)=[ 386.6,405.4,490.2,555.6,2500.0,806.5,2500.0 ]
      CASE ('NO        ')
        rsurf(7:13,n)=[ 2319.6,2432.4,2941.2,3333.3,15000.0,4838.7,15000.0 ]
      CASE ('HNO3      ','HONO2     ','B2ndry    ','A2ndry    ','N2O5      ')
        rsurf(7:13,n)=[ 12.0,59.1,65.4,12.8,13.9,16.0,19.4 ]
      CASE ('HCl       ','HOCl      ','HBr       ','HOBr      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:13,n)=[ 12.0,59.1,65.4,12.8,13.9,16.0,19.4 ]
        END IF
      CASE ('H2SO4     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:13,n)=[ 120.0, 118.1, 130.7, 128.5, 138.6, 160.4, 194.2 ]
        END IF
      CASE ('HNO4      ','HO2NO2    ')
        rsurf(7:13,n)=[ 24.0,23.6,26.1,25.7,27.7,32.1,38.8 ]
      CASE ('HONO      ')
        rsurf(7:13,n)=[ 60.0,59.1,65.4,64.2,69.3,80.2,97.1 ]
      CASE ('H2O2      ','HOOH      ')
        rsurf(7:13,n)=[ 120.0,118.1,130.7,128.5,138.6,160.4,194.2 ]
        ! CRI PAN-type species
      CASE ('PAN       ','PPAN      ',                                         &
            'PHAN      ','RU12PAN   ','RTN26PAN  ')
        rsurf(7:13,n)=[ 515.5,540.5,653.6,740.7,3333.3,1075.3,3333.3 ]
      CASE ('MPAN      ')
        rsurf(7:13,n)=[ 1030.9,1081.1,1307.2,1481.5,6666.7,2150.5,6666.7 ]
        ! CRI organic-nitrates, copying ORGANIT rates
      CASE ('OnitU     ','ISON      ','ORGNIT    ',                            &
            'HOC2H4NO3 ','RTX24NO3  ','RN9NO3    ',                            &
            'RN12NO3   ','RN15NO3   ','RN18NO3   ','RU14NO3   ',               &
            'RTN28NO3  ','RTN25NO3  ','RTX28NO3  ','RTX22NO3  ',               &
            'RA22NO3   ','RA25NO3   ','RTN23NO3  ','RU12NO3   ',               &
            'RU10NO3   ','RA13NO3   ','RA16NO3   ','RA19NO3   ')
        rsurf(7:13,n)=[ 618.6,648.6,784.3,888.9,4000.0,1290.3,4000.0 ]
      CASE ('NH3       ')
        rsurf(7:13,n)=[ 196.1,185.8,196.1,180.7,148.9,213.5,215.1 ]
      CASE ('CO        ')
        rsurf(7:13,n)=[ 1960.0,4550.0,4550.0,r_null,r_null,4550.0,r_null ]
          ! Shrub+bare soil set to C3 grass (guess)
      CASE ('HCHO      ')
        rsurf(7:13,n)=[ 211.6,203.1,217.9,205.1,182.7,246.5,261.8 ]
        ! CRI v2.2 adds IEPOX, HMML, HUCARB9 and DHCARB9 to this list
      CASE ('MeOH      ','EtOH      ','i-PrOH    ','n-PrOH    ',               &
            'AROH14    ','ARNOH14   ','AROH17    ','ARNOH17   ',               &
            'IEPOX     ','HMML      ','HUCARB9   ','DHCARB9   ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:13,n)=[ 295.6, 282.7, 301.7, 282.1, 245.1, 337.2, 352.2 ]
        END IF
        ! CRI Carbonyls, copying MeCHO rates
        !   Second generation nitrates deposit as NALD
      CASE ('MeCHO     ','EtCHO     ','MACR      ','NALD      ',               &
            'HOCH2CHO  ',                                                      &
            'CARB14    ','CARB17    ','CARB11A   ','UCARB10   ',               &
            'UCARB12   ','UDCARB8   ','UDCARB11  ',                            &
            'UDCARB14  ','TNCARB26  ','TNCARB10  ','TNCARB12  ',               &
            'TNCARB11  ','CCARB12   ','TNCARB15  ','TXCARB24  ',               &
            'TXCARB22  ','UDCARB17  ','NOA       ','NUCARB12  ',               &
            'ANHY      '                                                       &
            )
        rsurf(7:13,n)=[ 6185.6,6486.5,7843.1,8888.9,40000.0,12903.2,40000.0 ]
        ! CRI Glyoxal-type species, deposit as MGLY
        !   n.b. CARB6 ~ MGLOX == MGLY
      CASE ('MGLY      ',                                                      &
            'CARB3     ','CARB6     ','CARB9     ','CARB12    ',               &
            'CARB15    '                                                       &
            )
        rsurf(7:13,n)=[ 19607.8,18575.9,19607.8,18072.3,14888.3,21352.3,       &
                         21505.4 ]
      END SELECT
    END DO
  CASE (17, 27)
    ! Standard surface resistances (s m-1). Values are for 17/27 tiles in
    ! order: C3 Crop, C3 Pasture, C4 Grass, C4 Crop, C4 Pasture,
    !        Shrub deciduous, Shrub evergreen, Urban, Water, Bare Soil, Ice.
    DO n = 1, ndepd
      SELECT CASE (speci(nldepd(n)))
      CASE ('O3        ','O3S       ')
        rsurf(7:17,n)=[ 355.0,355.0,309.3, 309.3,309.3,                        &
                         324.3,392.2,444.4,2000.0,645.2,2000.0 ]
      CASE ('SO2       ')
        IF (ukca_config%l_ukca_dry_dep_so2wet) THEN
          ! For the Erisman, Pul and Wyers (1994) parameterization, C3
          ! and C4 crops are assumed to be irrigated/watered, and so their
          ! resistances are greatly reduced (although the choice of 30.0
          ! is an estimate). Additionally we reduce the surface resistance
          ! value of water to 1.0 s m-1.
          rsurf(7:17,n)=[  30.0,209.8,196.1, 30.0,196.1,                       &
                           185.8,196.1,180.7,  1.0,213.5,215.1 ]
        ELSE
          IF (l_fix_drydep_so2_water) THEN
            ! Setting resistance of water to 148.9 s m-1 was an error, so
            ! this change reduces it its previous value (10.0 s m-1)
            rsurf(7:17,n)=[ 209.8,209.8,196.1,196.1,196.1,                     &
                             185.8,196.1,180.7, 10.0,213.5,215.1 ]
          ELSE
            ! Retains previous suspect value of 148.9 s m-1 for water
            ! resistance (used in UKESM1.0)
            rsurf(7:17,n)=[ 209.8,209.8,196.1,196.1,196.1,                     &
                             185.8,196.1,180.7,148.9,213.5,215.1 ]
          END IF
        END IF
      CASE ('NO2       ','NO3       ')
        rsurf(7:17,n)=[ 443.8,443.8,386.6, 386.6,386.6,                        &
                         405.4,490.2,555.6,2500.0,806.5,2500.0 ]
      CASE ('NO        ')
        rsurf(7:17,n)=[ 3662.7,3662.7,2319.6, 2319.6,2319.6,                   &
                         2432.4,2941.2,3333.3,15000.0,4838.7,15000.0 ]
      CASE ('HNO3      ','HONO2     ','B2ndry    ','A2ndry    ','N2O5      ')
        rsurf(7:17,n)=[ 13.2,13.2,12.0,12.0,12.0,                              &
                         59.1,65.4,12.8,13.9,16.0,19.4 ]
      CASE ('HCl       ','HOCl      ','HBr       ','HOBr      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:17,n)=[ 13.2,13.2,12.0,12.0,12.0,                            &
                           59.1,65.4,12.8,13.9,16.0,19.4 ]
        END IF
      CASE ('H2SO4     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:17,n)=[ 131.9, 131.9, 120.0, 120.0, 120.0,                   &
                           118.1, 130.7, 128.5, 138.6, 160.4, 194.2 ]
        END IF
      CASE ('HNO4      ','HO2NO2    ')
        rsurf(7:17,n)=[ 26.4,26.4,24.0,24.0,24.0,                              &
                         23.6,26.1,25.7,27.7,32.1,38.8 ]
      CASE ('HONO      ')
        rsurf(7:17,n)=[ 65.9,65.9,60.0,60.0,60.0,                              &
                         59.1,65.4,64.2,69.3,80.2,97.1 ]
      CASE ('H2O2      ','HOOH      ')
        rsurf(7:17,n)=[ 131.9,131.9,120.0,120.0,120.0,                         &
                         118.1,130.7,128.5,138.6,160.4,194.2 ]
        ! CRI PAN-type species
      CASE ('PAN       ','PPAN      ',                                         &
            'PHAN      ','RU12PAN   ','RTN26PAN  ')
        rsurf(7:17,n)=[ 591.7,591.7,515.5, 515.5, 515.5,                       &
                         540.5,653.6,740.7,3333.3,1075.3,3333.3 ]
      CASE ('MPAN      ')
        rsurf(7:17,n)=[ 1183.4,1183.4,1030.9,1030.9,1030.9,                    &
                         1081.1,1307.2,1481.5,6666.7,2150.5,6666.7 ]
        ! CRI organic-nitrates, copying ORGANIT rates
      CASE ('OnitU     ','ISON      ','ORGNIT    ',                            &
            'HOC2H4NO3 ','RTX24NO3  ','RN9NO3    ',                            &
            'RN12NO3   ','RN15NO3   ','RN18NO3   ','RU14NO3   ',               &
            'RTN28NO3  ','RTN25NO3  ','RTX28NO3  ','RTX22NO3  ',               &
            'RA22NO3   ','RA25NO3   ','RTN23NO3  ','RU12NO3   ',               &
            'RU10NO3   ','RA13NO3   ','RA16NO3   ','RA19NO3   ')
        rsurf(7:17,n)=[ 710.1,710.1,618.6, 618.6, 618.6,                       &
                         648.6,784.3,888.9,4000.0,1290.3,4000.0 ]
      CASE ('NH3       ')
        rsurf(7:17,n)=[ 209.8,209.8,196.1,196.1,196.1,                         &
                         185.8,196.1,180.7,148.9,213.5,215.1 ]
      CASE ('CO        ')
        rsurf(7:17,n)=[ 4550.0,4550.0,1960.0,1960.0,1960.0,                    &
                         4550.0,4550.0,r_null,r_null,4550.0,r_null ]
          ! Shrub+bare soil set to C3 grass (guess)
      CASE ('HCHO      ')
        rsurf(7:17,n)=[ 228.5,228.5,211.6,211.6,211.6,                         &
                         203.1,217.9,205.1,182.7,246.5,261.8 ]
        ! CRI alcohols deposit as MeOH
      CASE ('MeOH      ','EtOH      ','i-PrOH    ','n-PrOH    ',               &
            'AROH14    ','ARNOH14   ','AROH17    ','ARNOH17   ',               &
            'IEPOX     ','HMML      ','HUCARB9   ','DHCARB9   ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:17,n)=[ 318.3, 318.3, 295.6, 295.6, 295.6,                   &
                           282.7, 301.7, 282.1, 245.1, 337.2, 352.2 ]
        END IF
        ! CRI Carbonyls, copying MeCHO rates
        !   Second generation nitrates deposit as NALD
      CASE ('MeCHO     ','EtCHO     ','MACR      ','NALD      ',               &
            'HOCH2CHO  ',                                                      &
            'CARB14    ','CARB17    ','CARB11A   ','UCARB10   ',               &
            'UCARB12   ','UDCARB8   ','UDCARB11  ',                            &
            'UDCARB14  ','TNCARB26  ','TNCARB10  ','TNCARB12  ',               &
            'TNCARB11  ','CCARB12   ','TNCARB15  ','TXCARB24  ',               &
            'TXCARB22  ','UDCARB17  ','NOA       ','NUCARB12  ',               &
            'ANHY      '                                                       &
            )
        rsurf(7:17,n)=[ 7100.6,7100.6,6185.6, 6185.6, 6185.6,                  &
                         6486.5,7843.1,8888.9,40000.0,12903.2,40000.0 ]
        ! CRI Glyoxal-type species, deposit as MGLY
        !   n.b. CARB6 ~ MGLOX == MGLY
      CASE ('MGLY      ',                                                      &
            'CARB3     ','CARB6     ','CARB9     ','CARB12    ',               &
            'CARB15    '                                                       &
            )
        rsurf(7:17,n)=[ 20979.0,20979.0,19607.8,19607.8,19607.8,               &
                         18575.9,19607.8,18072.3,14888.3,21352.3,21505.4 ]
      END SELECT
    END DO
  END SELECT

  IF (ntype ==27) THEN
    ! Standard surface resistances (s m-1). Values are for 27 tiles in
    ! order: Elev_Ice(1-10).
    DO n = 1, ndepd
      SELECT CASE (speci(nldepd(n)))
      CASE ('O3        ','O3S       ')
        rsurf(18:27,n)=[ 2000.0,2000.0,2000.0,2000.0,2000.0,                   &
                          2000.0,2000.0,2000.0,2000.0,2000.0 ]
      CASE ('SO2       ')
        rsurf(18:27,n)=[ 215.1,215.1,215.1,215.1,215.1,                        &
                          215.1,215.1,215.1,215.1,215.1 ]
      CASE ('NO2       ','NO3       ')
        rsurf(18:27,n)=[ 2500.0,2500.0,2500.0,2500.0,2500.0,                   &
                          2500.0,2500.0,2500.0,2500.0,2500.0 ]
      CASE ('NO        ')
        rsurf(18:27,n)=[ 15000.0,15000.0,15000.0,15000.0,15000.0,              &
                          15000.0,15000.0,15000.0,15000.0,15000.0 ]
      CASE ('HNO3      ','HONO2     ','B2ndry    ','A2ndry    ','N2O5      ')
        rsurf(18:27,n)=[ 19.4,19.4,19.4,19.4,19.4,                             &
                          19.4,19.4,19.4,19.4,19.4 ]
      CASE ('HCl       ','HOCl      ','HBr       ','HOBr      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(18:27,n)=[ 19.4,19.4,19.4,19.4,19.4,                           &
                            19.4,19.4,19.4,19.4,19.4 ]
        END IF
      CASE ('H2SO4     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(18:27,n)=[ 194.2, 194.2, 194.2, 194.2, 194.2,                  &
                            194.2, 194.2, 194.2, 194.2, 194.2 ]
        END IF
      CASE ('HNO4      ','HO2NO2    ')
        rsurf(18:27,n)=[ 38.8,38.8,38.8,38.8,38.8,                             &
                          38.8,38.8,38.8,38.8,38.8 ]
      CASE ('HONO      ')
        rsurf(18:27,n)=[ 97.1,97.1,97.1,97.1,97.1,                             &
                          97.1,97.1,97.1,97.1,97.1 ]
      CASE ('H2O2      ','HOOH      ')
        rsurf(18:27,n)=[ 194.2,194.2,194.2,194.2,194.2,                        &
                          194.2,194.2,194.2,194.2,194.2 ]
        ! CRI PAN-type species
      CASE ('PAN       ','PPAN      ',                                         &
            'PHAN      ','RU12PAN   ','RTN26PAN  ')
        rsurf(18:27,n)=[ 3333.3,3333.3,3333.3,3333.3,3333.3,                   &
                          3333.3,3333.3,3333.3,3333.3,3333.3 ]
      CASE ('MPAN      ')
        rsurf(18:27,n)=[ 6666.7,6666.7,6666.7,6666.7,6666.7,                   &
                          6666.7,6666.7,6666.7,6666.7,6666.7 ]
        ! CRI organic-nitrates, copying ORGANIT rates
      CASE ('OnitU     ','ISON      ','ORGNIT    ',                            &
            'HOC2H4NO3 ','RTX24NO3  ','RN9NO3    ',                            &
            'RN12NO3   ','RN15NO3   ','RN18NO3   ','RU14NO3   ',               &
            'RTN28NO3  ','RTN25NO3  ','RTX28NO3  ','RTX22NO3  ',               &
            'RA22NO3   ','RA25NO3   ','RTN23NO3  ','RU12NO3   ',               &
            'RU10NO3   ','RA13NO3   ','RA16NO3   ','RA19NO3   ')
        rsurf(18:27,n)=[ 4000.0,4000.0,4000.0,4000.0,4000.0,                   &
                          4000.0,4000.0,4000.0,4000.0,4000.0 ]
      CASE ('NH3       ')
        rsurf(18:27,n)=[ 215.1,215.1,215.1,215.1,215.1,                        &
                          215.1,215.1,215.1,215.1,215.1 ]
      CASE ('CO        ')
        rsurf(18:27,n)=[ r_null,r_null,r_null,r_null,r_null,                   &
                          r_null,r_null,r_null,r_null,r_null ]
          ! Shrub+bare soil set to C3 grass (guess)
      CASE ('HCHO      ')
        rsurf(18:27,n)=[ 261.8,261.8,261.8,261.8,261.8,                        &
                          261.8,261.8,261.8,261.8,261.8 ]
        ! CRI alcohols deposit as MeOH
      CASE ('MeOH      ','EtOH      ','i-PrOH    ','n-PrOH    ',               &
            'AROH14    ','ARNOH14   ','AROH17    ','ARNOH17   ',               &
            'IEPOX     ','HMML      ','HUCARB9   ','DHCARB9   ')
        IF (l_fix_improve_drydep) THEN
          rsurf(18:27,n)=[ 352.2, 352.2, 352.2, 352.2, 352.2,                  &
                            352.2, 352.2, 352.2, 352.2, 352.2 ]
        END IF
        ! CRI Carbonyls, copying MeCHO rates
        !   Second generation nitrates deposit as NALD
      CASE ('MeCHO     ','EtCHO     ','MACR      ','NALD      ',               &
            'HOCH2CHO  ',                                                      &
            'CARB14    ','CARB17    ','CARB11A   ','UCARB10   ',               &
            'UCARB12   ','UDCARB8   ','UDCARB11  ',                            &
            'UDCARB14  ','TNCARB26  ','TNCARB10  ','TNCARB12  ',               &
            'TNCARB11  ','CCARB12   ','TNCARB15  ','TXCARB24  ',               &
            'TXCARB22  ','UDCARB17  ','NOA       ','NUCARB12  ',               &
            'ANHY      '                                                       &
            )
        rsurf(18:27,n)=[ 40000.0,40000.0,40000.0,40000.0,40000.0,              &
                          40000.0,40000.0,40000.0,40000.0,40000.0 ]
        ! CRI Glyoxal-type species, deposit as MGLY
        !   n.b. CARB6 ~ MGLOX == MGLY
      CASE ('MGLY      ',                                                      &
            'CARB3     ','CARB6     ','CARB9     ','CARB12    ',               &
            'CARB15    '                                                       &
            )
        rsurf(18:27,n)=[ 21505.4,21505.4,21505.4,21505.4,21505.4,              &
                          21505.4,21505.4,21505.4,21505.4,21505.4 ]
      END SELECT
    END DO
  END IF
  first = .FALSE.
END IF !first

o3_stom_frac = 0.0

!     rco3global = 0.0
!     gsfglobal = 0.0


!     Set logical for surface types

DO n = 1, ntype
  DO k = 1, rows
    DO i = 1, row_length
      todo(i,k,n) = (gsf(i,k,n) > min_tile_frac)
    END DO
  END DO
END DO

! Set microb
microb(:,:) = (rh(:,:) > 0.4 .AND. t0(:,:) > 278.0)

!     Set surface resistances to standard values. rsurf is the
!     resistance of the soil, rock, water etc. Set all tiles to
!     standard values. These values will be modified below as
!     necessary. Extra terms for vegetated tiles (stomatal, cuticular)
!     will be added if required. Loop over all parts of array rc to
!     ensure all of it is assigned a value.

DO n = 1, ntype
  DO k = 1, rows
    DO i = 1, row_length
      p_surf1 = p0(i,k)
      t_surf1 = t0(i,k)
      ustar = u_s(i,k)
      sst   = t0(i,k) ! assume t0 is representative of SST (K)
      IF (todo(i,k,n)) THEN
        DO j = 1, ndepd
          rc(i,k,n,j) = rsurf(n,j)

          IF (ukca_config%l_ukca_ddepo3_ocean) THEN
            ! Call ozone dry deposition scheme for ocean surface, for both
            ! o3 and o3s. There are 4 ntype cases (i.e. 9, 13, 17 and 27)

            IF ((ntype == 9 .AND. n == 7)   .OR.                               &
                (ntype == 13 .AND. n == 11) .OR.                               &
                (ntype == 17 .AND. n == 15) .OR.                               &
                (ntype == 27 .AND. n == 15)) THEN

              IF (((speci(nldepd(j)) == 'O3        ') .OR.                     &
                   (speci(nldepd(j)) == 'O3S       ')) .AND.                   &
                   gsf(i,k,n) > 0.75) THEN
                ! gsf > 0.75 is to make sure that the bulk of the grid is
                ! water to invoke the new scheme
                CALL ukca_ddepo3_ocean(p_surf1,t_surf1,sst,ustar,rc_ocean)
                rc(i,k,n,j) = rc_ocean
              END IF
            ELSE IF ((ntype /= 9)  .AND. (ntype /= 13)  .AND.                  &
                     (ntype /= 17) .AND. (ntype /= 27)) THEN
              ! ntype must equal 9, 13, 17, 27
              WRITE(cmessage,'(A)')                                            &
                   'UKCA does not handle flexible tiles yet: ' //              &
                   'Dry deposition needs extending. '          //              &
                   'Please use standard tile configuration'
              errcode=1
              CALL ereport(RoutineName,errcode,cmessage)
            END IF
          END IF
        END DO
      ELSE
        DO j = 1, ndepd
          rc(i,k,n,j) = r_null
        END DO
      END IF
    END DO
  END DO
END DO

!     Calculate stomatal resistances

DO n = 1, npft
  DO k = 1, rows
    DO i = 1, row_length
      IF (todo(i,k,n) .AND. stcon(i,k,n) > glmin) THEN
        r_stom(i,k,n,r_stom_no2) = 1.5 * dif(r_stom_no2) / stcon(i,k,n) ! NO2
        r_stom(i,k,n,r_stom_o3)  =       dif(r_stom_o3)  / stcon(i,k,n) ! O3
        r_stom(i,k,n,r_stom_pan) =       dif(r_stom_pan) / stcon(i,k,n) ! PAN
        r_stom(i,k,n,r_stom_so2) =       dif(r_stom_so2) / stcon(i,k,n) ! SO2
        r_stom(i,k,n,r_stom_nh3) =       dif(r_stom_nh3) / stcon(i,k,n) ! NH3
      ELSE
        r_stom(i,k,n,r_stom_no2) = r_null
        r_stom(i,k,n,r_stom_o3)  = r_null
        r_stom(i,k,n,r_stom_pan) = r_null
        r_stom(i,k,n,r_stom_so2) = r_null
        r_stom(i,k,n,r_stom_nh3) = r_null
      END IF
    END DO
  END DO
END DO

!     Now begin assigning specific surface resistances.

DO j = 1, ndepd

  !       O3: Change land deposition values if surface is wet;
  !       soil moisture value of 0.3 fairly arbitrary.

  IF (speci(nldepd(j)) == 'O3        ') THEN
    DO k = 1, rows
      DO i = 1, row_length
        IF (smr(i,k) > soil_moist_thresh) THEN
          DO n = 1, npft
            IF (todo(i,k,n))   rc(i,k,n,j) = r_wet_o3
          END DO
          IF (todo(i,k,soil))  rc(i,k,soil,j) = r_wet_o3
        END IF
      END DO

      !           Change values for tundra regions
      ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
      ! When l_fix_ukca_h2dd_x is retired, this could become a CASE statement
      IF ( ( .NOT. l_fix_ukca_h2dd_x ) .OR. ( ntype == 9) ) THEN
        IF ( l_fix_ukca_h2dd_x ) THEN
          n = shrub
        ELSE
          n = npft
        END IF
        DO i = 1, row_length
          IF (sinlat(i,k) > tundra_s_limit) THEN
            IF (todo(i,k,n))     rc(i,k,n,j)     = r_tundra(r_tundra_o3)
          END IF
        END DO
      ELSE IF ( (ntype == 13) .OR. (ntype == 17) .OR. (ntype == 27) ) THEN
        DO n = shrub_dec, shrub_eg
          DO i = 1, row_length
            IF (sinlat(i,k) > tundra_s_limit) THEN
              IF (todo(i,k,n))   rc(i,k,n,j)     = r_tundra(r_tundra_o3)
            END IF
          END DO
        END DO
      ELSE
        ! ntype must equal 9, 13, 17, 27
        WRITE(cmessage,'(A)')                                                  &
             'UKCA does not handle flexible tiles yet: ' //                    &
             'Dry deposition needs extending. '          //                    &
             'Please use standard tile configuration'
        errcode=1
        CALL ereport(RoutineName,errcode,cmessage)
      END IF

      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,soil))    rc(i,k,soil,j)  = r_tundra(r_tundra_o3)
        END IF
      END DO

      !           Cuticular resistance for ozone.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n) .AND. lai_ft(i,k,n) > 0.0) THEN
            r_cut_o3(i,n) = cuticle_o3 / lai_ft(i,k,n)
          ELSE
            r_cut_o3(i,n) = r_null
          END IF
        END DO
      END DO

      !           Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,r_stom_o3)) +                               &
                 (1.0/r_cut_o3(i,n)) +                                         &
                 (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
            o3_stom_frac(i,k) = o3_stom_frac(i,k) +                            &
              gsf(i,k,n) / (rr * r_stom(i,k,n,r_stom_o3))
          END IF
        END DO
      END DO
    END DO

  ELSE IF (speci(nldepd(j)) == 'SO2       ') THEN

    ! The parametrization below is based on `Parametrization of surface
    ! resistance for the quantification of atmospheric deposition of
    ! acidifying pollutants and ozone' by Erisman, Pul and Wyers (1993)
    ! (hereafter EPW1994) and `The Elspeetsche Veld experiment on surface
    ! exchange of trace gases: Summary of results' by Erisman et al (1994).
    ! The surface resistance for S02 deposition is modified depending
    ! on if the surface is wet or dry or it is raining. If dry and no rain,
    ! make rc a function of surface relative humidity. Also calculate the
    ! cuticular resistance term for SO2 and include both r_cut_so2 and
    ! r_stom in final calculated value for SO2.
    ! Snow on veg/soil does not appear to be considered for any
    ! deposition term, so it is not included in the parameterization below.
    ! It is partly accounted for in the r_cut_so2 term via the use of t0tile.
    IF (ukca_config%l_ukca_dry_dep_so2wet) THEN

      ! Change to dry deposition for soil tile
      IF (soil > 0) THEN
        DO k = 1, rows
          DO i = 1, row_length
            IF (smr(i,k) > soil_moist_thresh .OR. surf_wetness(i,k) > tol) THEN
              ! Surface is wet
              IF (t0tile(i,k,soil) > minus1degc) THEN
                IF (todo(i,k,soil)) THEN
                  rc(i,k,soil,j) = r_wet_so2
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF

      ! Change to dry deposition for vegetation tiles
      DO n = 1, npft
        DO k = 1, rows
          DO i = 1, row_length
            IF (todo(i,k,n)) THEN

              ! Change land deposition values
              IF (smr(i,k) > soil_moist_thresh .OR.                            &
                  surf_wetness(i,k) > tol) THEN
                ! Surface is wet
                IF (t0tile(i,k,n) > minus1degc) THEN
                  rc(i,k,n,j) = r_wet_so2
                END IF
              END IF

              ! Cuticular resistance for SO2
              IF (lai_ft(i,k,n) > tol) THEN
                ! Determine the temperature range
                IF (t0tile(i,k,n) >= minus1degc) THEN
                  ! Temperature is greater than -1 degrees C
                  IF (smr(i,k) > soil_moist_thresh .OR.                        &
                      surf_wetness(i,k) > tol) THEN
                    ! Surface is wet
                    r_cut_so2 = r_wet_so2
                  ELSE
                    ! Equation 9 from EPW1994 is not applied to crops
                    IF ( (n == c3_crop) .OR. (n == c4_crop) ) THEN
                      r_cut_so2 = rc(i,k,n,j)
                    ELSE
                      IF (rh(i,k) < rh_thresh) THEN
                        r_cut_so2 = so2con1 * EXP(so2con3*rh(i,k))
                      ELSE
                        r_cut_so2 = so2con2 * EXP(so2con4*rh(i,k))
                      END IF
                    END IF

                    ! Ensure that r_cut_so2 stays within reasonable bounds
                    ! i.e. r_wet_so2 <= r_cut_so2 <= r_dry_so2
                    r_cut_so2 = MAX( r_wet_so2, MIN(r_dry_so2, r_cut_so2) )
                  END IF
                ELSE IF (t0tile(i,k,n) >= minus5degc) THEN
                  ! Temperature is between -5 degrees C and -1 degrees C
                  r_cut_so2 = r_cut_so2_5to_1
                ELSE
                  ! Temperature is < -5 degrees C
                  r_cut_so2 = r_cut_so2_5degc
                END IF
              ELSE
                r_cut_so2 = r_null
              END IF

              ! Include plant deposition and cuticular term in rc calculation
              rr = (1.0/r_stom(i,k,n,r_stom_so2)) + (1.0/r_cut_so2) +          &
                   (1.0/rc(i,k,n,j))
              rc(i,k,n,j) = 1.0 / rr

            END IF
          END DO
        END DO
      END DO

    END IF

    !       NO2

  ELSE IF (speci(nldepd(j)) == 'NO2       ') THEN

    DO k = 1, rows

      !           Change values for tundra regions
      ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
      ! When l_fix_ukca_h2dd_x is retired, this could become a CASE statement
      IF ( ( .NOT. l_fix_ukca_h2dd_x ) .OR. ( ntype == 9 ) ) THEN
        IF ( l_fix_ukca_h2dd_x ) THEN
          n = shrub
        ELSE
          n = npft
        END IF

        DO i = 1, row_length
          IF (sinlat(i,k) > tundra_s_limit) THEN
            IF (todo(i,k,n))    rc(i,k,n,j)    = r_tundra(r_tundra_no2)
          END IF
        END DO
      ELSE IF ( (ntype == 13) .OR. (ntype == 17) .OR. (ntype == 27) ) THEN
        DO n = shrub_dec, shrub_eg
          DO i = 1, row_length
            IF (sinlat(i,k) > tundra_s_limit) THEN
              IF (todo(i,k,n))    rc(i,k,n,j)    = r_tundra(r_tundra_no2)
            END IF
          END DO
        END DO
      ELSE
        ! ntype must equal 9, 13, 17, 27
        WRITE(cmessage,'(A)')                                                  &
             'UKCA does not handle flexible tiles yet: ' //                    &
             'Dry deposition needs extending. '          //                    &
             'Please use standard tile configuration'
        errcode=1
        CALL ereport(RoutineName,errcode,cmessage)
      END IF

      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,soil)) rc(i,k,soil,j) = r_tundra(r_tundra_no2)
        END IF
      END DO

      !           Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,r_stom_no2)) + (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
          END IF
        END DO
      END DO
    END DO

    !       PAN

  ELSE IF (speci(nldepd(j)) == 'PAN       ') THEN

    DO k = 1, rows

      !           Change values for tundra regions
      ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
      ! When l_fix_ukca_h2dd_x is retired, this could become a CASE statement
      IF ( ( .NOT. l_fix_ukca_h2dd_x ) .OR. ( ntype == 9 ) ) THEN
        IF ( l_fix_ukca_h2dd_x ) THEN
          n = shrub
        ELSE
          n = npft
        END IF

        DO i = 1, row_length
          IF (sinlat(i,k) > tundra_s_limit) THEN
            IF (todo(i,k,n))     rc(i,k,n,j)    = r_tundra(r_tundra_pan)
          END IF
        END DO
      ELSE IF ( (ntype == 13) .OR. (ntype == 17) .OR. (ntype == 27) ) THEN
        DO n = shrub_dec, shrub_eg
          DO i = 1, row_length
            IF (sinlat(i,k) > tundra_s_limit) THEN
              IF (todo(i,k,n))   rc(i,k,n,j)    = r_tundra(r_tundra_pan)
            END IF
          END DO
        END DO
      ELSE
        ! ntype must equal 9, 13, 17, 27
        WRITE(cmessage,'(A)')                                                  &
             'UKCA does not handle flexible tiles yet: ' //                    &
             'Dry deposition needs extending. '          //                    &
             'Please use standard tile configuration'
        errcode=1
        CALL ereport(RoutineName,errcode,cmessage)
      END IF

      ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,soil))  rc(i,k,soil,j) = r_tundra(r_tundra_pan)
        END IF
      END DO

      !           Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,r_stom_pan)) + (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
          END IF
        END DO
      END DO
    END DO

    !         PPAN

  ELSE IF (speci(nldepd(j)) == 'PPAN      ') THEN

    DO k = 1, rows

      !               Change values for tundra regions
      ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
      ! When l_fix_ukca_h2dd_x is retired, this could become a CASE statement
      IF ( ( .NOT. l_fix_ukca_h2dd_x ) .OR. ( ntype == 9 ) ) THEN
        IF ( l_fix_ukca_h2dd_x ) THEN
          n = shrub
        ELSE
          n = npft
        END IF
        DO i = 1, row_length
          IF (sinlat(i,k) > tundra_s_limit) THEN
            IF (todo(i,k,n))     rc(i,k,n,j)    = r_tundra(r_tundra_pan)
          END IF
        END DO
      ELSE IF ( (ntype == 13) .OR. (ntype == 17) .OR. (ntype == 27) ) THEN
        DO n = shrub_dec, shrub_eg
          DO i = 1, row_length
            IF (sinlat(i,k) > tundra_s_limit) THEN
              IF (todo(i,k,n))     rc(i,k,n,j)    = r_tundra(r_tundra_pan)
            END IF
          END DO
        END DO
      ELSE
        ! ntype must equal 9, 13, 17, 27
        WRITE(cmessage,'(A)')                                                  &
             'UKCA does not handle flexible tiles yet: ' //                    &
             'Dry deposition needs extending. '          //                    &
             'Please use standard tile configuration'
        errcode=1
        CALL ereport(RoutineName,errcode,cmessage)
      END IF

      ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,soil))  rc(i,k,soil,j) = r_tundra(r_tundra_pan)
        END IF
      END DO

      !               Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,r_stom_pan)) + (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
          END IF
        END DO
      END DO
    END DO

    !         MPAN

  ELSE IF (speci(nldepd(j)) == 'MPAN      ') THEN

    DO k = 1, rows

      !             Change values for tundra regions
      ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
      ! When l_fix_ukca_h2dd_x is retired, this could become a CASE statement
      IF ( ( .NOT. l_fix_ukca_h2dd_x ) .OR. ( ntype == 9 ) ) THEN
        IF ( l_fix_ukca_h2dd_x ) THEN
          n = shrub
        ELSE
          n = npft
        END IF

        DO i = 1, row_length
          IF (sinlat(i,k) > tundra_s_limit) THEN
            IF (todo(i,k,n))     rc(i,k,n,j)    = r_tundra(r_tundra_pan)
          END IF
        END DO
      ELSE IF ( (ntype == 13) .OR. (ntype == 17) .OR. (ntype == 27) ) THEN
        DO n = shrub_dec, shrub_eg
          DO i = 1, row_length
            IF (sinlat(i,k) > tundra_s_limit) THEN
              IF (todo(i,k,n))     rc(i,k,n,j)    = r_tundra(r_tundra_pan)
            END IF
          END DO
        END DO
      ELSE
        ! ntype must equal 9, 13, 17, 27
        WRITE(cmessage,'(A)')                                                  &
             'UKCA does not handle flexible tiles yet: ' //                    &
             'Dry deposition needs extending. '          //                    &
             'Please use standard tile configuration'
        errcode=1
        CALL ereport(RoutineName,errcode,cmessage)
      END IF

      ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,soil))  rc(i,k,soil,j) = r_tundra(r_tundra_pan)
        END IF
      END DO

      !             Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,r_stom_pan)) + (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
          END IF
        END DO
      END DO
    END DO

    !         ONITU

  ELSE IF (speci(nldepd(j)) == 'ONITU     ') THEN

    DO k = 1, rows

      !             Change values for tundra regions
      ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
      ! When l_fix_ukca_h2dd_x is retired, this could become a CASE statement
      IF ( ( .NOT. l_fix_ukca_h2dd_x ) .OR. ( ntype == 9 ) ) THEN
        IF ( l_fix_ukca_h2dd_x ) THEN
          n = shrub
        ELSE
          n = npft
        END IF

        DO i = 1, row_length
          IF (sinlat(i,k) > tundra_s_limit) THEN
            IF (todo(i,k,n))     rc(i,k,n,j)    = r_tundra(r_tundra_pan)
          END IF
        END DO
      ELSE IF ( (ntype == 13) .OR. (ntype == 17) .OR. (ntype == 27) ) THEN
        DO n = shrub_dec, shrub_eg
          DO i = 1, row_length
            IF (sinlat(i,k) > tundra_s_limit) THEN
              IF (todo(i,k,n))     rc(i,k,n,j)    = r_tundra(r_tundra_pan)
            END IF
          END DO
        END DO
      ELSE
        ! ntype must equal 9, 13, 17, 27
        WRITE(cmessage,'(A)')                                                  &
             'UKCA does not handle flexible tiles yet: ' //                    &
             'Dry deposition needs extending. '          //                    &
             'Please use standard tile configuration'
        errcode=1
        CALL ereport(RoutineName,errcode,cmessage)
      END IF

      ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,soil))  rc(i,k,soil,j) = r_tundra(r_tundra_pan)
        END IF
      END DO

      !             Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,r_stom_pan)) + (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
          END IF
        END DO
      END DO
    END DO

  ELSE IF (speci(nldepd(j)) == 'H2        ') THEN
    !       H2 dry dep vel has linear dependence on soil moisture
    !       Limit sm to avoid excessively high deposition velocities
    ! When l_fix_ukca_h2dd_x is retired, this could become a CASE statement
    IF ( ( .NOT. l_fix_ukca_h2dd_x ) .OR. ( ntype == 9 ) ) THEN
      DO k = 1, rows

        !  brd_leaf : H2 dry dep
        !  ndl_leaf : H2 dry dep
        !  c3_grass : H2 dry dep
        IF ( l_fix_ukca_h2dd_x ) THEN
          DO n = brd_leaf, c3_grass ! DO n = 1,3
            DO i = 1, row_length
              IF (todo(i,k,n) .AND. microb(i,k)) THEN
                sm = MAX(smr(i,k),0.1)
                rr = ( h2dd_m(n) * sm + h2dd_c(n) ) * 1.0e-4
                rc(i,k,n,j) = 1.0 / rr
              END IF
            END DO
          END DO
        ELSE
          DO n = 1, npft-2 ! DO n = 1,3 (9 tile)
            DO i = 1, row_length
              IF (todo(i,k,n) .AND. microb(i,k)) THEN
                sm = MAX(smr(i,k),0.1)
                rc(i,k,n,j) = 1.0 / (h2dd_m(n) * sm + h2dd_c(n))
              END IF
            END DO
          END DO
        END IF

        !  C4 grass: H2 dry dep has quadratic-log dependence on soil moisture
        IF ( l_fix_ukca_h2dd_x ) THEN
          n = c4_grass ! n = 4
          DO i = 1, row_length
            IF (todo(i,k,n) .AND. microb(i,k)) THEN
              sm = LOG(MAX(smr(i,k),0.1))
              rr = (          h2dd_c(c4_grass) +                               &
                     (sm * (  h2dd_m(c4_grass) +                               &
                     (sm *    h2dd_q(c4_grass)) )) ) * 1.0e-4
              IF (rr > 0.00131) rr = 0.00131 ! Conrad/Seiler Max value
              rc(i,k,n,j) = 1.0 / rr
            END IF
          END DO
        ELSE
          n = npft-1   ! n = 4 (9 tile)
          DO i = 1, row_length
            IF (todo(i,k,n) .AND. microb(i,k)) THEN
              sm = LOG(MAX(smr(i,k),0.1))
              rr = (          h2dd_c(c4_grass) +                               &
                     (sm * (  h2dd_m(c4_grass) +                               &
                     (sm *    h2dd_q(c4_grass)) )) ) * 1.0e-4
              IF (rr > 0.00131) rr = 0.00131 ! Conrad/Seiler Max value
              rc(i,k,n,j) = 1.0 / rr
            END IF
          END DO
        END IF


        ! Shrub: H2 dry dep velocity has no dependence on soil moisture
        ! Shrub and bare soil are assumed to be tundra if latitude > 60N;

        IF ( l_fix_ukca_h2dd_x ) THEN
          n = shrub ! n = 5
          DO i = 1, row_length
            IF (todo(i,k,n) .AND. microb(i,k)) THEN
              IF (sinlat(i,k) > tundra_s_limit) THEN
                ! treat shrub as shrub above tundra limit
                rr = h2dd_c(shrub) * 1.0e-4
                rc(i,k,n,j) = 1.0 / rr
              ELSE
                ! treat shrub as c3_grass below tundra limit
                sm = MAX(smr(i,k),0.1)
                rr = ( h2dd_m(c3_grass) * sm + h2dd_c(c3_grass) ) * 1.0e-4
                rc(i,k,n,j) = 1.0 / rr
              END IF
            END IF
          END DO
        ELSE
          n = npft ! n = 5 (9 tile)
          DO i = 1, row_length
            IF (todo(i,k,n) .AND. microb(i,k)) THEN
              IF (sinlat(i,k) > tundra_s_limit) THEN
                ! treat shrub as shrub above tundra limit
                rr = 1.0 / h2dd_c(n)
                rc(i,k,n,j) = rr
              ELSE
                ! treat shrub as c3_grass below tundra limit
                sm = MAX(smr(i,k),0.1)
                rc(i,k,n,j) = 1.0 / (h2dd_m(c3_grass) * sm + h2dd_c(c3_grass))
              END IF
            END IF
          END DO
        END IF

        ! Bare soil : H2 dry dep velocity has no dependence on soil moisture
        ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
        IF ( l_fix_ukca_h2dd_x ) THEN
          n = soil ! n = 8
          DO i = 1, row_length
            IF (todo(i,k,n) .AND. microb(i,k)) THEN
              IF (sinlat(i,k) > tundra_s_limit) THEN
                ! treat soil as shrub above tundra limit
                rr = 1.0 / h2dd_c(shrub) * 1.0e-4
                rc(i,k,n,j) = 1.0 / rr
              ELSE
                ! treat soil as c3_grass below tundra limit
                sm = MAX(smr(i,k),0.1)
                rr = ( h2dd_m(c3_grass)*sm + h2dd_c(c3_grass) ) * 1.0e-4
                rc(i,k,n,j) = 1.0 / rr
              END IF
            END IF
          END DO
        ELSE
          n = npft ! n = 5 (shrub not soil)
          DO i = 1, row_length
            IF (todo(i,k,npft) .AND. microb(i,k)) THEN  ! (shrub not soil)
              IF (sinlat(i,k) > tundra_s_limit) THEN
                ! treat soil as shrub above tundra limit
                rr = 1.0 / h2dd_c(npft) ! (shrub not soil)
                rc(i,k,soil,j) = rr
              ELSE
                ! treat soil as c3_grass below tundra limit
                sm = MAX(smr(i,k),0.1)
                rc(i,k,soil,j) = 1.0 / (h2dd_m(c3_grass)*sm + h2dd_c(c3_grass))
              END IF
            END IF
          END DO
        END IF

      END DO ! DO k = 1, rows

    ELSE IF ( ntype == 13 ) THEN

      DO k = 1, rows

        ! brd_leaf_dec      : H2 dry dep
        ! brd_leaf_eg_trop  : H2 dry dep
        ! brd_leaf_eg_temp  : H2 dry dep
        ! ndl_leaf_dec      : H2 dry dep
        ! ndl_leaf_eg       : H2 dry dep
        ! c3_grass          : H2 dry dep

        DO n = brd_leaf_dec, c3_grass ! DO n = 1,6
          DO i = 1, row_length
            IF (todo(i,k,n) .AND. microb(i,k)) THEN
              sm = MAX(smr(i,k),0.1)
              rr = ( h2dd_m(n) * sm + h2dd_c(n) ) * 1.0e-4
              rc(i,k,n,j) = 1.0 / rr
            END IF
          END DO
        END DO

        ! C4 grass: H2 dry dep has quadratic-log dependence on soil moisture

        n = c4_grass ! n = 7
        DO i = 1, row_length
          IF (todo(i,k,n) .AND. microb(i,k)) THEN
            sm = LOG(MAX(smr(i,k),0.1))
            rr = (          h2dd_c(c4_grass) +                                 &
                   (sm * (  h2dd_m(c4_grass) +                                 &
                   (sm *    h2dd_q(c4_grass)) )) ) * 1.0e-4
            IF (rr > 0.00131) rr = 0.00131 ! Conrad/Seiler Max value
            rc(i,k,n,j) = 1.0 / rr
          END IF
        END DO

        ! Shrub dec : H2 dry dep velocity has no dependence on soil moisture
        ! Shrub eg  : H2 dry dep velocity has no dependence on soil moisture
        ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
        DO n = shrub_dec , shrub_eg ! DO n = 8,9

          DO i = 1, row_length
            IF (todo(i,k,n) .AND. microb(i,k)) THEN
              IF (sinlat(i,k) > tundra_s_limit) THEN
                ! treat soil as shrub_dec above tundra limit
                rr = h2dd_c(shrub_dec) * 1.0e-4
                rc(i,k,n,j) = 1.0 / rr
              ELSE IF (sinlat(i,k) <= tundra_s_limit) THEN
                ! treat shrub as c3_grass below tundra limit
                sm = MAX(smr(i,k),0.1)
                rr = ( h2dd_m(c3_grass)*sm + h2dd_c(c3_grass) ) * 1.0e-4
                rc(i,k,n,j) = 1.0 / rr
              END IF
            END IF
          END DO
        END DO

        ! Bare soil :
        ! Shrub and bare soil are assumed to be tundra if latitude > 60N;

        n = soil ! n = 12

        DO i = 1, row_length
          IF (todo(i,k,n) .AND. microb(i,k)) THEN
            IF (sinlat(i,k) > tundra_s_limit) THEN
              ! treat soil as shrub_dec above tundra limit
              rr = h2dd_c(shrub_dec) * 1.0e-4
              rc(i,k,soil,j) = 1.0 / rr
            ELSE IF (sinlat(i,k) <= tundra_s_limit) THEN
              ! treat soil as c3_grass below tundra limit
              sm = MAX(smr(i,k),0.1)
              rr = ( h2dd_m(c3_grass) * sm + h2dd_c(c3_grass) ) * 1.0e-4
              rc(i,k,n,j) = 1.0 / rr
            END IF
          END IF
        END DO

      END DO ! DO k = 1, rows

    ELSE IF ( (ntype == 17) .OR. (ntype == 27) ) THEN

      DO k = 1, rows

        ! brd_leaf_dec      : H2 dry dep
        ! brd_leaf_eg_trop  : H2 dry dep
        ! brd_leaf_eg_temp  : H2 dry dep
        ! ndl_leaf_dec      : H2 dry dep
        ! ndl_leaf_eg       : H2 dry dep
        ! c3_grass          : H2 dry dep
        ! c3_crop           : H2 dry dep
        ! c3_pasture        : H2 dry dep

        DO n = brd_leaf_dec, c3_pasture ! DO n = 1,8
          DO i = 1, row_length
            IF (todo(i,k,n) .AND. microb(i,k)) THEN
              sm = MAX(smr(i,k),0.1)
              rr = ( h2dd_m(n) * sm + h2dd_c(n) ) * 1.0e-4
              rc(i,k,n,j) = 1.0 / rr
            END IF
          END DO
        END DO

        ! c4_grass   : H2 dry dep has quadratic-log dependence on soil moisture
        ! c4_crop    : H2 dry dep has quadratic-log dependence on soil moisture
        ! c4_pasture : H2 dry dep has quadratic-log dependence on soil moisture

        DO n = c4_grass , c4_pasture ! DO n = 9, 11
          DO i = 1, row_length
            IF (todo(i,k,n) .AND. microb(i,k)) THEN
              sm = LOG(MAX(smr(i,k),0.1))
              rr = (          h2dd_c(n) +                                      &
                     (sm * (  h2dd_m(n) +                                      &
                     (sm *    h2dd_q(n)) )) ) * 1.0e-4
              IF (rr > 0.00131) rr = 0.00131 ! Conrad/Seiler Max value
              rc(i,k,n,j) = 1.0 / rr
            END IF
          END DO
        END DO

        ! Shrub dec : H2 dry dep velocity has no dependence on soil moisture
        ! Shrub eg  : H2 dry dep velocity has no dependence on soil moisture
        ! Shrub and bare soil are assumed to be tundra if latitude > 60N;

        DO n = shrub_dec , shrub_eg !  DO n = 12, 13
          DO i = 1, row_length
            IF (todo(i,k,n) .AND. microb(i,k)) THEN
              IF (sinlat(i,k) > tundra_s_limit) THEN
                ! treat shrub as shrub_dec above tundra limit
                rr = h2dd_c(shrub_dec) * 1.0e-4
                rc(i,k,n,j) = 1.0 / rr
              ELSE IF (sinlat(i,k) <= tundra_s_limit) THEN
                ! treat shrub as c3_grass below tundra limit
                sm = MAX(smr(i,k),0.1)
                rr = ( h2dd_m(c3_grass)*sm + h2dd_c(c3_grass) ) * 1.0e-4
                rc(i,k,n,j) = 1.0 / rr
              END IF
            END IF
          END DO
        END DO

        !  Bare soil :
        ! Shrub and bare soil are assumed to be tundra if latitude > 60N;

        n = soil ! n = 16
        DO i = 1, row_length
          IF (todo(i,k,n) .AND. microb(i,k)) THEN
            IF (sinlat(i,k) > tundra_s_limit) THEN
              ! treat soil as shrub_dec above tundra limit
              rr = h2dd_c(shrub_dec) * 1.0e-4
              rc(i,k,soil,j) = 1.0 / rr
            ELSE IF (sinlat(i,k) <= tundra_s_limit) THEN
              ! treat soil as c3_grass below tundra limit
              sm = MAX(smr(i,k),0.1)
              rr = ( h2dd_m(c3_grass)*sm + h2dd_c(c3_grass) ) * 1.0e-4
              rc(i,k,n,j) = 1.0 / rr
            END IF
          END IF
        END DO

      END DO ! DO k = 1, rows

    ELSE
      ! ntype must equal 9, 13, 17, 27
      WRITE(cmessage,'(A)')                                                    &
         'UKCA does not handle flexible tiles yet: ' //                        &
         'Dry deposition needs extending. '          //                        &
         'Please use standard tile configuration'
      errcode=1
      CALL ereport(RoutineName,errcode,cmessage)
    END IF

  ELSE IF (speci(nldepd(j)) == 'CH4       ') THEN

    !        CH4: Calculate an uptake flux initially. The uptake
    !        flux depends on soil moisture, based on results of
    !        Reay et al. (2001). Do the PFTs Broadleaf, Needleleaf,
    !        C3 and C4 grasses first.

    DO k = 1, rows
      DO i = 1, row_length
        IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
          ! initialise to minimum uptake flux
          DO n=1,ntype
            rc(i,k,n,j) = 1.0/r_null ! uptake rate, NOT a resistance!
          END DO
        ELSE
          DO n=1,ntype
            rc(i,k,n,j) = r_null
          END DO
        END IF
        IF (microb(i,k)) THEN
          sm = smr(i,k)
          IF (sm < smfrac1) THEN
            f_ch4_uptk = sm / smfrac1
          ELSE IF (sm > smfrac2) THEN
            f_ch4_uptk = (smfrac3 - sm) / smfrac4
          ELSE
            f_ch4_uptk = 1.0
          END IF
          IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
            f_ch4_uptk = MIN(MAX(f_ch4_uptk,0.0),1.0)
          ELSE
            f_ch4_uptk = MAX(f_ch4_uptk,0.0)
          END IF
          ! When l_fix_ukca_h2dd_x is retired, this could become a CASE
          IF ( .NOT. l_fix_ukca_h2dd_x ) THEN
            IF (todo(i,k,1)) rc(i,k,1,j) = ch4_up_flux(1) * f_ch4_uptk
            IF (todo(i,k,2)) rc(i,k,2,j) = ch4_up_flux(2) * f_ch4_uptk
            IF (todo(i,k,3)) rc(i,k,3,j) = ch4_up_flux(3) * f_ch4_uptk
            IF (todo(i,k,4)) rc(i,k,4,j) = ch4_up_flux(4) * f_ch4_uptk
          ELSE IF ( ntype == 9 ) THEN
            ! brd_leaf
            ! ndl_leaf
            ! c3_grass
            ! c4_grass
            DO n = brd_leaf, c4_grass ! DO n = 1,4
              IF (todo(i,k,n)) rc(i,k,n,j) = ch4_up_flux(n) * f_ch4_uptk
            END DO
          ELSE IF ( ntype == 13 ) THEN
            ! brd_leaf_dec, brd_leaf_eg_trop, brd_leaf_eg_temp
            ! ndl_leaf_dec, ndl_leaf_eg
            ! c3_grass
            ! c4_grass
            DO n = brd_leaf_dec, c4_grass ! DO n = 1,7
              IF (todo(i,k,n)) rc(i,k,n,j) = ch4_up_flux(n) * f_ch4_uptk
            END DO
          ELSE IF ( (ntype == 17) .OR. (ntype == 27) ) THEN
             ! brd_leaf_dec, brd_leaf_eg_trop, brd_leaf_eg_temp
             ! ndl_leaf_dec, ndl_leaf_eg
             ! c3_grass, c3_crop, c3_pasture
             ! c4_grass, c4_crop, c4_pasture
            DO n = brd_leaf_dec, c4_pasture ! DO n = 1,11
              IF (todo(i,k,n)) rc(i,k,n,j) = ch4_up_flux(n) * f_ch4_uptk
            END DO
          ELSE
            ! ntype must equal 9, 13, 17, 27
            WRITE(cmessage,'(A)')                                              &
                 'UKCA does not handle flexible tiles yet: ' //                &
                 'Dry deposition needs extending. '          //                &
                 'Please use standard tile configuration'
            errcode=1
            CALL ereport(RoutineName,errcode,cmessage)
          END IF
        END IF   ! IF (microb(i,k))
      END DO
    END DO

    ! Now do shrub and bare soil, assumed to be tundra if latitude > 60N;
    !  Otherwise, calculate an uptake flux initially.
    !  The uptake flux depends on soil moisture,
    !  based on results of Reay et al. (2001).
    ! When l_fix_ukca_h2dd_x is retired, this could become a CASE statement
    IF ( ( .NOT. l_fix_ukca_h2dd_x ) .OR. ( ntype == 9 ) ) THEN

      IF ( l_fix_ukca_h2dd_x ) THEN
        n = shrub ! n = 5
      ELSE
        n = npft  ! n = 5 (9 tile)
      END IF

      DO k = 1, rows
        DO i = 1, row_length
          IF (microb(i,k)) THEN
            IF (sinlat(i,k) > tundra_s_limit) THEN
              ts = t0tile(i,k,n)
              rr =         ch4dd_tun(fourth_term) +                            &
                   (ts * ( ch4dd_tun(third_term)  +                            &
                   (ts * ( ch4dd_tun(second_term) +                            &
                   (ts *   ch4dd_tun(first_term)) )) ))
              rr = rr * 3600.0 ! Convert from s-1 to h-1
              IF (todo(i,k,n)) rc(i,k,n,j) = MAX(rr,0.0)
            ELSE
              sm = smr(i,k)
              IF (sm < smfrac1) THEN
                f_ch4_uptk = sm / smfrac1
              ELSE IF (sm > smfrac2) THEN
                f_ch4_uptk = (smfrac3 - sm) / smfrac4
              ELSE
                f_ch4_uptk = 1.0
              END IF
              IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
                f_ch4_uptk = MIN(MAX(f_ch4_uptk,0.0),1.0)
              ELSE
                f_ch4_uptk = MAX(f_ch4_uptk,0.0)
              END IF
              IF (todo(i,k,n)) rc(i,k,n,j) = ch4_up_flux(n) * f_ch4_uptk
            END IF
          ELSE
            IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
              rc(i,k,n,j) = 1.0/r_null ! uptake rate, NOT a resistance!
            ELSE
              rc(i,k,n,j) = r_null
            END IF
          END IF
        END DO
      END DO

    ELSE IF ( (ntype == 13) .OR. (ntype == 17) .OR. (ntype == 27) ) THEN

      DO n = shrub_dec, shrub_eg ! DO n=8,9 (ntype==13) ! DO n=12,13 otherwise
        DO k = 1, rows
          DO i = 1, row_length
            IF (microb(i,k)) THEN
              IF (sinlat(i,k) > tundra_s_limit) THEN
                ts = t0tile(i,k,n)
                rr =         ch4dd_tun(fourth_term) +                          &
                     (ts * ( ch4dd_tun(third_term)  +                          &
                     (ts * ( ch4dd_tun(second_term) +                          &
                     (ts *   ch4dd_tun(first_term)) )) ))
                rr = rr * 3600.0 ! Convert from s-1 to h-1
                IF (todo(i,k,n)) rc(i,k,n,j) = MAX(rr,0.0)
              ELSE
                sm = smr(i,k)
                IF (sm < smfrac1) THEN
                  f_ch4_uptk = sm / smfrac1
                ELSE IF (sm > smfrac2) THEN
                  f_ch4_uptk = (smfrac3 - sm) / smfrac4
                ELSE
                  f_ch4_uptk = 1.0
                END IF
                IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
                  f_ch4_uptk = MIN(MAX(f_ch4_uptk,0.0),1.0)
                ELSE
                  f_ch4_uptk = MAX(f_ch4_uptk,0.0)
                END IF
                IF (todo(i,k,n)) rc(i,k,n,j) = ch4_up_flux(n) * f_ch4_uptk
              END IF
            ELSE
              IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
                rc(i,k,n,j) = 1.0/r_null ! uptake rate, NOT a resistance!
              ELSE
                rc(i,k,n,j) = r_null
              END IF
            END IF
          END DO
        END DO
      END DO

    ELSE
      ! ntype must equal 9, 13, 17, 27
      WRITE(cmessage,'(A)')                                                    &
           'UKCA does not handle flexible tiles yet: ' //                      &
           'Dry deposition needs extending. '          //                      &
           'Please use standard tile configuration'
      errcode=1
      CALL ereport(RoutineName,errcode,cmessage)
    END IF

    ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
    IF ( l_fix_ukca_h2dd_x ) THEN
      n = soil
      DO k = 1, rows
        DO i = 1, row_length
          IF (microb(i,k)) THEN
            IF (sinlat(i,k) > tundra_s_limit) THEN
              ts = t0tile(i,k,soil)
              rr =         ch4dd_tun(fourth_term) +                            &
                   (ts * ( ch4dd_tun(third_term)  +                            &
                   (ts * ( ch4dd_tun(second_term) +                            &
                   (ts *   ch4dd_tun(first_term)) )) ))
              rr = rr * 3600.0 ! Convert from s-1 to h-1
              IF (todo(i,k,soil)) rc(i,k,soil,j) = MAX(rr,0.0)
            ELSE
              sm = smr(i,k)
              IF (sm < smfrac1) THEN
                f_ch4_uptk = sm / smfrac1
              ELSE IF (sm > smfrac2) THEN
                f_ch4_uptk = (smfrac3 - sm) / smfrac4
              ELSE
                f_ch4_uptk = 1.0
              END IF
              IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
                f_ch4_uptk = MIN(MAX(f_ch4_uptk,0.0),1.0)
              ELSE
                f_ch4_uptk = MAX(f_ch4_uptk,0.0)
              END IF
              IF (todo(i,k,soil)) rc(i,k,soil,j) = ch4_up_flux(soil) *         &
                                                   f_ch4_uptk
            END IF
          ELSE
            IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
              rc(i,k,n,j) = 1.0/r_null ! uptake rate, NOT a resistance!
            ELSE
              rc(i,k,n,j) = r_null
            END IF
          END IF
        END DO
      END DO
    ELSE
      n = npft ! ( shrub not soil )
      DO k = 1, rows
        DO i = 1, row_length
          IF (microb(i,k)) THEN
            IF (sinlat(i,k) > tundra_s_limit) THEN
              ts = t0tile(i,k,npft) ! ( shrub not soil )
              rr =         ch4dd_tun(fourth_term) +                            &
                   (ts * ( ch4dd_tun(third_term)  +                            &
                   (ts * ( ch4dd_tun(second_term) +                            &
                   (ts *   ch4dd_tun(first_term)) )) ))
              rr = rr * 3600.0 ! Convert from s-1 to h-1
              IF (todo(i,k,soil)) rc(i,k,soil,j) = MAX(rr,0.0)
            ELSE
              sm = smr(i,k)
              IF (sm < smfrac1) THEN
                f_ch4_uptk = sm / smfrac1
              ELSE IF (sm > smfrac2) THEN
                f_ch4_uptk = (smfrac3 - sm) / smfrac4
              ELSE
                f_ch4_uptk = 1.0
              END IF
              IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
                f_ch4_uptk = MIN(MAX(f_ch4_uptk,0.0),1.0)
              ELSE
                f_ch4_uptk = MAX(f_ch4_uptk,0.0)
              END IF
              IF (todo(i,k,soil)) rc(i,k,soil,j) = ch4_up_flux(soil) *         &
                                                   f_ch4_uptk
            END IF
          ELSE
            IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
              rc(i,k,n,j) = 1.0/r_null ! uptake rate, NOT a resistance!
            ELSE
              rc(i,k,n,j) = r_null
            END IF
          END IF
        END DO
      END DO
    END IF

    !         Convert CH4 uptake fluxes (ug m-2 h-1) to
    !         resistance (s m-1).

    IF (ukca_config%l_ukca_emsdrvn_ch4) THEN

      DO n = 1, ntype
        DO k = 1, rows
          DO i = 1, row_length
            IF (todo(i,k,n)) THEN ! global sfc fraction >0.0 for sfc type n
              IF (microb(i,k)) THEN ! there is microbial activity
                rr = rc(i,k,n,j)
                IF (rr > 0.0) THEN ! the default uptake flux is non-zero
                  IF (n == soil) THEN ! apply scaling factor for soil tiles
                    rc(i,k,n,j) = (p0(i,k)*mml)/(rmol*t0tile(i,k,n)*           &
                                   rr*TAR_scaling)
                  ELSE IF ((n == urban) .OR. (n == lake) .OR.                  &
                           (n >= ice)) THEN ! not in cities, lakes or on ice
                    rc(i,k,n,j) = r_null
                  ELSE
                    rc(i,k,n,j) = (p0(i,k)*mml)/(rmol*t0tile(i,k,n)*rr)
                  END IF
                ELSE ! default deposition resistance
                  rc(i,k,n,j) = r_null
                END IF
              ELSE ! no microbial activity in sfc fraction
                rc(i,k,n,j) = r_null
              END IF
              ! process non-PFT and non-soil tiles
              IF (n == urban) rc(i,k,n,j) = r_null
              IF (n == lake) rc(i,k,n,j) = r_null
              IF (n >= ice) rc(i,k,n,j) = r_null
            ELSE ! tile fraction below minimum for sfc type n
              rc(i,k,n,j) = r_null
            END IF
          END DO
        END DO
      END DO

    ELSE

      DO k = 1, rows
        DO n = 1, npft
          DO i = 1, row_length
            IF (todo(i,k,n) .AND. microb(i,k)) THEN
              rr = rc(i,k,n,j)
              IF (rr > 0.0) THEN
                rc(i,k,n,j) = p0(i,k) * mml / (rmol * t0tile(i,k,n) * rr)
              ELSE
                rc(i,k,n,j) = r_null
              END IF
            END IF
          END DO
        END DO

        n = soil
        DO i = 1, row_length
          IF (todo(i,k,n) .AND. microb(i,k)) THEN
            rr = rc(i,k,n,j)
            IF (rr > 0.0) THEN
              rc(i,k,n,j) = p0(i,k)*mml / (rmol*t0tile(i,k,n)*rr*TAR_scaling)
            ELSE
              rc(i,k,n,j) = r_null
            END IF
          END IF
        END DO
      END DO

    END IF

  ELSE IF (speci(nldepd(j)) == 'CO        ') THEN

    !         Only assign values for CO if microbes are active.
    ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
    ! When l_fix_ukca_h2dd_x is retired, this could become a CASE statement
    IF ( ( .NOT. l_fix_ukca_h2dd_x ) .OR. ( ntype == 9 ) ) THEN
      IF ( l_fix_ukca_h2dd_x ) THEN
        n = shrub
      ELSE
        n = npft
      END IF

      DO k = 1, rows
        DO i = 1, row_length
          IF (sinlat(i,k) > tundra_s_limit .AND. microb(i,k)) THEN
            IF (todo(i,k,n))     rc(i,k,n,j)    = r_tundra(r_tundra_co) ! CO
          END IF
        END DO
      END DO
    ELSE IF ( (ntype == 13) .OR. (ntype == 17) .OR. (ntype == 27) ) THEN

      DO n = shrub_dec, shrub_eg
        DO k = 1, rows
          DO i = 1, row_length
            IF (sinlat(i,k) > tundra_s_limit .AND. microb(i,k)) THEN
              IF (todo(i,k,n))     rc(i,k,n,j)    = r_tundra(r_tundra_co) ! CO
            END IF
          END DO
        END DO
      END DO
    ELSE
      ! ntype must equal 9, 13, 17, 27
      WRITE(cmessage,'(A)')                                                    &
           'UKCA does not handle flexible tiles yet: ' //                      &
           'Dry deposition needs extending. '          //                      &
           'Please use standard tile configuration'
      errcode=1
      CALL ereport(RoutineName,errcode,cmessage)
    END IF

    ! Shrub and bare soil are assumed to be tundra if latitude > 60N;
    DO k = 1, rows
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit .AND. microb(i,k)) THEN
          IF (todo(i,k,soil))  rc(i,k,soil,j) = r_tundra(r_tundra_co) ! CO
        END IF
      END DO
    END DO

    !       HONO2

    !       Calculate resistances for HONO2 deposition to ice, which
    !       depend on temperature. Ensure resistance for HONO2 does not fall
    !       below 10 s m-1.

  ELSE IF ((speci(nldepd(j)) == 'HNO3      ' .OR.                              &
            speci(nldepd(j)) == 'HONO2     ' .OR.                              &
            speci(nldepd(j)) == 'ISON      ') .OR.                             &
          ((speci(nldepd(j)) == 'HCl       ' .OR.                              &
            speci(nldepd(j)) == 'HOCl      ' .OR.                              &
            speci(nldepd(j)) == 'HBr       ' .OR.                              &
            speci(nldepd(j))== 'HOBr      ') .AND. (l_fix_improve_drydep))) THEN

    IF ( ( .NOT. l_fix_ukca_h2dd_x )                              .OR.         &
         ( (ntype == 9) .OR. (ntype == 13) .OR. (ntype == 17) ) ) THEN

      IF ( l_fix_ukca_h2dd_x ) THEN
        n = ice
      ELSE
        n = ntype
      END IF
      DO k = 1, rows
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN

            !               Limit temperature to a minimum of 252K. Curve used
            !               only used data between 255K and 273K.

            ts = MAX(t0tile(i,k,n), 252.0)
            rr =         hno3dd_ice(third_term)  +                             &
                 (ts * ( hno3dd_ice(second_term) +                             &
                 (ts *   hno3dd_ice(first_term)) ))
            rc(i,k,n,j) = MAX(rr,10.0)
          END IF
        END DO
      END DO
    ELSE IF ( ntype == 27 ) THEN
      DO n = ice, ice + n_elev_ice ! DO n = 17, 27
        DO k = 1, rows
          DO i = 1, row_length
            IF (todo(i,k,n)) THEN

              !               Limit temperature to a minimum of 252K. Curve used
              !               only used data between 255K and 273K.

              ts = MAX(t0tile(i,k,n), 252.0)
              rr =         hno3dd_ice(third_term)  +                           &
                   (ts * ( hno3dd_ice(second_term) +                           &
                   (ts *   hno3dd_ice(first_term)) ))
              rc(i,k,n,j) = MAX(rr,10.0)
            END IF
          END DO
        END DO
      END DO
    ELSE
      ! ntype must equal 9, 13, 17, 27
      WRITE(cmessage,'(A)')                                                    &
           'UKCA does not handle flexible tiles yet: ' //                      &
           'Dry deposition needs extending. '          //                      &
           'Please use standard tile configuration'
      errcode=1
      CALL ereport(RoutineName,errcode,cmessage)
    END IF

    !       ORGNIT is treated as an aerosol.
    !       Assume vd is valid for all land types and aerosol types.

  ELSE IF ((speci(nldepd(j)) == 'ORGNIT    '  .OR.                             &
            speci(nldepd(j)) == 'BSOA      '  .OR.                             &
            speci(nldepd(j)) == 'ASOA      '  .OR.                             &
            speci(nldepd(j)) == 'ISOSOA    ') .OR.                             &
           (speci(nldepd(j)) =='Sec_Org   ' .AND. (l_fix_improve_drydep)) .OR. &
           (speci(nldepd(j)) =='SEC_ORG_I ' .AND. (l_fix_improve_drydep))) THEN

    DO n = 1, urban     ! for all the functional plant types as
                        ! well as for surface type urban (n=npft+1)
      DO k = 1, rows
        DO i = 1, row_length
          IF (so4_vd(i,k) > 0.0 .AND. gsf(i,k,n) > 0.0) THEN
            rr          = 1.0 / so4_vd(i,k)
            rc(i,k,n,j) = rr
          END IF
        END DO
      END DO
    END DO

    n = soil            !for surface type soil (n=8)
    DO k = 1, rows
      DO i = 1, row_length
        IF (so4_vd(i,k) > 0.0 .AND. gsf(i,k,n) > 0.0) THEN
          rr          = 1.0 / so4_vd(i,k)
          rc(i,k,n,j) = rr
        END IF
      END DO
    END DO

  END IF ! End of IF (speci == species name)

END DO  ! End of DO j = 1, ndepd
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_surfddr

END MODULE ukca_surfddr_mod
