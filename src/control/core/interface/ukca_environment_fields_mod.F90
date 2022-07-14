! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module declaring UKCA's environmental driver fields. These are input
!   fields on the current UKCA model grid that may be varied by the parent
!   application during the run.
!
!   The following subroutines are provided for UKCA internal use
!     locate_land_points      - Determines number and location of land points
!                               in the land-sea mask
!     clear_land_only_fields  - Clears fields defined on land points only
!
!   The module also provides public arrays 'environ_field_info' and
!   'l_environ_field_available' for holding information about each required
!   field and its availability status respectively.
!
!   Public variables 'land_points' and 'land_index' are also provided, giving
!   the number of land points and their locations in the land-sea mask
!   respectively. These are defined when the land-sea mask environmental
!   driver is set and are required by UKCA for unpacking environment fields
!   which are defined at land points only.
!
!   Finally, a public array parameter 'drep' gives representative
!   particle diameters to be associated with input dust emissions in bins
!   'k1_dust_flux' to 'k2_dust_flux' and a public parameter 'ntphot' gives
!   the number of time points represented in photolysis data for the 2D
!   photolysis scheme.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds,
! University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_environment_fields_mod

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

! Public procedures
PUBLIC locate_land_points, clear_land_only_fields

! Public variables are declared below for all recognised environment fields.
! The notation [CF] against field names given in the comments below indicates
! that the preceeding name is a CF-compliant standard name for the field.
! The references {UM:nnnnn} give the STASH code of the equivalent field in the
! Unified Model (in which UKCA was originally developed).
! Some fields are required only by the UM routine for boundary layer mixing
! passed to 'ukca_setup' as 'bl_tracer_mix'. This callback routine is used for
! mixing of tracers after applying emissions when UKCA is coupled with UM
! physics. The fields are marked 'for bl_tracer_mix' and have arbitrary units
! with respect to UKCA.

! --- Scalar values of type real ---

REAL, PARAMETER, PUBLIC :: no_data_value = -999.0  ! Value when not set
REAL, SAVE, PUBLIC :: sin_declination = no_data_value
  ! sin of solar declination
REAL, SAVE, PUBLIC :: equation_of_time = no_data_value
  ! equation of time, specified as an hour angle (radians)
REAL, SAVE, PUBLIC :: atmospheric_ch4 = no_data_value
  ! mass_fraction_of_methane_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_co2 = no_data_value
  ! mass_fraction_of_carbon_dioxide_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_h2 = no_data_value
  ! mass_fraction_of_molecular_hydrogen_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_n2 = no_data_value
  ! mass fraction of molecular nitrogen in air
REAL, SAVE, PUBLIC :: atmospheric_o2 = no_data_value
  ! mass fraction of molecular oxygen in air
REAL, SAVE, PUBLIC :: atmospheric_n2o = no_data_value
  ! mass_fraction_of_nitrous_oxide_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_cfc11 = no_data_value
  ! mass_fraction_of_cfc11_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_cfc12 = no_data_value
  ! mass_fraction_of_cfc12_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_cfc113 = no_data_value
  ! mass_fraction_of_cfc113_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_hcfc22 = no_data_value
  ! mass_fraction_of_hcfc22_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_hfc125 = no_data_value
  ! mass_fraction_of_hfc125_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_hfc134a = no_data_value
  ! mass_fraction_of_hfc134a_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_mebr = no_data_value
  ! mass_fraction_of_methyl_bromide_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_mecl = no_data_value
  ! mass_fraction_of_methyl_chloride_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_ch2br2 = no_data_value
  ! mass fraction of dibromomethane in air
REAL, SAVE, PUBLIC :: atmospheric_chbr3 = no_data_value
  ! mass fraction of tribromomethane in air
REAL, SAVE, PUBLIC :: atmospheric_cfc114 = no_data_value
  ! mass_fraction_of_cfc114_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_cfc115 = no_data_value
  ! mass_fraction_of_cfc115_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_ccl4 = no_data_value
  ! mass_fraction_of_carbon_tetrachloride_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_meccl3 = no_data_value
  ! mass_fraction of methyl chloroform in air
REAL, SAVE, PUBLIC :: atmospheric_hcfc141b = no_data_value
  ! mass_fraction_of_hcfc141b_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_hcfc142b = no_data_value
  ! mass_fraction_of_hcfc142b_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_h1211 = no_data_value
  ! mass_fraction_of_halon1211_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_h1202 = no_data_value
  ! mass_fraction_of_halon1202_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_h1301 = no_data_value
  ! mass_fraction_of_halon1301_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_h2402 = no_data_value
  ! mass_fraction_of_halon2402_in_air [CF]
REAL, SAVE, PUBLIC :: atmospheric_cos = no_data_value
  ! mass fraction of carbonyl sulfide in air

! --- 1D fields of type real ---
REAL, ALLOCATABLE, SAVE, PUBLIC :: soil_moisture_layer1(:)
  ! moisture_content_of_soil_layer [CF] at surface (kg/m^2) {UM:00009}
REAL, ALLOCATABLE, SAVE, PUBLIC :: fland(:)
  ! land_area_fraction [CF] {UM:00505}
REAL, ALLOCATABLE, SAVE, PUBLIC :: ibvoc_isoprene(:)
  ! Gridbox mean isoprene emission flux (kgC/m2/s)
REAL, ALLOCATABLE, SAVE, PUBLIC :: ibvoc_terpene(:)
  ! Gridbox mean (mono-)terpene emission flux (kgC/m2/s)
REAL, ALLOCATABLE, SAVE, PUBLIC :: ibvoc_methanol(:)
  ! Gridbox mean methanol emission flux (kgC/m2/s)
REAL, ALLOCATABLE, SAVE, PUBLIC :: ibvoc_acetone(:)
  ! Gridbox mean acetone emission flux (kgC/m2/s)
REAL, ALLOCATABLE, SAVE, PUBLIC :: inferno_bc(:)
  ! Gridbox mean Black Carbon emissions from fires (kg/m2/s)
REAL, ALLOCATABLE, SAVE, PUBLIC :: inferno_ch4(:)
  ! Gridbox mean Methane emissions from fires (kg/m2/s)
REAL, ALLOCATABLE, SAVE, PUBLIC :: inferno_co(:)
  ! Gridbox mean Carbon monoxide emissions from fires (kg/m2/s)
REAL, ALLOCATABLE, SAVE, PUBLIC :: inferno_nox(:)
  ! Gridbox mean Nitric oxide emissions from fires (kg/m2/s)
REAL, ALLOCATABLE, SAVE, PUBLIC :: inferno_oc(:)
  ! Gridbox mean Organic Carbon emissions from fires (kg/m2/s)
REAL, ALLOCATABLE, SAVE, PUBLIC :: inferno_so2(:)
  ! Gridbox mean Sulphur dioxide emissions from fires (kg/m2/s)

! --- 2D fields of type real ---
REAL, ALLOCATABLE, SAVE, PUBLIC :: latitude(:,:)
  ! latitude [CF] (degrees north)
REAL, ALLOCATABLE, SAVE, PUBLIC :: longitude(:,:)
  ! longitude [CF] (degrees east, >=0, <360)
REAL, ALLOCATABLE, SAVE, PUBLIC :: sin_latitude(:,:)
  ! SIN(latitude [CF])
REAL, ALLOCATABLE, SAVE, PUBLIC :: cos_latitude(:,:)
  ! COS(latitude [CF])
REAL, ALLOCATABLE, SAVE, PUBLIC :: tan_latitude(:,:)
  ! TAN(latitude [CF])
REAL, ALLOCATABLE, SAVE, PUBLIC :: conv_cloud_lwp(:,:)
  ! convective cloud liquid water path (kg/m^2) {UM:00016}
REAL, ALLOCATABLE, SAVE, PUBLIC :: tstar(:,:)
  ! surface_temperature [CF] (K) {UM:00024}
REAL, ALLOCATABLE, SAVE, PUBLIC :: zbl(:,:)
  ! atmosphere_boundary_layer_thickness [CF] (m) {UM:00025}
REAL, ALLOCATABLE, SAVE, PUBLIC :: rough_length(:,:)
  ! surface_roughness_length [CF] (m)
REAL, ALLOCATABLE, SAVE, PUBLIC :: seaice_frac(:,:)
  ! sea_ice_area_fraction [CF] with respect to sea area {UM:00031}
REAL, ALLOCATABLE, SAVE, PUBLIC :: frac_types(:,:)
  ! fraction of surface type with respect to land area {UM:00216}
REAL, ALLOCATABLE, SAVE, PUBLIC :: laift_lp(:,:)
  ! leaf_area_index [CF] of plant functional type {UM:00217}
REAL, ALLOCATABLE, SAVE, PUBLIC :: canhtft_lp(:,:)
  ! canopy_height [CF] of plant functional type (m) {UM:00218}
REAL, ALLOCATABLE, SAVE, PUBLIC :: tstar_tile(:,:)
  ! surface_temperature [CF] on tiles (K) {UM:00233}
REAL, ALLOCATABLE, SAVE, PUBLIC :: z0tile_lp(:,:)
  ! surface_roughness_length [CF] on tiles (m) {UM:00234}
REAL, ALLOCATABLE, SAVE, PUBLIC :: pstar(:,:)
  ! surface_air_pressure [CF] (Pa) {UM:00409}
REAL, ALLOCATABLE, SAVE, PUBLIC :: surf_albedo(:,:)
  ! surface_albedo [CF]
REAL, ALLOCATABLE, SAVE, PUBLIC :: zhsc(:,:)
  ! height of top of decoupled stratocumulus layer, for bl_tracer_mix (L units)
  ! {UM:03073}
REAL, ALLOCATABLE, SAVE, PUBLIC :: u_scalar_10m(:,:)
  ! wind_speed [CF] at 10 m (m/s) {UM:03230}
REAL, ALLOCATABLE, SAVE, PUBLIC :: surf_hf(:,:)
  ! surface_upward_sensible_heat_flux [CF] (W/m^2) {UM:03217}
REAL, ALLOCATABLE, SAVE, PUBLIC :: u_s(:,:)
  ! explicit friction velocity (m/s) {UM:03465}
REAL, ALLOCATABLE, SAVE, PUBLIC :: ch4_wetl_emiss(:,:)
  ! CH4 wetland flux (ugC/m2/s) {UM:08242}
REAL, ALLOCATABLE, SAVE, PUBLIC :: dms_sea_conc(:,:)
  ! mole_concentration_of_dimethyl_sulphide_in_sea_water [CF] (nmol/l)
  ! {UM:00132}
REAL, ALLOCATABLE, SAVE, PUBLIC :: chloro_sea(:,:)
  ! ocean near-surface chlorophyll (kg/m^3) {UM:00096}
REAL, ALLOCATABLE, SAVE, PUBLIC :: dust_flux(:,:,:)
  ! dust emissions by CLASSIC size bin (kg/m^2/s) - combines multiple 2D fields
  ! {UM:03401-03406}
REAL, ALLOCATABLE, SAVE, PUBLIC :: surf_wetness(:,:)
  ! surface wetness {UM:00634}
REAL, ALLOCATABLE, SAVE, PUBLIC :: grid_surf_area(:,:)
  ! Gridbox surface area (m^2)
REAL, ALLOCATABLE, SAVE, PUBLIC :: ext_cg_flash(:,:)
  ! External cloud-to-ground flash rate [s-1]
REAL, ALLOCATABLE, SAVE, PUBLIC :: ext_ic_flash(:,:)
  ! External intracloud flash rate [s-1]

! --- 2D fields of type integer ---
INTEGER, ALLOCATABLE, SAVE, PUBLIC :: kent(:,:)
  ! grid level of surface mixed layer inversion, for bl_tracer_mix {UM:03065}
INTEGER, ALLOCATABLE, SAVE, PUBLIC :: kent_dsc(:,:)
  ! grid level of decoupled stratocumulus inversion, for bl_tracer_mix
  ! {UM:03069}
INTEGER, ALLOCATABLE, SAVE, PUBLIC :: conv_cloud_base(:,:)
  ! lowest convective cloud base level no. {UM:05218}
INTEGER, ALLOCATABLE, SAVE, PUBLIC :: conv_cloud_top(:,:)
  ! lowest convective cloud top level no. {UM:05219}
INTEGER, ALLOCATABLE, SAVE, PUBLIC :: lscat_zhang(:,:)
  ! surface type for gridbox as required by Zhang et al (Aero drydep) scheme

! --- 2D fields of type logical ---
LOGICAL, ALLOCATABLE, SAVE, PUBLIC :: land_sea_mask(:,:)
  ! land_binary_mask [CF] (land = .TRUE.) {UM:00030}
LOGICAL, ALLOCATABLE, SAVE, PUBLIC :: l_tile_active(:,:)
  ! active tile indicator (.TRUE. if tile is in use)

! --- 3D fields of type real ---
REAL, ALLOCATABLE, SAVE, PUBLIC :: theta(:,:,:)
  ! air_potential_temperature [CF] (K) {UM:00004}
REAL, ALLOCATABLE, SAVE, PUBLIC :: q(:,:,:)
  ! specific_humidity [CF] {UM:00010}
REAL, ALLOCATABLE, SAVE, PUBLIC :: qcf(:,:,:)
  ! mass_fraction_of_cloud_ice_in_air [CF] {UM:00012}
REAL, ALLOCATABLE, SAVE, PUBLIC :: conv_cloud_amount(:,:,:)
  ! fractional convective cloud amount with anvil {UM:00211}
REAL, ALLOCATABLE, SAVE, PUBLIC :: rho_r2(:,:,:)
  ! density*radius*radius (kg/m) {UM:00253}
REAL, ALLOCATABLE, SAVE, PUBLIC :: qcl(:,:,:)
  ! mass_fraction_of_cloud_liquid_water_in_air [CF] {UM:00254}
REAL, ALLOCATABLE, SAVE, PUBLIC :: exner_rho_levels(:,:,:)
  ! dimensionless_exner_function [CF] at rho levels {UM:00255}
REAL, ALLOCATABLE, SAVE, PUBLIC :: area_cloud_fraction(:,:,:)
  ! cloud_area_fraction_in_atmosphere_layer [CF] {UM:00265}
REAL, ALLOCATABLE, SAVE, PUBLIC :: cloud_frac(:,:,:)
  ! bulk cloud fraction {UM:00266}
REAL, ALLOCATABLE, SAVE, PUBLIC :: cloud_liq_frac(:,:,:)
  ! liquid cloud fraction {UM:00267}
REAL, ALLOCATABLE, SAVE, PUBLIC :: exner_theta_levels(:,:,:)
  ! dimensionless_exner_function [CF] at theta levels {UM:00406}
REAL, ALLOCATABLE, SAVE, PUBLIC :: p_rho_levels(:,:,:)
  ! air_pressure [CF] at rho levels (Pa) {UM:00407}
REAL, ALLOCATABLE, SAVE, PUBLIC :: p_theta_levels(:,:,:)
  ! air_pressure [CF] at theta levels (Pa) {UM:00408}
REAL, ALLOCATABLE, SAVE, PUBLIC :: rhokh_rdz(:,:,:)
  ! Mixing coefficient above surface:
  ! (scalar eddy diffusivity * density) / dz, for bl_tracer_mix
  ! (M/(L^2*T) units ) {UM:03060}
REAL, ALLOCATABLE, SAVE, PUBLIC :: dtrdz(:,:,:)
  ! dt/(density*radius*radius*dz) for scalar flux divergence, for bl_tracer_mix
  ! (T/M units) {UM:03064}
REAL, ALLOCATABLE, SAVE, PUBLIC :: we_lim(:,:,:)
  ! density * entrainment rate implied by placing of subsidence at surface mixed
  ! layer inversion, for bl_tracer_mix
  ! (M/L^2/T units) {UM:03066}
REAL, ALLOCATABLE, SAVE, PUBLIC :: t_frac(:,:,:)
  ! fraction of timestep surface mixed layer inversion is above level,
  ! for bl_tracer mix {UM:03067}
REAL, ALLOCATABLE, SAVE, PUBLIC :: zrzi(:,:,:)
  ! level height as fraction of surface mixed layer inversion height above ML
  ! base, for bl_tracer_mix {UM:03068}
REAL, ALLOCATABLE, SAVE, PUBLIC :: we_lim_dsc(:,:,:)
  ! density * entrainment rate implied by placing of subsidence at decoupled
  ! stratocumulus inversion, for bl_tracer_mix (M/L^2/T units) {UM:03070}
REAL, ALLOCATABLE, SAVE, PUBLIC :: t_frac_dsc(:,:,:)
  ! fraction of timestep decoupled stratocumulus inversion is above level,
  ! for bl_tracer_mix {UM:03071}
REAL, ALLOCATABLE, SAVE, PUBLIC :: zrzi_dsc(:,:,:)
  ! level height as fraction of decoupled stratocumulus inversion height above
  ! DSC ML base, for bl_tracer_mix
  ! {UM:03072}
REAL, ALLOCATABLE, SAVE, PUBLIC :: stcon(:,:,:)
  ! stomatal conductance on plant functional type (m/s) {UM:03462}
REAL, ALLOCATABLE, SAVE, PUBLIC :: ls_rain3d(:,:,:)
  ! large_scale_rainfall_flux [CF] out of model levels (kg/m^2/s) {UM:04222}
REAL, ALLOCATABLE, SAVE, PUBLIC :: ls_snow3d(:,:,:)
  ! large_scale_snowfall_flux [CF] out of model levels (kg/m^2/s) {UM:04223}
REAL, ALLOCATABLE, SAVE, PUBLIC :: autoconv(:,:,:)
  ! rain autoconversion rate, mass fraction (/s) {UM:04257}
REAL, ALLOCATABLE, SAVE, PUBLIC :: accretion(:,:,:)
  ! rain accretion rate, mass fraction (/s) {UM:04258}
REAL, ALLOCATABLE, SAVE, PUBLIC :: pv_on_theta_mlevs(:,:,:)
  ! ertel_potential_vorticity [CF] on theta levels (K m^2/(kg s) {UM:15218}
REAL, ALLOCATABLE, SAVE, PUBLIC :: conv_rain3d(:,:,:)
  ! convective_rainfall_flux [CF] (kg/m^2/s) {UM:05227}
REAL, ALLOCATABLE, SAVE, PUBLIC :: conv_snow3d(:,:,:)
  ! convective_snowfall_flux [CF] (kg/m^2/s) {UM:05228}
REAL, ALLOCATABLE, SAVE, PUBLIC :: so4_sa_clim(:,:,:)
  ! sulphate aerosol surface area density (cm^2/cm^3)
REAL, ALLOCATABLE, SAVE, PUBLIC :: so4_aitken(:,:,:)
  ! mass fraction of aitken mode sulfate dry aerosol expressed as sulfur in air
  ! {UM:00103}
REAL, ALLOCATABLE, SAVE, PUBLIC :: so4_accum(:,:,:)
  ! mass fraction of accumulation mode sulfate dry aerosol expressed as sulfur
  ! in air {UM:00104}
REAL, ALLOCATABLE, SAVE, PUBLIC :: soot_fresh(:,:,:)
  ! mass fraction of fresh black carbon dry aerosol in air {UM:00108}
REAL, ALLOCATABLE, SAVE, PUBLIC :: soot_aged(:,:,:)
  ! mass fraction of aged black carbon dry aerosol in air {UM:00109}
REAL, ALLOCATABLE, SAVE, PUBLIC :: ocff_fresh(:,:,:)
  ! mass fraction of fresh organic carbon from fossil fuel dry aerosol in air
  ! {UM:00114}
REAL, ALLOCATABLE, SAVE, PUBLIC :: ocff_aged(:,:,:)
  ! mass fraction of aged organic carbon from fossil fuel dry aerosol in air
  ! (kg/kg) {UM:00115}
REAL, ALLOCATABLE, SAVE, PUBLIC :: biogenic(:,:,:)
  ! mass fraction of biogenic secondary organic dry aerosol in air {UM:00351}
REAL, ALLOCATABLE, SAVE, PUBLIC :: dust_div1(:,:,:)
  ! mass fraction of dust division 1 dry aerosol in air {UM:00431}
REAL, ALLOCATABLE, SAVE, PUBLIC :: dust_div2(:,:,:)
  ! mass fraction of dust division 2 dry aerosol in air {UM:00432}
REAL, ALLOCATABLE, SAVE, PUBLIC :: dust_div3(:,:,:)
  ! mass fraction of dust division 3 dry aerosol in air {UM:00433}
REAL, ALLOCATABLE, SAVE, PUBLIC :: dust_div4(:,:,:)
  ! mass fraction of dust division 4 dry aerosol in air {UM:00434}
REAL, ALLOCATABLE, SAVE, PUBLIC :: dust_div5(:,:,:)
  ! mass fraction of dust division 5 dry aerosol in air {UM:00435}
REAL, ALLOCATABLE, SAVE, PUBLIC :: dust_div6(:,:,:)
  ! mass fraction of dust division 6 dry aerosol in air {UM:00436}
REAL, ALLOCATABLE, SAVE, PUBLIC :: sea_salt_film(:,:,:)
  ! atmosphere number concentration of film mode sea salt particles (/m^3)
  ! {UM:01247}
REAL, ALLOCATABLE, SAVE, PUBLIC :: sea_salt_jet(:,:,:)
  ! atmosphere number concentration of jet mode sea salt particles (/m^3)
  ! {UM:01248}
REAL, ALLOCATABLE, SAVE, PUBLIC :: co2_interactive(:,:,:)
  ! mass_fraction_of_carbon_dioxide_in_air [CF] {UM:00252}
REAL, ALLOCATABLE, SAVE, PUBLIC :: rim_cry(:,:,:)
  ! riming rate for ice crystals, mass fraction (/s)
  ! from_cloud_liquid [CF] (/s) {UM:04247}
REAL, ALLOCATABLE, SAVE, PUBLIC :: rim_agg(:,:,:)
  ! riming rate for ice aggregates, mass fraction (/s) {UM:04248}
REAL, ALLOCATABLE, SAVE, PUBLIC :: vertvel(:,:,:)
  ! upward_air_velocity [CF] (m/s) {UM:00150}
REAL, ALLOCATABLE, SAVE, PUBLIC :: bl_tke(:,:,:)
  ! turbulent kinetic energy (m2/s2) {UM:03473}
REAL, ALLOCATABLE, SAVE, PUBLIC :: interf_z(:,:,:)
  ! Altitude of grid-cell interfaces (m)

! --- Oxidant species for Offline Oxidants chemistry - 3-D of type Real ---
! Note: These fields can also be used for other chemistry schemes if needed
! e.g. o3_offline is used for overwriting stratospheric values in tropospheric
! chemistry schemes.
REAL, ALLOCATABLE, SAVE, PUBLIC :: h2o2_offline(:,:,:)
  ! Hydrogen peroxide as offline oxidant (kg/kg-air)
REAL, ALLOCATABLE, SAVE, PUBLIC :: ho2_offline(:,:,:)
  ! Hydroperoxyl radical as offline oxidant (kg/kg-air)
REAL, ALLOCATABLE, SAVE, PUBLIC :: no3_offline(:,:,:)
  ! Nitrate radical as offline oxidant (kg/kg-air)
REAL, ALLOCATABLE, SAVE, PUBLIC :: o3_offline(:,:,:)
  ! mass_fraction_of_ozone_in_air [CF] {UM:00060}
REAL, ALLOCATABLE, SAVE, PUBLIC :: oh_offline(:,:,:)
  ! Hydroxyl radical as offline oxidant (kg/kg-air)

! --- 4D fields of type real ---
REAL, ALLOCATABLE, SAVE, PUBLIC :: photol_rates(:,:,:,:)
  ! Interpolated photolysis rates
  ! Dimensions are rows, row_length, model_levels,
  ! and number of photolytic reactions (jppj)

! Public variables for unpacking land-only fields
INTEGER, SAVE, PUBLIC :: land_points = 0            ! Number of land points
INTEGER, ALLOCATABLE, SAVE, PUBLIC :: land_index(:) ! Location of land points
                                                    ! on the 2D grid

! Type used to hold information associated with an environment field
TYPE, PUBLIC :: env_field_info_type
  INTEGER :: group        ! Group that field belongs to
  INTEGER :: lbound_dim1  ! Lower bound of dimension 1
  INTEGER :: ubound_dim1  ! Upper bound of dimension 1
  INTEGER :: lbound_dim2  ! Lower bound of dimension 2
  INTEGER :: ubound_dim2  ! Upper bound of dimension 2
  INTEGER :: lbound_dim3  ! Lower bound of dimension 3
  INTEGER :: ubound_dim3  ! Upper bound of dimension 3
  INTEGER :: lbound_dim4  ! Lower bound of dimension 4
  INTEGER :: ubound_dim4  ! Upper bound of dimension 4
  LOGICAL :: l_land_only  ! True if field is defined on land points only
END TYPE env_field_info_type

! Missing data value for environment field array bound (for field bounds
! that cannot be fully defined in 'init_environment_req')
INTEGER, PARAMETER, PUBLIC :: no_bound_value = -999

! Bounds for 'dust_flux' array dimension 3. Elements correspond to size bins
! with boundaries at 0.0316, 0.1, 0.316, 1.0, 3.16, 10.0 and 31.6 um radius
INTEGER, PARAMETER, PUBLIC :: k1_dust_flux = 1
INTEGER, PARAMETER, PUBLIC :: k2_dust_flux = 6
! Representative particle diameters (m)
REAL, PARAMETER, PUBLIC :: drep(k1_dust_flux:k2_dust_flux) =                   &
  [0.112468e-06, 0.355656e-06, 0.112468e-05,                                   &
   0.355656e-05, 0.112468e-04, 0.355656e-04]

! Number of time points in 2D photolysis data
INTEGER, PARAMETER, PUBLIC :: ntphot = 3

! Meta-data corresponding to the list of environment fields required for the
! current UKCA configuration defined by 'environ_field_varnames' in
! module 'ukca_environment_req_mod'
LOGICAL, ALLOCATABLE, TARGET, SAVE, PUBLIC ::                                  &
  l_environ_field_available(:)     ! Availability status flag
TYPE(env_field_info_type), ALLOCATABLE, SAVE, PUBLIC ::                        &
  environ_field_info(:)            ! Field info

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_ENVIRONMENT_FIELDS_MOD'


CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE locate_land_points()
! ----------------------------------------------------------------------
! Description:
!   Set up index of land points for locating land-only environment fields
!   on a 2-D grid and update the bounds information for all such fields
!   to reflect the expected number of points.
!
! Method:
!   For each element set to true in the 2D field 'land_sea_mask', set an
!   index value giving the position of the element in an equivalent 1D
!   array and assign 'land_points' and 'land_index' to the number of
!   these indices and their values respectively.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Local variables

INTEGER :: i_land(SIZE(land_sea_mask))
INTEGER :: row_length
INTEGER :: rows
INTEGER :: i
INTEGER :: j

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'LOCATE_LAND_POINTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Determine land point indices
row_length = SIZE(land_sea_mask, DIM=1)
rows = SIZE(land_sea_mask, DIM=2)
land_points = 0
DO j =1, rows
  DO i = 1, row_length
    IF (land_sea_mask(i,j)) THEN
      land_points = land_points + 1
      i_land(land_points) = (j-1)*row_length + i
    END IF
  END DO
END DO

! Create array of land point indices from temporary array
IF (ALLOCATED(land_index)) DEALLOCATE(land_index)
ALLOCATE(land_index(land_points))
land_index = i_land(1:land_points)

! Update field info for land-only fields
DO i = 1, SIZE(environ_field_info)
  IF (environ_field_info(i)%l_land_only)                                       &
    environ_field_info(i)%ubound_dim1 = land_points
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE locate_land_points


! ----------------------------------------------------------------------
SUBROUTINE clear_land_only_fields()
! ----------------------------------------------------------------------
! Description:
!   Deallocates environment field arrays for fields defined on land
!   points only and updates their availability flags
! ----------------------------------------------------------------------

IMPLICIT NONE

! Local variables

INTEGER :: i

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CLEAR_LAND_ONLY_FIELDS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(soil_moisture_layer1)) DEALLOCATE(soil_moisture_layer1)
IF (ALLOCATED(fland)) DEALLOCATE(fland)
IF (ALLOCATED(l_tile_active)) DEALLOCATE(l_tile_active)
IF (ALLOCATED(frac_types)) DEALLOCATE(frac_types)
IF (ALLOCATED(laift_lp)) DEALLOCATE(laift_lp)
IF (ALLOCATED(canhtft_lp)) DEALLOCATE(canhtft_lp)
IF (ALLOCATED(tstar_tile)) DEALLOCATE(tstar_tile)
IF (ALLOCATED(z0tile_lp)) DEALLOCATE(z0tile_lp)
IF (ALLOCATED(inferno_so2)) DEALLOCATE(inferno_so2)
IF (ALLOCATED(inferno_oc)) DEALLOCATE(inferno_oc)
IF (ALLOCATED(inferno_nox)) DEALLOCATE(inferno_nox)
IF (ALLOCATED(inferno_co)) DEALLOCATE(inferno_co)
IF (ALLOCATED(inferno_ch4)) DEALLOCATE(inferno_ch4)
IF (ALLOCATED(inferno_bc)) DEALLOCATE(inferno_bc)
IF (ALLOCATED(ibvoc_acetone)) DEALLOCATE(ibvoc_acetone)
IF (ALLOCATED(ibvoc_methanol)) DEALLOCATE(ibvoc_methanol)
IF (ALLOCATED(ibvoc_terpene)) DEALLOCATE(ibvoc_terpene)
IF (ALLOCATED(ibvoc_isoprene)) DEALLOCATE(ibvoc_isoprene)

DO i = 1, SIZE(environ_field_info)
  IF (environ_field_info(i)%l_land_only) l_environ_field_available(i) = .FALSE.
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE clear_land_only_fields

END MODULE ukca_environment_fields_mod
