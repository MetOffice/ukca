! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module providing data structures to specify details of the active
!   UKCA model configuration:
!     ukca_config   - contains the UKCA configuration variables, excluding
!                     variables specific to GLOMAP-MODE
!     glomap_config - contains the config. variables specific to GLOMAP-MODE
!
!   The module also provides the following procedure for the UKCA_API
!     ukca_get_config - returns values for UKCA configuration variables.
!
!   The following additional public procedures are provided for use within UKCA
!     init_ukca_configuration  - initialises/resets all UKCA configuration
!                                data ready for a new UKCA configuration to
!                                be set up
!     copy_config_vector       - returns a copy of a vector-valued
!                                configuration variable
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_config_specification_mod

USE missing_data_mod, ONLY: imdi, rmdi

USE ukca_mode_setup,  ONLY: glomap_variables_type

IMPLICIT NONE
PUBLIC

! ---------------------------------------------------------------------------
! -- Type for holding UKCA configuration data (excludes GLOMAP-specifics) --
! ---------------------------------------------------------------------------

! Components in each section other than internal variables are set by the
! parent application via optional keyword arguments with matching names
! in 'ukca_setup'.
! In the event of any change to keyword argument names, old components may
! be retained in the structure as duplicates of the new input variables
! (to support retention of old keyword argument names for backward
! compatibility and/or to avoid the need to change variable names throughout
! the UKCA code), but should be re-categorised as internal variables.

TYPE :: ukca_config_spec_type

  ! -- Context information --
  ! Includes array sizes for input fields, indices for accessing specific
  ! items in the input arrays, model timestep etc
  INTEGER :: row_length                ! X dimension of UKCA domain (columns)
  INTEGER :: rows                      ! Y dimension of UKCA domain (rows)
  INTEGER :: model_levels              ! Z dimension of UKCA domain (levels)
  INTEGER :: bl_levels                 ! Number of boundary layer levels
  INTEGER :: nlev_ent_tr_mix           ! No. of grid levels for entrainment-
                                       ! related fields used by tr_mix
  INTEGER :: ntype                     ! Number of surface types considered in
                                       ! interactive dry deposition
  INTEGER :: npft                      ! Number of plant functional types
                                       ! considered in interactive dry
                                       ! deposition
  INTEGER :: i_brd_leaf                ! Index of surface type
                                       ! 'broad-leaf tree'
  INTEGER :: i_brd_leaf_dec            ! Index of surface type
                                       ! 'broad-leaf dec tree'
  INTEGER :: i_brd_leaf_eg_trop        ! Index of surface type
                                       ! 'broad-leaf eg trop tree'
  INTEGER :: i_brd_leaf_eg_temp        ! Index of surface type
                                       ! 'broad-leaf eg temp tree'
  INTEGER :: i_ndl_leaf                ! Index of surface type
                                       ! 'needle-leaf tree'
  INTEGER :: i_ndl_leaf_dec            ! Index of surface type
                                       ! 'needle-leaf dec tree'
  INTEGER :: i_ndl_leaf_eg             ! Index of surface type
                                       ! 'needle-leaf eg tree'
  INTEGER :: i_c3_grass                ! Index of surface type 'c3 grass'
  INTEGER :: i_c3_crop                 ! Index of surface type 'c3 crop'
  INTEGER :: i_c3_pasture              ! Index of surface type 'c3 pasture'
  INTEGER :: i_c4_grass                ! Index of surface type 'c4 grass'
  INTEGER :: i_c4_crop                 ! Index of surface type 'c4 crop'
  INTEGER :: i_c4_pasture              ! Index of surface type 'c4 pasture'
  INTEGER :: i_shrub                   ! Index of surface type 'shrub'
  INTEGER :: i_shrub_dec               ! Index of surface type 'shrub dec'
  INTEGER :: i_shrub_eg                ! Index of surface type 'shrub eg'
  INTEGER :: i_urban                   ! Index of surface type 'urban'
  INTEGER :: i_lake                    ! Index of surface type 'lake'
  INTEGER :: i_soil                    ! Index of surface type 'soil'
  INTEGER :: i_ice                     ! Index of surface type 'ice'
  INTEGER, ALLOCATABLE :: i_elev_ice(:)
                                       ! Indices of surface type 'elevated ice'
  REAL :: dzsoil_layer1                ! Thickness of surface soil layer (m)
  LOGICAL :: l_cal360                  ! True if UKCA is expected to use a
                                       ! 360-day calendar for applying a top
                                       ! boundary condition in stratospheric
                                       ! chemistry schemes and/or for
                                       ! validating external emissions.
                                       ! (Emissions registered with
                                       ! a day-of-week scaling requirement
                                       ! can't be used with a 360-day
                                       ! calendar).
                                       !!!! Also required for external
                                       !!!! photolysis code pending complete
                                       !!!! separation of photolysis from UKCA.
  REAL :: timestep                     ! Model timestep in seconds (an integer
                                       ! no of timesteps is expected in 1 hour)

  ! -- General UKCA configuration options --
  INTEGER :: i_ukca_chem               ! Chemistry scheme to use
  LOGICAL :: l_ukca_chem_aero          ! True for aerosol precursor chemistry
  LOGICAL :: l_ukca_mode               ! True for GLOMAP-mode aerosol scheme
  LOGICAL :: l_fix_tropopause_level    ! True to use an arbitrary fixed level
                                       ! for defining the tropopause height
                                       ! rather than diagnosing from PV and
                                       ! theta (Tropopause height affects
                                       ! heterogeneous polar stratospheric
                                       ! cloud chemistry, stratospheric updates
                                       ! for troposphere only schemes,
                                       ! interaction between nitrate and dust
                                       ! schemes and production of various
                                       ! diagnostics)
  INTEGER :: fixed_tropopause_level    ! Level to define a fixed tropopause
  LOGICAL :: l_ukca_ageair             ! True for age-of-air scheme
  INTEGER :: i_ageair_reset_method     ! Method for controlling age reset to 0
  INTEGER :: max_ageair_reset_level    ! Max. level for reset (Method 1)
  REAL :: max_ageair_reset_height      ! Max. height (m) for reset (Method 2)
  LOGICAL :: l_blankout_invalid_diags  ! True to set all non-valid diagnostic
                                       ! field data to missing data value
  LOGICAL :: l_enable_diag_um          ! True to enable diagnostic output via
                                       ! the Unified Model STASH system
  LOGICAL :: l_ukca_persist_off        ! True for no saving of horizontal
                                       ! arrays between timesteps
  LOGICAL :: l_timer                   ! True to use a timer routine for timing
                                       ! sub-models
  LOGICAL :: l_ukca_emissions_off      ! True to turn off emissions (and lower
                                       ! boundary conditions for stratospheric
                                       ! chemistry schemes)
  LOGICAL :: l_ukca_drydep_off         ! True to turn off dry deposition
  LOGICAL :: l_ukca_wetdep_off         ! True to turn off wet deposition
  INTEGER :: i_error_method            ! Error handling method to use
  LOGICAL :: l_ukca_scale_ppe          ! True to turn on PPE scaling parameters

  ! -- Chemistry configuration options --
  INTEGER :: i_ukca_chem_version       ! Chemical mechanism version identifier
  INTEGER :: nrsteps                   ! No. of N-R solver iterations before
                                       ! time step is halved
  INTEGER :: chem_timestep             ! Chemical timestep in seconds for N-R
                                       ! and offline oxidant schemes
  INTEGER :: dts0                      ! Backward Euler timestep (seconds)
  LOGICAL :: l_ukca_asad_columns       ! True to pass columns to ASAD solver
                                       ! rather than horizontal slices
  LOGICAL :: l_ukca_asad_full          ! True to pass the entire domain to ASAD
                                       ! solver rather than slices or columns
  LOGICAL :: l_ukca_debug_asad         ! Include additional print output
                                       ! specific to ASAD
  LOGICAL :: l_ukca_intdd              ! True for interactive dry deposition
  LOGICAL :: l_ukca_ddepo3_ocean       ! True to use Luhar et al. (2018)
                                       ! oceanic O3 dry-deposition scheme
  LOGICAL :: l_ukca_ddep_lev1          ! True to apply dry deposition losses
                                       ! only in level 1
  LOGICAL :: l_ukca_dry_dep_so2wet     ! True if considering the impact of
                                       ! surface wetness on dry deposition
  LOGICAL :: l_deposition_jules        ! True when using JULES-based interactive
                                       ! dry deposition routines
  INTEGER :: nit                       ! Number of iterations for B-E Solver
  LOGICAL :: l_ukca_quasinewton        ! True to use quasi-Newton (Broyden)
                                       ! method to reduce no. of iterations
                                       ! in N-R solver step
  INTEGER :: i_ukca_quasinewton_start  ! Iteration to start quasi-Newton step
                                       ! (2-50, 2 recommended)
  INTEGER :: i_ukca_quasinewton_end    ! Iteration to stop quasi-Newton step
                                       ! (2-50, 3 recommended)
  INTEGER :: ukca_chem_seg_size        ! Grid points per segment for the
                                       ! column-based chemical solver
  REAL :: max_z_for_offline_chem       ! Maximum height at which to integrate
                                       ! chemistry with the explicit B-E
                                       ! Offline Oxidants scheme
  INTEGER :: nlev_above_trop_o3_env    ! Number of levels above tropopause at
                                       ! at which to start overwriting O3
                                       ! and HNO3 with values derived from O3
                                       ! data in troposphere only schemes
  INTEGER :: nlev_ch4_stratloss        ! Number of top levels at which to
                                       ! apply stratospheric CH4 loss rate in
                                       ! troposhere only schemes
  LOGICAL :: l_tracer_lumping          ! True to handle lumped tracers in
                                       ! stratospheric chemistry schemes
  INTEGER :: i_ukca_topboundary        ! Method of treating top boundary for
                                       ! certain species in stratospheric
                                       ! schemes
  LOGICAL :: l_ukca_ro2_ntp            ! True to stop transport of peroxy
                                       ! radicals in StratTrop/CRI scheme
  LOGICAL :: l_ukca_ro2_perm           ! True for RO2-permutation chemistry in
                                       ! StratTrop/CRI scheme
  LOGICAL :: l_ukca_intph              ! True for calculating interactive cloud
                                       ! pH values (F for global value of 5)
  REAL :: ph_fit_coeff_a               ! value for fitting coefficient a of
                                       ! the interactive cloud pH code
  REAL :: ph_fit_coeff_b               ! value for fitting coefficient b of
                                       ! the interactive cloud pH code
  REAL :: ph_fit_intercept             ! y-intercept value of the relationship
                                       ! for the interactive cloud pH code
  LOGICAL :: l_ukca_scale_soa_yield_mt ! True to apply scaling factor to
                                       ! production of Secondary Organic
                                       ! Aerosol (SOA) from monoterpene
  REAL :: soa_yield_scaling_mt         ! Scaling factor for production of SOA
                                       ! from monoterpene
  LOGICAL :: l_ukca_scale_soa_yield_isop ! True to apply scaling factor to
                                       ! production of Secondary Organic
                                       ! Aerosol (SOA) from isoprene
  REAL :: soa_yield_scaling_isop       ! Scaling factor for production of SOA
                                       ! from isoprene
  REAL :: dry_depvel_so2_scaling       ! SO2 dry deposition velocity scaling
                                       ! factor

  ! -- Chemistry - Heterogeneous chemistry --

  LOGICAL :: l_ukca_het_psc            ! True for heterogeneous polar
                                       ! stratospheric cloud chemistry
  INTEGER :: i_ukca_hetconfig          ! 0 = default (5 reactions); 1 = JPL-15
                                       ! recommended coeff. 2 = JPL-15 + bromine
                                       ! reactions
  LOGICAL :: l_ukca_limit_nat          ! True to limit the formation of NAT in
                                       ! polar stratospheric clouds (type 1) to
                                       ! heights above 1 km
  LOGICAL :: l_ukca_sa_clim            ! True to use surface area density
                                       ! climatology for heterogeneous PSC chem.
  LOGICAL :: l_ukca_trophet            ! True for tropospheric het. chemistry
  LOGICAL :: l_ukca_classic_hetchem    ! True if accounting for heterogeneous
                                       ! chemistry on CLASSIC aerosols (RAQ
                                       ! chemistry scheme only)

  LOGICAL :: l_use_photolysis          ! True if UKCA is expecting to apply
                                       ! Photolysis rates in the solver and
                                       ! hence needs the rates from the parent
  !!!! Currently this will be set based on i_ukca_photol in the namelist !!!

  ! -- UKCA emissions configuration options --
  LOGICAL :: l_ukca_ibvoc              ! True for interactive bVOC emissions
  LOGICAL :: l_ukca_inferno            ! True for INFERNO fire emissions
  LOGICAL :: l_ukca_inferno_ch4        ! True for INFERNO CH4 fire emissions
  INTEGER :: i_inferno_emi             ! maximum INFERNO emission level
  LOGICAL :: l_ukca_so2ems_expvolc     ! True for interactive expvolc emissions
  LOGICAL :: l_ukca_so2ems_plumeria    ! True for expvolc emissions with plumeria
  LOGICAL :: l_ukca_qch4inter          ! True for interactive wetland CH4
                                       ! emissions
  LOGICAL :: l_ukca_emsdrvn_ch4        ! True when running UKCA in
                                       ! CH4 emissions-driven mode
  REAL :: mode_parfrac                 ! Fraction of SO2 emissions as aerosol(%)
  LOGICAL :: l_ukca_enable_seadms_ems  ! True to explicitly enable marine DMS
                                       ! emissions
                                       ! (deprecated redundant setting; if
                                       ! omitted from 'ukca_setup' call,
                                       ! DMS emissions are enabled or
                                       ! disabled according to the value of
                                       ! 'i_ukca_dms_flux')
  INTEGER :: i_ukca_dms_flux           ! Sea-air DMS exchange scheme to use
  LOGICAL :: l_ukca_scale_seadms_ems   ! True to apply scaling to marine DMS
                                       ! emissions
  REAL :: seadms_ems_scaling           ! Marine DMS emission scaling factor
  LOGICAL :: l_ukca_linox_scaling      ! True to use LOG(p) to distribute
                                       ! lightning NOx in the vertical
  REAL :: lightnox_scale_fac           ! Lightning NOx emission scale factor
  INTEGER :: i_ukca_light_param        ! Choice for flash-freq parameterisation
                                       ! (1 = Price & Rind, 2 = Luhar et al.,
                                       !  3 = external lightning scheme)
  LOGICAL :: l_support_ems_vertprof    ! True to support full range of vertical
                                       ! scaling options for 2D emission fields.
                                       ! (Required for SNAP sector, 'high_level'
                                       ! and 'single_level' emissions).
  LOGICAL :: l_support_ems_gridbox_units
                                       ! True to support the provision of
                                       ! offline emissions in grid-box units
  LOGICAL :: l_suppress_ems            ! True to suppress all tracer updates
                                       ! from emissions and boundary layer
                                       ! mixing (Emissions diagnostics
                                       ! are unaffected)
  REAL :: anth_so2_ems_scaling         ! Anthropogenic SO2 emissions scaling
                                       ! factor

  ! -- UKCA feedback configuration options --
  LOGICAL :: l_ukca_h2o_feedback       ! True for H2O feedback from chemistry
  LOGICAL :: l_ukca_conserve_h         ! True to include hydrogen conservation
                                       ! when 'l_ukca_h2o_feedback' is true

  ! -- UKCA environmental driver configuration options --
  LOGICAL :: l_param_conv              ! True if using precipitation diagnostics
                                       ! from a parameterized convection scheme
  LOGICAL :: l_ctile                   ! True if using partial land points at
                                       ! coasts
  LOGICAL :: l_zon_av_ozone            ! True if external O3, used for
                                       ! troposphere only schemes, is expected
                                       ! as a mean value for each domain row
                                       ! (i.e. zonal mean if lat const. in row)
  INTEGER :: i_strat_lbc_source        ! Source for gas MMR values for lower
                                       ! boundary conditions in stratospheric
                                       ! chemistry schemes
  LOGICAL :: l_chem_environ_gas_scalars
                                       ! True if using external values for
                                       ! CO2, CH4, O2, H2 and N2 in chemistry
  LOGICAL :: l_chem_environ_co2_fld    ! True if using an external CO2 field
                                       ! in chemistry
  LOGICAL :: l_ukca_prescribech4       ! True if prescribing surface CH4
  LOGICAL :: l_use_classic_so4         ! True to use CLASSIC SO4 for het. chem.
  LOGICAL :: l_use_classic_soot        ! True to use CLASSIC black carbon for
                                       ! het. chem. (RAQ scheme only)
  LOGICAL :: l_use_classic_ocff        ! True to use CLASSIC organic carbon from
                                       ! fossil fuels for het. chem. (RAQ only)
  LOGICAL :: l_use_classic_biogenic    ! True to use CLASSIC biogenic secondary
                                       ! organic aerosol for het. chem. (RAQ
                                       ! only)
  LOGICAL :: l_use_classic_seasalt     ! True to use CLASSIC sea salt for het.
                                       ! chem. (RAQ only)
  LOGICAL :: l_use_gridbox_volume      ! True to use gridbox volume in
                                       ! diagnostic calculations. (This is
                                       ! required for UM diagnostics and to
                                       ! enable use of ASAD framework
                                       ! diagnostics in non-UM applications.)
  LOGICAL :: l_use_gridbox_mass        ! True to use mass of air in grid box in
                                       ! prognostic and/or diagnostic
                                       ! calculations.
                                       ! (Note that GLOMAP-mode prognostics
                                       ! will run with or without this option
                                       ! giving equivalent but not
                                       ! bit-comparable results. The option is
                                       ! required for the stratospheric
                                       ! chemistry schemes and for UM
                                       ! diagnostics.)
  LOGICAL :: l_environ_z_top           ! True if using an external value for
                                       ! height at top of model (for bit-
                                       ! comparability with previous results
                                       ! when running in the UM)
  INTEGER :: env_log_step              ! Timestep number at which environment
                                       ! summary diagnostics should be printed
                                       ! to log file

  ! -- UKCA temporary logicals --
  LOGICAL :: l_fix_ukca_cloud_frac     ! True to fix cloud_frac offset bug
  LOGICAL :: l_fix_improve_drydep      ! True to fix dry deposition velocities
  LOGICAL :: l_fix_ukca_h2dd_x         ! True to fix H2 deposition to shrub/soil
  LOGICAL :: l_fix_drydep_so2_water    ! True to use correct surface resistance
                                       ! of water when calculating dry
                                       ! deposition of SO2
  LOGICAL :: l_fix_ukca_offox_h2o_fac  ! True to fix water vapour units in
                                       ! B-E Offline Oxidants scheme
  LOGICAL :: l_fix_ukca_h2so4_ystore   ! True to fix storage of H2SO4 in ASAD
                                       ! N-R schemes for updating in GLOMAP

  ! -- UKCA internal configuration variables (not modifiable by parent) --
  LOGICAL :: l_ukca_chem               ! True if chemistry is on
  LOGICAL :: l_ukca_trop               ! True for StdTrop chemistry (B-E)
  LOGICAL :: l_ukca_aerchem            ! True to add aerosol chemistry to
                                       ! StdTrop scheme
  LOGICAL :: l_ukca_raq                ! True for RAQ chemistry (B-E)
  LOGICAL :: l_ukca_raqaero            ! True for RAQ + aerosol chemistry (B-E)
  LOGICAL :: l_ukca_offline_be         ! True for B-E Offline Oxidants chemistry
  LOGICAL :: l_ukca_tropisop           ! True for Trop chem. + Isoprene (N-R)
  LOGICAL :: l_ukca_strattrop          ! True for Strat-trop chemistry (N-R)
  LOGICAL :: l_ukca_strat              ! True for Strat chemistry scheme (N-R)
  LOGICAL :: l_ukca_offline            ! True for N-R Offline Oxidants chemistry
                                       ! scheme
  LOGICAL :: l_ukca_cristrat           ! True for CRI-Strat chemistry (N-R)
  LOGICAL :: l_ukca_stratcfc           ! True for Strat-CFC chemistry (N-R)
  LOGICAL :: l_ukca_achem              ! True to add aerosol chemistry to an
                                       ! N-R scheme
  LOGICAL :: l_ukca_nr_aqchem          ! True when aqueous chem required for N-R
  LOGICAL :: l_ukca_advh2o             ! True if H2O treated as tracer by ASAD
  LOGICAL :: l_asad_chem_diags_support ! True if ASAD diagnostics support is
                                       ! included
  LOGICAL :: l_diurnal_isopems         ! True for applying diurnal cycle to
                                       ! isoprene emissions
  LOGICAL :: l_seawater_dms            ! True for seawater DMS emissions
  LOGICAL :: l_chem_environ_ch4_scalar ! True if external CH4 needed
  LOGICAL :: l_chem_environ_co2_scalar ! True if external CO2 needed
  LOGICAL :: l_chem_environ_h2_scalar  ! True if external H2 needed
  LOGICAL :: l_chem_environ_n2_scalar  ! True if external N2 needed
  LOGICAL :: l_chem_environ_o2_scalar  ! True if external O2 needed
  INTEGER :: ukca_int_method           ! Chemical integration method
  INTEGER :: timesteps_per_day         ! No. of model timesteps in a day
  INTEGER :: timesteps_per_hour        ! No. of model timesteps in an hour

END TYPE ukca_config_spec_type

! ---------------------------------------------------------------------------
! -- Type for holding GLOMAP-specific configuration data --
! ---------------------------------------------------------------------------

! Components in each section other than internal variables are set by the
! parent application via optional keyword arguments with matching names in
! 'ukca_setup'.
! In the event of any change to keyword argument names, old components may
! be retained in the structure as duplicates of the new input variables
! (to support retention of old keyword argument names for backward
! compatibility and/or to avoid the need to change variable names throughout
! the UKCA code), but should be re-categorised as internal variables.

TYPE :: glomap_config_spec_type

  ! -- General GLOMAP configuration options --
  INTEGER :: i_mode_nzts               ! No. of substeps for nucleation/
                                       ! sedimentation
  INTEGER :: ukca_mode_seg_size        ! Columns per segment for aerosol
  INTEGER :: i_mode_setup              ! MODE aerosol scheme to use
  LOGICAL :: l_mode_bhn_on             ! True for binary homogeneous sulphate
                                       ! nucleation
  LOGICAL :: l_mode_bln_on             ! True for boundary layer sulphate
                                       ! nucleation
  INTEGER :: i_mode_bln_param_method   ! Parameterisation method for boundary
                                       ! layer nucleation 1=activation;
                                       ! 2=kinetic; 3=organic mediated (as per
                                       ! Metzer et al, PNAS, 2010) - required
                                       ! when 'l_mode_bhn_on' is true.
  REAL :: mode_activation_dryr         ! Activation dry radius (nm)
  LOGICAL :: l_dust_ageing_on          ! True for dust being subjected
                                       ! to nuc scav, conden, coag,
                                       ! ageing, and activation
  REAL :: dry_depvel_acc_scaling       ! Scaling factor for the dry deposition
                                       ! velocity for the accumulation mode
  REAL :: acc_cor_scav_scaling         ! Scaling factor for scavenging
                                       ! parameters for the accumulation and
                                       ! coarse modes

  ! -- GLOMAP deposition configuration options --
  LOGICAL :: l_ddepaer                 ! True for aerosol dry deposition
  REAL :: mode_incld_so2_rfrac         ! Fraction of in-cloud oxidised SO2
                                       ! removed by precipitation
  LOGICAL :: l_rainout                 ! True for turning on nucleation
                                       ! scavenging (rainout)
  INTEGER :: i_mode_nucscav            ! Choice of nucl. scavenging co-effs:
                                       ! 1=original, 2=ECHAM5-HAM
                                       ! 3=as(1) but no scav of modes 6&7
  LOGICAL :: l_cv_rainout              ! True for convective rainout in UKCA
  LOGICAL :: l_impc_scav               ! True for turning on impaction
                                       ! scavenging
  LOGICAL :: l_dust_slinn_impc_scav    ! True for turning on the Slinn
                                       ! impaction scheme for dust

  ! -- GLOMAP emissions configuration options --
  LOGICAL :: l_ukca_primss             ! True for primary sea-salt emissions
  LOGICAL :: l_ukca_primsu             ! True for primary sulphate emissions
  LOGICAL :: l_ukca_primdu             ! True for primary dust emissions
  LOGICAL :: l_ukca_primbcoc           ! True for primary BC/OC emissions
  LOGICAL :: l_ukca_prim_moc           ! True for primary marine OC aerosol
                                       ! emissions
  LOGICAL :: l_bcoc_bf                 ! True for primary biofuel BC/OC emiss.
  LOGICAL :: l_bcoc_bm                 ! True for primary biomass BC/OC emiss.
  LOGICAL :: l_bcoc_ff                 ! True for primary fossil fuel BC/OC
                                       ! emissions
  LOGICAL :: l_ukca_scale_biom_aer_ems ! True to apply scaling factor to
                                       ! biomass burning BC/OC aerosol emissions
  LOGICAL :: l_ukca_scale_sea_salt_ems ! True to apply scaling factor to
                                       ! sea salt emissions
  LOGICAL :: l_ukca_scale_marine_pom_ems  ! True to apply scaling factor to
                                          ! marine particulate organice matter
  REAL :: biom_aer_ems_scaling         ! Biomass-burning emissions scaling
  LOGICAL :: l_ukca_fine_no3_prod      ! True for fine mode NO3/NH4 emissions
  LOGICAL :: l_ukca_coarse_no3_prod    ! True for coarse mode NO3 emissions
  LOGICAL :: l_no3_prod_in_aero_step   ! True for nitrate emissions in MODE
  REAL :: hno3_uptake_coeff            ! HNO3 uptake coefficient for nitrate
  REAL :: sea_salt_ems_scaling         ! Sea salt emission scaling factor
  REAL :: marine_pom_ems_scaling       ! Marine POM emission scaling factor
  INTEGER :: i_primss_method           ! Sea-salt emission options

  ! -- GLOMAP feedback configuration options --
  LOGICAL :: l_ukca_radaer             ! Provide output for calculating direct
                                       ! radiative effects of aerosols using
                                       ! RADAER
  INTEGER :: i_ukca_tune_bc            ! Options for tuning BC absorption
                                       ! via adjusting BC mass density and
                                       ! refractive index mixing approximations
  INTEGER :: i_ukca_activation_scheme  ! Activation scheme to use
  INTEGER :: i_ukca_nwbins             ! Value of nwbins in Activate (1-20)
                                       ! See Rosalind West paper for details
                                       ! doi:10.5194/acp-14-6369-2014
  REAL :: sigwmin                      ! Lower limit for std. dev. of updraft
                                       ! velocity pdf used in Activate
                                       ! (default: 0.01 m/s).
  LOGICAL :: l_ntpreq_n_activ_sum      ! True to include total number
                                       ! concentration of active aerosol in
                                       ! non-transported prognostics
                                       ! (e.g. for combination with an
                                       ! external cloud mask for calculating
                                       ! CDNC in the UM hybrid resolution model)
  LOGICAL :: l_ntpreq_dryd_nuc_sol     ! True to include to include the dry
                                       ! diameter for nucleation soluble mode.
  LOGICAL :: l_ukca_sfix               ! True for diagnosing UKCA CCN at
                                       ! fixed supersaturation
  REAL :: sigma_updraught_scaling      ! Scaling factor for standard deviation
                                       ! of updraught velocities

  ! -- GLOMAP temporary logicals --
  LOGICAL :: l_fix_neg_pvol_wat        ! True to check, trap and fix for
                                       ! negative values in 'pvol' and 'mdwat'
  LOGICAL :: l_fix_ukca_impscav        ! True to fix impaction scavenging bugs
  LOGICAL :: l_fix_nacl_density        ! True to correct sea-salt density
  LOGICAL :: l_improve_aero_drydep     ! if True the surface type used by the
                                       ! aerosol dry deposition is provided by
                                       ! the parent, rather than being inferred
                                       ! from roughness length in GLOMAP.
  LOGICAL :: l_fix_ukca_activate_pdf   ! True to fix probability density fn.
                                       ! of updraft velocities in Activate
  LOGICAL :: l_fix_ukca_activate_vert_rep
                                       ! True to fix CDNC vertical
                                       ! replication bug in Activate
  LOGICAL :: l_bug_repro_tke_index     ! True to reproduce effects of an old bug
                                       ! in the calculation of updraft pdf
                                       ! spread related to TKE indexing in the
                                       ! Unified Model
  LOGICAL :: l_fix_ukca_hygroscopicities
                                       ! True to update aerosol component
                                       ! hygroscopicities to kappa-Kohler
                                       ! values and fix bug in calculation of
                                       ! volume-weighted hygroscopicity value
  LOGICAL :: l_fix_ukca_water_content
                                       ! True to fix bugs in the calculation
                                       ! of aerosol water content values.

  ! -- GLOMAP internal configuration variables (not modifiable by parent) --
  INTEGER :: n_dust_emissions          ! Number of dust emission size ranges
  INTEGER :: i_dust_scheme             ! Dust scheme for NO3 second. prod.
  LOGICAL :: l_6bin_dust_no3           ! True for CLASSIC 6-bin dust scheme
  LOGICAL :: l_2bin_dust_no3           ! True for CLASSIC 2-bin dust scheme

END TYPE glomap_config_spec_type

! ---------------------------------------------------------------------------
! -- Option values for configuration variables --
! ---------------------------------------------------------------------------

! These option codes form part of the UKCA API and cannot be changed
! without potentially affecting backwards compatibility with UKCA input data

! Option codes for 'i_ukca_chem'
INTEGER, PARAMETER :: i_ukca_chem_off        = 0  ! No chemistry
! Backward Euler schemes
INTEGER, PARAMETER :: i_ukca_chem_trop       = 11 ! Std. Tropospheric scheme
INTEGER, PARAMETER :: i_ukca_chem_raq        = 13 ! Regional Air Quality scheme
INTEGER, PARAMETER :: i_ukca_chem_offline_be = 14 ! Offline Oxidants (B-E)
! Newton Raphson schemes
INTEGER, PARAMETER :: i_ukca_chem_tropisop   = 50 ! Trop-isoprene scheme
INTEGER, PARAMETER :: i_ukca_chem_strattrop  = 51 ! Strat-trop (CheST) scheme
INTEGER, PARAMETER :: i_ukca_chem_strat      = 52 ! Std. Stratospheric scheme
INTEGER, PARAMETER :: i_ukca_chem_offline    = 54 ! Offline Oxidants (N-R)
                                                  ! scheme
INTEGER, PARAMETER :: i_ukca_chem_cristrat   = 59 ! Common Representative
                                                  ! Intermidiates (CRI) scheme
                                                  ! with stratospheric chem.

! Option codes for 'i_mode_setup'
INTEGER, PARAMETER :: i_suss_4mode              = 1
INTEGER, PARAMETER :: i_sussbcoc_5mode          = 2
INTEGER, PARAMETER :: i_sussbcoc_4mode          = 3
INTEGER, PARAMETER :: i_sussbcocso_5mode        = 4
INTEGER, PARAMETER :: i_sussbcocso_4mode        = 5
INTEGER, PARAMETER :: i_du_2mode                = 6
! i_du_3mode = 7 has not been included yet
INTEGER, PARAMETER :: i_sussbcocdu_7mode        = 8
! i_sussbcocdu_4mode = 9 has not been included yet
INTEGER, PARAMETER :: i_sussbcocntnh_5mode_7cpt = 10
INTEGER, PARAMETER :: i_solinsol_6mode          = 11
INTEGER, PARAMETER :: i_sussbcocduntnh_8mode_8cpt = 12

! Option codes for 'i_ageair_reset_method', controlling how the near-surface
! values of the age-of-air tracer are reset to zero
INTEGER, PARAMETER :: i_age_reset_by_level = 1   ! Based on model level number
INTEGER, PARAMETER :: i_age_reset_by_height = 2  ! Based on height above ground

! Option codes for 'i_error_method'
INTEGER, PARAMETER :: i_error_method_abort = 1   ! Write error message and abort
INTEGER, PARAMETER :: i_error_method_return = 2  ! Return control to parent
INTEGER, PARAMETER :: i_error_method_warn_and_return = 3
                                                 ! Return control to parent
                                                 ! after printing error messagea
                                                 ! as a warning

! Option codes for 'i_strat_lbc_source'
INTEGER, PARAMETER :: i_strat_lbc_off = 0   ! LBCs off
INTEGER, PARAMETER :: i_strat_lbc_wmoa1 = 1 ! Use internal WMO A1 values
INTEGER, PARAMETER :: i_strat_lbc_env = 2   ! Use environment values
INTEGER, PARAMETER :: i_strat_lbc_rcp = 3   ! Use RCP file values

! Option codes for 'ukca_int_method'
INTEGER, PARAMETER :: int_method_impact = 1       ! IMPACT solver
INTEGER, PARAMETER :: int_method_nr = 3           ! ASAD Newton-Raphson solver
INTEGER, PARAMETER :: int_method_be = 5           ! ASAD Backward-Euler solver
INTEGER, PARAMETER :: int_method_be_explicit = 10 ! Explicit B-E solver

! Option codes for 'i_ukca_activation_scheme'
INTEGER, PARAMETER :: i_ukca_activation_off = 0   ! No activation scheme
INTEGER, PARAMETER :: i_ukca_activation_arg = 1   ! Abdul-Razzak & Ghan
INTEGER, PARAMETER :: i_ukca_activation_jones = 2 ! Jones

! Option codes for 'i_ukca_light_param'
INTEGER, PARAMETER :: i_light_param_off   = 0     ! No lightning scheme
INTEGER, PARAMETER :: i_light_param_pr    = 1     ! Original Price & Rind
INTEGER, PARAMETER :: i_light_param_luhar = 2     ! Updated Luhar et al.
INTEGER, PARAMETER :: i_light_param_ext   = 3     ! External lightning scheme

! Option codes for 'i_ukca_topboundary'
INTEGER, PARAMETER :: i_top_none    = 0   ! No overwriting
INTEGER, PARAMETER :: i_top_2levH2O = 1   ! Overwrite top 2 levels (except H2O)
INTEGER, PARAMETER :: i_top_1lev    = 2   ! Overwrite top level
INTEGER, PARAMETER :: i_top_BC      = 3   ! Boundary condition for NO, CO, & O3
INTEGER, PARAMETER :: i_top_BC_H2O  = 4   ! Boundary condition for NO, CO, O3
                                          ! and H2O

! Option codes for 'i_ukca_dms_flux'
INTEGER, PARAMETER :: i_dms_flux_off = 0    ! No marine DMS flux
INTEGER, PARAMETER :: i_liss_merlivat = 1   ! Liss & Merlivat (1986)
INTEGER, PARAMETER :: i_wanninkhof = 2      ! Wanninkhof (1992)
INTEGER, PARAMETER :: i_nightingale = 3     ! Nightingale et al. (2000)
INTEGER, PARAMETER :: i_blomquist = 4       ! Blomquist et al. (2017)

! Option codes for 'i_primss_method'
INTEGER, PARAMETER :: i_primss_method_smith = 1    ! Smith 1998
INTEGER, PARAMETER :: i_primss_method_monahan = 2  ! Gong (2003) and
                                                   ! Monahan (1986)
INTEGER, PARAMETER :: i_primss_method_combined = 3 ! Combination of above
INTEGER, PARAMETER :: i_primss_method_jaegle = 4   ! Jaegle (2011)

! -- Data structures specifying details of the active UKCA configuration --
! ---------------------------------------------------------------------------

TYPE(ukca_config_spec_type),   SAVE :: ukca_config
TYPE(glomap_config_spec_type), SAVE :: glomap_config

! ---------------------------------------------------------------------------
! -- Data structures specifying details of the glomap mode configuration --
! ---------------------------------------------------------------------------

TYPE(glomap_variables_type), SAVE, TARGET :: glomap_variables
TYPE(glomap_variables_type), SAVE, TARGET :: glomap_variables_climatology

! ---------------------------------------------------------------------------
! -- Configurable set of constants used in UKCA --
! ---------------------------------------------------------------------------
! These are constants that can potentially be overridden by a parent
! application (e.g. if required for consistency with values used elsewhere).
! Values are set in the routine 'init_ukca_configuration' defined below.
! Override functionality via the 'ukca_setup' argument list is not yet
! available. UM values are currently provided via 'ukca_um_legacy_mod'.
! These replace the default values from this module in UM builds.

!!!! Note that 'pi', 'zerodegc' and 'rmol' are currently declared as parameters
!!!! to avoid the need for adding them to OpenMP parallel sections in
!!!! 'ukca_abdulrazzak_ghan_mod', 'ukca_aero_ctl_mod' and
!!!! 'ukca_calc_drydiam_mod' and consequently having to make them variables
!!!! in the UM version of 'ukca_um_legacy_mod'. They will need to be changed to
!!!! variables and initialised in 'init_ukca_configuration' before they can
!!!! be overriden via 'ukca_setup'.

! Conversion Parameters
INTEGER :: isec_per_day = imdi    ! No. of seconds in a day (24 hours)
INTEGER :: isec_per_hour = imdi   ! No. of seconds in an hour
REAL :: rhour_per_day = rmdi      ! No. of hours in a day
REAL :: rsec_per_day = rmdi       ! No. of seconds in a day
REAL :: rsec_per_hour = rmdi      ! No. of seconds in an hour
REAL, PARAMETER :: pi = 3.14159265358979323846
                                  ! Pi
REAL :: pi_over_180 = rmdi        ! Pi / 180
REAL :: recip_pi_over_180 = rmdi  ! 180 / Pi
REAL, PARAMETER :: zerodegc = 273.15
                                  ! 0 degrees Celsius in Kelvin

! Molar universal gas constant (J/K/mol)
REAL, PARAMETER :: rmol = 8.314

! Planet Constants (Earth values by default)
REAL :: planet_radius = rmdi  ! Planet radius (m)
REAL :: g = rmdi              ! Mean acceleration due to gravity at
                              ! surface (m/s^2)
REAL :: r = rmdi              ! Gas constant for dry air (J/kg/K)
REAL :: cp = rmdi             ! Specific heat of dry air at constant
                              ! pressure (J/kg/K)
REAL :: pref = rmdi           ! Reference surface pressure (Pa)
REAL :: vkman = rmdi          ! Von Karman's constant
REAL :: repsilon = rmdi       ! Ratio of molecular weights of water and dry air
REAL :: rv = rmdi             ! Gas constant for water vapour (J/kg/K)
REAL :: kappa = rmdi          ! r/cp
REAL :: c_virtual = rmdi      ! (1/repsilon) - 1
REAL :: s2r = rmdi            ! Seconds-to-radians converter

! CLASSIC aerosol parameters for heterogeneous PSC chemistry
REAL :: rad_ait = rmdi  ! Mean radius of Aitken mode particles (m)
REAL :: rad_acc = rmdi  ! Meac radius of accumulation mode particles (m)
REAL :: chi = rmdi      ! Mole fraction of S in particle
REAL :: sigma = rmdi    ! Std. dev. of accumulation mode particle size
                        ! distribution

! ---------------------------------------------------------------------------
! -- Templates for parent callback procedures to be used in UKCA --
! ---------------------------------------------------------------------------
! These can be provided via the ukca_setup call to perform parent-specific
! processing.

ABSTRACT INTERFACE

  ! Subroutine to do boundary layer mixing for a tracer after applying
  ! emission(s). This does not have a UKCA default so a parent routine must
  ! be provided if emissions is on with tracer updates enabled.
  SUBROUTINE template_proc_bl_tracer_mix(row_length, rows, bl_levels,          &
                                         r_theta_levels, r_rho_levels,         &
                                         nlev_ent_tr_mix,                      &
                                         kent, kent_dsc, surf_em, zhnl, zhsc,  &
                                         we_lim, t_frac, zrzi,                 &
                                         we_lim_dsc, t_frac_dsc, zrzi_dsc,     &
                                         z_uv, rhokh_rdz, dtrdz, field)
  IMPLICIT NONE
  INTEGER :: row_length
  INTEGER :: rows
  INTEGER :: bl_levels
  INTEGER :: nlev_ent_tr_mix
  REAL, INTENT(IN) :: r_theta_levels(1:row_length,1:rows,0:bl_levels)
  ! Height of theta levels from Earth centre
  REAL, INTENT(IN) :: r_rho_levels(1:row_length,1:rows,bl_levels)
  ! Height of rho levels from Earth centre
  INTEGER, INTENT(IN) :: kent(row_length, rows)
    ! Grid level of surface mixed layer inversion
  INTEGER, INTENT(IN) :: kent_dsc(row_length, rows)
    ! Grid level of decoupled stratocumulus inversion
  REAL, INTENT(IN) :: surf_em(row_length, rows)
    ! Emission flux into surface level (kg/m^2/s)
  REAL, INTENT(IN) :: zhnl(row_length, rows)
    ! Atmosphere_boundary layer thickness (m) {UM stashcode:00025}
  REAL, INTENT(IN) :: zhsc(row_length, rows)
    ! Height of top of decoupled stratocumulus layer (m) {UM stashcode:03073}
  REAL, INTENT(IN) :: we_lim(row_length, rows, nlev_ent_tr_mix)
    ! Density * entrainment rate implied by placing of subsidence at surface
    ! mixed layer inversion (kg/m^2/s) {UM stashcode:03066}
  REAL, INTENT(IN) :: t_frac(row_length, rows, nlev_ent_tr_mix)
    ! Fraction of timestep surface mixed layer inversion is above level
    ! {UM stashcode:03067}
  REAL, INTENT(IN) :: zrzi(row_length, rows, nlev_ent_tr_mix)
    ! Level height as fraction of surface mixed layer inversion height above ML
    ! base {UM stashcode:03068}
  REAL, INTENT(IN) :: we_lim_dsc(row_length, rows, nlev_ent_tr_mix)
    ! Density * entrainment rate implied by placing of subsidence at decoupled
    ! stratocumulus inversion (kg/m^2/s) {UM stashcode:03070}
  REAL, INTENT(IN) :: t_frac_dsc(row_length, rows, nlev_ent_tr_mix)
    ! Fraction of timestep decoupled stratocumulus inversion is above level
    ! {UM stashcode:03071}
  REAL, INTENT(IN) :: zrzi_dsc(row_length, rows, nlev_ent_tr_mix)
    ! Level height as fraction of decoupled stratocumulus inversion height above
    ! DSC ML base {UM stashcode:03072}
  REAL, INTENT(IN) :: z_uv(row_length, rows, bl_levels)
    ! Height at rho levels (m)
  REAL, INTENT(IN) :: rhokh_rdz(row_length, rows, 2:bl_levels)
    ! Mixing coefficient above surface:
    ! (scalar eddy diffusivity * density) / dz (kg/m^2/s) {UM stashcode:03060}
  REAL, INTENT(IN) :: dtrdz(row_length, rows, bl_levels)
    ! Dt/(density*radius*radius*dz) for scalar flux divergence (s/kg)
    ! {UM stashcode:03064}
  REAL, INTENT(IN OUT) :: field(row_length, rows, bl_levels)
    ! Tracer mixing ratio (kg/kg)
  END SUBROUTINE template_proc_bl_tracer_mix

  ! Subroutine to calculate ozone column. This does not have a UKCA default
  ! so a parent routine must be provided.
  SUBROUTINE template_proc_calc_ozonecol(row_length, rows, model_levels,       &
    z_top_of_model, p_layer_boundaries, p_layer_centres, ozone_vmr, ozonecol)
  IMPLICIT NONE
  ! Model dimensions
  INTEGER, INTENT(IN) :: row_length
  INTEGER, INTENT(IN) :: rows
  INTEGER, INTENT(IN) :: model_levels

  REAL, INTENT(IN) :: z_top_of_model       ! model top (m)
  REAL, INTENT(IN) :: p_layer_boundaries(row_length, rows,                     &
    0:model_levels)
  REAL, INTENT(IN) :: p_layer_centres(row_length, rows, model_levels)
  REAL, INTENT(IN) :: ozone_vmr(row_length, rows, model_levels)

  REAL, INTENT(OUT) :: ozonecol(row_length, rows, model_levels)
  END SUBROUTINE template_proc_calc_ozonecol

  ! Subroutine to do parent-specific copy of 2D output for a named diagnostic
  ! (for direct copy to parent workspace other than array argument in API call)
  SUBROUTINE template_proc_diag2d_copy_out(diagname, field)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: diagname     ! Diagnostic name
  REAL, INTENT(IN) :: field(:,:)               ! 2D field for output
  END SUBROUTINE template_proc_diag2d_copy_out

  ! Subroutine to do parent-specific copy of 3D output for a named diagnostic
  ! (for direct copy to parent workspace other than array argument in API call)
  SUBROUTINE template_proc_diag3d_copy_out(diagname, field)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: diagname     ! Diagnostic name
  REAL, INTENT(IN) :: field(:,:,:)             ! 3D field for output
  END SUBROUTINE template_proc_diag3d_copy_out

END INTERFACE

PROCEDURE(template_proc_bl_tracer_mix), POINTER :: bl_tracer_mix
PROCEDURE(template_proc_calc_ozonecol), POINTER :: calc_ozonecol
PROCEDURE(template_proc_diag2d_copy_out), POINTER :: diag2d_copy_out
PROCEDURE(template_proc_diag3d_copy_out), POINTER :: diag3d_copy_out

! ---------------------------------------------------------------------------
! -- Flag to indicate whether a UKCA configuration is set up --
! ---------------------------------------------------------------------------
LOGICAL, SAVE :: l_ukca_config_available = .FALSE.

! ---------------------------------------------------------------------------
! -- Generic interface for vector copy subroutines --
! ---------------------------------------------------------------------------
INTERFACE copy_config_vector
  MODULE PROCEDURE copy_config_vector_integer
END INTERFACE copy_config_vector

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE init_ukca_configuration()
! ----------------------------------------------------------------------
! Description:
!   Initialises/resets the UKCA configuration data ready for a new
!   configuration to be set up.
!   In the UKCA and GLOMAP configuration structures, all logicals are
!   set to .FALSE., arrays are DEALLOCATED and other variables are set
!   to missing data values.
!   In addition, all configurable constants are set to their UKCA
!   default values.
! ----------------------------------------------------------------------

IMPLICIT NONE

l_ukca_config_available = .FALSE.

! -- Context information --
ukca_config%row_length = imdi
ukca_config%rows = imdi
ukca_config%model_levels = imdi
ukca_config%bl_levels = imdi
ukca_config%nlev_ent_tr_mix = imdi
ukca_config%ntype = imdi
ukca_config%npft = imdi
ukca_config%i_brd_leaf = imdi
ukca_config%i_brd_leaf_dec = imdi
ukca_config%i_brd_leaf_eg_trop = imdi
ukca_config%i_brd_leaf_eg_temp = imdi
ukca_config%i_ndl_leaf = imdi
ukca_config%i_ndl_leaf_dec = imdi
ukca_config%i_ndl_leaf_eg = imdi
ukca_config%i_c3_grass = imdi
ukca_config%i_c3_crop = imdi
ukca_config%i_c3_pasture = imdi
ukca_config%i_c4_grass = imdi
ukca_config%i_c4_crop = imdi
ukca_config%i_c4_pasture = imdi
ukca_config%i_shrub = imdi
ukca_config%i_shrub_dec = imdi
ukca_config%i_shrub_eg = imdi
ukca_config%i_urban = imdi
ukca_config%i_lake = imdi
ukca_config%i_soil = imdi
ukca_config%i_ice = imdi
IF (ALLOCATED(ukca_config%i_elev_ice)) DEALLOCATE(ukca_config%i_elev_ice)
ukca_config%dzsoil_layer1 = rmdi
ukca_config%l_cal360 = .FALSE.
ukca_config%timestep = rmdi

! -- General UKCA configuration options --
ukca_config%i_ukca_chem = imdi
ukca_config%l_ukca_chem_aero = .FALSE.
ukca_config%l_ukca_mode = .FALSE.
ukca_config%l_fix_tropopause_level = .FALSE.
ukca_config%fixed_tropopause_level = rmdi
ukca_config%l_ukca_ageair = .FALSE.
ukca_config%i_ageair_reset_method = imdi
ukca_config%max_ageair_reset_level = imdi
ukca_config%max_ageair_reset_height = rmdi
ukca_config%l_blankout_invalid_diags = .FALSE.
ukca_config%l_enable_diag_um = .FALSE.
ukca_config%l_ukca_persist_off = .FALSE.
ukca_config%l_timer = .FALSE.
ukca_config%l_ukca_emissions_off = .FALSE.
ukca_config%l_ukca_drydep_off = .FALSE.
ukca_config%l_ukca_wetdep_off = .FALSE.
ukca_config%i_error_method = imdi
ukca_config%l_ukca_scale_ppe = .FALSE.

! -- Chemistry configuration options --
ukca_config%i_ukca_chem_version = imdi
ukca_config%nrsteps = imdi
ukca_config%chem_timestep = imdi
ukca_config%dts0 = imdi
ukca_config%l_ukca_asad_columns = .FALSE.
ukca_config%l_ukca_asad_full = .FALSE.
ukca_config%l_ukca_debug_asad = .FALSE.
ukca_config%l_ukca_intdd = .FALSE.
ukca_config%l_ukca_ddepo3_ocean = .FALSE.
ukca_config%l_ukca_ddep_lev1 = .FALSE.
ukca_config%l_ukca_dry_dep_so2wet = .FALSE.
ukca_config%l_deposition_jules = .FALSE.
ukca_config%nit = imdi
ukca_config%l_ukca_quasinewton = .FALSE.
ukca_config%i_ukca_quasinewton_start = imdi
ukca_config%i_ukca_quasinewton_end = imdi
ukca_config%ukca_chem_seg_size = imdi
ukca_config%max_z_for_offline_chem = rmdi
ukca_config%nlev_above_trop_o3_env = imdi
ukca_config%nlev_ch4_stratloss = imdi
ukca_config%l_tracer_lumping = .FALSE.
ukca_config%i_ukca_topboundary = imdi
ukca_config%l_ukca_ro2_ntp = .FALSE.
ukca_config%l_ukca_ro2_perm = .FALSE.
ukca_config%l_ukca_intph = .FALSE.
ukca_config%ph_fit_coeff_a = rmdi
ukca_config%ph_fit_coeff_b = rmdi
ukca_config%ph_fit_intercept = rmdi
ukca_config%l_ukca_scale_soa_yield_mt = .FALSE.
ukca_config%soa_yield_scaling_mt = rmdi
ukca_config%l_ukca_scale_soa_yield_isop = .FALSE.
ukca_config%soa_yield_scaling_isop = rmdi
ukca_config%dry_depvel_so2_scaling = rmdi

! -- Chemistry - Heterogeneous chemistry --
ukca_config%l_ukca_het_psc = .FALSE.
ukca_config%i_ukca_hetconfig = imdi
ukca_config%l_ukca_limit_nat = .FALSE.
ukca_config%l_ukca_sa_clim = .FALSE.
ukca_config%l_ukca_trophet = .FALSE.
ukca_config%l_ukca_classic_hetchem = .FALSE.

! --- Photolysis requirement
ukca_config%l_use_photolysis = .FALSE.

! -- UKCA emissions configuration options --
ukca_config%l_ukca_ibvoc = .FALSE.
ukca_config%l_ukca_inferno = .FALSE.
ukca_config%l_ukca_inferno_ch4 = .FALSE.
ukca_config%i_inferno_emi = imdi
ukca_config%l_ukca_so2ems_expvolc = .FALSE.
ukca_config%l_ukca_so2ems_plumeria = .FALSE.
ukca_config%l_ukca_qch4inter = .FALSE.
ukca_config%mode_parfrac = rmdi
ukca_config%l_ukca_enable_seadms_ems = .FALSE.
ukca_config%l_ukca_emsdrvn_ch4 = .FALSE.
ukca_config%i_ukca_dms_flux = imdi
ukca_config%l_ukca_scale_seadms_ems = .FALSE.
ukca_config%seadms_ems_scaling = rmdi
ukca_config%l_ukca_linox_scaling = .FALSE.
ukca_config%lightnox_scale_fac = rmdi
ukca_config%i_ukca_light_param = imdi
ukca_config%l_support_ems_vertprof = .FALSE.
ukca_config%l_support_ems_gridbox_units = .FALSE.
ukca_config%l_suppress_ems = .FALSE.
ukca_config%anth_so2_ems_scaling = rmdi

! -- UKCA feedback configuration options --
ukca_config%l_ukca_h2o_feedback = .FALSE.
ukca_config%l_ukca_conserve_h = .FALSE.

! -- UKCA environmental driver configuration options --
ukca_config%l_param_conv = .FALSE.
ukca_config%l_ctile = .FALSE.
ukca_config%l_zon_av_ozone = .FALSE.
ukca_config%i_strat_lbc_source = imdi
ukca_config%l_chem_environ_gas_scalars = .FALSE.
ukca_config%l_chem_environ_co2_fld = .FALSE.
ukca_config%l_ukca_prescribech4 = .FALSE.
ukca_config%l_use_classic_so4 = .FALSE.
ukca_config%l_use_classic_soot = .FALSE.
ukca_config%l_use_classic_ocff = .FALSE.
ukca_config%l_use_classic_biogenic = .FALSE.
ukca_config%l_use_classic_seasalt = .FALSE.
ukca_config%l_use_gridbox_volume = .FALSE.
ukca_config%l_use_gridbox_mass = .FALSE.
ukca_config%l_environ_z_top = .FALSE.
ukca_config%env_log_step = imdi

! -- UKCA temporary logicals
ukca_config%l_fix_ukca_cloud_frac = .FALSE.
ukca_config%l_fix_improve_drydep = .FALSE.
ukca_config%l_fix_ukca_h2dd_x = .FALSE.
ukca_config%l_fix_drydep_so2_water = .FALSE.
ukca_config%l_fix_ukca_offox_h2o_fac = .FALSE.
ukca_config%l_fix_ukca_h2so4_ystore = .FALSE.

! -- UKCA internal configuration variables
ukca_config%l_ukca_chem = .FALSE.
ukca_config%l_ukca_trop = .FALSE.
ukca_config%l_ukca_aerchem = .FALSE.
ukca_config%l_ukca_raq = .FALSE.
ukca_config%l_ukca_raqaero = .FALSE.
ukca_config%l_ukca_offline_be = .FALSE.
ukca_config%l_ukca_tropisop = .FALSE.
ukca_config%l_ukca_strattrop = .FALSE.
ukca_config%l_ukca_strat = .FALSE.
ukca_config%l_ukca_offline = .FALSE.
ukca_config%l_ukca_cristrat = .FALSE.
ukca_config%l_ukca_stratcfc = .FALSE.
ukca_config%l_ukca_achem = .FALSE.
ukca_config%l_ukca_nr_aqchem = .FALSE.
ukca_config%l_ukca_advh2o = .FALSE.
ukca_config%l_asad_chem_diags_support = .FALSE.
ukca_config%l_diurnal_isopems = .FALSE.
ukca_config%l_seawater_dms = .FALSE.
ukca_config%l_chem_environ_ch4_scalar = .FALSE.
ukca_config%l_chem_environ_co2_scalar = .FALSE.
ukca_config%l_chem_environ_h2_scalar = .FALSE.
ukca_config%l_chem_environ_n2_scalar = .FALSE.
ukca_config%l_chem_environ_o2_scalar = .FALSE.
ukca_config%ukca_int_method = imdi
ukca_config%timesteps_per_day = imdi
ukca_config%timesteps_per_hour = imdi

! -- General GLOMAP configuration options --
glomap_config%i_mode_nzts = imdi
glomap_config%ukca_mode_seg_size = imdi
glomap_config%i_mode_setup = imdi
glomap_config%l_mode_bhn_on = .FALSE.
glomap_config%l_mode_bln_on = .FALSE.
glomap_config%i_mode_bln_param_method = imdi
glomap_config%mode_activation_dryr = rmdi
glomap_config%l_dust_ageing_on = .FALSE.
glomap_config%dry_depvel_acc_scaling = rmdi
glomap_config%acc_cor_scav_scaling = rmdi

! -- GLOMAP deposition configuration options --
glomap_config%l_ddepaer = .FALSE.
glomap_config%mode_incld_so2_rfrac = rmdi
glomap_config%l_rainout = .FALSE.
glomap_config%l_cv_rainout = .FALSE.
glomap_config%i_mode_nucscav = imdi
glomap_config%l_impc_scav = .FALSE.
glomap_config%l_dust_slinn_impc_scav = .FALSE.

! -- GLOMAP emissions configuration options --
glomap_config%l_ukca_primss = .FALSE.
glomap_config%l_ukca_primsu = .FALSE.
glomap_config%l_ukca_primdu = .FALSE.
glomap_config%l_ukca_primbcoc = .FALSE.
glomap_config%l_ukca_prim_moc = .FALSE.
glomap_config%l_bcoc_bf = .FALSE.
glomap_config%l_bcoc_bm = .FALSE.
glomap_config%l_bcoc_ff = .FALSE.
glomap_config%l_ukca_scale_biom_aer_ems = .FALSE.
glomap_config%biom_aer_ems_scaling = rmdi
glomap_config%l_ukca_fine_no3_prod = .FALSE.
glomap_config%l_ukca_coarse_no3_prod = .FALSE.
glomap_config%l_no3_prod_in_aero_step = .FALSE.
glomap_config%hno3_uptake_coeff = rmdi
glomap_config%l_ukca_scale_sea_salt_ems = .FALSE.
glomap_config%sea_salt_ems_scaling = rmdi
glomap_config%l_ukca_scale_marine_pom_ems = .FALSE.
glomap_config%marine_pom_ems_scaling = rmdi
glomap_config%i_primss_method = imdi

! -- GLOMAP feedback configuration options --
glomap_config%l_ukca_radaer = .FALSE.
glomap_config%i_ukca_tune_bc = imdi
glomap_config%i_ukca_activation_scheme = imdi
glomap_config%i_ukca_nwbins = imdi
glomap_config%sigwmin = rmdi
glomap_config%l_ntpreq_n_activ_sum = .FALSE.
glomap_config%l_ntpreq_dryd_nuc_sol = .FALSE.
glomap_config%l_ukca_sfix = .FALSE.
glomap_config%sigma_updraught_scaling = rmdi

! -- GLOMAP temporary logicals --
glomap_config%l_fix_neg_pvol_wat = .FALSE.
glomap_config%l_fix_ukca_impscav = .FALSE.
glomap_config%l_fix_nacl_density = .FALSE.
glomap_config%l_improve_aero_drydep = .FALSE.
glomap_config%l_fix_ukca_activate_pdf = .FALSE.
glomap_config%l_fix_ukca_activate_vert_rep = .FALSE.
glomap_config%l_bug_repro_tke_index = .FALSE.
glomap_config%l_fix_ukca_hygroscopicities = .FALSE.
glomap_config%l_fix_ukca_water_content = .FALSE.

! -- GLOMAP internal configuration variables --
glomap_config%n_dust_emissions = imdi
glomap_config%i_dust_scheme = imdi
glomap_config%l_6bin_dust_no3 = .FALSE.
glomap_config%l_2bin_dust_no3 = .FALSE.

! -- Configurable constants --
isec_per_day = 86400
isec_per_hour = 3600
rhour_per_day = 24.0
rsec_per_day = 86400.0
rsec_per_hour = 3600.0
pi_over_180 = pi / 180.0
recip_pi_over_180 = 180.0 / pi
planet_radius = 6371229.0
g = 9.80665
r = 287.05
cp = 1005.0
pref = 100000.0
vkman = 0.4
repsilon = 0.62198
rv = r / repsilon
kappa = r / cp
c_virtual = 1.0 / repsilon - 1.0
s2r = 2.0 * pi / rsec_per_day
rad_ait = 6.5e-9
rad_acc = 95.0e-9
chi = 32.0 / 132.0
sigma = 1.4

! -- Parent callback procedures --
NULLIFY(bl_tracer_mix)
NULLIFY(diag2d_copy_out)
NULLIFY(diag3d_copy_out)

RETURN
END SUBROUTINE init_ukca_configuration

! ----------------------------------------------------------------------
SUBROUTINE ukca_get_config(                                                    &
   row_length, rows, model_levels, bl_levels, nlev_ent_tr_mix,                 &
   ntype, npft,                                                                &
   i_brd_leaf, i_brd_leaf_dec, i_brd_leaf_eg_trop, i_brd_leaf_eg_temp,         &
   i_ndl_leaf, i_ndl_leaf_dec, i_ndl_leaf_eg,                                  &
   i_c3_grass, i_c3_crop, i_c3_pasture,                                        &
   i_c4_grass, i_c4_crop, i_c4_pasture,                                        &
   i_shrub, i_shrub_dec, i_shrub_eg,                                           &
   i_urban, i_lake, i_soil, i_ice, i_elev_ice,                                 &
   i_ukca_chem,                                                                &
   fixed_tropopause_level,                                                     &
   i_ageair_reset_method, max_ageair_reset_level,                              &
   i_error_method,                                                             &
   i_ukca_chem_version, nrsteps, chem_timestep,                                &
   dts0, nit,                                                                  &
   i_ukca_quasinewton_start, i_ukca_quasinewton_end,                           &
   ukca_chem_seg_size,                                                         &
   nlev_above_trop_o3_env,                                                     &
   nlev_ch4_stratloss,                                                         &
   i_ukca_topboundary,                                                         &
   i_ukca_hetconfig,                                                           &
   i_inferno_emi,                                                              &
   i_ukca_dms_flux,                                                            &
   i_ukca_light_param,                                                         &
   i_strat_lbc_source,                                                         &
   i_primss_method,                                                            &
   env_log_step,                                                               &
   ukca_int_method,                                                            &
   timesteps_per_day, timesteps_per_hour,                                      &
   i_mode_nzts, ukca_mode_seg_size,                                            &
   i_mode_setup, i_mode_bln_param_method, i_mode_nucscav,                      &
   i_ukca_tune_bc,                                                             &
   i_ukca_activation_scheme, i_ukca_nwbins,                                    &
   n_dust_emissions,                                                           &
   i_dust_scheme,                                                              &
   dzsoil_layer1,                                                              &
   timestep,                                                                   &
   max_ageair_reset_height,                                                    &
   max_z_for_offline_chem,                                                     &
   soa_yield_scaling_mt, soa_yield_scaling_isop,                               &
   dry_depvel_so2_scaling,                                                     &
   fastjx_prescutoff,                                                          &
   mode_parfrac,                                                               &
   seadms_ems_scaling,                                                         &
   sea_salt_ems_scaling, marine_pom_ems_scaling,                               &
   lightnox_scale_fac,                                                         &
   anth_so2_ems_scaling,                                                       &
   mode_activation_dryr,                                                       &
   l_dust_ageing_on,                                                           &
   dry_depvel_acc_scaling,                                                     &
   acc_cor_scav_scaling,                                                       &
   mode_incld_so2_rfrac,                                                       &
   biom_aer_ems_scaling,                                                       &
   hno3_uptake_coeff,                                                          &
   sigwmin,                                                                    &
   sigma_updraught_scaling,                                                    &
   l_cal360,                                                                   &
   l_ukca_chem_aero,                                                           &
   l_ukca_mode,                                                                &
   l_fix_tropopause_level,                                                     &
   l_ukca_ageair,                                                              &
   l_blankout_invalid_diags,                                                   &
   l_enable_diag_um,                                                           &
   l_ukca_persist_off,                                                         &
   l_timer,                                                                    &
   l_ukca_emissions_off,                                                       &
   l_ukca_drydep_off,                                                          &
   l_ukca_wetdep_off,                                                          &
   l_ukca_scale_ppe,                                                           &
   l_ukca_asad_columns,                                                        &
   l_ukca_asad_full,                                                           &
   l_ukca_debug_asad,                                                          &
   l_ukca_intdd, l_ukca_ddepo3_ocean, l_ukca_ddep_lev1, l_ukca_dry_dep_so2wet, &
   l_deposition_jules,                                                         &
   l_ukca_quasinewton,                                                         &
   l_tracer_lumping,                                                           &
   l_ukca_ro2_ntp, l_ukca_ro2_perm,                                            &
   l_ukca_intph,ph_fit_coeff_a,ph_fit_coeff_b,ph_fit_intercept,                &
   l_ukca_scale_soa_yield_mt,l_ukca_scale_soa_yield_isop,                      &
   l_ukca_het_psc,                                                             &
   l_ukca_limit_nat,                                                           &
   l_ukca_sa_clim,                                                             &
   l_ukca_trophet,                                                             &
   l_ukca_classic_hetchem,                                                     &
   l_use_photolysis,                                                           &
   l_ukca_ibvoc,                                                               &
   l_ukca_inferno, l_ukca_inferno_ch4,                                         &
   l_ukca_so2ems_expvolc,                                                      &
   l_ukca_so2ems_plumeria,                                                     &
   l_ukca_qch4inter, l_ukca_emsdrvn_ch4,                                       &
   l_ukca_enable_seadms_ems, l_ukca_scale_seadms_ems,                          &
   l_ukca_linox_scaling,                                                       &
   l_support_ems_vertprof,                                                     &
   l_support_ems_gridbox_units,                                                &
   l_suppress_ems,                                                             &
   l_ukca_h2o_feedback,                                                        &
   l_ukca_conserve_h,                                                          &
   l_param_conv,                                                               &
   l_ctile,                                                                    &
   l_zon_av_ozone,                                                             &
   l_chem_environ_gas_scalars, l_chem_environ_co2_fld, l_ukca_prescribech4,    &
   l_use_classic_so4, l_use_classic_soot, l_use_classic_ocff,                  &
   l_use_classic_biogenic, l_use_classic_seasalt,                              &
   l_use_gridbox_volume,                                                       &
   l_use_gridbox_mass,                                                         &
   l_environ_z_top,                                                            &
   l_fix_ukca_cloud_frac,                                                      &
   l_fix_improve_drydep,                                                       &
   l_fix_ukca_h2dd_x,                                                          &
   l_fix_drydep_so2_water,                                                     &
   l_fix_ukca_offox_h2o_fac,                                                   &
   l_fix_ukca_h2so4_ystore,                                                    &
   l_ukca_chem, l_ukca_trop, l_ukca_aerchem, l_ukca_raq, l_ukca_raqaero,       &
   l_ukca_offline_be, l_ukca_tropisop, l_ukca_strattrop, l_ukca_strat,         &
   l_ukca_offline, l_ukca_cristrat, l_ukca_stratcfc, l_ukca_achem,             &
   l_ukca_nr_aqchem,                                                           &
   l_ukca_advh2o,                                                              &
   l_asad_chem_diags_support,                                                  &
   l_diurnal_isopems,                                                          &
   l_seawater_dms,                                                             &
   l_chem_environ_ch4_scalar,                                                  &
   l_chem_environ_co2_scalar,                                                  &
   l_chem_environ_h2_scalar,                                                   &
   l_chem_environ_n2_scalar,                                                   &
   l_chem_environ_o2_scalar,                                                   &
   l_mode_bhn_on, l_mode_bln_on,                                               &
   l_ddepaer,                                                                  &
   l_rainout, l_cv_rainout,                                                    &
   l_impc_scav, l_dust_slinn_impc_scav,                                        &
   l_ukca_primss, l_ukca_primsu, l_ukca_primdu, l_ukca_primbcoc,               &
   l_ukca_prim_moc, l_bcoc_bf, l_bcoc_bm, l_bcoc_ff,                           &
   l_ukca_scale_biom_aer_ems,                                                  &
   l_ukca_fine_no3_prod,                                                       &
   l_ukca_coarse_no3_prod,                                                     &
   l_no3_prod_in_aero_step,                                                    &
   l_ukca_scale_sea_salt_ems,                                                  &
   l_ukca_scale_marine_pom_ems,                                                &
   l_ukca_radaer,                                                              &
   l_ntpreq_n_activ_sum,                                                       &
   l_ntpreq_dryd_nuc_sol,                                                      &
   l_ukca_sfix,                                                                &
   l_fix_neg_pvol_wat,                                                         &
   l_fix_ukca_impscav,                                                         &
   l_fix_nacl_density,                                                         &
   l_improve_aero_drydep,                                                      &
   l_fix_ukca_activate_pdf,                                                    &
   l_fix_ukca_activate_vert_rep,                                               &
   l_bug_repro_tke_index,                                                      &
   l_6bin_dust_no3,                                                            &
   l_2bin_dust_no3,                                                            &
   l_fix_ukca_hygroscopicities,                                                &
   l_fix_ukca_water_content,                                                   &
   l_config_available)
! ----------------------------------------------------------------------
! Description:
!   Returns values of all configuration variables that are supplied as
!   optional arguments.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
! (follow 'ukca_config_spec_type' & 'glomap_config_spec_type' order
!  within type groups)

INTEGER, OPTIONAL, INTENT(OUT) :: row_length
INTEGER, OPTIONAL, INTENT(OUT) :: rows
INTEGER, OPTIONAL, INTENT(OUT) :: model_levels
INTEGER, OPTIONAL, INTENT(OUT) :: bl_levels
INTEGER, OPTIONAL, INTENT(OUT) :: nlev_ent_tr_mix
INTEGER, OPTIONAL, INTENT(OUT) :: ntype
INTEGER, OPTIONAL, INTENT(OUT) :: npft
INTEGER, OPTIONAL, INTENT(OUT) :: i_brd_leaf
INTEGER, OPTIONAL, INTENT(OUT) :: i_brd_leaf_dec
INTEGER, OPTIONAL, INTENT(OUT) :: i_brd_leaf_eg_trop
INTEGER, OPTIONAL, INTENT(OUT) :: i_brd_leaf_eg_temp
INTEGER, OPTIONAL, INTENT(OUT) :: i_ndl_leaf
INTEGER, OPTIONAL, INTENT(OUT) :: i_ndl_leaf_dec
INTEGER, OPTIONAL, INTENT(OUT) :: i_ndl_leaf_eg
INTEGER, OPTIONAL, INTENT(OUT) :: i_c3_grass
INTEGER, OPTIONAL, INTENT(OUT) :: i_c3_crop
INTEGER, OPTIONAL, INTENT(OUT) :: i_c3_pasture
INTEGER, OPTIONAL, INTENT(OUT) :: i_c4_grass
INTEGER, OPTIONAL, INTENT(OUT) :: i_c4_crop
INTEGER, OPTIONAL, INTENT(OUT) :: i_c4_pasture
INTEGER, OPTIONAL, INTENT(OUT) :: i_shrub
INTEGER, OPTIONAL, INTENT(OUT) :: i_shrub_dec
INTEGER, OPTIONAL, INTENT(OUT) :: i_shrub_eg
INTEGER, OPTIONAL, INTENT(OUT) :: i_urban
INTEGER, OPTIONAL, INTENT(OUT) :: i_lake
INTEGER, OPTIONAL, INTENT(OUT) :: i_soil
INTEGER, OPTIONAL, INTENT(OUT) :: i_ice
INTEGER, ALLOCATABLE, OPTIONAL, INTENT(OUT) :: i_elev_ice(:)
INTEGER, OPTIONAL, INTENT(OUT) :: i_ukca_chem
INTEGER, OPTIONAL, INTENT(OUT) :: fixed_tropopause_level
INTEGER, OPTIONAL, INTENT(OUT) :: i_ageair_reset_method
INTEGER, OPTIONAL, INTENT(OUT) :: max_ageair_reset_level
INTEGER, OPTIONAL, INTENT(OUT) :: i_error_method
INTEGER, OPTIONAL, INTENT(OUT) :: i_ukca_chem_version
INTEGER, OPTIONAL, INTENT(OUT) :: nrsteps
INTEGER, OPTIONAL, INTENT(OUT) :: chem_timestep
INTEGER, OPTIONAL, INTENT(OUT) :: dts0
INTEGER, OPTIONAL, INTENT(OUT) :: nit
INTEGER, OPTIONAL, INTENT(OUT) :: i_ukca_quasinewton_start
INTEGER, OPTIONAL, INTENT(OUT) :: i_ukca_quasinewton_end
INTEGER, OPTIONAL, INTENT(OUT) :: ukca_chem_seg_size
INTEGER, OPTIONAL, INTENT(OUT) :: nlev_above_trop_o3_env
INTEGER, OPTIONAL, INTENT(OUT) :: nlev_ch4_stratloss
INTEGER, OPTIONAL, INTENT(OUT) :: i_ukca_topboundary
INTEGER, OPTIONAL, INTENT(OUT) :: i_ukca_hetconfig
INTEGER, OPTIONAL, INTENT(OUT) :: i_inferno_emi
INTEGER, OPTIONAL, INTENT(OUT) :: i_ukca_dms_flux
INTEGER, OPTIONAL, INTENT(OUT) :: i_ukca_light_param
INTEGER, OPTIONAL, INTENT(OUT) :: i_strat_lbc_source
INTEGER, OPTIONAL, INTENT(OUT) :: env_log_step
INTEGER, OPTIONAL, INTENT(OUT) :: ukca_int_method
INTEGER, OPTIONAL, INTENT(OUT) :: timesteps_per_day
INTEGER, OPTIONAL, INTENT(OUT) :: timesteps_per_hour
INTEGER, OPTIONAL, INTENT(OUT) :: i_mode_nzts
INTEGER, OPTIONAL, INTENT(OUT) :: ukca_mode_seg_size
INTEGER, OPTIONAL, INTENT(OUT) :: i_mode_setup
INTEGER, OPTIONAL, INTENT(OUT) :: i_mode_bln_param_method
INTEGER, OPTIONAL, INTENT(OUT) :: i_mode_nucscav
INTEGER, OPTIONAL, INTENT(OUT) :: i_ukca_tune_bc
INTEGER, OPTIONAL, INTENT(OUT) :: i_ukca_activation_scheme
INTEGER, OPTIONAL, INTENT(OUT) :: i_ukca_nwbins
INTEGER, OPTIONAL, INTENT(OUT) :: n_dust_emissions
INTEGER, OPTIONAL, INTENT(OUT) :: i_dust_scheme
INTEGER, OPTIONAL, INTENT(OUT) :: i_primss_method

REAL, OPTIONAL, INTENT(OUT) :: dzsoil_layer1
REAL, OPTIONAL, INTENT(OUT) :: timestep
REAL, OPTIONAL, INTENT(OUT) :: max_ageair_reset_height
REAL, OPTIONAL, INTENT(OUT) :: max_z_for_offline_chem
REAL, OPTIONAL, INTENT(OUT) :: soa_yield_scaling_mt
REAL, OPTIONAL, INTENT(OUT) :: soa_yield_scaling_isop
REAL, OPTIONAL, INTENT(OUT) :: dry_depvel_so2_scaling
REAL, OPTIONAL, INTENT(OUT) :: fastjx_prescutoff
REAL, OPTIONAL, INTENT(OUT) :: mode_parfrac
REAL, OPTIONAL, INTENT(OUT) :: seadms_ems_scaling
REAL, OPTIONAL, INTENT(OUT) :: sea_salt_ems_scaling
REAL, OPTIONAL, INTENT(OUT) :: marine_pom_ems_scaling
REAL, OPTIONAL, INTENT(OUT) :: lightnox_scale_fac
REAL, OPTIONAL, INTENT(OUT) :: anth_so2_ems_scaling
REAL, OPTIONAL, INTENT(OUT) :: mode_activation_dryr
REAL, OPTIONAL, INTENT(OUT) :: dry_depvel_acc_scaling
REAL, OPTIONAL, INTENT(OUT) :: acc_cor_scav_scaling
REAL, OPTIONAL, INTENT(OUT) :: mode_incld_so2_rfrac
REAL, OPTIONAL, INTENT(OUT) :: biom_aer_ems_scaling
REAL, OPTIONAL, INTENT(OUT) :: sigwmin
REAL, OPTIONAL, INTENT(OUT) :: sigma_updraught_scaling
REAL, OPTIONAL, INTENT(OUT) :: ph_fit_coeff_a
REAL, OPTIONAL, INTENT(OUT) :: ph_fit_coeff_b
REAL, OPTIONAL, INTENT(OUT) :: ph_fit_intercept
REAL, OPTIONAL, INTENT(OUT) :: hno3_uptake_coeff

LOGICAL, OPTIONAL, INTENT(OUT) :: l_cal360
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_chem_aero
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_mode
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_tropopause_level
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_ageair
LOGICAL, OPTIONAL, INTENT(OUT) :: l_blankout_invalid_diags
LOGICAL, OPTIONAL, INTENT(OUT) :: l_enable_diag_um
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_persist_off
LOGICAL, OPTIONAL, INTENT(OUT) :: l_timer
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_emissions_off
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_drydep_off
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_wetdep_off
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_scale_ppe
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_asad_columns
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_asad_full
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_debug_asad
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_intdd
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_ddepo3_ocean
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_ddep_lev1
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_dry_dep_so2wet
LOGICAL, OPTIONAL, INTENT(OUT) :: l_deposition_jules
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_quasinewton
LOGICAL, OPTIONAL, INTENT(OUT) :: l_tracer_lumping
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_ro2_ntp
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_ro2_perm
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_intph
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_scale_soa_yield_mt
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_scale_soa_yield_isop
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_het_psc
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_limit_nat
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_sa_clim
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_trophet
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_classic_hetchem
LOGICAL, OPTIONAL, INTENT(OUT) :: l_use_photolysis
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_ibvoc
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_inferno
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_inferno_ch4
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_so2ems_expvolc
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_so2ems_plumeria
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_qch4inter
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_emsdrvn_ch4
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_enable_seadms_ems
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_scale_seadms_ems
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_linox_scaling
LOGICAL, OPTIONAL, INTENT(OUT) :: l_support_ems_vertprof
LOGICAL, OPTIONAL, INTENT(OUT) :: l_support_ems_gridbox_units
LOGICAL, OPTIONAL, INTENT(OUT) :: l_suppress_ems
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_h2o_feedback
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_conserve_h
LOGICAL, OPTIONAL, INTENT(OUT) :: l_param_conv
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ctile
LOGICAL, OPTIONAL, INTENT(OUT) :: l_zon_av_ozone
LOGICAL, OPTIONAL, INTENT(OUT) :: l_chem_environ_gas_scalars
LOGICAL, OPTIONAL, INTENT(OUT) :: l_chem_environ_co2_fld
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_prescribech4
LOGICAL, OPTIONAL, INTENT(OUT) :: l_use_classic_so4
LOGICAL, OPTIONAL, INTENT(OUT) :: l_use_classic_soot
LOGICAL, OPTIONAL, INTENT(OUT) :: l_use_classic_ocff
LOGICAL, OPTIONAL, INTENT(OUT) :: l_use_classic_biogenic
LOGICAL, OPTIONAL, INTENT(OUT) :: l_use_classic_seasalt
LOGICAL, OPTIONAL, INTENT(OUT) :: l_use_gridbox_volume
LOGICAL, OPTIONAL, INTENT(OUT) :: l_use_gridbox_mass
LOGICAL, OPTIONAL, INTENT(OUT) :: l_environ_z_top
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_ukca_cloud_frac
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_improve_drydep
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_ukca_h2dd_x
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_drydep_so2_water
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_ukca_offox_h2o_fac
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_ukca_h2so4_ystore
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_chem
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_trop
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_aerchem
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_raq
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_raqaero
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_offline_be
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_tropisop
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_strattrop
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_strat
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_offline
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_cristrat
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_stratcfc
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_achem
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_nr_aqchem
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_advh2o
LOGICAL, OPTIONAL, INTENT(OUT) :: l_asad_chem_diags_support
LOGICAL, OPTIONAL, INTENT(OUT) :: l_diurnal_isopems
LOGICAL, OPTIONAL, INTENT(OUT) :: l_seawater_dms
LOGICAL, OPTIONAL, INTENT(OUT) :: l_chem_environ_ch4_scalar
LOGICAL, OPTIONAL, INTENT(OUT) :: l_chem_environ_co2_scalar
LOGICAL, OPTIONAL, INTENT(OUT) :: l_chem_environ_h2_scalar
LOGICAL, OPTIONAL, INTENT(OUT) :: l_chem_environ_n2_scalar
LOGICAL, OPTIONAL, INTENT(OUT) :: l_chem_environ_o2_scalar
LOGICAL, OPTIONAL, INTENT(OUT) :: l_mode_bhn_on
LOGICAL, OPTIONAL, INTENT(OUT) :: l_mode_bln_on
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ddepaer
LOGICAL, OPTIONAL, INTENT(OUT) :: l_rainout
LOGICAL, OPTIONAL, INTENT(OUT) :: l_cv_rainout
LOGICAL, OPTIONAL, INTENT(OUT) :: l_impc_scav
LOGICAL, OPTIONAL, INTENT(OUT) :: l_dust_slinn_impc_scav
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_primss
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_primsu
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_primdu
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_primbcoc
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_prim_moc
LOGICAL, OPTIONAL, INTENT(OUT) :: l_bcoc_bf
LOGICAL, OPTIONAL, INTENT(OUT) :: l_bcoc_bm
LOGICAL, OPTIONAL, INTENT(OUT) :: l_bcoc_ff
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_scale_biom_aer_ems
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_fine_no3_prod
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_coarse_no3_prod
LOGICAL, OPTIONAL, INTENT(OUT) :: l_no3_prod_in_aero_step
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_scale_sea_salt_ems
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_scale_marine_pom_ems
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_radaer
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ntpreq_n_activ_sum
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ntpreq_dryd_nuc_sol
LOGICAL, OPTIONAL, INTENT(OUT) :: l_ukca_sfix
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_neg_pvol_wat
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_ukca_impscav
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_nacl_density
LOGICAL, OPTIONAL, INTENT(OUT) :: l_improve_aero_drydep
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_ukca_activate_pdf
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_ukca_activate_vert_rep
LOGICAL, OPTIONAL, INTENT(OUT) :: l_bug_repro_tke_index
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_ukca_hygroscopicities
LOGICAL, OPTIONAL, INTENT(OUT) :: l_fix_ukca_water_content
LOGICAL, OPTIONAL, INTENT(OUT) :: l_6bin_dust_no3
LOGICAL, OPTIONAL, INTENT(OUT) :: l_2bin_dust_no3
LOGICAL, OPTIONAL, INTENT(OUT) :: l_config_available
LOGICAL, OPTIONAL, INTENT(OUT) :: l_dust_ageing_on

! -- Availability of a valid configuration --
IF (PRESENT(l_config_available)) l_config_available = l_ukca_config_available

! -- Context information --
IF (PRESENT(row_length)) row_length = ukca_config%row_length
IF (PRESENT(rows)) rows = ukca_config%rows
IF (PRESENT(model_levels)) model_levels = ukca_config%model_levels
IF (PRESENT(bl_levels)) bl_levels = ukca_config%bl_levels
IF (PRESENT(nlev_ent_tr_mix)) nlev_ent_tr_mix = ukca_config%nlev_ent_tr_mix
IF (PRESENT(ntype)) ntype = ukca_config%ntype
IF (PRESENT(npft)) npft = ukca_config%npft
IF (PRESENT(i_brd_leaf)) i_brd_leaf = ukca_config%i_brd_leaf
IF (PRESENT(i_brd_leaf_dec)) i_brd_leaf_dec = ukca_config%i_brd_leaf_dec
IF (PRESENT(i_brd_leaf_eg_trop))                                               &
  i_brd_leaf_eg_trop = ukca_config%i_brd_leaf_eg_trop
IF (PRESENT(i_brd_leaf_eg_temp))                                               &
  i_brd_leaf_eg_temp = ukca_config%i_brd_leaf_eg_temp
IF (PRESENT(i_ndl_leaf)) i_ndl_leaf = ukca_config%i_ndl_leaf
IF (PRESENT(i_ndl_leaf_dec)) i_ndl_leaf_dec = ukca_config%i_ndl_leaf_dec
IF (PRESENT(i_ndl_leaf_eg)) i_ndl_leaf_eg = ukca_config%i_ndl_leaf_eg
IF (PRESENT(i_c3_grass)) i_c3_grass = ukca_config%i_c3_grass
IF (PRESENT(i_c3_crop)) i_c3_crop = ukca_config%i_c3_crop
IF (PRESENT(i_c3_pasture)) i_c3_pasture = ukca_config%i_c3_pasture
IF (PRESENT(i_c4_grass)) i_c4_grass = ukca_config%i_c4_grass
IF (PRESENT(i_c4_crop)) i_c4_crop = ukca_config%i_c4_crop
IF (PRESENT(i_c4_pasture)) i_c4_pasture = ukca_config%i_c4_pasture
IF (PRESENT(i_shrub)) i_shrub = ukca_config%i_shrub
IF (PRESENT(i_shrub_dec)) i_shrub_dec = ukca_config%i_shrub_dec
IF (PRESENT(i_shrub_eg)) i_shrub_eg = ukca_config%i_shrub_eg
IF (PRESENT(i_urban)) i_urban = ukca_config%i_urban
IF (PRESENT(i_lake)) i_lake = ukca_config%i_lake
IF (PRESENT(i_soil)) i_soil = ukca_config%i_soil
IF (PRESENT(i_ice)) i_ice = ukca_config%i_ice
IF (PRESENT(i_elev_ice))                                                       &
  CALL copy_config_vector(ukca_config%i_elev_ice, i_elev_ice)
IF (PRESENT(dzsoil_layer1)) dzsoil_layer1 = ukca_config%dzsoil_layer1
IF (PRESENT(l_cal360)) l_cal360 = ukca_config%l_cal360
IF (PRESENT(timestep)) timestep = ukca_config%timestep

! -- General UKCA configuration options --
IF (PRESENT(i_ukca_chem)) i_ukca_chem = ukca_config%i_ukca_chem
IF (PRESENT(l_ukca_chem_aero)) l_ukca_chem_aero = ukca_config%l_ukca_chem_aero
IF (PRESENT(l_ukca_mode)) l_ukca_mode = ukca_config%l_ukca_mode
IF (PRESENT(l_fix_tropopause_level))                                           &
  l_fix_tropopause_level = ukca_config%l_fix_tropopause_level
IF (PRESENT(fixed_tropopause_level))                                           &
  fixed_tropopause_level = ukca_config%fixed_tropopause_level
IF (PRESENT(l_ukca_ageair)) l_ukca_ageair = ukca_config%l_ukca_ageair
IF (PRESENT(i_ageair_reset_method))                                            &
  i_ageair_reset_method = ukca_config%i_ageair_reset_method
IF (PRESENT(max_ageair_reset_level))                                           &
  max_ageair_reset_level = ukca_config%max_ageair_reset_level
IF (PRESENT(max_ageair_reset_height))                                          &
  max_ageair_reset_height = ukca_config%max_ageair_reset_height
IF (PRESENT(l_blankout_invalid_diags))                                         &
  l_blankout_invalid_diags = ukca_config%l_blankout_invalid_diags
IF (PRESENT(l_enable_diag_um)) l_enable_diag_um = ukca_config%l_enable_diag_um
IF (PRESENT(l_ukca_persist_off))                                               &
  l_ukca_persist_off = ukca_config%l_ukca_persist_off
IF (PRESENT(l_timer)) l_timer = ukca_config%l_timer
IF (PRESENT(l_ukca_emissions_off))                                             &
  l_ukca_emissions_off = ukca_config%l_ukca_emissions_off
IF (PRESENT(l_ukca_drydep_off))                                                &
  l_ukca_drydep_off = ukca_config%l_ukca_drydep_off
IF (PRESENT(l_ukca_wetdep_off))                                                &
  l_ukca_wetdep_off = ukca_config%l_ukca_wetdep_off
IF (PRESENT(i_error_method))                                                   &
  i_error_method = ukca_config%i_error_method
IF (PRESENT(l_ukca_scale_ppe))                                                 &
  l_ukca_scale_ppe = ukca_config%l_ukca_scale_ppe

! -- Chemistry configuration options --
IF (PRESENT(i_ukca_chem_version))                                              &
  i_ukca_chem_version = ukca_config%i_ukca_chem_version
IF (PRESENT(nrsteps)) nrsteps = ukca_config%nrsteps
IF (PRESENT(chem_timestep)) chem_timestep = ukca_config%chem_timestep
IF (PRESENT(dts0)) dts0 = ukca_config%dts0
IF (PRESENT(l_ukca_asad_columns))                                              &
  l_ukca_asad_columns = ukca_config%l_ukca_asad_columns
IF (PRESENT(l_ukca_asad_full))                                                 &
  l_ukca_asad_full = ukca_config%l_ukca_asad_full
IF (PRESENT(l_ukca_debug_asad))                                                &
  l_ukca_debug_asad = ukca_config%l_ukca_debug_asad
IF (PRESENT(l_ukca_intdd)) l_ukca_intdd = ukca_config%l_ukca_intdd
IF (PRESENT(l_ukca_ddepo3_ocean))                                              &
  l_ukca_ddepo3_ocean = ukca_config%l_ukca_ddepo3_ocean
IF (PRESENT(l_ukca_ddep_lev1)) l_ukca_ddep_lev1 = ukca_config%l_ukca_ddep_lev1
IF (PRESENT(l_ukca_dry_dep_so2wet))                                            &
  l_ukca_dry_dep_so2wet = ukca_config%l_ukca_dry_dep_so2wet
IF (PRESENT(l_deposition_jules))                                               &
  l_deposition_jules = ukca_config%l_deposition_jules
IF (PRESENT(nit)) nit = ukca_config%nit
IF (PRESENT(l_ukca_quasinewton))                                               &
  l_ukca_quasinewton = ukca_config%l_ukca_quasinewton
IF (PRESENT(i_ukca_quasinewton_start))                                         &
  i_ukca_quasinewton_start = ukca_config%i_ukca_quasinewton_start
IF (PRESENT(i_ukca_quasinewton_end))                                           &
  i_ukca_quasinewton_end = ukca_config%i_ukca_quasinewton_end
IF (PRESENT(ukca_chem_seg_size))                                               &
   ukca_chem_seg_size = ukca_config%ukca_chem_seg_size
IF (PRESENT(max_z_for_offline_chem))                                           &
  max_z_for_offline_chem = ukca_config%max_z_for_offline_chem
IF (PRESENT(nlev_above_trop_o3_env))                                           &
   nlev_above_trop_o3_env = ukca_config%nlev_above_trop_o3_env
IF (PRESENT(nlev_ch4_stratloss))                                               &
   nlev_ch4_stratloss = ukca_config%nlev_ch4_stratloss
IF (PRESENT(l_tracer_lumping))                                                 &
  l_tracer_lumping = ukca_config%l_tracer_lumping
IF (PRESENT(i_ukca_topboundary))                                               &
  i_ukca_topboundary = ukca_config%i_ukca_topboundary
IF (PRESENT(l_ukca_ro2_ntp)) l_ukca_ro2_ntp = ukca_config%l_ukca_ro2_ntp
IF (PRESENT(l_ukca_ro2_perm)) l_ukca_ro2_perm = ukca_config%l_ukca_ro2_perm
IF (PRESENT(l_ukca_intph)) l_ukca_intph = ukca_config%l_ukca_intph
IF (PRESENT(ph_fit_coeff_a)) ph_fit_coeff_a = ukca_config%ph_fit_coeff_a
IF (PRESENT(ph_fit_coeff_b)) ph_fit_coeff_b = ukca_config%ph_fit_coeff_b
IF (PRESENT(ph_fit_intercept)) ph_fit_intercept = ukca_config%ph_fit_intercept
IF (PRESENT(l_ukca_scale_soa_yield_mt))                                        &
  l_ukca_scale_soa_yield_mt = ukca_config%l_ukca_scale_soa_yield_mt
IF (PRESENT(soa_yield_scaling_mt))                                             &
  soa_yield_scaling_mt = ukca_config%soa_yield_scaling_mt
IF (PRESENT(l_ukca_scale_soa_yield_isop))                                      &
  l_ukca_scale_soa_yield_isop = ukca_config%l_ukca_scale_soa_yield_isop
IF (PRESENT(soa_yield_scaling_isop))                                           &
  soa_yield_scaling_isop = ukca_config%soa_yield_scaling_isop
IF (PRESENT(dry_depvel_so2_scaling))                                           &
  dry_depvel_so2_scaling = ukca_config%dry_depvel_so2_scaling

! -- Chemistry - Heterogeneous chemistry --
IF (PRESENT(l_ukca_het_psc)) l_ukca_het_psc = ukca_config%l_ukca_het_psc
IF (PRESENT(i_ukca_hetconfig)) i_ukca_hetconfig = ukca_config%i_ukca_hetconfig
IF (PRESENT(l_ukca_limit_nat)) l_ukca_limit_nat = ukca_config%l_ukca_limit_nat
IF (PRESENT(l_ukca_sa_clim)) l_ukca_sa_clim = ukca_config%l_ukca_sa_clim
IF (PRESENT(l_ukca_trophet)) l_ukca_trophet = ukca_config%l_ukca_trophet
IF (PRESENT(l_ukca_classic_hetchem))                                           &
  l_ukca_classic_hetchem = ukca_config%l_ukca_classic_hetchem

! -- Photolysis requirement --
IF (PRESENT(l_use_photolysis)) l_use_photolysis = ukca_config%l_use_photolysis
! -- UKCA emissions configuration options --
IF (PRESENT(l_ukca_ibvoc)) l_ukca_ibvoc = ukca_config%l_ukca_ibvoc
IF (PRESENT(l_ukca_inferno)) l_ukca_inferno = ukca_config%l_ukca_inferno
IF (PRESENT(l_ukca_inferno_ch4))                                               &
  l_ukca_inferno_ch4 = ukca_config%l_ukca_inferno_ch4
IF (PRESENT(i_inferno_emi)) i_inferno_emi = ukca_config%i_inferno_emi
IF (PRESENT(l_ukca_so2ems_expvolc))                                            &
  l_ukca_so2ems_expvolc = ukca_config%l_ukca_so2ems_expvolc
IF (PRESENT(l_ukca_so2ems_plumeria))                                           &
  l_ukca_so2ems_plumeria = ukca_config%l_ukca_so2ems_plumeria
IF (PRESENT(l_ukca_qch4inter)) l_ukca_qch4inter = ukca_config%l_ukca_qch4inter
IF (PRESENT(l_ukca_emsdrvn_ch4))                                               &
l_ukca_emsdrvn_ch4 = ukca_config%l_ukca_emsdrvn_ch4
IF (PRESENT(mode_parfrac)) mode_parfrac = ukca_config%mode_parfrac
IF (PRESENT(l_ukca_enable_seadms_ems))                                         &
  l_ukca_enable_seadms_ems = ukca_config%l_ukca_enable_seadms_ems
IF (PRESENT(i_ukca_dms_flux)) i_ukca_dms_flux = ukca_config%i_ukca_dms_flux
IF (PRESENT(l_ukca_scale_seadms_ems))                                          &
  l_ukca_scale_seadms_ems = ukca_config%l_ukca_scale_seadms_ems
IF (PRESENT(seadms_ems_scaling))                                               &
  seadms_ems_scaling = ukca_config%seadms_ems_scaling
IF (PRESENT(l_ukca_linox_scaling))                                             &
  l_ukca_linox_scaling = ukca_config%l_ukca_linox_scaling
IF (PRESENT(lightnox_scale_fac))                                               &
  lightnox_scale_fac = ukca_config%lightnox_scale_fac
IF (PRESENT(i_ukca_light_param))                                               &
  i_ukca_light_param = ukca_config%i_ukca_light_param
IF (PRESENT(l_support_ems_vertprof))                                           &
  l_support_ems_vertprof = ukca_config%l_support_ems_vertprof
IF (PRESENT(l_support_ems_gridbox_units))                                      &
  l_support_ems_gridbox_units = ukca_config%l_support_ems_gridbox_units
IF (PRESENT(l_suppress_ems))                                                   &
  l_suppress_ems = ukca_config%l_suppress_ems
IF (PRESENT(anth_so2_ems_scaling))                                             &
  anth_so2_ems_scaling = ukca_config%anth_so2_ems_scaling

! -- UKCA feedback configuration options --
IF (PRESENT(l_ukca_h2o_feedback))                                              &
  l_ukca_h2o_feedback = ukca_config%l_ukca_h2o_feedback
IF (PRESENT(l_ukca_conserve_h))                                                &
  l_ukca_conserve_h = ukca_config%l_ukca_conserve_h

! -- UKCA environmental driver configuration options --
IF (PRESENT(l_param_conv)) l_param_conv = ukca_config%l_param_conv
IF (PRESENT(l_ctile)) l_ctile = ukca_config%l_ctile
IF (PRESENT(l_zon_av_ozone)) l_zon_av_ozone = ukca_config%l_zon_av_ozone
IF (PRESENT(i_strat_lbc_source))                                               &
  i_strat_lbc_source = ukca_config%i_strat_lbc_source
IF (PRESENT(l_chem_environ_gas_scalars))                                       &
  l_chem_environ_gas_scalars = ukca_config%l_chem_environ_gas_scalars
IF (PRESENT(l_chem_environ_co2_fld))                                           &
  l_chem_environ_co2_fld = ukca_config%l_chem_environ_co2_fld
IF (PRESENT(l_ukca_prescribech4))                                              &
  l_ukca_prescribech4 = ukca_config%l_ukca_prescribech4
IF (PRESENT(l_use_classic_so4))                                                &
  l_use_classic_so4 = ukca_config%l_use_classic_so4
IF (PRESENT(l_use_classic_soot))                                               &
  l_use_classic_soot = ukca_config%l_use_classic_soot
IF (PRESENT(l_use_classic_ocff))                                               &
  l_use_classic_ocff = ukca_config%l_use_classic_ocff
IF (PRESENT(l_use_classic_biogenic))                                           &
  l_use_classic_biogenic = ukca_config%l_use_classic_biogenic
IF (PRESENT(l_use_classic_seasalt))                                            &
  l_use_classic_seasalt = ukca_config%l_use_classic_seasalt
IF (PRESENT(l_environ_z_top)) l_environ_z_top = ukca_config%l_environ_z_top
IF (PRESENT(l_use_gridbox_volume))                                             &
  l_use_gridbox_volume = ukca_config%l_use_gridbox_volume
IF (PRESENT(l_use_gridbox_mass))                                               &
  l_use_gridbox_mass = ukca_config%l_use_gridbox_mass
IF (PRESENT(env_log_step)) env_log_step = ukca_config%env_log_step

! -- UKCA temporary logicals
IF (PRESENT(l_fix_ukca_cloud_frac))                                            &
  l_fix_ukca_cloud_frac = ukca_config%l_fix_ukca_cloud_frac
IF (PRESENT(l_fix_improve_drydep))                                             &
  l_fix_improve_drydep = ukca_config%l_fix_improve_drydep
IF (PRESENT(l_fix_ukca_h2dd_x))                                                &
  l_fix_ukca_h2dd_x = ukca_config%l_fix_ukca_h2dd_x
IF (PRESENT(l_fix_drydep_so2_water))                                           &
  l_fix_drydep_so2_water = ukca_config%l_fix_drydep_so2_water
IF (PRESENT(l_fix_ukca_offox_h2o_fac))                                         &
  l_fix_ukca_offox_h2o_fac = ukca_config%l_fix_ukca_offox_h2o_fac
IF (PRESENT(l_fix_ukca_h2so4_ystore))                                          &
  l_fix_ukca_h2so4_ystore = ukca_config%l_fix_ukca_h2so4_ystore

! -- UKCA internal configuration variables
IF (PRESENT(l_ukca_chem)) l_ukca_chem = ukca_config%l_ukca_chem
IF (PRESENT(l_ukca_trop)) l_ukca_trop = ukca_config%l_ukca_trop
IF (PRESENT(l_ukca_aerchem)) l_ukca_aerchem = ukca_config%l_ukca_aerchem
IF (PRESENT(l_ukca_raq)) l_ukca_raq = ukca_config%l_ukca_raq
IF (PRESENT(l_ukca_raqaero)) l_ukca_raqaero = ukca_config%l_ukca_raqaero
IF (PRESENT(l_ukca_offline_be))                                                &
  l_ukca_offline_be = ukca_config%l_ukca_offline_be
IF (PRESENT(l_ukca_tropisop)) l_ukca_tropisop = ukca_config%l_ukca_tropisop
IF (PRESENT(l_ukca_strattrop)) l_ukca_strattrop = ukca_config%l_ukca_strattrop
IF (PRESENT(l_ukca_strat)) l_ukca_strat = ukca_config%l_ukca_strat
IF (PRESENT(l_ukca_offline)) l_ukca_offline = ukca_config%l_ukca_offline
IF (PRESENT(l_ukca_cristrat)) l_ukca_cristrat = ukca_config%l_ukca_cristrat
IF (PRESENT(l_ukca_stratcfc)) l_ukca_stratcfc = ukca_config%l_ukca_stratcfc
IF (PRESENT(l_ukca_achem)) l_ukca_achem = ukca_config%l_ukca_achem
IF (PRESENT(l_ukca_nr_aqchem)) l_ukca_nr_aqchem = ukca_config%l_ukca_nr_aqchem
IF (PRESENT(l_ukca_advh2o)) l_ukca_advh2o = ukca_config%l_ukca_advh2o
IF (PRESENT(l_asad_chem_diags_support))                                        &
  l_asad_chem_diags_support = ukca_config%l_asad_chem_diags_support
IF (PRESENT(l_diurnal_isopems))                                                &
  l_diurnal_isopems = ukca_config%l_diurnal_isopems
IF (PRESENT(l_seawater_dms)) l_seawater_dms = ukca_config%l_seawater_dms
IF (PRESENT(l_chem_environ_ch4_scalar))                                        &
  l_chem_environ_ch4_scalar = ukca_config%l_chem_environ_ch4_scalar
IF (PRESENT(l_chem_environ_co2_scalar))                                        &
  l_chem_environ_co2_scalar = ukca_config%l_chem_environ_co2_scalar
IF (PRESENT(l_chem_environ_h2_scalar))                                         &
  l_chem_environ_h2_scalar = ukca_config%l_chem_environ_h2_scalar
IF (PRESENT(l_chem_environ_n2_scalar))                                         &
  l_chem_environ_n2_scalar = ukca_config%l_chem_environ_n2_scalar
IF (PRESENT(l_chem_environ_o2_scalar))                                         &
  l_chem_environ_o2_scalar = ukca_config%l_chem_environ_o2_scalar
IF (PRESENT(ukca_int_method)) ukca_int_method = ukca_config%ukca_int_method
IF (PRESENT(timesteps_per_day))                                                &
  timesteps_per_day = ukca_config%timesteps_per_day
IF (PRESENT(timesteps_per_hour))                                               &
  timesteps_per_hour = ukca_config%timesteps_per_hour

! -- General GLOMAP configuration options --
IF (PRESENT(i_mode_nzts)) i_mode_nzts = glomap_config%i_mode_nzts
IF (PRESENT(ukca_mode_seg_size))                                               &
  ukca_mode_seg_size = glomap_config%ukca_mode_seg_size
IF (PRESENT(i_mode_setup)) i_mode_setup = glomap_config%i_mode_setup
IF (PRESENT(l_mode_bhn_on)) l_mode_bhn_on = glomap_config%l_mode_bhn_on
IF (PRESENT(l_mode_bln_on)) l_mode_bln_on = glomap_config%l_mode_bln_on
IF (PRESENT(i_mode_bln_param_method))                                          &
  i_mode_bln_param_method = glomap_config%i_mode_bln_param_method
IF (PRESENT(mode_activation_dryr))                                             &
  mode_activation_dryr = glomap_config%mode_activation_dryr
IF (PRESENT(l_dust_ageing_on))                                                 &
  l_dust_ageing_on = glomap_config%l_dust_ageing_on
IF (PRESENT(dry_depvel_acc_scaling))                                           &
  dry_depvel_acc_scaling = glomap_config%dry_depvel_acc_scaling
IF (PRESENT(acc_cor_scav_scaling))                                             &
  acc_cor_scav_scaling = glomap_config%acc_cor_scav_scaling

! -- GLOMAP deposition configuration optons --
IF (PRESENT(l_ddepaer))                                                        &
  l_ddepaer = glomap_config%l_ddepaer
IF (PRESENT(mode_incld_so2_rfrac))                                             &
  mode_incld_so2_rfrac = glomap_config%mode_incld_so2_rfrac
IF (PRESENT(l_rainout)) l_rainout = glomap_config%l_rainout
IF (PRESENT(l_cv_rainout)) l_cv_rainout = glomap_config%l_cv_rainout
IF (PRESENT(i_mode_nucscav)) i_mode_nucscav = glomap_config%i_mode_nucscav
IF (PRESENT(l_impc_scav)) l_impc_scav = glomap_config%l_impc_scav
IF (PRESENT(l_dust_slinn_impc_scav))                                           &
  l_dust_slinn_impc_scav = glomap_config%l_dust_slinn_impc_scav

! -- GLOMAP emissions configuration options --
IF (PRESENT(l_ukca_primss)) l_ukca_primss = glomap_config%l_ukca_primss
IF (PRESENT(l_ukca_primsu)) l_ukca_primsu = glomap_config%l_ukca_primsu
IF (PRESENT(l_ukca_primdu)) l_ukca_primdu = glomap_config%l_ukca_primdu
IF (PRESENT(l_ukca_primbcoc)) l_ukca_primbcoc = glomap_config%l_ukca_primbcoc
IF (PRESENT(l_ukca_prim_moc)) l_ukca_prim_moc = glomap_config%l_ukca_prim_moc
IF (PRESENT(l_bcoc_bf)) l_bcoc_bf = glomap_config%l_bcoc_bf
IF (PRESENT(l_bcoc_bm)) l_bcoc_bm = glomap_config%l_bcoc_bm
IF (PRESENT(l_bcoc_ff)) l_bcoc_ff = glomap_config%l_bcoc_ff
IF (PRESENT(l_ukca_scale_biom_aer_ems))                                        &
  l_ukca_scale_biom_aer_ems = glomap_config%l_ukca_scale_biom_aer_ems
IF (PRESENT(biom_aer_ems_scaling))                                             &
  biom_aer_ems_scaling = glomap_config%biom_aer_ems_scaling
IF (PRESENT(l_ukca_fine_no3_prod))                                             &
   l_ukca_fine_no3_prod = glomap_config%l_ukca_fine_no3_prod
IF (PRESENT(l_ukca_coarse_no3_prod))                                           &
   l_ukca_coarse_no3_prod = glomap_config%l_ukca_coarse_no3_prod
IF (PRESENT(l_no3_prod_in_aero_step))                                          &
   l_no3_prod_in_aero_step = glomap_config%l_no3_prod_in_aero_step
IF (PRESENT(hno3_uptake_coeff))                                                &
   hno3_uptake_coeff = glomap_config%hno3_uptake_coeff
IF (PRESENT(l_ukca_scale_sea_salt_ems))                                        &
  l_ukca_scale_sea_salt_ems = glomap_config%l_ukca_scale_sea_salt_ems
IF (PRESENT(sea_salt_ems_scaling))                                             &
  sea_salt_ems_scaling = glomap_config%sea_salt_ems_scaling
IF (PRESENT(l_ukca_scale_marine_pom_ems))                                      &
  l_ukca_scale_marine_pom_ems = glomap_config%l_ukca_scale_marine_pom_ems
IF (PRESENT(marine_pom_ems_scaling))                                           &
  marine_pom_ems_scaling = glomap_config%marine_pom_ems_scaling
IF (PRESENT(i_primss_method))                                                  &
  i_primss_method = glomap_config%i_primss_method

! -- GLOMAP feedback configuration options --
IF (PRESENT(l_ukca_radaer)) l_ukca_radaer = glomap_config%l_ukca_radaer
IF (PRESENT(i_ukca_tune_bc)) i_ukca_tune_bc = glomap_config%i_ukca_tune_bc
IF (PRESENT(i_ukca_activation_scheme))                                         &
  i_ukca_activation_scheme = glomap_config%i_ukca_activation_scheme
IF (PRESENT(i_ukca_nwbins)) i_ukca_nwbins = glomap_config%i_ukca_nwbins
IF (PRESENT(sigwmin)) sigwmin = glomap_config%sigwmin
IF (PRESENT(sigma_updraught_scaling))                                          &
  sigma_updraught_scaling = glomap_config%sigma_updraught_scaling
IF (PRESENT(l_ntpreq_n_activ_sum))                                             &
  l_ntpreq_n_activ_sum = glomap_config%l_ntpreq_n_activ_sum
IF (PRESENT(l_ntpreq_dryd_nuc_sol))                                            &
  l_ntpreq_dryd_nuc_sol = glomap_config%l_ntpreq_dryd_nuc_sol
IF (PRESENT(l_ukca_sfix)) l_ukca_sfix = glomap_config%l_ukca_sfix

! -- GLOMAP temporary logicals --
IF (PRESENT(l_fix_neg_pvol_wat))                                               &
  l_fix_neg_pvol_wat = glomap_config%l_fix_neg_pvol_wat
IF (PRESENT(l_fix_ukca_impscav))                                               &
  l_fix_ukca_impscav = glomap_config%l_fix_ukca_impscav
IF (PRESENT(l_fix_nacl_density))                                               &
  l_fix_nacl_density = glomap_config%l_fix_nacl_density
IF (PRESENT(l_improve_aero_drydep))                                            &
  l_improve_aero_drydep = glomap_config%l_improve_aero_drydep
IF (PRESENT(l_fix_ukca_activate_pdf))                                          &
  l_fix_ukca_activate_pdf = glomap_config%l_fix_ukca_activate_pdf
IF (PRESENT(l_fix_ukca_activate_vert_rep))                                     &
  l_fix_ukca_activate_vert_rep = glomap_config%l_fix_ukca_activate_vert_rep
IF (PRESENT(l_bug_repro_tke_index))                                            &
  l_bug_repro_tke_index = glomap_config%l_bug_repro_tke_index
IF (PRESENT(l_fix_ukca_hygroscopicities))                                      &
  l_fix_ukca_hygroscopicities = glomap_config%l_fix_ukca_hygroscopicities
IF (PRESENT(l_fix_ukca_water_content))                                         &
  l_fix_ukca_water_content = glomap_config%l_fix_ukca_water_content

! -- GLOMAP internal configuration variables --
IF (PRESENT(n_dust_emissions)) n_dust_emissions = glomap_config%n_dust_emissions
IF (PRESENT(i_dust_scheme)) i_dust_scheme = glomap_config%i_dust_scheme
IF (PRESENT(l_6bin_dust_no3)) l_6bin_dust_no3 = glomap_config%l_6bin_dust_no3
IF (PRESENT(l_2bin_dust_no3)) l_2bin_dust_no3 = glomap_config%l_2bin_dust_no3

RETURN
END SUBROUTINE ukca_get_config

! ----------------------------------------------------------------------
SUBROUTINE copy_config_vector_integer(vec_in, vec_out)
! ----------------------------------------------------------------------
! Description:
!   Returns a copy of an integer array valued configuration variable.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, ALLOCATABLE, INTENT(IN) :: vec_in(:)
INTEGER, ALLOCATABLE, INTENT(OUT) :: vec_out(:)

IF (ALLOCATED(vec_in)) THEN
  IF (ALLOCATED(vec_out)) THEN
    IF (SIZE(vec_out) /= SIZE(vec_in)) DEALLOCATE(vec_out)
  END IF
  IF (.NOT. ALLOCATED(vec_out)) ALLOCATE(vec_out(SIZE(vec_in)))
  vec_out = vec_in
ELSE
  IF (ALLOCATED(vec_out)) DEALLOCATE(vec_out)
END IF

RETURN
END SUBROUTINE copy_config_vector_integer

END MODULE ukca_config_specification_mod
