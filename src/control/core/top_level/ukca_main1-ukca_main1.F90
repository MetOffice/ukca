! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module containing the main subroutine for performing a UKCA time step
!    (called from the UKCA API routine 'ukca_step').
!
!  Method:
!  1) Tracers and non-transported prognostic fields are passed as arguments
!     to 'ukca_main1' with intent IN OUT.
!     Other input prognostic fields, referred to as environmental drivers are
!     provided as public variables in a separate UKCA module
!     'ukca_environment_fields_mod'.
!     Emissions fields (a special case of environmental drivers) are provided
!     in a data structure 'ncdf_emissions' in UKCA module
!     'ukca_emiss_struct_mod'.
!     Configuration variables are provided in data structures
!     'ukca_config' and 'glomap_config' in the UKCA module
!     'ukca_config_specification_mod'.
!  2) UKCA routines are called depending on the scheme selected:
!     - Emissions control routine
!     - Chemistry control routine
!     - Aerosol control routine for GLOMAP-mode
!     - Activate routine for CDNC calculation
!
! Part of the UKCA model, a community model supported by The Met Office
! and NCAS, with components initially provided by The University of
! Cambridge, University of Leeds and The Met Office. See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
MODULE ukca_main1_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: ukca_main1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_MAIN1_MOD'

CONTAINS
!
! Subroutine Interface:
! ----------------------------------------------------------------------
SUBROUTINE ukca_main1(error_code_ptr, timestep_number, current_time,           &
                      environ_ptrs, r_theta_levels, r_rho_levels,              &
                      diagnostics, all_tracers, all_ntp,                       &
                      previous_time, error_message, error_routine)
! ----------------------------------------------------------------------

USE missing_data_mod, ONLY: rmdi

USE ukca_config_specification_mod, ONLY: ukca_config, glomap_config,           &
                                         i_ukca_activation_arg,                &
                                         l_ukca_config_available,              &
                                         int_method_nr, int_method_be_explicit,&
                                         glomap_variables

USE ukca_environment_req_mod, ONLY: l_environ_req_available
USE ukca_environment_mod, ONLY: clear_environment_fields
USE ukca_environment_fields_mod, ONLY:                                         &
    env_field_ptrs_type,                                                       &
    land_points, land_index,                                                   &
    sin_declination,        equation_of_time,                                  &
    conv_cloud_base,        conv_cloud_top,       kent,                        &
    kent_dsc,               land_sea_mask,        dust_flux,                   &
    theta,                  soil_moisture_layer1, q,                           &
    qcf,                    Tstar,                                             &
    zbl,                    Rough_length,         seaice_frac,                 &
    o3_offline,             chloro_sea,                                        &
    so4_sa_clim,            so4_aitken,           so4_accum,                   &
    soot_fresh,             soot_aged,                                         &
    ocff_fresh,             ocff_aged,            dms_sea_conc,                &
    vertvel,                frac_types,                                        &
    l_tile_active,          laift_lp,             canhtft_lp,                  &
    Tstar_tile,             z0tile_lp,                                         &
    rho_r2,                 qcl,                  exner_rho_levels,            &
    cloud_frac,             cloud_liq_frac,                                    &
    biogenic,               fland,                                             &
    exner_theta_levels,     p_rho_levels,         p_theta_levels,              &
    pstar,                                                                     &
    sea_salt_film,          sea_salt_jet,                                      &
    rhokh_rdz,              dtrdz,                we_lim,                      &
    t_frac,                 zrzi,                 we_lim_dsc,                  &
    t_frac_dsc,             zrzi_dsc,             zhsc,                        &
    U_scalar_10m,           surf_hf,              stcon,                       &
    u_s,                    bl_tke,               h2o2_offline,                &
    ls_rain3d,              ls_snow3d,            rim_cry,                     &
    rim_agg,                autoconv,             accretion,                   &
    conv_rain3d,            conv_snow3d,          ch4_wetl_emiss,              &
    pv_on_theta_mlevs,      latitude,             longitude,                   &
    sin_latitude,           cos_latitude,         tan_latitude,                &
    dust_div1,              dust_div2,            dust_div3,                   &
    dust_div4,              dust_div5,            dust_div6,                   &
    interf_z,               grid_surf_area,       grid_area_fullht,            &
    grid_volume,            photol_rates,         ext_cg_flash,                &
    ext_ic_flash,           co2_interactive,      grid_airmass

USE ukca_pr_inputs_mod, ONLY: ukca_pr_inputs

USE ukca_diagnostics_type_mod, ONLY: diagnostics_type
USE ukca_diagnostics_output_mod, ONLY: update_skipped_diag_flags,              &
                                       blank_out_missing_diags

!!!! Note: LFRIC-specific pre-processor directives used in this module are
!!!! inappropriate in UKCA and should be removed but must be retained while
!!!! LFRic uses the UM version of ukca_um_legacy_mod which does not contain
!!!! dummy versions of timer and autotune routines and holds a .TRUE. value
!!!! of l_um_infrastruture which allows the possible use of rh_z_top_theta and
!!!! a_realhd which may not be initialised.

#if !defined(LFRIC)
USE ukca_um_legacy_mod, ONLY:                                                  &
    timer,                                                                     &
    l_autotune_segments,                                                       &
    autotune_type,                                                             &
    autotune_init,                                                             &
    autotune_entry,                                                            &
    autotune_return,                                                           &
    autotune_start_region,                                                     &
    autotune_stop_region,                                                      &
    rh_z_top_theta, a_realhd,                                                  &
    delta_lambda, delta_phi,                                                   &
    l_um_atmos_ukca_humidity, atmos_ukca_humidity, rhcrit, lsp_qclear,         &
    rad_ait, rad_acc, chi, sigma, pi,                                          &
    stashwork34, stashwork38, stashwork50,                                     &
    copydiag, copydiag_3d, stashcode_glomap_sec, len_stlist, stindex,          &
    stlist, num_stash_levels, stash_levels, si, sf, si_last,                   &
    l_ukca_stratflux, n_strat_fluxdiags,                                       &
    l_ukca_mode_diags, n_mode_diags,                                           &
    UkcaD1Codes, istrat_first,                                                 &
    imode_first, item1_mode_diags,                                             &
    UKCA_diag_sect,                                                            &
    item1_nitrate_diags, itemN_nitrate_diags,                                  &
    item1_dust3mode_diags, itemN_dust3mode_diags,                              &
    mype,                                                                      &
    gg => g, planet_radius
#else
USE ukca_um_legacy_mod, ONLY:                                                  &
    delta_lambda, delta_phi,                                                   &
    l_um_atmos_ukca_humidity, atmos_ukca_humidity, rhcrit, lsp_qclear,         &
    rad_ait, rad_acc, chi, sigma, pi,                                          &
    stashwork34, stashwork38, stashwork50,                                     &
    copydiag, copydiag_3d, stashcode_glomap_sec, len_stlist, stindex,          &
    stlist, num_stash_levels, stash_levels, si, sf, si_last,                   &
    L_ukca_stratflux, n_strat_fluxdiags,                                       &
    l_ukca_mode_diags, n_mode_diags,                                           &
    UkcaD1Codes, istrat_first,                                                 &
    imode_first, item1_mode_diags,                                             &
    UKCA_diag_sect,                                                            &
    item1_nitrate_diags, itemN_nitrate_diags,                                  &
    item1_dust3mode_diags, itemN_dust3mode_diags,                              &
    mype,                                                                      &
    gg => g, planet_radius
#endif

USE ukca_humidity_mod,      ONLY: ukca_vmrsat_liq, ukca_vmr_clear_sky

USE asad_mod,               ONLY: nnaf, interval, advt, jpctr, jpspec, jpdd,   &
                                  jpdw, jppj

USE ukca_tracer_vars,       ONLY: trmol_post_chem
USE ukca_cspecies,          ONLY: c_species, c_na_species, n_bro, n_h2o,       &
                                  n_hcl, n_no2, n_o3,                          &
                                  ukca_calc_cspecies, n_passive
USE chemistry_constants_mod, ONLY: boltzmann, avogadro, rho_so4
USE ukca_tropopause,        ONLY:                                              &
    p_tropopause,           theta_trop,            pv_trop,                    &
    tropopause_level,       l_stratosphere,        ukca_calc_tropopause

USE ukca_mode_verbose_mod,  ONLY: glob_verbose
USE ukca_chem_offline,      ONLY: o3_offline_diag, oh_offline_diag,            &
                                  no3_offline_diag, ho2_offline_diag
USE ukca_offline_oxidants_diags_mod, ONLY: ukca_offline_oxidants_diags

USE ukca_config_defs_mod,   ONLY:                                              &
    n_chem_tracers,         n_aero_tracers,         em_chem_spec,              &
    n_boundary_vals,        lbc_mmr,                n_mode_tracers,            &
    lbc_spec,               n_use_tracers,          nmax_mode_diags,           &
    nmax_strat_fluxdiags

USE ukca_transform_halogen_mod,                                                &
                            ONLY: ukca_transform_halogen
USE ukca_mode_tracer_maps_mod,                                                 &
                            ONLY: ukca_aero_tracer_init

USE asad_chem_flux_diags,   ONLY:                                              &
    asad_allocate_chemdiag,     asad_tendency_ste,                             &
    asad_mass_diagnostic,       asad_output_tracer,                            &
    calculate_STE,                                                             &
    calculate_tendency,         L_asad_use_chem_diags,                         &
    L_asad_use_STE,             L_asad_use_tendency,                           &
    L_asad_use_mass_diagnostic, L_asad_use_output_tracer,                      &
    L_asad_use_trop_mask,       asad_tropospheric_mask

USE ukca_diags_output_ctl_mod, ONLY: ukca_diags_output_ctl
USE asad_diags_output_ctl_mod, ONLY: asad_diags_output_ctl

USE ukca_diurnal_oxidant,   ONLY: ukca_int_cosz, dealloc_diurnal_oxidant
USE ukca_emiss_ctl_mod,     ONLY: ukca_emiss_ctl
USE ukca_chemistry_setup_mod, ONLY: ukca_chemistry_setup
USE ukca_chemistry_cleanup_mod, ONLY: ukca_chemistry_cleanup
USE ukca_chemistry_ctl_be_mod, ONLY: ukca_chemistry_ctl_be
USE ukca_chemistry_ctl_tropraq_mod, ONLY: ukca_chemistry_ctl_tropraq

USE ukca_scenario_ctl_mod,  ONLY: ukca_scenario_ctl
USE ukca_chem_diags_mod,    ONLY: ukca_chem_diags
USE ukca_chemistry_ctl_mod, ONLY: ukca_chemistry_ctl
USE ukca_chemistry_ctl_col_mod, ONLY: ukca_chemistry_ctl_col
USE ukca_chemistry_ctl_full_mod, ONLY: ukca_chemistry_ctl_full
USE ukca_activate_mod,      ONLY: ukca_activate
USE ukca_aero_ctl_mod,      ONLY: ukca_aero_ctl
USE ukca_mode_diags_mod,    ONLY: ukca_mode_diags_alloc, ukca_mode_diags,      &
                                  l_ukca_cmip6_diags,                          &
                                  l_ukca_pm_diags
USE ukca_age_air_mod,       ONLY: ukca_age_air
USE level_heights_mod,      ONLY: eta_theta_levels
USE ereport_mod,            ONLY: ereport
USE missing_data_mod,       ONLY: imdi
USE ukca_calc_cloud_ph_mod, ONLY: ukca_calc_cloud_ph
USE ukca_constants,         ONLY: H_plus

USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper,          &
                      PrStatus_Min

USE ukca_tracers_mod,       ONLY: all_tracers_names, n_tracers
USE ukca_ntp_mod,           ONLY: ntp_type, name2ntpindex

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim
USE ukca_iniasad_mod,       ONLY: ukca_iniasad, ukca_iniasad_spatial_vars,     &
                                  ukca_delasad_spatial_vars

USE ukca_solang_mod,        ONLY: ukca_solang

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_ukca_uninit, errcode_env_req_uninit,         &
                          errcode_env_field_missing,                           &
                          errcode_env_field_mismatch,                          &
                          errcode_value_missing

USE ukca_environment_check_mod, ONLY: check_environment

USE ukca_time_mod,   ONLY: i_year, i_month, i_day, i_hour, i_minute,           &
                           i_day_number, i_hour_previous,                      &
                           i_minute_previous, i_second_previous,               &
                           set_time, set_previous_time

IMPLICIT NONE

! ----------------------------------------------------------------------
! Declarations
! ----------------------------------------------------------------------

! Error code for status reporting
INTEGER, POINTER, INTENT(IN) :: error_code_ptr

! Model timestep number (counted from basis time at start of run)
INTEGER, INTENT(IN) :: timestep_number

! Current model time (year, month, day, hour, minute, second, day of year)
INTEGER, INTENT(IN) :: current_time(7)

! Environment field pointers
TYPE(env_field_ptrs_type), INTENT(IN) :: environ_ptrs(:)

! Height of theta and rho levels from Earth centre
REAL, INTENT(IN) :: r_theta_levels(:,:,0:), r_rho_levels(:,:,:)

! Diagnostic request info and pointers to parent arrays for diagnostic output
TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics

! UKCA tracers. Dimensions: X,Y,Z,N
! where X is row length of tracer field (= no. of columns)
!       Y is no. of rows in tracer field
!       Z is no. of levels in tracer fields
!       N is number of tracers
REAL, INTENT(IN OUT) :: all_tracers(:, :, :, :)

! Non-transported prognostics.
TYPE(ntp_type), INTENT(IN OUT) :: all_ntp(:)

! Model time at previous timestep (required for chemistry)
INTEGER, OPTIONAL, INTENT(IN) :: previous_time(7)

! Optional arguments for status reporting
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local scalars
INTEGER    :: row_length        ! Size of UKCA x dimension (columns)
INTEGER    :: rows              ! Size of UKCA y dimension (rows)
INTEGER    :: model_levels      ! Size of UKCA z dimension (levels)
INTEGER    :: theta_field_size  ! No. of points in horizontal plane
INTEGER    :: tot_n_pnts        ! No. of points in full domain
INTEGER    :: n_pnts            ! No. of points passed to ASAD
INTEGER    :: section           ! stash section
INTEGER    :: item              ! stash item
INTEGER    :: i                 ! loop variables
INTEGER    :: j                 ! loop variables
INTEGER    :: k,l,n             ! loop variables
INTEGER    :: kk                ! loop counter
INTEGER    :: icnt              ! counter
INTEGER    :: icode=0           ! local error status
INTEGER, PARAMETER :: im_index = 1  ! internal model index for STASH copy calls
INTEGER    :: n_fld_present     ! No. of required environment fields present
INTEGER    :: n_fld_missing     ! No. of required environment fields missing
INTEGER, SAVE :: wetox_in_aer   ! set for wet oxidation in MODE (=1)
                                ! or UKCA (=0)
INTEGER, SAVE :: uph2so4inaer   ! update H2SO4 tracer in MODE (=1)
                                ! or UKCA (=0) depending on chemistry
INTEGER, SAVE :: k_be_top       ! top level for offline oxidants (BE)
INTEGER    :: nmax_diags_inc_nitr ! number of nitrate diagnostics
INTEGER    :: n_nitrate_diags   ! number of nitrate diagnostics
INTEGER    :: n_sup_dust_diags  ! Number of diagnostics for dust 3rd mode

REAL       :: r_minute                  ! real equiv of i_minute
REAL       :: secondssincemidnight      ! day time
REAL       :: ssmn_incr                 ! sec. since midnight incr. counter

! Saved variables for when persistence is off
! Needed for preservation of int_zenith_angle
REAL, SAVE :: sindec_sav
REAL, SAVE :: secondssincemidnight_sav
REAL, SAVE :: eq_time_sav

LOGICAL       :: l_classic_aerosols     ! True if CLASSIC aerosols are modelled

! Local arrays
! Cache blocking feature : here sized to a maximum but within aero_ctl will be
! set to match the case in progress
INTEGER :: stride_seg   ! row_length*rows for striding through by column
INTEGER :: nseg         ! the number of segments on this MPI
INTEGER :: nbox_s(ukca_config%row_length * ukca_config%rows)
                        ! the number of boxes on each segment
INTEGER :: ncol_s(ukca_config%row_length * ukca_config%rows)
                        ! number of columns on each segment
INTEGER :: seg_rem      ! used to establish if more segments needed
INTEGER :: lbase(ukca_config%row_length * ukca_config%rows)
                        ! start point of segment in the 3d array
INTEGER :: nbox_max     ! largest chunk has highest number of boxes
INTEGER :: ik,lb        ! local loop iterators

! Tile index etc
INTEGER, ALLOCATABLE :: tile_index(:,:)    ! Indices of active tiles
INTEGER, ALLOCATABLE :: tile_pts(:)        ! No of tiles of each type

! Fields for humidity calculations when using UM routines
REAL, ALLOCATABLE :: qsatmr_ice_wat(:,:,:) ! qsat mixing ratio w.r.t. ice/water
REAL, ALLOCATABLE :: qsatmr_wat(:,:,:)     ! qsat mixing ratio w.r.t. water
REAL, ALLOCATABLE :: q1d(:)
REAL, ALLOCATABLE :: qsatmr_ice_wat1d(:)
REAL, ALLOCATABLE :: qsatmr_wat1d(:)
REAL, ALLOCATABLE :: qcf1d(:)
REAL, ALLOCATABLE :: cloud_liq_frac1d(:)
REAL, ALLOCATABLE :: cloud_frac1d(:)
REAL, ALLOCATABLE :: rhcrit1d(:)
REAL, ALLOCATABLE :: q_clr1d(:)
REAL, ALLOCATABLE :: rel_humid_frac_clr1d(:)

REAL, ALLOCATABLE :: rel_humid_frac_clr(:,:,:)
REAL, ALLOCATABLE :: ls_ppn3d(:,:,:)
REAL, ALLOCATABLE :: conv_ppn3d(:,:,:)
REAL, ALLOCATABLE :: delso2_wet_h2o2(:,:,:)
REAL, ALLOCATABLE :: delso2_wet_o3(:,:,:)
REAL, ALLOCATABLE :: delh2so4_chem(:,:,:)  ! Rate of net chem Pdn H2SO4
REAL, ALLOCATABLE :: delso2_drydep(:,:,:)
REAL, ALLOCATABLE :: delso2_wetdep(:,:,:)
REAL, ALLOCATABLE :: env_ozone3d(:,:,:)
REAL, ALLOCATABLE   :: v_latitude(:,:)     ! boundary lat
REAL, ALLOCATABLE   :: cos_zenith_angle(:,:)
REAL, ALLOCATABLE, SAVE   :: int_zenith_angle(:,:) ! integral over sza
REAL, ALLOCATABLE :: trmol_post_atmstep(:,:,:,:)
REAL, ALLOCATABLE :: mode_diags(:,:,:,:)
REAL, ALLOCATABLE :: strat_fluxdiags(:,:,:,:)
REAL, ALLOCATABLE :: totnodens(:,:,:)  ! density in molecs/m^3
REAL, ALLOCATABLE :: rel_humid_frac(:,:,:)
REAL, ALLOCATABLE :: water_vapour_mr_sat(:,:,:)
                              ! Water vapour saturation mixing ratio (kg/kg)
REAL, ALLOCATABLE :: water_vapour_mr(:,:,:)
                              ! Water vapour mixing ratio (kg/kg)
REAL, ALLOCATABLE :: water_vapour_mr_clr(:,:,:)
                              ! Clear sky water vapour mixing ratio (kg/kg)

! sat vap pressure with respect to liquid water irrespective of temperature
REAL, ALLOCATABLE :: qsvp(:,:,:)

REAL, SAVE, ALLOCATABLE :: z_half(:,:,:)
REAL, SAVE, ALLOCATABLE :: z_half_alllevs(:,:,:)

! 3D array for calulated pH values to use in chem_ctl
REAL, ALLOCATABLE :: H_plus_3d_arr(:,:,:)

! Plume height from explosive eruptions (m)
! Output from Plumeria model in ukca_volcanic_so2 subroutine
REAL, ALLOCATABLE :: plumeria_height(:,:)

! diagnostics from chemistry_ctl

! Nitric acid trihydrate (kg(nat)/kg(air))
REAL :: nat_psc(                                                               &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Trop CH4 burden in moles
REAL :: trop_ch4_mol(                                                          &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Trop O3 burden in moles
REAL :: trop_o3_mol(                                                           &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Trop OH burden in moles
REAL :: trop_oh_mol(                                                           &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Strat CH4 burden in moles
REAL :: strat_ch4_mol(                                                         &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Stratospheric CH4 loss (Moles/s)
REAL :: strat_ch4loss(                                                         &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Atmospheric Burden of CH4 in moles
REAL :: atm_ch4_mol(                                                           &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Atmospheric Burden of CO in moles
REAL :: atm_co_mol(                                                            &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Atmospheric Burden of Nitrous Oxide (N2O) in moles
REAL :: atm_n2o_mol(                                                           &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Atmospheric Burden of CFC-12 in moles
REAL :: atm_cf2cl2_mol(                                                        &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Atmospheric Burden of CFC-11 in moles
REAL :: atm_cfcl3_mol(                                                         &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Atmospheric Burden of CH3Br in moles
REAL :: atm_mebr_mol(                                                          &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Atmospheric Burden of H2 in moles
REAL :: atm_h2_mol(                                                            &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Derived data arrays, allocatable if necessary
REAL :: land_fraction(ukca_config%row_length, ukca_config%rows)

! Temperature on theta levels
REAL, ALLOCATABLE    :: t_theta_levels(:,:,:)

! thickness of BL layers in metres
REAL, ALLOCATABLE :: thick_bl_levels(:,:,:)

REAL, ALLOCATABLE :: p_layer_boundaries(:,:,:)
REAL, ALLOCATABLE :: so4_sa(:,:,:) ! aerosol surface area
REAL, ALLOCATABLE :: t_chem(:,:,:) ! Temperature for chemistry
REAL, ALLOCATABLE :: q_chem(:,:,:) ! Specific humidity for chemistry

REAL :: z_top_of_model      ! top of model
REAL :: z_top_of_model_ext  ! top of model (externally prescribed value)
REAL :: tol                 ! tolerance level for mismatch

! weighted cdnc = total cdnc * cldflg [m-3]
REAL :: cdncflag(                                                              &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Total activated aerosols [m-3]
REAL :: n_activ_sum(                                                           &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! weighted cdnc = total cdnc * liq_cloud_frac [m-3]
REAL :: cdncwt(                                                                &
          ukca_config%row_length,                                              &
          ukca_config%rows,                                                    &
          ukca_config%model_levels)

! Variables hoisted out of ukca_chemistry_ctl
INTEGER, SAVE :: istore_h2so4  ! location of H2SO4 in f array
REAL :: zdryrt(ukca_config%row_length,ukca_config%rows,jpdd)  ! dry dep rate
REAL :: zwetrt(ukca_config%row_length,ukca_config%rows,                        &
               ukca_config%model_levels,jpdw)                 ! wet dep rate
REAL :: shno3_3d(ukca_config%row_length,ukca_config%rows,                      &
                 ukca_config%model_levels)
INTEGER :: nlev_with_ddep(ukca_config%row_length,                              &
                          ukca_config%rows)  ! No. of levels in boundary layer

REAL, SAVE :: lambda_aitken, lambda_accum ! parameters for computation
                                          ! of surface area density
LOGICAL :: do_chemistry
LOGICAL, SAVE :: l_firstchem = .TRUE.     ! Logical for any operations
                   ! specific to first chemical timestep, which is usually
                   ! different to the first ukca_main call. This should
                   ! only be used within 'do_chemistry' sections

! required for ASAD Flux Diagnostics
INTEGER :: ierr, stashsize
REAL, ALLOCATABLE :: fluxdiag_all_tracers(:,:,:,:)

! Min value to assign SO4 used in FastJX if no CLASSIC
REAL, PARAMETER :: min_SO4_val = 1.0e-25

LOGICAL, SAVE :: l_first_call       = .TRUE.   ! true only on 1st call

LOGICAL :: l_store_value           ! T to store value, F to make difference

LOGICAL, PARAMETER :: unlump_species = .TRUE.
LOGICAL, PARAMETER :: lump_species = .FALSE.

! Bounds for relative humidity calculation
REAL, PARAMETER :: rh_min = 0.0  ! RH lower limit
REAL, PARAMETER :: rh_max = 1.0  ! RH upper limit

! Flags to track what humidity-related fields are required
LOGICAL :: l_using_rh      ! True if using relative humidity
LOGICAL :: l_using_rh_clr  ! True if using clear sky relative humidity
LOGICAL :: l_using_svp     ! True if using saturation vapour pressure of water

! Flags to track what humidity-related fields need to be calculated
LOGICAL :: l_calc_rh       ! True to calculate relative humidity
LOGICAL :: l_calc_rh_clr   ! True to calculate clear sky relative humidity
LOGICAL :: l_calc_svp      ! True to calculate saturation vapour pressure of
                           ! water

! Mask to limit formation of Nat below specified height
LOGICAL, SAVE, ALLOCATABLE  :: have_nat3d(:,:,:)
LOGICAL                     :: level_above_ht
REAL,PARAMETER :: nat_limit_ht = 1000.0  ! 1km
REAL           :: h_atmos   ! depth of atmosphere

CHARACTER (LEN=11) :: emiss_input  ! 'Ancillary  ' or  'NetCDF file'

REAL :: factor  ! Temporary working factor

! ErrorStatus
INTEGER                    :: errcode=0     ! Error flag (0 = OK)
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

#if !defined(LFRIC)
!Automatic segment size tuning
TYPE(autotune_type), ALLOCATABLE, SAVE :: autotune_state
#endif

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MAIN1'

!- End of header
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! 1. Initial set up
! ----------------------------------------------------------------------

! Set defaults for output arguments
error_code_ptr = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Create local copies of some frequently used configuration variables
row_length = ukca_config%row_length
rows = ukca_config%rows
model_levels = ukca_config%model_levels

theta_field_size = row_length * rows
tot_n_pnts = theta_field_size * model_levels
IF (ukca_config%l_ukca_asad_full) THEN
  n_pnts = tot_n_pnts
ELSE IF (ukca_config%l_ukca_asad_columns) THEN
  n_pnts = ukca_config%ukca_chem_seg_size
ELSE
  n_pnts = theta_field_size
END IF

! On first call, check for availability of UKCA configuration data
IF (l_first_call) THEN
  IF (.NOT. l_ukca_config_available) THEN
    error_code_ptr = errcode_ukca_uninit
    IF (PRESENT(error_message))                                                &
      error_message = 'No UKCA configuration has been set up'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
  IF (.NOT. l_environ_req_available) THEN
    error_code_ptr = errcode_env_req_uninit
    IF (PRESENT(error_message))                                                &
      error_message = 'No environment field requirement has been set up'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Set internal time values based on arguments

CALL set_time(current_time)

IF (ukca_config%l_ukca_chem) THEN
  IF (PRESENT(previous_time)) THEN
    CALL set_previous_time(previous_time)
  ELSE
    error_code_ptr = errcode_value_missing
    IF (PRESENT(error_message))                                                &
      error_message = 'Missing previous time argument required for chemistry'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Check whether all required environmental driver fields are present
! (None are required if the run does not use chemistry)
IF (ukca_config%l_ukca_chem) THEN
  CALL check_environment(n_fld_present, n_fld_missing)
  IF (n_fld_missing > 0) THEN
    error_code_ptr = errcode_env_field_missing
    IF (PRESENT(error_message))                                                &
      WRITE(error_message,'(A,I4,A,I4,A)') 'Missing ', n_fld_missing, ' of ',  &
      n_fld_present + n_fld_missing, ' environmental driver fields'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Ensure unsupported diagnostic output is suppressed if not using the
! UM STASH system
IF (.NOT. ukca_config%l_enable_diag_um) THEN
  l_ukca_stratflux = .FALSE.
  l_ukca_mode_diags = .FALSE.
  l_ukca_cmip6_diags = .FALSE.
  l_ukca_pm_diags = .FALSE.
  n_strat_fluxdiags = 0
  n_mode_diags = 0
  ! Allocate zero size STASH workspace for use in argument lists
  IF (.NOT. ALLOCATED(stashwork34)) ALLOCATE(stashwork34(0))
  IF (.NOT. ALLOCATED(stashwork38)) ALLOCATE(stashwork38(0))
  IF (.NOT. ALLOCATED(stashwork50)) ALLOCATE(stashwork50(0))
END IF

#if !defined(LFRIC)
!Set up automatic segment tuning
IF (l_autotune_segments) THEN
  IF (.NOT. ALLOCATED(autotune_state)) THEN
    ALLOCATE(autotune_state)
    CALL autotune_init(                                                        &
      autotune_state,                                                          &
      region_name    = 'ukca_mode_seg_size',                                   &
      tag            = 'UKCA-AERO',                                            &
      start_size     = glomap_config%ukca_mode_seg_size)
  END IF

  CALL autotune_entry(autotune_state, glomap_config%ukca_mode_seg_size)
END IF
#endif

! Calculate height above sea level at top of model.
z_top_of_model = r_theta_levels(1,1,model_levels) - planet_radius
#if !defined(LFRIC)
! If running in the UM, z_top_of_model may be overriden by an external
! UM-provided value to allow bit-comparability with previous results to be
! conserved. However, a tolerance check must then be passed to ensure
! that the UM value supplied is consistent with the value calculated above
! which UKCA takes to be the definitive value.
tol = 1e-10 * (r_theta_levels(1,1,model_levels) - r_theta_levels(1,1,1))
IF (ukca_config%l_environ_z_top) THEN
  z_top_of_model_ext = a_realhd(rh_z_top_theta)
  IF (ABS(z_top_of_model_ext - z_top_of_model) > tol) THEN
    error_code_ptr = errcode_env_field_mismatch
    IF (PRESENT(error_message))                                                &
      WRITE(error_message,'(A,E15.6,A)')                                       &
        'Difference between z_top_of_model prescribed and expected value = ',  &
        z_top_of_model_ext - z_top_of_model, ', exceeds tolerance'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
  z_top_of_model = z_top_of_model_ext
END IF
#endif

r_minute  = REAL(i_minute)     ! Real copy of i_minute
! Set Verbosity level for GLOMAP routines
glob_verbose = 0
IF ( PrintStatus > PrStatus_Oper ) glob_verbose = 2

IF (l_first_call) THEN

  IF (PrintStatus >= PrStatus_Oper) THEN
    WRITE(umMessage,'(A)') 'First call to ukca_main'
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E14.2)') 'timestep = ',ukca_config%timestep
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8)') 'timestep number = ',timestep_number
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  IF (ukca_config%l_ukca_chem) THEN

    CALL ukca_iniasad(n_pnts)

    ! Check that some CLASSIC aerosol types/processes are modelled if
    ! heterogeneous reactions on CLASSIC aerosols are used in UKCA.

    l_classic_aerosols = ukca_config%l_use_classic_so4 .OR.                    &
                         ukca_config%l_use_classic_soot .OR.                   &
                         ukca_config%l_use_classic_ocff .OR.                   &
                         ukca_config%l_use_classic_biogenic .OR.               &
                         ukca_config%l_use_classic_seasalt

    IF (ukca_config%l_ukca_classic_hetchem .AND. .NOT. l_classic_aerosols) THEN
      cmessage = 'Cannot use heterogeneous reactions on CLASSIC aerosols ' //  &
                 'if these are not being modelled'
      errcode  = 5
      CALL ereport (RoutineName, errcode, cmessage)
    END IF

    ! Heterogeneous reactions on CLASSIC aerosols are only allowed for the
    ! RAQ chemistry scheme, although this can be extended in the future.
    IF (ukca_config%l_ukca_classic_hetchem .AND. .NOT. ukca_config%l_ukca_raq) &
      THEN
      cmessage = 'Heterogeneous chemistry on CLASSIC aerosols is a '       //  &
                 'feature only available for the RAQ chemistry scheme'
      errcode  = 6
      CALL ereport (RoutineName, errcode, cmessage)
    END IF

    ALLOCATE(c_species(n_chem_tracers+n_aero_tracers))
    ALLOCATE(c_na_species(jpspec-jpctr))
    CALL ukca_calc_cspecies()

    errcode = 0

    DO i=1,jpctr        ! to allow for case when jpctr < n_chem_tracers...
      IF (c_species(i) < 0.001) THEN
        ! check that ratio of molec. wts. array has been initialised.
        errcode = errcode +1
        icode = -1*errcode
        cmessage=' c_species array is zero for '//advt(i)
        CALL ereport(RoutineName,icode,cmessage)
      END IF
    END DO

    IF (errcode > 0) THEN
      cmessage=' c_species array has zero values'//                            &
                   ', check UKCA_CSPECIES routine'
      CALL ereport(RoutineName,errcode,cmessage)
    END IF

    DO i=1,nnaf
      IF (c_na_species(i) < 0.001) THEN
        ! check that ratio of molec. wts. array has been initialised.
        cmessage=' c_na_species array contains zero value'
        errcode = i
        CALL ereport(RoutineName,errcode,cmessage)
      END IF
    END DO
    ! L_ukca_advh2o *MUST* be .TRUE. if using H2O as 'TR'
    IF ((n_h2o > 0) .AND. (.NOT. ukca_config%l_ukca_advh2o)) THEN
      cmessage=                                                                &
           ' H2O is defined as a tracer, but L_UKCA_ADVH2O=F'
      errcode = n_h2o
      CALL ereport(RoutineName,errcode,cmessage)
    END IF

    ! Ensure that particulate fraction of SO2 emissions is between 0 and 5%.
    IF (.NOT. ukca_config%l_ukca_emissions_off) THEN
      IF ( ( ANY(em_chem_spec == 'SO2_high  ') .OR.                            &
             ANY(em_chem_spec == 'SO2_nat   ') .OR.                            &
             ANY(em_chem_spec == 'SO2_low   ') ) .AND.                         &
           ( ukca_config%mode_parfrac < 0.0 .OR.                               &
             ukca_config%mode_parfrac > 5.0) ) THEN
        WRITE(cmessage,'(A,F16.8)') 'SO2 emissions are on but ' //             &
              'mode_parfrac is invalid (should be 0.0-5.0): ',                 &
              ukca_config%mode_parfrac
        errcode = 1
        CALL ereport(RoutineName,errcode, cmessage)
      END IF
    END IF

    ! set constants needed in sulphur chemistry
    IF (ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_strattrop .OR.        &
        ukca_config%l_ukca_stratcfc .OR. ukca_config%l_ukca_cristrat) THEN
      lambda_aitken = 3.0/rho_so4 / chi /                                      &
                     rad_ait * EXP(-2.5 * (LOG(sigma)**2))                     &
                     * 0.01 ! revert from m^-1 to cm^2/cm^3
      lambda_accum  = 3.0/rho_so4 / chi /                                      &
                     rad_acc * EXP(-2.5 * (LOG(sigma)**2))                     &
                     * 0.01 ! revert from m^-2 to cm^2/cm^3
    END IF

    ! Check that solver interval has been set correctly in ukca_init
    ! This is not relevant if running without chem i.e. Age-of-air mode

    icode = 0
    IF (ukca_config%ukca_int_method == int_method_nr .OR.                      &
        ukca_config%ukca_int_method == int_method_be_explicit) THEN
      ! interval can be greater than 1 for BE explicit and N-R solvers
      IF (interval == imdi) icode = 101
    ELSE
      ! only allowed to use interval = 1 for other solvers
      IF (interval == imdi .OR. interval /= 1) icode = 102
    END IF
    IF (icode > 0) THEN
      cmessage = 'Error for solver interval'
      WRITE(umMessage,'(A40,A12,I6)') cmessage,' Interval: ',interval
      CALL umPrint(umMessage,src=RoutineName)
      CALL ereport(RoutineName, icode, cmessage)
    END IF
  END IF      ! l_ukca_chem

END IF  ! l_first_call

! allocate flux array for ASAD chemical diagnostics
IF (L_asad_use_chem_diags) THEN
  CALL asad_allocate_chemdiag(row_length,rows)
END IF

! decide whether to do chemistry
IF ( ukca_config%l_ukca_chem ) THEN
  do_chemistry = (MOD(timestep_number, interval) == 0)
ELSE
  do_chemistry = .FALSE.
END IF

! Reallocate the spatial variables in ASAD if persistence is off
IF (ukca_config%l_ukca_chem .AND. ukca_config%l_ukca_persist_off .AND.         &
    (.NOT. l_first_call)) THEN
  CALL ukca_iniasad_spatial_vars(n_pnts)
END IF

! set up nmr and mmr indices on first call, only if using
! GLOMAP-mode
IF (ukca_config%l_ukca_mode .AND. l_first_call) CALL ukca_aero_tracer_init()

! ----------------------------------------------------------------------
! 3. Calculate any derived variables
! ----------------------------------------------------------------------

! Derived variables for chemistry
IF (ukca_config%l_ukca_chem) THEN
  ALLOCATE(t_theta_levels(row_length,rows,model_levels))
  ALLOCATE(Thick_bl_levels(row_length,rows,ukca_config%bl_levels))
  ALLOCATE(p_layer_boundaries(row_length,rows,0:model_levels))
  ALLOCATE(totnodens(row_length,rows,model_levels))
  ALLOCATE(ls_ppn3d(row_length,rows,model_levels))
  ALLOCATE(conv_ppn3d(row_length,rows,model_levels))
  ALLOCATE(tile_index(land_points,ukca_config%ntype))
  ALLOCATE(tile_pts(ukca_config%ntype))
  ALLOCATE(plumeria_height(1:row_length, 1:rows))
  ALLOCATE(so4_sa(row_length, rows, model_levels)) ! sulphate area density
  ! SO2 fluxes are required in chemistry_ctl
  ALLOCATE(delSO2_wet_h2o2(row_length,rows,model_levels))
  ALLOCATE(delSO2_wet_o3(row_length,rows,model_levels))
  ALLOCATE(delh2so4_chem(row_length,rows,model_levels))
  ALLOCATE(delSO2_drydep(row_length,rows,model_levels))
  ALLOCATE(delSO2_wetdep(row_length,rows,model_levels))

  ! Initialise fluxes to zero as may not be filled by chemistry
  delSO2_wet_h2o2(:,:,:)=0.0
  delSO2_wet_o3(:,:,:)=0.0
  delh2so4_chem(:,:,:)=0.0
  delSO2_drydep(:,:,:)=0.0
  delSO2_wetdep(:,:,:)=0.0
  ! END IF ukca_chem

  ! Derived variables for stratopheric flux diagnostics
  IF (.NOT. ALLOCATED(strat_fluxdiags)) THEN
    IF (.NOT. L_ukca_stratflux) n_strat_fluxdiags=0
    IF (n_strat_fluxdiags > 0) THEN
      ALLOCATE(strat_fluxdiags(row_length,rows,model_levels,                   &
           n_strat_fluxdiags))
      strat_fluxdiags=0.0
    END IF
  END IF

  ! Required in call to UKCA_MODE, but may be unallocated
  IF (.NOT. ALLOCATED(mode_diags)) THEN
    IF (.NOT. ukca_config%l_ukca_mode) n_mode_diags=0
    ALLOCATE(mode_diags(row_length,rows,model_levels,                          &
             n_mode_diags))
    mode_diags=0.0
  END IF

  ! Required in call to UKCA_EMISS_CTL, but may be unallocated
  IF (.NOT. ukca_config%l_ukca_emissions_off) THEN
    IF (.NOT. ALLOCATED(land_index)) ALLOCATE(land_index(0))
    IF (.NOT. ALLOCATED(latitude)) THEN
      ALLOCATE(latitude(row_length,rows))
      latitude(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(longitude)) THEN
      ALLOCATE(longitude(row_length,rows))
      longitude(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(sin_latitude)) THEN
      ALLOCATE(sin_latitude(row_length,rows))
      sin_latitude(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(cos_latitude)) THEN
      ALLOCATE(cos_latitude(row_length,rows))
      cos_latitude(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(tan_latitude)) THEN
      ALLOCATE(tan_latitude(row_length,rows))
      tan_latitude(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(zbl)) THEN
      ALLOCATE(zbl(row_length,rows))
      zbl(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(ch4_wetl_emiss)) THEN
      ALLOCATE(ch4_wetl_emiss(row_length,rows))
      ch4_wetl_emiss(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(seaice_frac)) THEN
      ALLOCATE(seaice_frac(row_length,rows))
      seaice_frac(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(u_scalar_10m)) THEN
      ALLOCATE(u_scalar_10m(row_length,rows))
      u_scalar_10m(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(tstar)) THEN
      ALLOCATE(tstar(row_length,rows))
      tstar(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(dms_sea_conc)) THEN
      ALLOCATE(dms_sea_conc(row_length,rows))
      dms_sea_conc(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(chloro_sea)) THEN
      ALLOCATE(chloro_sea(row_length,rows))
      chloro_sea(:,:) = 0.0
    END IF
    IF ( .NOT. ALLOCATED(dust_flux)) THEN
      ALLOCATE(dust_flux(row_length,rows,glomap_config%n_dust_emissions))
      dust_flux(:,:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(kent)) THEN
      ALLOCATE(kent(row_length,rows))
      kent = 0
    END IF
    IF (.NOT. ALLOCATED(kent_dsc)) THEN
      ALLOCATE(kent_dsc(row_length,rows))
      kent_dsc = 0
    END IF
    IF (.NOT. ALLOCATED(zhsc)) THEN
      ALLOCATE(zhsc(row_length,rows))
      zhsc = 0.0
    END IF
    IF (.NOT. ALLOCATED(rhokh_rdz)) THEN
      ALLOCATE(rhokh_rdz(row_length,rows,ukca_config%bl_levels))
      rhokh_rdz = 0.0
    END IF
    IF (.NOT. ALLOCATED(dtrdz)) THEN
      ALLOCATE(dtrdz(row_length,rows,ukca_config%bl_levels))
      dtrdz = 0.0
    END IF
    IF (.NOT. ALLOCATED(we_lim)) THEN
      ALLOCATE(we_lim(row_length,rows,ukca_config%nlev_ent_tr_mix))
      we_lim = 0.0
    END IF
    IF (.NOT. ALLOCATED(t_frac)) THEN
      ALLOCATE(t_frac(row_length,rows,ukca_config%nlev_ent_tr_mix))
      t_frac = 0.0
    END IF
    IF (.NOT. ALLOCATED(zrzi)) THEN
      ALLOCATE(zrzi(row_length,rows,ukca_config%nlev_ent_tr_mix))
      zrzi = 0.0
    END IF
    IF (.NOT. ALLOCATED(we_lim_dsc)) THEN
      ALLOCATE(we_lim_dsc(row_length,rows,ukca_config%nlev_ent_tr_mix))
      we_lim_dsc = 0.0
    END IF
    IF (.NOT. ALLOCATED(t_frac_dsc)) THEN
      ALLOCATE(t_frac_dsc(row_length,rows,ukca_config%nlev_ent_tr_mix))
      t_frac_dsc = 0.0
    END IF
    IF (.NOT. ALLOCATED(zrzi_dsc)) THEN
      ALLOCATE(zrzi_dsc(row_length,rows,ukca_config%nlev_ent_tr_mix))
      zrzi_dsc = 0.0
    END IF
    IF (.NOT. ALLOCATED(exner_rho_levels)) THEN
      ALLOCATE(exner_rho_levels(row_length,rows,model_levels+1))
      exner_rho_levels = 0.0
    END IF
    IF (.NOT. ALLOCATED(rho_r2)) THEN
      ALLOCATE(rho_r2(row_length,rows,model_levels))
      rho_r2 = 0.0
    END IF
    IF (.NOT. ALLOCATED(grid_surf_area)) THEN
      ALLOCATE(grid_surf_area(row_length,rows))
      grid_surf_area = 0.0
    END IF
    IF (.NOT. ALLOCATED(conv_cloud_base)) THEN
      ALLOCATE(conv_cloud_base(row_length,rows))
      conv_cloud_base(:,:)=0
    END IF
    IF (.NOT. ALLOCATED(conv_cloud_top)) THEN
      ALLOCATE(conv_cloud_top(row_length,rows))
      conv_cloud_top(:,:)=0
    END IF
    IF (.NOT. ALLOCATED(ext_cg_flash)) THEN
      ALLOCATE(ext_cg_flash(row_length, rows))
      ext_cg_flash = 0.0
    END IF
    IF (.NOT. ALLOCATED(ext_ic_flash)) THEN
      ALLOCATE(ext_ic_flash(row_length, rows))
      ext_ic_flash = 0.0
    END IF
    IF (.NOT. ALLOCATED(dust_div1)) THEN
      ALLOCATE(dust_div1(row_length,rows,model_levels))
      dust_div1(:,:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(dust_div2)) THEN
      ALLOCATE(dust_div2(row_length,rows,model_levels))
      dust_div2(:,:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(dust_div3)) THEN
      ALLOCATE(dust_div3(row_length,rows,model_levels))
      dust_div3(:,:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(dust_div4)) THEN
      ALLOCATE(dust_div4(row_length,rows,model_levels))
      dust_div4(:,:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(dust_div5)) THEN
      ALLOCATE(dust_div5(row_length,rows,model_levels))
      dust_div5(:,:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(dust_div6)) THEN
      ALLOCATE(dust_div6(row_length,rows,model_levels))
      dust_div6(:,:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(grid_area_fullht)) THEN
      ALLOCATE(grid_area_fullht(row_length,rows,model_levels))
      grid_area_fullht(:,:,:) = 0.0
    END IF
  END IF  ! .NOT. ukca_config%l_ukca_emissions_off

  ! Required in calls to chemistry control routines (UKCA_CHEMISTRY_CTL etc.)
  ! but may be unallocated
  IF (do_chemistry) THEN
    IF (.NOT. ALLOCATED(land_index)) ALLOCATE(land_index(0))
    IF (.NOT. ALLOCATED(latitude)) THEN
      ALLOCATE(latitude(row_length,rows))
      latitude(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(longitude)) THEN
      ALLOCATE(longitude(row_length,rows))
      longitude(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(sin_latitude)) THEN
      ALLOCATE(sin_latitude(row_length,rows))
      sin_latitude(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(tan_latitude)) THEN
      ALLOCATE(tan_latitude(row_length,rows))
      tan_latitude(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(frac_types)) THEN
      ALLOCATE(frac_types(land_points, ukca_config%ntype))
      frac_types(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(zbl)) THEN
      ALLOCATE(zbl(row_length,rows))
      zbl(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(tstar)) THEN
      ALLOCATE(tstar(row_length,rows))
      tstar(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(rough_length)) THEN
      ALLOCATE(rough_length(row_length, rows))
      rough_length(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(u_s)) THEN
      ALLOCATE(u_s(row_length, rows))
      u_s(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(surf_hf)) THEN
      ALLOCATE(surf_hf(row_length, rows))
      surf_hf(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(seaice_frac)) THEN
      ALLOCATE(seaice_frac(row_length,rows))
      seaice_frac(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(stcon)) THEN
      ALLOCATE(stcon(row_length, rows, ukca_config%npft))
      stcon(:,:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(soil_moisture_layer1)) THEN
      ALLOCATE(soil_moisture_layer1(land_points))
      soil_moisture_layer1(:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(laift_lp)) THEN
      ALLOCATE(laift_lp(land_points, ukca_config%npft))
      laift_lp(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(canhtft_lp)) THEN
      ALLOCATE(canhtft_lp(land_points, ukca_config%npft))
      canhtft_lp(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(z0tile_lp)) THEN
      ALLOCATE(z0tile_lp(land_points, ukca_config%ntype))
      z0tile_lp(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(tstar_tile)) THEN
      ALLOCATE(tstar_tile(land_points, ukca_config%ntype))
      tstar_tile(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(photol_rates)) THEN
      ALLOCATE(photol_rates(row_length, rows, model_levels, jppj))
      photol_rates(:,:,:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(grid_volume)) THEN
      ALLOCATE(grid_volume(row_length,rows,model_levels))
      grid_volume(:,:,:) = 0.0
    END IF

  END IF
  ! Air mass passed as argument for emissions, chemistry and aerosols
  IF ( ((.NOT. ukca_config%l_ukca_emissions_off) .OR. do_chemistry .OR.        &
       l_asad_use_mass_diagnostic) .AND. .NOT. ALLOCATED(grid_airmass) ) THEN
    ALLOCATE(grid_airmass(row_length,rows,model_levels))
    grid_airmass(:,:,:) = 0.0
  END IF

  ! Required in call to UKCA_AERO_CTL, but may be unallocated
  IF (do_chemistry .AND. ukca_config%l_ukca_mode) THEN
    IF (.NOT. ALLOCATED(zbl)) THEN
      ALLOCATE(zbl(row_length,rows))
      zbl(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(seaice_frac)) THEN
      ALLOCATE(seaice_frac(row_length,rows))
      seaice_frac(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(rough_length)) THEN
      ALLOCATE(rough_length(row_length, rows))
      rough_length(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(u_s)) THEN
      ALLOCATE(u_s(row_length, rows))
      u_s(:,:) = 0.0
    END IF
    IF (.NOT. ALLOCATED(ls_rain3d)) THEN
      ALLOCATE(ls_rain3d(row_length,rows,model_levels))
      ls_rain3d(:,:,:)=0.0
    END IF
    IF (.NOT. ALLOCATED(conv_rain3d)) THEN
      ALLOCATE(conv_rain3d(row_length,rows,model_levels))
      conv_rain3d(:,:,:)=0.0
    END IF
    IF (.NOT. ALLOCATED(ls_snow3d)) THEN
      ALLOCATE(ls_snow3d(row_length,rows,model_levels))
      ls_snow3d(:,:,:)=0.0
    END IF
    IF (.NOT. ALLOCATED(conv_snow3d)) THEN
      ALLOCATE(conv_snow3d(row_length,rows,model_levels))
      conv_snow3d(:,:,:)=0.0
    END IF
    IF (.NOT. ALLOCATED(autoconv)) THEN
      ALLOCATE(autoconv(row_length,rows,model_levels))
      autoconv(:,:,:)=0.0
    END IF
    IF (.NOT. ALLOCATED(accretion)) THEN
      ALLOCATE(accretion(row_length,rows,model_levels))
      accretion(:,:,:)=0.0
    END IF
    IF (.NOT. ALLOCATED(rim_agg)) THEN
      ALLOCATE(rim_agg(row_length,rows,model_levels))
      rim_agg(:,:,:)=0.0
    END IF
    IF (.NOT. ALLOCATED(rim_cry)) THEN
      ALLOCATE(rim_cry(row_length,rows,model_levels))
      rim_cry(:,:,:)=0.0
    END IF
  END IF

  ! initialise budget diagnostics at top of routine
  nat_psc = 0.0
  trop_ch4_mol = 0.0
  trop_o3_mol = 0.0
  trop_oh_mol = 0.0
  strat_ch4_mol = 0.0
  strat_ch4loss = 0.0
  atm_ch4_mol = 0.0
  atm_co_mol = 0.0
  atm_n2o_mol = 0.0
  atm_cf2cl2_mol = 0.0
  atm_cfcl3_mol = 0.0
  atm_mebr_mol = 0.0
  atm_h2_mol = 0.0

  IF ( l_first_call .OR. ukca_config%l_ukca_persist_off ) THEN

    IF ( .NOT. ALLOCATED(p_tropopause) )                                       &
      ALLOCATE(p_tropopause(row_length,rows))
    IF ( .NOT. ALLOCATED(tropopause_level) )                                   &
      ALLOCATE(tropopause_level(row_length,rows))
    IF ( .NOT. ALLOCATED(theta_trop) )                                         &
      ALLOCATE(theta_trop(row_length,rows))
    IF ( .NOT. ALLOCATED(pv_trop) )                                            &
      ALLOCATE(pv_trop(row_length,rows))
    IF ( .NOT. ALLOCATED(L_stratosphere) )                                     &
      ALLOCATE(L_stratosphere(row_length,rows,model_levels))

    IF ( .NOT. ALLOCATED(z_half) )                                             &
      ALLOCATE(z_half(row_length,rows,ukca_config%bl_levels))
    IF ( .NOT. ALLOCATED(z_half_alllevs) )                                     &
      ALLOCATE(z_half_alllevs(row_length,rows,model_levels))

    DO k = 1,model_levels
      DO j = 1, rows
        DO i= 1, row_length
          z_half_alllevs(i,j,k)=                                               &
                       r_rho_levels(i,j,k)-r_theta_levels(i,j,0)
        END DO
      END DO
    END DO
    DO k = 1,ukca_config%bl_levels
      z_half(:,:,k) = z_half_alllevs(:,:,k)
    END DO

    ! Calculate last model level for integration of offline oxidants (BE)
    k_be_top = model_levels        ! default
    IF (ukca_config%l_ukca_offline_be) THEN
      k_loop: DO k = 2,model_levels
        IF ((MINVAL(r_theta_levels(:,:,k)) - planet_radius) >                  &
             ukca_config%max_z_for_offline_chem) THEN
          k_be_top = k
          EXIT k_loop
        END IF
      END DO k_loop
      IF (PrintStatus >=  PrStatus_Oper) THEN
        WRITE(umMessage,'(A50,I5)') 'Top level for offline oxidants'//         &
                                    ' integration: ',k_be_top
        CALL umPrint(umMessage,src=RoutineName)
      END IF
    END IF    ! l_ukca_offline_be
  END IF     ! l_first_call or l_ukca_persist_off

    ! Read and interpolate aerosol surface area density.
    ! Depending on the choice of logicals, the array can contain
    ! surface area densities based on:
    ! - Climatological aerosol in whole domain
    ! - Climatology (above 12km)+ CLASSIC/AeroClim aerosol (below 12km)
    ! - Climatology (above 12km)+ GLOMAP aerosol  (below 12km)
    ! - GLOMAP aerosol in the whole domain
  so4_sa(:,:,:) = 0.0
  cmessage = ' '
  IF (do_chemistry .AND. ukca_config%l_ukca_het_psc) THEN
    IF (ukca_config%l_ukca_sa_clim) THEN
      cmessage =                                                               &
        'UKCA_Het_PSC: Using only climatological aerosol surface area'
      so4_sa(:,:,:) = so4_sa_clim(:,:,:)
      ! Below uses SO4 tracer from classic scheme for troposphere only
      IF (ukca_config%l_use_classic_so4 .AND.                                  &
          ALLOCATED(so4_aitken) .AND. ALLOCATED(so4_accum)) THEN
        cmessage =                                                             &
         'UKCA_Het_PSC: Using Climatology + CLASSIC/climatology'//             &
         ' aerosol surface area'
        ! Calculate surface area density fields from aerosol tracers
        !   and copy into cells below 12 km
        DO k=1,model_levels
          IF (eta_theta_levels(k) * z_top_of_model < 12000.0) THEN
            so4_sa(:,:,k) = (                                                  &
             lambda_aitken * so4_aitken(:,:,k)                                 &
            +lambda_accum  * so4_accum (:,:,k))                                &
            *rho_r2(:,:,k) / (r_theta_levels(1:row_length,1:rows,k)**2)
          END IF
        END DO
      ELSE IF ( ukca_config%l_ukca_mode ) THEN
        ! Use aerosol surface area from GLOMAP for cells below 12km
        cmessage =                                                             &
         'UKCA_Het_PSC: Using Climatology + GLOMAP aerosol surface area'
        i = name2ntpindex('surfarea  ')
        DO k=1,model_levels
          IF (eta_theta_levels(k) * z_top_of_model < 12000.0)                  &
              so4_sa(:,:,k) = all_ntp(i)%data_3d(:,:,k)
        END DO
      END IF    ! l_use_classic_so4/ Mode

    ELSE        ! Not ukca_sa_clim, take aerosol fully from GLOMAP

      IF (ukca_config%l_ukca_mode ) THEN
        cmessage='UKCA_Het_PSC: Using only GLOMAP aerosol surface area '
        i = name2ntpindex('surfarea  ')
        so4_sa = all_ntp(i)%data_3d
      ELSE
        cmessage='UKCA HET_PSC: Surface area undefined. At least '//           &
          'one from l_ukca_sa_clim and l_ukca_glomap has to be selected.'
        errcode = 1
        CALL ereport(RoutineName,errcode,cmessage)
      END IF     ! ukca_config%l_ukca_mode
    END IF       ! l_ukca_sa_clim

    WHERE (so4_sa(:,:,:) < 0.0)    ! Remove any negative values
      so4_sa(:,:,:) = 0.0
    END WHERE
    ! Inform users of the choice active for this run
    IF ( l_firstchem .AND. PrintStatus > PrStatus_Min )                        &
         CALL umPrint(cmessage,src=RoutineName)

  END IF         ! do_chemistry and l_ukca_het_psc

  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        t_theta_levels(i,j,k) = exner_theta_levels(i,j,k) * theta(i,j,k)
      END DO
    END DO
  END DO

  DO j = 1, rows
    DO i = 1, row_length
      Thick_bl_levels(i,j,1) = 2.0 * (r_theta_levels(i,j,1) -                  &
                                      r_theta_levels(i,j,0))
    END DO
  END DO
  DO k = 2, ukca_config%bl_levels
    DO j = 1, rows
      DO i = 1, row_length
        Thick_bl_levels(i,j,k) = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
      END DO
    END DO
  END DO
  DO j = 1, rows
    DO i = 1, row_length
      p_layer_boundaries(i,j,0)            = pstar(i,j)
      p_layer_boundaries(i,j,model_levels) = p_theta_levels(i,j,model_levels)
    END DO
  END DO
  DO k = 1, model_levels - 1
    DO j = 1, rows
      DO i = 1, row_length
        p_layer_boundaries(i,j,k) = p_rho_levels(i,j,k+1)
      END DO
    END DO
  END DO

  ! Mass:
  ! Where gridbox airmass is not provided but is required for GLOMAP a
  ! 'relative mass' is calculated here based on the hydrostatic assumption
  ! such that
  !     mass = -b*(solid_angle/(3*g))*(r_top^3 - r_bottom^3)
  ! where   b =                                                           &
  !    &    (p_layer_boundaries(:,:,k) - p_layer_boundaries(:,:,k-1))/    &
  !    &    (r_theta_levels(:,:,k)     - r_theta_levels(:,:,k-1))         &

  ! However, here the factor solid_angle is omitted and the mass of air will
  ! then be in arbitrary units for each column.
  ! The relative mass at each level within a column is preserved as
  ! required by GLOMAP sedimentation and dry deposition processing.

  IF (ukca_config%l_ukca_mode .AND. .NOT. ukca_config%l_use_gridbox_mass) THEN
    factor = 1.0
    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          grid_airmass(i,j,k) = (-factor / (3.0*gg)) *                         &
             	 	  ((p_layer_boundaries(i,j,k) - p_layer_boundaries(i,j,k-1))/  &
                   (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)) ) *       &
                  (r_theta_levels(i,j,k)**3.0 - r_theta_levels(i,j,k-1)**3.0)
        END DO
      END DO
    END DO
  END IF

  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        totnodens(i,j,k) = p_theta_levels(i,j,k)                               &
                    /(boltzmann*t_theta_levels(i,j,k))  ! molec/m3
      END DO
    END DO
  END DO

  ! Determine preciptation totals for wet deposition or conservation scheme
  IF ((.NOT. ukca_config%l_ukca_wetdep_off) .OR.                               &
      ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_stratcfc .OR.           &
      ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat) THEN
    ls_ppn3d(:,:,:)   = ls_rain3d  (:,:,:) + ls_snow3d  (:,:,:)
    ls_ppn3d(:,:,:)   = MAX(ls_ppn3d  (:,:,:), 0.0)
    IF (ukca_config%l_param_conv) THEN
      conv_ppn3d(:,:,:) = conv_rain3d(:,:,:) + conv_snow3d(:,:,:)
      conv_ppn3d(:,:,:) = MAX(conv_ppn3d(:,:,:), 0.0)
    ELSE
      conv_ppn3d(:,:,:) = 0.0
    END IF
  ELSE
    ls_ppn3d(:,:,:)   = 0.0
    conv_ppn3d(:,:,:) = 0.0
  END IF

  !-----------------------------------------------------------------------
  ! Calculate relative humidity (expressed as a fraction) if not provided.
  ! For GLOMAP-mode, also calculate clear sky relative humidity if not
  ! provided and saturation vapour pressure if needed for Activate.
  !-----------------------------------------------------------------------

  ! Set flags to indicate which fields are required depending on the
  ! science configuration and whether this is a chemistry time step.
  l_using_rh = do_chemistry
  ! If ukca_volcanic_so2 with plumeria call is TRUE, rel_humid_frac is required
  ! on UKCA timesteps
  IF (ukca_config%l_ukca_so2ems_expvolc .AND.                                  &
      ukca_config%l_ukca_so2ems_plumeria) THEN
    l_using_rh = ukca_config%l_ukca_mode
  END IF
  l_using_rh_clr = ukca_config%l_ukca_mode .AND. (do_chemistry .OR.            &
                   glomap_config%l_ukca_fine_no3_prod .OR.                     &
                   glomap_config%l_ukca_coarse_no3_prod)
  l_using_svp =                                                                &
    ukca_config%l_ukca_mode .AND.                                              &
    (glomap_config%i_ukca_activation_scheme == i_ukca_activation_arg) .AND.    &
    do_chemistry

  ! Allocate fields as necessary if not already available.

  ! rel_humid_frac must always be allocated on chemistry time steps for the
  ! chemistry and aerosol calls, and for ukca_volcanic_so2 subroutine
  ! with plumeria call
  IF (.NOT. ALLOCATED(rel_humid_frac)) THEN
    ALLOCATE(rel_humid_frac(row_length, rows, model_levels))
    rel_humid_frac = 0.0
  END IF
  ! rel_humid_frac_clr must always be allocated for the ukca_emiss_ctl call
  IF (.NOT. ALLOCATED(rel_humid_frac_clr)) THEN
    ALLOCATE(rel_humid_frac_clr(row_length, rows, model_levels))
    rel_humid_frac_clr = 0.0
  END IF
  ! qsvp need only be allocated when used
  IF (.NOT. ALLOCATED(qsvp) .AND. l_using_svp) THEN
    ALLOCATE(qsvp(row_length, rows, model_levels))
    qsvp = 0.0
  END IF

  IF (l_um_atmos_ukca_humidity) THEN

    ! Using external UM humidity routines for RH, clear sky RH and SVP.

    ! Use external UM atmos_ukca_humidity routine to calculate relative
    ! humidity, saturation mixing ratios (for clear sky relative humidity
    ! calculation) and saturation vapour pressure as required.
    ! Possible combinations required are:
    ! RH only, RH & sat MR, RH & sat MR & SVP (on chemistry time steps) or
    ! sat MR only (on non-chemistry time steps)

    IF (l_using_rh_clr) THEN
      ALLOCATE(qsatmr_wat(row_length,rows,model_levels))
      ALLOCATE(qsatmr_ice_wat(row_length,rows,model_levels))
    END IF

    IF (l_using_rh) THEN

      IF (l_using_rh_clr) THEN
        IF (l_using_svp) THEN
          ! RH & sat MR (for clear-sky RH) & SVP
          CALL atmos_ukca_humidity(row_length, rows, model_levels,             &
                                   t_theta_levels, p_theta_levels,             &
                                   q=q, rh=rel_humid_frac,                     &
                                   qsatmr_wat=qsatmr_wat,                      &
                                   qsatmr_ice_wat=qsatmr_ice_wat,              &
                                   svp=qsvp)
        ELSE
          ! RH & sat MR (for clear-sky RH)
          CALL atmos_ukca_humidity(row_length, rows, model_levels,             &
                                   t_theta_levels, p_theta_levels,             &
                                   q=q, rh=rel_humid_frac,                     &
                                   qsatmr_wat=qsatmr_wat,                      &
                                   qsatmr_ice_wat=qsatmr_ice_wat)

        END IF
      ELSE
        ! RH only
        CALL atmos_ukca_humidity(row_length, rows, model_levels,               &
                                 t_theta_levels, p_theta_levels,               &
                                 q=q, rh=rel_humid_frac)
      END IF

      ! Cap relative humidity values: 0.0 <= RH < 1.0
      !!!! Note that values of 1.0 or greater will be capped to 0.999 but
      !!!! values between 0.999 and 1.0 will not be modified. This lack of
      !!!! consistency is undesirable. It is noted but left unchanged at
      !!!! present to preserve bit-comparability of results.
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            IF (rel_humid_frac(i,j,k) < 0.0) THEN
              rel_humid_frac(i,j,k) = 0.0
            ELSE IF (rel_humid_frac(i,j,k) >= 1.0) THEN
              rel_humid_frac(i,j,k) = 0.999
            END IF
          END DO
        END DO
      END DO

    ELSE IF (l_using_rh_clr) THEN

      ! sat MR (for clear-sky RH) only
      CALL atmos_ukca_humidity(row_length, rows, model_levels,                 &
                               t_theta_levels, p_theta_levels,                 &
                               qsatmr_wat=qsatmr_wat,                          &
                               qsatmr_ice_wat=qsatmr_ice_wat)

    END IF

    ! Calculate clear sky RH using external UM lsp_qclear routine

    IF (l_using_rh_clr) THEN

      ALLOCATE(q1d(row_length*rows))
      ALLOCATE(qsatmr_ice_wat1d(row_length*rows))
      ALLOCATE(qsatmr_wat1d(row_length*rows))
      ALLOCATE(qcf1d(row_length*rows))
      ALLOCATE(cloud_liq_frac1d(row_length*rows))
      ALLOCATE(cloud_frac1d(row_length*rows))
      ALLOCATE(rhcrit1d(row_length*rows))
      ALLOCATE(q_clr1d(row_length*rows))
      ALLOCATE(rel_humid_frac_clr1d(row_length*rows))

      DO k = 1, model_levels

        q1d              = RESHAPE(q(:,:,k), [row_length * rows])
        qsatmr_ice_wat1d = RESHAPE(qsatmr_ice_wat(:,:,k), [row_length * rows])
        qsatmr_wat1d     = RESHAPE(qsatmr_wat(:,:,k), [row_length * rows])
        qcf1d            = RESHAPE(qcf(:,:,k), [row_length * rows])
        cloud_liq_frac1d = RESHAPE(cloud_liq_frac(:,:,k), [row_length * rows])
        cloud_frac1d     = RESHAPE(cloud_frac(:,:,k), [row_length * rows])

        rhcrit1d(:)      = rhcrit(k)

        !!!! NOTE - CORRECTION NEEDED: q_clr1d should be a mixing ratio for
        !!!! calculating RH which means that the qsp_clear input variables q1d,
        !!!! and qcf1d should be mixing ratios but they are are currently
        !!!! specific to moist air mass. This should be corrected
        !!!! at the same time as related issues noted in subroutine
        !!!! atmos_ukca_humidity. That routine currently returns a
        !!!! saturation specific humidity for qsatmr_ice_wat1d which should
        !!!! likewise be a mixing ratio.
        CALL lsp_qclear(q1d, qsatmr_ice_wat1d, qsatmr_wat1d,                   &
             qcf1d, cloud_liq_frac1d, cloud_frac1d,                            &
             rhcrit1d, q_clr1d, SIZE(q1d))

        rel_humid_frac_clr1d(:) = q_clr1d(:) / qsatmr_wat1d(:)

        rel_humid_frac_clr(:,:,k) = RESHAPE(rel_humid_frac_clr1d,              &
                                            [row_length, rows])

        ! Cap relative humidity values: 0.0 <= RH < 1.0
        !!!! Note that values of 1.0 or greater will be capped to 0.999 but
        !!!! values between 0.999 and 1.0 will not be modified. This lack of
        !!!! consistency is undesirable. It is noted but left unchanged at
        !!!! present to preserve bit-comparability of results.
        DO j = 1, rows
          DO i = 1, row_length
            IF (rel_humid_frac_clr(i,j,k) < 0.0) THEN
              rel_humid_frac_clr(i,j,k) = 0.0
            ELSE IF (rel_humid_frac_clr(i,j,k) >= 1.0) THEN
              rel_humid_frac_clr(i,j,k) = 0.999
            END IF
          END DO
        END DO

      END DO

      IF (ALLOCATED(rel_humid_frac_clr1d)) DEALLOCATE(rel_humid_frac_clr1d)
      IF (ALLOCATED(q_clr1d)) DEALLOCATE(q_clr1d)
      IF (ALLOCATED(rhcrit1d)) DEALLOCATE(rhcrit1d)
      IF (ALLOCATED(cloud_frac1d)) DEALLOCATE(cloud_frac1d)
      IF (ALLOCATED(cloud_liq_frac1d)) DEALLOCATE(cloud_liq_frac1d)
      IF (ALLOCATED(qcf1d)) DEALLOCATE(qcf1d)
      IF (ALLOCATED(qsatmr_wat1d)) DEALLOCATE(qsatmr_wat1d)
      IF (ALLOCATED(qsatmr_ice_wat1d)) DEALLOCATE(qsatmr_ice_wat1d)
      IF (ALLOCATED(q1d)) DEALLOCATE(q1d)
      IF (ALLOCATED(qsatmr_ice_wat)) DEALLOCATE(qsatmr_ice_wat)
      IF (ALLOCATED(qsatmr_wat)) DEALLOCATE(qsatmr_wat)

    END IF ! l_using_rh_clr

  ELSE

    ! Using UKCA humidity calculations for RH, clear sky RH and SVP.

    ! Set flags to indicate which fields need calculating. Currently, all of
    ! the used fields are calculated here pending the addition of functionality
    ! to supply any/all of the fields as environmental drivers, so the
    ! 'calc' flags just duplicate the 'using' flags. Subsequently,
    ! the used fields will only need to be calculated if not already supplied
    ! and the flags should be modified accordingly.

    l_calc_rh = l_using_rh
    l_calc_rh_clr = l_using_rh_clr
    l_calc_svp = l_using_svp

    ! Calculate saturation mixing ratio and/or saturation vapour pressure
    ! (with respect to liquid) as required
    IF (l_calc_rh .OR. l_calc_rh_clr) THEN
      ALLOCATE(water_vapour_mr_sat(row_length, rows, model_levels))
      IF (l_calc_svp) THEN
        ! sat MR & SVP
        CALL ukca_vmrsat_liq(row_length, rows, model_levels, t_theta_levels,   &
                             pres=p_theta_levels, vmr_sat=water_vapour_mr_sat, &
                             svp=qsvp)
      ELSE
        ! sat MR only
        CALL ukca_vmrsat_liq(row_length, rows, model_levels, t_theta_levels,   &
                             pres=p_theta_levels, vmr_sat=water_vapour_mr_sat)
      END IF
    ELSE IF (l_calc_svp) THEN
      ! SVP only
      CALL ukca_vmrsat_liq(row_length, rows, model_levels, t_theta_levels,     &
                           svp=qsvp)
    END IF

    ! Calculate water vapour mixing ratio for RH calculations.
    ! Determine mixing ratio from specific humidity taking into account
    ! mass fraction of cloud liquid and frozen water for the air parcel.
    IF (l_calc_rh .OR. l_calc_rh_clr) THEN
      ALLOCATE(water_vapour_mr(row_length, rows, model_levels))
      water_vapour_mr(:,:,:) = q(:,:,:) /                                      &
                               (1.0 - q(:,:,:) - qcl(:,:,:) - qcf(:,:,:))
    END IF

    ! Calculate RH if required
    IF (l_calc_rh) THEN
      rel_humid_frac(:,:,:) = MAX(MIN(                                         &
        water_vapour_mr(:,:,:) / water_vapour_mr_sat(:,:,:),                   &
        rh_max), rh_min)
    END IF

    ! Calculate clear sky RH if required
    IF (l_calc_rh_clr) THEN
      ALLOCATE(water_vapour_mr_clr(row_length, rows, model_levels))
      water_vapour_mr_clr = ukca_vmr_clear_sky(row_length, rows, model_levels, &
                            water_vapour_mr, water_vapour_mr_sat,              &
                            cloud_liq_frac)
      rel_humid_frac_clr(:,:,:) = MAX(MIN(                                     &
        water_vapour_mr_clr(:,:,:) / water_vapour_mr_sat(:,:,:),               &
        rh_max), rh_min)
      DEALLOCATE(water_vapour_mr_clr)
    END IF

    IF (ALLOCATED(water_vapour_mr)) DEALLOCATE(water_vapour_mr)
    IF (ALLOCATED(water_vapour_mr_sat)) DEALLOCATE(water_vapour_mr_sat)

  END IF  ! l_um_atmos_ukca_humidity

  !-----------------------------------------------------------------------

  ! If land fraction is not already allocated set to 1.
  ! This only happens if coastal tiling is off.
  ! When coastal tiling is off all land points have a land fraction of 1
  IF (.NOT. ALLOCATED(fland) ) THEN
    ALLOCATE(fland(land_points))
    fland(:)=1.0
  END IF

  ! Set up land fraction (required for MODE and DRY DEPN).
  ! Input data fland is a 1D packed array on land points only.
  ! We require a 2D data valid at all points.
  land_fraction(:,:)=0.0
  DO l = 1, land_points
    j = (land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    land_fraction(i,j) = fland(l)
  END DO

  ! ----------------------------------------------------------------------
  ! Print out values for debugging
  ! ----------------------------------------------------------------------

  IF (PrintStatus >= PrStatus_Oper) THEN
    CALL ukca_pr_inputs(error_code_ptr, timestep_number, environ_ptrs,         &
                        land_fraction, thick_bl_levels, t_theta_levels,        &
                        rel_humid_frac, z_half,                                &
                        error_message=error_message,                           &
                        error_routine=error_routine)
    IF (error_code_ptr > 0) THEN
      IF (lhook)                                                               &
        CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  END IF

  !       Warning statements about possible mismatch between interactive
  !       methane emissions and emissions ancillary
  IF (l_first_call .AND. .NOT. ukca_config%l_ukca_emissions_off) THEN
    emiss_input = 'NetCDF file'

    IF (ukca_config%l_ukca_qch4inter .AND.                                     &
        (.NOT. ukca_config%l_ukca_prescribech4)) THEN
      cmessage = 'CH4 WETLANDS EMS ARE ON - '  // TRIM(emiss_input) //         &
                 ' SHOULD NOT contain wetland ems'
      errcode=-1

      CALL ereport(RoutineName,errcode,cmessage)
    ELSE IF (.NOT. ukca_config%l_ukca_prescribech4 .AND.                       &
             ANY(advt == 'CH4       ')) THEN
      cmessage = 'CH4 WETLANDS EMS ARE OFF - ' // TRIM(emiss_input) //         &
                 ' SHOULD contain wetland ems'
      errcode=-2
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  END IF  ! End of IF (l_first_call)

END IF    ! End of IF (l_ukca_chem)

! Set up tile info
IF (ukca_config%l_ukca_chem) THEN
  tile_pts(:) = 0
  tile_index(:,:) = 0
  IF (ukca_config%l_ukca_intdd) THEN
    DO n = 1, ukca_config%ntype
      tile_index(:,n) = 0
      DO l = 1, land_points
        IF (l_tile_active(l,n)) THEN
          tile_pts(n) = tile_pts(n) + 1
          tile_index(tile_pts(n),n) = l
        END IF
      END DO
    END DO
  END IF
END IF

! ----------------------------------------------------------------------
! 4. Calls to science subroutines
! ----------------------------------------------------------------------
! 4.1 Age-of-air
! ----------------------------------------------------------------------

IF ( ukca_config%l_ukca_ageair ) THEN
  CALL ukca_age_air(row_length, rows, model_levels, ukca_config%timestep,      &
                    z_top_of_model, all_tracers_names, all_tracers)
END IF    ! ukca_config%l_ukca_ageair

! ----------------------------------------------------------------------
! 4.2 Chemistry related pre-processing
! ----------------------------------------------------------------------

IF (ukca_config%l_ukca_chem) THEN

#if !defined(LFRIC)
  IF (ukca_config%l_timer) CALL timer('UKCA CHEMISTRY MODEL',5)
#endif

  ! Transform tracers to ensure elemental conservation
  IF (ukca_config%l_tracer_lumping) THEN
    CALL ukca_transform_halogen(rows,row_length,model_levels,n_tracers,        &
                                unlump_species,all_tracers)
  END IF

  ! Locate the tropopause
  IF (ukca_config%l_fix_tropopause_level) THEN
    ! Fix tropopause at a specified model level
    k = ukca_config%fixed_tropopause_level
    tropopause_level(:,:) = k
    theta_trop(:,:) = rmdi
    pv_trop(:,:) = rmdi
    ! Catch in case using box model with 1 model level
    IF (model_levels > 1) THEN
      p_tropopause(:,:) = p_theta_levels(:,:,k)
      l_stratosphere(:,:,1:k) = .FALSE.
      l_stratosphere(:,:,k+1:model_levels) = .TRUE.
    ELSE
      ! Default to running in the troposphere if with one model level
      p_tropopause(:,:) = rmdi
      l_stratosphere(:,:,1) = .FALSE.
    END IF
  ELSE
    ! Calculate tropopause pressure using a combined
    ! theta and PV surface
    CALL ukca_calc_tropopause(row_length, rows, model_levels, r_theta_levels,  &
                              latitude, theta, pv_on_theta_mlevs,              &
                              p_layer_boundaries, p_theta_levels)
  END IF

  IF (L_asad_use_chem_diags .AND. L_asad_use_trop_mask)                        &
       CALL asad_tropospheric_mask(ierr)

  !       Calculate the difference in tracers per timestep in the
  !       troposphere (moles/s) due to transport using the
  !       trmol_post_atmstep array from this timestep and the
  !       trmol_post_chem array from the previous timestep.
  !       Do this only for those stratospheric flux diagnostics
  !       which are switched on.

  IF (l_first_call .AND. n_strat_fluxdiags > 0) THEN

    IF (.NOT. ALLOCATED(trmol_post_chem) )                                     &
      ALLOCATE(trmol_post_chem(row_length,rows,model_levels,                   &
                                 n_chem_tracers))  ! moles
    trmol_post_chem = 0.0

    DO l=1,n_strat_fluxdiags
      DO k=1,model_levels
        DO j=1,rows
          DO i=1,row_length
            strat_fluxdiags(i,j,k,l) = 0.0
          END DO
        END DO
      END DO
    END DO

  ELSE IF ((.NOT. l_first_call) .AND. n_strat_fluxdiags > 0) THEN

    IF (.NOT. ALLOCATED(trmol_post_atmstep) )                                  &
      ALLOCATE(trmol_post_atmstep(row_length,rows,model_levels,                &
                            n_chem_tracers))  ! moles
    trmol_post_atmstep = 0.0

    DO l=1,n_chem_tracers
      DO k=1,model_levels
        DO j=1,rows
          DO i=1,row_length
            trmol_post_atmstep(i,j,k,l) = all_tracers(i,j,k,l)                 &
                *totnodens(i,j,k)*grid_volume(i,j,k)/(c_species(l)*avogadro)
                                                                        ! moles
          END DO
        END DO
      END DO
    END DO

    icnt = 0
    DO l=1,n_chem_tracers
      IF (UkcaD1codes(istrat_first+l-1)%item /= imdi) THEN
        icnt = icnt + 1

        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              IF (L_stratosphere(i,j,k)) THEN                    ! stratosphere
                strat_fluxdiags(i,j,k,icnt) = 0.0
              ELSE                                               ! troposphere
                strat_fluxdiags(i,j,k,icnt) =                                  &
                 (trmol_post_atmstep(i,j,k,l)-                                 &
                  trmol_post_chem(i,j,k,l))/ukca_config%timestep ! moles/sec
              END IF
            END DO
          END DO
        END DO
      END IF
    END DO     ! n_chem_tracers

    IF (ALLOCATED(trmol_post_atmstep)) DEALLOCATE(trmol_post_atmstep)  ! moles

  END IF      ! End of IF l_first_call and n_strat_fluxdiags statement

  IF (.NOT. l_first_call) THEN ! DO NOT CALL ON FIRST TIMESTEP
     ! ASAD Flux Diags STE
    l_store_value = .FALSE.
    IF (L_asad_use_chem_diags .AND. L_asad_use_STE)                            &
         CALL asad_tendency_ste(                                               &
           row_length, rows, model_levels, n_chem_tracers,                     &
           all_tracers(:,:,:,1:n_chem_tracers),                                &
           totnodens, grid_volume, ukca_config%timestep,                       &
           calculate_STE, l_store_value, ierr)
  END IF

  l_store_value = .TRUE.
  IF (do_chemistry .AND. (L_asad_use_chem_diags                                &
       .AND. L_asad_use_tendency))                                             &
       CALL asad_tendency_ste(                                                 &
         row_length, rows, model_levels, n_chem_tracers,                       &
         all_tracers(:,:,:,1:n_chem_tracers),                                  &
         totnodens, grid_volume, ukca_config%timestep,                         &
       calculate_tendency,l_store_value,ierr)

  ! Setup ozone field for use in stratosphere with troposphere-only
  ! chemistry schemes
  IF (do_chemistry .AND. .NOT. ukca_config%l_ukca_offline_be) THEN
    IF (.NOT. ALLOCATED(env_ozone3d))                                          &
      ALLOCATE(env_ozone3d(row_length,rows,model_levels))
    IF ((ukca_config%l_ukca_trop .OR. ukca_config%l_ukca_raq .OR.              &
         ukca_config%l_ukca_raqaero .OR. ukca_config%l_ukca_tropisop) .AND.    &
         ukca_config%nlev_above_trop_o3_env < ukca_config%model_levels) THEN
      IF (SIZE(o3_offline, DIM=1) == 1) THEN ! Use zonal avg. (or single column)
        DO i = 1,row_length
          env_ozone3d(i,:,:) = o3_offline(1,:,:)
        END DO
      ELSE                                   ! 3-D ozone specifed
        env_ozone3d(:,:,:) = o3_offline(:,:,:)
      END IF
    ELSE
      env_ozone3d(:,:,:) = 0.0               ! Not used, set to zero
    END IF
  END IF

  ! ----------------------------------------------------------------------
  ! 4.3 Emissions
  ! ----------------------------------------------------------------------

  IF (.NOT. ukca_config%l_ukca_emissions_off) THEN

    ! Emission part. Before calling emissions,
    ! set up array of lower-boundary mixing ratios for stratospheric species.
    ! If for some species actual emissions (not prescribed lower boundary
    ! conditions) are required, remove these species from lbc_spec, e.g., for
    ! CH4.
    ! Mass mixing ratios are best estimates, taken from the WMO Scientific
    ! assessment of ozone depletion (2002), for 2000 conditions.
    ! H2 corresponds to 0.5 ppmv.
    ! Select lower boundary MMRs. Either from gui or predefined constants
    ! Note that for bromine compounds (MeBr, H1211, H1301) are increased
    ! by 24% to account for non-included species and HBr is set to 0.

    IF (ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_stratcfc .OR.         &
        ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat) THEN

      ! LBC code contained within UKCA_SCENARIO_CTL
      CALL ukca_scenario_ctl(n_boundary_vals, lbc_spec,                        &
                             i_year, i_day_number, lbc_mmr,                    &
                             .NOT. ukca_config%l_ukca_stratcfc)

    END IF

  END IF  ! .NOT. ukca_config%l_ukca_emissions_off

  !       Calculate solar zenith angle.
  !       Cos zenith angle and its integral over the day are required if
  !       applying a diurnal cycle to offline oxidants and/or to isoprene
  !       emissions.

  IF (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be .OR.       &
      (ukca_config%l_diurnal_isopems .AND.                                     &
       .NOT. ukca_config%l_ukca_emissions_off)) THEN

    IF (.NOT. ALLOCATED(cos_zenith_angle))                                     &
      ALLOCATE(cos_zenith_angle(row_length, rows))
    IF (.NOT. ALLOCATED(int_zenith_angle))                                     &
      ALLOCATE(int_zenith_angle(row_length, rows))

    secondssincemidnight = REAL(i_hour_previous * 3600                         &
                         +      i_minute_previous * 60                         &
                         +      i_second_previous)

    ! Compute COS(sza) integral only when day changes
    ! i.e. on first timestep of 'next' day, determined by checking
    ! if previous timestep was fully divisible by timesteps_per_day.
    IF ( l_first_call .OR.                                                     &
         MOD(timestep_number-1, ukca_config%timesteps_per_day) == 0 ) THEN
      int_zenith_angle(:,:) = 0.0
      secondssincemidnight_sav = secondssincemidnight
      sindec_sav = sin_declination
      eq_time_sav = equation_of_time
      DO kk = 1, (ukca_config%timesteps_per_hour*24)
        ssmn_incr = secondssincemidnight + ukca_config%timestep*(kk-1)
        CALL ukca_solang(sin_declination, ssmn_incr,                           &
                         ukca_config%timestep, equation_of_time,               &
                         sin_latitude,                                         &
                         cos_latitude,                                         &
                         longitude,                                            &
                         theta_field_size, cos_zenith_angle)

        WHERE (cos_zenith_angle > 0.0)
          int_zenith_angle(:,:) = int_zenith_angle(:,:) +                      &
            cos_zenith_angle(:,:) * ukca_config%timestep
        END WHERE
      END DO
    END IF

    ! Compute COS(sza) integral every timestep if persistence is off
    ! using saved values from when day changes or first timestep
    IF ( ukca_config%l_ukca_persist_off .AND. .NOT. (l_first_call .OR.         &
         MOD(timestep_number-1, ukca_config%timesteps_per_day) == 0) ) THEN
      int_zenith_angle(:,:) = 0.0
      DO kk = 1, (ukca_config%timesteps_per_hour*24)
        ssmn_incr = secondssincemidnight_sav + ukca_config%timestep*(kk-1)
        CALL ukca_solang(sindec_sav, ssmn_incr,                                &
                         ukca_config%timestep, eq_time_sav,                    &
                         sin_latitude,                                         &
                         cos_latitude,                                         &
                         longitude,                                            &
                         theta_field_size, cos_zenith_angle)

        WHERE (cos_zenith_angle > 0.0)
          int_zenith_angle(:,:) = int_zenith_angle(:,:) +                      &
            cos_zenith_angle(:,:) * ukca_config%timestep
        END WHERE
      END DO
    END IF

    ! This call needs to be after above loops to prevent errors with value in
    !  cos_zenith_angle
    CALL ukca_solang(sin_declination, secondssincemidnight,                    &
                     ukca_config%timestep, equation_of_time,                   &
                     sin_latitude,                                             &
                     cos_latitude,                                             &
                     longitude,                                                &
                     theta_field_size, cos_zenith_angle )

  END IF  ! ukca_config%l_ukca_offline etc.

  IF (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) THEN
    ! Compute COS(sza) integral and daylength, also stores cos(zenith) for
    !  Offline chemistry diurnal cycle
    CALL ukca_int_cosz(row_length, rows,                                       &
                       sin_latitude, cos_latitude, tan_latitude,               &
                       sin_declination, cos_zenith_angle, int_zenith_angle)
  END IF

  ! Emissions system (based on NetCDF emission input in the case of
  ! non-interactive emissions)

  ! Option to turn off emissions e.g. in UKCA box model
  IF (ukca_config%l_ukca_emissions_off) THEN

    IF (printstatus >= prstatus_oper) THEN
      WRITE(umMessage,'(A)')                                                   &
          'UKCA_MAIN1, NOT Calling ukca_emiss_ctl'
      CALL umPrint(umMessage,src=RoutineName)
    END IF

  ELSE

    IF (.NOT. ALLOCATED(cos_zenith_angle))                                     &
      ALLOCATE(cos_zenith_angle(row_length, rows))
    IF (.NOT. ALLOCATED(int_zenith_angle))                                     &
      ALLOCATE(int_zenith_angle(row_length, rows))

    ! Read emission fields and calculate online emissions, fill in
    ! the emissions structure, inject emissions, do tracer mixing and
    ! finally output emission diagnostics.
    CALL ukca_emiss_ctl (                                                      &
        row_length, rows, ukca_config%bl_levels, model_levels,                 &
        n_chem_tracers + n_aero_tracers + n_mode_tracers,                      &
        glomap_config%n_dust_emissions,                                        &
        latitude, longitude, sin_latitude, cos_latitude, tan_latitude,         &
        i_year, i_month, i_day, i_hour, ukca_config%timestep,                  &
        l_first_call, land_points, land_index,                                 &
        conv_cloud_base, conv_cloud_top, delta_lambda, delta_phi,              &
        r_theta_levels(1:row_length, 1:rows, 0:model_levels), grid_surf_area,  &
        cos_zenith_angle, int_zenith_angle, land_fraction, tropopause_level,   &
        r_rho_levels(1:row_length, 1:rows, 1:model_levels), t_theta_levels,    &
        p_theta_levels, p_layer_boundaries, rel_humid_frac_clr, grid_airmass,  &
        land_sea_mask, rel_humid_frac, plumeria_height,                        &
        theta, q, qcl, qcf, exner_rho_levels, rho_r2,                          &
        kent, kent_dsc, rhokh_rdz, dtrdz,                                      &
        we_lim, t_frac, zrzi, we_lim_dsc, t_frac_dsc, zrzi_dsc,                &
        zbl, zhsc, z_half, ch4_wetl_emiss, seaice_frac, grid_area_fullht,      &
        dust_flux, u_scalar_10m, tstar, dms_sea_conc, chloro_sea,              &
        dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6,      &
        all_tracers(:,:,:,1:n_chem_tracers+n_aero_tracers+n_mode_tracers),     &
        ext_cg_flash, ext_ic_flash,                                            &
        SIZE(stashwork38), stashwork38,                                        &
        SIZE(stashwork50), stashwork50)

    IF (error_code_ptr > 0) THEN
      IF (lhook) THEN
        CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      END IF
      RETURN
    END IF

  END IF ! l_ukca_emissions_off

  IF (ALLOCATED(cos_zenith_angle)) DEALLOCATE(cos_zenith_angle)

  IF (L_asad_use_chem_diags .AND. L_asad_use_mass_diagnostic)                  &
       CALL asad_mass_diagnostic(                                              &
       row_length,                                                             &
       rows,                                                                   &
       model_levels,                                                           &
       grid_airmass,                                                           &
       ierr)

    ! ----------------------------------------------------------------------
    ! 4.4 Call chemistry routines
    ! ----------------------------------------------------------------------
      ! Do chemistry calculation here (at chemistry timesteps)
  IF (do_chemistry) THEN

    IF (.NOT. ALLOCATED(t_chem)) ALLOCATE(t_chem(row_length,rows,model_levels))
    IF (.NOT. ALLOCATED(q_chem)) ALLOCATE(q_chem(row_length,rows,model_levels))
    t_chem = t_theta_levels
    q_chem = q

    !======================================
    ! PASSIVE OZONE - CheST/StratChem ONLY
    !  n_passive is set in ukca_calc_cspecies
             ! Do passive ozone
    IF (n_passive > 0) THEN
       ! assumes 360 day calendar
      IF ((i_day_number == 91) .OR. (i_day_number == 301)) THEN
        IF (i_hour == 0) THEN
          all_tracers(:,:,:,n_passive) = all_tracers(:,:,:,n_o3)
        END IF
      END IF
    END IF
    !======================================
    l_store_value = .TRUE.
    IF (L_asad_use_chem_diags .AND. L_asad_use_tendency)                       &
         CALL asad_tendency_ste(                                               &
           row_length, rows, model_levels, n_chem_tracers,                     &
           all_tracers(:,:,:,1:n_chem_tracers),                                &
           totnodens, grid_volume, ukca_config%timestep,                       &
           calculate_tendency, l_store_value, ierr)

    IF (ukca_config%l_ukca_mode) THEN
      IF (ukca_config%l_ukca_aerchem) THEN
        ! SO2 wet oxidation and H2SO4 updating done in UKCA
        wetox_in_aer = 0
        uph2so4inaer = 0
      ELSE IF (ukca_config%l_ukca_nr_aqchem .OR.                               &
               ukca_config%l_ukca_offline_be .OR.                              &
               ukca_config%l_ukca_raqaero) THEN
        ! SO2 wet oxidation in UKCA, and H2SO4 updating done in MODE
        wetox_in_aer = 0
        uph2so4inaer = 1
      ELSE                        ! No SO2 oxidation in UKCA chemistry
        wetox_in_aer = 1          ! Sulphur emissions are still needed
        uph2so4inaer = 1
      END IF
      IF (printstatus >= prstatus_oper) THEN
        WRITE(umMessage,'(A,I5)') 'uph2so4inaer is set to: ',                  &
              uph2so4inaer
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(umMessage,'(A,I5)') 'wetox_in_aer is set to: ',                  &
              wetox_in_aer
        CALL umPrint(umMessage,src=RoutineName)
      END IF
    ELSE
      wetox_in_aer = 0
      uph2so4inaer = 0
    END IF                 ! ukca_config%l_ukca_mode

    IF (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) THEN
      ! Allocate diagnostic arrays for offline chemistry
      IF (.NOT. ALLOCATED(o3_offline_diag))                                    &
        ALLOCATE(o3_offline_diag(theta_field_size,1:model_levels))
      IF (.NOT. ALLOCATED(oh_offline_diag))                                    &
        ALLOCATE(oh_offline_diag(theta_field_size,1:model_levels))
      IF (.NOT. ALLOCATED(no3_offline_diag))                                   &
        ALLOCATE(no3_offline_diag(theta_field_size,1:model_levels))
      IF (.NOT. ALLOCATED(ho2_offline_diag))                                   &
        ALLOCATE(ho2_offline_diag(theta_field_size,1:model_levels))
      ! Initialise to zero
      o3_offline_diag(:,:) = 0.0
      oh_offline_diag(:,:) = 0.0
      no3_offline_diag(:,:) = 0.0
      ho2_offline_diag(:,:) = 0.0
    END IF

    ! Persistence of spatial fields between timesteps will be turned off
    ! for LFRic
    ! Now setting up have_nat3d every timestep so saving will be unnecessary
    ! Only necessary if not l_ukca_offline_be
    IF (.NOT. ukca_config%l_ukca_offline_be) THEN
      IF (.NOT. ALLOCATED(have_nat3d))                                         &
         ALLOCATE(have_nat3d(row_length,rows,model_levels))

      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            ! True (allow NAT everywhere) if l_ukca_limit_nat=False,
            !      which is the default value for logical.
            ! False (supress everywhere) initially if logical is True
            have_nat3d(i,j,k) = (.NOT. ukca_config%l_ukca_limit_nat)
          END DO
        END DO
      END DO

      IF (ukca_config%l_ukca_limit_nat) THEN
        h_atmos = r_theta_levels(1,1,model_levels) - planet_radius
        DO k = 1, model_levels
          level_above_ht = ( (h_atmos * eta_theta_levels(k)) > nat_limit_ht )
          DO j = 1, rows
            DO i = 1, row_length
              have_nat3d(i,j,k) = level_above_ht
            END DO
          END DO
        END DO
      END IF
    END IF ! .NOT. ukca_config%l_ukca_offline_be

    ! Initially assign 3D pH array the default global value of 5.0
    IF (.NOT. ALLOCATED(H_plus_3d_arr))                                        &
        ALLOCATE(H_plus_3d_arr(row_length,rows,model_levels))
    H_plus_3d_arr(:,:,:) = H_plus

    ! Check whether to interactively calculate cloud pH values
    IF (ukca_config%l_ukca_intph) THEN ! Use logical set in Rose suite
      ! Calculate online pH values now to be used in chemistry control
      WRITE(umMessage,'(A)') 'Run Interactive Cloud pH routine'
      CALL umPrint(umMessage,src=RoutineName)
      CALL ukca_calc_cloud_ph(                                                 &
              all_tracers(:,:,:,1:n_chem_tracers+n_aero_tracers),              &
              row_length, rows, model_levels,                                  &
              p_theta_levels, t_chem, H_plus_3d_arr)

      WRITE(umMessage,'(A)') 'Finished Interactive Cloud pH routine'
      CALL umPrint(umMessage,src=RoutineName)
    END IF
    ! now pass H_plus array to routines that use it in chem ctl

    ! Prepare to run ASAD
    CALL ukca_chemistry_setup(                                                 &
         row_length, rows, model_levels,                                       &
         theta_field_size,                                                     &
         ukca_config%bl_levels,                                                &
         n_chem_tracers+n_aero_tracers,                                        &
         ukca_config%ntype,                                                    &
         ukca_config%npft,                                                     &
         SIZE(stashwork50),                                                    &
         i_month, i_day_number, i_hour,                                        &
         land_points, land_index, tile_pts, tile_index,                        &
         istore_h2so4,                                                         &
         nlev_with_ddep,                                                       &
         r_minute - ukca_config%timestep/60.0,                                 &
         REAL(ukca_config%chem_timestep),                                      &
         latitude, longitude,                                                  &
         sin_latitude, tan_latitude,                                           &
         t_chem,                                                               &
         p_theta_levels,                                                       &
         rel_humid_frac,                                                       &
         p_layer_boundaries,                                                   &
         Tstar,                                                                &
         Thick_bl_levels,                                                      &
         u_s,                                                                  &
         rough_length,                                                         &
         ls_ppn3d, conv_ppn3d,                                                 &
         frac_types,                                                           &
         seaice_frac, stcon, surf_hf,                                          &
         soil_moisture_layer1, fland,                                          &
         laift_lp, canhtft_lp,                                                 &
         z0tile_lp, tstar_tile,                                                &
         zbl,                                                                  &
         H_plus_3d_arr,                                                        &
         zdryrt, zwetrt,                                                       &
         all_tracers(:,:,:,1:n_chem_tracers+n_aero_tracers),                   &
         stashwork50,                                                          &
         l_firstchem                                                           &
         )

    IF (ukca_config%l_ukca_offline_be) THEN

      ! Offline chemistry with explicit backward-Euler solver
      CALL ukca_chemistry_ctl_be(                                              &
           row_length, rows, model_levels,                                     &
           theta_field_size, tot_n_pnts,                                       &
           n_chem_tracers+n_aero_tracers,                                      &
           ukca_config%chem_timestep,                                          &
           k_be_top,                                                           &
           p_theta_levels,                                                     &
           t_chem,                                                             &
           q_chem,                                                             &
           qcl,                                                                &
           all_tracers(:,:,:,1:n_chem_tracers+n_aero_tracers),                 &
           cloud_frac,                                                         &
           grid_volume,                                                        &
           ! Extra variables for new dry dep scheme
           uph2so4inaer,                                                       &
           delso2_wet_h2o2,                                                    &
           delso2_wet_o3,                                                      &
           delh2so4_chem,                                                      &
           delso2_drydep,                                                      &
           delso2_wetdep,                                                      &
           H_plus_3d_arr,                                                      &
           h2o2_offline,                                                       &
           zdryrt, zwetrt, nlev_with_ddep,                                     &
           l_firstchem                                                         &
           )

    ELSE IF (ukca_config%l_ukca_trop .OR. ukca_config%l_ukca_aerchem .OR.      &
             ukca_config%l_ukca_raq .OR. ukca_config%l_ukca_raqaero) THEN

      ! Aerosol mmr / numbers (from CLASSIC aerosol scheme) are only
      ! used by the RAQ chemistry scheme if heterogeneous chemistry
      ! is ON. Even if that functionality is OFF they need to be
      ! allocated before the call to UKCA_CHEMISTRY_CTL.

      IF (.NOT. ALLOCATED(so4_aitken)) THEN
        ALLOCATE(so4_aitken(row_length,rows,model_levels))
        so4_aitken(:,:,:) = min_SO4_val
      END IF

      IF (.NOT. ALLOCATED(so4_accum)) THEN
        ALLOCATE(so4_accum(row_length,rows,model_levels))
        so4_accum(:,:,:) = min_SO4_val
      END IF

      IF (.NOT. ALLOCATED(soot_fresh))                                         &
          ALLOCATE (soot_fresh (row_length, rows, model_levels))

      IF (.NOT. ALLOCATED(soot_aged))                                          &
          ALLOCATE (soot_aged  (row_length, rows, model_levels))

      IF (.NOT. ALLOCATED(ocff_fresh))                                         &
          ALLOCATE (ocff_fresh (row_length, rows, model_levels))

      IF (.NOT. ALLOCATED(ocff_aged))                                          &
          ALLOCATE (ocff_aged  (row_length, rows, model_levels))

      IF (.NOT. ALLOCATED(biogenic))                                           &
          ALLOCATE (biogenic   (row_length, rows, model_levels))

      IF (.NOT. ALLOCATED(sea_salt_film))                                      &
          ALLOCATE (sea_salt_film (row_length, rows, model_levels))

      IF (.NOT. ALLOCATED(sea_salt_jet))                                       &
          ALLOCATE (sea_salt_jet  (row_length, rows, model_levels))

      CALL ukca_chemistry_ctl_tropraq(                                         &
           row_length, rows, model_levels,                                     &
           theta_field_size, tot_n_pnts,                                       &
           n_chem_tracers+n_aero_tracers,                                      &
           REAL(ukca_config%chem_timestep),                                    &
           p_theta_levels,                                                     &
           t_chem,                                                             &
           q_chem,                                                             &
           qcf,                                                                &
           qcl,                                                                &
           rel_humid_frac,                                                     &
           all_tracers(:,:,:,1:n_chem_tracers+n_aero_tracers),                 &
           all_ntp,                                                            &
           cloud_frac,                                                         &
           photol_rates,                                                       &
           grid_volume,                                                        &
           so4_aitken,                                                         &
           so4_accum,                                                          &
           soot_fresh,                                                         &
           soot_aged,                                                          &
           ocff_fresh,                                                         &
           ocff_aged,                                                          &
           biogenic,                                                           &
           sea_salt_film,                                                      &
           sea_salt_jet,                                                       &
           uph2so4inaer,                                                       &
           delso2_wet_h2o2,                                                    &
           delso2_wet_o3,                                                      &
           delh2so4_chem,                                                      &
           delso2_drydep,                                                      &
           delso2_wetdep,                                                      &
           ! Diagnostics
           trop_ch4_mol,                                                       &
           trop_o3_mol,                                                        &
           trop_oh_mol,                                                        &
           strat_ch4_mol,                                                      &
           strat_ch4loss,                                                      &
           SIZE(stashwork50),                                                  &
           stashwork50,                                                        &
           H_plus_3d_arr,                                                      &
           zdryrt, zwetrt, nlev_with_ddep, L_stratosphere,                     &
           l_firstchem                                                         &
           )

    ELSE IF (ukca_config%l_ukca_asad_full) THEN

      CALL ukca_chemistry_ctl_full(                                            &
           row_length, rows, model_levels,                                     &
           theta_field_size, tot_n_pnts,                                       &
           n_chem_tracers+n_aero_tracers,                                      &
           istore_h2so4,                                                       &
           p_theta_levels,                                                     &
           t_chem,                                                             &
           q_chem,                                                             &
           qcf,                                                                &
           qcl,                                                                &
           all_tracers(:,:,:,1:n_chem_tracers+n_aero_tracers),                 &
           all_ntp,                                                            &
           cloud_frac,                                                         &
           photol_rates,                                                       &
           shno3_3d,                                                           &
           grid_volume,                                                        &
           have_nat3d,                                                         &
           uph2so4inaer,                                                       &
           delso2_wet_h2o2,                                                    &
           delso2_wet_o3,                                                      &
           delh2so4_chem,                                                      &
           so4_sa,                                                             &
           ! Diagnostics
           atm_ch4_mol,                                                        &
           atm_co_mol,                                                         &
           atm_n2o_mol,                                                        &
           atm_cf2cl2_mol,                                                     &
           atm_cfcl3_mol,                                                      &
           atm_mebr_mol,                                                       &
           atm_h2_mol,                                                         &
           H_plus_3d_arr,                                                      &
           zdryrt, zwetrt, nlev_with_ddep, L_stratosphere, co2_interactive,    &
           l_firstchem                                                         &
           )

    ELSE IF (ukca_config%l_ukca_asad_columns) THEN

      CALL ukca_chemistry_ctl_col(                                             &
           row_length, rows, model_levels,                                     &
           theta_field_size,                                                   &
           n_chem_tracers+n_aero_tracers,                                      &
           istore_h2so4,                                                       &
           p_theta_levels,                                                     &
           t_chem,                                                             &
           q_chem,                                                             &
           qcf,                                                                &
           qcl,                                                                &
           all_tracers(:,:,:,1:n_chem_tracers+n_aero_tracers),                 &
           all_ntp,                                                            &
           cloud_frac,                                                         &
           photol_rates,                                                       &
           shno3_3d,                                                           &
           grid_volume,                                                        &
           ! Extra variables for new dry dep scheme
           have_nat3d,                                                         &
           uph2so4inaer,                                                       &
           delso2_wet_h2o2,                                                    &
           delso2_wet_o3,                                                      &
           delh2so4_chem,                                                      &
           so4_sa,                                                             &
           ! Diagnostics
           atm_ch4_mol,                                                        &
           atm_co_mol,                                                         &
           atm_n2o_mol,                                                        &
           atm_cf2cl2_mol,                                                     &
           atm_cfcl3_mol,                                                      &
           atm_mebr_mol,                                                       &
           atm_h2_mol,                                                         &
           H_plus_3d_arr,                                                      &
           zdryrt, zwetrt, nlev_with_ddep                                      &
           )

    ELSE

      CALL ukca_chemistry_ctl(                                                 &
           row_length, rows, model_levels,                                     &
           theta_field_size, tot_n_pnts,                                       &
           n_chem_tracers+n_aero_tracers,                                      &
           istore_h2so4,                                                       &
           p_theta_levels,                                                     &
           t_chem,                                                             &
           q_chem,                                                             &
           qcf,                                                                &
           qcl,                                                                &
           all_tracers(:,:,:,1:n_chem_tracers+n_aero_tracers),                 &
           all_ntp,                                                            &
           cloud_frac,                                                         &
           photol_rates,                                                       &
           shno3_3d,                                                           &
           grid_volume,                                                        &
           have_nat3d,                                                         &
           uph2so4inaer,                                                       &
           delso2_wet_h2o2,                                                    &
           delso2_wet_o3,                                                      &
           delh2so4_chem,                                                      &
           so4_sa,                                                             &
           ! Diagnostics
           atm_ch4_mol,                                                        &
           atm_co_mol,                                                         &
           atm_n2o_mol,                                                        &
           atm_cf2cl2_mol,                                                     &
           atm_cfcl3_mol,                                                      &
           atm_mebr_mol,                                                       &
           atm_h2_mol,                                                         &
           H_plus_3d_arr,                                                      &
           zdryrt, zwetrt, nlev_with_ddep, co2_interactive, L_stratosphere,    &
           l_firstchem                                                         &
           )
    END IF

    ! ASAD post-processing
    IF (.NOT. ukca_config%l_ukca_offline_be) THEN
      CALL ukca_chemistry_cleanup(row_length, rows, model_levels,              &
             n_chem_tracers+n_aero_tracers,                                    &
             REAL(ukca_config%chem_timestep),                                  &
             p_theta_levels, ls_ppn3d, conv_ppn3d,                             &
             latitude, env_ozone3d, qcf,                                       &
             r_theta_levels(1:row_length,1:rows,:),                            &
             grid_airmass, z_top_of_model, shno3_3d, nat_psc,                  &
             all_tracers(:,:,:,1:n_chem_tracers+n_aero_tracers))
    END IF

    IF (ALLOCATED(t_chem)) DEALLOCATE(t_chem)
    IF (ALLOCATED(q_chem)) DEALLOCATE(q_chem)

    ! Add fields to section 50 diagnostics if requested
    IF (ukca_config%l_enable_diag_um .AND.                                     &
        (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be)) THEN
      CALL ukca_offline_oxidants_diags(row_length, rows, model_levels,         &
                                       totnodens, o3_offline_diag,             &
                                       oh_offline_diag, no3_offline_diag,      &
                                       ho2_offline_diag, h2o2_offline,         &
                                       SIZE(stashwork50), stashwork50)
      DEALLOCATE(o3_offline_diag)
      DEALLOCATE(oh_offline_diag)
      DEALLOCATE(no3_offline_diag)
      DEALLOCATE(ho2_offline_diag)
    END IF

  ELSE

    ! The current time step is not a chemistry time step so set the
    ! return status of any requests for diagnostics that are only output
    ! on chemistry time steps as skipped
    CALL update_skipped_diag_flags(diagnostics)

  END IF ! do_chemistry

  IF (n_strat_fluxdiags > 0) THEN
    DO l=1,n_chem_tracers
      DO k=1,model_levels
        DO j=1,rows
          DO i=1,row_length
            trmol_post_chem(i,j,k,l) = all_tracers(i,j,k,l)                    &
                                *totnodens(i,j,k)*grid_volume(i,j,k)           &
                                /(c_species(l)*avogadro)    ! moles
          END DO
        END DO
      END DO
    END DO
  END IF   ! End of IF n_strat_fluxdiags > 0 statement

  !  Save all_tracers array for STE calculation in ASAD Flux Diagnostics
  IF (.NOT. ALLOCATED(fluxdiag_all_tracers) .AND.                              &
       (L_asad_use_chem_diags .AND. L_asad_use_STE))                           &
       THEN
    ALLOCATE(fluxdiag_all_tracers(row_length,rows,model_levels,n_chem_tracers))
    fluxdiag_all_tracers = all_tracers(:,:,:,1:n_chem_tracers)
  END IF

  l_store_value=.FALSE.
  IF (do_chemistry .AND. (L_asad_use_chem_diags                                &
       .AND. L_asad_use_tendency))                                             &
    CALL asad_tendency_ste(row_length,rows,model_levels,                       &
                           n_chem_tracers,                                     &
                           all_tracers(:,:,:,1:n_chem_tracers),                &
                           totnodens,grid_volume,ukca_config%timestep,         &
                           calculate_tendency,l_store_value,ierr)

  ! Transform halogen/nitrogen/hydrogen species back
  IF (ukca_config%l_tracer_lumping) THEN

    ! first, copy unlumped NO2, BrO and HCl fields into all_ntp structure

    ! NO2
    i = name2ntpindex('NO2       ')
    all_ntp(i)%data_3d(:,:,:) = all_tracers(:,:,:,n_no2)

    ! BrO
    i = name2ntpindex('BrO       ')
    all_ntp(i)%data_3d(:,:,:) = all_tracers(:,:,:,n_bro)

    ! HCl
    i = name2ntpindex('HCl       ')
    all_ntp(i)%data_3d(:,:,:) = all_tracers(:,:,:,n_hcl)

    ! Now call ukca_transform halogen to do the lumping. After this,
    ! the all_tracers array contains the lumped species not NO2, BrO and HCl
    CALL ukca_transform_halogen(rows,row_length,model_levels,n_tracers,        &
                                lump_species,all_tracers)
  END IF

#if !defined(LFRIC)
  IF (ukca_config%l_timer) CALL timer('UKCA CHEMISTRY MODEL',6)
#endif
END IF    ! End of IF (l_ukca_chem) statement

! ----------------------------------------------------------------------
! 4.6 GLOMAP-mode aerosol scheme
! ----------------------------------------------------------------------
IF (ukca_config%l_ukca_mode .AND. do_chemistry) THEN

  ! Allocate space for copies of fields required for later diagnostic
  ! calculations
  IF (l_ukca_cmip6_diags .OR. l_ukca_pm_diags)                                 &
    CALL ukca_mode_diags_alloc(theta_field_size*model_levels)

  ! Cache blocking scheme
  stride_seg = row_length*rows    ! i.e. number of columns
#if !defined(LFRIC)
  IF (l_autotune_segments) THEN
    glomap_config%ukca_mode_seg_size =                                         &
      MIN(glomap_config%ukca_mode_seg_size, stride_seg)
  END IF
#endif
  nbox_max = glomap_config%ukca_mode_seg_size*model_levels
                                                     ! requested cols by levels
  IF ( (glomap_config%ukca_mode_seg_size < 1) .OR.                             &
       (glomap_config%ukca_mode_seg_size > stride_seg) ) THEN
    ! unrealistic value for ukca_mode_seg_size more than columns on MPI task
    WRITE(cmessage,'(A,i8,A,i8)') 'Select ukca_mode_seg_size as a factor of'// &
      ' stride_seg. Current values are: ukca_mode_seg_size=',                  &
      glomap_config%ukca_mode_seg_size, ' and stride_seg=',stride_seg
    errcode = 6
    CALL ereport(RoutineName,errcode,cmessage)
    ! code will abort (no way to test any earlier)
  ELSE
    ! divide number of total columns by the ncpc
    ! this allows segments to cross rows.
    nseg = stride_seg/glomap_config%ukca_mode_seg_size ! number of segments
                                                       ! (check remainder)
                                                       ! all uniform initially
    nbox_s(1:nseg) = glomap_config%ukca_mode_seg_size*model_levels
    nbox_s(nseg+1:stride_seg) = 0
    ncol_s(1:nseg) = glomap_config%ukca_mode_seg_size
    ncol_s(nseg+1:stride_seg) = 0

    seg_rem = MOD(stride_seg,glomap_config%ukca_mode_seg_size)
                                                       ! Only allow whole
                                                       ! factors of cols
    IF (seg_rem > 0) THEN
      ! account for final additional segment
      nseg = nseg+1
      ncol_s(nseg) = seg_rem    ! a smaller number of columns
      nbox_s(nseg) = seg_rem*model_levels  ! fewer boxes on final segment
      IF ( l_first_call ) THEN
        WRITE(cmessage,'(A,i8,A,i8,A)')                                        &
          'ukca_mode_seg_size is not a factor of columns on this PE ',         &
          stride_seg,'. Segments will be non-uniform with final segment of ',  &
          seg_rem,' columns.'
        errcode = -6
        CALL ereport(RoutineName,errcode,cmessage)
      END IF  ! l_first_call
    END IF

    ! determine start of each segment in terms of column base
    lb = 1
    DO ik = 1, nseg
      lbase(ik) = lb
      lb = lb + ncol_s(ik)
    END DO

  END IF

#if !defined(LFRIC)
  IF (ukca_config%l_timer) CALL timer('UKCA AEROSOL MODEL  ',5)

  IF (l_autotune_segments) THEN
    CALL autotune_start_region(autotune_state, stride_seg)
  END IF
#endif

  !!!! In the call below, note that cloud_frac field may have lower bound of 0
  !!!! hence the subrange 1:model_levels is specified. This relates to a
  !!!! bug controlled by the temporary logical l_fix_ukca_cloud_frac and can be
  !!!! removed once the bug fix is no longer conditional on this logical
  CALL ukca_aero_ctl(i_month, i_day_number, i_hour,                            &
       INT(r_minute - ukca_config%timestep/60.0),                              &
       REAL(ukca_config%chem_timestep),                                        &
       row_length, rows, model_levels,                                         &
       n_chem_tracers+n_aero_tracers,                                          &
       n_mode_tracers,                                                         &
       p_theta_levels,                                                         &
       t_theta_levels,                                                         &
       q,                                                                      &
       rel_humid_frac,                                                         &
       rel_humid_frac_clr,                                                     &
       p_layer_boundaries,                                                     &
       all_tracers_names(1:n_chem_tracers+n_aero_tracers),                     &
       all_tracers_names(                                                      &
            n_chem_tracers+n_aero_tracers+1:                                   &
            n_chem_tracers+n_aero_tracers+                                     &
            n_mode_tracers),                                                   &
       all_tracers(:,:,:,1:n_chem_tracers+n_aero_tracers),                     &
       all_tracers(:,:,:,                                                      &
            n_chem_tracers+n_aero_tracers+1:                                   &
            n_chem_tracers+n_aero_tracers+                                     &
            n_mode_tracers),                                                   &
       seaice_frac,                                                            &
       rough_length,                                                           &
       u_s,                                                                    &
       ls_rain3d,                                                              &
       conv_rain3d,                                                            &
       ls_snow3d,                                                              &
       conv_snow3d,                                                            &
       autoconv,                                                               &
       accretion,                                                              &
       rim_agg,                                                                &
       rim_cry,                                                                &
       land_fraction,                                                          &
       nbox_max,                                                               &
       delso2_wet_h2o2,                                                        &
       delso2_wet_o3,                                                          &
       delh2so4_chem,                                                          &
       mode_diags,                                                             &
       cloud_frac(:,:,1:model_levels),                                         &
       cloud_liq_frac,                                                         &
       qcl,                                                                    &
       z_half_alllevs,                                                         &
       grid_airmass,zbl,                                                       &
       uph2so4inaer,                                                           &
       wetox_in_aer,                                                           &
       all_ntp,                                                                &
       nseg, nbox_s, ncol_s, lbase, stride_seg                                 &
       )

#if !defined(LFRIC)
  !If autotuning is active, decide what to do with the
  !segment size and report the current status.
  IF (l_autotune_segments) THEN
    CALL autotune_stop_region(autotune_state)
  END IF

  IF (ukca_config%l_timer) CALL timer('UKCA AEROSOL MODEL  ',6)
#endif

  ! Call activation scheme if switched on
  IF ( glomap_config%i_ukca_activation_scheme == i_ukca_activation_arg ) THEN

    CALL ukca_activate(                                                        &
        row_length, rows,                                                      &
        model_levels,                                                          &
        ukca_config%bl_levels,                                                 &
        theta_field_size,                                                      &
        n_mode_tracers,                                                        &
        n_mode_diags,                                                          &
        glomap_config%i_ukca_nwbins,                                           &
        p_theta_levels,                                                        &
        t_theta_levels,                                                        &
        q,                                                                     &
        qsvp,                                                                  &
        bl_tke,                                                                &
        vertvel,                                                               &
        cloud_liq_frac,                                                        &
        qcl,                                                                   &
        cdncflag,                                                              &
        n_activ_sum,                                                           &
        cdncwt,                                                                &
        all_tracers(:,:,:,                                                     &
                    n_chem_tracers+n_aero_tracers+1:                           &
                    n_chem_tracers+n_aero_tracers+n_mode_tracers),             &
        mode_diags,                                                            &
        glomap_variables                                                       &
        )

    ! Write CDNC to all_ntp structure
    i = name2ntpindex('cdnc      ')
    all_ntp(i)%data_3d(:,:,:)=cdncflag(:,:,:)
    IF (glomap_config%l_ntpreq_n_activ_sum) THEN
      ! Write the total activated aerosol if required
      i = name2ntpindex('n_activ_sum   ')
      all_ntp(i)%data_3d(:,:,:)=n_activ_sum(:,:,:)
    END IF

  END IF   ! ( i_ukca_activation_scheme == i_ukca_activation_arg )

END IF   ! ukca_config%l_ukca_mode

! ----------------------------------------------------------------------
! 5. Output prognostics and diagnostics
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 5.1 Copy prognostics
! ----------------------------------------------------------------------
! Take a copy of tracer fields if requested for diagnostic purposes
IF (L_asad_use_chem_diags .AND. L_asad_use_output_tracer)                      &
     CALL asad_output_tracer(row_length, rows, model_levels,                   &
                             n_chem_tracers, all_tracers, ierr)

! ----------------------------------------------------------------------
! 5.2 Service diagnostic requests.
! ----------------------------------------------------------------------
! Requests processed in this section can be UKCA diagnostics requests
! received via the UKCA API or legacy-style requests received via the
! UM STASH system. The latter are only supported when UKCA is coupled
! with the UM parent model and are serviced by copying the requested
! data to the STASH work arrays rather than passing them back via the
! argument list.
! ----------------------------------------------------------------------

IF (ukca_config%l_enable_diag_um .AND. ukca_config%l_ukca_chem) THEN

  ! ----------------------------------------------------------------------
  ! 5.2.1 MODE diagnostics [UM legacy-style requests]
  !       (items 38,201 - 38,212 now done in ukca_mode_ems_um)
  !       Now adding nitrate diagnostics. 584-645 (item1_nitrate_diags to
  !       itemN_nitrate_diags) are the fluxes and partial volumes 646-668
  !       are the CMIP diagnostics which are not accounted for here.
  !       Now adding diagnostics for the super-coarse insoluble mode (dust)
  !       which was added after the other modes.
  ! ----------------------------------------------------------------------

  icnt=0
  DO l=1,nmax_mode_diags
    IF (UkcaD1codes(imode_first+l-1)%item /= imdi) THEN
      icnt = icnt + 1
      item = UkcaD1codes(imode_first+l-1)%item
      section = stashcode_glomap_sec
      IF (sf(item,section) .AND. item > item1_mode_diags+12 .AND.              &
          item < item1_mode_diags+284 .AND. item /= 388 .AND.                  &
          item /= 389) THEN
        CALL copydiag_3d (stashwork38(si(item,section,im_index):               &
          si_last(item,section,im_index)),                                     &
          mode_diags(:,:,:,icnt),                                              &
          row_length,rows,model_levels,                                        &
          stlist(:,stindex(1,item,section,im_index)),len_stlist,               &
          stash_levels,num_stash_levels+1)
      END IF
    END IF
  END DO       ! 1,nmax_mode_diags

  n_nitrate_diags = itemN_nitrate_diags -  item1_nitrate_diags + 1

  DO l=1,n_nitrate_diags
    IF (UkcaD1codes(imode_first+nmax_mode_diags+l-1)%item /= imdi) THEN
      icnt = icnt + 1
      item = UkcaD1codes(imode_first+nmax_mode_diags+l-1)%item
      section = stashcode_glomap_sec
      IF (sf(item,section) .AND. item >= item1_nitrate_diags  .AND.            &
          item <= item1_nitrate_diags+61) THEN
        CALL copydiag_3d (stashwork38(si(item,section,im_index):               &
             si_last(item,section,im_index)),                                  &
          mode_diags(:,:,:,icnt),                                              &
          row_length,rows,model_levels,                                        &
          stlist(:,stindex(1,item,section,im_index)),len_stlist,               &
          stash_levels,num_stash_levels+1)
      END IF
    END IF
  END DO       ! 1,n_nitrate_diags

  nmax_diags_inc_nitr = nmax_mode_diags + n_nitrate_diags
  n_sup_dust_diags = itemN_dust3mode_diags - item1_dust3mode_diags + 1

  DO l=1,n_sup_dust_diags
    IF (UkcaD1codes(imode_first+nmax_diags_inc_nitr+l-1)%item /= imdi) THEN
      icnt = icnt + 1
      item = UkcaD1codes(imode_first+nmax_diags_inc_nitr+l-1)%item
      section = stashcode_glomap_sec
      IF (sf(item,section) .AND. item >= item1_dust3mode_diags+1  .AND.        &
          item < item1_dust3mode_diags+19) THEN
        CALL copydiag_3d (stashwork38(si(item,section,im_index):               &
             si_last(item,section,im_index)),                                  &
          mode_diags(:,:,:,icnt),                                              &
          row_length,rows,model_levels,                                        &
          stlist(:,stindex(1,item,section,im_index)),len_stlist,               &
          stash_levels,num_stash_levels+1)
      END IF
    END IF
  END DO       ! 1,n_sup_dust_diags

  ! Copy CMIP6 diagnostics and/or PM diagnostics into STASHwork array
  IF (ukca_config%l_ukca_mode .AND. do_chemistry .AND.                         &
      (l_ukca_cmip6_diags .OR. l_ukca_pm_diags)) THEN
    CALL ukca_mode_diags(row_length, rows, model_levels,                       &
                         theta_field_size*model_levels,                        &
                         n_mode_tracers,                                       &
                         p_theta_levels,                                       &
                         t_theta_levels,                                       &
                         all_tracers(:,:,:,                                    &
                                     n_chem_tracers+n_aero_tracers+1:          &
                                     n_chem_tracers+n_aero_tracers+            &
                                     n_mode_tracers),                          &
                         interf_z,                                             &
                         SIZE(stashwork38),                                    &
                         stashwork38)
  END IF

  ! ----------------------------------------------------------------------
  ! 5.2.2 Stratospheric flux diagnostics [UM legacy-style requests]
  ! ----------------------------------------------------------------------
  icnt = 0
  DO l=1,nmax_strat_fluxdiags
    IF (UkcaD1codes(istrat_first+l-1)%item /= imdi) THEN
      icnt = icnt + 1
      item = UkcaD1codes(istrat_first+l-1)%item
      section = stashcode_glomap_sec
      IF (sf(item,section)) THEN
        CALL copydiag_3d (stashwork38(si(item,section,im_index):               &
          si_last(item,section,im_index)),                                     &
          strat_fluxdiags(:,:,:,icnt),                                         &
          row_length,rows,model_levels,                                        &
          stlist(:,stindex(1,item,section,im_index)),len_stlist,               &
          stash_levels,num_stash_levels+1)
      ELSE
        cmessage=' Strat flux item not found by STASH flag array'
        icode = section*1000+item
        CALL ereport(RoutineName,icode,cmessage)
      END IF
    END IF
  END DO       ! 1,nmax_strat_fluxdiags

  ! ----------------------------------------------------------------------
  ! 5.2.3 Diagnostics calculated from arrays passed back from
  ! chemistry_ctl [UM legacy-style requests]
  ! ----------------------------------------------------------------------

  IF (do_chemistry)                                                            &
    CALL ukca_chem_diags(                                                      &
      row_length, rows, model_levels,                                          &
      nat_psc,                                                                 &
      trop_ch4_mol,                                                            &
      trop_o3_mol,                                                             &
      trop_oh_mol,                                                             &
      strat_ch4_mol,                                                           &
      strat_ch4loss,                                                           &
      atm_ch4_mol,                                                             &
      atm_co_mol,                                                              &
      atm_n2o_mol,                                                             &
      atm_cf2cl2_mol,                                                          &
      atm_cfcl3_mol,                                                           &
      atm_mebr_mol,                                                            &
      atm_h2_mol,                                                              &
      so4_sa,                                                                  &
      H_plus_3d_arr,                                                           &
      SIZE(stashwork50), stashwork50)

END IF ! ukca_config%l_enable_diag_um .AND. ukca_config%l_ukca_chem

! ---------------------------------------------------------------------
! 5.2.4 Service any requests for available non-ASAD diagnostics
!       recognised by UKCA's diagnostic handling system.
!       [UKCA API requests]
! ---------------------------------------------------------------------

IF (ukca_config%l_ukca_chem) THEN

  CALL ukca_diags_output_ctl(error_code_ptr, row_length, rows, model_levels,   &
                             n_use_tracers, z_top_of_model,                    &
                             do_chemistry, p_tropopause, p_layer_boundaries,   &
                             p_theta_levels, plumeria_height, all_tracers,     &
                             photol_rates, diagnostics,                        &
                             error_message=error_message,                      &
                             error_routine=error_routine)

  IF (error_code_ptr > 0) THEN
    IF (lhook) THEN
      CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    END IF
    RETURN
  END IF

END IF

! ----------------------------------------------------------------------
! 5.2.5 ASAD flux diagnostics [UKCA API and/or UM legacy-style requests]
!        - needs to be done BEFORE ste calc below
! ----------------------------------------------------------------------

IF (L_asad_use_chem_diags) THEN

  stashsize = SIZE(stashwork50)

  CALL asad_diags_output_ctl(error_code_ptr,                                   &
                             row_length, rows, model_levels, stashsize,        &
                             diagnostics, stashwork50,                         &
                             error_message=error_message,                      &
                             error_routine=error_routine)

  IF (error_code_ptr > 0) THEN
    IF (lhook)                                                                 &
      CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
    RETURN
  END IF

END IF

IF (ukca_config%l_enable_diag_um .AND. ukca_config%l_ukca_chem) THEN

  ! ----------------------------------------------------------------------
  ! 5.2.6 ASAD Flux Diags Strat-Trop Exchange [UM legacy-style requests]
  !        - done AFTER asad_diags_output_ctl
  ! ----------------------------------------------------------------------
  l_store_value = .TRUE.
  IF (L_asad_use_chem_diags .AND. L_asad_use_STE) THEN
    CALL asad_tendency_ste(                                                    &
                 row_length, rows, model_levels, n_chem_tracers,               &
                 fluxdiag_all_tracers, totnodens, grid_volume,                 &
                 ukca_config%timestep, calculate_STE, l_store_value, ierr)
  END IF

  ! ----------------------------------------------------------------------
  ! 5.2.7 Grid-cell volume [UM legacy style request]
  ! ----------------------------------------------------------------------
  section = ukca_diag_sect
  item = 255
  IF (sf(item,ukca_diag_sect) .AND. ALLOCATED(grid_volume)) THEN
    CALL copydiag_3d (stashwork50(si(item,section,im_index):                   &
          si_last(item,section,im_index)),                                     &
          grid_volume,                                                         &
          row_length,rows,model_levels,                                        &
          stlist(:,stindex(1,item,section,im_index)),len_stlist,               &
          stash_levels,num_stash_levels+1)
  END IF
END IF

! ----------------------------------------------------------------------
! 6. Finally, blank out any missing diagnostics and deallocate arrays
! ----------------------------------------------------------------------

! Overwrite any non-valid diagnostic output if configured to do so
IF (ukca_config%l_blankout_invalid_diags) THEN
  CALL blank_out_missing_diags(error_code_ptr, diagnostics,                    &
                               error_message=error_message,                    &
                               error_routine=error_routine)
  IF (error_code_ptr > 0) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Reset/deallocate environmental driver fields
CALL clear_environment_fields()

! NB: have_nat3d is only allocated if chemistry is on and should be
!     deallocated at the end of every timestep regardless whether
!     persistence is turned off
IF (ALLOCATED(have_nat3d)) DEALLOCATE(have_nat3d)

! Deallocate global 3D pH array
IF (ALLOCATED(H_plus_3d_arr)) DEALLOCATE(H_plus_3d_arr)

! The following deallocations relate to the removal of persistent spatial
! fields in UKCA
IF (ukca_config%l_ukca_persist_off) THEN
  ! Deallocate zenith angle arrays
  CALL dealloc_diurnal_oxidant()
  ! Deallocate arrays allocated in ukca_main
  IF (ALLOCATED(int_zenith_angle)) DEALLOCATE(int_zenith_angle)
  IF (ALLOCATED(z_half_alllevs))   DEALLOCATE(z_half_alllevs)
  IF (ALLOCATED(z_half))           DEALLOCATE(z_half)
  IF (ALLOCATED(L_stratosphere))   DEALLOCATE(L_stratosphere)
  IF (ALLOCATED(pv_trop))          DEALLOCATE(pv_trop)
  IF (ALLOCATED(theta_trop))       DEALLOCATE(theta_trop)
  IF (ALLOCATED(tropopause_level)) DEALLOCATE(tropopause_level)
  IF (ALLOCATED(p_tropopause))     DEALLOCATE(p_tropopause)
END IF

! Deallocate variables local to ukca_main

IF (ALLOCATED(qsvp)) DEALLOCATE(qsvp)
IF (ALLOCATED(rel_humid_frac_clr)) DEALLOCATE(rel_humid_frac_clr)
IF (ALLOCATED(rel_humid_frac)) DEALLOCATE(rel_humid_frac)

IF (ALLOCATED(totnodens)) DEALLOCATE(totnodens)
IF (ALLOCATED(strat_fluxdiags)) DEALLOCATE(strat_fluxdiags)
IF (ALLOCATED(mode_diags)) DEALLOCATE(mode_diags)
IF (ALLOCATED(cos_zenith_angle)) DEALLOCATE(cos_zenith_angle)
IF (ALLOCATED(t_chem)) DEALLOCATE(t_chem)
IF (ALLOCATED(q_chem)) DEALLOCATE(q_chem)
IF (ALLOCATED(fluxdiag_all_tracers)) DEALLOCATE(fluxdiag_all_tracers)
IF (ALLOCATED(env_ozone3d)) DEALLOCATE(env_ozone3d)
IF (ALLOCATED(ls_ppn3d)) DEALLOCATE(ls_ppn3d)
IF (ALLOCATED(conv_ppn3d)) DEALLOCATE(conv_ppn3d)
IF (ALLOCATED(Thick_bl_levels)) DEALLOCATE(Thick_bl_levels)
IF (ALLOCATED(p_layer_boundaries)) DEALLOCATE(p_layer_boundaries)
IF (ALLOCATED(plumeria_height)) DEALLOCATE(plumeria_height)
IF (ALLOCATED(t_theta_levels)) DEALLOCATE(t_theta_levels)
IF (ALLOCATED(delSO2_wet_h2o2)) DEALLOCATE(delSO2_wet_h2o2)
IF (ALLOCATED(delSO2_wet_o3)) DEALLOCATE(delSO2_wet_o3)
IF (ALLOCATED(delh2so4_chem)) DEALLOCATE(delh2so4_chem)
IF (ALLOCATED(delSO2_drydep)) DEALLOCATE(delSO2_drydep)
IF (ALLOCATED(delSO2_wetdep)) DEALLOCATE(delSO2_wetdep)
IF (ALLOCATED(so4_sa)) DEALLOCATE(so4_sa)

! Deallocate ASAD arrays if persistence is off
IF (ukca_config%l_ukca_chem .AND. ukca_config%l_ukca_persist_off) THEN
  CALL ukca_delasad_spatial_vars()
END IF

l_first_call=.FALSE.
IF ( do_chemistry ) l_firstchem = .FALSE.

#if !defined(LFRIC)
IF (l_autotune_segments) THEN
  CALL autotune_return(autotune_state)
END IF
#endif

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_main1

END MODULE ukca_main1_mod
