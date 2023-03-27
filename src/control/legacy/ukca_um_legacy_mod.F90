! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module to satisfy legacy UM dependencies in UKCA.
!
! Method:
!
!   This module provides UKCA replacements for UM variables, parameters and
!   procedures that are required by UKCA when running in the UM and are
!   not currently provided via the UKCA API.
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

MODULE ukca_um_legacy_mod

! Used items below are declared elsewhere but accessed via this module to
! allow them to be replaced by UM versions in UM builds (by using a
! UM-specific version of this module). This may be required for consistency
! or simply for bit-reproducibility of results.
! In future, it should be possible to provide parent-specific values of
! these items via the 'ukca_setup' argument list whenever UKCA defaults are
! unsuitable.
! An item's use via this module can then be discontinued with a set of linked
! UM-LFRic-UKCA changes that includes the provision of any non-default values
! required using the new 'ukca_setup' functionality.

USE ukca_config_specification_mod, ONLY:                                       &
  isec_per_day,                                                                &
  isec_per_hour,                                                               &
  rhour_per_day,                                                               &
  rsec_per_day,                                                                &
  rsec_per_hour,                                                               &
  pi,                                                                          &
  pi_over_180,                                                                 &
  recip_pi_over_180,                                                           &
  zerodegc,                                                                    &
  rmol,                                                                        &
  planet_radius,                                                               &
  g,                                                                           &
  r,                                                                           &
  cp,                                                                          &
  pref,                                                                        &
  vkman,                                                                       &
  repsilon,                                                                    &
  rv,                                                                          &
  kappa,                                                                       &
  c_virtual,                                                                   &
  s2r,                                                                         &
  rad_ait,                                                                     &
  rad_acc,                                                                     &
  chi,                                                                         &
  sigma,                                                                       &
  i_liss_merlivat,                                                             &
  i_wanninkhof,                                                                &
  i_nightingale

USE ukca_dms_flux_mod, ONLY: dms_flux_4a => ukca_dms_flux

IMPLICIT NONE

PUBLIC

! Flags to indicate UM code availability
LOGICAL, PARAMETER :: l_um_infrastructure = .FALSE.
  ! UM infrastructure code including STASH support, grid parameters
LOGICAL, PARAMETER :: l_um_atmos_ukca_humidity = .FALSE.
  ! atmos_ukca_humidity subroutine
LOGICAL, PARAMETER :: l_um_emissions_updates = .FALSE.
  ! trsrc subroutine for emissions updating of tracers
LOGICAL, PARAMETER :: l_um_calc_surf_area = .FALSE.
  ! calc_surf_area subroutine

! ----------------------------------------------------------------------
! -- UM Parameters and Variables --
! ----------------------------------------------------------------------
! These items relate to UM functionality but are still required for
! compiling UKCA in non-UM applications. The values set here are either
! to dummy values or to fixed values appropriate for non-UM applications.

! Required for ukcad1codes (in this module) and by asad_chem_flux_diags
TYPE :: code
  INTEGER :: section
  INTEGER :: item
  INTEGER :: n_levels
  INTEGER :: address
  INTEGER :: length
  INTEGER :: halo_type
  INTEGER :: grid_type
  INTEGER :: field_type
  INTEGER :: len_dim1
  INTEGER :: len_dim2
  INTEGER :: len_dim3
  LOGICAL :: prognostic
  LOGICAL :: required
END TYPE code

! Required by ukca_environment_req_mod and ukca_ddcalc_mod
TYPE :: dummy_array_dims_type
  INTEGER :: i_end = 0   ! Dummy value
  INTEGER :: j_end = 0   ! Dummy value
  INTEGER :: k_start = 1
END TYPE
TYPE(dummy_array_dims_type) :: pdims, tdims

! Required by ukca_main1_mod

INTEGER, PARAMETER :: rh_z_top_theta = 0     ! Dummy value
INTEGER, PARAMETER :: global_row_length = 0  ! Dummy value

INTEGER, PARAMETER :: datastart(2) = [0,0]   ! Dummy value

REAL, PARAMETER :: delta_lambda = 0.0        ! Dummy value
REAL, PARAMETER :: delta_phi = 0.0           ! Dummy value
REAL, PARAMETER :: base_phi = 0.0            ! Dummy value

REAL, ALLOCATABLE :: a_realhd(:)
REAL, ALLOCATABLE :: fv_cos_theta_latitude(:,:)
REAL :: rhcrit(1) = [0]                      ! Dummy value

! Required by ukca_main1_mod and ukca_chemistry_ctl_col_mod
LOGICAL, PARAMETER :: l_autotune_segments = .FALSE.
TYPE :: autotune_type
  INTEGER :: dummy_value
END TYPE

! Required by ukca_add_emiss_mod
LOGICAL, PARAMETER :: l_ukca_so2ems_expvolc = .FALSE.

! Required by ukca_aer_no3_mod, ukca_prod_no3_mod and/or ukca_setup_mod;
! relate to use of CLASSIC dust

REAL, PARAMETER :: rhop = 1.0  ! Dummy value (must be non-zero)

REAL, ALLOCATABLE :: drep(:)

LOGICAL, PARAMETER :: l_dust = .FALSE.
LOGICAL, PARAMETER :: l_twobin_dust = .FALSE.

! Required by asad_inrats_mod
INTEGER, PARAMETER :: mype = 0  ! Serial case

! Required by ukca_activate_mod, ukca_mode_setup, ukca_vapour_mod and/or
! ukca_volume_mode_mod

INTEGER, PARAMETER :: i_gc_activation_arg = 0              ! Dummy value
INTEGER, PARAMETER :: i_glomap_clim_activation_scheme = 0  ! Dummy value
INTEGER, PARAMETER :: i_glomap_clim_setup = 0              ! Dummy value

LOGICAL, PARAMETER :: l_glomap_clim_radaer = .FALSE.
INTEGER, PARAMETER :: i_glomap_clim_tune_bc = 0

! Required by ukca_activate
REAL, POINTER :: gc_acc_sol_bc(:,:,:)
REAL, POINTER :: gc_acc_sol_oc(:,:,:)
REAL, POINTER :: gc_acc_sol_ss(:,:,:)
REAL, POINTER :: gc_acc_sol_su(:,:,:)
REAL, POINTER :: gc_ait_ins_bc(:,:,:)
REAL, POINTER :: gc_ait_ins_oc(:,:,:)
REAL, POINTER :: gc_ait_sol_bc(:,:,:)
REAL, POINTER :: gc_ait_sol_oc(:,:,:)
REAL, POINTER :: gc_ait_sol_su(:,:,:)
REAL, POINTER :: gc_cor_sol_bc(:,:,:)
REAL, POINTER :: gc_cor_sol_oc(:,:,:)
REAL, POINTER :: gc_cor_sol_ss(:,:,:)
REAL, POINTER :: gc_cor_sol_su(:,:,:)
REAL, POINTER :: gc_nd_acc_sol(:,:,:)
REAL, POINTER :: gc_nd_ait_ins(:,:,:)
REAL, POINTER :: gc_nd_ait_sol(:,:,:)
REAL, POINTER :: gc_nd_cor_sol(:,:,:)

! Required by asad_chem_flux_diags

INTEGER, PARAMETER :: st_levels_model_theta = 0  ! Dummy value
INTEGER, PARAMETER :: st_levels_single = 0       ! Dummy value
INTEGER, PARAMETER :: submodel_for_sm(1) = [0]   ! Dummy value
INTEGER, PARAMETER :: atmos_im = 1               ! Index for submodel_for_sm
INTEGER, PARAMETER :: ndiag = 0

INTEGER :: idom_b(1) = [0]  ! Dummy value
INTEGER :: iopl_d(1) = [0]  ! Dummy value
INTEGER :: isec_b(1) = [0]  ! Dummy value
INTEGER :: item_b(1) = [0]  ! Dummy value
INTEGER :: modl_b(1) = [0]  ! Dummy value

! Required by various modules for diagnostic processing

INTEGER            :: nukca_d1items = 0
INTEGER, PARAMETER :: istrat_first = 0         ! Dummy value
INTEGER, PARAMETER :: imode_first = 0          ! Dummy value
INTEGER, PARAMETER :: ukca_diag_sect = 0       ! Dummy value
INTEGER, PARAMETER :: item1_mode_diags = 0     ! Dummy value
INTEGER, PARAMETER :: item1_nitrate_diags = 0  ! Dummy value
INTEGER, PARAMETER :: itemN_nitrate_diags = 0  ! Dummy value
INTEGER            :: n_strat_fluxdiags = 0
INTEGER            :: n_mode_diags = 0
INTEGER, PARAMETER :: len_stlist = 0           ! Dummy value
INTEGER, PARAMETER :: num_stash_levels = 0     ! Dummy value
INTEGER, PARAMETER :: num_stash_pseudo = 0     ! Dummy value

INTEGER, ALLOCATABLE :: stindex(:,:,:,:)
INTEGER, ALLOCATABLE :: stlist(:,:)
INTEGER, ALLOCATABLE :: si(:,:,:)
INTEGER, ALLOCATABLE :: si_last(:,:,:)
INTEGER, ALLOCATABLE :: stash_levels(:,:)
INTEGER, ALLOCATABLE :: stash_pseudo_levels(:,:)

LOGICAL :: l_ukca_stratflux = .FALSE.
LOGICAL :: l_ukca_mode_diags = .FALSE.

LOGICAL, ALLOCATABLE :: sf(:,:)

TYPE(code), ALLOCATABLE, SAVE :: ukcad1codes(:)

! STASH codes, required by various modules (mainly diagnostics)
! All have dummy values here
INTEGER, PARAMETER :: stashcode_bc_acc_sol = 0
INTEGER, PARAMETER :: stashcode_bc_acc_sol_load = 0
INTEGER, PARAMETER :: stashcode_bc_ait_insol = 0
INTEGER, PARAMETER :: stashcode_bc_ait_insol_load = 0
INTEGER, PARAMETER :: stashcode_bc_ait_sol = 0
INTEGER, PARAMETER :: stashcode_bc_ait_sol_load = 0
INTEGER, PARAMETER :: stashcode_bc_cor_sol = 0
INTEGER, PARAMETER :: stashcode_bc_cor_sol_load = 0
INTEGER, PARAMETER :: stashcode_bc_total_load = 0
INTEGER, PARAMETER :: stashcode_du_acc_insol = 0
INTEGER, PARAMETER :: stashcode_du_acc_insol_load = 0
INTEGER, PARAMETER :: stashcode_du_acc_sol = 0
INTEGER, PARAMETER :: stashcode_du_acc_sol_load = 0
INTEGER, PARAMETER :: stashcode_du_cor_insol = 0
INTEGER, PARAMETER :: stashcode_du_cor_insol_load = 0
INTEGER, PARAMETER :: stashcode_du_cor_sol = 0
INTEGER, PARAMETER :: stashcode_du_cor_sol_load = 0
INTEGER, PARAMETER :: stashcode_du_total_load = 0
INTEGER, PARAMETER :: stashcode_glomap_sec = 0
INTEGER, PARAMETER :: stashcode_h2o_acc_sol = 0
INTEGER, PARAMETER :: stashcode_h2o_acc_sol_load = 0
INTEGER, PARAMETER :: stashcode_h2o_ait_sol = 0
INTEGER, PARAMETER :: stashcode_h2o_ait_sol_load = 0
INTEGER, PARAMETER :: stashcode_h2o_cor_sol = 0
INTEGER, PARAMETER :: stashcode_h2o_cor_sol_load = 0
INTEGER, PARAMETER :: stashcode_h2o_mmr = 0
INTEGER, PARAMETER :: stashcode_h2o_nuc_sol = 0
INTEGER, PARAMETER :: stashcode_h2o_nuc_sol_load = 0
INTEGER, PARAMETER :: stashcode_h2o_total = 0
INTEGER, PARAMETER :: stashcode_h2o_total_load = 0
INTEGER, PARAMETER :: stashcode_n_acc_insol = 0
INTEGER, PARAMETER :: stashcode_n_acc_sol = 0
INTEGER, PARAMETER :: stashcode_n_ait_insol = 0
INTEGER, PARAMETER :: stashcode_n_ait_sol = 0
INTEGER, PARAMETER :: stashcode_n_cor_insol = 0
INTEGER, PARAMETER :: stashcode_n_cor_sol = 0
INTEGER, PARAMETER :: stashcode_n_nuc_sol = 0
INTEGER, PARAMETER :: stashcode_nh4_acc_sol = 0
INTEGER, PARAMETER :: stashcode_nh4_acc_sol_load = 0
INTEGER, PARAMETER :: stashcode_nh4_ait_sol = 0
INTEGER, PARAMETER :: stashcode_nh4_ait_sol_load = 0
INTEGER, PARAMETER :: stashcode_nh4_cor_sol = 0
INTEGER, PARAMETER :: stashcode_nh4_cor_sol_load = 0
INTEGER, PARAMETER :: stashcode_nh4_total_load = 0
INTEGER, PARAMETER :: stashcode_nn_acc_sol = 0
INTEGER, PARAMETER :: stashcode_nn_acc_sol_load = 0
INTEGER, PARAMETER :: stashcode_nn_cor_sol = 0
INTEGER, PARAMETER :: stashcode_nn_cor_sol_load = 0
INTEGER, PARAMETER :: stashcode_nn_total_load = 0
INTEGER, PARAMETER :: stashcode_no3_acc_sol = 0
INTEGER, PARAMETER :: stashcode_no3_acc_sol_load = 0
INTEGER, PARAMETER :: stashcode_no3_ait_sol = 0
INTEGER, PARAMETER :: stashcode_no3_ait_sol_load = 0
INTEGER, PARAMETER :: stashcode_no3_cor_sol = 0
INTEGER, PARAMETER :: stashcode_no3_cor_sol_load = 0
INTEGER, PARAMETER :: stashcode_no3_total_load = 0
INTEGER, PARAMETER :: stashcode_oc_acc_sol = 0
INTEGER, PARAMETER :: stashcode_oc_acc_sol_load = 0
INTEGER, PARAMETER :: stashcode_oc_ait_insol = 0
INTEGER, PARAMETER :: stashcode_oc_ait_insol_load = 0
INTEGER, PARAMETER :: stashcode_oc_ait_sol = 0
INTEGER, PARAMETER :: stashcode_oc_ait_sol_load = 0
INTEGER, PARAMETER :: stashcode_oc_cor_sol = 0
INTEGER, PARAMETER :: stashcode_oc_cor_sol_load = 0
INTEGER, PARAMETER :: stashcode_oc_nuc_sol = 0
INTEGER, PARAMETER :: stashcode_oc_nuc_sol_load = 0
INTEGER, PARAMETER :: stashcode_oc_total_load = 0
INTEGER, PARAMETER :: stashcode_pm10_bc = 0
INTEGER, PARAMETER :: stashcode_pm10_dry = 0
INTEGER, PARAMETER :: stashcode_pm10_du = 0
INTEGER, PARAMETER :: stashcode_pm10_nh4 = 0
INTEGER, PARAMETER :: stashcode_pm10_nn = 0
INTEGER, PARAMETER :: stashcode_pm10_no3 = 0
INTEGER, PARAMETER :: stashcode_pm10_oc = 0
INTEGER, PARAMETER :: stashcode_pm10_so4 = 0
INTEGER, PARAMETER :: stashcode_pm10_ss = 0
INTEGER, PARAMETER :: stashcode_pm10_wet = 0
INTEGER, PARAMETER :: stashcode_pm2p5_bc = 0
INTEGER, PARAMETER :: stashcode_pm2p5_dry = 0
INTEGER, PARAMETER :: stashcode_pm2p5_du = 0
INTEGER, PARAMETER :: stashcode_pm2p5_nh4 = 0
INTEGER, PARAMETER :: stashcode_pm2p5_nn = 0
INTEGER, PARAMETER :: stashcode_pm2p5_no3 = 0
INTEGER, PARAMETER :: stashcode_pm2p5_oc = 0
INTEGER, PARAMETER :: stashcode_pm2p5_so4 = 0
INTEGER, PARAMETER :: stashcode_pm2p5_ss = 0
INTEGER, PARAMETER :: stashcode_pm2p5_wet = 0
INTEGER, PARAMETER :: stashcode_so4_acc_sol = 0
INTEGER, PARAMETER :: stashcode_so4_acc_sol_load = 0
INTEGER, PARAMETER :: stashcode_so4_ait_sol = 0
INTEGER, PARAMETER :: stashcode_so4_ait_sol_load = 0
INTEGER, PARAMETER :: stashcode_so4_cor_sol = 0
INTEGER, PARAMETER :: stashcode_so4_cor_sol_load = 0
INTEGER, PARAMETER :: stashcode_so4_nuc_sol = 0
INTEGER, PARAMETER :: stashcode_so4_nuc_sol_load = 0
INTEGER, PARAMETER :: stashcode_so4_total_load = 0
INTEGER, PARAMETER :: stashcode_ss_acc_sol = 0
INTEGER, PARAMETER :: stashcode_ss_acc_sol_load = 0
INTEGER, PARAMETER :: stashcode_ss_cor_sol = 0
INTEGER, PARAMETER :: stashcode_ss_cor_sol_load = 0
INTEGER, PARAMETER :: stashcode_ss_total_load = 0
INTEGER, PARAMETER :: stashcode_ukca_atmos_cfc11 = 0
INTEGER, PARAMETER :: stashcode_ukca_atmos_cfc12 = 0
INTEGER, PARAMETER :: stashcode_ukca_atmos_ch3br = 0
INTEGER, PARAMETER :: stashcode_ukca_atmos_ch4 = 0
INTEGER, PARAMETER :: stashcode_ukca_atmos_co = 0
INTEGER, PARAMETER :: stashcode_ukca_atmos_h2 = 0
INTEGER, PARAMETER :: stashcode_ukca_atmos_n2o = 0
INTEGER, PARAMETER :: stashcode_ukca_chem_diag = 0
INTEGER, PARAMETER :: stashcode_ukca_h_plus = 0
INTEGER, PARAMETER :: stashcode_ukca_nat = 0
INTEGER, PARAMETER :: stashcode_ukca_so4_sad = 0
INTEGER, PARAMETER :: stashcode_ukca_strat_ch4 = 0
INTEGER, PARAMETER :: stashcode_ukca_strt_ch4_lss = 0
INTEGER, PARAMETER :: stashcode_ukca_trop_ch4 = 0
INTEGER, PARAMETER :: stashcode_ukca_trop_o3 = 0
INTEGER, PARAMETER :: stashcode_ukca_trop_oh = 0

! ----------------------------------------------------------------------
! -- Interfaces for Dummy Procedures --
! ----------------------------------------------------------------------

INTERFACE copydiag
  MODULE PROCEDURE copydiag_real
  MODULE PROCEDURE copydiag_real_vector
  MODULE PROCEDURE copydiag_integer
  MODULE PROCEDURE copydiag_logical
END INTERFACE
INTERFACE copydiag_3d
  MODULE PROCEDURE copydiag_3d_real
  MODULE PROCEDURE copydiag_3d_logical
  MODULE PROCEDURE copydiag_3d_integer
END INTERFACE

CONTAINS

! ----------------------------------------------------------------------
! -- UKCA Replacements for UM Procedures in Use --
! ----------------------------------------------------------------------
! These are Vectlib replacement subroutines required by various modules

SUBROUTINE exp_v(n, x, y)

IMPLICIT NONE

! Sets y(i) to the exponential function of x(i), for i=1,..,n

INTEGER, INTENT(IN) :: n
REAL, INTENT(IN) :: x(n)
REAL, INTENT(OUT) :: y(n)

INTEGER :: i

DO i = 1, n
  y(i) = EXP(x(i))
END DO

RETURN
END SUBROUTINE exp_v

! ----------------------------------------------------------------------

SUBROUTINE powr_v(n, x, y, z)

IMPLICIT NONE

! Sets z(i) to x(i) raised to the power y, for i=1,..,n

INTEGER, INTENT(IN) :: n
REAL, INTENT(IN) :: x(n)
REAL, INTENT(IN) :: y
REAL, INTENT(OUT) :: z(n)

INTEGER :: i

DO i = 1, n
  z(i) = x(i) ** y
END DO

RETURN
END SUBROUTINE powr_v

! ----------------------------------------------------------------------

SUBROUTINE cubrt_v(n, x, y)

IMPLICIT NONE

! Sets y(i) to the cubed root of x(i), for i=1,..,n

INTEGER, INTENT(IN) :: n
REAL, INTENT(IN) :: x(n)
REAL, INTENT(OUT) :: y(n)

INTEGER :: i

DO i = 1, n
  y(i) = x(i) ** (1.0 / 3.0)
END DO

RETURN
END SUBROUTINE cubrt_v

! ----------------------------------------------------------------------

SUBROUTINE oneover_v(n, x, y)

IMPLICIT NONE

! Sets y(i) to the reciprocal of x(i), for i=1,..,n

INTEGER, INTENT(IN) :: n
REAL, INTENT(IN) :: x(n)
REAL, INTENT(OUT) :: y(n)

INTEGER :: i

DO i = 1, n
  y(i) = 1 / x(i)
END DO

RETURN
END SUBROUTINE oneover_v

! ----------------------------------------------------------------------

SUBROUTINE log_v(n, x, y)

IMPLICIT NONE

! Sets y(i) to the natural logarithm of x(i), for i=1,..,n

INTEGER, INTENT(IN) :: n
REAL, INTENT(IN) :: x(n)
REAL, INTENT(OUT) :: y(n)

INTEGER :: i

DO i = 1, n
  y(i) = LOG(x(i))
END DO

RETURN
END SUBROUTINE log_v

! ----------------------------------------------------------------------
! -- Unused Procedures (dummy versions) --
! ----------------------------------------------------------------------

! Subroutines for humidity-related calculations required by ukca_main1_mod

SUBROUTINE atmos_ukca_humidity(row_length, rows, model_levels,                 &
                               t_theta_levels, p_theta_levels, q, rh,          &
                               qsatmr_wat, qsatmr_ice_wat, svp)
IMPLICIT NONE
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
REAL, INTENT(IN) :: t_theta_levels(:,:,:)
REAL, INTENT(IN) :: p_theta_levels(:,:,:)
REAL, OPTIONAL, INTENT(IN) :: q(:,:,:)
REAL, OPTIONAL, INTENT(OUT) :: rh(:,:,:)
REAL, OPTIONAL, INTENT(OUT) :: qsatmr_wat(:,:,:)
REAL, OPTIONAL, INTENT(OUT) :: qsatmr_ice_wat(:,:,:)
REAL, OPTIONAL, INTENT(OUT) :: svp(:,:,:)
IF (PRESENT(rh)) rh = 0.0
IF (PRESENT(qsatmr_wat)) qsatmr_wat = 0.0
IF (PRESENT(qsatmr_ice_wat)) qsatmr_ice_wat = 0.0
IF (PRESENT(svp)) svp = 0.0
END SUBROUTINE atmos_ukca_humidity

SUBROUTINE lsp_qclear(q, qsmr, qsmr_wat, qcf, cloud_liq_frac, cloud_frac,      &
                      rhcrit, q_clear, npnts)
IMPLICIT NONE
REAL, INTENT(IN) :: q(:)
REAL, INTENT(IN) :: qsmr(:)
REAL, INTENT(IN) :: qsmr_wat(:)
REAL, INTENT(IN) :: qcf(:)
REAL, INTENT(IN) :: cloud_liq_frac(:)
REAL, INTENT(IN) :: cloud_frac(:)
REAL, INTENT(IN) :: rhcrit(:)
REAL, INTENT(OUT) :: q_clear(:)
INTEGER, INTENT(IN) :: npnts
q_clear = 0.0
RETURN
END SUBROUTINE lsp_qclear

! Tracer update subroutine required by ukca_add_emiss_mod

SUBROUTINE trsrce(rows, row_length, offx, offy, halo_i, halo_j,                &
                  r_theta_levels, r_rho_levels,                                &
                  theta, q, qcl, qcf, exner, rho, tracer, srce,                &
                  level, timestep, i_hour, i_minute, amp)
IMPLICIT NONE
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: offx
INTEGER, INTENT(IN) :: offy
INTEGER, INTENT(IN) :: halo_i
INTEGER, INTENT(IN) :: halo_j
REAL, INTENT(IN) :: r_theta_levels (:,:,0:)
REAL, INTENT(IN) :: r_rho_levels   (:,:,:)
REAL, INTENT(IN) :: theta(:,:,:)
REAL, INTENT(IN) :: q(:,:,:)
REAL, INTENT(IN) :: qcl(:,:,:)
REAL, INTENT(IN) :: qcf(:,:,:)
REAL, INTENT(IN) :: exner(:,:,:)
REAL, INTENT(IN) :: rho(:,:,:)
REAL, INTENT(IN OUT) :: tracer(:,:)
REAL, INTENT(IN) :: srce(:,:)
INTEGER, INTENT(IN) :: level
REAL, INTENT(IN) :: timestep
INTEGER, INTENT(IN) :: i_hour
INTEGER, INTENT(IN) :: i_minute
REAL, INTENT(IN) :: amp
END SUBROUTINE trsrce

! Subroutine for JULES-based dry deposition required by ukca_chemistry_ctl,
! ukca_chemistry_ctl_BE_mod & ukca_chemistry_ctl_col_mod

SUBROUTINE deposition_from_ukca_chemistry(                                     &
             secs_per_step, bl_levels, row_length, rows, ntype, npft,          &
             jpspec_ukca, ndepd_ukca, nldepd_ukca, speci_ukca,                 &
             land_points, land_index, tile_pts, tile_index,                    &
             seaice_frac, fland, sinlat,                                       &
             p_surf, rh, t_surf, surf_hf, surf_wetness,                        &
             z0tile_lp, stcon, laift_lp, canhtft_lp, t0tile_lp,                &
             soilmc_lp, zbl, dzl, frac_types, u_s,                             &
             dep_loss_rate_ij, nlev_with_ddep, len_stashwork50, stashwork50)

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
  row_length, rows, bl_levels,                                                 &
   ! size of UKCA x, y dimensions and no. of atmospheric boundary layer levels
  ntype, npft,                                                                 &
   ! no. of surface types and no. of plant functional types
  jpspec_ukca,                                                                 &
   ! no. of chemical species in the UKCA mechanism
  ndepd_ukca,                                                                  &
   ! no. of chemical species that are deposited
   ! Equal to and used interchangeably with jpdd
  land_points

INTEGER, INTENT(IN) ::                                                         &
  land_index(land_points),                                                     &
  tile_pts(ntype),                                                             &
  tile_index(land_points,ntype),                                               &
  nldepd_ukca(jpspec_ukca)
    ! Holds array elements of speci, identifying those chemical species
    ! in the chemical mechanism that are deposited

REAL, INTENT(IN) ::                                                            &
  secs_per_step,                                                               &
   ! time step (in s)
  sinlat(row_length, rows),                                                    &
   ! sin(latitude)
  p_surf(row_length, rows),                                                    &
   ! surface pressure
  dzl(row_length, rows, bl_levels),                                            &
   ! separation of boundary-layer levels
  u_s(row_length, rows),                                                       &
   ! surface friction velocity (m s-1)
  t_surf(row_length, rows),                                                    &
   ! surface temperature
  rh(row_length, rows),                                                        &
   ! relative humidity (-)
  surf_wetness(row_length, rows),                                              &
   ! surface wetness
  frac_types(land_points,ntype),                                               &
   ! surface tile fractions (-)
  zbl(row_length,rows),                                                        &
   ! boundary-layer height (m)
  surf_hf(row_length,rows),                                                    &
   ! sensible heat flux (W m-2)
  seaice_frac(row_length,rows),                                                &
   ! grid-cell sea ice fraction
  stcon(row_length,rows,npft),                                                 &
   ! stomatal conductance (s m-1)
  soilmc_lp(land_points),                                                      &
   ! soil moisture
  fland(land_points),                                                          &
   ! grid-cell land fraction
  laift_lp(land_points,npft),                                                  &
   ! leaf area index (m2 m-2)
  canhtft_lp(land_points,npft),                                                &
   ! canopy height (m)
  z0tile_lp(land_points,ntype),                                                &
   ! roughness length for heat and moisture (m) by surface type
  t0tile_lp(land_points,ntype)
   ! surface temperature (K) by surface type

CHARACTER(LEN=10), INTENT(IN) :: speci_ukca(jpspec_ukca)
                                   ! Names of all the chemical species
                                   ! in the UKCA chemical mechanism

! IN JULES-based code, INTENT(IN OUT) but make INRTENT(IN) here
INTEGER, INTENT(IN)  :: len_stashwork50
REAL, INTENT(IN)     :: stashwork50 (len_stashwork50)
                          ! Diagnostics array (UKCA stashwork50)

! Output variables
INTEGER, INTENT(OUT) ::                                                        &
  nlev_with_ddep(row_length, rows)
    ! Number of levs with deposition in boundary layer

REAL, INTENT(OUT) ::                                                           &
  dep_loss_rate_ij(row_length, rows, ndepd_ukca)
    ! dry deposition loss rate (s-1)

! Initialise output variables
nlev_with_ddep(:,:) = 0
dep_loss_rate_ij(:,:,:) = 0.0

RETURN

END SUBROUTINE deposition_from_ukca_chemistry

! Ozone column calculation subroutine required by ukca_chem_diags_allts_mod

SUBROUTINE ukca_calc_ozonecol(row_length, rows, model_levels, z_top_of_model,  &
                              p_layer_boundaries, p_layer_centres, ozone_vmr,  &
                              ozonecol)
IMPLICIT NONE
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
REAL, INTENT(IN) :: z_top_of_model
REAL, INTENT(IN) :: p_layer_boundaries(:,:,:)
REAL, INTENT(IN) :: p_layer_centres(:,:,:)
REAL, INTENT(IN) :: ozone_vmr(:,:,:)
REAL, INTENT(OUT) :: ozonecol(:,:,:)
ozonecol = 0.0
END SUBROUTINE ukca_calc_ozonecol

! Particle surface area calculation subroutine required by ukca_chemco_raq_mod

SUBROUTINE calc_surf_area (n_pnts, rho_air, rh_frac,                           &
                           so4_aitken,    so4_accum,                           &
                           soot_fresh,    soot_aged,                           &
                           ocff_fresh,    ocff_aged, biogenic,                 &
                           sea_salt_film, sea_salt_jet,                        &
                           sa_so4_ait,    sa_so4_acc,                          &
                           sa_bc_fresh,   sa_bc_aged,                          &
                           sa_ocff_fresh, sa_ocff_aged, sa_soa,                &
                           sa_ss_film,    sa_ss_jet,                           &
                           wr_so4_ait,    wr_so4_acc,                          &
                           wr_bc_fresh,   wr_bc_aged,                          &
                           wr_ocff_fresh, wr_ocff_aged, wr_soa,                &
                           wr_ss_film,    wr_ss_jet)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n_pnts
REAL, INTENT(IN) :: rho_air(:)
REAL, INTENT(IN) :: rh_frac(:)
REAL, INTENT(IN) :: so4_aitken(:)
REAL, INTENT(IN) :: so4_accum(:)
REAL, INTENT(IN) :: soot_fresh(:)
REAL, INTENT(IN) :: soot_aged(:)
REAL, INTENT(IN) :: ocff_fresh(:)
REAL, INTENT(IN) :: ocff_aged(:)
REAL, INTENT(IN) :: biogenic(:)
REAL, INTENT(IN) :: sea_salt_film(:)
REAL, INTENT(IN) :: sea_salt_jet(:)
REAL, INTENT(OUT) :: sa_so4_ait(:)
REAL, INTENT(OUT) :: sa_so4_acc(:)
REAL, INTENT(OUT) :: sa_bc_fresh(:)
REAL, INTENT(OUT) :: sa_bc_aged(:)
REAL, INTENT(OUT) :: sa_ocff_fresh(:)
REAL, INTENT(OUT) :: sa_ocff_aged(:)
REAL, INTENT(OUT) :: sa_soa(:)
REAL, INTENT(OUT) :: sa_ss_film(:)
REAL, INTENT(OUT) :: sa_ss_jet(:)
REAL, INTENT(OUT) :: wr_so4_ait(:)
REAL, INTENT(OUT) :: wr_so4_acc(:)
REAL, INTENT(OUT) :: wr_bc_fresh(:)
REAL, INTENT(OUT) :: wr_bc_aged(:)
REAL, INTENT(OUT) :: wr_ocff_fresh(:)
REAL, INTENT(OUT) :: wr_ocff_aged(:)
REAL, INTENT(OUT) :: wr_soa(:)
REAL, INTENT(OUT) :: wr_ss_film(:)
REAL, INTENT(OUT) :: wr_ss_jet(:)
sa_so4_ait = 0.0
sa_so4_acc = 0.0
sa_bc_fresh = 0.0
sa_bc_aged = 0.0
sa_ocff_fresh = 0.0
sa_ocff_aged = 0.0
sa_soa = 0.0
sa_ss_film = 0.0
sa_ss_jet = 0.0
wr_so4_ait = 0.0
wr_so4_acc = 0.0
wr_bc_fresh = 0.0
wr_bc_aged = 0.0
wr_ocff_fresh = 0.0
wr_ocff_aged = 0.0
wr_soa = 0.0
wr_ss_film = 0.0
wr_ss_jet = 0.0
END SUBROUTINE calc_surf_area

! SO2 volcanic emissions subroutine required by ukca_add_emiss_mod

SUBROUTINE ukca_volcanic_so2 (so2_mmr, mass, row_length, rows, model_levels,   &
           year, timestep, r_theta_levels)
IMPLICIT NONE
REAL, INTENT(IN OUT) :: so2_mmr(:,:,:)
REAL, INTENT(IN) :: mass(:,:,:)
INTEGER, INTENT(IN) :: row_length, rows, model_levels
INTEGER, INTENT(IN) :: year
REAL, INTENT(IN) :: timestep
REAL, INTENT(IN) :: r_theta_levels(:,:,:)
so2_mmr = 0.0
END SUBROUTINE ukca_volcanic_so2

! Diagnostics utility function required by ukca_emiss_diags_mod and
! ukca_update_emdiagstruct_mod

FUNCTION get_emdiag_stash(diag_name)
IMPLICIT NONE
CHARACTER (LEN=*), INTENT(IN) :: diag_name
INTEGER :: get_emdiag_stash
get_emdiag_stash = 0
END FUNCTION get_emdiag_stash

! Diagnostics utility function required by ukca_ddcalc_mod

SUBROUTINE set_pseudo_list(n_levels, len_stlist, stlist, pseudo_list,          &
                           stash_pseudo_levels, num_stash_pseudo)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n_levels
INTEGER, INTENT(IN) :: len_stlist
INTEGER, INTENT(IN) :: stlist(:)
LOGICAL, INTENT(OUT) :: pseudo_list(:)
INTEGER, INTENT(IN) :: stash_pseudo_levels(:,:)
INTEGER, INTENT(IN) :: num_stash_pseudo
pseudo_list = .FALSE.
END SUBROUTINE set_pseudo_list

! Routines for copying diagnostics, required by various modules

SUBROUTINE copydiag_real(diagout, diagin, row_length, rows)
IMPLICIT NONE
REAL, INTENT(OUT) :: diagout(:)
REAL, INTENT(IN)  :: diagin(:,:)
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
diagout = 0.0
END SUBROUTINE copydiag_real

SUBROUTINE copydiag_real_vector(diagout, diagin, length)
IMPLICIT NONE
REAL, INTENT(OUT) :: diagout(:)
REAL, INTENT(IN)  :: diagin(:)
INTEGER, INTENT(IN) :: length
diagout = 0.0
END SUBROUTINE copydiag_real_vector

SUBROUTINE copydiag_integer(diagout, diagin, row_length, rows)
IMPLICIT NONE
REAL, INTENT(OUT) :: diagout(:)
INTEGER, INTENT(IN) :: diagin(:,:)
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
diagout = 0.0
END SUBROUTINE copydiag_integer

SUBROUTINE copydiag_logical(diagout, diagin, row_length, rows)
IMPLICIT NONE
REAL, INTENT(OUT) :: diagout(:)
LOGICAL, INTENT(IN) :: diagin(:,:)
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
diagout = 0.0
END SUBROUTINE copydiag_logical

SUBROUTINE copydiag_3d_real(diagout, diagin, row_length, rows, levels,         &
                            stlist, len_stlist, stash_levels, len_stashlevels)
IMPLICIT NONE
REAL, INTENT(OUT) :: diagout(:)
REAL, INTENT(IN) :: diagin(:,:,:)
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: levels
INTEGER, INTENT(IN) :: stlist(:)
INTEGER, INTENT(IN) :: len_stlist
INTEGER, INTENT(IN) :: stash_levels(:,:)
INTEGER, INTENT(IN) :: len_stashlevels
diagout = 0.0
END SUBROUTINE copydiag_3d_real

SUBROUTINE copydiag_3d_integer(diagout, diagin, row_length, rows, levels,      &
                               stlist, len_stlist, stash_levels,               &
                               len_stashlevels)
IMPLICIT NONE
REAL, INTENT(OUT) :: diagout(:)
INTEGER, INTENT(IN) :: diagin(:,:,:)
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: levels
INTEGER, INTENT(IN) :: stlist(:)
INTEGER, INTENT(IN) :: len_stlist
INTEGER, INTENT(IN) :: stash_levels(:,:)
INTEGER, INTENT(IN) :: len_stashlevels
diagout = 0.0
END SUBROUTINE copydiag_3d_integer

SUBROUTINE copydiag_3d_logical(diagout, diagin, row_length, rows, levels,      &
                               stlist, len_stlist, stash_levels,               &
                               len_stashlevels)
IMPLICIT NONE
REAL, INTENT(OUT) :: diagout(:)
LOGICAL, INTENT(IN) :: diagin(:,:,:)
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: levels
INTEGER, INTENT(IN) :: stlist(:)
INTEGER, INTENT(IN) :: len_stlist
INTEGER, INTENT(IN) :: stash_levels(:,:)
INTEGER, INTENT(IN) :: len_stashlevels
diagout = 0.0
END SUBROUTINE copydiag_3d_logical

! Autotune subroutines required by ukca_main1_mod and ukca_chemistry_ctl_col_mod

SUBROUTINE autotune_entry(this, segment_size, first)
IMPLICIT NONE
TYPE(autotune_type), INTENT(IN OUT) :: this
INTEGER, INTENT(OUT) :: segment_size
LOGICAL, OPTIONAL, INTENT(IN) :: first
segment_size = 0
END SUBROUTINE autotune_entry

SUBROUTINE autotune_init(this, region_name, tag, start_size)
IMPLICIT NONE
TYPE(autotune_type), INTENT(OUT) :: this
CHARACTER(LEN=*),    INTENT(IN)  :: region_name
CHARACTER(LEN=*),    INTENT(IN)  :: tag
INTEGER,             INTENT(IN)  :: start_size
this%dummy_value = 0
END SUBROUTINE autotune_init

SUBROUTINE autotune_return(this, last)
IMPLICIT NONE
TYPE(autotune_type), INTENT(IN OUT) :: this
LOGICAL, OPTIONAL, INTENT(IN) :: last
END SUBROUTINE autotune_return

SUBROUTINE autotune_start_region(this, local_points)
IMPLICIT NONE
TYPE(autotune_type), INTENT(IN OUT) :: this
INTEGER, INTENT(IN) :: local_points
END SUBROUTINE autotune_start_region

SUBROUTINE autotune_stop_region(this)
IMPLICIT NONE
TYPE(autotune_type), INTENT(IN OUT) :: this
END SUBROUTINE autotune_stop_region

! Timer subroutine required by ukca_main1_mod

SUBROUTINE timer(sub, i_arg, starttime)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: sub
INTEGER, INTENT(IN)          :: i_arg
REAL, OPTIONAL, INTENT(IN)   :: starttime
END SUBROUTINE timer

END MODULE ukca_um_legacy_mod
