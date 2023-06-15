! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Top level module for the UKCA emission system, called from UKCA_MAIN1
!     at each timestep.
!
!  Method:
!  1) Call ukca_onl_emiss_init to check for and set up the 'online' emissions
!  2) Copy attributes and data for NetCDF/ Offline emissions provided by the
!     parent
!  3) Calculate online emissions (lightning NOx, DMS-flux, aerosols)
!  4) Remap JULES emissions (passed by parent) to 2-D lat/long grid
!  5) Call ukca_add_emiss to inject emissions and do tracer mixing
!  6) Call ukca_emiss_diags to produce emission diagnostics
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! --------------------------------------------------------------------------
!
MODULE ukca_emiss_ctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EMISS_CTL_MOD'

CONTAINS

SUBROUTINE ukca_emiss_ctl (                                                    &
  row_length, rows, bl_levels, model_levels, n_tracers, ndiv,                  &
  latitude, longitude, sin_latitude, cos_latitude, tan_latitude,               &
  iyear, imonth, iday, ihour,                                                  &
  timestep, l_first, land_points, land_index,                                  &
  conv_cloud_base, conv_cloud_top,                                             &
  delta_lambda, delta_phi, r_theta_levels, surf_area,                          &
  cos_zenith_angle, int_zenith_angle, land_fraction, tropopause_level,         &
  r_rho_levels, t_theta_levels,                                                &
  p_theta_levels, p_layer_boundaries, rel_humid_frac_clr, mass,                &
  ls_mask,  theta, q, qcl, qcf,                                                &
  exner_rho_levels, rho_r2, kent, kent_dsc, rhokh_rdz, dtrdz,                  &
  we_lim, t_frac, zrzi, we_lim_dsc, t_frac_dsc, zrzi_dsc, ml_depth,            &
  zhsc, z_half, ch4_wetl_emiss, seaice_frac, area, dust_flux, u_scalar_10m,    &
  tstar, dms_sea_conc, chloro_sea,                                             &
  dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6, tracers,   &
  ext_cg_flash, ext_ic_flash,                                                  &
  len_stashwork38, stashwork38, len_stashwork50, stashwork50)

USE ukca_emiss_mod,       ONLY: num_em_flds, num_cdf_em_flds, emissions,       &
                                ukca_onl_emiss_init,                           &
                                ukca_emiss_spatial_vars_init,                  &
                                inox_light, ich4_wetl,                         &
                                ic5h8_ibvoc, ic10h16_ibvoc,                    &
                                ich3oh_ibvoc, ich3coch3_ibvoc,                 &
                                ico_inferno, ich4_inferno,                     &
                                inox_inferno, iso2_inferno,                    &
                                ioc_inferno, ibc_inferno,                      &
                                ukca_emiss_update_mode, ukca_emiss_mode_map,   &
                                idms_seaflux, iseasalt_first,                  &
                                ipmoc_first, idust_first,                      &
                                ino3_first, inh4_first, iseasalt_hno3,         &
                                inano3_cond, idust_hno3
USE ukca_emiss_struct_mod, ONLY: ncdf_emissions, ukca_copy_emiss_struct
USE ukca_emiss_mode_mod,  ONLY: aero_ems_species

USE ukca_emiss_factors,   ONLY: vertical_emiss_factors
USE ukca_add_emiss_mod,   ONLY: ukca_add_emiss
USE ukca_um_legacy_mod,   ONLY: sf, stashcode_glomap_sec, ukca_diag_sect,      &
                                rgas => r, dms_flux_4a
USE ukca_config_defs_mod, ONLY: n_chem_emissions, n_3d_emissions, em_chem_spec,&
                                n_chem_tracers, n_aero_tracers, n_mode_tracers
USE ukca_emiss_diags_mod, ONLY: ukca_emiss_diags
USE ukca_emiss_diags_mode_mod, ONLY: ukca_emiss_diags_mode

USE ukca_mode_setup,      ONLY: nmodes, cp_cl, cp_du,                          &
                                cp_no3, cp_nh4, cp_nn, cp_su

USE ukca_prod_no3_mod,    ONLY: ukca_prod_no3_fine, ukca_prod_no3_coarse,      &
                                ukca_no3_check_values
USE ukca_calc_rho_mod,    ONLY: ukca_calc_rho
USE ukca_cspecies,        ONLY: n_nh3, n_hono2
USE ukca_mode_verbose_mod, ONLY: verbose => glob_verbose
USE ukca_prim_du_mod,     ONLY: ukca_prim_du
USE ukca_environment_fields_mod, ONLY: ibvoc_isoprene,  ibvoc_terpene,         &
                                 ibvoc_methanol, ibvoc_acetone,                &
                                 inferno_bc, inferno_ch4, inferno_co,          &
                                 inferno_nox, inferno_oc, inferno_so2
USE ukca_constants,       ONLY: m_c, m_ch4, m_dms,                             &
                                m_n, m_no, m_no2, m_s, m_c5h8, m_monoterp,     &
                                m_ch3oh, m_me2co
USE chemistry_constants_mod, ONLY: boltzmann, avogadro

USE ukca_config_specification_mod, ONLY: ukca_config, glomap_config,           &
                                         i_light_param_off, glomap_variables,  &
                                         i_solinsol_6mode

USE asad_mod,             ONLY: advt, jpctr
USE asad_chem_flux_diags, ONLY: L_asad_use_chem_diags, L_asad_use_light_ems,   &
                                asad_3d_emissions_diagnostics,                 &
                                lightning_emissions
USE ereport_mod,          ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1,             ONLY: jpim, jprb      ! DrHook
USE yomhook,              ONLY: lhook, dr_hook  ! DrHook

USE ukca_diurnal_isop_ems_mod, ONLY: ukca_diurnal_isop_ems
USE ukca_light_ctl_mod, ONLY: ukca_light_ctl
USE ukca_prim_ss_mod, ONLY: ukca_prim_ss
USE ukca_prim_moc_mod, ONLY: ukca_prim_moc

IMPLICIT NONE

! Subroutine arguments

! Iput arguments with info on model dimensions
INTEGER,          INTENT(IN)    :: row_length
INTEGER,          INTENT(IN)    :: rows
INTEGER,          INTENT(IN)    :: bl_levels
INTEGER,          INTENT(IN)    :: model_levels

! Input arguments to get nr tracers, nr dust divisions & current model time
INTEGER, INTENT(IN) :: n_tracers       ! nr tracers
INTEGER, INTENT(IN) :: ndiv            ! nr dust divisions
INTEGER, INTENT(IN) :: iyear           ! current yr, mon, day, hr
INTEGER, INTENT(IN) :: imonth
INTEGER, INTENT(IN) :: iday
INTEGER, INTENT(IN) :: ihour

REAL,    INTENT(IN) :: timestep        ! timestep length (sec)
LOGICAL, INTENT(IN) :: l_first         ! T if first call to UKCA

! Input arguments used if iBVOC emissions from JULES
INTEGER, INTENT(IN) :: land_points
INTEGER, INTENT(IN) :: land_index(land_points)

! Input arguments to apply diurnal cycle to isoprene emissions,
! get SO2 emiss from explosive volcanic eruptions, calculate
! lightning emissions of NOx, get vertical profiles or calculate
! surf_area
INTEGER, INTENT(IN) :: conv_cloud_base    (1:row_length, 1:rows)
INTEGER, INTENT(IN) :: conv_cloud_top     (1:row_length, 1:rows)

REAL, INTENT(IN) :: latitude              (1:row_length, 1:rows)  ! degrees N
REAL, INTENT(IN) :: longitude             (1:row_length, 1:rows)  ! degrees E
REAL, INTENT(IN) :: sin_latitude          (1:row_length, 1:rows)
REAL, INTENT(IN) :: cos_latitude          (1:row_length, 1:rows)
REAL, INTENT(IN) :: tan_latitude          (1:row_length, 1:rows)
REAL, INTENT(IN) :: delta_lambda
REAL, INTENT(IN) :: delta_phi
REAL, INTENT(IN) :: r_theta_levels        (1:row_length, 1:rows, 0:model_levels)
REAL, INTENT(IN) :: surf_area             (1:row_length, 1:rows)
REAL, INTENT(IN) :: cos_zenith_angle      (1:row_length, 1:rows)
REAL, INTENT(IN) :: int_zenith_angle      (1:row_length, 1:rows)
REAL, INTENT(IN) :: land_fraction         (1:row_length, 1:rows)
REAL, INTENT(IN) :: r_rho_levels          (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: t_theta_levels        (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: p_theta_levels        (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: p_layer_boundaries    (1:row_length, 1:rows, 0:model_levels)
REAL, INTENT(IN) :: mass                  (1:row_length, 1:rows, 1:model_levels)
INTEGER, INTENT(IN) :: tropopause_level   (1:row_length, 1:rows)
REAL, INTENT(IN) :: rel_humid_frac_clr    (1:row_length, 1:rows, 1:model_levels)

LOGICAL, INTENT(IN) :: ls_mask            (1:row_length, 1:rows)

! Input arguments needed by UKCA_ADD_EMISS to call TRSRCE
REAL, INTENT(IN) :: theta             (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: q                 (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: qcl               (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: qcf               (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: exner_rho_levels  (1:row_length, 1:rows, 1:model_levels+1)
REAL, INTENT(IN) :: rho_r2            (1:row_length, 1:rows, 1:model_levels)

! Input arguments needed by UKCA_ADD_EMISS to call TR_MIX
INTEGER, INTENT(IN) :: kent       (1:row_length, 1:rows)
INTEGER, INTENT(IN) :: kent_dsc   (1:row_length, 1:rows)

! Input arguments needed by UKCA_ADD_EMISS to call UKCA_PRIM_SS, UKCA_PRIM_MOC,
! UKCA_PRIM_DU and DMS_FLUX_4A.
REAL, INTENT(IN) :: seaice_frac    (1:row_length, 1:rows)
REAL, INTENT(IN) :: dust_flux(row_length,rows,ndiv) ! dust emiss (kg m-2 s-1)
REAL, INTENT(IN) :: u_scalar_10m   (:,:)
REAL, INTENT(IN) :: tstar          (1:row_length, 1:rows)
REAL, INTENT(IN) :: dms_sea_conc   (1:row_length, 1:rows)
REAL, INTENT(IN) :: chloro_sea     (1:row_length, 1:rows) ! surface chl-a kg/m3

! Area of grid cell, required by ukca_emiss_diags_mode
REAL, INTENT(IN) :: area(row_length, rows, model_levels)

REAL, INTENT(IN) :: rhokh_rdz  (1:row_length, 1:rows, 2:bl_levels)
REAL, INTENT(IN) :: dtrdz      (1:row_length, 1:rows, 1:bl_levels)

REAL, INTENT(IN) :: we_lim     (1:row_length, 1:rows,                          &
                                1:ukca_config%nlev_ent_tr_mix)
REAL, INTENT(IN) :: t_frac     (1:row_length, 1:rows,                          &
                                1:ukca_config%nlev_ent_tr_mix)
REAL, INTENT(IN) :: zrzi       (1:row_length, 1:rows,                          &
                                1:ukca_config%nlev_ent_tr_mix)

REAL, INTENT(IN) :: we_lim_dsc (1:row_length, 1:rows,                          &
                                1:ukca_config%nlev_ent_tr_mix)
REAL, INTENT(IN) :: t_frac_dsc (1:row_length, 1:rows,                          &
                                1:ukca_config%nlev_ent_tr_mix)
REAL, INTENT(IN) :: zrzi_dsc   (1:row_length, 1:rows,                          &
                                1:ukca_config%nlev_ent_tr_mix)

REAL, INTENT(IN) :: ml_depth   (1:row_length, 1:rows)
REAL, INTENT(IN) :: zhsc       (1:row_length, 1:rows)
REAL, INTENT(IN) :: z_half     (1:row_length, 1:rows,1:bl_levels)

! CLASSIC dust mass mixing ratios - will become obsolete
! when UKCA dust supercedes CLASSIC dust
REAL, INTENT(IN) :: dust_div1 (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: dust_div2 (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: dust_div3 (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: dust_div4 (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: dust_div5 (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: dust_div6 (1:row_length, 1:rows, 1:model_levels)

! External cloud-to-ground and intracloud lightning flashes
REAL, INTENT(IN) :: ext_cg_flash (1:row_length, 1:rows)
REAL, INTENT(IN) :: ext_ic_flash (1:row_length, 1:rows)

! Length of diagnostics arrays
INTEGER, INTENT(IN) :: len_stashwork38, len_stashwork50

! Wetland emissions of methane - scaled by land_fraction while adding
! to emissions structure. Any negative values are removed.
REAL, INTENT(IN OUT) :: ch4_wetl_emiss (1:row_length,1:rows)

! Tracer mass mixing ratios
REAL, INTENT(IN OUT) :: tracers (1:row_length, 1:rows, 1:model_levels,         &
                                1:n_tracers)
! Diagnostics arrays
REAL, INTENT(IN OUT) :: stashwork38(len_stashwork38)
REAL, INTENT(IN OUT) :: stashwork50(len_stashwork50)

! Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: mm (:)


INTEGER                 :: i, j, k, l, ilev
INTEGER                 :: section           ! stash section
REAL                    :: base_scaling      ! factor to convert emiss units
REAL                    :: vert_fact_3d (1:row_length, 1:rows, 1:model_levels)
                                     ! 3-D vertical profiles of emissions
REAL                    :: aer_emmas (row_length, rows, model_levels, nmodes)
                                     ! aerosol mass emissions by mode
REAL                    :: aer_emnum (row_length, rows, model_levels, nmodes)
                                     ! aerosol number emissions by mode
REAL                    :: chla(row_length, rows)
                                     ! chlorophyll-a, mg.m-3
REAL                    :: mass_pmoc (row_length, rows, 1)
                                     ! primary marine organic carbon emission
INTEGER                 :: emiss_levs     ! num levs for emissions
INTEGER                 :: icp            ! component
LOGICAL                 :: lmode_emiss(nmodes)   ! modes which are emitted into
                                                 ! for a given emission field
REAL                    :: f_dms_sea(row_length, rows)  ! DMS flux over sea

! Nitrate arrays
REAL, ALLOCATABLE :: aer_emmas_nh4(:,:,:,:)
REAL, ALLOCATABLE :: aer_emnum_nh4(:,:,:,:)
REAL, ALLOCATABLE :: aer_emmas_dust(:,:,:,:)
REAL, ALLOCATABLE :: aer_emmas_nacl(:,:,:,:)
REAL, ALLOCATABLE :: air_density(:,:,:) ! Air density (kg m-3)
REAL, ALLOCATABLE :: air_burden(:,:,:) ! Air burden (kg m-3)

! INFERNO emission profile weight
REAL :: inferno_vert_fact (1:model_levels)

! 2-D field to hold biogenic isoprene emissions, with base scaling
! applied (in case units are not kg m-2 s-1), but without diurnal cycle
REAL, ALLOCATABLE, SAVE :: biogenic_isop (:,:)
REAL, ALLOCATABLE :: tmp_in_em_field  (:,:)  ! Climatological and daily varying
REAL :: tmp_out_em_field (row_length, rows)  ! isoprene emission field
LOGICAL, SAVE    :: testdcycl = .FALSE.      ! True only for debugging
                                             ! UKCA_DIURNAL_ISOP_EMS

LOGICAL  :: l_3dsource = .FALSE.             ! source emissions are 3D

INTEGER :: ils_mask (row_length, rows)  ! Land/sea mask (1/0) used
                                        ! for lightning NOx

! Arrays of NOx lightning emiss with different units:
! kg(N)/grid box/s and kg(NO2)/kg(air)/s
REAL :: lightningem_n_gridbox  (1:row_length, 1:rows, 1:model_levels)
REAL :: lightningem_no2_to_air (1:row_length, 1:rows, 1:model_levels)

!  Molar mass of dry air (kg/mol) =avogadro*boltzmann/rgas
REAL    :: mm_da

! Aerosol mode setup 11 (SOL/INSOL) has a single aerosol component
! represented by the H2SO4 index. Other aerosol species (OC, BC, SS)
! are emitted into cp_so4. The variables this_cp_* are used to determine
! the active component which this component is emitted into
INTEGER :: this_cp_cl

INTEGER                            :: errcode    ! Error code for ereport
CHARACTER (LEN=errormessagelength) :: cmessage   ! Error message
INTEGER            :: ierr       ! Req. for ASAD diags (do not initialise here)
INTEGER, SAVE      :: inox       ! Index of NOx tracer for ASAD 3D diags

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_CTL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mm          => glomap_variables%mm

! Initialisation
errcode  = 0

inferno_vert_fact (:) = 0.0  ! INFERNO emission profile

IF ( l_first ) THEN

  ! Determine the number of 'online' (model-generated + coupling) emission
  ! fields and set up the data structures. Emissions from NetCDF files are read
  ! beforehand
  CALL ukca_onl_emiss_init (row_length, rows, model_levels)

  ! Copy the file emissions data into the combined structure, to avoid having to
  ! loop over both arrays twice. This needs to be repeated at each timestep.
  ! Check that each registered emission also has the 'values' array set, which
  ! should have been done by the parent before calling UKCA.
  IF (num_cdf_em_flds > 0) THEN
    DO k = 1, num_cdf_em_flds
      IF ( .NOT. ASSOCIATED(ncdf_emissions(k)%values) ) THEN
        errcode  = k
        cmessage = 'Emission field has been registered but no values set: '//  &
                   TRIM(ncdf_emissions(k)%var_name)
        CALL ereport (ModuleName//':'//RoutineName, errcode, cmessage)
      END IF
      ! Explicitly copy the required components
      CALL ukca_copy_emiss_struct(ncdf_emissions(k),emissions(k))
    END DO
  END IF

  ! Make sure that there are no missing fields in the emission structure,
  ! i.e. that all names in em_chem_spec are present in emissions structure
  ! An exception is that organic matter/BC emissions can be missing if
  ! l_ukca_primbcoc is false. The same exception doesn't apply to SO2 with
  ! l_ukca_primsu, since SO2 is emitted into the chemistry.
  DO k = 1, n_chem_emissions + n_3d_emissions
    IF (.NOT. ANY(emissions(:)%tracer_name == em_chem_spec (k)) ) THEN
      IF (.NOT. (                                                              &
          ((em_chem_spec(k)(1:2) == 'BC') .AND.                                &
           .NOT. glomap_config%l_ukca_primbcoc) .OR.                           &
          ((em_chem_spec(k)(1:2) == 'OM') .AND.                                &
           .NOT. glomap_config%l_ukca_primbcoc))) THEN
        cmessage = TRIM(em_chem_spec(k))//' missing from supplied emissions '
        errcode = k
        CALL ereport (ModuleName//':'//RoutineName, errcode, cmessage)
      END IF
    END IF
  END DO

  ! Get index for NOx tracer (needed only if ASAD diagnostics
  ! for NOx lightning emiss are required)
  DO k = 1, jpctr
    SELECT CASE (advt(k))
    CASE ('NOx       ')
      inox = k
    CASE ('NO        ')
      inox = k
    END SELECT
  END DO

  IF (inox == -99 .AND.                                                        &
      ukca_config%i_ukca_light_param /= i_light_param_off) THEN
    errcode  = 1
    cmessage = 'Did not find NO or NOx tracer'
    CALL ereport (ModuleName//':'//RoutineName, errcode, cmessage)
  END IF

  ! Calculate scaling factors indicating how to to spread emissions over
  ! different vertical levels - except for file based emissions, where the
  ! array can be received from the parent (vert_scaling_3d is
  ! already allocated)
  DO l = 1, num_em_flds
    IF ( .NOT. ASSOCIATED(emissions(l)%vert_scaling_3d) ) THEN
      ALLOCATE (emissions(l)%vert_scaling_3d (row_length, rows, model_levels))

      CALL vertical_emiss_factors (row_length, rows, model_levels,             &
             emissions(l)%lowest_lev, emissions(l)%highest_lev,                &
             emissions(l)%vert_fact, emissions(l)%var_name, vert_fact_3d )

      emissions(l)%vert_scaling_3d (:,:,:) = vert_fact_3d (:,:,:)
    END IF
  END DO

ELSE

  ! If persistence is off then spatial arrays in the emissions
  ! structure need to be reinitialised each timestep
  IF (ukca_config%l_ukca_persist_off) THEN

    ! Analogous to ukca_onl_emiss_init()
    CALL ukca_emiss_spatial_vars_init (row_length, rows, model_levels)

    ! Copy file emissions data into the combined structure, to avoid having
    ! to loop over both arrays twice. This needs to be repeated each timestep
    ! Check that each registered emission also has the 'values' array, which
    ! should have been done by the parent before calling UKCA.
    IF (num_cdf_em_flds > 0) THEN
      DO k = 1, num_cdf_em_flds
        IF ( .NOT. ASSOCIATED(ncdf_emissions(k)%values) ) THEN
          errcode  = k
          cmessage = 'Emission field has been registered but no values set: '//&
                     TRIM(ncdf_emissions(k)%var_name)
          CALL ereport (ModuleName//':'//RoutineName, errcode, cmessage)
        END IF
        ! Explicitly copy the required components
        CALL ukca_copy_emiss_struct(ncdf_emissions(k),emissions(k))
      END DO
    END IF

    ! Calculate scaling factors indicating how to to spread emissions over
    ! different vertical levels - except for file based emissions, where the
    ! array can be received from the parent (vert_scaling_3d is
    ! already allocated)
    DO l = 1, num_em_flds
      IF ( .NOT. ASSOCIATED(emissions(l)%vert_scaling_3d) ) THEN
        ALLOCATE (emissions(l)%vert_scaling_3d(row_length, rows, model_levels))

        CALL vertical_emiss_factors (row_length, rows, model_levels,           &
               emissions(l)%lowest_lev, emissions(l)%highest_lev,              &
               emissions(l)%vert_fact, emissions(l)%var_name, vert_fact_3d )

        emissions(l)%vert_scaling_3d (:,:,:) = vert_fact_3d (:,:,:)
      END IF
    END DO

  ELSE

    ! For other timesteps, only need to copy values (occassionally updated from
    ! file), vertical scaling array and update flag from ncdf_emissions part
    DO k = 1, num_cdf_em_flds
      IF ( .NOT. emissions(k)%l_online ) THEN
        IF ( ncdf_emissions(k)%l_update ) THEN
          emissions(k)%values(:,:,:) = ncdf_emissions(k)%values(:,:,:)
          emissions(k)%l_update = ncdf_emissions(k)%l_update
          emissions(k)%vert_scaling_3d(:,:,:) =                                &
                                       ncdf_emissions(k)%vert_scaling_3d(:,:,:)
        ELSE   ! Reset internal flag to avoid double processing
          emissions(k)%l_update = .FALSE.
        END IF
      END IF
    END DO

  END IF ! ukca_config%l_ukca_persist_off

END IF  ! IF l_first

! If the emissions structure contains an offline biogenic isoprene field
! that needs a diurnal cycle to be applied online (via the routine
! UKCA_DIURNAL_ISOP_EMS) then create a 2-D field which will hold
! such emissions before applying the diurnal cycle.
! That field will be updated only when emission values are
! updated, depending on the updating frequency, but will not
! include the diurnal correction to avoid re-applying it

! Also avoid double-counting emissions: Do not allow the use of
! iBVOC emissions for isoprene (from JULES) when an offline
! biogenic emission field for the same species has been read
! from the NetCDF files.

IF (num_cdf_em_flds > 0) THEN

  IF (ANY (ncdf_emissions(:)%tracer_name == 'C5H8      '  .AND.                &
    ncdf_emissions(:)%hourly_fact == 'diurnal_isopems') .AND.                  &
    ukca_config%l_diurnal_isopems .AND. .NOT. ukca_config%l_ukca_ibvoc .AND.   &
    ( .NOT. ALLOCATED(biogenic_isop) ) .AND.                                   &
    ( .NOT. ukca_config%l_ukca_persist_off) ) THEN

    ALLOCATE (biogenic_isop(row_length, rows))
    biogenic_isop(:,:) = 0.0

  END IF  !  biogenic_isop

END IF

ncdf_emissions(:)%l_update = .FALSE. ! Reset flag after use. This will be set
                                     ! by parent at appropriate timestep.

! ----------------------------------------------------------------------
! Deal with online emissions, which are not read from NetCDF files and
! are always updated at each time step. For the moment only NO2 from
! lightning, CH4 from wetlands and sea-air DMS flux are considered here.
!
IF (ukca_config%i_ukca_light_param /= i_light_param_off) THEN
  ! Call the routine to diagnose NO2 lightning emissions.
  lightningem_n_gridbox  = 0.0
  lightningem_no2_to_air = 0.0

  ! Set up integer land/sea mask (needed for lightning emissions)
  DO j = 1, rows
    DO i = 1, row_length
      IF (ls_mask (i, j)) THEN
        ils_mask (i, j) = 1
      ELSE
        ils_mask (i, j) = 0
      END IF
    END DO
  END DO

  CALL ukca_light_ctl (                                                        &
        rows, row_length, delta_lambda, delta_phi, model_levels,               &
        conv_cloud_base,                                                       &
        conv_cloud_top,                                                        &
        ils_mask,                                                              &
        latitude,                                                              &
        surf_area,                                                             &
        r_theta_levels,                                                        &
        r_rho_levels,                                                          &
        p_theta_levels,                                                        &
        p_layer_boundaries,                                                    &
        ext_cg_flash,                                                          &
        ext_ic_flash,                                                          &
        lightningem_n_gridbox,                                                 &
        lightningem_no2_to_air)

  ! ASAD diagnostics for lightning NOx emissions, in kg(NO2)/kg(air)/s. Note:
  !  lightningem_no2_to_air ==> lightning NOx emissions, in kg(NO2)/kg(air)/s
  !  lightning_emissions    ==> type flag used in ASAD diagnostic routines
  IF (L_asad_use_chem_diags .AND. L_asad_use_light_ems) THEN
    CALL asad_3d_emissions_diagnostics (                                       &
                row_length,                                                    &
                rows,                                                          &
                model_levels,                                                  &
                inox,                                                          &
                lightningem_no2_to_air,                                        &
                surf_area,                                                     &
                mass,                                                          &
                m_no2,                                                         &
                timestep,                                                      &
                lightning_emissions,                                           &
                ierr)
  END IF

  ! Convert the field 'lightningem_n_gridbox' from values calculated as
  ! 'kg(N) gridbox-1 s-1' to kg(NO) m^-2 s^-1 as expected for actual emissions
  DO ilev = 1, model_levels
    emissions(inox_light)%values (:,:,ilev) = lightningem_n_gridbox (:,:,ilev) &
                                              * (m_no / m_n)
    emissions(inox_light)%values (:,:,ilev) =                                  &
                       emissions(inox_light)%values (:,:,ilev)/ surf_area (:,:)
  END DO

END IF

! Set wetland methane emissions to zero over non-land surfaces
WHERE (ch4_wetl_emiss < 0.0) ch4_wetl_emiss = 0.0

! Fill in CH4 wetland emissions (if used) at surface level in the
! emissions structure and indicate that they are updated.
IF (ukca_config%l_ukca_qch4inter .AND.                                         &
    (.NOT. ukca_config%l_ukca_prescribech4)) THEN
  ! Scale wetland emission rates by the land fraction and convert from
  ! ug(C) (as received from JULES) to kg(CH4)
  emissions(ich4_wetl)%values (:,:,1) = ch4_wetl_emiss (:,:) *                 &
                                         land_fraction(:,:)  *                 &
                                         1.0e-09 * m_ch4 / m_c
END IF

! Interactive BVOC emissions
IF (ukca_config%l_ukca_ibvoc) THEN

  ! Regrid iBVOC emissions from landpoints (in JULES) to 2D-grid.
  ! If field is active, convert from kg(C) to kg(species) and store as
  ! surface emission

  IF (ic5h8_ibvoc > 0) THEN
    base_scaling = m_c5h8 / (5.0 * m_c)
    CALL regrid_jules_emi(row_length, rows, land_points, land_index,           &
                         ibvoc_isoprene, land_fraction,                        &
                         emissions(ic5h8_ibvoc)%values (:,:,1))  ! surface emiss
    emissions(ic5h8_ibvoc)%values (:,:,1) =                                    &
                  emissions(ic5h8_ibvoc)%values (:,:,1) * base_scaling
  END IF

  IF (ic10h16_ibvoc > 0) THEN
    base_scaling = m_monoterp / (10.0 * m_c)
    emissions(ic10h16_ibvoc)%values (:,:,:) = 0.0
    CALL regrid_jules_emi(row_length, rows, land_points, land_index,           &
                         ibvoc_terpene, land_fraction,                         &
                         emissions(ic10h16_ibvoc)%values (:,:,1))
    emissions(ic10h16_ibvoc)%values (:,:,1) =                                  &
             emissions(ic10h16_ibvoc)%values (:,:,1) * base_scaling
  END IF

  IF (ich3oh_ibvoc > 0) THEN
    base_scaling = m_ch3oh / m_c
    CALL regrid_jules_emi(row_length, rows, land_points, land_index,           &
                         ibvoc_methanol, land_fraction,                        &
                         emissions(ich3oh_ibvoc)%values (:,:,1))
    emissions(ich3oh_ibvoc)%values (:,:,1) =                                   &
           emissions(ich3oh_ibvoc)%values (:,:,1) * base_scaling
  END IF

  IF (ich3coch3_ibvoc > 0) THEN
    base_scaling = m_me2co / (3.0 * m_c)
    CALL regrid_jules_emi(row_length, rows, land_points, land_index,           &
                         ibvoc_acetone, land_fraction,                         &
                         emissions(ich3coch3_ibvoc)%values (:,:,1))
    emissions(ich3coch3_ibvoc)%values (:,:,1) =                                &
              emissions(ich3coch3_ibvoc)%values (:,:,1) * base_scaling
  END IF
END IF  ! If ukca_ibvoc

! Interactive INFERNO fire emissions
IF (ukca_config%l_ukca_inferno) THEN

  ! Linear increasing emission profile with model number, this results
  ! in higher emissions at top levels following a exponential increase
  ! due to hibrid levels coordinate system
  DO ilev = 1, model_levels
    IF (ilev <= ukca_config%i_inferno_emi) THEN
      inferno_vert_fact(ilev) = ilev *                                         &
        (2.0 /(ukca_config%i_inferno_emi * (ukca_config%i_inferno_emi + 1.0 )))
    ELSE
      inferno_vert_fact(ilev) = 0.0
    END IF
  END DO

  ! Gas-phase emissions (aerosols follow further down)
  IF (ico_inferno > 0) THEN
    CALL regrid_jules_emi(row_length, rows, land_points, land_index,           &
                         inferno_co, land_fraction,                            &
                         emissions(ico_inferno)%values (:,:,1))  ! surface emiss

  END IF

  IF (ich4_inferno > 0) THEN
    ! CH4 emissions can either be prescribed using the logical
    ! l_ukca_prescribech4 or interactive - e.g fire and wetland CH4
    ! emissions. l_ukca_inferno_ch4 makes sure that we only add
    ! interactive fire CH4 if l_ukca_inferno_ch4 is TRUE
    IF (ukca_config%l_ukca_inferno_ch4) THEN
      CALL regrid_jules_emi(row_length, rows, land_points, land_index,         &
                        inferno_ch4, land_fraction,                            &
                        emissions(ich4_inferno)%values (:,:,1))  ! surface emiss
    ELSE
      emissions(ich4_inferno)%values (:,:,1) = 0.0
    END IF

  END IF

  IF (inox_inferno > 0) THEN
    CALL regrid_jules_emi(row_length, rows, land_points, land_index,           &
                      inferno_nox, land_fraction,                              &
                      emissions(inox_inferno)%values (:,:,1))  ! surface emiss

  END IF

END IF ! (l_ukca_inferno)

! Call DMS flux routine if marine DMS emissions are required
IF (ukca_config%l_seawater_dms) THEN
  CALL dms_flux_4a(row_length, rows,                                           &
           u_scalar_10m, tstar, land_fraction, dms_sea_conc,                   &
           ukca_config%i_ukca_dms_flux, f_dms_sea)

  IF (ukca_config%l_ukca_scale_seadms_ems) THEN
    ! Multiplication by sea DMS emission factor.
    f_dms_sea = ukca_config%seadms_ems_scaling*f_dms_sea
  END IF
  ! Scale the emissions by non-land, non-seaice area fraction and
  ! convert from kg(S) to kg(DMS). Need to separate the operations as combining
  ! seems to give a different result
  emissions(idms_seaflux)%values(:,:,1) = f_dms_sea(:,:) *                     &
                           (1.0 - land_fraction(:,:)) * (1.0 - seaice_frac(:,:))
  emissions(idms_seaflux)%values(:,:,1) =                                      &
                           emissions(idms_seaflux)%values(:,:,1) * (m_dms / m_s)
END IF     ! If l_seawater_dms

DO l = 1, num_em_flds

  ! ----------------------------------------------------------------------
  ! Update isoprene emissions if they are diurnally varying.
  ! This is only applied to biogenically emitted isoprene (currently with hourly
  ! scaling = 'diurnal_isopems') but not to isoprene emiss from other sources.

  ! Need to store the updated emission values in another variable to avoid
  ! reapplying the hourly scaling with the call to UKCA_DIURNAL_ISOP_EMS

  IF (emissions(l)%tracer_name == 'C5H8      '      .AND.                      &
      emissions(l)%hourly_fact == 'diurnal_isopems' .AND.                      &
      ukca_config%l_diurnal_isopems                 .AND.                      &
      .NOT. ukca_config%l_ukca_ibvoc) THEN

    ! Only use the biogenic_isop fix if the persistence of 3D arrays is on
    ! otherwise assume that emissions are provided at every timestep
    ! and biogenic_isop is not required
    IF (ukca_config%l_ukca_persist_off) THEN

      CALL ukca_diurnal_isop_ems (row_length, rows,                            &
           emissions(l)%values(:,:,1),  cos_zenith_angle, int_zenith_angle,    &
           sin_latitude, cos_latitude, tan_latitude,                           &
           timestep, tmp_out_em_field, testdcycl)

    ELSE

      IF ( emissions(l)%l_update ) THEN
        biogenic_isop = emissions(l)%values (:,:,1)  ! always 2D field
      END IF

      ALLOCATE (tmp_in_em_field (row_length, rows))
      tmp_in_em_field (:,:) = biogenic_isop (:,:)

      CALL ukca_diurnal_isop_ems (row_length, rows,                            &
           tmp_in_em_field,  cos_zenith_angle, int_zenith_angle,               &
           sin_latitude, cos_latitude, tan_latitude,                           &
           timestep, tmp_out_em_field, testdcycl)

      DEALLOCATE (tmp_in_em_field)

    END IF  ! ukca_config%l_ukca_persist_off

    ! Update the emission field
    emissions(l)%values(:,:,1) = tmp_out_em_field (:,:)
  END IF

END DO  ! loop through number of emission fields

! --------------------------------
! Aerosol emissions into modes:
! Online emissions first (fire, sea salt, dust, marine OC) which return number
! and mass for each mode, and these are mapped onto a seperate emissions
! structure for each mode/moment.
! Then offline emissions, which are provided as total mass and are projected
! onto number and mass for each mode by ukca_emiss_update_mode, before being
! mapped to emissions structures for each mode/moment in the same way as the
! online emissions. This has to be done after the scaling factors are applied
! to the offline emission structures above.

! Online fire emissions from INFERNO
IF (ukca_config%l_ukca_inferno) THEN

  IF (iso2_inferno > 0) THEN
    DO ilev = 1, model_levels
      CALL regrid_jules_emi(row_length, rows, land_points, land_index,         &
                            inferno_so2, land_fraction,                        &
                            emissions(iso2_inferno)%values (:,:,ilev))

      emissions(iso2_inferno)%values(:,:,ilev) =                               &
        emissions(iso2_inferno)%values(:,:,ilev) * inferno_vert_fact(ilev)
    END DO

    aer_emmas(:,:,:,:) = 0.0
    aer_emnum(:,:,:,:) = 0.0

    CALL ukca_emiss_update_mode(row_length, rows, model_levels,                &
                                emissions(iso2_inferno)%tracer_name,           &
                                emissions(iso2_inferno)%three_dim,             &
                                emissions(iso2_inferno)%values(:,:,:),         &
                                aer_emmas, aer_emnum, lmode_emiss,             &
                                emiss_levs, icp)

    emissions(iso2_inferno)%values (:,:,:) = 0.0
    CALL ukca_emiss_mode_map(row_length, rows, model_levels,                   &
                             aer_emmas, aer_emnum, emiss_levs, icp,            &
                             lmode_emiss, iso2_inferno)
  END IF

  IF (ioc_inferno > 0) THEN
    DO ilev = 1, model_levels
      CALL regrid_jules_emi(row_length, rows, land_points, land_index,         &
                            inferno_oc, land_fraction,                         &
                            emissions(ioc_inferno)%values (:,:,ilev))

      emissions(ioc_inferno)%values(:,:,ilev) =                                &
        emissions(ioc_inferno)%values(:,:,ilev) * inferno_vert_fact(ilev)
    END DO

    aer_emmas(:,:,:,:) = 0.0
    aer_emnum(:,:,:,:) = 0.0

    CALL ukca_emiss_update_mode(row_length, rows, model_levels,                &
                                emissions(ioc_inferno)%tracer_name,            &
                                emissions(ioc_inferno)%three_dim,              &
                                emissions(ioc_inferno)%values(:,:,:),          &
                                aer_emmas, aer_emnum, lmode_emiss,             &
                                emiss_levs, icp)

    emissions(ioc_inferno)%values (:,:,:) = 0.0
    CALL ukca_emiss_mode_map(row_length, rows, model_levels,                   &
                             aer_emmas, aer_emnum, emiss_levs, icp,            &
                             lmode_emiss, ioc_inferno)
  END IF

  IF (ibc_inferno > 0) THEN
    DO ilev = 1, model_levels
      CALL regrid_jules_emi(row_length, rows, land_points, land_index,         &
                            inferno_bc, land_fraction,                         &
                            emissions(ibc_inferno)%values (:,:,ilev))

      emissions(ibc_inferno)%values(:,:,ilev) =                                &
        emissions(ibc_inferno)%values(:,:,ilev) * inferno_vert_fact(ilev)
    END DO

    aer_emmas(:,:,:,:) = 0.0
    aer_emnum(:,:,:,:) = 0.0

    CALL ukca_emiss_update_mode(row_length, rows, model_levels,                &
                                emissions(ibc_inferno)%tracer_name,            &
                                emissions(ibc_inferno)%three_dim,              &
                                emissions(ibc_inferno)%values(:,:,:),          &
                                aer_emmas, aer_emnum, lmode_emiss,             &
                                emiss_levs, icp)

    emissions(ibc_inferno)%values (:,:,:) = 0.0
    CALL ukca_emiss_mode_map(row_length, rows, model_levels,                   &
                             aer_emmas, aer_emnum, emiss_levs, icp,            &
                             lmode_emiss, ibc_inferno)
  END IF

END IF  ! (l_ukca_inferno)

! Online emissions for sea salt
IF (glomap_config%l_ukca_primss) THEN
  !
  ! In the soluble/insoluble configuration of GLOMAP, seasalt aerosols
  ! are emitted into the soluble component, which is hosted by sulfate
  !
  IF ( glomap_config%i_mode_setup == i_solinsol_6mode ) THEN
    this_cp_cl = cp_su
  ELSE
    this_cp_cl = cp_cl
  END IF
  aer_emmas(:,:,:,:) = 0.0
  aer_emnum(:,:,:,:) = 0.0
  CALL ukca_prim_ss(row_length, rows, model_levels,                            &
                    verbose, land_fraction, seaice_frac, u_scalar_10m,         &
                    aer_emmas, aer_emnum)

  ! Convert mass emission to kg(component)/m2/s
  mm_da = avogadro*boltzmann/rgas  ! Molar mass of dry air (kg/mol)
  aer_emmas = aer_emmas * mm(cp_cl) / mm_da
  lmode_emiss(:) = component(:,this_cp_cl) ! Seasalt emiss covers all ss modes.
  emiss_levs = 1                           ! Surface emissions.

  ! Map num and mass arrays from ukca_prim_ss to corresponding emissions
  ! structures.
  CALL ukca_emiss_mode_map(row_length, rows, model_levels,                     &
                           aer_emmas, aer_emnum, emiss_levs, this_cp_cl,       &
                           lmode_emiss, iseasalt_first)
END IF

! Online emissions for primary marine organic carbon.
! Requires sea-salt ems to be turned on
IF (glomap_config%l_ukca_prim_moc .AND. glomap_config%l_ukca_primbcoc .AND.    &
    glomap_config%l_ukca_primss) THEN
  ! Use the seasalt emissions and chlorophyll concentration to calculate
  ! marine organic primary emissions. aer_emass/num are the input seasalt
  ! emissions for each mode, while mass_pmoc is the total mass of marine
  ! organic emission. This is then distributed across modes by
  ! ukca_emission_update_mode, and mapped into the emissions super-array by
  ! ukca_emiss_mode_map. It is done this way to take advantage of
  ! ukca_emiss_update_mode to do the vol/number calculation and use the
  ! frac/size configuration information in ukca_emiss_mode_mod.

    ! Note that chloro_sea is not masked and contains land values.
    ! Land masking is accounted for implicitly here as a result of using
    ! sea-salt emissions as input to ukca_prim_moc
  chla(:,:) = chloro_sea(:,:) * 1e6 ! chlorophyll-a concentration, mg.m-3

  CALL ukca_prim_moc(row_length, rows, model_levels,                           &
                      aer_emmas, aer_emnum, chla, u_scalar_10m,                &
                      mass_pmoc)

  ! re-initialisation to zero after calculation of PMOC flux
  aer_emmas(:,:,:,:) = 0.0
  aer_emnum(:,:,:,:) = 0.0

  CALL ukca_emiss_update_mode(row_length, rows, model_levels,                  &
                              "PMOC      ", l_3dsource, mass_pmoc,             &
                              aer_emmas, aer_emnum, lmode_emiss,               &
                              emiss_levs,icp)

  CALL ukca_emiss_mode_map(row_length, rows, model_levels,                     &
                           aer_emmas, aer_emnum, emiss_levs, icp,              &
                           lmode_emiss, ipmoc_first)

ELSE IF (glomap_config%l_ukca_prim_moc .AND.                                   &
         (.NOT. glomap_config%l_ukca_primss)) THEN
  cmessage = 'Prim. marine OC requires sea-salt emissions to be turned on'
  errcode = 1
  CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)

END IF

! Online emissions for dust
IF (glomap_config%l_ukca_primdu) THEN
  aer_emmas(:,:,:,:) = 0.0
  aer_emnum(:,:,:,:) = 0.0
  CALL ukca_prim_du(row_length, rows, model_levels, ndiv,                      &
                    verbose, dust_flux, aer_emnum, aer_emmas)
  lmode_emiss(:) = component(:,cp_du)  ! Dust emiss covers all du modes.
  emiss_levs = 1                       ! Surface emissions.

  ! Map num and mass arrays from ukca_prim_du to corresponding emissions
  ! structures.
  CALL ukca_emiss_mode_map(row_length, rows, model_levels,                     &
                           aer_emmas, aer_emnum, emiss_levs, cp_du,            &
                           lmode_emiss, idust_first)
END IF ! l_ukca_primdu


! Online emissions for nitrate and ammonium
IF ((glomap_config%l_ukca_fine_no3_prod) .OR.                                  &
    (glomap_config%l_ukca_coarse_no3_prod)) THEN

  ! Allocate arrays needed for nitrate scheme and initialise to zero
  ! This should probably be moved to a subroutine
  ! e.g. initialise_nitrate_arrays()
  ALLOCATE (air_density(1:row_length,1:rows,1:model_levels))
  ALLOCATE (air_burden(1:row_length,1:rows,1:model_levels))
  air_density(:,:,:) = 0.0
  air_burden(:,:,:) = 0.0

  ! Calculate air density and mass burden using trsrce method
  CALL ukca_calc_rho( row_length, rows, model_levels,                          &
                      theta, q, qcl, qcf, exner_rho_levels, rho_r2,            &
                      r_theta_levels, r_rho_levels, air_density, air_burden )
END IF

IF (glomap_config%l_ukca_fine_no3_prod) THEN
  ALLOCATE (aer_emmas_nh4(1:row_length,1:rows,1:model_levels,1:nmodes))
  ALLOCATE (aer_emnum_nh4(1:row_length,1:rows,1:model_levels,1:nmodes))
  aer_emmas(:,:,:,:) = 0.0
  aer_emnum(:,:,:,:) = 0.0
  aer_emmas_nh4(:,:,:,:) = 0.0
  aer_emnum_nh4(:,:,:,:) = 0.0

  ! Fine mode NO3 / NH4 production from NH3 and HONO2
  ! NH3 and HONO2 MMRs are updated in situ
  ! Call UKCA_PROD_NO3_FINE
  CALL ukca_prod_no3_fine(                                                     &
                    row_length,rows,model_levels,timestep,                     &
                    t_theta_levels, rel_humid_frac_clr, air_density,           &
                    air_burden, tracers(:,:,:,                                 &
                                 n_chem_tracers+n_aero_tracers+1:              &
                                 n_chem_tracers+n_aero_tracers+n_mode_tracers),&
                    tracers(:,:,:,n_nh3), tracers(:,:,:,n_hono2),              &
                    dust_div1, dust_div2, dust_div3,                           &
                    aer_emmas, aer_emnum, aer_emmas_nh4, aer_emnum_nh4)

  lmode_emiss(:) = component(:,cp_no3) ! NO3/NH4 emiss covers all sol modes
  emiss_levs = model_levels            ! 3D emissions.

  ! Map num and mass arrays to corresponding emissions structures
  CALL ukca_emiss_mode_map(row_length, rows, model_levels,                     &
                           aer_emmas, aer_emnum, emiss_levs, cp_no3,           &
                           lmode_emiss, ino3_first)

  CALL ukca_emiss_mode_map(row_length, rows, model_levels,                     &
                           aer_emmas_nh4, aer_emnum_nh4, emiss_levs, cp_nh4,   &
                           lmode_emiss, inh4_first)

END IF  ! l_ukca_fine_no3_prod

IF (glomap_config%l_ukca_coarse_no3_prod) THEN
  ALLOCATE (aer_emmas_nacl(1:row_length,1:rows,1:model_levels,1:nmodes))
  ALLOCATE (aer_emmas_dust(1:row_length,1:rows,1:model_levels,1:nmodes))
  aer_emmas(:,:,:,:) = 0.0
  aer_emnum(:,:,:,:) = 0.0
  aer_emmas_nacl(:,:,:,:) = 0.0
  aer_emmas_dust(:,:,:,:) = 0.0

  ! Call UKCA_PROD_NO3_COARSE
  CALL ukca_prod_no3_coarse(                                                   &
                    row_length,rows,model_levels,timestep,                     &
                    t_theta_levels, rel_humid_frac_clr, air_density,           &
                    air_burden, tropopause_level,                              &
                    tracers(:,:,:,n_chem_tracers+n_aero_tracers+1:             &
                                 n_chem_tracers+n_aero_tracers+n_mode_tracers),&
                    tracers(:,:,:,n_hono2), dust_div1, dust_div2,              &
                    dust_div3, dust_div4, dust_div5, dust_div6,                &
                    aer_emmas, aer_emmas_dust, aer_emmas_nacl)

  ! Update NO3 / NACL
  lmode_emiss(:) = component(:,cp_nn)  ! NaNO3 covers acc/coa sol modes
  emiss_levs = model_levels            ! 3D emissions.

  CALL ukca_emiss_mode_map(row_length, rows, model_levels,                     &
                           aer_emmas, aer_emnum, emiss_levs, cp_nn,            &
                           lmode_emiss, inano3_cond)

  lmode_emiss(:) = component(:,cp_cl) ! Seasalt emiss covers all sol modes
                                      ! 3D emissions.

  CALL ukca_emiss_mode_map(row_length, rows, model_levels,                     &
                           aer_emmas_nacl, aer_emnum, emiss_levs, cp_cl,       &
                           lmode_emiss, iseasalt_hno3)

  IF (glomap_config%l_ukca_primdu) THEN
    ! Update DUST
    lmode_emiss(:) = component(:,cp_cl) ! Dust emiss covers accum/coarse

    CALL ukca_emiss_mode_map(row_length, rows, model_levels,                   &
                             aer_emmas_dust, aer_emnum, emiss_levs, cp_du,     &
                             lmode_emiss, idust_hno3)
  END IF

END IF ! l_ukca_coarse_no3_prod

! Deallocate nitrate arrays
IF (ALLOCATED(air_density))    DEALLOCATE (air_density)
IF (ALLOCATED(air_burden))     DEALLOCATE (air_burden)
IF (ALLOCATED(aer_emmas_nh4))  DEALLOCATE (aer_emmas_nh4)
IF (ALLOCATED(aer_emnum_nh4))  DEALLOCATE (aer_emnum_nh4)
IF (ALLOCATED(aer_emmas_dust)) DEALLOCATE (aer_emmas_dust)
IF (ALLOCATED(aer_emmas_nacl)) DEALLOCATE (aer_emmas_nacl)

! Map offline aerosol emissions into mode emissions structures (see block
! comment above).
DO l = 1, num_em_flds
  IF (ANY(aero_ems_species == emissions(l)%tracer_name) .AND.                  &
      (.NOT. emissions(l)%var_name(1:7) == "inferno") .AND.                    &
      ukca_config%l_ukca_mode) THEN

    aer_emmas(:,:,:,:) = 0.0
    aer_emnum(:,:,:,:) = 0.0
    CALL ukca_emiss_update_mode(row_length, rows, model_levels,                &
                                emissions(l)%tracer_name,                      &
                                emissions(l)%three_dim,                        &
                                emissions(l)%values(:,:,:),                    &
                                aer_emmas, aer_emnum, lmode_emiss, emiss_levs, &
                                icp)

    CALL ukca_emiss_mode_map(row_length, rows, model_levels,                   &
                             aer_emmas, aer_emnum, emiss_levs, icp,            &
                             lmode_emiss, l)
  END IF
END DO
! End of aerosol emission into modes
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Inject emissions and do tracer mixing
CALL ukca_add_emiss (                                                          &
  row_length, rows, model_levels, bl_levels,                                   &
  n_tracers, iyear, imonth, iday, ihour, timestep,                             &
  longitude,                                                                   &
  r_theta_levels, r_rho_levels,                                                &
  theta,                                                                       &
  q,                                                                           &
  qcl,                                                                         &
  qcf,                                                                         &
  exner_rho_levels,                                                            &
  rho_r2,                                                                      &
  kent,                                                                        &
  kent_dsc,                                                                    &
  rhokh_rdz,                                                                   &
  dtrdz,                                                                       &
  we_lim,                                                                      &
  t_frac,                                                                      &
  zrzi,                                                                        &
  we_lim_dsc,                                                                  &
  t_frac_dsc,                                                                  &
  zrzi_dsc,                                                                    &
  ml_depth,                                                                    &
  zhsc,                                                                        &
  z_half,                                                                      &
  surf_area,                                                                   &
  mass,                                                                        &
  tracers,                                                                     &
  len_stashwork50, stashwork50)

! ----------------------------------------------------------------------
! If nitrate is on then check aerosols are > 0.
IF ((glomap_config%l_ukca_fine_no3_prod) .OR.                                  &
    (glomap_config%l_ukca_coarse_no3_prod)) THEN
  CALL ukca_no3_check_values(row_length, rows, model_levels,                   &
                             tracers(:,:,:,                                    &
                                 n_chem_tracers+n_aero_tracers+1:              &
                                 n_chem_tracers+n_aero_tracers+n_mode_tracers))
END IF

! ----------------------------------------------------------------------
! Call the emission diagnostics code if any of the diagnostics present
! in the routine GET_EMDIAG_STASH has been selected via stash.
! Chemical diagnostics first, then aerosol diagnostics.

IF (ukca_config%l_enable_diag_um) THEN

  section = UKCA_diag_sect

  IF (sf(156,section) .OR. sf(157,section) .OR. sf(158,section) .OR.           &
      sf(159,section) .OR. sf(160,section) .OR. sf(161,section) .OR.           &
      sf(162,section) .OR. sf(163,section) .OR. sf(164,section) .OR.           &
      sf(165,section) .OR. sf(166,section) .OR. sf(167,section) .OR.           &
      sf(168,section) .OR. sf(169,section) .OR. sf(170,section) .OR.           &
      sf(171,section) .OR. sf(172,section) .OR. sf(211,section) .OR.           &
      sf(212,section) .OR. sf(213,section) .OR. sf(214,section) .OR.           &
      sf(215,section) .OR. sf(216,section) .OR. sf(217,section) .OR.           &
      sf(304,section) .OR. sf(305,section) .OR. sf(306,section) .OR.           &
      sf(307,section) .OR. sf(308,section) .OR. sf(309,section) .OR.           &
      sf(310,section) .OR. sf(311,section) .OR. sf(312,section) .OR.           &
      sf(313,section) .OR. sf(314,section) ) THEN

    CALL ukca_emiss_diags (row_length, rows, model_levels,                     &
                           len_stashwork50, stashwork50)
  END IF

  ! Aerosol diagnostics
  section = stashcode_glomap_sec

  IF (sf(201,section) .OR. sf(202,section) .OR. sf(203,section) .OR.           &
      sf(204,section) .OR. sf(205,section) .OR. sf(206,section) .OR.           &
      sf(207,section) .OR. sf(208,section) .OR. sf(209,section) .OR.           &
      sf(210,section) .OR. sf(211,section) .OR. sf(212,section) .OR.           &
      sf(213,section) .OR. sf(388,section) .OR. sf(389,section) .OR.           &
      sf(575,section) .OR. sf(576,section) .OR.                                &
      sf(577,section) .OR. sf(579,section) .OR. sf(580,section) .OR.           &
      sf(581,section) .OR. sf(582,section) .OR. sf(583,section)) THEN
    CALL ukca_emiss_diags_mode (row_length, rows, model_levels, area,          &
                                len_stashwork38, stashwork38)
  END IF

END IF

! New option to turn the persistence of spatial variables off
! Deallocated 3D arrays in emissions structure
! Arrays should be deallocated in reverse order that they were allocated
IF (ukca_config%l_ukca_persist_off) THEN

  ! Deallocate arrays in emissions structure
  DO k = num_em_flds, 1, -1
    IF (ASSOCIATED(emissions(k)%diags))                                        &
      DEALLOCATE (emissions(k)%diags)
    NULLIFY (emissions(k)%diags)
    IF (ASSOCIATED(emissions(k)%vert_scaling_3d))                              &
      DEALLOCATE (emissions(k)%vert_scaling_3d)
    NULLIFY (emissions(k)%vert_scaling_3d)
    IF (ASSOCIATED(emissions(k)%values))                                       &
      DEALLOCATE (emissions(k)%values)
    NULLIFY (emissions(k)%values)
  END DO

  ! Deallocate arrays in ncdf_emissions structure
  IF (num_cdf_em_flds > 0) THEN
    DO k = num_cdf_em_flds, 1, -1
      IF (ASSOCIATED(ncdf_emissions(k)%diags))                                 &
        DEALLOCATE (ncdf_emissions(k)%diags)
      NULLIFY (ncdf_emissions(k)%diags)
      IF (ASSOCIATED(ncdf_emissions(k)%vert_scaling_3d))                       &
        DEALLOCATE (ncdf_emissions(k)%vert_scaling_3d)
      NULLIFY (ncdf_emissions(k)%vert_scaling_3d)
      IF (ASSOCIATED(ncdf_emissions(k)%values))                                &
        DEALLOCATE (ncdf_emissions(k)%values)
      NULLIFY (ncdf_emissions(k)%values)
    END DO
  END IF

END IF  !  ukca_config%l_ukca_persist_off

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_emiss_ctl

SUBROUTINE regrid_jules_emi(                                                   &
  row_length, rows, land_points, land_index,                                   &
  emission, land_fraction, emission_2D)

USE parkind1,             ONLY: jpim, jprb      ! DrHook
USE yomhook,              ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE
!
! Description:
!     Regrids JULES emission vectors into UM grid

! Subroutine arguments

! Iput arguments with info on model dimensions
INTEGER,          INTENT(IN)    :: row_length
INTEGER,          INTENT(IN)    :: rows
INTEGER,          INTENT(IN)    :: land_points
INTEGER,          INTENT(IN)    :: land_index(land_points)


! Input arguments of JULES emission (1D) and land fraction
REAL,             INTENT(IN)    :: emission       (land_points)
REAL,             INTENT(IN)    :: land_fraction  (1:row_length, 1:rows)

! Outpit argument of regrided (2D) emission
REAL,             INTENT(OUT)   :: emission_2D    (row_length, rows)

! Local variables
INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

INTEGER                         :: ii, jj, ll

CHARACTER(LEN=*), PARAMETER :: RoutineName='REGRID_JULES_EMI'
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialise emission_2D
DO jj = 1, rows
  DO ii = 1, row_length
    emission_2D (ii,jj) = 0.0
  END DO
END DO

! Regrid JULES emission into UM grid
DO ll = 1, land_points

  jj = (land_index (ll) - 1)     / row_length + 1
  ii =  land_index (ll) - (jj-1) * row_length

  emission_2D (ii,jj) = emission (ll)

  ! Apply coastal tile correction
  emission_2D (ii,jj) = emission_2D (ii,jj) * land_fraction(ii,jj)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE regrid_jules_emi

END MODULE ukca_emiss_ctl_mod
