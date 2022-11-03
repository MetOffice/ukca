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
! Description:
!  Main driver routine for chemistry
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_MAIN1.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_chemistry_ctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CHEMISTRY_CTL_MOD'

CONTAINS

SUBROUTINE ukca_chemistry_ctl(                                                 &
                row_length, rows, model_levels, bl_levels,                     &
                theta_field_size,                                              &
                ntracers, ntype, npft,                                         &
                i_month, i_day_number, i_hour,                                 &
                r_minute, secs_per_step,                                       &
                latitude,                                                      &
                longitude,                                                     &
                sinlat,                                                        &
                tanlat,                                                        &
                pres, temp, q,                                                 &
                qcf, qcl, rh,                                                  &
                p_layer_boundaries,                                            &
                r_theta_levels,                                                &
                z_top_of_model,                                                &
                tracer,                                                        &
                all_ntp,                                                       &
                t_surf, dzl, z0m, u_s,                                         &
                drain, crain,                                                  &
                cloud_frac,                                                    &
                photol_rates,                                                  &
                volume, mass,                                                  &
                so4_aitken, so4_accum, soot_fresh, soot_aged,                  &
                ocff_fresh, ocff_aged, biogenic,                               &
                sea_salt_film, sea_salt_jet,                                   &
                land_points, land_index,                                       &
                tile_pts, tile_index, frac_types,                              &
                zbl, surf_hf, seaice_frac, stcon,                              &
                soilmc_lp, fland, laift_lp, canhtft_lp,                        &
                z0tile_lp, t0tile_lp,                                          &
                have_nat3d,                                                    &
                env_ozone3d,                                                   &
                uph2so4inaer,                                                  &
                delso2_wet_h2o2,                                               &
                delso2_wet_o3,                                                 &
                delh2so4_chem,                                                 &
                delso2_drydep,                                                 &
                delso2_wetdep,                                                 &
                so4_sa,                                                        &
                nat_psc,                                                       &
                trop_ch4_mol,                                                  &
                trop_o3_mol,                                                   &
                trop_oh_mol,                                                   &
                strat_ch4_mol,                                                 &
                strat_ch4loss,                                                 &
                atm_ch4_mol,                                                   &
                atm_co_mol,                                                    &
                atm_n2o_mol,                                                   &
                atm_cf2cl2_mol,                                                &
                atm_cfcl3_mol,                                                 &
                atm_mebr_mol,                                                  &
                atm_h2_mol,                                                    &
                len_stashwork,                                                 &
                stashwork,                                                     &
                H_plus_3d_arr                                                  &
                )

USE ukca_um_legacy_mod,   ONLY: pi_over_180, rgas => r,                        &
                                ! For JULES-based atmospheric deposition
                                deposition_from_ukca_chemistry
USE asad_mod,             ONLY: advt, cdt, ihso3_h2o2, ihso3_o3,               &
                                ih2so4_hv, iso2_oh, iso3_o3,                   &
                                jpctr, jpcspf, jpro2,                          &
                                jpdd, jpdw, jpeq,                              &
                                jphk, jppj, jpspec, jpnr,                      &
                                jpspj, jpspt, jptk,                            &
                                ldepd, ldepw,                                  &
                                nadvt, ndepd, ndepw, nnaf,                     &
                                nprkx, ntrkx,                                  &
                                o1d_in_ss, o3p_in_ss, rk,                      &
                                speci, sph2o, sphno3, spj, spt,                &
                                spro2, ctype,                                  &
                                tnd, y, za, specf, nlnaro2, nldepd
USE asad_chem_flux_diags, ONLY: l_asad_use_chem_diags,                         &
                                l_asad_use_drydep,                             &
                                l_asad_use_flux_rxns,                          &
                                l_asad_use_psc_diagnostic,                     &
                                l_asad_use_rxn_rates,                          &
                                l_asad_use_wetdep,                             &
                                asad_chemical_diagnostics,                     &
                                asad_psc_diagnostic
USE ukca_config_defs_mod, ONLY: nr_therm, nr_phot
USE ukca_cspecies,        ONLY: c_species, n_cf2cl2, n_cfcl3,                  &
                                n_ch4, n_co, n_h2, n_h2o2, n_n2o,              &
                                n_hono2, n_mebr, n_n2o, n_o3,                  &
                                nn_ch4, nn_cl, nn_h2o2, nn_h2so4,              &
                                nn_o1d, nn_o3, nn_o3p, nn_oh,                  &
                                n_h2o, nn_so2, c_na_species,                   &
                                n_h2so4
USE asad_findreaction_mod, ONLY: asad_findreaction
USE UKCA_tropopause,      ONLY: L_troposphere
USE ukca_conserve_mod,    ONLY: ukca_conserve
USE ukca_constants,       ONLY: c_h2o, c_hono2, c_o1d, c_o3p,                  &
                                c_co2, fxb, fxc
USE chemistry_constants_mod, ONLY: avogadro, boltzmann

USE ukca_config_specification_mod, ONLY: ukca_config, i_top_2levH2O,           &
                                         i_top_1lev, i_top_BC

USE ukca_raq_diags_mod, ONLY: ukca_raq_diags
USE ukca_ntp_mod,       ONLY: ntp_type, dim_ntp, name2ntpindex
! Modules for RAQ and RAQ-AERO
USE ukca_chemco_raq_mod,    ONLY: ukca_chemco_raq
USE ukca_deriv_raqaero_mod, ONLY: ukca_deriv_raqaero

USE ukca_environment_fields_mod, ONLY: co2_interactive, h2o2_offline,          &
                                       surf_wetness

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper

USE missing_data_mod,      ONLY: rmdi

USE errormessagelength_mod, ONLY: errormessagelength

USE asad_cdrive_mod, ONLY: asad_cdrive
USE ukca_be_drydep_mod, ONLY: ukca_be_drydep
USE ukca_be_wetdep_mod, ONLY: ukca_be_wetdep
USE ukca_ch4_stratloss_mod, ONLY: ukca_ch4_stratloss
USE ukca_chemco_mod, ONLY: ukca_chemco
USE ukca_ddepctl_mod, ONLY: ukca_ddepctl
USE ukca_ddeprt_mod, ONLY: ukca_ddeprt
USE ukca_deriv_mod, ONLY: ukca_deriv
USE ukca_deriv_aero_mod, ONLY: ukca_deriv_aero
USE ukca_deriv_raq_mod, ONLY: ukca_deriv_raq
USE ukca_fracdiss_mod, ONLY: ukca_fracdiss
USE ukca_sediment_mod, ONLY: ukca_sediment
USE ukca_stratf_mod, ONLY: ukca_stratf
USE ukca_wdeprt_mod, ONLY: ukca_wdeprt
USE ukca_topboundary_mod, ONLY: ukca_top_boundary

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length        ! size of UKCA x dimension
INTEGER, INTENT(IN) :: rows              ! size of UKCA y dimension
INTEGER, INTENT(IN) :: model_levels      ! size of UKCA z dimension
INTEGER, INTENT(IN) :: bl_levels         ! no. of boundary layer levels
INTEGER, INTENT(IN) :: theta_field_size  ! no. of points in horizontal
INTEGER, INTENT(IN) :: ntracers          ! no. of tracers
INTEGER, INTENT(IN) :: ntype             ! no. of surface types
INTEGER, INTENT(IN) :: npft              ! no. of plant functional types
INTEGER, INTENT(IN) :: i_month           ! month
INTEGER, INTENT(IN) :: i_day_number      ! day
INTEGER, INTENT(IN) :: i_hour            ! hour
INTEGER, INTENT(IN) :: uph2so4inaer      ! flag for H2SO4 updating

!       Variables for interactive dry deposition scheme
INTEGER, INTENT(IN) :: land_points
INTEGER, INTENT(IN) :: land_index(land_points)
INTEGER, INTENT(IN) :: tile_pts(ntype)
INTEGER, INTENT(IN) :: tile_index(land_points,ntype)

REAL, INTENT(IN) :: r_minute                           ! minute
REAL, INTENT(IN) :: secs_per_step                      ! time step
REAL, INTENT(IN) :: z_top_of_model                     ! top of model (m)
REAL, INTENT(IN) :: latitude(row_length,rows)          ! latitude (degrees)
REAL, INTENT(IN) :: longitude(row_length,rows)         ! longitude (degrees)
REAL, INTENT(IN) :: sinlat(row_length, rows)           ! sin(latitude)
REAL, INTENT(IN) :: tanlat(row_length, rows)           ! tan(latitude)
REAL, INTENT(IN) :: pres(row_length,rows,model_levels) ! pressure
REAL, INTENT(IN) :: p_layer_boundaries(row_length,rows,                        &
                                       0:model_levels) ! pressure
REAL, INTENT(IN) :: r_theta_levels(row_length,rows,                            &
                                   0:model_levels)
REAL, INTENT(IN) :: temp(row_length,rows,model_levels) ! actual temp
REAL, INTENT(IN) :: dzl(row_length, rows, bl_levels)   ! thickness
REAL, INTENT(IN) :: u_s(row_length, rows)              ! ustar
REAL, INTENT(IN) :: z0m(row_length, rows)              ! roughness
REAL, INTENT(IN) :: t_surf(row_length, rows)           ! surface temp
REAL, INTENT(IN) :: drain(row_length,rows,model_levels) ! 3-D LS rain
REAL, INTENT(IN) :: crain(row_length,rows,model_levels) ! 3-D convec
REAL, INTENT(IN) :: volume(row_length,rows,model_levels) ! cell vol.
REAL, INTENT(IN) :: mass(row_length, rows, model_levels) ! cell mass
! 3D pH array
REAL, INTENT(IN) :: H_plus_3d_arr(row_length, rows, model_levels)

! Aerosol MMR / numbers from CLASSIC, used for calculation
! of surface area if heterogeneous chemistry is ON.
! Aerosol MMR for most aerosol types (kg kg-1)
REAL, INTENT(IN) :: so4_aitken    (row_length, rows, model_levels)
REAL, INTENT(IN) :: so4_accum     (row_length, rows, model_levels)
REAL, INTENT(IN) :: soot_fresh    (row_length, rows, model_levels)
REAL, INTENT(IN) :: soot_aged     (row_length, rows, model_levels)
REAL, INTENT(IN) :: ocff_fresh    (row_length, rows, model_levels)
REAL, INTENT(IN) :: ocff_aged     (row_length, rows, model_levels)
REAL, INTENT(IN) :: biogenic      (row_length, rows, model_levels)
! Aerosol numbers for sea-salt (m-3)
REAL, INTENT(IN) :: sea_salt_film (row_length, rows, model_levels)
REAL, INTENT(IN) :: sea_salt_jet  (row_length, rows, model_levels)

REAL, INTENT(IN) :: env_ozone3d(row_length, rows, model_levels) ! O3
REAL, INTENT(IN) :: qcf(row_length, rows, model_levels)  ! qcf
REAL, INTENT(IN) :: qcl(row_length, rows, model_levels)  ! qcl
REAL, INTENT(IN) :: rh(row_length, rows, model_levels)   ! RH frac
REAL, INTENT(IN) :: cloud_frac(row_length, rows, model_levels)
REAL, INTENT(IN) :: so4_sa(row_length,rows,model_levels) ! aerosol surface area

!       Variables for interactive dry deposition scheme

REAL, INTENT(IN) :: frac_types(land_points,ntype)
REAL, INTENT(IN) :: zbl(row_length,rows)
REAL, INTENT(IN) :: surf_hf(row_length,rows)
REAL, INTENT(IN) :: seaice_frac(row_length,rows)
REAL, INTENT(IN) :: stcon(row_length,rows,npft)
REAL, INTENT(IN) :: soilmc_lp(land_points)
REAL, INTENT(IN) :: fland(land_points)
REAL, INTENT(IN) :: laift_lp(land_points,npft)
REAL, INTENT(IN) :: canhtft_lp(land_points,npft)
REAL, INTENT(IN) :: z0tile_lp(land_points,ntype)
REAL, INTENT(IN) :: t0tile_lp(land_points,ntype)

! Mask to limit formation of Nat below specified height
LOGICAL, INTENT(IN) :: have_nat3d(row_length, rows, model_levels)

REAL, INTENT(IN)     :: photol_rates(row_length, rows, model_levels, jppj)
REAL, INTENT(IN OUT) :: q(row_length,rows,model_levels)   ! water vapour
REAL, INTENT(IN OUT) :: tracer(row_length,rows,                                &
                              model_levels,ntracers)     ! tracer MMR

! Non transported prognostics
TYPE(ntp_type), INTENT(IN OUT) :: all_ntp(dim_ntp)

! SO2 increments
REAL, INTENT(IN OUT) :: delSO2_wet_H2O2(row_length,rows,                       &
                                       model_levels)
REAL, INTENT(IN OUT) :: delSO2_wet_O3(row_length,rows,model_levels)
REAL, INTENT(IN OUT) :: delh2so4_chem(row_length,rows,model_levels)
REAL, INTENT(IN OUT) :: delSO2_drydep(row_length,rows,model_levels)
REAL, INTENT(IN OUT) :: delSO2_wetdep(row_length,rows,model_levels)

! Nitric acid trihydrate (kg(nat)/kg(air))
REAL, INTENT(IN OUT) :: nat_psc(row_length,rows,model_levels)

! Trop CH4 burden (moles)
REAL, INTENT(IN OUT) :: trop_ch4_mol(row_length,rows,model_levels)

! Trop O3 burden (moles)
REAL, INTENT(IN OUT) :: trop_o3_mol(row_length,rows,model_levels)

! Trop OH burden (moles)
REAL, INTENT(IN OUT) :: trop_oh_mol(row_length,rows,model_levels)

! Strat CH4 burden (moles)
REAL, INTENT(IN OUT) :: strat_ch4_mol(row_length,rows,model_levels)

! Strat CH4 loss (Moles/s)
REAL, INTENT(IN OUT) :: strat_ch4loss(row_length,rows,model_levels)

! Atmospheric Burden of CH4
REAL, INTENT(IN OUT) :: atm_ch4_mol(row_length,rows,model_levels)

! Atmospheric Burden of CO
REAL, INTENT(IN OUT) :: atm_co_mol(row_length,rows,model_levels)

! Atmospheric Burden of Nitrous Oxide (N2O)
REAL, INTENT(IN OUT) :: atm_n2o_mol(row_length,rows,model_levels)

! Atmospheric Burden of CFC-12
REAL, INTENT(IN OUT) :: atm_cf2cl2_mol(row_length,rows,model_levels)

! Atmospheric Burden of CFC-11
REAL, INTENT(IN OUT) :: atm_cfcl3_mol(row_length,rows,model_levels)

! Atmospheric Burden of CH3Br
REAL, INTENT(IN OUT) :: atm_mebr_mol(row_length,rows,model_levels)

! Atmospheric Burden of H2
REAL, INTENT(IN OUT) :: atm_h2_mol(row_length,rows,model_levels)

! Diagnostics array
INTEGER, INTENT(IN) :: len_stashwork

REAL, INTENT(IN OUT) :: stashwork (len_stashwork)

! Local variables
INTEGER :: nlev_with_ddep(row_length, rows)     ! No levs in bl
INTEGER :: nlev_with_ddep2(theta_field_size)    ! No levs in bl

INTEGER, SAVE :: nr            ! no of rxns for BE
INTEGER, SAVE :: n_be_calls    ! no of call to BE solver

INTEGER :: ix                  ! dummy variable
INTEGER :: jy                  ! dummy variable
INTEGER :: i             ! Loop variable
INTEGER :: j             ! loop variable
INTEGER :: js            ! loop variable
INTEGER :: jtr           ! loop variable - transported tracers
INTEGER :: jro2          ! loop variable - NTP RO2 species
INTEGER :: jro2_copy     ! copy of loop variable jro2 for error reporting
INTEGER :: jna           ! loop variable, non-advected species
INTEGER :: jspf          ! loop variable - all active chemical species in f
INTEGER :: k             ! loop variable
INTEGER :: l             ! loop variable
INTEGER :: n_pnts        ! no. of pts in 2D passed to CDRIVE

INTEGER, SAVE :: istore_h2so4  ! location of H2SO4 in f array

INTEGER           :: ierr                     ! Error code: asad diags routines
INTEGER           :: errcode                  ! Error code: ereport
CHARACTER(LEN=errormessagelength) :: cmessage         ! Error message
CHARACTER(LEN=10)      :: prods(2)                 ! Products
CHARACTER(LEN=10)      :: prods3(3)                ! Products
LOGICAL           :: ddmask(theta_field_size) ! mask

REAL, SAVE      :: dts                       ! B. Euler timestep

REAL :: tgmt               ! GMT time (decimal represent
REAL :: tan_declin         ! TAN(declination)

!Dummy variables for compatability with column-call approach
REAL :: dpd_dummy(model_levels,jpspec)
REAL :: dpw_dummy(model_levels,jpspec)
REAL :: prk_dummy(model_levels,jpnr)
REAL :: y_dummy(model_levels,jpspec)
REAL :: fpsc1_dummy(model_levels)
REAL :: fpsc2_dummy(model_levels)

! SO2 increments in molecules/cm^3
REAL :: SO2_wetox_H2O2(theta_field_size)
REAL :: SO2_wetox_O3(theta_field_size)
REAL :: SO2_dryox_OH(theta_field_size)

REAL, ALLOCATABLE :: BE_rc(:,:)       ! 1-D Rate coeff array
REAL, ALLOCATABLE :: BE_hrc(:,:)      ! 1-D Rate coeff array (heterog reactions)
REAL, ALLOCATABLE :: zfnatr(:,:)      ! 1-D array of non-transported tracers
REAL, ALLOCATABLE :: ystore(:)        ! array for H2SO4 when updated in MODE
REAL :: zftr(theta_field_size,jpcspf) ! 1-D array of chemically active species
                                      !   including RO2 species, in VMR
REAL :: zp  (theta_field_size)        ! 1-D pressure
REAL :: zt  (theta_field_size)        ! 1-D temperature
REAL :: zclw(theta_field_size)        ! 1-D cloud liquid water
REAL :: zfcloud(theta_field_size)     ! 1-D cloud fraction
REAL :: cdot(theta_field_size,jpcspf)  ! 1-D chem. tendency
REAL :: zq(theta_field_size)          ! 1-D water vapour vmr
REAL :: co2_1d(theta_field_size)      ! 1-D CO2
REAL :: zprt1d(theta_field_size,jppj) ! 1-D photolysis rates for ASAD
REAL :: tloc  (row_length, rows)      ! local time
REAL :: daylen(row_length, rows)      ! local daylength
REAL :: cs_hour_ang(row_length, rows) ! cosine hour angle
REAL :: zdryrt(row_length, rows, jpdd)                ! dry dep rate
REAL :: zdryrt2(theta_field_size, jpdd)               ! dry dep rate
REAL :: zwetrt(row_length, rows, model_levels, jpdw)  ! wet dep rate
REAL :: zwetrt2(theta_field_size, jpdw)               ! wet dep rat
REAL :: zfrdiss2(theta_field_size,jpdw,jpeq+1)        ! dissolved fraction
REAL :: zfrdiss(row_length, rows, model_levels, jpdw, jpeq+1)
REAL :: rc_het(theta_field_size,2)                ! heterog rates for trop chem
REAL :: kp_nh(row_length, rows, model_levels)     ! Dissociation const
REAL :: BE_tnd(theta_field_size)                  ! total no density, molec cm-3
REAL :: BE_h2o(theta_field_size)                  ! water vapour concn
REAL :: BE_o2(theta_field_size)                   ! oxygen concn
REAL :: BE_vol(theta_field_size)                  ! gridbox volume
REAL :: BE_rho(theta_field_size)                  ! air density (kg m-3)
REAL :: BE_rh_frac(theta_field_size)              ! RH (fraction: 0.000-0.999)
REAL :: H_plus_2d_arr(theta_field_size,model_levels) ! 2-D pH array to use in
                                                     ! wet deposition
REAL :: H_plus_1d_arr(theta_field_size)  ! 1-D pH array to use in chemical
                                         ! solvers

! Aerosol mmr/numbers from CLASSIC
REAL :: BE_so4_aitken    (theta_field_size) ! MMR (kg kg-1)
REAL :: BE_so4_accum     (theta_field_size)
REAL :: BE_soot_fresh    (theta_field_size)
REAL :: BE_soot_aged     (theta_field_size)
REAL :: BE_ocff_fresh    (theta_field_size)
REAL :: BE_ocff_aged     (theta_field_size)
REAL :: BE_biogenic      (theta_field_size)
REAL :: BE_sea_salt_film (theta_field_size) ! number (m-3)
REAL :: BE_sea_salt_jet  (theta_field_size)

REAL :: BE_wetrt(theta_field_size,jpspec)         ! wet dep rates (s-1)
REAL :: BE_dryrt(theta_field_size,jpspec)         ! dry dep rates (s-1)
REAL :: BE_frdiss(theta_field_size,jpspec,jpeq+1) ! dissolved fraction
! concentrations for backward euler solver in volume mixing ratio
REAL :: BE_y (theta_field_size,jpspec)
! Local stratospheric CH4 loss rate
REAL :: strat_ch4loss_2d(theta_field_size,model_levels)
REAL :: k_dms(theta_field_size,5)                 ! dms rate coeffs

!     Dry and wet deposition fluxes (mol s-1)
REAL :: ddflux(theta_field_size, jpdd) ! dry deposition flux
REAL :: wdflux(theta_field_size, jpdw) ! wet deposition flux
REAL :: dry_dep_3d(row_length, rows, model_levels, jpdd) ! 3d dry dep
REAL :: wet_dep_3d(row_length, rows, model_levels, jpdw) ! 3d wet dep

LOGICAL, SAVE :: firstcall = .TRUE.


! The calls to ukca_conserve require a logical to be set.
! ukca_conserve calculates and conserves total chlorine, bromine, and
! hydrogen. For these elements closed chemistry should be prescribed.
! Called before chemistry, before_chem, it calculates
! total bromine, chlorine, and hydrogen as 3-D fields. Called afer
! chemistry, after_chem, it rescales the chlorine, bromine
! and hydrogen containing compounds so that total chlorine, bromine
! and hydrogen are conserved under chemistry.
LOGICAL, PARAMETER :: before_chem = .TRUE.
LOGICAL, PARAMETER :: after_chem = .FALSE.

! Variables for heterogeneous chemistry
REAL, ALLOCATABLE :: shno3_3d(:,:,:)
! 1-D masks for troposphere and NAT height limitation
LOGICAL :: stratflag(theta_field_size)
LOGICAL :: have_nat1d(theta_field_size)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMISTRY_CTL'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
n_pnts = rows * row_length
nr = nr_therm + nr_phot

! dummy variables for compatability with column call
ix=0
jy=0
dpd_dummy=0.0
dpw_dummy=0.0
fpsc1_dummy=0.0
fpsc2_dummy=0.0
prk_dummy=0.0
y_dummy=0.0

IF (firstcall) THEN

  !         Check that theta_field_size = n_pnts

  IF (theta_field_size /= n_pnts) THEN
    cmessage='theta_field_size not equal to n_pnts'
    CALL ereport('UKCA_CHEMISTRY_CTL',n_pnts,cmessage)
  END IF

  ! Backward Euler timestep variables

  n_be_calls = INT(secs_per_step/ukca_config%dts0)
  ! Ensure we call BE at least once
  IF (n_be_calls == 0) n_be_calls = 1
  ! Calculate the BE Timestep
  dts = secs_per_step/n_be_calls

  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'(A,I0,E12.4)') 'n_be_calls, dts= ',n_be_calls, dts
    CALL umPrint(umMessage,src='ukca_chemistry_ctl')
  END IF

  !         Check whether water vapour is advective tracer. Then,
  !         check whether UM and ASAD advective tracers correspond
  !         to each other.

  IF ((ukca_config%l_ukca_advh2o) .AND. (n_h2o == 0)) THEN
    cmessage='No tracer for advected water vapour'
    errcode = 4
    CALL ereport('UKCA_CHEMISTRY_CTL',errcode,cmessage)
  END IF

  ! Identify the SO2+OH rate coeff, the products are alternatives depending
  !  on whether the H2SO4 tracer updating is to be done in ASAD or in MODE.
  iso2_oh = 0
  ih2so4_hv = 0
  istore_h2so4 = 0
  IF (ukca_config%l_ukca_nr_aqchem) THEN
    ! find location of H2SO4 in zftr array
    DO jspf = 1, jpcspf
      IF (specf(jspf) == advt(n_h2so4)) THEN
        istore_h2so4 = jspf
      END IF
    END DO

    IF (ukca_config%l_ukca_offline) THEN
      prods = ['H2SO4     ','          ']
    ELSE
      prods = ['H2SO4     ','HO2       ']
    END IF
    iso2_oh = asad_findreaction( 'SO2       ', 'OH        ',                   &
                           prods, 2, spt, ntrkx, jptk+1, jpspt )

    IF (iso2_oh == 0) THEN   ! check for stratospheric sulphur chemistry
      ! product should be HO2 not H2O
      ! not considered when l_fix_ukca_h2so4_ystore=T
      prods = ['SO3       ','HO2       ']
      iso2_oh = asad_findreaction( 'SO2       ', 'OH        ',                 &
                           prods, 2, spt, ntrkx, jptk+1, jpspt )
      IF (ukca_config%i_ukca_chem_version >= 121) THEN
         ! additional product at version=121
        prods3 = ['SO3       ','OH        ','H         ']
        ih2so4_hv = asad_findreaction( 'H2SO4     ', 'PHOTON    ',             &
                             prods3, 3, spj, nprkx, jppj+1, jpspj )
      ELSE
        prods = ['SO3       ','OH        ']
        ih2so4_hv = asad_findreaction( 'H2SO4     ', 'PHOTON    ',             &
                             prods, 2, spj, nprkx, jppj+1, jpspj )
      END IF

    END IF

    IF (iso2_oh == 0 .AND. ih2so4_hv == 0) THEN
      cmessage=' Sulphur chemistry reactions not found'
      WRITE(umMessage,'(A)') cmessage
      CALL umPrint(umMessage,src='ukca_chemistry_ctl')
      WRITE(umMessage,'(A,I0,A,I0)') 'iso2_oh: ',iso2_oh,' ih2so4_hv: ',       &
                                     ih2so4_hv
      CALL umPrint(umMessage,src='ukca_chemistry_ctl')
      errcode = 1
      CALL ereport('UKCA_CHEMISTRY_CTL',errcode,cmessage)
    END IF
  END IF   ! l_ukca_nr_aqchem

END IF  ! of initialization of chemistry subroutine (firstcall)


!       Calculate local time as function of longitude

tgmt = REAL(i_hour) + r_minute/60.0                                            &
                    + secs_per_step * 0.5 / 3600.0
!        IF (tgmt < 0.) tgmt = tgmt + 24.

tloc = tgmt + 24.0 * longitude/360.0
WHERE (tloc > 24.0) tloc = tloc - 24.0

!       Calculate Declination Angle and Daylength for each row for
!       current day of the year.
!       Ensure COS of HOUR ANGLE does not exceed + or - 1
daylen = 0.0
tan_declin = TAN(fxb * SIN(pi_over_180*(266.0+i_day_number)))
IF (ABS(tan_declin) < EPSILON(0.0)) THEN
  cs_hour_ang = 0.0
ELSE
  cs_hour_ang = MAX(-1.0, MIN(1.0, -tanlat * tan_declin))
END IF
daylen = fxc * ACOS(cs_hour_ang)

! Call routine to calculate dry deposition rates.
zdryrt  = 0.0
zdryrt2 = 0.0
IF (ndepd /= 0) THEN

  IF (ukca_config%l_ukca_intdd) THEN           ! Call interactive dry dep

    IF (ukca_config%l_deposition_jules) THEN   ! Use JULES-based routines

      CALL deposition_from_ukca_chemistry(                                     &
        secs_per_step, bl_levels, row_length, rows, ntype, npft,               &
        jpspec, ndepd, nldepd, speci,                                          &
        land_points, land_index, tile_pts, tile_index,                         &
        seaice_frac, fland, sinlat,                                            &
        p_layer_boundaries(:,:,0), rh(:,:,1), t_surf, surf_hf, surf_wetness,   &
        z0tile_lp, stcon, laift_lp, canhtft_lp, t0tile_lp,                     &
        soilmc_lp, zbl, dzl, frac_types, u_s,                                  &
        zdryrt, nlev_with_ddep, len_stashwork, stashwork)

    ELSE                                       ! Use exising UKCA routines

      CALL ukca_ddepctl(row_length, rows, bl_levels, ntype, npft,              &
        land_points, land_index, tile_pts, tile_index,                         &
        secs_per_step, sinlat, frac_types, t_surf,                             &
        p_layer_boundaries(:,:,0), dzl, zbl, surf_hf, u_s,                     &
        rh, stcon, soilmc_lp, fland, seaice_frac, laift_lp,                    &
        canhtft_lp, z0tile_lp, t0tile_lp,                                      &
        nlev_with_ddep, zdryrt, len_stashwork, stashwork)

    END IF

  ELSE                             ! Call prescribed dry dep

    CALL ukca_ddeprt(daylen, tloc, n_pnts, dzl, bl_levels,                     &
                             z0m, u_s, t_surf,                                 &
                             latitude, i_month,                                &
                             1, n_pnts,                                        &
                             zdryrt)

  END IF
END IF

!       Call routine to calculate wet deposition rates.

zwetrt  = 0.0
zwetrt2 = 0.0
IF (ndepw /= 0) THEN

  ! Reshape 3D H_plus array to 2D to use in calculating Wet Deposition rates
  DO k=1,model_levels
    H_plus_2d_arr(:,k) = RESHAPE(H_plus_3d_arr(:,:,k),[n_pnts])
  END DO

  CALL ukca_wdeprt(n_pnts, model_levels, drain, crain, temp,                   &
                   latitude, secs_per_step, zwetrt, H_plus_2d_arr)
END IF

! Calculate dissolved fraction (only used in online chem with B-E Solver)
IF (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) THEN
  ! Send 3D H_plus array to calculate fraction dissolved
  CALL ukca_fracdiss(row_length, rows, model_levels,                           &
                     temp, pres, rh, qcl, zfrdiss, kp_nh, H_plus_3d_arr)
END IF

IF (ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_stratcfc .OR.             &
    ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat) THEN

  ! Calculate total chlorine and total bromine before chemistry

  IF (nn_cl > 0) THEN
    CALL ukca_conserve(row_length, rows, model_levels, ntracers,               &
         tracer, pres, drain, crain, before_chem)
  END IF

END IF    ! l_ukca_strat etc

! if heterogeneous chemistry is selected, allocate solid HNO3 array
IF (ukca_config%l_ukca_het_psc) THEN
  IF (.NOT. ALLOCATED(shno3_3d))                                               &
       ALLOCATE(shno3_3d(row_length, rows, model_levels))
  shno3_3d = 0.0
END IF

!       Initialize budget variables
strat_ch4loss_2d = 0.0

! Need this line here for ASAD interactive DD.
nlev_with_ddep2(:) = RESHAPE(nlev_with_ddep(:,:),[theta_field_size])

!$OMP PARALLEL DEFAULT(NONE)                                                   &
!$OMP PRIVATE(BE_biogenic, BE_dryrt, BE_frdiss, BE_h2o, BE_hrc,                &
!$OMP         BE_o2, BE_ocff_aged, BE_ocff_fresh,                              &
!$OMP         BE_rc, BE_rh_frac, BE_rho,                                       &
!$OMP         BE_sea_salt_film, BE_sea_salt_jet,                               &
!$OMP         BE_so4_aitken, BE_so4_accum,                                     &
!$OMP         BE_soot_aged, BE_soot_fresh,                                     &
!$OMP         BE_tnd, BE_wetrt, BE_vol, BE_y,                                  &
!$OMP         cdot, ddflux, ddmask, have_nat1d, ierr,                          &
!$OMP         jna, jro2, jro2_copy, jspf, k, k_dms, l, rc_het,                 &
!$OMP         SO2_dryox_OH, SO2_wetox_H2O2, SO2_wetox_O3,                      &
!$OMP         stratflag, wdflux, ystore,                                       &
!$OMP         zclw, zdryrt2, zfcloud, zfnatr, zfrdiss2, zftr,                  &
!$OMP         zp, zprt1d, zq, zt, co2_1d, zwetrt2,  H_plus_1d_arr)             &
!$OMP SHARED(advt, all_ntp, atm_cf2cl2_mol, atm_cfcl3_mol, atm_ch4_mol,        &
!$OMP        atm_co_mol, atm_h2_mol, atm_mebr_mol, atm_n2o_mol,                &
!$OMP        biogenic, c_species, c_na_species, cloud_frac,                    &
!$OMP        delh2so4_chem, delSO2_drydep, delSO2_wet_H2O2,                    &
!$OMP        delSO2_wet_O3, delSO2_wetdep, dpd_dummy, dpw_dummy,               &
!$OMP        dry_dep_3d, dts, fpsc1_dummy, fpsc2_dummy, h2o2_offline,          &
!$OMP        have_nat3d, prk_dummy, speci, specf, y_dummy, co2_interactive,    &
!$OMP        ih2so4_hv, ihso3_h2o2, ihso3_o3, iso2_oh, iso3_o3,                &
!$OMP        ix, jpctr, jpdd, jpdw, jphk, jppj, jpspec, jpro2, jy, kp_nh,      &
!$OMP        jpcspf, spro2, nlnaro2, ctype, istore_h2so4,                      &
!$OMP        l_asad_use_chem_diags, l_asad_use_drydep,                         &
!$OMP        l_asad_use_flux_rxns, l_asad_use_psc_diagnostic,                  &
!$OMP        l_asad_use_rxn_rates, l_asad_use_wetdep, l_troposphere,           &
!$OMP        ukca_config,                                                      &
!$OMP        ldepd, ldepw, model_levels,                                       &
!$OMP        n_be_calls, n_cf2cl2, n_cfcl3, n_ch4, n_co, n_h2, n_h2o2,         &
!$OMP        n_mebr, n_n2o, n_o3, n_pnts, nadvt, nlev_with_ddep2,              &
!$OMP        nn_ch4, nn_h2o2, nn_h2so4, nn_o1d, nn_o3, nn_o3p, nn_oh,          &
!$OMP        nn_so2, nnaf, nr, nr_therm,                                       &
!$OMP        o1d_in_ss, o3p_in_ss, ocff_aged, ocff_fresh, photol_rates,        &
!$OMP        pres, q, qcf, qcl, rgas, rh, row_length, rows,                    &
!$OMP        sea_salt_film, sea_salt_jet, shno3_3d,                            &
!$OMP        so4_accum, so4_aitken, so4_sa, soot_aged, soot_fresh,             &
!$OMP        strat_ch4_mol, strat_ch4loss, strat_ch4loss_2d,                   &
!$OMP        temp, theta_field_size, tracer,                                   &
!$OMP        trop_ch4_mol, trop_o3_mol, trop_oh_mol,                           &
!$OMP        uph2so4inaer, volume, wet_dep_3d,                                 &
!$OMP        zdryrt, zfrdiss, zwetrt, H_plus_3d_arr, cmessage)

IF (.NOT. ALLOCATED(ystore) .AND. uph2so4inaer == 1)                           &
                             ALLOCATE(ystore(theta_field_size))

!$OMP DO SCHEDULE(DYNAMIC)
DO k=1,model_levels

  ! Copy water vapour and ice field into 1-D arrays
  IF (ukca_config%l_ukca_het_psc) THEN
    IF (k <= model_levels) THEN
      sph2o(:) = RESHAPE(qcf(:,:,k),[theta_field_size])/c_h2o
    ELSE
      sph2o(:) = 0.0
    END IF
  END IF

  zdryrt2(:,:) = 0.0e0
  IF (ukca_config%l_ukca_intdd) THEN
     ! Interactive scheme extracts from levels in boundary layer
    ddmask(:) = (k <= nlev_with_ddep2(:))
    DO l=1,jpdd
      WHERE (ddmask(:))
        zdryrt2(:,l) =                                                         &
             RESHAPE(zdryrt(:,:,l),[theta_field_size])
      END WHERE
    END DO
  ELSE    ! non-interactive
    IF (k == 1) THEN
      DO l=1,jpdd
        zdryrt2(:,l) =                                                         &
             RESHAPE(zdryrt(:,:,l),[theta_field_size])
      END DO
    END IF
  END IF

  !       Put pressure, temperature and tracer mmr into 1-D arrays
  !       for use in ASAD chemical solver

  zp(:) = RESHAPE(pres(:,:,k),[theta_field_size])
  zt(:) = RESHAPE(temp(:,:,k),[theta_field_size])
  zq(:) = RESHAPE(q(:,:,k),[theta_field_size])/c_h2o

  IF (ANY(speci(:) == 'CO2       ')) THEN
    !  Copy the CO2 concentration into the asad module as VMR
    IF (ukca_config%l_chem_environ_co2_fld) THEN
      co2_1d(:) = RESHAPE(co2_interactive(:,:,k),[theta_field_size])/c_co2
    ELSE
      co2_1d(:) = rmdi
    END IF

  END IF       ! CO2 as species

  IF (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) THEN
    IF (k <= model_levels) THEN
      zclw(:) = RESHAPE(qcl(:,:,k),[theta_field_size])
      zfcloud(:) = RESHAPE(cloud_frac(:,:,k),[theta_field_size])
    ELSE
      zclw(:) = 0.0
      zfcloud(:) = 0.0
    END IF
  END IF

  ! Reduce over-prediction of H2O2 using ancillary value.
  IF (ukca_config%l_ukca_offline) THEN
    WHERE (tracer(:,:,k,n_h2o2) > h2o2_offline(:,:,k))                         &
       tracer(:,:,k,n_h2o2) = h2o2_offline(:,:,k)
  END IF

  ! Convert mmr into vmr for tracers. Pass data from the tracer 3D array,
  ! unwrap it and pass into the 1D zftr array before calling ASAD_CDRIVE.
  ! If running with nontransport RO2 species, data for the RO2 species
  ! needs to be passed from the all_ntp 3D array as well.

  ! First set counter for all chemical species in f array
  jspf = 0
  ! Loop through all species
  DO js = 1,jpspec
    ! First try to map with transported tracer
    DO jtr=1,jpctr
      IF (advt(jtr) == speci(js)) THEN
        jspf = jspf+1
        ! Map data from tracer 3D array into zftr
        zftr(:,jspf) = RESHAPE(tracer(:,:,k,jtr),                              &
                      [theta_field_size])/c_species(jtr)
      END IF
    END DO ! Close advected tracer do loop

    ! If RO2 species are not being transported, search in list of RO2 species
    IF (ukca_config%l_ukca_ro2_ntp) THEN
      DO jro2 = 1, jpro2
        IF (spro2(jro2) == speci(js) .AND. ctype(js) == 'OO') THEN
          ! Get index of speci(js) in NTP array.
          ! If not found this is a fatal error.
          l = name2ntpindex(spro2(jro2))

          ! Find location of species in nadvt array
          jna = nlnaro2(jro2)

          ! Map data from all_ntp 3D array into zftr
          IF (nadvt(jna) == spro2(jro2)) THEN
            jspf = jspf+1
            zftr(:,jspf) = RESHAPE(all_ntp(l)%data_3d(:,:,k),[n_pnts]) /       &
                          c_na_species(jna)
          ELSE
            jro2_copy = jro2
            WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl'
            CALL umPrint(umMessage,src='ukca_chemistry_ctl')
            cmessage='ERROR: Indices for RO2 species do not match with nadvt'
            CALL ereport('UKCA_CHEMISTRY_CTL', jro2_copy, cmessage)
          END IF ! Close IF nadvt and spro2 match

        END IF   ! Close IF RO2 species
      END DO     ! Close loop through RO2 species
    END IF       ! Close IF RO2_NTP
  END DO         ! Close loop through all species

  ! Check we have the correct number of active chemical species
  IF (ukca_config%l_ukca_ro2_ntp) THEN
    IF (jspf /= jpro2+jpctr) THEN
      WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl'
      CALL umPrint(umMessage,src='ukca_chemistry_ctl')
      cmessage = 'ERROR: Number of chemical active species /= jpro2+jpctr'
      CALL ereport('UKCA_CHEMISTRY_CTL', jspf, cmessage)
    END IF
  ELSE
    IF (jspf /= jpctr) THEN
      WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl'
      CALL umPrint(umMessage,src='ukca_chemistry_ctl')
      cmessage = 'ERROR: Number of chemical active species /= jpctr'
      CALL ereport('UKCA_CHEMISTRY_CTL', jspf, cmessage)
    END IF
  END IF


  ! Map photolysis rates onto 1-D array.
  IF (ukca_config%l_ukca_offline) THEN
    ! Offline chemistry has no photolysis
    zprt1d(:,:) = 0.0
  ELSE
    DO l=1,jppj
      zprt1d(:,l) = RESHAPE(photol_rates(:,:,k,l),[theta_field_size])
    END DO
  END IF

  !       Call ASAD routines to do chemistry integration
  !       In lowest levels choose half the dynamical timestep for
  !       chemistry. If dynamical timestep > 20 min, use half and
  !       quarter of dynamical timestep for chemistry.

  IF (.NOT. (ukca_config%l_ukca_trop .OR. ukca_config%l_ukca_aerchem .OR.      &
             ukca_config%l_ukca_raq .OR. ukca_config%l_ukca_raqaero)) THEN
    ! Not B-E

    ! retrieve tropospheric heterogeneous rates from previous time step
    ! for this model level (index k)
    IF (ukca_config%l_ukca_trophet) THEN
      ! N2O5
      l = name2ntpindex('het_n2o5  ')
      rc_het(:,1) = RESHAPE(all_ntp(l)%data_3d(:,:,k),                         &
                            [theta_field_size])
      ! HO2+HO2
      l = name2ntpindex('het_ho2   ')
      rc_het(:,2) = RESHAPE(all_ntp(l)%data_3d(:,:,k),                         &
                            [theta_field_size])
    ELSE
      rc_het(:,:) = 0.0
    END IF

    ! fill stratospheric flag indicator and SO4 surface area
    stratflag(:) = (.NOT. RESHAPE(L_troposphere(:,:,k),                        &
                              [theta_field_size]) )
    za(:) = RESHAPE(so4_sa(:,:,k),[theta_field_size])

    ! Reshape masked array for natpsc formation
    have_nat1d(:) = RESHAPE(have_nat3d(:,:,k),[theta_field_size])

    ! pass 2D wet-dep field to cdrive
    DO l=1,jpdw
      zwetrt2(:,l) = RESHAPE(zwetrt(:,:,k,l),[theta_field_size])
    END DO

    IF (uph2so4inaer == 1) THEN
      ! H2SO4 will be updated in MODE, so store old value here
      IF (ukca_config%l_fix_ukca_h2so4_ystore) THEN
         ! primary array passed is zftr, so save this, NOT y
        ystore(:) = zftr(:,istore_h2so4)
      ELSE
        ystore(:) = y(:,nn_h2so4)
      END IF
    END IF

    ! collate 1-D H_plus values to use in asad NR solver
    H_plus_1d_arr(:) = RESHAPE(H_plus_3d_arr(:,:,k),[n_pnts])

    CALL asad_cdrive(cdot, zftr, zp, zt, zq, co2_1d,                           &
                     RESHAPE(cloud_frac(:,:,k),[n_pnts]),                      &
                     RESHAPE(qcl(:,:,k),[n_pnts]),                             &
                     ix, jy, k,zdryrt2, zwetrt2, rc_het,                       &
                     zprt1d, n_pnts, have_nat1d, stratflag, H_plus_1d_arr)

    IF (ukca_config%l_ukca_het_psc) THEN
      ! Save MMR of NAT PSC particles into 3-D array for PSC sedimentation.
      ! Note that sphno3 is NAT in number density of HNO3.
      IF (ANY(sphno3(1:theta_field_size) > 0.0)) THEN
        shno3_3d(:,:,k) = RESHAPE(sphno3(:)/tnd(:),                            &
                             [row_length,rows])*c_hono2
      ELSE
        shno3_3d(:,:,k) = 0.0
      END IF
    END IF

    IF (ukca_config%l_ukca_chem .AND. ukca_config%l_ukca_nr_aqchem) THEN
      ! Calculate chemical fluxes for MODE
      IF (ihso3_h2o2 > 0) delSO2_wet_H2O2(:,:,k) =                             &
        delSO2_wet_H2O2(:,:,k) + RESHAPE(rk(:,ihso3_h2o2)*                     &
        y(:,nn_so2)*y(:,nn_h2o2),[row_length,rows])*cdt
      IF (ihso3_o3 > 0) delSO2_wet_O3(:,:,k) =                                 &
        delSO2_wet_O3(:,:,k) + RESHAPE(rk(:,ihso3_o3)*                         &
        y(:,nn_so2)*y(:,nn_o3),[row_length,rows])*cdt
      IF (iso3_o3 > 0) delSO2_wet_O3(:,:,k) =                                  &
        delSO2_wet_O3(:,:,k) + RESHAPE(rk(:,iso3_o3)*                          &
        y(:,nn_so2)*y(:,nn_o3),[row_length,rows])*cdt
      ! net H2SO4 production - note that this is affected by
      ! l_fix_ukca_h2so4_ystore above. Y value is concentration
      ! from chemistry prior to zftr being over-written below
      IF (iso2_oh > 0 .AND. ih2so4_hv > 0) THEN
        delh2so4_chem(:,:,k) = delh2so4_chem(:,:,k) +                          &
         (RESHAPE(rk(:,iso2_oh)*y(:,nn_so2)*y(:,nn_oh),                        &
         [row_length,rows]) - RESHAPE(rk(:,ih2so4_hv)*                         &
          y(:,nn_h2so4),[row_length,rows]))*cdt
      ELSE IF (iso2_oh > 0) THEN
        delh2so4_chem(:,:,k) = delh2so4_chem(:,:,k) +                          &
         RESHAPE(rk(:,iso2_oh)*y(:,nn_so2)*y(:,nn_oh),                         &
         [row_length,rows])*cdt
      END IF

      IF (uph2so4inaer == 1) THEN
         ! Restore H2SO4 tracer as it will be updated in MODE
         ! using delh2so4_chem
        IF (ukca_config%l_fix_ukca_h2so4_ystore) THEN
           ! calculate delh2so4_chem as the difference in H2SO4 over chemistry
           ! zftr is already in VMR, so divide by CDT to give as vmr/s
          delh2so4_chem(:,:,k) = RESHAPE((zftr(:,istore_h2so4) - ystore(:)),   &
                                         [row_length,rows]) / cdt
           ! primary array passed is zftr, so copy back to this, NOT y
          zftr(:,istore_h2so4) = ystore(:)
        ELSE
          y(:,nn_h2so4) = ystore(:)
        END IF
      END IF
    END IF

    ! 3D flux diagnostics
    IF (L_asad_use_chem_diags .AND.                                            &
         ((L_asad_use_flux_rxns .OR. L_asad_use_rxn_rates) .OR.                &
         (L_asad_use_wetdep .OR. L_asad_use_drydep)))                          &
         CALL asad_chemical_diagnostics(row_length,rows,                       &
         model_levels,dpd_dummy,dpw_dummy,prk_dummy,y_dummy,                   &
         ix,jy,k,volume,ierr)

    ! PSC diagnostics
    IF (L_asad_use_chem_diags .AND. L_asad_use_psc_diagnostic)                 &
         CALL asad_psc_diagnostic(row_length,rows,model_levels,                &
         fpsc1_dummy,fpsc2_dummy,ix,jy,k,ierr)

    ! Bring results back from vmr to mmr.
    ! Also bring back results for nontransported RO2 to all_ntp
    jspf = 0 ! reset counter for species in f array
    DO js = 1,jpspec
      DO jtr=1,jpctr
        IF (advt(jtr) == speci(js)) THEN
          jspf = jspf+1
          tracer(:,:,k,jtr) = RESHAPE(zftr(:,jspf),[row_length,rows])          &
                                  *c_species(jtr)
        END IF
      END DO

      ! If RO2 species are not being transported, map RO2 species
      ! back to the all_ntp 3D array
      IF (ukca_config%l_ukca_ro2_ntp) THEN
        DO jro2 = 1, jpro2
          IF (spro2(jro2) == speci(js) .AND. ctype(js) == 'OO') THEN
            jspf = jspf+1
            ! Get index of speci(js) in NTP array,
            ! If not found this is a fatal error.
            l = name2ntpindex(spro2(jro2))

            ! Find location of species in nadvt array
            jna = nlnaro2(jro2)
            all_ntp(l)%data_3d(:,:,k) = RESHAPE(zftr(:,jspf),                  &
                     [row_length,rows]) * c_na_species(jna)
          END IF ! Close IF RO2 species
        END DO   ! Close loop through RO2 species
      END IF     ! Close IF RO2_NTP
    END DO       ! Close loop through all species

    ! Set SS species concentrations for output (stratospheric configurations)

    ! O1D mmr
    IF (O1D_in_ss) THEN
      l = name2ntpindex('O(1D)     ')
      all_ntp(l)%data_3d(:,:,k) =                                              &
        RESHAPE(y(:,nn_o1d)/tnd(:),[row_length,rows]) * c_o1d
    END IF

    ! O3P mmr
    IF (O3P_in_ss) THEN
      l = name2ntpindex('O(3P)     ')
      all_ntp(l)%data_3d(:,:,k) =                                              &
        RESHAPE(y(:,nn_o3p)/tnd(:),[row_length,rows]) * c_o3p
    END IF

    ! First copy the concentrations from the zftr array to the
    ! diag arrays, reshape and convert to moles.
    ! The indices (n_ch4, n_n2o etc.) don't refer to the location
    ! in the zftr array if RO2_NTP is true, but instead to the
    ! location in the tracer array.
    DO jspf = 1, jpcspf

      IF (n_ch4 > 0) THEN
        IF (specf(jspf) == advt(n_ch4)) THEN
          atm_ch4_mol(:,:,k) = RESHAPE(zftr(:,jspf)*tnd(:),                    &
                           [row_length,rows])*volume(:,:,k)*                   &
                           1.0e6/avogadro
        END IF
      END IF

      ! CO
      IF (n_co > 0) THEN
        IF (specf(jspf) == advt(n_co)) THEN
          atm_co_mol(:,:,k) = RESHAPE(zftr(:,jspf)*tnd(:),                     &
             [row_length,rows])*volume(:,:,k)* 1.0e6/avogadro
        END IF
      END IF

      ! N2O
      IF (n_n2o > 0) THEN
        IF (specf(jspf) == advt(n_n2o)) THEN
          atm_n2o_mol(:,:,k) = RESHAPE(zftr(:,jspf)*tnd(:),                    &
            [row_length,rows])*volume(:,:,k)* 1.0e6/avogadro
        END IF
      END IF

      ! CFC-12
      IF (n_cf2cl2 > 0) THEN
        IF (specf(jspf) == advt(n_cf2cl2)) THEN
          atm_cf2cl2_mol(:,:,k) = RESHAPE(zftr(:,jspf)*tnd(:),                 &
            [row_length,rows])*volume(:,:,k)* 1.0e6/avogadro
        END IF
      END IF

      ! CFC-11
      IF (n_cfcl3 > 0) THEN
        IF (specf(jspf) == advt(n_cfcl3)) THEN
          atm_cfcl3_mol(:,:,k) = RESHAPE(zftr(:,jspf)*tnd(:),                  &
            [row_length,rows])*volume(:,:,k)* 1.0e6/avogadro
        END IF
      END IF

      ! CH3Br
      IF (n_mebr > 0) THEN
        IF (specf(jspf) == advt(n_mebr)) THEN
          atm_mebr_mol(:,:,k) = RESHAPE(zftr(:,jspf)*tnd(:),                   &
            [row_length,rows])*volume(:,:,k)* 1.0e6/avogadro
        END IF
      END IF

      ! H2
      IF (n_h2 > 0) THEN
        IF (specf(jspf) == advt(n_h2)) THEN
          atm_h2_mol(:,:,k) = RESHAPE(zftr(:,jspf)*tnd(:),                     &
              [row_length,rows])*volume(:,:,k)* 1.0e6/avogadro
        END IF
      END IF

    END DO ! End loop through species in zftr array

  ELSE  ! Backward Euler with non-families

    !         Calculate total number  density, o2, h2o, and tracer
    !         concentrations for Backward Euler solver

    BE_tnd(:) = zp(:) / (Boltzmann * 1.0e6 * zt(:))
    BE_o2(:)  = 0.2095 * BE_tnd(:)
    BE_h2o(:) = RESHAPE(q(:,:,k),[theta_field_size])*BE_tnd(:)
    BE_vol(:) = RESHAPE(volume(:,:,k),[theta_field_size])*1.0e6 !  m3->cm3

    ! Air density in kg m-3 (note that in this formula rgas is
    ! the gas constant for dry air:  287.05 J kg-1 K-1)
    BE_rho(:) = zp(:) / ( rgas * zt(:) )

    DO l=1,jpdw
      zwetrt2(:,l) = RESHAPE(zwetrt(:,:,k,l),[theta_field_size])
    END DO
    DO l=1,jpdd
      zdryrt2(:,l) = RESHAPE(zdryrt(:,:,l),[theta_field_size])
    END DO

    DO l=1,jpctr
      zftr(:,l) = zftr(:,l) * BE_tnd(:)
    END DO

    ! Update non-advected species (only for B-E solvers) from
    ! non transported prognostics
    IF (.NOT. ALLOCATED(zfnatr)) ALLOCATE(zfnatr(n_pnts,nnaf))
    DO js = 1, nnaf
      ! Get index of nadvt(js) in NTP array.
      ! If not found this is a fatal error.
      l = name2ntpindex(nadvt(js))
      zfnatr(:,js) = RESHAPE(all_ntp(l)%data_3d(:,:,k),[n_pnts]) /             &
                      c_na_species(js)
    END DO


    DO l=1,nnaf
      zfnatr(:,l) = zfnatr(:,l) * BE_tnd(:)
    END DO

    IF (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) THEN
      DO l=1,jpdw
        DO j=1,jpeq+1
          zfrdiss2(:,l,j) = RESHAPE(zfrdiss(:,:,k,l,j),                        &
                            [theta_field_size])
        END DO
      END DO
    END IF

    !         Assign wet and dry deposition rates to species

    CALL ukca_be_wetdep(n_pnts, zwetrt2, be_wetrt)

    CALL ukca_be_drydep(k, n_pnts, nlev_with_ddep2, zdryrt2, be_dryrt)

    ! Assign fractional dissociation
    IF (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) THEN
      DO i=1,jpeq+1
        CALL ukca_be_wetdep(n_pnts, zfrdiss2(:,:,i),                           &
                            be_frdiss(:,:,i))
      END DO
    END IF

    ! collate 1-D pH values to use in BE Solver
    H_plus_1d_arr(:) = RESHAPE(H_plus_3d_arr(:,:,k),[theta_field_size])

    !         Calculate reaction rate coefficients

    IF (.NOT. ALLOCATED(BE_rc))                                                &
         ALLOCATE(BE_rc(theta_field_size,nr_therm))

    IF (ukca_config%l_ukca_raq .OR. ukca_config%l_ukca_raqaero) THEN
      ! Set RH (fraction: 0 - 0.999), aerosol number (m-3) for the sea-salt
      ! modes and aerosol mmr (kg kg-1) for the other CLASSIC aerosol types,
      ! as well as heterogeneous rates (s-1) before the call to UKCA_CHEMCO_RAQ.
      !
      BE_rh_frac      (:) = RESHAPE (rh           (:,:,k), [theta_field_size])
      BE_so4_aitken   (:) = RESHAPE (so4_aitken   (:,:,k), [theta_field_size])
      BE_so4_accum    (:) = RESHAPE (so4_accum    (:,:,k), [theta_field_size])
      BE_soot_fresh   (:) = RESHAPE (soot_fresh   (:,:,k), [theta_field_size])
      BE_soot_aged    (:) = RESHAPE (soot_aged    (:,:,k), [theta_field_size])
      BE_ocff_fresh   (:) = RESHAPE (ocff_fresh   (:,:,k), [theta_field_size])
      BE_ocff_aged    (:) = RESHAPE (ocff_aged    (:,:,k), [theta_field_size])
      BE_biogenic     (:) = RESHAPE (biogenic     (:,:,k), [theta_field_size])
      BE_sea_salt_film(:) = RESHAPE (sea_salt_film(:,:,k), [theta_field_size])
      BE_sea_salt_jet (:) = RESHAPE (sea_salt_jet (:,:,k), [theta_field_size])

      IF (.NOT. ALLOCATED(BE_hrc))                                             &
        ALLOCATE (BE_hrc (theta_field_size, jphk))

      ! Fill in stratospheric flag indicator, which will be passed
      ! to UKCA_CHEMCO_RAQ in order to set tropospheric heterogeneous
      ! reactions (if present) to zero in the stratosphere
      stratflag (:) = (.NOT. RESHAPE ( L_troposphere (:,:,k),                  &
                                      [theta_field_size]) )

      CALL ukca_chemco_raq(nr_therm, n_pnts, zt(1:n_pnts),                     &
                      BE_tnd(1:n_pnts), BE_h2o(1:n_pnts),                      &
                      BE_o2(1:n_pnts),  zclw(1:n_pnts),                        &
                      zfcloud(1:n_pnts), be_frdiss(1:n_pnts,:,:),              &
                      BE_rho(1:n_pnts),                                        &
                      BE_rh_frac       (1:n_pnts),                             &
                      BE_so4_aitken    (1:n_pnts), BE_so4_accum    (1:n_pnts), &
                      BE_soot_fresh    (1:n_pnts), BE_soot_aged    (1:n_pnts), &
                      BE_ocff_fresh    (1:n_pnts), BE_ocff_aged    (1:n_pnts), &
                      BE_biogenic      (1:n_pnts),                             &
                      BE_sea_salt_film (1:n_pnts), BE_sea_salt_jet (1:n_pnts), &
                      stratflag        (1:n_pnts),                             &
                      BE_rc (1:n_pnts,:), BE_hrc (1:n_pnts,:),                 &
                      H_plus_1d_arr(1:n_pnts))
    ELSE

      CALL ukca_chemco(nr_therm, n_pnts, zt(1:n_pnts),                         &
                       BE_tnd(1:n_pnts), BE_h2o(1:n_pnts),                     &
                       BE_o2(1:n_pnts), zclw(1:n_pnts),                        &
                       zfcloud(1:n_pnts), BE_frdiss(1:n_pnts,:,:),             &
                       k_dms(1:n_pnts,:), BE_rc(1:n_pnts,:),                   &
                       H_plus_1d_arr(1:n_pnts))
    END IF

    !         Assign tracer concentrations to species concentrations

    BE_y(1:n_pnts,1:jpspec) = 0.0
    DO i = 1,jpspec
      j_loop1: DO j = 1,jpctr
        IF (speci(i) == advt(j)) THEN
          BE_y(1:n_pnts,i) = zftr(1:n_pnts,j)
          EXIT j_loop1
        END IF
      END DO j_loop1
    END DO

    ! Assign non-advected concentrations to species concentrations
    DO i = 1,jpspec
      j_loop2: DO j = 1,nnaf
        IF (speci(i) == nadvt(j)) THEN
          BE_y(1:n_pnts,i) = zfnatr(1:n_pnts,j)
          EXIT j_loop2
        END IF
      END DO j_loop2
    END DO

    IF (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) THEN
      SO2_wetox_H2O2(:) = 0.0
      SO2_dryox_OH(:)   = 0.0
      SO2_wetox_O3(:)   = 0.0
    ELSE
      delSO2_wet_H2O2(:,:,k) = 0.0
      delSO2_wet_O3(:,:,k)   = 0.0
      delh2so4_chem(:,:,k)   = 0.0
      delSO2_drydep(:,:,k)   = 0.0
      delSO2_wetdep(:,:,k)   = 0.0
    END IF

    !         Call Backward Euler solver
    !         N.B. Emissions already added, via call to TR_MIX from
    !         UKCA_EMISSION_CTL

    IF (ukca_config%l_ukca_aerchem) THEN

      CALL ukca_deriv_aero(nr_therm, n_be_calls, n_pnts,                       &
                    BE_rc(1:n_pnts,:),                                         &
                    BE_wetrt(1:n_pnts,:), BE_dryrt(1:n_pnts,:),                &
                    zprt1d(1:n_pnts,:), k_dms(1:n_pnts,:),                     &
                    BE_h2o(1:n_pnts), BE_tnd(1:n_pnts),                        &
                    BE_o2(1:n_pnts),                                           &
                    dts, BE_y(1:n_pnts,:),                                     &
                    SO2_wetox_H2O2(1:n_pnts),                                  &
                    SO2_wetox_O3(1:n_pnts),                                    &
                    SO2_dryox_OH(1:n_pnts) )
    ELSE IF (ukca_config%l_ukca_raq) THEN

      CALL ukca_deriv_raq(nr_therm, n_be_calls,                                &
                 n_pnts, dts,                                                  &
                 BE_h2o(1:n_pnts),                                             &
                 BE_vol(1:n_pnts),                                             &
                 BE_rc(1:n_pnts,:), BE_hrc(1:n_pnts,:),                        &
                 BE_dryrt(1:n_pnts,:),                                         &
                 BE_wetrt(1:n_pnts,:), zprt1d(1:n_pnts,:),                     &
                 ldepd, ldepw,                                                 &
                 BE_y(1:n_pnts,:), ddflux(1:n_pnts,:),                         &
                 wdflux(1:n_pnts,:))

      DO i=1,jpdd
        dry_dep_3d(:,:,k,i) = RESHAPE(ddflux(:,i),                             &
                             [row_length,rows])
      END DO
      DO i=1,jpdw
        wet_dep_3d(:,:,k,i) = RESHAPE(wdflux(:,i),                             &
                             [row_length,rows])
      END DO
    ELSE IF (ukca_config%l_ukca_raqaero) THEN

      ! Store H2SO4 tracer if it will be updated in MODE using delh2so4_chem
      IF (uph2so4inaer == 1) ystore(:) = BE_y(:,nn_h2so4)

      CALL ukca_deriv_raqaero(nr_therm, n_be_calls,                            &
                 n_pnts, dts,                                                  &
                 be_rc(1:n_pnts,:), be_wetrt(1:n_pnts,:),                      &
                 be_dryrt(1:n_pnts,:), zprt1d(1:n_pnts,:),                     &
                 be_h2o(1:n_pnts), BE_y(1:n_pnts,:),                           &
                 so2_wetox_H2O2(1:n_pnts),                                     &
                 so2_wetox_O3(1:n_pnts),                                       &
                 so2_dryox_OH(1:n_pnts) )

      ! Restore H2SO4 tracer as it will be updated in MODE using delh2so4_chem
      IF (uph2so4inaer == 1) BE_y(:,nn_h2so4) = ystore(:)

    ELSE

      CALL ukca_deriv(nr, n_be_calls, n_pnts,                                  &
                 BE_rc(1:n_pnts,:), BE_wetrt(1:n_pnts,:),                      &
                 BE_dryrt(1:n_pnts,:), zprt1d(1:n_pnts,:),                     &
                 BE_h2o(1:n_pnts), BE_tnd(1:n_pnts),                           &
                 BE_o2(1:n_pnts),                                              &
                 dts, BE_y(1:n_pnts,:) )
    END IF
    IF (ALLOCATED(BE_hrc)) DEALLOCATE(BE_hrc)
    IF (ALLOCATED(BE_rc))  DEALLOCATE(BE_rc)

    IF (k == model_levels .OR. k == model_levels-1 .OR.                        &
        k == model_levels - 2) THEN
      CALL ukca_ch4_stratloss(n_be_calls, n_pnts,                              &
                 BE_vol(1:n_pnts), dts,                                        &
                 BE_y(1:n_pnts,nn_ch4), strat_ch4loss_2d(1:n_pnts,k))
    END IF

    DO j = 1,jpctr
      i_loop1: DO i = 1,jpspec
        IF (advt(j) == speci(i)) THEN
          zftr(1:n_pnts,j) = BE_y(1:n_pnts,i)
          EXIT i_loop1
        END IF
      END DO i_loop1
    END DO

    DO js = 1, jpctr
      tracer(:,:,k,js) = RESHAPE(zftr(:,js)/BE_tnd(:),                         &
                             [row_length,rows])*c_species(js)
    END DO

    ! Convert non-advected tracers back to vmr
    DO j = 1,nnaf
      i_loop2: DO i = 1,jpspec
        IF (nadvt(j) == speci(i)) THEN
          zfnatr(:,j) = BE_y(1:n_pnts,i)/BE_tnd(1:n_pnts)
          EXIT i_loop2
        END IF
      END DO i_loop2
    END DO

    ! Set the value of all non-transported prognostics here
    DO js = 1, nnaf
      ! get index of nadvt(js) in NTP array
      l = name2ntpindex(nadvt(js))
      ! set data in appropriate entry in all_ntp
      all_ntp(l)%data_3d(:,:,k) = RESHAPE(zfnatr(:,js),                        &
        [row_length,rows])*c_na_species(js)
    END DO

    ! First copy the concentrations from the zftr array to the
    ! diag arrays, reshape and convert to moles.
    ! The values in the troposphere/stratosphere are masked
    ! off below

    ! CH4 burden in moles. Copy tropospheric values to stratospheric array
    ! for later masking using tropospheric mask.
    trop_ch4_mol(:,:,k) = RESHAPE(zftr(:,n_ch4),                               &
                          [row_length,rows])*volume(:,:,k)*                    &
                            1.0e6/avogadro
    strat_ch4_mol(:,:,k) = trop_ch4_mol(:,:,k)

    ! O3 burden in moles
    trop_o3_mol(:,:,k) = RESHAPE(zftr(:,n_o3),                                 &
                          [row_length,rows])*volume(:,:,k)*                    &
                            1.0e6/avogadro

    ! OH burden in moles
    trop_oh_mol(:,:,k) = RESHAPE(BE_y(:,nn_oh),                                &
                          [row_length,rows])*volume(:,:,k)*                    &
                            1.0e6/avogadro

    ! Stratospheric CH4 loss rate
    strat_ch4loss(:,:,k) = RESHAPE(strat_ch4loss_2d(:,k),                      &
                          [row_length,rows])

    IF (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) THEN
      delSO2_wet_H2O2(:,:,k)=RESHAPE(SO2_wetox_H2O2(:),                        &
                             [row_length,rows])
      delSO2_wet_O3(:,:,k)=RESHAPE(SO2_wetox_O3(:),                            &
                             [row_length,rows])
      delh2so4_chem(:,:,k)=RESHAPE(SO2_dryox_OH(:),                            &
                             [row_length,rows])
      delSO2_drydep(:,:,k) = RESHAPE(BE_dryrt(:,nn_so2)*                       &
                  BE_y(:,nn_so2),[row_length,rows]) * dts
      delSO2_wetdep(:,:,k) = RESHAPE(BE_wetrt(:,nn_so2)*                       &
                  BE_y(:,nn_so2),[row_length,rows]) * dts
    END IF

  END IF ! End of BE solver section
END DO ! level loop (k)
!$OMP END DO

IF (ALLOCATED(zfnatr)) DEALLOCATE(zfnatr)
IF (ALLOCATED(ystore)) DEALLOCATE(ystore)

!$OMP END PARALLEL

! Now mask off stratospheric and tropospheric diagnostics
WHERE (l_troposphere(:,:,:))
  strat_ch4_mol(:,:,:) = 0.0
ELSE WHERE
  trop_ch4_mol(:,:,:) = 0.0
  trop_o3_mol(:,:,:) = 0.0
  trop_oh_mol(:,:,:) = 0.0
END WHERE

! Call raq diagnostic routine
IF (ukca_config%l_enable_diag_um .AND. ukca_config%l_ukca_raq) THEN
  CALL ukca_raq_diags (row_length, rows, model_levels, ntracers,               &
       dry_dep_3d, wet_dep_3d, tracer, pres, temp,                             &
       len_stashwork, stashwork)
END IF

! Rescale bromine and chlorine tracers to guarantee conservation of total
! chlorine, bromine, and hydrogen over timestep. Only makes sense if at least
! chlorine chemistry is present.

IF (nn_cl > 0) THEN
  CALL ukca_conserve(row_length, rows, model_levels, ntracers,                 &
       tracer, pres, drain, crain, after_chem)
END IF

IF (ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_stratcfc .OR.             &
    ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat) THEN

  IF (ukca_config%l_ukca_het_psc) THEN
    ! Do NAT PSC sedimentation

    ! take NAT out of gasphase again
    tracer(:,:,:,n_hono2) = tracer(:,:,:,n_hono2) - shno3_3d

    CALL ukca_sediment(rows, row_length, model_levels, shno3_3d,               &
             qcf, r_theta_levels, mass, secs_per_step,                         &
             L_troposphere(:,:,:) )

    ! add solid-phase HNO3 back to gasphase HNO3
    tracer(:,:,:,n_hono2) = tracer(:,:,:,n_hono2) + shno3_3d
  END IF

  ! i_ukca_topboundary==0 (i_top_none) corresponds to no overwriting of
  ! top level(s) or any top boundary condition
  IF (ukca_config%i_ukca_topboundary == i_top_2levH2O) THEN
    ! Tracer overwrites required to stop accumulation of tracer mass
    ! in the uppermost layers.  Exclude water vapour.
    DO i=1,rows
      DO j=1,row_length
        DO k=1,ntracers
          IF (k /= n_h2o) THEN
            tracer(j,i,model_levels  ,k) = tracer(j,i,model_levels-2,k)
            tracer(j,i,model_levels-1,k) = tracer(j,i,model_levels-2,k)
          END IF
        END DO
      END DO
    END DO
  ELSE IF (ukca_config%i_ukca_topboundary == i_top_1lev) THEN
    ! over-write top level for all tracers with 2nd-highest level
    DO i=1,rows
      DO j=1,row_length
        DO k=1,ntracers
          tracer(j,i,model_levels,k) = tracer(j,i,model_levels-1,k)
        END DO
      END DO
    END DO
  ELSE IF (ukca_config%i_ukca_topboundary >= i_top_BC) THEN
    ! Apply top boundary condition for NO, CO, O3, optionally H2O, using
    ! ACE-FTS climatologies (assumes constant latitude on each row)
    ! no tracer over-writing at top levels
    IF ((z_top_of_model > 85500.0) .OR. (z_top_of_model < 79000.0)) THEN
      cmessage='Can only impose a top boundary condition at 85 km.'
      errcode=25
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
    IF (.NOT. ALL(ABS(latitude(row_length,:) - latitude(1,:))                  &
                                                          < EPSILON(0.0))) THEN
      cmessage=                                                                &
        'Can only impose top boundary condition if latitude is constant on rows'
      errcode=1
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
    CALL ukca_top_boundary(row_length, rows, model_levels, ntracers,           &
                           latitude(1,:), tracer)
  END IF

  ! Copy NAT MMR into user_diagostics
  IF (ukca_config%l_ukca_het_psc) THEN
    nat_psc(:,:,:)=shno3_3d(:,:,:)
  END IF

ELSE IF (.NOT. ukca_config%l_ukca_offline) THEN   ! tropospheric chemistry
  ! Call routine to overwrite O3 and HNO3 species once per day
  ! above tropopause. Only for tropospheric chemistry

  CALL ukca_stratf(row_length,rows, model_levels,                              &
                   jpctr, env_ozone3d,                                         &
                   tracer(1:row_length,1:rows,1:model_levels,1:jpctr))

END IF     ! l_ukca_strat etc


IF (ukca_config%l_ukca_het_psc .AND. ALLOCATED(shno3_3d)) DEALLOCATE(shno3_3d)

IF (firstcall) firstcall = .FALSE.


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chemistry_ctl

END MODULE ukca_chemistry_ctl_mod
