! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!     UKCA-MODE aerosol code: interface routine called from
!     UKCA_MAIN1 to perform a 1-timestep integration.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE ukca_aero_ctl_mod

USE ukca_um_legacy_mod, ONLY: stashcode_glomap_sec, pi,                        &
                              n_mode_diags,                                    &
                              nukca_d1items,                                   &
                              ukcaD1codes,                                     &
                              item1_mode_diags, L_ukca_mode_diags,             &
                              item1_nitrate_diags, itemN_nitrate_diags,        &
                              log_v, rgas => r
USE ukca_types_mod, ONLY: log_small, integer_32
USE parkind1,   ONLY: jprb, jpim
USE yomhook,    ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_AERO_CTL_MOD'

!-----------------------------------------------------------------------
! Derived types
!-----------------------------------------------------------------------

TYPE :: segment_data_type

  INTEGER           :: nbox_this_seg
  ! Size of nbox_s for this segment
  INTEGER           :: nchemg
  ! Number of gas phase chemistry tracers
  INTEGER           :: nhet
  ! Number of heterogeneous reaction rates
  INTEGER           :: nbudaer
  ! Number of aerosol budget fields
  INTEGER           :: nadvg
  ! Number of advected gas phase tracers

  REAL, POINTER ::delso2  (:)
  !  S(IV) --> S(VI) by H2O2 (molecules per cc) [input if WETOX_IN_AER=0]
  REAL, POINTER ::delso2_2(:)
  !  S(IV) --> S(VI) by O3   (molecules per cc) [input if WETOX_IN_AER=0]

  REAL, POINTER :: fac(:)
  ! Conversion factor for budget diagnostics
  REAL, POINTER :: fac_mmrconv(:)
  ! Conversion factor to convert MMRSO2/s to molecules/cc/DTC

  REAL, POINTER :: height(:)
  ! Mid-level height of gridbox
  REAL, POINTER :: htpblg(:)
  ! Height of boundary-layer in gridbox vertical-column

  REAL, POINTER :: mdtfixflag(:,:)
  REAL, POINTER :: mdtfixsink(:,:)
  ! Arrays storing info about how much mass is removed when ND>0 for o-o-r MDT
  ! mdtfixflag gets to set = 100.0 when fix applied so that means show the
  !         percentage of timesteps on which the fix is applied.
  ! mdtfixsink stores the amount of *total mass" (over all components) that is
  !         removed when the fix is applied, to support mass budget
  !         calculation. Units are moles/s.

  REAL, POINTER :: het_rates(:,:)
  ! Diagnostic to hold heterogeneous rates for tropospheric chemistry

  LOGICAL (KIND=log_small), POINTER :: mask1(:)
  INTEGER (KIND=integer_32), POINTER :: nbadmdt(:,:)
  INTEGER (KIND=integer_32), POINTER :: n_merge_1d(:,:)


  REAL, POINTER :: sarea(:,:)
  ! Surface area concentration for each mode (cm^2 / cm^3)
  REAL, POINTER :: vconc(:,:)
  ! Volume concentration for each mode

  REAL, POINTER ::aird(:)
  ! Grid box mass of air (kg)
  REAL, POINTER ::airdm3(:)
  ! Number density of air (per cm3)

  REAL, POINTER :: autoconv1d(:)
  ! Autoconversion rate (including accretion) together
  ! with snow and ice melting rate (kg.kg^-1.s^-1)

  REAL , POINTER :: bud_aer_mas(:,:)
  ! Outputs for budget calculations

  REAL, POINTER :: cn_3nm(:)
  ! CN concentration (dry diameter > 3nm)
  REAL, POINTER :: ccn_1(:)
  ! CCN concentration (acc-sol + cor-sol)
  REAL, POINTER :: ccn_2(:)
  ! CCN concentration (acc-sol + cor-sol + Aitsol>25nm drydp)
  REAL, POINTER :: ccn_3(:)
  ! CCN concentration (acc-sol + cor-sol + Aitsol>35nm drydp)
  REAL, POINTER :: cdn(:)
  ! CDN concentration

  REAL, POINTER :: erfterm(:)
  REAL, POINTER :: erf_arg(:)           ! Error Fn argument
  ! Items for CN,CCN,CDN calculation

  INTEGER, POINTER  :: lday(:)
  ! Switch for day/night (1/0)
  INTEGER, POINTER  :: jlabove(:)
  ! Index of box directly above this grid box
  INTEGER, POINTER  :: ilscat(:)
  ! Land surface category (based on 9 UM landsurf types)
  INTEGER, POINTER  :: iarr(:)
  ! Index of LAT for grid box
  INTEGER, POINTER  :: karr(:)
  ! Index of LON for grid box
  INTEGER, POINTER  :: larr(:)
  ! Index of vertical level for grid box
  REAL, POINTER :: nd(:,:)
  ! Aerosol ptcl number density for mode (cm^-3)
  REAL, POINTER :: mdt(:,:)
  ! Avg tot mass of aerosol ptcl in mode (particle^-1)
  REAL, POINTER :: mdwat(:,:)
  ! Molecular concentration of water (molecules per particle)
  REAL, POINTER :: drydp(:,:)
  ! Geometric mean dry diameter of particles in each mode (m)
  REAL, POINTER :: wetdp(:,:)
  ! Geometric mean wet diameter of particles in each mode (m)
  REAL, POINTER :: rhopar(:,:)
  ! Total particle density [incl. H2O & insoluble cpts] (kgm^-3)
  REAL, POINTER :: dvol(:,:)
  ! Geometric mean dry volume of particles in each mode (m^3)
  REAL, POINTER :: wvol(:,:)
  ! Geometric mean wet volume of particles in each mode (m^3)
  REAL, POINTER :: md(:,:,:)
  ! Avg cpt mass of aerosol particle in mode (particle^-1)
  REAL, POINTER :: s0(:,:)
  ! Partial masses of gas phase species (kg per gridbox)
  REAL, POINTER :: s0_dot_condensable(:,:)
  ! ASAD tendencies for condensable gas phase species (vmr per s)
  REAL, POINTER :: sm(:)     ! Mass of air in gridbox (kg)
  ! Number density of air (per m3)
  REAL, POINTER ::rhoa(:)
  ! Air density (kg/m3)
  REAL, POINTER :: vba(:)
  ! Mean free speed of air molecules (m/s)
  REAL, POINTER :: tsqrt(:)
  ! Square-root of centre level temperature (K)
  REAL, POINTER :: dvisc(:)
  ! Dynamic viscosity of air (kg m^-1 s^-1)
  REAL, POINTER :: mfpa(:)
  ! Mean free path of air (m)
  REAL, POINTER :: t(:)
  ! Air temperature at mid-point (K)
  REAL, POINTER :: rh(:)
  ! Relative humidity (fraction) - gridbox mean
  REAL, POINTER :: rh_clr(:)
  ! Relative humidity (fraction) - clear sky portion
  REAL, POINTER :: s(:)
  ! Specific humidity (kg/kg)
  REAL, POINTER :: pmid(:)
  ! Air pressure at mid-point (Pa)
  REAL, POINTER :: pupper(:)
  ! Air pressure at upper interface (Pa)
  REAL, POINTER :: plower(:)
  ! Air pressure at lower interface (Pa)
  REAL, POINTER :: zo3(:)
  ! Background vmr of O3 (dimensionless)
  REAL, POINTER ::zho2(:)
  ! Background conc. of HO2 (molecules per cc)
  REAL, POINTER ::zh2o2(:)
  ! Background conc. of H2O2 (molecules per cc)
  REAL, POINTER :: ustr(:)
  ! Surface friction velocity (m/s)
  REAL, POINTER :: znotg(:)
  ! Roughness length (m)
  REAL, POINTER :: surtp(:)
  ! Surface type: 0=seasurf,1=landsurf,2=above-seasurf,3=above-landsurf
  REAL, POINTER :: land_frac(:)
  ! Fraction of horizontal gridbox area covered by land
  REAL, POINTER :: seaice(:)
  ! Fraction of horizontal gridbox area containing seaice
  REAL, POINTER :: craing(:)
  ! Rain rate for conv precip. in box (kgm^-2s^-1)
  REAL, POINTER :: draing(:)
  ! Rain rate for dyn. precip. in box (kgm^-2s^-1)
  REAL, POINTER :: csnowg(:)
  ! Rain rate for conv snow. in box (kgm^-2s^-1)
  REAL, POINTER :: dsnowg(:)
  ! Rain rate for dyn. precip. in box (kgm^-2s^-1)
  REAL, POINTER :: craing_up(:)
  ! Rain rate for conv precip. in box above (kgm^-2s^-1)
  REAL, POINTER ::fconv_conv(:)
  ! Fraction of box condensate --> rain in 6 hours (conv)
  REAL, POINTER :: lowcloud(:)
  ! Horizontal low cloud fraction
  REAL, POINTER :: vfac(:)
  ! Vertical low cloud fraction
  REAL, POINTER :: lwc(:)
  ! Cloud liquid water content [kg/m3]
  REAL, POINTER :: clwc(:)
  ! Cloud liquid water content [kg/kg]
  REAL, POINTER :: clf(:)
  ! Liquid Cloud fraction
  REAL, POINTER :: pvol(:,:,:)
  ! Aerosol partial volume of each cpt in each mode
  REAL, POINTER :: pvol_wat(:,:)
  ! Aerosol partial volume of water in each mode
  REAL, POINTER :: tr_rs(:)
  ! Local variable to hold re-shaped aerosol tracers

  REAL, POINTER :: v1d_tmp(:) ! local temporary staging for multiplications

END TYPE segment_data_type

CONTAINS

! Subroutine Interface:
SUBROUTINE ukca_aero_ctl(i_month, i_day_number,                                &
                i_hour,i_minute, dtc,                                          &
                row_length, rows, model_levels,                                &
                n_chemistry_tracers,                                           &
                n_mode_tracers,                                                &
                pres,                                                          &
                temp,                                                          &
                q,                                                             &
                rh3d, rh3d_clr,                                                &
                p_bdrs,                                                        &
                chemistry_tracer_names,                                        &
                mode_tracer_names,                                             &
                chemistry_tracers,                                             &
                mode_tracers,                                                  &
                sea_ice_frac,                                                  &
                z0m,                                                           &
                u_s,                                                           &
                drain, crain,                                                  &
                dsnow, csnow,                                                  &
                autoconv, accretion,                                           &
                rim_agg, rim_cry,                                              &
                land_fraction,                                                 &
                nbox,                                                          &
                delso2_wet_h2o2,                                               &
                delso2_wet_o3,                                                 &
                delso2_dry_oh,                                                 &
                mode_diags,                                                    &
                cloud_frac,                                                    &
                cloud_liq_frac,                                                &
                cloud_liq_wat,                                                 &
                z_half_alllevs,                                                &
                mass,zbl,                                                      &
                dryox_in_aer,                                                  &
                wetox_in_aer,                                                  &
                all_ntp,                                                       &
                nseg, nbox_s, ncol_s, lbase, stride_s                          &
                )

USE ukca_drydiam_field_mod, ONLY: drydiam
USE ukca_config_defs_mod, ONLY: nmax_mode_diags

USE ukca_mode_setup,  ONLY: nmodes,                                            &
                            mode_nuc_sol,mode_ait_sol,                         &
                            mode_acc_sol, mode_cor_sol,                        &
                            mode_ait_insol, mode_acc_insol,                    &
                            mode_cor_insol, cp_su, cp_bc,                      &
                            cp_oc, cp_so, cp_cl, cp_du,                        &
                            cp_no3, cp_nh4, cp_nn

USE ukca_mode_tracer_maps_mod, ONLY:  nmr_index, mmr_index
USE ukca_mode_verbose_mod, ONLY: verbose => glob_verbose
USE chemistry_constants_mod, ONLY: avogadro, boltzmann
USE ukca_mode_check_artefacts_mod, ONLY:                                       &
                            ukca_mode_check_artefacts
USE ukca_cspecies,    ONLY: n_h2so4,n_h2o2,n_so2,n_o3,n_sec_org, n_sec_org_i
USE ukca_mode_diags_mod,  ONLY: l_ukca_cmip6_diags,                            &
                            l_ukca_pm_diags, mdwat_diag,                       &
                            wetdp_diag
USE ukca_impc_scav_mod, ONLY: ukca_mode_imscavcoff
USE ukca_impc_scav_dust_mod, ONLY: ukca_impc_scav_dust_init,                   &
                                   ukca_impc_scav_dust_dealloc
USE asad_mod, ONLY: jpctr
USE ukca_config_specification_mod, ONLY:                                       &
                            ukca_config, glomap_config,                        &
                            i_ukca_activation_arg,                             &
                            i_ukca_activation_jones,                           &
                            glomap_variables


USE ukca_ntp_mod,     ONLY: ntp_type, dim_ntp, name2ntpindex

USE ukca_setup_indices, ONLY: ntraer, nbudaer, mh2o2f, mh2so4,                 &
                              mm_gas, mox, msec_org, msec_orgi, msotwo,        &
                              nadvg, nchemg,                                   &
                              nmasagedsuintr52, nmasagedbcintr52,              &
                              nmasagedocintr52, nmasagedsointr52,              &
                              nmasclprsuaitsol1,nmasclprsuaccsol1,             &
                              nmasclprsucorsol1,nmasclprsuaitsol2,             &
                              nmasclprsuaccsol2,nmasclprsucorsol2,             &
                              nmascoagsuintr12,                                &
                              nmascoagsuintr13, nmascoagsuintr14,              &
                              nmascoagsuintr15, nmascoagsuintr16,              &
                              nmascoagsuintr17, nmascoagocintr12,              &
                              nmascoagocintr13, nmascoagocintr14,              &
                              nmascoagocintr15, nmascoagocintr16,              &
                              nmascoagocintr17, nmascoagsointr12,              &
                              nmascoagsointr13, nmascoagsointr14,              &
                              nmascoagsointr15, nmascoagsointr16,              &
                              nmascoagsointr17, nmascoagsuintr23,              &
                              nmascoagsuintr24, nmascoagbcintr23,              &
                              nmascoagbcintr24, nmascoagocintr23,              &
                              nmascoagocintr24, nmascoagsointr23,              &
                              nmascoagsointr24, nmascoagsuintr34,              &
                              nmascoagbcintr34, nmascoagocintr34,              &
                              nmascoagssintr34, nmascoagsointr34,              &
                              nmascoagduintr34, nmascoagbcintr53,              &
                              nmascoagocintr53, nmascoagbcintr54,              &
                              nmascoagocintr54, nmascoagduintr64,              &
                              nmascondsunucsol,                                &
                              nmascondsuaitsol, nmascondsuaccsol,              &
                              nmascondsucorsol, nmascondsuaitins,              &
                              nmascondsuaccins, nmascondsucorins,              &
                              nmascondocnucsol, nmascondocaitsol,              &
                              nmascondocaccsol, nmascondoccorsol,              &
                              nmascondocaitins, nmascondocaccins,              &
                              nmascondoccorins, nmascondocinucsol,             &
                              nmascondociaitsol, nmascondociaccsol,            &
                              nmascondocicorsol, nmascondociaitins,            &
                              nmascondociaccins, nmascondocicorins,            &
                              nmascondsonucsol,                                &
                              nmascondsoaitsol, nmascondsoaccsol,              &
                              nmascondsocorsol, nmascondsoaitins,              &
                              nmascondsoaccins, nmascondsocorins,              &
                              nmasddepsunucsol, nmasddepsuaitsol,              &
                              nmasddepsuaccsol, nmasddepsucorsol,              &
                              nmasddepssaccsol, nmasddepsscorsol,              &
                              nmasddepbcaitsol, nmasddepbcaccsol,              &
                              nmasddepbccorsol, nmasddepbcaitins,              &
                              nmasddepocnucsol, nmasddepocaitsol,              &
                              nmasddepocaccsol, nmasddepoccorsol,              &
                              nmasddepocaitins, nmasddepsonucsol,              &
                              nmasddepsoaitsol, nmasddepsoaccsol,              &
                              nmasddepsocorsol, nmasddepduaccsol,              &
                              nmasddepducorsol, nmasddepduaccins,              &
                              nmasddepducorins,                                &
                              nmasimscssaccsol, nmasimscocaitsol,              &
                              nmasimscsscorsol, nmasimscbcaitsol,              &
                              nmasimscbcaccsol, nmasimscsunucsol,              &
                              nmasimscsuaitsol, nmasimscsuaccsol,              &
                              nmasimscsucorsol, nmasimscbccorsol,              &
                              nmasimscbcaitins, nmasimscocnucsol,              &
                              nmasimscocaccsol, nmasimscoccorsol,              &
                              nmasimscocaitins, nmasimscsonucsol,              &
                              nmasimscsoaitsol, nmasimscsoaccsol,              &
                              nmasimscsocorsol, nmasimscduaccsol,              &
                              nmasimscducorsol, nmasimscduaccins,              &
                              nmasimscducorins,                                &
                              nmasmergsuintr12, nmasmergocintr12,              &
                              nmasmergsointr12, nmasmergsuintr23,              &
                              nmasmergbcintr23, nmasmergocintr23,              &
                              nmasmergsointr23, nmasmergsuintr34,              &
                              nmasmergbcintr34, nmasmergocintr34,              &
                              nmasmergssintr34, nmasmergduintr34,              &
                              nmasmergsointr34,                                &
                              nmasnuclsunucsol,                                &
                              nmasnuscsunucsol, nmasnuscsuaitsol,              &
                              nmasnuscsuaccsol, nmasnuscsucorsol,              &
                              nmasnuscssaccsol, nmasnuscsscorsol,              &
                              nmasnuscbcaitsol, nmasnuscbcaccsol,              &
                              nmasnuscbccorsol, nmasnuscbcaitins,              &
                              nmasnuscocnucsol, nmasnuscocaitsol,              &
                              nmasnuscocaccsol, nmasnuscoccorsol,              &
                              nmasnuscocaitins,                                &
                              nmasnuscsonucsol, nmasnuscsoaitsol,              &
                              nmasnuscsoaccsol, nmasnuscsocorsol,              &
                              nmasnuscduaccsol, nmasnuscducorsol,              &
                              nmasnuscduaccins, nmasnuscducorins,              &
                              nmasprocsuintr23, nmasprocbcintr23,              &
                              nmasprococintr23, nmasprocsointr23,              &
                              ! Nitrate indices
                              nmasddepntaitsol, nmasddepntaccsol,              &
                              nmasddepntcorsol, nmasddepnhaitsol,              &
                              nmasddepnhaccsol, nmasddepnhcorsol,              &
                              nmasddepnnaccsol, nmasddepnncorsol,              &
                              nmasnuscntaitsol, nmasnuscntaccsol,              &
                              nmasnuscntcorsol, nmasnuscnhaitsol,              &
                              nmasnuscnhaccsol, nmasnuscnhcorsol,              &
                              nmasnuscnnaccsol, nmasnuscnncorsol,              &
                              nmasimscntaitsol, nmasimscntaccsol,              &
                              nmasimscntcorsol, nmasimscnhaitsol,              &
                              nmasimscnhaccsol, nmasimscnhcorsol,              &
                              nmasimscnnaccsol, nmasimscnncorsol,              &
                              nmascoagntintr23, nmascoagnhintr23,              &
                              nmascoagntintr24, nmascoagnhintr24,              &
                              nmascoagntintr34, nmascoagnhintr34,              &
                              nmascoagnnintr34,                                &
                              nmasmergntintr23, nmasmergnhintr23,              &
                              nmasmergntintr34, nmasmergnhintr34,              &
                              nmasmergnnintr34,                                &
                              nmasprocntintr23, nmasprocnhintr23

USE ukca_trop_hetchem_mod, ONLY: ukca_trop_hetchem,                            &
                                 ihet_n2o5, ihet_ho2_ho2

USE ukca_fieldname_mod, ONLY: maxlen_fieldname
USE ukca_environment_fields_mod, ONLY: lscat_zhang

USE ereport_mod,      ONLY: ereport
USE missing_data_mod,      ONLY: rmdi, imdi

USE umErf_mod, ONLY: umErf

USE errormessagelength_mod, ONLY: errormessagelength

USE ukca_calc_drydiam_mod, ONLY: ukca_calc_drydiam
USE ukca_aero_step_mod, ONLY: ukca_aero_step
USE ukca_check_radaer_coupling_mod, ONLY: ukca_check_radaer_coupling
USE umPrintMgr, ONLY: umMessage, umPrint

USE ukca_cdnc_jones_mod, ONLY: ukca_cdnc_jones
USE ukca_ddepaer_coeff_mod, ONLY: set_ddepaer_coeff

!$ USE omp_lib, ONLY: omp_get_thread_num, omp_get_num_threads

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: i_month           ! month
INTEGER, INTENT(IN) :: i_day_number      ! day
INTEGER, INTENT(IN) :: i_hour            ! hour
INTEGER, INTENT(IN) :: i_minute          ! minute
INTEGER, INTENT(IN) :: row_length        ! # of pts in a patch row
INTEGER, INTENT(IN) :: rows              ! # of rows in patch
INTEGER, INTENT(IN) :: model_levels      ! # of model levels
INTEGER, INTENT(IN) :: n_chemistry_tracers ! # of chemistry tracers
INTEGER, INTENT(IN) :: n_mode_tracers    ! # of mode tracers
! nbox means the max. num. of boxes in a chunk and allows for possible variation
INTEGER, INTENT(IN) :: nbox              ! dimension of slice

! switch for updating of condensables in MODE or in UKCA-CHEMISTRY
INTEGER, INTENT(IN) :: dryox_in_aer      ! 0 external, 1 internal
! switch for doing aqueous SO4 production in MODE or in UKCA-CHEMISTRY
INTEGER, INTENT(IN) :: wetox_in_aer      ! 0 external, 1 internal

REAL, INTENT(IN) :: dtc                                ! timestep(s)
REAL, INTENT(IN) :: pres(row_length,rows,model_levels) ! pressure
REAL, INTENT(IN) :: temp(row_length,rows,model_levels) ! temperature
REAL, INTENT(IN) :: q(row_length,rows,model_levels)    ! sp humidity
REAL, INTENT(IN) :: rh3d(row_length,rows,model_levels) ! rh (frac)
REAL, INTENT(IN) :: rh3d_clr(row_length,rows,model_levels)
! rh (frac) - clear sky portion
REAL, INTENT(IN) :: p_bdrs(row_length,rows,0:model_levels)
! pressure on interfaces
REAL, INTENT(IN) :: sea_ice_frac(row_length, rows)     ! sea ice
REAL, INTENT(IN) :: u_s(row_length, rows)              ! friction velocity
REAL, INTENT(IN) :: z0m(row_length, rows)              ! roughness length
REAL, INTENT(IN) :: drain(row_length,rows, model_levels) ! 3-D LS rain rate
REAL, INTENT(IN) :: crain(row_length,rows, model_levels) ! 3-D conv rain
REAL, INTENT(IN) :: dsnow(row_length,rows, model_levels)
! 3-D LS snowfall rate
REAL, INTENT(IN) :: csnow(row_length,rows, model_levels)
! 3-D conv snow rate
REAL, INTENT(IN) :: autoconv(row_length,rows, model_levels)
! Autoconversion rate kg/kg/s
REAL, INTENT(IN) :: accretion(row_length,rows, model_levels)
! Accretion rate kg/kg/s
REAL, INTENT(IN) :: rim_agg(row_length,rows, model_levels)
! Riming rate of aggregates kg/kg/s
REAL, INTENT(IN) :: rim_cry(row_length,rows, model_levels)
! Riming rate of ice crystals kg/kg/s
REAL, INTENT(IN) :: land_fraction(row_length,rows)     ! land_fraction
! in-cloud oxidation rates (molecules/cc/DTC) from h2o2 & o3 (UKCA):
REAL, INTENT(IN) :: delso2_wet_h2o2(row_length,rows,model_levels)
REAL, INTENT(IN) :: delso2_wet_o3  (row_length,rows,model_levels)
! in-air   oxidation rate  (molecules/cc/DTC) from oh        (UKCA):
! Note: When l_fix_ukca_h2so4_ystore=T this in (vmr/s).
!       For ASAD-based chemical schemes (e.g. StratTrop)
!       this is the change in H2SO4 from chemistry.
REAL, INTENT(IN) :: delso2_dry_oh   (row_length,rows,model_levels)
! in-cloud oxidation rates (kgS/kgair/s     ) from h2o2 & o3 (CLASSIC):
!      REAL, INTENT(IN) :: delso2_wet_h2o2C(row_length,rows,model_levels)
!      REAL, INTENT(IN) :: delso2_wet_o3C  (row_length,rows,model_levels)
! in-air   oxidation rate  (kgS/kgair/s     ) from oh        (CLASSIC):
!      REAL, INTENT(IN) :: delso2_dry_ohC  (row_length,rows,model_levels)
! cloud fraction
REAL, INTENT(IN) :: cloud_frac(row_length, rows, model_levels)
REAL, INTENT(IN) :: cloud_liq_frac(row_length, rows, model_levels)
REAL, INTENT(IN) :: cloud_liq_wat(row_length, rows, model_levels)
REAL, INTENT(IN) :: mass(row_length,rows, model_levels)
REAL, INTENT(IN) :: zbl(row_length, rows)  ! BL height

! names of chemistry tracers
CHARACTER(LEN=maxlen_fieldname), INTENT(IN) :: chemistry_tracer_names(         &
                                                   n_chemistry_tracers)

! names of aerosol tracers
CHARACTER(LEN=maxlen_fieldname), INTENT(IN) :: mode_tracer_names(n_mode_tracers)

! chemistry tracer mass mixing ratios
REAL, INTENT(IN OUT) :: chemistry_tracers(row_length,rows,                     &
                            model_levels,n_chemistry_tracers)
! aerosol tracer mass mixing ratio
REAL, INTENT(IN OUT) :: mode_tracers(row_length,rows,                          &
                            model_levels,n_mode_tracers)
! 3-D diagnostic array
REAL, INTENT(IN OUT) :: mode_diags(row_length,rows,                            &
                            model_levels,n_mode_diags)

REAL, INTENT(IN) :: z_half_alllevs(1:row_length,1:rows,                        &
                                   1:model_levels)

! Non transported prognostics
TYPE(ntp_type), INTENT(IN OUT) :: all_ntp(dim_ntp)

! Segmentation
INTEGER, INTENT(IN) :: nseg, stride_s
INTEGER, INTENT(IN) :: lbase(nseg)    ! the box where the seg starts in 3d array
INTEGER, INTENT(IN) :: ncol_s(nseg)   ! the number of columns in each segment
INTEGER, INTENT(IN) :: nbox_s(nseg)   ! the size of each segment (grid-boxes)

! Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: ddplim0(:)
REAL,    POINTER :: ddplim1(:)
REAL,    POINTER :: mfrac_0(:,:)
REAL,    POINTER :: mlo(:)
REAL,    POINTER :: mm(:)
REAL,    POINTER :: mmid(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: modesol(:)
INTEGER, POINTER :: mode_choice(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: num_eps(:)
REAL,    POINTER :: sigmag(:)

! Relate to segmentation
INTEGER :: lb,ncs,nbs     ! short hand for element of lbase, ncol_s and nbox_s
INTEGER :: ik,ic          ! loop iterators for segments and columns in a segment
INTEGER :: tid_omp        ! thread id for parallel region

REAL :: a3d_tmp(row_length,rows,model_levels) ! temporary 3D array re-used


INTEGER, PARAMETER :: nhet = 2
! Number of heterogeneous reaction rates

INTEGER, PARAMETER :: nmts=1
! No. of microphysical sub-steps per DTC
INTEGER :: nzts
! No. of condensation-nucleation competition sub-steps per DTM
INTEGER :: iextra_checks
! Level of protection against bad MDT values following advection
INTEGER :: rainout_on
! Switch for whether rainout (nucl. scav.) is on/off
INTEGER :: imscav_on
! Switch for whether impaction scavenging is on/off
INTEGER, PARAMETER :: wetox_on=1
! Switch for whether wet oxidation (cloud processing) is on/off
INTEGER :: ddepaer_on
! Switch for whether aerosol dry deposition is on/off
INTEGER, PARAMETER :: sedi_on=1
! Switch for whether aerosol sedimentation is on/off
INTEGER, PARAMETER :: iso2wetoxbyo3=1
! Switch for whether SO2 wet oxidation by ozone is on/off
! Note that this switch is only used if WETOX_IN_AER=1
! When code used in UM    , WETOX_IN_AER is always set to 0
! When code used in TOMCAT, WETOX_IN_AER is always set to 1
INTEGER, PARAMETER :: cond_on=1
! Switch for whether vapour condensation is  on/off
INTEGER :: nucl_on
! Switch for whether binary nucleation is on/off
INTEGER :: bln_on
! Switch for whether binary BL nucleation is on/off
INTEGER, PARAMETER :: coag_on=1
! Switch for whether coagulation is on/off
INTEGER, PARAMETER :: icoag=1
! Switch for KIJ method (1:GLOMAP, 2: M7, 3: UMorig, 4:UMorig MFPP)
!   =3 Cunnigham scheme as in UM, =4 as in UM but computing values)
INTEGER, PARAMETER :: imerge=2
! Switch to use mid-pts (=1), edges (2) or dynamic (=3) in remode
INTEGER, PARAMETER :: ifuchs=2
! Switch for Fuchs(1964) (=1) or Fuchs-Sutugin(1971) for CC (=2)
INTEGER, PARAMETER :: idcmfp=2
! Switch for vapour-diffusion-method (1=as bin v1, 2=as bin v1.1)
INTEGER, PARAMETER :: icondiam=2
! Switch for what diameter to use for CONDEN (1=g.m.diam, 2=conden-diam)
INTEGER, PARAMETER :: i_nuc_method=2
!  I_NUC_METHOD: Switch for nucleation (how to combine BHN and BLN)
! (1=initial Pandis94 approach (no BLN even if switched on) -- Do not use!!
! (2=binary homogeneous nucleation applying BLN to BL only if switched on)
!   note there is an additional switch i_bhn_method (local to CALCNUCRATE)
!   to switch between using Kulmala98 or Vehkamaki02 for BHN rate
! (3=use same rate at all levels either activation(IBLN=1), kinetic(IBLN=2),
! PNAS(IBLN=3)
!  note that if I_NUC_METHOD=3 and IBLN=3 then also add on BHN rate as in PNAS.
INTEGER :: ibln
! Switch for BLN parametrisation rate calc (1=activation,2=kinetic,3=PNAS)
INTEGER :: iactmethod
! Switch for activation method (0=off,1=fixed ract,2=NSO3 scheme)
INTEGER :: inucscav
! Switch for nucl scav method (1=as GLOMAP Spr05, 2=use M7 scav coeffs,
!                              3=as (1) but no nucl scav of modes 6 & 7)
INTEGER, PARAMETER :: iddepaer=2
! Switch for dry dep method (1=as GLOMAP Spr05, 2=incl. sedi)
!      INTEGER, PARAMETER :: IDDEPAER=1
! Switch for dry dep method (1=as GLOMAP Spr05, 2=incl. sedi)
! INTEGER :: VERBOSE: copy of glob_verbose from UKCA_MODE_SETUP
! Switch to determine level of debug output (0=none, 1, 2)
! Now set based on value of PRINT_STATUS in UKCA_MAIN
INTEGER :: verbose_local ! local copy of glob_verbose to set independently
! from value of PRINT_STATUS in UKCA_MAIN (for GLOMAP de-bug
! uses this value within GLOMAP routines (other print statements
! in UKCA_AERO_CTL still controlled by PRINT_STATUS)
INTEGER, PARAMETER :: checkmd_nd=1
! Switch for whether to check for bad values of MD and ND
INTEGER, PARAMETER :: intraoff=0
! Switch to turn off intra-modal coagulation
INTEGER, PARAMETER :: interoff=0
! Switch to turn off inter-modal coagulation
INTEGER, PARAMETER :: idustems=0
! Switch for using Pringle scheme (=1) or AEROCOMdaily (=2)
INTEGER, PARAMETER :: iagecoagnucl67=0
! Switch to enable(1)/disable(0) ageing of modes 6&7 (insol dust only) and
! coagulation & nucleation involving those modes

REAL :: dp0                               ! Diam (nm)
!
REAL :: dtm
! Microphysics time step (s)
REAL :: dtz
! Competition (cond/nucl) time step (s)

REAL :: y(nmodes)
! EXP(2*LN(SIGMA)*LN(SIGMA)) for each mode (for surf area conc calc)

! Molar mass of dry air (kg/mol)
REAL :: mm_da  ! =avogadro*boltzmann/rgas

! Items for MDT too low/hi check
REAL :: mdtmin(nmodes,nseg)

INTEGER :: field_size              ! size of 2D field
INTEGER :: field_size3d            ! size of 3D field

INTEGER :: i,j,k,l,n,jl,imode,icp
INTEGER :: n_reqd_tracers           ! No of tracers required
INTEGER :: itra

REAL :: scale_delso2
! Scaling factor for the in-cloud production rates delso2,delso2_2
! to account for removal by precipitation before the cloud
! evaporates.
REAL, PARAMETER :: ma=4.78e-26 ! mass of air molecule (kg)
REAL :: act                    ! radius for activation (m)

LOGICAL :: lcvrainout
! Switch for convective rainout (.FALSE. if done with convective transport)

LOGICAL :: l_dust_slinn_impc_scav
! Switch to turn on the new impaction scavenging scheme for dust

! used for debug output
LOGICAL (KIND=log_small), ALLOCATABLE, SAVE :: mode_tracer_debug(:)

CHARACTER(LEN=errormessagelength) :: cmessage     ! Error message

LOGICAL :: logic  ! for use in storing of aerosol budget terms
LOGICAL :: logic1 ! for use in storing of aerosol budget terms
LOGICAL :: logic2 ! for use in storing of aerosol budget terms
LOGICAL, SAVE :: firstcall=.TRUE.
! counter: mode-merges applied
INTEGER (KIND=integer_32) :: n_merge_3d(row_length,rows,model_levels,nmodes)
INTEGER (KIND=integer_32) :: sum_nbadmdt(nmodes)
INTEGER (KIND=integer_32) :: thread_sum_nbadmdt
INTEGER (KIND=integer_32) :: nbadmdt_3d(row_length,rows,model_levels,nmodes)
INTEGER :: jv
INTEGER :: ifirst    ! index of first mode tracer in nmr_index, mmr_index

LOGICAL :: l_ukca_segment_uniform   ! are the cache-blocking segments uniform
! This taken out of run_ukca as set here for now
INTEGER :: i_mode_act_method

REAL :: root2                ! square root of 2
REAL :: log_sigmag(nmodes)

INTEGER :: thread_min
INTEGER :: thread_max
INTEGER :: num_of_threads

TYPE(segment_data_type) :: seg

INTEGER :: nseg_per_thread
INTEGER :: seg_remainder
INTEGER :: thread_nseg

INTEGER :: errcode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_AERO_CTL'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
ddplim0     => glomap_variables%ddplim0
ddplim1     => glomap_variables%ddplim1
mfrac_0     => glomap_variables%mfrac_0
mlo         => glomap_variables%mlo
mm          => glomap_variables%mm
mmid        => glomap_variables%mmid
mode        => glomap_variables%mode
modesol     => glomap_variables%modesol
mode_choice => glomap_variables%mode_choice
ncp         => glomap_variables%ncp
num_eps     => glomap_variables%num_eps
sigmag      => glomap_variables%sigmag

root2 = SQRT(2.0)

CALL segment_data_nullify(seg)

IF (glomap_config%l_ukca_radaer) THEN
  CALL ukca_check_radaer_coupling(all_ntp)
END IF

! Molar mass of dry air (kg/mol)
mm_da = avogadro*boltzmann/rgas

!-------------------------------------
! As well as VERBOSE being set from PrintStatus, also have local version
! of VERBOSE which can be set independently for debugging purposes
! Call this "VERBLOC" for local version of VERBOSE -- this is what gets
! used within GLOMAP routines. The reason it needed to be changed is that
! this GLOMAP switch is designed for de-bugging and not meant to be mapped
! to a top-level switch in UM --- prefer "devolved verbosity control" here.

!------------------------------------------
verbose_local = 0
! verbose_local = 0 is default setting for when running model experiments
! verbose_local = 2 is investigative setting for print statements to check
!             evolution of size-resolved aerosol properties
!             (number, mass, size) after each process
!------------------------------------------
! Set BL nucleation parametrisation rate calculation method
ibln = glomap_config%i_mode_bln_param_method

! Set scavenging coefficients
CALL ukca_mode_imscavcoff(verbose)
IF (verbose >= 2) THEN
  WRITE(umMessage,'(A45,2I6)') 'Set up aerosol etc., NTRAER,NBUDAER=',         &
                                   ntraer,nbudaer
  CALL umPrint(umMessage,src='ukca_aero_ctl')

  DO i=1,nmodes
    WRITE(umMessage,'(A45,I5,L7,3E12.3)')                                      &
              'I,MODE(I),DDPLIM0(I),DDPLIM1(I),SIGMAG(I)=',                    &
               i,mode(i),ddplim0(i),ddplim1(i),sigmag(i)
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    DO j=1,ncp
      WRITE(umMessage,'(A35,2I5,L7,E12.3)')                                    &
                 'I,J,COMPONENT(I,J),MFRAC_0(I,J)=',                           &
                  i,j,component(i,j),mfrac_0(i,j)
      CALL umPrint(umMessage,src='ukca_aero_ctl')
    END DO
  END DO
END IF

! below are the input configuration parameters set for MODE

IF (firstcall .AND. verbose >0 ) THEN
  WRITE(umMessage,'(A25,I6)') 'i_mode_setup =      ',                          &
                              glomap_config%i_mode_setup
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A25,F6.2)') 'mode_parfrac =      ',                        &
                                ukca_config%mode_parfrac
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A25,I6)') 'i_mode_nucscav =    ',                          &
                              glomap_config%i_mode_nucscav
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A25,I6)') 'i_mode_nzts =       ',                          &
                              glomap_config%i_mode_nzts
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A25,L7)') 'l_mode_bhn_on =     ',                          &
                              glomap_config%l_mode_bhn_on
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A25,L7)') 'l_mode_bln_on =     ',                          &
                              glomap_config%l_mode_bln_on
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A25,I6)') 'i_mode_bln_param_method',                       &
                              glomap_config%i_mode_bln_param_method
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF

! Set whether or not dry deposition of aerosols is on
IF (glomap_config%l_ddepaer) THEN
  ddepaer_on=1
ELSE
  ddepaer_on=0
END IF
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage,'(A22,L6,I5)') 'L_DDEPAER,DDEPAER_ON=',                      &
                                 glomap_config%l_ddepaer, ddepaer_on
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF

! Set whether or not wet deposition (nucleation and impaction scavenging) of
! aerosols is on
IF (glomap_config%l_rainout) THEN
  rainout_on=1
ELSE
  rainout_on=0
END IF
IF (glomap_config%l_impc_scav) THEN
  imscav_on=1
ELSE
  imscav_on=0
END IF
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage,'(A22,L6,I5)') 'L_RAINOUT,RAINOUT_ON=',                      &
                                 glomap_config%l_rainout, rainout_on
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A24,L6,I5)') 'L_IMPC_SCAV,IMSCAV_ON=',                     &
                                 glomap_config%l_impc_scav, imscav_on
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF

! set INUCSCAV according to i_mode_nucscav setting
inucscav = glomap_config%i_mode_nucscav
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage,'(A25,2I5)') 'I_MODE_NUCSCAV,INUCSCAV=',                     &
                                glomap_config%i_mode_nucscav,inucscav
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF

! set LCVRAINOUT according to l_cv_rainout setting
lcvrainout = glomap_config%l_cv_rainout
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage,'(A25,2L6)') 'L_CV_RAINOUT,LCVRAINOUT=',                     &
                                glomap_config%l_cv_rainout, lcvrainout
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF

! If new Slinn impaction scavenging scheme is on for dust then turn on flag
! (input to UKCA_AERO_STEP) and initialise scavenging arrays
l_dust_slinn_impc_scav = glomap_config%l_dust_slinn_impc_scav
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage,'(A24,L7)') 'L_DUST_SLINN_IMPC_SCAV=',                       &
                               glomap_config%l_dust_slinn_impc_scav
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF
IF (l_dust_slinn_impc_scav) CALL ukca_impc_scav_dust_init(verbose)

! Set IEXTRA_CHECKS to control checking for unacceptable MDT values
! With IEXTRA_CHECKS = 0, code in UKCA_AERO_CTL checks that mdt is above
! a minimum value after advection, and then resets nd (aerosol number density)
! and md (component mass per particle) to their default values where mdt is
! too low.
! With IEXTRA_CHECKS = 1, the routine UKCA_CHECK_ARTEFACTS is called
! from UKCA_AERO_CTL and mdt is checked for values above and below established
! limits for the mode, then nd and md reset to their default values for the
! exceptions.
! For IEXTRA_CHECKS = 2 (the recommended value), post advection checks are
! carried out as for IEXTRA_CHECKS = 1, and the routine  UKCA_CHECK_MDT
! is called from UKCA_COAGWITHNUCL to check that mdt remains with upper and
! lower limits with nd and md reset to their defaults for exceptions.
iextra_checks = 2
! set NZTS according to i_mode_nzts setting
nzts = glomap_config%i_mode_nzts
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage,'(A15,2I5)') 'I_MODE_NZTS,NZTS=',                            &
                               glomap_config%i_mode_nzts,nzts
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF
!
i_mode_act_method=1
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage,'(A44)') 'setting I_MODE_ACT_METHOD to default value'
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF
!
! set IACTMETHOD according to i_mode_act_method setting
iactmethod=i_mode_act_method
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage,'(A30,2I6)') 'I_MODE_ACT_METHOD,IACTMETHOD=',                &
      i_mode_act_method,iactmethod
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF
!
! set ACT according to mode_activation_dryr setting
IF (ABS(glomap_config%mode_activation_dryr - rmdi) < EPSILON(0.0)) THEN
  cmessage = ' mode_activation_dryr has not been set'
  errcode = 1
  WRITE(umMessage,'(A40)') cmessage
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
ELSE
  act = glomap_config%mode_activation_dryr*1.0e-9 ! convert nm to m
END IF

! Setting removed fraction of oxidised SO2 to a default of 0
IF (ABS(glomap_config%mode_incld_so2_rfrac) < EPSILON(0.0)) THEN
  WRITE(umMessage,'(A57)') 'MODE_INCLD_SO2_RFRAC has been set' //              &
                           ' to default value of 0.0'
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF

! Calculate the production rate scaling factor
scale_delso2=1.0-glomap_config%mode_incld_so2_rfrac

IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage,'(A26,2E12.3)') 'MODE_ACTIVATION_DRYR,ACT=',                 &
                                  glomap_config%mode_activation_dryr, act
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A22,E12.3)') 'MODE_INCLD_SO2_RFRAC=',                      &
                                 glomap_config%mode_incld_so2_rfrac
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF
!
! set NUCL_ON according to L_MODE_BHN_ON & L_MODE_BLN_ON settings
IF (glomap_config%l_mode_bhn_on) THEN
  nucl_on = 1 ! BHN
ELSE
  nucl_on = 0 ! BHN
END IF

IF (glomap_config%l_mode_bln_on) THEN
  bln_on = 1 ! BLN
ELSE
  bln_on = 0 ! BLN
END IF
!
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage,'(A16,2I6)') 'NUCL_ON,BLN_ON=',nucl_on,bln_on
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A10,I6)') 'IBLN=',ibln
  CALL umPrint(umMessage,src='ukca_aero_ctl')
END IF
!
dtm=dtc/REAL(nmts)
dtz=dtm/REAL(nzts)
field_size=row_length*rows
field_size3d=field_size*model_levels

IF (firstcall) THEN
  ALLOCATE(mode_tracer_debug(n_mode_tracers))
  mode_tracer_debug(:)=.TRUE.     ! all tracers with debug o/p
END IF

IF (verbose > 1) THEN

  WRITE(umMessage,'(A32,3I6)') 'nbox,field_size,field_size3d=',                &
      nbox,field_size,field_size3d
  CALL umPrint(umMessage,src='ukca_aero_ctl')

  WRITE(umMessage,'(A)') 'UKCA_MODE INPUT SETTINGS : '
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A14,I6)') 'i_mode_setup=',glomap_config%i_mode_setup
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A14,F8.2)') 'mode_parfrac=',ukca_config%mode_parfrac
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A16,I6)') 'i_mode_nucscav=',glomap_config%i_mode_nucscav
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A16,I6)') 'i_mode_nzts=',glomap_config%i_mode_nzts
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A16,L1)') 'l_mode_bhn_on=',glomap_config%l_mode_bhn_on
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A16,L1)') 'l_mode_bln_on=',glomap_config%l_mode_bln_on
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A26,I6)') 'i_mode_bln_param_method',                       &
                              glomap_config%i_mode_bln_param_method
  CALL umPrint(umMessage,src='ukca_aero_ctl')

  WRITE(umMessage,'(A16,I6)') 'i_month: ',i_month
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A16,I6)') 'i_day_number: ',i_day_number
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A16,I6)') 'i_hour: ',i_hour
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A16,I6)') 'i_minute: ',i_minute
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A16,F8.2)') 'DTC: ',dtc
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A16,I6)') 'model_levels: ',model_levels
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A16,I6)') 'rows: ',rows
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A16,I6)') 'row_length: ',row_length
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A22,I6)') 'n_chemistry_tracers: ',n_chemistry_tracers
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A18,I6)') 'n_mode_tracers: ',n_mode_tracers
  CALL umPrint(umMessage,src='ukca_aero_ctl')

  WRITE(umMessage,'(A40)') 'Array:     MIN        MAX         MEAN'
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  l=0
  WRITE(umMessage,'(A9,I6)') 'Level: ',l
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  WRITE(umMessage,'(A8,3E12.3)') 'p_bdrs: ',MINVAL(p_bdrs(:,:,l)),             &
                        MAXVAL(p_bdrs(:,:,l)),                                 &
                        SUM(p_bdrs(:,:,l))/REAL(SIZE(p_bdrs(:,:,l)))
  CALL umPrint(umMessage,src='ukca_aero_ctl')

  ! No model level 2 in UKCA box model - diagnostic print statements
  DO l=1,MIN(2, model_levels)            ! model_levels
    WRITE(umMessage,'(A9,I6)') 'Level: ',l
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A8,3E12.3)') 'pres: ',MINVAL(pres(:,:,l)),               &
               MAXVAL(pres(:,:,l)),                                            &
               SUM(pres(:,:,l))/REAL(SIZE(pres(:,:,l)))
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A8,3E12.3)') 'temp: ',MINVAL(temp(:,:,l)),               &
               MAXVAL(temp(:,:,l)),                                            &
               SUM(temp(:,:,l))/REAL(SIZE(temp(:,:,l)))
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A8,3E12.3)') 'q: ',MINVAL(q(:,:,l)),MAXVAL(q(:,:,l)),    &
               SUM(q(:,:,l))/REAL(SIZE(q(:,:,l)))
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A8,3E12.3)') 'rh3d: ',MINVAL(rh3d(:,:,l)),               &
               MAXVAL(rh3d(:,:,l)),                                            &
               SUM(rh3d(:,:,l))/REAL(SIZE(rh3d(:,:,l)))
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A8,3E12.3)') 'p_bdrs: ',MINVAL(p_bdrs(:,:,l)),           &
               MAXVAL(p_bdrs(:,:,l)),                                          &
               SUM(p_bdrs(:,:,l))/REAL(SIZE(p_bdrs(:,:,l)))
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A18,3E12.3)') 'delso2_wet_h2o2: ',                       &
               MINVAL(delso2_wet_h2o2(:,:,l)),                                 &
               SUM(delso2_wet_h2o2(:,:,l))/                                    &
               REAL(SIZE(delso2_wet_h2o2(:,:,l)))
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A18,3E12.3)') 'delso2_wet_o3  : ',                       &
               MINVAL(delso2_wet_o3(:,:,l)),                                   &
               SUM(delso2_wet_o3(:,:,l))/REAL(SIZE(delso2_wet_o3(:,:,l)))
  END DO
  IF (model_levels > 7) THEN
    l = 8
    WRITE(umMessage,'(A18,I4,3E12.3)') 'delso2_wet_h2o2: ',l,                  &
               MINVAL(delso2_wet_h2o2(:,:,l)),                                 &
               MAXVAL(delso2_wet_h2o2(:,:,l)),                                 &
               SUM(delso2_wet_h2o2(:,:,l))/REAL(SIZE(delso2_wet_h2o2(:,:,l)))
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A18,I4,3E12.3)') 'delso2_wet_o3  : ',l,                  &
               MINVAL(delso2_wet_o3(:,:,l)),                                   &
               MAXVAL(delso2_wet_o3(:,:,l)),                                   &
               SUM(delso2_wet_o3(:,:,l))/REAL(SIZE(delso2_wet_o3(:,:,l)))
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A18,I4,3E12.3)') 'delso2_dry_oh  : ',l,                  &
               MINVAL(delso2_dry_oh(:,:,l)),                                   &
               MAXVAL(delso2_dry_oh(:,:,l)),                                   &
               SUM(delso2_dry_oh(:,:,l))/SIZE(delso2_dry_oh(:,:,l))
    CALL umPrint(umMessage,src='ukca_aero_ctl')
  END IF

  DO l=1,MIN(2, model_levels)            ! model_levels
    DO j=1,n_chemistry_tracers
      WRITE(umMessage,'(A8,I4,A10,I4,A10)') 'Level: ',l,' Tracer: ',j,         &
           chemistry_tracer_names(j)
      CALL umPrint(umMessage,src='ukca_aero_ctl')
      WRITE(umMessage,'(A20,3E12.3)') 'chemistry_tracers: ',                   &
           MINVAL(chemistry_tracers(:,:,l,j)),                                 &
           MAXVAL(chemistry_tracers(:,:,l,j)),                                 &
           SUM(chemistry_tracers(:,:,l,j))/                                    &
           REAL(SIZE(chemistry_tracers(:,:,l,j)))
      CALL umPrint(umMessage,src='ukca_aero_ctl')
    END DO
    DO j=1,n_mode_tracers
      IF (mode_tracer_debug(j)) THEN
        WRITE(umMessage,'(A8,I4,A10,I4,A10)') 'Level: ',l,' Tracer: ',j,       &
             mode_tracer_names(j)
        CALL umPrint(umMessage,src='ukca_aero_ctl')
        WRITE(umMessage,'(A18,3E12.3)') 'mode_tracers: ',                      &
                MINVAL(mode_tracers(:,:,l,j)),                                 &
                MAXVAL(mode_tracers(:,:,l,j)),                                 &
                SUM(mode_tracers(:,:,l,j))/REAL(SIZE(mode_tracers(:,:,l,j)))
        CALL umPrint(umMessage,src='ukca_aero_ctl')
      END IF
    END DO
  END DO     ! model_levels

END IF ! IF (verbose > 1)

! Calculate number of aerosol tracers required for components and number
n_reqd_tracers = 0
DO imode=1,nmodes
  DO icp=1,ncp
    IF (component(imode,icp)) n_reqd_tracers = n_reqd_tracers + 1
  END DO
END DO
n_reqd_tracers = n_reqd_tracers + SUM(mode_choice)

IF (firstcall) THEN

  ! .. Check the number of tracers, warn if too many, stop if too few
  IF (n_mode_tracers > n_reqd_tracers) THEN
    errcode=-1
    cmessage=' Too many tracers input'
    WRITE(umMessage,'(A50,2I5)') cmessage,n_mode_tracers,n_reqd_tracers
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
  END IF
  IF (n_mode_tracers < n_reqd_tracers) THEN
    errcode=1
    cmessage=' Too few advected aerosol tracers input'
    WRITE(umMessage,'(A50,2I5)') cmessage,n_mode_tracers,n_reqd_tracers
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
  END IF

  ! Check that all tracer addresses are in range  (nmr_index and mmr_index are
  ! for all UKCA tracers, so subtract (ifirst-1) to index the mode tracer array)

  ifirst = jpctr + 1
  IF (verbose > 0) THEN
    WRITE(umMessage,'(A43)')'Checking MODE tracer addresses are in range:'
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A40)') 'Description, imode, [icp], ifirst, itra'
    CALL umPrint(umMessage,src='ukca_aero_ctl')
  END IF
  DO imode=1,nmodes
    IF (mode(imode)) THEN
      itra = nmr_index(imode) - ifirst + 1
      icp = 0
      IF (verbose > 0) THEN
        WRITE(umMessage,'(A10,4I6)') 'Number:  ',imode,icp,ifirst,itra
        CALL umPrint(umMessage,src='ukca_aero_ctl')
      END IF
      IF (itra <= 0 .OR. itra > n_mode_tracers) THEN
        errcode=1
        cmessage='Tracer address out of range for number'
        WRITE(umMessage,'(A72,2(A8,I6))') cmessage,' mode: ',imode,            &
                                  ' index: ',nmr_index(imode)
        CALL umPrint(umMessage,src='ukca_aero_ctl')
        CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
      END IF
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          itra = mmr_index(imode,icp) - ifirst + 1
          IF (verbose > 0) THEN
            WRITE(umMessage,'(A10,4I6)')'Mass MR: ',imode, icp, ifirst,        &
                 itra
            CALL umPrint(umMessage,src='ukca_aero_ctl')
          END IF
          IF (itra < 0 .OR. itra > n_mode_tracers) THEN
            errcode=1
            cmessage='Tracer address out of range for component'
            WRITE(umMessage,'(A72,3(A12,I6))') cmessage,' mode: ',imode,       &
                    ' component: ',icp,' index: ',nmr_index(imode)
            CALL umPrint(umMessage,src='ukca_aero_ctl')
            CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
          END IF
        END IF
      END DO     ! icp
    END IF     ! mode(imode)
  END DO     ! imode

  ! Allocate the coefficients required for aerosol dry deposition.
  ! The values used are chosen based on l_improve_aero_drydep
  ! The arrays are independent of location, so can be initialised once.
  CALL set_ddepaer_coeff(glomap_config%l_improve_aero_drydep)

END IF   ! firstcall

! the all_ntp array is at end of segment loop in this case to accumulate SA
i = name2ntpindex('surfarea  ')
all_ntp(i)%data_3d(:,:,:) = 0.0
n_merge_3d(:,:,:,:)=0
! temporary clumping of autoconv variables for later segment extraction
a3d_tmp(:,:,:) = autoconv(:,:,:)+accretion(:,:,:)+rim_agg(:,:,:)+rim_cry(:,:,:)

! this is used to track the number of cells where MDT was modified due to
! unphysical tracer values
sum_nbadmdt(1:nmodes) = 0

! other precomputations
CALL log_v(nmodes,sigmag,log_sigmag)

! Allocate the drydiam array (input to UKCA_ACTIVATE) if required.
IF ( ( glomap_config%i_ukca_activation_scheme == i_ukca_activation_arg )       &
     .AND. .NOT. ALLOCATED(drydiam) ) THEN
  ALLOCATE(drydiam(row_length,rows,model_levels,nmodes))
  drydiam(:,:,:,:) = 0.0
END IF

!
! OpenMP parallel region starts here (loop later around IK segments)
!
!$OMP PARALLEL  DEFAULT(NONE)                                                  &
!$OMP          SHARED(a3d_tmp, act, all_ntp, bln_on,                           &
!$OMP chemistry_tracers, cloud_frac,                                           &
!$OMP cloud_liq_frac, cloud_liq_wat, component, crain, csnow, delso2_dry_oh,   &
!$OMP delso2_wet_h2o2, delso2_wet_o3, drain, dryox_in_aer, dsnow,              &
!$OMP dtc, dtm, dtz, firstcall, iactmethod, ibln, inucscav,                    &
!$OMP jpctr, lbase,                                                            &
!$OMP ukca_config, glomap_config, glomap_variables,                            &
!$OMP l_ukca_mode_diags, l_ukca_cmip6_diags,                                   &
!$OMP l_ukca_pm_diags,                                                         &
!$OMP land_fraction, lcvrainout,  mass, mfrac_0, mh2o2f, mh2so4,               &
!$OMP mlo, mm, mm_da, mm_gas, mmid, mmr_index, mode, mode_diags, mdwat_diag,   &
!$OMP wetdp_diag, l_dust_slinn_impc_scav,                                      &
!$OMP mode_tracers, model_levels, modesol, mox, msec_org, msec_orgi, msotwo,   &
!$OMP n_h2o2, n_h2so4, n_merge_3d,                                             &
!$OMP n_o3, n_sec_org, n_sec_org_i, n_so2, nadvg, nbadmdt_3d, nbox, nbox_s,    &
!$OMP nbudaer, nchemg, ncol_s, ncp,                                            &
!$OMP nmasddepsunucsol, nmasddepsuaitsol, nmasddepsuaccsol, nmasddepsucorsol,  &
!$OMP nmasddepssaccsol, nmasddepsscorsol, nmasddepbcaitsol, nmasddepbcaccsol,  &
!$OMP nmasddepbccorsol, nmasddepbcaitins, nmasddepocnucsol, nmasddepocaitsol,  &
!$OMP nmasddepocaccsol, nmasddepoccorsol, nmasddepocaitins, nmasddepsonucsol,  &
!$OMP nmasddepsoaitsol, nmasddepsoaccsol, nmasddepsocorsol, nmasddepduaccsol,  &
!$OMP nmasddepducorsol, nmasddepduaccins, nmasddepducorins,                    &
!$OMP nmasnuscsunucsol, nmasnuscsuaitsol, nmasnuscsuaccsol, nmasnuscsucorsol,  &
!$OMP nmasnuscssaccsol, nmasnuscsscorsol, nmasnuscbcaitsol, nmasnuscbcaccsol,  &
!$OMP nmasnuscbccorsol, nmasnuscbcaitins, nmasnuscocnucsol, nmasnuscocaitsol,  &
!$OMP nmasnuscocaccsol, nmasnuscoccorsol, nmasnuscocaitins,                    &
!$OMP nmasnuscsonucsol, nmasnuscsoaitsol, nmasnuscsoaccsol, nmasnuscsocorsol,  &
!$OMP nmasnuscduaccsol, nmasnuscducorsol, nmasnuscduaccins, nmasnuscducorins,  &
!$OMP nmasimscsunucsol, nmasimscsuaitsol, nmasimscsuaccsol, nmasimscsucorsol,  &
!$OMP nmasimscssaccsol, nmasimscsscorsol, nmasimscbcaitsol, nmasimscbcaccsol,  &
!$OMP nmasimscbccorsol, nmasimscbcaitins, nmasimscocnucsol, nmasimscocaitsol,  &
!$OMP nmasimscocaccsol, nmasimscoccorsol, nmasimscocaitins, nmasimscsonucsol,  &
!$OMP nmasimscsoaitsol, nmasimscsoaccsol, nmasimscsocorsol, nmasimscduaccsol,  &
!$OMP nmasimscducorsol, nmasimscduaccins, nmasimscducorins, nmasclprsuaitsol1, &
!$OMP nmasclprsuaccsol1,nmasclprsucorsol1,nmasclprsuaitsol2,nmasclprsuaccsol2, &
!$OMP nmasclprsucorsol2, nmasprocsuintr23, nmasprocbcintr23, nmasprococintr23, &
!$OMP nmasprocsointr23, nmascondsunucsol, nmascondsuaitsol, nmascondsuaccsol,  &
!$OMP nmascondsucorsol, nmascondsuaitins, nmascondsuaccins, nmascondsucorins,  &
!$OMP nmascondocnucsol, nmascondocaitsol, nmascondocaccsol, nmascondoccorsol,  &
!$OMP nmascondocaitins, nmascondocaccins, nmascondoccorins, nmascondocinucsol, &
!$OMP nmascondociaitsol,nmascondociaccsol,nmascondocicorsol,nmascondociaitins, &
!$OMP nmascondociaccins,nmascondocicorins,nmascondsonucsol,                    &
!$OMP nmascondsoaitsol, nmascondsoaccsol, nmascondsocorsol, nmascondsoaitins,  &
!$OMP nmascondsoaccins, nmascondsocorins, nmasnuclsunucsol, nmascoagsuintr12,  &
!$OMP nmascoagsuintr13, nmascoagsuintr14, nmascoagsuintr15, nmascoagsuintr16,  &
!$OMP nmascoagsuintr17, nmascoagocintr12, nmascoagocintr13, nmascoagocintr14,  &
!$OMP nmascoagocintr15, nmascoagocintr16, nmascoagocintr17, nmascoagsointr12,  &
!$OMP nmascoagsointr13, nmascoagsointr14, nmascoagsointr15, nmascoagsointr16,  &
!$OMP nmascoagsointr17, nmascoagsuintr23, nmascoagsuintr24, nmascoagbcintr23,  &
!$OMP nmascoagbcintr24, nmascoagocintr23, nmascoagocintr24, nmascoagsointr23,  &
!$OMP nmascoagsointr24, nmascoagsuintr34, nmascoagbcintr34, nmascoagocintr34,  &
!$OMP nmascoagssintr34, nmascoagsointr34, nmascoagduintr34, nmascoagbcintr53,  &
!$OMP nmascoagocintr53, nmascoagbcintr54, nmascoagocintr54, nmascoagduintr64,  &
!$OMP nmasagedsuintr52, nmasagedbcintr52, nmasagedocintr52, nmasagedsointr52,  &
!$OMP nmasmergsuintr12, nmasmergocintr12, nmasmergsointr12, nmasmergsuintr23,  &
!$OMP nmasmergbcintr23, nmasmergocintr23, nmasmergsointr23, nmasmergsuintr34,  &
!$OMP nmasmergbcintr34, nmasmergocintr34, nmasmergssintr34, nmasmergduintr34,  &
!$OMP nmasmergsointr34,                                                        &
!$OMP nmasddepntaitsol, nmasddepntaccsol,                                      &
!$OMP nmasddepntcorsol, nmasddepnhaitsol,nmasddepnhaccsol, nmasddepnhcorsol,   &
!$OMP nmasnuscntaitsol, nmasnuscntaccsol,nmasnuscntcorsol, nmasnuscnhaitsol,   &
!$OMP nmasnuscnhaccsol, nmasnuscnhcorsol,nmasimscntaitsol, nmasimscntaccsol,   &
!$OMP nmasimscntcorsol, nmasimscnhaitsol,nmasimscnhaccsol, nmasimscnhcorsol,   &
!$OMP nmascoagntintr23, nmascoagnhintr23,nmascoagntintr24, nmascoagnhintr24,   &
!$OMP nmascoagntintr34, nmascoagnhintr34,nmasmergntintr23, nmasmergnhintr23,   &
!$OMP nmasmergntintr34, nmasmergnhintr34,nmasprocntintr23, nmasprocnhintr23,   &
!$OMP nmasddepnnaccsol, nmasddepnncorsol,                                      &
!$OMP nmasnuscnnaccsol, nmasnuscnncorsol,nmasimscnnaccsol, nmasimscnncorsol,   &
!$OMP nmascoagnnintr34, nmasmergnnintr34,                                      &
!$OMP nmax_mode_diags,nmr_index,nucl_on, nseg,nukca_d1items,                   &
!$OMP num_eps, nzts, p_bdrs, pres, q, rgas, rh3d, rh3d_clr,                    &
!$OMP root2, sum_nbadmdt, mdtmin,                                              &
!$OMP row_length, rows, sea_ice_frac, scale_delso2, sigmag, log_sigmag,        &
!$OMP stride_s, temp, u_s, ukcaD1codes, verbose, wetox_in_aer, z0m,            &
!$OMP z_half_alllevs, zbl,  num_of_threads, nseg_per_thread, seg_remainder,    &
!$OMP drydiam, iextra_checks, ddepaer_on, rainout_on, imscav_on,               &
!$OMP verbose_local, lscat_zhang)                                              &
!$OMP     PRIVATE(errcode, cmessage, ifirst, dp0,  itra,                       &
!$OMP jl,  k, l_ukca_segment_uniform, l,  lb,                                  &
!$OMP logic, logic1, logic2, ncs, nbs,                                         &
!$OMP tid_omp, y, i, ic, ik, jv, n, imode, icp, j,                             &
!$OMP thread_min, thread_max, thread_nseg, thread_sum_nbadmdt, seg)

!Divide the number of segments by the number of available threads. Also compute
!the number of remainder segments to catch the case where the division is not
!equal.
!$OMP SINGLE
num_of_threads = 1
!$ num_of_threads = omp_get_num_threads()
nseg_per_thread = nseg / num_of_threads
seg_remainder = MOD(nseg, num_of_threads)
!$OMP END SINGLE

!Local thread ID
tid_omp = 0
!$ tid_omp = omp_get_thread_num()

!If the number of segments does not divide equally between threads, then low
!thread IDs get one additional segment to mop up the remainders.  Note that
!thread_nseg is specific to this thread. nseg_per_thread is shared, and is the
!minimum number of segments that all threads have.
IF (tid_omp < MOD(nseg, num_of_threads)) THEN
  thread_nseg = nseg_per_thread + 1
ELSE
  thread_nseg = nseg_per_thread
END IF

!Compute the first and last segment to be done on this specific thread.  The
!third term for thread_min --- the MIN() --- offsets thread_min upwards,
!accounting for the remainder segments that may be assigned to lower thread IDs.
thread_min = 1 + (tid_omp*nseg_per_thread) + MIN(tid_omp, seg_remainder)
thread_max = thread_min + thread_nseg - 1

! dealing with non-uniform segments
l_ukca_segment_uniform = .TRUE.
IF ( nseg > 1 ) THEN
  IF ( ncol_s(nseg) /= ncol_s(nseg-1) ) THEN

    WRITE(ummessage,'(A)') 'AERO_CTL: Segments are not uniform'
    CALL umprint(ummessage,src='ukca_aero_ctl',pe=0)
    WRITE(ummessage,'(i6,1x,i6,1x,i6)') ncol_s(nseg),ncol_s(nseg-1), nseg
    CALL umprint(ummessage,src='ukca_aero_ctl',pe=0)

    l_ukca_segment_uniform = .FALSE.  ! indicate the segments are non-uniform
  END IF  ! last segment is smaller
END IF  ! nseg > 1
! These allocated assuming uniform nbox elements in all segments
! therefore allocation ahead of ik loop

CALL segment_data_allocate(seg, ncp, nbox, nchemg, nhet, nbudaer, nadvg)

!Each thread acts on its allotted segments. Note that a simple OMP DO
!SCHEDULE(STATIC) has been found to change answers in certain tests, for as yet
!unknown reasons. Therefore, we avoid its use here and apply a custom
!worksharing algorithm.
DO ik = thread_min, thread_max     ! the segments on this MPI task

  ! use local alias because they are used in  many places within  loop
  lb = lbase(ik)      ! base location on this segment
  ncs = ncol_s(ik)    ! The number of columns on this segment
  nbs = nbox_s(ik)    ! The number of boxes on this segment
  IF (firstcall .AND. verbose > 2 ) THEN
    WRITE(umMessage,'(A28,5(1x,i8))') 'AERO_CTL:lb,ncs,nbs,nseg,ik',           &
                                                            lb,ncs,nbs,nseg,ik
    CALL umPrint(umMessage,src='ukca_aero_ctl')
  END IF

  IF ( (.NOT. l_ukca_segment_uniform) .AND. (ik == nseg)) THEN
    ! Ensure first dimension matches true number of boxes on this final segment.
    CALL segment_data_deallocate (seg)
    CALL segment_data_allocate(seg, ncp, nbs, nchemg, nhet, nbudaer, nadvg)

  END IF                                 ! if last segment smaller
  !
  ! now work out which grid box is above each grid box
  jl = 0
  DO ic = 1, ncs                        ! loop over the columns in this segment
    DO l=1,(model_levels-1)                  ! let top level have JLABOVE=-1
      jl         =jl + 1                     ! the ID of the current box
      seg%jlabove(jl)=jl + 1                 ! the ID of the box above this one
      seg%karr(jl)=MOD(lb,row_length) + ic-1 ! this is longitude
      seg%iarr(jl)=((lb+ic-1)/rows) + 1      ! this is latitude
      seg%larr(jl)=l                         ! this is altitude (level)
    END DO                                   ! l, model_levels
    jl = jl +1                               ! the box at top of column
    seg%jlabove(jl) = -1                     ! there is no box above this one
    seg%karr(jl)=MOD(lb,row_length) + ic -1  ! this is longitude
    seg%iarr(jl)=lb/rows                     ! this is latitude
    seg%larr(jl)=l                           ! this is altitude (level)
  END DO                                     ! for each column

  seg%n_merge_1d(:,:)=0

  ! Reshape input quantities
  ! ========================
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                   temp(1,1,1) ,seg%t)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                   pres(1,1,1), seg%pmid)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                   p_bdrs(1,1,1), seg%pupper)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                   p_bdrs(1,1,0), seg%plower)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                   crain(1,1,1), seg%craing)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                   drain(1,1,1), seg%draing)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                   csnow(1,1,1), seg%csnowg)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                   dsnow(1,1,1), seg%dsnowg)

  ! .. add in-cloud accretion, ice and snow melt to the autoconversion rate
  ! NOTE previously accumulated a3d_tmp = autoconv+accretion+rim_agg+rim_cry
  CALL extract_seg (lb,ncs,seg%nbox_this_seg,stride_s,model_levels,            &
                      a3d_tmp(1,1,1), seg%autoconv1d)

  ! .. set CRAING_UP using JLABOVE as calculated above
  DO jl=1,seg%nbox_this_seg   ! limit is number of boxes on this segment
    IF (seg%jlabove(jl) > 0) THEN
      seg%craing_up(jl)=seg%craing( seg%jlabove(jl))
    ELSE
      seg%craing_up(jl)=0.0
    END IF
  END DO
  !
  ! .. currently set FCONV_CONV=0.99 -- need to change to take as input
  seg%fconv_conv(:)=0.99 ! fraction of condensate-->rain in 6 hrs
  ! .. weakened rainout -- FCONV_CONV=0.5
  !!      FCONV_CONV(:)=0.50 ! fraction of condensate-->rain in 6 hrs
  !
  ! calculate molecular concentration of air
  seg%aird(:)=seg%pmid(:)/(seg%t(:)*boltzmann*1.0e6)  ! no conc of air (/cm3)

  ! copy from delso2_wet_xxx arrays as output from UKCA_CHEMISTRY_CTL
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                     delso2_wet_h2o2(1,1,1), seg%delso2)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                     delso2_wet_o3(1,1,1), seg%delso2_2)

  ! Scale in-cloud production rates to account for lack of
  ! in-cloud wet removal of oxidised SO2
  seg%delso2  (:) = seg%delso2(:)   * scale_delso2
  seg%delso2_2(:) = seg%delso2_2(:) * scale_delso2

  !
  ! set these to zero as will not be used in UKCA_AERO_STEP
  seg%zo3(:)=0.0     ! currently do wet ox separately in UM
  seg%zho2(:)=0.0    ! currently do wet ox separately in UM
  seg%zh2o2(:)=0.0   ! currently do wet ox separately in UM
  seg%lday(:)=0      ! currently do wet ox separately in UM
  !
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                     rh3d(1,1,1), seg%rh)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                   rh3d_clr(1,1,1),seg%rh_clr)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                     q(1,1,1), seg%s)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                     cloud_liq_wat(1,1,1),seg%lwc)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                     cloud_liq_wat(1,1,1),seg%clwc)

  seg%lowcloud(:) = 0.0
  seg%vfac(:)     = 0.0
  seg%clf(:)      = 0.0

  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                     cloud_liq_frac(1,1,1), seg%clf)
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                     cloud_frac(1,1,1), seg%lowcloud)

  WHERE (seg%lowcloud > 0.0)
    seg%vfac = 1.0 ! set to 1 so that VFAC*LOWCLOUD=cloud_frac
  END WHERE

  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                     z_half_alllevs(1,1,1), seg%height)

  ! GM added code here to only set cloud fraction > 0 if in low cloud
  !    here low cloud is defined as being cloud with p>=680hPa
  seg%mask1(:)=(seg%pmid(:) < 680.0e2) ! PMID is in Pa
  WHERE (seg%mask1(:))
    seg%lowcloud(:)=0.0
    seg%lwc(:)=0.0
  END WHERE
  !
  ! 2D surface -> whole segment
  CALL surface_to_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,          &
                        u_s(1,1),seg%ustr)
  CALL surface_to_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,          &
                        z0m(1,1), seg%znotg )
  CALL surface_to_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,          &
                        sea_ice_frac(1,1),seg%seaice)
  ! fraction of land at surface
  CALL surface_to_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,          &
                        land_fraction(1,1),seg%land_frac )
  CALL surface_to_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,          &
                        zbl(1,1),seg%htpblg)

  jl = 0
  DO ic = 1, ncs                 ! for each column on segment
    jl = jl +1                   ! at surface layer jl = 1
    IF ( seg%land_frac(jl) < 0.5 ) THEN
      seg%surtp(jl) = 0.0
    ELSE
      seg%surtp(jl) = 1.0
    END IF
    DO l = 2, model_levels       ! above the surface layer
      jl = jl +1
      IF ( seg%land_frac(jl) < 0.5 ) THEN
        seg%surtp(jl) = 2.0
      ELSE
        seg%surtp(jl) = 3.0
      END IF
    END DO
  END DO

  seg%ilscat(:) = 0

  ! Corrected surface types array, using dominant land_surface type rather
  ! than roughness length
  IF ( glomap_config%l_improve_aero_drydep ) THEN
    CALL surface_to_seg_int(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,    &
                        lscat_zhang(1,1),seg%ilscat)
  ELSE
    ! Old method
    !---------------------------------------------------------------
    ! Put in section here to set land-surface category ILSCAT based on ZNOTG
    ! ILSCAT=1-9 based on 9 UM landsurf types but use existing approach
    ! to set 4 possible types according to roughness length

    seg%mask1(:)=(seg%znotg(:) < 1.0e-3)  ! water/sea - z0<0.001m
    WHERE (seg%mask1(:)) seg%ilscat(:)=7

    ! forests - z0>0.1m
    seg%mask1(:)=(seg%znotg(:) > 1.0e-1) ! forest
    WHERE (seg%mask1(:)) seg%ilscat(:)=1

    ! all other lands, grass 0.001<z0<0.1m
    seg%mask1(:)=((seg%znotg(:) >= 1.0e-3)                                     &
            .AND. (seg%znotg(:) <= 1.0e-1)) ! grass
    WHERE (seg%mask1(:)) seg%ilscat(:)=3

    ! If sea ice covers > 50% of sea surface, treat as sea ice
    seg%mask1(:)=(seg%seaice(:) > 0.5) !seaice
    WHERE (seg%mask1(:)) seg%ilscat(:)=9

  END IF
  !---------------------------------------------------------------

  ! Derived quantities
  ! ==================
  seg%airdm3(:)=seg%aird(:)*1.0e6              ! no conc of air (/m3)
  seg%rhoa(:)=seg%pmid(:)/(seg%t(:)*rgas)
  seg%vba(:)=SQRT(8.0*boltzmann*seg%t(:)/(pi*ma))
  seg%tsqrt(:)=SQRT(seg%t(:))
  seg%dvisc(:)=1.83e-5*(416.16/(seg%t(:)+120.0))*(SQRT(seg%t(:)/296.16)**3)
  seg%mfpa(:)=2.0*seg%dvisc(:)/(seg%rhoa(:)*seg%vba(:))
  CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,             &
                     mass(1,1,1), seg%sm)   ! mass air (kg/box)

  ! Gas-phase tracers required in aerosol code
  ! ==========================================
  !  S0 array that is passed in to aerosol code uses MH2SO4, msec_org, etc
  !     to index tracers (set in UKCA_SETUP_INDICES module procedures)

  seg%s0(:,:)=0.0                 ! set gas phase tracer masses to zero

  IF (dryox_in_aer == 0) THEN
    seg%s0_dot_condensable(:,:)=0.0 ! condensable tracer tendencies->0
    ! DRYOX_IN_AER=0 -> update of condensables done in UKCA_CHEMISTRY_CTL
  END IF ! if DRYOX_IN_AER=0

  IF (dryox_in_aer == 1) THEN
    ! Initialise S0_DOT_CONDENSABLE to 0
    seg%s0_dot_condensable(:,:)=0.0 ! condensable tracer tendencies->0
    ! DRYOX_IN_AER=1 -> update of condensables done in UKCA_AERO_STEP
    ! condensable tracer tendencies need to be set here to pass in as input

    CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,           &
                       delso2_dry_oh(1,1,1),seg%v1d_tmp )

    IF (ukca_config%l_fix_ukca_h2so4_ystore) THEN
      seg%s0_dot_condensable(:,mh2so4) = seg%v1d_tmp(:)
      ! .. delso2_dry_oh  is in units of vmr/s, correct for
      ! .. S0_DOT_CONDENSABLE
    ELSE
      seg%s0_dot_condensable(:,mh2so4) = seg%v1d_tmp(:)/(seg%aird(:)*dtc)
      ! .. delso2_dry_oh  is in units of molecules/cc/DTC
      ! .. need S0_DOT_CONDENSABLE to be in units of vmr/s
      ! .. so need to divide by (AIRD(:)*DTC)
    END IF

    IF (msec_org > 0) THEN
      seg%s0_dot_condensable(:,msec_org)=0.0
      ! for L_classSO2_inAer=T or F set Sec_Org prodn--> 0 (done in UKCA)
    END IF

    IF (msec_orgi > 0) THEN
      seg%s0_dot_condensable(:,msec_orgi)=0.0
      ! for L_classSO2_inAer=T or F set Sec_Org prodn--> 0 (done in UKCA)
    END IF
    !
  END IF      ! dryox_in_aer
  ! Set H2O2, O3 and SO2 for aqueous oxidation
  IF (wetox_in_aer == 1) THEN
    ! .. set h2o2 when required
    IF (mh2o2f  > 0 .AND. mm_gas(mh2o2f) > 1e-3 ) THEN
      CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,         &
                         chemistry_tracers(1,1,1,n_h2o2), seg%v1d_tmp )
      seg%s0(:,mh2o2f  )= seg%sm(:)*(mm_da/mm_gas(mh2o2f ))*seg%v1d_tmp
    ELSE
      cmessage=' H2O2 needs updating, but MH2O2F'//                            &
                 'or MM_GAS(MH2O2F) is wrong'
      WRITE(umMessage,'(A60,I6,E12.3)') cmessage,mh2o2f,mm_gas(mh2o2f)
      CALL umPrint(umMessage,src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
    END IF

    ! .. set O3
    IF (mox > 0 .AND. mm_gas(mox) > 1e-3 ) THEN
      CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,         &
                         chemistry_tracers(1,1,1,n_o3), seg%v1d_tmp )
      seg%zo3(:) = seg%sm(:)*(mm_da/mm_gas(mox))*seg%v1d_tmp
    ELSE
      cmessage=' O3 needs updating, but MOX or MM_GAS(MOX) is wrong'
      WRITE(umMessage,'(A52,I6,E12.3)') cmessage,' MOX = ',mox,mm_gas(mox)
      CALL umPrint(umMessage,src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
    END IF

    ! .. SO2 mmr in kg[SO2]/kg[dryair]
    IF (msotwo > 0 .AND. mm_gas(msotwo) > 1e-3) THEN
      CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,         &
                         chemistry_tracers(1,1,1,n_so2), seg%v1d_tmp )
      seg%s0(:,msotwo) = seg%sm(:)*(mm_da/mm_gas(msotwo))*seg%v1d_tmp
    ELSE
      cmessage=' SO2 needs updating, but MSOTWO or MM_GAS(MSOTWO)'//           &
               'is wrong'
      WRITE(umMessage,'(A20,I6,E12.3)') cmessage,' MSOTWO = ',msotwo,          &
            mm_gas(msotwo)
      CALL umPrint(umMessage,src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
    END IF

  END IF      ! wetox_in_aer=1

  IF (mh2so4   > 0) THEN
    CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,           &
                       chemistry_tracers(1,1,1,n_h2so4), seg%v1d_tmp )
    seg%s0(:,mh2so4  ) = seg%sm(:)*(mm_da/mm_gas(mh2so4  ))*seg%v1d_tmp(:)
  ELSE
    errcode=1
    cmessage='MH2SO4 <= 0'
    WRITE(umMessage,'(A20,I6)') cmessage,' MH2SO4 = ',mh2so4
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
  END IF

  ! .. Secondary Organic tracer mmr in kg[Sec_Org]/kg[dryair]
  IF (msec_org > 0) THEN
    CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,           &
                       chemistry_tracers(1,1,1,n_sec_org), seg%v1d_tmp )
    seg%s0(:,msec_org) = seg%sm(:)*(mm_da/mm_gas(msec_org))*seg%v1d_tmp
  ELSE
     ! may not have Sec_Org tracer - e.g. StratAer
    errcode=-1
    cmessage='msec_org <= 0'
    WRITE(umMessage,'(A20,I6)') cmessage,' msec_org = ',msec_org
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
  END IF

  ! .. Secondary Organic tracer from isoprene mmr in kg[SEC_ORG I]/kg[dryair]
  IF (msec_orgi > 0) THEN
    CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,           &
                       chemistry_tracers(1,1,1,n_sec_org_i), seg%v1d_tmp )
    seg%s0(:,msec_orgi) = seg%sm(:)*(mm_da/mm_gas(msec_orgi))*seg%v1d_tmp
  END IF


  ! Aerosol Tracers
  ! ===============
  !  Find index of 1st mode tracer, as nmr_index and
  !   mmr_index index all ukca tracers
  ifirst = jpctr + 1
  DO imode=1,nmodes
    seg%nbadmdt(:,imode) = 0
    seg%mdtfixflag(:,imode) = 0.0
    seg%mdtfixsink(:,imode) = 0.0
    IF (mode(imode)) THEN
      itra = nmr_index(imode) - ifirst + 1
      CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,         &
                         mode_tracers(1,1,1,itra), seg%tr_rs)
      seg%mask1(:)=(seg%tr_rs(:) < 0.0)
      WHERE (seg%mask1(:))
        seg%tr_rs(:)=0.0
      END WHERE
      ! .. above sets tr_rs to zero if negative
      seg%nd(:,imode)=seg%tr_rs(:)*seg%aird(:)
      ! .. above sets ND (particles per cc) from advected number-mixing-ratio

      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          itra = mmr_index(imode,icp) - ifirst + 1

          CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,     &
                             mode_tracers(1,1,1,itra), seg%tr_rs)
          seg%mask1(:)=(seg%tr_rs(:) < 0.0)
          WHERE (seg%mask1(:))
            seg%tr_rs(:)=0.0
          END WHERE
          ! .. above sets tr_rs to zero if negative
          seg%mask1(:)=(seg%nd(:,imode) > num_eps(imode))
          WHERE (seg%mask1(:))
            seg%md(:,imode,icp)=(mm_da/mm(icp))                                &
                               *seg%aird(:)*seg%tr_rs(:)/seg%nd(:,imode)
          ELSE WHERE
            seg%md(:,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
          END WHERE
          ! above sets MD (molecules / particle) from advected mass-mixing-ratio
          ! note that only "trusts" advected values where ND>NUM_EPS
        ELSE
          seg%md(:,imode,icp)=0.0
        END IF
      END DO ! loop over cpts
      !
      ! Set total mass array MDT from sum over individual component MDs
      seg%mdt(:,imode)=0.0
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          seg%mdt(:,imode)=seg%mdt(:,imode)+seg%md(:,imode,icp)
        END IF
      END DO
      !
      ! below checks if MDT is coming out too low after advection (af BLMIX)
        !
      IF (iextra_checks == 0) THEN
        mdtmin(imode,ik)=mlo(imode)*0.001 ! set equiv. to DPLIM0*0.1
        !
        DO icp=1,ncp
          IF (component(imode,icp)) THEN
            ! where MDT too low after advection set ND to zero and &
            ! set default MD
            WHERE (seg%mdt(:,imode) < mdtmin(imode,ik))
              seg%md(:,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
            END WHERE
          END IF
        END DO

        ! Count occurrences (NBADMDT) and store percent occurrence and mass-sink
        seg%fac(:)=seg%sm(:)/seg%aird(:)
        ! FAC converts aerosol mass fluxes from kg(dryair)/box/tstep to
        ! moles/gridbox/s

        seg%nbadmdt(:,imode)=0
        seg%mdtfixflag(:,imode)=0.0
        seg%mdtfixsink(:,imode)=0.0
        WHERE (seg%mdt(:,imode) < mdtmin(imode,ik))
          seg%nbadmdt(:,imode)=1
          seg%mdtfixflag(:,imode)=100.0
          ! mdtfixflag enables to track proportion of timesteps that fix is
          ! applied. units of mdtfixsink are moles/s
          seg%mdtfixsink(:,imode)=seg%nd(:,imode)*seg%mdt(:,imode)*            &
                                  seg%fac(:)/mm_da/dtc
          ! mdtfixsink stores total mass removed when fix is applied

        END WHERE
        ! put count of bad mdt cells into a 3D structure for later diagnostic
        CALL int_ins_seg(lb,ncs,seg%nbox_this_seg,stride_s,                    &
                   model_levels,seg%nbadmdt(:,imode), nbadmdt_3d(1,1,1,imode) )

        thread_sum_nbadmdt = SUM(seg%nbadmdt(:,imode))
!$OMP ATOMIC
        sum_nbadmdt(imode) = sum_nbadmdt(imode) + thread_sum_nbadmdt

        ! Set ND->0 where MDT too low (& set MDT->MMID) & count occurrences
        WHERE (seg%mdt(:,imode) < mdtmin(imode,ik))
          seg%nd(:,imode) = 0.0
          seg%mdt(:,imode) = mmid(imode)
          ! where MDT too low after advection (but +ve) set to MMID
        END WHERE
        ! moved reporting nbadmdt to after ik loop over segments
      END IF ! Iextra_checks = 0
    ELSE
      seg%nbadmdt(:,imode)=0
      seg%nd(:,imode)=0.0
      seg%mdt(:,imode)=mmid(imode)
      ! only one thread should set sum_nbadmdt
!$OMP MASTER
      sum_nbadmdt(imode)=0
!$OMP END MASTER
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          seg%md(:,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
        ELSE
          seg%md(:,imode,icp)=0.0
        END IF
      END DO
    END IF
  END DO ! loop over modes
  !
  IF (verbose >= 2) THEN

    DO itra=1,seg%nadvg
      WRITE(umMessage,'(A10,I4,2E12.3)') 'S0 : ',itra,                         &
            MINVAL(seg%s0(:,itra)),                                            &
            MAXVAL(seg%s0(:,itra))
      CALL umPrint(umMessage,src='ukca_aero_ctl')
    END DO

    DO imode=1,nmodes
      IF (mode(imode)) THEN
        WRITE(umMessage,'(A10,I4,2E12.3)') 'ND : ',imode,                      &
            MINVAL(seg%nd (:,imode)),                                          &
            MAXVAL(seg%nd (:,imode))
        CALL umPrint(umMessage,src='ukca_aero_ctl')
        WRITE(umMessage,'(A10,I4,2E12.3)') 'MDT: ',imode,                      &
            MINVAL(seg%mdt(:,imode)),                                          &
            MAXVAL(seg%mdt(:,imode))
        CALL umPrint(umMessage,src='ukca_aero_ctl')
        DO icp=1,ncp
          IF (component(imode,icp)) THEN
            WRITE(umMessage,'(A10,I4,2E12.3)') 'MD : ',imode,                  &
              MINVAL(seg%md(:,imode,icp)),                                     &
              MAXVAL(seg%md(:,imode,icp))
            CALL umPrint(umMessage,src='ukca_aero_ctl')
          END IF
        END DO
      END IF
    END DO
  END IF  ! verbose > 2
  !
  ! .. Apply extra checks for advection artefacts
  IF (iextra_checks > 0) THEN
    CALL ukca_mode_check_artefacts(verbose_local, seg%nbox_this_seg, seg%nd,   &
                  seg%md,                                                      &
                  seg%mdt, seg%sm, seg%aird, mm_da, dtc, seg%mdtfixflag,       &
                  seg%mdtfixsink)
  END IF

  ! .. zero aerosol budget terms before calling UKCA_AERO_STEP
  seg%bud_aer_mas(:,:)=0.0

  ! .. below is call to UKCA_AERO_STEP as at ukca_mode_v1_gm1.f90

  IF (firstcall .AND. verbose > 0) THEN
    WRITE(umMessage,'(A42)') 'Values of input variables passed to mode:'
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'RAINOUT_ON=',rainout_on
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'IMSCAV_ON=',imscav_on
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'WETOX_ON=',wetox_on
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'DDEPAER_ON=',ddepaer_on
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'SEDI_ON=',sedi_on
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'DRYOX_IN AER=',dryox_in_aer
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'WETOX_IN AER=',wetox_in_aer
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'COND_ON=',cond_on
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'NUCL_ON=',nucl_on
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'BLN_ON=',bln_on
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'COAG_ON=',coag_on
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'ICOAG=',icoag
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'IMERGE=',imerge
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'IFUCHS=',ifuchs
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'IDCMFP=',idcmfp
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'ICONDIAM=',icondiam
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'IBLN=',ibln
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'IACTMETHOD=',iactmethod
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'I_NUC_METHOD=',i_nuc_method
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'IDDEPAER=',iddepaer
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'INUCSCAV=',inucscav
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A12,L7)') 'lcvrainout=',lcvrainout
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A23,L7)') 'l_dust_slinn_impc_scav=',l_dust_slinn_impc_scav
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'VERBOSE=',verbose
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'VERBOSE_LOCAL=',verbose_local
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'CHECKMD_ND=',checkmd_nd
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'INTRAOFF=',intraoff
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'INTEROFF=',interoff
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    WRITE(umMessage,'(A15,I4)') 'IDUSTEMS=',idustems
    CALL umPrint(umMessage,src='ukca_aero_ctl')
  END IF
  !
  CALL ukca_aero_step(seg%nbox_this_seg,seg%nchemg,seg%nadvg,seg%nbudaer,      &
   seg%nd,seg%mdt,seg%md,seg%mdwat,seg%s0,seg%drydp,seg%wetdp,                 &
   seg%rhopar,seg%dvol,seg%wvol,seg%sm,seg%aird,seg%airdm3,seg%rhoa,           &
   seg%mfpa,seg%dvisc,seg%t,seg%tsqrt,seg%rh,seg%rh_clr,seg%s,                 &
   seg%pmid,seg%pupper,seg%plower,                                             &
   seg%zo3,seg%zho2,seg%zh2o2,seg%ustr,seg%znotg,                              &
   seg%surtp,seg%craing,seg%draing,seg%craing_up,                              &
   seg%csnowg,seg%dsnowg,seg%fconv_conv,seg%lowcloud,                          &
   seg%vfac,seg%clf,seg%autoconv1d,                                            &
   dtc,dtz,nmts,nzts,seg%lday,act,seg%bud_aer_mas,                             &
   rainout_on,iextra_checks,                                                   &
   imscav_on,wetox_on,ddepaer_on,sedi_on,iso2wetoxbyo3,                        &
   dryox_in_aer,wetox_in_aer,seg%delso2,seg%delso2_2,                          &
   cond_on,nucl_on,coag_on,bln_on,icoag,imerge,ifuchs,                         &
   idcmfp,icondiam,ibln,i_nuc_method,iagecoagnucl67,                           &
   iactmethod,iddepaer,inucscav,lcvrainout,l_dust_slinn_impc_scav,             &
   verbose_local,checkmd_nd,intraoff,                                          &
   interoff,seg%s0_dot_condensable,seg%lwc,seg%clwc,seg%pvol,                  &
   seg%pvol_wat,seg%jlabove,seg%ilscat,seg%n_merge_1d,seg%height,              &
   seg%htpblg)

  !
  ! Update tracers
  ! ==============

  IF (wetox_in_aer > 0) THEN
    IF (mh2o2f > 0) THEN
      ! .. update gas phase H2O2 mmr following SO2 aqueous phase oxidation
      seg%v1d_tmp(:) = (mm_gas(mh2o2f)/mm_da)*seg%s0(:,mh2o2f  )/seg%sm(:)
      CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,          &
                        seg%v1d_tmp(:), chemistry_tracers(1,1,1,n_h2o2))
    END IF
    IF (msotwo > 0) THEN
      ! .. update gas phase SO2 mmr following aqueous phase oxidation
      seg%v1d_tmp(:) = (seg%s0(:,msotwo)/seg%sm(:))*(mm_gas(msotwo)/mm_da)
      CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,          &
                        seg%v1d_tmp(:), chemistry_tracers(1,1,1,n_so2))
    END IF
  END IF

  IF (mh2so4   > 0) THEN
    ! .. update gas phase H2SO4 mmr following H2SO4 condensation/nucleation
    seg%v1d_tmp(:) = (seg%s0(:,mh2so4  )/seg%sm(:))*(mm_gas(mh2so4  )/mm_da)
    CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,            &
                      seg%v1d_tmp(:), chemistry_tracers(1,1,1,n_h2so4))
  END IF

  IF (msec_org > 0) THEN
    ! .. update gas phase Sec_Org mmr following condensation/nucleation
    seg%v1d_tmp(:) = (seg%s0(:,msec_org)/seg%sm(:))*(mm_gas(msec_org)/mm_da)
    CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,            &
                      seg%v1d_tmp(:), chemistry_tracers(1,1,1,n_sec_org))
  END IF

  IF (msec_orgi > 0) THEN
    ! .. update gas phase Sec_Org mmr following condensation/nucleation
    seg%v1d_tmp(:) = (seg%s0(:,msec_orgi)/seg%sm(:))*(mm_gas(msec_orgi)/mm_da)
    CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,            &
                      seg%v1d_tmp(:), chemistry_tracers(1,1,1,n_sec_org_i))
  END IF

  ! Set H2O2, O3 and SO2 for aqueous oxidation
  IF (wetox_in_aer == 1) THEN
    ! .. set h2o2 from h2o2_tracer when required
    IF (mh2o2f  > 0 .AND. mm_gas(mh2o2f) > 1e-3 ) THEN
      CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,         &
                         chemistry_tracers(1,1,1,n_h2o2), seg%v1d_tmp(:) )
      seg%s0(:,mh2o2f  )=seg%sm(:)*(mm_da/mm_gas(mh2o2f ))*seg%v1d_tmp(:)
    ELSE
      cmessage=' H2O2 needs updating, but MH2O2F'//                            &
                 'or MM_GAS(MH2O2F) is wrong'
      WRITE(umMessage,'(A10,I4,E12.3)') cmessage,mh2o2f,mm_gas(mh2o2f)
      CALL umPrint(umMessage,src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
    END IF

    ! .. set O3 from O3_tracer
    IF (mox > 0 .AND. mm_gas(mox) > 1e-3 ) THEN
      CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,        &
                         chemistry_tracers(1,1,1,n_o3),seg%v1d_tmp(:) )
      seg%zo3(:) = seg%sm(:)*(mm_da/mm_gas(mox))*seg%v1d_tmp(:)
    ELSE
      cmessage=' O3 needs updating, but MOX or MM_GAS(MOX) is wrong'
      WRITE(umMessage,'(A10,I4,E12.3)') cmessage,' MOX = ',mox,mm_gas(mox)
      CALL umPrint(umMessage,src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
    END IF

    IF (msotwo > 0 .AND. mm_gas(msotwo) > 1e-3) THEN
      CALL extract_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,         &
                         chemistry_tracers(1,1,1,n_so2),seg%v1d_tmp(:) )
      seg%s0(:,msotwo)=seg%sm(:)*(mm_da/mm_gas(msotwo))*seg%v1d_tmp(:)
    ELSE
      cmessage=' SO2 needs updating, but MSOTWO or MM_GAS(MSOTWO)'//           &
               'is wrong'
      WRITE(umMessage,'(A12,I4,E12.3)') cmessage,' MSOTWO = ',msotwo,          &
            mm_gas(msotwo)
      CALL umPrint(umMessage,src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
    END IF
  END IF      ! wetox_in_aer

  seg%sarea(:,:) = 0.0

  DO imode=1,nmodes
    IF (mode(imode)) THEN
      itra = nmr_index(imode) - ifirst + 1
      ! .. update aerosol no. conc. following aerosol microphysics
      seg%v1d_tmp(:) = (seg%nd(:,imode)/seg%aird(:))
      CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,          &
                        seg%v1d_tmp(:), mode_tracers(1,1,1,itra) )
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          itra = mmr_index(imode,icp) - ifirst + 1
          WHERE (seg%nd(:,imode) <= num_eps(imode))
            seg%md(:,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
          END WHERE
          ! .. update aerosol mmr following aerosol microphysics
          seg%v1d_tmp(:) = (mm(icp)/mm_da)                                     &
                         * (seg%md(:,imode,icp)*seg%nd(:,imode)/seg%aird(:))
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                            seg%v1d_tmp(:),mode_tracers(1,1,1,itra) )
        END IF
      END DO ! loop over cpts

      CALL int_ins_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,         &
                       seg%n_merge_1d(:,imode), n_merge_3d(1,1,1,imode) )
      seg%vconc(:,imode)=seg%wvol(:,imode)*seg%nd(:,imode)
      y(imode)=EXP(2.0*(log_sigmag(imode)**2))

      ! calculate surface area in units of units of cm^2 / cm^3.
      ! nd is aerosol ptcl number density for mode (cm^-3)
      ! wetdp is wet diameter in m
      ! So multiply by 1.0e+4 to convert from m^2 / cm^3 to cm^2 / cm^3
      seg%sarea(:,imode)=pi*1.0e+4*seg%nd(:,imode)                             &
                         *(seg%wetdp(:,imode)**2)*y(imode)

    END IF     ! mode(imode)
  END DO       ! imode=1,nmodes

  !--------------------------------------------------------
  ! below sets CN,CCN,CDN diagnostics

  ! Initialise to zero
  seg%cn_3nm(:)=0.0
  seg%ccn_2 (:)=0.0
  seg%ccn_3 (:)=0.0

  DO imode=1,nmodes
    IF (mode(imode)) THEN
      !
      dp0=3.0e-9
      seg%erf_arg(:)=LOG(dp0/seg%drydp(:,imode))/(root2*log_sigmag(imode))
      seg%erfterm(:)=0.5*seg%nd(:,imode)*(1.0-umErf(seg%erf_arg(:)))
      seg%cn_3nm(:)=seg%cn_3nm(:)+seg%erfterm(:)
      !
      IF (modesol(imode) == 1) THEN
        !
        ! for CCN_2 take CCN for particles > 25nm dry radius
        dp0=50.0e-9
        seg%erf_arg(:)=LOG(dp0/seg%drydp(:,imode))/(root2*log_sigmag(imode))
        seg%erfterm(:)=0.5*seg%nd(:,imode)*(1.0-umErf(seg%erf_arg(:)))
        seg%ccn_2(:)=seg%ccn_2(:)+seg%erfterm(:)
        !
        ! for CCN_3 take CCN for particles > 35nm dry radius
        dp0=70.0e-9
        seg%erf_arg(:)=LOG(dp0/seg%drydp(:,imode))/(root2*log_sigmag(imode))
        seg%erfterm(:)=0.5*seg%nd(:,imode)*(1.0-umErf(seg%erf_arg(:)))
        seg%ccn_3(:)=seg%ccn_3(:)+seg%erfterm(:)
        !
      END IF
    END IF
  END DO

  !----------------------------------------------------
  CALL ukca_cdnc_jones( seg%nbox_this_seg, act, seg%drydp, seg%nd,             &
                        glomap_variables, seg%ccn_1, seg%cdn )
  !----------------------------------------------------

  !
  ! Update DRYDP and DVOL after AERO_STEP and updation of MD
  CALL ukca_calc_drydiam( seg%nbox_this_seg, glomap_variables,                 &
                          seg%nd, seg%md, seg%mdt, seg%drydp, seg%dvol )

  IF (verbose > 1) THEN

    WRITE(umMessage,'(A30)') 'AFTER CALL TO UKCA_AERO_STEP'
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    DO itra=1,seg%nadvg
      WRITE(umMessage,'(A10,I4,2E12.3)') 'S0 : ',itra,                         &
            MINVAL(seg%s0(:,itra)),MAXVAL(seg%s0(:,itra))
      CALL umPrint(umMessage,src='ukca_aero_ctl')
    END DO
    DO imode=1,nmodes
      IF (mode(imode)) THEN
        WRITE(umMessage,'(A10,I4,2E12.3)') 'ND : ',imode,                      &
              MINVAL(seg%nd (:,imode)), MAXVAL(seg%nd (:,imode))
        CALL umPrint(umMessage,src='ukca_aero_ctl')
        WRITE(umMessage,'(A10,I4,2E12.3)') 'MDT: ',imode,                      &
              MINVAL(seg%mdt(:,imode)), MAXVAL(seg%mdt(:,imode))
        CALL umPrint(umMessage,src='ukca_aero_ctl')
        DO icp=1,ncp
          IF (component(imode,icp)) THEN
            WRITE(umMessage,'(A10,I4,2E12.3)') 'MD : ',imode,                  &
                  MINVAL(seg%md(:,imode,icp)), MAXVAL(seg%md(:,imode,icp))
            CALL umPrint(umMessage,src='ukca_aero_ctl')
          END IF
        END DO
      END IF
    END DO

    DO imode=1,nmodes
      IF (mode(imode)) THEN
        WRITE(umMessage,'(A10,I4,2E12.3)') 'DRYDP: ',imode,                    &
               MINVAL(seg%drydp(:,imode)), MAXVAL(seg%drydp(:,imode))
        CALL umPrint(umMessage,src='ukca_aero_ctl')
      END IF
    END DO

  END IF ! if VERBOSE > 1

  ! Calculate heterogeneous rate coeffs for tropospheric chemistry
  IF (ukca_config%l_ukca_trophet) THEN
    CALL ukca_trop_hetchem(seg%nbox_this_seg,seg%nhet,seg%t,seg%rh,seg%aird,   &
                           seg%pvol,seg%wetdp,seg%sarea,seg%het_rates)
    ! Now copy the het_rates into the all_ntp array
    i = name2ntpindex('het_n2o5  ')
    CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,            &
                      seg%het_rates(:,ihet_n2o5), all_ntp(i)%data_3d(1,1,1) )
    i = name2ntpindex('het_ho2   ')
    CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,            &
                      seg%het_rates(:,ihet_ho2_ho2), all_ntp(i)%data_3d(1,1,1))
  END IF

  ! Copy MODE diagnostics from BUD_AER_MAS to BUD_AER_MAS1D
  ! and convert all budget variables to be in kg/DTC

  seg%fac(:)=seg%sm(:)/seg%aird(:)
  ! FAC converts aerosol mass flux from kg(dryair)/box/tstep to moles/gridbox/s
  DO jv=1,seg%nbudaer
    ! NMASCLPR variables are in kg(dryair)/DTC, others in molecules/cc/DTC
    logic1=(jv <  nmasclprsuaitsol1)
    logic2=(jv >  nmasclprsucorsol2)
    logic=logic1 .OR. logic2
    IF (logic) THEN
      ! when passed out from AERO_STEP, BUD_AER_MAS is in molecules/cc/DTC
      seg%bud_aer_mas(:,jv)=seg%bud_aer_mas(:,jv)*seg%fac(:)/mm_da/dtc ! moles/s
    ELSE
      ! when passed out from AERO_STEP, BUD_AER_MAS (NMASCLPR) is
      !                                                   in kg(dryair)/box/DTC
      seg%bud_aer_mas(:,jv)=seg%bud_aer_mas(:,jv)/mm_da/dtc  ! moles/s
    END IF
  END DO

  ! Write 3_D diagnostics to mode_diags array
  !  N.B. L_ukca_mode_diags is set whenever a STASH request for a relevant item
  !  is found, and the ukcaD1codes(N)%item is then set, otherwise it is IMDI

  IF (L_ukca_mode_diags) THEN  ! fill 3D array
    k=0
    DO n=1,nukca_d1items
      IF (ukcaD1codes(n)%section == stashcode_glomap_sec .AND.                 &
          ((ukcaD1codes(n)%item >= item1_mode_diags .AND.                      &
            ukcaD1codes(n)%item <= item1_mode_diags+nmax_mode_diags-1) .OR.    &
           (ukcaD1codes(n)%item >= item1_nitrate_diags .AND.                   &
            ukcaD1codes(n)%item <= itemN_nitrate_diags)) .AND.                 &
          ukcaD1codes(n)%item /= imdi) THEN
        k=k+1

        ! number for user psm  mode_fluxdiagsv6.6_gm3_SUSSBCOC_reduced_gm1
        ! item1_mode_diags=201
        ! prim SU -- 201-203        ! Note that primary emission diagnostics
        ! prim SS -- 204-205        ! are now in routine UKCA_MODE_EMS_UM.
        ! prim BC -- 206-207        !                "
        ! prim OC -- 208-209        !                "
        ! prim DU -- 210-213        !                "
        ! ddep SU -- 214-217
        ! ddep SS -- 218-219
        ! ddep BC -- 220-223
        ! ddep OC -- 224-228
        ! ddep SO -- 229-232
        ! ddep DU -- 233-236
        ! nusc SU -- 237-240        ! These are fluxes from ukca_rainout
        ! nusc SS -- 241-242        !                "
        ! nusc BC -- 243-246        !                "
        ! nusc OC -- 247-251        !                "
        ! nusc SO -- 252-256        !                "
        ! nusc DU -- 257-260        !                "
        ! imsc SU -- 261-264
        ! imsc SS -- 265-266
        ! imsc BC -- 267-270
        ! imsc OC -- 271-275
        ! imsc SO -- 276-279
        ! imsc DU -- 280-283
        ! clpr SU -- 284-289 (this is wet oxidation of SO2)
        ! proc SU,BC,OC,SO -- 290-293 (this is processing Aitsol-->accsol)
        ! cond SU -- 294-300
        ! cond OC -- 301-307
        ! cond SO -- 308-314
        ! htox SU -- 315-318
        ! nucl SU -- 319
        ! coag SU,SS,BC,OC,SO,DU -- 320-370
        ! aged SU,SS,BC,OC,SO,DU -- 371-374
        ! merg SU,SS,BC,OC,SO,DU -- 375-387
        ! drydp1-7 - 401-407
        ! wetdp1-4 - 408-411
        ! mdwat1-4 - 412-415
        ! sarea1-7 - 416-422
        ! vconc1-7 - 423-429
        ! rhop 1-7 - 430-436
        ! cnccncdn - 437-441
        ! pvol,wat - 442-468
        ! ACTIVATE - 469-485
        ! plumescav - 486-496
        ! mdt fix mass loss -- 546-552
        ! mdt fix frequency -- 553-559
        ! ****************** NITRATE :
        ! prim NH -- 574-577
        ! prim NT -- 578-581
        ! cond NN -- 582-583
        ! ddep NH -- 584-587
        ! ddep NT -- 588-591
        ! ddep NN -- 592-593
        ! nusc NH -- 594-597
        ! nusc NT -- 598-601
        ! nusc NN -- 602-603
        ! imsc NH -- 604-607
        ! imsc NT -- 608-611
        ! imsc NN -- 612-613
        ! proc NH -- 614
        ! proc NT -- 615
        ! coag NH -- 616-618,622-623,626
        ! coag NT -- 619-621,624-625,627
        ! coag NN -- 628
        ! merg NH -- 629-631
        ! merg NT -- 632-634
        ! merg NN -- 635
        ! pvol NH -- 636-639
        ! pvol NT -- 640-643
        ! pvol NN -- 644-645

        IF (firstcall .AND. verbose > 0 ) THEN
          WRITE(umMessage,'(A45,3I6)') 'About to set mode diagnostics for'//   &
               ' N,k,item=',n,k,ukcaD1codes(n)%item
          CALL umPrint(umMessage,src='ukca_aero_ctl')
        END IF
        !
        SELECT CASE(UkcaD1codes(n)%item)
        CASE (item1_mode_diags:item1_mode_diags+12)
          !           Do nothing, primary emissions not handled here now
        CASE (item1_mode_diags+13)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepsunucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+14)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepsuaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+15)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepsuaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+16)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepsucorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+17)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepssaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+18)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepsscorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+19)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepbcaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+20)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepbcaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+21)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepbccorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+22)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepbcaitins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+23)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepocnucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+24)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepocaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+25)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepocaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+26)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepoccorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+27)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepocaitins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+28)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepsonucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+29)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepsoaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+30)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepsoaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+31)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepsocorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+32)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepduaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+33)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepducorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+34)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepduaccins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+35)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepducorins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+36)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscsunucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+37)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscsuaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+38)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscsuaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+39)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscsucorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+40)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscssaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+41)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscsscorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+42)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscbcaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+43)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscbcaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+44)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscbccorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+45)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscbcaitins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+46)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscocnucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+47)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscocaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+48)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscocaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+49)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscoccorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+50)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscocaitins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+51)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscsonucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+52)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscsoaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+53)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscsoaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+54)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscsocorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+55)
          mode_diags(:,:,:,k)=0.0
          !! .. there is no MASNUSCSOAITINS --- erroneously included in
          !! .. UKCA_mode stash section so set it to zero here
        CASE (item1_mode_diags+56)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscduaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+57)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscducorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+58)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscduaccins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+59)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscducorins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+60)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscsunucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+61)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscsuaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+62)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscsuaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+63)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscsucorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+64)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscssaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+65)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscsscorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+66)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscbcaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+67)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscbcaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+68)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscbccorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+69)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscbcaitins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+70)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscocnucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+71)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscocaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+72)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscocaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+73)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscoccorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+74)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscocaitins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+75)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscsonucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+76)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscsoaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+77)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscsoaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+78)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscsocorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+79)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscduaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+80)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscducorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+81)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscduaccins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+82)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscducorins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+83)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                     seg%bud_aer_mas(:,nmasclprsuaitsol1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+84)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                     seg%bud_aer_mas(:,nmasclprsuaccsol1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+85)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                     seg%bud_aer_mas(:,nmasclprsucorsol1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+86)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                     seg%bud_aer_mas(:,nmasclprsuaitsol2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+87)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                     seg%bud_aer_mas(:,nmasclprsuaccsol2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+88)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                     seg%bud_aer_mas(:,nmasclprsucorsol2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+89)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasprocsuintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+90)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasprocbcintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+91)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasprococintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+92)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasprocsointr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+93)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsunucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+94)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsuaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+95)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsuaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+96)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsucorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+97)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsuaitins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+98)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsuaccins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+99)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsucorins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+100)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondocnucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+101)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondocaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+102)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondocaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+103)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondoccorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+104)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondocaitins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+105)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondocaccins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+106)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondoccorins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+107)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsonucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+108)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsoaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+109)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsoaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+110)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsocorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+111)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsoaitins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+112)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsoaccins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+113)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondsocorins), mode_diags(1,1,1,k) )
          !!
          !! Code for heterogeneous oxidation of SO2 --> SO4 on dust not yet in
          !!
        CASE (item1_mode_diags+114)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+115)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+116)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+117)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+118)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuclsunucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+119)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsuintr12), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+120)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsuintr13), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+121)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsuintr14), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+122)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsuintr15), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+123)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsuintr16), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+124)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsuintr17), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+125)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagocintr12), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+126)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagocintr13), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+127)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagocintr14), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+128)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagocintr15), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+129)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagocintr16), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+130)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagocintr17), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+131)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsointr12), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+132)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsointr13), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+133)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsointr14), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+134)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsointr15), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+135)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsointr16), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+136)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsointr17), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+137)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsuintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+138)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsuintr24), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+139)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+140)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+141)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagbcintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+142)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagbcintr24), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+143)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+144)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+145)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagocintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+146)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagocintr24), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+147)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+148)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+149)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsointr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+150)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsointr24), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+151)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+152)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+153)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsuintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+154)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+155)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagbcintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+156)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+157)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagocintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+158)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+159)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagssintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+160)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+161)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagsointr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+162)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+163)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagduintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+164)
          mode_diags(:,:,:,k)=0.0
        CASE (item1_mode_diags+165)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagbcintr53), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+166)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagocintr53), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+167)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagbcintr54), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+168)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagocintr54), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+169)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagduintr64), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+170)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasagedsuintr52), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+171)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasagedbcintr52), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+172)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasagedocintr52), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+173)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasagedsointr52), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+174)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergsuintr12), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+175)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergocintr12), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+176)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergsointr12), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+177)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergsuintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+178)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergbcintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+179)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergocintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+180)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergsointr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+181)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergsuintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+182)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergssintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+183)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergbcintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+184)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergocintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+185)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergduintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+186)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergsointr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+187)
          ! Do nothing here, used for marine OM emissions
        CASE (item1_mode_diags+188)
          ! Do nothing here, used for marine OM emissions
          ! SEC_ORG_I condensation
        CASE (item1_mode_diags+189)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondocinucsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+190)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondociaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+191)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondociaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+192)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascondocicorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+193)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                     seg%bud_aer_mas(:,nmascondociaitins), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+200)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%drydp(:,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+201)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%drydp(:,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+202)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%drydp(:,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+203)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%drydp(:,4), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+204)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%drydp(:,5), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+205)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%drydp(:,6), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+206)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%drydp(:,7), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+207)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%wetdp(:,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+208)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%wetdp(:,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+209)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%wetdp(:,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+210)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%wetdp(:,4), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+211)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                  seg%mdwat(:,1)/avogadro, mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+212)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                  seg%mdwat(:,2)/avogadro, mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+213)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                  seg%mdwat(:,3)/avogadro, mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+214)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                  seg%mdwat(:,4)/avogadro, mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+215)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%sarea(:,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+216)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%sarea(:,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+217)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%sarea(:,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+218)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%sarea(:,4), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+219)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%sarea(:,5), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+220)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%sarea(:,6), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+221)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%sarea(:,7), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+222)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%vconc(:,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+223)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%vconc(:,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+224)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%vconc(:,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+225)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%vconc(:,4), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+226)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%vconc(:,5), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+227)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%vconc(:,6), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+228)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                           seg%vconc(:,7), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+229)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%rhopar(:,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+230)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%rhopar(:,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+231)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%rhopar(:,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+232)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%rhopar(:,4), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+233)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%rhopar(:,5), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+234)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%rhopar(:,6), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+235)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%rhopar(:,7), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+236)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                            seg%cn_3nm(:), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+237)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                             seg%ccn_1(:), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+238)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                             seg%ccn_2(:), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+239)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                             seg%ccn_3(:), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+240)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                               seg%cdn(:), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+241)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,1,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+242)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,1,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+243)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,1,6), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+244)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                        seg%pvol_wat(:,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+245)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,2,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+246)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,2,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+247)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,2,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+248)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,2,6), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+249)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                        seg%pvol_wat(:,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+250)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,3,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+251)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,3,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+252)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,3,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+253)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,3,4), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+254)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,3,5), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+255)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,3,6), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+256)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                         seg%pvol_wat(:,3), mode_diags(1,1,1,k))
        CASE (item1_mode_diags+257)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,4,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+258)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,4,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+259)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,4,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+260)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,4,4), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+261)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,4,5), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+262)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,4,6), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+263)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                         seg%pvol_wat(:,4), mode_diags(1,1,1,k))
        CASE (item1_mode_diags+264)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,5,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+265)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,5,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+266)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,6,5), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+267)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                          seg%pvol(:,7,5), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+268:item1_mode_diags+284)
          ! Do nothing - these are used by ukca_activate
        CASE (item1_mode_diags+345)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixsink(:,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+346)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixsink(:,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+347)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixsink(:,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+348)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixsink(:,4), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+349)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixsink(:,5), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+350)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixsink(:,6), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+351)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixsink(:,7), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+352)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixflag(:,1), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+353)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixflag(:,2), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+354)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixflag(:,3), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+355)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixflag(:,4), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+356)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixflag(:,5), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+357)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixflag(:,6), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+358)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                      seg%mdtfixflag(:,7), mode_diags(1,1,1,k) )
          !!
          !!  Nitrate diags
          !!
        CASE (item1_mode_diags+373:item1_mode_diags+380)
          !           Do nothing, NO3/NH4 emissions not handled here now
        CASE (item1_mode_diags+381:item1_mode_diags+382)
          !           Do nothing, NaNO3 emissions not handled here now
        CASE (item1_mode_diags+383)
          !           Do nothing, nucleation mode not on for NH4/NO3
          !           Functionality added as placeholder
        CASE (item1_mode_diags+384)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepnhaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+385)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepnhaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+386)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepnhcorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+387)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+388)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepntaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+389)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepntaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+390)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepntcorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+391)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepnnaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+392)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasddepnncorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+393)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+394)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscnhaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+395)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscnhaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+396)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscnhcorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+397)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+398)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscntaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+399)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscntaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+400)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscntcorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+401)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscnnaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+402)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasnuscnncorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+403)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+404)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscnhaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+405)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscnhaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+406)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscnhcorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+407)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+408)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscntaitsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+409)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscntaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+410)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscntcorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+411)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscnnaccsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+412)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasimscnncorsol), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+413)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasprocnhintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+414)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasprocntintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+415:item1_mode_diags+420)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+421)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagnhintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+422)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagnhintr24), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+423)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagntintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+424)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagntintr24), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+425)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagnhintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+426)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagntintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+427)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmascoagnnintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+428)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+429)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergnhintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+430)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergnhintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+431)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+432)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergntintr23), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+433)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergntintr34), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+434)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                      seg%bud_aer_mas(:,nmasmergnnintr34), mode_diags(1,1,1,k) )
          !!! Partial volumes for RADAER
        CASE (item1_mode_diags+435)
          !           Do nothing, nucleation mode not on for NH4
        CASE (item1_mode_diags+436)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,2,9), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+437)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,3,9), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+438)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,4,9), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+439)
          !           Do nothing, nucleation mode not on for NO3
        CASE (item1_mode_diags+440)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,2,7), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+441)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,3,7), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+442)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,4,7), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+443)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,3,8), mode_diags(1,1,1,k) )
        CASE (item1_mode_diags+444)
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,      &
                                          seg%pvol(:,4,8), mode_diags(1,1,1,k) )
        CASE DEFAULT
          cmessage=' Item not found in CASE statement'
          CALL ereport('UKCA_AERO_CTL',UkcaD1codes(n)%item,cmessage)
        END SELECT
        IF (verbose == 2) THEN
          WRITE(umMessage,'(A16,3i4,3e14.4)') 'UKCA_MODE diag: ',              &
                        k,n,ukcaD1codes(n)%item,                               &
                        SUM(mode_diags(:,:,:,k)),                              &
                     MAXVAL(mode_diags(:,:,:,k)),                              &
                     MINVAL(mode_diags(:,:,:,k))
          CALL umPrint(umMessage,src='ukca_aero_ctl')
        END IF
      END IF    ! section == MODE_diag_sec etc
    END DO      ! nukca_d1items
  END IF        ! L_ukca_mode_diags

  SELECT CASE (glomap_config%i_ukca_activation_scheme)
  CASE ( i_ukca_activation_arg )
    ! Fill the drydiam array as input to UKCA_ACTIVATE
    DO imode=1,nmodes
      IF ( mode(imode) ) THEN
        CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,        &
                        seg%drydp(:,imode), drydiam(1,1,1,imode) )
      END IF
    END DO
  CASE ( i_ukca_activation_jones )
    ! set ukca_cdnc prognostic without activation scheme
    ! change units from cm^-3 to m^-3
    i = name2ntpindex('cdnc      ')
    seg%v1d_tmp = 1.0e+6*seg%cdn
    CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels,            &
                    seg%v1d_tmp, all_ntp(i)%data_3d(1,1,1) )
  END SELECT

  ! Calculate total aerosol surface area for use in UKCA chemistry for
  ! heterogeneous reactions by summing across all soluble modes and
  ! store result in all_ntp structure
  seg%v1d_tmp = 0.0
  DO imode=mode_nuc_sol,mode_cor_sol ! only take the soluble modes
    ! SAREA was converted to cm^2/cm^3 when it was calculated (see above)
    seg%v1d_tmp = seg%v1d_tmp + seg%sarea(:,imode)
  END DO
  i = name2ntpindex('surfarea  ')
  CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,model_levels, seg%v1d_tmp, &
                                             all_ntp(i)%data_3d(1,1,1) )

  ! Fill the non-transported prognostic fields for RADAER coupling,
  ! depending on selected modes and components. The nucleation mode
  ! is not considered.
  ! --------------------------------------------------------------

  IF (glomap_config%l_ukca_radaer) THEN
    errcode = 0
    DO imode = mode_ait_sol,nmodes
      IF (mode(imode)) THEN
        SELECT CASE(imode)
        CASE (mode_ait_sol)
          i = name2ntpindex('drydiam_ait_sol     ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                               seg%drydp(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('wetdiam_ait_sol     ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                               seg%wetdp(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('aerdens_ait_sol     ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                              seg%rhopar(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('pvol_h2o_ait_sol    ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                            seg%pvol_wat(:,imode), all_ntp(i)%data_3d(1,1,1) )
        CASE (mode_acc_sol)
          i = name2ntpindex('drydiam_acc_sol     ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                               seg%drydp(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('wetdiam_acc_sol     ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                               seg%wetdp(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('aerdens_acc_sol     ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                              seg%rhopar(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('pvol_h2o_acc_sol    ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                            seg%pvol_wat(:,imode), all_ntp(i)%data_3d(1,1,1) )
        CASE (mode_cor_sol)
          i = name2ntpindex('drydiam_cor_sol     ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                               seg%drydp(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('wetdiam_cor_sol     ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                               seg%wetdp(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('aerdens_cor_sol     ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                              seg%rhopar(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('pvol_h2o_cor_sol    ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                            seg%pvol_wat(:,imode), all_ntp(i)%data_3d(1,1,1) )
        CASE (mode_ait_insol)
          i = name2ntpindex('drydiam_ait_insol   ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                               seg%drydp(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('aerdens_ait_insol   ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                              seg%rhopar(:,imode), all_ntp(i)%data_3d(1,1,1) )
        CASE (mode_acc_insol)
          i = name2ntpindex('drydiam_acc_insol   ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                               seg%drydp(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('aerdens_acc_insol   ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                              seg%rhopar(:,imode), all_ntp(i)%data_3d(1,1,1) )
        CASE (mode_cor_insol)
          i = name2ntpindex('drydiam_cor_insol   ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                               seg%drydp(:,imode), all_ntp(i)%data_3d(1,1,1) )
          i = name2ntpindex('aerdens_cor_insol   ')
          CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,     &
                              seg%rhopar(:,imode), all_ntp(i)%data_3d(1,1,1) )
        CASE DEFAULT
          cmessage=' Mode not found in RADAER coupling CASE statement'
          errcode = ABS(imode)
          CALL ereport(RoutineName,errcode,cmessage)
        END SELECT
        DO icp=1,ncp
          IF (component(imode,icp)) THEN
            IF (imode == mode_ait_sol) THEN
              SELECT CASE (icp)
              CASE (cp_su)
                i = name2ntpindex('pvol_su_ait_sol    ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                      model_levels,                                            &
                             seg%pvol(:,imode,icp),all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_bc)
                i = name2ntpindex('pvol_bc_ait_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_oc)
                i = name2ntpindex('pvol_oc_ait_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_so)
                i = name2ntpindex('pvol_so_ait_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_no3)
                i = name2ntpindex('pvol_no3_ait_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_nh4)
                i = name2ntpindex('pvol_nh4_ait_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE DEFAULT
                cmessage = ' Component not found in RADAER coupling CASE'//    &
                           ' statement'
                errcode = ABS(imode*100) + ABS(icp)
                CALL ereport(RoutineName,errcode,cmessage)
              END SELECT
            ELSE IF (imode == mode_acc_sol) THEN
              SELECT CASE (icp)
              CASE (cp_su)
                i = name2ntpindex('pvol_su_acc_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_bc)
                i = name2ntpindex('pvol_bc_acc_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_oc)
                i = name2ntpindex('pvol_oc_acc_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_cl)
                i = name2ntpindex('pvol_ss_acc_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_no3)
                i = name2ntpindex('pvol_no3_acc_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_nh4)
                i = name2ntpindex('pvol_nh4_acc_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_nn)
                i = name2ntpindex('pvol_nn_acc_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_du)
                i = name2ntpindex('pvol_du_acc_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_so)
                i = name2ntpindex('pvol_so_acc_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE DEFAULT
                cmessage = ' Component not found in RADAER coupling CASE'//    &
                           ' statement'
                errcode = ABS(imode*100) + ABS(icp)
                CALL ereport(RoutineName,errcode,cmessage)
              END SELECT
            ELSE IF (imode == mode_cor_sol) THEN
              SELECT CASE (icp)
              CASE (cp_su)
                i = name2ntpindex('pvol_su_cor_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                             seg%pvol(:,imode,icp),all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_bc)
                i = name2ntpindex('pvol_bc_cor_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_oc)
                i = name2ntpindex('pvol_oc_cor_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                             seg%pvol(:,imode,icp),all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_cl)
                i = name2ntpindex('pvol_ss_cor_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                             seg%pvol(:,imode,icp),all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_no3)
                i = name2ntpindex('pvol_no3_cor_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                             seg%pvol(:,imode,icp),all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_nh4)
                i = name2ntpindex('pvol_nh4_cor_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                             seg%pvol(:,imode,icp),all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_nn)
                i = name2ntpindex('pvol_nn_cor_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                            seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_du)
                i = name2ntpindex('pvol_du_cor_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                             seg%pvol(:,imode,icp),all_ntp(i)%data_3d(1,1,1) )
              CASE (cp_so)
                i = name2ntpindex('pvol_so_cor_sol     ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                             seg%pvol(:,imode,icp),all_ntp(i)%data_3d(1,1,1) )
              CASE DEFAULT
                cmessage = ' Component not found in RADAER coupling CASE'//    &
                           ' statement'
                errcode = ABS(imode*100) + ABS(icp)
                CALL ereport(RoutineName,errcode,cmessage)
              END SELECT
            ELSE IF (imode == mode_ait_insol) THEN
              SELECT CASE (icp)
              CASE (cp_bc)
                i = name2ntpindex('pvol_bc_ait_insol   ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                             seg% pvol(:,imode,icp),all_ntp(i)%data_3d(1,1,1))
              CASE (cp_oc)
                i = name2ntpindex('pvol_oc_ait_insol   ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                             seg%pvol(:,imode,icp),all_ntp(i)%data_3d(1,1,1) )
              CASE DEFAULT
                cmessage = ' Component not found in RADAER coupling CASE'//    &
                           ' statement'
                errcode = ABS(imode*100) + ABS(icp)
                CALL ereport(RoutineName,errcode,cmessage)
              END SELECT
            ELSE IF (imode == mode_acc_insol) THEN
              SELECT CASE (icp)
              CASE (cp_du)
                i = name2ntpindex('pvol_du_acc_insol   ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                             seg%pvol(:,imode,icp),all_ntp(i)%data_3d(1,1,1) )
              CASE DEFAULT
                cmessage = ' Component not found in RADAER coupling CASE'//    &
                           ' statement'
                errcode = ABS(imode*100) + ABS(icp)
                CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
              END SELECT
            ELSE IF (imode == mode_cor_insol) THEN
              SELECT CASE (icp)
              CASE (cp_du)
                i = name2ntpindex('pvol_du_cor_insol   ')
                CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s,             &
                     model_levels,                                             &
                             seg%pvol(:,imode,icp), all_ntp(i)%data_3d(1,1,1) )
              CASE DEFAULT
                cmessage = ' Component not found in RADAER coupling CASE'//    &
                           ' statement'
                errcode = ABS(imode*100) + ABS(icp)
                CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
              END SELECT
            ELSE
              cmessage = ' mode out of range in RADAER coupling IF clause'
              errcode = ABS(imode)
              CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
            END IF        ! imode == ?
          END IF         ! component
        END DO    ! icp

      END IF    ! mode(imode)
    END DO  ! imode

    IF ( errcode /= 0) THEN
      cmessage = ' Element of all_ntp array uninitialised'
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  END IF   ! l_ukca_radaer

  ! Dry diameter for nucleations-soluble mode is needed when running
  ! ACTIVATE in the Senior component of the hybrid resolution model
  ! (but is not needed by RADAER).
  IF (glomap_config%l_ntpreq_dryd_nuc_sol) THEN
    i = name2ntpindex('drydiam_nuc_sol     ')
    CALL insert_seg(lb,ncs,seg%nbox_this_seg,stride_s, model_levels,           &
                    seg%drydp(:,mode_nuc_sol), all_ntp(i)%data_3d(1,1,1) )

  END IF

  IF (l_ukca_cmip6_diags .OR. l_ukca_pm_diags) THEN
    ! Fill mdwat_diag array
    DO imode=1,nmodes
      IF (mode(imode)) THEN
        CALL insert_seg(lb, ncs, seg%nbox_this_seg, stride_s, model_levels,    &
                        seg%mdwat(:,imode), mdwat_diag(1,imode))
      END IF
    END DO
  END IF

  IF (l_ukca_pm_diags) THEN
    ! Copy wet diameter for calculating PM10 and PM2.5 diagnostics
    DO imode=1,nmodes
      IF (mode(imode)) THEN
        CALL insert_seg(lb, ncs, seg%nbox_this_seg, stride_s, model_levels,    &
                        seg%wetdp(:,imode), wetdp_diag(1,imode))
      END IF
    END DO
  END IF

END DO  ! ik loop over segments

CALL segment_data_deallocate(seg)

!$OMP END PARALLEL
!
! End of OpenMP

! Deallocate lookup tables in dust impaction scavenging routine
IF (l_dust_slinn_impc_scav) CALL ukca_impc_scav_dust_dealloc()

! postponed reporting of nbadmdt from before aero_step
DO imode = 1, nmodes
  IF (mode(imode) ) THEN
    IF ( verbose > 0) THEN
      IF (sum_nbadmdt(imode) > 0 ) THEN
        ! Below print out total occurrences if > 0
        WRITE(umMessage,'(A55)')'MDT<MDTMIN, ND=0:IMODE,MDTMIN,NBADMDT'//      &
                                ' (after bl mix)'
        CALL umPrint(umMessage,src='ukca_aero_ctl')
        IF (verbose > 1) THEN
          WRITE(umMessage,'(I6,E12.3,I12)') imode,MINVAL(mdtmin(imode,:)),     &
                                                             sum_nbadmdt(imode)
          CALL umPrint(umMessage,src='ukca_aero_ctl')
          ! Per level sum of nbadmdt in 4 groups
          jl = model_levels / 4
          j = 0
          DO i = 1, 4
            l = MIN(j+1+jl,model_levels)
            WRITE(umMessage,'((10000I6))')                                     &
                          (SUM(nbadmdt_3d(:,:,k,imode)),k=j+1,l)
            ! Enclosing format in brackets makes it repeat for all items
            CALL umPrint(umMessage,src='ukca_aero_ctl')
            j = j + jl
          END DO
        END IF   ! Verbose > 1
      END IF     ! sum(nbadmdt > 0)
    END IF       ! Verbose > 0
  END IF         ! mode is active
END DO           ! over modes

IF (verbose > 1) THEN

  WRITE(umMessage,'(A30)') ' Tracers at end of UKCA_MODE:'
  CALL umPrint(umMessage,src='ukca_aero_ctl')

  ! Only one model level in UKCA box model
  DO i=1,MIN(2, model_levels)           !model_levels
    DO j=1,n_chemistry_tracers
      WRITE(umMessage,'(A10,I6,A10,I6)') 'Level: ',i,' Tracer: ',j
      CALL umPrint(umMessage,src='ukca_aero_ctl')
      WRITE(umMessage,'(A20,3E12.3)') 'chemistry_tracers:',                    &
         MINVAL(chemistry_tracers(:,:,i,j)),                                   &
         MAXVAL(chemistry_tracers(:,:,i,j)),                                   &
         SUM(chemistry_tracers(:,:,i,j))/REAL(SIZE(chemistry_tracers(:,:,i,j)))
      CALL umPrint(umMessage,src='ukca_aero_ctl')
    END DO
    DO j=1,n_mode_tracers
      IF (mode_tracer_debug(j)) THEN
        WRITE(umMessage,'(A10,I4,A10,I4,A10)') 'Level: ',i,' Tracer: ',        &
               j,mode_tracer_names(j)
        CALL umPrint(umMessage,src='ukca_aero_ctl')
        WRITE(umMessage,'(A14,3E12.3)') 'mode_tracers: ',                      &
               MINVAL(mode_tracers(:,:,i,j)),                                  &
               MAXVAL(mode_tracers(:,:,i,j)),                                  &
               SUM(mode_tracers(:,:,i,j))/REAL(SIZE(mode_tracers(:,:,i,j)))
        CALL umPrint(umMessage,src='ukca_aero_ctl')
      END IF
    END DO
    WRITE(umMessage,'(A30,I6)') 'Number of merges for Level: ',i
    CALL umPrint(umMessage,src='ukca_aero_ctl')
    DO j=1,nmodes
      WRITE(umMessage,'(2I6)') j,SUM(n_merge_3d(:,:,i,j))
      CALL umPrint(umMessage,src='ukca_aero_ctl')
    END DO
  END DO      ! i

  WRITE(umMessage,'(A24,I9,E12.3)')                                            &
       'Total Number of merges=:',                                             &
        SUM(n_merge_3d),REAL(SUM(n_merge_3d))/REAL(SIZE(n_merge_3d))
  CALL umPrint(umMessage,src='ukca_aero_ctl')
  DO j=1,nmodes
    WRITE(umMessage,'(2I9,E12.3)') j,                                          &
         SUM(n_merge_3d(:,:,:,j)),                                             &
    REAL(SUM(n_merge_3d(:,:,:,j)))/REAL(SIZE(n_merge_3d(:,:,:,j)))
    CALL umPrint(umMessage,src='ukca_aero_ctl')
  END DO
END IF      ! verbose

IF (firstcall) firstcall=.FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_aero_ctl
!
SUBROUTINE segment_data_nullify(seg)
! Nullifies pointer arrays in segment.
IMPLICIT NONE

TYPE(segment_data_type), INTENT(IN OUT) :: seg

NULLIFY(seg%aird)
NULLIFY(seg%airdm3)
NULLIFY(seg%autoconv1d)
NULLIFY(seg%bud_aer_mas)
NULLIFY(seg%ccn_1)
NULLIFY(seg%ccn_2)
NULLIFY(seg%ccn_3)
NULLIFY(seg%cdn)
NULLIFY(seg%clf)
NULLIFY(seg%clwc)
NULLIFY(seg%cn_3nm)
NULLIFY(seg%craing)
NULLIFY(seg%craing_up)
NULLIFY(seg%csnowg)
NULLIFY(seg%delso2)
NULLIFY(seg%delso2_2)
NULLIFY(seg%draing)
NULLIFY(seg%drydp)
NULLIFY(seg%dsnowg)
NULLIFY(seg%dvisc)
NULLIFY(seg%dvol)
NULLIFY(seg%erf_arg)
NULLIFY(seg%erfterm)
NULLIFY(seg%fac)
NULLIFY(seg%fac_mmrconv)
NULLIFY(seg%fconv_conv)
NULLIFY(seg%height)
NULLIFY(seg%het_rates)
NULLIFY(seg%htpblg)
NULLIFY(seg%iarr)
NULLIFY(seg%ilscat)
NULLIFY(seg%jlabove)
NULLIFY(seg%karr)
NULLIFY(seg%land_frac)
NULLIFY(seg%larr)
NULLIFY(seg%lday)
NULLIFY(seg%lowcloud)
NULLIFY(seg%lwc)
NULLIFY(seg%mask1)
NULLIFY(seg%md)
NULLIFY(seg%mdt)
NULLIFY(seg%mdtfixflag)
NULLIFY(seg%mdtfixsink)
NULLIFY(seg%mdwat)
NULLIFY(seg%mfpa)
NULLIFY(seg%n_merge_1d)
NULLIFY(seg%nbadmdt)
NULLIFY(seg%nd)
NULLIFY(seg%plower)
NULLIFY(seg%pmid)
NULLIFY(seg%pupper)
NULLIFY(seg%pvol)
NULLIFY(seg%pvol_wat)
NULLIFY(seg%rh)
NULLIFY(seg%rh_clr)
NULLIFY(seg%rhoa)
NULLIFY(seg%rhopar)
NULLIFY(seg%s)
NULLIFY(seg%s0)
NULLIFY(seg%s0_dot_condensable)
NULLIFY(seg%sarea)
NULLIFY(seg%seaice)
NULLIFY(seg%sm)
NULLIFY(seg%surtp)
NULLIFY(seg%t)
NULLIFY(seg%tr_rs)
NULLIFY(seg%tsqrt)
NULLIFY(seg%ustr)
NULLIFY(seg%v1d_tmp)
NULLIFY(seg%vba)
NULLIFY(seg%vconc)
NULLIFY(seg%vfac)
NULLIFY(seg%wetdp)
NULLIFY(seg%wvol)
NULLIFY(seg%zh2o2)
NULLIFY(seg%zho2)
NULLIFY(seg%znotg)
NULLIFY(seg%zo3)

RETURN
END SUBROUTINE segment_data_nullify

SUBROUTINE segment_data_allocate(seg, ncp, nbox, nchemg, nhet, nbudaer, nadvg)
! Allocates segment data, according to the passed array sizes.

USE ukca_mode_setup,  ONLY: nmodes
IMPLICIT NONE

TYPE(segment_data_type), INTENT(IN OUT) :: seg
INTEGER,                 INTENT(IN)    :: ncp
INTEGER,                 INTENT(IN)    :: nbox
INTEGER,                 INTENT(IN)    :: nchemg
INTEGER,                 INTENT(IN)    :: nhet
INTEGER,                 INTENT(IN)    :: nbudaer
INTEGER,                 INTENT(IN)    :: nadvg

seg%nbox_this_seg = nbox
seg%nchemg        = nchemg
seg%nhet          = nhet
seg%nbudaer       = nbudaer
seg%nadvg         = nadvg

ALLOCATE(seg%aird(nbox))
ALLOCATE(seg%airdm3(nbox))
ALLOCATE(seg%autoconv1d(nbox))
ALLOCATE(seg%bud_aer_mas(nbox,0:nbudaer))
ALLOCATE(seg%ccn_1(nbox))
ALLOCATE(seg%ccn_2(nbox))
ALLOCATE(seg%ccn_3(nbox))
ALLOCATE(seg%cdn(nbox))
ALLOCATE(seg%clf(nbox))
ALLOCATE(seg%clwc(nbox))
ALLOCATE(seg%cn_3nm(nbox))
ALLOCATE(seg%craing(nbox))
ALLOCATE(seg%craing_up(nbox))
ALLOCATE(seg%csnowg(nbox))
ALLOCATE(seg%delso2(nbox))
ALLOCATE(seg%delso2_2(nbox))
ALLOCATE(seg%draing(nbox))
ALLOCATE(seg%drydp(nbox,nmodes))
ALLOCATE(seg%dsnowg(nbox))
ALLOCATE(seg%dvisc(nbox))
ALLOCATE(seg%dvol(nbox,nmodes))
ALLOCATE(seg%erf_arg(nbox))
ALLOCATE(seg%erfterm(nbox))
ALLOCATE(seg%fac(nbox))
ALLOCATE(seg%fac_mmrconv(nbox))
ALLOCATE(seg%fconv_conv(nbox))
ALLOCATE(seg%height(nbox))
ALLOCATE(seg%het_rates(nbox,nhet))
ALLOCATE(seg%htpblg(nbox))
ALLOCATE(seg%iarr(nbox))
ALLOCATE(seg%ilscat(nbox))
ALLOCATE(seg%jlabove(nbox))
ALLOCATE(seg%karr(nbox))
ALLOCATE(seg%land_frac(nbox))
ALLOCATE(seg%larr(nbox))
ALLOCATE(seg%lday(nbox))
ALLOCATE(seg%lowcloud(nbox))
ALLOCATE(seg%lwc(nbox))
ALLOCATE(seg%mask1(nbox))
ALLOCATE(seg%md(nbox,nmodes,ncp))
ALLOCATE(seg%mdt(nbox,nmodes))
ALLOCATE(seg%mdtfixflag(nbox,nmodes))
ALLOCATE(seg%mdtfixsink(nbox,nmodes))
ALLOCATE(seg%mdwat(nbox,nmodes))
ALLOCATE(seg%mfpa(nbox))
ALLOCATE(seg%n_merge_1d(nbox,nmodes))
ALLOCATE(seg%nbadmdt(nbox,nmodes))
ALLOCATE(seg%nd(nbox,nmodes))
ALLOCATE(seg%plower(nbox))
ALLOCATE(seg%pmid(nbox))
ALLOCATE(seg%pupper(nbox))
ALLOCATE(seg%pvol(nbox,nmodes,ncp))
ALLOCATE(seg%pvol_wat(nbox,nmodes))
ALLOCATE(seg%rh(nbox))
ALLOCATE(seg%rh_clr(nbox))
ALLOCATE(seg%rhoa(nbox))
ALLOCATE(seg%rhopar(nbox,nmodes))
ALLOCATE(seg%s(nbox))
ALLOCATE(seg%s0(nbox,nadvg))
ALLOCATE(seg%s0_dot_condensable(nbox,nchemg))
ALLOCATE(seg%sarea(nbox,nmodes))
ALLOCATE(seg%seaice(nbox))
ALLOCATE(seg%sm(nbox))
ALLOCATE(seg%surtp(nbox))
ALLOCATE(seg%t(nbox))
ALLOCATE(seg%tr_rs(nbox))
ALLOCATE(seg%tsqrt(nbox))
ALLOCATE(seg%ustr(nbox))
ALLOCATE(seg%v1d_tmp(nbox))
ALLOCATE(seg%vba(nbox))
ALLOCATE(seg%vconc(nbox,nmodes))
ALLOCATE(seg%vfac(nbox))
ALLOCATE(seg%wetdp(nbox,nmodes))
ALLOCATE(seg%wvol(nbox,nmodes))
ALLOCATE(seg%zh2o2(nbox))
ALLOCATE(seg%zho2(nbox))
ALLOCATE(seg%znotg(nbox))
ALLOCATE(seg%zo3(nbox))

RETURN
END SUBROUTINE segment_data_allocate

SUBROUTINE segment_data_deallocate(seg)
! Deallocates segment data, in reverse order to the allocations in
! segment_data_allocate(...)
IMPLICIT NONE

TYPE(segment_data_type), INTENT(IN OUT) :: seg

IF ( ASSOCIATED(seg%zo3) )               DEALLOCATE(seg%zo3 )
IF ( ASSOCIATED(seg%znotg) )             DEALLOCATE(seg%znotg )
IF ( ASSOCIATED(seg%zho2) )              DEALLOCATE(seg%zho2 )
IF ( ASSOCIATED(seg%zh2o2 ) )            DEALLOCATE(seg%zh2o2 )
IF ( ASSOCIATED(seg%wvol) )              DEALLOCATE(seg%wvol )
IF ( ASSOCIATED(seg%wetdp) )             DEALLOCATE(seg%wetdp )
IF ( ASSOCIATED(seg%vfac) )              DEALLOCATE(seg%vfac )
IF ( ASSOCIATED(seg%vconc) )             DEALLOCATE(seg%vconc )
IF ( ASSOCIATED(seg%vba) )               DEALLOCATE(seg%vba )
IF ( ASSOCIATED(seg%v1d_tmp) )           DEALLOCATE(seg%v1d_tmp )
IF ( ASSOCIATED(seg%ustr) )              DEALLOCATE(seg%ustr )
IF ( ASSOCIATED(seg%tsqrt) )             DEALLOCATE(seg%tsqrt )
IF ( ASSOCIATED(seg%tr_rs) )             DEALLOCATE(seg%tr_rs )
IF ( ASSOCIATED(seg%t) )                 DEALLOCATE(seg%t )
IF ( ASSOCIATED(seg%surtp) )             DEALLOCATE(seg%surtp )
IF ( ASSOCIATED(seg%sm) )                DEALLOCATE(seg%sm )
IF ( ASSOCIATED(seg%seaice) )            DEALLOCATE(seg%seaice )
IF ( ASSOCIATED(seg%sarea) )             DEALLOCATE(seg%sarea )
IF ( ASSOCIATED(seg%s0_dot_condensable)) DEALLOCATE(seg%s0_dot_condensable )
IF ( ASSOCIATED(seg%s0) )                DEALLOCATE(seg%s0 )
IF ( ASSOCIATED(seg%s) )                 DEALLOCATE(seg%s )
IF ( ASSOCIATED(seg%rhopar) )            DEALLOCATE(seg%rhopar )
IF ( ASSOCIATED(seg%rhoa) )              DEALLOCATE(seg%rhoa )
IF ( ASSOCIATED(seg%rh_clr) )            DEALLOCATE(seg%rh_clr )
IF ( ASSOCIATED(seg%rh ) )               DEALLOCATE(seg%rh )
IF ( ASSOCIATED(seg%pvol_wat) )          DEALLOCATE(seg%pvol_wat )
IF ( ASSOCIATED(seg%pvol) )              DEALLOCATE(seg%pvol )
IF ( ASSOCIATED(seg%pupper) )            DEALLOCATE(seg%pupper )
IF ( ASSOCIATED(seg%pmid) )              DEALLOCATE(seg%pmid )
IF ( ASSOCIATED(seg%plower) )            DEALLOCATE(seg%plower )
IF ( ASSOCIATED(seg%nd) )                DEALLOCATE(seg%nd )
IF ( ASSOCIATED(seg%nbadmdt) )           DEALLOCATE(seg%nbadmdt )
IF ( ASSOCIATED(seg%n_merge_1d) )        DEALLOCATE(seg%n_merge_1d )
IF ( ASSOCIATED(seg%mfpa) )              DEALLOCATE(seg%mfpa )
IF ( ASSOCIATED(seg%mdwat) )             DEALLOCATE(seg%mdwat )
IF ( ASSOCIATED(seg%mdtfixsink) )        DEALLOCATE(seg%mdtfixsink )
IF ( ASSOCIATED(seg%mdtfixflag) )        DEALLOCATE(seg%mdtfixflag )
IF ( ASSOCIATED(seg%mdt) )               DEALLOCATE(seg%mdt )
IF ( ASSOCIATED(seg%md) )                DEALLOCATE(seg%md )
IF ( ASSOCIATED(seg%mask1) )             DEALLOCATE(seg%mask1 )
IF ( ASSOCIATED(seg%lwc) )               DEALLOCATE(seg%lwc )
IF ( ASSOCIATED(seg%lowcloud) )          DEALLOCATE(seg%lowcloud )
IF ( ASSOCIATED(seg%lday) )              DEALLOCATE(seg%lday )
IF ( ASSOCIATED(seg%larr) )              DEALLOCATE(seg%larr )
IF ( ASSOCIATED(seg%land_frac) )         DEALLOCATE(seg%land_frac )
IF ( ASSOCIATED(seg%karr) )              DEALLOCATE(seg%karr )
IF ( ASSOCIATED(seg%jlabove) )           DEALLOCATE(seg%jlabove )
IF ( ASSOCIATED(seg%ilscat) )            DEALLOCATE(seg%ilscat )
IF ( ASSOCIATED(seg%iarr) )              DEALLOCATE(seg%iarr )
IF ( ASSOCIATED(seg%htpblg) )            DEALLOCATE(seg%htpblg )
IF ( ASSOCIATED(seg%het_rates) )         DEALLOCATE(seg%het_rates )
IF ( ASSOCIATED(seg%height) )            DEALLOCATE(seg%height )
IF ( ASSOCIATED(seg%fconv_conv) )        DEALLOCATE(seg%fconv_conv )
IF ( ASSOCIATED(seg%fac_mmrconv) )       DEALLOCATE(seg%fac_mmrconv )
IF ( ASSOCIATED(seg%fac) )               DEALLOCATE(seg%fac )
IF ( ASSOCIATED(seg%erfterm) )           DEALLOCATE(seg%erfterm )
IF ( ASSOCIATED(seg%erf_arg) )           DEALLOCATE(seg%erf_arg )
IF ( ASSOCIATED(seg%dvol) )              DEALLOCATE(seg%dvol )
IF ( ASSOCIATED(seg%dvisc) )             DEALLOCATE(seg%dvisc )
IF ( ASSOCIATED(seg%dsnowg) )            DEALLOCATE(seg%dsnowg )
IF ( ASSOCIATED(seg%drydp) )             DEALLOCATE(seg%drydp )
IF ( ASSOCIATED(seg%draing) )            DEALLOCATE(seg%draing )
IF ( ASSOCIATED(seg%delso2_2) )          DEALLOCATE(seg%delso2_2 )
IF ( ASSOCIATED(seg%delso2) )            DEALLOCATE(seg%delso2 )
IF ( ASSOCIATED(seg%csnowg) )            DEALLOCATE(seg%csnowg )
IF ( ASSOCIATED(seg%craing_up) )         DEALLOCATE(seg%craing_up )
IF ( ASSOCIATED(seg%craing) )            DEALLOCATE(seg%craing )
IF ( ASSOCIATED(seg%cn_3nm) )            DEALLOCATE(seg%cn_3nm )
IF ( ASSOCIATED(seg%clwc) )              DEALLOCATE(seg%clwc )
IF ( ASSOCIATED(seg%clf) )               DEALLOCATE(seg%clf )
IF ( ASSOCIATED(seg%cdn) )               DEALLOCATE(seg%cdn )
IF ( ASSOCIATED(seg%ccn_3) )             DEALLOCATE(seg%ccn_3 )
IF ( ASSOCIATED(seg%ccn_2) )             DEALLOCATE(seg%ccn_2 )
IF ( ASSOCIATED(seg%ccn_1) )             DEALLOCATE(seg%ccn_1 )
IF ( ASSOCIATED(seg%bud_aer_mas) )       DEALLOCATE(seg%bud_aer_mas )
IF ( ASSOCIATED(seg%autoconv1d) )        DEALLOCATE(seg%autoconv1d )
IF ( ASSOCIATED(seg%airdm3) )            DEALLOCATE(seg%airdm3 )
IF ( ASSOCIATED(seg%aird) )              DEALLOCATE(seg%aird )

CALL segment_data_nullify(seg)

RETURN
END SUBROUTINE segment_data_deallocate

SUBROUTINE surface_to_seg (lbase,ncol,nb,stride,nl, a, b)
! This replicates the surface value up a column
IMPLICIT NONE
INTEGER,INTENT(IN) :: lbase   ! Address of base of first column in segment
INTEGER,INTENT(IN) :: ncol    ! Number of columns in this segment
INTEGER,INTENT(IN) :: nb      ! The number of boxes in this segment
INTEGER,INTENT(IN) :: stride  ! The number of columns on the MPI task
INTEGER,INTENT(IN) :: nl      ! The number of model levels
REAL,INTENT(IN) :: a(stride)  ! The 2D surface data from which a column is taken
REAL,INTENT(IN OUT) :: b(nb)   ! Segment made up from columns

! local loop iterators
INTEGER :: ic, l, ia, ib

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SURFACE_TO_SEG'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
ia = lbase       ! the address of base of 1st column on segment
ib = 1

DO ic = 1, ncol           ! loop over columns on the segment
  DO l = 1, nl            ! climb a column
    b(ib) = a(ia)         ! copy the data into the shorter vector
    ib = ib + 1           ! next location in segment
  END DO
  ia = lbase+ic           ! start of base of next column
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE surface_to_seg
!
SUBROUTINE surface_to_seg_int (lbase,ncol,nb,stride,nl, a, b)
! This replicates the surface value up a column for integer arrays

IMPLICIT NONE
INTEGER,INTENT(IN) :: lbase   ! Address of base of first column in segment
INTEGER,INTENT(IN) :: ncol    ! Number of columns in this segment
INTEGER,INTENT(IN) :: nb      ! The number of boxes in this segment
INTEGER,INTENT(IN) :: stride  ! The number of columns on the MPI task
INTEGER,INTENT(IN) :: nl      ! The number of model levels
INTEGER,INTENT(IN) :: a(stride)  ! The 2D surface var from which column is taken
INTEGER,INTENT(IN OUT) :: b(nb)  ! Segment made up from columns

! local loop iterators
INTEGER :: ic, l, ia, ib

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SURFACE_TO_SEG_INT'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
ia = lbase       ! the address of base of 1st column on segment
ib = 1

DO ic = 1, ncol           ! loop over columns on the segment
  DO l = 1, nl            ! climb a column
    b(ib) = a(ia)         ! copy the data into the shorter vector
    ib = ib + 1           ! next location in segment
  END DO
  ia = lbase+ic           ! start of base of next column
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE surface_to_seg_int
!
SUBROUTINE extract_seg (lbase,ncol,nb,stride,nl, a, b)
! This puts the elements of a 3D array and packs ito a segment (of columns)
IMPLICIT NONE
INTEGER,INTENT(IN) :: lbase   ! Address of base of first column in segment
INTEGER,INTENT(IN) :: ncol    ! Number of columns in this segment
INTEGER,INTENT(IN) :: nb      ! The number of boxes in this segment
INTEGER,INTENT(IN) :: stride  ! The number of columns on the MPI task
INTEGER,INTENT(IN) :: nl      ! The number of model levels
REAL,INTENT(IN) :: a(stride*nl)    ! the 3D data from which a column is taken
REAL,INTENT(IN OUT) :: b(nb)        ! segment made up from columns

! local loop iterators
INTEGER :: ic, l, ia, ib

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXTRACT_SEG'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
ia = lbase                ! the address of base of 1st column on segment
ib = 1

DO ic = 1, ncol           ! loop over columns on the segment
  DO l = 1, nl            ! climb a column
    b(ib) = a(ia)         ! copy the data into the shorter vector
    ib = ib + 1           ! next location in segment
    ia = ia + stride      ! next location in unrolled 3D array
  END DO
  ia = lbase+ic           ! start of base of next column
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE extract_seg
!
SUBROUTINE insert_seg (lbase,ncol,nb,stride,nl, b, a)
! This subroutine puts the elements of the segment vector back into the 3D array
IMPLICIT NONE
INTEGER,INTENT(IN) :: lbase   ! Address of base of first column in segment
INTEGER,INTENT(IN) :: ncol    ! Number of columns in this segment
INTEGER,INTENT(IN) :: nb      ! The number of boxes in this segment
INTEGER,INTENT(IN) :: stride  ! The number of columns on the MPI task
INTEGER,INTENT(IN) :: nl      ! The number of model levels
REAL,INTENT(IN)    :: b(nb)   ! the incoming segment
REAL,INTENT(IN OUT) :: a(stride*nl) ! 3D array to be modified by b()

! local loop iterators
INTEGER :: ic, l, ia, ib

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INSERT_SEG'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


ia = lbase                ! the address of base of 1st column on segment
ib = 1                    ! the box id in the segment

DO ic = 1, ncol           ! loop over columns on the segment
  DO l = 1, nl            ! climb a column
    a(ia) = b(ib)         ! copy (spreading) the data into the 3d array
    ib = ib + 1           ! next location in segment
    ia = ia + stride      ! next location in unrolled 3D array
  END DO
  ia = lbase+ic           ! start of base of next column
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE insert_seg
!
SUBROUTINE int_ins_seg (lbase,ncol,nb,stride,nl, b, a)
! This subroutine puts the elements of the segment vector back into the 3D array
IMPLICIT NONE
INTEGER,INTENT(IN) :: lbase   ! Address of base of first column in segment
INTEGER,INTENT(IN) :: ncol    ! Number of columns in this segment
INTEGER,INTENT(IN) :: nb      ! The number of boxes in this segment
INTEGER,INTENT(IN) :: stride  ! The number of columns on the MPI task
INTEGER,INTENT(IN) :: nl      ! The number of model levels
INTEGER(KIND=integer_32) ,INTENT(IN)    :: b(nb)        ! incoming segment
INTEGER(KIND=integer_32) ,INTENT(IN OUT) :: a(stride*nl) ! 3D array unwound

! local loop iterators
INTEGER :: ic, l, ia, ib
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INT_INS_SEG'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ia = lbase                ! the address of base of 1st column on segment
ib = 1                    ! the box id in the segment

DO ic = 1, ncol           ! loop over columns on the segment
  DO l = 1, nl            ! climb a column
    a(ia) = b(ib)         ! copy (spreading) the data into the 3d array
    ib = ib + 1           ! next location in segment
    ia = ia + stride      ! next location in unrolled 3D array
  END DO
  ia = lbase+ic           ! start of base of next column
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE int_ins_seg
!
END MODULE ukca_aero_ctl_mod
