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
                              item1_dust3mode_diags, itemN_dust3mode_diags,    &
                              log_v, rgas => r
USE ukca_types_mod, ONLY: log_small, integer_32
USE parkind1,   ONLY: jprb, jpim
USE yomhook,    ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_AERO_CTL_MOD'

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
                all_tracer_names,                                              &
                all_tracers,                                                   &
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
                            mode_cor_insol, mode_sup_insol,                    &
                            cp_su, cp_bc, cp_oc, cp_so, cp_cl, cp_du,          &
                            cp_no3, cp_nh4, cp_nn

USE ukca_mode_tracer_maps_mod, ONLY:  nmr_index, mmr_index
USE ukca_mode_verbose_mod, ONLY: verbose => glob_verbose
USE chemistry_constants_mod, ONLY: avogadro, boltzmann
USE ukca_mode_check_artefacts_mod, ONLY:                                       &
                            ukca_mode_check_artefacts
USE ukca_cspecies,    ONLY: n_h2so4,n_h2o2,n_so2,n_o3,n_sec_org, n_sec_org_i,  &
                            n_hono2,n_nh3
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
                              mhno3, mnh3, nadvg, nchemg, ichem,               &
                              nmasagedsuintr52, nmasagedbcintr52,              &
                              nmasagedocintr52, nmasagedsointr52,              &
                              nmasagedduintr63, nmasagedduintr74,              &
                              nmasagedduintr84,                                &
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
                              nmasprimntaitsol, nmasprimntaccsol,              &
                              nmasprimntcorsol, nmasprimnhaitsol,              &
                              nmasprimnhaccsol, nmasprimnhcorsol,              &
                              nmascondnnaccsol, nmascondnncorsol,              &
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
                              nmasprocntintr23, nmasprocnhintr23,              &
                              ! 3-mode dust indices
                              nmasddepdusupins, nmasnuscdusupins,              &
                              nmasimscdusupins, nmascondsusupins,              &
                              nmascondocsupins, nmascondsosupins,              &
                              nmascoagsuintr18, nmascoagocintr18,              &
                              nmascoagsointr18

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

REAL, INTENT(IN) :: dtc                                  ! timestep(s)
REAL, INTENT(IN) :: pres(row_length, rows, model_levels) ! pressure
REAL, INTENT(IN) :: temp(row_length, rows, model_levels) ! temperature
REAL, INTENT(IN) :: q(row_length, rows, model_levels)    ! sp humidity
REAL, INTENT(IN) :: rh3d(row_length, rows, model_levels) ! rh (frac)
REAL, INTENT(IN) :: rh3d_clr(row_length, rows, model_levels)
! rh (frac) - clear sky portion
REAL, INTENT(IN) :: p_bdrs(row_length, rows, 0:model_levels)
! pressure on interfaces
REAL, INTENT(IN) :: sea_ice_frac(row_length, rows)     ! sea ice
REAL, INTENT(IN) :: u_s(row_length, rows)              ! friction velocity
REAL, INTENT(IN) :: z0m(row_length, rows)              ! roughness length
REAL, INTENT(IN) :: drain(row_length, rows, model_levels) ! 3-D LS rain rate
REAL, INTENT(IN) :: crain(row_length, rows, model_levels) ! 3-D conv rain
REAL, INTENT(IN) :: dsnow(row_length, rows, model_levels)
! 3-D LS snowfall rate
REAL, INTENT(IN) :: csnow(row_length, rows, model_levels)
! 3-D conv snow rate
REAL, INTENT(IN) :: autoconv(row_length, rows, model_levels)
! Autoconversion rate kg/kg/s
REAL, INTENT(IN) :: accretion(row_length, rows, model_levels)
! Accretion rate kg/kg/s
REAL, INTENT(IN) :: rim_agg(row_length, rows, model_levels)
! Riming rate of aggregates kg/kg/s
REAL, INTENT(IN) :: rim_cry(row_length, rows, model_levels)
! Riming rate of ice crystals kg/kg/s
REAL, INTENT(IN) :: land_fraction(row_length, rows)     ! land_fraction
! in-cloud oxidation rates (molecules/cc/DTC) from h2o2 & o3 (UKCA):
REAL, INTENT(IN) :: delso2_wet_h2o2(row_length, rows, model_levels)
REAL, INTENT(IN) :: delso2_wet_o3  (row_length, rows, model_levels)
! in-air   oxidation rate  (molecules/cc/DTC) from oh        (UKCA):
! Note: When l_fix_ukca_h2so4_ystore=T this in (vmr/s).
!       For ASAD-based chemical schemes (e.g. StratTrop)
!       this is the change in H2SO4 from chemistry.
REAL, INTENT(IN) :: delso2_dry_oh   (row_length, rows, model_levels)
! in-cloud oxidation rates (kgS/kgair/s     ) from h2o2 & o3 (CLASSIC):
!      REAL, INTENT(IN) :: delso2_wet_h2o2C(row_length, rows, model_levels)
!      REAL, INTENT(IN) :: delso2_wet_o3C  (row_length, rows, model_levels)
! in-air   oxidation rate  (kgS/kgair/s     ) from oh        (CLASSIC):
!      REAL, INTENT(IN) :: delso2_dry_ohC  (row_length, rows, model_levels)
! cloud fraction
REAL, INTENT(IN) :: cloud_frac(row_length, rows, model_levels)
REAL, INTENT(IN) :: cloud_liq_frac(row_length, rows, model_levels)
REAL, INTENT(IN) :: cloud_liq_wat(row_length, rows, model_levels)
REAL, INTENT(IN) :: mass(row_length, rows, model_levels)
REAL, INTENT(IN) :: zbl(row_length, rows)  ! BL height

! names of tracers
CHARACTER(LEN=maxlen_fieldname), INTENT(IN) ::                                 &
                          all_tracer_names(n_chemistry_tracers+ n_mode_tracers)

! tracer mass mixing ratios
REAL, INTENT(IN OUT) :: all_tracers(row_length, rows, model_levels,            &
                                    n_chemistry_tracers+ n_mode_tracers)

! 3-D diagnostic array
REAL, INTENT(IN OUT) :: mode_diags(row_length, rows,                           &
                                   model_levels, n_mode_diags)

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

REAL :: a3d_tmp(row_length, rows, model_levels) ! temporary 3D array re-used

INTEGER, PARAMETER :: nhet = 2
! Number of heterogeneous reaction rates

INTEGER, PARAMETER :: nmts = 1
! No. of microphysical sub-steps per DTC
INTEGER :: nzts
! No. of condensation-nucleation competition sub-steps per DTM
INTEGER :: iextra_checks
! Level of protection against bad MDT values following advection
INTEGER :: rainout_on
! Switch for whether rainout (nucl. scav.) is on/off
INTEGER :: imscav_on
! Switch for whether impaction scavenging is on/off
INTEGER :: wetox_on
! Switch for whether wet oxidation (cloud processing) is on/off
INTEGER :: ddepaer_on
! Switch for whether aerosol dry deposition is on/off
INTEGER, PARAMETER :: sedi_on = 1
! Switch for whether aerosol sedimentation is on/off
INTEGER, PARAMETER :: iso2wetoxbyo3 = 1
! Switch for whether SO2 wet oxidation by ozone is on/off
! Note that this switch is only used if WETOX_IN_AER=1
! When code used in UM    , WETOX_IN_AER is always set to 0
! When code used in TOMCAT, WETOX_IN_AER is always set to 1
INTEGER :: cond_on
! Switch for whether vapour condensation is  on/off
INTEGER :: nucl_on
! Switch for whether binary nucleation is on/off
INTEGER :: bln_on
! Switch for whether binary BL nucleation is on/off
INTEGER :: fine_no3_prod_on
! Switch to determine whether fine NO3 production is on/off
INTEGER :: coarse_no3_prod_on
! Switch to determine whether coarse NO3 production is on/off
INTEGER :: coag_on
! Switch for whether coagulation is on/off
INTEGER, PARAMETER :: icoag = 1
! Switch for KIJ method (1:GLOMAP, 2: M7, 3: UMorig, 4:UMorig MFPP)
!   =3 Cunnigham scheme as in UM, =4 as in UM but computing values)
INTEGER, PARAMETER :: imerge = 2
! Switch to use mid-pts (=1), edges (2) or dynamic (=3) in remode
INTEGER, PARAMETER :: ifuchs = 2
! Switch for Fuchs(1964) (=1) or Fuchs-Sutugin(1971) for CC (=2)
INTEGER, PARAMETER :: idcmfp = 2
! Switch for vapour-diffusion-method (1=as bin v1, 2=as bin v1.1)
INTEGER, PARAMETER :: icondiam = 2
! Switch for what diameter to use for CONDEN (1=g.m.diam, 2=conden-diam)
INTEGER, PARAMETER :: i_nuc_method = 2
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
INTEGER, PARAMETER :: iddepaer = 2
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
INTEGER, PARAMETER :: checkmd_nd = 1
! Switch for whether to check for bad values of MD and ND
INTEGER, PARAMETER :: intraoff = 0
! Switch to turn off intra-modal coagulation
INTEGER, PARAMETER :: interoff = 0
! Switch to turn off inter-modal coagulation
INTEGER, PARAMETER :: idustems = 0
! Switch for using Pringle scheme (=1) or AEROCOMdaily (=2)

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

! HNO3 uptake coefficient
REAL :: hno3_uptake_coeff

! Items for MDT too low/hi check
REAL :: mdtmin(nmodes, nseg)

INTEGER :: field_size              ! size of 2D field
INTEGER :: field_size3d            ! size of 3D field

INTEGER :: i, j, k, l, n, jl, jl2, imode, icp
! No of tracers required
INTEGER :: n_reqd_tracers
INTEGER :: itra, iaer

! Scaling factor for the in-cloud production rates delso2,delso2_2
! to account for removal by precipitation before the cloud
! evaporates.
REAL :: scale_delso2

! Mass of air molecule (kg)
REAL, PARAMETER :: ma = 4.78e-26
! radius for activation (m)
REAL :: act


! Switch for convective rainout (.FALSE. if done with convective transport)
LOGICAL :: lcvrainout
! Switch to turn on the new impaction scavenging scheme for dust
LOGICAL :: l_dust_slinn_impc_scav

! Used for debug output
LOGICAL (KIND=log_small), ALLOCATABLE, SAVE :: mode_tracer_debug(:)

! Error message
CHARACTER(LEN=errormessagelength) :: cmessage

! For use in storing of aerosol budget terms
LOGICAL :: logic
LOGICAL :: logic1
LOGICAL :: logic2

LOGICAL, SAVE :: firstcall=.TRUE.
! counter: mode-merges applied
INTEGER (KIND=integer_32) :: n_merge_3d(row_length, rows, model_levels, nmodes)
INTEGER (KIND=integer_32) :: sum_nbadmdt(nmodes)
INTEGER (KIND=integer_32) :: nbadmdt_3d(row_length, rows, model_levels, nmodes)
INTEGER :: jv
INTEGER :: ifirst    ! index of first mode tracer in nmr_index, mmr_index

! This taken out of run_ukca as set here for now
INTEGER :: i_mode_act_method

REAL :: root2                ! square root of 2
REAL :: log_sigmag(nmodes)

! Start and end positions of segments.
INTEGER :: i_start, i_start_cp, i_end, i_end_cp
! Used for indexing
INTEGER :: nbs_index(0:nmodes)

! INTEGER segments used in the OpenMP region.
! Index of LAT for grid box
INTEGER :: seg_iarr(nbox)
! Land surface category (based on 9 UM landsurf types)
INTEGER :: seg_ilscat(nbox)
! Index of box directly above this grid box
INTEGER :: seg_jlabove(nbox)
! Index of LON for grid box
INTEGER :: seg_karr(nbox)
! Index of vertical level for grid box
INTEGER :: seg_larr(nbox)
! Switch for day/night (1/0)
INTEGER :: seg_lday(nbox)
INTEGER (KIND=integer_32) :: seg_n_merge_1d(nbox * nmodes)
INTEGER (KIND=integer_32) :: seg_nbadmdt(nbox * nmodes)

! REAL segments used in the OpenMP region
! Grid box mass of air (kg)
REAL :: seg_aird(nbox)
! Number density of air (per cm3)
REAL :: seg_airdm3(nbox)
! Autoconversion rate (including accretion) together
! with snow and ice melting rate (kg.kg^-1.s^-1)
REAL :: seg_autoconv1d(nbox)
! Outputs for budget calculations
REAL :: seg_bud_aer_mas(nbox * (nbudaer + 1))
! CCN concentration (acc-sol + cor-sol)
REAL :: seg_ccn_1(nbox)
! CCN concentration (acc-sol + cor-sol + Aitsol>25nm drydp)
REAL :: seg_ccn_2(nbox)
! CCN concentration (acc-sol + cor-sol + Aitsol>35nm drydp)
REAL :: seg_ccn_3(nbox)
! CCN concentration (acc-sol + cor-sol + Aitsol/ins>15nm drydp)
REAL :: seg_ccn_4(nbox)
! CCN concentration (acc-sol + cor-sol + Aitsol/ins>25nm drydp)
REAL :: seg_ccn_5(nbox)
! CDN concentration
REAL :: seg_cdn(nbox)
! Liquid Cloud fraction
REAL :: seg_clf(nbox)
! Cloud liquid water content [kg/kg]
REAL :: seg_clwc(nbox)
REAL :: seg_cn_3nm(nbox)
! Rain rate for conv precip. in box (kgm^-2s^-1)
REAL :: seg_craing(nbox)
! Rain rate for conv precip. in box above (kgm^-2s^-1)
REAL :: seg_craing_up(nbox)
! Rain rate for conv snow. in box (kgm^-2s^-1)
REAL :: seg_csnowg(nbox)
! S(IV) --> S(VI) by H2O2 (molecules per cc) [input if WETOX_IN_AER=0]
REAL :: seg_delso2(nbox)
! S(IV) --> S(VI) by O3   (molecules per cc) [input if WETOX_IN_AER=0]
REAL :: seg_delso2_2(nbox)
! Rain rate for dyn. precip. in box (kgm^-2s^-1)
REAL :: seg_draing(nbox)
! Geometric mean dry diameter of particles in each mode (m)
REAL :: seg_drydp(nbox * nmodes)
! Rain rate for dyn. precip. in box (kgm^-2s^-1)
REAL :: seg_dsnowg(nbox)
! Dynamic viscosity of air (kg m^-1 s^-1)
REAL :: seg_dvisc(nbox)
! Geometric mean dry volume of particles in each mode (m^3)
REAL :: seg_dvol(nbox * nmodes)

! Items for CN,CCN,CDN calculation
REAL :: seg_erf_arg(nbox)
! Error Fn argument
REAL :: seg_erfterm(nbox)

! Conversion factor for budget diagnostics
REAL :: seg_fac(nbox)
! Fraction of box condensate --> rain in 6 hours (conv)
REAL :: seg_fconv_conv(nbox)
! Mid-level height of gridbox
REAL :: seg_height(nbox)
! Diagnostic to hold heterogeneous rates for tropospheric chemistry
REAL :: seg_het_rates(nbox * nhet)
! Height of boundary-layer in gridbox vertical-column
REAL :: seg_htpblg(nbox)
! Fraction of horizontal gridbox area covered by land
REAL :: seg_land_frac(nbox)
! Horizontal low cloud fraction
REAL :: seg_lowcloud(nbox)
! Cloud liquid water content [kg/m3]
REAL :: seg_lwc(nbox)
! Avg cpt mass of aerosol particle in mode (particle^-1)
REAL :: seg_md(nbox * nmodes * glomap_variables%ncp)
! Avg tot mass of aerosol ptcl in mode (particle^-1)
REAL :: seg_mdt(nbox * nmodes)

! Arrays storing info about how much mass is removed when ND>0 for o-o-r MDT
! mdtfixflag gets to set = 100.0 when fix applied so that means show the
! percentage of timesteps on which the fix is applied.
REAL :: seg_mdtfixflag(nbox * nmodes)
! mdtfixsink stores the amount of *total mass" (over all components) that is
! removed when the fix is applied, to support mass budget calculation. Units
! are moles/s.
REAL :: seg_mdtfixsink(nbox * nmodes)

! Molecular concentration of water (molecules per particle)
REAL :: seg_mdwat(nbox * nmodes)
! Mean free path of air (m)
REAL :: seg_mfpa(nbox)
! Aerosol ptcl number density for mode (cm^-3)
REAL :: seg_nd(nbox * nmodes)
! Air pressure at lower interface (Pa)
REAL :: seg_plower(nbox)
! Air pressure at mid-point (Pa)
REAL :: seg_pmid(nbox)
! Air pressure at upper interface (Pa)
REAL :: seg_pupper(nbox)
! Aerosol partial volume of each cpt in each mode
REAL :: seg_pvol(nbox * nmodes * glomap_variables%ncp)
! Aerosol partial volume of water in each mode
REAL :: seg_pvol_wat(nbox * nmodes)
! Relative humidity (fraction) - gridbox mean
REAL :: seg_rh(nbox)
! Relative humidity (fraction) - clear sky portion
REAL :: seg_rh_clr(nbox)
! Air density (kg/m3)
REAL :: seg_rhoa(nbox)
! Total particle density [incl. H2O & insoluble cpts] (kgm^-3)
REAL :: seg_rhopar(nbox * nmodes)
! Specific humidity (kg/kg)
REAL :: seg_s(nbox)
! Partial masses of gas phase species (kg per gridbox)
REAL :: seg_s0(nbox * nadvg)
! ASAD tendencies for condensable gas phase species (vmr per s)
REAL :: seg_s0_dot_condensable(nbox * nchemg)
! Surface area concentration for each mode (cm^2 / cm^3)
REAL :: seg_sarea(nbox * nmodes)
! Fraction of horizontal gridbox area containing seaice
REAL :: seg_seaice(nbox)
! Mass of air in gridbox (kg)
REAL :: seg_sm(nbox)
! Surface type: 0=seasurf, 1=landsurf, 2=above-seasurf, 3=above-landsurf
REAL :: seg_surtp(nbox)
! Air temperature at mid-point (K)
REAL :: seg_t(nbox)
! Local variable to hold re-shaped aerosol tracers
REAL :: seg_tr_rs(nbox)
! Square-root of centre level temperature (K)
REAL :: seg_tsqrt(nbox)
! Surface friction velocity (m/s)
REAL :: seg_ustr(nbox)
! Local temporary staging for multiplications
REAL :: seg_v1d_tmp(nbox)
! Mean free speed of air molecules (m/s)
REAL :: seg_vba(nbox)
! Volume concentration for each mode
REAL :: seg_vconc(nbox * nmodes)
! Vertical low cloud fraction
REAL :: seg_vfac(nbox)
! Geometric mean wet diameter of particles in each mode (m)
REAL :: seg_wetdp(nbox * nmodes)
! Geometric mean wet volume of particles in each mode (m^3)
REAL :: seg_wvol(nbox * nmodes)
! Background conc. of H2O2 (molecules per cc)
REAL :: seg_zh2o2(nbox)
! Background conc. of HO2 (molecules per cc)
REAL :: seg_zho2(nbox)
! Roughness length (m)
REAL :: seg_znotg(nbox)
! Background vmr of O3 (dimensionless)
REAL :: seg_zo3(nbox)

INTEGER :: errcode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_AERO_CTL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

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

IF (glomap_config%l_ukca_radaer) THEN
  CALL ukca_check_radaer_coupling(all_ntp)
END IF

! Molar mass of dry air (kg/mol)
mm_da = avogadro * boltzmann/rgas

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
  WRITE(umMessage, '(A45,2I6)') 'Set up aerosol etc., NTRAER,NBUDAER=',        &
                                ntraer, nbudaer
  CALL umPrint(umMessage, src='ukca_aero_ctl')

  DO i = 1, nmodes
    WRITE(umMessage, '(A45,I5,L7,3E12.3)')                                     &
              'I,MODE(I),DDPLIM0(I),DDPLIM1(I),SIGMAG(I)=',                    &
               i,mode(i),ddplim0(i),ddplim1(i),sigmag(i)
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    DO j = 1, ncp
      WRITE(umMessage, '(A35,2I5,L7,E12.3)')                                   &
                 'I,J,COMPONENT(I,J),MFRAC_0(I,J)=',                           &
                  i,j,component(i,j),mfrac_0(i,j)
      CALL umPrint(umMessage, src='ukca_aero_ctl')
    END DO
  END DO
END IF

! below are the input configuration parameters set for MODE

IF (firstcall .AND. verbose >0 ) THEN
  WRITE(umMessage, '(A25,I6)') 'i_mode_setup =      ',                         &
                               glomap_config%i_mode_setup
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A25,F6.2)') 'mode_parfrac =      ',                       &
                                 ukca_config%mode_parfrac
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A25,I6)') 'i_mode_nucscav =    ',                         &
                               glomap_config%i_mode_nucscav
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A25,I6)') 'i_mode_nzts =       ',                         &
                               glomap_config%i_mode_nzts
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A25,L7)') 'l_mode_bhn_on =     ',                         &
                               glomap_config%l_mode_bhn_on
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A25,L7)') 'l_mode_bln_on =     ',                         &
                               glomap_config%l_mode_bln_on
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A25,I6)') 'i_mode_bln_param_method',                      &
                               glomap_config%i_mode_bln_param_method
  CALL umPrint(umMessage, src='ukca_aero_ctl')
END IF

! Set whether or not dry deposition of aerosols is on
IF (glomap_config%l_ddepaer) THEN
  ddepaer_on = 1
ELSE
  ddepaer_on = 0
END IF
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage, '(A22,L6,I5)') 'L_DDEPAER,DDEPAER_ON=',                     &
                                  glomap_config%l_ddepaer, ddepaer_on
  CALL umPrint(umMessage, src='ukca_aero_ctl')
END IF

! Set whether or not wet deposition (nucleation and impaction scavenging) of
! aerosols is on
IF (glomap_config%l_rainout) THEN
  rainout_on = 1
ELSE
  rainout_on = 0
END IF
IF (glomap_config%l_impc_scav) THEN
  imscav_on = 1
ELSE
  imscav_on = 0
END IF
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage, '(A22,L6,I5)') 'L_RAINOUT,RAINOUT_ON=',                     &
                                 glomap_config%l_rainout, rainout_on
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A24,L6,I5)') 'L_IMPC_SCAV,IMSCAV_ON=',                    &
                                 glomap_config%l_impc_scav, imscav_on
  CALL umPrint(umMessage, src='ukca_aero_ctl')
END IF

! set INUCSCAV according to i_mode_nucscav setting
inucscav = glomap_config%i_mode_nucscav
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage, '(A25,2I5)') 'I_MODE_NUCSCAV,INUCSCAV=',                    &
                                glomap_config%i_mode_nucscav,inucscav
  CALL umPrint(umMessage, src='ukca_aero_ctl')
END IF

! set LCVRAINOUT according to l_cv_rainout setting
lcvrainout = glomap_config%l_cv_rainout
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage, '(A25,2L6)') 'L_CV_RAINOUT,LCVRAINOUT=',                    &
                                glomap_config%l_cv_rainout, lcvrainout
  CALL umPrint(umMessage, src='ukca_aero_ctl')
END IF

! If new Slinn impaction scavenging scheme is on for dust then turn on flag
! (input to UKCA_AERO_STEP) and initialise scavenging arrays
l_dust_slinn_impc_scav = glomap_config%l_dust_slinn_impc_scav
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage, '(A24,L7)') 'L_DUST_SLINN_IMPC_SCAV=',                      &
                               glomap_config%l_dust_slinn_impc_scav
  CALL umPrint(umMessage, src='ukca_aero_ctl')
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
  WRITE(umMessage, '(A15,2I5)') 'I_MODE_NZTS,NZTS=',                           &
                               glomap_config%i_mode_nzts,nzts
  CALL umPrint(umMessage, src='ukca_aero_ctl')
END IF
!
i_mode_act_method = 1
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage, '(A44)') 'setting I_MODE_ACT_METHOD to default value'
  CALL umPrint(umMessage, src='ukca_aero_ctl')
END IF
!
! set IACTMETHOD according to i_mode_act_method setting
iactmethod = i_mode_act_method
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage, '(A30,2I6)') 'I_MODE_ACT_METHOD,IACTMETHOD=',               &
                                i_mode_act_method, iactmethod
  CALL umPrint(umMessage, src='ukca_aero_ctl')
END IF
!
! set ACT according to mode_activation_dryr setting
act = rmdi
IF (ABS(glomap_config%mode_activation_dryr - rmdi) < EPSILON(0.0)) THEN
  IF (ichem == 1) THEN
    cmessage = ' mode_activation_dryr has not been set'
    errcode = 1
    WRITE(umMessage, '(A40)') cmessage
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
  END IF
ELSE
  act = glomap_config%mode_activation_dryr*1.0e-9 ! convert nm to m
END IF

! Setting removed fraction of oxidised SO2 to a default of 0
IF (ABS(glomap_config%mode_incld_so2_rfrac) < EPSILON(0.0)) THEN
  WRITE(umMessage, '(A57)') 'MODE_INCLD_SO2_RFRAC has been set' //             &
                           ' to default value of 0.0'
  CALL umPrint(umMessage, src='ukca_aero_ctl')
END IF

! Calculate the production rate scaling factor
scale_delso2 = 1.0 - glomap_config%mode_incld_so2_rfrac

IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage, '(A26,2E12.3)') 'MODE_ACTIVATION_DRYR,ACT=',                &
                                   glomap_config%mode_activation_dryr, act
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A22,E12.3)') 'MODE_INCLD_SO2_RFRAC=',                     &
                                  glomap_config%mode_incld_so2_rfrac
  CALL umPrint(umMessage, src='ukca_aero_ctl')
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

! Turn off wet oxidation, condensation and coagulation
! if chemistry is off reflecting that the only aerosol
! scheme that uses this configuration is 2-mode passive
! dust
IF (ichem == 1) THEN
  wetox_on = 1
  cond_on = 1
  coag_on = 1
ELSE
  wetox_on = 0
  cond_on = 0
  coag_on = 0
END IF
!
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage, '(A16,2I6)') 'NUCL_ON,BLN_ON=', nucl_on,bln_on
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A10,I6)') 'IBLN=', ibln
  CALL umPrint(umMessage, src='ukca_aero_ctl')
END IF

! Set whether nitrate emissions are handled here
hno3_uptake_coeff = glomap_config%hno3_uptake_coeff
IF (glomap_config%l_no3_prod_in_aero_step .AND.                                &
    glomap_config%l_ukca_fine_no3_prod) THEN
  fine_no3_prod_on = 1
ELSE
  fine_no3_prod_on = 0
END IF
IF (glomap_config%l_no3_prod_in_aero_step .AND.                                &
    glomap_config%l_ukca_coarse_no3_prod) THEN
  coarse_no3_prod_on = 1
ELSE
  coarse_no3_prod_on = 0
END IF
IF (firstcall .AND. verbose > 0) THEN
  WRITE(umMessage, '(A22,I6)') 'FINE_NO3_PROD_ON   = ', fine_no3_prod_on
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A22,I6)') 'COARSE_NO3_PROD_ON = ', coarse_no3_prod_on
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A22,E15.4)') 'HNO3_UPTAKE_COEFF = ', hno3_uptake_coeff
  CALL umPrint(umMessage, src='ukca_aero_ctl')
END IF

!
dtm = dtc/REAL(nmts)
dtz = dtm/REAL(nzts)
field_size = row_length*rows
field_size3d = field_size*model_levels

IF (firstcall) THEN
  ALLOCATE(mode_tracer_debug(n_mode_tracers))
  mode_tracer_debug(:) = .TRUE.     ! all tracers with debug o/p
END IF

IF (verbose > 1) THEN

  WRITE(umMessage, '(A,3I6)') 'nbox, field_size, field_size3d=',               &
                               nbox, field_size, field_size3d
  CALL umPrint(umMessage, src='ukca_aero_ctl')

  WRITE(umMessage, '(A)') 'UKCA_MODE INPUT SETTINGS : '
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A14,I6)') 'i_mode_setup=', glomap_config%i_mode_setup
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A14,F8.2)') 'mode_parfrac=', ukca_config%mode_parfrac
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A16,I6)') 'i_mode_nucscav=', glomap_config%i_mode_nucscav
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A16,I6)') 'i_mode_nzts=', glomap_config%i_mode_nzts
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A16,L1)') 'l_mode_bhn_on=', glomap_config%l_mode_bhn_on
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A16,L1)') 'l_mode_bln_on=', glomap_config%l_mode_bln_on
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A26,I6)') 'i_mode_bln_param_method',                      &
                               glomap_config%i_mode_bln_param_method
  CALL umPrint(umMessage, src='ukca_aero_ctl')

  WRITE(umMessage, '(A16,I6)') 'i_month: ', i_month
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A16,I6)') 'i_day_number: ', i_day_number
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A16,I6)') 'i_hour: ', i_hour
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A16,I6)') 'i_minute: ', i_minute
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A16,F8.2)') 'DTC: ', dtc
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A16,I6)') 'model_levels: ', model_levels
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A16,I6)') 'rows: ', rows
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A16,I6)') 'row_length: ', row_length
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A22,I6)') 'n_chemistry_tracers: ', n_chemistry_tracers
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A18,I6)') 'n_mode_tracers: ', n_mode_tracers
  CALL umPrint(umMessage, src='ukca_aero_ctl')

  WRITE(umMessage, '(A40)') 'Array:     MIN        MAX         MEAN'
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  l = 0
  WRITE(umMessage, '(A9,I6)') 'Level: ', l
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  WRITE(umMessage, '(A8,3E12.3)') 'p_bdrs: ', MINVAL(p_bdrs(:, :, l)),         &
                               MAXVAL(p_bdrs(:, :, l)),                        &
                               SUM(p_bdrs(:, :, l))/REAL(SIZE(p_bdrs(:, :, l)))
  CALL umPrint(umMessage, src='ukca_aero_ctl')

  ! No model level 2 in UKCA box model - diagnostic print statements
  DO l = 1, MIN(2, model_levels)            ! model_levels
    WRITE(umMessage, '(A9,I6)') 'Level: ',l
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A8,3E12.3)') 'pres: ',MINVAL(pres(:, :, l)),            &
               MAXVAL(pres(:, :, l)),                                          &
               SUM(pres(:, :, l))/REAL(SIZE(pres(:, :, l)))
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A8,3E12.3)') 'temp: ',MINVAL(temp(:, :, l)),            &
               MAXVAL(temp(:, :, l)),                                          &
               SUM(temp(:, :, l))/REAL(SIZE(temp(:, :, l)))
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A8,3E12.3)') 'q: ', MINVAL(q(:, :, l)),                 &
               MAXVAL(q(:, :, l)),                                             &
               SUM(q(:, :, l))/REAL(SIZE(q(:, :, l)))
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A8,3E12.3)') 'rh3d: ',MINVAL(rh3d(:, :, l)),            &
               MAXVAL(rh3d(:, :, l)),                                          &
               SUM(rh3d(:, :, l))/REAL(SIZE(rh3d(:, :, l)))
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A8,3E12.3)') 'p_bdrs: ',MINVAL(p_bdrs(:, :, l)),        &
               MAXVAL(p_bdrs(:, :, l)),                                        &
               SUM(p_bdrs(:, :, l))/REAL(SIZE(p_bdrs(:, :, l)))
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A18,3E12.3)') 'delso2_wet_h2o2: ',                      &
               MINVAL(delso2_wet_h2o2(:, :, l)),                               &
               SUM(delso2_wet_h2o2(:, :, l))/                                  &
               REAL(SIZE(delso2_wet_h2o2(:, :, l)))
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A18,3E12.3)') 'delso2_wet_o3  : ',                      &
               MINVAL(delso2_wet_o3(:, :, l)),                                 &
               SUM(delso2_wet_o3(:, :, l))/REAL(SIZE(delso2_wet_o3(:, :, l)))
  END DO
  IF (model_levels > 7) THEN
    l = 8
    WRITE(umMessage, '(A18,I4,3E12.3)') 'delso2_wet_h2o2: ',l,                 &
          MINVAL(delso2_wet_h2o2(:, :, l)),                                    &
          MAXVAL(delso2_wet_h2o2(:, :, l)),                                    &
          SUM(delso2_wet_h2o2(:, :, l))/REAL(SIZE(delso2_wet_h2o2(:, :, l)))
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A18,I4,3E12.3)') 'delso2_wet_o3  : ',l,                 &
               MINVAL(delso2_wet_o3(:, :, l)),                                 &
               MAXVAL(delso2_wet_o3(:, :, l)),                                 &
               SUM(delso2_wet_o3(:, :, l))/REAL(SIZE(delso2_wet_o3(:, :, l)))
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A18,I4,3E12.3)') 'delso2_dry_oh  : ',l,                 &
                        MINVAL(delso2_dry_oh(:, :, l)),                        &
                        MAXVAL(delso2_dry_oh(:, :, l)),                        &
                        SUM(delso2_dry_oh(:, :, l))/SIZE(delso2_dry_oh(:, :, l))
    CALL umPrint(umMessage, src='ukca_aero_ctl')
  END IF

  IF (ichem == 1) THEN
    DO l = 1, MIN(2, model_levels)            ! model_levels
      DO j = 1, n_chemistry_tracers
        WRITE(umMessage, '(A8,I4,A10,I4,A10)') 'Level: ', l, ' Tracer: ', j,   &
                                                             all_tracer_names(j)
        CALL umPrint(umMessage, src='ukca_aero_ctl')
        WRITE(umMessage, '(A20,3E12.3)') 'chemistry_tracers: ',                &
                                    MINVAL(all_tracers(:, :, l, j)),           &
                                    MAXVAL(all_tracers(:, :, l, j)),           &
                                    SUM(all_tracers(:, :, l, j))/              &
                                    REAL(SIZE(all_tracers(:, :, l, j)))
        CALL umPrint(umMessage, src='ukca_aero_ctl')
      END DO
    END DO
  END IF
  DO l = 1, MIN(2, model_levels)            ! model_levels
    DO j = 1, n_mode_tracers
      IF (mode_tracer_debug(j)) THEN
        iaer = n_chemistry_tracers + j
        WRITE(umMessage, '(A8,I4,A10,I4,A10)') 'Level: ', l,' Tracer: ', j,    &
                                                          all_tracer_names(iaer)
        CALL umPrint(umMessage, src='ukca_aero_ctl')
        WRITE(umMessage, '(A18,3E12.3)') 'mode_tracers: ',                     &
                MINVAL(all_tracers(:, :, l, iaer)),                            &
                MAXVAL(all_tracers(:, :, l, iaer)),                            &
                SUM(all_tracers(:, :, l, iaer)) /                              &
                REAL(SIZE(all_tracers(:, :, l, iaer)))
        CALL umPrint(umMessage, src='ukca_aero_ctl')
      END IF
    END DO
  END DO     ! model_levels

END IF ! IF (verbose > 1)

! Calculate number of aerosol tracers required for components and number
n_reqd_tracers = 0
DO imode = 1, nmodes
  DO icp = 1, ncp
    IF (component(imode,icp)) n_reqd_tracers = n_reqd_tracers + 1
  END DO
END DO
n_reqd_tracers = n_reqd_tracers + SUM(mode_choice)

IF (firstcall) THEN

  ! .. Check the number of tracers, warn if too many, stop if too few
  IF (n_mode_tracers > n_reqd_tracers) THEN
    errcode = -1
    cmessage = ' Too many tracers input'
    WRITE(umMessage, '(A50,2I5)') cmessage,n_mode_tracers,n_reqd_tracers
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
  END IF
  IF (n_mode_tracers < n_reqd_tracers) THEN
    errcode = 1
    cmessage = ' Too few advected aerosol tracers input'
    WRITE(umMessage, '(A50,2I5)') cmessage, n_mode_tracers, n_reqd_tracers
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
  END IF

  ! Check that all tracer addresses are in range  (nmr_index and mmr_index are
  ! for all UKCA tracers, so subtract (ifirst-1) to index the mode tracer array)

  ifirst = jpctr + 1
  IF (verbose > 0) THEN
    WRITE(umMessage, '(A43)') 'Checking MODE tracer addresses are in range:'
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A40)') 'Description, imode, [icp], ifirst, itra'
    CALL umPrint(umMessage, src='ukca_aero_ctl')
  END IF
  DO imode = 1, nmodes
    IF (mode(imode)) THEN
      itra = nmr_index(imode) - ifirst + 1
      icp = 0
      IF (verbose > 0) THEN
        WRITE(umMessage, '(A10,4I6)') 'Number:  ', imode, icp, ifirst, itra
        CALL umPrint(umMessage, src='ukca_aero_ctl')
      END IF
      IF (itra <= 0 .OR. itra > n_mode_tracers) THEN
        errcode = 1
        cmessage = 'Tracer address out of range for number'
        WRITE(umMessage, '(A72,2(A8,I6))') cmessage, ' mode: ', imode,         &
                                           ' index: ', nmr_index(imode)
        CALL umPrint(umMessage, src='ukca_aero_ctl')
        CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
      END IF
      DO icp = 1, ncp
        IF (component(imode,icp)) THEN
          itra = mmr_index(imode,icp) - ifirst + 1
          IF (verbose > 0) THEN
            WRITE(umMessage, '(A10,4I6)') 'Mass MR: ', imode, icp, ifirst,     &
                                          itra
            CALL umPrint(umMessage, src='ukca_aero_ctl')
          END IF
          IF (itra < 0 .OR. itra > n_mode_tracers) THEN
            errcode = 1
            cmessage = 'Tracer address out of range for component'
            WRITE(umMessage, '(A72,3(A12,I6))') cmessage, ' mode: ', imode,    &
                    ' component: ', icp, ' index: ', nmr_index(imode)
            CALL umPrint(umMessage, src='ukca_aero_ctl')
            CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
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
all_ntp(i)%data_3d(:, :, :) = 0.0
n_merge_3d(:, :, :, :) = 0
! temporary clumping of autoconv variables for later segment extraction
a3d_tmp(:, :, :) = autoconv(:, :, :) + accretion(:, :, :) +                    &
                   rim_agg(:, :, :) + rim_cry(:, :, :)

! this is used to track the number of cells where MDT was modified due to
! unphysical tracer values
nbadmdt_3d(:, :, :, :) = 0

! other precomputations
CALL log_v(nmodes, sigmag, log_sigmag)

! Allocate the drydiam array (input to UKCA_ACTIVATE) if required.
IF ( ( glomap_config%i_ukca_activation_scheme == i_ukca_activation_arg )       &
     .AND. .NOT. ALLOCATED(drydiam) ) THEN
  ALLOCATE(drydiam(row_length, rows, model_levels, nmodes))
  drydiam(:, :, :, :) = 0.0
END IF

!
! OpenMP parallel region starts here (loop later around IK segments)
!
!$OMP PARALLEL  DEFAULT(NONE)                                                  &
!$OMP          SHARED(a3d_tmp, act, all_ntp, bln_on,                           &
!$OMP all_tracers, cloud_frac, cloud_liq_frac, cloud_liq_wat,                  &
!$OMP n_chemistry_tracers, component, crain, csnow,                            &
!$OMP ddepaer_on, delso2_dry_oh, delso2_wet_h2o2, delso2_wet_o3,               &
!$OMP drain, drydiam, dryox_in_aer, dsnow, dtc, dtm, dtz,                      &
!$OMP firstcall, glomap_config, glomap_variables,                              &
!$OMP iactmethod, ibln, iextra_checks, imscav_on, inucscav,                    &
!$OMP jpctr, l_dust_slinn_impc_scav, l_ukca_cmip6_diags, l_ukca_mode_diags,    &
!$OMP l_ukca_pm_diags, land_fraction, lbase, lcvrainout,                       &
!$OMP log_sigmag, lscat_zhang,                                                 &
!$OMP mass, mdtmin, mdwat_diag, mfrac_0, mh2o2f, mh2so4, mlo,                  &
!$OMP mm, mm_da, mm_gas, mmid, mmr_index,                                      &
!$OMP mode, mode_diags, model_levels, modesol, mox,                            &
!$OMP msec_org, msec_orgi, msotwo,                                             &
!$OMP n_h2o2, n_h2so4, n_merge_3d, n_o3, n_sec_org, n_sec_org_i, n_so2,        &
!$OMP n_hono2, n_nh3, mhno3, mnh3,                                             &
!$OMP nadvg, nbadmdt_3d, nbox, nbox_s, nbudaer, nchemg, ichem, ncol_s, ncp,    &
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
!$OMP nmasagedduintr63, nmasagedduintr74, nmasagedduintr84,                    &
!$OMP nmasmergsuintr12, nmasmergocintr12, nmasmergsointr12, nmasmergsuintr23,  &
!$OMP nmasmergbcintr23, nmasmergocintr23, nmasmergsointr23, nmasmergsuintr34,  &
!$OMP nmasmergbcintr34, nmasmergocintr34, nmasmergssintr34, nmasmergduintr34,  &
!$OMP nmasddepdusupins, nmasnuscdusupins, nmasimscdusupins, nmascondsusupins,  &
!$OMP nmascondocsupins, nmascondsosupins, nmascoagsuintr18, nmascoagocintr18,  &
!$OMP nmascoagsointr18, nmasmergsointr34,                                      &
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
!$OMP nmasprimntaitsol, nmasprimntaccsol,nmasprimntcorsol, nmasprimnhaitsol,   &
!$OMP nmasprimnhaccsol, nmasprimnhcorsol,nmascondnnaccsol, nmascondnncorsol,   &
!$OMP nmax_mode_diags,nmr_index,nucl_on, nseg,nukca_d1items,                   &
!$OMP wetox_on, cond_on, coag_on,                                              &
!$OMP fine_no3_prod_on, coarse_no3_prod_on, hno3_uptake_coeff,                 &
!$OMP num_eps, nzts, p_bdrs, pres, q,                                          &
!$OMP rainout_on, rgas, rh3d, rh3d_clr, root2, row_length, rows,               &
!$OMP sea_ice_frac, scale_delso2, sigmag, stride_s, temp, u_s,                 &
!$OMP ukca_config, ukcaD1codes, verbose, verbose_local,                        &
!$OMP wetdp_diag, wetox_in_aer, z0m, z_half_alllevs, zbl)                      &
!$OMP     PRIVATE(cmessage, dp0, errcode,                                      &
!$OMP             ifirst, i, i_end, i_end_cp, i_start, i_start_cp,             &
!$OMP             ic, icp, ik, itra, iaer, imode, j, jl, jl2, jv, k,           &
!$OMP             l, lb, logic, logic1, logic2, n, ncs, nbs, nbs_index,        &
!$OMP             seg_aird, seg_airdm3, seg_autoconv1d, seg_bud_aer_mas,       &
!$OMP             seg_ccn_1, seg_ccn_2, seg_ccn_3, seg_ccn_4, seg_ccn_5,       &
!$OMP             seg_cdn,seg_clf, seg_clwc, seg_cn_3nm, seg_craing,           &
!$OMP             seg_craing_up,seg_csnowg, seg_delso2, seg_delso2_2,          &
!$OMP             seg_draing,seg_drydp, seg_dsnowg, seg_dvisc, seg_dvol,       &
!$OMP             seg_erf_arg,seg_erfterm, seg_fac, seg_fconv_conv,            &
!$OMP             seg_height,seg_het_rates, seg_htpblg, seg_iarr, seg_ilscat,  &
!$OMP             seg_jlabove, seg_karr, seg_land_frac, seg_larr, seg_lday,    &
!$OMP             seg_lowcloud, seg_lwc, seg_md, seg_mdt, seg_mdtfixflag,      &
!$OMP             seg_mdtfixsink, seg_mdwat, seg_mfpa, seg_n_merge_1d,         &
!$OMP             seg_nbadmdt, seg_nd, seg_plower, seg_pmid, seg_pupper,       &
!$OMP             seg_pvol, seg_pvol_wat, seg_rh, seg_rh_clr, seg_rhoa,        &
!$OMP             seg_rhopar, seg_s, seg_s0, seg_s0_dot_condensable,           &
!$OMP             seg_sarea, seg_seaice, seg_sm, seg_surtp, seg_t,             &
!$OMP             seg_tr_rs, seg_tsqrt, seg_ustr, seg_v1d_tmp, seg_vba,        &
!$OMP             seg_vconc, seg_vfac, seg_wetdp, seg_wvol, seg_zh2o2,         &
!$OMP             seg_zho2, seg_znotg, seg_zo3,                                &
!$OMP             y)

! dealing with non-uniform segments
IF ( nseg > 1 ) THEN
  IF ( ncol_s(nseg) /= ncol_s(nseg-1) ) THEN
    WRITE(ummessage,'(A)') 'AERO_CTL: Segments are not uniform'
    CALL umprint(ummessage, src='ukca_aero_ctl', pe=0)
    WRITE(ummessage,'(i6,1x,i6,1x,i6)') ncol_s(nseg), ncol_s(nseg-1), nseg
    CALL umprint(ummessage, src='ukca_aero_ctl', pe=0)
  END IF  ! last segment is smaller
END IF  ! nseg > 1

!$OMP DO SCHEDULE(DYNAMIC)
DO ik = 1, nseg
  ! use local alias because they are used in  many places within  loop
  lb = lbase(ik)      ! base location on this segment
  ncs = ncol_s(ik)    ! The number of columns on this segment
  nbs = nbox_s(ik)    ! The number of boxes on this segment
  IF (firstcall .AND. verbose > 2 ) THEN
    WRITE(umMessage, '(A28,5(1x,i8))') 'AERO_CTL:lb,ncs,nbs,nseg,ik',          &
                                      lb, ncs, nbs, nseg, ik
    CALL umPrint(umMessage, src='ukca_aero_ctl')
  END IF

  ! Now work out which grid box is above each grid box
  jl = 0
  DO ic = 1, ncs                        ! loop over the columns in this segment
    DO l = 1, (model_levels-1)               ! let top level have JLABOVE=-1
      jl = jl + 1                            ! the ID of the current box
      seg_jlabove(jl) = jl + 1               ! the ID of the box above this one
      seg_karr(jl)=MOD(lb,row_length) + ic-1 ! this is longitude
      seg_iarr(jl) = ((lb+ic-1)/rows) + 1    ! this is latitude
      seg_larr(jl) = l                       ! this is altitude (level)
    END DO                                   ! l, model_levels
    jl = jl +1                               ! the box at top of column
    seg_jlabove(jl) = -1                     ! there is no box above this one
    seg_karr(jl)=MOD(lb,row_length) + ic -1  ! this is longitude
    seg_iarr(jl) = lb/rows                   ! this is latitude
    seg_larr(jl) = l                         ! this is altitude (level)
  END DO                                     ! for each column

  DO i = 0, nmodes
    nbs_index(i) = i * nbs
  END DO

  DO jl = 1, nbs * nmodes
    seg_n_merge_1d(jl) = 0
  END DO

  ! Reshape input quantities
  ! ========================
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   temp(1, 1, 1) ,seg_t(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   pres(1, 1, 1), seg_pmid(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   p_bdrs(1, 1, 1), seg_pupper(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   p_bdrs(1, 1, 0), seg_plower(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   crain(1, 1, 1), seg_craing(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   drain(1, 1, 1), seg_draing(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   csnow(1, 1, 1), seg_csnowg(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   dsnow(1, 1, 1), seg_dsnowg(1))

  ! .. add in-cloud accretion, ice and snow melt to the autoconversion rate
  ! NOTE previously accumulated a3d_tmp = autoconv+accretion+rim_agg+rim_cry
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   a3d_tmp(1, 1, 1), seg_autoconv1d(1))

  ! .. set CRAING_UP using JLABOVE as calculated above
  DO jl = 1, nbs   ! limit is number of boxes on this segment
    IF (seg_jlabove(jl) > 0) THEN
      seg_craing_up(jl) = seg_craing( seg_jlabove(jl) )
    ELSE
      seg_craing_up(jl) = 0.0
    END IF
  END DO
  !
  DO jl = 1, nbs
    ! .. currently set FCONV_CONV=0.99 -- need to change to take as input
    ! fraction of condensate-->rain in 6 hrs
    seg_fconv_conv(jl) = 0.99
    ! .. weakened rainout -- FCONV_CONV=0.5
    !!      FCONV_CONV(:)=0.50 ! fraction of condensate-->rain in 6 hrs
    !
    ! Calculate molecular concentration of air. No conc of air (/cm3)
    seg_aird(jl) = seg_pmid(jl) / (seg_t(jl) * boltzmann * 1.0e6)
  END DO

  ! copy from delso2_wet_xxx arrays as output from UKCA_CHEMISTRY_CTL
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   delso2_wet_h2o2(1, 1, 1), seg_delso2(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   delso2_wet_o3(1, 1, 1), seg_delso2_2(1))

  ! Scale in-cloud production rates to account for lack of
  ! in-cloud wet removal of oxidised SO2
  DO jl = 1, nbs
    seg_delso2  (jl) = seg_delso2(jl)   * scale_delso2
    seg_delso2_2(jl) = seg_delso2_2(jl) * scale_delso2
  END DO

  !
  ! Set these to zero as will not be used in UKCA_AERO_STEP
  ! Currently do wet ox separately in UM
  DO jl = 1, nbs
    seg_zo3(jl) = 0.0
    seg_zho2(jl) = 0.0
    seg_zh2o2(jl) = 0.0
    seg_lday(jl) = 0
  END DO
  !
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   rh3d(1, 1, 1), seg_rh(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   rh3d_clr(1, 1, 1), seg_rh_clr(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   q(1, 1, 1), seg_s(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   cloud_liq_wat(1, 1, 1), seg_lwc(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   cloud_liq_wat(1, 1, 1), seg_clwc(1))

  DO jl = 1, nbs
    seg_lowcloud(jl) = 0.0
    seg_vfac(jl)     = 0.0
    seg_clf(jl)      = 0.0
  END DO

  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   cloud_liq_frac(1, 1, 1), seg_clf(1))
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   cloud_frac(1, 1, 1), seg_lowcloud(1))

  DO jl = 1, nbs
    IF (seg_lowcloud(jl) > 0.0) THEN
      ! set to 1 so that VFAC*LOWCLOUD=cloud_frac
      seg_vfac(jl) = 1.0
    END IF
  END DO

  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   z_half_alllevs(1, 1, 1), seg_height(1))

  ! GM added code here to only set cloud fraction > 0 if in low cloud
  !    here low cloud is defined as being cloud with p>=680hPa
  ! PMID is in Pa
  DO jl = 1, nbs
    IF (seg_pmid(jl) < 680.0e2) THEN
      seg_lowcloud(jl) = 0.0
      seg_lwc(jl) = 0.0
    END IF
  END DO
  !
  ! 2D surface -> whole segment
  CALL surface_to_seg(lb, ncs, nbs, stride_s, model_levels,                    &
                      u_s(1, 1), seg_ustr(1))
  CALL surface_to_seg(lb, ncs, nbs, stride_s, model_levels,                    &
                      z0m(1, 1), seg_znotg(1))
  CALL surface_to_seg(lb, ncs, nbs, stride_s, model_levels,                    &
                      sea_ice_frac(1, 1), seg_seaice(1))
  ! fraction of land at surface
  CALL surface_to_seg(lb, ncs, nbs, stride_s, model_levels,                    &
                      land_fraction(1, 1),seg_land_frac(1))
  CALL surface_to_seg(lb, ncs, nbs, stride_s, model_levels,                    &
                      zbl(1, 1), seg_htpblg(1))

  jl = 0
  DO ic = 1, ncs                 ! for each column on segment
    jl = jl +1                   ! at surface layer jl = 1
    IF ( seg_land_frac(jl) < 0.5 ) THEN
      seg_surtp(jl) = 0.0
    ELSE
      seg_surtp(jl) = 1.0
    END IF
    DO l = 2, model_levels       ! above the surface layer
      jl = jl +1
      IF ( seg_land_frac(jl) < 0.5 ) THEN
        seg_surtp(jl) = 2.0
      ELSE
        seg_surtp(jl) = 3.0
      END IF
    END DO
  END DO

  DO jl = 1, nbs
    seg_ilscat(jl) = 0
  END DO

  ! Corrected surface types array, using dominant land_surface type rather
  ! than roughness length
  IF ( glomap_config%l_improve_aero_drydep ) THEN
    CALL surface_to_seg_int(lb, ncs, nbs, stride_s, model_levels,              &
                            lscat_zhang(1, 1), seg_ilscat(1))
  ELSE
    ! Old method
    !---------------------------------------------------------------
    ! Put in section here to set land-surface category ILSCAT based on ZNOTG
    ! ILSCAT=1-9 based on 9 UM landsurf types but use existing approach
    ! to set 4 possible types according to roughness length
    DO jl = 1, nbs
      IF (seg_znotg(jl) < 1.0e-3) THEN
        ! Water/sea - z0<0.001m
        seg_ilscat(jl) = 7
      ELSE IF (seg_znotg(jl) > 1.0e-1) THEN
        ! Forest
        seg_ilscat(jl) = 1
      ELSE
        ! All other lands, grass 0.001<z0<0.1m
        seg_ilscat(jl) = 3
      END IF
    END DO

    ! If sea ice covers > 50% of sea surface, treat as sea ice
    DO jl = 1, nbs
      IF (seg_seaice(jl) > 0.5) THEN
        ! Sea ice
        seg_ilscat(jl) = 9
      END IF
    END DO

  END IF
  !---------------------------------------------------------------

  ! Derived quantities
  ! ==================
  DO jl = 1, nbs
    ! no conc of air (/m3)
    seg_airdm3(jl) = seg_aird(jl) * 1.0e6
    seg_rhoa(jl) = seg_pmid(jl) / (seg_t(jl) * rgas)
    seg_vba(jl) = SQRT( 8.0 * boltzmann * seg_t(jl) / (pi * ma) )
    seg_tsqrt(jl) = SQRT(seg_t(jl))
    seg_dvisc(jl) = 1.83e-5 * (416.16 / (seg_t(jl) + 120.0)) *                 &
                    (SQRT(seg_t(jl) / 296.16)**3)
    seg_mfpa(jl) = 2.0 * seg_dvisc(jl)/(seg_rhoa(jl) * seg_vba(jl))
  END DO

  ! Mass air (kg/box)
  CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                       &
                   mass(1, 1, 1), seg_sm(1))

  ! Gas-phase tracers required in aerosol code
  ! ==========================================
  !  S0 array that is passed in to aerosol code uses MH2SO4, msec_org, etc
  !     to index tracers (set in UKCA_SETUP_INDICES module procedures)

  ! Set gas phase tracer masses to zero
  DO itra = 1, nadvg
    i_start = nbs * (itra - 1)
    DO jl = 1, nbs
      seg_s0(i_start + jl) = 0.0
    END DO
  END DO

  IF (dryox_in_aer == 0) THEN
    ! Condensable tracer tendencies->0
    DO i = 1, nchemg
      i_start = nbs * (i - 1)
      DO jl = 1, nbs
        seg_s0_dot_condensable(i_start + jl) = 0.0
      END DO
    END DO
    ! DRYOX_IN_AER=0 -> update of condensables done in UKCA_CHEMISTRY_CTL
  END IF ! if DRYOX_IN_AER=0

  IF (dryox_in_aer == 1) THEN
    ! Initialise S0_DOT_CONDENSABLE to 0
    DO i = 1, nchemg
      i_start = nbs * (i - 1)
      DO jl = 1, nbs
        seg_s0_dot_condensable(i_start + jl) = 0.0
      END DO
    END DO
    ! DRYOX_IN_AER=1 -> update of condensables done in UKCA_AERO_STEP
    ! condensable tracer tendencies need to be set here to pass in as input

    CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                     &
                     delso2_dry_oh(1, 1, 1), seg_v1d_tmp(1))

    IF (ukca_config%l_fix_ukca_h2so4_ystore) THEN
      i_start = nbs * (mh2so4 - 1)
      DO jl = 1, nbs
        seg_s0_dot_condensable(i_start + jl) = seg_v1d_tmp(jl)
      END DO
      ! .. delso2_dry_oh  is in units of vmr/s, correct for
      ! .. S0_DOT_CONDENSABLE
    ELSE
      i_start = nbs * (mh2so4 - 1)
      DO jl = 1, nbs
        seg_s0_dot_condensable(i_start + jl) = seg_v1d_tmp(jl) /               &
                                               (seg_aird(jl) * dtc)
      END DO
      ! .. delso2_dry_oh  is in units of molecules/cc/DTC
      ! .. need S0_DOT_CONDENSABLE to be in units of vmr/s
      ! .. so need to divide by (AIRD(:)*DTC)
    END IF

    IF (msec_org > 0) THEN
      ! For L_classSO2_inAer=T or F set Sec_Org prodn--> 0 (done in UKCA)
      i_start = nbs * (msec_org - 1)
      DO jl = 1, nbs
        seg_s0_dot_condensable(i_start + jl) = 0.0
      END DO
    END IF

    IF (msec_orgi > 0) THEN
      ! For L_classSO2_inAer=T or F set Sec_Org prodn--> 0 (done in UKCA)
      i_start = nbs * (msec_orgi - 1)
      DO jl = 1, nbs
        seg_s0_dot_condensable(i_start + jl) = 0.0
      END DO
    END IF
    !
  END IF      ! dryox_in_aer
  ! Set H2O2, O3 and SO2 for aqueous oxidation
  IF (wetox_in_aer == 1) THEN
    ! .. set h2o2 when required
    IF (mh2o2f  > 0 .AND. mm_gas(mh2o2f) > 1e-3 ) THEN
      CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                   &
                       all_tracers(1, 1, 1, n_h2o2), seg_v1d_tmp(1))
      i_start = nbs * (mh2o2f - 1)
      DO jl = 1, nbs
        seg_s0(i_start + jl) = seg_sm(jl) * (mm_da / mm_gas(mh2o2f)) *         &
                               seg_v1d_tmp(jl)
      END DO
    ELSE
      cmessage = ' H2O2 needs updating, but MH2O2F' //                         &
                 'or MM_GAS(MH2O2F) is wrong'
      WRITE(umMessage, '(A60,I6,E12.3)') cmessage, mh2o2f, mm_gas(mh2o2f)
      CALL umPrint(umMessage, src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
    END IF

    ! .. set O3
    IF (mox > 0 .AND. mm_gas(mox) > 1e-3 ) THEN
      CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                   &
                       all_tracers(1, 1, 1, n_o3), seg_v1d_tmp(1))
      DO jl = 1, nbs
        seg_zo3(jl) = seg_sm(jl) * (mm_da / mm_gas(mox)) * seg_v1d_tmp(jl)
      END DO
    ELSE
      cmessage = ' O3 needs updating, but MOX or MM_GAS(MOX) is wrong'
      WRITE(umMessage, '(A52,I6,E12.3)') cmessage, ' MOX = ', mox,mm_gas(mox)
      CALL umPrint(umMessage, src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
    END IF

    ! .. SO2 mmr in kg[SO2]/kg[dryair]
    IF (msotwo > 0 .AND. mm_gas(msotwo) > 1e-3) THEN
      CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                   &
                       all_tracers(1, 1, 1, n_so2), seg_v1d_tmp(1))
      i_start = nbs * (msotwo - 1)
      DO jl = 1, nbs
        seg_s0(i_start + jl) = seg_sm(jl) * (mm_da / mm_gas(msotwo)) *         &
                               seg_v1d_tmp(jl)
      END DO
    ELSE
      cmessage = ' SO2 needs updating, but MSOTWO or MM_GAS(MSOTWO)'//         &
                 'is wrong'
      WRITE(umMessage, '(A20,I6,E12.3)') cmessage, ' MSOTWO = ', msotwo,       &
            mm_gas(msotwo)
      CALL umPrint(umMessage, src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
    END IF

  END IF      ! wetox_in_aer=1

  IF (mh2so4   > 0) THEN
    CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                     &
                     all_tracers(1, 1, 1, n_h2so4), seg_v1d_tmp(1))
    i_start = nbs * (mh2so4 - 1)
    DO jl = 1, nbs
      seg_s0(i_start + jl) = seg_sm(jl) * (mm_da / mm_gas(mh2so4)) *           &
                             seg_v1d_tmp(jl)
    END DO
  END IF

  ! .. Secondary Organic tracer mmr in kg[Sec_Org]/kg[dryair]
  IF (msec_org > 0) THEN
    CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                     &
                     all_tracers(1, 1, 1, n_sec_org), seg_v1d_tmp(1))
    i_start = nbs * (msec_org - 1)
    DO jl = 1, nbs
      seg_s0(i_start + jl) = seg_sm(jl) * (mm_da / mm_gas(msec_org)) *         &
                             seg_v1d_tmp(jl)
    END DO
  END IF

  ! .. Secondary Organic tracer from isoprene mmr in kg[SEC_ORG I]/kg[dryair]
  IF (msec_orgi > 0) THEN
    CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                     &
                     all_tracers(1, 1, 1, n_sec_org_i), seg_v1d_tmp(1))
    i_start = nbs * (msec_orgi - 1)
    DO jl = 1, nbs
      seg_s0(i_start + jl) = seg_sm(jl) * (mm_da / mm_gas(msec_orgi)) *        &
                             seg_v1d_tmp(jl)
    END DO
  END IF

  ! .. HNO3 tracer - initialise even if fine NO3 off via fine_no3_prod_on and
  !                  coarse NO3 off via coarse_no3_prod_on
  IF (mhno3 > 0) THEN
    CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                     &
                     all_tracers(1, 1, 1, n_hono2), seg_v1d_tmp(1))
    i_start = nbs * (mhno3 - 1)
    DO jl = 1, nbs
      seg_s0(i_start + jl) = seg_sm(jl) * (mm_da / mm_gas(mhno3)) *            &
                             seg_v1d_tmp(jl)
    END DO
  END IF

  IF (mnh3 > 0) THEN
    CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                     &
                     all_tracers(1, 1, 1, n_nh3), seg_v1d_tmp(1))
    i_start = nbs * (mnh3 - 1)
    DO jl = 1, nbs
      seg_s0(i_start + jl) = seg_sm(jl) * (mm_da / mm_gas(mnh3)) *             &
                             seg_v1d_tmp(jl)
    END DO
  END IF


  ! Aerosol Tracers
  ! ===============
  !  Find index of 1st mode tracer, as nmr_index and
  !   mmr_index index all ukca tracers
  ifirst = jpctr + 1
  DO imode = 1, nmodes
    DO jl = 1, nbs
      jl2 = nbs_index(imode - 1) + jl
      seg_nbadmdt(jl2) = 0
      seg_mdtfixflag(jl2) = 0.0
      seg_mdtfixsink(jl2) = 0.0
    END DO
    IF (mode(imode)) THEN
      iaer = n_chemistry_tracers + nmr_index(imode) - ifirst + 1
      CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                   &
                       all_tracers(1, 1, 1, iaer), seg_tr_rs(1))
      DO jl = 1, nbs
        ! Set tr_rs to zero if negative
        seg_tr_rs(jl) = MAX(0.0, seg_tr_rs(jl))
        ! Set ND (particles per cc) from advected number-mixing-ratio
        seg_nd(nbs_index(imode - 1) + jl) = seg_tr_rs(jl) * seg_aird(jl)
      END DO

      DO icp = 1, ncp
        IF (component(imode, icp)) THEN
          i_start_cp = (icp - 1) * nmodes * nbs
          iaer = n_chemistry_tracers + mmr_index(imode, icp) - ifirst + 1
          CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,               &
                           all_tracers(1, 1, 1, iaer), seg_tr_rs(1))

          DO jl = 1, nbs
            jl2 = nbs_index(imode - 1) + jl
            ! Set tr_rs to zero if negative
            seg_tr_rs(jl) = MAX(0.0, seg_tr_rs(jl))

            ! Set MD (molecules / particle) from advected mass-mixing-ratio
            ! note that only "trusts" advected values where ND>NUM_EPS
            IF (seg_nd(jl2) > num_eps(imode)) THEN
              seg_md(i_start_cp + jl2) = (mm_da / mm(icp)) *                   &
                                        seg_aird(jl) * seg_tr_rs(jl) /         &
                                        seg_nd(jl2)
            ELSE
              seg_md(i_start_cp + jl2) = mmid(imode) * mfrac_0(imode, icp)
            END IF
          END DO
        ELSE
          i_start_cp = nbs * (imode - 1 + (icp - 1) * nmodes)
          DO jl = 1, nbs
            seg_md(i_start_cp + jl) = 0.0
          END DO
        END IF
      END DO ! loop over cpts
      !
      ! Set total mass array MDT from sum over individual component MDs
      DO jl = 1, nbs
        seg_mdt(nbs_index(imode - 1) + jl) = 0.0
      END DO
      DO icp = 1, ncp
        IF (component(imode, icp)) THEN
          i_start_cp = (icp - 1) * nmodes * nbs
          DO jl = 1, nbs
            jl2 = nbs_index(imode - 1) + jl
            seg_mdt(jl2) = seg_mdt(jl2) + seg_md(i_start_cp + jl2)
          END DO
        END IF
      END DO
      !
      ! below checks if MDT is coming out too low after advection (af BLMIX)
      !
      IF (iextra_checks == 0) THEN
        mdtmin(imode, ik) = mlo(imode) * 0.001 ! set equiv. to DPLIM0*0.1
        !
        DO icp = 1, ncp
          IF (component(imode, icp)) THEN
            ! Where MDT too low after advection set ND to zero and &
            ! set default MD
            i_start_cp = (icp - 1) * nmodes * nbs
            DO jl = 1, nbs
              jl2 = nbs_index(imode - 1) + jl
              IF (seg_mdt(jl2) < mdtmin(imode, ik)) THEN
                seg_md(i_start_cp + jl2) = mmid(imode) * mfrac_0(imode, icp)
              END IF
            END DO
          END IF
        END DO

        ! Count occurrences (NBADMDT) and store percent occurrence and mass-sink
        DO jl = 1, nbs
          seg_fac(jl) = seg_sm(jl) / seg_aird(jl)
        END DO
        ! FAC converts aerosol mass fluxes from kg(dryair)/box/tstep to
        ! moles/gridbox/s

        DO jl = 1, nbs
          jl2 = nbs_index(imode - 1) + jl
          IF (seg_mdt(jl2) < mdtmin(imode, ik)) THEN
            seg_nbadmdt(jl2) = 1
            ! mdtfixflag enables to track proportion of timesteps that fix is
            ! applied. units of mdtfixsink are moles/s
            seg_mdtfixflag(jl2) = 100.0
            ! mdtfixsink stores total mass removed when fix is applied
            seg_mdtfixsink(jl2) = seg_nd(jl2) * seg_mdt(jl2) *                 &
                                  seg_fac(jl) / mm_da / dtc
          ELSE
            seg_nbadmdt(jl2) = 0
            seg_mdtfixflag(jl2) = 0.0
            seg_mdtfixsink(jl2) = 0.0
          END IF
        END DO

        ! Put count of bad mdt cells into a 3D structure for later diagnostic
        CALL int_ins_seg(lb, ncs, nbs, stride_s, model_levels,                 &
                         seg_nbadmdt(nbs_index(imode-1)+1:nbs_index(imode)),   &
                         nbadmdt_3d(1, 1, 1, imode))

        ! Set ND->0 where MDT too low (& set MDT->MMID) & count occurrences
        DO jl = 1, nbs
          jl2 = nbs_index(imode - 1) + jl
          IF (seg_mdt(jl2) < mdtmin(imode, ik)) THEN
            seg_nd(jl2) = 0.0
            ! Where MDT too low after advection (but +ve) set to MMID
            seg_mdt(jl2) = mmid(imode)
          END IF
        END DO
      END IF ! Iextra_checks = 0
    ELSE
      DO jl = 1, nbs
        jl2 = nbs_index(imode - 1) + jl
        seg_nd(jl2) = 0.0
        seg_mdt(jl2) = mmid(imode)
      END DO
      DO icp = 1, ncp
        i_start_cp = nbs * (imode - 1 + (icp - 1) * nmodes)
        IF (component(imode, icp)) THEN
          DO jl = 1, nbs
            seg_md(i_start_cp + jl) = mmid(imode) * mfrac_0(imode, icp)
          END DO
        ELSE
          DO jl = 1, nbs
            seg_md(i_start_cp + jl) = 0.0
          END DO
        END IF
      END DO
    END IF
  END DO ! loop over modes
  !
  IF (verbose >= 2) THEN

    IF (ichem == 1) THEN
      DO itra = 1, nadvg
        CALL select_array_segment(nbs, 1, itra, i_start, i_end)
        WRITE(umMessage, '(A10,I4,2E12.3)') 'S0 : ',itra,                      &
                                            MINVAL(seg_s0(i_start:i_end)),     &
                                            MAXVAL(seg_s0(i_start:i_end))
        CALL umPrint(umMessage, src='ukca_aero_ctl')
      END DO
    END IF

    DO imode = 1, nmodes
      IF (mode(imode)) THEN
        CALL select_array_segment(nbs, 1, imode, i_start, i_end)
        WRITE(umMessage, '(A10,I4,2E12.3)') 'ND : ',imode,                     &
                                            MINVAL(seg_nd(i_start:i_end)),     &
                                            MAXVAL(seg_nd(i_start:i_end))
        CALL umPrint(umMessage, src='ukca_aero_ctl')
        WRITE(umMessage, '(A10,I4,2E12.3)') 'MDT: ',imode,                     &
                                            MINVAL(seg_mdt(i_start:i_end)),    &
                                            MAXVAL(seg_mdt(i_start:i_end))
        CALL umPrint(umMessage, src='ukca_aero_ctl')
        DO icp = 1, ncp
          IF (component(imode, icp)) THEN
            CALL select_array_segment(nbs, 1, imode, i_start_cp, i_end_cp,     &
                                      dim2_len = nmodes, dim3_index = icp)
            WRITE(umMessage, '(A10,I4,2E12.3)') 'MD : ',imode,                 &
                                      MINVAL(seg_md(i_start_cp:i_end_cp)),     &
                                      MAXVAL(seg_md(i_start_cp:i_end_cp))
            CALL umPrint(umMessage, src='ukca_aero_ctl')
          END IF
        END DO
      END IF
    END DO
  END IF  ! verbose > 2
  !
  ! .. Apply extra checks for advection artefacts
  IF (iextra_checks > 0) THEN
    CALL ukca_mode_check_artefacts(verbose_local, nbs, seg_nd,  seg_md,        &
                                   seg_mdt, seg_sm, seg_aird, mm_da, dtc,      &
                                   seg_mdtfixflag, seg_mdtfixsink)
  END IF

  ! .. zero aerosol budget terms before calling UKCA_AERO_STEP
  DO jv = 0, nbudaer
    i_start = jv * nbs
    DO jl = 1, nbs
      seg_bud_aer_mas(i_start + jl) = 0.0
    END DO
  END DO

  ! .. below is call to UKCA_AERO_STEP as at ukca_mode_v1_gm1.f90

  IF (firstcall .AND. verbose > 0) THEN
    WRITE(umMessage, '(A42)') 'Values of input variables passed to mode:'
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'RAINOUT_ON=', rainout_on
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'IMSCAV_ON=', imscav_on
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'WETOX_ON=', wetox_on
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'DDEPAER_ON=', ddepaer_on
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'SEDI_ON=', sedi_on
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'DRYOX_IN AER=', dryox_in_aer
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'WETOX_IN AER=', wetox_in_aer
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'COND_ON=', cond_on
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'NUCL_ON=', nucl_on
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'BLN_ON=', bln_on
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'FINE_NO3_ON=', fine_no3_prod_on
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'COARSE_NO3_ON=', coarse_no3_prod_on
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A19,E15.4)') 'HNO3_UPTAKE_COEFF=', hno3_uptake_coeff
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'COAG_ON=', coag_on
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'ICOAG=', icoag
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'IMERGE=', imerge
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'IFUCHS=', ifuchs
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'IDCMFP=', idcmfp
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'ICONDIAM=', icondiam
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'IBLN=', ibln
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'IACTMETHOD=', iactmethod
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'I_NUC_METHOD=', i_nuc_method
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'IDDEPAER=', iddepaer
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'INUCSCAV=', inucscav
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A12,L7)') 'lcvrainout=', lcvrainout
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A23,L7)') 'l_dust_slinn_impc_scav=',                    &
                                 l_dust_slinn_impc_scav
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'VERBOSE=', verbose
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'VERBOSE_LOCAL=', verbose_local
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'CHECKMD_ND=', checkmd_nd
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'INTRAOFF=', intraoff
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'INTEROFF=', interoff
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    WRITE(umMessage, '(A15,I4)') 'IDUSTEMS=', idustems
    CALL umPrint(umMessage, src='ukca_aero_ctl')
  END IF
  !
  CALL ukca_aero_step(nbs, nchemg, nadvg, nbudaer,                             &
                      seg_nd, seg_mdt, seg_md, seg_mdwat, seg_s0, seg_drydp,   &
                      seg_wetdp, seg_rhopar, seg_dvol, seg_wvol, seg_sm,       &
                      seg_aird, seg_airdm3, seg_rhoa, seg_mfpa, seg_dvisc,     &
                      seg_t, seg_tsqrt, seg_rh, seg_rh_clr, seg_s, seg_pmid,   &
                      seg_pupper, seg_plower, seg_zo3, seg_zho2, seg_zh2o2,    &
                      seg_ustr, seg_znotg, seg_surtp, seg_craing, seg_draing,  &
                      seg_craing_up, seg_csnowg, seg_dsnowg, seg_fconv_conv,   &
                      seg_lowcloud, seg_vfac, seg_clf, seg_autoconv1d,         &
                      dtc, dtz, nmts, nzts, seg_lday, act, seg_bud_aer_mas,    &
                      rainout_on, iextra_checks, imscav_on, wetox_on,          &
                      ddepaer_on, sedi_on, iso2wetoxbyo3, dryox_in_aer,        &
                      wetox_in_aer, seg_delso2, seg_delso2_2,                  &
                      cond_on, nucl_on, coag_on, bln_on, icoag, imerge,        &
                      fine_no3_prod_on, coarse_no3_prod_on, hno3_uptake_coeff, &
                      ifuchs, idcmfp, icondiam, ibln, i_nuc_method,            &
                      iactmethod, iddepaer, inucscav, ichem,                   &
                      lcvrainout, l_dust_slinn_impc_scav, verbose_local,       &
                      checkmd_nd, intraoff, interoff,                          &
                      seg_s0_dot_condensable, seg_lwc,seg_clwc, seg_pvol,      &
                      seg_pvol_wat, seg_jlabove, seg_ilscat, seg_n_merge_1d,   &
                      seg_height, seg_htpblg)

  !
  ! Update tracers
  ! ==============

  IF (wetox_in_aer > 0) THEN
    IF (mh2o2f > 0) THEN
      ! .. update gas phase H2O2 mmr following SO2 aqueous phase oxidation
      i_start = nbs * (mh2o2f - 1)
      DO jl = 1, nbs
        seg_v1d_tmp(jl) = (mm_gas(mh2o2f) / mm_da) * seg_s0(i_start + jl) /    &
                          seg_sm(jl)
      END DO
      CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                    &
                      seg_v1d_tmp(1:nbs), all_tracers(1, 1, 1, n_h2o2))
    END IF
    IF (msotwo > 0) THEN
      ! .. update gas phase SO2 mmr following aqueous phase oxidation
      i_start = nbs * (msotwo - 1)
      DO jl = 1, nbs
        seg_v1d_tmp(jl) = (seg_s0(i_start + jl) / seg_sm(jl)) *                &
                          (mm_gas(msotwo) / mm_da)
      END DO
      CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                    &
                      seg_v1d_tmp(1:nbs), all_tracers(1, 1, 1, n_so2))
    END IF
  END IF

  IF (mh2so4   > 0) THEN
    ! .. update gas phase H2SO4 mmr following H2SO4 condensation/nucleation
    i_start = nbs * (mh2so4 - 1)
    DO jl = 1, nbs
      seg_v1d_tmp(jl) = (seg_s0(i_start + jl) / seg_sm(jl)) *                  &
                        (mm_gas(mh2so4) / mm_da)
    END DO
    CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                      &
                    seg_v1d_tmp(1:nbs), all_tracers(1, 1, 1, n_h2so4))
  END IF

  IF (msec_org > 0) THEN
    ! .. update gas phase Sec_Org mmr following condensation/nucleation
    i_start = nbs * (msec_org - 1)
    DO jl = 1, nbs
      seg_v1d_tmp(jl) = (seg_s0(i_start + jl) / seg_sm(jl)) *                  &
                        (mm_gas(msec_org)/mm_da)
    END DO
    CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                      &
                    seg_v1d_tmp(1:nbs), all_tracers(1, 1, 1, n_sec_org))
  END IF

  IF (msec_orgi > 0) THEN
    ! .. update gas phase Sec_Org mmr following condensation/nucleation
    i_start = nbs * (msec_orgi - 1)
    DO jl = 1, nbs
      seg_v1d_tmp(jl) = (seg_s0(i_start + jl) / seg_sm(jl)) *                  &
                        (mm_gas(msec_orgi)/mm_da)
    END DO
    CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                      &
                    seg_v1d_tmp(1:nbs),                                        &
                    all_tracers(1, 1, 1, n_sec_org_i))
  END IF

  IF (mhno3   > 0) THEN
    ! .. update gas phase HNO3 mmr following nitrate production
    i_start = nbs * (mhno3 - 1)
    DO jl = 1, nbs
      seg_v1d_tmp(jl) = (seg_s0(i_start + jl) / seg_sm(jl)) *                  &
                        (mm_gas(mhno3) / mm_da)
    END DO
    CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                      &
                    seg_v1d_tmp(1:nbs), all_tracers(1, 1, 1, n_hono2))
  END IF

  IF (mnh3   > 0) THEN
    ! .. update gas phase NH3 mmr following ammonium nitrate production
    i_start = nbs * (mnh3 - 1)
    DO jl = 1, nbs
      seg_v1d_tmp(jl) = (seg_s0(i_start + jl) / seg_sm(jl)) *                  &
                        (mm_gas(mnh3) / mm_da)
    END DO
    CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                      &
                    seg_v1d_tmp(1:nbs), all_tracers(1, 1, 1, n_nh3))
  END IF

  ! Set H2O2, O3 and SO2 for aqueous oxidation
  IF (wetox_in_aer == 1) THEN
    ! .. set h2o2 from h2o2_tracer when required
    IF (mh2o2f  > 0 .AND. mm_gas(mh2o2f) > 1e-3 ) THEN
      CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                   &
                       all_tracers(1, 1, 1, n_h2o2), seg_v1d_tmp(1))
      i_start = nbs * (mh2o2f - 1)
      DO jl = 1, nbs
        seg_s0(i_start + jl) = seg_sm(jl) * (mm_da / mm_gas(mh2o2f)) *         &
                               seg_v1d_tmp(jl)
      END DO
    ELSE
      cmessage = ' H2O2 needs updating, but MH2O2F'//                          &
                 'or MM_GAS(MH2O2F) is wrong'
      WRITE(umMessage, '(A10,I4,E12.3)') cmessage, mh2o2f, mm_gas(mh2o2f)
      CALL umPrint(umMessage, src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
    END IF

    ! .. set O3 from O3_tracer
    IF (mox > 0 .AND. mm_gas(mox) > 1e-3 ) THEN
      CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                   &
                       all_tracers(1, 1, 1, n_o3), seg_v1d_tmp(1))
      DO jl = 1, nbs
        seg_zo3(jl) = seg_sm(jl) * (mm_da / mm_gas(mox)) * seg_v1d_tmp(jl)
      END DO
    ELSE
      cmessage = ' O3 needs updating, but MOX or MM_GAS(MOX) is wrong'
      WRITE(umMessage, '(A10,I4,E12.3)') cmessage, ' MOX = ', mox, mm_gas(mox)
      CALL umPrint(umMessage, src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
    END IF

    IF (msotwo > 0 .AND. mm_gas(msotwo) > 1e-3) THEN
      CALL extract_seg(lb, ncs, nbs, stride_s, model_levels,                   &
                       all_tracers(1, 1, 1, n_so2), seg_v1d_tmp(1))
      i_start = nbs * (msotwo - 1)
      DO jl = 1, nbs
        seg_s0(i_start + jl) = seg_sm(jl) * (mm_da / mm_gas(msotwo)) *         &
                               seg_v1d_tmp(jl)
      END DO
    ELSE
      cmessage = ' SO2 needs updating, but MSOTWO or MM_GAS(MSOTWO)'//         &
                 'is wrong'
      WRITE(umMessage, '(A12,I4,E12.3)') cmessage, ' MSOTWO = ', msotwo,       &
                                         mm_gas(msotwo)
      CALL umPrint(umMessage, src='ukca_aero_ctl')
      errcode = 1
      CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
    END IF
  END IF      ! wetox_in_aer

  DO imode = 1, nmodes
    DO jl = 1, nbs
      seg_sarea(nbs_index(imode-1) + jl) = 0.0
    END DO
  END DO

  DO imode = 1, nmodes
    IF (mode(imode)) THEN
      iaer = n_chemistry_tracers + nmr_index(imode) - ifirst + 1
      ! .. update aerosol no. conc. following aerosol microphysics
      DO jl = 1, nbs
        seg_v1d_tmp(jl) = (seg_nd(nbs_index(imode - 1) + jl) / seg_aird(jl))
      END DO
      CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                    &
                      seg_v1d_tmp(1:nbs), all_tracers(1, 1, 1, iaer) )
      DO icp = 1, ncp
        IF (component(imode,icp)) THEN
          i_start_cp = (icp - 1) * nmodes * nbs
          DO jl = 1, nbs
            jl2 = nbs_index(imode - 1) + jl
            IF (seg_nd(jl2) <= num_eps(imode)) THEN
              seg_md(i_start_cp + jl2) = mmid(imode) * mfrac_0(imode, icp)
            END IF
            ! .. update aerosol mmr following aerosol microphysics
            seg_v1d_tmp(jl) = (mm(icp) / mm_da) *                              &
                              ( seg_md(i_start_cp + jl2) * seg_nd(jl2) /       &
                                seg_aird(jl) )
          END DO

          iaer = n_chemistry_tracers + mmr_index(imode,icp) - ifirst + 1
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_v1d_tmp(1:nbs), all_tracers(1, 1, 1, iaer) )
        END IF
      END DO ! loop over cpts

      CALL int_ins_seg(lb, ncs, nbs, stride_s, model_levels,                   &
                       seg_n_merge_1d(nbs_index(imode-1)+1:nbs_index(imode)),  &
                       n_merge_3d(1, 1, 1, imode) )
      y(imode) = EXP(2.0 * (log_sigmag(imode)**2))
      DO jl = 1, nbs
        jl2 = nbs_index(imode-1) + jl
        seg_vconc(jl2) = seg_wvol(jl2) * seg_nd(jl2)

        ! calculate surface area in units of units of cm^2 / cm^3.
        ! nd is aerosol ptcl number density for mode (cm^-3)
        ! wetdp is wet diameter in m
        ! So multiply by 1.0e+4 to convert from m^2 / cm^3 to cm^2 / cm^3
        seg_sarea(jl2) = pi * 1.0e+4 * seg_nd(jl2) * (seg_wetdp(jl2)**2) *     &
                         y(imode)
      END DO
    END IF     ! mode(imode)
  END DO       ! imode = 1, nmodes

  !--------------------------------------------------------
  ! below sets CN,CCN,CDN diagnostics

  ! Initialise to zero
  DO jl = 1, nbs
    seg_cn_3nm(jl) = 0.0
    seg_ccn_2(jl)  = 0.0
    seg_ccn_3(jl)  = 0.0
    seg_ccn_4(jl)  = 0.0
    seg_ccn_5(jl)  = 0.0
  END DO

  DO imode = 1, nmodes
    IF (mode(imode)) THEN
      DO jl = 1, nbs
        jl2 = nbs_index(imode - 1) + jl
        !
        dp0 = 3.0e-9
        seg_erf_arg(jl) = LOG(dp0 / seg_drydp(jl2)) /                          &
                          (root2 * log_sigmag(imode))
        seg_erfterm(jl) = 0.5 * seg_nd(jl2) * (1.0 - umErf(seg_erf_arg(jl)))
        seg_cn_3nm(jl) = seg_cn_3nm(jl) + seg_erfterm(jl)
        !
      END DO
      IF (modesol(imode) == 1) THEN
        DO jl = 1, nbs
          jl2 = nbs_index(imode - 1) + jl
          !
          ! For CCN_2 take CCN for particles > 25nm dry radius
          dp0 = 50.0e-9
          seg_erf_arg(jl) = LOG(dp0 / seg_drydp(jl2)) /                        &
                            (root2 * log_sigmag(imode))
          seg_erfterm(jl) = 0.5 * seg_nd(jl2) * (1.0 - umErf(seg_erf_arg(jl)))
          seg_ccn_2(jl) = seg_ccn_2(jl) + seg_erfterm(jl)
          !
          ! For CCN_3 take CCN for particles > 35nm dry radius
          dp0 = 70.0e-9
          seg_erf_arg(jl) = LOG(dp0 / seg_drydp(jl2)) /                        &
                            (root2 * log_sigmag(imode))
          seg_erfterm(jl) = 0.5 * seg_nd(jl2) * (1.0 - umErf(seg_erf_arg(jl)))
          seg_ccn_3(jl) = seg_ccn_3(jl) + seg_erfterm(jl)
          !
        END DO
      END IF
      DO jl = 1, nbs
        jl2 = nbs_index(imode - 1) + jl
        ! For CCN_4 take CCN for particles > 15nm dry radius, including insoluble
        dp0 = 30.0e-9
        seg_erf_arg(jl) = LOG(dp0 / seg_drydp(jl2)) /                          &
                          (root2 * log_sigmag(imode))
        seg_erfterm(jl) = 0.5 * seg_nd(jl2) * (1.0 - umErf(seg_erf_arg(jl)))
        seg_ccn_4(jl) = seg_ccn_4(jl) + seg_erfterm(jl)
        !
        ! For CCN_5 take CCN for particles > 25nm dry radius, including insoluble
        dp0 = 50.0e-9
        seg_erf_arg(jl) = LOG(dp0 / seg_drydp(jl2)) /                          &
                          (root2 * log_sigmag(imode))
        seg_erfterm(jl) = 0.5 * seg_nd(jl2) * (1.0 - umErf(seg_erf_arg(jl)))
        seg_ccn_5(jl) = seg_ccn_5(jl) + seg_erfterm(jl)
        !
      END DO
    END IF
  END DO

  !----------------------------------------------------
  CALL ukca_cdnc_jones(nbs, act, seg_drydp, seg_nd,                            &
                       glomap_variables, seg_ccn_1, seg_cdn)
  !----------------------------------------------------

  !
  ! Update DRYDP and DVOL after AERO_STEP and updation of MD
  CALL ukca_calc_drydiam(nbs, glomap_variables,                                &
                         seg_nd, seg_md, seg_mdt, seg_drydp, seg_dvol)

  IF (verbose > 1) THEN

    WRITE(umMessage, '(A30)') 'AFTER CALL TO UKCA_AERO_STEP'
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    IF (ichem == 1) THEN
      DO itra = 1, nadvg
        CALL select_array_segment(nbs, 1, itra, i_start, i_end)
        WRITE(umMessage, '(A10,I4,2E12.3)') 'S0 : ', itra,                     &
                                            MINVAL(seg_s0(i_start:i_end)),     &
                                            MAXVAL(seg_s0(i_start:i_end))
        CALL umPrint(umMessage, src='ukca_aero_ctl')
      END DO
    END IF
    DO imode = 1, nmodes
      IF (mode(imode)) THEN
        CALL select_array_segment(nbs, 1, imode, i_start, i_end)
        WRITE(umMessage, '(A10,I4,2E12.3)') 'ND : ', imode,                    &
                                            MINVAL(seg_nd(i_start:i_end)),     &
                                            MAXVAL(seg_nd(i_start:i_end))
        CALL umPrint(umMessage, src='ukca_aero_ctl')
        WRITE(umMessage, '(A10,I4,2E12.3)') 'MDT: ', imode,                    &
                                            MINVAL(seg_mdt(i_start:i_end)),    &
                                            MAXVAL(seg_mdt(i_start:i_end))
        CALL umPrint(umMessage, src='ukca_aero_ctl')
        DO icp = 1, ncp
          IF (component(imode,icp)) THEN
            CALL select_array_segment(nbs, 1, imode, i_start_cp, i_end_cp,     &
                                      dim2_len = nmodes, dim3_index = icp)
            WRITE(umMessage, '(A10,I4,2E12.3)') 'MD : ', imode,                &
                                      MINVAL(seg_md(i_start_cp:i_end_cp)),     &
                                      MAXVAL(seg_md(i_start_cp:i_end_cp))
            CALL umPrint(umMessage, src='ukca_aero_ctl')
          END IF
        END DO
      END IF
    END DO

    DO imode = 1, nmodes
      IF (mode(imode)) THEN
        WRITE(umMessage, '(A10,I4,2E12.3)') 'DRYDP: ', imode,                  &
                  MINVAL(seg_drydp(nbs_index(imode-1)+1:nbs_index(imode))),    &
                  MAXVAL(seg_drydp(nbs_index(imode-1)+1:nbs_index(imode)))
        CALL umPrint(umMessage, src='ukca_aero_ctl')
      END IF
    END DO

  END IF ! if VERBOSE > 1

  ! Calculate heterogeneous rate coeffs for tropospheric chemistry
  IF (ukca_config%l_ukca_trophet) THEN
    CALL ukca_trop_hetchem(nbs,nhet, seg_t,seg_rh, seg_aird,                   &
                           seg_pvol, seg_wetdp, seg_sarea, seg_het_rates)
    ! Now copy the het_rates into the all_ntp array
    CALL select_array_segment(nbs, 1,  ihet_n2o5, i_start, i_end)
    i = name2ntpindex('het_n2o5  ')
    CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                      &
                    seg_het_rates(i_start:i_end), all_ntp(i)%data_3d(1, 1, 1))
    CALL select_array_segment(nbs, 1,  ihet_ho2_ho2, i_start, i_end)
    i = name2ntpindex('het_ho2   ')
    CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                      &
                    seg_het_rates(i_start:i_end), all_ntp(i)%data_3d(1, 1, 1))
  END IF

  ! Copy MODE diagnostics from BUD_AER_MAS to BUD_AER_MAS1D
  ! and convert all budget variables to be in kg/DTC

  DO jl = 1, nbs
    seg_fac(jl) = seg_sm(jl) / seg_aird(jl)
  END DO
  ! FAC converts aerosol mass flux from kg(dryair)/box/tstep to moles/gridbox/s
  DO jv = 1, nbudaer
    ! NMASCLPR variables are in kg(dryair)/DTC, others in molecules/cc/DTC
    logic1=(jv <  nmasclprsuaitsol1)
    logic2=(jv >  nmasclprsucorsol2)
    logic=logic1 .OR. logic2
    IF (logic) THEN
      ! When passed out from AERO_STEP, BUD_AER_MAS is in molecules/cc/DTC
      ! Moles/s
      DO jl = 1, nbs
        seg_bud_aer_mas(jv * nbs + jl) = seg_bud_aer_mas(jv * nbs + jl) *      &
                                         seg_fac(jl) / mm_da / dtc
      END DO
    ELSE
      ! When passed out from AERO_STEP, BUD_AER_MAS (NMASCLPR) is
      ! in kg(dryair)/box/DTC
      ! Moles/s
      DO jl = 1, nbs
        seg_bud_aer_mas(jv * nbs + jl) = seg_bud_aer_mas(jv * nbs + jl) /      &
                                         mm_da / dtc
      END DO
    END IF
  END DO

  ! Write 3_D diagnostics to mode_diags array
  !  N.B. L_ukca_mode_diags is set whenever a STASH request for a relevant item
  !  is found, and the ukcaD1codes(N)%item is then set, otherwise it is IMDI

  IF (L_ukca_mode_diags) THEN  ! fill 3D array
    k = 0
    DO n = 1, nukca_d1items
      IF (ukcaD1codes(n)%section == stashcode_glomap_sec .AND.                 &
          ((ukcaD1codes(n)%item >= item1_mode_diags .AND.                      &
            ukcaD1codes(n)%item <= item1_mode_diags+nmax_mode_diags-1) .OR.    &
           (ukcaD1codes(n)%item >= item1_nitrate_diags .AND.                   &
            ukcaD1codes(n)%item <= itemN_nitrate_diags) .OR.                   &
           (ukcaD1codes(n)%item >= item1_dust3mode_diags .AND.                 &
            ukcaD1codes(n)%item <= itemN_dust3mode_diags)) .AND.               &
          ukcaD1codes(n)%item /= imdi) THEN
        k = k + 1

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
        ! cnccncdn - 437-441, 700-701
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
        ! ****************** DUST 3RD MODE :
        ! prim DU  -- 675
        ! ddep DU  -- 676
        ! nusc DU  -- 677
        ! imsc DU  -- 678
        ! cond SU,OC,SO -- 679-681
        ! htox SU  -- 682
        ! coag SU,OC,SO -- 683-685
        ! aged DU  -- 686-688
        ! drydp8   -- 689
        ! sarea8   -- 690
        ! vconc8   -- 691
        ! rhop8    -- 692
        ! pvol8 DU -- 693
        ! mdt fix mass loss -- 694
        ! mdt fix frequency -- 695

        IF (firstcall .AND. verbose > 0 ) THEN
          WRITE(umMessage, '(A,3I6)') 'About to set mode diagnostics for'//    &
               ' N, k, item=', n, k, ukcaD1codes(n)%item
          CALL umPrint(umMessage, src='ukca_aero_ctl')
        END IF
        !
        SELECT CASE(UkcaD1codes(n)%item)
        CASE (item1_mode_diags:item1_mode_diags+12)
          !           Do nothing, primary emissions not handled here now
        CASE (item1_mode_diags+13)
          CALL select_array_segment(nbs, 0, nmasddepsunucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+14)
          CALL select_array_segment(nbs, 0, nmasddepsuaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+15)
          CALL select_array_segment(nbs, 0, nmasddepsuaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+16)
          CALL select_array_segment(nbs, 0, nmasddepsucorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+17)
          CALL select_array_segment(nbs, 0, nmasddepssaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+18)
          CALL select_array_segment(nbs, 0, nmasddepsscorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+19)
          CALL select_array_segment(nbs, 0, nmasddepbcaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+20)
          CALL select_array_segment(nbs, 0, nmasddepbcaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+21)
          CALL select_array_segment(nbs, 0, nmasddepbccorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+22)
          CALL select_array_segment(nbs, 0, nmasddepbcaitins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+23)
          CALL select_array_segment(nbs, 0, nmasddepocnucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+24)
          CALL select_array_segment(nbs, 0, nmasddepocaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+25)
          CALL select_array_segment(nbs, 0, nmasddepocaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+26)
          CALL select_array_segment(nbs, 0, nmasddepoccorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+27)
          CALL select_array_segment(nbs, 0, nmasddepocaitins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+28)
          CALL select_array_segment(nbs, 0, nmasddepsonucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+29)
          CALL select_array_segment(nbs, 0, nmasddepsoaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+30)
          CALL select_array_segment(nbs, 0, nmasddepsoaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+31)
          CALL select_array_segment(nbs, 0, nmasddepsocorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+32)
          CALL select_array_segment(nbs, 0, nmasddepduaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+33)
          CALL select_array_segment(nbs, 0, nmasddepducorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+34)
          CALL select_array_segment(nbs, 0, nmasddepduaccins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+35)
          CALL select_array_segment(nbs, 0, nmasddepducorins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+36)
          CALL select_array_segment(nbs, 0, nmasnuscsunucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+37)
          CALL select_array_segment(nbs, 0, nmasnuscsuaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+38)
          CALL select_array_segment(nbs, 0, nmasnuscsuaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+39)
          CALL select_array_segment(nbs, 0, nmasnuscsucorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+40)
          CALL select_array_segment(nbs, 0, nmasnuscssaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+41)
          CALL select_array_segment(nbs, 0, nmasnuscsscorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+42)
          CALL select_array_segment(nbs, 0, nmasnuscbcaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+43)
          CALL select_array_segment(nbs, 0, nmasnuscbcaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+44)
          CALL select_array_segment(nbs, 0, nmasnuscbccorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+45)
          CALL select_array_segment(nbs, 0, nmasnuscbcaitins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+46)
          CALL select_array_segment(nbs, 0, nmasnuscocnucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+47)
          CALL select_array_segment(nbs, 0, nmasnuscocaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+48)
          CALL select_array_segment(nbs, 0, nmasnuscocaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+49)
          CALL select_array_segment(nbs, 0, nmasnuscoccorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+50)
          CALL select_array_segment(nbs, 0, nmasnuscocaitins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+51)
          CALL select_array_segment(nbs, 0, nmasnuscsonucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+52)
          CALL select_array_segment(nbs, 0, nmasnuscsoaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+53)
          CALL select_array_segment(nbs, 0, nmasnuscsoaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+54)
          CALL select_array_segment(nbs, 0, nmasnuscsocorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+55)
          mode_diags(:, :, :,k)=0.0
          !! .. there is no MASNUSCSOAITINS --- erroneously included in
          !! .. UKCA_mode stash section so set it to zero here
        CASE (item1_mode_diags+56)
          CALL select_array_segment(nbs, 0, nmasnuscduaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+57)
          CALL select_array_segment(nbs, 0, nmasnuscducorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+58)
          CALL select_array_segment(nbs, 0, nmasnuscduaccins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+59)
          CALL select_array_segment(nbs, 0, nmasnuscducorins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+60)
          CALL select_array_segment(nbs, 0, nmasimscsunucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+61)
          CALL select_array_segment(nbs, 0, nmasimscsuaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+62)
          CALL select_array_segment(nbs, 0, nmasimscsuaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+63)
          CALL select_array_segment(nbs, 0, nmasimscsucorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+64)
          CALL select_array_segment(nbs, 0, nmasimscssaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+65)
          CALL select_array_segment(nbs, 0, nmasimscsscorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+66)
          CALL select_array_segment(nbs, 0, nmasimscbcaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+67)
          CALL select_array_segment(nbs, 0, nmasimscbcaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+68)
          CALL select_array_segment(nbs, 0, nmasimscbccorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+69)
          CALL select_array_segment(nbs, 0, nmasimscbcaitins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+70)
          CALL select_array_segment(nbs, 0, nmasimscocnucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+71)
          CALL select_array_segment(nbs, 0, nmasimscocaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+72)
          CALL select_array_segment(nbs, 0, nmasimscocaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+73)
          CALL select_array_segment(nbs, 0, nmasimscoccorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+74)
          CALL select_array_segment(nbs, 0, nmasimscocaitins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+75)
          CALL select_array_segment(nbs, 0, nmasimscsonucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+76)
          CALL select_array_segment(nbs, 0, nmasimscsoaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+77)
          CALL select_array_segment(nbs, 0, nmasimscsoaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+78)
          CALL select_array_segment(nbs, 0, nmasimscsocorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+79)
          CALL select_array_segment(nbs, 0, nmasimscduaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+80)
          CALL select_array_segment(nbs, 0, nmasimscducorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+81)
          CALL select_array_segment(nbs, 0, nmasimscduaccins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+82)
          CALL select_array_segment(nbs, 0, nmasimscducorins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+83)
          CALL select_array_segment(nbs, 0, nmasclprsuaitsol1, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+84)
          CALL select_array_segment(nbs, 0, nmasclprsuaccsol1, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+85)
          CALL select_array_segment(nbs, 0, nmasclprsucorsol1, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+86)
          CALL select_array_segment(nbs, 0, nmasclprsuaitsol2, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+87)
          CALL select_array_segment(nbs, 0, nmasclprsuaccsol2, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+88)
          CALL select_array_segment(nbs, 0, nmasclprsucorsol2, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+89)
          CALL select_array_segment(nbs, 0, nmasprocsuintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+90)
          CALL select_array_segment(nbs, 0, nmasprocbcintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+91)
          CALL select_array_segment(nbs, 0, nmasprococintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+92)
          CALL select_array_segment(nbs, 0, nmasprocsointr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+93)
          CALL select_array_segment(nbs, 0, nmascondsunucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+94)
          CALL select_array_segment(nbs, 0, nmascondsuaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+95)
          CALL select_array_segment(nbs, 0, nmascondsuaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+96)
          CALL select_array_segment(nbs, 0, nmascondsucorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+97)
          CALL select_array_segment(nbs, 0, nmascondsuaitins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+98)
          CALL select_array_segment(nbs, 0, nmascondsuaccins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+99)
          CALL select_array_segment(nbs, 0, nmascondsucorins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+100)
          CALL select_array_segment(nbs, 0, nmascondocnucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+101)
          CALL select_array_segment(nbs, 0, nmascondocaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+102)
          CALL select_array_segment(nbs, 0, nmascondocaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+103)
          CALL select_array_segment(nbs, 0, nmascondoccorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+104)
          CALL select_array_segment(nbs, 0, nmascondocaitins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+105)
          CALL select_array_segment(nbs, 0, nmascondocaccins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+106)
          CALL select_array_segment(nbs, 0, nmascondoccorins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+107)
          CALL select_array_segment(nbs, 0, nmascondsonucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+108)
          CALL select_array_segment(nbs, 0, nmascondsoaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+109)
          CALL select_array_segment(nbs, 0, nmascondsoaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+110)
          CALL select_array_segment(nbs, 0, nmascondsocorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+111)
          CALL select_array_segment(nbs, 0, nmascondsoaitins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+112)
          CALL select_array_segment(nbs, 0, nmascondsoaccins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+113)
          CALL select_array_segment(nbs, 0, nmascondsocorins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
          !!
          !! Code for heterogeneous oxidation of SO2 --> SO4 on dust not yet in
          !!
        CASE (item1_mode_diags+114)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+115)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+116)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+117)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+118)
          CALL select_array_segment(nbs, 0, nmasnuclsunucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+119)
          CALL select_array_segment(nbs, 0, nmascoagsuintr12, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+120)
          CALL select_array_segment(nbs, 0, nmascoagsuintr13, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+121)
          CALL select_array_segment(nbs, 0, nmascoagsuintr14, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+122)
          CALL select_array_segment(nbs, 0, nmascoagsuintr15, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+123)
          CALL select_array_segment(nbs, 0, nmascoagsuintr16, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+124)
          CALL select_array_segment(nbs, 0, nmascoagsuintr17, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+125)
          CALL select_array_segment(nbs, 0, nmascoagocintr12, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+126)
          CALL select_array_segment(nbs, 0, nmascoagocintr13, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+127)
          CALL select_array_segment(nbs, 0, nmascoagocintr14, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+128)
          CALL select_array_segment(nbs, 0, nmascoagocintr15, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+129)
          CALL select_array_segment(nbs, 0, nmascoagocintr16, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+130)
          CALL select_array_segment(nbs, 0, nmascoagocintr17, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+131)
          CALL select_array_segment(nbs, 0, nmascoagsointr12, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+132)
          CALL select_array_segment(nbs, 0, nmascoagsointr13, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+133)
          CALL select_array_segment(nbs, 0, nmascoagsointr14, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+134)
          CALL select_array_segment(nbs, 0, nmascoagsointr15, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+135)
          CALL select_array_segment(nbs, 0, nmascoagsointr16, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+136)
          CALL select_array_segment(nbs, 0, nmascoagsointr17, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+137)
          CALL select_array_segment(nbs, 0, nmascoagsuintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+138)
          CALL select_array_segment(nbs, 0, nmascoagsuintr24, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+139)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+140)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+141)
          CALL select_array_segment(nbs, 0, nmascoagbcintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+142)
          CALL select_array_segment(nbs, 0, nmascoagbcintr24, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+143)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+144)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+145)
          CALL select_array_segment(nbs, 0, nmascoagocintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+146)
          CALL select_array_segment(nbs, 0, nmascoagocintr24, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+147)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+148)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+149)
          CALL select_array_segment(nbs, 0, nmascoagsointr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+150)
          CALL select_array_segment(nbs, 0, nmascoagsointr24, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+151)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+152)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+153)
          CALL select_array_segment(nbs, 0, nmascoagsuintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+154)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+155)
          CALL select_array_segment(nbs, 0, nmascoagbcintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+156)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+157)
          CALL select_array_segment(nbs, 0, nmascoagocintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+158)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+159)
          CALL select_array_segment(nbs, 0, nmascoagssintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+160)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+161)
          CALL select_array_segment(nbs, 0, nmascoagsointr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+162)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+163)
          CALL select_array_segment(nbs, 0, nmascoagduintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+164)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+165)
          CALL select_array_segment(nbs, 0, nmascoagbcintr53, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+166)
          CALL select_array_segment(nbs, 0, nmascoagocintr53, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+167)
          CALL select_array_segment(nbs, 0, nmascoagbcintr54, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+168)
          CALL select_array_segment(nbs, 0, nmascoagocintr54, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+169)
          CALL select_array_segment(nbs, 0, nmascoagduintr64, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+170)
          CALL select_array_segment(nbs, 0, nmasagedsuintr52, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+171)
          CALL select_array_segment(nbs, 0, nmasagedbcintr52, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+172)
          CALL select_array_segment(nbs, 0, nmasagedocintr52, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+173)
          CALL select_array_segment(nbs, 0, nmasagedsointr52, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+174)
          CALL select_array_segment(nbs, 0, nmasmergsuintr12, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+175)
          CALL select_array_segment(nbs, 0, nmasmergocintr12, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+176)
          CALL select_array_segment(nbs, 0, nmasmergsointr12, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+177)
          CALL select_array_segment(nbs, 0, nmasmergsuintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+178)
          CALL select_array_segment(nbs, 0, nmasmergbcintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+179)
          CALL select_array_segment(nbs, 0, nmasmergocintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+180)
          CALL select_array_segment(nbs, 0, nmasmergsointr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+181)
          CALL select_array_segment(nbs, 0, nmasmergsuintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+182)
          CALL select_array_segment(nbs, 0, nmasmergssintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+183)
          CALL select_array_segment(nbs, 0, nmasmergbcintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+184)
          CALL select_array_segment(nbs, 0, nmasmergocintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+185)
          CALL select_array_segment(nbs, 0, nmasmergduintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+186)
          CALL select_array_segment(nbs, 0, nmasmergsointr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+187)
          ! Do nothing here, used for marine OM emissions
        CASE (item1_mode_diags+188)
          ! Do nothing here, used for marine OM emissions
          ! SEC_ORG_I condensation
        CASE (item1_mode_diags+189)
          CALL select_array_segment(nbs, 0, nmascondocinucsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+190)
          CALL select_array_segment(nbs, 0, nmascondociaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+191)
          CALL select_array_segment(nbs, 0, nmascondociaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+192)
          CALL select_array_segment(nbs, 0, nmascondocicorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+193)
          CALL select_array_segment(nbs, 0, nmascondociaitins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+200)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(1:nbs),                                    &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+201)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(nbs_index(1)+1:nbs_index(2)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+202)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(nbs_index(2)+1:nbs_index(3)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+203)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(nbs_index(3)+1:nbs_index(4)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+204)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(nbs_index(4)+1:nbs_index(5)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+205)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(nbs_index(5)+1:nbs_index(6)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+206)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(nbs_index(6)+1:nbs_index(7)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+207)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_wetdp(1:nbs),                                    &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+208)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_wetdp(nbs_index(1)+1:nbs_index(2)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+209)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_wetdp(nbs_index(2)+1:nbs_index(3)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+210)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_wetdp(nbs_index(3)+1:nbs_index(4)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+211)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdwat(nbs_index(0)+1:nbs_index(1))/avogadro,     &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+212)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdwat(nbs_index(1)+1:nbs_index(2))/avogadro,     &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+213)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdwat(nbs_index(2)+1:nbs_index(3))/avogadro,     &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+214)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdwat(nbs_index(3)+1:nbs_index(4))/avogadro,     &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+215)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_sarea(1:nbs),                                    &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+216)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_sarea(nbs_index(1)+1:nbs_index(2)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+217)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_sarea(nbs_index(2)+1:nbs_index(3)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+218)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_sarea(nbs_index(3)+1:nbs_index(4)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+219)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_sarea(nbs_index(4)+1:nbs_index(5)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+220)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_sarea(nbs_index(5)+1:nbs_index(6)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+221)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_sarea(nbs_index(6)+1:nbs_index(7)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+222)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_vconc(1:nbs),                                    &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+223)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_vconc(nbs_index(1)+1:nbs_index(2)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+224)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_vconc(nbs_index(2)+1:nbs_index(3)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+225)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_vconc(nbs_index(3)+1:nbs_index(4)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+226)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_vconc(nbs_index(4)+1:nbs_index(5)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+227)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_vconc(nbs_index(5)+1:nbs_index(6)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+228)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_vconc(nbs_index(6)+1:nbs_index(7)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+229)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(1:nbs), mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+230)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(nbs_index(1)+1:nbs_index(2)),             &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+231)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(nbs_index(2)+1:nbs_index(3)),             &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+232)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(nbs_index(3)+1:nbs_index(4)),             &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+233)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(nbs_index(4)+1:nbs_index(5)),             &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+234)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(nbs_index(5)+1:nbs_index(6)),             &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+235)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(nbs_index(6)+1:nbs_index(7)),             &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+236)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_cn_3nm(1:nbs), mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+237)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_ccn_1(1:nbs), mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+238)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_ccn_2(1:nbs), mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+239)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_ccn_3(1:nbs), mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+240)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_cdn(1:nbs), mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+241)
          CALL select_array_segment(nbs, 1, 1, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 1)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+242)
          CALL select_array_segment(nbs, 1, 1, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 3)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+243)
          CALL select_array_segment(nbs, 1, 1, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 6)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+244)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol_wat(1:nbs),                                 &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+245)
          CALL select_array_segment(nbs, 1, 2, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 1)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+246)
          CALL select_array_segment(nbs, 1, 2, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 2)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+247)
          CALL select_array_segment(nbs, 1, 2, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 3)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+248)
          CALL select_array_segment(nbs, 1, 2, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 6)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+249)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol_wat(nbs_index(1)+1:nbs_index(2)),           &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+250)
          CALL select_array_segment(nbs, 1, 3, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 1)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+251)
          CALL select_array_segment(nbs, 1, 3, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 2)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+252)
          CALL select_array_segment(nbs, 1, 3, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 3)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+253)
          CALL select_array_segment(nbs, 1, 3, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 4)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+254)
          CALL select_array_segment(nbs, 1, 3, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 5)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+255)
          CALL select_array_segment(nbs, 1, 3, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 6)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+256)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol_wat(nbs_index(2)+1:nbs_index(3)),           &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+257)
          CALL select_array_segment(nbs, 1, 4, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 1)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+258)
          CALL select_array_segment(nbs, 1, 4, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 2)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+259)
          CALL select_array_segment(nbs, 1, 4, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 3)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+260)
          CALL select_array_segment(nbs, 1, 4, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 4)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+261)
          CALL select_array_segment(nbs, 1, 4, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 5)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+262)
          CALL select_array_segment(nbs, 1, 4, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 6)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+263)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol_wat(nbs_index(3)+1:nbs_index(4)),           &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+264)
          CALL select_array_segment(nbs, 1, 5, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 2)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+265)
          CALL select_array_segment(nbs, 1, 5, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 3)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+266)
          CALL select_array_segment(nbs, 1, 6, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 5)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+267)
          CALL select_array_segment(nbs, 1, 7, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 5)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+268:item1_mode_diags+284)
          ! Do nothing - these are used by ukca_activate
        CASE (item1_mode_diags+345)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixsink(1:nbs),                               &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+346)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixsink(nbs_index(1)+1:nbs_index(2)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+347)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixsink(nbs_index(2)+1:nbs_index(3)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+348)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixsink(nbs_index(3)+1:nbs_index(4)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+349)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixsink(nbs_index(4)+1:nbs_index(5)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+350)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixsink(nbs_index(5)+1:nbs_index(6)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+351)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixsink(nbs_index(6)+1:nbs_index(7)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+352)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixflag(1:nbs),                               &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+353)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixflag(nbs_index(1)+1:nbs_index(2)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+354)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixflag(nbs_index(2)+1:nbs_index(3)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+355)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixflag(nbs_index(3)+1:nbs_index(4)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+356)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixflag(nbs_index(4)+1:nbs_index(5)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+357)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixflag(nbs_index(5)+1:nbs_index(6)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+358)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixflag(nbs_index(6)+1:nbs_index(7)),         &
                          mode_diags(1, 1, 1, k))
          !!
          !!  Nitrate diags
          !!
        CASE (item1_mode_diags+373)
          !           Do nothing, nucleation-mode NH4 emissions not included
        CASE (item1_mode_diags+374)
          CALL select_array_segment(nbs, 0, nmasprimnhaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+375)
          CALL select_array_segment(nbs, 0, nmasprimnhaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+376)
          CALL select_array_segment(nbs, 0, nmasprimnhcorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+377)
          !           Do nothing, nucleation-mode NO3 emissions not included
        CASE (item1_mode_diags+378)
          CALL select_array_segment(nbs, 0, nmasprimntaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+379)
          CALL select_array_segment(nbs, 0, nmasprimntaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+380)
          CALL select_array_segment(nbs, 0, nmasprimntcorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+381)
          CALL select_array_segment(nbs, 0, nmascondnnaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+382)
          CALL select_array_segment(nbs, 0, nmascondnncorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+383)
          !           Do nothing, nucleation mode not on for NH4/NO3
          !           Functionality added as placeholder
        CASE (item1_mode_diags+384)
          CALL select_array_segment(nbs, 0, nmasddepnhaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+385)
          CALL select_array_segment(nbs, 0, nmasddepnhaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+386)
          CALL select_array_segment(nbs, 0, nmasddepnhcorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+387)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+388)
          CALL select_array_segment(nbs, 0, nmasddepntaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+389)
          CALL select_array_segment(nbs, 0, nmasddepntaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+390)
          CALL select_array_segment(nbs, 0, nmasddepntcorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+391)
          CALL select_array_segment(nbs, 0, nmasddepnnaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+392)
          CALL select_array_segment(nbs, 0, nmasddepnncorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                           mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+393)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+394)
          CALL select_array_segment(nbs, 0, nmasnuscnhaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+395)
          CALL select_array_segment(nbs, 0, nmasnuscnhaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+396)
          CALL select_array_segment(nbs, 0, nmasnuscnhcorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+397)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+398)
          CALL select_array_segment(nbs, 0, nmasnuscntaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+399)
          CALL select_array_segment(nbs, 0, nmasnuscntaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+400)
          CALL select_array_segment(nbs, 0, nmasnuscntcorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+401)
          CALL select_array_segment(nbs, 0, nmasnuscnnaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+402)
          CALL select_array_segment(nbs, 0, nmasnuscnncorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+403)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+404)
          CALL select_array_segment(nbs, 0, nmasimscnhaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+405)
          CALL select_array_segment(nbs, 0, nmasimscnhaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+406)
          CALL select_array_segment(nbs, 0, nmasimscnhcorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+407)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+408)
          CALL select_array_segment(nbs, 0, nmasimscntaitsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+409)
          CALL select_array_segment(nbs, 0, nmasimscntaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+410)
          CALL select_array_segment(nbs, 0, nmasimscntcorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+411)
          CALL select_array_segment(nbs, 0, nmasimscnnaccsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+412)
          CALL select_array_segment(nbs, 0, nmasimscnncorsol, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+413)
          CALL select_array_segment(nbs, 0, nmasprocnhintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+414)
          CALL select_array_segment(nbs, 0, nmasprocntintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+415:item1_mode_diags+420)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+421)
          CALL select_array_segment(nbs, 0, nmascoagnhintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+422)
          CALL select_array_segment(nbs, 0, nmascoagnhintr24, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+423)
          CALL select_array_segment(nbs, 0, nmascoagntintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+424)
          CALL select_array_segment(nbs, 0, nmascoagntintr24, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+425)
          CALL select_array_segment(nbs, 0, nmascoagnhintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+426)
          CALL select_array_segment(nbs, 0, nmascoagntintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+427)
          CALL select_array_segment(nbs, 0, nmascoagnnintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+428)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+429)
          CALL select_array_segment(nbs, 0, nmasmergnhintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+430)
          CALL select_array_segment(nbs, 0, nmasmergnhintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+431)
          !           Do nothing, nucleation mode not on for NH4/NO3
        CASE (item1_mode_diags+432)
          CALL select_array_segment(nbs, 0, nmasmergntintr23, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+433)
          CALL select_array_segment(nbs, 0, nmasmergntintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+434)
          CALL select_array_segment(nbs, 0, nmasmergnnintr34, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
          !!! Partial volumes for RADAER
        CASE (item1_mode_diags+435)
          !           Do nothing, nucleation mode not on for NH4
        CASE (item1_mode_diags+436)
          CALL select_array_segment(nbs, 1, 2, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 9)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+437)
          CALL select_array_segment(nbs, 1, 3, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 9)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+438)
          CALL select_array_segment(nbs, 1, 4, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 9)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+439)
          !           Do nothing, nucleation mode not on for NO3
        CASE (item1_mode_diags+440)
          CALL select_array_segment(nbs, 1, 2, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 7)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+441)
          CALL select_array_segment(nbs, 1, 3, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 7)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+442)
          CALL select_array_segment(nbs, 1, 4, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 7)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+443)
          CALL select_array_segment(nbs, 1, 3, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 8)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+444)
          CALL select_array_segment(nbs, 1, 4, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 8)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
          !!
          !!  Dust 3rd insol mode diags
          !!
        CASE (item1_mode_diags+474)
          !           Do nothing, emissions not handled here now
        CASE (item1_mode_diags+475)
          CALL select_array_segment(nbs, 0, nmasddepdusupins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+476)
          CALL select_array_segment(nbs, 0, nmasnuscdusupins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+477)
          CALL select_array_segment(nbs, 0, nmasimscdusupins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+478)
          CALL select_array_segment(nbs, 0, nmascondsusupins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+479)
          CALL select_array_segment(nbs, 0, nmascondocsupins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+480)
          CALL select_array_segment(nbs, 0, nmascondsosupins, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
          !! Code for heterogeneous oxidation of SO2 --> SO4 on dust not yet in
        CASE (item1_mode_diags+481)
          mode_diags(:, :, :, k) = 0.0
        CASE (item1_mode_diags+482)
          CALL select_array_segment(nbs, 0, nmascoagsuintr18, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+483)
          CALL select_array_segment(nbs, 0, nmascoagocintr18, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+484)
          CALL select_array_segment(nbs, 0, nmascoagsointr18, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+485)
          CALL select_array_segment(nbs, 0, nmasagedduintr63, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+486)
          CALL select_array_segment(nbs, 0, nmasagedduintr74, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+487)
          CALL select_array_segment(nbs, 0, nmasagedduintr84, i_start, i_end)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_bud_aer_mas(i_start:i_end),                      &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+488)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(nbs_index(7)+1:nbs_index(8)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+489)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_sarea(nbs_index(7)+1:nbs_index(8)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+490)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_vconc(nbs_index(7)+1:nbs_index(8)),              &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+491)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(nbs_index(7)+1:nbs_index(8)),             &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+492)
          CALL select_array_segment(nbs, 1, 8, i_start_cp, i_end_cp,           &
                                    dim2_len = nmodes, dim3_index = 5)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol(i_start_cp:i_end_cp),                       &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+493)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixsink(nbs_index(7)+1:nbs_index(8)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+494)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_mdtfixflag(nbs_index(7)+1:nbs_index(8)),         &
                          mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+499)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_ccn_4(1:nbs), mode_diags(1, 1, 1, k))
        CASE (item1_mode_diags+500)
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_ccn_5(1:nbs), mode_diags(1, 1, 1, k))
        CASE DEFAULT
          cmessage = ' Item not found in CASE statement'
          CALL ereport('UKCA_AERO_CTL',UkcaD1codes(n)%item,cmessage)
        END SELECT
        IF (verbose == 2) THEN
          WRITE(umMessage, '(A16,3i4,3e14.4)') 'UKCA_MODE diag: ',             &
                                               k, n, ukcaD1codes(n)%item,      &
                                               SUM(mode_diags(:, :, :, k)),    &
                                               MAXVAL(mode_diags(:, :, :, k)), &
                                               MINVAL(mode_diags(:, :, :, k))
          CALL umPrint(umMessage, src='ukca_aero_ctl')
        END IF
      END IF    ! section == MODE_diag_sec etc
    END DO      ! nukca_d1items
  END IF        ! L_ukca_mode_diags

  SELECT CASE (glomap_config%i_ukca_activation_scheme)
  CASE ( i_ukca_activation_arg )
    ! Fill the drydiam array as input to UKCA_ACTIVATE
    DO imode = 1, nmodes
      IF ( mode(imode) ) THEN
        CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                  &
                        seg_drydp(nbs_index(imode-1)+1:nbs_index(imode)),      &
                        drydiam(1, 1, 1, imode) )
      END IF
    END DO
  CASE ( i_ukca_activation_jones )
    ! set ukca_cdnc prognostic without activation scheme
    ! change units from cm^-3 to m^-3
    i = name2ntpindex('cdnc      ')
    DO jl = 1, nbs
      seg_v1d_tmp(jl) = 1.0e+6 * seg_cdn(jl)
    END DO
    CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                      &
                    seg_v1d_tmp, all_ntp(i)%data_3d(1, 1, 1))
  END SELECT

  ! Calculate total aerosol surface area for use in UKCA chemistry for
  ! heterogeneous reactions by summing across all soluble modes and
  ! store result in all_ntp structure
  DO jl = 1, nbs
    seg_v1d_tmp(jl) = 0.0
  END DO
  DO imode = mode_nuc_sol, mode_cor_sol ! only take the soluble modes
    ! SAREA was converted to cm^2/cm^3 when it was calculated (see above)
    DO jl = 1, nbs
      seg_v1d_tmp(jl) = seg_v1d_tmp(jl) + seg_sarea(nbs_index(imode-1) + jl)
    END DO
  END DO
  i = name2ntpindex('surfarea  ')
  CALL insert_seg(lb, ncs, nbs, stride_s, model_levels, seg_v1d_tmp,           &
                  all_ntp(i)%data_3d(1, 1, 1))

  ! Fill the non-transported prognostic fields for RADAER coupling,
  ! depending on selected modes and components. The nucleation mode
  ! is not considered.
  ! --------------------------------------------------------------

  IF (glomap_config%l_ukca_radaer) THEN
    errcode = 0
    DO imode = mode_ait_sol,nmodes
      i_start = nbs_index(imode - 1) + 1
      CALL select_array_segment(nbs, 1, imode, i_start, i_end)
      IF (mode(imode)) THEN
        SELECT CASE(imode)
        CASE (mode_ait_sol)
          i = name2ntpindex('drydiam_ait_sol     ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(i_start:i_end),                            &
                          all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('wetdiam_ait_sol     ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_wetdp(i_start:i_end),                            &
                          all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('aerdens_ait_sol     ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(i_start:i_end),                           &
                          all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('pvol_h2o_ait_sol    ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol_wat(i_start:i_end),                         &
                          all_ntp(i)%data_3d(1, 1, 1))
        CASE (mode_acc_sol)
          i = name2ntpindex('drydiam_acc_sol     ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(i_start:i_end),                            &
                          all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('wetdiam_acc_sol     ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_wetdp(i_start:i_end),                            &
                          all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('aerdens_acc_sol     ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(i_start:i_end),                           &
                          all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('pvol_h2o_acc_sol    ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol_wat(i_start:i_end),                         &
                          all_ntp(i)%data_3d(1, 1, 1))
        CASE (mode_cor_sol)
          i = name2ntpindex('drydiam_cor_sol     ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(i_start:i_end),                            &
                          all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('wetdiam_cor_sol     ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_wetdp(i_start:i_end),                            &
                          all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('aerdens_cor_sol     ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(i_start:i_end),                           &
                          all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('pvol_h2o_cor_sol    ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_pvol_wat(i_start:i_end),                         &
                          all_ntp(i)%data_3d(1, 1, 1))
        CASE (mode_ait_insol)
          i = name2ntpindex('drydiam_ait_insol   ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(i_start:i_end),                            &
                          all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('aerdens_ait_insol   ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(i_start:i_end),                           &
                          all_ntp(i)%data_3d(1, 1, 1))
        CASE (mode_acc_insol)
          i = name2ntpindex('drydiam_acc_insol   ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                           seg_drydp(i_start:i_end),                           &
                           all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('aerdens_acc_insol   ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(i_start:i_end),                           &
                          all_ntp(i)%data_3d(1, 1, 1))
        CASE (mode_cor_insol)
          i = name2ntpindex('drydiam_cor_insol   ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_drydp(i_start:i_end),                            &
                          all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('aerdens_cor_insol   ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(i_start:i_end),                           &
                          all_ntp(i)%data_3d(1, 1, 1))
        CASE (mode_sup_insol)
          i = name2ntpindex('drydiam_sup_insol   ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                           seg_drydp(i_start:i_end),                           &
                           all_ntp(i)%data_3d(1, 1, 1))
          i = name2ntpindex('aerdens_sup_insol   ')
          CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                &
                          seg_rhopar(i_start:i_end),                           &
                          all_ntp(i)%data_3d(1, 1, 1))
        CASE DEFAULT
          cmessage = ' Mode not found in RADAER coupling CASE statement'
          errcode = ABS(imode)
          CALL ereport(RoutineName,errcode,cmessage)
        END SELECT
        DO icp = 1, ncp
          IF (component(imode, icp)) THEN
            CALL select_array_segment(nbs, 1, imode, i_start_cp, i_end_cp,     &
                                      dim2_len = nmodes, dim3_index = icp)
            IF (imode == mode_ait_sol) THEN
              SELECT CASE (icp)
              CASE (cp_su)
                i = name2ntpindex('pvol_su_ait_sol    ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_bc)
                i = name2ntpindex('pvol_bc_ait_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_oc)
                i = name2ntpindex('pvol_oc_ait_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_so)
                i = name2ntpindex('pvol_so_ait_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_no3)
                i = name2ntpindex('pvol_no3_ait_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_nh4)
                i = name2ntpindex('pvol_nh4_ait_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
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
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_bc)
                i = name2ntpindex('pvol_bc_acc_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_oc)
                i = name2ntpindex('pvol_oc_acc_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_cl)
                i = name2ntpindex('pvol_ss_acc_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_no3)
                i = name2ntpindex('pvol_no3_acc_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_nh4)
                i = name2ntpindex('pvol_nh4_acc_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_nn)
                i = name2ntpindex('pvol_nn_acc_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_du)
                i = name2ntpindex('pvol_du_acc_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_so)
                i = name2ntpindex('pvol_so_acc_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
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
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_bc)
                i = name2ntpindex('pvol_bc_cor_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_oc)
                i = name2ntpindex('pvol_oc_cor_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_cl)
                i = name2ntpindex('pvol_ss_cor_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_no3)
                i = name2ntpindex('pvol_no3_cor_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_nh4)
                i = name2ntpindex('pvol_nh4_cor_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_nn)
                i = name2ntpindex('pvol_nn_cor_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_du)
                i = name2ntpindex('pvol_du_cor_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_so)
                i = name2ntpindex('pvol_so_cor_sol     ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
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
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE (cp_oc)
                i = name2ntpindex('pvol_oc_ait_insol   ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
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
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE DEFAULT
                cmessage = ' Component not found in RADAER coupling CASE'//    &
                           ' statement'
                errcode = ABS(imode*100) + ABS(icp)
                CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
              END SELECT
            ELSE IF (imode == mode_cor_insol) THEN
              SELECT CASE (icp)
              CASE (cp_du)
                i = name2ntpindex('pvol_du_cor_insol   ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE DEFAULT
                cmessage = ' Component not found in RADAER coupling CASE'//    &
                           ' statement'
                errcode = ABS(imode*100) + ABS(icp)
                CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
              END SELECT
            ELSE IF (imode == mode_sup_insol) THEN
              SELECT CASE (icp)
              CASE (cp_du)
                i = name2ntpindex('pvol_du_sup_insol   ')
                CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,          &
                                seg_pvol(i_start_cp:i_end_cp),                 &
                                all_ntp(i)%data_3d(1, 1, 1))
              CASE DEFAULT
                cmessage = ' Component not found in RADAER coupling CASE'//    &
                           ' statement'
                errcode = ABS(imode*100) + ABS(icp)
                CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
              END SELECT
            ELSE
              cmessage = ' mode out of range in RADAER coupling IF clause'
              errcode = ABS(imode)
              CALL ereport('UKCA_AERO_CTL', errcode, cmessage)
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
    CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                      &
                    seg_drydp(nbs_index(mode_nuc_sol-1)+1:                     &
                              nbs_index(mode_nuc_sol)),                        &
                    all_ntp(i)%data_3d(1, 1, 1))

  END IF

  IF (l_ukca_cmip6_diags .OR. l_ukca_pm_diags) THEN
    ! Fill mdwat_diag array
    DO imode = 1, nmodes
      IF (mode(imode)) THEN
        CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                  &
                        seg_mdwat(nbs_index(imode-1)+1:nbs_index(imode)),      &
                        mdwat_diag(1,imode))
      END IF
    END DO
  END IF

  IF (l_ukca_pm_diags) THEN
    ! Copy wet diameter for calculating PM10 and PM2.5 diagnostics
    DO imode = 1, nmodes
      IF (mode(imode)) THEN
        CALL insert_seg(lb, ncs, nbs, stride_s, model_levels,                  &
                        seg_wetdp(nbs_index(imode-1)+1:nbs_index(imode)),      &
                        wetdp_diag(1,imode))
      END IF
    END DO
  END IF

END DO  ! ik loop over segments

!$OMP END PARALLEL
!
! End of OpenMP

! Deallocate lookup tables in dust impaction scavenging routine
IF (l_dust_slinn_impc_scav) CALL ukca_impc_scav_dust_dealloc()

! postponed reporting of nbadmdt from before aero_step
DO imode = 1, nmodes
  IF (mode(imode) ) THEN
    IF ( verbose > 0) THEN
      sum_nbadmdt(imode) = SUM(nbadmdt_3d(:, :, :, imode))
      IF (sum_nbadmdt(imode) > 0 ) THEN
        ! Below print out total occurrences if > 0
        WRITE(umMessage, '(A55)') 'MDT<MDTMIN, ND=0:IMODE,MDTMIN,NBADMDT'//    &
                                  ' (after bl mix)'
        CALL umPrint(umMessage, src='ukca_aero_ctl')
        IF (verbose > 1) THEN
          WRITE(umMessage, '(I6,E12.3,I12)') imode, MINVAL(mdtmin(imode, :)),  &
                                             sum_nbadmdt(imode)
          CALL umPrint(umMessage, src='ukca_aero_ctl')
          ! Per level sum of nbadmdt in 4 groups
          jl = model_levels / 4
          j = 0
          DO i = 1, 4
            l = MIN(j+1+jl,model_levels)
            WRITE(umMessage, '((10000I6))')                                    &
                          (SUM(nbadmdt_3d(:, :, k, imode)), k=j+1,l)
            ! Enclosing format in brackets makes it repeat for all items
            CALL umPrint(umMessage, src='ukca_aero_ctl')
            j = j + jl
          END DO
        END IF   ! Verbose > 1
      END IF     ! sum(nbadmdt > 0)
    END IF       ! Verbose > 0
  END IF         ! mode is active
END DO           ! over modes

IF (verbose > 1) THEN

  WRITE(umMessage, '(A30)') ' Tracers at end of UKCA_MODE:'
  CALL umPrint(umMessage, src='ukca_aero_ctl')

  ! Only one model level in UKCA box model
  IF (ichem == 1) THEN
    DO i = 1, MIN(2, model_levels)           !model_levels
      DO j = 1, n_chemistry_tracers
        WRITE(umMessage, '(A10,I6,A10,I6)') 'Level: ',i, ' Tracer: ', j
        CALL umPrint(umMessage, src='ukca_aero_ctl')
        WRITE(umMessage, '(A20,3E12.3)') 'chemistry_tracers:',                 &
           MINVAL(all_tracers(: , :, i, j)),                                   &
           MAXVAL(all_tracers(:, :, i, j)),                                    &
           SUM(all_tracers(:, :, i, j)) /                                      &
               REAL(SIZE(all_tracers(:, :, i, j)))
        CALL umPrint(umMessage, src='ukca_aero_ctl')
      END DO
    END DO
  END IF
  DO i = 1, MIN(2, model_levels)           !model_levels
    DO j = 1, n_mode_tracers
      IF (mode_tracer_debug(j)) THEN
        iaer = n_chemistry_tracers + j
        WRITE(umMessage, '(A10,I4,A10,I4,A10)') 'Level: ', i, ' Tracer: ',     &
               j, all_tracer_names(iaer)
        CALL umPrint(umMessage, src='ukca_aero_ctl')
        WRITE(umMessage, '(A14,3E12.3)') 'mode_tracers: ',                     &
               MINVAL(all_tracers(:, :, i, iaer)),                             &
               MAXVAL(all_tracers(:, :, i, iaer)),                             &
               SUM(all_tracers(:, :, i, iaer))/                                &
               REAL(SIZE(all_tracers(:, :, i, iaer)))
        CALL umPrint(umMessage, src='ukca_aero_ctl')
      END IF
    END DO
    WRITE(umMessage, '(A30,I6)') 'Number of merges for Level: ', i
    CALL umPrint(umMessage, src='ukca_aero_ctl')
    DO j = 1, nmodes
      WRITE(umMessage, '(2I6)') j, SUM(n_merge_3d(:, :, i, j))
      CALL umPrint(umMessage, src='ukca_aero_ctl')
    END DO
  END DO      ! i

  WRITE(umMessage, '(A24,I9,E12.3)')                                           &
       'Total Number of merges=:',                                             &
       SUM(n_merge_3d), REAL(SUM(n_merge_3d))/REAL(SIZE(n_merge_3d))
  CALL umPrint(umMessage, src='ukca_aero_ctl')
  DO j = 1, nmodes
    WRITE(umMessage, '(2I9,E12.3)') j,                                         &
         SUM(n_merge_3d(:, :, :, j)),                                          &
         REAL(SUM(n_merge_3d(:, :, :, j)))/REAL(SIZE(n_merge_3d(:, :, :, j)))
    CALL umPrint(umMessage, src='ukca_aero_ctl')
  END DO
END IF      ! verbose

IF (firstcall) firstcall = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE ukca_aero_ctl
!
SUBROUTINE surface_to_seg(lbase, ncol, nb, stride,nl, a, b)
! This replicates the surface value up a column
IMPLICIT NONE
INTEGER,INTENT(IN) :: lbase   ! Address of base of first column in segment
INTEGER,INTENT(IN) :: ncol    ! Number of columns in this segment
INTEGER,INTENT(IN) :: nb      ! The number of boxes in this segment
INTEGER,INTENT(IN) :: stride  ! The number of columns on the MPI task
INTEGER,INTENT(IN) :: nl      ! The number of model levels
REAL,INTENT(IN) :: a(stride)  ! The 2D surface data from which a column is taken
REAL,INTENT(IN OUT) :: b(nb)  ! Segment made up from columns

! local loop iterators
INTEGER :: ic, l, ia, ib

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SURFACE_TO_SEG'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
ia = lbase       ! the address of base of 1st column on segment
ib = 1

DO ic = 1, ncol           ! loop over columns on the segment
  DO l = 1, nl            ! climb a column
    b(ib) = a(ia)         ! copy the data into the shorter vector
    ib = ib + 1           ! next location in segment
  END DO
  ia = lbase+ic           ! start of base of next column
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
END SUBROUTINE surface_to_seg
!
SUBROUTINE surface_to_seg_int(lbase, ncol, nb, stride,nl, a, b)
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

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SURFACE_TO_SEG_INT'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
ia = lbase       ! the address of base of 1st column on segment
ib = 1

DO ic = 1, ncol           ! loop over columns on the segment
  DO l = 1, nl            ! climb a column
    b(ib) = a(ia)         ! copy the data into the shorter vector
    ib = ib + 1           ! next location in segment
  END DO
  ia = lbase+ic           ! start of base of next column
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
END SUBROUTINE surface_to_seg_int
!
SUBROUTINE extract_seg(lbase, ncol, nb, stride, nl, a, b)
! This puts the elements of a 3D array and packs ito a segment (of columns)
IMPLICIT NONE
INTEGER,INTENT(IN) :: lbase   ! Address of base of first column in segment
INTEGER,INTENT(IN) :: ncol    ! Number of columns in this segment
INTEGER,INTENT(IN) :: nb      ! The number of boxes in this segment
INTEGER,INTENT(IN) :: stride  ! The number of columns on the MPI task
INTEGER,INTENT(IN) :: nl      ! The number of model levels
REAL,INTENT(IN) :: a(stride*nl)    ! the 3D data from which a column is taken
REAL,INTENT(IN OUT) :: b(nb)       ! segment made up from columns

! local loop iterators
INTEGER :: ic, l, ia, ib

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'EXTRACT_SEG'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
ia = lbase                ! the address of base of 1st column on segment
ib = 1

DO ic = 1, ncol           ! loop over columns on the segment
  DO l = 1, nl            ! climb a column
    b(ib) = a(ia)         ! copy the data into the shorter vector
    ib = ib + 1           ! next location in segment
    ia = ia + stride      ! next location in unrolled 3D array
  END DO
  ia = lbase + ic         ! start of base of next column
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
END SUBROUTINE extract_seg
!
SUBROUTINE insert_seg(lbase, ncol, nb, stride, nl, b, a)
! This subroutine puts the elements of the segment vector back into the 3D array
IMPLICIT NONE
INTEGER,INTENT(IN)  :: lbase   ! Address of base of first column in segment
INTEGER,INTENT(IN)  :: ncol    ! Number of columns in this segment
INTEGER,INTENT(IN)  :: nb      ! The number of boxes in this segment
INTEGER,INTENT(IN)  :: stride  ! The number of columns on the MPI task
INTEGER,INTENT(IN)  :: nl      ! The number of model levels
REAL,INTENT(IN)     :: b(nb)   ! the incoming segment
REAL,INTENT(IN OUT) :: a(stride*nl) ! 3D array to be modified by b()

! local loop iterators
INTEGER :: ic, l, ia, ib

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

END SUBROUTINE insert_seg
!
SUBROUTINE int_ins_seg(lbase, ncol, nb, stride, nl, b, a)
! This subroutine puts the elements of the segment vector back into the 3D array
IMPLICIT NONE
INTEGER,INTENT(IN) :: lbase   ! Address of base of first column in segment
INTEGER,INTENT(IN) :: ncol    ! Number of columns in this segment
INTEGER,INTENT(IN) :: nb      ! The number of boxes in this segment
INTEGER,INTENT(IN) :: stride  ! The number of columns on the MPI task
INTEGER,INTENT(IN) :: nl      ! The number of model levels
INTEGER(KIND=integer_32) ,INTENT(IN)     :: b(nb)        ! incoming segment
INTEGER(KIND=integer_32) ,INTENT(IN OUT) :: a(stride*nl) ! 3D array unwound

! local loop iterators
INTEGER :: ic, l, ia, ib
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'INT_INS_SEG'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ia = lbase                ! the address of base of 1st column on segment
ib = 1                    ! the box id in the segment

DO ic = 1, ncol           ! loop over columns on the segment
  DO l = 1, nl            ! climb a column
    a(ia) = b(ib)         ! copy (spreading) the data into the 3d array
    ib = ib + 1           ! next location in segment
    ia = ia + stride      ! next location in unrolled 3D array
  END DO
  ia = lbase + ic         ! start of base of next column
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
END SUBROUTINE int_ins_seg
!
SUBROUTINE select_array_segment(dim1_len, dim2_start, dim2_index,              &
                                i_start, i_end, dim2_len, dim3_index)
! To accommodate a chunking size (dim1_len) which is often smaller for the
! final remainder chunk, 2D and 3D arrays are put into a 1D array so they
! can form contiguous data even when the chunking size (dim1_len) is reduced
! in size.
! For an example, consider a 2D array of size (nbox, nmodes), where the only
! part of the array which is used is (1:nbs, 1:nmodes). When nbs=nbox the
! data is contiguous, with no gaps. However when nbs < nbox, there are gaps
! in the data at (nbs+1:nbox, :). This a problem when passing to subroutines
! like ukca_aero_step which aren't expecting these gaps.
! When the data is put into a 1D array of size (nbox * nmodes) the data
! can be sequeezed together so that it is contiguous, and the only gap,
! (nbs*nmodes+1:nbox*nmodes), is at the end of the array and is simply not
! used.
! When wanting to access (:, i) from the 2D array, this routine indicates
! the i_start and i_end in the form (i_start:i_end) which gives the same
! data in the 1D array. And it does something similar for 3D arrays.

IMPLICIT NONE

! Subroutine arguments
! Length of first dimension (nbs in ukca_aero_ctl, which is less than or
! equal to nbox)
INTEGER, INTENT(IN)           :: dim1_len
! 0 if 2nd dimension starts at 0, and 1 if 2nd dimension starts at 1
INTEGER, INTENT(IN)           :: dim2_start
! Index of 2nd dimension
INTEGER, INTENT(IN)           :: dim2_index
! Start point in the 1D array
INTEGER, INTENT(OUT)          :: i_start
! End point in the 1D array
INTEGER, INTENT(OUT)          :: i_end
! Size of 2nd dimension (needed if the original array was 3D)
INTEGER, INTENT(IN), OPTIONAL :: dim2_len
! Index of 3rd dimension
INTEGER, INTENT(IN), OPTIONAL :: dim3_index

! In the case of a 3rd dimension
IF (PRESENT(dim2_len) .AND. PRESENT(dim3_index)) THEN
  i_end = dim1_len * dim2_len * (dim3_index - 1)
ELSE
  i_end = 0
END IF

! The end of the data
i_end = i_end + (dim1_len * (dim2_index - dim2_start + 1))

! Starting point of data
i_start = i_end - dim1_len + 1

END SUBROUTINE select_array_segment
!
END MODULE ukca_aero_ctl_mod
