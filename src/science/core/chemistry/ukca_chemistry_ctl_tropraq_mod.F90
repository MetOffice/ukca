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
!  Main driver routine for chemistry under tropospheric and regional air
!  quality configurations.
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
MODULE ukca_chemistry_ctl_tropraq_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName =                           &
                                        'UKCA_CHEMISTRY_CTL_TROPRAQ_MOD'

CONTAINS

SUBROUTINE ukca_chemistry_ctl_tropraq(                                         &
                row_length, rows, model_levels, theta_field_size, tot_n_pnts,  &
                ntracers,                                                      &
                secs_per_step,                                                 &
                pres, temp, q,                                                 &
                qcf, qcl, rh,                                                  &
                tracer,                                                        &
                all_ntp,                                                       &
                cloud_frac,                                                    &
                photol_rates,                                                  &
                volume,                                                        &
                so4_aitken, so4_accum, soot_fresh, soot_aged,                  &
                ocff_fresh, ocff_aged, biogenic,                               &
                sea_salt_film, sea_salt_jet,                                   &
                uph2so4inaer,                                                  &
                delso2_wet_h2o2,                                               &
                delso2_wet_o3,                                                 &
                delh2so4_chem,                                                 &
                delso2_drydep,                                                 &
                delso2_wetdep,                                                 &
                trop_ch4_mol,                                                  &
                trop_o3_mol,                                                   &
                trop_oh_mol,                                                   &
                strat_ch4_mol,                                                 &
                strat_ch4loss,                                                 &
                len_stashwork,                                                 &
                stashwork,                                                     &
                H_plus,                                                        &
                zdryrt, zwetrt, nlev_with_ddep, L_stratosphere,                &
                firstcall                                                      &
                )

USE ukca_um_legacy_mod,   ONLY: rgas => r
USE asad_mod,             ONLY: advt, jpctr, jpcspf, jpro2, jpdd, jpdw,        &
                                jpeq, jphk, jppj, jpspec, jpnr, jpspj, jpspt,  &
                                jptk, ldepd, ldepw, nadvt, nnaf, nprkx, ntrkx, &
                                speci, sph2o, spj, spt, spro2, ctype, y, nlnaro2
USE ukca_config_defs_mod, ONLY: nr_therm, nr_phot
USE ukca_cspecies,        ONLY: c_species, n_ch4, n_hono2, n_o3,               &
                                nn_ch4, nn_cl, nn_h2o2, nn_h2so4,              &
                                nn_o1d, nn_o3, nn_o3p, nn_oh,                  &
                                n_h2o, nn_so2, c_na_species,                   &
                                n_h2so4
USE ukca_constants,       ONLY: c_h2o, c_hono2
USE ukca_config_constants_mod, ONLY: avogadro, boltzmann
USE ukca_config_specification_mod, ONLY: ukca_config

USE ukca_raq_diags_mod,   ONLY: ukca_raq_diags
USE ukca_ntp_mod,         ONLY: ntp_type, dim_ntp, name2ntpindex
USE ukca_chemco_raq_mod,  ONLY: ukca_chemco_raq
USE ukca_deriv_raqaero_mod, ONLY: ukca_deriv_raqaero

USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
USE ereport_mod,          ONLY: ereport
USE umPrintMgr,           ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper

USE ukca_missing_data_mod, ONLY: rmdi

USE errormessagelength_mod, ONLY: errormessagelength

USE ukca_be_drydep_mod,   ONLY: ukca_be_drydep
USE ukca_be_wetdep_mod,   ONLY: ukca_be_wetdep
USE ukca_ch4_stratloss_mod, ONLY: ukca_ch4_stratloss
USE ukca_chemco_mod,      ONLY: ukca_chemco
USE ukca_deriv_mod,       ONLY: ukca_deriv
USE ukca_deriv_aero_mod,  ONLY: ukca_deriv_aero
USE ukca_deriv_raq_mod,   ONLY: ukca_deriv_raq
USE ukca_fracdiss_mod,    ONLY: ukca_fracdiss

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length        ! size of UKCA x dimension
INTEGER, INTENT(IN) :: rows              ! size of UKCA y dimension
INTEGER, INTENT(IN) :: model_levels      ! size of UKCA z dimension
INTEGER, INTENT(IN) :: theta_field_size  ! no. of points in horizontal
INTEGER, INTENT(IN) :: tot_n_pnts        ! no. of points in full domain
INTEGER, INTENT(IN) :: ntracers          ! no. of tracers
INTEGER, INTENT(IN) :: uph2so4inaer      ! flag for H2SO4 updating
INTEGER, INTENT(IN) :: nlev_with_ddep(theta_field_size) ! No levs in bl

REAL, INTENT(IN) :: secs_per_step        ! time step
REAL, INTENT(IN) :: pres(tot_n_pnts)     ! pressure
REAL, INTENT(IN) :: temp(tot_n_pnts)     ! actual temperature
REAL, INTENT(IN) :: volume(tot_n_pnts)   ! cell volume
REAL, INTENT(IN) :: H_plus(tot_n_pnts)   ! pH array

! Aerosol MMR / numbers from CLASSIC, used for calculation
! of surface area if heterogeneous chemistry is ON.
! Aerosol MMR for most aerosol types (kg kg-1)
REAL, INTENT(IN) :: so4_aitken    (tot_n_pnts)
REAL, INTENT(IN) :: so4_accum     (tot_n_pnts)
REAL, INTENT(IN) :: soot_fresh    (tot_n_pnts)
REAL, INTENT(IN) :: soot_aged     (tot_n_pnts)
REAL, INTENT(IN) :: ocff_fresh    (tot_n_pnts)
REAL, INTENT(IN) :: ocff_aged     (tot_n_pnts)
REAL, INTENT(IN) :: biogenic      (tot_n_pnts)
! Aerosol numbers for sea-salt (m-3)
REAL, INTENT(IN) :: sea_salt_film (tot_n_pnts)
REAL, INTENT(IN) :: sea_salt_jet  (tot_n_pnts)

REAL, INTENT(IN) :: qcf(tot_n_pnts)
REAL, INTENT(IN) :: qcl(tot_n_pnts)
REAL, INTENT(IN) :: rh(tot_n_pnts)                  ! relative humidity frac
REAL, INTENT(IN) :: cloud_frac(tot_n_pnts)
REAL, INTENT(IN) :: zdryrt(theta_field_size,jpdd)              ! dry dep rate
REAL, INTENT(IN) :: zwetrt(theta_field_size,model_levels,jpdw) ! wet dep rate
REAL, INTENT(IN) :: photol_rates(tot_n_pnts,jppj)
REAL, INTENT(IN OUT) :: q(tot_n_pnts)               ! water vapour
REAL, INTENT(IN OUT) :: tracer(tot_n_pnts,ntracers) ! tracer MMR

! SO2 increments
REAL, INTENT(IN OUT) :: delSO2_wet_H2O2(tot_n_pnts)
REAL, INTENT(IN OUT) :: delSO2_wet_O3(tot_n_pnts)
REAL, INTENT(IN OUT) :: delh2so4_chem(tot_n_pnts)
REAL, INTENT(IN OUT) :: delSO2_drydep(tot_n_pnts)
REAL, INTENT(IN OUT) :: delSO2_wetdep(tot_n_pnts)

! Trop CH4 burden (moles)
REAL, INTENT(IN OUT) :: trop_ch4_mol(tot_n_pnts)

! Trop O3 burden (moles)
REAL, INTENT(IN OUT) :: trop_o3_mol(tot_n_pnts)

! Trop OH burden (moles)
REAL, INTENT(IN OUT) :: trop_oh_mol(tot_n_pnts)

! Strat CH4 burden (moles)
REAL, INTENT(IN OUT) :: strat_ch4_mol(tot_n_pnts)

! Strat CH4 loss (Moles/s)
REAL, INTENT(IN OUT) :: strat_ch4loss(tot_n_pnts)

! Non transported prognostics
TYPE(ntp_type), INTENT(IN OUT) :: all_ntp(dim_ntp)

! Diagnostics array
INTEGER, INTENT(IN) :: len_stashwork
REAL, INTENT(IN OUT) :: stashwork(len_stashwork)

! Stratosphere mask
LOGICAL, INTENT(IN) :: L_stratosphere(tot_n_pnts)

! Flag for determining if this is the first chemistry call
LOGICAL, INTENT(IN) :: firstcall

! Local variables
INTEGER :: nlev_ch4_stratloss  ! No top levs for CH4 stratospheric loss
INTEGER, SAVE :: nr            ! no of rxns for BE
INTEGER, SAVE :: n_be_calls    ! no of call to BE solver

INTEGER :: i             ! loop variable
INTEGER :: j             ! loop variable
INTEGER :: js            ! loop variable
INTEGER :: jtr           ! loop variable - transported tracers
INTEGER :: jro2          ! loop variable - NTP RO2 species
INTEGER :: jna           ! loop variable, non-advected species
INTEGER :: jspf          ! loop variable - all active chemical species in f
INTEGER :: k             ! loop variable
INTEGER :: l             ! loop variable

INTEGER :: kcs           ! start index of current model level
INTEGER :: kce           ! end index of current model level

INTEGER :: errcode                             ! Error code: ereport
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message

REAL, SAVE :: dts        ! B. Euler timestep

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
REAL :: zfrdiss(theta_field_size, model_levels, jpdw, jpeq+1)
REAL :: kp_nh(row_length, rows, model_levels)     ! Dissociation const
REAL :: BE_tnd(theta_field_size)                  ! total no density, molec cm-3
REAL :: BE_h2o(theta_field_size)                  ! water vapour concn
REAL :: BE_o2(theta_field_size)                   ! oxygen concn
REAL :: BE_vol(theta_field_size)                  ! gridbox volume
REAL :: BE_rho(theta_field_size)                  ! air density (kg m-3)
REAL :: BE_rh_frac(theta_field_size)              ! RH (fraction: 0.000-0.999)

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
REAL :: k_dms(theta_field_size,5)                 ! dms rate coeffs

!     Dry and wet deposition fluxes (mol s-1)
REAL :: dry_dep_3d(theta_field_size, model_levels, jpdd) ! 3d dry dep
REAL :: wet_dep_3d(theta_field_size, model_levels, jpdw) ! 3d wet dep

! Full ntp array
REAL :: ntp_data(tot_n_pnts,dim_ntp)

! 1-D masks for troposphere and NAT height limitation
LOGICAL :: stratflag(theta_field_size)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMISTRY_CTL_TROPRAQ'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
nr = nr_therm + nr_phot

IF (firstcall) THEN

  ! Backward Euler timestep variables

  n_be_calls = INT(secs_per_step/ukca_config%dts0)
  ! Ensure we call BE at least once
  IF (n_be_calls == 0) n_be_calls = 1
  ! Calculate the BE Timestep
  dts = secs_per_step/n_be_calls

  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'(A,I0,E12.4)') 'n_be_calls, dts= ',n_be_calls, dts
    CALL umPrint(umMessage,src='ukca_chemistry_ctl_tropraq')
  END IF

END IF  ! of initialization of chemistry subroutine (firstcall)

! Calculate dissolved fraction (only used in online chem with B-E Solver)
IF (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) THEN
  ! Send H_plus array to calculate fraction dissolved
  CALL ukca_fracdiss(row_length, rows, model_levels,                           &
                     temp, pres, rh, qcl, zfrdiss, kp_nh, H_plus)
END IF

! A stratospheric CH4 loss rate may be applied at a number of top levels.
! If the number of top levels is not set explicitly it is set to 3 as
! required by the UM. This is deprecated functionality. In future, the UM
! should set the number of levels explicitly via the ukca_setup call
! making this functionality redundant.
IF (ukca_config%nlev_ch4_stratloss >= 0) THEN
  nlev_ch4_stratloss = ukca_config%nlev_ch4_stratloss
ELSE
  nlev_ch4_stratloss = 3
END IF

DO l = 1, dim_ntp
  IF (all_ntp(l)%l_required) THEN
    ntp_data(:,l) = RESHAPE(all_ntp(l)%data_3d(:,:,:),[tot_n_pnts])
  END IF
END DO

!$OMP PARALLEL DEFAULT(NONE)                                                   &
!$OMP PRIVATE(BE_biogenic, BE_dryrt, BE_frdiss, BE_h2o, BE_hrc,                &
!$OMP         BE_o2, BE_ocff_aged, BE_ocff_fresh,                              &
!$OMP         BE_rc, BE_rh_frac, BE_rho,                                       &
!$OMP         BE_sea_salt_film, BE_sea_salt_jet,                               &
!$OMP         BE_so4_aitken, BE_so4_accum,                                     &
!$OMP         BE_soot_aged, BE_soot_fresh,                                     &
!$OMP         BE_tnd, BE_wetrt, BE_vol, BE_y,                                  &
!$OMP         errcode, jna, jro2, jspf, k, kcs, kce, k_dms, l,                 &
!$OMP         SO2_dryox_OH, SO2_wetox_H2O2, SO2_wetox_O3,                      &
!$OMP         stratflag, ystore,                                               &
!$OMP         zfnatr, zftr)                                                    &
!$OMP SHARED(advt, avogadro, biogenic, boltzmann, c_species, c_na_species,     &
!$OMP        cloud_frac, delh2so4_chem, delSO2_drydep, delSO2_wet_H2O2,        &
!$OMP        delSO2_wet_O3, delSO2_wetdep, dry_dep_3d, dts, speci,             &
!$OMP        jpctr, jpdd, jpdw, jphk, jppj, jpspec, jpro2,                     &
!$OMP        jpcspf, spro2, nlnaro2, ctype, l_stratosphere,                    &
!$OMP        ukca_config,                                                      &
!$OMP        ldepd, ldepw, model_levels,                                       &
!$OMP        n_be_calls, n_ch4, n_o3, nadvt, nlev_with_ddep, ntp_data,         &
!$OMP        nlev_ch4_stratloss, nn_ch4, nn_h2o2, nn_h2so4, nn_o1d, nn_o3,     &
!$OMP        nn_o3p, nn_oh, nn_so2, nnaf, nr, nr_therm,                        &
!$OMP        ocff_aged, ocff_fresh, photol_rates, pres, q, qcf, qcl, rgas, rh, &
!$OMP        row_length, rows, sea_salt_film, sea_salt_jet, so4_accum,         &
!$OMP        so4_aitken, soot_aged, soot_fresh, strat_ch4_mol, strat_ch4loss,  &
!$OMP        temp, theta_field_size, tracer, trop_ch4_mol,                     &
!$OMP        trop_o3_mol, trop_oh_mol, uph2so4inaer, volume, wet_dep_3d,       &
!$OMP        zdryrt, zfrdiss, zwetrt, H_plus, cmessage)

IF (.NOT. ALLOCATED(ystore) .AND. uph2so4inaer == 1)                           &
                             ALLOCATE(ystore(theta_field_size))

!$OMP DO SCHEDULE(DYNAMIC)
DO k=1,model_levels

  ! Determine the indices for extracting the slice of data which corresponds to
  ! the k-th model level
  kcs = (k - 1) * theta_field_size + 1
  kce = k * theta_field_size

  ! Copy water vapour and ice field into 1-D arrays
  IF (ukca_config%l_ukca_het_psc) THEN
    IF (k <= model_levels) THEN
      sph2o(:) = qcf(kcs:kce)/c_h2o
    ELSE
      sph2o(:) = 0.0
    END IF
  END IF

  ! Convert mmr into vmr for tracers. Pass data from the tracer 3D array,
  ! unwrap it and pass into the 1D zftr array before calling the chemistry
  ! solver. If running with nontransport RO2 species, data for the RO2 species
  ! needs to be passed from the ntp_data array as well.

  ! First set counter for all chemical species in f array
  jspf = 0
  ! Loop through all species
  DO js = 1,jpspec
    ! First try to map with transported tracer
    DO jtr=1,jpctr
      IF (advt(jtr) == speci(js)) THEN
        jspf = jspf+1
        ! Map data from tracer 3D array into zftr
        zftr(:,jspf) = tracer(kcs:kce,jtr)/c_species(jtr)
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
            zftr(:,jspf) = ntp_data(kcs:kce,l) / c_na_species(jna)
          ELSE
            errcode = jro2
            WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl_tropraq'
            CALL umPrint(umMessage,src='ukca_chemistry_ctl_tropraq')
            cmessage='ERROR: Indices for RO2 species do not match with nadvt'
            CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
          END IF ! Close IF nadvt and spro2 match

        END IF   ! Close IF RO2 species
      END DO     ! Close loop through RO2 species
    END IF       ! Close IF RO2_NTP
  END DO         ! Close loop through all species

  ! Check we have the correct number of active chemical species
  IF (ukca_config%l_ukca_ro2_ntp) THEN
    IF (jspf /= jpro2+jpctr) THEN
      errcode = jspf
      WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl_tropraq'
      CALL umPrint(umMessage,src='ukca_chemistry_ctl_tropraq')
      cmessage = 'ERROR: Number of chemical active species /= jpro2+jpctr'
      CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
    END IF
  ELSE
    IF (jspf /= jpctr) THEN
      errcode = jspf
      WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl_tropraq'
      CALL umPrint(umMessage,src='ukca_chemistry_ctl_tropraq')
      cmessage = 'ERROR: Number of chemical active species /= jpctr'
      CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
    END IF
  END IF

  !       Call ASAD routines to do chemistry integration

  !         Calculate total number  density, o2, h2o, and tracer
  !         concentrations for Backward Euler solver

  BE_tnd(:) = pres(kcs:kce) / (boltzmann * 1.0e6 * temp(kcs:kce))
  BE_o2(:)  = 0.2095 * BE_tnd(:)
  BE_h2o(:) = q(kcs:kce)*BE_tnd(:)
  BE_vol(:) = volume(kcs:kce)*1.0e6 !  m3->cm3

  ! Air density in kg m-3 (note that in this formula rgas is
  ! the gas constant for dry air:  287.05 J kg-1 K-1)
  BE_rho(:) = pres(kcs:kce) / ( rgas * temp(kcs:kce) )

  DO l=1,jpctr
    zftr(:,l) = zftr(:,l) * BE_tnd(:)
  END DO

  ! Update non-advected species (only for B-E solvers) from
  ! non transported prognostics
  IF (.NOT. ALLOCATED(zfnatr)) ALLOCATE(zfnatr(theta_field_size,nnaf))
  DO js = 1, nnaf
    ! Get index of nadvt(js) in NTP array.
    ! If not found this is a fatal error.
    l = name2ntpindex(nadvt(js))
    zfnatr(:,js) = ntp_data(kcs:kce,l) / c_na_species(js)
  END DO


  DO l=1,nnaf
    zfnatr(:,l) = zfnatr(:,l) * BE_tnd(:)
  END DO

  !         Assign wet and dry deposition rates to species

  CALL ukca_be_wetdep(theta_field_size, zwetrt(:,k,:), BE_wetrt)

  CALL ukca_be_drydep(k, theta_field_size, nlev_with_ddep, zdryrt, BE_dryrt)

  ! Assign fractional dissociation
  IF (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) THEN
    DO i=1,jpeq+1
      CALL ukca_be_wetdep(theta_field_size, zfrdiss(:,k,:,i), BE_frdiss(:,:,i))
    END DO
  END IF

  !         Calculate reaction rate coefficients

  IF (.NOT. ALLOCATED(BE_rc))                                                  &
       ALLOCATE(BE_rc(theta_field_size,nr_therm))

  IF (ukca_config%l_ukca_raq .OR. ukca_config%l_ukca_raqaero) THEN
    ! Set RH (fraction: 0 - 0.999), aerosol number (m-3) for the sea-salt
    ! modes and aerosol mmr (kg kg-1) for the other CLASSIC aerosol types,
    ! as well as heterogeneous rates (s-1) before the call to UKCA_CHEMCO_RAQ  .
    !
    BE_rh_frac      (:) = rh           (kcs:kce)
    BE_so4_aitken   (:) = so4_aitken   (kcs:kce)
    BE_so4_accum    (:) = so4_accum    (kcs:kce)
    BE_soot_fresh   (:) = soot_fresh   (kcs:kce)
    BE_soot_aged    (:) = soot_aged    (kcs:kce)
    BE_ocff_fresh   (:) = ocff_fresh   (kcs:kce)
    BE_ocff_aged    (:) = ocff_aged    (kcs:kce)
    BE_biogenic     (:) = biogenic     (kcs:kce)
    BE_sea_salt_film(:) = sea_salt_film(kcs:kce)
    BE_sea_salt_jet (:) = sea_salt_jet (kcs:kce)

    IF (.NOT. ALLOCATED(BE_hrc))                                               &
      ALLOCATE (BE_hrc (theta_field_size, jphk))

    ! Fill in stratospheric flag indicator, which will be passed
    ! to UKCA_CHEMCO_RAQ in order to set tropospheric heterogeneous
    ! reactions (if present) to zero in the stratosphere
    stratflag (:) = L_stratosphere(kcs:kce)

    CALL ukca_chemco_raq(nr_therm, theta_field_size, temp(kcs:kce), BE_tnd,    &
                         BE_h2o, BE_o2, qcl(kcs:kce), cloud_frac(kcs:kce),     &
                         BE_frdiss, BE_rho, BE_rh_frac, BE_so4_aitken,         &
                         BE_so4_accum, BE_soot_fresh, BE_soot_aged,            &
                         BE_ocff_fresh, BE_ocff_aged, BE_biogenic,             &
                         BE_sea_salt_film, BE_sea_salt_jet, stratflag,         &
                         BE_rc, BE_hrc, H_plus(kcs:kce))
  ELSE

    CALL ukca_chemco(nr_therm, theta_field_size, temp(kcs:kce), BE_tnd,        &
                     BE_h2o, BE_o2, qcl(kcs:kce), cloud_frac(kcs:kce),         &
                     BE_frdiss, k_dms, BE_rc, H_plus(kcs:kce))
  END IF

  !         Assign tracer concentrations to species concentrations

  BE_y(:,:) = 0.0
  DO i = 1,jpspec
    DO j = 1,jpctr
      IF (speci(i) == advt(j)) THEN
        BE_y(:,i) = zftr(:,j)
        EXIT
      END IF
    END DO
  END DO

  ! Assign non-advected concentrations to species concentrations
  DO i = 1,jpspec
    DO j = 1,nnaf
      IF (speci(i) == nadvt(j)) THEN
        BE_y(:,i) = zfnatr(:,j)
        EXIT
      END IF
    END DO
  END DO

  IF (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) THEN
    SO2_wetox_H2O2(:) = 0.0
    SO2_dryox_OH(:)   = 0.0
    SO2_wetox_O3(:)   = 0.0
  ELSE
    delSO2_wet_H2O2(kcs:kce) = 0.0
    delSO2_wet_O3(kcs:kce)   = 0.0
    delh2so4_chem(kcs:kce)   = 0.0
    delSO2_drydep(kcs:kce)   = 0.0
    delSO2_wetdep(kcs:kce)   = 0.0
  END IF

  !         Call Backward Euler solver
  !         N.B. Emissions already added, via call to TR_MIX from
  !         UKCA_EMISSION_CTL

  IF (ukca_config%l_ukca_aerchem) THEN

    CALL ukca_deriv_aero(nr_therm, n_be_calls, theta_field_size, BE_rc,        &
                         BE_wetrt, BE_dryrt, photol_rates(kcs:kce,:), k_dms,   &
                         BE_h2o, BE_tnd, BE_o2, dts, BE_y, SO2_wetox_H2O2,     &
                         SO2_wetox_O3, SO2_dryox_OH)
  ELSE IF (ukca_config%l_ukca_raq) THEN

    CALL ukca_deriv_raq(nr_therm, n_be_calls, theta_field_size, dts, BE_h2o,   &
                        BE_vol, BE_rc, BE_hrc, BE_dryrt, BE_wetrt,             &
                        photol_rates(kcs:kce,:), ldepd, ldepw, BE_y,           &
                        dry_dep_3d(:,k,:), wet_dep_3d(:,k,:))

  ELSE IF (ukca_config%l_ukca_raqaero) THEN

    ! Store H2SO4 tracer if it will be updated in MODE using delh2so4_chem
    IF (uph2so4inaer == 1) ystore(:) = BE_y(:,nn_h2so4)

    CALL ukca_deriv_raqaero(nr_therm, n_be_calls, theta_field_size, dts,       &
                            BE_rc, BE_wetrt, BE_dryrt,                         &
                            photol_rates(kcs:kce,:), BE_h2o, BE_y,             &
                            so2_wetox_H2O2, so2_wetox_O3, so2_dryox_OH)

    ! Restore H2SO4 tracer as it will be updated in MODE using delh2so4_chem
    IF (uph2so4inaer == 1) BE_y(:,nn_h2so4) = ystore(:)

  ELSE

    CALL ukca_deriv(nr, n_be_calls, theta_field_size, BE_rc, BE_wetrt,         &
                    BE_dryrt, photol_rates(kcs:kce,:), BE_h2o, BE_tnd,         &
                    BE_o2, dts, BE_y)
  END IF
  IF (ALLOCATED(BE_hrc)) DEALLOCATE(BE_hrc)
  IF (ALLOCATED(BE_rc))  DEALLOCATE(BE_rc)

  ! Apply stratospheric CH4 loss rate at top levels.
  IF (k >= model_levels + 1 - nlev_ch4_stratloss) THEN
    CALL ukca_ch4_stratloss(n_be_calls, theta_field_size, BE_vol, dts,         &
                            BE_y(:,nn_ch4), strat_ch4loss(kcs:kce))
  END IF

  DO j = 1,jpctr
    DO i = 1,jpspec
      IF (advt(j) == speci(i)) THEN
        zftr(:,j) = BE_y(:,i)
        EXIT
      END IF
    END DO
  END DO

  DO js = 1, jpctr
    tracer(kcs:kce,js) = zftr(:,js)/BE_tnd(:)*c_species(js)
  END DO

  ! Convert non-advected tracers back to vmr
  DO j = 1,nnaf
    DO i = 1,jpspec
      IF (nadvt(j) == speci(i)) THEN
        zfnatr(:,j) = BE_y(:,i)/BE_tnd(:)
        EXIT
      END IF
    END DO
  END DO

  ! Set the value of all non-transported prognostics here
  DO js = 1, nnaf
    ! get index of nadvt(js) in NTP array
    l = name2ntpindex(nadvt(js))
    ! set data in appropriate entry in all_ntp
    ntp_data(kcs:kce,l) = zfnatr(:,js)*c_na_species(js)
  END DO

  ! First copy the concentrations from the zftr array to the
  ! diag arrays, reshape and convert to moles.
  ! The values in the troposphere/stratosphere are masked
  ! off below

  ! CH4 burden in moles. Copy tropospheric values to stratospheric array
  ! for later masking using tropospheric mask.
  trop_ch4_mol(kcs:kce) = zftr(:,n_ch4)*volume(kcs:kce)*1.0e6/avogadro
  strat_ch4_mol(kcs:kce) = trop_ch4_mol(kcs:kce)

  ! O3 burden in moles
  trop_o3_mol(kcs:kce) = zftr(:,n_o3)*volume(kcs:kce)*1.0e6/avogadro

  ! OH burden in moles
  trop_oh_mol(kcs:kce) = BE_y(:,nn_oh)*volume(kcs:kce)*1.0e6/avogadro

  IF (ukca_config%l_ukca_aerchem .OR. ukca_config%l_ukca_raqaero) THEN
    delSO2_wet_H2O2(kcs:kce) = SO2_wetox_H2O2(:)
    delSO2_wet_O3(kcs:kce) = SO2_wetox_O3(:)
    delh2so4_chem(kcs:kce) = SO2_dryox_OH(:)
    delSO2_drydep(kcs:kce) = BE_dryrt(:,nn_so2)*BE_y(:,nn_so2)*dts
    delSO2_wetdep(kcs:kce) = BE_wetrt(:,nn_so2)*BE_y(:,nn_so2)*dts
  END IF
END DO ! level loop (k)
!$OMP END DO

IF (ALLOCATED(zfnatr)) DEALLOCATE(zfnatr)
IF (ALLOCATED(ystore)) DEALLOCATE(ystore)

!$OMP END PARALLEL

! Map ntp_data back into 3D array
DO l = 1, dim_ntp
  IF (all_ntp(l)%l_required) THEN
    all_ntp(l)%data_3d(:,:,:) = RESHAPE(ntp_data(:,l),                         &
                                        [row_length,rows,model_levels])
  END IF
END DO

! Now mask off stratospheric and tropospheric diagnostics
WHERE (L_stratosphere(:))
  trop_ch4_mol(:) = 0.0
  trop_o3_mol(:) = 0.0
  trop_oh_mol(:) = 0.0
ELSE WHERE
  strat_ch4_mol(:) = 0.0
END WHERE

! Call raq diagnostic routine
IF (ukca_config%l_enable_diag_um .AND. ukca_config%l_ukca_raq) THEN
  CALL ukca_raq_diags (row_length, rows, model_levels, ntracers,               &
       dry_dep_3d, wet_dep_3d, tracer, pres, temp,                             &
       len_stashwork, stashwork)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chemistry_ctl_tropraq

END MODULE ukca_chemistry_ctl_tropraq_mod
