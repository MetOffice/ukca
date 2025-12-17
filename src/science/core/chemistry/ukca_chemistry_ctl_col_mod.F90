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
!  Main driver routine for chemistry using column-wise mode.
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
MODULE ukca_chemistry_ctl_col_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName=                            &
                                     'UKCA_CHEMISTRY_CTL_COL_MOD'

CONTAINS

SUBROUTINE ukca_chemistry_ctl_col(                                             &
                row_length, rows, model_levels, theta_field_size, ntracers,    &
                istore_h2so4,                                                  &
                pres, temp, q,                                                 &
                qcf, qcl,                                                      &
                tracer,                                                        &
                all_ntp,                                                       &
                cloud_frac,                                                    &
                photol_rates,                                                  &
                shno3_3d,                                                      &
                volume,                                                        &
                have_nat3d,                                                    &
                uph2so4inaer,                                                  &
                delso2_wet_h2o2,                                               &
                delso2_wet_o3,                                                 &
                delh2so4_chem,                                                 &
                so4_sa,                                                        &
                atm_ch4_mol,                                                   &
                atm_co_mol,                                                    &
                atm_n2o_mol,                                                   &
                atm_cf2cl2_mol,                                                &
                atm_cfcl3_mol,                                                 &
                atm_mebr_mol,                                                  &
                atm_h2_mol,                                                    &
                H_plus_3d_arr,                                                 &
                zdryrt, zwetrt, nlev_with_ddep                                 &
                )

USE asad_mod,             ONLY: advt, cdt, ctype,                              &
                                dpd, dpw, fpsc1, fpsc2,                        &
                                ihso3_h2o2, ihso3_o3, ih2so4_hv, iso2_oh,      &
                                iso3_o3, jpctr, jpcspf, jpdd, jpdw, jpnr,      &
                                jppj, jpro2, jpspec, nadvt, nlnaro2, nprkx,    &
                                o1d_in_ss, o3p_in_ss, prk, rk,                 &
                                specf, speci, sph2o, sphno3, spro2, tnd, y, za
USE asad_chem_flux_diags, ONLY: l_asad_use_chem_diags,                         &
                                l_asad_use_drydep,                             &
                                l_asad_use_flux_rxns,                          &
                                l_asad_use_psc_diagnostic,                     &
                                l_asad_use_rxn_rates,                          &
                                l_asad_use_wetdep,                             &
                                asad_psc_diagnostic,                           &
                                asad_chemical_diagnostics
USE ukca_cspecies,        ONLY: c_species, c_na_species, n_cf2cl2, n_cfcl3,    &
                                n_ch4, n_co, n_h2, n_h2so4, n_mebr, n_n2o,     &
                                nn_h2o2, nn_h2so4, nn_o1d, nn_o3, nn_o3p,      &
                                nn_oh, nn_so2
USE ukca_tropopause,      ONLY: L_stratosphere
USE ukca_constants,       ONLY: c_h2o, c_hono2, c_o1d, c_o3p, c_co2

USE ukca_config_constants_mod, ONLY: avogadro
USE ukca_config_specification_mod, ONLY: ukca_config

USE ukca_ntp_mod,       ONLY: ntp_type, dim_ntp, name2ntpindex
USE ukca_environment_fields_mod, ONLY: co2_interactive

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint

USE ukca_missing_data_mod, ONLY: rmdi

USE errormessagelength_mod, ONLY: errormessagelength

USE asad_cdrive_mod, ONLY: asad_cdrive

!!!! Note: LFRIC-specific pre-processor directives used in this module are
!!!! inappropriate in UKCA and should be removed but must be retained while
!!!! LFRic uses the UM version of ukca_um_legacy_mod which does not contain
!!!! dummy versions of autotune routines.

#if !defined(LFRIC)
USE ukca_um_legacy_mod, ONLY:                                                  &
    l_autotune_segments,                                                       &
    autotune_type,                                                             &
    autotune_init,                                                             &
    autotune_entry,                                                            &
    autotune_return,                                                           &
    autotune_start_region,                                                     &
    autotune_stop_region
#endif

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length        ! size of UKCA x dimension
INTEGER, INTENT(IN) :: rows              ! size of UKCA y dimension
INTEGER, INTENT(IN) :: model_levels      ! size of UKCA z dimension
INTEGER, INTENT(IN) :: theta_field_size  ! no. of points in horizontal
INTEGER, INTENT(IN) :: ntracers          ! no. of tracers
INTEGER, INTENT(IN) :: uph2so4inaer      ! flag for H2SO4 updating
INTEGER, INTENT(IN) :: istore_h2so4      ! location of H2SO4 in f array
INTEGER, INTENT(IN) :: nlev_with_ddep(row_length,rows) ! No levs in bl

REAL, INTENT(IN) :: pres(row_length,rows,model_levels) ! pressure
REAL, INTENT(IN) :: temp(row_length,rows,model_levels) ! actual temp
REAL, INTENT(IN) :: volume(row_length,rows,model_levels) ! cell vol.
REAL, INTENT(IN) :: H_plus_3d_arr(row_length,rows,model_levels) ! 3D pH array
REAL, INTENT(IN) :: qcf(row_length,rows,model_levels)  ! qcf
REAL, INTENT(IN) :: qcl(row_length,rows,model_levels)  ! qcl
REAL, INTENT(IN) :: cloud_frac(row_length,rows,model_levels)
REAL, INTENT(IN) :: so4_sa(row_length,rows,model_levels) ! aerosol surface area
REAL, INTENT(IN) :: zdryrt(row_length,rows,jpdd)              ! dry dep rate
REAL, INTENT(IN) :: zwetrt(row_length,rows,model_levels,jpdw) ! wet dep rate
REAL, INTENT(IN) :: photol_rates(row_length,rows,model_levels,jppj)
REAL, INTENT(OUT)    :: shno3_3d(row_length,rows,model_levels)
REAL, INTENT(IN OUT) :: q(row_length,rows,model_levels)  ! water vapour
REAL, INTENT(IN OUT) :: tracer(row_length,rows,model_levels,                   &
                               ntracers)                 ! tracer MMR

! SO2 increments
REAL, INTENT(IN OUT) :: delSO2_wet_H2O2(row_length,rows,model_levels)
REAL, INTENT(IN OUT) :: delSO2_wet_O3(row_length,rows,model_levels)
REAL, INTENT(IN OUT) :: delh2so4_chem(row_length,rows,model_levels)

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

! Non transported prognostics
TYPE(ntp_type), INTENT(IN OUT) :: all_ntp(dim_ntp)

! Mask to limit formation of Nat below specified height
LOGICAL, INTENT(IN) :: have_nat3d(row_length,rows,model_levels)

! Local variables
INTEGER :: i             ! Loop variable
INTEGER :: j             ! loop variable
INTEGER :: js            ! loop variable
INTEGER :: jtr           ! loop variable - transported tracers
INTEGER :: jro2          ! loop variable - NTP RO2 species
INTEGER :: jna           ! loop variable, non-advected species
INTEGER :: jspf          ! loop variable - all active chemical species in f
INTEGER :: k             ! loop variable
INTEGER :: klevel        ! dummy variable
INTEGER :: l             ! loop variable

INTEGER :: kcs           ! loop variable, start level of current segment/chunk
INTEGER :: kce           ! loop variable, end level of current segment/chunk
INTEGER :: chunk_size    ! Number of points within the current segment/chunk

INTEGER           :: ierr                     ! Error code: asad diags routines
INTEGER           :: errcode                  ! Error code: ereport
CHARACTER(LEN=errormessagelength) :: cmessage         ! Error message
CHARACTER(LEN=10)      :: prods(2)                 ! Products
CHARACTER(LEN=10)      :: prods3(3)                ! Products

! array to store H2SO4 when updated in MODE
REAL, ALLOCATABLE :: ystore(:)
REAL :: zftr(model_levels,jpcspf)  ! 1-D array of tracers, including RO2
REAL :: zp  (model_levels)        ! 1-D pressure
REAL :: zt  (model_levels)        ! 1-D temperature
REAL :: zclw(model_levels)        ! 1-D cloud liquid water
REAL :: zfcloud(model_levels)     ! 1-D cloud fraction
REAL :: zq(model_levels)          ! 1-D water vapour vmr
REAL :: co2_1d(model_levels)      ! 1-D CO2 vmr
REAL :: zprt1d(model_levels,jppj) ! 1-D photolysis rates for ASAD
REAL :: zdryrt2(model_levels, jpdd)               ! dry dep rate
REAL :: zwetrt2(model_levels, jpdw)               ! wet dep rat
REAL :: rc_het(model_levels,2)                ! heterog rates for trop chem
REAL :: H_plus_2d_arr(theta_field_size,model_levels) ! 2-D pH array to use in
                                                     ! wet deposition
REAL :: H_plus_1d_arr(model_levels)  ! 1-D pH array to use in chemical solvers

! Local arrays to store full column values of selected asad_mod variables
! to pass to chemical_diagnostics and psc_diagnotics subroutines
REAL :: dpd_full(model_levels,jpspec)
REAL :: dpw_full(model_levels,jpspec)
REAL :: fpsc1_full(model_levels)
REAL :: fpsc2_full(model_levels)
REAL :: prk_full(model_levels,jpnr)
REAL :: y_full(model_levels,jpspec)

LOGICAL :: l_autotune_local
LOGICAL :: stratflag(model_levels)
LOGICAL :: have_nat1d(model_levels)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMISTRY_CTL_COL'

#if !defined(LFRIC)
TYPE(autotune_type), ALLOCATABLE, SAVE :: autotune_state
#endif

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if !defined(LFRIC)
! Set up automatic segment size tuning
IF (l_autotune_segments) THEN
  IF (.NOT. ALLOCATED(autotune_state)) THEN
    ALLOCATE(autotune_state)
    CALL autotune_init(                                                        &
      autotune_state,                                                          &
      region_name    = 'ukca_chem_seg_size',                                   &
      tag            = 'UKCA-CHEM',                                            &
      start_size     = ukca_config%ukca_chem_seg_size)
  END IF

  CALL autotune_entry(autotune_state, ukca_config%ukca_chem_seg_size)
END IF
#endif

! Dummy variables to satisfy expected numbers of arguments for ASAD_CDRIVE,
! ASAD_CHEMICAL_DIAGNOSTICS and ASAD_PSC_DIAGNOSTIC
klevel=0

! if heterogeneous chemistry is selected, allocate solid HNO3 array
IF (ukca_config%l_ukca_het_psc) THEN
  shno3_3d = 0.0
END IF

#if defined(LFRIC)
! No autotuning in LFRic
l_autotune_local = .FALSE.
#else
l_autotune_local = l_autotune_segments
IF (l_autotune_segments) THEN
  CALL autotune_start_region(autotune_state, row_length*rows)
END IF
#endif

! Model levels loop
!$OMP PARALLEL DEFAULT(NONE)                                                   &
!$OMP PRIVATE(cmessage, errcode, i, ierr, j, js, l, rc_het, stratflag,         &
!$OMP         ystore, zclw, zdryrt2, zfcloud, zftr, have_nat1d,                &
!$OMP         zp, zprt1d, zq, zt, co2_1d, zwetrt2,                             &
!$OMP         kcs, kce, chunk_size, dpd_full, dpw_full,                        &
!$OMP         fpsc1_full, fpsc2_full, prk_full, y_full, jspf, jna,             &
!$OMP         H_plus_1d_arr)                                                   &
!$OMP SHARED(advt, all_ntp, atm_cf2cl2_mol, atm_cfcl3_mol, atm_ch4_mol,        &
!$OMP        atm_co_mol, atm_h2_mol, atm_mebr_mol, atm_n2o_mol, avogadro,      &
!$OMP        c_species, c_na_species, cloud_frac,                              &
!$OMP        delh2so4_chem, delSO2_wet_H2O2, delSO2_wet_O3, have_nat3d,        &
!$OMP        ih2so4_hv, ihso3_h2o2, ihso3_o3, iso2_oh, iso3_o3,                &
!$OMP        jpctr, jpdd, jpspec, jpcspf, jpro2, klevel, speci, specf,         &
!$OMP        spro2, nlnaro2, ctype, co2_interactive, istore_h2so4,             &
!$OMP        l_asad_use_chem_diags, l_asad_use_drydep,                         &
!$OMP        l_asad_use_flux_rxns, l_asad_use_psc_diagnostic,                  &
!$OMP        l_asad_use_rxn_rates, l_asad_use_wetdep, l_stratosphere,          &
!$OMP        ukca_config,                                                      &
!$OMP        l_autotune_local,                                                 &
!$OMP        model_levels, n_cf2cl2, n_cfcl3, n_ch4, n_co, n_h2,               &
!$OMP        n_mebr, n_n2o, nadvt, nlev_with_ddep,                             &
!$OMP        nn_h2o2, nn_h2so4, nn_o1d, nn_o3, nn_o3p, nn_oh, nn_so2,          &
!$OMP        o1d_in_ss, o3p_in_ss, photol_rates, pres, q, qcf, qcl,            &
!$OMP        row_length, rows, so4_sa, shno3_3d,                               &
!$OMP        temp, tracer, uph2so4inaer, volume, zdryrt, zwetrt,               &
!$OMP        H_plus_3d_arr)

IF (.NOT. ALLOCATED(ystore) .AND. uph2so4inaer == 1)                           &
                             ALLOCATE(ystore(model_levels))

! With autotuning, the segment size will - in general - have changed since
! the previous call. Since this routine allocates THREADPRIVATE arrays,
! we need to reallocate inside the parallel region.
IF (l_autotune_local) THEN
  CALL ukca_reallocate_asad_arrays(ukca_config%ukca_chem_seg_size)
END IF

!$OMP DO SCHEDULE(STATIC)
DO i=1,rows
  DO j=1,row_length

    zdryrt2 = 0.0
    IF (ukca_config%l_ukca_intdd) THEN
      ! Interactive scheme extracts from levels in boundary layer
      DO l=1,jpdd
        zdryrt2(1:nlev_with_ddep(j,i),l) = zdryrt(j,i,l)
      END DO
    ELSE    ! non-interactive
      zdryrt2(1,:) = zdryrt(j,i,:)
    END IF

    !       Put pressure, temperature and tracer mmr into 1-D arrays
    !       for use in ASAD chemical solver

    zp(:) = pres(j,i,:)
    zt(:) = temp(j,i,:)
    zq(:) = q(j,i,:)/c_h2o

    ! collate 1-D H_plus values to use in asad NR solver
    H_plus_1d_arr(:) = H_plus_3d_arr(j,i,:)

    IF (ANY(speci(:) == 'CO2       ')) THEN
      ! Copy the CO2 concentration into the asad module as VMR
      IF (ukca_config%l_chem_environ_co2_fld) THEN
        co2_1d(:) = co2_interactive(j,i,:)/c_co2
      ELSE
        co2_1d(:) = rmdi
      END IF
    END IF

    ! Put cloud liquid water and cloud fraction into 1-D arrays
    ! for use in ASAD chemical solver
    zclw(:) = qcl(j,i,:)
    zfcloud(:) = cloud_frac(j,i,:)

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
          zftr(:,jspf) = tracer(j,i,:,jtr)/c_species(jtr)
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
              zftr(:,jspf) = all_ntp(l)%data_3d(j,i,:) / c_na_species(jna)
            ELSE
              errcode = jro2
              WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl'
              CALL umPrint(umMessage,src=RoutineName)
              cmessage='ERROR: Indices for RO2 species do not match w/ nadvt'
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
        WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl'
        CALL umPrint(umMessage,src=RoutineName)
        cmessage = 'ERROR: Number of chemical active species /= jpro2+jpctr'
        CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
      END IF
    ELSE
      IF (jspf /= jpctr) THEN
        errcode = jspf
        WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl'
        CALL umPrint(umMessage,src=RoutineName)
        cmessage = 'ERROR: Number of chemical active species /= jpctr'
        CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
      END IF
    END IF


    ! Map photolysis rates onto 1-D array.
    IF (ukca_config%l_ukca_offline) THEN
      ! Offline chemistry has no photolysis
      zprt1d(:,:) = 0.0
    ELSE
      ! take column+species info from 4D photolysis array
      zprt1d(:,:) = photol_rates(j,i,:,:)
    END IF

    !       Call ASAD routines to do chemistry integration

    IF (.NOT. (ukca_config%l_ukca_trop .OR. ukca_config%l_ukca_aerchem .OR.    &
               ukca_config%l_ukca_raq .OR. ukca_config%l_ukca_raqaero)) THEN
                                                                      ! Not B-E

      ! retrieve tropospheric heterogeneous rates from previous time step
      ! for this model level (index k)
      IF (ukca_config%l_ukca_trophet) THEN
        ! N2O5
        l = name2ntpindex('het_n2o5  ')
        rc_het(:,1) = all_ntp(l)%data_3d(j,i,:)
        ! HO2+HO2
        l = name2ntpindex('het_ho2   ')
        rc_het(:,2) = all_ntp(l)%data_3d(j,i,:)
      ELSE
        rc_het(:,:) = 0.0
      END IF

      have_nat1d(:) = have_nat3d(j,i,:)

      ! fill stratospheric flag indicator
      stratflag(:) = L_stratosphere(j,i,:)

      ! pass 2D wet-dep field to cdrive
      zwetrt2(:,:) = zwetrt(j,i,:,:)

      ! Strided loop to segment the column before passing to asad_cdrive
      DO kcs = 1, model_levels, ukca_config%ukca_chem_seg_size

        kce = MIN(kcs+(ukca_config%ukca_chem_seg_size-1),model_levels)

        ! Get current chunk_size (not necessarily equal to ukca_chem_seg_size)
        chunk_size  = (kce+1)-kcs

        ! For unequal chunk sizes, reallocate arrays
        ! to ensure array conformance
        IF (chunk_size /= ukca_config%ukca_chem_seg_size) THEN
          CALL ukca_reallocate_asad_arrays(chunk_size)
        END IF

        IF (uph2so4inaer == 1) THEN
          ! H2SO4 will be updated in MODE, so store old value here
          IF (ukca_config%l_fix_ukca_h2so4_ystore) THEN
            ! primary array passed is zftr, so save this, NOT y
            ystore(kcs:kce) = zftr(kcs:kce,istore_h2so4)
          ELSE
            ystore(kcs:kce) = y(1:chunk_size,nn_h2so4)
          END IF
        END IF

        ! Initialise za for current chunk, following reallocation
        za(:)=0.0
        za(1:chunk_size) = so4_sa(j,i,kcs:kce)

        ! Initialise sph2o for current chunk
        sph2o(:) = 0.0
        IF (ukca_config%l_ukca_het_psc) THEN
          sph2o(1:chunk_size) = qcf(j,i,kcs:kce)/c_h2o
        END IF

        ! Call asad_cdrive with segmented arrays
        CALL asad_cdrive(zftr(kcs:kce,:),                                      &
                         zp(kcs:kce),                                          &
                         zt(kcs:kce),                                          &
                         zq(kcs:kce),                                          &
                         co2_1d(kcs:kce),                                      &
                         zfcloud(kcs:kce),                                     &
                         zclw(kcs:kce),                                        &
                         j,i,klevel,                                           &
                         zdryrt2(kcs:kce,:),                                   &
                         zwetrt2(kcs:kce,:),                                   &
                         rc_het(kcs:kce,:),                                    &
                         zprt1d(kcs:kce,:),                                    &
                         chunk_size,                                           &
                         have_nat1d(kcs:kce),                                  &
                         stratflag(kcs:kce),                                   &
                         H_plus_1d_arr(kcs:kce))

        ! Store the full column values of dpd, dpw, fpsc1,
        ! fpsc2, prk and y - these are needed later on
        ! for the calculation of 3D flux diagnostics
        ! outside the chunking loop
        dpd_full(kcs:kce,:)=dpd(1:chunk_size,:)
        dpw_full(kcs:kce,:)=dpw(1:chunk_size,:)
        fpsc1_full(kcs:kce)=fpsc1(1:chunk_size)
        fpsc2_full(kcs:kce)=fpsc2(1:chunk_size)
        prk_full(kcs:kce,:)=prk(1:chunk_size,:)
        y_full(kcs:kce,:)=y(1:chunk_size,:)

        IF (ukca_config%l_ukca_het_psc) THEN
          ! Save MMR of NAT PSC particles into 3-D array for PSC sedimentation.
          ! Note that sphno3 is NAT in number density of HNO3.
          IF (ANY(sphno3(:) > 0.0)) THEN
            shno3_3d(j,i,kcs:kce) = (sphno3(:)/tnd(:))*c_hono2
          ELSE
            shno3_3d(j,i,kcs:kce) = 0.0
          END IF
        END IF

        IF (ukca_config%l_ukca_chem .AND. ukca_config%l_ukca_nr_aqchem) THEN
          ! Calculate chemical fluxes for MODE
          IF (ihso3_h2o2 > 0) delSO2_wet_H2O2(j,i,kcs:kce) =                   &
            delSO2_wet_H2O2(j,i,kcs:kce) + (rk(:,ihso3_h2o2)*                  &
            y(:,nn_so2)*y(:,nn_h2o2))*cdt
          IF (ihso3_o3 > 0) delSO2_wet_O3(j,i,kcs:kce) =                       &
            delSO2_wet_O3(j,i,kcs:kce) + (rk(:,ihso3_o3)*                      &
            y(:,nn_so2)*y(:,nn_o3))*cdt
          IF (iso3_o3 > 0) delSO2_wet_O3(j,i,kcs:kce) =                        &
            delSO2_wet_O3(j,i,kcs:kce) + (rk(:,iso3_o3)*                       &
            y(:,nn_so2)*y(:,nn_o3))*cdt
          ! net H2SO4 production - note that this is affected by
          ! l_fix_ukca_h2so4_ystore above. Y value is concentration
          ! from chemistry prior to zftr being over-written below
          IF (iso2_oh > 0 .AND. ih2so4_hv > 0) THEN
            delh2so4_chem(j,i,kcs:kce) = delh2so4_chem(j,i,kcs:kce) +          &
             ((rk(:,iso2_oh)*y(:,nn_so2)*y(:,nn_oh)) -                         &
              (rk(:,ih2so4_hv)*y(:,nn_h2so4)))*cdt
          ELSE IF (iso2_oh > 0) THEN
            delh2so4_chem(j,i,kcs:kce) = delh2so4_chem(j,i,kcs:kce) +          &
              (rk(:,iso2_oh)*y(:,nn_so2)*y(:,nn_oh))*cdt
          END IF

          IF (uph2so4inaer == 1) THEN
            ! Restore H2SO4 tracer as it will be updated in MODE
            ! using delh2so4_chem
            IF (ukca_config%l_fix_ukca_h2so4_ystore) THEN
              ! calculate delh2so4_chem as the difference in H2SO4 over
              ! chemistry
              ! zftr is already in VMR, so divide by CDT to give as vmr/s
              delh2so4_chem(j,i,kcs:kce) = (zftr(kcs:kce,istore_h2so4)         &
                                             - ystore(kcs:kce)) / cdt
              ! primary array passed is zftr, so copy back to this, NOT y
              zftr(kcs:kce,istore_h2so4) = ystore(kcs:kce)
            ELSE
              y(:,nn_h2so4) = ystore(kcs:kce)
            END IF
          END IF
        END IF

        ! Bring results back from vmr to mmr.
        ! Also bring back results for nontransported RO2 to all_ntp
        jspf = 0 ! reset counter for all chemical species in f array
        DO js = 1,jpspec
          DO jtr=1,jpctr
            IF (advt(jtr) == speci(js)) THEN

              jspf = jspf+1
              tracer(j,i,kcs:kce,jtr) = zftr(kcs:kce,jspf) * c_species(jtr)

            END IF
          END DO

          ! If RO2 species are not being transported, map RO2 species
          ! back to the all_ntp 3D array
          IF (ukca_config%l_ukca_ro2_ntp) THEN
            DO jro2 = 1, jpro2
              IF (spro2(jro2) == speci(js) .AND. ctype(js) == 'OO') THEN
                jspf = jspf+1
                ! Get index of speci(js) in NTP array.
                ! If not found this is a fatal error.
                l = name2ntpindex(spro2(jro2))

                ! Find location of species in nadvt array
                jna = nlnaro2(jro2)
                all_ntp(l)%data_3d(j,i,kcs:kce) = zftr(kcs:kce,jspf) *         &
                                                  c_na_species(jna)

              END IF ! Close IF RO2 species
            END DO   ! Close loop through RO2 species
          END IF     ! Close IF RO2_NTP
        END DO       ! Close loop through all species

        ! Set SS species concentrations for output (stratospheric
        ! configurations)

        ! O1D mmr
        IF (o1d_in_ss) THEN
          l = name2ntpindex('O(1D)     ')
          all_ntp(l)%data_3d(j,i,kcs:kce) = ( y(:,nn_o1d)/tnd(:) ) * c_o1d
        END IF

        ! O3P mmr
        IF (o3p_in_ss) THEN
          l = name2ntpindex('O(3P)     ')
          all_ntp(l)%data_3d(j,i,kcs:kce) = ( y(:,nn_o3p)/tnd(:) ) * c_o3p
        END IF

        ! First copy the concentrations from the zftr array to the
        ! diag arrays, reshape and convert to moles.
        ! The indices (n_ch4, n_n2o etc.) don't refer to the location
        ! in the zftr array if RO2_NTP is true, but instead to the
        ! location in the tracer array.
        DO jspf = 1, jpcspf
          IF (n_ch4 > 0) THEN
            IF (specf(jspf) == advt(n_ch4)) THEN
              atm_ch4_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf)*tnd(:)*            &
                                  volume(j,i,kcs:kce)*1.0e6/avogadro
            END IF
          END IF

          ! CO
          IF (n_co > 0) THEN
            IF (specf(jspf) == advt(n_co)) THEN
              atm_co_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf)*tnd(:)*             &
                                  volume(j,i,kcs:kce)*1.0e6/avogadro
            END IF
          END IF

          ! N2O
          IF (n_n2o > 0) THEN
            IF (specf(jspf) == advt(n_n2o)) THEN
              atm_n2o_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf)*tnd(:)*            &
                                   volume(j,i,kcs:kce)*1.0e6/avogadro
            END IF
          END IF

          ! CFC-12
          IF (n_cf2cl2 > 0) THEN
            IF (specf(jspf) == advt(n_cf2cl2)) THEN
              atm_cf2cl2_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf)*tnd(:)*         &
                                   volume(j,i,kcs:kce)* 1.0e6/avogadro
            END IF
          END IF

          ! CFC-11
          IF (n_cfcl3 > 0) THEN
            IF (specf(jspf) == advt(n_cfcl3)) THEN
              atm_cfcl3_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf)*tnd(:)*          &
                                   volume(j,i,kcs:kce)* 1.0e6/avogadro
            END IF
          END IF

          ! CH3Br
          IF (n_mebr > 0) THEN
            IF (specf(jspf) == advt(n_mebr)) THEN
              atm_mebr_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf)*tnd(:)*           &
                                   volume(j,i,kcs:kce)* 1.0e6/avogadro
            END IF
          END IF

          ! H2
          IF (n_h2 > 0) THEN
            IF (specf(jspf) == advt(n_h2)) THEN
              atm_h2_mol(j,i,kcs:kce) = zftr(kcs:kce,jspf)*tnd(:)*             &
                                   volume(j,i,kcs:kce)* 1.0e6/avogadro
            END IF
          END IF

        END DO ! End loop through species in zftr array

      END DO ! end chunking loop

      ! If current chunk size is not equal to ukca_chem_seg_size,
      !reallocate asad arrays ready for next column
      IF (chunk_size /= ukca_config%ukca_chem_seg_size) THEN
        CALL ukca_reallocate_asad_arrays(ukca_config%ukca_chem_seg_size)
      END IF

      ! 3D flux diagnostics
      IF (L_asad_use_chem_diags .AND.                                          &
         ((L_asad_use_flux_rxns .OR. L_asad_use_rxn_rates) .OR.                &
         (L_asad_use_wetdep .OR. L_asad_use_drydep)))                          &
         CALL asad_chemical_diagnostics(row_length,rows,model_levels,          &
            model_levels,dpd_full,dpw_full,prk_full,y_full,                    &
            j,i,klevel,volume,ierr)

      ! PSC diagnostics
      IF (L_asad_use_chem_diags .AND. L_asad_use_psc_diagnostic)               &
         CALL asad_psc_diagnostic(row_length,rows,model_levels,model_levels,   &
                             fpsc1_full,fpsc2_full,                            &
                             j,i,klevel,ierr)
    ELSE
      cmessage='Column call is not available for Backward Euler schemes'
      errcode = 5
      CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
    END IF
  END DO
END DO ! loop (j,i)
!$OMP END DO

IF (ALLOCATED(ystore)) DEALLOCATE(ystore)

!$OMP END PARALLEL

#if !defined(LFRIC)
IF (l_autotune_segments) THEN
  CALL autotune_stop_region(autotune_state)
END IF
#endif

#if !defined(LFRIC)
! If autotuning is active, decide what to do with the
! trial segment size and report the current state.
IF (l_autotune_segments) THEN
  CALL autotune_return(autotune_state)
END IF
#endif

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chemistry_ctl_col

SUBROUTINE ukca_reallocate_asad_arrays(n_pnts)

USE asad_mod, ONLY: prod, slos, pd, co3, deriv, dpd, dpw,                      &
                     ej, emr, f, fdot, fj, fpsc1, fpsc2,                       &
                     ftilde, ipa, lati, linfam, method, p,                     &
                     pmintnd, prk, qa, ratio, rk, sh2o, shno3,                 &
                     sph2o, sphno3, t, t300, tnd, wp, co2, y,                  &
                     ydot, za, spfj, spfjsize_max,                             &
                     jpspec, jpcspf, jpnr
USE ukca_config_specification_mod, ONLY: ukca_config, int_method_nr

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_pnts

! Set integration method (1 = IMPACT; 3 = N-R solver; 5 = Backward-Euler)
method = ukca_config%ukca_int_method

! Deallocate asad mod variables (in reverse order
! to which they are initially allocated in ukca_mod.F90)...

IF (method == int_method_NR) THEN ! sparse_vars
  IF (ALLOCATED(spfj)) DEALLOCATE(spfj)
END IF

IF (ALLOCATED(za)) DEALLOCATE(za)
IF (ALLOCATED(ydot)) DEALLOCATE(ydot)
IF (ALLOCATED(y)) DEALLOCATE(y)
IF (ALLOCATED(co2)) DEALLOCATE(co2)
IF (ALLOCATED(wp)) DEALLOCATE(wp)
IF (ALLOCATED(tnd)) DEALLOCATE(tnd)
IF (ALLOCATED(t300)) DEALLOCATE(t300)
IF (ALLOCATED(t)) DEALLOCATE(t)
IF (ALLOCATED(sphno3)) DEALLOCATE(sphno3)
IF (ALLOCATED(sph2o)) DEALLOCATE(sph2o)
IF (ALLOCATED(shno3)) DEALLOCATE(shno3)
IF (ALLOCATED(sh2o)) DEALLOCATE(sh2o)
IF (ALLOCATED(rk)) DEALLOCATE(rk)
IF (ALLOCATED(ratio)) DEALLOCATE(ratio)
IF (ALLOCATED(qa)) DEALLOCATE(qa)
IF (ALLOCATED(prk)) DEALLOCATE(prk)
IF (ALLOCATED(pmintnd)) DEALLOCATE(pmintnd)
IF (ALLOCATED(p)) DEALLOCATE(p)
IF (ALLOCATED(linfam)) DEALLOCATE(linfam)
IF (ALLOCATED(lati)) DEALLOCATE(lati)
IF (ALLOCATED(ipa)) DEALLOCATE(ipa)
IF (ALLOCATED(ftilde)) DEALLOCATE(ftilde)
IF (ALLOCATED(fpsc2)) DEALLOCATE(fpsc2)
IF (ALLOCATED(fpsc1))  DEALLOCATE(fpsc1)
IF (ALLOCATED(fj)) DEALLOCATE(fj)
IF (ALLOCATED(fdot)) DEALLOCATE(fdot)
IF (ALLOCATED(f)) DEALLOCATE(f)
IF (ALLOCATED(emr)) DEALLOCATE(emr)
IF (ALLOCATED(ej)) DEALLOCATE(ej)
IF (ALLOCATED(dpw)) DEALLOCATE(dpw)
IF (ALLOCATED(dpd)) DEALLOCATE(dpd)
IF (ALLOCATED(deriv)) DEALLOCATE(deriv)
IF (ALLOCATED(co3)) DEALLOCATE(co3)
IF (ALLOCATED(pd)) DEALLOCATE(pd)

!...and re-allocate based on value of n_pnts
IF (.NOT. ALLOCATED(pd)) ALLOCATE(pd(n_pnts,2*jpspec))
IF (.NOT. ALLOCATED(co3)) ALLOCATE(co3(n_pnts))
IF (.NOT. ALLOCATED(deriv)) ALLOCATE(deriv(n_pnts,4,4))
IF (.NOT. ALLOCATED(dpd)) ALLOCATE(dpd(n_pnts,jpspec))
IF (.NOT. ALLOCATED(dpw)) ALLOCATE(dpw(n_pnts,jpspec))
IF (.NOT. ALLOCATED(ej)) ALLOCATE(ej(n_pnts,jpcspf))
IF (.NOT. ALLOCATED(emr)) ALLOCATE(emr(n_pnts,jpspec))
IF (.NOT. ALLOCATED(f)) ALLOCATE(f(n_pnts,jpcspf))
IF (.NOT. ALLOCATED(fdot)) ALLOCATE(fdot(n_pnts,jpcspf))
IF (.NOT. ALLOCATED(fj)) ALLOCATE(fj(n_pnts,jpcspf,jpcspf))
IF (.NOT. ALLOCATED(fpsc1)) ALLOCATE(fpsc1(n_pnts))
IF (.NOT. ALLOCATED(fpsc2)) ALLOCATE(fpsc2(n_pnts))
IF (.NOT. ALLOCATED(ftilde)) ALLOCATE(ftilde(n_pnts, jpcspf))
IF (.NOT. ALLOCATED(ipa)) ALLOCATE(ipa(n_pnts,jpcspf))
IF (.NOT. ALLOCATED(lati)) ALLOCATE(lati(n_pnts))
IF (.NOT. ALLOCATED(linfam)) ALLOCATE(linfam(n_pnts,0:jpcspf))
IF (.NOT. ALLOCATED(p)) ALLOCATE(p(n_pnts))
IF (.NOT. ALLOCATED(pmintnd)) ALLOCATE(pmintnd(n_pnts))
IF (.NOT. ALLOCATED(prk)) ALLOCATE(prk(n_pnts,jpnr))
IF (.NOT. ALLOCATED(qa)) ALLOCATE(qa(n_pnts,jpspec))
IF (.NOT. ALLOCATED(ratio)) ALLOCATE(ratio(n_pnts,jpspec))
IF (.NOT. ALLOCATED(rk)) ALLOCATE(rk(n_pnts,jpnr))
IF (.NOT. ALLOCATED(sh2o)) ALLOCATE(sh2o(n_pnts))
IF (.NOT. ALLOCATED(shno3)) ALLOCATE(shno3(n_pnts))
IF (.NOT. ALLOCATED(sph2o)) ALLOCATE(sph2o(n_pnts))
IF (.NOT. ALLOCATED(sphno3)) ALLOCATE(sphno3(n_pnts))
IF (.NOT. ALLOCATED(t)) ALLOCATE(t(n_pnts))
IF (.NOT. ALLOCATED(t300)) ALLOCATE(t300(n_pnts))
IF (.NOT. ALLOCATED(tnd)) ALLOCATE(tnd(n_pnts))
IF (.NOT. ALLOCATED(wp)) ALLOCATE(wp(n_pnts))
IF (.NOT. ALLOCATED(co2)) ALLOCATE(co2(n_pnts))
IF (.NOT. ALLOCATED(y)) ALLOCATE(y(n_pnts,jpspec))
IF (.NOT. ALLOCATED(ydot)) ALLOCATE(ydot(n_pnts,jpspec))
IF (.NOT. ALLOCATED(za)) ALLOCATE(za(n_pnts))

IF (method == int_method_NR) THEN ! sparse vars
  IF (.NOT. ALLOCATED(spfj)) ALLOCATE(spfj(n_pnts,spfjsize_max))
END IF

NULLIFY(prod)
NULLIFY(slos)
prod => pd(:,1:jpspec)
slos => pd(:,jpspec+1:2*jpspec)

! (re-)initialise DERIV array to 1.0 before each call to ASAD_CDRIVE
! to ensure bit-comparability when changing domain decomposition
deriv(:,:,:) = 1.0

!     Clear the species arrays
f(:,:)      = 0.0
fdot(:,:)   = 0.0
ej(:,:)     = 0.0
linfam(:,:) = .FALSE.

y(:,:)    = 0.0
ydot(:,:) = 0.0
prod(:,:) = 0.0
slos(:,:) = 0.0
dpd(:,:)  = 0.0
dpw(:,:)  = 0.0
emr(:,:)  = 0.0

!     Clear the rates and index arrays
rk(:,:)   = 0.0
prk(:,:)  = 0.0

RETURN
END SUBROUTINE ukca_reallocate_asad_arrays

END MODULE ukca_chemistry_ctl_col_mod
