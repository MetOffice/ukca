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
!  Main driver routine for chemistry using level-wise mode.
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
                row_length, rows, model_levels, theta_field_size, tot_n_pnts,  &
                ntracers,                                                      &
                istore_h2so4,                                                  &
                pres, temp, q,                                                 &
                qcf, qcl,                                                      &
                tracer,                                                        &
                all_ntp,                                                       &
                cloud_frac,                                                    &
                photol_rates,                                                  &
                shno3_3d,                                                      &
                volume,                                                        &
                have_nat,                                                      &
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
                H_plus,                                                        &
                zdryrt, zwetrt, nlev_with_ddep, co2_interactive,               &
                L_stratosphere, firstcall                                      &
                )

USE asad_mod,             ONLY: advt, cdt, ctype,                              &
                                ihso3_h2o2, ihso3_o3, ih2so4_hv, iso2_oh,      &
                                iso3_o3, jpctr, jpcspf, jpdd, jpdw, jpnr,      &
                                jppj, jpro2, jpspec, nadvt, nlnaro2, nprkx,    &
                                o1d_in_ss, o3p_in_ss, rk,                      &
                                specf, speci, sph2o, sphno3, spro2, tnd, y, za
USE asad_chem_flux_diags, ONLY: l_asad_use_chem_diags,                         &
                                l_asad_use_drydep,                             &
                                l_asad_use_flux_rxns,                          &
                                l_asad_use_psc_diagnostic,                     &
                                l_asad_use_rxn_rates,                          &
                                l_asad_use_wetdep,                             &
                                asad_psc_diagnostic,                           &
                                asad_chemical_diagnostics
USE asad_cdrive_mod,      ONLY: asad_cdrive
USE ukca_cspecies,        ONLY: c_species, c_na_species, n_cf2cl2, n_cfcl3,    &
                                n_ch4, n_co, n_h2, n_h2so4, n_mebr, n_n2o,     &
                                nn_h2o2, nn_h2so4, nn_o1d, nn_o3, nn_o3p,      &
                                nn_oh, nn_so2
USE ukca_constants,       ONLY: c_h2o, c_hono2, c_o1d, c_o3p, c_co2
USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_ntp_mod,         ONLY: ntp_type, dim_ntp, name2ntpindex
USE chemistry_constants_mod, ONLY: avogadro
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
USE ereport_mod,          ONLY: ereport
USE umPrintMgr,           ONLY: umMessage, umPrint
USE missing_data_mod,     ONLY: rmdi
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length        ! size of UKCA x dimension
INTEGER, INTENT(IN) :: rows              ! size of UKCA y dimension
INTEGER, INTENT(IN) :: model_levels      ! size of UKCA z dimension
INTEGER, INTENT(IN) :: theta_field_size  ! no. of points in horizontal
INTEGER, INTENT(IN) :: tot_n_pnts        ! no. of points in full domain
INTEGER, INTENT(IN) :: ntracers          ! no. of tracers
INTEGER, INTENT(IN) :: uph2so4inaer      ! flag for H2SO4 updating
INTEGER, INTENT(IN) :: istore_h2so4      ! location of H2SO4 in f array
INTEGER, INTENT(IN) :: nlev_with_ddep(theta_field_size) ! No levs in bl

REAL, INTENT(IN) :: pres(tot_n_pnts)                ! pressure
REAL, INTENT(IN) :: temp(tot_n_pnts)                ! actual temperature
REAL, INTENT(IN) :: volume(tot_n_pnts)              ! cell volume
REAL, INTENT(IN) :: H_plus(tot_n_pnts)              ! pH array
REAL, INTENT(IN) :: qcf(tot_n_pnts)
REAL, INTENT(IN) :: qcl(tot_n_pnts)
REAL, INTENT(IN) :: cloud_frac(tot_n_pnts)
REAL, INTENT(IN) :: so4_sa(tot_n_pnts)              ! aerosol surface area
REAL, INTENT(IN) :: zdryrt(theta_field_size,jpdd)   ! dry dep rate
REAL, INTENT(IN) :: zwetrt(tot_n_pnts,jpdw)         ! wet dep rate
REAL, INTENT(IN) :: photol_rates(tot_n_pnts,jppj)
REAL, INTENT(IN) :: co2_interactive(tot_n_pnts)
REAL, INTENT(OUT) :: shno3_3d(tot_n_pnts)
REAL, INTENT(IN OUT) :: q(tot_n_pnts)               ! water vapour
REAL, INTENT(IN OUT) :: tracer(tot_n_pnts,ntracers) ! tracer MMR

! SO2 increments
REAL, INTENT(IN OUT) :: delSO2_wet_H2O2(tot_n_pnts)
REAL, INTENT(IN OUT) :: delSO2_wet_O3(tot_n_pnts)
REAL, INTENT(IN OUT) :: delh2so4_chem(tot_n_pnts)

! Atmospheric Burden of CH4
REAL, INTENT(IN OUT) :: atm_ch4_mol(tot_n_pnts)

! Atmospheric Burden of CO
REAL, INTENT(IN OUT) :: atm_co_mol(tot_n_pnts)

! Atmospheric Burden of Nitrous Oxide (N2O)
REAL, INTENT(IN OUT) :: atm_n2o_mol(tot_n_pnts)

! Atmospheric Burden of CFC-12
REAL, INTENT(IN OUT) :: atm_cf2cl2_mol(tot_n_pnts)

! Atmospheric Burden of CFC-11
REAL, INTENT(IN OUT) :: atm_cfcl3_mol(tot_n_pnts)

! Atmospheric Burden of CH3Br
REAL, INTENT(IN OUT) :: atm_mebr_mol(tot_n_pnts)

! Atmospheric Burden of H2
REAL, INTENT(IN OUT) :: atm_h2_mol(tot_n_pnts)

! Non transported prognostics
TYPE(ntp_type), INTENT(IN OUT) :: all_ntp(dim_ntp)

! Mask to limit formation of Nat below specified height
LOGICAL, INTENT(IN) :: have_nat(tot_n_pnts)

! Stratosphere mask
LOGICAL, INTENT(IN) :: L_stratosphere(tot_n_pnts)

! Flag for determining if this is the first chemistry call
LOGICAL, INTENT(IN) :: firstcall

! Local variables
INTEGER :: ix            ! dummy variable
INTEGER :: jy            ! dummy variable
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

INTEGER :: ierr                               ! Error code: asad diags routines
INTEGER :: errcode                            ! Error code: ereport
CHARACTER(LEN=errormessagelength) :: cmessage ! Error message
CHARACTER(LEN=10) :: prods(2)                 ! Products
CHARACTER(LEN=10) :: prods3(3)                ! Products
LOGICAL :: ddmask(theta_field_size)           ! mask

! Dummy variables to satisfy expected numbers of arguments for ASAD_CDRIVE,
! ASAD_CHEMICAL_DIAGNOSTICS and ASAD_PSC_DIAGNOSTIC
REAL :: dpd_dummy(model_levels,jpspec)
REAL :: dpw_dummy(model_levels,jpspec)
REAL :: prk_dummy(model_levels,jpnr)
REAL :: y_dummy(model_levels,jpspec)
REAL :: fpsc1_dummy(model_levels)
REAL :: fpsc2_dummy(model_levels)

REAL, ALLOCATABLE :: ystore(:)          ! array for H2SO4 when updated in MODE
REAL :: zftr(theta_field_size,jpcspf)   ! 1-D array of chemically active species
                                        !   including RO2 species, in VMR
REAL :: cdot(theta_field_size,jpcspf)   ! 1-D chem. tendency
REAL :: zq(theta_field_size)            ! 1-D water vapour vmr
REAL :: co2_1d(theta_field_size)        ! 1-D CO2
REAL :: zprt1d(theta_field_size,jppj)   ! 1-D photolysis rates for ASAD
REAL :: zdryrt2(theta_field_size, jpdd) ! dry dep rate
REAL :: rc_het(theta_field_size,2)      ! heterog rates for trop chem

! Full ntp array
REAL :: ntp_data(tot_n_pnts,dim_ntp)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMISTRY_CTL'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Dummy variables to satisfy expected numbers of arguments for ASAD_CDRIVE,
! ASAD_CHEMICAL_DIAGNOSTICS and ASAD_PSC_DIAGNOSTIC
ix=0
jy=0
dpd_dummy=0.0
dpw_dummy=0.0
fpsc1_dummy=0.0
fpsc2_dummy=0.0
prk_dummy=0.0
y_dummy=0.0

! if heterogeneous chemistry is selected, allocate solid HNO3 array
IF (ukca_config%l_ukca_het_psc) THEN
  shno3_3d = 0.0
END IF

! Fill stratospheric flag indicator and ntp_data array
DO l = 1, dim_ntp
  IF (all_ntp(l)%l_required) THEN
    ntp_data(:,l) = RESHAPE(all_ntp(l)%data_3d(:,:,:),[tot_n_pnts])
  END IF
END DO

!$OMP PARALLEL DEFAULT(NONE)                                                   &
!$OMP PRIVATE(cdot, ddmask, errcode, ierr, jna, jro2, js, jspf, jtr, k, kcs,   &
!$OMP         kce, l, rc_het, ystore, zdryrt2, zftr, zprt1d, zq, co2_1d)       &
!$OMP SHARED(advt, atm_cf2cl2_mol, atm_cfcl3_mol, atm_ch4_mol,                 &
!$OMP        atm_co_mol, atm_h2_mol, atm_mebr_mol, atm_n2o_mol,                &
!$OMP        c_species, c_na_species, cloud_frac, cmessage,                    &
!$OMP        delh2so4_chem, delSO2_wet_H2O2, delSO2_wet_O3,                    &
!$OMP        dpd_dummy, dpw_dummy, fpsc1_dummy, fpsc2_dummy, H_plus,           &
!$OMP        have_nat, prk_dummy, speci, specf, y_dummy, co2_interactive,      &
!$OMP        ih2so4_hv, ihso3_h2o2, ihso3_o3, iso2_oh, iso3_o3,                &
!$OMP        ix, jpctr, jpdd, jpdw, jppj, jpspec, jpro2, jy,                   &
!$OMP        jpcspf, spro2, nlnaro2, ctype, istore_h2so4,                      &
!$OMP        l_asad_use_chem_diags, l_asad_use_drydep,                         &
!$OMP        l_asad_use_flux_rxns, l_asad_use_psc_diagnostic,                  &
!$OMP        l_asad_use_rxn_rates, l_asad_use_wetdep, l_stratosphere,          &
!$OMP        ukca_config,                                                      &
!$OMP        model_levels,                                                     &
!$OMP        n_cf2cl2, n_cfcl3, n_ch4, n_co, n_h2,                             &
!$OMP        n_mebr, n_n2o, nadvt, nlev_with_ddep,                             &
!$OMP        nn_h2o2, nn_h2so4, nn_o1d, nn_o3, nn_o3p, nn_oh, nn_so2,          &
!$OMP        ntp_data, o1d_in_ss, o3p_in_ss, photol_rates,                     &
!$OMP        pres, q, qcf, qcl, row_length, rows,                              &
!$OMP        shno3_3d, so4_sa,                                                 &
!$OMP        temp, theta_field_size, tracer,                                   &
!$OMP        uph2so4inaer, volume,                                             &
!$OMP        zdryrt, zwetrt)

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

  zdryrt2(:,:) = 0.0e0
  IF (ukca_config%l_ukca_intdd) THEN
     ! Interactive scheme extracts from levels in boundary layer
    ddmask(:) = (k <= nlev_with_ddep(:))
    DO l=1,jpdd
      WHERE (ddmask(:))
        zdryrt2(:,l) = zdryrt(:,l)
      END WHERE
    END DO
  ELSE    ! non-interactive
    IF (k == 1) THEN
      zdryrt2(:,:) = zdryrt(:,:)
    END IF
  END IF

  ! Put tracer mmr into 1-D array for use in ASAD chemical solver
  zq(:) = q(kcs:kce)/c_h2o

  IF (ANY(speci(:) == 'CO2       ')) THEN
    !  Copy the CO2 concentration into the asad module as VMR
    IF (ukca_config%l_chem_environ_co2_fld) THEN
      co2_1d(:) = co2_interactive(kcs:kce)/c_co2
    ELSE
      co2_1d(:) = rmdi
    END IF

  END IF       ! CO2 as species

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

          ! Map data from ntp_data array into zftr
          IF (nadvt(jna) == spro2(jro2)) THEN
            jspf = jspf+1
            zftr(:,jspf) = ntp_data(kcs:kce,l) / c_na_species(jna)
          ELSE
            errcode = jro2
            WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl'
            CALL umPrint(umMessage,src='ukca_chemistry_ctl')
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
      WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl'
      CALL umPrint(umMessage,src='ukca_chemistry_ctl')
      cmessage = 'ERROR: Number of chemical active species /= jpro2+jpctr'
      CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
    END IF
  ELSE
    IF (jspf /= jpctr) THEN
      errcode = jspf
      WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl'
      CALL umPrint(umMessage,src='ukca_chemistry_ctl')
      cmessage = 'ERROR: Number of chemical active species /= jpctr'
      CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
    END IF
  END IF

  ! Map photolysis rates onto 1-D array.
  IF (ukca_config%l_ukca_offline) THEN
    ! Offline chemistry has no photolysis
    zprt1d(:,:) = 0.0
  ELSE
    zprt1d(:,:) = photol_rates(kcs:kce,:)
  END IF

  !       Call ASAD routines to do chemistry integration

  ! retrieve tropospheric heterogeneous rates from previous time step
  ! for this model level (index k)
  IF (ukca_config%l_ukca_trophet) THEN
    ! N2O5
    l = name2ntpindex('het_n2o5  ')
    rc_het(:,1) = ntp_data(kcs:kce,l)
    ! HO2+HO2
    l = name2ntpindex('het_ho2   ')
    rc_het(:,2) = ntp_data(kcs:kce,l)
  ELSE
    rc_het(:,:) = 0.0
  END IF

  za(:) = so4_sa(kcs:kce)

  IF (uph2so4inaer == 1) THEN
    ! H2SO4 will be updated in MODE, so store old value here
    IF (ukca_config%l_fix_ukca_h2so4_ystore) THEN
       ! primary array passed is zftr, so save this, NOT y
      ystore(:) = zftr(:,istore_h2so4)
    ELSE
      ystore(:) = y(:,nn_h2so4)
    END IF
  END IF

  CALL asad_cdrive(cdot,                                                       &
                   zftr,                                                       &
                   pres(kcs:kce),                                              &
                   temp(kcs:kce),                                              &
                   zq,                                                         &
                   co2_1d,                                                     &
                   cloud_frac(kcs:kce),                                        &
                   qcl(kcs:kce),                                               &
                   ix, jy, k,                                                  &
                   zdryrt2,                                                    &
                   zwetrt(kcs:kce,:),                                          &
                   rc_het,                                                     &
                   zprt1d,                                                     &
                   theta_field_size,                                           &
                   have_nat(kcs:kce),                                          &
                   L_stratosphere(kcs:kce),                                    &
                   H_plus(kcs:kce))

  IF (ukca_config%l_ukca_het_psc) THEN
    ! Save MMR of NAT PSC particles into 3-D array for PSC sedimentation.
    ! Note that sphno3 is NAT in number density of HNO3.
    IF (ANY(sphno3(1:theta_field_size) > 0.0)) THEN
      shno3_3d(kcs:kce) = sphno3(:)/tnd(:)*c_hono2
    ELSE
      shno3_3d(kcs:kce) = 0.0
    END IF
  END IF

  IF (ukca_config%l_ukca_chem .AND. ukca_config%l_ukca_nr_aqchem) THEN
    ! Calculate chemical fluxes for MODE
    IF (ihso3_h2o2 > 0) THEN
      delSO2_wet_H2O2(kcs:kce) = delSO2_wet_H2O2(kcs:kce) +                    &
        rk(:,ihso3_h2o2)*y(:,nn_so2)*y(:,nn_h2o2)*cdt
    END IF
    IF (ihso3_o3 > 0) THEN
      delSO2_wet_O3(kcs:kce) = delSO2_wet_O3(kcs:kce) +                        &
        rk(:,ihso3_o3)*y(:,nn_so2)*y(:,nn_o3)*cdt
    END IF
    IF (iso3_o3 > 0) THEN
      delSO2_wet_O3(kcs:kce) = delSO2_wet_O3(kcs:kce) +                        &
        rk(:,iso3_o3)*y(:,nn_so2)*y(:,nn_o3)*cdt
    END IF
    ! net H2SO4 production - note that this is affected by
    ! l_fix_ukca_h2so4_ystore above. Y value is concentration
    ! from chemistry prior to zftr being over-written below
    IF (iso2_oh > 0 .AND. ih2so4_hv > 0) THEN
      delh2so4_chem(kcs:kce) = delh2so4_chem(kcs:kce) +                        &
       (rk(:,iso2_oh)*y(:,nn_so2)*y(:,nn_oh) -                                 &
        rk(:,ih2so4_hv)*y(:,nn_h2so4))*cdt
    ELSE IF (iso2_oh > 0) THEN
      delh2so4_chem(kcs:kce) = delh2so4_chem(kcs:kce) +                        &
       rk(:,iso2_oh)*y(:,nn_so2)*y(:,nn_oh)*cdt
    END IF

    IF (uph2so4inaer == 1) THEN
       ! Restore H2SO4 tracer as it will be updated in MODE
       ! using delh2so4_chem
      IF (ukca_config%l_fix_ukca_h2so4_ystore) THEN
         ! calculate delh2so4_chem as the difference in H2SO4 over chemistry
         ! zftr is already in VMR, so divide by CDT to give as vmr/s
        delh2so4_chem(kcs:kce) = (zftr(:,istore_h2so4) - ystore(:)) / cdt
         ! primary array passed is zftr, so copy back to this, NOT y
        zftr(:,istore_h2so4) = ystore(:)
      ELSE
        y(:,nn_h2so4) = ystore(:)
      END IF
    END IF
  END IF

  ! 3D flux diagnostics
  IF (L_asad_use_chem_diags .AND.                                              &
       ((L_asad_use_flux_rxns .OR. L_asad_use_rxn_rates) .OR.                  &
       (L_asad_use_wetdep .OR. L_asad_use_drydep)))                            &
       CALL asad_chemical_diagnostics(row_length,rows,model_levels,            &
       theta_field_size,dpd_dummy,dpw_dummy,prk_dummy,y_dummy,                 &
       ix,jy,k,volume,ierr)

  ! PSC diagnostics
  IF (L_asad_use_chem_diags .AND. L_asad_use_psc_diagnostic)                   &
       CALL asad_psc_diagnostic(row_length,rows,model_levels,theta_field_size, &
       fpsc1_dummy,fpsc2_dummy,ix,jy,k,ierr)

  ! Bring results back from vmr to mmr.
  ! Also bring back results for nontransported RO2 to all_ntp
  jspf = 0 ! reset counter for species in f array
  DO js = 1,jpspec
    DO jtr=1,jpctr
      IF (advt(jtr) == speci(js)) THEN
        jspf = jspf+1
        tracer(kcs:kce,jtr) = zftr(:,jspf)*c_species(jtr)
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
          ntp_data(kcs:kce,l) = zftr(:,jspf) * c_na_species(jna)
        END IF ! Close IF RO2 species
      END DO   ! Close loop through RO2 species
    END IF     ! Close IF RO2_NTP
  END DO       ! Close loop through all species

  ! Set SS species concentrations for output (stratospheric configurations)

  ! O1D mmr
  IF (O1D_in_ss) THEN
    l = name2ntpindex('O(1D)     ')
    ntp_data(kcs:kce,l) = y(:,nn_o1d)/tnd(:)*c_o1d
  END IF

  ! O3P mmr
  IF (O3P_in_ss) THEN
    l = name2ntpindex('O(3P)     ')
    ntp_data(kcs:kce,l) = y(:,nn_o3p)/tnd(:)*c_o3p
  END IF

  ! First copy the concentrations from the zftr array to the
  ! diag arrays, reshape and convert to moles.
  ! The indices (n_ch4, n_n2o etc.) don't refer to the location
  ! in the zftr array if RO2_NTP is true, but instead to the
  ! location in the tracer array.
  DO jspf = 1, jpcspf

    IF (n_ch4 > 0) THEN
      IF (specf(jspf) == advt(n_ch4)) THEN
        atm_ch4_mol(kcs:kce) = zftr(:,jspf)*tnd(:)*volume(kcs:kce)*            &
                               1.0e6/avogadro
      END IF
    END IF

    ! CO
    IF (n_co > 0) THEN
      IF (specf(jspf) == advt(n_co)) THEN
        atm_co_mol(kcs:kce) = zftr(:,jspf)*tnd(:)*volume(kcs:kce)*             &
                              1.0e6/avogadro
      END IF
    END IF

    ! N2O
    IF (n_n2o > 0) THEN
      IF (specf(jspf) == advt(n_n2o)) THEN
        atm_n2o_mol(kcs:kce) = zftr(:,jspf)*tnd(:)*volume(kcs:kce)*            &
                               1.0e6/avogadro
      END IF
    END IF

    ! CFC-12
    IF (n_cf2cl2 > 0) THEN
      IF (specf(jspf) == advt(n_cf2cl2)) THEN
        atm_cf2cl2_mol(kcs:kce) = zftr(:,jspf)*tnd(:)*volume(kcs:kce)*         &
                                  1.0e6/avogadro
      END IF
    END IF

    ! CFC-11
    IF (n_cfcl3 > 0) THEN
      IF (specf(jspf) == advt(n_cfcl3)) THEN
        atm_cfcl3_mol(kcs:kce) = zftr(:,jspf)*tnd(:)*volume(kcs:kce)*          &
                                 1.0e6/avogadro
      END IF
    END IF

    ! CH3Br
    IF (n_mebr > 0) THEN
      IF (specf(jspf) == advt(n_mebr)) THEN
        atm_mebr_mol(kcs:kce) = zftr(:,jspf)*tnd(:)*volume(kcs:kce)*           &
                                1.0e6/avogadro
      END IF
    END IF

    ! H2
    IF (n_h2 > 0) THEN
      IF (specf(jspf) == advt(n_h2)) THEN
        atm_h2_mol(kcs:kce) = zftr(:,jspf)*tnd(:)*volume(kcs:kce)*             &
                              1.0e6/avogadro
      END IF
    END IF

  END DO ! End loop through species in zftr array

END DO ! level loop (k)
!$OMP END DO

IF (ALLOCATED(ystore)) DEALLOCATE(ystore)

!$OMP END PARALLEL

! Map ntp_data back into 3D array
DO l = 1, dim_ntp
  IF (all_ntp(l)%l_required) THEN
    all_ntp(l)%data_3d(:,:,:) = RESHAPE(ntp_data(:,l),                         &
                                        [row_length,rows,model_levels])
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chemistry_ctl

END MODULE ukca_chemistry_ctl_mod
