! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2024]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Description:
!  Main driver routine for chemistry using full domain mode.
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
MODULE ukca_chemistry_ctl_full_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CHEMISTRY_CTL_FULL_MOD'

CONTAINS

SUBROUTINE ukca_chemistry_ctl_full(                                            &
                row_length, rows, model_levels, theta_field_size, tot_n_pnts,  &
                chunk_n_x, chunk_n_y, chunk_n_z,                               &
                ntracers,                                                      &
                istore_h2so4,                                                  &
                pres, temp, q,                                                 &
                qcf, qcl,                                                      &
                tracer,                                                        &
                all_ntp,                                                       &
                cloud_frac,                                                    &
                photol_rates,                                                  &
                shno3,                                                         &
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
                zdryrt, zwetrt, nlev_with_ddep, L_stratosphere,                &
                co2_interactive, firstcall                                     &
                )

USE asad_mod,             ONLY: advt, cdt_diag, ctype,                         &
                                ihso3_h2o2, ihso3_o3, ih2so4_hv, iso2_oh,      &
                                iso3_o3, jpctr, jpcspf, jpdd, jpdw, jpnr,      &
                                jppj, jpro2, jpspec, nadvt, nlnaro2, nprkx,    &
                                o1d_in_ss, o3p_in_ss, rk,                      &
                                specf, speci, sph2o, sphno3, spro2, tnd, za,   &
                                dpd, dpw, prk, y, fpsc1, fpsc2
USE asad_cdrive_mod,      ONLY: asad_cdrive
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
USE ukca_constants,       ONLY: c_h2o, c_hono2, c_o1d, c_o3p, c_co2
USE ukca_config_constants_mod, ONLY: avogadro
USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_ntp_mod,         ONLY: ntp_type, dim_ntp, name2ntpindex

USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
USE ereport_mod,          ONLY: ereport
USE umPrintMgr,           ONLY: umMessage, umPrint

USE ukca_missing_data_mod, ONLY: rmdi

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length        ! size of UKCA x dimension
INTEGER, INTENT(IN) :: rows              ! size of UKCA y dimension
INTEGER, INTENT(IN) :: model_levels      ! size of UKCA z dimension
INTEGER, INTENT(IN) :: theta_field_size  ! no. of points in horizontal
INTEGER, INTENT(IN) :: tot_n_pnts        ! no. of points in full domain
INTEGER, INTENT(IN) :: chunk_n_x         ! chunk size in x dimension
INTEGER, INTENT(IN) :: chunk_n_y         ! chunk size in y dimension
INTEGER, INTENT(IN) :: chunk_n_z         ! chunk size in z dimension
INTEGER, INTENT(IN) :: ntracers          ! no. of tracers
INTEGER, INTENT(IN) :: uph2so4inaer      ! flag for H2SO4 updating
INTEGER, INTENT(IN) :: istore_h2so4      ! location of H2SO4 in f array
INTEGER, INTENT(IN) :: nlev_with_ddep(row_length,rows) ! No levs in bl

REAL, INTENT(IN) :: pres(row_length,rows,model_levels) ! pressure
REAL, INTENT(IN) :: temp(row_length,rows,model_levels) ! actual temperature
REAL, INTENT(IN) :: volume(row_length,rows,model_levels) ! cell volume
REAL, INTENT(IN) :: H_plus(row_length,rows,model_levels) ! pH array
REAL, INTENT(IN) :: qcf(row_length,rows,model_levels)
REAL, INTENT(IN) :: qcl(row_length,rows,model_levels)
REAL, INTENT(IN) :: cloud_frac(row_length,rows,model_levels)
REAL, INTENT(IN) :: so4_sa(row_length,rows,model_levels) ! aerosol surface area
REAL, INTENT(IN) :: zdryrt(row_length,rows,jpdd) ! dry dep rate
REAL, INTENT(IN) :: zwetrt(row_length,rows,model_levels,jpdw) ! wet dep rate
REAL, INTENT(IN) :: photol_rates(row_length,rows,model_levels,jppj)

! must be allocatable as passed unallocated from main if l_chem_environ_co2_fld
! is false
REAL, INTENT(IN), ALLOCATABLE :: co2_interactive(:,:,:)

REAL, INTENT(OUT) :: shno3(row_length,rows,model_levels)
REAL, INTENT(IN OUT) :: q(row_length,rows,model_levels) ! water vapour
REAL, INTENT(IN OUT) :: tracer(row_length,rows,model_levels,                   &
                               ntracers) ! tracer MMR

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
LOGICAL, INTENT(IN) :: have_nat(row_length,rows,model_levels)
! Stratospheric mask
LOGICAL, INTENT(IN) :: L_stratosphere(row_length,rows,model_levels)

! Flag for determining if this is the first chemistry call
LOGICAL, INTENT(IN) :: firstcall

! Local variables
INTEGER :: ix            ! dummy variable
INTEGER :: jy            ! dummy variable
INTEGER :: klevel        ! dummy variable
INTEGER :: js            ! loop variable
INTEGER :: jtr           ! loop variable - transported tracers
INTEGER :: jro2          ! loop variable - NTP RO2 species
INTEGER :: jna           ! loop variable, non-advected species
INTEGER :: jspf          ! loop variable - all active chemical species in f
INTEGER :: k             ! loop variable
INTEGER :: l             ! loop variable

INTEGER :: kcs           ! loop variable, start of current segment/chunk
INTEGER :: kce           ! loop variable, end of current segment/chunk

INTEGER :: ierr                               ! Error code: asad diags routines
INTEGER :: errcode                            ! Error code: ereport
CHARACTER(LEN=errormessagelength) :: cmessage ! Error message
CHARACTER(LEN=10) :: prods(2)                 ! Products
CHARACTER(LEN=10) :: prods3(3)                ! Products
LOGICAL :: ddmask(row_length,rows)            ! mask

REAL, ALLOCATABLE :: ystore(:,:,:) ! array for H2SO4 when updated in MODE
REAL :: zftr(row_length,rows,model_levels,jpcspf)
                                  ! array of chemically active species
                                  !   including RO2 species, in VMR
REAL :: zq(row_length,rows,model_levels) ! water vapour vmr
REAL :: co2_tmp(row_length,rows,model_levels) ! CO2
REAL :: zprt(row_length,rows,model_levels,jppj) ! photolysis rates for ASAD
REAL :: zdryrt2(row_length,rows,model_levels,jpdd) ! dry dep rate
REAL :: rc_het(row_length,rows,model_levels,2) ! heterog rates for trop chem

! 1-D mask for stratosphere
LOGICAL :: stratflag(row_length,rows,model_levels)

! Full ntp array
REAL :: ntp_data(row_length,rows,model_levels,dim_ntp)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Full-domain versions of ASAD modules variables
REAL :: full_dpd(row_length,rows,model_levels,jpspec)
REAL :: full_dpw(row_length,rows,model_levels,jpspec)
REAL :: full_fpsc1(row_length,rows,model_levels)
REAL :: full_fpsc2(row_length,rows,model_levels)
REAL :: full_prk(row_length,rows,model_levels,jpnr)
REAL :: full_y(row_length,rows,model_levels,jpspec)
REAL :: full_sph2o(row_length,rows,model_levels)
REAL :: full_za(row_length,rows,model_levels)
REAL :: full_rk(row_length,rows,model_levels,jpnr)
REAL :: full_tnd(row_length,rows,model_levels)
REAL :: full_sphno3(row_length,rows,model_levels)
REAL :: full_zftr(row_length,rows,model_levels,jpcspf)

! Chunked versions of ASAD modules variables
REAL :: chunk_dpd(chunk_n_x,chunk_n_y,chunk_n_z,jpspec)
REAL :: chunk_dpw(chunk_n_x,chunk_n_y,chunk_n_z,jpspec)
REAL :: chunk_fpsc1(chunk_n_x,chunk_n_y,chunk_n_z)
REAL :: chunk_fpsc2(chunk_n_x,chunk_n_y,chunk_n_z)
REAL :: chunk_prk(chunk_n_x,chunk_n_y,chunk_n_z,jpnr)
REAL :: chunk_y(chunk_n_x,chunk_n_y,chunk_n_z,jpspec)
REAL :: chunk_rk(chunk_n_x,chunk_n_y,chunk_n_z,jpnr)
REAL :: chunk_tnd(chunk_n_x,chunk_n_y,chunk_n_z)
REAL :: chunk_sphno3(chunk_n_x,chunk_n_y,chunk_n_z)
REAL :: chunk_zftr(chunk_n_x,chunk_n_y,chunk_n_z,jpcspf)

! Chunk start, end, iterator, and offset in each dimension
INTEGER :: xs, xe, xi, xo
INTEGER :: ys, ye, yi, yo
INTEGER :: zs, ze, zi, zo
INTEGER :: chunk_n_pnts
LOGICAL :: chunking_enabled

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMISTRY_CTL_FULL'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Dummy variables to satisfy expected numbers of arguments for ASAD_CDRIVE,
! ASAD_CHEMICAL_DIAGNOSTICS and ASAD_PSC_DIAGNOSTIC
ix=0
jy=0
klevel=0

! if heterogeneous chemistry is selected, allocate solid HNO3 array
IF (ukca_config%l_ukca_het_psc) THEN
  shno3 = 0.0
END IF

! Fill stratospheric flag indicator and ntp_data array
stratflag(:,:,:) = L_stratosphere(:,:,:)
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(l)                                        &
!$OMP SHARED(all_ntp, ntp_data, tot_n_pnts)
!$OMP DO SCHEDULE(DYNAMIC)
DO l = 1, dim_ntp
  IF (all_ntp(l)%l_required) THEN
    ntp_data(:,:,:,l) = all_ntp(l)%data_3d(:,:,:)
  END IF
END DO
!$OMP END DO
!$OMP END PARALLEL

! Store tracer mmr into for use in ASAD chemical solver
zq(:,:,:) = q(:,:,:)/c_h2o

! Store photolysis rates
IF (ukca_config%l_ukca_offline) THEN
  ! Offline chemistry has no photolysis
  zprt(:,:,:,:) = 0.0
ELSE
  zprt(:,:,:,:) = photol_rates(:,:,:,:)
END IF

! CO2 as species
IF (ANY(speci(:) == 'CO2       ')) THEN
  !  Copy the CO2 concentration into the asad module as VMR
  IF (ukca_config%l_chem_environ_co2_fld) THEN

    ! co2_interactive should be allocated if config option on
    IF (.NOT. ALLOCATED(co2_interactive)) THEN
      errcode = 1
      CALL ereport(ModuleName//':'//RoutineName, errcode,                      &
        'ERROR: co2_interactive array not allocated')
    END IF

    co2_tmp(:,:,:) = co2_interactive(:,:,:)/c_co2
  ELSE
    co2_tmp(:,:,:) = rmdi
  END IF
END IF

! Retrieve tropospheric heterogeneous rates from previous time step
IF (ukca_config%l_ukca_trophet) THEN
  ! N2O5
  l = name2ntpindex('het_n2o5  ')
  rc_het(:,:,:,1) = ntp_data(:,:,:,l)
  ! HO2+HO2
  l = name2ntpindex('het_ho2   ')
  rc_het(:,:,:,2) = ntp_data(:,:,:,l)
ELSE
  rc_het(:,:,:,:) = 0.0
END IF

zdryrt2(:,:,:,:) = 0.0
IF (ukca_config%l_ukca_intdd) THEN
  ! Interactive scheme extracts from levels in boundary layer
  DO k = 1,model_levels
    DO l=1,jpdd
      ddmask(:,:) = k <= nlev_with_ddep(:,:)
      WHERE (ddmask(:,:))
        zdryrt2(:,:,k,l) = zdryrt(:,:,l)
      END WHERE
    END DO
  END DO
ELSE    ! non-interactive
  zdryrt2(:,:,1,:) = zdryrt(:,:,:)
END IF

! Convert mmr into vmr for tracers.
! If running with nontransport RO2 species, data for the RO2 species
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
      full_zftr(:,:,:,jspf) = tracer(:,:,:,jtr)/c_species(jtr)
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
          full_zftr(:,:,:,jspf) = ntp_data(:,:,:,l) / c_na_species(jna)
        ELSE
          errcode = jro2
          WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl_full'
          CALL umPrint(umMessage,src='ukca_chemistry_ctl_full')
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
    WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl_full'
    CALL umPrint(umMessage,src='ukca_chemistry_ctl_full')
    cmessage = 'ERROR: Number of chemical active species /= jpro2+jpctr'
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
ELSE
  IF (jspf /= jpctr) THEN
    errcode = jspf
    WRITE(umMessage,'(A)') '** ERROR in ukca_chemistry_ctl_full'
    CALL umPrint(umMessage,src='ukca_chemistry_ctl_full')
    cmessage = 'ERROR: Number of chemical active species /= jpctr'
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
END IF

IF (.NOT. ALLOCATED(ystore) .AND. uph2so4inaer == 1)                           &
                             ALLOCATE(ystore(row_length,rows,model_levels))

! Copy water vapour and ice field into 1-D arrays
full_sph2o(:,:,:) = 0.0
IF (ukca_config%l_ukca_het_psc) THEN
  full_sph2o(:,:,:) = qcf(:,:,:) / c_h2o
END IF

! Fill SO4 surface area
full_za(:,:,:) = so4_sa(:,:,:)

! Number of grid points in every chunk
chunk_n_pnts = chunk_n_x * chunk_n_y * chunk_n_z
chunking_enabled = chunk_n_pnts /= tot_n_pnts

IF (uph2so4inaer == 1) THEN
  ! H2SO4 will be updated in MODE, so store old value here
  IF (ukca_config%l_fix_ukca_h2so4_ystore) THEN
     ! primary array passed is zftr, so save this, NOT y
    ystore(:,:,:) = full_zftr(:,:,:,istore_h2so4)
  ELSE IF (.not. chunking_enabled) THEN
    ! Preserve non-fixed behaviour when chunking is disabled
    ystore(:,:,:) = reshape(y(:,nn_h2so4), [row_length,rows,model_levels])
  END IF
END IF

! This next section of code calls ASAD routines to do chemistry integration
! with 3D chunking.  Unlike column-mode chunking, the chunk size is
! parameterised in all three dimensions (for maximum flexibility) and
! we do not resize ASAD arrays when the chunk size does not divide evenly
! into the full-domain size. Instead, we retract the starting point of the
! final chunk in any dimension to obtain a full chunk's worth of data. This
! leads to some redundant computation (some grid points form part of two
! different chunks) but, for small chunk sizes, this has a neglible impact
! on performance. As different chunks can be processed in parallel, care is
! taken to avoid data races where two different threads are computing the
! same grid point at the same time.

! 3D chunking loop
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) COLLAPSE(3)                  &
!$OMP PRIVATE(zi, ze, zs, zo, yi, ye, ys, yo, xi, xe, xs, xo,                  &
!$OMP         chunk_zftr, chunk_dpd, chunk_dpw, chunk_fpsc1, chunk_fpsc2,      &
!$OMP         chunk_prk, chunk_y, chunk_sphno3, chunk_tnd, chunk_rk)           &
!$OMP SHARED(chunk_n_pnts, chunk_n_x, chunk_n_y, chunk_n_z, model_levels,      &
!$OMP        rows, row_length, jpspec, full_sph2o, full_za, full_y, full_dpd,  &
!$OMP        full_dpw, full_fpsc1, full_fpsc2, full_prk, full_sphno3,          &
!$OMP        full_tnd, full_rk, pres, temp, zq, co2_tmp, cloud_frac, qcl,      &
!$OMP        zdryrt2, zwetrt, zprt, rc_het, have_nat, stratflag, H_plus,       &
!$OMP        zftr, full_zftr, ix, jy, k)
DO zi = 1, model_levels, chunk_n_z
  DO yi = 1, rows, chunk_n_y
    DO xi = 1, row_length, chunk_n_x
      ! Chunk start, end, and offset in Z dimension
      ze = min(model_levels, zi + (chunk_n_z - 1))
      zs = ze - chunk_n_z + 1
      zo = zi - zs + 1

      ! Chunk start, end, and offset in Y dimension
      ye = min(rows, yi + (chunk_n_y - 1))
      ys = ye - chunk_n_y + 1
      yo = yi - ys + 1

      ! Chunk start, end, and offset in X dimension
      xe = min(row_length, xi + (chunk_n_x - 1))
      xs = xe - chunk_n_x + 1
      xo = xi - xs + 1

      ! Copy from full-domain arrays to ASAD module variables
      sph2o(:) = reshape(full_sph2o(xs:xe,ys:ye,zs:ze), shape(sph2o))
      za(:) = reshape(full_za(xs:xe,ys:ye,zs:ze), shape(za))
      chunk_zftr(:,:,:,:) = full_zftr(xs:xe,ys:ye,zs:ze,:)

      CALL asad_cdrive(                                                        &
        chunk_zftr,                                                            &
        pres(xs:xe,ys:ye,zs:ze),                                               &
        temp(xs:xe,ys:ye,zs:ze),                                               &
        zq(xs:xe,ys:ye,zs:ze),                                                 &
        co2_tmp(xs:xe,ys:ye,zs:ze),                                            &
        cloud_frac(xs:xe,ys:ye,zs:ze),                                         &
        qcl(xs:xe,ys:ye,zs:ze),                                                &
        ix,                                                                    &
        jy,                                                                    &
        k,                                                                     &
        zdryrt2(xs:xe,ys:ye,zs:ze,:),                                          &
        zwetrt(xs:xe,ys:ye,zs:ze,:),                                           &
        rc_het(xs:xe,ys:ye,zs:ze,:),                                           &
        zprt(xs:xe,ys:ye,zs:ze,:),                                             &
        chunk_n_pnts,                                                          &
        have_nat(xs:xe,ys:ye,zs:ze),                                           &
        stratflag(xs:xe,ys:ye,zs:ze),                                          &
        H_plus(xs:xe,ys:ye,zs:ze))

      ! Copy from ASAD module variables to full-domain arrays
      zftr(xi:xe,yi:ye,zi:ze,:) = chunk_zftr(xo:,yo:,zo:,:)

      chunk_dpd(:,:,:,:) = reshape(dpd, shape(chunk_dpd))
      full_dpd(xi:xe,yi:ye,zi:ze,:) = chunk_dpd(xo:,yo:,zo:,:)

      chunk_dpw(:,:,:,:) = reshape(dpw, shape(chunk_dpw))
      full_dpw(xi:xe,yi:ye,zi:ze,:) = chunk_dpw(xo:,yo:,zo:,:)

      chunk_fpsc1(:,:,:) = reshape(fpsc1, shape(chunk_fpsc1))
      full_fpsc1(xi:xe,yi:ye,zi:ze) = chunk_fpsc1(xo:,yo:,zo:)

      chunk_fpsc2(:,:,:) = reshape(fpsc2, shape(chunk_fpsc2))
      full_fpsc2(xi:xe,yi:ye,zi:ze) = chunk_fpsc2(xo:,yo:,zo:)

      chunk_prk(:,:,:,:) = reshape(prk, shape(chunk_prk))
      full_prk(xi:xe,yi:ye,zi:ze,:) = chunk_prk(xo:,yo:,zo:,:)

      chunk_y(:,:,:,:) = reshape(y, shape(chunk_y))
      full_y(xi:xe,yi:ye,zi:ze,:) = chunk_y(xo:,yo:,zo:,:)

      chunk_sphno3(:,:,:) = reshape(sphno3, shape(chunk_sphno3))
      full_sphno3(xi:xe,yi:ye,zi:ze) = chunk_sphno3(xo:,yo:,zo:)

      chunk_tnd(:,:,:) = reshape(tnd, shape(chunk_tnd))
      full_tnd(xi:xe,yi:ye,zi:ze) = chunk_tnd(xo:,yo:,zo:)

      chunk_rk(:,:,:,:) = reshape(rk, shape(chunk_rk))
      full_rk(xi:xe,yi:ye,zi:ze,:) = chunk_rk(xo:,yo:,zo:,:)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (ukca_config%l_ukca_het_psc) THEN
  ! Save MMR of NAT PSC particles into 3-D array for PSC sedimentation.
  ! Note that sphno3 is NAT in number density of HNO3.
  IF (ANY(full_sphno3(:,:,:) > 0.0)) THEN
    shno3(:,:,:) = full_sphno3(:,:,:)/full_tnd(:,:,:)*c_hono2
  ELSE
    shno3(:,:,:) = 0.0
  END IF
END IF

IF (ukca_config%l_ukca_chem .AND. ukca_config%l_ukca_nr_aqchem) THEN
  ! Calculate chemical fluxes for MODE
  IF (ihso3_h2o2 > 0) THEN
    delSO2_wet_H2O2(:,:,:) = delSO2_wet_H2O2(:,:,:) +                          &
      full_rk(:,:,:,ihso3_h2o2) * full_y(:,:,:,nn_so2) *                       &
      full_y(:,:,:,nn_h2o2) * cdt_diag
  END IF
  IF (ihso3_o3 > 0) THEN
    delSO2_wet_O3(:,:,:) = delSO2_wet_O3(:,:,:) +                              &
      full_rk(:,:,:,ihso3_o3)*full_y(:,:,:,nn_so2)*full_y(:,:,:,nn_o3)*cdt_diag
  END IF
  IF (iso3_o3 > 0) THEN
    delSO2_wet_O3(:,:,:) = delSO2_wet_O3(:,:,:) +                              &
      full_rk(:,:,:,iso3_o3)*full_y(:,:,:,nn_so2)*full_y(:,:,:,nn_o3)*cdt_diag
  END IF
  ! net H2SO4 production - note that this is affected by
  ! l_fix_ukca_h2so4_ystore above. Y value is concentration
  ! from chemistry prior to zftr being over-written below
  IF (iso2_oh > 0 .AND. ih2so4_hv > 0) THEN
    delh2so4_chem(:,:,:) = delh2so4_chem(:,:,:) +                              &
      (full_rk(:,:,:,iso2_oh) * full_y(:,:,:,nn_so2) * full_y(:,:,:,nn_oh) -   &
       full_rk(:,:,:,ih2so4_hv) * full_y(:,:,:,nn_h2so4)) * cdt_diag
  ELSE IF (iso2_oh > 0) THEN
    delh2so4_chem(:,:,:) = delh2so4_chem(:,:,:) + full_rk(:,:,:,iso2_oh) *     &
      full_y(:,:,:,nn_so2) * full_y(:,:,:,nn_oh) * cdt_diag
  END IF

  IF (uph2so4inaer == 1) THEN
    ! Restore H2SO4 tracer as it will be updated in MODE
    ! using delh2so4_chem
    IF (ukca_config%l_fix_ukca_h2so4_ystore) THEN
      ! calculate delh2so4_chem as the difference in H2SO4 over chemistry
      ! zftr is already in VMR, so divide by diagnostic chemistry timestep to
      ! give as vmr/s
      delh2so4_chem(:,:,:) = (zftr(:,:,:,istore_h2so4) - ystore(:,:,:)) /      &
                               cdt_diag
      ! primary array passed is zftr, so copy back to this, NOT y
      zftr(:,:,:,istore_h2so4) = ystore(:,:,:)
    ELSE
      full_y(:,:,:,nn_h2so4) = ystore(:,:,:)
    END IF
  END IF
END IF

! 3D flux diagnostics
IF (L_asad_use_chem_diags .AND.                                                &
  ((L_asad_use_flux_rxns .OR. L_asad_use_rxn_rates) .OR.                       &
  (L_asad_use_wetdep .OR. L_asad_use_drydep))) THEN
  CALL asad_chemical_diagnostics(row_length,rows,model_levels,tot_n_pnts,      &
    full_dpd,full_dpw,full_prk,full_y,jy,ix,klevel,volume,ierr)
END IF

! PSC diagnostics
IF (L_asad_use_chem_diags .AND. L_asad_use_psc_diagnostic) THEN
  CALL asad_psc_diagnostic(row_length,rows,model_levels,tot_n_pnts,            &
                           full_fpsc1,full_fpsc2,jy,ix,klevel,ierr)
END IF

! Set SS species concentrations for output (stratospheric configurations)

! O1D mmr
IF (O1D_in_ss) THEN
  l = name2ntpindex('O(1D)     ')
  ntp_data(:,:,:,l) = full_y(:,:,:,nn_o1d)/full_tnd(:,:,:)*c_o1d
END IF

! O3P mmr
IF (O3P_in_ss) THEN
  l = name2ntpindex('O(3P)     ')
  ntp_data(:,:,:,l) = full_y(:,:,:,nn_o3p)/full_tnd(:,:,:)*c_o3p
END IF

! First copy the concentrations from the zftr array to the
! diag arrays and convert to moles.
! The indices (n_ch4, n_n2o etc.) don't refer to the location
! in the zftr array if RO2_NTP is true, but instead to the
! location in the tracer array.
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(jspf)                                     &
!$OMP SHARED(advt, atm_cf2cl2_mol, atm_cfcl3_mol, atm_ch4_mol, atm_co_mol,     &
!$OMP        atm_h2_mol, atm_mebr_mol, atm_n2o_mol, avogadro, jpcspf,          &
!$OMP        n_cf2cl2, n_cfcl3, n_ch4, n_co, n_h2, n_mebr, n_n2o, specf,       &
!$OMP        volume, zftr, full_tnd)
!$OMP DO SCHEDULE(DYNAMIC)
DO jspf = 1, jpcspf

  IF (n_ch4 > 0) THEN
    IF (specf(jspf) == advt(n_ch4)) THEN
      atm_ch4_mol(:,:,:) = zftr(:,:,:,jspf) * full_tnd(:,:,:) *                &
                             volume(:,:,:) * 1.0e6/avogadro
    END IF
  END IF

  ! CO
  IF (n_co > 0) THEN
    IF (specf(jspf) == advt(n_co)) THEN
      atm_co_mol(:,:,:) = zftr(:,:,:,jspf) * full_tnd(:,:,:) *                 &
                            volume(:,:,:) * 1.0e6/avogadro
    END IF
  END IF

  ! N2O
  IF (n_n2o > 0) THEN
    IF (specf(jspf) == advt(n_n2o)) THEN
      atm_n2o_mol(:,:,:) = zftr(:,:,:,jspf) * full_tnd(:,:,:) *                &
                             volume(:,:,:) * 1.0e6/avogadro
    END IF
  END IF

  ! CFC-12
  IF (n_cf2cl2 > 0) THEN
    IF (specf(jspf) == advt(n_cf2cl2)) THEN
      atm_cf2cl2_mol(:,:,:) = zftr(:,:,:,jspf) * full_tnd(:,:,:) *             &
                                volume(:,:,:) * 1.0e6/avogadro
    END IF
  END IF

  ! CFC-11
  IF (n_cfcl3 > 0) THEN
    IF (specf(jspf) == advt(n_cfcl3)) THEN
      atm_cfcl3_mol(:,:,:) = zftr(:,:,:,jspf) * full_tnd(:,:,:) *              &
                               volume(:,:,:) * 1.0e6/avogadro
    END IF
  END IF

  ! CH3Br
  IF (n_mebr > 0) THEN
    IF (specf(jspf) == advt(n_mebr)) THEN
      atm_mebr_mol(:,:,:) = zftr(:,:,:,jspf) * full_tnd(:,:,:) *               &
                              volume(:,:,:) * 1.0e6/avogadro
    END IF
  END IF

  ! H2
  IF (n_h2 > 0) THEN
    IF (specf(jspf) == advt(n_h2)) THEN
      atm_h2_mol(:,:,:) = zftr(:,:,:,jspf) * full_tnd(:,:,:) *                 &
                            volume(:,:,:) * 1.0e6/avogadro
    END IF
  END IF

END DO ! End loop through species in zftr array
!$OMP END DO
!$OMP END PARALLEL

IF (ALLOCATED(ystore)) DEALLOCATE(ystore)

! Bring results back from vmr to mmr.
! Also bring back results for nontransported RO2 to all_ntp
jspf = 0 ! reset counter for species in f array
DO js = 1,jpspec
  DO jtr=1,jpctr
    IF (advt(jtr) == speci(js)) THEN
      jspf = jspf+1
      tracer(:,:,:,jtr) = zftr(:,:,:,jspf) * c_species(jtr)
    END IF
  END DO

  ! If RO2 species are not being transported, map RO2 species
  ! back to the ntp_data array
  IF (ukca_config%l_ukca_ro2_ntp) THEN
    DO jro2 = 1, jpro2
      IF (spro2(jro2) == speci(js) .AND. ctype(js) == 'OO') THEN
        jspf = jspf+1
        ! Get index of speci(js) in NTP array,
        ! If not found this is a fatal error.
        l = name2ntpindex(spro2(jro2))

        ! Find location of species in nadvt array
        jna = nlnaro2(jro2)
        ntp_data(:,:,:,l) = zftr(:,:,:,jspf) * c_na_species(jna)
      END IF ! Close IF RO2 species
    END DO   ! Close loop through RO2 species
  END IF     ! Close IF RO2_NTP
END DO       ! Close loop through all species

! Map ntp_data back into 3D array
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(l)                                        &
!$OMP SHARED(all_ntp, model_levels, ntp_data, row_length, rows)
!$OMP DO SCHEDULE(DYNAMIC)
DO l = 1, dim_ntp
  IF (all_ntp(l)%l_required) THEN
    all_ntp(l)%data_3d(:,:,:) = ntp_data(:,:,:,l)
  END IF
END DO
!$OMP END DO
!$OMP END PARALLEL

! Preserve non-fixed behaviour when chunking is disabled
IF (.not. ukca_config%l_fix_ukca_h2so4_ystore .and.                            &
    .not. chunking_enabled) THEN
  y(:,:) = reshape(full_y, [tot_n_pnts,jpspec])
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chemistry_ctl_full

END MODULE ukca_chemistry_ctl_full_mod
