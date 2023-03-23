! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module for offline oxidants chemistry with the Backward-Euler solver.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds,
!  University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Contained subroutines:
!   ukca_chemistry_ctl_be
!   ukca_chemco_be_offline
!   ukca_deriv_be_offline
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_chemistry_ctl_be_mod

IMPLICIT NONE
PRIVATE

REAL, ALLOCATABLE :: rc(:,:)                ! 1-D Rate coefficient array

INTEGER, ALLOCATABLE, SAVE :: ibimol(:)     ! index array for bimolecular rxns
INTEGER, ALLOCATABLE, SAVE :: itrimol(:)    ! index array for termolecular rxns
INTEGER, ALLOCATABLE, SAVE :: ihetero(:)    ! index array for heterogeneous rxns

REAL, SAVE    :: frac_dms_so2      ! product yield DMS + OH => SO2 + DMSO
REAL, SAVE    :: frac_dms_dmso     ! product yield DMS + OH => SO2 + DMSO
REAL, SAVE    :: frac_dmso_so2     ! product yield DMSO + OH => SO2
REAL, SAVE    :: frac_monoterp_sec_org(3)
                                   ! product yield Monoterp + X => Sec_Org

PUBLIC ukca_chemistry_ctl_be

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CHEMISTRY_CTL_BE_MOD'

CONTAINS

! ##############################################################################

SUBROUTINE ukca_chemistry_ctl_be(                                              &
  row_length, rows, model_levels, bl_levels,                                   &
  theta_field_size,                                                            &
  ntracers, ntype, npft,                                                       &
  i_month, i_day_number, i_hour,                                               &
  r_minute, chem_timestep,                                                     &
  k_be_top,                                                                    &
  latitude,                                                                    &
  longitude,                                                                   &
  sinlat,                                                                      &
  tanlat,                                                                      &
  p_theta_levels, temp, q,                                                     &
  qcl, rh,                                                                     &
  p_layer_boundaries,                                                          &
  tracer,                                                                      &
  t_surf, dzl, z0m, u_s,                                                       &
  drain, crain,                                                                &
  cloud_frac,                                                                  &
  volume,                                                                      &
  land_points, land_index,                                                     &
  tile_pts, tile_index, frac_types,                                            &
  zbl, surf_hf, seaice_frac, stcon,                                            &
  soilmc_lp, fland, laift_lp, canhtft_lp,                                      &
  z0tile_lp, t0tile_lp,                                                        &
  uph2so4inaer,                                                                &
  delso2_wet_h2o2,                                                             &
  delso2_wet_o3,                                                               &
  delh2so4_chem,                                                               &
  delso2_drydep,                                                               &
  delso2_wetdep,                                                               &
  len_stashwork,                                                               &
  stashwork,                                                                   &
  H_plus_3d_arr                                                                &
  )

USE ukca_um_legacy_mod,   ONLY: deposition_from_ukca_chemistry
USE asad_mod,             ONLY: speci, y, tnd, ndepd, ndepw, advt, prk,        &
                                nitfg, lvmr, f, jpctr, jpspec,                 &
                                jpbk, jptk, jphk, jpnr, jpdd, jpdw, nldepd
USE asad_chem_flux_diags, ONLY: asad_chemical_diagnostics,                     &
                                l_asad_use_chem_diags, l_asad_use_flux_rxns,   &
                                l_asad_use_wetdep, l_asad_use_drydep
USE asad_ftoy_mod,        ONLY: asad_ftoy
USE ukca_cspecies,        ONLY: c_species, nn_h2so4, nn_so2, n_h2o2
USE ukca_constants,       ONLY: c_h2o
USE ukca_environment_fields_mod,   ONLY: h2o2_offline, surf_wetness
USE ukca_config_specification_mod, ONLY: ukca_config

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper

USE ukca_be_drydep_mod, ONLY: ukca_be_drydep
USE ukca_be_wetdep_mod, ONLY: ukca_be_wetdep
USE ukca_ddepctl_mod, ONLY: ukca_ddepctl
USE ukca_ddeprt_mod, ONLY: ukca_ddeprt
USE ukca_wdeprt_mod, ONLY: ukca_wdeprt
USE ukca_drydep_mod, ONLY: ukca_drydep
USE ukca_wetdep_mod, ONLY: ukca_wetdep

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
INTEGER, INTENT(IN) :: k_be_top          ! top level for integration

!       Variables for interactive dry deposition scheme
INTEGER, INTENT(IN) :: land_points
INTEGER, INTENT(IN) :: land_index(land_points)
INTEGER, INTENT(IN) :: tile_pts(ntype)
INTEGER, INTENT(IN) :: tile_index(land_points,ntype)

INTEGER, INTENT(IN) :: chem_timestep   ! chemical timestep

REAL, INTENT(IN) :: r_minute                           ! minute
REAL, INTENT(IN) :: latitude(row_length,rows)          ! lautitude (degrees)
REAL, INTENT(IN) :: longitude(row_length,rows)         ! longitude (degrees)
REAL, INTENT(IN) :: sinlat(row_length, rows)           ! sin(latitude)
REAL, INTENT(IN) :: tanlat(row_length, rows)           ! tan(latitude)
REAL, INTENT(IN) :: p_theta_levels(row_length,rows,model_levels) ! pressure
REAL, INTENT(IN) :: p_layer_boundaries(row_length,rows,0:model_levels)
                                                       ! pressure
REAL, INTENT(IN) :: temp(row_length,rows,model_levels) ! actual temp
REAL, INTENT(IN) :: dzl(row_length, rows, bl_levels)   ! thickness
REAL, INTENT(IN) :: u_s(row_length, rows)              ! ustar
REAL, INTENT(IN) :: z0m(row_length, rows)              ! roughness
REAL, INTENT(IN) :: t_surf(row_length, rows)           ! surface temp
REAL, INTENT(IN) :: drain(row_length,rows,model_levels)   ! 3-D LS rain
REAL, INTENT(IN) :: crain(row_length,rows,model_levels)   ! 3-D convec
REAL, INTENT(IN) :: volume(row_length,rows,model_levels)  ! cell vol.
REAL, INTENT(IN) :: qcl(row_length, rows, model_levels)   ! qcl
REAL, INTENT(IN) :: rh(row_length, rows, model_levels)    ! RH frac
REAL, INTENT(IN) :: cloud_frac(row_length, rows, model_levels)

!  Variables for interactive dry deposition scheme

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

REAL, INTENT(IN OUT) :: q(row_length,rows,model_levels)   ! water vapour
REAL, INTENT(IN OUT) :: tracer(row_length,rows,model_levels,ntracers)
                                                         ! tracer MMR

! SO2 increments
REAL, INTENT(IN OUT) :: delso2_wet_h2o2(row_length,rows,                       &
  model_levels)
REAL, INTENT(IN OUT) :: delso2_wet_o3(row_length,rows,model_levels)
REAL, INTENT(IN OUT) :: delh2so4_chem(row_length,rows,model_levels)
REAL, INTENT(IN OUT) :: delso2_drydep(row_length,rows,model_levels)
REAL, INTENT(IN OUT) :: delso2_wetdep(row_length,rows,model_levels)

! Diagnostics array
INTEGER, INTENT(IN) :: len_stashwork

REAL, INTENT(IN OUT) :: stashwork (len_stashwork)

! 3D pH array
REAL, INTENT(IN) :: H_plus_3d_arr(row_length, rows, model_levels)

! Chemical fluxes
LOGICAL             :: lflux             ! true when flux requests found
REAL, ALLOCATABLE   :: dflux(:,:)        ! chemical fluxes on each level
                                         ! (molecules / (cm3.s) )


! Local variables
INTEGER :: nlev_in_bl(row_length, rows)     ! No levs in bl
INTEGER :: nlev_in_bl2(theta_field_size)    ! No levs in bl

INTEGER, SAVE :: nr            ! no of rxns for BE
INTEGER, SAVE :: n_be_calls    ! no of call to BE solver

INTEGER :: ix                  ! dummy variable
INTEGER :: jy                  ! dummy variable
INTEGER :: i                   ! Loop variable
INTEGER :: j                   ! loop variable
INTEGER :: js                  ! loop variable
INTEGER :: k                   ! loop variable for model levels
INTEGER :: l                   ! loop variable
INTEGER :: n_pnts              ! number of points (= theta_field_size)
INTEGER :: jit                 ! dummy iteration count

INTEGER           :: ierr                     ! Error code: asad diags routines
INTEGER           :: errcode                  ! Error code: ereport
CHARACTER(LEN=72) :: cmessage                 ! Error message
LOGICAL           :: blmask(theta_field_size) ! mask
LOGICAL, SAVE     :: ofirst = .TRUE.          ! True for first call of asad_ftoy
                                              ! in current timestep
REAL, SAVE      :: dts                        ! Backward Euler timestep

REAL :: secs_per_step                         ! chemical time step

! Dummy variables for compatability with column call approach
REAL :: dpd_dummy(model_levels,jpspec)
REAL :: dpw_dummy(model_levels,jpspec)
REAL :: prk_dummy(model_levels,jpnr)
REAL :: y_dummy(model_levels,jpspec)

! SO2 increments in molecules/cm^3
REAL :: so2_wetox_h2o2(theta_field_size)
REAL :: so2_wetox_o3(theta_field_size)
REAL :: so2_dryox_oh(theta_field_size)

REAL, ALLOCATABLE :: ystore(:)        ! to store H2SO4 when updated in MODE
REAL :: zftr(theta_field_size,jpctr)  ! 1-D array of tracers
REAL :: zp(theta_field_size)          ! 1-D pressure
REAL :: zt(theta_field_size)          ! 1-D temperature
REAL :: zq(theta_field_size)          ! 1-D water vapour
REAL :: zclw(theta_field_size)        ! 1-D cloud liquid water
REAL :: zfcloud(theta_field_size)     ! 1-D cloud fraction
REAL :: zdryrt(row_length, rows, jpdd)                ! dry dep rate
REAL :: zdryrt2(theta_field_size, jpdd)               ! dry dep rate
REAL :: zwetrt(row_length, rows, model_levels, jpdw)  ! wet dep rate
REAL :: zwetrt2(theta_field_size, jpdw)               ! wet dep rat
REAL :: wetrt(theta_field_size,jpspec)         ! wet dep rates (s-1)
REAL :: dryrt(theta_field_size,jpspec)         ! dry dep rates (s-1)
REAL :: H_plus_2d_arr(theta_field_size,model_levels) ! 2-D pH array to use in
                                                     ! wet deposition
REAL :: H_plus_1d_arr(theta_field_size)  ! 1-D pH array to use in BE solver

LOGICAL, SAVE :: firstcall = .TRUE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMISTRY_CTL_BE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

n_pnts = theta_field_size
nr = jpbk + jptk + jphk
secs_per_step = REAL(chem_timestep)

IF (.NOT. ALLOCATED(dflux)) ALLOCATE(dflux(theta_field_size,nr))

! ystore is required when H2SO4 (gas-phase) is updated in MODE
IF (.NOT. ALLOCATED(ystore) .AND. uph2so4inaer == 1)                           &
  ALLOCATE(ystore(theta_field_size))

! dummy variables required for consistency with column call approach
ix=0
jy=0
dpd_dummy=0.0
dpw_dummy=0.0
prk_dummy=0.0
y_dummy=0.0

IF (firstcall) THEN

  ! Backward Euler timestep variables

  n_be_calls = INT(secs_per_step/ukca_config%dts0)
  ! Ensure we call BE at least once
  IF (n_be_calls == 0) n_be_calls = 1
  ! Calculate the BE Timestep
  dts = secs_per_step/n_be_calls

  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'(A12,I7,A6,E10.1)') 'n_be_calls: ',n_be_calls,' dts: ',dts
    CALL umPrint(umMessage,src='ukca_chemistry_ctl_be')
  END IF

  IF (.NOT. lvmr) THEN
    errcode = 1
    cmessage = ' lvmr must be set to true for BE offline chemistry'
    CALL ereport('UKCA_CHEMCO_BE_OFFLINE',errcode,cmessage)
  END IF

END IF  ! of initialization of chemistry subroutine (firstcall)

zdryrt  = 0.0
zdryrt2 = 0.0
IF (ndepd /= 0 .AND. (.NOT. ukca_config%l_ukca_drydep_off)) THEN

  IF (ukca_config%l_ukca_intdd) THEN          ! Call interactive dry deposition

    IF (ukca_config%l_deposition_jules) THEN   ! Use JULES-based routines

      CALL deposition_from_ukca_chemistry(                                     &
        secs_per_step, bl_levels, row_length, rows, ntype, npft,               &
        jpspec, ndepd, nldepd, speci,                                          &
        land_points, land_index, tile_pts, tile_index,                         &
        seaice_frac, fland, sinlat,                                            &
        p_layer_boundaries(:,:,0), rh(:,:,1), t_surf, surf_hf, surf_wetness,   &
        z0tile_lp, stcon, laift_lp, canhtft_lp, t0tile_lp,                     &
        soilmc_lp, zbl, dzl, frac_types, u_s,                                  &
        zdryrt, nlev_in_bl, len_stashwork, stashwork)

    ELSE                                      ! Use exising UKCA routines

      CALL ukca_ddepctl(row_length, rows, bl_levels, ntype, npft,              &
        land_points, land_index, tile_pts, tile_index,                         &
        secs_per_step, sinlat, frac_types, t_surf,                             &
        p_layer_boundaries(:,:,0), dzl, zbl, surf_hf, u_s,                     &
        rh, stcon, soilmc_lp, fland, seaice_frac, laift_lp,                    &
        canhtft_lp, z0tile_lp, t0tile_lp,                                      &
        nlev_in_bl, zdryrt, len_stashwork, stashwork)

    END IF

  ELSE                             ! Call prescribed dry deposition

    CALL ukca_ddeprt(n_pnts, bl_levels, i_day_number, i_month, i_hour,         &
                     r_minute, secs_per_step, longitude, latitude, tanlat,     &
                     dzl, z0m, u_s, t_surf, zdryrt)

  END IF
END IF

! Call routine to calculate wet deposition rates.

zwetrt  = 0.0
zwetrt2 = 0.0
IF (ndepw /= 0 .AND. .NOT. ukca_config%l_ukca_wetdep_off) THEN

  ! Reshape 3D H_plus array to 2D to use in calculating Wet Deposition rates
  DO k=1,model_levels
    H_plus_2d_arr(:,k) = RESHAPE(H_plus_3d_arr(:,:,k),[theta_field_size])
  END DO

  CALL ukca_wdeprt(n_pnts, model_levels, drain, crain, temp,                   &
                   latitude, secs_per_step, zwetrt, H_plus_2d_arr)
END IF

! Need this line here for interactive dry deposition.
nlev_in_bl2(:) = RESHAPE(nlev_in_bl(:,:),[theta_field_size])

IF (.NOT. ALLOCATED(rc)) ALLOCATE(rc(theta_field_size,nr))

! Model levels loop
DO k=1,k_be_top

  zdryrt2(:,:) = 0.0e0
  IF (ukca_config%l_ukca_intdd) THEN
    ! Interactive scheme extracts from levels in boundary layer
    blmask(:) = (k <= nlev_in_bl2(:))
    DO l=1,jpdd
      WHERE (blmask(:))
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

  ! Put pressure, temperature and tracer mmr into 1-D arrays
  ! for use in ASAD chemical solver

  zp(:) = RESHAPE(p_theta_levels(:,:,k),[theta_field_size])
  zt(:) = RESHAPE(temp(:,:,k),[theta_field_size])

  IF (ukca_config%l_fix_ukca_offox_h2o_fac) THEN
    ! convert q from mmr to vmr as required by chemical solver
    zq(:) = RESHAPE(q(:,:,k),[theta_field_size])/c_h2o
  ELSE
    zq(:) = RESHAPE(q(:,:,k),[theta_field_size])
  END IF

  IF (k <= model_levels) THEN
    zclw(:) = RESHAPE(qcl(:,:,k),[theta_field_size])
    zfcloud(:) = RESHAPE(cloud_frac(:,:,k),[theta_field_size])
  ELSE
    zclw(:) = 0.0
    zfcloud(:) = 0.0
  END IF

  ! Convert mmr into vmr for tracers and set f array
  DO js=1,jpctr
    zftr(:,js) = RESHAPE(tracer(:,:,k,js),[theta_field_size])/c_species(js)
  END DO

  ! reshape 3D H_plus values to 1-D for use in BE solver
  H_plus_1d_arr(:) = RESHAPE(H_plus_3d_arr(:,:,k),[theta_field_size])

  ! Backward Euler with non-families

  ! Calculate total number  density, o2, h2o, and tracer
  !  concentrations for Backward Euler solver

  DO l=1,jpdw
    zwetrt2(:,l) = RESHAPE(zwetrt(:,:,k,l),[theta_field_size])
  END DO

  ! Assign wet and dry deposition rates to species

  CALL ukca_be_wetdep(n_pnts, zwetrt2, wetrt)

  CALL ukca_be_drydep(k, n_pnts, nlev_in_bl2, zdryrt2, dryrt)

  ! Fill the asad arrays for wet and dry deposition diagnostics
  IF ( ndepw /= 0 ) CALL ukca_wetdep(zwetrt2, n_pnts)

  IF ( ndepd /= 0 ) CALL ukca_drydep(k, zdryrt2, n_pnts)


  ! Calculate reaction rate coefficients (rc)

  CALL ukca_chemco_be_offline(theta_field_size, zt, zp, zq, zfcloud, zclw,     &
                              H_plus_1d_arr)

  ! Convert tracers to concentration
  DO js=1,jpctr
    f(:,js) = zftr(:,js) * tnd(:)
  END DO

  ! Initialise y array, including the offline oxidants
  jit = 0
  CALL asad_ftoy(ofirst, nitfg, jit, n_pnts, ix, jy, k)

  !  Call Backward Euler solver
  !   N.B. Emissions already added, via call to TR_MIX from UKCA_EMISSION_CTL

  ! Set the flux logical according to the STASH requests found previously
  lflux = (l_asad_use_chem_diags .AND. (l_asad_use_flux_rxns  .OR.             &
           l_asad_use_wetdep .OR. l_asad_use_drydep))

  ! Store H2SO4 tracer if it will be updated in MODE using delh2so4_chem
  IF (uph2so4inaer == 1) ystore(:) = y(:,nn_h2so4)

  CALL ukca_deriv_offline(nr, n_be_calls, n_pnts, dts, y, dryrt, wetrt,        &
                          so2_wetox_h2o2, so2_wetox_o3, so2_dryox_oh, lflux,   &
                          dflux)

  ! Restore H2SO4 tracer as it will be updated in MODE using delh2so4_chem
  IF (uph2so4inaer == 1) y(:,nn_h2so4) = ystore(:)

  ! Retrieve tracer concentrations
  DO j = 1,jpctr
    i_loop: DO i = 1,jpspec
      IF (advt(j) == speci(i)) THEN
        tracer(:,:,k,j) = RESHAPE(y(:,i)/tnd(:),[row_length,rows])*            &
                               c_species(j)
        EXIT i_loop
      END IF
    END DO i_loop
  END DO

  ! Update the 3-D SO2 flux arrays (molecules cm-3 per timestep)
  delso2_wet_h2o2(:,:,k) = RESHAPE(so2_wetox_h2o2(:),[row_length,rows])
  delso2_wet_o3(:,:,k)   = RESHAPE(so2_wetox_o3(:),[row_length,rows])
  delh2so4_chem(:,:,k)   = RESHAPE(so2_dryox_oh(:),[row_length,rows])
  delso2_drydep(:,:,k)   = RESHAPE(dryrt(:,nn_so2)*                            &
                                 y(:,nn_so2),[row_length,rows]) * dts
  delso2_wetdep(:,:,k)   = RESHAPE(wetrt(:,nn_so2)*                            &
                                 y(:,nn_so2),[row_length,rows]) * dts

  ! Fill the prk array from the flux array
  IF (lflux) THEN
    DO i = 1,jpbk
      prk(:,ibimol(i)) = dflux(:,i)
    END DO

    j = jpbk
    DO i = 1,jptk
      prk(:,itrimol(i)) = dflux(:,i+j)
    END DO

    j = jpbk + jptk
    DO i = 1,jphk
      prk(:,ihetero(i)) = dflux(:,i+j)
    END DO

    ! 3D flux diagnostics
    CALL asad_chemical_diagnostics(row_length, rows, model_levels,             &
    dpd_dummy, dpw_dummy, prk_dummy, y_dummy, ix, jy, k,                       &
                                   volume, ierr)

  END IF    ! lflux

END DO    ! level loop (k)

IF (ALLOCATED(rc)) DEALLOCATE(rc)
IF (ALLOCATED(dflux)) DEALLOCATE(dflux)
IF (ALLOCATED(ystore)) DEALLOCATE(ystore)

! Reduce over-prediction of H2O2 using ancillary value.
WHERE (tracer(:,:,:,n_h2o2) > h2o2_offline(:,:,:))                             &
       tracer(:,:,:,n_h2o2) = h2o2_offline(:,:,:)

IF (firstcall) firstcall = .FALSE.


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chemistry_ctl_be

! ##############################################################################

SUBROUTINE ukca_chemco_be_offline(theta_field_size, zt, zp, zq, zfcloud,       &
                                  zclw, H_plus_1d_arr)

USE asad_bimol_mod,      ONLY: asad_bimol
USE asad_findreaction_mod, ONLY: asad_findreaction
USE asad_hetero_mod,     ONLY: asad_hetero
USE asad_totnud_mod,     ONLY: asad_totnud
USE asad_trimol_mod,     ONLY: asad_trimol
USE ukca_chem_defs_mod,  ONLY: ratb_defs, ratt_defs, rath_defs
USE asad_mod,            ONLY: p, t, wp, spb, spt, sph, nbrkx, ntrkx, nhrkx,   &
                               rk, jpspb, jpspt, jpsph, jpbk, jptk, jphk
USE missing_data_mod,    ONLY: rmdi, imdi
USE umPrintMgr,          ONLY: umPrint, umMessage
USE ereport_mod,         ONLY: ereport
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN) :: theta_field_size  ! No of points

REAL, INTENT(IN)  :: zt(1:theta_field_size)         ! Temperature (K)
REAL, INTENT(IN)  :: zp(1:theta_field_size)         ! Pressure (Pa)
REAL, INTENT(IN)  :: zq(1:theta_field_size)         ! Water vapour
REAL, INTENT(IN)  :: zfcloud(1:theta_field_size)    ! Cloud fraction
REAL, INTENT(IN)  :: zclw(1:theta_field_size)       ! Cloud liquid water (kg/kg)
! 1-D pH array to be passed to asad_hetero
REAL, INTENT(IN) :: H_plus_1d_arr(1:theta_field_size)

! Local variables

INTEGER :: i                             ! counter
INTEGER :: j                             ! counter
INTEGER :: npnts                         ! no of points
INTEGER :: errcode                       ! error code

REAL :: dummy(theta_field_size,2)        ! dummy variable for asad_hetero call

CHARACTER(LEN=10) :: r1                  ! 1st reactant name
CHARACTER(LEN=10) :: r2                  ! 2nd reactant name
CHARACTER(LEN=10) :: p1                  ! 1st product name
CHARACTER(LEN=10) :: p2                  ! 2nd product name
CHARACTER(LEN=10) :: prods(2)            ! product names

CHARACTER(LEN=72) :: cmessage            ! Error message

LOGICAL, SAVE     :: first = .TRUE.      ! True for first call to routine

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMCO_BE_OFFLINE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Copy pressure, temperature and water vapour into module variables
npnts = theta_field_size

p(1:npnts) = zp(1:npnts)
t(1:npnts) = zt(1:npnts)
wp(1:npnts) = zq(1:npnts)

! Calculate total number density (tnd)

CALL asad_totnud(npnts)

! Calculate rate coefficients

CALL asad_bimol(npnts)

CALL asad_trimol(npnts)

CALL asad_hetero(npnts, zfcloud, zclw, dummy, H_plus_1d_arr)

IF (first) THEN

  IF (.NOT. ALLOCATED(ibimol)) ALLOCATE(ibimol(jpbk))
  IF (.NOT. ALLOCATED(itrimol)) ALLOCATE(itrimol(jptk))
  IF (.NOT. ALLOCATED(ihetero)) ALLOCATE(ihetero(jphk))

  ibimol(:) = imdi
  itrimol(:) = imdi
  ihetero(:) = imdi

  frac_dms_so2 = rmdi
  frac_dms_dmso = rmdi
  frac_dmso_so2 = rmdi
  frac_monoterp_sec_org(3) = rmdi
  rc(:,:) = rmdi

  ! Locate reactions and read product yields
  errcode = 0
  DO i = 1,jpbk
    r1 = ratb_defs(i)%react1
    r2 = ratb_defs(i)%react2
    p1 = ratb_defs(i)%prod1
    p2 = ratb_defs(i)%prod2
    prods = [p1,p2]
    ibimol(i) = asad_findreaction(r1, r2, prods, 2, spb, nbrkx, jpbk+1, jpspb)

    IF (ibimol(i) == imdi) errcode = i

    IF (p1(1:4) == 'SO2 ' .AND. p2(1:4) == 'DMSO') THEN
      ! DMS + OH => SO2 + DMSO
      frac_dms_so2 = ratb_defs(i)%pyield1
      frac_dms_dmso = ratb_defs(i)%pyield2
    END IF
    IF (r1(1:4) == 'DMSO') THEN
      frac_dmso_so2 = ratb_defs(i)%pyield1
    END IF
    IF (r1(1:8) == 'Monoterp' .AND. r2(1:4) == 'OH  ') THEN
      frac_monoterp_sec_org(1) = ratb_defs(i)%pyield1
    END IF
    IF (r1(1:8) == 'Monoterp' .AND. r2(1:4) == 'O3  ') THEN
      frac_monoterp_sec_org(2) = ratb_defs(i)%pyield1
    END IF
    IF (r1(1:8) == 'Monoterp' .AND. r2(1:4) == 'NO3 ') THEN
      frac_monoterp_sec_org(3) = ratb_defs(i)%pyield1
    END IF
  END DO

  IF (errcode /= 0) THEN
    cmessage = 'Unidentified reaction: '//ratb_defs(errcode)%react1//          &
               ratb_defs(errcode)%react2
    CALL ereport('UKCA_CHEMCO_BE_OFFLINE',errcode,cmessage)
  END IF

  IF (ABS(frac_dms_so2 - rmdi) < EPSILON(0.0) .OR.                             &
      ABS(frac_dms_dmso - rmdi) < EPSILON(0.0) .OR.                            &
      ABS(frac_dmso_so2 - rmdi) < EPSILON(0.0) .OR.                            &
      ANY(ABS(frac_monoterp_sec_org - rmdi) < EPSILON(0.0))) THEN
    errcode = 1
    cmessage = ' One or more product yields are undefined'
    CALL ereport('UKCA_CHEMCO_BE_OFFLINE',errcode,cmessage)
  END IF

  errcode = 0
  DO i = 1,jptk
    r1 = ratt_defs(i)%react1
    r2 = ratt_defs(i)%react2
    p1 = ratt_defs(i)%prod1
    p2 = ratt_defs(i)%prod2
    prods = [p1,p2]
    itrimol(i) = asad_findreaction(r1, r2, prods, 2, spt, ntrkx, jptk+1, jpspt)

    IF (itrimol(i) == imdi) errcode = i
  END DO

  IF (errcode /= 0) THEN
    cmessage = 'Unidentified reaction: '//ratt_defs(errcode)%react1//          &
               ratt_defs(errcode)%react2
    CALL ereport('UKCA_CHEMCO_BE_OFFLINE',errcode,cmessage)
  END IF

  errcode = 0
  DO i = 1,jphk
    r1 = rath_defs(i)%react1
    r2 = rath_defs(i)%react2
    p1 = rath_defs(i)%prod1
    p2 = rath_defs(i)%prod2
    prods = [p1,p2]
    ihetero(i) = asad_findreaction(r1, r2, prods, 2, sph, nhrkx, jphk+1, jpsph)

    IF (ihetero(i) == imdi) errcode = i
  END DO

  IF (errcode /= 0) THEN
    cmessage = 'Unidentified reaction: '//rath_defs(errcode)%react1//          &
               rath_defs(errcode)%react2
    CALL ereport('UKCA_CHEMCO_BE_OFFLINE',errcode,cmessage)
  END IF

END IF     ! first

! Fill the backward-Euler rate coefficient array from rk, using the
!  addressing arrays

DO i = 1,jpbk
  rc(:,i) = rk(:,ibimol(i))
END DO

j = jpbk
DO i = 1,jptk
  rc(:,i+j) = rk(:,itrimol(i))
END DO

j = jpbk + jptk
DO i = 1,jphk
  rc(:,i+j) = rk(:,ihetero(i))
END DO

IF (first .AND. ANY(ABS(rc - rmdi) < EPSILON(0.0))) THEN
  cmessage = 'Missing rate coefficient data in array rc'
  WRITE(umMessage,'(A55,I2)') cmessage//' at location: ',MAXLOC(rc)
  CALL umPrint(umMessage,src='ukca_chemco_be_offline')
  errcode = 1
  CALL ereport('UKCA_CHEMCO_BE_OFFLINE',errcode,cmessage)
END IF

IF (first) THEN
  IF (frac_dms_so2 > 1.0 .OR. frac_dms_dmso > 1.0 .OR.                         &
      frac_dmso_so2 > 1.0 .OR. ANY(frac_monoterp_sec_org(:) > 1.0)) THEN
    cmessage = 'A yield fraction > 1.0 was encountered'
    errcode = 1
    CALL ereport('UKCA_CHEMCO_BE_OFFLINE',errcode,cmessage)
  END IF
END IF

first = .FALSE.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chemco_be_offline

! ##############################################################################

! Description:
!  To perform a chemical integration using the Backward Euler solver
!   for the offline chemistry mechanism
!
! Method:
!  Backward-Euler integration (no iterations are required)
!
! ----------------------------------------------------------------------

SUBROUTINE ukca_deriv_offline(nr, n_be_calls, theta_field_size, dts, y,        &
                              dryrt, wetrt, so2_wetox_h2o2, so2_wetox_o3,      &
                              so2_dryox_oh, lflux, dflux)

! This solver uses the following equations from module UKCA_CHEM_OFFLINE:

! rc(1) : DMS + OH => SO2
! rc(2) : DMS + OH => f1*SO2 + f2*DMSO
! rc(3) : DMS + NO3 => SO2
! rc(4) : DMSO + OH => f3*SO2
! rc(5) : Monoterpene + OH => Sec_Org
! rc(6) : Monoterpene + O3 => Sec_Org
! rc(7) : Monoterpene + NO3 => Sec_Org
! rc(8) : HO2 + HO2 => H2O2
! rc(9) : H2O2 + OH => H2O
! rc(10) : SO2 + OH => H2SO4
! rc(11) : SO2 + H2O2(aq) => NULL0  (HSO3- + H2O2(aq) => SO4)
! rc(12) : SO2 + O3(aq) => NULL1  (HSO3- + O3(aq) => SO4)
! rc(13) : SO2 + O3(aq) => NULL2  (SO3= + O3(aq) => SO4)

! where f1 = frac_dms_so2
!       f2 = frac_dms_dmso
!       f3 = frac_dmso_so2

USE asad_mod, ONLY: jpspec
USE ukca_um_legacy_mod, ONLY: oneover_v
USE yomhook,         ONLY: lhook, dr_hook
USE parkind1,        ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN) :: nr                            ! No. reactions
INTEGER, INTENT(IN) :: n_be_calls                    ! No. chemical steps
INTEGER, INTENT(IN) :: theta_field_size              ! No. of points

REAL, INTENT(IN)    :: dts                           ! timestep
REAL, INTENT(IN)    :: dryrt(theta_field_size,jpspec) ! dry deposition rates
REAL, INTENT(IN)    :: wetrt(theta_field_size,jpspec) ! dry deposition rates
LOGICAL, INTENT(IN) :: lflux                         ! true for chemical fluxes
! SO2 increments in molecules/cm^3
REAL, INTENT(OUT)   :: so2_wetox_h2o2(theta_field_size)
REAL, INTENT(OUT)   :: so2_wetox_o3(theta_field_size)
REAL, INTENT(OUT)   :: so2_dryox_oh(theta_field_size)
REAL, INTENT(OUT)   :: dflux(theta_field_size,nr)    ! reaction fluxes in
                                                     ! molecules/(cm3.s)

REAL, INTENT(IN OUT) :: y(theta_field_size,jpspec)    ! species concentrations

!     Local variables

INTEGER :: j                               ! loop counter
INTEGER :: n                               ! loop counter for chemical timestep
INTEGER :: npnts                           ! number of points

REAL :: p(theta_field_size)                ! production rate
REAL :: l(theta_field_size)                ! loss rate
REAL :: yp(theta_field_size,jpspec)        ! concentrations at iteration start
REAL :: tmp1(theta_field_size)             ! work array for vector reciprocal
REAL :: tmp2(theta_field_size)             ! work array for vector reciprocal

REAL :: f_inv                              ! inverse of n_be_calls

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DERIV_OFFLINE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

npnts = theta_field_size

! Initialise the flux terms
so2_wetox_h2o2(:) = 0.0
so2_wetox_o3(:) = 0.0
so2_dryox_oh(:) = 0.0
dflux(:,:) = 0.0

DO n = 1, n_be_calls           ! loop over chemical timesteps

  DO j = 1,jpspec
    yp(:,j) = y(:,j)
  END DO


  ! Note that there is no need for iterations in this mechanism

  ! DMS          Y( 2)
  p(:) = 0.0
  l(:) = dryrt(:,2) + wetrt(:,2) +                                             &
         ((rc(:,1) + rc(:,2)) * y(:,9)) + (rc(:,3) * y(:,10))
  y(:, 3) = (yp(:, 3)+dts*p(:))/(1.0+dts*l(:))

  tmp2(:) = 1.0+dts*l(:)
  CALL oneover_v(npnts, tmp2, tmp1)

  y(:, 2) = (yp(:, 2)+dts*p(:))*tmp1

  ! DMSO         Y( 5)
  p(:) = (rc(:,2) * y(:,2) * y(:,9) * frac_dms_dmso)
  l(:) = dryrt(:,5) + wetrt(:,5) +                                             &
         (rc(:,4) * y(:,9))

  tmp2(:) = 1.0+dts*l(:)
  CALL oneover_v(npnts, tmp2, tmp1)

  y(:, 5) = (yp(:, 5)+dts*p(:))*tmp1

  ! SO2          Y( 3)
  p(:) = ((rc(:,1) + rc(:,2) * frac_dms_so2) * y(:,2) * y(:,9)) +              &
         (rc(:,3) * y(:,2) * y(:,10)) +                                        &
         (rc(:,4) * y(:,5) * y(:,9) * frac_dmso_so2)
  l(:) = dryrt(:,3) + wetrt(:,3) +                                             &
         (rc(:,10) * y(:,9)) + (rc(:,11) * y(:,1)) +                           &
         ((rc(:,12) + rc(:,13)) * y(:,8))

  tmp2(:) = 1.0+dts*l(:)
  CALL oneover_v(npnts, tmp2, tmp1)

  y(:, 3) = (yp(:, 3)+dts*p(:))*tmp1

  ! H2SO4        Y( 4)
  p(:) = (rc(:,10) * y(:,3) * y(:,9))
  l(:) = dryrt(:,4) + wetrt(:,4)

  tmp2(:) = 1.0+dts*l(:)
  CALL oneover_v(npnts, tmp2, tmp1)

  y(:, 4) = (yp(:, 4)+dts*p(:))*tmp1

  ! Monoterp     Y( 6)
  p(:) = 0.0
  l(:) = dryrt(:,6) + wetrt(:,6) +                                             &
         (rc(:,5) * y(:,9)) + (rc(:,6) * y(:,8)) + (rc(:,7) * y(:,10))

  tmp2(:) = 1.0+dts*l(:)
  CALL oneover_v(npnts, tmp2, tmp1)

  y(:, 6) = (yp(:, 6)+dts*p(:))*tmp1

  ! Sec_org      Y( 7)
  p(:) = (rc(:,5) * y(:,6) * y(:,9) * frac_monoterp_sec_org(1)) +              &
         (rc(:,6) * y(:,6) * y(:,8) * frac_monoterp_sec_org(2)) +              &
         (rc(:,7) * y(:,6) * y(:,10) * frac_monoterp_sec_org(3))
  l(:) = dryrt(:,7) + wetrt(:,7)

  tmp2(:) = 1.0+dts*l(:)
  CALL oneover_v(npnts, tmp2, tmp1)

  y(:, 7) = (yp(:, 7)+dts*p(:))*tmp1

  ! H2O2         Y( 1)
  p(:) = rc(:,8) * y(:,11) * y(:,11)
  l(:) = dryrt(:,1) + wetrt(:,1) + (rc(:,9) * y(:,9)) + (rc(:,11) * y(:,3))

  tmp2(:) = 1.0+dts*l(:)
  CALL oneover_v(npnts, tmp2, tmp1)

  y(:, 1) = (yp(:, 1)+dts*p(:))*tmp1


  ! Calculate flux terms at end of iteration.

  ! Fluxes to aqueous sulphate (always required for MODE):
  so2_wetox_h2o2(:)= so2_wetox_h2o2(:) +                                       &
                                 rc(:,11) * y(:,3) * y(:,1) * dts
  so2_wetox_o3(:)  = so2_wetox_o3(:) +                                         &
                               ((rc(:,12) + rc(:,13)) * y(:,3) * y(:,8)) * dts
  so2_dryox_oh(:)  = so2_dryox_oh(:) +                                         &
                                 rc(:,10) * y(:,3) * y(:,9) * dts

  ! Chemical fluxes (only required if flux diagnostics are requested)
  IF (lflux) THEN
    f_inv = 1.0/n_be_calls
    dflux(:,1) = dflux(:,1) + rc(:,1)*y(:,2)*y(:,9)*f_inv      ! DMS + OH
    dflux(:,2) = dflux(:,1) + rc(:,2)*y(:,2)*y(:,9)*f_inv      ! DMS + OH
    dflux(:,3) = dflux(:,3) + rc(:,3)*y(:,2)*y(:,10)*f_inv     ! DMS + NO3
    dflux(:,4) = dflux(:,4) + rc(:,4)*y(:,5)*y(:,9)*f_inv      ! DMSO + OH
    dflux(:,5) = dflux(:,5) + rc(:,5)*y(:,6)*y(:,9)*f_inv      ! MonoTerp + OH
    dflux(:,6) = dflux(:,6) + rc(:,6)*y(:,6)*y(:,8)*f_inv      ! MonoTerp + O3
    dflux(:,7) = dflux(:,7) + rc(:,7)*y(:,6)*y(:,10)*f_inv     ! MonoTerp + NO3
    dflux(:,8) = dflux(:,8) + rc(:,8)*y(:,11)*y(:,11)*f_inv    ! HO2 + HO2
    dflux(:,9) = dflux(:,9) + rc(:,9)*y(:,1)*y(:,9)*f_inv      ! H2O2 + OH
    dflux(:,10) = dflux(:,10) + rc(:,10)*y(:,3)*y(:,9)*f_inv   ! SO2 + OH
    dflux(:,11) = dflux(:,11) + rc(:,11)*y(:,3)*y(:,1)*f_inv   ! SO2 + H2O2(aq)
    dflux(:,12) = dflux(:,12) + rc(:,12)*y(:,3)*y(:,8)*f_inv   ! SO2 + O3(aq)
    dflux(:,13) = dflux(:,13) + rc(:,13)*y(:,3)*y(:,8)*f_inv   ! SO2 + O3(aq)

  END IF    ! lflux

END DO  ! n_be_calls


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_deriv_offline

END MODULE ukca_chemistry_ctl_be_mod
