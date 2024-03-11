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
!  Setup routine for chemistry
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_CHEMISTRY_CTL, UKCA_CHEMISTRY_CTL_COL and
!   UKCA_CHEMISTRY_CTL_TROPRAQ.
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
MODULE ukca_chemistry_setup_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CHEMISTRY_SETUP_MOD'

CONTAINS

SUBROUTINE ukca_chemistry_setup(row_length, rows, model_levels,                &
                                theta_field_size, bl_levels, ntracers, ntype,  &
                                npft, len_stashwork, i_month, i_day_number,    &
                                i_hour, land_points, land_index, tile_pts,     &
                                tile_index, istore_h2so4, nlev_with_ddep,      &
                                r_minute, secs_per_step, latitude, longitude,  &
                                sinlat, tanlat, temp, pres, rh,                &
                                p_layer_boundaries, t_surf, dzl, u_s, z0m,     &
                                drain, crain, frac_types, seaice_frac, stcon,  &
                                surf_hf, soilmc_lp, fland, laift_lp,           &
                                canhtft_lp, z0tile_lp, t0tile_lp, zbl, H_plus, &
                                zdryrt, zwetrt, tracer, stashwork, firstcall)

USE asad_findreaction_mod,         ONLY: asad_findreaction
USE asad_mod,                      ONLY: advt, ih2so4_hv, iso2_oh, jpcspf,     &
                                         jpdd, jpdw, jppj, jpspec, jpspj,      &
                                         jpspt, jptk, ndepd, ndepw, nldepd,    &
                                         nprkx, ntrkx, specf, speci, spj, spt

USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_conserve_mod,             ONLY: ukca_conserve
USE ukca_cspecies,                 ONLY: n_h2o, n_h2o2, n_h2so4, nn_cl
USE ukca_ddepctl_mod,              ONLY: ukca_ddepctl
USE ukca_ddeprt_mod,               ONLY: ukca_ddeprt
USE ukca_environment_fields_mod,   ONLY: surf_wetness, h2o2_offline
USE ukca_wdeprt_mod,               ONLY: ukca_wdeprt
USE ukca_um_legacy_mod,            ONLY: deposition_from_ukca_chemistry

USE ereport_mod,                   ONLY: ereport
USE errormessagelength_mod,        ONLY: errormessagelength
USE parkind1,                      ONLY: jprb, jpim
USE umPrintMgr,                    ONLY: umMessage, umPrint
USE yomhook,                       ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN)     :: row_length       ! size of UKCA x dimension
INTEGER, INTENT(IN)     :: rows             ! size of UKCA y dimension
INTEGER, INTENT(IN)     :: model_levels     ! size of UKCA z dimension
INTEGER, INTENT(IN)     :: theta_field_size ! no. of points in horizontal
INTEGER, INTENT(IN)     :: bl_levels        ! no. of boundary layer levels
INTEGER, INTENT(IN)     :: ntracers         ! no. of tracers
INTEGER, INTENT(IN)     :: ntype            ! no. of surface types
INTEGER, INTENT(IN)     :: npft             ! no. of plant functional types
INTEGER, INTENT(IN)     :: len_stashwork
INTEGER, INTENT(IN)     :: i_month
INTEGER, INTENT(IN)     :: i_day_number
INTEGER, INTENT(IN)     :: i_hour
INTEGER, INTENT(IN)     :: land_points
INTEGER, INTENT(IN)     :: land_index(land_points)
INTEGER, INTENT(IN)     :: tile_pts(ntype)
INTEGER, INTENT(IN)     :: tile_index(land_points,ntype)
INTEGER, INTENT(OUT)    :: istore_h2so4  ! location of H2SO4 in f array
INTEGER, INTENT(OUT)    :: nlev_with_ddep(row_length,rows)

REAL,    INTENT(IN)     :: r_minute
REAL,    INTENT(IN)     :: secs_per_step ! timestep in seconds
REAL,    INTENT(IN)     :: latitude(row_length,rows)
REAL,    INTENT(IN)     :: longitude(row_length,rows)
REAL,    INTENT(IN)     :: sinlat(row_length,rows)             ! sin(latitude)
REAL,    INTENT(IN)     :: tanlat(row_length,rows)             ! tan(latitude)
REAL,    INTENT(IN)     :: temp(row_length,rows,model_levels)  ! temperature
REAL,    INTENT(IN)     :: pres(row_length,rows,model_levels)  ! pressure
REAL,    INTENT(IN)     :: rh(row_length,rows,model_levels)    ! RH frac
REAL,    INTENT(IN)     :: p_layer_boundaries(row_length,rows,                 &
                                              0:model_levels)  ! pressure
REAL,    INTENT(IN)     :: t_surf(row_length,rows)             ! surface temp
REAL,    INTENT(IN)     :: dzl(row_length,rows,bl_levels)      ! thickness
REAL,    INTENT(IN)     :: u_s(row_length,rows)                ! u-star
REAL,    INTENT(IN)     :: z0m(row_length,rows)                ! roughness
REAL,    INTENT(IN)     :: drain(row_length,rows,model_levels) ! 3-D LS rain
REAL,    INTENT(IN)     :: crain(row_length,rows,model_levels) ! 3-D convec
REAL,    INTENT(IN)     :: frac_types(land_points,ntype)
REAL,    INTENT(IN)     :: seaice_frac(row_length,rows)
REAL,    INTENT(IN)     :: stcon(row_length,rows,npft)
REAL,    INTENT(IN)     :: surf_hf(row_length,rows)
REAL,    INTENT(IN)     :: soilmc_lp(land_points)
REAL,    INTENT(IN)     :: fland(land_points)
REAL,    INTENT(IN)     :: laift_lp(land_points,npft)
REAL,    INTENT(IN)     :: canhtft_lp(land_points,npft)
REAL,    INTENT(IN)     :: z0tile_lp(land_points,ntype)
REAL,    INTENT(IN)     :: t0tile_lp(land_points,ntype)
REAL,    INTENT(IN)     :: zbl(row_length,rows)
REAL,    INTENT(IN)     :: H_plus(row_length,rows,model_levels) ! 3-D pH array
REAL,    INTENT(OUT)    :: zdryrt(row_length,rows,jpdd)         ! dry dep rate
REAL,    INTENT(OUT)    :: zwetrt(row_length,rows,model_levels,                &
                                  jpdw)                         ! wet dep rate
REAL,    INTENT(IN OUT) :: tracer(row_length,rows,model_levels,                &
                                  ntracers)                     ! tracer MMR
REAL,    INTENT(IN OUT) :: stashwork(len_stashwork) ! diagnostics array

LOGICAL, INTENT(IN)     :: firstcall ! is this the first chemistry call?

! Local variables

INTEGER :: jspf

LOGICAL :: tropraq

! The call to ukca_conserve requires a logical to be set.
! ukca_conserve calculates and conserves total chlorine, bromine, and
! hydrogen. For these elements closed chemistry should be prescribed.
! Called before chemistry, before_chem, it calculates
! total bromine, chlorine, and hydrogen as 3-D fields.
LOGICAL, PARAMETER :: before_chem = .TRUE.

INTEGER                           :: errcode   ! Error code: ereport
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message
CHARACTER(LEN=10)                 :: prods(2)  ! Pair of products
CHARACTER(LEN=10)                 :: prods3(3) ! Triplet of products

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'UKCA_CHEMISTRY_SETUP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

tropraq = (ukca_config%l_ukca_trop .OR. ukca_config%l_ukca_aerchem .OR.        &
           ukca_config%l_ukca_raq .OR. ukca_config%l_ukca_raqaero)

IF (firstcall .AND. (.NOT. ukca_config%l_ukca_offline_be)) THEN

  !         Check whether water vapour is advective tracer. Then,
  !         check whether UM and ASAD advective tracers correspond
  !         to each other.

  IF ((ukca_config%l_ukca_advh2o) .AND. (n_h2o == 0)) THEN
    errcode = 4
    cmessage='No tracer for advected water vapour'
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF

  IF (.NOT. tropraq) THEN

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
      iso2_oh = asad_findreaction('SO2       ', 'OH        ',                  &
                                  prods, 2, spt, ntrkx, jptk+1, jpspt)

      IF (iso2_oh == 0) THEN   ! check for stratospheric sulphur chemistry
        ! product should be HO2 not H2O
        ! not considered when l_fix_ukca_h2so4_ystore=T
        prods = ['SO3       ','HO2       ']
        iso2_oh = asad_findreaction('SO2       ', 'OH        ',                &
                                    prods, 2, spt, ntrkx, jptk+1, jpspt)
        IF (ukca_config%i_ukca_chem_version >= 121) THEN
           ! additional product at version=121
          prods3 = ['SO3       ','OH        ','H         ']
          ih2so4_hv = asad_findreaction('H2SO4     ', 'PHOTON    ',            &
                                        prods3, 3, spj, nprkx, jppj+1, jpspj)
        ELSE
          prods = ['SO3       ','OH        ']
          ih2so4_hv = asad_findreaction('H2SO4     ', 'PHOTON    ',            &
                                        prods, 2, spj, nprkx, jppj+1, jpspj)
        END IF

      END IF

      IF (iso2_oh == 0 .AND. ih2so4_hv == 0) THEN
        errcode = 1
        cmessage=' Sulphur chemistry reactions not found'
        WRITE(umMessage,'(A)') cmessage
        CALL umPrint(umMessage,src='ukca_chemistry_ctl')
        WRITE(umMessage,'(A,I0,A,I0)') 'iso2_oh: ',iso2_oh,' ih2so4_hv: ',     &
                                       ih2so4_hv
        CALL umPrint(umMessage,src='ukca_chemistry_ctl')
        CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
      END IF
    END IF   ! l_ukca_nr_aqchem

  END IF

END IF  ! of initialization of chemistry subroutine (firstcall)

! Call routine to calculate dry deposition rates.
zdryrt  = 0.0
IF (ndepd /= 0 .AND. .NOT. ukca_config%l_ukca_drydep_off) THEN

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

    ELSE IF (tropraq .AND. (.NOT. ukca_config%l_ukca_offline_be)) THEN

      CALL ukca_ddepctl(row_length, rows, bl_levels, ntype, npft,              &
        land_points, land_index, tile_pts, tile_index,                         &
        secs_per_step, sinlat, frac_types, t_surf,                             &
        p_layer_boundaries(:,:,0), dzl, zbl, surf_hf, u_s,                     &
        rh, stcon, soilmc_lp, fland, seaice_frac, laift_lp,                    &
        canhtft_lp, z0tile_lp, t0tile_lp,                                      &
        nlev_with_ddep, zdryrt, len_stashwork, stashwork)

    ELSE

      CALL ukca_ddepctl(row_length, rows, bl_levels, ntype, npft,              &
        land_points, land_index, tile_pts, tile_index,                         &
        secs_per_step, sinlat, frac_types, t_surf,                             &
        p_layer_boundaries(:,:,0), dzl, zbl, surf_hf, u_s,                     &
        rh(:,:,1), stcon, soilmc_lp, fland, seaice_frac, laift_lp,             &
        canhtft_lp, z0tile_lp, t0tile_lp,                                      &
        nlev_with_ddep, zdryrt, len_stashwork, stashwork)

    END IF

  ELSE                             ! Call prescribed dry dep

    CALL ukca_ddeprt(theta_field_size, bl_levels, i_day_number, i_month,       &
                     i_hour, r_minute, secs_per_step, longitude, latitude,     &
                     tanlat, dzl, z0m, u_s, t_surf, zdryrt)

  END IF
END IF

! Call routine to calculate wet deposition rates.
zwetrt = 0.0
IF (ndepw /= 0 .AND. .NOT. ukca_config%l_ukca_wetdep_off) THEN
  CALL ukca_wdeprt(theta_field_size, model_levels, drain, crain, temp,         &
                   latitude, secs_per_step, zwetrt, H_plus)
END IF

IF (ukca_config%l_ukca_offline_be) RETURN

IF (ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_stratcfc .OR.             &
    ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat) THEN

  ! Calculate total chlorine and total bromine before chemistry
  IF (nn_cl > 0) THEN
    CALL ukca_conserve(row_length, rows, model_levels, ntracers,               &
                       tracer, pres, drain, crain, before_chem)
  END IF

END IF    ! l_ukca_strat etc

! Reduce over-prediction of H2O2 using ancillary value.
IF (ukca_config%l_ukca_offline) THEN
  WHERE (tracer(:,:,:,n_h2o2) > h2o2_offline(:,:,:))                           &
    tracer(:,:,:,n_h2o2) = h2o2_offline(:,:,:)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_chemistry_setup

END MODULE ukca_chemistry_setup_mod
