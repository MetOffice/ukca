! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!  Print inputs to UKCA to help with debugging
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
!    Language:  Fortran 2003
!
! ----------------------------------------------------------------------
MODULE ukca_pr_inputs_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: ukca_pr_inputs

CONTAINS

SUBROUTINE ukca_pr_inputs(l_first_call)

USE ukca_um_legacy_mod,   ONLY: mype
USE umPrintMgr,           ONLY: umMessage, UmPrint, PrintStatus,               &
                                PrStatus_Diag, PrStatus_Oper
USE ukca_config_specification_mod, ONLY: ukca_config, glomap_config,           &
                                         i_ukca_fastjx, i_light_param_pr,      &
                                         i_light_param_luhar,                  &
                                         i_light_param_ext

USE ukca_environment_fields_mod,   ONLY:                                       &
  conv_cloud_base, conv_cloud_top, kent, kent_dsc, land_sea_mask,              &
  soil_moisture_layer1, q, qcf, zbl, rough_length, seaice_frac,                &
  so4_aitken, so4_accum, co2_interactive, soot_fresh, soot_aged, ocff_fresh,   &
  ocff_aged, frac_types, tstar_tile, rho_r2, qcl,                              &
  exner_rho_levels, biogenic, fland, p_rho_levels, p_theta_levels, pstar,      &
  tstar, sea_salt_film, sea_salt_jet, rhokh_rdz, dtrdz,                        &
  we_lim, t_frac, zrzi, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc, u_s,           &
  ls_rain3d, ls_snow3d, rim_cry, rim_agg, autoconv, accretion, conv_rain3d,    &
  conv_snow3d, pv_on_theta_mlevs, latitude, longitude, sin_latitude,           &
  cos_latitude, tan_latitude, ext_cg_flash, ext_ic_flash,                      &
  dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6

USE ukca_ntp_mod,         ONLY: print_all_ntp


IMPLICIT NONE

LOGICAL, INTENT(IN) :: l_first_call

INTEGER :: k ! loop counter

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PR_INPUTS'

! Debug Non-transported prognostics every TS, depending on PrintStatus

IF (PrintStatus >= PrStatus_Diag) THEN

  ! print out info on NTP
  CALL print_all_ntp

END IF

IF (l_first_call .AND. PrintStatus >= PrStatus_Oper) THEN

  WRITE(umMessage,'(A)') ' ==========================='
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A)') '  MAX and MIN of UKCA INPUTS'
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A)') ' ==========================='
  CALL umPrint(umMessage,src=RoutineName)

  IF (ukca_config%l_ukca_intdd) THEN
    WRITE(umMessage,'(A,I8, 2E12.4)') 'soil_moisture_layer1: ',mype,           &
      MAXVAL(soil_moisture_layer1),MINVAL(soil_moisture_layer1)
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, 2E12.4)') 'frac_types: ',mype,                     &
      MAXVAL(frac_types),MINVAL(frac_types)
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, 2E12.4)') 'tstar_tile: ',mype,                     &
      MAXVAL(tstar_tile),MINVAL(tstar_tile)
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  WRITE(umMessage,'(A,I8, 2E12.4)') 'rough_length: ',mype,                     &
    MAXVAL(rough_length),MINVAL(rough_length)
  CALL umPrint(umMessage,src=RoutineName)

  IF (ukca_config%i_ukca_light_param == i_light_param_pr .OR.                  &
      ukca_config%i_ukca_light_param == i_light_param_luhar) THEN
    WRITE(umMessage,'(A,I8, 2I8)') 'conv_cloud_base: ',mype,                   &
      MAXVAL(conv_cloud_base),MINVAL(conv_cloud_base)
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, 2I8)') 'conv_cloud_top: ',mype,                    &
      MAXVAL(conv_cloud_top),MINVAL(conv_cloud_top)
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  WRITE(umMessage,'(A,L5)') 'land_sea_mask(1,1): ',land_sea_mask(1,1)
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I8, 2E12.4)') 'latitude: ',mype,                         &
    MAXVAL(latitude),MINVAL(latitude)
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I8, 2E12.4)') 'longitude: ',mype,                        &
    MAXVAL(longitude),MINVAL(longitude)
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I8, 2E12.4)') 'sin_latitude: ',mype,                     &
    MAXVAL(sin_latitude),MINVAL(sin_latitude)
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I8, 2E12.4)') 'cos_latitude: ',mype,                     &
    MAXVAL(cos_latitude),MINVAL(cos_latitude)
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I8, 2E12.4)') 'tan_latitude: ',mype,                     &
    MAXVAL(tan_latitude),MINVAL(tan_latitude)
  CALL umPrint(umMessage,src=RoutineName)

  IF (ukca_config%i_ukca_light_param == i_light_param_ext) THEN
    WRITE(umMessage,'(A,I8, 2E12.4)') 'ext_cg_flash: ',mype,                   &
      MAXVAL(ext_cg_flash),MINVAL(ext_cg_flash)
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, 2E12.4)') 'ext_ic_flash: ',mype,                   &
      MAXVAL(ext_ic_flash),MINVAL(ext_ic_flash)
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  DO k=1,ukca_config%model_levels

    WRITE(umMessage,'(A,I6)') 'LEVEL: ',k
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'p_rho_levels:   ',mype, k,          &
      MAXVAL(p_rho_levels(:,:,k)),MINVAL(p_rho_levels(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'p_theta_levels: ',mype,k,           &
      MAXVAL(p_theta_levels(:,:,k)),MINVAL(p_theta_levels(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)

    IF ((.NOT. ukca_config%l_suppress_ems) .OR.                                &
      glomap_config%l_ukca_fine_no3_prod .OR.                                  &
      glomap_config%l_ukca_coarse_no3_prod .OR.                                &
      (ukca_config%l_ukca_het_psc .AND. ukca_config%l_ukca_sa_clim .AND.       &
       ukca_config%l_use_classic_so4)) THEN
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'rho_r2:         ',mype,k,         &
        MAXVAL(rho_r2(:,:,k)),MINVAL(rho_r2(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    IF (ALLOCATED(co2_interactive)) THEN
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'co2_interactive: ',mype,k,        &
      MAXVAL(co2_interactive(:,:,k)),MINVAL(co2_interactive(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END IF
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'ls_rain3d:      ',mype,k,           &
      MAXVAL(ls_rain3d(:,:,k)),MINVAL(ls_rain3d(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'ls_snow3d:      ',mype,k,           &
      MAXVAL(ls_snow3d(:,:,k)),MINVAL(ls_snow3d(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'conv_rain3d:      ',mype,k,         &
      MAXVAL(conv_rain3d(:,:,k)),MINVAL(conv_rain3d(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'conv_snow3d:      ',mype,k,         &
      MAXVAL(conv_snow3d(:,:,k)),MINVAL(conv_snow3d(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    IF (.NOT. ukca_config%l_fix_tropopause_level) THEN
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'pv_on_theta_mlevs:      ',mype,k, &
        MAXVAL(pv_on_theta_mlevs(:,:,k)),MINVAL(pv_on_theta_mlevs(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    IF (glomap_config%l_6bin_dust_no3) THEN
      WRITE(umMessage,'(A30,2I5,3E12.4)')                                      &
        'CLASSIC dust bin 1 MMR: ',mype,k,                                     &
        MAXVAL(dust_div1 (:,:,k)),MINVAL(dust_div1(:,:,k)),                    &
        SUM(dust_div1 (:,:,k)) / SIZE(dust_div1(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A30,2I5,3E12.4)')                                      &
        'CLASSIC dust bin 2 MMR: ',mype,k,                                     &
        MAXVAL(dust_div2 (:,:,k)),MINVAL(dust_div2(:,:,k)),                    &
        SUM(dust_div2 (:,:,k)) / SIZE(dust_div2(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A30,2I5,3E12.4)')                                      &
        'CLASSIC dust bin 3 MMR: ',mype,k,                                     &
        MAXVAL(dust_div3 (:,:,k)),MINVAL(dust_div3(:,:,k)),                    &
        SUM(dust_div3 (:,:,k)) / SIZE(dust_div3(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A30,2I5,3E12.4)')                                      &
        'CLASSIC dust bin 4 MMR: ',mype,k,                                     &
        MAXVAL(dust_div4 (:,:,k)),MINVAL(dust_div4(:,:,k)),                    &
        SUM(dust_div4 (:,:,k)) / SIZE(dust_div4(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A30,2I5,3E12.4)')                                      &
        'CLASSIC dust bin 5 MMR: ',mype,k,                                     &
        MAXVAL(dust_div5 (:,:,k)),MINVAL(dust_div5(:,:,k)),                    &
        SUM(dust_div5 (:,:,k)) / SIZE(dust_div5(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A30,2I5,3E12.4)')                                      &
        'CLASSIC dust bin 6 MMR: ',mype,k,                                     &
        MAXVAL(dust_div6 (:,:,k)),MINVAL(dust_div6(:,:,k)),                    &
        SUM(dust_div6 (:,:,k)) / SIZE(dust_div6(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    IF (glomap_config%l_2bin_dust_no3) THEN
      WRITE(umMessage,'(A30,2I5,3E12.4)')                                      &
        'CLASSIC dust bin 1 MMR: ',mype,k,                                     &
        MAXVAL(dust_div1 (:,:,k)),MINVAL(dust_div1(:,:,k)),                    &
        SUM(dust_div1 (:,:,k)) / SIZE(dust_div1(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A30,2I5,3E12.4)')                                      &
        'CLASSIC dust bin 2 MMR: ',mype,k,                                     &
        MAXVAL(dust_div2 (:,:,k)),MINVAL(dust_div2(:,:,k)),                    &
        SUM(dust_div2 (:,:,k)) / SIZE(dust_div2(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    ! If heterogeneous chemistry on CLASSIC aerosols is used then
    ! print their corresponding MMRs or aerosol numbers.
    IF (ukca_config%l_ukca_classic_hetchem) THEN

      IF (ukca_config%l_use_classic_soot) THEN
        WRITE(umMessage,'(A30,2I5,3E12.4)')                                    &
          'CLASSIC soot fresh MMR: ',mype,k,                                   &
          MAXVAL(soot_fresh(:,:,k)),MINVAL(soot_fresh(:,:,k)),                 &
          SUM(soot_fresh(:,:,k)) / SIZE(soot_fresh(:,:,k))
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(umMessage,'(A30,2I5,3E12.4)')                                    &
          'CLASSIC soot aged MMR: ', mype, k,                                  &
          MAXVAL(soot_aged(:,:,k)),MINVAL(soot_aged(:,:,k)),                   &
          SUM(soot_aged(:,:,k)) / SIZE(soot_aged(:,:,k))
        CALL umPrint(umMessage,src=RoutineName)
      END IF

      IF (ukca_config%l_use_classic_ocff) THEN
        WRITE(umMessage,'(A30,2I5,3E12.4)')                                    &
          'CLASSIC OCFF fresh MMR: ',mype,k,                                   &
          MAXVAL(ocff_fresh(:,:,k)),MINVAL(ocff_fresh(:,:,k)),                 &
          SUM(ocff_fresh(:,:,k)) / SIZE(ocff_fresh(:,:,k))
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(umMessage,'(A30,2I5,3E12.4)')                                    &
          'CLASSIC OCFF aged MMR: ',mype,k,                                    &
          MAXVAL(ocff_aged(:,:,k)),MINVAL(ocff_aged(:,:,k)),                   &
          SUM(ocff_aged(:,:,k)) / SIZE(ocff_aged(:,:,k))
        CALL umPrint(umMessage,src=RoutineName)
      END IF

      IF (ukca_config%l_use_classic_biogenic) THEN
        WRITE(umMessage,'(A30,2I5,3E12.4)')                                    &
          'CLASSIC biogenic aerosol MMR: ',mype,k,                             &
          MAXVAL(biogenic(:,:,k)),MINVAL(biogenic(:,:,k)),                     &
          SUM(biogenic(:,:,k)) / SIZE(biogenic(:,:,k))
        CALL umPrint(umMessage,src=RoutineName)
      END IF

      IF (ukca_config%l_use_classic_seasalt) THEN
        WRITE(umMessage,'(A30,2I5,3E12.4)')                                    &
          'CLASSIC sea salt film no.: ',mype,k,                                &
          MAXVAL(sea_salt_film (:,:,k)),MINVAL(sea_salt_film(:,:,k)),          &
          SUM(sea_salt_film (:,:,k)) / SIZE(sea_salt_film(:,:,k))
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(umMessage,'(A30,2I5,3E12.4)')                                    &
          'CLASSIC sea salt jet no.: ',mype,k,                                 &
          MAXVAL(sea_salt_jet (:,:,k)),MINVAL(sea_salt_jet(:,:,k)),            &
          SUM(sea_salt_jet (:,:,k)) / SIZE(sea_salt_jet(:,:,k))
        CALL umPrint(umMessage,src=RoutineName)
      END IF
    END IF

  END DO   ! k

  ! wet variables
  DO k=1,ukca_config%model_levels
    WRITE(umMessage,'(A,I6)') 'WET LEVEL: ',k
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'q:              ',mype,k,           &
      MAXVAL(q(:,:,k)),MINVAL(q(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'qcl:            ',mype,k,           &
      MAXVAL(qcl(:,:,k)),MINVAL(qcl(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'qcf:            ',mype,k,           &
      MAXVAL(qcf(:,:,k)),MINVAL(qcf(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A12,2I5,3E12.4)') 'autoconv: ',mype,k,                   &
      MAXVAL(autoconv(:,:,k)),MINVAL(autoconv(:,:,k)),                         &
      SUM(autoconv(:,:,k)) / SIZE(autoconv(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A12,2I5,3E12.4)') 'accretion: ',mype,k,                  &
      MAXVAL(accretion(:,:,k)),MINVAL(accretion(:,:,k)),                       &
      SUM(accretion(:,:,k)) / SIZE(accretion(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A12,2I5,3E12.4)') 'rim_agg: ',mype,k,                    &
      MAXVAL(rim_agg(:,:,k)),MINVAL(rim_agg(:,:,k)),                           &
      SUM(rim_agg(:,:,k)) / SIZE(rim_agg(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A12,2I5,3E12.4)') 'rim_cry: ',mype,k,                    &
      MAXVAL(rim_cry(:,:,k)),MINVAL(rim_cry(:,:,k)),                           &
      SUM(rim_cry(:,:,k)) / SIZE(rim_cry(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A12,2I5,3E12.4)') 'ls_rain3d: ',mype,k,                  &
      MAXVAL(ls_rain3d(:,:,k)),MINVAL(ls_rain3d(:,:,k)),                       &
      SUM(ls_rain3d(:,:,k)) / SIZE(ls_rain3d(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A12,2I5,3E12.4)') 'ls_snow : ',mype,k,                   &
      MAXVAL(ls_snow3d(:,:,k)),MINVAL(ls_snow3d(:,:,k)),                       &
      SUM(ls_snow3d(:,:,k)) / SIZE(ls_snow3d(:,:,k))
    CALL umPrint(umMessage,src=RoutineName)
  END DO

  ! Exner - additional level at top
  IF ((.NOT. ukca_config%l_suppress_ems) .OR.                                  &
      glomap_config%l_ukca_fine_no3_prod .OR.                                  &
      glomap_config%l_ukca_coarse_no3_prod) THEN
    DO k=1,ukca_config%model_levels+1
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'exner_rho_levels: ',mype,k,       &
        MAXVAL(exner_rho_levels(:,:,k)),MINVAL(exner_rho_levels(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF

  ! boundary layer variables
  IF (.NOT. ukca_config%l_suppress_ems) THEN
    DO k=1,ukca_config%bl_levels
      WRITE(umMessage,'(A,I6)') 'BL LEVEL: ',k
      IF (k > 1) THEN
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'rhokh_rdz: ',mype,k,            &
          MAXVAL(rhokh_rdz(:,:,k)),MINVAL(rhokh_rdz(:,:,k))
        CALL umPrint(umMessage,src=RoutineName)
      END IF
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'dtrdz: ',mype,k,                  &
        MAXVAL(dtrdz(:,:,k)),MINVAL(dtrdz(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF

  ! Sulfate aerosol inputs to photolysis schemes
  IF (ukca_config%l_use_classic_so4 .AND.                                      &
      ukca_config%i_ukca_photol == i_ukca_fastjx) THEN
    DO k=1,ukca_config%model_levels
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'so4_aitken: ',mype,k,             &
        MAXVAL(so4_aitken(:,:,k)),MINVAL(so4_aitken(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'so4_accum: ',mype,k,              &
        MAXVAL(so4_accum(:,:,k)),MINVAL(so4_accum(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF

  IF (.NOT. ukca_config%l_suppress_ems) THEN
    WRITE(umMessage,'(A,I8,2I8)') 'kent:     ',mype,                           &
      MAXVAL(kent),MINVAL(kent)
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8,2I8)') 'kent_dsc: ',mype,                           &
      MAXVAL(kent_dsc),MINVAL(kent_dsc)
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,I8,2E12.4)') 'zhsc:     ',mype,                        &
      MAXVAL(zhsc),MINVAL(zhsc)
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  WRITE(umMessage,'(A,I8,2E12.4)') 'zbl: ',mype,                               &
    MAXVAL(zbl),MINVAL(zbl)
  CALL umPrint(umMessage,src=RoutineName)

  ! Loop over vegetation types considered in boundary layer mixing (tr_mix)
  IF (.NOT. ukca_config%l_suppress_ems) THEN
    DO k=1,ukca_config%nlev_ent_tr_mix
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'we_lim: ',mype,k,                 &
        MAXVAL(we_lim(:,:,k)),MINVAL(we_lim(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 't_frac: ',mype,k,                 &
        MAXVAL(t_frac(:,:,k)),MINVAL(t_frac(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'zrzi: ',mype,k,                   &
        MAXVAL(zrzi(:,:,k)),MINVAL(zrzi(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'we_lim_dsc: ',mype,k,             &
        MAXVAL(we_lim_dsc(:,:,k)),MINVAL(we_lim_dsc(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 't_frac_dsc: ',mype,k,             &
        MAXVAL(t_frac_dsc(:,:,k)),MINVAL(t_frac_dsc(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'zrzi_dsc: ',mype,k,               &
        MAXVAL(zrzi_dsc(:,:,k)),MINVAL(zrzi_dsc(:,:,k))
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF

  WRITE(umMessage,'(A,I8,2E12.4)') 'pstar: ',mype,MAXVAL(pstar),MINVAL(pstar)
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I8,2E12.4)') 'tstar: ',mype,MAXVAL(tstar),MINVAL(tstar)
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I8,2E12.4)') 'fland: ',mype,MAXVAL(fland),MINVAL(fland)
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I8,2E12.4)') 'u_s  : ',mype,MAXVAL(u_s),MINVAL(u_s)
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A16,I5,2E12.4)') 'seaice_frac: ',mype,                     &
    MAXVAL(seaice_frac), MINVAL(seaice_frac)
  CALL umPrint(umMessage,src=RoutineName)

END IF  ! l_first_call .AND. PrintStatus >= PrStatus_Oper

END SUBROUTINE ukca_pr_inputs

END MODULE ukca_pr_inputs_mod
