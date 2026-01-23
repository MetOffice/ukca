! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Calculate the complex refractive index of an aerosol mixture in a mode
!  and for a given wavelength
!
! Method:
!  By default the effective refractive index of the aerosol mixture is
!  assumed to be the volume-weighted mean of the refractive index of the
!  chemical components present. However, if the Maxwell-Garnett approximation
!  is to be used (l_ukca_radaer_mg_mix = .TRUE.), then the mixing is done
!  in two stages. Firstly, the non-BC components being mixed together first
!  using volume-weighting. Secondly, the BC is mixed with the non-BC mixture
!  via the Maxwell-Garnett mixing rule.
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
MODULE ukca_radaer_ri_calc_mod

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook


IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName = 'UKCA_RADAER_RI_CALC_MOD'

CONTAINS

SUBROUTINE ukca_radaer_ri_calc(                                                &
     ! From the structure ukca_radaer for UKCA/radiation interaction
     nmodes, ncp_max, ncp_max_x_nmodes,                                        &
     i_cpnt_index, i_cpnt_type, n_cpnt_in_mode,                                &
     l_nitrate, l_soluble, l_sustrat,                                          &
     ! Refractive index
     refr_real, refr_imag,                                                     &
     ! Modal properties
     ukca_cpnt_volume, ukca_modal_volume, ukca_water_volume,                   &
     ! Indicies for arrays
     i_mode, n_ukca_cpnt,                                                      &
     ! Stratospheric aerosol treated as sulphuric acid?
     l_in_stratosphere,                                                        &
     ! Integer control switches
     i_ukca_tune_bc, i_glomap_clim_tune_bc, i_ukca_radaer_prescribe_ssa,       &
     ! Output refractive index real and imag parts
     re_m, im_m )

USE ukca_radaer_precalc,    ONLY:                                              &
     npd_ukca_maxcomptype

USE ukca_radaer_struct_mod, ONLY:                                              &
     ip_ukca_h2so4,                                                            &
     ip_ukca_water

USE ukca_mode_setup,        ONLY:                                              &
     cp_su, cp_bc

USE ukca_option_mod,       ONLY:                                               &
     do_not_prescribe

!
! Arguments
!

! From ukca_radaer Structure for UKCA/radiation interaction
INTEGER, INTENT(IN) :: nmodes
INTEGER, INTENT(IN) :: ncp_max
INTEGER, INTENT(IN) :: ncp_max_x_nmodes
INTEGER, INTENT(IN) :: i_cpnt_index( ncp_max, nmodes )
INTEGER, INTENT(IN) :: i_cpnt_type( ncp_max_x_nmodes )
INTEGER, INTENT(IN) :: n_cpnt_in_mode( nmodes )
LOGICAL, INTENT(IN) :: l_nitrate
LOGICAL, INTENT(IN) :: l_soluble( nmodes )
LOGICAL, INTENT(IN) :: l_sustrat

! Array dimensions
INTEGER, INTENT(IN) ::  n_ukca_cpnt
!
! Array index
INTEGER, INTENT(IN) :: i_mode
!
! Refractive index for the various aerosol components
REAL, INTENT(IN) :: refr_real (npd_ukca_maxcomptype)
REAL, INTENT(IN) :: refr_imag (npd_ukca_maxcomptype)

! Component volumes
REAL,    INTENT(IN) :: ukca_cpnt_volume (n_ukca_cpnt)
!
! Modal volumes
REAL, INTENT(IN) :: ukca_modal_volume
!
! Volume of water in modes
REAL, INTENT(IN) :: ukca_water_volume
!
!
! Stratospheric aerosol treated as sulphuric acid?
LOGICAL, INTENT(IN) :: l_in_stratosphere
!

INTEGER, INTENT(IN) :: i_glomap_clim_tune_bc
INTEGER, INTENT(IN) :: i_ukca_tune_bc
!
! When > 0, use a prescribed single scattering albedo field
!
INTEGER, INTENT(IN) :: i_ukca_radaer_prescribe_ssa

REAL, INTENT(OUT) :: re_m
REAL, INTENT(OUT) :: im_m

! Local variables

!
! Modal volume of BC
REAL     :: ukca_modal_bc_vol

!
! BC density tuned plus Maxwell-Garnet mixing method
INTEGER, PARAMETER :: i_ukca_bc_mg_mix        = 2

!
! Complex refractive index after MG mixing
COMPLEX :: refr_mix
!
! Loop integer
INTEGER :: i_cmpt
!
! Index for arrays
INTEGER :: this_cpnt, this_cpnt_type
!
! Switch controlling if MG mixing is required
LOGICAL :: l_mg_mix
!
! Dr Hook is not used to caliper this routine because the overheads
! are too large.

! Initialize refractive index
re_m = 0.0e+00
im_m = 0.0e+00
ukca_modal_bc_vol = 0.0e+00
l_mg_mix = .FALSE.

! If single-scattering albedo is prescribed, then only calculate the
! real refractive index and do not account for MG mixing
IF (i_ukca_radaer_prescribe_ssa /= do_not_prescribe) THEN

  DO i_cmpt = 1, n_cpnt_in_mode(i_mode)

    this_cpnt = i_cpnt_index(i_cmpt, i_mode)

    !
    ! If requested, switch the refractive index of the
    ! sulphate component to that for sulphuric acid
    ! for levels above the tropopause.
    !
    IF ( l_sustrat .AND.                                                       &
         ( i_cpnt_type(this_cpnt) == cp_su ) .AND.                             &
         l_in_stratosphere .AND.                                               &
         ( .NOT. l_nitrate ) ) THEN

      this_cpnt_type = ip_ukca_h2so4

    ELSE

      this_cpnt_type = i_cpnt_type(this_cpnt)

    END IF

    re_m = re_m + ukca_cpnt_volume(this_cpnt) * refr_real(this_cpnt_type)

  END DO ! i_cmpt

  IF (l_soluble(i_mode)) THEN

    !
    ! Account for refractive index of water
    !
    re_m = re_m + ukca_water_volume * refr_real(ip_ukca_water)

  END IF ! l_soluble

ELSE

  DO i_cmpt = 1, n_cpnt_in_mode(i_mode)

    this_cpnt = i_cpnt_index(i_cmpt, i_mode)

    !
    ! If requested, switch the refractive index of the
    ! sulphate component to that for sulphuric acid
    ! for levels above the tropopause.
    !
    IF ( l_sustrat .AND.                                                       &
         ( i_cpnt_type(this_cpnt) == cp_su ) .AND.                             &
         l_in_stratosphere .AND.                                               &
         ( .NOT. l_nitrate ) ) THEN

      this_cpnt_type = ip_ukca_h2so4

    ELSE

      this_cpnt_type = i_cpnt_type(this_cpnt)

    END IF

    !
    ! Work out if Maxwell-Garnett mixing approach will be required
    ! The decision is based on the integer value of
    ! i_ukca_tune_bc or i_ukca_glomap_clim_tune_bc being set to
    ! i_ukca_bc_mg_mix, and this component being BC with non-zero mass
    !
    IF ((i_ukca_tune_bc == i_ukca_bc_mg_mix .OR.                               &
         i_glomap_clim_tune_bc == i_ukca_bc_mg_mix) .AND.                      &
         this_cpnt_type == cp_bc .AND.                                         &
         ukca_cpnt_volume(this_cpnt) > 0.0) THEN

       ! If yes set the BC volume and the logical switch to later use MG mixing
       !
      l_mg_mix = .TRUE.
      ukca_modal_bc_vol = ukca_cpnt_volume(this_cpnt)
    ELSE
      ! If component is not BC, or MG mixing is not requested then
      ! sum up the RI, weighting by component volume
      !
      re_m = re_m + ukca_cpnt_volume(this_cpnt) * refr_real(this_cpnt_type)
      im_m = im_m + ukca_cpnt_volume(this_cpnt) * refr_imag(this_cpnt_type)
    END IF

  END DO ! i_cmpt

  IF ( l_soluble(i_mode) ) THEN

    !
    ! Account for refractive index of water
    !
    re_m = re_m + ukca_water_volume * refr_real(ip_ukca_water)
    im_m = im_m + ukca_water_volume * refr_imag(ip_ukca_water)

  END IF ! l_soluble

  !
  ! Mix in the BC via Maxwell-Garnett?
  ! Only if the BC component type is allowed in this mode
  ! and is present with non-zero volume, and i_ukca_tune_bc==2
  ! or i_glomap_clim_tune_bc==2
  !
  IF (l_mg_mix) THEN

    ! There is the potential for unphysical values or divide by zero errors in
    ! refract_mix_mg if ukca_modal_volume =<  ukca_modal_bc_vol
    ! or the RI assigned to the medium (re_m) <= zero. This would only occur
    ! if the mode contained only pure BC. Therefore, if either those conditions
    ! is true assign re_m and im_m to that of BC and do not call refract_mix_mg
    !
    IF ((ukca_modal_volume  <=  ukca_modal_bc_vol) .OR. (re_m  <=  0.0e+00)) THEN

      re_m = refr_real(cp_bc)
      im_m = refr_imag(cp_bc)

    ELSE

      re_m = re_m / (ukca_modal_volume - ukca_modal_bc_vol)
      im_m = im_m / (ukca_modal_volume - ukca_modal_bc_vol)

      refr_mix = refract_mix_mg(re_m, im_m,                                    &
           refr_real(cp_bc), refr_imag(cp_bc),                                 &
           ukca_modal_volume, ukca_modal_bc_vol)

      re_m = REAL(refr_mix)
      im_m = AIMAG(refr_mix)

    END IF

  ELSE

    re_m = re_m / ukca_modal_volume
    im_m = im_m / ukca_modal_volume

  END IF ! l_mg_mix
END IF ! i_ukca_radaer_prescribe_ssa /= do_not_prescribe

RETURN
END SUBROUTINE ukca_radaer_ri_calc

! ===============================================================
!       MAXWELL-GARNET refractive index mixing function
! ===============================================================

COMPLEX FUNCTION refract_mix_mg(re_m, im_m, re_bc, im_bc,                      &
                        ukca_modal_volume, ukca_modal_bc_vol)

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

! Arguments with intent(in)
REAL, INTENT(IN) :: re_m   ! refract index real part for homogenous medium
REAL, INTENT(IN) :: im_m   ! refract index imag part for homogenous medium
REAL, INTENT(IN) :: re_bc  ! refract index real part for black carbon
REAL, INTENT(IN) :: im_bc  ! refract index imag part for black carbon

!  Scalar values for a specific gridbox and aerosol mode
REAL, INTENT(IN) :: ukca_modal_volume   ! Total wet volume in mode
REAL, INTENT(IN) :: ukca_modal_bc_vol   ! Volume of BC in mode

! Local variables
COMPLEX :: eff_dicnst_m    ! Effective dielectric constant for the medium
COMPLEX :: eff_dicnst_bc   ! Effective dielectric constant for BC
COMPLEX :: eff_dicnst_mix  ! Effective dielectric constant for mixture
COMPLEX :: a               ! Terms in the Maxwell-Garnett expression
REAL    :: vol_frac_bc     ! Volume fraction of BC

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='REFRACT_MIX_MG'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

eff_dicnst_m  = CMPLX(re_m, im_m)**2.0
eff_dicnst_bc = CMPLX(re_bc, im_bc)**2.0
vol_frac_bc   = ukca_modal_bc_vol / ukca_modal_volume

a = vol_frac_bc * (eff_dicnst_bc - eff_dicnst_m)                               &
     / (eff_dicnst_bc + (2.0 * eff_dicnst_m))

eff_dicnst_mix = eff_dicnst_m * (1.0 + (3.0*a)/(1.0-a))
refract_mix_mg    = SQRT(eff_dicnst_mix)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END FUNCTION refract_mix_mg

END MODULE ukca_radaer_ri_calc_mod

