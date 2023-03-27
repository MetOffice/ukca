! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Calculate emissions of primary marine organic carbon (PMOC) using the
!    Gantt et al., 2012 and 2015 parameterization. Emissions are calculated
!    with respect to chlorophyll-a concentration and wind speed, and
!    reference to the just-calculated sea salt flux.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk. This routine was provided by CSIRO.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
! Subroutine Interface:
MODULE ukca_prim_moc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_PRIM_MOC_MOD'

CONTAINS

SUBROUTINE ukca_prim_moc(row_length, rows, model_levels,                       &
              aer_emmas, aer_emnum, chla, u_scalar_10m,                        &
              mass_pmoc)

!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculate emissions of primary marine organic carbon (PMOC) using the
! Gantt et al., 2012 and 2015 parameterization. Emissions are calculated
! with respect to chlorophyll-a concentration and wind speed, and
! reference to the just-calculated sea salt flux.
!
! Inputs
! ------
! row_length    : No of columns
! rows          : No of rows
! model levels  : No of model levels
! aer_emmas     : sea-salt mass emissions by mode, kg/m2/s
! aer_emnum     : sea-salt number emissions by mode, equiv-kg/m2/s
! chla          : chlorophyll-a concentration (mg/m3)
! u_scalar_10m  : Scalar wind at 10m (ms-1)
!
! Outputs
! -------
! mass_pmoc     : mass flux of primary marine organic carbon, kg(POM) / m2 / s
!
!
! Local Variables
! ---------------
! imode      : index for looping over modes
! mm_da      : Molar mass of dry air
! ssvol      : Volume of sea salt emissions into a given mode (um3/m2/s)
! ssnum      : Number of sea salt particles into a given mode (equiv-kg/m2/s)
! ss_dp      : Sea salt particle dry diameter as in Gantt et al., 2012 (um)
! factor     : Converts from molecls or partcls per m2/s to kg/m2/s
! gantt1     : Fragment of Gantt et al., 2012 Eq.1 broken down for clarity
! gantt2     : Fragment of Gantt et al., 2012 Eq.1 broken down for clarity
! frac_om_ssa : Organic mass fraction of sea-spray aerosol
! dens_ssa   : Apparent density of sea-spray aerosol (g/cm3)
!
! Inputted by module UKCA_UM_LEGACY_MOD, CHEMISTRY_CONSTANTS
! ----------------------------------------------------------
! avogadro   : Avogadro's constant (molecules per mole)
! boltzmann  : Boltzmann constant
! rgas       : Specific gas constant for dry air
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! glomap_variables%mode      : logical for mode on/off
! nmodes    : Number of possible aerosol modes
! glomap_variables%sigmag    : Geometric standard deviation of modes
! glomap_variables%rhocomp   : Mass density of each of the aerosol components
!                              (kgm^-3)
! cp_cl     : index of cpt in which sea-salt mass is stored
! cp_oc     : index of cpt in which organic carbon mass is stored
!
!
! References
! ----------
! Gantt et al., 2012. ACP 12: 8553-8566
! Gantt et al., 2015. GMD 8: 619-629
!
!--------------------------------------------------------------------

USE ukca_mode_setup,                ONLY: nmodes, cp_cl, cp_oc

USE ukca_um_legacy_mod,             ONLY: pi, rgas => r
USE chemistry_constants_mod,        ONLY: avogadro, boltzmann

USE ukca_config_specification_mod,  ONLY: glomap_config, glomap_variables

USE yomhook,                        ONLY: lhook, dr_hook
USE parkind1,                       ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: row_length     ! No of columns
INTEGER, INTENT(IN) :: rows           ! No of rows
INTEGER, INTENT(IN) :: model_levels   ! No of model levels

REAL, INTENT(IN) :: aer_emmas (row_length, rows, model_levels, nmodes)
    ! sea-salt mass emissions by mode, kg/m2/s
REAL, INTENT(IN) :: aer_emnum (row_length, rows, model_levels, nmodes)
    ! sea-salt number emissions by mode, equiv-kg/m2/s
REAL, INTENT(IN) :: chla(row_length,rows)
    ! chlorophyll-a concentration, mg/m-3
REAL, INTENT(IN) :: u_scalar_10m(row_length,rows)
    ! Scalar 10m wind

REAL, INTENT(OUT) :: mass_pmoc(row_length,rows,1)
    ! mass flux of PMOC, returned as kg(POM)/m2/s

! Local variables

INTEGER :: imode
REAL :: mm_da
    ! Molar mass of dry air (kg/mol) = avogadro*boltzmann/rgas
REAL :: lsigmag
    ! log of sigmag for given mode
REAL :: ssvol(row_length,rows)
    ! Volume of sea salt emitted into a given mode (um3/m2/s)
REAL :: ssnum(row_length,rows)
    ! Number of sea salt particles emitted into a given mode (equiv-kg/m2/s)
REAL :: ss_dp(row_length,rows)
    ! Sea salt particle dry diameter as in Gantt et al., 2012, 2015 (um)
REAL :: factor
    ! converts from molecules or particles per m2/s to kg/m2/s
REAL :: gantt1(row_length,rows)
    ! Fragment of Gantt et al., 2012 Eq.1 broken down for clarity
REAL :: gantt2(row_length,rows)
    ! Fragment of Gantt et al., 2012 Eq.1 broken down for clarity
REAL :: frac_om_ssa(row_length,rows)
    ! Organic mass fraction of sea-spray aerosol
REAL :: dens_ssa(row_length,rows)
    ! Apparent density of sea-spray aerosol (g/cm3)
REAL :: scale_emiss
    ! Scaling factor for the marine POM emission

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PRIM_MOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Molar mass of dry air (kg/mol)
mm_da = avogadro*boltzmann/rgas

factor=mm_da/avogadro ! converts from molecls or partcls per m2/s to kg/m2/s

! Determine the scaling factor for emission of marine POM
IF (glomap_config%l_ukca_scale_marine_pom_ems) THEN
  scale_emiss = glomap_config%marine_pom_ems_scaling
ELSE
  ! No scaling.
  scale_emiss = 1.0
END IF

mass_pmoc(:,:,:) = 0.0 ! accumulated over modes, initialise here
DO imode=1,nmodes

  IF (glomap_variables%mode(imode)) THEN

    lsigmag = LOG(glomap_variables%sigmag(imode))

    WHERE ( aer_emmas(:,:,1,imode) > 0.0) ! where sea salt is being emitted

      ! Calculate sea salt particle dry diameter
      ssvol(:,:) = 1e18 * aer_emmas(:,:,1,imode) /                             &
                   glomap_variables%rhocomp(cp_cl) ! um3/m2/s
      ssnum(:,:) = aer_emnum(:,:,1,imode) ! equiv-kg / m2 / s
      ss_dp(:,:) = (( ( factor * ssvol(:,:) / ssnum(:,:) ) / (pi/6.0) )        &
                         / (EXP(4.5*lsigmag*lsigmag))) ** (1.0/3.0) ! um

      ! Gantt et al, 2012 OM_SSA eqn split into separate calculations for
      ! clarity
      gantt1(:,:) = 1.0 + EXP( (3.0*(-2.63*chla(:,:)))                         &
                       + (3.0*0.18*u_scalar_10m(:,:) ) )
      gantt2(:,:) = 1.0 + (0.03 * EXP(6.81 * ss_dp(:,:)))
      frac_om_ssa(:,:) = ((1.0 / gantt1(:,:)) / gantt2(:,:))                   &
                         + (0.03 / gantt1(:,:))

      ! Calculate apparent density of emitted aerosol
      dens_ssa(:,:) = ((frac_om_ssa(:,:) * glomap_variables%rhocomp(cp_oc)) +  &
                      ((1.0 - frac_om_ssa(:,:)) *                              &
                       glomap_variables%rhocomp(cp_cl))) *1e-3 ! g.cm-3

      ! Total emission flux of marine POM
      ! Parameterization is derived as POM fraction, so no need for OC*1.4
      ! Note that Gantt et al., 2012 include additional factor of six
      ! .. scaling to match obs - not included here
      ! *1e-12 to convert volume flux from um3/m2/s to cm3/m2/s
      ! Final units for mass_pmoc is kg(POM)/m2/s.
      mass_pmoc(:,:,1) = mass_pmoc(:,:,1) +                                    &
                         (scale_emiss * (ssvol(:,:)*1e-12) *                   &
                          frac_om_ssa(:,:) * dens_ssa(:,:) * 1e-3)

    END WHERE ! where sea salt is emitted

  END IF  ! if glomap_variables%mode is defined

END DO ! loop over nmodes

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE ukca_prim_moc

END MODULE ukca_prim_moc_mod
