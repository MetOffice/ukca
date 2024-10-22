! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Module containing subroutine for calculating particulate matter diagnostics
!  (PMx where x is upper size threshold in um) and contributions to these from
!  aerosol components.
!  The diagnostics include both PM concentration at ambient conditions and the
!  corresponding dry mass concentrations. All diagnostic concentrations are in
!  ug/m3.
!
! Method:
!
!  Determine the concentrations of the required PMx diagnostics for each mode
!  represented in the model configuration separately as described below and
!  nd sum the concentrations over all modes.
!
!  For each mode, first calculate the mass concentrations for each component
!  present and the total dry and wet mass concentrations (the latter being
!  the concentration at ambient conditions). The PMx mass concentrations are
!  then determined using the cumulative distribution function for the
!  log-normal associated with each mode. For a particular size cut-off
!  d_cutoff (d_cutoff = x), the PM diagnostic is the fraction of the relevant
!  mass concentration in the size range 0 - d_cutoff, given by
!
!    lowfrac = 0.5 * (1 + erf(ln(d_cutoff/dbar)/(sqrt(2)*ln(sigmag)))
!
!  where dbar and sigmag are the geometric mean and standard deviation
!  parameters for the log-normal distribution of particle volume. This is
!  equivalent to the normalized cumulative distribution function given by
!  Eq (8.46) in Seinfeld and Pandis (2016). The particle volume is that at
!  ambient conditions.
!
!  At ambient conditions, the modes's particle number distribution is described
!  by the log-normal with geometric mean diameter wetdp and geometric standard
!  deviation sigmag. The volume distribution parameter dbar is related to the
!  number distribution parameter wetdp by
!
!    dbar = wetdp * exp (3 * (ln(sigmag))^2)
!
!  (Eq (8.53) in Seinfeld and Pandis, 2016).
!
!  Reference: John H. Seinfeld and Spyros N. Pandis. Atmospheric Chemistry
!  and Physics: From Air Pollution to Climate Change, Third Edition.
!  John Wiley & Sons, 2016, 1152pp.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_pm_diags_mod

USE ukca_config_specification_mod, ONLY: glomap_variables

USE ukca_mode_setup,               ONLY: nmodes

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_PM_DIAGS_MOD'

! Structure to hold request flags indicating which PM diagnostics are
! to be calculated
TYPE, PUBLIC :: pm_request_struct
  LOGICAL, ALLOCATABLE :: l_total_dry(:)   ! Total PM dry mass by size category
  LOGICAL, ALLOCATABLE :: l_total_wet(:)   ! Total PM wet mass by size category
  LOGICAL, ALLOCATABLE :: l_component(:,:) ! Component contributions to PM by
                                           ! component number and size category
END TYPE pm_request_struct

PUBLIC :: ukca_pm_diags

CONTAINS


SUBROUTINE ukca_pm_diags(nbox,nd,md,mdwat_diag,wetdp_diag,                     &
                         d_cutoff,pm_request,pm_dry,pm_wet,pm_component)

USE ukca_config_constants_mod, ONLY: avogadro
USE ukca_constants, ONLY: mmw
USE ukca_um_legacy_mod, ONLY: umErf

USE parkind1,               ONLY: jpim, jprb      ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Subroutine arguments

! Number of elements
INTEGER, INTENT(IN) :: nbox

! Aerosol particle number density for each mode (cm^-3)
REAL, INTENT(IN) :: nd(:,:)

! Average component concentrations of aerosol particle in each mode
! (molecules per particle) by component number
REAL, INTENT(IN) :: md(:,:,:)

! Molecular concentration of water in each mode (molecules per particle)
REAL, INTENT(IN) :: mdwat_diag(:,:)

! Geometric mean wet diameter of particles in each mode (m)
REAL, INTENT(IN) :: wetdp_diag(:,:)

! Size limits for particulate matter - upper threshold (m)
REAL, INTENT(IN) :: d_cutoff(:)

! Request flags for controlling PM diagnostic calculation
TYPE(pm_request_struct), INTENT(IN) :: pm_request

! PM diagnostics (ug m^-3)
REAL, INTENT(OUT) :: pm_dry(:,:)          ! Total PM dry mass by size category
REAL, INTENT(OUT) :: pm_wet(:,:)          ! Total PM wet mass by size category
REAL, INTENT(OUT) :: pm_component(:,:,:)  ! Component contributions to PM by
                                          ! component number and size category

! Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: mm(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: sigmag(:)

INTEGER :: i_size_cat       ! loop counter for PM size category
INTEGER :: imode            ! loop counter for modes
INTEGER :: icp              ! loop counter for components

REAL, PARAMETER :: kgpcm3_to_ugpm3 = 1.0e15
                            ! Conversion from kg cm^-3 to ug m^-3

REAL :: mass_contrib(nbox,glomap_variables%ncp)
                            ! Mass contributions from each component in mode

REAL :: mass_dry(nbox)      ! Sum of mass contributions for mode excluding H2O
REAL :: mass_wet(nbox)      ! Sum of mass contributions for mode including H20
REAL :: dbar(nbox)          ! Geometric mean diameter of particles in mode
REAL :: erf_arg(nbox)       ! Argument vector for error function
REAL :: lowfrac(nbox)       ! Fraction of mass below size cutoff

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PM_DIAGS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mm          => glomap_variables%mm
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
sigmag      => glomap_variables%sigmag

! Initialise PM arrays
pm_dry(:,:) = 0.0
pm_wet(:,:) = 0.0
pm_component(:,:,:) = 0.0

! Add PM contributions from all valid modes

DO imode = 1,nmodes
  IF (mode(imode)) THEN

    ! Determine mass concentration in this mode for each valid component and
    ! sum over components to get dry mass for mode
    mass_dry(:) = 0.0
    DO icp = 1,ncp
      IF (component(imode,icp) .AND.                                           &
          (ANY(pm_request%l_total_dry(:)) .OR.                                 &
           ANY(pm_request%l_total_wet(:)) .OR.                                 &
           ANY(pm_request%l_component(icp,:)))) THEN
        mass_contrib(:,icp) =                                                  &
          kgpcm3_to_ugpm3 * nd(:,imode) * md(:,imode,icp) * mm(icp) / avogadro
        mass_dry(:) = mass_dry(:) + mass_contrib(:,icp)
      END IF
    END DO

    ! Determine wet mass concentration for mode
    IF (ANY(pm_request%l_total_wet(:)))                                        &
      mass_wet(:) = mass_dry(:) +                                              &
        kgpcm3_to_ugpm3 * nd(:,imode) * mdwat_diag(:,imode) * mmw / avogadro

    ! Determine volume geometric mean diameters for the mode to define
    ! log-normal cumulative distribution function for the volume/mass
    ! distributions
    dbar(:) = wetdp_diag(:,imode) *                                            &
      EXP(3.0 * LOG(sigmag(imode)) * LOG(sigmag(imode)))

    ! Do the PM calculations for each particle size cutoff as required
    DO i_size_cat = 1,SIZE(d_cutoff)
      IF (pm_request%l_total_dry(i_size_cat) .OR.                              &
          pm_request%l_total_wet(i_size_cat) .OR.                              &
          ANY(pm_request%l_component(:,i_size_cat))) THEN

        ! Determine fraction of mass below size cutoff using c.d.f.
        erf_arg(:) = LOG(d_cutoff(i_size_cat)/dbar(:)) /                       &
                     (SQRT(2.0)*LOG(sigmag(imode)))
        lowfrac(:) = 0.5 * (1.0 + umErf(erf_arg(:)))

        ! Add PM contribution of dry mass conc. in this mode to total dry PM
        IF (pm_request%l_total_dry(i_size_cat)) THEN
          pm_dry(:,i_size_cat) = pm_dry(:,i_size_cat) +                        &
            lowfrac(:) * mass_dry(:)
        END IF

        ! Add PM contribution of wet mass conc. in this mode to total wet PM
        IF (pm_request%l_total_wet(i_size_cat)) THEN
          pm_wet(:,i_size_cat) = pm_wet(:,i_size_cat) +                        &
            lowfrac(:) * mass_wet(:)
        END IF

        ! Add PM contributions of components in this mode to component PMs
        IF (ANY(pm_request%l_component(:,i_size_cat))) THEN
          DO icp = 1,SIZE(pm_component,2)
            IF (component(imode,icp) .AND.                                     &
                pm_request%l_component(icp,i_size_cat)) THEN
              pm_component(:,icp,i_size_cat) =                                 &
                pm_component(:,icp,i_size_cat) +                               &
                lowfrac(:) * mass_contrib(:,icp)
            END IF
          END DO
        END IF

      END IF
    END DO

  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE ukca_pm_diags

END MODULE ukca_pm_diags_mod
