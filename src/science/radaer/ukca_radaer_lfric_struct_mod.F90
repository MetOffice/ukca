! -----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! -----------------------------------------------------------------------------
!
! Purpose:
! Defines maximum dimensions
! Defines type ukca_radaer_lfric_struct, the structure used by UKCA_RADAER
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Contained subroutines:
!      allocate_radaer_lfric_struct
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE ukca_radaer_lfric_struct_mod

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='UKCA_RADAER_LFRIC_STRUCT_MOD'

INTEGER, SAVE :: npd_ukca_cpnt       ! nmodes*ncp

INTEGER, SAVE :: ncp_max_x_nmodes    ! nmodes * ncp_max

! Thresholds on the modal mass-mixing ratio and modal number
! concentrations above which aerosol optical properties are to be
! computed. Placed here as used by multiple subroutines.
REAL, PARAMETER :: threshold_mmr = 1.0e-12 ! kg/kg
! Corresponds to burden of 0.01 mg/m2 if mmr=1.e-12 everywhere.

REAL, PARAMETER :: threshold_vol = 1.0e-25 ! m3/m3
! Corresponds to particle diameter < 10nm

REAL, PARAMETER :: threshold_nbr = 1.0e+00 ! m-3
! A single coarse-mode particle with d=10um per would give
! a mixing ratio of ~1.e-12kg/kg in the lower troposphere

! Main structure holding all the variables needed for
! interacting UKCA aerosols with radiation.

TYPE :: ukca_radaer_lfric_struct

  ! Information about UKCA aerosol modes
  ! ====================================

  ! Actual number of modes, with a default value for
  ! minimising array dimensions.
  INTEGER :: n_mode = 1

  ! Type of mode (i.e. nucleation, Aitken, accum, or coarse)
  INTEGER, ALLOCATABLE :: i_mode_type(:)

  ! Solubility of mode (soluble if true)
  LOGICAL, ALLOCATABLE :: l_soluble(:)

  ! Lower and upper limits on the geometric mean diameter (m)
  ! in each mode
  REAL, ALLOCATABLE :: d0low(:)
  REAL, ALLOCATABLE :: d0up(:)

  ! Geometric standard deviation in this mode
  REAL, ALLOCATABLE :: sigma(:)

  ! Number of components in each mode and index of each
  ! component in array ukca_cpnt_info
  INTEGER, ALLOCATABLE :: n_cpnt_in_mode(:)
  INTEGER, ALLOCATABLE :: i_cpnt_index(:,:)

  ! Modal diameter of the dry aerosol (m)
  REAL,    ALLOCATABLE :: dry_diam(:, :, :, :)

  ! Modal diameter of the wet aerosol (m)
  REAL,    ALLOCATABLE :: wet_diam(:, :, :, :)

  ! Modal densities (kg/m3)
  REAL,    ALLOCATABLE :: modal_rho(:, :, :, :)

  ! Modal volumes (including water for soluble modes)
  REAL,    ALLOCATABLE :: modal_vol(:, :, :, :)

  ! Fractional volume of water in each mode
  REAL,    ALLOCATABLE :: modal_wtv(:, :, :, :)

  ! Modal number concentrations
  REAL,    ALLOCATABLE :: modal_nbr(:, :, :, :)

  ! Information about UKCA aerosol components
  ! =========================================

  ! Actual number of components, with a default value for
  ! minimising array dimensions.
  INTEGER :: n_cpnt = 1

  ! Size of ncp_max
  INTEGER :: ncp_max

  ! Size of ncp
  INTEGER :: ncp

  ! Size of nmodes
  INTEGER :: nmodes

  ! Size of ncp_max_x_nmodes
  INTEGER :: ncp_max_x_nmodes

  ! Type of component (e.g. sulphate)
  INTEGER, ALLOCATABLE :: i_cpnt_type(:)

  ! Mass density of each component (kg/m3)
  REAL,    ALLOCATABLE :: density(:)

  ! Array index of the mode this component belongs to
  INTEGER, ALLOCATABLE :: i_mode(:)

  ! Array index for translating component and mode
  ! to ukca_radaer index
  INTEGER, ALLOCATABLE :: idx_cpnt_mode(:,:)

   ! Component mass-mixing ratio (kg/kg)
  REAL,    ALLOCATABLE :: mix_ratio(:, :, :, :)

  ! Component volumes
  REAL,    ALLOCATABLE :: comp_vol(:, :, :, :)

  ! Switch: if true, use sulphuric acid optical properties in the
  ! stratosphere, instead of ammonium sulphate (if false).
  ! Has the same default as run_ukca:l_ukca_radaer_sustrat or
  ! run_glomap_aeroclim:l_glomap_clim_radaer_sustrat
  ! from which is it assigned.
  LOGICAL :: l_sustrat = .FALSE.

  ! Switch: if true, nitrate scheme is on and sodium nitrate
  ! refractive indices are included in the pcalc file
  LOGICAL :: l_nitrate = .FALSE.

  ! Switch: if true, use the narrow coarse mode LUTs for the
  ! coarse insoluble mode. This is the case if the super coarse
  ! insoluble mode is used
  LOGICAL :: l_cornarrow_ins = .FALSE.

END TYPE ukca_radaer_lfric_struct

PUBLIC :: allocate_radaer_lfric_struct

CONTAINS

! #############################################################################

SUBROUTINE allocate_radaer_lfric_struct(ukca_radaer, glomap_variables)

! To allocate arrays in the structure ukca_radaer using the number of modes
! and number of components configured in the GLOMAP setup routine.

USE ukca_mode_setup, ONLY: nmodes, ncp_max, cp_no3, glomap_variables_type
USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook
IMPLICIT NONE

TYPE(ukca_radaer_lfric_struct), INTENT(IN OUT) :: ukca_radaer
TYPE(glomap_variables_type),    INTENT(IN)     :: glomap_variables

INTEGER :: ncp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_RADAER_LFRIC_STRUCT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


ncp                          = glomap_variables%ncp
npd_ukca_cpnt                = ncp * nmodes
ncp_max_x_nmodes             = ncp_max * nmodes

ukca_radaer%ncp              = ncp
ukca_radaer%ncp_max          = ncp_max
ukca_radaer%nmodes           = nmodes
ukca_radaer%ncp_max_x_nmodes = ncp_max_x_nmodes

! Amend the indices of H2SO4 and water for sodium nitrate inclusion
IF (ncp >= cp_no3) THEN
  IF (ANY(glomap_variables%component (:,cp_no3))) THEN
    ukca_radaer%l_nitrate = .TRUE.
  END IF
END IF

IF (.NOT. ALLOCATED(ukca_radaer%i_mode_type))                                  &
      ALLOCATE(ukca_radaer%i_mode_type(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%l_soluble))                                    &
      ALLOCATE(ukca_radaer%l_soluble(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d0low))                                        &
      ALLOCATE(ukca_radaer%d0low(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d0up))                                         &
      ALLOCATE(ukca_radaer%d0up(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%sigma))                                        &
      ALLOCATE(ukca_radaer%sigma(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%n_cpnt_in_mode))                               &
      ALLOCATE(ukca_radaer%n_cpnt_in_mode(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%i_cpnt_index))                                 &
      ALLOCATE(ukca_radaer%i_cpnt_index(ncp_max,nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%i_cpnt_type))                                  &
      ALLOCATE(ukca_radaer%i_cpnt_type(ncp_max_x_nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%density))                                      &
      ALLOCATE(ukca_radaer%density(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%i_mode))                                       &
      ALLOCATE(ukca_radaer%i_mode(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%idx_cpnt_mode))                                &
      ALLOCATE(ukca_radaer%idx_cpnt_mode(nmodes,ncp_max))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE allocate_radaer_lfric_struct

END MODULE ukca_radaer_lfric_struct_mod
