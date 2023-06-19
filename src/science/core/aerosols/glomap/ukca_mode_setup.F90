! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module to store MODE setup arrays
!    Contains public subroutines:
!      UKCA_MODE_ALLCP_4MODE
!      UKCA_MODE_SUSS_4MODE
!      UKCA_MODE_SUSSBCOC_4MODE
!      UKCA_MODE_SUSSBCOC_5MODE
!      UKCA_MODE_SUSSBCOCSO_5MODE
!      UKCA_MODE_SUSSBCOCSO_4MODE
!      UKCA_MODE_DUonly_2MODE
!      UKCA_MODE_DUonly_3MODE (needs to be added at some point)
!      UKCA_MODE_SUSSBCOCDU_7MODE
!      UKCA_MODE_SUSSBCOCDUNTNH_8MODE_8CPT
!      UKCA_MODE_SUSSBCOCDU_4MODE
!      UKCA_MODE_SUSSBCOCNTNH_5MODE_7CPT
!      UKCA_MODE_SOLINSOL_6MODE
!    which define modes and components for different components/modes setup.
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
!    Language:  Fortran
!
! ######################################################################
!
! Subroutine Interface:
MODULE ukca_mode_setup
! ---------------------------------------------------------------------|
!  Module to contain modes and components
!
! Description:
! To allow use throughout UM, module stores MODE setup arrays
!
! Note: currently code is hard-coded so that ordering of modes must
! 1) nucln, 2)   soluble Aitken, 3)   soluble accum, 4)   soluble coarse
!           5) insoluble Aitken, 6) insoluble accum, 7) insoluble coarse
!
! ---------------------------------------------------------------------|

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umPrint, umMessage
USE chemistry_constants_mod, ONLY: avogadro, rho_so4

IMPLICIT NONE

PUBLIC
SAVE

INTEGER, PARAMETER :: nmodes=8       ! No of modes
INTEGER, PARAMETER :: nmodes_sol = 4 ! No of soluble modes
INTEGER, PARAMETER :: nmodes_ins = 4 ! No of insoluble modes
INTEGER, PARAMETER :: ncp_max=9      ! No of components
INTEGER, PARAMETER :: ncation=3      ! No possible cation species
INTEGER, PARAMETER :: nanion =4      ! No possible anion species

INTEGER, PARAMETER :: cp_su=1   ! Index to store SO4    cpt
INTEGER, PARAMETER :: cp_bc=2   ! Index to store BC     cpt
INTEGER, PARAMETER :: cp_oc=3   ! Index to store 1st OC cpt
INTEGER, PARAMETER :: cp_cl=4   ! Index to store NaCl   cpt
INTEGER, PARAMETER :: cp_du=5   ! Index to store dust   cpt
INTEGER, PARAMETER :: cp_so=6   ! Index to store 2nd OC cpt
INTEGER, PARAMETER :: cp_no3=7  ! Index to store NO3    cpt
INTEGER, PARAMETER :: cp_nn=8   ! Index to store NaNO3  cpt
INTEGER, PARAMETER :: cp_nh4=9  ! Index to store NH4    cpt

INTEGER, PARAMETER :: mode_nuc_sol  =1 ! Index of nucleation sol mode
INTEGER, PARAMETER :: mode_ait_sol  =2 ! Index of Aitken sol mode
INTEGER, PARAMETER :: mode_acc_sol  =3 ! Index of accumulation sol
INTEGER, PARAMETER :: mode_cor_sol  =4 ! Index of coarse sol mode
INTEGER, PARAMETER :: mode_ait_insol=5 ! Index of Aitken insol mode
INTEGER, PARAMETER :: mode_acc_insol=6 ! Index of accumulation insol
INTEGER, PARAMETER :: mode_cor_insol=7 ! Index of coarse insol mode
INTEGER, PARAMETER :: mode_sup_insol=8 ! Index of super coarse insol mode

INTEGER, PARAMETER :: moment_number =0 ! Index of number moment
INTEGER, PARAMETER :: moment_mass   =3 ! Index of mass moment

! Internal IDs for the types of UKCA aerosol modes used by RADAER.
INTEGER, PARAMETER :: ip_ukca_mode_nucleation  = 0
INTEGER, PARAMETER :: ip_ukca_mode_aitken      = 1
INTEGER, PARAMETER :: ip_ukca_mode_accum       = 2
INTEGER, PARAMETER :: ip_ukca_mode_coarse      = 3
INTEGER, PARAMETER :: ip_ukca_mode_supercoarse = 4

INTEGER, PARAMETER :: i_ukca_bc_tuned         = 1 ! BC density tuned
INTEGER, PARAMETER :: i_ukca_bc_mg_mix        = 2 ! BC density tuned, plus
                                                  ! Maxwell-Garnet mixing method

! =============================================================================
! ukca_mode_sussbcoc_5mode specific settings
!
! This is to allow DO loops over modes and components without IF statements

! Number of modes used by sussbcoc_5mode
INTEGER, PARAMETER :: nmodes_list_sussbcoc_5mode = 5
! List of modes used by sussbcoc_5mode
INTEGER :: mode_list_sussbcoc_5mode(nmodes_list_sussbcoc_5mode) =              &
   [ mode_nuc_sol, mode_ait_sol, mode_acc_sol, mode_cor_sol, mode_ait_insol ]

! Number of components used by sussbcoc_5mode
INTEGER, PARAMETER :: ncp_list_sussbcoc_5mode = 15
! List of components used by sussbcoc_5mode (by mode)
INTEGER :: component_list_by_mode_sussbcoc_5mode(ncp_list_sussbcoc_5mode) =    &
            [ mode_nuc_sol,   mode_nuc_sol,                                    &
              mode_ait_sol,   mode_ait_sol,   mode_ait_sol,                    &
              mode_acc_sol,   mode_acc_sol,   mode_acc_sol,   mode_acc_sol,    &
              mode_cor_sol,   mode_cor_sol,   mode_cor_sol,   mode_cor_sol,    &
              mode_ait_insol, mode_ait_insol ]
! List of components used by sussbcoc_5mode (by component)
INTEGER :: component_list_by_cp_sussbcoc_5mode(ncp_list_sussbcoc_5mode) =      &
            [ cp_su,          cp_oc,                                           &
              cp_su,          cp_bc,          cp_oc,                           &
              cp_su,          cp_bc,          cp_oc,          cp_cl,           &
              cp_su,          cp_bc,          cp_oc,          cp_cl,           &
              cp_bc,          cp_oc ]

! =============================================================================
! ukca_mode_sussbcocdu_7mode specific settings
!
! This is to allow DO loops over modes and components without IF statements

! Number of modes used by sussbcocdu_7mode
INTEGER, PARAMETER :: nmodes_list_sussbcocdu_7mode = 7
! List of modes used by sussbcocdu_7mode
INTEGER :: mode_list_sussbcocdu_7mode(nmodes_list_sussbcocdu_7mode) =          &
    [ mode_nuc_sol, mode_ait_sol,   mode_acc_sol,   mode_cor_sol,              &
                    mode_ait_insol, mode_acc_insol, mode_cor_insol ]

! Number of components used by sussbcocdu_7mode
INTEGER, PARAMETER :: ncp_list_sussbcocdu_7mode = 19
! List of components used by sussbcocdu_7mode (by mode)
INTEGER :: component_list_by_mode_sussbcocdu_7mode(ncp_list_sussbcocdu_7mode)= &
   [ mode_nuc_sol,   mode_nuc_sol,                                             &
     mode_ait_sol,   mode_ait_sol,  mode_ait_sol,                              &
     mode_acc_sol,   mode_acc_sol,  mode_acc_sol, mode_acc_sol, mode_acc_sol,  &
     mode_cor_sol,   mode_cor_sol,  mode_cor_sol, mode_cor_sol, mode_cor_sol,  &
     mode_ait_insol, mode_ait_insol,                                           &
     mode_acc_insol,                                                           &
     mode_cor_insol ]
! List of components used by sussbcocdu_7mode (by component)
INTEGER :: component_list_by_cp_sussbcocdu_7mode(ncp_list_sussbcocdu_7mode) =  &
  [  cp_su,          cp_oc,                                                    &
     cp_su,          cp_bc,         cp_oc,                                     &
     cp_su,          cp_bc,         cp_oc,        cp_cl,        cp_du,         &
     cp_su,          cp_bc,         cp_oc,        cp_cl,        cp_du,         &
     cp_bc,          cp_oc,                                                    &
     cp_du,                                                                    &
     cp_du ]

! =============================================================================
! mode coagulation table
! Modes resulting when two modes coagulate...

INTEGER, PARAMETER :: coag_mode( nmodes, nmodes ) = RESHAPE( [                 &
                                                      1, 2, 3, 4, 2, 3, 4, 4,  &
                                                      2, 2, 3, 4, 2, 3, 4, 4,  &
                                                      3, 3, 3, 4, 3, 3, 4, 4,  &
                                                      4, 4, 4, 4, 4, 4, 4, 4,  &
                                                      2, 2, 3, 4, 5, 6, 7, 8,  &
                                                      3, 3, 3, 4, 6, 6, 7, 8,  &
                                                      4, 4, 4, 4, 7, 7, 7, 8,  &
                                                      4, 4, 4, 4, 8, 8, 8, 8], &
                                                           [ nmodes, nmodes ] )

! =============================================================================

REAL, PARAMETER :: rho_nacl = 2165.0       ! Correct NaCl density (kg m^-3)
!
! The default value for BC set in rhocomp is 1500.0 Kg m^-3 but two
! alternative estimates can be used as tuning options to reduce BC absorption
! efficiency
!
! Estimate for BC density within the range given by Bond and Bergstrom (2006)
! tuned for use with the Maxwell-Garnet mixing approximation in RADAER
REAL, PARAMETER :: rho_bc_mg_mix = 1800.0  ! (Kg m^-3)

! High estimate for BC density based on Bond and Bergstrom (2006)
! tuned for use with the volume-mixing approximation in RADAER
REAL, PARAMETER :: rho_bc_tuned = 1900.0  ! (Kg m^-3)

! Mode names - the same for all MODE set ups, different ones are
! on for different set ups - see mode_choice
! Set case to be consistent with names in ukca_fieldname_mod for easier pattern
! matching
CHARACTER(LEN=7), PARAMETER :: mode_names(nmodes) =                            &
    ['Nuc_SOL','Ait_SOL','Acc_SOL','Cor_SOL',                                  &
               'Ait_INS','Acc_INS','Cor_INS','Sup_INS']

! =============================================================================
! -- Type for holding GLOMAP_mode variables --

TYPE :: glomap_variables_type
  ! No. of components
  INTEGER :: ncp

  ! Mode switches (1=on, 0=0ff)
  INTEGER :: mode_choice(nmodes)

  ! Specify which modes are soluble
  INTEGER :: modesol(nmodes)

  ! Fraction of bc ems to go into each mode
  REAL :: fracbcem(nmodes)

  ! Fraction of om ems to go into each mode
  REAL :: fracocem(nmodes)

  ! Lower size limits of geometric mean diameter for each mode
  REAL :: ddplim0(nmodes)

  ! Mid-point of size mode (m)
  REAL :: ddpmid(nmodes)

  ! Upper size limits of geometric mean diameter for each mode
  REAL :: ddplim1(nmodes)

  ! Mid-point masses for initial radius grid
  REAL :: mmid(nmodes)

  ! Lo-interf masses for initial radius grid
  REAL :: mlo(nmodes)

  ! Hi-interf masses for initial radius grid
  REAL :: mhi(nmodes)

  ! Threshold for number in mode to carry out calculations
  REAL :: num_eps(nmodes)

  ! Fixed geometric standard deviation for each mode
  REAL :: sigmag(nmodes)

  ! Modes (T/F)
  ! Note this must be declared as .FALSE.
  ! otherwise it will be unitialised when required by for example AQUM suites
  LOGICAL :: mode(nmodes) = .FALSE.

  ! Component switches (1=on, 0=off)
  INTEGER, ALLOCATABLE :: component_choice(:)

  ! Components that are soluble
  INTEGER, ALLOCATABLE :: soluble_choice(:)

  ! Components allowed in each mode (must be consistent with coag_mode)
  INTEGER, ALLOCATABLE :: component_mode(:,:)

  ! Initial fractions of mass in each mode among components
  REAL, ALLOCATABLE :: mfrac_0(:,:)

  ! Molar masses of components (kg mol-1)
  REAL, ALLOCATABLE :: mm(:)

  ! Mass density of components (kg m^-3)
  REAL, ALLOCATABLE :: rhocomp(:)

  ! Number of dissociating ions in soluble components
  REAL, ALLOCATABLE :: no_ions(:)

  ! Component (T/F)
  LOGICAL, ALLOCATABLE :: component(:,:)

  ! Components that are soluble
  LOGICAL, ALLOCATABLE :: soluble(:)

  ! Component names
  CHARACTER(LEN=7), ALLOCATABLE :: component_names(:)

  ! EXP((9/2)*LOG^2(SIGMA_G))
  REAL :: x(nmodes)

END TYPE glomap_variables_type


! =============================================================================

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_MODE_SETUP'

CONTAINS


! Subroutine Interface:

SUBROUTINE ukca_mode_allcp_4mode ( glomap_variables_local,                     &
                                   l_radaer_in,                                &
                                   i_tune_bc_in,                               &
                                   l_fix_nacl_density_in,                      &
                                   pi )

! ---------------------------------------------------------------------|
!  Subroutine to define modes and components with all components
!  switched on but only 4 modes used.
! ---------------------------------------------------------------------|
IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp


! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_ALLCP_4MODE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

glomap_variables_local%ncp = 6
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[1,1,1,1,0,0,0,0]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[1,1,1,1,1,0]
! ***n.b. in above have kept all cpts on (not SO) for UM test***
! Components that are soluble
glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
glomap_variables_local%component_mode(1,1:ncp)=[1,0,1,0,0,1] !allowed in nuc_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,1,1,0,0,1] !allowed in ait_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,1,1,1,1,1] !allowed in acc_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,1,1,1,1,1] !allowed in cor_sol
glomap_variables_local%component_mode(5,1:ncp)=[0,1,1,0,0,0] !allowed in ait_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0] !allowed in acc_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0] !allowed in cor_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0] !allowed in sup_ins

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,1.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0,1.8]

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Nuc soluble
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Ait soluble
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Acc soluble
glomap_variables_local%mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !Cor soluble
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Ait insoluble
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Acc insoluble
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Cor insoluble
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Sup insoluble

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!                                 h2so4 bc    oc     nacl    dust  so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
                                    [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO

! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,2.0,0.0,0.0]
!
glomap_variables_local%fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).

! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_allcp_4mode

SUBROUTINE ukca_mode_suss_4mode( glomap_variables_local,                       &
                                 l_radaer_in,                                  &
                                 i_tune_bc_in,                                 &
                                 l_fix_nacl_density_in,                        &
                                 pi )
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate and sea-salt only in 4 modes.
!  Uses 10 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp

! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSS_4MODE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

glomap_variables_local%ncp = 6
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[1,1,1,1,0,0,0,0]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[1,0,0,1,0,0]
! *** n.b. only have h2so4 and nacl cpts on for this setup ***
! Components that are soluble
glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
glomap_variables_local%component_mode(1,1:ncp)=[1,0,1,0,0,1] !allowed in nuc_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,1,1,0,0,1] !allowed in ait_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,1,1,1,1,1] !allowed in acc_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,1,1,1,1,1] !allowed in cor_sol
glomap_variables_local%component_mode(5,1:ncp)=[0,1,1,0,0,0] !allowed in ait_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0] !allowed in acc_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0] !allowed in cor_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0] !allowed in sup_ins

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,1.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0,1.8]

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Nuc soluble
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Ait soluble
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Acc soluble
glomap_variables_local%mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !Cor soluble
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Ait insoluble
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Acc insoluble
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Cor insoluble
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Sup insoluble

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!                                 h2so4 bc    oc     nacl    dust  so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
                                    [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO

! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,2.0,0.0,0.0]
!
glomap_variables_local%fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).

! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_suss_4mode

SUBROUTINE ukca_mode_sussbcocdu_4mode( glomap_variables_local,                 &
                                       l_radaer_in,                            &
                                       i_tune_bc_in,                           &
                                       l_fix_nacl_density_in,                  &
                                       pi )
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  SO4, sea-salt, bc, oc (secondary & primary combined) & du in 4 modes.
!  Uses 19 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp

! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOCDU_4MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

glomap_variables_local%ncp = 6
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[1,1,1,1,0,0,0,0]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[1,1,1,1,1,0]
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
glomap_variables_local%component_mode(1,1:ncp)=[1,0,1,0,0,1] !allowed in nuc_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,1,1,0,0,1] !allowed in ait_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,1,1,1,1,1] !allowed in acc_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,1,1,1,1,1] !allowed in cor_sol
glomap_variables_local%component_mode(5,1:ncp)=[0,1,1,0,0,0] !allowed in ait_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0] !allowed in acc_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0] !allowed in cor_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0] !allowed in sup_ins

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,1.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0,1.8]

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Nuc soluble
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Ait soluble
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Acc soluble
glomap_variables_local%mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !Cor soluble
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Ait insoluble
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Acc insoluble
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Cor insoluble
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Sup insoluble

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!                                 h2so4 bc    oc     nacl    dust  so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
                                    [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO

! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,2.0,0.0,0.0]
!
glomap_variables_local%fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).


! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcocdu_4mode

SUBROUTINE ukca_mode_sussbcocdu_7mode( glomap_variables_local,                 &
                                       l_radaer_in,                            &
                                       i_tune_bc_in,                           &
                                       l_fix_nacl_density_in,                  &
                                       pi )
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  SO4, sea-salt, bc, oc (secondary & primary combined) & du in 7 modes.
!  Uses 26 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp

! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOCDU_7MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

glomap_variables_local%ncp = 6
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[1,1,1,1,1,1,1,0]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[1,1,1,1,1,0]
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
glomap_variables_local%component_mode(1,1:ncp)=[1,0,1,0,0,1] !allowed in nuc_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,1,1,0,0,1] !allowed in ait_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,1,1,1,1,1] !allowed in acc_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,1,1,1,1,1] !allowed in cor_sol
glomap_variables_local%component_mode(5,1:ncp)=[0,1,1,0,0,0] !allowed in ait_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0] !allowed in acc_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0] !allowed in cor_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0] !allowed in sup_ins

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,1.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0,1.8]

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Nuc soluble
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Ait soluble
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Acc soluble
glomap_variables_local%mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !Cor soluble
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Ait insoluble
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Acc insoluble
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Cor insoluble
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Sup insoluble

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!                                 h2so4 bc    oc     nacl    dust  so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
                                    [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO

! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,2.0,0.0,0.0]
!
glomap_variables_local%fracbcem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).

! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcocdu_7mode

SUBROUTINE ukca_mode_sussbcoc_4mode( glomap_variables_local,                   &
                                     l_radaer_in,                              &
                                     i_tune_bc_in,                             &
                                     l_fix_nacl_density_in,                    &
                                     pi )
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate, sea-salt, bc & oc (secondary & primary combined) in 4 modes.
!  Uses 17 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp

! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOC_4MODE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

glomap_variables_local%ncp = 6
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[1,1,1,1,0,0,0,0]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[1,1,1,1,0,0]
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
glomap_variables_local%component_mode(1,1:ncp)=[1,0,1,0,0,1] !allowed in nuc_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,1,1,0,0,1] !allowed in ait_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,1,1,1,1,1] !allowed in acc_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,1,1,1,1,1] !allowed in cor_sol
glomap_variables_local%component_mode(5,1:ncp)=[0,1,1,0,0,0] !allowed in ait_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0] !allowed in acc_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0] !allowed in cor_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0] !allowed in sup_ins

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,1.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0,1.8]

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Nuc soluble
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Ait soluble
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Acc soluble
glomap_variables_local%mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !Cor soluble
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Ait insoluble
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Acc insoluble
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Cor insoluble
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Sup insoluble

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!                                 h2so4 bc    oc     nacl    dust  so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
                                    [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO

! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,2.0,0.0,0.0]
!
glomap_variables_local%fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).

! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcoc_4mode

SUBROUTINE ukca_mode_sussbcoc_5mode( glomap_variables_local,                   &
                                     l_radaer_in,                              &
                                     i_tune_bc_in,                             &
                                     l_fix_nacl_density_in,                    &
                                     pi )
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate, sea-salt, bc & oc (secondary & primary combined) in 5 modes.
!  Uses 20 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp

! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOC_5MODE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

glomap_variables_local%ncp = 6
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[1,1,1,1,1,0,0,0]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[1,1,1,1,0,0]
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
glomap_variables_local%component_mode(1,1:ncp)=[1,0,1,0,0,1] !allowed in nuc_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,1,1,0,0,1] !allowed in ait_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,1,1,1,1,1] !allowed in acc_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,1,1,1,1,1] !allowed in cor_sol
glomap_variables_local%component_mode(5,1:ncp)=[0,1,1,0,0,0] !allowed in ait_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0] !allowed in acc_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0] !allowed in cor_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0] !allowed in sup_ins

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,1.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0,1.8]

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Nuc soluble
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Ait soluble
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Acc soluble
glomap_variables_local%mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !Cor soluble
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Ait insoluble
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Acc insoluble
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Cor insoluble
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Sup insoluble

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!                                 h2so4 bc    oc     nacl    dust  so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
                                    [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO

! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,2.0,0.0,0.0]
!
glomap_variables_local%fracbcem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).

! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcoc_5mode

SUBROUTINE ukca_mode_sussbcocso_4mode( glomap_variables_local,                 &
                                       l_radaer_in,                            &
                                       i_tune_bc_in,                           &
                                       l_fix_nacl_density_in,                  &
                                       pi )
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate, sea-salt, bc, primary oc & secondary oc cpts in 5 modes.
!  Uses 20 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp

! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOCSO_4MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

glomap_variables_local%ncp = 6
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[1,1,1,1,0,0,0,0]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[1,1,1,1,0,1]
! ***all cpts on except dust***
! Components that are soluble
glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
glomap_variables_local%component_mode(1,1:ncp)=[1,0,1,0,0,1] !allowed in nuc_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,1,1,0,0,1] !allowed in ait_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,1,1,1,1,1] !allowed in acc_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,1,1,1,1,1] !allowed in cor_sol
glomap_variables_local%component_mode(5,1:ncp)=[0,1,1,0,0,0] !allowed in ait_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0] !allowed in acc_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0] !allowed in cor_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0] !allowed in sup_ins

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,1.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0,1.8]

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Nuc soluble
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Ait soluble
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Acc soluble
glomap_variables_local%mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !Cor soluble
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Ait insoluble
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Acc insoluble
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Cor insoluble
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Sup insoluble

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!                                 h2so4 bc    oc     nacl    dust  so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
                                    [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO

! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,2.0,0.0,0.0]
!
glomap_variables_local%fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into   soluble Aitken for this setup).

! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcocso_4mode

SUBROUTINE ukca_mode_sussbcocso_5mode( glomap_variables_local,                 &
                                       l_radaer_in,                            &
                                       i_tune_bc_in,                           &
                                       l_fix_nacl_density_in,                  &
                                       pi )
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate, sea-salt, bc, primary oc & secondary oc cpts in 5 modes.
!  Uses 23 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp

! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOCSO_5MODE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

glomap_variables_local%ncp = 6
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[1,1,1,1,1,0,0,0]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[1,1,1,1,0,1]
! ***all cpts on except dust***
! Components that are soluble
glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
glomap_variables_local%component_mode(1,1:ncp)=[1,0,1,0,0,1] !allowed in nuc_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,1,1,0,0,1] !allowed in ait_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,1,1,1,1,1] !allowed in acc_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,1,1,1,1,1] !allowed in cor_sol
glomap_variables_local%component_mode(5,1:ncp)=[0,1,1,0,0,0] !allowed in ait_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0] !allowed in acc_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0] !allowed in cor_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0] !allowed in sup_ins

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,1.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0,1.8]

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Nuc soluble
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Ait soluble
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Acc soluble
glomap_variables_local%mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !Cor soluble
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Ait insoluble
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Acc insoluble
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Cor insoluble
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Sup insoluble

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!                                 h2so4 bc    oc     nacl    dust  so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
                                    [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO

! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,2.0,0.0,0.0]
!
glomap_variables_local%fracbcem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).

! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcocso_5mode

SUBROUTINE ukca_mode_duonly_2mode( glomap_variables_local,                     &
                                   l_radaer_in,                                &
                                   i_tune_bc_in,                               &
                                   l_fix_nacl_density_in,                      &
                                   pi )
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  only du cpt in 2 (insoluble) modes.
!  Uses  4 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp

! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_DUONLY_2MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

glomap_variables_local%ncp = 6
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[0,0,0,0,0,1,1,0]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[0,0,0,0,1,0]! ***only dust on***
! Components that are soluble
glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
glomap_variables_local%component_mode(1,1:ncp)=[1,0,1,0,0,1] !allowed in nuc_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,1,1,0,0,1] !allowed in ait_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,1,1,1,1,1] !allowed in acc_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,1,1,1,1,1] !allowed in cor_sol
glomap_variables_local%component_mode(5,1:ncp)=[0,1,1,0,0,0] !allowed in ait_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0] !allowed in acc_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0] !allowed in cor_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0] !allowed in sup_ins

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,1.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0,1.8]

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Nuc soluble
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Ait soluble
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Acc soluble
glomap_variables_local%mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !Cor soluble
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Ait insoluble
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Acc insoluble
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Cor insoluble
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Sup insoluble

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!                                 h2so4 bc    oc     nacl    dust  so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
                                    [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO

! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,2.0,0.0,0.0]
!
glomap_variables_local%fracbcem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).

! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_duonly_2mode

SUBROUTINE ukca_mode_sussbcocntnh_5mode_7cpt( glomap_variables_local,          &
                                              l_radaer_in,                     &
                                              i_tune_bc_in,                    &
                                              l_fix_nacl_density_in,           &
                                              pi )
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate, sea-salt, bc, oc (secondary & primary combined),
!  NH4, and NO3 in 5 modes, and 7 components.
!  Also adding NaNO3 - a coarse nitrate tracer in ACC/COR
!  soluble modes - 7 components
!  Uses 28 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp

! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOCNTNH_5MODE_7CPT'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! No of components
glomap_variables_local%ncp = 9
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org',          &
          'no3    ','nano3  ','nh4    ']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[1,1,1,1,1,0,0,0]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[1,1,1,1,0,0,1,1,1]
! *** n.b. only have h2so4, bc, oc, nacl, no3, nh4 cpts on for this setup ***
! Components that are soluble
glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,1,0,0,1,1,1]
! Components allowed in each mode (must be consistent with coag_mode)
!allowed nuc_sol
glomap_variables_local%component_mode(1,1:ncp)=[1,0,1,0,0,1,0,0,0]
!allowed ait_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,1,1,0,0,1,1,0,1]
!allowed acc_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,1,1,1,1,1,1,1,1]
!allowed cor_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,1,1,1,1,1,1,1,1]
!allowed ait_ins
glomap_variables_local%component_mode(5,1:ncp)=[0,1,1,0,0,0,0,0,0]
!allowed acc_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0,0,0,0]
!allowed cor_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0,0,0,0]
!allowed in sup_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0,0,0,0]

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,1.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0,1.8]

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
!NucSol
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
!AitSol
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
!AccSol
glomap_variables_local%mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0]
!CorSol
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0]
!AitIns
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0]
!AccIns
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0]
!CorIns
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0]
!SupIns

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=                                              &
            [0.098, 0.012, 0.0168, 0.05844, 0.100, 0.0168, 0.062, 0.084, 0.018]
!            h2so4  bc     oc      nacl     dust   so      no3    nano3  nh4
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
       [1769.0, 1500.0, 1500.0, 1600.0, 2650.0, 1500.0, 1500.0, 1600.0, 1769.0]
!       h2so4   bc      oc      nacl    dust    so      no3     nano3   nh4

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO

! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,2.0,0.0,0.0,2.0,2.0,2.0]

! Fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
glomap_variables_local%fracbcem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]

! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcocntnh_5mode_7cpt

SUBROUTINE  ukca_mode_solinsol_6mode( glomap_variables_local,                  &
                                      l_radaer_in,                             &
                                      i_tune_bc_in,                            &
                                      l_fix_nacl_density_in,                   &
                                      pi )
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  soluble and insoluble components in 6 modes.
!  Uses 12 aerosol tracers
! ---------------------------------------------------------------------|

IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp

! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SOLINSOL_6MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! No of components
glomap_variables_local%ncp = 6
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[1,1,1,1,0,1,1,0]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[1,0,0,0,1,0]
! *** n.b. only have soluble (as h2so4) and insoluble (as dust) in this setup
! Components that are soluble

glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,0,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
!allowed nuc_sol
glomap_variables_local%component_mode(1,1:ncp)=[1,0,0,0,0,0]
!allowed ait_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,0,0,0,0,0]
!allowed acc_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,0,0,0,0,0]
!allowed cor_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,0,0,0,0,0]
!allowed ait_ins
glomap_variables_local%component_mode(5,1:ncp)=[0,0,0,0,0,0]
!allowed acc_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0]
!allowed cor_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0]
!allowed sup_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0]


! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,1.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0,1.8]
! to match M7, but sigacc=1.4

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Nuc soluble
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Ait soluble
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Acc soluble
glomap_variables_local%mfrac_0(4,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Cor soluble
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Ait insoluble
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Acc insoluble
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Cor insoluble
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Sup insoluble

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=                                              &
           [ 0.098, 0.012, 0.0168, 0.05844, 0.100, 0.0168]
!            h2so4  bc     oc      nacl     dust   so
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
                                    [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO


! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,0.0,0.0,0.0]
!
glomap_variables_local%fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).

! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_solinsol_6mode

SUBROUTINE ukca_mode_sussbcocduntnh_8mode_8cpt( glomap_variables_local,        &
                                                l_radaer_in,                   &
                                                i_tune_bc_in,                  &
                                                l_fix_nacl_density_in,         &
                                                pi )
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate, sea-salt, bc, oc, du (secondary & primary combined),
!  NH4, and NO3 in 5 modes, and 7 components.
!  Also adding NaNO3 - a coarse nitrate tracer in ACC/COR
!  soluble modes - 8 components
!  Uses 34 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

! Arguments

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local
LOGICAL,                     INTENT(IN)     :: l_radaer_in
INTEGER,                     INTENT(IN)     :: i_tune_bc_in
LOGICAL,                     INTENT(IN)     :: l_fix_nacl_density_in
REAL,                        INTENT(IN)     :: pi

! Local variables

INTEGER :: imode
INTEGER :: icp
INTEGER :: ncp

! specifies average (rho/mm) for default composition given by mfrac_0
REAL :: rhommav

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER ::                                                 &
  RoutineName='UKCA_MODE_SUSSBCOCDUNTNH_8MODE_8CPT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! No of components
glomap_variables_local%ncp = 9
ncp = glomap_variables_local%ncp

CALL ukca_mode_allocate_ctl_vars ( ncp , glomap_variables_local )

glomap_variables_local%mode(:) = .FALSE.

! Component names
glomap_variables_local%component_names(1:ncp) =                                &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org',          &
          'no3    ','nano3  ','nh4    ']

! Mode switches (1=on, 0=0ff)
glomap_variables_local%mode_choice=[1,1,1,1,1,1,1,1]

! Specify which modes are soluble
glomap_variables_local%modesol=[1,1,1,1,0,0,0,0]

! Component switches (1=on, 0=off)
glomap_variables_local%component_choice(1:ncp)=[1,1,1,1,1,0,1,1,1]
! ** n.b. only have h2so4, bc, oc, du, nacl, no3, nh4 cpts on for this setup **
! Components that are soluble
glomap_variables_local%soluble_choice(1:ncp)=[1,0,0,1,0,0,1,1,1]
! Components allowed in each mode (must be consistent with coag_mode)
!allowed nuc_sol
glomap_variables_local%component_mode(1,1:ncp)=[1,0,1,0,0,1,0,0,0]
!allowed ait_sol
glomap_variables_local%component_mode(2,1:ncp)=[1,1,1,0,0,1,1,0,1]
!allowed acc_sol
glomap_variables_local%component_mode(3,1:ncp)=[1,1,1,1,1,1,1,1,1]
!allowed cor_sol
glomap_variables_local%component_mode(4,1:ncp)=[1,1,1,1,1,1,1,1,1]
!allowed ait_ins
glomap_variables_local%component_mode(5,1:ncp)=[0,1,1,0,0,0,0,0,0]
!allowed acc_ins
glomap_variables_local%component_mode(6,1:ncp)=[0,0,0,0,1,0,0,0,0]
!allowed cor_ins
glomap_variables_local%component_mode(7,1:ncp)=[0,0,0,0,1,0,0,0,0]
!allowed in sup_ins
glomap_variables_local%component_mode(8,1:ncp)=[0,0,0,0,1,0,0,0,0]

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
glomap_variables_local%ddplim0=                                                &
                      [1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6,5.0e-6]
glomap_variables_local%ddplim1=                                                &
                      [1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,5.0e-6,5.0e-5]

! Specify fixed geometric standard deviation for each mode
glomap_variables_local%sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,1.59,1.8]

DO imode=1,nmodes
  glomap_variables_local%x(imode)=EXP(4.5 *                                    &
                                 LOG( glomap_variables_local%sigmag(imode) ) * &
                                 LOG( glomap_variables_local%sigmag(imode) ) )
END DO

! Specify threshold for ND (per cc) below which don't do calculations
glomap_variables_local%num_eps =                                               &
                  [1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14,1.0e-20]

! Initial fractions of mass in each mode among components
glomap_variables_local%mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
!NucSol
glomap_variables_local%mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
!AitSol
glomap_variables_local%mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
!AccSol
glomap_variables_local%mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0]
!CorSol
glomap_variables_local%mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0]
!AitIns
glomap_variables_local%mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0]
!AccIns
glomap_variables_local%mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0]
!CorIns
glomap_variables_local%mfrac_0(8,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0]
!SupIns

! Molar masses of components (kg mol-1)
glomap_variables_local%mm(1:ncp)=                                              &
            [0.098, 0.012, 0.0168, 0.05844, 0.100, 0.0168, 0.062, 0.084, 0.018]
!            h2so4  bc     oc      nacl     dust   so      no3    nano3  nh4
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
glomap_variables_local%rhocomp(1:ncp) =                                        &
       [1769.0, 1500.0, 1500.0, 1600.0, 2650.0, 1500.0, 1500.0, 1600.0, 1769.0]
!       h2so4   bc      oc      nacl    dust    so      no3     nano3   nh4

IF ( l_radaer_in ) THEN
  SELECT CASE (i_tune_bc_in)
  CASE (i_ukca_bc_tuned)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_tuned
  CASE (i_ukca_bc_mg_mix)
    glomap_variables_local%rhocomp(cp_bc) = rho_bc_mg_mix
  END SELECT
END IF

IF ( l_fix_nacl_density_in ) THEN
  glomap_variables_local%rhocomp(cp_cl) = rho_nacl
END IF

DO imode=1,nmodes
  glomap_variables_local%ddpmid(imode) = EXP( 0.5 *                            &
                               ( LOG(glomap_variables_local%ddplim0(imode) ) + &
                                 LOG(glomap_variables_local%ddplim1(imode) ) ) )
END DO

DO imode=1,nmodes
  rhommav=0.0
  DO icp=1,ncp
    rhommav = rhommav +                                                        &
              glomap_variables_local%mfrac_0(imode,icp) *                      &
              ( glomap_variables_local%rhocomp(icp) /                          &
                glomap_variables_local%mm(icp) )
  END DO

  glomap_variables_local%mmid(imode) = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddpmid(imode)**3 ) * &
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mlo(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim0(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

  glomap_variables_local%mhi(imode)  = ( pi / 6.0 ) *                          &
                                 ( glomap_variables_local%ddplim1(imode)**3 ) *&
                                 ( rhommav * avogadro ) *                      &
                                 glomap_variables_local%x(imode)

END DO

! number of dissociating ions in soluble components
glomap_variables_local%no_ions(1:ncp)=[3.0,0.0,0.0,2.0,0.0,0.0,2.0,2.0,2.0]

! Fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
glomap_variables_local%fracbcem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
glomap_variables_local%fracocem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]

! Set logical variables
glomap_variables_local%mode      = ( glomap_variables_local%mode_choice > 0 )
glomap_variables_local%component = .FALSE.
glomap_variables_local%soluble   = .FALSE.
DO imode=1,nmodes
  DO icp=1,ncp

    IF ( ( ( glomap_variables_local%component_mode( imode, icp ) == 1 ) .AND.  &
           ( glomap_variables_local%component_choice( icp ) == 1 ) )    .AND.  &
         ( glomap_variables_local%mode_choice(imode) == 1 ) ) THEN
      glomap_variables_local%component( imode, icp ) = .TRUE.
    END IF

    IF ( glomap_variables_local%soluble_choice(icp) == 1 ) THEN
      glomap_variables_local%soluble(icp) = .TRUE.
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcocduntnh_8mode_8cpt

! ############################################################################

SUBROUTINE ukca_mode_allocate_ctl_vars( ncp , glomap_variables_local )
! ---------------------------------------------------------------------|
!  Subroutine to allocate arrays in this module
! ---------------------------------------------------------------------|

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: ncp

TYPE(glomap_variables_type), INTENT(IN OUT) :: glomap_variables_local

! Local variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_ALLOCATE_CTL_VARS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(glomap_variables_local%component_choice))                  &
          ALLOCATE( glomap_variables_local%component_choice(ncp))
IF (.NOT. ALLOCATED(glomap_variables_local%soluble_choice))                    &
          ALLOCATE( glomap_variables_local%soluble_choice(ncp))
IF (.NOT. ALLOCATED(glomap_variables_local%mm))                                &
          ALLOCATE( glomap_variables_local%mm(ncp))
IF (.NOT. ALLOCATED(glomap_variables_local%rhocomp))                           &
          ALLOCATE( glomap_variables_local%rhocomp(ncp))
IF (.NOT. ALLOCATED(glomap_variables_local%no_ions))                           &
          ALLOCATE( glomap_variables_local%no_ions(ncp))
IF (.NOT. ALLOCATED(glomap_variables_local%component_names))                   &
          ALLOCATE( glomap_variables_local%component_names(ncp))
IF (.NOT. ALLOCATED(glomap_variables_local%soluble))                           &
          ALLOCATE( glomap_variables_local%soluble(ncp))
IF (.NOT. ALLOCATED(glomap_variables_local%component_mode))                    &
          ALLOCATE( glomap_variables_local%component_mode(nmodes,ncp))
IF (.NOT. ALLOCATED(glomap_variables_local%mfrac_0))                           &
          ALLOCATE( glomap_variables_local%mfrac_0(nmodes,ncp))
IF (.NOT. ALLOCATED(glomap_variables_local%component))                         &
          ALLOCATE( glomap_variables_local%component(nmodes,ncp))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_allocate_ctl_vars

END MODULE ukca_mode_setup
