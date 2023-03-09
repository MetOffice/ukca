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
!      UKCA_MODE_SUSSBCOCDU_4MODE
!      UKCA_MODE_SUSSBCOCNTNH_5MODE_7CPT
!      UKCA_MODE_SOLINSOL_6MODE
!    which define modes and components for different components/modes setup.
!
!  The internal routine: UKCA_MODE_ALLOCATE_CTL_VARS
!  is used to allocate arrays after the No. of components (ncp) is defined.
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
USE ukca_config_specification_mod, ONLY: glomap_config, i_ukca_bc_tuned, i_ukca_bc_mg_mix
USE chemistry_constants_mod, ONLY: avogadro, rho_so4
USE ukca_um_legacy_mod,  ONLY: pi, l_glomap_clim_radaer, i_glomap_clim_tune_bc
IMPLICIT NONE

PUBLIC
SAVE

INTEGER, PARAMETER :: nmodes=7    ! No of modes
INTEGER, PARAMETER :: ncp_max=9   ! No of components
INTEGER, PARAMETER :: ncation=3   ! No possible cation species
INTEGER, PARAMETER :: nanion =4   ! No possible anion species

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

INTEGER, PARAMETER :: moment_number =0 ! Index of number moment
INTEGER, PARAMETER :: moment_mass   =3 ! Index of mass moment

! Internal IDs for the types of UKCA aerosol modes used by RADAER.
INTEGER, PARAMETER :: ip_ukca_mode_nucleation = 0
INTEGER, PARAMETER :: ip_ukca_mode_aitken     = 1
INTEGER, PARAMETER :: ip_ukca_mode_accum      = 2
INTEGER, PARAMETER :: ip_ukca_mode_coarse     = 3

! No. of components
INTEGER :: ncp
! Mode switches (1=on, 0=0ff)
INTEGER :: mode_choice(nmodes)
! Modes resulting when two modes coagulate
INTEGER :: coag_mode(nmodes,nmodes)
! Specify which modes are soluble
INTEGER :: modesol(nmodes)

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
               mode_ait_sol,   mode_ait_sol,   mode_ait_sol,                   &
               mode_acc_sol,   mode_acc_sol,   mode_acc_sol,   mode_acc_sol,   &
               mode_cor_sol,   mode_cor_sol,   mode_cor_sol,   mode_cor_sol,   &
               mode_ait_insol, mode_ait_insol ]
! List of components used by sussbcoc_5mode (by component)
INTEGER :: component_list_by_cp_sussbcoc_5mode(ncp_list_sussbcoc_5mode) =      &
            [ cp_su,          cp_oc,                                           &
               cp_su,          cp_bc,          cp_oc,                          &
               cp_su,          cp_bc,          cp_oc,          cp_cl,          &
               cp_su,          cp_bc,          cp_oc,          cp_cl,          &
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

INTEGER :: i_tune_bc_loc ! local copy of i_tune_ukca_bc or i_glomap_clim_tune_bc

! Fraction of bc ems to go into each mode
REAL :: fracbcem(nmodes)
! Fraction of om ems to go into each mode
REAL :: fracocem(nmodes)
! Lower size limits of geometric mean diameter for each mode
REAL :: ddplim0(nmodes)
! Upper size limits of geometric mean diameter for each mode
REAL :: ddplim1(nmodes)
! Mid-point of size mode (m)
REAL :: ddpmid(nmodes)
! Mid-point masses for initial radius grid
REAL :: mmid(nmodes)
! Lo-interf masses for initial radius grid
REAL :: mlo(nmodes)
! Hi-interf masses for initial radius grid
REAL :: mhi(nmodes)
! Fixed geometric standard deviation for each mode
REAL :: sigmag(nmodes)
! EXP((9/2)*LOG^2(SIGMA_G))
REAL :: x(nmodes)
! Threshold for number in mode to carry out calculations
REAL :: num_eps(nmodes)

REAL :: rhommav
! specifies average (rho/mm) for default composition given by mfrac_0

! Mode names - the same for all MODE set ups, different ones are
! on for different set ups - see mode_choice
! Set case to be consistent with names in ukca_fieldname_mod for easier pattern
! matching
CHARACTER(LEN=7), PARAMETER :: mode_names(nmodes) =                            &
    ['Nuc_SOL','Ait_SOL','Acc_SOL','Cor_SOL',                                  &
                'Ait_INS','Acc_INS','Cor_INS']

! Modes (T/F)
LOGICAL     :: mode(nmodes) = .FALSE.

! Allocatable arrays:
! -------------------

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
LOGICAL, ALLOCATABLE :: component(:,:)
! Components that are soluble
LOGICAL, ALLOCATABLE :: soluble(:)
! Component names
CHARACTER(LEN=7),ALLOCATABLE :: component_names(:)

PRIVATE :: ukca_mode_allocate_ctl_vars

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_MODE_SETUP'

CONTAINS


! Subroutine Interface:

SUBROUTINE ukca_mode_allcp_4mode
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components with all components
!  switched on but only 4 modes used.
! ---------------------------------------------------------------------|
IMPLICIT NONE

INTEGER :: imode
INTEGER :: icp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_ALLCP_4MODE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncp = 6
CALL ukca_mode_allocate_ctl_vars()
! Component names
component_names =                                                              &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']
! Mode switches (1=on, 0=0ff)
mode_choice=[1,1,1,1,0,0,0]
! Specify which modes are soluble
modesol=[1,1,1,1,0,0,0]
! Component switches (1=on, 0=off)
component_choice=[1,1,1,1,1,0]
! ***n.b. in above have kept all cpts on (not SO) for UM test***
! Components that are soluble
soluble_choice=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
component_mode(1,1:ncp)=[1,0,1,0,0,1]       !allowed in nuc_sol
component_mode(2,1:ncp)=[1,1,1,0,0,1]       !allowed in ait_sol
component_mode(3,1:ncp)=[1,1,1,1,1,1]       !allowed in acc_sol
component_mode(4,1:ncp)=[1,1,1,1,1,1]       !allowed in cor_sol
component_mode(5,1:ncp)=[0,1,1,0,0,0]       !allowed in ait_ins
component_mode(6,1:ncp)=[0,0,0,0,1,0]       !allowed in acc_ins
component_mode(7,1:ncp)=[0,0,0,0,1,0]       !allowed in cor_ins

! Specify size limits of geometric mean diameter for each mode
! Set dlim34 here to be 500nm to agree with bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]

! Specify fixed geometric standard deviation for each mode
sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0] ! to match M7, but sigacc=1.4

DO imode=1,nmodes
  x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
num_eps=[1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14]

! Initial fractions of mass in each mode among components
mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !nucln. soluble
mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Aitken soluble
mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !accum. soluble
mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !coarse soluble
mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Aitken insoluble
mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !accum. insoluble
mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !coarse insoluble

! Molar masses of components (kg mol-1)
mm=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)
! Mass density of components (kg m^-3)
rhocomp = [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

! Adjust the density of BC?
IF (glomap_config%l_ukca_radaer) THEN
   i_tune_bc_loc = glomap_config%i_ukca_tune_bc
ELSE IF (l_glomap_clim_radaer) THEN
   i_tune_bc_loc =  i_glomap_clim_tune_bc
ELSE
   ! Do not tune density of BC
   i_tune_bc_loc = 0 
END IF
SELECT CASE (i_tune_bc_loc)
   CASE (i_ukca_bc_tuned)
      rhocomp(cp_bc) = rho_bc_tuned
   CASE (i_ukca_bc_mg_mix)
      rhocomp(cp_bc) = rho_bc_mg_mix
END SELECT

! Adjust the density of NaCl (sea-salt)
IF (glomap_config%l_fix_nacl_density) rhocomp(cp_cl) = rho_nacl

DO imode=1,nmodes
  ddpmid(imode)=EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
  rhommav=0.0
  DO icp=1,ncp
    rhommav=rhommav+mfrac_0(imode,icp)*(rhocomp(icp)/mm(icp))
  END DO
  mmid(imode)=(pi/6.0)*(ddpmid (imode)**3)*(rhommav*avogadro)*x(imode)
  mlo (imode)=(pi/6.0)*(ddplim0(imode)**3)*(rhommav*avogadro)*x(imode)
  mhi (imode)=(pi/6.0)*(ddplim1(imode)**3)*(rhommav*avogadro)*x(imode)
END DO

! Modes resulting when two modes coagulate
coag_mode(1,1:nmodes)=[1,2,3,4,2,3,4]
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
coag_mode(2,1:nmodes)=[2,2,3,4,2,3,4]
! 2 coagulating with:    1,2,3,4,5,6,7
coag_mode(3,1:nmodes)=[3,3,3,4,3,3,4]
! 3 coagulating with:    1,2,3,4,5,6,7
coag_mode(4,1:nmodes)=[4,4,4,4,4,4,4]
! 4 coagulating with:    1,2,3,4,5,6,7
coag_mode(5,1:nmodes)=[2,2,3,4,5,6,7]
! 5 coagulating with:    1,2,3,4,5,6,7
coag_mode(6,1:nmodes)=[3,3,3,4,6,6,7]
! 6 coagulating with:    1,2,3,4,5,6,7
coag_mode(7,1:nmodes)=[4,4,4,4,7,7,7]
! 7 coagulating with:    1,2,3,4,5,6,7

! number of dissociating ions in soluble components
no_ions=[3.0,0.0,0.0,2.0,0.0,0.0]
!
fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
mode=(mode_choice > 0)
component=.FALSE.
soluble=.FALSE.
DO imode=1,nmodes
  DO icp=1,ncp
    IF (((component_mode(imode,icp) == 1) .AND.                                &
          (component_choice(icp) == 1)) .AND.                                  &
          (mode_choice(imode) == 1)) THEN
      component(imode,icp)=.TRUE.
    END IF
    IF (soluble_choice(icp) == 1) soluble(icp)=.TRUE.
  END DO
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_allcp_4mode
SUBROUTINE ukca_mode_suss_4mode
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate and sea-salt only in 4 modes.
!  Uses 10 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

INTEGER :: imode
INTEGER :: icp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSS_4MODE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncp = 6
CALL ukca_mode_allocate_ctl_vars()
! Component names
component_names =                                                              &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']
! Mode switches (1=on, 0=0ff)
mode_choice=[1,1,1,1,0,0,0]
! Specify which modes are soluble
modesol=[1,1,1,1,0,0,0]
! Component switches (1=on, 0=off)
component_choice=[1,0,0,1,0,0]
! *** n.b. only have h2so4 and nacl cpts on for this setup ***
! Components that are soluble
soluble_choice=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
component_mode(1,1:ncp)=[1,0,1,0,0,1]       !allowed in nuc_sol
component_mode(2,1:ncp)=[1,1,1,0,0,1]       !allowed in ait_sol
component_mode(3,1:ncp)=[1,1,1,1,1,1]       !allowed in acc_sol
component_mode(4,1:ncp)=[1,1,1,1,1,1]       !allowed in cor_sol
component_mode(5,1:ncp)=[0,1,1,0,0,0]       !allowed in ait_ins
component_mode(6,1:ncp)=[0,0,0,0,1,0]       !allowed in acc_ins
component_mode(7,1:ncp)=[0,0,0,0,1,0]       !allowed in cor_ins

! Specify size limits of geometric mean diameter for each mode
!! set dplim34 here to be 500nm to match value found from bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]

! Specify fixed geometric standard deviation for each mode
sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0] ! to match M7, but sigacc=1.4

DO imode=1,nmodes
  x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
num_eps=[1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14]

! Initial fractions of mass in each mode among components
mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !nucln. soluble
mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Aitken soluble
mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !accum. soluble
mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !coarse soluble
mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Aitken insoluble
mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !accum. insoluble
mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !coarse insoluble

! Molar masses of components (kg mol-1)
mm=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)
! Mass density of components (kg m^-3)
rhocomp = [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

! Adjust the density of BC?
IF (glomap_config%l_ukca_radaer) THEN
   i_tune_bc_loc = glomap_config%i_ukca_tune_bc
ELSE IF (l_glomap_clim_radaer) THEN
   i_tune_bc_loc =  i_glomap_clim_tune_bc
ELSE
   ! Do not tune density of BC
   i_tune_bc_loc = 0 
END IF
SELECT CASE (i_tune_bc_loc)
   CASE (i_ukca_bc_tuned)
      rhocomp(cp_bc) = rho_bc_tuned
   CASE (i_ukca_bc_mg_mix)
      rhocomp(cp_bc) = rho_bc_mg_mix
END SELECT

! Adjust the density of NaCl (sea-salt)
IF (glomap_config%l_fix_nacl_density) rhocomp(cp_cl) = rho_nacl

DO imode=1,nmodes
  ddpmid(imode)=EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
  rhommav=0.0
  DO icp=1,ncp
    rhommav=rhommav+mfrac_0(imode,icp)*(rhocomp(icp)/mm(icp))
  END DO
  mmid(imode)=(pi/6.0)*(ddpmid (imode)**3)*(rhommav*avogadro)*x(imode)
  mlo (imode)=(pi/6.0)*(ddplim0(imode)**3)*(rhommav*avogadro)*x(imode)
  mhi (imode)=(pi/6.0)*(ddplim1(imode)**3)*(rhommav*avogadro)*x(imode)
END DO

! Modes resulting when two modes coagulate
coag_mode(1,1:nmodes)=[1,2,3,4,2,3,4]
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
coag_mode(2,1:nmodes)=[2,2,3,4,2,3,4]
! 2 coagulating with:    1,2,3,4,5,6,7
coag_mode(3,1:nmodes)=[3,3,3,4,3,3,4]
! 3 coagulating with:    1,2,3,4,5,6,7
coag_mode(4,1:nmodes)=[4,4,4,4,4,4,4]
! 4 coagulating with:    1,2,3,4,5,6,7
coag_mode(5,1:nmodes)=[2,2,3,4,5,6,7]
! 5 coagulating with:    1,2,3,4,5,6,7
coag_mode(6,1:nmodes)=[3,3,3,4,6,6,7]
! 6 coagulating with:    1,2,3,4,5,6,7
coag_mode(7,1:nmodes)=[4,4,4,4,7,7,7]
! 7 coagulating with:    1,2,3,4,5,6,7

! number of dissociating ions in soluble components
no_ions=[3.0,0.0,0.0,2.0,0.0,0.0]
!
fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
mode=(mode_choice > 0)
component=.FALSE.
soluble=.FALSE.
DO imode=1,nmodes
  DO icp=1,ncp
    IF (((component_mode(imode,icp) == 1) .AND.                                &
          (component_choice(icp) == 1)) .AND.                                  &
          (mode_choice(imode) == 1)) THEN
      component(imode,icp)=.TRUE.
    END IF
    IF (soluble_choice(icp) == 1) soluble(icp)=.TRUE.
  END DO
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_suss_4mode
SUBROUTINE ukca_mode_sussbcocdu_4mode
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  SO4, sea-salt, bc, oc (secondary & primary combined) & du in 4 modes.
!  Uses 19 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

INTEGER :: imode
INTEGER :: icp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOCDU_4MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncp = 6
CALL ukca_mode_allocate_ctl_vars()
! Component names
component_names =                                                              &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']
! Mode switches (1=on, 0=0ff)
mode_choice=[1,1,1,1,0,0,0]
! Specify which modes are soluble
modesol=[1,1,1,1,0,0,0]
! Component switches (1=on, 0=off)
component_choice=[1,1,1,1,1,0]
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
soluble_choice=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
component_mode(1,1:ncp)=[1,0,1,0,0,1]       !allowed in nuc_sol
component_mode(2,1:ncp)=[1,1,1,0,0,1]       !allowed in ait_sol
component_mode(3,1:ncp)=[1,1,1,1,1,1]       !allowed in acc_sol
component_mode(4,1:ncp)=[1,1,1,1,1,1]       !allowed in cor_sol
component_mode(5,1:ncp)=[0,1,1,0,0,0]       !allowed in ait_ins
component_mode(6,1:ncp)=[0,0,0,0,1,0]       !allowed in acc_ins
component_mode(7,1:ncp)=[0,0,0,0,1,0]       !allowed in cor_ins

! Specify size limits of geometric mean diameter for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]

! Specify fixed geometric standard deviation for each mode
sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0] ! to match M7, but sigacc=1.4

DO imode=1,nmodes
  x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
num_eps=[1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14]

! Initial fractions of mass in each mode among components
mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !nucln. soluble
mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Aitken soluble
mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !accum. soluble
mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !coarse soluble
mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Aitken insoluble
mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !accum. insoluble
mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !coarse insoluble

! Molar masses of components (kg mol-1)
mm=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)
! Mass density of components (kg m^-3)
rhocomp = [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

! Adjust the density of BC?
IF (glomap_config%l_ukca_radaer) THEN
   i_tune_bc_loc = glomap_config%i_ukca_tune_bc
ELSE IF (l_glomap_clim_radaer) THEN
   i_tune_bc_loc =  i_glomap_clim_tune_bc
ELSE
   ! Do not tune BC density
   i_tune_bc_loc = 0
END IF
SELECT CASE (i_tune_bc_loc)
   CASE (i_ukca_bc_tuned)
      rhocomp(cp_bc) = rho_bc_tuned
   CASE (i_ukca_bc_mg_mix)
      rhocomp(cp_bc) = rho_bc_mg_mix
END SELECT

! Adjust the density of NaCl (sea-salt)
IF (glomap_config%l_fix_nacl_density) rhocomp(cp_cl) = rho_nacl

DO imode=1,nmodes
  ddpmid(imode)=EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
  rhommav=0.0
  DO icp=1,ncp
    rhommav=rhommav+mfrac_0(imode,icp)*(rhocomp(icp)/mm(icp))
  END DO
  mmid(imode)=(pi/6.0)*(ddpmid (imode)**3)*(rhommav*avogadro)*x(imode)
  mlo (imode)=(pi/6.0)*(ddplim0(imode)**3)*(rhommav*avogadro)*x(imode)
  mhi (imode)=(pi/6.0)*(ddplim1(imode)**3)*(rhommav*avogadro)*x(imode)
END DO

! Modes resulting when two modes coagulate
coag_mode(1,1:nmodes)=[1,2,3,4,2,3,4]
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
coag_mode(2,1:nmodes)=[2,2,3,4,2,3,4]
! 2 coagulating with:    1,2,3,4,5,6,7
coag_mode(3,1:nmodes)=[3,3,3,4,3,3,4]
! 3 coagulating with:    1,2,3,4,5,6,7
coag_mode(4,1:nmodes)=[4,4,4,4,4,4,4]
! 4 coagulating with:    1,2,3,4,5,6,7
coag_mode(5,1:nmodes)=[2,2,3,4,5,6,7]
! 5 coagulating with:    1,2,3,4,5,6,7
coag_mode(6,1:nmodes)=[3,3,3,4,6,6,7]
! 6 coagulating with:    1,2,3,4,5,6,7
coag_mode(7,1:nmodes)=[4,4,4,4,7,7,7]
! 7 coagulating with:    1,2,3,4,5,6,7

! number of dissociating ions in soluble components
no_ions=[3.0,0.0,0.0,2.0,0.0,0.0]
!
fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
mode=(mode_choice > 0)
component=.FALSE.
soluble=.FALSE.
DO imode=1,nmodes
  DO icp=1,ncp
    IF (((component_mode(imode,icp) == 1) .AND.                                &
          (component_choice(icp) == 1)) .AND.                                  &
          (mode_choice(imode) == 1)) THEN
      component(imode,icp)=.TRUE.
    END IF
    IF (soluble_choice(icp) == 1) soluble(icp)=.TRUE.
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcocdu_4mode
SUBROUTINE ukca_mode_sussbcocdu_7mode
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  SO4, sea-salt, bc, oc (secondary & primary combined) & du in 7 modes.
!  Uses 26 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

INTEGER :: imode
INTEGER :: icp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOCDU_7MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncp = 6
CALL ukca_mode_allocate_ctl_vars()
! Component names
component_names =                                                              &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
mode_choice=[1,1,1,1,1,1,1]
! Specify which modes are soluble
modesol=[1,1,1,1,0,0,0]
! Component switches (1=on, 0=off)
component_choice=[1,1,1,1,1,0]
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
soluble_choice=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
component_mode(1,1:ncp)=[1,0,1,0,0,1]       !allowed in nuc_sol
component_mode(2,1:ncp)=[1,1,1,0,0,1]       !allowed in ait_sol
component_mode(3,1:ncp)=[1,1,1,1,1,1]       !allowed in acc_sol
component_mode(4,1:ncp)=[1,1,1,1,1,1]       !allowed in cor_sol
component_mode(5,1:ncp)=[0,1,1,0,0,0]       !allowed in ait_ins
component_mode(6,1:ncp)=[0,0,0,0,1,0]       !allowed in acc_ins
component_mode(7,1:ncp)=[0,0,0,0,1,0]       !allowed in cor_ins

! Specify size limits of geometric mean diameter for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]

! Specify fixed geometric standard deviation for each mode
sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0] ! to match M7, but sigacc=1.4

DO imode=1,nmodes
  x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
num_eps=[1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14]

! Initial fractions of mass in each mode among components
mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !nucln. soluble
mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Aitken soluble
mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !accum. soluble
mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !coarse soluble
mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Aitken insoluble
mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !accum. insoluble
mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !coarse insoluble

! Molar masses of components (kg mol-1)
mm=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)
! Mass density of components (kg m^-3)
rhocomp = [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

! Adjust the density of BC?
IF (glomap_config%l_ukca_radaer) THEN
   i_tune_bc_loc = glomap_config%i_ukca_tune_bc
ELSE IF (l_glomap_clim_radaer) THEN
   i_tune_bc_loc =  i_glomap_clim_tune_bc
ELSE
   ! Do not tune BC density
   i_tune_bc_loc = 0 
END IF
SELECT CASE (i_tune_bc_loc)
   CASE (i_ukca_bc_tuned)
      rhocomp(cp_bc) = rho_bc_tuned
   CASE (i_ukca_bc_mg_mix)
      rhocomp(cp_bc) = rho_bc_mg_mix
END SELECT

! Adjust the density of NaCl (sea-salt)
IF (glomap_config%l_fix_nacl_density) rhocomp(cp_cl) = rho_nacl

DO imode=1,nmodes
  ddpmid(imode)=EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
  rhommav=0.0
  DO icp=1,ncp
    rhommav=rhommav+mfrac_0(imode,icp)*(rhocomp(icp)/mm(icp))
  END DO
  mmid(imode)=(pi/6.0)*(ddpmid (imode)**3)*(rhommav*avogadro)*x(imode)
  mlo (imode)=(pi/6.0)*(ddplim0(imode)**3)*(rhommav*avogadro)*x(imode)
  mhi (imode)=(pi/6.0)*(ddplim1(imode)**3)*(rhommav*avogadro)*x(imode)
END DO

! Modes resulting when two modes coagulate
coag_mode(1,1:nmodes)=[1,2,3,4,2,3,4]
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
coag_mode(2,1:nmodes)=[2,2,3,4,2,3,4]
! 2 coagulating with:    1,2,3,4,5,6,7
coag_mode(3,1:nmodes)=[3,3,3,4,3,3,4]
! 3 coagulating with:    1,2,3,4,5,6,7
coag_mode(4,1:nmodes)=[4,4,4,4,4,4,4]
! 4 coagulating with:    1,2,3,4,5,6,7
coag_mode(5,1:nmodes)=[2,2,3,4,5,6,7]
! 5 coagulating with:    1,2,3,4,5,6,7
coag_mode(6,1:nmodes)=[3,3,3,4,6,6,7]
! 6 coagulating with:    1,2,3,4,5,6,7
coag_mode(7,1:nmodes)=[4,4,4,4,7,7,7]
! 7 coagulating with:    1,2,3,4,5,6,7

! number of dissociating ions in soluble components
no_ions=[3.0,0.0,0.0,2.0,0.0,0.0]
!
fracbcem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0]
fracocem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
mode=(mode_choice > 0)
component=.FALSE.
soluble=.FALSE.
DO imode=1,nmodes
  DO icp=1,ncp
    IF (((component_mode(imode,icp) == 1) .AND.                                &
          (component_choice(icp) == 1)) .AND.                                  &
          (mode_choice(imode) == 1)) THEN
      component(imode,icp)=.TRUE.
    END IF
    IF (soluble_choice(icp) == 1) soluble(icp)=.TRUE.
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcocdu_7mode
SUBROUTINE ukca_mode_sussbcoc_4mode
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate, sea-salt, bc & oc (secondary & primary combined) in 4 modes.
!  Uses 17 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

INTEGER :: imode
INTEGER :: icp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOC_4MODE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncp = 6
CALL ukca_mode_allocate_ctl_vars()
! Component names
component_names =                                                              &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']
! Mode switches (1=on, 0=0ff)
mode_choice=[1,1,1,1,0,0,0]
! Specify which modes are soluble
modesol=[1,1,1,1,0,0,0]
! Component switches (1=on, 0=off)
component_choice=[1,1,1,1,0,0]
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
soluble_choice=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
component_mode(1,1:ncp)=[1,0,1,0,0,1]       !allowed in nuc_sol
component_mode(2,1:ncp)=[1,1,1,0,0,1]       !allowed in ait_sol
component_mode(3,1:ncp)=[1,1,1,1,1,1]       !allowed in acc_sol
component_mode(4,1:ncp)=[1,1,1,1,1,1]       !allowed in cor_sol
component_mode(5,1:ncp)=[0,1,1,0,0,0]       !allowed in ait_ins
component_mode(6,1:ncp)=[0,0,0,0,1,0]       !allowed in acc_ins
component_mode(7,1:ncp)=[0,0,0,0,1,0]       !allowed in cor_ins

! Specify size limits of geometric mean diameter for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]

! Specify fixed geometric standard deviation for each mode
sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0] ! to match M7, but sigacc=1.4

DO imode=1,nmodes
  x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
num_eps=[1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14]

! Initial fractions of mass in each mode among components
mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !nucln. soluble
mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Aitken soluble
mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !accum. soluble
mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !coarse soluble
mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Aitken insoluble
mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !accum. insoluble
mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !coarse insoluble

! Molar masses of components (kg mol-1)
mm=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)
! Mass density of components (kg m^-3)
rhocomp = [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

! Adjust the density of BC?
IF (glomap_config%l_ukca_radaer) THEN
   i_tune_bc_loc = glomap_config%i_ukca_tune_bc
ELSE IF (l_glomap_clim_radaer) THEN
   i_tune_bc_loc =  i_glomap_clim_tune_bc
ELSE
   ! Do not tune BC density
   i_tune_bc_loc = 0 
END IF
SELECT CASE (i_tune_bc_loc)
   CASE (i_ukca_bc_tuned)
      rhocomp(cp_bc) = rho_bc_tuned
   CASE (i_ukca_bc_mg_mix)
      rhocomp(cp_bc) = rho_bc_mg_mix
END SELECT

! Adjust the density of NaCl (sea-salt)
IF (glomap_config%l_fix_nacl_density) rhocomp(cp_cl) = rho_nacl

DO imode=1,nmodes
  ddpmid(imode)=EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
  rhommav=0.0
  DO icp=1,ncp
    rhommav=rhommav+mfrac_0(imode,icp)*(rhocomp(icp)/mm(icp))
  END DO
  mmid(imode)=(pi/6.0)*(ddpmid (imode)**3)*(rhommav*avogadro)*x(imode)
  mlo (imode)=(pi/6.0)*(ddplim0(imode)**3)*(rhommav*avogadro)*x(imode)
  mhi (imode)=(pi/6.0)*(ddplim1(imode)**3)*(rhommav*avogadro)*x(imode)
END DO

! Modes resulting when two modes coagulate
coag_mode(1,1:nmodes)=[1,2,3,4,2,3,4]
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
coag_mode(2,1:nmodes)=[2,2,3,4,2,3,4]
! 2 coagulating with:    1,2,3,4,5,6,7
coag_mode(3,1:nmodes)=[3,3,3,4,3,3,4]
! 3 coagulating with:    1,2,3,4,5,6,7
coag_mode(4,1:nmodes)=[4,4,4,4,4,4,4]
! 4 coagulating with:    1,2,3,4,5,6,7
coag_mode(5,1:nmodes)=[2,2,3,4,5,6,7]
! 5 coagulating with:    1,2,3,4,5,6,7
coag_mode(6,1:nmodes)=[3,3,3,4,6,6,7]
! 6 coagulating with:    1,2,3,4,5,6,7
coag_mode(7,1:nmodes)=[4,4,4,4,7,7,7]
! 7 coagulating with:    1,2,3,4,5,6,7

! number of dissociating ions in soluble components
no_ions=[3.0,0.0,0.0,2.0,0.0,0.0]
!
fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
mode=(mode_choice > 0)
component=.FALSE.
soluble=.FALSE.
DO imode=1,nmodes
  DO icp=1,ncp
    IF (((component_mode(imode,icp) == 1) .AND.                                &
          (component_choice(icp) == 1)) .AND.                                  &
          (mode_choice(imode) == 1)) THEN
      component(imode,icp)=.TRUE.
    END IF
    IF (soluble_choice(icp) == 1) soluble(icp)=.TRUE.
  END DO
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcoc_4mode
SUBROUTINE ukca_mode_sussbcoc_5mode
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate, sea-salt, bc & oc (secondary & primary combined) in 5 modes.
!  Uses 20 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

INTEGER :: imode
INTEGER :: icp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOC_5MODE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncp = 6
CALL ukca_mode_allocate_ctl_vars()
! Component names
component_names =                                                              &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']
! Mode switches (1=on, 0=0ff)
mode_choice=[1,1,1,1,1,0,0]
! Specify which modes are soluble
modesol=[1,1,1,1,0,0,0]
! Component switches (1=on, 0=off)
component_choice=[1,1,1,1,0,0]
! *** n.b. only have h2so4,bc,oc,nacl cpts on for this setup ***
! Components that are soluble
soluble_choice=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
component_mode(1,1:ncp)=[1,0,1,0,0,1]       !allowed in nuc_sol
component_mode(2,1:ncp)=[1,1,1,0,0,1]       !allowed in ait_sol
component_mode(3,1:ncp)=[1,1,1,1,1,1]       !allowed in acc_sol
component_mode(4,1:ncp)=[1,1,1,1,1,1]       !allowed in cor_sol
component_mode(5,1:ncp)=[0,1,1,0,0,0]       !allowed in ait_ins
component_mode(6,1:ncp)=[0,0,0,0,1,0]       !allowed in acc_ins
component_mode(7,1:ncp)=[0,0,0,0,1,0]       !allowed in cor_ins

! Specify size limits of geometric mean diameter for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]

! Specify fixed geometric standard deviation for each mode
sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0] ! to match M7, but sigacc=1.4

DO imode=1,nmodes
  x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
num_eps=[1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14]

! Initial fractions of mass in each mode among components
mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !nucln. soluble
mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Aitken soluble
mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !accum. soluble
mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !coarse soluble
mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Aitken insoluble
mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !accum. insoluble
mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !coarse insoluble

! Molar masses of components (kg mol-1)
mm=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)
! Mass density of components (kg m^-3)
rhocomp = [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

! Adjust the density of BC?
IF (glomap_config%l_ukca_radaer) THEN
   i_tune_bc_loc = glomap_config%i_ukca_tune_bc
ELSE IF (l_glomap_clim_radaer) THEN
   i_tune_bc_loc =  i_glomap_clim_tune_bc
ELSE
   ! Do not tune BC density
   i_tune_bc_loc = 0 
END IF
SELECT CASE (i_tune_bc_loc)
   CASE (i_ukca_bc_tuned)
      rhocomp(cp_bc) = rho_bc_tuned
   CASE (i_ukca_bc_mg_mix)
      rhocomp(cp_bc) = rho_bc_mg_mix
END SELECT

! Adjust the density of NaCl (sea-salt)
IF (glomap_config%l_fix_nacl_density) rhocomp(cp_cl) = rho_nacl

DO imode=1,nmodes
  ddpmid(imode)=EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
  rhommav=0.0
  DO icp=1,ncp
    rhommav=rhommav+mfrac_0(imode,icp)*(rhocomp(icp)/mm(icp))
  END DO
  mmid(imode)=(pi/6.0)*(ddpmid (imode)**3)*(rhommav*avogadro)*x(imode)
  mlo (imode)=(pi/6.0)*(ddplim0(imode)**3)*(rhommav*avogadro)*x(imode)
  mhi (imode)=(pi/6.0)*(ddplim1(imode)**3)*(rhommav*avogadro)*x(imode)
END DO

! Modes resulting when two modes coagulate
coag_mode(1,1:nmodes)=[1,2,3,4,2,3,4]
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
coag_mode(2,1:nmodes)=[2,2,3,4,2,3,4]
! 2 coagulating with:    1,2,3,4,5,6,7
coag_mode(3,1:nmodes)=[3,3,3,4,3,3,4]
! 3 coagulating with:    1,2,3,4,5,6,7
coag_mode(4,1:nmodes)=[4,4,4,4,4,4,4]
! 4 coagulating with:    1,2,3,4,5,6,7
coag_mode(5,1:nmodes)=[2,2,3,4,5,6,7]
! 5 coagulating with:    1,2,3,4,5,6,7
coag_mode(6,1:nmodes)=[3,3,3,4,6,6,7]
! 6 coagulating with:    1,2,3,4,5,6,7
coag_mode(7,1:nmodes)=[4,4,4,4,7,7,7]
! 7 coagulating with:    1,2,3,4,5,6,7

! number of dissociating ions in soluble components
no_ions=[3.0,0.0,0.0,2.0,0.0,0.0]
!
fracbcem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0]
fracocem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
mode=(mode_choice > 0)
component=.FALSE.
soluble=.FALSE.
DO imode=1,nmodes
  DO icp=1,ncp
    IF (((component_mode(imode,icp) == 1) .AND.                                &
          (component_choice(icp) == 1)) .AND.                                  &
          (mode_choice(imode) == 1)) THEN
      component(imode,icp)=.TRUE.
    END IF
    IF (soluble_choice(icp) == 1) soluble(icp)=.TRUE.
  END DO
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcoc_5mode
SUBROUTINE ukca_mode_sussbcocso_4mode
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate, sea-salt, bc, primary oc & secondary oc cpts in 5 modes.
!  Uses 20 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

INTEGER :: imode
INTEGER :: icp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOCSO_4MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncp = 6
CALL ukca_mode_allocate_ctl_vars()
! Component names
component_names =                                                              &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']
! Mode switches (1=on, 0=0ff)
mode_choice=[1,1,1,1,0,0,0]
! Specify which modes are soluble
modesol=[1,1,1,1,0,0,0]
! Component switches (1=on, 0=off)
component_choice=[1,1,1,1,0,1] ! ***all cpts on except dust***
! Components that are soluble
soluble_choice=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
component_mode(1,1:ncp)=[1,0,1,0,0,1]       !allowed in nuc_sol
component_mode(2,1:ncp)=[1,1,1,0,0,1]       !allowed in ait_sol
component_mode(3,1:ncp)=[1,1,1,1,1,1]       !allowed in acc_sol
component_mode(4,1:ncp)=[1,1,1,1,1,1]       !allowed in cor_sol
component_mode(5,1:ncp)=[0,1,1,0,0,0]       !allowed in ait_ins
component_mode(6,1:ncp)=[0,0,0,0,1,0]       !allowed in acc_ins
component_mode(7,1:ncp)=[0,0,0,0,1,0]       !allowed in cor_ins

! Specify size limits of geometric mean diameter for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]

! Specify fixed geometric standard deviation for each mode
sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0] ! to match M7, but sigacc=1.4

DO imode=1,nmodes
  x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
num_eps=[1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14]

! Initial fractions of mass in each mode among components
mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !nucln. soluble
mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Aitken soluble
mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !accum. soluble
mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !coarse soluble
mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Aitken insoluble
mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !accum. insoluble
mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !coarse insoluble

! Molar masses of components (kg mol-1)
mm=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)
! Mass density of components (kg m^-3)
rhocomp = [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

! Adjust the density of BC?
IF (glomap_config%l_ukca_radaer) THEN
   i_tune_bc_loc = glomap_config%i_ukca_tune_bc
ELSE IF (l_glomap_clim_radaer) THEN
   i_tune_bc_loc =  i_glomap_clim_tune_bc
ELSE
   ! Do not tune BC density
   i_tune_bc_loc = 0 
END IF
SELECT CASE (i_tune_bc_loc)
   CASE (i_ukca_bc_tuned)
      rhocomp(cp_bc) = rho_bc_tuned
   CASE (i_ukca_bc_mg_mix)
      rhocomp(cp_bc) = rho_bc_mg_mix
END SELECT

! Adjust the density of NaCl (sea-salt)
IF (glomap_config%l_fix_nacl_density) rhocomp(cp_cl) = rho_nacl

DO imode=1,nmodes
  ddpmid(imode)=EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
  rhommav=0.0
  DO icp=1,ncp
    rhommav=rhommav+mfrac_0(imode,icp)*(rhocomp(icp)/mm(icp))
  END DO
  mmid(imode)=(pi/6.0)*(ddpmid (imode)**3)*(rhommav*avogadro)*x(imode)
  mlo (imode)=(pi/6.0)*(ddplim0(imode)**3)*(rhommav*avogadro)*x(imode)
  mhi (imode)=(pi/6.0)*(ddplim1(imode)**3)*(rhommav*avogadro)*x(imode)
END DO

! Modes resulting when two modes coagulate
coag_mode(1,1:nmodes)=[1,2,3,4,2,3,4]
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
coag_mode(2,1:nmodes)=[2,2,3,4,2,3,4]
! 2 coagulating with:    1,2,3,4,5,6,7
coag_mode(3,1:nmodes)=[3,3,3,4,3,3,4]
! 3 coagulating with:    1,2,3,4,5,6,7
coag_mode(4,1:nmodes)=[4,4,4,4,4,4,4]
! 4 coagulating with:    1,2,3,4,5,6,7
coag_mode(5,1:nmodes)=[2,2,3,4,5,6,7]
! 5 coagulating with:    1,2,3,4,5,6,7
coag_mode(6,1:nmodes)=[3,3,3,4,6,6,7]
! 6 coagulating with:    1,2,3,4,5,6,7
coag_mode(7,1:nmodes)=[4,4,4,4,7,7,7]
! 7 coagulating with:    1,2,3,4,5,6,7

! number of dissociating ions in soluble components
no_ions=[3.0,0.0,0.0,2.0,0.0,0.0]
!
fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into   soluble Aitken for this setup).
!
! Set logical variables
mode=(mode_choice > 0)
component=.FALSE.
soluble=.FALSE.
DO imode=1,nmodes
  DO icp=1,ncp
    IF (((component_mode(imode,icp) == 1) .AND.                                &
          (component_choice(icp) == 1)) .AND.                                  &
          (mode_choice(imode) == 1)) THEN
      component(imode,icp)=.TRUE.
    END IF
    IF (soluble_choice(icp) == 1) soluble(icp)=.TRUE.
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcocso_4mode
SUBROUTINE ukca_mode_sussbcocso_5mode
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate, sea-salt, bc, primary oc & secondary oc cpts in 5 modes.
!  Uses 23 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

INTEGER :: imode
INTEGER :: icp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOCSO_5MODE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncp = 6
CALL ukca_mode_allocate_ctl_vars()
! Component names
component_names =                                                              &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']
! Mode switches (1=on, 0=0ff)
mode_choice=[1,1,1,1,1,0,0]
! Specify which modes are soluble
modesol=[1,1,1,1,0,0,0]
! Component switches (1=on, 0=off)
component_choice=[1,1,1,1,0,1] ! ***all cpts on except dust***
! Components that are soluble
soluble_choice=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
component_mode(1,1:ncp)=[1,0,1,0,0,1]       !allowed in nuc_sol
component_mode(2,1:ncp)=[1,1,1,0,0,1]       !allowed in ait_sol
component_mode(3,1:ncp)=[1,1,1,1,1,1]       !allowed in acc_sol
component_mode(4,1:ncp)=[1,1,1,1,1,1]       !allowed in cor_sol
component_mode(5,1:ncp)=[0,1,1,0,0,0]       !allowed in ait_ins
component_mode(6,1:ncp)=[0,0,0,0,1,0]       !allowed in acc_ins
component_mode(7,1:ncp)=[0,0,0,0,1,0]       !allowed in cor_ins

! Specify size limits of geometric mean diameter for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]

! Specify fixed geometric standard deviation for each mode
sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0] ! to match M7, but sigacc=1.4

DO imode=1,nmodes
  x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
num_eps=[1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14]

! Initial fractions of mass in each mode among components
mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !nucln. soluble
mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Aitken soluble
mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !accum. soluble
mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !coarse soluble
mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Aitken insoluble
mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !accum. insoluble
mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !coarse insoluble

! Molar masses of components (kg mol-1)
mm=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)
! Mass density of components (kg m^-3)
rhocomp = [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

! Adjust the density of BC?
IF (glomap_config%l_ukca_radaer) THEN
   i_tune_bc_loc = glomap_config%i_ukca_tune_bc
ELSE IF (l_glomap_clim_radaer) THEN
   i_tune_bc_loc =  i_glomap_clim_tune_bc
ELSE
   ! Do not tune BC density
   i_tune_bc_loc = 0 
END IF
SELECT CASE (i_tune_bc_loc)
   CASE (i_ukca_bc_tuned)
      rhocomp(cp_bc) = rho_bc_tuned
   CASE (i_ukca_bc_mg_mix)
      rhocomp(cp_bc) = rho_bc_mg_mix
END SELECT

! Adjust the density of NaCl (sea-salt)
IF (glomap_config%l_fix_nacl_density) rhocomp(cp_cl) = rho_nacl

DO imode=1,nmodes
  ddpmid(imode)=EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
  rhommav=0.0
  DO icp=1,ncp
    rhommav=rhommav+mfrac_0(imode,icp)*(rhocomp(icp)/mm(icp))
  END DO
  mmid(imode)=(pi/6.0)*(ddpmid (imode)**3)*(rhommav*avogadro)*x(imode)
  mlo (imode)=(pi/6.0)*(ddplim0(imode)**3)*(rhommav*avogadro)*x(imode)
  mhi (imode)=(pi/6.0)*(ddplim1(imode)**3)*(rhommav*avogadro)*x(imode)
END DO

! Modes resulting when two modes coagulate
coag_mode(1,1:nmodes)=[1,2,3,4,2,3,4]
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
coag_mode(2,1:nmodes)=[2,2,3,4,2,3,4]
! 2 coagulating with:    1,2,3,4,5,6,7
coag_mode(3,1:nmodes)=[3,3,3,4,3,3,4]
! 3 coagulating with:    1,2,3,4,5,6,7
coag_mode(4,1:nmodes)=[4,4,4,4,4,4,4]
! 4 coagulating with:    1,2,3,4,5,6,7
coag_mode(5,1:nmodes)=[2,2,3,4,5,6,7]
! 5 coagulating with:    1,2,3,4,5,6,7
coag_mode(6,1:nmodes)=[3,3,3,4,6,6,7]
! 6 coagulating with:    1,2,3,4,5,6,7
coag_mode(7,1:nmodes)=[4,4,4,4,7,7,7]
! 7 coagulating with:    1,2,3,4,5,6,7

! number of dissociating ions in soluble components
no_ions=[3.0,0.0,0.0,2.0,0.0,0.0]
!
fracbcem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0]
fracocem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
mode=(mode_choice > 0)
component=.FALSE.
soluble=.FALSE.
DO imode=1,nmodes
  DO icp=1,ncp
    IF (((component_mode(imode,icp) == 1) .AND.                                &
          (component_choice(icp) == 1)) .AND.                                  &
          (mode_choice(imode) == 1)) THEN
      component(imode,icp)=.TRUE.
    END IF
    IF (soluble_choice(icp) == 1) soluble(icp)=.TRUE.
  END DO
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcocso_5mode
SUBROUTINE UKCA_MODE_DUonly_2MODE
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  only du cpt in 2 (insoluble) modes.
!  Uses  4 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

INTEGER :: imode
INTEGER :: icp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_DUONLY_2MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncp = 6
CALL ukca_mode_allocate_ctl_vars()
! Component names
component_names =                                                              &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
mode_choice=[0,0,0,0,0,1,1]
! Specify which modes are soluble
modesol=[1,1,1,1,0,0,0]
! Component switches (1=on, 0=off)
component_choice=[0,0,0,0,1,0] ! ***only dust on***
! Components that are soluble
soluble_choice=[1,0,0,1,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
component_mode(1,1:ncp)=[1,0,1,0,0,1]       !allowed in nuc_sol
component_mode(2,1:ncp)=[1,1,1,0,0,1]       !allowed in ait_sol
component_mode(3,1:ncp)=[1,1,1,1,1,1]       !allowed in acc_sol
component_mode(4,1:ncp)=[1,1,1,1,1,1]       !allowed in cor_sol
component_mode(5,1:ncp)=[0,1,1,0,0,0]       !allowed in ait_ins
component_mode(6,1:ncp)=[0,0,0,0,1,0]       !allowed in acc_ins
component_mode(7,1:ncp)=[0,0,0,0,1,0]       !allowed in cor_ins

! Specify size limits of geometric mean diameter for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]

! Specify fixed geometric standard deviation for each mode
sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0] ! to match M7, but sigacc=1.4

DO imode=1,nmodes
  x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
num_eps=[1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14]

! Initial fractions of mass in each mode among components
mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !nucln. soluble
mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Aitken soluble
mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !accum. soluble
mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0] !coarse soluble
mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0] !Aitken insoluble
mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !accum. insoluble
mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !coarse insoluble

! Molar masses of components (kg mol-1)
mm=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)
! Mass density of components (kg m^-3)
rhocomp = [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

! Adjust the density of BC?
IF (glomap_config%l_ukca_radaer) THEN
   i_tune_bc_loc = glomap_config%i_ukca_tune_bc
ELSE IF (l_glomap_clim_radaer) THEN
   i_tune_bc_loc =  i_glomap_clim_tune_bc
ELSE
   ! Do not tune BC density
   i_tune_bc_loc = 0 
END IF
SELECT CASE (i_tune_bc_loc)
   CASE (i_ukca_bc_tuned)
      rhocomp(cp_bc) = rho_bc_tuned
   CASE (i_ukca_bc_mg_mix)
      rhocomp(cp_bc) = rho_bc_mg_mix
END SELECT

! Adjust the density of NaCl (sea-salt)
IF (glomap_config%l_fix_nacl_density) rhocomp(cp_cl) = rho_nacl

DO imode=1,nmodes
  ddpmid(imode)=EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
  rhommav=0.0
  DO icp=1,ncp
    rhommav=rhommav+mfrac_0(imode,icp)*(rhocomp(icp)/mm(icp))
  END DO
  mmid(imode)=(pi/6.0)*(ddpmid (imode)**3)*(rhommav*avogadro)*x(imode)
  mlo (imode)=(pi/6.0)*(ddplim0(imode)**3)*(rhommav*avogadro)*x(imode)
  mhi (imode)=(pi/6.0)*(ddplim1(imode)**3)*(rhommav*avogadro)*x(imode)
END DO

! Modes resulting when two modes coagulate
coag_mode(1,1:nmodes)=[1,2,3,4,2,3,4]
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
coag_mode(2,1:nmodes)=[2,2,3,4,2,3,4]
! 2 coagulating with:    1,2,3,4,5,6,7
coag_mode(3,1:nmodes)=[3,3,3,4,3,3,4]
! 3 coagulating with:    1,2,3,4,5,6,7
coag_mode(4,1:nmodes)=[4,4,4,4,4,4,4]
! 4 coagulating with:    1,2,3,4,5,6,7
coag_mode(5,1:nmodes)=[2,2,3,4,5,6,7]
! 5 coagulating with:    1,2,3,4,5,6,7
coag_mode(6,1:nmodes)=[3,3,3,4,6,6,7]
! 6 coagulating with:    1,2,3,4,5,6,7
coag_mode(7,1:nmodes)=[4,4,4,4,7,7,7]
! 7 coagulating with:    1,2,3,4,5,6,7

! number of dissociating ions in soluble components
no_ions=[3.0,0.0,0.0,2.0,0.0,0.0]
!
fracbcem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0]
fracocem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
!
! Set logical variables
mode=(mode_choice > 0)
component=.FALSE.
soluble=.FALSE.
DO imode=1,nmodes
  DO icp=1,ncp
    IF (((component_mode(imode,icp) == 1) .AND.                                &
          (component_choice(icp) == 1)) .AND.                                  &
          (mode_choice(imode) == 1)) THEN
      component(imode,icp)=.TRUE.
    END IF
    IF (soluble_choice(icp) == 1) soluble(icp)=.TRUE.
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE UKCA_MODE_DUonly_2MODE

! ######################################################################

SUBROUTINE ukca_mode_sussbcocntnh_5mode_7cpt
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  sulfate, sea-salt, bc, oc (secondary & primary combined),
!  NH4, and NO3 in 5 modes, and 7 components.
!  Also adding NaNO3 - a coarse nitrate tracer in ACC/COR
!  soluble modes - 7 components
!  Uses 28 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

INTEGER :: imode
INTEGER :: icp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SUSSBCOCNTNH_5MODE_7CPT'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! No of components
ncp = 9
CALL ukca_mode_allocate_ctl_vars()
! Component names
component_names =                                                              &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org',          &
          'no3    ','nano3  ','nh4    ']

! Mode switches (1=on, 0=0ff)
mode_choice=[1,1,1,1,1,0,0]
! Specify which modes are soluble
modesol=[1,1,1,1,0,0,0]
! Component switches (1=on, 0=off)
component_choice=[1,1,1,1,0,0,1,1,1]
! *** n.b. only have h2so4, bc, oc, nacl, no3, nh4 cpts on for this setup ***
! Components that are soluble
soluble_choice=[1,0,0,1,0,0,1,1,1]
! Components allowed in each mode (must be consistent with coag_mode)
component_mode(1,1:ncp)=[1,0,1,0,0,1,0,0,0]       !allowed in nuc_sol
component_mode(2,1:ncp)=[1,1,1,0,0,1,1,0,1]       !allowed in ait_sol
component_mode(3,1:ncp)=[1,1,1,1,1,1,1,1,1]       !allowed in acc_sol
component_mode(4,1:ncp)=[1,1,1,1,1,1,1,1,1]       !allowed in cor_sol
component_mode(5,1:ncp)=[0,1,1,0,0,0,0,0,0]       !allowed in ait_ins
component_mode(6,1:ncp)=[0,0,0,0,1,0,0,0,0]       !allowed in acc_ins
component_mode(7,1:ncp)=[0,0,0,0,1,0,0,0,0]       !allowed in cor_ins

! Specify size limits of geometric mean diameter for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]

! Specify fixed geometric standard deviation for each mode
sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0] ! to match M7, but sigacc=1.4

DO imode=1,nmodes
  x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
END DO

! Specify threshold for ND (per cc) below which don't do calculations
num_eps=[1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14]

! Initial fractions of mass in each mode among components
mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] !nucln. soluble
mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] !Aitken soluble
mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] !accum. soluble
mfrac_0(4,1:ncp)=[0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0] !coarse soluble
mfrac_0(5,1:ncp)=[0.0,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0] !Aitken insoluble
mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0] !accum. insoluble
mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0] !coarse insoluble

! Molar masses of components (kg mol-1)
mm=[0.098, 0.012, 0.0168, 0.05844, 0.100, 0.0168, 0.062, 0.084, 0.018]
!    h2so4  bc     oc      nacl     dust   so      no3     nano3   nh4
! n.b. mm_bc=0.012, mm_oc=0.012*1.4=0.168 (1.4 POM:OC ratio)

! Mass density of components (kg m^-3)
rhocomp=[1769.0, 1500.0, 1500.0, 1600.0, 2650.0, 1500.0,                       &
          1500.0, 1600.0, 1769.0]
!         h2so4   bc      oc      nacl    dust    so
!         no3     nano3   nh4

! Adjust the density of BC?
IF (glomap_config%l_ukca_radaer) THEN
   i_tune_bc_loc = glomap_config%i_ukca_tune_bc
ELSE IF (l_glomap_clim_radaer) THEN
   i_tune_bc_loc =  i_glomap_clim_tune_bc
ELSE
   ! Do not tune BC density
   i_tune_bc_loc = 0 
END IF
SELECT CASE (i_tune_bc_loc)
   CASE (i_ukca_bc_tuned)
      rhocomp(cp_bc) = rho_bc_tuned
   CASE (i_ukca_bc_mg_mix)
      rhocomp(cp_bc) = rho_bc_mg_mix
END SELECT

! Adjust the density of NaCl (sea-salt)
IF (glomap_config%l_fix_nacl_density) THEN
  rhocomp(cp_cl) = rho_nacl
  rhocomp(cp_nn) = rho_nacl
END IF

DO imode=1,nmodes
  ddpmid(imode)=EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
  rhommav=0.0
  DO icp=1,ncp
    rhommav=rhommav+mfrac_0(imode,icp)*(rhocomp(icp)/mm(icp))
  END DO
  mmid(imode)=(pi/6.0)*(ddpmid (imode)**3)*(rhommav*avogadro)*x(imode)
  mlo (imode)=(pi/6.0)*(ddplim0(imode)**3)*(rhommav*avogadro)*x(imode)
  mhi (imode)=(pi/6.0)*(ddplim1(imode)**3)*(rhommav*avogadro)*x(imode)
END DO

! Modes resulting when two modes coagulate
coag_mode(1,1:nmodes)=[1,2,3,4,2,3,4]
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
coag_mode(2,1:nmodes)=[2,2,3,4,2,3,4]
! 2 coagulating with:    1,2,3,4,5,6,7
coag_mode(3,1:nmodes)=[3,3,3,4,3,3,4]
! 3 coagulating with:    1,2,3,4,5,6,7
coag_mode(4,1:nmodes)=[4,4,4,4,4,4,4]
! 4 coagulating with:    1,2,3,4,5,6,7
coag_mode(5,1:nmodes)=[2,2,3,4,5,6,7]
! 5 coagulating with:    1,2,3,4,5,6,7
coag_mode(6,1:nmodes)=[3,3,3,4,6,6,7]
! 6 coagulating with:    1,2,3,4,5,6,7
coag_mode(7,1:nmodes)=[4,4,4,4,7,7,7]
! 7 coagulating with:    1,2,3,4,5,6,7

! number of dissociating ions in soluble components
no_ions=[3.0,0.0,0.0,2.0,0.0,0.0,2.0,2.0,2.0]

! Fractions of primary BC/POM emissions to go to each mode at emission
! (emit into insoluble Aitken for this setup).
fracbcem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0]
fracocem=[0.0,0.0,0.0,0.0,1.0,0.0,0.0]

! Set logical variables
mode=(mode_choice > 0)
component=.FALSE.
soluble=.FALSE.
DO imode=1,nmodes
  DO icp=1,ncp
    IF (((component_mode(imode,icp) == 1) .AND.                                &
          (component_choice(icp) == 1)) .AND.                                  &
          (mode_choice(imode) == 1)) THEN
      component(imode,icp)=.TRUE.
    END IF
    IF (soluble_choice(icp) == 1) soluble(icp)=.TRUE.
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_sussbcocntnh_5mode_7cpt

SUBROUTINE ukca_mode_solinsol_6mode
! ---------------------------------------------------------------------|
!  Subroutine to define modes and components for version with
!  soluble and insoluble components in 6 modes.
!  Uses 12 aerosol tracers
! ---------------------------------------------------------------------|
IMPLICIT NONE

INTEGER :: imode
INTEGER :: icp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SOLINSOL_6MODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncp = 6
CALL ukca_mode_allocate_ctl_vars()
! Component names
component_names =                                                              &
        ['h2so4  ','bcarbon','ocarbon','nacl   ','dust   ','sec_org']

! Mode switches (1=on, 0=0ff)
mode_choice=[1,1,1,1,0,1,1]
! Specify which modes are soluble
modesol=[1,1,1,1,0,0,0]
! Component switches (1=on, 0=off)
component_choice=[1,0,0,0,1,0]
! *** n.b. only have soluble (as h2so4) and insoluble (as dust) in this setup
! Components that are soluble
soluble_choice=[1,0,0,0,0,0]
! Components allowed in each mode (must be consistent with coag_mode)
component_mode(1,1:ncp)=[1,0,0,0,0,0]       !allowed in nuc_sol
component_mode(2,1:ncp)=[1,0,0,0,0,0]       !allowed in ait_sol
component_mode(3,1:ncp)=[1,0,0,0,0,0]       !allowed in acc_sol
component_mode(4,1:ncp)=[1,0,0,0,0,0]       !allowed in cor_sol
component_mode(5,1:ncp)=[0,0,0,0,0,0]       !allowed in ait_ins
component_mode(6,1:ncp)=[0,0,0,0,1,0]       !allowed in acc_ins
component_mode(7,1:ncp)=[0,0,0,0,1,0]       !allowed in cor_ins

! Specify size limits of geometric mean diameter for each mode
!  set ddplim34 here to be 500nm to match value found by bin-mode comparison
ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]

! Specify fixed geometric standard deviation for each mode
sigmag=[1.59,1.59,1.40,2.0,1.59,1.59,2.0] ! to match M7, but sigacc=1.4

DO imode=1,nmodes
  x(imode)=EXP(4.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
END DO
!
! Specify threshold for ND (per cc) below which don't do calculations
num_eps=[1.0e-8,1.0e-8,1.0e-8,1.0e-14,1.0e-8,1.0e-14,1.0e-14]

! Initial fractions of mass in each mode among components
mfrac_0(1,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !nucln. soluble
mfrac_0(2,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !Aitken soluble
mfrac_0(3,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !accum. soluble
mfrac_0(4,1:ncp)=[1.0,0.0,0.0,0.0,0.0,0.0] !coarse soluble
mfrac_0(5,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !Aitken insoluble
mfrac_0(6,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !accum. insoluble
mfrac_0(7,1:ncp)=[0.0,0.0,0.0,0.0,1.0,0.0] !coarse insoluble

! Molar masses of components (kg mol-1)
mm=[0.098,0.012,0.0168,0.05844,0.100,0.0168]
!          h2so4  bc     oc    nacl   dust    so
! n.b. mm_bc=0.012, mm_oc=mm_so=0.012*1.4=0.168 (1.4 POM:OC ratio)
! Mass density of components (kg m^-3)
rhocomp = [1769.0,1500.0,1500.0,1600.0,2650.0,1500.0]

IF (glomap_config%l_fix_nacl_density) rhocomp(cp_cl) = rho_nacl

DO imode=1,nmodes
  ddpmid(imode)=EXP(0.5*(LOG(ddplim0(imode))+LOG(ddplim1(imode))))
  rhommav=0.0
  DO icp=1,ncp
   rhommav=rhommav+mfrac_0(imode,icp)*(rhocomp(icp)/mm(icp))
  END DO
  mmid(imode)=(pi/6.0)*(ddpmid (imode)**3)*(rhommav*avogadro)*x(imode)
  mlo (imode)=(pi/6.0)*(ddplim0(imode)**3)*(rhommav*avogadro)*x(imode)
  mhi (imode)=(pi/6.0)*(ddplim1(imode)**3)*(rhommav*avogadro)*x(imode)
END DO

! Modes resulting when two modes coagulate
coag_mode(1,1:nmodes)=[1,2,3,4,2,3,4]
! 1 coagulating with:    1,2,3,4,5,6,7 produces...
coag_mode(2,1:nmodes)=[2,2,3,4,2,3,4]
! 2 coagulating with:    1,2,3,4,5,6,7
coag_mode(3,1:nmodes)=[3,3,3,4,3,3,4]
! 3 coagulating with:    1,2,3,4,5,6,7
coag_mode(4,1:nmodes)=[4,4,4,4,4,4,4]
! 4 coagulating with:    1,2,3,4,5,6,7
coag_mode(5,1:nmodes)=[2,2,3,4,5,6,7]
! 5 coagulating with:    1,2,3,4,5,6,7
coag_mode(6,1:nmodes)=[6,6,6,6,6,6,7]
! 6 coagulating with:    1,2,3,4,5,6,7
coag_mode(7,1:nmodes)=[7,7,7,7,7,7,7]
! 7 coagulating with:    1,2,3,4,5,6,7

! number of dissociating ions in soluble components
no_ions=[3.0,0.0,0.0,0.0,0.0,0.0]
!
fracbcem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
fracocem=[0.0,1.0,0.0,0.0,0.0,0.0,0.0]
! fractions of primary BC/POM emissions to go to each mode at emission
! (emit into soluble Aitken for this setup).
!
! Set logical variables
mode=(mode_choice > 0)
component=.FALSE.
soluble=.FALSE.
DO imode=1,nmodes
  DO icp=1,ncp
    IF (((component_mode(imode,icp) == 1) .AND.                                &
          (component_choice(icp) == 1)) .AND.                                  &
          (mode_choice(imode) == 1)) THEN
      component(imode,icp)=.TRUE.
    END IF
    IF (soluble_choice(icp) == 1) soluble(icp)=.TRUE.
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_solinsol_6mode

! ############################################################################

SUBROUTINE ukca_mode_allocate_ctl_vars
! ---------------------------------------------------------------------|
!  Subroutine to allocate arrays in this module
! ---------------------------------------------------------------------|

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_ALLOCATE_CTL_VARS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(component_choice)) ALLOCATE(component_choice(ncp))
IF (.NOT. ALLOCATED(soluble_choice)) ALLOCATE(soluble_choice(ncp))
IF (.NOT. ALLOCATED(mm)) ALLOCATE(mm(ncp))
IF (.NOT. ALLOCATED(rhocomp)) ALLOCATE(rhocomp(ncp))
IF (.NOT. ALLOCATED(no_ions)) ALLOCATE(no_ions(ncp))
IF (.NOT. ALLOCATED(component_names)) ALLOCATE(component_names(ncp))
IF (.NOT. ALLOCATED(soluble)) ALLOCATE(soluble(ncp))
IF (.NOT. ALLOCATED(component_mode)) ALLOCATE(component_mode(nmodes,ncp))
IF (.NOT. ALLOCATED(mfrac_0)) ALLOCATE(mfrac_0(nmodes,ncp))
IF (.NOT. ALLOCATED(component)) ALLOCATE(component(nmodes,ncp))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_allocate_ctl_vars

END MODULE ukca_mode_setup
