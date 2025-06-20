! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module providing subroutines to put GLOMAP-mode diagnostics into
!   the STASHwork array
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds,
! University of Oxford, and the Met. Office.
! See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v8 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_mode_diags_mod

USE ukca_config_constants_mod, ONLY: avogadro, boltzmann
USE ukca_constants,         ONLY: mmw, m_air
USE ukca_mode_setup,        ONLY: nmodes, ncp_max,                             &
                                  cp_su, cp_bc, cp_oc, cp_cl, cp_du,           &
                                  cp_no3, cp_nh4, cp_nn, cp_mp
USE ukca_config_specification_mod, ONLY: glomap_variables
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod,            ONLY: ereport
USE umPrintMgr,             ONLY: umPrint, umMessage
USE parkind1,               ONLY: jprb, jpim      ! for Dr Hook tracing
USE yomhook,                ONLY: lhook, dr_hook  ! for Dr Hook tracing

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_MODE_DIAGS_MOD'

LOGICAL, SAVE, PUBLIC :: l_ukca_cmip6_diags=.FALSE.
! Set to true to enable CMIP6 diagnostics. This is done by the parent via
! direct use of this module pending the addition of full diagnostic support to
! the UKCA API

LOGICAL, SAVE, PUBLIC :: l_ukca_pm_diags=.FALSE.
! Set to true to enable PM10 or PM2.5 38 diagnostics. This is done by the
! parent via direct use of this module pending the addition of full diagnostic
! support to the UKCA API

REAL, ALLOCATABLE, PUBLIC :: mdwat_diag(:,:)
! Molecular concentration of water in each mode (molecules per particle)

REAL, ALLOCATABLE, PUBLIC :: wetdp_diag(:,:)
! Geometric mean wet diameter of particles in each mode (m)

PUBLIC :: ukca_mode_diags_alloc
PUBLIC :: ukca_mode_diags

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE ukca_mode_diags_alloc(nbox)
! Description:
!   Allocate arrays for copies of water content and wet diameter fields
!   that are provided by the aerosol scheme and required as input to
!   diagnostic calculations in UKCA_MODE_DIAGS.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Number of elements
INTEGER, INTENT(IN) :: nbox

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_DIAGS_ALLOC'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate diagnostic field for copy of aerosol water array if required
IF (.NOT. ALLOCATED(mdwat_diag) .AND.                                          &
    (l_ukca_cmip6_diags .OR. l_ukca_pm_diags))                                 &
  ALLOCATE(mdwat_diag(nbox,nmodes))
! Allocate diagnostic field for copy of wet particle diameter array if required
IF (.NOT. ALLOCATED(wetdp_diag) .AND. l_ukca_pm_diags)                         &
  ALLOCATE(wetdp_diag(nbox,nmodes))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_mode_diags_alloc


! ----------------------------------------------------------------------
SUBROUTINE ukca_mode_diags(row_length, rows, model_levels,                     &
                           nbox, n_mode_tracers,                               &
                           p_theta_levels,                                     &
                           t_theta_levels, mode_tracers,                       &
                           interf_z, len_stashwork38, stashwork38)
! Description:
!   Obtain number densities and component material concentrations for
!   each mode to use in calculating diagnostics, derive the required
!   diagnostics and put them in the STASHwork array.
! ----------------------------------------------------------------------

USE ukca_mode_tracer_maps_mod, ONLY: nmr_index, mmr_index
USE asad_mod,                  ONLY: jpctr
USE ukca_types_mod,            ONLY: log_small

IMPLICIT NONE

! UKCA domain dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! Number of elements
INTEGER, INTENT(IN) :: nbox

! Number of MODE tracers
INTEGER, INTENT(IN) :: n_mode_tracers

! Pressure on theta levels
REAL, INTENT(IN)    :: p_theta_levels(row_length, rows, model_levels)

! Temperature on theta levels
REAL, INTENT(IN)    :: t_theta_levels(row_length, rows, model_levels)

! MODE tracer array
REAL, INTENT(IN)    :: mode_tracers(row_length, rows, model_levels,            &
                                    n_mode_tracers)

! Height of interface levels above surface (m)
REAL, INTENT(IN)    :: interf_z(row_length, rows, 0:model_levels)

! Length of diagnostics array
INTEGER, INTENT(IN) :: len_stashwork38

! Work array for STASH
REAL, INTENT(IN OUT) :: stashwork38(len_stashwork38)

! Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: mfrac_0(:,:)
REAL,    POINTER :: mm(:)
REAL,    POINTER :: mmid(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp

INTEGER :: ifirst         ! index of first mode tracer in nmr_index, mmr_index
INTEGER :: imode          ! loop counter for modes
INTEGER :: icp            ! loop counter for components
INTEGER :: itra           ! tracer number
INTEGER :: icode          ! error code
REAL    :: pmid(nbox)     ! Air pressure at mid-point (Pa)
REAL    :: tmid(nbox)     ! Temperature at mid-point (K)
REAL    :: tr_rs(nbox)    ! Local variable to hold re-shaped aerosol tracers
REAL    :: aird(nbox)     ! Number density of air (cm^-3)
REAL    :: nd(nbox,nmodes)! Aerosol particle number density for mode (cm^-3)
REAL    :: md(nbox,nmodes,glomap_variables%ncp)
                          ! Average component concentration of aerosol particle
                          ! in mode (molecules.particle^-1)
LOGICAL (KIND=log_small) :: mask1(nbox)              ! To mask negatives
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_DIAGS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mfrac_0     => glomap_variables%mfrac_0
mm          => glomap_variables%mm
mmid        => glomap_variables%mmid
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp

icode = 0

! -------------------------------------------------------------
! Check for availability of diagnostic fields required as input
! -------------------------------------------------------------

IF (.NOT. ALLOCATED(mdwat_diag)) THEN
  icode = 1
  cmessage = ' Aerosol water content diagnostics not found'
  WRITE(umMessage,'(A70)') cmessage
  CALL umPrint(umMessage,src=RoutineName)
  CALL ereport(RoutineName,icode,cmessage)
END IF

IF (.NOT. ALLOCATED(wetdp_diag) .AND. l_ukca_pm_diags) THEN
  icode = 3
  cmessage = ' Wet particle diameter diagnostics not found'
  WRITE(umMessage,'(A70)') cmessage
  CALL umPrint(umMessage,src=RoutineName)
  CALL ereport(RoutineName,icode,cmessage)
END IF

! ----------------------------------------------
! Obtain ND and MD arrays from mode tracer array
! ----------------------------------------------

tmid(:)=RESHAPE(t_theta_levels(:,:,:),[nbox])
pmid(:)=RESHAPE(p_theta_levels(:,:,:),[nbox])

! calculate molecular concentration of air (/cm3)
aird(:)=pmid(:)/(tmid(:)*boltzmann*1.0e6)

!  Find index of 1st mode tracer, as nmr_index and
!   mmr_index index all ukca tracers
ifirst = jpctr + 1
DO imode=1,nmodes
  IF (mode(imode)) THEN
    itra = nmr_index(imode) - ifirst + 1
    tr_rs(:)=RESHAPE(mode_tracers(:,:,:,itra),[nbox])
    mask1(:)=(tr_rs(:) < 0.0)
    WHERE (mask1(:))
      tr_rs(:)=0.0
    END WHERE
    ! .. above sets tr_rs to zero if negative
    nd(:,imode)=tr_rs(:)*aird(:)
    ! .. above sets ND (particles per cc) from advected number-mixing-ratio

    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        itra = mmr_index(imode,icp) - ifirst + 1
        tr_rs(:)=RESHAPE(mode_tracers(:,:,:,itra),[nbox])
        mask1(:)=(tr_rs(:) < 0.0)
        WHERE (mask1(:))
          tr_rs(:)=0.0
        END WHERE
        ! .. above sets tr_rs to zero if negative
        !        mask1(:)=(nd(:,imode) > num_eps(imode))
        mask1(:)=(nd(:,imode) > 1e-30)
        WHERE (mask1(:))
          md(:,imode,icp)=(m_air/mm(icp))*aird(:)*tr_rs(:)/nd(:,imode)
        ELSE WHERE
          md(:,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
        END WHERE
        ! above sets MD (molecules per particle) from advected mass-mix-ratio
        ! .. note that only "trusts" values where ND>NUM_EPS
      ELSE
        md(:,imode,icp)=0.0
      END IF
    END DO    ! loop over cpts

  END IF    ! mode

END DO  ! loop over modes

! --------------------------------------------------------------------
! Calculate required diagnostics and copy them to the STASHwork array
! --------------------------------------------------------------------

IF (l_ukca_cmip6_diags) THEN
  CALL mode_diags_cmip6(row_length, rows, model_levels,                        &
                        nbox, aird, nd, md, interf_z,                          &
                        len_stashwork38, stashwork38)
END IF

IF (l_ukca_pm_diags) THEN
  CALL mode_diags_pm(row_length, rows, model_levels, nbox, nd, md,             &
                     len_stashwork38, stashwork38)
END IF

! --------------------------------------------------------------------
! Deallocate arrays
! --------------------------------------------------------------------
IF (ALLOCATED(wetdp_diag)) DEALLOCATE(wetdp_diag)
IF (ALLOCATED(mdwat_diag)) DEALLOCATE(mdwat_diag)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_mode_diags


! ----------------------------------------------------------------------
SUBROUTINE mode_diags_cmip6(row_length, rows, model_levels, nbox, aird,        &
                            nd, md, interf_z, len_stashwork38, stashwork38)
! Description:
!   Calculate the required CMIP6 diagnostics and copy them to the
!   STASHwork array.
! ----------------------------------------------------------------------

USE ukca_um_legacy_mod, ONLY:                                                  &
  stashcode_so4_nuc_sol,  stashcode_so4_ait_sol,  stashcode_so4_acc_sol,       &
  stashcode_so4_cor_sol,  stashcode_bc_ait_sol,   stashcode_bc_acc_sol,        &
  stashcode_bc_cor_sol,   stashcode_bc_ait_insol, stashcode_oc_nuc_sol,        &
  stashcode_oc_nuc_sol,   stashcode_oc_ait_sol,   stashcode_oc_acc_sol,        &
  stashcode_oc_cor_sol,   stashcode_oc_ait_insol, stashcode_ss_acc_sol,        &
  stashcode_ss_cor_sol,   stashcode_du_acc_sol,   stashcode_du_cor_sol,        &
  stashcode_du_acc_insol, stashcode_du_cor_insol, stashcode_du_sup_insol,      &
  stashcode_n_nuc_sol,    stashcode_n_ait_sol,    stashcode_n_acc_sol,         &
  stashcode_n_cor_sol,    stashcode_n_ait_insol,  stashcode_n_acc_insol,       &
  stashcode_n_cor_insol,  stashcode_n_sup_insol,                               &
  stashcode_h2o_nuc_sol,  stashcode_h2o_ait_sol,  stashcode_h2o_acc_sol,       &
  stashcode_h2o_cor_sol,  stashcode_h2o_total,    stashcode_so4_nuc_sol_load,  &
  stashcode_so4_ait_sol_load,  stashcode_so4_acc_sol_load,                     &
  stashcode_so4_cor_sol_load,  stashcode_so4_total_load,                       &
  stashcode_bc_ait_sol_load,   stashcode_bc_acc_sol_load,                      &
  stashcode_bc_cor_sol_load,   stashcode_bc_ait_insol_load,                    &
  stashcode_bc_total_load,     stashcode_oc_nuc_sol_load,                      &
  stashcode_oc_ait_sol_load,   stashcode_oc_acc_sol_load,                      &
  stashcode_oc_cor_sol_load,   stashcode_oc_ait_insol_load,                    &
  stashcode_oc_total_load,     stashcode_ss_acc_sol_load,                      &
  stashcode_ss_cor_sol_load,   stashcode_ss_total_load,                        &
  stashcode_du_acc_sol_load,   stashcode_du_cor_sol_load,                      &
  stashcode_du_acc_insol_load, stashcode_du_cor_insol_load,                    &
  stashcode_du_sup_insol_load, stashcode_du_total_load,                        &
  stashcode_h2o_nuc_sol_load,  stashcode_h2o_ait_sol_load,                     &
  stashcode_h2o_acc_sol_load,  stashcode_h2o_cor_sol_load,                     &
  stashcode_h2o_total_load,    stashcode_h2o_mmr,                              &
  stashcode_nh4_ait_sol, stashcode_nh4_acc_sol,                                &
  stashcode_nh4_cor_sol, stashcode_no3_ait_sol,                                &
  stashcode_no3_acc_sol, stashcode_no3_cor_sol,                                &
  stashcode_nn_acc_sol , stashcode_nn_cor_sol ,                                &
  stashcode_nh4_ait_sol_load,                                                  &
  stashcode_nh4_acc_sol_load, stashcode_nh4_cor_sol_load,                      &
  stashcode_no3_ait_sol_load,                                                  &
  stashcode_no3_acc_sol_load, stashcode_no3_cor_sol_load,                      &
  stashcode_nn_acc_sol_load , stashcode_nn_cor_sol_load ,                      &
  stashcode_nh4_total_load, stashcode_no3_total_load, stashcode_nn_total_load, &
  stashcode_mp_ait_sol,   stashcode_mp_acc_sol,   stashcode_mp_cor_sol,        &
  stashcode_mp_ait_insol, stashcode_mp_acc_insol, stashcode_mp_cor_insol,      &
  stashcode_mp_sup_insol, stashcode_mp_sup_insol_load,                         &
  stashcode_mp_ait_sol_load,   stashcode_mp_acc_sol_load,                      &
  stashcode_mp_cor_sol_load,   stashcode_mp_ait_insol_load,                    &
  stashcode_mp_acc_insol_load, stashcode_mp_cor_insol_load,                    &
  stashcode_mp_total_load,                                                     &
  len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf,         &
  si_last, stashcode_glomap_sec,                                               &
  copydiag, copydiag_3d

IMPLICIT NONE

! UKCA domain dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! Number of elements
INTEGER, INTENT(IN) :: nbox

! Number density of air (cm^-3)
REAL, INTENT(IN)    :: aird(nbox)

! Aerosol particle number density for mode (cm^-3)
REAL, INTENT(IN)    :: nd(nbox,nmodes)

! Average component concentration of aerosol particle in mode
! (molecules.particle^-1)
REAL, INTENT(IN)    :: md(nbox,nmodes,glomap_variables%ncp)

! Height of interface levels above surface (m)
REAL, INTENT(IN)    :: interf_z(row_length, rows, 0:model_levels)

! Diagnostics array
INTEGER, INTENT(IN) :: len_stashwork38

! Work array for STASH
REAL, INTENT(IN OUT) :: stashwork38(len_stashwork38)

! Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: mm (:)
LOGICAL, POINTER :: mode (:)
INTEGER, POINTER :: modesol(:)
INTEGER, POINTER :: ncp

! Table to locate STASH item numbers from mode and component
! Final components are unused (SO, NO3, NH4, Na, Cl) currently
INTEGER, PARAMETER :: item_component_cmip6(nmodes,ncp_max) =                   &
  RESHAPE( [                                                                   &
! SO4
  stashcode_so4_nuc_sol, stashcode_so4_ait_sol, stashcode_so4_acc_sol,         &
  stashcode_so4_cor_sol, -1,                    -1,                            &
  -1                   , -1,                                                   &
! BC
  -1,                    stashcode_bc_ait_sol,  stashcode_bc_acc_sol,          &
  stashcode_bc_cor_sol,  stashcode_bc_ait_insol,-1,                            &
  -1                   , -1,                                                   &
! OC
  stashcode_oc_nuc_sol,  stashcode_oc_ait_sol,  stashcode_oc_acc_sol,          &
  stashcode_oc_cor_sol,  stashcode_oc_ait_insol,-1,                            &
  -1                  ,  -1,                                                   &
! NaCl
  -1,                    -1,                    stashcode_ss_acc_sol,          &
  stashcode_ss_cor_sol,  -1,                    -1,                            &
  -1                  ,  -1,                                                   &
! Dust
  -1,                    -1,                    stashcode_du_acc_sol,          &
  stashcode_du_cor_sol,  -1,                    stashcode_du_acc_insol,        &
  stashcode_du_cor_insol, stashcode_du_sup_insol,                              &
! Secondary organic
  -1,                    -1,                    -1,                            &
  -1,                    -1,                    -1,                            &
  -1,                    -1,                                                   &
! NO3
  -1,                    stashcode_no3_ait_sol, stashcode_no3_acc_sol,         &
  stashcode_no3_cor_sol,                    -1,                    -1,         &
  -1                   , -1,                                                   &
! NaNO3
  -1,                                       -1,  stashcode_nn_acc_sol,         &
   stashcode_nn_cor_sol,                    -1,                    -1,         &
  -1                   , -1,                                                   &
! NH4
  -1,                    stashcode_nh4_ait_sol, stashcode_nh4_acc_sol,         &
  stashcode_nh4_cor_sol,                    -1,                    -1,         &
  -1                   , -1,                                                   &
  ! Microplastics
  -1,                    stashcode_mp_ait_sol,  stashcode_mp_acc_sol,          &
  stashcode_mp_cor_sol,  stashcode_mp_ait_insol,stashcode_mp_acc_insol,        &
  stashcode_mp_cor_insol,stashcode_mp_sup_insol                                &
  ], [nmodes, ncp_max])

! Table to locate STASH item numbers for mode and component loads
! Final components are unused (SO,NO3,NH4,Na,Cl) currently
INTEGER, PARAMETER :: item_load_cmip6(nmodes,ncp_max) =                        &
  RESHAPE( [                                                                   &
! SO4
  stashcode_so4_nuc_sol_load, stashcode_so4_ait_sol_load,                      &
  stashcode_so4_acc_sol_load, stashcode_so4_cor_sol_load,                      &
  -1,                         -1,                                              &
  -1,                         -1,                                              &
! BC
  -1,                          stashcode_bc_ait_sol_load,                      &
  stashcode_bc_acc_sol_load,   stashcode_bc_cor_sol_load,                      &
  stashcode_bc_ait_insol_load, -1,                                             &
  -1,                          -1,                                             &
! OC
  stashcode_oc_nuc_sol_load  , stashcode_oc_ait_sol_load,                      &
  stashcode_oc_acc_sol_load  , stashcode_oc_cor_sol_load,                      &
  stashcode_oc_ait_insol_load, -1,                                             &
  -1,                          -1,                                             &
! NaCl
  -1,                          -1,                                             &
  stashcode_ss_acc_sol_load,   stashcode_ss_cor_sol_load,                      &
  -1,                          -1,                                             &
  -1,                          -1,                                             &
! Dust
  -1,                          -1,                                             &
  stashcode_du_acc_sol_load,   stashcode_du_cor_sol_load,                      &
  -1,                          stashcode_du_acc_insol_load,                    &
  stashcode_du_cor_insol_load, stashcode_du_sup_insol_load,                    &
! SO
  -1,                         -1,                                              &
  -1,                         -1,                                              &
  -1,                         -1,                                              &
  -1,                         -1,                                              &
!  NO3
  -1,                         stashcode_no3_ait_sol_load,                      &
  stashcode_no3_acc_sol_load, stashcode_no3_cor_sol_load,                      &
  -1,                         -1,                                              &
  -1,                         -1,                                              &
!  NaNO3
  -1,                                                 -1,                      &
   stashcode_nn_acc_sol_load,  stashcode_nn_cor_sol_load,                      &
  -1,                         -1,                                              &
  -1,                         -1,                                              &
!  NH4
  -1,                         stashcode_nh4_ait_sol_load,                      &
  stashcode_nh4_acc_sol_load, stashcode_nh4_cor_sol_load,                      &
  -1,                         -1,                                              &
  -1,                         -1,                                              &
! Microplastics
  -1,                          stashcode_mp_ait_sol_load,                      &
  stashcode_mp_acc_sol_load,   stashcode_mp_cor_sol_load,                      &
  stashcode_mp_ait_insol_load, stashcode_mp_acc_insol_load,                    &
  stashcode_mp_cor_insol_load, stashcode_mp_sup_insol_load                     &
  ], [nmodes, ncp_max])

! Table to locate STASH item numbers for aerosol water density
INTEGER :: item_water_cmip6(nmodes) = [                                        &
  stashcode_h2o_nuc_sol, stashcode_h2o_ait_sol, stashcode_h2o_acc_sol,         &
  stashcode_h2o_cor_sol, -1,                    -1,                            &
  -1,                    -1  ]

! Table to locate STASH item numbers for aerosol water loads
INTEGER :: item_water_load_cmip6(nmodes) = [                                   &
  stashcode_h2o_nuc_sol_load, stashcode_h2o_ait_sol_load,                      &
  stashcode_h2o_acc_sol_load, stashcode_h2o_cor_sol_load,                      &
  -1,                         -1,                                              &
  -1,                         -1  ]

! Table to locate STASH item numbers for number density
INTEGER :: item_number_cmip6(nmodes) =                                         &
 [stashcode_n_nuc_sol,   stashcode_n_ait_sol,   stashcode_n_acc_sol,           &
   stashcode_n_cor_sol,   stashcode_n_ait_insol, stashcode_n_acc_insol,        &
   stashcode_n_cor_insol, stashcode_n_sup_insol]

! Table to locate STASH item numbers for total aerosol loads
! Positions for SO, NO3, and NH4 unused currently
INTEGER, PARAMETER :: item_aerosol_load_cmip6(ncp_max+1) = [                   &
  stashcode_so4_total_load, stashcode_bc_total_load,                           &
  stashcode_oc_total_load,  stashcode_ss_total_load,                           &
  stashcode_du_total_load,  -1,                                                &
  stashcode_no3_total_load, stashcode_nh4_total_load,                          &
  stashcode_nn_total_load, stashcode_mp_total_load,                            &
  stashcode_h2o_total_load  ]

INTEGER :: n_soluble      ! No of soluble modes
INTEGER :: section        ! stash section
INTEGER :: tsection       ! stash section * 1000
INTEGER :: item           ! stash item
INTEGER :: imode          ! loop counter for modes
INTEGER :: icp            ! loop counter for components
INTEGER :: icp2           ! loop counter for H2O load
INTEGER :: icode          ! error code
INTEGER :: im_index       ! internal model index
INTEGER :: k              ! loop counter
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message

REAL, PARAMETER :: cm3_per_m3 = 1.0e6       ! cm^3 per m^3

REAL    :: field(nbox)                      ! Output field (1D)
REAL    :: field3d(row_length,rows,model_levels)   ! Output field (3D)
REAL    :: dz(row_length,rows,model_levels) ! Depth of each layer (m)

REAL    :: aerosol_component_density(row_length,rows,model_levels,             &
                                     nmodes,glomap_variables%ncp+1)
                                            ! aerosol cpnt. density (kg.m^-3)
REAL    :: aerosol_number_density(row_length,rows,model_levels,nmodes)
                                            ! aerosol no. density (m^-3)
REAL    :: aerosol_component_load(row_length,rows,nmodes,                      &
                                            glomap_variables%ncp+1)
                                            ! integrated aerosol load (kg.m^-2)
REAL    :: aerosol_total_load(row_length,rows,                                 &
                                            glomap_variables%ncp+1)
                                            ! total aerosol load of each
                                            ! component + h2o (kg.m^-2)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MODE_DIAGS_CMIP6'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mm          => glomap_variables%mm
mode        => glomap_variables%mode
modesol     => glomap_variables%modesol
ncp         => glomap_variables%ncp

im_index = 1
icode = 0

!  -------------------------------------------------------
!  CMIP6 Diagnostics for aerosol mass and number densities
!  -------------------------------------------------------

section = stashcode_glomap_sec
tsection = stashcode_glomap_sec*1000
n_soluble = SUM(modesol)

! -----------------------------------------------------
! Calculate aerosol component and number density arrays
! -----------------------------------------------------

aerosol_component_density = 0.0
aerosol_number_density = 0.0
field = 0.0

DO k=1,model_levels
  dz(:,:,k) = interf_z(:,:,k) - interf_z(:,:,k-1)
END DO

DO imode=1,nmodes
  IF (mode(imode)) THEN

    ! Number densities
    field(:) = nd(:,imode)*cm3_per_m3   ! Number per m^3
    aerosol_number_density(:,:,:,imode) = RESHAPE(field(:),                    &
                                             [row_length,rows,model_levels])
    ! Component densities
    DO icp = 1,ncp
      IF (component(imode,icp)) THEN
        field(:) = mm(icp)*md(:,imode,icp)*nd(:,imode)*cm3_per_m3/avogadro
        aerosol_component_density(:,:,:,imode,icp) = RESHAPE(field(:),         &
                                             [row_length,rows,model_levels])

      END IF
    END DO     ! icp
  END IF     ! mode

  !  H2O component density
  IF (modesol(imode) == 1) THEN
    field(:) = mmw*mdwat_diag(:,imode)*nd(:,imode)*cm3_per_m3/avogadro
    aerosol_component_density(:,:,:,imode,ncp+1) = RESHAPE(field(:),           &
                                           [row_length,rows,model_levels])
  END IF
END DO  ! imode


! Two-dimensional aerosol loads (kg.m^-2)
! ---------------------------------------

aerosol_component_load = 0.0
DO imode=1,nmodes
  IF (mode(imode)) THEN
    DO icp = 1,ncp
      IF (component(imode,icp)) THEN
        item = item_load_cmip6(imode,icp)   ! find correct stash item no.
        item = item - tsection
        IF (item < 0) THEN
          icode = 11
          cmessage = ' No valid item number identified for this mode'//        &
                     ' and component'
          WRITE(umMessage,'(A70,2I6)') cmessage,imode,icp
          CALL umPrint(umMessage,src=RoutineName)
          CALL ereport(RoutineName,icode,cmessage)
        END IF
        IF (sf(item,section)) THEN
          DO k = 1,model_levels
            aerosol_component_load(:,:,imode,icp) =                            &
                  aerosol_component_load(:,:,imode,icp) +                      &
                  aerosol_component_density(:,:,k,imode,icp) * dz(:,:,k)
          END DO
          CALL copydiag(stashwork38(si(item,section,im_index):                 &
               si_last(item,section,im_index)),                                &
               aerosol_component_load(:,:,imode,icp),row_length,rows)
        END IF
      END IF
    END DO     ! icp
    ! H2O
    IF (modesol(imode) == 1) THEN
      item = item_water_load_cmip6(imode)   ! find correct stash item no.
      item = item - tsection
      IF (item < 0) THEN
        icode = 12
        cmessage = ' No valid item number identified for this mode'//          &
                   ' and component'
        WRITE(umMessage,'(A70,2I6)') cmessage,imode,icp
        CALL umPrint(umMessage,src=RoutineName)
        CALL ereport(RoutineName,icode,cmessage)
      END IF
      IF (sf(item,section)) THEN
        DO k = 1,model_levels
          aerosol_component_load(:,:,imode,ncp+1) =                            &
                aerosol_component_load(:,:,imode,ncp+1) +                      &
                aerosol_component_density(:,:,k,imode,ncp+1) * dz(:,:,k)
        END DO
        CALL copydiag(stashwork38(si(item,section,im_index):                   &
             si_last(item,section,im_index)),                                  &
             aerosol_component_load(:,:,imode,icp),row_length,rows)
      END IF
    END IF     ! mode
  END IF     ! modesol
END DO    ! imode

! Total aerosol loads in each component
! -------------------------------------

aerosol_total_load = 0.0
DO imode=1,nmodes
  IF (mode(imode)) THEN
    DO icp = 1,ncp
      IF (component(imode,icp)) THEN
        DO k = 1,model_levels
          aerosol_total_load(:,:,icp) = aerosol_total_load(:,:,icp) +          &
            aerosol_component_density(:,:,k,imode,icp) * dz(:,:,k)
        END DO
      END IF
    END DO  ! icp
    ! H2O
    DO k = 1,model_levels
      aerosol_total_load(:,:,ncp+1) = aerosol_total_load(:,:,ncp+1) +          &
        aerosol_component_density(:,:,k,imode,ncp+1) * dz(:,:,k)
    END DO
  END IF
END DO   ! imode

! Total aerosol loads of components
DO icp=1,ncp
  IF (ANY(component(:,icp))) THEN
    item = item_aerosol_load_cmip6(icp)
    item = item - tsection
    IF (item < 0) THEN
      icode = 13
      cmessage = ' No valid item number identified for this component'
      WRITE(umMessage,'(A70,I6)') cmessage,icp
      CALL umPrint(umMessage,src=RoutineName)
      CALL ereport(RoutineName,icode,cmessage)
    END IF
    IF (sf(item,section)) THEN
      CALL copydiag(stashwork38(si(item,section,im_index):                     &
           si_last(item,section,im_index)),                                    &
           aerosol_total_load(:,:,icp),row_length,rows)
    END IF
  END IF
END DO   ! icp

! H2O load
icp = ncp_max + 1
item = item_aerosol_load_cmip6(icp)
item = item - tsection
IF (item < 0) THEN
  icode = 14
  cmessage = ' No valid item number identified for this component'
  WRITE(umMessage,'(A70,I6)') cmessage,icp
  CALL umPrint(umMessage,src=RoutineName)
  CALL ereport(RoutineName,icode,cmessage)
END IF

icp2 = ncp+1
IF (sf(item,section)) THEN
  CALL copydiag(stashwork38(si(item,section,im_index):                         &
       si_last(item,section,im_index)),                                        &
       aerosol_total_load(:,:,icp2),row_length,rows)
END IF

! Copy items into STASHwork array:
! --------------------------------

! 1) aerosol components in kg/m^3
! -------------------------------
DO imode=1,nmodes
  IF (mode(imode)) THEN
    DO icp = 1,ncp
      IF (component(imode,icp)) THEN
        item = item_component_cmip6(imode,icp)   ! find correct stash item no.
        item = item - tsection
        IF (item < 0) THEN
          icode = 15
          cmessage = ' No valid item number identified for this mode'//        &
                     ' and component'
          WRITE(umMessage,'(A70,2I6)') cmessage,imode,icp
          CALL umPrint(umMessage,src=RoutineName)
          CALL ereport(RoutineName,icode,cmessage)
        END IF

        IF (sf(item,section)) THEN
          CALL copydiag_3d(stashwork38(si(item,section,im_index):              &
            si_last(item,section,im_index)),                                   &
            aerosol_component_density(:,:,:,imode,icp),                        &
            row_length,rows,model_levels,                                      &
            stlist(:,stindex(1,item,section,im_index)),len_stlist,             &
            stash_levels,num_stash_levels+1)
        END IF       ! sf
      END IF        ! component
    END DO         ! icp
  END IF          ! mode(imode)
END DO           ! imode

! 2) aerosol water in kg/m^3
! --------------------------
DO imode=1,nmodes
  IF (mode(imode) .AND. modesol(imode) == 1) THEN
    item = item_water_cmip6(imode)   ! find correct stash item no.
    item = item - tsection
    IF (item < 0) THEN
      icode = 16
      cmessage = ' No valid item number identified for this mode'
      WRITE(umMessage,'(A70,2I6)') cmessage,imode,ncp+1
      CALL umPrint(umMessage,src=RoutineName)
      CALL ereport(RoutineName,icode,cmessage)
    END IF

    IF (sf(item,section)) THEN
      CALL copydiag_3d(stashwork38(si(item,section,im_index):                  &
        si_last(item,section,im_index)),                                       &
        aerosol_component_density(:,:,:,imode,ncp+1),                          &
        row_length,rows,model_levels,                                          &
        stlist(:,stindex(1,item,section,im_index)),len_stlist,                 &
        stash_levels,num_stash_levels+1)
    END IF       ! sf
  END IF       ! mode(imode)
END DO       ! imode

! Total aerosol water from all soluble modes
! ------------------------------------------
item = stashcode_h2o_total - tsection
IF (sf(item,section)) THEN
  field3d(:,:,:) = 0.0
  DO imode=1,nmodes
    IF (mode(imode) .AND. modesol(imode) == 1) THEN
      field3d(:,:,:) = field3d(:,:,:) +                                        &
                       aerosol_component_density(:,:,:,imode,ncp+1)
    END IF
  END DO

  CALL copydiag_3d(stashwork38(si(item,section,im_index):                      &
    si_last(item,section,im_index)),                                           &
    field3d(:,:,:),                                                            &
    row_length,rows,model_levels,                                              &
    stlist(:,stindex(1,item,section,im_index)),len_stlist,                     &
    stash_levels,num_stash_levels+1)
END IF   ! sf

! 3) aerosol number density in m^-3
! ---------------------------------
DO imode=1,nmodes
  IF (mode(imode)) THEN
    item = item_number_cmip6(imode)
    item = item - tsection
    IF (item < 0) THEN
      icode = 17
      cmessage = ' No valid item number identified for this mode'
      WRITE(umMessage,'(A70,I6)') cmessage,imode
      CALL umPrint(umMessage,src=RoutineName)
      CALL ereport(RoutineName,icode,cmessage)
    END IF

    IF (sf(item,section)) THEN
      CALL copydiag_3d(stashwork38(si(item,section,im_index):                  &
        si_last(item,section,im_index)),                                       &
        aerosol_number_density(:,:,:,imode),                                   &
        row_length,rows,model_levels,                                          &
        stlist(:,stindex(1,item,section,im_index)),len_stlist,                 &
        stash_levels,num_stash_levels+1)
    END IF       ! sf
  END IF       ! mode(imode)
END DO       ! imode

! 4) Aerosol water mass mixing ratio
! ----------------------------------

item = stashcode_h2o_mmr - tsection
IF (sf(item,section)) THEN
  field(:) = 0.0
  DO imode = 1,n_soluble
    ! aerosol water over all soluble modes
    field(:) = field(:) + mdwat_diag(:,imode)*nd(:,imode)
  END DO
  field(:) = field(:)*mmw/(aird(:)*m_air)  ! kg/kg
  field3d = RESHAPE(field,[row_length,rows,model_levels])

  CALL copydiag_3d(stashwork38(si(item,section,im_index):                      &
    si_last(item,section,im_index)),                                           &
    field3d,                                                                   &
    row_length,rows,model_levels,                                              &
    stlist(:,stindex(1,item,section,im_index)),len_stlist,                     &
    stash_levels,num_stash_levels+1)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE mode_diags_cmip6


! ----------------------------------------------------------------------
SUBROUTINE mode_diags_pm(row_length, rows, model_levels, nbox, nd, md,         &
                         len_stashwork38, stashwork38)
! Description:
!   Calculate the required PM10 and PM2.5 diagnostics and copy them to
!   the STASHwork array.
! ----------------------------------------------------------------------

USE ukca_pm_diags_mod,    ONLY: ukca_pm_diags, pm_request_struct
USE ukca_um_legacy_mod,   ONLY: stashcode_glomap_sec,                          &
                                stashcode_pm10_dry, stashcode_pm2p5_dry,       &
                                stashcode_pm10_wet, stashcode_pm2p5_wet,       &
                                stashcode_pm10_so4, stashcode_pm10_bc,         &
                                stashcode_pm10_oc, stashcode_pm10_ss,          &
                                stashcode_pm10_du,                             &
                                stashcode_pm2p5_so4, stashcode_pm2p5_bc,       &
                                stashcode_pm2p5_oc, stashcode_pm2p5_ss,        &
                                stashcode_pm2p5_du,                            &
                                stashcode_pm10_nh4, stashcode_pm2p5_nh4,       &
                                stashcode_pm10_no3, stashcode_pm2p5_no3,       &
                                stashcode_pm10_nn , stashcode_pm2p5_nn,        &
                                stashcode_pm10_mp , stashcode_pm2p5_mp,        &
                                len_stlist, stindex,                           &
                                stlist, num_stash_levels,                      &
                                stash_levels, si, sf, si_last,                 &
                                copydiag_3d

IMPLICIT NONE

! UKCA domain dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! Number of elements
INTEGER, INTENT(IN) :: nbox

! Aerosol particle number density for mode (cm^-3)
REAL, INTENT(IN)    :: nd(nbox,nmodes)

! Average component concentration of aerosol particle in mode
! (molecules.particle^-1)
REAL, INTENT(IN)    :: md(nbox,nmodes,glomap_variables%ncp)

! Length of diagnostics array
INTEGER, INTENT(IN) :: len_stashwork38

! Work array for STASH
REAL, INTENT(IN OUT) :: stashwork38(len_stashwork38)

! Local variables

INTEGER, PARAMETER :: n_size_cat = 2 ! Number of PM size categories
INTEGER, PARAMETER :: ndiagcp = 9    ! Number of component diagnostics
                                     ! available for each PM size category
INTEGER :: im_index       ! Internal model index
INTEGER :: section        ! Stash section
INTEGER :: tsection       ! Stash section * 1000
INTEGER :: item           ! Stash item
INTEGER :: i_size_cat     ! Loop counter for PM size category
INTEGER :: i_size_cat_max ! Highest PM size category index to be used
INTEGER :: icp            ! Loop counter for components
INTEGER :: icp_max        ! Highest component index to be used for component
                          ! contributions
INTEGER :: i              ! Loop counter
INTEGER :: icode          ! Error code
REAL, PARAMETER :: d_cutoff(n_size_cat) = [10.0e-6, 2.5e-6]
                          ! Size limits for particulate matter (m)
REAL    :: field3d(row_length,rows,model_levels)   ! Output field (3D)

! Specification for PM STASH items

TYPE :: pm_item_struct
  INTEGER :: stcode_total_dry(n_size_cat) ! STASH code of total PM dry mass
                                          !   (1=PM10, 2=PM2.5)
  INTEGER :: stcode_total_wet(n_size_cat) ! STASH code of total PM wet mass
  INTEGER :: stcode_component(ndiagcp,n_size_cat)
                                          ! STASH code of contribution to PM
  INTEGER :: i_ref_component(ndiagcp)     ! Reference component number
END TYPE pm_item_struct

TYPE(pm_item_struct), PARAMETER :: pm_item = pm_item_struct(                   &
  [stashcode_pm10_dry, stashcode_pm2p5_dry],                                   &
  [stashcode_pm10_wet, stashcode_pm2p5_wet],                                   &
  RESHAPE([stashcode_pm10_so4, stashcode_pm10_bc, stashcode_pm10_oc,           &
            stashcode_pm10_ss, stashcode_pm10_du, stashcode_pm10_no3,          &
            stashcode_pm10_nn, stashcode_pm10_nh4, stashcode_pm10_mp,          &
            stashcode_pm2p5_so4, stashcode_pm2p5_bc, stashcode_pm2p5_oc,       &
            stashcode_pm2p5_ss, stashcode_pm2p5_du, stashcode_pm2p5_no3,       &
            stashcode_pm2p5_nn, stashcode_pm2p5_nh4, stashcode_pm2p5_mp],      &
            [ndiagcp,2]),                                                      &
  [cp_su, cp_bc, cp_oc, cp_cl, cp_du, cp_no3, cp_nn, cp_nh4, cp_mp])

! Request flags indicating the PM diagnostics to be calculated
! in 'ukca_pm_diags'
TYPE(pm_request_struct) :: pm_request

! PM diagnostics from 'ukca_pm_diags' (ug m-3)
REAL, ALLOCATABLE :: pm_dry(:,:)         ! Total PM dry mass by size category
REAL, ALLOCATABLE :: pm_wet(:,:)         ! Total PM wet mass by size category
REAL, ALLOCATABLE :: pm_component(:,:,:) ! Component contributions to PM by
                                         ! component number and size category

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MODE_DIAGS_PM'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

section = stashcode_glomap_sec
tsection = stashcode_glomap_sec*1000

! Initialize PM request flags

ALLOCATE(pm_request%l_total_dry(n_size_cat))
ALLOCATE(pm_request%l_total_wet(n_size_cat))
ALLOCATE(pm_request%l_component(glomap_variables%ncp,n_size_cat))
pm_request%l_total_dry(:) = .FALSE.
pm_request%l_total_wet(:) = .FALSE.
pm_request%l_component(:,:) = .FALSE.

! Determine which PM diagnostics must be calculated based on STASH requests,
! set PM request flags accordingly and allocate storage (don't allocate above
! maximum index to be used)

! Total PM dry mass
i_size_cat_max = 0
DO i_size_cat = 1,n_size_cat
  item = pm_item%stcode_total_dry(i_size_cat) - tsection
  IF (sf(item,section)) THEN
    pm_request%l_total_dry(i_size_cat) = .TRUE.
    i_size_cat_max = i_size_cat
  END IF
END DO
IF (i_size_cat_max > 0) THEN
  ALLOCATE(pm_dry(nbox,i_size_cat_max))
ELSE
  ALLOCATE(pm_dry(1,1))
END IF

! Total PM wet mass
i_size_cat_max = 0
DO i_size_cat = 1,n_size_cat
  item = pm_item%stcode_total_wet(i_size_cat) - tsection
  IF (sf(item,section)) THEN
    pm_request%l_total_wet(i_size_cat) = .TRUE.
    i_size_cat_max = i_size_cat
  END IF
END DO
IF (i_size_cat_max > 0) THEN
  ALLOCATE(pm_wet(nbox,i_size_cat_max))
ELSE
  ALLOCATE(pm_wet(1,1))
END IF

! Component contribution requests.
i_size_cat_max = 0
icp_max = 0
DO i_size_cat = 1,n_size_cat
  DO i = 1,ndiagcp
    item = pm_item%stcode_component(i,i_size_cat) - tsection
    IF (sf(item,section)) THEN
      icp = pm_item%i_ref_component(i)
      pm_request%l_component(icp,i_size_cat) = .TRUE.
      i_size_cat_max = i_size_cat
      IF (icp > icp_max) icp_max = icp
    END IF
  END DO
END DO
IF (icp > 0) THEN
  ALLOCATE(pm_component(nbox,icp_max,i_size_cat_max))
ELSE
  ALLOCATE(pm_component(1,1,1))
END IF

! Derive the PM diagnostics needed for producing the required STASH items

CALL ukca_pm_diags(nbox,nd,md,mdwat_diag,wetdp_diag,d_cutoff,pm_request,       &
                   pm_dry,pm_wet,pm_component)

! Copy required items to STASH work array

im_index = 1
icode = 0

DO i_size_cat = 1,n_size_cat

  ! Total dry mass for PM size category
  item = pm_item%stcode_total_dry(i_size_cat) - tsection
  IF (sf(item,section)) THEN
    field3d = RESHAPE(pm_dry(:,i_size_cat),                                    &
                      [row_length,rows,model_levels])

    CALL copydiag_3d(stashwork38(si(item,section,im_index):                    &
      si_last(item,section,im_index)),                                         &
      field3d,                                                                 &
      row_length,rows,model_levels,                                            &
      stlist(:,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1)
  END IF

  ! Total wet mass for PM size category
  item = pm_item%stcode_total_wet(i_size_cat) - tsection
  IF (sf(item,section)) THEN
    field3d = RESHAPE(pm_wet(:,i_size_cat),                                    &
                      [row_length,rows,model_levels])

    CALL copydiag_3d(stashwork38(si(item,section,im_index):                    &
      si_last(item,section,im_index)),                                         &
      field3d,                                                                 &
      row_length,rows,model_levels,                                            &
      stlist(:,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1)
  END IF

  ! Component contributions to PM size category
  DO i = 1,ndiagcp
    item = pm_item%stcode_component(i,i_size_cat) - tsection
    IF (sf(item,section)) THEN
      icp = pm_item%i_ref_component(i)
      field3d = RESHAPE(pm_component(:,icp,i_size_cat),                        &
                        [row_length,rows,model_levels])

      CALL copydiag_3d(stashwork38(si(item,section,im_index):                  &
        si_last(item,section,im_index)),                                       &
        field3d,                                                               &
        row_length,rows,model_levels,                                          &
        stlist(:,stindex(1,item,section,im_index)),len_stlist,                 &
        stash_levels,num_stash_levels+1)
    END IF
  END DO

END DO

DEALLOCATE(pm_request%l_total_dry)
DEALLOCATE(pm_request%l_total_wet)
DEALLOCATE(pm_request%l_component)
DEALLOCATE(pm_dry)
DEALLOCATE(pm_wet)
DEALLOCATE(pm_component)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE mode_diags_pm

END MODULE ukca_mode_diags_mod
