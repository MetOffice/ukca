! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module holding procedures mainly for UKCA online emissions processing
!
! Method:
!   Declares and sets up the emissions data structure for online emissions.
!   Description of each subroutine contained here:
!    1. UKCA_ONL_EMISS_INIT: Identify online emissions (including those from
!         coupling components and initialise emissions structure accordingly.
!    2. UKCA_EM_STRUCT_INIT: Initialise emissions structure to
!         default values.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_emiss_mod
!

USE ukca_emiss_struct_mod, ONLY: ukca_em_struct, ukca_em_struct_init,          &
                                 ncdf_emissions
USE ukca_config_specification_mod, ONLY: ukca_config, glomap_config,           &
                                         i_light_param_off
USE ukca_config_defs_mod,  ONLY: em_chem_spec
USE ukca_emiss_mode_mod,   ONLY: ukca_def_mode_emiss
USE ukca_mode_setup,       ONLY: nmodes, rhocomp, mode, component,             &
                                 cp_cl, cp_du, cp_su, cp_oc, cp_bc,            &
                                 cp_no3, cp_nh4, cp_nn,                        &
                                 moment_number, moment_mass

! UM specific modules
USE ukca_um_legacy_mod,    ONLY: pi, rgas => r
USE ereport_mod,           ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,            ONLY: umPrint, umMessage, PrintStatus,              &
                                 PrStatus_Normal, newline
USE chemistry_constants_mod, ONLY: avogadro, boltzmann

USE parkind1,              ONLY: jpim, jprb      ! DrHook
USE yomhook,               ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Default private
PRIVATE

! Emission Namelist Items
INTEGER         :: num_onln_em_flds  ! Number of online emission fields
INTEGER, PUBLIC :: num_em_flds       ! Number of total emission fields
INTEGER, PUBLIC :: num_cdf_em_flds   ! Number of file-based emission fields
INTEGER         :: num_seasalt_modes ! Number of modes used by seasalt
INTEGER         :: num_pmoc_modes    ! Number of modes used by primary marine OC
INTEGER         :: num_dust_modes    ! Number of modes used by dust
INTEGER         :: num_iso2_modes    ! Number of modes used by interactive SO2
INTEGER         :: num_ioc_modes     ! Number of modes used by interactive OC
INTEGER         :: num_ibc_modes     ! Number of modes used by interactive BC
INTEGER         :: num_no3_modes     ! Number of modes used by nitrate
INTEGER         :: num_nh4_modes     ! Number of modes used by ammonium
INTEGER         :: num_nn_modes      ! Number of modes used by NaNO3

! Indices for any online emissions which are not input as NetCDF
INTEGER, PUBLIC :: inox_light          ! Index for lightning emiss of NOx
INTEGER, PUBLIC :: ich4_wetl           ! Index for wetland emiss of CH4
INTEGER, PUBLIC :: ic5h8_ibvoc         ! Index for iBVOC emiss of isoprene
INTEGER, PUBLIC :: ic10h16_ibvoc       ! Index for iBVOC emiss of monoterpenes
INTEGER, PUBLIC :: ich3oh_ibvoc        ! Index for iBVOC emiss of methanol
INTEGER, PUBLIC :: ich3coch3_ibvoc     ! Index for iBVOC emiss of acetone
INTEGER, PUBLIC :: idms_seaflux        ! Index for marine DMS emissions
INTEGER, PUBLIC :: iseasalt_first      ! Index for seasalt emiss 1st mode
INTEGER, PUBLIC :: ipmoc_first         ! Index for PMOC emiss 1st mode
INTEGER, PUBLIC :: idust_first         ! Index for dust emiss 1st mode
INTEGER, PUBLIC :: ino3_first          ! Index for nitrate emiss 1st mode
INTEGER, PUBLIC :: inh4_first          ! Index for ammonium emiss 1st mode
INTEGER, PUBLIC :: iseasalt_hno3       ! Index for seasalt lost to nitrate
INTEGER, PUBLIC :: idust_hno3          ! Index for dust lost to nitrate
INTEGER, PUBLIC :: inano3_cond         ! Index for nitrate condensing on ss
INTEGER, PUBLIC :: ico_inferno         ! Index for INFERNO emiss of CO
INTEGER, PUBLIC :: ich4_inferno        ! Index for INFERNO emiss of CH4
INTEGER, PUBLIC :: inox_inferno        ! Index for INFERNO emiss of NOx
INTEGER, PUBLIC :: iso2_inferno        ! Index for INFERNO emiss of SO2
INTEGER, PUBLIC :: ioc_inferno         ! Index for INFERNO emiss of OC
INTEGER, PUBLIC :: ibc_inferno         ! Index for INFERNO emiss of BC

! Super array of emissions
TYPE (ukca_em_struct), ALLOCATABLE, PUBLIC :: emissions (:)

! Subroutines available outside this module
PUBLIC :: ukca_onl_emiss_init, ukca_emiss_update_mode, ukca_emiss_mode_map,    &
          ukca_emiss_spatial_vars_init

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EMISS_MOD'

CONTAINS

! ---------------------------------------------------------------------
! Description:
!  Check if any emissions fields are being calculated online, or obtained from
!  coupling components and initialise the emissions structure.
!
! ---------------------------------------------------------------------
SUBROUTINE ukca_onl_emiss_init (row_length, rows,  model_levels)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT(IN) :: row_length       ! model dimensions
INTEGER,           INTENT(IN) :: rows
INTEGER,           INTENT(IN) :: model_levels

! Local variables
INTEGER :: ecount          ! Index in emissions array

CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER            :: ierror    ! Error code

! GLOMAP related variables
INTEGER :: moment_step = moment_mass - moment_number
INTEGER :: imode
INTEGER :: imoment

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_ONL_EMISS_INIT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierror   = 0

! Set this from size of ncdf_emissions -assuming that structure has been filled
! before call to this routine.
num_cdf_em_flds = MAX(0, SIZE(ncdf_emissions))

num_onln_em_flds = 0
! Initialise indices for interative BVOC emissions
ic5h8_ibvoc     = -99
ic10h16_ibvoc   = -99
ich3oh_ibvoc    = -99
ich3coch3_ibvoc = -99

! Initialise indices for INFERNO emissions
ico_inferno     = -99
ich4_inferno    = -99
inox_inferno    = -99
iso2_inferno    = -99
ioc_inferno     = -99
ibc_inferno     = -99

! -----------------------------------------------------------------------------
! The total nr of emission fields will be equal to the number of fields
! found in the NetCDF files (i.e. num_cdf_em_flds) plus any other fields
! that account for online emissions.
!
! Add first an index for online NOx emiss from lightning if present.
IF (ukca_config%i_ukca_light_param /= i_light_param_off) THEN
  inox_light = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 1
END IF
!
! Add index for online CH4 from wetland emissions if these are
! active and prescribed surface CH4 concentrations are not used.
! Also increment the index for total number of emissions
IF (ukca_config%l_ukca_qch4inter .AND.                                         &
    (.NOT. ukca_config%l_ukca_prescribech4)) THEN
  ich4_wetl        = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 1
ELSE
  ich4_wetl = -99
END IF

! Interactive emissions of biogenic volatile organic compounds (iBVOC)
IF (ukca_config%l_ukca_ibvoc) THEN

  IF (ANY (em_chem_spec == 'C5H8      ')) THEN
    ic5h8_ibvoc      = num_cdf_em_flds + num_onln_em_flds + 1
    num_onln_em_flds = num_onln_em_flds + 1
  ELSE
    ! Report a warning to indicate that it will not be used
    ! because species not emitted in the chemistry scheme
    ierror      = -2
    cmessage    = 'Biogenic emission field of C5H8 '        //                 &
                  'ignored for this chemistry scheme'
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF

  IF (ANY (em_chem_spec == 'Monoterp  ')) THEN
    ic10h16_ibvoc    = num_cdf_em_flds  + num_onln_em_flds + 1
    num_onln_em_flds = num_onln_em_flds + 1
  ELSE
    ierror        = -3
    cmessage      = 'Biogenic emission field of Monoterp '  //                 &
                    'ignored for this chemistry scheme'
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF

  ! Methanol and Acetone emissions, although available from JULES are
  ! currently DEACTIVATED (index set to -99). This is due to a lack
  ! of reliable checks for double-counting (prescribed via files as
  ! well as online). It is also estimated that the contribution from
  ! vegetation (i.e. online) is small compared to other e.g. industrial
  ! sources (prescribed via files).

  ! Note that methanol can be referred to with two different names
  ! depending on the chemistry scheme.
  IF (ANY (em_chem_spec == 'MeOH      ') .OR.                                  &
      ANY (em_chem_spec == 'CH3OH     ')) THEN
    !     ich3oh_ibvoc     = num_cdf_em_flds  + num_onln_em_flds + 1
    !     num_onln_em_flds = num_onln_em_flds + 1
    ich3oh_ibvoc  = -99
  ELSE
    ierror        = -4
    cmessage      = 'Biogenic emission field of methanol '  //                 &
                    'ignored for this chemistry scheme'
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF

  IF (ANY (em_chem_spec == 'Me2CO     ')) THEN
    !     ich3coch3_ibvoc  = num_cdf_em_flds  + num_onln_em_flds + 1
    !     num_onln_em_flds = num_onln_em_flds + 1
    ich3coch3_ibvoc  = -99
  ELSE
    ierror          = -5
    cmessage        = 'Biogenic emission field of Me2CO '   //                 &
                      'ignored for this chemistry scheme'
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF
END IF

! Interactive emissions of fire emissions from INFERNO
IF (ukca_config%l_ukca_inferno) THEN
  cmessage='WARNING: INFERNO fire interactive emissions have'//       newline//&
           'been requested. Please ensure fire emissions are not'//   newline//&
           'prescribed through ancillary files.'
  ierror = -1
  CALL ereport(RoutineName, ierror, cmessage)

  IF (ANY (em_chem_spec == 'CO        ')) THEN
    ico_inferno      = num_cdf_em_flds  + num_onln_em_flds + 1
    num_onln_em_flds = num_onln_em_flds + 1
  ELSE
    ! Report a warning to indicate that it will not be used
    ! because species not emitted in the chemistry scheme
    ierror      = -2
    cmessage    = 'INFERNO fire emissions of CO '         //                   &
                  'ignored for this chemistry scheme'
    CALL ereport(RoutineName, ierror, cmessage)
  END IF

  IF (ukca_config%l_ukca_inferno_ch4 .AND.                                     &
      (.NOT. ukca_config%l_ukca_prescribech4)) THEN
    IF (ANY (em_chem_spec == 'CH4       ')) THEN
      ich4_inferno     = num_cdf_em_flds  + num_onln_em_flds + 1
      num_onln_em_flds = num_onln_em_flds + 1
    ELSE
      ! Report a warning to indicate that it will not be used
      ! because species not emitted in the chemistry scheme
      ierror      = -2
      cmessage    = 'INFERNO fire emissions of CH4 '        //                 &
                    'ignored for this chemistry scheme'
      CALL ereport(RoutineName, ierror, cmessage)
    END IF
  END IF

  IF (ANY (em_chem_spec == 'NO        ')) THEN
    inox_inferno     = num_cdf_em_flds  + num_onln_em_flds + 1
    num_onln_em_flds = num_onln_em_flds + 1
  ELSE
    ! Report a warning to indicate that it will not be used
    ! because species not emitted in the chemistry scheme
    ierror      = -2
    cmessage    = 'INFERNO fire emissions of NOx '        //                   &
                  'ignored for this chemistry scheme'
    CALL ereport(RoutineName, ierror, cmessage)
  END IF

  IF (ANY (em_chem_spec == 'SO2_nat   ')) THEN
    num_iso2_modes = COUNT(component(:,cp_su))
    iso2_inferno = num_cdf_em_flds + num_onln_em_flds + 1
    num_onln_em_flds = num_onln_em_flds + 2*num_iso2_modes
  ELSE
    ! Report a warning to indicate that it will not be used
    ! because species not emitted in the chemistry scheme
    ierror      = -2
    cmessage    = 'INFERNO fire emissions of SO2 '        //                   &
                  'ignored for this chemistry scheme'
    CALL ereport(RoutineName, ierror, cmessage)
  END IF

  IF (ANY (em_chem_spec == 'OM_biomass')) THEN
    num_ioc_modes = COUNT(component(:,cp_oc))
    ioc_inferno = num_cdf_em_flds + num_onln_em_flds + 1
    num_onln_em_flds = num_onln_em_flds + 2*num_ioc_modes
  ELSE
    ! Report a warning to indicate that it will not be used
    ! because species not emitted in the chemistry scheme
    ierror      = -2
    cmessage    = 'INFERNO fire emissions of OC '         //                   &
                  'ignored for this chemistry scheme'
    CALL ereport(RoutineName, ierror, cmessage)
  END IF

  IF (ANY (em_chem_spec == 'BC_biomass')) THEN
    num_ibc_modes = COUNT(component(:,cp_bc))
    ibc_inferno = num_cdf_em_flds + num_onln_em_flds + 1
    num_onln_em_flds = num_onln_em_flds + 2*num_ibc_modes
  ELSE
    ! Report a warning to indicate that it will not be used
    ! because species not emitted in the chemistry scheme
    ierror      = -2
    cmessage    = 'INFERNO fire emissions of BC '         //                   &
                  'ignored for this chemistry scheme'
    CALL ereport(RoutineName, ierror, cmessage)
  END IF

END IF !(l_ukca_inferno)

! Add index for online marine DMS emissions
IF (ukca_config%l_seawater_dms) THEN
  idms_seaflux = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 1
ELSE
  idms_seaflux = -99
END IF

! Add index for online seasalt emissions
! Make space for mass and number for every seasalt mode
IF (glomap_config%l_ukca_primss) THEN
  num_seasalt_modes = COUNT(component(:,cp_cl))
  iseasalt_first = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 2*num_seasalt_modes
ELSE
  iseasalt_first = -99
END IF
! Add index for online marine primary organic emissions
! Make space for number and mass for all OC modes
IF (glomap_config%l_ukca_prim_moc) THEN
  IF (.NOT. glomap_config%l_ukca_primss) THEN
    cmessage = 'Prim. marine OC requires sea-salt'          //                 &
               'emissions to be turned on'
    ierror = 1
    CALL ereport('UKCA_EMISS_INIT',ierror,cmessage)
  END IF
  num_pmoc_modes = COUNT(component(:,cp_oc))
  ipmoc_first = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 2*num_pmoc_modes
ELSE
  ipmoc_first = -99
END IF
! Add index for online dust emissions
IF (glomap_config%l_ukca_primdu) THEN
  num_dust_modes = COUNT(component(:,cp_du))
  idust_first = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 2*num_dust_modes
ELSE
  idust_first = -99
END IF

! Add indices for online nitrate/ammonium emissions
IF (glomap_config%l_ukca_fine_no3_prod) THEN
  num_no3_modes = COUNT(component(:,cp_no3))
  ino3_first = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 2*num_no3_modes
  num_nh4_modes = COUNT(component(:,cp_nh4))
  inh4_first = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 2*num_nh4_modes
ELSE
  ino3_first = -99
  inh4_first = -99
END IF
! Add indices for online coarse nitrate emissions
IF (glomap_config%l_ukca_coarse_no3_prod) THEN
  num_seasalt_modes = COUNT(component(:,cp_cl))
  iseasalt_hno3 = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 2*num_seasalt_modes
  num_nn_modes = COUNT(component(:,cp_nn))
  inano3_cond = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 2*num_nn_modes
ELSE
  iseasalt_hno3 = -99
  inano3_cond = -99
END IF

IF ((glomap_config%l_ukca_coarse_no3_prod) .AND.                               &
    (glomap_config%l_ukca_primdu)) THEN
  num_dust_modes = COUNT(component(:,cp_du))
  idust_hno3 = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 2*num_dust_modes
ELSE
  idust_hno3 = -99
END IF

num_em_flds = num_cdf_em_flds + num_onln_em_flds
IF (PrintStatus >= PrStatus_Normal) THEN
  WRITE (umMessage,'(A,1X,I3)') 'Number of online emission fields:',           &
                                 num_onln_em_flds
  CALL umPrint(umMessage,src=RoutineName)
  WRITE (umMessage,'(A,1X,I3)') 'Total number of emission fields:',            &
                                 num_em_flds
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Allocate the emission structure to hold NetCDF as well as Online emissions
ALLOCATE (emissions(num_em_flds))

CALL ukca_em_struct_init(emissions)

! Initialise fields for lightning NOx in the emissions structure
! (Note that there is no need to update file_name, update_freq
! and update_type for online emissions)
IF (ukca_config%i_ukca_light_param /= i_light_param_off) THEN
  emissions(inox_light)%var_name    = 'NOx_lightning_emissions'
  emissions(inox_light)%tracer_name = 'NO_lightng'
  emissions(inox_light)%l_online    = .TRUE.

  ! Units should be consistent with the emission field we will take from
  ! the call to ukca_light_ctl within ukca_new_emiss_ctl, i.e.
  ! 'kg(N) gridbox-1 s-1'. However, to follow the same approach as
  ! for NetCDF emission fields, we do not indicate 'N' in the units
  ! but use the substring 'expressed as nitrogen' in the long_name.
  emissions(inox_light)%units       = 'kg gridbox-1 s-1'
  emissions(inox_light)%lng_name    = 'tendency of atmosphere mass content' // &
    ' of nitrogen monoxide expressed as nitrogen due to emission from lightning'

  emissions(inox_light)%hourly_fact = 'none'         ! No need to apply
  emissions(inox_light)%daily_fact  = 'none'         ! time profiles

  emissions(inox_light)%vert_fact   = 'all_levels'   ! 3D emissions
  emissions(inox_light)%three_dim   = .TRUE.
  emissions(inox_light)%lowest_lev  = 1
  emissions(inox_light)%highest_lev = model_levels

  ! Emission values and diagostics allocated and initialised here.
  ! They will be updated in ukca_new_emiss_ctl and ukca_add_emiss,
  ! respectively
  ALLOCATE (emissions(inox_light)%values(row_length, rows, model_levels))
  ALLOCATE (emissions(inox_light)%diags (row_length, rows, model_levels))
  emissions(inox_light)%values(:,:,:) = 0.0
  emissions(inox_light)%diags (:,:,:) = 0.0
END IF

! Initialise fields for online wetland emiss of CH4 in emissions structure
! (done in a similar way as done above for lightning NOx, but only if
! these emissions exist)
IF (ukca_config%l_ukca_qch4inter .AND.                                         &
    (.NOT. ukca_config%l_ukca_prescribech4)) THEN
  emissions(ich4_wetl)%var_name    = 'CH4_wetland_emissions'
  emissions(ich4_wetl)%tracer_name = 'CH4_wetlnd'
  emissions(ich4_wetl)%l_online    = .TRUE.

  !  Units should be consistent with CH4 wetland emiss from JULES,
  !  given in 'kg(C) m-2 s-1'. Similarly as seen above for 'NO_lightng',
  !  we do not indicate 'C' in the units but use the substring
  !  'expressed as carbon' in the long_name.
  emissions(ich4_wetl)%units      = 'kg m-2 s-1'
  emissions(ich4_wetl)%lng_name   = 'tendency of atmosphere mass content' //   &
     ' of methane expressed as carbon due to emission from wetlands'

  emissions(ich4_wetl)%hourly_fact = 'none'      ! No need to apply
  emissions(ich4_wetl)%daily_fact  = 'none'      ! time profiles

  emissions(ich4_wetl)%vert_fact   = 'surface'   ! surface emissions
  emissions(ich4_wetl)%three_dim   = .FALSE.
  emissions(ich4_wetl)%lowest_lev  = 1
  emissions(ich4_wetl)%highest_lev = 1

  ALLOCATE (emissions(ich4_wetl)%values(row_length, rows, 1))
  ALLOCATE (emissions(ich4_wetl)%diags (row_length, rows, 1))
  emissions(ich4_wetl)%values(:,:,:) = 0.0
  emissions(ich4_wetl)%diags (:,:,:) = 0.0
END IF

! Initialise fields for iBVOC in emissions structure
IF (ukca_config%l_ukca_ibvoc) THEN
  IF (ic5h8_ibvoc > 0) THEN
    emissions(ic5h8_ibvoc)%var_name    = 'interactive_C5H8'
    emissions(ic5h8_ibvoc)%tracer_name = 'C5H8'
    emissions(ic5h8_ibvoc)%l_online    = .TRUE.

    ! iBVOC emissions from JULES are in 'kg(C) m-2 s-1'. Don't indicate this
    ! in %units but use substring 'expressed as carbon' in %lng_name.
    emissions(ic5h8_ibvoc)%units      = 'kg m-2 s-1'
    emissions(ic5h8_ibvoc)%lng_name   =                                        &
     'tendency of atmosphere mass content'                               //    &
     ' of isoprene expressed as carbon due to emission from vegetation'

    emissions(ic5h8_ibvoc)%hourly_fact = 'none'      ! No need to apply
    emissions(ic5h8_ibvoc)%daily_fact  = 'none'      ! time profiles

    emissions(ic5h8_ibvoc)%vert_fact   = 'surface'   ! surface emissions
    emissions(ic5h8_ibvoc)%three_dim   = .FALSE.
    emissions(ic5h8_ibvoc)%lowest_lev  = 1
    emissions(ic5h8_ibvoc)%highest_lev = 1

    ALLOCATE (emissions(ic5h8_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ic5h8_ibvoc)%diags (row_length, rows, 1))
    emissions(ic5h8_ibvoc)%values(:,:,:) = 0.0
    emissions(ic5h8_ibvoc)%diags (:,:,:) = 0.0
  END IF

  IF (ic10h16_ibvoc > 0) THEN
    emissions(ic10h16_ibvoc)%var_name    = 'interactive_terpene'
    emissions(ic10h16_ibvoc)%tracer_name = 'Monoterp'
    emissions(ic10h16_ibvoc)%l_online    = .TRUE.

    ! iBVOC emissions from JULES are in 'kg(C) m-2 s-1'. Don't indicate this
    ! in %units but use substring 'expressed as carbon' in %lng_name.
    emissions(ic10h16_ibvoc)%units      = 'kg m-2 s-1'
    emissions(ic10h16_ibvoc)%lng_name   =                                      &
     'tendency of atmosphere mass content'                               //    &
     ' of terpenes expressed as carbon due to emission from vegetation'

    emissions(ic10h16_ibvoc)%hourly_fact = 'none'      ! No need to apply
    emissions(ic10h16_ibvoc)%daily_fact  = 'none'      ! time profiles

    emissions(ic10h16_ibvoc)%vert_fact   = 'surface'   ! surface emissions
    emissions(ic10h16_ibvoc)%three_dim   = .FALSE.
    emissions(ic10h16_ibvoc)%lowest_lev  = 1
    emissions(ic10h16_ibvoc)%highest_lev = 1

    ALLOCATE (emissions(ic10h16_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ic10h16_ibvoc)%diags (row_length, rows, 1))
    emissions(ic10h16_ibvoc)%values(:,:,:) = 0.0
    emissions(ic10h16_ibvoc)%diags (:,:,:) = 0.0
  END IF

  IF (ich3oh_ibvoc > 0) THEN
    emissions(ich3oh_ibvoc)%var_name    = 'interactive_CH3OH'

    ! Note that methanol can be referred to with two different names
    ! depending on the chemistry scheme.
    IF (ANY (em_chem_spec == 'MeOH      ')) THEN
      emissions(ich3oh_ibvoc)%tracer_name = 'MeOH'
    ELSE IF (ANY (em_chem_spec == 'CH3OH     ')) THEN
      emissions(ich3oh_ibvoc)%tracer_name = 'CH3OH'
    ELSE
      ierror   = 1
      cmessage = 'iBVOC field for methanol found in emission structure ' //    &
                 ' but not found in em_chem_spec'
      CALL ereport ('UKCA_EMISS_INIT', ierror, cmessage)
    END IF

    emissions(ich3oh_ibvoc)%l_online   = .TRUE.

    ! iBVOC emissions from JULES are in 'kg(C) m-2 s-1'. Don't indicate this
    ! in %units but use substring 'expressed as carbon' in %lng_name.
    emissions(ich3oh_ibvoc)%units      = 'kg m-2 s-1'
    emissions(ich3oh_ibvoc)%lng_name   =                                       &
     'tendency of atmosphere mass content'                               //    &
     ' of methanol expressed as carbon due to emission from vegetation'

    emissions(ich3oh_ibvoc)%hourly_fact = 'none'      ! No need to apply
    emissions(ich3oh_ibvoc)%daily_fact  = 'none'      ! time profiles

    emissions(ich3oh_ibvoc)%vert_fact   = 'surface'   ! surface emissions
    emissions(ich3oh_ibvoc)%three_dim   = .FALSE.
    emissions(ich3oh_ibvoc)%lowest_lev  = 1
    emissions(ich3oh_ibvoc)%highest_lev = 1

    ALLOCATE (emissions(ich3oh_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ich3oh_ibvoc)%diags (row_length, rows, 1))
    emissions(ich3oh_ibvoc)%values(:,:,:) = 0.0
    emissions(ich3oh_ibvoc)%diags (:,:,:) = 0.0
  END IF

  IF (ich3coch3_ibvoc > 0) THEN
    emissions(ich3coch3_ibvoc)%var_name    = 'interactive_Me2CO'
    emissions(ich3coch3_ibvoc)%tracer_name = 'Me2CO'
    emissions(ich3coch3_ibvoc)%l_online    = .TRUE.

    ! iBVOC emissions from JULES are in 'kg(C) m-2 s-1'. Don't indicate this
    ! in %units but use substring 'expressed as carbon' in %lng_name.
    emissions(ich3coch3_ibvoc)%units      = 'kg m-2 s-1'
    emissions(ich3coch3_ibvoc)%lng_name   =                                    &
     'tendency of atmosphere mass content'                               //    &
     ' of acetone expressed as carbon due to emission from vegetation'

    emissions(ich3coch3_ibvoc)%hourly_fact = 'none'      ! No need to apply
    emissions(ich3coch3_ibvoc)%daily_fact  = 'none'      ! time profiles

    emissions(ich3coch3_ibvoc)%vert_fact   = 'surface'   ! surface emissions
    emissions(ich3coch3_ibvoc)%three_dim   = .FALSE.
    emissions(ich3coch3_ibvoc)%lowest_lev  = 1
    emissions(ich3coch3_ibvoc)%highest_lev = 1

    ALLOCATE (emissions(ich3coch3_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ich3coch3_ibvoc)%diags (row_length, rows, 1))
    emissions(ich3coch3_ibvoc)%values(:,:,:) = 0.0
    emissions(ich3coch3_ibvoc)%diags (:,:,:) = 0.0
  END IF
END IF

! Initialise fields for INFERNO fire emissions in emissions structure
IF (ukca_config%l_ukca_inferno) THEN

  IF (ico_inferno > 0) THEN
    emissions(ico_inferno)%var_name    = 'inferno_CO'
    emissions(ico_inferno)%tracer_name = 'CO        '
    emissions(ico_inferno)%l_online    = .TRUE.
    emissions(ico_inferno)%units       = 'kg m-2 s-1'
    emissions(ico_inferno)%lng_name    =                                       &
     'tendency of atmosphere mass content'                               //    &
     ' of carbon monoxide due to emission from fires'

    emissions(ico_inferno)%hourly_fact = 'none'      ! No need to apply
    emissions(ico_inferno)%daily_fact  = 'none'      ! time profiles

    emissions(ico_inferno)%vert_fact   = 'surface'   ! surface emissions
    emissions(ico_inferno)%three_dim   = .FALSE.
    emissions(ico_inferno)%lowest_lev  = 1
    emissions(ico_inferno)%highest_lev = 1

    ALLOCATE (emissions(ico_inferno)%values(row_length, rows, 1))
    ALLOCATE (emissions(ico_inferno)%diags (row_length, rows, 1))
    emissions(ico_inferno)%values(:,:,:) = 0.0
    emissions(ico_inferno)%diags (:,:,:) = 0.0
  END IF

  IF (ukca_config%l_ukca_inferno_ch4 .AND.                                     &
      (.NOT. ukca_config%l_ukca_prescribech4)) THEN
    IF (ich4_inferno > 0) THEN
      emissions(ich4_inferno)%var_name    = 'inferno_CH4'
      emissions(ich4_inferno)%tracer_name = 'CH4       '
      emissions(ich4_inferno)%l_online    = .TRUE.
      emissions(ich4_inferno)%units       = 'kg m-2 s-1'
      emissions(ich4_inferno)%lng_name    =                                    &
       'tendency of atmosphere mass content'                                // &
       ' of methane due to emission from fires'

      emissions(ich4_inferno)%hourly_fact = 'none'      ! No need to apply
      emissions(ich4_inferno)%daily_fact  = 'none'      ! time profiles

      emissions(ich4_inferno)%vert_fact   = 'surface'   ! surface emissions
      emissions(ich4_inferno)%three_dim   = .FALSE.
      emissions(ich4_inferno)%lowest_lev  = 1
      emissions(ich4_inferno)%highest_lev = 1

      ALLOCATE (emissions(ich4_inferno)%values(row_length, rows, 1))
      ALLOCATE (emissions(ich4_inferno)%diags (row_length, rows, 1))
      emissions(ich4_inferno)%values(:,:,:) = 0.0
      emissions(ich4_inferno)%diags (:,:,:) = 0.0
    END IF
  END IF  !(l_ukca_inferno_ch4) .AND. (.NOT. l_ukca_prescribech4)

  IF (inox_inferno > 0) THEN
    emissions(inox_inferno)%var_name    = 'inferno_NO'
    emissions(inox_inferno)%tracer_name = 'NO       '
    emissions(inox_inferno)%l_online    = .TRUE.
    emissions(inox_inferno)%units       = 'kg m-2 s-1'
    emissions(inox_inferno)%lng_name    =                                      &
     'tendency of atmosphere mass content'                               //    &
     ' of nitrogen oxide due to emission from fires'

    emissions(inox_inferno)%hourly_fact = 'none'      ! No need to apply
    emissions(inox_inferno)%daily_fact  = 'none'      ! time profiles

    emissions(inox_inferno)%vert_fact   = 'surface'   ! surface emissions
    emissions(inox_inferno)%three_dim   = .FALSE.
    emissions(inox_inferno)%lowest_lev  = 1
    emissions(inox_inferno)%highest_lev = 1

    ALLOCATE (emissions(inox_inferno)%values(row_length, rows, 1))
    ALLOCATE (emissions(inox_inferno)%diags (row_length, rows, 1))
    emissions(inox_inferno)%values(:,:,:) = 0.0
    emissions(inox_inferno)%diags (:,:,:) = 0.0
  END IF

  IF (iso2_inferno > 0) THEN
    ecount = iso2_inferno - 1
    DO imode = 1, nmodes
      DO imoment = moment_number, moment_mass, moment_step
        IF (mode(imode) .AND. component(imode,cp_su)) THEN
          ecount = ecount + 1
          emissions(ecount)%var_name    = 'inferno_SO2_nat'
          emissions(ecount)%tracer_name = 'SO2_nat   '
          emissions(ecount)%from_emiss  = iso2_inferno

          emissions(ecount)%l_mode     = .TRUE.
          emissions(ecount)%l_mode_so2 = .TRUE.
          emissions(ecount)%l_online   = .TRUE.
          emissions(ecount)%mode       = imode
          emissions(ecount)%moment     = imoment
          emissions(ecount)%component  = cp_su
          emissions(ecount)%units      = 'kg m-2 s-1'

          emissions(ecount)%lng_name   =                                       &
           'tendency of atmosphere mass content'                         //    &
           ' of sulfur dioxide due to emission from fires'

          emissions(ecount)%hourly_fact = 'none'      ! No need to apply
          emissions(ecount)%daily_fact  = 'none'      ! time profiles

          ! SO2_Nat is a 3D emission
          emissions(ecount)%vert_fact   = 'all_levels'   ! 3D emissions
          emissions(ecount)%three_dim   = .TRUE.
          emissions(ecount)%lowest_lev  = 1
          emissions(ecount)%highest_lev = model_levels

          ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
          ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
          emissions(ecount)%values(:,:,:) = 0.0
          emissions(ecount)%diags (:,:,:) = 0.0
        END IF  ! mode and component in use
      END DO  ! imoment
    END DO  ! imode
  END IF  ! (iso2_inferno > 0)

  IF (ioc_inferno > 0) THEN
    ecount = ioc_inferno - 1
    DO imode = 1, nmodes
      DO imoment = moment_number, moment_mass, moment_step
        IF (mode(imode) .AND. component(imode,cp_oc)) THEN
          ecount = ecount + 1
          emissions(ecount)%var_name    = 'inferno_OM_biomass'
          emissions(ecount)%tracer_name = 'OM_biomass'
          emissions(ecount)%from_emiss  = ioc_inferno

          emissions(ecount)%l_mode    = .TRUE.
          emissions(ecount)%l_online  = .TRUE.
          emissions(ecount)%mode      = imode
          emissions(ecount)%moment    = imoment
          emissions(ecount)%component = cp_oc
          emissions(ecount)%units      = 'kg m-2 s-1'

          emissions(ecount)%lng_name   =                                       &
           'tendency of atmosphere mass content'                          //   &
           ' of organic matter due to emission from fires'

          emissions(ecount)%hourly_fact = 'none'      ! No need to apply
          emissions(ecount)%daily_fact  = 'none'      ! time profiles

          emissions(ecount)%vert_fact   = 'all_levels'   ! 3D emissions
          emissions(ecount)%three_dim   = .TRUE.
          emissions(ecount)%lowest_lev  = 1
          emissions(ecount)%highest_lev = model_levels

          ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
          ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))

          emissions(ecount)%values(:,:,:) = 0.0
          emissions(ecount)%diags (:,:,:) = 0.0
        END IF  ! mode and component in use
      END DO  ! imoment
    END DO  ! imode
  END IF  ! (ioc_inferno > 0)

  IF (ibc_inferno > 0) THEN
    ecount = ibc_inferno - 1
    DO imode = 1, nmodes
      DO imoment = moment_number, moment_mass, moment_step
        IF (mode(imode) .AND. component(imode,cp_bc)) THEN
          ecount = ecount + 1
          emissions(ecount)%var_name    = 'inferno_BC_biomass'
          emissions(ecount)%tracer_name = 'BC_biomass'
          emissions(ecount)%from_emiss  = ibc_inferno

          emissions(ecount)%l_mode    = .TRUE.
          emissions(ecount)%l_online  = .TRUE.
          emissions(ecount)%mode      = imode
          emissions(ecount)%moment    = imoment
          emissions(ecount)%component = cp_bc
          emissions(ecount)%units      = 'kg m-2 s-1'

          emissions(ecount)%lng_name   =                                       &
           'tendency of atmosphere mass content'                          //   &
           ' of black carbon due to emission from fires'

          emissions(ecount)%hourly_fact = 'none'      ! No need to apply
          emissions(ecount)%daily_fact  = 'none'      ! time profiles

          emissions(ecount)%vert_fact   = 'all_levels'   ! 3D emissions
          emissions(ecount)%three_dim   = .TRUE.
          emissions(ecount)%lowest_lev  = 1
          emissions(ecount)%highest_lev = model_levels

          ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
          ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))

          emissions(ecount)%values(:,:,:) = 0.0
          emissions(ecount)%diags (:,:,:) = 0.0
        END IF  ! mode and component in use
      END DO  ! imoment
    END DO  ! imode
  END IF  ! (ibc_inferno > 0)

END IF ! (l_ukca_inferno)

! Initialise online marine DMS emissions
IF (ukca_config%l_seawater_dms) THEN
  emissions(idms_seaflux)%var_name    = 'DMS_marine_emissions'
  emissions(idms_seaflux)%tracer_name = 'DMS       '
  emissions(idms_seaflux)%l_online  = .TRUE.

  !  Units should be consistent with DMS emiss from DMS_FLUX_4A,
  !  given in 'kg(S) m-2 s-1'. Similarly as seen above for 'NO_lightng',
  !  we do not indicate 'S' in the units but use the substring
  !  'expressed as sulfur' in the long_name.
  emissions(idms_seaflux)%units      = 'kg m-2 s-1'
  emissions(idms_seaflux)%lng_name   = 'tendency of atmosphere mass ' //       &
     'content of dimethyl sulfide expressed as sulfur due to emission from sea'

  emissions(idms_seaflux)%hourly_fact = 'none'      ! No need to apply
  emissions(idms_seaflux)%daily_fact  = 'none'      ! time profiles

  emissions(idms_seaflux)%vert_fact   = 'surface'   ! surface emissions
  emissions(idms_seaflux)%three_dim   = .FALSE.
  emissions(idms_seaflux)%lowest_lev  = 1
  emissions(idms_seaflux)%highest_lev = 1

  ALLOCATE (emissions(idms_seaflux)%values(row_length, rows, 1))
  ALLOCATE (emissions(idms_seaflux)%diags (row_length, rows, 1))
  emissions(idms_seaflux)%values(:,:,:) = 0.0
  emissions(idms_seaflux)%diags (:,:,:) = 0.0
END IF

! Initialise fields for seasalt emissions.
! There is a 2D mass and number field for every mode, consistent with the
! shape of the arrays returned by ukca_prim_ss.
! Number flux is multiplied by MM_DA/AVC to give units of kg(air) m-2 s-1
! Set from_emiss=iseasalt_first to identify as seasalt for ukca_emiss_mode_map
IF (glomap_config%l_ukca_primss) THEN
  ecount = iseasalt_first - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_cl)) THEN
        ecount = ecount + 1
        emissions(ecount)%var_name    = 'seasalt_online_emissions'
        emissions(ecount)%tracer_name = 'mode_emiss'
        emissions(ecount)%from_emiss  = iseasalt_first

        emissions(ecount)%l_mode    = .TRUE.
        emissions(ecount)%l_online  = .TRUE.
        emissions(ecount)%mode      = imode
        emissions(ecount)%moment    = imoment
        emissions(ecount)%component = cp_cl

        IF (imoment == moment_number) THEN
          emissions(ecount)%units      = 'kg(air) m-2 s-1'
        ELSE
          emissions(ecount)%units      = 'kg m-2 s-1'
        END IF
        emissions(ecount)%std_name  = 'tendency_of_atmosphere_mass_' //        &
                     'content_of_seasalt_dry_aerosol_due_to_emission'

        emissions(ecount)%hourly_fact = 'none'      ! No need to apply
        emissions(ecount)%daily_fact  = 'none'      ! time profiles

        emissions(ecount)%vert_fact   = 'surface'   ! surface emissions
        emissions(ecount)%three_dim   = .FALSE.
        emissions(ecount)%lowest_lev  = 1
        emissions(ecount)%highest_lev = 1

        ALLOCATE (emissions(ecount)%values(row_length, rows, 1))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, 1))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF  ! l_ukca_primss

! Initialise fields for primary marine organic carbon emissions (PMOC).
! Allow emission into any mode for flexibility. Emissions are 2D and emit into
! OC component.
! Number flux is multiplied by MM_DA/AVC to give units of kg(air) m-2 s-1
! Set from_emiss=ipmao_first to identify as pmao for ukca_emiss_mode_map
IF (glomap_config%l_ukca_prim_moc) THEN
  ecount = ipmoc_first - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_oc)) THEN
        ecount = ecount + 1
        emissions(ecount)%var_name    = 'pmoc_online_emissions'
        emissions(ecount)%tracer_name = 'mode_emiss'
        emissions(ecount)%from_emiss  = ipmoc_first

        emissions(ecount)%l_mode    = .TRUE.
        emissions(ecount)%l_online  = .TRUE.
        emissions(ecount)%mode      = imode
        emissions(ecount)%moment    = imoment
        emissions(ecount)%component = cp_oc

        IF (imoment == moment_number) THEN
          emissions(ecount)%units      = 'kg(air) m-2 s-1'
        ELSE
          emissions(ecount)%units      = 'kg m-2 s-1'
        END IF
        emissions(ecount)%std_name  = 'tendency_of_atmosphere_mass_' //        &
           'content_of_particulate_organic_matter_dry_aerosol_due_to_emission'

        emissions(ecount)%hourly_fact = 'none'      ! No need to apply
        emissions(ecount)%daily_fact  = 'none'      ! time profiles

        emissions(ecount)%vert_fact   = 'surface'   ! surface emissions
        emissions(ecount)%three_dim   = .FALSE.
        emissions(ecount)%lowest_lev  = 1
        emissions(ecount)%highest_lev = 1

        ALLOCATE (emissions(ecount)%values(row_length, rows, 1))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, 1))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF  ! l_ukca_prim_moc

! Initialise fields for dust emissions.
! There is a 2D mass and number field for every mode, consistent with the
! shape of the arrays returned by ukca_prim_du.
! Number flux is multiplied by MM_DA/AVC to give units of kg(air) m-2 s-1
! Set from_emiss=idust_first to identify as dust for ukca_emiss_mode_map
IF (glomap_config%l_ukca_primdu) THEN
  ecount = idust_first - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_du)) THEN
        ecount = ecount + 1
        emissions(ecount)%var_name    = 'dust_online_emissions'
        emissions(ecount)%tracer_name = 'mode_emiss'
        emissions(ecount)%from_emiss  = idust_first

        emissions(ecount)%l_mode    = .TRUE.
        emissions(ecount)%l_online  = .TRUE.
        emissions(ecount)%mode      = imode
        emissions(ecount)%moment    = imoment
        emissions(ecount)%component = cp_du

        IF (imoment == moment_number) THEN
          emissions(ecount)%units      = 'kg(air) m-2 s-1'
        ELSE
          emissions(ecount)%units      = 'kg m-2 s-1'
        END IF
        emissions(ecount)%std_name  = 'tendency_of_atmosphere_mass_' //        &
                     'content_of_dust_dry_aerosol_due_to_emission'

        emissions(ecount)%hourly_fact = 'none'      ! No need to apply
        emissions(ecount)%daily_fact  = 'none'      ! time profiles

        emissions(ecount)%vert_fact   = 'surface'   ! surface emissions
        emissions(ecount)%three_dim   = .FALSE.
        emissions(ecount)%lowest_lev  = 1
        emissions(ecount)%highest_lev = 1

        ALLOCATE (emissions(ecount)%values(row_length, rows, 1))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, 1))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF  ! l_ukca_primdu


! Initialise fields for nitrate and ammonium emissions.
! There is a 2D mass and number field for every mode, consistent with the
! shape of the arrays returned by ukca_prim_du.
! Number flux is multiplied by MM_DA/AVC to give units of kg(air) m-2 s-1
! Set from_emiss=ino3_first and inh4_first to identify as
! no3 or nh4 for ukca_emiss_mode_map
IF (glomap_config%l_ukca_fine_no3_prod) THEN
  ecount = ino3_first - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_no3)) THEN
        ecount = ecount + 1
        emissions(ecount)%var_name    = 'no3_online_emissions'
        emissions(ecount)%tracer_name = 'mode_emiss'
        emissions(ecount)%from_emiss  = ino3_first

        emissions(ecount)%l_mode     = .TRUE.
        emissions(ecount)%l_online   = .TRUE.
        emissions(ecount)%mode       = imode
        emissions(ecount)%moment     = imoment
        emissions(ecount)%component  = cp_no3

        IF (imoment == moment_number) THEN
          emissions(ecount)%units      = 'kg(air) m-2 s-1'
        ELSE
          emissions(ecount)%units      = 'kg m-2 s-1'
        END IF

        emissions(ecount)%lng_name   =                                         &
         'tendency of atmosphere mass content'                         //      &
         ' of nitrate aerosol from atmos chemistry'

        emissions(ecount)%hourly_fact = 'none'      ! No need to apply
        emissions(ecount)%daily_fact  = 'none'      ! time profiles

        ! SO2_Nat is a 3D emission
        emissions(ecount)%vert_fact   = 'all_levels'   ! 3D emissions
        emissions(ecount)%three_dim   = .TRUE.
        emissions(ecount)%lowest_lev  = 1
        emissions(ecount)%highest_lev = model_levels

        ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
  ecount = inh4_first - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_nh4)) THEN
        ecount = ecount + 1
        emissions(ecount)%var_name    = 'nh4_online_emissions'
        emissions(ecount)%tracer_name = 'mode_emiss'
        emissions(ecount)%from_emiss  = inh4_first

        emissions(ecount)%l_mode     = .TRUE.
        emissions(ecount)%l_online   = .TRUE.
        emissions(ecount)%mode       = imode
        emissions(ecount)%moment     = imoment
        emissions(ecount)%component  = cp_nh4

        IF (imoment == moment_number) THEN
          emissions(ecount)%units      = 'kg(air) m-2 s-1'
        ELSE
          emissions(ecount)%units      = 'kg m-2 s-1'
        END IF

        emissions(ecount)%lng_name   =                                         &
         'tendency of atmosphere mass content'                         //      &
         ' of ammonium aerosol from atmos chemistry'

        emissions(ecount)%hourly_fact = 'none'      ! No need to apply
        emissions(ecount)%daily_fact  = 'none'      ! time profiles

        ! SO2_Nat is a 3D emission
        emissions(ecount)%vert_fact   = 'all_levels'   ! 3D emissions
        emissions(ecount)%three_dim   = .TRUE.
        emissions(ecount)%lowest_lev  = 1
        emissions(ecount)%highest_lev = model_levels

        ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF

IF (glomap_config%l_ukca_coarse_no3_prod) THEN
  ecount = iseasalt_hno3 - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_cl)) THEN
        ecount = ecount + 1
        emissions(ecount)%var_name    = 'seasalt_cnvtd_to_nitrate'
        emissions(ecount)%tracer_name = 'mode_emiss'
        emissions(ecount)%from_emiss  = iseasalt_hno3

        emissions(ecount)%l_mode     = .TRUE.
        emissions(ecount)%l_online   = .TRUE.
        emissions(ecount)%mode       = imode
        emissions(ecount)%moment     = imoment
        emissions(ecount)%component  = cp_cl

        IF (imoment == moment_number) THEN
          emissions(ecount)%units      = 'kg(air) m-2 s-1'
        ELSE
          emissions(ecount)%units      = 'kg m-2 s-1'
        END IF

        emissions(ecount)%lng_name   =                                         &
         'tendency of atmosphere mass content'                         //      &
         ' of seasalt from nitrate chemistry'

        emissions(ecount)%hourly_fact = 'none'      ! No need to apply
        emissions(ecount)%daily_fact  = 'none'      ! time profiles

        ! SO2_Nat is a 3D emission
        emissions(ecount)%vert_fact   = 'all_levels'   ! 3D emissions
        emissions(ecount)%three_dim   = .TRUE.
        emissions(ecount)%lowest_lev  = 1
        emissions(ecount)%highest_lev = model_levels

        ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
  ecount = inano3_cond - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_nn)) THEN
        ecount = ecount + 1
        emissions(ecount)%var_name    = 'nano3_condensed_on_coarse_aero'
        emissions(ecount)%tracer_name = 'mode_emiss'
        emissions(ecount)%from_emiss  = inano3_cond

        emissions(ecount)%l_mode     = .TRUE.
        emissions(ecount)%l_online   = .TRUE.
        emissions(ecount)%mode       = imode
        emissions(ecount)%moment     = imoment
        emissions(ecount)%component  = cp_nn

        IF (imoment == moment_number) THEN
          emissions(ecount)%units      = 'kg(air) m-2 s-1'
        ELSE
          emissions(ecount)%units      = 'kg m-2 s-1'
        END IF

        emissions(ecount)%lng_name   =                                         &
         'tendency of atmosphere mass content'                         //      &
         ' of nano3 aerosol from condensation'

        emissions(ecount)%hourly_fact = 'none'      ! No need to apply
        emissions(ecount)%daily_fact  = 'none'      ! time profiles

        ! SO2_Nat is a 3D emission
        emissions(ecount)%vert_fact   = 'all_levels'   ! 3D emissions
        emissions(ecount)%three_dim   = .TRUE.
        emissions(ecount)%lowest_lev  = 1
        emissions(ecount)%highest_lev = model_levels

        ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF

IF ((glomap_config%l_ukca_coarse_no3_prod) .AND.                               &
    (glomap_config%l_ukca_primdu)) THEN
  ecount = idust_hno3 - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_du)) THEN
        ecount = ecount + 1
        emissions(ecount)%var_name    = 'dust_cnvtd_to_nitrate'
        emissions(ecount)%tracer_name = 'mode_emiss'
        emissions(ecount)%from_emiss  = idust_hno3

        emissions(ecount)%l_mode     = .TRUE.
        emissions(ecount)%l_online   = .TRUE.
        emissions(ecount)%mode       = imode
        emissions(ecount)%moment     = imoment
        emissions(ecount)%component  = cp_du

        IF (imoment == moment_number) THEN
          emissions(ecount)%units      = 'kg(air) m-2 s-1'
        ELSE
          emissions(ecount)%units      = 'kg m-2 s-1'
        END IF

        emissions(ecount)%lng_name   =                                         &
         'tendency of atmosphere mass content'                         //      &
         ' of dust from nitrate chemistry'

        emissions(ecount)%hourly_fact = 'none'      ! No need to apply
        emissions(ecount)%daily_fact  = 'none'      ! time profiles

        ! SO2_Nat is a 3D emission
        emissions(ecount)%vert_fact   = 'all_levels'   ! 3D emissions
        emissions(ecount)%three_dim   = .TRUE.
        emissions(ecount)%lowest_lev  = 1
        emissions(ecount)%highest_lev = model_levels

        ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF

emissions(num_cdf_em_flds+1:num_em_flds)%l_update = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN

END SUBROUTINE ukca_onl_emiss_init

! ---------------------------------------------------------------------
! Description:
!   Initialise the 3D/2D arrays in the emissions structure.
!
! ---------------------------------------------------------------------
SUBROUTINE ukca_emiss_spatial_vars_init (row_length, rows,  model_levels)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT(IN) :: row_length       ! model dimensions
INTEGER,           INTENT(IN) :: rows
INTEGER,           INTENT(IN) :: model_levels

! Local variables
INTEGER :: ecount          ! Index in emissions array

! GLOMAP related variables
INTEGER :: moment_step = moment_mass - moment_number
INTEGER :: imode
INTEGER :: imoment

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_SPATIAL_VARS_INIT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialise fields for lightning NOx in the emissions structure
! (Note that there is no need to update file_name, update_freq
! and update_type for online emissions)
IF (ukca_config%i_ukca_light_param /= i_light_param_off) THEN
  ALLOCATE (emissions(inox_light)%values(row_length, rows, model_levels))
  ALLOCATE (emissions(inox_light)%diags (row_length, rows, model_levels))
  emissions(inox_light)%values(:,:,:) = 0.0
  emissions(inox_light)%diags (:,:,:) = 0.0
END IF

! Initialise fields for online wetland emiss of CH4 in emissions structure
! (done in a similar way as done above for lightning NOx, but only if
! these emissions exist)
IF (ukca_config%l_ukca_qch4inter .AND.                                         &
    (.NOT. ukca_config%l_ukca_prescribech4)) THEN
  ALLOCATE (emissions(ich4_wetl)%values(row_length, rows, 1))
  ALLOCATE (emissions(ich4_wetl)%diags (row_length, rows, 1))
  emissions(ich4_wetl)%values(:,:,:) = 0.0
  emissions(ich4_wetl)%diags (:,:,:) = 0.0
END IF

! Initialise fields for iBVOC in emissions structure
IF (ukca_config%l_ukca_ibvoc) THEN
  IF (ic5h8_ibvoc > 0) THEN
    ALLOCATE (emissions(ic5h8_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ic5h8_ibvoc)%diags (row_length, rows, 1))
    emissions(ic5h8_ibvoc)%values(:,:,:) = 0.0
    emissions(ic5h8_ibvoc)%diags (:,:,:) = 0.0
  END IF

  IF (ic10h16_ibvoc > 0) THEN
    ALLOCATE (emissions(ic10h16_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ic10h16_ibvoc)%diags (row_length, rows, 1))
    emissions(ic10h16_ibvoc)%values(:,:,:) = 0.0
    emissions(ic10h16_ibvoc)%diags (:,:,:) = 0.0
  END IF

  IF (ich3oh_ibvoc > 0) THEN
    ALLOCATE (emissions(ich3oh_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ich3oh_ibvoc)%diags (row_length, rows, 1))
    emissions(ich3oh_ibvoc)%values(:,:,:) = 0.0
    emissions(ich3oh_ibvoc)%diags (:,:,:) = 0.0
  END IF

  IF (ich3coch3_ibvoc > 0) THEN
    ALLOCATE (emissions(ich3coch3_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ich3coch3_ibvoc)%diags (row_length, rows, 1))
    emissions(ich3coch3_ibvoc)%values(:,:,:) = 0.0
    emissions(ich3coch3_ibvoc)%diags (:,:,:) = 0.0
  END IF
END IF

! Initialise fields for INFERNO fire emissions in emissions structure
IF (ukca_config%l_ukca_inferno) THEN

  IF (ico_inferno > 0) THEN
    ALLOCATE (emissions(ico_inferno)%values(row_length, rows, 1))
    ALLOCATE (emissions(ico_inferno)%diags (row_length, rows, 1))
    emissions(ico_inferno)%values(:,:,:) = 0.0
    emissions(ico_inferno)%diags (:,:,:) = 0.0
  END IF

  IF (ukca_config%l_ukca_inferno_ch4 .AND.                                     &
      (.NOT. ukca_config%l_ukca_prescribech4)) THEN
    IF (ich4_inferno > 0) THEN
      ALLOCATE (emissions(ich4_inferno)%values(row_length, rows, 1))
      ALLOCATE (emissions(ich4_inferno)%diags (row_length, rows, 1))
      emissions(ich4_inferno)%values(:,:,:) = 0.0
      emissions(ich4_inferno)%diags (:,:,:) = 0.0
    END IF
  END IF  !(l_ukca_inferno_ch4) .AND. (.NOT. l_ukca_prescribech4)

  IF (inox_inferno > 0) THEN
    ALLOCATE (emissions(inox_inferno)%values(row_length, rows, 1))
    ALLOCATE (emissions(inox_inferno)%diags (row_length, rows, 1))
    emissions(inox_inferno)%values(:,:,:) = 0.0
    emissions(inox_inferno)%diags (:,:,:) = 0.0
  END IF

  IF (iso2_inferno > 0) THEN
    ecount = iso2_inferno - 1
    DO imode = 1, nmodes
      DO imoment = moment_number, moment_mass, moment_step
        IF (mode(imode) .AND. component(imode,cp_su)) THEN
          ecount = ecount + 1
          ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
          ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
          emissions(ecount)%values(:,:,:) = 0.0
          emissions(ecount)%diags (:,:,:) = 0.0
        END IF  ! mode and component in use
      END DO  ! imoment
    END DO  ! imode
  END IF  ! (iso2_inferno > 0)

  IF (ioc_inferno > 0) THEN
    ecount = ioc_inferno - 1
    DO imode = 1, nmodes
      DO imoment = moment_number, moment_mass, moment_step
        IF (mode(imode) .AND. component(imode,cp_oc)) THEN
          ecount = ecount + 1
          ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
          ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
          emissions(ecount)%values(:,:,:) = 0.0
          emissions(ecount)%diags (:,:,:) = 0.0
        END IF  ! mode and component in use
      END DO  ! imoment
    END DO  ! imode
  END IF  ! (ioc_inferno > 0)

  IF (ibc_inferno > 0) THEN
    ecount = ibc_inferno - 1
    DO imode = 1, nmodes
      DO imoment = moment_number, moment_mass, moment_step
        IF (mode(imode) .AND. component(imode,cp_bc)) THEN
          ecount = ecount + 1
          ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
          ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
          emissions(ecount)%values(:,:,:) = 0.0
          emissions(ecount)%diags (:,:,:) = 0.0
        END IF  ! mode and component in use
      END DO  ! imoment
    END DO  ! imode
  END IF  ! (ibc_inferno > 0)

END IF ! (l_ukca_inferno)

! Initialise online marine DMS emissions
IF (ukca_config%l_seawater_dms) THEN
  ALLOCATE (emissions(idms_seaflux)%values(row_length, rows, 1))
  ALLOCATE (emissions(idms_seaflux)%diags (row_length, rows, 1))
  emissions(idms_seaflux)%values(:,:,:) = 0.0
  emissions(idms_seaflux)%diags (:,:,:) = 0.0
END IF

! Initialise fields for seasalt emissions.
! There is a 2D mass and number field for every mode, consistent with the
! shape of the arrays returned by ukca_prim_ss.
! Number flux is multiplied by MM_DA/AVC to give units of kg(air) m-2 s-1
! Set from_emiss=iseasalt_first to identify as seasalt for ukca_emiss_mode_map
IF (glomap_config%l_ukca_primss) THEN
  ecount = iseasalt_first - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_cl)) THEN
        ecount = ecount + 1
        ALLOCATE (emissions(ecount)%values(row_length, rows, 1))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, 1))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF  ! l_ukca_primss

! Initialise fields for primary marine organic carbon emissions (PMOC).
! Allow emission into any mode for flexibility. Emissions are 2D and emit into
! OC component.
! Number flux is multiplied by MM_DA/AVC to give units of kg(air) m-2 s-1
! Set from_emiss=ipmao_first to identify as pmao for ukca_emiss_mode_map
IF (glomap_config%l_ukca_prim_moc) THEN
  ecount = ipmoc_first - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_oc)) THEN
        ecount = ecount + 1
        ALLOCATE (emissions(ecount)%values(row_length, rows, 1))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, 1))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF  ! l_ukca_prim_moc

! Initialise fields for dust emissions.
! There is a 2D mass and number field for every mode, consistent with the
! shape of the arrays returned by ukca_prim_du.
! Number flux is multiplied by MM_DA/AVC to give units of kg(air) m-2 s-1
! Set from_emiss=idust_first to identify as dust for ukca_emiss_mode_map
IF (glomap_config%l_ukca_primdu) THEN
  ecount = idust_first - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_du)) THEN
        ecount = ecount + 1
        ALLOCATE (emissions(ecount)%values(row_length, rows, 1))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, 1))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF  ! l_ukca_primdu


! Initialise fields for nitrate and ammonium emissions.
! There is a 2D mass and number field for every mode, consistent with the
! shape of the arrays returned by ukca_prim_du.
! Number flux is multiplied by MM_DA/AVC to give units of kg(air) m-2 s-1
! Set from_emiss=ino3_first and inh4_first to identify as
! no3 or nh4 for ukca_emiss_mode_map
IF (glomap_config%l_ukca_fine_no3_prod) THEN
  ecount = ino3_first - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_no3)) THEN
        ecount = ecount + 1
        ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
  ecount = inh4_first - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_nh4)) THEN
        ecount = ecount + 1
        ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF

IF (glomap_config%l_ukca_coarse_no3_prod) THEN
  ecount = iseasalt_hno3 - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_cl)) THEN
        ecount = ecount + 1
        ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
  ecount = inano3_cond - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_nn)) THEN
        ecount = ecount + 1
        ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF

IF ((glomap_config%l_ukca_coarse_no3_prod) .AND.                               &
    (glomap_config%l_ukca_primdu)) THEN
  ecount = idust_hno3 - 1
  DO imode = 1, nmodes
    DO imoment = moment_number, moment_mass, moment_step
      IF (mode(imode) .AND. component(imode,cp_du)) THEN
        ecount = ecount + 1
        ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF

emissions(num_cdf_em_flds+1:num_em_flds)%l_update = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN

END SUBROUTINE ukca_emiss_spatial_vars_init

!---------------------------------------------------------------------------
! Description:
!   Convert total aerosol emission mass into mass and number for each mode.
!
! Method:
!   For each mode and component this emission field maps to, compute number and
!   mass, and add to appropriate emissions structures.
!   This routine doesn't care whether the emissions are 2d or 3d: it just
!   converts the source mass into number and mass for each mode. The (level and
!   time) attributes of the modal emissions are defined in ukca_emiss_init.
!   The source emission is normally read from a netcdf file (offline emission),
!   but could be calculated online by the model.
!---------------------------------------------------------------------------
SUBROUTINE ukca_emiss_update_mode(row_length, rows, model_levels,              &
                                  tracer_name, three_dim, emiss_values,        &
                                  emmas, emnum, lmode_emiss,                   &
                                  nlev, icp)

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: row_length    ! Model dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

CHARACTER (LEN=10), INTENT(IN) :: tracer_name  ! name of emission tracer
LOGICAL, INTENT(IN) :: three_dim  ! source emissions are 3D

! Source emission mass (kg(component)/m2/s) which needs converting to modal
! emissions. Assumed-shape array because 3rd dim can be 1 or model_levels.
REAL, INTENT(IN) :: emiss_values(:,:,:)

! 3D number (equivalent kg/m2/s) and mass (kg(component)/m2/s) emissions for
! each mode.
REAL, INTENT(OUT) :: emnum(row_length, rows, model_levels, nmodes)
REAL, INTENT(OUT) :: emmas(row_length, rows, model_levels, nmodes)

! Which modes receive some emission for this species (TRUE for those that do)
LOGICAL, INTENT(OUT) :: lmode_emiss(nmodes)

INTEGER, INTENT(OUT) :: nlev  ! Number of levels emission is applied to
INTEGER, INTENT(OUT) :: icp   ! Aerosol component emission is applied to

! Mass (kg/m2/s) and particle volume (nm3/m2/s) emissions for an individual mode
REAL :: modemass(row_length, rows, model_levels)
REAL :: modevol(row_length, rows, model_levels)

INTEGER :: imode            ! Loop index for mode number.
REAL :: mode_frac(nmodes)   ! Fraction of emission in each mode (0->1).
REAL :: mode_diam(nmodes)   ! Geometric mean diameter for each mode emiss (nm)
REAL :: mode_stdev(nmodes)  ! Geometric stdev for each mode emiss (units=1)

REAL    :: mm_da  !  Molar mass of dry air (kg/mol) =avc*zboltz/ra
REAL    :: factor  !  converts from partcls per m2/s to kg/m2/s
REAL    :: lstdev  !  log(stdev)

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_UPDATE_MODE'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Call ukca_def_mode_emiss to obtain emission mask, mode fraction, diameter
! and stdev for each mode, as well as component emissions applied to (=icp).
CALL ukca_def_mode_emiss(tracer_name, lmode_emiss, icp,                        &
                         mode_frac, mode_diam, mode_stdev)

! emiss%values might only have 1 level, so ensure consistency when assiging to
! modemass. All other arrays used below have all levels, so don't need to use
! nlev after that.
IF (three_dim) THEN
  nlev = model_levels
ELSE
  nlev = 1
END IF

mm_da = avogadro*boltzmann/rgas  ! Molar mass of dry air (kg/mol)
factor = mm_da / avogadro   ! Converts from partcls per m2/s to kg/m2/s
modemass(:,:,:) = 0.0
modevol(:,:,:) = 0.0
emnum(:,:,:,:) = 0.0
emmas(:,:,:,:) = 0.0

! For each mode, calculate mass, volume and hence number emissions and save into
! emmas and emnum arrays.
DO imode = 1, nmodes
  IF (lmode_emiss(imode)) THEN
    lstdev = LOG(mode_stdev(imode))

    ! Particulate mass emissions (kg of component /m2/s) in mode.
    ! Note that units conversions (C to OC, S to SO2) are done in
    ! ukca_new_emiss_ctl using base_scaling, and scaling for particulate
    ! fraction and conversion from SO2 to H2SO4 is done in ukca_add_emiss using
    ! aero_scaling.
    modemass(:,:,1:nlev) = emiss_values(:,:,1:nlev) * mode_frac(imode)
    emmas(:,:,:,imode) = modemass(:,:,:)

    ! Total particle volume emission rate (nm3 per m2 per s) in this mode is
    ! the rate of mass emission (kg/m2/s) divided by the density of this mode
    ! (kg/m3) multiplied by 1e27 (1e9^3) to go from m3/m2/s to nm3/m2/s.
    modevol(:,:,:) = 1e27*modemass(:,:,:)/rhocomp(icp)

    ! Total particle number (particles per m2 per s) * factor.
    ! factor converts from ptcls/m2/s to equivalent kg/m2/s as dry air.
    emnum(:,:,:,imode) = modevol(:,:,:) * factor /                             &
      ( (pi/6.0) * (mode_diam(imode)**3) * EXP(4.5*lstdev*lstdev) )

  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_emiss_update_mode


!---------------------------------------------------------------------------
! Description:
!   Map number and mass emissions arrays to appropriate mode emissions
!   structures
!
! Method:
!   For a given emission field (from file or online), take the provided mass and
!   number and find the emissions structures with the corresponding moment,
!   mode, component and source species. Copy the emissions into
!   emissions(l)%values.
!---------------------------------------------------------------------------
SUBROUTINE ukca_emiss_mode_map(row_length, rows, model_levels, emmas, emnum,   &
                               emiss_levs, icp, lmode_emiss, from_emiss)

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: row_length    ! Model dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! 3D number (equivalent kg/m2/s) and mass (kg(component)/m2/s) emissions for
! each mode.
REAL, INTENT(IN) :: emmas(row_length, rows, model_levels, nmodes)
REAL, INTENT(IN) :: emnum(row_length, rows, model_levels, nmodes)

! Number levels of emission field. To ensure that emnum/emmas slice is correct
! shape for emissions(l)%values.
INTEGER, INTENT(IN) :: emiss_levs

INTEGER, INTENT(IN) :: icp         ! Aerosol component for this emission
INTEGER, INTENT(IN) :: from_emiss  ! Index of source emission
LOGICAL, INTENT(IN) :: lmode_emiss(nmodes)   ! Which modes receive emission?

! Local variables
INTEGER :: nlev            ! Short name from emiss_levs argument
INTEGER :: imode, imoment  ! Loop indices for mode and moment number
INTEGER :: l               ! Loop index for emissions structure array
LOGICAL :: lfound          ! Flag for whether matching emission found
INTEGER :: moment_step = moment_mass - moment_number

CHARACTER (LEN=errormessagelength) :: cmessage    ! Error message
INTEGER                            :: ierror      ! Error code

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_MODE_MAP'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

nlev = emiss_levs

! For each mode and for both number and mass, loop through all emissions
! structures until we find a mode emission with matching mode, component,
! moment and name of source tracer. Copy into the %values tag of the matching
! structure.
DO imode = 1, nmodes
  IF (lmode_emiss(imode)) THEN
    DO imoment = moment_number, moment_mass, moment_step
      lfound = .FALSE.
      l_loop: DO l = 1, num_em_flds
        ! Find the corresponding mode emissions structures for mass and number
        IF (emissions(l)%l_mode .AND.                                          &
            emissions(l)%mode == imode .AND.                                   &
            emissions(l)%component  == icp .AND.                               &
            emissions(l)%moment == imoment .AND.                               &
            emissions(l)%from_emiss == from_emiss) THEN
          SELECT CASE(imoment)
          CASE (moment_number)
            emissions(l)%values(:,:,1:nlev) = emnum(:,:,1:nlev,imode)
          CASE (moment_mass)
            emissions(l)%values(:,:,1:nlev) = emmas(:,:,1:nlev,imode)
          END SELECT

          lfound = .TRUE.
          EXIT l_loop
        END IF  ! mode==imode and component==icp
      END DO l_loop  ! num_em_flds
    END DO  ! imoment

    IF (.NOT. lfound) THEN
      ! If this happens something has gone wrong in ukca_emiss_init or
      ! ukca_def_mode_emiss
      WRITE(cmessage,'(2A,3I3)') 'Missing mode emission structure: ',          &
        from_emiss, imode, icp, imoment
      ierror = icp
      CALL ereport('UKCA_EMISS_MODE_MAP', ierror, cmessage)
    END IF
  END IF  ! lmode_emiss(imode) and component(imode,icp)
END DO  ! nmodes

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_emiss_mode_map

END MODULE ukca_emiss_mod
