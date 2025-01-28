! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module with procedures for interfacing between emission file data and UKCA
!
! Method:
!   Provides the list of required emission fields to the parent model, lets the
!   parent define the emission attributes for UKCA and passess to UKCA the
!   emission values on its grid and relevant for current timestep.
!    1. UKCA_GET_EMISSION_VARLIST: Returns the list of active emission species
!           and the number of expected MODE emissions.
!    2. UKCA_REGISTER_EMISSION: Setup a single emission field by populating
!           its slot in the emissions data structure on UKCA side
!    3. UKCA_SET_EMISSION: Pass emission value received from parent into the
!           corresponding slot in the emissions data structure.
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
MODULE ukca_emiss_api_mod
!
USE ukca_emiss_struct_mod, ONLY: ncdf_emissions, ukca_em_struct_init
USE ukca_config_defs_mod,  ONLY: em_chem_spec
USE ukca_emiss_mode_mod,   ONLY: aero_ems_species, ukca_def_mode_emiss
USE ukca_mode_setup,       ONLY: nmodes, moment_number, moment_mass
USE ukca_emiss_factors,    ONLY: vertical_emiss_factors, base_emiss_factors
USE ukca_config_specification_mod, ONLY: ukca_config, l_ukca_config_available
USE ukca_error_mod,        ONLY: maxlen_message, maxlen_procname,              &
                                 errcode_ukca_uninit

! UM specific modules
USE ereport_mod,           ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,              ONLY: jpim, jprb      ! DrHook
USE yomhook,               ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Default private
PRIVATE

INTEGER, PUBLIC :: n_ems_registered = 0
                            ! Number of emissions fields registered,
                            ! to be returned to the parent as emiss_id so that
                            ! each field has a unique identifier
INTEGER :: ecount = 0       ! Index of actual position in data array
                            ! after considering aerosol slots
INTEGER :: n_cdf_emiss      ! Total number of emissions from files
INTEGER, ALLOCATABLE :: emiss_map(:)
                            ! Mapping of source emission to actual slot
                            ! in the data array

! Subroutines available outside this module
PUBLIC :: ukca_get_emission_varlist, ukca_register_emission,                   &
          ukca_set_emission, get_registered_ems_info

! Dr Hook parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EMISS_API_MOD'

CONTAINS

! ---------------------------------------------------------------------
! Description:
!  Return the list of active emission species in this configuration, alongwith
!  the number of MODE emission slots for data read from files
!  Need to use pointer arguments (and targets) as the number of species
!  is not known to the calling routine, so is unable to fix the size of
!  the array arguments before calling.
! ---------------------------------------------------------------------
SUBROUTINE ukca_get_emission_varlist (emiss_species, num_per_species,          &
                                      error_code, error_message, error_routine)

USE ukca_emiss_struct_mod, ONLY: maxlen_emiss_tracer_name

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=maxlen_emiss_tracer_name), POINTER, INTENT(OUT) ::               &
  emiss_species(:)                                    ! Names of species
INTEGER, POINTER, INTENT(OUT)  :: num_per_species(:)  ! No. of emission slots
                               ! for each species. Will be > 1 for aerosol emiss

! Error code for status reporting
INTEGER, INTENT(OUT) :: error_code

! Further arguments for status reporting
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

CHARACTER(LEN=maxlen_emiss_tracer_name), ALLOCATABLE, SAVE, TARGET :: em_spec(:)
INTEGER, ALLOCATABLE, SAVE, TARGET :: n_per_species(:)
INTEGER :: n, n_species, icp
LOGICAL :: lmode_emiss(nmodes)

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_GET_EMISSION_VARLIST'

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Nullify pointer before use
emiss_species => NULL()
num_per_species => NULL()

! Set defaults for output arguments
error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check if UKCA configuration data has been initialised - needed to obtain
! the emissions information
IF (.NOT. l_ukca_config_available) THEN
  error_code = errcode_ukca_uninit
  IF (PRESENT(error_message))                                                  &
    error_message = 'No UKCA configuration has been set up'
  IF (PRESENT(error_routine)) error_routine = RoutineName

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Store active emission species and names
! (The list will be empty if emissions are turned off)

IF (ukca_config%l_ukca_emissions_off .OR. (.NOT. ukca_config%l_ukca_chem)) THEN
  n_species = 0
ELSE
  n_species = SIZE(em_chem_spec)
END IF

ALLOCATE(em_spec(n_species))
ALLOCATE(n_per_species(n_species))

IF (n_species > 0) THEN
  em_spec(:) = ''
  n_per_species(:) = 0
END IF

! Store names of emission fields, including those feeding into aerosols.
! At the same time calculate the expected slots required in the emissions
! structure for each (i.e. the no. of modes that a given aerosol species will
! get added to). The latter is obtained by calling UKCA_DEF_MODE_EMISS which
! returns this number in the lmode_emiss array. Note that in reality, two slots
! are required for each mode (to hold mass and number)
DO n = 1, n_species
  em_spec(n) = em_chem_spec(n)
  IF ( ukca_config%l_ukca_mode .AND. ANY(aero_ems_species == em_spec(n)) ) THEN
    lmode_emiss(:) = .FALSE.
    CALL ukca_def_mode_emiss(em_spec(n), lmode_emiss, icp)
    n_per_species(n) = 1 + 2*COUNT(lmode_emiss)
  ELSE
    n_per_species(n) = 1   ! One entry for each non-aerosol species
  END IF
END DO      ! loop over all species

emiss_species => em_spec
num_per_species => n_per_species

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE ukca_get_emission_varlist

! ----------------------------------------------------------------------------
! Description:
!  Sets up an entry in the UKCA emissions structure, based on component values
!  received as arguments and returns an unique emission identifier.
!  On the first call, expects the total number of slots to be able to allocate
!  the structure.
! --------------------------------------------------------------------------
SUBROUTINE ukca_register_emission (num_tot, varname, tracer_name, units,       &
                three_dim, emiss_id, l_omit_std_name, std_name, long_name,     &
                vert_fact, lowest_lev, highest_lev, hourly_fact, daily_fact)

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)          :: num_tot     ! Total number of emission slots
CHARACTER(LEN=*), INTENT(IN) :: varname     ! Variable name for referring to
                                            ! field in log/error messages
CHARACTER(LEN=*), INTENT(IN) :: tracer_name ! Name of emissions species
CHARACTER(LEN=*), INTENT(IN) :: units       ! Emission units
LOGICAL, INTENT(IN)          :: three_dim   ! T if the field is full height
INTEGER, INTENT(OUT)         :: emiss_id    ! Return: Position in fields list

! Emissions must be in kg m-2 s-1. A mandatory units argument avoids ambiguity
! and must contain an equivalent string.
! If emissions are not expressed in kg m-2 s-1 of the emitted species a
! descriptive name must be included (a CF standard name std_name or an
! alternative long_name).
! To force conversion, the descriptive name must contain the string
! 'expressed_as_<element>' (or 'expressed as <element>' in the case of
! long_name) where <element> is 'nitrogen', 'carbon' or 'sulfur'.
! If no conversion is required the need for a descriptive name can be
! optionally overriden by setting l_omit_std_name = T.
LOGICAL, OPTIONAL, INTENT(IN) :: l_omit_std_name
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: std_name
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: long_name

! Keywords relating to vertical and temporal scaling factors
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: vert_fact
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: hourly_fact
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: daily_fact
INTEGER, OPTIONAL, INTENT(IN)          :: lowest_lev
INTEGER, OPTIONAL, INTENT(IN)          :: highest_lev

! Local variables

INTEGER :: icp         ! Indices for aerosol component and mode into which a
INTEGER :: imode       ! source emission field is directly apportioned
INTEGER :: imoment     ! Indicator for mass (0) or number (3) of aerosol field

INTEGER :: from_emiss  ! Emissions index from which a mode emission is derived
LOGICAL :: lmode_emiss(nmodes) ! True for modes which a given source emission
                               ! emits into
INTEGER :: moment_step = moment_mass - moment_number

INTEGER :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_REGISTER_EMISSION'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Check for availability of UKCA configuration data
IF (.NOT. l_ukca_config_available) THEN
  errcode = 1
  cmessage = 'No UKCA configuration has been set up'
  CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
END IF

IF (ukca_config%l_ukca_emissions_off .OR. (.NOT. ukca_config%l_ukca_chem)) THEN
  errcode = 1
  cmessage = 'Not expecting any emissions to be registered ' //                &
             'for the current UKCA configuration'
  CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
END IF

! Check that the emissions species name is valid for the UKCA configuration
IF (.NOT. ANY(tracer_name == em_chem_spec)) THEN
  errcode = 1
  cmessage = 'The name ' // TRIM(tracer_name) // ' is not an expected ' //     &
             'emission species for the current UKCA configuration'
  CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
END IF

! If the emissions structure is not allocated, allocate using total number
! received from the parent
IF ( .NOT. ALLOCATED(ncdf_emissions) ) THEN
  IF ( num_tot < 0 ) THEN
    errcode = 1
    WRITE(cmessage,'(A,I0)') 'Total number of emissions expected, received ',  &
                             num_tot
    CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
  END IF

  ALLOCATE(ncdf_emissions(num_tot))
  CALL ukca_em_struct_init(ncdf_emissions)
  n_cdf_emiss = num_tot    ! Store for downstream use.
END IF

IF ( .NOT. ALLOCATED(emiss_map) ) THEN
  ALLOCATE(emiss_map(num_tot))
  emiss_map(:) = -99
END IF

n_ems_registered = n_ems_registered + 1  ! Returned index, always in
                                         ! contiguous order
ecount = ecount + 1                      ! Actual position in data structure,
                                         ! next slot
IF ( ecount > n_cdf_emiss ) THEN
  errcode = 2
  WRITE(cmessage,'(A,2(1x,I0))') 'Number of emission slots exceeds maximum',   &
                                 n_cdf_emiss, ecount
  CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
END IF

emiss_map(n_ems_registered) = ecount  ! Store mapping between the two indices
emiss_id = n_ems_registered           ! Index value to be returned to parent

! Assign received values to components in emission structure. Defaults for
! optional components are set in ukca_em_struct_init
ncdf_emissions(ecount)%var_name    = varname
ncdf_emissions(ecount)%tracer_name = tracer_name
ncdf_emissions(ecount)%units       = units
ncdf_emissions(ecount)%three_dim   = three_dim
ncdf_emissions(ecount)%l_update    = .FALSE.

IF (PRESENT(l_omit_std_name))                                                  &
  ncdf_emissions(ecount)%l_omit_std_name = l_omit_std_name
IF (PRESENT(std_name)) ncdf_emissions(ecount)%std_name = std_name
IF (PRESENT(long_name)) ncdf_emissions(ecount)%lng_name = long_name
IF (PRESENT(vert_fact)) ncdf_emissions(ecount)%vert_fact = vert_fact
IF (PRESENT(lowest_lev)) ncdf_emissions(ecount)%lowest_lev = lowest_lev
IF (PRESENT(highest_lev)) ncdf_emissions(ecount)%highest_lev = highest_lev

IF (PRESENT(hourly_fact)) ncdf_emissions(ecount)%hourly_fact = hourly_fact
IF (PRESENT(daily_fact)) ncdf_emissions(ecount)%daily_fact = daily_fact

! If this emission field feeds into aerosols, set up additional slots for these
!  (with same tracer name) in the emissions structure
IF ( ANY(aero_ems_species == tracer_name) ) THEN
  CALL ukca_def_mode_emiss(tracer_name, lmode_emiss, icp)
  from_emiss = ecount   ! Source of emission for the aerosol modes
  DO imode = 1, nmodes
    IF (lmode_emiss(imode)) THEN
      DO imoment = moment_number,moment_mass,moment_step
        ecount = ecount + 1
        ! Set MODE specific attributes, first copying default values from the
        ! 'source' emission; tracer_name is fixed.
        ! For mass emissions units are same as in file, but for number
        ! they need overwritten. std_name/lng_name will be wrong for
        ! number but these are not used in the code.
        ! Number flux is multiplied by MM_DA/AVC to give units of
        ! kg(air) m-2 s-1
        ! If tracer_name contains substring SO2 or biomass, set
        ! logicals which control scaling in ukca_add_emiss
        ncdf_emissions(ecount) = ncdf_emissions(from_emiss)
        ncdf_emissions(ecount)%l_mode      = .TRUE.
        ncdf_emissions(ecount)%l_online    = .TRUE.
        ncdf_emissions(ecount)%mode        = imode
        ncdf_emissions(ecount)%moment      = imoment
        ncdf_emissions(ecount)%component   = icp
        ncdf_emissions(ecount)%tracer_name = 'mode_emiss'
        ncdf_emissions(ecount)%from_emiss  = from_emiss
        ncdf_emissions(ecount)%l_update    = .FALSE.

        IF (imoment == moment_number) THEN
          ncdf_emissions(ecount)%units     = 'kg(air) m-2 s-1'
        ELSE
          ncdf_emissions(ecount)%units     = 'kg m-2 s-1'
        END IF

        IF (INDEX(tracer_name,'SO2') /= 0) THEN
          ncdf_emissions(ecount)%l_mode_so2 = .TRUE.
        ELSE IF (INDEX(tracer_name,'biomass') /= 0 ) THEN
          ncdf_emissions(ecount)%l_mode_biom = .TRUE.
        END IF
      END DO  ! imoment loop
    END IF    ! IF lmode_emiss
  END DO      ! Loop over modes
END IF        ! MODE emission

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_register_emission

! ----------------------------------------------------------------------------
! Description:
!   Returns information about the registered emissions: number of 2D
!   emissions, number of 3D emissions and whether these emissions have
!   been registered in dimension order.
! ----------------------------------------------------------------------------
SUBROUTINE get_registered_ems_info(n_2d_ems, n_3d_ems, l_ndim_order)

IMPLICIT NONE

! Arguments

INTEGER, INTENT(OUT) :: n_2d_ems      ! Number of 2D emissions registered
INTEGER, INTENT(OUT) :: n_3d_ems      ! Number of 3D emissions registered

LOGICAL, INTENT(OUT) :: l_ndim_order  ! T if emissions are registerd in order
                                      ! of their dimensionality

! Local variables

INTEGER :: i
INTEGER :: ecount

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_REGISTERED_EMS_INFO'

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

n_2d_ems = 0
n_3d_ems = 0
l_ndim_order = .TRUE.

DO i = 1, n_ems_registered
  ecount = emiss_map(i)
  IF (ecount > 0) THEN
    IF (ncdf_emissions(ecount)%three_dim) THEN
      n_3d_ems = n_3d_ems + 1
    ELSE
      n_2d_ems = n_2d_ems + 1
      IF (n_3d_ems > 0) l_ndim_order = .FALSE.
    END IF
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE get_registered_ems_info

! ----------------------------------------------------------------------------
! Description:
!  Recives the emission value for a given species, valid for the current
!  timestep and enters this into the emissions data structure
!  Also populates the array of 3-D vertical scaling factors for that species.
!  This is achieved by one of the following methods:
!  1. IF vert_scaling array is Present and Allocated use directly
!  2. IF vert_scaling array is Present and Not Allocated, calculate scaling and
!     pass it back to parent
!  3. If vert_scaling array is absent, the calculation will be done inside UKCA
!
! ---------------------------------------------------------------------------
SUBROUTINE ukca_set_emission (emiss_id, emiss_value, vert_scaling)

USE ukca_fieldname_mod, ONLY: fldname_grid_surf_area
USE ukca_environment_req_mod, ONLY: environ_field_available
USE ukca_environment_fields_mod, ONLY: grid_surf_area

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN) :: emiss_id       ! Index to map to actual position of
                                      ! emission field in struct
REAL,    INTENT(IN) :: emiss_value(:,:,:)  ! Emission value

REAL, ALLOCATABLE, OPTIONAL, INTENT(IN OUT) :: vert_scaling(:,:,:)  ! 3-D
                                                 ! vertical scaling factor

! Local variables

INTEGER :: idx, l, ilev
INTEGER :: em_rowlen, em_rows, em_levs  ! incoming array dimensions

REAL    :: base_scaling    ! factor to convert emiss units
LOGICAL :: gridbox_emiss   ! True for emissions reported per grid box

LOGICAL :: vert_scaling_present_allocated

INTEGER :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage

! Dr Hook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_EMISSION'

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Check that the index given to identify the emission is in range

IF (emiss_id < 1) THEN
  errcode = 1
  WRITE(cmessage,'(A,I0,A)') 'The given emission index ', emiss_id,            &
                             ' is not a positive integer'
  CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
END IF

IF (emiss_id > n_ems_registered) THEN
  errcode = 1
  WRITE(cmessage,'(A,I0,A,I0)')                                                &
    'The given emission index ', emiss_id,                                     &
    ' exceeds the number of registered emissions ', n_ems_registered
  CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
END IF

! Get the emission slot corresponding to this emiss_id
idx = emiss_map(emiss_id)

IF ( idx < 1 .OR. idx > n_cdf_emiss ) THEN
  errcode = 1
  WRITE(cmessage,'(A,I0,A,I0,A,I0,A)')                                         &
    'Bad emissions structure index idx = ', idx,                               &
    ' (emiss_id = ', emiss_id, ', n_cdf_emiss = ', n_cdf_emiss, ')'
  CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
END IF

! Ensure that the particular emission slot has been registered - by checking
! one of the mandatory parameters from the data structure
IF ( SIZE(ncdf_emissions) > 0 .AND. ncdf_emissions(idx)%var_name == ' ' ) THEN
  errcode = 2
  WRITE(cmessage,'(A,I4)')'This slot in emission structure is not registered ',&
                           idx
  CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
END IF

! Update values in the emission structure
em_rowlen = SIZE(emiss_value,DIM=1)
em_rows = SIZE(emiss_value,DIM=2)
em_levs = SIZE(emiss_value,DIM=3)
! Check horizontal dimensions - (will need change if zonal emissions added)
IF (em_rowlen /= ukca_config%row_length .OR. em_rows /= ukca_config%rows) THEN
  errcode = 3
  WRITE(cmessage,'(A,A,4(1x,I0))') 'Array dimensions do not match: ',          &
         ncdf_emissions(idx)%var_name, em_rowlen, ukca_config%row_length,      &
         em_rows, ukca_config%rows
  CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
END IF
IF ( .NOT. ASSOCIATED(ncdf_emissions(idx)%values) ) THEN
  ALLOCATE(ncdf_emissions(idx)%values(em_rowlen, em_rows, em_levs))
END IF

ncdf_emissions(idx)%values(:,:,:) = emiss_value(:,:,:)
ncdf_emissions(idx)%l_update = .TRUE.

! Get conversion factors for the emission field. A base_scaling different
! from 1 will be applied if what is being emitted does not exactly coincide
! with the tracer the emission field is mapped to (e.g. if an NO emission
! field is expressed as kg of nitrogen, or if a VOC emission is expressed as
! kg of carbon).
CALL base_emiss_factors (ncdf_emissions(idx)%l_omit_std_name,                  &
           ncdf_emissions(idx)%tracer_name,                                    &
           ncdf_emissions(idx)%units, ncdf_emissions(idx)%var_name,            &
           ncdf_emissions(idx)%std_name, ncdf_emissions(idx)%lng_name,         &
           base_scaling, gridbox_emiss)

! Get emissions in kg of emitted tracer
ncdf_emissions(idx)%values (:,:,:) = ncdf_emissions(idx)%values (:,:,:) *      &
                                     base_scaling

! If emissions are in kg gridbox-1 s-1 then convert to kg m-2 s-1
IF (gridbox_emiss) THEN
  IF (.NOT. ukca_config%l_support_ems_gridbox_units) THEN
    errcode = 4
    WRITE(cmessage,'(A,I0,A)') 'Offline emission ', emiss_id,                  &
      ' has grid-box units. These are not supported in this configuration.'
    CALL ereport (ModuleName//':'//RoutineName, errcode, cmessage)
  END IF
  IF (.NOT. environ_field_available(fldname_grid_surf_area)) THEN
    errcode = 5
    cmessage = 'Required field ' //fldname_grid_surf_area//                    &
               ' is not available in environment'
    CALL ereport (ModuleName//':'//RoutineName, errcode, cmessage)
  END IF
  IF (ncdf_emissions(idx)%three_dim) THEN      ! 3D emiss
    DO ilev = 1, ukca_config%model_levels
      ncdf_emissions(idx)%values (:,:,ilev) =                                  &
                   ncdf_emissions(idx)%values (:,:,ilev) / grid_surf_area (:,:)
    END DO
  ELSE                                  ! surface or single-level emiss
    ncdf_emissions(idx)%values (:,:,1) =                                       &
                      ncdf_emissions(idx)%values (:,:,1) / grid_surf_area (:,:)
  END IF
END IF

! If this emission feeds into aerosols, allocate the values array for the
! associated slots, where from_emiss = idx. The array will be populated later
! in ukca_emiss_mode_map.
DO l = 1, n_cdf_emiss
  IF ( ncdf_emissions(l)%l_mode .AND. ncdf_emissions(l)%from_emiss == idx      &
      .AND. .NOT. ASSOCIATED(ncdf_emissions(l)%values) ) THEN
    ALLOCATE( ncdf_emissions(l)%values(ukca_config%row_length,                 &
                                       ukca_config%rows, em_levs))
    ncdf_emissions(l)%values(:,:,:) = 0.0
  END IF
END DO

! Vertical scaling factor calculation. If vert_scaling array is present and
! allocated, value is coming from parent so use it directly. If not allocated,
! parent is still supposed to hold the values, so calculate this pass it back
!
! If vert_scaling is not an argument, then this is supposed to be calculated
! and stored internally. However, this will happen in the main emissions routine
! inside UKCA, so cannot be checked here (on the first timestep)

vert_scaling_present_allocated = .FALSE.
IF (PRESENT(vert_scaling)) THEN
  IF (ALLOCATED(vert_scaling)) THEN
    vert_scaling_present_allocated = .TRUE.
  END IF
END IF

IF ( vert_scaling_present_allocated ) THEN
  IF ( MINVAL(vert_scaling) < 0.0 ) THEN
    errcode = 6
    WRITE(cmessage,'(A)') 'Vertical scaling array passed as argument not' //   &
        ' correctly set for ' // ncdf_emissions(idx)%var_name
    CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
  END IF

  IF ( .NOT. ASSOCIATED(ncdf_emissions(idx)%vert_scaling_3d) ) THEN
    ALLOCATE(ncdf_emissions(idx)%vert_scaling_3d(ukca_config%row_length,       &
                                  ukca_config%rows, ukca_config%model_levels))
  END IF
  ncdf_emissions(idx)%vert_scaling_3d(:,:,:) = vert_scaling(:,:,:)

ELSE  ! vert_scaling array not present or not allocated
  ! If internal array is also not allocated calculate scaling factor
  IF ( .NOT. ASSOCIATED(ncdf_emissions(idx)%vert_scaling_3d) ) THEN
    ALLOCATE(ncdf_emissions(idx)%vert_scaling_3d(ukca_config%row_length,       &
                                  ukca_config%rows, ukca_config%model_levels))

    CALL vertical_emiss_factors (ukca_config%row_length, ukca_config%rows,     &
           ukca_config%model_levels, ncdf_emissions(idx)%lowest_lev,           &
           ncdf_emissions(idx)%highest_lev, ncdf_emissions(idx)%vert_fact,     &
           ncdf_emissions(idx)%tracer_name, ncdf_emissions(idx)%vert_scaling_3d)
  END IF
END IF  ! vert_scaling Present and Allocated

IF ( PRESENT(vert_scaling) ) THEN
  IF ( .NOT. ALLOCATED(vert_scaling) ) THEN
    ! Pass back internal or freshly calculated value
    ALLOCATE( vert_scaling(ukca_config%row_length, ukca_config%rows,           &
        ukca_config%model_levels) )
    vert_scaling(:,:,:) = ncdf_emissions(idx)%vert_scaling_3d(:,:,:)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE ukca_set_emission

!---------------------------------------------------------------------
END MODULE ukca_emiss_api_mod
