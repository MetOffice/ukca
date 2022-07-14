! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module holding the data structure and procedures for initialising emissions
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
MODULE ukca_emiss_struct_mod
!

IMPLICIT NONE

! Emission Data structure

INTEGER, PARAMETER, PUBLIC :: maxlen_emiss_var_name = 80
INTEGER, PARAMETER, PUBLIC :: maxlen_emiss_tracer_name = 10
INTEGER, PARAMETER, PUBLIC :: maxlen_emiss_std_name = 256
INTEGER, PARAMETER, PUBLIC :: maxlen_emiss_long_name = 256
INTEGER, PARAMETER, PUBLIC :: maxlen_emiss_units = 30
INTEGER, PARAMETER, PUBLIC :: maxlen_emiss_hourly_fact = 20
INTEGER, PARAMETER, PUBLIC :: maxlen_emiss_daily_fact = 20
INTEGER, PARAMETER, PUBLIC :: maxlen_emiss_vert_fact = 30

TYPE, PUBLIC :: ukca_em_struct

  CHARACTER (LEN=maxlen_emiss_var_name)                                        &
                       :: var_name     ! Name of variable in file
  CHARACTER (LEN=maxlen_emiss_tracer_name)                                     &
                       :: tracer_name  ! Emitted species

  LOGICAL              :: l_omit_std_name
                                       ! True to allow omission of any
                                       ! descriptive name (implies no need
                                       ! for name-based unit conversion)

  CHARACTER (LEN=maxlen_emiss_std_name)                                        &
                       :: std_name     ! Standard name
  CHARACTER (LEN=maxlen_emiss_long_name)                                       &
                       :: lng_name     ! Long name (alternative to std_name)
  CHARACTER (LEN=maxlen_emiss_units)                                           &
                       :: units        ! Units of the emission field

  LOGICAL              :: l_update     ! True if field updated in a given tstep
  ! (to indicate if conversion factors need to be applied in UKCA_NEW_EMISS_CTL)
  LOGICAL              :: l_online     ! True if field is calculated online
                                       ! (not read from file)
  LOGICAL              :: three_dim    ! True if 3D emiss field

  REAL                 :: base_fact    ! Base conv factor

  CHARACTER (LEN=maxlen_emiss_hourly_fact)                                     &
                       :: hourly_fact  ! hourly_scaling: traffic, none,
                                       ! diurnal_isopems
  CHARACTER (LEN=maxlen_emiss_daily_fact)                                      &
                       :: daily_fact   ! daily_scaling:  traffic, none
  CHARACTER (LEN=maxlen_emiss_var_name)                                        &
                       :: vert_fact    ! vertical_scaling: surface,
                                       ! all_levels, etc.
  INTEGER              :: lowest_lev   ! Lowest and highest level where
  INTEGER              :: highest_lev  ! emiss can be injected

  ! Variables for aerosol emissions...
  LOGICAL              :: l_mode         ! Field is a num/mass mode emiss into
                                         ! which offline emissions are mapped.
  LOGICAL              :: l_mode_so2     ! Mode emission from SO2
  LOGICAL              :: l_mode_biom    ! Mode emission from biomass burning

  ! If l_mode = .TRUE. the following define the moment, mode, component the
  ! field applies to, and the name of the tracer which supplied the source
  ! emissions field.
  INTEGER              :: moment         ! Modal emission moment (0 or 3)
  INTEGER              :: mode           ! Modal emission mode
  INTEGER              :: component      ! Modal emission component
  INTEGER              :: from_emiss     ! Modal source emission index

  ! Allocatable fields of derived types are not allowed in F95 (they are
  ! an extension of Standard F95). As a consequence we use pointers for
  ! variables that need to be allocated within the emissions structure.
  REAL, POINTER:: values (:,:,:) => NULL()  ! emission data
  REAL, POINTER:: diags  (:,:,:) => NULL()  ! emission diagnostics
  REAL, POINTER:: vert_scaling_3d (:,:,:) => NULL()! Vertical conv fact

END TYPE ukca_em_struct

! Super array of emissions
!   File-based, data received from parent
TYPE (ukca_em_struct), ALLOCATABLE, PUBLIC :: ncdf_emissions (:)

! Subroutines available outside this module
PUBLIC :: ukca_em_struct_init, ukca_copy_emiss_struct

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EMISS_STRUCT_MOD'

CONTAINS

!---------------------------------------------------------------------------
! Description:
!   Initialisation of emissions structure to default values.
! Method:
!   Initialise elements of the emissions structure, which has been
!   previously allocated with upper bound equal to 'num_em_flds'.
!   Note that as indicated below some elements cannot be
!   allocated yet within this subroutine.
!---------------------------------------------------------------------------
SUBROUTINE ukca_em_struct_init(emissions)

USE parkind1,             ONLY: jpim, jprb      ! DrHook
USE yomhook,              ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Emissions structure input, as automatic array
TYPE (ukca_em_struct), INTENT(IN OUT) :: emissions (:)

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! For Dr Hook
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EM_STRUCT_INIT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

emissions(:)%var_name    = REPEAT(' ',80)
emissions(:)%tracer_name = REPEAT(' ',10)
emissions(:)%l_omit_std_name = .FALSE.
emissions(:)%std_name    = REPEAT(' ',256)
emissions(:)%lng_name    = REPEAT(' ',256)
emissions(:)%units       = REPEAT(' ',30)

emissions(:)%l_update    = .FALSE.
emissions(:)%l_online    = .FALSE.
emissions(:)%three_dim   = .FALSE.
emissions(:)%vert_fact   = REPEAT(' ',30)
emissions(:)%lowest_lev  = 1   ! Initially assume surface emiss (changed later
emissions(:)%highest_lev = 1   ! in the code for both 3D and high level emiss)

emissions(:)%base_fact   = 1.0 ! by default do not apply any conversion factor

emissions(:)%hourly_fact = REPEAT(' ',20)
emissions(:)%daily_fact  = REPEAT(' ',20)

emissions(:)%l_mode      = .FALSE.
emissions(:)%l_mode_so2  = .FALSE.
emissions(:)%l_mode_biom = .FALSE.
emissions(:)%moment      = -1
emissions(:)%mode        = -1
emissions(:)%component   = -1
emissions(:)%from_emiss  = -1


! Note that the pointers emissions(:)%values and emissions(:)%diags
! cannot be initialised yet. Before that they need to be allocated
! independently for each element of emissions(:). This is done in
! ukca_emiss_init after calling this routine.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_em_struct_init


!-------------------------------------------------------------------------
!
! Copy components of structure from one instance to another
! For multi-dimension fields, allocates the arrays in target instance first
!
! ------------------------------------------------------------------------
SUBROUTINE ukca_copy_emiss_struct(srce, targ)

USE parkind1,             ONLY: jpim, jprb      ! DrHook
USE yomhook,              ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Emissions structure input and output, as automatic arrays
TYPE (ukca_em_struct), INTENT(IN)  :: srce
TYPE (ukca_em_struct), INTENT(OUT) :: targ

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! For Dr Hook
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_COPY_EMISS_STRUCT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

targ%var_name    = srce%var_name
targ%tracer_name = srce%tracer_name
targ%std_name    = srce%std_name
targ%lng_name    = srce%lng_name
targ%units       = srce%units

targ%l_update    = srce%l_update
targ%l_online    = srce%l_online
targ%three_dim   = srce%three_dim
targ%vert_fact   = srce%vert_fact
targ%lowest_lev  = srce%lowest_lev
targ%highest_lev = srce%highest_lev

targ%base_fact   = srce%base_fact

targ%hourly_fact = srce%hourly_fact
targ%daily_fact  = srce%daily_fact

targ%l_mode      = srce%l_mode
targ%l_mode_so2  = srce%l_mode_so2
targ%l_mode_biom = srce%l_mode_biom
targ%moment      = srce%moment
targ%mode        = srce%mode
targ%component   = srce%component
targ%from_emiss  = srce%from_emiss

! For array/ pointer components, clear target arrays and re-allocate to sizes
!  in source before copying
IF (ASSOCIATED(srce%values)) THEN
  IF (ASSOCIATED(targ%values)) DEALLOCATE(targ%values)
  NULLIFY(targ%values)
  ALLOCATE( targ%values(SIZE(srce%values,DIM=1), SIZE(srce%values,DIM=2),      &
                        SIZE(srce%values,DIM=3)) )
  targ%values(:,:,:) = srce%values(:,:,:)
END IF

IF (ASSOCIATED(srce%vert_scaling_3d)) THEN
  IF (ASSOCIATED(targ%vert_scaling_3d)) DEALLOCATE(targ%vert_scaling_3d)
  NULLIFY(targ%vert_scaling_3d)
  ALLOCATE( targ%vert_scaling_3d( SIZE(srce%vert_scaling_3d,DIM=1),            &
                                  SIZE(srce%vert_scaling_3d,DIM=2),            &
                                  SIZE(srce%vert_scaling_3d,DIM=3)) )
  targ%vert_scaling_3d(:,:,:) = srce%vert_scaling_3d(:,:,:)
END IF

IF (ASSOCIATED(srce%diags)) THEN
  IF (ASSOCIATED(targ%diags)) DEALLOCATE(targ%diags)
  NULLIFY(targ%diags)
  ALLOCATE( targ%diags(SIZE(srce%diags,DIM=1), SIZE(srce%diags,DIM=2),         &
                       SIZE(srce%diags,DIM=3)) )
  targ%diags(:,:,:) = srce%diags(:,:,:)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_copy_emiss_struct

END MODULE ukca_emiss_struct_mod
