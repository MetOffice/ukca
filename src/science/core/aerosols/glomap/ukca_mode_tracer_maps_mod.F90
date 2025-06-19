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
!    Module to store indices of MODE data in UM and UKCA tracer arrays

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  Fortran
!
! ######################################################################
!
! Subroutine Interface:
MODULE ukca_mode_tracer_maps_mod

IMPLICIT NONE

! Index of tracers for mode number, allocated size is the number of modes.
! Holds the tracer number in the all_tracers array (internal to UKCA).
! Defaults to zero if the mode is not found.
INTEGER, ALLOCATABLE, SAVE :: nmr_index(:)

! Index of tracers for mass mixing ratio of component, allocated size is
! no. of modes * no. of components.
! Holds the appropriate tracer number in the all_tracers array.
! Defaults to -1 if the component is not found.
INTEGER, ALLOCATABLE, SAVE :: mmr_index(:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: modulename='UKCA_MODE_TRACER_MAPS_MOD'


CONTAINS


SUBROUTINE ukca_aero_tracer_init()

!   Initialise arrays used to identify mode tracers.

!   Procedure:  The arrays nmr_index(nmodes) and mmr_index(nmodes,ncp) index
!   the number and component tracers in the all_tracers array. All
!   UKCA tracers are checked in turn to determine if they are a number
!   or mass tracer and if so mmr_index and nmr_index are then filled with the
!   tracer number for the enabled tracers in the ukca tracer array.
!   So mmr_index(2,3) gives the index of the tracer for the second mode and
!   third component.

USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength
USE parkind1,                  ONLY: jprb, jpim
USE ukca_tracers_mod,          ONLY: all_tracers_names, n_tracers

USE ukca_config_specification_mod, ONLY: glomap_variables

USE ukca_mode_setup,           ONLY: nmodes, mode_names, cp_su, cp_cl,         &
                                     cp_bc, cp_oc, cp_du, cp_so, cp_nh4,       &
                                     cp_no3, cp_nn

USE ukca_um_legacy_mod,        ONLY: mype
USE umprintmgr,                ONLY: umprint, ummessage, printstatus,          &
                                     prstatus_diag
USE yomhook,                   ONLY: lhook, dr_hook
USE parkind1,                  ONLY: jprb, jpim

IMPLICIT NONE

! Local variables

INTEGER :: icp                  ! component index
INTEGER :: imode                ! mode index
INTEGER :: icode                ! error code
INTEGER :: i, j                 ! loop counters
#if defined(NAG_FORTRAN) && (NAG_FORTRAN == 7000000)
LOGICAL :: tracer_belongs
#endif
CHARACTER(LEN=errormessagelength) :: cmessage   ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: routinename='UKCA_AERO_TRACER_INIT'

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_in,zhook_handle)

! Allocate the nmr_index and initialise to
! zero. Zero inidicates that a mode
! is not active
IF (.NOT. ALLOCATED(nmr_index)) THEN
  ALLOCATE(nmr_index(nmodes))
  nmr_index(:) = 0
END IF

! Allocate mmr_index and initialise to -1
! -1 indicates that the component is not active
! It is only used if > 0
IF (.NOT. ALLOCATED(mmr_index)) THEN
  ALLOCATE(mmr_index(nmodes,glomap_variables%ncp))
  mmr_index(:,:) = -1
END IF

! Construct the number and component mass index arrays
! by looping over all entries in the all_tracers array
DO i=1,n_tracers
  ! (Slow) Workaround for NAG Fortran vn7.0 internal compiler error.
  ! Fixed in vn7.1
#if defined(NAG_FORTRAN) && (NAG_FORTRAN == 7000000)
  tracer_belongs = .FALSE.
  DO j=1,nmodes
    IF ( mode_names(j) == all_tracers_names(i)(1:7) ) THEN
      tracer_belongs = .TRUE.
    END IF
  END DO
  IF ( tracer_belongs ) THEN
#else
  IF ( ANY(mode_names(:) == all_tracers_names(i)(1:7)) ) THEN
#endif
      ! Tracer belongs to UKCA_MODE aerosol scheme

      ! first 7 characters of tracer name are the mode name
      ! imode is the position in the mode_names array
    imode = arrayloc(mode_names, all_tracers_names(i)(1:7))

    ! last 2 characters are the component name
    ! or 'N ' for number mixing ratio
    SELECT CASE (all_tracers_names(i)(9:10))
    CASE ('N ')
      icp = 0
    CASE ('SU')
      icp = cp_su
    CASE ('SS')
      icp = cp_cl
    CASE ('BC')
      icp = cp_bc
    CASE ('OM')
      icp = cp_oc
    CASE ('DU')
      icp = cp_du
    CASE ('SO')
      icp = cp_so
    CASE ('NH')
      icp = cp_nh4
    CASE ('NT')
      icp = cp_no3
    CASE ('NN')
      icp = cp_nn
    CASE DEFAULT
      icode = 1
      cmessage = ' Unknown component in CASE statement'
      WRITE(umMessage,'(A40,A10)') cmessage,all_tracers_names(i)
      CALL umPrint(umMessage,src='ukca_mode_setup')
      CALL ereport(routinename,icode,cmessage)
    END SELECT

    ! now use this information to set up the indices in
    ! all tracers for number and mass
    IF (icp == 0) THEN
      ! This is a number density - set index in nmr_index
      nmr_index(imode) = i

    ELSE IF (icp > 0) THEN
      ! This is a component mass - set index in mmr_index
      mmr_index(imode,icp) = i

    END IF
  END IF
END DO

! print out nmr and mmr indices
IF (printstatus >= prstatus_diag .AND. mype == 0) THEN
  WRITE(ummessage,'(A)') 'number tracer maps'
  CALL umprint(ummessage,src=routinename)
  DO i=1, nmodes
    WRITE(ummessage,'(I6,1X,I6)') i, nmr_index(i)
    CALL umPrint(ummessage,src=routinename)
  END DO
END IF

IF (printstatus >= prstatus_diag .AND. mype == 0) THEN
  WRITE(umMessage,'(A)') 'mass tracer maps'
  CALL umprint(ummessage,src=routinename)
  DO j=1, glomap_variables%ncp
    DO i=1,nmodes
      WRITE(ummessage,'(I6,1X,I6,1X,I6)') i, j, mmr_index(i, j)
      CALL umprint(ummessage,src=routinename)
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_aero_tracer_init


INTEGER FUNCTION arrayloc(ch_arr,ch_test)

! Takes as input a 1D CHARACTER arrray and
! a CHARACTER variable
! If the value is in the array, return the index,
! if it isn't return -99
! Note that the match is case sensitive but ignores
! all leading and trailing white space
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

INTEGER :: i
CHARACTER(LEN=*),INTENT(IN) :: ch_arr(:)
CHARACTER(LEN=*),INTENT(IN) :: ch_test

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: routinename='ARRAYLOC'

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_in,zhook_handle)

! Start with default value
arrayloc = -99

find_arrayloc : DO i = 1, SIZE (ch_arr)
  IF (TRIM(ADJUSTL(ch_arr(i))) == TRIM(ADJUSTL(ch_test))) THEN
    arrayloc = i
    EXIT find_arrayloc
  END IF
END DO find_arrayloc

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_out,zhook_handle)

END FUNCTION arrayloc

END MODULE ukca_mode_tracer_maps_mod
