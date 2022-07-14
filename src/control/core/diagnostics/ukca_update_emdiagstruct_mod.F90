! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT***

!  Description:
!    Update emission diagnostic values and flags in the "emdiags" structure.
!
!  Method:
!    In the first time step initialise all emission diagnostic flags
!    to FALSE in the "emdiags" structure.
!    Identify the given emission diagnostic by the argument diag_name.
!    When it is the  first time we look for it then allocate the
!    corresponding field in the "emdiag" structure and set the flag
!    to TRUE. Then store diagnostic values in that field.
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

MODULE ukca_update_emdiagstruct_mod

USE ukca_emdiags_struct_mod, ONLY: emdiags_struct
USE ukca_emiss_mode_mod,     ONLY: aero_ems_species

IMPLICIT NONE

! Declaration of new emiss diags
TYPE (emdiags_struct) :: emdiags

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName='UKCA_UPDATE_EMDIAGSTRUCT_MOD'

CONTAINS
! --------------------------------------------------------------------------

SUBROUTINE update_emdiagstruct (                                               &
           row_length, rows, model_levels,                                     &
           em_diags, diag_name)

USE ukca_um_legacy_mod,     ONLY: get_emdiag_stash, sf, UKCA_diag_sect
USE ereport_mod,            ONLY: ereport
USE parkind1,               ONLY: jpim,  jprb     ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Subroutine arguments
INTEGER, INTENT(IN)    :: row_length        ! Model dimensions
INTEGER, INTENT(IN)    :: rows
INTEGER, INTENT(IN)    :: model_levels

REAL,    INTENT(IN)    :: em_diags (:,:,:)  ! values of the emiss diagnostic

CHARACTER (LEN=10), INTENT(IN) :: diag_name ! name of the emiss diagnostic

! Local variables
INTEGER                        :: ierr
INTEGER                        :: item, section
CHARACTER (LEN=errormessagelength)            :: cmessage
LOGICAL, SAVE                  :: l_first = .TRUE.

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_EMDIAGSTRUCT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierr    = 0

! Initialise all emission diagnostic flags in 'emdiags' structure
IF (l_first) THEN
  emdiags%l_em_no     = .FALSE.
  emdiags%l_em_ch4    = .FALSE.
  emdiags%l_em_co     = .FALSE.
  emdiags%l_em_hcho   = .FALSE.
  emdiags%l_em_c2h6   = .FALSE.
  emdiags%l_em_c3h8   = .FALSE.
  emdiags%l_em_me2co  = .FALSE.
  emdiags%l_em_mecho  = .FALSE.
  emdiags%l_em_c5h8   = .FALSE.
  emdiags%l_em_c4h10  = .FALSE.
  emdiags%l_em_c2h4   = .FALSE.
  emdiags%l_em_c3h6   = .FALSE.
  emdiags%l_em_tol    = .FALSE.
  emdiags%l_em_oxyl   = .FALSE.
  emdiags%l_em_ch3oh  = .FALSE.
  emdiags%l_em_h2     = .FALSE.
  emdiags%l_em_no_air = .FALSE.
  emdiags%l_em_montrp = .FALSE.
  emdiags%l_em_meoh   = .FALSE.
  emdiags%l_em_nh3    = .FALSE.
  emdiags%l_em_dms    = .FALSE.
  emdiags%l_em_so2low = .FALSE.
  emdiags%l_em_so2hi  = .FALSE.
  emdiags%l_em_so2nat = .FALSE.
  ! CRI emission logicals
  emdiags%l_em_etoh    = .FALSE.
  emdiags%l_em_c2h2    = .FALSE.
  emdiags%l_em_tbut2ene= .FALSE.
  emdiags%l_em_apinene = .FALSE.
  emdiags%l_em_bpinene = .FALSE.
  emdiags%l_em_benzene = .FALSE.
  emdiags%l_em_hcooh   = .FALSE.
  emdiags%l_em_meco2h  = .FALSE.
  emdiags%l_em_etcho   = .FALSE.
  emdiags%l_em_hoch2cho= .FALSE.
  emdiags%l_em_mek     = .FALSE.

  l_first             = .FALSE.
END IF

! -----------------------------------------------------------------------
! Store emiss diagnostics in the corresponding field of the 'emdiags'
! structure. Note that most fields are 2D. However 3D fields are also
! allowed as it is the case for 'NO_aircrft'.
! If it is the first time we look for a given diagnostics then also
! set the emission diagnostic flag to TRUE. The detailed comments
! for the first diagnostics ('NO        ') in the CASE statement
! below are valid for any other emission diagnostic.
section = UKCA_diag_sect
IF (.NOT. ANY(aero_ems_species == diag_name)) THEN
  item = get_emdiag_stash (diag_name)
ELSE
  item = 0  ! Dummy value; Item No. is not required for aerosol emissions here
END IF

SELECT CASE ( diag_name )
CASE ('NO        ')
  IF (emdiags%l_em_no) THEN
    ! If the emission diagnostic flag (emdiags%l_em_no) is set to TRUE
    ! then the array with emission diagnostics (emdiags%em_no) has already
    ! been allocated. We only need to fill it in with the values 'em_diags'
    ! passed to this routine as an INTEN(IN) argument.
    emdiags%em_no (:,:) = em_diags (:,:,1)

  ELSE
    IF (sf(item,section)) THEN
      ! If the NO diagnostic is requested by STASH for the first time, then
      ! set the flag emdiags%l_em_no to TRUE and allocate the array
      ! emdiags%em_no before filling it with the values 'em_diags'.
      ALLOCATE (emdiags%em_no (row_length, rows))
      emdiags%em_no (:,:) = em_diags (:,:,1)
      emdiags%l_em_no     = .TRUE.
    END IF
  END IF

CASE ('CH4       ')
  IF (emdiags%l_em_ch4) THEN
    emdiags%em_ch4 (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_ch4 (row_length, rows))
      emdiags%em_ch4 (:,:) = em_diags (:,:,1)
      emdiags%l_em_ch4     = .TRUE.
    END IF
  END IF

CASE ('CO        ')
  IF (emdiags%l_em_co) THEN
    emdiags%em_co (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_co (row_length, rows))
      emdiags%em_co (:,:) = em_diags (:,:,1)
      emdiags%l_em_co     = .TRUE.
    END IF
  END IF

CASE ('HCHO      ')
  IF (emdiags%l_em_hcho) THEN
    emdiags%em_hcho (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_hcho (row_length, rows))
      emdiags%em_hcho (:,:) = em_diags (:,:,1)
      emdiags%l_em_hcho     = .TRUE.
    END IF
  END IF

CASE ('C2H6      ')
  IF (emdiags%l_em_c2h6) THEN
    emdiags%em_c2h6 (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_c2h6 (row_length, rows))
      emdiags%em_c2h6 (:,:) = em_diags (:,:,1)
      emdiags%l_em_c2h6     = .TRUE.
    END IF
  END IF

CASE ('C3H8      ')
  IF (emdiags%l_em_c3h8) THEN
    emdiags%em_c3h8 (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_c3h8 (row_length, rows))
      emdiags%em_c3h8 (:,:) = em_diags (:,:,1)
      emdiags%l_em_c3h8     = .TRUE.
    END IF
  END IF

CASE ('Me2CO     ')
  IF (emdiags%l_em_me2co) THEN
    emdiags%em_me2co (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_me2co (row_length, rows))
      emdiags%em_me2co (:,:) = em_diags (:,:,1)
      emdiags%l_em_me2co     = .TRUE.
    END IF
  END IF

CASE ('MeCHO     ')
  IF (emdiags%l_em_mecho) THEN
    emdiags%em_mecho (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_mecho (row_length, rows))
      emdiags%em_mecho (:,:) = em_diags (:,:,1)
      emdiags%l_em_mecho     = .TRUE.
    END IF
  END IF

CASE ('C5H8      ')
  IF (emdiags%l_em_c5h8) THEN
    emdiags%em_c5h8 (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_c5h8 (row_length, rows))
      emdiags%em_c5h8 (:,:) = em_diags (:,:,1)
      emdiags%l_em_c5h8     = .TRUE.
    END IF
  END IF

CASE ('C4H10     ')
  IF (emdiags%l_em_c4h10) THEN
    emdiags%em_c4h10 (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_c4h10 (row_length, rows))
      emdiags%em_c4h10 (:,:) = em_diags (:,:,1)
      emdiags%l_em_c4h10     = .TRUE.
    END IF
  END IF

CASE ('C2H4      ')
  IF (emdiags%l_em_c2h4) THEN
    emdiags%em_c2h4 (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_c2h4 (row_length, rows))
      emdiags%em_c2h4 (:,:) = em_diags (:,:,1)
      emdiags%l_em_c2h4     = .TRUE.
    END IF
  END IF

CASE ('C3H6      ')
  IF (emdiags%l_em_c3h6) THEN
    emdiags%em_c3h6 (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_c3h6 (row_length, rows))
      emdiags%em_c3h6 (:,:) = em_diags (:,:,1)
      emdiags%l_em_c3h6     = .TRUE.
    END IF
  END IF

CASE ('TOLUENE   ')
  IF (emdiags%l_em_tol) THEN
    emdiags%em_tol (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_tol (row_length, rows))
      emdiags%em_tol (:,:) = em_diags (:,:,1)
      emdiags%l_em_tol     = .TRUE.
    END IF
  END IF

CASE ('oXYLENE   ')
  IF (emdiags%l_em_oxyl) THEN
    emdiags%em_oxyl (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_oxyl (row_length, rows))
      emdiags%em_oxyl (:,:) = em_diags (:,:,1)
      emdiags%l_em_oxyl     = .TRUE.
    END IF
  END IF

CASE ('CH3OH     ')
  IF (emdiags%l_em_ch3oh) THEN
    emdiags%em_ch3oh (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_ch3oh (row_length, rows))
      emdiags%em_ch3oh (:,:) = em_diags (:,:,1)
      emdiags%l_em_ch3oh     = .TRUE.
    END IF
  END IF

CASE ('H2        ')
  IF (emdiags%l_em_h2) THEN
    emdiags%em_h2 (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_h2 (row_length, rows))
      emdiags%em_h2 (:,:) = em_diags (:,:,1)
      emdiags%l_em_h2     = .TRUE.
    END IF
  END IF

CASE ('NO_aircrft')
  IF (emdiags%l_em_no_air) THEN
    emdiags%em_no_air (:,:,:) = em_diags (:,:,:)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_no_air (row_length, rows, model_levels))
      emdiags%em_no_air (:,:,:) = em_diags (:,:,:)
      emdiags%l_em_no_air       = .TRUE.
    END IF
  END IF

CASE ('Monoterp  ')
  IF (emdiags%l_em_montrp) THEN
    emdiags%em_montrp (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_montrp (row_length, rows))
      emdiags%em_montrp (:,:) = em_diags (:,:,1)
      emdiags%l_em_montrp       = .TRUE.
    END IF
  END IF

CASE ('NVOC      ','MeOH      ')
  IF (emdiags%l_em_meoh) THEN
    emdiags%em_meoh (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_meoh (row_length, rows))
      emdiags%em_meoh (:,:) = em_diags (:,:,1)
      emdiags%l_em_meoh       = .TRUE.
    END IF
  END IF

CASE ('NH3       ')
  IF (emdiags%l_em_nh3) THEN
    emdiags%em_nh3 (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_nh3 (row_length, rows))
      emdiags%em_nh3 (:,:) = em_diags (:,:,1)
      emdiags%l_em_nh3       = .TRUE.
    END IF
  END IF

CASE ('DMS       ')
  IF (emdiags%l_em_dms) THEN
    emdiags%em_dms (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_dms (row_length, rows))
      emdiags%em_dms (:,:) = em_diags (:,:,1)
      emdiags%l_em_dms       = .TRUE.
    END IF
  END IF

CASE ('SO2_low   ')
  IF (emdiags%l_em_so2low) THEN
    emdiags%em_so2low (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_so2low (row_length, rows))
      emdiags%em_so2low (:,:) = em_diags (:,:,1)
      emdiags%l_em_so2low       = .TRUE.
    END IF
  END IF

CASE ('SO2_high   ')
  IF (emdiags%l_em_so2hi) THEN
    emdiags%em_so2hi (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_so2hi (row_length, rows))
      emdiags%em_so2hi (:,:) = em_diags (:,:,1)
      emdiags%l_em_so2hi       = .TRUE.
    END IF
  END IF

CASE ('SO2_nat    ')
  IF (emdiags%l_em_so2nat) THEN
    emdiags%em_so2nat (:,:,:) = em_diags (:,:,:)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_so2nat (row_length, rows, model_levels))
      emdiags%em_so2nat (:,:,:) = em_diags (:,:,:)
      emdiags%l_em_so2nat       = .TRUE.
    END IF
  END IF

  ! CASE selects for CRI emission diagnostics
CASE ('EtOH       ')
  IF (emdiags%l_em_etoh) THEN
    emdiags%em_etoh (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_etoh (row_length, rows))
      emdiags%em_etoh (:,:)  = em_diags (:,:,1)
      emdiags%l_em_etoh      = .TRUE.
    END IF
  END IF

CASE ('C2H2       ')
  IF (emdiags%l_em_c2h2) THEN
    emdiags%em_c2h2 (:,:) = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_c2h2 (row_length, rows))
      emdiags%em_c2h2 (:,:)  = em_diags (:,:,1)
      emdiags%l_em_c2h2      = .TRUE.
    END IF
  END IF

CASE ('TBUT2ENE   ')
  IF (emdiags%l_em_tbut2ene) THEN
    emdiags%em_tbut2ene (:,:)   = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_tbut2ene (row_length, rows))
      emdiags%em_tbut2ene (:,:) = em_diags (:,:,1)
      emdiags%l_em_tbut2ene     = .TRUE.
    END IF
  END IF

CASE ('APINENE    ')
  IF (emdiags%l_em_apinene) THEN
    emdiags%em_apinene (:,:)    = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_apinene (row_length, rows))
      emdiags%em_apinene (:,:)  = em_diags (:,:,1)
      emdiags%l_em_apinene      = .TRUE.
    END IF
  END IF

CASE ('BPINENE    ')
  IF (emdiags%l_em_bpinene) THEN
    emdiags%em_bpinene (:,:)    = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_bpinene (row_length, rows))
      emdiags%em_bpinene (:,:)  = em_diags (:,:,1)
      emdiags%l_em_bpinene      = .TRUE.
    END IF
  END IF

CASE ('BENZENE    ')
  IF (emdiags%l_em_benzene) THEN
    emdiags%em_benzene (:,:)    = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_benzene (row_length, rows))
      emdiags%em_benzene (:,:)  = em_diags (:,:,1)
      emdiags%l_em_benzene      = .TRUE.
    END IF
  END IF

CASE ('HCOOH      ')
  IF (emdiags%l_em_hcooh) THEN
    emdiags%em_hcooh (:,:)      = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_hcooh (row_length, rows))
      emdiags%em_hcooh (:,:)    = em_diags (:,:,1)
      emdiags%l_em_hcooh        = .TRUE.
    END IF
  END IF

CASE ('MeCO2H     ')
  IF (emdiags%l_em_meco2h) THEN
    emdiags%em_meco2h (:,:)     = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_meco2h (row_length, rows))
      emdiags%em_meco2h (:,:)   = em_diags (:,:,1)
      emdiags%l_em_meco2h       = .TRUE.
    END IF
  END IF

CASE ('EtCHO      ')
  IF (emdiags%l_em_etcho) THEN
    emdiags%em_etcho (:,:)      = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_etcho (row_length, rows))
      emdiags%em_etcho (:,:)    = em_diags (:,:,1)
      emdiags%l_em_etcho        = .TRUE.
    END IF
  END IF

CASE ('HOCH2CHO   ')
  IF (emdiags%l_em_hoch2cho) THEN
    emdiags%em_hoch2cho (:,:)   = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_hoch2cho (row_length, rows))
      emdiags%em_hoch2cho (:,:) = em_diags (:,:,1)
      emdiags%l_em_hoch2cho     = .TRUE.
    END IF
  END IF

CASE ('MEK        ')
  IF (emdiags%l_em_mek) THEN
    emdiags%em_mek   (:,:)      = em_diags (:,:,1)
  ELSE
    IF (sf(item,section)) THEN
      ALLOCATE (emdiags%em_mek (row_length, rows))
      emdiags%em_mek   (:,:)    = em_diags (:,:,1)
      emdiags%l_em_mek          = .TRUE.
    END IF
  END IF


CASE DEFAULT
    ! Report error unless this is an aerosol emission (BC_fossil:, etc)
  IF (.NOT. ANY(aero_ems_species == diag_name))  ierr = 1
END SELECT

! Report error if diagnostics not found
IF (ierr == 1) THEN
  cmessage = 'Unexpected diagnostic in S-50: ' // TRIM (diag_name)
  CALL ereport ('UPDATE_EMDIAGSTRUCT', ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE update_emdiagstruct

END MODULE ukca_update_emdiagstruct_mod
