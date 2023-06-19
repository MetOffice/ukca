! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE ukca_emiss_diags_mode_mod

USE parkind1, ONLY: jprb,  jpim

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
PRIVATE
PUBLIC :: ukca_emiss_diags_mode

!  Description:
!    Produce emission diagnostics for GLOMAP-mode emissions under
!    the new UKCA emission system (based on NetCDF input emissions).
!
!  Method:
!    Information about all available diagnostics is contained in the
!    mode_diag array, which links each STASH item code to mode and component.
!    This array is filled by ukca_emiss_diags_mode_init on the first timestep.
!    The routine ukca_emiss_diags_mode loops over each diagnostic, finds
!    corresponding mass emissions and copies the diags elements of this into a
!    local array em_diags. em_diags is then copied to stash in copydiag_3d.
!    Marine OC emissions are separated by name.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

TYPE :: type_mode_diag_struct
  INTEGER :: item  ! stash item number
  INTEGER :: mode  ! diag's mode
  INTEGER :: component  ! diag's component
  CHARACTER(LEN=30) :: vname  ! name of variable in file
END TYPE type_mode_diag_struct

INTEGER, PARAMETER :: num_mode_diags = 24  ! 13 + 8 for nitrate scheme
                                           !    + 2 for marine OC
                                           !    + 1 for supins dust
TYPE(type_mode_diag_struct) :: mode_diag(num_mode_diags)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EMISS_DIAGS_MODE_MOD'

CONTAINS


SUBROUTINE ukca_emiss_diags_mode_init

USE ukca_mode_setup,  ONLY:                                                    &
    cp_su,                                                                     &
    cp_bc,                                                                     &
    cp_oc,                                                                     &
    cp_cl,                                                                     &
    cp_du,                                                                     &
    cp_no3,                                                                    &
    cp_nh4,                                                                    &
    cp_nn,                                                                     &
    mode_ait_sol,                                                              &
    mode_acc_sol,                                                              &
    mode_cor_sol,                                                              &
    mode_ait_insol,                                                            &
    mode_acc_insol,                                                            &
    mode_cor_insol,                                                            &
    mode_sup_insol

USE ukca_emiss_mod,  ONLY: marine_oc_online
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


IMPLICIT NONE

INTEGER :: icount
CHARACTER (LEN=errormessagelength) :: cmessage

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_DIAGS_MODE_INIT'

INTEGER(KIND=jpim), PARAMETER :: zhook_in=0
INTEGER(KIND=jpim), PARAMETER :: zhook_out=1
REAL   (KIND=jprb)            :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icount = 0
mode_diag(:)%item = -1
mode_diag(:)%mode = -1
mode_diag(:)%component = -1
mode_diag(:)%vname = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'

! SO4 to aitken-sol
icount = icount + 1
mode_diag(icount)%item = 201
mode_diag(icount)%mode = mode_ait_sol
mode_diag(icount)%component = cp_su

! SO4 to accum-sol
icount = icount + 1
mode_diag(icount)%item = 202
mode_diag(icount)%mode = mode_acc_sol
mode_diag(icount)%component = cp_su

! SO4 to coarse-sol
icount = icount + 1
mode_diag(icount)%item = 203
mode_diag(icount)%mode = mode_cor_sol
mode_diag(icount)%component = cp_su

! sea-salt to accum-sol
icount = icount + 1
mode_diag(icount)%item = 204
mode_diag(icount)%mode = mode_acc_sol
mode_diag(icount)%component = cp_cl

! sea-salt to coarse-sol
icount = icount + 1
mode_diag(icount)%item = 205
mode_diag(icount)%mode = mode_cor_sol
mode_diag(icount)%component = cp_cl

! black carbon to aitken-sol
icount = icount + 1
mode_diag(icount)%item = 206
mode_diag(icount)%mode = mode_ait_sol
mode_diag(icount)%component = cp_bc

! black carbon to aitken-ins
icount = icount + 1
mode_diag(icount)%item = 207
mode_diag(icount)%mode = mode_ait_insol
mode_diag(icount)%component = cp_bc

! organic carbon to aitken-sol
icount = icount + 1
mode_diag(icount)%item = 208
mode_diag(icount)%mode = mode_ait_sol
mode_diag(icount)%component = cp_oc

! organic carbon to aitken-ins
icount = icount + 1
mode_diag(icount)%item = 209
mode_diag(icount)%mode = mode_ait_insol
mode_diag(icount)%component = cp_oc

! marine OC to Aitken-sol
icount = icount + 1
mode_diag(icount)%item = 388
mode_diag(icount)%mode = mode_ait_sol
mode_diag(icount)%component = cp_oc
mode_diag(icount)%vname = marine_oc_online

! marine OC to Aitken-insol
icount = icount + 1
mode_diag(icount)%item = 389
mode_diag(icount)%mode = mode_ait_insol
mode_diag(icount)%component = cp_oc
mode_diag(icount)%vname = marine_oc_online

! dust to accum-sol
icount = icount + 1
mode_diag(icount)%item = 210
mode_diag(icount)%mode = mode_acc_sol
mode_diag(icount)%component = cp_du

! dust to accum-ins
icount = icount + 1
mode_diag(icount)%item = 211
mode_diag(icount)%mode = mode_acc_insol
mode_diag(icount)%component = cp_du

! dust to coarse-sol
icount = icount + 1
mode_diag(icount)%item = 212
mode_diag(icount)%mode = mode_cor_sol
mode_diag(icount)%component = cp_du

! dust to coarse-ins
icount = icount + 1
mode_diag(icount)%item = 213
mode_diag(icount)%mode = mode_cor_insol
mode_diag(icount)%component = cp_du

! dust to sup-ins
icount = icount + 1
mode_diag(icount)%item = 675
mode_diag(icount)%mode = mode_sup_insol
mode_diag(icount)%component = cp_du

! ammonium to aitken-sol
icount = icount + 1
mode_diag(icount)%item = 575
mode_diag(icount)%mode = mode_ait_sol
mode_diag(icount)%component = cp_nh4

! ammonium to accum-sol
icount = icount + 1
mode_diag(icount)%item = 576
mode_diag(icount)%mode = mode_acc_sol
mode_diag(icount)%component = cp_nh4

! ammonium to coarse-sol
icount = icount + 1
mode_diag(icount)%item = 577
mode_diag(icount)%mode = mode_cor_sol
mode_diag(icount)%component = cp_nh4

! nitrate to aitken-sol
icount = icount + 1
mode_diag(icount)%item = 579
mode_diag(icount)%mode = mode_ait_sol
mode_diag(icount)%component = cp_no3

! nitrate to accum-sol
icount = icount + 1
mode_diag(icount)%item = 580
mode_diag(icount)%mode = mode_acc_sol
mode_diag(icount)%component = cp_no3

! nitrate to coarse-sol
icount = icount + 1
mode_diag(icount)%item = 581
mode_diag(icount)%mode = mode_cor_sol
mode_diag(icount)%component = cp_no3

! sodium nitrate to accum-sol
icount = icount + 1
mode_diag(icount)%item = 582
mode_diag(icount)%mode = mode_acc_sol
mode_diag(icount)%component = cp_nn

! sodium nitrate to coarse-sol
icount = icount + 1
mode_diag(icount)%item = 583
mode_diag(icount)%mode = mode_cor_sol
mode_diag(icount)%component = cp_nn

! check that number of diagnostics matches the length of the array
IF (icount /= num_mode_diags) THEN
  cmessage = "Number of mode diagnostics wrong: icount /= num_mode_diags"
  CALL ereport ('UKCA_EMISS_DIAGS_MODE_INIT', icount, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_emiss_diags_mode_init


SUBROUTINE ukca_emiss_diags_mode (row_length, rows, model_levels, area,        &
    len_stashwork, stashwork)

USE ukca_config_specification_mod, ONLY: glomap_variables

USE ukca_mode_setup,    ONLY: moment_mass

USE ukca_emiss_mod,     ONLY: emissions, num_em_flds, marine_oc_online
USE parkind1,           ONLY: jprb, jpim

USE yomhook,            ONLY: lhook, dr_hook
USE copydiag_mod,       ONLY: copydiag

USE ukca_um_legacy_mod, ONLY: len_stlist, stindex, stlist, num_stash_levels,   &
                              stash_levels, si, sf, si_last,                   &
                              stashcode_glomap_sec, copydiag_3d

IMPLICIT NONE


! Subroutine arguments
INTEGER, INTENT(IN)    :: row_length        ! Model dimensions
INTEGER, INTENT(IN)    :: rows
INTEGER, INTENT(IN)    :: model_levels
REAL, INTENT(IN) :: area(row_length, rows, model_levels)  ! Area of grid cell

INTEGER, INTENT(IN)    :: len_stashwork     ! Length of diagnostics array
REAL,    INTENT(IN OUT) :: stashwork (len_stashwork) ! Diagnostics array

! Local variables

INTEGER    :: k, l         ! counters / indices
INTEGER    :: section      ! stash section
INTEGER    :: item         ! stash item
INTEGER    :: im_index     ! internal model index
INTEGER    :: imode
INTEGER    :: icp
INTEGER    :: ilev
INTEGER    :: err_code

REAL :: em_diags (row_length, rows, model_levels)

LOGICAL, SAVE :: lfirst = .TRUE.  ! Indicator of first call to this routine
LOGICAL       :: lfound = .TRUE.  ! Has corresponding emission been found

CHARACTER(LEN=errormessagelength) :: cmessage  ! Error return message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=30) :: varname                   ! Name of emission
CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_DIAGS_MODE'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (lfirst) THEN
  CALL ukca_emiss_diags_mode_init
  lfirst = .FALSE.
END IF

im_index = 1
section = stashcode_glomap_sec

! ------------------------------------------------------------------------
! For each diagnostic add contributions from each source emission, then
! copy to corresponding stash item.
DO k = 1, num_mode_diags
  imode = mode_diag(k)%mode
  icp = mode_diag(k)%component
  item = mode_diag(k)%item
  varname = mode_diag(k)%vname

  IF (sf(item, section)) THEN
    ! Loop through all emissions to store the diagnostics for this mode/cpt
    em_diags(:,:,:) = 0.0
    lfound = .FALSE.
    DO l = 1, num_em_flds
      IF (varname == marine_oc_online) THEN
        ! Match emissions to name for online marine OC to separate
        IF ((emissions(l)%tracer_name == 'mode_emiss') .AND.                   &
            (emissions(l)%moment == moment_mass)       .AND.                   &
            (emissions(l)%mode == imode)               .AND.                   &
            (emissions(l)%var_name == 'pmoc_online_emission') .AND.            &
            (emissions(l)%component == icp)) THEN

          lfound = .TRUE.

          IF (emissions(l)%three_dim) THEN
            em_diags(:,:,:) = em_diags(:,:,:) + emissions(l)%diags(:,:,:)
          ELSE
            em_diags(:,:,1) = em_diags(:,:,1) + emissions(l)%diags(:,:,1)
            ! For 2D emissions, diags stored in level 1, vert_scaling unused
          END IF
        END IF  ! if this emission matches the diagnostic
      ELSE
        ! Not marine OC emissions
        IF ((emissions(l)%tracer_name == 'mode_emiss') .AND.                   &
            (emissions(l)%moment == moment_mass)       .AND.                   &
            (emissions(l)%mode == imode)               .AND.                   &
            (emissions(l)%component == icp)) THEN

          lfound = .TRUE.
          IF (emissions(l)%three_dim) THEN
            em_diags(:,:,:) = em_diags(:,:,:) + emissions(l)%diags(:,:,:)
          ELSE
            ! For 2D emissions, diags stored in level 1 but vertical scaling
            ! stored and we use it to re-project onto levels.  Column total of
            ! vert_scaling_3d is always equal to 1.
            DO ilev = 1, model_levels
              em_diags(:,:,ilev) = em_diags(:,:,ilev) +                        &
                             emissions(l)%diags(:,:,1) *                       &
                             emissions(l)%vert_scaling_3d(:,:,ilev)
            END DO
          END IF  ! three_dim
        END IF  ! if this emission matches the diagnostic

      END IF    ! varname

    END DO  ! loop over emissions

    ! Raise a warning (and skip to next diag) if no emissions for this diag
    IF (.NOT. lfound) THEN
      cmessage = 'No mode emissions corresponding to section 38 ' //           &
                 'diagnostic request'
      err_code = -item
      CALL ereport ('UKCA_EMISS_DIAGS_MODE', err_code, cmessage)
      CYCLE
    END IF

    ! Convert from kg/m2/s to mol/gridbox/s
    em_diags(:,:,:) = em_diags(:,:,:) * area(:,:,:) / glomap_variables%mm(icp)

    IF (varname == marine_oc_online) THEN
      ! Use first layer of em_diags only
      CALL copydiag(stashwork (si(item,section,im_index):                      &
                 si_last(item,section,im_index)),                              &
                 em_diags(:,:,1), row_length, rows)
    ELSE
      ! everything else

      CALL copydiag_3d (stashwork (si(item,section,im_index):                  &
           si_last(item,section,im_index)),                                    &
               em_diags (:,:,:),                                               &
               row_length, rows,model_levels,                                  &
               stlist(:,stindex(1,item,section,im_index)), len_stlist,         &
               stash_levels, num_stash_levels+1)
    END IF  ! varname
  END IF  ! if diagnostic requested
END DO  ! loop over diags
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_emiss_diags_mode

END MODULE ukca_emiss_diags_mode_mod
