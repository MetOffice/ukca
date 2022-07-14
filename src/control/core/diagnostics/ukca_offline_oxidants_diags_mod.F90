! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to contain code output the offline oxidant fields as diagnostics.
!
! Method:
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and
!  The Met Office. See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! Contained subroutines in this module:
!   ukca_offline_oxidants_diags        (Public)
!
! --------------------------------------------------------------------------

MODULE ukca_offline_oxidants_diags_mod

USE parkind1,         ONLY: jpim, jprb      ! DrHook
USE yomhook,          ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Default private
PRIVATE

PUBLIC :: ukca_offline_oxidants_diags

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName='UKCA_OFFLINE_OXIDANTS_DIAGS_MOD'

CONTAINS

SUBROUTINE ukca_offline_oxidants_diags(row_length, rows, model_levels,         &
  totnodens, o3_offline_diag, oh_offline_diag, no3_offline_diag,               &
  ho2_offline_diag, h2o2_offline_diag,                                         &
  len_stashwork,STASHwork50)
! ---------------------------------------------------------------------
! Description:
!  To place offline oxidants fields into the section 50 STASHwork array
! Method:
!  For the oxidants which are constant fields, the diagnostic fields are
!  allocated and filled in routine UKCA_SET_DIURNAL_OX.  For the H2O2
!  limiting field, there is no requirement for a seperate diagnostic field
!  as it is not given a diurnal cycle.
! ---------------------------------------------------------------------

USE ukca_constants,      ONLY: c_o3, c_oh, c_no3, c_ho2
USE ukca_um_legacy_mod,  ONLY: len_stlist, stindex, stlist, num_stash_levels,  &
                               stash_levels, si, sf, si_last, copydiag_3d,     &
                               UKCA_diag_sect


IMPLICIT NONE

! Input arguments
INTEGER, INTENT(IN) :: row_length         ! No. of rows
INTEGER, INTENT(IN) :: rows               ! No. of columns
INTEGER, INTENT(IN) :: model_levels       ! No. of levels
INTEGER, INTENT(IN) :: len_stashwork      ! Length of stashwork array
! Density in molecules/m^3
REAL, INTENT(IN)    :: totnodens(row_length,rows,model_levels)
! Diagnostic fields (mass mixing ratio)
REAL, INTENT(IN) :: o3_offline_diag(row_length*rows,model_levels)   ! O3 diag
REAL, INTENT(IN) :: oh_offline_diag(row_length*rows,model_levels)   ! OH diag
REAL, INTENT(IN) :: no3_offline_diag(row_length*rows,model_levels)  ! NO3 diag
REAL, INTENT(IN) :: ho2_offline_diag(row_length*rows,model_levels)  ! HO2 diag
REAL, INTENT(IN) :: h2o2_offline_diag(row_length,rows,model_levels) ! H2O2 diag

REAL, INTENT(IN OUT) :: STASHwork50(len_stashwork)  ! Work array

! Local variables

INTEGER :: section                        ! Stash section
INTEGER :: item                           ! Stash item
INTEGER :: im_index                       ! internal model index
INTEGER :: k                              ! loop index

REAL :: field3d(row_length,rows,model_levels)  ! To contain 3d diagnostic

REAL (KIND=jprb)   :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_OFFLINE_OXIDANTS_DIAGS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

im_index = 1

! Offline oxidant diagnostics

section = ukca_diag_sect
DO item = 206,210
  IF (sf(item,section) .AND. item == 206) THEN
    DO k=1,model_levels
      field3d(:,:,k) = RESHAPE(o3_offline_diag(:,k),[row_length,rows])*        &
        c_o3/(totnodens(:,:,k)/1.0e6)
    END DO

    CALL copydiag_3d (stashwork50(si(item,section,im_index):                   &
         si_last(item,section,im_index)),                                      &
      field3d,                                                                 &
      row_length,rows,model_levels,                                            &
      stlist(:,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1)
  ELSE IF (sf(item,section) .AND. item == 207) THEN
    DO k=1,model_levels
      field3d(:,:,k) = RESHAPE(oh_offline_diag(:,k),[row_length,rows])*        &
        c_oh/(totnodens(:,:,k)/1.0e6)
    END DO

    CALL copydiag_3d (stashwork50(si(item,section,im_index):                   &
         si_last(item,section,im_index)),                                      &
      field3d,                                                                 &
      row_length,rows,model_levels,                                            &
      stlist(:,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1)
  ELSE IF (sf(item,section) .AND. item == 208) THEN
    DO k=1,model_levels
      field3d(:,:,k) = RESHAPE(no3_offline_diag(:,k),[row_length,rows])*       &
        c_no3/(totnodens(:,:,k)/1.0e6)
    END DO

    CALL copydiag_3d (stashwork50(si(item,section,im_index):                   &
         si_last(item,section,im_index)),                                      &
      field3d,                                                                 &
      row_length,rows,model_levels,                                            &
      stlist(:,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1)
  ELSE IF (sf(item,section) .AND. item == 209) THEN
    DO k=1,model_levels
      field3d(:,:,k) = RESHAPE(ho2_offline_diag(:,k),[row_length,rows])*       &
        c_ho2/(totnodens(:,:,k)/1.0e6)
    END DO

    CALL copydiag_3d (stashwork50(si(item,section,im_index):                   &
         si_last(item,section,im_index)),                                      &
      field3d,                                                                 &
      row_length,rows,model_levels,                                            &
      stlist(:,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1)
  ELSE IF (sf(item,section) .AND. item == 210) THEN

    CALL copydiag_3d (stashwork50(si(item,section,im_index):                   &
         si_last(item,section,im_index)),                                      &
      h2o2_offline_diag(:,:,:),                                                &
      row_length,rows,model_levels,                                            &
      stlist(:,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1)
  END IF

END DO       ! item

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE ukca_offline_oxidants_diags

END MODULE ukca_offline_oxidants_diags_mod
