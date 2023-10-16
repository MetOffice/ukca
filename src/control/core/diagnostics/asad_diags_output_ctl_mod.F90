! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module containing top-level subroutine for output of ASAD diagnostics
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds,
! University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE asad_diags_output_ctl_mod

IMPLICIT NONE

PRIVATE

PUBLIC asad_diags_output_ctl

CHARACTER(LEN=*), PARAMETER :: ModuleName='ASAD_DIAGS_OUTPUT_CTL_MOD'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE asad_diags_output_ctl(error_code_ptr,                               &
                                 row_length, rows, model_levels,               &
                                 stashwork_size, diagnostics, stashwork,       &
                                 error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
!   Output ASAD flux diagnostics by summing the contributions to each of
!   the requested fields and copying the results to arrays supplied by
!   the parent application. If the parent is the Unified Model, the
!   diagnostic requests may include UM legacy requests that are output
!   by copying to the relevant STASH work array instead.
! ----------------------------------------------------------------------

USE ukca_um_legacy_mod, ONLY: copydiag, copydiag_3d,                           &
                              len_stlist, stindex, stlist,                     &
                              num_stash_levels, stash_levels, si, sf, si_last
USE asad_chem_flux_diags, ONLY: stash_handling, asad_chemdiags, n_asad_diags,  &
                                n_chemdiags
USE ukca_diagnostics_type_mod, ONLY: diagnostics_type
USE ukca_diagnostics_output_mod, ONLY: update_diagnostics_2d_real,             &
                                       update_diagnostics_3d_real
USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine arguments

! Return code for status reporting
INTEGER, POINTER, INTENT(IN) :: error_code_ptr

! Array sizes
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: stashwork_size

! Diagnostic request info and pointers to parent arrays for diagnostic output
TYPE(diagnostics_type), INTENT(IN OUT) :: diagnostics

! STASH array for UM diagnostic output (legacy requests only)
REAL, INTENT(IN OUT) :: stashwork(1:stashwork_size)

! Optional arguments for status reporting
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables

! Flux arrays
REAL, ALLOCATABLE :: upload_array_3D(:,:,:)
REAL, ALLOCATABLE :: upload_array_2D(:,:)

INTEGER       :: item                             ! stash item
INTEGER       :: section                          ! stash section
INTEGER       :: im_index                         ! atmos model index
INTEGER       :: i,j,k,l                          ! Counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_DIAGS_OUTPUT_CTL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

im_index = 1

ALLOCATE(upload_array_2D(1:row_length,1:rows))
ALLOCATE(upload_array_3D(1:row_length,1:rows,1:model_levels))

DO l=1,n_asad_diags

  item = stash_handling(l)%stash_item
  section = stash_handling(l)%stash_section
  upload_array_2D(:,:)   = 0.0
  upload_array_3D(:,:,:) = 0.0

  IF (stash_handling(l)%len_dim3 == 1) THEN                   ! 2D field

    DO k=1,stash_handling(l)%number_of_fields
         ! will only be summing surface fields in this case
      upload_array_2D(:,:) = upload_array_2D(:,:) +                            &
              asad_chemdiags(                                                  &
              stash_handling(l)%chemdiags_location(k)                          &
              )%throughput(:,:,1)
    END DO

    ! Output diagnostic field
    IF (stash_handling(l)%l_um_legacy_request) THEN
      ! UM legacy request: copy to STASHwork if required on this time step
      IF (sf(item,section)) THEN
        CALL copydiag(stashwork(si(item,section,im_index):                     &
                      si_last(item,section,im_index)),                         &
                      upload_array_2D(:,:),row_length,rows)
      END IF
    ELSE
      ! Copy to diagnostic array supplied by parent model if required on this
      ! time step (i.e. not deactivated by parent or skipped by UKCA on a
      ! non-chemistry time step)
      CALL update_diagnostics_2d_real(error_code_ptr,                          &
                                      stash_handling(l)%diagname,              &
                                      upload_array_2D, diagnostics,            &
                                      error_message=error_message,             &
                                      error_routine=error_routine)
      IF (error_code_ptr > 0) THEN
        IF (lhook)                                                             &
          CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
        RETURN
      END IF
    END IF

  ELSE ! len_dim3 /= 1

    DO k=1,stash_handling(l)%number_of_fields
        ! may be summing both 3D and surface fields into the same 3D array
      DO j=1,asad_chemdiags(                                                   &
        stash_handling(l)%chemdiags_location(k))%num_levs
        upload_array_3D(:,:,j) = upload_array_3D(:,:,j) +                      &
                asad_chemdiags(                                                &
                stash_handling(l)%chemdiags_location(k)                        &
                )%throughput(:,:,j)
      END DO
    END DO

    ! Output diagnostic field
    IF (stash_handling(l)%l_um_legacy_request) THEN
      ! UM legacy request: copy to STASHwork if required on this time step
      IF (sf(item,section)) THEN
        CALL copydiag_3d(stashwork(si(item,section,im_index):                  &
                         si_last(item,section,im_index)),                      &
                         upload_array_3D(:,:,:),                               &
                         row_length,rows,model_levels,                         &
                         stlist(:,stindex(1,item,section,im_index)),           &
                         len_stlist,stash_levels,num_stash_levels+1)
      END IF
    ELSE
      ! Copy to diagnostic array(s) supplied by parent model if required on
      ! this time step (i.e. to satisfy request(s) not deactivated by parent or
      ! skipped by UKCA on a non-chemistry time step)
      CALL update_diagnostics_3d_real(error_code_ptr,                          &
                                      stash_handling(l)%diagname,              &
                                      upload_array_3D, diagnostics,            &
                                      error_message=error_message,             &
                                      error_routine=error_routine)
      IF (error_code_ptr > 0) THEN
        IF (lhook)                                                             &
          CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
        RETURN
      END IF
    END IF

  END IF ! len_dim3

END DO       ! l,n_asad_diags

! Deallocate arrays, if are able to
DO i=1,n_chemdiags
  IF (asad_chemdiags(i)%can_deallocate)                                        &
       DEALLOCATE(asad_chemdiags(i)%throughput)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_diags_output_ctl

END MODULE asad_diags_output_ctl_mod
