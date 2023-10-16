! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module containing type definitions and parameters used in UKCA's
!   diagnostic handling system.
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

MODULE ukca_diagnostics_type_mod

USE ukca_fieldname_mod, ONLY: maxlen_diagname

IMPLICIT NONE

PUBLIC

! --- Parameter definitions ---

! Field group codes for collating diagnostic fields in arrays
! 'flat' group comprises fields defined on a flat spatial grid with 2D
! representation internally (output can have reduced dimension)
! 'fullht' group comprises fields defined on a full height spatial grid with 3D
! representation internally (output can have reduced dimension)
INTEGER, PARAMETER :: dgroup_flat_real = 1     ! 2D spatial real
INTEGER, PARAMETER :: dgroup_fullht_real = 2   ! 3D spatial real

! Number of diagnostic groups
INTEGER, PARAMETER :: n_diag_group = 2

! --- Status flag values for requested diagnostics ---
! The parent can set these to indicate that individual requests should be
! treated as active or inactive.
! Requests are accepted as such if the diagnostic is available in the current
! UKCA configuration. Otherwise they are marked unavailable.
! The status of each active request is updated when serviced or when
! skipped due to non-availability on the current time step.
INTEGER, PARAMETER :: diag_status_inactive = 0
                                ! Request set/accepted as inactive
INTEGER, PARAMETER :: diag_status_requested = 1
                                ! Request set/accepted as active
INTEGER, PARAMETER :: diag_status_valid = 2
                                ! Data provided
INTEGER, PARAMETER :: diag_status_skipped = 3
                                ! Data unavailable on current time step
INTEGER, PARAMETER :: diag_status_unavailable = 4
                                ! Diagnostic unavailable in UKCA configuration
                                ! (request denied)

! --- Type definitions ---

! Type to hold the array bounds associated with a diagnostic
TYPE :: bounds_type
  INTEGER :: dim1_lower  ! Lower bound of dimension 1
  INTEGER :: dim1_upper  ! Upper bound of dimension 1
  INTEGER :: dim2_lower  ! Lower bound of dimension 2
  INTEGER :: dim2_upper  ! Upper bound of dimension 2
  INTEGER :: dim3_lower  ! Lower bound of dimension 3
  INTEGER :: dim3_upper  ! Upper bound of dimension 3
END TYPE bounds_type

! Type to hold an entry in the master diagnostics list
TYPE :: diag_entry_type
  CHARACTER(LEN=maxlen_diagname) :: varname  ! Diagnostic name
  INTEGER :: group            ! Diagnostic group that field belongs to
  INTEGER :: group_alt        ! Alternative diagnostic group that subfield
                              ! belongs to
  TYPE(bounds_type) :: bound  ! Array bounds of diagnostic field
  LOGICAL :: l_available      ! Availability in current UKCA configuration
  LOGICAL :: l_chem_timestep  ! True if only output on chemistry time steps
  INTEGER :: asad_id          ! ASAD diagnostic identifier (currently
                              ! assumed to be a UM Section-50 STASH code)
END TYPE diag_entry_type

! Type to hold the set of diagnostic requests for a single group
TYPE :: diag_requests_type
  CHARACTER(LEN=maxlen_diagname), ALLOCATABLE :: varnames(:)
                                           ! Names of requested diagnostics
  INTEGER, ALLOCATABLE :: status_flags(:)  ! Request status flags
  INTEGER, ALLOCATABLE :: i_master(:)      ! Indices in master diagnostics list
END TYPE diag_requests_type

! Type to hold diagnostic status flags for a single group
TYPE :: diag_status_type
  INTEGER, ALLOCATABLE :: status_flags(:)  ! Diagnostic status flags
END TYPE diag_status_type

! Type to hold data for servicing diagnostic requests
TYPE :: diagnostics_type
  ! Pointer to access diagnostic requests in each group
  ! Note: flat real group requests for full height group diagnostics
  ! are allowed (level 1 subfield is supplied)
  TYPE(diag_requests_type), POINTER :: requests_ptr(:)
  ! Number of requests for each group
  INTEGER :: n_request(n_diag_group)
  ! Dimension of output array provided by parent for each group
  INTEGER :: dimension_out(n_diag_group)
  ! Pointers to access output arrays provided by parent
  REAL, POINTER :: value_0d_real_ptr(:)
  REAL, POINTER :: value_1d_real_ptr(:,:)
  REAL, POINTER :: value_2d_real_ptr(:,:,:)
  REAL, POINTER :: value_3d_real_ptr(:,:,:,:)
  ! Status flags for diagnostic output
  TYPE(diag_status_type) :: outvalue_status(n_diag_group)
END TYPE diagnostics_type

END MODULE ukca_diagnostics_type_mod
