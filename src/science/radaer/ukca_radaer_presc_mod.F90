! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module structure and procedures for UKCA RADAER prescriptions
!
! Method:
!   Declares the UKCA RADAER prescription data structure where the prescribed
!   distributions and other relevant data are stored. UKCA RADAER prescriptions
!   distributions are read from netCDF files.
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
MODULE ukca_radaer_presc_mod

USE filenamelength_mod, ONLY:                                                  &
    filenamelength

IMPLICIT NONE

! All variables and functions private by default
PRIVATE
PUBLIC :: set_presc_ssa_field

! UKCA RADAER prescription data structure
TYPE, PUBLIC :: ukca_radaer_presc_struct
  CHARACTER (LEN=filenamelength)  :: file_name    ! Name of source file
  CHARACTER (LEN=80)   :: var_name     ! Name of variable in file
  INTEGER              :: varid        ! ID of variable in file
  CHARACTER (LEN=256)  :: std_name     ! standard_name attrib in NetCDF files
  CHARACTER (LEN=256)  :: long_name    ! long_name     attrib in NetCDF files
  CHARACTER (LEN=30)   :: units        ! Units of the field

  INTEGER              :: update_freq  ! Update frequency (hours)
  INTEGER              :: update_type  ! 1 serial, 2 periodic, ...
  INTEGER              :: last_update  ! Num anc update intervals
                                       ! at last update
  LOGICAL              :: l_update     ! True if field updated in a given tstep

  INTEGER              :: ndims        ! Number of dimensions (excluding time)
                                       ! (default: 3)

  INTEGER              :: n_specbands  ! Additional dimension for optical data
                                       ! that is defined on spectral bands
                                       ! (default: 0 meaning none)

  REAL, ALLOCATABLE:: values_3d(:,:,:)   ! UKCA RADAER 3D data
  REAL, ALLOCATABLE:: values_4d(:,:,:,:) ! UKCA RADAER 4D data

END TYPE ukca_radaer_presc_struct

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_RADAER_PRESC_MOD'

CONTAINS

SUBROUTINE set_presc_ssa_field(                                                &
      ! Fixed array dimensions
      npd_profile, npd_layer, npd_band,                                        &
      ! Actual array dimensions
      n_profile, n_layer, n_band,                                              &
      ! Full grid indices
      col_list, row_list,                                                      &
      ! Input and output prescribed extinction properties
      ukca_radaer_presc, ukca_radaer_presc_ssa)
  !
  ! Modules used
  !

USE def_ukca_radaer_presc, ONLY: t_ukca_radaer_presc

USE rad_input_mod,   ONLY: l_extra_top

USE um_types, ONLY: real_umphys

USE yomhook,         ONLY: lhook, dr_hook
USE parkind1,        ONLY: jprb, jpim

IMPLICIT NONE

!
! Fixed array dimensions
!
INTEGER, INTENT(IN) :: npd_profile,                                            &
                       npd_layer,                                              &
                       npd_band

!
! Actual array dimensions
!
INTEGER, INTENT(IN) :: n_profile,                                              &
                       n_layer,                                                &
                       n_band


! List of column and row indices on the full grid
INTEGER, INTENT(IN) :: col_list(npd_profile)
INTEGER, INTENT(IN) :: row_list(npd_profile)

! EasyAerosol optical properties
TYPE (t_ukca_radaer_presc), INTENT(IN) :: ukca_radaer_presc

! Output SSA field
REAL(KIND=real_umphys), INTENT(OUT) ::                                         &
    ukca_radaer_presc_ssa(npd_profile, npd_layer, npd_band)

!
! Local variables
!
INTEGER :: i, j, k
INTEGER :: i_top_copy
INTEGER :: this_layr
REAL(KIND=real_umphys) :: this_ssa

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SET_PRESC_SSA_FIELD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                          &
                        zhook_in, zhook_handle)

IF (l_extra_top) THEN
  i_top_copy = 2
ELSE
  i_top_copy = 1
END IF

! Calculate SSA using the Easy Aerosol extinctions
DO k = 1, n_band

  DO j = i_top_copy, n_layer

    this_layr = n_layer+1-j

    DO i = 1, n_profile

      this_ssa = 1.0 - (ukca_radaer_presc%absorption(col_list(i),              &
                                                     row_list(i),              &
                                                     this_layr,                &
                                                     k) /                      &
                        ukca_radaer_presc%extinction(col_list(i),              &
                                                     row_list(i),              &
                                                     this_layr,                &
                                                     k))

      ukca_radaer_presc_ssa(i, j, k) = MIN(1.0, MAX(0.0, this_ssa))

    END DO ! i

  END DO ! j

END DO ! k

!
! If using an extra top layer, set its prescribed SSA properties
! to one.
!
IF (l_extra_top) THEN

  DO k = 1, n_band

    DO i = 1, n_profile

      ukca_radaer_presc_ssa(i, 1, k) = 1.0

    END DO ! i

  END DO ! k

END IF ! l_extra_top

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                          &
                        zhook_out, zhook_handle)

RETURN
END SUBROUTINE set_presc_ssa_field

END MODULE ukca_radaer_presc_mod
