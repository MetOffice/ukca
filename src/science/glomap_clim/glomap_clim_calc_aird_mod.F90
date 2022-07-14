! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!   Calculate the density of dry air using the ideal gas law.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE glomap_clim_calc_aird_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='GLOMAP_CLIM_CALC_AIRD_MOD'

CONTAINS

SUBROUTINE glomap_clim_calc_aird ( n_points, p_theta_levels_1d,                &
                                   t_theta_levels_1d, aird )

USE chemistry_constants_mod,         ONLY:                                     &
    boltzmann

USE parkind1,                        ONLY:                                     &
    jprb,                                                                      &
    jpim

USE yomhook,                         ONLY:                                     &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: n_points

!                                   Pressure on theta levels (in 1D)
REAL(KIND=real_umphys),INTENT(IN) :: p_theta_levels_1d( n_points )

!                                   Temperature on theta levels (in 1D)
REAL(KIND=real_umphys),INTENT(IN) :: t_theta_levels_1d( n_points )

REAL(KIND=real_umphys),INTENT(OUT):: aird(n_points)

! Local variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName='GLOMAP_CLIM_CALC_AIRD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! The density of dry air is calculated using the ideal gas law
aird(:) = p_theta_levels_1d(:) /                                               &
          ( t_theta_levels_1d(:) * ( boltzmann * 1000000.0 ) )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE glomap_clim_calc_aird

END MODULE glomap_clim_calc_aird_mod
