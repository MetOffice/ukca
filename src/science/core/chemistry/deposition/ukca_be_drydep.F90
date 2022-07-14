! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  To assign dry deposition rates in s-1 to array odryrt
!
!           Called from UKCA_CHEMISTRY_CTL
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_be_drydep_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_BE_DRYDEP_MOD'

CONTAINS

SUBROUTINE ukca_be_drydep(k, n_pnts, nlev_with_ddep2, dryrt, odryrt)

USE asad_mod,          ONLY:  ndepd, nldepd, jpspec, jpdd
USE ukca_config_specification_mod, ONLY: ukca_config
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: k                  ! Level number
INTEGER, INTENT(IN) :: n_pnts             ! No of spatial points

! Number of levels in dry depositon applied
INTEGER, INTENT(IN) :: nlev_with_ddep2(n_pnts)

REAL, INTENT(IN) :: dryrt(n_pnts,jpdd)    ! Dry dep rates in s-1
REAL, INTENT(OUT):: odryrt(n_pnts,jpspec) ! Dry dep rates in s-1

!     Local variables

INTEGER, PARAMETER :: bottom_level = 1

INTEGER :: i                              ! Loop variable
INTEGER :: js                             ! Loop variable
INTEGER :: nspec                          ! Pointer for species

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_BE_DRYDEP'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO js = 1, jpspec
  DO i = 1, n_pnts
    odryrt(i,js) = 0.0
  END DO
END DO

IF (ukca_config%l_ukca_intdd) THEN

  DO i = 1,n_pnts
    IF (k <= nlev_with_ddep2(i)) THEN
      DO js = 1,ndepd
        nspec = nldepd(js)
        odryrt(i,nspec) = dryrt(i,js)
      END DO
    END IF
  END DO

ELSE

  IF (k == bottom_level) THEN

    DO js = 1,ndepd
      nspec = nldepd(js)
      DO i = 1,n_pnts
        odryrt(i,nspec) = dryrt(i,js)
      END DO
    END DO

  END IF

END IF   ! End of IF (l_ukca_intdd) statement

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_be_drydep
END MODULE ukca_be_drydep_mod
