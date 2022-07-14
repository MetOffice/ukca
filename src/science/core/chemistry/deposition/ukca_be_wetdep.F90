! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Subroutine to assign wet deposition rates in s-1 to
!  array owetrt.
!
!  Called from UKCA routine UKCA_CHEMISTRY_CTL.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_be_wetdep_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_BE_WETDEP_MOD'

CONTAINS

SUBROUTINE ukca_be_wetdep(n_pnts, wetrt, owetrt)

USE asad_mod,      ONLY: ndepw, nldepw, jpspec, jpdw
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_pnts             ! No of spatial pts
REAL, INTENT(IN) :: wetrt(n_pnts,jpdw)    ! Wet dep rates (s-1)
REAL, INTENT(OUT):: owetrt(n_pnts,jpspec) ! Wet dep rates (s-1)

!       Local variables

INTEGER :: i                         ! Loop variable
INTEGER :: js                        ! Loop variable
INTEGER :: nspec                     ! Pointer for species

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_BE_WETDEP'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO js = 1, jpspec
  DO i = 1, n_pnts
    owetrt(i,js) = 0.0
  END DO
END DO

DO js = 1,ndepw
  nspec = nldepw(js)
  DO i = 1,n_pnts
    owetrt(i,nspec)= wetrt(i,js)
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_be_wetdep
END MODULE ukca_be_wetdep_mod
