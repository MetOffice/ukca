! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Purpose: Subroutine to to assign wet deposition rates in s-1 to
!          array dpw in module asad_mod. This is passed for use in
!          ASAD chemistry
!          Based on routine from Cambridge TOMCAT model
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD routine CDRIVE.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_wetdep_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_WETDEP_MOD'

CONTAINS

SUBROUTINE ukca_wetdep(wetrt,n_points)

USE asad_mod,         ONLY: ndepw, nldepw, dpw, jpdw
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook


IMPLICIT NONE


INTEGER, INTENT(IN) :: n_points                  ! No of spatial

REAL, INTENT(IN) :: wetrt(n_points,jpdw) ! Wet dep rates

!       Local variables

INTEGER :: i                         ! Loop variable
INTEGER :: js                        ! Loop variable
INTEGER :: nspec                     ! Pointer for species

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_WETDEP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO js = 1,ndepw
  nspec = nldepw(js)
  DO i = 1,n_points
    dpw(i,nspec)= wetrt(i,js)
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_wetdep

END MODULE ukca_wetdep_mod
