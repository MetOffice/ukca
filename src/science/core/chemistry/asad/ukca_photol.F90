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
! Purpose: Subroutine to to assign photolysis rates in s-1 to
!          array rk in module asad_mod. This is passed for use in
!          asad chemistry.
!          Based on photol.F from Cambridge TOMCAT model
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
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_photol_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_PHOTOL_MOD'

CONTAINS

SUBROUTINE ukca_photol(prt,n_points)

USE asad_mod,        ONLY:  rk, nprkx, jppj
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) :: n_points         ! No of spatial points

REAL, INTENT(IN) :: prt(n_points,jppj)  ! Photolysis rates

! Local variables

INTEGER :: jl                           ! Loop variable
INTEGER :: jr                           ! Loop variable

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PHOTOL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO jl = 1, n_points
  DO jr = 1, jppj
    rk(jl,nprkx(jr)) = prt(jl,jr)
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_photol
END MODULE ukca_photol_mod
