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
! Purpose: Routine to read in dry deposition velocities at 1 metre.
!          Original version taken from the Cambridge TOMCAT model.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD routine ASAD_CINIT.
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
MODULE ukca_inddep_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_INDDEP_MOD'

CONTAINS

SUBROUTINE ukca_inddep

USE asad_mod,           ONLY: depvel, jddepc, jddept, jpdd
USE ukca_chem_defs_mod, ONLY: depvel_defs
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport


USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER :: errcode                ! Variable passed to ereport
INTEGER :: ns                   ! Loop variable
INTEGER :: nt                   ! Loop variable
INTEGER :: nc                   ! Loop variable

CHARACTER(LEN=errormessagelength)  :: cmessage  ! String for error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INDDEP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!       Reading dry deposition velocities from module
!       Dataset is for 6 time points and 5 land categories.

IF (SIZE(depvel_defs) /= SIZE(depvel)) THEN
  cmessage='sizes of depvel_defs and depvel are inconsistent'
  errcode=1
  CALL ereport('UKCA_INDDEP',errcode,cmessage)
END IF

DO ns=1,jpdd
  DO nc=1,jddepc
    DO nt=1,jddept
      depvel(nt,nc,ns)=depvel_defs(nt,nc,ns)
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_inddep
END MODULE ukca_inddep_mod
