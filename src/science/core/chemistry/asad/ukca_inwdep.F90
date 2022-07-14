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
! Purpose: Subroutine to read in coefficients for calculating
!          effective Henry's Law coefficients. Original version
!          taken from the Cambridge TOMCAT model.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD routine ASAD_cinit.
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
MODULE ukca_inwdep_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_INWDEP_MOD'

CONTAINS

SUBROUTINE ukca_inwdep

USE asad_mod,             ONLY: kd298, k298, ddhr, dhr,                        &
                                ct_k298, ct_kd298, ct_ddhr,                    &
                                ct_dhr, jpdw
USE ukca_chem_defs_mod,   ONLY: henry_defs
USE ukca_chem_offline,    ONLY: henry_defs_const, nwet_constant
USE ukca_config_specification_mod, ONLY: ukca_config
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER :: errcode                ! Variable passed to ereport
INTEGER            :: ns       ! Loop variable

CHARACTER (LEN=errormessagelength) :: cmessage ! String for error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INWDEP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!       Use module to define Henry constants

IF (SIZE(henry_defs) /= jpdw*6) THEN
  cmessage='jpdw and henry_defs are inconsistent'
  errcode=1
  CALL ereport('UKCA_INWDEP',errcode,cmessage)
END IF

DO ns=1,jpdw
  k298(ns)    = henry_defs(1,ns)
  dhr(ns)     = henry_defs(2,ns)
  kd298(ns,1) = henry_defs(3,ns)
  ddhr(ns,1)  = henry_defs(4,ns)
  kd298(ns,2) = henry_defs(5,ns)
  ddhr(ns,2)  = henry_defs(6,ns)
END DO

! Dissociation of offline ozone
IF (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) THEN
  IF (nwet_constant > 0) THEN
    DO ns=1,nwet_constant
      ct_k298(ns)    = henry_defs_const(1,ns)
      ct_dhr(ns)     = henry_defs_const(2,ns)
      ct_kd298(ns,1) = henry_defs_const(3,ns)
      ct_ddhr(ns,1)  = henry_defs_const(4,ns)
      ct_kd298(ns,2) = henry_defs_const(5,ns)
      ct_ddhr(ns,2)  = henry_defs_const(6,ns)
    END DO
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_inwdep
END MODULE ukca_inwdep_mod
