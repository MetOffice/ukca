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
! Purpose: This is a dummy version of a user-supplied subroutine.
!     The role of this routine is to perform any post-processing
!     or house keeping required by the user's heterogeneous
!     chemistry scheme.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE
!
!     Example
!     -------
!     As an example of what this routine might be used for,
!     in our stratospheric chemistry scheme, we have a
!     heterogeneous chemistry scheme that takes into account
!     the amount of water and nitic acid in the solid phase.
!     Therefore, the routine hetero in our chemistry partitions
!     the gaseous HNO3 into a solid and gaseous form. The amount
!     in the solid phase is stored in a module. This routine
!     is responsible for adding that solid phase amount of
!     HNO3 back into the gaseous phase at the end of the timestep.
!     ie. routine hetero computes the amount of solid HNO3, subtracts
!     that from the array y(), put it in common, where this routine
!     adds it back to y() before ASAD is exited.
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
MODULE asad_posthet_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_POSTHET_MOD'

CONTAINS

SUBROUTINE asad_posthet


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER                       ::  icode
CHARACTER (LEN=errormessagelength)            ::  cmessage
CHARACTER (LEN=* ), PARAMETER ::  RoutineName='ASAD_POSTHET'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
cmessage = 'Routine should not be callable'
icode = 1

CALL ereport(RoutineName,icode,cmessage)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_posthet
END MODULE asad_posthet_mod
