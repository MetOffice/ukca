! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Dummy module to replace the DrHook library.
!
! Part of the UKCA model, a community model supported by The Met Office
! and NCAS, with components initially provided by The University of
! Cambridge, University of Leeds and The Met Office. See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE yomhook

USE parkind1, ONLY: jpim, jprb

IMPLICIT NONE

LOGICAL, PARAMETER :: lhook = .FALSE.

CONTAINS

SUBROUTINE dr_hook(routine_name, code, handle)

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*),   INTENT(IN)    :: routine_name
INTEGER(KIND=jpim), INTENT(IN)    :: code
REAL(KIND=jprb),    INTENT(IN OUT) :: handle

RETURN
END SUBROUTINE dr_hook

END MODULE yomhook
