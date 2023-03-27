! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  UM legacy interface to UKCA error reporting.
!
! Method:
!
!  This is a simple prototype version.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ereport_mod

IMPLICIT NONE

PUBLIC ereport

CONTAINS

SUBROUTINE ereport(routine_name, error_status, message)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: routine_name
INTEGER, INTENT(IN OUT) :: error_status
CHARACTER(LEN=*), INTENT(IN) :: message

! Take action depending on UM error status convention
IF (error_status > 0) THEN
  WRITE(*,*) 'UKCA ERROR in ', routine_name, ':', message
  STOP
ELSE IF (error_status < 0) THEN
  WRITE(*,*) 'UKCA WARNING in', routine_name, ':', message
ELSE IF (error_status == 0) THEN
  WRITE(*,*) 'UKCA INFO in', routine_name, ':', message
END IF

! Reset error_status
error_status = 0

END SUBROUTINE ereport

END MODULE ereport_mod
