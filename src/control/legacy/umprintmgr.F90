! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description:
!
!   UM legacy interface to UKCA printing.
!
! Method:
!
!   This is a simple prototype version.
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

MODULE umprintmgr

IMPLICIT NONE

PUBLIC

! The following public variables/parameters are all used by
! the legacy UM code in UKCA, so their names/values should not change

! ========================
! PARAMETERS
! ========================

INTEGER, PARAMETER            :: PrStatus_Min    = 1  ! Minimum output
INTEGER, PARAMETER            :: PrStatus_Normal = 2  ! Short info output
INTEGER, PARAMETER            :: PrStatus_Oper   = 3  ! Full info output
INTEGER, PARAMETER            :: PrStatus_Diag   = 4  ! Diagnostic output

! parameterisation of storage/buffering ammounts
INTEGER, PARAMETER            :: maxLineLen=1024

! Declare newline character
CHARACTER(LEN=1), PARAMETER   :: newline=NEW_LINE('a')

! ========================
! A buffer that clients
! can use for list
! directed writes
! ========================
CHARACTER (LEN=maxLineLen)    :: umMessage
!$OMP THREADPRIVATE (umMessage)

! ========================
! Runtime variables and
! initial defaults.
! ========================
INTEGER                       :: PrintStatus = PrStatus_Diag

CONTAINS

SUBROUTINE umprint(line, level, pe, src, model,                                &
                   UsrPrefix, HangIndent, stdErrorToo)

IMPLICIT NONE

CHARACTER(LEN=*)            :: line
INTEGER, OPTIONAL           :: level
INTEGER, OPTIONAL           :: pe
CHARACTER(LEN=*), OPTIONAL  :: src
CHARACTER(LEN=*), OPTIONAL  :: model
CHARACTER(LEN=*), OPTIONAL  :: UsrPrefix       ! User prefix for each line
CHARACTER(LEN=*), OPTIONAL  :: HangIndent      ! Hanging indent for new lines
LOGICAL, OPTIONAL           :: stdErrorToo

WRITE(*,'(A,A)') 'UKCA INFO:', TRIM(line)

END SUBROUTINE umprint

! Dummy routine
SUBROUTINE umprintflush
IMPLICIT NONE
END SUBROUTINE umprintflush

END MODULE umprintmgr
