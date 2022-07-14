! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module to provide availability checking of UKCA environmental driver
!   inputs
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds,
! University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_environment_check_mod

IMPLICIT NONE

PUBLIC


CONTAINS


! ----------------------------------------------------------------------
SUBROUTINE check_environment(n_fld_present, n_fld_missing)
! ----------------------------------------------------------------------
! Description:
!   Check availability of required UKCA environmental driver fields
!   and print details if any fields are missing
! ----------------------------------------------------------------------

USE ukca_environment_req_mod, ONLY: check_environment_availability,            &
                                    ukca_get_environment_varlist
USE ukca_fieldname_mod, ONLY: maxlen_fieldname

USE umPrintMgr, ONLY: umMessage, umPrint

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(OUT) :: n_fld_present   ! No. of required fields present
INTEGER, INTENT(OUT) :: n_fld_missing   ! No. of required fields missing

! Local variables
LOGICAL, POINTER :: available(:) ! Availability flags associated with
                                 ! list of required fields
CHARACTER(LEN=maxlen_fieldname), POINTER :: fieldnames(:)
                                 ! List of required fields by name
INTEGER :: i                     ! Field index
CHARACTER(LEN=7) :: i_char       ! Character variable for field index
CHARACTER(LEN=8) :: status_label ! Label to indicate presence/absence of field
INTEGER :: error_code            ! Error code (error if > 0)

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CHECK_ENVIRONMENT'

! Get field availability record
CALL check_environment_availability(n_fld_present, n_fld_missing,              &
                                    availability_ptr=available)

! If there are missing fields then print a list showing all required fields
! with availability status

IF (n_fld_missing > 0) THEN

  WRITE(umMessage,'(A)')                                                       &
    '*** ONE OR MORE REQUIRED UKCA ENVIRONMENTAL DRIVER FIELDS ARE UNSET ***'
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I4)') 'Number of required fields present: ', n_fld_present
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I4)') 'Number of required fields missing: ', n_fld_missing
  CALL umPrint(umMessage,src=RoutineName)

  CALL ukca_get_environment_varlist(fieldnames, error_code)
  IF (error_code > 0) THEN
    WRITE(umMessage,'(A)') 'The list of required fields is unavailable'
    CALL umPrint(umMessage,src=RoutineName)
  ELSE
    WRITE(umMessage,'(A7,1X,A8,1X,A)') '-INDEX-','-STATUS-','-FIELDNAME-'
    CALL umPrint(umMessage,src=RoutineName)
    DO i = 1, SIZE(fieldnames)
      IF (available(i)) THEN
        status_label = 'PRESENT'
      ELSE
        status_label = 'MISSING'
      END IF
      WRITE(i_char,'(I7)') i
      WRITE(umMessage,'(A7,1X,A8,1X,A)') ADJUSTL(i_char), status_label,        &
        fieldnames(i)
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF

END IF

RETURN
END SUBROUTINE check_environment

END MODULE ukca_environment_check_mod
