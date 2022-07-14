! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module providing parameters for UKCA error handling.
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

MODULE ukca_error_mod

IMPLICIT NONE
PUBLIC

INTEGER, PARAMETER :: maxlen_message = 133
INTEGER, PARAMETER :: maxlen_procname = 80
INTEGER, PARAMETER :: errcode_ukca_uninit = 1
INTEGER, PARAMETER :: errcode_tracer_req_uninit = 2
INTEGER, PARAMETER :: errcode_tracer_mismatch = 3
INTEGER, PARAMETER :: errcode_ntp_uninit = 4
INTEGER, PARAMETER :: errcode_ntp_mismatch = 5
INTEGER, PARAMETER :: errcode_env_req_uninit = 6
INTEGER, PARAMETER :: errcode_env_field_unknown = 7
INTEGER, PARAMETER :: errcode_env_field_mismatch = 8
INTEGER, PARAMETER :: errcode_env_field_missing = 9
INTEGER, PARAMETER :: errcode_array_unallocated = 10
INTEGER, PARAMETER :: errcode_value_missing = 11

END MODULE ukca_error_mod
