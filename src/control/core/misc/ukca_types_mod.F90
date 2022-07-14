! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module containing kind definitions.
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

MODULE ukca_types_mod

IMPLICIT NONE

PRIVATE

! The definitions below are replicated from the UM module um_types
! (at UM vn12.2)

! Range for 32 bit integer
INTEGER, PARAMETER :: irange32 = 9
! Kind for 32 bit integer
INTEGER, PARAMETER, PUBLIC :: integer32 = SELECTED_INT_KIND(irange32)
! Range for small logicals
INTEGER, PARAMETER :: lrange1=1
! Kind for 32 bit logicals. Note that there is no
! "selected_logical_kind", but using the equivalent integer kind is a
! workaround that works on every platform tested.
INTEGER, PARAMETER, PUBLIC :: logical32 = integer32
! Kind for small logicals
INTEGER, PARAMETER, PUBLIC :: log_small = SELECTED_INT_KIND(lrange1)

END MODULE ukca_types_mod
