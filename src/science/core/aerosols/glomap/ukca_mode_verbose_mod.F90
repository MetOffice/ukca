! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module to store  Global variable for debug_level
!    passed to various GLOMAP routines

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  Fortran
!
! ######################################################################
MODULE ukca_mode_verbose_mod

IMPLICIT NONE

INTEGER, SAVE :: glob_verbose     ! Global variable for debug_level
                                  ! passed to various GLOMAP routines

END MODULE ukca_mode_verbose_mod
