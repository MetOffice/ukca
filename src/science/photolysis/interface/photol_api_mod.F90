! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Application program interface (API) module for UKCA Photolysis.
!
! Method:
!
!  This module provides access to subroutines and parameters
!  required by a parent or component model for running Photolysis.
!  It acts as a collation point for components of the API defined
!  in other Photolysis modules rather than including any definitions
!  itself.
!
!  Note that 'photol_api_mod' should be the only Photolysis module used
!  by
!  a parent application.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA/Photolysis
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE photol_api_mod

! Procedures and parameters made available below constitute the formal
! Photolysis API. All names made available begin with 'photol_' for
! clarity.
! This is important to avoid pollution of a parent application's
! namespace.

USE photol_calc_ozonecol_mod, ONLY: photol_calc_ozonecol

IMPLICIT NONE

END MODULE photol_api_mod
