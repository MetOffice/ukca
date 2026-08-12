! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!   To save structure ukca_radaer_lfric
!
! ---------------------------------------------------------------------

MODULE ukca_radaer_lfric_saved_mod

USE ukca_radaer_lfric_struct_mod, ONLY:                                        &
    ukca_radaer_lfric_struct

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_RADAER_LFRIC_SAVED_MOD'

! Structure for UKCA/radiation interaction
TYPE (ukca_radaer_lfric_struct), SAVE :: ukca_radaer_lfric

END MODULE ukca_radaer_lfric_saved_mod
