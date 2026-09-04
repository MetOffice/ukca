! -----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! -----------------------------------------------------------------------------
!
! Purpose:
!   To save structure ukca_radaer_lfric
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
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
