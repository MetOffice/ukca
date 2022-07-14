! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Field to hold dry diameter
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: UKCA
!
! Code Description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.

MODULE ukca_drydiam_field_mod

IMPLICIT NONE

REAL, ALLOCATABLE :: drydiam(:,:,:,:)
! Geometric mean dry diameter of particles in each mode (m)

END MODULE ukca_drydiam_field_mod
