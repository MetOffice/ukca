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
!    This subroutine sediments solid NAT particles (treated as solid HNO3).
!    The fall velocity is relatively slow for NAT and much faster for NAT/ice,
!    i.e., at lower temperatures. Method and parameters are from SLIMCAT.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE ukca_sediment_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_SEDIMENT_MOD'

CONTAINS

SUBROUTINE ukca_sediment(rows, row_length, model_levels,                       &
                  shno3, qcf, r_theta_levels, mass, secs_per_step,             &
                  L_stratosphere)
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Subroutine interface

INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: model_levels
REAL   , INTENT(IN) :: secs_per_step

REAL   , INTENT(IN) :: r_theta_levels(row_length, rows,                        &
                                   0:model_levels)
REAL   , INTENT(IN) :: mass(row_length, rows, model_levels)
REAL   , INTENT(IN) :: qcf(row_length, rows, model_levels)
LOGICAL, INTENT(IN) :: L_stratosphere(row_length, rows, model_levels)
REAL, INTENT(IN OUT) :: shno3(row_length, rows, model_levels)

! Local variables
! Slow sedimentation in the presence of NAT only, fast in the presence of ice
! fall velocity for NAT PSCs, from SLIMCAT:
REAL, PARAMETER :: fv1 = 0.000463  ! m/s
! fall velocity for NAT/ICE PSCs, from SLIMCAT
REAL, PARAMETER :: fv2 = 0.01736   ! m/s

INTEGER :: i,j,k

REAL :: thick(row_length, rows)
REAL :: sed_mmr(row_length, rows)
REAL :: frac(row_length, rows)
REAL :: psc1(row_length, rows)
REAL :: psc2(row_length, rows)
REAL :: sh2o(row_length, rows, model_levels)
! Local copy of ice cloud fraction- masked out in troposphere.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SEDIMENT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! do sedimentation. Don't do it at bottom level.
! any limiting of NAT formation will already have been considered in the
! amounts of shno3

! set sh2o=0.0 in the troposphere, elsewhere to frozen cloud fraction
DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      sh2o(i,j,k) = qcf(i,j,k)
    END DO
  END DO
END DO
WHERE (.NOT. L_stratosphere) sh2o = 0.0

DO k=2,model_levels
  IF (MAXVAL(shno3(:,:,k)) > 0.0) THEN

    ! set presence of PSC1/PSC2 indicators
    psc1 = 0.0
    psc2 = 0.0
    WHERE (shno3(:,:,k) > 0.0) psc1 = 1.0
    WHERE (sh2o (:,:,k) > 0.0) psc2 = 1.0

    ! calculate model level thickness
    thick = r_theta_levels(:,:,k) - r_theta_levels(:,:,k-1)

    ! calculate relative fraction of HNO3S to be moved into lower level
    frac = (psc1*fv1 + psc2*fv2) * secs_per_step / thick

    ! limit fraction to 1.
    WHERE (frac > 1.0) frac = 1.0

    ! calculate fraction in MMR or VMR to be sedimented
    sed_mmr = shno3(:,:,k) * frac

    ! subtract from layer under consideration
    shno3(:,:,k  ) = shno3(:,:,k  ) - sed_mmr

    ! add to next layer below.
    shno3(:,:,k-1) = shno3(:,:,k-1) +                                          &
                     sed_mmr * mass(:,:,k) / mass(:,:,k-1)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_sediment
END MODULE ukca_sediment_mod
