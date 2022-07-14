! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  To perform an explicit loss of methane in the top two model
!  levels when tropospheric chemistry is switched on
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_ch4_stratloss_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CH4_STRATLOSS_MOD'

CONTAINS

SUBROUTINE ukca_ch4_stratloss(n_be_calls, n_pnts,                              &
                    vol, dts,                                                  &
                    y, stratloss)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE chemistry_constants_mod, ONLY: avogadro
IMPLICIT NONE

INTEGER, INTENT(IN) :: n_be_calls ! No. chemical steps
INTEGER, INTENT(IN) :: n_pnts     ! Actual no. calculations

REAL, INTENT(IN)    :: dts
REAL, INTENT(IN)    :: vol(n_pnts)

REAL, INTENT(IN OUT) :: y(n_pnts)

REAL, INTENT(IN OUT) :: stratloss(n_pnts)

!     Local variables

INTEGER :: n

REAL, PARAMETER :: Loss_rate = 2.0e-7   ! per second
REAL :: p(n_pnts)
REAL :: l(n_pnts)
REAL :: yp(n_pnts)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CH4_STRATLOSS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO n = 1, n_be_calls

  yp(:) = y(:)
  p(:)  = 0.0
  l(:)  = Loss_rate
  y(:)  = (yp(:)+dts*p(:))/(1.0+dts*l(:))

  !       Calculate flux terms.

  stratloss(:) = stratloss(:)                                                  &
                + Loss_rate*y(:)*dts*Vol(:)/avogadro

END DO  ! n_be_calls

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_ch4_stratloss
END MODULE ukca_ch4_stratloss_mod
