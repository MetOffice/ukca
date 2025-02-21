! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Calculates the fraction of each species in each
!   dissolved state.
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
MODULE ukca_fdiss_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_FDISS_MOD'

CONTAINS

SUBROUTINE ukca_fdiss(n_points, qcl_min, t, p, qcl, fdiss, H_plus_1d_arr)


USE asad_mod,            ONLY: ddhr, dhr, kd298, k298, jpeq, jpdw
USE ukca_constants,      ONLY: m_air
USE ukca_config_constants_mod,  ONLY: rmol
USE parkind1,            ONLY: jprb, jpim
USE yomhook,             ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_points                 ! No of points

REAL, INTENT(IN) :: qcl_min                     ! do calcs when qcl > qcl_min
REAL, INTENT(IN) :: t(n_points)                 ! Temperature (K)
REAL, INTENT(IN) :: p(n_points)                 ! Pressure (Pa)
REAL, INTENT(IN) :: qcl(n_points)               ! Cloud liquid water (kg/kg)

REAL, INTENT(OUT) :: fdiss(n_points, jpdw, jpeq+1)

! 1-D pH array to be used in ukca_fdiss
REAL, INTENT(IN) :: H_plus_1d_arr(n_points)

! Local variables
INTEGER :: ns                   ! Loop variable

REAL, PARAMETER :: p0      = 1.0135e5   ! P0 in Pa
REAL, PARAMETER :: inv_298 = 1.0/298.0  ! 1/298

REAL :: tmp1(n_points)                  ! Temporary variable
REAL :: tmp2(n_points)                  ! Temporary variable
REAL :: Henry(n_points)                 ! Henry's Law consts
REAL :: Henry_eff(n_points)             ! Effective  ---"---
REAL :: qcla(n_points)                  ! qcl in kg/m^3
LOGICAL :: todo(n_points)               ! defines qcl > qcl_min

! Aqueous equilibrium constants:
REAL :: eq_con(n_points,jpeq)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_FDISS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


fdiss(:,:,:) = 0.0
todo(:) = .FALSE.
WHERE (qcl(:) > qcl_min)
  todo(:) = .TRUE.
  qcla(:) = qcl(:)*M_air*p(:)/(Rmol*t(:))
  tmp1(:) = (1.0/t(:)) - inv_298
  tmp2(:) = Rmol*qcla(:)*t(:)
END WHERE

DO ns=1,jpdw
  IF (kd298(ns,2) > 1e-12) THEN      ! 2nd dissociation constant
    WHERE (todo(:))
      Henry(:) = k298(ns) * EXP(dhr(ns)*tmp1(:))
      eq_con(:,1) = kd298(ns,1) * EXP(ddhr(ns,1)*tmp1(:))
      eq_con(:,2) = kd298(ns,2) * EXP(ddhr(ns,2)*tmp1(:))
      ! Replace Global pH with 1D values from H_plus array
      Henry_eff(:) = Henry(:)*(1.0+eq_con(:,1)/H_plus_1d_arr(:))
      fdiss(:,ns,1) = 1.0/(1.0 + (p0/(tmp2(:)*Henry(:))))
      fdiss(:,ns,2) = 1.0/(1.0 + (p0/(tmp2(:)*Henry_eff(:))))
      fdiss(:,ns,2) = fdiss(:,ns,2) - fdiss(:,ns,1)
      ! Replace Global pH with 1D values from H_plus array
      Henry_eff(:) = Henry(:) * (1.0 + eq_con(:,1)/                            &
               H_plus_1d_arr(:) + eq_con(:,1)*eq_con(:,2)/H_plus_1d_arr(:)**2)
      fdiss(:,ns,3) = 1.0/(1.0 + (p0/(tmp2(:)*Henry_eff(:))))
      fdiss(:,ns,3) = fdiss(:,ns,3)- fdiss(:,ns,2) - fdiss(:,ns,1)
    END WHERE
  ELSE IF (kd298(ns,1) > 1e-12) THEN  ! one dissociation
    WHERE (todo(:))
      Henry(:) = k298(ns) * EXP(dhr(ns)*tmp1(:))
      eq_con(:,1) = kd298(ns,1) * EXP(ddhr(ns,1)*tmp1(:))
      ! Replace Global pH with 1D values from H_plus array
      Henry_eff(:)=Henry(:)*(1.0+eq_con(:,1)/H_plus_1d_arr(:))
      tmp2(:) = Rmol*qcla(:)*t(:)
      fdiss(:,ns,1) = 1.0/(1.0 + (p0/(tmp2(:)*Henry(:))))
      fdiss(:,ns,2) = 1.0/(1.0 + (p0/(tmp2(:)*Henry_eff(:))))
      fdiss(:,ns,2) = fdiss(:,ns,2) - fdiss(:,ns,1)
    END WHERE
  ELSE                                  ! No dissociation
    IF (k298(ns) > 1e-12) THEN          ! avoid glitch with BrONO2
      WHERE (todo(:))
        Henry(:) = k298(ns) * EXP(dhr(ns)*tmp1(:))
        fdiss(:,ns,1) = 1.0/(1.0 + (p0/(tmp2(:)*Henry(:))))
      END WHERE
    END IF
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_fdiss
END MODULE ukca_fdiss_mod
