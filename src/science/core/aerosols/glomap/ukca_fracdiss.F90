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
MODULE ukca_fracdiss_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_FRACDISS_MOD'

CONTAINS

SUBROUTINE ukca_fracdiss(row_length, rows, model_levels, temp, p_theta_levels, &
                         rh, qcl, frac_diss,kp_nh, H_plus_3d_arr)


USE asad_mod,             ONLY: ddhr, dhr, kd298, k298, jpeq, jpdw
USE ukca_um_legacy_mod,   ONLY: rmol
USE chemistry_constants_mod, ONLY: avogadro
USE ukca_constants,       ONLY: m_air
USE parkind1,             ONLY: jprb, jpim
USE yomhook,              ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length     ! No of points per row
INTEGER, INTENT(IN) :: rows           ! No of rows
INTEGER, INTENT(IN) :: model_levels   ! No of levels

REAL, INTENT(IN) :: temp(row_length,rows,model_levels)           ! Temperature
REAL, INTENT(IN) :: p_theta_levels(row_length,rows,model_levels) ! Pressure (Pa)
REAL, INTENT(IN) :: rh(row_length,rows,model_levels)  ! Relative Humidity
REAL, INTENT(IN) :: qcl(row_length,rows,model_levels) ! Cloud liquid water
                                                      !   (kg/kg)
! 3D pH array
REAL, INTENT(IN) :: H_plus_3d_arr(row_length, rows, model_levels)

REAL, INTENT(OUT) :: frac_diss(row_length,rows,model_levels,                   &
                               jpdw,jpeq+1)
REAL, INTENT(OUT) :: kp_nh(row_length,rows,model_levels)

! Local variables
INTEGER :: k                    ! Loop variable
INTEGER :: ns                   ! Loop variable

REAL, PARAMETER :: qcl_min=1.0e-12  ! do calcs when qcl > qcl_min

REAL :: inv_298                        ! 1/298
REAL :: tmp1(row_length,rows,model_levels)        ! Temporary variable
REAL :: tmp2(row_length,rows,model_levels)        ! Temporary variable
REAL :: Henry(row_length,rows,model_levels)       ! Henry's Law consts
REAL :: Henry_eff(row_length,rows,model_levels)   ! Effective  ---"---
REAL :: tlev(row_length,rows,model_levels)        ! Temperature on model levels
REAL :: qcla(row_length,rows,model_levels)        ! qcl in kg/m^3
LOGICAL :: todo(row_length,rows,model_levels)     ! defines qcl > qcl_min

! Aqueous equilibrium constants:
REAL :: eq_con(row_length,rows,model_levels,jpeq)

! For ammonium nitrate dissociation
REAL :: rhd(row_length,rows,model_levels)      ! rhd deliquesance RH
REAL :: kp(row_length,rows)                    ! kp
REAL :: kp1(row_length,rows)                   ! kp1
REAL :: p1(row_length,rows)                    ! p1
REAL :: p2(row_length,rows)                    ! p2
REAL :: p3(row_length,rows)                    ! p3

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_FRACDISS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

inv_298 = 1.0/298.0

frac_diss(:,:,:,:,:)=0.0
todo(:,:,:) = .FALSE.
DO k=1,model_levels
  WHERE (qcl(:,:,k) > qcl_min)
    todo(:,:,k) = .TRUE.
    tlev(:,:,k) = temp(:,:,k)
    qcla(:,:,k)=qcl(:,:,k)*M_air*p_theta_levels(:,:,k)/                        &
                (Rmol*temp(:,:,k))
  END WHERE
END DO

DO ns=1,jpdw
  IF (kd298(ns,2) > 1e-12) THEN      ! 2nd dissociation constant
    WHERE (todo(:,:,:))
      tmp1(:,:,:) = (1.0/tlev(:,:,:)) - inv_298
      Henry(:,:,:) = k298(ns) * EXP(dhr(ns)*tmp1(:,:,:))
      eq_con(:,:,:,1) = kd298(ns,1) * EXP(ddhr(ns,1)*tmp1(:,:,:))
      eq_con(:,:,:,2) = kd298(ns,2) * EXP(ddhr(ns,2)*tmp1(:,:,:))
      ! Replace Global pH with 3D values from H_plus array
      Henry_eff(:,:,:)=Henry(:,:,:)*(1.0+eq_con(:,:,:,1)/H_plus_3d_arr(:,:,:))
      tmp2(:,:,:) = Rmol*qcla(:,:,:)*tlev(:,:,:)
      frac_diss(:,:,:,ns,1) = 1.0/(1.0 + (1.013e5/                             &
                              (tmp2(:,:,:)*Henry(:,:,:))))
      frac_diss(:,:,:,ns,2) = 1.0/(1.0 + (1.013e5/                             &
                             (tmp2(:,:,:)*Henry_eff(:,:,:))))
      frac_diss(:,:,:,ns,2) = frac_diss(:,:,:,ns,2)-                           &
                              frac_diss(:,:,:,ns,1)
      ! Replace Global pH with 3D values from H_plus array
      Henry_eff(:,:,:) = Henry(:,:,:) * (1.0 + eq_con(:,:,:,1)/                &
               H_plus_3d_arr(:,:,:) + eq_con(:,:,:,1)*eq_con(:,:,:,2)/         &
               H_plus_3d_arr(:,:,:)**2)
      frac_diss(:,:,:,ns,3) = 1.0/(1.0 + (1.013e5/                             &
                             (tmp2(:,:,:)*Henry_eff(:,:,:))))
      frac_diss(:,:,:,ns,3)=frac_diss(:,:,:,ns,3)-                             &
               frac_diss(:,:,:,ns,2) - frac_diss(:,:,:,ns,1)
    END WHERE
  ELSE IF (kd298(ns,1) > 1e-12) THEN  ! one dissociation
    WHERE (todo(:,:,:))
      tmp1(:,:,:) = (1.0/tlev(:,:,:)) - inv_298
      Henry(:,:,:) = k298(ns) * EXP(dhr(ns)*tmp1(:,:,:))
      eq_con(:,:,:,1) = kd298(ns,1) * EXP(ddhr(ns,1)*tmp1(:,:,:))
      ! Replace Global pH with 3D values from H_plus array
      Henry_eff(:,:,:)=Henry(:,:,:)*(1.0+eq_con(:,:,:,1)/H_plus_3d_arr(:,:,:))
      tmp2(:,:,:) = Rmol*qcla(:,:,:)*tlev(:,:,:)
      frac_diss(:,:,:,ns,1) = 1.0/(1.0 + (1.013e5/                             &
                              (tmp2(:,:,:)*Henry(:,:,:))))
      frac_diss(:,:,:,ns,2) = 1.0/(1.0 + (1.013e5/                             &
                              (tmp2(:,:,:)*Henry_eff(:,:,:))))
      frac_diss(:,:,:,ns,2) = frac_diss(:,:,:,ns,2)-                           &
                              frac_diss(:,:,:,ns,1)
    END WHERE
  ELSE                                  ! No dissociation
    WHERE (todo(:,:,:))
      tmp1(:,:,:) = (1.0/tlev(:,:,:)) - inv_298
      Henry(:,:,:) = k298(ns) * EXP(dhr(ns)*tmp1(:,:,:))
      tmp2(:,:,:)=qcla(:,:,:)*Henry(:,:,:)*Rmol*tlev(:,:,:)
      frac_diss(:,:,:,ns,1) = 1.0/(1.0 + (1.013e5/tmp2(:,:,:)))
    END WHERE
  END IF
END DO

! Calculate NH4NO3 dissociation constant
DO k=1,model_levels
  rhd(:,:,k) = EXP((618.3/temp(:,:,k))-2.551)
END DO

kp_nh(:,:,:)=0.0
todo(:,:,:)=.FALSE.
todo = (rh < rhd)
DO k=1,model_levels
  WHERE (todo(:,:,k))
    kp(:,:)=(temp(:,:,k)**(-6.025))*EXP(118.87-24084.0/                        &
             temp(:,:,k))
  ELSE WHERE
    p1(:,:)=EXP(-135.94+(8763/temp(:,:,k)))*temp(:,:,k)**19.12
    p2(:,:)=EXP(-122.65+(9969/temp(:,:,k)))*temp(:,:,k)**16.22
    p3(:,:)=EXP(-182.61+(13875/temp(:,:,k)))*temp(:,:,k)**24.46
    kp1(:,:)=(temp(:,:,k)**(-6.025))*EXP(118.87-24084.0/                       &
              temp(:,:,k))
    kp(:,:)=(p1(:,:)-p2(:,:)*(1.0-rh(:,:,k))+p3(:,:)*                          &
            (1.0-rh(:,:,k))**2)*((1.0-rh(:,:,k))**1.75)*kp1(:,:)
  END WHERE
  kp_nh(:,:,k)=kp(:,:)*1.0e-8/((rmol*1.0e6*temp(:,:,k)/                        &
               avogadro)**2)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_fracdiss
END MODULE ukca_fracdiss_mod
