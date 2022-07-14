! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module containing subroutine to calculate the pressure of
!          the 2.0pvu surface and the pressure of the 380K surface.
!          Routine combines them to calculate the pressure of the
!          tropopause and set L_troposphere to .true. for those
!          gridboxes in the troposphere.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_tropopause

USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Diag, newline
USE ereport_mod, ONLY: ereport
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength
USE ukca_um_legacy_mod, ONLY: planet_radius
USE level_heights_mod, ONLY: r_theta_levels
USE missing_data_mod, ONLY: rmdi, imdi

IMPLICIT NONE
SAVE
PRIVATE

REAL, PARAMETER :: tpv = 2.0e-6        ! tropopause PV (pvu)
REAL, PARAMETER :: tpt = 380.0         ! tropopause theta (K)

INTEGER,ALLOCATABLE,PUBLIC  :: tropopause_level(:,:)

!     Pressures of the theta, pv, and combined tropopause

REAL,ALLOCATABLE,PUBLIC  :: p_tropopause(:,:)
REAL,ALLOCATABLE,PUBLIC  :: theta_trop(:,:)
REAL,ALLOCATABLE,PUBLIC  :: pv_trop(:,:)

!     Logical set to true for gridpoints within the troposphere

LOGICAL, ALLOCATABLE, PUBLIC :: L_troposphere(:,:,:)

CHARACTER(LEN=errormessagelength) :: cmessage           ! Error message
INTEGER           :: ierr               !   "   code

PUBLIC ukca_calc_tropopause

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_TROPOPAUSE'

CONTAINS

SUBROUTINE ukca_calc_tropopause(                                               &
  row_length, rows, model_levels,                                              &
  latitude, theta, pv, pr_boundaries, pr_levels)

!      Description:
!       Subroutine to calculate p_tropopause. This is a weighted
!       average of 2 tropopause definitions. In the extratropics
!       (>= 28 deg), p_tropopause is the pressure (in Pa) of the
!       2.0 PVU surface. In the tropics (<= 13.0 deg), it is the
!       pressure of the 380K isentropic surface. Between the two
!       latitudes, the weighting function is
!       W = A * sech (lat) + B * lat**2 + C and then
!       p_tropopause = W *(PV tropopause) + (1-W) *(380K tropopause)
!
!       This function is from Hoerling, Monthly Weather Review,
!       Vol 121 162-172. "A Global Analysis of STE during Northern
!       Winter."
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

REAL, PARAMETER :: lat_etropics = 28.0    ! define extratropics
REAL, PARAMETER :: lat_tropics  = 13.0    ! define tropics
REAL, PARAMETER :: fixed_pres   = 40000.0 ! default tropopause

REAL, INTENT(IN):: latitude(row_length,rows)           ! latitude (degrees N)
REAL, INTENT(IN):: theta(row_length,rows,model_levels) ! theta
REAL, INTENT(IN):: pv(row_length,rows,model_levels)    ! pv on model levs
REAL, INTENT(IN):: pr_boundaries(row_length,rows,0:model_levels)
                          ! pressure at layer boundaries
REAL, INTENT(IN):: pr_levels(row_length,rows,1:model_levels)
                          ! pressure at theta levels

!     Local variables

INTEGER                          :: i,j,l           ! Loop counters
INTEGER                          :: jll,jlu         ! Level indices

REAL                             :: thalph

REAL :: wt(row_length,rows)           ! weighting fn
REAL :: wth(model_levels)             ! theta profile
REAL :: wpl(model_levels)             ! pressure profile
REAL :: wpv(model_levels)             ! PV profile

REAL, PARAMETER :: max_z_for_tropopause_height = 30000.0
                   ! Max height (m) upto which to search for tropopause
INTEGER, SAVE   :: max_trop_level
                   ! level corresponding to max_z_for_tropopause_height

LOGICAL, SAVE   :: l_first = .TRUE.
LOGICAL         :: l_level_found

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALC_TROPOPAUSE'


!     Calculate weighting function from Hoerling, 1993 paper
!     W = a sech (lat) + B (lat)**2 + C

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise values for level indices
jll = 0
jlu = 0

! Initialise missing values
p_tropopause(:,:)     = rmdi
theta_trop(:,:)       = rmdi
pv_trop(:,:)          = rmdi
L_troposphere(:,:,:)  = .FALSE.
tropopause_level(:,:) = imdi

! If running with single level tropopause cannot be calculated
IF (model_levels > 1) THEN
  !!Calculate level corresponding to maximum tropopause height
  IF ( l_first ) THEN
    max_trop_level = model_levels-2        ! default
    LOOP_trop_1: DO l = 2,model_levels
      IF ((MINVAL(r_theta_levels(:,:,l)) - planet_radius) >                    &
         max_z_for_tropopause_height) THEN
        max_trop_level = l
        EXIT LOOP_trop_1
      END IF
    END DO LOOP_trop_1
    l_first = .FALSE.
  END IF

  DO j = 1,rows
    DO i = 1,row_length
      IF (ABS(latitude(i,j)) >= lat_etropics) THEN
        wt(i,j) = 1.0                                   ! extratropcs
      ELSE IF (ABS(latitude(i,j)) <= lat_tropics) THEN
        wt(i,j) = 0.0                                   ! tropics
      ELSE
        wt(i,j) = 5560.74*2                                                    &
              /(EXP(latitude(i,j))+EXP(-1.0*latitude(i,j)))                    &
              + 1.67e-3*(latitude(i,j))**2 - 0.307    ! sub-tropics
      END IF
    END DO
  END DO

  !     Calculate theta and pv tropopauses

  DO j = 1,rows
    DO i = 1,row_length
      DO l=1,model_levels
        wth(l) = theta(i,j,l)         ! theta profile
        wpl(l) = pr_levels(i,j,l)     ! pressure profile
        wpv(l) = pv(i,j,l)            ! PV profile
      END DO

      !         Find theta levels which straddle tpt

      l_level_found = .FALSE.
      LOOP_tropopause_down_1: DO l = max_trop_level,2,-1
        IF (wth(l) <= tpt .AND. wth(l+1) >= tpt) THEN
          jll = l
          jlu = l+1
          l_level_found = .TRUE.
          EXIT LOOP_tropopause_down_1
        END IF
      END DO LOOP_tropopause_down_1

      !         Calculate pressure of theta tropopause

      thalph = (tpt-wth(jll))/(wth(jlu)-wth(jll))
      theta_trop(i,j) = (1.0-thalph)*LOG(wpl(jll))                             &
                      + thalph*LOG(wpl(jlu))
      IF (wpl(jll) < 0.0 .OR. wpl(jlu) < 0.0 .OR. thalph < 0.0                 &
        .OR. thalph > 1.0) THEN
        IF (PrintStatus >= PrStatus_Diag) THEN
          DO l=1,model_levels
            WRITE(umMessage,'(3I6,3E12.3)') i,j,l,theta(i,j,l),                &
                       pr_levels(i,j,l), pv(i,j,l)
            CALL umPrint(umMessage,src='ukca_tropopause')
          END DO
          WRITE(umMessage,'(2I6,2E12.3)') jll,jlu,thalph,theta_trop(i,j)
          CALL umPrint(umMessage,src='ukca_tropopause')
        END IF
      END IF
      theta_trop(i,j) = EXP(theta_trop(i,j))

      IF ((.NOT. l_level_found) .AND. PrintStatus >= PrStatus_Diag) THEN
        WRITE(umMessage,'(A,E12.3,A,E12.3)')                                   &
          'Level scan for theta tropopause failed.' // newline //              &
          'Diagnosed tropopause pressure (', theta_trop(i,j),                  &
          'Pa) may be invalid.' // newline //                                  &
          'Weighting for theta tropopause is ', 1.0 - wt(i,j)
        CALL umPrint(umMessage,src=RoutineName)
      END IF

      !         Find potential vorticity levels which straddle tpv

      l_level_found = .FALSE.
      LOOP_tropopause_down_2: DO l = max_trop_level,2,-1
        IF (ABS(wpv(l)) <= tpv .AND. ABS(wpv(l+1)) >= tpv) THEN
          jll = l
          jlu = l+1
          l_level_found = .TRUE.
          EXIT LOOP_tropopause_down_2
        END IF
      END DO LOOP_tropopause_down_2

      !         Calculate pressure of pv tropopause

      thalph = (tpv-ABS(wpv(jll)))/(ABS(wpv(jlu))-ABS(wpv(jll)))
      pv_trop(i,j) = (1.0-thalph)*LOG(wpl(jll))                                &
                   + thalph*LOG(wpl(jlu))
      IF (wpl(jll) < 0.0 .OR. wpl(jlu) < 0.0 .OR. thalph < 0.0                 &
        .OR. thalph > 1.0) THEN
        IF (PrintStatus >= PrStatus_Diag) THEN
          WRITE(umMessage,'(A20)') 'UKCA_TROPOPAUSE:'
          CALL umPrint(umMessage,src='ukca_tropopause')
          DO l=1,model_levels
            WRITE(umMessage,'(3I6,3E12.3)') i,j,l,theta(i,j,l),                &
                       pr_levels(i,j,l), pv(i,j,l)
            CALL umPrint(umMessage,src='ukca_tropopause')
          END DO
        END IF
        cmessage = ' Difficulty diagnosing pv tropopause, '//                  &
                   'Reverting to default tropopause pressure'
        ierr = -1
        CALL ereport(RoutineName,ierr,cmessage)

        pv_trop(i,j) = LOG(fixed_pres)
      END IF

      pv_trop(i,j) = EXP(pv_trop(i,j))

      IF ((.NOT. l_level_found) .AND. PrintStatus >= PrStatus_Diag) THEN
        WRITE(umMessage,'(A,E12.3,A,E12.3)')                                   &
          'Level scan for PV tropopause failed.' // newline //                 &
          'Tropopause pressure (', pv_trop(i,j),                               &
          'Pa) may be invalid or default.' // newline //                       &
          'Weighting for PV tropopause is ', wt(i,j)
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    END DO
  END DO

  !     Calculate combined tropopause based on weighting function

  L_troposphere = .FALSE.
  DO j = 1,rows
    DO i = 1,row_length
      p_tropopause(i,j) = (wt(i,j)*pv_trop(i,j))                               &
                        + ((1.0-wt(i,j))*theta_trop(i,j))

      !         Check for whole gridboxes below tropopause and find
      !         the model level which contains the tropopause

      LOOP_trop_2: DO l=1,model_levels
        IF (pr_boundaries(i,j,l) >= p_tropopause(i,j)) THEN
          L_troposphere(i,j,l)  = .TRUE.
        ELSE
          tropopause_level(i,j) = l
          EXIT LOOP_trop_2
        END IF
      END DO LOOP_trop_2

    END DO
  END DO

ELSE
  ! Case where only one model level, default to being in troposphere
  L_troposphere(:,:,:) = .TRUE.
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calc_tropopause

END MODULE ukca_tropopause
