! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module containing the photol_curve subroutine called in
!  photol_ctl.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE photol_curve_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PHOTOL_CURVE_MOD'

CONTAINS

SUBROUTINE photol_curve(pjinda, tloc, dayl, tot_p_rows, row_length, jppj, wks)
!
! Purpose: Subroutine to interpolate tropospheric photolysis rates
!          in time. Based on curve.F from Cambridge TOMCAT model.
!
!          Called from PHOTOL_CTL.
!
! ---------------------------------------------------------------------
!
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim
USE umPrintMgr,          ONLY: umPrint, umMessage
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

INTEGER, INTENT(IN) :: tot_p_rows                ! no of rows
INTEGER, INTENT(IN) :: row_length                ! no of cols
INTEGER, INTENT(IN) :: jppj                      ! number of chemical species

REAL, INTENT(IN) :: dayl(:)                      ! day length
REAL, INTENT(IN) :: tloc(:)                      ! local time
REAL, INTENT(IN) :: pjinda(:,:,:)                ! 2D photolys

REAL, INTENT(OUT) :: wks(:,:)                    ! interpolated

! Local variables

INTEGER :: i                                       ! loop variab
INTEGER :: j                                       ! loop variables
INTEGER :: k                                       ! loop variables
INTEGER :: jr                                      ! loop variables
INTEGER :: errcode                                 ! Variable passed to ereport

REAL, PARAMETER :: tfrac1 = 0.04691008             ! determines
REAL, PARAMETER :: tfrac2 = 0.23076534             ! 2D photolys

REAL :: dawn                ! time of dawn
REAL :: dusk                ! time of dusk
REAL :: timel               ! local time
REAL :: slope               ! slope used in linear interpolation
REAL :: const               ! intercept used in linear interpola
REAL :: fgmt(7)             ! times at which photol rates are va


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHOTOL_CURVE'

! Initialise wks

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
wks = 0.0

! Calculate rates using a simple linear interpolation.

DO j = 1,tot_p_rows
  DO i = 1,row_length
    k = i+(j-1)*row_length

    ! Non-Polar night

    IF (dayl(k) > 0.0) THEN
      dawn = 12.00 - (dayl(k)/2.0)
      dusk = 12.00 + (dayl(k)/2.0)

      fgmt(1) = dawn
      fgmt(2) = dawn + tfrac1*dayl(k)
      fgmt(3) = dawn + tfrac2*dayl(k)
      fgmt(4) = 12.00
      fgmt(5) = dawn + (1.0-tfrac2)*dayl(k)
      fgmt(6) = dawn + (1.0-tfrac1)*dayl(k)
      fgmt(7) = dusk

      timel = tloc(k)

      IF (timel > 24.0) timel = timel-24.0

      timel = MIN(timel,24.0)
      timel = MAX(timel, 0.0)

      ! Local Night-time

      IF (timel < dawn .OR. timel > dusk) THEN

        wks(k,1:jppj) = 0.0

        ! For the time between dawn and PJIN(1) or PJIN(5) and dusk

      ELSE IF ((timel >= dawn   .AND. timel < fgmt(2))                         &
           .OR. (timel > fgmt(6) .AND. timel <= dusk)) THEN

        IF (timel > fgmt(6)) timel = 24.00 - timel

        ! trap for -ve (timel-fgmt(1))
        IF ((fgmt(1) - timel) < 1.0e-6) timel = fgmt(1)

        DO jr = 1, jppj
          slope = pjinda(j,1,jr)/(fgmt(2) - fgmt(1))
          wks(k,jr) = slope*(timel - fgmt(1))

          IF (wks(k,jr) < 0.0) THEN
            WRITE(umMessage,'(A,F12.6)')                                       &
            'negative wks in photol_curve 1',wks(k,jr)
            CALL umPrint(umMessage,src='ukca_phot2d')
            WRITE(umMessage,'(I6,I6,I6,I6,F12.6,F12.6,F12.6,F12.6)')           &
            i,j,k,jr,slope,fgmt(1),fgmt(2),timel
            CALL umPrint(umMessage,src='ukca_phot2d')
            WRITE(umMessage,'(F12.6,F12.6)') pjinda(j,1,jr),fgmt(1)-timel
            CALL umPrint(umMessage,src='ukca_phot2d')
            WRITE(umMessage,'(F12.6,F12.6,F12.6)') fgmt,dawn,dusk
            CALL umPrint(umMessage,src='ukca_phot2d')
            errcode = jr
            CALL ereport('PHOTOL_CURVE',errcode,' Negative photolysis')
          END IF
        END DO

        ! For the time between PJIN(1) and PJIN(2) or PJIN(4) and PJIN(5)

      ELSE IF ((timel >= fgmt(2) .AND. timel < fgmt(3))                        &
           .OR. (timel >  fgmt(5) .AND. timel <= fgmt(6))) THEN

        IF (timel > fgmt(5)) timel = 24.00 - timel

        DO jr = 1, jppj
          slope = (pjinda(j,2,jr)-pjinda(j,1,jr))/                             &
               (fgmt(3) - fgmt(2))
          const = pjinda(j,1,jr)- slope* fgmt(2)
          wks(k,jr)= slope*timel + const
        END DO

        ! For the time between PJIN(2), PJIN(3) and PJIN(4)

      ELSE IF (timel >= fgmt(3) .AND. timel <= fgmt(5)) THEN

        IF (timel > fgmt(4)) timel = 24.00 - timel

        DO jr = 1, jppj
          slope = (pjinda(j,3,jr)-pjinda(j,2,jr))/                             &
              (fgmt(4) - fgmt(3))
          const = pjinda(j,2,jr)- slope* fgmt(3)
          wks(k,jr)= slope*timel + const
        END DO

      END IF    ! end of IF (timel < dawn .OR. timel > dusk)

      ! End of the condition on the non polar night

    ELSE

      wks(k,1:jppj) = 0.0

    END IF     ! end of IF (dayl(k) > 0.0)

  END DO       ! end of looping over row length
END DO         ! end of looping over rows

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE photol_curve

END MODULE photol_curve_mod