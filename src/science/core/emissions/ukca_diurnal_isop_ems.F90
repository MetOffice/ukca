! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
! Description:
!  Main driver routine for chemistry
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_MAIN1.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Adapted from original code written by Guang Zeng
!
!------------------------------------------------------------------
!
MODULE ukca_diurnal_isop_ems_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_DIURNAL_ISOP_EMS_MOD'

CONTAINS

SUBROUTINE ukca_diurnal_isop_ems(row_length,                                   &
                                 rows,                                         &
                                 emi_in,                                       &
                                 cosza_in,                                     &
                                 intza_in,                                     &
                                 sinlat,                                       &
                                 coslat,                                       &
                                 tanlat,                                       &
                                 timestep,                                     &
                                 emi_out,                                      &
                                 testdcycl)

USE ukca_constants,  ONLY: rhour_per_day, rsec_per_hour,  pi, pi_over_180,     &
                           recip_pi_over_180, fxb, fxc
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Diag

USE ukca_time_mod, ONLY: i_day_number

IMPLICIT NONE


INTEGER,                              INTENT(IN)  :: row_length
INTEGER,                              INTENT(IN)  :: rows
LOGICAL,                              INTENT(IN)  :: testdcycl

REAL,                                 INTENT(IN)  :: timestep
REAL, INTENT(IN)  :: cosza_in(1:row_length,1:rows)    ! COS (zenith angle)
REAL, INTENT(IN)  :: intza_in(1:row_length,1:rows)    ! INT(COS(sza))
REAL, INTENT(IN)  :: sinlat(1:row_length,1:rows)      ! sin (latitude)
REAL, INTENT(IN)  :: coslat(1:row_length,1:rows)      ! cos (latitude)
REAL, INTENT(IN)  :: tanlat(1:row_length,1:rows)      ! tan (latitude)
REAL, INTENT(IN)  :: emi_in(1:row_length,1:rows)      ! IN isoprene emission
REAL, INTENT(OUT) :: emi_out(1:row_length,1:rows)
                                    ! OUT diurnally varying isoprene emission

!     Local variables

INTEGER :: i,j

REAL :: declin                   ! Solar declination angle
REAL :: tan_declin               ! TAN(solar declination angle)
REAL :: int_a,b,int_h,sza_int    ! SZA integration variables
REAL :: emit_day                 ! C5H8 emission in one day
REAL :: trise                    ! Time of the sunrise
REAL :: mpd, mod_mpd

REAL  :: cosza(1:row_length,1:rows)       ! COS (zenith angle)
REAL  :: intza (1:row_length,1:rows)      ! INT(COS(sza))
REAL  :: daylen(1:row_length,1:rows)      ! day length for curent day
REAL  :: cs_hour_ang(1:row_length,1:rows) ! cosine hour angle


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DIURNAL_ISOP_EMS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Calculate Declination Angle and Daylength for each row for
!     curent day of the year.
!     Ensure COS of HOUR ANGLE does not exceed + or - 1, and set DAY
!     LENGTH of 1st & last rows to be the same as that of adjacent rows
!     to avoid possible problems at the poles (tan(90)=infinity).

daylen(:,:)      = 0.0
cs_hour_ang(:,:) = 0.0
emi_out(:,:)     = 0.0

cosza(:,:)    = cosza_in(:,:)
intza(:,:)    = intza_in(:,:)

declin = fxb * SIN(pi_over_180*(266.0+i_day_number))
tan_declin = TAN(declin)
DO j = 1,rows
  DO i = 1,row_length

    IF (ABS(tan_declin) < EPSILON(0.0)) THEN
      cs_hour_ang(i,j) = 0.0
    ELSE
      cs_hour_ang(i,j) = MAX(-1.0, MIN(1.0, -1.0*tanlat(i,j) * tan_declin))
    END IF

    mpd     = (fxc * ACOS(cs_hour_ang(i,j)))*60.0 ! compute minutes of
                                                  ! sunshine per day
    mod_mpd = MOD(mpd,(timestep/60.0)) ! compute residual mins
                                       ! (in excess to mins per timestep)
    daylen(i,j) = (mpd - mod_mpd)/60.0 ! compute sunshine hours
                                       ! a multiple of timestep

    trise = 12.0 - 0.5*daylen(i,j)

    int_a   = sinlat(i,j)*SIN(declin)*daylen(i,j)
    b       = coslat(i,j)*COS(declin)
    int_h   = (24.0/pi)*SIN(pi*trise/12.0)
    sza_int = int_a + b*int_h

    !         adjust factor to model time step units

    sza_int = sza_int*rsec_per_hour

    !         calculate emission(day)

    emit_day = emi_in(i,j)*rsec_per_hour*rhour_per_day

    !         now scale the emissions

    IF ((cosza(i,j) > 0.0) .AND. (sza_int > 1.0e-1)) THEN
      emi_out(i,j) = emit_day*(cosza(i,j)/intza(i,j))
      IF (emi_out(i,j) < 0.0) THEN
        emi_out(i,j) = 0.0
      END IF
    ELSE
      emi_out(i,j) = 0.0
    END IF

    IF ((j == 1) .AND. (i == 1) .AND.                                          &
        (emi_out(i,j) > 0.0)     .AND.                                         &
        PrintStatus >= PrStatus_Diag .AND.                                     &
        (testdcycl)) THEN

      WRITE(umMessage,'(A8,2A5,6A15,3A12)')                                    &
              'UKCA_DIURNAL_ISOP_EMS:',                                        &
              'i', 'j',                                                        &
              'cont_daylen',                                                   &
              'disc_daylen',                                                   &
              'deg_lat',                                                       &
              'SZA',                                                           &
              'INT(SZA)',                                                      &
              'INT(SZA) prec',                                                 &
              'emi_in',                                                        &
              'emi_out',                                                       &
              'emi_day'
      CALL umPrint(umMessage,src='ukca_diurnal_isop_ems')

      WRITE(umMessage,'(A8,2I5,6F15.1,3E12.4)')                                &
              i, j,                                                            &
              (fxc * ACOS(cs_hour_ang(i,j))),                                  &
              daylen(i,j),                                                     &
              ACOS(coslat(i,j))*recip_pi_over_180,                             &
              ACOS(cosza(i,j))*recip_pi_over_180,                              &
              sza_int,                                                         &
              intza(i,j),                                                      &
              emi_in(i,j),                                                     &
              emi_out(i,j),                                                    &
              emit_day
      CALL umPrint(umMessage,src='ukca_diurnal_isop_ems')
    END IF
  END DO                  !end domain loop
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_diurnal_isop_ems
END MODULE ukca_diurnal_isop_ems_mod
