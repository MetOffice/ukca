! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Purpose: To calculate dry deposition rates for use in UKCA chemistry.
!           Original version from Cambridge TOMCAT model.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!           Called from UKCA_CHEMISTRY_CTL.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_ddeprt_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_DDEPRT_MOD'

CONTAINS

SUBROUTINE ukca_ddeprt(p_fieldda, bl_levels, i_day_number, i_mon, i_hour,      &
                       r_minute, secs_per_step, lon, lat, tanlat, dzl, z0m,    &
                       u_s, temp, dryrt)

USE ukca_um_legacy_mod, ONLY: vkman
USE ukca_constants,     ONLY: pi_over_180, fxb, fxc
USE asad_mod, ONLY: depvel, jpdd
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: p_fieldda             ! No of spatial points
INTEGER, INTENT(IN) :: bl_levels             ! No of boundary layer levels
INTEGER, INTENT(IN) :: i_day_number          ! Day of year
INTEGER, INTENT(IN) :: i_mon                 ! Month of year
INTEGER, INTENT(IN) :: i_hour                ! Hour of day
REAL, INTENT(IN) :: r_minute                 ! Minute of hour
REAL, INTENT(IN) :: secs_per_step            ! Time step (s)
REAL, INTENT(IN) :: lon(p_fieldda)           ! Longitude (degrees)
REAL, INTENT(IN) :: lat(p_fieldda)           ! Latitude (degrees)
REAL, INTENT(IN) :: tanlat(p_fieldda)        ! tan(latitude)
REAL, INTENT(IN) :: dzl(p_fieldda,bl_levels) ! Boundary layer thickness (m)
REAL, INTENT(IN) :: z0m(p_fieldda)           ! Roughness length (m)
REAL, INTENT(IN) :: u_s(p_fieldda)           ! Surface friction velocity (m/s)
REAL, INTENT(IN) :: temp(p_fieldda)          ! Surface temperature (K)

REAL, INTENT(OUT) :: dryrt(p_fieldda,jpdd)   ! Dry deposition rate (/s)

!       Local variables

INTEGER :: i                                 ! Loop variable
INTEGER :: ns                                ! Loop variable
INTEGER :: iday                              ! Integer indicatin
INTEGER :: isum                              ! Integer indicatin
INTEGER :: icat                              ! Integer indicatin

REAL, PARAMETER :: midday = 12.0             ! Midday time
REAL, PARAMETER :: temp_summer = 275.0       ! Temp threshold fo
REAL, PARAMETER :: z0m_min = 1.0e-3          ! Min z0m used to g
REAL, PARAMETER :: z0m_max = 1.0e-1          ! Max z0m used to g

REAL :: dawn                                 ! Time of dawn
REAL :: dusk                                 ! Time of dusk
REAL :: timel                                ! Local time
REAL :: vdep1                                ! Depn velocity at
REAL :: vdeph                                ! Depn velocity at

REAL :: tgmt                   ! GMT time (decimal representation)
REAL :: tan_declin             ! TAN(declination)
REAL :: cs_hour_ang            ! cosine hour angle
REAL :: tloc  (p_fieldda)      ! local time
REAL :: daylen(p_fieldda)      ! local daylength

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DDEPRT'


!       Initialise dryrt to zero

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
dryrt(:,:) = 0.0

DO i = 1,p_fieldda  ! Loop over spatial scale

  ! Calculate local time as function of longitude
  tgmt = REAL(i_hour) + (r_minute / 60.0) + (secs_per_step * 0.5 / 3600.0)
  tloc(i) = tgmt + (24.0 * lon(i) / 360.0)
  IF (tloc(i) > 24.0) tloc(i) = tloc(i) - 24.0

  ! Calculate declination angle and daylength for current day of year
  ! Ensure cos of hour angle does not exceed + or - 1
  tan_declin = TAN(fxb * SIN(pi_over_180 * (266.0 + i_day_number)))
  IF (ABS(tan_declin) < EPSILON(0.0)) THEN
    cs_hour_ang = 0.0
  ELSE
    cs_hour_ang = MAX(-1.0, MIN(1.0, -tanlat(i) * tan_declin))
  END IF
  daylen(i) = fxc * ACOS(cs_hour_ang)

END DO

!       Set up deposition arrays for each species

DO ns = 1,jpdd                     ! Loop over deposited species
  DO i = 1,p_fieldda               ! Loop over spatial scale

    !           Calculate time of dawn and dusk for particular day.

    dawn  = midday - (daylen(i)/2.0)
    dusk  = midday + (daylen(i)/2.0)
    timel = tloc(i)

    !           Local night-time/day-time

    IF (timel < dawn .OR. timel > dusk) THEN
      iday = 0
    ELSE
      iday = 1
    END IF

    !           Find out if grid point is summer/winter depending on tempera

    IF (temp(i) > temp_summer) THEN
      isum = 1
    ELSE IF (temp(i) <= temp_summer) THEN
      isum = 0
    END IF

    !           Find out what category (water,forest,grass,desert)
    !           desert not used at present

    IF (z0m(i) < z0m_min) THEN
      icat = 1                      ! water/sea - 0.001

    ELSE IF (z0m(i) > z0m_max) THEN
      icat = 2                      ! forests - 0.1

    ELSE IF (z0m(i) <= z0m_max .OR. z0m(i) >= z0m_min) THEN
      icat = 3                      ! all other lands, grass

    END IF

    IF (lat(i) >= 0.0) THEN         ! Northern hemisphere

      !             Overwrite switch if there is ice cover - code from 2-D mod

      IF ((i_mon == 12 .OR. i_mon <= 3) .AND. lat(i) >= 45.0)                  &
        icat = 5
      IF ((i_mon == 11 .OR. i_mon == 4) .AND. lat(i) >= 55.0)                  &
        icat = 5
      IF ((i_mon == 10 .OR. i_mon == 5 .OR. i_mon == 6)                        &
      .AND. lat(i)>=65.0)                                                      &
        icat = 5
      IF ((i_mon >= 7 .AND. i_mon <= 9) .AND. lat(i) >= 70.0)                  &
        icat = 5

    ELSE IF (lat(i) < 0.0) THEN      ! Southern hemisphere

      IF ((i_mon == 12 .OR. i_mon <= 4) .AND. lat(i) <= -75.0)                 &
        icat = 5
      IF ((i_mon >= 5 .AND. i_mon <= 11) .AND. lat(i) <= -70.0)                &
        icat = 5

    END IF

    !           Select appropriate dry deposition velocity - 1m values
    !           and convert from cm/s to m/s

    IF (iday == 1 .AND. isum == 1) THEN       ! Summer day
      vdep1 = depvel(1,icat,ns)/100.0

    ELSE IF (iday == 0 .AND. isum == 1) THEN  ! Summer night
      vdep1 = depvel(2,icat,ns)/100.0

    ELSE IF (iday == 1 .AND. isum == 0) THEN  ! Winter day
      vdep1 = depvel(4,icat,ns)/100.0

    ELSE IF (iday == 0 .AND. isum == 0) THEN  ! Winter night
      vdep1 = depvel(5,icat,ns)/100.0

    END IF

    !           Extrapolate to the middle of lowest model layer

    vdeph = vdep1/(1.0 + vdep1*LOG(dzl(i,1)/2.0)/(vkman*u_s(i)))
    dryrt(i,ns) = vdeph/dzl(i,1)

  END DO                     ! End of loop over spatial points
END DO                       ! End of loop over deposited species

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_ddeprt
END MODULE ukca_ddeprt_mod
