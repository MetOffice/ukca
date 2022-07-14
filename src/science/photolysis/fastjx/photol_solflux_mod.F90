! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine to calculate solar flux, using average solar flux, the shape
!  of the solar cycle, and the time series of the solar cycle.
!
! Method: A singular value decomposition has been performed on the spectrally
!  resolved solar irradiance data. The data is very well described by the
!  product of a spectrally varying term and a time series. The spectral
!  composition of the solar cycle has been fed through the FAST-JX binning
!  algorithm. Here we compute the product of the spectral signature
!  and the solar cycle time series to infer the spectral variation.
!
!  If i_ukca_solcyc is set to (1) the observational data are used for the time
!  period for which they exist, before and after this time period
!  an average cycle is repeated. The period of this repeated average cycle
!  is 10 years 8 months.
!  If i_ukca_solcyc is set to (2) this average cycle is used at all times.
!
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE photol_solflux_mod
USE ukca_config_specification_mod, ONLY: ukca_config
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PHOTOL_SOLFLUX_MOD'

CONTAINS

SUBROUTINE photol_solflux(current_time, tau, lookup, solcyc, x)

USE fastjx_data, ONLY: solcyc_av, solcyc_ts, n_solcyc_av, n_solcyc_ts
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: current_time(7) ! current model time
REAL, INTENT(IN)    :: tau         ! FastJX time (in hours, FJX only)

LOGICAL, INTENT(IN) :: lookup      ! Are we calculating for:
                                   !   lookup tables (T)
                                   !   FastJX (F)

REAL, INTENT(IN)    :: solcyc(:)   ! spectral component of sol cyc var

REAL, INTENT(OUT)   :: x(:)        ! fl or quanta mod modification


! Local variables

REAL :: realtime        ! time in years
REAL :: sol_ph          ! phase of periodic cycle
INTEGER :: obs_idx
INTEGER :: cyc_idx
REAL :: frac

! Dependent on calendar and number of days in the year:
! number of days in each month (Jan-Dec) dependent
INTEGER :: days_in_month(12)
! other local time variables
INTEGER :: i_year
INTEGER :: i_month
INTEGER :: i_day

REAL :: reftime        ! reference time for periodic cycle
REAL :: init_time      ! initial time for the solar cycle
REAL :: obs_end_year   ! end year for the observations

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHOTOL_SOLFLUX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! initialise to zero
x = 0.0

! set local time variables
i_year = current_time(1)
i_month = current_time(2)
i_day = current_time(3)

! January of the start year of the observations
init_time = REAL(ukca_config%i_ukca_solcyc_start_year) + (0.5/12.0)

! reftime for the periodic cycle is given as December 2004 (t=0)
reftime  = 2004.0 + (11.5/12.0)

! end year for the observations, n_solcyc_ts is in months
obs_end_year = REAL(ukca_config%i_ukca_solcyc_start_year) +                    &
     (REAL(n_solcyc_ts)/12.0)


! get the correct number of days in the months
IF (ukca_config%l_cal360) THEN
  days_in_month = [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]
ELSE
  IF (days_in_year(i_year) < 366) THEN
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  ELSE
    days_in_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  END IF
END IF


IF (lookup) THEN

   ! doing for lookup tables - need to do from the middle of the month ONLY
   ! works for both 360-day and Gregorian calendars

   ! calculate real time in years to the middle of the month only for
   ! interpolation. As this is only refreshed every month must
   ! only interpolate to/from middle of months
  realtime = REAL(i_year) + ((REAL(i_month-1) + 0.5)/12.0)

ELSE

   ! Doing for FJX, so can interpolate for each chemical timestep
   ! need to calculate real time in years.

   ! Real time in years done using the current year, the current day,
   ! and the current time used in FJX calculations, tau (in hours,
   ! calculated in UKCA_FASTJX)
   ! Works for both 360-day and Gregorian calendars, where the progress
   ! the month (i.e. from XXXX.00 -> XXXX.(1/12) etc.) is referenced to
   ! the middle of that month. This means that for a 360-day calendar 12
   ! midnight on the 16th of the month would be equivalent to 0.5 through the
   ! month, but for a Gregorian calendar, e.g. a 31-day month 12 noon on the
   ! 16th would be equivalent to 0.5 through the month.
   ! This is to ensure that time interpolation works correctly in all calendars
   ! and when considering the periodic/average cycle.
  realtime = REAL(i_year) + ((REAL(i_month-1) +                                &
      (REAL(i_day-1)/REAL(days_in_month(i_month))))/12.0) +                    &
      (tau/(24.0*REAL(days_in_year(i_year))))

END IF


! if outside of obs or using the periodic/average cycle, interpolate
! between end of obs and avg cyc.
! Calculate phase of cycle, phase 0 set to Dec 2004 (reftime),
! cycle period = 10y8m (128m). Need to -1.0 in the modulo to ensure start
! time is correct.
sol_ph = MODULO(((realtime-reftime)*12.0)-1.0,REAL(n_solcyc_av))
cyc_idx = FLOOR(sol_ph) + 1

! modify solar flux with solar cycle
IF ((ukca_config%i_ukca_solcyc==1) .AND.                                       &
    realtime>=REAL(ukca_config%i_ukca_solcyc_start_year)                       &
    .AND. realtime<=obs_end_year) THEN
  obs_idx = FLOOR((realtime - init_time)*12.0) + 1
  frac = (realtime - init_time)*12.0 - (obs_idx-1)
  IF (obs_idx < 1) THEN
    ! if at start of obs, interpolate between avg. cycle and beginning
    ! of obs
    x = solcyc * ((1.0 - frac) * solcyc_av(cyc_idx) + frac * solcyc_ts(1))
  ELSE IF (obs_idx >= n_solcyc_ts) THEN
    x = solcyc * ((1.0 - frac) * solcyc_ts(n_solcyc_ts) +                      &
         frac * solcyc_av(cyc_idx))
  ELSE
    x = solcyc * ((1.0 - frac) * solcyc_ts(obs_idx) +                          &
         frac * solcyc_ts(obs_idx+1))
  END IF
ELSE ! using periodic avg. cycle
  frac = sol_ph - (cyc_idx-1)
  IF (cyc_idx < n_solcyc_av) THEN
    x = solcyc * ((1.0 - frac) * solcyc_av(cyc_idx) +                          &
         frac * solcyc_av(cyc_idx+1))
  ELSE
    x = solcyc * ((1.0 - frac) * solcyc_av(cyc_idx) + frac * solcyc_av(1))
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE photol_solflux

INTEGER FUNCTION days_in_year(year)
! ----------------------------------------------------------------------
! Description:
!   Returns the year length in days based on the model calendar.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Function argument
INTEGER, INTENT(IN) :: year

IF (ukca_config%l_cal360) THEN
  ! 360-day calendar
  days_in_year = 360
ELSE
  ! Gregorian calendar; look out for leap years
  IF (MOD(year,4) /= 0) THEN
    days_in_year = 365
  ELSE IF (MOD(year,100) /= 0) THEN
    days_in_year = 366
  ELSE IF (MOD(year,400) /= 0) THEN
    days_in_year = 365
  ELSE
    days_in_year = 366
  END IF
END IF

RETURN
END FUNCTION days_in_year

END MODULE photol_solflux_mod
