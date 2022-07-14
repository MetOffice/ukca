! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module for handling UKCA time variables.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds,
! University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_time_mod

USE ukca_config_specification_mod, ONLY: ukca_config

IMPLICIT NONE

PUBLIC

INTEGER :: i_year                 ! Current model time (year)
INTEGER :: i_month                ! Current model time (month)
INTEGER :: i_day                  ! Current model time (day)
INTEGER :: i_hour                 ! Current model time (hour)
INTEGER :: i_minute               ! Current model time (minute)
INTEGER :: i_second               ! Current model time (second)
INTEGER :: i_day_number           ! Current model time (day of year)

! Time variables at previous time step
INTEGER :: i_year_previous        ! Previous model time (year)
INTEGER :: i_hour_previous        ! Previous model time (hour)
INTEGER :: i_minute_previous      ! Previous model time (minute)
INTEGER :: i_second_previous      ! Previous model time (second)
INTEGER :: i_day_number_previous  ! Previous model time (day of year)


CONTAINS


! ----------------------------------------------------------------------
SUBROUTINE set_time(time_data)
! ----------------------------------------------------------------------
! Description:
!   Sets the values of the time variables for the current time step
!   to the values in the input time array having 7 elements:
!   year, month, day, hour, minute, second, day of year
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine argument
INTEGER, INTENT(IN) :: time_data(7)

i_year = time_data(1)
i_month = time_data(2)
i_day = time_data(3)
i_hour = time_data(4)
i_minute = time_data(5)
i_second = time_data(6)
i_day_number = time_data(7)

RETURN
END SUBROUTINE set_time


! ----------------------------------------------------------------------
SUBROUTINE set_previous_time(time_data)
! ----------------------------------------------------------------------
! Description:
!   Sets the values of the time variables for the previous time step
!   to the required values from the input time array having 7 elements:
!   year, month, day, hour, minute, second, day of year
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine argument
INTEGER, INTENT(IN) :: time_data(7)

i_year_previous = time_data(1)
i_hour_previous = time_data(4)
i_minute_previous = time_data(5)
i_second_previous = time_data(6)
i_day_number_previous = time_data(7)

RETURN
END SUBROUTINE set_previous_time


! ----------------------------------------------------------------------
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

END MODULE ukca_time_mod
