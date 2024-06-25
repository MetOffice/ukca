! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE setzen_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     This routine calculates the zenith angle grid used, also max zenith
!     angle and sets up temperature grid and O3 factor grid.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SETZEN_MOD'

CONTAINS
SUBROUTINE setzen(tabang,zenmax,altc,tabt,tabo3)

USE photol_constants_mod, ONLY: pi => const_pi,                                &
                                planet_radius => const_planet_radius

USE ukca_parpho_mod, ONLY: jplev, jpchi, jps90, jpchin,                        &
                           jptem, jpo3p, jps90,                                &
                           szamax, tmin, tmax, o3min, o3max
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)  :: altc(jplev)    ! Altitude of each level centre in km

REAL, INTENT(OUT) :: tabang(jpchi)  ! Solar zenith angle in radians
REAL, INTENT(OUT) :: tabt(jptem)    ! Temperature
REAL, INTENT(OUT) :: tabo3(jpo3p)   ! O3 factor

! Mazimum zenith angle at each altitude at which direct sunlight is
! received.
REAL, INTENT(OUT) :: zenmax(jplev)

! Local variables
REAL :: re  ! = planet_radius/1.0e3
REAL :: cosint
REAL :: radint
REAL :: tdif
REAL :: o3dif

INTEGER :: j
INTEGER :: jc
INTEGER :: jt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETZEN'



!     Spacing of angles in lookup table
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
cosint=1.0/REAL(jpchi-1-jps90)
radint=(szamax - 90.0)*pi/(180.0*REAL(jps90))

!     First angle = 0 degrees
tabang(1) = 0.0

!     From 0 to 90 degrees
DO j=2,jpchin-1
  jc = ABS(j-jpchin)
  tabang(j) = ACOS(cosint*REAL(jc))
END DO

!     Angle 90 degrees
tabang(jpchin) = 0.5*pi

DO j = jpchin+1,jpchi
  tabang(j) = 0.5*pi + REAL(j-jpchin)*radint
END DO

!     Maximum zenith angle
re = planet_radius/1.0e3
DO j=1,jplev
  zenmax(j) = (pi - ASIN(re/(re + altc(j))))
END DO

!     Temperature grid
IF (jptem == 1) THEN
  tabt(    1)=0.5*(tmin + tmax)
ELSE
  tabt(    1)=tmin
  tabt(jptem)=tmax
  tdif=(tmax-tmin)/REAL(jptem-1)
  DO jt=2,jptem-1
    tabt(jt) = tmin + tdif*REAL(jt-1)
  END DO
END IF

!     O3 grid
IF (jpo3p == 1) THEN
  tabo3(   1)=1.0
ELSE
  tabo3(    1)=o3min
  tabo3(jpo3p)=o3max
  o3dif=(o3max-o3min)/REAL(jpo3p-1)
  DO jt=2,jpo3p-1
    tabo3(jt) = o3min + o3dif*REAL(jt-1)
  END DO
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE setzen
END MODULE setzen_mod
