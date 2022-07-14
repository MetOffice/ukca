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
!  Description:
!   Calculate solar zenith angle for photolysis applications.
!   Adapted from UM SOLANG routine, but cos(SZA) is not
!   assumed to be 0 at night.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE photol_solang_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'PHOTOL_SOLANG_MOD'

CONTAINS

SUBROUTINE photol_solang(sindec, t, dt, eqt, sinlat, coslat, longit,           &
                        k, cosz)

USE planet_constants_mod, ONLY: s2r
USE conversions_mod, ONLY: pi, pi_over_180
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE


INTEGER, INTENT(IN) ::  k         ! Number of points

REAL, INTENT(IN)    :: sindec     ! Sin(solar declination)
REAL, INTENT(IN)    :: t          ! Start time (GMT)
REAL, INTENT(IN)    :: dt         ! timestep
REAL, INTENT(IN)    :: eqt        ! Eqn of time

REAL, INTENT(IN)    :: sinlat(k)  ! sin(latitude)
REAL, INTENT(IN)    :: coslat(k)  ! cos(latitude)
REAL, INTENT(IN)    :: longit(k)  ! longitude (degrees)

REAL, INTENT(OUT)   :: cosz(k)    ! Mean Cos(sza)

!     Local variables

INTEGER          :: j         ! Loop counter over points

!!!      REAL, PARAMETER  :: s2r = pi/43200.0  ! sec-to-rads converter

REAL             :: sinsin    ! Products of the sines and of the cosines
REAL             :: coscos    ! of solar declination and of latitude.
REAL             :: hat       ! Local hour angle at the start time.
REAL             :: difsin    ! A difference-of-sines intermediate value
REAL             :: trad      ! Start and length of timestep (T & DT)
REAL             :: dtrad     ! converted to radians after midday GMT

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHOTOL_SOLANG'

REAL :: cosdec

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
trad  =  t * s2r - pi
dtrad = dt * s2r
cosdec = SQRT(1.0 - sindec**2)

DO j = 1, k
  sinsin  = sindec * sinlat(j)
  coscos  = cosdec * coslat(j)
  hat     = pi_over_180 * longit(j) + trad + eqt
  difsin  = SIN(hat + dtrad) - SIN(hat)
  cosz(j) = difsin*coscos/dtrad + sinsin
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE photol_solang
END MODULE photol_solang_mod
