! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module containing all routines relating to 2D photolysis
!  used in UKCA sub-model.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_Photolysis
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_phot2d

USE ukca_flupj_mod,      ONLY: ukca_flupj
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

USE photol_constants_mod, ONLY: pi_over_180 => const_pi_over_180

IMPLICIT NONE

PRIVATE

PUBLIC ukca_photin

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_PHOT2D'

CONTAINS

!-----------------------------------------------------------------------

SUBROUTINE ukca_photin(row_lengthda, tot_p_rows, p_levelsda,                   &
                       jppj, nolat, nolev, nlphot, ntphot,                     &
                       sinlat, pl, pjin2d, photol_rates_2d)

! Purpose: Subroutine to interpolate 2D photolysis rates on
!          3-d latitudes, heights. reconstruct daily curves (every 5 day
!          plus code to account for hour angle (longitude)
!
!          3 values are stored for each level and each latitude
!          symmetrical distribution plus zero values at dawn and
!          dusk gives 7 points altogether.
!          Based on PHOTIN.F from Cambridge TOMCAT model.

! dataset: 51(levels)x19(latitudes)x3(points)x74(every 5days)
!
! ---------------------------------------------------------------------
!

IMPLICIT NONE


INTEGER, INTENT(IN) :: row_lengthda
INTEGER, INTENT(IN) :: tot_p_rows
INTEGER, INTENT(IN) :: p_levelsda
INTEGER, INTENT(IN) :: jppj
! File dimensions
INTEGER, INTENT(IN) :: nolat
INTEGER, INTENT(IN) :: nolev
INTEGER, INTENT(IN) :: nlphot
INTEGER, INTENT(IN) :: ntphot

REAL, INTENT(IN) :: sinlat(:)
REAL, INTENT(IN) :: pl(:,:,:)
REAL, INTENT(IN) :: pjin2d(:,:,:,:) ! 2D photol rates

REAL, INTENT(OUT) :: photol_rates_2d(:,:,:,:,:)

! Local variables
REAL :: pr2d(nolev)     ! 2D model level pressures
REAL :: pr2dj(nlphot)   ! 2D photolysis level pressures

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PHOTIN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up 2-D pressure arrays

CALL ukca_inpr2d(nolev,nlphot,pr2d,pr2dj)

! Interpolate 2D photolysis rates onto 3-D levels and
! latitudes. Longitude comes later.

CALL ukca_interpj(pjin2d, pr2dj, pl, row_lengthda,                             &
             tot_p_rows, p_levelsda, jppj, nolat, nlphot, ntphot,              &
             sinlat, pi_over_180, photol_rates_2d)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_photin

!------------------------------------------------------------

SUBROUTINE ukca_inpr2d(nolev, nlphot, pr2d, pr2dj)
!
! Purpose: Subroutine to calculate the pressure levels for the
!          2D photolysis rates. Original version taken from the
!          Cambridge TOMCAT model.
!
!          Called from UKCA_PHOTIN.
!
! ---------------------------------------------------------------------
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: nolev
INTEGER, INTENT(IN) :: nlphot
REAL, INTENT(OUT)   :: pr2d(nolev)       ! Pressures on 2D model levels
REAL, INTENT(OUT)   :: pr2dj(nlphot)     ! Pressures on 2D photolysis level

!       Local variables

INTEGER, PARAMETER :: maxlev = 30

INTEGER :: j                        ! Loop variable
INTEGER :: ij                       ! Loop variable

REAL, PARAMETER :: fac  = 1.284025417  ! Used to ensure that pre
REAL, PARAMETER :: psur = 1.0e5        ! Surface pressure in Pas
REAL, PARAMETER :: ares = 3.0          ! Factor related to verti
REAL, PARAMETER :: eps  = 1.0e-10      ! Factor used to calculat

REAL :: ee                       ! Temporary store
REAL :: fj                       ! Factor related to vertical re
REAL :: za                       ! Factor related to vertical re
REAL :: pp(nlphot+1)
REAL :: pres(maxlev)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INPR2D'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO j = 1, nolev
  ee = EXP((j-1)/2.0)
  pr2d(j) = psur/(ee * fac)
END DO

!       2D pressure levels - normal 2-D order - up to level 30
!       for photolysis

DO j = 1,maxlev-1
  ee = EXP((j-1) / 2.0)
  pres(j) = psur / (ee * fac)
END DO
pres(maxlev) = pres(maxlev-1)

fj    = 2.0/ares
pp(1) = (1.0-fj)*LOG(psur)+fj*LOG(pres(1))
pp(1) = EXP(pp(1))

DO ij = 2,nlphot+1
  za     = ij/(2.0*ares)
  fj     = 2.0*za+0.5+eps
  j      = INT(fj)
  fj     = fj-j-eps
  j      = j+1
  pp(ij) = (1.0-fj)*LOG(pres(j-1))+ fj*LOG(pres(j))
  pp(ij) = EXP(pp(ij))
END DO

pr2dj(1)=(psur+pp(1))*0.5

DO ij = 2,nlphot
  pr2dj(ij) = (pp(ij)+pp(ij-1))*0.5
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_inpr2d

!----------------------------------------------------------------

SUBROUTINE ukca_interpj(pjin2d, pr2dj, pl, lon, lat, lev, jppj,                &
                        nolat, nlphot, ntphot, sinlat, degrad, pjin)
!
! Purpose: Subroutine to interpolate photolysis rates
!
!          Called from UKCA_PHOTIN.
!
! ---------------------------------------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN) :: lon           ! No of longitudes
INTEGER, INTENT(IN) :: lat           ! No of latitudes
INTEGER, INTENT(IN) :: lev           ! No of levels
INTEGER, INTENT(IN) :: jppj          ! No of species
! File dimensions
INTEGER, INTENT(IN) :: nolat
INTEGER, INTENT(IN) :: nlphot
INTEGER, INTENT(IN) :: ntphot

REAL, INTENT(IN)       :: pjin2d(:,:,:,:)  ! 2D photol rates
REAL, INTENT(IN)       :: pl(:,:,:)  ! local press
REAL, INTENT(IN)       :: pr2dj(:)   ! 2D photo
REAL, INTENT(IN)       :: sinlat(:)  ! Sine (3D
REAL, INTENT(IN)                     :: degrad  ! To conve

REAL, INTENT(OUT)      :: pjin(:,:,:,:,:)

!       Local variables

INTEGER :: i                     ! Loop variable
INTEGER :: ii                    ! Loop variable
INTEGER :: j                     ! Loop variable
INTEGER :: jr                    ! Loop variable
INTEGER :: k                     ! Loop variable
INTEGER :: kk                    ! Loop variable
INTEGER :: l                     ! Loop variable

REAL, PARAMETER :: Npole = 90.0

REAL :: p2d(nlphot)
REAL :: lat2d(nolat)          ! 2D model latitudes
REAL :: wks1(lat,nlphot)      ! Working array
REAL :: wks2(nolat)           ! Working array
REAL :: wks3(nlphot)          ! Working array
REAL :: wks4(lat,lev)         ! Working array
REAL :: lati
REAL :: press
REAL :: delphi

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INTERPJ'


!       Set up 2D latitudes. lat=1 pole nord
!       LAT2D() is the latitude in the centre of the 2D box in radians

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
delphi = Npole/nolat
DO i = 1,nolat
  lat2d(i)=SIN((90.0 - (2*i-1)*delphi)*degrad)
END DO

DO jr = 1,jppj              ! Loop over photolysis reactions

  !         Interpolate linearly in Sin(lat) (KK is the point through the

  DO kk = 1,ntphot          ! Loop over times of day

    DO j = 1,nlphot
      DO i = 1,nolat
        wks2(i) = pjin2d(i,j,kk,jr)
      END DO

      DO k = 1,lat
        lati = sinlat((k-1)*lon+1)
        lati = MAX(lati, lat2d(nolat))
        lati = MIN(lati, lat2d(    1))
        wks1(k,j) = ukca_flupj(lati,lat2d,wks2,nolat)
        wks1(k,j) = MAX(wks1(k,j), 0.0)
      END DO
    END DO

    !           Interpolate linearly in log(P)
    DO ii = 1,lon
      DO k = 1,lat
        DO j = 1,nlphot
          wks3(j) = wks1(k,j)
          p2d(j)  = LOG(pr2dj(j))
        END DO

        DO l = 1,lev
          press = LOG(pl(ii,k,l))
          press = MAX(press, p2d(nlphot))
          press = MIN(press, p2d(1))
          wks4(k,l) = ukca_flupj(press,p2d,wks3,nlphot)
          wks4(k,l) = MAX(wks4(k,l), 0.0)
        END DO
      END DO


      DO l = 1,lev
        DO k = 1,lat
          pjin(ii,k,l,kk,jr) = wks4(k,l)
        END DO
      END DO
    END DO


  END DO     ! End of loop over times of day
END DO       ! End of loop over photolysis reactions

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_interpj

END MODULE ukca_phot2d
