! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE settab_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate photolysis factors

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SETTAB_MOD'

CONTAINS
SUBROUTINE settab(alt,altc,do2c,do3c,tempc,                                    &
  drsc,tabs,tabang,albedo,lscat,                                               &
  tspo2,scs,ao2,ao2sr,ao3,dalt)
USE photol_constants_mod, ONLY: planet_radius => const_planet_radius
USE ukca_parpho_mod, ONLY: jplevp1, jplev, jpchi, jps90, jpchin,               &
                           jpwav, jplo, jphi, jps90
! Module procedures
USE cso2o3_mod,   ONLY: cso2o3
USE invert_mod,   ONLY: invert
USE ei2_mod,      ONLY: ei2
USE ei3_mod,      ONLY: ei3
USE isrchfgt_mod, ONLY: isrchfgt

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
LOGICAL, INTENT(IN) :: lscat
REAL, INTENT(IN) :: albedo             !     Ground albedo.
REAL, INTENT(IN) :: alt(jplevp1)
REAL, INTENT(IN) :: altc(jplev)
REAL, INTENT(IN) :: do2c(jplev)
REAL, INTENT(IN) :: do3c (jplev)
REAL, INTENT(IN) :: tempc(jplev)
REAL, INTENT(IN) :: drsc(jplev)
REAL, INTENT(OUT):: tabs (jplev,jpchi,jpwav)
REAL, INTENT(IN) :: tabang(jpchi)
REAL, INTENT(OUT):: tspo2(jplev,jpchi)
REAL, INTENT(IN) :: scs(jpwav)
REAL, INTENT(IN) :: ao2(jpwav)
! Olaf check, ao2sr and ao3 are intent(inout) in CSO2O3
REAL, INTENT(IN OUT) :: ao2sr(jpwav)
REAL, INTENT(IN OUT) :: ao3(jpwav)
REAL, INTENT(IN) :: dalt (jplev)

! Local variables
INTEGER :: i
INTEGER :: j
INTEGER :: jc
INTEGER :: jci
INTEGER :: jk
INTEGER :: jn
INTEGER :: jt
INTEGER :: jw
REAL :: re  ! = planet_radius/1.0e3
REAL :: alpha
REAL :: alta
REAL :: altacm
REAL :: altcur
REAL :: altb
REAL :: altbcm
REAL :: alt90
REAL :: arg
REAL :: beta
REAL :: spl
REAL :: spl2
REAL :: tauk
REAL :: taun
REAL :: tauv
REAL :: teak
REAL :: tean
REAL :: teav
REAL :: twoamu
REAL :: teanp
REAL :: taunp

INTEGER :: ijt(jps90,jplev)
INTEGER :: indx(jplev)

!     Slant path between centre of levels.
REAL :: tablen(jplev,jplev,jpchin)
REAL :: deptha(jplev,jpchi)
REAL :: depths(jplev,jpchi)
REAL :: tablenoa(jplev,jplev,jps90)
REAL :: tablenob(jplev,jplev,jps90)


!     Vertical optical depths etc.
REAL :: teac(jplevp1)
REAL :: tauc(jplevp1)
REAL :: jac (jplev)
REAL :: tea (jplevp1)
REAL :: tau (jplevp1)

!     Enhancement factor table.
REAL :: tabs0(jplev,jpchi,jpwav)

REAL :: b(jplev,jplev)
REAL :: binv (jplev,jplev)
REAL :: bcopy(jplev,jplev)
REAL :: a(jplev,jplev)
REAL :: delta(jplev,jplev)

REAL :: tmp1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETTAB'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

re = planet_radius/1.0e3

!     Set up path lengths (cm) up to and including 90 degrees.
!$OMP PARALLEL DEFAULT(NONE)                                                   &
!$OMP PRIVATE(jc,i,alpha,j,alta,altb,beta,altacm,altbcm,spl2,spl)              &
!$OMP PRIVATE(jci,alt90,jt,altcur)                                             &
!$OMP SHARED(tabang,altc,alt,tablen,ijt,tablenoa,tablenob,tspo2,do2c,re)
!$OMP DO SCHEDULE(STATIC)
DO jc = 1, jpchin
  DO i = 1, jplev

    !           Zenith angle at the starting point.
    alpha = tabang(jc)

    !           Follow the path back in spherical geometry.
    DO j = i, jplev

      !              If at the starting level
      IF (j == i) THEN
        alta = altc(j)
      ELSE
        !                For all other levels.
        alta = alt (j)
      END IF

      altb = alt(j+1)

      !              Calculate zenith angle for next level up.
      beta = ASIN((alta+re)*SIN(alpha)/(altb+re))

      !              Calculate the slant path within this level in cm
      altacm = (alta+re)*1.0e5
      altbcm = (altb+re)*1.0e5
      spl2 = altacm*altacm + altbcm*altbcm -                                   &
        2.0*altacm*altbcm*COS(alpha-beta)
      spl  = SQRT(MAX(0.0,spl2))

      !              Distance in cm within level.
      tablen(i,j,jc) = spl

      !              Reinitialise variables.
      alpha = beta

    END DO
  END DO
END DO
!$OMP END DO


!     Find the level index for each zenith angle and each level for
!     which the zenith angle is 90 degrees.
!$OMP DO SCHEDULE(STATIC)
DO jc = jpchin + 1, jpchi
  jci = jc - jpchin
  DO j = 1, jplev
    alt90 = (altc(j)+re)*SIN(tabang(jc)) - re
    !           Altitudes monotonically increasing
    jt=isrchfgt(jplevp1,alt,alt90) - 1
    ijt(jci,j) = jt
  END DO
END DO
!$OMP END DO


!     When the zenith angle TABANG(JC) > 90.0 degrees.
!$OMP DO SCHEDULE(STATIC)
DO jc = jpchin + 1, jpchi

  !        A zenith angle index.
  jci = jc - jpchin

  DO j = 1, jplev

    !           Initialise to zero.
    spl = 0.0
    DO i = 1, jplev
      tablenoa(i,j,jci) = 0.0
      tablenob(i,j,jci) = 0.0
    END DO

    !           Level, JT, which contains tangent point.
    jt = ijt(jci,j)

    !           Check if light can get there.
    IF (jt > 0 .AND. j >= jt) THEN

      !              Start by calculating path lengths just to the left of
      !              the tangent point.

      !              Check if this part of the calculation is needed.
      IF (j >= jt+1) THEN

        !              Altitude where zenith angle is 90 degrees.
        alt90 = (altc(j)+re)*SIN(tabang(jc))

        !              Altitude of interface JT+1.
        altcur = alt(jt+1) + re

        !              Set zenith angle at current point.
        alpha = ASIN(alt90/altcur)

        !              Go from current level
        !              to the level above JT summing the path lengths.
        DO i = jt + 1, j - 1

          alta = alt(i  )
          altb = alt(i+1)

          ! Calc. zenith angle at the centre of the next level down.
          beta = ASIN((alta+re)*SIN(alpha)/(altb+re))

          ! Calculate the slant path between the centre of
          ! this level and the centre of the next level down in cm.
          altacm = (alta+re)*1.0e5
          altbcm = (altb+re)*1.0e5
          spl2 = altacm*altacm + altbcm*altbcm -                               &
            2.0*altacm*altbcm*COS(alpha-beta)
          spl  = SQRT(MAX(0.0,spl2))

          tablenoa(i,j,jci) = spl

          !                 Reintialise variables.
          !                 adjust zenith angle for spherical geometry.
          alpha = beta

        END DO

        !              End of check if this part of the calculation is needed.
      END IF

      !              Path length within level containing tangent point
      !              Case (A) J >= JT +1
      IF (j >= jt+1) THEN

        alt90 = (altc(j)+re)*SIN(tabang(jc))
        altcur = alt(jt+1) + re
        spl = 2.0e5*SQRT(altcur*altcur - alt90*alt90)

        !              Case (B) J = JT. Tangent height within level J
      ELSE IF (j == jt) THEN

        alt90 = (altc(j)+re)*SIN(tabang(jc))

        altcur = alt(jt+1) + re
        spl =       1.0e5*SQRT(altcur*altcur - alt90*alt90)

        altcur = altc(j) + re
        spl = spl + 1.0e5*SQRT(altcur*altcur - alt90*alt90)

      END IF

      tablenoa(jt,j,jci) = spl

      ! Now calculate path lengths to the right of the tangent point.

      ! Set zenith angle at the current point.
      alpha = ASIN(alt90/altcur)

      ! Now go from JT to the top of the atmosphere.
      DO i = jt + 1, jplev

        ! Initialise variables.
        alta = alt(i  )
        altb = alt(i+1)

        ! Calculate zenith angle for the next level up.
        beta = ASIN((alta+re)*SIN(alpha)/(altb+re))

        ! Calculate the slant path between this level and the
        ! next level down in cm.
        altacm = (alta+re)*1.0e5
        altbcm = (altb+re)*1.0e5
        spl2 = altacm*altacm + altbcm*altbcm -                                 &
          2.0*altacm*altbcm*COS(alpha-beta)
        spl  = SQRT(MAX(0.0,spl2))

        tablenob(i,j,jci) = spl

        ! Reintialise variables.
        ! Adjust zenith angle for spherical geometry.
        alpha = beta

      END DO
    END IF
  END DO
END DO
!$OMP END DO


!     Set up slant path O2 column. Used in AO2SR and ANO parameterisations.

!     Set up path lengths up to and including 90 degrees.
!$OMP DO SCHEDULE(STATIC)
DO jc = 1, jpchin
  DO i = 1, jplev

    tspo2(i,jc) = 0.0

    !           Follow the path back in spherical geometry.
    DO j = i, jplev
      tspo2(i,jc) = tspo2(i,jc) + tablen(i,j,jc)*do2c(j)
    END DO

  END DO
END DO
!$OMP END DO

!     When the zenith angle TABANG(JC) > 90.0 degrees.
!$OMP DO SCHEDULE(STATIC)
DO jc = jpchin + 1, jpchi

  !        A zenith angle index.
  jci = jc - jpchin

  DO j = 1, jplev

    !           Set O2 column to small value to prevent LOG problems in ACSSRW
    tspo2(j,jc) = 1.0e10

    !           Find level, JT, just below where ALPHA=90 degrees.
    jt = ijt(jci,j)

    !           Check if light can get there.
    IF (jt > 0 .AND. j >= jt) THEN

      ! Having found this level first go from the current level
      ! down to the level above JT summing the path lengths.
      DO i = jt + 1, j - 1
        tspo2(j,jc)=tspo2(j,jc) + tablenoa(i,j,jci)*do2c(i)
      END DO

      ! Now calculate the path length between the level
      ! JT + 1 as it stradles the zenith angle of 90 degrees.

      tspo2(j,jc) = tspo2(j,jc) + tablenoa(jt,j,jci)*do2c(jt)

      ! Now calculate path lengths to the right of the tangent point.

      ! Now go from JT to the top of the atmosphere.
      DO i = jt + 1, jplev
        tspo2(j,jc) = tspo2(j,jc) + tablenob(i,j,jci)*do2c(i)
      END DO

    END IF
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

!     Set up optical depths

!     Wavelength loop.
DO jw = jplo, jphi

  !        Up to and including 90 degrees.
  DO jc = 1, jpchin
    DO i = 1, jplev
      deptha(i,jc) = 0.0
      depths(i,jc) = 0.0
    END DO
  END DO

  DO jc = 1, jpchin
    DO j = 1, jplev

      ! Reset the O2 and O3 absorption cross sections.
      CALL cso2o3(ao2sr,ao3,tempc(j),tspo2(j,jc),jw,jpwav)
      tmp1= (do3c(j)*ao3(jw) + do2c(j)*(ao2(jw)+ao2sr(jw)))
      DO i = 1,j
        deptha(i,jc) = deptha(i,jc) + tablen(i,j,jc)* tmp1
        depths(i,jc) = depths(i,jc) + tablen(i,j,jc)* drsc(j)*scs(jw)
      END DO
    END DO
  END DO


  ! Greater than 90 degrees.
  DO jc = jpchin + 1, jpchi

    ! A zenith angle index.
    jci = jc - jpchin

    DO j = 1, jplev

      ! Initialise DEPTHA
      deptha(j,jc) = 0.0
      depths(j,jc) = 0.0

      ! Recall what JT was
      jt = ijt(jci,j)

      ! Having found this level first go from the current level
      ! down to the level above JT summing the path lengths.

      ! If the calculation is required:
      IF (jt > 0 .AND. j >= jt) THEN
        DO i = jt + 1, j - 1

          ! Reset the O2 and O3 absorption cross sections.
          CALL cso2o3(ao2sr,ao3,tempc(i),tspo2(i,jc),                          &
            jw,jpwav)

          deptha(j,jc) = deptha(j,jc) + tablenoa(i,j,jci)*                     &
            (do3c(i)*ao3(jw) + do2c(i)*(ao2(jw)+ao2sr(jw)))
          depths(j,jc) = depths(j,jc) + tablenoa(i,j,jci)*                     &
            drsc(i)*scs(jw)
        END DO

        ! Reset the O2 and O3 absorption cross sections.
        CALL cso2o3(ao2sr,ao3,tempc(jt),tspo2(jt,jc),                          &
          jw,jpwav)

        deptha(j,jc) = deptha(j,jc) + tablenoa(jt,j,jci)*                      &
          (do3c(jt)*ao3(jw) + do2c(jt)*(ao2(jw)+ao2sr(jw)))
        depths(j,jc) = depths(j,jc) + tablenoa(jt,j,jci)*                      &
          drsc(jt)*scs(jw)

        ! Now go from JT+1 all the way to the top of the atmosphere.
        DO i = jt + 1, jplev

          ! Reset the O2 and O3 absorption cross sections.
          CALL cso2o3(ao2sr,ao3,tempc(i),tspo2(i,jc),                          &
            jw,jpwav)

          deptha(j,jc) = deptha(j,jc) + tablenob(i,j,jci)*                     &
            (do3c(i)*ao3(jw) + do2c(i)*(ao2(jw)+ao2sr(jw)))
          depths(j,jc) = depths(j,jc) + tablenob(i,j,jci)*                     &
            drsc(i)*scs(jw)
        END DO
      END IF
    END DO
  END DO


  ! Set up the vertical optical depths etc.
  DO j = 1, jplev

    teac(j) = deptha(j,1)
    tauc(j) = depths(j,1)

  END DO
  teac(jplevp1) = 0.0
  tauc(jplevp1) = 0.0


  ! Assign the total vertical optical depth at the ground;

  tauv = 0.0
  teav = 0.0

  DO j = 1, jplev

    ! Reset O2 and O3 absorption cross sections.
    CALL cso2o3(ao2sr,ao3,tempc(j),tspo2(j,1),                                 &
      jw,jpwav)

    ! Scattering.
    tauv = tauv + dalt(j)*1.0e5* drsc(j)*scs(jw)

    ! Absorption.
    teav = teav + dalt(j)*1.0e5*(do3c(j)*ao3(jw) +                             &
      do2c(j)*(ao2(jw)+ao2sr(jw)))

  END DO


  tau(1) = tauv
  tea(1) = teav
  tau(jplevp1) = 0.0
  tea(jplevp1) = 0.0

  DO j = jplev, 2, -1
    !           Reset O2 and O3 absorption cross sections.
    CALL cso2o3(ao2sr,ao3,tempc(j),tspo2(j,1),                                 &
      jw,jpwav)

    tea(j) = tea(j+1) + dalt(j)*1.0e5*(do3c(j)*ao3(jw) +                       &
      do2c(j)*(ao2(jw)+ao2sr(jw)))
    tau(j) = tau(j+1) + dalt(j)*1.0e5* drsc(j)*scs(jw)
  END DO

  !        Calculate d tau / d delta tau

  DO i = 1, jplev

    !           Reset O2 and O3 absorption cross sections.
    CALL cso2o3(ao2sr,ao3,tempc(i),tspo2(i,1),                                 &
      jw,jpwav)

    jac(i) = 1.0/(1.0 + (do3c(i)*ao3(jw) + do2c(i)*                            &
      (ao2(jw)+ao2sr(jw)))/(drsc(i)*scs(jw)))

  END DO


  !        Calculate the initial scattering rate.

  !        Zenith angle loop.
  DO jc = 1, jpchi

    !          Ground reflected component.
    IF (jc  <=  jpchin) THEN
      arg = deptha(1,jc) + depths(1,jc)
      twoamu = 2.0*albedo*COS(tabang(jc))*EXP(-arg)
    ELSE
      twoamu = 0.0
    END IF

    !          Level loop.
    DO jk = 1, jplev

      arg = deptha(jk,jc) + depths(jk,jc)
      tabs0(jk,jc,jw) = EXP(-arg)

      IF (albedo > 0.0) THEN
        teak = teac(jk)
        tauk = tauc(jk)
        tabs0(jk,jc,jw) = tabs0(jk,jc,jw) + twoamu*ei2(                        &
          MAX(0.0,tauv - tauk + teav - teak))
      END IF

      !          End of level loop.
    END DO

    !        End of zenith angle loop.
  END DO


  IF (lscat) THEN

    !           Set up identity matrix
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(jn,jk) SHARED(delta)
!$OMP DO SCHEDULE(STATIC)
    DO jn = 1, jplev
      DO jk = 1, jplev
        delta(jn,jk)=0.0
      END DO
    END DO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
    DO jk = 1, jplev
      delta(jk,jk)=1.0
    END DO
!$OMP END DO
!$OMP END PARALLEL

        !           Two nested level loops.
!$OMP PARALLEL DO SCHEDULE(STATIC)  DEFAULT(NONE)                              &
!$OMP PRIVATE(jn,tean,taun,teanp,taunp,jk,teak,tauk)                           &
!$OMP SHARED(tea,tau,teac,tauc,a,jac,albedo,b,delta,tauv,teav)
    DO jn = 1, jplev

      tean = tea(jn)
      taun = tau(jn)
      teanp = tea(jn+1)
      taunp = tau(jn+1)

      DO jk = 1, jplev

        teak = teac(jk)
        tauk = tauc(jk)

        IF (jn /= jk) THEN
          a(jn,jk) = 0.5*ABS(                                                  &
            ei2(ABS(tauk - taun)  + ABS(teak - tean))                          &
            - ei2(ABS(tauk - taunp) + ABS(teak - teanp)))*jac(jn)
        ELSE
          a(jn,jk) = 0.5*(2.0 -                                                &
            ei2(ABS(tauk - taun)  + ABS(teak - tean))                          &
            - ei2(ABS(tauk - taunp) + ABS(teak - teanp)))*jac(jn)
        END IF

        IF (albedo > 0.0) a(jn,jk) = a(jn,jk) + albedo*                        &
          ei2(MAX(0.0,tauv - tauk  + teav - teak))*                            &
          ABS(ei3(MAX(0.0,tauv - taun  + teav - tean))                         &
          - ei3(MAX(0.0,tauv - taunp + teav - teanp)))*                        &
          jac(jn)

        !               Set up B matrix
        b(jn,jk) = delta(jn,jk) - a(jn,jk)

        !             End of nested level loops.
      END DO
    END DO
!$OMP END PARALLEL DO


        !           Invert the matrix B.
    CALL invert(b,binv,bcopy,indx,jplev)


    !           For greater than 90 degrees make sure where there is no direct
    !           sun that the initial scattering rate is zero.
!$OMP PARALLEL  DEFAULT(NONE)                                                  &
!$OMP PRIVATE(jc,jn,jk,jci,j,jt) SHARED(tabs,tabs0,binv,jw,ijt)
!$OMP DO SCHEDULE(STATIC)
    DO jc = jpchin + 1, jpchi

      !              A zenith angle index.
      jci = jc - jpchin

      DO j = 1, jplev

        !                 Recall what JT was.
        jt = ijt(jci,j)
        IF (jt  <=  0) tabs0(j,jc,jw) = 0.0
      END DO
    END DO
!$OMP END DO


        !           Zenith angle loop.
!$OMP DO SCHEDULE(STATIC)
    DO jc = 1, jpchi

      DO jn = 1, jplev

        !                 Initialise S.
        tabs(jn,jc,jw) = 0.0

        !                 Sum: S(JN)=S0(JK)*BINV(JK,JN)
        DO jk = 1, jplev
          tabs(jn,jc,jw) = tabs(jn,jc,jw) + tabs0(jk,jc,jw)*                   &
            binv(jk,jn)
        END DO
      END DO

      !           End of zenith angle loop.
    END DO
!$OMP END DO


        !        For greater than 90 degrees make sure where there is no direct
        !        sun that the scattering rate is zero.
!$OMP DO SCHEDULE(STATIC)
    DO jc = jpchin + 1, jpchi

      !              A zenith angle index.
      jci = jc - jpchin

      DO j = 1, jplev

        !                 Recall what JT was.
        jt = ijt(jci,j)
        IF (jt  <=  0) tabs(j,jc,jw) = 0.0
      END DO
    END DO
!$OMP END DO


        !           Check values are not negative.
!$OMP DO SCHEDULE(STATIC)
    DO jc = 1, jpchi
      DO j = 1, jplev
        tabs(j,jc,jw) = MAX(0.0,tabs(j,jc,jw))
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL


        !        Else if no scattering.
  ELSE


    !           Zenith angle loop.
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                               &
!$OMP PRIVATE(jc,jn) SHARED(tabs,jw,tabs0)
    DO jc = 1, jpchi

      DO jn = 1, jplev

        !                 Just initial scattering rate.
        tabs(jn,jc,jw) = tabs0(jn,jc,jw)

      END DO

    END DO
!$OMP END PARALLEL DO

  END IF


  !     End of wavelength loop.
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE settab
END MODULE settab_mod
