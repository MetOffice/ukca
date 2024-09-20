! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE calcjs_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE


! Description:
!     Look up photolysis rates from  j tables

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALCJS_MOD'

CONTAINS
SUBROUTINE calcjs(ihmin,ihmax,f,ppa,t,o3col,lnt)

USE ukca_um_dissoc_mod, ONLY: ajhno3, ajpna,  ajh2o2, aj2a,  aj2b, aj3,        &
       aj3a, ajcnita,ajcnitb,ajbrno3,ajbrcl, ajoclo, ajcl2o2,                  &
       ajhocl, ajno,  ajno2,  ajn2o5, ajno31, ajno32, ajbro,                   &
       ajhcl,  ajn2o,  ajhobr, ajf11,  ajf12,  ajh2o,  ajccl4,                 &
       ajf113, ajf22,  ajch3cl, ajc2oa, ajc2ob, ajmhp,  ajch3br,               &
       ajmcfm,ajch4,  ajf12b1,ajf13b1, ajcof2, ajcofcl, ajco2,                 &
       ajcos,  ajhono, ajmena, ajchbr3, ajdbrm, ajcs2, ajh2so4,                &
       ajso3
USE photol_constants_mod, ONLY: pi => const_pi

USE ukca_tbjs_mod,  ONLY:    tabj2a   , tabj2b   , tabj3    ,                  &
       tabj3a   , tabjno   , tabjno31 , tabjno32 , tabjno2  ,                  &
       tabjn2o5 , tabjhno3 , tabjh2o  , tabjh2o2 , tabjf11  ,                  &
       tabjf12  , tabjf22  , tabjf113 , tabjch3cl, tabjccl4 ,                  &
       tabjcnita, tabjcnitb, tabjhcl  , tabjhocl , tabjpna  ,                  &
       tabjcl2o2, tabjn2o  , tabjbrcl , tabjbrno3, tabjhobr ,                  &
       tabjbro  , tabjoclo , tabjc2oa , tabjc2ob , tabjmhp  ,                  &
       tabjmcfm , tabjch3br, tabjf12b1, tabjf13b1, tabjcof2 ,                  &
       tabjcofcl, tabjch4  , tabjcos  , tabjco2  ,                             &
       tabjhono , tabjmena , tabjchbr3, tabjdbrm , tabjcs2  ,                  &
       tabjh2so4, tabjso3  , tabpres,   tabang,    tabt     ,                  &
       tabo3, sao3c
USE ukca_parpho_mod, ONLY: jplev, jpchi, jps90, jpchin,                        &
                           jptem, jpo3p, jps90,                                &
                           szamax, tmin, tmax, o3min, o3max
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: ihmin        ! first point
INTEGER, INTENT(IN) :: ihmax        ! last point
INTEGER, INTENT(IN) :: lnt          ! No of points

REAL, INTENT(IN) :: f(lnt)          ! Cos of zenith angle
REAL, INTENT(IN) :: t(lnt)          ! Temperature
REAL, INTENT(IN) :: o3col(lnt)      ! Ozone column
REAL, INTENT(IN) :: ppa(lnt)        ! Pressure

! Local variables

INTEGER :: ih

REAL :: csup
REAL :: dlnp
REAL :: tdif
REAL :: o3dif
REAL :: radint
REAL :: cosint

INTEGER :: jx(lnt)
INTEGER :: jxp1(lnt)
INTEGER :: jy(lnt)
INTEGER :: jyp1(lnt)
INTEGER :: jz(lnt)
INTEGER :: jzp1(lnt)
INTEGER :: jo(lnt)
INTEGER :: jop1(lnt)

! Scalar copies of individual elements above, to reduce compile times:
INTEGER :: jx_0
INTEGER :: jxp1_0
INTEGER :: jy_0
INTEGER :: jyp1_0
INTEGER :: jz_0
INTEGER :: jzp1_0
INTEGER :: jo_0
INTEGER :: jop1_0

REAL :: angle(lnt)
REAL :: o3fac(lnt)
REAL :: p(lnt)
REAL :: tu(lnt)
REAL :: zo (lnt)
REAL :: zx (lnt)
REAL :: zy (lnt)
REAL :: zz (lnt)
REAL :: co (lnt)
REAL :: cx (lnt)
REAL :: cy (lnt)
REAL :: cz (lnt)
REAL :: w1 (lnt)
REAL :: w2 (lnt)
REAL :: w3 (lnt)
REAL :: w4 (lnt)
REAL :: w5 (lnt)
REAL :: w6 (lnt)
REAL :: w7 (lnt)
REAL :: w8 (lnt)
REAL :: w9 (lnt)
REAL :: w10(lnt)
REAL :: w11(lnt)
REAL :: w12(lnt)
REAL :: w13(lnt)
REAL :: w14(lnt)
REAL :: w15(lnt)
REAL :: w16(lnt)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALCJS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Spacing of angles in lookup table
cosint=1.0/REAL(jpchi-1-jps90)
radint=(szamax - 90.0)*pi/(180.0*REAL(jps90))

!     Model pressure levels in hPa. SZA in radians. T within range.
DO ih=ihmin,ihmax
  p (ih)=ppa(ih)*0.01
  p (ih)=MIN(p (ih), tabpres(    1))
  p (ih)=MAX(p (ih), tabpres(jplev))
  tu(ih)=t(ih)
  tu(ih)=MIN(tu(ih), tmax)
  tu(ih)=MAX(tu(ih), tmin)
  angle(ih)=ACOS(f(ih))
END DO

!     Find the location of the pressure levels,
!     zenith angle, temperature and O3 in the photolysis look up tables.

!     i)  X. Pressure levels equally spaced in log(p)
!     ii) Y. Zenith angle
dlnp=(LOG(tabpres(1))-LOG(tabpres(jplev)))/REAL(jplev-1)

DO ih = ihmin, ihmax
  jx(ih)=INT((LOG(tabpres(1))-LOG(p(ih)))/dlnp) + 1

  IF (angle(ih) < tabang(jpchin)) THEN
    jy(ih)=INT((1.0 - f(ih))/cosint) + 1
  ELSE
    jy(ih)=jpchin + INT((angle(ih)-0.5*pi)/radint)
  END IF
END DO

!     iii) Z  temperature evenly spaced
IF (jptem == 1) THEN
  DO ih = ihmin, ihmax
    jz(ih)=1
  END DO
ELSE
  tdif=(tmax-tmin)/REAL(jptem-1)
  DO ih = ihmin, ihmax
    jz(ih)=INT((tu(ih)-tmin)/tdif) + 1
  END DO
END IF

!     iv) O3 factor
IF (jpo3p == 1) THEN
  DO ih = ihmin, ihmax
    jo(ih)=1
  END DO
ELSE
  !       Find O3 factor using tabulated O3 column
  o3dif=(o3max-o3min)/REAL(jpo3p-1)
  DO ih=ihmin,ihmax
    csup=(p(ih)            -tabpres(jx(ih)))/                                  &
      (tabpres(jx(ih)+1)-tabpres(jx(ih)))

    !         Ratio O3 to standard atmosphere O3 above model pressure.
    o3fac(ih)=o3col(ih)/                                                       &
      ((1.0-csup)*sao3c(jx(ih))+csup*sao3c(jx(ih)+1))

    !         Check that O3 column ratio is within range
    o3fac(ih)=MIN(o3fac(ih), o3max)
    o3fac(ih)=MAX(o3fac(ih), o3min)

    jo(ih)=INT((o3fac(ih)-o3min)/o3dif) + 1
  END DO
END IF

!     If out of range set pointers accordingly.
DO ih=ihmin,ihmax
  jx(ih) = MAX(jx(ih),      1)
  jx(ih) = MIN(jx(ih),jplev-1)
  jy(ih) = MAX(jy(ih),      1)
  jy(ih) = MIN(jy(ih),jpchi-1)
  jz(ih) = MAX(jz(ih),      1)
  jz(ih) = MIN(jz(ih),jptem-1)
  jo(ih) = MAX(jo(ih),      1)
  jo(ih) = MIN(jo(ih),jpo3p-1)

  jz(ih) = MAX(jz(ih),1)
  jo(ih) = MAX(jo(ih),1)

  jxp1(ih) =      jx(ih) + 1
  jyp1(ih) =      jy(ih) + 1
  jzp1(ih) = MIN(jz(ih) + 1, jptem)
  jop1(ih) = MIN(jo(ih) + 1, jpo3p)
END DO

!     Find points used in interpolation
DO ih=ihmin,ihmax
  !        i) X: pressure
  zx(ih)=(LOG(p      (ih))       - LOG(tabpres(jx(ih))))/                      &
    (LOG(tabpres(jxp1(ih))) - LOG(tabpres(jx(ih))))

  !        ii) Y: zenith angle
  zy(ih)=(angle(ih)         - tabang (jy(ih)))/                                &
    (tabang (jyp1(ih)) - tabang (jy(ih)))

  !        iii) Z: temperature
  IF (jptem == 1) THEN
    zz(ih)=1.0
  ELSE
    zz(ih)=(tu   (ih)       - tabt (jz(ih)))/                                  &
      (tabt (jzp1(ih)) - tabt (jz(ih)))
  END IF

  !        iv) O: O3 profile
  IF (jpo3p == 1) THEN
    zo(ih)=1.0
  ELSE
    zo(ih)=(o3fac(ih)       - tabo3(jo(ih)))/                                  &
      (tabo3(jop1(ih)) - tabo3(jo(ih)))
  END IF

  cx(ih)=1.0 - zx(ih)
  cy(ih)=1.0 - zy(ih)
  cz(ih)=1.0 - zz(ih)
  co(ih)=1.0 - zo(ih)
END DO

!     Calculate weights for interpolation
DO ih=ihmin,ihmax
  w1 (ih)=cx(ih)*cy(ih)*cz(ih)*co(ih)
  w2 (ih)=zx(ih)*cy(ih)*cz(ih)*co(ih)
  w3 (ih)=zx(ih)*zy(ih)*cz(ih)*co(ih)
  w4 (ih)=cx(ih)*zy(ih)*cz(ih)*co(ih)
  w5 (ih)=cx(ih)*cy(ih)*zz(ih)*co(ih)
  w6 (ih)=zx(ih)*cy(ih)*zz(ih)*co(ih)
  w7 (ih)=zx(ih)*zy(ih)*zz(ih)*co(ih)
  w8 (ih)=cx(ih)*zy(ih)*zz(ih)*co(ih)
  w9 (ih)=cx(ih)*cy(ih)*cz(ih)*zo(ih)
  w10(ih)=zx(ih)*cy(ih)*cz(ih)*zo(ih)
  w11(ih)=zx(ih)*zy(ih)*cz(ih)*zo(ih)
  w12(ih)=cx(ih)*zy(ih)*cz(ih)*zo(ih)
  w13(ih)=cx(ih)*cy(ih)*zz(ih)*zo(ih)
  w14(ih)=zx(ih)*cy(ih)*zz(ih)*zo(ih)
  w15(ih)=zx(ih)*zy(ih)*zz(ih)*zo(ih)
  w16(ih)=cx(ih)*zy(ih)*zz(ih)*zo(ih)
END DO

!     Look up photolysis rates. 38 rates.
DO ih=ihmin,ihmax

  !     If angle in range
  IF (angle(ih) <= tabang(jpchi)) THEN
    jx_0   = jx(ih)
    jy_0   = jy(ih)
    jz_0   = jz(ih)
    jo_0   = jo(ih)
    jxp1_0 = jxp1(ih)
    jyp1_0 = jyp1(ih)
    jzp1_0 = jzp1(ih)
    jop1_0 = jop1(ih)

    !1
    aj2a  (ih)=                                                                &
      w1 (ih)*tabj2a   (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabj2a   (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabj2a   (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabj2a   (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabj2a   (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabj2a   (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabj2a   (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabj2a   (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabj2a   (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabj2a   (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabj2a   (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabj2a   (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabj2a   (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabj2a   (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabj2a   (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabj2a   (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !1B
    aj2b  (ih)=                                                                &
      w1 (ih)*tabj2b   (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabj2b   (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabj2b   (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabj2b   (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabj2b   (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabj2b   (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabj2b   (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabj2b   (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabj2b   (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabj2b   (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabj2b   (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabj2b   (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabj2b   (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabj2b   (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabj2b   (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabj2b   (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !2
    aj3   (ih)=                                                                &
      w1 (ih)*tabj3    (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabj3    (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabj3    (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabj3    (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabj3    (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabj3    (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabj3    (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabj3    (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabj3    (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabj3    (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabj3    (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabj3    (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabj3    (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabj3    (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabj3    (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabj3    (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !3
    aj3a  (ih)=                                                                &
      w1 (ih)*tabj3a   (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabj3a   (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabj3a   (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabj3a   (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabj3a   (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabj3a   (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabj3a   (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabj3a   (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabj3a   (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabj3a   (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabj3a   (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabj3a   (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabj3a   (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabj3a   (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabj3a   (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabj3a   (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !4
    ajno  (ih)=                                                                &
      w1 (ih)*tabjno   (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjno   (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjno   (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjno   (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjno   (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjno   (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjno   (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjno   (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjno   (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjno   (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjno   (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjno   (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjno   (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjno   (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjno   (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjno   (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !5
    ajno2 (ih)=                                                                &
      w1 (ih)*tabjno2  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjno2  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjno2  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjno2  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjno2  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjno2  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjno2  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjno2  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjno2  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjno2  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjno2  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjno2  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjno2  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjno2  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjno2  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjno2  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !6
    ajno31(ih)=                                                                &
      w1 (ih)*tabjno31 (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjno31 (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjno31 (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjno31 (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjno31 (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjno31 (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjno31 (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjno31 (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjno31 (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjno31 (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjno31 (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjno31 (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjno31 (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjno31 (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjno31 (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjno31 (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !7
    ajno32(ih)=                                                                &
      w1 (ih)*tabjno32 (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjno32 (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjno32 (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjno32 (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjno32 (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjno32 (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjno32 (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjno32 (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjno32 (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjno32 (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjno32 (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjno32 (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjno32 (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjno32 (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjno32 (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjno32 (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !8
    ajn2o (ih)=                                                                &
      w1 (ih)*tabjn2o  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjn2o  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjn2o  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjn2o  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjn2o  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjn2o  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjn2o  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjn2o  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjn2o  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjn2o  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjn2o  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjn2o  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjn2o  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjn2o  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjn2o  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjn2o  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !9
    ajn2o5(ih)=                                                                &
      w1 (ih)*tabjn2o5 (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjn2o5 (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjn2o5 (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjn2o5 (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjn2o5 (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjn2o5 (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjn2o5 (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjn2o5 (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjn2o5 (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjn2o5 (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjn2o5 (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjn2o5 (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjn2o5 (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjn2o5 (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjn2o5 (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjn2o5 (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !10
    ajhno3(ih)=                                                                &
      w1 (ih)*tabjhno3 (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjhno3 (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjhno3 (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjhno3 (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjhno3 (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjhno3 (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjhno3 (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjhno3 (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjhno3 (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjhno3 (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjhno3 (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjhno3 (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjhno3 (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjhno3 (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjhno3 (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjhno3 (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !11
    ajcnita(ih)=                                                               &
      w1 (ih)*tabjcnita(jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjcnita(jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjcnita(jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjcnita(jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjcnita(jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjcnita(jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjcnita(jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjcnita(jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjcnita(jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjcnita(jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjcnita(jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjcnita(jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjcnita(jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjcnita(jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjcnita(jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjcnita(jx_0  ,jyp1_0,jzp1_0,jop1_0)
    ajcnitb(ih)=                                                               &
      w1 (ih)*tabjcnitb(jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjcnitb(jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjcnitb(jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjcnitb(jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjcnitb(jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjcnitb(jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjcnitb(jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjcnitb(jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjcnitb(jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjcnitb(jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjcnitb(jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjcnitb(jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjcnitb(jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjcnitb(jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjcnitb(jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjcnitb(jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !12
    ajpna (ih)=                                                                &
      w1 (ih)*tabjpna  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjpna  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjpna  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjpna  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjpna  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjpna  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjpna  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjpna  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjpna  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjpna  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjpna  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjpna  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjpna  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjpna  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjpna  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjpna  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !13
    ajh2o2(ih)=                                                                &
      w1 (ih)*tabjh2o2 (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjh2o2 (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjh2o2 (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjh2o2 (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjh2o2 (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjh2o2 (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjh2o2 (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjh2o2 (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjh2o2 (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjh2o2 (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjh2o2 (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjh2o2 (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjh2o2 (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjh2o2 (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjh2o2 (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjh2o2 (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !14
    ajh2o (ih)=                                                                &
      w1 (ih)*tabjh2o  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjh2o  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjh2o  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjh2o  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjh2o  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjh2o  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjh2o  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjh2o  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjh2o  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjh2o  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjh2o  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjh2o  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjh2o  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjh2o  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjh2o  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjh2o  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !15
    ajhocl(ih)=                                                                &
      w1 (ih)*tabjhocl (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjhocl (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjhocl (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjhocl (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjhocl (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjhocl (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjhocl (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjhocl (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjhocl (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjhocl (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjhocl (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjhocl (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjhocl (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjhocl (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjhocl (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjhocl (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !16
    ajhcl (ih)=                                                                &
      w1 (ih)*tabjhcl  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjhcl  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjhcl  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjhcl  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjhcl  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjhcl  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjhcl  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjhcl  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjhcl  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjhcl  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjhcl  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjhcl  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjhcl  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjhcl  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjhcl  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjhcl  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !17
    ajcl2o2(ih)=                                                               &
      w1 (ih)*tabjcl2o2(jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjcl2o2(jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjcl2o2(jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjcl2o2(jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjcl2o2(jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjcl2o2(jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjcl2o2(jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjcl2o2(jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjcl2o2(jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjcl2o2(jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjcl2o2(jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjcl2o2(jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjcl2o2(jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjcl2o2(jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjcl2o2(jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjcl2o2(jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !18
    ajbro (ih)=                                                                &
      w1 (ih)*tabjbro  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjbro  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjbro  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjbro  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjbro  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjbro  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjbro  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjbro  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjbro  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjbro  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjbro  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjbro  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjbro  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjbro  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjbro  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjbro  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !19
    ajbrcl(ih)=                                                                &
      w1 (ih)*tabjbrcl (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjbrcl (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjbrcl (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjbrcl (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjbrcl (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjbrcl (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjbrcl (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjbrcl (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjbrcl (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjbrcl (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjbrcl (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjbrcl (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjbrcl (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjbrcl (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjbrcl (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjbrcl (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !20
    ajco2(ih)=                                                                 &
      w1 (ih)*tabjco2  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjco2  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjco2  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjco2  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjco2  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjco2  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjco2  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjco2  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjco2  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjco2  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjco2  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjco2  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjco2  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjco2  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjco2  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjco2  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    ! COS
    ajcos(ih)=                                                                 &
      w1 (ih)*tabjcos  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjcos  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjcos  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjcos  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjcos  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjcos  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjcos  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjcos  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjcos  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjcos  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjcos  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjcos  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjcos  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjcos  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjcos  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjcos  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    ! HONO
    ajhono(ih) =                                                               &
      w1 (ih)*tabjhono (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjhono (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjhono (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjhono (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjhono (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjhono (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjhono (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjhono (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjhono (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjhono (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjhono (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjhono (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjhono (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjhono (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjhono (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjhono (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    ! MeONO2
    ajmena(ih) =                                                               &
      w1 (ih)*tabjmena (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjmena (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjmena (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjmena (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjmena (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjmena (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjmena (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjmena (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjmena (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjmena (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjmena (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjmena (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjmena (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjmena (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjmena (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjmena (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    ! CHBr3
    ajchbr3(ih) =                                                              &
      w1 (ih)*tabjchbr3(jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjchbr3(jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjchbr3(jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjchbr3(jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjchbr3(jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjchbr3(jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjchbr3(jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjchbr3(jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjchbr3(jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjchbr3(jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjchbr3(jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjchbr3(jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjchbr3(jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjchbr3(jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjchbr3(jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjchbr3(jx_0  ,jyp1_0,jzp1_0,jop1_0)
    ! CH2Br2
    ajdbrm (ih) =                                                              &
      w1 (ih)*tabjdbrm (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjdbrm (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjdbrm (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjdbrm (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjdbrm (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjdbrm (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjdbrm (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjdbrm (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjdbrm (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjdbrm (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjdbrm (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjdbrm (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjdbrm (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjdbrm (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjdbrm (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjdbrm (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    ! CS2
    ajcs2 (ih) =                                                               &
      w1 (ih)*tabjcs2 (jx_0  ,jy_0  ,jz_0  ,jo_0  )                            &
      + w2 (ih)*tabjcs2 (jxp1_0,jy_0  ,jz_0  ,jo_0  )                          &
      + w3 (ih)*tabjcs2 (jxp1_0,jyp1_0,jz_0  ,jo_0  )                          &
      + w4 (ih)*tabjcs2 (jx_0  ,jyp1_0,jz_0  ,jo_0  )                          &
      + w5 (ih)*tabjcs2 (jx_0  ,jy_0  ,jzp1_0,jo_0  )                          &
      + w6 (ih)*tabjcs2 (jxp1_0,jy_0  ,jzp1_0,jo_0  )                          &
      + w7 (ih)*tabjcs2 (jxp1_0,jyp1_0,jzp1_0,jo_0  )                          &
      + w8 (ih)*tabjcs2 (jx_0  ,jyp1_0,jzp1_0,jo_0  )                          &
      + w9 (ih)*tabjcs2 (jx_0  ,jy_0  ,jz_0  ,jop1_0)                          &
      + w10(ih)*tabjcs2 (jxp1_0,jy_0  ,jz_0  ,jop1_0)                          &
      + w11(ih)*tabjcs2 (jxp1_0,jyp1_0,jz_0  ,jop1_0)                          &
      + w12(ih)*tabjcs2 (jx_0  ,jyp1_0,jz_0  ,jop1_0)                          &
      + w13(ih)*tabjcs2 (jx_0  ,jy_0  ,jzp1_0,jop1_0)                          &
      + w14(ih)*tabjcs2 (jxp1_0,jy_0  ,jzp1_0,jop1_0)                          &
      + w15(ih)*tabjcs2 (jxp1_0,jyp1_0,jzp1_0,jop1_0)                          &
      + w16(ih)*tabjcs2 (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    ! H2SO4
    ajh2so4 (ih) =                                                             &
      w1 (ih)*tabjh2so4 (jx_0  ,jy_0  ,jz_0  ,jo_0  )                          &
      + w2 (ih)*tabjh2so4 (jxp1_0,jy_0  ,jz_0  ,jo_0  )                        &
      + w3 (ih)*tabjh2so4 (jxp1_0,jyp1_0,jz_0  ,jo_0  )                        &
      + w4 (ih)*tabjh2so4 (jx_0  ,jyp1_0,jz_0  ,jo_0  )                        &
      + w5 (ih)*tabjh2so4 (jx_0  ,jy_0  ,jzp1_0,jo_0  )                        &
      + w6 (ih)*tabjh2so4 (jxp1_0,jy_0  ,jzp1_0,jo_0  )                        &
      + w7 (ih)*tabjh2so4 (jxp1_0,jyp1_0,jzp1_0,jo_0  )                        &
      + w8 (ih)*tabjh2so4 (jx_0  ,jyp1_0,jzp1_0,jo_0  )                        &
      + w9 (ih)*tabjh2so4 (jx_0  ,jy_0  ,jz_0  ,jop1_0)                        &
      + w10(ih)*tabjh2so4 (jxp1_0,jy_0  ,jz_0  ,jop1_0)                        &
      + w11(ih)*tabjh2so4 (jxp1_0,jyp1_0,jz_0  ,jop1_0)                        &
      + w12(ih)*tabjh2so4 (jx_0  ,jyp1_0,jz_0  ,jop1_0)                        &
      + w13(ih)*tabjh2so4 (jx_0  ,jy_0  ,jzp1_0,jop1_0)                        &
      + w14(ih)*tabjh2so4 (jxp1_0,jy_0  ,jzp1_0,jop1_0)                        &
      + w15(ih)*tabjh2so4 (jxp1_0,jyp1_0,jzp1_0,jop1_0)                        &
      + w16(ih)*tabjh2so4 (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    ! SO3
    ajso3 (ih) =                                                               &
      w1 (ih)*tabjso3 (jx_0  ,jy_0  ,jz_0  ,jo_0  )                            &
      + w2 (ih)*tabjso3 (jxp1_0,jy_0  ,jz_0  ,jo_0  )                          &
      + w3 (ih)*tabjso3 (jxp1_0,jyp1_0,jz_0  ,jo_0  )                          &
      + w4 (ih)*tabjso3 (jx_0  ,jyp1_0,jz_0  ,jo_0  )                          &
      + w5 (ih)*tabjso3 (jx_0  ,jy_0  ,jzp1_0,jo_0  )                          &
      + w6 (ih)*tabjso3 (jxp1_0,jy_0  ,jzp1_0,jo_0  )                          &
      + w7 (ih)*tabjso3 (jxp1_0,jyp1_0,jzp1_0,jo_0  )                          &
      + w8 (ih)*tabjso3 (jx_0  ,jyp1_0,jzp1_0,jo_0  )                          &
      + w9 (ih)*tabjso3 (jx_0  ,jy_0  ,jz_0  ,jop1_0)                          &
      + w10(ih)*tabjso3 (jxp1_0,jy_0  ,jz_0  ,jop1_0)                          &
      + w11(ih)*tabjso3 (jxp1_0,jyp1_0,jz_0  ,jop1_0)                          &
      + w12(ih)*tabjso3 (jx_0  ,jyp1_0,jz_0  ,jop1_0)                          &
      + w13(ih)*tabjso3 (jx_0  ,jy_0  ,jzp1_0,jop1_0)                          &
      + w14(ih)*tabjso3 (jxp1_0,jy_0  ,jzp1_0,jop1_0)                          &
      + w15(ih)*tabjso3 (jxp1_0,jyp1_0,jzp1_0,jop1_0)                          &
      + w16(ih)*tabjso3 (jx_0  ,jyp1_0,jzp1_0,jop1_0)

    !     Else darkness
  ELSE
    aj2a   (ih) = 0.0
    aj2b   (ih) = 0.0
    aj3    (ih) = 0.0
    aj3a   (ih) = 0.0
    ajno   (ih) = 0.0
    ajno2  (ih) = 0.0
    ajno31 (ih) = 0.0
    ajno32 (ih) = 0.0
    ajn2o  (ih) = 0.0
    ajn2o5 (ih) = 0.0
    ajhno3 (ih) = 0.0
    ajcnita(ih) = 0.0
    ajcnitb(ih) = 0.0
    ajpna  (ih) = 0.0
    ajh2o2 (ih) = 0.0
    ajh2o  (ih) = 0.0
    ajhocl (ih) = 0.0
    ajhcl  (ih) = 0.0
    ajcl2o2(ih) = 0.0
    ajbro  (ih) = 0.0
    ajbrcl (ih) = 0.0
    ajco2  (ih) = 0.0
    ajcos  (ih) = 0.0
    ajhono (ih) = 0.0
    ajmena (ih) = 0.0
    ajchbr3(ih) = 0.0
    ajdbrm (ih) = 0.0
    ajcs2  (ih) = 0.0
    ajh2so4(ih) = 0.0
    ajso3  (ih) = 0.0
  END IF

END DO

DO ih=ihmin,ihmax

  !     If angle in range
  IF (angle(ih) <= tabang(jpchi)) THEN
    jx_0   = jx(ih)
    jy_0   = jy(ih)
    jz_0   = jz(ih)
    jo_0   = jo(ih)
    jxp1_0 = jxp1(ih)
    jyp1_0 = jyp1(ih)
    jzp1_0 = jzp1(ih)
    jop1_0 = jop1(ih)

    !20
    ajbrno3(ih)=                                                               &
      w1 (ih)*tabjbrno3(jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjbrno3(jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjbrno3(jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjbrno3(jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjbrno3(jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjbrno3(jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjbrno3(jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjbrno3(jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjbrno3(jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjbrno3(jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjbrno3(jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjbrno3(jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjbrno3(jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjbrno3(jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjbrno3(jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjbrno3(jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !21
    ajhobr(ih)=                                                                &
      w1 (ih)*tabjhobr (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjhobr (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjhobr (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjhobr (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjhobr (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjhobr (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjhobr (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjhobr (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjhobr (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjhobr (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjhobr (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjhobr (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjhobr (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjhobr (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjhobr (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjhobr (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !22
    ajoclo(ih)=                                                                &
      w1 (ih)*tabjoclo (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjoclo (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjoclo (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjoclo (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjoclo (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjoclo (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjoclo (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjoclo (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjoclo (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjoclo (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjoclo (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjoclo (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjoclo (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjoclo (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjoclo (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjoclo (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !23
    ajc2oa(ih)=                                                                &
      w1 (ih)*tabjc2oa (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjc2oa (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjc2oa (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjc2oa (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjc2oa (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjc2oa (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjc2oa (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjc2oa (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjc2oa (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjc2oa (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjc2oa (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjc2oa (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjc2oa (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjc2oa (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjc2oa (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjc2oa (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !24
    ajc2ob(ih)=                                                                &
      w1 (ih)*tabjc2ob (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjc2ob (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjc2ob (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjc2ob (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjc2ob (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjc2ob (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjc2ob (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjc2ob (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjc2ob (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjc2ob (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjc2ob (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjc2ob (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjc2ob (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjc2ob (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjc2ob (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjc2ob (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !25
    ajmhp (ih)=                                                                &
      w1 (ih)*tabjmhp  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjmhp  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjmhp  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjmhp  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjmhp  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjmhp  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjmhp  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjmhp  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjmhp  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjmhp  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjmhp  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjmhp  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjmhp  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjmhp  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjmhp  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjmhp  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !26
    ajch3cl(ih)=                                                               &
      w1 (ih)*tabjch3cl(jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjch3cl(jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjch3cl(jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjch3cl(jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjch3cl(jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjch3cl(jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjch3cl(jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjch3cl(jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjch3cl(jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjch3cl(jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjch3cl(jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjch3cl(jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjch3cl(jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjch3cl(jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjch3cl(jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjch3cl(jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !27
    ajmcfm(ih)=                                                                &
      w1 (ih)*tabjmcfm (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjmcfm (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjmcfm (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjmcfm (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjmcfm (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjmcfm (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjmcfm (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjmcfm (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjmcfm (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjmcfm (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjmcfm (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjmcfm (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjmcfm (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjmcfm (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjmcfm (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjmcfm (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !28
    ajf11 (ih)=                                                                &
      w1 (ih)*tabjf11  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjf11  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjf11  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjf11  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjf11  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjf11  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjf11  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjf11  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjf11  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjf11  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjf11  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjf11  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjf11  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjf11  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjf11  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjf11  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !29
    ajf12 (ih)=                                                                &
      w1 (ih)*tabjf12  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjf12  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjf12  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjf12  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjf12  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjf12  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjf12  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjf12  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjf12  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjf12  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjf12  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjf12  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjf12  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjf12  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjf12  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjf12  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !30
    ajf22 (ih)=                                                                &
      w1 (ih)*tabjf22  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjf22  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjf22  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjf22  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjf22  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjf22  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjf22  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjf22  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjf22  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjf22  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjf22  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjf22  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjf22  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjf22  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjf22  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjf22  (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !31
    ajf113(ih)=                                                                &
      w1 (ih)*tabjf113 (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjf113 (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjf113 (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjf113 (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjf113 (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjf113 (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjf113 (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjf113 (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjf113 (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjf113 (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjf113 (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjf113 (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjf113 (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjf113 (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjf113 (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjf113 (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !32
    ajccl4(ih)=                                                                &
      w1 (ih)*tabjccl4 (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjccl4 (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjccl4 (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjccl4 (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjccl4 (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjccl4 (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjccl4 (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjccl4 (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjccl4 (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjccl4 (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjccl4 (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjccl4 (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjccl4 (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjccl4 (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjccl4 (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjccl4 (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !33
    ajch3br(ih)=                                                               &
      w1 (ih)*tabjch3br(jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjch3br(jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjch3br(jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjch3br(jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjch3br(jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjch3br(jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjch3br(jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjch3br(jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjch3br(jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjch3br(jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjch3br(jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjch3br(jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjch3br(jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjch3br(jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjch3br(jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjch3br(jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !34
    ajf12b1(ih)=                                                               &
      w1 (ih)*tabjf12b1(jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjf12b1(jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjf12b1(jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjf12b1(jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjf12b1(jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjf12b1(jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjf12b1(jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjf12b1(jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjf12b1(jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjf12b1(jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjf12b1(jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjf12b1(jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjf12b1(jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjf12b1(jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjf12b1(jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjf12b1(jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !35
    ajf13b1(ih)=                                                               &
      w1 (ih)*tabjf13b1(jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjf13b1(jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjf13b1(jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjf13b1(jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjf13b1(jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjf13b1(jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjf13b1(jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjf13b1(jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjf13b1(jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjf13b1(jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjf13b1(jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjf13b1(jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjf13b1(jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjf13b1(jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjf13b1(jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjf13b1(jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !36
    ajcof2(ih)=                                                                &
      w1 (ih)*tabjcof2 (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjcof2 (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjcof2 (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjcof2 (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjcof2 (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjcof2 (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjcof2 (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjcof2 (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjcof2 (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjcof2 (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjcof2 (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjcof2 (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjcof2 (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjcof2 (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjcof2 (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjcof2 (jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !37
    ajcofcl(ih)=                                                               &
      w1 (ih)*tabjcofcl(jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjcofcl(jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjcofcl(jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjcofcl(jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjcofcl(jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjcofcl(jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjcofcl(jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjcofcl(jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjcofcl(jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjcofcl(jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjcofcl(jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjcofcl(jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjcofcl(jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjcofcl(jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjcofcl(jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjcofcl(jx_0  ,jyp1_0,jzp1_0,jop1_0)
    !38
    ajch4 (ih)=                                                                &
      w1 (ih)*tabjch4  (jx_0  ,jy_0  ,jz_0  ,jo_0  )                           &
      + w2 (ih)*tabjch4  (jxp1_0,jy_0  ,jz_0  ,jo_0  )                         &
      + w3 (ih)*tabjch4  (jxp1_0,jyp1_0,jz_0  ,jo_0  )                         &
      + w4 (ih)*tabjch4  (jx_0  ,jyp1_0,jz_0  ,jo_0  )                         &
      + w5 (ih)*tabjch4  (jx_0  ,jy_0  ,jzp1_0,jo_0  )                         &
      + w6 (ih)*tabjch4  (jxp1_0,jy_0  ,jzp1_0,jo_0  )                         &
      + w7 (ih)*tabjch4  (jxp1_0,jyp1_0,jzp1_0,jo_0  )                         &
      + w8 (ih)*tabjch4  (jx_0  ,jyp1_0,jzp1_0,jo_0  )                         &
      + w9 (ih)*tabjch4  (jx_0  ,jy_0  ,jz_0  ,jop1_0)                         &
      + w10(ih)*tabjch4  (jxp1_0,jy_0  ,jz_0  ,jop1_0)                         &
      + w11(ih)*tabjch4  (jxp1_0,jyp1_0,jz_0  ,jop1_0)                         &
      + w12(ih)*tabjch4  (jx_0  ,jyp1_0,jz_0  ,jop1_0)                         &
      + w13(ih)*tabjch4  (jx_0  ,jy_0  ,jzp1_0,jop1_0)                         &
      + w14(ih)*tabjch4  (jxp1_0,jy_0  ,jzp1_0,jop1_0)                         &
      + w15(ih)*tabjch4  (jxp1_0,jyp1_0,jzp1_0,jop1_0)                         &
      + w16(ih)*tabjch4  (jx_0  ,jyp1_0,jzp1_0,jop1_0)

    !     Else darkness
  ELSE
    ajbrno3(ih) = 0.0
    ajhobr (ih) = 0.0
    ajoclo (ih) = 0.0
    ajc2oa (ih) = 0.0
    ajc2ob (ih) = 0.0
    ajmhp  (ih) = 0.0
    ajch3cl(ih) = 0.0
    ajmcfm (ih) = 0.0
    ajf11  (ih) = 0.0
    ajf12  (ih) = 0.0
    ajf22  (ih) = 0.0
    ajf113 (ih) = 0.0
    ajccl4 (ih) = 0.0
    ajch3br(ih) = 0.0
    ajf12b1(ih) = 0.0
    ajf13b1(ih) = 0.0
    ajcof2 (ih) = 0.0
    ajcofcl(ih) = 0.0
    ajch4  (ih) = 0.0
  END IF

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE calcjs
END MODULE calcjs_mod
