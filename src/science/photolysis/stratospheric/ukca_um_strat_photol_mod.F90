! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************

!  Description:
!   Module containing subroutine strat_photol.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM

!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################

MODULE ukca_um_strat_photol_mod



USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus,                         &
                      PrStatus_Normal, PrStatus_Oper
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_UM_STRAT_PHOTOL_MOD'

CONTAINS

!---------------------------------------------------------------------------
! Subroutine STRAT_PHOTOL
!------------------------------------------------------------------------

! This routine computes stratospheric photolysis rates and merges the
! rates, where necessary, with the tropospheric rates. This is done for
! one level at a time. The stratospheric photolysis routines are taken
! from SLIMCAT.


SUBROUTINE strat_photol(                                                       &
  row_length, rows,                                                            &
  jppj, spj,                                                                   &
  pressure,                                                                    &
  temp,                                                                        &
  ozonecol,                                                                    &
  cos_zenith_angle,                                                            &
  current_time,                                                                &
  photrates)

! SLIMCAT stratospheric photolysis routines
USE calcjs_mod, ONLY: calcjs
USE inijtab_mod, ONLY: inijtab

USE photol_config_specification_mod, ONLY: photol_config, i_scheme_fastjx,     &
                               fjx_mode_2Donly

! Holds photolysis rates
USE ukca_um_dissoc_mod,   ONLY: aj2a, aj2b, aj3, aj3a, ajbrcl, ajbrno3, ajbro, &
                                ajc2oa, ajc2ob, ajccl4, ajch3br, ajch3cl,      &
                                ajch4, ajchbr3, ajcl2o2, ajco2, ajcof2,        &
                                ajcofcl,ajcos, ajcnita, ajcnitb, ajcs2, ajdbrm,&
                                ajf11, ajf113, ajf12, ajf12b1, ajf13b1,        &
                                ajf22, ajh2o, ajh2o2, ajh2so4, ajhcl, ajhobr,  &
                                ajhocl, ajhno3, ajmcfm, ajmhp, ajn2o,          &
                                ajn2o5, ajno, ajno2, ajno31, ajno32, ajoclo,   &
                                ajpna, ajso3

USE um_parcore, ONLY: mype

IMPLICIT NONE


! Subroutine interface

! Model dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows

! Photolysis reaction data
INTEGER, INTENT(IN) :: jppj
CHARACTER(LEN=10), POINTER, INTENT(IN) :: spj(:,:)

REAL, INTENT(IN) :: pressure(row_length, rows)
REAL, INTENT(IN) :: temp(row_length, rows)
REAL, INTENT(IN) :: ozonecol(row_length, rows)
REAL, INTENT(IN) :: cos_zenith_angle(row_length, rows)
INTEGER, INTENT(IN)  :: current_time(7)
REAL, INTENT(IN OUT) :: photrates(row_length, rows, jppj)

! local variables
INTEGER :: theta_field_size
INTEGER :: i
INTEGER :: j
INTEGER :: l
INTEGER :: k

REAL :: frac
LOGICAL, SAVE :: firstcall = .TRUE.
LOGICAL, SAVE :: tables_filled = .FALSE.
INTEGER, SAVE :: current_month = 0

INTEGER, SAVE :: ih2o2=0
INTEGER, SAVE :: ihchoa=0
INTEGER, SAVE :: ihchob=0
INTEGER, SAVE :: iho2no2a=0
INTEGER, SAVE :: ihno3=0
INTEGER, SAVE :: imeooh=0
INTEGER, SAVE :: in2o5=0
INTEGER, SAVE :: ino2=0
INTEGER, SAVE :: ino3a=0
INTEGER, SAVE :: ino3b=0
INTEGER, SAVE :: io2a=0
INTEGER, SAVE :: io2b=0
INTEGER, SAVE :: io3a=0
INTEGER, SAVE :: io3b=0
INTEGER, SAVE :: io3sa=0
INTEGER, SAVE :: io3sb=0
INTEGER, SAVE :: ich4=0
INTEGER, SAVE :: ih2o=0
INTEGER, SAVE :: ino=0
INTEGER, SAVE :: in2o=0
INTEGER, SAVE :: if11=0
INTEGER, SAVE :: if12=0
INTEGER, SAVE :: iclono2a=0
INTEGER, SAVE :: iclono2b=0
INTEGER, SAVE :: ihcl=0
INTEGER, SAVE :: ihocl=0
INTEGER, SAVE :: ioclo=0
INTEGER, SAVE :: icl2o2=0
INTEGER, SAVE :: ibro=0
INTEGER, SAVE :: ihobr=0
INTEGER, SAVE :: ibrono2a=0
INTEGER, SAVE :: ibrcl=0
INTEGER, SAVE :: imebr=0
INTEGER, SAVE :: iccl4=0
INTEGER, SAVE :: if113=0
INTEGER, SAVE :: imecl=0
INTEGER, SAVE :: imcfm=0
INTEGER, SAVE :: if22=0
INTEGER, SAVE :: ih1211=0
INTEGER, SAVE :: ih1301=0
INTEGER, SAVE :: icof2=0
INTEGER, SAVE :: icofcl=0
INTEGER, SAVE :: ico2=0
INTEGER, SAVE :: ibrono2b=0
INTEGER, SAVE :: iho2no2b=0
INTEGER, SAVE :: icos=0
INTEGER, SAVE :: idbrm=0
INTEGER, SAVE :: ichbr3=0
INTEGER, SAVE :: ics2=0
INTEGER, SAVE :: ih2so4=0
INTEGER, SAVE :: iso3=0

! Local time variable
INTEGER :: i_month

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STRAT_PHOTOL'

! upon first entry initialize positions of photolysis reactions
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set local time variable
i_month = current_time(2)

IF (firstcall) THEN
  DO k=1,jppj
    ! Merge tropospheric and stratospheric photolysis rates
    SELECT CASE (spj(k,1))

      ! H2O2   + h nu -> 2 OH
    CASE ('H2O2      ')
      ih2o2 = k

      ! consider branching in case of HCHO
    CASE ('HCHO      ')
      ! HCHO   + h nu -> H2 + CO
      IF ((spj(k,3) == 'H2        ') .OR.                                      &
        (spj(k,4) == 'H2        ')) THEN
        ihchob = k
      ELSE
        ! HCHO   + h nu -> H + CHO -> CO + 2 HO2
        ihchoa = k
      END IF

      ! HO2NO2 + h nu -> NO2 + HO2
    CASE ('HO2NO2    ')
      IF ((spj(k,3) == 'HO2       ')                                           &
        .OR. (spj(k,4) == 'HO2       ')) THEN
        ! HO2NO2 + h nu -> NO2 + HO2
        iho2no2a = k
      ELSE
        ! HO2NO2 + h nu -> NO3 + OH
        iho2no2b = k
      END IF

      ! HNO3   + h nu -> NO2 + OH
    CASE ('HONO2     ')
      ihno3 = k

      ! MeOOH  + h nu -> MeO + OH -> HCHO + HO2 + OH
    CASE ('MeOOH     ')
      imeooh = k

      ! N2O5   + h nu -> NO2 + NO3
    CASE ('N2O5      ')
      in2o5 = k

      ! NO2    + h nu -> NO  + O
    CASE ('NO2       ')
      ino2 = k

      ! consider branching for NO3
    CASE ('NO3       ')
      IF ((spj(k,3) == 'O(3P)     ')                                           &
        .OR. (spj(k,4) == 'O(3P)     ')) THEN
        ! NO3    + h nu -> NO2 + O
        ino3b = k
      ELSE
        ! NO3    + h nu -> NO  + O2
        ino3a = k
      END IF

    CASE ('O2        ')
      ! O2     + h nu -> 2 O / O + O(1D)
      IF ((spj(k,3) == 'O(1D)     ')                                           &
        .OR. (spj(k,4) == 'O(1D)     ')) THEN
        io2b = k
      ELSE
        io2a = k
      END IF

      ! consider branching for O3
    CASE ('O3        ')
      IF ((spj(k,3) == 'O(1D)     ')                                           &
        .OR. (spj(k,4) == 'O(1D)     ')) THEN
        ! O3     + h nu -> O2  + O(1D)
        io3a = k
      ELSE
        ! O3     + h nu -> O2  + O
        io3b = k
      END IF

      ! consider branching for O3S
    CASE ('O3S       ')
      IF ((spj(k,3) == 'O(1D)S    ')                                           &
        .OR. (spj(k,4) == 'O(1D)S    ')) THEN
        ! O3S    + h nu -> O2  + O(1D)S
        io3sa = k
      ELSE
        ! O3S    + h nu -> O2  + O(3P)S
        io3sb = k
      END IF

      ! CH4    + h nu -> CH3   + H -> MeOO + HO2
    CASE ('CH4       ')
      ich4 = k

      ! H2O    + h nu -> OH    + H -> OH   + HO2
    CASE ('H2O       ','H2OS      ')
      ih2o = k

      ! NO     + h nu -> O     + N -> 2O    + NO
    CASE ('NO        ')
      ino = k

      ! N2O    + h nu -> O(1D) + N2
    CASE ('N2O       ')
      in2o = k

      ! F11    + h nu -> 3Cl
    CASE ('CFCl3     ')
      if11 = k

      ! F12    + h nu -> 2Cl
    CASE ('CF2Cl2    ')
      if12 = k

      ! ClONO2 + h nu -> Cl   + NO3 / ClO + NO2
    CASE ('ClONO2    ')
      IF ((spj(k,3) == 'Cl        ')                                           &
        .OR. (spj(k,4) == 'Cl        ')) THEN
        ! ClONO2 + h nu -> Cl + NO3
        iclono2a = k
      ELSE
        ! ClONO2 + h nu -> ClO + NO2
        iclono2b = k
      END IF

      ! HCl    + h nu -> H    + Cl
    CASE ('HCl       ')
      ihcl = k

      ! HOCl   + h nu -> OH   + Cl
    CASE ('HOCl      ')
      ihocl = k

      ! OClO   + h nu -> O    + ClO
    CASE ('OClO      ')
      ioclo = k

      ! Cl2O2  + h nu -> 2Cl  + O2
    CASE ('Cl2O2     ')
      icl2o2 = k

      ! BrO    + h nu -> Br   + O
    CASE ('BrO       ')
      ibro = k

      ! HOBr   + h nu -> Br   + OH
    CASE ('HOBr      ')
      ihobr = k

      ! BrONO2 + h nu -> Br   + NO3 / BrO + NO2
    CASE ('BrONO2    ')
      IF ((spj(k,3) == 'Br        ')                                           &
        .OR. (spj(k,4) == 'Br        ')) THEN
        ! BrONO2 + h nu -> Br + NO3
        ibrono2a = k
      ELSE
        ! BrONO2 + h nu -> BrO + NO2
        ibrono2b = k
      END IF

      ! BrCl   + h nu -> Br   + Cl
    CASE ('BrCl      ')
      ibrcl = k

      ! MeBr   + h nu -> Br
    CASE ('MeBr      ')
      imebr = k

      ! CCl4   + h nu -> 4 Cl
    CASE ('CCl4      ')
      iccl4 = k

      ! CF2ClCFCl2+h nu->3 Cl
    CASE ('CF2ClCFCl2')
      if113 = k

      ! CH3Cl +h nu   -> Cl
    CASE ('MeCl      ')
      imecl = k

      ! CH3CCl3 +h nu -> 3 Cl
    CASE ('MeCCl3    ')
      imcfm = k

      ! CHF2Cl +h nu -> Cl
    CASE ('CHF2Cl    ')
      if22 = k

      ! CBrClF2 +h nu -> Cl + Br
    CASE ('CF2ClBr   ')
      ih1211 = k

      ! CBrF3 +h nu   -> Br
    CASE ('CF3Br     ')
      ih1301 = k

      ! COF2 +h nu    -> 2F + CO
    CASE ('COF2      ')
      icof2 = k

      ! COFCl +h nu   -> F + Cl
    CASE ('COFCl     ')
      icofcl = k

      ! CO2 +h nu     -> CO + O(3P)
    CASE ('CO2       ')
      ico2 = k

      ! COS + h nu    -> CO + S
    CASE ('COS       ')
      icos = k

      ! HONO + h nu   -> OH + NO
    CASE ('HONO      ')
      ! ihono is not used

      ! MeONO2 + h nu -> HO2 + HCHO + NO2
    CASE ('MeONO2    ')
      ! imeono2 is not used

      ! CHBr3 + h nu  -> HBr + 2Br
    CASE ('CHBr3     ')
      ichbr3 = k

      ! CH2Br2 + h nu -> H2O + 2Br
    CASE ('CH2Br2    ')
      idbrm = k

      ! CS2 + h nu -> COS + SO2
    CASE ('CS2     ')
      ics2 = k
      ! H2SO4 + h nu -> SO3 + OH
    CASE ('H2SO4     ')
      ih2so4 = k
      ! SO3 + h nu -> SO2 + O(3P)
    CASE ('SO3     ')
      iso3 = k

    END SELECT
  END DO

  ! needed for solar cycle below
  current_month = i_month

  firstcall = .FALSE.
END IF ! firstcall

! Work out logic of whether to fill the tables. If at the beginning of the
! month and the tables haven't been filled, then fill the tables.
IF (tables_filled .AND. (current_month /= i_month)) THEN
  ! Reset tables filled indicator so that the tables will be re-filled
  tables_filled = .FALSE.
  current_month = i_month
END IF

IF (.NOT. tables_filled) THEN
  CALL inijtab(current_time,                                                   &
    ((photol_config%i_photol_scheme == i_scheme_fastjx) .AND.                  &
     (photol_config%fastjx_mode /= fjx_mode_2Donly) ) )
  tables_filled = .TRUE.
END IF

! CALCJS fills the stratospheric photolysis arrays AJxyz with sensible
! values.

theta_field_size = row_length * rows

CALL calcjs(1, theta_field_size,                                               &
  RESHAPE(cos_zenith_angle,[theta_field_size]),                                &
  RESHAPE(pressure,[theta_field_size]),                                        &
  RESHAPE(temp,[theta_field_size]),                                            &
  RESHAPE(ozonecol,[theta_field_size]),                                        &
  theta_field_size)

! here: use only existing photolysis reactions where pressure is less than
! 300 hPa, with a linear transition into stratospheric rates

DO i=1,rows
  DO j=1,row_length
    l=(i-1) * row_length + j

    IF (pressure(j,i) < 30000.0) THEN
      IF (pressure(j,i) < 20000.0) THEN
        frac = 1.0
      ELSE
        frac = (30000.0 - pressure(j,i))/10000.0
      END IF

      ! Merge tropospheric and stratospheric photolysis rates
      ! H2O2   + h nu -> 2 OH
      photrates(j,i,ih2o2) = frac  * ajh2o2(l)                                 &
        + (1.0-frac) * photrates(j,i,ih2o2)

      ! consider branching in case of HCHO
      ! HCHO   + h nu -> H2 + CO
      photrates(j,i,ihchob) = frac  * ajc2ob(l)                                &
        + (1.0-frac) * photrates(j,i,ihchob)

      ! HCHO   + h nu -> H + CHO -> CO + 2 HO2
      photrates(j,i,ihchoa) = frac  * ajc2oa(l)                                &
        + (1.0-frac) * photrates(j,i,ihchoa)

      ! HO2NO2 + h nu. Branching ratio from JPL (2002)
      IF (iho2no2a > 0) THEN
        IF (iho2no2b > 0) THEN
          ! HO2NO2 + h nu -> NO2 + HO2
          photrates(j,i,iho2no2a) = 0.667 * frac  * ajpna(l)                   &
            + (1.0-frac) * photrates(j,i,iho2no2a)
          ! HO2NO2 + h nu -> NO3 + OH
          photrates(j,i,iho2no2b) = 0.333 * frac  * ajpna(l)                   &
            + (1.0-frac) * photrates(j,i,iho2no2b)
        ELSE
          photrates(j,i,iho2no2a) = frac  * ajpna(l)                           &
            + (1.0-frac) * photrates(j,i,iho2no2a)
        END IF
      END IF

      ! HNO3   + h nu -> NO2 + OH
      photrates(j,i,ihno3) = frac  * ajhno3(l)                                 &
        + (1.0-frac) * photrates(j,i,ihno3)

      ! MeOOH  + h nu -> MeO + OH -> HCHO + HO2 + OH
      photrates(j,i,imeooh) = frac  * ajmhp(l)                                 &
        + (1.0-frac) * photrates(j,i,imeooh)

      ! N2O5   + h nu -> NO2 + NO3
      photrates(j,i,in2o5) = frac  * ajn2o5(l)                                 &
        + (1.0-frac) * photrates(j,i,in2o5)

      ! NO2    + h nu -> NO  + O
      photrates(j,i,ino2) = frac  * ajno2(l)                                   &
        + (1.0-frac) * photrates(j,i,ino2)

      ! consider branching for NO3
      ! NO3    + h nu -> NO2 + O
      photrates(j,i,ino3b) = frac  * ajno32(l)                                 &
        + (1.0-frac) * photrates(j,i,ino3b)

      ! NO3    + h nu -> NO  + O2
      photrates(j,i,ino3a) = frac  * ajno31(l)                                 &
        + (1.0-frac) * photrates(j,i,ino3a)

      ! O2     + h nu -> 2 O
      photrates(j,i,io2a) = frac  * aj2a(l)                                    &
        + (1.0-frac) * photrates(j,i,io2a)

      ! consider branching for O3
      ! O3     + h nu -> O2  + O(1D)
      photrates(j,i,io3a) = frac  * aj3a(l)                                    &
        + (1.0-frac) * photrates(j,i,io3a)
      ! O3     + h nu -> O2  + O
      photrates(j,i,io3b) = frac  * aj3(l)                                     &
        + (1.0-frac) * photrates(j,i,io3b)

      ! consider branching for O3S
      ! O3S    + h nu -> O2  + O(1D)S
      IF (io3sa > 0)                                                           &
        photrates(j,i,io3sa) = frac  * aj3a(l)                                 &
        + (1.0-frac) * photrates(j,i,io3sa)
      ! O3S    + h nu -> O2  + O(3P)S
      IF (io3sb > 0)                                                           &
        photrates(j,i,io3sb) = frac  * aj3(l)                                  &
        + (1.0-frac) * photrates(j,i,io3sb)

    END IF    ! pressure < 3000

    ! purely stratospheric photolysis rates
    ! SLIMCAT specific photolysis reactions

    ! O2     + h nu -> O + O(1D)
    IF (io2b > 0) THEN
      photrates(j,i,io2b) = aj2b(l)
    ELSE
      photrates(j,i,io2a) = photrates(j,i,io2a) + aj2b(l)
    END IF
    ! CH4    + h nu -> CH3   + H -> MeOO + HO2
    IF (ich4 > 0)    photrates(j,i,ich4) = ajch4(l)

    ! H2O    + h nu -> OH    + H -> OH   + HO2
    IF (ih2o > 0)    photrates(j,i,ih2o) = ajh2o(l)

    ! NO     + h nu -> O     + N -> 2O    + NO
    IF (ino > 0)     photrates(j,i,ino)  = ajno(l)

    ! N2O    + h nu -> O(1D) + N2
    IF (in2o > 0)    photrates(j,i,in2o) = ajn2o(l)

    ! F11    + h nu -> 3Cl
    IF (if11 > 0)    photrates(j,i,if11) = ajf11(l)

    ! F12    + h nu -> 2Cl
    IF (if12 > 0)    photrates(j,i,if12) = ajf12(l)

    ! ClONO2 + h nu -> Cl   + NO3
    IF (iclono2a > 0) THEN
      photrates(j,i,iclono2a) = ajcnita(l)
      IF (iclono2b == 0)                                                       &
        photrates(j,i,iclono2a) = photrates(j,i,iclono2a)                      &
        + ajcnitb(l)
    END IF

    ! ClONO2 + h nu -> ClO  + NO2
    IF (iclono2b > 0) THEN
      photrates(j,i,iclono2b) = ajcnitb(l)
      IF (iclono2a == 0)                                                       &
        photrates(j,i,iclono2b) = photrates(j,i,iclono2b)                      &
        + ajcnita(l)
    END IF

    ! HCl    + h nu -> H    + Cl
    IF (ihcl > 0)    photrates(j,i,ihcl) = ajhcl(l)

    ! HOCl   + h nu -> OH   + Cl
    IF (ihocl > 0)   photrates(j,i,ihocl) = ajhocl(l)

    ! OClO   + h nu -> O    + ClO
    IF (ioclo > 0)   photrates(j,i,ioclo) = ajoclo(l)

    ! Cl2O2  + h nu -> 2Cl  + O2
    IF (icl2o2 > 0)  photrates(j,i,icl2o2) = ajcl2o2(l)

    ! BrO    + h nu -> Br   + O
    IF (ibro > 0)    photrates(j,i,ibro) = ajbro(l)

    ! HOBr   + h nu -> Br   + OH
    IF (ihobr > 0)   photrates(j,i,ihobr) = ajhobr(l)

    ! BrONO2 + h nu -> Br   + NO3 / BrO + NO2
    ! Consider branching ratio (JPL, 2002)
    IF (ibrono2a*ibrono2b > 0) THEN
      photrates(j,i,ibrono2a) = 0.29*ajbrno3(l) ! Br + NO3 channel
      photrates(j,i,ibrono2b) = 0.71*ajbrno3(l) ! BrO + NO2 channel
    ELSE IF (ibrono2a > 0) THEN ! no branching
      photrates(j,i,ibrono2a) = ajbrno3(l)
    ELSE IF (ibrono2b > 0) THEN
      photrates(j,i,ibrono2b) = ajbrno3(l)
    END IF

    ! BrCl   + h nu -> Br   + Cl
    IF (ibrcl > 0)   photrates(j,i,ibrcl) = ajbrcl(l)

    ! MeBr   + h nu -> Br
    IF (imebr > 0)   photrates(j,i,imebr) = ajch3br(l)

    ! CCl4   + h nu -> 4 Cl
    IF (iccl4 > 0)   photrates(j,i,iccl4) = ajccl4(l)

    ! CF2ClCFCl2+h nu->3 Cl
    IF (if113 > 0)   photrates(j,i,if113) = ajf113(l)

    ! CH3Cl +h nu   -> Cl
    IF (imecl > 0)   photrates(j,i,imecl) = ajch3cl(l)

    ! CH3CCl3 +h nu -> 3 Cl
    IF (imcfm > 0)   photrates(j,i,imcfm) = ajmcfm(l)

    ! CHF2Cl +h nu -> Cl
    IF (if22 > 0)    photrates(j,i,if22) = ajf22(l)

    ! CBrClF2 +h nu -> Cl + Br
    IF (ih1211 > 0)  photrates(j,i,ih1211) = ajf12b1(l)

    ! CBrF3 +h nu   -> Br
    IF (ih1301 > 0)  photrates(j,i,ih1301) = ajf13b1(l)

    ! COF2 +h nu    -> 2F + CO
    IF (icof2 > 0)   photrates(j,i,icof2) = ajcof2(l)

    ! COFCl +h nu   -> F + Cl
    IF (icofcl > 0)  photrates(j,i,icofcl) = ajcofcl(l)

    ! CO2 +h nu     -> CO + O(3P)
    IF (ico2 > 0)    photrates(j,i,ico2) = ajco2(l)

    ! COS +h nu     -> CO + S
    IF (icos > 0)    photrates(j,i,icos) = ajcos(l)

    ! CHBr3 +h nu   -> HBr + 2Br
    IF (ichbr3 > 0)  photrates(j,i,ichbr3) = ajchbr3(l)

    ! CH2Br2 +h nu  -> H2O + 2Br
    IF (idbrm > 0)   photrates(j,i,idbrm) = ajdbrm(l)

    ! CS2 +h nu     -> COS + SO2
    IF (ics2 > 0)    photrates(j,i,ics2) = ajcs2(l)
    ! H2SO4 +h nu     -> SO3 + OH
    IF (ih2so4 > 0)  photrates(j,i,ih2so4) = ajh2so4(l)
    ! SO3 +h nu     -> SO2 + O(3P)
    IF (iso3 > 0)    photrates(j,i,iso3) = ajso3(l)

  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE strat_photol

END MODULE ukca_um_strat_photol_mod
