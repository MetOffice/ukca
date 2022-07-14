! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************

MODULE ukca_binapara_mod

! Module containing the subroutine UKCA_BINAPARA which calculates
! parametrized values of nucleation rate for tropospheric and
! stratospheric conditions

IMPLICIT NONE

!  Description:
!    Calculates parametrized values of nucleation rate,
!    mole fraction of sulphuric acid
!    total number of particles, and the radius of the critical cluster
!    in H2O-H2SO4 system if temperature, saturation ratio of water and
!    sulfuric acid concentration are given.
!
!  Reference:
!    See H. Vehkamaki, et al., An improved parameterisation for
!    sulfuric-acid-water nucleation rates for tropospheric and
!    stratospheric conditions,
!    J. Geophys. Res., 107, (D22), 4622, doi:10.1029/2002JD002184, 2002.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford and The Met Office.
!  See:  www.ukca.ac.uk
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_BINAPARA_MOD'

CONTAINS

SUBROUTINE ukca_binapara(nbox,t_in,rh_in,h2so4_in,jveh,rc)

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ukca_um_legacy_mod, ONLY: exp_v, log_v

IMPLICIT NONE

INTEGER, INTENT(IN) :: nbox            ! No of points

REAL, INTENT(IN)    :: t_in(nbox)      ! temperature in K (190.15-300.15K)
REAL, INTENT(IN)    :: rh_in(nbox)     ! satur'n ratio of water (0.0001-1)
REAL, INTENT(IN)    :: h2so4_in(nbox)  ! sulfuric acid concentration in
!                                            ! 1/cm3 (10^4-10^11 1/cm3)

REAL, INTENT(OUT)   :: jveh(nbox)      ! nucleation rate in 1/cm3s
!                                              (10^-7-10^10 1/cm3s)
REAL, INTENT(OUT)   :: rc(nbox)        ! radius of the critical cluster
!                                            ! in nm

INTEGER :: jl                          ! counter

REAL :: t(nbox)                        ! from t_in, modified
REAL :: rh(nbox)                       ! from rh_in, modified
REAL :: tdegk                          ! temperature (K)
REAL :: tdegk2                         ! temperature (K) ^2
REAL :: tdegk3                         ! temperature (K) ^3
REAL :: logrh                          ! LOG(rh)
REAL :: logrh2                         ! LOG(rh) ^2
REAL :: logrh3                         ! LOG(rh) ^3
REAL :: logh2so4                       ! LOG(h2so4)
REAL :: logh2so42                      ! LOG(h2so4) ^2
REAL :: logh2so43                      ! LOG(h2so4) ^3
REAL :: h2so4(nbox)                    ! from h2so4_in, modified
REAL :: termx(nbox)                    ! molefraction of H2SO4 in the
!                                            ! critical cluster
REAL :: ntot(nbox)                     ! total number of molecules in the
!                                            ! critical cluster (ntot>4)
REAL :: ntot_out(nbox)                 ! temporary array
REAL :: jveh_out(nbox)                 ! temporary array
REAL :: rc_out(nbox)                   ! temporary array

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_BINAPARA'

REAL :: logrha(nbox)                   ! LOG(rh)
REAL :: logh2so4a(nbox)                ! LOG(h2so4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Perform check on whether Temp, RH and H2SO4 are out of bounds

t    (:) = t_in    (:)
rh   (:) = rh_in   (:)
h2so4(:) = h2so4_in(:)

WHERE (t(:) < 190.15) t(:)=190.15

WHERE (t(:) > 300.15) t(:)=300.15

WHERE (rh(:) < 0.0001) rh(:)=0.0001

WHERE (rh(:) > 1.0   ) rh(:)=1.0

WHERE (h2so4(:) < 1.0e4 ) h2so4(:)=1.0e4

WHERE (h2so4(:) > 1.0e11) h2so4(:)=1.0e11


CALL log_v(nbox, rh, logrha)
CALL log_v(nbox, h2so4, logh2so4a)

DO jl=1,nbox

  logrh = logrha(jl)
  logrh2 = logrh*logrh
  logrh3 = logrh2*logrh
  logh2so4 = logh2so4a(jl)
  logh2so42 = logh2so4*logh2so4
  logh2so43 = logh2so42*logh2so4
  tdegk = t(jl)
  tdegk2 = tdegk*tdegk
  tdegk3 = tdegk2*tdegk


  termx(jl)=  0.7409967177282139 - 0.002663785665140117*tdegk                  &
      + 0.002010478847383187*logrh                                             &
      - 0.0001832894131464668*tdegk*logrh                                      &
      + 0.001574072538464286*logrh2                                            &
      - 0.00001790589121766952*tdegk*logrh2                                    &
      + 0.0001844027436573778*logrh3                                           &
      -  1.503452308794887e-6*tdegk*logrh3                                     &
      - 0.003499978417957668*logh2so4                                          &
      + 0.0000504021689382576*tdegk*logh2so4


  jveh(jl)= 0.1430901615568665 + 2.219563673425199*tdegk -                     &
      0.02739106114964264*tdegk2 +                                             &
      0.00007228107239317088*tdegk3 + 5.91822263375044/termx(jl)+              &
      0.1174886643003278*logrh + 0.4625315047693772*tdegk*logrh -              &
      0.01180591129059253*tdegk2*logrh +                                       &
      0.0000404196487152575*tdegk3*logrh +                                     &
      (15.79628615047088*logrh)/termx(jl) -                                    &
      0.215553951893509*logrh2 -                                               &
      0.0810269192332194*tdegk*logrh2 +                                        &
      0.001435808434184642*tdegk2*logrh2 -                                     &
      4.775796947178588e-6*tdegk3*logrh2 -                                     &
      (2.912974063702185*logrh2)/termx(jl) -                                   &
      3.588557942822751*logrh3 +                                               &
      0.04950795302831703*tdegk*logrh3 -                                       &
      0.0002138195118737068*tdegk2*logrh3 +                                    &
      3.108005107949533e-7*tdegk3*logrh3 -                                     &
      (0.02933332747098296*logrh3)/termx(jl) +                                 &
      1.145983818561277*logh2so4 -                                             &
      0.6007956227856778*tdegk*logh2so4 +                                      &
      0.00864244733283759*tdegk2*logh2so4 -                                    &
      0.00002289467254710888*tdegk3*logh2so4 -                                 &
      (8.44984513869014*logh2so4)/termx(jl) +                                  &
      2.158548369286559*logrh*logh2so4 +                                       &
      0.0808121412840917*tdegk*logrh*logh2so4 -                                &
      0.0004073815255395214*tdegk2*logrh*logh2so4 -                            &
      4.019572560156515e-7*tdegk3*logrh*logh2so4 +                             &
      (0.7213255852557236*logrh*logh2so4)/termx(jl) +                          &
      1.62409850488771*logrh2*logh2so4 -                                       &
      0.01601062035325362*tdegk*logrh2*logh2so4 +                              &
      0.00003771238979714162*tdegk2*logrh2*logh2so4+                           &
      3.217942606371182e-8*tdegk3*logrh2*logh2so4 -                            &
      (0.01132550810022116*logrh2*logh2so4)/termx(jl)+                         &
      9.71681713056504*logh2so42 -                                             &
      0.1150478558347306*tdegk*logh2so42 +                                     &
      0.0001570982486038294*tdegk2*logh2so42 +                                 &
      4.009144680125015e-7*tdegk3*logh2so42 +                                  &
      (0.7118597859976135*logh2so42)/termx(jl) -                               &
      1.056105824379897*logrh*logh2so42 +                                      &
      0.00903377584628419*tdegk*logrh*logh2so42 -                              &
      0.00001984167387090606*tdegk2*logrh*logh2so42+                           &
      2.460478196482179e-8*tdegk3*logrh*logh2so42 -                            &
      (0.05790872906645181*logrh*logh2so42)/termx(jl)-                         &
      0.1487119673397459*logh2so43 +                                           &
      0.002835082097822667*tdegk*logh2so43 -                                   &
      9.24618825471694e-6*tdegk2*logh2so43 +                                   &
      5.004267665960894e-9*tdegk3*logh2so43 -                                  &
      (0.01270805101481648*logh2so43)/termx(jl)
  !     1/(cm3s)


  ntot(jl) =-0.002954125078716302 - 0.0976834264241286*tdegk +                 &
      0.001024847927067835*tdegk2-2.186459697726116e-6*tdegk3 -                &
      0.1017165718716887/termx(jl) - 0.002050640345231486*logrh -              &
      0.007585041382707174*tdegk*logrh +                                       &
      0.0001926539658089536*tdegk2*logrh -                                     &
      6.70429719683894e-7*tdegk3*logrh -                                       &
      (0.2557744774673163*logrh)/termx(jl) +                                   &
      0.003223076552477191*logrh2 +                                            &
      0.000852636632240633*tdegk*logrh2 -                                      &
      0.00001547571354871789*tdegk2*logrh2 +                                   &
      5.666608424980593e-8*tdegk3*logrh2 +                                     &
      (0.03384437400744206*logrh2)/termx(jl) +                                 &
      0.04743226764572505*logrh3 -                                             &
      0.0006251042204583412*tdegk*logrh3 +                                     &
      2.650663328519478e-6*tdegk2*logrh3 -                                     &
      3.674710848763778e-9*tdegk3*logrh3 -                                     &
      (0.0002672510825259393*logrh3)/termx(jl) -                               &
      0.01252108546759328*logh2so4 +                                           &
      0.005806550506277202*tdegk*logh2so4 -                                    &
      0.0001016735312443444*tdegk2*logh2so4 +                                  &
      2.881946187214505e-7*tdegk3*logh2so4 +                                   &
      (0.0942243379396279*logh2so4)/termx(jl) -                                &
      0.0385459592773097*logrh*logh2so4 -                                      &
      0.0006723156277391984*tdegk*logrh*logh2so4 +                             &
      2.602884877659698e-6*tdegk2*logrh*logh2so4 +                             &
      1.194163699688297e-8*tdegk3*logrh*logh2so4 -                             &
      (0.00851515345806281*logrh*logh2so4)/termx(jl) -                         &
      0.01837488495738111*logrh2*logh2so4 +                                    &
      0.0001720723574407498*tdegk*logrh2*logh2so4 -                            &
      3.717657974086814e-7*tdegk2*logrh2*logh2so4 -                            &
      5.148746022615196e-10*tdegk3*logrh2*logh2so4 +                           &
      (0.0002686602132926594*logrh2*logh2so4)/termx(jl) -                      &
      0.06199739728812199*logh2so42 +                                          &
      0.000906958053583576*tdegk*logh2so42 -                                   &
      9.11727926129757e-7*tdegk2*logh2so42 -                                   &
      5.367963396508457e-9*tdegk3*logh2so42 -                                  &
      (0.007742343393937707*logh2so42)/termx(jl) +                             &
      0.0121827103101659*logrh*logh2so42 -                                     &
      0.0001066499571188091*tdegk*logrh*logh2so42 +                            &
      2.534598655067518e-7*tdegk2*logrh*logh2so42 -                            &
      3.635186504599571e-10*tdegk3*logrh*logh2so42 +                           &
      (0.0006100650851863252*logrh*logh2so42)/termx(jl) +                      &
      0.0003201836700403512*logh2so43 -                                        &
      0.0000174761713262546*tdegk*logh2so43 +                                  &
      6.065037668052182e-8*tdegk2*logh2so43 -                                  &
      1.421771723004557e-11*tdegk3*logh2so43 +                                 &
      (0.0001357509859501723*logh2so43)/termx(jl)
  rc(jl)=-1.6524245+0.42316402*termx(jl)+0.33466487*                           &
               ntot(jl)

END DO       ! jl

CALL exp_v(nbox, ntot, ntot_out)
CALL exp_v(nbox, jveh, jveh_out)
CALL exp_v(nbox, rc  , rc_out  )
jveh=jveh_out
rc=rc_out
DO jl=1,nbox

  IF (ntot_out(jl) < 4.0 ) THEN
    IF ( t(jl) < 195.15) THEN
      !            print *,'Warning: critical cluster < 4 molecules',
      !     &              ' and T<195.15 : rough estimate J=1e5/(cm3s)'
      jveh(jl) = 1.0e5
    ELSE
      !            print *,'Warning: critical cluster < 4 molecules:',
      !     &              ' accuracy of the parameterisation reduced,',
      !     &              ' J_theor/J_para up to 50'
    END IF
  END IF


  IF (jveh(jl) < 1.0e-7) THEN
    !         print *,'Warning: nucleation rate < 1e-7/cm3s,',
    !     &           ' using 0.0/cm3s'
    jveh(jl) = 0.0
  END IF

  IF (jveh(jl) > 1.0e10) THEN
    !         print *,'Warning: nucleation rate > 1e10/cm3s,',
    !     &           ' using 1e10/cm3s'
    jveh(jl) = 1.0e10
  END IF

END DO       ! jl

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_binapara

END MODULE ukca_binapara_mod
