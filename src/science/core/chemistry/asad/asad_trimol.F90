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
! Purpose: Calculates trimolecular rate coefficients
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     Method
!     ------
!     See the IUPAC reference material on their website for details on
!     calculation of termolecular rates.
!     http://www.iupac-kinetic.ch.cam.ac.uk/
!
!     Local variables
!     ---------------
!     zo         Low pressure limit to rate*density
!     zi         High pressure limit to rate
!     zr         Ratio of zo/zi
!     iho2       Reaction index for HO2+HO2+M
!     ih2o       Array index for advected tracer H2O
!     in2o5      Reaction index for N2O5+M
!     ino2no3    Reaction index for NO2+NO3+M
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_trimol_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_TRIMOL_MOD'

CONTAINS

SUBROUTINE asad_trimol(n_points)

USE asad_mod,        ONLY: rk, at, ntrkx, spt, t300, t, peps,                  &
                           tnd, f, wp, specf, jpcspf, jptk
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ukca_um_legacy_mod, ONLY: exp_v, powr_v, oneover_v

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_points

!       Local variables

INTEGER, SAVE :: ih2o              ! Index for h2o in tracer array
INTEGER       :: iho2              ! Index for ho2+ho2 in rk array
INTEGER       :: in2o5             ! Reaction index for N2O5+M
INTEGER       :: ino2no3           ! Reaction index for NO2+NO3+M
INTEGER       :: idmsoho2_a        ! Reaction index for DMS+OH+O2
INTEGER       :: idmsoho2_b        ! Reaction index for DMS+OH+O2

INTEGER       :: j                 ! Loop variable
INTEGER       :: jl                ! Loop variable
INTEGER       :: jtr               ! Loop variable
INTEGER       :: jr                ! Index

REAL :: zo                         ! k_0
REAL :: zi                         ! k_infinity
REAL :: zfc                        ! F_c
REAL :: zr                         ! k_0/k_infinity
REAL, ALLOCATABLE, SAVE :: nf(:)   ! component of broadening factor
REAL :: nf2                        ! Covers a small change to nf

LOGICAL, SAVE :: first = .TRUE.
LOGICAL, SAVE :: first_pass = .TRUE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_TRIMOL'

REAL :: tmp(1:n_points)
REAL :: tmp_out(1:n_points)
REAL :: inv_t_power4(1:n_points)
REAL :: inv_t_power7(1:n_points)
REAL :: inv_t_power4_out(1:n_points)
REAL :: inv_t_power7_out(1:n_points)
REAL :: t_power3(1:n_points)
REAL :: t_power6(1:n_points)
REAL :: inv_t(1:n_points)


!       1.  Calculate trimolecular rate coefficients
!           --------- ------------ ---- ------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
iho2  = 0
in2o5   = 0
ino2no3 = 0
idmsoho2_a = 0
idmsoho2_b = 0

nf2 = 0.0
zfc = 0.0
zr =  0.0
zi = 0.0
zo = 0.0

CALL oneover_v(n_points, t, inv_t)

! OMP CRITICAL will only allow one thread through this code at a time,
! while the other threads are held until completion.
!$OMP CRITICAL (asad_trimol_init)
IF (first_pass) THEN
  IF (first) THEN
    ! Calculate the N factor to be used in broadening factor
    ALLOCATE(nf(jptk))
    nf(:) = 1.0
    WHERE (at(:,1) > 1e-3) nf(:) = 0.75 - 1.27*LOG10(at(:,1))

    ! Check if H2O is an advected tracer
    ih2o = 0
    DO jtr = 1, jpcspf
      IF ( specf(jtr)  ==  'H2O       ' ) ih2o = jtr
    END DO
    first = .FALSE.
  END IF
  first_pass = .FALSE.
END IF
!$OMP END CRITICAL (asad_trimol_init)

DO j = 1, jptk
  jr = ntrkx(j)

  IF ( spt(j,1) == 'HO2    ' .AND. spt(j,2) == 'HO2    ' )                     &
       iho2 = jr

  !         N2O5 + M
  !         reset to zero each step

  in2o5 = 0
  IF (spt(j,1) == 'N2O5      ') in2o5 = jr

  !         NO2 + NO3 + M
  !         reset to zero each step

  ino2no3 = 0
  IF ((spt(j,1) == 'NO2       ' .AND. spt(j,2) == 'NO3       ')                &
       .OR. (spt(j,2) == 'NO2       ' .AND.                                    &
       spt(j,1) == 'NO3       ')) ino2no3 = jr

  ! Catch for DMS + OH + O2 reaction
  idmsoho2_a = 0
  IF (spt(j,1) == 'DMS       ' .AND. spt(j,2) == 'OH        '                  &
      .AND. spt(j,3) == 'DMSO      ' .AND.                                     &
      spt(j,4) == 'HO2       ')  idmsoho2_a = jr

  ! Other DMS +OH + O2 branch in StratTrop
  idmsoho2_b = 0
  IF (spt(j,1) == 'DMS       ' .AND. spt(j,2) == 'OH        '                  &
      .AND. spt(j,3) == 'SO2       ' .AND.                                     &
      spt(j,4) == 'MeOO      ')  idmsoho2_b = jr


  CALL powr_v(n_points,t300,at(j,3),t_power3)
  CALL powr_v(n_points,t300,at(j,6),t_power6)
  inv_t_power4(1:n_points) = -at(j,4)*inv_t(1:n_points)
  inv_t_power7(1:n_points) = -at(j,7)*inv_t(1:n_points)
  CALL exp_v(n_points,inv_t_power4,inv_t_power4_out)
  CALL exp_v(n_points,inv_t_power7,inv_t_power7_out)
  DO jl = 1, n_points
    zo = at(j,2) * t_power3(jl) *                                              &
                   inv_t_power4_out(jl)  * tnd(jl)
    zi = at(j,5) * t_power6(jl) * inv_t_power7_out(jl)
    IF ( zo < peps ) THEN
      rk(jl,jr) = zi
    ELSE IF ( zi < peps ) THEN
      rk(jl,jr) = zo
    ELSE

      ! Special case for DMS+OH+O2 -> DMSO + HO2 / SO2 + MeOO
      ! k = (9.5e-39 * [O2] * exp(5270/T))/(1+7.5e-29[O2]exp(5610/T))
      ! IUPAC 1997+
      ! http://iupac.pole-ether.fr/htdocs/datasheets/pdf/SOx22_HO_CH3SCH3.pdf
      IF (idmsoho2_a /= 0 .OR. idmsoho2_b /= 0) THEN
        ! Multiply by 0.21 to give in terms of [O2] instead of [M]
        zo = zo * 0.2095
        zi = zi * 0.2095 * tnd(jl)
        rk(jl,jr) = zo/(1 + zi)
      ELSE
        nf2 = nf(j)
        IF ( at(j,1) <= 1.0 ) THEN
          zfc = at(j,1)
        ELSE IF ( in2o5 /= 0 .OR. ino2no3 /= 0 ) THEN ! dependent
                                                      ! reactions
          zfc = 2.5*EXP(-1950.0/t(jl))+0.9*EXP(-t(jl)/at(j,1))
        ELSE              ! temperature dependent Fc
          zfc = EXP( -t(jl)/at(j,1) )
          nf2 = 0.75 - 1.27*LOG10(zfc)
        END IF
        zr = zo / zi
        rk(jl,jr) = (zo/(1.0+zr)) *                                            &
                     zfc**(1.0/(1.0 + (LOG10(zr)/nf2)**2))
      END IF
    END IF
  END DO
END DO              ! end of loop over jptk

!       2. Dependent reactions.
!          --------- ----------

! HO2 + HO2 [+ M]
IF (ih2o /= 0 .AND. iho2 /= 0 ) THEN
  !         h2o is an advected tracer
  tmp(1:n_points) = 2200.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  DO jl = 1, n_points
    rk(jl,iho2) = rk(jl,iho2) *                                                &
    ( 1.0 + 1.4e-21*f(jl,ih2o)*tmp_out(jl) )
  END DO
ELSE IF (ih2o == 0 .AND. iho2 /= 0) THEN
  !         use modelled water concentration
  tmp(1:n_points) = 2200.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  DO jl = 1, n_points
    rk(jl,iho2) = rk(jl,iho2) *                                                &
    ( 1.0 + 1.4e-21*wp(jl)*tnd(jl)*tmp_out(jl) )
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_trimol
END MODULE asad_trimol_mod
