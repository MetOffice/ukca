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
! Purpose: Controlling routine for the IMPACT integration scheme
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
!     Interface
!     ---------
!     Called from chemistry driver routine cdrive.
!
!     Arguments:
!          lphot - Logical array. One entry per gridpoint. If set .T.
!                  then photolysis is occurring and the scheme will
!                  adjust the accuracy of the Jacobian. See method.
!
!     All other variables are obtained via modules.
!
!     Method
!     ------
!     The IMPACT scheme used the Euler backward implicit timescheme
!     as its starting point. This is a first order conservative
!     scheme. A predictor-corrector approach is used to solve this
!     by employing a Newton-Raphson type iteration. However, rather
!     than compute and invert a near full Jacobian matrix, this scheme
!     only considers the diagonal terms only. For mildly stiff chemistry
!     this should work fine. However, if the scheme is not converging
!     rapidly enough during daylight, this is a symptom that the
!     off-diagonal production terms, which are usually ignored, are
!     important. In this case, set the switch ljacx .T. and make sure
!     that the lphot array is set .true. where photolysis is on.
!     Since, only approximate forms of the Jacobian are considered, if
!     switched on, the method will add in the terms in the Jacobian
!     due to photolysis. It does this in an approximate way to avoid
!     an matrix inversion but this greatly improves the performance of
!     the scheme under stiffer conditions.
!
!     The IMPACT scheme is fully described in the forthcoming paper;
!     Carver and Stott, 1997, 'IMPACT: an implicit time integration
!     scheme for chemical species and families', to be submitted to
!     Annales Geophysicae. Preprints available.
!
!     Externals
!     ---------
!     ftoy         Converts f concentrations to y (partitions
!                  families).
!     diffun       Computes tendencies due to chemistry and
!                  main diagonal Jacobian elements.
!
!     Local variables
!     ---------------
!     gconv         .true. if convergence in Newton-Rhapson has been
!                   achieved.
!     ifi           Number of ftoy iterations.
!     zf            Values of f at start of chemistry step.
!     zprf          Value of f on previous Newton-Rhapson iteration.
!     binv          Inverse of the approximate (main diagonal only)
!                   of the matrix ( I - dt.J) where J is Jacobian.
!     dely          Represents the change in variable, y, produced by
!                   multiplying binv to the function we are trying to
!                   find the root of.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_impact_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_IMPACT_MOD'

CONTAINS

SUBROUTINE asad_impact(n_points, ix, jy, nlev)

USE asad_mod,        ONLY: cdt, ctype, ej, f, fdot, jpif,                      &
                           linfam, ljacx, madvtr, majors, moffam,              &
                           nlpdv, nltrim, nitnr, nprkx,                        &
                           nspi, peps, pmintnd, ptol, rk, y, jpcspf
USE ukca_config_specification_mod, ONLY: ukca_config
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook


USE asad_diffun_mod, ONLY: asad_diffun
USE asad_inimpct_mod, ONLY: asad_inimpct
USE asad_jac_mod, ONLY: asad_jac
USE asad_ftoy_mod, ONLY: asad_ftoy
IMPLICIT NONE


INTEGER, INTENT(IN) :: n_points   ! No of spatial points
INTEGER, INTENT(IN) :: ix         ! i counter
INTEGER, INTENT(IN) :: jy         ! j counter
INTEGER, INTENT(IN) :: nlev       ! Model level

!       Local variables

INTEGER :: lphot(n_points)
INTEGER :: inl
INTEGER :: nl
INTEGER :: j                        ! Loop variable
INTEGER :: jit                      ! Loop variable
INTEGER :: jl                       ! Loop variable
INTEGER :: jr                       ! Loop variable
INTEGER :: jtr                      ! Loop variable
INTEGER :: jt1                      ! Loop variable
INTEGER :: isp                      ! Index
INTEGER :: ifi
INTEGER :: ntr1
INTEGER :: ntr2
INTEGER :: ipos0
INTEGER :: ipos
INTEGER :: ir
INTEGER :: irk
INTEGER :: njr
INTEGER :: ireac
INTEGER :: iprod

REAL :: zf(n_points,jpcspf)
REAL :: zprf(n_points,jpcspf)
REAL :: binv(n_points,jpcspf)
REAL :: dely(n_points,jpcspf)
REAL :: dd(n_points)
REAL :: corrn(n_points)

LOGICAL :: gconv
LOGICAL :: gfam
LOGICAL, SAVE :: gfirst = .TRUE.
LOGICAL, SAVE :: first_pass = .TRUE.
LOGICAL :: not_first_call = .FALSE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_IMPACT'


!       1.  Predictor step (first guess)
!           --------- ---- ------ ------

!       Crude way of determining if photolysis is on or not.
!       Needs improving.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
inl = 0
DO jl = 1, n_points
  IF ( rk(jl,nprkx(1)) > peps ) THEN
    inl = inl + 1
    lphot(inl) = jl
  END IF
END DO

IF ( first_pass ) THEN
  ! OMP CRITICAL will only allow one thread through this code at a time,
  ! while the other threads are held until completion.
!$OMP CRITICAL (asad_impact_init)
  IF ( gfirst ) THEN
    gfirst = .FALSE.
    CALL asad_inimpct
  END IF
  first_pass = .FALSE.
!$OMP END CRITICAL (asad_impact_init)
END IF

!       1.1  Do the linearised first guess to give first approx.
!            solution at y(n+1).

nl = n_points
DO jtr = 1, jpcspf
  isp = majors(jtr)
  DO jl = 1, n_points
    zf(jl,jtr) = f(jl,jtr)

    !           1.2  Test Jacobian has not gone positive. Possible since
    !                we only use approximate form. If it has, use
    !                an explicit step.

    IF ( ej(jl,jtr) > 0.0 ) THEN
      f(jl,jtr) = zf(jl,jtr) + cdt*fdot(jl,jtr)
    ELSE
      f(jl,jtr) = zf(jl,jtr) + ( cdt*fdot(jl,jtr) )                            &
                       / ( 1.0 - cdt*ej(jl,jtr) )
    END IF
    IF ( linfam(jl,jtr) ) f(jl,jtr) = y(jl,isp)

  END DO
END DO

!       2.  Newton-Rhapson iteration (corrector step).
!           -------------- --------- ---------- ------

DO jit = 1, ukca_config%nrsteps

  !         2.1  Decide if we are recomputing ratios or not.

  IF ( jit < 4 ) THEN
    ifi = 0
  ELSE
    ifi = nitnr
  END IF

  !         2.2  Work out rates of change and main diagonal of J.

  ! Call to asad_ftoy includes iteration count
  CALL asad_ftoy(not_first_call, ifi, jit, n_points, ix, jy, nlev)
  CALL asad_diffun( nl )
  CALL asad_jac(n_points)

  !         2.3  Do the normal Newton-Raphson iteration with
  !              Jacobian approximated to main diagonal.

  DO jtr = 1, jpcspf
    isp = majors(jtr)
    DO jl = 1, n_points
      zprf(jl,jtr) = f(jl,jtr)
      binv(jl,jtr) = 1.0 / ( 1.0 - cdt*ej(jl,jtr) )
      dely(jl,jtr) = binv(jl,jtr) *                                            &
             ( (zf(jl,jtr) - f(jl,jtr)) + cdt*fdot(jl,jtr) )
      f(jl,jtr) = f(jl,jtr) + dely(jl,jtr)

      IF ( linfam(jl,jtr) ) f(jl,jtr) = y(jl,isp)
    END DO
  END DO

  !         2.4  If requested for all points where photolysis is
  !              occurring, include the terms in the Jacobian
  !              arising from just the photolysis terms. Note
  !              that we only do this for points where photolysis
  !              is actually turned on.

  IF ( ljacx .AND. inl > 0 ) THEN
    ntr1 = nltrim(0,1)

    !           2.4.1  For each tracer that needs correcting (ie.
    !                  loop over the rows in the Jacobian matrix).

    DO jt1 = 1, ntr1

      !             2.4.2  Compute the contribution from non-zero
      !                    elements in the Jacobian due to photolysis.
      !                    ie. work our way along the columns. We first
      !                    compute the contribution from ordinary tracers
      !                    'TR' and then from families.

      ntr2  = nltrim(jt1,2)
      ipos0 = nltrim(jt1,3)
      ipos  = ipos0
      DO jl = 1, n_points
        corrn(jl) = 0.0
      END DO

      undef_loop: DO

        ir   = nlpdv(ipos,1)
        irk  = nprkx(ir)
        isp  = nspi(irk,1)
        njr  = nlpdv(ipos,2)
        ireac = madvtr(isp)
        gfam = ireac  ==  0
        DO j = 1, inl
          dd(j) = 0.0
        END DO

        IF ( gfam ) THEN
          ireac = moffam(isp)
          DO jr = 1, njr
            ir  = nlpdv(ipos,1)
            irk = nprkx(ir)
            isp = nspi(irk,1)
            DO j = 1, inl
              jl    = lphot(j)
              dd(j) = dd(j) + rk(jl,irk)*y(jl,isp)
            END DO
            ipos = ipos + 1
          END DO
          DO j = 1, inl
            jl    = lphot(j)
            dd(j) = dd(j) / zprf(jl,ireac)
          END DO

        ELSE

          !             2.4.3.  Species of type 'FT' will come here.
          !                     Check, at each point, whether it's gone
          !                     into the family or not.

          DO jr = 1, njr
            ir  = nlpdv(ipos,1)
            irk = nprkx(ir)
            isp = nspi(irk,1)
            DO j = 1, inl
              jl = lphot(j)
              IF ( linfam(jl,ireac) ) THEN
                dd(j) = dd(j)+rk(jl,irk)*y(jl,isp)/zprf(jl,ireac)
              ELSE
                dd(j) = dd(j)+rk(jl,irk)
              END IF
            END DO
            ipos = ipos + 1
          END DO
        END IF

        DO j = 1, inl
          jl = lphot(j)
          corrn(j) = corrn(j) + dd(j)*dely(jl,ireac)
        END DO

        IF ( ipos - ipos0 >= ntr2 ) EXIT undef_loop
      END DO undef_loop

      !             2.4.3  Now add correction to the tracer (rows). If the
      !                    tracer is of type 'FT', then for each pt, we must
      !                    check whether the tracer has been put into the f
      !                    or not. If it has, we add the correction to the
      !                    family and not to the tracer.

      iprod = nltrim(jt1,1)
      isp   = majors(iprod)
      IF ( ctype(isp) == jpif ) THEN
        DO j = 1, inl
          jl    = lphot(j)
          iprod = nltrim(jt1,1)
          IF ( linfam(jl,iprod) ) iprod = moffam(isp)
          f(jl,iprod) = f(jl,iprod)+cdt*binv(jl,iprod)*corrn(j)
        END DO
      ELSE
        DO j = 1, inl
          jl = lphot(j)
          f(jl,iprod) = f(jl,iprod)+cdt*binv(jl,iprod)*corrn(j)
        END DO
      END IF

    END DO
  END IF      ! End of IF ( ljacx .AND. inl > 0 ) statement

  !         9. Check for convergence
  !            ----- --- -----------

  gconv = .TRUE.
  jtr_loop: DO jtr = 1,jpcspf
    DO jl = 1, n_points
      IF ( ABS(f(jl,jtr)-zprf(jl,jtr)) >  ptol*f(jl,jtr)                       &
      .AND. f(jl,jtr) >  pmintnd(jl) ) gconv=.FALSE.
    END DO
    IF ( .NOT. gconv ) EXIT jtr_loop
  END DO jtr_loop
  IF (gconv) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END DO           ! End of NR (jit) loop

!       9.1  Convergence achieved or max. no. of iterations reached

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_impact
END MODULE asad_impact_mod
