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
! Purpose: To compute the main diagonal of the Jacobian matrix df/dy
!     where f is dy/dt. The IMPACT time integration scheme uses
!     the main diagonal.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_IMPACT
!
!     Interface
!     ---------
!        Called by impact.
!
!        n_points - set on entry
!
!        On exit, the array ej will have been updated.
!
!     Method
!     ------
!     The production term resulting from  a particular reaction
!     is given by P = k*Y1*Y2, or P = J*Y1 for photolysis and
!     P = K*Y1 for thermal decomposition. It is the same for all
!     products. The loss term for reactant 1 can be either
!     L1 = k*Y2 or L1 = J for photolysis and L1 = K for thermal
!     decomposition. The loss rate for reactant 2 , if present,
!     is L2 = K*Y1.  These terms are stored in the array prk.
!
!     The main diagonal Jacobian elements are calculated for
!     each family/tracer by differentiating d(fdot)/df.
!     Considering the term due to each individual reaction:
!     If neither of the reactants are
!     members of the family/tracer, the term for the reaction is
!     zero. However, if one of reactants is a family member, a
!     term  - N*k*Y1*(Y2/F) should be added to the Jacobian element.
!     N is the net number of odd atoms of that family/tracer type
!     resulting from the reaction (calculated in inodd). Since
!     P = k*Y1*Y2 has already been calculated, we add N*P to the
!     running sum. The division by F, the family concentration is
!     done at the end of the subroutine.
!     If both reactants are members of the same family/tracer, the
!     term due to the reaction is  - 2*N*k*Y1*(Y2/F) where the
!     factor 2 results from the differentiation of the F^2 term,
!     and so 2*N*P is added to the running sum.
!
!     Deposition terms:
!     Since emissions are generally atmospheric concentration
!     independent, they do not enter the expression for the main
!     diagonal Jacobian elements. However, deposition is a first order
!     loss process and so the Jacobian elements must also be updated.
!     The loss rate of a family due to deposition of one of its
!     members is given by -N*dp*(Y/F)*F where N is the number of odd
!     atoms in the species of the family type and dp is the loss rate.
!     The Jacobian element is given by dfdot/df, and so a term =
!     -N*dp*(Y/F) must be added to the Jacobian element. Since the
!     final division by F is done at the end, the term added to the
!     sum kept in memory is -N*dp*Y.
!
!     Externals
!     ---------
!          None.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_jac_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_JAC_MOD'

CONTAINS

SUBROUTINE asad_jac(n_points)

USE asad_mod, ONLY: ctype, dpd, dpw, ej, f, jpcspf, jpif, jpfm,                &
                    jpmsp, jpnr, jpspec, linfam, madvtr, moffam,               &
                    ndepd, ndepw, njacx1, njacx2, njacx3, njcgrp,              &
                    nltr3, nltrf, nmpjac, nodd, npjac1, nspi,                  &
                    ntr3, ntrf, peps, prk, y
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) :: n_points

!       Local variables

INTEGER :: ifam(jpmsp)
INTEGER :: itr(jpmsp)
INTEGER :: j            ! Loop variable
INTEGER :: jc           ! Loop variable
INTEGER :: j1           ! Loop variable
INTEGER :: j2           ! Loop variable
INTEGER :: j3           ! Loop variable
INTEGER :: jl           ! Loop variable
INTEGER :: jp           ! Loop variable
INTEGER :: jr           ! Loop variable
INTEGER :: js           ! Loop variable
INTEGER :: jtr          ! Loop variable
INTEGER :: ip1          ! Index
INTEGER :: ip2          ! Index
INTEGER :: ip3          ! Index
INTEGER :: i1           ! Index
INTEGER :: i2           ! Index
INTEGER :: i3           ! Index
INTEGER :: ir1          ! Index
INTEGER :: ir2          ! Index
INTEGER :: itrcr        ! Index
INTEGER :: ks
INTEGER :: ifamd
INTEGER :: itrd

LOGICAL :: gif(jpmsp)
LOGICAL :: gtype(jpmsp)
LOGICAL :: gir1
LOGICAL :: gir2
LOGICAL :: ginclude2
LOGICAL :: ginclude3
LOGICAL :: ginclude4
LOGICAL :: ginclude5

CHARACTER :: ityped*2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_JAC'


!       1.  Initialise local variables.
!           ---------- ----- ----------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO jtr = 1, jpcspf
  DO jl = 1, n_points
    ej(jl,jtr) = 0.0
  END DO
END DO

!       2.  Calc. main diagonal of Jacobian matrix.
!           ----- ---- -------- -- -------- -------

!           For those model families and non-family species
!           that don't have species that swap in and out of
!           the family we can unroll the Jacobian summation.

DO jc = 1, ntrf
  itrcr = nltrf(jc)
  ip3   = njcgrp(itrcr,3)
  ip2   = njcgrp(itrcr,2)
  ip1   = njcgrp(itrcr,1)

  !         2.1  Terms in groups of 3 (negative contrib.)

  DO j3 = 1, ip3
    i1 = njacx3(1,j3,itrcr)
    i2 = njacx3(2,j3,itrcr)
    i3 = njacx3(3,j3,itrcr)

    DO jl = 1, n_points
      ej(jl,itrcr) = ej(jl,itrcr)                                              &
                   - prk(jl,i1) - prk(jl,i2) - prk(jl,i3)
    END DO
  END DO

  !         2.2 Terms in groups of 2; this loop is either
  !             0 or 1.

  DO j2 = 1, ip2
    i1 = njacx2(1,itrcr)
    i2 = njacx2(2,itrcr)
    DO jl = 1, n_points
      ej(jl,itrcr) = ej(jl,itrcr) - prk(jl,i1) - prk(jl,i2)
    END DO
  END DO

  !         2.3  Terms in groups of 3; either 0 or 1.

  DO j1 = 1, ip1
    i1 = njacx1(itrcr)
    DO jl = 1, n_points
      ej(jl,itrcr) = ej(jl,itrcr) - prk(jl,i1)
    END DO
  END DO

  !         2.4  Terms that make a positive contribution to the
  !              Jacobian. This isn't unrolled as there are not
  !              likely to be many of them.

  DO jp = 1, nmpjac(itrcr)
    i1 = npjac1(jp,itrcr)
    DO jl = 1, n_points
      ej(jl,itrcr) = ej(jl,itrcr) + prk(jl,i1)
    END DO
  END DO

END DO

!       3.  For families with in/out species.
!           --- -------- ---- ------ --------

!           For families with in/out species, we have to recompute
!           which reactions contribute every dynamical step since
!           this now varies with gridpoint. Duplicates the code in
!           inijac.f

IF ( ntr3 /= 0 ) THEN
  DO jr = 1, jpnr
    ir1 = nspi(jr,1)
    ir2 = nspi(jr,2)
    ip1 = nspi(jr,3)
    ip2 = nspi(jr,4)
    ip3 = nspi(jr,5)

    DO j = 1, jpmsp
      gif(j)   = .FALSE.
      gtype(j) = .FALSE.
      ifam(j)  = 0
      itr(j)   = 0
    END DO
    IF ( ir1 /= 0 ) THEN
      ifam(1) = moffam(ir1)
      itr(1)  = madvtr(ir1)
    END IF
    IF ( ir2 /= 0 ) THEN
      ifam(2) = moffam(ir2)
      itr(2)  = madvtr(ir2)
    END IF

    gir1 = .FALSE.
    gir2 = .FALSE.
    DO j3 = 1, ntr3
      gir1 = gir1 .OR. ifam(1) == nltr3(j3)
      gir2 = gir2 .OR. ifam(2) == nltr3(j3)
    END DO

    IF ( gir1 .OR. gir2 ) THEN

      !             3.1  Set indices to determine species type.

      IF ( ir1 /= 0 ) THEN
        gtype(1) = ctype(ir1)  ==  jpfm
        gif(1) = ctype(ir1)  ==  jpif
      END IF
      IF ( ir2 /= 0 ) THEN
        gtype(2) = ctype(ir2)  ==  jpfm
        gif(2) = ctype(ir2)  ==  jpif
      END IF
      IF ( ip1 /= 0 ) THEN
        ifam(3) = moffam(ip1)
        gtype(3) = ctype(ip1)  ==  jpfm
        gif(3) = ctype(ip1)  ==  jpif
        itr(3) = madvtr(ip1)
      END IF
      IF ( ip2 /= 0 ) THEN
        ifam(4) = moffam(ip2)
        gtype(4) = ctype(ip2)  ==  jpfm
        gif(4) = ctype(ip2)  ==  jpif
        itr(4) = madvtr(ip2)
      END IF
      IF ( ip3 /= 0 ) THEN
        ifam(5) = moffam(ip3)
        gtype(5) = ctype(ip3)  ==  jpfm
        gif(5) = ctype(ip3)  ==  jpif
        itr(5) = madvtr(ip3)
      END IF
    END IF

    !           3.2  If first reactant is a family member.

    IF ( gir1 ) THEN
      DO jl = 1, n_points
        ginclude2 = (gif(2) .AND. linfam(jl,itr(2)) .OR. gtype(2))
        ginclude3 = (gif(3) .AND. linfam(jl,itr(3)) .OR. gtype(3))
        ginclude4 = (gif(4) .AND. linfam(jl,itr(4)) .OR. gtype(4))
        ginclude5 = (gif(5) .AND. linfam(jl,itr(5)) .OR. gtype(5))
        IF ( gif(1) .AND. .NOT. linfam(jl,itr(1)) ) THEN
          ks = 0
        ELSE
          ks = -nodd(ir1)
          IF (ifam(2)==ifam(1) .AND. ginclude2) ks = ks-nodd(ir2)
          IF (ifam(3)==ifam(1) .AND. ginclude3) ks = ks+nodd(ip1)
          IF (ifam(4)==ifam(1) .AND. ginclude4) ks = ks+nodd(ip2)
          IF (ifam(5)==ifam(1) .AND. ginclude5) ks = ks+nodd(ip3)
        END IF
        IF ( ks  /=  0 )                                                       &
          ej(jl,ifam(1)) = ej(jl,ifam(1)) + prk(jl,jr) * ks
      END DO
    END IF

    !           3.3  If second reactant is a family member.

    IF ( gir2 ) THEN
      DO jl = 1, n_points
        ginclude2 = (gif(2) .AND. linfam(jl,itr(2)) .OR. gtype(2))
        ginclude3 = (gif(3) .AND. linfam(jl,itr(3)) .OR. gtype(3))
        ginclude4 = (gif(4) .AND. linfam(jl,itr(4)) .OR. gtype(4))
        ginclude5 = (gif(5) .AND. linfam(jl,itr(5)) .OR. gtype(5))
        IF ( gif(2) .AND. .NOT. linfam(jl,itr(2)) ) THEN
          ks = 0
        ELSE
          ks = -nodd(ir2)
          IF (ifam(1)==ifam(2) .AND. ginclude2) ks = ks-nodd(ir1)
          IF (ifam(3)==ifam(2) .AND. ginclude3) ks = ks+nodd(ip1)
          IF (ifam(4)==ifam(2) .AND. ginclude4) ks = ks+nodd(ip2)
          IF (ifam(5)==ifam(2) .AND. ginclude5) ks = ks+nodd(ip3)
        END IF
        IF ( ks  /=  0 )                                                       &
          ej(jl,ifam(2)) = ej(jl,ifam(2)) + prk(jl,jr) * ks
      END DO
    END IF

  END DO

END IF           ! End of IF ( ntr3 /= 0 ) statement

!       4.  Add deposition terms to Jacobian diagonal.
!           --- ---------- ----- -- -------- ---------

IF ( ndepw /= 0 .OR. ndepd /= 0 ) THEN
  DO js = 1, jpspec
    ifamd = moffam(js)
    itrd = madvtr(js)
    ityped = ctype(js)

    IF ( ifamd /= 0 ) THEN
      DO jl = 1, n_points
        IF ((ityped == jpfm) .OR.                                              &
           (ityped == jpif .AND. linfam(jl,itrd)))                             &
          ej(jl,ifamd) = ej(jl,ifamd) - nodd(js) *                             &
                    ( dpd(jl,js)+dpw(jl,js)) * y(jl,js)
      END DO
    END IF
    IF ( itrd /= 0 ) THEN
      DO jl = 1, n_points
        ej(jl,itrd) = ej(jl,itrd)                                              &
                    - (dpd(jl,js)+dpw(jl,js))*y(jl,js)
      END DO
    END IF
  END DO
END IF

!       5.  Jacobian elements in final form
!           -------- -------- -- ----- ----

DO jtr = 1, jpcspf
  DO jl = 1, n_points
    IF ( f(jl,jtr)  >   peps ) THEN
      ej(jl,jtr) = ej(jl,jtr) / f(jl,jtr)
    ELSE
      ej(jl,jtr) = 0.0
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_jac
END MODULE asad_jac_mod
