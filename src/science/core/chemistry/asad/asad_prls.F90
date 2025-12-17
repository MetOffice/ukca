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
! Purpose: To compute the rate of change from each reaction and to
!     compute the production and loss for each species.
!
!     This routine can be regarded as the main engine to ASAD and
!     accounts for a significant part of the total computational
!     effort.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_DIFFUN and ASAD_FTOY
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     Interface
!     ---------
!     Called by ftoy and diffun.
!
!        kl       - No. of spatial points to compute for. Upper limit
!                   to inner loop over gridpoints. Usually this is set
!                   to jpnl. However, some stiff integrators only comput
!                   a single gridpoint at a time, in which case, this
!                   is set to 1. ***NOTE**** Setting kl to 1 on a vector
!                   computer will severely hamper performance and
!                   effectively cause this routine (which accounts for
!                   most of the computation in ASAD) to run in scalar
!                   mode.
!        knspec   - Number of species to compute product and loss terms
!                   for. Normally, this routine is called to compute
!                   prod. & loss terms for steady state species or non-
!                   steady state species (rather than all the species).
!        kspec    - Array of size knspec which holds the species number
!                   to compute the product and loss terms for.
!        ldepem   - If true, deposition and emission terms are considere
!                   although the user may not have turned them on for
!                   individual species. If .false., then dry & wet
!                   deposition and emissions are ignored regardless of
!                   what the user has set. This option is useful when
!                   the product and loss are required without these
!                   effects.
!
!     On exit from this routine, the arrays; prod, slos and prk will hav
!     been assigned.
!
!     Method
!     ------
!
!     First the reaction rates of change are computed. These are
!     divided into reactions with one reactant and those with two
!     reactants.
!
!     If we write dy/dt = P - Ly then prod holds P and slos holds Ly,
!     so that prod and slos are always positive. The same code is used
!     for both prod and slos which means we equivalence these arrays
!     to one which spans both. This allows the summation of reaction
!     rates to be unrolled so that we can compute the sums in groups of
!     three. This has proved to be much more efficient on the Cray
!     computer architecture than summing the terms one by one. Tests hav
!     shown there is little to be gained by unrolling over more than 3
!     terms. Therefore when indexing the pd array we loop over it twice
!     with an interval of jpspec. The first pass updates the products an
!     the second updates the losses.  [Equivalence replaced by creating
!     prod and slos as pointers in UM code - see routine asad_mod_init.]
!
!     For some reactions, the products may be fractional e.g. the actual
!     rate of production is a percentage of the production given by the
!     rate of reaction. All fractional products are treated separated as
!     they involve another multiplication. This loop is not unrolled as
!     it is assumed that the no. of reactions with fractional products
!     is small compared to the total number of reactions.
!
!     The individual species loss rates are updated by adding the
!     terms due to wet & dry deposition, using previously calculated
!     rates. Similarly, the production rates are updated by adding
!     previously calculated emission rates. Note that wet and dry
!     deposition are treated are first order loss rates
!     ie. dy/dt = -Wy - Dy where W is the
!     wet and D is the dry deposition rates respectively.
!
!     Local variables
!     ---------------
!     ip1, ip2, ip3   : no. of groups of 1, 2 & 3 term sums for each
!                       species
!     ngrp  : for each species, holds the no. of 3 sum, 2 sum and single
!             sum terms and therefore controls the loop limits
!     nprdx3, nprdx2, nprdx1  : arrays holds the species indicies for
!                               3 term, 2 term and single term sums.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_prls_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_PRLS_MOD'

CONTAINS

SUBROUTINE asad_prls( kl, knspec, kspec, ldepem )

USE asad_mod,            ONLY: prod, slos, y, rk, prk,                         &
                               dpd, dpw, pd,                                   &
                               nuni, nspi, ngrp, nprdx1,                       &
                               nprdx2, nprdx3, ntabfp, frpx,                   &
                               nldepx, nnfrp, jpspec, jpnr
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) :: kl
INTEGER, INTENT(IN) :: knspec
INTEGER, INTENT(IN) :: kspec(knspec)

LOGICAL, INTENT(IN) ::  ldepem

!       Local variables

INTEGER :: j        ! Loop variable
INTEGER :: j1       ! Loop variable
INTEGER :: j2       ! Loop variable
INTEGER :: j3       ! Loop variable
INTEGER :: jl       ! Loop variable
INTEGER :: jr       ! Loop variable
INTEGER :: js       ! Loop variable
INTEGER :: jskip    ! Loop variable
INTEGER :: i1       ! Index
INTEGER :: i2       ! Index
INTEGER :: i3       ! Index
INTEGER :: ir1
INTEGER :: ir2
INTEGER :: ip1
INTEGER :: ip2
INTEGER :: ip3
INTEGER :: iunip1
INTEGER :: iss
INTEGER :: ir
INTEGER :: ix       ! index for frpx array
INTEGER :: istart
INTEGER :: iend

REAL :: fr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_PRLS'


!       1.  Initialise variables.
!           ---------- ----------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO js = 1, jpspec
  DO jl = 1, kl
    prod(jl,js) = 0.0
    slos(jl,js) = 0.0
  END DO
END DO

!       2.  Compute reaction rate of change (kA or kAB).
!           ------- -------- ---- -- ------ --- -- -----

!       2.1  Reactions with one reactant.


DO jr = 1, nuni
  ir1 = nspi(jr,1)

  DO jl = 1, kl
    prk(jl,jr) = rk(jl,jr) * y(jl,ir1)
  END DO
END DO

!       2.2  Reactions with two reactants.

iunip1 = nuni + 1
DO jr = iunip1, jpnr
  ir1 = nspi(jr,1)
  ir2 = nspi(jr,2)
  DO jl = 1, kl
    prk(jl,jr) = rk(jl,jr) * y(jl,ir1) * y(jl,ir2)
  END DO
END DO

!       3.  Compute production & loss term; see method for details.
!           ------- ---------- - ---- ----- --- ------ --- --------

DO jskip = 0, jpspec, jpspec

  DO j = 1, knspec
    js  = kspec(j) + jskip
    ip3 = ngrp(js,3)
    ip2 = ngrp(js,2)
    ip1 = ngrp(js,1)

    !           3.1  Terms in groups of 3

    DO j3 = 1, ip3
      i1 = nprdx3(1,j3,js)
      i2 = nprdx3(2,j3,js)
      i3 = nprdx3(3,j3,js)

      DO jl = 1, kl
        pd(jl,js) = pd(jl,js) + prk(jl,i1) + prk(jl,i2) + prk(jl,i3)
      END DO
    END DO

    !           3.2  Terms in groups of 2; this loop will only ever be
    !                ip2=0 or 1.

    DO j2 = 1, ip2
      i1 = nprdx2(1,js)
      i2 = nprdx2(2,js)
      DO jl = 1, kl
        pd(jl,js) = pd(jl,js) + prk(jl,i1) + prk(jl,i2)
      END DO
    END DO

    !           3.3  Single term left over; ip1 will only be 0 or 1.

    DO j1 = 1, ip1
      i1 = nprdx1(js)
      DO jl = 1, kl
        pd(jl,js) = pd(jl,js) + prk(jl,i1)
      END DO
    END DO

  END DO
END DO

!       3.4  Sum contribution from fractional products.

DO jr = 1, nnfrp
  iss = ntabfp(jr,1)
  ir = ntabfp(jr,2)
  ix = ntabfp(jr,3)
  fr = frpx(ix)
  DO jl = 1, kl
    pd(jl,iss) = pd(jl,iss) + fr * prk(jl,ir)
  END DO
END DO

!       4.  Add contribution from depositions and emissions.
!           --- ------------ ---- ----------- --- ----------

IF ( ldepem ) THEN

  !         4.1  Loss due to both dry and wet deposition.

  istart = nldepx(1)
  iend   = nldepx(2)
  DO j = istart, iend
    js = nldepx(j)
    DO  jl = 1, kl
      slos(jl,js) = slos(jl,js) + ( dpd(jl,js) + dpw(jl,js) ) * y(jl,js)
    END DO
  END DO

  !         4.2  Loss due to dry deposition only.

  istart = nldepx(3)
  iend   = nldepx(4)
  DO j = istart, iend
    js = nldepx(j)
    DO jl = 1, kl
      slos(jl,js) = slos(jl,js) + dpd(jl,js) * y(jl,js)
    END DO
  END DO

  !         4.3  Loss due to wet deposition

  istart = nldepx(5)
  iend   = nldepx(6)
  DO j = istart, iend
    js = nldepx(j)
    DO jl = 1, kl
      slos(jl,js) = slos(jl,js) + dpw(jl,js) * y(jl,js)
    END DO
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_prls
END MODULE asad_prls_mod
