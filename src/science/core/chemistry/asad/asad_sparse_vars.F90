! *****************************COPYRIGHT*******************************
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution.
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module holding variables and routines for sparse algebra.
!    Part of the ASAD chemical solver. Contains the following
!    routines:
!      setup_spfuljac
!      spfuljac
!      splinslv2
!      spresolv2
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
!     Method
!     ------
!     Sparse algebra works in the same way as dense algebra, namely by LU
!     decomposition (Gaussian elimination) of the linear system. In contrast
!     to dense algebra, here we keep track of non-zero matrix elements
!     and hence cut out algebraic manipulations whose result would be zero.
!     To make this efficient, species need to be reorder, such that those
!     with the fewest non-zero matrix elements (reactions) associated with
!     them, occur first in the list, and those with most (generally OH)
!     occur last.
!
! ######################################################################

MODULE asad_sparse_vars

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ASAD_SPARSE_VARS'

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE setup_spfuljac()
!
! This routine is divided into 4 parts.
!
! 1. Create nonzero_map(1:jpcspf,1:jpcspf) and base_tracer(1:spfjsize_max)
!
!    nonzero_map(1:jpcspf,1:jpcspf) is an integer indexing array such that, for
!    a given species index jl, spfj(jl,nonzero_map(i,j)) = A(i,j) where A is
!    the dense Jacobian matrix for species jl and spfj is the compressed
!    (sparse) storage of the dense Jacobian matrix.
!
! 2. Create permuted_nonzero_map(1:jpcspf,1:jpcspf)
!
!    permuted_nonzero_map(1:jpcspf,1:jpcspf) is an integer indexing array such
!    that spfj(jl,permuted_nonzero_map(i,j)) = (PAP')(i,j) where (PAP')(i,j) is
!    the (i,j) entry of the matrix matmul(P,matmul(A,transpose(P))) where P is
!    the permutation matrix and A is the dense Jacobian.
!
! 3. Check number of nonzero entries in LU factorization of the dense Jacobian
!
!    Gaussian elimination is used to solve a linear equation in splinslv2 and
!    the array spfj is used to hold the LU factorization. This section checks
!    that when this LU factorization is done the spfj array is sufficiently
!    large enough to hold the nonzero entries of the LU factorization.
!
! 4. Calculate product and loss indexing arrays for Jacobian
!
!    posterms, negterms and fracterms are used in calculating the values for
!    the Jacobian matrix.

USE asad_mod, ONLY: specf, frpx, jpcspf, jpfrpx, jpmsp, jpspec,                &
                    madvtr, modified_map, ndepd, ndepw, nfrpx, njcoth, nltrf,  &
                    nmsjac, nmzjac, nonzero_map, permuted_nonzero_map,         &
                    npdfr, nsjac1, nstst, ntabpd, ntrf, ntro3, nzjac1,         &
                    reorder, spfjsize_max, maxterms, maxfterms,                &
                    nposterms, nnegterms, nfracterms, posterms, negterms,      &
                    fracterms, base_tracer, ffrac, ztabpd, total
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY:                                                          &
    PrintStatus, PrStatus_Oper, PrStatus_Normal, umMessage, umPrint
USE errormessagelength_mod, ONLY: errormessagelength
USE ukca_um_legacy_mod, ONLY: mype

IMPLICIT NONE

! Local variables

INTEGER :: total1            ! total number of nonzero entries in Jacobian
INTEGER :: errcode           ! Variable passed to ereport
INTEGER :: errcodes(spfjsize_max,3)  ! An array for recording error codes

INTEGER :: irj
INTEGER :: itrd
INTEGER :: i
INTEGER :: j
INTEGER :: jc
INTEGER :: itrcr
INTEGER :: j3
INTEGER :: jn
INTEGER :: js
INTEGER :: i1
INTEGER :: i2
INTEGER :: isx
INTEGER :: kr
INTEGER :: ikr
INTEGER :: krj
INTEGER :: ij
INTEGER :: itemp1
INTEGER :: activity(jpcspf)

INTEGER, ALLOCATABLE :: permute(:,:)  ! permutation matrix
INTEGER, ALLOCATABLE :: map(:,:)      ! a matrix of 1's indicating where the
                                      ! dense Jacobian A(i,j) has nonzero
                                      ! entries

CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETUP_SPFULJAC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ------------------------------------------------------------------------------
! Section 1: Create nonzero_map (index array for Jacobian matrix) and
!            base_tracer
! ------------------------------------------------------------------------------

IF (.NOT. ALLOCATED(map)) ALLOCATE(map(jpcspf, jpcspf))

map(:,:) = 0
DO i = 1, jpcspf
  map(i,i) = 1
END DO

DO jc = 1, ntrf
  itrcr = nltrf(jc)
  DO j3 = 1, nmzjac(itrcr)
    irj = nzjac1(j3,itrcr)
    !
    DO jn = 1, jpmsp
      ij = njcoth(irj,jn)
      IF (ij /= 0) map(ij,itrcr) = 1
    END DO
    !
    i1 = npdfr(irj,1)
    IF (i1 /= 0) THEN
      DO jn = i1, npdfr(irj,2)
        isx = ntabpd(jn,1)
        map(isx,itrcr) = 1
      END DO
    END IF
    !
  END DO
END DO
!
!  Go through the steady state additions to the Jacobian; currently
!  assume that only O(1D) and O(3P) are modelled, and that both are
!  required for the O3 loss rate.
!
DO jc = 1, nstst
  DO j3 = 1, nmsjac(jc)
    irj = nsjac1(j3,jc)
    DO jn = 1, jpmsp
      ij = njcoth(irj,jn)
      IF (ij /= 0) map(ij,ntro3) = 1
    END DO
  END DO
END DO
!
!     -----------------------------------------------------------------
!              Add deposition terms to Jacobian diagonal.
!     -----------------------------------------------------------------
!
IF ((ndepw /= 0) .OR. (ndepd /= 0) ) THEN
  DO js = 1, jpspec
    itrd = madvtr(js)
    IF (itrd /= 0) map(itrd,itrd) = 1
  END DO
END IF

total = SUM(map)
IF (mype == 0 .AND. printstatus >= prstatus_oper) THEN
  WRITE(umMessage,'(A,I6)')                                                    &
  'TOTAL NUMBER OF NONZERO ENTRIES IN JACOBIAN: ',total
  CALL umPrint(umMessage,src='asad_sparse_vars')
END IF

isx = 0
nonzero_map(:,:) = 0
DO i = 1, jpcspf
  DO j = 1, jpcspf
    ! calculate forward and backward pointers for sparse representation
    IF (map(i,j) == 1) THEN
      isx = isx + 1
      ! backward pointer
      nonzero_map(i,j) = isx
    END IF
  END DO
END DO

base_tracer(:) = 0
DO j = 1, jpcspf
  DO i = 1, jpcspf
    i1 = nonzero_map(i,j)
    IF (i1 > 0) base_tracer(i1) = j
  END DO
END DO

! ------------------------------------------------------------------------------
! Section 2: Create permuted_nonzero_map (index array for PAP')
! ------------------------------------------------------------------------------

! Reorder species by their reactivity to minimize fill-in
DO i = 1, jpcspf
  activity(i) = SUM(map(:,i)) + SUM(map(i,:))
  reorder(i) = i
END DO

DO i = 1, jpcspf-1
  DO j = i+1, jpcspf
    IF (activity(i) > activity(j)) THEN
      ! exchange i and j tracers if i is more active than j.
      itemp1 = reorder(i)
      reorder(i) = reorder(j)
      reorder(j) = itemp1
      itemp1 = activity(i)
      activity(i) = activity(j)
      activity(j) = itemp1
    END IF
  END DO
END DO

DO i = 1, jpcspf
  IF (mype == 0 .AND. printstatus >= prstatus_oper) THEN
    WRITE(umMessage,'(A,I3,1X,I3,1X,A10,1X,I3)') 'IN ASAD SETUP_SPFULJAC: ',   &
         i, reorder(i), specf(reorder(i)), activity(i)
    CALL umPrint(umMessage,src='asad_sparse_vars')
  END IF
END DO

! Determine permutation matrix P
IF (.NOT. ALLOCATED(permute)) ALLOCATE(permute(jpcspf,jpcspf))
permute(:,:) = 0
DO i = 1, jpcspf
  permute(i,reorder(i)) = 1
END DO

! Calculate the index array permuted_nonzero_map such that
!   spfj(jl,permuted_nonzero_map(i,j)) = P*A*P'(i,j)
! for each chemical species jl where P is the permutation matrix and A is the
! dense Jacobian matrix with spfj the compressed storage (sparse) array for
! matrix A.
permuted_nonzero_map = MATMUL(MATMUL(permute, nonzero_map), TRANSPOSE(permute))

IF (ALLOCATED(permute)) DEALLOCATE(permute)
IF (ALLOCATED(map))     DEALLOCATE(map)

! ------------------------------------------------------------------------------
! Section 3: Check the number of nonzero elements in LU factorization
! ------------------------------------------------------------------------------

! Calculate the number of nonzero matrix elements in the LU factorization of the
! array PAP' and check that it is less than spfjsize_max.
total1 = total
modified_map(:,:) = permuted_nonzero_map(:,:)
DO kr = 1, jpcspf
  DO i = kr+1, jpcspf
    ikr = modified_map(i,kr)
    IF (ikr > 0) THEN
      DO j = kr+1, jpcspf
        krj = modified_map(kr,j)
        IF (krj > 0) THEN
          ij = modified_map(i,j)
          ! If the matrix entry is zero then create new matrix element.
          IF (ij <= 0) THEN
            total1 = total1 + 1
            modified_map(i,j) = total1
          END IF
        END IF
      END DO
    END IF
  END DO
END DO

! Perform error check outside of the loop to better suit GPU runs
IF (total1 > spfjsize_max) THEN
  errcode = total1
  WRITE(umMessage,'(A,2I4)') 'Total1 exceeded spfjsize_max: ', total1,         &
    spfjsize_max
  CALL umPrint(umMessage,src='asad_sparse_vars')
  CALL ereport('SETUP_SPFULJAC',errcode,'Increase spfjsize_max')
END IF

IF (mype == 0 .AND. printstatus >= prstatus_normal) THEN
  WRITE(umMessage,'(A,2I4)') 'Initial and final fill-in: ',total, total1
  CALL umPrint(umMessage,src='asad_sparse_vars')
  WRITE(umMessage,'(A,I4)') 'spfjsize_max is currently set to: ',spfjsize_max
  CALL umPrint(umMessage,src='asad_sparse_vars')
  WRITE(umMessage,'(A)') 'Reduce spfjsize_max if much greater than total1 '
  CALL umPrint(umMessage,src='asad_sparse_vars')
END IF

! ------------------------------------------------------------------------------
! Section 4: Calculate production and loss terms
! ------------------------------------------------------------------------------

nposterms(:)   = 0
nnegterms(:)   = 0
nfracterms(:)  = 0
posterms(:,:)  = 0
negterms(:,:)  = 0
fracterms(:,:) = 0
ffrac(:,:)     = 0.0

! Note that posterms, negterms and fracterms are particularly large arrays,
! typically with size 1000 x 160, less than 0.5% of the array contain nonzero
! values so there is potential here to use compressed storage format.

errcodes(:,:) = 0
DO jc = 1, ntrf
  itrcr = nltrf(jc)
  !
  DO j3 = 1, nmzjac(itrcr)
    irj = nzjac1(j3,itrcr)
    !
    DO jn = 1, jpmsp
      ij = njcoth(irj,jn)
      IF (ij /= 0) THEN
        i = nonzero_map(ij,itrcr)
        IF (jn < 3) THEN
          nnegterms(i) = nnegterms(i) + 1
          IF (nnegterms(i) > maxterms) THEN
            errcodes(i,1) = 1
          ELSE IF (ALL(errcodes(i,:) == 0)) THEN
            negterms(i,nnegterms(i)) = irj
          END IF
        ELSE
          IF (nfrpx(irj) == 0) THEN
            nposterms(i) = nposterms(i) + 1
            IF (nposterms(i) > maxterms) THEN
              errcodes(i,1) = 2
            ELSE IF (ALL(errcodes(i,:) == 0)) THEN
              posterms(i,nposterms(i)) = irj
            END IF
          ELSE
            nfracterms(i) = nfracterms(i) + 1
            IF (nfracterms(i) > maxfterms) THEN
              errcodes(i,1) = 3
            END IF
            IF (nfrpx(irj)+jn-3 > jpfrpx) THEN
              errcodes(i,2) = i
              errcodes(i,3) = jn
            ELSE IF (ALL(errcodes(i,:) == 0)) THEN
              fracterms(i,nfracterms(i)) = irj
              ffrac(i,nfracterms(i))     = frpx(nfrpx(irj)+jn-3)
            END IF
          END IF
        END IF
      END IF
    END DO
    !
    IF (npdfr(irj,1) /= 0) THEN
      i1 = npdfr(irj,1)
      i2 = npdfr(irj,2)
      DO jn = i1, i2
        isx = ntabpd(jn,1)
        i = nonzero_map(isx,itrcr)
        nfracterms(i) = nfracterms(i) + 1
        IF (nfracterms(i) > maxfterms) THEN
          errcodes(i,1) = 3
        ELSE IF (.NOT. ANY(errcodes(i,:) /= 0)) THEN
          fracterms(i,nfracterms(i)) = irj
          ffrac(i,nfracterms(i))     = ztabpd(jn,1)
        END IF
      END DO
    END IF
    !
  END DO
END DO

! Perform error checks outside of the loop to better suit GPU runs
IF (ANY(errcodes(:,:) /= 0)) THEN
  IF (ANY(errcodes(:,1) /= 0)) THEN
    errcode = MINVAL(errcodes(:,1),mask=(errcodes(:,1) /= 0))
    cmessage = ' Increase maxterms'
    CALL ereport('SETUP_SPFULJAC',errcode,cmessage)
  ELSE
    errcode = 4
    i = MINVAL(errcodes(:,2),mask=(errcodes(:,2) /= 0))
    irj = negterms(i,nnegterms(i))
    jn = errcodes(i,3)
    cmessage = ' frpx array index > jpfrpx'
    WRITE(umMessage,'(A)') cmessage
    CALL umPrint(umMessage,src='asad_sparse_vars')
    WRITE(umMessage,'(A,I4,A,I4)')'irj: ',irj,' nfrpx(irj): ',nfrpx(irj)
    CALL umPrint(umMessage,src='asad_sparse_vars')
    WRITE(umMessage,'(A,I4,A,I4)') 'i: ',i,' jn: ',jn
    CALL umPrint(umMessage,src='asad_sparse_vars')
    CALL ereport('SETUP_SPFULJAC',errcode,cmessage)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE setup_spfuljac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE spfuljac(n_points, cdt, min_pivot, nonzero_map, spfj)
!
!  Routine to calculate the Jacobian in sparse format
!
USE asad_mod, ONLY: ctype, deriv, dpd, dpw, f, jpcspf, jpfm, jpif, jpmsp,      &
                    jpspec, linfam, madvtr, moffam, ndepd, ndepw, njcoth,      &
                    nmsjac, nodd, nsjac1, nstst, ntro3, prk, spfjsize_max,     &
                    nposterms, nnegterms, nfracterms, posterms, negterms,      &
                    fracterms, base_tracer, ffrac, y, total
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: n_points
REAL, INTENT(IN)    :: cdt
REAL, INTENT(IN)    :: min_pivot
INTEGER, INTENT(IN) :: nonzero_map(jpcspf,jpcspf)
REAL, INTENT(OUT)   :: spfj(n_points,spfjsize_max)

! Local variables
REAL :: deltt

CHARACTER (LEN=2) :: ityped

INTEGER :: p
INTEGER :: irj
INTEGER :: ifamd
INTEGER :: itrd
INTEGER :: i
INTEGER :: jc
INTEGER :: j3
INTEGER :: jn
INTEGER :: js
INTEGER :: jl
INTEGER :: ij

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SPFULJAC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! At the bottom of this routine we divide by f, so ensure that
! f is not too small.
DO i = 1, jpcspf
  DO jl = 1, n_points
    IF (f(jl, i) < min_pivot) THEN
      f(jl, i) = min_pivot
    END IF
  END DO
END DO

deltt = 1.0 / cdt

! Calculate diagonal element of Jacobian
spfj(:,:) = 0.0
DO i = 1, jpcspf
  spfj(:,nonzero_map(i,i)) = -deltt*f(:,i)
END DO
!
!
!     -----------------------------------------------------------------
!           2.  Calc. full Jacobian matrix.
!               ----- ---- -------- -------
!
! Sum up positive non-fractional, negative, and fractional terms

DO p = 1, total
  DO i = nposterms(p), 1, -1
    spfj(:,p) = spfj(:,p) + prk(:,posterms(p,i))
  END DO
  DO i = nnegterms(p), 1, -1
    spfj(:,p) = spfj(:,p) - prk(:,negterms(p,i))
  END DO
  DO i = nfracterms(p), 1, -1
    spfj(:,p) = spfj(:,p) + ffrac(p,i)*prk(:,fracterms(p,i))
  END DO
END DO

!  Go through the steady state additions to the Jacobian; currently
!  assume that only O(1D) [and O(3P)] are modelled, and that [both] are
!  required for the O3 loss rate.
!
DO jc = 1, nstst
  DO j3 = 1, nmsjac(jc)
    irj = nsjac1(j3,jc)
    DO jn = 1, jpmsp
      ij = njcoth(irj,jn)
      IF (ij /= 0) THEN
        p = nonzero_map(ij,ntro3)
        IF (jn < 3) THEN
          spfj(:,p) = spfj(:,p) - prk(:,irj)*deriv(:,jc,1)
        ELSE
          spfj(:,p) = spfj(:,p) + prk(:,irj)*deriv(:,jc,1)
        END IF
      END IF
    END DO
  END DO
END DO
!
!     -----------------------------------------------------------------
!          4.  Add deposition terms to Jacobian diagonal.
!              --- ---------- ----- -- -------- ---------
!
IF ((ndepw /= 0) .OR. (ndepd /= 0)) THEN
  DO js = 1, jpspec
    ifamd = moffam(js)
    itrd = madvtr(js)
    !
    IF (ifamd /= 0) THEN
      ityped = ctype(js)
      p = nonzero_map(ifamd,ifamd)
      DO jl = 1, n_points
        IF ((ityped == jpfm) .OR. ((ityped == jpif) .AND. linfam(jl,itrd))) THEN
          spfj(jl,p) = spfj(jl,p) - nodd(js)*(dpd(jl,js) + dpw(jl,js))*y(jl,js)
        END IF
      END DO
    END IF
    IF (itrd /= 0) THEN
      p = nonzero_map(itrd,itrd)
      spfj(:,p) = spfj(:,p) - (dpd(:,js) + dpw(:,js))*y(:,js)
    END IF
  END DO
END IF
!
!     -------------------------------------------------------------
!          5.  Jacobian elements in final form
!              -------- -------- -- ----- ----
!
DO p = 1, total
  spfj(:,p) = spfj(:,p) / f(:,base_tracer(p))   ! filter f earlier!
END DO

!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE spfuljac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE splinslv2(n_points, bb, xx, min_pivot, max_val,                     &
                     permuted_nonzero_map, modified_map, spfj)

USE asad_mod, ONLY: jpcspf, spfjsize_max, total
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE
!
! linear_solve solves the equation A x = b where A and b are known.
!
! In UKCA A x = b is re-written as solving P A P' z = P b where P' z = x
! Here P' denotes transpose(P) and P is the permutation matrix calculated in
! setup_jacobian()
! PAP' is written in terms of L (lower triangular) and U (upper triangular)
! factors such that P A P' = LU with L containing diagonal entries of 1's and
! U containing the pivot values.
! linear_solve is divided into 2 sections:
!
! Section 1: calculate L and U such that P A P' = L U
! Here P denotes the permutation matrix calculated in setup_jacobian()
! The L and U matrices are stored in the original sparse Jacobian matrix to
! reuse the array.
!
! Section 2: Solve L U z = P b with P' z = x
! This is done in four parts:
! Part (a) calculate right-hand side P b
! Part (b) forward-substitution: find w where L w = P b
! Part (c) back-substitution: find z where U z = w
! Part (d) determine x, apply transpose(P) to z, P'z = x

INTEGER, INTENT(IN)     :: n_points
REAL, INTENT(IN OUT)    :: bb(n_points,jpcspf)
REAL, INTENT(OUT)       :: xx(n_points,jpcspf)
REAL, INTENT(IN)        :: min_pivot
REAL, INTENT(IN)        :: max_val
INTEGER, INTENT(IN)     :: permuted_nonzero_map(jpcspf,jpcspf)
INTEGER, INTENT(IN OUT) :: modified_map(jpcspf,jpcspf)
REAL, INTENT(IN OUT)    :: spfj(1:n_points,1:spfjsize_max)

! Local variables
INTEGER :: kr
INTEGER :: jl
INTEGER :: i
INTEGER :: j
INTEGER :: ikr
INTEGER :: krj
INTEGER :: ij

REAL :: multiplier(n_points)
REAL :: pivot(n_points)
REAL :: lower(n_points)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SPLINSLV2'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! NOTE: splinslv2 is always preceded by spfuljac, which computes spfj to be a
!       sparse (compressed) representation of the Jacobian matrix for each
!       species.

! Filter sparse Jacobian
#if defined(IBM_XL_FORTRAN)
! Version optimised for IBM by using the fsel IBM-only intrinsic
DO j = 1, total
  DO jl = 1, n_points
    tmp = fsel(spfj(jl,j)+max_val, spfj(jl,j), -max_val)
    spfj(jl,j) = fsel(tmp-max_val, max_val, tmp)
  END DO
END DO
#else
DO j = 1, total
  DO jl = 1, n_points
    spfj(jl,j) = MIN(MAX(spfj(jl,j), -max_val), max_val)
  END DO
END DO
#endif


! Section 1: Determine L U factors such that L U = P A P'
!
! The L U factors are overwritten onto the original sparse Jacobian array. That
! is, for each species index jl, the upper triangular part of the matrix
! represented by spfj(jl,:) holds the U factor and the strictly lower triangular
! part holds the L factor below its unit diagonal.

! Loop over rows for pivoting
DO kr = 1, jpcspf

  ! The multipliers are the diagonal entries for each species
  multiplier(:) = spfj(:,modified_map(kr,kr))

  ! Determine the pivot for each species using a threshold
  WHERE (ABS(multiplier) > min_pivot)
    pivot(:) = 1.0 / multiplier
  ELSE WHERE
    pivot(:) = max_val
  END WHERE

  ! Loop over non-zero entries in the lower-diagonal part of column kr
  DO i = kr+1, jpcspf
    ikr = modified_map(i,kr)
    IF (ikr > 0) THEN

      ! Compute contribution to the lower diagonal part and store it in spfj
      lower(:) = spfj(:,ikr)*pivot
      spfj(:,ikr) = lower

      ! Loop over non-zero entries in the upper-diagonal part of row kr
      DO j = kr+1, jpcspf
        krj = modified_map(kr,j)
        IF (krj > 0) THEN

          ! Compute the contribution towards the upper triangular part and store
          ! it in the upper-diagonal part of spfj
          ! NOTE: We can rely on the entry being nonzero because the sparsity
          ! pattern was pre-computed.
          ij = modified_map(i,j)
          spfj(:,ij) = spfj(:,ij) - lower*spfj(:,krj)

        END IF
      END DO
    END IF
  END DO
END DO

! Filter sparse Jacobian
#if defined(IBM_XL_FORTRAN)
! Version optimised for IBM by using the fsel IBM-only intrinsic
DO j = 1, total
  DO jl = 1, n_points
    tmp = fsel(spfj(jl,j)+max_val, spfj(jl,j), -max_val)
    spfj(jl,j) = fsel(tmp-max_val, max_val, tmp)
  END DO
END DO
#else
DO j = 1, total
  DO jl = 1, n_points
    spfj(jl,j) = MIN(MAX(spfj(jl,j), -max_val), max_val)
  END DO
END DO
#endif


! Section 2: Solve P A P' z = P b with P'z = x using L U z = P b
CALL spresolv2(n_points, bb, xx, min_pivot, modified_map, spfj, max_val)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE splinslv2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE spresolv2(n_points, bb, xx, min_pivot, modified_map, spfj, max_val)

! This subroutine determines x where L U z = P b with P' z = x
! The L U factors are supplied to this routine and contained with spfj array.
!
! This is done in four parts:
! Part (a) calculate the right-hand side P b
! Part (b) forward-substitution: find w where L w = P b
! Part (c) back-substitution: find z where U z = w
! Part (d) determine x, apply transpose(P) to z, P'z = x

USE asad_mod, ONLY: jpcspf, reorder, spfjsize_max
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN)  :: n_points
REAL,    INTENT(IN)  :: bb(n_points,jpcspf)
REAL,    INTENT(OUT) :: xx(n_points,jpcspf)
REAL,    INTENT(IN)  :: min_pivot
INTEGER, INTENT(IN)  :: modified_map(jpcspf,jpcspf)
REAL,    INTENT(IN)  :: spfj(n_points,spfjsize_max)

! Maximum tolerated value for use in filtering step. Unused if negative
REAL,    INTENT(IN)  :: max_val


INTEGER :: kr
INTEGER :: i
INTEGER :: j
INTEGER :: krj
INTEGER :: ikr

REAL :: bb1(n_points, jpcspf)
REAL :: xx1(n_points, jpcspf)
REAL :: pivot(n_points)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SPRESOLV2'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Part (a) calculate right-hand side P b
DO i = 1, jpcspf
  bb1(:,i) = bb(:,reorder(i))
END DO

! Part (b) forward-substitution: find w where L w = P b
DO kr = 1, jpcspf-1
  DO i = kr+1, jpcspf
    ikr = modified_map(i,kr)
    IF (ikr > 0) bb1(:,i) = bb1(:,i) - spfj(:,ikr)*bb1(:,kr)
  END DO
END DO

! Part (c) back-substitution: find z where U z = w
DO kr = jpcspf, 1, -1
  xx1(:,kr) = bb1(:,kr)
  DO j = kr+1, jpcspf
    krj = modified_map(kr,j)
    IF (krj > 0) xx1(:,kr) = xx1(:,kr) - spfj(:,krj)*xx1(:,j)
  END DO
  pivot(:) = spfj(:,modified_map(kr,kr))
  WHERE (ABS(pivot) < min_pivot) pivot = min_pivot
  IF (max_val > 0) THEN
    xx1(:,kr) = MIN(MAX(xx1(:,kr), -max_val), max_val) / pivot
  ELSE
    xx1(:,kr) = xx1(:,kr) / pivot
  END IF
END DO

! Part (d) determine x, apply transpose(P) to z, P'z = x
xx(:,:) = 0.0
DO i = 1, jpcspf
  xx(:,reorder(i)) = xx1(:,i)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE spresolv2

END MODULE asad_sparse_vars
