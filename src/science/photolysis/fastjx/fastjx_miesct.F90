! *****************************COPYRIGHT*******************************
! (c) [University of California] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
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
!    Fast-JX routine for calculating online photolysis rates
!
!-----------------------------------------------------------------------
!   this is an adaption of the prather radiative transfer code, (mjp, 10/95)
!     Prather, 1974, astrophys. j. 192, 787-792.
!         sol_n of inhomogeneous rayleigh scattering atmosphere.
!         (original rayleigh w/ polarization)
!     Cochran and Trafton, 1978, ap.j., 219, 756-762.
!         raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, j.geophys.res., 94, 12975-13002.
!         chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!
!    takes atmospheric structure and source terms from std j-code
!    also limited to 4 gauss points, only calculates mean field! (m=1)
!
!    subroutine contains all other core scattering routines
!-----------------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
!-----------------------------------------------------------------------

MODULE fastjx_miesct_mod

IMPLICIT NONE

! Default private
PRIVATE

! Subroutines available outside this module
PUBLIC :: fastjx_miesct

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FASTJX_MIESCT_MOD'

CONTAINS

SUBROUTINE fastjx_miesct(fj,fjt,fjb, pomega,fz,ztau,zflux,rfl,u0,nd)

USE fastjx_data, ONLY: w_, m_, m2_, n_, emu, wt
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN)  :: nd                 ! no of depth pts for scattering
                                           ! calculation
REAL,    INTENT(IN)  :: pomega(m2_,n_,w_)  ! scattering phase fs on levels
                                           ! used for scattering calcs
REAL,    INTENT(IN)  :: fz(n_,w_)          ! attenuated incoming radiation
REAL,    INTENT(IN)  :: ztau(n_,w_)        ! optical depth on levels used for
                                           ! scattering
REAL,    INTENT(IN)  :: rfl(w_)            ! surface albedo
REAL,    INTENT(IN)  :: u0                 ! cosine (sza)
REAL,    INTENT(IN)  :: zflux(w_)          ! upward flux from Earth's surface

REAL,    INTENT(OUT) :: fj(n_,w_)          ! total 'scattered' radiation flux
REAL,    INTENT(OUT) :: fjt(w_)            ! 'upward' bound 'scattered'
                                           ! radiation flux
REAL,    INTENT(OUT) :: fjb(w_)            ! 'downward' bound 'scattered'
                                           ! radiation flux


! End of I/O

REAL :: pm(m_,m2_)
REAL :: pm0(m2_)

INTEGER :: i, im                            ! Loop variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_MIESCT'

! End of Header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i = 1,m_
  CALL legnd0 (emu(i),pm0,m2_)
  DO im = 1,m2_
    pm(i,im) = pm0(im)
  END DO
END DO

CALL legnd0 (-u0,pm0,m2_)
DO im=1,m2_
  pm0(im) = 0.25e0*pm0(im)
END DO

! ---blkslv now CALLed with all the wavelength arrays (k=1:w_)
CALL blkslv(fj,pomega,fz,ztau,zflux,rfl,pm,pm0,fjt,fjb, nd)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

!###########################################################################

CONTAINS

  !###########################################################################
  !-----------------------------------------------------------------------
  !---calculates ordinary legENDre fns of x (real)
  !---   from p[0] = pl(1) = 1,  p[1] = x, .... p[n-1] = pl(n)
SUBROUTINE legnd0 (x,pl,n)

IMPLICIT NONE

INTEGER, INTENT(IN)  :: n
REAL,    INTENT(IN)  :: x
REAL,    INTENT(OUT) :: pl(n)

INTEGER              :: i
REAL                 :: den
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LEGND0'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---always does pl(2) = p[1]
pl(1) = 1.0e0
pl(2) = x

DO i = 3,n
  den = (i-1)
  pl(i) = pl(i-1)*x*(2.0e0-1.0/den) - pl(i-2)*(1.0e0-1.0e0/den)
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE legnd0

!###########################################################################

!-----------------------------------------------------------------------
!  sets up and solves the block tri-diagonal system:
!               a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = h(i)
!  this goes back to the old, dumb, fast version 5.3
!-----------------------------------------------------------------------
SUBROUTINE blkslv                                                              &
     (fj,pomega,fz,ztau,zflux,rfl,pm,pm0,fjtop,fjbot,nd)

IMPLICIT NONE

INTEGER, INTENT(IN)  :: nd
REAL,    INTENT(IN)  :: pomega(m2_,n_,w_)
REAL,    INTENT(IN)  :: fz(n_,w_)
REAL,    INTENT(IN)  :: ztau(n_,w_)
REAL,    INTENT(IN)  :: pm(m_,m2_)
REAL,    INTENT(IN)  :: pm0(m2_)
REAL,    INTENT(IN)  :: rfl(w_)
REAL,    INTENT(IN)  :: zflux(w_)

REAL,    INTENT(OUT) :: fj(n_,w_)
REAL,    INTENT(OUT) :: fjtop(w_)
REAL,    INTENT(OUT) :: fjbot(w_)
! End of I/O

REAL                 :: a(m_,n_,w_)
REAL                 :: c(m_,n_,w_)
REAL                 :: h(m_,n_,w_)
REAL                 :: rr(m_,n_,w_)
REAL                 :: b(m_,m_,n_,w_)
REAL                 :: aa(m_,m_,n_,w_)
REAL                 :: cc(m_,m_,n_,w_)
REAL                 :: dd(m_,m_,n_,w_)
REAL                 :: e(m_,m_,w_)
REAL                 :: f(m_,m_,w_)

REAL                 :: sumb
REAL                 :: sumbx
REAL                 :: sumt

INTEGER              :: i, j, k, l
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BLKSLV'

! End of Header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO k=1,w_
  CALL gen_id (pomega(1,1,k),fz(1,k),ztau(1,k),zflux(k),rfl(k),                &
       pm,pm0, b(1,1,1,k),cc(1,1,1,k),aa(1,1,1,k),                             &
       a(1,1,k),h(1,1,k),c(1,1,k), nd)
END DO

!-----------upper boundary l=1
DO k = 1,w_
  DO j = 1,m_
    DO i = 1,m_
      f(i,j,k) = b(i,j,1,k)
    END DO
  END DO
END DO

CALL matinw (f,e)

DO k = 1,w_
  DO j = 1,m_
    DO i = 1,m_
      dd(i,j,1,k) = -e(i,1,k)*cc(1,j,1,k)-e(i,2,k)*cc(2,j,1,k)                 &
                    -e(i,3,k)*cc(3,j,1,k)-e(i,4,k)*cc(4,j,1,k)
    END DO
    rr(j,1,k) = e(j,1,k)*h(1,1,k)+e(j,2,k)*h(2,1,k)                            &
                + e(j,3,k)*h(3,1,k)+e(j,4,k)*h(4,1,k)
  END DO
END DO

! --------continue through all depth points id=2 to id=nd-1
DO l = 2,nd-1
  DO k = 1,w_
    DO j = 1,m_
      DO i = 1,m_
        b(i,j,l,k) = b(i,j,l,k) + a(i,l,k)*dd(i,j,l-1,k)
      END DO
      h(j,l,k) = h(j,l,k) - a(j,l,k)*rr(j,l-1,k)
    END DO

    DO j = 1,m_
      DO i = 1,m_
        f(i,j,k) = b(i,j,l,k)
      END DO
    END DO
  END DO

  CALL matinw (f,e)

  DO k = 1,w_
    DO j = 1,m_
      DO i = 1,m_
        dd(i,j,l,k) = - e(i,j,k)*c(j,l,k)
      END DO
      rr(j,l,k) = e(j,1,k)*h(1,l,k)+e(j,2,k)*h(2,l,k)                          &
              + e(j,3,k)*h(3,l,k)+e(j,4,k)*h(4,l,k)
    END DO
  END DO
END DO     ! l

!---------final depth point: l=nd
l = nd

DO k = 1,w_
  DO j = 1,m_
    DO i = 1,m_
      b(i,j,l,k) = b(i,j,l,k)                                                  &
              + aa(i,1,l,k)*dd(1,j,l-1,k)                                      &
              + aa(i,2,l,k)*dd(2,j,l-1,k)                                      &
              + aa(i,3,l,k)*dd(3,j,l-1,k)                                      &
              + aa(i,4,l,k)*dd(4,j,l-1,k)
    END DO
    h(j,l,k) = h(j,l,k)                                                        &
           - aa(j,1,l,k)*rr(1,l-1,k) - aa(j,2,l,k)*rr(2,l-1,k)                 &
           - aa(j,3,l,k)*rr(3,l-1,k) - aa(j,4,l,k)*rr(4,l-1,k)
  END DO
  DO j = 1,m_
    DO i = 1,m_
      f(i,j,k) = b(i,j,l,k)
    END DO
  END DO
END DO

CALL matinw (f,e)

DO k = 1,w_
  DO j = 1,m_
    rr(j,l,k) = e(j,1,k)*h(1,l,k)+e(j,2,k)*h(2,l,k)                            &
           + e(j,3,k)*h(3,l,k)+e(j,4,k)*h(4,l,k)
  END DO
END DO

! ----------back solution
DO l = nd-1,1,-1
  DO k = 1,w_
    DO j = 1,m_
      rr(j,l,k) = rr(j,l,k)                                                    &
           + dd(j,1,l,k)*rr(1,l+1,k)                                           &
           + dd(j,2,l,k)*rr(2,l+1,k)                                           &
           + dd(j,3,l,k)*rr(3,l+1,k)                                           &
           + dd(j,4,l,k)*rr(4,l+1,k)
    END DO
  END DO
END DO

! ---------mean j & h
fj(:,:) = 0.0e0
DO l = 1,nd,2
  DO k = 1,w_
    fj(l,k) = rr(1,l,k)*wt(1) + rr(2,l,k)*wt(2)                                &
         + rr(3,l,k)*wt(3) + rr(4,l,k)*wt(4)
  END DO
END DO

DO l = 2,nd,2
  DO k = 1,w_
    fj(l,k) = rr(1,l,k)*wt(1)*emu(1) + rr(2,l,k)*wt(2)*emu(2)                  &
         + rr(3,l,k)*wt(3)*emu(3) + rr(4,l,k)*wt(4)*emu(4)
  END DO
END DO

!---fjtop = scaled diffuse flux out top-of-atmosphere (limit = mu0)
!---fjbot = scaled diffuse flux onto surface:
!---zflux = reflect/(1 + reflect) * mu0 * fsolar(lower boundary)
!---sumbx = flux from lambert reflected i+
DO k = 1,w_
  sumt = rr(1, 1,k)*wt(1)*emu(1) + rr(2, 1,k)*wt(2)*emu(2)                     &
       + rr(3, 1,k)*wt(3)*emu(3) + rr(4, 1,k)*wt(4)*emu(4)
  sumb = rr(1,nd,k)*wt(1)*emu(1) + rr(2,nd,k)*wt(2)*emu(2)                     &
       + rr(3,nd,k)*wt(3)*emu(3) + rr(4,nd,k)*wt(4)*emu(4)
  sumbx = 4.0e0*sumb*rfl(k)/(1.0e0 + rfl(k)) + zflux(k)

  fjtop(k) = 4.0e0*sumt
  fjbot(k) = 4.0e0*sumb - sumbx
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE blkslv

!####################################################################
!
!  generates coefficient matrices for the block tri-diagonal system:
!               a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = h(i)
!
!-----------------------------------------------------------------------
SUBROUTINE gen_id(pomega,fz,ztau,zflux,rfl,pm,pm0                              &
     ,b,cc,aa,a,h,c,  nd)

IMPLICIT NONE

INTEGER, INTENT(IN)  :: nd
REAL,    INTENT(IN)  :: pomega(m2_,n_)
REAL,    INTENT(IN)  :: pm(m_,m2_)
REAL,    INTENT(IN)  :: pm0(m2_)
REAL,    INTENT(IN)  :: rfl
REAL,    INTENT(IN)  :: zflux
REAL,    INTENT(IN)  :: fz(n_)
REAL,    INTENT(IN)  :: ztau(n_)

REAL,    INTENT(OUT) :: b  (m_,m_,n_)
REAL,    INTENT(OUT) :: aa (m_,m_,n_)
REAL,    INTENT(OUT) :: cc (m_,m_,n_)
REAL,    INTENT(OUT) :: a  (m_,n_)
REAL,    INTENT(OUT) :: c  (m_,n_)
REAL,    INTENT(OUT) :: h  (m_,n_)

! End of I/O
REAL :: s(m_,m_)
REAL :: t(m_,m_)
REAL :: u(m_,m_)
REAL :: v(m_,m_)
REAL :: w(m_,m_)

REAL :: sum0
REAL :: sum1
REAL :: sum2
REAL :: sum3

REAL :: deltau
REAL :: d1
REAL :: d2
REAL :: surfac

INTEGER :: i, j, l1,l2,ll
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GEN_ID'


!---------------------------------------------
! End of Header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------upper boundary:  2nd-order terms
l1 = 1
l2 = 2

DO i = 1,m_
  sum0 =                                                                       &
       pomega(1,l1)*pm(i,1)*pm0(1) + pomega(3,l1)*pm(i,3)*pm0(3)               &
       + pomega(5,l1)*pm(i,5)*pm0(5) + pomega(7,l1)*pm(i,7)*pm0(7)

  sum2 =                                                                       &
       pomega(1,l2)*pm(i,1)*pm0(1) + pomega(3,l2)*pm(i,3)*pm0(3)               &
       + pomega(5,l2)*pm(i,5)*pm0(5) + pomega(7,l2)*pm(i,7)*pm0(7)

  sum1 =                                                                       &
       pomega(2,l1)*pm(i,2)*pm0(2) + pomega(4,l1)*pm(i,4)*pm0(4)               &
       + pomega(6,l1)*pm(i,6)*pm0(6) + pomega(8,l1)*pm(i,8)*pm0(8)

  sum3 =                                                                       &
       pomega(2,l2)*pm(i,2)*pm0(2) + pomega(4,l2)*pm(i,4)*pm0(4)               &
       + pomega(6,l2)*pm(i,6)*pm0(6) + pomega(8,l2)*pm(i,8)*pm0(8)

  h(i,l1) = 0.5e0*(sum0*fz(l1) + sum2*fz(l2))
  a(i,l1) = 0.5e0*(sum1*fz(l1) + sum3*fz(l2))

END DO

DO i = 1,m_
  DO j = 1,i
    sum0 =                                                                     &
         pomega(1,l1)*pm(i,1)*pm(j,1)                                          &
         + pomega(3,l1)*pm(i,3)*pm(j,3)                                        &
         + pomega(5,l1)*pm(i,5)*pm(j,5)                                        &
         + pomega(7,l1)*pm(i,7)*pm(j,7)

    sum2 =                                                                     &
         pomega(1,l2)*pm(i,1)*pm(j,1)                                          &
         + pomega(3,l2)*pm(i,3)*pm(j,3)                                        &
         + pomega(5,l2)*pm(i,5)*pm(j,5)                                        &
         + pomega(7,l2)*pm(i,7)*pm(j,7)

    sum1 =                                                                     &
         pomega(2,l1)*pm(i,2)*pm(j,2)                                          &
         + pomega(4,l1)*pm(i,4)*pm(j,4)                                        &
         + pomega(6,l1)*pm(i,6)*pm(j,6)                                        &
         + pomega(8,l1)*pm(i,8)*pm(j,8)

    sum3 =                                                                     &
         pomega(2,l2)*pm(i,2)*pm(j,2)                                          &
         + pomega(4,l2)*pm(i,4)*pm(j,4)                                        &
         + pomega(6,l2)*pm(i,6)*pm(j,6)                                        &
         + pomega(8,l2)*pm(i,8)*pm(j,8)

    s(i,j) = - sum2*wt(j)
    s(j,i) = - sum2*wt(i)
    t(i,j) = - sum1*wt(j)
    t(j,i) = - sum1*wt(i)
    v(i,j) = - sum3*wt(j)
    v(j,i) = - sum3*wt(i)
    b(i,j,l1) = - 0.5e0*(sum0 + sum2)*wt(j)
    b(j,i,l1) = - 0.5e0*(sum0 + sum2)*wt(i)
  END DO
END DO

DO i = 1,m_
  s(i,i)    = s(i,i)    + 1.0e0
  t(i,i)    = t(i,i)    + 1.0e0
  v(i,i)    = v(i,i)    + 1.0e0
  b(i,i,l1) = b(i,i,l1) + 1.0e0

  c(i,l1)= s(i,1)*a(1,l1)/emu(1) + s(i,2)*a(2,l1)/emu(2)                       &
       + s(i,3)*a(3,l1)/emu(3) + s(i,4)*a(4,l1)/emu(4)
END DO

DO i = 1,m_
  DO j = 1,m_

    w(j,i) = s(j,1)*t(1,i)/emu(1) + s(j,2)*t(2,i)/emu(2)                       &
         + s(j,3)*t(3,i)/emu(3) + s(j,4)*t(4,i)/emu(4)

    u(j,i) = s(j,1)*v(1,i)/emu(1) + s(j,2)*v(2,i)/emu(2)                       &
         + s(j,3)*v(3,i)/emu(3) + s(j,4)*v(4,i)/emu(4)
  END DO
END DO

!-------------upper boundary, 2nd-order, c-matrix is full (cc)
deltau = ztau(l2) - ztau(l1)
d2 = 0.25e0*deltau

DO i = 1,m_
  DO j = 1,m_
    b(i,j,l1) = b(i,j,l1) + d2*w(i,j)
    cc(i,j,l1) = d2*u(i,j)
  END DO
  h(i,l1) = h(i,l1) + 2.0e0*d2*c(i,l1)
  a(i,l1) = 0.0e0
END DO
DO i = 1,m_
  d1 = emu(i)/deltau
  b(i,i,l1)  = b(i,i,l1) + d1
  cc(i,i,l1) = cc(i,i,l1) - d1
END DO

!------------intermediate points:  can be even or odd, a & c diagonal
!---mid-layer h-points, legENDre terms 2,4,6,8
DO ll=2,nd-1,2
  deltau = ztau(ll+1) - ztau(ll-1)
  DO i = 1,m_
    a(i,ll) = emu(i)/deltau
    c(i,ll) = -a(i,ll)
    h(i,ll) = fz(ll)*(                                                         &
         pomega(2,ll)*pm(i,2)*pm0(2) + pomega(4,ll)*pm(i,4)*pm0(4)             &
         +pomega(6,ll)*pm(i,6)*pm0(6) +pomega(8,ll)*pm(i,8)*pm0(8))
  END DO

  DO i = 1,m_
    DO j=1,i
      sum0 =                                                                   &
           pomega(2,ll)*pm(i,2)*pm(j,2)                                        &
           + pomega(4,ll)*pm(i,4)*pm(j,4)                                      &
           + pomega(6,ll)*pm(i,6)*pm(j,6)                                      &
           + pomega(8,ll)*pm(i,8)*pm(j,8)

      b(i,j,ll) =  - sum0*wt(j)
      b(j,i,ll) =  - sum0*wt(i)
    END DO
  END DO
  DO i = 1,m_
    b(i,i,ll) = b(i,i,ll) + 1.0e0
  END DO
END DO

!---odd-layer j-points, legENDre terms 1,3,5,7
DO ll=3,nd-2,2
  deltau = ztau(ll+1) - ztau(ll-1)
  DO i = 1,m_
    a(i,ll) = emu(i)/deltau
    c(i,ll) = -a(i,ll)
    h(i,ll) = fz(ll)*(                                                         &
         pomega(1,ll)*pm(i,1)*pm0(1)                                           &
         + pomega(3,ll)*pm(i,3)*pm0(3)                                         &
         + pomega(5,ll)*pm(i,5)*pm0(5)                                         &
         + pomega(7,ll)*pm(i,7)*pm0(7))
  END DO

  DO i = 1,m_
    DO j=1,i
      sum0 =                                                                   &
           pomega(1,ll)*pm(i,1)*pm(j,1)                                        &
           + pomega(3,ll)*pm(i,3)*pm(j,3)                                      &
           + pomega(5,ll)*pm(i,5)*pm(j,5)                                      &
           + pomega(7,ll)*pm(i,7)*pm(j,7)

      b(i,j,ll) =  - sum0*wt(j)
      b(j,i,ll) =  - sum0*wt(i)
    END DO
  END DO

  DO i = 1,m_
    b(i,i,ll) = b(i,i,ll) + 1.0e0
  END DO
END DO

!---------lower boundary:  2nd-order terms
l1 = nd
l2 = nd-1

DO i = 1,m_
  sum0 =                                                                       &
       pomega(1,l1)*pm(i,1)*pm0(1) + pomega(3,l1)*pm(i,3)*pm0(3)               &
       + pomega(5,l1)*pm(i,5)*pm0(5) + pomega(7,l1)*pm(i,7)*pm0(7)

  sum2 =                                                                       &
       pomega(1,l2)*pm(i,1)*pm0(1) + pomega(3,l2)*pm(i,3)*pm0(3)               &
       + pomega(5,l2)*pm(i,5)*pm0(5) + pomega(7,l2)*pm(i,7)*pm0(7)

  sum1 =                                                                       &
       pomega(2,l1)*pm(i,2)*pm0(2) + pomega(4,l1)*pm(i,4)*pm0(4)               &
       + pomega(6,l1)*pm(i,6)*pm0(6) + pomega(8,l1)*pm(i,8)*pm0(8)

  sum3 =                                                                       &
       pomega(2,l2)*pm(i,2)*pm0(2) + pomega(4,l2)*pm(i,4)*pm0(4)               &
       + pomega(6,l2)*pm(i,6)*pm0(6) + pomega(8,l2)*pm(i,8)*pm0(8)

  h(i,l1) = 0.5e0*(sum0*fz(l1) + sum2*fz(l2))
  a(i,l1) = 0.5e0*(sum1*fz(l1) + sum3*fz(l2))

END DO

DO i = 1,m_
  DO j = 1,i
    sum0 =                                                                     &
         pomega(1,l1)*pm(i,1)*pm(j,1) + pomega(3,l1)*pm(i,3)*pm(j,3)           &
         + pomega(5,l1)*pm(i,5)*pm(j,5) + pomega(7,l1)*pm(i,7)*pm(j,7)

    sum2 =                                                                     &
         pomega(1,l2)*pm(i,1)*pm(j,1) + pomega(3,l2)*pm(i,3)*pm(j,3)           &
         + pomega(5,l2)*pm(i,5)*pm(j,5) + pomega(7,l2)*pm(i,7)*pm(j,7)

    sum1 =                                                                     &
         pomega(2,l1)*pm(i,2)*pm(j,2) + pomega(4,l1)*pm(i,4)*pm(j,4)           &
         + pomega(6,l1)*pm(i,6)*pm(j,6) + pomega(8,l1)*pm(i,8)*pm(j,8)

    sum3 =                                                                     &
         pomega(2,l2)*pm(i,2)*pm(j,2) + pomega(4,l2)*pm(i,4)*pm(j,4)           &
         + pomega(6,l2)*pm(i,6)*pm(j,6) + pomega(8,l2)*pm(i,8)*pm(j,8)

    s(i,j) = - sum2*wt(j)
    s(j,i) = - sum2*wt(i)
    t(i,j) = - sum1*wt(j)
    t(j,i) = - sum1*wt(i)
    v(i,j) = - sum3*wt(j)
    v(j,i) = - sum3*wt(i)
    b(i,j,l1) = - 0.5e0*(sum0 + sum2)*wt(j)
    b(j,i,l1) = - 0.5e0*(sum0 + sum2)*wt(i)

  END DO
END DO

DO i = 1,m_
  s(i,i)   = s(i,i)    + 1.0e0
  t(i,i)   = t(i,i)    + 1.0e0
  v(i,i)   = v(i,i)    + 1.0e0
  b(i,i,l1)= b(i,i,l1) + 1.0e0

  c(i,l1) = s(i,1)*a(1,l1)/emu(1) + s(i,2)*a(2,l1)/emu(2)                      &
       + s(i,3)*a(3,l1)/emu(3) + s(i,4)*a(4,l1)/emu(4)
END DO

DO i = 1,m_
  DO j = 1,m_
    w(j,i) = s(j,1)*t(1,i)/emu(1) + s(j,2)*t(2,i)/emu(2)                       &
         + s(j,3)*t(3,i)/emu(3) + s(j,4)*t(4,i)/emu(4)
    u(j,i) = s(j,1)*v(1,i)/emu(1) + s(j,2)*v(2,i)/emu(2)                       &
         + s(j,3)*v(3,i)/emu(3) + s(j,4)*v(4,i)/emu(4)
  END DO
END DO

!------------lower boundary, 2nd-order, a-matrix is full (aa)
deltau = ztau(l1) - ztau(l2)
d2 = 0.25e0*deltau
surfac = 4.0e0*rfl/(1.0e0 + rfl)

DO i = 1,m_

  d1 = emu(i)/deltau
  sum0 = d1 + d2*(w(i,1)+w(i,2)+w(i,3)+w(i,4))
  sum1 = surfac*sum0

  DO j = 1,m_

    aa(i,j,l1) = - d2*u(i,j)
    b(i,j,l1) = b(i,j,l1) + d2*w(i,j) - sum1*emu(j)*wt(j)

  END DO

  h(i,l1) = h(i,l1) - 2.0e0*d2*c(i,l1) + sum0*zflux

END DO

DO i = 1,m_

  d1 = emu(i)/deltau
  aa(i,i,l1) = aa(i,i,l1) + d1
  b(i,i,l1)  = b(i,i,l1) + d1
  c(i,l1) = 0.0e0

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE gen_id

!###############################################################################
!
!  invert 4x4 matrix b(4,4) with l-u decomposition, return as a(4,4) (mjp, old)

SUBROUTINE matinw (b,a)

IMPLICIT NONE

REAL, INTENT(IN)  :: b(4,4,w_)
REAL, INTENT(OUT) :: a(4,4,w_)

INTEGER :: k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MATINW'

! ****************************
! End of Header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

a(:,:,:) = b(:,:,:)

DO k=1,w_

  a(1,1,k) = 1.0e0/a(1,1,k)
  !---setup l and u
  a(2,1,k) = a(2,1,k)*a(1,1,k)
  a(2,2,k) = a(2,2,k)-a(2,1,k)*a(1,2,k)
  a(2,2,k) = 1.0e0/a(2,2,k)
  a(2,3,k) = a(2,3,k)-a(2,1,k)*a(1,3,k)
  a(2,4,k) = a(2,4,k)-a(2,1,k)*a(1,4,k)
  a(3,1,k) = a(3,1,k)*a(1,1,k)
  a(3,2,k) = (a(3,2,k)-a(3,1,k)*a(1,2,k))*a(2,2,k)
  a(3,3,k) = a(3,3,k)-a(3,1,k)*a(1,3,k)-a(3,2,k)*a(2,3,k)
  a(3,3,k) = 1.0e0/a(3,3,k)
  a(3,4,k) = a(3,4,k)-a(3,1,k)*a(1,4,k)-a(3,2,k)*a(2,4,k)
  a(4,1,k) = a(4,1,k)*a(1,1,k)
  a(4,2,k) = (a(4,2,k)-a(4,1,k)*a(1,2,k))*a(2,2,k)
  a(4,3,k) = (a(4,3,k)-a(4,1,k)*a(1,3,k)-a(4,2,k)*a(2,3,k))*a(3,3,k)
  a(4,4,k) =                                                                   &
     a(4,4,k)-a(4,1,k)*a(1,4,k)-a(4,2,k)*a(2,4,k)-a(4,3,k)*a(3,4,k)

  !---invert l
  a(4,3,k) = -a(4,3,k)
  a(4,2,k) = -a(4,2,k)-a(4,3,k)*a(3,2,k)
  a(4,1,k) = -a(4,1,k)-a(4,2,k)*a(2,1,k)-a(4,3,k)*a(3,1,k)
  a(3,2,k) = -a(3,2,k)
  a(3,1,k) = -a(3,1,k)-a(3,2,k)*a(2,1,k)
  a(2,1,k) = -a(2,1,k)

  !---invert u
  a(4,4,k) = 1.0e0/a(4,4,k)
  a(3,4,k) = -a(3,4,k)*a(4,4,k)*a(3,3,k)
  a(2,4,k) = -(a(2,3,k)*a(3,4,k)+a(2,4,k)*a(4,4,k))*a(2,2,k)
  a(2,3,k) = -a(2,3,k)*a(3,3,k)*a(2,2,k)
  a(1,4,k) =                                                                   &
     -(a(1,2,k)*a(2,4,k)+a(1,3,k)*a(3,4,k)+a(1,4,k)*a(4,4,k))*a(1,1,k)
  a(1,3,k) = -(a(1,2,k)*a(2,3,k)+a(1,3,k)*a(3,3,k))*a(1,1,k)
  a(1,2,k) = -a(1,2,k)*a(2,2,k)*a(1,1,k)

  !---multiply (k,u-inverse)*(k,l-inverse)
  a(1,1,k) =                                                                   &
     a(1,1,k)+a(1,2,k)*a(2,1,k)+a(1,3,k)*a(3,1,k)+a(1,4,k)*a(4,1,k)
  a(1,2,k) = a(1,2,k)+a(1,3,k)*a(3,2,k)+a(1,4,k)*a(4,2,k)
  a(1,3,k) = a(1,3,k)+a(1,4,k)*a(4,3,k)
  a(2,1,k) = a(2,2,k)*a(2,1,k)+a(2,3,k)*a(3,1,k)+a(2,4,k)*a(4,1,k)
  a(2,2,k) = a(2,2,k)+a(2,3,k)*a(3,2,k)+a(2,4,k)*a(4,2,k)
  a(2,3,k) = a(2,3,k)+a(2,4,k)*a(4,3,k)
  a(3,1,k) = a(3,3,k)*a(3,1,k)+a(3,4,k)*a(4,1,k)
  a(3,2,k) = a(3,3,k)*a(3,2,k)+a(3,4,k)*a(4,2,k)
  a(3,3,k) = a(3,3,k)+a(3,4,k)*a(4,3,k)
  a(4,1,k) = a(4,4,k)*a(4,1,k)
  a(4,2,k) = a(4,4,k)*a(4,2,k)
  a(4,3,k) = a(4,4,k)*a(4,3,k)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE matinw

!#############################################################################

END SUBROUTINE fastjx_miesct

END MODULE fastjx_miesct_mod
