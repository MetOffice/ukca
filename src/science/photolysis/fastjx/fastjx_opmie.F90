! *****************************COPYRIGHT*******************************
! (c) [University of California] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

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
!   Fast-jx routine for calculating online photolysis rates
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################

MODULE fastjx_opmie_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FASTJX_OPMIE_MOD'

CONTAINS

SUBROUTINE fastjx_opmie (dtaux,pomegax,u0,rfl,amf2,jxtra,                      &
  fjact,fjtop,fjbot,fsbot,fjflx,flxd,flxd0)
!
!      This routine is one of the central routines for the fast-jx scattering.
!      It refines various quantities such as the scattering phase functions to
!      feed into the routine that actually performs the scattering calculations.

USE fastjx_data,    ONLY: atau0, lpar, m2_, n_, w_
USE ereport_mod,    ONLY: ereport
USE yomhook,        ONLY: lhook, dr_hook
USE parkind1,       ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength
USE fastjx_miesct_mod, ONLY: fastjx_miesct
IMPLICIT NONE

REAL, INTENT(IN)    :: dtaux((lpar+1),w_)       ! optical depth on model levels
REAL, INTENT(IN)    :: pomegax(8,(lpar+1),w_) ! sctterng phase fns on model lvls
REAL, INTENT(IN)    :: amf2(2*(lpar+1)+1,2*(lpar+1)+1) ! air mass factors
REAL, INTENT(IN)    :: u0                              ! cos(solar zenith angle)
REAL, INTENT(IN)    :: rfl(w_)                         ! surface albedo
INTEGER, INTENT(IN) :: jxtra((2*(lpar+1))+1) ! no extra lvls for opt. dense lyrs

REAL, INTENT(OUT)   :: fjact(lpar,w_)    ! mean intensity on model levels
REAL, INTENT(OUT)   :: fjtop(w_)         ! scattered rad. flux at top of atm
REAL, INTENT(OUT)   :: fjbot(w_)         ! scattered rad. flux at bottom of atm
REAL, INTENT(OUT)   :: fsbot(w_)        ! reflected flux at bottom of atmosphere
REAL, INTENT(OUT)   :: fjflx(lpar,w_)    ! mean diffuse flux on model levels
REAL, INTENT(OUT)   :: flxd((lpar+1),w_) ! drct sol. flx deposited in model lvl
REAL, INTENT(OUT)   :: flxd0(w_)  ! integr. sol. flx deposited over whole column

! End of I/O
INTEGER :: jndlev(lpar)
INTEGER :: jnelev((lpar+1))
INTEGER :: jaddlv((2*(lpar+1))+1)
INTEGER :: jaddto((2*(lpar+1))+1),l2lev((2*(lpar+1))+1)

! Loop variables
INTEGER :: i, ii, k, l, ll
INTEGER :: l2, l2l, l22, lz, lzz, nd



REAL    :: dtau((lpar+1)+1,w_)
REAL    :: ttau((2*(lpar+1))+1,w_)
REAL    :: ftau2((2*(lpar+1))+1,w_)

REAL    :: ataua
REAL    :: atauz
REAL    :: xltau
REAL    :: taudn
REAL    :: tauup
REAL    :: dtauj
REAL    :: fjflx0

REAL    :: taubtm(w_)
REAL    :: tautop(w_)
REAL    :: fbtm(w_)
REAL    :: ftop(w_)
REAL    :: zflux(w_)

! Variables added to try and speed up code
REAL    :: dtau_k((lpar+1)+1)
REAL    :: amf2_ll(2*(lpar+1)+1)
REAL    :: pomega_up(m2_), pomega_down(m2_)
REAL    ::  pomega_k(m2_,n_), pomegaj_k(m2_,(2*(lpar+1))+1)
REAL    ::  pomegax_k(8,(lpar+1))
REAL    :: pomegab_k(m2_)
REAL, PARAMETER :: point_five=0.5e0
REAL, PARAMETER :: nada=0.0e0

!--- variables used in mie code-----------------------------------------
REAL    :: fjt(w_)
REAL    :: fjb(w_)
REAL    :: fj(n_,w_)
REAL    :: fz(n_,w_)
REAL    :: ztau(n_,w_)
REAL    :: pomega(m2_,n_,w_)
REAL    :: flxd2(2*(lpar+1),w_)

INTEGER                   :: errcode            ! error code
CHARACTER (LEN=errormessagelength)        :: cmessage           ! error message
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_OPMIE'


!-----------------------------------------------
! End of Header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO l2 = 1,(2*(lpar+1)),1
  jaddlv(l2) = jxtra(l2)
END DO
jaddto((2*(lpar+1))+1) = 0
DO l2 = (2*(lpar+1)),1,-1
  jaddto(l2) = jaddto(l2+1) + jaddlv(l2)
END DO

!---expanded grid now included ctm edge and mid layers plus expanded
!---    grid to allow for finer delta-tau at tops of clouds.
!---    dim of new grid = (2*(lpar+1)) + jaddto(1) + 1

!---l2lev(l2) = l2-index for old level l2 in expanded j-grid (w/jaddlv)
!     in absence of jaddlv, l2lev(l2) = l2
l2lev(1)  = 1
DO l2 = 2,(2*(lpar+1))+1
  l2lev(l2) = l2lev(l2-1) + 1 + jaddlv(l2-1)
END DO

!---jndlev(l=1:lpar) = l2-index in expanded grid for ctm mid-layer l
!---jnelev(l=1:lpar) = l2-index for top of layer l
DO l = 1,lpar
  jndlev(l) = l2lev(2*l)
  jnelev(l) = l2lev(2*l+1)
END DO
jnelev(lpar+1) = 0  ! need to set this to top-of-atmosphere

nd = 2*(2*(lpar+1)) + 2*jaddto(1) + 1

IF (nd > n_) THEN
  errcode = 1
  cmessage = ' Overflow of scatter arrays'
  CALL ereport('FASTJX_OPMIE',errcode,cmessage)
END IF

! ---------------Begin wavelength dependent set up------------------------------

!---reinitialize arrays
ztau(:,:)     = nada
fz(:,:)       = nada
!pomega(:,:,:) = nada

DO k=1,w_

  pomega_k(:,:) = nada
  pomegax_k(:,:) = pomegax(:,:,k)

  !---set up optical depth dtau(l)
  DO l = 1,(lpar+1)
    dtau(l,k) = dtaux(l,k)
  END DO
  dtau((lpar+1)+1,k) = 0.0e0

  ! Try and speed code up by slicing in k
  dtau_k(:) = point_five*dtau(:,k)

  !---define the total scattering phase fn for each ctm layer l=1:lpar+1
  !---   from a dtau-wt_d mix of aerosols, cloud & rayleigh
  !---no. of quadrature pts fixed at 4(m_), expansion of phase fn @ 8
  DO l = 1,(lpar+1)
    DO i = 1,m2_
      pomegaj_k(i,l) = pomegax_k(i,l)
    END DO
  END DO

  !-  calculate attenuated incident beam exp(-ttau/u0 = dtau * amf)
  !-      at the middle & edges of the ctm layers l=1:2*(lpar+1)+1
  !-  (lpar+1) is top-edge of ctm (ie, l=38 = 2 hpa) which has tau > 0
  !-  note that dtau((lpar+1)) is optical depth in the full ctm layer just above
  ftau2(:,:) = 0.0e0
  ftau2((2*(lpar+1))+1,:) = 1.0e0
  DO ll = 1,2*(lpar+1)+1
    l = (ll+1)/2

    IF (amf2(ll,ll) > 0.0e0) THEN
      ! Try and speed up code by slicing in ll
      amf2_ll(:) = amf2(:,ll)

      xltau = 0.0e0
      DO ii = 1,2*(lpar+1)+1
        i = ii+1
        i = i/2
        xltau = xltau + dtau_k(i)*amf2_ll(ii)
      END DO
      IF (xltau < 76.0e0) THEN   ! zero out flux at 1e-33
        ftau2(ll,k) = EXP(-xltau)
      END IF
    END IF
  END DO     ! ll

  !-- calculate direct solar flux deposited in each ctm half-layer:
  !-- l=1:(2*(lpar+1))
  !--     use fsbot for surface flux, cannot DO layer above ctm (lpar+1)
  flxd2(:,:) = 0.0e0
  DO ll = 1,2*(lpar+1)
    IF (amf2(ll,ll) > 0.0e0) THEN
      flxd2(ll,k) = (ftau2(ll+1,k) - ftau2(ll,k))/amf2(ll,ll)
    END IF
  END DO
  IF (amf2(1,1) > 0.0e0) THEN
    fsbot(k) = ftau2(1,k)/amf2(1,1)
  ELSE
    fsbot(k) = 0.0e0
  END IF

  DO ll = 2,2*(lpar+1),2
    l=ll/2
    flxd(l,k) = flxd2(ll,k)+flxd2(ll-1,k)
  END DO

  !--integrate solar flux depositied in ctm layers l=1:lpar, cannot DO top layer
  !--  note flxd0 .ne. (1.d0 - ftau(lpar+1))/amf(lpar+1,lpar+1)
  !--  with spherical atmos
  flxd0(k) = 0.0e0
  IF (amf2(2*(lpar+1),2*(lpar+1)) > 0.0e0) THEN
    DO l=1,(lpar+1)
      flxd0(k) = flxd0(k) + flxd(l,k)
    END DO
  END IF

  !------------------------------------------------------------------------
  !  take optical properties on ctm layers and convert to a photolysis
  !  level grid corresponding to layer centres and boundaries. this is
  !  required so that j-values can be calculated for the centre of ctm
  !  layers; the index of these layers is kept in the jndlev array.
  !------------------------------------------------------------------------
  !---now combine the ctm layer edges (1:lpar+2) with the ctm mid-layer
  !---    points (1:lpar) plus 1 for the mid point of added top layer.
  !---combine these edge- and mid-layer points into grid of size:
  !---              (2*(lpar+1))+1 = 2*(lpar+1)+1 = 2*lpar+3
  !---calculate column optical depths above each level, ttau(1:(2*(lpar+1))+1)
  !---      note that ttau((2*(lpar+1))+1)=0 and ttau(1)=total od

  ttau((2*(lpar+1))+1,k) = 0.0e0
  DO l2 = (2*(lpar+1)),1,-1
    l          = (l2+1)/2
    dtauj      = point_five * dtau(l,k)
    ttau(l2,k)   = ttau(l2+1,k) + dtauj
  END DO

  ! ---solar flux incident on lower boundary & lambertian reflect factor:
  IF (fsbot(k) > 0.0e0) THEN
    zflux(k) = fsbot(k)*rfl(k)/(1.0e0+rfl(k))
  ELSE
    zflux(k) = 0.0e0
  END IF

  ! --Calculate scattering properties, level centres then level boundaries
  ! --be careful of order, we are overwriting/shifting the 'pomegaj' up in index
  DO l2 = (2*(lpar+1)),2,-2
    l   = l2/2
    DO i = 1,m2_
      pomegaj_k(i,l2) = pomegaj_k(i,l)
    END DO
  END DO

  ! --lower boundary value is set (pomegaj(i,1), but set upper:
  DO i = 1,m2_
    pomegaj_k(i,(2*(lpar+1))+1) = pomegaj_k(i,(2*(lpar+1)))
  END DO

  ! --now have pomegaj filled at even points from l2=3:(2*(lpar+1))-1
  ! --use inverse interpolation for correct tau-weighted values at edges
  DO l2 = 3,(2*(lpar+1))-1,2
    taudn = ttau(l2-1,k)-ttau(l2,k)
    tauup = ttau(l2,k)-ttau(l2+1,k)
    DO i = 1,m2_
      pomegaj_k(i,l2) = (pomegaj_k(i,l2-1)*taudn +                             &
          pomegaj_k(i,l2+1)*tauup) / (taudn+tauup)
    END DO
  END DO

  ! --at this point ftau2(1:(2*(lpar+1))+1) and
  !    pomeagj(1:8, 1:(2*(lpar+1))+1) where ftau2((2*(lpar+1))+1) = 1.0
  !    = top-of-atmos, ftau2(1) = surface

          ! l2 = index of ctm edge- and mid-layers
  DO l2 = 1,(2*(lpar+1))+1
    ! l2l = index for l2 in expanded scale(jadd)
    l2l = l2lev(l2)
    ! lz = index for l2 in scatt arrays
    lz  = nd + 2 - 2*l2l
    ztau(lz,k) = ttau(l2,k)
    fz(lz,k)   = ftau2(l2,k)
    DO i=1,m2_
      pomega_k(i,lz) = pomegaj_k(i,l2)
    END DO
  END DO

  ! now go thru the pairs of l2 levels to see if we need jadd levels
  ! l2 = index of ctm edge- and mid-layers
  DO l2 = 1,(2*(lpar+1))
    ! l2l = index for l2 in expanded scale(jadd)
    l2l = l2lev(l2)
    ! lz = index for l2 in scatt arrays
    lz  = nd + 2 - 2*l2l
    ! l22 = 0 if no added levels
    l22 = l2lev(l2+1) - l2lev(l2) - 1

    IF (l22 > 0) THEN
      taubtm(k) = ttau(l2,k)
      tautop(k) = ttau(l2+1,k)
      fbtm(k)   = ftau2(l2,k)
      ftop(k)   = ftau2(l2+1,k)
      DO i = 1,m2_
        pomegab_k(i) = pomegaj_k(i,l2)
      END DO

      !---to fit l22 new layers between taubot > tautop,
      ! calculate new 1/atau factor such that tau(just above tau-btm)
      ! equals atuaz * taubtm < taubtm
      atauz = EXP(-LOG(taubtm(k)/MAX(tautop(k),atau0))                         &
            /REAL(l22+1))

      ! add odd levels between l2lev(l2) &l2lev(l2+1)
      DO l = 1,l22

        ! lzz = index(odd) of added level in scatt arrays
        lzz = lz - 2*l
        ztau(lzz,k) = taubtm(k) * atauz

        !---fraction from taubtm=>tautop
        ataua=(taubtm(k)-ztau(lzz,k))/(taubtm(k)-tautop(k))
        !---solar flux at interp-levels: use EXP(tau/u0) if
        ! u0>0.02 (89 deg), else scale by tau
        IF (u0 > 0.02e0) THEN
          fz(lzz,k) = ftop(k) * EXP((tautop(k)-ztau(lzz,k))/u0)
        ELSE
          IF (fbtm(k) < 1.0e-32) THEN
            fz(lzz,k) = 0.0e0
          ELSE
            fz(lzz,k) = fbtm(k) * (ftop(k)/fbtm(k))**ataua
          END IF
        END IF
        DO i = 1,m2_
          pomega_k(i,lzz) = pomegab_k(i) +                                     &
                ataua*(pomegaj_k(i,l2+1)-pomegab_k(i))
        END DO
        taubtm(k)    = ztau(lzz,k)
        fbtm(k)      = fz(lzz,k)
        DO i = 1,m2_
          pomegab_k(i) = pomega_k(i,lzz)
        END DO
      END DO
    END IF
  END DO

  ! Fill in the even points with simple interp. in scatter arrays:
  DO lz = 2,nd-1,2
    ztau(lz,k) = point_five*(ztau(lz-1,k)+ztau(lz+1,k))
    fz(lz,k)   = SQRT(fz(lz-1,k)*fz(lz+1,k))

    pomega_up(:) = pomega_k(:,lz+1)
    pomega_down(:) = pomega_k(:,lz-1)

    DO i=1,m2_
      pomega_k(i,lz) = point_five*(pomega_down(i)+pomega_up(i))
    END DO
  END DO

  pomega(:,:,k) = pomega_k(:,:)

END DO  ! wavelength loop!

!----------------------------------------------------------------------
CALL fastjx_miesct(fj,fjt,fjb,pomega,fz,ztau,zflux,rfl,u0,nd)


!---Move mean intensity from scatter array fj(lz=1:nd)
!--- to ctm mid-level array fjact(l=1:lpar)

DO k=1,w_

  !---mean intensity:  4*<i> + solar at mid-layer
  DO l = 1,lpar
    l2l = jndlev(l)
    lz  = nd+2 - 2*l2l
    fjact(l,k) = 4.0e0*fj(lz,k) + fz(lz,k)
  END DO

  !---mean diffuse flux:  4<i*mu> (not solar) at top of layer l
  !---average (tau-wtd) the h's just above and below the l-edge
  DO l = 1,lpar
    l2l = jnelev(l)
    lz  = nd+2 - 2*l2l
    fjflx0 = (ztau(lz+1,k)-ztau(lz,k))/(ztau(lz+1,k)-ztau(lz-1,k))
    fjflx(l,k)=4.0e0*(fj(lz-1,k)*fjflx0 + fj(lz+1,k)*(1.0e0-fjflx0))

  END DO

  !---dIFfuse fluxes reflected at top, incident at bottom
  fjtop(k) = fjt(k)
  fjbot(k) = fjb(k)

END DO  ! wavelength loop!

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE fastjx_opmie
END MODULE fastjx_opmie_mod
