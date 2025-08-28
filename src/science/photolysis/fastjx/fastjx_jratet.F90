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
!  Description:
!   Fast-JX routine for calculating online photolysis rates
!
!  Calculate and print J-values. Note that the loop in this routine
!  only covers the jpcl levels actually needed by the CTM.
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
MODULE fastjx_jratet_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FASTJX_JRATET_MOD'

CONTAINS

SUBROUTINE fastjx_jratet(ppj,ttj,fff, valjl)
!
!     Description: Routine which takes pressure profiles at edges of the grid
!     boxes (ppj), temperatures at mid-level (ttj),and the mean actinic flux
!     (fff) as inputs to calculate photolysis rates as a function of model
!     levels and photolysed species
!
! #######################################################################
USE fastjx_data, ONLY: lpar, jpcl, njval, q1d, qqq,                            &
                       qo2, qo3, tqq, w_, x_, titlej
USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim
IMPLICIT NONE

! Interface
REAL,   INTENT(IN)  ::  ppj((lpar+1)+1)    ! pressure profile at edges
REAL,   INTENT(IN)  ::  ttj((lpar+1))      ! temperatures at mid-level
REAL,   INTENT(IN)  ::  fff(w_,jpcl)       ! mean actinic flux
REAL,   INTENT(OUT) ::  valjl(jpcl,njval)  ! photolysis rates

! Local variables
REAL :: valj(x_)          ! temp for call j's at one l
REAL :: qo2tot(w_)
REAL :: qo3tot(w_)
REAL :: qo31dy(w_)
REAL :: qo31d
REAL :: qo3p
REAL :: qqqt
REAL :: tfact
REAL :: tt                ! temperature at mid layer (K)
REAL :: pp                ! pressure at mid layer
REAL :: dd                ! density at mid layer
REAL :: tt200
REAL :: tfaca
REAL :: tfac0
REAL :: tfac1
REAL :: tfac2
REAL :: qqqa
REAL :: qq2
REAL :: qq1a
REAL :: qq1b

! Loop variables
INTEGER ::  j,k,l, iv

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_JRATET'

! *************************************************************
! End of Header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO l = 1,jpcl    ! master loop over layer = l

  !Need temperature, pressure, and dens. at mid-layer (for some quantum yields):
  tt   = ttj(l)
  IF (l  ==  1) THEN
    pp = ppj(1)
  ELSE
    pp  = (ppj(l)+ppj(l+1))*0.5e0
  END IF
  dd = 7.24e18*pp/tt

  ! Initialise to 0
  DO j = 1,njval
    valj(j) = 0.0e0
  END DO

  ! Interpolate temperatures
  IF (tt <= tqq(2,1)) THEN
    IF (tt <= tqq(1,1)) THEN
      tfact = 0.0e0
    ELSE
      tfact = (tt -tqq(1,1))/(tqq(2,1) -tqq(1,1))
    END IF
    DO k = 1,w_
      qo2tot(k) = qo2(k,1) + (qo2(k,2) - qo2(k,1))*tfact
    END DO
  ELSE
    IF (tt >= tqq(3,1)) THEN
      tfact = 1.0e0
    ELSE
      tfact = (tt -tqq(2,1))/(tqq(3,1) -tqq(2,1))
    END IF
    DO k = 1,w_
      qo2tot(k) = qo2(k,2) + (qo2(k,3) - qo2(k,2))*tfact
    END DO
  END IF

  IF (tt <= tqq(2,2)) THEN
    IF (tt <= tqq(1,2)) THEN
      tfact = 0.0e0
    ELSE
      tfact = (tt -tqq(1,2))/(tqq(2,2) -tqq(1,2))
    END IF
    DO k = 1,w_
      qo3tot(k) = qo3(k,1) + (qo3(k,2) - qo3(k,1))*tfact
    END DO
  ELSE
    IF (tt >= tqq(3,2)) THEN
      tfact = 1.0e0
    ELSE
      tfact = (tt -tqq(2,2))/(tqq(3,2) -tqq(2,2))
    END IF
    DO k = 1,w_
      qo3tot(k) = qo3(k,2) + (qo3(k,3) - qo3(k,2))*tfact
    END DO
  END IF

  IF (tt <= tqq(2,3)) THEN
    IF (tt <= tqq(1,3)) THEN
      tfact = 0.0e0
    ELSE
      tfact = (tt -tqq(1,3))/(tqq(2,3) -tqq(1,3))
    END IF
    DO k = 1,w_
      qo31dy(k) = q1d(k,1) + (q1d(k,2) - q1d(k,1))*tfact
    END DO
  ELSE
    IF (tt >= tqq(3,3)) THEN
      tfact = 1.0e0
    ELSE
      tfact = (tt -tqq(2,3))/(tqq(3,3) -tqq(2,3))
    END IF
    DO k = 1,w_
      qo31dy(k) = q1d(k,2) + (q1d(k,3) - q1d(k,2))*tfact
    END DO
  END IF

  ! Loop over wavelengths
  DO k = 1,w_
    qo31d  = qo31dy(k)*qo3tot(k)
    qo3p  = (1-qo31dy(k))*qo3tot(k)
    valj(1) = valj(1) + qo2tot(k)*fff(k,l)   ! O2 -> O(3P) + O(3P)
    valj(2) = valj(2) + qo3p*fff(k,l)        ! O3 -> O2 + O(3P)
    valj(3) = valj(3) + qo31d*fff(k,l)       ! O3 -> O2 + O(1D)
  END DO

  DO j = 4,njval

    IF (tqq(2,j) > tqq(1,j)) THEN
      tfact = MAX(0.0,MIN(1.0,(tt-tqq(1,j))/(tqq(2,j)-tqq(1,j))))
    ELSE
      tfact = 0.0e0
    END IF

    DO k = 1,w_
      qqqt    = qqq(k,1,j) + (qqq(k,2,j) - qqq(k,1,j))*tfact
      valj(j) = valj(j) + qqqt*fff(k,l)
    END DO

    ! ********Special case for methyl glyoxal***********************************

    ! The cross section for methyl glyoxal is taken from the IUPAC assessments
    ! based on the work of Chen et al. (Y. Chen, W. Wang and L. Zhu, J. Phys.
    ! Chem. A, 104 11126 (2000).). This modifies the cross section with a
    ! wavelength dependent quantum yield:
    ! 1/j =1/j_0 + k*P(N_2 (torrs))*exp(-(5639/lambda(nm)))
    ! Here k is a constant 19300, P is the pressure of nitrogen in torrs and
    ! lambda is the wavelength in nm. The value in the code (4.12*10-3) is a
    ! weighted average of this expression over all wavelengths used by fast-jx
    ! (weighting by incident solar flux....not perfect but seems to work).
    ! This is an approximation used to avoid the complexity of having a
    ! dependence of the j rate on wavelength online. It has been shown to give
    ! a good description of the rate in the PHOTOCOMP comparison.

    IF (titlej(j) == 'jmkal') THEN
      valj(j) = 1/((1/valj(j))+4.12e-3*pp)
    END IF
  END DO

  !--------------J-ref v8.3 includes Blitz ACETONE q-yields--------------
  !---Acetone is a special case:   (as per Blitz et al GRL, 2004)
  !---     61 = NJVAL-1 = J1(acetone-a) ==> CH3CO + CH3
  !---     62 = NJVAL   = J2(acetone-b) ==> CH3 + CO + CH3
  ! Note we only use the first channel which means the original fast-jx code
  ! has been modified.

  IF (titlej(njval) == 'jaceto') THEN

    valj(njval) = 0.0

    !---IV=NJVAL = Xsect (total abs) for Acetone - pre-calc Temp interp factors
    iv    = njval
    tfaca = (tt-tqq(1,iv))/(tqq(2,iv)-tqq(1,iv))
    tfaca = MAX(0.0, MIN(1.0, tfaca))

    !--NJVAL+1 = Q2 for Acetone=>(2), specifically designed for quadratic interp
    iv   = njval+1
    tfac0 = ( (tt-tqq(1,iv))/(tqq(2,iv)-tqq(1,iv)) )**2
    IF (tt < tqq(1,iv)) THEN
      tfac0 = (tt - 210.0e0)/(tqq(1,iv)-210.0e0)
    END IF
    tfac0 = MAX(0.0, MIN(1.0, tfac0))

    !---IV=NJVAL+2 = Q1A for Acetone => (1), allow full range of T = 200K-300K
    iv    = njval+2
    tt200 = MIN(300.0, MAX(200.0, tt))
    tfac1 = (tt200-tqq(1,iv))/(tqq(2,iv)-tqq(1,iv))

    !---IV=NJVAL+3 = Q1B for Acetone => (1)
    iv    = njval+3
    tfac2 = (tt200-tqq(1,iv))/(tqq(2,iv)-tqq(1,iv))

    !---NJVAL = Xsect (total abs) for Acetone - integrate over wavelengths
    DO k = 1,w_
      iv   = njval
      qqqa = qqq(k,1,iv) + (qqq(k,2,iv)-qqq(k,1,iv))*tfaca

      !NJVAL+1 = Q2 for Acetone=>(2), specifically designed for quadratic interp
      iv   = njval
      qq2  = qqq(k,1,iv) + (qqq(k,2,iv)-qqq(k,1,iv))*tfac0
      IF (tt < tqq(1,iv)) THEN
        qq2 = qqq(k,1,iv)*tfac0
      END IF

      !---NJVAL+2 = Q1A for Acetone => (1), allow full range of T = 200K-300K
      iv   = njval+2
      qq1a = qqq(k,1,iv) + (qqq(k,2,iv)-qqq(k,1,iv))*tfac1

      !---NJVAL+3 = Q1B for Acetone => (1)   ! scaled to [M]=2.5e19
      iv   = njval+3
      qq1b = qqq(k,1,iv) + (qqq(k,2,iv)-qqq(k,1,iv))*tfac2
      qq1b = qq1b*4.0e-20

      !---J(61)
      valj(njval) = valj(njval)+ fff(k,l)*qqqa*(1.0e0-qq2)/(qq1a + qq1b*dd)
    END DO

  END IF

  !----load array of j-values in native order, need to be indexed/scaled
  !    by asad-related code later: zpj(l,jj) = valjl(l,jind(jj))*jfacta(jj)
  DO j=1,njval
    valjl(l,j) = valj(j)
  END DO

END DO    ! master loop over l=1,jpcl

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE fastjx_jratet
END MODULE fastjx_jratet_mod
