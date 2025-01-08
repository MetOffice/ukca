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
!  Routine taken from fast-jx. Only changes tidying up of format.
!  Calculation of spherical geometry; derive tangent heights, slant path
!  lengths and air mass factor for each layer. Not called when
!  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
!  beam (where tangent height is below altitude J-value desired at).
!
!     GMU     MU, cos(solar zenith angle)
!     RZ      Distance from centre of Earth to each point (cm)
!     RQ      Square of radius ratios
!     TANHT   Tangent height for the current SZA
!     XL      Slant path between points
!     AMF     Air mass factor for slab between level and level above
!     U0      cos(solar zenith angle*pi/180)
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
!
MODULE fastjx_sphere_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FASTJX_SPHERE_MOD'

CONTAINS

SUBROUTINE fastjx_sphere(gmu,rz,amf2,l1_)

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN)  :: l1_                   ! Number of levels
REAL,    INTENT(IN)  :: gmu                   ! cosine sza
REAL,    INTENT(IN)  :: rz(l1_+1)             ! Distance from centre of Earth
                                              !  to each point (cm)
REAL,    INTENT(OUT) :: amf2(2*l1_+1,2*l1_+1) ! return air mass fraction

! End of I/O

REAL                 :: xmu1
REAL                 :: xmu2
REAL                 :: xl                     ! Slant path between points
REAL                 :: diff
REAL                 :: shadht                 ! Shadow height for current zsa
REAL                 :: rz2(2*l1_+1)           ! Distance on split levels
REAL                 :: rq2(2*l1_+1)           ! Ratio of distances

! Loop variables
INTEGER  ::   i, j, ii, l2
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_SPHERE'

!---------------------------------------------------
! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate heights for edges of split ctm-layers
l2 = 2*l1_
DO ii = 2,l2,2
  i = ii/2
  rz2(ii-1) = rz(i)
  rz2(ii) = 0.5e0*(rz(i)+rz(i+1))
END DO
rz2(l2+1) = rz(l1_+1)
DO ii = 1,l2
  rq2(ii) = (rz2(ii)/rz2(ii+1))**2
END DO


! Shadow height for sza > 90
IF (gmu < 0.0e0) THEN
  shadht = rz2(1)/SQRT(1.0e0 - gmu**2)
ELSE
  shadht = 0.0e0
END IF

! Up from the surface calculating the slant paths between each level
!---  and the level above, and deriving the appropriate air mass factor
amf2(:,:) = 0.0e0

DO j = 1,2*l1_+1
  ! air mass factors all zero if below the tangent height
  IF (rz2(j) < shadht) CYCLE

  ! ascend from layer j calculating amf2s
  xmu1 = ABS(gmu)
  DO i = j,2*l1_
    xmu2      =  SQRT(1.0e0 - rq2(i)*(1.0e0-xmu1**2))
    xl        =  rz2(i+1)*xmu2 - rz2(i)*xmu1
    amf2(i,j) =  xl / (rz2(i+1)-rz2(i))
    xmu1      =  xmu2
  END DO

  !--fix above top-of-atmos (l=l1_+1), must set dtau(l1_+1)=0
  amf2(2*l1_+1,j) = 1.0e0

  !  twilight case - emergent beam, calc air mass factors below layer
  IF (gmu >= 0.0e0) CYCLE

  !  Descend from layer j
  xmu1       = ABS(gmu)
  DO ii = j-1,1,-1
    diff        = rz2(ii+1)*SQRT(1.0e0-xmu1**2)-rz2(ii)
    IF (ii == 1)  diff = MAX(diff,0.0e0)   ! filter

    !  tangent height below current level - beam passes through twice
    IF (diff < 0.0e0) THEN
      xmu2       =  SQRT(1.0e0 - (1.0e0-xmu1**2)/rq2(ii))
      xl         =  ABS(rz2(ii+1)*xmu1-rz2(ii)*xmu2)
      amf2(ii,j) =  2.0e0*xl/(rz2(ii+1)-rz2(ii))
      xmu1       =  xmu2

      !  lowest level intersected by emergent beam
    ELSE
      xl        = rz2(ii+1)*xmu1*2.0e0
      amf2(ii,j) = xl/(rz2(ii+1)-rz2(ii))
      EXIT
    END IF
  END DO    ! ii
END DO    ! j

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_sphere
END MODULE fastjx_sphere_mod
