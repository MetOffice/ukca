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
!    Part of Fast-jx, a routine for calculating online photolysis rates
!    Calculates where, due to optical depth, need to add extra layers
!
!    new version 6.1, add sub-layers (jxtra) to thick cloud/aerosol layers
!    this version sets up log-spaced sub-layers of increasing thickness atau
!
!     dtaux(l=1:l1x) = optical depth in layer l (generally 600 nm od)
!        this can be just cloud or cloud+aerosol, it is used only to set
!        the number in levels to insert in each layer l
!        set for log-spacing of tau levels, increasing top-down.
!
!     n.b. the ttau, etc calculated here are not used elsewhere
!
!---the log-spacing parameters have been tested for convergence and chosen
!---  to be within 0.5% for ranges od=1-500, rflect=0-100%, mu0=0.1-1.0
!---  use of atau = 1.18 and min = 0.01, gives at most +135 pts for od=100
!---  atau = 1.12 now recommended for more -accurate heating rates (not j's)
!-----------------------------------------------------------------------
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
MODULE fastjx_extral_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FASTJX_EXTRAL_MOD'

CONTAINS

SUBROUTINE fastjx_extral(dtaux,l1x,l2x,nx,jtaumx,atau,atau0, jxtra)

USE ereport_mod,          ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! index of cloud/aerosol
INTEGER, INTENT(IN)  ::  jtaumx,l1x,l2x
! mie scattering array size
INTEGER, INTENT(IN)  ::  nx
! cloud+3aerosol od in each layer
REAL,    INTENT(IN)  ::  dtaux(l1x)
!
REAL,    INTENT(IN)  ::  atau,atau0
! number of sub-layers to be added
INTEGER, INTENT(OUT) ::  jxtra(l2x+1)

! End of I/O

! Loop variables
INTEGER              :: jtotl,l,l2

REAL                 :: ttau(l2x+1)
REAL                 :: dtauj
REAL                 :: atauln
REAL                 :: ataum
REAL                 :: ataun1
INTEGER              :: errcode
CHARACTER (LEN=errormessagelength)   :: cmessage      ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_EXTRAL'


! **************************************************************
! End of Header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set arrays to 0
ttau(:)  = 0.0e0
jxtra(:) = 0

! Combine these edge- and mid-layer points into grid of size:
!  l2x+1 = 2*l1x+1 = 2*l_+3
! Calculate column optical depths above each level, ttau(1:l2x+1)
!  note that ttau(l2x+1)=0 and ttau(1)=total od

! Divide thick layers to achieve better accuracy in the scattering code
!  in the original fast-j, equal sub-layers were chosen, this is wasteful
!  and this new code (ver 5.3) uses log-scale:
!  each succesive layer (down) increase thickness by atau > 1
!  e.g., if atau = 2, a layer with od = 15 could be divided into
!  4 sub-layers with ods = 1 - 2 - 4 - 8
!  the key parameters are:
!    atau = factor increase from one layer to the next
!    ataumn = the smallest od layer desired
!    jtaumx = maximum number of divisions (i.e., may not get to ataumn)
!  these are hardwired below, can be changed, but have been tested/optimized

atauln = LOG(atau)
ttau(l2x+1)  = 0.0e0
DO l2 = l2x,1,-1
  l         = (l2+1)/2
  dtauj     = 0.5e0 * dtaux(l)
  ttau(l2)  = ttau(l2+1) + dtauj

  ! now compute the number of log-spaced sub-layers to be added in
  ! the interval ttau(l2) > ttau(l2+1)
  ! the objective is to have successive tau-layers increasing by
  ! factor atau >1 the number of sub-layers + 1
  IF (ttau(l2) < atau0) THEN
    jxtra(l2) = 0
  ELSE
    ataum     = MAX(atau0, ttau(l2+1))
    ataun1    = LOG(ttau(l2)/ataum) / atauln
    jxtra(l2) = MIN(jtaumx, MAX(0, INT(ataun1 - 0.5e0)))
  END IF
END DO

! Check on overflow of arrays, cut off jxtra at lower l if too many levels
jtotl    = l2x + 2
l2_loop: DO l2 = l2x,1,-1
  jtotl  = jtotl + jxtra(l2)
  IF (jtotl > nx/2) THEN
    cmessage = ' Too many levels, cutting off'
    errcode  = -1*jtotl
    IF (printstatus > Prstatus_oper) THEN
      WRITE(umMessage,'(A24,3I6)') 'n_/l2_/l2-cutoff jxtra: ',nx,l2x,l2
      CALL umPrint(umMessage,src='fastjx_extral')
    END IF
    CALL ereport('FASTJX_EXTRAL',errcode,cmessage)
    DO l = l2,1,-1
      jxtra(l) = 0
    END DO
    EXIT l2_loop
  END IF
END DO l2_loop

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_extral
END MODULE fastjx_extral_mod
