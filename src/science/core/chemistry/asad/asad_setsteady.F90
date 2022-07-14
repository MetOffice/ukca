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
!    Sets up ASAD variables needed for the calculation of steady-state
!    species.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! *** *set_steady* - chemistry driver routine
!
!     ASAD: set_steady             Version:  set_steady.f 4.1 01/15/97
!
!     Purpose
!     -------
!
!     Sets up ASAD variables needed for the calculation of SS species.
!
!     Interface
!     ---------
!     Called from cinit.
!
!     No arguments
!
!     Method
!     ------
!
!     The routine goes through the arrays defining chemistry and finds
!     SS species and associated source and sink reactions.
!
!     Externals
!     ---------
!     ereport -- report errors in setup.
!
!     Local variables
!     ---------------
!     indss  -- index of steady-state species
!     indsx  -- index of reactions of steady-state species
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE asad_setsteady_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_SETSTEADY_MOD'

CONTAINS

SUBROUTINE asad_setsteady

USE asad_mod, ONLY: ctype, h_in_ss,                                            &
                    jpmsp, jpnr, jpspec, jpss, jpssr,                          &
                    n_in_ss, nspi,                                             &
                    nssi, nsspi, nsspt, nssri, nssrt, nssrx, nsst,             &
                    o1d_in_ss, o3p_in_ss, speci
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper

USE ukca_um_legacy_mod, ONLY: mype
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

! Local variables
INTEGER :: js
INTEGER :: jr
INTEGER :: jp
INTEGER :: ix
INTEGER :: i
INTEGER :: j
INTEGER :: k
INTEGER :: issr
INTEGER :: indss(jpss)
INTEGER :: issx(jpssr)
CHARACTER(LEN=errormessagelength) :: cmessage        ! error message
INTEGER :: errcode                   ! variable passed to ereport

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_SETSTEADY'


! Zero all arrays
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
nsst=0
nssi=0
indss=0
nssrt=0
nssrx=0
nsspt=0
nssri=0
nsspi=0
ix=0

! Select steady state species - O(1D) before O(3P)
DO js = 1,jpspec
  IF (speci(js) == 'O(1D)     ' .AND. o1d_in_ss) THEN
    indss(1)=js
    ix = ix + 1
  ELSE IF (speci(js) == 'O(3P)     ' .AND. o3p_in_ss) THEN
    indss(2)=js
    ix = ix + 1
  ELSE IF (speci(js) == 'N         ' .AND. n_in_ss) THEN
    indss(3)=js
    ix = ix + 1
  ELSE IF (speci(js) == 'H         ' .AND. h_in_ss) THEN
    indss(4)=js
    ix = ix + 1
  ELSE IF (ctype(js) == 'SS') THEN
    ix=ix+1
    indss(ix)=js
  END IF
END DO
!
!  Check that we want them in steady state
DO i=1,jpss
  IF (indss(i) /= 0) THEN
    IF (ctype(indss(i)) == 'SS') THEN
      nsst=nsst+1
      nssi(nsst)=indss(i)
    END IF
  END IF
END DO
IF (nsst == 0) GO TO 9999
!
! Catch array bounds problems (IF necessary)
IF (nsst > jpss) THEN
  errcode=1
  cmessage='Too many steady state species - increase jpss in asad_mod'
  WRITE(umMessage,'(A,I0)') 'INCREASE JPSS TO: ',nsst
  CALL umPrint(umMessage,src='asad_setsteady')

  CALL ereport('ASAD_SETSTEADY',errcode,cmessage)
END IF
!
!  Loop through reactions collecting relevant ones
DO jr = 1,jpnr
  DO jp = 1,jpmsp
    DO ix = 1,nsst
      IF (nspi(jr,jp) == nssi(ix)) THEN
        IF (jp <= 2) THEN
          nssrt(ix) = nssrt(ix)+1
          nssri(ix,nssrt(ix)) = jr
          nssrx(ix,nssrt(ix)) = 3-jp
        ELSE
          nsspt(ix) = nsspt(ix)+1
          nsspi(ix,nsspt(ix)) = jr
        END IF
      END IF
    END DO
  END DO
END DO
!
! Check array bounds
DO ix=1,nsst
  i = MAX(nssrt(ix),nsspt(ix))
  IF (i > jpssr) THEN
    errcode=2
    cmessage=' Too many reactions for steady state '//speci(nssi(ix))//        &
        '- increase jpssr in asad_mod'
    WRITE(umMessage,'(A,I0)') 'INCREASE JPSSR TO: ',i
    CALL umPrint(umMessage,src='asad_setsteady')

    CALL ereport('ASAD_SETSTEADY',errcode,cmessage)
  END IF
END DO
!
!  Check for self-reactions - can't deal with quadratics at the moment
DO ix=1,nsst
  issr=0
  DO i=1,jpssr
    issx(i)=0
  END DO
  DO i=1,nssrt(ix)
    IF (nspi(nssri(ix,i),nssrx(ix,i)) == nssi(ix)) THEN
      IF (issr == 0) THEN
        cmessage=' Steady state '//speci(nssi(ix))//                           &
                 ' has self reaction - disallowed'
        errcode=-1

        CALL ereport('ASAD_SETSTEADY',errcode,cmessage)
      END IF
      issr=issr+1
      issx(issr)=i
    END IF
  END DO
  !  Complain - THEN remove reaction from list
  IF (issr /= 0) THEN
    nssrt(ix)=nssrt(ix)-issr
    DO i=1,nssrt(ix)
      k=i
      DO j=1,issr
        IF (k >= issx(j)) k=k+1
      END DO
      nssri(ix,i) = nssri(ix,k)
      nssrx(ix,i) = nssrx(ix,k)
    END DO
  END IF
END DO
!
!  Output details
DO ix=1,nsst
  IF (mype == 0 .AND. printstatus >= PrStatus_Oper) THEN
    WRITE(umMessage,'(3A)') ' Putting ',speci(nssi(ix)),'in steady state:'
    CALL umPrint(umMessage,src='asad_setsteady')
  END IF
END DO

9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_setsteady
END MODULE asad_setsteady_mod
