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
!   Set aerosol/cloud types
!
!     MX       Number of different types of aerosol to be considered
!     MIEDX    Index of aerosol types in jv_spec.dat - hardwire in here
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
MODULE fastjx_set_aer_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FASTJX_SET_AER_MOD'

CONTAINS

SUBROUTINE fastjx_set_aer

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper
USE fastjx_data, ONLY: mx, miedx, naa
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

INTEGER                       :: i                ! Loop variable
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_SET_AER'
CHARACTER(LEN=errormessagelength) :: cmessage   ! Error message


! ********************************
! EOH
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise aerosol index

DO i = 1,mx
  miedx(i) = 0
END DO

! Select Aerosol/Cloud types to be used - define types here
miedx(1) = 9     !  Water Cloud (Deirmenjian 8 micron)
miedx(2) = 13    !  Irregular Ice Cloud (Mishchenko)
miedx(3) = 16    !  UT sulphate (CHECK: doesn't exactly correspond to fastjx)

! Loop over mx types
DO i = 1,mx
  IF (printstatus >= prstatus_oper) THEN
    WRITE(umMessage,'(A,I0,A,I0)') 'Mie scattering type ', i, ' ', miedx(i)
    CALL umPrint(umMessage,src='fastjx_set_aer')
  END IF

  IF (miedx(i) > naa .OR. miedx(i) <= 0) THEN
    WRITE(umMessage,'(I0,A,I0)') miedx(i), ' ', naa
    CALL umPrint(umMessage,src='fastjx_set_aer')
    cmessage= 'MIEDX(i) is negative or less than naa'
    CALL ereport('FASTJX_SET_AER',naa,cmessage)
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_set_aer
END MODULE fastjx_set_aer_mod
