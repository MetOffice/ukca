! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Used to calculate and print out max,min,mean of
!    each condensable gas phase concentration for
!    error checking in box model.
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
! ######################################################################
!
! Subroutine Interface:
MODULE ukca_calcminmaxgc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CALCMINMAXGC_MOD'

CONTAINS

SUBROUTINE ukca_calcminmaxgc(str_at,nbox,nchemg,gc)
!----------------------------------------------------------------------
!
!     Purpose
!     -------
!     Used to calculate and print out max,min,mean of
!     each condensable gas phase concentration for
!     error checking in box model
!
!     Inputs
!     ------
!     STR_AT     : String indicating which process the code is up to
!     NBOX       : Number of grid boxes
!     nchemg     : Number of gas phase chemistry tracers
!     GC         : Condensable cpt number density (molecules cm-3)
!
!     Outputs
!     -------
!     None
!
!     Local variables
!     ---------------
!     GCMIN,GCMAX,GCMEAN : min,max,mean of gas phase condensable conc.
!
!--------------------------------------------------------------------
USE ukca_setup_indices, ONLY: condensable
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Diag
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Arguments
CHARACTER(LEN=30), INTENT(IN) :: str_at
INTEGER,           INTENT(IN) :: nbox
INTEGER,           INTENT(IN) :: nchemg
REAL,              INTENT(IN) :: gc(nbox,nchemg)


! Local variables
INTEGER :: jl
INTEGER :: jv
REAL    :: gcmin
REAL    :: gcmax
REAL    :: gcmean

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALCMINMAXGC'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO jv=1,nchemg
  IF (condensable(jv)) THEN
    IF (PrintStatus >= PrStatus_Diag) THEN
      CALL umPrint( '**********************************',                      &
          src='ukca_calcminmaxgc')
      CALL umPrint( str_at,src='ukca_calcminmaxgc')
    END IF
    gcmin=1.0e9
    gcmax=-1.0e9
    gcmean=0.0e0
    DO jl=1,nbox
      gcmin=MIN(gc(jl,jv),gcmin)
      gcmax=MAX(gc(jl,jv),gcmax)
      gcmean=gcmean+gc(jl,jv)/REAL(nbox)
    END DO
    IF (PrintStatus >= PrStatus_Diag) THEN
      WRITE(umMessage,'(A20,I6,3E12.3)') 'GC:JV,Min,max,mean=',jv,gcmin,       &
                                 gcmax,gcmean
      CALL umPrint(umMessage,src='ukca_calcminmaxgc')
      CALL umPrint( '**********************************',                      &
          src='ukca_calcminmaxgc')
    END IF
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calcminmaxgc
END MODULE ukca_calcminmaxgc_mod
