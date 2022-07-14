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
!    number conc. (ND) & total mass/ptcl (MDT) for each aerosol mode
!    Used for error checking in box model.
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
MODULE ukca_calcminmaxndmdt_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CALCMINMAXNDMDT_MOD'

CONTAINS

SUBROUTINE ukca_calcminmaxndmdt(nbox, nd, mdt, verbose)
!----------------------------------------------------------------------
!
!     Purpose
!     -------
!     Used to calculate and print out max,min,mean of
!     number conc. (ND) & total mass/ptcl (MDT) for each aerosol mode
!     Used for error checking in box model
!
!     Inputs
!     ------
!     NBOX       : Number of grid boxes
!     ND         : Aerosol ptcl number density for mode (cm^-3)
!     MDT        : Avg tot mass of aerosol ptcl in mode (particle^-1)
!     verbose    : Switch for printing out test print statements
!
!     Outputs
!     -------
!     None
!
!     Local variables
!     ---------------
!     XMINN,XMAXN,XMEANN : min,max,mean of particle number conc.
!     XMINM,XMAXM,XMEANM : min,max,mean of avg total mass per particle
!
!--------------------------------------------------------------------
USE ukca_mode_setup,   ONLY: mode, nmodes
USE parkind1,          ONLY: jprb, jpim
USE yomhook,           ONLY: lhook, dr_hook
USE umPrintMgr, ONLY:                                                          &
    umPrint,                                                                   &
    umMessage
IMPLICIT NONE


INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: verbose               ! Verbose level
REAL, INTENT(IN)    :: nd(nbox,nmodes)
REAL, INTENT(IN)    :: mdt(nbox,nmodes)

! Local variables
INTEGER :: imode
INTEGER :: jl
REAL    :: xminn
REAL    :: xmaxn
REAL    :: xmeann
REAL    :: xminm
REAL    :: xmaxm
REAL    :: xmeanm

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALCMINMAXNDMDT'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (verbose > 0) THEN
  !  Write header
  WRITE(umMessage,'(A25,A65)') 'UKCA_CALCMINMAXNDMDT:',                        &
         'imode min(nd) max(nd) mean(nd) min(mdt) max(mdt) mean(mdt)'
  CALL umPrint(umMessage,src='ukca_calcminmaxndmdt')
END IF

DO imode=1,nmodes
  IF (mode(imode)) THEN
    xminn=1.0e16
    xmaxn=-1.0e16
    xmeann=0.0e0
    xminm=1.0e16
    xmaxm=-1.0e16
    xmeanm=0.0
    DO jl=1,nbox
      xminn=MIN(nd(jl,imode),xminn)
      xmaxn=MAX(nd(jl,imode),xmaxn)
      xmeann=xmeann+nd(jl,imode)/REAL(nbox)
      xminm=MIN(mdt(jl,imode),xminm)
      xmaxm=MAX(mdt(jl,imode),xmaxm)
      xmeanm=xmeanm+mdt(jl,imode)/REAL(nbox)
    END DO

    IF (verbose > 0) THEN
      WRITE(umMessage,'(1i4,6e15.6)') imode,xminn,xmaxn,xmeann,xminm,          &
                                      xmaxm,xmeanm
      CALL umPrint(umMessage,src='ukca_calcminmaxndmdt')
    END IF
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calcminmaxndmdt
END MODULE ukca_calcminmaxndmdt_mod
