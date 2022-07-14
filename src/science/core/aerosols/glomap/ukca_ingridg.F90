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
!    Generate sectional particle mass grid in molecules per particle.
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
MODULE ukca_ingridg_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_INGRIDG_MOD'

CONTAINS

SUBROUTINE ukca_ingridg(nbins,mbsmall,mblarge,mbmid,mblo,mbhi)
!-----------------------------------------------------------------------
!
!     Purpose
!     -------
!     Generate sectional particle mass grid (in molecules per particle)
!
!     Inputs
!     ------
!     NBINS   : Number of aerosol ptcl size bins
!     MBSMALL : Number of molecules for ptcl at smallest bin mid-pt
!     MBLARGE : Number of molecules for ptcl at largest  bin mid-pt
!
!     Outputs
!     -------
!     MBLO, MBMID, MBHI : Mass of lower, mid and upper bin edge
!
!     Local variables
!     ---------------
!     LGMRANGE : Diff. in logarithms of masses of largest/smallest bins
!     LGMGRID  : Diff. in logarithms of masses of limits of each bin
!-----------------------------------------------------------------------

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN)  :: nbins
REAL, INTENT(IN)     :: mbsmall
REAL, INTENT(IN)     :: mblarge
REAL, INTENT(OUT)    :: mbmid(nbins)
REAL, INTENT(OUT)    :: mblo(nbins)
REAL, INTENT(OUT)    :: mbhi(nbins)

! Local variables
INTEGER :: jv
REAL    :: lgmrange
REAL    :: lgmgrid

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INGRIDG'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

lgmrange=LOG(mblarge)-LOG(mbsmall)
DO jv=1,nbins+1
  lgmgrid=LOG(mbsmall)+lgmrange*REAL(jv-1)/REAL(nbins)
  IF (jv < (nbins+1)) mblo(jv)=EXP(lgmgrid)
  IF (jv > 1) mbhi(jv-1)=EXP(lgmgrid)
END DO
DO jv=1,nbins
  mbmid(jv)=EXP(0.5*LOG(mblo(jv)*mbhi(jv))) ! geometric mean
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_ingridg
END MODULE ukca_ingridg_mod
