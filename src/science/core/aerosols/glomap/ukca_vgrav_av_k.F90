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
!    Calculates the average (over kth moment) gravitational
!    settling velocity for a log-normally distributed aerosol
!    population with geometric mean diameter DPG and geometric
!    mean standard deviation SIGMA.
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
MODULE ukca_vgrav_av_k_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_VGRAV_AV_K_MOD'

CONTAINS

SUBROUTINE ukca_vgrav_av_k(nbox,k,dp,sigma,dvisc,mfpa,rhop,                    &
                           vgrav_av)
!------------------------------------------------------------------
!
! Purpose
! -------
! Calculates the average (over kth moment) gravitational
! settling velocity for a log-normally distributed aerosol
! population with geometric mean diameter DP and geometric
! mean standard deviation SIGMA following method in Regional
! Particulate Model as described by Binkowski & Shankar (1995).
!
! Parameters
! ----------
! None
!
! Inputs
! ------
! NBOX       : Number of boxes in domain
! K          : Index of moment for calculation
! DP         : Geometric mean particle diameter for mode (m)
! SIGMA      : Geometric standard deviation for mode
! DVISC      : Dynamic viscosity of air (kg m-1 s-1)
! MFPA       : Mean free path of air (m)
! GG         : Gravitational acceleration = (ms^-2)
! RHOP       : Density of aerosol particle (kgm^-3)
!
! Outputs
! -------
! VGRAV_AV   : Avg. grav. settling velocity (m s-1)
!
! Local variables
! ---------------
! KNG        : Knudsen number for geo. mean sized particle
! PREF       : Prefactor term to expression
! LNSQSG     : ln(SIGMA)*ln(SIGMA)
! TERM1,TERM2: Terms in average diff coeff expression
!
! Inputted by module
! ------------------
! GG        : Gravitational acceleration (ms^-2)
!
!--------------------------------------------------------------------
USE ukca_um_legacy_mod,  ONLY: gg => g
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: k
INTEGER, INTENT(IN) :: nbox
REAL,    INTENT(IN) :: dp(nbox)
REAL,    INTENT(IN) :: sigma
REAL,    INTENT(IN) :: dvisc(nbox)
REAL,    INTENT(IN) :: mfpa(nbox)
REAL,    INTENT(IN) :: rhop(nbox)

REAL,    INTENT(OUT) :: vgrav_av(nbox)

! Local variables
REAL :: kng(nbox)
REAL :: pref(nbox)
REAL :: lnsqsg
REAL :: term1
REAL :: term2
REAL :: term3
REAL :: term4

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_VGRAV_AV_K'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

lnsqsg=LOG(sigma)*LOG(sigma)
term1=(4.0*REAL(k)+4.0)/2.0
term2=(2.0*REAL(k)+1.0)/2.0
term3=EXP(term1*lnsqsg)
term4=1.246*EXP(term2*lnsqsg)

kng(:)=2.0*mfpa(:)/dp(:)
pref(:)=rhop(:)*dp(:)*dp(:)*gg/(18.0*dvisc(:))
vgrav_av(:)=pref(:)*(term3+term4*kng(:))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_vgrav_av_k
END MODULE ukca_vgrav_av_k_mod
