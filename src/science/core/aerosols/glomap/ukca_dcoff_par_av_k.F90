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
!    Calculates the average (over kth moment) diffusion coefficient
!    for a log-normally distributed aerosol population with
!    geometric mean diameter DG, geometric mean standard
!    deviation SSIGMA.
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
MODULE ukca_dcoff_par_av_k_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_DCOFF_PAR_AV_K_MOD'

CONTAINS

SUBROUTINE ukca_dcoff_par_av_k(nbox,k,dp,sigma,t,dvisc,mfpa,                   &
                               dcoff_par_av_k)
!--------------------------------------------------------------------
!
! Purpose
! -------
! Calculates the average (over kth moment) diffusion coefficient
! for a log-normally distributed aerosol population with
! geometric mean diameter DG, geometric mean standard
! deviation SIGMA following method in Regional Particulate
! Model as described by Binkowski & Shankar (1995).
! Also follows expression in Seinfeld & Pandis pg 474
! (Stokes-Einstein relation with slip-flow correction)
!
! Parameters
! ----------
! None
!
! Inputs
! ------
! NBOX           : Number of boxes in domain
! K              : Index of moment for calculation
! DP             : Geometric mean particle diameter for mode (m)
! SIGMA          : Geometric standard deviation for mode
! T              : Temperature of air (K)
! DVISC          : Dynamic viscosity of air (kg m-1 s-1)
! MFPA           : Mean free path of air (m)
!
! Outputs
! -------
! DCOFF_PAR_AV_K : Avg. ptcl diffusion coefficient (m^2 s-1)
!
! Local variables
! ---------------
! KNG            : Knudsen number for geo. mean sized particle
! LNSQSG         : ln(SIGMA)*ln(SIGMA)
! PREF           : Prefactor term to expression
! TERM1,TERM2    : Terms in average diff coeff expression
!
! Inputted by module CHEMISTRY_CONSTANTS
! ---------------------------------
! BOLTZMANN      : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
!
!--------------------------------------------------------------------
USE ukca_um_legacy_mod,  ONLY: pi
USE chemistry_constants_mod, ONLY: boltzmann
USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim
IMPLICIT NONE
!
! Subroutine interface
INTEGER, INTENT(IN) :: k
INTEGER, INTENT(IN) :: nbox
REAL, INTENT(IN)    :: dp(nbox)
REAL, INTENT(IN)    :: sigma
REAL, INTENT(IN)    :: t(nbox)
REAL, INTENT(IN)    :: dvisc(nbox)
REAL, INTENT(IN)    :: mfpa(nbox)
REAL, INTENT(OUT)   :: dcoff_par_av_k(nbox)

! Local variables
REAL    :: lnsqsg
REAL    :: term1
REAL    :: term2
REAL    :: term3
REAL    :: term4
REAL    :: kng(nbox)
REAL    :: pref(nbox)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DCOFF_PAR_AV_K'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

lnsqsg=LOG(sigma)*LOG(sigma)
term1=(-2.0*REAL(k)+1.0)/2.0
term2=(-4.0*REAL(k)+4.0)/2.0
term3=EXP(term1*lnsqsg)
term4=1.246*EXP(term2*lnsqsg)

kng(:)=2.0*mfpa(:)/dp(:)
pref(:)=boltzmann*t(:)/(3.0*pi*dvisc(:)*dp(:))
dcoff_par_av_k(:)=pref(:)*(term3+term4*kng(:))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_dcoff_par_av_k
END MODULE ukca_dcoff_par_av_k_mod
