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
!    Solve the coagulation/nucleation equation.
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
MODULE ukca_solvecoagnucl_v_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_SOLVECOAGNUCL_V_MOD'

CONTAINS

SUBROUTINE ukca_solvecoagnucl_v(nv,mask,a,b,c,nd,dtz,deln)
!---------------------------------------------------------------
!
! Purpose
! -------
! Solve the coagulation/nucleation equation which is of the form
! dN/dt=A*N^2 + B*N + C
!
! Note that the A term refers to intra-modal coagulation
!           the B term refers to inter-modal coagulation
!           the C term refers to new particle formation
!
! Rearrange this to be in the form of an indefinite integral
!
! int_{N_0}^{N} dx/X = int_{t_0}^{t_0+deltat} dt = deltat
!
! where X= A*x^2 + B*x + C
!
! Use analytical solutions to solve this indefinite integral for
! different cases below following pg 28 of Bronshtein & Semendyayev,
! "Handbook of Mathematics" (1997, 3rd edition) Equation 40.
!
! D = -discriminant = 4*A*C - B^2
!
! 1) [LOGIC1A is true] : A /= 0, and D < 0
!
! int dx/X = {2/sqrt(D)} * arctan{ (2*A*x+B)/sqrt(D) }
!
! 2) [LOGIC1B is true] : A /= 0, and D > 0
!
! int dx/X = {1/sqrt(-D)}*log{ (2*A*x+B-sqrt(-D)) / (2*A*x+B+sqrt(-D) }
!
! 3) [LOGIC1CA is true] : A /= 0, D = 0 and B/=0 (causes error message)
!
! In this case, B^2 = 4AC (B=/0) and there is no solution given.
!
! 4) [LOGIC1CB is true] : A /= 0, D = 0 and B=0
!
! In this case, B^2 = 4AC, and B=0, so C must also = 0
! in which case we have dN/dt = A*N^2 which has solution
!
!       N=1/(1/N0-3*A*deltat)
!
! 5) [LOGIC2A is true]  : A == 0 and B /= 0
!
!    In this case dN/dt = B*N + C which has the solution
!
!       N={(B*N0+C)*exp(B*deltat)-C}/B
!
!     n.b. if no nucleation (C=0) then N=N0*exp{B*deltat}
!
! 6) [LOGIC2B is true]  : A == 0 and B == 0
!
!    In this case dN/dt = C and N=N0+C*deltat.
!
! Parameters
! ----------
! None
!
! Inputs:
! ------
! A      : constant in coag/nucl equation
! B      : constant in coag/nucl equation
! C      : constant in coag/nucl equation
! ND     : initial aerosol particle no. conc. (ptcls/cc)
! NV     : No of points
! DTZ    : nucl/coag/cond time step (s)
!
! Outputs:
! -------
! DELN   : change in ND due to nucl & coag (ptcls/cc)
!
!
! Local variables:
! ---------------
! D      : 4*A*C - B^2
! SQD    : SQRT(D) if D>0 or SQRT(-D) if D<0
! NDNEW  : Number concentration after combined coag-nucl.
! LOGIC1,LOGIC2,LOGIC3,etc. : Logicals for various cases (see above)
! IERR   : stores which boxes had no solution (should not occur).
!
!--------------------------------------------------------------------

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE ukca_types_mod, ONLY: integer_32, logical_32
USE umPrintMgr, ONLY:                                                          &
    umPrint,                                                                   &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


!     Subroutine interface
INTEGER, INTENT(IN) :: nv
REAL, INTENT(IN)    :: a(nv)
REAL, INTENT(IN)    :: b(nv)
REAL, INTENT(IN)    :: c(nv)
REAL, INTENT(IN)    :: nd(nv)
REAL, INTENT(IN)    :: dtz
LOGICAL (KIND=logical_32), INTENT(IN) :: mask(nv)
REAL, INTENT(IN OUT) :: deln(nv)

!     Local variables
INTEGER :: i
INTEGER (KIND=integer_32) :: ierr(nv)
LOGICAL (KIND=logical_32) :: logic1(nv)
LOGICAL (KIND=logical_32) :: logic2(nv)
LOGICAL (KIND=logical_32) :: logic3(nv)
LOGICAL (KIND=logical_32) :: logic1a(nv)
LOGICAL (KIND=logical_32) :: logic1b(nv)
LOGICAL (KIND=logical_32) :: logic1c(nv)
LOGICAL (KIND=logical_32) :: logic1ca(nv)
LOGICAL (KIND=logical_32) :: logic1cb(nv)
LOGICAL (KIND=logical_32) :: logic2a(nv)
LOGICAL (KIND=logical_32) :: logic2b(nv)
LOGICAL (KIND=logical_32) :: logic1aok(nv)
LOGICAL (KIND=logical_32) :: logic1aok2(nv) ! logical for extra robustness
REAL    :: ndnew(nv)
REAL    :: d(nv)
REAL    :: sqd(nv)
REAL    :: sqd_tms_dtz(nv) ! intermediate array within extra robustness code
REAL    :: term1(nv)
REAL    :: term2(nv)
REAL    :: term3(nv)
REAL    :: term4(nv)
REAL    :: eps_ab
REAL    :: eps_d
CHARACTER(LEN=errormessagelength) :: cmessage        ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SOLVECOAGNUCL_V'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

eps_ab=1.0e-20
eps_d=eps_ab*eps_ab

ierr(:) = 0

ndnew(:)=nd(:)  ! initialise NDNEW to ND rather than 0
d(:)=4.0*a(:)*c(:)-b(:)*b(:)

logic1(:)=mask(:) .AND. ((a(:) > eps_ab) .OR. (a(:) < -eps_ab))
logic2(:)=mask(:) .AND. ((a(:) < eps_ab) .AND. (a(:) > -eps_ab))
logic3(:)=((b(:) <= eps_ab) .AND. (b(:) >= -eps_ab))

logic1a(:)=logic1(:) .AND. (d(:) < -eps_d)
logic1b(:)=logic1(:) .AND. (d(:) >  eps_d)
logic1c(:)=logic1(:) .AND. ((d(:) <= eps_d) .AND. (d(:) >= -eps_d))

logic1ca(:)=logic1c(:) .AND. (.NOT. logic3(:))
logic1cb(:)=logic1c(:) .AND. (logic3(:))

logic2a(:)=logic2(:) .AND. (.NOT. logic3(:))
logic2b(:)=logic2(:) .AND. logic3(:)

! Below is for case where A /= 0 & D < 0
WHERE (logic1a(:))
  sqd(:)=SQRT(-d(:))
  term1(:)=2.0*a(:)*nd(:)+b(:)+sqd(:)
  term2(:)=2.0*a(:)*nd(:)+b(:)-sqd(:)
  term3(:)=term1(:)/term2(:)
ELSE WHERE
  term3(:) = 0.0e0
END WHERE

! lines below make sure we're not dividing by tiny TERM3
logic1aok(:)=logic1a(:) .AND.                                                  &
             ((term3(:) > eps_ab) .OR. (term3(:) < -eps_ab))

!    Added code here to give additional robustness of code to
!    avoid potential numerical problems where EXP(sqd(:)*dtz).
!    Code limits sqd(:)*dtz so that it does not exceed 50.
!    Done by introducing a new local array sqd_tms_dtz(nv)

DO i=1,nv
  IF (logic1aok(i)) THEN
    sqd_tms_dtz(i)=sqd(i)*dtz
    IF (sqd_tms_dtz(i) > 50.0) THEN
      sqd_tms_dtz(i)=50.0 ! limit to 50
      !                     (EXP(X) is INFTY for X>=50)
    END IF
    term4(i)=1.0-EXP(sqd_tms_dtz(i))/term3(i)
  END IF
END DO

! lines below make sure we're not dividing by tiny TERM4 (as for TERM3)
DO i=1,nv
  IF (logic1aok(i)) THEN
    IF ((term4(i) > eps_ab) .OR. (term4(i) < -eps_ab)) THEN
      logic1aok2(i) = .TRUE.
    ELSE
      logic1aok2(i) = .FALSE.
    END IF
  ELSE
    logic1aok2(i) = .FALSE.
  END IF
END DO


WHERE (logic1aok2(:))
  ndnew(:)=(2.0*sqd(:)/term4(:)-b(:)-sqd(:))/2.0/a(:)
END WHERE

! Below is for case where A /= 0 & D > 0
WHERE (logic1b(:))
  sqd(:)=SQRT(d(:))
  term1(:)=(2.0*a(:)*nd(:)+b(:))/sqd(:)
  term2(:)=ATAN(term1(:))+sqd(:)*dtz/2.0
  ndnew(:)=(sqd(:)*TAN(term2(:))-b(:))/2.0/a(:)
END WHERE

! Below is for case where A /= 0, D = 0 and B/=0 (error)
WHERE (logic1ca(:)) ierr(:)=1 ! ! A /= 0, D=0 & B=0
!
! Below is for case where A /= 0, D = 0 and B=0
WHERE (logic1cb(:)) ndnew(:)=1.0/(1.0/nd(:)-3.0*a(:)*dtz)

! Below is for case where A == 0, B /= 0 (==> D<0 )
WHERE (logic2a(:))
  term1(:)=EXP(b(:)*dtz)
  term2(:)=b(:)*nd(:)+c(:)
  ndnew(:)=(term1(:)*term2(:)-c(:))/b(:)
END WHERE

! Below is for case where A == 0, B == 0 (==> D=0 )
WHERE (logic2b(:)) ndnew(:)=nd(:)+c(:)*dtz

! Calculate increment to ND due to combined coag & nucl., only set DELN to
! calculated value using MASK criterion
deln(:)=0.0
WHERE (mask(:)) deln(:)=ndnew(:)-nd(:)

IF (SUM(ierr(:)) > 0) THEN
  DO i=1,nv
    IF (ierr(i) == 0) CYCLE
    d(i)=4.0*a(i)*c(i)-b(i)*b(i)
    WRITE(umMessage,'(A)') 'D=0, A not = 0, B not = 0, C not = 0'
    CALL umPrint(umMessage,src='ukca_solvecoagnucl_v')
    WRITE(umMessage,'(A,F0.5,A,F0.5,A,F0.5,A,F0.5)')                           &
                                'so B^2=4AC>0',a(i),' ',b(i),' ',c(i),' ',d(i)
    CALL umPrint(umMessage,src='ukca_solvecoagnucl_v')
    ! Note that if C=0 & D=0 then B=0
  END DO
  cmessage=' IERR > ZERO, D = 0'
  WRITE(umMessage,'(A)') cmessage
  CALL umPrint(umMessage,src='ukca_solvecoagnucl_v')

  CALL ereport('UKCA_SOLVECOAGNUCL_V',i,'IERR>0')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_solvecoagnucl_v
END MODULE ukca_solvecoagnucl_v_mod
