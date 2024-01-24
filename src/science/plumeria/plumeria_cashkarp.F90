! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module containing the subroutines CASHKARPQS and CASHKARP.
!
!  Methods:
!    Description of each subroutine contained here:
!    1. PLUMERIA_CASHKARPQS: Integrates properties with elevation using
!       the Cash-Karp
!       Method: Cash, J.R. and Karp, A.H., A variable order Runge-Kutta method
!               for initial value problems with rapidly varying right-hand
!               sides, ACM Transactions on Mathematical Software,
!               v. 16, no. 3, pp. 201-222.
!    2. PLUMERIA_CASHKARP: Cash-carp ode solver, uses the cash-karp method of
!       solving ode's using the fifth-order Runge-Kutta solution with
!       embedded fourth-order solution:
!       y_out(i) = y_in(i) + c(1)* dydx(i) + c(2)*dydx2(i) + c(3)*dydx3(i)
!                  + c(4)*dydx4(i) + c(5)*dydx5(i) + c(6)*dydx6(i)
!       where c(i) is a constant with values given in Cash and Karp,
!       and dydx, dydx2 etc. are the slopes of y with x evaluated
!       at the following intermediate values of x: x+h/5, x+3h/10,
!       x+3h/5, x+h, 1+7h/8.

!       This is compared with the fourth-order Runge-Kutta solution:
!       ystar_out(i) = y_in(i) + cstar(1)* dydx(i) + cstar(2)*dydx2(i) +
!                      cstar(3)*dydx3(i) + cstar(4)*dydx4(i) +
!                      cstar(5)*dydx5(i) + cstar(6)*dydx6(i)
!       where cstar(i) are coefficients with different values.
!
!       The differences between the fourth-order and fifth-order solutions
!       are used in the quality step adjustment in cashkarpqs.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE plumeria_cashkarp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PLUMERIA_CASHKARP_MOD'

CONTAINS

!=======================================================================!

SUBROUTINE plumeria_cashkarpqs(y, dydx, x, x_out, h_in, yscal, h_out, hnext)
  ! ---------------------------------------------------------------------
  ! Description:
  !   Integrates properties with elevation using the Cash-Karp.
  ! ---------------------------------------------------------------------

USE plumeria_derivs_mod,       ONLY: plumeria_derivs
USE parkind1,                  ONLY: jpim, jprb      ! DrHook
USE yomhook,                   ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Arguments
REAL, INTENT(IN) :: dydx(13)  ! Input data for integration
REAL, INTENT(IN) :: x         ! Initial x value for integration
REAL, INTENT(IN) :: h_in      ! Initial step size for integration
REAL, INTENT(IN) :: yscal(13) ! Scaling factors for variables
REAL, INTENT(IN OUT) :: y(13) ! Output data after integration
REAL, INTENT(OUT) :: x_out    ! Updated value of x after integration
REAL, INTENT(OUT) :: h_out    ! Stepsize used in current integration step
REAL, INTENT(OUT) :: hnext    ! Stepsize for next integration step

! Return values from calling plumeria_cashkarp
REAL :: yerr(13)       ! Store estimated errors at each integration step
REAL :: ytemp(13)      ! Temporary storage of values during integration

! Local variables
! if yerr/yscal exceeds the number below, reduce step size.
REAL, PARAMETER   :: eps_small = 5.0e-4
REAL :: alpha, max_error, h, xnew

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PLUMERIA_CASHKARPQS'
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

h = h_in
CALL plumeria_cashkarp(y, dydx, x, h, ytemp, yerr)
max_error = 0.0
!Find maximum error, scaled to scale factors
max_error = MAXVAL(ABS(yerr / yscal))
max_error = max_error / eps_small

check_maxerror: DO
  IF (max_error > 1.0) THEN
    !adjustment following Dormund, J.R.
    alpha = max_error ** 0.25
    h = 0.9 * h / MIN(alpha, 10.0)
    xnew = x + h
    CALL plumeria_cashkarp(y, dydx, x, h, ytemp, yerr)
    max_error = 0.0
    !Find maximum error, scaled to scale factors
    max_error = MAXVAL(ABS(yerr / yscal))
    max_error = max_error / eps_small

  ELSE
    EXIT check_maxerror
  END IF
END DO check_maxerror

!Increase h, but not by more than 5 times
hnext = MIN(0.9 * h / (max_error ** 0.2), 5.0 * h)
h_out = h
x_out = x + h
y = ytemp
RETURN

END SUBROUTINE plumeria_cashkarpqs

!=======================================================================!

SUBROUTINE plumeria_cashkarp(y_in, dydx, x, h, yout, yerr)
  ! ---------------------------------------------------------------------
  ! Description:
  !   Cash-carp ode solver.
  ! ---------------------------------------------------------------------

USE plumeria_derivs_mod,   ONLY: plumeria_derivs
USE parkind1,              ONLY: jpim, jprb      ! DrHook
USE yomhook,               ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Arguments
REAL, INTENT(IN) :: y_in(13)       !input values
REAL, INTENT(IN) :: x              !initial x value
REAL, INTENT(IN) :: h              !step size
REAL, INTENT(OUT):: yout(13)       !output y
REAL, INTENT(OUT):: yerr(13)       !Difference between fourth-order
                                   !and fifth-order solutions

! Local variables for cashcarp ode solver
REAL          :: dydx(13), dydx2(13), dydx3(13), dydx4(13),                    &
                 dydx5(13), dydx6(13), ytemp(13), ystar_out(13)

!Constants
!Constants taken from p. 206 of Cash and Karp (ACM Transactions of on
!Mathematical Software, v. 16, pp. 201-222)
REAL, PARAMETER  :: b_21    = .2
REAL, PARAMETER  :: b_31    = 3.0/40.0
REAL, PARAMETER  :: b_32    = 9.0/40.0
REAL, PARAMETER  :: b_41    = .3
REAL, PARAMETER  :: b_42    = -.9
REAL, PARAMETER  :: b_43    = 1.2
REAL, PARAMETER  :: b_51    = -11.0/54.0
REAL, PARAMETER  :: b_52    = 2.5
REAL, PARAMETER  :: b_53    = -70.0/27.0
REAL, PARAMETER  :: b_54    = 35.0/27.0
REAL, PARAMETER  :: b_61    = 1631.0/55296.0
REAL, PARAMETER  :: b_62    = 175.0/512.0
REAL, PARAMETER  :: b_63    = 575.0/13824.0
REAL, PARAMETER  :: b_64    = 44275.0/110592.0
REAL, PARAMETER  :: b_65    = 253.0/4096.0
REAL, PARAMETER  :: c_1     = 37.0/378.0
REAL, PARAMETER  :: c_3     = 250.0/621.0
REAL, PARAMETER  :: c_4     = 125.0/594.0
REAL, PARAMETER  :: c_6     = 512.0/1771.0
REAL, PARAMETER  :: cstar_1 = 2825.0/27648.0
REAL, PARAMETER  :: cstar_3 = 18575.0/48384.0
REAL, PARAMETER  :: cstar_4 = 13525.0/55296.0
REAL, PARAMETER  :: cstar_5 = -277.0/14336.0
REAL, PARAMETER  :: cstar_6 = 1.0/4.0
REAL, PARAMETER  :: dc_1    = c_1-2825.0/27648.0
REAL, PARAMETER  :: dc_3    = c_3-18575.0/48384.0
REAL, PARAMETER  :: dc_4    = c_4-13525.0/55296.0
REAL, PARAMETER  :: dc_5    = 277.0/14336.0
REAL, PARAMETER  :: dc_6    = c_6-.25

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PLUMERIA_CASHKARP'
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!Calculate y and dydx at intermediate values of x
ytemp = y_in + b_21 * h * dydx

!at h/5
CALL plumeria_derivs(x + h / 5.0, ytemp, dydx2)
ytemp = y_in + h *(b_31 * dydx + b_32 * dydx2)

!at 3h/10
CALL plumeria_derivs(x + 3.0 * h / 10.0, ytemp, dydx3)
ytemp = y_in + h * (b_41 * dydx + b_42 * dydx2 + b_43 * dydx3)

!at 6h/10
CALL plumeria_derivs(x + 6.0 * h / 10.0, ytemp, dydx4)
ytemp = y_in + h * (b_51 * dydx + b_52 * dydx2 + b_53 * dydx3 + b_54 * dydx4)

!at h
CALL plumeria_derivs(x + h, ytemp, dydx5)
ytemp = y_in + h * (b_61 * dydx + b_62 * dydx2 + b_63 * dydx3 + b_64 * dydx4   &
        + b_65 * dydx5)

! at 7h/8
CALL plumeria_derivs(x + 0.875 * h, ytemp, dydx6)

!Fifth order Runge-Kutta formula
yout = y_in + h * (c_1 * dydx + c_3 * dydx3 + c_4 * dydx4 + c_6 * dydx6)

!Embedded fourth-order Runge-Kutta formula
ystar_out = y_in + h * (cstar_1 * dydx + cstar_3 * dydx3 +                     &
            cstar_4 * dydx4 + cstar_5 * dydx5 + cstar_6 * dydx6)

!Difference between fourth-order and fifth-order solutions
yerr = yout - ystar_out

RETURN

END SUBROUTINE plumeria_cashkarp
END MODULE plumeria_cashkarp_mod
