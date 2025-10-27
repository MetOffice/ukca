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
!    Uses the Newton-Raphson method to solve the Backward Euler equations to
!    find the species concentration at the next timestep.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!   Called from asad_spmjpdriv
!
!     Purpose.
!     --------
!     To solve the chemical rate equations using the backward Euler method.
!
!     Interface
!     ---------
!     Called from chemistry driver routine via the mjpdriv driver.
!
!     Method.
!     -------
!     This routine solves the ODE equation:
!     df / dt = fdot where f is the species concentrations.
!
!     This equation is rearranged and solved using the Backward-Euler method:
!     (f(t=n+1) - f(t=n))/dt - fdot = 0
!     We let function G(f) = (f(t=n+1) - f(t=n))/dt - fdot and use
!     Newton-Raphson to find f at t=n+1 such that G(f)=0.
!     An iterative method is used to find f with the first guess f_0 found using
!     the forward Euler method.
!     Subsequent estimates, f_1, f_2, f_3,...,f_k,... are found for f(t=n+1).
!
!     The solver converges in two possible ways:
!     1. The residual relative error |G(f)| is less than a specified tolerance
!     2. The relative error |f_k+1 - f_k| is less than a specified tolerance.
!
!     The solver may exit the routine without a solution if
!     1. the species become too large or
!     2. convergence has not been achieved
!     In each of these cases the timestep is halved.
!
!     Currently a limiter is applied to the first 7 iterations to prevent very
!     rapid initial changes leading to divergence.
!     A limiter is also applied at all iterations to ensure that the chemical
!     species do not become smaller than f_min.
!
!     A quasi-Newton solver is typically done on iterations 2 and 3 of the NR
!     solver. This uses a cheaper estimate of the Jacobian matrix.
!
!     Bit-comparability: The error norms used in this routine are a
!     function of both the number of grid boxes (n_points) and the number of
!     chemical species (jpcspf).
!     Changing n_points (or jpcspf) will change the results of the error norm
!     and hence results will not bit compare for different values of n_points.
!     This is seen when changing either the PE domain decomposition when
!     using horizontal slices or the chunk size if using vertical columns.
!
!     Global variables
!     ----------------
!     RelTol_residual_error - tolerance - set to  1.0E-10
!     rafmin - limit for first few iterations, set to 0.1
!     rafmax - limit for first few iterations, set to 1.0E+04
!     RelTol_error - tolerance - set to 1.0E-04
!     f_max - maximum concentration, above which divergence is assumed
!     f_min - smallest non-zero concentration
!
!     Local variables
!     ---------------
!     ifi           Number of ftoy iterations.
!     f_initial     Values of f at start of chemistry step.
!     damp1         Damping factor to apply to the first iteration
!     deltt         Reciprocal time step length  (1./cdt)
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE asad_spimpmjp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_SPIMPMJP_MOD'

CONTAINS

! *********************************************************************

SUBROUTINE forward_euler(n_points, f, f_initial, f_min, nonzero_map, spfj)

USE asad_mod,            ONLY: cdt, jpcspf, fdot, spfjsize_max
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN)  :: n_points
REAL, INTENT(OUT)    :: f(1:n_points,1:jpcspf)
REAL, INTENT(IN)     :: f_initial(1:n_points,1:jpcspf)
REAL, INTENT(IN)     :: f_min
INTEGER, INTENT(IN)  :: nonzero_map(1:jpcspf,1:jpcspf)
REAL, INTENT(IN)     :: spfj(1:n_points,spfjsize_max)

INTEGER :: jl, jtr, ip

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*),   PARAMETER :: RoutineName='FORWARD_EULER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Solve for all species in f array, including non-transported RO2 species
DO jtr=1,jpcspf
  ip = nonzero_map(jtr,jtr)
  DO jl=1,n_points
    IF (spfj(jl,ip) > 0.0) THEN
      f(jl,jtr) = f_initial(jl,jtr) + cdt*fdot(jl,jtr)
    ELSE
      f(jl,jtr) = f_initial(jl,jtr) + (cdt*fdot(jl,jtr))/(1.0-cdt*spfj(jl,ip))
    END IF
    IF (f(jl,jtr) < f_min) f(jl,jtr) = f_min
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE forward_euler

! *********************************************************************

SUBROUTINE calc_residual_error(n_points,residual_error,G_f,f_min)

USE asad_mod,            ONLY: jpcspf, nlf
USE asad_mod,            ONLY: prod, slos
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_points
REAL, INTENT(OUT) :: residual_error
REAL, INTENT(IN) :: G_f(1:n_points,1:jpcspf)
REAL, INTENT(IN) :: f_min

REAL :: tmprc(1:n_points,1:jpcspf)
INTEGER :: jl, jtr, j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='CALC_RESIDUAL_ERROR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

residual_error = 0.0

!  Temporary prod+loss array used to test for convergence
!  by calculating residual_error further down
DO jtr=1,jpcspf
  j=nlf(jtr)
  tmprc(1:n_points,jtr) = prod(1:n_points,j) + slos(1:n_points,j)
END DO

DO jtr=1,jpcspf
  DO jl=1,n_points
    IF (ABS(tmprc(jl,jtr)) > f_min) THEN
      residual_error=MAX(residual_error,ABS(G_f(jl,jtr)/tmprc(jl,jtr)))
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_residual_error

! *********************************************************************

SUBROUTINE calc_error_norm(n_points,error_norm,f,f_incr,f_min)

USE asad_mod,            ONLY: jpcspf
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_points
REAL, INTENT(OUT)   :: error_norm
REAL, INTENT(IN)    :: f(1:n_points,1:jpcspf)
REAL, INTENT(IN)    :: f_incr(1:n_points,1:jpcspf)
REAL, INTENT(IN)    :: f_min

INTEGER :: jl
INTEGER :: jtr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='CALC_ERROR_NORM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_norm = 0.0
DO jtr=1,jpcspf
  DO jl=1,n_points
    IF (ABS(f_incr(jl,jtr)) > 1.0e-16) THEN
      error_norm = MAX(error_norm,ABS(f_incr(jl,jtr)/MAX(f(jl,jtr),f_min)))
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_error_norm

! *********************************************************************

SUBROUTINE asad_spimpmjp(exit_code, ix, jy, nlev, n_points, location,          &
                         solver_iter)

USE asad_mod,           ONLY: ptol, peps, cdt, f, fdot, nitnr, nstst, y,       &
                              fj, nonzero_map, ltrig, jpcspf, spfj,            &
                              modified_map, nonzero_map_unordered
USE asad_sparse_vars,   ONLY: setup_spfuljac, spfuljac, spresolv2, splinslv2
USE ukca_config_specification_mod, ONLY: ukca_config
USE yomhook,            ONLY: lhook, dr_hook
USE parkind1,           ONLY: jprb, jpim
USE umPrintMgr,         ONLY: printstatus, PrStatus_Oper, PrStatus_Diag,       &
                              umMessage, umPrint
USE errormessagelength_mod, ONLY: errormessagelength
USE ukca_um_legacy_mod, ONLY: mype
USE asad_diffun_mod,    ONLY: asad_diffun
USE asad_fuljac_mod,    ONLY: asad_fuljac
USE asad_steady_mod,    ONLY: asad_steady
USE asad_ftoy_mod,      ONLY: asad_ftoy

IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: n_points
INTEGER, INTENT(IN) :: ix
INTEGER, INTENT(IN) :: jy
INTEGER, INTENT(IN) :: nlev
INTEGER, INTENT(IN) :: location
INTEGER, INTENT(OUT):: exit_code
INTEGER, INTENT(OUT):: solver_iter ! No. of iterations

! Local variables
INTEGER, PARAMETER :: maxneg=2000     ! Max No. negatives allowed
INTEGER :: iter
INTEGER :: jl
INTEGER :: jtr
INTEGER :: ifi
INTEGER :: i
INTEGER :: itr
INTEGER :: count_negatives
REAL :: ztmp
! The maximum concentration allowed was previously f_max = 1.0/f_min
REAL, PARAMETER :: f_max = 1.0e30
REAL :: f_min
! Negative threshold to turn off filtering in quasi-Newton mode
REAL, PARAMETER :: max_val = -1.0
REAL :: RelTol_residual_error
REAL :: RelTol_error
REAL :: rafmin
REAL :: rafmax
REAL :: deltt
REAL :: damp1
REAL :: error_norm
REAL :: residual_error

LOGICAL :: not_first_call = .FALSE.

REAL :: f_initial(n_points,jpcspf) ! Concentration from previous timestep,f(t=n)
REAL :: f_incr(n_points,jpcspf)  ! Increment at end of NR step that is added to
                                 ! the previous estimate of the chemical species
                                 ! f(t=n+1) at the next timestep.
REAL :: G_f(n_points,jpcspf)   ! Function G_f = (f(t=n+1) - f(t=n))/dt - fdot
                               ! Newton-Raphson finds f such that G_f(f) = 0

INTEGER, PARAMETER :: ltrig_iter=51   ! Set to nrsteps if want LTRIG
INTEGER, PARAMETER :: incr_limiter=7  ! The maximum number of iterations for
                                      !  which the increment limiter is applied.

!  Variables required for possible quasi-Newton step:
REAL :: G_f_old(n_points,jpcspf)
REAL :: G_ftmp(n_points,jpcspf)
REAL :: delta_G(n_points,jpcspf)
REAL :: coeff

REAL :: zsum(jpcspf) ! Diagnostic output, sum of f_initial across all points

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*),   PARAMETER :: RoutineName='ASAD_SPIMPMJP'

CHARACTER(LEN=errormessagelength) :: cmessage1 = "(1x,i2,20(1x,1pG12.4))"
CHARACTER(LEN=errormessagelength) :: cmessage2 = "(1x,a3,20(1x,1pG12.4))"

LOGICAL, SAVE :: first = .TRUE.
LOGICAL, SAVE :: first_pass = .TRUE.


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

f_min = SQRT(peps) ! Minimum species concentration

RelTol_residual_error = 1.0e-10 ! Relative tolerance for |G(f)|
RelTol_error = 10*ptol          ! Relative tolerance for |f_incr|

! Size of increment limiter
rafmin = 1.0e-01
rafmax = 1.0e+04

exit_code = 0
error_norm = 1.0
deltt = 1.0/cdt
damp1 = 0.5

!  Save values of f at start of step
f_initial = f
WHERE (f<f_min) f = f_min

! Call ASAD_STEADY at start of step to initialise deriv properly
IF (nstst /= 0)  CALL asad_steady( n_points )

! OMP CRITICAL will only allow one thread through this code at a time,
! while the other threads are held until completion.
!$OMP CRITICAL (setup_jacobian_init)
IF (first_pass) THEN
  IF (first) THEN
    ! Determine number and positions of nonzero elements in sparse
    ! full Jacobian
    CALL setup_spfuljac()
    first = .FALSE.
  END IF
  first_pass = .FALSE.
END IF
!$OMP END CRITICAL (setup_jacobian_init)

CALL spfuljac(n_points,cdt,f_min,nonzero_map,spfj)

! Call forward Euler to make first guess, f_0 for f(t=n+1) (next timestep)
CALL forward_euler(n_points, f, f_initial, f_min, nonzero_map, spfj)

! Start Newton-Raphson loop to generate estimates f_1, f_2, f_3,... for f(t=n+1)
DO iter=1,ukca_config%nrsteps

  ifi = 0
  IF (iter == 1) ifi=nitnr
  ! Asad_ftoy needs iteration count to calculate concentration of total RO2
  ! if l_ro2_perm_chem = TRUE
  CALL asad_ftoy(not_first_call, ifi, iter, n_points, ix, jy, nlev)

  IF (nstst /= 0 .AND. ifi ==0) CALL asad_steady( n_points )

  IF (ltrig .AND. printstatus >= prstatus_oper) THEN
    DO jl=1,n_points
      WRITE(umMessage,"('Point: ',i4)") jl
      CALL umPrint(umMessage,src='asad_spimpmjp')
      IF (iter == 1) THEN
        WRITE(umMessage,cmessage1) iter-1,(f_initial(jl,jtr),jtr=1,jpcspf)
        CALL umPrint(umMessage,src='asad_spimpmjp')
        WRITE(umMessage,cmessage1) iter-1,                                     &
                      ((f_initial(jl,jtr)+cdt*fdot(jl,jtr)),jtr=1,jpcspf)
        CALL umPrint(umMessage,src='asad_spimpmjp')
      END IF
      WRITE(umMessage,cmessage1) iter-1, (f(jl,jtr),jtr=1,jpcspf),             &
                                                                 (y(jl,i),i=1,2)
      CALL umPrint(umMessage,src='asad_spimpmjp')
    END DO
  END IF

  CALL asad_diffun( n_points ) ! calculates fdot

  IF (error_norm < RelTol_error) THEN
    exit_code = 0           ! Successful exit
    solver_iter = iter - 1  ! Technically converged on previous iteration
    GO TO 9999
  END IF

  IF (iter == ukca_config%nrsteps) THEN
    exit_code = 2
    GO TO 9999
  END IF

  ! Calculate G_f = (f(t=n+1) - f(t=n))/dt - fdot.
  ! Find f(t=n+1) such that G_f = 0.
  G_f = (f - f_initial)*deltt - fdot

  ! Calculate residual error (the relative magnitude of G_f)
  CALL calc_residual_error(n_points,residual_error,G_f,f_min)
  IF (residual_error < RelTol_residual_error) THEN
    exit_code = 0 ! Successful exit
    solver_iter = iter
    GO TO 9999
  END IF

  DO jl=1,n_points
    DO jtr=1,jpcspf
      IF (f(jl,jtr) > f_max) THEN
        ! Exit solver if chemical species is larger than f_max
        exit_code = 3
      END IF
    END DO
  END DO

  ! moved outside of loop to allow parallelisation
  IF (exit_code ==3) THEN
    GO TO 9999
  END IF

  CALL spfuljac(n_points,cdt,f_min,nonzero_map,spfj)

  IF (ltrig .AND. printstatus == PrStatus_Diag) THEN
    WRITE(umMessage,"('Iteration ',i4)") iter
    CALL umPrint(umMessage,src='asad_spimpmjp')
    DO jl=1,n_points
      WRITE(umMessage,"('Point: ',i4)") jl
      CALL umPrint(umMessage,src='asad_spimpmjp')
      fj(jl,:,:) = 0.0
      DO jtr=1,jpcspf
        DO itr=1,jpcspf
          IF (nonzero_map(jtr,itr) > 0)                                        &
            fj(jl,jtr,itr) = spfj(jl,nonzero_map(jtr,itr))
        END DO
      END DO
      DO jtr=1,jpcspf
        WRITE(umMessage,cmessage1) jtr, (fj(jl,jtr,itr),itr=1,jpcspf)
        CALL umPrint(umMessage,src='asad_spimpmjp')
      END DO
    END DO
  END IF

  CALL splinslv2(n_points,G_f,f_incr,f_min,f_max,nonzero_map_unordered,        &
                    modified_map,spfj)

  IF (ltrig .AND. printstatus == PrStatus_Diag) THEN
    DO jl=1,n_points
      WRITE(umMessage,"('Point: ',i4)") jl
      CALL umPrint(umMessage,src='asad_spimpmjp')
      WRITE(umMessage,cmessage2) 'G_f',(G_f(jl,jtr),jtr=1,jpcspf)
      CALL umPrint(umMessage,src='asad_spimpmjp')
      WRITE(umMessage,cmessage2) 'fdt',(fdot(jl,jtr),jtr=1,jpcspf)
      CALL umPrint(umMessage,src='asad_spimpmjp')
      WRITE(umMessage,cmessage2) 'del',                                        &
                            ((f(jl,jtr)-f_initial(jl,jtr))*deltt,jtr=1,jpcspf)
      CALL umPrint(umMessage,src='asad_spimpmjp')
      WRITE(umMessage,cmessage2) 'f  ',(f(jl,jtr),jtr=1,jpcspf)
      CALL umPrint(umMessage,src='asad_spimpmjp')
      WRITE(umMessage,cmessage2) 'f_incr',(f_incr(jl,jtr),jtr=1,jpcspf)
      CALL umPrint(umMessage,src='asad_spimpmjp')
    END DO
    CALL asad_fuljac(n_points)
    WRITE(umMessage,"('Iteration ',i4)") iter
    CALL umPrint(umMessage,src='asad_spimpmjp')
    DO jl=1,n_points
      DO jtr=1,jpcspf
        zsum(jtr) = 0.0
        DO itr=1,jpcspf
          zsum(jtr)=zsum(jtr)+fj(jl,jtr,itr)*f_incr(jl,itr)
        END DO
        WRITE(umMessage,cmessage1) jtr,                                        &
                                  (fj(jl,jtr,itr)*f_incr(jl,itr),itr=1,jpcspf)
        CALL umPrint(umMessage,src='asad_spimpmjp')
      END DO
      WRITE(umMessage,cmessage2) 'sum', (zsum(jtr),jtr=1,jpcspf)
      CALL umPrint(umMessage,src='asad_spimpmjp')
    END DO
  END IF

  ! Damp increment on first Newton-Raphson iteration
  IF (iter == 1) f_incr = damp1*f_incr

  !  Filter increments
  f_incr = MIN(MAX(f_incr,-f_max),f_max)
  CALL calc_error_norm(n_points,error_norm,f,f_incr,f_min)

  ! Apply increment f_k+1 = f_k + f_incr
  count_negatives = 0
  DO jtr=1,jpcspf
    DO jl=1,n_points
      !  New mixing ratios
      ztmp = f(jl,jtr) + f_incr(jl,jtr)
      !  Put limit on increment for first few iterations
      IF (iter < incr_limiter) THEN
        ztmp = MAX(rafmin*f(jl,jtr),MIN(rafmax*f(jl,jtr),ztmp))
      END IF

            !  Filter negatives and zeros
      IF (ztmp == 0.0) ztmp = f_min
      IF (ztmp < 0.0) THEN
        ztmp = f_min
        count_negatives = count_negatives + 1
      END IF
      !  Final mixing ratios
      f(jl,jtr) = ztmp
    END DO
  END DO

  IF (count_negatives > maxneg) THEN
    ! Exit if number of negative values exceeds threshold
    exit_code = 2
    GO TO 9999
  END IF

  ! Perform quasi-Newton (Broyden) Method to reduce number of iterations
  ! This is done on iterations 2 <= iter <= 50, and recommended on steps 2 & 3
  ! This step will not be done if error_norm < RelTol_error, i.e. the values are
  ! converged and the routine is about to exit (this is actually tested
  ! at the top of the next loop).
  IF (ukca_config%l_ukca_quasinewton .AND. (error_norm >= RelTol_error)) THEN
    IF ((iter >= ukca_config%i_ukca_quasinewton_start) .AND.                   &
        (iter <= ukca_config%i_ukca_quasinewton_end)) THEN

      CALL asad_ftoy(not_first_call,ifi, iter, n_points, ix, jy, nlev)
      CALL asad_steady(n_points)
      CALL asad_diffun(n_points) ! updates fdot

      G_f_old = G_f
      G_f     = (f - f_initial)*deltt - fdot
      delta_G = G_f - G_f_old

      DO jl=1,n_points
        coeff = DOT_PRODUCT(G_f(jl,:), delta_G(jl,:))                          &
                    /DOT_PRODUCT(delta_G(jl,:),delta_G(jl,:))
        G_ftmp(jl,:) = G_f(jl,:)*(1.0 - coeff)
      END DO

      CALL spresolv2(n_points,G_ftmp,f_incr,f_min,modified_map,spfj,max_val)

      f = f + f_incr
      ! remove negative values. Does not need to be done in
      ! as intelligent a way as above, as we are not exiting
      ! the routine directly after this step.
      f = ABS(f)
    END IF
  END IF ! l_ukca_quasinewton

END DO

9999 CONTINUE

IF (exit_code /= 0) THEN
  ! Solver has not found a solution.
  f = f_initial
  solver_iter = iter

  IF (count_negatives > maxneg) THEN
    WRITE(umMessage,"('Negatives - exceeds maxneg')")
    CALL umPrint(umMessage,src='asad_spimpmjp')
    WRITE(umMessage,                                                           &
      "(1x,'Too many negatives (>',i4,') in spimpmjp (iter',i3,"//             &
      "')  lon=',i3,'  lat=',i3,'; halving step')")                            &
      maxneg,iter, location, mype
    CALL umPrint(umMessage,src='asad_spimpmjp')
  END IF

  IF (iter >= ltrig_iter) THEN
    exit_code = 4 ! exit with debug option for use in asad_spmjpdriv
    WRITE(umMessage,                                                           &
"('Convergence problems (',i3,1x,'iter) at location=',i3,' pe=',i3)")          &
    iter, location, mype
    CALL umPrint(umMessage,src='asad_spimpmjp')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_spimpmjp

END MODULE asad_spimpmjp_mod
