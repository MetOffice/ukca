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
!    Driver for fully implicit ODE integrator.
!    Part of the ASAD chemical solver.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!   Called from ASAD_CDRIVE
!
!
!     MJPDRIV  - Driver for MJP fully implicit ODE integrator.
!
!     Michael Prather            Earth System Science
!     Oliver Wild                University of California, Irvine
!
!     ASAD: mjpdriv              Version: mjpdriv.f 1.0 04/17/97
!
!     Purpose.
!     --------
!     To organise the integration of the chemical rate equations using
!     the MJP implicit integrator.
!
!     Interface
!     ---------
!     Called from chemistry driver routine *cdrive*.
!
!     This routine assumes that all the bi-,tri-, phot- and het-
!     reaction rates have been computed prior to this routine.
!     It is also assumed that the species array, y, has been set
!     from the array, f, passed to ASAD, and that constant species
!     have been set. This can be done by calling the routine fyinit.
!
!     Method.
!     -------
!     This routine calls the MJP integrator once for each gridpoint
!     of the one-dimensional arrays passed to it.
!     If convergence isn't achieved, the time step is halved, and the
!     integrator called again - this is continued until either
!     convergence is achieved or the minimum time step length is
!     encountered (currently 1.E-05 seconds).
!
!     Local variables
!     ---------------
!     ncst    -  Stores number of basic chemical steps 'ncsteps'
!     ctrd    -  Stores basic chemical time step length 'ctd'
!     f_initial      -  Stores family concentrations at beginning of call
!     exit_code  -  Error code from the integrator:
!                   0 = successful return
!                   1 = negatives encountered
!                   2 = convergence failure after 'nrsteps' iterations
!                   3 = convergence failure due to divergence - 'NaN's
!                   4 = convergence failure (as '2') but set debugging
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE asad_spmjpdriv_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_SPMJPDRIV_MOD'

CONTAINS

SUBROUTINE asad_spmjpdriv(ix,jy,nlev,n_points,num_iter)

USE asad_mod, ONLY: cdt, f, jpcspf, jpspec, ltrig,                             &
                    ncsteps, nitfg, speci, y
USE ukca_config_specification_mod, ONLY: ukca_config
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper

USE ukca_um_legacy_mod,  ONLY: mype

USE errormessagelength_mod, ONLY: errormessagelength

USE asad_diffun_mod, ONLY: asad_diffun
USE asad_spimpmjp_mod, ONLY: asad_spimpmjp
USE asad_ftoy_mod, ONLY: asad_ftoy

IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: n_points
INTEGER, INTENT(IN) :: ix                    ! i counter
INTEGER, INTENT(IN) :: jy                    ! j counter
INTEGER, INTENT(IN) :: nlev
INTEGER, INTENT(OUT) :: num_iter ! added iteration counter

! Local variables
INTEGER, PARAMETER :: max_redo=128      ! Max times for halving TS, was 16
INTEGER :: exit_code                        ! Convergence exit code
INTEGER :: ncst
INTEGER :: iredo
INTEGER :: jl
INTEGER :: iter
INTEGER :: location
INTEGER :: solver_iter
INTEGER :: jit

INTEGER :: errcode                ! Variable passed to ereport
LOGICAL :: not_first_call = .FALSE.

CHARACTER(LEN=errormessagelength) :: cmessage

REAL :: ctrd
REAL :: f_initial(n_points,jpcspf) ! Saved f array from previous solver call

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_SPMJPDRIV'

!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
ncst = ncsteps
ctrd = cdt
ltrig=.FALSE.
!
CALL asad_diffun( n_points )
!
iredo = 1
f_initial(1:n_points,:)=f(1:n_points,:)
!
IF (ukca_config%l_ukca_asad_columns) THEN
   ! mapping to theta_field
  location = ix + ((jy - 1)*ukca_config%row_length)
ELSE
  location = nlev
END IF
! Start iterations here.
num_iter = 0
iter = 1
DO WHILE (iter <= iredo)
  CALL asad_spimpmjp(exit_code, ix, jy, nlev, n_points, location, solver_iter)
  num_iter = num_iter + solver_iter

  !  Debug slow convergence systems - switch this on in 'spimpmjp'
  IF (exit_code == 4) THEN
    IF (ltrig) THEN
      errcode=1
      cmessage='Slow-converging system, Set printstatus for Jacobian debug'
      DO jl=1,jpspec
        WRITE(umMessage,'(a4,i6,a12,2e14.5,i12)') 'y: ',jl,speci(jl),          &
            MAXVAL(y(:,jl)), MINVAL(y(:,jl)),SIZE(y(:,jl))
        CALL umPrint(umMessage,src='asad_spmjpdriv')
      END DO
      CALL ereport('ASAD_SPMJPDRIV',errcode,cmessage)
    END IF
    ltrig=.TRUE.
    f(1:n_points,:)=f_initial(1:n_points,:)
    ! Call asad_ftoy with jit = 0 to reinitialise y array
    jit = 0
    CALL asad_ftoy( not_first_call, nitfg, jit, n_points, ix, jy, nlev )
    CALL asad_diffun( n_points )
    iter = 1
  ELSE
    !
    !  Reset for failed convergence
    IF (exit_code > 1) THEN
      ncsteps = ncsteps*2
      cdt = cdt/2.0
      iredo = iredo*2
      IF (ukca_config%l_ukca_debug_asad) THEN
        ! Added extra print statements here for verbosity
        WRITE(umMessage,                                                       &
             "('ASAD: failed to converge at location  = ',I0)") location
        CALL umPrint(umMessage,src='asad_spmjpdriv')
        WRITE(umMessage,'(A,I0,A,I0,A,E18.8)')                                 &
              'ASAD: halving timestep: ncsteps = ', ncsteps,                   &
              ' iredo = ', iredo, ' cdt = ', cdt
        CALL umPrint(umMessage,src='asad_spmjpdriv')
      END IF
      IF (cdt < 1.0e-05) THEN
        errcode=2
        cmessage=' Time step now too short'
        CALL ereport('ASAD_SPMJPDRIV',errcode,cmessage)
      END IF
      f(1:n_points,:)=f_initial(1:n_points,:)
      ! Drop out at some point - if 3 successive halvings fail
      IF (iredo >= max_redo) THEN
        IF (printstatus >= prstatus_oper) THEN
          WRITE(umMessage,"(' Resetting array after',i4,' iterations')")       &
              iredo
          CALL umPrint(umMessage,src='asad_spmjpdriv')
          WRITE(umMessage,"('NO CONVERGENCE location: ',i4,' pe: ',i4)")       &
              location,mype
          CALL umPrint(umMessage,src='asad_spmjpdriv')
        END IF
        ncsteps = ncst
        cdt = ctrd
        GO TO 9999
      END IF
      ! Call asad_ftoy with jit = 0 to reinitialise y array
      jit = 0
      CALL asad_ftoy( not_first_call, nitfg, jit, n_points, ix, jy, nlev )
      CALL asad_diffun( n_points )
      iter = 1
    ELSE
      iter = iter + 1
    END IF
  END IF
END DO

IF (iredo > 1) THEN
  IF (iredo > 2) THEN
    WRITE(umMessage,"('   No. iterations =',i2)") iredo
    CALL umPrint(umMessage,src='asad_spmjpdriv')
  END IF
  ncsteps = ncst
  cdt = ctrd
END IF

9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_spmjpdriv


END MODULE asad_spmjpdriv_mod
