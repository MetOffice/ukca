! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Purpose: Main chemistry driver routine
!
!     Chemistry driver routine. If necessary, model concentrations are
!     converted from vmr to number density. If the chemistry is to be
!     "process-split" from the transport, the chemistry tendencies are
!     integrated by the chosen method and an average chemistry tendency
!     over the model timestep is returned. If the chemistry is not to
!     be integrated, instantaneous chemistry tendencies are returned.
!     The tracer tendencies returned to the calling routine are the
!     final values held in the chemistry.
!
!     Note: It is important for conservation that the concentrations
!     passed to this routine are non-negative.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_CHEMISTRY_CTL
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     Arguments:
!        cdot        - Tracer tendencies due to chemistry.
!        ftr         - Tracer concentrations.
!        pp          - Pressure (Nm-2).
!        pt          - Temperature (K).
!        pq          - Water vapor field (vmr).
!        nlev        - Model level
!        dryrt       - Dry deposition rates (s-1)
!        wetrt       - Wet deposition rates (s-1)
!        n_points    - No. of points calculations be done.
!
!     Method
!     ------
!     Since the heterogeneous rates may depend on species
!     concentrations, the call to hetero is made inside the
!     loop over chemistry sub-steps.
!
!     Photolysis rates may be computed at frequency set by the
!     variable nfphot and so the call to photol is also inside
!     the chemical sub-step loop.
!
!     Externals
!     ---------
!     asad_bimol  - Calculates bimolecular rate coefficients.
!     asad_trimol - Calculates trimolecular rate coefficients.
!     ukca_photol - Calculates photolysis rate coefficients.
!     ukca_drydep - Calculates dry deposition rates.
!     ukca_wetdep - Calculates wet deposition rates.
!     asad_totnud - Calculates total number densities.
!     asad_impact - IMPACT time integration scheme.
!     asad_hetero - Calculates heterogeneous rate coefficients.
!     asad_posthet- Performs housekeeping after heterogeneous chemistry.
!     asad_ftoy   - Partitions families.
!     asad_diffun - Calculates chemistry tendencies.
!
!     Local variables
!     ---------------
!     ifam       Family index of in/out species.
!     itr        Tracer index of in/out species.
!     iodd       Number of odd atoms in in/out species.
!     gfirst     .true. when the species need to be
!                initialised on the first chemical step.
!     gphot      .true. if photol needs to be called.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_cdrive_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_CDRIVE_MOD'

CONTAINS

SUBROUTINE asad_cdrive(cdot, ftr, pp, pt, pq, co2_1d, cld_f, cld_l,            &
                       ix, jy, nlev, dryrt, wetrt, rc_het, prt,                &
                       n_points, have_nat, stratflag, H_plus_1d_arr)

USE asad_mod,        ONLY: cdt, ctype, fdot, f,                                &
                           jpspec, jpcspf, jppj,                               &
                           jpdd, jpdw, jpif, jsubs,                            &
                           linfam, lvmr,                                       &
                           madvtr, method, moffam,                             &
                           ncsteps, ndepd, ndepw,                              &
                           nfphot, nit0, nitfg, nodd,                          &
                           p, t, tnd, wp, co2
USE ukca_hetero_mod, ONLY: ukca_hetero, ukca_solidphase
USE ukca_config_specification_mod, ONLY: ukca_config

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport

USE umPrintMgr, ONLY: umMessage, umPrint
USE errormessagelength_mod, ONLY: errormessagelength

USE asad_bedriv_mod, ONLY: asad_bedriv
USE asad_bimol_mod, ONLY: asad_bimol
USE asad_diffun_mod, ONLY: asad_diffun
USE asad_hetero_mod, ONLY: asad_hetero
USE asad_impact_mod, ONLY: asad_impact
USE asad_jac_mod, ONLY: asad_jac
USE asad_spmjpdriv_mod, ONLY: asad_spmjpdriv
USE asad_totnud_mod, ONLY: asad_totnud
USE asad_trimol_mod, ONLY: asad_trimol
USE ukca_drydep_mod, ONLY: ukca_drydep
USE ukca_wetdep_mod, ONLY: ukca_wetdep
USE ukca_photol_mod, ONLY: ukca_photol
USE asad_posthet_mod, ONLY: asad_posthet
USE asad_ftoy_mod, ONLY: asad_ftoy
IMPLICIT NONE


! Subroutine interface
INTEGER, INTENT(IN) :: n_points              ! No of points
INTEGER, INTENT(IN) :: ix                    ! i counter
INTEGER, INTENT(IN) :: jy                    ! j counter
INTEGER, INTENT(IN) :: nlev                  ! Model level

REAL, INTENT(IN) :: prt(n_points,jppj)       ! Photolysis rates
REAL, INTENT(IN) :: dryrt(n_points,jpdd)     ! Dry dep rates
REAL, INTENT(IN) :: wetrt(n_points,jpdw)     ! Wet dep rates
REAL, INTENT(IN) :: rc_het(n_points,2)       ! Hetero. Chemistry rates
REAL, INTENT(IN) :: pp(n_points)             ! Pressure
REAL, INTENT(IN) :: pt(n_points)             ! Temperature
REAL, INTENT(IN) :: pq(n_points)             ! Water vapour
REAL, INTENT(IN) :: co2_1d(n_points)         ! CO2
REAL, INTENT(IN) :: cld_f(n_points)          ! Cloud fraction
REAL, INTENT(IN) :: cld_l(n_points)          ! Cloud liquid water (kg/kg)
LOGICAL, INTENT(IN) :: have_nat(n_points)    ! Mask for natpsc
LOGICAL, INTENT(IN) :: stratflag(n_points)   ! Strat indicator

! 1-D pH array to be passed to asad_hetero
REAL, INTENT(IN) :: H_plus_1d_arr(n_points)

! Size of arrays adapted to number of chemically active species,
! rather than number of tracers
REAL, INTENT(IN OUT) :: ftr(n_points,jpcspf)   ! Tracer concs

REAL, INTENT(OUT)   :: cdot(n_points,jpcspf)  ! Tracer tendencies

!       Local variables

INTEGER :: errcode               ! Variable passed to ereport

INTEGER :: jtr                                ! Loop variable
INTEGER :: jl                                 ! Loop variable
INTEGER :: js                                 ! Loop variable
INTEGER :: nl
INTEGER :: ifam
INTEGER :: itr
INTEGER :: iodd

INTEGER :: num_iter ! To store no.of iterations by chem solver

LOGICAL :: gfirst
LOGICAL :: gphot

LOGICAL :: first_call = .TRUE.

CHARACTER(LEN=errormessagelength) :: cmessage          ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_CDRIVE'


!       1.  Initialise variables and arrays

!       1.1   Clear tendencies to avoid contributions from levels
!             on which no chemistry is performed

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO jtr = 1, jpcspf
  DO jl = 1, n_points
    cdot(jl,jtr) = 0.0
  END DO
END DO

!       1.2  Copy pressure and temperature to asad_mod

DO jl = 1, n_points
  p(jl) = pp(jl)
  t(jl) = pt(jl)
END DO

!       1.2.1 Copy water vapor and co2 to asad_mod

DO jl = 1, n_points
  wp(jl) = pq(jl)
  co2(jl) = co2_1d(jl)
END DO

!       2.  Calculate total number densities

CALL asad_totnud(n_points)

!       3.  Read model tracer concentrations into working array,
!           and if necessary, convert vmr to number densities

IF ( lvmr ) THEN
  DO jtr = 1, jpcspf
    DO jl  = 1, n_points
      ftr(jl,jtr) = ftr(jl,jtr) * tnd(jl)
      f(jl,jtr)   = ftr(jl,jtr)
    END DO
  END DO
ELSE
  DO jtr = 1, jpcspf
    DO jl  = 1, n_points
      f(jl,jtr)   = ftr(jl,jtr)
    END DO
  END DO
END IF

!       4.  Calculate reaction rate coefficients
!           --------- -------- ---- ------------

CALL asad_bimol (n_points)
CALL asad_trimol(n_points)

! Calculate aqueous-phase SO2 oxdn. and tropospheric heterogeneous rates
IF (ukca_config%l_ukca_nr_aqchem .OR. ukca_config%l_ukca_trophet)              &
  CALL asad_hetero(n_points, cld_f, cld_l, rc_het, H_plus_1d_arr)

!       5.  Calculate deposition and emission rates
!           --------- ---------- --- -------- -----

IF ( ndepw /= 0 ) CALL ukca_wetdep(wetrt, n_points)
IF ( ndepd /= 0 ) CALL ukca_drydep(nlev, dryrt, n_points)

!       6.  Integrate chemistry by chosen method. Otherwise,
!           simply calculate tendencies due to chemistry
!           ------ --------- ---------- --- -- ---------

gphot = .TRUE.
IF ( method /= 0 ) THEN

  DO js = 1, ncsteps
    jsubs = js

    gfirst = jsubs  ==  1
    gphot  = gfirst
    IF ( nfphot /= 0 .AND. .NOT. gfirst )                                      &
                  gphot = MOD(jsubs-1,nfphot) == 0

    !           ---------------------------------------------------
    !           NON-STIFF integrators take values in the range 1-9.
    !           ---------------------------------------------------

    !           6.1  IMPACT integration: first compute heterogeneous
    !                and photolysis rates, species and tendencies.
    !                ===============================================

    nl = n_points

    ! pass num_iter to asad_spmjpdriv and asad_ftoy
    num_iter=0

    IF ( method == 1 ) THEN
      IF ( gphot ) CALL ukca_photol(prt,n_points)
      IF (ukca_config%l_ukca_het_psc)                                          &
          CALL ukca_hetero(n_points,have_nat,stratflag)
      CALL asad_ftoy( gfirst, nitfg, num_iter, n_points, ix, jy, nlev )
      CALL asad_diffun( nl )
      CALL asad_jac( n_points )
      CALL asad_impact( n_points, ix, jy, nlev )

      !             6.2.  Quasi-steady state scheme.
      !             ================================

    ELSE IF ( method == 2 ) THEN
      cmessage='QSSA not in UM6.5 build'
      errcode=1
      CALL ereport('ASAD_CDRIVE',errcode,cmessage)

      !              6.3   Sparse Newton-Raphson solver
      !              ==================================

    ELSE IF ( method == 3 ) THEN
      IF ( gphot ) CALL ukca_photol(prt,n_points)
      IF (ukca_config%l_ukca_het_psc)                                          &
          CALL ukca_hetero(n_points,have_nat,stratflag)
      CALL asad_ftoy( gfirst, nitfg, num_iter, n_points, ix, jy, nlev )

      CALL asad_spmjpdriv(ix,jy,nlev,n_points,num_iter)

      IF (ukca_config%l_ukca_debug_asad) THEN
        ! Select which print statement to write depending
        ! on whether asad solver is passed vertical columns
        ! or horizontal slices
        IF (ukca_config%l_ukca_asad_columns) THEN
          WRITE(umMessage,                                                     &
          "('Iterations in spmjpdriv = ',I0,' ix = ',I0,' jy = ',I0)")         &
          num_iter, ix, jy
          CALL umPrint(umMessage,src='asad_cdrive')
        ELSE
          WRITE(umMessage,"('Iterations in spmjpdriv = ',I0,' k = ',I0)")      &
          num_iter, nlev
          CALL umPrint(umMessage,src='asad_cdrive')
        END IF
      END IF

      !              6.5   Backward Euler solver
      !              ===========================

    ELSE IF ( method == 5 ) THEN
      IF ( gphot ) CALL ukca_photol(prt,n_points)
      IF (ukca_config%l_ukca_het_psc)                                          &
           CALL ukca_hetero(n_points,have_nat,stratflag)
      CALL asad_ftoy( gfirst, nitfg, num_iter, n_points, ix, jy, nlev )
      CALL asad_bedriv(ix, jy, n_points, nlev)

      !           -------------------------------------------------
      !           STIFF integrators take values in the range 10-19.
      !           -------------------------------------------------

      !           6.10  NAG BDF stiff integrator.

    ELSE IF ( method == 10 ) THEN
      cmessage='NAG BDF not in UM6.5 build'
      errcode=1
      CALL ereport('ASAD_CDRIVE',errcode,cmessage)

      !             6.11  SVODE ODE stiff integrator from NETLIB.

    ELSE IF ( method == 11 ) THEN
      cmessage='SVODE not in UM6.5 build'
      errcode=1
      CALL ereport('ASAD_CDRIVE',errcode,cmessage)

    END IF

    !           6.12  Do any final work for the heterogeneous
    !                 chemistry before the end of the time loop.
    !                 =========================================

    IF ( ukca_config%l_ukca_het_psc ) CALL ukca_solidphase(n_points)

  END DO        ! End of looping over chemical timesteps

ELSE           ! Method is equal to zero

  !       6.99  Not integrating: just compute tendencies.

  method = 0
  IF ( gphot ) CALL ukca_photol(prt,n_points)
  IF ( ukca_config%l_ukca_het_psc )                                            &
       CALL ukca_hetero(n_points,have_nat,stratflag)
  CALL asad_ftoy( first_call, nit0, num_iter, n_points, ix, jy, nlev )
  CALL asad_diffun( nl )
  IF ( ukca_config%l_ukca_het_psc ) CALL ukca_solidphase( n_points )
  IF ( ukca_config%l_ukca_het_psc ) CALL asad_posthet()

END IF        ! End of IF statement for method


!       7.  Determine concentrations and tendencies to be returned to
!           the model depending on whether or not the chemistry has
!           been integrated  -- -- ------- -- --- --- --------- ---
!           ---- ----------

!       7.1  Obtain model family concentrations by subtracting
!            concentrations of in/out species from ASAD families

DO js = 1, jpspec
  IF ( ctype(js) == jpif ) THEN
    ifam = moffam(js)
    itr  = madvtr(js)
    iodd = nodd(js)
    DO jl = 1, n_points
      IF ( linfam(jl,itr) ) f(jl,ifam) =                                       &
                            f(jl,ifam) - iodd*f(jl,itr)
    END DO
  END IF
END DO

!       7.2  Returned values of concentration and chemical tendency

DO jtr = 1, jpcspf
  IF ( method /= 0 ) THEN
    DO jl = 1, n_points
      cdot(jl,jtr) = ( f(jl,jtr)-ftr(jl,jtr)) / (cdt*ncsteps)
      ftr(jl,jtr)  = f(jl,jtr)
    END DO
  ELSE
    DO jl = 1, n_points
      cdot(jl,jtr) = fdot(jl,jtr)
      ftr(jl,jtr)  = f(jl,jtr)
    END DO
  END IF
END DO


!       8.  If necessary, convert from number densities back to vmr
!           -- ---------- ------- ---- ------ --------- ---- -- ---

IF ( lvmr ) THEN
  DO jtr = 1, jpcspf
    DO jl = 1, n_points
      ftr(jl,jtr)  = ftr(jl,jtr)  / tnd(jl)
      cdot(jl,jtr) = cdot(jl,jtr) / tnd(jl)
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_cdrive
END MODULE asad_cdrive_mod
