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
! Purpose: Sets the species concentrations prior to calculating the
!          chemistry tendencies. This includes partitioning families.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE and ASAD_IMPACT
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     Interface
!     ---------
!     Must always be called before diffun.
!
!     Arguments:
!        ofirst     Should be .true. if this is first call to ftoy
!                   during the current model/dynamical timestep.
!        iter       Maximum number of iterations.
!
!     Method
!     ------
!    After initialisation, the family member concentrations and
!    steady state species are iterated until they converge to a
!    solution. Firstly, the ratios for family members and steady
!    state concentrations for the non-family steady state species
!    are computed. Then, the family member concentrations are
!    evaluated using the ratios. A convergence test is performed
!    for all species at all spatial points. If convergence is not
!    achieved for any species at any point, the iterations continue
!    until either convergence is achieved or the maximum number of
!    iterations is reached.
!
!    In/out species are determined to be in or out of the family
!    on the first call to ftoy during a model timestep according
!    to their lifetime compared with some threshold value. If the
!    species is found to be 'in' the family, its concentration
!    is added to the family f concentration before the sum
!    is partitioned amongst the members. The in/out test is
!    performed at all spatial points. This species then stays in
!    the family for the rest of the call to cdrive. The species is
!    subtracted from the family at the end of cdrive.
!
!    If the species self reacts, the steady state concentration
!    is found by solving a quadratic. Otherwise, it is simply
!    found by dividing the production rate by the loss rate.
!    The self reacting term, qa, has previously been calculated
!    in fyself. Note qa is +ve.
!
!    The ratio of a minor family member to the major member is
!    determined by dividing the steady state concentration of the
!    minor member from fystst by the family member concentration:
!                 R1m = Y1/Ym, R2m = Y2/Ym ....
!    The ratio of the major member to the family as a whole is
!    determined from the minor member ratios weighted by the number
!    of odd atoms in each species:
!                 Rmf = 1 / ( Nm + N1*R1m +N2*R2m + .....)
!    The concentrations are then found from:
!                 Ym = Rmf*Z, Y1 = R1m*Ym, Y2 = R2m*Ym, .....
!    where Z is the remaining family concentration after the
!    appropriate in/out species concentrations have been subtracted.
!
!    Although ASAD supports the user of deposition and emissions,
!    we do not include these terms in the calculation of the ratios
!    even though they represent a loss or production from or in to
!    the family. This is because, family members are assumed to be
!    in steady state, generally taken to be a 'chemical' steady
!    state, one that is controlled by the concentrations of other
!    species. If for instance, a strong emission source was present,
!    for NO, this would give a high bogus value for the NO ratio,
!    which would be wrong because by assuming steady state we assume
!    that the extra NO has already reacted to reach a new equilibrium
!    with its family members.
!
!    Externals
!    ---------
!    fyinit         Sets initial species concentrations.
!    fyself         Calculates self-reacting terms.
!    fyfixr         Calculates family member concentrations
!                   using ratios computed on a previous call.
!    fystst         Calculates concentration of a species
!                   in steady state.
!    prls           Calculates production & loss terms.
!
!    Local variables
!    ---------------
!    ifam           Index of family to which species belongs.
!    imaj           Index of major member of family to which
!                   species belongs.
!    itr            Index of model tracer to which species
!                   corresponds.
!    gconv          .true. if convergence has been achieved.
!    gonce          .true. for only the first call to this routine.
!    zthresh        Threshold value of loss rate to determine
!                   whether in/out species are in or out of family.
!    zy             Value of y on previous iteration.
!    ilstmin        List of species which are either steady state,
!                   'SS', family members ('FM' and 'FT') but EXCLUDING
!                   major family species.
!    istmin         No. of entries in ilstmin.
!    ilft           List of species of type 'FT'.
!    ift            No. of entries in ilft.
!    zb, zc, zd     Variables used in the quadratic. Note that
!                   zb is really '-b' in the quadratic.
!    zb             The loss rate less self reacting terms.
!    zc             The production rate.
!    zd             "b^2 - 4*a*c"
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_ftoy_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_FTOY_MOD'

CONTAINS

SUBROUTINE asad_ftoy(ofirst,iter, num_iter, n_points, ix, jy, nlev)

USE asad_mod,            ONLY: f, y, prod, slos, ratio, qa,                    &
                               peps, cdt, pmintnd, nstst,                      &
                               jpfm, jpif, jpna, moffam,                       &
                               majors, ilstmin, ilft, nodd,                    &
                               nlmajmin, madvtr, linfam,                       &
                               nlstst, ctype, ftol,                            &
                               jsro2, nlfro2, jpspec, jpro2
USE ukca_config_specification_mod, ONLY: ukca_config
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport

USE umPrintMgr, ONLY: umPrint

USE errormessagelength_mod, ONLY: errormessagelength

USE asad_fyfixr_mod, ONLY: asad_fyfixr
USE asad_fyinit_mod, ONLY: asad_fyinit
USE asad_prls_mod, ONLY: asad_prls
USE asad_fyself_mod, ONLY: asad_fyself

IMPLICIT NONE



INTEGER, INTENT(IN)    :: n_points  ! No of spatial points
INTEGER, INTENT(IN)    :: ix        ! i counter
INTEGER, INTENT(IN)    :: jy        ! j counter
INTEGER, INTENT(IN)    :: nlev      ! Model level
INTEGER, INTENT(IN OUT) :: iter      ! Max no of iterations
INTEGER, INTENT(IN)    :: num_iter  ! Current iteration number

LOGICAL, INTENT(IN) :: ofirst

!       Local variables

! Variables for when using RO2-permutation chemistry
REAL          :: fro2(n_points)     ! Total RO2 concentration (molecules/cm3)
REAL          :: fro2_old(n_points) ! Total RO2 concentration previous iteration
REAL          :: damp_ro2           ! Dampening factor for RO2 calculation

INTEGER       :: j           ! Loop variable
INTEGER       :: jit         ! Loop variable
INTEGER       :: jl          ! Loop variable
INTEGER       :: js          ! Loop variable
INTEGER       :: ifam
INTEGER       :: imaj
INTEGER       :: nl
INTEGER       :: istart
INTEGER       :: iend
INTEGER       :: iodd
INTEGER       :: itr
INTEGER       :: icode       ! Error code
INTEGER       :: iro2        ! Counter for RO2 species


INTEGER, SAVE :: istmin
INTEGER, SAVE :: ift
!        INTEGER       :: ilstmin(jpspec)      ! Now in ASAD_MOD
!        INTEGER       :: ilft(jpspec)         ! Now in ASAD_MOD

CHARACTER (LEN=errormessagelength) :: cmessage     ! Error message
INTEGER :: errcode          ! Variable passed to ereport
INTEGER :: errcodes(nstst)  ! Array for recording error codes

REAL          :: zthresh
REAL          :: sl
REAL          :: zy(n_points,jpspec)
REAL          :: zb(n_points)
REAL          :: zc(n_points)
REAL          :: zd(n_points)

LOGICAL       :: gconv
LOGICAL, SAVE :: gonce  = .TRUE.
LOGICAL, SAVE :: first_pass = .TRUE.
LOGICAL, SAVE :: gdepem = .FALSE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_FTOY'


!       1.  Initialise.
!           -----------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                          &
                        zhook_in,zhook_handle)
! OMP CRITICAL will only allow one thread through this code at a time,
! while the other threads are held until completion.
!$OMP CRITICAL (asad_ftoy_init)
errcodes(:) = 0
IF ( first_pass ) THEN
  IF ( gonce ) THEN
    gonce  = .FALSE.
    istmin = 0
    ift    = 0
    DO j = 1, jpspec
      ilstmin(j) = 0
      ilft(j) = 0
    END DO

    DO j = 1, nstst
      js = nlstst(j)
      IF ( ctype(js) == jpfm ) THEN
        ifam = moffam(js)
        imaj = majors(ifam)
        IF ( imaj /= js ) THEN
          istmin          = istmin + 1
          ilstmin(istmin) = js
        END IF
      ELSE IF ( ctype(js) == jpif ) THEN
        istmin          = istmin + 1
        ilstmin(istmin) = js
        ift             = ift + 1
        ilft(ift)       = js
      ELSE IF ( ctype(js) == jpna ) THEN
        istmin          = istmin + 1
        ilstmin(istmin) = js
      ELSE
        errcodes(j) = js
      END IF
    END DO

  END IF   ! End of IF (gonce) statement
  first_pass = .FALSE.
END IF     ! End of IF (first_pass) statement
!$OMP END CRITICAL (asad_ftoy_init)

! Perform error check outside of the loop to better suit GPU runs
IF (ANY(errcodes(:) /= 0)) THEN
  errcode = MINVAL(errcodes,mask=(errcodes(:) /= 0))
  CALL umPrint( '**** ASAD ERROR in ASAD_FTOY!! ',src='asad_ftoy')
  CALL umPrint( 'ASAD_FTOY found an unexpected species type',src='asad_ftoy')
  CALL umPrint( 'in the species list nlstst ',src='asad_ftoy')
  CALL ereport('ASAD_FTOY',errcode,'Found unexpected species type')
END IF

!       1.1 Set concentrations/initialise species

CALL asad_fyinit(ofirst,n_points, ix, jy, nlev)

! If computing RO2-permutation reactions, update total RO2 concentration
IF ( ukca_config%l_ukca_ro2_perm ) THEN


  ! Save previous calculation of total RO2
  fro2_old(1:n_points) = y(1:n_points, jsro2)

  ! Calculate sum of all RO2 species on current iteration
  fro2(1:n_points)     = 0.0

  DO j = 1, jpro2

    ! Get index location of each RO2 species and sum
    iro2             = nlfro2(j)
    fro2(1:n_points) = fro2(1:n_points) + f(1:n_points, iro2)

  END DO   ! End iteration over RO2 species

  ! Apply a dampening factor on the first couple of iterations
  ! to prevent problem solution oscillations. Note the values below
  ! have been found to perform well in tests, but there may be room
  ! for further optimisation.
  IF (num_iter <= 0) THEN
    damp_ro2 = 1.0    ! Don't use a dampening factor if initialising RO2
  ELSE IF (num_iter == 1) THEN
    damp_ro2 = 0.5   ! Midway between previous and new value on 1st iteration
  ELSE IF (num_iter == 2) THEN
    damp_ro2 = 0.9   ! 90% towards new value on 2nd iteration
  ELSE
    damp_ro2 = 0.99  ! Thereafter set to (almost) the new value
  END IF

  ! Calculate the new concentration of RO2
  y(1:n_points, jsro2) = fro2_old(1:n_points) +                                &
      damp_ro2*(fro2(1:n_points) - fro2_old(1:n_points))

END IF

!       1.2 If there are no family or steady species then exit

IF ( nstst == 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                        &
                          zhook_out,zhook_handle)
  RETURN
END IF

!       1.3 Initialise local variables and do sanity checks.

nl = n_points

DO j = 1, nstst
  js = nlstst(j)
  DO jl = 1, n_points
    zy(jl,js) = 0.0
  END DO
END DO

IF (ofirst .AND. iter < 5) THEN
  iter = 5
  icode = -1
  cmessage = 'iter too low on first call, resetting to 5'
  CALL ereport('ASAD_FTOY',icode,cmessage)
END IF

zthresh = 2.0 / cdt

!       2.  Calculate self-reacting terms
!           --------- ------------- -----

IF ( ofirst ) CALL asad_fyself(n_points)

!       3.  Calculate family members using previous ratios
!           --------- ------ ------- ----- -------- ------

IF (iter == 0) THEN
  CALL asad_fyfixr(n_points)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                        &
                          zhook_out,zhook_handle)
  RETURN
END IF

!       4.  Iterate family members (including in/out species)
!           ------- ------ ------- --------------------------
!           and non-model steady-state species
!           --- --------- ------------ -------

DO jit = 1, iter

  !         check to see if a species has gone zero with a nonzero
  !         production.

  IF ( jit > 1 ) THEN
    DO j = 1, nstst
      js = nlstst(j)
      DO jl = 1, n_points
        IF (ABS(y(jl,js)) < peps .AND.                                         &
            prod(jl,js) > peps ) THEN
          y(jl,js) = peps
        END IF
      END DO
    END DO
  END IF

  !         4.1 Calculate production and loss terms.
  !             Note that the effect of deposition (wet & dry) and
  !             emissions are not included in the ratios calculation.
  !             Regardless of whether the user has them turned on or not.
  !             See method above.

  CALL asad_prls( nl, istmin, ilstmin, gdepem )

  !         4.2 Initialise ratios

  DO j = 1, nstst
    js = nlstst(j)
    DO jl = 1, n_points
      ratio(jl,js) = 0.0
    END DO
  END DO

  !         4.3  Compute species values for species in steady state
  !              (minor members of family, 'SS' and 'FT' species)
  !              N.B. because of the way asad computes slos, the y
  !              value can go slightly negative during the quadratic
  !              due to loss of precision. We forcibly fix it here.
  !              Also note that we cannot permit y=0.0 during the ftoy
  !              iteration because of the need to get 'sl'.

  DO j = 1, istmin
    js = ilstmin(j)
    DO jl = 1, n_points
      IF ( y(jl,js) < peps ) THEN
        sl = 0.0
      ELSE
        sl = slos(jl,js) / y(jl,js)
      END IF
      IF ( qa(jl,js) > peps ) THEN
        zb(jl) = sl - qa(jl,js) * y(jl,js)
        zc(jl) = prod(jl,js)
        zd(jl) = zb(jl)*zb(jl) + 4.0*qa(jl,js)*zc(jl)
        IF ( zd(jl) > 0.0 ) THEN
          y(jl,js) = (zb(jl) - SQRT(zd(jl)))/(-2.0*qa(jl,js))
          IF ( y(jl,js)  <   0.0 ) y(jl,js) = 10.0*peps
        ELSE
          y(jl,js) = zb(jl) / ( -2.0 * qa(jl,js) )
        END IF
      ELSE IF (qa(jl,js) <= peps .AND. sl > peps ) THEN
        y(jl,js) = prod(jl,js) / sl
      ELSE
        y(jl,js) = 0.0
      END IF
    END DO
  END DO

  !         4.4  Now compute ratios for minor species members.
  !              ** could possibly take out 1/y(imaj) and convert to '*'

  istart = nlmajmin(3)
  iend   = nlmajmin(4)
  DO j = nlmajmin(3), nlmajmin(4)
    js = nlmajmin(j)
    iodd = nodd(js)
    ifam = moffam(js)
    imaj = 0
    IF ( ifam /= 0 ) imaj = majors(ifam)
    DO jl = 1, n_points
      IF ( y(jl,imaj) > peps ) THEN
        ratio(jl,js)   = y(jl,js) / y(jl,imaj)
        ratio(jl,imaj) = ratio(jl,imaj) +                                      &
                         iodd * ratio(jl,js)
      ELSE
        ratio(jl,js) = 0.0
      END IF
    END DO
  END DO

  !         4.5  Now compute ratios from species of type 'FT' if
  !              applicable. If this is the first call, we also set
  !              whether the species is in or out of the family.

  DO j = 1, ift
    js = ilft(j)
    iodd = nodd(js)
    ifam = moffam(js)
    itr  = madvtr(js)
    imaj = majors(ifam)
    IF ( ofirst .AND. jit == 1 ) THEN
      DO jl = 1, n_points
        IF ( y(jl,js) > peps ) THEN
          sl = slos(jl,js) / y(jl,js)
        ELSE
          sl = 0.0
        END IF
        linfam(jl,itr) = sl  >   zthresh
        IF ( linfam(jl,itr) )                                                  &
             f(jl,ifam) = f(jl,ifam) + iodd*f(jl,itr)
      END DO
    END IF

    DO jl = 1, n_points
      IF ( linfam(jl,itr) ) THEN
        IF ( y(jl,imaj)  >   peps ) THEN
          ratio(jl,js)   = y(jl,js) / y(jl,imaj)
          ratio(jl,imaj) = ratio(jl,imaj) +                                    &
                           iodd * ratio(jl,js)
        ELSE
          ratio(jl,js) = 0.0
        END IF
      ELSE
        y(jl,js) = f(jl,itr)
      END IF
    END DO
  END DO

  !         4.6  Finally compute ratio of major species to family and
  !              hence set the minor family members.

  istart = nlmajmin(1)
  iend   = nlmajmin(2)
  DO j = istart, iend
    js = nlmajmin(j)
    iodd = nodd(js)
    ifam = moffam(js)
    DO jl = 1, n_points
      ratio(jl,js) = 1.0 / (iodd + ratio(jl,js))
      y(jl,js)     = f(jl,ifam) * ratio(jl,js)
    END DO
  END DO

  istart = nlmajmin(3)
  iend   = nlmajmin(4)
  DO j = istart, iend
    js = nlmajmin(j)
    ifam = moffam(js)
    imaj = majors(ifam)
    IF ( ctype(js) /= jpif ) THEN
      DO jl = 1, n_points
        y(jl,js) = y(jl,imaj) * ratio(jl,js)
      END DO
    ELSE
      itr = madvtr(js)
      DO jl = 1, n_points
        IF ( linfam(jl,itr) ) y(jl,js) =                                       &
             y(jl,imaj) * ratio(jl,js)
      END DO
    END IF
  END DO

  !         4.7 Test for convergence.

  gconv = .TRUE.
  DO j = 1,nstst
    js = nlstst(j)
    DO jl = 1, n_points
      IF ( ABS(y(jl,js)-zy(jl,js)) >  ftol*y(jl,js)                            &
      .AND. y(jl,js) >  pmintnd(jl) ) gconv=.FALSE.
    END DO
    IF ( .NOT. gconv ) EXIT
  END DO
  IF (gconv) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                      &
                            zhook_out,zhook_handle)
    RETURN
  END IF

  DO j = 1, nstst
    js = nlstst(j)
    DO jl = 1, n_points
      zy(jl,js) = y(jl,js)
    END DO
  END DO

END DO    ! End of iterations


!       9. Convergence achieved or max. iterations reached.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                          &
                        zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_ftoy
END MODULE asad_ftoy_mod
