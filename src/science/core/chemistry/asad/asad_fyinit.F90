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
! Purpose: sets or initialises species concentrations
!
!     Sets the individual species concentrations y from those of
!     the model tracers/families f. If the species concentrations
!     are later to be iterated in ftoy, then they are merely
!     initialised in this routine.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_FTOY
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     Method
!     ------
!     If this is the first call to ftoy for the chemical time loop
!     ie. we're starting a new model/dynamical timestep, then we need
!     to initialise the species array, y. For subsequent calls, all we
!     need to do is copy the species of type 'TR' and 'FT' from the f
!     (family/tracer) array to their places in the species array, y.
!
!     When ofirst is .true., family are initialised by setting the
!     major family member to the family concentration and the minor
!     family members to a small number. Steady state species are also
!     set to a small number. Constant species (species type 'CT')
!     are also only set when .ofirst. is true, since these species are
!     never integrated (see diffun.f). Constant field species (type
!     'CF') are allowed to change every chemical time step however.
!
!     ofirst should be .true. if this is first call to FTOY
!     for the chemistry time loop ie. ofirst should be
!     true only when starting a new dynamical timestep.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_fyinit_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_FYINIT_MOD'

CONTAINS

SUBROUTINE asad_fyinit(ofirst, n_points, ix, jy, nlev)

USE asad_mod,             ONLY: ctype, f,                                      &
                                ilcf, ilct, ilftr, ilss, jpoo,                 &
                                jpcf, jpco, jpif, jpna, jpsp,                  &
                                madvtr, moffam, nlmajmin, peps,                &
                                speci, tnd, y, jsro2, jpspec
USE ukca_diurnal_oxidant, ONLY: ukca_set_diurnal_ox, ukca_set_diurnal_ox_col
USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_environment_fields_mod, ONLY: atmospheric_co2,                        &
                                       atmospheric_h2,                         &
                                       atmospheric_n2,                         &
                                       atmospheric_o2,                         &
                                       atmospheric_ch4
USE ukca_constants, ONLY: c_co2, c_h2, c_n2, c_o2, c_ch4
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE asad_inicnt_mod, ONLY: asad_inicnt
USE asad_inicnt_col_mod, ONLY: asad_inicnt_col
IMPLICIT NONE

INTEGER, INTENT(IN) :: n_points   ! No of spatial points
INTEGER, INTENT(IN) :: ix         ! i counter
INTEGER, INTENT(IN) :: jy         ! j counter
INTEGER, INTENT(IN) :: nlev       ! model level

LOGICAL, INTENT(IN) :: ofirst     ! True on first call

!       Local variables

INTEGER, SAVE :: iss
INTEGER, SAVE :: ict
INTEGER, SAVE :: iftr
INTEGER, SAVE :: icf

INTEGER       :: j           ! Loop variable
INTEGER       :: jl          ! Loop variable
INTEGER       :: js          ! Index
INTEGER       :: istart
INTEGER       :: iend
INTEGER       :: ifam
INTEGER       :: itr

INTEGER       :: errcode    ! variable passed to ereport

REAL          :: peps10
REAL          :: fch4,fco2,fh2,fn2,fo2
                            ! volume mixing ratios of fixed species (mol mol-1)
REAL          :: zcnst

! Default volume mixing ratios (mol mol-1) for fixed species
! ('present day' values)
REAL, PARAMETER :: fco2_default = 350.0e-6
REAL, PARAMETER :: fh2_default  = 5.0e-7
REAL, PARAMETER :: fn2_default  = 0.78084
REAL, PARAMETER :: fo2_default  = 0.20945
REAL, PARAMETER :: fch4_default = 1.76e-6

LOGICAL, SAVE :: gonce = .TRUE.
LOGICAL, SAVE :: first_pass = .TRUE.
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_FYINIT'


!       1.  Species initialised only at start of dynamical step.
!           ------- ----------- ---- -- ----- -- --------- -----

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! OMP CRITICAL will only allow one thread through this code at a time,
! while the other threads are held until completion.
!$OMP CRITICAL (asad_fyinit_init)
IF ( first_pass ) THEN
  IF ( gonce ) THEN

    !         Build lists for this routine. This lot should really be a
    !         common variable - next version!

    gonce = .FALSE.
    iss = 0
    ict = 0
    icf = 0
    iftr = 0
    DO js = 1, jpspec
      IF ( ctype(js) == jpna ) THEN
        iss       = iss + 1
        ilss(iss) = js
      END IF
      IF ( ctype(js) == jpco ) THEN
        ict       = ict + 1
        ilct(ict) = js
      END IF
      IF ( ctype(js) == jpcf ) THEN
        icf       = icf + 1
        ilcf(icf) = js
        ! Save index of RO2
        IF ( speci(js) == 'RO2       ') jsro2 = js
      END IF
      ! Family members, normal tracers and RO2-type species mapped
      ! between f and y arrays
      IF (ctype(js) == jpif .OR.                                               &
          ctype(js) == jpsp .OR.                                               &
          ctype(js) == jpoo ) THEN
        iftr        = iftr + 1
        ilftr(iftr) = js
      END IF
    END DO

  END IF
  first_pass = .FALSE.
END IF           ! End of IF (gonce) statement
!$OMP END CRITICAL (asad_fyinit_init)

IF ( ofirst ) THEN
  peps10 = 10.0 * peps

  !         1.1 Ordinary family members (major then minor )

  istart = nlmajmin(1)
  iend   = nlmajmin(2)
  DO j = istart, iend
    js = nlmajmin(j)
    ifam = moffam(js)
    DO jl = 1, n_points
      y(jl,js) = f(jl,ifam)
    END DO
  END DO

  istart = nlmajmin(3)
  iend   = nlmajmin(4)
  DO j = istart, iend
    js = nlmajmin(j)
    DO jl = 1, n_points
      y(jl,js) = peps10
    END DO
  END DO

  !         1.2 Non-model species in steady state

  DO j = 1, iss
    js = ilss(j)
    DO jl = 1, n_points
      y(jl,js) = peps10
    END DO
  END DO

  !         1.3 Species set to constants (unaffected by remainder of ftoy)

  ! ----------------------------------------------------------------------
  ! Set volume mixing ratios for fixed species:
  ! Base these on external mixing ratios if requested, otherwise use
  ! default values.
  ! ----------------------------------------------------------------------

  fco2 = fco2_default
  fch4 = fch4_default
  fh2 = fh2_default
  fo2 = fo2_default
  fn2 = fn2_default

  ! Replace defaults with environmental values if prescribed
  ! (only available for species in the selected chemistry scheme)

  IF (ukca_config%l_chem_environ_co2_scalar) THEN
    fco2 = atmospheric_co2/c_co2
    IF (ukca_config%l_chem_environ_co2_fld)                                    &
      ! Set to zero as fco2 should not be used if prescribing a full field
      fco2 = 0.0
  END IF

  IF (ukca_config%l_chem_environ_h2_scalar)  fh2 = atmospheric_h2/c_h2
  IF (ukca_config%l_chem_environ_n2_scalar)  fn2 = atmospheric_n2/c_n2
  IF (ukca_config%l_chem_environ_o2_scalar)  fo2 = atmospheric_o2/c_o2
  IF (ukca_config%l_chem_environ_ch4_scalar) fch4 = atmospheric_ch4/c_ch4

  DO j = 1, ict
    js = ilct(j)

    IF ( speci(js) == 'N2        ' ) THEN
      zcnst = fn2
    ELSE IF ( speci(js) == 'O2        ' ) THEN
      zcnst = fo2
    ELSE IF ( speci(js) == 'CO2       ' ) THEN
      zcnst = fco2
    ELSE IF ( speci(js) == 'H2        ' ) THEN
      zcnst = fh2
    ELSE IF ( speci(js) == 'CH4       ' ) THEN
      zcnst = fch4
    ELSE
      cmessage=' Value not supplied for '//speci(js)//                         &
               ',setting value to zero'
      errcode=-1

      CALL ereport('ASAD_FYINIT',errcode,cmessage)
      zcnst = 0.0
    END IF

    DO jl = 1, n_points
      y(jl,js) = zcnst * tnd(jl)
    END DO
  END DO

END IF       ! ofirst

!       2.   Initialise rest of species (potentially each call).
!            ---------- ---- -- ------- ------------ ---- ------

!       2.1  Species of type 'TR' and 'FT' are copied from tracers
!            on every call to this routine.
!            Species of type 'OO' also copied here

DO j = 1, iftr
  js = ilftr(j)
  itr = madvtr(js)
  DO jl = 1, n_points
    y(jl,js) = f(jl,itr)
  END DO
END DO

!       2.2  Constant but set by user, type 'CF' are
!            initialised at every chemical step.
!            Normally only for H2O, but the offline oxidants code
!            initialises O3, OH, HO2, NO3, and H2O2.
!            With RO2-permutation reactions, RO2 is also a CF
IF (ofirst) THEN
  DO j = 1, icf
    js = ilcf(j)
    IF (ukca_config%l_ukca_asad_columns) THEN
      CALL asad_inicnt_col(speci(js), y(1,js), n_points, ix, jy)
    ELSE
      CALL asad_inicnt(speci(js), y(1,js), n_points, nlev)
    END IF
  END DO
END IF

! Impose a diurnal profile on offline oxidants read in as time mean
IF (ofirst .AND.                                                               &
    (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be)) THEN
  DO j = 1, icf
    js = ilcf(j)
    IF (ukca_config%l_ukca_asad_columns) THEN
      CALL ukca_set_diurnal_ox_col(speci(js),  y(:,js), n_points, ix, jy)
    ELSE
      CALL ukca_set_diurnal_ox(speci(js),  y(:,js), n_points, nlev)
    END IF
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_fyinit
END MODULE asad_fyinit_mod
