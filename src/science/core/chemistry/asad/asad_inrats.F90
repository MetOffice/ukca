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
! Purpose: Reads species and reaction data. Combines reactions into one
!          array and reorders them to put single reactant reactions first
!          to improve code in the prls routine.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CINIT
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     Interface
!     ---------
!     Reads the information:
!                  - Chosen chemistry , contains information
!                    on species involved in the chemistry and
!                    the families/tracers to which they belong.
!                  - Data for bimolecular reactions.
!                  - Data for trimolecular reactions.
!                  - Data for photolysis reactions.
!                  - Data for heterogeneous reactions.
!     from ukca_chem1 module
!
!     Method
!     ------
!     The file specifies the species types using 2 letter
!     codes for easier reading.
!
!             ctype         Meaning
!             'FM'          Family member
!             'FT'          Tracer but will be put into a family
!                           if lifetime becomes short.
!             'TR'          Tracer, advected by calling model.
!             'SS'          Steady state species.
!             'CT'          Constant species.
!             'OO'          Peroxy-radical (RO2) species (MeOO, EtOO etc.) seen
!                           as a normal tracer within ASAD. All summed to give
!                           total RO2 concentration.
!                           May or may not be not transported outside of ASAD.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_inrats_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_INRATS_MOD'

INTEGER, SAVE :: iadv               ! Counter for advected species
INTEGER, SAVE :: ispf               ! Counter for chemical species in f array

CONTAINS

SUBROUTINE asad_inrats_set_sp_lists()

! Sets up lists of tracer (advt) and non transported prognostic (nadvt) species
! and switches to indicate whether certain species are in steady state (for the
! N-R solver). Also initialises some other related ASAD module variables.

USE asad_mod,             ONLY: o1d_in_ss, o3p_in_ss, n_in_ss, h_in_ss,        &
                                advt, specf, ctype, family, jpfm, jpif,        &
                                jpna, jpoo, jpsp, nadvt, nltr3, nltrf,         &
                                nnaf, nlfro2, nlnaro2, nodd, ntr3,             &
                                ntrf, nro2, speci, spro2,                      &
                                jpctr, jpspec, jpro2
USE ukca_chem_defs_mod,   ONLY: chch_defs
USE ukca_config_specification_mod, ONLY: ukca_config
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!       Local variables

INTEGER :: errcode               ! Variable passed to ereport

INTEGER :: inadv                 ! Counter for non-advected species
INTEGER :: iro2                  ! Counter for tracer-steady type RO2 species
INTEGER :: jadv                  ! Loop variable
INTEGER :: js                    ! Loop variable
INTEGER :: k                     ! Loop variable

CHARACTER (LEN=errormessagelength) :: cmessage        ! Error message

LOGICAL :: l_fa

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_INRATS_SET_SP_LISTS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Find out which species are in steady state (for N-R solver)

o1d_in_ss = .FALSE.
o3p_in_ss = .FALSE.
n_in_ss   = .FALSE.
h_in_ss   = .FALSE.
DO k=1,jpspec
  IF (chch_defs(k)%speci=='O(1D)     ' .AND.                                   &
      chch_defs(k)%ctype(1:2)==jpna) o1d_in_ss=.TRUE.
  IF (chch_defs(k)%speci=='O(3P)     ' .AND.                                   &
      chch_defs(k)%ctype(1:2)==jpna) o3p_in_ss=.TRUE.
  IF (chch_defs(k)%speci=='N         ' .AND.                                   &
      chch_defs(k)%ctype(1:2)==jpna) n_in_ss=.TRUE.
  IF (chch_defs(k)%speci=='H         ' .AND.                                   &
      chch_defs(k)%ctype(1:2)==jpna) h_in_ss=.TRUE.
END DO

! Currently N-R solver assumes that O(1D) is in steady-state
IF (.NOT. o1d_in_ss .AND. (ukca_config%l_ukca_strat .OR.                       &
    ukca_config%l_ukca_stratcfc .OR. ukca_config%l_ukca_strattrop .OR.         &
    ukca_config%l_ukca_cristrat)) THEN
  cmessage=' O(1D) is not a Steady-State species'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Setup species lists and find which species are in families.

nltrf(:)    = 0
nltr3(:)    = 0

IF (SIZE(chch_defs) /= jpspec) THEN
  errcode=1
  cmessage=' jpspec and chch_defs are inconsistent'
  WRITE(umMessage,'(A)') cmessage
  CALL umPrint(umMessage,src=RoutineName)

  CALL ereport(RoutineName,errcode,cmessage)
END IF
DO k=1,jpspec
  speci(k)  = chch_defs(k)%speci
  nodd(k)   = chch_defs(k)%nodd
  ctype(k)  = chch_defs(k)%ctype(1:2)
  family(k) = chch_defs(k)%family
END DO

iadv = 0
inadv = 0
ntrf = 0
nnaf = 0
ntr3 = 0
! Set counters for number RO2 species
iro2 = 0
nro2 = 0
! Add counter for total number of chemically active species
ispf = 0

DO js = 1, jpspec
  IF ( ctype(js) /= jpfm .AND. ctype(js) /= jpif )                             &
    family(js)='          '
  IF ( ctype(js) == jpsp ) THEN               ! tracers
    iadv = iadv + 1
    ntrf = ntrf + 1
    ispf = ispf + 1
    IF ( iadv > jpctr ) THEN
      WRITE(umMessage,'(A)') '** ASAD ERROR in subroutine ' // RoutineName
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I3)') '** Parameter jpctr is too low; found',iadv
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I3)') '** tracers so far with ',jpspec-js
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A)') '** species to check.'
      CALL umPrint(umMessage,src=RoutineName)
      cmessage = 'ASAD ERROR: jpctr is too low'
      CALL ereport(RoutineName,iadv,cmessage)
    END IF
    advt(iadv)  = speci(js)

    ! Tracers added to lists of all species treated in f array,
    ! which includes non-transported RO2 species if option is set
    nltrf(ntrf) = ispf ! ispf is always == ntrf
    specf(ntrf) = speci(js)

    ! Second block to separate tracer-steady type RO2 species from other tracers
  ELSE IF (ctype(js) == jpoo) THEN
    ! First add to RO2 counters and arrays
    iro2 = iro2 + 1
    nro2 = nro2 + 1
    spro2(iro2) = speci(js)
    IF ( iro2 > jpro2 ) THEN
      WRITE(umMessage,'(A)') '** ASAD ERROR in subroutine ' // RoutineName
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I3)') '** Parameter jpro2 is too low; found',iro2
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I3)') '** tracers so far with ',jpspec-js
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A)') '** species to check.'
      CALL umPrint(umMessage,src=RoutineName)
      cmessage = 'ASAD ERROR: jpro2 is too low'
      CALL ereport(RoutineName,iadv,cmessage)
    END IF

    ! Add RO2 species to the list of species treated by ASAD (ntrf)
    ntrf        = ntrf + 1
    ispf        = ispf + 1
    nltrf(ntrf) = ispf ! ispf is always == ntrf
    specf(ntrf) = speci(js)

    ! Create mapping indices to position in f-array
    nlfro2(iro2) = ispf

    ! See whether or not RO2 species are transported
    IF (ukca_config%l_ukca_ro2_ntp) THEN

      ! Add RO2 species to non-transported counters
      inadv        = inadv + 1
      nnaf         = nnaf + 1
      nadvt(inadv) = speci(js)
      ! Create mapping indices to position in non-advected species list
      nlnaro2(iro2) = inadv
    ELSE

      ! Otherwise add RO2 species to list of transported tracers
      iadv = iadv + 1
      advt(iadv)  = speci(js)
    END IF

  ELSE IF ( ctype(js) == jpna ) THEN       ! Steady-state species
    inadv       = inadv + 1
    nnaf        = nnaf + 1
    nadvt(inadv)= speci(js)

  ELSE IF ( ctype(js) == jpfm .OR. ctype(js) == jpif ) THEN
    cmessage = ' Family chemistry not available in this version'
    errcode = js
    CALL ereport('ASAD_INRATS',errcode,cmessage)
    IF ( ctype(js) == jpif ) THEN
      iadv = iadv + 1
      ntr3 = ntr3 + 1
      IF ( iadv  >   jpctr ) THEN
        WRITE(umMessage,'(A)') '** ASAD ERROR in subroutine ' // RoutineName
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(umMessage,'(A,I3)') '** Parameter jpctr is too low; found',iadv
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(umMessage,'(A,I3)') '** tracers so far with ',jpspec-js
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(umMessage,'(A)') '** species to check.'
        CALL umPrint(umMessage,src=RoutineName)
        cmessage = 'ERROR in jpctr'
        errcode = js
        CALL ereport(RoutineName,errcode,cmessage)
      END IF
      advt(iadv)  = speci(js)
      nltr3(ntr3) = iadv
    END IF
    l_fa=.TRUE.
    jadv_loop: DO jadv = 1, iadv
      IF ( family(js) == advt(jadv) ) THEN
        l_fa=.FALSE.
        EXIT jadv_loop
      END IF
    END DO jadv_loop
    IF (l_fa) THEN
      iadv = iadv + 1
      ntrf = ntrf + 1
      ispf = ispf + 1
      IF ( iadv  >   jpctr ) THEN
        WRITE(umMessage,'(A)') '***** ASAD ERROR in subroutine ' // RoutineName
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(umMessage,'(A,I3)') '***** Param jpctr is too low; found',iadv
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(umMessage,'(A,I3)') '***** tracers so far with ',jpspec-js
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(umMessage,'(A)') '***** species to check.'
        CALL umPrint(umMessage,src=RoutineName)
        cmessage = 'INRATS ERROR : jpctr is too low'
        CALL ereport(RoutineName,iadv,cmessage)
      END IF
      advt(iadv) = family(js)
      nltrf(ntrf) = ispf
    END IF      ! l_fa
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE asad_inrats_set_sp_lists


SUBROUTINE asad_inrats()

USE asad_mod,             ONLY: ab, advt, specf, at, ctype,                    &
                                family, frpb, frph, frpj, frpt, frpx,          &
                                jpfrpb, jpfrph, jpfrpj, jpfrpt, jpif,          &
                                jpspb, jpsph, jpspj, jpspt,                    &
                                madvtr, majors, moffam,                        &
                                nadvt, nbrkx, nfrpx, nhrkx,                    &
                                nlmajmin, nnaf, nuni, nprkx, nspi,             &
                                ntrf, ntrkx,                                   &
                                spb, speci, sph, spj, spt,                     &
                                jpctr, jpspec, jpcspf, jpbk, jptk,             &
                                jppj, jphk, jpnr

USE ukca_chem_defs_mod,   ONLY: ratb_defs, ratt_defs,                          &
                                ratj_defs, rath_defs

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper
USE ukca_um_legacy_mod, ONLY: mype

USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE


!       Local variables

INTEGER :: errcode                ! Variable passed to ereport

INTEGER :: ispb(jpbk+1,jpspb)
INTEGER :: ispt(jptk+1,jpspt)
INTEGER :: ispj(jppj+1,jpspj)
INTEGER :: isph(jphk+1,jpsph)
INTEGER :: ifrpbx(jpbk+1)
INTEGER :: ifrpjx(jppj+1)
INTEGER :: ifrptx(jptk+1)
INTEGER :: ifrphx(jphk+1)
INTEGER :: ilmin(jpspec)
INTEGER :: ilmaj(jpspec)
INTEGER :: imajor                ! Counter
INTEGER :: iminor                ! Counter
INTEGER :: ix                    ! Counter
INTEGER :: icount                ! Counter
INTEGER :: ifam                  ! Index
INTEGER :: imaj                  ! Index
INTEGER :: iflag                 ! Used to test family order
INTEGER :: j                     ! Loop variable
INTEGER :: jf                    ! Loop variable
INTEGER :: jb                    ! Loop variable
INTEGER :: jctr                  ! Loop variable
INTEGER :: jcspf                 ! Loop variable
INTEGER :: jh                    ! Loop variable
INTEGER :: jj                    ! Loop variable
INTEGER :: jp                    ! Loop variable
INTEGER :: jr                    ! Loop variable
INTEGER :: js                    ! Loop variable
INTEGER :: jspb                  ! Loop variable
INTEGER :: jsph                  ! Loop variable
INTEGER :: jspj                  ! Loop variable
INTEGER :: jspt                  ! Loop variable
INTEGER :: jt                    ! Loop variable
INTEGER :: k                     ! Loop variable
INTEGER :: ind                   ! Loop index

CHARACTER (LEN=10), PARAMETER :: nullx='          '
CHARACTER (LEN=errormessagelength) :: cmessage        ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_INRATS'


!       1.  Determine chemistry
!           --------- ---------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise local counters, see asad_cinit for (e.g.) spb and frpb
ispb(:,:) = 0
ifrpbx(:) = 0
ispt(:,:) = 0
ifrptx(:) = 0
ispj(:,:) = 0
ifrpjx(:) = 0
isph(:,:) = 0
ifrphx(:) = 0

!       1.1 Find major species of families

! Find major species of families including arrays with RO2 species
DO jcspf = 1, ntrf
  DO js = 1, jpspec
    IF ( family(js) == specf(jcspf) .AND. ctype(js) /= jpif ) THEN
      majors(jcspf) = js
    END IF

    IF ( speci(js) == specf(jcspf) ) THEN
      majors(jcspf) = js
    END IF
  END DO
END DO

!       1.2 Allocate families to species

! Also including RO2 species in loop
DO js = 1, jpspec
  moffam(js) = 0
  madvtr(js) = 0
  DO jcspf = 1, ntrf
    IF (family(js) == specf(jcspf) ) THEN
      moffam(js) = jcspf
    END IF

    IF (speci(js)  == specf(jcspf) ) THEN
      madvtr(js) = jcspf
    END IF
    IF (family(js) == specf(jcspf) .AND. js > majors(jcspf)) THEN
      WRITE(umMessage,'(A)') '** ASAD ERROR: '
      CALL umPrint(umMessage,src='asad_inrats')
      WRITE(umMessage,'(A)') 'RE-ORDER SPECIES FILE SO THAT THE MAJOR '
      CALL umPrint(umMessage,src='asad_inrats')
      WRITE(umMessage,'(A)') 'SPECIES OF A FAMILY OCCURS AFTER THE OTHERS'
      CALL umPrint(umMessage,src='asad_inrats')
      cmessage = 'INRATS ERROR : Order of species is incorrect'
      errcode = jcspf
      CALL ereport('ASAD_INRATS',errcode,cmessage)
    END IF
  END DO
END DO

!       1.3  Build the list of major and minor species

nlmajmin(1) = 5
imajor      = 0
iminor      = 0
DO js = 1, jpspec
  ifam = moffam(js)
  imaj = 0
  IF ( ifam /= 0 ) THEN
    imaj = majors(ifam)
    IF ( imaj /= js ) THEN
      iminor = iminor + 1
      ilmin(iminor) = js
    ELSE
      imajor = imajor + 1
      ilmaj(imajor) = js
    END IF
  END IF
END DO
nlmajmin(2) = nlmajmin(1) + imajor - 1
nlmajmin(3) = nlmajmin(1) + imajor
nlmajmin(4) = nlmajmin(3) + iminor - 1
DO j = nlmajmin(1), nlmajmin(2)
  nlmajmin(j) = ilmaj(j-nlmajmin(1)+1)
END DO
DO j = nlmajmin(3), nlmajmin(4)
  nlmajmin(j) = ilmin(j-nlmajmin(3)+1)
END DO

IF ( iadv /= jpctr ) THEN
  WRITE(umMessage,'(2A)') '** ASAD ERROR: Number of advected tracers',         &
   ' specified in chch_defs does not match jpctr'
  CALL umPrint(umMessage,src='asad_inrats')
  WRITE(umMessage,'(A,I3,A,I3)') 'Found ',iadv,' but expected ',jpctr
  CALL umPrint(umMessage,src='asad_inrats')
  cmessage = 'INRATS ERROR : iadv and jpctr do not match'

  CALL ereport('ASAD_INRATS',iadv,cmessage)
END IF

! Check number of chemical species within f array
! Should be == iadv if RO2 species are transported,
! should be equal to transported tracers + RO2 species if not.
IF ( ispf /= jpcspf ) THEN
  WRITE(umMessage,'(2A)') '** ASAD ERROR: Number of species treated by asad',  &
   ' specified in chch_defs does not match jpcspf'
  CALL umPrint(umMessage,src='asad_inrats')
  WRITE(umMessage,'(A,I3,A,I3)') 'Found ',ispf,' but expected ',jpcspf
  CALL umPrint(umMessage,src='asad_inrats')
  cmessage = 'INRATS ERROR : ispf and jpcspf do not match'

  CALL ereport('ASAD_INRATS',ispf,cmessage)
END IF


!       2.  Write details of chemistry selection to log file
!           ----- ------- -- --------- --------- -- --- ----

IF (mype == 0 .AND. printstatus >= prstatus_oper) THEN
  CALL umPrint('',src='asad_inrats')
  WRITE(umMessage,'(A)')'  ***  CHEMISTRY INFORMATION  ***'
  CALL umPrint(umMessage,src='asad_inrats')
  CALL umPrint('',src='asad_inrats')
  WRITE(umMessage,'(A)') 'ASAD IS TREATING ADVECTED TRACERS IN THE ORDER:'
  CALL umPrint(umMessage,src='asad_inrats')
  CALL umPrint('',src='asad_inrats')
  DO jctr=1,jpctr
    WRITE(umMessage,'(A,I3,A,A10)')'  ',jctr,' ',advt(jctr)
    CALL umPrint(umMessage,src='asad_inrats')
  END DO
  CALL umPrint('',src='asad_inrats')
  WRITE(umMessage,'(A)')                                                       &
      'ASAD IS TREATING NON-ADVECTED TRACERS IN THE ORDER:'
  CALL umPrint(umMessage,src='asad_inrats')
  CALL umPrint('',src='asad_inrats')
  DO jctr=1,nnaf
    WRITE(umMessage,'(A,I2,A,A10)')'  ',jctr,' ',nadvt(jctr)
    CALL umPrint(umMessage,src='asad_inrats')
  END DO
  CALL umPrint('',src='asad_inrats')

  ! Extra print statement for chemically active species in f array
  WRITE(umMessage,'(A)')                                                       &
      'ASAD IS TREATING CHEMICALLY ACTIVE SPECIES IN THE ORDER:'
  CALL umPrint(umMessage,src='asad_inrats')
  CALL umPrint('',src='asad_inrats')
  DO jcspf=1,ntrf
    WRITE(umMessage,'(A,I2,A,A10)')'  ',jcspf,' ',specf(jcspf)
    CALL umPrint(umMessage,src='asad_inrats')
  END DO
  CALL umPrint('',src='asad_inrats')

  WRITE(umMessage,'(A)')'IF THE TRACERS WERE NOT INITIALISED IN THIS '
  CALL umPrint(umMessage,src='asad_inrats')
  WRITE(umMessage,'(A)')'ORDER THEN THE MODEL RESULTS ARE WORTHLESS '
  CALL umPrint(umMessage,src='asad_inrats')
  CALL umPrint('',src='asad_inrats')

  iflag = 0
  DO jcspf = 1, jpcspf
    IF (specf(jcspf) /= speci(majors(jcspf)) ) THEN
      IF (iflag == 0 ) THEN
        WRITE(umMessage,'(A)')'THE MAJOR MEMBER OF EACH OF THE FAMILIES'
        CALL umPrint(umMessage,src='asad_inrats')
        WRITE(umMessage,'(A)')'IS GIVEN BELOW. IF THIS IS NOT ACCEPTABLE,'
        CALL umPrint(umMessage,src='asad_inrats')
        WRITE(umMessage,'(A)')'THEN YOU MUST REORDER THE SPECIES IN'//         &
                  ' chch_defs'
        CALL umPrint(umMessage,src='asad_inrats')
        WRITE(umMessage,'(A)')'SO THE MAJOR SPECIES FOLLOWS THE OTHERS.'
        CALL umPrint(umMessage,src='asad_inrats')
        WRITE(umMessage,'(A)')
        CALL umPrint(umMessage,src='asad_inrats')
        iflag = 1
      END IF
      WRITE(umMessage,'(A10,1X,A10)') specf(jcspf), speci(majors(jcspf))
      CALL umPrint(umMessage,src='asad_inrats')
    END IF
  END DO
END IF     ! End of IF mype statement

!       3.  Bimolecular ratefile
!           ----------- --------

!       Get bimolecular rates from module

IF (SIZE(ratb_defs) /= jpbk) THEN
  errcode=1
  cmessage='size of ratb_defs is inconsistent with jpbk'

  CALL ereport('ASAD_INRATS',errcode,cmessage)
END IF
icount=1
DO k=1,jpbk
  spb(k,1) = ratb_defs(k)%react1
  spb(k,2) = ratb_defs(k)%react2
  spb(k,3) = ratb_defs(k)%prod1
  spb(k,4) = ratb_defs(k)%prod2
  spb(k,5) = ratb_defs(k)%prod3
  spb(k,6) = ratb_defs(k)%prod4
  ab(k,1)  = ratb_defs(k)%k0
  ab(k,2)  = ratb_defs(k)%alpha
  ab(k,3)  = ratb_defs(k)%beta
  IF (ratb_defs(k)%pyield1 > 1e-18) THEN
    ifrpbx(k)     = icount
    frpb(icount)  = ratb_defs(k)%pyield1
    IF (spb(k,4) /= nullx) frpb(icount+1) = ratb_defs(k)%pyield2
    IF (spb(k,5) /= nullx) frpb(icount+2) = ratb_defs(k)%pyield3
    IF (spb(k,6) /= nullx) frpb(icount+3) = ratb_defs(k)%pyield4
    icount = icount + 4
  END IF
END DO

DO jb = 1, jpbk
  DO js = 1, jpspec
    DO jspb = 1, jpspb
      IF ( speci(js) == spb(jb,jspb) ) ispb(jb,jspb) = js
    END DO
  END DO
END DO

! Load in bimol fractional prod coefs to frpx array
DO jf = 1, jpfrpb
  frpx(jf) = frpb(jf)
END DO

!       4.  Trimolecular ratefile
!           ------------ --------

!       Get trimolecular rates from module for UM version

IF (SIZE(ratt_defs) /= jptk) THEN
  errcode=1
  cmessage='size of ratt_defs is inconsistent with jptk'

  CALL ereport('ASAD_INRATS',errcode,cmessage)
END IF
icount=1
DO k=1,jptk
  spt(k,1) = ratt_defs(k)%react1
  spt(k,2) = ratt_defs(k)%react2
  spt(k,3) = ratt_defs(k)%prod1
  spt(k,4) = ratt_defs(k)%prod2
  at(k,1)  = ratt_defs(k)%f
  at(k,2)  = ratt_defs(k)%k1
  at(k,3)  = ratt_defs(k)%alpha1
  at(k,4)  = ratt_defs(k)%beta1
  at(k,5)  = ratt_defs(k)%k2
  at(k,6)  = ratt_defs(k)%alpha2
  at(k,7)  = ratt_defs(k)%beta2
  IF (ratt_defs(k)%pyield1 > 1e-18) THEN
    ifrptx(k)      = icount
    frpt(icount)   = ratt_defs(k)%pyield1
    IF (spt(k,4) /= nullx) frpt(icount+1) = ratt_defs(k)%pyield2
    icount = icount + 2
  END IF
END DO

DO jt = 1, jptk
  DO js = 1, jpspec
    DO jspt = 1, jpspt
      IF (speci(js) == spt(jt,jspt) ) ispt(jt,jspt) = js
    END DO
  END DO
END DO

! Load in trimol fractional prod coefs to frpx array
DO jf = 1, jpfrpt
  ind = jf + jpfrpb
  frpx(ind) = frpt(jf)
END DO


!       5.  Photolysis ratefile
!           ---------- --------

!       use module to get spj

IF (SIZE(ratj_defs) /= jppj) THEN
  errcode=1
  cmessage='size of ratj_defs is not equal to jppj'

  CALL ereport('ASAD_INRATS',errcode,cmessage)
END IF
icount=1
DO k=1,jppj
  spj(k,1) = ratj_defs(k)%react1
  spj(k,2) = ratj_defs(k)%react2
  spj(k,3) = ratj_defs(k)%prod1
  spj(k,4) = ratj_defs(k)%prod2
  spj(k,5) = ratj_defs(k)%prod3
  spj(k,6) = ratj_defs(k)%prod4
  IF (ratj_defs(k)%pyield1 > 1e-18) THEN
    ifrpjx(k)     = icount
    frpj(icount)  = ratj_defs(k)%pyield1
    IF (spj(k,4) /= nullx) frpj(icount+1) = ratj_defs(k)%pyield2
    IF (spj(k,5) /= nullx) frpj(icount+2) = ratj_defs(k)%pyield3
    IF (spj(k,6) /= nullx) frpj(icount+3) = ratj_defs(k)%pyield4
    icount = icount + 4
  END IF
END DO

DO jj = 1, jppj
  DO js = 1, jpspec
    DO jspj = 1, jpspj
      IF (speci(js) == spj(jj,jspj) ) ispj(jj,jspj) = js
    END DO
  END DO
END DO

! Load in photol fractional prod coefs to frpx array
DO jf = 1, jpfrpj
  ind = jf + jpfrpb + jpfrpt
  frpx(ind) = frpj(jf)
END DO

!       6.  Heterogeneous ratefile
!           ------------- --------

IF (jphk > 0) THEN

  !         use module to get sph

  IF (SIZE(rath_defs) /= jphk) THEN
    errcode=1
    cmessage='size of rath_defs is not equal to jphk'

    CALL ereport('ASAD_INRATS',errcode,cmessage)
  END IF
  icount=1
  DO k=1,jphk
    sph(k,1) = rath_defs(k)%react1
    sph(k,2) = rath_defs(k)%react2
    sph(k,3) = rath_defs(k)%prod1
    sph(k,4) = rath_defs(k)%prod2
    sph(k,5) = rath_defs(k)%prod3
    sph(k,6) = rath_defs(k)%prod4
    IF (rath_defs(k)%pyield1 > 1e-18) THEN
      ifrphx(k)     = icount
      frph(icount)  = rath_defs(k)%pyield1
      IF (sph(k,4) /= nullx) frph(icount+1)=rath_defs(k)%pyield2
      IF (sph(k,5) /= nullx) frph(icount+2)=rath_defs(k)%pyield3
      IF (sph(k,6) /= nullx) frph(icount+3)=rath_defs(k)%pyield4
      icount = icount + 4
    END IF
  END DO

  DO jh = 1, jphk
    DO js = 1, jpspec
      DO jsph = 1, jpsph
        IF (speci(js) == sph(jh,jsph) ) isph(jh,jsph) = js
      END DO
    END DO
  END DO

  ! Load in het fractional prod coefs to frpx array
  DO jf = 1, jpfrph
    ind = jf + jpfrpb + jpfrpt + jpfrpj
    frpx(ind) = frph(jf)
  END DO

END IF       ! jphk > 0

!       7.  Reorder reactions, putting single reactants first.
!           ------- ---------- ------- ------ --------- ------

nuni = 0

!       7.1  Single reactants; scan ratefiles in turn.

!
DO jr = 1, jpbk
  IF (ispb(jr,2) == 0 ) THEN
    nuni = nuni + 1
    nbrkx(jr) = nuni
    DO jp = 1, jpspb
      nspi(nuni,jp) = ispb(jr,jp)
    END DO
    IF ( ifrpbx(jr) /= 0 ) nfrpx(nuni) = ifrpbx(jr)
  END IF
END DO

DO jr = 1, jptk
  IF ( ispt(jr,2) == 0 ) THEN
    nuni = nuni + 1
    ntrkx(jr) = nuni
    DO jp = 1, jpspt
      nspi(nuni,jp) = ispt(jr,jp)
    END DO
    IF (ifrptx(jr) /= 0) nfrpx(nuni) = ifrptx(jr)+jpfrpb
  END IF
END DO

DO jr = 1, jppj
  IF ( ispj(jr,2) == 0 ) THEN
    nuni = nuni + 1
    nprkx(jr) = nuni
    DO jp = 1, jpspj
      nspi(nuni,jp) = ispj(jr,jp)
    END DO
    IF (ifrpjx(jr) /= 0) nfrpx(nuni) = ifrpjx(jr)+jpfrpb+jpfrpt
  END IF
END DO

IF ( jphk > 0 ) THEN
  DO jr = 1, jphk
    IF ( isph(jr,2) == 0 ) THEN
      nuni = nuni + 1
      nhrkx(jr) = nuni
      DO jp = 1, jpsph
        nspi(nuni,jp) = isph(jr,jp)
      END DO
      IF (ifrphx(jr) /= 0) nfrpx(nuni) = ifrphx(jr)+jpfrpb+                    &
                                         jpfrpt+jpfrpj
    END IF
  END DO
END IF

!       7.2  Two reactants; copy remaining reactions

ix = nuni
DO jr = 1, jpbk
  IF ( ispb(jr,2) /= 0 ) THEN
    ix = ix + 1
    nbrkx(jr) = ix
    DO jp = 1, jpspb
      nspi(ix,jp) = ispb(jr,jp)
    END DO

    IF ( ifrpbx(jr) /= 0 ) nfrpx(ix) = ifrpbx(jr)
  END IF
END DO

DO jr = 1, jptk
  IF ( ispt(jr,2) /= 0 ) THEN
    ix = ix + 1
    ntrkx(jr) = ix
    DO jp = 1, jpspt
      nspi(ix,jp) = ispt(jr,jp)
    END DO
    IF (ifrptx(jr) /= 0) nfrpx(ix) = ifrptx(jr)+jpfrpb
  END IF
END DO

DO jr = 1, jppj
  IF (ispj(jr,2) /= 0 ) THEN
    ix = ix + 1
    nprkx(jr) = ix
    DO jp = 1, jpspj
      nspi(ix,jp) = ispj(jr,jp)
    END DO
    IF (ifrpjx(jr) /= 0) nfrpx(ix) = ifrpjx(jr)+jpfrpb+jpfrpt
  END IF
END DO

IF ( jphk > 0 ) THEN
  DO jr = 1, jphk
    IF ( isph(jr,2) /= 0 ) THEN
      ix = ix + 1
      nhrkx(jr) = ix
      DO jp = 1, jpsph
        nspi(ix,jp) = isph(jr,jp)
      END DO
      IF (ifrphx(jr) /= 0) nfrpx(ix) = ifrphx(jr)+jpfrpb+                      &
                                       jpfrpt+jpfrpj
    END IF
  END DO
END IF

IF ( ix /= jpnr ) THEN
  WRITE(umMessage,'(2A)') '*** INTERNAL ASAD ERROR: Number of reactions',      &
              ' placed in nspi array does not equal jpnr. '
  CALL umPrint(umMessage,src='asad_inrats')
  WRITE(umMessage,'(2A)') '                         Check that reaction',      &
              ' files and value of jpnr in UKCA namelist are'
  CALL umPrint(umMessage,src='asad_inrats')
  WRITE(umMessage,'(2A)') '                         consistent. Found: ',      &
              ix,' jpnr: ',jpnr
  CALL umPrint(umMessage,src='asad_inrats')
  cmessage = 'No of rxns in nspi array is not equal to jpnr'

  CALL ereport('ASAD_INRATS',ix,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_inrats
END MODULE asad_inrats_mod
