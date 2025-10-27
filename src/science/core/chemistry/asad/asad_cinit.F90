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
! Purpose: To initialize variables used in the chemistry
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_INIASAD
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     Method
!     ------
!     Input arguments are checked and copied to asad_mod. Species
!     and reaction data are read in. Other variables used in the
!     chemistry are initialised.
!
!     -- Calculation of peps
!
!     At several places in the ASAD code, we needed a min. value
!     to use to guard against zero divides. We have to compute
!     this to allow for all possible computer hardwares and precisions
!     that ASAD might be run at.
!
!     Externals
!     ---------
!     inrats      - Reads species and reaction data.
!     inphot      - Initialises photolysis scheme.
!     ukca_inwdep - Initialises wet deposition data
!                   (user supplied routine)
!     ukca_inddep - Initialises dry deposition data
!                   (user supplied routine)
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_cinit_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_CINIT_MOD'

CONTAINS

SUBROUTINE asad_cinit()

USE asad_mod,        ONLY: cdt, ctype, dpd, dpw, ej, emr,                      &
                           f, fdot, frpb, frph, frpj, frpt, frpx,              &
                           jpif, jpfm, jpfrpx, jppjac, kfphot,                 &
                           ldepd, ldepw, linfam, ljacx,                        &
                           method, nbrkx, ncsteps, ndepd, ndepw,               &
                           nfphot, nfrpx, ngrp, nhrkx,                         &
                           nit0, nitfg, nitnr,                                 &
                           njacx1, njacx2, njacx3, njcgrp,                     &
                           nlall, nldepd, nldepw, nldepx,                      &
                           nlf, nlmajmin, nlstst,                              &
                           nmpjac, nmsjac, nmzjac,                             &
                           npdfr, npjac1,                                      &
                           nprdx1, nprdx2, nprdx3, nprkx,                      &
                           nsjac1, nspi,                                       &
                           ntabfp, ntabpd, ntrkx, nzjac1,                      &
                           peps, prk, prod, rk,                                &
                           slos, spb, sph, spj, spt,                           &
                           y, ydot, ztabpd, nrsteps_max,                       &
                           jpspec, jpcspf, jpnr
USE ukca_config_specification_mod, ONLY: ukca_config, int_method_nr,           &
                                         int_method_be_explicit

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint
USE errormessagelength_mod, ONLY: errormessagelength

USE asad_inijac_mod, ONLY: asad_inijac
USE asad_inix_mod, ONLY: asad_inix
USE asad_inrats_mod, ONLY: asad_inrats
USE asad_setsteady_mod, ONLY: asad_setsteady
USE ukca_inddep_mod, ONLY: ukca_inddep
USE ukca_inwdep_mod, ONLY: ukca_inwdep
IMPLICIT NONE

!       Local variables

INTEGER :: j                      ! Loop variable
INTEGER :: jc                     ! Loop variable
INTEGER :: jp                     ! Loop variable
INTEGER :: js                     ! Loop variable
INTEGER :: jx                     ! Loop variable
INTEGER :: jpnpx3                 ! Loop variable
INTEGER :: errcode                ! Variable passed to ereport
INTEGER :: errcodes(jpspec,2)     ! Array for recording error codes
INTEGER, PARAMETER :: nit0_max    = 50   ! max
INTEGER, PARAMETER :: nitfg_max   = 50   ! max
INTEGER, PARAMETER :: nitnr_max   = 50   ! max

REAL :: sfmin

CHARACTER (LEN=10), PARAMETER :: nullx='          '
CHARACTER (LEN=errormessagelength) :: cmessage    ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_CINIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
ljacx = .TRUE.

jpnpx3 = (jpnr/(3*3))+3*3


! 1.  Check the default number of iterations chosen for the solver.
!     These are only relevant to non-Explicit Backward Euler schemes
! ---------------------------------------------------------------

IF ( method /= int_method_be_explicit ) THEN
  IF (ukca_config%nrsteps < 0 .OR. ukca_config%nrsteps > nrsteps_max) THEN
    WRITE(cmessage,'(2(A,1x,I0))') 'NRSTEPS is not in range - 0:',nrsteps_max, &
          '. NRSTEPS=',ukca_config%nrsteps
    errcode=1
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF

  IF (nit0 < 0 .OR. nit0 > nit0_max) THEN
    WRITE(cmessage,'(2(A,1x,I0))') 'NIT0 is not in range - 0:',nit0_max,       &
          '. NIT0=',nit0
    errcode=2
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF

  IF (nitfg < 0 .OR. nitfg > nitfg_max) THEN
    WRITE(cmessage,'(2(A,1x,I0))') 'NITFG is not in range - 0:',nitfg_max,     &
          '. NITFG=',nitfg
    errcode=3
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF

  IF (nitnr < 0 .OR. nitnr > nitnr_max ) THEN
    WRITE(cmessage,'(2(A,1x,I0))') 'NITNR is not in range - 0:',nitnr_max,     &
          '. NITNR=',nitnr
    errcode=4
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
END IF   ! if method not BE-explicit

!       2.1  Set photolysis frequency.

IF (kfphot < 0 .AND. ABS(kfphot) > INT(ukca_config%timestep)) THEN
  WRITE(umMessage,'(A,I0,A,A)') '**CINIT WARNING: VALUE OF KFPHOT ',kfphot,    &
  'EXCEEDS THE MODEL TIMESTEP. ROUTINE PHOTOL WILL ONLY',                      &
  ' BE CALLED ONCE.'
  CALL umPrint(umMessage,src=RoutineName)
  nfphot = 0
ELSE IF ( kfphot > 0 .AND. kfphot > ncsteps ) THEN
  WRITE(umMessage,'(A,I0,A,A)') '**CINIT WARNING: FREQUENCY KFPHOT ',kfphot,   &
   ' EXCEEDS THE TOTAL NUMBER OF CHEMICAL SUBSTEPS. ROUTINE ',                 &
   ' PHOTOL WILL BE CALLED ONCE ONLY.'
  CALL umPrint(umMessage,src=RoutineName)
  nfphot = 0
ELSE IF (kfphot < 0) THEN
  nfphot = INT( ABS(kfphot)/cdt )
ELSE
  nfphot = kfphot
END IF

!       2.2  Compute minimum safe value (see Method above)

sfmin = TINY(1.0)
sfmin = 10.0**(INT(LOG10(sfmin))+1)
peps  = 1.0e19 * sfmin

!       3.  Set fixed vmrs (Now done in UKCA_MAIN1)

! The arrays below are stored on each thread
!$OMP PARALLEL

!       4.  Clear the species arrays

f      = 0.0
fdot   = 0.0
ej     = 0.0
linfam = .FALSE.

y    = 0.0
ydot = 0.0
prod = 0.0
slos = 0.0
dpd  = 0.0
dpw  = 0.0
emr  = 0.0

!       5.   Clear the rates and index arrays.
!            ----- --- ----- --- ----- -------


rk   = 0.0
prk  = 0.0
!$OMP END PARALLEL

nspi = 0

DO js = 1, jpspec
  ngrp(js,1)          = 0
  ngrp(js,2)          = 0
  ngrp(js,3)          = 0
  nprdx2(1,js)        = 0
  nprdx2(2,js)        = 0
  nprdx1(js)          = 0
  ngrp(js+jpspec,1)   = 0
  ngrp(js+jpspec,2)   = 0
  ngrp(js+jpspec,3)   = 0
  nprdx2(1,js+jpspec) = 0
  nprdx2(2,js+jpspec) = 0
  nprdx1(js+jpspec)   = 0
  nlall(js)           = 0
  nlstst(js)          = 0
  nlf(js)             = 0
  nlmajmin(js)        = 0
  nldepd(js)          = 0
  nldepw(js)          = 0
END DO

DO js = 1, jpspec+6
  nldepx(js)          = 0
END DO

DO js = 1, 2*jpspec
  DO jx = 1, jpnpx3
    nprdx3(1,jx,js) = 0
    nprdx3(2,jx,js) = 0
    nprdx3(3,jx,js) = 0
  END DO
END DO

nbrkx = 0
ntrkx = 0
nprkx = 0
nhrkx = 0

njacx3(1,:,:) = 0
njacx3(2,:,:) = 0
njacx3(3,:,:) = 0

njcgrp(:,1) = 0
njcgrp(:,2) = 0
njcgrp(:,3) = 0
njacx2(1,:) = 0
njacx2(2,:) = 0
njacx1(:)   = 0

! Loop through all species in F array
DO jc = 1, jpcspf
  nmpjac(jc) = 0
  DO jp = 1, jppjac
    npjac1(jp,jc) = 0
  END DO
END DO

! Initialise the character arrays
spb(:,:) = nullx
spt(:,:) = nullx
spj(:,:) = nullx
sph(:,:) = nullx

! Initialise the fractional product arrays
frpb(:)  = 0.0
frpt(:)  = 0.0
frpj(:)  = 0.0
frph(:)  = 0.0
frpx(:)  = 0.0
nfrpx    = 0

ntabfp(1:jpfrpx,1) = 0
ntabfp(1:jpfrpx,2) = 0
ntabfp(1:jpfrpx,3) = 0
nmzjac = 0
nmsjac = 0
nzjac1 = 0
nsjac1 = 0
ntabpd = 0
ztabpd = 0.0
npdfr  = 0

!       6.  Read chemistry data
!           ---- --------- ----

CALL asad_inrats

! Check that deposition and emission is not on for constant species
errcodes(:,:) = 0
DO js = 1, jpspec
  IF ( ldepd(js) .AND. ctype(js)(1:1)  ==  'C' ) THEN
    errcodes(js,1) = js
  END IF
  IF ( ldepw(js) .AND. ctype(js)(1:1)  ==  'C' ) THEN
    errcodes(js,2) = js
  END IF
END DO

! Perform error check outside of the loop to better suit GPU runs
IF (ANY(errcodes(:,:) /= 0)) THEN
  IF (ANY(errcodes(:,1) /= 0)) THEN
    cmessage='Dry deposition turned on for constant species'
    errcode = MINVAL(errcodes(:,1),mask=(errcodes(:,1) /= 0))
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  ELSE
    cmessage='Wet deposition turned on for constant species'
    errcode = MINVAL(errcodes(:,2),mask=(errcodes(:,2) /= 0))
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
END IF


IF ( method == int_method_nr) THEN  ! For Newton-Raphson solver only
  CALL asad_setsteady    ! Initialize steady-state species
END IF

IF ( method >= 10 ) THEN
  DO j = 1, jpspec
    IF ( ctype(j)  ==  jpfm .OR. ctype(j)  ==  jpif ) THEN
      WRITE(umMessage,'(A,A,A,I0,A,I0)')                                       &
    '*** ASAD ERROR: You cannot use families ',                                &
    ' with one of the stiff integrators. If method  >=  10 ',                  &
    ' you cannot have species specified as ',jpfm,' or ',jpif
      CALL umPrint(umMessage,src=RoutineName)
      cmessage = 'ASAD ABORTED'
      errcode = j
      CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
    END IF
  END DO
END IF

!       7.  Set up the index arrays.
!           --- -- --- ----- -------

CALL asad_inix
CALL asad_inijac

!       8.  Read deposition and emission data
!           ---- ---------- --- -------- ----

IF ( ndepw /= 0 ) CALL ukca_inwdep()
IF ( ndepd /= 0 ) CALL ukca_inddep()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_cinit
END MODULE asad_cinit_mod
