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
! Purpose: To set up the arrays holding the indicies of reactions,
!     species etc used in the main parts of the code e.g. prls.
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
!     On exit from this routine, the following arrays will have
!     been set: ngrp, nprdx3, nprdx2, nprdx1, nlall, nlstst, nlf,
!     nldepd and nldepw.
!
!     IMPORTANT! This routine MUST be called after inrats since it
!     expects the nspi array to be initialised.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_inix_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_INIX_MOD'

CONTAINS

SUBROUTINE asad_inix

USE asad_mod,        ONLY: ctype, jpif, jpfm, jpfrpx, jpmsp, jpoo,             &
                           jpna, jpsp, jpspb, ldepd, ldepw,                    &
                           ndepd, ndepw,        nf, nfrpx, ngrp,               &
                           nlall, nldepd, nldepw, nldepx,                      &
                           nlf, nlstst, nnfrp,                                 &
                           nprdx1, nprdx2, nprdx3,                             &
                           nspi, nstst, ntabfp,                                &
                           jpspec, jpnr
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umPrint, umMessage
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE



!       Local variables

INTEGER :: j                    ! Loop variable
INTEGER :: jr                   ! Loop variable
INTEGER :: js                   ! Loop variable
INTEGER :: jp                   ! Loop variable
INTEGER :: iprd
INTEGER :: igrp3
INTEGER :: iloss
INTEGER :: ngmax
INTEGER :: iss
INTEGER :: ispec
INTEGER :: jdry
INTEGER :: jwet
INTEGER :: jpnpx3

INTEGER :: igrp(3)
INTEGER :: idepd(jpspec)
INTEGER :: idepw(jpspec)

LOGICAL :: errcodes(jpnr,jpmsp)
CHARACTER (LEN=2)  :: itype
CHARACTER (LEN=errormessagelength) :: cmessage  ! Error message

LOGICAL :: gdry
LOGICAL :: gwet

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_INIX'


!       1.  Initialise the types of specie arrays.
!           ---------- --- ----- -- ------ -------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
jpnpx3 = (jpnr/(3*3))+3*3
nstst = 0
nf = 0
DO js = 1, jpspec
  itype = ctype(js)
  nlall(js) = js
  IF (itype==jpfm .OR. itype==jpif .OR. itype==jpna) THEN
    nstst = nstst + 1
    nlstst(nstst) = js
  END IF
  ! 'OO' type peroxy radicals are treated as chemically active species in ASAD
  IF ( itype==jpfm .OR. itype==jpif .OR. itype==jpsp .OR. itype==jpoo ) THEN
    nf = nf + 1
    ! Track indexes between f and y arrays
    nlf(nf) = js
  END IF
END DO

!       2.  Fill in the 3, 2 & 1 term product index arrays.
!           ---- -- --- -- - - - ---- ------- ----- -------

DO js = 1, jpspec

  !         2.1  For each species, scan the reactions for groups
  !              of 3, 2 & 1 product terms and loss terms.

  iprd  = 0
  igrp3 = 0
  DO jr = 1, jpnr
    IF ( nfrpx(jr) == 0 ) THEN

      !             If reaction does NOT have fractional products.

      DO jp = 3, jpmsp
        IF ( nspi(jr,jp) == js ) THEN
          iprd = iprd + 1
          igrp(iprd) = jr
          IF ( iprd == 3 ) THEN
            iprd = 0
            igrp3 = igrp3 + 1
            DO j = 1, 3
              nprdx3(j,igrp3,js) = igrp(j)
            END DO
          END IF
        END IF
      END DO
    END IF
  END DO
  ngrp(js,3) = igrp3
  IF ( iprd == 2 ) THEN
    ngrp(js,2) = 1
    nprdx2(1,js) = igrp(1)
    nprdx2(2,js) = igrp(2)
  ELSE IF ( iprd == 1 ) THEN
    ngrp(js,1) = 1
    nprdx1(js) = igrp(1)
  END IF

  iloss = 0
  igrp3 = 0
  DO jr = 1, jpnr
    DO jp = 1, 2
      IF ( nspi(jr,jp) == js ) THEN
        iloss = iloss + 1
        igrp(iloss) = jr
        IF ( iloss == 3 ) THEN
          iloss = 0
          igrp3 = igrp3 + 1
          DO j = 1, 3
            nprdx3(j,igrp3,jpspec+js) = igrp(j)
          END DO
        END IF
      END IF
    END DO
  END DO

  ngrp(jpspec+js,3) = igrp3
  IF ( iloss == 2 ) THEN
    ngrp(js+jpspec,2) = 1
    nprdx2(1,jpspec+js) = igrp(1)
    nprdx2(2,jpspec+js) = igrp(2)
  ELSE IF ( iloss == 1 ) THEN
    ngrp(js+jpspec,1) = 1
    nprdx1(jpspec+js) = igrp(1)
  END IF

END DO

!       2.2  Check array size is ok.

ngmax = ngrp(1,3)
DO js = 2, jpspec
  ngmax = MAX( ngmax, ngrp(js,3) )
END DO
IF ( ngmax > jpnpx3 ) THEN
  WRITE(umMessage,'(2A,I0)') '** ASAD ERROR: The parameter jpnpx3 is too',     &
' small. Please change it not less than ',ngmax
  CALL umPrint(umMessage,src='asad_inix')
  cmessage = ' Parameter jpnpx3 is too small'

  CALL ereport('ASAD_INIX',ngmax,cmessage)
END IF

!       2.3  Build table for fractional products

errcodes(:,:) = .FALSE.
nnfrp = 0
DO jr = 1, jpnr
  IF ( nfrpx(jr) /= 0 ) THEN
    DO jp = 3, jpspb
      iss = nspi(jr,jp)
      IF ( iss /= 0 ) THEN
        nnfrp = nnfrp + 1
        IF (nnfrp > jpfrpx) THEN
          errcodes(jr,jp) = .TRUE.
        ELSE
          ntabfp(nnfrp,1) = iss
          ntabfp(nnfrp,2) = jr
          ntabfp(nnfrp,3) = nfrpx(jr) + (jp-3)
        END IF
      END IF
    END DO
  END IF
END DO

! Perform error check outside of the loop to better suit GPU runs
IF (ANY(errcodes(:,:))) THEN
  WRITE(umMessage,'(2A)') '*ASAD ERROR: The parameter jpfrpx is',              &
  ' too small. Please increase it.'
  CALL umPrint(umMessage,src='asad_inix')
  cmessage=' Parameter jpfrpx is too small for no. of fractional products'
  CALL ereport('ASAD_INIX',nnfrp,cmessage)
END IF

!       3.   Set up lists for deposition and emissions.
!            --- -- ----- --- ---------- --- ----------

ndepd = 0
ndepw = 0
DO  js = 1, jpspec
  IF ( ldepd(js) ) THEN
    ndepd = ndepd + 1
    nldepd(ndepd) = js
  END IF
  IF ( ldepw(js) ) THEN
    ndepw = ndepw + 1
    nldepw(ndepw) = js
  END IF
END DO

!       3.1.  Now form the list used in prls.f

DO js = 1, jpspec
  idepd(js) = nldepd(js)
  idepw(js) = nldepw(js)
END DO

!       First look for species with both dry and wet on.

nldepx(1) = 7
nldepx(2) = nldepx(1) - 1

DO ispec = 1,jpspec
  gdry = .FALSE.
  gwet = .FALSE.

  DO j = 1, jpspec
    IF ( idepd(j) == ispec ) THEN
      gdry = .TRUE.
      jdry = j
    END IF
    IF ( idepw(j) == ispec ) THEN
      gwet = .TRUE.
      jwet = j
    END IF
  END DO

  IF ( gdry .AND. gwet ) THEN
    nldepx(2)         = nldepx(2) + 1
    nldepx(nldepx(2)) = ispec
    idepd(jdry)       = 0
    idepw(jwet)       = 0
  END IF
END DO   ! ispec loop

!       Now get the remaining species with dry dep. only

nldepx(3) = nldepx(2) + 1
nldepx(4) = nldepx(3) - 1
DO j = 1, jpspec
  IF ( idepd(j) /= 0 ) THEN
    nldepx(4) = nldepx(4) + 1
    nldepx(nldepx(4)) = idepd(j)
  END IF
END DO

!       Remaining species with wet deposition on.

nldepx(5) = nldepx(4) + 1
nldepx(6) = nldepx(5) - 1
DO j = 1, jpspec
  IF ( idepw(j) /= 0 ) THEN
    nldepx(6) = nldepx(6) + 1
    nldepx(nldepx(6)) = idepw(j)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_inix
END MODULE asad_inix_mod
