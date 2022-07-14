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
! Purpose: Integer function to find reaction number
!   Given a reaction, this function looks for a matching entry in the
!   rate file. The reactants and products given to this routine aren't
!   order dependent. This function will look for the reactants in any
!   order, ie. react1 / react2 or react2 / react1.  And likewise for
!   the products. If a reaction only two products, then specify the
!   third product as a blank string.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_BIMOL etc.
!
!     Interface
!     ---------
!     react1 & react2 - character strings holding the two reactants of
!                       the reaction to look for.
!     products     - character array holding the products of the
!                    reaction to look for.
!     np           - no. of products to search for.
!     reactions    - character array holding the array to search throu
!                    Dimensioned as (nr,ns). nr is the no. of
!                    reactions. ns will be the total no. of reactants
!                    and products (e.g. 2 reactants + 3 products = 5 )
!     index        - this is an integer array dimensioned as (nreacts)
!                    reorders the reactions from the external ratefile
!                    computational efficiency so the order they are st
!                    is not the same order they appear in the ratefile
!                    the appropriate indexing array in order that the
!                    position of this reaction in the ascii ratefile c
!                    returned. e.g. use one of the arrays, nbrkx, ntrk
!
!     NOTE!!  Because of the way the string matching is done in this
!     routine
!     **MUST** pass the strings with trailing spaces. e.g. if the
!     string length allowed for species names is 6, then the argument
!     list must look like
!
!     i = findreaction( 'OH    ', 'CO2   ', 'HOCO2 ', '      ', i
!                       '      ', .
!
!     If this function finds a matching reaction, it will return the
!     index to that reaction in the array reactions. If it doesn't
!     find a match it returns 0.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_findreaction_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_FINDREACTION_MOD'

CONTAINS
INTEGER FUNCTION asad_findreaction( react1,react2,products,                    &
                                    np,reactions,idx_in,nr,ns )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN)    :: nr, idx_in(nr), ns, np

CHARACTER(LEN=10), INTENT(IN) :: react1, react2
CHARACTER(LEN=10), INTENT(IN) :: products(np)
CHARACTER(LEN=10), INTENT(IN) :: reactions(nr,ns)

!     Local variables

INTEGER           :: j, jp, jp2, jr, nreacts
INTEGER           :: nprods, iprods, ipargs

CHARACTER(LEN=10) :: prods(ns-2), inprods(ns-2), empty_ten

LOGICAL           :: found(ns-2), thelast

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_FINDREACTION'


!     -----------------------------------------------------------------
!          Search for reactants first.
!          ------ --- --------- ------
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
asad_findreaction = 0
empty_ten             = '         '

nreacts = SIZE(reactions,1)
nprods  = SIZE(reactions,2) - 2    ! always 2 reactants

!     Get the non-blank products that we're looking for.

ipargs = 0
inprods = empty_ten
DO j = 1, np
  IF ( products(j) /= empty_ten ) THEN
    ipargs = ipargs + 1
    inprods(ipargs) = products(j)
  END IF
END DO

j = 1
DO
  jr = idx_in(j)             ! loop over reactions

  IF ((reactions(j,1)==react1 .AND. reactions(j,2)==react2 ) .OR.              &
     ( reactions(j,1)==react2 .AND. reactions(j,2)==react1 )) THEN

    !         Found a reaction with matching reactants, now check all the
    !         possible cases that the products could match. Need to allow
    !         for rate files which have varying no. of products.

    prods = empty_ten
    found = .FALSE.
    iprods = 0
    DO jp = 1, nprods   ! only copy the non-blank products
      IF ( reactions(j,2+jp) /= empty_ten ) THEN
        iprods = iprods + 1
        prods(iprods) = reactions(j,2+jp)
      END IF
    END DO

    !         If no. of nonblank products doesn't match then
    !         try next reaction.

    IF ( iprods /= ipargs ) THEN
      j = j + 1
      IF ( j > nreacts ) GO TO 9999
      CYCLE
    END IF

    !         Otherwise check the names of the products

    DO jp = 1, iprods
      DO jp2 = 1, iprods
        IF (.NOT. found(jp) .AND. inprods(jp)==prods(jp2)) THEN
          found(jp) = .TRUE.
          prods(jp2) = empty_ten
        END IF
      END DO
    END DO

    !         check to see if we have found all the products.
    !         If we have then exit

    thelast = .TRUE.
    DO jp = 1, iprods
      thelast = thelast .AND. found(jp)
    END DO
    IF ( thelast ) THEN
      asad_findreaction = jr
      GO TO 9999
    END IF
  END IF

  !       next reaction

  j = j + 1
  IF ( j > nreacts ) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END DO

9999 CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION asad_findreaction
END MODULE asad_findreaction_mod
