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
! Purpose: Subroutine to initialise ASAD and fill ldepd and ldepw arrays
!          Adapted from original version written by Olaf Morgenstern.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_MAIN1.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_iniasad_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_INIASAD_MOD'

CONTAINS

SUBROUTINE ukca_iniasad(npoints)

USE ukca_chem_defs_mod,   ONLY: chch_defs
USE asad_mod,             ONLY: ldepd, ldepw, ih_o3, ih_h2o2,                  &
                                ih_hno3, ih_so2, asad_mod_init,                &
                                jpctr, jpspec, jpro2, jpcspf, jpbk,            &
                                jptk, jppj, jphk, jpnr, jpdd, jpdw
USE parkind1,             ONLY: jprb, jpim
USE yomhook,              ONLY: lhook, dr_hook
USE ereport_mod,          ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper

USE errormessagelength_mod, ONLY: errormessagelength
USE asad_cinit_mod, ONLY: asad_cinit
IMPLICIT NONE

INTEGER, INTENT(IN) :: npoints   ! no of spatial points

!       Local variables

INTEGER :: errcode              ! Variable passed to ereport
INTEGER :: k                    ! Loop variable
INTEGER :: iw                   ! Loop variable
INTEGER :: jerr(jpspec)         ! For error analysis

CHARACTER (LEN=errormessagelength) :: cmessage  ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INIASAD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (printstatus >= prstatus_oper) THEN
  WRITE(umMessage,'(A)') 'ASAD initialised using:'
  CALL umPrint(umMessage,src='ukca_iniasad')
  WRITE(umMessage,'(A,I3,A,I3,A,I3)')                                          &
      'jpctr: ',jpctr, ' jpspec: ',jpspec, ' jpnr: ',jpnr
  CALL umPrint(umMessage,src='ukca_iniasad')
  WRITE(umMessage,'(A,I3,A,I3,A,I3)')                                          &
      'jpbk: ', jpbk,  ' jptk: ',  jptk,   ' jppj: ',jppj
  CALL umPrint(umMessage,src='ukca_iniasad')
  WRITE(umMessage,'(A,I3,A,I3,A,I3)')                                          &
      'jphk: ', jphk,  ' jpdd: ',  jpdd,   ' jpdw: ',jpdw
  CALL umPrint(umMessage,src='ukca_iniasad')
  WRITE(umMessage,'(A,I3,A,I3)') 'jpro2: ', jpro2, 'jpcspf: ', jpcspf
  CALL umPrint(umMessage,src='ukca_iniasad')
END IF

!$OMP PARALLEL
CALL asad_mod_init(npoints)
!$OMP END PARALLEL

! Set up dry and wet deposition logicals using module switches
ldepd(:) = .FALSE.
ldepw(:) = .FALSE.
DO k=1,jpspec
  ldepd(k) = (chch_defs(k)%switch1 == 1)
  ldepw(k) = (chch_defs(k)%switch2 == 1)
END DO

! Indentify index of henry_defs array for SO2 etc to use in asad_hetero
iw = 0
ih_o3 = 0
ih_h2o2 = 0
ih_hno3 = 0
ih_so2 = 0
DO k=1,jpspec
  IF (ldepw(k)) iw = iw + 1
  IF (chch_defs(k)%speci == 'O3        ' .AND. ldepw(k)) ih_o3 = iw
  IF (chch_defs(k)%speci == 'H2O2      ' .AND. ldepw(k)) ih_h2o2 = iw
  IF ((chch_defs(k)%speci == 'HONO2     ' .OR.                                 &
         chch_defs(k)%speci == 'HNO3      ') .AND. ldepw(k)) ih_hno3 = iw
  IF (chch_defs(k)%speci == 'SO2       ' .AND. ldepw(k)) ih_so2 = iw
END DO

! Check if module sizes are compatible with namelist values

IF (SIZE(chch_defs) /= jpspec) THEN
  cmessage='size of chch_defs inconsistent with jpspec'
  WRITE(umMessage,'(A,I0,I0)') cmessage, SIZE(chch_defs), jpspec
  CALL umPrint(umMessage,src='ukca_iniasad')
  errcode = 1
  CALL ereport('UKCA_INIASAD',errcode,cmessage)
END IF

jerr(:)=0
WHERE (ldepd) jerr=1
IF (SUM(jerr) /= jpdd) THEN
  cmessage='chch_defs%switch1 values inconsistent with jpdd'
  WRITE(umMessage,'(A,I0,I0)') cmessage, SUM(jerr), jpdd
  CALL umPrint(umMessage,src=RoutineName)
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

jerr(:)=0
WHERE (ldepw) jerr=1
IF (SUM(jerr) /= jpdw) THEN
  cmessage='chch_defs%switch2 values inconsistent with jpdw'
  WRITE(umMessage,'(A,I0,I0)') cmessage, SUM(jerr), jpdw
  CALL umPrint(umMessage,src=RoutineName)
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

CALL asad_cinit()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_iniasad

! ######################################################################

SUBROUTINE ukca_iniasad_spatial_vars(npoints)

USE asad_mod,             ONLY: asad_mod_init_spatial_vars
USE parkind1,             ONLY: jprb, jpim
USE yomhook,              ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: npoints   ! no of spatial points

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INIASAD_SPATIAL_VARS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL
CALL asad_mod_init_spatial_vars(npoints)
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_iniasad_spatial_vars

! ######################################################################

SUBROUTINE ukca_delasad_spatial_vars()

USE asad_mod,             ONLY: asad_mod_dealloc_spatial_vars
USE parkind1,             ONLY: jprb, jpim
USE yomhook,              ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DELASAD_SPATIAL_VARS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL
CALL asad_mod_dealloc_spatial_vars()
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_delasad_spatial_vars

END MODULE ukca_iniasad_mod
