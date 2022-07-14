! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
! Description:
!  Contains routines to perform operations that are common to all 'scenarios'
!  i.e. methods of specifying Lower Boundary Conditions (LBCs) in UKCA.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_scenario_common_mod

USE umPrintMgr, ONLY: umPrint, umMessage, PrStatus_Diag, PrintStatus

USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod,       ONLY: ereport
USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

! Integers holding location of possible species that could be passed-in into
! the lbc_mmr array
INTEGER, SAVE, PUBLIC :: icfcl3=0, icf2cl2=0, icf2clcfcl2=0, icf2clcf2cl=0
INTEGER, SAVE, PUBLIC :: icf2clcf3=0, iccl4=0, imeccl3=0, ichf2cl=0
INTEGER, SAVE, PUBLIC :: imecfcl2=0, imecf2cl=0, icf2clbr=0, icf2br2=0
INTEGER, SAVE, PUBLIC :: icf3br=0, icf2brcf2br=0, imecl=0, imebr=0
INTEGER, SAVE, PUBLIC :: ich2br2=0, icos=0, in2o=0, ich4=0, ico2=0, ih2=0
INTEGER, SAVE, PUBLIC :: in2=0, icf3chf2=0, ich2fcf3=0, ichbr3=0

! Integers holding location of species in the full_lbc array used for lumping
! and to populate the final lbc_mmr array. Location of each species is fixed
INTEGER, PARAMETER, PUBLIC :: i_cfcl3=1,      i_cf2cl2=2,      i_cf2clcfcl2=3
INTEGER, PARAMETER, PUBLIC :: i_cf2clcf2cl=4, i_cf2clcf3=5,    i_ccl4=6
INTEGER, PARAMETER, PUBLIC :: i_meccl3=7,     i_chf2cl=8,      i_mecfcl2=9
INTEGER, PARAMETER, PUBLIC :: i_mecf2cl=10,   i_cf2clbr=11,    i_cf2br2=12
INTEGER, PARAMETER, PUBLIC :: i_cf3br=13,     i_cf2brcf2br=14, i_mecl=15
INTEGER, PARAMETER, PUBLIC :: i_mebr=16,      i_ch2br2=17,     i_cos=18
INTEGER, PARAMETER, PUBLIC :: i_n2o=19,       i_ch4=20,        i_co2=21
INTEGER, PARAMETER, PUBLIC :: i_h2=22,        i_n2=23,         i_cf3chf2=24
INTEGER, PARAMETER, PUBLIC :: i_ch2fcf3=25,   i_chbr3=26

! Array holding all lower boundary species (in all scenarios) during processing
INTEGER, PARAMETER, PUBLIC :: n_full_spc = 26
REAL, PUBLIC :: full_lbc(n_full_spc)

! Routines available from this module
PUBLIC :: ukca_set_lbc_index, ukca_lump_lbcs, ukca_set_lbc_mmr

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_SCENARIO_COMMON_MOD'

CONTAINS

SUBROUTINE ukca_set_lbc_index(n_lbc, lbc_spec)

! ----------------------------------------------------------------------
! Description
!  Set up indices for species that are used as lower boundary conditions
!  in this configuration
! ----------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN)  :: n_lbc                    ! Size of LBC-names array
CHARACTER(LEN=*), INTENT(IN) :: lbc_spec(n_lbc) ! Names of Lower BC species

INTEGER :: i
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_LBC_INDEX'

! Body of subroutine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise full_lbc array (assuming this routine is only called once)
DO i = 1, n_full_spc
  full_lbc(i) = 0.0
END DO

DO i= 1, n_lbc
  SELECT CASE (lbc_spec(i))
  CASE ('CFCl3     ')
    icfcl3=i
  CASE ('CF2Cl2    ')
    icf2cl2=i
  CASE ('CF2ClCFCl2')
    icf2clcfcl2=i
  CASE ('CF2ClCF2Cl')
    icf2clcf2cl=i
  CASE ('CF2ClCF3  ')
    icf2clcf3=i
  CASE ('CCl4      ')
    iccl4=i
  CASE ('MeCCl3    ')
    imeccl3=i
  CASE ('CHF2Cl    ')
    ichf2cl=i
  CASE ('MeCF2Cl   ')
    imecf2cl=i
  CASE ('MeCFCl2   ')
    imecfcl2=i
  CASE ('CF2ClBr   ')
    icf2clbr=i
  CASE ('CF2Br2    ')
    icf2br2=i
  CASE ('CF3Br     ')
    icf3br=i
  CASE ('CF2BrCF2Br')
    icf2brcf2br=i
  CASE ('MeCl      ')
    imecl=i
  CASE ('MeBr      ')
    imebr=i
  CASE ('CH2Br2    ')
    ich2br2=i
  CASE ('COS       ')
    icos=i
  CASE ('N2O       ')
    in2o=i
  CASE ('CH4       ')
    ich4=i
  CASE ('CO2       ')
    ico2=i
  CASE ('H2        ')
    ih2=i
  CASE ('N2        ')
    in2=i
  CASE ('CF3CHF2   ')
    icf3chf2=i
  CASE ('CH2FCF3   ')
    ich2fcf3=i
  CASE ('CHBr3     ')
    ichbr3=i
  END SELECT
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_lbc_index


SUBROUTINE ukca_lump_lbcs(n_spc,full_lbc)

! ----------------------------------------------------------------------
! Description
!  If perform_lumping is .TRUE. then the concentrations of several Cl
!  species will be added to either CFC11 or CFC12 (which goes where
!  is loosly judged by the lifetime of the species), and the
!  concetrations of several Br species will be added to MeBr (CH3Br).
!
!  It should be noted that this may cause issues if radiative feedback
!  from UKCA CFC fields is requested.
! ----------------------------------------------------------------------

USE ukca_constants,    ONLY: c_cf2cl2,c_cf2clcfcl2, c_cf2clcf2cl,              &
    c_cf2clcf3, c_cfcl3, c_ccl4, c_meccl3, c_chf2cl, c_mecfcl2,                &
    c_mecf2cl, c_cf2clbr, c_mecl, c_mebr, c_cf2br2, c_cf3br,                   &
    c_cf2brcf2br

IMPLICIT NONE

INTEGER, INTENT(IN)   :: n_spc           ! Size of LBC array
REAL, INTENT(IN OUT)  :: full_lbc(n_spc) ! Array of all LBC species. Larger than
                                         ! no. of spc involved in lumping but
                                         ! position of each species in the array
                                         ! is assumed to be fixed

! Max. number of species (source, target) involved in lumping of LBCs. This is
! expected to be less than the total LBCs. The order is same for all 'scenarios'
! Currently there are 13 species lumped into 3 major ones (13+3 = total 16 CFC)
! The names of the 13 species are defined in array 'cfc_lumped', while the
! 3 'active' species are part of 'lbc_spec' in module: ukca_config_defs_mod
INTEGER, PARAMETER :: n_spc_cfc = 16

! List that specifies which species replaces which one. 1 = CFC11, 2 = CFC12
! 16 = CH3Br
INTEGER, SAVE :: replace_Cl(n_spc_cfc) =                                       &
                                  [0,0,2,2,                                    &
                                    2,1,1,2,                                   &
                                    1,1,1,0,                                   &
                                    0,0,1,0]

INTEGER, SAVE :: replace_Br(n_spc_cfc) =                                       &
                                  [0,0,0,0,                                    &
                                    0,0,0,0,                                   &
                                    0,0,16,16,                                 &
                                    16,16,0,0]

! conversions factors to use when lumping, filled at first call
REAL, SAVE :: convfac_Cl(n_spc_cfc)    ! Chlorine species
REAL, SAVE :: convfac_Br(n_spc_cfc)    ! Bromine

LOGICAL, SAVE :: first = .TRUE.

INTEGER :: i   ! loop variable

! error handling
INTEGER :: ierr
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_LUMP_LBCS'

! Body of subroutine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first) THEN

  ! Check that total number of LBCs is not less than no. of lumped species
  IF ( n_spc < n_spc_cfc ) THEN
    ierr = 1
    WRITE(cmessage,'(2(A,I0))')                                                &
     ' Total no. of LBCs cannot be less than species involved in lumping',     &
      n_spc, ' < ', n_spc_cfc
    CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  END IF

  ! set-up conversion factor arrays
  convfac_Cl = [                                                               &
       0.0, 0.0, 1.5*c_cf2cl2/c_cf2clcfcl2,                                    &
       c_cf2cl2/c_cf2clcf2cl,                                                  &
       0.5*c_cf2cl2/c_cf2clcf3,                                                &
       1.3333*c_cfcl3/c_ccl4, c_cfcl3/c_meccl3,                                &
       0.5*c_cf2cl2/c_chf2cl,                                                  &
       0.6667*c_cfcl3/c_mecfcl2, 0.3333*c_cfcl3/c_mecf2cl,                     &
       0.3333*c_cfcl3/c_cf2clbr, 0.0,                                          &
       0.0,0.0,0.3333*c_cfcl3/c_mecl,0.0 ]

  convfac_Br = [                                                               &
       0.0, 0.0, 0.0, 0.0,                                                     &
       0.0, 0.0, 0.0, 0.0,                                                     &
       0.0, 0.0, c_mebr/c_cf2clbr, 2.0*c_mebr/c_cf2br2,                        &
       c_mebr/c_cf3br,2.0*c_mebr/c_cf2brcf2br,                                 &
       0.0, 0.0 ]

  !  first = .FALSE.
END IF ! first

! Perform lumping
DO i = 1, n_spc_cfc
  IF (replace_Cl(i) > 0) THEN
    full_lbc(replace_Cl(i)) = full_lbc(replace_Cl(i)) +                        &
                              convfac_Cl(i) * full_lbc(i)
  END IF
  IF (replace_Br(i) > 0) THEN
    full_lbc(replace_Br(i)) = full_lbc(replace_Br(i)) +                        &
                              convfac_Br(i) * full_lbc(i)
  END IF
END DO

first = .FALSE.
! diagnostic output of the LBC values at this timestep (if requested)
IF (PrintStatus >= Prstatus_Diag) THEN
  DO i=1,n_full_spc
    WRITE(umMessage,'(A,I0,1X,E18.8)') ' UKCA Full_LBC: ', i, full_lbc(i)
    CALL umPrint(umMessage,src=RoutineName)
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_lump_lbcs

SUBROUTINE ukca_set_lbc_mmr(n_bdy_vals, lbc_mmr)

! ----------------------------------------------------------------------
! Description
!  Populate the array of lower boundary values from the full_lbc array
!  set up by different scenarios.
!  This involves all the species across the current scenarios, which
!  should not be an issue since index for absent species will be 0 and
!  its value in the full_lbc array will also be set to 0.0
! ----------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN) :: n_bdy_vals          ! Size of LBC-MMR array
REAL, INTENT(IN OUT) :: lbc_mmr(n_bdy_vals) ! Final values of Lower BC species

INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_LBC_MMR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise lbc array
DO i = 1, n_bdy_vals
  lbc_mmr(i) = 0.0
END DO

IF (icfcl3      > 0) lbc_mmr( icfcl3    ) = full_lbc(i_cfcl3) ! F11
IF (icf2cl2     > 0) lbc_mmr(icf2cl2    ) = full_lbc(i_cf2cl2) ! F12
IF (icf2clcfcl2 > 0) lbc_mmr(icf2clcfcl2) = full_lbc(i_cf2clcfcl2) ! F113
IF (icf2clcf2cl > 0) lbc_mmr(icf2clcf2cl) = full_lbc(i_cf2clcf2cl) ! F114
IF (icf2clcf3   > 0) lbc_mmr(icf2clcf3  ) = full_lbc(i_cf2clcf3) ! F115
IF (iccl4       > 0) lbc_mmr(iccl4      ) = full_lbc(i_ccl4) ! CCl4
IF (imeccl3     > 0) lbc_mmr(imeccl3    ) = full_lbc(i_meccl3) ! MeCCl3
IF (ichf2cl     > 0) lbc_mmr(ichf2cl    ) = full_lbc(i_chf2cl) ! HFCF-22
IF (imecfcl2    > 0) lbc_mmr(imecfcl2   ) = full_lbc(i_mecfcl2) ! HCFC-141b
IF (imecf2cl    > 0) lbc_mmr(imecf2cl   ) = full_lbc(i_mecf2cl) ! HCFC-142b
IF (icf2clbr    > 0) lbc_mmr(icf2clbr   ) = full_lbc(i_cf2clbr) ! H-1211
IF (icf2br2     > 0) lbc_mmr(icf2br2    ) = full_lbc(i_cf2br2) ! H-1202
IF (icf3br      > 0) lbc_mmr(icf3br     ) = full_lbc(i_cf3br) ! H-1301
IF (icf2brcf2br > 0) lbc_mmr(icf2brcf2br) = full_lbc(i_cf2brcf2br) ! H-2402
IF (imecl       > 0) lbc_mmr(imecl      ) = full_lbc(i_mecl) ! MeCl
IF (imebr       > 0) lbc_mmr(imebr      ) = full_lbc(i_mebr) ! MeBr
IF (ich2br2     > 0) lbc_mmr(ich2br2    ) = full_lbc(i_ch2br2) ! CH2Br2

IF (icos        > 0) lbc_mmr(icos       ) = full_lbc(i_cos) ! COS
IF (in2o        > 0) lbc_mmr(in2o       ) = full_lbc(i_n2o) ! N2O
IF (ich4        > 0) lbc_mmr(ich4       ) = full_lbc(i_ch4) ! CH4
IF (ico2        > 0) lbc_mmr(ico2       ) = full_lbc(i_co2) ! CO2
IF (ih2         > 0) lbc_mmr(ih2        ) = full_lbc(i_h2)  ! H2
IF (in2         > 0) lbc_mmr(in2        ) = full_lbc(i_n2)  ! N2

! Currently only used by RCP scenario
IF (icf3chf2    > 0) lbc_mmr(icf3chf2   ) = full_lbc(i_cf3chf2) ! HFC125
IF (ich2fcf3    > 0) lbc_mmr(ich2fcf3   ) = full_lbc(i_ch2fcf3) ! HFC134a
IF (ichbr3      > 0) lbc_mmr(ichbr3     ) = full_lbc(i_chbr3) ! CHBr3

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_lbc_mmr

END MODULE ukca_scenario_common_mod
