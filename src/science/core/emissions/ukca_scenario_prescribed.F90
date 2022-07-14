! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
! Description:
!  Specify surface boundary conditions for chemical tracers
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  This subroutine specifies surface MMRs for halogenated source gases
!  following the A1 scenario of WMO (2006). It also lumps halogen species
!  to account for total chlorine and bromine, where only a restricted
!  number of species is modelled. Trace gas MMRs are taken from scalar
!  environment field values set via UKCA API calls.
!
!  Called from UKCA_MAIN1.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_scenario_prescribed_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName = 'UKCA_SCENARIO_PRESCRIBED_MOD'

CONTAINS

SUBROUTINE ukca_scenario_prescribed(n_lbc_specs, lbc_specs,                    &
                                    lbc_mmr, perform_lumping)

USE ukca_environment_fields_mod, ONLY: atmospheric_ccl4,                       &
                       atmospheric_cfc113, atmospheric_cfc114,                 &
                       atmospheric_cfc115, atmospheric_cfc11,                  &
                       atmospheric_cfc12, atmospheric_ch2br2,                  &
                       atmospheric_ch4, atmospheric_co2,                       &
                       atmospheric_cos, atmospheric_h1202,                     &
                       atmospheric_h1211, atmospheric_h1301,                   &
                       atmospheric_h2, atmospheric_h2402,                      &
                       atmospheric_hfc125, atmospheric_hfc134a,                &
                       atmospheric_hcfc141b, atmospheric_hcfc142b,             &
                       atmospheric_hcfc22, atmospheric_mebr,                   &
                       atmospheric_meccl3, atmospheric_mecl,                   &
                       atmospheric_n2o, atmospheric_chbr3,                     &
                       atmospheric_n2
USE ukca_constants,    ONLY: c_mebr

USE ukca_scenario_common_mod, ONLY: ich2br2, i_cfcl3, i_cf2cl2, i_cf2clcfcl2,  &
                       i_cf2clcf2cl, i_cf2clcf3, i_ccl4, i_meccl3, i_chf2cl,   &
                       i_mecfcl2, i_mecf2cl, i_cf2clbr, i_cf2br2, i_cf3br,     &
                       i_chbr3, i_cf3chf2, i_ch2fcf3,                          &
                       i_cf2brcf2br, i_mecl, i_mebr, i_ch2br2, i_cos, i_n2o,   &
                       i_ch4, i_co2, i_h2, i_n2, n_full_spc, full_lbc,         &
                       ukca_set_lbc_index, ukca_lump_lbcs, ukca_set_lbc_mmr

USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper
USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN)      :: n_lbc_specs            ! Number of lower BC species
REAL, INTENT(IN OUT)      :: lbc_mmr(n_lbc_specs)  ! Lower BC mass mixing ratios
CHARACTER(LEN=10), INTENT(IN) :: lbc_specs(n_lbc_specs)    ! LBC species
LOGICAL, INTENT(IN)      :: perform_lumping        ! T for lumping of lbc's

LOGICAL, SAVE :: first = .TRUE.

! local variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SCENARIO_PRESCRIBED'


! Body of subroutine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first) THEN
  ! output UKCA environment values from parent model
  IF (PrintStatus >= PrStatus_Oper) THEN
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED CH4      = ',atmospheric_ch4
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED N2O      = ',atmospheric_n2o
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED CO2      = ',atmospheric_co2
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED CFC11    = ',atmospheric_cfc11
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED CFC12    = ',atmospheric_cfc12
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED CFC113   = ',atmospheric_cfc113
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED CFC114   = ',atmospheric_cfc114
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED CFC115   = ',atmospheric_cfc115
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED CCl4     = ',atmospheric_ccl4
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED MeCCl3   = ',atmospheric_meccl3
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED HCFC22   = ',atmospheric_hcfc22
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED HCFC141b = ',atmospheric_hcfc141b
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED HCFC142b = ',atmospheric_hcfc142b
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED H1211    = ',atmospheric_h1211
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED H1202    = ',atmospheric_h1202
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED H1301    = ',atmospheric_h1301
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED H2402    = ',atmospheric_h2402
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED MeBr     = ',atmospheric_mebr
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED MeCl     = ',atmospheric_mecl
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED CH2Br2   = ',atmospheric_ch2br2
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED H2       = ',atmospheric_h2
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(A,E12.3)') 'PRESCRIBED COS      = ',atmospheric_cos
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  CALL ukca_set_lbc_index(n_lbc_specs, lbc_specs)
  first = .FALSE.

END IF

! F11/ CFCl3
full_lbc(i_cfcl3) = atmospheric_cfc11
! F12/ CF2Cl2
full_lbc(i_cf2cl2) = atmospheric_cfc12
! F113/ CF2ClCFCl2
full_lbc(i_cf2clcfcl2) = atmospheric_cfc113
! F114/ CF2ClCF2Cl
full_lbc(i_cf2clcf2cl) = atmospheric_cfc114
! F115/ CF2ClCF3
full_lbc(i_cf2clcf3) = atmospheric_cfc115
! CCl4
full_lbc(i_ccl4) = atmospheric_ccl4
! MeCCl3
full_lbc(i_meccl3) = atmospheric_meccl3
! HCFC-22/ CHF2Cl
full_lbc(i_chf2cl) = atmospheric_hcfc22
! HCFC-141b/ MeCFCl2
full_lbc(i_mecfcl2) = atmospheric_hcfc141b
! HCFC-142b/ MeCF2Cl
full_lbc(i_mecf2cl)= atmospheric_hcfc142b
! H-1211/ CF2ClBr
full_lbc(i_cf2clbr)= atmospheric_h1211
! H-1202/ CF2Br2
full_lbc(i_cf2br2)= atmospheric_h1202
! H-1301/ CF3Br
full_lbc(i_cf3br)= atmospheric_h1301
! H-2402/ CF2BrCF2Br
full_lbc(i_cf2brcf2br)= atmospheric_h2402
! MeCl
full_lbc(i_mecl)= atmospheric_mecl
! MeBr
full_lbc(i_mebr)= atmospheric_mebr
! Add 6 pptv of bromine to MeBr to account for the effect of very
! short-lived bromine gases - !! SHOULD WE STILL DO THIS PRE-INDUSTRIAL?
IF (ich2br2 == 0) THEN
  full_lbc(i_mebr) = full_lbc(i_mebr) + 5.0e-12*c_mebr
END IF

! CH2Br2
full_lbc(i_ch2br2)= atmospheric_ch2br2
! COS
full_lbc(i_cos)= atmospheric_cos
! N2O
full_lbc(i_n2o) = atmospheric_n2o
! CH4
full_lbc(i_ch4) = atmospheric_ch4
! CO2
full_lbc(i_co2) = atmospheric_co2
! H2
full_lbc(i_h2)  = atmospheric_h2

! Species available only for RCP scenario (set to 0 or default value otherwise)
! N2
full_lbc(i_n2)  = atmospheric_n2
! CHBr3
full_lbc(i_chbr3)= atmospheric_chbr3
! HFC125 / CF3CHF2
full_lbc(i_cf3chf2)= atmospheric_hfc125
! HFC134a / CH2FCF3
full_lbc(i_ch2fcf3)= atmospheric_hfc134a

! Perform lumping
IF (perform_lumping) THEN
  CALL ukca_lump_lbcs(n_full_spc, full_lbc)
END IF

! Assign values to LBC array
CALL ukca_set_lbc_mmr(n_lbc_specs, lbc_mmr)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_scenario_prescribed
END MODULE ukca_scenario_prescribed_mod
