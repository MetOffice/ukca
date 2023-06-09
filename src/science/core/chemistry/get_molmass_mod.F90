! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module containing the function GET_MOLMASS to
!   get the molar mass of a given chemical species
!
! Method:
!   Use values in the module ukca_constants to return the molar
!   mass of a given chemical species. The name of the species
!   (species_name) is an input argument that has been passed by
!   another routine where the species is used as a tracer,
!   emission field, lower boundary condition, etc.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE get_molmass_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GET_MOLMASS_MOD'

CONTAINS

FUNCTION get_molmass (species_name)

USE ukca_constants,  ONLY: m_br, m_brcl, m_bro, m_brono2, m_c, m_c2h4,         &
                           m_c2h6, m_c3h6, m_c3h8, m_c4h10, m_ccl4,            &
                           m_cf2cl2, m_cf2clbr, m_cf2clcfcl2, m_cfcl3,         &
                           m_ch2br2, m_ch4, m_chf2cl, m_cl, m_cl2o2,           &
                           m_clo, m_clono2, m_co, m_cs2, m_dms, m_h2,          &
                           m_h2s, m_hbr, m_hcho, m_hcl, m_hobr, m_hocl,        &
                           m_isop, m_me2co, m_me2s, m_mebr, m_meccl3,          &
                           m_mecho, m_mecl, m_meoh, m_monoterp, m_n2o,         &
                           m_nh3, m_no, m_no2, m_oclo, m_ocs, m_oxylene,       &
                           m_so2, m_so4, m_toluene, m_c2h2, m_tbut2ene,        &
                           m_benzene, m_etoh, m_etcho, m_mek, m_hcooh,         &
                           m_meco2h, m_hoch2cho

USE ereport_mod,     ONLY: ereport
USE parkind1,        ONLY: jpim,  jprb     ! DrHook
USE yomhook,         ONLY: lhook, dr_hook  ! DrHook
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Function arguments
CHARACTER (LEN=10), INTENT(IN) :: species_name ! name of a chemical species,
                                               ! emission field, etc
! Function return value
REAL :: get_molmass

! Local variables
INTEGER                        :: ierr
CHARACTER (LEN=errormessagelength)             :: cmessage

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_MOLMASS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierr = 0

! Use values from the module ukca_constants to return the
! molar mass of a given chemical species
SELECT CASE ( species_name )

CASE ('NO2       ', 'NOx       ')
  !  Fields reported in kg of NO2
  get_molmass = m_no2

CASE ('NO        ')
  get_molmass = m_no

CASE ('H2        ', 'no_H2     ')
  get_molmass = m_h2

CASE ('NH3       ')
  get_molmass = m_nh3

CASE ('N2O       ', 'no_N2O    ')
  get_molmass = m_n2o

CASE ('BC_fossil ', 'BC_biofuel', 'BC_biomass',                                &
      'OM_fossil ', 'OM_biofuel', 'OM_biomass')
  get_molmass = m_c

  ! -----------------------------------------
  ! Sulphur containing
CASE ('SO2       ', 'SO2_low   ', 'SO2_high  ', 'SO2_nat   ')
  ! Some of these are emission fields
  get_molmass = m_so2

CASE ('DMS       ')
  get_molmass = m_dms

CASE ('Me2S      ')
  get_molmass = m_me2s

CASE ('SO4       ')
  get_molmass = m_so4

CASE ('COS       ')
  get_molmass = m_ocs

CASE ('H2S       ')
  get_molmass = m_h2s

CASE ('CS2       ')
  get_molmass = m_cs2

  ! -----------------------------------------
  ! CO, CH4 and NMVOCs
CASE ('CO        ')
  get_molmass = m_co

CASE ('CH4       ')
  get_molmass = m_ch4

CASE ('C2H6      ')
  get_molmass = m_c2h6

CASE ('C2H4      ')
  get_molmass = m_c2h4

CASE ('C3H8      ')
  get_molmass = m_c3h8

CASE ('C3H6      ')
  get_molmass = m_c3h6

CASE ('C4H10     ')
  get_molmass = m_c4h10

CASE ('MeOH      ', 'CH3OH     ')
  get_molmass = m_meoh

CASE ('Me2CO     ')
  get_molmass = m_me2co

CASE ('HCHO      ')
  get_molmass = m_hcho

CASE ('MeCHO     ')
  get_molmass = m_mecho

CASE ('C5H8      ')
  get_molmass = m_isop

CASE ('Monoterp  ', 'APINENE   ', 'BPINENE   ')
  get_molmass = m_monoterp

CASE ('TOLUENE   ')
  get_molmass = m_toluene

CASE ('oXYLENE   ')
  get_molmass = m_oxylene

  ! -----------------------------------------
  ! Long-lived halocarbons containing chlorine
CASE ('MeCl      ')
  get_molmass = m_mecl

CASE ('MeCCl3    ')
  get_molmass = m_meccl3

CASE ('CCl4      ')
  get_molmass = m_ccl4

CASE ('CF2Cl2    ')
  get_molmass = m_cf2cl2

CASE ('CFCl3     ')
  get_molmass = m_cfcl3

CASE ('CF2ClCFCl2')
  get_molmass = m_cf2clcfcl2

CASE ('CF2ClBr   ')
  get_molmass = m_cf2clbr

CASE ('CHF2Cl    ')
  get_molmass = m_chf2cl

  ! -----------------------------------------
  ! Reactive and total chlorine
CASE ('Cl        ', 'TOT_Cl    ')
  get_molmass = m_cl

CASE ('ClO       ', 'Clx       ')
  get_molmass = m_clo

CASE ('Cl2O2     ')
  get_molmass = m_cl2o2

CASE ('HCl       ')
  get_molmass = m_hcl

CASE ('HOCl      ')
  get_molmass = m_hocl

CASE ('ClONO2    ')
  get_molmass = m_clono2

CASE ('OClO      ')
  get_molmass = m_oclo

  ! -----------------------------------------
  !   Bromocarbons, reactive bromine and total bromine
CASE ('MeBr      ')
  get_molmass = m_mebr

CASE ('CH2Br2    ')
  get_molmass = m_ch2br2

CASE ('Br        ', 'TOT_Br    ')
  get_molmass = m_br

CASE ('BrO       ', 'Brx       ')
  get_molmass = m_bro

CASE ('HBr       ')
  get_molmass = m_hbr

CASE ('HOBr      ')
  get_molmass = m_hobr

CASE ('BrONO2    ')
  get_molmass = m_brono2

CASE ('BrCl      ')
  get_molmass = m_brcl

  ! -----------------------------------------
  !   Emitted CRI species
CASE ('C2H2      ')
  get_molmass = m_c2h2

CASE ('TBUT2ENE  ')
  get_molmass = m_tbut2ene

CASE ('BENZENE   ')
  get_molmass = m_benzene

CASE ('EtOH      ')
  get_molmass = m_etoh

CASE ('EtCHO     ')
  get_molmass = m_etcho

CASE ('MEK       ')
  get_molmass = m_mek

CASE ('HCOOH     ')
  get_molmass = m_hcooh

CASE ('MeCO2H    ')
  get_molmass = m_meco2h

CASE ('HOCH2CHO  ')
  get_molmass = m_hoch2cho

  ! -----------------------------------------
  !   Others (report warning)
CASE ('AGE       ')
  get_molmass =  1.0
  ierr        = -1

CASE ('AGE OF AIR')
  get_molmass =  1.0
  ierr        = -1

CASE ('PASSIVE O3')
  get_molmass =  1.0
  ierr        = -1

CASE ('XXX       ')
  get_molmass =  1.0   ! unused species
  ierr        = -1

CASE DEFAULT  ! report error
  ierr     = 1
END SELECT

! -------------------------------------------
! Report warning when the mass is set to 1
! and error if the species is not found
IF (ierr /= 0) THEN
  SELECT CASE (ierr)
  CASE (-1)
    cmessage = 'Molecular mass of ' // TRIM(species_name) // ' set to 1'
  CASE (1)
    cmessage = 'Add molecular mass of ' // TRIM(species_name) //               &
             ' to this routine and check if it is present in ukca_constants'
  END SELECT

  CALL ereport ('GET_MOLMASS', ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN

END FUNCTION get_molmass

END MODULE get_molmass_mod
