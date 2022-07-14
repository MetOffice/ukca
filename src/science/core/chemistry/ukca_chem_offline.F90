! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module providing UKCA with inputs describing the chemistry for the
!    off line oxidants version of the sulphur and terpene chemistry for the
!    ASAD system. The subroutine UKCA_INIT_OFFLINE constructs the final
!    arrays from the components.
!
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and
!  The Met Office. See www.ukca.ac.uk
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
!  Contained subroutines:
!   ukca_init_offline
!
! ----------------------------------------------------------------------

MODULE ukca_chem_offline

USE ukca_chem_defs_mod, ONLY: chch_t, ratb_t, rath_t, ratj_t, ratt_t

IMPLICIT NONE
PRIVATE

! Chemistry for sulphur and monoterpene species using offline oxidants
! =============================================================================

TYPE(chch_t), ALLOCATABLE, PUBLIC :: chch_defs_offline(:)
TYPE(ratb_t), ALLOCATABLE, PUBLIC :: ratb_defs_offline(:)
TYPE(ratt_t), ALLOCATABLE, PUBLIC :: ratt_defs_offline(:)
TYPE(ratj_t), ALLOCATABLE, PUBLIC :: ratj_defs_offline(:)
TYPE(rath_t), ALLOCATABLE, PUBLIC :: rath_defs_offline(:)
REAL,         ALLOCATABLE, PUBLIC :: depvel_defs_offline(:,:,:)
REAL,         ALLOCATABLE, PUBLIC :: henry_defs_offline(:,:)

! No of heterogeneous processes
INTEGER, PARAMETER :: naq_offline   = 3        ! Aqueous phase

! No of dry deposited tracers
INTEGER, PARAMETER, PUBLIC :: ndry_offline  = 6

! No of wet deposited species
INTEGER, PARAMETER, PUBLIC :: nwet_offline  = 4        ! tracers

! No of constant species in aqueous phase
INTEGER, PARAMETER, PUBLIC :: nwet_constant = 1        ! O3 field

INTEGER :: n_tracers

! Diagnostic fields (mass mixing ratio)
REAL, PUBLIC, ALLOCATABLE, SAVE :: o3_offline_diag(:,:)   ! O3 diagnostic
REAL, PUBLIC, ALLOCATABLE, SAVE :: oh_offline_diag(:,:)   ! OH diagnostic
REAL, PUBLIC, ALLOCATABLE, SAVE :: no3_offline_diag(:,:)  ! NO3 diagnostic
REAL, PUBLIC, ALLOCATABLE, SAVE :: ho2_offline_diag(:,:)  ! HO2 diagnostic

! The CHCH definitions array lists the species used by the chemical solver.
! The elements of the chch_t type are:
! idum (reference integer)
! name (species name)
! nodd (unused)
! ctype (either TR - tracer or CF - constant field here)
! family (unused)
! switch1 (1 if species dry deposits)
! switch2 (1 if species wet deposits)

TYPE(chch_t) :: chch_defs_a(11) = [                                            &
  chch_t( 1,'H2O2      ',  1,'TR        ','          ',  1,  1),               &
  chch_t( 2,'DMS       ',  1,'TR        ','          ',  0,  0),               &
  chch_t( 3,'SO2       ',  1,'TR        ','          ',  1,  1),               &
  chch_t( 4,'H2SO4     ',  1,'TR        ','          ',  1,  0),               &
  chch_t( 5,'DMSO      ',  1,'TR        ','          ',  1,  1),               &
  chch_t( 6,'Monoterp  ',  1,'TR        ','          ',  1,  0),               &
  chch_t( 7,'Sec_Org   ',  1,'TR        ','          ',  1,  1),               &
  chch_t( 8,'O3        ',  1,'CF        ','          ',  0,  0),               &
  chch_t( 9,'OH        ',  1,'CF        ','          ',  0,  0),               &
  chch_t(10,'NO3       ',  1,'CF        ','          ',  0,  0),               &
  chch_t(11,'HO2       ',  1,'CF        ','          ',  0,  0)                &
  ]


! Bimolecular reactions
TYPE(ratb_t), PARAMETER :: ratb_defs_a(9) = [                                  &
  ratb_t('DMS       ','OH        ','SO2       ','          ','          ',     &
  '          ',9.60e-12,  0.00,    240.00, 0.000, 0.000, 0.000, 0.000),        &
  ratb_t('DMS       ','OH        ','SO2       ','DMSO      ','          ',     &
  '          ',3.04e-12,  0.00,   -350.00, 0.600, 0.400, 0.000, 0.000),        &
  ratb_t('DMS       ','NO3       ','SO2       ','          ','          ',     &
  '          ',1.90e-13,  0.00,   -500.00, 0.000, 0.000, 0.000, 0.000),        &
  ratb_t('DMSO      ','OH        ','SO2       ','          ','          ',     &
  '          ',5.80e-11,  0.00,      0.00, 0.600, 0.000, 0.000, 0.000),        &
  ratb_t('Monoterp  ','OH        ','Sec_Org   ','          ','          ',     &
  '          ',1.20e-11,  0.00,   -444.00, 0.130, 0.000, 0.000, 0.000),        &
  ratb_t('Monoterp  ','O3        ','Sec_Org   ','          ','          ',     &
  '          ',1.01e-15,  0.00,   732.00, 0.130, 0.000, 0.000, 0.000),         &
  ratb_t('Monoterp  ','NO3       ','Sec_Org   ','          ','          ',     &
  '          ',1.19e-12,  0.00,   -925.00, 0.130, 0.000, 0.000, 0.000),        &
  ratb_t('HO2       ','HO2       ','H2O2      ','          ','          ',     &
  '          ',2.20e-13,  0.00,   -600.00, 0.000, 0.000, 0.000, 0.000),        &
  ratb_t('OH        ','H2O2      ','H2O       ','          ','          ',     &
  '          ',2.90e-12,  0.00,    160.00, 0.000, 0.000, 0.000, 0.000)         &
  ]

! The rate coefficients for these reactions are
! defined in the asad_hetero routine. The 'NULL' products serve to identify the
! reactions - these reactions have no gas-phase products.
TYPE(rath_t) :: rath_defs_a(naq_offline) = [                                   &
  rath_t('SO2       ','H2O2      ','NULL0     ','          ','          ',     &
  '          ', 0.000, 0.000, 0.000, 0.000),                                   &
  ! HSO3+H2O2(aq)
  rath_t('SO2       ','O3        ','NULL1     ','          ','          ',     &
  '          ', 0.000, 0.000, 0.000, 0.000),                                   &
  ! HSO3+O3(aq)
  rath_t('SO2       ','O3        ','NULL2     ','          ','          ',     &
  '          ', 0.000, 0.000, 0.000, 0.000)                                    &
  ! SO3+O3(aq)
  ]


! Termolecular reactions
TYPE(ratt_t) :: ratt_defs_a(1) = [                                             &
  ratt_t('SO2       ','OH        ','H2SO4     ','          ',                  &
  0.6, 3.00e-31, -3.30,     0.0, 1.50e-12, 0.00, 0.0, 0.000, 0.000)            &
  ]

! Dry deposition array:
! NOTE: This array is only used when the interactives scheme is not selected
! Dry deposition rates are written in a table format, with the columns
! corresponding to times and the rows corresponding to surface type:
!
!           SUMMER                       WINTER
!       day  night  24hr    day  night  24hr
!         x     x     x       x     x     x     Water
!         x     x     x       x     x     x     Forest
!         x     x     x       x     x     x     Grass/Shrub
!         x     x     x       x     x     x     Desert
!         x     x     x       x     x     x     Snow/Ice
!
! '24 hr' corresponds to the 24 hour average deposition velocity,
!  currently not used. The units of x are cm s-1

! These have to be in the same order as the species which dry deposit in the
! chch array.


REAL :: depvel_defs_a(6,5,ndry_offline) = RESHAPE([                            &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & ! H2O2
  1.25,  0.16,  0.71,  0.28,  0.12,  0.20,                                     &
  1.25,  0.53,  0.89,  0.83,  0.78,  0.81,                                     &
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,                                     &
  0.32,  0.32,  0.32,  0.32,  0.32,  0.32,                                     &
  0.80,  0.80,  0.80,  0.80,  0.80,  0.80,  & ! SO2
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,                                     &
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,                                     &
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,                                     &
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,                                     &
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  & ! H2SO4
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,                                     &
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,                                     &
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,                                     &
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,                                     &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & ! DMSO
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                                     &
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  & ! Monoterp
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,                                     &
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,                                     &
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,                                     &
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,                                     &
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  & ! Sec_Org
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,                                     &
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,                                     &
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,                                     &
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10                                      &
  ],[  6,  5,  ndry_offline] )


! The following formula is used to calculate the effective Henry's Law coef,
! which takes the affects of dissociation and complex formation on a species'
! solubility (see Giannakopoulos, 1998)
!
!       H(eff) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
!
! The data in columns 1 and 2 above give the data for gas-aqueous transfer,
!       Column 1 = K(298) [M/atm]
!       Column 2 = -deltaH/R [K-1]
!
! If the species dissociates in the aqueous phase the above term is multiplied
! by another factor of 1+{K(aq)/[H+]}, where
!       K(aq) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
! The data in columns 3 and 4 give the data for this aqueous-phase dissociation,
!       Column 3 = K(298) [M]
!       Column 4 = -deltaH/R [K-1]
! The data in columns 5 and 6 give the data for a second dissociation,
! e.g for SO2, HSO3^{-}, and SO3^{2-}
!       Column 5 = K(298) [M]
!       Column 6 = -deltaH/R [K-1]

! Tracers (must be in the same order as in chch_defs array):
!   H2O2, SO2, DMSO, Sec_Org
REAL :: henry_defs_a(6,nwet_offline) = RESHAPE([                               &
  0.830e+05, 0.740e+04, 0.240e-11,-0.373e+04, 0.000e+00, 0.000e+00, & ! H2O2
  0.123e+01, 0.302e+04, 0.123e-01, 0.201e+04, 0.600e-07, 0.112e+04, & ! SO2
  0.500e+05, 0.6425e+04,0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, & ! DMSO
  0.100e+06, 0.120e+02, 0.000e-00, 0.000e+00, 0.000e+00, 0.000e+00  & ! Sec_Org
  ],[  6,  nwet_offline] )

! Constant fields: O3
REAL, PUBLIC :: henry_defs_const(6,nwet_constant) = RESHAPE([                  &
     0.113e-01, 0.230e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00          &
     ],[  6,  nwet_constant] )

PUBLIC ukca_init_offline

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CHEM_OFFLINE'

CONTAINS

! ######################################################################

SUBROUTINE ukca_init_offline

! To initialise chemistry files from specified species and rate arrays
! using the defined logicals.

USE asad_mod,    ONLY: jddepc, jddept, ih_o3_const,                            &
                       ct_k298, ct_dhr, ct_kd298, ct_ddhr, jpeq, jpsp,         &
                       jpctr, jpspec, jpbk, jppj, jptk, jphk, jpdd, jpdw
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr,  ONLY: umPrint, umMessage
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CHARACTER (LEN=errormessagelength) :: cmessage           ! Error message
INTEGER            :: errcode            ! Error code
INTEGER            :: i                  ! Loop counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INIT_OFFLINE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Define species index for the constant aqueous species
ih_o3_const = 1

n_tracers = 0
! Count no of tracers in chch_defs array
DO i = 1,SIZE(chch_defs_a)
  IF (chch_defs_a(i)%ctype == jpsp) n_tracers = n_tracers + 1
END DO

! Check that number of tracers declared is equal to jpctr
IF (n_tracers /= jpctr) THEN
  cmessage = 'No. of tracers found not equal to jpctr'
  errcode = 1
  WRITE(umMessage,'(A40,A12,I5,A9,I5)') cmessage,' n_tracers:',n_tracers,      &
    ' jpctr: ',jpctr
  CALL umPrint(umMessage,src='UKCA_INIT_OFFLINE')
  CALL ereport('UKCA_INIT_OFFLINE',errcode,cmessage)
END IF

ALLOCATE(chch_defs_offline(jpspec))
ALLOCATE(ratb_defs_offline(jpbk))
ALLOCATE(ratj_defs_offline(jppj))
ALLOCATE(ratt_defs_offline(jptk))
ALLOCATE(rath_defs_offline(jphk))
ALLOCATE(depvel_defs_offline(jddept,jddepc,jpdd))
ALLOCATE(henry_defs_offline(6,jpdw))

! Set number of constant species in aqueous phase
IF (nwet_constant > 0) THEN
  ! allow for dissociation for constant species as well as tracers
  IF (.NOT. ALLOCATED(ct_k298))  ALLOCATE(ct_k298(nwet_constant))
  IF (.NOT. ALLOCATED(ct_dhr))   ALLOCATE(ct_dhr(nwet_constant))
  IF (.NOT. ALLOCATED(ct_kd298)) ALLOCATE(ct_kd298(nwet_constant,jpeq))
  IF (.NOT. ALLOCATED(ct_ddhr))  ALLOCATE(ct_ddhr(nwet_constant,jpeq))
END IF

! Offline chemistry
chch_defs_offline = [chch_defs_a]
ratb_defs_offline = [ratb_defs_a]
ratt_defs_offline = [ratt_defs_a]
rath_defs_offline = [rath_defs_a]
depvel_defs_offline = depvel_defs_a
henry_defs_offline  = henry_defs_a

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_init_offline

END MODULE ukca_chem_offline
