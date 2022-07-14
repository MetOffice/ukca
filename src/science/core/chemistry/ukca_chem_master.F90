! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module providing UKCA with inputs describing chemistry
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_chem_master_mod

IMPLICIT NONE

PRIVATE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CHEM_MASTER_MOD'

PUBLIC  :: ukca_chem_master

CONTAINS

SUBROUTINE ukca_chem_master

USE ukca_chem_defs_mod, ONLY: chch_t1, ratb_t1, rath_t1, ratj_t1, ratt_t1,     &
                              depvel_t, wetdep,                                &
                              chch_defs_new, ratb_defs_new, rath_defs_new,     &
                              ratj_defs_new, ratt_defs_new, depvel_defs_new,   &
                              henry_defs_new
USE ukca_config_specification_mod, ONLY: ukca_config
USE asad_mod,           ONLY: jddepc, jddept, jpctr, jpspec, jpro2, jpcspf,    &
                              jpbk, jptk, jppj, jphk, jpnr, jpdd, jpdw
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,  ONLY: umPrint, umMessage

IMPLICIT NONE


! In the below items, item_no identifies different versions of the same entry
! (e.g., updates of the same reaction). If e.g. in different schemes, the
! branching ratios or products are different, these should be identified
! through different settings of chemistry. If multiple versions are defined,
! the highest version <= the selected version is used. Versions of items need
! to be listed together.


! Define chemistry scheme types
! Entries will only be selected if the active chemistry scheme flag matches
! with one of the values in the "scheme" column.
! Only one chemistry scheme flag can be active at a time.
INTEGER, PARAMETER :: st  = 1   ! stratosphere-troposphere scheme
INTEGER, PARAMETER :: t   = 2   ! troposphere scheme
INTEGER, PARAMETER :: s   = 4   ! stratosphere scheme
INTEGER, PARAMETER :: r   = 8   ! RAQ scheme
INTEGER, PARAMETER :: ol  = 16  ! offline scheme
INTEGER, PARAMETER :: ti  = 32  ! troposphere-isoprene scheme
INTEGER, PARAMETER :: cs  = 64  ! Common Representative Intermediates +
                                ! Stratosphere scheme

! CRI-Strat 2 (CS2) is an update to CRI-Strat (CS) which is turned on if
! CRI-Strat is selected and i_ukca_chem_version >= 119

! Define qualifiers
! Entries with the corresponding qualifier flag in their 'qualifier' column
! will be selected, whilst those with the corresponding flag in their
! 'disqualifier' column will be deselected. Any entry with no value in
! the qualifier and disqualifier columns will be automatically selected.
! Multiple qualifier flags can be used at the same time. However, a match
! in the disqualifier column will take precedent over a match in the
! qualifier column.
INTEGER, PARAMETER :: a  = 1      ! aerosol chemistry
INTEGER, PARAMETER :: th = 2      ! tropospheric heterogeneous reactions
INTEGER, PARAMETER :: hp = 4      ! heterogeneous PSC chemistry
INTEGER, PARAMETER :: es = 8      ! extended stratospheric reactions
INTEGER, PARAMETER :: ci = 16     ! 3-D CO2
INTEGER, PARAMETER :: eh = 32     ! extended heterogeneous chemistry
INTEGER, PARAMETER :: rn = 64     ! RO2 Non-transported prognostics (ST only)
INTEGER, PARAMETER :: rp = 128    ! RO2-Permutation chemistry (ST only)
INTEGER, PARAMETER :: cs2 = 256   ! Used to deselect species which are
                                  ! present in CS but not in CS2
INTEGER, PARAMETER :: ce = 512    ! CH4 emission-driven chemistry
INTEGER, PARAMETER :: rm = 1024   ! Used to remove reactions that existed at
                                  ! previous chemistry versions but not at
                                  ! later ones

! Define size of master chemistry
INTEGER, PARAMETER :: n_chch_master = 356 ! number of known species
INTEGER, PARAMETER :: n_het_master  =  18 ! number of heterogeneous reactions
INTEGER, PARAMETER :: n_dry_master  = 161 ! number of dry deposition reactions
INTEGER, PARAMETER :: n_wet_master  = 159 ! number of wet deposition reactions
INTEGER, PARAMETER :: n_bimol_master = 1209 ! number of bimolecular reactions
INTEGER, PARAMETER :: n_ratj_master = 183 ! number of photolysis reactions
INTEGER, PARAMETER :: n_ratt_master = 116 ! number of termolecular reactions

! Start and end ranges for each section of arrays. Makes adding extra reactions
! to multi-table arrays easier
INTEGER :: n_ratb_s, n_ratb_e


! Arrays split to avoid continuations exceeded error
! All declarations moved up to here to avoid type-declaration
! out of order errors
TYPE(chch_t1)  :: chch_defs_master(1:n_chch_master)
TYPE(rath_t1)  :: rath_defs_master(1:n_het_master)
TYPE(ratj_t1)  :: ratj_defs_master(1:n_ratj_master)
TYPE(ratt_t1)  :: ratt_defs_master(1:n_ratt_master)
TYPE(wetdep)   :: henry_defs_master(1:n_wet_master)
TYPE(depvel_t) :: depvel_defs_master(1:n_dry_master)
TYPE(ratb_t1)  :: ratb_defs_master(1:n_bimol_master)

! To initialise chemistry files from specified species and rate arrays
! using the defined logicals. Composes chemistry definition arrays from
! components using logical definitions from UKCA configuration data.
INTEGER            :: ierr
CHARACTER (LEN=errormessagelength) :: cmessage = ''

INTEGER :: chem_scheme
INTEGER :: qualifier
INTEGER :: i
INTEGER :: j
INTEGER :: k
LOGICAL :: found
LOGICAL :: take_this
TYPE(chch_t1) :: last_chch
TYPE(ratb_t1) :: last_bimol
TYPE(ratt_t1) :: last_termol
TYPE(ratj_t1) :: last_photol
TYPE(rath_t1) :: last_het
TYPE(depvel_t):: last_depvel
TYPE(wetdep)  :: last_wetdep

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEM_MASTER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Define below the chemistry schemes

! ST scheme follows ATA NLA CheST Chemistry v1.2
! These reactions generally have comments identifying where they come from
! (e.g., "T001" is the first termolecular reaction. Reactions not labelled are
! taken from the previous "trop", "tropisop", "strat", or offline schemes,
! where they were not labelled (hence I'm not sure about the origins of the
! data). I suggest that eventually these reactions are superseded by the ST
! scheme.


! Species Definitions
! The columns take the following meanings:
! item number, species name, (unused), species classification, family (unused),
! dry deposition flag, wet deposition flag, scheme,
! selection qualifier, deselection / disqualifier, (UM) version when entered
! into the master mechanism.
chch_defs_master(1:149)=[                                                      &
!   1
chch_t1( 1,'O(3P)     ',1,'TR        ','Ox        ',0,0,s+st+cs,0,0,107),      &
!   1
chch_t1( 1,'O(3P)     ',1,'SS        ','Ox        ',0,0,ti+t+r,0,0,107),       &
!   2
chch_t1( 2,'O(1D)     ',1,'SS        ','Ox        ',0,0,ti+s+t+st+r+cs,        &
                                                                0,0,107),      &
!   3 DD:  1,
chch_t1( 3,'O3        ',1,'TR        ','Ox        ',1,0,ti+s+t+st+r+cs,        &
                                                                0,a,107),      &
!   3        WD:  1,
!      Wet deposition of ozone was introduced for aerosol chemistry
!      The GLOMAP-mode aerosol scheme needs to know about dissolved O3 in the
!      aqueous phase to be able to perform SO2 in-cloud oxidation.
!      Wet deposition is a negligible loss process for O3.
chch_t1( 3,'O3        ',1,'TR        ','Ox        ',1,1,ti+st+s+t+r+cs,        &
                                                                a,0,107),      &
!   3
chch_t1( 3,'O3        ',1,'CF        ','Ox        ',0,0,ol,0,0,107),           &
!   4
chch_t1( 4,'N         ',1,'TR        ','NOx       ',0,0,s+st+cs,0,0,107),      &
! NOx is emitted as NO2 for all schemes except RAQ which causes ozone
! depletion (NO + O3 -> NO2 + O2).
! This could be made consistent. Also RAQ scheme does not do dry deposition of
! NO.
!   5 DD:  2,
chch_t1( 5,'NO        ',1,'TR        ','NOx       ',1,0,s+st+cs,0,0,107),      &
chch_t1( 5,'NO        ',1,'TR        ','NOx       ',1,0,ti+t,0,0,107),         &
!   5                 (10)
chch_t1( 5,'NO        ',1,'TR        ','NOx       ',0,0,r,0,0,107),            &
!   6 DD:  3,WD:  2,
! No dry deposition for NO3 in RAQ/S/T schemes.
chch_t1( 6,'NO3       ',1,'TR        ','NOx       ',0,1,r+s+t,0,0,107),        &
chch_t1( 6,'NO3       ',1,'TR        ','NOx       ',1,1,ti+st+cs,0,0,107),     &
chch_t1( 6,'NO3       ',1,'CF        ','NOx       ',0,0,ol,0,0,107),           &
!   7 DD:  4,
chch_t1( 7,'NO2       ',1,'TR        ','NOx       ',1,0,ti+s+st+cs,0,0,107),   &
chch_t1( 7,'NO2       ',1,'TR        ','NOx       ',1,0,t+r,0,0,107),          &
! No dry deposition for N2O5 in R/S/T schemes
!   8 DD:  5,WD:  3,
chch_t1( 8,'N2O5      ',2,'TR        ','          ',0,1,r+s+t,0,0,107),        &
chch_t1( 8,'N2O5      ',2,'TR        ','          ',1,1,ti+st+cs,0,0,107),     &
!   9 DD:  6,WD:  4,
! No dry deposition for HO2NO2 in R and T schemes
chch_t1( 9,'HO2NO2    ',1,'TR        ','          ',0,1,r+t,0,0,107),          &
chch_t1( 9,'HO2NO2    ',1,'TR        ','          ',1,1,ti+s+st+cs,0,0,107),   &
!  10 DD:  7,WD:  5,   (20)
chch_t1(10,'HONO2     ',1,'TR        ','          ',1,1,ti+s+t+st+r+cs,        &
                                                                0,0,107),      &
!  11 DD:  8,WD:  6,
chch_t1(11,'H2O2      ',1,'TR        ','          ',1,1,ti+s+t+st+ol+r+cs,     &
                                                                   0,0,107),   &
!  12
! Dry deposition of CH4 not enabled in S and ST schemes.
chch_t1(12,'CH4       ',1,'TR        ','          ',0,0,s+st+cs,0,ce,107),     &
!  12 DD:  9,
chch_t1(12,'CH4       ',1,'TR        ','          ',1,0,ti+t+r,0,ce,107),      &
!  12
! Dry deposition of CH4 not enabled in S and CS schemes.
chch_t1(12,'CH4       ',1,'TR        ','          ',0,0,s+cs,ce,0,107),        &
!  12 DD:  9,
! Dry deposition of CH4 enabled in TI, T, ST and R schemes.
chch_t1(12,'CH4       ',1,'TR        ','          ',1,0,ti+t+st+r,ce,0,107),   &
!  13 DD: 10,
chch_t1(13,'CO        ',1,'TR        ','          ',1,0,ti+s+t+st+r+cs,        &
                                                                0,0,107),      &
!  14 DD: 11,WD:  7,
! No dry deposition for HCHO in S/T/R schemes
chch_t1(14,'HCHO      ',1,'TR        ','          ',0,1,s+t+r,0,0,107),        &
chch_t1(14,'HCHO      ',1,'TR        ','          ',1,1,ti+st+cs,0,0,107),     &
!  15        WD:  8,
chch_t1(15,'MeOO      ',1,'TR        ','          ',0,1,ti+s,0,0,107),         &
chch_t1(15,'MeOO      ',1,'SS        ','          ',0,1,t+r,0,0,107),          &
chch_t1(15,'MeOO      ',1,'TR        ','          ',0,1,st,0,rn+rp,107),       &
! Select 'OO' type MeOO in ST mechanism if RN or RP == T
!                     (30)
chch_t1(15,'MeOO      ',1,'OO        ','          ',0,1,st+cs,rn+rp,0,107),    &
!  16 DD: 12,WD:  9,
chch_t1(16,'MeOOH     ',1,'TR        ','          ',1,1,ti+s+t+st+r+cs,        &
    0,0,107),                                                                  &
!  17
chch_t1(17,'H         ',1,'TR        ','HOx       ',0,0,s+st+cs,0,0,107),      &
!  18
! No water in RAQ scheme? I presume H2O is implicit.
chch_t1(18,'H2O       ',1,'TR        ','          ',0,0,s+st+cs,0,0,107),      &
!  18
chch_t1(18,'H2O       ',1,'CF        ','          ',0,0,ti+t+r,0,0,107),       &
!  19
chch_t1(19,'OH        ',1,'TR        ','HOx       ',0,0,ti+s+st+cs,0,0,107),   &
chch_t1(19,'OH        ',1,'SS        ','HOx       ',0,0,t+r,0,0,107),          &
chch_t1(19,'OH        ',1,'CF        ','HOx       ',0,0,ol,0,0,107),           &
!  20        WD: 10,
chch_t1(20,'HO2       ',1,'TR        ','HOx       ',0,1,ti+s+st+cs,0,0,107),   &
chch_t1(20,'HO2       ',1,'SS        ','HOx       ',0,1,t+r,0,0,107),          &
chch_t1(20,'HO2       ',1,'CF        ','HOx       ',0,0,ol,0,0,107),           &
!  21
chch_t1(21,'Cl        ',1,'TR        ','Clx       ',0,0,s+st+cs,0,0,107),      &
!  22
chch_t1(22,'Cl2O2     ',1,'TR        ','Clx       ',0,0,s+st+cs,0,0,107),      &
!  23
chch_t1(23,'ClO       ',1,'TR        ','Clx       ',0,0,s+st+cs,0,0,107),      &
!  24
chch_t1(24,'OClO      ',1,'TR        ','          ',0,0,s+st+cs,0,0,107),      &
!  25
chch_t1(25,'Br        ',1,'TR        ','Brx       ',0,0,s+st+cs,0,0,107),      &
!  26
chch_t1(26,'BrO       ',1,'TR        ','Brx       ',0,0,s+st+cs,0,0,107),      &
!  27
chch_t1(27,'BrCl      ',1,'TR        ','          ',0,0,s+st+cs,0,0,107),      &
!  28        WD: 11,
chch_t1(28,'BrONO2    ',1,'TR        ','          ',0,1,s+st+cs,0,0,107),      &
!  29
chch_t1(29,'N2O       ',1,'TR        ','          ',0,0,s+st+cs,0,0,107),      &
!  30 DD: 13,WD: 12,   (50)
chch_t1(30,'HCl       ',1,'TR        ','          ',1,1,s+st+cs,0,0,107),      &
!  31 DD: 14,WD: 13,
chch_t1(31,'HOCl      ',1,'TR        ','          ',1,1,s+st+cs,0,0,107),      &
!  32 DD: 15,WD: 14,
chch_t1(32,'HBr       ',1,'TR        ','          ',1,1,s+st+cs,0,0,107),      &
!  33 DD: 16,WD: 15,
chch_t1(33,'HOBr      ',1,'TR        ','          ',1,1,s+st+cs,0,0,107),      &
!  34        WD: 16,
chch_t1(34,'ClONO2    ',1,'TR        ','          ',0,1,s+st+cs,0,0,107),      &
!  35
chch_t1(35,'CFCl3     ',1,'TR        ','          ',0,0,s+st+cs,0,0,107),      &
!  36
chch_t1(36,'CF2Cl2    ',1,'TR        ','          ',0,0,s+st+cs,0,0,107),      &
!  37
chch_t1(37,'MeBr      ',1,'TR        ','          ',0,0,s+st+cs,0,0,107),      &
!  38 DD: 17,WD: 17,
chch_t1(38,'HONO      ',1,'TR        ','          ',1,1,ti+t+st+cs,0,0,107),   &
!  39
chch_t1(39,'C2H6      ',1,'TR        ','          ',0,0,ti+t+st+r+cs,          &
                                                              0,0,107),        &
!  40
chch_t1(40,'EtOO      ',1,'TR        ','          ',0,0,ti,0,0,107),           &
!  40
chch_t1(40,'EtOO      ',1,'SS        ','          ',0,0,t+r,0,0,107),          &
!  40 ; Select 'OO' type EtOO in ST mechanism if RN or RP == T
chch_t1(40,'EtOO      ',1,'TR        ','          ',0,0,st,0,rn+rp,107),       &
chch_t1(40,'EtOO      ',1,'OO        ','          ',0,0,st+cs,rn+rp,0,107),    &
!  41 DD: 18,WD: 18,
chch_t1(41,'EtOOH     ',1,'TR        ','          ',1,1,ti+t+st+r+cs,          &
                                                              0,0,107),        &
!  42 DD: 19,
! No dry deposition for MeCHO in R and T schemes
chch_t1(42,'MeCHO     ',1,'TR        ','          ',0,0,r+t,0,0,107),          &
chch_t1(42,'MeCHO     ',1,'TR        ','          ',1,0,ti+st+cs,0,0,107),     &
!  43
chch_t1(43,'MeCO3     ',1,'TR        ','          ',0,0,ti,0,0,107),           &
!  43
chch_t1(43,'MeCO3     ',1,'SS        ','          ',0,0,t+r,0,0,107),          &
!  43 ; Select 'OO' type MeCO3 in ST mechanism if RN or RP == T
chch_t1(43,'MeCO3     ',1,'TR        ','          ',0,0,st,0,rn+rp,107),       &
!  43                 (70)
chch_t1(43,'MeCO3     ',1,'OO        ','          ',0,0,st+cs,rn+rp,0,107),    &
!  44 DD: 20,
chch_t1(44,'PAN       ',1,'TR        ','          ',1,0,ti+t+st+r+cs,          &
                                                              0,0,107),        &
!  45
chch_t1(45,'C3H8      ',1,'TR        ','          ',0,0,ti+t+st+r+cs,          &
                                                              0,0,107),        &
!  46
chch_t1(46,'n-PrOO    ',1,'TR        ','          ',0,0,ti,0,0,107),           &
!  46
chch_t1(46,'n-PrOO    ',1,'SS        ','          ',0,0,t,0,0,107),            &
!  46 ; Select 'OO' type n-PrOO in ST mechanism if RN or RP == T
chch_t1(46,'n-PrOO    ',1,'TR        ','          ',0,0,st,0,rn+rp,107),       &
chch_t1(46,'n-PrOO    ',1,'OO        ','          ',0,0,st,rn+rp,0,107),       &
!  47
chch_t1(47,'i-PrOO    ',1,'TR        ','          ',0,0,ti,0,0,107),           &
!  47
chch_t1(47,'i-PrOO    ',1,'SS        ','          ',0,0,t+r,0,0,107),          &
!  47 ; Select 'OO' type i-PrOO in ST mechanism if RN or RP == T
chch_t1(47,'i-PrOO    ',1,'TR        ','          ',0,0,st,0,rn+rp,107),       &
!                     (80)
chch_t1(47,'i-PrOO    ',1,'OO        ','          ',0,0,st+cs,rn+rp,0,107),    &
!  48 DD: 21,WD: 19,
chch_t1(48,'n-PrOOH   ',1,'TR        ','          ',1,1,ti+t+st,0,0,107),      &
!  49 DD: 22,WD: 20,
chch_t1(49,'i-PrOOH   ',1,'TR        ','          ',1,1,ti+t+st+r+cs,          &
                                                              0,0,107),        &
!  50 DD: 23,
chch_t1(50,'EtCHO     ',1,'TR        ','          ',1,0,ti+t+st,0,0,107),      &
!  50
chch_t1(50,'EtCHO     ',1,'TR        ','          ',1,0,cs,0,0,107),           &
!  51
chch_t1(51,'EtCO3     ',1,'TR        ','          ',0,0,ti,0,0,107),           &
!  51
chch_t1(51,'EtCO3     ',1,'SS        ','          ',0,0,t,0,0,107),            &
!  51 ; Select 'OO' type EtCO3 in ST mechanism if RN or RP == T
chch_t1(51,'EtCO3     ',1,'TR        ','          ',0,0,st,0,rn+rp,107),       &
chch_t1(51,'EtCO3     ',1,'OO        ','          ',0,0,st+cs,rn+rp,0,107),    &
!  52
chch_t1(52,'Me2CO     ',1,'TR        ','          ',0,0,ti+t+st+r+cs,          &
                                                              0,0,107),        &
!  53                 (90)
chch_t1(53,'MeCOCH2OO ',1,'TR        ','          ',0,0,ti,0,0,107),           &
!  53
chch_t1(53,'MeCOCH2OO ',1,'SS        ','          ',0,0,t+r,0,0,107),          &
!  53 ; Select 'OO' type MeCOCH2OO in ST mechanism if RN or RP == T
chch_t1(53,'MeCOCH2OO ',1,'TR        ','          ',0,0,st,0,rn+rp,107),       &
chch_t1(53,'MeCOCH2OO ',1,'OO        ','          ',0,0,st,rn+rp,0,107),       &
!  54 DD: 24,WD: 21,
chch_t1(54,'MeCOCH2OOH',1,'TR        ','          ',1,1,ti+t+st,0,0,107),      &
!  55 DD: 25,
chch_t1(55,'PPAN      ',1,'TR        ','          ',1,0,ti+t+st+cs,0,0,107),   &
!  56
chch_t1(56,'MeONO2    ',1,'TR        ','          ',0,0,ti+t+st+cs,0,0,107),   &
!  57
chch_t1(57,'C5H8      ',1,'TR        ','          ',0,0,ti+st+r+cs,0,0,107),   &
!  58
chch_t1(58,'ISO2      ',1,'TR        ','          ',0,0,ti,0,0,107),           &
!  58 ; Select 'OO' type ISO2 in ST mechanism if RN or RP == T
chch_t1(58,'ISO2      ',1,'TR        ','          ',0,0,st,0,rn+rp,107),       &
!  58                 (100)
chch_t1(58,'ISO2      ',1,'OO        ','          ',0,0,st,rn+rp,0,107),       &
!  59 DD: 26,WD: 22,
chch_t1(59,'ISOOH     ',1,'TR        ','          ',1,1,ti+st+r,0,0,107),      &
!  60 DD: 27,WD: 23,
! No DD for ISON in R scheme
chch_t1(60,'ISON      ',1,'TR        ','          ',0,1,r,0,0,107),            &
chch_t1(60,'ISON      ',1,'TR        ','          ',1,1,ti+st,0,0,107),        &
!  61 DD: 28,
chch_t1(61,'MACR      ',1,'TR        ','          ',1,0,ti+st,0,0,107),        &
!  62
chch_t1(62,'MACRO2    ',1,'TR        ','          ',0,0,ti,0,0,107),           &
!  62 ; Select 'OO' type MACRO2 in ST mechanism if RN or RP == T
chch_t1(62,'MACRO2    ',1,'TR        ','          ',0,0,st,0,rn+rp,107),       &
chch_t1(62,'MACRO2    ',1,'OO        ','          ',0,0,st,rn+rp,0,107),       &
!  63 DD: 29,WD: 24,
chch_t1(63,'MACROOH   ',1,'TR        ','          ',1,1,ti+st,0,0,107),        &
!  64 DD: 30,
chch_t1(64,'MPAN      ',1,'TR        ','          ',1,0,ti+st+cs,0,0,107),     &
!  65 DD: 31,WD: 25,   (110)
chch_t1(65,'HACET     ',1,'TR        ','          ',1,1,ti+st,0,0,107),        &
!  66 DD: 32,WD: 26,
! No dry deposition for methyl glyoxal in R scheme
chch_t1(66,'MGLY      ',1,'TR        ','          ',0,1,r,0,0,107),            &
chch_t1(66,'MGLY      ',1,'TR        ','          ',1,1,ti+st,0,0,107),        &
!  67 DD: 33,
chch_t1(67,'NALD      ',1,'TR        ','          ',1,0,ti+st,0,0,107),        &
!  68 DD: 34,WD: 27,
chch_t1(68,'HCOOH     ',1,'TR        ','          ',1,1,ti+st,0,0,107),        &
!  68
chch_t1(68,'HCOOH     ',1,'TR        ','          ',1,1,cs,0,0,107),           &
!  69 DD: 35,WD: 28,
chch_t1(69,'MeCO3H    ',1,'TR        ','          ',1,1,ti+st+cs,0,0,107),     &
!  70 DD: 36,WD: 29,
chch_t1(70,'MeCO2H    ',1,'TR        ','          ',1,1,ti+st,0,0,107),        &
!  70 MeCO2H emitted in CRI
chch_t1(70,'MeCO2H    ',1,'TR        ','          ',1,1,cs,0,0,107),           &
!  71 DD: 37,
! Dry deposition for H2 in R scheme but not in other schemes
! H2 is prescribed at the surface. Could it become interactive?
chch_t1(71,'H2        ',1,'TR        ','          ',0,0,s+st+cs,0,0,107),      &
!  71                 (120)
chch_t1(71,'H2        ',1,'TR        ','          ',1,0,r,0,0,107),            &
!  71
chch_t1(71,'H2        ',1,'CT        ','          ',0,0,ti+t,0,0,107),         &
!  72 DD: 38,WD: 30,
! No dry deposition for MeOH in R scheme. No emissions of MeOH in other schemes.
! Note MeOH is emitted as NVOC in ST with aerosols.
chch_t1(72,'MeOH      ',1,'TR        ','          ',1,1,ti+st,0,0,107),        &
chch_t1(72,'MeOH      ',1,'TR        ','          ',0,1,r,0,0,107),            &
!  72
chch_t1(72,'MeOH      ',1,'TR        ','          ',1,1,cs,0,0,107),           &
!  73
chch_t1(73,'CO2       ',1,'CT        ','          ',0,0,ti+t+s+st+cs,          &
                                                             0,ci,107),        &
!  73
! Stratospheric chemical schemes are enabled to take 3D CO2 as input
chch_t1(73,'CO2       ',1,'CF        ','          ',0,0,s+st+cs,ci,0,111),     &
!  74
chch_t1(74,'O2        ',1,'CT        ','          ',0,0,ti+s+t+st+cs,          &
                                                              0,0,107),        &
!  75
chch_t1(75,'N2        ',1,'CT        ','          ',0,0,ti+s+t+st+cs,          &
                                                              0,0,107),        &
!  76
chch_t1(76,'DMS       ',1,'TR        ','          ',0,0,ti+s+st+ol+r+cs,       &
                                                                 a,0,107),     &
!  77 DD: 39,WD: 31,
chch_t1(77,'SO2       ',1,'TR        ','          ',1,1,ti+s+st+ol+r+cs,       &
                                                                 a,0,107),     &
!  78 DD: 40,
chch_t1(78,'H2SO4     ',1,'TR        ','          ',1,0,ti+s+st+ol+r+cs,       &
                                                                 a,0,107),     &
!  79 DD: 41,
chch_t1(79,'MSA       ',1,'TR        ','          ',0,0,s+st+r+cs,a,0,107),    &
chch_t1(79,'MSA       ',1,'TR        ','          ',1,0,ti,a,0,107),           &
!  80 DD: 42,WD: 32,
! No dry deposition for DMSO in R scheme.
chch_t1(80,'DMSO      ',1,'TR        ','          ',0,1,r,a,0,107),            &
chch_t1(80,'DMSO      ',1,'TR        ','          ',1,1,ti+st+ol+cs,a,0,107),  &
!  81 DD: 43,WD: 33,
chch_t1(81,'NH3       ',1,'TR        ','          ',1,1,ti+st+r+cs,a,0,107),   &
!  82
chch_t1(82,'CS2       ',1,'TR        ','          ',0,0,ti+s+st+cs,a,0,107),   &
!  83
chch_t1(83,'COS       ',1,'TR        ','          ',0,0,ti+s+st+cs,a,0,107),   &
!  84
chch_t1(84,'H2S       ',1,'TR        ','          ',0,0,ti+s+st+cs,a,0,107),   &
!  85                 (140)
chch_t1(85,'SO3       ',1,'TR        ','          ',0,0,s+st+cs,a,0,107),      &
!  86 DD: 44,
chch_t1(86,'Monoterp  ',1,'TR        ','          ',1,0,ti+st+ol+r,a,0,107),   &
!  87 DD: 45,WD: 34,
chch_t1(87,'Sec_Org   ',1,'TR        ','          ',1,1,ti+st+ol+r+cs,         &
                                                               a,0,107),       &
chch_t1(88,'O(3P)S    ',1,'SS        ','          ',0,0,t,0,0,107),            &
chch_t1(89,'O(1D)S    ',1,'SS        ','          ',0,0,t,0,0,107),            &
! O3S does not have wet deposition but O3 does for A chemistry!
!  90 DD: 46,
chch_t1(90,'O3S       ',1,'TR        ','          ',1,0,t+r,0,0,107),          &
!  91
chch_t1(91,'OHS       ',1,'SS        ','          ',0,0,t,0,0,107),            &
!  92        WD: 35,
chch_t1(92,'HO2S      ',1,'SS        ','          ',0,1,t,0,0,107)             &
]


! Species definitions continued
chch_defs_master(150:244)=[                                                    &
!  93
chch_t1(93,'s-BuOO    ',1,'SS        ','          ',0,0,r,0,0,107),            &
!  94
chch_t1(94,'MVK       ',1,'TR        ','          ',0,0,r,0,0,107),            &
!  95 DD: 47,WD: 36,
chch_t1(95,'MVKOOH    ',1,'TR        ','          ',1,1,r,0,0,107),            &
!  96
chch_t1(96,'MEKO2     ',1,'SS        ','          ',0,0,r,0,0,107),            &
!  97
chch_t1(97,'HOC2H4O2  ',1,'SS        ','          ',0,0,r,0,0,107),            &
!  98 DD: 48,WD: 37,
chch_t1(98,'ORGNIT    ',1,'TR        ','          ',1,1,r,0,0,107),            &
!  99
chch_t1(99,'HOC3H6O2  ',1,'SS        ','          ',0,0,r,0,0,107),            &
! 100
chch_t1(100,'OXYL1     ',1,'SS        ','          ',0,0,r,0,0,107),           &
! 101
chch_t1(101,'MEMALD1   ',1,'SS        ','          ',0,0,r,0,0,107),           &
! 102
chch_t1(102,'RNC2H4    ',1,'TR        ','          ',0,0,r,0,0,107),           &
! 103
chch_t1(103,'HOIPO2    ',1,'SS        ','          ',0,0,r,0,0,107),           &
! 104
chch_t1(104,'RNC3H6    ',1,'TR        ','          ',0,0,r,0,0,107),           &
! 105
chch_t1(105,'HOMVKO2   ',1,'SS        ','          ',0,0,r,0,0,107),           &
! 106
chch_t1(106,'C2H4      ',1,'TR        ','          ',0,0,r+cs,0,0,107),        &
! 107
chch_t1(107,'C3H6      ',1,'TR        ','          ',0,0,r+cs,0,0,107),        &
! 108
chch_t1(108,'C4H10     ',1,'TR        ','          ',0,0,r+cs,0,0,107),        &
! 109 DD: 49,WD: 38,
chch_t1(109,'s-BuOOH   ',1,'TR        ','          ',1,1,r,0,0,107),           &
! 110
chch_t1(110,'MEK       ',1,'TR        ','          ',0,0,r,0,0,107),           &
! 110
chch_t1(110,'MEK       ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 111
chch_t1(111,'TOLUENE   ',1,'TR        ','          ',0,0,r+cs,0,0,107),        &
! 112
chch_t1(112,'TOLP1     ',1,'SS        ','          ',0,0,r,0,0,107),           &
! 113
chch_t1(113,'MEMALD    ',1,'TR        ','          ',0,0,r,0,0,107),           &
! 114        WD: 39,
chch_t1(114,'GLY       ',1,'TR        ','          ',0,1,r,0,0,107),           &
! 115
chch_t1(115,'oXYLENE   ',1,'TR        ','          ',0,0,r+cs,0,0,107),        &
! 116 - Lumped peroxy radical field for use with RO2-permutation chemistry
chch_t1(116,'RO2       ',1,'CF        ','          ',0,0,st+cs,rp,0,107),      &
! From here on out, all new CRI species
! 117
chch_t1(117,'TBUT2ENE  ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 118 DD: 50,
chch_t1(118,'APINENE   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 119 DD: 51,
chch_t1(119,'BPINENE   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 120
chch_t1(120,'C2H2      ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 121
chch_t1(121,'BENZENE   ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 122 DD: 52,WD: 40,
chch_t1(122,'EtOH      ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 123 DD: 53,WD: 41,
chch_t1(123,'i-PrOH    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 124 DD: 54,WD: 42,    (180)
chch_t1(124,'n-PrOH    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 125 DD: 55,WD: 43,
chch_t1(125,'HOCH2CHO  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 126
chch_t1(126,'HOCH2CH2O2',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 127 DD: 56,WD: 44,
chch_t1(127,'HOC2H4OOH ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 128
chch_t1(128,'HOCH2CO3  ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 129 DD: 57,WD: 45,
chch_t1(129,'EtCO3H    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 130 DD: 58,WD: 46,
chch_t1(130,'HOCH2CO3H ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 131 DD: 59,WD: 47,
chch_t1(131,'NOA       ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 132
chch_t1(132,'EtONO2    ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 133
chch_t1(133,'i-PrONO2  ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 134                  (190)
chch_t1(134,'MeO2NO2   ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 135 DD: 60,
chch_t1(135,'HOC2H4NO3 ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! WD + DD in CS2
chch_t1(135,'HOC2H4NO3 ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 136 DD: 61,
chch_t1(136,'PHAN      ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! WD + DD in CS2
chch_t1(136,'PHAN      ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 137
chch_t1(137,'MeSCH2OO  ',1,'TR        ','          ',0,0,cs,a,0,107),          &
! 138
chch_t1(138,'MeS       ',1,'TR        ','          ',0,0,cs,a,0,107),          &
! 139
chch_t1(139,'MeSO      ',1,'TR        ','          ',0,0,cs,a,0,107),          &
! 140
chch_t1(140,'MeSO2     ',1,'TR        ','          ',0,0,cs,a,0,107),          &
! 141
chch_t1(141,'MeSO3     ',1,'TR        ','          ',0,0,cs,a,0,107),          &
! 142
chch_t1(142,'MSIA      ',1,'TR        ','          ',0,0,cs,a,0,107),          &
! 143 DD: 62,
chch_t1(143,'CARB14    ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 144 DD: 63,          (200)
chch_t1(144,'CARB17    ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 145 DD: 64,
chch_t1(145,'CARB11A   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 146 DD: 65,WD: 48, CARB7  ~ HACET
chch_t1(146,'CARB7     ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 147 DD: 66,WD: 49, CARB10 ~ MACR
chch_t1(147,'CARB10    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 148 DD: 67,WD: 50,
chch_t1(148,'CARB13    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 149 DD: 68,WD: 51,
chch_t1(149,'CARB16    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 150 DD: 69,
chch_t1(150,'UCARB10   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 151 DD: 70,WD: 52, CARB3  ~ GLY
chch_t1(151,'CARB3     ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 152 DD: 71,WD: 53, CARB6  ~ MGLY
chch_t1(152,'CARB6     ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 153 DD: 72,WD: 54,
chch_t1(153,'CARB9     ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 154 DD: 73,WD: 55,    (210)
chch_t1(154,'CARB12    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 155 DD: 74,WD: 56,
chch_t1(155,'CARB15    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 156 DD: 75,
chch_t1(156,'UCARB12   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 157 DD: 76,WD: 57,
chch_t1(157,'NUCARB12  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 158 DD: 77,
chch_t1(158,'UDCARB8   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 159 DD: 78,
chch_t1(159,'UDCARB11  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 160 DD: 79,
chch_t1(160,'UDCARB14  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 161 DD: 80,
chch_t1(161,'TNCARB26  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 162 DD: 81,
chch_t1(162,'TNCARB10  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 163 RN10NO3 ~ n-PrONO2
chch_t1(163,'RN10NO3   ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 164                  (220)
chch_t1(164,'RN13NO3   ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 165
chch_t1(165,'RN16NO3   ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 166
chch_t1(166,'RN19NO3   ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 167
chch_t1(167,'RA13NO3   ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! WD + DD in CS2
chch_t1(167,'RA13NO3   ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 168
chch_t1(168,'RA16NO3   ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! WD + DD in CS2
chch_t1(168,'RA16NO3   ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 169
chch_t1(169,'RA19NO3   ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! WD + DD in CS2
chch_t1(169,'RA19NO3   ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 170 DD: 82,
chch_t1(170,'RTX24NO3  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 171 DD: 83,WD: 58, RN10OOH ~ n-PrOOH
chch_t1(171,'RN10OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 172 DD: 84,WD: 59,
chch_t1(172,'RN13OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 173 DD: 85,WD: 60,
chch_t1(173,'RN16OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 174 DD: 86,WD: 61,    (230)
chch_t1(174,'RN19OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 175 DD: 87,WD: 62, RN8OOH ~ ALKAOOH (from acetone degredation)
chch_t1(175,'RN8OOH    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 176 DD: 88,WD: 63, RN11OOH ~ MEKOOH
chch_t1(176,'RN11OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 177 DD: 89,WD: 64,
chch_t1(177,'RN14OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 178 DD: 90,WD: 65,
chch_t1(178,'RN17OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 179 DD: 91,WD: 66, RU14OOH ~ ISOOH
chch_t1(179,'RU14OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 180 DD: 92,WD: 67,
chch_t1(180,'RU12OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 181 DD: 93,WD: 68, RU10OOH ~ MACROOH
chch_t1(181,'RU10OOH   ',1,'TR        ','          ',1,1,cs,0,0,107)           &
]


chch_defs_master(245:n_chch_master)=[                                          &
! 182 DD: 94,WD: 69,
chch_t1(182,'NRU14OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 183 DD: 95,WD: 70,
chch_t1(183,'NRU12OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 184 DD: 96,WD: 71, RN9OOH from C3H6 degredation (240)
chch_t1(184,'RN9OOH    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 185 DD: 97,WD: 72, RN9OOH from TBUT2ENE degredation
chch_t1(185,'RN12OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 186 DD: 98,WD: 73,
chch_t1(186,'RN15OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 187 DD: 99,WD: 74,
chch_t1(187,'RN18OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 188 DD:100,WD: 75,
chch_t1(188,'NRN6OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 189 DD:101,WD: 76,
chch_t1(189,'NRN9OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 190 DD:102,WD: 77,
chch_t1(190,'NRN12OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 191 DD:103,WD: 78,
chch_t1(191,'RA13OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 192 DD:104,WD: 79,
chch_t1(192,'RA16OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 193 DD:105,WD: 80,
chch_t1(193,'RA19OOH   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 194 DD:106,WD: 81 ,   (250)
chch_t1(194,'RTN28OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 195 DD:107,WD: 82,
chch_t1(195,'NRTN28OOH ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 196 DD:108,WD: 83,
chch_t1(196,'RTN26OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 197 DD:109,WD: 84,
chch_t1(197,'RTN25OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 198 DD:110,WD: 85,
chch_t1(198,'RTN24OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 199 DD:111,WD: 86,
chch_t1(199,'RTN23OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 200 DD:112,WD: 87,
chch_t1(200,'RTN14OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 201 DD:113,WD: 88,
chch_t1(201,'RTN10OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 202 DD:114,WD: 89,
chch_t1(202,'RTX28OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 203 DD:115,WD: 90,
chch_t1(203,'RTX24OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 204 DD:116,WD: 91,    (260)
chch_t1(204,'RTX22OOH  ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 205 DD:117,WD: 92,
chch_t1(205,'NRTX28OOH ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 206 RN10O2 ~ n-PrOO
chch_t1(206,'RN10O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 207
chch_t1(207,'RN13O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 208
chch_t1(208,'RN16O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 209
chch_t1(209,'RN19O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 210
chch_t1(210,'RN13AO2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 211
chch_t1(211,'RN16AO2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 212
chch_t1(212,'RA13O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 213
chch_t1(213,'RA16O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 214                  (270)
chch_t1(214,'RA19AO2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 215
chch_t1(215,'RA19CO2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 216
chch_t1(216,'RN9O2     ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 217
chch_t1(217,'RN12O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 218
chch_t1(218,'RN15O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 219
chch_t1(219,'RN18O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 220
chch_t1(220,'RN15AO2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 221
chch_t1(221,'RN18AO2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 222
chch_t1(222,'RN8O2     ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 223 RN11O2 ~ MEKO2
chch_t1(223,'RN11O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 224                  (280)
chch_t1(224,'RN14O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 225
chch_t1(225,'RN17O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 226 RU14O2 ~ ISO2
chch_t1(226,'RU14O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 227
chch_t1(227,'RU12O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 228 RU14O2 ~ MACRO2
chch_t1(228,'RU10O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 229
chch_t1(229,'NRN6O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 230
chch_t1(230,'NRN9O2    ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 231
chch_t1(231,'NRN12O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 232
chch_t1(232,'NRU14O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 233
chch_t1(233,'NRU12O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 234 RTN products - from APINENE degredation (290)
chch_t1(234,'RTN28O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 235
chch_t1(235,'NRTN28O2  ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 236
chch_t1(236,'RTN26O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 237
chch_t1(237,'RTN25O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 238
chch_t1(238,'RTN24O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 239
chch_t1(239,'RTN23O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 240
chch_t1(240,'RTN14O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 241
chch_t1(241,'RTN10O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 242; RTX products - from BPINENE degredation
chch_t1(242,'RTX28O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 243
chch_t1(243,'NRTX28O2  ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 244                  (300)
chch_t1(244,'RTX24O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 245
chch_t1(245,'RTX22O2   ',1,'OO        ','          ',0,0,cs,rp,0,107),         &
! 246
chch_t1(246,'RAROH14   ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 247
chch_t1(247,'RAROH17   ',1,'TR        ','          ',0,0,cs,0,0,107),          &
! 248 DD:118, Species not used in CS2
chch_t1(248,'RU12PAN   ',1,'TR        ','          ',1,0,cs,0,cs2,107),        &
! 249 DD:119,
chch_t1(249,'RTN26PAN  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! WD + DD in CS2
chch_t1(249,'RTN26PAN  ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 250 DD:120, Species not used in CS2
chch_t1(250,'TNCARB12  ',1,'TR        ','          ',1,0,cs,0,cs2,107),        &
! 251 DD:121, Species not used in CS2
chch_t1(251,'TNCARB11  ',1,'TR        ','          ',1,0,cs,0,cs2,107),        &
! 252 DD:122, Species not used in CS2
chch_t1(252,'RTN23NO3  ',1,'TR        ','          ',1,0,cs,0,cs2,107),        &
! 253 DD:123,
chch_t1(253,'CCARB12   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 254 DD:124,          (310)
chch_t1(254,'TNCARB15  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 255 DD:125,
chch_t1(255,'RCOOH25   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 256 DD:126,
chch_t1(256,'TXCARB24  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 257 DD:127,
chch_t1(257,'TXCARB22  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 258 DD:128,
chch_t1(258,'RN9NO3    ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! DD and WD in CS2
chch_t1(258,'RN9NO3    ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 259 DD:129,
chch_t1(259,'RN12NO3   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! DD and WD in CS2
chch_t1(259,'RN12NO3   ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 260 DD:130,
chch_t1(260,'RN15NO3   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! WD and DD in CS2
chch_t1(260,'RN15NO3   ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 261 DD:131,
chch_t1(261,'RN18NO3   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! DD and WD in CS2
chch_t1(261,'RN18NO3   ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 262 DD:132,
chch_t1(262,'RU14NO3   ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! WD + DD in CS2
chch_t1(262,'RU14NO3   ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 263 DD:133,
chch_t1(263,'RTN28NO3  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! WD + DD in CS2
chch_t1(263,'RTN28NO3  ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 264 DD:134,          (320)
chch_t1(264,'RTN25NO3  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! WD + DD in CS2
chch_t1(264,'RTN25NO3  ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 265 DD:135,
chch_t1(265,'RTX28NO3  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! WD + DD in CS2
chch_t1(265,'RTX28NO3  ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 266 DD:136,
chch_t1(266,'RTX22NO3  ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! WD + DD in CS2
chch_t1(266,'RTX22NO3  ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 267 DD:137,WD: 93, AROH14 ~ PHENOL;  from BENZENE oxidation
chch_t1(267,'AROH14    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 268 DD:138,WD: 94,
chch_t1(268,'ARNOH14   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 269 DD:139,WD: 95, AROH17 ~ CRESOL;  from TOLUENE oxidation
chch_t1(269,'AROH17    ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 270 DD:140,WD: 96,
chch_t1(270,'ARNOH17   ',1,'TR        ','          ',1,1,cs,0,0,107),          &
! 271 DD:141, ANHY ~ MALANHY (327)
chch_t1(271,'ANHY      ',1,'TR        ','          ',1,0,cs,0,0,107),          &
! 272
chch_t1(272,'IEPOX     ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 273
chch_t1(273,'HMML      ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 274
chch_t1(274,'HUCARB9   ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 275
chch_t1(275,'HPUCARB12 ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 276
chch_t1(276,'DHPCARB9  ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 277
chch_t1(277,'DHPR12O2  ',1,'OO        ','          ',0,0,cs,rp,0,119),         &
! 278
chch_t1(278,'DHPR12OOH ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 279
chch_t1(279,'RU10AO2   ',1,'OO        ','          ',0,0,cs,rp,0,119),         &
! 280
chch_t1(280,'DHCARB9   ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 281
chch_t1(281,'RU12NO3   ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 282
chch_t1(282,'RU10NO3   ',1,'TR        ','          ',1,1,cs,0,0,119),          &
! 283
chch_t1(283,'MACO3     ',1,'OO        ','          ',0,0,cs,rp,0,119)          &
]



! ----------------------------------------------------------------------


! Heterogeneous chemistry
! Columns take the following meanings:
! Item number, reactant1, reactant2, product1, product2, product3, product4,
! product yield x 4, chemistry scheme, qualifier, disqualifier, version
rath_defs_master(1:n_het_master)=[                                             &
rath_t1(1,'ClONO2    ','H2O       ','HOCl      ','HONO2     ','          ',    &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,hp,0,107),                   &
rath_t1(2,'ClONO2    ','HCl       ','Cl        ','Cl        ','HONO2     ',    &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,hp,0,107),                   &
rath_t1(3,'HOCl      ','HCl       ','Cl        ','Cl        ','H2O       ',    &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,hp,0,107),                   &
rath_t1(4,'N2O5      ','H2O       ','HONO2     ','HONO2     ','          ',    &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,hp,0,107),                   &
rath_t1(5,'N2O5      ','HCl       ','Cl        ','NO2       ','HONO2     ',    &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,hp,0,107),                   &
rath_t1(6,'ClONO2    ','HBr       ','BrCl      ','HONO2     ','          ',    &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,eh,0,111),                   &
rath_t1(7,'HOCl      ','HBr       ','BrCl      ','H2O       ','          ',    &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,eh,0,111),                   &
rath_t1(8,'HOBr      ','HCl       ','BrCl      ','H2O       ','          ',    &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,eh,0,111),                   &
rath_t1(9,'BrONO2    ','HCl       ','BrCl      ','HONO2     ','          ',    &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,eh,0,111),                   &
rath_t1(10,'BrONO2    ','H2O       ','HOBr      ','HONO2     ','          ',   &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,eh,0,111),                   &
rath_t1(11,'HOBr      ','HBr       ','Br        ','Br        ','H2O       ',   &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,eh,0,111),                   &
rath_t1(12,'BrONO2    ','HBr       ','Br        ','Br        ','HONO2     ',   &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,eh,0,111),                   &
rath_t1(13,'N2O5      ','HBr       ','Br        ','NO2       ','HONO2     ',   &
'          ', 0.000, 0.000, 0.000, 0.000, s+st+cs,eh,0,111),                   &

! Aerosol chemistry: there are no gas phase products, the 'NULLx' products
!  identify the reactions in asad_hetero

!HSO3+H2O2(aq)
rath_t1(14,'SO2       ','H2O2      ','NULL0     ','          ','          ',   &
'          ', 0.000, 0.000, 0.000, 0.000, ti+s+st+ol+r+cs,a,0,107),            &
!HSO3+O3(aq)
rath_t1(15,'SO2       ','O3        ','NULL1     ','          ','          ',   &
'          ', 0.000, 0.000, 0.000, 0.000, ti+s+st+ol+r+cs,a,0,107),            &
!HSO3+O3(aq)
rath_t1(16,'SO2       ','O3        ','NULL2     ','          ','          ',   &
'          ', 0.000, 0.000, 0.000, 0.000, ti+s+st+ol+r+cs,a,0,107),            &

! Tropospheric heterogeneous reactions
rath_t1(17,'N2O5      ','          ','HONO2     ','          ','          ',   &
'          ', 2.000, 0.000, 0.000, 0.000, ti+st+r+cs,th,0,107),                &
! Heterogeneous
rath_t1(18,'HO2       ','          ','H2O2      ','          ','          ',   &
'          ', 0.500, 0.000, 0.000, 0.000, ti+st,th,0,107)  ]

! ----------------------------------------------------------------------


! Photolysis reactions
! Item number, reactant, 'PHOTON', product x4,
! fractional production coefficient x4, quantum yield (%),
! photolysis rate label, chemistry scheme, qualifier, disqualifier, version
ratj_defs_master(1:100)=[                                                      &
! 1
ratj_t1(1,'EtOOH     ','PHOTON    ','MeCHO     ','HO2       ','OH        ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',ti+t+st+r+cs,0,0,107),&
! 2
ratj_t1(2,'H2O2      ','PHOTON    ','OH        ','OH        ','          ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.0,'jh2o2     ',ti+s+t+st+r+cs,0,0,107),&
! 3
! This should produce H+ CHO -> H + HO2 + CO in ST scheme.
ratj_t1(3,'HCHO      ','PHOTON    ','HO2       ','HO2       ','CO        ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhchoa    ',ti+t+st+r+cs,0,0,107),&
! CS2
ratj_t1(3,'HCHO      ','PHOTON    ','HO2       ','HO2       ','CO        ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jforma    ',cs,0,0,119),          &
ratj_t1(3,'HCHO      ','PHOTON    ','H         ','CO        ','HO2       ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhchoa    ',s,0,0,107) ,          &
! 4
ratj_t1(4,'HCHO      ','PHOTON    ','H2        ','CO        ','          ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.0,'jhchob    ',ti+s+t+st+r+cs,0,0,107),&
! CS2
ratj_t1(4,'HCHO      ','PHOTON    ','H2        ','CO        ','          ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.0,'jformb    ',cs,0,0,119),            &
! 5
! In S and ST, branching is assumed following JPL 2011. Update T, TI, and R?
! Added photolysis of HO2NO2 for CRI-Strat.
ratj_t1(5,'HO2NO2    ','PHOTON    ','HO2       ','NO2       ','          ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna67    ',st+cs,0,0,107) ,      &
ratj_t1(5,'HO2NO2    ','PHOTON    ','NO2       ','HO2       ','          ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna67    ',s,0,0,107) ,          &
ratj_t1(5,'HO2NO2    ','PHOTON    ','HO2       ','NO2       ','          ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna      ',ti+t+r,0,0,107) ,     &
! 6
ratj_t1(6,'HONO2     ','PHOTON    ','OH        ','NO2       ','          ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhono2    ',ti+t+st+r+cs,0,0,107),&
ratj_t1(6,'HONO2     ','PHOTON    ','NO2       ','OH        ','          ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhono2    ',s,0,0,107) ,          &
! 7
ratj_t1(7,'MeCHO     ','PHOTON    ','MeOO      ','HO2       ','CO        ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jaceta    ',ti+t+st+r+cs,0,0,107),&
! 8
! Missing in R (And CRI, Added for CRI-Strat).
ratj_t1(8,'MeCHO     ','PHOTON    ','CH4       ','CO        ','          ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jacetb    ',ti+t+st+cs,0,0,107) , &
! 9
! THIS REACTION IS MISSING IN THE TROPOSPHERIC CHEMISTRY SCHEME.
ratj_t1(9,'MeOOH     ','PHOTON    ','HO2       ','HCHO      ','OH        ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',ti+st+r+cs,0,0,107) , &
ratj_t1(9,'MeOOH     ','PHOTON    ','HCHO      ','HO2       ','OH        ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',s,0,0,107) ,          &
! 10
ratj_t1(10,'N2O5      ','PHOTON    ','NO3       ','NO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jn2o5     ',ti+t+st+r+cs,0,0,107),&
ratj_t1(10,'N2O5      ','PHOTON    ','NO2       ','NO3       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jn2o5     ',s,0,0,107) ,          &
! 11
ratj_t1(11,'NO2       ','PHOTON    ','NO        ','O(3P)     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.0,'jno2      ',ti+s+t+st+r+cs,0,0,107),&
! 12
ratj_t1(12,'NO3       ','PHOTON    ','NO        ','O2        ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.0,'jno3a     ',ti+s+t+st+r+cs,0,0,107),&
! 13
ratj_t1(13,'NO3       ','PHOTON    ','NO2       ','O(3P)     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.0,'jno3b     ',ti+s+t+st+r+cs,0,0,107),&
! 14
! Missing in R, Added to CRI-Strat
ratj_t1(14,'O2        ','PHOTON    ','O(3P)     ','O(3P)     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.0,'jo2       ',ti+s+t+st+cs,0,0,107),  &
! 15
ratj_t1(15,'O3        ','PHOTON    ','O2        ','O(1D)     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.0,'jo3a      ',ti+s+t+st+r+cs,0,0,107),&
! 16
ratj_t1(16,'O3        ','PHOTON    ','O2        ','O(3P)     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.0,'jo3b      ',ti+s+t+st+r+cs,0,0,107),&
! 17
ratj_t1(17,'PAN       ','PHOTON    ','MeCO3     ','NO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.0,'jpan      ',ti+t+st+r+cs,0,0,107) , &
! 18
! HONO not modelled in R
ratj_t1(18,'HONO      ','PHOTON    ','OH        ','NO        ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhono     ',ti+t+st+cs,0,0,107) , &
! 19
ratj_t1(19,'EtCHO     ','PHOTON    ','EtOO      ','HO2       ','CO        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jetcho    ',ti+t+st+cs,0,0,107) , &
! CS2
ratj_t1(19,'EtCHO     ','PHOTON    ','EtOO      ','HO2       ','CO        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jprpal    ',cs,0,0,119) ,         &
! 20
ratj_t1(20,'Me2CO     ','PHOTON    ','MeCO3     ','MeOO      ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jaceto    ',ti+st+r+cs,0,0,107) , &
! The naming of 'jacetone' is probably inconsistent with FAST-JX. Error?
ratj_t1(20,'Me2CO     ','PHOTON    ','MeCO3     ','MeOO      ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jacetone  ',t,0,0,107) ,          &
! 21 missing in A
ratj_t1(21,'n-PrOOH   ','PHOTON    ','EtCHO     ','HO2       ','OH        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',t+st,0,0,107) ,       &
ratj_t1(21,'n-PrOOH   ','PHOTON    ','HO2       ','EtCHO     ','OH        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',ti,0,0,107) ,         &
! 22 missing in A
ratj_t1(22,'i-PrOOH   ','PHOTON    ','Me2CO     ','HO2       ','OH        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',t+st+r+cs,0,0,107) ,  &
ratj_t1(22,'i-PrOOH   ','PHOTON    ','HO2       ','Me2CO     ','OH        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',ti,0,0,107) ,         &
! 23
! reaction missing in tropospheric chemistry scheme
ratj_t1(23,'MeCOCH2OOH','PHOTON    ','MeCO3     ','HCHO      ','OH        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',st,0,0,107) ,         &
ratj_t1(23,'MeCOCH2OOH','PHOTON    ','HCHO      ','MeCO3     ','OH        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',ti,0,0,107) ,         &
! 24 - Note no CRI photolysis of PPAN-> Added for FT
ratj_t1(24,'PPAN      ','PHOTON    ','EtCO3     ','NO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpan      ',ti+t+st+cs,0,0,107) , &
! 25
ratj_t1(25,'MeONO2    ','PHOTON    ','HO2       ','HCHO      ','NO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmena     ',ti+t+st+cs,0,0,107) , &
! 26 Different products in TI and ST versus R.
ratj_t1(26,'ISOOH     ','PHOTON    ','OH        ','MACR      ','HCHO      ',   &
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',ti+st,0,0,107) ,      &
ratj_t1(26,'ISOOH     ','PHOTON    ','OH        ','MVK      ','HCHO      ',    &
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',r,0,0,107) ,          &
! 27
ratj_t1(27,'ISON      ','PHOTON    ','NO2       ','MACR      ','HCHO      ',   &
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jiprn     ',ti+st,0,0,107) ,      &
! 28
ratj_t1(28,'MACR      ','PHOTON    ','MeCO3     ','HCHO      ','CO        ',   &
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jmacr     ',ti+st,0,0,107) ,      &
! 29
ratj_t1(29,'MPAN      ','PHOTON    ','MACRO2    ','NO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpan      ',ti+st,0,0,107) ,      &
! 30
ratj_t1(30,'MACROOH   ','PHOTON    ','OH        ','HO2       ','OH        ',   &
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jmacro    ',ti+st,0,0,107) ,      &
! 31
ratj_t1(31,'MACROOH   ','PHOTON    ','HACET     ','CO        ','MGLY      ',   &
     'HCHO      ', 0.0,0.0,0.0,0.0, 100.000,'jmacro    ',ti+st,0,0,107) ,      &
! 32
ratj_t1(32,'HACET     ','PHOTON    ','MeCO3     ','HCHO      ','HO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhacet    ',ti+st,0,0,107) ,      &
! 33
ratj_t1(33,'MGLY      ','PHOTON    ','MeCO3     ','CO        ','HO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmkal     ',ti+st+r,0,0,107) ,    &
! 34
ratj_t1(34,'NALD      ','PHOTON    ','HCHO      ','CO        ','NO2       ',   &
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jaceta    ',ti+st,0,0,107) ,      &
! 35
ratj_t1(35,'MeCO3H    ','PHOTON    ','MeOO      ','OH        ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmeco3h   ',ti+st,0,0,107) ,      &
! 35 - Same as MeOOH in CRI
ratj_t1(35,'MeCO3H    ','PHOTON    ','MeOO      ','OH        ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 36
ratj_t1(36,'BrCl      ','PHOTON    ','Br        ','Cl        ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jbrcl     ',s+st+cs,0,0,107) ,    &
! 37
ratj_t1(37,'BrO       ','PHOTON    ','Br        ','O(3P)     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jbro      ',s+st+cs,0,0,107) ,    &
! 38
ratj_t1(38,'BrONO2    ','PHOTON    ','Br        ','NO3       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jbrna     ',s+st+cs,0,0,107) ,    &
! 39
ratj_t1(39,'BrONO2    ','PHOTON    ','BrO       ','NO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jbrnb     ',s+st+cs,0,0,107) ,    &
! 40
ratj_t1(40,'O2        ','PHOTON    ','O(3P)     ','O(1D)     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jo2b      ',s+st+cs,0,0,107) ,    &
! 41
ratj_t1(41,'OClO      ','PHOTON    ','O(3P)     ','ClO       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'joclo     ',s+st+cs,0,0,107) ,    &
! 42
ratj_t1(42,'NO        ','PHOTON    ','N         ','O(3P)     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jno       ',s+st+cs,0,0,107) ,    &
! 43
ratj_t1(43,'HOBr      ','PHOTON    ','OH        ','Br        ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhobr     ',s+st+cs,0,0,107) ,    &
! 44
ratj_t1(44,'N2O       ','PHOTON    ','N2        ','O(1D)     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jn2o      ',s+st+cs,0,0,107) ,    &
! 45
ratj_t1(45,'H2O       ','PHOTON    ','OH        ','H         ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jh2o      ',s+st+cs,0,0,107) ,    &
! 46
ratj_t1(46,'ClONO2    ','PHOTON    ','Cl        ','NO3       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jclna     ',s+st+cs,0,0,107) ,    &
! 47
ratj_t1(47,'ClONO2    ','PHOTON    ','ClO       ','NO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jclnb     ',s+st+cs,0,0,107) ,    &
! 48
ratj_t1(48,'HCl       ','PHOTON    ','H         ','Cl        ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhcl      ',s+st+cs,0,0,107) ,    &
! 49
ratj_t1(49,'HOCl      ','PHOTON    ','OH        ','Cl        ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhocl     ',s+st+cs,0,0,107) ,    &
! 50
ratj_t1(50,'Cl2O2     ','PHOTON    ','Cl        ','Cl        ','O2        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jcl2o2    ',s+st+cs,0,0,107) ,    &
! 51
ratj_t1(51,'CFCl3     ','PHOTON    ','Cl        ','Cl        ','Cl        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jcfcl3    ',s+st+cs,0,0,107) ,    &
! 52
ratj_t1(52,'CF2Cl2    ','PHOTON    ','Cl        ','Cl        ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jcfc2     ',s+st+cs,0,0,107) ,    &
! 53
ratj_t1(53,'MeBr      ','PHOTON    ','Br        ','H         ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmebr     ',s+st+cs,0,0,107) ,    &
! 54
ratj_t1(54,'CH4       ','PHOTON    ','MeOO      ','H         ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jch4      ',s+st+cs,0,0,107) ,    &
! 55
ratj_t1(55,'CO2       ','PHOTON    ','CO        ','O(3P)     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jco2      ',s+st+cs,0,0,107),     &
! 56
! Reaction missing in the R scheme.
ratj_t1(56,'HO2NO2    ','PHOTON    ','OH        ','NO3       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna33    ',st+cs,0,0,107),       &
ratj_t1(56,'HO2NO2    ','PHOTON    ','NO3       ','OH       ','          ',    &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna33    ',s,0,0,107),           &
! 57
ratj_t1(57,'CS2       ','PHOTON    ','COS       ','SO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jcs2      ',s+st+cs,a,0,107),     &
! 58
ratj_t1(58,'COS       ','PHOTON    ','CO        ','SO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jcos      ',s+st+cs,a,0,107),     &
! 59
ratj_t1(59,'H2SO4     ','PHOTON    ','SO3       ','OH        ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jh2so4    ',s+st+cs,a,0,107),     &
! update for 121
ratj_t1(59,'H2SO4     ','PHOTON    ','SO3       ','OH        ','H         ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jh2so4    ',s+st+cs,a,0,121),     &
! 60
ratj_t1(60,'SO3       ','PHOTON    ','SO2       ','O(3P)     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jso3      ',s+st+cs,a,0,107),     &
! 61
ratj_t1(61,'MEK       ','PHOTON    ','MeCO3     ','EtOO      ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jaceto    ',r,0,0,107),           &
! 61 - approximate rate relative to jno2 in CRI
ratj_t1(61,'MEK       ','PHOTON    ','MeCO3     ','EtOO      ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 5.497e-02,'jno2      ',cs,0,0,107),        &
! 62
ratj_t1(62,'GLY       ','PHOTON    ','CO        ','CO        ','HO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmkal     ',r,0,0,107),           &
! 63
ratj_t1(63,'s-BuOOH   ','PHOTON    ','MEK       ','HO2       ','OH        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',r,0,0,107),           &
! 64
ratj_t1(64,'MVKOOH    ','PHOTON    ','MGLY      ','HCHO      ','OH        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',r,0,0,107),           &
! From here, photolysis reactions for CRI. Approximate photolysis rates
! used where equivalent cross-sections in UKCA not available.
! See http://mcm.york.ac.uk/parameters/photolysis.htt for details
! 65
ratj_t1(65,'CARB14    ','PHOTON    ','MeCO3     ','RN10O2    ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 2.605e-01 ,'jno2      ',cs,0,0,107),       &
! 66
ratj_t1(66,'CARB17    ','PHOTON    ','RN8O2     ','RN10O2    ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 7.311e-02 ,'jno2      ',cs,0,0,107),       &
! 67
ratj_t1(67,'CARB11A   ','PHOTON    ','MeCO3     ','EtOO      ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 5.497e-02 ,'jno2      ',cs,0,0,107),       &
! 68
ratj_t1(68,'CARB7     ','PHOTON    ','MeCO3     ','HCHO      ','HO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 5.497e-02 ,'jno2      ',cs,0,0,107),       &
! 69
ratj_t1(69,'CARB10    ','PHOTON    ','MeCO3     ','MeCHO     ','HO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 5.497e-02 ,'jno2      ',cs,0,0,107),       &
! 70
ratj_t1(70,'CARB13    ','PHOTON    ','RN8O2     ','MeCHO     ','HO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 1.649e-01 ,'jno2      ',cs,0,0,107),       &
! 71
ratj_t1(71,'CARB16    ','PHOTON    ','RN8O2     ','EtCHO     ','HO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 1.841e-01 ,'jno2      ',cs,0,0,107),       &
! 72
ratj_t1(72,'HOCH2CHO  ','PHOTON    ','HCHO      ','CO        ','HO2       ',   &
     'HO2       ', 0.0,0.0,0.0,0.0, 9.640e+02 ,'jhchoa    ',cs,0,0,107),       &
! 73
ratj_t1(73,'UCARB10   ','PHOTON    ','MeCO3     ','HCHO      ','HO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 200.000,'jmacr     ',cs,0,0,107),          &
! CS2
ratj_t1(73,'UCARB10   ','PHOTON    ','MeCO3     ','HCHO      ','HO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 200.000,'jmacr2    ',cs,0,0,119),          &
! 74
ratj_t1(74,'CARB3     ','PHOTON    ','CO        ','CO        ','HO2       ',   &
     'HO2       ', 0.0,0.0,0.0,0.0, 9.640e+02 ,'jhchoa    ',cs,0,0,107),       &
! 74 CS2
ratj_t1(74,'CARB3     ','PHOTON    ','CO        ','CO        ','HO2       ',   &
     'HO2       ', 0.0,0.0,0.0,0.0, 100.00 ,'jcb3c     ',cs,0,0,119),          &
! 75
ratj_t1(75,'CARB6     ','PHOTON    ','MeCO3     ','CO        ','HO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.00 ,'jmkal     ',cs,0,0,107),          &
! 76
ratj_t1(76,'CARB9     ','PHOTON    ','MeCO3     ','MeCO3     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 3.224e+00 ,'jno2      ',cs,0,0,107),       &
! 77
ratj_t1(77,'CARB12    ','PHOTON    ','MeCO3     ','RN8O2     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 3.224e+00 ,'jno2      ',cs,0,0,107),       &
! 78
ratj_t1(78,'CARB15    ','PHOTON    ','RN8O2     ','RN8O2     ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 3.224e+00 ,'jno2      ',cs,0,0,107),       &
! 79
ratj_t1(79,'UCARB12   ','PHOTON    ','MeCO3     ','HOCH2CHO  ','CO        ',   &
     'HO2       ', 0.0,0.0,0.0,0.0, 200.000,'jmacr     ',cs,0,0,107),          &
! CS2
ratj_t1(79,'UCARB12   ','PHOTON    ','MeCO3     ','HOCH2CHO  ','CO        ',   &
     'HO2       ', 0.0,0.0,0.0,0.0, 25.000,'jmacr2    ',cs,0,0,119)            &
]


! Photolysis reactions continued
ratj_defs_master(101:n_ratj_master)=[                                          &
! 80
ratj_t1(80,'NUCARB12  ','PHOTON    ','NOA       ','CO        ','HO2       ',   &
     '          ', 1.0,2.0,2.0,0.0, 100.000,'jmacr     ',cs,0,0,107),          &
! 80
ratj_t1(80,'NUCARB12  ','PHOTON    ','HUCARB9   ','CO        ','NO2       ',   &
     'OH        ', 0.0,0.0,0.0,0.0, 0.6066, 'jno2     ',cs,0,0,119),           &
! 81
ratj_t1(81,'NOA       ','PHOTON    ','MeCO3     ','HCHO      ','NO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 7.582e-02 ,'jno2      ',cs,0,0,107),       &
! 82
ratj_t1(82,'UDCARB8   ','PHOTON    ','EtOO      ','HO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 1.280e+00 ,'jno2      ',cs,0,0,107),       &
! 83
ratj_t1(83,'UDCARB11  ','PHOTON    ','RN10O2    ','HO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 1.100e+00 ,'jno2      ',cs,0,0,107),       &
! 84
ratj_t1(84,'UDCARB14  ','PHOTON    ','RN13O2    ','HO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 1.100e+00 ,'jno2      ',cs,0,0,107),       &
! 85
ratj_t1(85,'TNCARB26  ','PHOTON    ','RTN26O2   ','HO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 9.640e+02 ,'jhchoa    ',cs,0,0,107),       &
! 86
ratj_t1(86,'TNCARB10  ','PHOTON    ','MeCO3     ','MeCO3     ','CO        ',   &
     '          ', 0.0,0.0,0.0,0.0, 1.612e+00 ,'jno2      ',cs,0,0,107),       &
! 87
ratj_t1(87,'EtONO2    ','PHOTON    ','MeCHO     ','HO2       ','NO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 1.813e-02 ,'jno2      ',cs,0,0,107),       &
! 88
ratj_t1(88,'RN10NO3   ','PHOTON    ','EtCHO     ','HO2       ','NO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 2.402e-02 ,'jno2      ',cs,0,0,107),       &
! 89
ratj_t1(89,'i-PrONO2  ','PHOTON    ','Me2CO     ','HO2       ','NO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 4.078e-02 ,'jno2      ',cs,0,0,107),       &
! 90
ratj_t1(90,'RN13NO3   ','PHOTON    ','MeCHO     ','EtOO      ','NO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 9.558e-03 ,'jno2      ',cs,0,0,107),       &
! 91
ratj_t1(91,'RN13NO3   ','PHOTON    ','CARB11A   ','HO2       ','NO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 1.446e-02 ,'jno2      ',cs,0,0,107),       &
! 92
ratj_t1(92,'RN16NO3   ','PHOTON    ','RN15O2    ','NO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 2.402e-02 ,'jno2      ',cs,0,0,107),       &
! 93
ratj_t1(93,'RN19NO3   ','PHOTON    ','RN18O2    ','NO2       ','          ',   &
     '          ', 0.0,0.0,0.0,0.0, 2.402e-02 ,'jno2      ',cs,0,0,107),       &
! 94
ratj_t1(94,'RA13NO3   ','PHOTON    ','CARB3     ','UDCARB8   ','HO2       ',   &
     'NO2       ', 0.0,0.0,0.0,0.0, 4.078e-02 ,'jno2      ',cs,0,0,107),       &
! 95
ratj_t1(95,'RA16NO3   ','PHOTON    ','CARB3     ','UDCARB11  ','HO2       ',   &
     'NO2       ', 0.0,0.0,0.0,0.0, 4.078e-02 ,'jno2      ',cs,0,0,107),       &
! 96
ratj_t1(96,'RA19NO3   ','PHOTON    ','CARB6     ','UDCARB11  ','HO2       ',   &
     'NO2       ', 0.0,0.0,0.0,0.0, 4.078e-02 ,'jno2      ',cs,0,0,107),       &
! 97
ratj_t1(97,'RTX24NO3  ','PHOTON    ','TXCARB22  ','HO2       ','NO2       ',   &
     '          ', 0.0,0.0,0.0,0.0, 4.078e-02 ,'jno2      ',cs,0,0,107),       &
! 98
ratj_t1(98,'RN10OOH   ','PHOTON    ','EtCHO     ','HO2       ','OH        ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 99
ratj_t1(99,'RN13OOH   ','PHOTON    ','MeCHO     ','EtOO      ','OH        ',   &
     '          ', 0.0,0.0,0.0,0.0,  39.800,'jmhp      ',cs,0,0,107),          &
! 100
ratj_t1(100,'RN13OOH   ','PHOTON    ','CARB11A   ','HO2       ','OH        ',  &
     '          ', 0.0,0.0,0.0,0.0,  60.200,'jmhp      ',cs,0,0,107),          &
! 101
ratj_t1(101,'RN16OOH   ','PHOTON    ','RN15AO2   ','OH        ','          ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 102
ratj_t1(102,'RN19OOH   ','PHOTON    ','RN18AO2   ','OH        ','          ',  &
     '          ', 0.0,0.0,0.0,0.0,  100.000,'jmhp      ',cs,0,0,107),         &
! 103
ratj_t1(103,'EtCO3H    ','PHOTON    ','EtOO      ','OH        ','          ',  &
     '          ', 0.0,0.0,0.0,0.0,  100.000,'jmhp      ',cs,0,0,107),         &
! 104
ratj_t1(104,'HOCH2CO3H ','PHOTON    ','HCHO      ','HO2       ','OH        ',  &
     '          ', 0.0,0.0,0.0,0.0,  100.000,'jmhp      ',cs,0,0,107),         &
! 105
ratj_t1(105,'RN8OOH    ','PHOTON    ','EtOO      ','OH        ','          ',  &
     '          ', 0.0,0.0,0.0,0.0,  100.000,'jmhp      ',cs,0,0,107),         &
! 106
ratj_t1(106,'RN11OOH   ','PHOTON    ','RN10O2    ','OH        ','          ',  &
     '          ', 0.0,0.0,0.0,0.0,  100.000,'jmhp      ',cs,0,0,107),         &
! 107
ratj_t1(107,'RN14OOH   ','PHOTON    ','RN13O2    ','OH        ','          ',  &
     '          ', 0.0,0.0,0.0,0.0,  100.000,'jmhp      ',cs,0,0,107),         &
! 108
ratj_t1(108,'RN17OOH   ','PHOTON    ','RN16O2    ','OH        ','          ',  &
     '          ', 0.0,0.0,0.0,0.0,  100.000,'jmhp      ',cs,0,0,107),         &
! 109
ratj_t1(109,'RU14OOH   ','PHOTON    ','UCARB12   ','HO2       ','OH        ',  &
     '          ', 0.0,0.0,0.0,0.0,  25.200,'jmhp      ',cs,0,0,107),          &
! 110
ratj_t1(110,'RU14OOH   ','PHOTON    ','UCARB10   ','HCHO      ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0,  74.800,'jmhp      ',cs,0,0,107),          &
! 111
ratj_t1(111,'RU12OOH   ','PHOTON    ','CARB6     ','HOCH2CHO  ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 112
ratj_t1(112,'RU10OOH   ','PHOTON    ','MeCO3     ','HOCH2CHO  ','OH        ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 113
ratj_t1(113,'NRU14OOH  ','PHOTON    ','NUCARB12  ','HO2       ','OH        ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 114
ratj_t1(114,'NRU12OOH  ','PHOTON    ','NOA       ','CO        ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 114
ratj_t1(114,'NRU12OOH  ','PHOTON    ','NOA       ','CARB3     ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,119),          &
! 115
ratj_t1(115,'HOC2H4OOH ','PHOTON    ','HCHO      ','HCHO      ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 116
ratj_t1(116,'RN9OOH    ','PHOTON    ','MeCHO     ','HCHO      ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 117
ratj_t1(117,'RN12OOH   ','PHOTON    ','MeCHO     ','MeCHO     ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 118
ratj_t1(118,'RN15OOH   ','PHOTON    ','EtCHO     ','MeCHO     ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 119
ratj_t1(119,'RN18OOH   ','PHOTON    ','EtCHO     ','EtCHO     ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 120
ratj_t1(120,'NRN6OOH   ','PHOTON    ','HCHO      ','HCHO      ','NO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 121
ratj_t1(121,'NRN9OOH   ','PHOTON    ','MeCHO     ','HCHO      ','NO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 122
ratj_t1(122,'NRN12OOH  ','PHOTON    ','MeCHO     ','MeCHO     ','NO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 123
ratj_t1(123,'RA13OOH   ','PHOTON    ','CARB3     ','UDCARB8   ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 124
ratj_t1(124,'RA16OOH   ','PHOTON    ','CARB3     ','UDCARB11  ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 125
ratj_t1(125,'RA19OOH   ','PHOTON    ','CARB6     ','UDCARB11  ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 126
ratj_t1(126,'RTN28OOH  ','PHOTON    ','TNCARB26  ','HO2       ','OH        ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 127
ratj_t1(127,'NRTN28OOH ','PHOTON    ','TNCARB26  ','NO2       ','OH        ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 128
ratj_t1(128,'RTN26OOH  ','PHOTON    ','RTN25O2   ','OH        ','          ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 129
ratj_t1(129,'RTN25OOH  ','PHOTON    ','RTN24O2   ','OH        ','          ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 130
ratj_t1(130,'RTN24OOH  ','PHOTON    ','RTN23O2   ','OH        ','          ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 131
ratj_t1(131,'RTN23OOH  ','PHOTON    ','Me2CO     ','RTN14O2   ','OH        ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 132
ratj_t1(132,'RTN14OOH  ','PHOTON    ','TNCARB10  ','HCHO      ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 133
ratj_t1(133,'RTN10OOH  ','PHOTON    ','RN8O2     ','CO        ','OH        ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 134
ratj_t1(134,'RTX28OOH  ','PHOTON    ','TXCARB24  ','HCHO      ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 135
ratj_t1(135,'RTX24OOH  ','PHOTON    ','TXCARB22  ','HO2       ','OH        ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 136
ratj_t1(136,'RTX22OOH  ','PHOTON    ','Me2CO     ','RN13O2    ','OH        ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 137
ratj_t1(137,'NRTX28OOH ','PHOTON    ','TXCARB24  ','HCHO      ','NO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',cs,0,0,107),          &
! 138
ratj_t1(138,'UDCARB8   ','PHOTON    ','ANHY      ','HO2       ','HO2       ',  &
     '          ', 0.0,0.0,0.0,0.0, 7.200e-01 ,'jno2      ',cs,0,0,107),       &
! 139
ratj_t1(139,'UDCARB11  ','PHOTON    ','ANHY      ','HO2       ','MeOO      ',  &
     '          ', 0.0,0.0,0.0,0.0, 9.000e-01 ,'jno2      ',cs,0,0,107),       &
! 140
ratj_t1(140,'UDCARB14  ','PHOTON    ','ANHY      ','HO2       ','EtOO      ',  &
     '          ', 0.0,0.0,0.0,0.0, 9.000e-01 ,'jno2      ',cs,0,0,107),       &
! 141
ratj_t1(141,'TNCARB12  ','PHOTON    ','RN9O2     ','HOCH2CO3  ','          ',  &
     '          ', 0.0,0.0,0.0,0.0, 5.497e-02 ,'jno2      ',cs,0,0,107),       &
! 141 CS2
ratj_t1(141,'CARB3     ','PHOTON    ','HCHO      ','CO        ','         ',   &
     '          ', 0.0,0.0,0.0,0.0, 100.00,'jcb3a     ',cs,0,0,119),           &
! 142
ratj_t1(142,'TNCARB11  ','PHOTON    ','RTN10O2   ','CO        ','HO2       ',  &
     '          ', 0.0,0.0,0.0,0.0, 1.500e+00 ,'jno2      ',cs,0,0,107),       &
! 142
ratj_t1(142,'RU14NO3   ','PHOTON    ','UCARB10   ','HCHO      ','HO2       ',  &
     'NO2       ', 0.0,0.0,0.0,0.0, 100.000,'jiprn     ',cs,0,0,119),          &
! Added MeO2NO2 photolysis for accurate UT NOx
! 143
ratj_t1(143,'MeO2NO2   ','PHOTON    ','MeOO      ','NO2       ','          ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna67    ',cs,0,0,107) ,         &
! 144
ratj_t1(144,'MeO2NO2   ','PHOTON    ','HCHO      ','HO2       ','NO3       ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna33    ',cs,0,0,107),          &
! 145
ratj_t1(145,'DHPCARB9  ','PHOTON    ','RN8OOH    ','CO        ','HO2       ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 9.64e02,'jforma    ',cs,0,0,119),          &
! 146
ratj_t1(146,'DHPCARB9  ','PHOTON    ','CARB6     ','HCHO      ','OH        ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 200 ,'jmhp      ',cs,0,0,119),             &
! 147
ratj_t1(147,'HPUCARB12 ','PHOTON    ','HUCARB9   ','CO        ','OH        ',  &
     'OH        ', 1.0,1.0,1.0,1.0,183.3*50 ,'jmacr2    ',cs,0,0,119),         &
! 148
ratj_t1(148,'HPUCARB12 ','PHOTON    ','CARB7     ','CO        ','OH        ',  &
     'HO2       ', 1.0,2.0,1.0,1.0,183.3*50 ,'jmacr2    ',cs,0,0,119),         &
! 149
ratj_t1(149,'HUCARB9   ','PHOTON    ','CARB6     ','CO        ','OH        ',  &
     'HO2       ', 0.0,0.0,0.0,0.0,183.3*50 ,'jmacr2    ',cs,0,0,119),         &
! 150
ratj_t1(150,'DHPR12OOH ','PHOTON    ','DHPCARB9  ','CO        ','OH        ',  &
     'HO2       ', 0.0,0.0,0.0,0.0, 100 ,'jprpal    ',cs,0,0,119),             &
! 151
ratj_t1(151,'DHPR12OOH ','PHOTON    ','CARB3     ','RN8OOH    ','OH        ',  &
     'OH        ', 0.0,0.0,0.0,0.0, 300 ,'jmhp      ',cs,0,0,119),             &
! 152
ratj_t1(152,'DHCARB9   ','PHOTON    ','CARB7     ','CO        ','HO2       ',  &
     'HO2       ', 0.0,0.0,0.0,0.0, 100 ,'jetcho    ',cs,0,0,119),             &
! 153
ratj_t1(153,'RU12NO3   ','PHOTON    ','CARB6     ','HOCH2CHO  ','HO2       ',  &
     'NO2       ', 0.0,0.0,0.0,0.0,7.583e-01,'jno2    ',cs,0,0,119),           &
! 154
ratj_t1(154,'RU10NO3   ','PHOTON    ','MeCO3     ','HOCH2CHO  ','NO2       ',  &
     '          ', 0.0,0.0,0.0,0.0, 1.213e-01 ,'jno2    ',cs,0,0,119),         &
! 155
ratj_t1(155,'NUCARB12  ','PHOTON    ','CARB7     ','CO        ','HO2       ',  &
     'NO2       ', 1.0,2.0,1.0,1.0, 6.066e-01 ,'jno2    ',cs,0,0,119),         &
! 156
ratj_t1(156,'UCARB12   ','PHOTON    ','RU12O2    ','HO2       ','          ',  &
     '          ', 0.0,0.0,0.0,0.0, 50.000,'jmacr2    ',cs,0,0,119),           &
! 157
ratj_t1(157,'UCARB12   ','PHOTON    ','CARB7     ','CO        ','HO2       ',  &
     '          ', 1.0,2.0,2.0,0.0, 25.000,'jmacr2    ',cs,0,0,119),           &
! 158 CS2
ratj_t1(158,'CARB3     ','PHOTON    ','CO        ','CO        ','H2        ',  &
     '          ', 0.0,0.0,0.0,0.0, 100.00,'jcb3b     ',cs,0,0,119)            &
]


! ----------------------------------------------------------------------


! Termolecular reactions
! Rates taken from  file "UKCA_Reaction_Rates - Termolecular Reactions.csv"
!on 11:41 2012/05/01

! Item number, reactant x2, product x2, broadening coefficient parameter F,
! low-pressure limit coefficients (k, alpha and beta), infinite-pressure limit
! coefficients (k, alpha and beta), product yields x 2, chemistry scheme,
! qualifier, disqualifier, version
ratt_defs_master(1:n_ratt_master)=[                                            &
ratt_t1(1,'O(3P)     ','O2        ','O3        ','m         ',     0.0,        &
  5.70e-34, -2.60, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, ti+t,0,0,107),        &
! T001 JPL 2011
ratt_t1(1,'O(3P)     ','O2        ','O3        ','m         ',     0.0,        &
  6.00e-34, -2.50, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, st+r,0,0,107),        &
! update for 121
ratt_t1(1,'O(3P)     ','O2        ','O3        ','m         ',     0.0,        &
  6.10e-34, -2.40, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, st,0,0,121),          &
! C
ratt_t1(1,'O(3P)     ','O2        ','O3        ','m         ',     0.0,        &
  5.60e-34, -2.60, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, s,0,0,107) ,          &
! (Slightly) different rate in CRI
ratt_t1(1,'O(3P)     ','O2        ','O3        ','          ',     0.0,        &
  6.00e-34, -2.60, 0.00,  0.00e+00,  0.0, 0.00, 0.0,0.000,cs,0,0,107),         &
! CS2 --> corrected for weighted average
ratt_t1(1,'O(3P)     ','O2        ','O3        ','          ',     0.0,        &
  5.68e-34, -2.60, 0.00,  0.00e+00,  0.0, 0.00, 0.0,0.000,cs,0,0,119),         &
! T002 JPL 2011
! Missing in T and TI schemes. Note the discrepancy in low-pres rate
! coefficient
ratt_t1(2,'O(3P)     ','NO        ','NO2       ','m         ',     0.6,        &
  9.00e-32, -1.50, 0.00,  3.00e-11,  0.0, 0.0, 0.0, 0.0, st+r+cs,0,0,107),     &
! do not do for cs as covered by cs2 below (119) - update for 121
ratt_t1(2,'O(3P)     ','NO        ','NO2       ','m         ',     0.6,        &
  9.10e-32, -1.50, 0.00,  3.00e-11,  0.0, 0.0, 0.0, 0.0, st,0,0,121),          &
ratt_t1(2,'O(3P)     ','NO        ','NO2       ','m         ',     0.6,        &
  9.00e-31, -1.50, 0.00,  3.00e-11,  0.0, 0.0, 0.0, 0.0, s,0,0,107),           &
! Rate change in CS2
ratt_t1(2,'O(3P)     ','NO        ','NO2       ','m         ',     0.85,       &
  1.00e-31, -1.60, 0.00,  5.00e-11,  -0.3, 0.0, 0.0, 0.0, cs,0,0,119),         &
! T003 JPL 2011, Missing in T and TI scheme
ratt_t1(3,'O(3P)     ','NO2       ','NO3       ','m         ',     0.6,        &
  2.50e-31, -1.80, 0.00,  2.20e-11, -0.7, 0.0, 0.0, 0.0, st,0,0,107),          &
! update for 121 - cs treated differently below
ratt_t1(3,'O(3P)     ','NO2       ','NO3       ','m         ',     0.6,        &
  3.40e-31, -1.60, 0.00,  2.30e-11, -0.2, 0.0, 0.0, 0.0, st,0,0,121),          &
ratt_t1(3,'O(3P)     ','NO2       ','NO3       ','m         ',     0.6,        &
  1.30e-31, -1.50, 0.00,  2.30e-11, 0.24, 0.0, 0.0, 0.0, s,0,0,107) ,          &
!  T003 - Different rate in CRI
ratt_t1(3,'O(3P)     ','NO2       ','NO3       ','          ',    0.60,        &
  9.00e-32, -2.00, 0.00, 2.20e-11,  0.0, 0.00, 0.0,0.000,cs,0,0,107),          &
! Rate change in CS2
ratt_t1(3,'O(3P)     ','NO2       ','NO3       ','          ',    0.60,        &
  1.30e-31, -1.50, 0.00, 2.30e-11,  0.24, 0.00, 0.0,0.000,cs,0,0,119),         &
! T004 JPL 2011
ratt_t1(4,'O(1D)     ','N2        ','N2O       ','m         ',     0.0,        &
  2.80e-36, -0.90, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, st+cs,0,0,107),       &
! update for 121
ratt_t1(4,'O(1D)     ','N2        ','N2O       ','m         ',    0.60,        &
  2.80e-36, -0.90, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, st+cs,0,0,121),       &
! JPL 2003
ratt_t1(4,'O(1D)     ','N2        ','N2O       ','m         ',     0.0,        &
  3.50e-37, -0.60, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, s,0,0,107) ,          &
! T005 JPL 2011
ratt_t1(5,'BrO       ','NO2       ','BrONO2    ','m         ',     0.6,        &
  5.20e-31, -3.20, 0.00,  6.90e-12,  0.0, 0.0, 0.0, 0.0, st+cs,0,0,107),       &
! update for 121
ratt_t1(5,'BrO       ','NO2       ','BrONO2    ','m         ',     0.6,        &
  5.50e-31, -3.10, 0.00,  6.60e-12, -2.9, 0.0, 0.0, 0.0, st+cs,0,0,121),       &
ratt_t1(5,'BrO       ','NO2       ','BrONO2    ','m         ',   327.0,        &
  4.70e-31, -3.10, 0.00,  1.40e-11, -1.2, 0.0, 0.0, 0.0, s,0,0,107) ,          &
! T006 JPL 2011
ratt_t1(6,'ClO       ','ClO       ','Cl2O2     ','m         ',     0.6,        &
  1.60e-32, -4.50, 0.00,  3.00e-12, -2.0, 0.0, 0.0, 0.0, st+cs,0,0,107),       &
! update for 121
ratt_t1(6,'ClO       ','ClO       ','Cl2O2     ','m         ',     0.6,        &
  1.90e-32, -3.60, 0.00,  3.70e-12, -1.6, 0.0, 0.0, 0.0, st+cs,0,0,121),       &
ratt_t1(6,'ClO       ','ClO       ','Cl2O2     ','m         ',     0.6,        &
  1.70e-32, -4.00, 0.00,  5.40e-12,  0.0, 0.0, 0.0, 0.0, s,0,0,107) ,          &
! T007 IUPAC 2005
ratt_t1(7,'Cl2O2     ','m         ','ClO       ','ClO       ',     0.5,        &
  3.70e-07,  0.00, 7690.0, 1.80e+14,  0.0,7690.0, 0.0, 0.0, st+cs,0,0,107),    &
! update for 121
ratt_t1(7,'Cl2O2     ','m         ','ClO       ','ClO       ',    0.45,        &
  3.70e-07,  0.00, 7690.0, 1.80e+14,  0.0,7690.0, 0.0, 0.0, st+cs,0,0,121),    &
ratt_t1(7,'Cl2O2     ','m         ','ClO       ','ClO       ',     0.6,        &
  1.00e-06,  1.00, 8000.0, 4.80e+15,  0.00,8820.0, 0.0, 0.0, s,0,0,107) ,      &
! T008 JPL 2011
ratt_t1(8,'ClO       ','NO2       ','ClONO2    ','m         ',     0.6,        &
  1.80e-31, -3.40, 0.00,  1.50e-11,  0.0, 0.0, 0.0, 0.0, st+cs,0,0,107),       &
! update for 121
ratt_t1(8,'ClO       ','NO2       ','ClONO2    ','m         ',     0.6,        &
  1.80e-31, -2.00, 0.00,  1.00e-10, -1.0, 0.0, 0.0, 0.0, st+cs,0,0,121),       &
ratt_t1(8,'ClO       ','NO2       ','ClONO2    ','m         ',   430.0,        &
  1.60e-31, -3.40, 0.00,  1.50e-11,  0.0, 0.0, 0.0, 0.0, s,0,0,107) ,          &
! T009 JPL 2011 22nd
ratt_t1(9,'H         ','O2        ','HO2       ','m         ',     0.6,        &
  4.40e-32, -1.30, 0.00,  7.50e-11,  0.0, 0.0, 0.0, 0.0, st+cs,0,0,107),       &
! JPL2003
ratt_t1(9,'H         ','O2        ','HO2       ','m         ',     0.6,        &
  5.70e-32, -1.60, 0.00,  7.50e-11,  0.0, 0.0, 0.0, 0.0, s,0,0,107),           &
! update for 121, include for s
ratt_t1(9,'H         ','O2        ','HO2       ','m         ',     0.6,        &
  5.30e-32, -1.80, 0.00,  9.50e-11,  0.4, 0.0, 0.0, 0.0, s+st+cs,0,0,121),     &
! T010 JPL 2011 see also asad_trimol.F90
ratt_t1(10,'HO2       ','HO2       ','H2O2      ','O2        ',     0.0,       &
  2.10e-33,  0.00,-920.0,  0.00e+00,  0.0, 0.0, 0.0, 0.0, st+r,0,0,107),       &
ratt_t1(10,'HO2       ','HO2       ','H2O2      ','O2        ',     0.0,       &
  1.90e-33,  0.00,-980.0,  0.00e+00,  0.0, 0.0, 0.0, 0.0, ti+s+t+cs,0,0,107),  &
! T011 JPL 2011
ratt_t1(11,'HO2       ','NO2       ','HO2NO2    ','m         ',     0.6,       &
  2.00e-31, -3.40, 0.00,  2.90e-12,  0.0, 0.0, 0.0, 0.0, st+r,0,0,107),        &
ratt_t1(11,'HO2       ','NO2       ','HO2NO2    ','m         ',     0.6,       &
  1.80e-31, -3.20, 0.00,  4.70e-12,  0.0, 0.0, 0.0, 0.0, ti+s+t+cs,0,0,107),   &
! Rate change in CS2
ratt_t1(11,'HO2       ','NO2       ','HO2NO2    ','m         ',     0.4,       &
  1.40e-31, -3.1, 0.00,  4.00e-12,  0.0, 0.0, 0.0, 0.0, cs,0,0,119),           &
! updates for 121
ratt_t1(11,'HO2       ','NO2       ','HO2NO2    ','m         ',     0.6,       &
  1.90e-31, -3.40, 0.00,  4.00e-12, -0.3, 0.0, 0.0, 0.0, st,0,0,121),          &
! T012 IUPAC 2001
ratt_t1(12,'HO2NO2    ','m         ','HO2       ','NO2       ',     0.5,       &
  4.10e-05,  0.00,10650.0, 4.80e+15,  0.0,11170.0,0.0, 0.0, t+st+r,0,0,107),   &
ratt_t1(12,'HO2NO2    ','m         ','HO2       ','NO2       ',     0.6,       &
  4.10e-05,  0.00,10650.0, 4.80e+15,  0.0,11170.0,0.0, 0.0, ti,0,0,107),       &
! Note the different rate
ratt_t1(12,'HO2NO2    ','m         ','HO2       ','NO2       ',     0.6,       &
  4.10e-06,  0.00,10650.0, 4.80e+15,  0.0,11170.0,0.0, 0.0, s,0,0,107) ,       &
!  T012 - Different rates in CRI
ratt_t1(12,'HO2NO2    ','m         ','NO2       ','HO2       ',    0.40,       &
  4.10e-05,  0.00,10650.0, 6.00e+15,  0.0,11170.0, 0.0, 0.0,cs,0,0,107),       &
! update for 121
ratt_t1(12,'HO2NO2    ','m         ','NO2       ','HO2       ',    0.40,       &
  4.10e-05,  0.00,10650.0, 6.00e+15,  0.0,11170.0, 0.0, 0.0,st,0,0,121),       &
!  T013
ratt_t1(13,'OH        ','NO        ','HONO      ','m         ',  1420.0,       &
  7.40e-31, -2.40, 0.00,  3.30e-11, -0.3, 0.0, 0.0, 0.0, ti+t+cs,0,0,107),     &
! CS2 rate changed
ratt_t1(13,'OH        ','NO        ','HONO      ','m         ',   0.81 ,       &
  7.40e-31, -2.40, 0.00,  3.30e-11, -0.3, 0.0, 0.0, 0.0,      cs,0,0,119),     &
! T013 JPL 2011
ratt_t1(13,'OH        ','NO        ','HONO      ','m         ',     0.6,       &
  7.00e-31, -2.60, 0.00,  3.60e-11, -0.10, 0.0, 0.0, 0.0, st,0,0,107),         &
! T014 JPL 2011
ratt_t1(14,'OH        ','NO2       ','HONO2     ','m         ',     0.6,       &
  1.80e-30, -3.00, 0.00,  2.80e-11,  0.0, 0.0, 0.0, 0.0, st+r,0,0,107),        &
ratt_t1(14,'OH        ','NO2       ','HONO2     ','m         ',     0.4,       &
  3.30e-30, -3.00, 0.00,  4.10e-11,  0.0, 0.0, 0.0, 0.0, ti+s+t,0,0,107),      &
!  T014 - Different rates in CRI
ratt_t1(14,'OH        ','NO2       ','HONO2     ','m         ',    0.60,       &
  2.60e-30, -3.20, 0.00, 2.40e-11, -1.30, 0.00, 0.0,0.000,cs,0,0,107),         &
!  T014 CS2 rate changed
ratt_t1(14,'OH        ','NO2       ','HONO2     ','m         ',    0.41,       &
  3.20e-30, -4.50, 0.00, 3.00e-11, 0.00, 0.00, 0.0,0.000,cs,0,0,119),          &
! T015 JPL 2011 - Added to CS scheme for Stratosphere
ratt_t1(15,'OH        ','OH        ','H2O2      ','m         ',     0.6,       &
  6.90e-31, -1.00, 0.00,  2.60e-11,  0.0, 0.0, 0.0, 0.0, st+cs,0,0,107),       &
ratt_t1(15,'OH        ','OH        ','H2O2      ','m         ',     0.5,       &
  6.90e-31, -0.80, 0.00,  2.60e-11,  0.0, 0.0, 0.0, 0.0, ti+s+t,0,0,107),      &
! T016 MCMv3.2
ratt_t1(16,'MeCO3     ','NO2       ','PAN       ','m         ',     0.3,       &
  2.70e-28, -7.10, 0.00,  1.20e-11, -0.9, 0.0, 0.0, 0.0, ti+t+st+r,0,0,107),   &
!  T016 - CRI (note still branched from MCM3.1...)
ratt_t1(16,'MeCO3     ','NO2       ','PAN       ','m         ',    0.60,       &
  8.50e-29, -6.50, 0.00, 1.10e-11, -1.00, 0.00, 0.0,0.000,cs,0,0,107),         &
! T016 CS2 rate changed
ratt_t1(16,'MeCO3     ','NO2       ','PAN       ','m         ',    0.30,       &
  3.28e-28, -6.87, 0.00, 1.125e-11, -1.105, 0.00, 0.0,0.000,cs,0,0,119),       &
! T017 MCMv3.2
ratt_t1(17,'PAN       ','m         ','MeCO3     ','NO2       ',     0.3,       &
  4.90e-03,  0.00,12100.0, 5.40e+16,  0.0,13830.0,0.0, 0.0, ti+t+st+r,0,0,107),&
!  T017 - Different rates in CRI
ratt_t1(17,'PAN       ','m         ','MeCO3     ','NO2       ',    0.30,       &
  1.10e-05,  0.0, 10100, 1.90e+17,  0.0, 14100, 0.0, 0.0,cs,0,0,107),          &
! T018 MCMv3.2
ratt_t1(18,'EtCO3     ','NO2       ','PPAN      ','m         ',     0.3,       &
  2.70e-28, -7.10, 0.00,  1.20e-11, -0.9, 0.0, 0.0, 0.0, ti+t+st,0,0,107),     &
!  T018 - Different rates in CRI
ratt_t1(18,'EtCO3     ','NO2       ','PPAN      ','m         ',    0.60,       &
  8.50e-29, -6.50, 0.00, 1.10e-11, -1.00, 0.00, 0.0,0.000,cs,0,0,107),         &
!  T018 CS2 rate changed
ratt_t1(18,'EtCO3     ','NO2       ','PPAN      ','m         ',    0.30,       &
  3.28e-28, -6.87, 0.00, 1.125e-11, -1.105, 0.00, 0.0,0.000,cs,0,0,119),       &
!  T019
ratt_t1(19,'PPAN      ','m         ','EtCO3     ','NO2       ',     0.4,       &
  1.70e-03,  0.00,11280.0, 8.30e+16,  0.0,13940.0,0.0, 0.0, ti+t,0,0,107),     &
! T019 MCMv3.2
ratt_t1(19,'PPAN      ','m         ','EtCO3     ','NO2       ',     0.3,       &
  4.90e-03,  0.00,12100.0, 5.40e+16,  0.0,13830.0,0.0, 0.0, st,0,0,107),       &
! T019 - Different rates in CRI
ratt_t1(19,'PPAN      ','m         ','EtCO3     ','NO2       ',    0.30,       &
  1.10e-05,  0.0, 10100, 1.90e+17,  0.0, 14100, 0.0, 0.0,cs,0,0,107),          &
! T020 MVMv3.2
ratt_t1(20,'MACRO2    ','NO2       ','MPAN      ','m         ',     0.3,       &
  2.70e-28, -7.10, 0.00,  1.20e-11, -0.9, 0.0, 0.0, 0.0, ti+st,0,0,107),       &
! updates for 121
ratt_t1(20,'MACRO2    ','NO2       ','MPAN      ','m         ',     0.3,       &
  3.28e-28, -6.87, 0.00,  1.13e-11, -1.105, 0.0, 0.0, 0.0, st,0,0,121),        &
! T021 MCMv3.2
ratt_t1(21,'MPAN      ','m         ','MACRO2    ','NO2       ',     0.3,       &
  4.90e-03,  0.00,12100.0, 5.40e+16,  0.0,13830.0,0.0, 0.0, ti+st,0,0,107),    &
! T022 IUPAC 2002
! Missing in S schemes
ratt_t1(22,'NO2       ','NO3       ','N2O5      ','m         ',     0.3,       &
  3.60e-30, -4.10, 0.00,  1.90e-12,  0.2, 0.0, 0.0, 0.0, ti+t+s+st+r,0,0,107), &
! updates for 121 - also do for s
ratt_t1(22,'NO2       ','NO3       ','N2O5      ','m         ',     0.6,       &
  2.40e-30, -3.00, 0.00,  1.60e-12,  0.1, 0.0, 0.0, 0.0, s+st,0,0,121),        &
!  T022 - Different rates in CRI
ratt_t1(22,'NO2       ','NO3       ','N2O5      ','m         ',    0.60,       &
  2.20e-30, -3.90, 0.00, 1.50e-12, -0.70, 0.00, 0.0,0.000,cs,0,0,107),         &
!  T022 !! CS2 rate changed
ratt_t1(22,'NO2       ','NO3       ','N2O5      ','m         ',    0.35,       &
  3.60e-30, -4.10, 0.00, 1.90e-12, 0.20, 0.00, 0.0,0.000,cs,0,0,119),          &
! T023 IUPAC 2002 - There is a rounding error in ffac (should be 0.35)!
ratt_t1(23,'N2O5      ','m         ','NO2       ','NO3       ',     0.3,       &
  1.30e-03, -3.50,11000.0, 9.70e+14,0.1,11080.0,0.0, 0.0,ti+s+t+st+r,0,0,107), &
! updates for 121 - now 0.35 for st
ratt_t1(23,'N2O5      ','m         ','NO2       ','NO3       ',    0.35,       &
  1.30e-03, -3.50,11000.0, 9.70e+14, 0.10, 11080.0,0.0, 0.0,st,0,0,121),       &
!  T023 - CRI
ratt_t1(23,'N2O5      ','m         ','NO3       ','NO2       ',    0.35,       &
  1.30e-03, -3.50, 11000, 9.70e+14,  0.10, 11080, 0.0, 0.0,cs,0,0,107),        &
! T024 IUPAC 2001
! not in TI/TI scheme
ratt_t1(24,'NO        ','NO        ','NO2       ','NO2       ',     0.0,       &
  3.30e-39,  0.00,-530.0,  0.00e+00,  0.0, 0.0, 0.0, 0.0, t+st,0,0,107),       &
! reaction with [O2] NOT [M] -> Should be factor 0.21 slower than ST rate
ratt_t1(24,'NO        ','NO        ','NO2       ','NO2       ',     0.0,       &
  6.93e-40,  0.00,-530.0,  0.00e+00,  0.0, 0.0, 0.0, 0.0, s+cs,0,0,107),       &
! T025
ratt_t1(25,'SO2       ','OH        ','SO3       ','HO2       ',     0.6,       &
  3.00e-31, -3.30, 0.00,  1.50e-12,  0.0, 0.0, 0.0, 0.0, st+s,a,0,107),        &
! updates for 121 - do for st & s
ratt_t1(25,'SO2       ','OH        ','SO3       ','HO2       ',     0.6,       &
  2.90e-31, -4.10, 0.00,  1.70e-12,  0.2, 0.0, 0.0, 0.0, st+s,a,0,121),        &
ratt_t1(25,'SO2       ','OH        ','HO2       ','H2SO4     ',     0.6,       &
  3.00e-31, -3.30, 0.00,  1.50e-12,  0.0, 0.0, 0.0, 0.0, ti,a,0,107),          &
ratt_t1(25,'SO2       ','OH        ','H2SO4     ','          ',     0.6,       &
  3.00e-31, -3.30, 0.00,  1.50e-12,  0.0, 0.0, 0.0, 0.0, ol,a,0,107),          &
!  T025 - Different rates in CRI
ratt_t1(25,'SO2       ','OH        ','SO3       ','HO2       ',    0.53,       &
  2.50e-31, -2.60, 0.00, 2.00e-12,  0.0, 0.00, 0.0, 0.0,cs,a,0,107),           &
! T026
ratt_t1(26,'O(3P)S    ','O2        ','O3S       ','m         ',     0.0,       &
  5.70e-34, -2.60, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, t,0,0,107),           &
! T027
ratt_t1(27,'OH        ','C2H4      ','HOC2H4O2  ','          ',     0.0,       &
  8.60e-29, -3.10, 0.00,  9.00e-12, -0.85, 0.0, 0.0, 0.0, r,0,0,107),          &
!  T027 - different ffac in CRI
ratt_t1(27,'OH        ','C2H4      ','HOCH2CH2O2','m         ',    0.48,       &
  8.60e-29, -3.10, 0.00, 9.00e-12, -0.85, 0.00, 0.0,0.000,cs,0,0,107),         &
! T028
ratt_t1(28,'OH        ','C3H6      ','HOC3H6O2  ','          ',     0.0,       &
  8.00e-27, -3.50, 0.00,  3.00e-11, -1.00, 0.0, 0.0, 0.0, r,0,0,107),          &
!  T028 - different ffac in CRI
ratt_t1(28,'OH        ','C3H6      ','RN9O2     ','m         ',    0.50,       &
  8.00e-27, -3.50, 0.00, 3.00e-11, -1.00, 0.00, 0.0,0.000,cs,0,0,107),         &
! From here on, CRI reactions
!  T029
ratt_t1(29,'MeO2NO2   ','m         ','MeOO      ','NO2       ',    0.36,       &
  9.00e-05,  0.0,  9690, 1.10e+16,  0.0, 10560, 0.0, 0.0,cs,0,0,107),          &
!  T030 - Not sure why this is in CRI and not the others...
ratt_t1(30,'O(3P)     ','SO2       ','SO3       ','m         ',     0.0,       &
  4.00e-32,  0.0,  1000,  0.00e+00,  0.0, 0.00, 0.0,0.000,cs,a,0,107),         &
!  T031
ratt_t1(31,'C2H2      ','OH        ','HCOOH     ','CO        ',    0.37,       &
  1.82e-30, -1.50, 0.00, 3.64e-13,  0.0, 0.00, 0.0, 0.0,cs,0,0,107),           &
!  T032
ratt_t1(32,'C2H2      ','OH        ','CARB3     ','OH        ',    0.37,       &
  3.18e-30, -1.50, 0.00, 6.36e-13,  0.0, 0.00, 0.0, 0.0,cs,0,0,107),           &
!  T033
ratt_t1(33,'HOCH2CO3  ','NO2       ','PHAN      ','m         ',    0.60,       &
  8.50e-29, -6.50, 0.00, 1.10e-11, -1.00, 0.00, 0.0,0.000,cs,0,0,107),         &
!  T033 CS2 rate changed
ratt_t1(33,'HOCH2CO3  ','NO2       ','PHAN      ','m         ',    0.30,       &
  3.28e-28, -6.87, 0.00, 1.125e-11, -1.105, 0.00, 0.0,0.000,cs,0,0,119),       &
!  T034
ratt_t1(34,'PHAN      ','m         ','HOCH2CO3  ','NO2       ',    0.30,       &
  1.10e-05,  0.0, 10100, 1.90e+17,  0.0, 14100, 0.0, 0.0,cs,0,0,107),          &
!  T035
ratt_t1(35,'RU12O2    ','NO2       ','RU12PAN   ','m         ',    0.60,       &
  5.19e-30, -6.50, 0.00, 6.71e-13, -1.00, 0.00, 0.0,0.000,cs,0,0,107),         &
!  T035
ratt_t1(35,'RU12O2    ','m         ','DHCARB9   ','CO        ',    0.0,        &
  0.0, 0.0, 0.00, 1.20e5,0.0 , 5300.0 , 2.0,2.0,cs,0,0,119),                   &
!  T036
ratt_t1(36,'RU12PAN   ','m         ','RU12O2    ','NO2       ',    0.30,       &
  1.10e-05,  0.0, 10100, 1.90e+17,  0.0, 14100, 0.0, 0.0,cs,0,0,107),          &
! T036
ratt_t1(36,'RU12O2    ','m         ','OH        ','m         ',     0.0,       &
  0.0,  0.0, 0.0, 1.20e5, 0.0, 5300, 2.00, 0.00,cs,0,0,119),                   &
!  T037
ratt_t1(37,'RU10O2    ','NO2       ','MPAN      ','m         ',    0.60,       &
  3.48e-30, -6.50, 0.00, 4.51e-13, -1.00, 0.00, 0.0,0.000,cs,0,0,107),         &
! T037
ratt_t1(37,'MACO3     ','NO2       ','MPAN      ','m         ',    0.30,       &
  3.28e-28, -6.87, 0.00, 1.125e-11, -1.105, 0.00, 0.0,0.000,cs,0,0,119),       &
!  T038
ratt_t1(38,'MPAN      ','m         ','RU10O2    ','NO2       ',    0.30,       &
  1.10e-05,  0.0, 10100, 1.90e+17,  0.0, 14100, 0.0, 0.0,cs,0,0,107),          &
! T038
ratt_t1(38,'MPAN      ','m         ','MACO3     ','NO2       ',    0.0,        &
  0.0,  0.0,0.0, 1.60e16,  0.0, 13500.0, 0.0, 0.0,cs,0,0,119),                 &
!  T039
ratt_t1(39,'RTN26O2   ','NO2       ','RTN26PAN  ','m         ',    0.60,       &
  6.14e-29, -6.50, 0.00, 7.94e-12, -1.00, 0.00, 0.0,0.000,cs,0,0,107),         &
! T039
ratt_t1(39,'RTN26O2   ','NO2       ','RTN26PAN  ','m         ',    0.30,       &
  2.368e-28, -6.87, 0.00, 8.123e-12, -1.105, 0.00, 0.0,0.000,cs,0,0,119),      &
!  T040
ratt_t1(40,'RTN26PAN  ','m         ','RTN26O2   ','NO2       ',    0.30,       &
  1.10e-05,  0.0, 10100, 1.90e+17,  0.0, 14100, 0.0, 0.0,cs,0,0,107),          &
!  T041
ratt_t1(41,'MeOO      ','NO2       ','MeO2NO2   ','m         ',    0.36,       &
  2.50e-30, -5.50, 0.00, 1.80e-11,  0.0, 0.00, 0.0,0.000,cs,0,0,107),          &
!  T042 - DMS + OH + O2, IUPAC 2002, See asad_trimol
ratt_t1(42,'DMS       ','OH        ','DMSO      ','HO2       ',     0.0,       &
  9.50e-39,  0.00, -5270, 7.50e-29,  0.00, -5610, 0.00, 0.00,cs,a,0,107),      &
! updates for 121 - include for st - see also T051, replaces B201
ratt_t1(42,'DMS       ','OH        ','DMSO      ','HO2       ',     0.0,       &
  3.80e-39,  0.00, -5270, 3.00e-29,  0.00, -5610, 0.00, 0.00,st,a,0,121),      &
!  T043
ratt_t1(43,'RU14O2    ','m         ','UCARB10   ','HCHO      ',     0.0,       &
  0.0,  0.0, 0.0, 6.20e10, 0.0, 9750, 2.00, 2.00,cs,0,0,119),                  &
!  T044
ratt_t1(44,'RU14O2    ','m         ','OH        ','m         ',     0.0,       &
  0.0,  0.0, 0.0, 6.20e10, 0.0, 9750, 2.00, 0.00,cs,0,0,119),                  &
!  T045
ratt_t1(45,'RU14O2    ','m         ','HPUCARB12 ','HO2       ',     0.0,       &
  0.0,  0.0,0.0,1.38e7 , 0.0,6579 ,0.00, 0.00,cs,0,0,119),                     &
!  T046
ratt_t1(46,'RU14O2    ','m         ','DHPR12O2  ','m         ',     0.0,       &
  0.0,  0.0,0.0,1.38e7 , 0.0, 6759,0.00, 0.00,cs,0,0,119),                     &
!  T047
ratt_t1(47,'DHPR12O2  ','m         ','DHPCARB9  ','CO        ',     0.0,       &
  0.0,  0.0,0.0,1.50e7 ,0.0, 5300.00,2.00, 2.00,cs,0,0,119),                   &
!  T048
ratt_t1(48,'DHPR12O2  ','m         ','OH        ','m         ',     0.0,       &
  0.0,  0.0,0.0,1.50e7 ,0.0, 5300.00,2.00, 2.00,cs,0,0,119),                   &
!  T049
ratt_t1(49,'RU10AO2   ','m         ','CARB7     ','CO        ',     0.0,       &
  0.0,  0.0,0.0,1.50e7 ,0.0, 5300.00,2.00, 2.00,cs,0,0,119),                   &
!  T050
ratt_t1(50,'RU10AO2   ','m         ','OH        ','m         ',     0.0,       &
  0.0,  0.0,0.0,1.50e7 ,0.0, 5300.00,2.00, 2.00,cs,0,0,119),                   &
!  T051
! updates for 121 - include for st - see also T042, replaces B201
ratt_t1(51,'DMS       ','OH        ','SO2       ','MeOO      ',     0.0,       &
  5.70e-39,  0.00, -5270, 4.50e-29,  0.00, -5610, 1.00, 2.00,st,a,0,121)       &
]

!----------------------------------------------------------------------
! NOTES: CheST Termolecular Reactions
!----------------------------------------------------------------------
! T001 O(3P)+O2 -> O3 m JPL 2011
! T001 IUPAC 2002 recommend k = 5.623E-34*(T/300)^-2.6 (based on
! T001 weighted mean of kO2+kN2)
!----------------------------------------------------------------------
! T002 O(3P)+NO -> NO2 m JPL 2011
! T002 IUPAC 2002 recommend k0 = 1.0E-3*(T/300)-1.6*[N2]
!----------------------------------------------------------------------
! T003 O(3P)+NO2 -> NO3 m JPL 2011
! T003 IUPAC 2002 recommend k0 = 1.3E-31*(T/300)^-1.5 kinf =
! T003 2.3E-11*(T/300)^0.24 Fc = 0.6
!----------------------------------------------------------------------
! T004 O(1D)+N2 -> N2O m JPL 2011
! T004 IUPAC 2002 k0=2.8E-36*[N2]
!----------------------------------------------------------------------
! T005 BrO+NO2 -> BrONO2 m JPL 2011
! T005 IUPAC Fc = 0.55 k0 = 4.2E-31*exp(T/300)^-2.4*[N2] kinf=2.7E-11
!----------------------------------------------------------------------
! T006 ClO+ClO -> Cl2O2 m JPL 2011
! T006 IUPAC recommend k0 = 2.0E-32*(T/300)^-4*[N2] kinf = 1.0E-11.
! T006 Fc=0.45 In general this is a problematic rate
!----------------------------------------------------------------------
! T007 Cl2O2+m -> ClO ClO IUPAC 2005
! T007 No JPL data
!----------------------------------------------------------------------
! T008 ClO+NO2 -> ClONO2 m JPL 2011
! T008 IUPAC recommend k0 = 1.6E-31*(T/300)^-3.4*[N2] kinf = 7.0E-11
! T008 Fc=0.4
!----------------------------------------------------------------------
! T009 H+O2 -> HO2 m JPL 2011
! T009 IUPAC 2009 recommend k0 = 4.3E-32*(T/300)^-1.2*[N2] kinf =
! T009 9.6E-11 and FC(ent) = 0.5
!----------------------------------------------------------------------
! T010 HO2+HO2 -> H2O2 O2 JPL 2011 see also asad_trimol.F90
! T010 IUPAC (2001) k = 1.9E-33*exp(980/T)*[N2]. Note that this reaction
! T010 is special and in the presence of H2O the rate constant needs to
! T010 be adjusted by: {1 + 1.4E-21*[H2O]*exp(2200/T)} (same H2O
! T010 expression for JPL). JPL also include the formation of HO2-H2O as
! T010 a separate species.
!----------------------------------------------------------------------
! T011 HO2+NO2 -> HO2NO2 m JPL 2011
! T011 IUPAC 2009 k0 = 1.4E-31*(T/300)^-3.1*[N2] kinf = 4.0E-12 Fc = 0.4
!----------------------------------------------------------------------
! T012 HO2NO2+m -> HO2 NO2 IUPAC 2001
! T012 No JPL data.. NOTE should multiply the expression by [N2] NOT [M]
!----------------------------------------------------------------------
! T013 OH+NO -> HONO m JPL 2011
! T013 IUPAC 2002 recommend k0 = 7.40E-31*(T/300)^-2.4 kinf =
! T013 3.3E-11*(T/300)^-0.3 Fc = 0.81
!----------------------------------------------------------------------
! T014 OH+NO2 -> HONO2 m JPL 2011
! T014 IUPCA 2009 recommend k0 = 3.3E-30*(T/300)^-3.0 kinf = 6.0E-11 Fc
! T014 = 0.4
!----------------------------------------------------------------------
! T015 OH+OH -> H2O2 m JPL 2011
! T015 IUPAC 2009 recommend k0 = 6.9E-31*(T/300)^-0.8 kinf =
! T015 3.9E-11*(T/300)^-0.47
!----------------------------------------------------------------------
! T016 MeCO3+NO2 -> PAN m MCMv3.2
! T016 Based on IUPAC 2003. JPL recommend k0 = 9.7E-29*(T/300)^-5.6 kinf
! T016 = 9.3E-12*(T/300)^-1.5 Fc = 0.6
!----------------------------------------------------------------------
! T017 PAN+m -> MeCO3 NO2 MCMv3.2
! T017 IUPAC 2003
!----------------------------------------------------------------------
! T018 EtCO3+NO2 -> PPAN m MCMv3.2
! T018 JPL recommend k0 = 9.0E-28*(T/300)^-8.9 kinf =
! T018 7.7E-12*(T/300)^-0.2 Fc = 0.6. IUPAC 2003 k=
! T018 1.70E-3*exp(-11280/T) kinf = 8.3E16*exp(-13940/T) Fc = 0.36
!----------------------------------------------------------------------
! T019 PPAN+m -> EtCO3 NO2 MCMv3.2
! T019 IUPAC 2003 k0 = 1.70E-03*exp(-11280/T) kinf =
! T019 8.30E+16*exp(-13940/T) Fc = 0.36
!----------------------------------------------------------------------
! T020 MACRO2+NO2 -> MPAN m MVMv3.2
! T020 -
!----------------------------------------------------------------------
! T021 MPAN+m -> MACRO2 NO2 MCMv3.2
! T021 IUPAC recommend k = 1.6E16*exp(-13500/T). No JPL data
!----------------------------------------------------------------------
! T022 NO2+NO3 -> N2O5 m IUPAC 2002
! T022 JPL recommend k0 = 2.0E-30*(T/300)^-4.4 kinf =
! T022 1.4E-12*(T/300)^-0.7 Fc = 0.6
!----------------------------------------------------------------------
! T023 N2O5+m -> NO2 NO3 IUPAC 2002
! T023 No JPL data. NOTE there is code in asad_trimol.F90 which looks
! T023 for this reaction and modifies the rate by an extra factor. I can
! T023 NOT find why it does this in the literature.
!----------------------------------------------------------------------
! T024 NO+NO -> NO2 NO2 IUPAC 2001
! T024 No JPL data..
! ---------------------------------------------------------------------


! Wet deposition.

! The following formula is used to calculate the effective Henry's Law
! coefficient, which takes the effects of dissociation and complex formation
! on a species' solubility into account: (see Giannakopoulos, 1998)
!
!       H(eff) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
!
! The data in columns 1 & 2 above give the data for this gas-aqueous transfer,
!       Column 1 = K(298) [M/atm]
!       Column 2 = -deltaH/R [K-1]
!
! If the species dissociates in the aqueous phase, the above term is multiplied
! by 1+{K(aq)/[H+]}, where
!       K(aq) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
! The data in columns 3 & 4 give the data for this aqueous-phase dissociation,
!       Column 3 = K(298) [M]
!       Column 4 = -deltaH/R [K-1]
! The data in columns 5 and 6 give the data for a second dissociation,
! e.g for SO2, HSO3^{-}, and SO3^{2-}
!       Column 5 = K(298) [M]
!       Column 6 = -deltaH/R [K-1]

! Item number, species name, Henry's coefficients x6, chemistry scheme,
! qualifier, disqualifier, version

henry_defs_master(1:100)=[                                                     &
! WD: 1
wetdep(1,'O3        ',                                                         &
[0.113e-01,0.23e+04,0.0e+00,0.0e+00,0.0e+00,0.0e+00],ti+s+st+ol+cs,0,0,107),   &
! WD: 2
wetdep(2,'NO3       ',                                                         &
[0.20e+01,0.20e+04,0.00e+00,0.00e+00,0.00e+00,0.0e+00],ti+t+st+r+cs,0,0,107),  &
wetdep(2,'NO3       ',                                                         &
[0.60e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.0e+00],s,0,0,107),             &
! WD: 3
wetdep(3,'N2O5      ',                                                         &
[0.21e+06,0.87e+04,0.20e+02,0.00e+00,0.00e+00,0.0e+00],ti+t+st+r+cs,0,0,107),  &
wetdep(3,'N2O5      ',                                                         &
[0.21e+06,0.87e+04,0.157e+02,0.00e+00,0.00e+00,0.0e+00],s,0,0,107),            &
! WD: 4
wetdep(4,'HO2NO2    ',                                                         &
[0.13e+05,0.69e+04,0.10e-04,0.00e+00,0.00e+00,0.0e+00],ti+t+st+r+cs,0,0,107),  &
wetdep(4,'HO2NO2    ',                                                         &
[0.20e+05,0.00e+00,0.10e-04,0.00e+00,0.00e+00,0.00e+00],s,0,0,107),            &
! WD: 5
wetdep(5,'HONO2     ',                                                         &
[0.21e+06,0.87e+04,0.20e+02,0.00e+00,0.00e+00,0.0e+00],ti+t+st+r+cs,0,0,107),  &
wetdep(5,'HONO2     ',                                                         &
[0.21e+06,0.87e+04,0.157e+02,0.00e+00,0.00e+00,0.00e+00],s,0,0,107),           &
! WD: 6
wetdep(6,'H2O2      ',                                                         &
[0.83e+05,0.74e+04,0.24e-11,-0.373e+04,0.00e+00,0.00e+00],ti+t+st+ol+r+cs,     &
    0,0,107),                                                                  &
wetdep(6,'H2O2      ',                                                         &
[0.83e+05,0.74e+04,0.22e-11,-0.373e+04,0.00e+00,0.00e+00],s,0,0,107),          &
! WD: 7
wetdep(7,'HCHO      ',                                                         &
[0.33e+04,0.65e+04,0.00e+00,0.00e+00,0.00e+00,0.0e+00],ti+t+st+r+cs,0,0,107),  &
wetdep(7,'HCHO      ',                                                         &
[0.30e+04,0.72e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],s,0,0,107),            &
! WD: 8
wetdep(8,'MeOO      ',                                                         &
[0.20e+04,0.66e+04,0.00e+00,0.00e+00,0.00e+00,0.0e+00],ti+t+st+r+cs,0,0,107),  &
wetdep(8,'MeOO      ',                                                         &
[0.20e+04,0.564e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],s,0,0,107),           &
! WD: 9
wetdep(9,'MeOOH     ',                                                         &
[0.31e+03,0.50e+04,0.00e+00,0.00e+00,0.00e+00,0.0e+00],ti+t+st+r+cs,0,0,107),  &
wetdep(9,'MeOOH     ',                                                         &
[0.31e+03,0.52e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],s,0,0,107),            &
! WD: 10
wetdep(10,'HO2       ',                                                        &
[0.40e+04,0.59e+04,0.20e-04,0.00e+00,0.00e+00,0.0e+00],ti+t+st+r+cs,0,0,107),  &
wetdep(10,'HO2       ',                                                        &
[0.40e+04,0.59e+04,0.16e-04,0.00e+00,0.00e+00,0.00e+00],s,0,0,107),            &
! WD: 11
wetdep(11,'BrONO2    ',                                                        &
[0.21e+06,0.87e+04,0.157e+03,0.00e+00,0.00e+00,0.00e+00],st+cs,0,0,107),       &
wetdep(11,'BrONO2    ',                                                        &
[0.00e+05,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00],s,0,0,107),            &
! WD: 12
wetdep(12,'HCl       ',                                                        &
[0.19e+02,0.60e+03,0.10e+05,0.00e+00,0.00e+00,0.00e+00],s+st+cs,0,0,107),      &
! WD: 13
wetdep(13,'HOCl      ',                                                        &
[0.93e+03,0.59e+04,0.32e-07,0.00e+00,0.00e+00,0.00e+00],s+st+cs,0,0,107),      &
! WD: 14
wetdep(14,'HBr       ',                                                        &
[0.13e+01,0.102e+05,0.10e+10,0.00e+00,0.00e+00,0.00e+00],s+st+cs,0,0,107),     &
! WD: 15
wetdep(15,'HOBr      ',                                                        &
[0.61e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00],s+st+cs,0,0,107),      &
! WD: 16
wetdep(16,'ClONO2    ',                                                        &
[0.21e+06,0.87e+04,0.157e+02,0.00e+00,0.00e+00,0.00e+00],s+st+cs,0,0,107),     &
! WD: 17
wetdep(17,'HONO      ',                                                        &
[0.50e+02,0.49e+04,0.56e-03,-0.126e+04,0.00e+00,0.0e+00],ti+t+st+cs,0,0,107),  &
! WD: 18
wetdep(18,'EtOOH     ',                                                        &
[0.34e+03,0.57e+04,0.00e+00,0.00e+00,0.00e+00,0.0e+00],ti+t+st+r+cs,0,0,107),  &
! WD: 19
wetdep(19,'n-PrOOH   ',                                                        &
[0.34e+03,0.57e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],ti+t+st,0,0,107),      &
! WD: 20
wetdep(20,'i-PrOOH   ',                                                        &
[0.34e+03,0.57e+04,0.00e+00,0.00e+00,0.00e+00,0.0e+00],ti+t+st+r+cs,0,0,107),  &
! WD: 21
wetdep(21,'MeCOCH2OOH',                                                        &
[0.34e+03,0.57e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],ti+t+st,0,0,107),      &
! WD: 22
wetdep(22,'ISOOH     ',                                                        &
[0.17e+07,0.97e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],ti+st+r,0,0,107),      &
! WD: 23
wetdep(23,'ISON      ',                                                        &
[0.30e+04,0.74e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],ti+st+r,0,0,107),      &
! WD: 24
wetdep(24,'MACROOH   ',                                                        &
[0.17e+07,0.97e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],ti+st,0,0,107),        &
! WD: 25
wetdep(25,'HACET     ',                                                        &
[0.14e+03,0.72e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],ti+st,0,0,107),        &
! WD: 26
wetdep(26,'MGLY      ',                                                        &
[0.35e+04,0.72e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],ti+st+r,0,0,107),      &
! WD: 27
wetdep(27,'HCOOH     ',                                                        &
[0.69e+04,0.56e+04,0.18e-03,-0.151e+04,0.00e+00,0.00e+00],ti+st+cs,0,0,107),   &
! WD: 28
wetdep(28,'MeCO3H    ',                                                        &
[0.75e+03,0.53e+04,0.63e-08,0.00e+00,0.00e+00,0.00e+00],ti+st+cs,0,0,107),     &
! WD: 29
wetdep(29,'MeCO2H    ',                                                        &
[0.47e+04,0.60e+04,0.18e-04,0.00e+00,0.00e+00,0.00e+00],ti+st+cs,0,0,107),     &
! WD: 30
wetdep(30,'MeOH      ',                                                        &
[0.23e+03,0.49e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],ti+st+r+cs,0,0,107),   &
! WD: 31
wetdep(31,'SO2       ',                                                        &
[0.123e+01,0.302e+04,0.123e-01,0.201e+04,0.60e-07,0.112e+04],ti+s+st+ol+cs,    &
    a,0,107),                                                                  &
! WD: 32
wetdep(32,'DMSO      ',                                                        &
[0.50e+05,0.6425e+04,0.00e+00,0.00e+00,0.0e+00,0.0e+00],ti+st+ol+cs,a,0,107),  &
! WD: 33
wetdep(33,'NH3       ',                                                        &
[0.10e+07,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00],ti+st+cs,a,0,107),     &
! WD: 34
wetdep(34,'Sec_Org   ',                                                        &
[0.10e+06,0.12e+02,0.00e-00,0.00e+00,0.0e+00,0.0e+00],ti+st+ol+cs,a,0,107),    &
! WD: 35
wetdep(35,'HO2S      ',                                                        &
[0.40e+04,0.59e+04,0.20e-04,0.00e+00,0.00e+00,0.00e+00],t,0,0,107),            &
! WD: 36
wetdep(36,'MVKOOH    ',                                                        &
[0.17e+07,0.97e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],r,0,0,107),            &
! WD: 37
wetdep(37,'ORGNIT    ',                                                        &
[0.13e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00],r,0,0,107),            &
! WD: 38
wetdep(38,'s-BuOOH   ',                                                        &
[0.34e+03,0.57e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],r,0,0,107),            &
! WD: 39
wetdep(39,'GLY       ',                                                        &
[0.36e+06,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00],r,0,0,107),            &
! WD 40
wetdep( 40, 'EtOH      ',                                                      &
[2.30e+02,4.90e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 41. Using MeOH Henry's Law Coeffs
wetdep( 41, 'i-PrOH    ',                                                      &
[2.30e+02,4.90e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 42. Using MeOH Henry's Law Coeffs
wetdep( 42, 'n-PrOH    ',                                                      &
[2.30e+02,4.90e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 43. Using HOCH2CHO Henry's Law Coeffs
wetdep( 43, 'HOCH2CHO  ',                                                      &
[4.15e+04,4.60e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 44. Using EtOOH Henry's Law Coeffs
wetdep( 44, 'HOC2H4OOH ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020))
wetdep( 44, 'HOC2H4OOH ',                                                      &
[1.90e+06,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 45. Using MeCO3H Henry's Law Coeffs
wetdep( 45, 'EtCO3H    ',                                                      &
[7.50e+02,5.30e+03,6.30e-09,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 46. Using MeCO3H Henry's Law Coeffs
wetdep( 46, 'HOCH2CO3H ',                                                      &
[7.50e+02,5.30e+03,6.30e-09,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 47. Using NALD Henry's Law Coeffs
wetdep( 47, 'NOA       ',                                                      &
[1.00e+03,7.40e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 48. Using HACET Henry's Law Coeffs
wetdep( 48, 'CARB7     ',                                                      &
[1.40e+02,7.20e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep
wetdep( 48, 'CARB7     ',                                                      &
[1.46e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 49. Using HACET Henry's Law Coeffs
wetdep( 49, 'CARB10    ',                                                      &
[1.40e+02,7.20e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep 62nd
wetdep( 49, 'CARB10    ',                                                      &
[1.46e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 50. Using HACET Henry's Law Coeffs
wetdep( 50, 'CARB13    ',                                                      &
[1.40e+02,7.20e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep
wetdep( 50, 'CARB13    ',                                                      &
[1.46e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 51. Using HACET Henry's Law Coeffs
wetdep( 51, 'CARB16    ',                                                      &
[1.40e+02,7.20e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep
wetdep( 51, 'CARB16    ',                                                      &
[1.46e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 52. Using MGLY Henry's Law Coeffs
wetdep( 52, 'CARB3     ',                                                      &
[3.50e+03,7.20e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep
wetdep( 52, 'CARB3     ',                                                      &
[4.19e+05,7.48e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 53. Using MGLY Henry's Law Coeffs
wetdep( 53, 'CARB6     ',                                                      &
[3.50e+03,7.20e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 54. Using MGLY Henry's Law Coeffs
wetdep( 54, 'CARB9     ',                                                      &
[3.50e+03,7.20e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 55. Using MGLY Henry's Law Coeffs
wetdep( 55, 'CARB12    ',                                                      &
[3.50e+03,7.20e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 56. Using MGLY Henry's Law Coeffs
wetdep( 56, 'CARB15    ',                                                      &
[3.50e+03,7.20e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 57. Using NALD Henry's Law Coeffs
wetdep( 57, 'NUCARB12  ',                                                      &
[1.00e+03,7.40e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 58. Using EtOOH Henry's Law Coeffs
wetdep( 58, 'RN10OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 59. Using EtOOH Henry's Law Coeffs
wetdep( 59, 'RN13OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 60. Using EtOOH Henry's Law Coeffs
wetdep( 60, 'RN16OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 61. Using EtOOH Henry's Law Coeffs
wetdep( 61, 'RN19OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 62. Using EtOOH Henry's Law Coeffs
wetdep( 62, 'RN8OOH    ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) HYPERACT)
wetdep( 62, 'RN8OOH    ',                                                      &
[1.16e+04,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 63. Using EtOOH Henry's Law Coeffs
wetdep( 63, 'RN11OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) HYPERACT)
wetdep( 63, 'RN11OOH   ',                                                      &
[1.16e+04,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 64. Using EtOOH Henry's Law Coeffs
wetdep( 64, 'RN14OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) HYPERACT)
wetdep( 64, 'RN14OOH   ',                                                      &
[1.16e+04,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 65. Using EtOOH Henry's Law Coeffs
wetdep( 65, 'RN17OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) HYPERACT)
wetdep( 65, 'RN17OOH   ',                                                      &
[1.16e+04,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 66. Using EtOOH Henry's Law Coeffs
wetdep( 66, 'RU14OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) ISOPOOH)
wetdep( 66, 'RU14OOH   ',                                                      &
[3.50e+06,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 67. Using EtOOH Henry's Law Coeffs
wetdep( 67, 'RU12OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) MVKOOH)
wetdep( 67, 'RU12OOH   ',                                                      &
[1.24e+06,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 68. Using EtOOH Henry's Law Coeffs
wetdep( 68, 'RU10OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) MVKOOH)
wetdep( 68, 'RU10OOH   ',                                                      &
[1.24e+06,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 69. Using EtOOH Henry's Law Coeffs
wetdep( 69, 'NRU14OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) ISOPNOOH)
wetdep( 69, 'NRU14OOH  ',                                                      &
[8.75e+04,6.00e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 70. Using EtOOH Henry's Law Coeffs
wetdep( 70, 'NRU12OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) MVKN)
wetdep( 70, 'NRU12OOH  ',                                                      &
[1.84e+05,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 71. Using EtOOH Henry's Law Coeffs
wetdep( 71, 'RN9OOH    ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep
wetdep( 71, 'RN9OOH    ',                                                      &
[1.50e+06,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 72. Using EtOOH Henry's Law Coeffs
wetdep( 72, 'RN12OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 73. Using EtOOH Henry's Law Coeffs
wetdep( 73, 'RN15OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 74. Using EtOOH Henry's Law Coeffs
wetdep( 74, 'RN18OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107)            &
]



henry_defs_master(101:n_wet_master)=[                                          &
! WD: 75. Using EtOOH Henry's Law Coeffs
wetdep( 75, 'NRN6OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 76. Using EtOOH Henry's Law Coeffs
wetdep( 76, 'NRN9OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 77. Using EtOOH Henry's Law Coeffs
wetdep( 77, 'NRN12OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 78. Using EtOOH Henry's Law Coeffs
wetdep( 78, 'RA13OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) BENZOOH)
wetdep( 78, 'RA13OOH   ',                                                      &
[2.30e+03,6.00e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 79. Using EtOOH Henry's Law Coeffs
wetdep( 79, 'RA16OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) BENZOOH)
wetdep( 79, 'RA16OOH   ',                                                      &
[2.30e+03,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 80. Using EtOOH Henry's Law Coeffs
wetdep( 80, 'RA19OOH   ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 81. Using EtOOH Henry's Law Coeffs
wetdep( 81, 'RTN28OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep  (Schwantes (2020) hydroperoxy acetone)
wetdep( 81, 'RTN28OOH  ',                                                      &
[1.16e+04,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 82. Using EtOOH Henry's Law Coeffs
wetdep( 82, 'NRTN28OOH ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 83. Using EtOOH Henry's Law Coeffs
wetdep( 83, 'RTN26OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) hydroperoxy acetone)
wetdep( 83, 'RTN26OOH  ',                                                      &
[1.16e+04,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 84. Using EtOOH Henry's Law Coeffs
wetdep( 84, 'RTN25OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) hydroperoxy acetone)
wetdep( 84, 'RTN25OOH  ',                                                      &
[1.16e+04,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 85. Using EtOOH Henry's Law Coeffs
wetdep( 85, 'RTN24OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 86. Using EtOOH Henry's Law Coeffs
wetdep( 86, 'RTN23OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 87. Using EtOOH Henry's Law Coeffs
wetdep( 87, 'RTN14OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) hydroperoxy acetone)
wetdep( 87, 'RTN14OOH  ',                                                      &
[1.16e+04,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 88. Using EtOOH Henry's Law Coeffs 120th
wetdep( 88, 'RTN10OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 89. Using EtOOH Henry's Law Coeffs
wetdep( 89, 'RTX28OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) hydroperoxy acetone)
wetdep( 89, 'RTX28OOH  ',                                                      &
[1.16e+04,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 90. Using EtOOH Henry's Law Coeffs
wetdep( 90, 'RTX24OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) hydroperoxy acetone)
wetdep( 90, 'RTX24OOH  ',                                                      &
[1.16e+04,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 91. Using EtOOH Henry's Law Coeffs 25th
wetdep( 91, 'RTX22OOH  ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) hydroperoxy acetone)
wetdep( 91, 'RTX22OOH  ',                                                      &
[1.16e+04,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 92. Using EtOOH Henry's Law Coeffs
wetdep( 92, 'NRTX28OOH ',                                                      &
[3.40e+02,5.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 93. Using MeOH Henry's Law Coeffs
wetdep( 93, 'AROH14    ',                                                      &
[2.30e+02,4.90e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) PHENOL)
wetdep( 93, 'AROH14    ',                                                      &
[2.84e+03,2.70e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 94. Using MeOH Henry's Law Coeffs
wetdep( 94, 'ARNOH14   ',                                                      &
[2.30e+02,4.90e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Guo and Brimblecombe (2007) - NITROPHENOL)
wetdep( 94, 'ARNOH14   ',                                                      &
[8.5e+01,6.270e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 95. Using MeOH Henry's Law Coeffs
wetdep( 95, 'AROH17    ',                                                      &
[2.30e+02,4.90e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Schwantes (2020) CRESOL)
wetdep( 95, 'AROH17    ',                                                      &
[5.67e+02,5.80e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 96. Using MeOH Henry's Law Coeffs
wetdep( 96, 'ARNOH17   ',                                                      &
[2.30e+02,4.90e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! Updated wet dep (Guo and Brimblecombe (2007) - NITROPHENOL)
wetdep( 96, 'ARNOH17   ',                                                      &
[8.50e+01,6.27e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 97. (Schwantes (2020) TERPA2PAN)
wetdep(  97, 'RTN26PAN  ',                                                     &
[9.59e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 98 (Schwantes (2020) HONITR)
wetdep( 98, 'RN9NO3    ',                                                      &
[2.64e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 99 (Schwantes (2020) HONITR)
wetdep( 99, 'RN12NO3   ',                                                      &
[2.64e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 100 (Schwantes (2020) HONITR)
wetdep(100, 'RN15NO3   ',                                                      &
[2.64e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 101 (Schwantes (2020) HONITR)
wetdep(101, 'RN18NO3   ',                                                      &
[2.64e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 102 (Schwantes (2020) beta isoprene hydroxy nitrate)
wetdep(102, 'RU14NO3   ',                                                      &
[8.34e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 103 (Schwantes (2020) beta isoprene hydroxy nitrate)
wetdep(103, 'RTN28NO3  ',                                                      &
[8.34e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 104 (Schwantes (2020) beta isoprene hydroxy nitrate)
wetdep(104, 'RTN25NO3  ',                                                      &
[8.34e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 105 (Schwantes (2020) beta isoprene hydroxy nitrate)
wetdep(105, 'RTX28NO3  ',                                                      &
[8.34e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 106 (Schwantes (2020) beta isoprene hydroxy nitrate)
wetdep(106, 'RTX22NO3 ',                                                       &
[8.34e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 107. (Schwantes (2020) beta isoprene hydroxy nitrate)
wetdep(107, 'RA13NO3   ',                                                      &
[8.34e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 108. (Schwantes (2020) beta isoprene hydroxy nitrate)
wetdep(108, 'RA16NO3   ',                                                      &
[8.34e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 109. (Schwantes (2020) beta isoprene hydroxy nitrate)
wetdep(109, 'RA19NO3   ',                                                      &
[8.34e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 110 (Schwantes (2020) HONITR)
wetdep(110, 'HOC2H4NO3 ',                                                      &
[2.64e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,107),           &
! WD: 111 (Schwantes (2020) HONITR)
wetdep(111, 'PHAN      ',                                                      &
[2.64e+03,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 112 (Schwantes (2020) IEPOX)
wetdep(112, 'IEPOX     ',                                                      &
[3.0e+07,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),            &
! WD: 113 (Schwantes (2020) ICHE)
wetdep(113, 'HMML      ',                                                      &
[2.09e+06,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 114
wetdep(114, 'HUCARB9   ',                                                      &
[1.1e+05,6.00e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),            &
! WD: 115 (Schwantes (2020) HPALD)
wetdep(115, 'HPUCARB12',                                                       &
[2.30e+05,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 116 (Schwantes (2020) DHPMPAL)
wetdep(116, 'DHPCARB9 ',                                                       &
[9.37e+07,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 117 (Schwantes (2020) DHPMPAL)
wetdep(117, 'DHPR12OOH ',                                                      &
[9.37e+07,6.01e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 118 (Schwantes (2020) hydroxy acetoneL)
wetdep(118, 'DHCARB9   ',                                                      &
[1.46e+03,0.49e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 119 (Schwantes (2020) MVKN)
wetdep(119, 'RU12NO3   ',                                                      &
[1.84e+05,6.00e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119),           &
! WD: 120  (Schwantes (2020) MVKN)
wetdep(120, 'RU10NO3   ',                                                      &
[1.84e+05,6.00e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00],cs,0,0,119)            &
]

! ---------------------------------------------------------------------

! Bimolecular reactions

! Rates taken from  file "UKCA_Reaction_Rates - Bimolecular Reactions.csv"
! on 12:08 2012/05/01

! Too many to define here in one statement.
! Start and end bounds for 1st section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = 1
n_ratb_e = n_ratb_s+55

! Item number, reactants x2, products x4, Arrhenius coefficients (k0, alpha
! and beta), fractional production rates x4, chemistry scheme, qualifier,
! disqualifier, version
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B001 JPL2011
ratb_t1(1,'Br        ','Cl2O2     ','BrCl      ','Cl        ','O2        ',    &
'          ',5.90e-12,  0.00,  170.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(1,'Br        ','Cl2O2     ','BrCl      ','Cl        ','O2        ',    &
'          ',3.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! 3 B002 JPL2011
ratb_t1(2,'Br        ','HCHO      ','HBr       ','CO        ','HO2       ',    &
'          ',1.70e-11,  0.00,  800.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B003 JPL2011
ratb_t1(3,'Br        ','HO2       ','HBr       ','O2        ','          ',    &
'          ',4.80e-12,  0.00,  310.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(3,'Br        ','HO2       ','HBr       ','O2        ','          ',    &
'          ',1.40e-11,  0.00,  590.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B004 JPL2011
ratb_t1(4,'Br        ','O3        ','BrO       ','O2        ','          ',    &
'          ',1.60e-11,  0.00, 780.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),   &
ratb_t1(4,'Br        ','O3        ','BrO       ','O2        ','          ',    &
'          ',1.70e-11,  0.00,  800.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B005 JPL2011
ratb_t1(5,'Br        ','OClO      ','BrO       ','ClO       ','          ',    &
'          ',2.60e-11,  0.00, 1300.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(5,'Br        ','OClO      ','BrO       ','ClO       ','          ',    &
'          ',2.70e-11,  0.00, 1300.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B006 JPL2011
ratb_t1(6,'BrO       ','BrO       ','Br        ','Br        ','O2        ',    &
'          ',2.40e-12,  0.00,  -40.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - do for s, st, & cs
ratb_t1(6,'BrO       ','BrO       ','Br        ','Br        ','O2        ',    &
'          ',1.50e-12,  0.00, -230.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B007 JPL2011
ratb_t1(7,'BrO       ','ClO       ','Br        ','Cl        ','O2        ',    &
'          ',2.30e-12,  0.00, -260.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(7,'BrO       ','ClO       ','Br        ','Cl        ','O2        ',    &
'          ',2.90e-12,  0.00, -220.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B008 JPL2011
ratb_t1(8,'BrO       ','ClO       ','Br        ','OClO      ','          ',    &
'          ',9.50e-13,  0.00, -550.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(8,'BrO       ','ClO       ','Br        ','OClO      ','          ',    &
'          ',1.60e-12,  0.00, -430.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B009 JPL2011
ratb_t1(9,'BrO       ','ClO       ','BrCl      ','O2        ','          ',    &
'          ',4.10e-13,  0.00, -290.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(9,'BrO       ','ClO       ','BrCl      ','O2        ','          ',    &
'          ',5.80e-13,  0.00, -170.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B010 JPL2011
! Reaction essentially removed (see comments) by setting k=0.
! updates for 121 - REMOVE for st & cs
ratb_t1(10,'BrO       ','HO2       ','HBr       ','O3        ','          ',   &
'          ',0.00e+00,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,rm,107), &
ratb_t1(10,'BrO       ','HO2       ','HBr       ','O3        ','          ',   &
'          ',7.00e-14,  0.00, -540.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B011 JPL2011
ratb_t1(11,'BrO       ','HO2       ','HOBr      ','O2        ','          ',   &
'          ',4.50e-12,  0.00, -460.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(11,'BrO       ','HO2       ','HOBr      ','O2        ','          ',   &
'          ',3.30e-12,  0.00, -540.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B012 JPL2011
ratb_t1(12,'BrO       ','NO        ','Br        ','NO2       ','          ',   &
'          ',8.80e-12,  0.00, -260.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(12,'BrO       ','NO        ','Br        ','NO2       ','          ',   &
'          ',8.70e-12,  0.00, -260.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B013 JPL2011
ratb_t1(13,'BrO       ','OH        ','Br        ','HO2       ','          ',   &
'          ',1.70e-11,  0.00, -250.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(13,'BrO       ','OH        ','Br        ','HO2       ','          ',   &
'          ',7.50e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B014 JPL2011
ratb_t1(14,'CF2Cl2    ','O(1D)     ','Cl        ','ClO       ','          ',   &
'          ',1.40e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - do for s, st, & cs
ratb_t1(14,'CF2Cl2    ','O(1D)     ','Cl        ','ClO       ','          ',   &
'          ',1.20e-10,  0.00,  -25.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B015 JPL2011
ratb_t1(15,'CFCl3     ','O(1D)     ','Cl        ','Cl        ','ClO       ',   &
'          ',2.30e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - do for s, st, & cs
ratb_t1(15,'CFCl3     ','O(1D)     ','Cl        ','Cl        ','ClO       ',   &
'          ',2.07e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B016 JPL2011
ratb_t1(16,'Cl        ','CH4       ','HCl       ','MeOO      ','          ',   &
'          ',7.30e-12,  0.00, 1280.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(16,'Cl        ','CH4       ','HCl       ','MeOO      ','          ',   &
'          ',6.60e-12,  0.00, 1240.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! updates for 121 - do for s, st, & cs
ratb_t1(16,'Cl        ','CH4       ','HCl       ','MeOO      ','          ',   &
'          ',7.10e-12,  0.00, 1270.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B017 IUPAC2006
ratb_t1(17,'Cl        ','Cl2O2     ','Cl        ','Cl        ','Cl        ',   &
'          ',7.60e-11,  0.00,  -65.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(17,'Cl        ','Cl2O2     ','Cl        ','Cl        ','Cl        ',   &
'          ',1.00e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B018 JPL2011
ratb_t1(18,'Cl        ','ClONO2    ','Cl        ','Cl        ','NO3       ',   &
'          ',6.50e-12,  0.00, -135.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B019 JPL2011
ratb_t1(19,'Cl        ','H2        ','HCl       ','H         ','          ',   &
'          ',3.05e-11,  0.00, 2270.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(19,'Cl        ','H2        ','HCl       ','H         ','          ',   &
'          ',3.90e-11,  0.00, 2310.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! updates for 121 - align s to st & cs
ratb_t1(19,'Cl        ','H2        ','HCl       ','H         ','          ',   &
'          ',3.05e-11,  0.00, 2270.00, 0.00, 0.00, 0.00, 0.00,s,0,0,121),      &
! B020 JPL2011
ratb_t1(20,'Cl        ','H2O2      ','HCl       ','HO2       ','          ',   &
'          ',1.10e-11,  0.00,  980.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B021 JPL2011
ratb_t1(21,'Cl        ','HCHO      ','HCl       ','CO        ','HO2       ',   &
'          ',8.10e-11,  0.00,   30.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(21,'Cl        ','HCHO      ','HCl       ','CO        ','HO2       ',   &
'          ',8.20e-11,  0.00,   34.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B022 JPL2011
ratb_t1(22,'Cl        ','HO2       ','ClO       ','OH        ','          ',   &
'          ',3.65e-11,  0.00,  375.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(22,'Cl        ','HO2       ','ClO       ','OH        ','          ',   &
'          ',4.10e-11,  0.00,  450.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! updates for 121 - do for s, st, & cs
ratb_t1(22,'Cl        ','HO2       ','ClO       ','OH        ','          ',   &
'          ',3.60e-11,  0.00,  375.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B023 JPL2011
ratb_t1(23,'Cl        ','HO2       ','HCl       ','O2        ','          ',   &
'          ',1.40e-11,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(23,'Cl        ','HO2       ','HCl       ','O2        ','          ',   &
'          ',1.80e-11,  0.00, -170.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B024 JPL2011
ratb_t1(24,'Cl        ','HOCl      ','Cl        ','Cl        ','OH        ',   &
'          ',3.40e-12,  0.00,  130.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(24,'Cl        ','HOCl      ','Cl        ','Cl        ','OH        ',   &
'          ',2.50e-12,  0.00,  130.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B025 JPL2011
ratb_t1(25,'Cl        ','MeOOH     ','HCl       ','MeOO      ','          ',   &
'          ',5.70e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B026 JPL2011
ratb_t1(26,'Cl        ','NO3       ','ClO       ','NO2       ','          ',   &
'          ',2.40e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B027 JPL2011
ratb_t1(27,'Cl        ','O3        ','ClO       ','O2        ','          ',   &
'          ',2.30e-11,  0.00,  200.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(27,'Cl        ','O3        ','ClO       ','O2        ','          ',   &
'          ',2.90e-11,  0.00,  260.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B028 JPL2011
ratb_t1(28,'Cl        ','OClO      ','ClO       ','ClO       ','          ',   &
'          ',3.40e-11,  0.00, -160.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(28,'Cl        ','OClO      ','ClO       ','ClO       ','          ',   &
'          ',3.20e-11,  0.00, -170.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B029 JPL2011
ratb_t1(29,'ClO       ','ClO       ','Cl        ','Cl        ','O2        ',   &
'          ',1.00e-12,  0.00, 1590.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B030 JPL2011
ratb_t1(30,'ClO       ','ClO       ','Cl        ','Cl        ','O2        ',   &
'          ',3.00e-11,  0.00, 2450.00, 0.00, 0.00, 0.00,0.00,s+st+cs,0,0,107)  &
]

! Start and end bounds for 2nd section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+67
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B031 JPL2011
ratb_t1(31,'ClO       ','ClO       ','Cl        ','OClO      ','          ',   &
'          ',3.50e-13,  0.00, 1370.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B032 JPL2011
! Reaction not active - k=0 (see notes)
! updates for 121 - REMOVE for s, st, & cs
ratb_t1(32,'ClO       ','HO2       ','HCl       ','O3        ','          ',   &
'          ',0.00e+00,  0.00,   0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,rm,107),&
ratb_t1(33,'ClO       ','HO2       ','HOCl      ','O2        ','          ',   &
'          ',4.60e-13,  0.00, -710.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B033 JPL2011
ratb_t1(33,'ClO       ','HO2       ','HOCl      ','O2        ','          ',   &
'          ',2.60e-12,  0.00, -290.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
! B034 JPL2011
ratb_t1(34,'ClO       ','MeOO      ','Cl        ','HCHO      ','HO2       ',   &
'          ',3.30e-12,  0.00,  115.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - do for s, st, & cs
ratb_t1(34,'ClO       ','MeOO      ','Cl        ','HCHO      ','HO2       ',   &
'          ',1.80e-12,  0.00,  600.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B035 JPL2011
ratb_t1(35,'ClO       ','NO        ','Cl        ','NO2       ','          ',   &
'          ',6.40e-12,  0.00, -290.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(35,'ClO       ','NO        ','Cl        ','NO2       ','          ',   &
'          ',6.20e-12,  0.00, -295.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B036 JPL2011
ratb_t1(36,'ClO       ','NO3       ','Cl        ','O2        ','NO2       ',   &
'          ',4.60e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - do for s, st, & cs
ratb_t1(36,'ClO       ','NO3       ','Cl        ','O2        ','NO2       ',   &
'          ',4.70e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B037 JPL2011
! Reaction not active - k=0 (see notes
! updates for 121 - REMOVE for s, st, & cs
ratb_t1(37,'ClO       ','NO3       ','OClO      ','NO2       ','          ',   &
'          ',0.00e+00,  0.00,   0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,rm,107),&
! B038 IUPAC2002
ratb_t1(38,'EtCO3     ','NO        ','EtOO      ','CO2       ','NO2       ',   &
'          ',6.70e-12,  0.00, -340.00, 0.00, 0.00, 0.00, 0.0,ti+t+st,0,0,107), &
! CS2
ratb_t1(38,'EtCO3     ','NO        ','EtOO      ','CO2       ','NO2       ',   &
'          ',7.50e-12,  0.00, -290.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! - Different reaction rates in CRI
ratb_t1(38,'EtCO3     ','NO        ','EtOO      ','CO2       ','NO2       ',   &
'          ',8.10e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B039 MCMv3.2
ratb_t1(39,'EtCO3     ','NO3       ','EtOO      ','CO2       ','NO2       ',   &
'          ',4.00e-12,  0.00,    0.00, 0.00, 0.00,0.0, 0.0,ti+t+st+cs,0,0,107),&
! CS2
ratb_t1(39,'EtCO3     ','NO3       ','EtOO      ','CO2       ','NO2       ',   &
'          ',3.68e-12,  0.00,    0.00, 0.00, 0.00,0.0, 0.0,cs,0,0,119),        &
! B040 IUPAC2002
ratb_t1(40,'EtOO      ','MeCO3     ','MeCHO     ','HO2       ','MeOO      ',   &
'          ',4.40e-13,  0.00,-1070.00, 0.00, 0.00, 0.0, 0.0,ti+t+st,0,rp,107), &
! B041 IUPAC2005
ratb_t1(41,'EtOO      ','NO        ','MeCHO     ','HO2       ','NO2       ',   &
'          ',2.55e-12,  0.00, -380.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
! updates for 121 - st only
ratb_t1(41,'EtOO      ','NO        ','MeCHO     ','HO2       ','NO2       ',   &
'          ',2.60e-12,  0.00, -365.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(41,'EtOO      ','NO        ','MeCHO     ','HO2       ','NO2       ',   &
'          ',2.60e-12,  0.00, -380.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
! (Slightly) different reaction rates in CRI
!    Remaining fraction goes to EtONO2
ratb_t1(41,'EtOO      ','NO        ','MeCHO     ','HO2       ','NO2       ',   &
'          ',2.58e-12,  0.00, -365.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(41,'EtOO      ','NO        ','MeCHO     ','HO2       ','NO2       ',   &
'          ',2.53e-12,  0.00, -380.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B042 IUPAC2008
ratb_t1(42,'EtOO      ','NO3       ','MeCHO     ','HO2       ','NO2       ',   &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.0,ti+t+st,0,0,107), &
! (Slightly) different reaction rates in CRI
ratb_t1(42,'EtOO      ','NO3       ','MeCHO     ','HO2       ','NO2       ',   &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
ratb_t1(42,'EtOO      ','NO3       ','MeCHO     ','HO2       ','NO2       ',   &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B043 JPL2011
ratb_t1(43,'H         ','HO2       ','H2        ','O2        ','          ',   &
'          ',6.90e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(43,'H         ','HO2       ','H2        ','O2        ','          ',   &
'          ',2.35e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B044 JPL2011
ratb_t1(44,'H         ','HO2       ','O(3P)     ','H2O       ','          ',   &
'          ',1.62e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B045 JPL2011
ratb_t1(45,'H         ','HO2       ','OH        ','OH        ','          ',   &
'          ',7.20e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(45,'H         ','HO2       ','OH        ','OH        ','          ',   &
'          ',5.59e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B046 JPL2011
ratb_t1(46,'H         ','NO2       ','OH        ','NO        ','          ',   &
'          ',4.00e-10,  0.00,  340.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - do for s, st, & cs
ratb_t1(46,'H         ','NO2       ','OH        ','NO        ','          ',   &
'          ',1.35e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B047 JPL2011
ratb_t1(47,'H         ','O3        ','OH        ','O2        ','          ',   &
'          ',1.40e-10,  0.00,  470.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B048
ratb_t1(48,'HO2       ','EtCO3     ','O2        ','EtCO3H    ','          ',   &
'          ',3.05e-13,  0.00,  -1040.00, 0.00, 0.00, 0.00, 0.00,ti,0,0,107),   &
! B048 MCMv3.2*
ratb_t1(48,'HO2       ','EtCO3     ','O2        ','EtCO3H    ','          ',   &
'          ',4.40e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,st+t,0,0,107),   &
! updates for 121 - st only
ratb_t1(48,'HO2       ','EtCO3     ','O2        ','EtCO3H    ','          ',   &
'          ',1.76e-13,  0.00, -1040.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),    &
! - Different reaction rates, single branch in CRI
ratb_t1(48,'HO2       ','EtCO3     ','O2        ','EtCO3H    ','          ',   &
'          ',4.30e-13,  0.00,-1040.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(48,'HO2       ','EtCO3     ','OH        ','EtCO3H    ','EtOO      ',   &
'          ',5.20e-13,  0.00,-980.00, 0.44, 0.56, 0.44, 0.00,cs,0,0,119),      &
! B049 - NOTE EtCO2H IS *NOT* CURRENTLY A SPECIES CONSIDERED BY ASAD
ratb_t1(49,'HO2       ','EtCO3     ','O3        ','EtCO2H    ','          ',   &
'          ',1.25e-13,  0.00,  -1040.00, 0.00, 0.00, 0.00, 0.00,ti,0,0,107),   &
! B049 MCMv3.2
ratb_t1(49,'HO2       ','EtCO3     ','O3        ','EtCO2H    ','          ',   &
'          ',7.80e-14,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,t+st,0,0,107),   &
! updates for 121 - st only
ratb_t1(49,'HO2       ','EtCO3     ','O3        ','EtCO2H    ','          ',   &
'          ',6.45e-14,  0.00,  -1040.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),   &
! B050 IUPAC2011
ratb_t1(50,'HO2       ','EtOO      ','EtOOH     ','          ','          ',   &
'          ',6.40e-13,  0.00, -710.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
ratb_t1(50,'HO2       ','EtOO      ','EtOOH     ','          ','          ',   &
'          ',3.80e-13,  0.00, -900.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
! - Different reaction rates in CRI
ratb_t1(50,'HO2       ','EtOO      ','EtOOH     ','          ','          ',   &
'          ',7.50e-13,  0.00, -700.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(50,'HO2       ','EtOO      ','EtOOH     ','          ','          ',   &
'          ',4.30e-13,  0.00, -870.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B051 JPL2011 see also asad_bimol
ratb_t1(51,'HO2       ','HO2       ','H2O2      ','          ','          ',   &
'          ',3.00e-13,  0.00, -460.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
ratb_t1(51,'HO2       ','HO2       ','H2O2      ','          ','          ',   &
'          ',2.20e-13,  0.00, -600.00, 0.00,0.00,0.00,0.00,ti+t+ol+cs,0,0,107),&
ratb_t1(51,'HO2       ','HO2       ','H2O2      ','O2        ','          ',   &
'          ',2.20e-13,  0.00, -600.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B052 Poschl00
ratb_t1(52,'HO2       ','ISO2      ','ISOOH     ','          ','          ',   &
'          ',2.05e-13,  0.00,  -1300.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,0,107),&
! B053 Poschl00
ratb_t1(53,'HO2       ','MACRO2    ','MACROOH   ','          ','          ',   &
'          ',1.82e-13,  0.00,  -1300.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,0,107),&
! B054 IUPAC2009
ratb_t1(54,'HO2       ','MeCO3     ','MeCO2H    ','O3        ','          ',   &
'          ',7.80e-14,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
! updates for 121
ratb_t1(54,'HO2       ','MeCO3     ','MeCO2H    ','O3        ','          ',   &
'          ',2.25e-13,  0.00, -730.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(54,'HO2       ','MeCO3     ','MeCO2H    ','O3        ','          ',   &
'          ',1.04e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
! B055 IUPAC2009
ratb_t1(55,'HO2       ','MeCO3     ','MeCO3H    ','          ','          ',   &
'          ',2.13e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
! updates for 121
ratb_t1(55,'HO2       ','MeCO3     ','MeCO3H    ','          ','          ',   &
'          ',6.40e-13,  0.00, -730.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(55,'HO2       ','MeCO3     ','MeCO3H    ','          ','          ',   &
'          ',2.08e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
! - Different rates, single branch in CRI
ratb_t1(55,'HO2       ','MeCO3     ','MeCO3H    ','          ','          ',   &
'          ',4.30e-13,  0.00,-1040.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(55,'HO2       ','MeCO3     ','MeCO3H    ','          ','          ',   &
'          ',2.91e-13,  0.00,-980.00, 0.56, 0.44, 0.44, 0.00,cs,0,0,119),      &
! B056 IUPAC2009
ratb_t1(56,'HO2       ','MeCO3     ','OH        ','MeOO      ','          ',   &
'          ',2.29e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,t+st,0,0,107),   &
! updates for 121 - st only
ratb_t1(56,'HO2       ','MeCO3     ','OH        ','MeOO      ','          ',   &
'          ',8.65e-13,  0.00, -730.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(56,'HO2       ','MeCO3     ','OH        ','MeOO      ','          ',   &
'          ',2.08e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,ti,0,0,107),     &
ratb_t1(56,'HO2       ','MeCO3     ','O3        ','MeOO      ','MeCO3     ',   &
'          ',5.20e-13,  0.00, -980.00, 0.300, 0.800, 0.200, 0.00,r,0,0,107),   &
! CS2
ratb_t1(56,'HO2       ','MeCO3     ','OH        ','MeOO      ','          ',   &
'          ',2.29e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B057 IUPAC2009
ratb_t1(57,'HO2       ','MeCOCH2OO ','MeCOCH2OOH','          ','          ',   &
'          ',9.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
! updates for 121
ratb_t1(57,'HO2       ','MeCOCH2OO ','MeCOCH2OOH','          ','          ',   &
'          ',8.60e-13,  0.00, -700.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(57,'HO2       ','MeCOCH2OO ','MeCOCH2OOH','          ','          ',   &
'          ',1.36e-13,  0.00,  -1250.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107), &
! B058 IUPAC2009 see also asad_bimol
ratb_t1(58,'HO2       ','MeOO      ','HCHO      ','          ','          ',   &
'          ',3.80e-13,  0.00, -780.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! B059 IUPAC2009 see also asad_bimol
! Is branching not handled consistently? No branching in R.
! n.b. also no branching for CRI. Will only use T-dependent
! branching-Ratio if asad_bimol detects both reactions
ratb_t1(59,'HO2       ','MeOO      ','MeOOH     ','          ','          ',   &
'          ',3.80e-13,  0.00, -780.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r+cs,0,0,107)&
]


! Start and end bounds for 3rd section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+71
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
ratb_t1(59,'HO2       ','MeOO      ','O2        ','MeOOH     ','          ',   &
'          ',3.80e-13,  0.00, -780.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B060 JPL2011
ratb_t1(60,'HO2       ','NO        ','OH        ','NO2       ','          ',   &
'          ',3.30e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
! updates for 121 - st only
ratb_t1(60,'HO2       ','NO        ','OH        ','NO2       ','          ',   &
'          ',3.44e-12,  0.00, -260.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(60,'HO2       ','NO        ','OH        ','NO2       ','          ',   &
'          ',3.50e-12,  0.00, -250.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
ratb_t1(60,'HO2       ','NO        ','OH        ','NO2       ','          ',   &
'          ',3.60e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,ti+t+cs,0,0,107),&
! CS2
ratb_t1(60,'HO2       ','NO        ','OH        ','NO2       ','          ',   &
'          ',3.45e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B060a added Alex
! Not in TI, S, ST and R schemes
ratb_t1(60,'HO2       ','NO        ','HONO2     ','          ','          ',   &
'          ',3.60e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,t,0,0,107),      &
! B061 JPL2011
ratb_t1(61,'HO2       ','NO3       ','OH        ','NO2       ','          ',   &
'          ',3.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
ratb_t1(61,'HO2       ','NO3       ','OH        ','NO2       ','          ',   &
'          ',4.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t+cs,0,0,107),&
ratb_t1(61,'HO2       ','NO3       ','OH        ','NO2       ','O2        ',   &
'          ',2.15e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B062 IUPAC2001
ratb_t1(62,'HO2       ','O3        ','OH        ','O2        ','          ',   &
'          ',2.03e-16,  4.57, -693.00, 0.0, 0.0,0.0,0.0,ti+t+st+r+cs,0,0,107), &
! updates for 121 - st only
ratb_t1(62,'HO2       ','O3        ','OH        ','O2        ','          ',   &
'          ',1.00e-14,  0.00,  490.00, 0.0, 0.0,0.0,0.0,st,0,0,121),           &
ratb_t1(62,'HO2       ','O3        ','OH        ','O2        ','O2        ',   &
'          ',2.03e-16,  4.57, -693.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B063 MCMv3.2
ratb_t1(63,'HO2       ','i-PrOO    ','i-PrOOH   ','          ','          ',   &
'          ',1.51e-13,  0.00,-1300.00, 0.0, 0.0, 0.0,0.0,ti+t+st+r+cs,0,0,107),&
! B064 MCMv3.2
ratb_t1(64,'HO2       ','n-PrOO    ','n-PrOOH   ','          ','          ',   &
'          ',1.51e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.0,ti+t+st,0,0,107), &
! B065 Poschl00
ratb_t1(65,'ISO2      ','ISO2      ','MACR      ','MACR      ','HCHO      ',   &
'HO2       ',2.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,rp,107), &
! B066 Poschl00
ratb_t1(66,'MACRO2    ','MACRO2    ','HACET     ','MGLY      ','HCHO      ',   &
'CO        ',1.00e-12,  0.00,    0.00, 2.00, 2.00, 1.00, 1.00,ti+st,0,rp,107), &
! B067 Poschl00
ratb_t1(67,'MACRO2    ','MACRO2    ','HO2       ','          ','          ',   &
'          ',1.00e-12,  0.00,    0.00, 2.00, 0.00, 0.00, 0.00,ti+st,0,rp,107), &
! B068 JPL2011
ratb_t1(68,'MeBr      ','Cl        ','Br        ','HCl       ','          ',   &
'          ',1.40e-11,  0.00, 1030.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
! updates for 121 - st & cs
ratb_t1(68,'MeBr      ','Cl        ','Br        ','HCl       ','          ',   &
'          ',1.46e-11,  0.00, 1040.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,121),  &
ratb_t1(68,'MeBr      ','Cl        ','Br        ','HCl       ','          ',   &
'          ',1.70e-11,  0.00, 1080.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B069 JPL2011
ratb_t1(69,'MeBr      ','O(1D)     ','Br        ','OH        ','          ',   &
'          ',1.80e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B070 JPL2011
ratb_t1(70,'MeBr      ','OH        ','Br        ','H2O       ','          ',   &
'          ',2.35e-12,  0.00, 1300.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(70,'MeBr      ','OH        ','Br        ','H2O       ','          ',   &
'          ',1.42e-12,  0.00, 1150.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B071 IUPAC2002
ratb_t1(71,'MeCO3     ','NO        ','MeOO      ','CO2       ','NO2       ',   &
'          ',7.50e-12,  0.00, -290.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r,0,0,107),  &
! updates for 121 - st only
ratb_t1(71,'MeCO3     ','NO        ','MeOO      ','CO2       ','NO2       ',   &
'          ',8.10e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
! (Slightly) different reaction rates in CRI
ratb_t1(71,'MeCO3     ','NO        ','MeOO      ','CO2       ','NO2       ',   &
'          ',8.10e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(71,'MeCO3     ','NO        ','MeOO      ','CO2       ','NO2       ',   &
'          ',7.50e-12,  0.00, -290.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B072 IUPAC2008
! Missing in R
ratb_t1(72,'MeCO3     ','NO3       ','MeOO      ','CO2       ','NO2       ',   &
'          ',4.00e-12,  0.00,    0.00, 0.00,0.00,0.00,0.00,ti+t+st+cs,0,0,107),&
! CS2
ratb_t1(72,'MeCO3     ','NO3       ','MeOO      ','CO2       ','NO2       ',   &
'          ',3.68e-12,  0.00,    0.00, 0.00,0.00,0.00,0.00,     cs,0,0,119),   &
! B073 MCMv3.2
ratb_t1(73,'MeCOCH2OO ','NO        ','MeCO3     ','HCHO      ','NO2       ',   &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
! update to 121 - st only
ratb_t1(73,'MeCOCH2OO ','NO        ','MeCO3     ','HCHO      ','NO2       ',   &
'          ',2.90e-12,  0.00, -300.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(73,'MeCOCH2OO ','NO        ','MeCO3     ','HCHO      ','NO2       ',   &
'          ',2.80e-12,  0.00, -300.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
! B074 MCMv3.2
ratb_t1(74,'MeCOCH2OO ','NO3       ','MeCO3     ','HCHO      ','NO2       ',   &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
ratb_t1(74,'MeCOCH2OO ','NO3       ','MeCO3     ','HCHO      ','NO2       ',   &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
! B075 IUPAC2002
ratb_t1(75,'MeOO      ','MeCO3     ','HO2       ','HCHO      ','MeOO      ',   &
'          ',1.80e-12,  0.00, -500.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r,0,rp,107), &
! B076 IUPAC2002
ratb_t1(76,'MeOO      ','MeCO3     ','MeCO2H    ','HCHO      ','          ',   &
'          ',2.00e-13,  0.00, -500.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r,0,rp,107), &
! B077 IUPAC2002 see also asad_bimol
ratb_t1(77,'MeOO      ','MeOO      ','HO2       ','HO2       ','HCHO      ',   &
'HCHO      ',1.03e-13,  0.00, -365.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r,0,rp,107), &
! updates for 121 - st only
ratb_t1(77,'MeOO      ','MeOO      ','HO2       ','HO2       ','HCHO      ',   &
'HCHO      ',1.03e-13,  0.00, -390.00, 0.0, 0.0, 0.0, 0.0,st,0,rp,121),        &
ratb_t1(77,'MeOO      ','MeOO      ','HCHO      ','HO2       ','O2        ',   &
'          ',9.50e-14,  0.00, -390.00, 2.000, 2.000, 1.000, 0.00,s,0,0,107),   &
! B078 IUPAC2002 see also asad_bimol
ratb_t1(78,'MeOO      ','MeOO      ','MeOH      ','HCHO      ','          ',   &
'          ',1.03e-13,  0.00, -365.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r,0,rp,107), &
! updates for 121 - st only
ratb_t1(78,'MeOO      ','MeOO      ','MeOH      ','HCHO      ','          ',   &
'          ',1.03e-13,  0.00, -390.00, 0.0, 0.0, 0.0, 0.0,st,0,rp,121),        &
! B079 IUPAC2005
ratb_t1(79,'MeOO      ','NO        ','HO2       ','HCHO      ','NO2       ',   &
'          ',2.30e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
! updates for 121 - st only
ratb_t1(79,'MeOO      ','NO        ','HO2       ','HCHO      ','NO2       ',   &
'          ',2.80e-12,  0.00, -300.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(79,'MeOO      ','NO        ','HO2       ','HCHO      ','NO2       ',   &
'          ',2.95e-12,  0.00, -285.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
ratb_t1(79,'MeOO      ','NO        ','HCHO      ','HO2       ','NO2       ',   &
'          ',2.95e-12,  0.00, -285.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! (Slightly) different reaction rates in CRI
ratb_t1(79,'MeOO      ','NO        ','HCHO      ','HO2       ','NO2       ',   &
'          ',3.00e-12,  0.00, -280.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(79,'MeOO      ','NO        ','HCHO      ','HO2       ','NO2       ',   &
'          ',2.298e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B080 IUPAC2005
ratb_t1(80,'MeOO      ','NO        ','MeONO2    ','          ','          ',   &
'          ',2.30e-15,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
! updates for 121
ratb_t1(80,'MeOO      ','NO        ','MeONO2    ','          ','          ',   &
'          ',8.97e-15,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(80,'MeOO      ','NO        ','MeONO2    ','          ','          ',   &
'          ',2.95e-15,  0.00, -285.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
! (Slightly) different reaction rates in CRI
ratb_t1(80,'MeOO      ','NO        ','MeONO2    ','          ','          ',   &
'          ',3.00e-15,  0.00, -280.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B080 CS2
ratb_t1(80,'MeOO      ','NO        ','MeONO2    ','          ','          ',   &
'          ',2.30e-15,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B081 MCMv3.2
ratb_t1(81,'MeOO      ','NO3       ','HO2       ','HCHO      ','NO2       ',   &
'          ',  1.20e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),   &
ratb_t1(81,'MeOO      ','NO3       ','HO2       ','HCHO      ','NO2       ',   &
'          ',1.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
! - Different reaction rates in CRI
ratb_t1(81,'MeOO      ','NO3       ','HO2       ','HCHO      ','NO2       ',   &
'          ',1.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(81,'MeOO      ','NO3       ','HO2       ','HCHO      ','NO2       ',   &
'          ',1.20e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B082 JPL2011
ratb_t1(82,'N         ','NO        ','N2        ','O(3P)     ','          ',   &
'          ',2.10e-11,  0.00, -100.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B083 JPL2011
ratb_t1(83,'N         ','NO2       ','N2O       ','O(3P)     ','          ',   &
'          ',5.80e-12,  0.00, -220.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B084 JPL2011
ratb_t1(84,'N         ','O2        ','NO        ','O(3P)     ','          ',   &
'          ',1.50e-11,  0.00, 3600.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(84,'N         ','O2        ','NO        ','O(3P)     ','          ',   &
'          ',3.30e-12,  0.00, 3150.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B085 ????
! Not in strat scheme
! Added N2O5+H2O gasphase reaction to CRI-Strat temporarily to account for loss
! of N2O5, at least until tropospheric heterogeneous chemistry is working...
ratb_t1(85,'N2O5      ','H2O       ','HONO2     ','HONO2     ','          ',   &
 '          ',2.50e-22, 0.00,0.00,0.00,0.00,0.00,0.00,ti+t+st+cs,0,0,107),     &
! B086 Poschl00
ratb_t1(86,'NO        ','ISO2      ','ISON      ','          ','          ',   &
'          ',1.12e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,0,107),  &
! updates for 121 - st only
ratb_t1(86,'NO        ','ISO2      ','ISON      ','          ','          ',   &
'          ',9.68e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
! B087 Poschl00
ratb_t1(87,'NO        ','ISO2      ','NO2       ','MACR      ','HCHO      ',   &
'HO2       ',2.43e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,0,107),  &
! updates for 121 - st only
ratb_t1(87,'NO        ','ISO2      ','NO2       ','MACR      ','HCHO      ',   &
'HO2       ',7.83e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
! B088 Poschl00
ratb_t1(88,'NO        ','MACRO2    ','MGLY      ','HCHO      ','HO2       ',   &
'          ',1.27e-12,  0.00, -360.00, 1.00, 1.50, 1.50, 0.00,ti+st,0,0,107),  &
! updates for 121 - st only
ratb_t1(88,'NO        ','MACRO2    ','MGLY      ','HCHO      ','HO2       ',   &
'          ',1.35e-12,  0.00, -360.00, 1.00, 1.50, 1.50, 0.00,st,0,0,121),     &
! B089 Poschl00
ratb_t1(89,'NO        ','MACRO2    ','NO2       ','MeCO3     ','HACET     ',   &
'CO        ',1.27e-12,  0.00, -360.00, 2.00, 0.50, 0.50, 0.50,ti+st,0,0,107),  &
! updates for 121 - st only
ratb_t1(89,'NO        ','MACRO2    ','NO2       ','MeCO3     ','HACET     ',   &
'CO        ',1.25e-12,  0.00, -360.00, 2.00, 0.50, 0.50, 0.50,st,0,0,121),     &
! B090 JPL2011
ratb_t1(90,'NO        ','NO3       ','NO2       ','NO2       ','          ',   &
'          ',1.50e-11,  0.00, -170.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
! updates for 121 - st only
ratb_t1(90,'NO        ','NO3       ','NO2       ','NO2       ','          ',   &
'          ',1.70e-11,  0.00, -125.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121)      &
]


! Start and end bounds for 4th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+68
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
ratb_t1(90,'NO        ','NO3       ','NO2       ','NO2       ','          ',   &
'          ',1.80e-11,  0.00, -110.00, 0.00, 0.00,0.00,0.00,ti+t+s+cs,0,0,107),&
! B091 JPL2011
ratb_t1(91,'NO        ','O3        ','NO2       ','          ','          ',   &
'          ',3.00e-12,  0.00, 1500.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
ratb_t1(91,'NO        ','O3        ','NO2       ','O2        ','          ',   &
'          ',1.40e-12,  0.00, 1310.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
ratb_t1(91,'NO        ','O3        ','NO2       ','          ','          ',   &
'          ',1.40e-12,  0.00, 1310.00, 0.00, 0.00, 0.00, 0.00,ti+t+cs,0,0,107),&
! B092 JPL2011
ratb_t1(92,'NO2       ','NO3       ','NO        ','NO2       ','O2        ',   &
'          ',4.50e-14,  0.00, 1260.00, 0.0, 0.0, 0.0, 0.0,s+t+st+r+cs,0,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(92,'NO2       ','NO3       ','NO        ','NO2       ','O2        ',   &
'          ',4.35e-14,  0.00, 1340.00, 0.0, 0.0, 0.0, 0.0,s+st+cs,0,0,121),    &
! B093 JPL2011
ratb_t1(93,'NO2       ','O3        ','NO3       ','          ','          ',   &
'          ',1.20e-13,  0.00, 2450.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
ratb_t1(93,'NO2       ','O3        ','NO3       ','O2        ','          ',   &
'          ',1.40e-13,  0.00, 2470.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
ratb_t1(93,'NO2       ','O3        ','NO3       ','          ','          ',   &
'          ',1.40e-13,  0.00, 2470.00, 0.00, 0.00, 0.00, 0.00,ti+t+cs,0,0,107),&
! B094 JPL2011
ratb_t1(94,'NO3       ','Br        ','BrO       ','NO2       ','          ',   &
'          ',1.60e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B095 IUPAC2007
ratb_t1(95,'NO3       ','C5H8      ','ISON      ','          ','          ',   &
'          ',3.15e-12,  0.00,  450.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,0,107),  &
! updates for 121 - st only
ratb_t1(95,'NO3       ','C5H8      ','ISON      ','          ','          ',   &
'          ',2.95e-12,  0.00,  450.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
! '(NO3)C4H6CHO' is too long for the ratb_t1 constructor (12 characters
! instead of 10). We truncate `prod1` to be '(NO3)C4H6C'.
ratb_t1(95,'NO3       ','C5H8      ','(NO3)C4H6C','HO2       ','          ',   &
'          ',3.03e-12,  0.00,  446.00, 0.00, 0.00, 0.00, 0.00,r,0,0,107),      &
! B095 - Different rates and prods in CRI
ratb_t1(95,'NO3       ','C5H8      ','NRU14O2   ','          ','          ',   &
'          ',3.03e-12,  0.00,  446.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(95,'NO3       ','C5H8      ','NRU14O2   ','          ','          ',   &
'          ',3.15e-12,  0.00,  450.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B096 IUPAC2007
ratb_t1(96,'NO3       ','EtCHO     ','HONO2     ','EtCO3     ','          ',   &
'          ',6.30e-15,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
ratb_t1(96,'NO3       ','EtCHO     ','HONO2     ','EtCO3     ','          ',   &
'          ',3.46e-12,  0.00, 1862.00, 0.00, 0.00, 0.00, 0.00,ti+t+cs,0,0,107),&
! CS2
ratb_t1(96,'NO3       ','EtCHO     ','HONO2     ','EtCO3     ','          ',   &
'          ',3.24e-12,  0.00, 1860.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B097 IUPAC2007
ratb_t1(97,'NO3       ','HCHO      ','HONO2     ','HO2       ','CO        ',   &
'          ',2.00e-12,  0.00, 2440.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r,0,0,107),  &
! updates for 121 - st only
ratb_t1(97,'NO3       ','HCHO      ','HONO2     ','HO2       ','CO        ',   &
'          ',5.80e-16,  0.00,    0.00, 0.0, 0.0, 0.0, 0.0,st,0,0,121),         &
ratb_t1(97,'NO3       ','HCHO      ','HONO2     ','CO        ','HO2       ',   &
'          ',5.60e-16,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! Different rate in CRI
ratb_t1(97,'NO3       ','HCHO      ','HONO2     ','HO2       ','CO        ',   &
'          ',5.80e-16,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(97,'NO3       ','HCHO      ','HONO2     ','HO2       ','CO        ',   &
'          ',5.50e-16,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B098 MCMv3.2
ratb_t1(98,'NO3       ','MGLY      ','MeCO3     ','CO        ','HONO2     ',   &
'          ',3.36e-12,  0.00, 1860.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
! updates for 121
ratb_t1(98,'NO3       ','MGLY      ','MeCO3     ','CO        ','HONO2     ',   &
'          ',5.00e-16,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(98,'NO3       ','MGLY      ','MeCO3     ','CO        ','HONO2     ',   &
'          ',3.46e-12,  0.00, 1860.00, 0.00, 0.00, 0.00, 0.00,ti,0,0,107),     &
! B099 IUPAC2007
ratb_t1(99,'NO3       ','Me2CO     ','HONO2     ','MeCOCH2OO ','          ',   &
'          ',3.00e-17,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! B100 IUPAC2007
ratb_t1(100,'NO3       ','MeCHO     ','HONO2     ','MeCO3     ','          ',  &
'          ',1.40e-12,  0.00, 1860.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r,0,0,107),  &
! (Slightly) different rates in CRI - could be rounding error?
ratb_t1(100,'NO3       ','MeCHO     ','HONO2     ','MeCO3     ','          ',  &
'          ',1.44e-12,  0.00, 1862.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(100,'NO3       ','MeCHO     ','HONO2     ','MeCO3     ','          ',  &
 '          ',1.40e-12,  0.00, 1860.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B101 JPL2011
ratb_t1(101,'O(1D)     ','CH4       ','HCHO      ','H2        ','          ',  &
'          ',9.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
! updates for 121 - st & cs
ratb_t1(101,'O(1D)     ','CH4       ','HCHO      ','H2        ','          ',  &
'          ',8.75e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,121),  &
ratb_t1(101,'O(1D)     ','CH4       ','HCHO      ','H2        ','          ',  &
'          ',7.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
ratb_t1(101,'O(1D)     ','CH4       ','HCHO      ','H2        ','          ',  &
'          ',1.50e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B102 JPL2011
ratb_t1(102,'O(1D)     ','CH4       ','HCHO      ','HO2       ','HO2       ',  &
'          ',3.45e-11,  0.00,    0.00, 0.00,0.00,0.00,0.00,ti+t+st+cs,0,0,107),&
! updates for 121 - st & cs only
ratb_t1(102,'O(1D)     ','CH4       ','HCHO      ','HO2       ','HO2       ',  &
'          ',3.50e-11,  0.00,    0.00, 0.00,0.00,0.00,0.00,st+cs,0,0,121),     &
! B103 JPL2011
ratb_t1(103,'O(1D)     ','CH4       ','OH        ','MeOO      ','          ',  &
'          ',1.31e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(103,'O(1D)     ','CH4       ','OH        ','MeOO      ','          ',  &
'          ',1.35e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
ratb_t1(103,'O(1D)     ','CH4       ','OH        ','MeOO      ','          ',  &
'          ',1.05e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
! B104 JPL2011
! Not in TI scheme
ratb_t1(104,'O(1D)     ','CO2       ','O(3P)     ','CO2       ','          ',  &
'          ',7.50e-11,  0.00, -115.00, 0.00, 0.00, 0.00, 0.00,t+st+cs,0,0,107),&
ratb_t1(104,'O(1D)     ','CO2       ','O(3P)     ','CO2       ','          ',  &
'          ',7.40e-11,  0.00, -120.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B105 IUPAC2008
! Not in TI scheme
ratb_t1(105,'O(1D)     ','H2        ','OH        ','H         ','          ',  &
'          ',1.20e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,t+st+cs,0,0,107),&
ratb_t1(105,'O(1D)     ','H2        ','OH        ','H         ','          ',  &
'          ',1.10e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B106 JPL2011
ratb_t1(106,'O(1D)     ','H2O       ','OH        ','OH        ','          ',  &
'          ',1.63e-10,  0.00,  -60.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
ratb_t1(106,'O(1D)     ','H2O       ','OH        ','OH        ','          ',  &
'          ',2.20e-10,  0.00,    0.00, 0.00, 0.00,0.00,0.00,ti+s+t+cs,0,0,107),&
! CS2
ratb_t1(106,'O(1D)     ','H2O       ','OH        ','OH        ','          ',  &
'          ',2.14e-10,  0.00,    0.00, 0.00, 0.00,0.00,0.00,cs,0,0,119),       &
! B107 JPL2011
ratb_t1(107,'O(1D)     ','HBr       ','HBr       ','O(3P)     ','          ',  &
'          ',3.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B108 JPL2011
ratb_t1(108,'O(1D)     ','HBr       ','OH        ','Br        ','          ',  &
'          ',1.20e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(108,'O(1D)     ','HBr       ','OH        ','Br        ','          ',  &
'          ',9.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B109 JPL2011
ratb_t1(109,'O(1D)     ','HCl       ','H         ','ClO       ','          ',  &
'          ',3.60e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(109,'O(1D)     ','HCl       ','H         ','ClO       ','          ',  &
'          ',3.30e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B110 JPL2011
ratb_t1(110,'O(1D)     ','HCl       ','O(3P)     ','HCl       ','          ',  &
'          ',1.35e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(110,'O(1D)     ','HCl       ','O(3P)     ','HCl       ','          ',  &
'          ',1.80e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B111 JPL2011
ratb_t1(111,'O(1D)     ','HCl       ','OH        ','Cl        ','          ',  &
'          ',1.01e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(111,'O(1D)     ','HCl       ','OH        ','Cl        ','          ',  &
'          ',9.90e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,121),&
! B112 JPL2011
ratb_t1(112,'O(1D)     ','N2        ','O(3P)     ','N2        ','          ',  &
'          ',2.15e-11,  0.00, -110.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
ratb_t1(112,'O(1D)     ','N2        ','O(3P)     ','N2        ','          ',  &
'          ',2.10e-11,  0.00, -115.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
ratb_t1(112,'O(1D)     ','N2        ','O(3P)     ','N2        ','          ',  &
'          ',1.80e-11,  0.00, -110.00, 0.00, 0.00, 0.00, 0.00,s+cs,0,0,107),   &
! CS2
ratb_t1(112,'O(1D)     ','N2        ','O(3P)     ','N2        ','          ',  &
'          ',2.00e-11,  0.00, -130.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B113 JPL2011
ratb_t1(113,'O(1D)     ','N2O       ','N2        ','O2        ','          ',  &
'          ',4.60e-11,  0.00,  -20.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
! updates for 121 - st & cs
ratb_t1(113,'O(1D)     ','N2O       ','N2        ','O2        ','          ',  &
'          ',4.64e-11,  0.00,  -20.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,121),  &
ratb_t1(113,'O(1D)     ','N2O       ','N2        ','O2        ','          ',  &
'          ',4.90e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B114 JPL2011
ratb_t1(114,'O(1D)     ','N2O       ','NO        ','NO        ','          ',  &
'          ',7.30e-11,  0.00,  -20.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
! updates for 121 - st & cs
ratb_t1(114,'O(1D)     ','N2O       ','NO        ','NO        ','          ',  &
'          ',7.26e-11,  0.00,  -20.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,121),  &
ratb_t1(114,'O(1D)     ','N2O       ','NO        ','NO        ','          ',  &
'          ',6.70e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B115 JPL2011
ratb_t1(115,'O(1D)     ','O2        ','O(3P)     ','O2        ','          ',  &
'          ',3.30e-11,  0.00,  -55.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
ratb_t1(115,'O(1D)     ','O2        ','O(3P)     ','O2        ','          ',  &
'          ',3.20e-11,  0.00,  -67.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
ratb_t1(115,'O(1D)     ','O2        ','O(3P)     ','O2        ','          ',  &
'          ',3.20e-11,  0.00,  -70.00, 0.00, 0.00, 0.00, 0.00,s+cs,0,0,107),   &
! B115 CS2
ratb_t1(115,'O(1D)     ','O2        ','O(3P)     ','O2        ','          ',  &
'          ',3.20e-11,  0.00,  -67.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119)      &
]


! Start and end bounds for 5th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+66
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B116 JPL2011
! Not in TI scheme (or CRI)
ratb_t1(116,'O(1D)     ','O3        ','O2        ','O(3P)     ','O(3P)     ',  &
'          ',1.20e-10,  0.00,    0.00, 0.00, 0.00,0.00,0.00,s+t+st+cs,0,0,107),&
! B117 JPL2011
! Not in TI scheme
ratb_t1(117,'O(1D)     ','O3        ','O2        ','O2        ','          ',  &
'          ',1.20e-10,  0.00,    0.00, 0.00, 0.00,0.00,0.00,s+t+st+cs,0,0,107),&
! B118 JPL2011
ratb_t1(118,'O(3P)     ','BrO       ','O2        ','Br        ','          ',  &
'          ',1.90e-11,  0.00, -230.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B119 JPL2011
ratb_t1(119,'O(3P)     ','ClO       ','Cl        ','O2        ','          ',  &
'          ',2.80e-11,  0.00,  -85.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(119,'O(3P)     ','ClO       ','Cl        ','O2        ','          ',  &
'          ',3.80e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B120 JPL2011
ratb_t1(120,'O(3P)     ','ClONO2    ','ClO       ','NO3       ','          ',  &
'          ',3.60e-12,  0.00,  840.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(120,'O(3P)     ','ClONO2    ','ClO       ','NO3       ','          ',  &
'          ',2.90e-12,  0.00,  800.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B121 ????
! Not in TI scheme
ratb_t1(121,'O(3P)     ','H2        ','OH        ','HO2       ','          ',  &
'          ',  9.00e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,t,0,0,107),    &
! updates for 121 - REMOVE for s, st, & cs
ratb_t1(121,'O(3P)     ','H2        ','OH        ','H         ','          ',  &
'          ',9.00e-18,  0.00,   0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,rm,107),&
! B122 JPL2011
! Not in T/TI/R schemes
ratb_t1(122,'O(3P)     ','H2O2      ','OH        ','HO2       ','          ',  &
'          ',1.40e-12,  0.00, 2000.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B123 JPL2011
ratb_t1(123,'O(3P)     ','HBr       ','OH        ','Br        ','          ',  &
'          ',5.80e-12,  0.00, 1500.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B124 JPL2011
! Not in TI/R schemes
ratb_t1(124,'O(3P)     ','HCHO      ','OH        ','CO        ','HO2       ',  &
'          ',3.40e-11,  0.00, 1600.00, 0.00, 0.00,0.00,0.00,t,0,0,107),        &
! updates for 121 - REMOVE for s, st, & cs
ratb_t1(124,'O(3P)     ','HCHO      ','OH        ','CO        ','HO2       ',  &
'          ',3.40e-11,  0.00, 1600.00, 0.00, 0.00,0.00,0.00,s+st+cs,0,rm,107), &
! B125 JPL2011
ratb_t1(125,'O(3P)     ','HCl       ','OH        ','Cl        ','          ',  &
'          ',1.00e-11,  0.00, 3300.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B126 IUPAC2001
! Not in TI/R schemes
ratb_t1(126,'O(3P)     ','HO2       ','OH        ','O2        ','          ',  &
'          ',2.70e-11,  0.00, -224.00, 0.00, 0.00,0.00,0.00,s+t+st+cs,0,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(126,'O(3P)     ','HO2       ','OH        ','O2        ','          ',  &
'          ',3.00e-11,  0.00, -200.00, 0.00, 0.00,0.00,0.00,s+st+cs,0,0,121),  &
! B127 JPL2011
ratb_t1(127,'O(3P)     ','HOCl      ','OH        ','ClO       ','          ',  &
'          ',1.70e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B128 JPL2011
ratb_t1(128,'O(3P)     ','NO2       ','NO        ','O2        ','          ',  &
'          ',5.10e-12,  0.00, -210.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
! updates for 121 - st only
ratb_t1(128,'O(3P)     ','NO2       ','NO        ','O2        ','          ',  &
'          ',5.30e-12,  0.00, -200.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(128,'O(3P)     ','NO2       ','O2        ','NO        ','          ',  &
'          ',5.60e-12,  0.00, -180.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
ratb_t1(128,'O(3P)     ','NO2       ','NO        ','O2        ','          ',  &
'          ',5.50e-12,  0.00, -188.00, 0.00, 0.00, 0.00, 0.00,ti+t+cs,0,0,107),&
! B129 IUPAC2009
! Not in TI/R schemes
ratb_t1(129,'O(3P)     ','NO3       ','O2        ','NO2       ','          ',  &
'          ',1.70e-11,  0.00,    0.00, 0.00, 0.00,0.00,0.00,s+t+st+cs,0,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(129,'O(3P)     ','NO3       ','O2        ','NO2       ','          ',  &
'          ',1.30e-11,  0.00,    0.00, 0.00, 0.00,0.00,0.00,s+st+cs,0,0,121),  &
! B130 JPL2011
ratb_t1(130,'O(3P)     ','O3        ','O2        ','O2        ','          ',  &
'          ',8.00e-12,  0.00, 2060.00, 0.0, 0.0, 0.0,0.0,ti+s+t+st+cs,0,0,107),&
! B131 JPL2011
ratb_t1(131,'O(3P)     ','OClO      ','O2        ','ClO       ','          ',  &
'          ',2.40e-12,  0.00,  960.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B132 JPL2011
! Not in TI/R scheme
ratb_t1(132,'O(3P)     ','OH        ','O2        ','HO2       ','          ',  &
'          ',1.80e-11,  0.00, -180.00, 0.00, 0.00, 0.00, 0.00,t,0,0,107),      &
ratb_t1(132,'O(3P)     ','OH        ','O2        ','H         ','          ',  &
'          ',1.80e-11,  0.00, -180.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(132,'O(3P)     ','OH        ','O2        ','H         ','          ',  &
'          ',2.40e-11,  0.00, -110.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B133 IUPAC2007*
ratb_t1(133,'O3        ','C5H8      ','HO2       ','OH        ','          ',  &
'          ',3.33e-15,  0.00, 1995.00, 0.75, 0.75, 0.00, 0.00,ti+st,0,0,107),  &
! updates for 121 - st only
ratb_t1(133,'O3        ','C5H8      ','HO2       ','OH        ','          ',  &
'          ',3.30e-15,  0.00, 2000.00, 0.75, 0.75, 0.00, 0.00,st,0,0,121),     &
! Number mislabelled - renamed to 133!
ratb_t1(133,'O3        ','C5H8      ','OH        ','          ','          ',  &
'          ',3.93e-15,  0.00, 1913.00, 0.54, 0.00, 0.00, 0.00,r,0,0,107),      &
! B134 IUPAC2007*
ratb_t1(134,'O3        ','C5H8      ','MACR      ','HCHO      ','MACRO2    ',  &
'MeCO3     ',3.33e-15,  0.00, 1995.00, 1.95, 1.74, 0.30, 0.30,ti+st,0,0,107),  &
! updates for 121 - st only
ratb_t1(134,'O3        ','C5H8      ','MACR      ','HCHO      ','MACRO2    ',  &
'MeCO3     ',3.30e-15,  0.00, 2000.00, 1.95, 1.74, 0.30, 0.30,st,0,0,121),     &
! B134 - Different rates and prods in CRI
ratb_t1(134,'O3        ','C5H8      ','UCARB10   ','CO        ','HO2       ',  &
'OH        ',2.12e-15,  0.00, 1913.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(134,'O3        ','C5H8      ','UCARB10   ','CO        ','HO2       ',  &
'OH        ',1.288e-15,  0.00, 1995.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B135 IUPAC2007*
ratb_t1(135,'O3        ','C5H8      ','MeOO      ','HCOOH     ','CO        ',  &
'H2O2      ',3.33e-15,  0.00, 1995.00, 0.24, 0.84, 0.42, 0.27,ti+st,0,0,107),  &
! updates for 121 - st only
ratb_t1(135,'O3        ','C5H8      ','MeOO      ','HCOOH     ','CO        ',  &
'H2O2      ',3.30e-15,  0.00, 2000.00, 0.24, 0.84, 0.42, 0.27,st,0,0,121),     &
! Don't know how to map CH2O2.
ratb_t1(135,'O3        ','C5H8      ','MVK       ','CO        ','CH2O2     ',  &
'HO2       ',3.93e-15,  0.00, 1913.00, 2.00, 1.56, 0.44, 0.54,r,0,0,107),      &
! B135 - Different rates and prods in CRI
ratb_t1(135,'O3        ','C5H8      ','UCARB10   ','HCOOH     ','          ',  &
'          ',5.74e-15,  0.00, 1913.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(135,'O3        ','C5H8      ','HCHO      ','MeOO      ','CO        ',  &
'HO2       ',9.79e-16,  0.00, 1995, 2.00, 1.00,1.00 ,1.00, cs,0,0,119),        &
! B136 IUPAC2007* - IDENTICAL REACTANTS/PRODUCTS TO B137
ratb_t1(136,'O3        ','MACR      ','MGLY      ','HCOOH     ','HO2       ',  &
'CO        ',2.13e-16,  0.00, 1520.00, 1.80, 0.90, 0.64, 0.44,ti+st,0,0,107),  &
! updates for 121 - st only
ratb_t1(136,'O3        ','MACR      ','MGLY      ','HCOOH     ','HO2       ',  &
'CO        ',1.23e-15,  0.00, 2100.00, 1.80, 0.90, 0.64, 0.44,st,0,0,121),     &
! B137 IUPAC2007* - IDENTICAL REACTANTS/PRODUCTS TO B136
ratb_t1(137,'O3        ','MACR      ','MGLY      ','HCOOH     ','HO2       ',  &
'CO        ',3.50e-16,  0.00, 2100.00, 1.80, 0.90, 0.64, 0.44,ti+st,0,0,107),  &
! B138 IUPAC2007*
ratb_t1(138,'O3        ','MACR      ','OH        ','MeCO3     ','          ',  &
'          ',2.13e-16,  0.00, 1520.00, 0.38, 0.20, 0.00, 0.00,ti+st,0,0,107),  &
! updates for 121 - st only
ratb_t1(138,'O3        ','MACR      ','OH        ','MeCO3     ','          ',  &
'          ',1.68e-16,  0.00, 2100.00, 0.38, 0.20, 0.00, 0.00,st,0,0,121),     &
! B139 IUPAC2007*
ratb_t1(139,'O3        ','MACR      ','OH        ','MeCO3     ','          ',  &
'          ',3.50e-16,  0.00, 2100.00, 0.38, 0.20, 0.00, 0.00,ti+st,0,0,107),  &
! B140 JPL2011
ratb_t1(140,'OClO      ','NO        ','NO2       ','ClO       ','          ',  &
'          ',2.50e-12,  0.00,  600.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(140,'OClO      ','NO        ','NO2       ','ClO       ','          ',  &
'          ',3.40e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B141 IUPAC2007
ratb_t1(141,'OH        ','C2H6      ','H2O       ','EtOO      ','          ',  &
'          ',6.90e-12,  0.00, 1000.00, 0.0, 0.0, 0.0,0.0,ti+t+st+r+cs,0,0,107),&
! updates for 121 - st & cs
ratb_t1(141,'OH        ','C2H6      ','H2O       ','EtOO      ','          ',  &
'          ',7.66e-12,  0.00, 1020.00, 0.0, 0.0, 0.0,0.0,st+cs,0,0,121),       &
! B142 IUPAC2007 see also asad_bimol
ratb_t1(142,'OH        ','C3H8      ','i-PrOO    ','H2O       ','          ',  &
'          ',7.60e-12,  0.00,  585.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r,0,0,107),  &
! updates for 121 - st only
ratb_t1(142,'OH        ','C3H8      ','i-PrOO    ','H2O       ','          ',  &
'          ',6.76e-12,  0.00,  630.00, 0.0, 0.0, 0.0, 0.0,st,0,0,121),         &
! RAQ does not do n-PrOO production; this would get only half the loss rate
! for C3H8.
! B142 - Different rates in CRI. Note no n-PrOO
ratb_t1(142,'OH        ','C3H8      ','i-PrOO    ','          ','          ',  &
'          ',5.59e-12,  0.00,  585.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B143 IUPAC2007 see also asad_bimol
ratb_t1(143,'OH        ','C3H8      ','n-PrOO    ','H2O       ','          ',  &
'          ',7.60e-12,  0.00,  585.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! updates for 121 - st only
ratb_t1(143,'OH        ','C3H8      ','n-PrOO    ','H2O       ','          ',  &
'          ',2.43e-12,  0.00,  630.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
! B143 - Different rates in CRI, and different product (equiv. to n-PrOO)
ratb_t1(143,'OH        ','C3H8      ','RN10O2    ','          ','          ',  &
'          ',2.01e-12,  0.00,  585.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B144 IUPAC2009
ratb_t1(144,'OH        ','C5H8      ','ISO2      ','          ','          ',  &
'          ',2.70e-11,  0.00, -390.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,0,107),  &
! updates for 121 - st only
ratb_t1(144,'OH        ','C5H8      ','ISO2      ','          ','          ',  &
'          ',3.00e-11,  0.00, -390.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(144,'OH        ','C5H8      ','HOIPO2    ','H2O       ','          ',  &
'          ',2.54e-11,  0.00, -410.00, 0.00, 0.00, 0.00, 0.00,r,0,0,107),      &
! B144 - Different rates and prods in CRI
ratb_t1(144,'OH        ','C5H8      ','RU14O2    ','          ','          ',  &
'          ',2.54e-11,  0.00, -410.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(144,'OH        ','C5H8      ','RU14O2    ','          ','          ',  &
'          ',2.70e-11,  0.00, -390.00, 1.00, 0.015, 0.00, 0.00,cs,0,0,119),    &
! B145 JPL2011
ratb_t1(145,'OH        ','CH4       ','H2O       ','MeOO      ','          ',  &
'          ',2.45e-12,  0.00, 1775.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
ratb_t1(145,'OH        ','CH4       ','H2O       ','MeOO      ','          ',  &
'          ',1.85e-12,  0.00, 1690.00, 0.00, 0.00,0.00,0.00,ti+s+t+cs,0,0,107),&
! B146 IUPAC2005 see also asad_bimol
ratb_t1(146,'OH        ','CO        ','HO2       ','          ','          ',  &
'          ',1.44e-13,  0.00,    0.00, 0.0, 0.0,0.0,0.0,ti+t+st+r+cs,0,0,107), &
! updates for 121 - st & cs - NOTE CHANGE OF PRODUCTS
ratb_t1(146,'OH        ','CO        ','H         ','CO2       ','          ',  &
'          ',1.85e-13,  0.00,   65.00, 0.0, 0.0,0.0,0.0,st+cs,0,0,121),        &
ratb_t1(146,'OH        ','CO        ','H         ','CO2       ','          ',  &
'          ',1.30e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B147 JPL2011
ratb_t1(147,'OH        ','ClO       ','HCl       ','O2        ','          ',  &
'          ',6.00e-13,  0.00, -230.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107) &
]


! Start and end bounds for 6th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+61
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B148 JPL2011
ratb_t1(148,'OH        ','ClO       ','HO2       ','Cl        ','          ',  &
'          ',7.40e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B149 JPL2011
ratb_t1(149,'OH        ','ClONO2    ','HOCl      ','NO3       ','          ',  &
'          ',1.20e-12,  0.00,  330.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B150 IUPAC2007
ratb_t1(150,'OH        ','EtCHO     ','H2O       ','EtCO3     ','          ',  &
'          ',4.90e-12,  0.00, -405.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
ratb_t1(150,'OH        ','EtCHO     ','H2O       ','EtCO3     ','          ',  &
'          ',5.10e-12,  0.00, -405.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
! - Different rates in CRI
ratb_t1(150,'OH        ','EtCHO     ','H2O       ','EtCO3     ','          ',  &
'          ',1.96e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(150,'OH        ','EtCHO     ','H2O       ','EtCO3     ','          ',  &
'          ',4.9e-12,  0.00,  -405.0, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),      &
! B151 MCMv3.2
ratb_t1(151,'OH        ','EtOOH     ','H2O       ','EtOO      ','          ',  &
'          ',1.90e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! B152 MCMv3.2
ratb_t1(152,'OH        ','EtOOH     ','H2O       ','MeCHO     ','OH        ',  &
'          ',8.01e-12,  0.00,    0.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r,0,0,107),  &
! Different reaction rates in CRI, still using MCM3.1 rates
ratb_t1(152,'OH        ','EtOOH     ','H2O       ','MeCHO     ','OH        ',  &
'          ',1.36e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B153 JPL2011
! Products should be H2O + H not H2O + HO2 in ST (and CS) chemistry.
ratb_t1(153,'OH        ','H2        ','H2O       ','HO2       ','          ',  &
'          ',2.80e-12,  0.00, 1800.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
! updates for 121 - st change products
ratb_t1(153,'OH        ','H2        ','H2O       ','H         ','          ',  &
'          ',2.80e-12,  0.00, 1800.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
! Products should be H2O + H not H2O + HO2 in ST (and CS) chemistry.
ratb_t1(153,'OH        ','H2        ','H2O       ','HO2       ','          ',  &
'          ',7.70e-12,  0.00, 2100.00, 0.00, 0.00,0.00,0.0,ti+t+r+cs,0,0,107), &
! updates for 121 - cs change products (& keep cs rates)
ratb_t1(153,'OH        ','H2        ','H2O       ','H         ','          ',  &
'          ',7.70e-12,  0.00, 2100.00, 0.00, 0.00,0.00,0.0,cs,0,0,121),        &
ratb_t1(153,'OH        ','H2        ','H2O       ','H         ','          ',  &
'          ',7.70e-12,  0.00, 2100.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B154 IUPAC2001
ratb_t1(154,'OH        ','H2O2      ','H2O       ','HO2       ','          ',  &
'          ',2.90e-12,  0.00,  160.00, 0.0,0.0,0.0,0.0,ti+s+t+st+r+cs,0,0,107),&
ratb_t1(154,'OH        ','H2O2      ','H2O       ','          ','          ',  &
'          ',2.90e-12,  0.00,  160.00, 0.00, 0.00, 0.00, 0.00,ol,0,0,107),     &
! B155 IUPAC2007
ratb_t1(155,'OH        ','HACET     ','MGLY      ','HO2       ','          ',  &
'          ',1.60e-12,  0.00, -305.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
! updates for 121
ratb_t1(155,'OH        ','HACET     ','MGLY      ','HO2       ','          ',  &
'          ',2.00e-12,  0.00, -320.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(155,'OH        ','HACET     ','MGLY      ','HO2       ','          ',  &
'          ',3.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti,0,0,107),     &
! B156 JPL2011
ratb_t1(156,'OH        ','HBr       ','H2O       ','Br        ','          ',  &
'          ',5.50e-12,  0.00, -200.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(156,'OH        ','HBr       ','H2O       ','Br        ','          ',  &
'          ',1.10e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B157 IUPAC2007
ratb_t1(157,'OH        ','HCHO      ','H2O       ','HO2       ','CO        ',  &
'          ',5.40e-12,  0.00, -135.00, 0.0, 0.0,0.0,0.0,ti+t+st+r+cs,0,0,107), &
! updates for 121 - st & cs
ratb_t1(157,'OH        ','HCHO      ','H2O       ','HO2       ','CO        ',  &
'          ',5.50e-12,  0.00, -125.00, 0.0, 0.0,0.0,0.0,st+cs,0,0,121),        &
ratb_t1(157,'OH        ','HCHO      ','H2O       ','CO        ','HO2       ',  &
'          ',8.20e-12,  0.00,  -40.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B158 IUPAC2007
! more products needed here: OH + HCOOH -> H2O + CO2 + H
ratb_t1(158,'OH        ','HCOOH     ','HO2       ','          ','          ',  &
'          ',4.50e-13,  0.00,    0.00, 0.00, 0.00,0.00,0.00,ti+st+cs,0,0,107), &
! updates for 121 - st & cs only NOTE CHANGE OF PRODUCTS
ratb_t1(158,'OH        ','HCOOH     ','CO2       ','H2O       ','H         ',  &
'          ',4.00e-13,  0.00,    0.00, 0.00, 0.00,0.00,0.00,st+cs,0,0,121),    &
! B159 JPL2011
ratb_t1(159,'OH        ','HCl       ','H2O       ','Cl        ','          ',  &
'          ',1.80e-12,  0.00,  250.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(159,'OH        ','HCl       ','H2O       ','Cl        ','          ',  &
'          ',1.80e-12,  0.00,  240.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B160 JPL2011
ratb_t1(160,'OH        ','HO2       ','H2O       ','          ','          ',  &
'          ',4.80e-11,  0.00, -250.00, 0.0, 0.0,0.0,0.0,ti+t+st+r+cs,0,0,107), &
ratb_t1(160,'OH        ','HO2       ','H2O       ','O2        ','          ',  &
'          ',4.80e-11,  0.00, -250.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B161 IUPAC2007
ratb_t1(161,'OH        ','HO2NO2    ','H2O       ','NO2       ','O2        ',  &
'          ',3.20e-13,  0.00, -690.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
! updates for 121 - st only
ratb_t1(161,'OH        ','HO2NO2    ','H2O       ','NO2       ','O2        ',  &
'          ',4.50e-13,  0.00, -610.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(161,'OH        ','HO2NO2    ','H2O       ','NO2       ','O2        ',  &
'          ',1.30e-12,  0.00, -380.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
ratb_t1(161,'OH        ','HO2NO2    ','H2O       ','NO2       ','          ',  &
'          ',1.90e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,ti+t+cs,0,0,107),&
! CS2
! B161 !!! CRI v2.2 rate change *
ratb_t1(161,'OH        ','HO2NO2    ','H2O       ','NO2       ','O2        ',  &
'          ',3.20e-13,  0.00,-690.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),      &
! B162 JPL2011
ratb_t1(162,'OH        ','HOCl      ','ClO       ','H2O       ','          ',  &
'          ',3.00e-12,  0.00,  500.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,0,0,107),&
! B163 IUPAC2004
ratb_t1(163,'OH        ','HONO      ','H2O       ','NO2       ','          ',  &
'          ',2.50e-12,  0.00, -260.00, 0.00,0.00,0.00,0.00,ti+t+st+cs,0,0,107),&
! updates for 121 - st & cs only
ratb_t1(163,'OH        ','HONO      ','H2O       ','NO2       ','          ',  &
'          ',3.00e-12,  0.00, -250.00, 0.00,0.00,0.00,0.00,st+cs,0,0,121),     &
! B164 IUPAC2004 see also asad_bimol
ratb_t1(164,'OH        ','HONO2     ','H2O       ','NO3       ','          ',  &
'          ',2.40e-14,  0.00, -460.00, 0.00, 0.00, 0.00, 0.00,st+r+cs,0,0,107),&
! updates for 121 - st & cs only
ratb_t1(164,'OH        ','HONO2     ','H2O       ','NO3       ','          ',  &
'          ',3.70e-14,  0.00, -240.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,121),  &
ratb_t1(164,'OH        ','HONO2     ','H2O       ','NO3       ','          ',  &
'          ',1.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t+s,0,0,107), &
! B165 Poschl00
ratb_t1(165,'OH        ','ISON      ','HACET     ','NALD      ','          ',  &
'          ',1.30e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,0,107),  &
! B166 Poschl00
ratb_t1(166,'OH        ','ISOOH     ','MACR      ','OH        ','          ',  &
'          ',1.00e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,0,107),  &
! Very different rates & products
ratb_t1(166,'OH        ','ISOOH     ','MVK       ','HCHO      ','OH        ',  &
'          ',4.20e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,r,0,0,107),      &
! B167 IUPAC2007
ratb_t1(167,'OH        ','MACR      ','MACRO2    ','          ','          ',  &
'          ',1.30e-12,  0.00, -610.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,0,107),  &
! B168 IUPAC2007
ratb_t1(168,'OH        ','MACR      ','MACRO2    ','          ','          ',  &
'          ',4.00e-12,  0.00, -380.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,0,107),  &
! B169 MCMv3.2
ratb_t1(169,'OH        ','MACROOH   ','MACRO2    ','          ','          ',  &
'          ',3.77e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
ratb_t1(169,'OH        ','MACROOH   ','MACRO2    ','          ','          ',  &
'          ',3.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti,0,0,107),     &
! B170 IUPAC2008
ratb_t1(170,'OH        ','MGLY      ','MeCO3     ','CO        ','          ',  &
'          ',1.90e-12,  0.00, -575.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
ratb_t1(170,'OH        ','MGLY      ','MeCO3     ','CO        ','          ',  &
'          ',1.50e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti,0,0,107),     &
! Note the different rates
ratb_t1(170,'OH        ','MGLY      ','MeCO3     ','CO        ','          ',  &
'          ',1.72e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,r,0,0,107),      &
! B171 IUPAC2006.
ratb_t1(171,'OH        ','MPAN      ','HACET     ','NO2       ','          ',  &
'          ',2.90e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+st,0,0,107),  &
! Different products and rates in CRI
ratb_t1(171,'OH        ','MPAN      ','CARB7     ','CO        ','NO2       ',  &
'          ',3.60e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(171,'OH        ','MPAN      ','CARB7     ','CO        ','NO2       ',  &
'          ',6.38e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B172 IUPAC2007
ratb_t1(172,'OH        ','Me2CO     ','H2O       ','MeCOCH2OO ','          ',  &
'          ',1.70e-14,  0.00, -423.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
ratb_t1(172,'OH        ','Me2CO     ','H2O       ','MeCOCH2OO ','          ',  &
'          ',1.70e-14,  0.00, -420.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107),   &
! - Different products in CRI
ratb_t1(172,'OH        ','Me2CO     ','H2O       ','RN8O2     ','          ',  &
'          ',1.70e-14,  0.00, -423.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B173 IUPAC2007
ratb_t1(173,'OH        ','Me2CO     ','H2O       ','MeCOCH2OO ','          ',  &
'          ',8.80e-12,  0.00, 1320.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r,0,0,107),  &
! - Different products in CRI
ratb_t1(173,'OH        ','Me2CO     ','H2O       ','RN8O2     ','          ',  &
'          ',8.80e-12,  0.00, 1320.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B174 IUPAC2009
ratb_t1(174,'OH        ','MeCHO     ','H2O       ','MeCO3     ','          ',  &
'          ',4.70e-12,  0.00, -345.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
! updates for 121 - st only
ratb_t1(174,'OH        ','MeCHO     ','H2O       ','MeCO3     ','          ',  &
'          ',4.63e-12,  0.00, -350.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(174,'OH        ','MeCHO     ','H2O       ','MeCO3     ','          ',  &
'          ',4.40e-12,  0.00, -365.00, 0.00, 0.00, 0.00, 0.00,ti+t,0,0,107)    &
]

! Start and end bounds for 7th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+68
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B174 - Different rates in CRI
ratb_t1(174,'OH        ','MeCHO     ','H2O       ','MeCO3     ','          ',  &
'          ',5.55e-12,  0.00, -311.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(174,'OH        ','MeCHO     ','H2O       ','MeCO3     ','          ',  &
'          ',4.70e-12,  0.00, -345.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B175 MCMv3.2
ratb_t1(175,'OH        ','MeCO2H    ','MeOO      ','          ','          ',  &
'          ',8.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
! updates for 121 - st only
ratb_t1(175,'OH        ','MeCO2H    ','MeOO      ','          ','          ',  &
'          ',3.15e-14,  0.00, -920.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(175,'OH        ','MeCO2H    ','MeOO      ','          ','          ',  &
'          ',4.00e-13,  0.00, -200.00, 0.00, 0.00, 0.00, 0.00,ti,0,0,107),     &
! B176 MCMv3.2
ratb_t1(176,'OH        ','MeCO3H    ','MeCO3     ','          ','          ',  &
'          ',3.70e-12,  0.00,    0.00, 0.00,0.00,0.00,0.00,ti+st+cs,0,0,107),  &
! B177 MCMv3.2
ratb_t1(177,'OH        ','MeCOCH2OOH','H2O       ','MeCOCH2OO ','          ',  &
'          ',1.90e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! B178 MCMv3.2
ratb_t1(178,'OH        ','MeCOCH2OOH','OH        ','MGLY      ','          ',  &
'          ',8.39e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! B179 IUPAC2006
ratb_t1(179,'OH        ','MeOH      ','HO2       ','HCHO      ','          ',  &
'          ',2.85e-12,  0.00,  345.00, 0.00,0.00,0.00,0.00,ti+st+r+cs,0,0,107),&
! updates for 121 - st & cs
ratb_t1(179,'OH        ','MeOH      ','HO2       ','HCHO      ','          ',  &
'          ',2.90e-12,  0.00,  345.00, 0.00,0.00,0.00,0.00,st+cs,0,0,121),     &
! B180 IUPAC2006
ratb_t1(180,'OH        ','MeONO2    ','HCHO      ','NO2       ','H2O       ',  &
'          ',4.00e-13,  0.00,  845.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! updates for 121 - st only
ratb_t1(180,'OH        ','MeONO2    ','HCHO      ','NO2       ','H2O       ',  &
'          ',8.00e-13,  0.00, 1000.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
! - (Very!) Different reaction rates in CRI
ratb_t1(180,'OH        ','MeONO2    ','HCHO      ','NO2       ','H2O       ',  &
'          ',1.00e-14,  0.00,-1060.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(180,'OH        ','MeONO2    ','HCHO      ','NO2       ','H2O       ',  &
'          ',4.00e-13,  0.00, 845.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),      &
! B181 IUPAC2007
! branching not in S scheme
ratb_t1(181,'OH        ','MeOOH     ','H2O       ','HCHO      ','OH        ',  &
'          ',2.12e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,st+r,0,0,107),   &
! updates for 121 - st only
ratb_t1(181,'OH        ','MeOOH     ','H2O       ','HCHO      ','OH        ',  &
'          ',1.10e-12,  0.00, -200.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
! Note the different rates
ratb_t1(181,'OH        ','MeOOH     ','H2O       ','HCHO      ','OH        ',  &
'          ',1.02e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,ti+t+cs,0,0,107),&
! CS2
ratb_t1(181,'OH        ','MeOOH     ','H2O       ','HCHO      ','OH        ',  &
'          ',2.12e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B182 IUPAC2007
ratb_t1(182,'OH        ','MeOOH     ','H2O       ','MeOO      ','          ',  &
'          ',1.89e-12,  0.00, -190.00, 0.0, 0.0, 0.0,0.0,ti+t+st+r+cs,0,0,107),&
! updates for 121 - st only
ratb_t1(182,'OH        ','MeOOH     ','H2O       ','MeOO      ','          ',  &
'          ',2.70e-12,  0.00, -200.00, 0.0, 0.0, 0.0,0.0,st,0,0,121),          &
ratb_t1(182,'OH        ','MeOOH     ','H2O       ','MeOO      ','          ',  &
'          ',2.70e-12,  0.00, -200.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! CS2
ratb_t1(182,'OH        ','MeOOH     ','H2O       ','MeOO      ','          ',  &
'          ',3.18e-12,  0.00, -190.00, 0.0, 0.0, 0.0,0.0,cs,0,0,119),          &
! B183 IUPAC2009
ratb_t1(183,'OH        ','NALD      ','HCHO      ','CO        ','NO2       ',  &
'          ',4.70e-12,  0.00, -345.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
ratb_t1(183,'OH        ','NALD      ','HCHO      ','CO        ','NO2       ',  &
'          ',4.40e-12,  0.00, -365.00, 0.00, 0.00, 0.00, 0.00,ti,0,0,107),     &
! B184 JPL2011
ratb_t1(184,'OH        ','NO3       ','HO2       ','NO2       ','          ',  &
'          ',2.20e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st,0,0,107),   &
ratb_t1(184,'OH        ','NO3       ','HO2       ','NO2       ','          ',  &
! updates for 121 - s & st
'          ',2.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s+st,0,0,121),   &
ratb_t1(184,'OH        ','NO3       ','HO2       ','NO2       ','          ',  &
'          ',2.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t+cs,0,0,107),&
! B185 JPL2011
ratb_t1(185,'OH        ','O3        ','HO2       ','O2        ','          ',  &
'          ',1.70e-12,  0.00,  940.00, 0.0,0.0,0.0,0.0,ti+t+s+st+r+cs,0,0,107),&
! B186 JPL2011
ratb_t1(186,'OH        ','OClO      ','HOCl      ','O2        ','          ',  &
'          ',1.40e-12,  0.00, -600.00, 0.00, 0.00, 0.00, 0.00,st+cs,0,0,107),  &
ratb_t1(186,'OH        ','OClO      ','HOCl      ','O2        ','          ',  &
'          ',4.50e-13,  0.00, -800.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B187 IUPAC2001
! Note the different signs of temperature dependence
! Added to CRI-Strat in case important in stratosphere
ratb_t1(187,'OH        ','OH        ','H2O       ','O(3P)     ','          ',  &
'          ',6.31e-14,  2.60, -945.00, 0.00,0.00,0.00,0.00,ti+t+st+cs,0,0,107),&
! updates for 121 - st & cs
ratb_t1(187,'OH        ','OH        ','H2O       ','O(3P)     ','          ',  &
'          ',1.80e-12,  0.00,    0.00, 0.00,0.00,0.00,0.00,st+cs,0,0,121),     &
ratb_t1(187,'OH        ','OH        ','H2O       ','O(3P)     ','          ',  &
'          ',4.20e-12,  0.00,  240.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B188 MCMv3.2
ratb_t1(188,'OH        ','PAN       ','HCHO      ','NO2       ','H2O       ',  &
'          ',3.00e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
ratb_t1(188,'OH        ','PAN       ','HCHO      ','NO3       ','          ',  &
'          ',9.50e-13,  0.00,  650.00, 0.00, 0.00, 0.00, 0.00,r,0,0,107),      &
! - Different prods in CRI
ratb_t1(188,'OH        ','PAN       ','HCHO      ','CO        ','NO2       ',  &
'H2O       ',9.50e-13,  0.00,  650.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(188,'OH        ','PAN       ','HCHO      ','CO        ','NO2       ',  &
'H2O       ',3.0e-14,  0.00,  0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),        &
! B189 MCMv3.2
ratb_t1(189,'OH        ','PPAN      ','MeCHO     ','NO2       ','H2O       ',  &
'          ',1.27e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! Different prods in CRI
! http://mcm.leeds.ac.uk/MCMv3.2/browse.htt?species=PAN
! PPN + OH -> CH3CHO + CO + NO2 (+ H2O)
ratb_t1(189,'OH        ','PPAN      ','MeCHO     ','CO        ','NO2       ',  &
'H2O       ',1.27e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B190 MCMv3.2
ratb_t1(190,'OH        ','i-PrOOH   ','Me2CO     ','OH        ','          ',  &
'          ',1.66e-11,  0.00,    0.00, 0.00, 0.0, 0.0, 0.0,ti+t+st+r,0,0,107), &
! Different reaction rate, no branch in CRI
ratb_t1(190,'OH        ','i-PrOOH   ','Me2CO     ','OH        ','          ',  &
'          ',2.78e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B191 MCMv3.2
ratb_t1(191,'OH        ','i-PrOOH   ','i-PrOO    ','H2O       ','          ',  &
'          ',1.90e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! B192 MCMv3.2
ratb_t1(192,'OH        ','n-PrOOH   ','EtCHO     ','H2O       ','OH        ',  &
'          ',1.10e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! B193 MCMv3.2
ratb_t1(193,'OH        ','n-PrOOH   ','n-PrOO    ','H2O       ','          ',  &
'          ',1.90e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! B194 IUPAC2005
ratb_t1(194,'i-PrOO    ','NO        ','Me2CO     ','HO2       ','NO2       ',  &
'          ',2.70e-12,  0.00, -360.00, 0.0, 0.0, 0.0, 0.0,ti+t+st+r,0,0,107),  &
! - (Slightly) different reaction rates in CRI - rest forms i-PrONO2
ratb_t1(194,'i-PrOO    ','NO        ','Me2CO     ','HO2       ','NO2       ',  &
'          ',2.59e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B195 MCMv3.2
ratb_t1(195,'i-PrOO    ','NO3       ','Me2CO     ','HO2       ','NO2       ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
! updates for 121 - st only
ratb_t1(195,'i-PrOO    ','NO3       ','Me2CO     ','HO2       ','NO2       ',  &
'          ',2.20e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(195,'i-PrOO    ','NO3       ','Me2CO     ','HO2       ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ti+t+cs,0,0,107),&
! CS2
ratb_t1(195,'i-PrOO    ','NO3       ','Me2CO     ','HO2       ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B196 IUPAC2005
ratb_t1(196,'n-PrOO    ','NO        ','EtCHO     ','HO2       ','NO2       ',  &
'          ',2.90e-12,  0.00, -350.00, 0.00, 0.00, 0.00, 0.00,ti+t+st,0,0,107),&
! B197 MCM3.2
ratb_t1(197,'n-PrOO    ','NO3       ','EtCHO     ','HO2       ','NO2       ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,st,0,0,107),     &
! updates for 121
ratb_t1(197,'n-PrOO    ','NO3       ','EtCHO     ','HO2       ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121),     &
ratb_t1(197,'n-PrOO    ','NO3       ','EtCHO     ','HO2       ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,t+ti,0,0,107),   &
! B198
! Missing in TI scheme
ratb_t1(198,'CS2       ','O(3P)     ','COS       ','SO2       ','CO        ',  &
'          ',3.20e-11,  0.00,  650.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,a,0,107),&
! updates for 121 - s, st & cs
ratb_t1(198,'CS2       ','O(3P)     ','CO        ','SO2       ','SO2       ',  &
'          ',3.30e-11,  0.00,  650.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,a,0,121),&
! B199
ratb_t1(199,'CS2       ','OH        ','COS       ','SO2       ','          ',  &
'          ',1.25e-16,  0.00,-4550.00, 0.00, 0.00,0.00,0.00,s+st+r+cs,a,0,107),&
! updates for 121 - s, st, & cs only
ratb_t1(199,'CS2       ','OH        ','COS       ','SO2       ','          ',  &
'          ',2.00e-15,  0.00,    0.00, 0.00, 0.00,0.00,0.00,s+st+cs,a,0,121),  &
ratb_t1(199,'CS2       ','OH        ','SO2       ','COS       ','          ',  &
'          ',8.80e-16,  0.00,  -2300.00, 0.00, 0.00, 0.00, 0.00, ti,a,0,107),  &
! B200
ratb_t1(200,'DMS       ','OH        ','SO2       ','          ','          ',  &
'          ',1.20e-11,  0.00,  260.00, 0.00, 0.00, 0.00, 0.00, s+st+r,a,0,107),&
! updates for 121 - s & st only
ratb_t1(200,'DMS       ','OH        ','SO2       ','          ','          ',  &
'          ',1.20e-11,  0.00,  280.00, 0.00, 0.00, 0.00, 0.00, s+st,a,0,121),  &
! Fractional production of MeOO can't be 0
ratb_t1(200,'DMS       ','OH        ','SO2       ','DMSO      ','MeOO      ',  &
'          ',3.04e-12,  0.00, -350.00, 0.60, 0.40, 0.00, 0.00, ti,a,0,107),    &
ratb_t1(200,'DMS       ','OH        ','SO2       ','          ','          ',  &
'          ',9.60e-12,  0.00,  240.00, 0.00, 0.00, 0.00, 0.00, ol,a,0,107),    &
! B201
ratb_t1(201,'DMS       ','OH        ','SO2       ','MeOO      ','HCHO      ',  &
'          ',9.60e-12,  0.00,  240.00, 0.00, 0.00, 0.00, 0.00, ti,a,0,107),    &
! Difference in fractional products. Inconsistent?
ratb_t1(201,'DMS       ','OH        ','MSA       ','SO2       ','          ',  &
'          ',3.04e-12,  0.00, -350.00, 0.00, 0.00, 0.00, 0.00, s+r,a,0,107),   &
! updates for 121 - REMOVE for st as covered by T042
ratb_t1(201,'DMS       ','OH        ','MSA       ','SO2       ','          ',  &
'          ',3.04e-12,  0.00, -350.00, 0.00, 0.00, 0.00, 0.00, st,a,rm,107),   &
ratb_t1(201,'DMS       ','OH        ','SO2       ','DMSO      ','MeOO      ',  &
'          ',3.04e-12,  0.00, -350.00, 0.60, 0.40, 1.00, 0.00, st,a,rm,117),   &
ratb_t1(201,'DMS       ','OH        ','SO2       ','DMSO      ','          ',  &
'          ',3.04e-12,  0.00, -350.00, 0.60, 0.40, 0.00, 0.00, ol,a,0,107),    &
! B202
! Incomplete products
ratb_t1(202,'DMS       ','NO3       ','SO2       ','          ','          ',  &
'          ',1.90e-13,  0.00, -500.00, 0.0, 0.0, 0.0, 0.0, s+st+ol+r,a,0,107)  &
]


! Start and end bounds for 8th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+60
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B202
! updates for 121 - s & st only
ratb_t1(202,'DMS       ','NO3       ','SO2       ','          ','          ',  &
'          ',1.90e-13,  0.00, -530.00, 0.0, 0.0, 0.0, 0.0, s+st,a,0,121),      &
ratb_t1(202,'DMS       ','NO3       ','SO2       ','HONO2     ','MeOO      ',  &
'HCHO      ',1.90e-13,  0.00, -500.00, 0.00, 0.00, 0.00, 0.00, ti,a,0,107),    &
! B203
! Not in TIA scheme
! updates for 121 - REMOVE for s & st
ratb_t1(203,'DMS       ','O(3P)     ','SO2       ','          ','          ',  &
'          ',1.30e-11,  0.00, -410.00, 0.00, 0.00, 0.00, 0.00, s+st,a,rm,107), &
! B204
ratb_t1(204,'DMSO      ','OH        ','SO2       ','          ','          ',  &
'          ',5.80e-11,  0.00,    0.00, 0.60, 0.00, 0.00, 0.00, ol+r,a,0,107),  &
ratb_t1(204,'DMSO      ','OH        ','SO2       ','MSA       ','          ',  &
'          ',5.80e-11,  0.00,    0.00, 0.60, 0.40, 0.00, 0.00, ti,a,0,107),    &
ratb_t1(204,'DMSO      ','OH        ','SO2       ','MSA       ','          ',  &
'          ',5.80e-11,  0.00,    0.00, 0.60, 0.40, 0.00, 0.00, st,a,0,117),    &
! updates for 121
ratb_t1(204,'DMSO      ','OH        ','SO2       ','MSA       ','          ',  &
'          ',6.10e-12,  0.00, -800.00, 0.60, 0.40, 0.00, 0.00, st,a,0,121),    &
! B205
! Not in TIA scheme. Incomplete products
ratb_t1(205,'H2S       ','O(3P)     ','OH        ','SO2       ','          ',  &
'          ',9.20e-12,  0.00, 1800.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,a,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(205,'H2S       ','O(3P)     ','OH        ','SO2       ','OH        ',  &
'O(3P)     ',9.20e-12,  0.00, 1800.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,a,0,121),&
! B206
ratb_t1(206,'H2S       ','OH        ','SO2       ','H2O       ','          ',  &
'          ',6.00e-12,  0.00,   75.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,a,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(206,'H2S       ','OH        ','SO2       ','H2O       ','          ',  &
'          ',3.30e-12,  0.00, -100.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,a,0,121),&
ratb_t1(206,'H2S       ','OH        ','SO2       ','          ','          ',  &
'          ',6.00e-12,  0.00,   75.00, 0.00, 0.00, 0.00, 0.00, ti,a,0,107),    &
! B207
! Not in TIA scheme
ratb_t1(207,'COS       ','O(3P)     ','CO        ','SO2       ','          ',  &
'          ',2.10e-11,  0.00, 2200.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,a,0,107),&
! B208
ratb_t1(208,'COS       ','OH        ','CO2       ','SO2       ','          ',  &
'          ',1.10e-13,  0.00, 1200.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,a,0,107),&
! updates for 121 - s, st, & cs
ratb_t1(208,'COS       ','OH        ','CO2       ','SO2       ','          ',  &
'          ',7.20e-14,  0.00, 1070.00, 0.00, 0.00, 0.00, 0.00,s+st+cs,a,0,121),&
ratb_t1(208,'COS       ','OH        ','SO2       ','          ','          ',  &
'          ',1.10e-13,  0.00, 1200.00, 0.00, 0.00, 0.00, 0.00, ti,a,0,107),    &
! B209
! Not in TI scheme
ratb_t1(209,'SO2       ','O3        ','SO3       ','          ','          ',  &
'          ',3.00e-12,  0.00, 7000.00, 0.00, 0.00, 0.00, 0.00, s+st,a,0,107),  &
! B210
! Not in TI scheme
ratb_t1(210,'SO3       ','H2O       ','H2SO4     ','H2O       ','          ',  &
'          ',8.50e-41,  0.00,-6540.00, 0.00, 0.00, 0.00, 0.00, s+st,a,0,107),  &
! B210 - Different rate in CRI
ratb_t1(210,'SO3       ','H2O       ','H2SO4     ','H2O       ','          ',  &
'          ',3.90e-41,  0.00,-6830.60, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B211
ratb_t1(211,'Monoterp  ','OH        ','Sec_Org   ','          ','          ',  &
'          ',1.20e-11,  0.00, -444.0, 0.13, 0.0, 0.0, 0.0, ti+st+ol+r,a,0,107),&
! updates for 121 - st only
ratb_t1(211,'Monoterp  ','OH        ','Sec_Org   ','          ','          ',  &
'          ',1.34e-11,  0.00, -410.0, 0.13, 0.0, 0.0, 0.0, st,a,0,121),        &
! B212
ratb_t1(212,'Monoterp  ','O3        ','Sec_Org   ','          ','          ',  &
'          ',1.01e-15,  0.00,  732.00, 0.13, 0.0, 0.0, 0.0, ti+st+ol,a,0,107), &
! updates for 121 - st only
ratb_t1(212,'Monoterp  ','O3        ','Sec_Org   ','          ','          ',  &
'          ',8.22e-16,  0.00,  640.00, 0.13, 0.0, 0.0, 0.0, st,a,0,121),       &
! B213
ratb_t1(213,'Monoterp  ','NO3       ','Sec_Org   ','          ','          ',  &
'          ',1.19e-12,  0.00, -925.00, 0.13, 0.0, 0.0, 0.0, ti+st+ol,a,0,107), &
! updates for 121 - st only
ratb_t1(213,'Monoterp  ','NO3       ','Sec_Org   ','          ','          ',  &
'          ',1.20e-12,  0.00, -490.00, 0.13, 0.0, 0.0, 0.0, st,a,0,121),       &
! B214
! Make sure the below five reactions are in sync with the corresponding
! reactions listed above.
ratb_t1(214,'HO2S      ','O3S       ','HO2       ','O2        ','          ',  &
'          ',2.03e-16,  4.57, -693.00, 0.00, 0.00, 0.00, 0.00, t,0,0,107),     &
! B215
ratb_t1(215,'OHS       ','O3S       ','OH        ','O2        ','          ',  &
'          ',1.70e-12,  0.00,  940.00, 0.00, 0.00, 0.00, 0.00, t,0,0,107),     &
! B216
ratb_t1(216,'O(1D)S    ','H2O       ','H2O       ','          ','          ',  &
'          ',2.20e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, t,0,0,107),     &
! B217
ratb_t1(217,'O(1D)S    ','N2        ','O(3P)S    ','N2        ','          ',  &
'          ',2.10e-11,  0.00, -115.00, 0.00, 0.00, 0.00, 0.00, t,0,0,107),     &
! B218
ratb_t1(218,'O(1D)S    ','O2        ','O(3P)S    ','O2        ','          ',  &
'          ',3.20e-11,  0.00,  -67.00, 0.00, 0.00, 0.00, 0.00, t,0,0,107),     &
! B219
! These are reactions that are in R but not in the other mechs. Should they be?
ratb_t1(219,'NO3       ','NO3       ','NO2       ','NO2       ','O2        ',  &
'          ',8.50e-13,  0.00, 2450.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B220
ratb_t1(220,'OH        ','NH3       ','NH2       ','H2O       ','          ',  &
'          ',3.50e-12,  0.00,  925.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B221
ratb_t1(221,'NO3       ','HO2       ','HONO2     ','O2        ','          ',  &
'          ',4.20e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B222
ratb_t1(222,'EtOO      ','MeOO      ','MeCHO     ','HO2       ','HO2       ',  &
'HCHO      ',2.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B223
ratb_t1(223,'C4H10     ','OH        ','s-BuOO    ','H2O       ','          ',  &
'          ',7.90e-13,  2.00, -300.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B223 - Different rates and prods in CRI
ratb_t1(223,'OH        ','C4H10     ','RN13O2    ','H2O       ','          ',  &
'          ',9.80e-12,  0.00, -425.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B224
ratb_t1(224,'s-BuOO    ','NO        ','MEK       ','HO2       ','NO2       ',  &
'          ',2.54e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B225
ratb_t1(225,'s-BuOO    ','MeOO      ','MEK       ','HO2       ','HO2       ',  &
'HCHO      ',2.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B226
ratb_t1(226,'OH        ','MEK       ','MEKO2     ','H2O       ','          ',  &
'          ',1.30e-12,  0.00,   25.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! Different rate and products in CRI
ratb_t1(226,'MEK       ','OH        ','RN11O2    ','H2O       ','          ',  &
'          ',1.50e-12,  0.00,   90.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B227
! don't know how to map the products onto the RAC species
ratb_t1(227,'EtOO      ','EtOO      ','EtO       ','EtO       ','O2        ',  &
'          ',6.40e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B228
ratb_t1(228,'MeCO3     ','MeCO3     ','MeOO      ','MeOO      ','CO2       ',  &
'CO2       ',2.90e-12,  0.00, -500.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B229
ratb_t1(229,'MeOO      ','MeCOCH2OO ','MeCO3     ','HCHO      ','HCHO      ',  &
'HO2       ',3.80e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B230
ratb_t1(230,'i-PrOO    ','MeOO      ','Me2CO     ','HCHO      ','HO2       ',  &
'HO2       ',4.00e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B231
ratb_t1(231,'HO2       ','s-BuOO    ','s-BuOOH   ','O2        ','          ',  &
'          ',1.82e-13,  0.00,  -1300.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),   &
! B232
ratb_t1(232,'OH        ','s-BuOOH   ','MEK       ','OH        ','          ',  &
'          ',2.15e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B233
! how does MeCOCH2OOMe map onto the species?
ratb_t1(233,'NO        ','MeCOCHO2Me','MeCO3     ','MeCHO     ','NO2       ',  &
'          ',2.54e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B234
ratb_t1(234,'MeOO      ','MeCOCHO2Me','HCHO      ','HO2       ','MeCO3     ',  &
'MeCHO     ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B235
ratb_t1(235,'HOC2H4O2  ','NO        ','HCHO      ','HCHO      ','HO2       ',  &
'NO2       ',9.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B236
ratb_t1(236,'MeOO      ','HOC2H4O2  ','HCHO      ','HO2       ','O2        ',  &
'          ',2.00e-12,  0.00,    0.00, 3.00, 2.00, 1.00, 0.00, r,0,0,107),     &
! B237
! products: HCHO, 0.47 CH2O2, 0.31 CO, 0.22 CO2, 0.31 H2O, 0.13 H2 and 0.2 HO2.
! How does CH2O2 map onto the RAQ species? If it is dropped, and CO2 and H2O
! are dropped, this can be combined into one.
ratb_t1(237,'O3        ','C2H4      ','HCHO      ','CH2O2     ','CO        ',  &
'CO2       ',6.00e-15,  0.00, 2630.00, 2.00, 0.94, 0.62, 0.44, r,0,0,107),     &
! B237 - Different rates and prods in CRI
ratb_t1(237,'O3        ','C2H4      ','HCHO      ','CO        ','HO2       ',  &
'OH        ',1.19e-15,  0.00, 2580.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(237,'O3        ','C2H4      ','HCHO      ','CO        ','HO2       ',  &
'OH        ',1.18e-15,  0.00, 2580.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B238
ratb_t1(238,'O3        ','C2H4      ','H2O       ','H2        ','HO2       ',  &
'CO2       ',6.00e-15,  0.00, 2630.00, 0.62, 0.26, 0.40, 0.00, r,0,0,107),     &
! B238 - Different rates and prods in CRI
ratb_t1(238,'O3        ','C2H4      ','HCHO      ','HCOOH     ','          ',  &
'          ',7.95e-15,  0.00, 2580.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(238,'O3        ','C2H4      ','HCHO      ','HCOOH     ','          ',  &
'          ',7.91e-15,  0.00, 2580.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B239
ratb_t1(239,'O3        ','C3H6      ','HCHO      ','CH4       ','CO        ',  &
'CO2       ',1.375e-15, 0.00, 1878.00, 2.00, 0.60, 0.80, 1.20, r,0,0,107),     &
! B239 - Different rates and prods in CRI
ratb_t1(239,'O3        ','C3H6      ','HCHO      ','CO        ','MeOO      ',  &
'OH        ',1.98e-15,  0.00, 1878.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B240
ratb_t1(240,'O3        ','C3H6      ','OH        ','MeOH      ','HO2       ',  &
'MeOO      ',1.375e-15, 0.00, 1878.00, 0.56, 0.24, 0.60, 1.16, r,0,0,107),     &
! B240 - Different rates and prods in CRI
ratb_t1(240,'O3        ','C3H6      ','HCHO      ','MeCO2H    ','          ',  &
'          ',3.53e-15,  0.00, 1878.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B241
! products 0.420 CO2 and 0.580 H2O are dropped here
ratb_t1(241,'O3        ','C3H6      ','MeCHO     ','H2        ','CO        ',  &
'HO2       ',2.75e-15,  0.00, 1878.00, 1.00, 0.24, 1.16, 0.18, r,0,0,107)      &
]

! Start and end bounds for 9th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+52
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B242
ratb_t1(242,'HOC3H6O2  ','NO        ','HCHO      ','HO2       ','MeCHO     ',  &
'NO2       ',2.54e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B243
ratb_t1(243,'MeOO      ','HOC3H6O2  ','HCHO      ','HO2       ','MeCHO     ',  &
'O2        ',6.00e-13,  0.00,    0.00, 2.00, 2.00, 1.00, 1.00, r,0,0,107),     &
! B244
ratb_t1(244,'O3        ','MVK       ','MGLY      ','CO        ','CH2O2     ',  &
'HO2       ',3.78e-15,  0.00,  -1521.00, 2.00, 1.52, 0.48, 0.72, r,0,0,107),   &
! B245
ratb_t1(245,'O3        ','MVK       ','OH        ','          ','          ',  &
'          ',3.78e-15,  0.00,  -1521.00, 0.72, 0.00, 0.00, 0.00, r,0,0,107),   &
! B246
ratb_t1(246,'HOIPO2    ','MeOO      ','HO2       ','HO2       ','HCHO      ',  &
'MVK       ',5.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B247
ratb_t1(247,'HOMVKO2   ','MeOO      ','HO2       ','HO2       ','HCHO      ',  &
'MGLYOX    ',2.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B248
ratb_t1(248,'HOIPO2    ','HO2       ','ISOOH     ','O2        ','          ',  &
'          ',2.45e-13,  0.00,  -1250.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),   &
! B249
ratb_t1(249,'HOMVKO2   ','HO2       ','MVKOOH    ','O2        ','          ',  &
'          ',2.23e-13,  0.00,  -1250.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),   &
! B250
ratb_t1(250,'MVKOOH    ','OH        ','MGLY      ','HCHO      ','OH        ',  &
'          ',5.77e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B251
ratb_t1(251,'GLY       ','OH        ','HO2       ','CO        ','CO        ',  &
'          ',1.14e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B252
ratb_t1(252,'HOC5H8O2  ','NO        ','MVK       ','HO2       ','HCHO      ',  &
'NO2       ',2.08e-12,  0.00, -180.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B253
ratb_t1(253,'OH        ','MVK       ','HOMVKO2   ','H2O       ','          ',  &
'          ',4.13e-12,  0.00, -452.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B254
! don't know how to map MeCOCHO
ratb_t1(254,'HOMVKO2   ','NO        ','MGLY      ','HCHO      ','HO2       ',  &
'NO2       ',2.50e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B255
ratb_t1(255,'NO3       ','C2H6      ','EtOO      ','HONO2     ','          ',  &
'          ',5.70e-12,  0.00, 4426.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B256
ratb_t1(256,'NO3       ','C4H10     ','s-BuOO    ','HONO2     ','          ',  &
'          ',2.80e-12,  0.00, 3280.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B257
! 'CH2(NO3)CHO' is too long for the ratb_t1 constructor (11 characters
! instead of 10). We truncate `prod1` to be 'CH2(NO3)CH'.
ratb_t1(257,'NO3       ','C2H4      ','CH2(NO3)CH','HO2       ','          ',  &
'          ',3.30e-12,  2.00, 2880.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B257 - Different rates and prods in CRI
ratb_t1(257,'NO3       ','C2H4      ','NRN6O2    ','          ','          ',  &
'          ',2.10e-16,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(257,'NO3       ','C2H4      ','NRN6O2    ','          ','          ',  &
'          ',3.3e-12,  0.00,   2880.0, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B258
! 'CH2(NO3)CHO' is too long for the ratb_t1 constructor (11 characters
! instead of 10). We truncate `react1` to be 'CH2(NO3)CH'.
ratb_t1(258,'CH2(NO3)CH','OH        ','HCHO      ','NO2       ','CO2       ',  &
'          ',4.95e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B259
! 'CH3CH(NO3)CHO' is too long for the ratb_t1 constructor (13 characters
! instead of 10). We truncate `prod1` to be 'CH3CH(NO3)'.
ratb_t1(259,'NO3       ','C3H6      ','CH3CH(NO3)','HO2       ',               &
            '          ',                                                      &
'          ',4.59e-13,  0.00, 1156.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B259 - Different rates and prods in CRI
ratb_t1(259,'NO3       ','C3H6      ','NRN9O2    ','          ','          ',  &
'          ',9.40e-15,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B260
! 'CH3CH(NO3)CHO' is too long for the ratb_t1 constructor (13 characters
! instead of 10). We truncate `react1` to be 'CH3CH(NO3)'.
ratb_t1(260,'CH3CH(NO3)','OH        ','MeCHO     ','NO2       ',               &
            'CO2       ',                                                      &
'          ',5.25e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B261
! '(NO3)C4H6CHO' is too long for the ratb_t1 constructor (12 characters
! instead of 10). We truncate `prod1` to be '(NO3)C4H6C'.
ratb_t1(261,'(NO3)C4H6C','OH        ','MVK       ','NO2       ','          ',  &
'          ',4.16e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B262
ratb_t1(262,'OH        ','oXYLENE   ','HO2       ','MEMALD    ','MGLY      ',  &
'          ',1.36e-11,  0.00,    0.00, 1.00, 0.80, 0.80, 0.00, r,0,0,107),     &
! B263
ratb_t1(263,'OXYL1     ','NO2       ','ORGNIT    ','          ','          ',  &
'          ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B264
ratb_t1(264,'OH        ','MEMALD    ','MEMALD1   ','          ','          ',  &
'          ',5.60e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B265
! MGLY = MeCOCHO; GLY =  CHOCHO
ratb_t1(265,'MEMALD1   ','NO        ','HO2       ','MGLY      ','GLY       ',  &
'NO2       ',2.54e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B266
ratb_t1(266,'OH        ','TOLUENE   ','MEMALD    ','GLY       ','HO2       ',  &
'          ',1.18e-12,  0.00, -338.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B266 - Different rates and prods in CRI
ratb_t1(266,'OH        ','TOLUENE   ','AROH17    ','HO2       ','          ',  &
'          ',3.26e-13,  0.00, -338.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(266,'OH        ','TOLUENE   ','AROH17    ','HO2       ','          ',  &
'          ',3.24e-13,  0.00, -340.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B267
ratb_t1(267,'OH        ','TOLUENE   ','TOLP1     ','          ','          ',  &
'          ',3.60e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B267 - Different rates and prods in CRI
ratb_t1(267,'OH        ','TOLUENE   ','RA16O2    ','          ','          ',  &
'          ',1.48e-12,  0.00, -338.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(267,'OH        ','TOLUENE   ','RA16O2    ','          ','          ',  &
'          ',1.48e-12,  0.00, -340.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B268
ratb_t1(268,'OH        ','oXYLENE   ','OXYL1     ','          ','          ',  &
'          ',1.36e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B269
ratb_t1(269,'OH        ','ORGNIT    ','MEMALD    ','GLY       ','NO2       ',  &
'          ',2.70e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B270
ratb_t1(270,'NO3       ','ORGNIT    ','MEMALD    ','GLY       ','NO2       ',  &
'NO2       ',7.00e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B271
ratb_t1(271,'OXYL1     ','HO2       ','MGLY      ','MEMALD    ','          ',  &
'          ',2.50e-13,  0.00,  -1300.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),   &
! B272
ratb_t1(272,'MEMALD1   ','MeOO      ','HO2       ','HCHO      ','MGLY      ',  &
'GLY       ',1.00e-13,  0.00,    0.00, 2.00, 1.00, 1.00, 1.00, r,0,0,107),     &
! B273
ratb_t1(273,'HO2       ','TOLP1     ','MEMALD    ','GLY       ','OH        ',  &
'          ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B274
ratb_t1(274,'NO2       ','TOLP1     ','ORGNIT    ','          ','          ',  &
'          ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,0,0,107),     &
! B275
! not sure about the rate or products
ratb_t1(275,'SO2       ','OH        ','          ','          ','          ',  &
'          ',0.00e+00,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, r,a,0,107),     &
! B276
! The below reaction is missing in ST.
ratb_t1(276,'MACRO2    ','MeOO      ','MGLY      ','HACET     ','MeCO3     ',  &
'HCHO      ',1.00e-12,  0.00,    0.00, 1.00, 0.75, 0.25, 2.75, ti,0,0,107),    &
! B277
ratb_t1(277,'MACRO2    ','MeOO      ','HO2       ','CO        ','          ',  &
'          ',1.00e-12,  0.00,    0.00, 1.17, 0.25, 0.00, 0.00, ti,0,0,107),    &
! RO2+RO2 permutation reactions, only activated if RP = .TRUE. and
! using StratTrop chemical mechanism
! B278 IUPAC2002 see also asad_bimol
ratb_t1(278,'MeOO      ','RO2       ','HO2       ','HCHO      ','          ',  &
'          ',2.06e-13, 0.00, -365.00, 0.00, 0.00, 0.00, 0.00, st,rp,0,107),    &
! B279 IUPAC2002 see also asad_bimol
ratb_t1(279,'MeOO      ','RO2       ','MeOH      ','HCHO      ','          ',  &
'          ',2.06e-13, 0.00, -365.00, 0.50, 0.50, 0.00, 0.00,st,rp,0,107),     &
! B280 MCM 3.3.1
ratb_t1(280,'EtOO      ','RO2       ','HO2       ','MeCHO     ','          ',  &
'          ',9.74e-14, 0.00, -182.50, 0.00, 0.00, 0.00, 0.00,st,rp,0,107),     &
! B281 MCM 3.3.1; 2nd & 3rd pathways merged into one reaction
! with MeCHO as product for branch 3, as EtOH is not a species in ST
ratb_t1(281,'EtOO      ','RO2       ','MeCHO     ','          ','          ',  &
'          ',6.50e-14, 0.00, -182.50, 0.0,  0.00, 0.00, 0.00,st,rp,0,107),     &
! B282 MCM 3.3.1
ratb_t1(282,'n-PrOO    ','RO2       ','HO2       ','EtCHO     ','          ',  &
'          ',3.89e-13, 0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,rp,0,107),     &
! B283 MCM 3.3.1; 2nd & 3rd pathways merged into one reaction
! EtCHO used instead of n-PrOH for branch 3
ratb_t1(283,'n-PrOO    ','RO2       ','EtCHO     ','          ','          ',  &
'          ',2.59e-13, 0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,rp,0,107),     &
! B284 MCM 3.3.1
ratb_t1(284,'i-PrOO    ','RO2       ','HO2       ','Me2CO     ','          ',  &
'          ',4.87e-13, 0.00,  917.50, 0.00, 0.00, 0.00, 0.00,st,rp,0,107),     &
! B285 MCM 3.3.1; 2nd & 3rd pathways merged into one reaction
! Me2CO used instead of i-PrOH for branch 3
ratb_t1(285,'i-PrOO    ','RO2       ','Me2CO     ','          ','          ',  &
'          ',3.25e-13, 0.00, 917.50, 0.0,  0.00,  0.00, 0.00,st,rp,0,107),     &
! B286 MCM 3.3.1; Two pathways merged into one reaction
ratb_t1(286,'MeCO3     ','RO2       ','MeCO2H    ','MeOO      ','          ',  &
'          ',1.00e-11, 0.00,    0.00, 0.30, 0.70, 0.00, 0.00,st,rp,0,107),     &
! B287 MCM 3.3.1; Two pathways merged into one reaction
!  - NOTE EtCO2H IS *NOT* CURRENTLY A SPECIES CONSIDERED BY ASAD
ratb_t1(287,'EtCO3     ','RO2       ','EtCO2H    ','EtOO      ','          ',  &
'          ',1.00e-11,  0.00,    0.00, 0.30, 0.70, 0.00, 0.00,st,rp,0,107)     &
]


! Start and end bounds for 10th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+80
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B288 MCM 3.3.1
ratb_t1(288,'MeCOCH2OO ','RO2       ','HCHO      ','MeCO3     ','          ',  &
'          ',2.01e-12, 0.00,    0.00, 0.00, 0.00, 0.00, 0.00,st,rp,0,107),     &
! B289 MCM 3.3.1; Two pathways merged into one reaction
ratb_t1(289,'MeCOCH2OO ','RO2       ','HACET     ','MGLY      ','          ',  &
'          ',1.34e-12, 0.00,    0.00, 0.50, 0.50, 0.00, 0.00,st,rp,0,107),     &
! B290 Poschl00; reformulated to RO2+RO2 style
ratb_t1(290,'ISO2      ','RO2       ','MACR      ','HCHO      ','HO2       ',  &
'          ',1.67e-12, 0.00,    0.00, 1.00, 0.50, 0.50, 0.00,st,rp,0,107),     &
! updates for 121
ratb_t1(290,'ISO2      ','RO2       ','MACR      ','HCHO      ','HO2       ',  &
'          ',1.26e-12, 0.00,    0.00, 1.00, 0.50, 0.50, 0.00,st,rp,0,121),     &
! B291 Poschl00; reformulated to RO2+RO2 style
ratb_t1(291,'MACRO2    ','RO2       ','HACET     ','MGLY      ','HCHO      ',  &
'CO        ',8.37e-13, 0.00,    0.00, 1.00, 1.00, 0.50, 0.50,st,rp,0,107),     &
! updates for 121
ratb_t1(291,'MACRO2    ','RO2       ','HACET     ','MGLY      ','HCHO      ',  &
'CO        ',2.10e-13, 0.00,    0.00, 1.00, 1.00, 0.50, 0.50,st,rp,0,121),     &
! B292 Poschl00; reformulated to RO2+RO2 style
ratb_t1(292,'MACRO2    ','RO2       ','HO2       ','          ','          ',  &
'          ',8.37e-13,  0.00,    0.00, 1.00, 0.00, 0.00, 0.00,st,rp,0,107),    &
! updates for 121
ratb_t1(292,'MACRO2    ','RO2       ','HO2       ','          ','          ',  &
'          ',2.10e-13,  0.00,    0.00, 1.00, 0.00, 0.00, 0.00,st,rp,0,121),    &
! From here on out, new CRI reactions!
! B313
ratb_t1(313,'OH        ','TBUT2ENE  ','RN12O2    ','          ','          ',  &
'          ',1.01e-11,  0.00, -550.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B314
ratb_t1(314,'NO3       ','TBUT2ENE  ','NRN12O2   ','          ','          ',  &
'          ',3.90e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B315
ratb_t1(315,'O3        ','TBUT2ENE  ','MeCHO     ','CO        ','MeOO      ',  &
'OH        ',4.58e-15,  0.00, 1059.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B316
ratb_t1(316,'O3        ','TBUT2ENE  ','MeCHO     ','MeCO2H    ','          ',  &
'          ',2.06e-15,  0.00, 1059.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B317
! APINENE/BPINENE+Ox reactions produce a fraction of Sec-Org, equivalent to the
! Monoterp+Ox reactions. In the future, a more realistic depiction of SOA-
! formation should be implemented with different yields for different reaction
! parthways.
ratb_t1(317,'APINENE   ','OH        ','RTN28O2   ','          ','          ',  &
'          ',1.20e-11,  0.00, -444.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,107),     &
ratb_t1(317,'APINENE   ','OH        ','RTN28O2   ','Sec_Org   ','          ',  &
'          ',1.20e-11,  0.00, -444.00, 1.00, 0.13, 0.00, 0.00,cs,a,0,107),     &
! B318
ratb_t1(318,'APINENE   ','NO3       ','NRTN28O2  ','          ','          ',  &
'          ',1.19e-12,  0.00, -490.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,107),     &
! CS2
ratb_t1(318,'APINENE   ','NO3       ','NRTN28O2  ','          ','          ',  &
'          ',1.20e-12,  0.00, -490.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,119),     &
ratb_t1(318,'APINENE   ','NO3       ','NRTN28O2  ','Sec_Org   ','          ',  &
'          ',1.19e-12,  0.00, -490.00, 1.00, 0.13, 0.00, 0.00,cs,a,0,107),     &
! CS2
ratb_t1(318,'APINENE   ','NO3       ','NRTN28O2  ','Sec_Org   ','          ',  &
'          ',1.20e-12,  0.00, -490.00, 1.00, 0.13, 0.00, 0.00,cs,a,0,119),     &
! B319
ratb_t1(319,'APINENE   ','O3        ','OH        ','Me2CO     ','RN18AO2   ',  &
'          ',8.08e-16,  0.00,  732.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,107),     &
! CS2
ratb_t1(319,'APINENE   ','O3        ','OH        ','Me2CO     ','RN18AO2   ',  &
'          ',6.44e-16,  0.00,  640.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,119),     &
ratb_t1(319,'APINENE   ','O3        ','OH        ','Me2CO     ','RN18AO2   ',  &
'Sec_Org   ',8.08e-16,  0.00,  732.00, 1.00, 1.00, 1.00, 0.13,cs,a,0,107),     &
! CS2
ratb_t1(319,'APINENE   ','O3        ','OH        ','Me2CO     ','RN18AO2   ',  &
'Sec_Org   ',6.44e-16,  0.00,  640.00, 1.00, 1.00, 1.00, 0.13,cs,a,0,119),     &
! B320
ratb_t1(320,'APINENE   ','O3        ','TNCARB26  ','H2O2      ','          ',  &
'          ',7.57e-17,  0.00,  732.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,107),     &
! CS2
ratb_t1(320,'APINENE   ','O3        ','TNCARB26  ','H2O2      ','          ',  &
'          ',1.41e-16,  0.00,  640.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,119),     &
ratb_t1(320,'APINENE   ','O3        ','TNCARB26  ','H2O2      ','Sec_Org   ',  &
'          ',7.57e-17,  0.00,  732.00, 1.00, 1.00, 0.13, 0.00,cs,a,0,107),     &
! CS2
ratb_t1(320,'APINENE   ','O3        ','TNCARB26  ','H2O2      ','Sec_Org   ',  &
'          ',1.41e-16,  0.00,  640.00, 1.00, 1.00, 0.13, 0.00,cs,a,0,119),     &
! B321
ratb_t1(321,'APINENE   ','O3        ','RCOOH25   ','          ','          ',  &
'          ',1.26e-16,  0.00,  732.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,107),     &
! CS2
ratb_t1(321,'APINENE   ','O3        ','RCOOH25   ','          ','          ',  &
'          ',2.01e-17,  0.00,  640.00, 1.00, 0.0, 0.00, 0.00,cs,0,a,119),      &
ratb_t1(321,'APINENE   ','O3        ','RCOOH25   ','Sec_Org   ','          ',  &
'          ',1.26e-16,  0.00,  732.00, 1.00, 0.13, 0.00, 0.00,cs,a,0,107),     &
! CS2
ratb_t1(321,'APINENE   ','O3        ','RCOOH25   ','Sec_Org   ','          ',  &
'          ',2.01e-17,  0.00,  640.00, 1.00, 0.13, 0.00, 0.00,cs,a,0,119),     &
! B322
ratb_t1(322,'BPINENE   ','OH        ','RTX28O2   ','          ','          ',  &
'          ',2.38e-11,  0.00, -357.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,107),     &
ratb_t1(322,'BPINENE   ','OH        ','RTX28O2   ','Sec_Org   ','          ',  &
'          ',2.38e-11,  0.00, -357.00, 1.00, 0.13, 0.00, 0.00,cs,a,0,107),     &
! B323
ratb_t1(323,'BPINENE   ','NO3       ','NRTX28O2  ','          ','          ',  &
'          ',2.51e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,107),     &
ratb_t1(323,'BPINENE   ','NO3       ','NRTX28O2  ','Sec_Org   ','          ',  &
'          ',2.51e-12,  0.00,    0.00, 1.00, 0.13, 0.00, 0.00,cs,a,0,107),     &
! B324
ratb_t1(324,'BPINENE   ','O3        ','RTX24O2   ','OH        ','          ',  &
'          ',5.25e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,107),     &
! CS2
ratb_t1(324,'BPINENE   ','O3        ','RTX24O2   ','OH        ','          ',  &
'          ',4.73e-16,  0.00, 1270.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,119),     &
ratb_t1(324,'BPINENE   ','O3        ','RTX24O2   ','OH        ','Sec_Org   ',  &
'          ',5.25e-18,  0.00,    0.00, 1.00, 1.00, 0.13, 0.00,cs,a,0,107),     &
! CS2
ratb_t1(324,'BPINENE   ','O3        ','RTX24O2   ','OH        ','Sec_Org   ',  &
'          ',4.73e-16,  0.00, 1270.00, 1.00, 1.00, 0.13, 0.00,cs,a,0,119),     &
! B325
ratb_t1(325,'BPINENE   ','O3        ','HCHO      ','TXCARB24  ','H2O2      ',  &
'          ',3.00e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,107),     &
! CS2
ratb_t1(325,'BPINENE   ','O3        ','HCHO      ','TXCARB24  ','H2O2      ',  &
'          ',2.70e-16,  0.00,  1270.0, 0.00, 0.00, 0.00, 0.00,cs,0,a,119),     &
ratb_t1(325,'BPINENE   ','O3        ','HCHO      ','TXCARB24  ','H2O2      ',  &
'Sec_Org   ',3.00e-18,  0.00,    0.00, 1.00, 1.00, 1.00, 0.13,cs,a,0,107),     &
! CS2
ratb_t1(325,'BPINENE   ','O3        ','HCHO      ','TXCARB24  ','H2O2      ',  &
'Sec_Org   ',2.70e-16,  0.00,  1270.0, 1.00, 1.00, 1.00, 0.13,cs,a,0,119),     &
! B326
ratb_t1(326,'BPINENE   ','O3        ','HCHO      ','TXCARB22  ','          ',  &
'          ',3.75e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,107),     &
! CS2
ratb_t1(326,'BPINENE   ','O3        ','HCHO      ','TXCARB22  ','          ',  &
'          ',3.38e-16,  0.00,  1270.0, 0.00, 0.00, 0.00, 0.00,cs,0,a,119),     &
ratb_t1(326,'BPINENE   ','O3        ','HCHO      ','TXCARB22  ','Sec_Org   ',  &
'          ',3.75e-18,  0.00,    0.00, 1.00, 1.00, 0.13, 0.00,cs,a,0,107),     &
! CS2
ratb_t1(326,'BPINENE   ','O3        ','HCHO      ','TXCARB22  ','Sec_Org   ',  &
'          ',3.38e-16,  0.00, 1270.0, 1.00, 1.00, 0.13, 0.00,cs,a,0,119),      &
! B327
ratb_t1(327,'BPINENE   ','O3        ','TXCARB24  ','CO        ','          ',  &
'          ',3.00e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,107),     &
! CS2
ratb_t1(327,'BPINENE   ','O3        ','TXCARB24  ','CO        ','          ',  &
'          ',2.70e-16,  0.00 ,1270.00, 0.00, 0.00, 0.00, 0.00,cs,0,a,119),     &
ratb_t1(327,'BPINENE   ','O3        ','TXCARB24  ','CO        ','Sec_Org   ',  &
'          ',3.00e-18,  0.00,    0.00, 1.00, 1.00, 0.13, 0.00,cs,a,0,107),     &
! CS2
ratb_t1(327,'BPINENE   ','O3        ','TXCARB24  ','CO        ','Sec_Org   ',  &
'          ',2.70e-16,  0.00, 1270.00, 1.00, 1.00, 0.13, 0.00,cs,a,0,119),     &
! B328
ratb_t1(328,'BENZENE   ','OH        ','RA13O2    ','          ','          ',  &
'          ',1.10e-12,  0.00,  193.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(328,'BENZENE   ','OH        ','RA13O2    ','          ','          ',  &
'          ',1.08e-12,  0.00,  190.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B329
ratb_t1(329,'BENZENE   ','OH        ','AROH14    ','HO2       ','          ',  &
'          ',1.23e-12,  0.00,  193.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(329,'BENZENE   ','OH        ','AROH14    ','HO2       ','          ',  &
'          ',1.22e-12,  0.00,  190.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B330 - OH+oXYLENE reaction added here, two branches rather than one.
!        Could be merged into a single reaction...
ratb_t1(330,'oXYLENE   ','OH        ','RA19AO2   ','          ','          ',  &
'          ',9.52e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B331
ratb_t1(331,'oXYLENE   ','OH        ','RA19CO2   ','          ','          ',  &
'          ',4.08e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B332
ratb_t1(332,'OH        ','EtOH      ','MeCHO     ','HO2       ','          ',  &
'          ',2.85e-12,  0.00,  -20.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B333
ratb_t1(333,'OH        ','EtOH      ','HOCH2CH2O2','          ','          ',  &
'          ',1.50e-13,  0.00,  -20.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B334
ratb_t1(334,'n-PrOH    ','OH        ','EtCHO     ','HO2       ','          ',  &
'          ',2.71e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(334,'n-PrOH    ','OH        ','EtCHO     ','HO2       ','          ',  &
'          ',2.27e-12,  0.00,  -70.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B335
ratb_t1(335,'n-PrOH    ','OH        ','RN9O2     ','          ','          ',  &
'          ',2.82e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(335,'n-PrOH    ','OH        ','RN9O2     ','          ','          ',  &
'          ',2.33e-12,  0.00,  -70.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B336
ratb_t1(336,'OH        ','i-PrOH    ','Me2CO     ','HO2       ','          ',  &
'          ',2.24e-12,  0.00, -200.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B337
ratb_t1(337,'OH        ','i-PrOH    ','RN9O2     ','          ','          ',  &
'          ',3.61e-13,  0.00, -200.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B338
ratb_t1(338,'RN10O2    ','NO        ','EtCHO     ','HO2       ','NO2       ',  &
'          ',2.74e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B339
ratb_t1(339,'RN13O2    ','NO        ','MeCHO     ','EtOO      ','NO2       ',  &
'          ',8.76e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(339,'RN13O2    ','NO        ','MeCHO     ','EtOO      ','NO2       ',  &
'          ',9.854e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B340
ratb_t1(340,'RN13O2    ','NO        ','CARB11A   ','HO2       ','NO2       ',  &
'          ',1.32e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(340,'RN13O2    ','NO        ','CARB11A   ','HO2       ','NO2       ',  &
'          ',1.491e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B341
ratb_t1(341,'RN16O2    ','NO        ','RN15AO2   ','NO2       ','          ',  &
'          ',2.10e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(341,'RN16O2    ','NO        ','RN15AO2   ','NO2       ','          ',  &
'          ',2.368e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B342
ratb_t1(342,'RN19O2    ','NO        ','RN18AO2   ','NO2       ','          ',  &
'          ',1.89e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(342,'RN19O2    ','NO        ','RN18AO2   ','NO2       ','          ',  &
'          ',2.395e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B343
ratb_t1(343,'RN13AO2   ','NO        ','RN12O2    ','NO2       ','          ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(343,'RN13AO2   ','NO        ','RN12O2    ','NO2       ','          ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B344
ratb_t1(344,'RN16AO2   ','NO        ','RN15O2    ','NO2       ','          ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(344,'RN16AO2   ','NO        ','RN15O2    ','NO2       ','          ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B345
ratb_t1(345,'RA13O2    ','NO        ','CARB3     ','UDCARB8   ','HO2       ',  &
'NO2       ',2.20e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(345,'RA13O2    ','NO        ','CARB3     ','UDCARB8   ','HO2       ',  &
'NO2       ',2.48e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B346
ratb_t1(346,'RA16O2    ','NO        ','CARB3     ','UDCARB11  ','HO2       ',  &
'NO2       ',1.49e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(346,'RA16O2    ','NO        ','CARB3     ','UDCARB11  ','HO2       ',  &
'NO2       ',1.68e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119)      &
]



! Start and end bounds for 11th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+96
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B347
ratb_t1(347,'RA16O2    ','NO        ','CARB6     ','UDCARB8   ','HO2       ',  &
'NO2       ',6.40e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(347,'RA16O2    ','NO        ','CARB6     ','UDCARB8   ','HO2       ',  &
'NO2       ',7.20e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B348
ratb_t1(348,'RA19AO2   ','NO        ','CARB3     ','UDCARB14  ','HO2       ',  &
'NO2       ',2.07e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(348,'RA19AO2   ','NO        ','CARB3     ','UDCARB14  ','HO2       ',  &
'NO2       ',2.327e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B349
ratb_t1(349,'RA19CO2   ','NO        ','CARB9     ','UDCARB8   ','HO2       ',  &
'NO2       ',2.07e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(349,'RA19CO2   ','NO        ','CARB9     ','UDCARB8   ','HO2       ',  &
'NO2       ',2.327e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B350
ratb_t1(350,'HOCH2CH2O2','NO        ','HCHO      ','HCHO      ','HO2       ',  &
'NO2       ',1.85e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(350,'HOCH2CH2O2','NO        ','HCHO      ','HCHO      ','HO2       ',  &
'NO2       ',2.084e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B351
ratb_t1(351,'HOCH2CH2O2','NO        ','HOCH2CHO  ','HO2       ','NO2       ',  &
'          ',5.35e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(351,'HOCH2CH2O2','NO        ','HOCH2CHO  ','HO2       ','NO2       ',  &
'          ',6.018e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B352
ratb_t1(352,'RN9O2     ','NO        ','MeCHO     ','HCHO      ','HO2       ',  &
'NO2       ',2.35e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(352,'RN9O2     ','NO        ','MeCHO     ','HCHO      ','HO2       ',  &
'NO2       ',2.643e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B353
ratb_t1(353,'RN12O2    ','NO        ','MeCHO     ','MeCHO     ','HO2       ',  &
'NO2       ',2.30e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(353,'RN12O2    ','NO        ','MeCHO     ','MeCHO     ','HO2       ',  &
'NO2       ',2.589e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B354
ratb_t1(354,'RN15O2    ','NO        ','EtCHO     ','MeCHO     ','HO2       ',  &
'NO2       ',2.25e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(354,'RN15O2    ','NO        ','EtCHO     ','MeCHO     ','HO2       ',  &
'NO2       ',2.527e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B355
ratb_t1(355,'RN18O2    ','NO        ','EtCHO     ','EtCHO     ','HO2       ',  &
'NO2       ',2.17e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(355,'RN18O2    ','NO        ','EtCHO     ','EtCHO     ','HO2       ',  &
'NO2       ',2.438e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B356
ratb_t1(356,'RN15AO2   ','NO        ','CARB13    ','HO2       ','NO2       ',  &
'          ',2.34e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(356,'RN15AO2   ','NO        ','CARB13    ','HO2       ','NO2       ',  &
'          ',2.633e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B357
ratb_t1(357,'RN18AO2   ','NO        ','CARB16    ','HO2       ','NO2       ',  &
'          ',2.27e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(357,'RN18AO2   ','NO        ','CARB16    ','HO2       ','NO2       ',  &
'          ',2.554e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B358
ratb_t1(358,'HOCH2CO3  ','NO        ','HO2       ','HCHO      ','NO2       ',  &
'          ',8.10e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(358,'HOCH2CO3  ','NO        ','HO2       ','HCHO      ','NO2       ',  &
'          ',7.50e-12,  0.00, -290.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B359
ratb_t1(359,'RN8O2     ','NO        ','MeCO3     ','HCHO      ','NO2       ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(359,'RN8O2     ','NO        ','MeCO3     ','HCHO      ','NO2       ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B360
ratb_t1(360,'RN11O2    ','NO        ','MeCO3     ','MeCHO     ','NO2       ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(360,'RN11O2    ','NO        ','MeCO3     ','MeCHO     ','NO2       ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B361
ratb_t1(361,'RN14O2    ','NO        ','EtCO3     ','MeCHO     ','NO2       ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(361,'RN14O2    ','NO        ','EtCO3     ','MeCHO     ','NO2       ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B362
ratb_t1(362,'RN17O2    ','NO        ','RN16AO2   ','NO2       ','          ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(362,'RN17O2    ','NO        ','RN16AO2   ','NO2       ','          ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B363
ratb_t1(363,'RU14O2    ','NO        ','UCARB12   ','HO2       ','NO2       ',  &
'          ',5.44e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(363,'RU14O2    ','NO        ','UCARB12   ','HO2       ','NO2       ',  &
'          ',7.776e-14,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B364
ratb_t1(364,'RU14O2    ','NO        ','UCARB10   ','HCHO      ','HO2       ',  &
'NO2       ',1.62e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(364,'RU14O2    ','NO        ','UCARB10   ','HCHO      ','HO2       ',  &
'NO2       ',2.352e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B365
ratb_t1(365,'RU12O2    ','NO        ','MeCO3     ','HOCH2CHO  ','NO2       ',  &
'          ',1.68e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(365,'RU12O2    ','NO        ','CARB6     ','HOCH2CHO  ','NO2       ',  &
'HO2       ',1.099e-12, 0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B366
ratb_t1(366,'RU12O2    ','NO        ','CARB7     ','CO        ','HO2       ',  &
'NO2       ',7.20e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(366,'RU12O2    ','NO        ','CARB7     ','CARB3     ','NO2       ',  &
'HO2       ',1.048e-12,  0.00, -360.00, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),        &
! B367
ratb_t1(367,'RU10O2    ','NO        ','MeCO3     ','HOCH2CHO  ','NO2       ',  &
'          ',1.20e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(367,'RU10O2    ','NO        ','MeCO3     ','HOCH2CHO  ','NO2       ',  &
'          ',1.809e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B368
ratb_t1(368,'RU10O2    ','NO        ','CARB6     ','HCHO      ','HO2       ',  &
'NO2       ',7.20e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(368,'RU10O2    ','NO        ','CARB6     ','HCHO      ','HO2       ',  &
'NO2       ',7.965e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B369
ratb_t1(369,'RU10O2    ','NO        ','CARB7     ','HCHO      ','HO2       ',  &
'NO2       ',4.80e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(369,'RU10O2    ','NO        ','RU10NO3   ','          ','          ',  &
'          ',9.45e-14,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B370
ratb_t1(370,'NRN6O2    ','NO        ','HCHO      ','HCHO      ','NO2       ',  &
'NO2       ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(370,'NRN6O2    ','NO        ','HCHO      ','HCHO      ','NO2       ',  &
'NO2       ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B371
ratb_t1(371,'NRN9O2    ','NO        ','MeCHO     ','HCHO      ','NO2       ',  &
'NO2       ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(371,'NRN9O2    ','NO        ','MeCHO     ','HCHO      ','NO2       ',  &
'NO2       ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B372
ratb_t1(372,'NRN12O2   ','NO        ','MeCHO     ','MeCHO     ','NO2       ',  &
'NO2       ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(372,'NRN12O2   ','NO        ','MeCHO     ','MeCHO     ','NO2       ',  &
'NO2       ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B373
ratb_t1(373,'NRU14O2   ','NO        ','NUCARB12  ','HO2       ','NO2       ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(373,'NRU14O2   ','NO        ','NUCARB12  ','HO2       ','NO2       ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B374
ratb_t1(374,'NRU12O2   ','NO        ','NOA       ','CO        ','HO2       ',  &
'NO2       ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(374,'NRU12O2   ','NO        ','NOA       ','CO        ','HO2       ',  &
'NO2       ',1.35e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B375
ratb_t1(375,'RTN28O2   ','NO        ','TNCARB26  ','HO2       ','NO2       ',  &
'          ',1.68e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(375,'RTN28O2   ','NO        ','TNCARB26  ','HO2       ','NO2       ',  &
'          ',1.89e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B376
ratb_t1(376,'RTN28O2   ','NO        ','Me2CO     ','RN19O2    ','NO2       ',  &
'          ',1.56e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(376,'RTN28O2   ','NO        ','Me2CO     ','RN19O2    ','NO2       ',  &
'          ',1.76e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B377
ratb_t1(377,'NRTN28O2  ','NO        ','TNCARB26  ','NO2       ','NO2       ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(377,'NRTN28O2  ','NO        ','TNCARB26  ','NO2       ','NO2       ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B378
ratb_t1(378,'RTN26O2   ','NO        ','RTN25O2   ','NO2       ','          ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(378,'RTN26O2   ','NO        ','RTN25O2   ','NO2       ','          ',  &
'          ',7.50e-12,  0.00, -290.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B379
ratb_t1(379,'RTN25O2   ','NO        ','RTN24O2   ','NO2       ','          ',  &
'          ',2.02e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(379,'RTN25O2   ','NO        ','RTN24O2   ','NO2       ','          ',  &
'          ',2.27e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B380
ratb_t1(380,'RTN24O2   ','NO        ','RTN23O2   ','NO2       ','          ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(380,'RTN24O2   ','NO        ','RTN23O2   ','NO2       ','          ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B381
ratb_t1(381,'RTN23O2   ','NO        ','Me2CO     ','RTN14O2   ','NO2       ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(381,'RTN23O2   ','NO        ','Me2CO     ','RTN14O2   ','NO2       ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B382
ratb_t1(382,'RTN14O2   ','NO        ','HCHO      ','TNCARB10  ','HO2       ',  &
'NO2       ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(382,'RTN14O2   ','NO        ','HCHO      ','TNCARB10  ','HO2       ',  &
'NO2       ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B383
ratb_t1(383,'RTN10O2   ','NO        ','RN8O2     ','CO        ','NO2       ',  &
'          ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(383,'RTN10O2   ','NO        ','RN8O2     ','CO        ','NO2       ',  &
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B384
ratb_t1(384,'RTX28O2   ','NO        ','TXCARB24  ','HCHO      ','HO2       ',  &
'NO2       ',1.68e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(384,'RTX28O2   ','NO        ','TXCARB24  ','HCHO      ','HO2       ',  &
'NO2       ',1.89e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B385
ratb_t1(385,'RTX28O2   ','NO        ','Me2CO     ','RN19O2    ','NO2       ',  &
'          ',1.56e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B386
ratb_t1(386,'NRTX28O2  ','NO        ','TXCARB24  ','HCHO      ','NO2       ',  &
'NO2       ',2.40e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(386,'NRTX28O2  ','NO        ','TXCARB24  ','HCHO      ','NO2       ',  &
'NO2       ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B387
ratb_t1(387,'RTX24O2   ','NO        ','TXCARB22  ','HO2       ','NO2       ',  &
'          ',1.21e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(387,'RTX24O2   ','NO        ','TXCARB22  ','HO2       ','NO2       ',  &
'          ',1.37e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B388
ratb_t1(388,'RTX24O2   ','NO        ','Me2CO     ','RN13AO2   ','HCHO      ',  &
'NO2       ',8.09e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(388,'RTX24O2   ','NO        ','Me2CO     ','RN13AO2   ','HCHO      ',  &
'NO2       ',9.10e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B389
ratb_t1(389,'RTX22O2   ','NO        ','Me2CO     ','RN13O2    ','NO2       ',  &
'          ',1.68e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(389,'RTX22O2   ','NO        ','Me2CO     ','RN13O2    ','NO2       ',  &
'          ',1.89e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B390 - Organinitrate products not in other schemes
ratb_t1(390,'EtOO      ','NO        ','EtONO2    ','          ','          ',  &
'          ',2.34e-14,  0.00, -365.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(390,'EtOO      ','NO        ','EtONO2    ','          ','          ',  &
'          ',2.30e-14,  0.00, -380.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B391
ratb_t1(391,'RN10O2    ','NO        ','RN10NO3   ','          ','          ',  &
'          ',5.60e-14,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B392
ratb_t1(392,'i-PrOO    ','NO        ','i-PrONO2  ','          ','          ',  &
'          ',1.13e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B393
ratb_t1(393,'RN13O2    ','NO        ','RN13NO3   ','          ','          ',  &
'          ',1.99e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(393,'RN13O2    ','NO        ','RN13NO3   ','          ','          ',  &
'          ',2.241e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B394
ratb_t1(394,'RN16O2    ','NO        ','RN16NO3   ','          ','          ',  &
'          ',2.95e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(394,'RN16O2    ','NO        ','RN16NO3   ','          ','          ',  &
'          ',3.321e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B395
ratb_t1(395,'RN19O2    ','NO        ','RN19NO3   ','          ','          ',  &
'          ',5.09e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(395,'RN19O2    ','NO        ','RN19NO3   ','          ','          ',  &
'          ',5.724e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B396
ratb_t1(396,'HOCH2CH2O2','NO        ','HOC2H4NO3 ','          ','          ',  &
'          ',1.20e-14,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(396,'HOCH2CH2O2','NO        ','HOC2H4NO3 ','          ','          ',  &
'          ',1.35e-14,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119)      &
]


! Start and end bounds for 12th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+99
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B397
ratb_t1(397,'RN9O2     ','NO        ','RN9NO3    ','          ','          ',  &
'          ',5.04e-14,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(397,'RN9O2     ','NO        ','RN9NO3    ','          ','          ',  &
'          ',5.67e-14,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B398
ratb_t1(398,'RN12O2    ','NO        ','RN12NO3   ','          ','          ',  &
'          ',9.84e-14,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(398,'RN12O2    ','NO        ','RN12NO3   ','          ','          ',  &
'          ',1.107e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B399
ratb_t1(399,'RN15O2    ','NO        ','RN15NO3   ','          ','          ',  &
'          ',1.54e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(399,'RN15O2    ','NO        ','RN15NO3   ','          ','          ',  &
'          ',1.728e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B400
ratb_t1(400,'RN18O2    ','NO        ','RN18NO3   ','          ','          ',  &
'          ',2.33e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(400,'RN18O2    ','NO        ','RN18NO3   ','          ','          ',  &
'          ',2.619e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B401
ratb_t1(401,'RN15AO2   ','NO        ','RN15NO3   ','          ','          ',  &
'          ',6.00e-14,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(401,'RN15AO2   ','NO        ','RN15NO3   ','          ','          ',  &
'          ',6.75e-14,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B402
ratb_t1(402,'RN18AO2   ','NO        ','RN18NO3   ','          ','          ',  &
'          ',1.30e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(402,'RN18AO2   ','NO        ','RN18NO3   ','          ','          ',  &
'          ',1.458e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B403
ratb_t1(403,'RU14O2    ','NO        ','RU14NO3   ','          ','          ',  &
'          ',2.40e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(403,'RU14O2    ','NO        ','RU14NO3   ','          ','          ',  &
'          ',2.70e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B404
ratb_t1(404,'RA13O2    ','NO        ','RA13NO3   ','          ','          ',  &
'          ',1.97e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(404,'RA13O2    ','NO        ','RA13NO3   ','          ','          ',  &
'          ',2.214e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B405
ratb_t1(405,'RA16O2    ','NO        ','RA16NO3   ','          ','          ',  &
'          ',2.66e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(405,'RA16O2    ','NO        ','RA16NO3   ','          ','          ',  &
'          ',2.997e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B406
ratb_t1(406,'RA19AO2   ','NO        ','RA19NO3   ','          ','          ',  &
'          ',3.31e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2 20th
ratb_t1(406,'RA19AO2   ','NO        ','RA19NO3   ','          ','          ',  &
'          ',3.726e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B407
ratb_t1(407,'RA19CO2   ','NO        ','RA19NO3   ','          ','          ',  &
'          ',3.31e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(407,'RA19CO2   ','NO        ','RA19NO3   ','          ','          ',  &
'          ',3.726e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B408
ratb_t1(408,'RTN28O2   ','NO        ','RTN28NO3  ','          ','          ',  &
'          ',5.59e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(408,'RTN28O2   ','NO        ','RTN28NO3  ','          ','          ',  &
'          ',6.29e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B409
ratb_t1(409,'RTN25O2   ','NO        ','RTN25NO3  ','          ','          ',  &
'          ',3.84e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(409,'RTN25O2   ','NO        ','RTN25NO3  ','          ','          ',  &
'          ',4.32e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B410
ratb_t1(410,'RTX28O2   ','NO        ','RTX28NO3  ','          ','          ',  &
'          ',5.59e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(410,'RTX28O2   ','NO        ','RTX28NO3  ','          ','          ',  &
'          ',6.29e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B411
ratb_t1(411,'RTX24O2   ','NO        ','RTX24NO3  ','          ','          ',  &
'          ',3.77e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(411,'RTX24O2   ','NO        ','RTX24NO3  ','          ','          ',  &
'          ',4.24e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B412
ratb_t1(412,'RTX22O2   ','NO        ','RTX22NO3  ','          ','          ',  &
'          ',7.20e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(412,'RTX22O2   ','NO        ','RTX22NO3  ','          ','          ',  &
'          ',8.10e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B413
ratb_t1(413,'RN10O2    ','NO3       ','EtCHO     ','HO2       ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(413,'RN10O2    ','NO3       ','EtCHO     ','HO2       ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B414
ratb_t1(414,'RN13O2    ','NO3       ','MeCHO     ','EtOO      ','NO2       ',  &
'          ',9.95e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(414,'RN13O2    ','NO3       ','MeCHO     ','EtOO      ','NO2       ',  &
'          ',9.15e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B415
ratb_t1(415,'RN13O2    ','NO3       ','CARB11A   ','HO2       ','NO2       ',  &
'          ',1.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(415,'RN13O2    ','NO3       ','CARB11A   ','HO2       ','NO2       ',  &
'          ',1.38e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B416
ratb_t1(416,'RN16O2    ','NO3       ','RN15AO2   ','NO2       ','          ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(416,'RN16O2    ','NO3       ','RN15AO2   ','NO2       ','          ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B417
ratb_t1(417,'RN19O2    ','NO3       ','RN18AO2   ','NO2       ','          ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(417,'RN19O2    ','NO3       ','RN18AO2   ','NO2       ','          ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B418
ratb_t1(418,'RN13AO2   ','NO3       ','RN12O2    ','NO2       ','          ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(418,'RN13AO2   ','NO3       ','RN12O2    ','NO2       ','          ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B419
ratb_t1(419,'RN16AO2   ','NO3       ','RN15O2    ','NO2       ','          ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(419,'RN16AO2   ','NO3       ','RN15O2    ','NO2       ','          ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B420
ratb_t1(420,'RA13O2    ','NO3       ','CARB3     ','UDCARB8   ','HO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(420,'RA13O2    ','NO3       ','CARB3     ','UDCARB8   ','HO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B421
ratb_t1(421,'RA16O2    ','NO3       ','CARB3     ','UDCARB11  ','HO2       ',  &
'NO2       ',1.75e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(421,'RA16O2    ','NO3       ','CARB3     ','UDCARB11  ','HO2       ',  &
'NO2       ',1.61e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B422
ratb_t1(422,'RA16O2    ','NO3       ','CARB6     ','UDCARB8   ','HO2       ',  &
'NO2       ',7.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(422,'RA16O2    ','NO3       ','CARB6     ','UDCARB8   ','HO2       ',  &
'NO2       ',6.90e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B423
ratb_t1(423,'RA19AO2   ','NO3       ','CARB3     ','UDCARB14  ','HO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(423,'RA19AO2   ','NO3       ','CARB3     ','UDCARB14  ','HO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B424
ratb_t1(424,'RA19CO2   ','NO3       ','CARB9     ','UDCARB8   ','HO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(424,'RA19CO2   ','NO3       ','CARB9     ','UDCARB8   ','HO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B425
ratb_t1(425,'HOCH2CH2O2','NO3       ','HCHO      ','HCHO      ','HO2       ',  &
'NO2       ',1.94e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(425,'HOCH2CH2O2','NO3       ','HCHO      ','HCHO      ','HO2       ',  &
'NO2       ',1.78e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B426
ratb_t1(426,'HOCH2CH2O2','NO3       ','HOCH2CHO  ','HO2       ','NO2       ',  &
'          ',5.60e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(426,'HOCH2CH2O2','NO3       ','HOCH2CHO  ','HO2       ','NO2       ',  &
'          ',5.15e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B427
ratb_t1(427,'RN9O2     ','NO3       ','MeCHO     ','HCHO      ','HO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(427,'RN9O2     ','NO3       ','MeCHO     ','HCHO      ','HO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B428
ratb_t1(428,'RN12O2    ','NO3       ','MeCHO     ','MeCHO     ','HO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(428,'RN12O2    ','NO3       ','MeCHO     ','MeCHO     ','HO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B429
ratb_t1(429,'RN15O2    ','NO3       ','EtCHO     ','MeCHO     ','HO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(429,'RN15O2    ','NO3       ','EtCHO     ','MeCHO     ','HO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B430
ratb_t1(430,'RN18O2    ','NO3       ','EtCHO     ','EtCHO     ','HO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(430,'RN18O2    ','NO3       ','EtCHO     ','EtCHO     ','HO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B431
ratb_t1(431,'RN15AO2   ','NO3       ','CARB13    ','HO2       ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(431,'RN15AO2   ','NO3       ','CARB13    ','HO2       ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B432
ratb_t1(432,'RN18AO2   ','NO3       ','CARB16    ','HO2       ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(432,'RN18AO2   ','NO3       ','CARB16    ','HO2       ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B433
ratb_t1(433,'HOCH2CO3  ','NO3       ','HO2       ','HCHO      ','NO2       ',  &
'          ',4.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(433,'HOCH2CO3  ','NO3       ','HO2       ','HCHO      ','NO2       ',  &
'          ',3.68e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B434
ratb_t1(434,'RN8O2     ','NO3       ','MeCO3     ','HCHO      ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(434,'RN8O2     ','NO3       ','MeCO3     ','HCHO      ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B435
ratb_t1(435,'RN11O2    ','NO3       ','MeCO3     ','MeCHO     ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(435,'RN11O2    ','NO3       ','MeCO3     ','MeCHO     ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B436
ratb_t1(436,'RN14O2    ','NO3       ','EtCO3     ','MeCHO     ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(436,'RN14O2    ','NO3       ','EtCO3     ','MeCHO     ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B437
ratb_t1(437,'RN17O2    ','NO3       ','RN16AO2   ','NO2       ','          ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2 82
ratb_t1(437,'RN17O2    ','NO3       ','RN16AO2   ','NO2       ','          ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B438
ratb_t1(438,'RU14O2    ','NO3       ','UCARB12   ','HO2       ','NO2       ',  &
'          ',6.30e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(438,'RU14O2    ','NO3       ','UCARB12   ','HO2       ','NO2       ',  &
'          ',7.36e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B439
ratb_t1(439,'RU14O2    ','NO3       ','UCARB10   ','HCHO      ','HO2       ',  &
'NO2       ',1.87e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(439,'RU14O2    ','NO3       ','UCARB10   ','HCHO      ','HO2       ',  &
'NO2       ',2.23e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B440
ratb_t1(440,'RU12O2    ','NO3       ','MeCO3     ','HOCH2CHO  ','NO2       ',  &
'          ',1.75e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(440,'RU12O2    ','NO3       ','CARB6     ','HOCH2CHO  ','NO2       ',  &
'HO2       ',9.71e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B441
ratb_t1(441,'RU12O2    ','NO3       ','CARB7     ','CO        ','HO2       ',  &
'NO2       ',7.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(441,'RU12O2    ','NO3       ','CARB7     ','CARB3     ','HO2       ',  &
'NO2       ',9.246e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B442
ratb_t1(442,'RU10O2    ','NO3       ','MeCO3     ','HOCH2CHO  ','NO2       ',  &
'          ',1.25e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(442,'RU10O2    ','NO3       ','MeCO3     ','HOCH2CHO  ','NO2       ',  &
'          ',1.61e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B443
ratb_t1(443,'RU10O2    ','NO3       ','CARB6     ','HCHO      ','HO2       ',  &
'NO2       ',7.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(443,'RU10O2    ','NO3       ','CARB6     ','HCHO      ','HO2       ',  &
'NO2       ',6.9e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),      &
! B444
ratb_t1(444,'RU10O2    ','NO3       ','CARB7     ','HCHO      ','HO2       ',  &
'NO2       ',5.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(444,'RU10AO2   ','NO3       ','CARB7     ','HO2       ','CO        ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B445
ratb_t1(445,'NRN6O2    ','NO3       ','HCHO      ','HCHO      ','NO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
ratb_t1(445,'NRN6O2    ','NO3       ','HCHO      ','HCHO      ','NO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B446
ratb_t1(446,'NRN9O2    ','NO3       ','MeCHO     ','HCHO      ','NO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(446,'NRN9O2    ','NO3       ','MeCHO     ','HCHO      ','NO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119)      &
!
]


! Start and end bounds for 13th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+69
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B447
ratb_t1(447,'NRN12O2   ','NO3       ','MeCHO     ','MeCHO     ','NO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(447,'NRN12O2   ','NO3       ','MeCHO     ','MeCHO     ','NO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B448
ratb_t1(448,'NRU14O2   ','NO3       ','NUCARB12  ','HO2       ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(448,'NRU14O2   ','NO3       ','NUCARB12  ','HO2       ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B449
ratb_t1(449,'NRU12O2   ','NO3       ','NOA       ','CO        ','HO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(449,'NRU12O2   ','NO3       ','NOA       ','CO        ','HO2       ',  &
'NO2       ',1.15e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B450
ratb_t1(450,'RTN28O2   ','NO3       ','TNCARB26  ','HO2       ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(450,'RTN28O2   ','NO3       ','TNCARB26  ','HO2       ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B451
ratb_t1(451,'NRTN28O2  ','NO3       ','TNCARB26  ','NO2       ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(451,'NRTN28O2  ','NO3       ','TNCARB26  ','NO2       ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B452
ratb_t1(452,'RTN26O2   ','NO3       ','RTN25O2   ','NO2       ','          ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(452,'RTN26O2   ','NO3       ','RTN25O2   ','NO2       ','          ',  &
'          ',3.68e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B453
ratb_t1(453,'RTN25O2   ','NO3       ','RTN24O2   ','NO2       ','          ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(453,'RTN25O2   ','NO3       ','RTN24O2   ','NO2       ','          ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B454
ratb_t1(454,'RTN24O2   ','NO3       ','RTN23O2   ','NO2       ','          ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(454,'RTN24O2   ','NO3       ','RTN23O2   ','NO2       ','          ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B455
ratb_t1(455,'RTN23O2   ','NO3       ','Me2CO     ','RTN14O2   ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(455,'RTN23O2   ','NO3       ','Me2CO     ','RTN14O2   ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B456
ratb_t1(456,'RTN14O2   ','NO3       ','HCHO      ','TNCARB10  ','HO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(456,'RTN14O2   ','NO3       ','HCHO      ','TNCARB10  ','HO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B457
ratb_t1(457,'RTN10O2   ','NO3       ','RN8O2     ','CO        ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(457,'RTN10O2   ','NO3       ','RN8O2     ','CO        ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B458
ratb_t1(458,'RTX28O2   ','NO3       ','TXCARB24  ','HCHO      ','HO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(458,'RTX28O2   ','NO3       ','TXCARB24  ','HCHO      ','HO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B459
ratb_t1(459,'RTX24O2   ','NO3       ','TXCARB22  ','HO2       ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(459,'RTX24O2   ','NO3       ','TXCARB22  ','HO2       ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B460
ratb_t1(460,'RTX22O2   ','NO3       ','Me2CO     ','RN13O2    ','NO2       ',  &
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(460,'RTX22O2   ','NO3       ','Me2CO     ','RN13O2    ','NO2       ',  &
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B461
ratb_t1(461,'NRTX28O2  ','NO3       ','TXCARB24  ','HCHO      ','NO2       ',  &
'NO2       ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(461,'NRTX28O2  ','NO3       ','TXCARB24  ','HCHO      ','NO2       ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B462
ratb_t1(462,'RN10O2    ','HO2       ','RN10OOH   ','          ','          ',  &
'          ',1.51e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B463
ratb_t1(463,'RN13O2    ','HO2       ','RN13OOH   ','          ','          ',  &
'          ',1.82e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B464
ratb_t1(464,'RN16O2    ','HO2       ','RN16OOH   ','          ','          ',  &
'          ',2.05e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B465
ratb_t1(465,'RN19O2    ','HO2       ','RN19OOH   ','          ','          ',  &
'          ',2.24e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B466
ratb_t1(466,'RN13AO2   ','HO2       ','RN13OOH   ','          ','          ',  &
'          ',1.82e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B467
ratb_t1(467,'RN16AO2   ','HO2       ','RN16OOH   ','          ','          ',  &
'          ',2.05e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B468
ratb_t1(468,'RA13O2    ','HO2       ','RA13OOH   ','          ','          ',  &
'          ',2.24e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B469
ratb_t1(469,'RA16O2    ','HO2       ','RA16OOH   ','          ','          ',  &
'          ',2.39e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B470
ratb_t1(470,'RA19AO2   ','HO2       ','RA19OOH   ','          ','          ',  &
'          ',2.50e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B471
ratb_t1(471,'RA19CO2   ','HO2       ','RA19OOH   ','          ','          ',  &
'          ',2.50e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B472
ratb_t1(472,'HOCH2CH2O2','HO2       ','HOC2H4OOH ','          ','          ',  &
'          ',2.03e-13,  0.00,-1250.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B473
ratb_t1(473,'RN9O2     ','HO2       ','RN9OOH    ','          ','          ',  &
'          ',1.51e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B474
ratb_t1(474,'RN12O2    ','HO2       ','RN12OOH   ','          ','          ',  &
'          ',1.82e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B475
ratb_t1(475,'RN15O2    ','HO2       ','RN15OOH   ','          ','          ',  &
'          ',2.05e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B476
ratb_t1(476,'RN18O2    ','HO2       ','RN18OOH   ','          ','          ',  &
'          ',2.24e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B477
ratb_t1(477,'RN15AO2   ','HO2       ','RN15OOH   ','          ','          ',  &
'          ',2.05e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B478
ratb_t1(478,'RN18AO2   ','HO2       ','RN18OOH   ','          ','          ',  &
'          ',2.24e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B479
ratb_t1(479,'HOCH2CO3  ','HO2       ','HOCH2CO3H ','          ','          ',  &
'          ',4.30e-13,  0.00,-1040.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(479,'HOCH2CO3  ','HO2       ','HOCH2CO3H ','          ','          ',  &
'          ',2.912e-13,  0.00,-980.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B480
ratb_t1(480,'RN8O2     ','HO2       ','RN8OOH    ','          ','          ',  &
'          ',1.51e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B481
ratb_t1(481,'RN11O2    ','HO2       ','RN11OOH   ','          ','          ',  &
'          ',1.82e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B482
ratb_t1(482,'RN14O2    ','HO2       ','RN14OOH   ','          ','          ',  &
'          ',2.05e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B483
ratb_t1(483,'RN17O2    ','HO2       ','RN17OOH   ','          ','          ',  &
'          ',2.24e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B484
ratb_t1(484,'RU14O2    ','HO2       ','RU14OOH   ','          ','          ',  &
'          ',2.24e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(484,'RU14O2    ','HO2       ','RU14OOH   ','          ','          ',  &
'          ',2.05e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B485
ratb_t1(485,'RU12O2    ','HO2       ','RU12OOH   ','          ','          ',  &
'          ',2.05e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B486
ratb_t1(486,'RU10O2    ','HO2       ','RU10OOH   ','          ','          ',  &
'          ',1.82e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B487
ratb_t1(487,'NRN6O2    ','HO2       ','NRN6OOH   ','          ','          ',  &
'          ',1.13e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B488
ratb_t1(488,'NRN9O2    ','HO2       ','NRN9OOH   ','          ','          ',  &
'          ',1.51e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B489
ratb_t1(489,'NRN12O2   ','HO2       ','NRN12OOH  ','          ','          ',  &
'          ',1.82e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B490
ratb_t1(490,'NRU14O2   ','HO2       ','NRU14OOH  ','          ','          ',  &
'          ',2.24e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(490,'NRU14O2   ','HO2       ','NRU14OOH  ','          ','          ',  &
'          ',2.05e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B491
ratb_t1(491,'NRU12O2   ','HO2       ','NRU12OOH  ','          ','          ',  &
'          ',1.82e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(491,'NRU12O2   ','HO2       ','NRU12OOH  ','          ','          ',  &
'          ',2.05e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B492
ratb_t1(492,'RTN28O2   ','HO2       ','RTN28OOH  ','          ','          ',  &
'          ',2.66e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B493
ratb_t1(493,'NRTN28O2  ','HO2       ','NRTN28OOH ','          ','          ',  &
'          ',2.66e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B494
ratb_t1(494,'RTN26O2   ','HO2       ','RTN26OOH  ','          ','          ',  &
'          ',2.66e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(494,'RTN26O2   ','HO2       ','RTN26OOH  ','RTN25O2   ','OH        ',  &
'          ',5.20e-13,  0.00,-980.00, 0.56, 0.44, 0.44, 0.00,cs,0,0,119),      &
! B495
ratb_t1(495,'RTN25O2   ','HO2       ','RTN25OOH  ','          ','          ',  &
'          ',2.59e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B496
ratb_t1(496,'RTN24O2   ','HO2       ','RTN24OOH  ','          ','          ',  &
'          ',2.59e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107)      &
]

! Start and end bounds for 14th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+57
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B497
ratb_t1(497,'RTN23O2   ','HO2       ','RTN23OOH  ','          ','          ',  &
'          ',2.59e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B498
ratb_t1(498,'RTN14O2   ','HO2       ','RTN14OOH  ','          ','          ',  &
'          ',2.24e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B499
ratb_t1(499,'RTN10O2   ','HO2       ','RTN10OOH  ','          ','          ',  &
'          ',2.05e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B500
ratb_t1(500,'RTX28O2   ','HO2       ','RTX28OOH  ','          ','          ',  &
'          ',2.66e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B501
ratb_t1(501,'RTX24O2   ','HO2       ','RTX24OOH  ','          ','          ',  &
'          ',2.59e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B502
ratb_t1(502,'RTX22O2   ','HO2       ','RTX22OOH  ','          ','          ',  &
'          ',2.59e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B503
ratb_t1(503,'NRTX28O2  ','HO2       ','NRTX28OOH ','          ','          ',  &
'          ',2.66e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CRI RO2+RO2 reactions
! B504
ratb_t1(504,'MeOO      ','RO2       ','HCHO      ','HO2       ','          ',  &
'          ',6.01e-14,  0.00, -416.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(504,'MeOO      ','RO2       ','HCHO      ','HO2       ','          ',  &
'          ',2.06e-13, 0.00, -365.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,119),     &
! B505
ratb_t1(505,'MeOO      ','RO2       ','HCHO      ','          ','          ',  &
'          ',6.10e-14,  0.00, -416.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(505,'MeOO      ','RO2       ','MeOH      ','HCHO      ','          ',  &
'          ',2.06e-13, 0.00, -365.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,119),     &
! B506
ratb_t1(506,'MeOO      ','RO2       ','MeOH      ','          ','          ',  &
'          ',6.10e-14,  0.00, -416.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(506,'HOCH2CO3  ','HO2       ','HCHO      ','HO2       ','OH        ',  &
'          ',2.288e-13,  0.00,-980.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B507
ratb_t1(507,'EtOO      ','RO2       ','MeCHO     ','HO2       ','          ',  &
'          ',1.86e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B508
ratb_t1(508,'EtOO      ','RO2       ','MeCHO     ','          ','          ',  &
'          ',6.20e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B509
ratb_t1(509,'EtOO      ','RO2       ','EtOH      ','          ','          ',  &
'          ',6.20e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B510
ratb_t1(510,'RN10O2    ','RO2       ','EtCHO     ','HO2       ','          ',  &
'          ',3.60e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B511
ratb_t1(511,'RN10O2    ','RO2       ','EtCHO     ','          ','          ',  &
'          ',1.20e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B512
ratb_t1(512,'RN10O2    ','RO2       ','n-PrOH    ','          ','          ',  &
'          ',1.20e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B513
ratb_t1(513,'i-PrOO    ','RO2       ','Me2CO     ','HO2       ','          ',  &
'          ',2.40e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B514
ratb_t1(514,'i-PrOO    ','RO2       ','Me2CO     ','          ','          ',  &
'          ',8.00e-15,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B515
ratb_t1(515,'i-PrOO    ','RO2       ','i-PrOH    ','          ','          ',  &
'          ',8.00e-15,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B516
ratb_t1(516,'RN13O2    ','RO2       ','MeCHO     ','EtOO      ','          ',  &
'          ',9.95e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B517
ratb_t1(517,'RN13O2    ','RO2       ','CARB11A   ','HO2       ','          ',  &
'          ',1.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B518
ratb_t1(518,'RN13AO2   ','RO2       ','RN12O2    ','          ','          ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B519
ratb_t1(519,'RN16AO2   ','RO2       ','RN15O2    ','          ','          ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B520
ratb_t1(520,'RA13O2    ','RO2       ','CARB3     ','UDCARB8   ','HO2       ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B521
ratb_t1(521,'RA16O2    ','RO2       ','CARB3     ','UDCARB11  ','HO2       ',  &
'          ',6.16e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B522
ratb_t1(522,'RA16O2    ','RO2       ','CARB6     ','UDCARB8   ','HO2       ',  &
'          ',2.64e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B523
ratb_t1(523,'RA19AO2   ','RO2       ','CARB3     ','UDCARB14  ','HO2       ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B524
ratb_t1(524,'RA19CO2   ','RO2       ','CARB3     ','UDCARB14  ','HO2       ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B525
ratb_t1(525,'RN16O2    ','RO2       ','RN15AO2   ','          ','          ',  &
'          ',2.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B526
ratb_t1(526,'RN19O2    ','RO2       ','RN18AO2   ','          ','          ',  &
'          ',2.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B527
ratb_t1(527,'HOCH2CH2O2','RO2       ','HCHO      ','HCHO      ','HO2       ',  &
'          ',1.55e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B528
ratb_t1(528,'HOCH2CH2O2','RO2       ','HOCH2CHO  ','HO2       ','          ',  &
'          ',4.48e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B529
ratb_t1(529,'RN9O2     ','RO2       ','MeCHO     ','HCHO      ','HO2       ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B530
ratb_t1(530,'RN12O2    ','RO2       ','MeCHO     ','MeCHO     ','HO2       ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B531
ratb_t1(531,'RN15O2    ','RO2       ','EtCHO     ','MeCHO     ','HO2       ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B532
ratb_t1(532,'RN18O2    ','RO2       ','EtCHO     ','EtCHO     ','HO2       ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B533
ratb_t1(533,'RN15AO2   ','RO2       ','CARB13    ','HO2       ','          ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B534
ratb_t1(534,'RN18AO2   ','RO2       ','CARB16    ','HO2       ','          ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B535
ratb_t1(535,'MeCO3     ','RO2       ','MeOO      ','          ','          ',  &
'          ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B536
ratb_t1(536,'EtCO3     ','RO2       ','EtOO      ','          ','          ',  &
'          ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B537
ratb_t1(537,'HOCH2CO3  ','RO2       ','HCHO      ','HO2       ','          ',  &
'          ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B538
ratb_t1(538,'RN8O2     ','RO2       ','MeCO3     ','HCHO      ','          ',  &
'          ',1.40e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B539
ratb_t1(539,'RN11O2    ','RO2       ','MeCO3     ','MeCHO     ','          ',  &
'          ',1.40e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B540
ratb_t1(540,'RN14O2    ','RO2       ','EtCO3     ','MeCHO     ','          ',  &
'          ',1.40e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B541
ratb_t1(541,'RN17O2    ','RO2       ','RN16AO2   ','          ','          ',  &
'          ',1.40e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B542
ratb_t1(542,'RU14O2    ','RO2       ','UCARB12   ','HO2       ','          ',  &
'          ',4.31e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(542,'RU14O2    ','RO2       ','UCARB12   ','HO2       ','          ',  &
'          ',1.26e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,119),    &
! B543
ratb_t1(543,'RU14O2    ','RO2       ','UCARB10   ','HCHO      ','HO2       ',  &
'          ',1.28e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(543,'RU14O2    ','RO2       ','UCARB10   ','HCHO      ','HO2       ',  &
'          ',1.13e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,119),    &
! B544
ratb_t1(544,'RU12O2    ','RO2       ','MeCO3     ','HOCH2CHO  ','          ',  &
'          ',1.40e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(544,'RU12O2    ','RO2       ','MeCO3     ','HOCH2CHO  ','          ',  &
'          ',1.772e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,119),   &
! B545
ratb_t1(545,'RU12O2    ','RO2       ','CARB7     ','HOCH2CHO  ','HO2       ',  &
'          ',6.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(545,'RU12O2    ','RO2       ','HO2      ','CARB7     ','CARB3     ',   &
'          ',1.688e-13,  0.00,    0.00, 0.0, 0.0, 0.0, 0.0,cs,rp,0,119),       &
! B546
ratb_t1(546,'RU10O2    ','RO2       ','MeCO3     ','HOCH2CHO  ','          ',  &
'          ',1.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(546,'RU10O2    ','RO2       ','MeCO3     ','HOCH2CHO  ','          ',  &
'          ',1.28e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,119)     &
]

! Start and end bounds for 15th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+65
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B547
ratb_t1(547,'RU10O2    ','RO2       ','CARB6     ','HCHO      ','HO2       ',  &
'          ',6.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(547,'RU10O2    ','RO2       ','CARB6     ','HCHO      ','HO2       ',  &
'          ',5.49e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,119),    &
! B548
ratb_t1(548,'RU10O2    ','RO2       ','CARB7     ','HCHO      ','HO2       ',  &
'          ',4.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(548,'RU10AO2   ','RO2      ','CARB7     ','CO        ','HO2       ',   &
'          ',3.60e-13,  0.00,   0.00, 0.0, 0.0, 0.0, 0.0,cs,rp,0,119),         &
! B549
ratb_t1(549,'NRN6O2    ','RO2       ','HCHO      ','HCHO      ','NO2       ',  &
'          ',6.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B550
ratb_t1(550,'NRN9O2    ','RO2       ','MeCHO     ','HCHO      ','NO2       ',  &
'          ',2.30e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B551
ratb_t1(551,'NRN12O2   ','RO2       ','MeCHO     ','MeCHO     ','NO2       ',  &
'          ',2.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B552
ratb_t1(552,'NRU14O2   ','RO2       ','NUCARB12  ','HO2       ','          ',  &
'          ',1.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B553
ratb_t1(553,'NRU12O2   ','RO2       ','NOA       ','CO        ','HO2       ',  &
'          ',9.60e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(553,'NRU12O2   ','RO2       ','NOA       ','CO        ','HO2       ',  &
'          ',4.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00,0.00,cs,rp,0,119),     &
! B554
ratb_t1(554,'RTN28O2   ','RO2       ','TNCARB26  ','HO2       ','          ',  &
'          ',2.85e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B555
ratb_t1(555,'NRTN28O2  ','RO2       ','TNCARB26  ','NO2       ','          ',  &
'          ',1.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B556
ratb_t1(556,'RTN26O2   ','RO2       ','RTN25O2   ','          ','          ',  &
'          ',2.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! CS2
ratb_t1(556,'RTN26O2   ','RO2       ','RTN25O2   ','          ','          ',  &
'          ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,119),    &
! B557
ratb_t1(557,'RTN25O2   ','RO2       ','RTN24O2   ','          ','          ',  &
'          ',1.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B558
ratb_t1(558,'RTN24O2   ','RO2       ','RTN23O2   ','          ','          ',  &
'          ',6.70e-15,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B559
ratb_t1(559,'RTN23O2   ','RO2       ','Me2CO     ','RTN14O2   ','          ',  &
'          ',6.70e-15,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B560
ratb_t1(560,'RTN14O2   ','RO2       ','HCHO      ','TNCARB10  ','HO2       ',  &
'          ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B561
ratb_t1(561,'RTN10O2   ','RO2       ','RN8O2     ','CO        ','          ',  &
'          ',2.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B562
ratb_t1(562,'RTX28O2   ','RO2       ','TXCARB24  ','HCHO      ','HO2       ',  &
'          ',2.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B563
ratb_t1(563,'RTX24O2   ','RO2       ','TXCARB22  ','HO2       ','          ',  &
'          ',2.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B564
ratb_t1(564,'RTX22O2   ','RO2       ','Me2CO     ','RN13O2    ','          ',  &
'          ',2.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B565
ratb_t1(565,'NRTX28O2  ','RO2       ','TXCARB24  ','HCHO      ','NO2       ',  &
'          ',9.20e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,107),    &
! B566
ratb_t1(566,'OH        ','CARB14    ','RN14O2    ','          ','          ',  &
'          ',1.87e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B567
ratb_t1(567,'OH        ','CARB17    ','RN17O2    ','          ','          ',  &
'          ',4.36e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B568
ratb_t1(568,'OH        ','CARB11A   ','RN11O2    ','          ','          ',  &
'          ',1.50e-12,  0.00,   90.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B569
ratb_t1(569,'OH        ','CARB7     ','CARB6     ','HO2       ','          ',  &
'          ',3.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B570
ratb_t1(570,'OH        ','CARB10    ','CARB9     ','HO2       ','          ',  &
'          ',5.86e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B571
ratb_t1(571,'OH        ','CARB13    ','RN13O2    ','          ','          ',  &
'          ',1.65e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B572
ratb_t1(572,'OH        ','CARB16    ','RN16O2    ','          ','          ',  &
'          ',1.25e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B573
ratb_t1(573,'OH        ','UCARB10   ','RU10O2    ','          ','          ',  &
'          ',2.50e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(573,'OH        ','UCARB10   ','RU10O2    ','RU10AO2   ','MACO3     ',  &
'          ',3.84e-12,  0.00, -533.00, 0.693, 0.157, 0.15, 0.00,cs,0,0,119),   &
! B574
ratb_t1(574,'NO3       ','UCARB10   ','RU10O2    ','HONO2     ','          ',  &
'          ',1.44e-12,  0.00, 1862.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(574,'NO3       ','UCARB10   ','RU10O2    ','HONO2     ','          ',  &
'          ',5.98e-13,  0.00, 1862, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),        &
! B575
ratb_t1(575,'O3        ','UCARB10   ','HCHO      ','MeCO3     ','CO        ',  &
'OH        ',1.68e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(575,'O3        ','UCARB10   ','HCHO      ','MeCO3     ','CO        ',  &
'OH        ',3.84e-16,  0.00,  1710.0, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B576
ratb_t1(576,'O3        ','UCARB10   ','HCHO      ','CARB6     ','H2O2      ',  &
'          ',1.17e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(576,'O3        ','UCARB10   ','HCHO      ','CARB6     ','          ',  &
'          ',8.16e-16,  0.00,  1710.0, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B577
ratb_t1(577,'OH        ','HOCH2CHO  ','HOCH2CO3  ','          ','          ',  &
'          ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B578
ratb_t1(578,'NO3       ','HOCH2CHO  ','HOCH2CO3  ','HONO2     ','          ',  &
'          ',1.44e-12,  0.00, 1862.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(578,'NO3       ','HOCH2CHO  ','HOCH2CO3  ','HONO2     ','          ',  &
'          ',1.44e-12,  0.00, 1862.0, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),      &
! B579
ratb_t1(579,'OH        ','CARB3     ','CO        ','CO        ','HO2       ',  &
'          ',1.14e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(579,'OH        ','CARB3     ','CO        ','CO        ','HO2       ',  &
'          ',2.48e-12,  0.00,   -340.0, 0.0, 0.0, 0.0, 0.00,cs,0,0,119),       &
! B580
ratb_t1(580,'OH        ','CARB6     ','MeCO3     ','CO        ','          ',  &
'          ',1.72e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(580,'OH        ','CARB6     ','MeCO3     ','CO        ','          ',  &
'          ',1.9e-12,  0.00,   -575.0, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B581
ratb_t1(581,'OH        ','CARB9     ','RN9O2     ','          ','          ',  &
'          ',2.40e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B582
ratb_t1(582,'OH        ','CARB12    ','RN12O2    ','          ','          ',  &
'          ',1.38e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B583
ratb_t1(583,'OH        ','CARB15    ','RN15O2    ','          ','          ',  &
'          ',4.81e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B584
ratb_t1(584,'OH        ','CCARB12   ','RN12O2    ','          ','          ',  &
'          ',4.79e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B585
ratb_t1(585,'OH        ','UCARB12   ','RU12O2    ','          ','          ',  &
'          ',4.52e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(585,'OH        ','UCARB12   ','RU12O2    ','          ','          ',  &
'          ',6.42e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B586
ratb_t1(586,'NO3       ','UCARB12   ','RU12O2    ','HONO2     ','          ',  &
'          ',6.12e-12,  0.00, 1862.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(586,'NO3       ','UCARB12   ','RU12O2    ','HONO2     ','          ',  &
'          ',6.12e-12,  0.00, 1862.0, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),      &
! B587
ratb_t1(587,'O3        ','UCARB12   ','HOCH2CHO  ','MeCO3     ','CO        ',  &
'OH        ',2.14e-17,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(587,'O3        ','UCARB12   ','HOCH2CHO  ','MeCO3     ','CARB3     ',  &
'OH        ',6.00e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B588
ratb_t1(588,'O3        ','UCARB12   ','HOCH2CHO  ','CARB6     ','H2O2      ',  &
'          ',2.64e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(588,'O3        ','UCARB12   ','CARB6     ','CO        ','OH        ',  &
'HO2       ',1.20e-17,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B589
ratb_t1(589,'OH        ','NUCARB12  ','NRU12O2   ','          ','          ',  &
'          ',4.16e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B590
ratb_t1(590,'OH        ','NOA       ','CARB6     ','NO2       ','          ',  &
'          ',1.30e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(590,'OH        ','NOA       ','CARB6     ','NO2       ','          ',  &
'          ',6.70e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B591
ratb_t1(591,'OH        ','UDCARB8   ','EtOO      ','          ','          ',  &
'          ',2.60e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B592
ratb_t1(592,'OH        ','UDCARB11  ','RN10O2    ','          ','          ',  &
'          ',3.07e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B593
ratb_t1(593,'OH        ','UDCARB14  ','RN13O2    ','          ','          ',  &
'          ',3.85e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B594
ratb_t1(594,'OH        ','TNCARB26  ','RTN26O2   ','          ','          ',  &
'          ',4.20e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B595
ratb_t1(595,'OH        ','TNCARB15  ','RN15AO2   ','          ','          ',  &
'          ',1.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B596
ratb_t1(596,'OH        ','TNCARB10  ','RTN10O2   ','          ','          ',  &
'          ',1.00e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107)      &
]

! Start and end bounds for 16th section of ratb_defs.
! If you add extra reactions to this section, increment n_ratb_e
n_ratb_s = n_ratb_e+1
n_ratb_e = n_ratb_s+54
ratb_defs_master(n_ratb_s:n_ratb_e) = [                                        &
! B597
ratb_t1(597,'NO3       ','TNCARB26  ','RTN26O2   ','HONO2     ','          ',  &
'          ',3.80e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B598
ratb_t1(598,'NO3       ','TNCARB10  ','RTN10O2   ','HONO2     ','          ',  &
'          ',7.92e-12,  0.00, 1862.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B599
ratb_t1(599,'OH        ','RCOOH25   ','RTN25O2   ','          ','          ',  &
'          ',6.65e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B600
ratb_t1(600,'OH        ','TXCARB24  ','RTX24O2   ','          ','          ',  &
'          ',1.55e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B601
ratb_t1(601,'OH        ','TXCARB22  ','RTX22O2   ','          ','          ',  &
'          ',4.55e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B602
ratb_t1(602,'OH        ','EtONO2    ','MeCHO     ','NO2       ','          ',  &
'          ',4.40e-14,  0.00, -720.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B603
ratb_t1(603,'OH        ','RN10NO3   ','EtCHO     ','NO2       ','          ',  &
'          ',7.30e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B604
ratb_t1(604,'OH        ','i-PrONO2  ','Me2CO     ','NO2       ','          ',  &
'          ',4.90e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B605
ratb_t1(605,'OH        ','RN13NO3   ','CARB11A   ','NO2       ','          ',  &
'          ',9.20e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B606
ratb_t1(606,'OH        ','RN16NO3   ','CARB14    ','NO2       ','          ',  &
'          ',1.85e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B607
ratb_t1(607,'OH        ','RN19NO3   ','CARB17    ','NO2       ','          ',  &
'          ',3.02e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B608
ratb_t1(608,'OH        ','HOC2H4NO3 ','HOCH2CHO  ','NO2       ','          ',  &
'          ',1.09e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B609
ratb_t1(609,'OH        ','RN9NO3    ','CARB7     ','NO2       ','          ',  &
'          ',1.31e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B610
ratb_t1(610,'OH        ','RN12NO3   ','CARB10    ','NO2       ','          ',  &
'          ',1.79e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B611
ratb_t1(611,'OH        ','RN15NO3   ','CARB13    ','NO2       ','          ',  &
'          ',1.03e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B612
ratb_t1(612,'OH        ','RN18NO3   ','CARB16    ','NO2       ','          ',  &
'          ',1.34e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B613
ratb_t1(613,'OH        ','RU14NO3   ','UCARB12   ','NO2       ','          ',  &
'          ',5.55e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(613,'OH        ','RU14NO3   ','UCARB12   ','NO2       ','          ',  &
'          ',1.02e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B614
ratb_t1(614,'OH        ','RA13NO3   ','CARB3     ','UDCARB8   ','NO2       ',  &
'          ',7.30e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B615
ratb_t1(615,'OH        ','RA16NO3   ','CARB3     ','UDCARB11  ','NO2       ',  &
'          ',7.16e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B616
ratb_t1(616,'OH        ','RA19NO3   ','CARB6     ','UDCARB11  ','NO2       ',  &
'          ',8.31e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B617
ratb_t1(617,'OH        ','RTN28NO3  ','TNCARB26  ','NO2       ','          ',  &
'          ',4.35e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B618
ratb_t1(618,'OH        ','RTN25NO3  ','Me2CO     ','TNCARB15  ','NO2       ',  &
'          ',2.88e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B619
ratb_t1(619,'OH        ','RTX28NO3  ','TXCARB24  ','HCHO      ','NO2       ',  &
'          ',3.53e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B620
ratb_t1(620,'OH        ','RTX24NO3  ','TXCARB22  ','NO2       ','          ',  &
'          ',6.48e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B621
ratb_t1(621,'OH        ','RTX22NO3  ','Me2CO     ','CCARB12   ','NO2       ',  &
'          ',4.74e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B622
ratb_t1(622,'OH        ','AROH14    ','RAROH14   ','          ','          ',  &
'          ',2.63e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B623
ratb_t1(623,'NO3       ','AROH14    ','RAROH14   ','HONO2     ','          ',  &
'          ',3.78e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B624
ratb_t1(624,'RAROH14   ','NO2       ','ARNOH14   ','          ','          ',  &
'          ',2.08e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B625
ratb_t1(625,'OH        ','ARNOH14   ','CARB13    ','NO2       ','          ',  &
'          ',9.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B626
ratb_t1(626,'NO3       ','ARNOH14   ','CARB13    ','NO2       ','HONO2     ',  &
'          ',9.00e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B627
ratb_t1(627,'OH        ','AROH17    ','RAROH17   ','          ','          ',  &
'          ',4.65e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B628
ratb_t1(628,'NO3       ','AROH17    ','RAROH17   ','HONO2     ','          ',  &
'          ',1.25e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B629
ratb_t1(629,'RAROH17   ','NO2       ','ARNOH17   ','          ','          ',  &
'          ',2.08e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B630
ratb_t1(630,'OH        ','ARNOH17   ','CARB16    ','NO2       ','          ',  &
'          ',1.53e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B631
ratb_t1(631,'NO3       ','ARNOH17   ','CARB16    ','NO2       ','HONO2     ',  &
'          ',3.13e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B632
ratb_t1(632,'OH        ','RN10OOH   ','EtCHO     ','OH        ','          ',  &
'          ',1.89e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B633
ratb_t1(633,'OH        ','RN13OOH   ','CARB11A   ','OH        ','          ',  &
'          ',3.57e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B634
ratb_t1(634,'OH        ','RN16OOH   ','CARB14    ','OH        ','          ',  &
'          ',4.21e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B635
ratb_t1(635,'OH        ','RN19OOH   ','CARB17    ','OH        ','          ',  &
'          ',4.71e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B636
ratb_t1(636,'OH        ','EtCO3H    ','EtCO3     ','          ','          ',  &
'          ',4.42e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B637
ratb_t1(637,'OH        ','HOCH2CO3H ','HOCH2CO3  ','          ','          ',  &
'          ',6.19e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B638
ratb_t1(638,'OH        ','RN8OOH    ','CARB6     ','OH        ','          ',  &
'          ',4.42e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(638,'OH        ','RN8OOH    ','CARB6     ','OH        ','          ',  &
'          ',1.20e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B639
ratb_t1(639,'OH        ','RN11OOH   ','CARB9     ','OH        ','          ',  &
'          ',2.50e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B640
ratb_t1(640,'OH        ','RN14OOH   ','CARB12    ','OH        ','          ',  &
'          ',3.20e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B641
ratb_t1(641,'OH        ','RN17OOH   ','CARB15    ','OH        ','          ',  &
'          ',3.35e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B642
ratb_t1(642,'OH        ','RU14OOH   ','UCARB12   ','OH        ','          ',  &
'          ',7.51e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(642,'OH        ','RU14OOH   ','UCARB12   ','OH        ','          ',  &
'          ',6.66e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B643
ratb_t1(643,'OH        ','RU12OOH   ','RU12O2    ','          ','          ',  &
'          ',3.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(643,'OH        ','RU12OOH   ','RU10OOH   ','CO        ','HO2       ',  &
'          ',3.50e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B644
ratb_t1(644,'OH        ','RU10OOH   ','RU10O2    ','          ','          ',  &
'          ',3.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(644,'OH        ','RU10OOH   ','CARB7    ','CO        ','OH        ',   &
'          ',3.84e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B645
ratb_t1(645,'OH        ','NRU14OOH  ','NUCARB12  ','OH        ','          ',  &
'          ',1.03e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B646
ratb_t1(646,'OH        ','NRU12OOH  ','NOA       ','CO        ','OH        ',  &
'          ',2.65e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107)      &
]

! Start and end bounds for FINAL section of ratb_defs.
! This one must end on n_bimol_master!
! Remember to increment n_bimol_master if you add a reaction in any ratb array
n_ratb_s = n_ratb_e+1
ratb_defs_master(n_ratb_s:n_bimol_master) = [                                  &
! B647
ratb_t1(647,'OH        ','HOC2H4OOH ','HOCH2CHO  ','OH        ','          ',  &
'          ',2.13e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B648
ratb_t1(648,'OH        ','RN9OOH    ','CARB7     ','OH        ','          ',  &
'          ',2.50e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B649
ratb_t1(649,'OH        ','RN12OOH   ','CARB10    ','OH        ','          ',  &
'          ',3.25e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B650
ratb_t1(650,'OH        ','RN15OOH   ','CARB13    ','OH        ','          ',  &
'          ',3.74e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B651
ratb_t1(651,'OH        ','RN18OOH   ','CARB16    ','OH        ','          ',  &
'          ',3.83e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B652
ratb_t1(652,'OH        ','NRN6OOH   ','HCHO      ','HCHO      ','NO2       ',  &
'OH        ',5.22e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B653
ratb_t1(653,'OH        ','NRN9OOH   ','MeCHO     ','HCHO      ','NO2       ',  &
'OH        ',6.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B654
ratb_t1(654,'OH        ','NRN12OOH  ','MeCHO     ','MeCHO     ','NO2       ',  &
'OH        ',7.15e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B655
ratb_t1(655,'OH        ','RA13OOH   ','CARB3     ','UDCARB8   ','OH        ',  &
'          ',9.77e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B656
ratb_t1(656,'OH        ','RA16OOH   ','CARB3     ','UDCARB11  ','OH        ',  &
'          ',9.64e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B657
ratb_t1(657,'OH        ','RA19OOH   ','CARB6     ','UDCARB11  ','OH        ',  &
'          ',1.12e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B658
ratb_t1(658,'OH        ','RTN28OOH  ','TNCARB26  ','OH        ','          ',  &
'          ',2.38e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B659
ratb_t1(659,'OH        ','RTN26OOH  ','RTN26O2   ','          ','          ',  &
'          ',1.20e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B660
ratb_t1(660,'OH        ','NRTN28OOH ','TNCARB26  ','NO2       ','OH        ',  &
'          ',9.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B661
ratb_t1(661,'OH        ','RTN25OOH  ','RTN25O2   ','          ','          ',  &
'          ',1.66e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B662
ratb_t1(662,'OH        ','RTN24OOH  ','RTN24O2   ','          ','          ',  &
'          ',1.05e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B663
ratb_t1(663,'OH        ','RTN23OOH  ','RTN23O2   ','          ','          ',  &
'          ',2.05e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B664
ratb_t1(664,'OH        ','RTN14OOH  ','RTN14O2   ','          ','          ',  &
'          ',8.69e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B665
ratb_t1(665,'OH        ','RTN10OOH  ','RTN10O2   ','          ','          ',  &
'          ',4.23e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B666
ratb_t1(666,'OH        ','RTX28OOH  ','RTX28O2   ','          ','          ',  &
'          ',2.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B667
ratb_t1(667,'OH        ','RTX24OOH  ','TXCARB22  ','OH        ','          ',  &
'          ',8.59e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B668
ratb_t1(668,'OH        ','RTX22OOH  ','Me2CO     ','CCARB12   ','OH        ',  &
'          ',7.50e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B669
ratb_t1(669,'OH        ','NRTX28OOH ','NRTX28O2  ','          ','          ',  &
'          ',9.58e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B670
ratb_t1(670,'OH        ','PHAN      ','HCHO      ','CO        ','NO2       ',  &
'          ',1.12e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B671
ratb_t1(671,'OH        ','RU12PAN   ','UCARB10   ','NO2       ','          ',  &
'          ',2.52e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(671,'DHPCARB9  ','OH        ','RN8OOH    ','CO        ','OH        ',  &
'          ',3.64e-11,  0.00,   0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),      &
! B672
ratb_t1(672,'OH        ','RTN26PAN  ','Me2CO     ','CARB16    ','NO2       ',  &
'          ',3.66e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B673
ratb_t1(673,'OH        ','ANHY      ','HOCH2CH2O2','          ','          ',  &
'          ',1.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B674
ratb_t1(674,'OH        ','UDCARB8   ','ANHY      ','HO2       ','          ',  &
'          ',2.60e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B675
ratb_t1(675,'OH        ','UDCARB11  ','ANHY      ','MeOO      ','          ',  &
'          ',2.51e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B676
ratb_t1(676,'OH        ','UDCARB14  ','ANHY      ','EtOO      ','          ',  &
'          ',3.15e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! B677
ratb_t1(677,'RTN23O2   ','NO        ','RTN23NO3  ','          ','          ',  &
'          ',2.83e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(677,'RU12O2    ','NO        ','RU12NO3   ','         ','          ',   &
'          ',9.45e-14,  0.00, -360, 0.00, 0.00, 0.00 ,0.00, cs,0,0,119),       &
! B678
ratb_t1(678,'RTN23NO3  ','OH        ','Me2CO     ','TNCARB12  ','NO2       ',  &
'          ',5.37e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(678,'RU14NO3   ','OH        ','RU10NO3   ','HCHO      ','HO2       ',  &
'          ',1.44e-11,  0.00,   0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),      &
! B679
ratb_t1(679,'TNCARB12  ','OH        ','TNCARB11  ','HO2       ','          ',  &
'          ',3.22e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS
ratb_t1(679,'RU12NO3   ','OH        ','CARB7     ','CARB3     ','NO2       ',  &
'          ',2.50e-12,  0.00, 0.00, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),            &
! B680
ratb_t1(680,'TNCARB11  ','OH        ','RTN10O2   ','CO        ','          ',  &
'          ',1.33e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(680,'IEPOX     ','OH        ','RU12O2    ','          ','          ',  &
'          ',1.16e-11,  0.00,00.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),        &
! B681
ratb_t1(681,'TNCARB11  ','NO3       ','RTN10O2   ','CO        ','HONO2     ',  &
'          ',7.92e-12,  0.00, 1862.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,107),     &
! CS2
ratb_t1(681,'HMML      ','OH        ','CARB6     ','OH        ','HCOOH     ',  &
'MeCO3     ',4.33e-12,  0.00,00.00, 0.70, 0.70, 0.30, 0.30,cs,0,0,119),        &
! DMS scheme in CRI is quite different to rest of UKCA, so putting here
! B682
ratb_t1(682,'DMS       ','OH        ','MeSCH2OO  ','H2O       ','          ',  &
'          ',1.12e-11,  0.00,  250.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B683
ratb_t1(683,'DMS       ','NO3       ','MeSCH2OO  ','HONO2     ','          ',  &
'          ',1.90e-13,  0.00, -520.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B684
ratb_t1(684,'MeSCH2OO  ','NO        ','HCHO      ','MeS       ','NO2       ',  &
'          ',4.90e-12,  0.00, -263.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B685
ratb_t1(685,'MeSCH2OO  ','MeSCH2OO  ','HCHO      ','HCHO      ','MeS       ',  &
'MeS       ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B686
ratb_t1(686,'MeS       ','O3        ','MeSO      ','          ','          ',  &
'          ',1.15e-12,  0.00, -432.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B687
ratb_t1(687,'MeS       ','NO2       ','MeSO      ','NO        ','          ',  &
'          ',3.00e-11,  0.00, -210.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B688
ratb_t1(688,'MeSO      ','NO2       ','MeSO2     ','SO2       ','MeOO      ',  &
'NO        ',1.20e-11,  0.00,    0.00, 0.82, 0.18, 0.18, 1.00,cs,a,0,107),     &
! B689
ratb_t1(689,'MeSO      ','O3        ','MeSO2     ','          ','          ',  &
'          ',6.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B690
ratb_t1(690,'MeSO2     ','          ','SO2       ','MeOO      ','          ',  &
'          ',5.00e+13,  0.00, 9673.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B691
ratb_t1(691,'MeSO2     ','NO2       ','MeSO3     ','NO        ','          ',  &
'          ',2.20e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B692
ratb_t1(692,'MeSO2     ','O3        ','MeSO3     ','          ','          ',  &
'          ',3.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B693
ratb_t1(693,'MeSO3     ','HO2       ','MSA       ','          ','          ',  &
'          ',5.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B694
ratb_t1(694,'MeSO3     ','          ','MeOO      ','H2SO4     ','          ',  &
'          ',1.36e+14,  0.00,11071.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B695
ratb_t1(695,'DMSO      ','OH        ','MSIA      ','MeOO      ','          ',  &
'          ',8.70e-11,  0.00,    0.00, 0.95, 0.95, 0.00, 0.00,cs,a,0,107),     &
! B696
ratb_t1(696,'MSIA      ','OH        ','MeSO2     ','MSA       ','HO2       ',  &
'H2O       ',9.00e-11,  0.00,    0.00, 0.95, 0.05, 0.05, 1.00,cs,a,0,107),     &
! B697
ratb_t1(697,'MSIA      ','NO3       ','MeSO2     ','HONO2     ','          ',  &
'          ',1.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,a,0,107),     &
! B698
ratb_t1(698,'DHPR12O2  ','NO        ','CARB3     ','RN8OOH    ','OH        ',  &
'NO2       ',2.70e-12,  0.00,    -360, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B699
ratb_t1(699,'DHPR12O2  ','NO3       ','CARB3     ','RN8OOH    ','OH        ',  &
'NO2       ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B700
ratb_t1(700,'DHPR12O2  ','HO2       ','DHPR12OOH ','          ','          ',  &
'          ',2.054e-13,  0.00,   -1300, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),    &
! B701
ratb_t1(701,'DHPR12O2  ','RO2       ','CARB3     ','RN8OOH    ','OH        ',  &
'          ',7.60e-13,  0.00,   0.00, 0.00, 0.00, 0.00, 0.00,cs,rp,0,119),     &
! B702
ratb_t1(702,'DHPR12OOH ','OH        ','DHPCARB9  ','CO        ','OH        ',  &
'          ',5.64e-11,  0.00,   0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),      &
! B703
ratb_t1(703,'HPUCARB12 ','OH        ','HUCARB9   ','CO        ','OH        ',  &
'          ',5.20e-11,  0.00,   0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),      &
! B704
ratb_t1(704,'HPUCARB12 ','O3        ','CARB3     ','CARB6     ','OH        ',  &
'OH        ',2.4e-17,  0.00,   0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),       &
! B705
ratb_t1(705,'HPUCARB12 ','NO3       ','HUCARB9   ','CO        ','OH        ',  &
'HONO2     ',6.12e-12,  0.00,1862.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),      &
! B706
ratb_t1(706,'HUCARB9   ','OH        ','CARB6     ','CO        ','HO2       ',  &
'          ',5.78e-11,  0.00,00.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),        &
! B707
ratb_t1(707,'MPAN      ','OH        ','HMML      ','NO3       ','          ',  &
'          ',2.262e-11,  0.00,00.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),       &
! B708
ratb_t1(708,'RU10AO2   ','NO        ','RU10NO3   ','          ','          ',  &
'          ',3.51e-14,  0.00,-360.00, 0.013, 0.0, 0.0, 0.0,cs,0,0,119),        &
! B709
ratb_t1(709,'RU10AO2   ','NO        ','CARB7     ','HO2       ','CO        ',  &
'NO2       ',2.665e-12,  0.00,-360.00, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),         &
! B710
ratb_t1(710,'RU10AO2   ','HO2       ','RU10OOH   ','          ','          ',  &
'          ',1.819e-13,  0.00,   -1300.00, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),     &
! B711
ratb_t1(711,'OH        ','CARB3     ','CO        ','OH        ','          ',  &
'          ',6.2e-13, 0.0, -340.0, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),             &
! B712
ratb_t1(712,'DHCARB9   ','OH        ','CARB6     ','HO2       ','          ',  &
'          ',3.42e-11,  0.00,   0.00, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),          &
! B713
ratb_t1(713,'RU12O2    ','NO        ','CARB7     ','HOCH2CO3  ','NO2       ',  &
'          ',4.59e-13,  0.00, -360.00, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),         &
! B714
ratb_t1(714,'RU10NO3   ','OH        ','CARB7     ','CO        ','NO2       ',  &
'          ',5.26e-13,  0.00, 0.00, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),            &
! B715
ratb_t1(715,'RU12O2    ','NO3       ','CARB7     ','HOCH2CO3  ','NO2       ',  &
'          ',4.05e-13,  0.00, 0.00, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),            &
! B716
ratb_t1(716,'MACO3     ','NO        ','MeOO      ','CO        ','HCHO      ',  &
'NO2       ',2.438e-12, 0, -290, 2.0, 2.0, 2.0, 2.0,cs,0,0,119),               &
! B717
ratb_t1(717,'MACO3     ','NO        ','HO2       ','          ','          ',  &
'          ',2.438e-12, 0.0,-290,2.0, 0.0, 0.0, 0.0,cs,0,0,119),               &
! B718
ratb_t1(718,'MACO3     ','NO        ','MeCO3     ','HCHO      ','HO2       ',  &
'NO2       ',2.625e-12, 0.0, -290, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),             &
! B719
ratb_t1(719,'MACO3     ','NO3       ','MeOO      ','HCHO      ','HO2       ',  &
'CO        ',1.3e-12, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0,cs,0,0,119),                &
! B720
ratb_t1(720,'MACO3     ','NO3       ','NO2       ','          ','          ',  &
'          ',1.3e-12, 0.0, 0.00, 2.0, 0.0, 0.0, 0.0,cs,0,0,119),               &
! B721
ratb_t1(721,'MACO3     ','NO3       ','MeCO3     ','HCHO      ','HO2       ',  &
'NO2       ',1.4e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),                &
! B722
ratb_t1(722,'MACO3     ','HO2       ','RU10OOH   ','          ','          ',  &
'          ',2.91e-13, 0.0, -980, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),              &
! B723
ratb_t1(723,'MACO3     ','HO2       ','MeOO      ','CO        ','HCHO      ',  &
'OH        ',2.29e-13, 0.0, -980, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),              &
! B724
ratb_t1(724,'MACO3     ','RO2       ','MeOO      ','CO        ','HCHO      ',  &
'OH        ',6.5e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),                &
! B725
ratb_t1(725,'MACO3     ','RO2       ','MeCO3     ','HCHO      ','HO2       ',  &
'          ',3.5e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,cs,0,0,119),                &
! B726
ratb_t1(726,'NRU12O2   ','NO        ','NOA       ','CARB3     ','HO2       ',  &
'NO2       ',1.35e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B727
ratb_t1(727,'NRU12O2   ','NO3       ','NOA       ','CARB3     ','HO2       ',  &
'NO2       ',1.15e-12,  0.00, 0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),        &
! B728
ratb_t1(728,'O3        ','UCARB12   ','CARB3     ','CARB7     ','CO        ',  &
'          ',6.00e-18,  0.00, 0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),        &
! B729
ratb_t1(729,'O3        ','C5H8      ','UCARB10   ','HCOOH     ','         ',   &
'          ',1.80e-15,  0.00, 1995, 0.00, 0.00,0.00,0.00,cs,0,0,119),          &
! B730
ratb_t1(730,'O3        ','C5H8      ','UCARB10   ','HCHO     ','H2O2      ',   &
'          ',3.97e-15,  0.00, 1995, 0.00, 0.00,0.00 ,0.00,cs,0,0,119),         &
! B731
ratb_t1(731,'O3        ','C5H8      ','MeCO3     ','HCHO      ','CO        ',  &
'OH        ',1.29e-15,  0.00, 1995, 1.00, 2.00,1.00 ,1.00, cs,0,0,119),        &
! B732
ratb_t1(732,'RU14O2    ','NO        ','HPUCARB12 ','HO2       ','DHPR12O2  ',  &
'NO        ',1.399e-15,  0.00,  -1668.00, 0.50, 0.50, 0.50, 1.00,cs,0,0,119),  &
! B733
ratb_t1(733,'RU14O2    ','NO3       ','HPUCARB12 ','HO2       ','DHPR12O2  ',  &
'NO3       ',1.191e-15,  0.00,  -1308.00, 0.50, 0.50, 0.50, 1.00,cs,0,0,119),  &
! B734
ratb_t1(734,'RU14O2    ','HO2       ','HPUCARB12 ','HO2       ','DHPR12O2  ',  &
'HO2       ',1.064e-16,  0.00, -2608.00, 0.5, 0.5,0.5 ,1.0, cs,0,0,119),       &
! B735
ratb_t1(735,'RU14O2    ','RO2       ','HPUCARB12 ','HO2       ','DHPR12O2  ',  &
'RO2       ',6.527e-16,  0.00,  -1308.00, 0.50, 0.50, 0.50, 1.00,cs,rp,0,119), &
! B736
ratb_t1(736,'NRU12O2   ','RO2       ','NOA       ','CARB3     ','HO2       ',  &
'          ',4.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00,0.00,cs,rp,0,119),     &
! B737
ratb_t1(737,'HOCH2CO3  ','HO2       ','HCHO      ','HO2       ','OH        ',  &
'          ',2.288e-13,  0.00,-980.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B738
ratb_t1(738,'O3        ','C5H8      ','UCARB10   ','CO        ','          ',  &
'          ',9.785e-16,  0.00, 1995.00, 0.00, 0.00, 0.00, 0.00 ,cs,0,0,119),   &
! B739
ratb_t1(739,'RU14NO3   ','OH       ','RU12NO3   ','HO2       ','          ',   &
'          ',5.4e-12,  0.00,   0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),       &
! B740
ratb_t1(740,'RU12O2    ','RO2       ','CARB7    ','HOCH2CO3  ','          ',   &
'          ',7.392e-14,  0.00,    0.00, 0.0, 0.0, 0.0, 0.0,cs,rp,0,119),       &
! B741
ratb_t1(741,'OH        ','RU14OOH   ','RU14O2    ','          ','          ',  &
'          ',4.44e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B742
ratb_t1(742,'OH        ','RU14OOH   ','IEPOX     ','OH        ','          ',  &
'          ',6.29e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,cs,0,0,119),     &
! B743 c.f. ~B061
ratb_t1(743,'HO2       ','NO3       ','O2        ','HONO2     ','          ',  &
'          ',2.15e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,s,0,0,107),      &
! B744 c.f. ~B135
ratb_t1(744,'O3        ','C5H8      ','OH        ','          ','          ',  &
'          ',3.93e-15,  0.00, 1913.00, 0.54, 0.00, 0.00, 0.00,r,0,0,107),      &
! B745 - updates for 121
ratb_t1(745,'HO2       ','EtCO3     ','OH        ','EtOO      ','CO2       ',  &
'          ',1.89e-13,  0.00,  -1040.00, 0.00, 0.00, 0.00, 0.00,st,0,0,121)    &
]

!----------------------------------------------------------------------
! NOTES: CheST Bimolecular Reactions
!----------------------------------------------------------------------
! B001 Br+Cl2O2 -> BrCl Cl O2 JPL2011
! B001 k(298)=3.3E-12. Both JPL2011 and IUPAC2006 agree. Reaction forms
! B001 ClOO rather than Cl and O2.
!----------------------------------------------------------------------
! B002 Br+HCHO -> HBr CO HO2 JPL2011
! B002 k(298)=1.1E-12. Unchanged since JPL2009. No reference by IUPAC.
! B002 JPL cite two studies of which this value is derived from by a
! B002 weighted mean.
!----------------------------------------------------------------------
! B003 Br+HO2 -> HBr O2  JPL2011
! B003 k(298)=1.7E-12. Unchanged since JPL2006. IUPAC2003 recommend
! B003 k=7.7E-12*exp(-450/temp).
!----------------------------------------------------------------------
! B004 Br+O3 -> BrO O2  JPL2011
! B004 k(298)=1.2E-12. Unchanged since JPL2009. IUPAC2003 recommend k=1
! B004 .7E-11exp(-800/temp)
!----------------------------------------------------------------------
! B005 Br+OClO -> BrO ClO  JPL2011
! B005 k(298)=3.4E-13. Originally from JPL1990
!----------------------------------------------------------------------
! B006 BrO+BrO -> Br Br O2 JPL2011
! B006 Originally from JPL1997
!----------------------------------------------------------------------
! B007 BrO+ClO -> Br Cl O2 JPL2011
! B007 k(298)=5.5E-12. Unchanged since JPL2009. IUPAC2004 recommend
! B007 k=2.9E-12*exp(220/temp) with k(298)=6.1E-12. Both JPL and IUPAC
! B007 recommended products are Br + ClOO.
!----------------------------------------------------------------------
! B008 BrO+ClO -> Br OClO  JPL2011
! B008 k(298)=6.0E-12. Unchanged since JPL2009. IUPAC2004 recommend
! B008 k=1.60E-13*exp(430/temp) with k(298)=6.8E-12 (at 220 K k =
! B008 1.13E-12 c.f 1.16E-12). This channel accounts for around 60% of
! B008 the total for this reaction.
!----------------------------------------------------------------------
! B009 BrO+ClO -> BrCl O2  JPL2011
! B009 k(298)=1.1E-12. Unchanged since JPL2009. IUPAC2004 recommend
! B009 k=5.8E-13*exp(170/temp) with k(298)=1.0E-12. Branching ratios put
! B009 this channel at ~8%.
!----------------------------------------------------------------------
! B010 BrO+HO2 -> HBr O3  JPL2011
! B010 Remove this reaction since both JPL2006 and IUPAC2003 suggest
! B010 that the branching ratio for this channel is 0.0
!----------------------------------------------------------------------
! B011 BrO+HO2 -> HOBr O2  JPL2011
! B011 k(298)=2.1E-11. Unchanged since JPL2006. IUPAC2003 recommend k=4
! B011 .5E-12exp(500/temp) with k(298)=2.4E-11. Both JPL and IUPAC
! B011 suggest that essentially this is the only channel for this
! B011 reaction.
!----------------------------------------------------------------------
! B012 BrO+NO -> Br NO2  JPL2011
! B012 k(298)=2.1E-11. Originally from JPL1982
!----------------------------------------------------------------------
! B013 BrO+OH -> Br HO2  JPL2011
! B013 k(298)=3.9E-11. Unchanged since JPL2006. IUPAC2003 recommend
! B013 k=1.8E-11*exp(250/temp) with k(298)=4.1E-11. The temperature
! B013 dependence comes from the study by Bedjanian et al.
!----------------------------------------------------------------------
! B014 CF2Cl2+O(1D) -> Cl ClO  JPL2011
! B014 k(298)=1.4E-10. Unchanged since JPL1992. Both JPL1992 and
! B014 IUPAC2011 agree. The ClO channel is ~87%.
!----------------------------------------------------------------------
! B015 CFCl3+O(1D) -> Cl Cl ClO JPL2011
! B015 k(298)=2.3E-10. Unchanged since JPL1992. Both JPL1992 and
! B015 IUPAC2011 agree. The ClO channel is ~80%.
!----------------------------------------------------------------------
! B016 Cl+CH4 -> HCl MeOO  JPL2011
! B016 k(298)=1.0E-13. IUPAC recommend k=6.6E-12*exp(-1240/temp).
!----------------------------------------------------------------------
! B017 Cl+Cl2O2 -> Cl Cl Cl IUPAC2006
! B017 k(298)=1.0E-10. Probable error in JPL2011 gives beta=9.4E-11!
! B017 IUPAC2006 recommend k = 7.6E-11*exp(65/temp) based on study by
! B017 Ingham et al.
!----------------------------------------------------------------------
! B018 Cl+ClONO2 -> Cl Cl NO3 JPL2011
! B018 IUPAC similar (6.2 and -145). Actual reaction forms Cl2 rather
! B018 than 2Cl.
!----------------------------------------------------------------------
! B019 Cl+H2 -> HCl H  JPL2011
! B019 IUPAC similar (3.9 and 2310)
!----------------------------------------------------------------------
! B020 Cl+H2O2 -> HCl HO2  JPL2011
! B020 IUPAC identical.
!----------------------------------------------------------------------
! B021 Cl+HCHO -> HCl CO HO2 JPL2011
! B021 IUPAC similar (8.1 and 34).
!----------------------------------------------------------------------
! B022 Cl+HO2 -> ClO OH  JPL2011
! B022 IUPAC gives total rate + k2/k. Hard to compare
!----------------------------------------------------------------------
! B023 Cl+HO2 -> HCl O2  JPL2011
! B023 IUPAC gives total rate + k2/k. Hard to compare
!----------------------------------------------------------------------
! B024 Cl+HOCl -> Cl Cl OH JPL2011
! B024 No IUPAC rate
!----------------------------------------------------------------------
! B025 Cl+MeOOH -> HCl MeOO  JPL2011
! B025 IUPAC rate similar (5.9E-11)
!----------------------------------------------------------------------
! B026 Cl+NO3 -> ClO NO2  JPL2011
! B026 IUPAC rates identical
!----------------------------------------------------------------------
! B027 Cl+O3 -> ClO O2  JPL2011
! B027 IUPAC rates similar (2.8 and 250)
!----------------------------------------------------------------------
! B028 Cl+OClO -> ClO ClO  JPL2011
! B028 IUPAC rates similar (3.2 and -170)
!----------------------------------------------------------------------
! B029 ClO+ClO -> Cl Cl O2 JPL2011
! B029 IUPAC rates identical
!----------------------------------------------------------------------
! B030 ClO+ClO -> Cl Cl O2 JPL2011
! B030 IUPAC rates identical
!----------------------------------------------------------------------
! B031 ClO+ClO -> Cl OClO  JPL2011
! B031 IUPAC rates identical
!----------------------------------------------------------------------
! B032 ClO+HO2 -> HCl O3  JPL2011
! B032 No rate given in IUPAC or JPL. Reaction should probably be
! B032 removed.
!----------------------------------------------------------------------
! B033 ClO+HO2 -> HOCl O2  JPL2011
! B033 IUPAC rates similar (2.2 and -340)
!----------------------------------------------------------------------
! B034 ClO+MeOO -> Cl HCHO HO2 JPL2011
! B034 IUPAC provides a more detailed breakdown of products. Would
! B034 require a modification of asad_bimol to treat properly. Plus
! B034 retain policy of using JPL for ClOx reactions.
!----------------------------------------------------------------------
! B035 ClO+NO -> Cl NO2  JPL2011
! B035 IUPAC rates very similar (6.2 and -295)
!----------------------------------------------------------------------
! B036 ClO+NO3 -> Cl O2 NO2 JPL2011
! B036 Assume all ClO + NO3 forms ClOO. Note IUPAC still has 1.2/4.6
! B036 going to other channel
!----------------------------------------------------------------------
! B037 ClO+NO3 -> OClO NO2  JPL2011
! B037 Assume all ClO + NO3 forms ClOO. Note IUPAC still has 1.2/4.6
! B037 going to other channel. If use JPL approach reaction could be
! B037 removed.
!----------------------------------------------------------------------
! B038 EtCO3+NO -> EtOO CO2 NO2 IUPAC2002
! B038 No JPL rate
!----------------------------------------------------------------------
! B039 EtCO3+NO3 -> EtOO CO2 NO2 MCMv3.2
! B039 No IUPAC/JPL rate
!----------------------------------------------------------------------
! B040 EtOO+MeCO3 -> MeCHO HO2 MeOO IUPAC2002
! B040 No JPL rate
!----------------------------------------------------------------------
! B041 EtOO+NO -> MeCHO HO2 NO2 IUPAC2005
! B041 No JPL rate
!----------------------------------------------------------------------
! B042 EtOO+NO3 -> MeCHO HO2 NO2 IUPAC2008
! B042 No JPL rate
!----------------------------------------------------------------------
! B043 H+HO2 -> H2 O2  JPL2011
! B043 No IUPAC rate
!----------------------------------------------------------------------
! B044 H+HO2 -> O(3P) H2O  JPL2011
! B044 No IUPAC rate
!----------------------------------------------------------------------
! B045 H+HO2 -> OH OH  JPL2011
! B045 No IUPAC rate
!----------------------------------------------------------------------
! B046 H+NO2 -> OH NO  JPL2011
! B046 No IUPAC rate
!----------------------------------------------------------------------
! B047 H+O3 -> OH O2  JPL2011
! B047 No IUPAC rate
!----------------------------------------------------------------------
! B048 HO2+EtCO3 -> O2 EtCO3H  MCMv3.2*
! B048 No IUPAC/JPL rate. No exact correspondence in MCM. Take total
! B048 rate and subtract O3 channel.
!----------------------------------------------------------------------
! B049 HO2+EtCO3 -> O3 EtCO2H  MVMv3.2
! - NOTE EtCO2H IS *NOT* CURRENTLY A SPECIES CONSIDERED BY ASAD
! B049 No IUPAC/JPL rate
!----------------------------------------------------------------------
! B050 HO2+EtOO -> EtOOH   IUPAC2011
! B050 -
!----------------------------------------------------------------------
! B051 HO2+HO2 -> H2O2   JPL2011 see also asad_bimol
! B051 IUPAC slightly different (2.2 & -600) . Rate modified in the
! B051 presence of H2O in asad bimol.
!----------------------------------------------------------------------
! B052 HO2+ISO2 -> ISOOH   Poschl00
! B052 -
!----------------------------------------------------------------------
! B053 HO2+MACRO2 -> MACROOH   Poschl00
! B053 -
!----------------------------------------------------------------------
! B054 HO2+MeCO3 -> MeCO2H O3  IUPAC2009
! B054 -
!----------------------------------------------------------------------
! B055 HO2+MeCO3 -> MeCO3H   IUPAC2009
! B055 -
!----------------------------------------------------------------------
! B056 HO2+MeCO3 -> OH MeOO  IUPAC2009
! B056 -
!----------------------------------------------------------------------
! B057 HO2+MeCOCH2OO -> MeCOCH2OOH   IUPAC2009
! B057 Use measured value rather than the expression from the MCM (which
! B057 does include T dependence)
!----------------------------------------------------------------------
! B058 HO2+MeOO -> HCHO H2O  IUPAC2009 see also asad_bimol
! B058 Total rate given. Split between channels controlled by asad_bimol
!----------------------------------------------------------------------
! B059 HO2+MeOO -> MeOOH O2  IUPAC2009 see also asad_bimol
! B059 Total rate given. Split between channels controlled by asad_bimol
!----------------------------------------------------------------------
! B060 HO2+NO -> OH NO2  JPL2011
! B060 IUPAC values similar (3.45 and 270)
!----------------------------------------------------------------------
! B061 HO2+NO3 -> OH NO2 O2 JPL2011
! B061 IUPAC value slightly larger (4 cf 3.5)
!----------------------------------------------------------------------
! B062 HO2+O3 -> OH O2  IUPAC2001
! B062 Use more sophisticated expression of IUPAC. JPL have simple term
! B062 1*10^-14 and 490.
!----------------------------------------------------------------------
! B063 HO2+i-PrOO -> i-PrOOH   MCMv3.2
! B063 No rate in IUPAC or JPL
!----------------------------------------------------------------------
! B064 HO2+n-PrOO -> n-PrOOH   MCMv3.2
! B064 No rate in IUPAC or JPL
!----------------------------------------------------------------------
! B065 ISO2+ISO2 -> MACR MACR HCHOHO2 Poschl00
! B065 -
!----------------------------------------------------------------------
! B066 MACRO2+MACRO2 -> HACET MGLY HCHOCO Poschl00
! B066 -
!----------------------------------------------------------------------
! B067 MACRO2+MACRO2 -> HO2   Poschl00
! B067 -
!----------------------------------------------------------------------
! B068 MeBr+Cl -> Br HCl  JPL2011
! B068 -
!----------------------------------------------------------------------
! B069 MeBr+O(1D) -> Br OH  JPL2011
! B069 -
!----------------------------------------------------------------------
! B070 MeBr+OH -> Br H2O  JPL2011
! B070 -
!----------------------------------------------------------------------
! B071 MeCO3+NO -> MeOO CO2 NO2 IUPAC2002
! B071 -
!----------------------------------------------------------------------
! B072 MeCO3+NO3 -> MeOO CO2 NO2 IUPAC2008
! B072 -
!----------------------------------------------------------------------
! B073 MeCOCH2OO+NO -> MeCO3 HCHO NO2 MCMv3.2
! B073 No rate in IUPAC or JPL
!----------------------------------------------------------------------
! B074 MeCOCH2OO+NO3 -> MeCO3 HCHO NO2 MCMv3.2
! B074 No rate in IUPAC or JPL
!----------------------------------------------------------------------
! B075 MeOO+MeCO3 -> HO2 HCHO MeOO IUPAC2002
! B075 JPL only gives total rates
!----------------------------------------------------------------------
! B076 MeOO+MeCO3 -> MeCO2H HCHO  IUPAC2002
! B076 JPL only gives total rates
!----------------------------------------------------------------------
! B077 MeOO+MeOO -> HO2 HO2 HCHOHCHO IUPAC2002 see also asad_bimol
! B077 IUPAC gives total k and k2 as f(T) resulting in complicated
! B077 expression evaluated in asad_bimol. JPL only gives a total rate
! B077 (9.5 and -390).
!----------------------------------------------------------------------
! B078 MeOO+MeOO -> MeOH HCHO  IUPAC2002 see also asad_bimol
! B078 IUPAC gives total k and k2 as f(T) resulting in complicated
! B078 expression evaluated in asad_bimol. JPL only gives a total rate
! B078 (9.5 and -390).
!----------------------------------------------------------------------
! B079 MeOO+NO -> HO2 HCHO NO2 IUPAC2005
! B079 IUPAC gives rates to CH3O. Assume 99.9% of the CH3O goes to HCHO
! B079 and HO2 and 0.1% goes to MeONO2. Where does this partition come
! B079 from?
!----------------------------------------------------------------------
! B080 MeOO+NO -> MeONO2   IUPAC2005
! B080 IUPAC gives rates to CH3O. Assume 99.9% of the CH3O goes to HCHO
! B080 and HO2 and 0.1% goes to MeONO2. Where does this partition come
! B080 from?
!----------------------------------------------------------------------
! B081 MeOO+NO3 -> HO2 HCHO NO2 MCMv3.2
! B081 No JPL value.
!----------------------------------------------------------------------
! B082 N+NO -> N2 O(3P)  JPL2011
! B082 No IUPAC value.
!----------------------------------------------------------------------
! B083 N+NO2 -> N2O O(3P)  JPL2011
! B083 No IUPAC value.
!----------------------------------------------------------------------
! B084 N+O2 -> NO O(3P)  JPL2011
! B084 IUPAC  values slightly different (5.1 and -198). Different
! B084 choices of results and fitting.
!----------------------------------------------------------------------
! B085 N2O5+H2O -> HONO2 HONO2  ????
! B085 Where is this from? Both JPL and IUPAC only give upper limits
! B085 (which are lower than this value)
!----------------------------------------------------------------------
! B086 NO+ISO2 -> ISON   Poschl00
! B086 -
!----------------------------------------------------------------------
! B087 NO+ISO2 -> NO2 MACR HCHOHO2 Poschl00
! B087 -
!----------------------------------------------------------------------
! B088 NO+MACRO2 -> MGLY HCHO HO2 Poschl00
! B088 Note to accommodate all the products split into two reactions and
! B088 half the rates.
!----------------------------------------------------------------------
! B089 NO+MACRO2 -> NO2 MeCO3 HACETCO Poschl00
! B089 Note to accommodate all the products split into two reactions and
! B089 half the rates.
!----------------------------------------------------------------------
! B090 NO+NO3 -> NO2 NO2  JPL2011
! B090 Discrepancy with IUPAC (1.8E-11 & beta =170) related to
! B090 differences in averaging.
!----------------------------------------------------------------------
! B091 NO+O3 -> NO2   JPL2011
! B091 Considerable discrepancy with IUPAC (1.4E-12 & beta=1310)
! B091 resulting in differences ~10% between the two. Use JPL as include
! B091 the extra (room temperature only though) of Stedman & Niki and
! B091 Bemand et al. Believe difference arises from fitting & averaging
! B091 vs averaging & fitting.
!----------------------------------------------------------------------
! B092 NO2+NO3 -> NO NO2 O2 JPL2011
! B092 This is a termolecular reaction. Why is it here? IUPAC don't
! B092 recognise the reaction and JPL note its existence is not
! B092 definitively proven. However they do use studies that infer the
! B092 rates from thermal decomposition of N2O5.
!----------------------------------------------------------------------
! B093 NO2+O3 -> NO3   JPL2011
! B093 Slight disagreement with IUPAC rates (1.4E-11 & -2470). Use JPL
! B093 by default.
!----------------------------------------------------------------------
! B094 NO3+Br -> BrO NO2  JPL2011
! B094 Both IUPAC & JPL recommend approach of Mellouki et al (1989)
!----------------------------------------------------------------------
! B095 NO3+C5H8 -> ISON   IUPAC2007
! B095 -
!----------------------------------------------------------------------
! B096 NO3+EtCHO -> HONO2 EtCO3  IUPAC2007
! B096 No JPL rates. No T dependence given buy IUPAC.
!----------------------------------------------------------------------
! B097 NO3+HCHO -> HONO2 HO2 CO IUPAC2007
! B097 No direct measurements of T dependence. Infer from T dependence
! B097 of MeCHO + NO3.  JPL values at 298K same as IUPAC.
!----------------------------------------------------------------------
! B098 NO3+MGLY -> MeCO3 CO HONO2 MCMv3.2
! B098 No rate in JPL/IUPAC. 2.4*IUPAC MeCHO + NO3.
!----------------------------------------------------------------------
! B099 NO3+Me2CO -> HONO2 MeCOCH2OO  IUPAC2007
! B099 Number is actually an upper limit. Use it as no other number
! B099 available (JPL don't provide one either)
!----------------------------------------------------------------------
! B100 NO3+MeCHO -> HONO2 MeCO3  IUPAC2007
! B100 -
!----------------------------------------------------------------------
! B101 O(1D)+CH4 -> HCHO H2  JPL2011
! B101 IUPAC values slightly lower (7.50E-12)
!----------------------------------------------------------------------
! B102 O(1D)+CH4 -> HCHO HO2 HO2 JPL2011
! B102 IUPAC values same
!----------------------------------------------------------------------
! B103 O(1D)+CH4 -> OH MeOO  JPL2011
! B103 IUPAC values slightly lower (1.05E-10)
!----------------------------------------------------------------------
! B104 O(1D)+CO2 -> O(3P) CO2  JPL2011
! B104 -
!----------------------------------------------------------------------
! B105 O(1D)+H2 -> OH H  IUPAC2008
! B105 JPL values identical
!----------------------------------------------------------------------
! B106 O(1D)+H2O -> OH OH  JPL2011
! B106 JPL value dominated by the work of Dunlea & Ravishankra. Slightly
! B106 lower than the IUPAC values (5% at room temperature)
!----------------------------------------------------------------------
! B107 O(1D)+HBr -> HBr O(3P)  JPL2011
! B107 Use Wine al for BR to O3P
!----------------------------------------------------------------------
! B108 O(1D)+HBr -> OH Br  JPL2011
! B108 Use Wine et al O3P BR and assume rest of HBR gos to OH + Br
!----------------------------------------------------------------------
! B109 O(1D)+HCl -> H ClO  JPL2011
! B109 Use Wine et al for branching ratios.
!----------------------------------------------------------------------
! B110 O(1D)+HCl -> O(3P) HCl  JPL2011
! B110 Use Wine et al for branching ratios.
!----------------------------------------------------------------------
! B111 O(1D)+HCl -> OH Cl  JPL2011
! B111 Use Wine et al for branching ratios.
!----------------------------------------------------------------------
! B112 O(1D)+N2 -> O(3P) N2  JPL2011
! B112 IUPAC identical
!----------------------------------------------------------------------
! B113 O(1D)+N2O -> N2 O2  JPL2011
! B113 Use JPL values (slightly higher than IUPAC 4.3e-11) as slightly
! B113 newer and include T dependence
!----------------------------------------------------------------------
! B114 O(1D)+N2O -> NO NO  JPL2011
! B114 Use JPL values (slightly lower than IUPAC 7.6e-11) as slightly
! B114 newer and include T dependence
!----------------------------------------------------------------------
! B115 O(1D)+O2 -> O(3P) O2  JPL2011
! B115 IUPAC values slightly different (3.2 and -67)
!----------------------------------------------------------------------
! B116 O(1D)+O3 -> O2 O(3P) O(3P) JPL2011
! B116 IUPAC values identical
!----------------------------------------------------------------------
! B117 O(1D)+O3 -> O2 O2  JPL2011
! B117 IUPAC values identical
!----------------------------------------------------------------------
! B118 O(3P)+BrO -> O2 Br  JPL2011
! B118 -
!----------------------------------------------------------------------
! B119 O(3P)+ClO -> Cl O2  JPL2011
! B119 Updated JPL 10-6
!----------------------------------------------------------------------
! B120 O(3P)+ClONO2 -> ClO NO3  JPL2011
! B120 Updated JPL 10-6
!----------------------------------------------------------------------
! B121 O(3P)+H2 -> OH H  ????
! B121 Not defined in IUPAC or JPL. Where does this come from?
!----------------------------------------------------------------------
! B122 O(3P)+H2O2 -> OH HO2  JPL2011
! B122 IUPAC values identical
!----------------------------------------------------------------------
! B123 O(3P)+HBr -> OH Br  JPL2011
! B123 -
!----------------------------------------------------------------------
! B124 O(3P)+HCHO -> OH CO HO2 JPL2011
! B124 Value not given in IUPAC (are we sure?)
!----------------------------------------------------------------------
! B125 O(3P)+HCl -> OH Cl  JPL2011
! B125 -
!----------------------------------------------------------------------
! B126 O(3P)+HO2 -> OH O2  IUPAC2001
! B126 JPL are the same to one significant figure.
!----------------------------------------------------------------------
! B127 O(3P)+HOCl -> OH ClO  JPL2011
! B127 -
!----------------------------------------------------------------------
! B128 O(3P)+NO2 -> NO O2  JPL2011
! B128 Note T dependence in IUPAC is slightly different (198 cf 210)
!----------------------------------------------------------------------
! B129 O(3P)+NO3 -> O2 NO2  IUPAC2009
! B129 Rate smaller in JPL2011 (1.1e-11) who prefer the results of
! B129 Graham & Johnston (1978) to Canosa-Mas et al. (1989). Though
! B129 IUPAC note that the results agree within uncertainties.
!----------------------------------------------------------------------
! B130 O(3P)+O3 -> O2 O2  JPL2011
! B130 IUPAC identical
!----------------------------------------------------------------------
! B131 O(3P)+OClO -> O2 ClO  JPL2011
! B131 -
!----------------------------------------------------------------------
! B132 O(3P)+OH -> O2 H  JPL2011
! B132 IUPAC slightly different (2.4 & -110). JPL include slightly more
! B132 results and newer results. However main difference probably
! B132 derives from fitting.
!----------------------------------------------------------------------
! B133 O3+C5H8 -> HO2 OH  IUPAC2007*
! B133 Total rates from IUPAC. Arbitrarily split into three reactions.
! B133 Fractional products given by Poschl et al (2000).
!----------------------------------------------------------------------
! B134 O3+C5H8 -> MACR HCHO MACRO2MeCO3 IUPAC2007*
! B134 Total rates from IUPAC. Arbitrarily split into three reactions.
! B134 Fractional products given by Poschl et al (2000).
!----------------------------------------------------------------------
! B135 O3+C5H8 -> MeOO HCOOH COH2O2 IUPAC2007*
! B135 Total rates from IUPAC. Arbitrarily split into three reactions.
! B135 Fractional products given by Poschl et al (2000).
!----------------------------------------------------------------------
! B136 O3+MACR -> MGLY HCOOH HO2CO IUPAC2007*
! B136 Complicated expression. Rate is average of IUPAC MACR + O3 and
! B136 MVK + O3. Split into further two reactions to allow all products
! B136 to be included. Hence rates are IUPAC/4. Fractional Products from
! B136 Poschl et al (2000).
!----------------------------------------------------------------------
! B137 O3+MACR -> MGLY HCOOH HO2CO IUPAC2007*
! B137 Complicated expression. Rate is average of IUPAC MACR + O3 and
! B137 MVK + O3. Split into further two reactions to allow all products
! B137 to be included. Hence rates are IUPAC/4. Fractional Products from
! B137 Poschl et al (2000).
!----------------------------------------------------------------------
! B138 O3+MACR -> OH MeCO3  IUPAC2007*
! B138 Complicated expression. Rate is average of IUPAC MACR + O3 and
! B138 MVK + O3. Split into further two reactions to allow all products
! B138 to be included. Hence rates are IUPAC/4. Fractional Products from
! B138 Poschl et al (2000).
!----------------------------------------------------------------------
! B139 O3+MACR -> OH MeCO3  IUPAC2007*
! B139 Complicated expression. Rate is average of IUPAC MACR + O3 and
! B139 MVK + O3. Split into further two reactions to allow all products
! B139 to be included. Hence rates are IUPAC/4. Fractional Products from
! B139 Poschl et al (2000).
!----------------------------------------------------------------------
! B140 OClO+NO -> NO2 ClO  JPL2011
! B140 Introduce T dependence of Bemand et al. New measurements by Li et
! B140 al in good agreement with these
!----------------------------------------------------------------------
! B141 OH+C2H6 -> H2O EtOO  IUPAC2007
! B141 JPL slightly different (7.22 & 1020). Difference arises from
! B141 fitting. Use IUPAC by default for organics.
!----------------------------------------------------------------------
! B142 OH+C3H8 -> i-PrOO H2O  IUPAC2007 see also asad_bimol
! B142 Total k rate (given) from IUPAC2007. Split between the two
! B142 channels is temperature dependent using measurements of Droege
! B142 and Tully (1986) and covered by asad_bimol. Note JPL values
! B142 slightly higher
!----------------------------------------------------------------------
! B143 OH+C3H8 -> n-PrOO H2O  IUPAC2007 see also asad_bimol
! B143 Total k rate (given) from IUPAC2007. Split between the two
! B143 channels is temperature dependent using measurements of Droege
! B143 and Tully (1986) and covered by asad_bimol. Note JPL values
! B143 slightly higher
!----------------------------------------------------------------------
! B144 OH+C5H8 -> ISO2   IUPAC2009
! B144 No JPL rate.
!----------------------------------------------------------------------
! B145 OH+CH4 -> H2O MeOO  JPL2011
! B145 Small difference with IUPAC results (1.85 and 1690)
!----------------------------------------------------------------------
! B146 OH+CO -> HO2   IUPAC2005 see also asad_bimol
! B146 Base rate modified by density dependence term in asad_bimol. JPL
! B146 treat as a termolecular reaction.
!----------------------------------------------------------------------
! B147 OH+ClO -> HCl O2  JPL2011
! B147 -
!----------------------------------------------------------------------
! B148 OH+ClO -> HO2 Cl  JPL2011
! B148 -
!----------------------------------------------------------------------
! B149 OH+ClONO2 -> HOCl NO3  JPL2011
! B149 -
!----------------------------------------------------------------------
! B150 OH+EtCHO -> H2O EtCO3  IUPAC2007
! B150 Latest value includes measurements of Le Crane et al (2005)>
!----------------------------------------------------------------------
! B151 OH+EtOOH -> H2O EtOO  MCMv3.2
! B151 IUPAC provides a limit
!----------------------------------------------------------------------
! B152 OH+EtOOH -> H2O MeCHO OH MCMv3.2
! B152 IUPAC provides a limit
!----------------------------------------------------------------------
! B153 OH+H2 -> H2O HO2  JPL2011
! B153 IUPAC different parameterisation (7.7 and 2100). Same at 298
!----------------------------------------------------------------------
! B154 OH+H2O2 -> H2O HO2  IUPAC2001
! B154 -
!----------------------------------------------------------------------
! B155 OH+HACET -> MGLY HO2  IUPAC2007
! B155 Includes temperature dependence of Dillon et al (2006)
!----------------------------------------------------------------------
! B156 OH+HBr -> H2O Br  JPL2011
! B156 -
!----------------------------------------------------------------------
! B157 OH+HCHO -> H2O HO2 CO IUPAC2007
! B157 -
!----------------------------------------------------------------------
! B158 OH+HCOOH -> HO2   IUPAC2007
! B158 -
!----------------------------------------------------------------------
! B159 OH+HCl -> H2O Cl  JPL2011
! B159 -
!----------------------------------------------------------------------
! B160 OH+HO2 -> H2O   JPL2011
! B160 IUPAC identical
!----------------------------------------------------------------------
! B161 OH+HO2NO2 -> H2O NO2  IUPAC2007
! B161 The preferred values are based on the recent and extensive
! B161 absolute rate study of Jimenez et al. (2004)
!----------------------------------------------------------------------
! B162 OH+HOCl -> ClO H2O  JPL2011
! B162 -
!----------------------------------------------------------------------
! B163 OH+HONO -> H2O NO2  IUPAC2004
! B163 IUPAC includes the more modern results of Burkholder. JPL values
! B163 slightly different (1.8 and 390)
!----------------------------------------------------------------------
! B164 OH+HONO2 -> H2O NO3  IUPAC2004 see also asad_bimol
! B164 Include rate with no density dependence. Density dependence
! B164 calculated using asad_bimol
!----------------------------------------------------------------------
! B165 OH+ISON -> HACET NALD  Poschl00
! B165 Use original MIM rate of Poschl et al (2000)
!----------------------------------------------------------------------
! B166 OH+ISOOH -> MACR OH  Poschl00
! B166 Use original MIM rate of Poschl et al (2000)
!----------------------------------------------------------------------
! B167 OH+MACR -> MACRO2   IUPAC2007
! B167 MACR rate is average of MACR+OH and MVK+OH. This is the IUPAC MVK
! B167 rate.
!----------------------------------------------------------------------
! B168 OH+MACR -> MACRO2   IUPAC2007
! B168 MACR rate is average of MACR+OH and MVK+OH. This is the IUPAC
! B168 MACR rate.
!----------------------------------------------------------------------
! B169 OH+MACROOH -> MACRO2   MCMv3.2
! B169  No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B170 OH+MGLY -> MeCO3 CO  IUPAC2008
! B170 -
!----------------------------------------------------------------------
! B171 OH+MPAN -> HACET NO2  IUPAC2006
! B171 -
!----------------------------------------------------------------------
! B172 OH+Me2CO -> H2O MeCOCH2OO  IUPAC2007
! B172 Two part reaction due to IUPAC recommendation of the form k1+k2
!----------------------------------------------------------------------
! B173 OH+Me2CO -> H2O MeCOCH2OO  IUPAC2007
! B173  Two part reaction due to IUPAC recommendation of the form k1+k2
!----------------------------------------------------------------------
! B174 OH+MeCHO -> H2O MeCO3  IUPAC2009
! B174 Use total cross section.
!----------------------------------------------------------------------
! B175 OH+MeCO2H -> MeOO   MCMv3.2
! B175 No data in IUPAC or JPL. Rate k is based on MCM v3.2. Can't find
! B175 the data used for previous versions of the code.
!----------------------------------------------------------------------
! B176 OH+MeCO3H -> MeCO3   MCMv3.2
! B176 No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B177 OH+MeCOCH2OOH -> H2O MeCOCH2OO  MCMv3.2
! B177  No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B178 OH+MeCOCH2OOH -> OH MGLY  MCMv3.2
! B178  No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B179 OH+MeOH -> HO2 HCHO  IUPAC2006
! B179 JPL rates same.
!----------------------------------------------------------------------
! B180 OH+MeONO2 -> HCHO NO2 H2O IUPAC2006
! B180 JPL rates slightly different 8.0E-13 and 1000.  Large
! B180 uncertainties on measurement.
!----------------------------------------------------------------------
! B181 OH+MeOOH -> H2O HCHO OH IUPAC2007
! B181 Use total rate from IUPAC * 0.4 (IUPAC BR to HCHO). Note JPL rate
! B181 is considerably lower. There are discrepancies between measured
! B181 values of this cross section which JPL and IUPAC reconcile
! B181 differently. We employ the IUPAC values.
!----------------------------------------------------------------------
! B182 OH+MeOOH -> H2O MeOO  IUPAC2007
! B182 Use total rate from IUPAC * 0.6 (IUPAC BR to MeOO). Note JPL rate
! B182 is considerably lower.
!----------------------------------------------------------------------
! B183 OH+NALD -> HCHO CO NO2 IUPAC2009
! B183 As in Poschl et al use MeCHO + OH rates (but updated).
!----------------------------------------------------------------------
! B184 OH+NO3 -> HO2 NO2  JPL2011
! B184 Note small discrepancy with IUPAC values (JPL rates 10% higher).
! B184 Same measurements used
!----------------------------------------------------------------------
! B185 OH+O3 -> HO2 O2  JPL2011
! B185 Same values in IUPAC.
!----------------------------------------------------------------------
! B186 OH+OClO -> HOCl O2  JPL2011
! B186 Updated based on the recommended value reported by Gierczak et
! B186 al.
!----------------------------------------------------------------------
! B187 OH+OH -> H2O O(3P)  IUPAC2001
! B187 No T dependence given by JPL
!----------------------------------------------------------------------
! B188 OH+PAN -> HCHO NO2 H2O MCMv3.2
! B188 No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B189 OH+PPAN -> MeCHO NO2 H2O MCMv3.2
! B189 No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B190 OH+i-PrOOH -> Me2CO OH  MCMv3.2
! B190 No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B191 OH+i-PrOOH -> i-PrOO H2O  MCMv3.2
! B191 No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B192 OH+n-PrOOH -> EtCHO H2O OH MCMv3.2
! B192 Rate k is based on MCM v3.2. IUPAC only provides total cross
! B192 sections.
!----------------------------------------------------------------------
! B193 OH+n-PrOOH -> n-PrOO H2O  MCMv3.2
! B193 Rate k is based on MCM v3.2. IUPAC only provides total cross
! B193 sections.
!----------------------------------------------------------------------
! B194 i-PrOO+NO -> Me2CO HO2 NO2 IUPAC2005
! B194 -
!----------------------------------------------------------------------
! B195 i-PrOO+NO3 -> Me2CO HO2 NO2 MCMv3.2
! B195 No data in IUPAC or JPL. Take rates from MCM v3.2.
!----------------------------------------------------------------------
! B196 n-PrOO+NO -> EtCHO HO2 NO2 IUPAC2005
! B196 The recommendation accepts the Arrhenius expression of Eberhard
! B196 and Howard (1996).  Neglect n-propyl nitrate formation.
!----------------------------------------------------------------------
! B197 n-PrOO+NO3 -> EtCHO HO2 NO2 MCM3.2
! B197 No data in IUPAC or JPL. Take rates from MCM v3.2.
!----------------------------------------------------------------------
! Notes on RO2+RO2 permutation reactions
! B278 MeOO+RO2 -> HO2+HCHO
! Equivalent to IUPAC2002 MeOO+MeOO->2HO2+2HCHO reaction, with k rate doubled
! To give same consumption of MeOO when in MeOO+RO2 form
! See asad_bimol for temperature dependent branching ratio
!----------------------------------------------------------------------
! B279 MeOO+RO2 -> 0.5*HO2+0.5*HCHO
! Equivalent to IUPAC2002 MeOO+MeOO->2MeOH+2HCHO reaction
! with k rate doubled and products halved to give same net
! consumption of MeOO when in MeOO+RO2 form.
! See asad_bimol for temperature dependent branching ratio
!----------------------------------------------------------------------
! B280 EtOO+RO2->HO2+EtCHO
! From MCM 3.3.1
!----------------------------------------------------------------------
! B281 EtOO+RO2->0.5*EtOH+0.5*MeCHO
! From MCM3.3.1, pathways 2&3 merged into single reaction, EtOH=>MeCHO
!----------------------------------------------------------------------
! B282 n-PrOO+RO2->HO2+EtCHO
! From MCM 3.3.1
!----------------------------------------------------------------------
! B283 n-PrOO+RO2->0.5*n-PrOH+0.5*EtCHO
! From MCM3.3.1, pathways 2&3 merged into single reaction, n-PrOH=>EtCHO
!----------------------------------------------------------------------
! B284 i-PrOO+RO2->HO2+Me2CO
! From MCM 3.3.1
!----------------------------------------------------------------------
! B285 i-PrOO+RO2->0.5*i-PrOH+0.5*Me2CO
! From MCM3.3.1, pathways 2&3 merged into single reaction, i-PrOH=>Me2CO
!----------------------------------------------------------------------
! B286 MeCO3+RO2->0.3*MeCO2H+0.7*MeOO
! From MCM3.3.1, 2 pathways merged into single reaction
!----------------------------------------------------------------------
! B287 EtCO3+RO2->0.3*EtCO2H+0.7*EtOO
! - NOTE EtCO2H IS *NOT* CURRENTLY A SPECIES CONSIDERED BY ASAD
! From MCM3.3.1, 2 pathways merged into single reaction
!----------------------------------------------------------------------
! B288 MeCOCH2OO+RO2->HCHO+MeCO3
! From MCM3.3.1
!----------------------------------------------------------------------
! B289 MeCOCH2OO+RO2->0.5*HACET+0.5*MGLY
! From MCM3.3.1, 2 pathways merged into single reaction
!----------------------------------------------------------------------
! B290 ISO2+RO2->MACR+0.5*HCHO+0.5HO2
! Poschl00, reformulated to RO2+RO2 style by doubling k rate and
! halving products
! Merged reaction rate calculated as:
! k0 = 2*(k_iso2+iso2 * k_MeOO+MeOO)^0.5
! Where:
! k_iso2+iso2 = 2E-12; k_MeOO+MeOO = 3.5E-13
! See MCM website for details of methodology
! mcm.leeds.ac.uk/MCM/categories/saunders-2003-4_6_5-gen-master.htt?rxnId=57
!----------------------------------------------------------------------
! B291 MACRO2+RO2->HACET+MGLY+0.5HCHO+0.5*CO
! Poschl00, reformulated to RO2+RO2 style by doubling k rate and
! halving products, following same methodology as ISO2
!----------------------------------------------------------------------
! B292 MACRO2+RO2->HO2
! Poschl00, reformulated to RO2+RO2 style by doubling k rate and
! halving products, following same methodology as ISO2
!----------------------------------------------------------------------
! B313-B697 Common Representative Intermediate V2-R5 reactions
! Please see CRI website for details
! http://mcm.leeds.ac.uk/CRI/
!----------------------------------------------------------------------
! CRI-Strat 2 (CS2) includes updates to some existing reactions in
! the range B313-B697 and new reactions in the range B698-B743
! Please see CRI website for details
! http://cri.york.ac.uk
!----------------------------------------------------------------------


! Dry deposition coefficients. Only used for off-line dry-deposition scheme.
! Item number, species name, dry dep velocities x30, chemistry scheme,
! qualifier, disqualifier, version

! Split to avoid warning "Too many continuation lines"
! 1
! R and T are at older revision than S and ST. Make consistent
depvel_defs_master(1:29) = [                                                   &
depvel_t(1,'O3        ',& ! (Ganzeveld & Lelieveld (1995) note 1)
                          ! (modified to be the same as Guang version)
[0.05,  0.05,  0.05,  0.05,  0.05,  0.05,&   ! DD:  1.1
  0.85,  0.30,  0.65,  0.65,  0.25,  0.45,&   ! DD:  1.2
  0.65,  0.25,  0.45,  0.65,  0.25,  0.45,&   ! DD:  1.3
  0.18,  0.18,  0.18,  0.18,  0.18,  0.18,&   ! DD:  1.4
  0.05,  0.05,  0.05,  0.05,  0.05,  0.05],& ! DD:  1.5
  ti+st+cs,0,0,107),                                                           &
! 2*
depvel_t(1,'O3        ',                                                       &
[0.05,  0.05,  0.05,  0.05,  0.05,  0.05,                                      &
  1.00,  0.11,  0.56,  0.26,  0.11,  0.19,                                     &
  1.00,  0.37,  0.69,  0.59,  0.46,  0.53,                                     &
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,                                     &
  0.05,  0.05,  0.05,  0.05,  0.05,  0.05],                                    &
  t+r,0,0,107),                                                                &
! 3**
! O3 (Ganzeveld& Lelieveld (1995) - note 1)
depvel_t(1,'O3        ',                                                       &
[0.07,  0.07,  0.07,  0.07,  0.07,  0.07,                                      &
  1.00,  0.11,  0.56,  0.26,  0.11,  0.19,                                     &
  1.00,  0.37,  0.69,  0.59,  0.46,  0.53,                                     &
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,                                     &
  0.07,  0.07,  0.07,  0.07,  0.07,  0.07],                                    &
  s,0,0,107),                                                                  &
! 4
! No DD of NO in R scheme
depvel_t(4,'NO        ',& ! (inferred from NO2 - see Giannakopoulos (1998))
[0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:  2.1
  0.14,  0.01,  0.07,  0.01,  0.01,  0.01,&   ! DD:  2.2
  0.10,  0.01,  0.06,  0.01,  0.01,  0.01,&   ! DD:  2.3
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,&   ! DD:  2.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:  2.5
  ti+t+s+st+cs,0,0,107),                                                       &
! 5
! No DD of NO3 in R scheme
depvel_t(5,'NO3     ',& ! (as NO2)
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:  3.1
  0.83,  0.04,  0.44,  0.06,  0.04,  0.05,&   ! DD:  3.2
  0.63,  0.06,  0.35,  0.07,  0.06,  0.07,&   ! DD:  3.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:  3.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:  3.5
  ti+t+s+st+cs,0,0,107),                                                       &
! 6
depvel_t(6,'NO2       ',                                                       &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:  4.1
  0.83,  0.04,  0.44,  0.06,  0.04,  0.05,&   ! DD:  4.2
  0.63,  0.06,  0.35,  0.07,  0.06,  0.07,&   ! DD:  4.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:  4.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:  4.5
  ti+t+s+st+r+cs,0,0,107),                                                     &
! 7
! No DD of N2O5 in R scheme
depvel_t(7,'N2O5      ',& ! (as HNO3)
[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD:  5.1
  4.00,  3.00,  3.50,  2.00,  1.00,  1.50,&   ! DD:  5.2
  2.50,  1.50,  2.00,  1.00,  1.00,  1.00,&   ! DD:  5.3
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD:  5.4
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50],& ! DD:  5.5
  ti+t+s+st+cs,0,0,107),                                                       &
! 8
! No DD of HO2NO2 in R scheme
depvel_t(8,'HO2NO2    ',& ! (as HNO3)
[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD:  6.1
  3.20,  1.40,  2.30,  2.00,  1.00,  1.50,&   ! DD:  6.2
  1.80,  0.80,  1.30,  1.00,  1.00,  1.00,&   ! DD:  6.3
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD:  6.4
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50],& ! DD:  6.5
  ti+st+cs,0,0,107),                                                           &
! 9*
depvel_t(8,'HO2NO2    ',& ! (as HNO3)
[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                                      &
  4.00,  3.00,  3.50,  2.00,  1.00,  1.50,                                     &
  2.50,  1.50,  2.00,  1.00,  1.00,  1.00,                                     &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                                     &
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50],                                    &
  t+s,0,0,107),                                                                &
! 10
! T and R are outdated versus S, ST
depvel_t(10,'HONO2     ',& ! (Zhang et al., 2003)
[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD:  7.1
  3.20,  1.40,  2.30,  2.00,  1.00,  1.50,&   ! DD:  7.2
  1.80,  0.80,  1.30,  1.00,  1.00,  1.00,&   ! DD:  7.3
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD:  7.4
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50],& ! DD:  7.5
  ti+st+cs,0,0,107),                                                           &
! 11*
depvel_t(10,'HONO2     ',                                                      &
[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                                      &
  4.00,  3.00,  3.50,  2.00,  1.00,  1.50,                                     &
  2.50,  1.50,  2.00,  1.00,  1.00,  1.00,                                     &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                                     &
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50],                                    &
  s+t+r,0,0,107),                                                              &
! 12
depvel_t(12,'H2O2      ',                                                      &
[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD:  8.1
  1.25,  0.16,  0.71,  0.28,  0.12,  0.20,&   ! DD:  8.2
  1.25,  0.53,  0.89,  0.83,  0.78,  0.81,&   ! DD:  8.3
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,&   ! DD:  8.4
  0.32,  0.32,  0.32,  0.32,  0.32,  0.32],& ! DD:  8.5
  ti+t+s+st+ol+r+cs,0,0,107),                                                  &
! 13
depvel_t(13,'CH4       ',                                                      &
[0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:  9.1
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:  9.2
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:  9.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:  9.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:  9.5
  ti+t+st+r,0,0,107),                                                          &
! 14
depvel_t(14,'CO        ',& ! (see Giannokopoulos (1998))
[0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 10.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 10.2
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 10.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 10.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 10.5
  ti+t+s+st+r+cs,0,0,107),                                                     &
! 15
! No DD for HCHO in R scheme
depvel_t(15,'HCHO      ',& ! (Zhang et al., 2003)
[1.00,  0.60,  0.80,  1.00,  0.60,  0.80,&   ! DD: 11.1
  1.60,  0.40,  1.00,  0.01,  0.01,  0.01,&   ! DD: 11.2
  0.71,  0.20,  0.45,  0.03,  0.03,  0.03,&   ! DD: 11.3
  0.40,  0.40,  0.40,  0.00,  0.00,  0.00,&   ! DD: 11.4
  0.30,  0.30,  0.30,  0.30,  0.30,  0.30],& ! DD: 11.5
  ti+st+cs,0,0,107),                                                           &
! 16*
depvel_t(15,'HCHO      ',                                                      &
[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                                      &
  1.00,  0.01,  0.50,  0.01,  0.01,  0.01,                                     &
  0.71,  0.03,  0.37,  0.03,  0.03,  0.03,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                                     &
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04],                                    &
  s+t,0,0,107),                                                                &
! 17
depvel_t(17,'MeOOH     ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 12.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 12.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 12.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 12.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 12.5
  ti+t+s+st+r+cs,0,0,107),                                                     &
! 18
depvel_t(18,'HCl       ',& ! (same as HBr)
[0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 13.1
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD: 13.2
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD: 13.3
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 13.4
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20],& ! DD: 13.5
  s+st+cs,0,0,107),                                                            &
! 19
depvel_t(19,'HOCl      ',& ! (same as HOBr)
[0.20,  0.20,  0.20,  0.20,  0.20,  0.20,&   ! DD: 14.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 14.2
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 14.3
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,&   ! DD: 14.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10],& ! DD: 14.5
  s+st+cs,0,0,107),                                                            &
! 20
depvel_t(20,'HBr       ',& ! (note 3)
[0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 15.1
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD: 15.2
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD: 15.3
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 15.4
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20],& ! DD: 15.5
  s+st+cs,0,0,107),                                                            &
! 21
depvel_t(21,'HOBr      ',& ! (note 3)
[0.20,  0.20,  0.20,  0.20,  0.20,  0.20,&   ! DD: 16.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 16.2
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 16.3
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,&   ! DD: 16.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10],& ! DD: 16.5
  s+st+cs,0,0,107),                                                            &
! 22
depvel_t(22,'HONO      ',                                                      &
[0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 17.1
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 17.2
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 17.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 17.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 17.5
  ti+t+st+cs,0,0,107),                                                         &
! 23
depvel_t(23,'EtOOH     ',& !(as MeOOH)
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 18.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 18.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 18.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 18.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 18.5
  ti+t+st+r+cs,0,0,107),                                                       &
! 24
! No DD in R scheme for MeCHO
depvel_t(24,'MeCHO     ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 19.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 19.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 19.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 19.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 19.5
  ti+t+st+cs,0,0,107),                                                         &
! 25
! Make consistent
depvel_t(25,'PAN       ',& ! (Zhang et al., 2003)
[0.01,  0.01,  0.01,  0.01,  0.01,  0.01,&   ! DD: 20.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,&   ! DD: 20.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,&   ! DD: 20.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 20.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00],& ! DD: 20.5
  ti+st+cs,0,0,107),                                                           &
! 26*
depvel_t(25,'PAN       ',                                                      &
[0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                                      &
  0.53,  0.04,  0.29,  0.06,  0.05,  0.06,                                     &
  0.42,  0.06,  0.24,  0.07,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00],                                    &
  t+r,0,0,107),                                                                &
! 27
depvel_t(27,'n-PrOOH   ',& ! (as MeOOH)
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 21.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 21.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 21.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 21.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 21.5
  ti+t+st,0,0,107),                                                            &
! 28
depvel_t(28,'i-PrOOH   ',& ! (as MeOOH)
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 22.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 22.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 22.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 22.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 22.5
  ti+t+st+r+cs,0,0,107),                                                       &
! 29
depvel_t(29,'EtCHO     ',& ! (as MeCHO)
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 23.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 23.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 23.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 23.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 23.5
  ti+t+st+cs,0,0,107) ]

! Split to avoid warning "Too many continuation lines"
depvel_defs_master(30:59) = [                                                  &
! 30
depvel_t(30,'MeCOCH2OOH',& !(as MeCO3H - where do the desert data come from!?)
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 24.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 24.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 24.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 24.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 24.5
  ti+t+st,0,0,107),                                                            &
! 31
depvel_t(31,'PPAN      ',& ! (as PAN)
[0.01,  0.01,  0.01,  0.01,  0.01,  0.01,&   ! DD: 25.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,&   ! DD: 25.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,&   ! DD: 25.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 25.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00],& ! DD: 25.5
  ti+st+cs,0,0,107),                                                           &
! 32*
depvel_t(31,'PPAN      ',                                                      &
[0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                                      &
  0.53,  0.04,  0.29,  0.06,  0.05,  0.06,                                     &
  0.42,  0.06,  0.24,  0.07,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00],                                    &
  t,0,0,107),                                                                  &
! 33
! Inconsistent DD for ISOOH
depvel_t(33,'ISOOH     ',& ! (as MeCO3H) (MIM)
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 26.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 26.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 26.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 26.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 26.5
  ti+st,0,0,107),                                                              &
! 34*
depvel_t(33,'ISOOH     ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,& !ISOOH
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],                                    &
  r,0,0,107),                                                                  &
! 35
depvel_t(35,'ISON      ',& ! (as MeCO2H) (MIM) - note 2
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD: 27.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD: 27.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD: 27.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD: 27.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD: 27.5
  ti+st,0,0,107),                                                              &
! 36
depvel_t(36,'MACR      ',& ! (as MeCHO) (MIM)
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 28.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 28.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 28.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 28.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 28.5
  ti+st,0,0,107),                                                              &
! 37
depvel_t(37,'MACROOH   ',& ! (as MeCO3H) (MIM)
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 29.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 29.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 29.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 29.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 29.5
  ti+st,0,0,107),                                                              &
! 38
depvel_t(38,'MPAN      ',& ! (as PAN) (MIM)
[0.01,  0.01,  0.01,  0.01,  0.01,  0.01,&   ! DD: 30.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,&   ! DD: 30.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,&   ! DD: 30.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 30.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00],& ! DD: 30.5
  ti+st+cs,0,0,107),                                                           &
! 39
depvel_t(39,'HACET     ',& ! (as MeCO3H) (MIM)
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 31.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 31.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 31.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 31.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 31.5
  ti+st,0,0,107),                                                              &
! 40
depvel_t(40,'MGLY      ',& ! (as MeCHO) (MIM)
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 32.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 32.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 32.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 32.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 32.5
  ti+st,0,0,107),                                                              &
! 41
depvel_t(41,'NALD      ',& ! (as MeCHO) (MIM)
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 33.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 33.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 33.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 33.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 33.5
  ti+st,0,0,107),                                                              &
! 42
depvel_t(42,'HCOOH     ',& ! (from Sander & Crutzen (1996) - note 3) (MIM)
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD: 34.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD: 34.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD: 34.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD: 34.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD: 34.5
  ti+st+cs,0,0,107),                                                           &
! 43
depvel_t(43,'MeCO3H    ',& ! (from Giannakopoulos (1998)) (MIM)
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 35.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 35.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 35.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 35.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 35.5
  ti+st+cs,0,0,107),                                                           &
! 44
depvel_t(44,'MeCO2H    ',& ! (as HCOOH) (MIM)
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD: 36.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD: 36.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD: 36.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD: 36.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD: 36.5
  ti+st+cs,0,0,107),                                                           &
! 45
depvel_t(45,'H2        ',                                                      &
[0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 37.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 37.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 37.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 37.1
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 37.1
  r,0,0,107),                                                                  &
! 46
depvel_t(46,'MeOH      ',& ! (as MeOOH)
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 38.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 38.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 38.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 38.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 38.5
  ti+st+cs,0,0,107),                                                           &
! 47
depvel_t(47,'SO2       ',                                                      &
[0.80,  0.80,  0.80,  0.80,  0.80,  0.80,&   ! DD: 39.1
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,&   ! DD: 39.2
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,&   ! DD: 39.3
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,&   ! DD: 39.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10],& ! DD: 39.5
  ti+s+st+ol+r+cs,a,0,107),                                                    &
! 48
depvel_t(48,'H2SO4     ',                                                      &
[0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 40.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 40.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 40.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 40.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50],& ! DD: 40.1
  ti+ol+s+st+r+cs,a,0,107),                                                    &
! 49
depvel_t(49,'MSA       ',                                                      &
[0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 41.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 41.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 41.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,&   ! DD: 41.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50],& ! DD: 41.1
  ti,a,0,107),                                                                 &
! 50
! no DD of DMSO in R scheme
depvel_t(50,'DMSO      ',                                                      &
[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&   ! DD: 42.1
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 42.2
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 42.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 42.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 42.5
  ti+st+ol+cs,a,0,107),                                                        &
! 51
depvel_t(51,'NH3       ',                                                      &
[0.80,  0.80,  0.80,  0.80,  0.80,  0.80,&   ! DD: 43.1
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,&   ! DD: 43.2
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,&   ! DD: 43.3
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,&   ! DD: 43.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10],& ! DD: 43.5
  ti+st+r+cs,a,0,107),                                                         &
! 52
! No DD of Monoterp in R scheme
depvel_t(52,'Monoterp  ',                                                      &
[0.20,  0.20,  0.20,  0.20,  0.20,  0.20,&   ! DD: 44.1
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 44.2
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 44.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 44.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10],& ! DD: 44.5
  ti+st+ol,a,0,107),                                                           &
! 53
! No DD of Sec_Org in R scheme
depvel_t(53,'Sec_Org   ',                                                      &
[0.20,  0.20,  0.20,  0.20,  0.20,  0.20,&   ! DD: 45.1
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 45.2
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 45.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 45.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10],& ! DD: 45.5
  ti+st+ol+cs,a,0,107),                                                        &
! 54
depvel_t(54,'O3S       ',                                                      &
[0.07,  0.07,  0.07,  0.07,  0.07,  0.07,&   ! DD: 46.1
  1.00,  0.11,  0.56,  0.26,  0.11,  0.19,&   ! DD: 46.2
  1.00,  0.37,  0.69,  0.59,  0.46,  0.53,&   ! DD: 46.3
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,&   ! DD: 46.4
  0.07,  0.07,  0.07,  0.07,  0.07,  0.07],& ! DD: 46.1
  t+r,0,0,107),                                                                &
! 55
depvel_t(55,'MVKOOH    ',                                                      &
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 47.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 47.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 47.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 47.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 47.1
  r,0,0,107),                                                                  &
! 56
depvel_t(56,'ORGNIT    ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 48.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 48.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 48.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 48.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 48.1
  r,0,0,107),                                                                  &
! 57
depvel_t(57,'s-BuOOH   ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 49.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 49.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 49.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 49.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 49.1
  r,0,0,107),                                                                  &
! CRI Dry-dep species from here!
!  58 APINENE  Using monoterp DD rates
depvel_t(58,'APINENE   ',                                                      &
[0.20,  0.20,  0.20,  0.20,  0.20,  0.20,&   ! DD: 50.1
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 50.2
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 50.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 50.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10],& ! DD: 50.5
  cs,0,0,107),                                                                 &
!  59 BPINENE  Using monoterp DD rates
depvel_t(59,'BPINENE   ',                                                      &
[0.20,  0.20,  0.20,  0.20,  0.20,  0.20,&   ! DD: 51.1
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 51.2
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 51.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 51.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10],& ! DD: 51.5
  cs,0,0,107)                                                                  &
]

depvel_defs_master(60:89) = [                                                  &
!  60 EtOH  Using MeOH DD rates
depvel_t(60,'EtOH      ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 52.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 52.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 52.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 52.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 52.5
  cs,0,0,107),                                                                 &
!  61 i-PrOH  Using MeOH DD rates
depvel_t(61,'i-PrOH    ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 53.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 53.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 53.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 53.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 53.5
  cs,0,0,107),                                                                 &
!  62 n-PrOH  Using MeOH DD rates
depvel_t(62,'n-PrOH    ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 54.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 54.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 54.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 54.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 54.5
  cs,0,0,107),                                                                 &
!  63 HOCH2CHO  Using MeCHO DD rates
depvel_t(63,'HOCH2CHO  ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 55.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 55.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 55.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 55.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 55.5
  cs,0,0,107),                                                                 &
!  64 HOC2H4OOH  Using EtOOH DD rates
depvel_t(64,'HOC2H4OOH ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 56.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 56.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 56.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 56.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 56.5
  cs,0,0,107),                                                                 &
!  65 EtCO3H  Using MeCO3H DD rates
depvel_t(65,'EtCO3H    ',                                                      &
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 57.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 57.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 57.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 57.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 57.5
  cs,0,0,107),                                                                 &
!  66 HOCH2CO3H  Using MeCO3H DD rates
depvel_t(66,'HOCH2CO3H ',                                                      &
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 58.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 58.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 58.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 58.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 58.5
  cs,0,0,107),                                                                 &
!  67 NOA  Using NALD DD rates
depvel_t(67,'NOA       ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 59.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 59.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 59.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 59.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 59.5
  cs,0,0,107),                                                                 &
!  68 HOC2H4NO3  Using ISON DD rates
depvel_t(68,'HOC2H4NO3 ',                                                      &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD: 60.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD: 60.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD: 60.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD: 60.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD: 60.5
  cs,0,0,107),                                                                 &
!  69 PHAN  Using PAN DD rates
depvel_t(69,'PHAN      ',                                                      &
[0.01,  0.01,  0.01,  0.01,  0.01,  0.01,&   ! DD: 61.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,&   ! DD: 61.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,&   ! DD: 61.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD: 61.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00],& ! DD: 61.5
  cs,0,0,107),                                                                 &
!  70 CARB14  Using MeCHO DD rates
depvel_t(70,'CARB14    ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 62.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 62.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 62.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 62.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 62.5
  cs,0,0,107),                                                                 &
!  71 CARB17  Using MeCHO DD rates
depvel_t(71,'CARB17    ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 63.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 63.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 63.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 63.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 60.5
  cs,0,0,107),                                                                 &
!  72 CARB11A  Using MeCHO DD rates
depvel_t(72,'CARB11A   ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 64.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 64.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 64.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 64.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 64.5
  cs,0,0,107),                                                                 &
!  73 CARB7  Using HACET DD rates
depvel_t(73,'CARB7     ',                                                      &
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 65.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 65.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 65.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 65.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 65.5
  cs,0,0,107),                                                                 &
!  74 CARB10  Using HACET DD rates
depvel_t(74,'CARB10    ',                                                      &
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 66.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 66.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 66.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 66.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 66.5
  cs,0,0,107),                                                                 &
!  75 CARB13  Using HACET DD rates
depvel_t(75,'CARB13    ',                                                      &
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 67.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 67.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 67.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 67.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 67.5
  cs,0,0,107),                                                                 &
!  76 CARB16  Using HACET DD rates
depvel_t(76,'CARB16    ',                                                      &
[0.36,  0.36,  0.36,  0.19,  0.19,  0.19,&   ! DD: 68.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,&   ! DD: 68.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,&   ! DD: 68.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 68.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02],& ! DD: 68.5
  cs,0,0,107),                                                                 &
!  77 UCARB10  Using MeCHO DD rates
depvel_t(77,'UCARB10   ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 69.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 69.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 69.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 69.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 69.5
  cs,0,0,107),                                                                 &
!  78 CARB3  Using MGLY DD rates
depvel_t(78,'CARB3     ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 70.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 70.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 70.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 70.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 70.5
  cs,0,0,107),                                                                 &
!  79 CARB6  Using MGLY DD rates
depvel_t(79,'CARB6     ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 71.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 71.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 71.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 71.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 71.5
  cs,0,0,107),                                                                 &
!  80 CARB9  Using MGLY DD rates
depvel_t(80,'CARB9     ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 72.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 72.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 72.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 72.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 72.5
  cs,0,0,107),                                                                 &
!  81 CARB12  Using MGLY DD rates
depvel_t(81,'CARB12    ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 73.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 73.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 73.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 73.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 73.5
  cs,0,0,107),                                                                 &
!  82 CARB15  Using MGLY DD rates
depvel_t(82,'CARB15    ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 74.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 74.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 74.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 74.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 74.5
  cs,0,0,107),                                                                 &
!  83 UCARB12  Using MeCHO DD rates
depvel_t(83,'UCARB12   ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 75.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 75.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 75.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 75.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 75.5
  cs,0,0,107),                                                                 &
!  84 NUCARB12  Using MeCHO DD rates
depvel_t(84,'NUCARB12  ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 76.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 76.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 76.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 76.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 76.5
  cs,0,0,107),                                                                 &
!  85 UDCARB8  Using MeCHO DD rates
depvel_t(85,'UDCARB8   ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 77.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 77.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 77.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 77.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 77.5
  cs,0,0,107),                                                                 &
!  86 UDCARB11  Using MeCHO DD rates
depvel_t(86,'UDCARB11  ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 78.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 78.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 78.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 78.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 78.5
  cs,0,0,107),                                                                 &
!  87 UDCARB14  Using MeCHO DD rates
depvel_t(87,'UDCARB14  ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 79.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 79.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 79.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 79.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 79.5
  cs,0,0,107),                                                                 &
!  88 TNCARB26  Using MeCHO DD rates
depvel_t(88,'TNCARB26  ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 80.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 80.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 80.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 80.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 80.5
  cs,0,0,107),                                                                 &
!  89 TNCARB10  Using MeCHO DD rates
depvel_t(89,'TNCARB10  ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 81.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 81.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 81.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 81.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 81.5
  cs,0,0,107)                                                                  &
]

! Dry deposition continued
depvel_defs_master(90:122) = [                                                 &
!  90 RA13NO3  Using ISON DD rates
depvel_t(90,'RA13NO3   ',                                                      &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD:128.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD:128.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD:128.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD:128.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD:128.5
 cs,0,0,119),                                                                  &
!  91 RA16NO3  Using ISON DD rates
depvel_t(91,'RA16NO3   ',                                                      &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD:128.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD:128.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD:128.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD:128.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD:128.5
 cs,0,0,119),                                                                  &
!  92 RA19NO3  Using ISON DD rates
depvel_t(92,'RA19NO3   ',                                                      &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD:128.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD:128.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD:128.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD:128.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD:128.5
 cs,0,0,119),                                                                  &
!  93 RTX24NO3  Using NALD DD rates
depvel_t(93,'RTX24NO3  ',                                                      &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD: 82.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD: 82.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD: 82.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD: 82.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD: 82.5
  cs,0,0,107),                                                                 &
!  94 RN10OOH  Using EtOOH DD rates
depvel_t(94,'RN10OOH   ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 83.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 83.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 83.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 83.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 83.5
  cs,0,0,107),                                                                 &
!  95 RN13OOH  Using EtOOH DD rates
depvel_t(95,'RN13OOH   ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 84.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 84.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 84.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 84.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 84.5
  cs,0,0,107),                                                                 &
!  96 RN16OOH  Using EtOOH DD rates
depvel_t(96,'RN16OOH   ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 85.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 85.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 85.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 85.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 85.5
  cs,0,0,107),                                                                 &
!  97 RN19OOH  Using EtOOH DD rates
depvel_t(97,'RN19OOH   ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 86.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 86.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 86.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 86.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 86.5
  cs,0,0,107),                                                                 &
!  98 RN8OOH  Using EtOOH DD rates
depvel_t(98,'RN8OOH    ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 87.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 87.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 87.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 87.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 87.5
  cs,0,0,107),                                                                 &
!  99 RN11OOH  Using EtOOH DD rates
depvel_t(99,'RN11OOH   ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 88.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 88.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 88.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 88.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 88.5
  cs,0,0,107),                                                                 &
! 100 RN14OOH  Using EtOOH DD rates
depvel_t(100,'RN14OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 89.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 89.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 89.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 89.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 89.5
  cs,0,0,107),                                                                 &
! 101 RN17OOH  Using EtOOH DD rates
depvel_t(101,'RN17OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 90.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 90.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 90.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 90.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 90.5
  cs,0,0,107),                                                                 &
! 102 RU14OOH  Using EtOOH DD rates
depvel_t(102,'RU14OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 91.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 91.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 91.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 91.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 91.5
  cs,0,0,107),                                                                 &
!  103 RU12OOH  Using EtOOH DD rates
depvel_t(103,'RU12OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 92.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 92.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 92.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 92.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 92.5
  cs,0,0,107),                                                                 &
! 104 RU10OOH  Using EtOOH DD rates
depvel_t(104,'RU10OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 93.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 93.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 93.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 93.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 93.5
  cs,0,0,107),                                                                 &
! 105 NRU14OOH  Using EtOOH DD rates
depvel_t(105,'NRU14OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 94.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 94.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 94.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 94.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 94.5
  cs,0,0,107),                                                                 &
!  106 NRU12OOH  Using EtOOH DD rates
depvel_t(106,'NRU12OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 95.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 95.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 95.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 95.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 95.5
  cs,0,0,107),                                                                 &
!  107 RN9OOH  Using EtOOH DD rates
depvel_t(107,'RN9OOH    ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 96.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 96.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 96.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 96.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 96.5
  cs,0,0,107),                                                                 &
!  108 RN12OOH  Using EtOOH DD rates
depvel_t(108,'RN12OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 97.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 97.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 97.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 97.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 97.5
  cs,0,0,107),                                                                 &
!  109 RN15OOH  Using EtOOH DD rates
depvel_t(109,'RN15OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 98.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 98.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 98.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 98.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 98.5
  cs,0,0,107),                                                                 &
!  110 RN18OOH  Using EtOOH DD rates
depvel_t(110,'RN18OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD: 99.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD: 99.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD: 99.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD: 99.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD: 99.5
  cs,0,0,107),                                                                 &
!  111 NRN6OOH  Using EtOOH DD rates
depvel_t(111,'NRN6OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:100.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:100.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:100.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:100.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:100.5
  cs,0,0,107),                                                                 &
!  112 NRN9OOH  Using EtOOH DD rates
depvel_t(112,'NRN9OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:101.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:101.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:101.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:101.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:101.5
  cs,0,0,107),                                                                 &
!  113 NRN12OOH  Using EtOOH DD rates
depvel_t(113,'NRN12OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:102.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:102.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:102.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:102.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:102.5
  cs,0,0,107),                                                                 &
!  114 RA13OOH  Using EtOOH DD rates
depvel_t(114,'RA13OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:103.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:103.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:103.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:103.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:103.5
  cs,0,0,107),                                                                 &
!  115 RA16OOH  Using EtOOH DD rates
depvel_t(115,'RA16OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:104.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:104.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:104.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:104.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:104.5
  cs,0,0,107),                                                                 &
!  116 RA19OOH  Using EtOOH DD rates
depvel_t(116,'RA19OOH   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:105.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:105.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:105.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:105.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:105.5
  cs,0,0,107),                                                                 &
!  117 RTN28OOH  Using EtOOH DD rates
depvel_t(117,'RTN28OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:106.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:106.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:106.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:106.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:106.5
  cs,0,0,107),                                                                 &
!  118 NRTN28OOH  Using EtOOH DD rates
depvel_t(118,'NRTN28OOH ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:107.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:107.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:107.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:107.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:107.5
  cs,0,0,107),                                                                 &
!  119 RTN26OOH  Using EtOOH DD rates
depvel_t(119,'RTN26OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:108.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:108.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:108.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:108.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:108.5
  cs,0,0,107),                                                                 &
!  120 RTN25OOH  Using EtOOH DD rates
depvel_t(120,'RTN25OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:109.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:109.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:109.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:109.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:109.5
  cs,0,0,107),                                                                 &
!  121 RTN24OOH  Using EtOOH DD rates
depvel_t(121,'RTN24OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:110.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:110.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:110.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:110.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:110.5
  cs,0,0,107),                                                                 &
!  122 RTN23OOH  Using EtOOH DD rates
depvel_t(122,'RTN23OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:111.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:111.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:111.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:111.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:111.5
  cs,0,0,107)                                                                  &
]

! Dry deposition continued
depvel_defs_master(123:150) = [                                                &
!  123 RTN14OOH  Using EtOOH DD rates
depvel_t(123,'RTN14OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:112.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:112.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:112.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:112.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:112.5
  cs,0,0,107),                                                                 &
!  124 RTN10OOH  Using EtOOH DD rates
depvel_t(124,'RTN10OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:113.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:113.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:113.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:113.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:113.5
  cs,0,0,107),                                                                 &
!  125 RTX28OOH  Using EtOOH DD rates
depvel_t(125,'RTX28OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:114.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:114.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:114.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:114.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:114.5
  cs,0,0,107),                                                                 &
!  126 RTX24OOH  Using EtOOH DD rates
depvel_t(126,'RTX24OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:115.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:115.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:115.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:115.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:115.5
  cs,0,0,107),                                                                 &
!  127 RTX22OOH  Using EtOOH DD rates
depvel_t(127,'RTX22OOH  ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:116.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:116.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:116.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:116.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:116.5
  cs,0,0,107),                                                                 &
!  128 NRTX28OOH  Using EtOOH DD rates
depvel_t(128,'NRTX28OOH ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:117.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:117.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:117.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:117.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:117.5
  cs,0,0,107),                                                                 &
!  129 RU12PAN  Using PAN DD rates
depvel_t(129,'RU12PAN   ',                                                     &
[0.01,  0.01,  0.01,  0.01,  0.01,  0.01,&   ! DD:118.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,&   ! DD:118.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,&   ! DD:118.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD:118.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00],& ! DD:118.5
  cs,0,cs2,107),                                                               &
!  130 RTN26PAN  Using PAN DD rates
depvel_t(130,'RTN26PAN  ',                                                     &
[0.01,  0.01,  0.01,  0.01,  0.01,  0.01,&   ! DD:119.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,&   ! DD:119.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,&   ! DD:119.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,&   ! DD:119.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00],& ! DD:119.5
  cs,0,0,107),                                                                 &
!  131 TNCARB12  Using MeCHO DD rates
depvel_t(131,'TNCARB12  ',                                                     &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:120.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD:120.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD:120.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:120.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:120.5
  cs,0,cs2,107),                                                               &
!  132 TNCARB11  Using MeCHO DD rates
depvel_t(132,'TNCARB11  ',                                                     &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:121.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD:121.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD:121.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:121.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:121.5
  cs,0,cs2,107),                                                               &
!  133 RTN23NO3  Using NALD DD rates
depvel_t(133,'RTN23NO3  ',                                                     &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:122.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD:122.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD:122.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:122.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:122.5
  cs,0,cs2,107),                                                               &
!  134 CCARB12  Using MeCHO DD rates
depvel_t(134,'CCARB12   ',                                                     &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:123.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD:123.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD:123.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:123.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:123.5
  cs,0,0,107),                                                                 &
!  135 TNCARB15  Using MeCHO DD rates
depvel_t(135,'TNCARB15  ',                                                     &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:124.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD:124.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD:124.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:124.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:124.5
  cs,0,0,107),                                                                 &
!  136 RCOOH25  Using MeCO2H DD rates
depvel_t(136,'RCOOH25   ',                                                     &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD:125.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD:125.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD:125.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD:125.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD:125.5
  cs,0,0,107),                                                                 &
!  137 TXCARB24  Using MeCHO DD rates
depvel_t(137,'TXCARB24  ',                                                     &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:126.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD:126.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD:126.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:126.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:126.5
  cs,0,0,107),                                                                 &
!  138 TXCARB22  Using MeCHO DD rates
depvel_t(138,'TXCARB22  ',                                                     &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:127.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD:127.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD:127.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:127.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:127.5
  cs,0,0,107),                                                                 &
!  139 RN9NO3  Using ISON DD rates
depvel_t(139,'RN9NO3    ',                                                     &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD:128.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD:128.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD:128.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD:128.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD:128.5
  cs,0,0,107),                                                                 &
!  140 RN12NO3  Using ISON DD rates
depvel_t(140,'RN12NO3   ',                                                     &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD:129.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD:129.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD:129.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD:129.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD:129.5
  cs,0,0,107),                                                                 &
!  141 RN15NO3  Using ISON DD rates
depvel_t(141,'RN15NO3   ',                                                     &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD:130.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD:130.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD:130.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD:130.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD:130.5
  cs,0,0,107),                                                                 &
!  142 RN18NO3  Using ISON DD rates
depvel_t(142,'RN18NO3   ',                                                     &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD:131.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD:131.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD:131.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD:131.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD:131.5
  cs,0,0,107),                                                                 &
!  143 RU14NO3  Using ISON DD rates
depvel_t(143,'RU14NO3   ',                                                     &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD:132.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD:132.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD:132.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD:132.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD:132.5
  cs,0,0,107),                                                                 &
!  144 RTN28NO3  Using ISON DD rates
depvel_t(144,'RTN28NO3  ',                                                     &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD:133.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD:133.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD:133.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD:133.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD:133.5
  cs,0,0,107),                                                                 &
!  145 RTN25NO3  Using NALD DD rates
depvel_t(145,'RTN25NO3  ',                                                     &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:134.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD:134.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD:134.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:134.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:134.5
  cs,0,0,107),                                                                 &
!  146 RTX28NO3  Using ISON DD rates
depvel_t(146,'RTX28NO3  ',                                                     &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,&   ! DD:135.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,&   ! DD:135.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,&   ! DD:135.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,&   ! DD:135.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],& ! DD:135.5
  cs,0,0,107),                                                                 &
!  144 RTX22NO3  Using NALD DD rates
depvel_t(147,'RTX22NO3  ',                                                     &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:136.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD:136.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD:136.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:136.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:136.5
  cs,0,0,107),                                                                 &
!  148 AROH14  Using MeOH DD rates
depvel_t(148,'AROH14    ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:137.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:137.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:137.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:137.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:137.5
  cs,0,0,107),                                                                 &
!  149 ARNOH14  Using MeOH DD rates
depvel_t(149,'ARNOH14   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:138.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:138.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:138.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:138.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:138.5
  cs,0,0,107),                                                                 &
!  150 AROH17  Using MeOH DD rates
depvel_t(150,'AROH17    ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:139.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:139.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:139.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:139.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:139.5
  cs,0,0,107)                                                                  &
  ]

depvel_defs_master(151:n_dry_master) = [                                       &
!  151 ARNOH17  Using MeOH DD rates
depvel_t(151,'ARNOH17   ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,&   ! DD:140.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&   ! DD:140.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&   ! DD:140.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&   ! DD:140.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],& ! DD:140.5
  cs,0,0,107),                                                                 &
!  152 ANHY  Using MeCHO DD rates
depvel_t(152,'ANHY      ',                                                     &
[0.02,  0.02,  0.02,  0.02,  0.02,  0.02,&   ! DD:141.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,&   ! DD:141.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,&   ! DD:141.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&   ! DD:141.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],& ! DD:141.5
  cs,0,0,107),                                                                 &
  ! 153 DHPR12OOH CRIv2.2 addition Using EtOOH DD rates
depvel_t(153,'DHPR12OOH',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,                                      &
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],                                    &
  cs,0,0,119),                                                                 &
!  154 IEPOX CRIv2.2 addition Using MeOH DD rates - could be refined
depvel_t(154,'IEPOX     ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,                                      &
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],                                    &
  cs,0,0,119),                                                                 &
!  155 HMML CRIv2.2 addition Using MeOH DD rates - could be refined
depvel_t(155,'HMML      ',                                                     &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,                                      &
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],                                    &
  cs,0,0,119),                                                                 &
!  156 RU12NO3 CRIv2.2 addition  Using ISON DD rates
depvel_t(156,'RU12NO3   ',                                                     &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,                                      &
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,                                     &
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,                                     &
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,                                     &
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],                                    &
  cs,0,0,119),                                                                 &
!  157 RU10NO3 CRIv2.2 addition  Using ISON DD rates
depvel_t(157,'RU10NO3   ',                                                     &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,                                      &
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,                                     &
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,                                     &
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,                                     &
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],                                    &
  cs,0,0,119),                                                                 &
!  158 HPUCARB12  CRIv2.2 addition  Using EtOOH DD rates
depvel_t(158,'HPUCARB12 ',                                                     &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,                                      &
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,                                     &
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,                                     &
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,                                     &
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],                                    &
  cs,0,0,119),                                                                 &
!  159 HUCARB9 CRIv2.2 addition Using MeOH DD rates - could be refined
depvel_t(159,'HUCARB9  ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,                                      &
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],                                    &
  cs,0,0,119),                                                                 &
!  160 DHPCARB9  CRIv2.2 addition  Using EtOOH DD rates
depvel_t(160,'DHPCARB9 ',                                                      &
[0.50,  0.50,  0.50,  0.25,  0.25,  0.25,                                      &
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,                                     &
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,                                     &
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,                                     &
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03],                                    &
  cs,0,0,119),                                                                 &
!  161 DHCARB9 CRIv2.2 addition Using MeOH DD rates - could be refined
depvel_t(161,'DHCARB9  ',                                                      &
[0.25,  0.25,  0.25,  0.10,  0.10,  0.10,                                      &
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01],                                    &
  cs,0,0,119)                                                                  &
]

! ----------------------------------------------------------------------


! determine which chemistry is to be used. Test here that only one scheme is
! selected.
ierr = 0
chem_scheme = 0
i = 0
IF (ukca_config%l_ukca_strat) THEN
  chem_scheme = s
  i = i + 1
END IF
IF (ukca_config%l_ukca_strattrop) THEN
  chem_scheme = st
  i = i + 1
END IF
IF (ukca_config%l_ukca_trop) THEN
  chem_scheme = t
  i = i + 1
END IF
IF (ukca_config%l_ukca_tropisop) THEN
  chem_scheme = ti
  i = i + 1
END IF
IF (ukca_config%l_ukca_raq) THEN
  chem_scheme = r
  i = i + 1
END IF
IF (ukca_config%l_ukca_offline) THEN
  chem_scheme = ol
  i = i + 1
END IF
IF (ukca_config%l_ukca_cristrat) THEN
  chem_scheme = cs
  i = i + 1
END IF


! Determine which qualifiers to use, based on namelist options.
qualifier = 0
IF (ukca_config%l_ukca_achem .OR. ukca_config%l_ukca_offline) qualifier = a
IF (ukca_config%l_ukca_trophet) qualifier = qualifier + th
IF (ukca_config%l_ukca_het_psc) qualifier = qualifier + hp
IF (ukca_config%l_ukca_stratcfc) qualifier = qualifier + es
IF (ukca_config%l_chem_environ_co2_fld) qualifier = qualifier + ci
IF (ukca_config%i_ukca_hetconfig == 2) qualifier = qualifier + eh
IF (ukca_config%i_ukca_chem_version >= 119  .AND.                              &
    ukca_config%i_ukca_chem == 59) qualifier = qualifier + cs2

! Qualifier for RO2 Non-transport options
IF (ukca_config%l_ukca_ro2_ntp) qualifier = qualifier + rn
! Qualifier for RO2-permutation reactions
IF (ukca_config%l_ukca_ro2_perm) qualifier = qualifier + rp
! Qualifier for CH4 emissions driven StratTrop mechanism
IF (ukca_config%l_ukca_emsdrvn_ch4) qualifier = qualifier + ce
! (dis)qualifier to remove reactions used in previous chemistry versions
! introduced at version 121.
IF (ukca_config%i_ukca_chem_version >= 121) qualifier = qualifier + rm

WRITE (umMessage,'(A,1X,I9)') 'UKCA_CHEM_MASTER: qualifier = ',qualifier
CALL umPrint(umMessage,src=RoutineName)

! No chemistry scheme selected / chemistry scheme not implemented here
IF (chem_scheme == 0) THEN
  ierr = 11
  cmessage = 'No chemistry scheme selected'
END IF
IF (i > 1) THEN
  ierr = 12
  cmessage = "Can't cope with more than one chemistry scheme."
END IF

! calculate sizes of arrays to be extracted
jpspec = 0
jpbk = 0
jppj = 0
jptk = 0
jphk = 0
jpdd = 0
jpdw = 0

! select species and work out number of dry- and wet-deposited species.
IF (ANY(IAND(chch_defs_master%qual,chch_defs_master%disqual) > 0)) THEN
  ierr=30
  cmessage = 'Do not select and deselect the same species.'
END IF

last_chch = chch_t1( 0,'',  0,'','',  0,  0, 0,0,0,0)
DO i=1,n_chch_master
  ! only use the entry if the version is less than or equal to the selected
  ! version
  IF (chch_defs_master(i)%version <= ukca_config%i_ukca_chem_version) THEN
    ! only use if it is part of the chosen chemistry scheme
    IF (IAND(chch_defs_master(i)%scheme,chem_scheme) > 0) THEN
      IF (match(chch_defs_master(i)%qual,chch_defs_master(i)%disqual,          &
          qualifier)) THEN
        ! only use if it is of a higher version than the last chosen item,
        ! if that's for the same item
        take_this = .FALSE.
        IF (chch_defs_master(i)%item_no /= last_chch%item_no) THEN
          take_this = .TRUE.
        ELSE IF (chch_defs_master(i)%version == last_chch%version) THEN
          ! Can't have two entries for the same item number and version
          ierr = 1
          WRITE(cmessage,'(A68,I4)')                                           &
          'Two entries in chch_defs_master with same item and version numbers.'&
          ,chch_defs_master(i)%item_no
        ELSE IF (chch_defs_master(i)%version > last_chch%version) THEN
          ! This is a higher version than one previously selected.
          take_this = .TRUE.
        END IF
        IF (take_this) THEN
          ! Fresh item, no previous version had been selected
          IF (chch_defs_master(i)%item_no /= last_chch%item_no) THEN
            ! the case of this being the first entry of this item number.
            ! Increase counts jpspec, jpdd, and jpdw by 1 as needed.
            jpspec = jpspec + 1
            jpdd = jpdd + chch_defs_master(i)%switch1
            jpdw = jpdw + chch_defs_master(i)%switch2
          ELSE
            ! This is the case of this being a higher version than one
            ! previously taken.
            ! Do not increase jpspec (done that for the older version).
            ! Increase or decrease jpdd and jpdw if there are differences
            ! between the two versions.
            jpdd = jpdd + chch_defs_master(i)%switch1 - last_chch%switch1
            jpdw = jpdw + chch_defs_master(i)%switch2 - last_chch%switch2
          END IF
          ! memorize entry for next round
          last_chch = chch_defs_master(i)
        END IF
      END IF
    END IF
  END IF
END DO


! work out number of bimolecular reactions
last_bimol = ratb_t1(0,'','','','','','',0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0,0)

IF (ANY(IAND(ratb_defs_master%qual,ratb_defs_master%disqual) > 0)) THEN
  ierr=40
  cmessage = 'Do not select and deselect the same bimol reactions.'
END IF
DO i=1,n_bimol_master
  ! only count if part of selected chemistry scheme
  IF (IAND(ratb_defs_master(i)%scheme,chem_scheme) > 0) THEN
    IF (match(ratb_defs_master(i)%qual,ratb_defs_master(i)%disqual,qualifier)) &
       THEN
      ! only count if version is less than or equal to chosen global version
      IF (ratb_defs_master(i)%version <= ukca_config%i_ukca_chem_version) THEN
        ! if previously chosen item has same item number, ignore this
        take_this = .FALSE.
        IF (ratb_defs_master(i)%item_no /= last_bimol%item_no) THEN
          take_this = .TRUE.
        ELSE IF (ratb_defs_master(i)%version == last_bimol%version) THEN
          ! we can't have two entries with the same item number and version.
          ierr = 2
          WRITE (cmessage,'(A68,I4)')                                          &
   'Two entries in ratb_defs_master with same item and version numbers.',      &
            ratb_defs_master(i)%item_no
        ELSE IF (ratb_defs_master(i)%version > last_bimol%version) THEN
          ! higher version
          take_this = .TRUE.
        END IF
        IF (take_this) THEN
          IF (ratb_defs_master(i)%item_no /= last_bimol%item_no) THEN
            jpbk = jpbk + 1
          END IF
          last_bimol = ratb_defs_master(i)
        END IF
      END IF
    END IF
  END IF
END DO

! find number of termolecular reactions
last_termol = ratt_t1(0,'','','','',0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0,0)
IF (ANY(IAND(ratt_defs_master%qual,ratt_defs_master%disqual) > 0)) THEN
  ierr=50
  cmessage = 'Do not select and deselect the same termol reactions.'
END IF
DO i=1,n_ratt_master
  IF (IAND(ratt_defs_master(i)%scheme,chem_scheme) > 0) THEN
    IF (match(ratt_defs_master(i)%qual,ratt_defs_master(i)%disqual,qualifier)) &
    THEN
      IF (ratt_defs_master(i)%version <= ukca_config%i_ukca_chem_version) THEN
        take_this = .FALSE.
        IF (ratt_defs_master(i)%item_no /= last_termol%item_no) THEN
          take_this = .TRUE.
        ELSE IF (ratt_defs_master(i)%version == last_termol%version) THEN
          ierr = 3
          cmessage =                                                           &
          'Two items in ratt_defs_master with same item and version numbers.'
        ELSE IF (ratt_defs_master(i)%version > last_termol%version) THEN
          take_this = .TRUE.
        END IF
        IF (take_this) THEN
          IF (ratt_defs_master(i)%item_no /= last_termol%item_no) THEN
            jptk = jptk + 1
          END IF
          last_termol = ratt_defs_master(i)
        END IF
      END IF
    END IF
  END IF
END DO

! find number of photolysis reactions
last_photol = ratj_t1(0,'','','','','','',0.0,0.0,0.0,0.0,0.0,'',0,0,0,0)
IF (ANY(IAND(ratj_defs_master%qual,ratj_defs_master%disqual) > 0)) THEN
  ierr=60
  cmessage = 'Do not select and deselect the same photol reactions.'
END IF
DO i=1,n_ratj_master
  IF (IAND(ratj_defs_master(i)%scheme,chem_scheme) > 0) THEN
    IF (match(ratj_defs_master(i)%qual,ratj_defs_master(i)%disqual,qualifier)) &
    THEN
      IF (ratj_defs_master(i)%version <= ukca_config%i_ukca_chem_version) THEN
        take_this = .FALSE.
        IF (ratj_defs_master(i)%item_no /= last_photol%item_no) THEN
          take_this = .TRUE.
        ELSE IF (ratj_defs_master(i)%version == last_photol%version) THEN
          ierr = 4
          WRITE(cmessage,'(A58,I4)')                                           &
            'Two items in ratj_defs_master with same item and version ',       &
            ratj_defs_master(i)%item_no
        ELSE IF (ratj_defs_master(i)%version > last_photol%version) THEN
          take_this = .TRUE.
        END IF
        IF (take_this) THEN
          IF (ratj_defs_master(i)%item_no /= last_photol%item_no)              &
            jppj = jppj + 1
          last_photol = ratj_defs_master(i)
        END IF
      END IF
    END IF
  END IF
END DO

! find number of heterogeneous reactions
last_het = rath_t1(0,'','','','','','',0.0,0.0,0.0,0.0,0,0,0,0)
IF (ANY(IAND(rath_defs_master%qual,rath_defs_master%disqual) > 0)) THEN
  ierr=70
  cmessage = 'Do not select and deselect the same hetero reactions.'
END IF
DO i=1,n_het_master
  IF (IAND(rath_defs_master(i)%scheme,chem_scheme) > 0) THEN
    IF (match(rath_defs_master(i)%qual,rath_defs_master(i)%disqual,qualifier)) &
    THEN
      IF (rath_defs_master(i)%version <= ukca_config%i_ukca_chem_version) THEN
        take_this = .FALSE.
        IF (rath_defs_master(i)%item_no /= last_het%item_no) THEN
          take_this = .TRUE.
        ELSE IF (rath_defs_master(i)%version == last_het%version) THEN
          ierr = 5
          cmessage =                                                           &
          'Two items in rath_defs_master with same item and version numbers.'
        ELSE IF (rath_defs_master(i)%version > last_het%version) THEN
          take_this = .TRUE.
        END IF
        IF (take_this) THEN
          IF (rath_defs_master(i)%item_no /= last_het%item_no) THEN
            jphk = jphk + 1
          END IF
          last_het = rath_defs_master(i)
        END IF
      END IF
    END IF
  END IF
END DO

! calculate total number of reactions
jpnr = jpbk + jptk + jppj + jphk

! allocate final arrays for species and reactions
ALLOCATE(chch_defs_new(jpspec))
ALLOCATE(ratb_defs_new(jpbk))
ALLOCATE(ratj_defs_new(jppj))
ALLOCATE(ratt_defs_new(jptk))
ALLOCATE(rath_defs_new(jphk))
ALLOCATE(depvel_defs_new(jddept,jddepc,jpdd))
ALLOCATE(henry_defs_new(6,jpdw))

! pick out reactions/species needed
j = 0
last_chch = chch_t1( 0,'',  0,'','',  0,  0, 0,0,0,0)
DO i=1,n_chch_master
  ! only use the entry if the version is less than or equal to the selected
  ! version
  IF ((chch_defs_master(i)%version <= ukca_config%i_ukca_chem_version) .AND.   &
   (IAND(chch_defs_master(i)%scheme,chem_scheme) > 0) .AND.                    &
   (match(chch_defs_master(i)%qual,chch_defs_master(i)%disqual,qualifier)))    &
   THEN
    ! increase count only if fresh item not just different version of previous
    ! item. If not, overwrite previously chosen item.
    IF (chch_defs_master(i)%item_no /= last_chch%item_no) j = j + 1
    ! If new item or higher version of previously selected item, copy into array
    IF ((chch_defs_master(i)%item_no /= last_chch%item_no) .OR.                &
        (chch_defs_master(i)%version > last_chch%version)) THEN
      ! memorize entry for next round
      chch_defs_new(j) = chch_defs_master(i)
      last_chch = chch_defs_master(i)
    END IF
  END IF
END DO

WRITE (umMessage,'(A9,I4)') 'jpspec = ',jpspec
CALL umPrint(umMessage,src=RoutineName)
WRITE (umMessage,'(A4)') 'CHCH'
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jpspec
  WRITE (umMessage,'(I3,1X,A10,1X,A2,2(1X,I1),1X,I3)')                         &
         chch_defs_new(i)%item_no,                                             &
         chch_defs_new(i)%speci,chch_defs_new(i)%ctype,                        &
         chch_defs_new(i)%switch1,chch_defs_new(i)%switch2,                    &
         chch_defs_new(i)%version
  CALL umPrint(umMessage,src=RoutineName)
END DO

! Pick out bimolecular reactions needed
j=0
last_bimol = ratb_t1(0,'','','','','','',0.0,0.0,0.0, 0.0, 0.0, 0.0, 0.0,0,0,0,0)
DO i=1,n_bimol_master
  IF ((IAND(ratb_defs_master(i)%scheme,chem_scheme) > 0) .AND.                 &
    (ratb_defs_master(i)%version <= ukca_config%i_ukca_chem_version) .AND.     &
    (match(ratb_defs_master(i)%qual,ratb_defs_master(i)%disqual,qualifier)))   &
  THEN
    IF (ratb_defs_master(i)%item_no /= last_bimol%item_no) j = j + 1
    IF ((ratb_defs_master(i)%item_no /= last_bimol%item_no) .OR.               &
        (ratb_defs_master(i)%version > last_bimol%version)) THEN
      ratb_defs_new(j) = ratb_defs_master(i)
      last_bimol = ratb_defs_master(i)
    END IF
  END IF
END DO

! Pick out termolecular reactions needed
j=0
last_termol = ratt_t1(0,'','','','',0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, 0.0,0,0,0,0)
DO i=1,n_ratt_master
  IF ((IAND(ratt_defs_master(i)%scheme,chem_scheme) > 0) .AND.                 &
    (ratt_defs_master(i)%version <= ukca_config%i_ukca_chem_version) .AND.     &
    (match(ratt_defs_master(i)%qual,ratt_defs_master(i)%disqual,qualifier)))   &
  THEN
    IF (ratt_defs_master(i)%item_no /= last_termol%item_no) j = j + 1
    IF ((ratt_defs_master(i)%item_no /= last_termol%item_no) .OR.              &
      (ratt_defs_master(i)%version > last_termol%version)) THEN
      ratt_defs_new(j) = ratt_defs_master(i)
      last_termol = ratt_defs_master(i)
    END IF
  END IF
END DO

! Pick out photolysis reactions needed
j=0
last_photol = ratj_t1(0,'','','','','','',0.0,0.0,0.0,0.0,0.0,'',0,0,0,0)
DO i=1,n_ratj_master
  IF ((IAND(ratj_defs_master(i)%scheme,chem_scheme) > 0) .AND.                 &
    (ratj_defs_master(i)%version <= ukca_config%i_ukca_chem_version) .AND.     &
    (match(ratj_defs_master(i)%qual,ratj_defs_master(i)%disqual,qualifier)))   &
  THEN
    IF (ratj_defs_master(i)%item_no /= last_photol%item_no) j = j + 1
    IF ((ratj_defs_master(i)%item_no /= last_photol%item_no) .OR.              &
        (ratj_defs_master(i)%version > last_photol%version)) THEN
      ratj_defs_new(j) = ratj_defs_master(i)
      last_photol = ratj_defs_master(i)
    END IF
  END IF
END DO

! Pick out heterogeneous reactions needed
j=0
last_het = rath_t1(0,'','','','','','',0.0,0.0,0.0,0.0,0,0,0,0)
DO i=1,n_het_master
  IF ((IAND(rath_defs_master(i)%scheme,chem_scheme) > 0) .AND.                 &
    (rath_defs_master(i)%version <= ukca_config%i_ukca_chem_version) .AND.     &
    (match(rath_defs_master(i)%qual,rath_defs_master(i)%disqual,qualifier)))   &
  THEN
    IF (rath_defs_master(i)%item_no /= last_het%item_no) j = j + 1
    IF ((rath_defs_master(i)%item_no /= last_het%item_no) .OR.                 &
        (rath_defs_master(i)%version > last_het%version)) THEN
      rath_defs_new(j) = rath_defs_master(i)
      last_het = rath_defs_master(i)
    END IF
  END IF
END DO


! pick out species with dry deposition
j = 0
last_depvel = depvel_t(0,'',[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,          &
                              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
                              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],0,0,0,0)
DO i=1,jpspec
  IF (chch_defs_new(i)%switch1 == 1) THEN
    found = .FALSE.
    j = j + 1
    DO k=1,n_dry_master
      ! check for matching species, chemistry scheme, and version
      IF (depvel_defs_master(k)%speci == chch_defs_new(i)%speci) THEN
        IF (IAND(depvel_defs_master(k)%scheme,chem_scheme) > 0) THEN
          IF (depvel_defs_master(k)%version <=                                 &
              ukca_config%i_ukca_chem_version) THEN
            IF (match(depvel_defs_master(k)%qual,depvel_defs_master(k)%disqual,&
                qualifier)) THEN
              take_this = .FALSE.
              IF (depvel_defs_master(k)%item_no == last_depvel%item_no) THEN
                ! can't have two entries at the same version
                IF (depvel_defs_master(k)%version == last_depvel%version) THEN
                  ierr = 7
                  WRITE (cmessage,'(A70,I4)')                                  &
    'Two entries in depvel_defs_master with same item and version numbers.',   &
                    depvel_defs_master(k)%speci
                END IF
                IF (depvel_defs_master(k)%version > last_depvel%version)       &
                  take_this = .TRUE.
              ELSE
                take_this = .TRUE.
              END IF
              IF (take_this) THEN
                depvel_defs_new(:,:,j) = RESHAPE(                              &
                  depvel_defs_master(k)%velocity,[jddept,jddepc])
                last_depvel =  depvel_defs_master(k)
                found = .TRUE.
              END IF
            END IF
          END IF
        END IF
      END IF
    END DO
    ! check for whether dry deposition parameters have been located.
    IF (.NOT. (found)) THEN
      ierr = 8
      cmessage = 'Dry deposition selected for '//chch_defs_new(i)%speci//      &
        ' but dry deposition parameters missing.'
    END IF
  END IF
END DO

! pick out species with wet deposition
j = 0
last_wetdep = wetdep(0,'',[0.0,0.0,0.0,0.0,0.0,0.0],0,0,0,0)
DO i=1,jpspec
  IF (chch_defs_new(i)%switch2 == 1) THEN
    found = .FALSE.
    j = j + 1
    DO k=1,n_wet_master
      ! check for matching species, chemistry scheme, and version
      IF (henry_defs_master(k)%speci == chch_defs_new(i)%speci) THEN
        IF (IAND(henry_defs_master(k)%scheme,chem_scheme) > 0) THEN
          IF (henry_defs_master(k)%version <=                                  &
              ukca_config%i_ukca_chem_version) THEN
            IF (match(henry_defs_master(k)%qual,henry_defs_master(k)%disqual,  &
                qualifier)) THEN
              take_this = .FALSE.
              IF (henry_defs_master(k)%item_no == last_wetdep%item_no) THEN
                ! can't have two entries at the same version
                IF (henry_defs_master(k)%version == last_wetdep%version) THEN
                  ierr = 9
                  cmessage =                                                   &
    'Two entries in henry_defs_master with same item and version numbers.'
                END IF
                IF (henry_defs_master(k)%version > last_wetdep%version)        &
                  take_this = .TRUE.
              ELSE
                take_this = .TRUE.
              END IF
              IF (take_this) THEN
                henry_defs_new(:,j) = henry_defs_master(k)%params
                last_wetdep =  henry_defs_master(k)
                found = .TRUE.
              END IF
            END IF
          END IF
        END IF
      END IF
    END DO
    ! check for whether dry deposition parameters have been located
    IF (.NOT. (found)) THEN
      ierr = 10
      cmessage = 'Wet deposition selected for '//chch_defs_new(i)%speci//      &
        ' but wet deposition parameters missing.'
    END IF
  END IF
END DO

! Calculate total number of tracers
jpctr = 0  ! Total number of tracers ('TR')
jpro2 = 0  ! Total number of peroxy radicals ('OO')
jpcspf = 0 ! Total number of chemically active species (TR+OO)
DO i=1,jpspec
  IF (TRIM(chch_defs_new(i)%ctype) == 'TR') THEN
    jpctr = jpctr + 1
    jpcspf = jpcspf + 1
    ! Code to count the number of RO2 species
  ELSE IF (TRIM(chch_defs_new(i)%ctype) == 'OO') THEN
    jpro2 = jpro2 + 1
    jpcspf = jpcspf + 1
    ! If we are transporting RO2 species, also add to tracer array
    IF (.NOT. ukca_config%l_ukca_ro2_ntp) THEN
      jpctr = jpctr + 1
    END IF
  END IF
END DO

WRITE(umMessage,'(A8,I4)') 'jpctr  = ',jpctr
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A8,I4)') 'jpro2  = ',jpro2
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A8,I4)') 'jpcspf  = ',jpcspf
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A8,I4)') 'jpbk   = ',jpbk
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A4)') 'ratb'
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jpbk
  WRITE (umMessage,'(I3,1X,6(A10,1X))')                                        &
                       ratb_defs_new(i)%item_no,ratb_defs_new(i)%react1,       &
                       ratb_defs_new(i)%react2,ratb_defs_new(i)%prod1,         &
                       ratb_defs_new(i)%prod2,ratb_defs_new(i)%prod3,          &
                       ratb_defs_new(i)%prod4
  CALL umPrint(umMessage,src=RoutineName)
  WRITE (umMessage,'(I3,3(1X,E15.6))')                                         &
                       ratb_defs_new(i)%version,                               &
                       ratb_defs_new(i)%k0   ,ratb_defs_new(i)%alpha,          &
                       ratb_defs_new(i)%beta
  CALL umPrint(umMessage,src=RoutineName)
END DO
WRITE (umMessage,'(A8,I4)') 'jppj   = ',jppj
CALL umPrint(umMessage,src=RoutineName)
WRITE (umMessage,'(A4)') 'ratj'
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jppj
  WRITE (umMessage,'(I3,1X,7(A10,1X),F9.4,1X,I3)')ratj_defs_new(i)%item_no,    &
    ratj_defs_new(i)%react1,ratj_defs_new(i)%react2,ratj_defs_new(i)%prod1,    &
    ratj_defs_new(i)%prod2,ratj_defs_new(i)%prod3,ratj_defs_new(i)%prod4,      &
    ratj_defs_new(i)%fname,ratj_defs_new(i)%jfacta,                            &
    ratj_defs_new(i)%version
  CALL umPrint(umMessage,src='ukca_chem_master')
END DO
WRITE (umMessage,'(A8,I4)') 'jptk   = ',jptk
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A4)') 'ratt'
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jptk
  WRITE(umMessage,'(I3,1X,4(A10,1X))') ratt_defs_new(i)%item_no,               &
    ratt_defs_new(i)%react1,ratt_defs_new(i)%react2,ratt_defs_new(i)%prod1,    &
    ratt_defs_new(i)%prod2
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(2(E15.6,1X,F10.6,1X,F12.6,1X),I3)')                        &
    ratt_defs_new(i)%k1,                                                       &
    ratt_defs_new(i)%alpha1,ratt_defs_new(i)%beta1,ratt_defs_new(i)%k2,        &
    ratt_defs_new(i)%alpha2,ratt_defs_new(i)%beta2,ratt_defs_new(i)%version
  CALL umPrint(umMessage,src=RoutineName)
END DO
WRITE(umMessage,'(A8,I4)') 'jphk   = ',jphk
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A4)') 'rath'
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jphk
  WRITE(umMessage,'(I3,1X,6(A10,1X),I3)') rath_defs_new(i)%item_no,            &
    rath_defs_new(i)%react1,rath_defs_new(i)%react2,rath_defs_new(i)%prod1,    &
    rath_defs_new(i)%prod2,rath_defs_new(i)%prod3,rath_defs_new(i)%prod4,      &
    rath_defs_new(i)%version
  CALL umPrint(umMessage,src=RoutineName)
END DO

WRITE(umMessage,'(A8,I4)') 'jpdd   = ',jpdd
CALL umPrint(umMessage,src=RoutineName)
j=0
DO i=1,jpspec
  IF (chch_defs_new(i)%switch1 == 1) THEN
    j=j+1
    WRITE(umMessage,'(A10)') chch_defs_new(i)%speci
    CALL umPrint(umMessage,src=RoutineName)
    DO k=1,jddepc
      WRITE(umMessage,'(6(F10.6))') depvel_defs_new(:,k,j)
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF
END DO

j=0
WRITE(umMessage,'(A8,I4)') 'jpdw   = ',jpdw
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jpspec
  IF (chch_defs_new(i)%switch2 == 1) THEN
    j=j+1
    WRITE(umMessage,'(A10)') chch_defs_new(i)%speci
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(6(E15.6))') henry_defs_new(:,j)
    CALL umPrint(umMessage,src=RoutineName)
  END IF
END DO

IF (ierr > 1) THEN
  CALL ereport(RoutineName,ierr,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chem_master


! ######################################################################

LOGICAL FUNCTION match(qual, disqual, qualifier)

! If qual == 0 and disqual == 0, select unconditionally.
! If qual > 0, only select if the corresponding flag match is True.
! If disqual > 0, select if the corresponding flag is False,
!                 deselect if the corresponding flag is True.
IMPLICIT NONE

INTEGER, INTENT(IN) :: qual, disqual, qualifier

! initialise to False
match=.FALSE.

! set to True if qual==0, disqual value may change this later
match = (qual == 0)

IF (qual > 0) THEN
  ! if qual in qualifier, take reaction
  IF (IAND(qual,  qualifier) > 0) match = .TRUE.
END IF

IF (disqual > 0) THEN
  ! if disqual in qualifier, reject reaction
  IF (IAND(disqual,qualifier) > 0) match = .FALSE.
END IF

RETURN
END FUNCTION match


END MODULE ukca_chem_master_mod
