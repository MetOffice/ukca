! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module providing UKCA with inputs describing RAQ-AERO chemistry
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
!   Language:  Fortran 95
!   This code is written to UMDP3 standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_chem_raqaero_mod

USE ukca_chem_defs_mod, ONLY:  chch_t, ratb_t, rath_t, ratj_t, ratt_t
IMPLICIT NONE
PRIVATE

INTEGER, PARAMETER, PUBLIC :: ndry_raqaero = 19    ! No of dry deposited species
INTEGER, PARAMETER, PUBLIC :: nwet_raqaero = 24    ! No of wet deposited species


! Regional air quality (RAQ) chemistry, plus aerosol chem
! ========================================================
!
TYPE(chch_t), PUBLIC :: chch_defs_raqaero( 66)=[                               &
chch_t(  1,'O(3P)     ',  1,'SS        ','          ',  0,  0),                &
chch_t(  2,'O(1D)     ',  1,'SS        ','          ',  0,  0),                &
chch_t(  3,'OH        ',  1,'SS        ','          ',  0,  0),                &
chch_t(  4,'O3        ',  1,'TR        ','          ',  1,  1),                &
chch_t(  5,'NO        ',  1,'TR        ','          ',  0,  0),                &
chch_t(  6,'NO3       ',  1,'TR        ','          ',  0,  1),                &
chch_t(  7,'NO2       ',  1,'TR        ','          ',  1,  0),                &
chch_t(  8,'N2O5      ',  2,'TR        ','          ',  0,  1),                &
chch_t(  9,'HO2NO2    ',  1,'TR        ','          ',  0,  1),                &
chch_t( 10,'HONO2     ',  1,'TR        ','          ',  1,  1),                &
chch_t( 11,'H2O2      ',  1,'TR        ','          ',  1,  1),                &
chch_t( 12,'CH4       ',  1,'TR        ','          ',  1,  0),                &
chch_t( 13,'CO        ',  1,'TR        ','          ',  1,  0),                &
chch_t( 14,'HCHO      ',  1,'TR        ','          ',  0,  1),                &
chch_t( 15,'HO2       ',  1,'SS        ','          ',  0,  1),                &
chch_t( 16,'MeOO      ',  1,'SS        ','          ',  0,  1),                &
chch_t( 17,'MeOOH     ',  1,'TR        ','          ',  1,  1),                &
chch_t( 18,'EtOO      ',  1,'SS        ','          ',  0,  0),                &
chch_t( 19,'C2H6      ',  1,'TR        ','          ',  0,  0),                &
chch_t( 20,'MeCO3     ',  1,'SS        ','          ',  0,  0),                &
chch_t( 21,'EtOOH     ',  1,'TR        ','          ',  1,  1),                &
chch_t( 22,'MeCHO     ',  1,'TR        ','          ',  0,  0),                &
chch_t( 23,'PAN       ',  1,'TR        ','          ',  1,  0),                &
chch_t( 24,'s-BuOO    ',  1,'SS        ','          ',  0,  0),                &
chch_t( 25,'C3H8      ',  1,'TR        ','          ',  0,  0),                &
chch_t( 26,'i-PrOOH   ',  1,'TR        ','          ',  1,  1),                &
chch_t( 27,'Me2CO     ',  1,'TR        ','          ',  0,  0),                &
chch_t( 28,'O3S       ',  1,'TR        ','          ',  1,  0),                &
chch_t( 29,'C5H8      ',  1,'TR        ','          ',  0,  0),                &
chch_t( 30,'i-PrOO    ',  1,'SS        ','          ',  0,  0),                &
chch_t( 31,'ISOOH     ',  1,'TR        ','          ',  1,  1),                &
chch_t( 32,'ISON      ',  1,'TR        ','          ',  0,  1),                &
chch_t( 33,'MGLY      ',  1,'TR        ','          ',  0,  1),                &
chch_t( 34,'MVK       ',  1,'TR        ','          ',  0,  0),                &
chch_t( 35,'MVKOOH    ',  1,'TR        ','          ',  1,  1),                &
chch_t( 36,'MeCOCH2OO ',  1,'SS        ','          ',  0,  0),                &
chch_t( 37,'MEKO2     ',  1,'SS        ','          ',  0,  0),                &
chch_t( 38,'HOC2H4O2  ',  1,'SS        ','          ',  0,  0),                &
chch_t( 39,'ORGNIT    ',  1,'TR        ','          ',  1,  1),                &
chch_t( 40,'HOC3H6O2  ',  1,'SS        ','          ',  0,  0),                &
chch_t( 41,'CH3OH     ',  1,'TR        ','          ',  0,  1),                &
chch_t( 42,'OXYL1     ',  1,'SS        ','          ',  0,  0),                &
chch_t( 43,'H2        ',  1,'TR        ','          ',  1,  0),                &
chch_t( 44,'DMS       ',  1,'TR        ','          ',  0,  0),                &
chch_t( 45,'SO2       ',  1,'TR        ','          ',  1,  1),                &
chch_t( 46,'H2SO4     ',  1,'TR        ','          ',  1,  0),                &
chch_t( 47,'MSA       ',  1,'TR        ','          ',  0,  0),                &
chch_t( 48,'DMSO      ',  1,'TR        ','          ',  0,  1),                &
chch_t( 49,'NH3       ',  1,'TR        ','          ',  1,  1),                &
chch_t( 50,'RNC3H6    ',  1,'TR        ','          ',  0,  0),                &
chch_t( 51,'C2H4      ',  1,'TR        ','          ',  0,  0),                &
chch_t( 52,'MEMALD1   ',  1,'SS        ','          ',  0,  0),                &
chch_t( 53,'RNC2H4    ',  1,'TR        ','          ',  0,  0),                &
chch_t( 54,'HOIPO2    ',  1,'SS        ','          ',  0,  0),                &
chch_t( 55,'HOMVKO2   ',  1,'SS        ','          ',  0,  0),                &
chch_t( 56,'Monoterp  ',  1,'TR        ','          ',  0,  0),                &
chch_t( 57,'Sec_Org   ',  1,'TR        ','          ',  0,  1),                &
chch_t( 58,'C3H6      ',  1,'TR        ','          ',  0,  0),                &
chch_t( 59,'C4H10     ',  1,'TR        ','          ',  0,  0),                &
chch_t( 60,'s-BuOOH   ',  1,'TR        ','          ',  1,  1),                &
chch_t( 61,'MEK       ',  1,'TR        ','          ',  0,  0),                &
chch_t( 62,'TOLUENE   ',  1,'TR        ','          ',  0,  0),                &
chch_t( 63,'TOLP1     ',  1,'SS        ','          ',  0,  0),                &
chch_t( 64,'MEMALD    ',  1,'TR        ','          ',  0,  0),                &
chch_t( 65,'GLY       ',  1,'TR        ','          ',  0,  1),                &
chch_t( 66,'oXYLENE   ',  1,'TR        ','          ',  0,  0)                 &
 ]


! The RAQ-AERO Chemistry, only uses Backward Euler
! solver. Because of that ratb_defs_raqaero and ratt_defs_raqaero are
! initialised with empty names for the reactants/products and zero
! values for the rate coeffs and products' fractional yields.
! This could be modified in future UM releases if the
! use of ASAD routines is considered to do RAQ-AERO
! chemistry integration with the Newton-Raphson solver.
!

TYPE(ratb_t), PARAMETER, PUBLIC :: ratb_defs_raqaero( 123) =                   &
ratb_t('          ','          ','          ','          ','          ',       &
       '          ',0.00e+00,  0.00,      0.00, 0.000, 0.000, 0.000, 0.000)


TYPE(ratt_t), PARAMETER, PUBLIC :: ratt_defs_raqaero(  14) =                   &
ratt_t('          ','          ','          ','          ',                    &
       0.0, 0.00e+00, 0.00,    0.0, 0.00e+00,  0.00,  0.0, 0.000, 0.000)

TYPE(rath_t), PUBLIC :: rath_defs_raqaero(2)=                                  &
rath_t('          ','          ','          ','          ','          ',       &
       '          ', 0.0, 0.0, 0.0, 0.0)

!
TYPE(ratj_t), PUBLIC :: ratj_defs_raqaero(  23)=[                              &
ratj_t('O3        ','PHOTON    ','O2        ','O(3P)     ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jo3b      ') ,              &
ratj_t('O3        ','PHOTON    ','O2        ','O(1D)     ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jo3a      ') ,              &
ratj_t('NO2       ','PHOTON    ','NO        ','O(3P)     ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jno2      ') ,              &
ratj_t('H2O2      ','PHOTON    ','OH        ','OH        ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jh2o2     ') ,              &
ratj_t('HONO2     ','PHOTON    ','OH        ','NO2       ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jhno3     ') ,              &
ratj_t('HCHO      ','PHOTON    ','HO2       ','HO2       ','CO        ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jhchoa    ') ,              &
ratj_t('HCHO      ','PHOTON    ','H2        ','CO        ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jhchob    ') ,              &
ratj_t('MeCHO     ','PHOTON    ','MeOO      ','HO2       ','CO        ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jaceta    ') ,              &
ratj_t('MEK       ','PHOTON    ','MeCO3     ','EtOO      ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jaceto    ') ,              &
ratj_t('Me2CO     ','PHOTON    ','MeCO3     ','MeOO      ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jaceto    ') ,              &
ratj_t('HO2NO2    ','PHOTON    ','HO2       ','NO2       ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jpna      ') ,              &
ratj_t('MGLY      ','PHOTON    ','MeCO3     ','CO        ','HO2       ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmkal     ') ,              &
ratj_t('GLY       ','PHOTON    ','CO        ','CO        ','HO2       ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmkal     ') ,              &
ratj_t('NO3       ','PHOTON    ','NO        ','O2        ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jno3a     ') ,              &
ratj_t('NO3       ','PHOTON    ','NO2       ','O(3P)     ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jno3b     ') ,              &
ratj_t('N2O5      ','PHOTON    ','NO3       ','NO2       ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jn2o5     ') ,              &
ratj_t('MeOOH     ','PHOTON    ','HCHO      ','HO2       ','OH        ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmhp      ') ,              &
ratj_t('PAN       ','PHOTON    ','MeCO3     ','NO2       ','          ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jpan      ') ,              &
ratj_t('EtOOH     ','PHOTON    ','MeCHO     ','HO2       ','OH        ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmhp      ') ,              &
ratj_t('i-PrOOH   ','PHOTON    ','Me2CO     ','HO2       ','OH        ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmhp      ') ,              &
ratj_t('s-BuOOH   ','PHOTON    ','MEK       ','HO2       ','OH        ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmhp      ') ,              &
ratj_t('ISOOH     ','PHOTON    ','MVK       ','HCHO      ','OH        ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmhp      ') ,              &
ratj_t('MVKOOH    ','PHOTON    ','MGLY      ','HCHO      ','OH        ',       &
  '          ',   0.0,   0.0,   0.0,   0.0, 100.0,'jmhp      ')                &
  ]


! Dry deposition array:
! Dry deposition velocities are written in a table format, with the columns
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
! '24 hr' corresponds to the 24 hour average deposition velocity, currently not
! used.
! The units of x are cm s-1

REAL, PUBLIC :: depvel_defs_raqaero(  6,  5, ndry_raqaero)=RESHAPE([           &
  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  & !O3
  1.00,  0.11,  0.56,  0.26,  0.11,  0.19,                                     &
  1.00,  0.37,  0.69,  0.59,  0.46,  0.53,                                     &
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,                                     &
  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,                                     &
  0.02,  0.02,  0.02,  0.02,  0.02,  0.02,  & !NO2
  0.83,  0.04,  0.44,  0.06,  0.04,  0.05,                                     &
  0.63,  0.06,  0.35,  0.07,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                                     &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !HONO2
  4.00,  3.00,  3.50,  2.00,  1.00,  1.50,                                     &
  2.50,  1.50,  2.00,  1.00,  1.00,  1.00,                                     &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,                                     &
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,                                     &
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  & !H2O2
  1.25,  0.16,  0.71,  0.28,  0.12,  0.20,                                     &
  1.25,  0.53,  0.89,  0.83,  0.78,  0.81,                                     &
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,                                     &
  0.32,  0.32,  0.32,  0.32,  0.32,  0.32,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !CH4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !CO
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                                     &
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & !MeOOH
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                                     &
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & !EtOOH
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  & !PAN
  0.53,  0.04,  0.29,  0.06,  0.05,  0.06,                                     &
  0.42,  0.06,  0.24,  0.07,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00,                                     &
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & !i-PrOOH
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                                     &
  0.07,  0.07,  0.07,  0.07,  0.07,  0.07,  & !O3S
  1.00,  0.11,  0.56,  0.26,  0.11,  0.19,                                     &
  1.00,  0.37,  0.69,  0.59,  0.46,  0.53,                                     &
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,                                     &
  0.07,  0.07,  0.07,  0.07,  0.07,  0.07,                                     &
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & !ISOOH
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                                     &
  0.36,  0.36,  0.36,  0.19,  0.19,  0.19,  & !MVKOOH
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,                                     &
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02,                                     &
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & !ORGNIT
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  & !H2
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,                                     &
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
  0.80,  0.80,  0.80,  0.80,  0.80,  0.80,  & ! NH3
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,                                     &
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,                                     &
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,                                     &
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,                                     &
  0.25,  0.25,  0.25,  0.10,  0.10,  0.10,  & !s-BuOOH
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,                                     &
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,                                     &
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,                                     &
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01                                      &
  ],[  6,  5, ndry_raqaero] )


! The following formula is used to calculate the effective Henry's Law coef,
! which takes the effects of dissociation and complex formation on a species'
! solubility (see Giannakopoulos, 1998)
!
!       H(eff) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
!
! The data in columns 1 and 2 above give the data for this gas-aqueous transfer,
!       Column 1 = K(298) [M/atm]
!       Column 2 = -deltaH/R [K-1]
!
! If the species dissociates in the aqueous phase the above term is
! multiplied by another factor of 1+{K(aq)/[H+]}, where
!       K(aq) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
! The data in columns 3 and 4 give the data for this aqueous-phase dissociation,
!       Column 3 = K(298) [M]
!       Column 4 = -deltaH/R [K-1]
!
! The data in columns 5 and 6 give the data for a second dissociation,
! see e.g for SO2, HSO3^{-}, and SO3^{2-} in module UKCA_CHEM_AER
!       Column 5 = K(298) [M]
!       Column 6 = -deltaH/R [K-1]

REAL, PUBLIC :: henry_defs_raqaero(6, nwet_raqaero)=RESHAPE([                  &
 0.113e-01, 0.230e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! O3
 0.200e+01, 0.200e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! NO3
 0.210e+06, 0.870e+04, 0.200e+02, 0.000e+00, 0.000e+00, 0.000e+00,&  ! N2O5
 0.130e+05, 0.690e+04, 0.100e-04, 0.000e+00, 0.000e+00, 0.000e+00,&  ! HO2NO2
 0.210e+06, 0.870e+04, 0.200e+02, 0.000e+00, 0.000e+00, 0.000e+00,&  ! HONO2
 0.830e+05, 0.740e+04, 0.240e-11,-0.373e+04, 0.000e+00, 0.000e+00,&  ! H2O2
 0.330e+04, 0.650e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! HCHO
 0.400e+04, 0.590e+04, 0.200e-04, 0.000e+00, 0.000e+00, 0.000e+00,&  ! HO2
 0.200e+04, 0.660e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! MeOO
 0.310e+03, 0.500e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! MeOOH
 0.340e+03, 0.570e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! EtOOH
 0.340e+03, 0.570e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! i-PrOOH
 0.170e+07, 0.970e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! ISOOH
 0.300e+04, 0.740e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! ISON
 0.350e+04, 0.720e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! MGLY
 0.170e+07, 0.970e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! MVKOOH
 0.130e+03, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! ORGNIT
 0.220e+03, 0.520e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! CH3OH
 0.123e+01, 0.302e+04, 0.123e-01, 0.201e+04, 0.600e-07, 0.112e+04,&  ! SO2
 0.500e+05, 0.642e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! DMSO
 0.100e+07, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! NH3 (H*)
 0.100e+06, 0.120e+02, 0.000e-00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! Sec_Org
 0.340e+03, 0.570e+04, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,&  ! s-BuOOH
 0.360e+06, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00 &  ! GLY
  ],[  6, nwet_raqaero] )

END MODULE ukca_chem_raqaero_mod
