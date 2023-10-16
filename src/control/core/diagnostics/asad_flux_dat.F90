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
! Description:
!  Module containing the reactions which define the chemical reaction
!  fluxes for STASH output.
!
!  Part of the UKCA model. UKCA is a community model supported
!  by The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!    This code is written to UMDP3 v8 programming standards.
!
! ######################################################################
!
! TP       = RXN for reaction
!            DEP for deposition
!            EMS for emission
!            NET for a tendency of a particular species
!            STE for the Strat-Trop exchange across the 2PVU+380K Tropopause
!                for a particular speices
!            MAS for the mass of air per gridcell
!            PSC for location of type 1/2 PSCs
! STASH#   = First 2 digits are section number, last 3 are item number,
!            5 digits in total, e.g. 34412 etc.
! WHICH    = IF TP==RXN: B = bimolecular         3D
!                        H = heterogeneous       3D
!                        T = termolecular        3D
!                        J = photolysis          3D
!            IF TP==DEP: D = dry deposition      2D
!                        W = wet deposition      3D
!            IF TP==EMS: S = surface emissions   2D
!                        A = aircraft emissions  3D
!                        L = lightning emissions 3D
!                        V = volcanic emissions  3D
!                                            (only available if have SO2 tracer)
!                        S = 3D SO2   emissions  3D
!                                            (only available if have SO2 tracer)
!            IF TP==NET  X = not used
!            IF TP==STE  X = not used
!            IF TP==MAS  X = not used
!            IF TP==PSC  1/2 3D
! RXN#     = If there are multiple reactions with the same reactants AND
!            products then this number specifies which of those (in order) you
!            wish to budget.
!            Set to 0 otherwise.
! #SPECIES = number of species. IF TP==RXN this is total of reactants + products
!                               IF TP!=RXN this should be 1 ONLY
! R1       = IF TP!=RXN this is the species being deposited/emitted
!
!
! NOTE: If the same stash number is specified for multiple diagnostics then
!       these will be summed up before outputting into stash for processing.
!       Specifing the same reaction multiple times is allowed, but only works
!       for diagnostics of the same dimensionallity, i.e. 2D *or* 3D
!
! UNITS OF DIAGNOSTICS ARE: moles/gridcell/s
!
! ######################################################################

MODULE asad_flux_dat

USE missing_data_mod,      ONLY: imdi

IMPLICIT NONE

PRIVATE

! Describes the bimolecular reaction rates
TYPE :: asad_flux_defn
  CHARACTER(LEN=3)  :: diag_type         ! which type of flux:RXN,DEP,EMS
  INTEGER           :: stash_number      ! stash number, e.g. 50001 etc.
  CHARACTER(LEN=1)  :: rxn_type          ! which rxn type: B,H,T,J,D,W,S,A,L,V
  LOGICAL           :: tropospheric_mask ! T or F
  INTEGER           :: rxn_location      ! for multiple reacions with the same
                                         ! reactants and products, default=0
  INTEGER           :: num_species       ! number of reactants+products
  CHARACTER(LEN=10) :: reactants(2)      ! list of reactants
  CHARACTER(LEN=10) :: products(4)       ! list of products
END TYPE asad_flux_defn
PUBLIC asad_flux_defn

! Number of chemical fluxes defined in asad_load_default_fluxes
INTEGER :: n_chemical_fluxes = 0

! STASH section for ASAD diagnostics. This is hard-wired in the flux
! definitions below (as the first 2 digits of stash_number). It is expected
! to match the equivalent UM STASH section number to allow UM legacy style
! requests for ASAD diagnostics to be recognised when UKCA is coupled with
! the UM parent model. Otherwise such requests will be ignored.
INTEGER, PARAMETER, PUBLIC :: stashcode_ukca_chem_diag = 50

TYPE(asad_flux_defn), ALLOCATABLE, SAVE, PUBLIC :: asad_chemical_fluxes(:)

TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_prod(20) =              &
! Production of Ox
       [ asad_flux_defn('RXN',50001,'B',.TRUE.,0,4,                            &
       ['HO2       ','NO        '],                                            &
       ['OH        ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50002,'B',.TRUE.,0,5,                              &
       ['MeOO      ','NO        '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          ']),                 &

! NO + RO2 reactions: sum into STASH section 50 item 3
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['EtOO      ','NO        '],                                            &
       ['MeCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['MeCO3     ','NO        '],                                            &
       ['MeOO      ','CO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['n-PrOO    ','NO        '],                                            &
       ['EtCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['i-PrOO    ','NO        '],                                            &
       ['Me2CO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['EtCO3     ','NO        '],                                            &
       ['EtOO      ','CO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['MeCOCH2OO ','NO        '],                                            &
       ['MeCO3     ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['NO        ','ISO2      '],                                            &
       ['NO2       ','MACR      ','HCHO      ','HO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['NO        ','MACRO2    '],                                            &
       ['NO2       ','MeCO3     ','HACET     ','CO        ']),                 &

! OH + inorganic acid reactions: sum into STASH section 50 item 4
       asad_flux_defn('RXN',50004,'B',.TRUE.,0,4,                              &
       ['OH        ','HONO2     '],                                            &
       ['H2O       ','NO3       ','          ','          ']),                 &
       asad_flux_defn('RXN',50004,'B',.TRUE.,0,4,                              &
       ['OH        ','HONO      '],                                            &
       ['H2O       ','NO2       ','          ','          ']),                 &

! OH + organic nitrate reactions: sum into STASH section 50 item 5
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','MeONO2    '],                                            &
       ['HCHO      ','NO2       ','H2O       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','NALD      '],                                            &
       ['HCHO      ','CO        ','NO2       ','          ']),                 &

! Organic nitrate photolysis: sum into STASH section 50 item 6
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['MeONO2    ','PHOTON    '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                              &
       ['NALD      ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','NO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                              &
       ['ISON      ','PHOTON    '],                                            &
       ['NO2       ','MACR      ','HCHO      ','HO2       ']),                 &

! OH + PAN-type reactions: sum into STASH section 50 item 7
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,5,                              &
       ['OH        ','PAN       '],                                            &
       ['HCHO      ','NO2       ','H2O       ','          ']),                 &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,5,                              &
       ['OH        ','PPAN      '],                                            &
       ['MeCHO     ','NO2       ','H2O       ','          ']),                 &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,4,                              &
       ['OH        ','MPAN      '],                                            &
       ['HACET     ','NO2       ','          ','          '])                  &
       ]


! Different version of OX production for CRI chemistry
! Split into two to avoid continuation errors
TYPE(asad_flux_defn), PARAMETER :: cri_trop_ox_budget_prod01(57) =             &
! Production of Ox
       [ asad_flux_defn('RXN',50001,'B',.TRUE.,0,4,                            &
       ['HO2       ','NO        '],                                            &
       ['OH        ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50002,'B',.TRUE.,0,5,                              &
       ['MeOO      ','NO        '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          ']),                 &

! NO + RO2 reactions: sum into STASH section 50 item 3
! Total = 55
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['EtOO      ','NO        '],                                            &
       ['MeCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['MeCO3     ','NO        '],                                            &
       ['MeOO      ','CO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['i-PrOO    ','NO        '],                                            &
       ['Me2CO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['EtCO3     ','NO        '],                                            &
       ['EtOO      ','CO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN10O2    ','NO        '],                                            &
       ['EtCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN13O2    ','NO        '],                                            &
       ['MeCHO     ','EtOO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN13O2    ','NO        '],                                            &
       ['CARB11A   ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN16O2    ','NO        '],                                            &
       ['RN15AO2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN19O2    ','NO        '],                                            &
       ['RN18AO2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN13AO2   ','NO        '],                                            &
       ['RN12O2    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN16AO2   ','NO        '],                                            &
       ['RN15O2    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RA13O2    ','NO        '],                                            &
       ['CARB3     ','UDCARB8   ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RA16O2    ','NO        '],                                            &
       ['CARB3     ','UDCARB11  ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RA16O2    ','NO        '],                                            &
       ['CARB6     ','UDCARB8   ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RA19AO2   ','NO        '],                                            &
       ['CARB3     ','UDCARB14  ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RA19CO2   ','NO        '],                                            &
       ['CARB9     ','UDCARB8   ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['HOCH2CH2O2','NO        '],                                            &
       ['HCHO      ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['HOCH2CH2O2','NO        '],                                            &
       ['HOCH2CHO  ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RN9O2     ','NO        '],                                            &
       ['MeCHO     ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RN12O2    ','NO        '],                                            &
       ['MeCHO     ','MeCHO     ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RN15O2    ','NO        '],                                            &
       ['EtCHO     ','MeCHO     ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RN18O2    ','NO        '],                                            &
       ['EtCHO     ','EtCHO     ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN15AO2   ','NO        '],                                            &
       ['CARB13    ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN18AO2   ','NO        '],                                            &
       ['CARB16    ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['HOCH2CO3  ','NO        '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN8O2     ','NO        '],                                            &
       ['MeCO3     ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN11O2    ','NO        '],                                            &
       ['MeCO3     ','MeCHO     ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN14O2    ','NO        '],                                            &
       ['EtCO3     ','MeCHO     ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN17O2    ','NO        '],                                            &
       ['RN16AO2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RU14O2    ','NO        '],                                            &
       ['UCARB12   ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RU14O2    ','NO        '],                                            &
       ['UCARB10   ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RU12O2    ','NO        '],                                            &
       ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RU12O2    ','NO        '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RU10O2    ','NO        '],                                            &
       ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RU10O2    ','NO        '],                                            &
       ['CARB6     ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RU10O2    ','NO        '],                                            &
       ['CARB7     ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['NRN6O2    ','NO        '],                                            &
       ['HCHO      ','HCHO      ','NO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['NRN9O2    ','NO        '],                                            &
       ['MeCHO     ','HCHO      ','NO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['NRN12O2   ','NO        '],                                            &
       ['MeCHO     ','MeCHO     ','NO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['NRU14O2   ','NO        '],                                            &
       ['NUCARB12  ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTN28O2   ','NO        '],                                            &
       ['TNCARB26  ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTN28O2   ','NO        '],                                            &
       ['Me2CO     ','RN19O2    ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['NRTN28O2  ','NO        '],                                            &
       ['TNCARB26  ','NO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RTN26O2   ','NO        '],                                            &
       ['RTN25O2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RTN25O2   ','NO        '],                                            &
       ['RTN24O2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RTN24O2   ','NO        '],                                            &
       ['RTN23O2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTN23O2   ','NO        '],                                            &
       ['Me2CO     ','RTN14O2   ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RTN14O2   ','NO        '],                                            &
       ['HCHO      ','TNCARB10  ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO        '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RTX28O2   ','NO        '],                                            &
       ['TXCARB24  ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTX28O2   ','NO        '],                                            &
       ['Me2CO     ','RN19O2    ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['NRTX28O2  ','NO        '],                                            &
       ['TXCARB24  ','HCHO      ','NO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTX24O2   ','NO        '],                                            &
       ['TXCARB22  ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RTX24O2   ','NO        '],                                            &
       ['Me2CO     ','RN13AO2   ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTX22O2   ','NO        '],                                            &
       ['Me2CO     ','RN13O2    ','NO2       ','          '])                  &
]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_trop_ox_budget_prod01(64) =      &
! Production of Ox
       [ asad_flux_defn('RXN',50001,'B',.TRUE.,0,4,                            &
       ['HO2       ','NO        '],                                            &
       ['OH        ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50002,'B',.TRUE.,0,5,                              &
       ['MeOO      ','NO        '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['EtOO      ','NO        '],                                            &
       ['MeCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['MeCO3     ','NO        '],                                            &
       ['MeOO      ','CO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['i-PrOO    ','NO        '],                                            &
       ['Me2CO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['EtCO3     ','NO        '],                                            &
       ['EtOO      ','CO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN10O2    ','NO        '],                                            &
       ['EtCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN13O2    ','NO        '],                                            &
       ['MeCHO     ','EtOO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN13O2    ','NO        '],                                            &
       ['CARB11A   ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN16O2    ','NO        '],                                            &
       ['RN15AO2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN19O2    ','NO        '],                                            &
       ['RN18AO2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN13AO2   ','NO        '],                                            &
       ['RN12O2    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN16AO2   ','NO        '],                                            &
       ['RN15O2    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN13AO2   ','NO        '],                                            &
       ['RN12O2    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN16AO2   ','NO        '],                                            &
       ['RN15O2    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RA13O2    ','NO        '],                                            &
       ['CARB3     ','UDCARB8   ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RA16O2    ','NO        '],                                            &
       ['CARB3     ','UDCARB11  ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RA16O2    ','NO        '],                                            &
       ['CARB6     ','UDCARB8   ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RA19AO2   ','NO        '],                                            &
       ['CARB3     ','UDCARB14  ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RA19CO2   ','NO        '],                                            &
       ['CARB9     ','UDCARB8   ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['HOCH2CH2O2','NO        '],                                            &
       ['HCHO      ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['HOCH2CH2O2','NO        '],                                            &
       ['HOCH2CHO  ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RN9O2     ','NO        '],                                            &
       ['MeCHO     ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RN12O2    ','NO        '],                                            &
       ['MeCHO     ','MeCHO     ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RN15O2    ','NO        '],                                            &
       ['EtCHO     ','MeCHO     ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RN18O2    ','NO        '],                                            &
       ['EtCHO     ','EtCHO     ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN15AO2   ','NO        '],                                            &
       ['CARB13    ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN18AO2   ','NO        '],                                            &
       ['CARB16    ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['HOCH2CO3  ','NO        '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN8O2     ','NO        '],                                            &
       ['MeCO3     ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN11O2    ','NO        '],                                            &
       ['MeCO3     ','MeCHO     ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RN14O2    ','NO        '],                                            &
       ['EtCO3     ','MeCHO     ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RN17O2    ','NO        '],                                            &
       ['RN16AO2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RU14O2    ','NO        '],                                            &
       ['UCARB12   ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RU14O2    ','NO        '],                                            &
       ['UCARB10   ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RU12O2    ','NO        '],                                            &
       ['CARB6     ','HOCH2CHO  ','NO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RU12O2    ','NO        '],                                            &
       ['CARB7     ','CARB3     ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RU12O2    ','NO        '],                                            &
       ['CARB7     ','HOCH2CO3  ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RU10O2    ','NO        '],                                            &
       ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RU10O2    ','NO        '],                                            &
       ['CARB6     ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['NRN6O2    ','NO        '],                                            &
       ['HCHO      ','HCHO      ','NO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['NRN9O2    ','NO        '],                                            &
       ['MeCHO     ','HCHO      ','NO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['NRN12O2   ','NO        '],                                            &
       ['MeCHO     ','MeCHO     ','NO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['NRU14O2   ','NO        '],                                            &
       ['NUCARB12  ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTN28O2   ','NO        '],                                            &
       ['TNCARB26  ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTN28O2   ','NO        '],                                            &
       ['Me2CO     ','RN19O2    ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['NRTN28O2  ','NO        '],                                            &
       ['TNCARB26  ','NO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RTN26O2   ','NO        '],                                            &
       ['RTN25O2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RTN25O2   ','NO        '],                                            &
       ['RTN24O2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,4,                              &
       ['RTN24O2   ','NO        '],                                            &
       ['RTN23O2   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTN23O2   ','NO        '],                                            &
       ['Me2CO     ','RTN14O2   ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RTN14O2   ','NO        '],                                            &
       ['HCHO      ','TNCARB10  ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO        '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RTX28O2   ','NO        '],                                            &
       ['TXCARB24  ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTX28O2   ','NO        '],                                            &
       ['Me2CO     ','RN19O2    ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['NRTX28O2  ','NO        '],                                            &
       ['TXCARB24  ','HCHO      ','NO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTX24O2   ','NO        '],                                            &
       ['TXCARB22  ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RTX24O2   ','NO        '],                                            &
       ['Me2CO     ','RN13AO2   ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                              &
       ['RTX22O2   ','NO        '],                                            &
       ['Me2CO     ','RN13O2    ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['RU10AO2   ','NO        '],                                            &
       ['CARB7     ','HO2       ','CO        ','NO2       ']),                 &
       ! MACO3 new RO2 species
       ! Below 2 entries are single reaction split over 2 chem_master entries
       ! due to 5 products. First part added twice
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO        '],                                            &
       ['MeOO      ','CO        ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO        '],                                            &
       ['MeOO      ','CO        ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO        '],                                            &
       ['MeCO3     ','HCHO      ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                              &
       ['DHPR12O2  ','NO        '],                                            &
       ['CARB3     ','RN8OOH    ','OH        ','NO2       '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: cri_trop_ox_budget_prod02(41) =             &
! OH + inorganic acid reactions: sum into STASH section 50 item 4
       [ asad_flux_defn('RXN',50004,'B',.TRUE.,0,4,                            &
       ['OH        ','HONO2     '],                                            &
       ['H2O       ','NO3       ','          ','          ']),                 &
       asad_flux_defn('RXN',50004,'B',.TRUE.,0,4,                              &
       ['OH        ','HONO      '],                                            &
       ['H2O       ','NO2       ','          ','          ']),                 &
! OH + organic nitrate reactions: sum into STASH section 50 item 5
! Including CRI reactions, Total = 22
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','MeONO2    '],                                            &
       ['HCHO      ','NO2       ','H2O       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','EtONO2    '],                                            &
       ['MeCHO     ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN10NO3   '],                                            &
       ['EtCHO     ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','i-PrONO2  '],                                            &
       ['Me2CO     ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN13NO3   '],                                            &
       ['CARB11A   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN16NO3   '],                                            &
       ['CARB14    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN19NO3   '],                                            &
       ['CARB17    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','HOC2H4NO3 '],                                            &
       ['HOCH2CHO  ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN9NO3    '],                                            &
       ['CARB7     ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN12NO3   '],                                            &
       ['CARB10    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN15NO3   '],                                            &
       ['CARB13    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN18NO3   '],                                            &
       ['CARB16    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RU14NO3   '],                                            &
       ['UCARB12   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RA13NO3   '],                                            &
       ['CARB3     ','UDCARB8   ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RA16NO3   '],                                            &
       ['CARB3     ','UDCARB11  ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RA19NO3   '],                                            &
       ['CARB6     ','UDCARB11  ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RTN28NO3  '],                                            &
       ['TNCARB26  ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RTN25NO3  '],                                            &
       ['Me2CO     ','TNCARB15  ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RTX28NO3  '],                                            &
       ['TXCARB24  ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RTX24NO3  '],                                            &
       ['TXCARB22  ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RTX22NO3  '],                                            &
       ['Me2CO     ','CCARB12   ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['RTN23NO3  ','OH        '],                                            &
       ['Me2CO     ','TNCARB12  ','NO2       ','          ']),                 &

! Organic nitrate photolysis: sum into STASH section 50 item 6
! Total = 12
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['MeONO2    ','PHOTON    '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['EtONO2    ','PHOTON    '],                                            &
       ['MeCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['RN10NO3   ','PHOTON    '],                                            &
       ['EtCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['i-PrONO2  ','PHOTON    '],                                            &
       ['Me2CO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['RN13NO3   ','PHOTON    '],                                            &
       ['MeCHO     ','EtOO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['RN13NO3   ','PHOTON    '],                                            &
       ['CARB11A   ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,4,                              &
       ['RN16NO3   ','PHOTON    '],                                            &
       ['RN15O2    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,4,                              &
       ['RN19NO3   ','PHOTON    '],                                            &
       ['RN18O2    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                              &
       ['RA13NO3   ','PHOTON    '],                                            &
       ['CARB3     ','UDCARB8   ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                              &
       ['RA16NO3   ','PHOTON    '],                                            &
       ['CARB3     ','UDCARB11  ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                              &
       ['RA19NO3   ','PHOTON    '],                                            &
       ['CARB6     ','UDCARB11  ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['RTX24NO3  ','PHOTON    '],                                            &
       ['TXCARB22  ','HO2       ','NO2       ','          ']),                 &
! OH + PAN-type reactions:
! These do not actually produce net Ox, so should not be included,
! but are included here to compare like-for-like with StratTrop
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,6,                              &
       ['OH        ','PAN       '],                                            &
       ['HCHO      ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,6,                              &
       ['OH        ','PPAN      '],                                            &
       ['MeCHO     ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,5,                              &
       ['OH        ','MPAN      '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,4,                              &
       ['OH        ','RU12PAN   '],                                            &
       ['UCARB10   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,5,                              &
       ['OH        ','RTN26PAN  '],                                            &
       ['Me2CO     ','CARB16    ','NO2       ','          '])                  &
        ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_trop_ox_budget_prod02(44) =      &
! OH + inorganic acid reactions: sum into STASH section 50 item 4
       [ asad_flux_defn('RXN',50004,'B',.TRUE.,0,4,                            &
       ['OH        ','HONO2     '],                                            &
       ['H2O       ','NO3       ','          ','          ']),                 &
       asad_flux_defn('RXN',50004,'B',.TRUE.,0,4,                              &
       ['OH        ','HONO      '],                                            &
       ['H2O       ','NO2       ','          ','          ']),                 &
! OH + organic nitrate reactions: sum into STASH section 50 item 5
! Including CRI reactions, Total = 22 + 2 new = 24
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','MeONO2    '],                                            &
       ['HCHO      ','NO2       ','H2O       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','EtONO2    '],                                            &
       ['MeCHO     ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN10NO3   '],                                            &
       ['EtCHO     ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','i-PrONO2  '],                                            &
       ['Me2CO     ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN13NO3   '],                                            &
       ['CARB11A   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN16NO3   '],                                            &
       ['CARB14    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN19NO3   '],                                            &
       ['CARB17    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','HOC2H4NO3 '],                                            &
       ['HOCH2CHO  ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN9NO3    '],                                            &
       ['CARB7     ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN12NO3   '],                                            &
       ['CARB10    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN15NO3   '],                                            &
       ['CARB13    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RN18NO3   '],                                            &
       ['CARB16    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RU14NO3   '],                                            &
       ['UCARB12   ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RA13NO3   '],                                            &
       ['CARB3     ','UDCARB8   ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RA16NO3   '],                                            &
       ['CARB3     ','UDCARB11  ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RA19NO3   '],                                            &
       ['CARB6     ','UDCARB11  ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RTN28NO3  '],                                            &
       ['TNCARB26  ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RTN25NO3  '],                                            &
       ['Me2CO     ','TNCARB15  ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RTX28NO3  '],                                            &
       ['TXCARB24  ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,4,                              &
       ['OH        ','RTX24NO3  '],                                            &
       ['TXCARB22  ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['OH        ','RTX22NO3  '],                                            &
       ['Me2CO     ','CCARB12   ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['RU12NO3   ','OH        '],                                            &
       ['CARB7     ','CARB3     ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                              &
       ['RU10NO3   ','OH        '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['MeONO2    ','PHOTON    '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['EtONO2    ','PHOTON    '],                                            &
       ['MeCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['RN10NO3   ','PHOTON    '],                                            &
       ['EtCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['i-PrONO2  ','PHOTON    '],                                            &
       ['Me2CO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['RN13NO3   ','PHOTON    '],                                            &
       ['MeCHO     ','EtOO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['RN13NO3   ','PHOTON    '],                                            &
       ['CARB11A   ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,4,                              &
       ['RN16NO3   ','PHOTON    '],                                            &
       ['RN15O2    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,4,                              &
       ['RN19NO3   ','PHOTON    '],                                            &
       ['RN18O2    ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                              &
       ['RA13NO3   ','PHOTON    '],                                            &
       ['CARB3     ','UDCARB8   ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                              &
       ['RA16NO3   ','PHOTON    '],                                            &
       ['CARB3     ','UDCARB11  ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                              &
       ['RA19NO3   ','PHOTON    '],                                            &
       ['CARB6     ','UDCARB11  ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['RTX24NO3  ','PHOTON    '],                                            &
       ['TXCARB22  ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                              &
       ['RU12NO3   ','PHOTON    '],                                            &
       ['CARB6     ','HOCH2CHO  ','NO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                              &
       ['RU10NO3   ','PHOTON    '],                                            &
       ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                              &
       ['RU14NO3   ','PHOTON    '],                                            &
       ['UCARB10   ','HCHO      ','HO2       ','NO2       ']),                 &
! Skipping OH + PAN-type reactions:
! These do not produce net Ox, so should not be included
!++SAN - readding these just so that I can compare like-for-like with
! StratTrop (+5)
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,6,                              &
       ['OH        ','PAN       '],                                            &
       ['HCHO      ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,6,                              &
       ['OH        ','PPAN      '],                                            &
       ['MeCHO     ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,5,                              &
       ['OH        ','MPAN      '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,5,                              &
       ['OH        ','RTN26PAN  '],                                            &
       ['Me2CO     ','CARB16    ','NO2       ','          '])                  &
       ]



TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_loss01(23)              &
! Loss of Ox
       = [ asad_flux_defn('RXN',50011,'B',.TRUE.,0,4,                          &
       ['O(1D)     ','H2O       '],                                            &
       ['OH        ','OH        ','          ','          ']),                 &

! Minor reactions, should have negligble impact:
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['OH        ','MeOO      ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','H2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,5,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','HO2       ','HO2       ','          ']),                 &

! Include with above. Each is loss of 2xOx, so include twice to sum
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50013,'B',.TRUE.,0,4,                              &
       ['HO2       ','O3        '],                                            &
       ['OH        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50014,'B',.TRUE.,0,4,                              &
       ['OH        ','O3        '],                                            &
       ['HO2       ','O2        ','          ','          ']),                 &

! O3 + alkene:
! O3 + isoprene reactions
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MACR      ','HCHO      ','MACRO2    ','MeCO3     ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeOO      ','HCOOH     ','CO        ','H2O2      ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C5H8      '],                                            &
       ['HO2       ','OH        ','          ','          ']),                 &

! O3 + MACR reactions. ratb_defs specifies 2 different rates for each of these,
!  so need to select each one in turn
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,4,                              &
       ['O3        ','MACR      '],                                            &
       ['OH        ','MeCO3     ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,2,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,2,4,                              &
       ['O3        ','MACR      '],                                            &
       ['OH        ','MeCO3     ','          ','          ']),                 &

! N2O5 + H20 reaction
       asad_flux_defn('RXN',50016,'B',.TRUE.,0,4,                              &
       ['N2O5      ','H2O       '],                                            &
       ['HONO2     ','HONO2     ','          ','          ']),                 &
! Heterogenious reactions are not included yet
!RXN    50016     H       4       N2O5    H2O    HONO2    HONO2

! NO3 chemical loss
! sink of 2xOx:
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
! these are sinks of 1xOx
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['HO2       ','NO3       '],                                            &
       ['OH        ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['OH        ','NO3       '],                                            &
       ['HO2       ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['MeOO      ','NO3       '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER, PUBLIC :: asad_trop_ox_budget_loss01_132(23)  &
! Version for use when i_ukca_chem_version >= 132
! Loss of Ox
       = [ asad_flux_defn('RXN',50011,'B',.TRUE.,0,4,                          &
       ['O(1D)     ','H2O       '],                                            &
       ['OH        ','OH        ','          ','          ']),                 &

! Minor reactions, should have negligble impact:
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['OH        ','MeOO      ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','H2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,5,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','HO2       ','HO2       ','          ']),                 &

! Include with above. Each is loss of 2xOx, so include twice to sum
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50013,'B',.TRUE.,0,4,                              &
       ['HO2       ','O3        '],                                            &
       ['OH        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50014,'B',.TRUE.,0,4,                              &
       ['OH        ','O3        '],                                            &
       ['HO2       ','O2        ','          ','          ']),                 &

! O3 + alkene:
! O3 + isoprene reactions
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MACR      ','HCHO      ','MACRO2    ','MeCO3     ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeOO      ','HCOOH     ','CO        ','H2O2      ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['O3        ','C5H8      '],                                            &
       ['HO2       ','OH        ','SEC_ORG_I ','          ']),                 &

! O3 + MACR reactions. ratb_defs specifies 2 different rates for each of these,
!  so need to select each one in turn
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,4,                              &
       ['O3        ','MACR      '],                                            &
       ['OH        ','MeCO3     ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,2,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,2,4,                              &
       ['O3        ','MACR      '],                                            &
       ['OH        ','MeCO3     ','          ','          ']),                 &
! N2O5 + H20 reaction
       asad_flux_defn('RXN',50016,'B',.TRUE.,0,4,                              &
       ['N2O5      ','H2O       '],                                            &
       ['HONO2     ','HONO2     ','          ','          ']),                 &
! Heterogenious reactions are not included yet
!RXN    50016     H       4       N2O5    H2O    HONO2    HONO2

! NO3 chemical loss
! sink of 2xOx:
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
! these are sinks of 1xOx
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['HO2       ','NO3       '],                                            &
       ['OH        ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['OH        ','NO3       '],                                            &
       ['HO2       ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['MeOO      ','NO3       '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          '])                  &
       ]



! Whole atmosphere CH4 rxn-flux diagnostics
TYPE(asad_flux_defn), PARAMETER :: asad_atmos_ch4_budget_loss(6)               &
  ! Whole atmosphere CH4+OH rxn-flux
       = [ asad_flux_defn('RXN',50428,'B',.FALSE.,0,4,                         &
       ['OH        ','CH4       '],                                            &
       ['H2O       ','MeOO      ','          ','          ']),                 &
  ! Whole atmosphere CH4+O(1D) rxn-fluxs
       asad_flux_defn('RXN',50429,'B',.FALSE.,0,4,                             &
       ['O(1D)     ','CH4       '],                                            &
       ['OH        ','MeOO      ','          ','          ']),                 &
       asad_flux_defn('RXN',50429,'B',.FALSE.,0,4,                             &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','H2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50429,'B',.FALSE.,0,5,                             &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','HO2       ','HO2       ','          ']),                 &
  ! Whole atmosphere CH4+Cl rxn-flux
       asad_flux_defn('RXN',50430,'B',.FALSE.,0,4,                             &
       ['Cl        ','CH4       '],                                            &
       ['HCl       ','MeOO      ','          ','          ']),                 &
  ! Whole atmosphere CH4+hv rxn-flux
       asad_flux_defn('RXN',50431,'J',.FALSE.,0,4,                             &
       ['CH4       ','PHOTON    '],                                            &
       ['MeOO      ','H         ','          ','          '])                  &
       ]

! Methane dry deposition (3D) written into s50i439
! Dry deposition:
TYPE(asad_flux_defn), PARAMETER :: asad_ch4_drydep(1)                          &
       = [ asad_flux_defn('DEP',50438,'D',.FALSE.,0,1,                         &
       ['CH4       ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]

! Methane strat-trop exchange flux (3D) written into s50i438
! Strat-trop xchng:
TYPE(asad_flux_defn), PARAMETER :: asad_ch4_ste(1)                             &
       = [ asad_flux_defn('STE',50439,'X',.FALSE.,0,1,                         &
       ['CH4       ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]

! Equivalent ox budget losses for CRI mechanism
TYPE(asad_flux_defn), PARAMETER :: cri_trop_ox_budget_loss01(34)               &
! Loss of Ox
       = [ asad_flux_defn('RXN',50011,'B',.TRUE.,0,4,                          &
       ['O(1D)     ','H2O       '],                                            &
       ['OH        ','OH        ','          ','          ']),                 &
! Minor reactions, should have negligble impact:
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['OH        ','MeOO      ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','H2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,5,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','HO2       ','HO2       ','          ']),                 &

! Include with above. Each is loss of 2xOx, so include twice to sum
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50013,'B',.TRUE.,0,4,                              &
       ['HO2       ','O3        '],                                            &
       ['OH        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50014,'B',.TRUE.,0,4,                              &
       ['OH        ','O3        '],                                            &
       ['HO2       ','O2        ','          ','          ']),                 &
! O3 + alkene:
! From here is where CRI diverges from StratTrop. N = 18
! O3 + isoprene reactions
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','HCOOH     ','          ','          ']),                 &
! O3 + other alkenes
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','HCOOH     ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','MeCO2H    ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','MeCO2H    ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['APINENE   ','O3        '],                                            &
       ['OH        ','Me2CO     ','RN18AO2   ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['APINENE   ','O3        '],                                            &
       ['TNCARB26  ','H2O2      ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,3,                              &
       ['APINENE   ','O3        '],                                            &
       ['RCOOH25   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['BPINENE   ','O3        '],                                            &
       ['RTX24O2   ','OH        ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB24  ','H2O2      ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB22  ','          ','          ']),                 &
! O3 + unsaturated aldehyde reactions (O3 + MACR and similar)
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,5,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','CARB6     ','H2O2      ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,5,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','CARB6     ','H2O2      ','          ']),                 &
! N2O5 + H20 reaction - added to CRI mechanism
       asad_flux_defn('RXN',50016,'B',.TRUE.,0,4,                              &
       ['N2O5      ','H2O       '],                                            &
       ['HONO2     ','HONO2     ','          ','          ']),                 &
! Heterogenious reactions are not included yet
!RXN    50016     H       4       N2O5    H2O    HONO2    HONO2

! NO3 chemical loss
! sink of 2xOx:
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
! these are sinks of 1xOx
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['HO2       ','NO3       '],                                            &
       ['OH        ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['OH        ','NO3       '],                                            &
       ['HO2       ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['MeOO      ','NO3       '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_trop_ox_budget_loss01(41)        &
! Loss of Ox
       = [ asad_flux_defn('RXN',50011,'B',.TRUE.,0,4,                          &
       ['O(1D)     ','H2O       '],                                            &
       ['OH        ','OH        ','          ','          ']),                 &
! Minor reactions, should have negligble impact:
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['OH        ','MeOO      ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','H2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,5,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','HO2       ','HO2       ','          ']),                 &
! Include with above. Each is loss of 2xOx, so include twice to sum
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50013,'B',.TRUE.,0,4,                              &
       ['HO2       ','O3        '],                                            &
       ['OH        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50014,'B',.TRUE.,0,4,                              &
       ['OH        ','O3        '],                                            &
       ['HO2       ','O2        ','          ','          ']),                 &
! O3 + alkene:
! From here is where CRI diverges from StratTrop. N = 18
! O3 + isoprene reactions
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','HCOOH     ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeCO3     ','HCHO      ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['HCHO      ','MeOO      ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','HCHO      ','H2O2      ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','          ','          ']),                 &
! O3 + other alkenes
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','HCOOH     ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','MeCO2H    ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','MeCO2H    ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['APINENE   ','O3        '],                                            &
       ['OH        ','Me2CO     ','RN18AO2   ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['APINENE   ','O3        '],                                            &
       ['TNCARB26  ','H2O2      ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,3,                              &
       ['APINENE   ','O3        '],                                            &
       ['RCOOH25   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['BPINENE   ','O3        '],                                            &
       ['RTX24O2   ','OH        ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB24  ','H2O2      ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB22  ','          ','          ']),                 &
! O3 + unsaturated aldehyde reactions (O3 + MACR and similar)
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,4,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','CARB6     ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','MeCO3     ','CARB3     ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB6     ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,5,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB3     ','CARB7     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','HPUCARB12 '],                                            &
       ['CARB3     ','CARB6     ','OH        ','OH        ']),                 &
! N2O5 + H20 reaction - added to CRI mechanism
       asad_flux_defn('RXN',50016,'B',.TRUE.,0,4,                              &
       ['N2O5      ','H2O       '],                                            &
       ['HONO2     ','HONO2     ','          ','          ']),                 &
! Heterogenious reactions are not included yet
       asad_flux_defn('RXN',50016,'H',.TRUE.,0,4,                              &
       ['N2O5      ','H2O       '],                                            &
       ['HONO2     ','HONO2     ','          ','          ']),                 &
! NO3 chemical loss
! sink of 2xOx:
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
! these are sinks of 1xOx
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['HO2       ','NO3       '],                                            &
       ['OH        ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['OH        ','NO3       '],                                            &
       ['HO2       ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['MeOO      ','NO3       '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          '])                  &
       ]



! Equivalent ox budget losses for CRI mechanism with GLOMAP aerosols
! (Different Monoterp+Ox type reactions)
TYPE(asad_flux_defn), PARAMETER :: cri_aer_trop_ox_budget_loss01(34)           &
! Loss of Ox
       = [ asad_flux_defn('RXN',50011,'B',.TRUE.,0,4,                          &
       ['O(1D)     ','H2O       '],                                            &
       ['OH        ','OH        ','          ','          ']),                 &
! Minor reactions, should have negligble impact:
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['OH        ','MeOO      ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','H2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,5,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','HO2       ','HO2       ','          ']),                 &

! Include with above. Each is loss of 2xOx, so include twice to sum
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50013,'B',.TRUE.,0,4,                              &
       ['HO2       ','O3        '],                                            &
       ['OH        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50014,'B',.TRUE.,0,4,                              &
       ['OH        ','O3        '],                                            &
       ['HO2       ','O2        ','          ','          ']),                 &
! O3 + alkenes:
! Here is where CRI diverges from StratTrop. N = 18
! O3 + isoprene reactions
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','HCOOH     ','          ','          ']),                 &
! O3 + other alkenes
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','HCOOH     ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','MeCO2H    ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','MeCO2H    ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['APINENE   ','O3        '],                                            &
       ['OH        ','Me2CO     ','RN18AO2   ','Sec_Org   ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['APINENE   ','O3        '],                                            &
       ['TNCARB26  ','H2O2      ','Sec_Org   ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['APINENE   ','O3        '],                                            &
       ['RCOOH25   ','Sec_Org   ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['BPINENE   ','O3        '],                                            &
       ['RTX24O2   ','OH        ','Sec_Org   ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB24  ','H2O2      ','Sec_Org   ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB22  ','Sec_Org   ','          ']),                 &
! O3 + unsaturated aldehyde reactions (O3 + MACR and similar)
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,5,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','CARB6     ','H2O2      ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,5,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','CARB6     ','H2O2      ','          ']),                 &
! N2O5 + H20 reaction - added to CRI mechanism
       asad_flux_defn('RXN',50016,'B',.TRUE.,0,4,                              &
       ['N2O5      ','H2O       '],                                            &
       ['HONO2     ','HONO2     ','          ','          ']),                 &
! Heterogenious reactions are not included yet
!RXN    50016     H       4       N2O5    H2O    HONO2    HONO2

! NO3 chemical loss
! sink of 2xOx:
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
! these are sinks of 1xOx
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['HO2       ','NO3       '],                                            &
       ['OH        ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['OH        ','NO3       '],                                            &
       ['HO2       ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['MeOO      ','NO3       '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_aer_trop_ox_budget_loss01(41)    &
! Loss of Ox
       = [ asad_flux_defn('RXN',50011,'B',.TRUE.,0,4,                          &
       ['O(1D)     ','H2O       '],                                            &
       ['OH        ','OH        ','          ','          ']),                 &
! Minor reactions, should have negligble impact:
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['OH        ','MeOO      ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','H2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,5,                              &
       ['O(1D)     ','CH4       '],                                            &
       ['HCHO      ','HO2       ','HO2       ','          ']),                 &

! Include with above. Each is loss of 2xOx, so include twice to sum
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','O3        '],                                            &
       ['O2        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                              &
       ['O(3P)     ','NO2       '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50013,'B',.TRUE.,0,4,                              &
       ['HO2       ','O3        '],                                            &
       ['OH        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50014,'B',.TRUE.,0,4,                              &
       ['OH        ','O3        '],                                            &
       ['HO2       ','O2        ','          ','          ']),                 &
! O3 + alkenes:
! Here is where CRI diverges from StratTrop. N = 18
! O3 + isoprene reactions
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','HCOOH     ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeCO3     ','HCHO      ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['HCHO      ','MeOO      ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','HCHO      ','H2O2      ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','          ','          ']),                 &
! O3 + other alkenes
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','HCOOH     ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','MeCO2H    ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','MeCO2H    ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['APINENE   ','O3        '],                                            &
       ['OH        ','Me2CO     ','RN18AO2   ','Sec_Org   ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['APINENE   ','O3        '],                                            &
       ['TNCARB26  ','H2O2      ','Sec_Org   ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                              &
       ['APINENE   ','O3        '],                                            &
       ['RCOOH25   ','Sec_Org   ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,5,                              &
       ['BPINENE   ','O3        '],                                            &
       ['RTX24O2   ','OH        ','Sec_Org   ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB24  ','H2O2      ','Sec_Org   ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB22  ','Sec_Org   ','          ']),                 &
! O3 + unsaturated aldehyde reactions (O3 + MACR and similar)
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,4,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','CARB6     ','          ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','MeCO3     ','CARB3     ','OH        ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB6     ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,5,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB3     ','CARB7     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                              &
       ['O3        ','HPUCARB12 '],                                            &
       ['CARB3     ','CARB6     ','OH        ','OH        ']),                 &
! N2O5 + H20 reaction - added to CRI mechanism
       asad_flux_defn('RXN',50016,'B',.TRUE.,0,4,                              &
       ['N2O5      ','H2O       '],                                            &
       ['HONO2     ','HONO2     ','          ','          ']),                 &
! Heterogenious reactions are not included yet
       asad_flux_defn('RXN',50016,'H',.TRUE.,0,4,                              &
       ['N2O5      ','H2O       '],                                            &
       ['HONO2     ','HONO2     ','          ','          ']),                 &
! NO3 chemical loss
! sink of 2xOx:
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                              &
       ['NO3       ','PHOTON    '],                                            &
       ['NO        ','O2        ','          ','          ']),                 &
! these are sinks of 1xOx
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['HO2       ','NO3       '],                                            &
       ['OH        ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['OH        ','NO3       '],                                            &
       ['HO2       ','NO2       ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['MeOO      ','NO3       '],                                            &
       ['HO2       ','HCHO      ','NO2       ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_loss02(12)              &
     = [ asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                            &
       ['EtOO      ','NO3       '],                                            &
       ['MeCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['MeCO3     ','NO3       '],                                            &
       ['MeOO      ','CO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['n-PrOO    ','NO3       '],                                            &
       ['EtCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['i-PrOO    ','NO3       '],                                            &
       ['Me2CO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['EtCO3     ','NO3       '],                                            &
       ['EtOO      ','CO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['MeCOCH2OO ','NO3       '],                                            &
       ['MeCO3     ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['NO3       ','MeCHO     '],                                            &
       ['HONO2     ','MeCO3     ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['NO3       ','EtCHO     '],                                            &
       ['HONO2     ','EtCO3     ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['NO3       ','Me2CO     '],                                            &
       ['HONO2     ','MeCOCH2OO ','          ','          ']),                 &
      ! sink of 2xOx
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                              &
       ['NO3       ','C5H8      '],                                            &
       ['ISON      ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                              &
       ['NO3       ','C5H8      '],                                            &
       ['ISON      ','          ','          ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER, PUBLIC :: asad_trop_ox_budget_loss02_132(12)  &
! Version for use when i_ukca_chem_version >= 132
     = [ asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                            &
       ['EtOO      ','NO3       '],                                            &
       ['MeCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['MeCO3     ','NO3       '],                                            &
       ['MeOO      ','CO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['n-PrOO    ','NO3       '],                                            &
       ['EtCHO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['i-PrOO    ','NO3       '],                                            &
       ['Me2CO     ','HO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['EtCO3     ','NO3       '],                                            &
       ['EtOO      ','CO2       ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['MeCOCH2OO ','NO3       '],                                            &
       ['MeCO3     ','HCHO      ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['NO3       ','MeCHO     '],                                            &
       ['HONO2     ','MeCO3     ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['NO3       ','EtCHO     '],                                            &
       ['HONO2     ','EtCO3     ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['NO3       ','Me2CO     '],                                            &
       ['HONO2     ','MeCOCH2OO ','          ','          ']),                 &
      ! sink of 2xOx
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['NO3       ','C5H8      '],                                            &
       ['ISON      ','SEC_ORG_I ','          ','          ']),                 &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                              &
       ['NO3       ','C5H8      '],                                            &
       ['ISON      ','SEC_ORG_I ','          ','          '])                  &
       ]


! CRI version with O3 loss via NO3+org reactions
TYPE(asad_flux_defn), PARAMETER :: cri_trop_ox_budget_loss02(72)               &
      = [ asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                           &
        ['EtOO      ','NO3       '],                                           &
        ['MeCHO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['MeCO3     ','NO3       '],                                           &
        ['MeOO      ','CO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['i-PrOO    ','NO3       '],                                           &
        ['Me2CO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['EtCO3     ','NO3       '],                                           &
        ['EtOO      ','CO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','HCHO      '],                                           &
        ['HONO2     ','HO2       ','CO        ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','MeCHO     '],                                           &
        ['HONO2     ','MeCO3     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','EtCHO     '],                                           &
        ['HONO2     ','EtCO3     ','          ','          ']),                &
        ! sink of 1 x Ox => Considering NRU14O2 as a reservoir
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C5H8      '],                                           &
        ['NRU14O2   ','          ','          ','          ']),                &
        ! From here, new reactions in CRI mechanism
        ! NO3 + Alkene reactions => Net -1*Ox loss (5)
        !  NO3 counts for 2xOx, N-RO2 species as reservoir for one Ox
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C2H4      '],                                           &
        ['NRN6O2    ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C3H6      '],                                           &
        ['NRN9O2    ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','TBUT2ENE  '],                                           &
        ['NRN12O2   ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['APINENE   ','NO3       '],                                           &
        ['NRTN28O2  ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['BPINENE   ','NO3       '],                                           &
        ['NRTX28O2  ','          ','          ','          ']),                &
        ! NO3 + RO2 reactions - 1 Ox each (50)
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN10O2    ','NO3       '],                                           &
        ['EtCHO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN13O2    ','NO3       '],                                           &
        ['MeCHO     ','EtOO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN13O2    ','NO3       '],                                           &
        ['CARB11A   ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16O2    ','NO3       '],                                           &
        ['RN15AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN19O2    ','NO3       '],                                           &
        ['RN18AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN13AO2   ','NO3       '],                                           &
        ['RN12O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16AO2   ','NO3       '],                                           &
        ['RN15O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16AO2   ','NO3       '],                                           &
        ['RN15O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA13O2    ','NO3       '],                                           &
        ['CARB3     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA16O2    ','NO3       '],                                           &
        ['CARB3     ','UDCARB11  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA16O2    ','NO3       '],                                           &
        ['CARB6     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA19AO2   ','NO3       '],                                           &
        ['CARB3     ','UDCARB14  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA19CO2   ','NO3       '],                                           &
        ['CARB9     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['HOCH2CH2O2','NO3       '],                                           &
        ['HCHO      ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['HOCH2CH2O2','NO3       '],                                           &
        ['HOCH2CHO  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN9O2     ','NO3       '],                                           &
        ['MeCHO     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN12O2    ','NO3       '],                                           &
        ['MeCHO     ','MeCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN15O2    ','NO3       '],                                           &
        ['EtCHO     ','MeCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN18O2    ','NO3       '],                                           &
        ['EtCHO     ','EtCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN15AO2   ','NO3       '],                                           &
        ['CARB13    ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN18AO2   ','NO3       '],                                           &
        ['CARB16    ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['HOCH2CO3  ','NO3       '],                                           &
        ['HO2       ','HCHO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN8O2     ','NO3       '],                                           &
        ['MeCO3     ','HCHO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN11O2    ','NO3       '],                                           &
        ['MeCO3     ','MeCHO     ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN14O2    ','NO3       '],                                           &
        ['EtCO3     ','MeCHO     ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN17O2    ','NO3       '],                                           &
        ['RN16AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RU14O2    ','NO3       '],                                           &
        ['UCARB12   ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU14O2    ','NO3       '],                                           &
        ['UCARB10   ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RU12O2    ','NO3       '],                                           &
        ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU12O2    ','NO3       '],                                           &
        ['CARB7     ','CO        ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RU10O2    ','NO3       '],                                           &
        ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU10O2    ','NO3       '],                                           &
        ['CARB6     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU10O2    ','NO3       '],                                           &
        ['CARB7     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN6O2    ','NO3       '],                                           &
        ['HCHO      ','HCHO      ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN9O2    ','NO3       '],                                           &
        ['MeCHO     ','HCHO      ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN12O2   ','NO3       '],                                           &
        ['MeCHO     ','MeCHO     ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NRU14O2   ','NO3       '],                                           &
        ['NUCARB12  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRU12O2   ','NO3       '],                                           &
        ['NOA       ','CO        ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN28O2   ','NO3       '],                                           &
        ['TNCARB26  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NRTN28O2  ','NO3       '],                                           &
        ['TNCARB26  ','NO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN26O2   ','NO3       '],                                           &
        ['RTN25O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN25O2   ','NO3       '],                                           &
        ['RTN24O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN24O2   ','NO3       '],                                           &
        ['RTN23O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN23O2   ','NO3       '],                                           &
        ['Me2CO     ','RTN14O2   ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RTN14O2   ','NO3       '],                                           &
        ['HCHO      ','TNCARB10  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN10O2   ','NO3       '],                                           &
        ['RN8O2     ','CO        ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RTX28O2   ','NO3       '],                                           &
        ['TXCARB24  ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTX24O2   ','NO3       '],                                           &
        ['TXCARB22  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTX22O2   ','NO3       '],                                           &
        ['Me2CO     ','RN13O2    ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRTX28O2  ','NO3       '],                                           &
        ['TXCARB24  ','HCHO      ','NO2       ','NO2       ']),                &
        ! NO3 + carbonyls (7)
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','UCARB10   '],                                           &
        ['RU10O2    ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','HOCH2CHO  '],                                           &
        ['HOCH2CO3  ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','UCARB12   '],                                           &
        ['RU12O2    ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','TNCARB26  '],                                           &
        ['RTN26O2   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','TNCARB10  '],                                           &
        ['RTN10O2   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','AROH14    '],                                           &
        ['RAROH14   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','AROH14    '],                                           &
        ['RAROH14   ','HONO2     ','          ','          ']),                &
        ! Nitrophenols: Net -1 Ox loss when react with NO3 (2)
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','ARNOH14   '],                                           &
        ['CARB13    ','NO2       ','HONO2     ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','ARNOH17   '],                                           &
        ['CARB16    ','NO2       ','HONO2     ','          '])                 &
        ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_trop_ox_budget_loss02(77)        &
      = [ asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                           &
        ['EtOO      ','NO3       '],                                           &
        ['MeCHO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['MeCO3     ','NO3       '],                                           &
        ['MeOO      ','CO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['i-PrOO    ','NO3       '],                                           &
        ['Me2CO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['EtCO3     ','NO3       '],                                           &
        ['EtOO      ','CO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','HCHO      '],                                           &
        ['HONO2     ','HO2       ','CO        ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','MeCHO     '],                                           &
        ['HONO2     ','MeCO3     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','EtCHO     '],                                           &
        ['HONO2     ','EtCO3     ','          ','          ']),                &
        ! sink of 1 x Ox => Considering NRU14O2 as a reservoir
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C5H8      '],                                           &
        ['NRU14O2   ','          ','          ','          ']),                &
        ! From here, new reactions in CRI mechanism
        ! NO3 + Alkene reactions => Net -1*Ox loss (5)
        !  NO3 counts for 2xOx, N-RO2 species as reservoir for one Ox
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C2H4      '],                                           &
        ['NRN6O2    ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C3H6      '],                                           &
        ['NRN9O2    ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','TBUT2ENE  '],                                           &
        ['NRN12O2   ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['APINENE   ','NO3       '],                                           &
        ['NRTN28O2  ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['BPINENE   ','NO3       '],                                           &
        ['NRTX28O2  ','          ','          ','          ']),                &
        ! NO3 + RO2 reactions - 1 Ox each (50)
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN10O2    ','NO3       '],                                           &
        ['EtCHO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN13O2    ','NO3       '],                                           &
        ['MeCHO     ','EtOO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN13O2    ','NO3       '],                                           &
        ['CARB11A   ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16O2    ','NO3       '],                                           &
        ['RN15AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN19O2    ','NO3       '],                                           &
        ['RN18AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN13AO2   ','NO3       '],                                           &
        ['RN12O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16AO2   ','NO3       '],                                           &
        ['RN15O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16AO2   ','NO3       '],                                           &
        ['RN15O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA13O2    ','NO3       '],                                           &
        ['CARB3     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA16O2    ','NO3       '],                                           &
        ['CARB3     ','UDCARB11  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA16O2    ','NO3       '],                                           &
        ['CARB6     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA19AO2   ','NO3       '],                                           &
        ['CARB3     ','UDCARB14  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA19CO2   ','NO3       '],                                           &
        ['CARB9     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['HOCH2CH2O2','NO3       '],                                           &
        ['HCHO      ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['HOCH2CH2O2','NO3       '],                                           &
        ['HOCH2CHO  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN9O2     ','NO3       '],                                           &
        ['MeCHO     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN12O2    ','NO3       '],                                           &
        ['MeCHO     ','MeCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN15O2    ','NO3       '],                                           &
        ['EtCHO     ','MeCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN18O2    ','NO3       '],                                           &
        ['EtCHO     ','EtCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN15AO2   ','NO3       '],                                           &
        ['CARB13    ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN18AO2   ','NO3       '],                                           &
        ['CARB16    ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['HOCH2CO3  ','NO3       '],                                           &
        ['HO2       ','HCHO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN8O2     ','NO3       '],                                           &
        ['MeCO3     ','HCHO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN11O2    ','NO3       '],                                           &
        ['MeCO3     ','MeCHO     ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN14O2    ','NO3       '],                                           &
        ['EtCO3     ','MeCHO     ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN17O2    ','NO3       '],                                           &
        ['RN16AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RU14O2    ','NO3       '],                                           &
        ['UCARB12   ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU14O2    ','NO3       '],                                           &
        ['UCARB10   ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU12O2    ','NO3       '],                                           &
        ['CARB6     ','HOCH2CHO  ','NO2       ','HO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU12O2    ','NO3       '],                                           &
        ['CARB7     ','CARB3     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RU10O2    ','NO3       '],                                           &
        ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU10O2    ','NO3       '],                                           &
        ['CARB6     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN6O2    ','NO3       '],                                           &
        ['HCHO      ','HCHO      ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN9O2    ','NO3       '],                                           &
        ['MeCHO     ','HCHO      ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN12O2   ','NO3       '],                                           &
        ['MeCHO     ','MeCHO     ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NRU14O2   ','NO3       '],                                           &
        ['NUCARB12  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRU12O2   ','NO3       '],                                           &
        ['NOA       ','CO        ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN28O2   ','NO3       '],                                           &
        ['TNCARB26  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NRTN28O2  ','NO3       '],                                           &
        ['TNCARB26  ','NO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN26O2   ','NO3       '],                                           &
        ['RTN25O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN25O2   ','NO3       '],                                           &
        ['RTN24O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN24O2   ','NO3       '],                                           &
        ['RTN23O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN23O2   ','NO3       '],                                           &
        ['Me2CO     ','RTN14O2   ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RTN14O2   ','NO3       '],                                           &
        ['HCHO      ','TNCARB10  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN10O2   ','NO3       '],                                           &
        ['RN8O2     ','CO        ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RTX28O2   ','NO3       '],                                           &
        ['TXCARB24  ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTX24O2   ','NO3       '],                                           &
        ['TXCARB22  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTX22O2   ','NO3       '],                                           &
        ['Me2CO     ','RN13O2    ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRTX28O2  ','NO3       '],                                           &
        ['TXCARB24  ','HCHO      ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['DHPR12O2  ','NO3       '],                                           &
        ['CARB3     ','RN8OOH    ','OH        ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['MACO3     ','NO3       '],                                           &
        ['MeOO      ','HCHO      ','HO2       ','CO        ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['MACO3     ','NO3       '],                                           &
        ['MeOO      ','HCHO      ','HO2       ','CO        ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['MACO3     ','NO3       '],                                           &
        ['MeCO3     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU10AO2   ','NO3       '],                                           &
        ['CARB7     ','CO        ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NO3       ','HPUCARB12 '],                                           &
        ['HUCARB9   ','CO        ','OH        ','HONO2     ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','UCARB10   '],                                           &
        ['RU10O2    ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','HOCH2CHO  '],                                           &
        ['HOCH2CO3  ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','UCARB12   '],                                           &
        ['RU12O2    ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','TNCARB26  '],                                           &
        ['RTN26O2   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','TNCARB10  '],                                           &
        ['RTN10O2   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','AROH14    '],                                           &
        ['RAROH14   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','AROH14    '],                                           &
        ['RAROH14   ','HONO2     ','          ','          ']),                &
        ! Nitrophenols: Net -1 Ox loss when react with NO3 (2)
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','ARNOH14   '],                                           &
        ['CARB13    ','NO2       ','HONO2     ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','ARNOH17   '],                                           &
        ['CARB16    ','NO2       ','HONO2     ','          '])                 &
        ]





! CRI version with O3 loss via NO3+org reactions
! Different versions of monoterpene+NO3 reactions when running with aerosol
TYPE(asad_flux_defn), PARAMETER :: cri_aer_trop_ox_budget_loss02(72)           &
      = [ asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                           &
        ['EtOO      ','NO3       '],                                           &
        ['MeCHO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['MeCO3     ','NO3       '],                                           &
        ['MeOO      ','CO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['i-PrOO    ','NO3       '],                                           &
        ['Me2CO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['EtCO3     ','NO3       '],                                           &
        ['EtOO      ','CO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','HCHO      '],                                           &
        ['HONO2     ','HO2       ','CO        ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','MeCHO     '],                                           &
        ['HONO2     ','MeCO3     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','EtCHO     '],                                           &
        ['HONO2     ','EtCO3     ','          ','          ']),                &
        ! sink of 1 x Ox => Considering NRU14O2 as a reservoir
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C5H8      '],                                           &
        ['NRU14O2   ','          ','          ','          ']),                &
        ! From here, new reactions in CRI mechanism
        ! NO3 + Alkene reactions => Net -1*Ox loss (5)
        !  NO3 counts for 2xOx, N-RO2 species as reservoir for one Ox
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C2H4      '],                                           &
        ['NRN6O2    ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C3H6      '],                                           &
        ['NRN9O2    ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','TBUT2ENE  '],                                           &
        ['NRN12O2   ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['APINENE   ','NO3       '],                                           &
        ['NRTN28O2  ','Sec_Org   ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['BPINENE   ','NO3       '],                                           &
        ['NRTX28O2  ','Sec_Org   ','          ','          ']),                &
        ! NO3 + RO2 reactions - 1 Ox each (50)
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN10O2    ','NO3       '],                                           &
        ['EtCHO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN13O2    ','NO3       '],                                           &
        ['MeCHO     ','EtOO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN13O2    ','NO3       '],                                           &
        ['CARB11A   ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16O2    ','NO3       '],                                           &
        ['RN15AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN19O2    ','NO3       '],                                           &
        ['RN18AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN13AO2   ','NO3       '],                                           &
        ['RN12O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16AO2   ','NO3       '],                                           &
        ['RN15O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16AO2   ','NO3       '],                                           &
        ['RN15O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA13O2    ','NO3       '],                                           &
        ['CARB3     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA16O2    ','NO3       '],                                           &
        ['CARB3     ','UDCARB11  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA16O2    ','NO3       '],                                           &
        ['CARB6     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA19AO2   ','NO3       '],                                           &
        ['CARB3     ','UDCARB14  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA19CO2   ','NO3       '],                                           &
        ['CARB9     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['HOCH2CH2O2','NO3       '],                                           &
        ['HCHO      ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['HOCH2CH2O2','NO3       '],                                           &
        ['HOCH2CHO  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN9O2     ','NO3       '],                                           &
        ['MeCHO     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN12O2    ','NO3       '],                                           &
        ['MeCHO     ','MeCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN15O2    ','NO3       '],                                           &
        ['EtCHO     ','MeCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN18O2    ','NO3       '],                                           &
        ['EtCHO     ','EtCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN15AO2   ','NO3       '],                                           &
        ['CARB13    ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN18AO2   ','NO3       '],                                           &
        ['CARB16    ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['HOCH2CO3  ','NO3       '],                                           &
        ['HO2       ','HCHO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN8O2     ','NO3       '],                                           &
        ['MeCO3     ','HCHO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN11O2    ','NO3       '],                                           &
        ['MeCO3     ','MeCHO     ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN14O2    ','NO3       '],                                           &
        ['EtCO3     ','MeCHO     ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN17O2    ','NO3       '],                                           &
        ['RN16AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RU14O2    ','NO3       '],                                           &
        ['UCARB12   ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU14O2    ','NO3       '],                                           &
        ['UCARB10   ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RU12O2    ','NO3       '],                                           &
        ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU12O2    ','NO3       '],                                           &
        ['CARB7     ','CO        ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RU10O2    ','NO3       '],                                           &
        ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU10O2    ','NO3       '],                                           &
        ['CARB6     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU10O2    ','NO3       '],                                           &
        ['CARB7     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN6O2    ','NO3       '],                                           &
        ['HCHO      ','HCHO      ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN9O2    ','NO3       '],                                           &
        ['MeCHO     ','HCHO      ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN12O2   ','NO3       '],                                           &
        ['MeCHO     ','MeCHO     ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NRU14O2   ','NO3       '],                                           &
        ['NUCARB12  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRU12O2   ','NO3       '],                                           &
        ['NOA       ','CO        ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN28O2   ','NO3       '],                                           &
        ['TNCARB26  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NRTN28O2  ','NO3       '],                                           &
        ['TNCARB26  ','NO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN26O2   ','NO3       '],                                           &
        ['RTN25O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN25O2   ','NO3       '],                                           &
        ['RTN24O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN24O2   ','NO3       '],                                           &
        ['RTN23O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN23O2   ','NO3       '],                                           &
        ['Me2CO     ','RTN14O2   ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RTN14O2   ','NO3       '],                                           &
        ['HCHO      ','TNCARB10  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN10O2   ','NO3       '],                                           &
        ['RN8O2     ','CO        ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RTX28O2   ','NO3       '],                                           &
        ['TXCARB24  ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTX24O2   ','NO3       '],                                           &
        ['TXCARB22  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTX22O2   ','NO3       '],                                           &
        ['Me2CO     ','RN13O2    ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRTX28O2  ','NO3       '],                                           &
        ['TXCARB24  ','HCHO      ','NO2       ','NO2       ']),                &
        ! NO3 + carbonyls (7)
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','UCARB10   '],                                           &
        ['RU10O2    ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','HOCH2CHO  '],                                           &
        ['HOCH2CO3  ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','UCARB12   '],                                           &
        ['RU12O2    ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','TNCARB26  '],                                           &
        ['RTN26O2   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','TNCARB10  '],                                           &
        ['RTN10O2   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','AROH14    '],                                           &
        ['RAROH14   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','AROH14    '],                                           &
        ['RAROH14   ','HONO2     ','          ','          ']),                &
        ! Nitrophenols: Net -1 Ox loss when react with NO3 (2)
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','ARNOH14   '],                                           &
        ['CARB13    ','NO2       ','HONO2     ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','ARNOH17   '],                                           &
        ['CARB16    ','NO2       ','HONO2     ','          '])                 &
        ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_aer_trop_ox_budget_loss02(77)    &
      = [ asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                           &
        ['EtOO      ','NO3       '],                                           &
        ['MeCHO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['MeCO3     ','NO3       '],                                           &
        ['MeOO      ','CO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['i-PrOO    ','NO3       '],                                           &
        ['Me2CO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['EtCO3     ','NO3       '],                                           &
        ['EtOO      ','CO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','HCHO      '],                                           &
        ['HONO2     ','HO2       ','CO        ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','MeCHO     '],                                           &
        ['HONO2     ','MeCO3     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','EtCHO     '],                                           &
        ['HONO2     ','EtCO3     ','          ','          ']),                &
        ! sink of 1 x Ox => Considering NRU14O2 as a reservoir
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C5H8      '],                                           &
        ['NRU14O2   ','          ','          ','          ']),                &
        ! From here, new reactions in CRI mechanism
        ! NO3 + Alkene reactions => Net -1*Ox loss (5)
        !  NO3 counts for 2xOx, N-RO2 species as reservoir for one Ox
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C2H4      '],                                           &
        ['NRN6O2    ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','C3H6      '],                                           &
        ['NRN9O2    ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['NO3       ','TBUT2ENE  '],                                           &
        ['NRN12O2   ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['APINENE   ','NO3       '],                                           &
        ['NRTN28O2  ','Sec_Org   ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['BPINENE   ','NO3       '],                                           &
        ['NRTX28O2  ','Sec_Org   ','          ','          ']),                &
        ! NO3 + RO2 reactions - 1 Ox each (50)
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN10O2    ','NO3       '],                                           &
        ['EtCHO     ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN13O2    ','NO3       '],                                           &
        ['MeCHO     ','EtOO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN13O2    ','NO3       '],                                           &
        ['CARB11A   ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16O2    ','NO3       '],                                           &
        ['RN15AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN19O2    ','NO3       '],                                           &
        ['RN18AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN13AO2   ','NO3       '],                                           &
        ['RN12O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16AO2   ','NO3       '],                                           &
        ['RN15O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN16AO2   ','NO3       '],                                           &
        ['RN15O2    ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA13O2    ','NO3       '],                                           &
        ['CARB3     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA16O2    ','NO3       '],                                           &
        ['CARB3     ','UDCARB11  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA16O2    ','NO3       '],                                           &
        ['CARB6     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA19AO2   ','NO3       '],                                           &
        ['CARB3     ','UDCARB14  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RA19CO2   ','NO3       '],                                           &
        ['CARB9     ','UDCARB8   ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['HOCH2CH2O2','NO3       '],                                           &
        ['HCHO      ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['HOCH2CH2O2','NO3       '],                                           &
        ['HOCH2CHO  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN9O2     ','NO3       '],                                           &
        ['MeCHO     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN12O2    ','NO3       '],                                           &
        ['MeCHO     ','MeCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN15O2    ','NO3       '],                                           &
        ['EtCHO     ','MeCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RN18O2    ','NO3       '],                                           &
        ['EtCHO     ','EtCHO     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN15AO2   ','NO3       '],                                           &
        ['CARB13    ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN18AO2   ','NO3       '],                                           &
        ['CARB16    ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['HOCH2CO3  ','NO3       '],                                           &
        ['HO2       ','HCHO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN8O2     ','NO3       '],                                           &
        ['MeCO3     ','HCHO      ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN11O2    ','NO3       '],                                           &
        ['MeCO3     ','MeCHO     ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RN14O2    ','NO3       '],                                           &
        ['EtCO3     ','MeCHO     ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RN17O2    ','NO3       '],                                           &
        ['RN16AO2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RU14O2    ','NO3       '],                                           &
        ['UCARB12   ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU14O2    ','NO3       '],                                           &
        ['UCARB10   ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU12O2    ','NO3       '],                                           &
        ['CARB6     ','HOCH2CHO  ','NO2       ','HO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU12O2    ','NO3       '],                                           &
        ['CARB7     ','CARB3     ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RU10O2    ','NO3       '],                                           &
        ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU10O2    ','NO3       '],                                           &
        ['CARB6     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN6O2    ','NO3       '],                                           &
        ['HCHO      ','HCHO      ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN9O2    ','NO3       '],                                           &
        ['MeCHO     ','HCHO      ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRN12O2   ','NO3       '],                                           &
        ['MeCHO     ','MeCHO     ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NRU14O2   ','NO3       '],                                           &
        ['NUCARB12  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRU12O2   ','NO3       '],                                           &
        ['NOA       ','CO        ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN28O2   ','NO3       '],                                           &
        ['TNCARB26  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NRTN28O2  ','NO3       '],                                           &
        ['TNCARB26  ','NO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN26O2   ','NO3       '],                                           &
        ['RTN25O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN25O2   ','NO3       '],                                           &
        ['RTN24O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['RTN24O2   ','NO3       '],                                           &
        ['RTN23O2   ','NO2       ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN23O2   ','NO3       '],                                           &
        ['Me2CO     ','RTN14O2   ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RTN14O2   ','NO3       '],                                           &
        ['HCHO      ','TNCARB10  ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTN10O2   ','NO3       '],                                           &
        ['RN8O2     ','CO        ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RTX28O2   ','NO3       '],                                           &
        ['TXCARB24  ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTX24O2   ','NO3       '],                                           &
        ['TXCARB22  ','HO2       ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['RTX22O2   ','NO3       '],                                           &
        ['Me2CO     ','RN13O2    ','NO2       ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NRTX28O2  ','NO3       '],                                           &
        ['TXCARB24  ','HCHO      ','NO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['DHPR12O2  ','NO3       '],                                           &
        ['CARB3     ','RN8OOH    ','OH        ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['MACO3     ','NO3       '],                                           &
        ['MeOO      ','HCHO      ','HO2       ','CO        ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                             &
        ['MACO3     ','NO3       '],                                           &
        ['NO2       ','          ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['MACO3     ','NO3       '],                                           &
        ['MeCO3     ','HCHO      ','HO2       ','NO2       ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['RU10AO2   ','NO3       '],                                           &
        ['CARB7     ','CO        ','HO2       ','NO2       ']),                &
        ! NO3 + carbonyls
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,6,                             &
        ['NO3       ','HPUCARB12 '],                                           &
        ['HUCARB9   ','CO        ','OH        ','HONO2     ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','UCARB10   '],                                           &
        ['RU10O2    ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','HOCH2CHO  '],                                           &
        ['HOCH2CO3  ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','UCARB12   '],                                           &
        ['RU12O2    ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','TNCARB26  '],                                           &
        ['RTN26O2   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','TNCARB10  '],                                           &
        ['RTN10O2   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','AROH14    '],                                           &
        ['RAROH14   ','HONO2     ','          ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                             &
        ['NO3       ','AROH14    '],                                           &
        ['RAROH14   ','HONO2     ','          ','          ']),                &
        ! Nitrophenols: Net -1 Ox loss when react with NO3 (2)
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','ARNOH14   '],                                           &
        ['CARB13    ','NO2       ','HONO2     ','          ']),                &
        asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                             &
        ['NO3       ','ARNOH17   '],                                           &
        ['CARB16    ','NO2       ','HONO2     ','          '])                 &
        ]



TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_drydep(11)              &
! Dry deposition:
! O3 dry deposition
       = [asad_flux_defn('DEP',50021,'D',.TRUE.,0,1,                           &
       ['O3        ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! NOy dry deposition
! sink of 2xOx
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! sink on 3xOx
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['PAN       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['PPAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['MPAN      ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


! Trop Ox DD loss in CRI
TYPE(asad_flux_defn), PARAMETER :: cri_trop_ox_budget_drydep(26)               &
! Dry deposition:
! O3 dry deposition
       = [asad_flux_defn('DEP',50021,'D',.TRUE.,0,1,                           &
       ['O3        ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! NOy dry deposition
! sink of 2xOx
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! sink of 3xOx
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['PAN       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['PPAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['MPAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       ! Include other PAN-type species in CRI
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['PHAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['RU12PAN   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['RTN26PAN  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['RTN26PAN  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       ! Also include nitro-phenols in Ox budget
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['ARNOH14   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['ARNOH17   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       ! Also include products of NO3 + RH reactions
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NOA       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NUCARB12  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRU14OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRU12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRN6OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRN9OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRN12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRTN28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRTX28OOH ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_trop_ox_budget_drydep(24)        &
! Dry deposition:
! O3 dry deposition
       = [asad_flux_defn('DEP',50021,'D',.TRUE.,0,1,                           &
       ['O3        ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! NOy dry deposition
! sink of 2xOx
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! sink of 3xOx
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['PAN       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['PPAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['MPAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       ! Include other PAN-type species in CRI
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['PHAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['RTN26PAN  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       ! Also include nitro-phenols in Ox budget
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['ARNOH14   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['ARNOH17   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       ! Also include products of NO3 + RH reactions
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NOA       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NUCARB12  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRU14OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRU12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRN6OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRN9OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRN12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRTN28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                              &
       ['NRTX28OOH ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]



TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_wetdep(7)=              &
! Wet deposition:
! NOy wet deposition
! sink of 2xOx
       [ asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                            &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! sink of 3xOx
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER :: cri_trop_ox_budget_wetdep(18)=              &
! Wet deposition:
! O3 should wet deposit when running with aerosol (small)
! NOy wet deposition
! sink of 2xOx
       [ asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                            &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! sink of 3xOx
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! Include other NOy species in CRI
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['ARNOH14   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['ARNOH17   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NOA       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NUCARB12  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRU14OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRU12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRN6OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRN9OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRN12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRTN28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRTX28OOH ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_trop_ox_budget_wetdep(20)=       &
! Wet deposition:
! O3 should wet deposit when running with aerosol (small)
! NOy wet deposition
! sink of 2xOx
       [ asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                            &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! sink of 3xOx
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! Include other NOy species in CRI
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['ARNOH14   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['ARNOH17   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NOA       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NUCARB12  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRU14OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRU12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRN6OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRN9OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRN12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRTN28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['NRTX28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['RTN26PAN  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                              &
       ['PHAN      ','          '],                                            &
       ['          ','          ','          ','          '])                  &
]


! Total oxidised nitrogen 3D WET deposition into 50-241

TYPE(asad_flux_defn), PARAMETER :: asad_oxidN_wetdep(9)=                       &
! Wet deposition of NO3, N2O5, HONO3, HONO2, ISON, BrONO2, ClONO2,
! and HONO added together
       [ asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                           &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! 2 N per N2O5
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['ISON      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['BrONO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['ClONO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['HONO      ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: cri_oxidN_wetdep(19)=                       &
! Wet deposition of NO3, N2O5, HONO3, HONO2, ISON, BrONO2, ClONO2,
! and HONO added together
       [ asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                           &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! 2 N per N2O5
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['BrONO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['ClONO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['HONO      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       ! CRI organonitrates that wet deposit
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NOA       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &

       ['NUCARB12  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRU14OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRU12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRN6OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRN9OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRN12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRTN28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRTX28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['ARNOH14   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['ARNOH17   ','          '],                                            &

       ['          ','          ','          ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_oxidN_wetdep(36)=                &
! Wet deposition of NO3, N2O5, HONO3, HONO2, ISON, BrONO2, ClONO2,
! and HONO added together
       [ asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                           &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! 2 N per N2O5
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['BrONO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['ClONO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['HONO      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       ! CRI organonitrates that wet deposit
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NOA       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NUCARB12  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRU14OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRU12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRN6OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRN9OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRN12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRTN28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['NRTX28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['ARNOH14   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['ARNOH17   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RU14NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RU12NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RU10NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RA13NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RA16NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RA19NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RTN26PAN  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RN9NO3    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RN12NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RN15NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RN18NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RTN28NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RTN25NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RTX28NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['RTX22NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['PHAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50241,'W',.FALSE.,0,1,                             &
       ['HOC2H4NO3 ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]



! Total oxidised nitrogen 3D DRY deposition into 50-242

TYPE(asad_flux_defn), PARAMETER :: asad_oxidN_drydep(13)                       &
! Dry deposition: NO, NO2, NO3, HONO, HONO2, HONO3, N2O5, PAN, PPAN,
! MPAN, ISON, and NALD
       = [asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                          &
       ['NO        ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NO2       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['HONO      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! 2 N per N2O5
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['PAN       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['PPAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['MPAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['ISON      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NALD      ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: cri_oxidN_drydep(37)                        &
! Dry deposition: NO, NO2, NO3, HONO, HONO2, HONO3, N2O5, PAN, PPAN,
! MPAN, ISON, and NALD
       = [asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                          &
       ['NO        ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NO2       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['HONO      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! 2 N per N2O5
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['PAN       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['PPAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['MPAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NOA       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['HOC2H4NO3 ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['PHAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NUCARB12  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTX24NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRU14OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRU12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRN6OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRN9OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRN12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRTN28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRTX28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RU12PAN   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTN26PAN  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTN23NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RN9NO3    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RN12NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RN15NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RN18NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RU14NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTN28NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTN25NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTX28NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTX22NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['ARNOH14   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['ARNOH17   ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_oxidN_drydep(40)                 &
! Dry deposition: NO, NO2, NO3, HONO, HONO2, HONO3, N2O5, PAN, PPAN,
! MPAN, ISON, and NALD
       = [asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                          &
       ['NO        ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NO2       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NO3       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['HONO      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['HONO2     ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['HO2NO2    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! 2 N per N2O5
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['N2O5      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['PAN       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['PPAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['MPAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NOA       ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['HOC2H4NO3 ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['PHAN      ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NUCARB12  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTX24NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRU14OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRU12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRN6OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRN9OOH   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRN12OOH  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRTN28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['NRTX28OOH ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTN26PAN  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RN9NO3    ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RN12NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RN15NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RN18NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RU14NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTN28NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTN25NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTX28NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RTX22NO3  ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['ARNOH14   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['ARNOH17   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RU12NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RU10NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RA13NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RA16NO3   ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('DEP',50242,'D',.FALSE.,0,1,                             &
       ['RA19NO3   ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: asad_trop_other_fluxes(16) =                &
! Extra fluxes of interest
       [ asad_flux_defn('RXN',50042,'B',.TRUE.,0,3,                            &
       ['NO3       ','C5H8      '],                                            &
       ['ISON      ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50043,'B',.TRUE.,0,3,                              &
       ['NO        ','ISO2      '],                                            &
       ['ISON      ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'T',.TRUE.,0,4,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,3,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeOO      '],                                            &
       ['MeOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','EtOO      '],                                            &
       ['EtOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeCO3     '],                                            &
       ['MeCO3H    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                              &
       ['HO2       ','MeCO3     '],                                            &
       ['MeCO2H    ','O3        ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','n-PrOO    '],                                            &
       ['n-PrOOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','i-PrOO    '],                                            &
       ['i-PrOOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                              &
       ['HO2       ','EtCO3     '],                                            &
       ['O2        ','EtCO3H    ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                              &
       ['HO2       ','EtCO3     '],                                            &
       ['EtCO2H    ','O3        ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeCOCH2OO '],                                            &
       ['MeCOCH2OOH','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','ISO2      '],                                            &
       ['ISOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MACRO2    '],                                            &
       ['MACROOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50046,'T',.TRUE.,0,4,                              &
       ['OH        ','NO2       '],                                            &
       ['HONO2     ','m         ','          ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER, PUBLIC :: asad_trop_other_fluxes_132(16) =    &
! Extra fluxes of interest when i_ukca_chem_version >= 132
       [ asad_flux_defn('RXN',50042,'B',.TRUE.,0,4,                            &
       ['NO3       ','C5H8      '],                                            &
       ['ISON      ','SEC_ORG_I ','          ','          ']),                 &
       asad_flux_defn('RXN',50043,'B',.TRUE.,0,3,                              &
       ['NO        ','ISO2      '],                                            &
       ['ISON      ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'T',.TRUE.,0,4,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,3,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeOO      '],                                            &
       ['MeOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','EtOO      '],                                            &
       ['EtOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeCO3     '],                                            &
       ['MeCO3H    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                              &
       ['HO2       ','MeCO3     '],                                            &
       ['MeCO2H    ','O3        ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','n-PrOO    '],                                            &
       ['n-PrOOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','i-PrOO    '],                                            &
       ['i-PrOOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                              &
       ['HO2       ','EtCO3     '],                                            &
       ['O2        ','EtCO3H    ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                              &
       ['HO2       ','EtCO3     '],                                            &
       ['EtCO2H    ','O3        ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeCOCH2OO '],                                            &
       ['MeCOCH2OOH','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','ISO2      '],                                            &
       ['ISOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MACRO2    '],                                            &
       ['MACROOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50046,'T',.TRUE.,0,4,                              &
       ['OH        ','NO2       '],                                            &
       ['HONO2     ','m         ','          ','          '])                  &
       ]

! Equivalent fluxes in CRI mechanism
TYPE(asad_flux_defn), PARAMETER :: cri_trop_other_fluxes(56) =                 &
! Extra fluxes of interest
       [ asad_flux_defn('RXN',50042,'B',.TRUE.,0,3,                            &
       ['NO3       ','C5H8      '],                                            &
       ['NRU14O2   ','          ','          ','          ']),                 &
       ! Equiv. to NO+ISO2
       asad_flux_defn('RXN',50043,'B',.TRUE.,0,3,                              &
       ['NO        ','RU14O2    '],                                            &
       ['RU14NO3   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'T',.TRUE.,0,4,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,3,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','          ','          ','          ']),                 &
! Including production of H2O2 from O3+alkene reactions (4)
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,4,                              &
       ['APINENE   ','O3        '],                                            &
       ['TNCARB26  ','H2O2      ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,5,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB24  ','H2O2      ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,5,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','CARB6     ','H2O2      ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,5,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','CARB6     ','H2O2      ','          ']),                 &
! 50045 - ROOH production (5)
! Note not including acetic acid production - not a hydroperoxide!
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeOO      '],                                            &
       ['MeOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','EtOO      '],                                            &
       ['EtOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeCO3     '],                                            &
       ['MeCO3H    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','i-PrOO    '],                                            &
       ['i-PrOOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                              &
       ['HO2       ','EtCO3     '],                                            &
       ['O2        ','EtCO3H    ','          ','          ']),                 &
! New peroxide production in CRI (43):
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN10O2    ','HO2       '],                                            &
       ['RN10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN13O2    ','HO2       '],                                            &
       ['RN13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN16O2    ','HO2       '],                                            &
       ['RN16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN19O2    ','HO2       '],                                            &
       ['RN19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN13AO2   ','HO2       '],                                            &
       ['RN13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN16AO2   ','HO2       '],                                            &
       ['RN16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA13O2    ','HO2       '],                                            &
       ['RA13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA16O2    ','HO2       '],                                            &
       ['RA16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA19AO2   ','HO2       '],                                            &
       ['RA19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA19CO2   ','HO2       '],                                            &
       ['RA19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HOCH2CH2O2','HO2       '],                                            &
       ['HOC2H4OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN9O2     ','HO2       '],                                            &
       ['RN9OOH    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN12O2    ','HO2       '],                                            &
       ['RN12OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN15O2    ','HO2       '],                                            &
       ['RN15OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN18O2    ','HO2       '],                                            &
       ['RN18OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN15AO2   ','HO2       '],                                            &
       ['RN15OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN18AO2   ','HO2       '],                                            &
       ['RN18OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HOCH2CO3  ','HO2       '],                                            &
       ['HOCH2CO3H ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN8O2     ','HO2       '],                                            &
       ['RN8OOH    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN11O2    ','HO2       '],                                            &
       ['RN11OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN14O2    ','HO2       '],                                            &
       ['RN14OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN17O2    ','HO2       '],                                            &
       ['RN17OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU14O2    ','HO2       '],                                            &
       ['RU14OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU12O2    ','HO2       '],                                            &
       ['RU12OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU10O2    ','HO2       '],                                            &
       ['RU10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN6O2    ','HO2       '],                                            &
       ['NRN6OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN9O2    ','HO2       '],                                            &
       ['NRN9OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN12O2   ','HO2       '],                                            &
       ['NRN12OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRU14O2   ','HO2       '],                                            &
       ['NRU14OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRU12O2   ','HO2       '],                                            &
       ['NRU12OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN28O2   ','HO2       '],                                            &
       ['RTN28OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRTN28O2  ','HO2       '],                                            &
       ['NRTN28OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN26O2   ','HO2       '],                                            &
       ['RTN26OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN25O2   ','HO2       '],                                            &
       ['RTN25OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN24O2   ','HO2       '],                                            &
       ['RTN24OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN23O2   ','HO2       '],                                            &
       ['RTN23OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN14O2   ','HO2       '],                                            &
       ['RTN14OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN10O2   ','HO2       '],                                            &
       ['RTN10OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX28O2   ','HO2       '],                                            &
       ['RTX28OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX24O2   ','HO2       '],                                            &
       ['RTX24OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX22O2   ','HO2       '],                                            &
       ['RTX22OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRTX28O2  ','HO2       '],                                            &
       ['NRTX28OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50046,'T',.TRUE.,0,4,                              &
       ['OH        ','NO2       '],                                            &
       ['HONO2     ','m         ','          ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_trop_other_fluxes(64) =          &
! Extra fluxes of interest
       [ asad_flux_defn('RXN',50042,'B',.TRUE.,0,3,                            &
       ['NO3       ','C5H8      '],                                            &
       ['NRU14O2   ','          ','          ','          ']),                 &
       ! Equiv. to NO+ISO2
       asad_flux_defn('RXN',50043,'B',.TRUE.,0,3,                              &
       ['NO        ','RU14O2    '],                                            &
       ['RU14NO3   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'T',.TRUE.,0,4,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,3,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','          ','          ','          ']),                 &
! Including production of H2O2 from O3+alkene reactions (4)
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,4,                              &
       ['APINENE   ','O3        '],                                            &
       ['TNCARB26  ','H2O2      ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,5,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB24  ','H2O2      ','          ']),                 &
! CS2 O3 + C5H8 now makes H2O2
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,5,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','HCHO      ','H2O2      ','          ']),                 &
! 50045 - ROOH production (5)
! Note not including acetic acid production - not a hydroperoxide!
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeOO      '],                                            &
       ['MeOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','EtOO      '],                                            &
       ['EtOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeCO3     '],                                            &
       ['MeCO3H    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','i-PrOO    '],                                            &
       ['i-PrOOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,5,                              &
       ['HO2       ','EtCO3     '],                                            &
       ['OH        ','EtCO3H    ','EtOO      ','          ']),                 &
! New peroxide production in CRI (43):
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN10O2    ','HO2       '],                                            &
       ['RN10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN13O2    ','HO2       '],                                            &
       ['RN13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN16O2    ','HO2       '],                                            &
       ['RN16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN19O2    ','HO2       '],                                            &
       ['RN19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN13AO2   ','HO2       '],                                            &
       ['RN13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN16AO2   ','HO2       '],                                            &
       ['RN16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA13O2    ','HO2       '],                                            &
       ['RA13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA16O2    ','HO2       '],                                            &
       ['RA16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA19AO2   ','HO2       '],                                            &
       ['RA19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA19CO2   ','HO2       '],                                            &
       ['RA19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HOCH2CH2O2','HO2       '],                                            &
       ['HOC2H4OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN9O2     ','HO2       '],                                            &
       ['RN9OOH    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN12O2    ','HO2       '],                                            &
       ['RN12OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN15O2    ','HO2       '],                                            &
       ['RN15OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN18O2    ','HO2       '],                                            &
       ['RN18OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN15AO2   ','HO2       '],                                            &
       ['RN15OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN18AO2   ','HO2       '],                                            &
       ['RN18OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HOCH2CO3  ','HO2       '],                                            &
       ['HOCH2CO3H ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN8O2     ','HO2       '],                                            &
       ['RN8OOH    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN11O2    ','HO2       '],                                            &
       ['RN11OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN14O2    ','HO2       '],                                            &
       ['RN14OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN17O2    ','HO2       '],                                            &
       ['RN17OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU14O2    ','HO2       '],                                            &
       ['RU14OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU12O2    ','HO2       '],                                            &
       ['RU12OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU10O2    ','HO2       '],                                            &
       ['RU10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN6O2    ','HO2       '],                                            &
       ['NRN6OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN9O2    ','HO2       '],                                            &
       ['NRN9OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN12O2   ','HO2       '],                                            &
       ['NRN12OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRU14O2   ','HO2       '],                                            &
       ['NRU14OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRU12O2   ','HO2       '],                                            &
       ['NRU12OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN28O2   ','HO2       '],                                            &
       ['RTN28OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRTN28O2  ','HO2       '],                                            &
       ['NRTN28OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,5,                              &
       ['RTN26O2   ','HO2       '],                                            &
       ['RTN26OOH  ','RTN25O2   ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN25O2   ','HO2       '],                                            &
       ['RTN25OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN24O2   ','HO2       '],                                            &
       ['RTN24OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN23O2   ','HO2       '],                                            &
       ['RTN23OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN14O2   ','HO2       '],                                            &
       ['RTN14OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN10O2   ','HO2       '],                                            &
       ['RTN10OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX28O2   ','HO2       '],                                            &
       ['RTX28OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX24O2   ','HO2       '],                                            &
       ['RTX24OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX22O2   ','HO2       '],                                            &
       ['RTX22OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRTX28O2  ','HO2       '],                                            &
       ['NRTX28OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU10AO2   ','HO2       '],                                            &
       ['RU10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['DHPR12O2  ','HO2       '],                                            &
       ['DHPR12OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['MACO3     ','HO2       '],                                            &
       ['RU10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,6,                              &
       ['DHPR12O2  ','NO        '],                                            &
       ['CARB3     ','RN8OOH    ','OH        ','NO2       ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,6,                              &
       ['DHPR12O2  ','NO3       '],                                            &
       ['CARB3     ','RN8OOH    ','OH        ','NO2       ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,5,                              &
       ['DHPR12O2  ','RO2       '],                                            &
       ['CARB3     ','RN8OOH    ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,5,                              &
       ['DHPCARB9  ','OH        '],                                            &
       ['RN8OOH    ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50045,'J',.TRUE.,0,6,                              &
       ['DHPCARB9  ','PHOTON    '],                                            &
       ['RN8OOH    ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,5,                              &
       ['OH        ','RU12OOH   '],                                            &
       ['RU10OOH   ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50046,'T',.TRUE.,0,4,                              &
       ['OH        ','NO2       '],                                            &
       ['HONO2     ','m         ','          ','          '])                  &
       ]



! Equivalent fluxes in CRI mechanism
TYPE(asad_flux_defn), PARAMETER :: cri_aer_trop_other_fluxes(56) =             &
! Extra fluxes of interest
       [ asad_flux_defn('RXN',50042,'B',.TRUE.,0,3,                            &
       ['NO3       ','C5H8      '],                                            &
       ['NRU14O2   ','          ','          ','          ']),                 &
       ! Equiv. to NO+ISO2
       asad_flux_defn('RXN',50043,'B',.TRUE.,0,3,                              &
       ['NO        ','RU14O2    '],                                            &
       ['RU14NO3   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'T',.TRUE.,0,4,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,3,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','          ','          ','          ']),                 &
! Including production of H2O2 from O3+alkene reactions (4)
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,5,                              &
       ['APINENE   ','O3        '],                                            &
       ['TNCARB26  ','H2O2      ','Sec_Org   ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,6,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB24  ','H2O2      ','Sec_Org   ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,5,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','CARB6     ','H2O2      ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,5,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','CARB6     ','H2O2      ','          ']),                 &
! ROOH production (5)
! Note not including acetic acid production - not a hydro-peroxide!
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeOO      '],                                            &
       ['MeOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','EtOO      '],                                            &
       ['EtOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeCO3     '],                                            &
       ['MeCO3H    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','i-PrOO    '],                                            &
       ['i-PrOOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                              &
       ['HO2       ','EtCO3     '],                                            &
       ['O2        ','EtCO3H    ','          ','          ']),                 &
! New peroxide production in CRI (43):
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN10O2    ','HO2       '],                                            &
       ['RN10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN13O2    ','HO2       '],                                            &
       ['RN13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN16O2    ','HO2       '],                                            &
       ['RN16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN19O2    ','HO2       '],                                            &
       ['RN19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN13AO2   ','HO2       '],                                            &
       ['RN13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN16AO2   ','HO2       '],                                            &
       ['RN16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA13O2    ','HO2       '],                                            &
       ['RA13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA16O2    ','HO2       '],                                            &
       ['RA16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA19AO2   ','HO2       '],                                            &
       ['RA19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA19CO2   ','HO2       '],                                            &
       ['RA19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HOCH2CH2O2','HO2       '],                                            &
       ['HOC2H4OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN9O2     ','HO2       '],                                            &
       ['RN9OOH    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN12O2    ','HO2       '],                                            &
       ['RN12OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN15O2    ','HO2       '],                                            &
       ['RN15OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN18O2    ','HO2       '],                                            &
       ['RN18OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN15AO2   ','HO2       '],                                            &
       ['RN15OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN18AO2   ','HO2       '],                                            &
       ['RN18OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HOCH2CO3  ','HO2       '],                                            &
       ['HOCH2CO3H ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN8O2     ','HO2       '],                                            &
       ['RN8OOH    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN11O2    ','HO2       '],                                            &
       ['RN11OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN14O2    ','HO2       '],                                            &
       ['RN14OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN17O2    ','HO2       '],                                            &
       ['RN17OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU14O2    ','HO2       '],                                            &
       ['RU14OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU12O2    ','HO2       '],                                            &
       ['RU12OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU10O2    ','HO2       '],                                            &
       ['RU10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN6O2    ','HO2       '],                                            &
       ['NRN6OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN9O2    ','HO2       '],                                            &
       ['NRN9OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN12O2   ','HO2       '],                                            &
       ['NRN12OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRU14O2   ','HO2       '],                                            &
       ['NRU14OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRU12O2   ','HO2       '],                                            &
       ['NRU12OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN28O2   ','HO2       '],                                            &
       ['RTN28OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRTN28O2  ','HO2       '],                                            &
       ['NRTN28OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN26O2   ','HO2       '],                                            &
       ['RTN26OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN25O2   ','HO2       '],                                            &
       ['RTN25OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN24O2   ','HO2       '],                                            &
       ['RTN24OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN23O2   ','HO2       '],                                            &
       ['RTN23OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN14O2   ','HO2       '],                                            &
       ['RTN14OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN10O2   ','HO2       '],                                            &
       ['RTN10OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX28O2   ','HO2       '],                                            &
       ['RTX28OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX24O2   ','HO2       '],                                            &
       ['RTX24OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX22O2   ','HO2       '],                                            &
       ['RTX22OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRTX28O2  ','HO2       '],                                            &
       ['NRTX28OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50046,'T',.TRUE.,0,4,                              &
       ['OH        ','NO2       '],                                            &
       ['HONO2     ','m         ','          ','          '])                  &
       ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_aer_trop_other_fluxes(64) =      &
! Extra fluxes of interest
       [ asad_flux_defn('RXN',50042,'B',.TRUE.,0,3,                            &
       ['NO3       ','C5H8      '],                                            &
       ['NRU14O2   ','          ','          ','          ']),                 &
       ! Equiv. to NO+ISO2
       asad_flux_defn('RXN',50043,'B',.TRUE.,0,3,                              &
       ['NO        ','RU14O2    '],                                            &
       ['RU14NO3   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'T',.TRUE.,0,4,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','O2        ','          ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,3,                              &
       ['HO2       ','HO2       '],                                            &
       ['H2O2      ','          ','          ','          ']),                 &
! Including production of H2O2 from O3+alkene reactions (4)
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,5,                              &
       ['APINENE   ','O3        '],                                            &
       ['TNCARB26  ','H2O2      ','Sec_Org   ','          ']),                 &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,6,                              &
       ['BPINENE   ','O3        '],                                            &
       ['HCHO      ','TXCARB24  ','H2O2      ','Sec_Org   ']),                 &
!CS2 O3 + C5H8 now makes H2O2
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,5,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','HCHO      ','H2O2      ','          ']),                 &
! ROOH production (5)
! Note not including acetic acid production - not a hydro-peroxide!
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeOO      '],                                            &
       ['MeOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','EtOO      '],                                            &
       ['EtOOH     ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','MeCO3     '],                                            &
       ['MeCO3H    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HO2       ','i-PrOO    '],                                            &
       ['i-PrOOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,5,                              &
       ['HO2       ','EtCO3     '],                                            &
       ['OH        ','EtCO3H    ','EtOO      ','          ']),                 &
! New peroxide production in CRI (43):
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN10O2    ','HO2       '],                                            &
       ['RN10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN13O2    ','HO2       '],                                            &
       ['RN13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN16O2    ','HO2       '],                                            &
       ['RN16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN19O2    ','HO2       '],                                            &
       ['RN19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN13AO2   ','HO2       '],                                            &
       ['RN13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN16AO2   ','HO2       '],                                            &
       ['RN16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA13O2    ','HO2       '],                                            &
       ['RA13OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA16O2    ','HO2       '],                                            &
       ['RA16OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA19AO2   ','HO2       '],                                            &
       ['RA19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RA19CO2   ','HO2       '],                                            &
       ['RA19OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HOCH2CH2O2','HO2       '],                                            &
       ['HOC2H4OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN9O2     ','HO2       '],                                            &
       ['RN9OOH    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN12O2    ','HO2       '],                                            &
       ['RN12OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN15O2    ','HO2       '],                                            &
       ['RN15OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN18O2    ','HO2       '],                                            &
       ['RN18OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN15AO2   ','HO2       '],                                            &
       ['RN15OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN18AO2   ','HO2       '],                                            &
       ['RN18OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['HOCH2CO3  ','HO2       '],                                            &
       ['HOCH2CO3H ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN8O2     ','HO2       '],                                            &
       ['RN8OOH    ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN11O2    ','HO2       '],                                            &
       ['RN11OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN14O2    ','HO2       '],                                            &
       ['RN14OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RN17O2    ','HO2       '],                                            &
       ['RN17OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU14O2    ','HO2       '],                                            &
       ['RU14OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU12O2    ','HO2       '],                                            &
       ['RU12OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU10O2    ','HO2       '],                                            &
       ['RU10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN6O2    ','HO2       '],                                            &
       ['NRN6OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN9O2    ','HO2       '],                                            &
       ['NRN9OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRN12O2   ','HO2       '],                                            &
       ['NRN12OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRU14O2   ','HO2       '],                                            &
       ['NRU14OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRU12O2   ','HO2       '],                                            &
       ['NRU12OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN28O2   ','HO2       '],                                            &
       ['RTN28OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRTN28O2  ','HO2       '],                                            &
       ['NRTN28OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,5,                              &
       ['RTN26O2   ','HO2       '],                                            &
       ['RTN26OOH  ','RTN25O2   ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN25O2   ','HO2       '],                                            &
       ['RTN25OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN24O2   ','HO2       '],                                            &
       ['RTN24OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN23O2   ','HO2       '],                                            &
       ['RTN23OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN14O2   ','HO2       '],                                            &
       ['RTN14OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTN10O2   ','HO2       '],                                            &
       ['RTN10OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX28O2   ','HO2       '],                                            &
       ['RTX28OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX24O2   ','HO2       '],                                            &
       ['RTX24OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RTX22O2   ','HO2       '],                                            &
       ['RTX22OOH  ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['NRTX28O2  ','HO2       '],                                            &
       ['NRTX28OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['RU10AO2   ','HO2       '],                                            &
       ['RU10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['DHPR12O2  ','HO2       '],                                            &
       ['DHPR12OOH ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                              &
       ['MACO3     ','HO2       '],                                            &
       ['RU10OOH   ','          ','          ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,6,                              &
       ['DHPR12O2  ','NO        '],                                            &
       ['CARB3     ','RN8OOH    ','OH        ','NO2       ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,6,                              &
       ['DHPR12O2  ','NO3       '],                                            &
       ['CARB3     ','RN8OOH    ','OH        ','NO2       ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,5,                              &
       ['DHPR12O2  ','RO2       '],                                            &
       ['CARB3     ','RN8OOH    ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,5,                              &
       ['DHPCARB9  ','OH        '],                                            &
       ['RN8OOH    ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50045,'J',.TRUE.,0,6,                              &
       ['DHPCARB9  ','PHOTON    '],                                            &
       ['RN8OOH    ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,5,                              &
       ['OH        ','RU12OOH   '],                                            &
       ['RU10OOH   ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50046,'T',.TRUE.,0,4,                              &
       ['OH        ','NO2       '],                                            &
       ['HONO2     ','m         ','          ','          '])                  &
       ]








TYPE(asad_flux_defn), PARAMETER :: asad_general_interest(8) = [                &
! For CH4 lifetime
       asad_flux_defn('RXN',50041,'B',.TRUE.,0,4,                              &
       ['OH        ','CH4       '],                                            &
       ['H2O       ','MeOO      ','          ','          ']),                 &
! Ozone STE
       asad_flux_defn('STE',50051,'X',.TRUE.,0,1,                              &
       ['O3        ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! Ozone tendency in the troposphere
       asad_flux_defn('NET',50052,'X',.TRUE.,0,1,                              &
       ['O3        ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! Ozone in the troposphere
       asad_flux_defn('OUT',50053,'X',.TRUE.,0,1,                              &
       ['O3        ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! Ozone tendency - whole atmos
       asad_flux_defn('NET',50054,'X',.FALSE.,0,1,                             &
       ['O3        ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! Tropospheric mass
       asad_flux_defn('MAS',50061,'X',.TRUE.,0,1,                              &
       ['X         ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! Tropospheric mask
       asad_flux_defn('TPM',50062,'X',.TRUE.,0,1,                              &
       ['X         ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
! Total atmospheric mass
       asad_flux_defn('MAS',50063,'X',.FALSE.,0,1,                             &
       ['X         ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: asad_trop_co_budget(21) = [                 &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,3,                              &
       ['OH        ','CO        '],                                            &
       ['HO2       ','          ','          ','          ']),                 &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['OH        ','HCHO      '],                                            &
       ['H2O       ','HO2       ','CO        ','          ']),                 &
! MGLY + OH/NO3
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','MGLY      '],                                            &
       ['MeCO3     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['NO3       ','MGLY      '],                                            &
       ['MeCO3     ','CO        ','HONO2     ','          ']),                 &
! O3 + MACR/ISOP & OTHER FLUXES
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeOO      ','HCOOH     ','CO        ','H2O2      ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,1,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,2,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NO        ','MACRO2    '],                                            &
       ['NO2       ','MeCO3     ','HACET     ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACRO2    ','MACRO2    '],                                            &
       ['HACET     ','MGLY      ','HCHO      ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NALD      '],                                            &
       ['HCHO      ','CO        ','NO2       ','          ']),                 &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['HO2       ','HO2       ','CO        ','          ']),                 &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['H2        ','CO        ','          ','          ']),                 &
! MGLY photolysis
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['MGLY      ','PHOTON    '],                                            &
       ['MeCO3     ','CO        ','HO2       ','          ']),                 &
! OTHER CO PROD PHOTOLYSIS RXNS
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['MeOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,4,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['CH4       ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['EtCHO     ','PHOTON    '],                                            &
       ['EtOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['MACR      ','PHOTON    '],                                            &
       ['MeCO3     ','HCHO      ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['MACROOH   ','PHOTON    '],                                            &
       ['HACET     ','CO        ','MGLY      ','HCHO      ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NALD      ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','NO2       ','HO2       ']),                 &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                              &
       ['CO        ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: asad_trop_co_budget_121(21) = [             &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,4,                              &
       ['OH        ','CO        '],                                            &
       ['H         ','CO2       ','          ','          ']),                 &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['OH        ','HCHO      '],                                            &
       ['H2O       ','HO2       ','CO        ','          ']),                 &
! MGLY + OH/NO3
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','MGLY      '],                                            &
       ['MeCO3     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['NO3       ','MGLY      '],                                            &
       ['MeCO3     ','CO        ','HONO2     ','          ']),                 &
! O3 + MACR/ISOP & OTHER FLUXES
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeOO      ','HCOOH     ','CO        ','H2O2      ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,1,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,2,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NO        ','MACRO2    '],                                            &
       ['NO2       ','MeCO3     ','HACET     ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACRO2    ','MACRO2    '],                                            &
       ['HACET     ','MGLY      ','HCHO      ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NALD      '],                                            &
       ['HCHO      ','CO        ','NO2       ','          ']),                 &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['HO2       ','HO2       ','CO        ','          ']),                 &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['H2        ','CO        ','          ','          ']),                 &
! MGLY photolysis
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['MGLY      ','PHOTON    '],                                            &
       ['MeCO3     ','CO        ','HO2       ','          ']),                 &
! OTHER CO PROD PHOTOLYSIS RXNS
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['MeOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,4,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['CH4       ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['EtCHO     ','PHOTON    '],                                            &
       ['EtOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['MACR      ','PHOTON    '],                                            &
       ['MeCO3     ','HCHO      ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['MACROOH   ','PHOTON    '],                                            &
       ['HACET     ','CO        ','MGLY      ','HCHO      ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NALD      ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','NO2       ','HO2       ']),                 &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                              &
       ['CO        ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: asad_ro2perm_trop_co_budget(21) = [         &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,3,                              &
       ['OH        ','CO        '],                                            &
       ['HO2       ','          ','          ','          ']),                 &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['OH        ','HCHO      '],                                            &
       ['H2O       ','HO2       ','CO        ','          ']),                 &
! MGLY + OH/NO3
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','MGLY      '],                                            &
       ['MeCO3     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['NO3       ','MGLY      '],                                            &
       ['MeCO3     ','CO        ','HONO2     ','          ']),                 &
! O3 + MACR/ISOP & OTHER FLUXES
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeOO      ','HCOOH     ','CO        ','H2O2      ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,1,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,2,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NO        ','MACRO2    '],                                            &
       ['NO2       ','MeCO3     ','HACET     ','CO        ']),                 &
  ! Alternative form for MACRO2+RO2 reaction if l_ukca_ro2_perm = TRUE
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACRO2    ','RO2       '],                                            &
       ['HACET     ','MGLY      ','HCHO      ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NALD      '],                                            &
       ['HCHO      ','CO        ','NO2       ','          ']),                 &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['HO2       ','HO2       ','CO        ','          ']),                 &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['H2        ','CO        ','          ','          ']),                 &
! MGLY photolysis
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['MGLY      ','PHOTON    '],                                            &
       ['MeCO3     ','CO        ','HO2       ','          ']),                 &
! OTHER CO PROD PHOTOLYSIS RXNS
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['MeOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,4,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['CH4       ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['EtCHO     ','PHOTON    '],                                            &
       ['EtOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['MACR      ','PHOTON    '],                                            &
       ['MeCO3     ','HCHO      ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['MACROOH   ','PHOTON    '],                                            &
       ['HACET     ','CO        ','MGLY      ','HCHO      ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NALD      ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','NO2       ','HO2       ']),                 &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                              &
       ['CO        ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: asad_ro2perm_trop_co_budget_121(21) = [     &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,4,                              &
       ['OH        ','CO        '],                                            &
       ['H         ','CO2       ','          ','          ']),                 &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['OH        ','HCHO      '],                                            &
       ['H2O       ','HO2       ','CO        ','          ']),                 &
! MGLY + OH/NO3
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','MGLY      '],                                            &
       ['MeCO3     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['NO3       ','MGLY      '],                                            &
       ['MeCO3     ','CO        ','HONO2     ','          ']),                 &
! O3 + MACR/ISOP & OTHER FLUXES
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeOO      ','HCOOH     ','CO        ','H2O2      ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,1,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,2,6,                              &
       ['O3        ','MACR      '],                                            &
       ['MGLY      ','HCOOH     ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NO        ','MACRO2    '],                                            &
       ['NO2       ','MeCO3     ','HACET     ','CO        ']),                 &
  ! Alternative form for MACRO2+RO2 reaction if l_ukca_ro2_perm = TRUE
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACRO2    ','RO2       '],                                            &
       ['HACET     ','MGLY      ','HCHO      ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NALD      '],                                            &
       ['HCHO      ','CO        ','NO2       ','          ']),                 &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['HO2       ','HO2       ','CO        ','          ']),                 &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['H2        ','CO        ','          ','          ']),                 &
! MGLY photolysis
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['MGLY      ','PHOTON    '],                                            &
       ['MeCO3     ','CO        ','HO2       ','          ']),                 &
! OTHER CO PROD PHOTOLYSIS RXNS
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['MeOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,4,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['CH4       ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['EtCHO     ','PHOTON    '],                                            &
       ['EtOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['MACR      ','PHOTON    '],                                            &
       ['MeCO3     ','HCHO      ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['MACROOH   ','PHOTON    '],                                            &
       ['HACET     ','CO        ','MGLY      ','HCHO      ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NALD      ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','NO2       ','HO2       ']),                 &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                              &
       ['CO        ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: cri_trop_co_budget(47) = [                  &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,3,                              &
       ['OH        ','CO        '],                                            &
       ['HO2       ','          ','          ','          ']),                 &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['OH        ','HCHO      '],                                            &
       ['H2O       ','HO2       ','CO        ','          ']),                 &
! MGLY + OH/NO3
! in CRI, CARB3 ~ Glyoxal, CARB6 ~ MethylGlyoxal
! No CARB6+NO3 reaction in CRI.
! OH + CARB3 -> 2CO + HO2
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','CARB6     '],                                            &
       ['MeCO3     ','CO        ','          ','          ']),                 &
! O3 + MACR/ISOP & OTHER FLUXES
!   O3 + unsaturated VOCs in CRI (8)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['BPINENE   ','O3        '],                                            &
       ['TXCARB24  ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','MeCO3     ','CO        ','OH        ']),                 &
! Miscelleneous reactions in CRI
! (These could be split into separate diags)
! First, OH + C2H2 reaction in CRI
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['C2H2      ','OH        '],                                            &
       ['HCOOH     ','CO        ','          ','          ']),                 &
! OH+PAN reactions create CO in CRI (4)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','MPAN      '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PAN       '],                                            &
       ['HCHO      ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PPAN      '],                                            &
       ['MeCHO     ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PHAN      '],                                            &
       ['HCHO      ','CO        ','NO2       ','          ']),                 &
! Secondary production from NO/NO3/OH + CARB (10)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU12O2    ','NO        '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO        '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO        '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU12O2    ','NO3       '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO3       '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO3       '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['TNCARB11  ','OH        '],                                            &
       ['RTN10O2   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['TNCARB11  ','NO3       '],                                            &
       ['RTN10O2   ','CO        ','HONO2     ','          ']),                 &
! Finally, CO from RO2+RO2 reactions (2)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['NRU12O2   ','RO2       '],                                            &
       ['NOA       ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['RTN10O2   ','RO2       '],                                            &
       ['RN8O2     ','CO        ','          ','          ']),                 &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['HO2       ','HO2       ','CO        ','          ']),                 &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['H2        ','CO        ','          ','          ']),                 &
! MGLY photolysis
!   CARB3 ~ GLY; CARB3 + hv -> 2CO + 2HO2 (3)
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB6     ','PHOTON    '],                                            &
       ['MeCO3     ','CO        ','HO2       ','          ']),                 &
! OTHER CO PROD PHOTOLYSIS RXNS
! Only one MeCHO+hv reaction in CRI (2)
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['MeOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['EtCHO     ','PHOTON    '],                                            &
       ['EtOO      ','HO2       ','CO        ','          ']),                 &
! CRI species photolysis CO production (7)
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HOCH2CHO  ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['MeCO3     ','HOCH2CHO  ','CO        ','HO2       ']),                 &
! NUCARB12 + photon -> NOA + 2*CO + 2*HO2
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['NOA       ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['NOA       ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['TNCARB10  ','PHOTON    '],                                            &
       ['MeCO3     ','MeCO3     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NRU12OOH  ','PHOTON    '],                                            &
       ['NOA       ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['RTN10OOH  ','PHOTON    '],                                            &
       ['RN8O2     ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['TNCARB11  ','PHOTON    '],                                            &
       ['RTN10O2   ','CO        ','HO2       ','          ']),                 &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                              &
       ['CO        ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: cri_strat2_trop_co_budget(81) = [           &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,3,                              &
       ['OH        ','CO        '],                                            &
       ['HO2       ','          ','          ','          ']),                 &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['OH        ','HCHO      '],                                            &
       ['H2O       ','HO2       ','CO        ','          ']),                 &
! MGLY + OH/NO3
! in CRI, CARB3 ~ Glyoxal, CARB6 ~ MethylGlyoxal
! No CARB6/CARB3 +NO3 reaction in CRI.
! OH + CARB3 -> 2CO + HO2
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','OH        ','          ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','CARB6     '],                                            &
       ['MeCO3     ','CO        ','          ','          ']),                 &
! O3 + MACR/ISOP & OTHER FLUXES
!   O3 + unsaturated VOCs in CRI
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeCO3     ','HCHO      ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['HCHO      ','MeOO      ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['BPINENE   ','O3        '],                                            &
       ['TXCARB24  ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB6     ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB7     ','CARB3     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['DHPR12OOH ','OH        '],                                            &
       ['DHPCARB9  ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['HPUCARB12 ','OH        '],                                            &
       ['HUCARB9   ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['HPUCARB12 ','NO3       '],                                            &
       ['HUCARB9   ','CO        ','OH        ','HONO2     ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['HUCARB9   ','OH        '],                                            &
       ['CARB6     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU10AO2   ','NO        '],                                            &
       ['CARB7     ','HO2       ','CO        ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RU10NO3   ','OH        '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU10AO2   ','NO3       '],                                            &
       ['CARB7     ','HO2       ','CO        ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO3       '],                                            &
       ['MeOO      ','HCHO      ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO3       '],                                            &
       ['MeOO      ','HCHO      ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','HO2       '],                                            &
       ['MeOO      ','CO        ','HCHO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','RO2       '],                                            &
       ['MeOO      ','CO        ','HCHO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO        '],                                            &
       ['MeOO      ','CO        ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO        '],                                            &
       ['MeOO      ','CO        ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU12O2    ','m         '],                                            &
       ['DHCARB9   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU12O2    ','m         '],                                            &
       ['DHCARB9   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['DHPR12O2  ','m         '],                                            &
       ['DHPCARB9  ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['DHPR12O2  ','m         '],                                            &
       ['DHPCARB9  ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU10AO2   ','m         '],                                            &
       ['CARB7     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU10AO2   ','m         '],                                            &
       ['CARB7     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU10AO2   ','RO2       '],                                            &
       ['CARB7     ','CO        ','HO2       ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['DHPCARB9  ','OH        '],                                            &
       ['RN8OOH    ','CO        ','OH        ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU12OOH   ','OH        '],                                            &
       ['RU10OOH   ','CO        ','HO2       ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU10OOH   ','OH        '],                                            &
       ['CARB7     ','CO        ','OH        ','          '] ),                &


! Miscelleneous reactions in CRI
! (These could be split into separate diags)
! First, OH + C2H2 reaction in CRI
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['C2H2      ','OH        '],                                            &
       ['HCOOH     ','CO        ','          ','          ']),                 &
! OH+PAN reactions create CO in CRI (4)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','MPAN      '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PAN       '],                                            &
       ['HCHO      ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PPAN      '],                                            &
       ['MeCHO     ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PHAN      '],                                            &
       ['HCHO      ','CO        ','NO2       ','          ']),                 &
!! Secondary production from NO/NO3/OH + CARB (10)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO        '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO        '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO3       '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO3       '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
! Finally, CO from RO2+RO2 reactions (2)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['NRU12O2   ','RO2       '],                                            &
       ['NOA       ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['RTN10O2   ','RO2       '],                                            &
       ['RN8O2     ','CO        ','          ','          ']),                 &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['HO2       ','HO2       ','CO        ','          ']),                 &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['H2        ','CO        ','          ','          ']),                 &
! MGLY photolysis
!   CARB3 ~ GLY; CARB3 + hv -> 2CO + 2HO2 (3)
      asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                               &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB6     ','PHOTON    '],                                            &
       ['MeCO3     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','H2        ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','H2        ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,4,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','          ','          ']),                 &
! OTHER CO PROD PHOTOLYSIS RXNS
! Only one MeCHO+hv reaction in CRI (2)
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['MeOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['EtCHO     ','PHOTON    '],                                            &
       ['EtOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HOCH2CHO  ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['MeCO3     ','HOCH2CHO  ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['HUCARB9   ','CO        ','NO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['TNCARB10  ','PHOTON    '],                                            &
       ['MeCO3     ','MeCO3     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['RTN10OOH  ','PHOTON    '],                                            &
       ['RN8O2     ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHPCARB9  ','PHOTON    '],                                            &
       ['RN8OOH    ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['HUCARB9   ','CO        ','OH        ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HUCARB9   ','PHOTON    '],                                            &
       ['CARB6     ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHPR12OOH ','PHOTON    '],                                            &
       ['DHPCARB9  ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHCARB9   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                              &
       ['CO        ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: cri_strat2_trop_co_budget_121(81) = [       &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,4,                              &
       ['OH        ','CO        '],                                            &
       ['H         ','CO2       ','          ','          ']),                 &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['OH        ','HCHO      '],                                            &
       ['H2O       ','HO2       ','CO        ','          ']),                 &
! MGLY + OH/NO3
! in CRI, CARB3 ~ Glyoxal, CARB6 ~ MethylGlyoxal
! No CARB6/CARB3 +NO3 reaction in CRI.
! OH + CARB3 -> 2CO + HO2
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','OH        ','          ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','CARB6     '],                                            &
       ['MeCO3     ','CO        ','          ','          ']),                 &
! O3 + MACR/ISOP & OTHER FLUXES
!   O3 + unsaturated VOCs in CRI
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeCO3     ','HCHO      ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['HCHO      ','MeOO      ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['BPINENE   ','O3        '],                                            &
       ['TXCARB24  ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB6     ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB7     ','CARB3     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['DHPR12OOH ','OH        '],                                            &
       ['DHPCARB9  ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['HPUCARB12 ','OH        '],                                            &
       ['HUCARB9   ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['HPUCARB12 ','NO3       '],                                            &
       ['HUCARB9   ','CO        ','OH        ','HONO2     ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['HUCARB9   ','OH        '],                                            &
       ['CARB6     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU10AO2   ','NO        '],                                            &
       ['CARB7     ','HO2       ','CO        ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RU10NO3   ','OH        '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU10AO2   ','NO3       '],                                            &
       ['CARB7     ','HO2       ','CO        ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO3       '],                                            &
       ['MeOO      ','HCHO      ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO3       '],                                            &
       ['MeOO      ','HCHO      ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','HO2       '],                                            &
       ['MeOO      ','CO        ','HCHO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','RO2       '],                                            &
       ['MeOO      ','CO        ','HCHO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO        '],                                            &
       ['MeOO      ','CO        ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO        '],                                            &
       ['MeOO      ','CO        ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU12O2    ','m         '],                                            &
       ['DHCARB9   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU12O2    ','m         '],                                            &
       ['DHCARB9   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['DHPR12O2  ','m         '],                                            &
       ['DHPCARB9  ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['DHPR12O2  ','m         '],                                            &
       ['DHPCARB9  ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU10AO2   ','m         '],                                            &
       ['CARB7     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU10AO2   ','m         '],                                            &
       ['CARB7     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU10AO2   ','RO2       '],                                            &
       ['CARB7     ','CO        ','HO2       ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['DHPCARB9  ','OH        '],                                            &
       ['RN8OOH    ','CO        ','OH        ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU12OOH   ','OH        '],                                            &
       ['RU10OOH   ','CO        ','HO2       ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU10OOH   ','OH        '],                                            &
       ['CARB7     ','CO        ','OH        ','          '] ),                &


! Miscelleneous reactions in CRI
! (These could be split into separate diags)
! First, OH + C2H2 reaction in CRI
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['C2H2      ','OH        '],                                            &
       ['HCOOH     ','CO        ','          ','          ']),                 &
! OH+PAN reactions create CO in CRI (4)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','MPAN      '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PAN       '],                                            &
       ['HCHO      ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PPAN      '],                                            &
       ['MeCHO     ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PHAN      '],                                            &
       ['HCHO      ','CO        ','NO2       ','          ']),                 &
!! Secondary production from NO/NO3/OH + CARB (10)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO        '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO        '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO3       '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO3       '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
! Finally, CO from RO2+RO2 reactions (2)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['NRU12O2   ','RO2       '],                                            &
       ['NOA       ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['RTN10O2   ','RO2       '],                                            &
       ['RN8O2     ','CO        ','          ','          ']),                 &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['HO2       ','HO2       ','CO        ','          ']),                 &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['H2        ','CO        ','          ','          ']),                 &
! MGLY photolysis
!   CARB3 ~ GLY; CARB3 + hv -> 2CO + 2HO2 (3)
      asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                               &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB6     ','PHOTON    '],                                            &
       ['MeCO3     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','H2        ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','H2        ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,4,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','          ','          ']),                 &
! OTHER CO PROD PHOTOLYSIS RXNS
! Only one MeCHO+hv reaction in CRI (2)
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['MeOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['EtCHO     ','PHOTON    '],                                            &
       ['EtOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HOCH2CHO  ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['MeCO3     ','HOCH2CHO  ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['HUCARB9   ','CO        ','NO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['TNCARB10  ','PHOTON    '],                                            &
       ['MeCO3     ','MeCO3     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['RTN10OOH  ','PHOTON    '],                                            &
       ['RN8O2     ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHPCARB9  ','PHOTON    '],                                            &
       ['RN8OOH    ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['HUCARB9   ','CO        ','OH        ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HUCARB9   ','PHOTON    '],                                            &
       ['CARB6     ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHPR12OOH ','PHOTON    '],                                            &
       ['DHPCARB9  ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHCARB9   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                              &
       ['CO        ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: cri_aer_trop_co_budget(47) = [              &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,3,                              &
       ['OH        ','CO        '],                                            &
       ['HO2       ','          ','          ','          ']),                 &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['OH        ','HCHO      '],                                            &
       ['H2O       ','HO2       ','CO        ','          ']),                 &
! MGLY + OH/NO3
! in CRI, CARB3 ~ Glyoxal, CARB6 ~ MethylGlyoxal
! No CARB6+NO3 reaction in CRI.
! OH + CARB3 -> 2CO + HO2
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','CARB6     '],                                            &
       ['MeCO3     ','CO        ','          ','          ']),                 &
! O3 + MACR/ISOP & OTHER FLUXES
!   O3 + unsaturated VOCs in CRI (8)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['BPINENE   ','O3        '],                                            &
       ['TXCARB24  ','CO        ','Sec_Org   ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['HOCH2CHO  ','MeCO3     ','CO        ','OH        ']),                 &
! Miscelleneous reactions in CRI
! (These could be split into separate diags)
! First, OH + C2H2 reaction in CRI
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['C2H2      ','OH        '],                                            &
       ['HCOOH     ','CO        ','          ','          ']),                 &
! OH+PAN reactions create CO in CRI (4)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','MPAN      '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PAN       '],                                            &
       ['HCHO      ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PPAN      '],                                            &
       ['MeCHO     ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PHAN      '],                                            &
       ['HCHO      ','CO        ','NO2       ','          ']),                 &
! Secondary production from NO/NO3/OH + CARB (10)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU12O2    ','NO        '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO        '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO        '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU12O2    ','NO3       '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO3       '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO3       '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['TNCARB11  ','OH        '],                                            &
       ['RTN10O2   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['TNCARB11  ','NO3       '],                                            &
       ['RTN10O2   ','CO        ','HONO2     ','          ']),                 &
! Finally, CO from RO2+RO2 reactions (2)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['NRU12O2   ','RO2       '],                                            &
       ['NOA       ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['RTN10O2   ','RO2       '],                                            &
       ['RN8O2     ','CO        ','          ','          ']),                 &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['HO2       ','HO2       ','CO        ','          ']),                 &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['H2        ','CO        ','          ','          ']),                 &
! MGLY photolysis
!   CARB3 ~ GLY; CARB3 + hv -> 2CO + 2HO2 (3)
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB6     ','PHOTON    '],                                            &
       ['MeCO3     ','CO        ','HO2       ','          ']),                 &
! OTHER CO PROD PHOTOLYSIS RXNS
! Only one MeCHO+hv reaction in CRI (2)
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['MeOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['EtCHO     ','PHOTON    '],                                            &
       ['EtOO      ','HO2       ','CO        ','          ']),                 &
! CRI species photolysis CO production (7)
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HOCH2CHO  ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['MeCO3     ','HOCH2CHO  ','CO        ','HO2       ']),                 &
! NUCARB12 + photon -> NOA + 2*CO + 2*HO2
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['NOA       ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['NOA       ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['TNCARB10  ','PHOTON    '],                                            &
       ['MeCO3     ','MeCO3     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NRU12OOH  ','PHOTON    '],                                            &
       ['NOA       ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['RTN10OOH  ','PHOTON    '],                                            &
       ['RN8O2     ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['TNCARB11  ','PHOTON    '],                                            &
       ['RTN10O2   ','CO        ','HO2       ','          ']),                 &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                              &
       ['CO        ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: cri_strat2_aer_trop_co_budget(81) = [       &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,3,                              &
       ['OH        ','CO        '],                                            &
       ['HO2       ','          ','          ','          ']),                 &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['OH        ','HCHO      '],                                            &
       ['H2O       ','HO2       ','CO        ','          ']),                 &
! MGLY + OH/NO3
! in CRI, CARB3 ~ Glyoxal, CARB6 ~ MethylGlyoxal
! No CARB6+NO3 reaction in CRI.
! OH + CARB3 -> 2CO + HO2
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','OH        ','          ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','CARB6     '],                                            &
       ['MeCO3     ','CO        ','          ','          ']),                 &
! O3 + MACR/ISOP & OTHER FLUXES
!   O3 + unsaturated VOCs in CRI
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','HO2       ','OH        ']),                 &
       ! additional O3 + C5H8 --> CO reactions
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeCO3     ','HCHO      ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['HCHO      ','MeOO      ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['BPINENE   ','O3        '],                                            &
       ['TXCARB24  ','CO        ','Sec_Org   ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB6     ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB7     ','CARB3     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['DHPR12OOH ','OH        '],                                            &
       ['DHPCARB9  ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['HPUCARB12 ','OH        '],                                            &
       ['HUCARB9   ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['HPUCARB12 ','NO3       '],                                            &
       ['HUCARB9   ','CO        ','OH        ','HONO2     ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['HUCARB9   ','OH        '],                                            &
       ['CARB6     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU10AO2   ','NO        '],                                            &
       ['CARB7     ','HO2       ','CO        ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RU10NO3   ','OH        '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU10AO2   ','NO3       '],                                            &
       ['CARB7     ','HO2       ','CO        ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO3       '],                                            &
       ['MeOO      ','HCHO      ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO3       '],                                            &
       ['MeOO      ','HCHO      ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','HO2       '],                                            &
       ['MeOO      ','CO        ','HCHO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','RO2       '],                                            &
       ['MeOO      ','CO        ','HCHO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO        '],                                            &
       ['MeOO      ','CO        ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO        '],                                            &
       ['MeOO      ','CO        ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU12O2    ','m         '],                                            &
       ['DHCARB9   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU12O2    ','m         '],                                            &
       ['DHCARB9   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['DHPR12O2  ','m         '],                                            &
       ['DHPCARB9  ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['DHPR12O2  ','m         '],                                            &
       ['DHPCARB9  ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU10AO2   ','m         '],                                            &
       ['CARB7     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU10AO2   ','m         '],                                            &
       ['CARB7     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU10AO2   ','RO2       '],                                            &
       ['CARB7     ','CO        ','HO2       ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['DHPCARB9  ','OH        '],                                            &
       ['RN8OOH    ','CO        ','OH        ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU12OOH   ','OH        '],                                            &
       ['RU10OOH   ','CO        ','HO2       ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU10OOH   ','OH        '],                                            &
       ['CARB7     ','CO        ','OH        ','          '] ),                &
! Miscelleneous reactions in CRI
! (These could be split into separate diags)
! First, OH + C2H2 reaction in CRI
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['C2H2      ','OH        '],                                            &
       ['HCOOH     ','CO        ','          ','          ']),                 &
! OH+PAN reactions create CO in CRI (4)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','MPAN      '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PAN       '],                                            &
       ['HCHO      ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PPAN      '],                                            &
       ['MeCHO     ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PHAN      '],                                            &
       ['HCHO      ','CO        ','NO2       ','          ']),                 &
! Secondary production from NO/NO3/OH + CARB (10)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO        '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO        '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO3       '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO3       '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
! Finally, CO from RO2+RO2 reactions (2)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['NRU12O2   ','RO2       '],                                            &
       ['NOA       ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['RTN10O2   ','RO2       '],                                            &
       ['RN8O2     ','CO        ','          ','          ']),                 &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['HO2       ','HO2       ','CO        ','          ']),                 &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['H2        ','CO        ','          ','          ']),                 &
! MGLY photolysis
!   CARB3 ~ GLY; CARB3 + hv -> 2CO + 2HO2 (3)
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB6     ','PHOTON    '],                                            &
       ['MeCO3     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','H2        ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','H2        ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,4,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','          ','          ']),                 &
! OTHER CO PROD PHOTOLYSIS RXNS
! Only one MeCHO+hv reaction in CRI (2)
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['MeOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['EtCHO     ','PHOTON    '],                                            &
       ['EtOO      ','HO2       ','CO        ','          ']),                 &
! CRI species photolysis CO production (7)
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HOCH2CHO  ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['MeCO3     ','HOCH2CHO  ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['HUCARB9   ','CO        ','NO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['TNCARB10  ','PHOTON    '],                                            &
       ['MeCO3     ','MeCO3     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['RTN10OOH  ','PHOTON    '],                                            &
       ['RN8O2     ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHPCARB9  ','PHOTON    '],                                            &
       ['RN8OOH    ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['HUCARB9   ','CO        ','OH        ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HUCARB9   ','PHOTON    '],                                            &
       ['CARB6     ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHPR12OOH ','PHOTON    '],                                            &
       ['DHPCARB9  ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHCARB9   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                              &
       ['CO        ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: cri_strat2_aer_trop_co_budget_121(81) = [   &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,4,                              &
       ['OH        ','CO        '],                                            &
       ['H         ','CO2       ','          ','          ']),                 &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['NO3       ','HCHO      '],                                            &
       ['HONO2     ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                              &
       ['OH        ','HCHO      '],                                            &
       ['H2O       ','HO2       ','CO        ','          ']),                 &
! MGLY + OH/NO3
! in CRI, CARB3 ~ Glyoxal, CARB6 ~ MethylGlyoxal
! No CARB6+NO3 reaction in CRI.
! OH + CARB3 -> 2CO + HO2
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','CARB3     '],                                            &
       ['CO        ','OH        ','          ','          ']),                 &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                              &
       ['OH        ','CARB6     '],                                            &
       ['MeCO3     ','CO        ','          ','          ']),                 &
! O3 + MACR/ISOP & OTHER FLUXES
!   O3 + unsaturated VOCs in CRI
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','HO2       ','OH        ']),                 &
       ! additional O3 + C5H8 --> CO reactions
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['MeCO3     ','HCHO      ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C5H8      '],                                            &
       ['HCHO      ','MeOO      ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['O3        ','C5H8      '],                                            &
       ['UCARB10   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C2H4      '],                                            &
       ['HCHO      ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','C3H6      '],                                            &
       ['HCHO      ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','TBUT2ENE  '],                                            &
       ['MeCHO     ','CO        ','MeOO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['BPINENE   ','O3        '],                                            &
       ['TXCARB24  ','CO        ','Sec_Org   ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB10   '],                                            &
       ['HCHO      ','MeCO3     ','CO        ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB6     ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['O3        ','UCARB12   '],                                            &
       ['CARB7     ','CARB3     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['DHPR12OOH ','OH        '],                                            &
       ['DHPCARB9  ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['HPUCARB12 ','OH        '],                                            &
       ['HUCARB9   ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['HPUCARB12 ','NO3       '],                                            &
       ['HUCARB9   ','CO        ','OH        ','HONO2     ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['HUCARB9   ','OH        '],                                            &
       ['CARB6     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU10AO2   ','NO        '],                                            &
       ['CARB7     ','HO2       ','CO        ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RU10NO3   ','OH        '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['RU10AO2   ','NO3       '],                                            &
       ['CARB7     ','HO2       ','CO        ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO3       '],                                            &
       ['MeOO      ','HCHO      ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO3       '],                                            &
       ['MeOO      ','HCHO      ','HO2       ','CO        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','HO2       '],                                            &
       ['MeOO      ','CO        ','HCHO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','RO2       '],                                            &
       ['MeOO      ','CO        ','HCHO      ','OH        ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO        '],                                            &
       ['MeOO      ','CO        ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['MACO3     ','NO        '],                                            &
       ['MeOO      ','CO        ','HCHO      ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU12O2    ','m         '],                                            &
       ['DHCARB9   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU12O2    ','m         '],                                            &
       ['DHCARB9   ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['DHPR12O2  ','m         '],                                            &
       ['DHPCARB9  ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['DHPR12O2  ','m         '],                                            &
       ['DHPCARB9  ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU10AO2   ','m         '],                                            &
       ['CARB7     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['RU10AO2   ','m         '],                                            &
       ['CARB7     ','CO        ','          ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU10AO2   ','RO2       '],                                            &
       ['CARB7     ','CO        ','HO2       ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['DHPCARB9  ','OH        '],                                            &
       ['RN8OOH    ','CO        ','OH        ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU12OOH   ','OH        '],                                            &
       ['RU10OOH   ','CO        ','HO2       ','          '] ),                &
       asad_flux_defn('RXN',50074,'B',.FALSE.,0,5,                             &
       ['RU10OOH   ','OH        '],                                            &
       ['CARB7     ','CO        ','OH        ','          '] ),                &
! Miscelleneous reactions in CRI
! (These could be split into separate diags)
! First, OH + C2H2 reaction in CRI
       asad_flux_defn('RXN',50074,'T',.TRUE.,0,4,                              &
       ['C2H2      ','OH        '],                                            &
       ['HCOOH     ','CO        ','          ','          ']),                 &
! OH+PAN reactions create CO in CRI (4)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','MPAN      '],                                            &
       ['CARB7     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PAN       '],                                            &
       ['HCHO      ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PPAN      '],                                            &
       ['MeCHO     ','CO        ','NO2       ','H2O       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['OH        ','PHAN      '],                                            &
       ['HCHO      ','CO        ','NO2       ','          ']),                 &
! Secondary production from NO/NO3/OH + CARB (10)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO        '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO        '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                              &
       ['NRU12O2   ','NO3       '],                                            &
       ['NOA       ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['RTN10O2   ','NO3       '],                                            &
       ['RN8O2     ','CO        ','NO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['OH        ','NRU12OOH  '],                                            &
       ['NOA       ','CO        ','OH        ','          ']),                 &
! Finally, CO from RO2+RO2 reactions (2)
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                              &
       ['NRU12O2   ','RO2       '],                                            &
       ['NOA       ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,4,                              &
       ['RTN10O2   ','RO2       '],                                            &
       ['RN8O2     ','CO        ','          ','          ']),                 &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['HO2       ','HO2       ','CO        ','          ']),                 &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                              &
       ['HCHO      ','PHOTON    '],                                            &
       ['H2        ','CO        ','          ','          ']),                 &
! MGLY photolysis
!   CARB3 ~ GLY; CARB3 + hv -> 2CO + 2HO2 (3)
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,6,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB6     ','PHOTON    '],                                            &
       ['MeCO3     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','H2        ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['CO        ','CO        ','H2        ','          ']),                 &
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,4,                              &
       ['CARB3     ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','          ','          ']),                 &
! OTHER CO PROD PHOTOLYSIS RXNS
! Only one MeCHO+hv reaction in CRI (2)
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['MeCHO     ','PHOTON    '],                                            &
       ['MeOO      ','HO2       ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['EtCHO     ','PHOTON    '],                                            &
       ['EtOO      ','HO2       ','CO        ','          ']),                 &
! CRI species photolysis CO production (7)
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HOCH2CHO  ','PHOTON    '],                                            &
       ['HCHO      ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['MeCO3     ','HOCH2CHO  ','CO        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['UCARB12   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['HUCARB9   ','CO        ','NO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['TNCARB10  ','PHOTON    '],                                            &
       ['MeCO3     ','MeCO3     ','CO        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                              &
       ['RTN10OOH  ','PHOTON    '],                                            &
       ['RN8O2     ','CO        ','OH        ','          ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHPCARB9  ','PHOTON    '],                                            &
       ['RN8OOH    ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['HUCARB9   ','CO        ','OH        ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HPUCARB12 ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','OH        ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['HUCARB9   ','PHOTON    '],                                            &
       ['CARB6     ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHPR12OOH ','PHOTON    '],                                            &
       ['DHPCARB9  ','CO        ','OH        ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['DHCARB9   ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','HO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                              &
       ['NUCARB12  ','PHOTON    '],                                            &
       ['CARB7     ','CO        ','HO2       ','NO2       ']),                 &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                              &
       ['CO        ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER :: asad_lightning_diags(5) = [                 &
! WARNING: Lightning emissions are calculated as NO2, but emitted as NO
       asad_flux_defn('EMS',50081,'L',.FALSE.,0,1,                             &
       ['NO        ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('LGT',50082,'T',.FALSE.,0,1,                             &
       ['X         ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('LGT',50083,'G',.FALSE.,0,1,                             &
       ['X         ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('LGT',50084,'C',.FALSE.,0,1,                             &
       ['X         ','          '],                                            &
       ['          ','          ','          ','          ']),                 &
       asad_flux_defn('LGT',50085,'N',.FALSE.,0,1,                             &
       ['X         ','          '],                                            &
       ['          ','          ','          ','          '])                  &
       ]


TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                                     asad_strat_oh_prod(31) = [                &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,5,                                    &
['Cl        ','HOCl      '],                                                   &
['Cl        ','Cl        ','OH        ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(3P)     ','HBr       '],                                                  &
['OH        ','Br        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(3P)     ','HCl       '],                                                  &
['OH        ','Cl        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(3P)     ','HOCl      '],                                                  &
['OH        ','ClO       ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','HBr       '],                                                  &
['OH        ','Br        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','HCl       '],                                                  &
['OH        ','Cl        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','MeBr      '],                                                  &
['OH        ','Br        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['Cl        ','HO2       '],                                                  &
['ClO       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['Cl        ','HO2       '],                                                  &
['ClO       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['NO3       ','HO2       '],                                                  &
['NO2       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O3        ','HO2       '],                                                  &
['O2        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(3P)     ','HO2       '],                                                  &
['OH        ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','CH4       '],                                                  &
['OH        ','MeOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','H2        '],                                                  &
['OH        ','H         ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','H2O       '],                                                  &
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','H2O       '],                                                  &
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['H         ','HO2       '],                                                  &
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['H         ','HO2       '],                                                  &
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['H         ','O3        '],                                                  &
['OH        ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['H         ','NO2       '],                                                  &
['OH        ','NO        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(3P)     ','H2        '],                                                  &
['OH        ','H         ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(3P)     ','H2O2      '],                                                  &
['OH        ','HO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,5,                                    &
 ['O(3P)     ','HCHO      '],                                                  &
['OH        ','CO        ','HO2       ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['NO        ','HO2       '],                                                  &
['NO2       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['HOBr      ','PHOTON    '],                                                 &
['OH        ','Br        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['HOCl      ','PHOTON    '],                                                 &
['OH        ','Cl        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['HONO2     ','PHOTON    '],                                                 &
['NO2       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['HO2NO2    ','PHOTON    '],                                                 &
['NO3       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['H2O       ','PHOTON    '],                                                 &
['OH        ','H         ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['H2O2      ','PHOTON    '],                                                 &
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,5,                                    &
  ['MeOOH     ','PHOTON    '],                                                 &
['HCHO      ','HO2       ','OH        ','          '])                         &
]


TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                                     asad_strat_oh_prod_121(29) = [            &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,5,                                    &
['Cl        ','HOCl      '],                                                   &
['Cl        ','Cl        ','OH        ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(3P)     ','HBr       '],                                                  &
['OH        ','Br        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(3P)     ','HCl       '],                                                  &
['OH        ','Cl        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(3P)     ','HOCl      '],                                                  &
['OH        ','ClO       ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','HBr       '],                                                  &
['OH        ','Br        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','HCl       '],                                                  &
['OH        ','Cl        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','MeBr      '],                                                  &
['OH        ','Br        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['Cl        ','HO2       '],                                                  &
['ClO       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['Cl        ','HO2       '],                                                  &
['ClO       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['NO3       ','HO2       '],                                                  &
['NO2       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O3        ','HO2       '],                                                  &
['O2        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(3P)     ','HO2       '],                                                  &
['OH        ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','CH4       '],                                                  &
['OH        ','MeOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','H2        '],                                                  &
['OH        ','H         ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','H2O       '],                                                  &
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(1D)     ','H2O       '],                                                  &
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['H         ','HO2       '],                                                  &
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['H         ','HO2       '],                                                  &
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['H         ','O3        '],                                                  &
['OH        ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['H         ','NO2       '],                                                  &
['OH        ','NO        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['O(3P)     ','H2O2      '],                                                  &
['OH        ','HO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                                    &
 ['NO        ','HO2       '],                                                  &
['NO2       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['HOBr      ','PHOTON    '],                                                 &
['OH        ','Br        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['HOCl      ','PHOTON    '],                                                 &
['OH        ','Cl        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['HONO2     ','PHOTON    '],                                                 &
['NO2       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['HO2NO2    ','PHOTON    '],                                                 &
['NO3       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['H2O       ','PHOTON    '],                                                 &
['OH        ','H         ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                                    &
  ['H2O2      ','PHOTON    '],                                                 &
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,5,                                    &
  ['MeOOH     ','PHOTON    '],                                                 &
['HCHO      ','HO2       ','OH        ','          '])                         &
]


! OH loss reactions
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                                       asad_strat_oh_loss(26) = [              &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['BrO       ','OH        '],                                                 &
  ['Br        ','HO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['ClO       ','OH        '],                                                 &
  ['Cl        ','HO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','ClO       '],                                                 &
  ['HCl       ','O2        ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','HBr       '],                                                 &
  ['H2O       ','Br        ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','HCl       '],                                                 &
  ['H2O       ','Cl        ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','HOCl      '],                                                 &
  ['ClO       ','H2O       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','OClO      '],                                                 &
  ['HOCl      ','O2        ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','ClONO2    '],                                                 &
  ['HOCl      ','NO3       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','MeBr      '],                                                 &
  ['Br        ','H2O       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
   ['OH        ','O3        '],                                                &
  ['HO2       ','O2        ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['O(3P)     ','OH        '],                                                 &
  ['O2        ','H         ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','CH4       '],                                                 &
  ['H2O       ','MeOO      ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,3,                                  &
  ['OH        ','CO        '],                                                 &
  ['HO2       ','          ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
   ['OH        ','NO3       '],                                                &
  ['HO2       ','NO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                                  &
  ['OH        ','NO2       '],                                                 &
  ['HONO2     ','m         ','          ','          ']),                      &
  Asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','HONO2     '],                                                 &
  ['NO3       ','H2O       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,5,                                  &
  ['OH        ','HO2NO2    '],                                                 &
  ['H2O       ','NO2       ','O2        ','          ']),                      &
  asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                                  &
   ['OH        ','OH        '],                                                &
  ['H2O2      ','m         ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                                  &
  ['OH        ','OH        '],                                                 &
  ['H2O2      ','m         ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','OH        '],                                                 &
  ['H2O       ','O(3P)     ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','OH        '],                                                 &
  ['H2O       ','O(3P)     ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,3,                                  &
  ['OH        ','HO2       '],                                                 &
  ['H2O       ','          ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
   ['OH        ','H2O2      '],                                                &
  ['H2O       ','HO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,5,                                  &
   ['OH        ','HCHO      '],                                                &
  ['H2O       ','CO        ','HO2       ','          ']),                      &
  Asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','H2        '],                                                 &
  ['H2O       ','H02       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','MeOOH     '],                                                 &
  ['MeOO      ','H2O       ','          ','          '])                       &
  ]


! OH loss reactions
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                                       asad_strat_oh_loss_121(26) = [          &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['BrO       ','OH        '],                                                 &
  ['Br        ','HO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['ClO       ','OH        '],                                                 &
  ['Cl        ','HO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','ClO       '],                                                 &
  ['HCl       ','O2        ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','HBr       '],                                                 &
  ['H2O       ','Br        ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','HCl       '],                                                 &
  ['H2O       ','Cl        ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','HOCl      '],                                                 &
  ['ClO       ','H2O       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','OClO      '],                                                 &
  ['HOCl      ','O2        ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','ClONO2    '],                                                 &
  ['HOCl      ','NO3       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','MeBr      '],                                                 &
  ['Br        ','H2O       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
   ['OH        ','O3        '],                                                &
  ['HO2       ','O2        ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['O(3P)     ','OH        '],                                                 &
  ['O2        ','H         ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','CH4       '],                                                 &
  ['H2O       ','MeOO      ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','CO        '],                                                 &
  ['H         ','CO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
   ['OH        ','NO3       '],                                                &
  ['HO2       ','NO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                                  &
  ['OH        ','NO2       '],                                                 &
  ['HONO2     ','m         ','          ','          ']),                      &
  Asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','HONO2     '],                                                 &
  ['NO3       ','H2O       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,5,                                  &
  ['OH        ','HO2NO2    '],                                                 &
  ['H2O       ','NO2       ','O2        ','          ']),                      &
  asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                                  &
   ['OH        ','OH        '],                                                &
  ['H2O2      ','m         ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                                  &
  ['OH        ','OH        '],                                                 &
  ['H2O2      ','m         ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','OH        '],                                                 &
  ['H2O       ','O(3P)     ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','OH        '],                                                 &
  ['H2O       ','O(3P)     ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,3,                                  &
  ['OH        ','HO2       '],                                                 &
  ['H2O       ','          ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
   ['OH        ','H2O2      '],                                                &
  ['H2O       ','HO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,5,                                  &
   ['OH        ','HCHO      '],                                                &
  ['H2O       ','CO        ','HO2       ','          ']),                      &
  Asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','H2        '],                                                 &
  ['H2O       ','H         ','          ','          ']),                      &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                                  &
  ['OH        ','MeOOH     '],                                                 &
  ['MeOO      ','H2O       ','          ','          '])                       &
  ]


! Simple strat ozone budget
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                                   asad_strat_o3_budget(20) = [                &
! production
asad_flux_defn('RXN',50101,'J',.FALSE.,0,4,                                    &
['O2        ','PHOTON    '],                                                   &
['O(3P)     ','O(3P)     ','          ','          ']),                        &
asad_flux_defn('RXN',50101,'J',.FALSE.,0,4,                                    &
['O2        ','PHOTON    '],                                                   &
['O(3P)     ','O(1D)     ','          ','          ']),                        &
asad_flux_defn('RXN',50102,'B',.FALSE.,0,4,                                    &
['HO2       ','NO        '],                                                   &
['OH        ','NO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50103,'B',.FALSE.,0,5,                                    &
['MeOO      ','NO        '],                                                   &
['HO2       ','HCHO      ','NO2       ','          ']),                        &
asad_flux_defn('RXN',50104,'B',.FALSE.,0,4,                                    &
['OH        ','HONO2     '],                                                   &
['H2O       ','NO3       ','          ','          ']),                        &
! loss
asad_flux_defn('RXN',50111,'J',.FALSE.,0,5,                                    &
['Cl2O2     ','PHOTON    '],                                                   &
['Cl        ','Cl        ','O2        ','          ']),                        &
asad_flux_defn('RXN',50112,'B',.FALSE.,0,4,                                    &
['BrO       ','ClO       '],                                                   &
['Br        ','OClO      ','          ','          ']),                        &
asad_flux_defn('RXN',50113,'B',.FALSE.,0,4,                                    &
['HO2       ','O3        '],                                                   &
['OH        ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50114,'B',.FALSE.,0,4,                                    &
['ClO       ','HO2       '],                                                   &
['HOCl      ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50115,'B',.FALSE.,0,4,                                    &
['BrO       ','HO2       '],                                                   &
['HOBr      ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50116,'B',.FALSE.,0,4,                                    &
['O(3P)     ','ClO       '],                                                   &
['Cl        ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50117,'B',.FALSE.,0,4,                                    &
['O(3P)     ','NO2       '],                                                   &
['NO        ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50118,'B',.FALSE.,0,4,                                    &
['O(3P)     ','HO2       '],                                                   &
['OH        ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50119,'B',.FALSE.,0,4,                                    &
['O3        ','H         '],                                                   &
['OH        ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50120,'B',.FALSE.,0,4,                                    &
['O(3P)     ','O3        '],                                                   &
['O2        ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50121,'J',.FALSE.,0,4,                                    &
['NO3       ','PHOTON    '],                                                   &
['NO        ','O2        ','          ','          ']),                        &
asad_flux_defn('RXN',50122,'B',.FALSE.,0,4,                                    &
['O(1D)     ','H2O       '],                                                   &
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50123,'B',.FALSE.,0,5,                                    &
['HO2       ','NO3       '],                                                   &
['OH        ','NO2       ','O2        ','          ']),                        &
asad_flux_defn('RXN',50124,'B',.FALSE.,0,4,                                    &
['OH        ','NO3       '],                                                   &
['HO2       ','NO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50125,'B',.FALSE.,0,5,                                    &
['NO3       ','HCHO      '],                                                   &
['HONO2     ','HO2       ','CO        ','          '])                         &
]

TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                                      asad_strat_o3_misc(15) = [               &
 asad_flux_defn('DEP',50131,'D',.FALSE.,0,1,                                   &
 ['O3        ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                                   &
 ['NO3       ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                                   &
 ['NO3       ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                                   &
 ['N2O5      ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                                   &
 ['N2O5      ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                                   &
 ['N2O5      ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                                   &
 ['HO2NO2    ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                                   &
 ['HONO2     ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                                   &
 ['NO3       ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                                   &
 ['NO3       ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                                   &
 ['N2O5      ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                                   &
 ['N2O5      ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                                   &
 ['N2O5      ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                                   &
 ['HO2NO2    ','          '],                                                  &
 ['          ','          ','          ','          ']),                       &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                                   &
 ['HONO2     ','          '],                                                  &
 ['          ','          ','          ','          '])                        &
 ]

! RC(O)O2 + NO2 --> *PAN production reactions, output to 50-248
TYPE(asad_flux_defn), PARAMETER :: asad_rco2no2_pan_prod(3) = [                &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T016
['MeCO3     ','NO2       '],                                                   &
['PAN       ','m         ','          ','          '] ),                       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T018
['EtCO3     ','NO2       '],                                                   &
['PPAN      ','m         ','          ','          '] ),                       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T020
['MACRO2    ','NO2       '],                                                   &
['MPAN      ','m         ','          ','          '] )                        &
  ]


! CRI-PAN forming reactions
TYPE(asad_flux_defn), PARAMETER :: cri_rco2no2_pan_prod(6) = [                 &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T016
['MeCO3     ','NO2       '],                                                   &
['PAN       ','m         ','          ','          '] ),                       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T018
['EtCO3     ','NO2       '],                                                   &
['PPAN      ','m         ','          ','          '] ),                       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T020
['HOCH2CO3  ','NO2       '],                                                   &
['PHAN      ','m         ','          ','          '] ),                       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T020
['RU12O2    ','NO2       '],                                                   &
['RU12PAN   ','m         ','          ','          '] ),                       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T020
['RU10O2    ','NO2       '],                                                   &
['MPAN      ','m         ','          ','          '] ),                       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T020
['RTN26O2   ','NO2       '],                                                   &
['RTN26PAN  ','m         ','          ','          '] )                        &
  ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_rco2no2_pan_prod(5) = [          &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T016
['MeCO3     ','NO2       '],                                                   &
['PAN       ','m         ','          ','          '] ),                       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T018
['EtCO3     ','NO2       '],                                                   &
['PPAN      ','m         ','          ','          '] ),                       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T020
['HOCH2CO3  ','NO2       '],                                                   &
['PHAN      ','m         ','          ','          '] ),                       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T020
['MACO3     ','NO2       '],                                                   &
['MPAN      ','m         ','          ','          '] ),                       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T020
['RTN26O2   ','NO2       '],                                                   &
['RTN26PAN  ','m         ','          ','          '] )                        &
]

! RO2 + HO2 type reactions, output to 50-249
TYPE(asad_flux_defn), PARAMETER :: asad_ro2ho2_reacn(13) = [                   &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B048
['HO2       ','EtCO3     '],                                                   &
['O2        ','EtCO3H    ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B049
['HO2       ','EtCO3     '],                                                   &
['O3        ','EtCO2H    ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B050
['HO2       ','EtOO      '],                                                   &
['EtOOH     ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B052
['HO2       ','ISO2      '],                                                   &
['ISOOH     ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B053
['HO2       ','MACRO2    '],                                                   &
['MACROOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B054
['HO2       ','MeCO3     '],                                                   &
['MeCO2H    ','O3        ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B055
['HO2       ','MeCO3     '],                                                   &
['MeCO3H    ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B056
['HO2       ','MeCO3     '],                                                   &
['OH        ','MeOO      ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B057
['HO2       ','MeCOCH2OO '],                                                   &
['MeCOCH2OOH','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B058
['HO2       ','MeOO      '],                                                   &
['HCHO      ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B059
['HO2       ','MeOO      '],                                                   &
['MeOOH     ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B063
['HO2       ','i-PrOO    '],                                                   &
['i-PrOOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['HO2       ','n-PrOO    '],                                                   &
['n-PrOOH   ','          ','          ','          '] )                        &
  ]


! RO2 + HO2 type reactions, output to 50-249
TYPE(asad_flux_defn), PARAMETER :: cri_ro2ho2_reacn(47) = [                    &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                                    &
['HO2       ','EtCO3     '],                                                   &
['O2        ','EtCO3H    ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B050
['HO2       ','EtOO      '],                                                   &
['EtOOH     ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B055
['HO2       ','MeCO3     '],                                                   &
['MeCO3H    ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B059
['HO2       ','MeOO      '],                                                   &
['MeOOH     ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B063
['HO2       ','i-PrOO    '],                                                   &
['i-PrOOH   ','          ','          ','          '] ),                       &
! CRI ROOH forming reactions from here
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN10O2    ','HO2       '],                                                   &
['RN10OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN13O2    ','HO2       '],                                                   &
['RN13OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN16O2    ','HO2       '],                                                   &
['RN16OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN19O2    ','HO2       '],                                                   &
['RN19OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN13AO2   ','HO2       '],                                                   &
['RN13OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN16AO2   ','HO2       '],                                                   &
['RN16OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RA13O2    ','HO2       '],                                                   &
['RA13OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RA16O2    ','HO2       '],                                                   &
['RA16OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RA19AO2   ','HO2       '],                                                   &
['RA19OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RA19CO2   ','HO2       '],                                                   &
['RA19OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['HOCH2CH2O2','HO2       '],                                                   &
['HOC2H4OOH ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN9O2     ','HO2       '],                                                   &
['RN9OOH    ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN12O2    ','HO2       '],                                                   &
['RN12OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN15O2    ','HO2       '],                                                   &
['RN15OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN18O2    ','HO2       '],                                                   &
['RN18OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN15AO2   ','HO2       '],                                                   &
['RN15OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN18AO2   ','HO2       '],                                                   &
['RN18OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['HOCH2CO3  ','HO2       '],                                                   &
['HOCH2CO3H ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN8O2     ','HO2       '],                                                   &
['RN8OOH    ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN11O2    ','HO2       '],                                                   &
['RN11OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN14O2    ','HO2       '],                                                   &
['RN14OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN17O2    ','HO2       '],                                                   &
['RN17OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RU14O2    ','HO2       '],                                                   &
['RU14OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RU12O2    ','HO2       '],                                                   &
['RU12OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RU10O2    ','HO2       '],                                                   &
['RU10OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRN6O2    ','HO2       '],                                                   &
['NRN6OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRN9O2    ','HO2       '],                                                   &
['NRN9OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRN12O2   ','HO2       '],                                                   &
['NRN12OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRU14O2   ','HO2       '],                                                   &
['NRU14OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRU12O2   ','HO2       '],                                                   &
['NRU12OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN28O2   ','HO2       '],                                                   &
['RTN28OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRTN28O2  ','HO2       '],                                                   &
['NRTN28OOH ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN26O2   ','HO2       '],                                                   &
['RTN26OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN25O2   ','HO2       '],                                                   &
['RTN25OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN24O2   ','HO2       '],                                                   &
['RTN24OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN23O2   ','HO2       '],                                                   &
['RTN23OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN14O2   ','HO2       '],                                                   &
['RTN14OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN10O2   ','HO2       '],                                                   &
['RTN10OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTX28O2   ','HO2       '],                                                   &
['RTX28OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTX24O2   ','HO2       '],                                                   &
['RTX24OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTX22O2   ','HO2       '],                                                   &
['RTX22OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRTX28O2  ','HO2       '],                                                   &
['NRTX28OOH ','          ','          ','          '] )                        &
  ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_ro2ho2_reacn(52) = [             &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,5,                                    &
['HO2       ','EtCO3     '],                                                   &
['OH        ','EtCO3H    ','EtOO      ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B050
['HO2       ','EtOO      '],                                                   &
['EtOOH     ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B055
['HO2       ','MeCO3     '],                                                   &
['MeCO3H    ','MeOO      ','OH        ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B055
['HO2       ','MeCO3     '],                                                   &
['OH        ','MeOO      ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B059
['HO2       ','MeOO      '],                                                   &
['MeOOH     ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B063
['HO2       ','i-PrOO    '],                                                   &
['i-PrOOH   ','          ','          ','          '] ),                       &
! CRI ROOH forming reactions from here
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN10O2    ','HO2       '],                                                   &
['RN10OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN13O2    ','HO2       '],                                                   &
['RN13OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN16O2    ','HO2       '],                                                   &
['RN16OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN19O2    ','HO2       '],                                                   &
['RN19OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN13AO2   ','HO2       '],                                                   &
['RN13OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN16AO2   ','HO2       '],                                                   &
['RN16OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RA13O2    ','HO2       '],                                                   &
['RA13OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RA16O2    ','HO2       '],                                                   &
['RA16OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RA19AO2   ','HO2       '],                                                   &
['RA19OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RA19CO2   ','HO2       '],                                                   &
['RA19OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['HOCH2CH2O2','HO2       '],                                                   &
['HOC2H4OOH ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN9O2     ','HO2       '],                                                   &
['RN9OOH    ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
['RN12O2    ','HO2       '],                                                   &
['RN12OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN15O2    ','HO2       '],                                                   &
['RN15OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN18O2    ','HO2       '],                                                   &
['RN18OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN15AO2   ','HO2       '],                                                   &
['RN15OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN18AO2   ','HO2       '],                                                   &
['RN18OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['HOCH2CO3  ','HO2       '],                                                   &
['HOCH2CO3H ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN8O2     ','HO2       '],                                                   &
['RN8OOH    ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN11O2    ','HO2       '],                                                   &
['RN11OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN14O2    ','HO2       '],                                                   &
['RN14OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RN17O2    ','HO2       '],                                                   &
['RN17OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RU14O2    ','HO2       '],                                                   &
['RU14OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RU12O2    ','HO2       '],                                                   &
['RU12OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RU10O2    ','HO2       '],                                                   &
['RU10OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRN6O2    ','HO2       '],                                                   &
['NRN6OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRN9O2    ','HO2       '],                                                   &
['NRN9OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRN12O2   ','HO2       '],                                                   &
['NRN12OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRU14O2   ','HO2       '],                                                   &
['NRU14OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRU12O2   ','HO2       '],                                                   &
['NRU12OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN28O2   ','HO2       '],                                                   &
['RTN28OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRTN28O2  ','HO2       '],                                                   &
['NRTN28OOH ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,5,                                    &
['RTN26O2   ','HO2       '],                                                   &
['RTN26OOH  ','RTN25O2   ','OH        ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN25O2   ','HO2       '],                                                   &
['RTN25OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN24O2   ','HO2       '],                                                   &
['RTN24OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN23O2   ','HO2       '],                                                   &
['RTN23OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN14O2   ','HO2       '],                                                   &
['RTN14OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTN10O2   ','HO2       '],                                                   &
['RTN10OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTX28O2   ','HO2       '],                                                   &
['RTX28OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTX24O2   ','HO2       '],                                                   &
['RTX24OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RTX22O2   ','HO2       '],                                                   &
['RTX22OOH  ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['NRTX28O2  ','HO2       '],                                                   &
['NRTX28OOH ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['DHPR12O2  ','HO2       '],                                                   &
['DHPR12OOH ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['MACO3     ','HO2       '],                                                   &
['RU10OOH   ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,6,                                    &
['MACO3     ','HO2       '],                                                   &
['MeOO      ','CO        ','HCHO      ','OH        '] ),                       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                                    &
['RU10AO2   ','HO2       '],                                                   &
['RU10OOH   ','          ','          ','          '] )                        &
]



! RO2 + NO3 type reactions, output to 50-250
TYPE(asad_flux_defn), PARAMETER :: asad_ro2no3_reacn(7) = [                    &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B039
['EtCO3     ','NO3       '],                                                   &
['EtOO      ','CO2       ','NO2       ','          '] ),                       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B042
['EtOO      ','NO3       '],                                                   &
['MeCHO     ','HO2       ','NO2       ','          '] ),                       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B072
['MeCO3     ','NO3       '],                                                   &
['MeOO      ','CO2       ','NO2       ','          '] ),                       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B074
['MeCOCH2OO ','NO3       '],                                                   &
['MeCO3     ','HCHO      ','NO2       ','          '] ),                       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B081
['MeOO      ','NO3       '],                                                   &
['HO2       ','HCHO      ','NO2       ','          '] ),                       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B195
['i-PrOO    ','NO3       '],                                                   &
['Me2CO     ','HO2       ','NO2       ','          '] ),                       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B197
['n-PrOO    ','NO3       '],                                                   &
['EtCHO     ','HO2       ','NO2       ','          '] )                        &
  ]


! CRI RO2 + NO3 type reactions, output to 50-250
TYPE(asad_flux_defn), PARAMETER :: cri_ro2no3_reacn(54) = [                    &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B039
  ['EtCO3     ','NO3       '],                                                 &
  ['EtOO      ','CO2       ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B042
  ['EtOO      ','NO3       '],                                                 &
  ['MeCHO     ','HO2       ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B072
  ['MeCO3     ','NO3       '],                                                 &
  ['MeOO      ','CO2       ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B081
  ['MeOO      ','NO3       '],                                                 &
  ['HO2       ','HCHO      ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B195
  ['i-PrOO    ','NO3       '],                                                 &
  ['Me2CO     ','HO2       ','NO2       ','          '] ),                     &
  !CRI RO2 + NO3 (49)
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                                  &
  ['RN10O2    ','NO3       '],                                                 &
  ['EtCHO     ','HO2       ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                                  &
  ['RN13O2    ','NO3       '],                                                 &
  ['MeCHO     ','EtOO      ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                                  &
  ['RN13O2    ','NO3       '],                                                 &
  ['CARB11A   ','HO2       ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,4,                                  &
  ['RN16O2    ','NO3       '],                                                 &
  ['RN15AO2   ','NO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,4,                                  &
  ['RN19O2    ','NO3       '],                                                 &
  ['RN18AO2   ','NO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,4,                                  &
  ['RN13AO2   ','NO3       '],                                                 &
  ['RN12O2    ','NO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,4,                                  &
  ['RN16AO2   ','NO3       '],                                                 &
  ['RN15O2    ','NO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,6,                                  &
  ['RA13O2    ','NO3       '],                                                 &
  ['CARB3     ','UDCARB8   ','HO2       ','NO2       '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,6,                                  &
  ['RA16O2    ','NO3       '],                                                 &
  ['CARB3     ','UDCARB11  ','HO2       ','NO2       '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,6,                                  &
  ['RA16O2    ','NO3       '],                                                 &
  ['CARB6     ','UDCARB8   ','HO2       ','NO2       '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,6,                                  &
  ['RA19AO2   ','NO3       '],                                                 &
  ['CARB3     ','UDCARB14  ','HO2       ','NO2       '] ),                     &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RA19CO2   ','NO3       '],                                                 &
  ['CARB9     ','UDCARB8   ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['HOCH2CH2O2','NO3       '],                                                 &
  ['HCHO      ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['HOCH2CH2O2','NO3       '],                                                 &
  ['HOCH2CHO  ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RN9O2     ','NO3       '],                                                 &
  ['MeCHO     ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RN12O2    ','NO3       '],                                                 &
  ['MeCHO     ','MeCHO     ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RN15O2    ','NO3       '],                                                 &
  ['EtCHO     ','MeCHO     ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RN18O2    ','NO3       '],                                                 &
  ['EtCHO     ','EtCHO     ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RN15AO2   ','NO3       '],                                                 &
  ['CARB13    ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RN18AO2   ','NO3       '],                                                 &
  ['CARB16    ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['HOCH2CO3  ','NO3       '],                                                 &
  ['HO2       ','HCHO      ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RN8O2     ','NO3       '],                                                 &
  ['MeCO3     ','HCHO      ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RN11O2    ','NO3       '],                                                 &
  ['MeCO3     ','MeCHO     ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RN14O2    ','NO3       '],                                                 &
  ['EtCO3     ','MeCHO     ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,4,                                   &
  ['RN17O2    ','NO3       '],                                                 &
  ['RN16AO2   ','NO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RU14O2    ','NO3       '],                                                 &
  ['UCARB12   ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RU14O2    ','NO3       '],                                                 &
  ['UCARB10   ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RU12O2    ','NO3       '],                                                 &
  ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RU12O2    ','NO3       '],                                                 &
  ['CARB7     ','CO        ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RU10O2    ','NO3       '],                                                 &
  ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RU10O2    ','NO3       '],                                                 &
  ['CARB6     ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RU10O2    ','NO3       '],                                                 &
  ['CARB7     ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['NRN6O2    ','NO3       '],                                                 &
  ['HCHO      ','HCHO      ','NO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['NRN9O2    ','NO3       '],                                                 &
  ['MeCHO     ','HCHO      ','NO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['NRN12O2   ','NO3       '],                                                 &
  ['MeCHO     ','MeCHO     ','NO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['NRU14O2   ','NO3       '],                                                 &
  ['NUCARB12  ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['NRU12O2   ','NO3       '],                                                 &
  ['NOA       ','CO        ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RTN28O2   ','NO3       '],                                                 &
  ['TNCARB26  ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['NRTN28O2  ','NO3       '],                                                 &
  ['TNCARB26  ','NO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,4,                                   &
  ['RTN26O2   ','NO3       '],                                                 &
  ['RTN25O2   ','NO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,4,                                   &
  ['RTN25O2   ','NO3       '],                                                 &
  ['RTN24O2   ','NO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,4,                                   &
  ['RTN24O2   ','NO3       '],                                                 &
  ['RTN23O2   ','NO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RTN23O2   ','NO3       '],                                                 &
  ['Me2CO     ','RTN14O2   ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RTN14O2   ','NO3       '],                                                 &
  ['HCHO      ','TNCARB10  ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RTN10O2   ','NO3       '],                                                 &
  ['RN8O2     ','CO        ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RTX28O2   ','NO3       '],                                                 &
  ['TXCARB24  ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RTX24O2   ','NO3       '],                                                 &
  ['TXCARB22  ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RTX22O2   ','NO3       '],                                                 &
  ['Me2CO     ','RN13O2    ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['NRTX28O2  ','NO3       '],                                                 &
  ['TXCARB24  ','HCHO      ','NO2       ','NO2       '])                       &
  ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_ro2no3_reacn(59) = [             &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B039
  ['EtCO3     ','NO3       '],                                                 &
  ['EtOO      ','CO2       ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B042
  ['EtOO      ','NO3       '],                                                 &
  ['MeCHO     ','HO2       ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B072
  ['MeCO3     ','NO3       '],                                                 &
  ['MeOO      ','CO2       ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B081
  ['MeOO      ','NO3       '],                                                 &
  ['HO2       ','HCHO      ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B195
  ['i-PrOO    ','NO3       '],                                                 &
  ['Me2CO     ','HO2       ','NO2       ','          '] ),                     &
  !CRI RO2 + NO3 (49)
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                                  &
  ['RN10O2    ','NO3       '],                                                 &
  ['EtCHO     ','HO2       ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                                  &
  ['RN13O2    ','NO3       '],                                                 &
  ['MeCHO     ','EtOO      ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                                  &
  ['RN13O2    ','NO3       '],                                                 &
  ['CARB11A   ','HO2       ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,4,                                  &
  ['RN16O2    ','NO3       '],                                                 &
  ['RN15AO2   ','NO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,4,                                  &
  ['RN19O2    ','NO3       '],                                                 &
  ['RN18AO2   ','NO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,4,                                  &
  ['RN13AO2   ','NO3       '],                                                 &
  ['RN12O2    ','NO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,4,                                  &
  ['RN16AO2   ','NO3       '],                                                 &
  ['RN15O2    ','NO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,6,                                  &
  ['RA13O2    ','NO3       '],                                                 &
  ['CARB3     ','UDCARB8   ','HO2       ','NO2       '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,6,                                  &
  ['RA16O2    ','NO3       '],                                                 &
  ['CARB3     ','UDCARB11  ','HO2       ','NO2       '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,6,                                  &
  ['RA16O2    ','NO3       '],                                                 &
  ['CARB6     ','UDCARB8   ','HO2       ','NO2       '] ),                     &
  asad_flux_defn('RXN',50250,'B',.FALSE.,0,6,                                  &
  ['RA19AO2   ','NO3       '],                                                 &
  ['CARB3     ','UDCARB14  ','HO2       ','NO2       '] ),                     &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RA19CO2   ','NO3       '],                                                 &
  ['CARB9     ','UDCARB8   ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['HOCH2CH2O2','NO3       '],                                                 &
  ['HCHO      ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['HOCH2CH2O2','NO3       '],                                                 &
  ['HOCH2CHO  ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RN9O2     ','NO3       '],                                                 &
  ['MeCHO     ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RN12O2    ','NO3       '],                                                 &
  ['MeCHO     ','MeCHO     ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RN15O2    ','NO3       '],                                                 &
  ['EtCHO     ','MeCHO     ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RN18O2    ','NO3       '],                                                 &
  ['EtCHO     ','EtCHO     ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RN15AO2   ','NO3       '],                                                 &
  ['CARB13    ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RN18AO2   ','NO3       '],                                                 &
  ['CARB16    ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['HOCH2CO3  ','NO3       '],                                                 &
  ['HO2       ','HCHO      ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RN8O2     ','NO3       '],                                                 &
  ['MeCO3     ','HCHO      ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RN11O2    ','NO3       '],                                                 &
  ['MeCO3     ','MeCHO     ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RN14O2    ','NO3       '],                                                 &
  ['EtCO3     ','MeCHO     ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,4,                                   &
  ['RN17O2    ','NO3       '],                                                 &
  ['RN16AO2   ','NO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RU14O2    ','NO3       '],                                                 &
  ['UCARB12   ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RU14O2    ','NO3       '],                                                 &
  ['UCARB10   ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RU12O2    ','NO3       '],                                                 &
  ['CARB6     ','HOCH2CHO  ','NO2       ','HO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RU12O2    ','NO3       '],                                                 &
  ['CARB7     ','CARB3     ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RU12O2    ','NO3       '],                                                 &
  ['CARB7     ','HOCH2CO3  ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RU10O2    ','NO3       '],                                                 &
  ['MeCO3     ','HOCH2CHO  ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RU10O2    ','NO3       '],                                                 &
  ['CARB6     ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RU10AO2   ','NO3       '],                                                 &
  ['CARB7     ','HO2       ','CO        ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['NRN6O2    ','NO3       '],                                                 &
  ['HCHO      ','HCHO      ','NO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['NRN9O2    ','NO3       '],                                                 &
  ['MeCHO     ','HCHO      ','NO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['NRN12O2   ','NO3       '],                                                 &
  ['MeCHO     ','MeCHO     ','NO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['NRU14O2   ','NO3       '],                                                 &
  ['NUCARB12  ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['NRU12O2   ','NO3       '],                                                 &
  ['NOA       ','CO        ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RTN28O2   ','NO3       '],                                                 &
  ['TNCARB26  ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['NRTN28O2  ','NO3       '],                                                 &
  ['TNCARB26  ','NO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,4,                                   &
  ['RTN26O2   ','NO3       '],                                                 &
  ['RTN25O2   ','NO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,4,                                   &
  ['RTN25O2   ','NO3       '],                                                 &
  ['RTN24O2   ','NO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,4,                                   &
  ['RTN24O2   ','NO3       '],                                                 &
  ['RTN23O2   ','NO2       ','          ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RTN23O2   ','NO3       '],                                                 &
  ['Me2CO     ','RTN14O2   ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RTN14O2   ','NO3       '],                                                 &
  ['HCHO      ','TNCARB10  ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RTN10O2   ','NO3       '],                                                 &
  ['RN8O2     ','CO        ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['RTX28O2   ','NO3       '],                                                 &
  ['TXCARB24  ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RTX24O2   ','NO3       '],                                                 &
  ['TXCARB22  ','HO2       ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,5,                                   &
  ['RTX22O2   ','NO3       '],                                                 &
  ['Me2CO     ','RN13O2    ','NO2       ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['NRTX28O2  ','NO3       '],                                                 &
  ['TXCARB24  ','HCHO      ','NO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['MACO3     ','NO3       '],                                                 &
  ['MeOO      ','HCHO      ','HO2       ','CO        ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,3,                                   &
  ['MACO3     ','NO3       '],                                                 &
  ['NO2       ','          ','          ','          ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['MACO3     ','NO3       '],                                                 &
  ['MeCO3     ','HCHO      ','HO2       ','NO2       ']),                      &
  asad_flux_defn('RXN',50250,'B',.TRUE.,0,6,                                   &
  ['DHPR12O2  ','NO3       '],                                                 &
  ['CARB3     ','RN8OOH    ','OH        ','NO2       '])                       &
  ]


! RO2 + RO2 type reactions, output to 50-251
TYPE(asad_flux_defn), PARAMETER :: asad_ro2ro2_reacn(8) = [                    &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B040
['EtOO      ','MeCO3     '],                                                   &
['MeCHO     ','HO2       ','MeOO      ','          '] ),                       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,6,                      & ! B065
['ISO2      ','ISO2      '],                                                   &
['MACR      ','MACR      ','HCHO      ','HO2       '] ),                       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,6,                      & ! B066
['MACRO2    ','MACRO2    '],                                                   &
['HACET     ','MGLY      ','HCHO      ','CO        '] ),                       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B067
['MACRO2    ','MACRO2    '],                                                   &
['HO2       ','          ','          ','          '] ),                       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B075
['MeOO      ','MeCO3     '],                                                   &
['HO2       ','HCHO      ','MeOO      ','          '] ),                       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B076
['MeOO      ','MeCO3     '],                                                   &
['MeCO2H    ','HCHO      ','          ','          '] ),                       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,6,                      & ! B077
['MeOO      ','MeOO      '],                                                   &
['HO2       ','HO2       ','HCHO      ','HCHO      '] ),                       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B078
['MeOO      ','MeOO      '],                                                   &
['MeOH      ','HCHO      ','          ','          '] )                        &
  ]

! Make RO2+RO2 for StratTrop with RO2-permuation chem
! RO2 + RO2 type reactions, output to 50-251
TYPE(asad_flux_defn), PARAMETER :: asad_ro2perm_ro2ro2_reacn(15) = [           &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B278
  ['MeOO      ','RO2       '],                                                 &
  ['HO2       ','HCHO      ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B279
  ['MeOO      ','RO2       '],                                                 &
  ['MeOH      ','HCHO      ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B300
  ['EtOO      ','RO2       '],                                                 &
  ['HO2       ','MeCHO     ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B301
  ['EtOO      ','RO2       '],                                                 &
  ['MeCHO     ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B302
  ['n-PrOO    ','RO2       '],                                                 &
  ['HO2       ','EtCHO     ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B303
  ['n-PrOO    ','RO2       '],                                                 &
  ['EtCHO     ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B304
  ['i-PrOO    ','RO2       '],                                                 &
  ['HO2       ','Me2CO     ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B305
  ['i-PrOO    ','RO2       '],                                                 &
  ['Me2CO     ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B306
  ['MeCO3     ','RO2       '],                                                 &
  ['MeCO2H    ','MeOO      ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B307
  ['EtCO3     ','RO2       '],                                                 &
  ['EtCO2H    ','EtOO      ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B308
  ['MeCOCH2OO ','RO2       '],                                                 &
  ['HCHO      ','MeCO3     ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B309
  ['MeCOCH2OO ','RO2       '],                                                 &
  ['HACET     ','MGLY      ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B310
  ['ISO2      ','RO2       '],                                                 &
  ['MACR      ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,6,                      & ! B311
  ['MACRO2    ','RO2       '],                                                 &
  ['HACET     ','MGLY      ','HCHO      ','CO        '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B312
  ['MACRO2    ','RO2       '],                                                 &
  ['HO2       ','          ','          ','          '] )                      &
  ]


!--> Make RO2+RO2 for CRI chem
TYPE(asad_flux_defn), PARAMETER :: cri_ro2ro2_reacn(63) = [                    &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B504
  ['MeOO      ','RO2       '],                                                 &
  ['HCHO      ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B505
  ['MeOO      ','RO2       '],                                                 &
  ['HCHO      ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B506
  ['MeOO      ','RO2       '],                                                 &
  ['MeOH      ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B507
  ['EtOO      ','RO2       '],                                                 &
  ['MeCHO     ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B508
  ['EtOO      ','RO2       '],                                                 &
  ['MeCHO     ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B509
  ['EtOO      ','RO2       '],                                                 &
  ['EtOH      ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B510
  ['RN10O2    ','RO2       '],                                                 &
  ['EtCHO     ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B511
  ['RN10O2    ','RO2       '],                                                 &
  ['EtCHO     ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B512
  ['RN10O2    ','RO2       '],                                                 &
  ['n-PrOH    ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B513
  ['i-PrOO    ','RO2       '],                                                 &
  ['Me2CO     ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B514
  ['i-PrOO    ','RO2       '],                                                 &
  ['Me2CO     ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B515
  ['i-PrOO    ','RO2       '],                                                 &
  ['i-PrOH    ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B516
  ['RN13O2    ','RO2       '],                                                 &
  ['MeCHO     ','EtOO      ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B517
  ['RN13O2    ','RO2       '],                                                 &
  ['CARB11A   ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B518
  ['RN13AO2   ','RO2       '],                                                 &
  ['RN12O2    ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B519
  ['RN16AO2   ','RO2       '],                                                 &
  ['RN15O2    ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B520
  ['RA13O2    ','RO2       '],                                                 &
  ['CARB3     ','UDCARB8   ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B521
  ['RA16O2    ','RO2       '],                                                 &
  ['CARB3     ','UDCARB11  ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B522
  ['RA16O2    ','RO2       '],                                                 &
  ['CARB6     ','UDCARB8   ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B523
  ['RA19AO2   ','RO2       '],                                                 &
  ['CARB3     ','UDCARB14  ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B524
  ['RA19CO2   ','RO2       '],                                                 &
  ['CARB3     ','UDCARB14  ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B525
  ['RN16O2    ','RO2       '],                                                 &
  ['RN15AO2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B526
  ['RN19O2    ','RO2       '],                                                 &
  ['RN18AO2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B527
  ['HOCH2CH2O2','RO2       '],                                                 &
  ['HCHO      ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B528
  ['HOCH2CH2O2','RO2       '],                                                 &
  ['HOCH2CHO  ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B529
  ['RN9O2     ','RO2       '],                                                 &
  ['MeCHO     ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B530
  ['RN12O2    ','RO2       '],                                                 &
  ['MeCHO     ','MeCHO     ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B531
  ['RN15O2    ','RO2       '],                                                 &
  ['EtCHO     ','MeCHO     ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B532
  ['RN18O2    ','RO2       '],                                                 &
  ['EtCHO     ','EtCHO     ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B533
  ['RN15AO2   ','RO2       '],                                                 &
  ['CARB13    ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B534
  ['RN18AO2   ','RO2       '],                                                 &
  ['CARB16    ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B535
  ['MeCO3     ','RO2       '],                                                 &
  ['MeOO      ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B536
  ['EtCO3     ','RO2       '],                                                 &
  ['EtOO      ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B537
  ['HOCH2CO3  ','RO2       '],                                                 &
  ['HCHO      ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B538
  ['RN8O2     ','RO2       '],                                                 &
  ['MeCO3     ','HCHO      ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B539
  ['RN11O2    ','RO2       '],                                                 &
  ['MeCO3     ','MeCHO     ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B540
  ['RN14O2    ','RO2       '],                                                 &
  ['EtCO3     ','MeCHO     ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B541
  ['RN17O2    ','RO2       '],                                                 &
  ['RN16AO2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B542
  ['RU14O2    ','RO2       '],                                                 &
  ['UCARB12   ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B543
  ['RU14O2    ','RO2       '],                                                 &
  ['UCARB10   ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B544
  ['RU12O2    ','RO2       '],                                                 &
  ['MeCO3     ','HOCH2CHO  ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B545
  ['RU12O2    ','RO2       '],                                                 &
  ['CARB7     ','HOCH2CHO  ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B546
  ['RU10O2    ','RO2       '],                                                 &
  ['MeCO3     ','HOCH2CHO  ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B547
  ['RU10O2    ','RO2       '],                                                 &
  ['CARB6     ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B548
  ['RU10O2    ','RO2       '],                                                 &
  ['CARB7     ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B549
  ['NRN6O2    ','RO2       '],                                                 &
  ['HCHO      ','HCHO      ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B550
  ['NRN9O2    ','RO2       '],                                                 &
  ['MeCHO     ','HCHO      ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B551
  ['NRN12O2   ','RO2       '],                                                 &
  ['MeCHO     ','MeCHO     ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B552
  ['NRU14O2   ','RO2       '],                                                 &
  ['NUCARB12  ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B553
  ['NRU12O2   ','RO2       '],                                                 &
  ['NOA       ','CO        ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B554
  ['RTN28O2   ','RO2       '],                                                 &
  ['TNCARB26  ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B555
  ['NRTN28O2  ','RO2       '],                                                 &
  ['TNCARB26  ','NO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B556
  ['RTN26O2   ','RO2       '],                                                 &
  ['RTN25O2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B557
  ['RTN25O2   ','RO2       '],                                                 &
  ['RTN24O2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B558
  ['RTN24O2   ','RO2       '],                                                 &
  ['RTN23O2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B559
  ['RTN23O2   ','RO2       '],                                                 &
  ['Me2CO     ','RTN14O2   ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B560
  ['RTN14O2   ','RO2       '],                                                 &
  ['HCHO      ','TNCARB10  ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B561
  ['RTN10O2   ','RO2       '],                                                 &
  ['RN8O2     ','CO        ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B562
  ['RTX28O2   ','RO2       '],                                                 &
  ['TXCARB24  ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B563
  ['RTX24O2   ','RO2       '],                                                 &
  ['TXCARB22  ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B564
  ['RTX22O2   ','RO2       '],                                                 &
  ['Me2CO     ','RN13O2    ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B565
  ['NRTX28O2  ','RO2       '],                                                 &
  ['TXCARB24  ','HCHO      ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B566
  ['NRTX28O2  ','RO2       '],                                                 &
  ['TXCARB24  ','HCHO      ','NO2       ','          '] )                      &
  ]

TYPE(asad_flux_defn), PARAMETER :: cri_strat2_ro2ro2_reacn(66) = [             &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B504
  ['MeOO      ','RO2       '],                                                 &
  ['HCHO      ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B505
  ['MeOO      ','RO2       '],                                                 &
  ['HCHO      ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B506
  ['MeOO      ','RO2       '],                                                 &
  ['MeOH      ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B507
  ['EtOO      ','RO2       '],                                                 &
  ['MeCHO     ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B508
  ['EtOO      ','RO2       '],                                                 &
  ['MeCHO     ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B509
  ['EtOO      ','RO2       '],                                                 &
  ['EtOH      ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B510
  ['RN10O2    ','RO2       '],                                                 &
  ['EtCHO     ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B511
  ['RN10O2    ','RO2       '],                                                 &
  ['EtCHO     ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B512
  ['RN10O2    ','RO2       '],                                                 &
  ['n-PrOH    ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B513
  ['i-PrOO    ','RO2       '],                                                 &
  ['Me2CO     ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B514
  ['i-PrOO    ','RO2       '],                                                 &
  ['Me2CO     ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B515
  ['i-PrOO    ','RO2       '],                                                 &
  ['i-PrOH    ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B516
  ['RN13O2    ','RO2       '],                                                 &
  ['MeCHO     ','EtOO      ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B517
  ['RN13O2    ','RO2       '],                                                 &
  ['CARB11A   ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B518
  ['RN13AO2   ','RO2       '],                                                 &
  ['RN12O2    ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B519
  ['RN16AO2   ','RO2       '],                                                 &
  ['RN15O2    ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B520
  ['RA13O2    ','RO2       '],                                                 &
  ['CARB3     ','UDCARB8   ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B521
  ['RA16O2    ','RO2       '],                                                 &
  ['CARB3     ','UDCARB11  ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B522
  ['RA16O2    ','RO2       '],                                                 &
  ['CARB6     ','UDCARB8   ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B523
  ['RA19AO2   ','RO2       '],                                                 &
  ['CARB3     ','UDCARB14  ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B524
  ['RA19CO2   ','RO2       '],                                                 &
  ['CARB3     ','UDCARB14  ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B525
  ['RN16O2    ','RO2       '],                                                 &
  ['RN15AO2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B526
  ['RN19O2    ','RO2       '],                                                 &
  ['RN18AO2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B527
  ['HOCH2CH2O2','RO2       '],                                                 &
  ['HCHO      ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B528
  ['HOCH2CH2O2','RO2       '],                                                 &
  ['HOCH2CHO  ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B529
  ['RN9O2     ','RO2       '],                                                 &
  ['MeCHO     ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B530
  ['RN12O2    ','RO2       '],                                                 &
  ['MeCHO     ','MeCHO     ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B531
  ['RN15O2    ','RO2       '],                                                 &
  ['EtCHO     ','MeCHO     ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B532
  ['RN18O2    ','RO2       '],                                                 &
  ['EtCHO     ','EtCHO     ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B533
  ['RN15AO2   ','RO2       '],                                                 &
  ['CARB13    ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B534
  ['RN18AO2   ','RO2       '],                                                 &
  ['CARB16    ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B535
  ['MeCO3     ','RO2       '],                                                 &
  ['MeOO      ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B536
  ['EtCO3     ','RO2       '],                                                 &
  ['EtOO      ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B537
  ['HOCH2CO3  ','RO2       '],                                                 &
  ['HCHO      ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B538
  ['RN8O2     ','RO2       '],                                                 &
  ['MeCO3     ','HCHO      ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B539
  ['RN11O2    ','RO2       '],                                                 &
  ['MeCO3     ','MeCHO     ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B540
  ['RN14O2    ','RO2       '],                                                 &
  ['EtCO3     ','MeCHO     ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B541
  ['RN17O2    ','RO2       '],                                                 &
  ['RN16AO2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B542
  ['RU14O2    ','RO2       '],                                                 &
  ['UCARB12   ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B543
  ['RU14O2    ','RO2       '],                                                 &
  ['UCARB10   ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B544
  ['RU12O2    ','RO2       '],                                                 &
  ['MeCO3     ','HOCH2CHO  ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B545
  ['RU12O2    ','RO2       '],                                                 &
  ['CARB7     ','HOCH2CHO  ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B546
  ['RU10O2    ','RO2       '],                                                 &
  ['MeCO3     ','HOCH2CHO  ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B547
  ['RU10O2    ','RO2       '],                                                 &
  ['CARB6     ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B547
  ['RU10AO2   ','RO2       '],                                                 &
  ['CARB7     ','CO        ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B549
  ['NRN6O2    ','RO2       '],                                                 &
  ['HCHO      ','HCHO      ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B550
  ['NRN9O2    ','RO2       '],                                                 &
  ['MeCHO     ','HCHO      ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B551
  ['NRN12O2   ','RO2       '],                                                 &
  ['MeCHO     ','MeCHO     ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B552
  ['NRU14O2   ','RO2       '],                                                 &
  ['NUCARB12  ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B553
  ['NRU12O2   ','RO2       '],                                                 &
  ['NOA       ','CO        ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B554
  ['RTN28O2   ','RO2       '],                                                 &
  ['TNCARB26  ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B555
  ['NRTN28O2  ','RO2       '],                                                 &
  ['TNCARB26  ','NO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B556
  ['RTN26O2   ','RO2       '],                                                 &
  ['RTN25O2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B557
  ['RTN25O2   ','RO2       '],                                                 &
  ['RTN24O2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B558
  ['RTN24O2   ','RO2       '],                                                 &
  ['RTN23O2   ','          ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B559
  ['RTN23O2   ','RO2       '],                                                 &
  ['Me2CO     ','RTN14O2   ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B560
  ['RTN14O2   ','RO2       '],                                                 &
  ['HCHO      ','TNCARB10  ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B561
  ['RTN10O2   ','RO2       '],                                                 &
  ['RN8O2     ','CO        ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B562
  ['RTX28O2   ','RO2       '],                                                 &
  ['TXCARB24  ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B563
  ['RTX24O2   ','RO2       '],                                                 &
  ['TXCARB22  ','HO2       ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B564
  ['RTX22O2   ','RO2       '],                                                 &
  ['Me2CO     ','RN13O2    ','          ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B565
  ['NRTX28O2  ','RO2       '],                                                 &
  ['TXCARB24  ','HCHO      ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B566
  ['NRTX28O2  ','RO2       '],                                                 &
  ['TXCARB24  ','HCHO      ','NO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,6,                      & ! B565
  ['MACO3     ','RO2       '],                                                 &
  ['MeOO      ','CO        ','HCHO      ','OH        '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B566
  ['MACO3     ','RO2       '],                                                 &
  ['MeCO3     ','HCHO      ','HO2       ','          '] ),                     &
  asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B566
  ['DHPR12O2  ','RO2       '],                                                 &
  ['CARB3     ','RN8OOH    ','OH        ','          '] )                      &
  ]

! Methane Oxidation reactions, output to 50-247
TYPE(asad_flux_defn), PARAMETER :: asad_ch4_oxidn(6) = [                       &
asad_flux_defn('RXN',50247,'J',.FALSE.,0,4,                      & ! RATJ_T
['CH4       ','PHOTON    '],                                                   &
['MeOO      ','H         ','          ','          ']),                        &
asad_flux_defn('RXN',50247,'B',.FALSE.,0,4,                      & ! B016
['Cl        ','CH4       '],                                                   &
['HCl       ','MeOO      ','          ','          '] ),                       &
asad_flux_defn('RXN',50247,'B',.FALSE.,0,4,                      & ! B101
['O(1D)     ','CH4       '],                                                   &
['HCHO      ','H2        ','          ','          '] ),                       &
asad_flux_defn('RXN',50247,'B',.FALSE.,0,5,                      & ! B102
['O(1D)     ','CH4       '],                                                   &
['HCHO      ','HO2       ','HO2       ','          '] ),                       &
asad_flux_defn('RXN',50247,'B',.FALSE.,0,4,                      & ! B103
['O(1D)     ','CH4       '],                                                   &
['OH        ','MeOO      ','          ','          '] ),                       &
asad_flux_defn('RXN',50247,'B',.FALSE.,0,4,                      & ! B145
['OH        ','CH4       '],                                                   &
['H2O       ','MeOO      ','          ','          '] )                        &
  ]

! Chemical Production of O(1D), output to 50-254
TYPE(asad_flux_defn), PARAMETER :: asad_o1d_prod(3) = [                        &
asad_flux_defn('RXN',50254,'J',.FALSE.,0,4,                                    &
['O3        ','PHOTON    '],                                                   &
['O2        ','O(1D)     ','          ','          ']),                        &
asad_flux_defn('RXN',50254,'J',.FALSE.,0,4,                                    &
['O2        ','PHOTON    '],                                                   &
['O(3P)     ','O(1D)     ','          ','          ']),                        &
asad_flux_defn('RXN',50254,'J',.FALSE.,0,4,                                    &
['N2O       ','PHOTON    '],                                                   &
['N2        ','O(1D)     ','          ','          '])                         &
  ]

! Tropospheric sulphur chemistry for online oxidants
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                            asad_aerosol_chem_online(16) = [                   &
asad_flux_defn('RXN',50140,'B',.FALSE.,0,5,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','MeOO      ','HCHO      ','          ']),                        &
asad_flux_defn('RXN',50141,'B',.FALSE.,0,5,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','DMSO      ','MeOO      ','          ']),                        &
asad_flux_defn('RXN',50142,'B',.FALSE.,0,6,                                    &
['DMS       ','NO3       '],                                                   &
['SO2       ','HONO2     ','MeOO      ','HCHO      ']),                        &
asad_flux_defn('RXN',50143,'B',.FALSE.,0,4,                                    &
['DMSO      ','OH        '],                                                   &
['SO2       ','MSA       ','          ','          ']),                        &
asad_flux_defn('RXN',50144,'B',.FALSE.,0,4,                                    &
['CS2       ','OH        '],                                                   &
['SO2       ','COS       ','          ','          ']),                        &
asad_flux_defn('RXN',50145,'B',.FALSE.,0,3,                                    &
['H2S       ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50146,'B',.FALSE.,0,3,                                    &
['COS       ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50147,'B',.FALSE.,0,3,                                    &
['Monoterp  ','OH        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,3,                                    &
['Monoterp  ','O3        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50149,'B',.FALSE.,0,3,                                    &
['Monoterp  ','NO3       '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50150,'T',.FALSE.,0,4,                                    &
['SO2       ','OH        '],                                                   &
['HO2       ','H2SO4     ','          ','          ']),                        &
asad_flux_defn('RXN',50151,'H',.FALSE.,0,3,                                    &
['SO2       ','H2O2      '],                                                   &
['NULL0     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50152,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL1     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50153,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL2     ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50154,'D',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50155,'W',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          '])                         &
]

! Tropospheric sulphur chemistry for Offline oxidants chemistry
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                                   asad_aerosol_chem_offline(16)               &
                                                         = [                   &
asad_flux_defn('RXN',50140,'B',.FALSE.,0,3,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50141,'B',.FALSE.,0,4,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','DMSO      ','          ','          ']),                        &
asad_flux_defn('RXN',50142,'B',.FALSE.,0,3,                                    &
['DMS       ','NO3       '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50143,'B',.FALSE.,0,3,                                    &
['DMSO      ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50144,'B',.FALSE.,0,4,                                    &
['CS2       ','OH        '],                                                   &
['SO2       ','COS       ','          ','          ']),                        &
asad_flux_defn('RXN',50145,'B',.FALSE.,0,3,                                    &
['H2S       ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50146,'B',.FALSE.,0,3,                                    &
['COS       ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50147,'B',.FALSE.,0,3,                                    &
['Monoterp  ','OH        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,3,                                    &
['Monoterp  ','O3        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50149,'B',.FALSE.,0,3,                                    &
['Monoterp  ','NO3       '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50150,'T',.FALSE.,0,3,                                    &
['SO2       ','OH        '],                                                   &
['H2SO4     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50151,'H',.FALSE.,0,3,                                    &
['SO2       ','H2O2      '],                                                   &
['NULL0     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50152,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL1     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50153,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL2     ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50154,'D',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50155,'W',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          '])                         &
]

! Water production: sum into section 50, item 238
! Reaction identifiers from ukca_chem_strattrop are shown
TYPE(asad_flux_defn), PARAMETER :: asad_h2o_budget(38) = [                     &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['H         ','HO2       '],                           & ! B044
['O(3P)     ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['MeBr      ','OH        '],                           & ! B070
['Br        ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','C2H6      '],                           & ! B141
['H2O       ','EtOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','C3H8      '],                           & ! B142
['i-PrOO    ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','C3H8      '],                           & ! B143
['n-PrOO    ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','CH4       '],                           & ! B145
['H2O       ','MeOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','EtCHO     '],                           & ! B150
['H2O       ','EtCO3     ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','EtOOH     '],                           & ! B151
['H2O       ','EtOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','EtOOH     '],                           & ! B152
['H2O       ','MeCHO     ','OH        ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','H2        '],                           & ! B153
['H2O       ','HO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','H2O2      '],                           & ! B154
['H2O       ','HO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','HBr       '],                           & ! B156
['H2O       ','Br        ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','HCHO      '],                           & ! B157
['H2O       ','HO2       ','CO        ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','HCl       '],                           & ! B159
['H2O       ','Cl        ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,3,                                    &
['OH        ','HO2       '],                           & ! B160
['H2O       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','HONO2     '],                           & ! B164
['H2O       ','NO3       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','HO2NO2    '],                           & ! B161
['H2O       ','NO2       ','O2        ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','HOCl      '],                           & ! B162
['ClO       ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','HONO      '],                           & ! B163
['H2O       ','NO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,1,4,                                    &
['OH        ','Me2CO     '],                           & ! B172
['H2O       ','MeCOCH2OO ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,2,4,                                    &
['OH        ','Me2CO     '],                           & ! B173
['H2O       ','MeCOCH2OO ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','MeCHO     '],                           & ! B174
['H2O       ','MeCO3     ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','MeCOCH2OOH'],                           & ! B177
['H2O       ','MeCOCH2OO ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','MeONO2    '],                           & ! B180
['HCHO      ','NO2       ','H2O       ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','MeOOH     '],                           & ! B181
['H2O       ','HCHO      ','OH        ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','MeOOH     '],                           & ! B182
['H2O       ','MeOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','OH        '],                           & ! B187
['H2O       ','O(3P)     ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','PAN       '],                           & ! B188
['HCHO      ','NO2       ','H2O       ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','PPAN      '],                           & ! B189
['MeCHO     ','NO2       ','H2O       ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','i-PrOOH   '],                           & ! B191
['i-PrOO    ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','n-PrOOH   '],                           & ! B192
['n-PrOO    ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','n-PrOOH   '],                           & ! B193
['EtCHO     ','H2O       ','OH        ','          ']),                        &
asad_flux_defn('RXN',50238,'H',.FALSE.,0,5,                                    &
['HOCl      ','HCl       '],                           & ! Het PSC
['Cl        ','Cl        ','H2O       ','          ']),                        &

! Water loss: sum into section 50, item 239
asad_flux_defn('RXN',50239,'B',.FALSE.,0,4,                                    &
['O(1D)     ','H2O       '],                            & ! B106
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50239,'B',.FALSE.,0,4,                                    &
['N2O5      ','H2O       '],                            & ! B085
['HONO2     ','HONO2     ','          ','          ']),                        &
asad_flux_defn('RXN',50239,'J',.FALSE.,0,4,                                    &
['H2O       ','PHOTON    '],                            & ! Photol
['OH        ','H         ','          ','          ']),                        &
asad_flux_defn('RXN',50239,'H',.FALSE.,0,4,                                    &
['ClONO2    ','H2O       '],                            & ! Het PSC
['HOCl      ','HONO2     ','          ','          ']),                        &
asad_flux_defn('RXN',50239,'H',.FALSE.,0,4,                                    &
['N2O5      ','H2O       '],                            & ! Het PSC
['HONO2     ','HONO2     ','          ','          '])                         &
]


! Water production: sum into section 50, item 238
! Reaction identifiers from ukca_chem_strattrop are shown
TYPE(asad_flux_defn), PARAMETER :: asad_h2o_budget_121(39) = [                 &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['H         ','HO2       '],                           & ! B044
['O(3P)     ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['MeBr      ','OH        '],                           & ! B070
['Br        ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','C2H6      '],                           & ! B141
['H2O       ','EtOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','C3H8      '],                           & ! B142
['i-PrOO    ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','C3H8      '],                           & ! B143
['n-PrOO    ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','CH4       '],                           & ! B145
['H2O       ','MeOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','EtCHO     '],                           & ! B150
['H2O       ','EtCO3     ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','EtOOH     '],                           & ! B151
['H2O       ','EtOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','EtOOH     '],                           & ! B152
['H2O       ','MeCHO     ','OH        ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','H2        '],                           & ! B153
['H2O       ','H         ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','H2O2      '],                           & ! B154
['H2O       ','HO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','HBr       '],                           & ! B156
['H2O       ','Br        ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','HCHO      '],                           & ! B157
['H2O       ','HO2       ','CO        ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','HCOOH     '],                           & ! B158
['CO2       ','H2O       ','H         ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','HCl       '],                           & ! B159
['H2O       ','Cl        ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,3,                                    &
['OH        ','HO2       '],                           & ! B160
['H2O       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','HONO2     '],                           & ! B164
['H2O       ','NO3       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','HO2NO2    '],                           & ! B161
['H2O       ','NO2       ','O2        ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','HOCl      '],                           & ! B162
['ClO       ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','HONO      '],                           & ! B163
['H2O       ','NO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,1,4,                                    &
['OH        ','Me2CO     '],                           & ! B172
['H2O       ','MeCOCH2OO ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,2,4,                                    &
['OH        ','Me2CO     '],                           & ! B173
['H2O       ','MeCOCH2OO ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','MeCHO     '],                           & ! B174
['H2O       ','MeCO3     ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','MeCOCH2OOH'],                           & ! B177
['H2O       ','MeCOCH2OO ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','MeONO2    '],                           & ! B180
['HCHO      ','NO2       ','H2O       ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','MeOOH     '],                           & ! B181
['H2O       ','HCHO      ','OH        ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','MeOOH     '],                           & ! B182
['H2O       ','MeOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','OH        '],                           & ! B187
['H2O       ','O(3P)     ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','PAN       '],                           & ! B188
['HCHO      ','NO2       ','H2O       ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','PPAN      '],                           & ! B189
['MeCHO     ','NO2       ','H2O       ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','i-PrOOH   '],                           & ! B191
['i-PrOO    ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                                    &
['OH        ','n-PrOOH   '],                           & ! B192
['n-PrOO    ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                                    &
['OH        ','n-PrOOH   '],                           & ! B193
['EtCHO     ','H2O       ','OH        ','          ']),                        &
asad_flux_defn('RXN',50238,'H',.FALSE.,0,5,                                    &
['HOCl      ','HCl       '],                           & ! Het PSC
['Cl        ','Cl        ','H2O       ','          ']),                        &

! Water loss: sum into section 50, item 239
asad_flux_defn('RXN',50239,'B',.FALSE.,0,4,                                    &
['O(1D)     ','H2O       '],                            & ! B106
['OH        ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50239,'B',.FALSE.,0,4,                                    &
['N2O5      ','H2O       '],                            & ! B085
['HONO2     ','HONO2     ','          ','          ']),                        &
asad_flux_defn('RXN',50239,'J',.FALSE.,0,4,                                    &
['H2O       ','PHOTON    '],                            & ! Photol
['OH        ','H         ','          ','          ']),                        &
asad_flux_defn('RXN',50239,'H',.FALSE.,0,4,                                    &
['ClONO2    ','H2O       '],                            & ! Het PSC
['HOCl      ','HONO2     ','          ','          ']),                        &
asad_flux_defn('RXN',50239,'H',.FALSE.,0,4,                                    &
['N2O5      ','H2O       '],                            & ! Het PSC
['HONO2     ','HONO2     ','          ','          '])                         &
]


! Strat-Trop sulphur chemistry (contains explicit SO3)
! For i_ukca_chem_version < 117
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                         asad_aerosol_chem_strattrop(22) = [                   &
asad_flux_defn('RXN',50140,'B',.FALSE.,0,5,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50141,'B',.FALSE.,0,5,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','MSA       ','          ','          ']),                        &
asad_flux_defn('RXN',50142,'B',.FALSE.,0,6,                                    &
['DMS       ','NO3       '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50144,'B',.FALSE.,0,4,                                    &
['CS2       ','OH        '],                                                   &
['COS       ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50145,'B',.FALSE.,0,4,                                    &
['H2S       ','OH        '],                                                   &
['SO2       ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50146,'B',.FALSE.,0,4,                                    &
['COS       ','OH        '],                                                   &
['CO2       ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50147,'B',.FALSE.,0,3,                                    &
['Monoterp  ','OH        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,3,                                    &
['Monoterp  ','O3        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50149,'B',.FALSE.,0,3,                                    &
['Monoterp  ','NO3       '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50150,'T',.FALSE.,0,4,                                    &
['SO2       ','OH        '],                                                   &
['HO2       ','SO3       ','          ','          ']),                        &
asad_flux_defn('RXN',50151,'H',.FALSE.,0,3,                                    &
['SO2       ','H2O2      '],                                                   &
['NULL0     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50152,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL1     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50153,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL2     ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50154,'D',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50155,'W',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50331,'J',.TRUE.,0,4,                                     &
['COS       ','PHOTON    '],                                                   &
['CO        ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50332,'J',.TRUE.,0,4,                                     &
['H2SO4     ','PHOTON    '],                                                   &
['SO3       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50333,'J',.TRUE.,0,4,                                     &
['SO3       ','PHOTON    '],                                                   &
['SO2       ','O(3P)     ','          ','          ']),                        &
asad_flux_defn('RXN',50334,'B',.TRUE.,0,3,                                     &
['SO2       ','O3        '],                                                   &
['SO3       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50335,'B',.TRUE.,0,3,                                     &
['DMS       ','O(3P)     '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50336,'B',.TRUE.,0,4,                                     &
['COS       ','O(3P)     '],                                                   &
['CO        ','SO2       ','          ','          ']),                        &
! Retain H2O here as it is in original scheme
asad_flux_defn('RXN',50337,'B',.TRUE.,0,4,                                     &
['SO3       ','H2O       '],                                                   &
['H2SO4     ','H2O       ','          ','          '])                         &
]

! Strat-Trop sulphur chemistry (contains explicit SO3) for
!                                         i_ukca_chem_version >=117
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                         asad_aerosol_chem_strattrop_117(23) = [               &
asad_flux_defn('RXN',50140,'B',.FALSE.,0,3,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50142,'B',.FALSE.,0,3,                                    &
['DMS       ','NO3       '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50143,'B',.FALSE.,0,4,                                    &
['DMSO      ','OH        '],                                                   &
['SO2       ','MSA       ','          ','          ']),                        &
asad_flux_defn('RXN',50144,'B',.FALSE.,0,4,                                    &
['CS2       ','OH        '],                                                   &
['COS       ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50145,'B',.FALSE.,0,4,                                    &
['H2S       ','OH        '],                                                   &
['SO2       ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50146,'B',.FALSE.,0,4,                                    &
['COS       ','OH        '],                                                   &
['CO2       ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50147,'B',.FALSE.,0,3,                                    &
['Monoterp  ','OH        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,3,                                    &
['Monoterp  ','O3        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50149,'B',.FALSE.,0,3,                                    &
['Monoterp  ','NO3       '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50150,'T',.FALSE.,0,4,                                    &
['SO2       ','OH        '],                                                   &
['HO2       ','SO3       ','          ','          ']),                        &
asad_flux_defn('RXN',50151,'H',.FALSE.,0,3,                                    &
['SO2       ','H2O2      '],                                                   &
['NULL0     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50152,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL1     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50153,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL2     ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50154,'D',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50155,'W',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50331,'J',.FALSE.,0,4,                                    &
['COS       ','PHOTON    '],                                                   &
['CO        ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50332,'J',.FALSE.,0,4,                                    &
['H2SO4     ','PHOTON    '],                                                   &
['SO3       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50333,'J',.FALSE.,0,4,                                    &
['SO3       ','PHOTON    '],                                                   &
['SO2       ','O(3P)     ','          ','          ']),                        &
asad_flux_defn('RXN',50334,'B',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['SO3       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50335,'B',.FALSE.,0,3,                                    &
['DMS       ','O(3P)     '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50336,'B',.FALSE.,0,4,                                    &
['COS       ','O(3P)     '],                                                   &
['CO        ','SO2       ','          ','          ']),                        &
! retain H2O product here as in the original scheme
asad_flux_defn('RXN',50337,'B',.FALSE.,0,4,                                    &
['SO3       ','H2O       '],                                                   &
['H2SO4     ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50338,'B',.FALSE.,0,5,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','DMSO      ','MeOO      ','          '])                         &
]

! Strat-Trop sulphur chemistry (contains explicit SO3) for
!                                         i_ukca_chem_version >=121
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                         asad_aerosol_chem_strattrop_121(24) = [               &
asad_flux_defn('RXN',50140,'B',.FALSE.,0,3,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50142,'B',.FALSE.,0,3,                                    &
['DMS       ','NO3       '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50143,'B',.FALSE.,0,4,                                    &
['DMSO      ','OH        '],                                                   &
['SO2       ','MSA       ','          ','          ']),                        &
asad_flux_defn('RXN',50144,'B',.FALSE.,0,4,                                    &
['CS2       ','OH        '],                                                   &
['COS       ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50145,'B',.FALSE.,0,4,                                    &
['H2S       ','OH        '],                                                   &
['SO2       ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50146,'B',.FALSE.,0,4,                                    &
['COS       ','OH        '],                                                   &
['CO2       ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50147,'B',.FALSE.,0,3,                                    &
['Monoterp  ','OH        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,3,                                    &
['Monoterp  ','O3        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50149,'B',.FALSE.,0,3,                                    &
['Monoterp  ','NO3       '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50150,'T',.FALSE.,0,4,                                    &
['SO2       ','OH        '],                                                   &
['HO2       ','SO3       ','          ','          ']),                        &
asad_flux_defn('RXN',50151,'H',.FALSE.,0,3,                                    &
['SO2       ','H2O2      '],                                                   &
['NULL0     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50152,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL1     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50153,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL2     ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50154,'D',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50155,'W',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50331,'J',.FALSE.,0,4,                                    &
['COS       ','PHOTON    '],                                                   &
['CO        ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50332,'J',.FALSE.,0,4,                                    &
['H2SO4     ','PHOTON    '],                                                   &
['SO3       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50333,'J',.FALSE.,0,4,                                    &
['SO3       ','PHOTON    '],                                                   &
['SO2       ','O(3P)     ','          ','          ']),                        &
asad_flux_defn('RXN',50334,'B',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['SO3       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50335,'B',.FALSE.,0,3,                                    &
['DMS       ','O(3P)     '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50336,'B',.FALSE.,0,4,                                    &
['COS       ','O(3P)     '],                                                   &
['CO        ','SO2       ','          ','          ']),                        &
! retain H2O product here as in the original scheme
asad_flux_defn('RXN',50337,'B',.FALSE.,0,4,                                    &
['SO3       ','H2O       '],                                                   &
['H2SO4     ','H2O       ','          ','          ']),                        &
! replace 50338 with the following
asad_flux_defn('RXN',50339,'T',.FALSE.,0,4,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','MeOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50340,'T',.FALSE.,0,4,                                    &
['DMS       ','OH        '],                                                   &
['DMSO      ','HO2       ','          ','          '])                         &
]

! Strat-Trop sulphur chemistry (contains explicit SO3) for
!                                         i_ukca_chem_version >=132
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                         asad_aerosol_chem_strattrop_132(27) = [               &
asad_flux_defn('RXN',50140,'B',.FALSE.,0,3,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50142,'B',.FALSE.,0,3,                                    &
['DMS       ','NO3       '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50143,'B',.FALSE.,0,4,                                    &
['DMSO      ','OH        '],                                                   &
['SO2       ','MSA       ','          ','          ']),                        &
asad_flux_defn('RXN',50144,'B',.FALSE.,0,4,                                    &
['CS2       ','OH        '],                                                   &
['COS       ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50145,'B',.FALSE.,0,4,                                    &
['H2S       ','OH        '],                                                   &
['SO2       ','H2O       ','          ','          ']),                        &
asad_flux_defn('RXN',50146,'B',.FALSE.,0,4,                                    &
['COS       ','OH        '],                                                   &
['CO2       ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50147,'B',.FALSE.,0,3,                                    &
['Monoterp  ','OH        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,3,                                    &
['Monoterp  ','O3        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50149,'B',.FALSE.,0,3,                                    &
['Monoterp  ','NO3       '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50150,'T',.FALSE.,0,4,                                    &
['SO2       ','OH        '],                                                   &
['HO2       ','SO3       ','          ','          ']),                        &
asad_flux_defn('RXN',50151,'H',.FALSE.,0,3,                                    &
['SO2       ','H2O2      '],                                                   &
['NULL0     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50152,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL1     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50153,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL2     ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50154,'D',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50155,'W',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50331,'J',.FALSE.,0,4,                                    &
['COS       ','PHOTON    '],                                                   &
['CO        ','SO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50332,'J',.FALSE.,0,4,                                    &
['H2SO4     ','PHOTON    '],                                                   &
['SO3       ','OH        ','          ','          ']),                        &
asad_flux_defn('RXN',50333,'J',.FALSE.,0,4,                                    &
['SO3       ','PHOTON    '],                                                   &
['SO2       ','O(3P)     ','          ','          ']),                        &
asad_flux_defn('RXN',50334,'B',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['SO3       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50335,'B',.FALSE.,0,3,                                    &
['DMS       ','O(3P)     '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50336,'B',.FALSE.,0,4,                                    &
['COS       ','O(3P)     '],                                                   &
['CO        ','SO2       ','          ','          ']),                        &
! retain H2O product here as in the original scheme
asad_flux_defn('RXN',50337,'B',.FALSE.,0,4,                                    &
['SO3       ','H2O       '],                                                   &
['H2SO4     ','H2O       ','          ','          ']),                        &
! replace 50338 with the following
asad_flux_defn('RXN',50339,'T',.FALSE.,0,4,                                    &
['DMS       ','OH        '],                                                   &
['SO2       ','MeOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50340,'T',.FALSE.,0,4,                                    &
['DMS       ','OH        '],                                                   &
['DMSO      ','HO2       ','          ','          ']),                        &
asad_flux_defn('RXN',50042,'B',.TRUE.,0,4,                                     &
['NO3       ','C5H8      '],                                                   &
['ISON      ','SEC_ORG_I ','          ','          ']),                        &
asad_flux_defn('RXN',50440,'B',.TRUE.,0,4,                                     &
['OH        ','C5H8      '],                                                   &
['ISO2      ','SEC_ORG_I ','          ','          ']),                        &
asad_flux_defn('RXN',50441,'B',.TRUE.,0,5,                                     &
['O3        ','C5H8      '],                                                   &
['HO2       ','OH        ','SEC_ORG_I ','          '])                         &
]


! CRI-Strat with updated DMS and monoterp chemistry
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                                     &
                         asad_aerosol_chem_cri(24) = [                         &
asad_flux_defn('RXN',50140,'B',.FALSE.,0,5,                                    &
['DMS       ','OH        '],                                                   &
['MeSCH2OO  ','H2O       ','          ','          ']),                        &
! No equivalent of DMS+OH -> SO2+MSA reactions in CRI (50142)
asad_flux_defn('RXN',50142,'B',.FALSE.,0,6,                                    &
['DMS       ','NO3       '],                                                   &
['MeSCH2OO  ','HONO2     ','          ','          ']),                        &
asad_flux_defn('RXN',50143,'B',.FALSE.,0,4,                                    &
['DMSO      ','OH        '],                                                   &
['MSIA      ','MeOO      ','          ','          ']),                        &
asad_flux_defn('RXN',50144,'B',.FALSE.,0,4,                                    &
['CS2       ','OH        '],                                                   &
['SO2       ','COS       ','          ','          ']),                        &
asad_flux_defn('RXN',50145,'B',.FALSE.,0,3,                                    &
['H2S       ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50146,'B',.FALSE.,0,3,                                    &
['COS       ','OH        '],                                                   &
['SO2       ','          ','          ','          ']),                        &
! CRI equivalent to Monoterp+OH -> Sec_org
asad_flux_defn('RXN',50147,'B',.FALSE.,0,4,                                    &
['APINENE   ','OH        '],                                                   &
['RTN28O2   ','Sec_Org   ','          ','          ']),                        &
asad_flux_defn('RXN',50147,'B',.FALSE.,0,4,                                    &
['BPINENE   ','OH        '],                                                   &
['RTX28O2   ','Sec_Org   ','          ','          ']),                        &
! CRI equivalent to Monoterp+O3 -> sec_org
asad_flux_defn('RXN',50148,'B',.FALSE.,0,3,                                    &
['Monoterp  ','O3        '],                                                   &
['Sec_Org   ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,6,                                    &
['APINENE   ','O3        '],                                                   &
['OH        ','Me2CO     ','RN18AO2   ','Sec_Org   ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,5,                                    &
['APINENE   ','O3        '],                                                   &
['TNCARB26  ','H2O2      ','Sec_Org   ','          ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,4,                                    &
['APINENE   ','O3        '],                                                   &
['RCOOH25   ','Sec_Org   ','          ','          ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,5,                                    &
['BPINENE   ','O3        '],                                                   &
['RTX24O2   ','OH        ','Sec_Org   ','          ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,6,                                    &
['BPINENE   ','O3        '],                                                   &
['HCHO      ','TXCARB24  ','H2O2      ','Sec_Org   ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,5,                                    &
['BPINENE   ','O3        '],                                                   &
['HCHO      ','TXCARB22  ','Sec_Org   ','          ']),                        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,5,                                    &
['BPINENE   ','O3        '],                                                   &
['TXCARB24  ','CO        ','Sec_Org   ','          ']),                        &
! CRI equivalent to Monoterp+NO3 -> sec_org
asad_flux_defn('RXN',50149,'B',.FALSE.,0,4,                                    &
['APINENE   ','NO3       '],                                                   &
['NRTN28O2  ','Sec_Org   ','          ','          ']),                        &
asad_flux_defn('RXN',50149,'B',.FALSE.,0,4,                                    &
['BPINENE   ','NO3       '],                                                   &
['NRTX28O2  ','Sec_Org   ','          ','          ']),                        &
asad_flux_defn('RXN',50150,'T',.FALSE.,0,4,                                    &
['SO2       ','OH        '],                                                   &
['HO2       ','SO3       ','          ','          ']),                        &
asad_flux_defn('RXN',50151,'H',.FALSE.,0,3,                                    &
['SO2       ','H2O2      '],                                                   &
['NULL0     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50152,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL1     ','          ','          ','          ']),                        &
asad_flux_defn('RXN',50153,'H',.FALSE.,0,3,                                    &
['SO2       ','O3        '],                                                   &
['NULL2     ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50154,'D',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          ']),                        &
asad_flux_defn('DEP',50155,'W',.TRUE.,0,1,                                     &
['SO2       ','          '],                                                   &
['          ','          ','          ','          '])                         &
]



PUBLIC :: asad_load_default_fluxes

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ASAD_FLUX_DAT'

CONTAINS

! #####################################################################

SUBROUTINE asad_load_default_fluxes
! To load flux information into asad_chemical_fluxes array for use in
!  asad flux diagnostics scheme. Note that both the total number of
!  fluxes and the number in each section (currently: reaction fluxes;
!  photolytic fluxes; dry deposition; and wet deposition) are limited.
!  Reactions etc are either taken from the definition arrays above
!  or are specified in a user-STASHmaster file, and read in using
!  the chemical definition arrays.
!  NOT currently handling the Strat-trop exchange as this should be
!  in section 38.
!  Called from UKCA_MAIN1.

! Procedure:

! 1) The arrays defining the diagnostics used by asad are selected
!    from alternatives using the model configuration variables

! 2) Generic allocatable arrays are defined and filled with one of the above
!    according to the model configuration, these are used to fill the array
!    asad_chemical_fluxes

! 3) The size of the asad_chemical_fluxes array is calculated from the sizes
!    of the allocated generic arrays and the asad_chemical_fluxes array
!    can then be allocated

! 4) The asad_chemical fluxes array is filled with the contents of the
!    allocated generic arrays

! 5) The generic arrays are deallocated

USE ukca_config_specification_mod,  ONLY: ukca_config
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Normal
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! initialisation string
CHARACTER(LEN = 10), PARAMETER :: initstring = 'XXXXXXXXXX'

INTEGER, PARAMETER :: ichem_version_117 = 117
INTEGER, PARAMETER :: ichem_version_119 = 119
INTEGER, PARAMETER :: ichem_version_121 = 121
INTEGER, PARAMETER :: ichem_version_132 = 132

INTEGER :: i, j            ! counters

! Generic temporary allocatable arrays
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_trop_ox_budget_prod01(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_trop_ox_budget_prod02(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_trop_ox_budget_loss01(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_trop_ox_budget_loss02(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_trop_ox_budget_drydep(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_trop_ox_budget_wetdep(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_oxidN_drydep(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_oxidN_wetdep(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_trop_other_fluxes(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_general_interest(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_trop_co_budget(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_lightning_diags(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_strat_oh_prod(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_strat_oh_loss(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_strat_o3_budget(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_strat_o3_misc(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_rco2no2_pan_prod(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_ro2ho2_reacn(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_ro2no3_reacn(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_ro2ro2_reacn(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_ch4_oxidn(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_o1d_prod(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_aerosol_chem(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_h2o_budget(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_ch4_budget_loss(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_ch4_drydep(:)
TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: aa_ch4_ste(:)

INTEGER :: p1                   ! start position in asad_chemical_fluxes array
INTEGER :: p2                   ! end position in asad_chemical_fluxes array
INTEGER :: errcode              ! error code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_LOAD_DEFAULT_FLUXES'
CHARACTER(LEN=errormessagelength) :: cmessage     ! error message

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ukca_config%l_ukca_strattrop) THEN
  ! Select the asad diagnostics appropriate for StratTrop chemistry and
  ! allocate the generic flux arrays according to the model configuration
  ! ===============================================================

    ! aa_trop_ox_budget_prod01
  ALLOCATE(aa_trop_ox_budget_prod01(SIZE(asad_trop_ox_budget_prod)))
  aa_trop_ox_budget_prod01 = asad_trop_ox_budget_prod
  ! aa_trop_ox_budget_loss01
  IF (ukca_config%i_ukca_chem_version >= ichem_version_132) THEN
    ALLOCATE(aa_trop_ox_budget_loss01(SIZE(asad_trop_ox_budget_loss01_132)))
    aa_trop_ox_budget_loss01 = asad_trop_ox_budget_loss01_132
  ELSE
    ALLOCATE(aa_trop_ox_budget_loss01(SIZE(asad_trop_ox_budget_loss01)))
    aa_trop_ox_budget_loss01 = asad_trop_ox_budget_loss01
  END IF
  ! aa_trop_ox_budget_loss02
  IF (ukca_config%i_ukca_chem_version >= ichem_version_132) THEN
    ALLOCATE(aa_trop_ox_budget_loss02(SIZE(asad_trop_ox_budget_loss02_132)))
    aa_trop_ox_budget_loss02 = asad_trop_ox_budget_loss02_132
  ELSE
    ALLOCATE(aa_trop_ox_budget_loss02(SIZE(asad_trop_ox_budget_loss02)))
    aa_trop_ox_budget_loss02 = asad_trop_ox_budget_loss02
  END IF
  ! aa_trop_ox_budget_drydep
  ALLOCATE(aa_trop_ox_budget_drydep(SIZE(asad_trop_ox_budget_drydep)))
  aa_trop_ox_budget_drydep = asad_trop_ox_budget_drydep
  ! aa_trop_ox_budget_wetdep
  ALLOCATE(aa_trop_ox_budget_wetdep(SIZE(asad_trop_ox_budget_wetdep)))
  aa_trop_ox_budget_wetdep = asad_trop_ox_budget_wetdep
  ! aa_trop_other_fluxes
  IF (ukca_config%i_ukca_chem_version >= ichem_version_132) THEN
    ALLOCATE(aa_trop_other_fluxes(SIZE(asad_trop_other_fluxes_132)))
    aa_trop_other_fluxes = asad_trop_other_fluxes_132
  ELSE
    ALLOCATE(aa_trop_other_fluxes(SIZE(asad_trop_other_fluxes)))
    aa_trop_other_fluxes = asad_trop_other_fluxes
  END IF
  ! aa_general_interest
  ALLOCATE(aa_general_interest(SIZE(asad_general_interest)))
  aa_general_interest = asad_general_interest
  ! aa_trop_co_budget
  IF (ukca_config%l_ukca_ro2_perm) THEN
    IF (ukca_config%i_ukca_chem_version >= ichem_version_121) THEN
      ALLOCATE(aa_trop_co_budget(SIZE(asad_ro2perm_trop_co_budget_121)))
      aa_trop_co_budget = asad_ro2perm_trop_co_budget_121
    ELSE
      ALLOCATE(aa_trop_co_budget(SIZE(asad_ro2perm_trop_co_budget)))
      aa_trop_co_budget = asad_ro2perm_trop_co_budget
    END IF
  ELSE
    IF (ukca_config%i_ukca_chem_version >= ichem_version_121) THEN
      ALLOCATE(aa_trop_co_budget(SIZE(asad_trop_co_budget_121)))
      aa_trop_co_budget = asad_trop_co_budget_121
    ELSE
      ALLOCATE(aa_trop_co_budget(SIZE(asad_trop_co_budget)))
      aa_trop_co_budget = asad_trop_co_budget
    END IF
  END IF
  ! aa_lightning_diags
  ALLOCATE(aa_lightning_diags(SIZE(asad_lightning_diags)))
  aa_lightning_diags = asad_lightning_diags
  ! aa_strat_oh_prod
  IF (ukca_config%i_ukca_chem_version >= ichem_version_121) THEN
    ALLOCATE(aa_strat_oh_prod(SIZE(asad_strat_oh_prod_121)))
    aa_strat_oh_prod = asad_strat_oh_prod_121
  ELSE
    ALLOCATE(aa_strat_oh_prod(SIZE(asad_strat_oh_prod)))
    aa_strat_oh_prod = asad_strat_oh_prod
  END IF
  ! aa_strat_oh_loss
  IF (ukca_config%i_ukca_chem_version >= ichem_version_121) THEN
    ALLOCATE(aa_strat_oh_loss(SIZE(asad_strat_oh_loss_121)))
    aa_strat_oh_loss = asad_strat_oh_loss_121
  ELSE
    ALLOCATE(aa_strat_oh_loss(SIZE(asad_strat_oh_loss)))
    aa_strat_oh_loss = asad_strat_oh_loss
  END IF
  ! aa_strat_o3_budget
  ALLOCATE(aa_strat_o3_budget(SIZE(asad_strat_o3_budget)))
  aa_strat_o3_budget = asad_strat_o3_budget
  ! aa_strat_o3_misc
  ALLOCATE(aa_strat_o3_misc(SIZE(asad_strat_o3_misc)))
  aa_strat_o3_misc = asad_strat_o3_misc
  ! aa_aerosol_chem
  IF (ukca_config%l_ukca_chem_aero) THEN
    IF (ukca_config%i_ukca_chem_version >= ichem_version_132) THEN
      ALLOCATE(aa_aerosol_chem(SIZE(asad_aerosol_chem_strattrop_132)))
      aa_aerosol_chem = asad_aerosol_chem_strattrop_132
    ELSE IF (ukca_config%i_ukca_chem_version >= ichem_version_121) THEN
      ALLOCATE(aa_aerosol_chem(SIZE(asad_aerosol_chem_strattrop_121)))
      aa_aerosol_chem = asad_aerosol_chem_strattrop_121
    ELSE IF ((ukca_config%i_ukca_chem_version < ichem_version_121) .AND.       &
         (ukca_config%i_ukca_chem_version >= ichem_version_117)) THEN
      ALLOCATE(aa_aerosol_chem(SIZE(asad_aerosol_chem_strattrop_117)))
      aa_aerosol_chem = asad_aerosol_chem_strattrop_117
    ELSE
      ALLOCATE(aa_aerosol_chem(SIZE(asad_aerosol_chem_strattrop)))
      aa_aerosol_chem = asad_aerosol_chem_strattrop
    END IF
  END IF
  ! aa_rco2no2_pan_prod
  ALLOCATE(aa_rco2no2_pan_prod(SIZE(asad_rco2no2_pan_prod)))
  aa_rco2no2_pan_prod = asad_rco2no2_pan_prod
  ! aa_ro2ho2_reacn
  ALLOCATE(aa_ro2ho2_reacn(SIZE(asad_ro2ho2_reacn)))
  aa_ro2ho2_reacn = asad_ro2ho2_reacn
  ! aa_ro2no3_reacn
  ALLOCATE(aa_ro2no3_reacn(SIZE(asad_ro2no3_reacn)))
  aa_ro2no3_reacn = asad_ro2no3_reacn
  ! aa_ro2ro2_reacn
  IF (ukca_config%l_ukca_ro2_perm) THEN
    ALLOCATE(aa_ro2ro2_reacn(SIZE(asad_ro2perm_ro2ro2_reacn)))
    aa_ro2ro2_reacn = asad_ro2perm_ro2ro2_reacn
  ELSE
    ALLOCATE(aa_ro2ro2_reacn(SIZE(asad_ro2ro2_reacn)))
    aa_ro2ro2_reacn = asad_ro2ro2_reacn
  END IF
  ! aa_ch4_oxidn
  ALLOCATE(aa_ch4_oxidn(SIZE(asad_ch4_oxidn)))
  aa_ch4_oxidn = asad_ch4_oxidn
  ! aa_o1d_prod
  ALLOCATE(aa_o1d_prod(SIZE(asad_o1d_prod)))
  aa_o1d_prod = asad_o1d_prod
  ! aa_h2o_budget
  IF (ukca_config%i_ukca_chem_version >= ichem_version_121) THEN
    ALLOCATE(aa_h2o_budget(SIZE(asad_h2o_budget_121)))
    aa_h2o_budget = asad_h2o_budget_121
  ELSE
    ALLOCATE(aa_h2o_budget(SIZE(asad_h2o_budget)))
    aa_h2o_budget = asad_h2o_budget
  END IF
  ! aa_oxidN_drydep
  ALLOCATE(aa_oxidN_drydep(SIZE(asad_oxidN_drydep)))
  aa_oxidN_drydep = asad_oxidN_drydep
  ! aa_oxidN_wetdep
  ALLOCATE(aa_oxidN_wetdep(SIZE(asad_oxidN_wetdep)))
  aa_oxidN_wetdep = asad_oxidN_wetdep
  ! aa_ch4_budget_loss
  IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
    ALLOCATE(aa_ch4_budget_loss(SIZE(asad_atmos_ch4_budget_loss)))
    aa_ch4_budget_loss = asad_atmos_ch4_budget_loss
  END IF
  ! aa_ch4_drydep
  IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
    ALLOCATE(aa_ch4_drydep(SIZE(asad_ch4_drydep)))
    aa_ch4_drydep = asad_ch4_drydep
  END IF
  ! aa_ch4_ste
  IF (ukca_config%l_ukca_emsdrvn_ch4) THEN
    ALLOCATE(aa_ch4_ste(SIZE(asad_ch4_ste)))
    aa_ch4_ste = asad_ch4_ste
  END IF

ELSE IF (ukca_config%l_ukca_cristrat) THEN
  ! Select the asad diagnostics appropriate for CRI-Strat chemistry and
  ! allocate the generic flux arrays according to the model configuration
  ! ===============================================================

  ! aa_trop_ox_budget_prod01
  IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
    ALLOCATE(aa_trop_ox_budget_prod01(SIZE(cri_strat2_trop_ox_budget_prod01)))
    aa_trop_ox_budget_prod01 = cri_strat2_trop_ox_budget_prod01
  ELSE
    ALLOCATE(aa_trop_ox_budget_prod01(SIZE(cri_trop_ox_budget_prod01)))
    aa_trop_ox_budget_prod01 = cri_trop_ox_budget_prod01
  END IF
  ! aa_trop_ox_budget_prod02
  IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
    ALLOCATE(aa_trop_ox_budget_prod02(SIZE(cri_strat2_trop_ox_budget_prod02)))
    aa_trop_ox_budget_prod02 = cri_strat2_trop_ox_budget_prod02
  ELSE
    ALLOCATE(aa_trop_ox_budget_prod02(SIZE(cri_trop_ox_budget_prod02)))
    aa_trop_ox_budget_prod02 = cri_trop_ox_budget_prod02
  END IF
  ! aa_trop_ox_budget_loss01
  IF (ukca_config%l_ukca_chem_aero) THEN
    IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
      ALLOCATE(aa_trop_ox_budget_loss01(                                       &
                                  SIZE(cri_strat2_aer_trop_ox_budget_loss01)))
      aa_trop_ox_budget_loss01 = cri_strat2_aer_trop_ox_budget_loss01
    ELSE
      ALLOCATE(aa_trop_ox_budget_loss01(SIZE(cri_aer_trop_ox_budget_loss01)))
      aa_trop_ox_budget_loss01 = cri_aer_trop_ox_budget_loss01
    END IF
  ELSE
    IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
      ALLOCATE(aa_trop_ox_budget_loss01(SIZE(cri_strat2_trop_ox_budget_loss01)))
      aa_trop_ox_budget_loss01 = cri_strat2_trop_ox_budget_loss01
    ELSE
      ALLOCATE(aa_trop_ox_budget_loss01(SIZE(cri_trop_ox_budget_loss01)))
      aa_trop_ox_budget_loss01 = cri_trop_ox_budget_loss01
    END IF
  END IF
  ! aa_trop_ox_budget_loss02
  IF (ukca_config%l_ukca_chem_aero) THEN
    IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
      ALLOCATE(aa_trop_ox_budget_loss02(                                       &
                                  SIZE(cri_strat2_aer_trop_ox_budget_loss02)))
      aa_trop_ox_budget_loss02 = cri_strat2_aer_trop_ox_budget_loss02
    ELSE
      ALLOCATE(aa_trop_ox_budget_loss02(SIZE(cri_aer_trop_ox_budget_loss02)))
      aa_trop_ox_budget_loss02 = cri_aer_trop_ox_budget_loss02
    END IF
  ELSE
    IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
      ALLOCATE(aa_trop_ox_budget_loss02(SIZE(cri_strat2_trop_ox_budget_loss02)))
      aa_trop_ox_budget_loss02 = cri_strat2_trop_ox_budget_loss02
    ELSE
      ALLOCATE(aa_trop_ox_budget_loss02(SIZE(cri_trop_ox_budget_loss02)))
      aa_trop_ox_budget_loss02 = cri_trop_ox_budget_loss02
    END IF
  END IF
  ! aa_trop_ox_budget_drydep
  IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
    ALLOCATE(aa_trop_ox_budget_drydep(SIZE(cri_strat2_trop_ox_budget_drydep)))
    aa_trop_ox_budget_drydep = cri_strat2_trop_ox_budget_drydep
  ELSE
    ALLOCATE(aa_trop_ox_budget_drydep(SIZE(cri_trop_ox_budget_drydep)))
    aa_trop_ox_budget_drydep = cri_trop_ox_budget_drydep
  END IF
  ! aa_trop_ox_budget_wetdep
  IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
    ALLOCATE(aa_trop_ox_budget_wetdep(SIZE(cri_strat2_trop_ox_budget_wetdep)))
    aa_trop_ox_budget_wetdep = cri_strat2_trop_ox_budget_wetdep
  ELSE
    ALLOCATE(aa_trop_ox_budget_wetdep(SIZE(cri_trop_ox_budget_wetdep)))
    aa_trop_ox_budget_wetdep = cri_trop_ox_budget_wetdep
  END IF
  ! aa_trop_other_fluxes
  IF (ukca_config%l_ukca_chem_aero) THEN
    IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
      ALLOCATE(aa_trop_other_fluxes(SIZE(cri_strat2_aer_trop_other_fluxes)))
      aa_trop_other_fluxes = cri_strat2_aer_trop_other_fluxes
    ELSE
      ALLOCATE(aa_trop_other_fluxes(SIZE(cri_aer_trop_other_fluxes)))
      aa_trop_other_fluxes = cri_aer_trop_other_fluxes
    END IF
  ELSE
    IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
      ALLOCATE( aa_trop_other_fluxes(SIZE(cri_strat2_trop_other_fluxes)))
      aa_trop_other_fluxes = cri_strat2_trop_other_fluxes
    ELSE
      ALLOCATE( aa_trop_other_fluxes(SIZE(cri_trop_other_fluxes)))
      aa_trop_other_fluxes = cri_trop_other_fluxes
    END IF
  END IF
  ! aa_general_interest
  ALLOCATE(aa_general_interest(SIZE(asad_general_interest)))
  aa_general_interest = asad_general_interest
  ! aa_trop_co_budget
  IF (ukca_config%l_ukca_chem_aero) THEN
    IF (ukca_config%i_ukca_chem_version >= ichem_version_119 .AND.             &
         ukca_config%i_ukca_chem_version < ichem_version_121) THEN
      ALLOCATE(aa_trop_co_budget(SIZE(cri_strat2_aer_trop_co_budget)))
      aa_trop_co_budget = cri_strat2_aer_trop_co_budget
    ELSE IF (ukca_config%i_ukca_chem_version >= ichem_version_121) THEN
      ALLOCATE(aa_trop_co_budget(SIZE(cri_strat2_aer_trop_co_budget_121)))
      aa_trop_co_budget = cri_strat2_aer_trop_co_budget_121
    ELSE
      ALLOCATE(aa_trop_co_budget(SIZE(cri_aer_trop_co_budget)))
      aa_trop_co_budget = cri_aer_trop_co_budget
    END IF
  ELSE
    IF (ukca_config%i_ukca_chem_version >= ichem_version_119 .AND.             &
         ukca_config%i_ukca_chem_version < ichem_version_121) THEN
      ALLOCATE(aa_trop_co_budget(SIZE(cri_strat2_trop_co_budget)))
      aa_trop_co_budget = cri_strat2_trop_co_budget
    ELSE IF (ukca_config%i_ukca_chem_version >= ichem_version_121) THEN
      ALLOCATE(aa_trop_co_budget(SIZE(cri_strat2_trop_co_budget_121)))
      aa_trop_co_budget = cri_strat2_trop_co_budget_121
    ELSE
      ALLOCATE(aa_trop_co_budget(SIZE(cri_trop_co_budget)))
      aa_trop_co_budget = cri_trop_co_budget
    END IF
  END IF
  ! aa_lightning_diags
  ALLOCATE(aa_lightning_diags(SIZE(asad_lightning_diags)))
  aa_lightning_diags = asad_lightning_diags
  ! aa_strat_oh_prod
  IF (ukca_config%i_ukca_chem_version >= ichem_version_121) THEN
    ALLOCATE(aa_strat_oh_prod(SIZE(asad_strat_oh_prod_121)))
    aa_strat_oh_prod = asad_strat_oh_prod_121
  ELSE
    ALLOCATE(aa_strat_oh_prod(SIZE(asad_strat_oh_prod)))
    aa_strat_oh_prod = asad_strat_oh_prod
  END IF
  ! aa_strat_oh_loss
  IF (ukca_config%i_ukca_chem_version >= ichem_version_121) THEN
    ALLOCATE(aa_strat_oh_loss(SIZE(asad_strat_oh_loss_121)))
    aa_strat_oh_loss = asad_strat_oh_loss_121
  ELSE
    ALLOCATE(aa_strat_oh_loss(SIZE(asad_strat_oh_loss)))
    aa_strat_oh_loss = asad_strat_oh_loss
  END IF
  ! aa_strat_o3_budget
  ALLOCATE(aa_strat_o3_budget(SIZE(asad_strat_o3_budget)))
  aa_strat_o3_budget = asad_strat_o3_budget
  ! aa_strat_o3_misc
  ALLOCATE(aa_strat_o3_misc(SIZE(asad_strat_o3_misc)))
  aa_strat_o3_misc = asad_strat_o3_misc
  ! aa_aerosol_chem
  IF (ukca_config%l_ukca_chem_aero) THEN
    ALLOCATE(aa_aerosol_chem(SIZE(asad_aerosol_chem_cri)))
    aa_aerosol_chem = asad_aerosol_chem_cri
  END IF
  ! aa_rco2no2_pan_prod
  IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
    ALLOCATE(aa_rco2no2_pan_prod(SIZE(cri_strat2_rco2no2_pan_prod)))
    aa_rco2no2_pan_prod = cri_strat2_rco2no2_pan_prod
  ELSE
    ALLOCATE(aa_rco2no2_pan_prod(SIZE(cri_rco2no2_pan_prod)))
    aa_rco2no2_pan_prod = cri_rco2no2_pan_prod
  END IF
  ! aa_ro2ho2_reacn
  IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
    ALLOCATE(aa_ro2ho2_reacn(SIZE(cri_strat2_ro2ho2_reacn)))
    aa_ro2ho2_reacn = cri_strat2_ro2ho2_reacn
  ELSE
    ALLOCATE(aa_ro2ho2_reacn(SIZE(cri_ro2ho2_reacn)))
    aa_ro2ho2_reacn = cri_ro2ho2_reacn
  END IF
  ! aa_ro2no3_reacn
  IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
    ALLOCATE(aa_ro2no3_reacn(SIZE(cri_strat2_ro2no3_reacn)))
    aa_ro2no3_reacn = cri_strat2_ro2no3_reacn
  ELSE
    ALLOCATE(aa_ro2no3_reacn(SIZE(cri_ro2no3_reacn)))
    aa_ro2no3_reacn = cri_ro2no3_reacn
  END IF
  ! aa_ro2ro2_reacn
  IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
    ALLOCATE(aa_ro2ro2_reacn(SIZE(cri_strat2_ro2ro2_reacn)))
    aa_ro2ro2_reacn = cri_strat2_ro2ro2_reacn
  ELSE
    ALLOCATE(aa_ro2ro2_reacn(SIZE(cri_ro2ro2_reacn)))
    aa_ro2ro2_reacn = cri_ro2ro2_reacn
  END IF
  ! aa_ch4_oxidn
  ALLOCATE(aa_ch4_oxidn(SIZE(asad_ch4_oxidn)))
  aa_ch4_oxidn = asad_ch4_oxidn
  ! aa_o1d_prodn
  ALLOCATE(aa_o1d_prod(SIZE(asad_o1d_prod)))
  aa_o1d_prod = asad_o1d_prod
  ! aa_h2o_budget (nothing suitable yet)
  ! aa_oxidN_drydep
  IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
    ALLOCATE(aa_oxidN_drydep(SIZE(cri_strat2_oxidN_drydep)))
    aa_oxidN_drydep = cri_strat2_oxidN_drydep
  ELSE
    ALLOCATE(aa_oxidN_drydep(SIZE(cri_oxidN_drydep)))
    aa_oxidN_drydep = cri_oxidN_drydep
  END IF
  ! aa_oxidN_wetdep
  IF (ukca_config%i_ukca_chem_version >= ichem_version_119) THEN
    ALLOCATE(aa_oxidN_wetdep(SIZE(cri_strat2_oxidN_wetdep)))
    aa_oxidN_wetdep = cri_strat2_oxidN_wetdep
  ELSE
    ALLOCATE(aa_oxidN_wetdep(SIZE(cri_oxidN_wetdep)))
    aa_oxidN_wetdep = cri_oxidN_wetdep
  END IF


ELSE IF (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) THEN
  ! Select the asad diagnostics appropriate for Offline Oxidants chemistry and
  ! allocate the generic flux arrays according to the model configuration
  ! ==========================================================================

    ! aa_aerosol_chem
  ALLOCATE(aa_aerosol_chem(SIZE(asad_aerosol_chem_offline)))
  aa_aerosol_chem = asad_aerosol_chem_offline
ELSE
  cmessage = 'No rules to create asad diagnostics for this configuration'
  errcode = 1
  CALL ereport('asad_flux_dat',errcode,cmessage)
END IF

! Calculate the total number of diagnostics and allocate the
! asad_chemical_fluxes array
n_chemical_fluxes = 0
IF (ALLOCATED(aa_trop_ox_budget_prod01))                                       &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_trop_ox_budget_prod01)
IF (ALLOCATED(aa_trop_ox_budget_prod02))                                       &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_trop_ox_budget_prod02)
IF (ALLOCATED(aa_trop_ox_budget_loss01))                                       &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_trop_ox_budget_loss01)
IF (ALLOCATED(aa_trop_ox_budget_loss02))                                       &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_trop_ox_budget_loss02)
IF (ALLOCATED(aa_trop_ox_budget_drydep))                                       &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_trop_ox_budget_drydep)
IF (ALLOCATED(aa_trop_ox_budget_wetdep))                                       &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_trop_ox_budget_wetdep)
IF (ALLOCATED(aa_oxidN_drydep))                                                &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_oxidN_drydep)
IF (ALLOCATED(aa_oxidN_wetdep))                                                &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_oxidN_wetdep)
IF (ALLOCATED(aa_trop_other_fluxes))                                           &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_trop_other_fluxes)
IF (ALLOCATED(aa_general_interest))                                            &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_general_interest)
IF (ALLOCATED(aa_trop_co_budget))                                              &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_trop_co_budget)
IF (ALLOCATED(aa_lightning_diags))                                             &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_lightning_diags)
IF (ALLOCATED(aa_strat_oh_prod))                                               &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_strat_oh_prod)
IF (ALLOCATED(aa_strat_oh_loss))                                               &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_strat_oh_loss)
IF (ALLOCATED(aa_strat_o3_budget))                                             &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_strat_o3_budget)
IF (ALLOCATED(aa_strat_o3_misc))                                               &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_strat_o3_misc)
IF (ALLOCATED(aa_rco2no2_pan_prod))                                            &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_rco2no2_pan_prod)
IF (ALLOCATED(aa_ro2ho2_reacn))                                                &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_ro2ho2_reacn)
IF (ALLOCATED(aa_ro2no3_reacn))                                                &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_ro2no3_reacn)
IF (ALLOCATED(aa_ro2ro2_reacn))                                                &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_ro2ro2_reacn)
IF (ALLOCATED(aa_ch4_oxidn))                                                   &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_ch4_oxidn)
IF (ALLOCATED(aa_o1d_prod))                                                    &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_o1d_prod)
IF (ALLOCATED(aa_aerosol_chem))                                                &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_aerosol_chem)
IF (ALLOCATED(aa_h2o_budget))                                                  &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_h2o_budget)
IF (ALLOCATED(aa_ch4_budget_loss))                                             &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_ch4_budget_loss)
IF (ALLOCATED(aa_ch4_drydep))                                                  &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_ch4_drydep)
IF (ALLOCATED(aa_ch4_ste))                                                     &
   n_chemical_fluxes = n_chemical_fluxes + SIZE(aa_ch4_ste)

ALLOCATE(asad_chemical_fluxes(n_chemical_fluxes))

! Initialise for checking
asad_chemical_fluxes(:)%reactants(1) = initstring

! Fill the asad_chemical_fluxes array from the allocated generic arrays
p1 = 1
IF (ALLOCATED(aa_trop_ox_budget_prod01)) THEN
  p2 = p1 + SIZE(aa_trop_ox_budget_prod01) - 1
  asad_chemical_fluxes(p1:p2) = aa_trop_ox_budget_prod01
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_trop_ox_budget_prod02)) THEN
  p2 = p1 + SIZE(aa_trop_ox_budget_prod02) - 1
  asad_chemical_fluxes(p1:p2) = aa_trop_ox_budget_prod02(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_trop_ox_budget_loss01)) THEN
  p2 = p1 + SIZE(aa_trop_ox_budget_loss01) - 1
  asad_chemical_fluxes(p1:p2) = aa_trop_ox_budget_loss01(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_trop_ox_budget_loss02)) THEN
  p2 = p1 + SIZE(aa_trop_ox_budget_loss02) - 1
  asad_chemical_fluxes(p1:p2) = aa_trop_ox_budget_loss02(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_trop_ox_budget_drydep)) THEN
  p2 = p1 + SIZE(aa_trop_ox_budget_drydep) - 1
  asad_chemical_fluxes(p1:p2) = aa_trop_ox_budget_drydep(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_trop_ox_budget_wetdep)) THEN
  p2 = p1 + SIZE(aa_trop_ox_budget_wetdep) - 1
  asad_chemical_fluxes(p1:p2) = aa_trop_ox_budget_wetdep(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_oxidN_drydep)) THEN
  p2 = p1 + SIZE(aa_oxidN_drydep) - 1
  asad_chemical_fluxes(p1:p2) = aa_oxidN_drydep(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_oxidN_wetdep)) THEN
  p2 = p1 + SIZE(aa_oxidN_wetdep) - 1
  asad_chemical_fluxes(p1:p2) = aa_oxidN_wetdep(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_trop_other_fluxes)) THEN
  p2 = p1 + SIZE(aa_trop_other_fluxes) - 1
  asad_chemical_fluxes(p1:p2) = aa_trop_other_fluxes(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_general_interest)) THEN
  p2 = p1 + SIZE(aa_general_interest) - 1
  asad_chemical_fluxes(p1:p2) = aa_general_interest(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_trop_co_budget)) THEN
  p2 = p1 + SIZE(aa_trop_co_budget) - 1
  asad_chemical_fluxes(p1:p2) = aa_trop_co_budget(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_lightning_diags)) THEN
  p2 = p1 + SIZE(aa_lightning_diags) - 1
  asad_chemical_fluxes(p1:p2) = aa_lightning_diags(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_strat_oh_prod)) THEN
  p2 = p1 + SIZE(aa_strat_oh_prod) - 1
  asad_chemical_fluxes(p1:p2) = aa_strat_oh_prod(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_strat_oh_loss)) THEN
  p2 = p1 + SIZE(aa_strat_oh_loss) - 1
  asad_chemical_fluxes(p1:p2) = aa_strat_oh_loss(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_strat_o3_budget)) THEN
  p2 = p1 + SIZE(aa_strat_o3_budget) - 1
  asad_chemical_fluxes(p1:p2) = aa_strat_o3_budget(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_strat_o3_misc)) THEN
  p2 = p1 + SIZE(aa_strat_o3_misc) - 1
  asad_chemical_fluxes(p1:p2) = aa_strat_o3_misc(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_rco2no2_pan_prod)) THEN
  p2 = p1 + SIZE(aa_rco2no2_pan_prod) - 1
  asad_chemical_fluxes(p1:p2) = aa_rco2no2_pan_prod(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_ro2ho2_reacn)) THEN
  p2 = p1 + SIZE(aa_ro2ho2_reacn) - 1
  asad_chemical_fluxes(p1:p2) = aa_ro2ho2_reacn(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_ro2no3_reacn)) THEN
  p2 = p1 + SIZE(aa_ro2no3_reacn) - 1
  asad_chemical_fluxes(p1:p2) = aa_ro2no3_reacn(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_ro2ro2_reacn)) THEN
  p2 = p1 + SIZE(aa_ro2ro2_reacn) - 1
  asad_chemical_fluxes(p1:p2) = aa_ro2ro2_reacn(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_ch4_oxidn)) THEN
  p2 = p1 + SIZE(aa_ch4_oxidn) - 1
  asad_chemical_fluxes(p1:p2) = aa_ch4_oxidn(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_o1d_prod)) THEN
  p2 = p1 + SIZE(aa_o1d_prod) - 1
  asad_chemical_fluxes(p1:p2) = aa_o1d_prod(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_aerosol_chem)) THEN
  p2 = p1 + SIZE(aa_aerosol_chem) - 1
  asad_chemical_fluxes(p1:p2) = aa_aerosol_chem(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_h2o_budget)) THEN
  p2 = p1 + SIZE(aa_h2o_budget) - 1
  asad_chemical_fluxes(p1:p2) = aa_h2o_budget(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_ch4_budget_loss)) THEN
  p2 = p1 + SIZE(aa_ch4_budget_loss) - 1
  asad_chemical_fluxes(p1:p2) = aa_ch4_budget_loss(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_ch4_drydep)) THEN
  p2 = p1 + SIZE(aa_ch4_drydep) - 1
  asad_chemical_fluxes(p1:p2) = aa_ch4_drydep(:)
  p1 = p2 + 1
END IF
IF (ALLOCATED(aa_ch4_ste)) THEN
  p2 = p1 + SIZE(aa_ch4_ste) - 1
  asad_chemical_fluxes(p1:p2) = aa_ch4_ste(:)
  p1 = p2 + 1
END IF

IF (p2 /= n_chemical_fluxes) THEN
  cmessage = ' n_chemical_fluxes and p2 are different'
  WRITE(umMessage,'(A)') cmessage
  CALL umPrint(umMessage,src='asad_flux_dat')
  WRITE(umMessage,'(2I6)') n_chemical_fluxes,p2
  CALL umPrint(umMessage,src='asad_flux_dat')
  errcode = 1
  CALL ereport('asad_flux_dat',errcode,cmessage)
END IF

IF (ANY(asad_chemical_fluxes(:)%reactants(1) == initstring)) THEN
  ! Check that entire array is filled
  cmessage = ' Updating asad_chemical_fluxes failed '//                        &
             '- array not completely filled'
  errcode = 1
  CALL ereport('asad_flux_dat',errcode,cmessage)
END IF

IF (printstatus > PrStatus_Normal) THEN
  WRITE(umMessage,'(A,I5)') ' n_chemical_fluxes: ',n_chemical_fluxes
  CALL umPrint(umMessage,src='asad_flux_dat')
  WRITE(umMessage,'(A)') ' '
  CALL umPrint(umMessage,src='asad_flux_dat')
  WRITE(umMessage,'(A)') 'LIST OF DEFINED ASAD DIAGNOSTICS:'
  CALL umPrint(umMessage,src='asad_flux_dat')
  WRITE(umMessage,'(A)') '========================================='
  CALL umPrint(umMessage,src='asad_flux_dat')
  WRITE(umMessage,'(A)') ' '
  CALL umPrint(umMessage,src='asad_flux_dat')
  WRITE(umMessage,'(A)') 'N=STASH=TYPE=R1=======R2========P1========P2'//      &
             '========P3========P4========No.='
  CALL umPrint(umMessage,src='asad_flux_dat')
  j = 0
  DO i = 1, SIZE(asad_chemical_fluxes)
    IF (asad_chemical_fluxes(i)%stash_number /= imdi) THEN
      WRITE(umMessage,'(I3,1x,I5,1x,A3,1x,6A10,I3)') i,                        &
                  asad_chemical_fluxes(i)%stash_number,                        &
                  asad_chemical_fluxes(i)%diag_type,                           &
                  asad_chemical_fluxes(i)%reactants(1),                        &
                  asad_chemical_fluxes(i)%reactants(2),                        &
                  asad_chemical_fluxes(i)%products(1),                         &
                  asad_chemical_fluxes(i)%products(2),                         &
                  asad_chemical_fluxes(i)%products(3),                         &
                  asad_chemical_fluxes(i)%products(4),                         &
                  asad_chemical_fluxes(i)%num_species
      CALL umPrint(umMessage,src='asad_flux_dat')
      j = j + 1
    END IF
  END DO
END IF

! Deallocate the generic arrays
IF (ALLOCATED(aa_ch4_ste))               DEALLOCATE(aa_ch4_ste)
IF (ALLOCATED(aa_ch4_drydep))            DEALLOCATE(aa_ch4_drydep)
IF (ALLOCATED(aa_ch4_budget_loss))       DEALLOCATE(aa_ch4_budget_loss)
IF (ALLOCATED(aa_oxidN_wetdep))          DEALLOCATE(aa_oxidN_wetdep)
IF (ALLOCATED(aa_oxidN_drydep))          DEALLOCATE(aa_oxidN_drydep)
IF (ALLOCATED(aa_h2o_budget))            DEALLOCATE(aa_h2o_budget)
IF (ALLOCATED(aa_o1d_prod))              DEALLOCATE(aa_o1d_prod)
IF (ALLOCATED(aa_ch4_oxidn))             DEALLOCATE(aa_ch4_oxidn)
IF (ALLOCATED(aa_ro2ro2_reacn))          DEALLOCATE(aa_ro2ro2_reacn)
IF (ALLOCATED(aa_ro2no3_reacn))          DEALLOCATE(aa_ro2no3_reacn)
IF (ALLOCATED(aa_ro2ho2_reacn))          DEALLOCATE(aa_ro2ho2_reacn)
IF (ALLOCATED(aa_rco2no2_pan_prod))      DEALLOCATE(aa_rco2no2_pan_prod)
IF (ALLOCATED(aa_aerosol_chem))          DEALLOCATE(aa_aerosol_chem)
IF (ALLOCATED(aa_strat_o3_misc))         DEALLOCATE(aa_strat_o3_misc)
IF (ALLOCATED(aa_strat_o3_budget))       DEALLOCATE(aa_strat_o3_budget)
IF (ALLOCATED(aa_strat_oh_loss))         DEALLOCATE(aa_strat_oh_loss)
IF (ALLOCATED(aa_strat_oh_prod))         DEALLOCATE(aa_strat_oh_prod)
IF (ALLOCATED(aa_lightning_diags))       DEALLOCATE(aa_lightning_diags)
IF (ALLOCATED(aa_trop_co_budget))        DEALLOCATE(aa_trop_co_budget)
IF (ALLOCATED(aa_general_interest))      DEALLOCATE(aa_general_interest)
IF (ALLOCATED(aa_trop_other_fluxes))     DEALLOCATE(aa_trop_other_fluxes)
IF (ALLOCATED(aa_trop_ox_budget_wetdep)) DEALLOCATE(aa_trop_ox_budget_wetdep)
IF (ALLOCATED(aa_trop_ox_budget_drydep)) DEALLOCATE(aa_trop_ox_budget_drydep)
IF (ALLOCATED(aa_trop_ox_budget_loss02)) DEALLOCATE(aa_trop_ox_budget_loss02)
IF (ALLOCATED(aa_trop_ox_budget_loss01)) DEALLOCATE(aa_trop_ox_budget_loss01)
IF (ALLOCATED(aa_trop_ox_budget_prod02)) DEALLOCATE(aa_trop_ox_budget_prod02)
IF (ALLOCATED(aa_trop_ox_budget_prod01)) DEALLOCATE(aa_trop_ox_budget_prod01)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE asad_load_default_fluxes

END MODULE asad_flux_dat
