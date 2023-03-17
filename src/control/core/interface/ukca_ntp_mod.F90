! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module to hold derived type for non-transported prognostic variables
!   and procedures related to their use.
!
!   The module provides the following procedure for the UKCA API.
!
!     ukca_get_ntp_varlist - Return list of names of required NTP fields
!
!   The following additional public procedures are provided for use within UKCA.
!
!     ntp_init        - Initialise NTP structure containing all field
!                       names and an indication of which are required
!     ntp_copy_in     - Copy NTP fields from a 4D array to the NTP structure
!                       (overloaded for different dimension parent arrays)
!     ntp_copy_out    - Copy NTP fields from the NTP structure to a 4D array
!                       (overloaded for different dimension parent arrays)
!     print_all_ntp   - Print details of all fields in the NTP structure
!     ntp_dealloc     - Deallocate data space in NTP structure
!                       allocated by ntp_copy_in
!     name2ntpindex   - General purpose function for finding the index of
!                       an NTP field in the NTP structure given its name
!
!   Non-transported prognostics are variables which are part of the
!   model state but are not transported. For example short lived radicals
!   such as O(3P) where the value from the previous timestep is used
!   as an initial value to improve the solution accuracy.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds,
! University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_ntp_mod

! Standard UM modules used
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,             ONLY: umPrint, umMessage

USE ukca_fieldname_mod, ONLY: maxlen_fieldname, is_mode_ntp,                   &
  fldname_DHPR12O2,                                                            &
  fldname_RU10AO2,                                                             &
  fldname_MACO3,                                                               &
  fldname_RTX22O2,                                                             &
  fldname_RTX24O2,                                                             &
  fldname_NRTX28O2,                                                            &
  fldname_RTX28O2,                                                             &
  fldname_RTN10O2,                                                             &
  fldname_RTN14O2,                                                             &
  fldname_RTN23O2,                                                             &
  fldname_RTN24O2,                                                             &
  fldname_RTN25O2,                                                             &
  fldname_RTN26O2,                                                             &
  fldname_NRTN28O2,                                                            &
  fldname_RTN28O2,                                                             &
  fldname_NRU12O2,                                                             &
  fldname_NRU14O2,                                                             &
  fldname_NRN12O2,                                                             &
  fldname_NRN9O2,                                                              &
  fldname_NRN6O2,                                                              &
  fldname_RU10O2,                                                              &
  fldname_RU12O2,                                                              &
  fldname_RU14O2,                                                              &
  fldname_RN17O2,                                                              &
  fldname_RN14O2,                                                              &
  fldname_RN11O2,                                                              &
  fldname_RN8O2,                                                               &
  fldname_HOCH2CO3,                                                            &
  fldname_RN18AO2,                                                             &
  fldname_RN15AO2,                                                             &
  fldname_RN18O2,                                                              &
  fldname_RN15O2,                                                              &
  fldname_RN12O2,                                                              &
  fldname_RN9O2,                                                               &
  fldname_HOCH2CH2O2,                                                          &
  fldname_RA19CO2,                                                             &
  fldname_RA19AO2,                                                             &
  fldname_RA16O2,                                                              &
  fldname_RA13O2,                                                              &
  fldname_RN16AO2,                                                             &
  fldname_RN13AO2,                                                             &
  fldname_RN19O2,                                                              &
  fldname_RN16O2,                                                              &
  fldname_RN13O2,                                                              &
  fldname_RN10O2,                                                              &
  fldname_n_activ_sum,                                                         &
  fldname_MACRO2,                                                              &
  fldname_ISO2,                                                                &
  fldname_surfarea,                                                            &
  fldname_cdnc,                                                                &
  fldname_HO2S,                                                                &
  fldname_OHS,                                                                 &
  fldname_O1DS,                                                                &
  fldname_O3PS,                                                                &
  fldname_het_ho2,                                                             &
  fldname_het_n2o5,                                                            &
  fldname_TOLP1,                                                               &
  fldname_HOIPO2,                                                              &
  fldname_HOMVKO2,                                                             &
  fldname_MEMALD1,                                                             &
  fldname_OXYL1,                                                               &
  fldname_HOC3H6O2,                                                            &
  fldname_HOC2H4O2,                                                            &
  fldname_MEKO2,                                                               &
  fldname_MeCOCH2OO,                                                           &
  fldname_MeCOC2OO,                                                            &
  fldname_EtCO3,                                                               &
  fldname_i_PrOO,                                                              &
  fldname_s_BuOO,                                                              &
  fldname_n_PrOO,                                                              &
  fldname_MeCO3,                                                               &
  fldname_EtOO,                                                                &
  fldname_MeOO,                                                                &
  fldname_HCl,                                                                 &
  fldname_HO2,                                                                 &
  fldname_BrO,                                                                 &
  fldname_OH,                                                                  &
  fldname_NO2,                                                                 &
  fldname_O1D,                                                                 &
  fldname_O3P,                                                                 &
  fldname_drydiam_nuc_sol,                                                     &
  fldname_drydiam_ait_sol,                                                     &
  fldname_drydiam_acc_sol,                                                     &
  fldname_drydiam_cor_sol,                                                     &
  fldname_drydiam_ait_insol,                                                   &
  fldname_drydiam_acc_insol,                                                   &
  fldname_drydiam_cor_insol,                                                   &
  fldname_wetdiam_ait_sol,                                                     &
  fldname_wetdiam_acc_sol,                                                     &
  fldname_wetdiam_cor_sol,                                                     &
  fldname_aerdens_ait_sol,                                                     &
  fldname_aerdens_acc_sol,                                                     &
  fldname_aerdens_cor_sol,                                                     &
  fldname_aerdens_ait_insol,                                                   &
  fldname_aerdens_acc_insol,                                                   &
  fldname_aerdens_cor_insol,                                                   &
  fldname_pvol_su_ait_sol,                                                     &
  fldname_pvol_bc_ait_sol,                                                     &
  fldname_pvol_oc_ait_sol,                                                     &
  fldname_pvol_so_ait_sol,                                                     &
  fldname_pvol_no3_ait_sol,                                                    &
  fldname_pvol_nh4_ait_sol,                                                    &
  fldname_pvol_h2o_ait_sol,                                                    &
  fldname_pvol_su_acc_sol,                                                     &
  fldname_pvol_bc_acc_sol,                                                     &
  fldname_pvol_oc_acc_sol,                                                     &
  fldname_pvol_ss_acc_sol,                                                     &
  fldname_pvol_no3_acc_sol,                                                    &
  fldname_pvol_nh4_acc_sol,                                                    &
  fldname_pvol_du_acc_sol,                                                     &
  fldname_pvol_so_acc_sol,                                                     &
  fldname_pvol_h2o_acc_sol,                                                    &
  fldname_pvol_su_cor_sol,                                                     &
  fldname_pvol_bc_cor_sol,                                                     &
  fldname_pvol_oc_cor_sol,                                                     &
  fldname_pvol_ss_cor_sol,                                                     &
  fldname_pvol_no3_cor_sol,                                                    &
  fldname_pvol_nh4_cor_sol,                                                    &
  fldname_pvol_du_cor_sol,                                                     &
  fldname_pvol_so_cor_sol,                                                     &
  fldname_pvol_h2o_cor_sol,                                                    &
  fldname_pvol_bc_ait_insol,                                                   &
  fldname_pvol_oc_ait_insol,                                                   &
  fldname_pvol_du_acc_insol,                                                   &
  fldname_pvol_du_cor_insol,                                                   &
  fldname_pvol_nn_acc_sol,                                                     &
  fldname_pvol_nn_cor_sol

IMPLICIT NONE

! All variables and subroutines are private by default
PRIVATE

! Dr hook variables/parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

! subroutines/functions which are public
PUBLIC ntp_init, ukca_get_ntp_varlist, ntp_copy_in, ntp_copy_out,              &
       print_all_ntp, ntp_dealloc, name2ntpindex

! The size of the all_ntp array is defined here.
! If adding or removing entries remember to change
! the size of dim_ntp
INTEGER, PARAMETER, PUBLIC :: dim_ntp = 127

! Type used to hold all information for each non-transported prognostic.
! data_3d, l_required, name
TYPE, PUBLIC :: ntp_type
  REAL,ALLOCATABLE  :: data_3d(:,:,:)   ! only 3D data allowed
  LOGICAL           :: l_required       ! required by UKCA
  CHARACTER(LEN=maxlen_fieldname) :: varname
END TYPE ntp_type

! All non-transported prognostics
TYPE(ntp_type), SAVE, PUBLIC :: all_ntp(dim_ntp)

! Flag to indicate whether the NTP array has been initialised
LOGICAL, SAVE :: l_all_ntp_available = .FALSE.

! List of NTP variables required by UKCA
CHARACTER(LEN=maxlen_fieldname), TARGET, ALLOCATABLE, SAVE ::                  &
  ntp_varnames(:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_NTP_MOD'

! Generic interfaces for NTP copying routines - overloaded according to the
! dimensions of the data

INTERFACE ntp_copy_in
  MODULE PROCEDURE ntp_copy_in_1d
  MODULE PROCEDURE ntp_copy_in_3d
END INTERFACE ntp_copy_in

INTERFACE ntp_copy_out
  MODULE PROCEDURE ntp_copy_out_1d
  MODULE PROCEDURE ntp_copy_out_3d
END INTERFACE ntp_copy_out

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE ntp_init()
! ----------------------------------------------------------------------
! Description:
! initialise values of data in all_ntp structure
!
! Method:
! Variable name is set via a call to add_ntp_item which also determines
! whether a specific entry is required.
!
! To add an additional non-transported prognostic:
! 1. increment dim_ntp above
! 2. add a call to add_ntp_item in this subroutine to add it to the array
! 3. unless it is a chemical species in an explicit BE scheme, add
!    the logic to the function ntp_req to control when it is on or off.
!
! ----------------------------------------------------------------------

IMPLICIT NONE

! Local variables
INTEGER       :: n_reqvar              ! Count of NTP variables required
INTEGER       :: i                     ! Loop counter
INTEGER       :: j                     ! Variable name index
INTEGER       :: errcode               ! error code
CHARACTER(LEN=*), PARAMETER :: uninitialised_name = 'Not set!!!!!!!!!!!!!'
                                       ! Missing data string for NTP name
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NTP_INIT'
CHARACTER(LEN=errormessagelength)   :: cmessage   ! Error return message
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set dummy varname to do an error check at the end
all_ntp(:)%varname    = uninitialised_name
all_ntp(:)%l_required = .FALSE.

! Initialise metadata for all entries. Where these are chemical compounds
! which are not transported, but are stored in the dump, these names must
! be the same as those in the chch_defs arrays as these are used to check
! if a compound is on and where to put the data obtained from D1

! RTX22O2
CALL add_ntp_item(varname=fldname_RTX22O2)

! RTX24O2
CALL add_ntp_item(varname=fldname_RTX24O2)

! NRTX28O2
CALL add_ntp_item(varname=fldname_NRTX28O2)

! RTX28O2
CALL add_ntp_item(varname=fldname_RTX28O2)

! RTN10O2
CALL add_ntp_item(varname=fldname_RTN10O2)

! RTN14O2
CALL add_ntp_item(varname=fldname_RTN14O2)

! RTN23O2
CALL add_ntp_item(varname=fldname_RTN23O2)

! RTN24O2
CALL add_ntp_item(varname=fldname_RTN24O2)

! RTN25O2
CALL add_ntp_item(varname=fldname_RTN25O2)

! RTN26O2
CALL add_ntp_item(varname=fldname_RTN26O2)

! NRTN28O2
CALL add_ntp_item(varname=fldname_NRTN28O2)

! RTN28O2
CALL add_ntp_item(varname=fldname_RTN28O2)

! NRU12O2
CALL add_ntp_item(varname=fldname_NRU12O2)

! NRU14O2
CALL add_ntp_item(varname=fldname_NRU14O2)

! NRN12O2
CALL add_ntp_item(varname=fldname_NRN12O2)

! NRN9O2
CALL add_ntp_item(varname=fldname_NRN9O2)

! NRN6O2
CALL add_ntp_item(varname=fldname_NRN6O2)

! RU10O2
CALL add_ntp_item(varname=fldname_RU10O2)

! RU12O2
CALL add_ntp_item(varname=fldname_RU12O2)

! RU14O2
CALL add_ntp_item(varname=fldname_RU14O2)

! RN17O2
CALL add_ntp_item(varname=fldname_RN17O2)

! RN14O2
CALL add_ntp_item(varname=fldname_RN14O2)

! RN11O2
CALL add_ntp_item(varname=fldname_RN11O2)

! RN8O2
CALL add_ntp_item(varname=fldname_RN8O2)

! HOCH2CO3
CALL add_ntp_item(varname=fldname_HOCH2CO3)

! RN18AO2
CALL add_ntp_item(varname=fldname_RN18AO2)

! RN15AO2
CALL add_ntp_item(varname=fldname_RN15AO2)

! RN18O2
CALL add_ntp_item(varname=fldname_RN18O2)

! RN15O2
CALL add_ntp_item(varname=fldname_RN15O2)

! RN12O2
CALL add_ntp_item(varname=fldname_RN12O2)

! RN9O2
CALL add_ntp_item(varname=fldname_RN9O2)

! HOCH2CH2O2
CALL add_ntp_item(varname=fldname_HOCH2CH2O2)

! RA19CO2
CALL add_ntp_item(varname=fldname_RA19CO2)

! RA19AO2
CALL add_ntp_item(varname=fldname_RA19AO2)

! RA16O2
CALL add_ntp_item(varname=fldname_RA16O2)

! RA13O2
CALL add_ntp_item(varname=fldname_RA13O2)

! RN16AO2
CALL add_ntp_item(varname=fldname_RN16AO2)

! RN13AO2
CALL add_ntp_item(varname=fldname_RN13AO2)

! RN19O2
CALL add_ntp_item(varname=fldname_RN19O2)

! RN16O2
CALL add_ntp_item(varname=fldname_RN16O2)

! RN13O2
CALL add_ntp_item(varname=fldname_RN13O2)

! RN10O2
CALL add_ntp_item(varname=fldname_RN10O2)

! The total number concentration of activated particles
CALL add_ntp_item(varname=fldname_n_activ_sum)

! MACRO2
CALL add_ntp_item(varname=fldname_MACRO2)

! ISO2
CALL add_ntp_item(varname=fldname_ISO2)

! Aerosol surface area
CALL add_ntp_item(varname=fldname_surfarea)

! Cloud droplet number concentration
CALL add_ntp_item(varname=fldname_cdnc)

! Stratospheric HO2
CALL add_ntp_item(varname=fldname_HO2S)

! Stratospheric OH
CALL add_ntp_item(varname=fldname_OHS)

! Stratospheric O(1D)
CALL add_ntp_item(varname=fldname_O1DS)

! Stratospheric O(3P)
CALL add_ntp_item(varname=fldname_O3PS)

! Heterogeneous self reaction rate of HO2
CALL add_ntp_item(varname=fldname_het_ho2)

! Heterogeneous loss rate of N2O5
CALL add_ntp_item(varname=fldname_het_n2o5)

! TOLP1
CALL add_ntp_item(varname=fldname_TOLP1)

! HOIPO2
CALL add_ntp_item(varname=fldname_HOIPO2)

! HOMVKO2
CALL add_ntp_item(varname=fldname_HOMVKO2)

! MEMALD1
CALL add_ntp_item(varname=fldname_MEMALD1)

! OXYL1
CALL add_ntp_item(varname=fldname_OXYL1)

! HOC3H6O2
CALL add_ntp_item(varname=fldname_HOC3H6O2)

! HOC2H4O
CALL add_ntp_item(varname=fldname_HOC2H4O2)

! MEKO2
CALL add_ntp_item(varname=fldname_MEKO2)

! MeCOCH2OO
CALL add_ntp_item(varname=fldname_MeCOCH2OO)

! MeCOC2OO
CALL add_ntp_item(varname=fldname_MeCOC2OO)

! EtCO3
CALL add_ntp_item(varname=fldname_EtCO3)

! i-PrOO
CALL add_ntp_item(varname=fldname_i_PrOO)

! s-BuOO
CALL add_ntp_item(varname=fldname_s_BuOO)

! n-PrOO
CALL add_ntp_item(varname=fldname_n_PrOO)

! MeCO3
CALL add_ntp_item(varname=fldname_MeCO3)

! EtOO
CALL add_ntp_item(varname=fldname_EtOO)

! MeOO
CALL add_ntp_item(varname=fldname_MeOO)

! HCl
CALL add_ntp_item(varname=fldname_HCl)

! HO2
CALL add_ntp_item(varname=fldname_HO2)

! BrO
CALL add_ntp_item(varname=fldname_BrO)

! OH.
CALL add_ntp_item(varname=fldname_OH)

! NO2
CALL add_ntp_item(varname=fldname_NO2)

! O(1D)
CALL add_ntp_item(varname=fldname_O1D)

! O(3P)
CALL add_ntp_item(varname=fldname_O3P)

! DHPR12O2
CALL add_ntp_item(varname=fldname_DHPR12O2)
! RU10AO2
CALL add_ntp_item(varname=fldname_RU10AO2)
! MACO3
CALL add_ntp_item(varname=fldname_MACO3)

! items for use by RADAER and hybrid resolution model:
! Dry diameter
CALL add_ntp_item(varname=fldname_drydiam_nuc_sol)
CALL add_ntp_item(varname=fldname_drydiam_ait_sol)
CALL add_ntp_item(varname=fldname_drydiam_acc_sol)
CALL add_ntp_item(varname=fldname_drydiam_cor_sol)
CALL add_ntp_item(varname=fldname_drydiam_ait_insol)
CALL add_ntp_item(varname=fldname_drydiam_acc_insol)
CALL add_ntp_item(varname=fldname_drydiam_cor_insol)

! Wet diameter
CALL add_ntp_item(varname=fldname_wetdiam_ait_sol)
CALL add_ntp_item(varname=fldname_wetdiam_acc_sol)
CALL add_ntp_item(varname=fldname_wetdiam_cor_sol)

! Aerosol density
CALL add_ntp_item(varname=fldname_aerdens_ait_sol)
CALL add_ntp_item(varname=fldname_aerdens_acc_sol)
CALL add_ntp_item(varname=fldname_aerdens_cor_sol)
CALL add_ntp_item(varname=fldname_aerdens_ait_insol)
CALL add_ntp_item(varname=fldname_aerdens_acc_insol)
CALL add_ntp_item(varname=fldname_aerdens_cor_insol)

! Partial volume
CALL add_ntp_item(varname=fldname_pvol_su_ait_sol)
CALL add_ntp_item(varname=fldname_pvol_bc_ait_sol)
CALL add_ntp_item(varname=fldname_pvol_oc_ait_sol)
CALL add_ntp_item(varname=fldname_pvol_so_ait_sol)
CALL add_ntp_item(varname=fldname_pvol_no3_ait_sol)
CALL add_ntp_item(varname=fldname_pvol_nh4_ait_sol)
CALL add_ntp_item(varname=fldname_pvol_h2o_ait_sol)
CALL add_ntp_item(varname=fldname_pvol_su_acc_sol)
CALL add_ntp_item(varname=fldname_pvol_bc_acc_sol)
CALL add_ntp_item(varname=fldname_pvol_oc_acc_sol)
CALL add_ntp_item(varname=fldname_pvol_ss_acc_sol)
CALL add_ntp_item(varname=fldname_pvol_no3_acc_sol)
CALL add_ntp_item(varname=fldname_pvol_nh4_acc_sol)
CALL add_ntp_item(varname=fldname_pvol_nn_acc_sol)
CALL add_ntp_item(varname=fldname_pvol_du_acc_sol)
CALL add_ntp_item(varname=fldname_pvol_so_acc_sol)
CALL add_ntp_item(varname=fldname_pvol_h2o_acc_sol)
CALL add_ntp_item(varname=fldname_pvol_su_cor_sol)
CALL add_ntp_item(varname=fldname_pvol_bc_cor_sol)
CALL add_ntp_item(varname=fldname_pvol_oc_cor_sol)
CALL add_ntp_item(varname=fldname_pvol_ss_cor_sol)
CALL add_ntp_item(varname=fldname_pvol_no3_cor_sol)
CALL add_ntp_item(varname=fldname_pvol_nh4_cor_sol)
CALL add_ntp_item(varname=fldname_pvol_nn_cor_sol)
CALL add_ntp_item(varname=fldname_pvol_du_cor_sol)
CALL add_ntp_item(varname=fldname_pvol_so_cor_sol)
CALL add_ntp_item(varname=fldname_pvol_h2o_cor_sol)
CALL add_ntp_item(varname=fldname_pvol_bc_ait_insol)
CALL add_ntp_item(varname=fldname_pvol_oc_ait_insol)
CALL add_ntp_item(varname=fldname_pvol_du_acc_insol)
CALL add_ntp_item(varname=fldname_pvol_du_cor_insol)

! Finally, check metadata for all entries is set
IF (ANY(all_ntp(:)%varname == uninitialised_name)) THEN
  ! If some values are not set, write some useful messages and stop the model
  WRITE(umMessage,'(A)') 'one or more entries in all_ntp not set'
  CALL umPrint(umMessage,src=RoutineName)
  CALL print_all_ntp()
  WRITE(cmessage,'(A)') 'one or more entries in all_ntp not set'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Create list of required variables
n_reqvar = 0
DO i = 1, dim_ntp
  IF (all_ntp(i)%l_required) n_reqvar = n_reqvar + 1
END DO
ALLOCATE(ntp_varnames(n_reqvar))
j = 0
DO i = 1, dim_ntp
  IF (all_ntp(i)%l_required) THEN
    j = j + 1
    ntp_varnames(j) = all_ntp(i)%varname
  END IF
END DO

! Set flag to show availability of all_ntp array
l_all_ntp_available = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ntp_init

! ----------------------------------------------------------------------
SUBROUTINE add_ntp_item(varname)
! ----------------------------------------------------------------------
! Description:
! Add data to the all_ntp array. If try and add too many stop the model.
! 1. Increment the internal counter ntp_index
! 2. Check if the counter is within the the array size
! 3. Add the values passed in
! 4. Set the logical flag for whether this is required by calling
!    the function ntp_req
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine argument
CHARACTER(LEN=*), INTENT(IN) :: varname

! Local variables
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ADD_NTP_ITEM'
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message

! Counter for position in all_ntp array. Initialised to zero
! and incremented by one on each call to this routine.
INTEGER, SAVE :: ntp_index = 0
INTEGER       :: errcode                                   ! error code
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! increment the location we store information to
ntp_index = ntp_index + 1

! Test that the all_ntp array is big enough - if not, abort.
IF (ntp_index > dim_ntp) THEN
  WRITE(cmessage,'(A,I6)') 'dim_ntp too small in ukca_ntp_mod:', dim_ntp
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Set name for this entry and use the logical function ntp_req to test
! whether this is required for the current model run
all_ntp(ntp_index)%varname = varname
all_ntp(ntp_index)%l_required = ntp_req(varname)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE add_ntp_item

! ----------------------------------------------------------------------
LOGICAL FUNCTION ntp_req(varname)
! ----------------------------------------------------------------------
! Description:
! Use run time logic to test if varname is on in any specific model run.
! For chemistry, test nadvt. Need to treat others as special cases.
!
! ----------------------------------------------------------------------
USE asad_mod, ONLY: nadvt, O1D_in_ss, O3P_in_ss
USE ukca_config_specification_mod, ONLY: ukca_config, glomap_config,           &
                                         i_ukca_activation_arg,                &
                                         i_ukca_activation_jones,              &
                                         int_method_be_explicit,               &
                                         glomap_variables

USE ukca_mode_setup,       ONLY: mode_nuc_sol, mode_ait_sol, mode_acc_sol,     &
                                 mode_cor_sol, mode_ait_insol, mode_acc_insol, &
                                 mode_cor_insol, cp_su, cp_bc, cp_oc, cp_cl,   &
                                 cp_no3, cp_du, cp_so, cp_nh4, cp_nn

IMPLICIT NONE

! Function argument
CHARACTER(LEN=*), INTENT(IN) :: varname

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
LOGICAL, POINTER :: mode (:)

! Local variables
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NTP_REQ'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mode        => glomap_variables%mode

ntp_req = .FALSE.

! First handle the GLOMAP variables. To allow for the use of allocatable
! arrays for mode and component, handle these together

IF (is_mode_ntp(varname)) THEN

  IF (ukca_config%l_ukca_mode) THEN
    SELECT CASE (varname)

      ! Aerosol surface area stored in dump if using GLOMAP-mode.
      ! Only needed as prognostic under certain circumstances,
      ! however, we want to make it available as a diagnostic whenever
      ! GLOMAP-mode is being used.
    CASE (fldname_surfarea)
      ntp_req = .TRUE.

      ! RADAER items, check which modes and components are active
      ! Dry diameter
    CASE (fldname_drydiam_ait_sol)          ! Aitken-sol dry diameter
      ntp_req = mode(mode_ait_sol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_drydiam_acc_sol)          ! accumulation-sol dry diameter
      ntp_req = mode(mode_acc_sol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_drydiam_cor_sol)          ! coarse-sol dry diameter
      ntp_req = mode(mode_cor_sol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_drydiam_ait_insol)        ! Aitken-sol dry diameter
      ntp_req = mode(mode_ait_insol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_drydiam_acc_insol)        ! accumulation-sol dry diameter
      ntp_req = mode(mode_acc_insol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_drydiam_cor_insol)        ! coarse-sol dry diameter
      ntp_req = mode(mode_cor_insol) .AND. glomap_config%l_ukca_radaer

      ! Wet diameter
    CASE (fldname_wetdiam_ait_sol)          ! Aitken-sol wet diameter
      ntp_req = mode(mode_ait_sol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_wetdiam_acc_sol)          ! accumulation-sol wet diameter
      ntp_req = mode(mode_acc_sol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_wetdiam_cor_sol)          ! coarse-sol wet diameter
      ntp_req = mode(mode_cor_sol) .AND. glomap_config%l_ukca_radaer

      ! Aerosol density
    CASE (fldname_aerdens_ait_sol)          ! Aitken-sol aerosol density
      ntp_req = mode(mode_ait_sol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_aerdens_acc_sol)          ! accumulation-sol " density
      ntp_req = mode(mode_acc_sol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_aerdens_cor_sol)          ! coarse-sol aerosol density
      ntp_req = mode(mode_cor_sol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_aerdens_ait_insol)        ! Aitken-sol aerosol density
      ntp_req = mode(mode_ait_insol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_aerdens_acc_insol)        ! accumulation-sol " density
      ntp_req = mode(mode_acc_insol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_aerdens_cor_insol)        ! coarse-sol aerosol density
      ntp_req = mode(mode_cor_insol) .AND. glomap_config%l_ukca_radaer

      ! Partial volume
    CASE (fldname_pvol_su_ait_sol)          ! Aitken-sol sulphate
      ntp_req = component(mode_ait_sol,cp_su) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_bc_ait_sol)          ! Aitken-sol black carbon
      ntp_req = component(mode_ait_sol,cp_bc) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_oc_ait_sol)          ! Aitken-sol organic matter
      ntp_req = component(mode_ait_sol,cp_oc) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_so_ait_sol)          ! Aitken-sol secondary organic
      ntp_req = component(mode_ait_sol,cp_so) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_no3_ait_sol)         ! Aitken-sol nitrate
      IF (UBOUND(component,DIM=2) >= cp_no3)                                   &
      ntp_req = component(mode_ait_sol,cp_no3) .AND.                           &
                glomap_config%l_ukca_radaer
    CASE (fldname_pvol_nh4_ait_sol)         ! Aitken-sol ammonium
      IF (UBOUND(component,DIM=2) >= cp_nh4)                                   &
      ntp_req = component(mode_ait_sol,cp_nh4) .AND.                           &
                glomap_config%l_ukca_radaer
    CASE (fldname_pvol_h2o_ait_sol)         ! Aitken-sol H2O
      ntp_req = mode(mode_ait_sol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_su_acc_sol)          ! accumulation-sol sulphate
      ntp_req = component(mode_acc_sol,cp_su) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_bc_acc_sol)          ! accumulation-sol black carbon
      ntp_req = component(mode_acc_sol,cp_bc) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_oc_acc_sol)          ! accumulatn-sol organic matter
      ntp_req = component(mode_acc_sol,cp_oc) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_ss_acc_sol)          ! accumulation-sol sea-salt
      ntp_req = component(mode_acc_sol,cp_cl) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_no3_acc_sol)         ! accumulation-sol nitrate
      IF (UBOUND(component,DIM=2) >= cp_no3)                                   &
      ntp_req = component(mode_acc_sol,cp_no3) .AND.                           &
                glomap_config%l_ukca_radaer
    CASE (fldname_pvol_nh4_acc_sol)         ! accumulation-sol ammonium
      IF (UBOUND(component,DIM=2) >= cp_nh4)                                   &
      ntp_req = component(mode_acc_sol,cp_nh4) .AND.                           &
                glomap_config%l_ukca_radaer
    CASE (fldname_pvol_nn_acc_sol)         ! accum-sol sodium nitrate
      IF (UBOUND(component,DIM=2) >= cp_nn)                                    &
      ntp_req = component(mode_acc_sol,cp_nn) .AND.                            &
                glomap_config%l_ukca_radaer
    CASE (fldname_pvol_du_acc_sol)          ! accumulation-sol dust
      ntp_req = component(mode_acc_sol,cp_du) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_so_acc_sol)          ! accumulation-sol 2ndy organic
      ntp_req = component(mode_acc_sol,cp_so) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_h2o_acc_sol)         ! accumulation-sol H2O
      ntp_req = mode(mode_acc_sol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_su_cor_sol)          ! coarse-sol sulphate
      ntp_req = component(mode_cor_sol,cp_su) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_bc_cor_sol)          ! coarse-sol black carbon
      ntp_req = component(mode_cor_sol,cp_bc) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_oc_cor_sol)          ! coarse-sol organic matter
      ntp_req = component(mode_cor_sol,cp_oc) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_ss_cor_sol)          ! coarse-sol sea-salt
      ntp_req = component(mode_cor_sol,cp_cl) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_no3_cor_sol)         ! coarse-sol nitrate
      IF (UBOUND(component,DIM=2) >= cp_no3)                                   &
      ntp_req = component(mode_cor_sol,cp_no3) .AND.                           &
                glomap_config%l_ukca_radaer
    CASE (fldname_pvol_nh4_cor_sol)         ! coarse-sol ammonium
      IF (UBOUND(component,DIM=2) >= cp_nh4)                                   &
      ntp_req = component(mode_cor_sol,cp_nh4) .AND.                           &
                glomap_config%l_ukca_radaer
    CASE (fldname_pvol_nn_cor_sol)         ! coarse-sol sodium nitrate
      IF (UBOUND(component,DIM=2) >= cp_nn)                                    &
      ntp_req = component(mode_cor_sol,cp_nn) .AND.                            &
                glomap_config%l_ukca_radaer
    CASE (fldname_pvol_du_cor_sol)          ! coarse-sol dust
      ntp_req = component(mode_cor_sol,cp_du) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_so_cor_sol)          ! coarse-sol 2ndy organic
      ntp_req = component(mode_cor_sol,cp_so) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_h2o_cor_sol)         ! coarse-sol H2O
      ntp_req = mode(mode_cor_sol) .AND. glomap_config%l_ukca_radaer
    CASE (fldname_pvol_bc_ait_insol)        ! coarse-insol black carbon
      ntp_req = component(mode_ait_insol,cp_bc) .AND.                          &
                glomap_config%l_ukca_radaer
    CASE (fldname_pvol_oc_ait_insol)        ! coarse-insol organic matter
      ntp_req = component(mode_ait_insol,cp_oc) .AND.                          &
                glomap_config%l_ukca_radaer
    CASE (fldname_pvol_du_acc_insol)        ! accumulation-insol dust
      ntp_req = component(mode_acc_insol,cp_du) .AND.                          &
                glomap_config%l_ukca_radaer
    CASE (fldname_pvol_du_cor_insol)        ! coarse-insol dust
      ntp_req = component(mode_cor_insol,cp_du) .AND.                          &
                glomap_config%l_ukca_radaer

      ! In hybrid model it's possible to run ACTIVATE only in senior
      ! component, which requires the dry diameter for the nucleation
      ! soluble mode, which is the one mode not already needed for RADAER.
    CASE (fldname_drydiam_nuc_sol)
      ntp_req = mode(mode_nuc_sol) .AND. glomap_config%l_ntpreq_dryd_nuc_sol

    CASE DEFAULT
      ntp_req = .FALSE.
    END SELECT
  ELSE
    ntp_req = .FALSE.
  END IF

ELSE
  SELECT CASE (varname)

    ! CDNC prognostics - If CDNC is calculated
  CASE (fldname_cdnc)
    ntp_req =                                                                  &
      (glomap_config%i_ukca_activation_scheme == i_ukca_activation_jones) .OR. &
      (glomap_config%i_ukca_activation_scheme == i_ukca_activation_arg)

    ! The total number concentration of activated aerosol
  CASE (fldname_n_activ_sum)
    ntp_req =                                                                  &
      ((glomap_config%i_ukca_activation_scheme == i_ukca_activation_jones) .OR.&
       (glomap_config%i_ukca_activation_scheme == i_ukca_activation_arg)) .AND.&
      glomap_config%l_ntpreq_n_activ_sum

    ! The heterogeneous loss rates are on if l_ukca_trophet
    ! is true
  CASE (fldname_het_ho2, fldname_het_n2o5)
    ntp_req = ukca_config%l_ukca_trophet

    ! O1D and O3P - special case, used by both NR and BE schemes
    ! and on if they are in SS (some NR schemes have O3P as a tracer
    ! in which case it is stored in the tracer array not in the NTP
    ! structure)
  CASE (fldname_O1D)
    ntp_req = O1D_in_ss

  CASE (fldname_O3P)
    ntp_req = O3P_in_ss

    ! Lumped species - on for chemical schemes using lumping
  CASE (fldname_NO2, fldname_BrO, fldname_HCl)
    ntp_req = ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_strattrop .OR.  &
              ukca_config%l_ukca_stratcfc .OR. ukca_config%l_ukca_cristrat

    ! All others are checked against whether they are in
    ! the nadvt array (only for BE explicit schemes)
  CASE DEFAULT
    ! Logic expanded to include case with RO2 species not transported
    ! if StratTrop chemical mechanism is being used, or CRI-Strat used
    ntp_req = ANY(nadvt(:) == varname) .AND.                                   &
         ((ukca_config%ukca_int_method  == int_method_be_explicit) .OR.        &
          (ukca_config%l_ukca_ro2_ntp .AND. ukca_config%l_ukca_strattrop) .OR. &
          ukca_config%l_ukca_cristrat )

  END SELECT
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END FUNCTION ntp_req

! ----------------------------------------------------------------------
SUBROUTINE ukca_get_ntp_varlist(ntp_varnames_ptr, error_code,                  &
                                error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
! UKCA API routine to return pointer to the array of variable names of
! required NTP variables. This may have zero length if no NTP variable
! is marked as required. A non-zero error code is returned if the all_ntp
! array has not been initialised.
! ----------------------------------------------------------------------

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname, errcode_ntp_uninit

IMPLICIT NONE

!Subroutine arguments
CHARACTER(LEN=maxlen_fieldname), POINTER, INTENT(OUT) :: ntp_varnames_ptr(:)
INTEGER, INTENT(OUT) :: error_code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine

! Local variables
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_GET_NTP_VARLIST'

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

IF (.NOT. l_all_ntp_available) THEN
  error_code = errcode_ntp_uninit
  IF (PRESENT(error_message)) error_message =                                  &
    'Non-transported prognostics requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  NULLIFY(ntp_varnames_ptr)
  RETURN
END IF

ntp_varnames_ptr => ntp_varnames

END SUBROUTINE ukca_get_ntp_varlist

! ----------------------------------------------------------------------
SUBROUTINE ntp_copy_in_1d(ntp_data, model_levels,                              &
                          error_code, error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
! Variant of the generic procedure ntp_copy_in.
! Copy NTP fields from the given ntp_data array into the all_ntp data
! structure. Each input field is defined for a single column.
! Bounds for the vertical array dimension are checked for validity based
! on the UKCA extent i.e. 1:model_levels.
! The input fields must span the UKCA extent but may optionally extend
! beyond it in either direction. The vertical extent of the data copied
! is restricted to the UKCA bounds.
! ----------------------------------------------------------------------

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_ntp_uninit, errcode_ntp_mismatch

IMPLICIT NONE

! Subroutine arguments
REAL, ALLOCATABLE, INTENT(IN) :: ntp_data(:,:)  ! NTP data field array
INTEGER, INTENT(IN) :: model_levels        ! Size of UKCA z dimension
INTEGER, INTENT(OUT) :: error_code         ! Return status code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
                                           ! Return error message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
                                           ! Routine name where error trapped

! Local variables

INTEGER :: n_reqvar ! Count of NTP variables required
INTEGER :: i_ntp    ! Index of NTP variable in data field array
INTEGER :: i        ! Loop counter

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NTP_COPY_IN_1D'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check that the all_ntp structure has been initialised
IF (.NOT. l_all_ntp_available) THEN
  error_code = errcode_ntp_uninit
  IF (PRESENT(error_message))                                                  &
    error_message =                                                            &
      'Non-transported prognostics requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check the data field array contains the expected number of NTP fields
n_reqvar = 0
DO i = 1, dim_ntp
  IF (all_ntp(i)%l_required) n_reqvar = n_reqvar + 1
END DO
IF (SIZE(ntp_data,DIM=2) /= n_reqvar) THEN
  error_code = errcode_ntp_mismatch
  IF (PRESENT(error_message))                                                  &
    error_message =                                                            &
      'Number of non-transported prognostics fields does not match requirement'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (n_reqvar > 0) THEN

  ! Check array bounds for data fields are compatible with UKCA.
  ! The field data supplied must fill the required domain but may extend
  ! beyond it to avoid the need for pre-trimming (e.g. halo removal) by
  ! the parent model.
  IF (LBOUND(ntp_data,DIM=1) > 1 .OR. UBOUND(ntp_data,DIM=1) < model_levels)   &
  THEN
    error_code = errcode_ntp_mismatch
    IF (PRESENT(error_message))                                                &
      error_message =                                                          &
        'The non-transported prognostics fields have one or more invalid ' //  &
        'array bounds'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  ! Copy data to all_ntp structure.
  ! Any data outside the required bounds (e.g. halos) are discarded.
  i_ntp = 0
  DO i = 1, dim_ntp
    IF (all_ntp(i)%l_required) THEN
      i_ntp = i_ntp + 1
      ALLOCATE(all_ntp(i)%data_3d(1,1,model_levels))
      all_ntp(i)%data_3d(1,1,:) = ntp_data(1:model_levels,i_ntp)
    END IF
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ntp_copy_in_1d

! ----------------------------------------------------------------------
SUBROUTINE ntp_copy_in_3d(ntp_data, row_length, rows, model_levels,            &
                          error_code, error_message, error_routine)
! ----------------------------------------------------------------------
! Description:
! Variant of the generic procedure ntp_copy_in.
! Copy NTP fields from the given ntp_data array into the all_ntp data
! structure. Each field is defined on a 3D grid.
! Bounds are checked for validity based on the UKCA horizontal and
! vertical extents i.e. 1:row_length, 1:rows and 1:model_levels.
! The input fields must span the UKCA extents but may optionally extend
! beyond them. The horizontal and vertical extents of the data copied
! are restricted to the UKCA bounds.
! ----------------------------------------------------------------------

USE ukca_error_mod, ONLY: maxlen_message, maxlen_procname,                     &
                          errcode_ntp_uninit, errcode_ntp_mismatch

IMPLICIT NONE

! Subroutine arguments
REAL, ALLOCATABLE, INTENT(IN) :: ntp_data(:,:,:,:)  ! NTP data field array
INTEGER, INTENT(IN) :: row_length          ! Size of UKCA x dimension
INTEGER, INTENT(IN) :: rows                ! Size of UKCA y dimension
INTEGER, INTENT(IN) :: model_levels        ! Size of UKCA z dimension
INTEGER, INTENT(OUT) :: error_code         ! Return status code
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
                                           ! Return error message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
                                           ! Routine name where error trapped

! Local variables

INTEGER :: n_reqvar ! Count of NTP variables required
INTEGER :: i_ntp    ! Index of NTP variable in data field array
INTEGER :: i        ! Loop counter

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NTP_COPY_IN_3D'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

! Check that the all_ntp structure has been initialised
IF (.NOT. l_all_ntp_available) THEN
  error_code = errcode_ntp_uninit
  IF (PRESENT(error_message))                                                  &
    error_message =                                                            &
      'Non-transported prognostics requirement has not been initialised'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check the data field array contains the expected number of NTP fields
n_reqvar = 0
DO i = 1, dim_ntp
  IF (all_ntp(i)%l_required) n_reqvar = n_reqvar + 1
END DO
IF (SIZE(ntp_data,DIM=4) /= n_reqvar) THEN
  error_code = errcode_ntp_mismatch
  IF (PRESENT(error_message))                                                  &
    error_message =                                                            &
      'Number of non-transported prognostics fields does not match requirement'
  IF (PRESENT(error_routine)) error_routine = RoutineName
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (n_reqvar > 0) THEN

  ! Check array bounds for data fields are compatible with UKCA.
  ! The field data supplied must fill the required domain but may extend
  ! beyond it to avoid the need for pre-trimming (e.g. halo removal) by
  ! the parent model.
  IF (LBOUND(ntp_data,DIM=1) > 1 .OR. UBOUND(ntp_data,DIM=1) < row_length .OR. &
      LBOUND(ntp_data,DIM=2) > 1 .OR. UBOUND(ntp_data,DIM=2) < rows .OR.       &
      LBOUND(ntp_data,DIM=3) > 1 .OR. UBOUND(ntp_data,DIM=3) < model_levels)   &
  THEN
    error_code = errcode_ntp_mismatch
    IF (PRESENT(error_message))                                                &
      error_message =                                                          &
        'The non-transported prognostics fields have one or more invalid ' //  &
        'array bounds'
    IF (PRESENT(error_routine)) error_routine = RoutineName
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  ! Copy data to all_ntp structure.
  ! Any data outside the required bounds (e.g. halos) are discarded.
  i_ntp = 0
  DO i = 1, dim_ntp
    IF (all_ntp(i)%l_required) THEN
      i_ntp = i_ntp + 1
      ALLOCATE(all_ntp(i)%data_3d(row_length,rows,model_levels))
      all_ntp(i)%data_3d(:,:,:) = ntp_data(1:row_length,1:rows,1:model_levels, &
                                           i_ntp)
    END IF
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ntp_copy_in_3d

! ----------------------------------------------------------------------
SUBROUTINE ntp_copy_out_1d(model_levels, ntp_data)
! ----------------------------------------------------------------------
! Description:
! Variant of the generic procedure ntp_copy_out.
! Copy the NTP fields from the all_ntp data structure to the given
! ntp_data array. Each NTP field is defined on a single column.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: model_levels          ! Size of UKCA z dimension
REAL, ALLOCATABLE, INTENT(IN OUT) :: ntp_data(:,:)    ! NTP data fields

! Local variables

INTEGER :: i_ntp   ! Index of NTP variable in data field array
INTEGER :: i       ! Loop counter

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NTP_COPY_OUT_1D'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Copy data from all_ntp structure
i_ntp = 0
DO i = 1, dim_ntp
  IF (all_ntp(i)%l_required) THEN
    i_ntp = i_ntp + 1
    ntp_data(1:model_levels,i_ntp) = all_ntp(i)%data_3d(1,1,:)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ntp_copy_out_1d

! ----------------------------------------------------------------------
SUBROUTINE ntp_copy_out_3d(row_length, rows, model_levels, ntp_data)
! ----------------------------------------------------------------------
! Description:
! Variant of the generic procedure ntp_copy_out.
! Copy the NTP fields from the all_ntp data structure to the given
! ntp_data array. Each NTP field is defined on a 3D grid.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: row_length          ! Size of UKCA x dimension
INTEGER, INTENT(IN) :: rows                ! Size of UKCA y dimension
INTEGER, INTENT(IN) :: model_levels        ! Size of UKCA z dimension
REAL, ALLOCATABLE, INTENT(IN OUT) :: ntp_data(:,:,:,:)    ! NTP data fields

! Local variables

INTEGER :: i_ntp   ! Index of NTP variable in data field array
INTEGER :: i       ! Loop counter

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NTP_COPY_OUT_3D'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Copy data from all_ntp structure
i_ntp = 0
DO i = 1, dim_ntp
  IF (all_ntp(i)%l_required) THEN
    i_ntp = i_ntp + 1
    ntp_data(1:row_length,1:rows,1:model_levels,i_ntp) =                       &
      all_ntp(i)%data_3d(:,:,:)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ntp_copy_out_3d

! ----------------------------------------------------------------------
SUBROUTINE print_all_ntp()
! ----------------------------------------------------------------------
! Description:
! Print summary of all information in the all_ntp array.
!
! Method:
! Loop through each entry and print values
! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'PRINT_ALL_NTP'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Print header
WRITE(umMessage,'(A)') ' ==========================='
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)') ' Non-transported Prognostics'
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)') ' ==========================='
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)') '  i l_required varname'
CALL umPrint(umMessage,src=RoutineName)

! Loop over all items in ntp, printing metadata and min/max if allocated
DO i =1, dim_ntp
  WRITE(umMessage,'(I4,2X,L1,9X,A)')                                           &
      i, all_ntp(i)%l_required, all_ntp(i)%varname
  CALL umPrint(umMessage,src=RoutineName)
  IF (ALLOCATED(all_ntp(i)%data_3d)) THEN
    WRITE(umMessage,'(A,2E12.4)') 'Min/max: ',                                 &
      MINVAL(all_ntp(i)%data_3d),MAXVAL(all_ntp(i)%data_3d)
    CALL umPrint(umMessage,src=routinename)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_all_ntp

! ----------------------------------------------------------------------
SUBROUTINE ntp_dealloc()
! ----------------------------------------------------------------------
! Description:
! Deallocate all data arrays in the all_ntp array.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Local variables

INTEGER :: i   ! Loop counter

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NTP_DEALLOC'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i = 1, SIZE(all_ntp)
  IF (ALLOCATED(all_ntp(i)%data_3d)) DEALLOCATE(all_ntp(i)%data_3d)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ntp_dealloc

! ----------------------------------------------------------------------
INTEGER FUNCTION name2ntpindex(varname)
! ----------------------------------------------------------------------
! Description:
! Given a variable name look up where it is in the given ntp table.
! Abort here if failed to look up.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Function arguments
CHARACTER(LEN=*), INTENT(IN) :: varname

! Local variables
INTEGER :: i
INTEGER :: errcode                                   ! error code
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NAME2NTPINDEX'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! set a default value. If it is still this at the end, the search has failed.
name2ntpindex = -999

! Search all entries in ntp array to find the one we want
ntp_loop: DO i = 1, SIZE(all_ntp)
  IF (all_ntp(i)%varname == varname) THEN
    name2ntpindex = i
    EXIT ntp_loop
  END IF
END DO ntp_loop

! If name2ntpindex is -999 then call ereport
IF (name2ntpindex == -999) THEN
  WRITE(cmessage,'(A,A,A)') 'Failed to find: ', varname, ' in NTP structure.'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END FUNCTION name2ntpindex
! ----------------------------------------------------------------------

END MODULE ukca_ntp_mod

