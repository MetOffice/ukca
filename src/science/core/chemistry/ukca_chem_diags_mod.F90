! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module to contain code to place chemistry time step diagnostics
!   into the section 50 STASHwork array
!
! Method:
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and
!  The Met Office. See www.ukca.ac.uk.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

MODULE ukca_chem_diags_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CHEM_DIAGS_MOD'

CONTAINS

SUBROUTINE ukca_chem_diags (                                                   &
! IN model dimensions
  row_length, rows, model_levels,                                              &
! IN fields for diagnostics
  nat_psc, trop_ch4_mol, trop_o3_mol, trop_oh_mol,                             &
  strat_ch4_mol, strat_ch4loss,                                                &
  atm_ch4_mol, atm_co_mol, atm_n2o_mol,                                        &
  atm_cf2cl2_mol, atm_cfcl3_mol, atm_mebr_mol,                                 &
  atm_h2_mol, so4_sa, H_plus_3d_arr,                                           &
! INOUT stash workspace
  len_stashwork, stashwork)

! Description:
!   To place chemistry time step diagnostics into the section 50
!   STASHwork array

USE ukca_um_legacy_mod, ONLY: copydiag_3d, len_stlist, stindex, stlist,        &
                              num_stash_levels, stash_levels, si, sf, si_last, &
                              stashcode_ukca_chem_diag,                        &
                              stashcode_ukca_nat,                              &
                              stashcode_ukca_trop_ch4,                         &
                              stashcode_ukca_trop_o3,                          &
                              stashcode_ukca_trop_oh,                          &
                              stashcode_ukca_strat_ch4,                        &
                              stashcode_ukca_strt_ch4_lss,                     &
                              stashcode_ukca_atmos_ch4,                        &
                              stashcode_ukca_atmos_co,                         &
                              stashcode_ukca_atmos_n2o,                        &
                              stashcode_ukca_atmos_cfc12,                      &
                              stashcode_ukca_atmos_cfc11,                      &
                              stashcode_ukca_atmos_ch3br,                      &
                              stashcode_ukca_atmos_h2,                         &
                              stashcode_ukca_so4_sad,                          &
                              stashcode_ukca_h_plus,                           &
                              ukca_diag_sect

USE ereport_mod,  ONLY: ereport
USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN) :: row_length        ! Model dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! Diagnostic tracers
! Nitric acid trihydrate (kg(nat)/kg(air))
REAL, INTENT(IN) :: nat_psc(row_length,rows,model_levels)
! Trop CH4 burden (moles)
REAL, INTENT(IN) :: trop_ch4_mol(row_length,rows,model_levels)
! Trop O3 burden (moles)
REAL, INTENT(IN) :: trop_o3_mol(row_length,rows,model_levels)
! Trop OH burden (moles)
REAL, INTENT(IN) :: trop_oh_mol(row_length,rows,model_levels)
! Strat CH4 burden (moles)
REAL, INTENT(IN) :: strat_ch4_mol(row_length,rows,model_levels)
! Strat CH4 loss (Moles/s)
REAL, INTENT(IN) :: strat_ch4loss(row_length,rows,model_levels)
! Atmospheric Burden of CH4 in moles
REAL, INTENT(IN) :: atm_ch4_mol(row_length,rows,model_levels)
! Atmospheric Burden of CO in moles
REAL, INTENT(IN) :: atm_co_mol(row_length,rows,model_levels)
! Atmospheric Burden of Nitrous Oxide (N2O) in moles
REAL, INTENT(IN) :: atm_n2o_mol(row_length,rows,model_levels)
! Atmospheric Burden of CFC-12 in moles
REAL, INTENT(IN) :: atm_cf2cl2_mol(row_length,rows,model_levels)
! Atmospheric Burden of CFC-11 in moles
REAL, INTENT(IN) :: atm_cfcl3_mol(row_length,rows,model_levels)
! Atmospheric Burden of CH3Br in moles
REAL, INTENT(IN) :: atm_mebr_mol(row_length,rows,model_levels)
! Atmospheric Burden of H2 in moles
REAL, INTENT(IN) :: atm_h2_mol(row_length,rows,model_levels)
! Aerosol surface area used in chemistry
REAL, INTENT(IN) :: so4_sa(row_length,rows,model_levels)
! H+ concentrations used for pH calculation in chemistry
REAL, INTENT(IN) :: H_plus_3d_arr(row_length,rows,model_levels)

! Diagnostics info
INTEGER, INTENT(IN) :: len_stashwork ! Length of diagnostics array
REAL, INTENT(IN OUT) :: stashwork(len_stashwork)  ! STASH workspace

! Local variables
INTEGER :: item                                ! STASH item
INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0  ! DrHook tracing entry
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1 ! DrHook tracing exit
INTEGER :: icode                               ! error code for EReport
INTEGER :: im_index                            ! internal model index

REAL(KIND=jprb) :: zhook_handle ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEM_DIAGS' ! used for EReport


! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode = 0 ! Initialise error status
im_index = 1

! ----------------------------------------------------------------------
!   Copy diagnostic information to STASHwork for STASH processing
! ----------------------------------------------------------------------
! DIAG.50218 Nitric acid trihydrate
! ----------------------------------------------------------------------
item = stashcode_ukca_nat - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       nat_psc(:,:,:),                                                         &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50220 Trop CH4 burden in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_trop_ch4 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       trop_ch4_mol (:,:,:),                                                   &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50221 Trop O3 burden in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_trop_o3 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       trop_o3_mol(:,:,:),                                                     &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50222 Trop OH burden in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_trop_oh - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       trop_oh_mol(:,:,:),                                                     &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50223 Strat CH4 burden in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_strat_ch4 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       strat_ch4_mol(:,:,:),                                                   &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50226 Strat CH4 loss
! ----------------------------------------------------------------------
item = stashcode_ukca_strt_ch4_lss - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       strat_ch4loss(:,:,:),                                                   &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50231 Atmospheric Burden of CH4 in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_ch4 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       atm_ch4_mol(:,:,:),                                                     &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50232 Atmospheric Burden of CO in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_co - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       atm_co_mol(:,:,:),                                                      &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50233 Atmospheric Burden of N2O in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_n2o - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       atm_n2o_mol(:,:,:),                                                     &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50234 Atmospheric Burden of CFC-12 in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_cfc12 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       atm_cf2cl2_mol(:,:,:),                                                  &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50235 Atmospheric Burden of CFC-11 in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_cfc11 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       atm_cfcl3_mol(:,:,:),                                                   &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50236 Atmospheric Burden of CH3Br in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_ch3br - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       atm_mebr_mol(:,:,:),                                                    &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50237 Atmospheric Burden of H2 in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_h2 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       atm_h2_mol(:,:,:),                                                      &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! --------------------------------------------------------------------
! DIAG.50256 Aerosol surface area density as used in heteogeneous
!            chemistry & photolysis
! --------------------------------------------------------------------
item = stashcode_ukca_so4_sad - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       so4_sa(:,:,:),                                                          &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)


! --------------------------------------------------------------------
! DIAG.50442 [H+] concentrations used in the calculation of pH values
!            that are used in chemistry and deposition
! --------------------------------------------------------------------
item = stashcode_ukca_h_plus - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index):                &
       si_last(item,UKCA_diag_sect,im_index)),                                 &
       H_plus_3d_arr(:,:,:),                                                   &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1)
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chem_diags

END MODULE ukca_chem_diags_mod
