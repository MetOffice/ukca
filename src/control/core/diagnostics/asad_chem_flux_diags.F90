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
!  Module to deal with the chemical flux diagnostics from ASAD
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
!
! ######################################################################
!
MODULE asad_chem_flux_diags

USE ukca_um_legacy_mod, ONLY: modl_b, isec_b, item_b, iopl_d,                  &
                              idom_b, ndiag,                                   &
                              st_levels_model_theta, st_levels_single,         &
                              submodel_for_sm, atmos_im, mype
USE asad_flux_dat,    ONLY: asad_chemical_fluxes, stashcode_ukca_chem_diag
USE asad_mod,         ONLY: advt, dpd, dpw, fpsc1, fpsc2,                      &
                            jpspb, jpsph, jpspj, jpspt,                        &
                            ldepd, ldepw,                                      &
                            nbrkx, nhrkx, nprkx, ntrkx, prk,                   &
                            spb, speci, sph, spj, spt, y,                      &
                            jpspec, jpbk, jptk, jppj, jphk, jpnr, jpctr
USE ukca_tropopause,  ONLY: L_stratosphere
USE ukca_cspecies,    ONLY: c_species
USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook
USE ereport_mod,      ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus,                         &
                      PrStatus_Oper, PrStatus_Diag, umPrintFlush
USE missing_data_mod,      ONLY: imdi
USE errormessagelength_mod, ONLY: errormessagelength
USE ukca_fieldname_mod, ONLY: maxlen_diagname

IMPLICIT NONE

PRIVATE

! public types
PUBLIC :: chemdiag
! public variables
PUBLIC :: asad_chemdiags
PUBLIC :: stash_handling
PUBLIC :: n_chemdiags
PUBLIC :: n_asad_diags
PUBLIC :: L_asad_use_chem_diags
PUBLIC :: L_asad_use_air_ems
PUBLIC :: L_asad_use_light_ems
PUBLIC :: L_asad_use_light_diags_tot
PUBLIC :: L_asad_use_light_diags_c2g
PUBLIC :: L_asad_use_light_diags_c2c
PUBLIC :: L_asad_use_light_diags_N
PUBLIC :: L_asad_use_volc_ems
PUBLIC :: L_asad_use_sulp_ems
PUBLIC :: L_asad_use_surf_ems
PUBLIC :: L_asad_use_flux_rxns
PUBLIC :: L_asad_use_rxn_rates
PUBLIC :: L_asad_use_wetdep
PUBLIC :: L_asad_use_drydep
PUBLIC :: L_asad_use_tendency
PUBLIC :: L_asad_use_mass_diagnostic
PUBLIC :: L_asad_use_psc_diagnostic
PUBLIC :: L_asad_use_trop_mask
PUBLIC :: L_asad_use_STE
PUBLIC :: L_asad_use_output_tracer
PUBLIC :: aircraft_emissions
PUBLIC :: lightning_emissions
PUBLIC :: total_lightning_flashes
PUBLIC :: lightning_flashes_cloud2ground
PUBLIC :: lightning_flashes_cloud2cloud
PUBLIC :: total_N_2D
PUBLIC :: volcanic_emissions
PUBLIC :: so2_emchar
PUBLIC :: calculate_tendency
PUBLIC :: calculate_STE
! public routines
PUBLIC :: asad_chemical_diagnostics
PUBLIC :: asad_setstash_chemdiag
PUBLIC :: asad_init_chemdiag
PUBLIC :: asad_emissions_diagnostics
PUBLIC :: asad_3d_emissions_diagnostics
PUBLIC :: asad_tendency_ste
PUBLIC :: asad_tropospheric_mask
PUBLIC :: asad_mass_diagnostic
PUBLIC :: asad_allocate_chemdiag
PUBLIC :: asad_psc_diagnostic
PUBLIC :: asad_output_tracer
PUBLIC :: asad_lightning_diagnostics

TYPE :: chemdiag
  INTEGER :: location          ! location of reaction in
                               ! y/rk/prk/em_field/dpd/dpw array
  INTEGER :: stash_number      ! Stash number
  INTEGER :: find_rxn_loc      ! Reaction location
  INTEGER :: num_reactants     ! No of reactants
  INTEGER :: num_products      ! No of products
  LOGICAL :: tropospheric_mask ! are we only outputting the troposphere?
  REAL    :: molecular_mass    ! needed for emissions, tendency, STE
  REAL    :: c_vmr_to_mmr      ! conversion factor from vmr to mmr
  INTEGER :: num_levs          ! No of levels
  LOGICAL :: can_deallocate    ! can the throughput array be deallocated ?
  CHARACTER(LEN=3) :: diag_type ! Type: RXN,DDP,WDP,EMS,NET,STE,TPM
  CHARACTER(LEN=1) :: rxn_type  ! RXN only: B,J,H,T
                                ! (DDP==D,WDP==W,EMS==E),X
  CHARACTER(LEN=10) :: species  ! For DDP, WDP, EMS, NET only
  ! fix number of reactants at 2, but allow number of products to change
  CHARACTER(LEN=10) :: reactants(1:2) ! for use with RXN only
  ! for use with RXN only
  CHARACTER(LEN=10), ALLOCATABLE :: products(:)
  REAL, ALLOCATABLE :: throughput(:,:,:)
END TYPE chemdiag

! Type for holding STASH codes and dimension sizes and index in
! asad_chemdiags array for requested fields
TYPE :: stashsumdiag
  CHARACTER(LEN=maxlen_diagname) :: diagname  ! Field name
  INTEGER :: stash_value               ! section*1000 + item
  INTEGER :: stash_section             ! section
  INTEGER :: stash_item                ! item
  INTEGER :: number_of_fields          ! No of fields
  INTEGER :: len_dim1                  ! Array dimension
  INTEGER :: len_dim2                  ! Array dimension
  INTEGER :: len_dim3                  ! Array dimension
  INTEGER, ALLOCATABLE :: chemdiags_location(:)
  LOGICAL :: l_um_legacy_request       ! True if this is a UM legacy-style
                                       ! request
END TYPE stashsumdiag

CHARACTER(LEN=10) :: BlankStr = '          '  ! Represents undefined string

! Array to define diagnostic properties
TYPE(chemdiag), ALLOCATABLE, SAVE :: asad_chemdiags(:)

TYPE(stashsumdiag), ALLOCATABLE, SAVE :: stash_handling(:)

INTEGER, SAVE :: n_asad_diags = 0 ! No. of ASAD diagnostics actually
                                  ! requested in STASH (initialised to 0 here
                                  ! to show that no ASAD diagnostics have been
                                  ! set up)

INTEGER, SAVE :: n_chemdiags      ! No. of entries from 'asad_chemical_fluxes'
                                  ! contributing to requested diagnostics

! Switch to indicate whether any ASAD diagnostics are requested.
LOGICAL, SAVE :: L_asad_use_chem_diags=.FALSE.

! Initially set to these - if not required then will be turned off.
LOGICAL, SAVE :: L_asad_use_air_ems=.FALSE.   ! For aircraft emission
LOGICAL, SAVE :: L_asad_use_light_ems=.FALSE. ! For lightning emission
LOGICAL, SAVE :: L_asad_use_surf_ems=.FALSE.  ! For surface emission
LOGICAL, SAVE :: L_asad_use_volc_ems=.FALSE.  ! For volcanic emission
LOGICAL, SAVE :: L_asad_use_sulp_ems=.FALSE.  ! For 3D SO2 emission
LOGICAL, SAVE :: L_asad_use_flux_rxns=.FALSE. ! For reaction flux
LOGICAL, SAVE :: L_asad_use_rxn_rates=.FALSE.
LOGICAL, SAVE :: L_asad_use_wetdep=.FALSE.    ! For wet deposition
LOGICAL, SAVE :: L_asad_use_drydep=.FALSE.    ! For dry deposition
LOGICAL, SAVE :: L_asad_use_tendency=.FALSE.  ! For tendency
LOGICAL, SAVE :: L_asad_use_STE=.FALSE.       ! For strat-trop exchange

! For air mass diagnostic
LOGICAL, SAVE :: L_asad_use_mass_diagnostic=.FALSE.

! For psc diagnostic
LOGICAL, SAVE :: L_asad_use_psc_diagnostic=.FALSE.

! To output tropospheric mask
LOGICAL, SAVE :: L_asad_use_trop_mask_output=.FALSE.

! To output tracer
LOGICAL, SAVE :: L_asad_use_output_tracer=.FALSE.

! For o/p in troposphere only
LOGICAL, SAVE :: L_asad_use_trop_mask=.FALSE.

LOGICAL, SAVE :: L_asad_use_light_diags_tot=.FALSE.
LOGICAL, SAVE :: L_asad_use_light_diags_c2g=.FALSE.
LOGICAL, SAVE :: L_asad_use_light_diags_c2c=.FALSE.
LOGICAL, SAVE :: L_asad_use_light_diags_N=.FALSE.


CHARACTER(LEN=3), PARAMETER :: cdrxn='RXN'  ! chars for reaction
CHARACTER(LEN=3), PARAMETER :: cddep='DEP'  ! chars for deposition
CHARACTER(LEN=3), PARAMETER :: cdems='EMS'  ! chars for emission
CHARACTER(LEN=3), PARAMETER :: cdnet='NET'  ! chars for tendency
CHARACTER(LEN=3), PARAMETER :: cdste='STE'  ! chars for strattrop exchange
CHARACTER(LEN=3), PARAMETER :: cdmas='MAS'  ! chars for air mass
CHARACTER(LEN=3), PARAMETER :: cdpsc='PSC'  ! chars for polar strat cloud
CHARACTER(LEN=3), PARAMETER :: cdtpm='TPM'  ! chars for tropospheric mask
CHARACTER(LEN=3), PARAMETER :: cdout='OUT'  ! chars for tracer output
CHARACTER(LEN=1), PARAMETER :: cdbimol='B'  ! char for bimolecular rxn
CHARACTER(LEN=1), PARAMETER :: cdtermol='T' ! char for termolecular rxn
CHARACTER(LEN=1), PARAMETER :: cdphot='J'   ! char for photolytic rxn
CHARACTER(LEN=1), PARAMETER :: cdhetero='H' ! char for heterogeneous rxn
CHARACTER(LEN=1), PARAMETER :: cdair='A'    ! char for aircraft emission
CHARACTER(LEN=1), PARAMETER :: cdsurf='S'   ! char for surface emission
CHARACTER(LEN=1), PARAMETER :: cdsulp='T'   ! char for 3D SO2 emission
CHARACTER(LEN=1), PARAMETER :: cdvolc='V'   ! char for volcanic emission
CHARACTER(LEN=1), PARAMETER :: cdlight='L'  ! char for lightning emission
CHARACTER(LEN=1), PARAMETER :: cdwet='W'    ! char for wet deposition
CHARACTER(LEN=1), PARAMETER :: cddry='D'    ! char for dry deposition
CHARACTER(LEN=1), PARAMETER :: cdpsc_typ1='1' ! char for PSC type 1
CHARACTER(LEN=1), PARAMETER :: cdpsc_typ2='2' ! char for PSC type 2
CHARACTER(LEN=1), PARAMETER :: so2_emchar='S' !    "         "
CHARACTER(LEN=3), PARAMETER :: calculate_tendency=cdnet
CHARACTER(LEN=3), PARAMETER :: calculate_STE=cdste
CHARACTER(LEN=3), PARAMETER :: cdrte='RTE'
CHARACTER(LEN=3), PARAMETER :: cdlgt='LGT'
CHARACTER(LEN=1), PARAMETER :: cdlgttot='T'
CHARACTER(LEN=1), PARAMETER :: cdlgtc2g='G'
CHARACTER(LEN=1), PARAMETER :: cdlgtc2c='C'
CHARACTER(LEN=1), PARAMETER :: cdlgtN='N'
CHARACTER(LEN=1), PARAMETER :: cdignore='X'
CHARACTER(LEN=1), PARAMETER :: aircraft_emissions='A'
CHARACTER(LEN=1), PARAMETER :: lightning_emissions='L'
CHARACTER(LEN=1), PARAMETER :: volcanic_emissions='V'
CHARACTER(LEN=1), PARAMETER :: total_lightning_flashes=cdlgttot
CHARACTER(LEN=1), PARAMETER :: lightning_flashes_cloud2ground=cdlgtc2g
CHARACTER(LEN=1), PARAMETER :: lightning_flashes_cloud2cloud=cdlgtc2c
CHARACTER(LEN=1), PARAMETER :: total_N_2D=cdlgtN



INTEGER, PARAMETER :: number_of_reactants=2       ! No of reactants

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ASAD_CHEM_FLUX_DIAGS'

CONTAINS

! #####################################################################

SUBROUTINE asad_setstash_chemdiag(row_length, rows, model_levels,              &
                                  diag_requests, master_diag_list)
! Subroutine to collate some basic information about each requested ASAD
! diagnostic in the array 'stash_handling' to support the processing of
! the field for output by the routine 'asad_flux_put_stash'.
! Each ASAD diagnostic is identified by a unique code, corresponding to
! to the STASH code used for that diagnostic in the UM. A single
! diagnostic may be the sum of multiple fluxes having the same STASH
! code in its definition.
! Note that the number of individual fluxes contributing to each
! diagnostic and the location of metadata in the 'asad_chemdiags' array
! defining these fluxes will be added subsequently when that array is
! set up in the routine 'asad_init_chemdiag'.

USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_diagnostics_type_mod, ONLY: diag_entry_type, diag_requests_type,      &
                                     n_diag_group, diag_status_unavailable

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN) :: row_length     ! length of row
INTEGER, INTENT(IN) :: rows           ! number of rows
INTEGER, INTENT(IN) :: model_levels   ! number of levels

TYPE(diag_requests_type), INTENT(IN) :: diag_requests(:)
                                      ! Diagnostic request info by group

TYPE(diag_entry_type), INTENT(IN) :: master_diag_list(:)
                                      ! Master diagnostics list

! Type for temporary handling of STASH item code and other info specific to
! each requested field.
TYPE :: field_spec_type
  CHARACTER(LEN=maxlen_diagname) :: diagname
  INTEGER :: item
  INTEGER :: len_dim3
  LOGICAL :: l_um_legacy_request
END TYPE field_spec_type

! Local variables

INTEGER :: group
INTEGER :: n_requests
INTEGER :: i_master
INTEGER :: i
INTEGER :: ierr
INTEGER :: idiag
INTEGER :: stash_value
INTEGER :: item

LOGICAL :: l_done(999)   ! True if STASH item has already been processed

CHARACTER (LEN=errormessagelength) :: cmessage          ! Error message

! Array to define STASH codes and sizes of requested fields
TYPE(field_spec_type), ALLOCATABLE, SAVE :: chemcodes(:)

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_SETSTASH_CHEMDIAG'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr = -1

stash_value = 0
n_asad_diags = 0

! Allocate temporary array to total number of diagnostics defined in
! 'asad_chemical_fluxes'. This latter array defines all the chemical fluxes
! potentially contributing to ASAD diagnostics and contains duplicate STASH
! codes, so the number of ASAD diagnostics requested should always be smaller
! than the allocated size.
ALLOCATE(chemcodes(SIZE(asad_chemical_fluxes)))
chemcodes(:)%diagname = ''
chemcodes(:)%item = imdi
chemcodes(:)%len_dim3 = imdi
chemcodes(:)%l_um_legacy_request = .FALSE.

! Loop over the diagnostic requests below to determine number of ASAD
! diagnostics requested (i.e. those that are also defined in the
! 'asad_chemical_fluxes' array) and collate info for each.

! --- Process diagnostic requests made via the UKCA API ---

! Note that the same diagnostic may be requested in different groups but only
! one ASAD diagnostic will be set up for each diagnostic (as identified by
! name or stash value).

l_done(:) = .FALSE.

DO group = 1, n_diag_group

  IF (ALLOCATED(diag_requests(group)%varnames)) THEN
    n_requests = SIZE(diag_requests(group)%varnames)
  ELSE
    n_requests = 0
  END IF

  DO i = 1, n_requests
    IF (.NOT. diag_requests(group)%status_flags(i) ==                          &
              diag_status_unavailable) THEN
      ! Requested diagnostic is available in the present UKCA configuration
      ! (The request may be active/inactive)
      i_master = diag_requests(group)%i_master(i)
      stash_value = master_diag_list(i_master)%asad_id
      IF (stash_value > 0) THEN
        ! Diagnostic is an ASAD diagnostic so add to temporary array
        ! if not already done
        item = stash_value - (stashcode_ukca_chem_diag * 1000)
        IF (.NOT. l_done(item)) THEN
          l_done(item) = .TRUE.
          ! Check present in asad_chemical_fluxes array for safety before adding
          IF (ANY(asad_chemical_fluxes(:)%stash_number == stash_value)) THEN
            n_asad_diags = n_asad_diags + 1
            chemcodes(n_asad_diags)%diagname = diag_requests(group)%varnames(i)
            chemcodes(n_asad_diags)%item = item
            ! Note: all ASAD diagnostics currently supported by the API
            ! are 3D (otherwise we would need to check 'rxn_type' in the
            ! 'asad_chemical_fluxes' entry to determine dimensionality)
            chemcodes(n_asad_diags)%len_dim3 = model_levels
          ELSE
            ierr = stash_value
            cmessage = 'Expected STASH code not found in asad_chemical_fluxes'
            CALL ereport(ModuleName//':'//RoutineName,ierr,cmessage)
          END IF
        END IF
      END IF
    END IF
  END DO

END DO

! --- For UM parent, also process UM legacy-style requests ---

! Note that legacy-style requests are ignored for any STASH items that have
! already been requested via the API. Also each ASAD request is recorded just
! once, even though there may be multiple instances in STASH.

IF (ukca_config%l_enable_diag_um) THEN
  DO idiag = 1, ndiag
    IF (modl_b(idiag) == submodel_for_sm(atmos_im) .AND.                       &
        isec_b(idiag) == stashcode_ukca_chem_diag) THEN
      ! Requested diagnostic is in chemistry section so process here if not
      ! already done
      item = item_b(idiag)
      IF (.NOT. l_done(item)) THEN
        l_done(item) = .TRUE.
        stash_value = stashcode_ukca_chem_diag * 1000 + item
        IF (ANY(asad_chemical_fluxes(:)%stash_number == stash_value)) THEN
          ! Diagnostic is an available ASAD diagnostic so add to temporary array
          n_asad_diags = n_asad_diags + 1
          chemcodes(n_asad_diags)%item = item
          ! check on size of model levels - currently only
          ! contiguous theta levels or surface allowed
          IF (iopl_d(idom_b(idiag)) == st_levels_model_theta) THEN
            chemcodes(n_asad_diags)%len_dim3 = model_levels
          ELSE IF (iopl_d(idom_b(idiag)) == st_levels_single) THEN
            chemcodes(n_asad_diags)%len_dim3 = 1
          ELSE
            WRITE(cmessage,'(A,5(i6,1X),A)') 'ERROR: ',                        &
                  idiag, modl_b(idiag),                                        &
                  isec_b(idiag), item_b(idiag),                                &
                  iopl_d(idom_b(idiag)),                                       &
                  ' incorrect model levels'
            ierr = isec_b(idiag)*1000 + item_b(idiag)
            CALL ereport(ModuleName//':'//RoutineName,ierr,cmessage)
          END IF
          chemcodes(n_asad_diags)%l_um_legacy_request = .TRUE.
        END IF
      END IF
    END IF
  END DO
END IF

l_asad_use_chem_diags = (n_asad_diags > 0)

IF (l_asad_use_chem_diags) THEN

  ! Copy STASH codes and dimension sizes from temporary array to an array
  ! of the correct length

  ALLOCATE(stash_handling(1:n_asad_diags))
  stash_handling(:)%number_of_fields = 0

  DO i = 1, n_asad_diags
    stash_handling(i)%diagname = chemcodes(i)%diagname
    stash_handling(i)%stash_section = stashcode_ukca_chem_diag
    stash_handling(i)%stash_item = chemcodes(i)%item
    stash_handling(i)%stash_value =                                            &
          (stash_handling(i)%stash_section * 1000) +                           &
          stash_handling(i)%stash_item
    stash_handling(i)%len_dim1 = row_length
    stash_handling(i)%len_dim2 = rows
    stash_handling(i)%len_dim3 = chemcodes(i)%len_dim3
    stash_handling(i)%l_um_legacy_request = chemcodes(i)%l_um_legacy_request
  END DO

END IF

DEALLOCATE(chemcodes)

IF (PrintStatus >= Prstatus_Oper) THEN
  DO i=1,n_asad_diags
    WRITE(umMessage,'(A,6(I5,1x))') 'ASAD_FLUXES: ',i,                         &
          stash_handling(i)%stash_section, stash_handling(i)%stash_item,       &
          stash_handling(i)%len_dim1, stash_handling(i)%len_dim2,              &
          stash_handling(i)%len_dim3
    CALL umPrint(umMessage,src=RoutineName)
  END DO
END IF

ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_setstash_chemdiag

! #####################################################################

SUBROUTINE asad_zero_chemdiag(diagnostic)
IMPLICIT NONE


TYPE(chemdiag), INTENT(IN OUT) :: diagnostic

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_ZERO_CHEMDIAG'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     "zero" all fields
diagnostic%location             = 0
diagnostic%stash_number         = 0
diagnostic%find_rxn_loc         = 0
diagnostic%num_reactants        = 0
diagnostic%num_products         = 0
diagnostic%num_levs             = 0
diagnostic%can_deallocate       = .FALSE.
diagnostic%tropospheric_mask    = .FALSE.
diagnostic%molecular_mass       = 0.0
diagnostic%c_vmr_to_mmr         = 0.0
diagnostic%species              = BlankStr
diagnostic%diag_type            = BlankStr(1:3)
diagnostic%rxn_type             = BlankStr(1:1)
diagnostic%reactants(:)         = BlankStr
IF (ALLOCATED(diagnostic%products))                                            &
     diagnostic%products(:)          = BlankStr
IF (ALLOCATED(diagnostic%throughput))                                          &
     diagnostic%throughput(:,:,:)    = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_zero_chemdiag
! #####################################################################
SUBROUTINE asad_init_chemdiag()

IMPLICIT NONE

INTEGER :: i,j,k,l

! counters
INTEGER :: icount,jdiag,istash

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=72) :: outformat
INTEGER       :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_INIT_CHEMDIAG'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

n_chemdiags = 0
DO i = 1, n_asad_diags
  ! calculate how many diagnostics from flux_dat we are using
  DO j=1,SIZE(asad_chemical_fluxes)
    IF (stash_handling(i)%stash_value ==                                       &
               asad_chemical_fluxes(j)%stash_number) THEN
      n_chemdiags = n_chemdiags + 1
      ! NOTE: do not exit after first check, since we may have
      !       more than one diagnostic associated with each
      !       stash number.
    END IF
  END DO
END DO        ! i=1,n_asad_diags

! initially set to these - if not required then will be turned off.
! set up above, but also re-iterate here.
L_asad_use_air_ems         = .FALSE.
L_asad_use_light_ems       = .FALSE.
L_asad_use_surf_ems        = .FALSE.
L_asad_use_volc_ems        = .FALSE.
L_asad_use_sulp_ems        = .FALSE.
L_asad_use_flux_rxns       = .FALSE.
L_asad_use_wetdep          = .FALSE.
L_asad_use_drydep          = .FALSE.
L_asad_use_tendency        = .FALSE.
L_asad_use_mass_diagnostic = .FALSE.
L_asad_use_trop_mask       = .FALSE.

! allocate diagnostics array
ALLOCATE(asad_chemdiags(1:n_chemdiags))
DO l=1,n_chemdiags
  CALL asad_zero_chemdiag(asad_chemdiags(l))
END DO

! Now we do a rather complicated do-loop to set up the asad_chemdiags array
jdiag = 0
DO istash=1,n_asad_diags
  DO j=1,SIZE(asad_chemical_fluxes)
    ! Check for requested stash code whether
    ! we have a corresponding flux specified
    IF (stash_handling(istash)%stash_value ==                                  &
        asad_chemical_fluxes(j)%stash_number) THEN

      ! increment jdiag. We may have more than one
      ! diagnostic associated with each stash number
      jdiag = jdiag + 1
      ! STASH number
      asad_chemdiags(jdiag)%stash_number = asad_chemical_fluxes(j)%stash_number
      ! diagnostic type
      asad_chemdiags(jdiag)%diag_type =                                        &
        TRIM(asad_cd_uppercase(asad_chemical_fluxes(j)%diag_type))
      ! reaction (or otherwise) type
      asad_chemdiags(jdiag)%rxn_type =                                         &
        TRIM(asad_cd_uppercase(asad_chemical_fluxes(j)%rxn_type))
      ! do we use a tropospheric mask only output the
      ! troposphere?
      IF (asad_chemical_fluxes(j)%tropospheric_mask) THEN
        asad_chemdiags(jdiag)%tropospheric_mask = .TRUE.
        ! set logical for global value
        IF (.NOT. L_asad_use_trop_mask) L_asad_use_trop_mask = .TRUE.
      ELSE
        asad_chemdiags(jdiag)%tropospheric_mask = .FALSE.
      END IF
      ! if we have two reactions with the exactly the same
      ! reactants and products, but different rates,
      ! find_reaction will only find the first one, this
      ! flag over-rides that
      IF (asad_chemical_fluxes(j)%rxn_location <= 0) THEN
        asad_chemdiags(jdiag)%find_rxn_loc = 1
      ELSE
        asad_chemdiags(jdiag)%find_rxn_loc =                                   &
          asad_chemical_fluxes(j)%rxn_location
      END IF
      ! fill in number of reactants and products
      IF (asad_chemical_fluxes(j)%num_species == 1) THEN
        ! are not dealing with a reaction
        asad_chemdiags(jdiag)%num_reactants = 0
        asad_chemdiags(jdiag)%num_products = 0
      ELSE IF (asad_chemical_fluxes(j)%num_species >= 3) THEN
        ! we have a reaction (probably, or user error!)
        asad_chemdiags(jdiag)%num_reactants = number_of_reactants ! 2 currently
        asad_chemdiags(jdiag)%num_products =                                   &
           asad_chemical_fluxes(j)%num_species -                               &
           asad_chemdiags(jdiag)%num_reactants
      ELSE ! incorrect numbers specified!
        WRITE(cmessage,'(A,I0,A)')                                             &
             'WRONG NUMBER OF SPECIES: ',                                      &
             asad_chemical_fluxes(j)%num_species, ' FOUND'
        errcode=asad_chemdiags(jdiag)%stash_number
        WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',             &
             errcode,' PE: ',mype
        CALL umPrint(umMessage,src=RoutineName)
        CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
      END IF

      ! initialise reactants and products
      asad_chemdiags(jdiag)%reactants(:) = BlankStr
      IF (ALLOCATED(asad_chemdiags(jdiag)%products))                           &
        DEALLOCATE(asad_chemdiags(jdiag)%products)
      ALLOCATE(asad_chemdiags(jdiag)%products(1:                               &
               asad_chemdiags(jdiag)%num_products))
      asad_chemdiags(jdiag)%products(:) = BlankStr

      ! fill in from list - rxn specific
      SELECT CASE (asad_chemdiags(jdiag)%diag_type)
      CASE (cdrxn,cdrte)
        ! copy over reacants
        DO i=1,number_of_reactants
          asad_chemdiags(jdiag)%reactants(i) =                                 &
            TRIM(ADJUSTL(asad_chemical_fluxes(j)%reactants(i)))
        END DO
        ! copy over products
        DO i=1,asad_chemdiags(jdiag)%num_products
          asad_chemdiags(jdiag)%products(i) =                                  &
            TRIM(ADJUSTL(asad_chemical_fluxes(j)%products(i)))
        END DO
      CASE (cddep,cdems,cdnet,cdste,cdout)
        ! only one species in this case
        asad_chemdiags(jdiag)%species =                                        &
          TRIM(ADJUSTL(asad_chemical_fluxes(j)%reactants(1)))
      CASE (cdpsc,cdmas,cdtpm,cdlgt)
        ! no need to do anything - no species
      CASE DEFAULT
        cmessage='DIAGNOSTIC TYPE '// asad_chemdiags(jdiag)%diag_type //       &
                 ' NOT FOUND'
        errcode=asad_chemdiags(jdiag)%stash_number
        WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',             &
             errcode,' PE: ',mype
        CALL umPrint(umMessage,src=RoutineName)
        CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
      END SELECT

      ! Check on emissions, deposition etc, and set allocate/output logicals
      SELECT CASE (asad_chemdiags(jdiag)%diag_type)
      CASE (cdems)
         ! emissions only on chemical timesteps, and can
         ! deallocate at other times
        asad_chemdiags(jdiag)%can_deallocate = .TRUE.
        SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
        CASE (cdair) ! aircraft
          L_asad_use_air_ems=.TRUE.
        CASE (cdlight) ! lightning
          L_asad_use_light_ems=.TRUE.
        CASE (cdvolc) ! volcanic
          L_asad_use_volc_ems=.TRUE.
        CASE (cdsulp) ! volcanic
          L_asad_use_sulp_ems = .TRUE.
        CASE (cdsurf) ! surface
          L_asad_use_surf_ems=.TRUE.
        CASE DEFAULT
          cmessage='INCORRECT EMISSION TYPE: '//                               &
                   asad_chemdiags(jdiag)%diag_type//                           &
                   ':'//asad_chemdiags(jdiag)%rxn_type
          errcode=asad_chemdiags(jdiag)%stash_number
          WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',           &
               errcode,' PE: ',mype
          CALL umPrint(umMessage,src=RoutineName)
          CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
        END SELECT ! asad_chemdiags(jdiag)%rxn_type
      CASE (cdrxn)
        ! reactions only on chemical timesteps, and can
        ! deallocate at other times
        asad_chemdiags(jdiag)%can_deallocate = .TRUE.
        SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
        CASE (cdbimol,cdtermol,cdphot,cdhetero)
          L_asad_use_flux_rxns=.TRUE.
        CASE DEFAULT
          cmessage='INCORRECT REACTION TYPE: '//                               &
                   asad_chemdiags(jdiag)%diag_type//':'//                      &
                   asad_chemdiags(jdiag)%rxn_type
          errcode=asad_chemdiags(jdiag)%stash_number
          WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',           &
               errcode,' PE: ',mype
          CALL umPrint(umMessage,src=RoutineName)
          CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
        END SELECT ! asad_chemdiags(jdiag)%rxn_type
      CASE (cdrte)
         ! rates only on chemical timesteps, and can
         ! deallocate at other times
        asad_chemdiags(jdiag)%can_deallocate=.TRUE.
        SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
        CASE (cdbimol,cdtermol,cdphot,cdhetero)
          L_asad_use_rxn_rates=.TRUE.
        CASE DEFAULT
          cmessage='INCORRECT RATE TYPE: '//asad_chemdiags(jdiag)%diag_type//  &
                   ':'//asad_chemdiags(jdiag)%rxn_type
          errcode=asad_chemdiags(jdiag)%stash_number
          WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',           &
               errcode,' PE: ',mype
          CALL umPrint(umMessage,src=RoutineName)
          CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
        END SELECT ! asad_chemdiags(jdiag)%rxn_type
      CASE (cddep)
         ! deposition only on chemical timesteps, and can
         ! deallocate at other times
        asad_chemdiags(jdiag)%can_deallocate = .TRUE.
        SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
        CASE (cdwet)
          L_asad_use_wetdep=.TRUE.
        CASE (cddry)
          L_asad_use_drydep=.TRUE.
        CASE DEFAULT
          cmessage='INCORRECT DEPOSTITION TYPE: '//                            &
                   asad_chemdiags(jdiag)%diag_type//':'//                      &
                   asad_chemdiags(jdiag)%rxn_type
          errcode=asad_chemdiags(jdiag)%stash_number
          WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',           &
               errcode,' PE: ',mype
          CALL umPrint(umMessage,src=RoutineName)
          CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
        END SELECT ! asad_chemdiags(jdiag)%rxn_type
      CASE (cdnet)
        ! tendency only on chemical timesteps, and can
        ! deallocate at other times
        asad_chemdiags(jdiag)%can_deallocate = .TRUE.
        ! tendency
        L_asad_use_tendency = .TRUE.
      CASE (cdste)
        ! STE on all timesteps, and cannot be
        ! deallocated at any time
        asad_chemdiags(jdiag)%can_deallocate = .FALSE.
        ! tendency
        L_asad_use_STE = .TRUE.
      CASE (cdmas)
        ! air mass on all timesteps, and can
        ! deallocate at other times
        asad_chemdiags(jdiag)%can_deallocate = .TRUE.
        L_asad_use_mass_diagnostic=.TRUE.
      CASE (cdpsc)
        ! air mass on all timesteps, and can
        ! deallocate at other times
        asad_chemdiags(jdiag)%can_deallocate = .TRUE.
        L_asad_use_psc_diagnostic=.TRUE.
      CASE (cdtpm)
        ! Tropospheric mask on all timesteps, and can
        ! deallocate at other times
        asad_chemdiags(jdiag)%can_deallocate = .TRUE.
        ! set logical value to T for calculation of
        !tropospheric mask
        L_asad_use_trop_mask = .TRUE.
        ! set logical to true for outputting trop mask
        L_asad_use_trop_mask_output = .TRUE.
      CASE (cdout)
        ! tracer on all timesteps, and cannot be
        ! deallocated at any time
        asad_chemdiags(jdiag)%can_deallocate = .TRUE.
        ! tracer
        L_asad_use_output_tracer = .TRUE.
      CASE (cdlgt)
        ! lightning diagnostics on all timesteps, and cannot
        ! be deallocated at any time
        asad_chemdiags(jdiag)%can_deallocate = .TRUE.
        ! check which diagnostics
        SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
        CASE (cdlgttot)
          L_asad_use_light_diags_tot = .TRUE.
        CASE (cdlgtc2g)
          L_asad_use_light_diags_c2g = .TRUE.
        CASE (cdlgtc2c)
          L_asad_use_light_diags_c2c = .TRUE.
        CASE (cdlgtN)
          L_asad_use_light_diags_N = .TRUE.
        CASE DEFAULT
          cmessage='INCORRECT LIGHTNING DIAG TYPE: '//                         &
                   asad_chemdiags(jdiag)%diag_type//':'//                      &
                   asad_chemdiags(jdiag)%rxn_type
          errcode=asad_chemdiags(jdiag)%stash_number
          WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',           &
               errcode,' PE: ',mype
          CALL umPrint(umMessage,src=RoutineName)
          CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
        END SELECT
      END SELECT ! asad_chemdiags(jdiag)%diag_type

    END IF ! stash_handling(istash)%stash_value ==
           ! asad_chemical_fluxes(j)%stash_number
  END DO !j,SIZE(asad_chemical_fluxes)
END DO ! istash,n_asad_diags

! Check on summing options - simplifies STASH output
DO k=1,n_asad_diags
  DO l=1,n_chemdiags
     ! calcs number of fields with the same stash code
    IF (stash_handling(k)%stash_value ==                                       &
         asad_chemdiags(l)%stash_number) THEN
      stash_handling(k)%number_of_fields =                                     &
           stash_handling(k)%number_of_fields + 1
    END IF
  END DO
  ! creates an array of this size
  ALLOCATE(stash_handling(k)%chemdiags_location(1:                             &
       stash_handling(k)%number_of_fields))
  stash_handling(k)%chemdiags_location(:) = 0
  ! allocate %throughput array for each flux
  icount=0
  DO l=1,n_chemdiags
    IF (stash_handling(k)%stash_value ==                                       &
         asad_chemdiags(l)%stash_number) THEN
      icount=icount+1
      stash_handling(k)%chemdiags_location(icount) = l
      ! get number of levels (needed in putstash)
      asad_chemdiags(l)%num_levs =                                             &
           stash_handling(k)%len_dim3
    END IF
  END DO           ! l,n_chemdiags
END DO              ! k,n_asad_diags

! error check
DO k=1,n_asad_diags
  IF (stash_handling(k)%number_of_fields == 0) THEN
    cmessage=' STASH CODE NOT FOUND IN'//                                      &
             ' stash_handling ARRAY'
    errcode = -1*stash_handling(k)%stash_value
    WRITE(umMessage,'(2A,I0,A,I0,A,I0)') cmessage, ' Error code: ',            &
         errcode,' PE: ',mype,' K: ',k
    CALL umPrint(umMessage,src=RoutineName)
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
END DO

! print out diagnostics requested
IF ((PrintStatus >= PrStatus_Oper ) .AND. (mype == 0)) THEN
  WRITE(umMessage,'(A,A)') 'THE FOLLOWING UKCA FLUX DIAGNOSTICS ',             &
       'HAVE BEEN REQUESTED'
  CALL umPrint(umMessage,src=RoutineName)
  DO l=1,n_chemdiags
    IF (asad_chemdiags(l)%num_products > 0) THEN
       ! fancy formatting to cope with variable number of prods
      IF (asad_chemdiags(l)%num_products < 10) THEN
        WRITE(outformat,FMT='(A38,I1,A7)')                                     &
          '(I5,1X,A3,1X,A1,1X,A,1X,2(A,1X),A2,1X,',                            &
          asad_chemdiags(l)%num_products,                                      &
          '(A,1X))'
      ELSE
        WRITE(outformat,FMT='(A38,I2,A7)')                                     &
             '(I5,1X,A3,1X,A1,1X,A,1X,2(A,1X),A2,1X,',                         &
             asad_chemdiags(l)%num_products,                                   &
             '(A,1X))'
      END IF
      WRITE(umMessage,FMT=TRIM(ADJUSTL(outformat)))                            &
           asad_chemdiags(l)%stash_number,                                     &
           asad_chemdiags(l)%diag_type,                                        &
           asad_chemdiags(l)%rxn_type,                                         &
           TRIM(ADJUSTL(asad_chemdiags(l)%species)),                           &
           (TRIM(ADJUSTL(asad_chemdiags(l)%reactants(j))),                     &
           j=1,2),'->',                                                        &
           (TRIM(ADJUSTL(asad_chemdiags(l)%products(j))),                      &
           j=1,asad_chemdiags(l)%num_products)
      CALL umPrint(umMessage,src=RoutineName)
    ELSE
       ! not a reaction (or rather, has no products)
      WRITE(umMessage,FMT='(I5,1X,A3,1X,A1,1X,A,1X,2(A,1X))')                  &
           asad_chemdiags(l)%stash_number,                                     &
           asad_chemdiags(l)%diag_type,                                        &
           asad_chemdiags(l)%rxn_type,                                         &
           TRIM(ADJUSTL(asad_chemdiags(l)%species)),                           &
           (TRIM(ADJUSTL(asad_chemdiags(l)%reactants(j))),                     &
           j=1,2)
      CALL umPrint(umMessage,src=RoutineName)
    END IF
  END DO
END IF

! print out status of logicals
IF (PrintStatus == Prstatus_Diag) THEN
  WRITE(umMessage,'(A,L1)') 'L_asad_use_air_ems         = ',L_asad_use_air_ems
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,L1)') 'L_asad_use_light_ems       = ',L_asad_use_light_ems
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,L1)') 'L_asad_use_surf_ems        = ',L_asad_use_surf_ems
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,L1)') 'L_asad_use_volc_ems        = ',L_asad_use_volc_ems
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,L1)') 'L_asad_use_sulp_ems        = ',L_asad_use_sulp_ems
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,L1)') 'L_asad_use_flux_rxns       = ',L_asad_use_flux_rxns
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,L1)') 'L_asad_use_wetdep          = ',L_asad_use_wetdep
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,L1)') 'L_asad_use_drydep          = ',L_asad_use_drydep
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,L1)') 'L_asad_use_tendency        = ',L_asad_use_tendency
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,L1)') 'L_asad_use_mass_diagnostic = ',                   &
                            L_asad_use_mass_diagnostic
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,L1)') 'L_asad_use_trop_mask       = ',L_asad_use_trop_mask
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,L1)') 'L_asad_use_chem_diags      = ',                   &
                            L_asad_use_chem_diags
  CALL umPrint(umMessage,src=RoutineName)
END IF ! PrintStatus = Prstatus_diag

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_init_chemdiag

! #####################################################################
FUNCTION asad_cd_uppercase(string)
! To convert any string into uppercase letters

IMPLICIT NONE


CHARACTER(LEN=*), INTENT(IN) :: string
CHARACTER(LEN=30)            :: asad_cd_uppercase
INTEGER :: ia,iz,ishift
INTEGER :: i,ic

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_CD_UPPERCASE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! required for shift
ia     = ICHAR('a')
iz     = ICHAR('z')
ishift = ICHAR('A')-ia
! will need to strip leading dashes

asad_cd_uppercase = string

!uppercase the string a letter at a time using the shift values calculated above
DO i=1,30
  ic = ICHAR(asad_cd_uppercase(i:i))
  IF ((ic >= ia) .AND. (ic <= iz)) asad_cd_uppercase(i:i) =                    &
       CHAR(ishift+ic)
END DO

asad_cd_uppercase = ADJUSTL(asad_cd_uppercase)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION asad_cd_uppercase

! #####################################################################
SUBROUTINE asad_allocate_chemdiag(row_length,rows)
IMPLICIT NONE


INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows

INTEGER :: l


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_ALLOCATE_CHEMDIAG'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate flux array and initialise
DO l=1,n_chemdiags
  IF (.NOT. ALLOCATED(asad_chemdiags(l)%throughput)) THEN
    ALLOCATE(asad_chemdiags(l)%throughput(                                     &
             1:row_length,                                                     &
             1:rows,                                                           &
             1:asad_chemdiags(l)%num_levs))
    asad_chemdiags(l)%throughput(:,:,:) = 0.0
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_allocate_chemdiag

! #####################################################################
SUBROUTINE asad_chemical_diagnostics_init(ierr)

IMPLICIT NONE

INTEGER, INTENT(IN OUT) :: ierr          ! error code

INTEGER :: l,j,idep                     ! counters

INTEGER, ALLOCATABLE :: findrxn_tmp(:)   ! array to hold locations

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=72) :: outformat
INTEGER       :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_CHEMICAL_DIAGNOSTICS_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up locations of reactions using (an altered version of)
! asad_findreaction
DO l=1,n_chemdiags
  SELECT CASE (asad_chemdiags(l)%diag_type)
  CASE (cdrxn,cdrte)
    SELECT CASE (asad_chemdiags(l)%rxn_type)
    CASE (cdbimol) ! bimolecular reacions/rate
      ALLOCATE(findrxn_tmp(1:(jpbk+1)))
      findrxn_tmp = cd_findreaction(                                           &
           asad_chemdiags(l)%num_products,                                     &
           asad_chemdiags(l)%reactants,                                        &
           asad_chemdiags(l)%products,                                         &
           spb, nbrkx, (jpbk+1), jpspb )
      asad_chemdiags(l)%location =                                             &
           findrxn_tmp(asad_chemdiags(l)%find_rxn_loc)
      DEALLOCATE(findrxn_tmp)
      IF (asad_chemdiags(l)%location <= 0) THEN
        IF (asad_chemdiags(l)%num_products < 10) THEN
          WRITE(outformat,FMT='(A28,I1,A7)')                                   &
               '(A13,1X,A1,1X,2(A,1X),A2,1X,',                                 &
               asad_chemdiags(l)%num_products,                                 &
               '(A,1X))'
        ELSE
          WRITE(outformat,FMT='(A28,I2,A7)')                                   &
               '(A13,1X,A1,1X,2(A,1X),A2,1X,',                                 &
               asad_chemdiags(l)%num_products,                                 &
               '(A,1X))'
        END IF
        WRITE(umMessage,'(A)') outformat
        CALL umPrint(umMessage,src=RoutineName)
        WRITE(cmessage,FMT=TRIM(ADJUSTL(outformat)))                           &
          'RXN NOT FOUND', cdbimol,                                            &
          (TRIM(ADJUSTL(asad_chemdiags(l)%reactants(j))),                      &
           j=1,2),'->',                                                        &
          (TRIM(ADJUSTL(asad_chemdiags(l)%products(j))),                       &
           j=1,asad_chemdiags(l)%num_products)
        errcode=asad_chemdiags(l)%stash_number
        WRITE(umMessage,'(A,A13,I5,A5,I6)')                                    &
          cmessage, ' Error code: ',errcode, ' PE: ',mype
        CALL umPrint(umMessage,src=RoutineName)
        CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
      END IF
    CASE (cdhetero) ! Heterogeneous reactions/rates
      ALLOCATE(findrxn_tmp(1:(jphk+1)))
      findrxn_tmp = cd_findreaction(                                           &
           asad_chemdiags(l)%num_products,                                     &
           asad_chemdiags(l)%reactants,                                        &
           asad_chemdiags(l)%products,                                         &
           sph, nhrkx, (jphk+1), jpsph )
      asad_chemdiags(l)%location =                                             &
           findrxn_tmp(asad_chemdiags(l)%find_rxn_loc)
      DEALLOCATE(findrxn_tmp)
      IF (asad_chemdiags(l)%location <= 0) THEN
        IF (asad_chemdiags(l)%num_products < 10) THEN
          WRITE(outformat,FMT='(A28,I1,A7)')                                   &
               '(A13,1X,A1,1X,2(A,1X),A2,1X,',                                 &
               asad_chemdiags(l)%num_products,                                 &
               '(A,1X))'
        ELSE
          WRITE(outformat,FMT='(A28,I2,A7)')                                   &
               '(A13,1X,A1,1X,2(A,1X),A2,1X,',                                 &
               asad_chemdiags(l)%num_products,                                 &
               '(A,1X))'
        END IF
        WRITE(cmessage,FMT=TRIM(ADJUSTL(outformat)))                           &
           'RXN NOT FOUND', cdhetero,                                          &
           (TRIM(ADJUSTL(asad_chemdiags(l)%reactants(j))),                     &
           j=1,2),'->',                                                        &
           (TRIM(ADJUSTL(asad_chemdiags(l)%products(j))),                      &
           j=1,asad_chemdiags(l)%num_products)
        errcode=asad_chemdiags(l)%stash_number
        CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
      END IF
    CASE (cdtermol) ! termolecular reactions/rates
      ALLOCATE(findrxn_tmp(1:(jptk+1)))
      findrxn_tmp = cd_findreaction(                                           &
           asad_chemdiags(l)%num_products,                                     &
           asad_chemdiags(l)%reactants,                                        &
           asad_chemdiags(l)%products,                                         &
           spt, ntrkx, (jptk+1), jpspt )
      asad_chemdiags(l)%location =                                             &
           findrxn_tmp(asad_chemdiags(l)%find_rxn_loc)
      DEALLOCATE(findrxn_tmp)
      IF (asad_chemdiags(l)%location <= 0) THEN
        IF (asad_chemdiags(l)%num_products < 10) THEN
          WRITE(outformat,FMT='(A28,I1,A7)')                                   &
               '(A13,1X,A1,1X,2(A,1X),A2,1X,',                                 &
               asad_chemdiags(l)%num_products,                                 &
               '(A,1X))'
        ELSE
          WRITE(outformat,FMT='(A28,I2,A7)')                                   &
               '(A13,1X,A1,1X,2(A,1X),A2,1X,',                                 &
               asad_chemdiags(l)%num_products,                                 &
               '(A,1X))'
        END IF
        WRITE(cmessage,FMT=TRIM(ADJUSTL(outformat)))                           &
           'RXN NOT FOUND', cdtermol,                                          &
           (TRIM(ADJUSTL(asad_chemdiags(l)%reactants(j))),                     &
           j=1,2),'->',                                                        &
           (TRIM(ADJUSTL(asad_chemdiags(l)%products(j))),                      &
           j=1,asad_chemdiags(l)%num_products)
        errcode=asad_chemdiags(l)%stash_number
        WRITE(umMessage,'(A,A13,I5,A5,I6)')                                    &
           cmessage, ' Error code: ',errcode, ' PE: ',mype
        CALL umPrint(umMessage,src=RoutineName)
        CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
      END IF
    CASE (cdphot) ! photolysis reactions/rates
      ALLOCATE(findrxn_tmp(1:(jppj+1)))
      findrxn_tmp = cd_findreaction(                                           &
           asad_chemdiags(l)%num_products,                                     &
           asad_chemdiags(l)%reactants,                                        &
           asad_chemdiags(l)%products,                                         &
           spj, nprkx, (jppj+1), jpspj )
      asad_chemdiags(l)%location =                                             &
           findrxn_tmp(asad_chemdiags(l)%find_rxn_loc)
      DEALLOCATE(findrxn_tmp)
      IF (asad_chemdiags(l)%location <= 0) THEN
        IF (asad_chemdiags(l)%num_products < 10) THEN
          WRITE(outformat,FMT='(A28,I1,A7)')                                   &
               '(A13,1X,A1,1X,2(A,1X),A2,1X,',                                 &
               asad_chemdiags(l)%num_products,                                 &
               '(A,1X))'
        ELSE
          WRITE(outformat,FMT='(A28,I2,A7)')                                   &
               '(A13,1X,A1,1X,2(A,1X),A2,1X,',                                 &
               asad_chemdiags(l)%num_products,                                 &
               '(A,1X))'
        END IF
        WRITE(cmessage,FMT=TRIM(ADJUSTL(outformat)))                           &
           'RXN NOT FOUND', cdphot,                                            &
           (TRIM(ADJUSTL(asad_chemdiags(l)%reactants(j))),                     &
           j=1,2),'->',                                                        &
           (TRIM(ADJUSTL(asad_chemdiags(l)%products(j))),                      &
           j=1,asad_chemdiags(l)%num_products)
        errcode=asad_chemdiags(l)%stash_number
        WRITE(umMessage,'(A,A13,I5,A5,I6)')                                    &
           cmessage, ' Error code: ',errcode, ' PE: ',mype
        CALL umPrint(umMessage,src=RoutineName)
        CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
      END IF
    CASE DEFAULT
      cmessage='REACTION TYPE '//                                              &
           asad_chemdiags(l)%rxn_type//' NOT FOUND'
      errcode=asad_chemdiags(l)%stash_number
      WRITE(umMessage,'(A,A13,I5,A5,I6)')                                      &
         cmessage, ' Error code: ',errcode, ' PE: ',mype
      CALL umPrint(umMessage,src=RoutineName)
      CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
    END SELECT       ! asad_chemdiags(l)%rxn_type
  CASE (cddep)
    SELECT CASE (asad_chemdiags(l)%rxn_type)
    CASE (cddry) ! DRY DEP
      ! %species contains species to be dry deposited
      ! %location contains location in requisite array
      idep=1
      DO j=1,jpspec
        IF (ldepd(j) .AND. (asad_chemdiags(l)%species ==                       &
             speci(j))) THEN
            ! it is being deposited logically and in diagnostic sense
          asad_chemdiags(l)%location = j
          idep=idep+1
        ELSE IF (.NOT. ldepd(j) .AND.                                          &
           (asad_chemdiags(l)%species == speci(j))) THEN
          ! diagnostics of deposition requested, but species is not
          ! labelled as being dry deposited in the chch_defs file
          cmessage='SPECIES '//asad_chemdiags(l)%species                       &
                 //' NOT DRY DEPOSITED'
          errcode=asad_chemdiags(l)%stash_number
          WRITE(umMessage,'(A,A13,I5,A5,I6)')                                  &
             cmessage, ' Error code: ', errcode,' PE: ',mype
          CALL umPrint(umMessage,src=RoutineName)
          CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
        END IF
      END DO
    CASE (cdwet) ! WET DEP
      DO j=1,jpspec
        IF (ldepw(j) .AND. (asad_chemdiags(l)%species ==                       &
             speci(j))) THEN
          asad_chemdiags(l)%location = j
          idep=idep+1
        ELSE IF (.NOT. ldepw(j) .AND.                                          &
           (asad_chemdiags(l)%species == speci(j))) THEN
           ! diagnostics of deposition requested, but species is not
          ! labelled as being wet deposited in the chch_defs file
          cmessage='SPECIES '//asad_chemdiags(l)%species                       &
                  //' NOT WET DEPOSITED'
          errcode=asad_chemdiags(l)%stash_number
          WRITE(umMessage,'(A,A13,I5,A5,I6)')                                  &
             cmessage, ' Error code: ', errcode,' PE: ',mype
          CALL umPrint(umMessage,src=RoutineName)
          CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
        END IF
      END DO
      ! %species contains species to be wet deposited
      ! %location contains location in requisite array
    END SELECT    ! asad_chemdiags(l)%rxn_type
  END SELECT    ! asad_chemdiags(l)%diag_type


  ! assumes 2 reactants
  IF (PrintStatus == Prstatus_diag) THEN
    SELECT CASE (asad_chemdiags(l)%diag_type)
    CASE (cdrxn,cdrte)
      WRITE(umMessage,'(A5,I6,A10,I5,A1,A1,A3,A10,A3,A10,A7,I6)')              &
           'mype=',mype,' CD:FIND: ',                                          &
           asad_chemdiags(l)%stash_number,' ',                                 &
           asad_chemdiags(l)%rxn_type,' : ',                                   &
           asad_chemdiags(l)%reactants(1),                                     &
           ' + ',asad_chemdiags(l)%reactants(2),                               &
           ' rxn = ',                                                          &
           asad_chemdiags(l)%location
      CALL umPrint(umMessage,src=RoutineName)
    CASE (cddep)
      WRITE(umMessage,'(A5,I6,A10,I5,A1,A1,A3,A10,A13,I6)')                    &
           'mype=',mype,' CD:FIND: ',                                          &
           asad_chemdiags(l)%stash_number,' ',                                 &
           asad_chemdiags(l)%rxn_type,' : ',                                   &
           asad_chemdiags(l)%species,                                          &
           ' deposition: ',                                                    &
           asad_chemdiags(l)%location
      CALL umPrint(umMessage,src=RoutineName)
    END SELECT
  END IF           ! Printstatus
END DO ! i

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

CONTAINS

! #####################################################################

FUNCTION cd_findreaction( numprods, reactants, products,                       &
                         reactions, r_index, nr, ns )

!  cd_findreaction   - does a search for a specified reaction in the given
!                   ratefile. Note that this function is clever enough to
!                   look for all the various combinations in which the
!                   species might appear.

!  Glenn Carver, Centre for Atmospheric Science, University of Cambridge
!
!  ASAD: cd_findreaction            Version: cd_findreaction.f 1.2 12/21/01

!  Purpose
!  -------
!  Given a reaction, this function looks for a matching entry in the ratefile.
!  The reactants and products given to this routine are not order dependent.
!  This function will look for the reactants in any order, ie. react1 / react2
!  or react2 / react1. And likewise for the products. If a reaction only has
!  two products, then specify the third product as a blank string.
!
!  Interface
!  ---------
!  react1 & react2   - character strings holding the two reactants of the
!                      reaction to look for.
!  prod1, prod2 & prod3 - character strings holding the 3 products of the
!                      reaction to look for. If there are less than 3 products
!                      then give a blank string.
!  reactions         - character array holding the array to search through.
!                      Dimensioned as (Nreacts,Entries). Nreacts (nr)is the no.
!                      of reactions. Entries (ns) will be the total no. of
!                      reactants and products (e.g. 2 reactants + 3 products =5)
!  index             - this is an integer array dimensioned as (nreacts). ASAD
!                      reorders the reactions from the external ratefiles for
!                      computational efficiency so the order they are stored
!                      is not the same order they appear in the ratefile. Pass
!                      the appropriate indexing array in order that the correct
!                      position of this reaction in the ascii ratefile can be
!                      returned. e.g. use one of the arrays, nbrkx, ntrkx etc
!
!  NOTE!!  Because of the way the string matching is done in this routine, you
!  **MUST** pass the strings with trailing spaces. e.g. if the string length
!  allowed for species names is 6, then the argument list must look like:
!
!  i = cd_findreaction( 'OH    ', 'CO2   ', 'HOCO2 ', '      ', '      ', .... )
!
!  If this function finds a matching reaction, it will return the index
!  to that reaction in the array reactions. If it doesn't find a match it
!  returns 0.
!
!----------------------------------------------------------------------
!
IMPLICIT NONE

! number of products in product array
INTEGER, INTENT(IN)           :: numprods
INTEGER, INTENT(IN)           :: nr
INTEGER, INTENT(IN)           :: ns
CHARACTER(LEN=10), INTENT(IN) :: reactants(number_of_reactants)
CHARACTER(LEN=10), INTENT(IN) :: products(numprods)
! array containing reactions
CHARACTER(LEN=10), INTENT(IN) :: reactions(nr,ns)
INTEGER, INTENT(IN)           :: r_index(nr)

!Local
CHARACTER(LEN=10), ALLOCATABLE :: prods(:)     ! found products
CHARACTER(LEN=10), ALLOCATABLE :: inprods(:)   ! input products
LOGICAL, ALLOCATABLE           :: found(:)     ! T if prods==inprods
LOGICAL                        :: Match        ! T if all prods match inprods

INTEGER :: j, jp, jp2, jr
INTEGER :: nreacts, nprods, iprods, ipargs, jprod

INTEGER :: cd_findreaction(1:nr)
INTEGER :: find_count

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CD_FINDREACTION'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

cd_findreaction(:) = -1 ! set to -1 to see if haven't got it at end
find_count = 0

!     -----------------------------------------------------------------
!          Search for reactants first.
!          ------ --- --------- ------
!

nreacts = SIZE(reactions,DIM=1)
nprods  = SIZE(reactions,DIM=2) - number_of_reactants ! always 2 reactants

IF (nprods < numprods) THEN ! if number of products is not equal to number of
                            ! number of products in reactions array then exit
  cmessage='wrong number of products in find reaction'
  errcode=nprods
  WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',errcode,' PE: ',mype
  CALL umPrint(umMessage,src=RoutineName)
  CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
END IF

ALLOCATE(prods(nprods), inprods(nprods), found(nprods) )

!  Get the non-blank products that we're looking for.
prods(:) = BlankStr
inprods(:) = BlankStr
ipargs = 0
inprods = BlankStr
DO jprod=1,numprods ! cycle over number of products
  IF ( products(jprod) /= BlankStr ) THEN
    ipargs = ipargs + 1
    inprods(ipargs) = products(jprod)
  END IF
END DO

j = 1
find_loop: DO
  jr = r_index(j)            ! loop over reactions

  ! ASSUME 2 REACTANTS HERE!!
  IF ((reactions(j,1)==reactants(1) .AND.                                      &
       reactions(j,2)==reactants(2)) .OR.                                      &
       (reactions(j,1)==reactants(2) .AND.                                     &
       reactions(j,2)==reactants(1))) THEN

    !          Found a reaction with matching reactants, now check all the
    !          possible cases that the products could match. Need to allow for
    !          rate files which have varying no. of products.

    prods = BlankStr
    found = .FALSE.
    iprods = 0
    DO jp = 1, nprods   ! only copy the non-blank products
      IF ( reactions(j,2+jp) /= BlankStr ) THEN
        iprods = iprods + 1
        prods(iprods) = reactions(j,2+jp)
      END IF
    END DO

    ! If no. of nonblank products doesn't match then try next reaction.
    IF ( iprods /= ipargs ) THEN
      j = j + 1
      IF ( j > nreacts ) THEN
        IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,        &
                                zhook_handle)
        RETURN
      END IF
      CYCLE
    END IF

    ! Otherwise check the names of the products

    DO jp = 1, iprods
      DO jp2 = 1, iprods
        IF (.NOT. found(jp) .AND. inprods(jp)==prods(jp2))                     &
          THEN
          found(jp) = .TRUE.
          prods(jp2) = BlankStr
        END IF
      END DO
    END DO

    ! Check to see if we have found all the products. If we have then exit

    Match = .TRUE.
    DO jp = 1, iprods
      Match = Match .AND. found(jp)
    END DO

    IF ( Match ) THEN
      find_count = find_count+1
      cd_findreaction(find_count) = jr
      ! don't return from here since need to see if there is
      ! another reaction with the same reactants & products
    END IF
  END IF

  ! next reaction

  j = j + 1
  IF ( j > nreacts ) THEN
    EXIT find_loop
  END IF
END DO find_loop

IF (ALLOCATED(prods)) DEALLOCATE( prods )
IF (ALLOCATED(inprods)) DEALLOCATE( inprods )
IF (ALLOCATED(found)) DEALLOCATE(found)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION cd_findreaction

END SUBROUTINE asad_chemical_diagnostics_init

! #####################################################################
SUBROUTINE asad_chemical_diagnostics(row_length, rows, model_levels,           &
     chunk_size, dpd_full, dpw_full, prk_full, y_full,                         &
     ix, jy, klevel, volume, ierr)

USE chemistry_constants_mod,  ONLY: avogadro
USE ukca_config_specification_mod, ONLY: ukca_config

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length       ! length of row
INTEGER, INTENT(IN) :: rows             ! number of rows
INTEGER, INTENT(IN) :: model_levels     ! the level that we are at
INTEGER, INTENT(IN) :: chunk_size       ! size of input arrays
INTEGER, INTENT(IN) :: ix               ! i counter
INTEGER, INTENT(IN) :: jy               ! j counter
INTEGER, INTENT(IN) :: klevel           ! the level that we are at
REAL, INTENT(IN)    :: volume(row_length, rows, model_levels)   ! cell volume
INTEGER, INTENT(IN OUT) :: ierr          ! error code
REAL, INTENT(IN)    :: dpd_full(chunk_size,jpspec)
REAL, INTENT(IN)    :: dpw_full(chunk_size,jpspec)
REAL, INTENT(IN)    :: prk_full(chunk_size,jpnr)
REAL, INTENT(IN)    :: y_full(chunk_size,jpspec)

REAL, PARAMETER :: convfac = 1.0e6/avogadro ! conversion factor to convert
                                            ! volume into cm^3 and result in mol

LOGICAL, SAVE :: firstcall=.TRUE.
LOGICAL, SAVE :: first_pass=.TRUE.

INTEGER :: k                      ! index for vertical coordinate
INTEGER :: l                      ! counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER                           :: errcode

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_CHEMICAL_DIAGNOSTICS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ukca_config%l_ukca_asad_full) THEN
  IF (chunk_size /= row_length * rows * model_levels) THEN
    cmessage = 'Must pass full domain when l_ukca_asad_full is set to .true.'
    errcode = chunk_size
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
ELSE IF (ukca_config%l_ukca_asad_columns) THEN
  IF (chunk_size /= model_levels) THEN
    cmessage = 'Must pass column when l_ukca_asad_columns is set to .true.'
    errcode = chunk_size
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
END IF

ierr = -1
! OMP CRITICAL will only allow one thread through this code at a time,
! while the other threads are held until completion.
!$OMP CRITICAL (chemical_diagnostics_init)
IF (first_pass) THEN
  IF (firstcall) THEN
    CALL asad_chemical_diagnostics_init(ierr)
    firstcall=.FALSE.
  END IF   ! firstcall
  first_pass=.FALSE.
END IF     ! first_pass
!$OMP END CRITICAL (chemical_diagnostics_init)

IF (ukca_config%l_ukca_asad_full) THEN
  ! Go through and pick up fluxes from ASAD arrays
  ! prk is in units of molecules.cm^-3.s^-1
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(l)                                        &
!$OMP SHARED(asad_chemdiags, dpd_full, dpw_full, L_stratosphere,               &
!$OMP        model_levels, n_chemdiags, prk_full, row_length, rows, volume,    &
!$OMP        y_full)
!$OMP DO SCHEDULE(DYNAMIC)
  DO l=1,n_chemdiags
    ! in this case ASAD is being called over the full domain so we
    ! need to fill %throughput by reshaping the array
    SELECT CASE (asad_chemdiags(l)%diag_type)
    CASE (cdrxn)
      asad_chemdiags(l)%throughput(:,:,:) =                                    &
           RESHAPE(prk_full(:,asad_chemdiags(l)%location),                     &
           [row_length,rows,model_levels])*volume(:,:,:)*convfac
      IF (asad_chemdiags(l)%tropospheric_mask) THEN
        WHERE (L_stratosphere(:,:,:))
          asad_chemdiags(l)%throughput(:,:,:) = 0.0
        END WHERE
      END IF
    CASE (cddep)
      SELECT CASE (asad_chemdiags(l)%rxn_type)
      CASE (cddry) ! DRY DEP
        asad_chemdiags(l)%throughput(:,:,:)=                                   &
             RESHAPE(dpd_full(:,asad_chemdiags(l)%location),                   &
                     [row_length,rows,model_levels])*                          &
             RESHAPE(y_full(:,asad_chemdiags(l)%location),                     &
                     [row_length,rows,model_levels])*                          &
             volume(:,:,:)*convfac
        IF (asad_chemdiags(l)%tropospheric_mask) THEN
          WHERE (L_stratosphere(:,:,:))
            asad_chemdiags(l)%throughput(:,:,:) = 0.0
          END WHERE
        END IF
        ! Not needed (?)
      CASE (cdwet) ! WET DEP
        asad_chemdiags(l)%throughput(:,:,:)=                                   &
             RESHAPE(dpw_full(:,asad_chemdiags(l)%location),                   &
                     [row_length,rows,model_levels])*                          &
             RESHAPE(y_full(:,asad_chemdiags(l)%location),                     &
                     [row_length,rows,model_levels])*                          &
             volume(:,:,:)*convfac
        ! Not needed (?)
        IF (asad_chemdiags(l)%tropospheric_mask) THEN
          WHERE (L_stratosphere(:,:,:))
            asad_chemdiags(l)%throughput(:,:,:) = 0.0
          END WHERE
        END IF
      END SELECT
    END SELECT
  END DO ! l=1,n_chemdiags
!$OMP END DO
!$OMP END PARALLEL
ELSE IF (ukca_config%l_ukca_asad_columns) THEN
  ! Go through and pick up fluxes from ASAD arrays
  ! prk is in units of molecules.cm^-3.s^-1
  DO l=1,n_chemdiags
    ! in this case ASAD is being called column-by column so we need
    ! to fill %throughput by iterating over X & Y
    SELECT CASE (asad_chemdiags(l)%diag_type)
    CASE (cdrxn)
      asad_chemdiags(l)%throughput(ix,jy,:) =                                  &
           prk_full(:,asad_chemdiags(l)%location)*volume(ix,jy,:)*convfac
      IF (asad_chemdiags(l)%tropospheric_mask) THEN
        WHERE (L_stratosphere(ix,jy,:))
          asad_chemdiags(l)%throughput(ix,jy,:) = 0.0
        END WHERE
      END IF
    CASE (cddep)
      SELECT CASE (asad_chemdiags(l)%rxn_type)
      CASE (cddry) ! DRY DEP
         ! dpd array is over a column, no need to consider klevel
         ! NOTE: Current STASH settings have dry-deposition output
         !       as 3D. If this changes the code below will no
         !       longer work as the code expects an entire
         !       atmospheric column.
        asad_chemdiags(l)%throughput(ix,jy,:)=                                 &
             dpd_full(:,asad_chemdiags(l)%location)*                           &
             y_full(:,asad_chemdiags(l)%location)*                             &
             volume(ix,jy,:)*convfac
        ! Not needed (probably) but kept for consistency
        IF (asad_chemdiags(l)%tropospheric_mask) THEN
          WHERE (L_stratosphere(ix,jy,:))
            asad_chemdiags(l)%throughput(ix,jy,:) = 0.0
          END WHERE
        END IF
      CASE (cdwet) ! WET DEP
        asad_chemdiags(l)%throughput(ix,jy,:)=                                 &
             dpw_full(:,asad_chemdiags(l)%location)*                           &
             y_full(:,asad_chemdiags(l)%location)*                             &
             volume(ix,jy,:)*convfac
        ! Not needed (probably) but kept for consistency
        IF (asad_chemdiags(l)%tropospheric_mask) THEN
          WHERE (L_stratosphere(ix,jy,:))
            asad_chemdiags(l)%throughput(ix,jy,:) = 0.0
          END WHERE
        END IF
      END SELECT
    END SELECT
  END DO ! l=1,n_chemdiags
ELSE
  ! Go through and pick up fluxes from ASAD arrays
  ! prk is in units of molecules.cm^-3.s^-1
  DO l=1,n_chemdiags
    ! in this case ASAD is being called horizontally so we need
    ! to fill %throughput by iterating over Z
    SELECT CASE (asad_chemdiags(l)%diag_type)
    CASE (cdrxn)
      asad_chemdiags(l)%throughput(:,:,klevel) =                               &
           RESHAPE(prk(:,asad_chemdiags(l)%location),                          &
           [row_length,rows])*volume(:,:,klevel)*convfac
      IF (asad_chemdiags(l)%tropospheric_mask) THEN
        WHERE (L_stratosphere(:,:,klevel))
          asad_chemdiags(l)%throughput(:,:,klevel) = 0.0
        END WHERE
      END IF
    CASE (cddep)
      SELECT CASE (asad_chemdiags(l)%rxn_type)
      CASE (cddry) ! DRY DEP
        IF (klevel <=                                                          &
             SIZE(asad_chemdiags(l)%throughput(:,:,:),DIM=3)) THEN
           ! if 2D then only take lowest level, otherwise will be 3D
          asad_chemdiags(l)%throughput(:,:,klevel)=                            &
               RESHAPE(dpd(:,asad_chemdiags(l)%location),                      &
                       [row_length,rows])*                                     &
               RESHAPE(y(:,asad_chemdiags(l)%location),                        &
                       [row_length,rows])*                                     &
               volume(:,:,klevel)*convfac
          IF (asad_chemdiags(l)%tropospheric_mask) THEN
            WHERE (L_stratosphere(:,:,klevel))
              asad_chemdiags(l)%throughput(:,:,klevel) = 0.0
            END WHERE
          END IF
        END IF
        ! Not needed (?)
      CASE (cdwet) ! WET DEP
        asad_chemdiags(l)%throughput(:,:,klevel)=                              &
             RESHAPE(dpw(:,asad_chemdiags(l)%location),                        &
                     [row_length,rows])*                                       &
             RESHAPE(y(:,asad_chemdiags(l)%location),                          &
                     [row_length,rows])*                                       &
             volume(:,:,klevel)*convfac
        ! Not needed (?)
        IF (asad_chemdiags(l)%tropospheric_mask) THEN
          WHERE (L_stratosphere(:,:,klevel))
            asad_chemdiags(l)%throughput(:,:,klevel) = 0.0
          END WHERE
        END IF
      END SELECT
    END SELECT
  END DO ! l=1,n_chemdiags
END IF ! l_ukca_asad_columns / l_ukca_asad_full


ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE asad_chemical_diagnostics

! #####################################################################
SUBROUTINE asad_emissions_diagnostics(row_length, rows,                        &
                                 n_chem_tracers,                               &
                                 em_field, surf_area,                          &
                                 n_emissions,                                  &
                                 n_boundary_vals, em_spec,                     &
                                 lbc_spec, molmass,                            &
                                 lbc_molmass, ierr)

USE ukca_constants,    ONLY: c_no, c_no2
IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length                        ! array dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: n_chem_tracers                    ! number of tracers
REAL, INTENT(IN)    :: em_field(row_length,rows,n_chem_tracers)
REAL, INTENT(IN)    :: surf_area(row_length,rows)        ! cell area
INTEGER, INTENT(IN) :: n_emissions                       ! no of emitted species
INTEGER, INTENT(IN) :: n_boundary_vals
CHARACTER(LEN=10), INTENT(IN) :: em_spec(n_emissions)
CHARACTER(LEN=10), INTENT(IN) :: lbc_spec(n_boundary_vals)
REAL, INTENT(IN) :: molmass(n_emissions)
REAL, INTENT(IN) :: lbc_molmass(n_boundary_vals)

INTEGER, INTENT(OUT) :: ierr                               ! error code

INTEGER :: l,ntracer,mspecies                          ! loop variables

LOGICAL :: L_cd_emitted(jpctr)                 ! T for emission diags

CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER       :: errcode

LOGICAL, SAVE :: firstcall=.TRUE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_EMISSIONS_DIAGNOSTICS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr = -1

! Set up logicals and molmass for use in later calculations
IF (firstcall) THEN

  L_cd_emitted(:) =.FALSE.
  DO l=1,n_chemdiags
    IF (asad_chemdiags(l)%diag_type == cdems) THEN
      DO ntracer=1,jpctr
        DO mspecies=1,n_emissions
           ! work out which molmass to take!
          IF ((advt(ntracer) == em_spec(mspecies)) .AND.                       &
             (asad_chemdiags(l)%species==em_spec(mspecies))) THEN
            asad_chemdiags(l)%molecular_mass = molmass(mspecies)
            asad_chemdiags(l)%c_vmr_to_mmr = c_species(ntracer)
            L_cd_emitted(ntracer) = .TRUE.
            IF (asad_chemdiags(l)%species == 'NO        ')                     &
                       asad_chemdiags(l)%molecular_mass =                      &
                       asad_chemdiags(l)%molecular_mass*c_no/c_no2
          END IF
        END DO        ! mspecies,n_emissions
        DO mspecies=1,n_boundary_vals
           ! work out which molmass to take!, can overwrite with
           ! lbc_molmass, since emmissions_ctl does this anyway
          IF ((advt(ntracer) == lbc_spec(mspecies)) .AND.                      &
              (asad_chemdiags(l)%species == lbc_spec(mspecies)))               &
                     THEN
            asad_chemdiags(l)%molecular_mass=lbc_molmass(mspecies)
            asad_chemdiags(l)%c_vmr_to_mmr = c_species(ntracer)
            L_cd_emitted(ntracer) = .TRUE.
            IF (asad_chemdiags(l)%species == 'NO        ')                     &
                 asad_chemdiags(l)%molecular_mass =                            &
                 asad_chemdiags(l)%molecular_mass*c_no/c_no2
          END IF
        END DO        ! mspecies,n_boundary_vals

        IF (L_cd_emitted(ntracer)) THEN
          IF (advt(ntracer) == asad_chemdiags(l)%species)                      &
                    asad_chemdiags(l)%location = ntracer
        ELSE IF ((.NOT. L_cd_emitted(ntracer)) .AND.                           &
           (advt(ntracer) == asad_chemdiags(l)%species)) THEN
          cmessage='SPECIES '// asad_chemdiags(l)%species                      &
                //' NOT EMITTED'
          errcode=asad_chemdiags(l)%stash_number
          WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',           &
                errcode,' PE: ',mype
          CALL umPrint(umMessage,src=RoutineName)
          CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
        END IF
      END DO           ! ntracer,jpctr
    END IF              ! %diag_type == cdems
  END DO                 ! l,n_chemdiags
  firstcall=.FALSE.
END IF                    ! first

DO l=1,n_chemdiags
  IF (asad_chemdiags(l)%diag_type == cdems) THEN
    ! emissions are 2D and converted to (kg/m2/s)
    asad_chemdiags(l)%throughput(:,:,1) =                                      &
        em_field(:,:,asad_chemdiags(l)%location)*surf_area(:,:)*               &
        (1000.0/asad_chemdiags(l)%molecular_mass)
  END IF                 ! %diag_type
END DO                    ! l,n_chemdiags


ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_emissions_diagnostics

! #####################################################################
SUBROUTINE asad_3d_emissions_diagnostics( row_length, rows,                    &
     model_levels, i_spec_emiss, em_field_3d_in, surf_area,                    &
     mass, em_molmass, timestep,                                               &
     emission_type, ierr)

IMPLICIT NONE


INTEGER, INTENT(IN) :: row_length                      ! Array dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: i_spec_emiss
REAL, INTENT(IN) :: em_field_3d_in(1:row_length,1:rows,                        &
                                   1:model_levels)
REAL, INTENT(IN) :: surf_area(1:row_length,1:rows)
REAL, INTENT(IN) :: mass(1:row_length,1:rows,1:model_levels)
REAL, INTENT(IN) :: em_molmass
REAL, INTENT(IN) :: timestep
CHARACTER(LEN=1), INTENT(IN) :: emission_type
INTEGER, INTENT(OUT) :: ierr

REAL :: em_field_3d(1:row_length,1:rows,1:model_levels)


INTEGER :: i,j,l,klevel

CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER       :: errcode

LOGICAL, SAVE :: firstcall=.TRUE.
LOGICAL       :: spec_emitted

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_3D_EMISSIONS_DIAGNOSTICS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr = -1

IF (firstcall) THEN
  spec_emitted = .FALSE.
  DO l=1,n_chemdiags
    IF (asad_chemdiags(l)%diag_type == cdems .AND.                             &
        asad_chemdiags(l)%rxn_type == emission_type .AND.                      &
        asad_chemdiags(l)%species == advt(i_spec_emiss)) THEN
      asad_chemdiags(l)%location = i_spec_emiss
      asad_chemdiags(l)%molecular_mass = em_molmass
      asad_chemdiags(l)%c_vmr_to_mmr   = c_species(i_spec_emiss)
      spec_emitted = .TRUE.
    END IF         ! diag_type == cdems &etc
  END DO
  IF (.NOT. spec_emitted) THEN
    cmessage='SPECIES '//asad_chemdiags(l)%species                             &
         //' NOT EMITTED IN 3D'
    errcode=asad_chemdiags(l)%stash_number
    WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',errcode,         &
         ' PE: ',mype
    CALL umPrint(umMessage,src=RoutineName)
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
  firstcall = .FALSE.
END IF                    ! firstcall

! need to do some conversion because of the way that things are done in
! UKCA for 3D emissions
em_field_3d(:,:,:) = -1.0 ! i.e. if values are -ve, have a problem
! sanity check, just in case
!em_field_3d(:,:,:) = em_field_3d_in(:,:,:)
DO l=1,n_chemdiags
  SELECT CASE (asad_chemdiags(l)%diag_type)
  CASE (cdems)
    IF ((asad_chemdiags(l)%rxn_type == aircraft_emissions)                     &
         .AND. (emission_type == aircraft_emissions)) THEN
       ! these are in kg(NO2)/m2/s
       ! - need to multiply by the surface area and / air mass
      DO klevel=1,model_levels
        DO j=1,rows
          DO i=1,row_length
            em_field_3d(i,j,klevel) =                                          &
             (em_field_3d_in(i,j,klevel)/mass(i,j,klevel))*                    &
              surf_area(i,j)
          END DO ! i
        END DO ! j
      END DO ! klevel
    ELSE IF ((asad_chemdiags(l)%rxn_type ==                                    &
         lightning_emissions) .AND. (emission_type ==                          &
         lightning_emissions)) THEN
       ! do nothing - in kg(NO)/kg(air)/gridcell/s
      em_field_3d(:,:,:) = em_field_3d_in(:,:,:)
    ELSE IF ((asad_chemdiags(l)%rxn_type ==                                    &
         volcanic_emissions) .AND. (emission_type ==                           &
         volcanic_emissions)) THEN
       ! volcanic emissions are kg(SO2)/kg(air)/gridcell/timestep
       ! - so divide by timestep to get /s
      em_field_3d(:,:,:) = em_field_3d_in(:,:,:)/timestep
    ELSE
       ! do nothing - assume kg(Species)/kg(air)/gridcell/s
      em_field_3d(:,:,:) = em_field_3d_in(:,:,:)
    END IF
  END SELECT
END DO


! mass calculation method
! converts from kg(species)/kg(air)/s to mol/s
DO l=1,n_chemdiags
  IF (asad_chemdiags(l)%diag_type == cdems) THEN
    IF ((advt(i_spec_emiss) == asad_chemdiags(l)%species) .AND.                &
             (emission_type == asad_chemdiags(l)%rxn_type)) THEN
      asad_chemdiags(l)%throughput(:,:,:)= (em_field_3d(                       &
                     1:row_length,1:rows,1:model_levels)                       &
                     *mass(1:row_length,1:rows,1:model_levels))                &
                      *(1000.0/asad_chemdiags(l)%molecular_mass)
      IF (asad_chemdiags(l)%tropospheric_mask) THEN
        WHERE (L_stratosphere)
          asad_chemdiags(l)%throughput(:,:,:) = 0.0
        END WHERE
      END IF
    END IF              ! advt==%species
  END IF             ! %diag_type == cdems
END DO                    ! l,n_chemdiags

ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_3d_emissions_diagnostics

! #####################################################################
SUBROUTINE asad_tropospheric_mask(ierr)

USE ukca_tropopause , ONLY: L_stratosphere
IMPLICIT NONE

INTEGER, INTENT(OUT) :: ierr

INTEGER :: l

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_TROPOSPHERIC_MASK'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr = 1

IF (L_asad_use_trop_mask_output) THEN
  DO l=1,n_chemdiags
    SELECT CASE (asad_chemdiags(l)%diag_type)
    CASE (cdtpm)
      asad_chemdiags(l)%throughput(:,:,:) = 1.0e0
      WHERE (L_stratosphere(:,:,:))
        asad_chemdiags(l)%throughput(:,:,:) = 0.0e0
      END WHERE
    END SELECT ! %diag_type
  END DO ! l
END IF

ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_tropospheric_mask

! #####################################################################
SUBROUTINE asad_tendency_ste(row_length, rows, model_levels,                   &
     n_chem_tracers,all_tracers, total_number_density, volume,                 &
     timestep, which_diagnostic, L_store_value, ierr)

USE chemistry_constants_mod,   ONLY: avogadro
USE ukca_constants,   ONLY: m_air
IMPLICIT NONE


INTEGER, INTENT(IN)  :: row_length     ! length of row
INTEGER, INTENT(IN)  :: rows           ! number of rows
INTEGER, INTENT(IN)  :: model_levels   ! number of levels
INTEGER, INTENT(IN)  :: n_chem_tracers ! number of chemical tracers

! tracer array
REAL, INTENT(IN)    :: all_tracers(row_length, rows, model_levels,             &
                                   n_chem_tracers)
! number density array
REAL, INTENT(IN)    :: total_number_density(row_length, rows,                  &
                                            model_levels)
! cell volume array
REAL, INTENT(IN)    :: volume(row_length,rows,model_levels)
REAL, INTENT(IN) :: timestep

CHARACTER(LEN=3), INTENT(IN) :: which_diagnostic

LOGICAL, INTENT(IN) :: L_store_value    ! T to store value, F to make difference

INTEGER, INTENT(OUT) :: ierr            ! error code

INTEGER :: l,m                    ! counters

LOGICAL, SAVE :: firstcall_tendency=.TRUE.
LOGICAL, SAVE :: firstcall_STE=.TRUE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_TENDENCY_STE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


ierr = 1

! set up location and molecular mass from c_species array
IF (firstcall_tendency .OR. firstcall_STE) THEN

  DO m=1,n_chem_tracers
    chemdiags_loop_first: DO l=1,n_chemdiags
      IF (which_diagnostic /= asad_chemdiags(l)%diag_type)                     &
              CYCLE chemdiags_loop_first
      SELECT CASE (asad_chemdiags(l)%diag_type)
      CASE (cdnet,cdste)
        IF (advt(m) == asad_chemdiags(l)%species) THEN
          asad_chemdiags(l)%location = m
          asad_chemdiags(l)%molecular_mass = c_species(m)*m_air
          asad_chemdiags(l)%c_vmr_to_mmr = c_species(m)
        END IF
      END SELECT
    END DO chemdiags_loop_first
  END DO

  SELECT CASE (which_diagnostic)
  CASE (cdnet)
    firstcall_tendency=.FALSE.
  CASE (cdste)
    firstcall_STE=.FALSE.
  END SELECT
END IF         ! firstcall_tendency .OR. firstcall_STE

! calculate tendency/STE
chemdiags_loop: DO l=1,n_chemdiags
  IF (which_diagnostic /= asad_chemdiags(l)%diag_type)                         &
        CYCLE chemdiags_loop
  SELECT CASE (asad_chemdiags(l)%diag_type)
  CASE (cdnet,cdste)
    IF (L_store_value) THEN         ! store value in %throughput
      ! Volume calculation method, converts to mol
      asad_chemdiags(l)%throughput(:,:,:) =                                    &
            all_tracers(:,:,:,asad_chemdiags(l)%location)                      &
            *total_number_density(1:row_length,                                &
            1:rows,1:model_levels)*                                            &
            volume(1:row_length,1:rows,1:model_levels)                         &
                   /(avogadro*asad_chemdiags(l)%c_vmr_to_mmr)
    ELSE
      ! Take difference from previously stored value
      asad_chemdiags(l)%throughput(:,:,:) =                                    &
          ((all_tracers(:,:,:,asad_chemdiags(l)%location)                      &
                   *total_number_density(1:row_length,1:rows,                  &
                   1:model_levels)*                                            &
                   volume(1:row_length,1:rows,1:model_levels)                  &
                   /(avogadro *asad_chemdiags(l)%c_vmr_to_mmr))                &
                   -asad_chemdiags(l)%throughput(:,:,:))                       &
                   /timestep
      ! Mask out stratospheric values if requested
      IF (asad_chemdiags(l)%tropospheric_mask) THEN
        WHERE (L_stratosphere)
          asad_chemdiags(l)%throughput(:,:,:) = 0.0
        END WHERE
      END IF
    END IF              ! L_store_value
  END SELECT             ! asad_chemdiags(l)%diag_type
END DO chemdiags_loop     ! l=1,n_chemdiags

ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_tendency_ste

! #####################################################################
SUBROUTINE asad_mass_diagnostic(row_length, rows, model_levels,                &
                                mass, ierr)
IMPLICIT NONE


INTEGER, INTENT(IN)  :: row_length                           ! array dimension
INTEGER, INTENT(IN)  :: rows                                 ! array dimension
INTEGER, INTENT(IN)  :: model_levels                         ! array dimension
REAL,    INTENT(IN)  :: mass(row_length,rows,model_levels)   ! mass of cell
INTEGER, INTENT(OUT) :: ierr                                 ! error code

INTEGER :: l                      ! counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_MASS_DIAGNOSTIC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr = -1

DO l=1,n_chemdiags
  IF (asad_chemdiags(l)%diag_type == cdmas) THEN
    asad_chemdiags(l)%throughput(:,:,:) = mass(:,:,:)
    IF (asad_chemdiags(l)%tropospheric_mask) THEN
      WHERE (L_stratosphere)
        asad_chemdiags(l)%throughput(:,:,:) = 0.0
      END WHERE
    END IF
  END IF
END DO                    ! l,n_chemdiags

ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_mass_diagnostic
! #####################################################################
SUBROUTINE asad_psc_diagnostic(row_length, rows, model_levels, chunk_size,     &
                               fpsc1_full, fpsc2_full, ix, jy, klevel, ierr)

USE ukca_config_specification_mod, ONLY: ukca_config
IMPLICIT NONE


INTEGER, INTENT(IN)  :: row_length             ! array dimension
INTEGER, INTENT(IN)  :: rows                   ! array dimension
INTEGER, INTENT(IN)  :: model_levels           ! number of model levels
INTEGER, INTENT(IN)  :: chunk_size             ! size of input arrays
INTEGER, INTENT(IN)  :: ix                     ! i counter
INTEGER, INTENT(IN)  :: jy                     ! j counter
INTEGER, INTENT(IN)  :: klevel                 ! level number
INTEGER, INTENT(OUT) :: ierr                   ! error code
REAL, INTENT(IN) :: fpsc1_full(chunk_size)
REAL, INTENT(IN) :: fpsc2_full(chunk_size)

INTEGER :: l                                   ! counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER                           :: errcode

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_PSC_DIAGNOSTIC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
ierr = -1

IF (ukca_config%l_ukca_asad_full) THEN
  IF (chunk_size /= row_length * rows * model_levels) THEN
    cmessage = 'Must pass full domain when l_ukca_asad_full is set to .true.'
    errcode = chunk_size
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
ELSE IF (ukca_config%l_ukca_asad_columns) THEN
  IF (chunk_size /= model_levels) THEN
    cmessage = 'Must pass column when l_ukca_asad_columns is set to .true.'
    errcode = chunk_size
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
END IF

IF (ukca_config%l_ukca_asad_full) THEN
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(l)                                        &
!$OMP SHARED(asad_chemdiags, fpsc1_full, fpsc2_full, L_stratosphere,           &
!$OMP        model_levels, n_chemdiags, row_length, rows)
!$OMP DO SCHEDULE(DYNAMIC)
  DO l=1,n_chemdiags
    ! in this case ASAD is being called over the full domain so we
    ! need to fill %throughput by reshaping the array
    SELECT CASE (asad_chemdiags(l)%diag_type)
    CASE (cdpsc)
      SELECT CASE(asad_chemdiags(l)%rxn_type)
      CASE (cdpsc_typ1)
        asad_chemdiags(l)%throughput(:,:,:) =                                  &
             RESHAPE(fpsc1_full(:),[row_length,rows,model_levels])
        IF (asad_chemdiags(l)%tropospheric_mask) THEN
          WHERE (L_stratosphere(:,:,:))
            asad_chemdiags(l)%throughput(:,:,:) = 0.0
          END WHERE
        END IF
      CASE (cdpsc_typ2)
        asad_chemdiags(l)%throughput(:,:,:) =                                  &
             RESHAPE(fpsc2_full(:),[row_length,rows,model_levels])
        IF (asad_chemdiags(l)%tropospheric_mask) THEN
          WHERE (L_stratosphere(:,:,:))
            asad_chemdiags(l)%throughput(:,:,:) = 0.0
          END WHERE
        END IF
      END SELECT
    END SELECT
  END DO
!$OMP END DO
!$OMP END PARALLEL
ELSE IF (ukca_config%l_ukca_asad_columns) THEN
  DO l=1,n_chemdiags
    ! in this case ASAD is being called column-by column so we need
    ! to fill %throughput by iterating over X & Y
    SELECT CASE (asad_chemdiags(l)%diag_type)
    CASE (cdpsc)
      SELECT CASE (asad_chemdiags(l)%rxn_type)
      CASE (cdpsc_typ1)
        asad_chemdiags(l)%throughput(ix,jy,:) = fpsc1_full(:)
        IF (asad_chemdiags(l)%tropospheric_mask) THEN
          WHERE (L_stratosphere(ix,jy,:))
            asad_chemdiags(l)%throughput(ix,jy,:) = 0.0
          END WHERE
        END IF
      CASE (cdpsc_typ2)
        asad_chemdiags(l)%throughput(ix,jy,:) = fpsc2_full(:)
        IF (asad_chemdiags(l)%tropospheric_mask) THEN
          WHERE (L_stratosphere(ix,jy,:))
            asad_chemdiags(l)%throughput(ix,jy,:) = 0.0
          END WHERE
        END IF
      END SELECT
    END SELECT
  END DO
ELSE
  DO l=1,n_chemdiags
    ! in this case ASAD is being called horizontally so we need
    ! to fill %throughput by iterating over Z
    SELECT CASE (asad_chemdiags(l)%diag_type)
    CASE (cdpsc)
      SELECT CASE (asad_chemdiags(l)%rxn_type)
      CASE (cdpsc_typ1)
        asad_chemdiags(l)%throughput(:,:,klevel) =                             &
             RESHAPE(fpsc1(:),[row_length,rows])
        IF (asad_chemdiags(l)%tropospheric_mask) THEN
          WHERE (L_stratosphere(:,:,klevel))
            asad_chemdiags(l)%throughput(:,:,klevel) = 0.0
          END WHERE
        END IF
      CASE (cdpsc_typ2)
        asad_chemdiags(l)%throughput(:,:,klevel) =                             &
             RESHAPE(fpsc2(:),[row_length,rows])
        IF (asad_chemdiags(l)%tropospheric_mask) THEN
          WHERE (L_stratosphere(:,:,klevel))
            asad_chemdiags(l)%throughput(:,:,klevel) = 0.0
          END WHERE
        END IF
      END SELECT
    END SELECT
  END DO
END IF ! l_ukca_asad_columns / l_ukca_asad_full

ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_psc_diagnostic

! #####################################################################
SUBROUTINE asad_output_tracer(row_length, rows, model_levels,                  &
                              n_chem_tracers, all_tracers, ierr)

USE ukca_constants, ONLY: m_air
IMPLICIT NONE


INTEGER, INTENT(IN)  :: row_length     ! length of row
INTEGER, INTENT(IN)  :: rows           ! number of rows
INTEGER, INTENT(IN)  :: model_levels   ! number of levels
INTEGER, INTENT(IN)  :: n_chem_tracers ! number of chemical tracers

! Tracer array
REAL, INTENT(IN)     :: all_tracers(row_length, rows, model_levels,            &
                                    n_chem_tracers)

INTEGER, INTENT(OUT) :: ierr           ! error code
INTEGER :: l,m                   ! counters

LOGICAL, SAVE :: firstcall=.TRUE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_OUTPUT_TRACER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr = 1

! set up location and molecular mass from c_species array
IF (firstcall) THEN

  DO m=1,n_chem_tracers
    DO l=1,n_chemdiags
      SELECT CASE (asad_chemdiags(l)%diag_type)
      CASE (cdout)
        IF (advt(m) == asad_chemdiags(l)%species) THEN
          asad_chemdiags(l)%location = m
          asad_chemdiags(l)%molecular_mass = c_species(m)*m_air
          asad_chemdiags(l)%c_vmr_to_mmr = c_species(m)
        END IF
      END SELECT
    END DO
  END DO

  firstcall=.FALSE.
END IF

! take copy of tracer, and apply tropospheric mask if requested
DO l=1,n_chemdiags
  SELECT CASE (asad_chemdiags(l)%diag_type)
  CASE (cdout)
    asad_chemdiags(l)%throughput(:,:,:) =                                      &
         all_tracers(:,:,:,asad_chemdiags(l)%location)
    IF (asad_chemdiags(l)%tropospheric_mask) THEN
      WHERE (L_stratosphere(:,:,:))
        asad_chemdiags(l)%throughput(:,:,:) = 0.0
      END WHERE
    END IF
  END SELECT             ! asad_chemdiags(l)%diag_type
END DO      ! l=1,n_chemdiags

ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_output_tracer

! #####################################################################
SUBROUTINE asad_lightning_diagnostics( row_length, rows,                       &
     which_diagnostic, diagnostic_field, ierr)

IMPLICIT NONE


INTEGER, INTENT(IN) :: row_length                   ! Field dimension
INTEGER, INTENT(IN) :: rows                         ! Field dimension
CHARACTER(LEN=1), INTENT(IN) :: which_diagnostic
REAL, INTENT(IN) :: diagnostic_field(1:row_length,1:rows)
INTEGER, INTENT(OUT) :: ierr

INTEGER :: l

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_LIGHTNING_DIAGNOSTICS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr = -1

DO l=1,n_chemdiags
  SELECT CASE (asad_chemdiags(l)%diag_type)
  CASE (cdlgt)
    ! only copy into diagnostic requested
    IF (asad_chemdiags(l)%rxn_type == which_diagnostic) THEN
      asad_chemdiags(l)%throughput(:,:,1) =                                    &
                 diagnostic_field(:,:)

      ! never used
      IF (asad_chemdiags(l)%tropospheric_mask) THEN
        WHERE (L_stratosphere(:,:,1))
          asad_chemdiags(l)%throughput(:,:,1) = 0.0
        END WHERE
      END IF
    END IF
  END SELECT
END DO

ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_lightning_diagnostics

END MODULE asad_chem_flux_diags
