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

USE ukca_um_legacy_mod, ONLY: copydiag, copydiag_3d,                           &
                              modl_b, isec_b, item_b, iopl_d,                  &
                              idom_b, ndiag, len_stlist, stindex, stlist,      &
                              num_stash_levels, stash_levels, si, sf, si_last, &
                              st_levels_model_theta, st_levels_single,         &
                              stashcode_ukca_chem_diag,                        &
                              submodel_for_sm, atmos_im, code, mype
USE asad_flux_dat,    ONLY: asad_chemical_fluxes
USE asad_mod,         ONLY: advt, dpd, dpw, fpsc1, fpsc2,                      &
                            jpspb, jpsph, jpspj, jpspt,                        &
                            ldepd, ldepw,                                      &
                            nbrkx, nhrkx, nprkx, ntrkx, prk,                   &
                            spb, speci, sph, spj, spt, y,                      &
                            jpspec, jpbk, jptk, jppj, jphk, jpnr, jpctr
USE ukca_tropopause,  ONLY: L_troposphere
USE ukca_cspecies,    ONLY: c_species
USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook
USE ereport_mod,      ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus,                         &
                      PrStatus_Oper, PrStatus_Diag, umPrintFlush
USE missing_data_mod,      ONLY: imdi
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

PRIVATE

! public types
PUBLIC :: chemdiag
! public variables
PUBLIC :: asad_chemdiags
PUBLIC :: chemD1codes
PUBLIC :: n_chemdiags
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
PUBLIC :: asad_flux_put_stash
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

TYPE :: stashsumdiag
  INTEGER :: stash_value               ! section*1000 + item
  INTEGER :: stash_section             ! section
  INTEGER :: stash_item                ! item
  INTEGER :: chemd1codes_location      ! D1 location
  INTEGER :: number_of_fields          ! No of fields
  INTEGER :: len_dim1                  ! Array dimension
  INTEGER :: len_dim2                  ! Array dimension
  INTEGER :: len_dim3                  ! Array dimension
  INTEGER, ALLOCATABLE :: chemdiags_location(:)
END TYPE stashsumdiag

INTEGER :: n_stashsumdiag              ! No of requests

CHARACTER(LEN=10) :: BlankStr = '          '  ! Represents undefined string

! Array to define diagnostic properties
TYPE(chemdiag), ALLOCATABLE, SAVE :: asad_chemdiags(:)

! Array to define D1 lengths, address, etc.
TYPE(code), ALLOCATABLE, SAVE :: chemD1codes(:)

TYPE(stashsumdiag), ALLOCATABLE, SAVE :: stash_handling(:)

INTEGER, SAVE :: n_asad_diags ! No. of ASAD diagnostics actually
                             ! requested in STASH
INTEGER, ALLOCATABLE, SAVE :: st_asad_index(:)
                             ! Index of ASAD diags in full stash list

INTEGER, SAVE :: n_chemdiags

! This is set to .FALSE. initially, and changed in setstash
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

SUBROUTINE asad_setstash_chemdiag(row_length,rows,model_levels)
IMPLICIT NONE


INTEGER, INTENT(IN) :: row_length     ! length of row
INTEGER, INTENT(IN) :: rows           ! number of rows
INTEGER, INTENT(IN) :: model_levels   ! number of levels

INTEGER :: i,ierr,idiag,itm,ilast   ! counters

CHARACTER (LEN=errormessagelength) :: cmessage          ! Error message

LOGICAL, SAVE :: firsttime=.TRUE.     ! T for first call
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_SETSTASH_CHEMDIAG'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (firsttime) THEN

  itm=0
  ilast=0
  n_asad_diags = 0
  ! Allocate index array to total number of diagnostics defined in
  ! asad_chem_fluxes. The latter array contains duplication, so number of
  ! ASAD Stash requests should always be smaller than the allocated size.
  ALLOCATE( st_asad_index(SIZE(asad_chemical_fluxes)) )
  st_asad_index(:) = -99

  ! Loop over the Stash requests and determine number of ASAD diagnostics
  ! requested (i.e. Section-50 diags that are also defined in the
  ! asad_chem_fluxes array), and store their position in full Stash list.
  ! The number and index will be used for all further processing.
  DO idiag = 1, ndiag
    IF ( modl_b(idiag) == submodel_for_sm(atmos_im) .AND.                      &
          isec_b(idiag) == stashcode_ukca_chem_diag ) THEN
      itm = stashcode_ukca_chem_diag*1000+item_b(idiag)
      ! Each ASAD request is recorded just once, even though there
      ! may be multiple instances in the Stash, so compare against the
      ! previous item. This assumes requests are stored in serial order
      ! by STASH
      IF ( itm /= ilast .AND.                                                  &
         ANY(asad_chemical_fluxes(:)%stash_number == itm) ) THEN
        n_asad_diags = n_asad_diags + 1
        st_asad_index(n_asad_diags) = idiag   ! Store position in Stash
        ilast = itm
      END IF
    END IF    ! If isec_b = stashcode_ukca_chem_diag
  END DO      ! Loop over Stash list

  IF ( n_asad_diags > 0 ) THEN
    ALLOCATE(chemD1codes(n_asad_diags))

    ! These are initialised after being read in from D1
    chemD1codes(:)%section    = imdi
    chemD1codes(:)%item       = imdi
    chemD1codes(:)%n_levels   = imdi
    chemD1codes(:)%address    = imdi
    chemD1codes(:)%length     = imdi
    chemD1codes(:)%halo_type  = imdi
    chemD1codes(:)%grid_type  = imdi
    chemD1codes(:)%field_type = imdi
    ! These are standard
    chemD1codes(:)%len_dim1   = row_length
    chemD1codes(:)%len_dim2   = rows
    ! Dry dep is single level - change when assigned, this will be
    ! re-set in init_chemdiag
    chemD1codes(:)%len_dim3   = model_levels
    chemD1codes(:)%prognostic = .TRUE.
    ! Required set after D1 reading
    chemD1codes(:)%required   = .FALSE.
  END IF

  ! this should now set up the model to all the stash requests that are
  ! needed
  n_stashsumdiag = 0

  ! should have already been set above, but re-set to .false. here
  L_asad_use_chem_diags = .FALSE.

  DO i=1,n_asad_diags
    idiag = st_asad_index(i)
    Chemd1codes(i)%section    = stashcode_ukca_chem_diag
    Chemd1codes(i)%item       = item_b(idiag)
    Chemd1codes(i)%len_dim1   = row_length
    Chemd1codes(i)%len_dim2   = rows
    ! check on size of model levels - currently only
    ! contiguous theta levels or surface allowed
    IF (iopl_d(idom_b(idiag)) == st_levels_model_theta) THEN
      Chemd1codes(i)%len_dim3 = model_levels
    ELSE IF (iopl_d(idom_b(idiag)) == st_levels_single) THEN
      Chemd1codes(i)%len_dim3 = 1
    ELSE
      WRITE(cmessage,'(A,5(i6,1X),A)') 'ERROR: ',                              &
              idiag,modl_b(idiag),                                             &
              isec_b(idiag),item_b(idiag),                                     &
              iopl_d(idom_b(idiag)),                                           &
              ' incorrect model levels'
      ierr = isec_b(idiag)*1000 + item_b(idiag)
      CALL ereport('ASAD_SETSTASH_CHEMDIAG',ierr,                              &
                   cmessage)
    END IF
    ! required is set to true here, but in the context of a
    ! UKCA framework, it would be set to .false.
    Chemd1codes(i)%required   = .TRUE.
    IF (.NOT. L_asad_use_chem_diags)                                           &
        L_asad_use_chem_diags = .TRUE.
    Chemd1codes(i)%prognostic = .FALSE.
    n_stashsumdiag              = n_stashsumdiag + 1

    IF (PrintStatus >= Prstatus_Oper .AND.                                     &
             Chemd1codes(i)%section /= imdi) THEN
      WRITE(umMessage,'(A,6(I5,1x),L1)') 'ASAD_FLUXES: ',i,                    &
           Chemd1codes(i)%section, Chemd1codes(i)%item,                        &
           Chemd1codes(i)%len_dim1, Chemd1codes(i)%len_dim2,                   &
           Chemd1codes(i)%len_dim3, Chemd1codes(i)%required
      CALL umPrint(umMessage,src='asad_chem_flux_diags')
    END IF
  END DO       ! I=1,n_asad_diags

  firsttime = .FALSE.
END IF ! firstime

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
SUBROUTINE asad_init_chemdiag(ierr)

IMPLICIT NONE

INTEGER, INTENT(IN OUT) :: ierr       ! Error code

INTEGER :: i,j,k

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


ALLOCATE(stash_handling(1:n_stashsumdiag))
stash_handling(:)%stash_value   = 0
stash_handling(:)%stash_section = 0
stash_handling(:)%stash_item    = 0
stash_handling(:)%chemd1codes_location = 0
stash_handling(:)%len_dim1 = 0
stash_handling(:)%len_dim2 = 0
stash_handling(:)%len_dim3 = 0
stash_handling(:)%number_of_fields = 0

! Allocate diags array since we know how large it will
! need to be from STASH
! get stash item numbers - used to check/ch
! check to make sure they are required
n_chemdiags = 0
icount=0
DO i=1,n_asad_diags
  IF (chemD1codes(i)%required) THEN
    icount=icount+1
    stash_handling(icount)%stash_section =                                     &
              chemd1codes(i)%section
    stash_handling(icount)%stash_item =                                        &
              chemd1codes(i)%item
    stash_handling(icount)%stash_value =                                       &
              ((stash_handling(icount)%stash_section)*1000) +                  &
              stash_handling(icount)%stash_item
    stash_handling(icount)%chemd1codes_location = i
    stash_handling(icount)%len_dim1 = chemd1codes(i)%len_dim1
    stash_handling(icount)%len_dim2 = chemd1codes(i)%len_dim2
    stash_handling(icount)%len_dim3 =                                          &
              MAX(chemd1codes(i)%len_dim3,                                     &
             stash_handling(icount)%len_dim3)
    ! calculate how many diagnostics from flux_dat we are using
    DO j=1,SIZE(asad_chemical_fluxes)
      IF (stash_handling(icount)%stash_value ==                                &
                 asad_chemical_fluxes(j)%stash_number) THEN
        n_chemdiags = n_chemdiags + 1
        ! NOTE: do not exit after first check, since we may have
        !       more than one diagnostic associated with each
        !       stash number.
      END IF
    END DO
    ! debug output

  END IF      ! chemD1codes(i)%required
END DO        ! i=1,n_asad_diags

IF (L_asad_use_chem_diags) THEN

      ! allocate diagnostics array
  ALLOCATE(asad_chemdiags(1:n_chemdiags))
  DO i=1,n_chemdiags
    CALL asad_zero_chemdiag(asad_chemdiags(i))
  END DO

  ! Now we do a rather complicated do-loop to set up the asad_chemdiags array
  jdiag = 0
  icount=0
  DO istash=1,n_asad_diags
     ! check to see if stash is requested
    IF (chemD1codes(istash)%required) THEN
      icount=icount+1
      DO j=1,SIZE(asad_chemical_fluxes)
         ! check if stash code has been requested, and if
         ! we have a corresponding flux specified
        IF (stash_handling(icount)%stash_value ==                              &
             asad_chemical_fluxes(j)%stash_number) THEN
           ! increment jdiag. We may have more than one
           ! diagnostic associated with each stash number
          jdiag = jdiag + 1
          ! STASH number
          asad_chemdiags(jdiag)%stash_number=                                  &
               asad_chemical_fluxes(j)%stash_number
          ! diagnostic type
          asad_chemdiags(jdiag)%diag_type =                                    &
               TRIM(asad_cd_uppercase(                                         &
                    asad_chemical_fluxes(j)%diag_type))
          ! reaction (or otherwise) type
          asad_chemdiags(jdiag)%rxn_type =                                     &
               TRIM(asad_cd_uppercase(                                         &
                    asad_chemical_fluxes(j)%rxn_type))
          ! do we use a tropospheric mask only output the
          ! troposphere?
          IF (asad_chemical_fluxes(                                            &
               j)%tropospheric_mask) THEN
            asad_chemdiags(jdiag)%tropospheric_mask                            &
                 = .TRUE.
            ! set logical for global value
            IF (.NOT. L_asad_use_trop_mask)                                    &
                 L_asad_use_trop_mask = .TRUE.
          ELSE
            asad_chemdiags(jdiag)%tropospheric_mask                            &
                 = .FALSE.
          END IF
          ! if we have two reactions with the exactly the same
          ! reactants and products, but different rates,
          ! find_reaction will only find the first one, this
          ! flag over-rides that
          IF (asad_chemical_fluxes(j)%rxn_location                             &
               <= 0) THEN
            asad_chemdiags(jdiag)%find_rxn_loc = 1
          ELSE
            asad_chemdiags(jdiag)%find_rxn_loc =                               &
                 asad_chemical_fluxes(j)%rxn_location
          END IF
          ! fill in number of reactants and products
          IF (asad_chemical_fluxes(j)%num_species                              &
              == 1) THEN
             ! are not dealing with a reaction
            asad_chemdiags(jdiag)%num_reactants = 0
            asad_chemdiags(jdiag)%num_products  = 0
          ELSE IF (asad_chemical_fluxes(j)%num_species                         &
               >= 3) THEN
             ! we have a reaction (probably, or user error!)
            asad_chemdiags(jdiag)%num_reactants =                              &
                 number_of_reactants ! = 2 currently
            asad_chemdiags(jdiag)%num_products  =                              &
                 asad_chemical_fluxes(j)%num_species -                         &
                 asad_chemdiags(jdiag)%num_reactants
          ELSE ! incorrect numbers specified!
            WRITE(cmessage,'(A,I0,A)')                                         &
                 'WRONG NUMBER OF SPECIES: ',                                  &
                 asad_chemical_fluxes(j)%num_species,                          &
                 ' FOUND'
            errcode=asad_chemdiags(jdiag)%stash_number
            WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',         &
                 errcode,' PE: ',mype
            CALL umPrint(umMessage,src='asad_chem_flux_diags')
            CALL ereport('ASAD_CHEMICAL_DIAGNOSTICS',                          &
                 errcode,cmessage)
          END IF
          ! initialise reactants and products
          asad_chemdiags(jdiag)%reactants(:) = BlankStr
          IF (ALLOCATED(asad_chemdiags(                                        &
               jdiag)%products))                                               &
               DEALLOCATE(asad_chemdiags(                                      &
               jdiag)%products)
          ALLOCATE(asad_chemdiags(jdiag)%products(1:                           &
               asad_chemdiags(jdiag)%num_products))
          asad_chemdiags(jdiag)%products(:) = BlankStr
          ! fill in from list - rxn specific
          SELECT CASE (asad_chemdiags(jdiag)%diag_type)
          CASE (cdrxn,cdrte)
             ! copy over reacants
            DO i=1,number_of_reactants
              asad_chemdiags(jdiag)%reactants(i) =                             &
                   TRIM(ADJUSTL(                                               &
                   asad_chemical_fluxes(                                       &
                   j)%reactants(i)))
            END DO
            ! copy over products
            DO i=1,asad_chemdiags(jdiag)%num_products
              asad_chemdiags(jdiag)%products(i) =                              &
                   TRIM(ADJUSTL(                                               &
                   asad_chemical_fluxes(                                       &
                   j)%products(i)))
            END DO
          CASE (cddep,cdems,cdnet,cdste,cdout)
             ! only one species in this case
            asad_chemdiags(jdiag)%species =                                    &
                 TRIM(ADJUSTL(                                                 &
                 asad_chemical_fluxes(j)%reactants(1)))
          CASE (cdpsc,cdmas,cdtpm,cdlgt)
             ! no need to do anything - no species
          CASE DEFAULT
            cmessage='DIAGNOSTIC TYPE '//                                      &
                 asad_chemdiags(jdiag)%diag_type                               &
                 //' NOT FOUND'
            errcode=asad_chemdiags(jdiag)%stash_number
            WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',         &
                 errcode,' PE: ',mype
            CALL umPrint(umMessage,src='asad_chem_flux_diags')
            CALL ereport('ASAD_INIT_CHEMDIAG',                                 &
                    errcode,cmessage)
          END SELECT

          ! Check on emissions, deposition etc, and set allocate/output logicals
          SELECT CASE (asad_chemdiags(jdiag)%diag_type)
          CASE (cdems)
             ! emissions only on chemical timesteps, and can
             ! deallocate at other times
            asad_chemdiags(jdiag)%can_deallocate                               &
                 = .TRUE.
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
              cmessage='INCORRECT EMISSION TYPE: '//                           &
                   asad_chemdiags(jdiag)%diag_type//                           &
                  ':'//asad_chemdiags(jdiag)%rxn_type
              errcode=                                                         &
                   asad_chemdiags(jdiag)%stash_number
              WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',       &
                   errcode,' PE: ',mype
              CALL umPrint(umMessage,src='asad_chem_flux_diags')
              CALL ereport('ASAD_INIT_CHEMDIAG',                               &
                   errcode,cmessage)
            END SELECT ! asad_chemdiags(jdiag)%rxn_type
          CASE (cdrxn)
             ! reactions only on chemical timesteps, and can
             ! deallocate at other times
            asad_chemdiags(jdiag)%can_deallocate = .TRUE.
            SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
            CASE (cdbimol,cdtermol,cdphot,cdhetero)
              L_asad_use_flux_rxns=.TRUE.
            CASE DEFAULT
              cmessage='INCORRECT REACTION TYPE: '//                           &
                   asad_chemdiags(jdiag)%diag_type//                           &
                   ':'//                                                       &
                   asad_chemdiags(jdiag)%rxn_type
              errcode=asad_chemdiags(                                          &
                   jdiag)%stash_number
              WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',       &
                   errcode,' PE: ',mype
              CALL umPrint(umMessage,src='asad_chem_flux_diags')
              CALL ereport('ASAD_CHEMICAL_DIAGNOSTICS',                        &
                   errcode,cmessage)
            END SELECT ! asad_chemdiags(jdiag)%rxn_type
          CASE (cdrte)
             ! rates only on chemical timesteps, and can
             ! deallocate at other times
            asad_chemdiags(jdiag)%can_deallocate                               &
                 = .TRUE.
            SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
            CASE (cdbimol,cdtermol,cdphot,cdhetero)
              L_asad_use_rxn_rates=.TRUE.
            CASE DEFAULT
              cmessage='INCORRECT RATE TYPE: '//                               &
                   asad_chemdiags(jdiag)%diag_type//                           &
                   ':'//                                                       &
                   asad_chemdiags(jdiag)%rxn_type
              errcode=asad_chemdiags(                                          &
                   jdiag)%stash_number
              WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',       &
                   errcode,' PE: ',mype
              CALL umPrint(umMessage,src='asad_chem_flux_diags')
              CALL ereport('ASAD_INIT_CHEMDIAG',                               &
                   errcode,cmessage)
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
              cmessage='INCORRECT DEPOSTITION TYPE: '//                        &
                   asad_chemdiags(jdiag)%diag_type//                           &
                   ':'//                                                       &
                   asad_chemdiags(jdiag)%rxn_type
              errcode=asad_chemdiags(                                          &
                   jdiag)%stash_number
              WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',       &
                   errcode,' PE: ',mype
              CALL umPrint(umMessage,src='asad_chem_flux_diags')
              CALL ereport('ASAD_INIT_CHEMDIAG',                               &
                   errcode,cmessage)
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
            L_asad_use_trop_mask_output=.TRUE.
          CASE (cdout)
             ! tracer on all timesteps, and cannot be
             ! deallocated at any time
            asad_chemdiags(jdiag)%can_deallocate                               &
                 = .TRUE.
            ! tracer
            L_asad_use_output_tracer = .TRUE.
          CASE (cdlgt)
             ! lightning diagnostics on all timesteps, and cannot
             ! be deallocated at any time
            asad_chemdiags(jdiag)%can_deallocate                               &
                 = .TRUE.
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
              cmessage=                                                        &
                   'INCORRECT LIGHTNING DIAG TYPE: '//                         &
                   asad_chemdiags(jdiag)%diag_type//                           &
                   ':'//                                                       &
                   asad_chemdiags(jdiag)%rxn_type
              errcode=                                                         &
                   asad_chemdiags(jdiag)%stash_number
              WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',       &
                   errcode,' PE: ',mype
              CALL umPrint(umMessage,src='asad_chem_flux_diags')
              CALL ereport(                                                    &
                   'ASAD_CHEMICAL_DIAGNOSTICS',                                &
                   errcode,cmessage)
            END SELECT
          END SELECT ! asad_chemdiags(jdiag)%diag_type
        END IF ! stash_handling(icount)%stash_value ==
               ! asad_chemical_fluxes(j)%stash_number
      END DO !j,SIZE(asad_chemical_diagnostics)
    END IF ! chemd1codes(i)%required
  END DO ! i,n_asad_diags

  ! Check on summing options - simplifies STASH output
  DO k=1,n_stashsumdiag
    DO i=1,n_chemdiags
       ! calcs number of fields with the same stash code
      IF (stash_handling(k)%stash_value ==                                     &
           asad_chemdiags(i)%stash_number) THEN
        stash_handling(k)%number_of_fields =                                   &
             stash_handling(k)%number_of_fields + 1
      END IF
    END DO
    ! creates an array of this size
    ALLOCATE(stash_handling(k)%chemdiags_location(1:                           &
         stash_handling(k)%number_of_fields))
    stash_handling(k)%chemdiags_location(:) = 0
    ! allocate %throughput array for each flux
    icount=0
    DO i=1,n_chemdiags
      IF (stash_handling(k)%stash_value ==                                     &
           asad_chemdiags(i)%stash_number) THEN
        icount=icount+1
        stash_handling(k)%chemdiags_location(icount) = i
        ! get number of levels (needed in putstash)
        asad_chemdiags(i)%num_levs =                                           &
             stash_handling(k)%len_dim3
      END IF
    END DO           ! i,n_chemdiags
  END DO              ! k,n_stashsumdiag
  ! error check
  DO k=1,n_stashsumdiag
    IF (stash_handling(k)%number_of_fields == 0) THEN
      cmessage=' STASH CODE NOT FOUND IN'//                                    &
               ' stash_handling ARRAY'
      errcode = -1*stash_handling(k)%stash_value
      WRITE(umMessage,'(2A,I0,A,I0,A,I0)') cmessage, ' Error code: ',          &
           errcode,' PE: ',mype,' K: ',k
      CALL umPrint(umMessage,src='asad_chem_flux_diags')
      CALL ereport('ASAD_INIT_CHEMDIAG',errcode,cmessage)
    END IF
  END DO

END IF                 ! L_asad_use_chem_diags



! print out diagnostics requested
IF ((PrintStatus >= PrStatus_Oper ) .AND. (mype == 0)) THEN
  WRITE(umMessage,'(A,A)') 'THE FOLLOWING UKCA FLUX DIAGNOSTICS ',             &
       'HAVE BEEN REQUESTED'
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  DO i=1,n_chemdiags
    IF (asad_chemdiags(i)%num_products > 0) THEN
       ! fancy formatting to cope with variable number of prods
      IF (asad_chemdiags(i)%num_products < 10) THEN
        WRITE(outformat,FMT='(A38,I1,A7)')                                     &
          '(I5,1X,A3,1X,A1,1X,A,1X,2(A,1X),A2,1X,',                            &
          asad_chemdiags(i)%num_products,                                      &
          '(A,1X))'
      ELSE
        WRITE(outformat,FMT='(A38,I2,A7)')                                     &
             '(I5,1X,A3,1X,A1,1X,A,1X,2(A,1X),A2,1X,',                         &
             asad_chemdiags(i)%num_products,                                   &
             '(A,1X))'
      END IF
      WRITE(umMessage,FMT=TRIM(ADJUSTL(outformat)))                            &
           asad_chemdiags(i)%stash_number,                                     &
           asad_chemdiags(i)%diag_type,                                        &
           asad_chemdiags(i)%rxn_type,                                         &
           TRIM(ADJUSTL(asad_chemdiags(i)%species)),                           &
           (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))),                     &
           j=1,2),'->',                                                        &
           (TRIM(ADJUSTL(asad_chemdiags(i)%products(j))),                      &
           j=1,asad_chemdiags(i)%num_products)
      CALL umPrint(umMessage,src='asad_chem_flux_diags')
    ELSE
       ! not a reaction (or rather, has no products)
      WRITE(umMessage,FMT='(I5,1X,A3,1X,A1,1X,A,1X,2(A,1X))')                  &
           asad_chemdiags(i)%stash_number,                                     &
           asad_chemdiags(i)%diag_type,                                        &
           asad_chemdiags(i)%rxn_type,                                         &
           TRIM(ADJUSTL(asad_chemdiags(i)%species)),                           &
           (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))),                     &
           j=1,2)
      CALL umPrint(umMessage,src='asad_chem_flux_diags')
    END IF
  END DO
END IF

! print out status of logicals
IF (PrintStatus == Prstatus_Diag) THEN
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
       'L_asad_use_air_ems         = ',L_asad_use_air_ems
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
       'L_asad_use_light_ems       = ',L_asad_use_light_ems
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
       'L_asad_use_surf_ems        = ',L_asad_use_surf_ems
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
       'L_asad_use_volc_ems        = ',L_asad_use_volc_ems
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
       'L_asad_use_sulp_ems        = ',L_asad_use_sulp_ems
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
       'L_asad_use_flux_rxns       = ',L_asad_use_flux_rxns
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
       'L_asad_use_wetdep          = ',L_asad_use_wetdep
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
  'L_asad_use_drydep          = ',L_asad_use_drydep
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
       'L_asad_use_tendency        = ',L_asad_use_tendency
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
       'L_asad_use_mass_diagnostic = ',                                        &
       L_asad_use_mass_diagnostic
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
       'L_asad_use_trop_mask       = ',L_asad_use_trop_mask
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  WRITE(umMessage,'(A,A,L1)') ' asad_init_chemdiag: ',                         &
       'L_asad_use_chem_diags      = ',L_asad_use_chem_diags
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
END IF ! PrintStatus = Prstatus_diag

ierr = 0

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

INTEGER :: i


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_ALLOCATE_CHEMDIAG'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate flux array and initialise
DO i=1,n_chemdiags
  IF (.NOT. ALLOCATED(asad_chemdiags(i)%throughput)) THEN
    ALLOCATE(asad_chemdiags(i)%throughput(                                     &
             1:row_length,                                                     &
             1:rows,                                                           &
             1:asad_chemdiags(i)%num_levs))
    asad_chemdiags(i)%throughput(:,:,:) = 0.0
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_allocate_chemdiag

! #####################################################################
SUBROUTINE asad_chemical_diagnostics(row_length, rows,                         &
     model_levels, dpd_full, dpw_full, prk_full, y_full,                       &
     ix, jy, klevel, volume, ierr)

USE chemistry_constants_mod,  ONLY: avogadro
USE ukca_config_specification_mod, ONLY: ukca_config

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length       ! length of row
INTEGER, INTENT(IN) :: rows             ! number of rows
INTEGER, INTENT(IN) :: model_levels     ! the level that we are at
INTEGER, INTENT(IN) :: ix               ! i counter
INTEGER, INTENT(IN) :: jy               ! j counter
INTEGER, INTENT(IN) :: klevel           ! the level that we are at
REAL, INTENT(IN)    :: volume(row_length, rows, model_levels)   ! cell volume
INTEGER, INTENT(IN OUT) :: ierr          ! error code
REAL, INTENT(IN)    :: dpd_full(model_levels,jpspec)
REAL, INTENT(IN)    :: dpw_full(model_levels,jpspec)
REAL, INTENT(IN)    :: prk_full(model_levels,jpnr)
REAL, INTENT(IN)    :: y_full(model_levels,jpspec)

REAL, PARAMETER :: convfac = 1.0e6/avogadro ! conversion factor to convert
                                            ! volume into cm^3 and result in mol

LOGICAL, SAVE :: firstcall=.TRUE.
LOGICAL, SAVE :: first_pass=.TRUE.

INTEGER :: i,j,idep                     ! counters

INTEGER, ALLOCATABLE :: findrxn_tmp(:)   ! array to hold locations

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=72) :: outformat
INTEGER       :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_CHEMICAL_DIAGNOSTICS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr = -1

IF (first_pass) THEN
  ! OMP CRITICAL will only allow one thread through this code at a time,
  ! while the other threads are held until completion.
!$OMP CRITICAL (asad_chemical_diagnostics_init)
  IF (firstcall) THEN

    ! Set up locations of reactions using (an altered version of)
    ! asad_findreaction
    DO i=1,n_chemdiags
      SELECT CASE (asad_chemdiags(i)%diag_type)
      CASE (cdrxn,cdrte)
        SELECT CASE (asad_chemdiags(i)%rxn_type)
        CASE (cdbimol) ! bimolecular reacions/rate
          ALLOCATE(findrxn_tmp(1:(jpbk+1)))
          findrxn_tmp = cd_findreaction(                                       &
               asad_chemdiags(i)%num_products,                                 &
               asad_chemdiags(i)%reactants,                                    &
               asad_chemdiags(i)%products,                                     &
               spb, nbrkx, (jpbk+1), jpspb )
          asad_chemdiags(i)%location =                                         &
               findrxn_tmp(asad_chemdiags(i)%find_rxn_loc)
          DEALLOCATE(findrxn_tmp)
          IF (asad_chemdiags(i)%location <= 0) THEN
            IF (asad_chemdiags(i)%num_products < 10) THEN
              WRITE(outformat,FMT='(A28,I1,A7)')                               &
                   '(A13,1X,A1,1X,2(A,1X),A2,1X,',                             &
                   asad_chemdiags(i)%num_products,                             &
                   '(A,1X))'
            ELSE
              WRITE(outformat,FMT='(A28,I2,A7)')                               &
                   '(A13,1X,A1,1X,2(A,1X),A2,1X,',                             &
                   asad_chemdiags(i)%num_products,                             &
                   '(A,1X))'
            END IF
            WRITE(umMessage,'(A)') outformat
            CALL umPrint(umMessage,src='asad_chem_flux_diags')
            WRITE(cmessage,FMT=TRIM(ADJUSTL(outformat)))                       &
              'RXN NOT FOUND', cdbimol,                                        &
              (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))),                  &
               j=1,2),'->',                                                    &
              (TRIM(ADJUSTL(asad_chemdiags(i)%products(j))),                   &
               j=1,asad_chemdiags(i)%num_products)
            errcode=asad_chemdiags(i)%stash_number
            WRITE(umMessage,'(A,A13,I5,A5,I6)')                                &
              cmessage, ' Error code: ',errcode, ' PE: ',mype
            CALL umPrint(umMessage,src='asad_chem_flux_diags')
            CALL ereport('ASAD_CHEMICAL_DIAGNOSTICS',                          &
                 errcode,cmessage)
          END IF
        CASE (cdhetero) ! Heterogeneous reactions/rates
          ALLOCATE(findrxn_tmp(1:(jphk+1)))
          findrxn_tmp = cd_findreaction(                                       &
               asad_chemdiags(i)%num_products,                                 &
               asad_chemdiags(i)%reactants,                                    &
               asad_chemdiags(i)%products,                                     &
               sph, nhrkx, (jphk+1), jpsph )
          asad_chemdiags(i)%location =                                         &
               findrxn_tmp(asad_chemdiags(i)%find_rxn_loc)
          DEALLOCATE(findrxn_tmp)
          IF (asad_chemdiags(i)%location <= 0) THEN
            IF (asad_chemdiags(i)%num_products < 10) THEN
              WRITE(outformat,FMT='(A28,I1,A7)')                               &
                   '(A13,1X,A1,1X,2(A,1X),A2,1X,',                             &
                   asad_chemdiags(i)%num_products,                             &
                   '(A,1X))'
            ELSE
              WRITE(outformat,FMT='(A28,I2,A7)')                               &
                   '(A13,1X,A1,1X,2(A,1X),A2,1X,',                             &
                   asad_chemdiags(i)%num_products,                             &
                   '(A,1X))'
            END IF
            WRITE(cmessage,FMT=TRIM(ADJUSTL(outformat)))                       &
               'RXN NOT FOUND', cdhetero,                                      &
               (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))),                 &
               j=1,2),'->',                                                    &
               (TRIM(ADJUSTL(asad_chemdiags(i)%products(j))),                  &
               j=1,asad_chemdiags(i)%num_products)
            errcode=asad_chemdiags(i)%stash_number
            CALL ereport('ASAD_CHEMICAL_DIAGNOSTICS',                          &
                 errcode,cmessage)
          END IF
        CASE (cdtermol) ! termolecular reactions/rates
          ALLOCATE(findrxn_tmp(1:(jptk+1)))
          findrxn_tmp = cd_findreaction(                                       &
               asad_chemdiags(i)%num_products,                                 &
               asad_chemdiags(i)%reactants,                                    &
               asad_chemdiags(i)%products,                                     &
               spt, ntrkx, (jptk+1), jpspt )
          asad_chemdiags(i)%location =                                         &
               findrxn_tmp(asad_chemdiags(i)%find_rxn_loc)
          DEALLOCATE(findrxn_tmp)
          IF (asad_chemdiags(i)%location <= 0) THEN
            IF (asad_chemdiags(i)%num_products < 10) THEN
              WRITE(outformat,FMT='(A28,I1,A7)')                               &
                   '(A13,1X,A1,1X,2(A,1X),A2,1X,',                             &
                   asad_chemdiags(i)%num_products,                             &
                   '(A,1X))'
            ELSE
              WRITE(outformat,FMT='(A28,I2,A7)')                               &
                   '(A13,1X,A1,1X,2(A,1X),A2,1X,',                             &
                   asad_chemdiags(i)%num_products,                             &
                   '(A,1X))'
            END IF
            WRITE(cmessage,FMT=TRIM(ADJUSTL(outformat)))                       &
               'RXN NOT FOUND', cdtermol,                                      &
               (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))),                 &
               j=1,2),'->',                                                    &
               (TRIM(ADJUSTL(asad_chemdiags(i)%products(j))),                  &
               j=1,asad_chemdiags(i)%num_products)
            errcode=asad_chemdiags(i)%stash_number
            WRITE(umMessage,'(A,A13,I5,A5,I6)')                                &
               cmessage, ' Error code: ',errcode, ' PE: ',mype
            CALL umPrint(umMessage,src='asad_chem_flux_diags')
            CALL ereport('ASAD_CHEMICAL_DIAGNOSTICS',                          &
                 errcode,cmessage)
          END IF
        CASE (cdphot) ! photolysis reactions/rates
          ALLOCATE(findrxn_tmp(1:(jppj+1)))
          findrxn_tmp = cd_findreaction(                                       &
               asad_chemdiags(i)%num_products,                                 &
               asad_chemdiags(i)%reactants,                                    &
               asad_chemdiags(i)%products,                                     &
               spj, nprkx, (jppj+1), jpspj )
          asad_chemdiags(i)%location =                                         &
               findrxn_tmp(asad_chemdiags(i)%find_rxn_loc)
          DEALLOCATE(findrxn_tmp)
          IF (asad_chemdiags(i)%location <= 0) THEN
            IF (asad_chemdiags(i)%num_products < 10) THEN
              WRITE(outformat,FMT='(A28,I1,A7)')                               &
                   '(A13,1X,A1,1X,2(A,1X),A2,1X,',                             &
                   asad_chemdiags(i)%num_products,                             &
                   '(A,1X))'
            ELSE
              WRITE(outformat,FMT='(A28,I2,A7)')                               &
                   '(A13,1X,A1,1X,2(A,1X),A2,1X,',                             &
                   asad_chemdiags(i)%num_products,                             &
                   '(A,1X))'
            END IF
            WRITE(cmessage,FMT=TRIM(ADJUSTL(outformat)))                       &
               'RXN NOT FOUND', cdphot,                                        &
               (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))),                 &
               j=1,2),'->',                                                    &
               (TRIM(ADJUSTL(asad_chemdiags(i)%products(j))),                  &
               j=1,asad_chemdiags(i)%num_products)
            errcode=asad_chemdiags(i)%stash_number
            WRITE(umMessage,'(A,A13,I5,A5,I6)')                                &
               cmessage, ' Error code: ',errcode, ' PE: ',mype
            CALL umPrint(umMessage,src='asad_chem_flux_diags')
            CALL ereport('ASAD_CHEMICAL_DIAGNOSTICS',                          &
                 errcode,cmessage)
          END IF
        CASE DEFAULT
          cmessage='REACTION TYPE '//                                          &
               asad_chemdiags(i)%rxn_type//' NOT FOUND'
          errcode=asad_chemdiags(i)%stash_number
          WRITE(umMessage,'(A,A13,I5,A5,I6)')                                  &
             cmessage, ' Error code: ',errcode, ' PE: ',mype
          CALL umPrint(umMessage,src='asad_chem_flux_diags')
          CALL ereport('ASAD_CHEMICAL_DIAGNOSTICS',                            &
               errcode,cmessage)
        END SELECT       ! asad_chemdiags(i)%rxn_type
      CASE (cddep)
        SELECT CASE (asad_chemdiags(i)%rxn_type)
        CASE (cddry) ! DRY DEP
          ! %species contains species to be dry deposited
          ! %location contains location in requisite array
          idep=1
          DO j=1,jpspec
            IF (ldepd(j) .AND. (asad_chemdiags(i)%species ==                   &
                 speci(j))) THEN
                ! it is being deposited logically and in diagnostic sense
              asad_chemdiags(i)%location = j
              idep=idep+1
            ELSE IF (.NOT. ldepd(j) .AND.                                      &
               (asad_chemdiags(i)%species == speci(j))) THEN
              ! diagnostics of deposition requested, but species is not
              ! labelled as being dry deposited in the chch_defs file
              cmessage='SPECIES '//asad_chemdiags(i)%species                   &
                     //' NOT DRY DEPOSITED'
              errcode=asad_chemdiags(i)%stash_number
              WRITE(umMessage,'(A,A13,I5,A5,I6)')                              &
                 cmessage, ' Error code: ', errcode,' PE: ',mype
              CALL umPrint(umMessage,src='asad_chem_flux_diags')
              CALL ereport('ASAD_CHEMICAL_DIAGNOSTICS',                        &
                   errcode,cmessage)
            END IF
          END DO
        CASE (cdwet) ! WET DEP
          DO j=1,jpspec
            IF (ldepw(j) .AND. (asad_chemdiags(i)%species ==                   &
                 speci(j))) THEN
              asad_chemdiags(i)%location = j
              idep=idep+1
            ELSE IF (.NOT. ldepw(j) .AND.                                      &
               (asad_chemdiags(i)%species == speci(j))) THEN
               ! diagnostics of deposition requested, but species is not
              ! labelled as being wet deposited in the chch_defs file
              cmessage='SPECIES '//asad_chemdiags(i)%species                   &
                      //' NOT WET DEPOSITED'
              errcode=asad_chemdiags(i)%stash_number
              WRITE(umMessage,'(A,A13,I5,A5,I6)')                              &
                 cmessage, ' Error code: ', errcode,' PE: ',mype
              CALL umPrint(umMessage,src='asad_chem_flux_diags')
              CALL ereport('ASAD_CHEMICAL_DIAGNOSTICS',                        &
                   errcode,cmessage)
            END IF
          END DO
          ! %species contains species to be wet deposited
          ! %location contains location in requisite array
        END SELECT    ! asad_chemdiags(i)%rxn_type
      END SELECT    ! asad_chemdiags(i)%diag_type


      ! assumes 2 reactants
      IF (PrintStatus == Prstatus_diag) THEN
        SELECT CASE (asad_chemdiags(i)%diag_type)
        CASE (cdrxn,cdrte)
          WRITE(umMessage,'(A5,I6,A10,I5,A1,A1,A3,A10,A3,A10,A7,I6)')          &
               'mype=',mype,' CD:FIND: ',                                      &
               asad_chemdiags(i)%stash_number,' ',                             &
               asad_chemdiags(i)%rxn_type,' : ',                               &
               asad_chemdiags(i)%reactants(1),                                 &
               ' + ',asad_chemdiags(i)%reactants(2),                           &
               ' rxn = ',                                                      &
               asad_chemdiags(i)%location
          CALL umPrint(umMessage,src='asad_chem_flux_diags')
        CASE (cddep)
          WRITE(umMessage,'(A5,I6,A10,I5,A1,A1,A3,A10,A13,I6)')                &
               'mype=',mype,' CD:FIND: ',                                      &
               asad_chemdiags(i)%stash_number,' ',                             &
               asad_chemdiags(i)%rxn_type,' : ',                               &
               asad_chemdiags(i)%species,                                      &
               ' deposition: ',                                                &
               asad_chemdiags(i)%location
          CALL umPrint(umMessage,src='asad_chem_flux_diags')
        END SELECT
      END IF           ! Printstatus
    END DO ! i
    firstcall=.FALSE.
  END IF   ! firstcall
  first_pass=.FALSE.
!$OMP END CRITICAL (asad_chemical_diagnostics_init)
END IF     ! first_pass

! Go through and pick up fluxes from ASAD arrays
! prk is in units of molecules.cm^-3.s^-1
DO i=1,n_chemdiags

  IF (ukca_config%l_ukca_asad_columns) THEN
     ! in this case ASAD is being called column-by column so we need
     ! to fill %throughput by iterating over X & Y
    SELECT CASE (asad_chemdiags(i)%diag_type)
    CASE (cdrxn)
      asad_chemdiags(i)%throughput(ix,jy,:) =                                  &
           prk_full(:,asad_chemdiags(i)%location)*volume(ix,jy,:)*convfac
      IF (asad_chemdiags(i)%tropospheric_mask) THEN
        WHERE (.NOT. L_troposphere(ix,jy,:))
          asad_chemdiags(i)%throughput(ix,jy,:) = 0.0
        END WHERE
      END IF
    CASE (cddep)
      SELECT CASE (asad_chemdiags(i)%rxn_type)
      CASE (cddry) ! DRY DEP
         ! dpd array is over a column, no need to consider klevel
         ! NOTE: Current STASH settings have dry-deposition output
         !       as 3D. If this changes the code below will no
         !       longer work as the code expects an entire
         !       atmospheric column.
        asad_chemdiags(i)%throughput(ix,jy,:)=                                 &
             dpd_full(:,asad_chemdiags(i)%location)*                           &
             y_full(:,asad_chemdiags(i)%location)*                             &
             volume(ix,jy,:)*convfac
        ! Not needed (probably) but kept for consistency
        IF (asad_chemdiags(i)%tropospheric_mask) THEN
          WHERE (.NOT. L_troposphere(ix,jy,:))
            asad_chemdiags(i)%throughput(ix,jy,:) = 0.0
          END WHERE
        END IF
      CASE (cdwet) ! WET DEP
        asad_chemdiags(i)%throughput(ix,jy,:)=                                 &
             dpw_full(:,asad_chemdiags(i)%location)*                           &
             y_full(:,asad_chemdiags(i)%location)*                             &
             volume(ix,jy,:)*convfac
        ! Not needed (probably) but kept for consistency
        IF (asad_chemdiags(i)%tropospheric_mask) THEN
          WHERE (.NOT. L_troposphere(ix,jy,:))
            asad_chemdiags(i)%throughput(ix,jy,:) = 0.0
          END WHERE
        END IF
      END SELECT
    END SELECT
  ELSE
     ! in this case ASAD is being called horizontally so we need
     ! to fill %throughput by iterating over Z
    SELECT CASE (asad_chemdiags(i)%diag_type)
    CASE (cdrxn)
      asad_chemdiags(i)%throughput(:,:,klevel) =                               &
           RESHAPE(prk(:,asad_chemdiags(i)%location),                          &
           [row_length,rows])*volume(:,:,klevel)*convfac
      IF (asad_chemdiags(i)%tropospheric_mask) THEN
        WHERE (.NOT. L_troposphere(:,:,klevel))
          asad_chemdiags(i)%throughput(:,:,klevel) = 0.0
        END WHERE
      END IF
    CASE (cddep)
      SELECT CASE (asad_chemdiags(i)%rxn_type)
      CASE (cddry) ! DRY DEP
        IF (klevel <=                                                          &
             SIZE(asad_chemdiags(i)%throughput(:,:,:),DIM=3)) THEN
           ! if 2D then only take lowest level, otherwise will be 3D
          asad_chemdiags(i)%throughput(:,:,klevel)=                            &
               RESHAPE(dpd(:,asad_chemdiags(i)%location),                      &
                       [row_length,rows])*                                     &
               RESHAPE(y(:,asad_chemdiags(i)%location),                        &
                       [row_length,rows])*                                     &
               volume(:,:,klevel)*convfac
          IF (asad_chemdiags(i)%tropospheric_mask) THEN
            WHERE (.NOT. L_troposphere(:,:,klevel))
              asad_chemdiags(i)%throughput(:,:,klevel) = 0.0
            END WHERE
          END IF
        END IF
        ! Not needed (?)
      CASE (cdwet) ! WET DEP
        asad_chemdiags(i)%throughput(:,:,klevel)=                              &
             RESHAPE(dpw(:,asad_chemdiags(i)%location),                        &
                     [row_length,rows])*                                       &
             RESHAPE(y(:,asad_chemdiags(i)%location),                          &
                     [row_length,rows])*                                       &
             volume(:,:,klevel)*convfac
        ! Not needed (?)
        IF (asad_chemdiags(i)%tropospheric_mask) THEN
          WHERE (.NOT. L_troposphere(:,:,klevel))
            asad_chemdiags(i)%throughput(:,:,klevel) = 0.0
          END WHERE
        END IF
      END SELECT
    END SELECT
  END IF             ! l_ukca_asad_columns
END DO                ! i=1,n_chemdiags


ierr = 0

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
  CALL umPrint(umMessage,src='asad_chem_flux_diags')
  CALL ereport('cd_findreaction',errcode,cmessage)
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

INTEGER :: i,ntracer,mspecies                          ! loop variables

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
  DO i=1,n_chemdiags
    IF (asad_chemdiags(i)%diag_type == cdems) THEN
      DO ntracer=1,jpctr
        DO mspecies=1,n_emissions
           ! work out which molmass to take!
          IF ((advt(ntracer) == em_spec(mspecies)) .AND.                       &
             (asad_chemdiags(i)%species==em_spec(mspecies))) THEN
            asad_chemdiags(i)%molecular_mass = molmass(mspecies)
            asad_chemdiags(i)%c_vmr_to_mmr = c_species(ntracer)
            L_cd_emitted(ntracer) = .TRUE.
            IF (asad_chemdiags(i)%species == 'NO        ')                     &
                       asad_chemdiags(i)%molecular_mass =                      &
                       asad_chemdiags(i)%molecular_mass*c_no/c_no2
          END IF
        END DO        ! mspecies,n_emissions
        DO mspecies=1,n_boundary_vals
           ! work out which molmass to take!, can overwrite with
           ! lbc_molmass, since emmissions_ctl does this anyway
          IF ((advt(ntracer) == lbc_spec(mspecies)) .AND.                      &
              (asad_chemdiags(i)%species == lbc_spec(mspecies)))               &
                     THEN
            asad_chemdiags(i)%molecular_mass=lbc_molmass(mspecies)
            asad_chemdiags(i)%c_vmr_to_mmr = c_species(ntracer)
            L_cd_emitted(ntracer) = .TRUE.
            IF (asad_chemdiags(i)%species == 'NO        ')                     &
                 asad_chemdiags(i)%molecular_mass =                            &
                 asad_chemdiags(i)%molecular_mass*c_no/c_no2
          END IF
        END DO        ! mspecies,n_boundary_vals

        IF (L_cd_emitted(ntracer)) THEN
          IF (advt(ntracer) == asad_chemdiags(i)%species)                      &
                    asad_chemdiags(i)%location = ntracer
        ELSE IF ((.NOT. L_cd_emitted(ntracer)) .AND.                           &
           (advt(ntracer) == asad_chemdiags(i)%species)) THEN
          cmessage='SPECIES '// asad_chemdiags(i)%species                      &
                //' NOT EMITTED'
          errcode=asad_chemdiags(i)%stash_number
          WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',           &
                errcode,' PE: ',mype
          CALL umPrint(umMessage,src='asad_chem_flux_diags')
          CALL ereport('ASAD_EMISSIONS_DIGNOSTICS',                            &
                errcode,cmessage)
        END IF
      END DO           ! ntracer,jpctr
    END IF              ! %diag_type == cdems
  END DO                 ! i,n_chemdiags
  firstcall=.FALSE.
END IF                    ! first

DO i=1,n_chemdiags
  IF (asad_chemdiags(i)%diag_type == cdems) THEN
    ! emissions are 2D and converted to (kg/m2/s)
    asad_chemdiags(i)%throughput(:,:,1) =                                      &
        em_field(:,:,asad_chemdiags(i)%location)*surf_area(:,:)*               &
        (1000.0/asad_chemdiags(i)%molecular_mass)
  END IF                 ! %diag_type
END DO                    ! i,n_chemdiags


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
  DO i=1,n_chemdiags
    IF (asad_chemdiags(i)%diag_type == cdems .AND.                             &
        asad_chemdiags(i)%rxn_type == emission_type .AND.                      &
        asad_chemdiags(i)%species == advt(i_spec_emiss)) THEN
      asad_chemdiags(i)%location = i_spec_emiss
      asad_chemdiags(i)%molecular_mass = em_molmass
      asad_chemdiags(i)%c_vmr_to_mmr   = c_species(i_spec_emiss)
      spec_emitted = .TRUE.
    END IF         ! diag_type == cdems &etc
  END DO
  IF (.NOT. spec_emitted) THEN
    cmessage='SPECIES '//asad_chemdiags(i)%species                             &
         //' NOT EMITTED IN 3D'
    errcode=asad_chemdiags(i)%stash_number
    WRITE(umMessage,'(2A,I0,A,I0)') cmessage, ' Error code: ',errcode,         &
         ' PE: ',mype
    CALL umPrint(umMessage,src='asad_chem_flux_diags')
    CALL ereport('ASAD_3D_EMISSIONS_DIGNOSTICS',                               &
         errcode,cmessage)
  END IF
  firstcall = .FALSE.
END IF                    ! firstcall

! need to do some conversion because of the way that things are done in
! UKCA for 3D emissions
em_field_3d(:,:,:) = -1.0 ! i.e. if values are -ve, have a problem
! sanity check, just in case
!em_field_3d(:,:,:) = em_field_3d_in(:,:,:)
DO i=1,n_chemdiags
  SELECT CASE (asad_chemdiags(i)%diag_type)
  CASE (cdems)
    IF ((asad_chemdiags(i)%rxn_type == aircraft_emissions)                     &
         .AND. (emission_type == aircraft_emissions)) THEN
       ! these are in kg(NO2)/m2/s
       ! - need to multiply by the surface area and / air mass
      DO klevel=1,model_levels
        DO j=1,rows
          DO l=1,row_length
            em_field_3d(l,j,klevel) =                                          &
             (em_field_3d_in(l,j,klevel)/mass(l,j,klevel))*                    &
              surf_area(l,j)
          END DO ! l
        END DO ! j
      END DO ! klevel
    ELSE IF ((asad_chemdiags(i)%rxn_type ==                                    &
         lightning_emissions) .AND. (emission_type ==                          &
         lightning_emissions)) THEN
       ! do nothing - in kg(NO)/kg(air)/gridcell/s
      em_field_3d(:,:,:) = em_field_3d_in(:,:,:)
    ELSE IF ((asad_chemdiags(i)%rxn_type ==                                    &
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
DO i=1,n_chemdiags
  IF (asad_chemdiags(i)%diag_type == cdems) THEN
    IF ((advt(i_spec_emiss) == asad_chemdiags(i)%species) .AND.                &
             (emission_type == asad_chemdiags(i)%rxn_type)) THEN
      asad_chemdiags(i)%throughput(:,:,:)= (em_field_3d(                       &
                     1:row_length,1:rows,1:model_levels)                       &
                     *mass(1:row_length,1:rows,1:model_levels))                &
                      *(1000.0/asad_chemdiags(i)%molecular_mass)
      IF (asad_chemdiags(i)%tropospheric_mask) THEN
        WHERE (.NOT. L_troposphere)
          asad_chemdiags(i)%throughput(:,:,:) = 0.0
        END WHERE
      END IF
    END IF              ! advt==%species
  END IF             ! %diag_type == cdems
END DO                    ! i,n_chemdiags

ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_3d_emissions_diagnostics

! #####################################################################
SUBROUTINE asad_tropospheric_mask(ierr)

USE ukca_tropopause , ONLY: L_troposphere
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
      WHERE (.NOT. L_Troposphere(:,:,:))
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
        WHERE (.NOT. L_troposphere)
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

INTEGER :: i                      ! counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_MASS_DIAGNOSTIC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr = -1

DO i=1,n_chemdiags
  IF (asad_chemdiags(i)%diag_type == cdmas) THEN
    asad_chemdiags(i)%throughput(:,:,:) = mass(:,:,:)
    IF (asad_chemdiags(i)%tropospheric_mask) THEN
      WHERE (.NOT. L_Troposphere)
        asad_chemdiags(i)%throughput(:,:,:) = 0.0
      END WHERE
    END IF
  END IF
END DO                    ! i,n_chemdiags

ierr = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_mass_diagnostic
! #####################################################################
SUBROUTINE asad_psc_diagnostic(row_length, rows, model_levels,                 &
                               fpsc1_full, fpsc2_full, ix, jy, klevel, ierr)

USE ukca_config_specification_mod, ONLY: ukca_config
IMPLICIT NONE


INTEGER, INTENT(IN)  :: row_length             ! array dimension
INTEGER, INTENT(IN)  :: rows                   ! array dimension
INTEGER, INTENT(IN)  :: model_levels           ! number of model levels
INTEGER, INTENT(IN)  :: ix                     ! i counter
INTEGER, INTENT(IN)  :: jy                     ! j counter
INTEGER, INTENT(IN)  :: klevel                 ! level number
INTEGER, INTENT(OUT) :: ierr                   ! error code
REAL, INTENT(IN) :: fpsc1_full(model_levels)
REAL, INTENT(IN) :: fpsc2_full(model_levels)

INTEGER :: i                                   ! counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_PSC_DIAGNOSTIC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
ierr = -1

DO i=1,n_chemdiags

  IF (ukca_config%l_ukca_asad_columns) THEN
     ! in this case ASAD is being called column-by column so we need
     ! to fill %throughput by iterating over X & Y
    SELECT CASE (asad_chemdiags(i)%diag_type)
    CASE (cdpsc)
      SELECT CASE (asad_chemdiags(i)%rxn_type)
      CASE (cdpsc_typ1)
        asad_chemdiags(i)%throughput(ix,jy,:) = fpsc1_full(:)
        IF (asad_chemdiags(i)%tropospheric_mask) THEN
          WHERE (.NOT. L_Troposphere(ix,jy,:))
            asad_chemdiags(i)%throughput(ix,jy,:) = 0.0
          END WHERE
        END IF
      CASE (cdpsc_typ2)
        asad_chemdiags(i)%throughput(ix,jy,:) = fpsc2_full(:)
        IF (asad_chemdiags(i)%tropospheric_mask) THEN
          WHERE (.NOT. L_Troposphere(ix,jy,:))
            asad_chemdiags(i)%throughput(ix,jy,:) = 0.0
          END WHERE
        END IF
      END SELECT
    END SELECT
  ELSE
     ! in this case ASAD is being called horizontally so we need
     ! to fill %throughput by iterating over Z
    SELECT CASE (asad_chemdiags(i)%diag_type)
    CASE (cdpsc)
      SELECT CASE (asad_chemdiags(i)%rxn_type)
      CASE (cdpsc_typ1)
        asad_chemdiags(i)%throughput(:,:,klevel) =                             &
             RESHAPE(fpsc1(:),[row_length,rows])
        IF (asad_chemdiags(i)%tropospheric_mask) THEN
          WHERE (.NOT. L_Troposphere(:,:,klevel))
            asad_chemdiags(i)%throughput(:,:,klevel) = 0.0
          END WHERE
        END IF
      CASE (cdpsc_typ2)
        asad_chemdiags(i)%throughput(:,:,klevel) =                             &
             RESHAPE(fpsc2(:),[row_length,rows])
        IF (asad_chemdiags(i)%tropospheric_mask) THEN
          WHERE (.NOT. L_Troposphere(:,:,klevel))
            asad_chemdiags(i)%throughput(:,:,klevel) = 0.0
          END WHERE
        END IF
      END SELECT
    END SELECT
  END IF ! l_ukca_asad_columns
END DO

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
      WHERE (.NOT. L_Troposphere(:,:,:))
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
        WHERE (.NOT. L_Troposphere(:,:,1))
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


! #####################################################################
SUBROUTINE asad_flux_put_stash(                                                &
   stashwork_size,STASHwork,                                                   &
   row_length,rows,model_levels)

IMPLICIT NONE


INTEGER, INTENT(IN) :: stashwork_size
REAL, INTENT(IN OUT) :: STASHwork(1:stashwork_size)
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! Flux arrays
REAL, ALLOCATABLE :: upload_array_3D(:,:,:)
REAL, ALLOCATABLE :: upload_array_2D(:,:)

INTEGER       :: item                             ! stash item
INTEGER       :: section                          ! stash section
INTEGER       :: im_index                         ! atmos model index
INTEGER       :: i,j,k,l                          ! Counters
INTEGER       :: icode                            ! Error code
CHARACTER(LEN=errormessagelength) :: cmessage     ! Error return message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_FLUX_PUT_STASH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode = 0

im_index = 1

ALLOCATE(upload_array_2D(1:row_length,1:rows))
ALLOCATE(upload_array_3D(1:row_length,1:rows,1:model_levels))

! copy items into STASHwork array here, which is then uploaded with
! the call to STASH at the end of UKCA_MAIN.
DO l=1,n_stashsumdiag
  item = stash_handling(l)%stash_item
  section = stash_handling(l)%stash_section
  upload_array_2D(:,:)   = 0.0
  upload_array_3D(:,:,:) = 0.0

  ! Now check if can output diagnostic, since some of these
  !  will only be allowed on chemical timesteps
  IF (stash_handling(l)%len_dim3 == 1) THEN                   ! 2D field
    DO k=1,stash_handling(l)%number_of_fields
         ! will only be summing surface fields in this case
      upload_array_2D(:,:) = upload_array_2D(:,:) +                            &
              asad_chemdiags(                                                  &
              stash_handling(l)%chemdiags_location(k)                          &
              )%throughput(:,:,1)
    END DO
    IF (sf(item,section)) THEN

      CALL copydiag (stashwork(si(item,section,im_index):                      &
           si_last(item,section,im_index)),                                    &
           upload_array_2D(:,:), row_length,rows)
    END IF
  ELSE ! len_dim3 /= 1
    DO k=1,stash_handling(l)%number_of_fields
        ! may be summing both 3D and surface fields into the same 3D array
      DO j=1,asad_chemdiags(                                                   &
        stash_handling(l)%chemdiags_location(k))%num_levs
        upload_array_3D(:,:,j) = upload_array_3D(:,:,j) +                      &
                asad_chemdiags(                                                &
                stash_handling(l)%chemdiags_location(k)                        &
                )%throughput(:,:,j)
      END DO
    END DO

    IF (sf(item,section)) THEN
      CALL copydiag_3d (stashwork(si(item,section,im_index):                   &
            si_last(item,section,im_index)),                                   &
            upload_array_3D(:,:,:),                                            &
            row_length,rows,model_levels,                                      &
            stlist(:,stindex(1,item,section,im_index)),                        &
            len_stlist,                                                        &
            stash_levels,num_stash_levels+1)
    END IF
  END IF ! len_dim3
  IF (icode >  0) THEN
    WRITE(cmessage,'(A,A,5(I6,1X))') 'ASAD_FLUX_PUT_STASH: ',                  &
          'ERROR WITH copydiag ',                                              &
          item,section,im_index,atmos_im,icode
    WRITE(umMessage,'(A)') cmessage
    CALL umPrint(umMessage,src='asad_chem_flux_diags')
    CALL umPrintFlush()
    icode=(1000*section) + item
    CALL ereport('ASAD_FLUX_PUT_STASH',icode,cmessage)
  END IF
END DO       ! l,n_stashsumdiag

! Deallocate arrays, if are able to
DO i=1,n_chemdiags
  IF (asad_chemdiags(i)%can_deallocate)                                        &
       DEALLOCATE(asad_chemdiags(i)%throughput)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_flux_put_stash

END MODULE asad_chem_flux_diags
