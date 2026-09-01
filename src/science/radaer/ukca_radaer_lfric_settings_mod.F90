! -----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! -----------------------------------------------------------------------------
!
! Purpose: This could allow RADAER to have different GLOMAP settings
!          to UKCA (e.g. with option dust_and_clim)
!
! Procedure:
!   1) Allocate structure ukca_radaer_lfric
!   2) Obtain number of modes and components for specific GLOMAP setting
!   3) Populate structure ukca_radaer_lfric
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
MODULE ukca_radaer_lfric_settings_mod

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName =                             &
                                               'UKCA_RADAER_LFRIC_SETTINGS_MOD'

CONTAINS

SUBROUTINE ukca_radaer_lfric_settings( l_ukca_radaer_sustrat,                  &
                                       glomap_variables_radaer,                &
                                       ukca_radaer_lfric,                      &
                                       n_loc_mode,                             &
                                       n_loc_cpnt )

USE ereport_mod,                       ONLY:                                   &
    ereport

USE errormessagelength_mod,            ONLY:                                   &
    errormessagelength

USE ukca_mode_setup,                   ONLY:                                   &
    glomap_variables_type,                                                     &
    nmodes,                                                                    &
    mode_names,                                                                &
    ncp_max,                                                                   &
    mode_sup_insol,                                                            &
    ip_ukca_mode_nucleation,                                                   &
    ip_ukca_mode_aitken,                                                       &
    ip_ukca_mode_accum,                                                        &
    ip_ukca_mode_coarse,                                                       &
    ip_ukca_mode_supercoarse,                                                  &
    cp_su,                                                                     &
    cp_bc,                                                                     &
    cp_oc,                                                                     &
    cp_cl,                                                                     &
    cp_du,                                                                     &
    cp_so,                                                                     &
    cp_no3,                                                                    &
    cp_nh4,                                                                    &
    cp_nn,                                                                     &
    cp_mp

USE ukca_radaer_lfric_struct_mod,      ONLY:                                   &
    ukca_radaer_lfric_struct,                                                  &
    allocate_radaer_lfric_struct

USE parkind1,                          ONLY:                                   &
    jpim,                                                                      &
    jprb

USE yomhook,                           ONLY:                                   &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments
! =========

! Switch: Use sulphuric acid optical properties in the stratosphere,
! instead of ammonium sulphate.
LOGICAL, INTENT(IN) :: l_ukca_radaer_sustrat

TYPE(glomap_variables_type), TARGET, INTENT(IN OUT) :: glomap_variables_radaer
TYPE (ukca_radaer_lfric_struct),     INTENT(IN OUT) :: ukca_radaer_lfric

! Counters for number of modes and components
INTEGER, INTENT(IN OUT) :: n_loc_mode
INTEGER, INTENT(IN OUT) :: n_loc_cpnt

! Local variables
! ===============

LOGICAL,            POINTER :: component(:,:)
CHARACTER(LEN=7),   POINTER :: component_names(:)
REAL,               POINTER :: ddplim0(:)
REAL,               POINTER :: ddplim1(:)
LOGICAL,            POINTER :: mode(:)
INTEGER,            POINTER :: modesol(:)
INTEGER,            POINTER :: ncp
REAL,               POINTER :: rhocomp(:)
REAL,               POINTER :: sigmag(:)

! Loop counters
INTEGER :: i, j

! In-loop copy of mode names
CHARACTER(LEN=7) :: this_name

! In-loop mode type
INTEGER :: this_type

! Error message
INTEGER :: ierr
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*),   PARAMETER :: RoutineName='UKCA_RADAER_LFRIC_SETTINGS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component       => glomap_variables_radaer%component
component_names => glomap_variables_radaer%component_names
ddplim0         => glomap_variables_radaer%ddplim0
ddplim1         => glomap_variables_radaer%ddplim1
mode            => glomap_variables_radaer%mode
modesol         => glomap_variables_radaer%modesol
ncp             => glomap_variables_radaer%ncp
rhocomp         => glomap_variables_radaer%rhocomp
sigmag          => glomap_variables_radaer%sigmag

! Allocate elements that depend on ncp or nmodes
CALL allocate_radaer_lfric_struct( ukca_radaer_lfric, glomap_variables_radaer )

! Loop over all modes and components per mode,
! and only retain included modes and components.
ukca_radaer_lfric%idx_cpnt_mode(:,:)=-1
DO i = 1, nmodes

  IF (mode(i)) THEN

    this_name = mode_names(i)

    ! Get the mode type. Since there is no direct information,
    ! it is obtained from the mode names.

    SELECT CASE ( this_name(1:3) )

    CASE ('Nuc')
      this_type = ip_ukca_mode_nucleation

    CASE ('Ait')
      this_type = ip_ukca_mode_aitken

    CASE ('Acc')
      this_type = ip_ukca_mode_accum

    CASE ('Cor')
      this_type = ip_ukca_mode_coarse

    CASE ('Sup')
      this_type = ip_ukca_mode_supercoarse

    CASE DEFAULT
      ierr = 1
      cmessage = 'Unexpected mode name.' // this_name
      CALL ereport(RoutineName,ierr,cmessage)

    END SELECT

    ! Interaction of nucleation modes with radiation is neglected.
    IF ( this_type /= ip_ukca_mode_nucleation ) THEN

      ! Increase mode counter by one
      n_loc_mode = n_loc_mode + 1

      ukca_radaer_lfric%i_mode_type(n_loc_mode) = this_type
      ukca_radaer_lfric%l_soluble(n_loc_mode) = modesol(i) == 1
      ukca_radaer_lfric%d0low(n_loc_mode) = ddplim0(i)
      ukca_radaer_lfric%d0up(n_loc_mode) = ddplim1(i)
      ukca_radaer_lfric%sigma(n_loc_mode) = sigmag(i)
      ukca_radaer_lfric%n_cpnt_in_mode(n_loc_mode) = 0

      ! Loop on components within that mode.
      DO j = 1, ncp
        IF ( component(i, j) ) THEN

          ! Increase component counter by one
          n_loc_cpnt = n_loc_cpnt + 1

          ! Update the number of components in that mode and
          ! retain the array index of the current component.

          ukca_radaer_lfric%n_cpnt_in_mode(n_loc_mode) =                       &
                               ukca_radaer_lfric%n_cpnt_in_mode(n_loc_mode) + 1

          ukca_radaer_lfric%i_cpnt_index(                                      &
                    ukca_radaer_lfric%n_cpnt_in_mode(n_loc_mode),n_loc_mode) = &
                                          n_loc_cpnt

          ukca_radaer_lfric%idx_cpnt_mode(i,j) = n_loc_cpnt

          ukca_radaer_lfric%density(n_loc_cpnt) = rhocomp(j)

          ukca_radaer_lfric%i_mode(n_loc_cpnt) = n_loc_mode

          ! Get the component type. Since there is no direct
          ! information, it is obtained from the component
          ! names.

          SELECT CASE ( component_names(j) )

          CASE ('h2so4  ')
            ukca_radaer_lfric%i_cpnt_type(n_loc_cpnt) = cp_su

          CASE ('bcarbon')
            ukca_radaer_lfric%i_cpnt_type(n_loc_cpnt) = cp_bc

          CASE ('ocarbon')
            ukca_radaer_lfric%i_cpnt_type(n_loc_cpnt) = cp_oc

          CASE ('nacl   ')
            ukca_radaer_lfric%i_cpnt_type(n_loc_cpnt) = cp_cl

          CASE ('dust   ')
            ukca_radaer_lfric%i_cpnt_type(n_loc_cpnt) = cp_du

          CASE ('sec_org')
            ukca_radaer_lfric%i_cpnt_type(n_loc_cpnt) = cp_so

          CASE ('no3    ')
            ukca_radaer_lfric%i_cpnt_type(n_loc_cpnt) = cp_no3

          CASE ('nano3  ')
            ukca_radaer_lfric%i_cpnt_type(n_loc_cpnt) = cp_nn

          CASE ('nh4    ')
            ukca_radaer_lfric%i_cpnt_type(n_loc_cpnt) = cp_nh4

          CASE ('mp     ')
            ukca_radaer_lfric%i_cpnt_type(n_loc_cpnt) = cp_mp

          CASE DEFAULT
            ierr = 2
            cmessage = 'Unexpected component name: ' // component_names(j)
            CALL ereport(RoutineName,ierr,cmessage)

          END SELECT

        END IF
      END DO ! j

    END IF

  END IF
END DO ! i

ukca_radaer_lfric%n_mode = n_loc_mode
ukca_radaer_lfric%n_cpnt = n_loc_cpnt
ukca_radaer_lfric%l_sustrat = l_ukca_radaer_sustrat

! Ensure that the coarse narrow LUTs are used for the coarse insoluble mode
! if the super-coarse insoluble mode is selected
ukca_radaer_lfric%l_cornarrow_ins = mode(mode_sup_insol)

IF (ukca_radaer_lfric%n_mode == 0 .OR. ukca_radaer_lfric%n_cpnt == 0) THEN
  ierr = 3
  cmessage = 'Setup includes no UKCA aerosols.'
  CALL ereport(RoutineName,ierr,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE ukca_radaer_lfric_settings

END MODULE ukca_radaer_lfric_settings_mod
