! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module to apply vertical and temporal profiles to emissions, and do
!   some conversions depending on the units of the emission fields.
!
! Method:
!   Four main subroutines are contained in this module to correct the
!   emission fields:
!   1) vertical_emiss_factors: Get factors to spread emiss over model layers.
!   2) hourly_emiss_factors:   Set array with diurnal cycle to be applied.
!   3) daily_emiss_factors     Set array with day-to-day variability.
!   4) base_emiss_factors      Get correction factors based on units
!                              of emission fields.
!
! Part of the UKCA model, a community model supported by
! The Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------------
!
MODULE ukca_emiss_factors

USE ereport_mod,      ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,       ONLY: umPrint, umMessage
USE ukca_constants,   ONLY: m_c, m_c5h8, m_ch3oh, m_ch4, m_dms,                &
                            m_me2co, m_monoterp, m_n, m_nh3, m_no,             &
                            m_no2, m_s, m_so2
USE parkind1,         ONLY: jpim, jprb     ! DrHook
USE yomhook,          ONLY: lhook, dr_hook ! DrHook

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EMISS_FACTORS'

CONTAINS

!---------------------------------------------------------------------------
! Description:
!   Determine the factors needed to inject emissions in each model level.
!   This may be used to add emissions to all vertical levels or spread
!   them over a fixed number as specified in the netCDF attributes.
!
! Method:
!   Store vertical factors in the array vert_scaling_3d, which contains
!   one value for each grid cell. The sum of these factors for all cells
!   within a column is always one (except for 3D emisssions).
!---------------------------------------------------------------------------
SUBROUTINE vertical_emiss_factors (row_length, rows, model_levels,             &
               lowest_lev, highest_lev, vert_fact, field_name, vert_scaling_3d)

USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_fieldname_mod, ONLY: fldname_interf_z
USE ukca_environment_req_mod, ONLY: environ_field_available
USE ukca_environment_fields_mod, ONLY: interf_z

IMPLICIT NONE

! Subroutine arguments
INTEGER,            INTENT(IN) :: row_length   ! Model dimensions
INTEGER,            INTENT(IN) :: rows
INTEGER,            INTENT(IN) :: model_levels
INTEGER,            INTENT(IN) :: lowest_lev   ! Lowest and highest lev for
INTEGER,            INTENT(IN) :: highest_lev  ! point or high/single lev emiss

CHARACTER (LEN=*),  INTENT(IN) :: vert_fact    ! Surface, several model levs ...
CHARACTER (LEN=*),  INTENT(IN) :: field_name   ! Name of emiss field

REAL, INTENT(OUT) :: vert_scaling_3d (1:row_length, 1:rows, 1:model_levels)
! Emitted fraction for each grid cell (accounts for vertical profiles)

! Local variables
INTEGER     :: nlevs       ! Nr of model levels considered for
                           ! 'high_level' or 'single_level' emiss
INTEGER     :: nlayers_ref ! Number of layers in reference emiss profiles
INTEGER     :: ierror      ! Error code

! Altitude of the layer interfaces and vertical factors for average
! reference profiles used for specific source sectors:
REAL, ALLOCATABLE   :: interf_ref       (:)  ! dims --> 0:nlayers_ref
REAL, ALLOCATABLE   :: vert_scaling_ref (:)  ! dims --> 1:nlayers_ref

! Vertical thickness of model layers and total depth over a range of layers.
! Used for weighting vertical scaling factors by layer thickness
REAL :: zthick (1:row_length, 1:rows, 1:model_levels)
REAL :: total_depth (1:row_length, 1:rows)

CHARACTER (LEN=errormessagelength) :: cmessage  ! error message

! This logical is TRUE if we use avg ref vertic profs at specific
! altitudes for SNAP source sectors (typical for air quality).
! If other types of emission sectors (e.g. following UNECE-LRTAP
! or IPCC conventions) are needed in the future then this logical
! should be replaced by an integer variable which could adopt a
! different value (e.g. 1, 2, 3, ...) depending on the conventions.
LOGICAL :: l_ref_snap_prof

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='VERTICAL_EMISS_FACTORS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierror = 0

! Some of the profiles used below (e.g. those for SNAP source
! sectors, i.e. 'EMEP_modified_SNAP01' ... 'EMEP_modified_SNAP10',
! and 'Bieser_modified_SNAP01' ... 'Bieser_modified_SNAP10') consider
! the average percetage contribution of emissions over some fixed
! vertical layers. We use a logical variable to indicate that
! these cases need special treatment.
l_ref_snap_prof = ( INDEX (vert_fact, 'EMEP_modified')   > 0 .OR.              &
                    INDEX (vert_fact, 'Bieser_modified') > 0 )

! When emissions are given for independent SNAP source sectors we use
! average profiles within 7 vertical layers (note that profiles are
! applied at the mid-point of a layer and that one element more
! is needed for the interfaces).
IF ( l_ref_snap_prof ) THEN
  nlayers_ref = 7
  ALLOCATE (vert_scaling_ref (1:nlayers_ref))
  ALLOCATE (interf_ref       (0:nlayers_ref))

  ! Altitude of the interfaces in metres as reported by Mailler et al. (ACP,
  ! 2013). The vertical scaling is filled in within the CASE statement below
  interf_ref (0:nlayers_ref) =                                                 &
       [0.0, 20.0, 92.0, 184.0, 324.0, 522.0, 781.0, 1106.0]
END IF

! Ensure that the interf_z array is available as an environment field if needed
IF (l_ref_snap_prof .OR. vert_fact == 'high_level' .OR.                        &
    vert_fact == 'single_level') THEN
  IF (.NOT. ukca_config%l_support_ems_vertprof) THEN
    ierror = 1
    cmessage = 'Emission field ' // field_name // ' has vertical factor ' //   &
               vert_fact // ' that is not supported in this configuration'
    CALL ereport (ModuleName//':'//RoutineName, ierror, cmessage)
  END IF
  IF (.NOT. environ_field_available(fldname_interf_z)) THEN
    ierror = 1
    cmessage = 'Required field ' //fldname_interf_z//                          &
               ' is not available in environment'
    CALL ereport (ModuleName//':'//RoutineName, ierror, cmessage)
  END IF
END IF

!----------------------------------------------------------------------------
! Get a 3D vertical scaling factor for the given profile
SELECT CASE (TRIM (vert_fact))
  ! Surface emiss in lowest model lev unless something different specified
CASE ('surface', '')
  vert_scaling_3d (:,:,1)   = 1.0
  vert_scaling_3d (:,:,2:)  = 0.0

  ! 3D emissions over all model levels
CASE ('all_levels', '3D')
  vert_scaling_3d (:,:,:)   = 1.0

  ! Special vertical profile used in air quality simulations. This
  ! spreads surface emissions over some of the lowest model levels
  ! for a L38 vertical configuration of the UM: 80%, 15% and 5%
  ! of the total emission are spread over the 1st, 2nd and 3rd
  ! lowest model layers, respectively.
CASE ('step1')
  IF (model_levels == 38) THEN
    vert_scaling_3d (:,:,1)  = 0.80
    vert_scaling_3d (:,:,2)  = 0.15
    vert_scaling_3d (:,:,3)  = 0.05
    vert_scaling_3d (:,:,4:) = 0.00
  ELSE
    ierror = 1
  END IF

  ! Spread high_level emiss or any emiss specified as 'single_level'
  ! (including surface emissions) equally between the lowest
  ! and highest levels read from the NetCDF files. In order to have a uniform
  ! concentration (kg/m3) in height, the scaling factors must be weighted by
  ! the thickness of each model level, since the applied emission is total kg/m2
  ! for each model level. The scaling factors still total 1.0 in any given
  ! column. Although thickness varies by up to a factor of 2 horizontally,
  ! the weights themselves only very by a few percent between columns, since
  ! the proportional thicknesses are approximately constant.
  ! SPREAD is required to cast the 2D array total_depth to the same shape as
  ! the 3D slice of zthick.
CASE ('high_level', 'single_level')
  IF (highest_lev < 1 .OR. highest_lev > model_levels) THEN
    ierror = 2  ! Highest level out of range
  ELSE IF (lowest_lev < 1 .OR. lowest_lev > model_levels) THEN
    ierror = 3  ! Lowest level out of range
  ELSE
    nlevs = highest_lev - lowest_lev + 1
    IF (vert_fact == 'single_level' .AND. nlevs /= 1) THEN
      ierror = 4  ! Highest and lowest level mismatch for 'single_level' profile
    ELSE IF (nlevs <= 0 .OR. nlevs >= model_levels) THEN
      ierror = 5  ! Invalid number of levels
    ELSE
      zthick(:,:,:) = interf_z(:,:,1:model_levels)                             &
                    - interf_z(:,:,0:model_levels-1)
      total_depth(:,:) = interf_z(:,:,highest_lev) - interf_z(:,:,lowest_lev-1)

      vert_scaling_3d (:,:,:) = 0.0
      vert_scaling_3d (:,:,lowest_lev:highest_lev) =                           &
          zthick(:,:,lowest_lev:highest_lev) /                                 &
          SPREAD(total_depth(:,:), 3, nlevs)
    END IF
  END IF

  !--------------------------------------------------------------------
  !   Average vertical profiles for various SNAP source sectors.
  !   Similarly as done for the EMEP model, which considers 6
  !   vertical layers (e.g. Simpson et al., 2012, ACP; Bieser
  !   et al., 2011, Environ. Pollut.)
  !
  !   We add an additional 0-20 m layer for emissions close to surface
  !   as done by Mailler et al. (2013, ACP). Only UM configurations
  !   with a high enough vertical resolution close to surface will
  !   benefit from this extra layer. The altitudes of interfaces of
  !   the 7 resulting layers are stored in 'interf_ref' (see above).
  !
  !   Note that these are percentages and will be divided by 100 below.
  !   The array 'vert_scaling_3d' will also be filled later.
  !
CASE ('EMEP_modified_SNAP01')
  vert_scaling_ref (:) = [  0.0,  0.0,  0.0,  8.0, 46.0, 29.0, 17.0 ]

CASE ('EMEP_modified_SNAP02')
  vert_scaling_ref (:) = [ 11.0, 39.0, 50.0,  0.0,  0.0,  0.0,  0.0 ]

CASE ('EMEP_modified_SNAP03')
  vert_scaling_ref (:) = [  0.0,  0.0,  4.0, 19.0, 41.0, 30.0,  6.0 ]

CASE ('EMEP_modified_SNAP04')
  vert_scaling_ref (:) = [ 20.0, 70.0, 10.0,  0.0,  0.0,  0.0,  0.0 ]

CASE ('EMEP_modified_SNAP05')
  vert_scaling_ref (:) = [ 20.0, 70.0, 10.0,  0.0,  0.0,  0.0,  0.0 ]

CASE ('EMEP_modified_SNAP06')
  vert_scaling_ref (:) = [ 100.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 ]

CASE ('EMEP_modified_SNAP07')
  vert_scaling_ref (:) = [ 100.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 ]

CASE ('EMEP_modified_SNAP08')
  vert_scaling_ref (:) = [ 100.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 ]

CASE ('EMEP_modified_SNAP09')
  vert_scaling_ref (:) = [  2.0,  8.0, 15.0, 40.0, 35.0,  0.0,  0.0 ]

CASE ('EMEP_modified_SNAP10')
  vert_scaling_ref (:) = [ 100.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 ]

  !-----------------------------------------------------------------------
  !   Average vertical profiles for various SNAP source sectors.
  !   Similarly as done by Bieser et al. (2011, Environ. Pollut.)
  !   with the SMOKE model and including fugitive emissions.
  !
  !   As before, add additional 0-20 m layer for emissions close to surface
  !   as done by Mailler et al. (2013, ACP). These are percentages which will
  !   be divided by 100 below in the code.
  !
CASE ('Bieser_modified_SNAP01')
  vert_scaling_ref (:) = [  0.0,  0.0, 0.25, 51.0, 45.3, 3.25, 0.2 ]

CASE ('Bieser_modified_SNAP02')
  vert_scaling_ref (:) = [ 11.0, 89.0, 0.0,   0.0,  0.0,  0.0, 0.0 ]

CASE ('Bieser_modified_SNAP03')
  vert_scaling_ref (:) = [  0.0, 21.3, 75.4,  3.3,  0.0,  0.0, 0.0 ]

CASE ('Bieser_modified_SNAP04')
  vert_scaling_ref (:) = [ 20.0, 70.0, 7.0,   3.0,  0.0,  0.0, 0.0 ]

CASE ('Bieser_modified_SNAP05')
  vert_scaling_ref (:) = [ 20.0, 70.0, 7.0,   3.0,  0.0,  0.0, 0.0 ]

CASE ('Bieser_modified_SNAP06')
  vert_scaling_ref (:) = [ 100.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 ]

CASE ('Bieser_modified_SNAP07')
  vert_scaling_ref (:) = [ 100.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 ]

CASE ('Bieser_modified_SNAP08')
  vert_scaling_ref (:) = [ 100.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 ]

CASE ('Bieser_modified_SNAP09')
  vert_scaling_ref (:) = [   2.0,  8.0, 37.0,  51.0, 2.0,  0.0, 0.0 ]

CASE ('Bieser_modified_SNAP10')
  vert_scaling_ref (:) =  [ 100.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 ]

  ! Report error if something different specified
CASE DEFAULT
  ierror = 1

  ! The code might also need to be extended in the future to account
  ! for emissions which are BL height dependent
END SELECT

! Report error if an unexpected value of the attribute 'vert_scaling'
! was found in the NetCDF file
IF (ierror > 0) THEN
  SELECT CASE (ierror)
  CASE (1)
    cmessage = 'Unexpected vertical profile ('  // TRIM (vert_fact) //         &
               ') found for ' // TRIM (field_name)
  CASE (2)
    cmessage = 'Highest level specified is out of range'
  CASE (3)
    cmessage = 'Lowest level specified is out of range'
  CASE (4)
    cmessage = 'Non-matching highest and lowest levels specified for ' //      &
               '''single_level'' vertical profile'
  CASE (5)
    WRITE (umMessage,'(A,1X,A,1X,A,1X,I3)')                                    &
      'Number of vertical levels used to spread emissions of',                 &
      TRIM (field_name), 'is', nlevs
    CALL umPrint(umMessage,src='vertical_emiss_factors')
    cmessage = 'Number of vertical levels must be greater than zero '  //      &
               'and lower than the total number of model levels for '  //      &
               TRIM (vert_fact) // ' emissions'
  END SELECT

  CALL ereport (RoutineName, ierror, cmessage)
END IF

! When using avg. profiles given for specific altitudes we still need to
! calculate 'vert_scaling_3d' for all the model grid cells. This is done
! with the call to REGRID_VERTIC_PROFS
IF (l_ref_snap_prof) THEN
  vert_scaling_ref (:)     = vert_scaling_ref (:) / 100.0 ! now between 0 & 1
  vert_scaling_3d  (:,:,:) = 0

  CALL regrid_vertic_profs (                                                   &
    nlayers_ref, row_length, rows, model_levels,                               &
    interf_ref       (0:nlayers_ref),                                          &
    vert_scaling_ref (1:nlayers_ref),                                          &
    interf_z         (1:row_length, 1:rows, 0:model_levels),                   &
    field_name, vert_fact,                                                     &
    vert_scaling_3d  (1:row_length, 1:rows, 1:model_levels) )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE vertical_emiss_factors


!-----------------------------------------------------------------------
! Description:
!   Set up array with diurnal cycle to be applied to the emission fields
!
! Method:
!   Store values in the array hourly_scaling, with values for
!   hours 0, 1, 2, ..., 22, 23.
!-----------------------------------------------------------------------
SUBROUTINE hourly_emiss_factors (hourly_fact, field_name, hourly_scaling)

IMPLICIT NONE

! Subroutine arguments
CHARACTER (LEN=*), INTENT(IN)  :: hourly_fact     ! traffic, none,
                                                  ! diurnal_isopems
CHARACTER (LEN=*), INTENT(IN)  :: field_name      ! Name of emiss field

REAL,              INTENT(OUT) :: hourly_scaling (0:23)
                                   ! Factors for hour = 0, ... 23

! Local variables
INTEGER                            :: ierror
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HOURLY_EMISS_FACTORS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierror = 0

SELECT CASE (TRIM (hourly_fact))

  ! No hourly factors applied unless specified
CASE ('none', '')
  hourly_scaling =                                                             &
     [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                          &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                         &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]

  ! This is used for air quality studies over the UK
CASE ('traffic_uk')
  hourly_scaling =                                                             &
     [0.19, 0.09, 0.06, 0.05, 0.09, 0.22, 0.86, 1.84,                          &
       1.86, 1.41, 1.24, 1.20, 1.32, 1.44, 1.45, 1.59,                         &
       2.03, 2.08, 1.51, 1.06, 0.74, 0.62, 0.61, 0.44]

  ! -----------------------------------------------------------------------
  ! Hourly factors of emissions for Europe. Calculated by TNO for the
  ! MACC project.
  !
  ! Source:
  ! Hugo Denier van der Gon, Carlijn Hendriks, Jeroen Kuenen, Arjo Segers,
  ! Antoon Visschedijk: Description of current temporal emission patterns
  ! and sensitivity of predicted AQ for temporal emission patterns.
  ! EU FP7 MACC deliverable report D_D-EMIS_1.3, TNO report, Dec 2011.
  ! https://gmes-atmosphere.eu/documents/deliverables/d-emis/
  !         MACC_TNO_del_1_3_v2.pdf
  !
  ! The following hourly factors have been taken from Table on pg 22 of
  ! the report, which only includes anthropogenic emiss. An additional
  ! flat profile has been added here for SNAP 11 (emissions from natural
  ! and biogenic sources).
  !
  ! Profiles for SNAP 7 (road transport) are exactly like the 'traffic_uk'
  ! profile introduced above for air quality simulations over the UK.
  !
CASE ('TNO_MACC_EU_SNAP01')
  hourly_scaling =                                                             &
     [0.79, 0.72, 0.72, 0.71, 0.74, 0.80, 0.92, 1.08,                          &
       1.19, 1.22, 1.21, 1.21, 1.17, 1.15, 1.14, 1.13,                         &
       1.10, 1.07, 1.04, 1.02, 1.02, 1.01, 0.96, 0.88]

CASE ('TNO_MACC_EU_SNAP02')
  hourly_scaling =                                                             &
     [0.38, 0.36, 0.36, 0.36, 0.37, 0.50, 1.19, 1.53,                          &
       1.57, 1.56, 1.35, 1.16, 1.07, 1.06, 1.00, 0.98,                         &
       0.99, 1.12, 1.41, 1.52, 1.39, 1.35, 1.00, 0.42]

CASE ('TNO_MACC_EU_SNAP03')
  hourly_scaling =                                                             &
     [0.75, 0.75, 0.78, 0.82, 0.88, 0.95, 1.02, 1.09,                          &
       1.16, 1.22, 1.28, 1.30, 1.22, 1.24, 1.25, 1.16,                         &
       1.08, 1.01, 0.95, 0.90, 0.85, 0.81, 0.78, 0.75]

CASE ('TNO_MACC_EU_SNAP04')
  hourly_scaling =                                                             &
     [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                          &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                         &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]

CASE ('TNO_MACC_EU_SNAP05')
  hourly_scaling =                                                             &
     [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                          &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                         &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]

CASE ('TNO_MACC_EU_SNAP06')
  hourly_scaling =                                                             &
     [0.50, 0.35, 0.20, 0.10, 0.10, 0.20, 0.75, 1.25,                          &
       1.40, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,                         &
       1.50, 1.40, 1.25, 1.10, 1.00, 0.90, 0.80, 0.70]

CASE ('TNO_MACC_EU_SNAP07')
  hourly_scaling =                                                             &
     [0.19, 0.09, 0.06, 0.05, 0.09, 0.22, 0.86, 1.84,                          &
       1.86, 1.41, 1.24, 1.20, 1.32, 1.44, 1.45, 1.59,                         &
       2.03, 2.08, 1.51, 1.06, 0.74, 0.62, 0.61, 0.44]

CASE ('TNO_MACC_EU_SNAP08')
  hourly_scaling =                                                             &
     [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                          &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                         &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]

CASE ('TNO_MACC_EU_SNAP09')
  hourly_scaling =                                                             &
     [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                          &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                         &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]

CASE ('TNO_MACC_EU_SNAP10')
  hourly_scaling =                                                             &
     [0.60, 0.60, 0.60, 0.60, 0.60, 0.65, 0.75, 0.90,                          &
       1.10, 1.35, 1.45, 1.60, 1.65, 1.75, 1.70, 1.55,                         &
       1.35, 1.10, 0.90, 0.75, 0.65, 0.60, 0.60, 0.60]

CASE ('TNO_MACC_EU_SNAP11')
  hourly_scaling =                                                             &
     [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                          &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                         &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]

  ! ----------------------------------------------------------------
  ! The routine UKCA_DIURNAL_ISOP_EMS applies a diurnal cycle
  ! to isoprene emissions. Therefore a flat cycle is used here
  ! just in case the attribute 'diurnal_isopems' is specified
  ! in the NetCDF files.
CASE ('diurnal_isopems')
  hourly_scaling =                                                             &
     [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                          &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,                         &
       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]

  ! If something different then report error, just in case
  ! any other biogenic emissions with some diurnal cycle
  ! are specified in the NetCDF files.
CASE DEFAULT
  ierror = 1
END SELECT

! Call ereport if attribute 'hourly_scaling' of the NetCDF files
! contains information on a diurnal cycle not expected
IF (ierror > 0) THEN
  cmessage = 'Unexpected diurnal cycle ('  // TRIM (hourly_fact) //            &
             ') found for ' // TRIM (field_name)
  CALL ereport ('HOURLY_EMISS_FACTORS', ierror, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE hourly_emiss_factors


!-----------------------------------------------------------------------
! Description:
!   Apply weekly cycle to the emission fields.
!
! Method:
!   Store factos in the array daily_scaling, with values for
!   Sunday (1), Monday (2), ..., Friday (6), and Saturday (7).
!-----------------------------------------------------------------------
SUBROUTINE daily_emiss_factors (daily_fact, field_name, daily_scaling)

IMPLICIT NONE

! Subroutine arguments
CHARACTER (LEN=*), INTENT(IN)  :: daily_fact       ! traffic, none, ...
CHARACTER (LEN=*), INTENT(IN)  :: field_name       ! Name of emiss field

REAL,              INTENT(OUT) :: daily_scaling(7) ! Daily factors for day of
                                                   ! week: Sun, Mon, ..., Sat

! Local variables
INTEGER                            :: ierror
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DAILY_EMISS_FACTORS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierror = 0

SELECT CASE (TRIM (daily_fact))

  ! No daily factors applied unless specified
CASE ('none', '')
  daily_scaling = [1.000000, 1.000000, 1.000000, 1.000000,                     &
                    1.000000, 1.000000, 1.000000]

  ! Factors applied depending on day of the week (Sun, Mon, ..., Sat)
  ! This is used for air quality studies over the UK.
CASE ('traffic_uk')
  daily_scaling = [0.653264, 1.063470, 1.095050, 1.126510,                     &
                    1.126810, 1.200730, 0.734164]

  ! ----------------------------------------------------------------
  ! Daily factors of emissions derived for Great Britain, based on
  ! calculations by TNO for the MACC project (Denier van der Gon et
  ! al., 2011).
  !
  ! These profiles are generated for Great Britain, by averaging
  ! daily variation of VOC, CO and NOx emissions independently
  ! for each SNAP sector. They are not exactly like those given
  ! in the report, which were calculated for Europe as a whole.
  !
  ! A flat profile has been introduced for SNAP 11, as done above
  ! for the analogous daily factors.
  !
  ! The profile for SNAP 7 (road transport) is flatter than the
  ! 'traffic_uk' introduced above for air quality simulations
  ! over the UK.
  !
CASE ('TNO_MACC_GB_SNAP01')
  daily_scaling = [0.794600, 1.023600, 1.053500, 1.089400,                     &
                    1.093600, 1.065500, 0.879800]

CASE ('TNO_MACC_GB_SNAP02')
  daily_scaling = [0.941900, 0.992100, 1.028200, 1.080600,                     &
                    1.028800, 0.968800, 0.959600]

CASE ('TNO_MACC_GB_SNAP03')
  daily_scaling = [0.802200, 1.050800, 1.088800, 1.094500,                     &
                    1.084800, 1.075500, 0.803400]

CASE ('TNO_MACC_GB_SNAP04')
  daily_scaling = [0.803708, 1.058310, 1.086110, 1.087611,                     &
                    1.081810, 1.079810, 0.802641]

CASE ('TNO_MACC_GB_SNAP05')
  daily_scaling = [1.014486, 1.015886, 1.016985, 1.015885,                     &
                    0.962186, 0.962186, 1.012386]

CASE ('TNO_MACC_GB_SNAP06')
  daily_scaling = [0.267067, 1.209367, 1.210600, 1.212033,                     &
                    1.211733, 1.210967, 0.678233]

CASE ('TNO_MACC_GB_SNAP07')
  daily_scaling = [0.886409, 1.040343, 1.017710, 1.033410,                     &
                    1.037143, 1.105577, 0.879408]

CASE ('TNO_MACC_GB_SNAP08')
  daily_scaling = [0.997385, 1.000486, 1.002086, 1.003686,                     &
                    1.002586, 1.001086, 0.992685]

CASE ('TNO_MACC_GB_SNAP09')
  daily_scaling = [1.000000, 1.000000, 1.000000, 1.000000,                     &
                    1.000000, 1.000000, 1.000000]

CASE ('TNO_MACC_GB_SNAP10')
  daily_scaling = [1.007833, 0.999600, 0.991467, 0.991500,                     &
                    1.003967, 1.012100, 0.993533]

CASE ('TNO_MACC_GB_SNAP11')
  daily_scaling = [1.000000, 1.000000, 1.000000, 1.000000,                     &
                    1.000000, 1.000000, 1.000000]

  ! -----------------------------------------------------------
  ! If any other information then report error, just in case
  ! the user has specified a type of weekly cycle that is not
  ! accounted for in the present version of the code
CASE DEFAULT
  ierror = 1
END SELECT

! Call ereport if attribute 'daily_scaling' of the NetCDF files
! contains information on an unexpected weekly cycle
IF (ierror > 0) THEN
  cmessage = 'Unexpected weekly cycle ('  // TRIM (daily_fact) //              &
             ') found for ' // TRIM (field_name)
  CALL ereport ('DAILY_EMISS_FACTORS', ierror, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE daily_emiss_factors

!---------------------------------------------------------------------
! Description
!   Identify the necessary scaling factors needed to correct the
!   emission values so that they are consistent with the units used
!   for reporting emissions ('kg m-2 s-1' in most cases, following
!   NetCDF conventions).
!
! Method:
!   1) Check if units are 'kg m-2 s-1' as would be expected from NetCDF
!      files that comply with CF conventions. Different units are only
!      allowed for some online emissions (e.g. NOx from lightning) which
!      are calculated as totals per grid box; in that particular case the
!      logical 'gridbox_emiss' is set to TRUE to indicate that further
!      conversions need to be done in ukca_new_emiss_ctl.
!
!   2) Check if the standard_name attribute is present in the emiss field;
!      if so search for the substrings 'expressed_as_nitrogen',
!      'expressed_as_carbon' and 'expressed_as_sulfur' in it. If
!      standard_name is empty then look for the following substrings in
!      the long_name attribute: 'expressed as nitrogen',
!      'expressed as carbon' and 'expressed as sulfur'.
!      Note that standard_name contains underscores while long_name contains
!      spaces, and that IUPAC recommendations should be used (e.g. 'sulfur'
!      rather than 'sulphur').
!
!   3) If any of the mentioned substrings is found then call
!      GET_BASE_SCALING to calculate the appropriate scaling factors,
!      which will then be passed to UKCA_NEW_EMISS_CTL where the final
!      conversions will be applied.
!---------------------------------------------------------------------

SUBROUTINE base_emiss_factors (l_omit_std_name, tracer_name, units,            &
                               field_name, standard_name, long_name,           &
                               base_scaling, gridbox_emiss)

USE ukca_emiss_struct_mod, ONLY: maxlen_emiss_units

IMPLICIT NONE

! Subroutine arguments
LOGICAL,             INTENT (IN)   :: l_omit_std_name
                                        ! True to indicate that both standard
                                        ! and long names can be omitted
                                        ! (implies no need for name-based
                                        ! unit conversion)

CHARACTER (LEN=*),   INTENT (IN)   :: tracer_name    ! Name of emitted tracer
CHARACTER (LEN=*),   INTENT (IN)   :: units          ! Units of emiss field
CHARACTER (LEN=*),   INTENT (IN)   :: field_name     ! Name of emiss field
CHARACTER (LEN=*),   INTENT (IN)   :: standard_name
CHARACTER (LEN=*),   INTENT (IN)   :: long_name

REAL,                INTENT(IN OUT) :: base_scaling  ! Factor to apply to emiss
                                                     ! field, depending on units
LOGICAL,             INTENT(OUT)   :: gridbox_emiss  ! True if emiss are given
                                                     ! as kg gridbox-1 s-1

! Local variables
INTEGER                       :: ierror
CHARACTER(errormessagelength) :: cmessage          ! Error return message

LOGICAL             :: l_special_units  ! To indicate that special units are
                                        ! found in the standard_/long_name
CHARACTER (LEN=maxlen_emiss_units) :: special_units

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BASE_EMISS_FACTORS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierror = 0

! First make sure that emissions are given in kg m-2 s-1 , which are the only
! units accepted for NetCDF emission files (to comply with CF conventions).
IF (TRIM (units) == 'kg/m2/s'     .OR.                                         &
    TRIM (units) == 'kg m-2 s-1') THEN
  gridbox_emiss = .FALSE.  !  emiss not reported per grid cell

  ! However we allow special units for NOx lightning emissions, because they
  ! are calculated online and given per grid cell. As a consequence the full
  ! conversion cannot be made here. We therefore set gridbox_emiss to TRUE
  ! to make sure that such a correction is done UKCA_NEW_EMISS_CTL.
ELSE IF (TRIM (units) == 'kg/gridbox/s'       .OR.                             &
         TRIM (units) == 'kg gridbox-1 s-1'   .OR.                             &
         TRIM (units) == 'kg/gridcell/s'      .OR.                             &
         TRIM (units) == 'kg gridcell-1 s-1'  .OR.                             &
         TRIM (units) == 'kg/s'               .OR.                             &
         TRIM (units) == 'kg s-1') THEN
  IF (TRIM (tracer_name) == 'NO_lightng') THEN
    gridbox_emiss = .TRUE.
  ELSE
    ierror = 1
  END IF

ELSE
  ierror = 1
END IF

! As indicated above units in all NetCDF files are given as kg m-2 s-1.
! However units of some compounds can be expressed as kg of carbon,
! kg of nitrogen, kg of sulfur, etc. That needs to be specified in
! the standard_name (if it exists), otherwise in the long_name.
! The code below deals with that.

IF (standard_name /= ' ') THEN

  ! If the attribute standard_name is given for an emission field then
  ! check  if the substring 'expressed_as' is included in it.
  ! If that is the case then set l_special_units to TRUE.
  l_special_units = ( INDEX (standard_name, 'expressed_as') > 0 )

  ! If the substring was found then search again in the standard_name
  ! for the appearance of substrings indicating the special units
  ! for emissions (e.g expressed as nitrogen, carbon, sulfur, etc.).
  ! Then call GET_BASE_SCALING to find the appropriate scaling_factor
  ! for the emitted tracer.
  IF (l_special_units) THEN

    ! Check if 'expressed_as_nitrogen'
    IF ( INDEX (standard_name, 'expressed_as_nitrogen') > 0 ) THEN
      special_units = 'expressed_as_nitrogen'
      CALL get_base_scaling (tracer_name, special_units, base_scaling)

      ! Check if 'expressed_as_carbon'
    ELSE IF ( INDEX (standard_name, 'expressed_as_carbon') > 0 ) THEN
      special_units = 'expressed_as_carbon'
      CALL get_base_scaling (tracer_name, special_units, base_scaling)

      ! Check if 'expressed_as_sulfur'
    ELSE IF ( INDEX (standard_name, 'expressed_as_sulfur') > 0 ) THEN
      special_units = 'expressed_as_sulfur'
      CALL get_base_scaling (tracer_name, special_units, base_scaling)

    ELSE
      ! Report that unexpected units appear in the standard_name
      ierror = 2

    END IF

  ELSE

    ! No substrings indicating special units were found in standard_name.
    ! Therefore the emission field is given as kg m-2 s-1 of the emitted
    ! tracer and no conversion needed.
    base_scaling = 1.0

  END IF    ! check if special units used

ELSE IF (long_name /= ' ') THEN

  ! Check if the substring 'expressed as' is included in long_name.
  ! If that is the case then set l_special_units to TRUE.
  l_special_units = ( INDEX (long_name, 'expressed as') > 0 )

  ! If  the substring was found then fill in special_units as above,
  ! but this time using the long_name.
  IF (l_special_units) THEN

    ! Check if 'expressed as nitrogen'
    IF ( INDEX (long_name, 'expressed as nitrogen') > 0 ) THEN
      special_units = 'expressed as nitrogen'
      CALL get_base_scaling (tracer_name, special_units, base_scaling)

      ! Check if 'expressed as carbon'
    ELSE IF ( INDEX (long_name, 'expressed as carbon') > 0 ) THEN
      special_units = 'expressed as carbon'
      CALL get_base_scaling (tracer_name, special_units,  base_scaling)

      ! Check if 'expressed as sulfur'
    ELSE IF ( INDEX (long_name, 'expressed as sulfur') > 0 ) THEN
      special_units = 'expressed as sulfur'
      CALL get_base_scaling (tracer_name, special_units, base_scaling)

    ELSE
      ! Report that unexpected units appear in the long_name
      ierror = 4

    END IF

  ELSE

    ! No special units: emission should be given in kg m-2 s-1 of the
    ! emitted tracer and no conversion needed
    base_scaling = 1.0

  END IF    ! check if special units used

ELSE IF (l_omit_std_name) THEN
  ! Omission of descriptive name is explictly allowed.
  ! Emission should be given in kg m-2 s-1 of the emitted tracer
  ! and no conversion needed
  base_scaling = 1.0

ELSE
  ! Neither standard name or long name are available and omitted name
  ! is not allowed so report error
  ierror = 3

END IF        ! check if standard_name filled

! Stop the model with call to ereport ir errors found
IF (ierror > 0) THEN
  SELECT CASE (ierror)
  CASE (1)
    cmessage = TRIM (field_name)   // ' is expected to have '  //              &
               'other units than ' // TRIM (units)
  CASE (2)
    cmessage = 'Unexpected substring found in standard_name: ' //              &
                standard_name
  CASE (3)
    cmessage = 'long_name needs to be present in the absence ' //              &
               'of standard_name for ' // TRIM (field_name)
  CASE (4)
    cmessage = 'Unexpected substring found in long_name: '     //              &
                long_name
  END SELECT

  CALL ereport ('BASE_EMISS_FACTORS', ierror, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE base_emiss_factors

!---------------------------------------------------------------------
! Description
!   Calculate scaling factors to convert some emission fields.
!
! Method:
!   1) Use the variables tracer_name identifying the atmospheric tracer
!      to which an emission field is mapped (e.g. NO, CH4), and
!      special_units indicating what is currently being emitted
!      (e.g. kg of nitrogen instead of kg of NO, kg of carbon instead
!      of kg of CH4).
!
!   2) Use the molar masses from the module ukca_constants to make the
!      appropriate conversions.
!---------------------------------------------------------------------

SUBROUTINE get_base_scaling (tracer_name, special_units, base_scaling)

IMPLICIT NONE

! Subroutine arguments
CHARACTER (LEN=*),  INTENT(IN)  :: tracer_name
CHARACTER (LEN=*),  INTENT(IN)  :: special_units

REAL,               INTENT(OUT) :: base_scaling  ! scaling factor
                                                 ! for emission field
! Local variables
INTEGER                       :: ierror
CHARACTER(errormessagelength) :: cmessage    ! Error return message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_BASE_SCALING'

REAL, PARAMETER :: ocfact = 1.4    ! Factor to convert from C to POM

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierror = 0

IF (TRIM(special_units) == 'expressed_as_nitrogen'  .OR.                       &
    TRIM(special_units) == 'expressed as nitrogen') THEN
  SELECT CASE (TRIM (tracer_name))
  CASE ('NO', 'NO_aircrft', 'NO_lightng')
    ! Convert from kg(N)  to kg(NO)
    base_scaling  =  m_no / m_n
  CASE ('NH3')
    ! Convert from kg(N) to kg(NH3).
    base_scaling = m_nh3 / m_n
  CASE ('NO2')
    ! Convert from kg(N) to kg(NO2). Note that this is for the future
    ! since NO2 is currently not being emitted in any chemistry scheme.
    base_scaling = m_no2 / m_n
  CASE DEFAULT
    ! Report the need of conversion factor for the unexpected tracer
    ierror   = 1
  END SELECT

ELSE IF (TRIM(special_units) == 'expressed_as_carbon'  .OR.                    &
         TRIM(special_units) == 'expressed as carbon') THEN
  SELECT CASE (TRIM (tracer_name))
  CASE ('C5H8')
    ! Convert from kg(C) to kg(C5H8)
    base_scaling = m_c5h8 / (5.0 * m_c)
  CASE ('CH4')
    ! Convert from kg(C) to kg(CH4)
    base_scaling = m_ch4 / m_c
  CASE ('CH4_wetlnd')
    ! Convert from ug(C) to kg(CH4)
    base_scaling = 1.0e-09 * m_ch4 / m_c
  CASE ('Monoterp')
    ! Convert from kg(C) to kg(C10H16)
    base_scaling = m_monoterp / (10.0 * m_c)
  CASE ('MeOH', 'CH3OH' )
    ! Convert from kg(C) to kg(CH3OH)
    base_scaling = m_ch3oh / m_c
  CASE ('Me2CO')
    ! Convert from kg(C) to kg(CH3COCH3)
    base_scaling = m_me2co / (3.0 * m_c)
    ! In the future this might be revised to use factors such as
    ! me2co_factor=0.134. That might be needed to convert GEIA
    ! NVOC to emissions of the appropriate species.
  CASE ('OM_fossil', 'OM_biofuel', 'OM_biomass')
    ! Convert from mass of carbon to mass of POM
    base_scaling = ocfact
  CASE DEFAULT
    ! Report the need of conversion factor for the unexpected tracer
    ierror   = 1
  END SELECT

ELSE IF (TRIM(special_units) == 'expressed_as_sulfur'  .OR.                    &
         TRIM(special_units) == 'expressed as sulfur') THEN
  SELECT CASE (TRIM (tracer_name))
  CASE ('DMS')
    ! Convert from kg(S) to kg(DMS)
    base_scaling = m_dms / m_s
  CASE ('SO2_low', 'SO2_high', 'SO2_nat', 'SO2')
    ! Convert from kg(S) to kg(SO2)
    base_scaling = m_so2 / m_s
  CASE DEFAULT
    ! Report the need of conversion factor for the unexpected tracer
    ierror   = 1
  END SELECT
ELSE
  ! Report that unexpected special units found
  ierror = 2
END IF


! Stop the model with call to ereport ir errors found
IF (ierror > 0) THEN
  SELECT CASE (ierror)
  CASE (1)
    cmessage = 'Need new conversion factor if emission field '    //           &
                TRIM (tracer_name) // ' is ' // TRIM(special_units)
  CASE (2)
    cmessage = 'Unexpected special units found: '                 //           &
                TRIM (special_units)
  END SELECT

  CALL ereport ('GET_BASE_SCALING', ierror, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE get_base_scaling

!---------------------------------------------------------------------
! Description:
!   Use average reference profiles at given altitudes and create
!   equivalent 3-D profiles to indicate the fraction emitted within
!   each theta grid cell of the model.
!
! Method:
!   Loop over each model grid cell to find out if any layer
!   of the reference profile overlaps with it in the vertical.
!   If so then add the corresponding fraction of that profile
!   to 'vert_scaling_3d', which has been initialised to
!   zero before calling this routine.
!---------------------------------------------------------------------

SUBROUTINE regrid_vertic_profs (                                               &
              nlayers_ref, row_length, rows, model_levels,                     &
              interf_ref, vert_scaling_ref, interf_z,                          &
              field_name, vert_fact, vert_scaling_3d)

IMPLICIT NONE

! Subroutine arguments
!
INTEGER, INTENT(IN) :: nlayers_ref   ! Nr of layers in the ref profile
INTEGER, INTENT(IN) :: row_length    ! Model dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! Reference profiles
REAL, INTENT(IN) :: interf_ref       (0:nlayers_ref) ! altitude of interface
REAL, INTENT(IN) :: vert_scaling_ref (1:nlayers_ref) ! vertical scaling

! Profiles in the model grid: altitude of the interface
REAL, INTENT(IN)    :: interf_z      (1:row_length, 1:rows, 0:model_levels)

CHARACTER (LEN=*), INTENT(IN) :: field_name  ! name of emiss field
CHARACTER (LEN=*), INTENT(IN) :: vert_fact   ! name of the vertical profile

! Profiles in the model grid: vertical scaling applied at mid-theta levels
REAL, INTENT(IN OUT) :: vert_scaling_3d  (1:row_length, 1:rows, 1:model_levels)


! Local variables
INTEGER :: i, j, k, l        ! loop variables

REAL :: cell_lev1, cell_lev2 ! lowest & higest altitude (m) of a model grid cell
REAL :: ref_lev1,  ref_lev2  ! idem for a layer in a reference profile

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='REGRID_VERTIC_PROFS'

LOGICAL, PARAMETER :: test_vertprof = .FALSE.
! True only to test vertical profiles

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Loops over all model grid cells to calculate the emission contribution
! for each of them. Values stored in vert_scaling_3d.
DO k = 0, model_levels-1
  DO j = 1, rows
    DO i = 1, row_length
      ! Lower and upper altitude of a model grid cell
      cell_lev1  = interf_z (i, j, k)
      cell_lev2  = interf_z (i, j, k+1)

      ! Loop through all layers in the reference vertical profile. This
      ! identifies all cases where the model grid cell and the reference layer
      ! overlap. If there is overlap then the percentage of the emission
      ! that should be added to that model level is added to 'vert_scaling_3d'
      DO l = 0, nlayers_ref-1
        ! Lower and upper altitude of a layer in reference profile
        ref_lev1 = interf_ref (l)
        ref_lev2 = interf_ref (l+1)

        ! Check the 4 possible 'overlap options' and add the contributions:
        !  a1: Layer in ref. profile is completely covered by model grid cell
        !  a2: Layer in ref. profile covers upper part of model grid cell
        !  a3: Layer in ref. profile covers lowest part of model grid cell
        !  a4: Layer in ref. profile covers the whole model grid cell
        !
        IF (ref_lev1 >= cell_lev1 .AND. ref_lev1 < cell_lev2) THEN
          ! Bottom of layer in ref. profile is within current model grid cell ...
          IF (ref_lev2 < cell_lev2) THEN
            !and top of layer in ref. profile is also within model grid cell(a1)
            vert_scaling_3d (i,j,k+1) = vert_scaling_3d  (i,j,k+1) +           &
                                        vert_scaling_ref (l+1)
          ELSE
            !and top of layer in ref. profile is higher than model grid cell(a2)
            vert_scaling_3d (i,j,k+1) = vert_scaling_3d (i,j,k+1)  +           &
               vert_scaling_ref (l+1) * (cell_lev2 - ref_lev1) /               &
                                        (ref_lev2  - ref_lev1)
          END IF

        ELSE
          IF (ref_lev1 < cell_lev1 .AND. ref_lev2 > cell_lev1) THEN
            ! Bottom of layer in ref. profile is lower than model grid cell and
            IF (ref_lev2 < cell_lev2) THEN
              ! ... top of layer in ref. profile is within model grid cell (a3)
              vert_scaling_3d (i,j,k+1) = vert_scaling_3d (i,j,k+1) +          &
                 vert_scaling_ref (l+1) * (ref_lev2 - cell_lev1) /             &
                                          (ref_lev2 - ref_lev1)
            ELSE
              ! top of layer in ref. profile is higher than model grid cell (a4)
              vert_scaling_3d (i,j,k+1) = vert_scaling_3d (i,j,k+1) +          &
                 vert_scaling_ref (l+1) * (cell_lev2 - cell_lev1) /            &
                                          (ref_lev2  - ref_lev1)
            END IF
          END IF
        END IF
      END DO     ! nlayers_ref

    END DO      ! i
  END DO       ! j
END DO        ! k

! The switch test_vertprof can be set to TRUE to check that the code works as
! expected for any new profile or model configuration. The code below simply
! checks that the average contribution within a column should be around 1.00
IF (test_vertprof) THEN
  WRITE (umMessage,'(A,A,1X,A,1X,f12.9)') TRIM (field_name),                   &
    ': average vertical contribution for',                                     &
    vert_fact, SUM (vert_scaling_3d) / (rows * row_length)
  CALL umPrint(umMessage,src='regrid_vertic_profs')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE regrid_vertic_profs

END MODULE ukca_emiss_factors
