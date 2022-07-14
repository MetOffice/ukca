! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description:
!  The routine UKCA_SET_DIURNAL_OX distributes the mean oxidant concentrations
!  read in from external files to give a diurnally varying concentration.
!  Part of the offline oxidant chemistry scheme. A routine to calculate the time
!  integral of cos(zenith angle) is also included (UKCA_INT_COSZ).
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds, University of Oxford,
!  and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! Contained subroutines:
!   ukca_set_diurnal_ox
!   ukca_int_cosz
!
!-----------------------------------------------------------------------
!
MODULE ukca_diurnal_oxidant

USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim

IMPLICIT NONE
PRIVATE

! cos (zenith angle) and integral of cos(zenith angle) for diurnal variation
REAL, PUBLIC, ALLOCATABLE, SAVE :: cosza(:)          ! cos(za)
REAL, PUBLIC, ALLOCATABLE, SAVE :: intcosza(:)       ! integral cos(za)
REAL, PUBLIC, ALLOCATABLE, SAVE :: daylength1d(:)    ! radians

PUBLIC ukca_set_diurnal_ox                ! called from asad_fyinit
PUBLIC ukca_set_diurnal_ox_col            ! called from asad_fyinit
PUBLIC ukca_int_cosz                      ! called from ukca_main1
PUBLIC dealloc_diurnal_oxidant            ! called from ukca_main1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_DIURNAL_OXIDANT'

CONTAINS

! #############################################################################

SUBROUTINE ukca_set_diurnal_ox(species, field, n_points, nlev)

! This routine provides code to add a diurnal profile to oxidants provided
! as mean fields. This varies with each oxidant. In the case of OH and HO2
! the diurnal profile is proportional to cos(zenith angle).

USE ukca_chem_offline,        ONLY: o3_offline_diag, oh_offline_diag,          &
                                    no3_offline_diag, ho2_offline_diag
USE ukca_um_legacy_mod,       ONLY: rsec_per_day, pi
USE ereport_mod,              ONLY: ereport

USE errormessagelength_mod,   ONLY: errormessagelength
IMPLICIT NONE

INTEGER, INTENT(IN)  :: n_points                ! field dimension
INTEGER, INTENT(IN)  :: nlev                    ! model level
CHARACTER (LEN=10), INTENT(IN)  :: species      ! Species char strng

REAL,    INTENT(IN OUT) :: field(n_points)       ! field to be scaled

! Local variables

REAL, PARAMETER :: daylen_cutoff = 3600.0
REAL            :: daylen_secs(n_points)        ! daylength (s)

CHARACTER (LEN=errormessagelength) :: cmessage                  ! Error message
INTEGER            :: errcode                   ! Variable passed to ereport

REAL(KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_DIURNAL_OX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE(species(1:4))
CASE ('OH  ','HO2 ')

  ! Scale the field according to cos(zenith angle)/integral(cos(zenith_angle))
  ! Night time values are zero
  WHERE ((cosza > 0.0) .AND. (intcosza > 1.0e-1))
    field(:) = field(:)*rsec_per_day*(cosza(:)/intcosza(:))
  ELSE WHERE
    field(:) = 0.0
  END WHERE

CASE ('NO3 ')

  ! convert from radians to seconds
  daylen_secs = daylength1d*rsec_per_day/(2.0*pi)

  ! With short day length ensures no NO3 diurnal variation done
  WHERE (daylen_secs < daylen_cutoff)
    daylen_secs = rsec_per_day
  END WHERE

  ! Restrict NO3 to nightime
  WHERE ((cosza > 0.0) .AND. (intcosza > 1.0e-1) .AND.                         &
    (daylen_secs < rsec_per_day))
    field(:) = 0.0
  END WHERE

CASE ('O3  ')
  ! No diurnal profile for O3

CASE DEFAULT

  cmessage = 'ukca_set_diurnal_ox: Species: '//species//' not found'
  errcode = 1
  CALL ereport('ukca_set_diurnal_ox',errcode,cmessage)

END SELECT

! Get rid of any negative values
WHERE (field < 0.0)
  field = 0.0
END WHERE

! Set diagnostic fields
SELECT CASE(species(1:4))
CASE ('O3  ')
  o3_offline_diag(:,nlev) = field(:)
CASE ('OH  ')
  oh_offline_diag(:,nlev) = field(:)
CASE ('NO3 ')
  no3_offline_diag(:,nlev) = field(:)
CASE ('HO2 ')
  ho2_offline_diag(:,nlev) = field(:)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_diurnal_ox

! ######################################################################

SUBROUTINE ukca_set_diurnal_ox_col(species, field, n_points, ix, jy)

! This routine provides code to add a diurnal profile to oxidants provided
! as mean fields. This varies with each oxidant. In the case of OH and HO2
! the diurnal profile is proportional to cos(zenith angle).

USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_chem_offline,      ONLY: o3_offline_diag, oh_offline_diag,            &
                                  no3_offline_diag, ho2_offline_diag
USE ukca_um_legacy_mod,     ONLY: pi, rsec_per_day
USE ereport_mod,            ONLY: ereport

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN)  :: n_points                ! field dimension
INTEGER, INTENT(IN) :: ix                       ! i counter
INTEGER, INTENT(IN) :: jy                       ! j counter
CHARACTER (LEN=10), INTENT(IN)  :: species      ! Species char strng

REAL,    INTENT(IN OUT) :: field(n_points)       ! field to be scaled

! Local variables

REAL, PARAMETER :: daylen_cutoff = 3600.0
REAL            :: daylen_secs        ! daylength (s)

CHARACTER (LEN=errormessagelength) :: cmessage                  ! Error message
INTEGER            :: errcode                   ! Variable passed to ereport

INTEGER :: location ! mapping to theta_field

REAL(KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_DIURNAL_OX_COL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! get mapping from theta_field to proper x*y. Field is 1 column here!
location = ix + ((jy - 1)*ukca_config%row_length)

SELECT CASE(species(1:4))
CASE ('OH  ','HO2 ')

  ! Scale the field according to cos(zenith angle)/integral(cos(zenith_angle))
  ! Night time values are zero
  IF ((cosza(location) > 0.0) .AND. (intcosza(location) > 1.0e-1)) THEN
    field(:) = field(:)*rsec_per_day*(cosza(location) / intcosza(location))
  ELSE
    field(:) = 0.0
  END IF

CASE ('NO3 ')

  ! convert from radians to seconds
  daylen_secs = daylength1d(location)*rsec_per_day / (2.0*pi)

  ! With short day length ensures no NO3 diurnal variation done
  IF (daylen_secs < daylen_cutoff) daylen_secs = rsec_per_day

  ! Restrict NO3 to nightime
  IF (((cosza(location) > 0.0) .AND. (intcosza(location) > 1.0e-1)) .AND.      &
       (daylen_secs < rsec_per_day)) THEN
    field(:) = 0.0
  END IF

CASE ('O3  ')
  ! No diurnal profile for O3

CASE DEFAULT

  cmessage = 'ukca_set_diurnal_ox: Species: '//species//' not found'
  errcode = 1
  CALL ereport('ukca_set_diurnal_ox',errcode,cmessage)

END SELECT

! Get rid of any negative values
WHERE (field < 0.0)
  field = 0.0
END WHERE

! Set diagnostic fields
SELECT CASE(species(1:4))
CASE ('O3  ')
  o3_offline_diag(location,:) = field(:)
CASE ('OH  ')
  oh_offline_diag(location,:) = field(:)
CASE ('NO3 ')
  no3_offline_diag(location,:) = field(:)
CASE ('HO2 ')
  ho2_offline_diag(location,:) = field(:)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_diurnal_ox_col

! ######################################################################

SUBROUTINE ukca_int_cosz(row_length, rows,                                     &
                         sin_latitude, cos_latitude, tan_latitude, sindec,     &
                         cos_zenith_angle, int_zenith_angle)

!  Description:
!   Calculate integral of cos(solar zenith angle) dt over a period of
!   one day multiplied by seconds per day.
!   See Gates, D.M. (1980) Biophysical Ecology, Springer-Verlag, p103-4
!   Only depends on latitude and solar declination:
!   int_za = (daysecs/pi)*(hs*sin(phi)*sin(decls)+cos(phi)*cos(decls)*sin(hs)),
!   where hs is day length / 2.0.
!
! ######################################################################

USE ukca_um_legacy_mod,  ONLY: pi, rsec_per_day

IMPLICIT NONE

INTEGER, INTENT(IN) ::  row_length       ! field dimension
INTEGER, INTENT(IN) ::  rows             ! field dimension

REAL, INTENT(IN)    :: sin_latitude(row_length, rows)    ! SIN(latitude)
REAL, INTENT(IN)    :: cos_latitude(row_length, rows)    ! COS(latitude)
REAL, INTENT(IN)    :: tan_latitude(row_length, rows)    ! TAN(latitude)
REAL, INTENT(IN)    :: sindec            ! Sin(solar declination)
REAL, INTENT(IN)    :: cos_zenith_angle(row_length, rows)
! Cosine of zenith angle

REAL, INTENT(OUT)   :: int_zenith_angle(row_length, rows)
! Integral of zenith angle

! Local variables
REAL            :: daylength(row_length,rows)  ! Daylength
REAL            :: decls                       ! solar declination
REAL            :: tan_decls                   ! TAN(solar declination)
REAL            :: hs(row_length,rows)         ! half day length (radians)
REAL            :: cos_hs(row_length,rows)     ! cos of hs

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INT_COSZ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

decls = ASIN(sindec)

tan_decls = TAN(decls)
IF (ABS(tan_decls) < EPSILON(0.0)) THEN
  cos_hs = 0.0
ELSE
  cos_hs = -tan_latitude*tan_decls
END IF

! Polar summer, set cos_hs ~ -1
WHERE (ABS(cos_hs) > 1.0 .AND. sin_latitude*decls > 0.0)
  cos_hs = -0.99999999
END WHERE

WHERE (ABS(cos_hs) < 1.0)
  hs = ACOS(cos_hs)
  int_zenith_angle = (rsec_per_day/pi) *                                       &
    (hs * sin_latitude * sindec + cos_latitude * COS(decls) * SIN(hs))
  daylength = 2.0*hs
ELSE WHERE
  ! Polar night
  int_zenith_angle = 0.0
  daylength = 0.0
END WHERE

! Make 1-D arays for set_ukca_diurnal_ox routine
IF (.NOT. ALLOCATED(intcosza))    ALLOCATE(intcosza(row_length*rows))
IF (.NOT. ALLOCATED(cosza))       ALLOCATE(cosza(row_length*rows))
IF (.NOT. ALLOCATED(daylength1d)) ALLOCATE(daylength1d(row_length*rows))

intcosza = RESHAPE(int_zenith_angle,[row_length*rows])
cosza = RESHAPE(cos_zenith_angle,[row_length*rows])
daylength1d = RESHAPE(daylength,[row_length*rows])

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_int_cosz

! ######################################################################

SUBROUTINE dealloc_diurnal_oxidant()
! ----------------------------------------------------------------------
! Description:
! Deallocate the saved cosza, intcosza and daylength1d arrays
! ----------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'DEALLOC_DIURNAL_OXIDANT'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(daylength1d)) DEALLOCATE(daylength1d)
IF (ALLOCATED(cosza))       DEALLOCATE(cosza)
IF (ALLOCATED(intcosza))    DEALLOCATE(intcosza)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE dealloc_diurnal_oxidant
END MODULE ukca_diurnal_oxidant
