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
! Purpose: Species labelled 'CF': ASAD will treat the species as a constant
!     but will call this routine so that the user may set the
!     values differently at each gridpoint for example. Currently
!     used for setting the water vapour and CO2 concentrations, as well
!     as the offline oxidants.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_FYINIT
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     Interface
!     On entry, the following will be set:
!              species - character name of species to set.
!                        Will be the same as listed in chch.d file
!              klen    - length of array, y_out.
!
!     On exit, the following must be set:
!              y_out   - Array of points to set for the species.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_inicnt_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_INICNT_MOD'

CONTAINS

SUBROUTINE asad_inicnt( species, y_out, klen, nlev )

USE asad_mod,              ONLY: wp, co2, tnd, nlfro2, f, jpro2
USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_constants,        ONLY: c_oh, c_o3, c_no3, c_ho2
USE ukca_environment_fields_mod, ONLY: o3_offline, oh_offline,                 &
                           no3_offline, ho2_offline
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: klen      ! No of spatial points
INTEGER, INTENT(IN) :: nlev      ! Model level

CHARACTER (LEN=10), INTENT(IN)  :: species  ! Species char strng

REAL, INTENT(OUT)   :: y_out(klen)! Species concentration this may
                                  ! be in volumetric mixing ratio
                                  ! units (H2O) or as mass mixing
                                  ! ratio (offline oxidants)

!       Local variables

INTEGER :: errcode                ! Variable passed to ereport
INTEGER :: row_length             ! row_length for theta field
INTEGER :: rows                   ! rows for theta field
INTEGER :: iro2                   ! Counter for RO2 species
INTEGER :: j                      ! Loop variable

CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

REAL                          :: fro2(klen) ! Total RO2 concentration
                                            ! (molecules/cm3)

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_INICNT'


!     1.  Copy water, CO2, and offline oxidants (if required) into ASAD array.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

row_length = ukca_config%row_length
rows = ukca_config%rows

IF (species(1:4) == 'CO2 ' .AND. ukca_config%l_chem_environ_co2_fld) THEN
  !  The CO2 field is set to the UM prognostic if it is available and
  !  is used in the chemical scheme.
  y_out(:) = co2(:)
ELSE IF ( species(1:4) == 'H2O ' ) THEN
  ! Note that wp is in units of volumetric mixing ratio
  y_out(:) = wp(:)

  ! First call to calculate total RO2, using initial concentration
  ! of all RO2 species in mechanism
ELSE IF (species(1:4) == 'RO2 ') THEN
  IF (ukca_config%l_ukca_ro2_perm) THEN
    fro2(:) = 0.0
    ! Loop through all RO2 species
    DO j = 1, jpro2
      ! Get index location of each RO2 species and sum
      iro2    = nlfro2(j)
      fro2(:) = fro2(:) + f(:, iro2)
    END DO   ! End iteration over RO2 species
    ! Convert to VMR - fro2 will be in molecules/cm3
    y_out(:) = fro2(:)/tnd(:)
  ELSE
    errcode = 126
    cmessage = 'RO2 should only be a species if l_ukca_ro2_perm == T'
    CALL ereport('ASAD_INICNT',errcode,cmessage)
  END IF

ELSE IF (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) THEN
  ! These species are converted from mass mixing ratio to vmr
  IF (species(1:4) == 'OH  ') THEN
    y_out(:) = RESHAPE(oh_offline(1:row_length,1:rows,nlev),[klen])
    y_out(:) = y_out(:)/c_oh
  ELSE IF (species(1:4) == 'O3  ') THEN
    y_out(:) = RESHAPE(o3_offline(1:row_length,1:rows,nlev),[klen])
    y_out(:) = y_out(:)/c_o3
  ELSE IF (species(1:4) == 'NO3 ') THEN
    y_out(:) = RESHAPE(no3_offline(1:row_length,1:rows,nlev),[klen])
    y_out(:) = y_out(:)/c_no3
  ELSE IF (species(1:4) == 'HO2 ') THEN
    y_out(:) = RESHAPE(ho2_offline(1:row_length,1:rows,nlev),[klen])
    y_out(:) = y_out(:)/c_ho2
  ELSE
    errcode = 125
    cmessage = ' Species '//species//' is not treated by this routine'
    CALL ereport('ASAD_INICNT',errcode,cmessage)
  END IF
ELSE
  errcode=124
  cmessage= ' Species '//species//' not treated by this routine'
  CALL ereport('ASAD_INICNT',errcode,cmessage)
END IF

! Convert to molecules/cm^3 from vmr
y_out(:) = y_out(:)*tnd(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE asad_inicnt
END MODULE asad_inicnt_mod
