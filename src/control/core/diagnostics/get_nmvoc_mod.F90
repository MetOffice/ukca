! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module containing the routine GET_NMVOC to calculate the
!   concentration of total non-methane volatile organic compounds
!   (NMVOC) expressed as micro-g m-3 of carbon.
!
! Method:
!   Retrieve the mass mixing ratios (MMR) for all contributing species
!   in a given chemistry scheme and return the total, converting to mass
!   concentration of carbon.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE get_nmvoc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GET_NMVOC_MOD'

CONTAINS

SUBROUTINE get_nmvoc (row_length, rows, model_levels, ntracers, tracer,        &
                      pres, temp, nmvoc_3d)

USE ukca_um_legacy_mod,     ONLY: rgas => r
  ! Specific gas constant for dry air (J kg-1 K-1)

USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_constants,         ONLY: m_c,                                         &
                                  m_hcho, m_c2h6, m_mecho, m_c3h8,             &
                                  m_me2co, m_c5h8, m_mgly, m_mvk,              &
                                  m_ch3oh, m_c2h4, m_c3h6, m_c4h10,            &
                                  m_mek, m_toluene, m_gly, m_oxylene
USE ukca_cspecies,          ONLY: n_hcho, n_c2h6, n_mecho, n_c3h8,             &
                                  n_me2co, n_c5h8, n_mgly, n_mvk,              &
                                  n_ch3oh, n_c2h4, n_c3h6, n_c4h10,            &
                                  n_mek, n_toluene, n_gly, n_oxylene
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,               ONLY: jpim,  jprb     ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN)  :: row_length        ! Model dimensions
INTEGER, INTENT(IN)  :: rows
INTEGER, INTENT(IN)  :: model_levels
INTEGER, INTENT(IN)  :: ntracers          ! No. of tracers

REAL,    INTENT(IN)  :: tracer (row_length, rows, model_levels, ntracers)
! Tracers: mixing ratios, generally MMRs (kg kg-1)

REAL,    INTENT(IN)  :: pres (row_length, rows, model_levels)
! Pressure (Pa)

REAL,    INTENT(IN)  :: temp (row_length, rows, model_levels)
! Temperature (K)

REAL,    INTENT(OUT) :: nmvoc_3d (row_length, rows, model_levels)
! Total non-methane VOC, NMVOC (ug C m-3)

! Local variables

INTEGER                            :: ierr      ! Arguments for
CHARACTER (LEN=errormessagelength) :: cmessage  ! ereport

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_NMVOC'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierr = 0
nmvoc_3d (:,:,:) = 0.0

IF (ukca_config%l_ukca_raq .OR. ukca_config%l_ukca_raqaero) THEN

  ! The sum of the quotient carbon mass mixing ratio (MMR) / atomic mass of C
  ! is calculated for all relevant species present in the chemistry scheme
  ! (see array chch_defs). This quotient is the species MMR * no. of C atoms
  ! per molecule / molecular mass.
  ! The sum is then converted to NMVOC concentration by multiplying by the
  ! atomic mass of carbon and the density of air and scaled to the correct
  ! units.

  nmvoc_3d(:,:,:) = (                                                          &
                     tracer(:,:,:,n_hcho) / m_hcho +                           &
                     tracer(:,:,:,n_c2h6) * 2.0 / m_c2h6 +                     &
                     tracer(:,:,:,n_mecho) * 2.0 / m_mecho +                   &
                     tracer(:,:,:,n_c3h8) * 3.0 / m_c3h8 +                     &
                     tracer(:,:,:,n_me2co) * 3.0 / m_me2co +                   &
                     tracer(:,:,:,n_c5h8) * 5.0 / m_c5h8 +                     &
                     tracer(:,:,:,n_mgly) * 3.0 / m_mgly +                     &
                     tracer(:,:,:,n_mvk) * 4.0 / m_mvk +                       &
                     tracer(:,:,:,n_ch3oh) / m_ch3oh +                         &
                     tracer(:,:,:,n_c2h4) * 2.0 / m_c2h4 +                     &
                     tracer(:,:,:,n_c3h6) * 3.0 / m_c3h6 +                     &
                     tracer(:,:,:,n_c4h10) * 4.0 / m_c4h10 +                   &
                     tracer(:,:,:,n_mek) * 4.0 / m_mek +                       &
                     tracer(:,:,:,n_toluene) * 7.0 / m_toluene +               &
                     tracer(:,:,:,n_gly) * 2.0 / m_gly +                       &
                     tracer(:,:,:,n_oxylene) * 8.0 / m_oxylene                 &
                    ) *                                                        &
                    m_c * 1.0e9 * pres(:,:,:) / (rgas * temp(:,:,:))

ELSE
  ierr     = 1
  cmessage ='NMVOC diagnostic not supported for this chemistry scheme yet'
  CALL ereport ('GET_NMVOC', ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN

END SUBROUTINE get_nmvoc

END MODULE get_nmvoc_mod
