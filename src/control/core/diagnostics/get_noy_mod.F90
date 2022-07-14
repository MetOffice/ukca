! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module containing the function GET_NOY to calculate the volume mixing
!   ratio (VMR) of total reactive nitrogen (NOy).
!
! Method:
!   Look for all species contributing to NOy in a given chemistry
!   scheme, get the corresponding mass mixing ratios (MMR), convert
!   them to volume mixing ratios (VMR) and add them up to get
!   NOy (expressed as VMR).
!   The current code can be easily extended for more chemistry schemes
!   (one only needs to look for the species contributing to NOy
!   and make sure the right number of nitrogen atoms are included in
!   the calculation).
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
!
MODULE get_noy_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GET_NOY_MOD'

CONTAINS

SUBROUTINE get_noy (row_length, rows, model_levels, ntracers, tracer, noy_3d)

USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_cspecies,          ONLY: c_species, n_no, n_no3, n_no2, n_n2o5,       &
                                  n_ho2no2, n_hono2, n_pan, n_ison, n_orgnit,  &
                                  n_rnc2h4, n_rnc3h6
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,               ONLY: jpim,  jprb     ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN)  :: row_length        ! Model dimensions
INTEGER, INTENT(IN)  :: rows
INTEGER, INTENT(IN)  :: model_levels
INTEGER, INTENT(IN)  :: ntracers          ! no. of tracers

REAL,    INTENT(IN)  :: tracer (row_length, rows, model_levels, ntracers)
! tracers (MMR, i.e. g g-1 or kg kg-1)

REAL,    INTENT(OUT) :: noy_3d (row_length, rows, model_levels)
! Total reactive nitrogen: NOy (VMR, i.e. mol mol-1)


! Local variables

INTEGER                            :: ierr      ! Arguments for
CHARACTER (LEN=errormessagelength) :: cmessage  ! ereport

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_NOY'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierr           = 0
noy_3d (:,:,:) = 0

IF (ukca_config%l_ukca_raq) THEN
  ! NOy is calculated as the sum of the volume mixing ratios (VMR) of all
  ! nitrogen species present in the chemistry scheme (see array chch_defs),
  ! multiplied by the number of nitrogen atoms they contain.
  ! Since the array tracer is given in MMR it is converted here to VMR by
  ! dividing by c_species (= m_species / m_air).
  noy_3d (:,:,:) = tracer (:,:,:, n_no)         / c_species (n_no)     +       &
                   tracer (:,:,:, n_no3)        / c_species (n_no3)    +       &
                   tracer (:,:,:, n_no2)        / c_species (n_no2)    +       &
                   tracer (:,:,:, n_n2o5) * 2.0 / c_species (n_n2o5)   +       &
                   tracer (:,:,:, n_ho2no2)     / c_species (n_ho2no2) +       &
                   tracer (:,:,:, n_hono2 )     / c_species (n_hono2)  +       &
                   tracer (:,:,:, n_pan)        / c_species (n_pan)    +       &
                   tracer (:,:,:, n_ison)       / c_species (n_ison)   +       &
                   tracer (:,:,:, n_orgnit)     / c_species (n_orgnit) +       &
                   tracer (:,:,:, n_rnc2h4 )    / c_species (n_rnc2h4) +       &
                   tracer (:,:,:, n_rnc3h6 )    / c_species (n_rnc3h6)
ELSE
  ierr     = 1
  cmessage ='NOy diagnostic not supported for this chemistry scheme yet'
  CALL ereport ('GET_NOY', ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN

END SUBROUTINE get_noy

END MODULE get_noy_mod
