! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!   Module containing subroutine for calculating the sea-to-air flux of
!   DMS (as kg m-2 s-1 of sulphur) as a function of its concentration in
!   in seawater and windspeed.
!
! Method:
!
!   The sea-air exchange can be determined according to one of three
!   commonly-used parametrization schemes, those of Liss & Merlivat (1986),
!   Wanninkhof (1992) and Nightingale et al. (2000).
!
!   The Schmidt number for DMS is calculated as in Saltzman et al. (1993)
!   and used with the windspeed to determine the mass transfer (or
!   "piston") velocity according to the desired parametrization.
!   This is then used to determine the sea-air mass flux of DMS
!   as a function of sea-water DMS concentration. High surface
!   temperatures (caused by the land portion of a gridbox when
!   coastal tiling is not active) cause negative Sc values which
!   would give a floating-point error in the k_DMS calculation,
!   so the Tstar values are capped. This shouldn't be a problem
!   when coastal tiling is on as then the Tstar values passed in
!   are those for sea only.
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
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_dms_flux_mod

USE ukca_config_specification_mod, ONLY: i_liss_merlivat, i_wanninkhof,        &
                                         i_nightingale, zerodegc

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_DMS_FLUX_MOD'

CONTAINS

SUBROUTINE ukca_dms_flux(row_length, rows, wind_10m, tstar, land_fract,        &
                         dms_conc, i_dms_flux, f_dms)

USE water_constants_mod, ONLY: tfs

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umPrint, umMessage
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows

REAL, INTENT(IN) :: wind_10m(row_length, rows)   ! 10m windspeed (ms-1)
REAL, INTENT(IN) :: tstar(row_length, rows)      ! Surface temperature (K)
REAL, INTENT(IN) :: land_fract(row_length, rows) ! Fraction of land in gridbox
REAL, INTENT(IN) :: dms_conc(row_length, rows)   ! Concentraction of DMS in
                                                 ! seawater (nmol l-1)
INTEGER, INTENT(IN) :: i_dms_flux                ! Option code to select
                                                 ! parameterization scheme

REAL, INTENT(OUT) :: f_dms(row_length, rows)     ! Sea-air flux of DMS
                                                 ! (kg[S] m-2 s-1)

! Local variables:
INTEGER  :: i,j                     ! Loop counters
REAL     :: sc(row_length, rows)    ! Schmidt number
REAL     :: k_dms(row_length, rows) ! Piston velocity of DMS (cm h-1)
REAL     :: t_c                     ! Surface temperature (degrees Celsius)
REAL     :: k_600                   ! Piston vel. for gas with Schmidt no.= 600
REAL     :: k_660                   ! Piston vel. for gas with Schmidt no.= 660
REAL     :: n                       ! Schmidt number exponent
REAL, PARAMETER :: t_max=47.0       ! Max T to avoid breaking Sc fit (C)

INTEGER  :: errcode                              ! Error code
CHARACTER(LEN=errormessagelength) :: cmessage    ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DMS_FLUX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the Schmidt number (Sc):

DO j = 1, rows
  DO i = 1, row_length
    t_c = MIN((MAX(tstar(i, j), tfs) - zerodegc), t_max)
    sc(i, j) = 2674.0 - (147.12*t_c) + (3.726*t_c**2)                          &
      - (0.038*t_c**3)
  END DO
END DO

! Determine the mass transfer (or "piston") velocity (k_DMS) over sea
! according to the specified parametrization scheme:

SELECT CASE (i_dms_flux)

CASE (i_liss_merlivat)
  DO j = 1, rows
    DO i = 1, row_length
      IF (wind_10m(i, j)  <=  3.6) THEN
        k_600 = 0.17 * wind_10m(i, j)
        n = -2.0/3.0
      END IF
      IF (wind_10m(i, j)  >   3.6 .AND.                                        &
        wind_10m(i, j)  <=  13.0) THEN
        k_600 = (2.85*wind_10m(i, j)) - 9.65
        n = -0.5
      END IF
      IF (wind_10m(i, j)  >   13.0) THEN
        k_600 = (5.9*wind_10m(i, j)) - 49.3
        n = -0.5
      END IF
      IF (land_fract(i, j)  <   1.0) THEN
        k_dms(i, j) = k_600 * (sc(i, j)/600.0)**n
      ELSE
        k_dms(i, j) = 0.0
      END IF
    END DO
  END DO

CASE (i_wanninkhof)
  DO j = 1, rows
    DO i = 1, row_length
      k_660 = 0.31 * wind_10m(i, j)**2
      n = -0.5
      IF (land_fract(i, j)  <   1.0) THEN
        k_dms(i, j) = k_660 * (sc(i, j)/660.0)**n
      ELSE
        k_dms(i, j) = 0.0
      END IF
    END DO
  END DO

CASE (i_nightingale)
  DO j = 1, rows
    DO i = 1, row_length
      k_600 = (0.222*wind_10m(i, j)**2) + (0.333*wind_10m(i, j))
      n = -0.5
      IF (land_fract(i, j)  <   1.0) THEN
        k_dms(i, j) = k_600 * (sc(i, j)/600.0)**n
      ELSE
        k_dms(i, j) = 0.0
      END IF
    END DO
  END DO

CASE DEFAULT
  WRITE(cmessage,'(A,I4,A)') 'i_dms_flux = ', i_dms_flux,                      &
                              ' is not a valid choice.'
  errcode = 1
  CALL ereport(RoutineName, errcode, cmessage)

END SELECT

! Finally, calculate the sea-air flux of DMS as a function of k_DMS
! and dissolved DMS concentration. The former requires a conversion
! from cm hour-1 to ms-1, and the latter from nanomoles per litre to
! kg[S] m-3, to return the flux in kg[S] m-2 sec-1.

DO j = 1, rows
  DO i = 1, row_length
    f_dms(i, j) = (k_dms(i, j) / 3.6e5)                                        &
      * (dms_conc(i, j) * 32.0e-9)
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_dms_flux

END MODULE ukca_dms_flux_mod
