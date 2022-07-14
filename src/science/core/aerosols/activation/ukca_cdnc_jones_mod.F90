! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Calculate CDNC using Jones method doi:10.1038/370450a0
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE ukca_cdnc_jones_mod

!USE

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CDNC_JONES_MOD'

CONTAINS

SUBROUTINE ukca_cdnc_jones ( nbox,                                             &
                             act,                                              &
                             drydp,                                            &
                             nd,                                               &
                             ccn_1,                                            &
                             cdnc )

USE parkind1,         ONLY:                                                    &
    jpim,                                                                      &
    jprb

USE ukca_mode_setup,  ONLY:                                                    &
    mode,                                                                      &
    modesol,                                                                   &
    nmodes,                                                                    &
    sigmag

USE umErf_mod,        ONLY:                                                    &
    umErf

USE ukca_um_legacy_mod,      ONLY:                                             &
    exp_v,                                                                     &
    log_v

USE yomhook,          ONLY:                                                    &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN)  :: nbox
REAL,    INTENT(IN)  :: act                ! radius for activation (m)
REAL,    INTENT(IN)  :: drydp(nbox,nmodes) ! Dry diameter
REAL,    INTENT(IN)  :: nd(nbox,nmodes)    ! Number density
REAL,    INTENT(OUT) :: ccn_1(nbox)        ! CCN concentration
REAL,    INTENT(OUT) :: cdnc(nbox)         ! Cloud Droplet Number Concentration

! Local variables

INTEGER :: imode
INTEGER :: i

REAL    :: dp0                ! Diam (nm)
REAL    :: root2              ! 2^(1/2)
REAL    :: log_sigmag(nmodes) ! Log of fixed geometric standard dev of each mode
REAL    :: erf_arg            ! Error Fn argument
REAL    :: erfterm            !
REAL    :: ccn_1_exp(nbox)    ! exp of ccn_1

REAL, PARAMETER :: cdnmin = 5.0 ! Min CDN (no/cm^-3)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName='UKCA_CDNC_JONES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

! This is the Jones method for calculating CDNC, see doi:10.1038/370450a0

CALL log_v(nmodes,sigmag,log_sigmag)

root2 = SQRT(2.0)

! Initialise to zero
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i)                    &
!$OMP SHARED(nbox, ccn_1)
DO i = 1, nbox
  ccn_1(i) = 0.0
END DO
!$OMP END PARALLEL DO

DO imode=1,nmodes
  IF (mode(imode)) THEN
    IF (modesol(imode) == 1) THEN
      ! for CCN_1 take CCN for particles > ACT dry radius
      dp0=2.0*act
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,erf_arg, erfterm)   &
!$OMP SHARED(nbox, ccn_1, drydp,imode, log_sigmag, nd, dp0, root2)
      DO i = 1, nbox
        erf_arg = LOG(dp0/drydp(i,imode)) / (root2*log_sigmag(imode))
        erfterm = 0.5 * nd(i,imode) * (1.0-umErf(erf_arg))
        ccn_1(i)   = ccn_1(i) - (0.0025*erfterm)
      END DO
!$OMP END PARALLEL DO
    END IF
  END IF
END DO

CALL exp_v( nbox , ccn_1, ccn_1_exp )

! Jones method, see doi:10.1038/370450a0
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i)                    &
!$OMP SHARED(nbox, ccn_1_exp, cdnc)
DO i = 1, nbox
  cdnc(i)=375.0*(1.0-ccn_1_exp(i))
  IF (cdnc(i) < cdnmin ) cdnc(i) = cdnmin
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_cdnc_jones

END MODULE ukca_cdnc_jones_mod
