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
! Purpose: Subroutine to overwrite values at top of model
!          using interpolated 5-day fields from the 2-d model.
!          Based on STRATF.F from Cambridge TOMCAT model and
!          modified by Olaf Morgenstern to allow for flexible
!          positioning of the NOy species.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_CHEMISTRY_CTL.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_stratf_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_STRATF_MOD'

CONTAINS

SUBROUTINE ukca_stratf(row_length, rows,                                       &
                       model_levels,                                           &
                       ntracer,                                                &
                       env_ozone3d, tracer)

USE ukca_tropopause,   ONLY: tropopause_level
USE ukca_cspecies,     ONLY: n_o3, n_o3s, n_hono2
USE ukca_constants,    ONLY: c_n, c_hno3
USE parkind1,          ONLY: jprb, jpim
USE yomhook,           ONLY: lhook, dr_hook

USE ereport_mod,       ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Diag

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length         ! No of longitudes
INTEGER, INTENT(IN) :: rows               ! No of latitudes
INTEGER, INTENT(IN) :: model_levels       ! No of levels
INTEGER, INTENT(IN) :: ntracer            ! No of chemical tracers

REAL, INTENT(IN) :: env_ozone3d(row_length,rows,model_levels) ! O3

! Tracer concentrations in mass mixing ratio
REAL, INTENT(IN OUT) :: tracer(row_length,rows,model_levels,ntracer)

!     Local variables

INTEGER :: l            ! Loop variables

LOGICAL :: mask(row_length,rows,model_levels) ! mask to identify stratosphere

REAL, PARAMETER :: o3_hno3_ratio = 1.0/1000.0 ! kg[N]/kg[O3] from
                                              ! Murphy and Fahey 1994

REAL :: o33d(row_length,rows,model_levels)   ! 2D field interpolated onto 3D
REAL :: hno33d(row_length,rows,model_levels) ! 3D field from fixed o3:hno3 ratio

! Parameter to overwrite stratosphere (fixed no of levels above tropopause)
INTEGER, PARAMETER :: no_above_trop1 = 3        ! Suitable for L38/L60
INTEGER, PARAMETER :: no_above_trop2 = 10       ! Suitable for L63/L85
INTEGER :: no_above_trop

INTEGER           :: errcode                    ! Error code: ereport
CHARACTER(LEN=errormessagelength) :: cmessage                   ! Error message
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_STRATF'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set number of levels above tropopause depending on vertical resolution
IF (model_levels == 38 .OR. model_levels == 60) THEN
  no_above_trop = no_above_trop1
ELSE IF (model_levels == 63 .OR. model_levels == 70 .OR.                       &
         model_levels == 85) THEN
  no_above_trop = no_above_trop2
ELSE
  errcode = 1
  cmessage = 'Levels above tropopause not set at this resolution'
  CALL ereport('UKCA_STRATF',errcode,cmessage)
END IF

IF (printstatus == Prstatus_Diag) THEN
  WRITE(umMessage,'(A)') 'UKCA_STRATF: Logicals in use:'
  CALL umPrint(umMessage,src='ukca_stratf')
  WRITE(umMessage,'(A,I0)') 'no_above_trop: ',no_above_trop
  CALL umPrint(umMessage,src='ukca_stratf')
  WRITE(umMessage,'(A)') ' '
  CALL umPrint(umMessage,src='ukca_stratf')
END IF

o33d(:,:,:) = env_ozone3d(:,:,:)
hno33d(:,:,:) = o33d(:,:,:)*o3_hno3_ratio*c_hno3/c_n

! Overwrite o3, and hno3 at all gridboxes a fixed
!  number of model levels above the tropopause
!  O3   - UM ancillary
!  HNO3 - using fixed o3:hno3 ratio

mask(:,:,:) = .FALSE.
DO l = model_levels,1,-1
  mask(:,:,l) = tropopause_level(:,:)+no_above_trop <= l
END DO

WHERE (mask(:,:,:))
  tracer(:,:,:,n_o3)    = o33d(:,:,:)
  tracer(:,:,:,n_hono2) = hno33d(:,:,:)
END WHERE
IF (n_o3s > 0) THEN
  WHERE (mask(:,:,:))
    tracer(:,:,:,n_o3s) = o33d(:,:,:)
  END WHERE
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_stratf
END MODULE ukca_stratf_mod
