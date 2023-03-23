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
! Purpose: Subroutine to overwrite values at top of model using
!          external data
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
                       o33d, tracer)

USE ukca_config_specification_mod, ONLY: ukca_config

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

REAL, INTENT(IN) :: o33d(row_length,rows,model_levels) ! 3D O3 field

! Tracer concentrations in mass mixing ratio
REAL, INTENT(IN OUT) :: tracer(row_length,rows,model_levels,ntracer)

!     Local variables

INTEGER :: l            ! Loop variables

LOGICAL :: mask(row_length,rows,model_levels) ! mask to identify stratosphere

REAL, PARAMETER :: o3_hno3_ratio = 1.0/1000.0 ! kg[N]/kg[O3] from
                                              ! Murphy and Fahey 1994

REAL :: hno33d(row_length,rows,model_levels) ! 3D field from fixed o3:hno3 ratio

! Parameter to overwrite stratosphere (fixed no of levels above tropopause)
INTEGER, PARAMETER :: no_above_trop1 = 3        ! Suitable for L38/L60
INTEGER, PARAMETER :: no_above_trop2 = 10       ! Suitable for L63/L70/L85
INTEGER :: no_above_trop

INTEGER           :: errcode                    ! Error code: ereport
CHARACTER(LEN=errormessagelength) :: cmessage   ! Error message
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_STRATF'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set number of levels above tropopause at which to start overwriting values.
! If number of levels is not set explicitly then set it depending on vertical
! resolution. This is deprecated functionality designed to support specific UM
! grids. In future, the UM should set the number of levels explcitly via the
! ukca_setup call making this functionality redundant.
IF (ukca_config%nlev_above_trop_o3_env > 0) THEN
  no_above_trop = ukca_config%nlev_above_trop_o3_env
ELSE
  IF (model_levels == 38 .OR. model_levels == 60) THEN
    no_above_trop = no_above_trop1
  ELSE IF (model_levels == 63 .OR. model_levels == 70 .OR.                     &
           model_levels == 85) THEN
    no_above_trop = no_above_trop2
  ELSE
    errcode = 1
    cmessage =                                                                 &
      'Levels above tropopause must be set explicitly at this resolution'
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF

IF (printstatus == Prstatus_Diag) THEN
  WRITE(umMessage,'(A)') 'UKCA_STRATF:'
  CALL umPrint(umMessage,src='ukca_stratf')
  WRITE(umMessage,'(A,I0)') 'no_above_trop: ',no_above_trop
  CALL umPrint(umMessage,src='ukca_stratf')
  WRITE(umMessage,'(A)') ' '
  CALL umPrint(umMessage,src='ukca_stratf')
END IF

! Skip processing if there are no updates to do
IF (no_above_trop >= model_levels) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

hno33d(:,:,:) = o33d(:,:,:)*o3_hno3_ratio*c_hno3/c_n

! Overwrite o3, and hno3 at all gridboxes a fixed
!  number of model levels above the tropopause
!  O3   - external field
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
