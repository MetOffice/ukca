! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************

!  Description:
!   Module containing subroutine photol_calc_ozonecol.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################

MODULE photol_calc_ozonecol_mod

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus,                         &
                      PrStatus_Normal, PrStatus_Oper
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PHOTOL_CALC_OZONECOL_MOD'

CONTAINS

! ######################################################################
!------------------------------------------------------------------
! Subroutine CALC_OZONECOL
!------------------------------------------------------------------

! This routine calculates overhead ozone columns in molecules/cm^2 needed
! for photolysis calculations.

SUBROUTINE photol_calc_ozonecol(row_length, rows, model_levels,                &
  z_top_of_model, p_layer_boundaries,                                          &
  p_layer_centres, ozone_vmr, ozonecol)


USE ereport_mod,          ONLY: ereport

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine interface

! Model dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

REAL, INTENT(IN) :: z_top_of_model       ! model top (m)
REAL, INTENT(IN) :: p_layer_boundaries(row_length, rows,                       &
  0:model_levels)
REAL, INTENT(IN) :: p_layer_centres(row_length, rows, model_levels)
REAL, INTENT(IN) :: ozone_vmr(row_length, rows, model_levels)

REAL, INTENT(OUT) :: ozonecol(row_length, rows, model_levels)
! Ozone column above level in molecules/cm^2

! local variables

! Ozone column at model top. At 38 km, calculated from 60-level ozone clima
! tology. At 85 km, extrapolated (hence a wild guess but should not matter
! too much...).  In molecules/cm^2 units

REAL, PARAMETER :: ozcol_39km = 5.0e17
REAL, PARAMETER :: ozcol_41km = 4.8e17
REAL, PARAMETER :: ozcol_85km = 6.7e13

!REAL, PARAMETER :: colfac = 2.132e20 !! OLD VALUE
REAL, PARAMETER :: colfac = 2.117e20      ! Na/(g*M_air*1e4)

INTEGER :: l
INTEGER           :: errcode      ! Error code
CHARACTER(LEN=errormessagelength) :: cmessage     !   "   message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHOTOL_CALC_OZONECOL'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
ozonecol = 0.0
DO l=model_levels-1,1,-1

  ! compute the contributions from layers above l
  ozonecol(:,:,l) = ozonecol(:,:,l+1) +                                        &
    ozone_vmr(:,:,l+1) *                                                       &
    (p_layer_boundaries(:,:,l) - p_layer_boundaries(:,:,l+1))                  &
    * colfac
END DO

! add contribution within top of level L

DO l=1,model_levels-1
  ozonecol(:,:,l) = ozonecol(:,:,l) +                                          &
    ozone_vmr(:,:,model_levels) *                                              &
    (p_layer_centres(:,:,model_levels) -                                       &
    p_layer_boundaries(:,:,model_levels)) * colfac
END DO


! Add contribution above model top.
! If z_top_of_model lies in the given ranges, an approximate contribution
! of ozone column for above the model top (selected at 39, 41 or 85 km)
! is added to the total.

IF (z_top_of_model > 38000.0 .AND. z_top_of_model < 40500.0) THEN
  ozonecol = ozonecol + ozcol_39km
ELSE IF (z_top_of_model > 40500.0 .AND. z_top_of_model < 42000.0) THEN
  ozonecol = ozonecol + ozcol_41km
ELSE IF (z_top_of_model > 77000.0 .AND. z_top_of_model < 85500.0) THEN
  ozonecol = ozonecol + ozcol_85km
ELSE
  cmessage = ' Ozone column undefined for model top'
  errcode = 1
  CALL umPrint( cmessage,src='photol_calc_ozonecol_mod')
  WRITE(umMessage,'(A,E12.5)') 'z_top_of_model: ',z_top_of_model
  CALL umPrint(umMessage,src='photol_calc_ozonecol_mod')
  CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE photol_calc_ozonecol

END MODULE photol_calc_ozonecol_mod
