! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Checks component and total average masses MD, MDT [per ptcl] and
!    number concentrations ND for each mode for bad values.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
! Subroutine Interface:
MODULE ukca_check_md_nd_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CHECK_MD_ND_MOD'

CONTAINS

SUBROUTINE ukca_check_md_nd(nbox,process,nd,md,mdt)
!-----------------------------------------------------------
!
!   Purpose
!   -------
!   Checks cpt and total average masses MD, MDT [per ptcl] and
!   number concentrations ND for each mode for bad values.
!
!   Bad is taken to be:
!
!   MD : greater than 10^20 or less than zero.
!   MDT: greater than 10^20 or less than or equal to zero.
!   ND : greater than 10^9 per cc or less than zero.
!
!   Inputs
!   ------
!   NBOX      : Number of grid boxes
!   PROCESS   : Character string with process just completed
!   ND        : Aerosol ptcl number density for mode (cm^-3)
!   MD        : Avg cpt mass of aerosol ptcl in size mode (particle^-1)
!   MDT       : Avg tot mass of aerosol ptcl in size mode (particle^-1)
!
!   Outputs
!   -------
!   None
!
!   Inputted by module UKCA_MODE_SETUP
!   ----------------------------------
!   NMODES    : Number of possible aerosol modes
!   NCP       : Number of possible aerosol components
!   MODE      : Logical variable defining which modes are set.
!   COMPONENT : Logical variable defining which cpt are in which dsts
!   MLO       : Lo-interf masses for initial radius grid
!-----------------------------------------------------------

USE ukca_mode_setup,   ONLY: nmodes, ncp, mode, component, mlo
USE parkind1,          ONLY: jprb, jpim
USE yomhook,           ONLY: lhook, dr_hook
USE ereport_mod,       ONLY: ereport
USE ukca_types_mod,    ONLY: logical32
USE umPrintMgr,        ONLY:                                                   &
    umPrint,                                                                   &
    umMessage,                                                                 &
    PrintStatus,                                                               &
    PrStatus_Oper

USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: nbox
REAL, INTENT(IN)    :: nd(nbox,nmodes)
REAL, INTENT(IN)    :: mdt(nbox,nmodes)
REAL, INTENT(IN)    :: md(nbox,nmodes,ncp)
CHARACTER(LEN=30), INTENT(IN) :: process

! Local variables
INTEGER :: jl
INTEGER :: imode
INTEGER :: icp
INTEGER :: toterr
INTEGER :: errcode
LOGICAL (KIND=logical32) :: mask1(nbox)
LOGICAL (KIND=logical32) :: mask2(nbox)
LOGICAL (KIND=logical32) :: mask12(nbox)
LOGICAL (KIND=logical32) :: mask3(nbox)
LOGICAL (KIND=logical32) :: mask4(nbox)
LOGICAL (KIND=logical32) :: mask5(nbox)
REAL    :: summd(nbox,nmodes)
REAL    :: mdtmin
REAL    :: mdtmax
REAL    :: ndmax
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHECK_MD_ND'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! .. also checks if the sum of MD over defined cpts is zero
summd(:,:)=0.0
DO imode=1,nmodes
  IF (mode(imode)) THEN
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        summd(:,imode)=summd(:,imode)+md(:,imode,icp)
      END IF
    END DO
  END IF
END DO

toterr=0
mdtmax=1.0e20
ndmax =1.0e9
DO imode=1,nmodes
  IF (mode(imode)) THEN
    !
    ! ..    set minimum allowed value of MDT
    mdtmin=mlo(imode)*0.001 ! MDTMIN equiv. to DPLIM0*0.1
    !
    mask1(:)=((nd (:,imode) <  ndmax) .AND.                                    &
              (nd (:,imode) >= 0.0e0))
    mask2(:)=((mdt(:,imode) <  mdtmax) .AND.                                   &
              (mdt(:,imode) >= mdtmin*0.001)) ! 0.01*DPLIM0 here
    mask12(:)=((.NOT. mask1(:)) .OR. (.NOT. mask2(:)))
    mask4(:)=(summd(:,imode) <= mdtmin*0.001) ! 0.01*DPLIM0 here
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        mask3(:)=((md(:,imode,icp) < mdtmax) .AND.                             &
                  (md(:,imode,icp) >= 0.0))
        mask5(:)=((mask12(:) .OR. (.NOT. mask3(:))) .OR. mask4(:))
        DO jl=1,nbox
          IF (mask5(jl)) THEN
            IF (PrintStatus >= PrStatus_Oper) THEN
              WRITE(umMessage,'(2A35)') ' Bad value in CHECK_MD_ND:  ',        &
                       'jl, imode, icp, nd, ndt, md, nd*md'
              CALL umPrint( process,src='ukca_check_md_nd')
              WRITE(umMessage,'(3i6,4e13.5)') jl,imode,icp,                    &
                nd(jl,imode),mdt(jl,imode),md(jl,imode,icp),                   &
                (nd(jl,imode)*md(jl,imode,icp))
              CALL umPrint(umMessage,src='ukca_check_md_nd')
            END IF
            toterr=toterr+1
          END IF ! if bad value for either ND,MDT or MD
        END DO ! loop over boxes
      END IF ! COMPONENT(IMODE,ICP)
    END DO
  END IF ! MODE(IMODE)
END DO ! loop over modes

IF (toterr > 0) THEN
  cmessage='Extreme values of ND, MDT, or MD found '
  IF (PrintStatus >= PrStatus_Oper) THEN
    WRITE(umMessage,'(A40,A10,I10)') cmessage,' Toterr: ',toterr
    CALL umPrint(umMessage,src='ukca_check_md_nd')
    errcode=-1
    CALL ereport('UKCA_CHECK_MD_ND',errcode,cmessage)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_check_md_nd
END MODULE ukca_check_md_nd_mod
