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
!   Calculates geometric mean dry diameter for multi-component
!   aerosol population which is lognormally distributed with
!   number concentration ND in mode, component mass concentration
!   MD in mode, component density RHOCOMP and component molecular mass MM
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
MODULE ukca_calc_drydiam_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CALC_DRYDIAM_MOD'

CONTAINS

SUBROUTINE ukca_calc_drydiam( nbox, glomap_variables_local,                    &
                              nd, md, mdt, drydp, dvol )

!-------------------------------------------------------------
!
!     Calculates geometric mean dry diameter for multi-component
!     aerosol population which is lognormally distributed with
!     number concentration ND in mode,
!     cpt mass concentration MD in mode,
!     cpt density RHOCOMP and cpt molecular mass MM
!
!     Calculate dry volume per particle using composition info as:
!
!     DVOL = sum_cp { MD(ICP)*MM(ICP)/(AVOGADRO*RHOCOMP(ICP)) }
!
!     Where AVOGADRO is Avogadro's constant. Then, from Jacobsen,
!     "Fundamentals of Atmospheric Modeling", pp. 412, have
!
!     dry volume conc. = ND*(PI/6)*(Dp^3)*exp{9/2 log^2(sigma_g)}
!
!     i.e. DVOL  = (PI/6)*(Dp^3)*exp{9/2 log^2(sigma_g)}
!
!     where Dp is the number mean dry diameter,
!     and sigma_g is the geometric standard deviation.
!
!     Then calculate Dp as:
!
!     Dp=CBRT( DVOL*(6/PI)/EXP{9.2 log^2(sigma_g)} )
!
!     Inputs
!     ------
!     NBOX      : Number of grid boxes
!     ND        : Aerosol ptcl no. concentration (ptcls per cc)
!     MD        : Component median aerosol mass (molecules per ptcl)
!     MDT       : Total median aerosol mass (molecules per ptcl)
!
!     Outputs
!     -------
!     DRYDP     : Median particle dry diameter for each mode (m)
!     DVOL      : Median particle dry volume for each mode (m^3)
!     MD        : Component median aerosol mass (molecules per ptcl)
!     MDT       : Total median aerosol mass (molecules per ptcl)
!
!     Local Variables
!     ---------------
!     None
!
!     Inputted by modules CHEMISTRY_CONSTANT, UKCA_CONSTANTS
!     ---------------------------------
!     AVOGADRO  : Avogadro's constant (per mole)
!     MMSUL     : Molar mass of a pure H2SO4 aerosol (kg per mole)
!     RHO_SO4    : Mass density of a pure H2SO4 aerosol (kg per m^3)
!
!     Inputted via argument glomap_variables_local
!     ----------------------------------
!     NMODES    : Number of aerosol modes
!     NCP       : Number of aerosol components
!     MM        : Molecular mass of each component
!     RHOCOMP   : Densities (dry) of each component (kg/m^3)
!     MODE      : Which modes are being carried
!     COMPONENT : Which components are in each of modes
!     DDPLIM0   : Lower limits for dry diameter for each mode (m)
!     MFRAC_0   : Initial mass fraction to set when no particles.
!     X_EQUATION: EXP((9/2)*LOG^2(SIGMA_G))
!     NUM_EPS   : Value of NEWN below which do not carry out process
!     MMID      : Mid-point masses for initial radius grid
!     MLO       : Lo-interf masses for initial radius grid
!
!--------------------------------------------------------------------
USE ukca_config_constants_mod, ONLY: avogadro, rho_so4
USE ukca_constants,  ONLY: pi, mmsul
USE ukca_um_legacy_mod, ONLY: cubrt_v
USE ukca_types_mod, ONLY: log_small

USE ukca_mode_setup, ONLY:                                                     &
    glomap_variables_type,                                                     &
    nmodes,                                                                    &
    mode_nuc_sol,                                                              &
    mode_acc_sol

USE yomhook,         ONLY: lhook, dr_hook
USE parkind1,        ONLY: jprb, jpim
USE ereport_mod,     ONLY: ereport
USE umPrintMgr,      ONLY: umPrint, umMessage


USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Interface
INTEGER, INTENT(IN) :: nbox
TYPE(glomap_variables_type), TARGET, INTENT(IN) :: glomap_variables_local
REAL, INTENT(IN)    :: nd(nbox,nmodes)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables_local%ncp)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
REAL, INTENT(OUT)   :: drydp(nbox,nmodes)
REAL, INTENT(OUT)   :: dvol(nbox,nmodes)

! Local variables

! Caution - pointers to TYPE glomap_variables_local%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: ddplim0(:)
REAL,    POINTER :: mfrac_0(:,:)
REAL,    POINTER :: mlo(:)
REAL,    POINTER :: mm(:)
REAL,    POINTER :: mmid(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: num_eps(:)
REAL,    POINTER :: rhocomp(:)
REAL,    POINTER :: x_equation(:)

INTEGER :: jl
INTEGER :: imode
INTEGER :: icp
REAL    :: ddpcub(nbox,nmodes)
REAL    :: sixovrpix(nmodes)
LOGICAL (KIND=log_small) :: mask(nbox)
REAL    :: ratio1(glomap_variables_local%ncp)
REAL    :: ratio2
REAL    :: dp_thresh1
REAL    :: dp
REAL :: min_dvol, min_drydp
INTEGER :: i
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER           :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALC_DRYDIAM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables_local%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables_local%component
ddplim0     => glomap_variables_local%ddplim0
mfrac_0     => glomap_variables_local%mfrac_0
mlo         => glomap_variables_local%mlo
mm          => glomap_variables_local%mm
mmid        => glomap_variables_local%mmid
mode        => glomap_variables_local%mode
ncp         => glomap_variables_local%ncp
num_eps     => glomap_variables_local%num_eps
rhocomp     => glomap_variables_local%rhocomp
x_equation  => glomap_variables_local%x

! Below is over ncp
ratio1(:)=mm(:)/(avogadro*rhocomp(:))
!
! Below is over NMODES
!
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                               &
!$OMP PRIVATE(i, icp, imode, mask, ratio2)                                     &
!$OMP SHARED(avogadro, component, ddpcub, dvol, md, mmid, mode, nbox, ncp, nd, &
!$OMP num_eps, ratio1, rho_so4, sixovrpix, x_equation)
DO imode=1,nmodes
  IF (mode(imode)) THEN
    DO i=1,nbox
      mask(i) =(nd(i,imode) > num_eps(imode))
      IF (mask(i)) THEN
        dvol(i,imode)=0.0
      ELSE
        dvol(i,imode)=mmid(imode)*mmsul/(avogadro*rho_so4)
      END IF
    END DO
    !     calculate particle dry volume using composition info
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        DO i=1,nbox
          IF (mask(i)) THEN
            dvol(i,imode)=dvol(i,imode)+ratio1(icp)*md(i,imode,icp)
          END IF
        END DO
      END IF
    END DO
    !     DVOL calculates particle dry volume assuming pure H2SO4
  ELSE
    ratio2=mmsul*mmid(imode)/(avogadro*rho_so4)
    DO i=1,nbox
      dvol(i,imode)=ratio2
    END DO
  END IF
  sixovrpix(imode)=6.0/(pi*x_equation(imode))
  DO i=1,nbox
    ddpcub(i,imode)=sixovrpix(imode)*dvol(i,imode)
  END DO
END DO
!$OMP END PARALLEL DO

CALL cubrt_v(nmodes*nbox, ddpcub, drydp)

! .. also check whether mean diameter too low for mode
! only do check for solvent modes nuc,ait,acc
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                               &
!$OMP PRIVATE(dp, dp_thresh1, icp, imode, jl)                                  &
!$OMP SHARED(avogadro, component, ddplim0, drydp, dvol, md, mdt, mfrac_0,      &
!$OMP mlo, mode, nbox, ncp, rho_so4, sixovrpix)
DO imode=mode_nuc_sol,mode_acc_sol
  IF (mode(imode)) THEN
    DO jl=1,nbox
      dp    =drydp(jl,imode)
      dp_thresh1=ddplim0(imode)*0.1
      IF (dp < dp_thresh1) THEN
        DO icp=1,ncp
          IF (component(imode,icp)) THEN
            md(jl,imode,icp)=mlo(imode)*mfrac_0(imode,icp)
          END IF
        END DO
        mdt(jl,imode)=mlo(imode)
        dvol(jl,imode)=mlo(imode)*mmsul/(avogadro*rho_so4)
        drydp(jl,imode)=(sixovrpix(imode)*dvol(jl,imode))**(1.0/3.0)
      END IF
    END DO
  END IF
END DO
!$OMP END PARALLEL DO

DO imode=1,nmodes
  min_dvol=MINVAL(dvol (:,imode))
  min_drydp=MINVAL(drydp (:,imode))
  IF ((min_dvol               <= 0.0)  .OR.                                    &
      (min_drydp              <= 0.0)) THEN

    WRITE(umMessage,'(A,I0)') 'In calcdrydiam: drydp (min,max,sum) imode=',imode
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E20.8)') 'MINVAL(drydp(:,imode)) = ',                  &
                                  min_drydp
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E20.8)') 'MAXVAL(drydp(:,imode)) = ',                  &
                                  MAXVAL(drydp(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E20.8)') 'SUM(drydp(:,imode)) = ',                     &
                                  SUM(drydp(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'Location of min: ',MINLOC(drydp(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'Location of max: ',MAXLOC(drydp(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'In calcdrydiam: dvol (min,max,sum) imode=',imode
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E20.8)') 'MINVAL(dvol(:,imode)) = ',                   &
                                  min_dvol
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E20.8)') 'MAXVAL(dvol(:,imode)) = ',                   &
                                  MAXVAL(dvol(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,E20.8)') 'SUM(dvol(:,imode)) = ',                      &
                                  SUM(dvol(:,imode))
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'Location of min: ',MINLOC(drydp(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    WRITE(umMessage,'(A,I0)') 'Location of max: ',MAXLOC(drydp(:,imode), DIM=1)
    CALL umPrint(umMessage,src=RoutineName)

    cmessage = ' dvol or drydp <= 0'
    errcode=1
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calc_drydiam
END MODULE ukca_calc_drydiam_mod
