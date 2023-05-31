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
!    Calculate the change in number and mass concentration of
!    insoluble modes due to a combination of coagulation and nucleation.
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
MODULE ukca_coagwithnucl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_COAGWITHNUCL_MOD'

CONTAINS

SUBROUTINE ukca_coagwithnucl(nbox,nchemg,nbudaer,nd,md,mdt,delgc_nucl,dtz,     &
  ageterm2,intraoff,interoff,bud_aer_mas,kii_arr,kij_arr,                      &
  iagecoagnucl67,iextra_checks)

!---------------------------------------------------------------
!
! Purpose
! -------
! Calculate the change in number and mass concentration of
! insoluble modes due to a combination of coagulation and nucleation.
!
! Newly nucleated particles assumed to be 100 H2SO4 molecules
! per particle which equates to 5nm diameter.
!
! Coagulation equations are of the form dN/dt = a*N^2 + b*N + c
!
! Can rearrange this to be in the form of an indefinite integral
!
! int_{N_0}^{N} dx/X = int_{t_0}^{t_0+deltat} dt = deltat
!
! where X= A*x^2 + B*x + C
!
! Then solve this analytically (see header to ukca_solvecoagnucl_v.F90)
!
! Parameters
! ----------
! None
!
! Inputs
! ------
! NBOX    : Number of grid boxes
! nchemg  : Number of gas phase chemistry tracers
! nbudaer : Number of aerosol budget fields
! ND      : Initial no. concentration of aerosol mode (ptcls/cc)
! MD      : Initial avg cpt mass conc of aerosol mode (molecules/ptcl)
! MDT     : Initial avg total mass conc of aerosol mode (molecules/ptcl)
! delgc_nucl : Change in H2SO4 (g) due to nucleation (molecules/cc)
! DTZ     : nucl/coag/cond time step (s)
! INTRAOFF: Switch to turn off intra-modal coagulation
! INTEROFF: Switch to turn off intra-modal coagulation
! KII_ARR : Coag coeff for intra-modal coag (IMODE-IMODE) (cm^3/s)
! KIJ_ARR : Coag coeff for inter-modal coag (IMODE-JMODE) (cm^3/s)
! IAGECOAGNUCL67: Switch ageing,coag & nucle involving modes6&7 on (1)/off(0)
!
! Outputs
! -------
! ND      : Updated no. concentration of aerosol mode (ptcls/cc)
! MD      : Updated avg cpt mass conc of aerosol mode (molecules/ptcl)
! MDT     : Updated avg total mass conc of aerosol mode (molecules/ptcl)
! AGETERM2: Rate of accomodation of material to each insoluble mode
!           as a result of coagulation with smaller soluble modes
!           (in molecules cpt /cm3/DTZ)
!           Used for calculation of ageing rate in UKCA_AGEING.
! BUD_AER_MAS : Aerosol mass budgets
!
! Local variables
! ---------------
! RPI       : Geometric mean radius for mode IMODE (m)
! RPJ       : Geometric mean radius for mode JMODE (m)
! VPI       : Volume of particle of radius RPI (m^3)
! VPJ       : Volume of particle of radius RPJ (m^3)
! NDOLD     : Initial aerosol ptcl number conc (ptcls/cc)
! NDOLD_V   : Initial aerosol ptcl number conc (ptcls/cc) [1D version]
! MDOLD     : Initial aerosol cpt mass/ptcl (molecules/ptcl)
! MDCPOLD   : Initial aerosol cpt mass conc. (molecules/cm3)
! MDCPNEW   : New aerosol cpt mass conc. (molecules/cm3)
! MTRAN     : Stores transfer of cpt masses between modes by coag (/cm3)
! MTRANFMI  : Stores transfer of cpt masses from mode (/cm3)
! MTRANTOI  : Stores transfer of cpt masses to   mode (/cm3)
! MTRANNET  : Stores net transfer of cpt masses for mode (/cm3)
! DELN      : Change in ND (from NDOLD) due to coag & nucl (ptcls/cc)
! A         : Constant in nucl/coag equation
! B         : Constant in nucl/coag equation
! C         : Constant in nucl/coag equation
! BTERM     : Inter-modal coag rate term for single i-j
! KII       : Coag coeff for intra-modal coag (IMODE-IMODE) (cm^3/s)
! KIJ       : Coag coeff for inter-modal coag (IMODE-JMODE) (cm^3/s)
! XXX       : Term in exponential (1.0-EXP(-XXX(:)))
! XXX_EPS   : Tolerance for XXX below which don't evaluate exponential
! MASK1,MASK2,.. : Logicals to define domain regions for where loops
! TOPMODE   : Highest number mode for which coag & nucl is done
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! NMOL      : Number of molecules per particle at nucleation
! CONC_EPS  : Threshold for GC to calc nucl+coag (molecules per cc)
! DN_EPS    : Value of DELN below which do not carry out process
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! NCP       : Number of possible aerosol components
! MODE      : Defines which modes are set
! COMPONENT : Defines which cpts are allowed in each mode
! COAG_MODE : Defines which mode an IMODE-JMODE coagulation goes into
! MMID      : Mid-point masses for initial radius grid
! MFRAC_0   : Initial mass fraction to set when no particles.
! NUM_EPS   : Value of NEWN below which do not carry out process
! CP_SU     : Index of component in which sulfate is stored
! CP_BC     : Index of component in which BC is stored
! CP_OC     : Index of component in which OC is stored
! CP_CL     : Index of component in which NaCl is stored
! CP_DU     : Index of component in which dust is stored
! CP_SO     : Index of component in which condensible organic is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
USE ukca_constants,   ONLY: nmol, conc_eps, dn_eps

USE ukca_config_specification_mod, ONLY: glomap_variables

USE ukca_mode_setup,  ONLY: nmodes, coag_mode,                                 &
                            cp_su, cp_bc, cp_oc, cp_cl, cp_du, cp_so,          &
                            mode_nuc_sol, mode_ait_sol, mode_acc_sol,          &
                            mode_cor_sol, cp_no3, cp_nh4, cp_nn,               &
                            mode_ait_insol, mode_acc_insol, mode_cor_insol

USE ukca_setup_indices, ONLY: nmascoagsuintr12,                                &
         nmascoagocintr12, nmascoagsointr12, nmascoagsuintr13,                 &
         nmascoagocintr13, nmascoagsointr13, nmascoagsuintr14,                 &
         nmascoagocintr14, nmascoagsointr14, nmascoagsuintr15,                 &
         nmascoagocintr15, nmascoagsointr15, nmascoagsuintr16,                 &
         nmascoagocintr16, nmascoagsointr16, nmascoagsuintr17,                 &
         nmascoagocintr17, nmascoagsointr17, nmascoagsuintr23,                 &
         nmascoagbcintr23, nmascoagocintr23, nmascoagsointr23,                 &
         nmascoagsuintr24, nmascoagbcintr24, nmascoagocintr24,                 &
         nmascoagsointr24, nmascoagsuintr34, nmascoagbcintr34,                 &
         nmascoagocintr34, nmascoagssintr34, nmascoagsointr34,                 &
         nmascoagduintr34, nmascoagbcintr53, nmascoagocintr53,                 &
         nmascoagbcintr54, nmascoagocintr54, nmascoagduintr64,                 &
         nmascoagntintr23, nmascoagntintr24, nmascoagntintr34,                 &
         nmascoagnhintr23, nmascoagnhintr24, nmascoagnhintr34,                 &
         nmascoagnnintr34, mh2so4

USE ukca_mode_check_artefacts_mod, ONLY: ukca_mode_check_mdt
USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim
USE ukca_solvecoagnucl_v_mod, ONLY: ukca_solvecoagnucl_v
USE ukca_types_mod,   ONLY: logical_32

IMPLICIT NONE

! Subroutine interface:
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: nchemg
INTEGER, INTENT(IN) :: nbudaer
INTEGER, INTENT(IN) :: intraoff
INTEGER, INTENT(IN) :: interoff
INTEGER, INTENT(IN) :: iagecoagnucl67
INTEGER, INTENT(IN) :: iextra_checks

REAL, INTENT(IN)    :: delgc_nucl(nbox,nchemg)
REAL, INTENT(IN)    :: dtz
REAL, INTENT(IN)    :: kii_arr(nbox,nmodes)
REAL, INTENT(IN)    :: kij_arr(nbox,nmodes,nmodes)
REAL, INTENT(IN OUT) :: nd(nbox,nmodes)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables%ncp)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
REAL, INTENT(IN OUT) :: bud_aer_mas(nbox,0:nbudaer)
REAL, INTENT(OUT)   :: ageterm2(nbox,4,3,glomap_variables%ncp)

!  .. Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: mfrac_0(:,:)
REAL,    POINTER :: mmid(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: num_eps(:)


INTEGER :: imode
INTEGER :: jmode
INTEGER :: icp
INTEGER :: jcp
INTEGER :: topmode
LOGICAL (KIND=logical_32) :: mask1(nbox)
LOGICAL (KIND=logical_32) :: mask1a(nbox)
LOGICAL (KIND=logical_32) :: mask2(nbox)
LOGICAL (KIND=logical_32) :: mask3(nbox)
LOGICAL (KIND=logical_32) :: mask4(nbox)
REAL    :: ndold(nbox,nmodes)
REAL    :: ndold_v(nbox)
REAL    :: deln(nbox)
REAL    :: mdold(nbox,glomap_variables%ncp,nmodes)
REAL    :: mtran(nbox,glomap_variables%ncp,nmodes,nmodes)
REAL    :: mtranfmi(nbox,nmodes,glomap_variables%ncp)
REAL    :: mtrantoi(nbox,nmodes,glomap_variables%ncp)
REAL    :: mtrannet(nbox)
REAL    :: mdcpold(nbox)
REAL    :: mdcpnew(nbox)
REAL    :: a(nbox)
REAL    :: b(nbox)
REAL    :: c(nbox)
REAL    :: bterm(nbox)
REAL    :: kii(nbox)
REAL    :: kij(nbox)
REAL    :: xxx(nbox)
REAL    :: delh2so4_nucl(nbox)
REAL, PARAMETER :: xxx_eps=1.0e-3

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_COAGWITHNUCL'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mfrac_0     => glomap_variables%mfrac_0
mmid        => glomap_variables%mmid
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
num_eps     => glomap_variables%num_eps

ageterm2 = 0.0

!set limit of modes for which coagulation & nucleation occur
IF (iagecoagnucl67 == 1 ) THEN
  topmode = mode_cor_insol
ELSE
  topmode = mode_ait_insol
END IF


! Copy H2SO4 values from delgc_nucl to local variable
IF (mh2so4 > 0) THEN
  delh2so4_nucl(:)=delgc_nucl(:,mh2so4)
ELSE
  delh2so4_nucl(:)=0.0
END IF

bterm(:) = 0.0

DO imode=1,nmodes
  ndold(:,imode)=nd(:,imode)
  DO icp=1,ncp
    mdold(:,icp,imode)=md(:,imode,icp)
  END DO
END DO
mtran(:,:,:,:) = 0.0

DO imode = mode_nuc_sol, mode_cor_sol
  IF (mode(imode)) THEN

    kii(:)=kii_arr(:,imode) ! copy in pre-calculated KII

    a(:)=0.0
    b(:)=0.0
    c(:)=0.0

    ! Calculations are only done where NDOLD(:,IMODE) > NUM_EPS
    mask1(:) =(ndold(:,imode) > num_eps(imode))

    IF (intraoff /= 1) THEN
      WHERE (mask1(:)) a(:)=-0.5*kii(:)
    END IF

    DO jmode=(imode+1),mode_cor_sol ! inter-coag with larger soluble modes
      IF (mode(jmode)) THEN

        kij(:)=kij_arr(:,imode,jmode) ! copy in pre-calculated KIJ

        ! Calculations are only done where NDOLD(:,IMODE) > NUM_EPS
        ! and NDOLD(:,JMODE) > NUM_EPS(IMODE)
        mask2(:) =(ndold(:,imode) > num_eps(imode)) .AND.                      &
                  (ndold(:,jmode) > num_eps(jmode))

        IF (interoff /= 1) THEN
          WHERE (mask2(:))
            bterm(:)=-kij(:)*ndold(:,jmode)
            b(:)=b(:)+bterm(:)
          END WHERE
        END IF
        xxx(:)=-bterm(:)*dtz
        mask4(:)=ABS(xxx(:)) > xxx_eps
        ! .. above only evaluates exponential where it is "worth it"
        ! .. (in cases where XXX is larger than specified tolerance XXX_EPS)
        DO icp=1,ncp
          IF (component(jmode,icp)) THEN
            !
            WHERE (mask4(:) .AND. mask2(:))                                    &
              mtran(:,icp,imode,jmode)=                                        &
              mdold(:,icp,imode)*ndold(:,imode)*(1.0-EXP(-xxx(:)))
            !
            WHERE ((.NOT. mask4(:)) .AND. mask2(:))                            &
              mtran(:,icp,imode,jmode)=                                        &
              mdold(:,icp,imode)*ndold(:,imode)*xxx(:)
            !
            ! .. MTRAN(:,ICP,IMODE,JMODE) transfers from IMODE to JMODE
            !
          END IF
        END DO ! loop over components
      END IF ! if (MODE(JMODE))
    END DO ! end loop of JMODE over soluble modes

    DO jmode=(imode+4),topmode ! inter-coag with larger insoluble modes
      IF (mode(jmode)) THEN

        kij(:)=kij_arr(:,imode,jmode) ! copy in pre-calculated KIJ

        ! Calculations are only done where NDOLD(:,IMODE) > NUM_EPS
        ! and NDOLD(:,JMODE) > NUM_EPS
        mask2(:) =(ndold(:,imode) > num_eps(imode)) .AND.                      &
                  (ndold(:,jmode) > num_eps(jmode))

        IF (interoff /= 1) THEN
          WHERE (mask2(:))
            bterm(:)=-kij(:)*ndold(:,jmode)
            b(:)=b(:)+bterm(:)
          END WHERE
        END IF
        ! .. calculate MTRAN for sol-ins inter-modal coag to store in AGETERM2
        ! .. & carry out transfer of mass at end of the subroutine by MTRANNET.
        ! .. Also transfer number (include sol-ins term in B summation)
        xxx(:)=-bterm(:)*dtz
        mask4(:)=ABS(xxx(:)) > xxx_eps
        ! .. above only evaluates exponential where it is "worth it"
        ! .. (in cases where XXX is larger than specified tolerance XXX_EPS).
        DO icp=1,ncp
          WHERE (mask4(:) .AND. mask2(:))                                      &
      mtran(:,icp,imode,jmode)=                                                &
      mdold(:,icp,imode)*ndold(:,imode)*(1.0-EXP(-xxx(:)))
          WHERE ((.NOT. mask4(:)) .AND. mask2(:))                              &
      mtran(:,icp,imode,jmode)=                                                &
      mdold(:,icp,imode)*ndold(:,imode)*xxx(:)
        END DO
      END IF
    END DO ! end loop of JMODE over insoluble modes

    ! reset mass if insignificant # of particles
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        WHERE (.NOT. mask1(:))                                                 &
         md(:,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
      END IF
    END DO
    WHERE (.NOT. mask1(:)) mdt(:,imode)=mmid(imode)

    IF (imode == mode_nuc_sol) THEN

      WHERE (delh2so4_nucl(:) > conc_eps)
        c(:)=delh2so4_nucl(:)/dtz/nmol
      END WHERE ! if nucl mode & some particles are to be nucleated

    END IF

    mask1a(:)=mask1(:) .OR.                                                    &
              ((imode == mode_nuc_sol) .AND. (delh2so4_nucl(:) > conc_eps))

    ndold_v(:)=ndold(:,imode)
    CALL ukca_solvecoagnucl_v(nbox,mask1a,a,b,c,ndold_v,dtz,deln)

    WHERE (ABS(deln(:)) > dn_eps)
      nd(:,imode)=ndold(:,imode)+deln(:)
    END WHERE

  END IF ! if mode is defined
END DO ! end loop of IMODE over soluble modes
      ! (have updated ND here but not MD,MDT)
      ! (have stored transfer of mass from I->J in MTRAN(I,J)
      ! (updating of MD,MDT done at end of subroutine)
!
! .. Below is the section which has insoluble modes coagulating
! .. either with themselves (intra-modal) or inter-modally
! .. with larger soluble modes (inter-modal with ins neglected).
! .. Also the AGETERM2 is stored for the mass that coagulates
! .. with each insoluble mode from the smaller soluble modes.
!
! Solve for change in no. & mass due to intra and inter-modal coag here.
! Transfer of soluble mass from each soluble mode to each insoluble mode
! by inter-modal coag is stored here in AGETERM2 for use in UKCA_AGEING.
! Inter-modal coagulation with insoluble modes considered inefficient
! and neglected (as in M7).
!
! .. The section below stores the inter-modal coagulation between
! .. soluble modes and larger insoluble modes for calculation of ageing.
! .. Note that mass will be transferred to corresponding soluble mode
! .. (as prescribed by COAG_MODE) at end of coagulation subroutine
! .. because all soluble material is transferred in timestep by ageing
!
DO imode=mode_ait_insol,topmode
  IF (mode(imode)) THEN

    ageterm2(:,:,imode-4,:)=0.0

    DO jmode=mode_nuc_sol,mode_cor_sol ! loop over soluble modes
      IF (mode(jmode)) THEN
        DO jcp=1,ncp
          IF (component(jmode,jcp)) THEN
            ageterm2(:,jmode,imode-4,jcp)=mtran(:,jcp,jmode,imode)
          END IF
        END DO
      END IF
    END DO

    ! .. The section below calculates change in no. & mass due to
    ! .. intra-modal coag of insoluble modes and inter-modal coag with
    ! .. larger soluble modes (this does not contribute to AGETERM2 because
    ! .. it is automatically passed over when it coagulates with the larger
    ! .. soluble particles). Note that the number is updated here whereas
    ! .. the mass is updated at the end of the subroutine (stored in MTRAN).

    kii(:)=kii_arr(:,imode) ! copy in pre-calculated KII

    a(:)=0.0
    b(:)=0.0
    c(:)=0.0

    ! Calculations are only done where NDOLD(:,IMODE) > NUM_EPS
    mask1(:) =(ndold(:,imode) > num_eps(imode))

    IF (intraoff /= 1) THEN
      WHERE (mask1(:)) a(:)=-0.5*kii(:)
    END IF

    DO jmode=(imode-2),mode_cor_sol ! ins inter-coag with larger soluble modes
      IF (mode(jmode)) THEN

        kij(:)=kij_arr(:,imode,jmode) ! copy in pre-calculated KIJ

        ! Calculations are only done where NDOLD(:,IMODE) > NUM_EPS
        ! and NDOLD(:,JMODE) > NUM_EPS
        mask2(:) =(ndold(:,imode) > num_eps(imode)) .AND.                      &
                  (ndold(:,jmode) > num_eps(jmode))

        IF (interoff /= 1) THEN
          WHERE (mask2(:))
            bterm(:)=-kij(:)*ndold(:,jmode)
            b(:)=b(:)+bterm(:)
          END WHERE
        END IF
        ! .. calculate MTRAN for ins-sol inter-modal coag for carrying
        ! .. out transfer of mass at end of the subroutine by MTRANNET.
        ! .. Also transfer number (include ins-sol term in B summation)
        xxx(:)=-bterm(:)*dtz
        mask4(:)=ABS(xxx(:)) > xxx_eps
        ! .. above only evaluates exponential where it is "worth it"
        ! .. (in cases where XXX is larger than specified tolerance XXX_EPS)
        DO icp=1,ncp
          IF (component(imode,icp)) THEN
            WHERE (mask4(:) .AND. mask2(:))                                    &
             mtran(:,icp,imode,jmode)=                                         &
              mdold(:,icp,imode)*ndold(:,imode)*(1.0-EXP(-xxx(:)))
            WHERE ((.NOT. mask4(:)) .AND. mask2(:))                            &
             mtran(:,icp,imode,jmode)=                                         &
              mdold(:,icp,imode)*ndold(:,imode)*xxx(:)
          END IF
        END DO
      END IF
    END DO ! end loop of JMODE over larger soluble modes

    ndold_v(:)=ndold(:,imode)
    CALL ukca_solvecoagnucl_v(nbox,mask1,a,b,c,ndold_v,dtz,deln)

    WHERE (mask1(:) .AND. ABS(deln(:)) > dn_eps)
      nd(:,imode)=ndold(:,imode)+deln(:)
    END WHERE

    ! reset mass if insignificant # of particles
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        WHERE (.NOT. mask1(:))                                                 &
         md(:,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
      END IF
    END DO
    WHERE (.NOT. mask1(:)) mdt(:,imode)=mmid(imode)

  END IF ! end of IF(MODE(IMODE))
END DO ! end loop of IMODE over insoluble modes
!
DO icp=1,ncp
  DO imode=1,nmodes
    mtranfmi(:,imode,icp)=0.0
    mtrantoi(:,imode,icp)=0.0
  END DO
END DO
DO imode=1,nmodes
  IF (mode(imode)) THEN
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        DO jmode=1,nmodes
          IF (mode(jmode)) THEN
            mtranfmi(:,imode,icp)=mtranfmi(:,imode,icp)+                       &
                 mtran(:,icp,imode,jmode)
            mtrantoi(:,coag_mode(imode,jmode),icp)=                            &
                 mtrantoi(:,coag_mode(imode,jmode),icp)+                       &
                 mtran(:,icp,imode,jmode)
          END IF
        END DO
      END IF
    END DO
  END IF
END DO

DO imode=1,nmodes
  IF (mode(imode)) THEN

    ! Calculations are only done where ND(:,IMODE) > NUM_EPS
    mask1(:) =(nd(:,imode) > num_eps(imode))

    mdt(:,imode)=0.0
    mtrannet(:)=0.0
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        IF (imode == mode_nuc_sol) THEN
          IF (icp == cp_su) THEN ! add nucleated mass to sulfate cpt
            mtrannet(:)=mtrantoi(:,imode,icp)-                                 &
                        mtranfmi(:,imode,icp)+delh2so4_nucl(:)
          ELSE
            mtrannet(:)=mtrantoi(:,imode,icp)-mtranfmi(:,imode,icp)
          END IF
        ELSE
          mtrannet(:)=mtrantoi(:,imode,icp)-mtranfmi(:,imode,icp)
        END IF
        mdcpold(:)=ndold(:,imode)*mdold(:,icp,imode)
        mdcpnew(:)=mdcpold(:)+mtrannet(:)
        ! where MDCPNEW<0, set ND to zero (MDT,MD reset at end routine)
        WHERE (mask1(:) .AND. (mdcpnew(:) < 0.0))
          nd(:,imode)=0.0
          mask1(:)=.FALSE. ! set false so not used for other icp values
          ! n.b. MD and MDT will be re-set at end of routine
        END WHERE
        mask3(:)=mask1(:) .AND. (mdcpnew(:) >= 0.0)
        WHERE (mask3(:))
          md(:,imode,icp)=mdcpnew(:)/nd(:,imode)
          mdt(:,imode)=mdt(:,imode)+md(:,imode,icp)
        END WHERE
      END IF ! if COMPONENT(IMODE,ICP)
    END DO ! loop over ICP=1,NCP

    ! Apply extra checks on MDT being out of range after coagulation
    IF (iextra_checks > 1 ) THEN !do this only when ie_c = 2
      CALL ukca_mode_check_mdt(nbox, imode, mdt, md, nd, mask1)
    END IF

    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        mask3(:)=(mask1(:) .AND. (md(:,imode,icp) >= 0.0))
        DO jmode=1,nmodes
          IF (mode(jmode)) THEN
            IF ((imode == mode_nuc_sol) .AND. (jmode == mode_ait_sol)) THEN
              IF ((icp == cp_su) .AND. (nmascoagsuintr12 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsuintr12)=                             &
                  bud_aer_mas(:,nmascoagsuintr12)+mtran(:,cp_su,imode,jmode)
              END IF
              IF ((icp == cp_oc) .AND. (nmascoagocintr12 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagocintr12)=                             &
                  bud_aer_mas(:,nmascoagocintr12)+mtran(:,cp_oc,imode,jmode)
              END IF
              IF ((icp == cp_so) .AND. (nmascoagsointr12 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsointr12)=                             &
                  bud_aer_mas(:,nmascoagsointr12)+mtran(:,cp_so,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=1,2 (from mode 1 to mode 2)
            IF ((imode == mode_nuc_sol) .AND. (jmode == mode_acc_sol)) THEN
              IF ((icp == cp_su) .AND. (nmascoagsuintr13 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsuintr13)=                             &
                  bud_aer_mas(:,nmascoagsuintr13)+mtran(:,cp_su,imode,jmode)
              END IF
              IF ((icp == cp_oc) .AND. (nmascoagocintr13 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagocintr13)=                             &
                  bud_aer_mas(:,nmascoagocintr13)+mtran(:,cp_oc,imode,jmode)
              END IF
              IF ((icp == cp_so) .AND. (nmascoagsointr13 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsointr13)=                             &
                  bud_aer_mas(:,nmascoagsointr13)+mtran(:,cp_so,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=1,3 (from mode 1 to mode 3)
            IF ((imode == mode_nuc_sol) .AND. (jmode == mode_cor_sol)) THEN
              IF ((icp == cp_su) .AND. (nmascoagsuintr14 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsuintr14)=                             &
                  bud_aer_mas(:,nmascoagsuintr14)+mtran(:,cp_su,imode,jmode)
              END IF
              IF ((icp == cp_oc) .AND. (nmascoagocintr14 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagocintr14)=                             &
                  bud_aer_mas(:,nmascoagocintr14)+mtran(:,cp_oc,imode,jmode)
              END IF
              IF ((icp == cp_so) .AND. (nmascoagsointr14 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsointr14)=                             &
                  bud_aer_mas(:,nmascoagsointr14)+mtran(:,cp_so,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=1,4 (from mode 1 to mode 4)
            IF ((imode == mode_nuc_sol) .AND. (jmode == mode_ait_insol)) THEN
              IF ((icp == cp_su) .AND. (nmascoagsuintr15 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsuintr15)=                             &
                  bud_aer_mas(:,nmascoagsuintr15)+mtran(:,cp_su,imode,jmode)
              END IF
              IF ((icp == cp_oc) .AND. (nmascoagocintr15 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagocintr15)=                             &
                  bud_aer_mas(:,nmascoagocintr15)+mtran(:,cp_oc,imode,jmode)
              END IF
              IF ((icp == cp_so) .AND. (nmascoagsointr15 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsointr15)=                             &
                  bud_aer_mas(:,nmascoagsointr15)+mtran(:,cp_so,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=1,5 (from mode 1 to mode 5)
            IF ((imode == mode_nuc_sol) .AND. (jmode == mode_acc_insol)) THEN
              IF ((icp == cp_su) .AND. (nmascoagsuintr16 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsuintr16)=                             &
                  bud_aer_mas(:,nmascoagsuintr16)+mtran(:,cp_su,imode,jmode)
              END IF
              IF ((icp == cp_oc) .AND. (nmascoagocintr16 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagocintr16)=                             &
                  bud_aer_mas(:,nmascoagocintr16)+mtran(:,cp_oc,imode,jmode)
              END IF
              IF ((icp == cp_so) .AND. (nmascoagsointr16 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsointr16)=                             &
                  bud_aer_mas(:,nmascoagsointr16)+mtran(:,cp_so,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=1,6 (from mode 1 to mode 6)
            IF ((imode == mode_nuc_sol) .AND. (jmode == mode_cor_insol)) THEN
              IF ((icp == cp_su) .AND. (nmascoagsuintr17 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsuintr17)=                             &
                  bud_aer_mas(:,nmascoagsuintr17)+mtran(:,cp_su,imode,jmode)
              END IF
              IF ((icp == cp_oc) .AND. (nmascoagocintr17 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagocintr17)=                             &
                  bud_aer_mas(:,nmascoagocintr17)+mtran(:,cp_oc,imode,jmode)
              END IF
              IF ((icp == cp_so) .AND. (nmascoagsointr17 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsointr17)=                             &
                  bud_aer_mas(:,nmascoagsointr17)+mtran(:,cp_so,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=1,7 (from mode 1 to mode 7)
            IF ((imode == mode_ait_sol) .AND. (jmode == mode_acc_sol)) THEN
              IF ((icp == cp_su) .AND. (nmascoagsuintr23 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsuintr23)=                             &
                  bud_aer_mas(:,nmascoagsuintr23)+mtran(:,cp_su,imode,jmode)
              END IF
              IF ((icp == cp_bc) .AND. (nmascoagbcintr23 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagbcintr23)=                             &
                  bud_aer_mas(:,nmascoagbcintr23)+mtran(:,cp_bc,imode,jmode)
              END IF
              IF ((icp == cp_oc) .AND. (nmascoagocintr23 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagocintr23)=                             &
                  bud_aer_mas(:,nmascoagocintr23)+mtran(:,cp_oc,imode,jmode)
              END IF
              IF ((icp == cp_so) .AND. (nmascoagsointr23 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsointr23)=                             &
                  bud_aer_mas(:,nmascoagsointr23)+mtran(:,cp_so,imode,jmode)
              END IF
              IF ((icp == cp_no3) .AND. (nmascoagntintr23 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagntintr23)=                             &
                  bud_aer_mas(:,nmascoagntintr23)+mtran(:,cp_no3,imode,jmode)
              END IF
              IF ((icp == cp_nh4) .AND. (nmascoagnhintr23 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagnhintr23)=                             &
                  bud_aer_mas(:,nmascoagnhintr23)+mtran(:,cp_nh4,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=2,3 (from mode 2 to mode 3)
            IF ((imode == mode_ait_sol) .AND. (jmode == mode_cor_sol)) THEN
              IF ((icp == cp_su) .AND. (nmascoagsuintr24 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsuintr24)=                             &
                  bud_aer_mas(:,nmascoagsuintr24)+mtran(:,cp_su,imode,jmode)
              END IF
              IF ((icp == cp_bc) .AND. (nmascoagbcintr24 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagbcintr24)=                             &
                  bud_aer_mas(:,nmascoagbcintr24)+mtran(:,cp_bc,imode,jmode)
              END IF
              IF ((icp == cp_oc) .AND. (nmascoagocintr24 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagocintr24)=                             &
                  bud_aer_mas(:,nmascoagocintr24)+mtran(:,cp_oc,imode,jmode)
              END IF
              IF ((icp == cp_so) .AND. (nmascoagsointr24 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsointr24)=                             &
                  bud_aer_mas(:,nmascoagsointr24)+mtran(:,cp_so,imode,jmode)
              END IF
              IF ((icp == cp_no3) .AND. (nmascoagntintr24 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagntintr24)=                             &
                  bud_aer_mas(:,nmascoagntintr24)+mtran(:,cp_no3,imode,jmode)
              END IF
              IF ((icp == cp_nh4) .AND. (nmascoagnhintr24 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagnhintr24)=                             &
                  bud_aer_mas(:,nmascoagnhintr24)+mtran(:,cp_nh4,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=2,4 (from mode 2 to mode 4)
            IF ((imode == mode_acc_sol) .AND. (jmode == mode_cor_sol)) THEN
              IF ((icp == cp_su) .AND. (nmascoagsuintr34 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsuintr34)=                             &
                  bud_aer_mas(:,nmascoagsuintr34)+mtran(:,cp_su,imode,jmode)
              END IF
              IF ((icp == cp_bc) .AND. (nmascoagbcintr34 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagbcintr34)=                             &
                  bud_aer_mas(:,nmascoagbcintr34)+mtran(:,cp_bc,imode,jmode)
              END IF
              IF ((icp == cp_oc) .AND. (nmascoagocintr34 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagocintr34)=                             &
                  bud_aer_mas(:,nmascoagocintr34)+mtran(:,cp_oc,imode,jmode)
              END IF
              IF ((icp == cp_cl) .AND. (nmascoagssintr34 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagssintr34)=                             &
                  bud_aer_mas(:,nmascoagssintr34)+mtran(:,cp_cl,imode,jmode)
              END IF
              IF ((icp == cp_so) .AND. (nmascoagsointr34 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagsointr34)=                             &
                  bud_aer_mas(:,nmascoagsointr34)+mtran(:,cp_so,imode,jmode)
              END IF
              IF ((icp == cp_du) .AND. (nmascoagduintr34 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagduintr34)=                             &
                  bud_aer_mas(:,nmascoagduintr34)+mtran(:,cp_du,imode,jmode)
              END IF
              IF ((icp == cp_no3) .AND. (nmascoagntintr34 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagntintr34)=                             &
                  bud_aer_mas(:,nmascoagntintr34)+mtran(:,cp_no3,imode,jmode)
              END IF
              IF ((icp == cp_nh4) .AND. (nmascoagnhintr34 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagnhintr34)=                             &
                  bud_aer_mas(:,nmascoagnhintr34)+mtran(:,cp_nh4,imode,jmode)
              END IF
              IF ((icp == cp_nn) .AND. (nmascoagnnintr34 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagnnintr34)=                             &
                  bud_aer_mas(:,nmascoagnnintr34)+mtran(:,cp_nn,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=3,4 (from mode 3 to mode 4)
            IF ((imode == mode_ait_insol) .AND. (jmode == mode_acc_sol)) THEN
              IF ((icp == cp_bc) .AND. (nmascoagbcintr53 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagbcintr53)=                             &
                  bud_aer_mas(:,nmascoagbcintr53)+mtran(:,cp_bc,imode,jmode)
              END IF
              IF ((icp == cp_oc) .AND. (nmascoagocintr53 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagocintr53)=                             &
                  bud_aer_mas(:,nmascoagocintr53)+mtran(:,cp_oc,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=5,3 (from mode 5 to mode 3)
            IF ((imode == mode_ait_insol) .AND. (jmode == mode_cor_sol)) THEN
              IF ((icp == cp_bc) .AND. (nmascoagbcintr54 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagbcintr54)=                             &
                  bud_aer_mas(:,nmascoagbcintr54)+mtran(:,cp_bc,imode,jmode)
              END IF
              IF ((icp == cp_oc) .AND. (nmascoagocintr54 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagocintr54)=                             &
                  bud_aer_mas(:,nmascoagocintr54)+mtran(:,cp_oc,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=5,4 (from mode 5 to mode 4)
            IF ((imode == mode_acc_insol) .AND. (jmode == mode_cor_sol)) THEN
              IF ((icp == cp_du) .AND. (nmascoagduintr64 > 0)) THEN
                WHERE (mask3(:))                                               &
                  bud_aer_mas(:,nmascoagduintr64)=                             &
                  bud_aer_mas(:,nmascoagduintr64)+mtran(:,cp_du,imode,jmode)
              END IF
            END IF ! IF IMODE,JMODE=6,4 (from mode 6 to mode 4)
          END IF ! IF MODE(JMODE)
        END DO ! LOOP JMODE (those that may have been transferred
              !             to mode JMODE from mode IMODE)
      END IF ! IF COMPONENT(IMODE,ICP)
    END DO ! LOOP OVER COMPONENTS
    !
            ! reset mass if insignificant # of particles
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        WHERE (.NOT. mask1(:))                                                 &
         md(:,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
      END IF
    END DO
    WHERE (.NOT. mask1(:)) mdt(:,imode)=mmid(imode)
    !
  END IF ! IF(MODE(IMODE)
END DO ! LOOP IMODE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_coagwithnucl
END MODULE ukca_coagwithnucl_mod
