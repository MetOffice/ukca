! *****************************COPYRIGHT*******************************

! (c) [University of Leeds] [2013]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************

!  Description:
!    Module for routine UKCA_MODE_CHECK_ARTEFACTS

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and
!  The Met Office. See:  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  FORTRAN 90

! ######################################################################

MODULE ukca_mode_check_artefacts_mod
! ---------------------------------------------------------------------|
!+ Module to contain subroutine

! Description:
! To check and correct for artefacts caused by advection of modal tracers

! Note: currently code is hard-coded so that ordering of modes must
! 1) nucln, 2)   soluble Aitken, 3)   soluble accum, 4)   soluble coarse
!           5) insoluble Aitken, 6) insoluble accum, 7) insoluble coarse
!           8) insoluble super-coarse

! Contains subroutines:
!       1) ukca_mode_check_artefacts   .. called from ukca_aero_ctl
!       2) ukca_mode_check_mdt         .. called from ukca_coagwithnucl


! ---------------------------------------------------------------------|
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER,                                                   &
              PRIVATE :: ModuleName='UKCA_MODE_CHECK_ARTEFACTS_MOD'

REAL, PARAMETER :: mdtmin_par = 10 ! Absolute min MDT for a particle

CONTAINS

! #############################################################################

! Subroutine Interface:
SUBROUTINE ukca_mode_check_artefacts(verbose, nbox, nd, md, mdt,               &
                                  sm, aird, mm_da, dtc, mdtfixflag, mdtfixsink)

! Check to see if mdt is outside limits and reset nd and md to defaults if
! it is.
! Called from UKCA_AERO_CTL.

USE ukca_config_specification_mod, ONLY: glomap_variables

USE ukca_mode_setup,    ONLY:                                                  &
    nmodes,                                                                    &
    mode_nuc_sol,                                                              &
    mode_ait_sol,                                                              &
    mode_acc_sol,                                                              &
    mode_cor_sol,                                                              &
    mode_ait_insol,                                                            &
    mode_acc_insol,                                                            &
    mode_cor_insol,                                                            &
    mode_sup_insol

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umMessage

IMPLICIT NONE

INTEGER, INTENT(IN) :: verbose
! Switch to determine level of debug o/p (0=none, 1, 2)
INTEGER, INTENT(IN) :: nbox           ! No. of points

REAL, INTENT(IN)    :: sm(nbox)   ! Air mass (kg per gridbox)
REAL, INTENT(IN)    :: aird(nbox) ! Molecular air density (molecules per cm^3)
REAL, INTENT(IN)    :: mm_da      ! Molar mass of dry air (kg per mole)
REAL, INTENT(IN)    :: dtc        ! Time step on which GLOMAP is integrated

REAL, INTENT(IN OUT) :: mdtfixflag(nbox,nmodes)
! Arrays storing info about where ND->0 when MDT is out of range.
REAL, INTENT(IN OUT) :: mdtfixsink(nbox,nmodes)
! Arrays storing info about how much mass is removed when ND>0 for o-o-r MDT
! mdtfixflag gets set = to 100.0 when fix applied so that, by looking at the
!    climate-means see what average percentage of timesteps the fix is applied.
! mdtfixsink stores the amount of *total mass" (over all components) that is
!         removed when the fix is applied to enable to monitor that if required.
! Both diags stored in section 38 (temporarily overwrite available budget diags)

REAL, INTENT(IN OUT) :: nd(nbox,nmodes)
! Each mode's particle number density (ptcls cm^-3)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables%ncp)
! Each mode's average component mass per particle (molecules ptcle^-1)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
! Each mode's average total mass per particle (over all cpts) molecules ptcl^-1)

! Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: mfrac_0 (:,:)
REAL,    POINTER :: mhi (:)
REAL,    POINTER :: mlo (:)
REAL,    POINTER :: mmid (:)
LOGICAL, POINTER :: mode (:)
INTEGER, POINTER :: ncp


INTEGER :: nbadmdt(nbox,nmodes)           ! Count number of bad MDT points
INTEGER :: icp                            ! Counter for loop over components
INTEGER :: imode                          ! Counter for loop over modes

REAL    :: fac(nbox) ! Conversion factor for budget diagnostics

REAL    :: mdtmin                  ! min MDT value for each mode
REAL    :: mdtmax                  ! max MDT value for each mode

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_CHECK_ARTEFACTS'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mfrac_0     => glomap_variables%mfrac_0
mhi         => glomap_variables%mhi
mlo         => glomap_variables%mlo
mmid        => glomap_variables%mmid
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp

IF (verbose > 1) THEN
  !  Check min, max values of ND, MDT and MD
  WRITE(umMessage,'(A50)') 'BEFORE CALL TO UKCA_AERO_STEP (before fix)'
  DO imode=1,nmodes
    IF (mode(imode)) THEN
      WRITE(umMessage,'(A10,I4,2E12.4)') 'ND : ',imode,MINVAL(nd(:,imode)),    &
        MAXVAL(nd(:,imode))
      WRITE(umMessage,'(A10,I4,2E12.4)') 'MDT: ',imode,MINVAL(mdt(:,imode)),   &
        MAXVAL(mdt(:,imode))
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WRITE(umMessage,'(A10,I4,2E12.4)') 'MD : ',imode,                    &
            MINVAL(md(:,imode,icp)),MAXVAL(md(:,imode,icp))
        END IF
      END DO
    END IF
  END DO
END IF ! if VERBOSE > 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! "MDT checks" code for separate-moment tracer artefacts
!!
!! If outside MDTMIN/MDTMAX range then set ND-->0 -- in this case *must*
!! be artefact from separate advection of tracer moments and enforced
!! conservation or different re-scaling applied to one of other of the
!! moments (number or mass) which lead to values going out of valid range.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Original checks retained for reference
!!        MDTMIN=MLO(IMODE)*0.001  ! set equiv. to DPLIM0*0.1
!!        MDTMAX=MHI(IMODE)*1000.0 ! set equiv. to DPLIM1*10.0

DO imode=1,nmodes
  IF (mode(imode)) THEN
    ! below checks if MDT is coming out too low after advection
    ! 1) Min MDT criteria
    ! -------------------
    !  1.1 Soluble modes (all) : min MDT factor-10 lower in diam than
    !      min of mode
    IF ((imode == mode_nuc_sol) .OR. (imode == mode_ait_sol) .OR.              &
        (imode == mode_acc_sol) .OR. (imode == mode_cor_sol)) THEN
      mdtmin = mlo(imode)*0.001
    END IF
    !  1.2 Insoluble Aitken mode : min MDT factor-10 lower in diam than
    !      min of mode
    IF (imode == mode_ait_insol) mdtmin = mlo(imode)*0.001
    !  1.3 Insoluble accumn mode : only minimal constraint: OK down
    !      to mlo(1)*0.001
    !      only minimal constraint on mode 6
    !      (where contains mainly MSP can be v small)
    IF (imode == mode_acc_insol) mdtmin = mlo(1)*0.001
    !  1.4 Insoluble coarse mode : min MDT factor-10 lower in diam than
    !      min of mode
    IF (imode == mode_cor_insol) mdtmin = mlo(imode)*0.001
    !  1.5 Insoluble super-coarse mode: min MDT factor-10 lower in diam
    !      than min of mode
    IF (imode == mode_sup_insol) mdtmin = mlo(imode)*0.001

    ! 2) Max MDT criteria
    ! -------------------
    ! All modes have same maxMDT: max MDT factor-10 higher in diam than
    ! max of mode
    mdtmax = mhi(imode)*1000.0

    ! 3) Don't allow MDT to be below minimum possible particle size
    ! -------------------------------------------------------------
    IF (mdtmin < mdtmin_par) mdtmin = mdtmin_par

    ! 4) Where MDT outside (mdtmin,mdtmax) range re-set MD to default
    ! mid-pt value
    ! ---------------------------------------------------------------
    ! n.b.will also set ND-->0 later so these MD values are nominal
    ! (but needed).
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        ! where MDT too low after advection set ND to zero and set default MD
        WHERE (mdt(:,imode) < mdtmin)
          md(:,imode,icp) = mmid(imode)*mfrac_0(imode,icp)
        END WHERE
        ! where MDT too high after advection set ND-->zero & MD-->default mid-pt
        WHERE (mdt(:,imode) > mdtmax)
          md(:,imode,icp) = mmid(imode)*mfrac_0(imode,icp)
        END WHERE
      END IF
    END DO

    ! 5a) Initialise nbadmdt MDTFIX arrays to zero
    ! --------------------------------------------
    nbadmdt(:,imode)=0
    mdtfixflag(:,imode) = 0.0
    mdtfixsink(:,imode) = 0.0
    fac(:)=sm(:)/aird(:)
    ! FAC converts aerosol mass fluxes from kg(dryair)/box/tstep
    ! to moles/gridbox/s

    ! 5b) count occurrences of MDT<MDTMIN (NBADMDT) and set 100% where occurs
    ! -----------------------------------------------------------------------
    WHERE (mdt(:,imode) < mdtmin)
      nbadmdt(:,imode) = 1
      mdtfixflag(:,imode)=100.0
      ! mdtfixflag enables to track proportion of timesteps that fix is applied
      mdtfixsink(:,imode)=nd(:,imode)*mdt(:,imode)*fac(:)/mm_da/dtc ! moles/s
      ! mdtfixsink stores total mass removed when fix is applied
    END WHERE

    ! 6) print out number of MDT-too-low occurrences if > 0
    IF (SUM(nbadmdt(:,imode)) > 0 .AND. verbose > 0) THEN
      WRITE(umMessage,'(A39)') 'MDT<MDTMIN, ND=0: IMODE,MDTMIN,NBADMDT:'
      WRITE(umMessage,'(I6,1E12.3,I12)') imode,mdtmin,SUM(nbadmdt(:,imode))
    END IF

    ! 7) count occurrences of MDT>MDTMAX (add to NBADMDT) and set 100%
    !    where occurs
    WHERE (mdt(:,imode) > mdtmax)
      nbadmdt(:,imode) = 1
      mdtfixflag(:,imode)=100.0
      ! mdtfixflag enables to track proportion of timesteps that fix is applied
      mdtfixsink(:,imode)=nd(:,imode)*mdt(:,imode)*fac(:)/mm_da/dtc ! moles/s
      ! mdtfixsink stores total mass removed when fix is applied
    END WHERE

    ! 8) print out MDT too low or too high occurrences if > 0
    IF (SUM(nbadmdt(:,imode)) > 0 .AND. verbose > 0) THEN
      WRITE(umMessage,'(A53)') 'MDT<MDTMIN or MDT>MDTMAX, ND=0: '//            &
        'IMODE,MDTMAX,NBADMDT:'
      WRITE(umMessage,'(I6,1E12.3,I12)') imode,mdtmax,SUM(nbadmdt(:,imode))
    END IF

    ! 9) set ND->0 where MDT too low (& set MDT->MMID)
    WHERE (mdt(:,imode) < mdtmin)
      nd(:,imode) = 0.0
      mdt(:,imode) = mmid(imode)
    END WHERE

    ! 10) set ND->0 where MDT too high (& set MDT->MMID)
    WHERE (mdt(:,imode) > mdtmax)
      nd(:,imode) = 0.0
      mdt(:,imode) = mmid(imode)
    END WHERE
  END IF ! if mode is active in this aerosol configuration
END DO ! loop over modes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (verbose > 1) THEN
  !  Check min, max values of ND, MDT and MD after 1st check
  WRITE(umMessage,'(A50)') 'BEFORE CALL TO UKCA_AERO_STEP (after fix 1)'
  DO imode=1,nmodes
    IF (mode(imode)) THEN
      WRITE(umMessage,'(A9,I5,2E12.4)') 'ND : ',imode,MINVAL(nd(:,imode)),     &
        MAXVAL(nd(:,imode))
      WRITE(umMessage,'(A9,I5,2E12.4)') 'MDT: ',imode,MINVAL(mdt(:,imode)),    &
        MAXVAL(mdt(:,imode))
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WRITE(umMessage,'(A9,I5,2E12.4)') 'MD : ',imode,                     &
            MINVAL(md(:,imode,icp)),MAXVAL(md(:,imode,icp))
        END IF
      END DO
    END IF
  END DO
END IF ! if VERBOSE > 1


IF (verbose > 1) THEN
  !  Check min, max values of ND, MDT and MD after 2nd check
  WRITE(umMessage,'(A50)') 'BEFORE CALL TO UKCA_AERO_STEP (after fix 2)'
  DO imode=1,nmodes
    IF (mode(imode)) THEN
      WRITE(umMessage,'(A9,I4,2E12.4)') 'ND : ',imode,MINVAL(nd(:,imode)),     &
        MAXVAL(nd(:,imode))
      WRITE(umMessage,'(A9,I4,2E12.4)') 'MDT: ',imode,MINVAL(mdt(:,imode)),    &
        MAXVAL(mdt(:,imode))
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WRITE(umMessage,'(A9,I4,2E12.4)') 'MD : ',imode,                     &
            MINVAL(md(:,imode,icp)),MAXVAL(md(:,imode,icp))
        END IF
      END DO
    END IF
  END DO
END IF ! if VERBOSE > 1

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_check_artefacts

! #############################################################################

SUBROUTINE ukca_mode_check_mdt(nbox, imode, mdt, md, nd, mask1)

! Check to see if mdt is outside limits and reset nd and md to defaults if
! it is.
! Called from UKCA_COAGWITHNUCL.

USE ukca_config_specification_mod, ONLY:                                       &
    glomap_variables

USE ukca_mode_setup,    ONLY:                                                  &
    nmodes,                                                                    &
    mode_nuc_sol,                                                              &
    mode_ait_sol,                                                              &
    mode_acc_sol,                                                              &
    mode_cor_sol,                                                              &
    mode_ait_insol,                                                            &
    mode_acc_insol,                                                            &
    mode_cor_insol,                                                            &
    mode_sup_insol

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ukca_types_mod, ONLY: logical_32

IMPLICIT NONE

INTEGER, INTENT(IN) :: nbox                 ! No. of elements
INTEGER, INTENT(IN) :: imode                ! mode index

REAL, INTENT(IN OUT) :: nd(nbox,nmodes)
! Each mode's particle number density (ptcls cm^-3)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables%ncp)
! Each mode's average component mass per particle (molecules ptcle^-1)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
! Each mode's average total mass per particle (over all cpts) molecules ptcl^-1)

LOGICAL (KIND=logical_32), INTENT(IN OUT) :: mask1(nbox)
! Logical mask array -- gets set false where mdt < mdtmin or > mdtmax

! Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: mfrac_0 (:,:)
REAL,    POINTER :: mhi (:)
REAL,    POINTER :: mlo (:)
REAL,    POINTER :: mmid (:)
INTEGER, POINTER :: ncp

INTEGER :: icp                              ! Loop counter for components

REAL    :: mdtmin                  ! min MDT value for each mode
REAL    :: mdtmax                  ! max MDT value for each mode

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_CHECK_MDT'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mfrac_0     => glomap_variables%mfrac_0
mhi         => glomap_variables%mhi
mlo         => glomap_variables%mlo
mmid        => glomap_variables%mmid
ncp         => glomap_variables%ncp

! Below carries out key element of part A of MDT checks within GLOMAP code.
! Protects against occasional occurrences found in ukca_coagwithinucl.
! Essentially protects against strange values after coagulation.

!! Original checks retained for reference
!!        MDTMIN=MLO(IMODE)*0.001  ! set equiv. to DPLIM0*0.1
!!        MDTMAX=MHI(IMODE)*1000.0 ! set equiv. to DPLIM1*10.0

! 1) Min MDT criteria
! -------------------
! 1.1   Soluble modes (all) : min MDT factor-10 lower in diam than min of mode
IF ((imode == mode_nuc_sol) .OR. (imode == mode_ait_sol)  .OR.                 &
    (imode == mode_acc_sol) .OR. (imode == mode_cor_sol)) THEN
  mdtmin = mlo(imode)*0.001
END IF

! 1.2 Insoluble Aitken mode : min MDT factor-10 lower in diam than min of mode
IF (imode == mode_ait_insol) mdtmin = mlo(imode)*0.001

! 1.3 Insoluble accumn mode : only minimal constraint: OK down to mlo(1)*0.001
! only minimal constraint on mode 6
! (where contains mainly MSP can be v small)
IF (imode == mode_acc_insol) mdtmin = mlo(1)*0.001

! 1.4 Insoluble coarse mode : min MDT factor-10 lower in diam than min of mode
IF (imode == mode_cor_insol) mdtmin = mlo(imode)*0.001

! 1.5 Insoluble super-coarse mode : min MDT factor-10 lower in diam than
!     min of mode
IF (imode == mode_sup_insol) mdtmin = mlo(imode)*0.001

! 2) Max MDT criteria
! -------------------
! All modes have same maxMDT: max MDT factor-10 higher in diam than max of mode
mdtmax = mhi(imode)*1000.0

! below makes sure MDTMIN is never below MDTMIN_PAR (e.g. for mode 1)
IF (mdtmin < mdtmin_par) mdtmin = mdtmin_par

! where MDT<MDTMIN, set ND to zero (MDT,MD reset at end routine)
WHERE (mask1(:) .AND. (mdt(:,imode) < mdtmin))
  ! seem to very occasionally get very low MDT occurring which can have
  ! knock-on effects....  this repeats check applied in ukca_aero_ctl
  ! to make sure values are "sensible".

  nd(:,imode) = 0.0
  mask1(:) = .FALSE. ! set false so not used for other icp values
                     ! n.b. MD and MDT will be re-set at end of routine
END WHERE

! where MDT>MDTMAX, set ND to zero (MDT,MD reset at end routine)
! check also for when can get very high MDT occurring which can have
! knock-on effects....  this repeats check applied in ukca_aero_ctl
! to make sure values are "sensible".

WHERE (mask1(:) .AND. (mdt(:,imode) > mdtmax))
  nd(:,imode) = 0.0
  mask1(:) = .FALSE. ! set false so not used for other icp values
                     ! n.b. MD and MDT will be re-set at end of routine
END WHERE

DO icp=1,ncp
  IF (component(imode,icp)) THEN
    ! where MDT too low after advection set ND to zero and set default MD
    WHERE (mdt(:,imode) < mdtmin)
      md(:,imode,icp) = mmid(imode)*mfrac_0(imode,icp)
    END WHERE
    ! where MDT too high after advection set ND to zero and set default MD
    WHERE (mdt(:,imode) > mdtmax)
      md(:,imode,icp) = mmid(imode)*mfrac_0(imode,icp)
    END WHERE
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_check_mdt

END MODULE ukca_mode_check_artefacts_mod
