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
!   Carries out ageing of particles in insoluble mode.
!   Calculates the number of particles which will be coated
!   by the total soluble material and transfers number and
!   mass accordingly.
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
MODULE ukca_ageing_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_AGEING_MOD'

CONTAINS

SUBROUTINE ukca_ageing(nbox,nchemg,nbudaer,nd,md,mdt,                          &
 ageterm1,ageterm2,wetdp,bud_aer_mas,iagecoagnucl67)
!----------------------------------------------------------
!
! Purpose
! ------
! Carries out ageing of particles in insoluble mode.
! Calculates the number of particles which will be coated
! by the total soluble material and transfers number and
! mass accordingly.
! Number of ptcls given by assuming 1 particle ages with
! a 1 molecule layer thickness of particle of geometric
! mean size of mode.  i.e. 1 ptcl requires r_g^2/r_molec^2
! = 10^4,10^6,10^8 molecules for Aitken, accum, coarse modes
!
! Inputs
! ------
! NBOX      : Number of grid boxes
! nchemg    : Number of gas phase chemistry tracers
! nbudaer   : Number of aerosol budget fields
! ND        : Aerosol ptcl no. concentration (ptcls per cc)
! MD        : Component median aerosol mass (molecules per ptcl)
! MDT       : Total median aerosol mass (molecules per ptcl)
! AGETERM1  : Depletion rate of each component (molecules cpt/cc/DTZ)
!             from condensation onto the 3 insoluble modes.
! AGETERM2  : Rate of accomodation of material to each insoluble mode
!             as a result of coagulation with smaller soluble modes
!             (in molecules cpt /cm3/DTZ)
! WETDP     : Wet diameter corresponding to mean sized particle
! IAGECOAGNUCL67: Switch ageing,coag & nucl involving modes6&7 on(1)/off(0)
!
! Outputs
! -------
! ND        : Updated number concentration [each mode] (ptcls/cc)
! MD        : Updated avg cpt   mass conc. [each mode] (molecules/ptcl)
! MDT       : Updated avg total mass conc. [each mode] (molecules/ptcl)
! BUD_AER_MAS : Aerosol mass budgets
!
! Local variables
! ---------------
! AGE1PTCL  : Mass of soluble material needed to age 1 ptcl (molecules)
! NAGED     : # of insoluble ptcls aged to soluble mode [total] (/cc)
! NAGED_JV  : # of insoluble ptcls aged to soluble mode [by gas](/cc)
! TOTAGE_JV : Ageing flux to cpt by particular gas (molecules gas/cc)
! TOTAGE1   : Total ageing flux to cpt by conden.  (molecules gas/cc)
! TOTAGE2   : Total ageing flux to cpt by coaguln. (molecules gas/cc)
! TOTAGE    : Total ageing flux to cpt (coag+cond) (molecules gas/cc)
! NDINSNEW  : # in ins mode after reduction due to ageing (/cc)
! NDSOLNEW  : # in corresp. sol mode after reduction due to ageing (/cc)
! F_MM      : Ratio of molar masses of condensable gas to aerosol cpt
! CP_COAG_ADDED : Switch for whether added on ageing flux by
!                 coagulation to that cpt already
!                 (loop over jv --- need to make sure only count once)
! TOPMODE   : Highest number mode for which ageing occurs.
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! CONC_EPS  : Likewise, threshold for soluble material (molecules per cc)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! NCP       : Number of possible aerosol components
! MODE      : Which modes are being carried
! COMPONENT : Which components are in each of modes
! MM        : Molar masses of components (kg per mole)
! NUM_EPS   : Value of NEWN below which do not recalculate MD (per cc)
!             or carry out process
! CP_SU     : Component where sulfate is stored
! CP_BC     : Component in which black carbon is stored
! CP_OC     : Component in which primary organic carbon is stored
! CP_DU     : Component where dust is stored
! CP_SO     : Component where secondary organic carbon is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! MM_GAS    : Array of molar masses for gas phase species (kg/mol)
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY:                                       &
    glomap_variables

USE ukca_mode_setup,    ONLY:                                                  &
    nmodes,                                                                    &
    cp_su,                                                                     &
    cp_bc,                                                                     &
    cp_oc,                                                                     &
    cp_du,                                                                     &
    cp_so,                                                                     &
    mode_nuc_sol,                                                              &
    mode_cor_sol,                                                              &
    mode_ait_insol,                                                            &
    mode_acc_insol,                                                            &
    mode_cor_insol

USE ukca_setup_indices, ONLY: condensable,                                     &
        condensable_choice, mm_gas, dimen, nmasagedsuintr52,                   &
        nmasagedocintr52, nmasagedsointr52, nmasagedbcintr52,                  &
        nmasagedsuintr63, nmasagedocintr63, nmasagedsointr63,                  &
        nmasagedduintr63, nmasagedsuintr74, nmasagedocintr74,                  &
        nmasagedsointr74, nmasagedduintr74

USE yomhook,            ONLY: lhook, dr_hook
USE parkind1,           ONLY: jprb, jpim

IMPLICIT NONE
!
! .. Subroutine Interface
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: nchemg
INTEGER, INTENT(IN) :: nbudaer
REAL, INTENT(IN)    :: ageterm1(nbox,3,nchemg)
REAL, INTENT(IN)    :: ageterm2(nbox,4,3,glomap_variables%ncp)
REAL, INTENT(IN)    :: wetdp(nbox,nmodes)
REAL, INTENT(IN OUT) :: nd(nbox,nmodes)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables%ncp)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
REAL, INTENT(IN OUT) :: bud_aer_mas(nbox,0:nbudaer)
INTEGER, INTENT(IN) :: iagecoagnucl67
!
! .. Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: mm(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: num_eps(:)

INTEGER :: jl
INTEGER :: jv
INTEGER :: imode
INTEGER :: jmode
INTEGER :: icp
INTEGER :: cp_coag_added(glomap_variables%ncp)
INTEGER :: topmode
REAL    :: totage(glomap_variables%ncp)
REAL    :: totage_jv
REAL    :: totage1(glomap_variables%ncp)
REAL    :: totage2(glomap_variables%ncp)
REAL    :: age1ptcl
REAL    :: naged
REAL    :: naged_jv(nchemg)
REAL    :: ndinsnew
REAL    :: ndsolnew
REAL    :: f_mm(glomap_variables%ncp)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_AGEING'

!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mm          => glomap_variables%mm
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
num_eps     => glomap_variables%num_eps

!set limit of modes to be aged
IF (iagecoagnucl67 == 1 ) THEN
  topmode = mode_cor_insol
ELSE
  topmode = mode_ait_insol
END IF

DO imode=mode_ait_insol,topmode ! loop over insoluble modes
  IF (mode(imode)) THEN
    DO jl=1,nbox
      IF (nd(jl,imode) > num_eps(imode)) THEN
        cp_coag_added(:)=0
        totage(:)=0.0
        totage1(:)=0.0 ! from condensation
        totage2(:)=0.0 ! from coagulation
        f_mm(:)=0.0 ! need to initialise F_MM to 0 so it is defined
                    ! for components which don't condense
        naged=0.0
        DO jv=1,nchemg
          totage_jv=0.0
          naged_jv(jv)=0.0
          IF (condensable(jv)) THEN

            ! .. Below add on amount of soluble material taken up by insoluble
            !     modes
            ! .. as result of condensation of this gas onto insoluble modes
            ! .. This is stored in AGETERM1(:,3,NCP) for the 3 insoluble modes
            icp=condensable_choice(jv)
            f_mm(icp)=mm_gas(jv)/mm(icp)
            totage_jv=ageterm1(jl,imode-4,jv)/f_mm(icp)
            totage (icp)=totage (icp)+totage_jv
            totage1(icp)=totage1(icp)+totage_jv
            ! .. condensable material taken up by aerosol assumed 100% soluble.
            ! .. AGETERM1 is from condensation onto insoluble modes
            ! .. Divide by F_MM because AGETERM1 is in "molecules" of cpt,
            ! .. whereas TOTAGE needs to be in molecules of condensible gas
            ! .. (H2SO4/Sec_Org) since AGE1PTCL is in these units.

            ! .. Below add on amount of soluble material taken up by insoluble
            !     modes
            ! .. as result of coagulation with soluble modes
            !
            ! .. IF(CP_COAG_ADDED(ICP) == 0) checks if already added on this
            ! .. aerosol component's coagulated material (make sure don't
            !     doublecount)
            !
            IF (cp_coag_added(icp) == 0) THEN
              DO jmode=mode_nuc_sol,mode_cor_sol
                totage_jv=totage_jv+                                           &
                   ageterm2(jl,jmode,imode-4,icp)/f_mm(icp)
                ! .. add on amount coagulated
                totage (icp)=totage (icp)+                                     &
                   ageterm2(jl,jmode,imode-4,icp)/f_mm(icp)
                totage2(icp)=totage2(icp)+                                     &
                   ageterm2(jl,jmode,imode-4,icp)/f_mm(icp)
              END DO
              cp_coag_added(icp)=1
            END IF
            ! .. condensable material taken up by aerosol assumed 100% soluble.
            ! .. AGETERM2 is from coagn of soluble modes with larger insoluble
            ! .. modes Divide by F_MM because AGETERM2 is in "molecules" of cpt,
            ! .. whereas TOTAGE needs to be in molecules of condensible gas
            ! .. (H2SO4/Sec_Org) since AGE1PTCL is in these units.
            !
            age1ptcl=wetdp(jl,imode)*wetdp(jl,imode)/dimen(jv)/dimen(jv)
            ! above is number of molecules of condensable gas to age 1 particle
            naged_jv(jv)=totage_jv/age1ptcl/10.0
            !!! divide by 10 is to make particle age by 10 monolayers
            !    instead of 1
            !!          NAGED_JV(JV)=TOTAGE_JV/AGE1PTCL
            ! here keep as single monolayer ageing.
            naged=naged+naged_jv(jv)

          END IF ! if gas phase species is condensable
        END DO ! loop over gas phase species
        !
        IF (naged > num_eps(imode)) THEN ! if significant ptcls to age
          !
          IF (naged > nd(jl,imode)) THEN
            naged=nd(jl,imode) ! limit so no -ves
            totage(:)=totage(:)*nd(jl,imode)/naged
            ! above reduces ageing if limited by insoluble particles
          END IF
          !
          ndinsnew=nd(jl,imode)-naged
          ! set new insoluble mode no. (but don't update ND yet)
          ndsolnew=nd(jl,imode-3)+naged
          ! set new   soluble mode no. (but don't update ND yet)
          !

          IF (SUM(totage) > 0.0) THEN
            DO icp=1,ncp
              ! below transfers aged cpt masses from ins modes
              !  (doesn't include SU)
              IF (component(imode-3,icp)) THEN
                ! above if statement checks whether cpt is in corresponding
                !  soluble mode
                IF (imode == mode_ait_insol) THEN
                  IF ((icp == cp_su) .AND. (nmasagedsuintr52 > 0))             &
                bud_aer_mas(jl,nmasagedsuintr52)=                              &
                bud_aer_mas(jl,nmasagedsuintr52)+totage(icp)*f_mm(icp)
                  IF ((icp == cp_oc) .AND. (nmasagedocintr52 > 0))             &
                bud_aer_mas(jl,nmasagedocintr52)=                              &
                bud_aer_mas(jl,nmasagedocintr52)+totage(icp)*f_mm(icp)
                  IF ((icp == cp_so) .AND. (nmasagedsointr52 > 0))             &
                bud_aer_mas(jl,nmasagedsointr52)=                              &
                bud_aer_mas(jl,nmasagedsointr52)+totage(icp)*f_mm(icp)
                  !! above accounts for transfer of mass of condensed/coagulated
                  !! material multiply above by F_MM 'cos TOTAGE is in molecules
                  !! of gas phase species (H2SO4/Sec_Org) whereas needs to be
                  !! in "molecules" of aerosol component (CP_SU/CP_OC/CP_SU)
                  IF ((icp == cp_bc) .AND. (nmasagedbcintr52 > 0))             &
                bud_aer_mas(jl,nmasagedbcintr52)=                              &
                bud_aer_mas(jl,nmasagedbcintr52)+naged*md(jl,imode,icp)
                  IF ((icp == cp_oc) .AND. (nmasagedocintr52 > 0))             &
                bud_aer_mas(jl,nmasagedocintr52)=                              &
                bud_aer_mas(jl,nmasagedocintr52)+naged*md(jl,imode,icp)
                  !! .. above 2 lines account for transfer of mass due to aged
                  !      aerosol
                END IF ! if mode is Aitken-insoluble
                IF (imode == mode_acc_insol) THEN
                  IF ((icp == cp_su) .AND. (nmasagedsuintr63 > 0))             &
                bud_aer_mas(jl,nmasagedsuintr63)=                              &
                bud_aer_mas(jl,nmasagedsuintr63)+totage(icp)*f_mm(icp)
                  IF ((icp == cp_oc) .AND. (nmasagedocintr63 > 0))             &
                bud_aer_mas(jl,nmasagedocintr63)=                              &
                bud_aer_mas(jl,nmasagedocintr63)+totage(icp)*f_mm(icp)
                  IF ((icp == cp_so) .AND. (nmasagedsointr63 > 0))             &
                bud_aer_mas(jl,nmasagedsointr63)=                              &
                bud_aer_mas(jl,nmasagedsointr63)+totage(icp)*f_mm(icp)
                  !! above accounts for transfer of mass of condensed/coagulated
                  !! material multiply above by F_MM 'cos TOTAGE is in molecules
                  !! of gas phase species (H2SO4/Sec_Org) whereas needs to be in
                  !! "molecules" of aerosol component (CP_SU/CP_OC/CP_SU)
                  IF ((icp == cp_du) .AND. (nmasagedduintr63 > 0))             &
                bud_aer_mas(jl,nmasagedduintr63)=                              &
                bud_aer_mas(jl,nmasagedduintr63)+naged*md(jl,imode,icp)
                  !! .. above 1 line accounts for transfer of mass due to aged
                  !      aerosol
                END IF ! if mode is accum.-insoluble
                IF (imode == mode_cor_insol) THEN
                  IF ((icp == cp_su) .AND. (nmasagedsuintr74 > 0))             &
                bud_aer_mas(jl,nmasagedsuintr74)=                              &
                bud_aer_mas(jl,nmasagedsuintr74)+totage(icp)*f_mm(icp)
                  IF ((icp == cp_oc) .AND. (nmasagedocintr74 > 0))             &
                bud_aer_mas(jl,nmasagedocintr74)=                              &
                bud_aer_mas(jl,nmasagedocintr74)+totage(icp)*f_mm(icp)
                  IF ((icp == cp_so) .AND. (nmasagedsointr74 > 0))             &
                bud_aer_mas(jl,nmasagedsointr74)=                              &
                bud_aer_mas(jl,nmasagedsointr74)+totage(icp)*f_mm(icp)
                  !! above accounts for transfer of mass of condensed/coagulated
                  !! material multiply above by F_MM 'cos TOTAGE is in molecules
                  !! of gas phase species (H2SO4/Sec_Org) whereas needs to be in
                  !! "molecules" of aerosol component (CP_SU/CP_OC/CP_SU)
                  IF ((icp == cp_du) .AND. (nmasagedduintr74 > 0))             &
                bud_aer_mas(jl,nmasagedduintr74)=                              &
                bud_aer_mas(jl,nmasagedduintr74)+naged*md(jl,imode,icp)
                  !! .. above 1 line accounts for transfer of mass due to aged
                  !      aerosol
                END IF ! if mode is coarse-insoluble

                !! .. below calculates new cpt total masses in soluble modes
                !! .. due to transfer of mass from Ait-ins/acc-ins/cor-ins
                !! .. (n.b. insoluble avg. masses unchanged)
                !! .. (n.b. soluble mode mass not updated due to condensation
                !! ..  onto insoluble mode yet)

                IF (component(imode,icp)) THEN ! if in sol mode & ins mode
                  ! if sol. mode cpt is in insoluble mode then update
                  ! corresponding soluble mode cpt mass due to trans from
                  ! ins. mode & cond onto ins (n.b. coag already transferred
                  ! in coagulation routine using COAG_MODE)
                  md(jl,imode-3,icp)=(nd(jl,imode-3)*md(jl,imode-3,icp)        &
                   +naged*md(jl,imode,icp)+totage1(icp)*f_mm(icp))/ndsolnew
                ELSE ! if in sol mode but not in insoluble mode
                  ! if sol. mode component is not in insoluble mode then just
                  ! update soluble mode cpt mass due to cond onto ins
                  ! (& change in number)
                  ! (n.b. coag already transferred in coagulation routine
                  !  using COAG_MODE)
                  md(jl,imode-3,icp)=(nd(jl,imode-3)*md(jl,imode-3,icp)        &
                                          +totage1(icp)*f_mm(icp))/ndsolnew
                END IF

              END IF ! if component is in corresponding soluble mode
            END DO ! loop over components
          END IF ! if total amount of accomodated material > 0
          !
          !! update ND and MDT for ins mode
          nd(jl,imode  )=ndinsnew
          mdt(jl,imode)=0.0
          DO icp=1,ncp
            IF (component(imode,icp)) THEN
              mdt(jl,imode)=mdt(jl,imode)+md(jl,imode,icp)
            END IF
          END DO ! end loop over cpts
          !
          !! update ND and MDT for sol mode
          nd(jl,imode-3)=ndsolnew
          mdt(jl,imode-3)=0.0
          DO icp=1,ncp
            IF (component(imode-3,icp)) THEN
              mdt(jl,imode-3)=mdt(jl,imode-3)+md(jl,imode-3,icp)
            END IF
          END DO ! end loop over cpts
          !
        END IF ! if number of aged particles > NUM_EPS(IMODE)
      END IF ! if some particles in insoluble modes (ND>epsilon)
    END DO ! end loop over boxes
  END IF ! if insoluble mode is present (IMODE=5,7)
END DO ! Loop IMODE=5,topmode
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_ageing
END MODULE ukca_ageing_mod
