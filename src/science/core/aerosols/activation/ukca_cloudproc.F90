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
!    Cloud processing of aerosol.
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
MODULE ukca_cloudproc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CLOUDPROC_MOD'

CONTAINS

SUBROUTINE ukca_cloudproc(nbox,nbudaer,nd,md,mdt,drydp,                        &
                          lowcloud,vfac,act,iactmethod,bud_aer_mas)
!----------------------------------------------------------------------
!
!  Purpose
!  -------
!  If gridbox is in-cloud then accounts for the fact that
!  particles in the Aitken soluble mode will be activated
!  to form cloud droplets.
!  Calculates the fraction of the number & mass of Aitsol
!  particles which have dry radii larger than the activation
!  radius "ACT". Then transfers this mass and number
!  from Aitsol to accsol mode (updates number and mass of each)
!  so that Aitken soluble mode particles which were large
!  enough will then receive in-cloud produced S(VI) as
!  they will now be in the accsol mode.
!
!  Currently uses constant value for activation radius as in
!  first versions of GLOMAP (as published in Spracklen et al (2005)
!  This is IACTMETHOD=1.
!
!  Later want to use input updraft velocity to diagnose
!  maximum supersaturation and hence activation radius using
!  Nenes & Seinfeld (2003) parameterisation (IACTMETHOD=2).
!
!  Inputs
!  ------
!  NBOX        : Number of grid boxes
!  nbudaer     : Number of aerosol budget fields
!  ND          : Aerosol ptcl number density for mode (cm^-3)
!  MDT         : Avg tot mass of aerosol ptcl in mode (particle^-1)
!  MD          : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!  DRYDP       : Geometric mean dry diameter of particles in each mode (m)
!  LOWCLOUD    : Horizontal low cloud fraction
!  VFAC        : Vertical low cloud fraction
!  ACT         : Particle dry radius above which activation is assumed
!  IACTMETHOD  : Switch for activation method (0=off,1=fixed ract,2=NSO3 scheme)
!
!  Outputs
!  -------
!  ND          : Aerosol ptcl number density for mode (cm^-3)
!  MDT         : Avg tot mass of aerosol ptcl in mode (particle^-1)
!  MD          : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!  BUD_AER_MAS : Aerosol mass budgets
!
!  Local variables
!  ---------------
!  DP       : Number median diameter of mode (m)
!  LNRATN   : Log of ratio of threshold diameter to no. median diameter
!  ERFNUM   : LNRATN/(sqrt(2)*log(sigmag))
!  FRAC_N   : Fraction of ptcl number which is within bounds
!  DELN     : Number concentration to transfer due to cloud processing
!  LOG2SG   : log(sigmag)*log(sigmag)
!  DP2      : Volume median diameter of mode (m)
!  LNRATM   : Log of ratio of threshold diameter to mass median diameter
!  ERFMAS   : LNRATM/(sqrt(2)*log(sigmag))
!  FRAC_M   : Fraction of ptcl mass which is within bounds
!  DM       : Cpt mass concentration to transfer due to cloud processing
!  NEWN     : New particle number conc in mode MODE_AIT_SOL
!  NEWNP1   : New particle number conc in mode MODE_ACC_SOL
!
!  Inputted by module UKCA_MODE_SETUP
!  ----------------------------------
!  NMODES      : Number of possible aerosol modes
!  NCP         : Number of possible aerosol components
!  COMPONENT   : Logical variable defining which cpt are in which dsts
!  SIGMAG      : Geometric standard deviation for each mode
!  MFRAC_0     : Initial mass fraction to set when no particles.
!  MMID      : Mass of particle with dp=dpmed_g=exp(0.5*(lndp0+lndp1)) (ptcl^-1)
!  NUM_EPS     : Value of NEWN below which do not recalculate MD (per cc)
!              : or carry out process
!  CP_SU       : Component where sulfate is stored
!  CP_BC       : Component where black carbon is stored
!  CP_OC       : Component where organic carbon is stored
!  CP_SO       : Component where condensible organic species is stored
!  MODE_AIT_SOL : Parameter used to idenitfy Soluble Aitkin mode
!  MODE_ACC_SOL : Parameter used to idenitfy Soluble Accumulation mode
!
!  Inputted by module UKCA_SETUP_INDICES
!  -------------------------------------
!  Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: glomap_variables

USE ukca_mode_setup,    ONLY: nmodes,                                          &
                              cp_su, cp_bc, cp_oc, cp_so, cp_mp,               &
                              cp_no3, cp_nh4, mode_ait_sol, mode_acc_sol

USE ukca_setup_indices, ONLY: nmasprocsuintr23,                                &
                              nmasprocbcintr23, nmasprococintr23,              &
                              nmasprocsointr23,                                &
                              nmasprocntintr23, nmasprocnhintr23,              &
                              nmasprocmpintr23
USE yomhook,            ONLY: lhook, dr_hook
USE parkind1,           ONLY: jprb, jpim
USE ukca_um_legacy_mod, ONLY: umErf

IMPLICIT NONE

! .. Input/output variables
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: nbudaer
INTEGER, INTENT(IN) :: iactmethod
REAL, INTENT(IN)    :: drydp(nbox,nmodes)
REAL, INTENT(IN)    :: lowcloud(nbox),vfac(nbox)
REAL, INTENT(IN)    :: act
REAL, INTENT(IN OUT) :: nd(nbox,nmodes)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables%ncp)
REAL, INTENT(IN OUT) :: bud_aer_mas(nbox,0:nbudaer)

! .. Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: mfrac_0(:,:)
REAL,    POINTER :: mmid(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: num_eps(:)
REAL,    POINTER :: sigmag(:)

INTEGER :: jl
INTEGER :: icp
REAL    :: f
REAL    :: dp
REAL    :: lnratn
REAL    :: erfnum
REAL    :: frac_n
REAL    :: deln
REAL    :: log2sg
REAL    :: dp2
REAL    :: lnratm
REAL    :: erfmas
REAL    :: frac_m
REAL    :: dm(glomap_variables%ncp)
REAL    :: newn
REAL    :: newnp1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CLOUDPROC'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mfrac_0     => glomap_variables%mfrac_0
mmid        => glomap_variables%mmid
ncp         => glomap_variables%ncp
num_eps     => glomap_variables%num_eps
sigmag      => glomap_variables%sigmag

IF (iactmethod == 1) THEN
  DO jl=1,nbox
    f=lowcloud(jl)*vfac(jl)
    IF (f > 0.0) THEN ! if in cloud
      ! apply for Aitsol mode (transfer to accsol)
      ! .. calculate fraction of Aitsol mode no & mass with r>ACT
      dp=drydp(jl,mode_ait_sol)
      IF (nd(jl,mode_ait_sol) > num_eps(mode_ait_sol)) THEN
        lnratn=LOG(act*2.0/dp) ! use 2nd threshold to work out fraction
        erfnum=lnratn/SQRT(2.0)/LOG(sigmag(mode_ait_sol))
        frac_n=0.5*(1.0+umErf(erfnum)) ! fraction remaining in mode
        IF (frac_n < 0.5) frac_n=0.5 ! limit DELN to be max = half # of ptcls
        IF (frac_n < 0.99) THEN ! if more than 1% activated
          deln=nd(jl,mode_ait_sol)*(1.0-frac_n)
          deln=deln*f ! modify change in number taking into account gridbox
          !                        cloud fraction
          log2sg=LOG(sigmag(mode_ait_sol))*LOG(sigmag(mode_ait_sol))
          dp2=EXP(LOG(dp)+3.0*log2sg) ! volume median diameter
          lnratm=LOG(act*2.0/dp2)
          erfmas=lnratm/SQRT(2.0)/LOG(sigmag(mode_ait_sol))
          frac_m=0.5*(1.0+umErf(erfmas))
          newn=nd(jl,mode_ait_sol)-deln
          newnp1=nd(jl,mode_acc_sol)+deln
          IF (newn > num_eps(mode_ait_sol)) THEN
            DO icp=1,ncp
              IF (component(mode_ait_sol,icp)) THEN
                ! .. calculate amount of each component to transfer to
                !    next mode up (use old no/mass)
                dm(icp)=md(jl,mode_ait_sol,icp)*nd(jl,mode_ait_sol)*(1.0-frac_m)
                dm(icp)=dm(icp)*f ! modify change in mass taking into account
                !                                 gridbox cloud fraction
              ELSE
                dm(icp)=0.0
              END IF
              IF ((icp == cp_su) .AND. (nmasprocsuintr23 > 0))                 &
                           bud_aer_mas(jl,nmasprocsuintr23)=                   &
                           bud_aer_mas(jl,nmasprocsuintr23)+dm(icp)
              IF ((icp == cp_bc) .AND. (nmasprocbcintr23 > 0))                 &
                           bud_aer_mas(jl,nmasprocbcintr23)=                   &
                           bud_aer_mas(jl,nmasprocbcintr23)+dm(icp)
              IF ((icp == cp_oc) .AND. (nmasprococintr23 > 0))                 &
                           bud_aer_mas(jl,nmasprococintr23)=                   &
                           bud_aer_mas(jl,nmasprococintr23)+dm(icp)
              IF ((icp == cp_so) .AND. (nmasprocsointr23 > 0))                 &
                           bud_aer_mas(jl,nmasprocsointr23)=                   &
                           bud_aer_mas(jl,nmasprocsointr23)+dm(icp)
              IF ((icp == cp_no3) .AND. (nmasprocntintr23 > 0))                &
                           bud_aer_mas(jl,nmasprocntintr23)=                   &
                           bud_aer_mas(jl,nmasprocntintr23)+dm(icp)
              IF ((icp == cp_nh4) .AND. (nmasprocnhintr23 > 0))                &
                           bud_aer_mas(jl,nmasprocnhintr23)=                   &
                           bud_aer_mas(jl,nmasprocnhintr23)+dm(icp)
              IF ((icp == cp_mp) .AND. (nmasprocmpintr23 > 0))                 &
                           bud_aer_mas(jl,nmasprocmpintr23)=                   &
                           bud_aer_mas(jl,nmasprocmpintr23)+dm(icp)
            END DO
            !
            ! .. first remove mass to be transferred from mode MODE_AIT_SOL
            mdt(jl,mode_ait_sol)=0.0
            DO icp=1,ncp
              IF (component(mode_ait_sol,icp)) THEN
                md(jl,mode_ait_sol,icp) = ( nd(jl,mode_ait_sol) *              &
                                            md(jl,mode_ait_sol,icp) -          &
                                            dm(icp) ) / newn
                mdt(jl,mode_ait_sol) = mdt(jl,mode_ait_sol) +                  &
                                        md(jl,mode_ait_sol,icp)
              ELSE
                md(jl,mode_ait_sol,icp)=0.0
              END IF ! COMPONENT(MODE_AIT_SOL,ICP)
            END DO
            ! .. now set new number to mode MODE_AIT_SOL
            nd(jl,mode_ait_sol)=newn ! set particle number to new value
            !
            ! .. then add mass to be transferred to mode MODE_ACC_SOL
            mdt(jl,mode_acc_sol)=0.0
            DO icp=1,ncp
              IF (component(mode_acc_sol,icp)) THEN
                md(jl,mode_acc_sol,icp) = ( nd(jl,mode_acc_sol) *              &
                                            md(jl,mode_acc_sol,icp) +          &
                                            dm(icp) ) / newnp1
                mdt(jl,mode_acc_sol) = mdt(jl,mode_acc_sol) +                  &
                                        md(jl,mode_acc_sol,icp)
              ELSE
                md(jl,mode_acc_sol,icp)=0.0
              END IF ! COMPONENT(MODE_ACC_SOL,ICP)
            END DO
            ! .. now set new number to mode MODE_ACC_SOL
            nd(jl,mode_acc_sol)=newnp1
          END IF ! IF NEWN>0
          !
        END IF ! if FRAC_N<0.99 (if more than 1% activated)
      ELSE
        DO icp=1,ncp
          IF (component(mode_ait_sol,icp)) THEN
            md(jl,mode_ait_sol,icp) = mmid(mode_ait_sol) *                     &
                                      mfrac_0(mode_ait_sol,icp)
          END IF
        END DO
        mdt(jl,mode_ait_sol)=mmid(mode_ait_sol)
      END IF ! if significant number of particles in lower mode
    END IF ! if F>0 (if in cloud)
  END DO ! end loop over boxes
END IF ! if iactmethod=1

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_cloudproc
END MODULE ukca_cloudproc_mod
