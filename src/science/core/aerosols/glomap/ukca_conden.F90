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
!    Calculates condensation of condensable cpt vapours onto pre-existing
!    aerosol particles. Includes switch for using either Fuchs (1964) or
!    modified Fuchs and Sutugin (1971) calculation of CC.
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
MODULE ukca_conden_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CONDEN_MOD'

CONTAINS

SUBROUTINE ukca_conden(nbox,nchemg,nbudaer,ifuchs,idcmfp,icondiam,             &
                       nd,tsqrt,rhoa,airdm3,dtz,wetdp,pmid,t,                  &
                       md,mdt,gc,bud_aer_mas,                                  &
                       delgc_cond,ageterm1,s_cond_s)
!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculates condensation of condensable cpt vapours onto pre-existing
! aerosol particles. Includes switch for using either Fuchs (1964) or
! modified Fuchs and Sutugin (1971) calculation of CC.
!
! Parameters
! ----------
! SE_SOL : Sticking efficiency for   soluble modes [set to 1.0 as in M7]
! SE_INS : Sticking efficiency for insoluble modes [set to 0.3 as in M7]
!
! Inputs
! ------
! NBOX     : Number of grid boxes
! nchemg   : Number of gas phase chemistry tracers
! nbudaer  : Number of aerosol budget fields
! IFUCHS   : Switch for Fuchs (1964) or Fuchs-Sutugin (1971) for CC
! IDCMFP   : Switch : diffusion/mfp  (1=as Gbin v1, 2=as Gbin v1_1)
! ICONDIAM : Switch : wet diam in UKCA_CONDEN (1=g.mean,2=condiam.)
! ND       : Aerosol ptcl no. concentration (ptcls per cc)
! TSQRT    : Square-root of mid-level temperature (K)
! RHOA     : Air density (kg/m3)
! AIRDM3   : Number density of air (per m3)
! DTZ      : Time Step for nucl/cond competition (s)
! WETDP    : Avg wet diameter for each aerosol mode (m)
! PMID     : Centre level pressure (Pa)
! T        : Centre level temperature (K)
!
! InOut
! ------
! MD          : Component median aerosol mass (molecules per ptcl)
! MDT         : Total median aerosol mass (molecules per ptcl)
! GC          : Condensable cpt number density (molecules cm-3)
! BUD_AER_MAS : Aerosol mass fluxes (molecules/cc/DTC)
!
! Outputs:
! -------
! DELGC_COND : Change in vapour conc due to cond (molecules cpt/cm3)
! AGETERM1   : stores mass of soluble material which has condensed onto
!              each of the insoluble modes for use in UKCA_AGEING
!              (molecules per cc)
! S_COND_S   : Condensation sink
!
! Local variables:
! ---------------
! DMOL    : Molecular diameter of condensable cpt (m)
! CC      : Conden. coeff. for condensable cpt onto particle (m^3s^-1)
! RP      : Radius of aerosol particle (m)
! SE      : Sticking efficiency (accomodation coeff)
! SE_SOL  : Sticking efficiency (accomodation coeff) for soluble mode
! SE_INS  : Sticking efficiency (accomodation coeff) for insoluble mode
! MMCG    : Molar mass of condensing gas (kg/mole)
! NC      : Product of number conc and condensation coefficient
! SUMNC   : Sum of NC over all modes
! DELTAMS : Mass of condensing gas taken up by this   soluble mode
! DELTAMI : Mass of condensing gas taken up by this insoluble mode
!    n.b. DELTAMS,DELTAMI all in molecules per cc.
! MASK1-3 : Logical array to define regions of domain to work on
!
! References
! ----------
! Gong et al, JGR, 108(D1), 4007, doi:10.1029/2001JD002002, 2003.
! Raes et al, J. Aerosol Sci., 23 (7), pp. 759--771, 1992.
! Fuchs & Sutugin, Topics in aerosol research, 1971.
! Fuchs, "Mechanics of aerosols", Pergamon Press, 1964.
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! CONC_EPS : Threshold for condensable conc (molecules per cc)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES   : Number of possible aerosol modes
! NMODES_INS: Number of possible insoluble aerosol modes
! NCP      : Number of possible aerosol components
! MODE     : Defines which modes are set
! CONDENSABLE : Logical variable defining which cpts are condensable
! MODESOL  : Defines whether the mode is soluble or not (=1 or 0)
! DIMEN    : Molecular diamters of condensable components (m)
! NUM_EPS  : Value of NEWN below which do not recalculate MD (per cc)
!                                             or carry out process
! CP_SU    : Index of component in which sulfate is stored
! CP_OC    : Index of component in which 1st OC cpt is stored
! CP_SO    : Index of component in which 2nd OC cpt is stored
! SIGMAG   : Geometric standard deviation of mode
! TOPMODE   : Highest number mode for which coag & nucl is done
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! MH2SO4   : Index of MM_GAS, WTRATC and S0G for H2SO4
! MM_GAS   : Array of molar masses for gas phase species (kg/mol)
! DIMEN    : Molecular diamters of condensable components (m)
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
USE ukca_constants,       ONLY: conc_eps

USE ukca_config_specification_mod, ONLY: glomap_variables

USE ukca_mode_setup,      ONLY:                                                &
    nmodes,                                                                    &
    nmodes_ins,                                                                &
    cp_su,                                                                     &
    cp_oc,                                                                     &
    cp_so,                                                                     &
    mode_nuc_sol,                                                              &
    mode_ait_sol,                                                              &
    mode_acc_sol,                                                              &
    mode_cor_sol,                                                              &
    mode_ait_insol,                                                            &
    mode_acc_insol,                                                            &
    mode_cor_insol,                                                            &
    mode_sup_insol

USE ukca_setup_indices,   ONLY: mh2so4, msec_org, msec_orgi,                   &
    nmascondocaccins, nmascondocaccsol, nmascondocaitins,                      &
    nmascondocaitsol, nmascondoccorins, nmascondoccorsol,                      &
    nmascondocnucsol, nmascondocinucsol, nmascondociaitsol,                    &
    nmascondociaccsol, nmascondocicorsol, nmascondociaitins,                   &
    nmascondsoaccins, nmascondsoaccsol,                                        &
    nmascondsoaitins, nmascondsoaitsol, nmascondsocorins,                      &
    nmascondsocorsol, nmascondsonucsol, nmascondsuaccins,                      &
    nmascondsuaccsol, nmascondsuaitins, nmascondsuaitsol,                      &
    nmascondsucorins, nmascondsucorsol, nmascondsunucsol,                      &
    nmascondsusupins, nmascondocsupins, nmascondsosupins, dimen,               &
    condensable, condensable_choice, mm_gas

USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
USE ereport_mod,          ONLY: ereport
USE ukca_types_mod,       ONLY: log_small
USE umPrintMgr, ONLY: umMessage, umPrint
USE ukca_cond_coff_v_mod, ONLY: ukca_cond_coff_v

IMPLICIT NONE

! Subroutine interface:
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: nchemg
INTEGER, INTENT(IN) :: nbudaer
INTEGER, INTENT(IN) :: ifuchs
INTEGER, INTENT(IN) :: idcmfp
INTEGER, INTENT(IN) :: icondiam
REAL, INTENT(IN)    :: nd(nbox,nmodes)
REAL, INTENT(IN)    :: tsqrt(nbox)
REAL, INTENT(IN)    :: rhoa(nbox)
REAL, INTENT(IN)    :: airdm3(nbox)
REAL, INTENT(IN)    :: dtz
REAL, INTENT(IN)    :: wetdp(nbox,nmodes)
REAL, INTENT(IN)    :: pmid(nbox)
REAL, INTENT(IN)    :: t(nbox)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables%ncp)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
REAL, INTENT(IN OUT) :: gc(nbox,nchemg)
REAL, INTENT(IN OUT) :: bud_aer_mas(nbox,0:nbudaer)
REAL, INTENT(OUT)   :: delgc_cond(nbox,nchemg)
REAL, INTENT(OUT)   :: ageterm1(nbox,nmodes_ins,nchemg)
REAL, INTENT(OUT)   :: s_cond_s(nbox)

! Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: modesol(:)
REAL,    POINTER :: num_eps(:)
REAL,    POINTER :: sigmag(:)
INTEGER, POINTER :: topmode

INTEGER :: icp
INTEGER :: imode
INTEGER :: jv
INTEGER :: ierr        ! Error code
LOGICAL (KIND=log_small) :: mask1(nbox)
LOGICAL (KIND=log_small) :: mask2(nbox)
LOGICAL (KIND=log_small) :: mask3(nbox)
LOGICAL (KIND=log_small) :: mask3i(nbox)
LOGICAL (KIND=log_small) :: mask4i(nbox)
REAL    :: dmol
REAL    :: mmcg
REAL    :: cc(nbox)
REAL    :: rp(nbox)
REAL    :: sumnc(nbox)
REAL    :: nc(nbox,nmodes)
REAL    :: deltams(nbox)
REAL    :: deltami(nbox)
REAL    :: se
REAL    :: y2
REAL    :: aa
REAL    :: aa_modes(nmodes)   ! AA specific to each mode
REAL, PARAMETER :: se_sol=1.0
!!      REAL, PARAMETER :: SE_INS=0.3
REAL, PARAMETER :: se_ins=1.0
!
REAL :: sinkarr(nbox)
REAL :: difvol
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CONDEN'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
mode        => glomap_variables%mode
modesol     => glomap_variables%modesol
num_eps     => glomap_variables%num_eps
sigmag      => glomap_variables%sigmag
topmode     => glomap_variables%topmode

ageterm1(:,:,:)=0.0
s_cond_s(:)=0.0
sinkarr(:)=0.0
!
! Set value of AA for each mode according to ICONDIAM
! .. these values are taken from Figure 1 Lehtinen et al (2003)
SELECT CASE (icondiam)
CASE (1)
  aa_modes(:) = 0.0     ! as v1_gm4c, take g.m. number radius
CASE (2)
  aa_modes(1) = 2.0 ! continuum regime -- 2nd radial moment
  aa_modes(2) = 1.9
  aa_modes(3) = 1.5
  aa_modes(4) = 1.1 ! molecular regime -- 1st radial moment
  aa_modes(5) = 1.9
  aa_modes(6) = 1.5
  aa_modes(7) = 1.1
  aa_modes(8) = 1.1
CASE DEFAULT
  ierr = 1
  WRITE(umMessage,'(A,I5)')'Unexpected Value of ICONDIAM ',icondiam
  CALL umPrint(umMessage,src='ukca_conden')
  CALL ereport('UKCA_CONDEN',ierr,'Unexpected ICONDIAM value')
END SELECT

DO jv=1,nchemg
  IF (condensable(jv)) THEN

    ! .. Set component into which component will condense
    icp=condensable_choice(jv)

    dmol=dimen(jv)
    mmcg=mm_gas(jv)
    delgc_cond(:,jv)=0.0

    mask1(:) = gc(:,jv) > conc_eps

    sumnc(:)=0.0
    DO imode=1,topmode
      IF (mode(imode)) THEN

        aa = aa_modes(imode)

        !!          Y2=EXP(2.0*LOG(SIGMAG(IMODE))*LOG(SIGMAG(IMODE)))
        !! original runs tried just setting equivalent to AA=2.0 (above)
        y2=EXP(0.5*aa*aa*LOG(sigmag(imode))*LOG(sigmag(imode)))
        ! now use values of AA for each mode from Lehtinen et al (2003)

        nc(:,imode)=0.0
        mask2(:) = mask1(:) .AND. (nd(:,imode) > num_eps(imode))

        rp(:)=wetdp(:,imode)*0.5*y2 ! radius to use for condensation

        SELECT CASE (modesol(imode))
        CASE (1)
          se=se_sol
        CASE (0)
          se=se_ins
        END SELECT

        !!          IF(JV == MH2SO4  ) DIFVOL=DH2SO4 ! as H2SO4+H2O hydrate
        !!          IF(JV == msec_org) DIFVOL=DSECOR ! as OH-a-pinene radical
        !!          IF(JV == msec_orgi) DIFVOL=DSECOR ! as OH-a-pinene radical
        IF (jv == mh2so4) THEN
          difvol = 51.96          ! values from dist_data (BIN)
        ELSE IF (jv == msec_org .OR. jv == msec_orgi) THEN
          difvol = 204.14         ! values from dist_data (BIN)
        ELSE            ! trap if DIFVOL is being used w/o being defined
          ierr = 100+jv
          CALL ereport('UKCA_CONDEN',ierr,'DIFVOL remains undefined')
        END IF

        ! ..  Calculate change in condensable cpt conc (molecules cm^-3)
        CALL ukca_cond_coff_v(nbox,mask2,rp,tsqrt,airdm3,rhoa,                 &
               mmcg,se,dmol,ifuchs,cc,sinkarr,pmid,t,difvol,idcmfp)
        WHERE (mask2(:))
          nc(:,imode)=nd(:,imode)*cc(:)
          sumnc(:)=sumnc(:)+nc(:,imode)
        END WHERE
        !!          IF(CONDENSABLE_CHOICE(JV).EQ.1) THEN ! if H2SO4
        IF (jv == mh2so4) THEN ! if H2SO4
          s_cond_s(:)=s_cond_s(:)+nd(:,imode)*sinkarr(:)
        END IF

      END IF ! if mode is present
    END DO ! Over modes
    !
    WHERE (mask1(:))                                                           &
     delgc_cond(:,jv)=gc(:,jv)*(1.0-EXP(-sumnc(:)*dtz))

    ! .. Update condensable cpt concentration (molecules cm^-3)
    mask2(:) = mask1(:) .AND. (delgc_cond(:,jv) > conc_eps)
    WHERE (mask2(:) .AND. delgc_cond(:,jv) > gc(:,jv))
      delgc_cond(:,jv)=delgc_cond(:,jv)/gc(:,jv) ! make sure no -ves
    END WHERE
    WHERE (mask2(:)) gc(:,jv)=gc(:,jv)-delgc_cond(:,jv)

    ! loop over sol modes (do cond sol -> ins here too)
    DO imode=mode_nuc_sol,mode_cor_sol
      IF (mode(imode)) THEN

        !         Calculate increase in total & cpt masses in each soluble mode
        deltams(:)=0.0
        deltami(:)=0.0

        mask3 (:) = mask2(:) .AND. (nd(:,imode  ) > num_eps(imode))

        mask3i(:) = .FALSE.
        mask4i(:) = .FALSE.

        IF ( mode(mode_ait_insol) .AND. ( imode == mode_ait_sol ) ) THEN
          mask3i(:) = mask2(:) .AND. ( nd(:,mode_ait_insol) > num_eps(imode) )
        END IF

        IF ( mode(mode_acc_insol) .AND. ( imode == mode_acc_sol ) .AND.        &
             (topmode > mode_ait_insol) ) THEN
          mask3i(:) = mask2(:) .AND. ( nd(:,mode_acc_insol) > num_eps(imode) )
        END IF

        IF ( mode(mode_cor_insol) .AND. ( imode == mode_cor_sol ) .AND.        &
             (topmode > mode_ait_insol) ) THEN
          mask3i(:) = mask2(:) .AND. ( nd(:,mode_cor_insol) > num_eps(imode) )
        END IF

        IF ( mode(mode_sup_insol) .AND. ( imode == mode_cor_sol ) .AND.        &
             (topmode > mode_ait_insol) ) THEN
          mask4i(:) = mask2(:) .AND. ( nd(:,mode_sup_insol) > num_eps(imode) )
        END IF

        IF ( imode == mode_nuc_sol ) THEN

          IF ((icp == cp_su) .AND. (nmascondsunucsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsunucsol)=                                 &
              bud_aer_mas(:,nmascondsunucsol)+deltams(:)

            END WHERE
          END IF

          IF (msec_orgi > 0 .AND. jv == msec_orgi) THEN
            IF ((icp == cp_oc) .AND. (nmascondocinucsol > 0)) THEN
              ! condensation of sec_org_i to nucleation-sol
              WHERE (mask3(:))
                deltams(:) = delgc_cond(:,jv)*nc(:,imode)/sumnc(:)
                bud_aer_mas(:,nmascondocinucsol) =                             &
                bud_aer_mas(:,nmascondocinucsol) + deltams(:)
              END WHERE
            END IF
          ELSE
              ! condensation of sec_org to nucleation-sol
            IF ((icp == cp_oc) .AND. (nmascondocnucsol > 0)) THEN
              WHERE (mask3(:))
                deltams(:) = delgc_cond(:,jv)*nc(:,imode)/sumnc(:)
                bud_aer_mas(:,nmascondocnucsol) =                              &
                bud_aer_mas(:,nmascondocnucsol) + deltams(:)
              END WHERE
            END IF
          END IF

          IF ((icp == cp_so) .AND. (nmascondsonucsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsonucsol)=                                 &
              bud_aer_mas(:,nmascondsonucsol)+deltams(:)

            END WHERE
          END IF

        END IF ! ( imode == mode_nuc_sol )

        IF ( imode == mode_ait_sol ) THEN

          IF ((icp == cp_su) .AND. (nmascondsuaitsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsuaitsol)=                                 &
              bud_aer_mas(:,nmascondsuaitsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_su) .AND. (nmascondsuaitins > 0)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_ait_insol)/sumnc(:)

              bud_aer_mas(:,nmascondsuaitins)=                                 &
              bud_aer_mas(:,nmascondsuaitins)+deltami(:)

              ageterm1(:,mode_nuc_sol,jv)=deltami(:)

            END WHERE
          END IF

          IF (msec_orgi > 0 .AND. jv == msec_orgi) THEN
            IF ((icp == cp_oc) .AND. (nmascondociaitsol > 0)) THEN
              ! condensation of sec_org_i to Aitken-sol
              WHERE (mask3(:))
                deltams(:) = delgc_cond(:,jv)*nc(:,imode)/sumnc(:)
                bud_aer_mas(:,nmascondociaitsol) =                             &
                bud_aer_mas(:,nmascondociaitsol) + deltams(:)
              END WHERE
            END IF
          ELSE
            IF ((icp == cp_oc) .AND. (nmascondocaitsol > 0)) THEN
              ! condensation of sec_org to Aitken-sol
              WHERE (mask3(:))
                deltams(:) = delgc_cond(:,jv)*nc(:,imode)/sumnc(:)
                bud_aer_mas(:,nmascondocaitsol) =                              &
                bud_aer_mas(:,nmascondocaitsol) + deltams(:)
              END WHERE
            END IF
          END IF

          IF (msec_orgi > 0 .AND. jv == msec_orgi) THEN
            IF ((icp == cp_oc) .AND. (nmascondociaitins > 0)) THEN
              ! condensation of sec_org_i to Aitken-insol
              WHERE (mask3i(:))
                deltami(:) = delgc_cond(:,jv)*nc(:,mode_ait_insol)/sumnc(:)
                bud_aer_mas(:,nmascondociaitins) =                             &
                bud_aer_mas(:,nmascondociaitins) + deltami(:)
                ageterm1(:,mode_nuc_sol,jv) = deltami(:)
              END WHERE
            END IF
          ELSE
            IF ((icp == cp_oc) .AND. (nmascondocaitins > 0)) THEN
              ! condensation of sec_org to Aitken-insol
              WHERE (mask3i(:))
                deltami(:) = delgc_cond(:,jv)*nc(:,mode_ait_insol)/sumnc(:)
                bud_aer_mas(:,nmascondocaitins) =                              &
                bud_aer_mas(:,nmascondocaitins) + deltami(:)
                ageterm1(:,mode_nuc_sol,jv) = deltami(:)
              END WHERE
            END IF
          END IF

          IF ((icp == cp_so) .AND. (nmascondsoaitsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsoaitsol)=                                 &
              bud_aer_mas(:,nmascondsoaitsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsoaitins > 0)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_ait_insol)/sumnc(:)

              bud_aer_mas(:,nmascondsoaitins)=                                 &
              bud_aer_mas(:,nmascondsoaitins)+deltami(:)

              ageterm1(:,mode_nuc_sol,jv)=deltami(:)

            END WHERE
          END IF

        END IF ! IF ( imode == mode_ait_sol )

        IF ( imode == mode_acc_sol ) THEN

          IF ((icp == cp_su) .AND. (nmascondsuaccsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsuaccsol)=                                 &
              bud_aer_mas(:,nmascondsuaccsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_su) .AND. (nmascondsuaccins > 0) .AND.                &
              (topmode > mode_ait_insol)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_acc_insol)/sumnc(:)

              bud_aer_mas(:,nmascondsuaccins)=                                 &
              bud_aer_mas(:,nmascondsuaccins)+deltami(:)

              ageterm1(:,mode_ait_sol,jv)=deltami(:)

            END WHERE
          END IF

          IF (msec_orgi > 0 .AND. jv == msec_orgi) THEN
            IF ((icp == cp_oc) .AND. (nmascondociaccsol > 0)) THEN
              ! condensation of sec_org_i to accumulation-sol
              WHERE (mask3(:))
                deltams(:) = delgc_cond(:,jv)*nc(:,imode)/sumnc(:)
                bud_aer_mas(:,nmascondociaccsol) =                             &
                bud_aer_mas(:,nmascondociaccsol) + deltams(:)
              END WHERE
            END IF
          ELSE
            IF ((icp == cp_oc) .AND. (nmascondocaccsol > 0)) THEN
              ! condensation of sec_org to accumulation-sol
              WHERE (mask3(:))
                deltams(:) = delgc_cond(:,jv)*nc(:,imode)/sumnc(:)
                bud_aer_mas(:,nmascondocaccsol) =                              &
                bud_aer_mas(:,nmascondocaccsol) + deltams(:)
              END WHERE
            END IF
          END IF

          IF ((icp == cp_oc) .AND. (nmascondocaccins > 0) .AND.                &
              (topmode > mode_ait_insol)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_acc_insol)/sumnc(:)

              bud_aer_mas(:,nmascondocaccins)=                                 &
              bud_aer_mas(:,nmascondocaccins)+deltami(:)

              ageterm1(:,mode_ait_sol,jv)=deltami(:)

            END WHERE
          END IF


          IF ((icp == cp_oc) .AND. (nmascondocaccins > 0) .AND.                &
              (topmode > mode_ait_insol)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_acc_insol)/sumnc(:)

              bud_aer_mas(:,nmascondocaccins)=                                 &
              bud_aer_mas(:,nmascondocaccins)+deltami(:)

              ageterm1(:,mode_ait_sol,jv)=deltami(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsoaccsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsoaccsol)=                                 &
              bud_aer_mas(:,nmascondsoaccsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsoaccins > 0) .AND.                &
              (topmode > mode_ait_insol)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_acc_insol)/sumnc(:)

              bud_aer_mas(:,nmascondsoaccins)=                                 &
              bud_aer_mas(:,nmascondsoaccins)+deltami(:)

              ageterm1(:,mode_ait_sol,jv)=deltami(:)

            END WHERE
          END IF

        END IF ! IF ( imode == mode_acc_sol )

        IF ( imode == mode_cor_sol ) THEN

          IF ((icp == cp_su) .AND. (nmascondsucorsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsucorsol)=                                 &
              bud_aer_mas(:,nmascondsucorsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_su) .AND. (nmascondsucorins > 0) .AND.                &
              (topmode > mode_ait_insol)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_cor_insol)/sumnc(:)

              bud_aer_mas(:,nmascondsucorins)=                                 &
              bud_aer_mas(:,nmascondsucorins)+deltami(:)

              ageterm1(:,mode_acc_sol,jv)=deltami(:)

            END WHERE
          END IF

          IF ((icp == cp_su) .AND. (nmascondsusupins > 0) .AND.                &
              (topmode > mode_ait_insol)) THEN
            WHERE (mask4i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_sup_insol)/sumnc(:)

              bud_aer_mas(:,nmascondsusupins)=                                 &
              bud_aer_mas(:,nmascondsusupins)+deltami(:)

              ageterm1(:,mode_cor_sol,jv)=deltami(:)

            END WHERE
          END IF

          IF (msec_orgi > 0 .AND. jv == msec_orgi) THEN
            IF ((icp == cp_oc) .AND. (nmascondocicorsol > 0)) THEN
              ! condensation of sec_org_i to coarse-sol
              WHERE (mask3(:))
                deltams(:) = delgc_cond(:,jv)*nc(:,imode)/sumnc(:)
                bud_aer_mas(:,nmascondocicorsol) =                             &
                bud_aer_mas(:,nmascondocicorsol) + deltams(:)
              END WHERE
            END IF
          ELSE
            IF ((icp == cp_oc) .AND. (nmascondoccorsol > 0)) THEN
              ! condensation of sec_org to coarse-sol
              WHERE (mask3(:))
                deltams(:) = delgc_cond(:,jv)*nc(:,imode)/sumnc(:)
                bud_aer_mas(:,nmascondoccorsol) =                              &
                bud_aer_mas(:,nmascondoccorsol) + deltams(:)
              END WHERE
            END IF
          END IF

          IF ((icp == cp_oc) .AND. (nmascondoccorins > 0) .AND.                &
              (topmode > mode_ait_insol)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_cor_insol)/sumnc(:)

              bud_aer_mas(:,nmascondoccorins)=                                 &
              bud_aer_mas(:,nmascondoccorins)+deltami(:)

              ageterm1(:,mode_acc_sol,jv)=deltami(:)

            END WHERE
          END IF

          IF ((icp == cp_oc) .AND. (nmascondocsupins > 0) .AND.                &
              (topmode > mode_ait_insol)) THEN
            WHERE (mask4i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_sup_insol)/sumnc(:)

              bud_aer_mas(:,nmascondocsupins)=                                 &
              bud_aer_mas(:,nmascondocsupins)+deltami(:)

              ageterm1(:,mode_cor_sol,jv)=deltami(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsocorsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsocorsol)=                                 &
              bud_aer_mas(:,nmascondsocorsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsocorins > 0) .AND.                &
              (topmode > mode_ait_insol)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_cor_insol)/sumnc(:)

              bud_aer_mas(:,nmascondsocorins)=                                 &
              bud_aer_mas(:,nmascondsocorins)+deltami(:)

              ageterm1(:,mode_acc_sol,jv)=deltami(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsosupins > 0) .AND.                &
              (topmode > mode_ait_insol)) THEN
            WHERE (mask4i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,mode_sup_insol)/sumnc(:)

              bud_aer_mas(:,nmascondsosupins)=                                 &
              bud_aer_mas(:,nmascondsosupins)+deltami(:)

              ageterm1(:,mode_cor_sol,jv)=deltami(:)

            END WHERE
          END IF

        END IF ! ( imode == mode_cor_sol )

        WHERE (mask3(:))
          md(:,imode,icp)=                                                     &
           (md(:,imode,icp)*nd(:,imode)+deltams(:))/nd(:,imode)
          mdt(:,imode)=                                                        &
           (mdt(:,imode)*nd(:,imode)+deltams(:))/nd(:,imode)
        END WHERE

        !!
        !! **** HERE HAVE COMMENTED OUT THE UPDATING OF THE SOLUBLE MODE
        !! **** MD,MDT FOR THE CONDENSATION ONTO INSOLUBLE MODES BECAUSE IN
        !! **** THE CASE WHERE THERE ARE NO PARTICLES IN THE AITKEN SOLUBLE
        !! **** MODE BUT THERE ARE IN THE AITKEN INSOLUBLE MODE, THE UPDATING
        !! **** OF THE SOLUBLE MODE MD,MDT WILL NOT BE POSSIBLE UNTIL THE
        !! **** SOLUBLE MODE ND HAS BEEN UPDATED DUE TO THE TRANSFER OF AGED
        !! **** PARTICLES. THIS IS DONE IN THE UKCA_AGEING ROUTINE -- SO THIS IS
        !! **** WHERE THE UPDATE OF THE SOLUBLE MODE MD,MDT FOR THE CONDENSATION
        !! **** ONTO THE SHOULD BE DONE ALSO.
        !!
      END IF ! if mode is present

    END DO ! IMODE=1,4 (soluble modes)

  END IF ! IF CONDENSABLE(JV)

END DO ! DO jv=1,nchemg

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_conden
END MODULE ukca_conden_mod
