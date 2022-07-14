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
!    Calculates in-cloud aerosol wet deposition (nucleation-scavenging).
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and The Met Office.
!  See: www.ukca.ac.uk
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
MODULE ukca_rainout_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_RAINOUT_MOD'

CONTAINS

SUBROUTINE ukca_rainout(nbox,nbudaer,nd,md,mdt,fconv_conv,                     &
                        crain,crain_up,clwc,clf, autoconv1d,t,dtc,             &
                        bud_aer_mas,inucscav, lcvrainout,wetdp)
!---------------------------------------------------------------------------
!
!     Calculates in-cloud aerosl wet deposition (nucleation-scavenging)
!
!     Includes large- (dynamic) & small- scale (convective) precip.
!
!     Include conversion rates FCONV_CONV and FCONV_DYN which represent
!     fraction of condensate this is converted to rain in 6 hours.
!
!     Currently takes FCONV_CONV as input with FCONV_DYN set constant
!
!     Inputs
!     ------
!     NBOX       : Number of grid boxes
!     nbudaer    : Number of aerosol budget fields
!     ND         : Initial no. concentration of aer mode (ptcls/cm3)
!     MD         : Initial avg cpt mass conc of aer mode (molcules/ptcl)
!     MDT        : Avg tot mass of aerosol ptcl in mode (particle^-1)
!     FCONV_CONV : Fraction of box condensate converted to rain in
!                                                   6 hours (convective)
!     CRAIN      : Rain rate for conv. precip. in box (kgm^-2s^-1)
!     CRAIN_UP   : Rain rate for conv. precip. in box above (kgm^-2s^-1)
!     T          : Air temperature (K)
!     DTC        : Chemistry timestep (s)
!     INUCSCAV   : Switch for scheme for removal by nucl scav
!     LCVRAINOUT : Switch for convective rainout (logical) normally set to FALSE
!                : as plume scavenging is selected and convective scavenging is
!                : done within the convection scheme.
!     CLF        : Liquid cloud fraction
!     autoconv1d : Autoconversion rate including accretion (kg.kg^-1.s^-1)
!     CLWC       : Cloud liquid water content (kg/kg)
!     WETDP      : Geometric mean wet diameter for each mode (m)
!
!     Outputs
!     -------
!     ND           : Updated no. concentration in each mode (ptcls/cm3)
!     BUD_AER_MASS : Aerosol mass budgets
!
!     Local Variables
!     ---------------
!     SWITCH       : Rain created in box? (0=none,1=conv,2=dyn,3=both=3)
!     RSCAV        : Scavenging parameters for each mode
!     FCONV_DYN    : Fraction of box condensate converted to rain in
!                                                   6 hours (dynamic)
!     TAU_CONV_DYN : e-folding timescale for conversion of
!                                        condensate to rain (dynamic)
!     TAU_CONV_CONV: e-folding timescale for conversion of
!                                        condensate to rain (convective)
!     FBOX_CONV    : Gridbox fraction over which convective rain occurs
!     TICE         : Temperature below which insoluble aerosol can act
!                    as ice nucleii and hence be removed
!     DELN         : Change in number conc. due to nucleation-scavenging
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES       : Number of possible aerosol modes
!     NCP          : Number of possible aerosol components
!     MODE         : Defines which modes are set
!     COMPONENT    : Defines which components are set in each mode
!     SIGMAG       : Geometric standard deviation of mode
!     MODESOL      : Defines whether mode is soluble of not (integer)
!     NUM_EPS      : Value of NEWN below which don't recalculate MD or
!                                                    carry out process
!     CP_SU        : Index of component in which H2SO4  cpt is stored
!     CP_BC        : Index of component in which BC     cpt is stored
!     CP_OC        : Index of component in which 1st OC cpt is stored
!     CP_CL        : Index of component in which NaCl   cpt is stored
!     CP_DU        : Index of component in which dust   cpt is stored
!     CP_SO        : Index of component in which 2nd OC cpt is stored
!
!     Inputted by module UKCA_SETUP_INDICES
!     -------------------------------------
!     Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
USE ukca_mode_setup, ONLY:                                                     &
    nmodes,                                                                    &
    ncp,                                                                       &
    mode,                                                                      &
    component,                                                                 &
    sigmag,                                                                    &
    modesol,                                                                   &
    num_eps,                                                                   &
    cp_su,                                                                     &
    cp_bc,                                                                     &
    cp_oc,                                                                     &
    cp_cl,                                                                     &
    cp_du,                                                                     &
    cp_so,                                                                     &
    mode_nuc_sol,                                                              &
    mode_ait_sol,                                                              &
    mode_acc_sol,                                                              &
    mode_cor_sol,                                                              &
    mode_ait_insol,                                                            &
    mode_acc_insol,                                                            &
    mode_cor_insol,                                                            &
    cp_no3, cp_nh4, cp_nn

USE ukca_setup_indices,  ONLY: nmasnuscsunucsol,                               &
                               nmasnuscsuaitsol, nmasnuscsuaccsol,             &
                               nmasnuscsucorsol, nmasnuscbcaitsol,             &
                               nmasnuscbcaccsol, nmasnuscbccorsol,             &
                               nmasnuscbcaitins, nmasnuscocnucsol,             &
                               nmasnuscocaitsol, nmasnuscocaccsol,             &
                               nmasnuscoccorsol, nmasnuscocaitins,             &
                               nmasnuscssaccsol, nmasnuscsscorsol,             &
                               nmasnuscsonucsol, nmasnuscsoaitsol,             &
                               nmasnuscsoaccsol, nmasnuscsocorsol,             &
                               nmasnuscduaccsol, nmasnuscducorsol,             &
                               nmasnuscduaccins, nmasnuscducorins,             &
                               nmasnuscntaitsol, nmasnuscntaccsol,             &
                               nmasnuscntcorsol, nmasnuscnhaitsol,             &
                               nmasnuscnhaccsol, nmasnuscnhcorsol,             &
                               nmasnuscnnaccsol, nmasnuscnncorsol

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umErf_mod, ONLY: umErf

IMPLICIT NONE

! .. Subroutine interface
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: nbudaer
INTEGER, INTENT(IN) :: inucscav
REAL, INTENT(IN)    :: fconv_conv(nbox)
REAL, INTENT(IN)    :: crain(nbox)
REAL, INTENT(IN)    :: crain_up(nbox)
REAL, INTENT(IN)    :: clwc(nbox)
REAL, INTENT(IN)    :: clf(nbox)
REAL, INTENT(IN)    :: autoconv1d(nbox)
REAL, INTENT(IN)    :: t(nbox)
REAL, INTENT(IN)    :: dtc
REAL, INTENT(IN)    :: wetdp(nbox,nmodes)
LOGICAL, INTENT(IN) :: lcvrainout

REAL, INTENT(IN OUT) :: nd (nbox,nmodes)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
REAL, INTENT(IN OUT) :: md (nbox,nmodes,ncp)
REAL, INTENT(IN OUT) :: bud_aer_mas(nbox,0:nbudaer)

! .. Local variables
INTEGER, PARAMETER :: i_none = 0     ! defines precipitation type
INTEGER, PARAMETER :: i_conv = 1
INTEGER, PARAMETER :: i_dynm = 2
INTEGER, PARAMETER :: i_both = 3
INTEGER :: jl
INTEGER :: imode
INTEGER :: icp
INTEGER :: switch                    ! precipitation type
REAL    :: rscav(nmodes)
REAL    :: rscavm(nmodes)
REAL    :: deln
REAL    :: deln0
REAL    :: deln0_conv                 ! stores convective deln0
REAL    :: ndnew
REAL    :: tau_conv_conv
REAL    :: tau_conv_dyn
REAL    :: beta
REAL    :: craindiff
REAL, PARAMETER :: nscavact = 103.0e-9     ! Activation dry radius
REAL, PARAMETER :: fconv_dyn=0.9999
REAL, PARAMETER :: fbox_conv=0.3
REAL, PARAMETER :: rain_frac=0.3
REAL, PARAMETER :: tice=258.0

REAL    :: lnratn
REAL    :: erfnum
REAL    :: frac_n
REAL    :: log2sg
REAL    :: lnratm
REAL    :: erfmas
REAL    :: frac_m
REAL    :: dm(ncp)
REAL    :: scav1
REAL    :: scav2
REAL    :: scav1m
REAL    :: scav2m

REAL    :: dp
REAL    :: dp2
REAL    :: maxfracn

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RAINOUT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

maxfracn=0.90
!
! ..  Calculate timescale for conversion to dynamic rain assuming
! ..  a factor FCONV_DYN is converted to rain in 6 hours
tau_conv_dyn=(-6.0*3600.0)/(LOG(1.0-fconv_dyn))

!! now follow approach implemented at start of AEROS to better mimic
!! GLOMAP-bin nucleation-scavenging approach with size threshold.

DO jl=1,nbox
  !
  IF ((inucscav == 1) .OR. (inucscav == 3)) THEN


    imode = mode_ait_sol

    dp=wetdp(jl,imode)


    IF (nd(jl,imode) > num_eps(imode)) THEN
      lnratn=LOG(nscavact*2.0/dp)
      erfnum=lnratn/SQRT(2.0)/LOG(sigmag(imode))
      frac_n=0.5*(1.0+umErf(erfnum)) ! fraction remaining in mode

      IF (frac_n > maxfracn) THEN ! if less than 1% of number is
                                 ! larger than size then don't remove
        scav1 =0.0
        scav1m=0.0
      ELSE
        scav1 =1.0-frac_n

        log2sg=LOG(sigmag(imode))*LOG(sigmag(imode))
        dp2=EXP(LOG(dp)+3.0*log2sg) ! volume median diameter

        ! .. 2nd threshold determines fraction of number/mass
        lnratm=LOG(nscavact*2.0/dp2)
        erfmas=lnratm/SQRT(2.0)/LOG(sigmag(imode))
        frac_m=0.5*(1.0+umErf(erfmas))
        scav1m=1.0-frac_m
      END IF
    ELSE
      scav1 =0.0
      scav1m=0.0
    END IF  ! IF(ND(JL,IMODE) > NUM_EPS(IMODE))

    imode = mode_acc_sol

    dp=wetdp(jl,imode)


    IF (nd(jl,imode) > num_eps(imode)) THEN
      lnratn=LOG(nscavact*2.0/dp)
      erfnum=lnratn/SQRT(2.0)/LOG(sigmag(imode))
      frac_n=0.5*(1.0+umErf(erfnum)) ! fraction remaining in mode
      IF (frac_n > maxfracn) THEN ! if less than 10% of number is
                                 ! larger than size then don't remove
        scav2 =0.0
        scav2m=0.0
      ELSE
        scav2 =1.0-frac_n

        log2sg=LOG(sigmag(imode))*LOG(sigmag(imode))
        dp2=EXP(LOG(dp)+3.0*log2sg) ! volume median diameter

        ! .. 2nd threshold determines fraction of number/mass
        lnratm=LOG(nscavact*2.0/dp2)
        erfmas=lnratm/SQRT(2.0)/LOG(sigmag(imode))
        frac_m=0.5*(1.0+umErf(erfmas))
        ! .. limit DELM to be max 99.9%
        IF (frac_m < 0.001) frac_m=0.001

        scav2m=1.0-frac_m
      END IF

    ELSE
      scav2 =0.0
      scav2m=0.0
    END IF  ! IF(ND(JL,IMODE) > NUM_EPS(IMODE))

    rscav =[0.00,scav1 ,scav2 ,1.00,1.00,1.00,1.00]
    rscavm=[0.00,scav1m,scav2m,1.00,1.00,1.00,1.00]

    IF (inucscav == 3) THEN
      rscav(6)=0.0
      rscav(7)=0.0
      rscavm(6)=0.0
      rscavm(7)=0.0
    END IF  !inucscav=3

  END IF ! INUCSCAV=1 or 3

  IF (inucscav == 2) THEN ! set number and mass same as M7
    rscav =[0.10,0.25,0.85,0.99,0.20,0.40,0.40]
    rscavm=[0.10,0.25,0.85,0.99,0.20,0.40,0.40]
  END IF

  ! Determine precipitation type, and activate switch above threshold:
  !   (craindiff > 1e-10 for convective; autoconv1d > 1.0e-10 for dynamic)
  switch = i_none
  IF (lcvrainout) THEN
    craindiff = crain(jl) - crain_up(jl)
    IF (craindiff > 1.0e-10) switch = i_conv
  END IF
  IF (autoconv1d(jl) > 1.0e-10) THEN
    IF (switch == i_none) switch = i_dynm
    IF (switch == i_conv) switch = i_both
  END IF

  ! .. Switch > 0 only if rain is FORMED in that level
  IF (switch > i_none) THEN

    DO imode=1,nmodes
      IF (mode(imode)) THEN
        IF (nd(jl,imode) > num_eps(imode)) THEN

          ! .. Only apply to soluble modes or insoluble when T<TICE
          IF ((modesol(imode) == 1) .OR. (t(jl) < tice)) THEN

            !------------------------------------------------------------------
            !
            ! This section does removal by small-scale (conv) precipitation
            ! .. Convective rain : ND -> ND*(1-FCONV_CONV) over 6 hours
            !                      only apply over fraction FBOX_CONV
            !
            IF (switch == i_conv .OR. switch == i_both) THEN
              !  If convective rain
              IF (fconv_conv(jl) < 1.0) THEN
                tau_conv_conv=(-6.0*3600.0)/(LOG(1.0-fconv_conv(jl)))
                deln0=fbox_conv*nd(jl,imode)*(1.0-EXP(-dtc/tau_conv_conv))
              ELSE
                deln0=fbox_conv*nd(jl,imode)
              END IF
            END IF ! if small-scale (convective) rain

            ! When both dynamic and convective present, store deln0
            IF (switch == i_both) deln0_conv = deln0

            !------------------------------------------------------------------
            !
            ! This section does removal by large-scale (dyn.) precipitation
            !
            ! .. Dynamic rain    : ND -> ND*(1-FCONV_DYN) over 6 hours
            !                      apply over all of box
            !
            IF ((switch == i_dynm) .OR. (switch == i_both)) THEN
              ! if dynamic rain
              IF ( clwc(jl) < 1.0e-10 ) THEN
                ! Old method retained for low CLWC
                deln0=clf(jl)*nd(jl,imode)*(1.0-EXP(-dtc/tau_conv_dyn ))
              ELSE
                beta = autoconv1d(jl)/clwc(jl)
                ! Nucleation scavenging only occurs in liquid cloud and where
                ! rain is produced -> scale expression by CLF
                ! and an assumed raining fraction
                deln0 = rain_frac*clf(jl)*nd(jl,imode)*(1.0-EXP(-dtc*beta))
              END IF

            END IF ! if large-scale (dynamic) rain

            !------------------------------------------------------------------

            ! Add convective and dynamic contributions together
            IF (switch == i_both) THEN
              deln0 = deln0 + deln0_conv

              ! restrict deln0 so that a maximum of 90% of the aerosol number
              !  can be removed
              IF ((nd(jl,imode) - deln0*rscav(imode)) <                        &
                  0.1*nd(jl,imode)) deln0 = 0.9*nd(jl,imode)/rscav(imode)
            END IF       ! switch == i_both

            ! .. Multiply DELN0 by scavenging coefficient for delta for number
            !     (DELN)
            !    set DM for each aerosol component based on RSCAVM value
            !    Note, if INUCSCAV=1, then  RSCAV & RSCAVM are calculated based
            !     on size

            deln=deln0*rscav(imode)

            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                ! .. calculate cpt mass conc to transfer to next mode
                !     (use old no/mass)
                dm(icp)=deln0*md(jl,imode,icp)*rscavm(imode)
              END IF
            END DO

            !------------------------------------------------------------------

            ! .. calculate updated number concentration due to
            !     nucleation-scavenging
            ndnew=nd(jl,imode)-deln

            IF (ndnew > num_eps(imode)) THEN

              ! .. first remove mass to be transferred from mode IMODE
              mdt(jl,imode)=0.0
              DO icp=1,ncp
                IF (component(imode,icp)) THEN
                  md(jl,imode,icp)=                                            &
                      (nd(jl,imode)*md(jl,imode,icp)-dm(icp))/ndnew
                  mdt(jl,imode)=mdt(jl,imode)+md(jl,imode,icp)
                ELSE
                  md(jl,imode,icp)=0.0
                END IF ! COMPONENT(IMODE,ICP)
              END DO

              ! .. update number concentration following nucleation-scavening
              nd(jl,imode)=ndnew

              !----------------------------------------------------------------
              !
              ! .. This section stores removal of each cpt mass for budget
              !     calculations

              DO icp=1,ncp
                IF (component(imode,icp)) THEN
                  IF (icp == cp_su) THEN
                    IF ((imode == mode_nuc_sol) .AND. (nmasnuscsunucsol > 0))  &
                      bud_aer_mas(jl,nmasnuscsunucsol)=                        &
                      bud_aer_mas(jl,nmasnuscsunucsol)+dm(icp)
                    IF ((imode == mode_ait_sol) .AND. (nmasnuscsuaitsol > 0))  &
                      bud_aer_mas(jl,nmasnuscsuaitsol)=                        &
                      bud_aer_mas(jl,nmasnuscsuaitsol)+dm(icp)
                    IF ((imode == mode_acc_sol) .AND. (nmasnuscsuaccsol > 0))  &
                      bud_aer_mas(jl,nmasnuscsuaccsol)=                        &
                      bud_aer_mas(jl,nmasnuscsuaccsol)+dm(icp)
                    IF ((imode == mode_cor_sol) .AND. (nmasnuscsucorsol > 0))  &
                      bud_aer_mas(jl,nmasnuscsucorsol)=                        &
                      bud_aer_mas(jl,nmasnuscsucorsol)+dm(icp)
                  END IF
                  IF (icp == cp_bc) THEN
                    IF ((imode == mode_ait_sol) .AND. (nmasnuscbcaitsol > 0))  &
                      bud_aer_mas(jl,nmasnuscbcaitsol)=                        &
                      bud_aer_mas(jl,nmasnuscbcaitsol)+dm(icp)
                    IF ((imode == mode_acc_sol) .AND. (nmasnuscbcaccsol > 0))  &
                      bud_aer_mas(jl,nmasnuscbcaccsol)=                        &
                      bud_aer_mas(jl,nmasnuscbcaccsol)+dm(icp)
                    IF ((imode == mode_cor_sol) .AND. (nmasnuscbccorsol > 0))  &
                      bud_aer_mas(jl,nmasnuscbccorsol)=                        &
                      bud_aer_mas(jl,nmasnuscbccorsol)+dm(icp)
                    IF ((imode == mode_ait_insol) .AND. (nmasnuscbcaitins >0)) &
                      bud_aer_mas(jl,nmasnuscbcaitins)=                        &
                      bud_aer_mas(jl,nmasnuscbcaitins)+dm(icp)
                  END IF
                  IF (icp == cp_oc) THEN
                    IF ((imode == mode_nuc_sol) .AND. (nmasnuscocnucsol > 0))  &
                      bud_aer_mas(jl,nmasnuscocnucsol)=                        &
                      bud_aer_mas(jl,nmasnuscocnucsol)+dm(icp)
                    IF ((imode == mode_ait_sol) .AND. (nmasnuscocaitsol > 0))  &
                      bud_aer_mas(jl,nmasnuscocaitsol)=                        &
                      bud_aer_mas(jl,nmasnuscocaitsol)+dm(icp)
                    IF ((imode == mode_acc_sol) .AND. (nmasnuscocaccsol > 0))  &
                      bud_aer_mas(jl,nmasnuscocaccsol)=                        &
                      bud_aer_mas(jl,nmasnuscocaccsol)+dm(icp)
                    IF ((imode == mode_cor_sol) .AND. (nmasnuscoccorsol > 0))  &
                      bud_aer_mas(jl,nmasnuscoccorsol)=                        &
                      bud_aer_mas(jl,nmasnuscoccorsol)+dm(icp)
                    IF ((imode == mode_ait_insol) .AND. (nmasnuscocaitins >0)) &
                      bud_aer_mas(jl,nmasnuscocaitins)=                        &
                      bud_aer_mas(jl,nmasnuscocaitins)+dm(icp)
                  END IF
                  IF (icp == cp_cl) THEN
                    IF ((imode == mode_acc_sol) .AND. (nmasnuscssaccsol > 0))  &
                      bud_aer_mas(jl,nmasnuscssaccsol)=                        &
                      bud_aer_mas(jl,nmasnuscssaccsol)+dm(icp)
                    IF ((imode == mode_cor_sol) .AND. (nmasnuscsscorsol > 0))  &
                      bud_aer_mas(jl,nmasnuscsscorsol)=                        &
                      bud_aer_mas(jl,nmasnuscsscorsol)+dm(icp)
                  END IF
                  IF (icp == cp_so) THEN
                    IF ((imode == mode_nuc_sol) .AND. (nmasnuscsonucsol > 0))  &
                      bud_aer_mas(jl,nmasnuscsonucsol)=                        &
                      bud_aer_mas(jl,nmasnuscsonucsol)+dm(icp)
                    IF ((imode == mode_ait_sol) .AND. (nmasnuscsoaitsol > 0))  &
                      bud_aer_mas(jl,nmasnuscsoaitsol)=                        &
                      bud_aer_mas(jl,nmasnuscsoaitsol)+dm(icp)
                    IF ((imode == mode_acc_sol) .AND. (nmasnuscsoaccsol > 0))  &
                      bud_aer_mas(jl,nmasnuscsoaccsol)=                        &
                      bud_aer_mas(jl,nmasnuscsoaccsol)+dm(icp)
                    IF ((imode == mode_cor_sol) .AND. (nmasnuscsocorsol > 0))  &
                      bud_aer_mas(jl,nmasnuscsocorsol)=                        &
                      bud_aer_mas(jl,nmasnuscsocorsol)+dm(icp)
                  END IF
                  IF (icp == cp_du) THEN
                    IF ((imode == mode_acc_sol) .AND. (nmasnuscduaccsol > 0))  &
                      bud_aer_mas(jl,nmasnuscduaccsol)=                        &
                      bud_aer_mas(jl,nmasnuscduaccsol)+dm(icp)
                    IF ((imode == mode_cor_sol) .AND. (nmasnuscducorsol > 0))  &
                      bud_aer_mas(jl,nmasnuscducorsol)=                        &
                      bud_aer_mas(jl,nmasnuscducorsol)+dm(icp)
                    IF ((imode == mode_acc_insol) .AND. (nmasnuscduaccins >0)) &
                      bud_aer_mas(jl,nmasnuscduaccins)=                        &
                      bud_aer_mas(jl,nmasnuscduaccins)+dm(icp)
                    IF ((imode == mode_cor_insol) .AND. (nmasnuscducorins >0)) &
                      bud_aer_mas(jl,nmasnuscducorins)=                        &
                      bud_aer_mas(jl,nmasnuscducorins)+dm(icp)
                  END IF
                  IF (icp == cp_no3) THEN
                    IF ((imode == mode_ait_sol) .AND. (nmasnuscntaitsol > 0))  &
                      bud_aer_mas(jl,nmasnuscntaitsol)=                        &
                      bud_aer_mas(jl,nmasnuscntaitsol)+dm(icp)
                    IF ((imode == mode_acc_sol) .AND. (nmasnuscntaccsol > 0))  &
                      bud_aer_mas(jl,nmasnuscntaccsol)=                        &
                      bud_aer_mas(jl,nmasnuscntaccsol)+dm(icp)
                    IF ((imode == mode_cor_sol) .AND. (nmasnuscntcorsol > 0))  &
                      bud_aer_mas(jl,nmasnuscntcorsol)=                        &
                      bud_aer_mas(jl,nmasnuscntcorsol)+dm(icp)
                  END IF
                  IF (icp == cp_nh4) THEN
                    IF ((imode == mode_ait_sol) .AND. (nmasnuscnhaitsol > 0))  &
                      bud_aer_mas(jl,nmasnuscnhaitsol)=                        &
                      bud_aer_mas(jl,nmasnuscnhaitsol)+dm(icp)
                    IF ((imode == mode_acc_sol) .AND. (nmasnuscnhaccsol > 0))  &
                      bud_aer_mas(jl,nmasnuscnhaccsol)=                        &
                      bud_aer_mas(jl,nmasnuscnhaccsol)+dm(icp)
                    IF ((imode == mode_cor_sol) .AND. (nmasnuscnhcorsol > 0))  &
                      bud_aer_mas(jl,nmasnuscnhcorsol)=                        &
                      bud_aer_mas(jl,nmasnuscnhcorsol)+dm(icp)
                  END IF
                  IF (icp == cp_nn) THEN
                    IF ((imode == mode_acc_sol) .AND. (nmasnuscnnaccsol > 0))  &
                      bud_aer_mas(jl,nmasnuscnnaccsol)=                        &
                      bud_aer_mas(jl,nmasnuscnnaccsol)+dm(icp)
                    IF ((imode == mode_cor_sol) .AND. (nmasnuscnncorsol > 0))  &
                      bud_aer_mas(jl,nmasnuscnncorsol)=                        &
                      bud_aer_mas(jl,nmasnuscnncorsol)+dm(icp)
                  END IF

                END IF ! if component(imode,icp)
              END DO ! loop over components

            END IF ! if ND>NUM_EPS(IMODE)

            !------------------------------------------------------------------

          END IF ! if mode is soluble or T<TICE

        END IF ! IF ND>NUM_EPS
      END IF ! IF MODE is switched on
    END DO ! Loop over modes

  END IF ! IF RAIN IS PRODUCED IN THIS LEVEL

END DO ! Loop over gridboxes

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_rainout
END MODULE ukca_rainout_mod
