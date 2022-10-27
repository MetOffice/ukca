! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description:  West scheme for cloud droplet activation.
!                 doi:10.5194/acp-14-6369-2014
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge,
!  University of Leeds, University of Oxford, and The Met Office.
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
!
!  West scheme for cloud droplet activation.
!  based on *activ_box* activation box model, but no longer a box model!
!
!  Purpose: calls aero_activ_abdulrazzak_ghan which calculates
!  -------- number concentration of aerosol particles which become
!           "activated" into cloud droplets
!
! THINGS AERO_ACTIV_ABDULRAZZAK_GHAN NEEDS:
!
! kbdim
! model_levels
! zesw(kbdim,model_levels)                 saturation water vapour pressure
! zrho(kbdim,model_levels)                 air density[kg m-3]
! zn(kbdim,model_levels,nmodes)            aerosol number concentration for
!                                            each mode [m-3]
! zxtm1(kbdim,model_levels,nmodes,ncp)     tracer mass mixing ratio [kg kg-1]
! ztm1(kbdim,model_levels)                 temperature[T]
! zapm1(kbdim,model_levels)                pressure [K]
! zqm1(kbdim,model_levels)                 specific humidity[kg kg-1]
! zrdry(kbdim,model_levels,nmod)           dry count median radius (m)
! nwbins                                   number of 'bins' for calculating pdf
! zwarr(kbdim,model_levels,nwbins)         lin array of vert vel [m s-1]
! zwpdf(kbdim,model_levels,nwbins)         lin array of pdf of w
! zwbin(kbdim,model_levels,nwbins)         w bin width [m s-1]
! zsmax(kbdim,model_levels)                maximum supersaturation [fraction]
! zwchar(kbdim,model_levels)               calculated characteristic
!                                            updraught [m s-1]
! zcdncactm(kbdim,model_levels,nmodes)     number concentration of activate
! zcdncact(kbdim,model_levels,nmodes)      number concentration of
!                                            activated particles
!
! --------------------------------------------------------------------------
MODULE ukca_activate_mod

IMPLICIT NONE

PRIVATE

PUBLIC  :: ukca_activate

! Pre-calculated fraction two_thirds. Calculated as in the UM module
! bl_option_mod where it was originally defined to avoid changing results
! at the time of refactoring.
REAL, PARAMETER :: one_third = 1.0/3.0
REAL, PARAMETER :: two_thirds = 2.0 * one_third

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_ACTIVATE_MOD'

CONTAINS

SUBROUTINE ukca_activate ( row_length,                                         &
                           rows,                                               &
                           model_levels,                                       &
                           bl_levels,                                          &
                           kbdim,                                              &
                           n_mode_tracers,                                     &
                           n_mode_diags,                                       &
                           nwbins,                                             &
                           p_theta_levels,                                     &
                           t_theta_levels,                                     &
                           q,                                                  &
                           qsvp,                                               &
                           bl_tke,                                             &
                           vertvel,                                            &
                           liq_cloud_frac,                                     &
                           cloud_liq_water,                                    &
                           cdncflag,                                           &
                           n_activ_sum,                                        &
                           cdncwt,                                             &
                           mode_tracers,                                       &
                           mode_diags                                          &
                           )

USE ereport_mod,               ONLY:                                           &
    ereport

USE errormessagelength_mod,    ONLY:                                           &
    errormessagelength

USE missing_data_mod,          ONLY:                                           &
    imdi

USE parkind1,                  ONLY:                                           &
    jprb,                                                                      &
    jpim

USE ukca_abdulrazzak_ghan_mod, ONLY:                                           &
    ukca_abdulrazzak_ghan

USE ukca_activ_mod,            ONLY:                                           &
    activmklin,                                                                &
    activmkskew

USE chemistry_constants_mod,   ONLY:                                           &
    boltzmann                           ! Boltzmann's constant

USE ukca_um_legacy_mod,        ONLY:                                           &
    stashcode_glomap_sec,                                                      &
    ukcaD1codes,                                                               &
    item1_mode_diags,                                                          &
    Nukca_D1items,                                                             &
    l_ukca_mode_diags,                                                         &
    rgas => r,                                                                 &
    i_gc_activation_arg,                                                       &
    i_glomap_clim_activation_scheme

USE ukca_config_defs_mod,      ONLY:                                           &
    nmax_mode_diags

USE ukca_fixeds_mod,           ONLY:                                           &
    ukca_fixeds

USE ukca_mode_setup,           ONLY:                                           &
    nmodes,                                                                    &
    ncp,                                                                       &
    mode,                                                                      &
    modesol

USE ukca_config_specification_mod, ONLY:                                       &
    glomap_config,                                                             &
    i_ukca_activation_arg

USE umPrintMgr,                ONLY:                                           &
    umPrint,                                                                   &
    umMessage,                                                                 &
    PrintStatus,                                                               &
    PrStatus_Oper,                                                             &
    PrStatus_Diag,                                                             &
    newline

USE yomhook,                   ONLY:                                           &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Interface:

INTEGER, INTENT(IN) :: row_length        ! array dimension
INTEGER, INTENT(IN) :: rows              ! array dimension
INTEGER, INTENT(IN) :: model_levels      ! array dimension
INTEGER, INTENT(IN) :: bl_levels         ! array dimension
INTEGER, INTENT(IN) :: kbdim             ! = theta_field_size=row_length*rows
INTEGER, INTENT(IN) :: n_mode_tracers    ! # of mode tracers
INTEGER, INTENT(IN) :: n_mode_diags      ! # of mode diagnostics

! number of 'bins' or 'class intervals' for calculating pdf. Has to be >= 1
INTEGER, INTENT(IN)  :: nwbins

! Pressure aka p_theta_levels
REAL, INTENT(IN)  :: p_theta_levels(row_length,rows,model_levels)

! Temperature aka t_theta_levels
REAL, INTENT(IN)  :: t_theta_levels(row_length,rows,model_levels)

! Specific humidity aka q (e.g. from Section 0 Item 10 in UM parent model)
REAL, INTENT(IN)  :: q(row_length,rows,model_levels)

! Sat vap pressure with respect to liquid water irrespective of temperature
REAL, INTENT(IN)  :: qsvp(row_length,rows,model_levels)

! Turbulent kinetic energy from BL scheme (e.g. from Section 3 Item 473 in
! UM parent model). Note: TKE is assumed to be unavailable at the top level of
! the boundary layer as is the case in the UM.
REAL, INTENT(IN)  :: bl_tke(row_length,rows,bl_levels-1)

! w component of wind aka w (e.g. from Section 0 Item 150 in UM parent model)
REAL, INTENT(IN)  :: vertvel(row_length,rows,model_levels)

! liquid cloud fraction by volume aka cf_liquid (e.g. from Section 0 Item 267
! in UM parent model)
REAL, INTENT(IN)  :: liq_cloud_frac(row_length,rows,model_levels)

! cloud liquid water aka qcl (from Section 0 Item 254 in UM parent model)
REAL, INTENT(IN)  :: cloud_liq_water(row_length,rows,model_levels)

! weighted cdnc = total cdnc * cldflg [m-3]
REAL, INTENT(OUT) :: cdncflag(row_length,rows,model_levels)

! Total number concentration of activated aerosol [m-3]
REAL, INTENT(OUT) :: n_activ_sum(row_length,rows,model_levels)

! weighted cdnc = total cdnc * liq_cloud_frac [m-3]
REAL, INTENT(OUT) :: cdncwt(row_length,rows,model_levels)

! aerosol tracers mass mixing ratio
REAL, INTENT(IN)  :: mode_tracers(row_length,rows,model_levels,n_mode_tracers)

! MODE diagnostics array
REAL, INTENT(IN OUT) :: mode_diags(row_length,rows,model_levels,n_mode_diags)

! Local variables:

! Create an array of fixed supersaturations to run the code at:
INTEGER, PARAMETER :: nsfix=19           ! number of elements in fixed-S array

INTEGER, PARAMETER :: ccnlev=1 ! sets the level on which to call fixeds

INTEGER :: imode               ! loop index to modes
INTEGER :: k                   ! loop index
INTEGER :: i                   ! loop index
INTEGER :: j                   ! loop index
INTEGER :: m                   ! loop index
INTEGER :: n                   ! loop index
INTEGER :: isfix               ! for looping
INTEGER :: icode               ! Return code  =0 Normal exit  >1 Error

LOGICAL :: first=.TRUE.

REAL, PARAMETER :: skewness=0.0   ! define fixed value of skewness of updraft
                                  ! distribution
! RW variables for ukca_abdulrazzak_ghan
REAL :: ztm1(kbdim,model_levels)        ! temperature[K]
REAL :: zapm1(kbdim,model_levels)       ! pressure [Pa=N m-2]
REAL :: zrho(kbdim,model_levels)        ! air density (kg/m3)
REAL :: zaird(kbdim, model_levels)      ! number concentration of air
                                        ! molecules (/m3)
REAL :: zqm1(kbdim,model_levels)        ! specific humidity[kg kg-1]
REAL :: zesw(kbdim,model_levels)        ! saturation water vapour pressure
REAL :: zsmax(kbdim,model_levels)       ! maximum supersaturation [fraction]
REAL :: zcdncact(kbdim,model_levels)    ! no. concentration of activated
                                        ! particles> [m-3]
REAL :: zcdncactm(kbdim,model_levels,nmodes) ! number concentration of activated
                                             ! particles by mode [m-3]
REAL :: zrdry(kbdim,model_levels,nmodes)     ! dry count median radius [m]
REAL :: zn(kbdim,model_levels,nmodes)        ! aerosol no. concentration for
                                             ! modes [m-3]
REAL :: zxtm1(kbdim,model_levels,nmodes,ncp) ! component mass mixing
                                             ! ratio [kg kg-1]

! updraft velocity array for pdf
REAL :: zsigmaw(kbdim,model_levels) ! std dev of gaussian distn of w

REAL :: zvervel_min(kbdim,model_levels)   ! min limit of updraft vel pdf [m s-1]
REAL :: zvervel_max(kbdim,model_levels)   ! max limit of updraft vel pdf [m s-1]
REAL :: zwarr(kbdim,model_levels,nwbins)  ! lin array of vert vel [m s-1]
REAL :: zwpdf(kbdim,model_levels,nwbins)  ! lin array of pdf of w
REAL :: zwbin(kbdim,model_levels,nwbins)  ! w bin width [m s-1]
REAL :: zwchar(kbdim,model_levels)        ! calculated characteristic
                                          ! updraught [m s-1]
REAL :: zwbar(kbdim,model_levels)         ! reshaped large scale vertical
                                          ! velocity [m s-1]
REAL :: zwalpha(kbdim,model_levels)       ! skewness of vertical velocity
                                          ! distribution

! RW local reshaped fields for output
REAL :: n_activated(row_length,rows,model_levels,nmodes)
                                            ! number concentration of
                                            ! activated particles by mode[m-3]
REAL :: smaxpc(row_length,rows,model_levels)      ! maximum supersaturation [%]
REAL :: wchar(row_length,rows,model_levels)   ! calculated characteristic
                                              ! updraught [m s-1]
REAL :: sigw(row_length,rows,model_levels)  ! spread of pdf of updraught [m s-1]

INTEGER :: cldflag(row_length,rows,model_levels) ! cloud flag=1 if cloud in grid
                                              ! box, else 0

REAL :: cldbase(row_length,rows,model_levels) ! cloud base=1 and no cloud below,
                                              !       else 0. Not used as a flag
                                              ! (just for diag output)

REAL :: smaxflag(row_length,rows,model_levels)  ! weighted Smax = Smax * cldflg
                                                !                           [%]
REAL :: wcharflag(row_length,rows,model_levels) ! weighted wchar = wchar*cldflg
                                                !                        [m s-1]
REAL :: sigwflag(row_length,rows,model_levels)  ! weighted sigmaw = sigw*cldflg
                                                !                        [m s-1]
REAL :: tkeflag(row_length,rows,model_levels)   ! weighted BL TKE = sigw*cldflg
                                                !                        [m s-1]
REAL :: zrho_3d(row_length,rows,model_levels)   ! air density
                                                ! (reshaped from zrho)  [kg m-3]

! array of values of fixed-s (*100 to get in %)
REAL    :: zsfix(nsfix)=[ 0.0002, 0.0004, 0.0006, 0.0008,                      &
                           0.001,  0.0016, 0.002,  0.0023,                     &
                           0.003,  0.0033, 0.0038, 0.004,                      &
                           0.005,  0.006,  0.0075, 0.008,                      &
                           0.0085, 0.01,   0.012           ]

REAL    :: zccn(kbdim,nsfix)
REAL    :: n_ccn(row_length,rows,model_levels)
REAL    :: aird_fac !Factor to convert the CDNC to per kg from per m3 for
                    !cloud base level, and then to per m3 for current level

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_ACTIVATE'

CHARACTER(LEN=errormessagelength) :: cmessage

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode = 0 ! Initialise error status

IF ( first .AND. PrintStatus >= PrStatus_Oper ) THEN
  WRITE(umMessage,'(A)') 'INPUTS to ukca_activate:'
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2I4)') 'nmodes, ncp ', nmodes, ncp
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,7L1)') 'mode ', mode
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,7I2)') 'modesol',  modesol
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Derived quantities
! ==================
DO k=1,model_levels
  ztm1(:,k) = RESHAPE(t_theta_levels(:,:,k),[kbdim])
  zapm1(:,k) = RESHAPE(p_theta_levels(:,:,k),[kbdim])
  zrho(:,k) = zapm1(:,k)/(ztm1(:,k)*rgas) ! now in kg/m3
  !(RGAS is specific gas constant for dry air in J kg^-1 K^-1
  zaird(:,k) = zapm1(:,k)/(ztm1(:,k)*boltzmann) ! number air molecules /m3
  zrho_3d(:,:,k)=RESHAPE(zrho(:,k),[row_length,rows]) !Reshape this.
  zqm1(:,k) = RESHAPE(q(:,:,k),[kbdim])
  zesw(:,k) = RESHAPE(qsvp(:,:,k),[kbdim])
END DO

IF ( glomap_config%i_ukca_activation_scheme == i_ukca_activation_arg ) THEN
  ! Set No Density (particles per m3) & MMRs from aerosol tracer array
  CALL ukca_activate_calc_zn_zxtm1 ( row_length,                               &
                                     rows,                                     &
                                     model_levels,                             &
                                     kbdim,                                    &
                                     n_mode_tracers,                           &
                                     mode_tracers,                             &
                                     zaird,                                    &
                                     zn,                                       &
                                     zxtm1 )

  ! Fill zrdry from drydiam
  CALL ukca_activate_calc_zrdry ( row_length,                                  &
                                  rows,                                        &
                                  model_levels,                                &
                                  kbdim,                                       &
                                  zrdry )

  ! Calculate updraft velocity sigw
  CALL ukca_activate_calc_updraft_velocity ( row_length,                       &
                                             rows,                             &
                                             model_levels,                     &
                                             bl_levels,                        &
                                             bl_tke,                           &
                                             sigw )

ELSE IF ( i_glomap_clim_activation_scheme == i_gc_activation_arg ) THEN
  ! Set No Density (particles per m3) & MMRs from aerosol climatolgy
  ! Fill zrdry from ukca_radaer%dry_diam
  CALL gc_activate_calc_zn_zxtm1 ( row_length,                                 &
                                   rows,                                       &
                                   model_levels,                               &
                                   kbdim,                                      &
                                   zaird,                                      &
                                   zn,                                         &
                                   zxtm1 )

  ! Fill zrdry from drydiam
  CALL ukca_activate_calc_zrdry ( row_length,                                  &
                                  rows,                                        &
                                  model_levels,                                &
                                  kbdim,                                       &
                                  zrdry )

  ! Calculate updraft velocity sigw
  CALL gc_activate_calc_updraft_velocity ( row_length,                         &
                                           rows,                               &
                                           model_levels,                       &
                                           bl_levels,                          &
                                           bl_tke,                             &
                                           sigw )

ELSE
  WRITE(cmessage,'(A,I0)') 'i_glomap_clim_activation_scheme = ',               &
                            i_glomap_clim_activation_scheme
  CALL ereport(Modulename//':'//RoutineName,icode,cmessage)

  WRITE(cmessage,'(A,I0)') 'i_ukca_activation_scheme = ',                      &
                            glomap_config%i_ukca_activation_scheme
  CALL ereport(Modulename//':'//RoutineName,icode,cmessage)

  WRITE(umMessage,'(A)') 'It should not be possible to CALL ukca_activate'     &
                          //newline// 'with these options.'
  CALL umPrint(umMessage,src=RoutineName)

  cmessage ='It should not be possible to CALL ukca_activate with these options'
  icode = 1
  CALL ereport(Modulename//':'//RoutineName,icode,cmessage)
END IF

DO k=1,model_levels
  zsigmaw(:,k)=RESHAPE(sigw(:,:,k),[kbdim])
  zwbar(:,k)=RESHAPE(vertvel(:,:,k),[kbdim])
END DO

! Define integration limits for the PDF of updraft speeds
! The original ACTIVATE code had zvervel_min=0, zvervel_max=4sigmaw.
! At convection-permitting resolution, zwbar can be greater than 4*zsigmaw
! which led to incorrect or even undefined (zero divided by zero) behaviour
! In the code below, the PDF range is still 4*sigmaw, but the range
! floats with the mean updraft if the mean updraft > 2*sigma.
IF (glomap_config%l_fix_ukca_activate_pdf) THEN
  DO k=1,model_levels
    DO j=1,kbdim
      zvervel_min(j,k)=MAX(0.0,zwbar(j,k)-2.0*zsigmaw(j,k)) ! m/s
      zvervel_max(j,k)=MAX(4.0*zsigmaw(j,k),zwbar(j,k)+2.0*zsigmaw(j,k)) !m/s
    END DO
  END DO
ELSE
  zvervel_min(:,:)=0.0
  zvervel_max(:,:)=4.0*zsigmaw(:,:)
END IF
! Set shape parameter (skewness) of skew-normal updraft distribution
zwalpha(:,:)=skewness

! vertical velocity
IF (nwbins > 1) THEN
  ! generate pdf of vertical velocity
  CALL activmklin(kbdim, model_levels, nwbins, zvervel_min, zvervel_max,       &
                  zwbin, zwarr)
  CALL activmkskew(kbdim, model_levels, nwbins, zwarr, zsigmaw(:,:),           &
                   zwbar(:,:), zwalpha(:,:), zwpdf)
ELSE IF (nwbins == 1) THEN
  ! set fixed vertical velocity
  zwbin(:,:,:)=1.0
  zwpdf(:,:,:)=1.0
  zwarr(:,:,1)=zvervel_max(:,:)
ELSE
  WRITE(umMessage,'(A,I0,A)') 'Invalid value of NWBINS: ',                     &
                   nwbins,', should be >=1'
  CALL umPrint(umMessage,src=RoutineName)
  cmessage = 'Invalid value of NWBINS, should be >=1'
  icode = 1
  CALL ereport(Modulename//':'//RoutineName,icode,cmessage)
END IF

! Initialisation:
zcdncactm(:,:,:)=0.0e0     ! initialise to 0
zcdncact(:,:)=0.0e0        ! initialise to 0
n_activ_sum(:,:,:)=0.0e0   ! initialise to 0
n_activated(:,:,:,:)=0.0e0 ! initialise to 0
n_ccn(:,:,:)=0.0e0         ! initialise to 0
zccn(:,:)=0.0e0            ! initialise to 0
zsmax(:,:)=0.0e0           ! initialise to 0
smaxpc(:,:,:)=0.0e0        ! initialise to 0
zwchar(:,:)=0.0e0          ! initialise to 0
wchar(:,:,:)=0.0e0         ! initialise to 0
cdncwt(:,:,:)=0.0e0        ! initialise to 0
cldflag(:,:,:)=0           ! initialise to 0
cldbase(:,:,:)=0.0e0       ! initialise to 0
cdncflag(:,:,:)=0.0e0      ! initialise to 0
smaxflag(:,:,:)=0.0e0      ! initialise to 0
wcharflag(:,:,:)=0.0e0     ! initialise to 0
sigwflag(:,:,:)=0.0e0      ! initialise to 0
tkeflag(:,:,:)=0.0e0       ! initialise to 0

IF ( PrintStatus == PrStatus_Diag ) THEN
  WRITE(umMessage,'(A)') 'UKCA_ACTIVATE: Summary of data: '
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'bl_tke',MINVAL(bl_tke(:,:,:)),           &
              MAXVAL(bl_tke(:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'vertvel',MINVAL(vertvel(:,:,:)),         &
              MAXVAL(vertvel(:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'liq_cloud_frac',                         &
         MINVAL(liq_cloud_frac(:,:,:)), MAXVAL(liq_cloud_frac(:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'cloud_liq_water',                        &
        MINVAL(cloud_liq_water(:,:,:)), MAXVAL(cloud_liq_water(:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zcdncactm',                              &
       MINVAL(zcdncactm(:,:,:)), MAXVAL(zcdncactm(:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zcdncact',MINVAL(zcdncact(:,:)),         &
                                      MAXVAL(zcdncact(:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zwchar',MINVAL(zwchar(:,:)),             &
                                      MAXVAL(zwchar(:,:))
  CALL umPrint(umMessage,src=RoutineName)
  ! write out min and max values of aerosol number conc. and Drydiam
  ! for each mode
  DO imode=1,nmodes
    IF (mode(imode)) THEN
      WRITE(umMessage,'(A,I4,2(E12.4,1x))') 'zn, imode=',imode,                &
                           MINVAL(zn(:,:,imode)), MAXVAL(zn(:,:,imode))
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,I4,2(E12.4,1x))') 'zrdry, imode=', imode,            &
                      MINVAL(zrdry(:,:,imode)), MAXVAL(zrdry(:,:,imode))
      CALL umPrint(umMessage,src=RoutineName)
    END IF
  END DO

  WRITE(umMessage,'(A,2(E12.4,1x))') 'ztm1', MINVAL(ztm1(:,:)),                &
                                             MAXVAL(ztm1(:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))')'zapm1', MINVAL(zapm1(:,:)),               &
                                             MAXVAL(zapm1(:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zqm1', MINVAL(zqm1(:,:)),                &
                                             MAXVAL(zqm1(:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zvervel_max',                            &
               MINVAL(zvervel_max(:,:)), MAXVAL(zvervel_max(:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zwpdf', MINVAL(zwpdf(:,:,:)),            &
                                              MAXVAL(zwpdf(:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zwarr', MINVAL(zwarr(:,:,:)),            &
                                              MAXVAL(zwarr(:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zwbin', MINVAL(zwbin(:,:,:)),            &
                                              MAXVAL(zwbin(:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
END IF     ! Prstatus

CALL ukca_abdulrazzak_ghan(kbdim,   model_levels,                              &
                           zesw,                                               &
                           zn,      zxtm1,                                     &
                           ztm1,    zapm1,                                     &
                           zqm1,    zrdry,                                     &
                           nwbins,  zwarr,                                     &
                           zwpdf,   zwbin,                                     &
                           zsmax,   zwchar,                                    &
                           zcdncactm, zcdncact )

IF ( PrintStatus == PrStatus_Diag ) THEN
  WRITE(umMessage,'(A)') 'UKCA_ACTIVATE: Summary after ARG'
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zcdncactm', MINVAL(zcdncactm(:,:,:)),    &
                                                  MAXVAL(zcdncactm(:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zcdncact',  MINVAL(zcdncact(:,:)),       &
                                                  MAXVAL(zcdncact(:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zwchar',    MINVAL(zwchar(:,:)),         &
                                                  MAXVAL(zwchar(:,:))
  CALL umPrint(umMessage,src=RoutineName)
END IF

! reshape zcdncactm to n_activated
DO imode=1, nmodes
  IF (mode(imode)) THEN
    DO k=1,model_levels
      n_activated(:,:,k,imode)=RESHAPE(zcdncactm(:,k,imode),                   &
                                [row_length,rows])
    END DO
  END IF
END DO

!reshape zcdncact to n_activ_sum
DO k=1,model_levels
  n_activ_sum(:,:,k)=RESHAPE(zcdncact(:,k),[row_length,rows])
  !reshape zsmax to smaxpc and convert to %
  smaxpc(:,:,k)=100.0*(RESHAPE(zsmax(:,k),[row_length,rows]))
  !reshape zwchar to wchar
  wchar(:,:,k)=RESHAPE(zwchar(:,k),[row_length,rows])
END DO

IF ( PrintStatus == PrStatus_Diag ) THEN
  WRITE(umMessage,'(A,2(E12.4,1x))') ' n_activated',                           &
         MINVAL(n_activated(:,:,:,:)),MAXVAL(n_activated(:,:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'n_activ_sum',                            &
                 MINVAL(n_activ_sum(:,:,:)), MAXVAL(n_activ_sum(:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zsmax', MINVAL(zsmax(:,:)),              &
                                              MAXVAL(zsmax(:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'smaxpc', MINVAL(smaxpc(:,:,:)),          &
                                               MAXVAL(smaxpc(:,:,:))
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,2(E12.4,1x))') 'zsigmaw', MINVAL(zsigmaw(:,:)),          &
                                                MAXVAL(zsigmaw(:,:))
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Set CDNC minimum to the equivalent of 5 cm-3 for consistency with the
! minimum values assumed in the radiation and precipitation schemes.
WHERE (n_activ_sum < 5.0e6) n_activ_sum = 5.0e6

! For bottom level set CDNC if liquid cloud fraction > 0
! and LWC > 0 kg kg^-1
DO j=1, rows
  DO i=1, row_length
    IF (liq_cloud_frac(i,j,1) > 0.0 .AND.                                      &
        cloud_liq_water(i,j,1) > 0.0 ) THEN
      cdncwt(i,j,1)=liq_cloud_frac(i,j,1)*n_activ_sum(i,j,1)
      cldflag(i,j,1)=1
      cldbase(i,j,1)=1.0
      cdncflag(i,j,1)=n_activ_sum(i,j,1)
      smaxflag(i,j,1)=smaxpc(i,j,1)
      wcharflag(i,j,1)=wchar(i,j,1)
      sigwflag(i,j,1)=sigw(i,j,1)
      tkeflag(i,j,1)=bl_tke(i,j,1)
    END IF
  END DO
END DO

! Calculate weighted CDNC (total CDNC weighted by liquid cloud fraction)
! make a cloud flag diagnostic to tell us when and where there is > 0% cloud
! fraction. Also make an unweighted CDNC diag, zero outside cloud.
DO k=2, model_levels
  DO j=1, rows
    DO i=1, row_length
      ! Define cloudy if liquid cloud fraction > 0 and LWC > 0 kg kg^-1
      IF (liq_cloud_frac(i,j,k) > 0.0 .AND.                                    &
          cloud_liq_water(i,j,k) > 0.0 ) THEN
        cldflag(i,j,k)=1
        ! If layer below (i.e. k-1) is cloudy then set CDNC in this level (k)
        ! to number activated at cloud base
        IF ( cldflag(i,j,k-1) == 1 ) THEN
          ! Note that the bug fix below is automatically applied to the
          ! GLOMAP climatology aerosol scheme whereas it controlled by the
          ! temporary logical for interactive GLOMAP aerosols.
          IF ( glomap_config%l_fix_ukca_activate_vert_rep .OR.                 &
               ( i_glomap_clim_activation_scheme == i_gc_activation_arg ) ) THEN
            !cdncflag is in per m3, but actually want to vertically replicate
            !the per kg value. So, divide by air density at level below to
            !convert to per kg and then use this for the current level.
            !Multiply by air density at current level to get per m3 again for
            !output from this subroutine.
            aird_fac = zrho_3d(i,j,k) / zrho_3d(i,j,k-1)
            cdncwt(i,j,k)=liq_cloud_frac(i,j,k)*cdncflag(i,j,k-1) * aird_fac
            cdncflag(i,j,k)=cdncflag(i,j,k-1) * aird_fac
          ELSE
            cdncwt(i,j,k)=liq_cloud_frac(i,j,k)*cdncflag(i,j,k-1)
            cdncflag(i,j,k)=cdncflag(i,j,k-1)
          END IF
          smaxflag(i,j,k)=smaxflag(i,j,k-1)
          wcharflag(i,j,k)=wcharflag(i,j,k-1)
          sigwflag(i,j,k)=sigwflag(i,j,k-1)
          ! If layer below (k-1) is cloudy, then this layer(k) is not the
          ! cloud base, so
          cldbase(i,j,k)=0.0
          IF ( k <= bl_levels-1 ) tkeflag(i,j,k)=tkeflag(i,j,k-1)
        ELSE
          cdncwt(i,j,k)=liq_cloud_frac(i,j,k)*n_activ_sum(i,j,k)
          cdncflag(i,j,k)=n_activ_sum(i,j,k)
          smaxflag(i,j,k)=smaxpc(i,j,k)
          wcharflag(i,j,k)=wchar(i,j,k)
          sigwflag(i,j,k)=sigw(i,j,k)
          ! But if layer below(k-1) is not cloudy, then this layer(k) is the
          ! cloud base
          cldbase(i,j,k)=1.0
          IF ( k <= bl_levels-1 ) tkeflag(i,j,k)=bl_tke(i,j,k)
        END IF
      END IF    ! liq_cloud_frac > 0 etc
    END DO ! i
  END DO ! j
END DO ! k

IF ( i_glomap_clim_activation_scheme == i_gc_activation_arg ) THEN
  ! Return cdncflag and cdncwt
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Calculate CCN at fixed supersaturation:
IF (glomap_config%l_ukca_sfix) THEN

  CALL ukca_fixeds(kbdim,                                                      &
                   zccn,    zn(:,ccnlev,:),                                    &
                   zxtm1(:,ccnlev,:,:),                                        &
                   ztm1(:,ccnlev),                                             &
                   zrdry(:,ccnlev,:),                                          &
                   zsfix,   nsfix )

  !reshape zccn to n_ccn:
  DO isfix=1, nsfix
    n_ccn(:,:,isfix)=RESHAPE(zccn(:,isfix),                                    &
         [row_length,rows])
  END DO

END IF !l_ukca_sfix


!  Write 3_D diagnostics to mode_diags array
!  N.B. l_ukca_mode_diags is set whenever a STASH request for a relevant
!  item is found, and the ukcaD1codes(N)%item is then set, otherwise it is IMDI
IF (l_ukca_mode_diags) THEN ! fill 3D array
  m=0
  DO n=1,Nukca_D1items
    IF (ukcaD1codes(n)%section == stashcode_glomap_sec .AND.                   &
        ukcaD1codes(n)%item >= item1_mode_diags .AND.                          &
        ukcaD1codes(n)%item <= item1_mode_diags+                               &
                            nmax_mode_diags-1 .AND.                            &
        ukcaD1codes(n)%item /= imdi) THEN
      m=m+1

      SELECT CASE(UkcaD1codes(n)%item)
      CASE (item1_mode_diags+268) ! CDNC - NUC(sol)
        mode_diags(:,:,:,m)=n_activated(:,:,:,1)
      CASE (item1_mode_diags+269) ! CDNC - AIT(sol)
        mode_diags(:,:,:,m)=n_activated(:,:,:,2)
      CASE (item1_mode_diags+270)     ! CDNC - ACC(sol)
        mode_diags(:,:,:,m)=n_activated(:,:,:,3)
      CASE (item1_mode_diags+271)     ! CDNC - COR(sol)
        mode_diags(:,:,:,m)=n_activated(:,:,:,4)
      CASE (item1_mode_diags+272)     ! Smax %
        mode_diags(:,:,:,m)=smaxpc(:,:,:)
      CASE (item1_mode_diags+273)     ! cloud base
        mode_diags(:,:,:,m)=cldbase(:,:,:)
      CASE (item1_mode_diags+274)     ! sig_w [m/s]
        mode_diags(:,:,:,m)=sigw(:,:,:)
      CASE (item1_mode_diags+275)     ! liq cloud frac
        mode_diags(:,:,:,m)=liq_cloud_frac(:,:,:)
      CASE (item1_mode_diags+276)     ! total CDNC * liq cloud frac [m-3]
        mode_diags(:,:,:,m)=cdncwt(:,:,:)
      CASE (item1_mode_diags+277)     ! cloud flag
        mode_diags(:,:,:,m)=REAL(cldflag(:,:,:))
                       !Convert to REAL for averaging and diag output
      CASE (item1_mode_diags+278)     ! total CDNC * cloud flag [m-3]
        mode_diags(:,:,:,m)=cdncflag(:,:,:)
      CASE (item1_mode_diags+279)     ! Smaxpc * cloud flag [%]
        mode_diags(:,:,:,m)=smaxflag(:,:,:)
      CASE (item1_mode_diags+280)     ! Wchar * cloud flag [m s-1]
        mode_diags(:,:,:,m)=wcharflag(:,:,:)
      CASE (item1_mode_diags+281)     ! Sigw * cloud flag [m s-1]
        mode_diags(:,:,:,m)=sigwflag(:,:,:)
      CASE (item1_mode_diags+282)     ! BL TKE * cloud flag [m2 s-2]
        mode_diags(:,:,:,m)=tkeflag(:,:,:)
      CASE (item1_mode_diags+283)     ! CCN at fixed S [m-3]
        IF (glomap_config%l_ukca_sfix) THEN
          mode_diags(:,:,:,m)= n_ccn(:,:,:)
        ELSE
          mode_diags(:,:,:,m)=0.0
        END IF !l_ukca_sfix
      END SELECT

    END IF

  END DO
END IF       ! l_ukca_mode_diags

first=.FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE ukca_activate

SUBROUTINE ukca_activate_calc_zn_zxtm1 ( row_length,                           &
                                         rows,                                 &
                                         model_levels,                         &
                                         kbdim,                                &
                                         n_mode_tracers,                       &
                                         mode_tracers,                         &
                                         zaird,                                &
                                         zn,                                   &
                                         zxtm1 )

USE ereport_mod,               ONLY:                                           &
    ereport

USE errormessagelength_mod,    ONLY:                                           &
    errormessagelength

USE parkind1,                  ONLY:                                           &
    jprb,                                                                      &
    jpim

USE ukca_mode_setup,           ONLY:                                           &
    nmodes,                                                                    &
    mode,                                                                      &
    ncp,                                                                       &
    component

USE ukca_mode_tracer_maps_mod, ONLY:                                           &
    mmr_index,                                                                 &
    nmr_index

USE asad_mod, ONLY: jpctr

USE umPrintMgr,                ONLY:                                           &
    umPrint,                                                                   &
    umMessage

USE yomhook,                   ONLY:                                           &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: row_length     ! array dimension
INTEGER, INTENT(IN) :: rows           ! array dimension
INTEGER, INTENT(IN) :: model_levels   ! array dimension
INTEGER, INTENT(IN) :: kbdim          ! = theta_field_size=row_length*rows
INTEGER, INTENT(IN) :: n_mode_tracers ! # of mode tracers

! aerosol tracers mass mixing ratio
REAL, INTENT(IN)   :: mode_tracers(row_length,rows,model_levels,n_mode_tracers)

! number concentration of air molecules (/m3)
REAL, INTENT(IN)  :: zaird(kbdim,model_levels)

! aerosol no. concentration for modes [m-3]
REAL, INTENT(OUT) :: zn(kbdim,model_levels,nmodes)

! component mass mixing ratio [kg kg-1]
REAL, INTENT(OUT) :: zxtm1(kbdim,model_levels,nmodes,ncp)

! Local variables:
INTEGER :: icode               ! Return code  =0 Normal exit  >1 Error
INTEGER :: imode               ! loop index to modes
INTEGER :: icp                 ! loop index to components
INTEGER :: itra                ! tracer index
INTEGER :: ifirst              ! First tracer index in nmr_index, mmr_index
INTEGER :: k                   ! loop index

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='UKCA_ACTIVATE_CALC_ZN_ZXTM1'

INTEGER(KIND=jpim), PARAMETER     :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER     :: zhook_out = 1
REAL(KIND=jprb)                   :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! find the first tracer in nmr_index and mmr_index (these refer to all tracers)
ifirst = jpctr + 1

! initialise to zero
zn     = 0.0
zxtm1  = 0.0

! Aerosol Tracers
! ===============
! Loop through the modes
DO imode=1,nmodes
  IF (mode(imode)) THEN
    itra = nmr_index(imode) - ifirst + 1
    ! Set No Density (particles per m3) from aerosol tracer array
    != number particles per molecule air * number concentration of air molecules
    DO k=1,model_levels
      zn(:,k,imode)=RESHAPE(mode_tracers(:,:,k,itra),[kbdim]) * zaird(:,k)
    END DO

    ! Loop through the components
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        itra = mmr_index(imode,icp) - ifirst + 1
        IF (itra <= 0) THEN
          WRITE(umMessage,'(A)') '***** ERROR: MD_CP specified by COMPONENT'
          CALL umPrint(umMessage,src=RoutineName)

          WRITE(umMessage,'(A)') '*****        doesnt exist in MMR_INDEX'
          CALL umPrint(umMessage,src=RoutineName)

          WRITE(umMessage,'(A,2I5)') 'IMODE,ICP = ', imode, icp
          CALL umPrint(umMessage,src=RoutineName)

          WRITE(umMessage,'(A,7L2)') 'COMPONENT(IMODE,ICP)= ',                 &
                                      component(imode,icp)
          CALL umPrint(umMessage,src=RoutineName)

          WRITE(umMessage,'(A,I5)') 'mmr_index = ', mmr_index(imode, icp)
          CALL umPrint(umMessage,src=RoutineName)

          cmessage=' MD_CP specified by COMPONENT not in mmr_index'
          icode = 1
          CALL ereport(Modulename//':'//RoutineName,icode,cmessage)
        END IF
        ! Set mass mixing ratios from aerosol tracer array
        DO k=1,model_levels
          zxtm1(:,k,imode,icp) = RESHAPE(mode_tracers(:,:,k,itra),[kbdim])
        END DO
      END IF
    END DO
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE ukca_activate_calc_zn_zxtm1

SUBROUTINE gc_activate_calc_zn_zxtm1 ( row_length,                             &
                                       rows,                                   &
                                       model_levels,                           &
                                       kbdim,                                  &
                                       zaird,                                  &
                                       zn,                                     &
                                       zxtm1 )

USE ukca_um_legacy_mod,     ONLY:                                              &
    gc_nd_ait_sol,                                                             &
    gc_ait_sol_su,                                                             &
    gc_ait_sol_bc,                                                             &
    gc_ait_sol_oc,                                                             &
    gc_nd_acc_sol,                                                             &
    gc_acc_sol_su,                                                             &
    gc_acc_sol_bc,                                                             &
    gc_acc_sol_oc,                                                             &
    gc_acc_sol_ss,                                                             &
    gc_nd_cor_sol,                                                             &
    gc_cor_sol_su,                                                             &
    gc_cor_sol_bc,                                                             &
    gc_cor_sol_oc,                                                             &
    gc_cor_sol_ss,                                                             &
    gc_nd_ait_ins,                                                             &
    gc_ait_ins_bc,                                                             &
    gc_ait_ins_oc,                                                             &
    i_glomap_clim_setup,                                                       &
    i_gc_sussocbc_5mode

USE ereport_mod,            ONLY:                                              &
    ereport

USE errormessagelength_mod, ONLY:                                              &
    errormessagelength

USE parkind1,               ONLY:                                              &
    jprb,                                                                      &
    jpim

USE ukca_mode_setup,        ONLY:                                              &
    nmodes,                                                                    &
    mode,                                                                      &
    ncp,                                                                       &
    component,                                                                 &
    mode_ait_sol,                                                              &
    mode_acc_sol,                                                              &
    mode_cor_sol,                                                              &
    mode_ait_insol,                                                            &
    cp_su,                                                                     &
    cp_bc,                                                                     &
    cp_oc,                                                                     &
    cp_cl

USE umPrintMgr,             ONLY:                                              &
    newline

USE yomhook,                ONLY:                                              &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: row_length     ! array dimension
INTEGER, INTENT(IN) :: rows           ! array dimension
INTEGER, INTENT(IN) :: model_levels   ! array dimension
INTEGER, INTENT(IN) :: kbdim          ! = theta_field_size=row_length*rows

! number concentration of air molecules (/m3)
REAL, INTENT(IN)  :: zaird ( kbdim, model_levels )

! aerosol no. concentration for modes [m-3]
REAL, INTENT(OUT) :: zn ( kbdim, model_levels, nmodes )

! component mass mixing ratio [kg kg-1]
REAL, INTENT(OUT) :: zxtm1 ( kbdim, model_levels, nmodes, ncp )

! Local variables:

INTEGER :: icode   ! error code
INTEGER :: imode   ! mode counter
INTEGER :: icp     ! component counter
INTEGER :: k,j,i,l ! model_levels counter

REAL, ALLOCATABLE :: znmr3d  (:,:,:,:)
REAL, ALLOCATABLE :: zammr3d (:,:,:,:,:)

CHARACTER(LEN=errormessagelength) :: cmessage ! error message

CHARACTER(LEN=*), PARAMETER   :: RoutineName='GC_ACTIVATE_CALC_ZN_ZXTM1'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE(i_glomap_clim_setup)
CASE (i_gc_sussocbc_5mode)
  ! There are 4 modes
  ALLOCATE(znmr3d( 1:row_length, 1:rows, 1:model_levels, nmodes ) )
  znmr3d = 0.0
  znmr3d(:,:,1:model_levels,mode_ait_sol)   = gc_nd_ait_sol(:,:,1:model_levels)
  znmr3d(:,:,1:model_levels,mode_acc_sol)   = gc_nd_acc_sol(:,:,1:model_levels)
  znmr3d(:,:,1:model_levels,mode_cor_sol)   = gc_nd_cor_sol(:,:,1:model_levels)
  znmr3d(:,:,1:model_levels,mode_ait_insol) = gc_nd_ait_ins(:,:,1:model_levels)


  ! There are 13 components
  ALLOCATE(zammr3d( 1:row_length, 1:rows, 1:model_levels, nmodes, ncp ) )
  zammr3d = 0.0
  zammr3d(:,:,1:model_levels,mode_ait_sol,cp_su)   =                           &
                                              gc_ait_sol_su(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_ait_sol,cp_bc)   =                           &
                                              gc_ait_sol_bc(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_ait_sol,cp_oc)   =                           &
                                              gc_ait_sol_oc(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_acc_sol,cp_su)   =                           &
                                              gc_acc_sol_su(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_acc_sol,cp_bc)   =                           &
                                              gc_acc_sol_bc(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_acc_sol,cp_oc)   =                           &
                                              gc_acc_sol_oc(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_acc_sol,cp_cl)   =                           &
                                              gc_acc_sol_ss(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_cor_sol,cp_su)   =                           &
                                              gc_cor_sol_su(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_cor_sol,cp_bc)   =                           &
                                              gc_cor_sol_bc(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_cor_sol,cp_oc)   =                           &
                                              gc_cor_sol_oc(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_cor_sol,cp_cl)   =                           &
                                              gc_cor_sol_ss(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_ait_insol,cp_bc) =                           &
                                              gc_ait_ins_bc(:,:,1:model_levels)
  zammr3d(:,:,1:model_levels,mode_ait_insol,cp_oc) =                           &
                                              gc_ait_ins_oc(:,:,1:model_levels)

CASE DEFAULT
  WRITE(cmessage,'(A,I0)') 'Unrecognised glomap_clim setup' //newline//        &
                           'i_glomap_clim_setup = ', i_glomap_clim_setup
  icode = 1
  CALL ereport(Modulename//':'//RoutineName,icode,cmessage)
END SELECT

zn    = 0.0
zxtm1 = 0.0

! Loop through the modes
DO imode=1,nmodes
  IF (mode(imode)) THEN

    DO k=1,model_levels
      l = 0
      DO j = 1, rows
        DO i = 1, row_length

          !Calculate vector position
          l = l + 1

          !Set No Density (particles per m3)
          !number particles per molecule air * no concentration of air molecules
          zn(l,k,imode) = znmr3d(i,j,k,imode) * zaird(l,k)

        END DO
      END DO
    END DO

    ! Loop through the components
    DO icp=1,ncp
      IF (component(imode,icp)) THEN

        DO k=1,model_levels
          l = 0
          DO j = 1, rows
            DO i = 1, row_length

              !Calculate vector position
              l = l + 1

              ! Set mass mixing ratios from aerosol climatology fields
              zxtm1(l,k,imode,icp) = zammr3d(i,j,k,imode,icp)

            END DO
          END DO
        END DO

      END IF
    END DO

  END IF
END DO

DEALLOCATE(zammr3d)
DEALLOCATE(znmr3d)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE gc_activate_calc_zn_zxtm1

SUBROUTINE ukca_activate_calc_zrdry (row_length,rows,model_levels,kbdim,zrdry)

USE ereport_mod,            ONLY:                                              &
    ereport

USE errormessagelength_mod, ONLY:                                              &
    errormessagelength

USE parkind1,               ONLY:                                              &
    jprb,                                                                      &
    jpim

USE ukca_drydiam_field_mod, ONLY:                                              &
    drydiam

USE ukca_mode_setup,        ONLY:                                              &
    nmodes,                                                                    &
    mode

USE missing_data_mod,       ONLY:                                              &
    rmdi

USE umPrintMgr,             ONLY:                                              &
    newline

USE yomhook,                ONLY:                                              &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! = theta_field_size=row_length*rows
INTEGER, INTENT(IN) :: kbdim

! dry count median radius [m]
REAL, INTENT(OUT)   :: zrdry ( kbdim, model_levels, nmodes )

! Local variables:

INTEGER :: imode               ! loop index to modes
INTEGER :: k,j,i,l             ! loop index
INTEGER :: icode               ! Return code  =0 Normal exit  >1 Error

CHARACTER(LEN=errormessagelength) :: cmessage ! error message

INTEGER(KIND=jpim), PARAMETER     :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER     :: zhook_out = 1
REAL(KIND=jprb)                   :: zhook_handle

CHARACTER(LEN=*), PARAMETER       :: RoutineName='UKCA_ACTIVATE_CALC_ZRDRY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise zrdry
zrdry(:,:,:) = rmdi

! Fill zrdry from drydiam
IF (ALLOCATED(drydiam)) THEN

  DO imode=1,nmodes
    IF (mode(imode)) THEN

      DO k=1,model_levels
        l = 0
        DO j = 1, rows
          DO i = 1, row_length

            !Calculate vector position
            l = l + 1

            ! Convert values from diameter to radius as expected by ARG method
            zrdry(l,k,imode) = drydiam(i,j,k,imode) * 0.5

          END DO
        END DO
      END DO

    END IF
  END DO

  DEALLOCATE(drydiam)


ELSE
  cmessage = 'Aerosol Dry diameter not available from GLOMAP.' //newline//     &
             'Array drydiam is not allocated (in Aero_Ctl)'
  icode = 1
  CALL ereport(RoutineName,icode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE ukca_activate_calc_zrdry

SUBROUTINE ukca_activate_calc_updraft_velocity ( row_length,                   &
                                                 rows,                         &
                                                 model_levels,                 &
                                                 bl_levels,                    &
                                                 bl_tke,                       &
                                                 sigw )

USE ukca_config_specification_mod, ONLY: glomap_config

USE parkind1,             ONLY:                                                &
    jprb,                                                                      &
    jpim

USE yomhook,              ONLY:                                                &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: row_length        ! array dimension
INTEGER, INTENT(IN) :: rows              ! array dimension
INTEGER, INTENT(IN) :: model_levels      ! array dimension
INTEGER, INTENT(IN) :: bl_levels         ! array dimension

! Turbulent kinetic energy from BL scheme
! Note: TKE is assumed to be unavailable at the top level of the boundary layer
REAL, INTENT(IN)    :: bl_tke ( row_length, rows, bl_levels-1 )

! spread of pdf of updraught [m s-1]
REAL, INTENT(OUT)   :: sigw ( row_length, rows, model_levels )

! Local variables:

INTEGER :: k                             ! loop index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_ACTIVATE_CALC_UPDRAFT_VELOCITY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Updraught velocity pdf std. dev. (m/s)
! ======================================

! Define width of updraught velocity pdf based on BL TKE
! (Factor of 2/3 owing to the assumption that TKE is isotropic)
IF (glomap_config%l_bug_repro_tke_index) THEN
  ! Reproduce effects of an old bug related to TKE indexing in the UM
  sigw(:,:,1)=glomap_config%sigwmin
  DO k=2, bl_levels
    sigw(:,:,k)=MAX(glomap_config%sigwmin,(two_thirds*bl_tke(:,:,k-1))**0.5)
  END DO
  DO k = bl_levels + 1, model_levels
    sigw(:,:,k)=glomap_config%sigwmin
  END DO
ELSE
  DO k=1, bl_levels - 1
    sigw(:,:,k)=MAX(glomap_config%sigwmin,(two_thirds*bl_tke(:,:,k))**0.5)
  END DO
  DO k = bl_levels, model_levels
    sigw(:,:,k)=glomap_config%sigwmin
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE ukca_activate_calc_updraft_velocity

SUBROUTINE gc_activate_calc_updraft_velocity ( row_length,                     &
                                               rows,                           &
                                               model_levels,                   &
                                               bl_levels,                      &
                                               bl_tke,                         &
                                               sigw )

USE parkind1,             ONLY:                                                &
    jprb,                                                                      &
    jpim

USE yomhook,              ONLY:                                                &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: row_length        ! array dimension
INTEGER, INTENT(IN) :: rows              ! array dimension
INTEGER, INTENT(IN) :: model_levels      ! array dimension
INTEGER, INTENT(IN) :: bl_levels         ! array dimension

! Turbulent kinetic energy from BL scheme
! Note: TKE is assumed to be unavailable at the top level of the boundary layer
REAL, INTENT(IN)    :: bl_tke ( row_length, rows, bl_levels - 1 )

! spread of pdf of updraught [m s-1]
REAL, INTENT(OUT)   :: sigw ( row_length, rows, model_levels )

! Local variables:

INTEGER :: i,j,k                   ! loop index

REAL, PARAMETER :: sigwmin = 0.01  ! minimum value of sigmaw

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GC_ACTIVATE_CALC_UPDRAFT_VELOCITY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

sigw(:,:,:)=0.0e0         ! initialise to 0

! Define width of updraught velocity pdf based on BL TKE
! (Factor of 2/3 owing to the assumption that TKE is isotropic)

DO k = 1, bl_levels-1
  DO j = 1, rows
    DO i = 1, row_length
      ! Fix zsigmaw (std dev of updraught vel pdf in m/s) :
      sigw(i,j,k) = MAX( sigwmin, (two_thirds*bl_tke(i,j,k))**0.5 )
    END DO
  END DO
END DO

DO k = bl_levels, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      sigw(i,j,k) = sigwmin
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE gc_activate_calc_updraft_velocity

END MODULE ukca_activate_mod

