! *****************************COPYRIGHT*******************************
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
! Description:
!    Module used to calculate first order reaction coefficients for
!    the loss of gas phase species (N2O5 and HO2) to aerosol surfaces.
!
! Method:
!    The module contains the following functions and subroutines for
!    the calculation of the first order reaction rates:
!
!    * ukca_trop_hetchem: Subroutine including some calls to cal_n2o5,
!        cal_ho2 and cal_loss1k in order to calculate the first order
!        rate coefficients for each aerosol mode and add them up to
!        get a final rate.
!
!    * cal_n2o5:   Function to compute the GAMMA sticking factor for
!                  the heterogeneous hydrolysis of N2O5.
!
!    * cal_ho2:    Function to compute the gamma sticking factor for
!                  the HO2 self reaction.
!
!    * cal_loss1k: Function to calculate the first order loss rate of
!                  species on wet aerosol surface.
!
! UKCA is a community model supported by The Met Office and
! NCAS, with components initially provided by the University of
! Cambridge, University of Leeds and the Met Office. See
! www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 2003.
!   This code is written to UMDP3 programming standards.
!
! ######################################################################
!
MODULE ukca_trop_hetchem_mod

USE ukca_config_constants_mod, ONLY: rmol, avogadro
USE ukca_constants,       ONLY: pi, m_ho2, m_n2o5, m_air

USE ukca_mode_setup,      ONLY: nmodes,                                        &
                                cp_su, cp_bc, cp_oc, cp_cl,                    &
                                cp_du, cp_so, cp_no3, cp_nh4, cp_nn

USE ukca_config_specification_mod, ONLY: ukca_config, glomap_variables
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
USE ereport_mod,          ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Diag
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! All variables and functions private by default
PRIVATE
PUBLIC :: ukca_trop_hetchem, cal_n2o5, cal_loss1k

! Indices for the location of each rate in the returned array.
! 1. Index for heterogeneous hydrolysis of N2O5
INTEGER, PARAMETER, PUBLIC :: ihet_n2o5    = 1
! 2. Index for self reaction of HO2 on surfaces
INTEGER, PARAMETER, PUBLIC :: ihet_ho2_ho2 = 2

INTEGER    :: errcode            ! Error code

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_TROP_HETCHEM_MOD'

CONTAINS

SUBROUTINE ukca_trop_hetchem(nbox, nhet, temp, rh, aird, pvol,                 &
                             wetdp, sarea, het_rates)

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN) :: nbox                 ! No of points
INTEGER, INTENT(IN) :: nhet                 ! No of heterogeneous reaction rates

REAL, INTENT(IN)    :: temp(nbox)           ! temperature [K]
REAL, INTENT(IN)    :: rh(nbox)             ! Relative humidity [fraction]
! Total number density [molecules/cm^3]
REAL, INTENT(IN)    :: aird(nbox)
REAL, INTENT(IN)    :: pvol(nbox,nmodes,glomap_variables%ncp)
! .. Partial volume of soluble components as fraction of
! .. solution volume according to fraction of total mass (including water)
REAL, INTENT(IN)    :: wetdp(nbox,nmodes)   ! Mean wet radius for the mode [m]
REAL, INTENT(IN)    :: sarea(nbox,nmodes)   ! Surface area concentn. [cm^2/cm^3]
REAL, INTENT(OUT)   :: het_rates(nbox,nhet) ! rate coefficients

! Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
LOGICAL, POINTER :: mode (:)
INTEGER, POINTER :: ncp

INTEGER :: i                            ! Loop variable
INTEGER :: imode                        ! mode loop counter
INTEGER :: icp                          ! component loop counter

REAL    :: molec_weight(nhet)           ! Molecular weights of species

! Gamma calculated by CAL_* [unitless]
REAL    :: sticking_cf(nbox,nhet,glomap_variables%ncp)

REAL    :: total_vol(nbox,nmodes)       ! Total volume of each mode

! pvol/total_vol for each mode, unitless
REAL    :: frac_vol(nbox,nmodes,glomap_variables%ncp)

REAL    :: mode_gamma(nbox,nhet,nmodes) ! Weighted gammas for each mode
REAL    :: mode_rate(nbox,nhet,nmodes)  ! First order loss rate coefficient
REAL    :: total_sa(nbox)               ! Total aerosol surface area

CHARACTER(LEN=errormessagelength) :: cmessage           ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_TROP_HETCHEM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp

IF (PrintStatus >= PrStatus_Diag) THEN
  WRITE(umMessage,'(A)') 'UKCA_TROP_HETCHEM Routine:'
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
  WRITE(umMessage,'(A,I5)') 'nbox: ',nbox
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
  WRITE(umMessage,'(A,3E15.6)')                                                &
             'temp: ',MAXVAL(temp),MINVAL(temp),SUM(temp)/SIZE(temp)
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
  WRITE(umMessage,'(A,3E15.6)')                                                &
             'rh  : ',MAXVAL(rh),MINVAL(rh),SUM(rh)/SIZE(rh)
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
  WRITE(umMessage,'(A,3E15.6)')                                                &
             'aird: ',MAXVAL(aird),MINVAL(aird),SUM(aird)/SIZE(aird)
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
  DO imode=1,nmodes
    IF (mode(imode)) THEN
      WRITE(umMessage,'(A,I0,3E15.6)')                                         &
                  'sarea: ',imode,MAXVAL(sarea(:,imode)),                      &
                  MINVAL(sarea(:,imode)),                                      &
                  SUM(sarea(:,imode))/SIZE(sarea(:,1))
      CALL umPrint(umMessage,src='ukca_trop_hetchem')
      WRITE(umMessage,'(A,I0,3E15.6)')                                         &
                  'wetdp: ',imode,MAXVAL(wetdp(:,imode)),                      &
                  MINVAL(wetdp(:,imode)),                                      &
                  SUM(wetdp(:,imode))/SIZE(wetdp(:,1))
      CALL umPrint(umMessage,src='ukca_trop_hetchem')
    END IF
  END DO

  DO imode=1,nmodes
    IF (mode(imode)) THEN
      DO icp=1,ncp
        IF ( (component(imode,icp)) .AND. (icp /= cp_nh4) ) THEN
          WRITE(umMessage,'(A,2I0,3E15.6)')                                    &
                      'pvol: ',imode,icp,MAXVAL(pvol(:,imode,icp)),            &
                      MINVAL(pvol(:,imode,icp)),                               &
                      SUM(pvol(:,imode,icp))/SIZE(pvol(:,1,1))
          CALL umPrint(umMessage,src='ukca_trop_hetchem')
        END IF
      END DO
    END IF
  END DO
END IF

molec_weight(:) = 9.0e5
molec_weight(ihet_n2o5)    = m_n2o5
molec_weight(ihet_ho2_ho2) = m_ho2
IF (ANY(molec_weight > 1.0e5)) THEN
  cmessage = ' Molecular weights not defined for all cases'
  errcode = 1
  CALL ereport('UKCA_TROP_HETCHEM',errcode,cmessage)
END IF

! ========================================================
!  1. Calculate gamma on each aerosol type for 1,nhet
! ========================================================

! Get <gamma> for N2O5 or HO2 reaction, which is a function of aerosol type,
! temperature, and RH.

IF (PrintStatus >= PrStatus_Diag) THEN
  CALL umPrint( 'Sticking_cf:')
END IF
sticking_cf(:,:,:) = 0.0e0
DO i=1,nhet
  IF ( i == ihet_n2o5 ) THEN
    DO icp=1,ncp      ! No of components
      IF (icp /= cp_nh4) THEN
        sticking_cf(:,i,icp) = cal_n2o5(icp, nbox, temp(:), rh(:))
        IF (PrintStatus >= PrStatus_Diag) THEN
          WRITE(umMessage,'(2I0,3E15.6)')                                      &
              i,icp,MAXVAL(sticking_cf(:,i,icp)),                              &
              MINVAL(sticking_cf(:,i,icp)),                                    &
              SUM(sticking_cf(:,i,icp))/SIZE(sticking_cf(:,i,icp))
          CALL umPrint(umMessage,src='ukca_trop_hetchem')
        END IF
      END IF
    END DO
  ELSE IF ( i == ihet_ho2_ho2 ) THEN
    DO icp=1,ncp
      IF (icp /= cp_nh4) THEN
        sticking_cf(:,i,icp) = cal_ho2(icp, nbox, temp(:), rh(:))
        IF (PrintStatus >= PrStatus_Diag) THEN
          WRITE(umMessage,'(2I0,E15.6)') i,icp,MAXVAL(sticking_cf(:,i,icp))
          CALL umPrint(umMessage,src='ukca_trop_hetchem')
        END IF
      END IF
    END DO
  ELSE
    cmessage = ' Error in calculating sticking_cf'
    errcode = 1
    CALL ereport('UKCA_TROP_HETCHEM',errcode,cmessage)
  END IF
END DO   ! end of loop over nhet

! So now we have an array sticking_cf(nbox,nhet,aerotypes)
!  with gamma for n2o5 and ho2 for each aerosol type.

! =================================================================
!  2. Calculate mean gamma for each modes based on partial volumes
! =================================================================

      ! Calc fractional volumes

total_vol(:,:) = 0.0e0
DO imode=1,nmodes
  IF (mode(imode)) THEN
    DO icp=1,ncp
      IF ( (component(imode,icp)) .AND. (icp /= cp_nh4) ) THEN
        total_vol(:,imode) = total_vol(:,imode) +                              &
                            pvol(:,imode,icp)
      END IF
    END DO  ! end of loop over components
  END IF
END DO  ! end of loop over modes

IF (PrintStatus >= PrStatus_Diag) THEN
  CALL umPrint( 'frac_vol:' )
END IF
frac_vol(:,:,:) = 0.0e0
DO imode=1,nmodes
  IF (mode(imode)) THEN
    DO icp=1,ncp
      IF ( (component(imode,icp)) .AND. (icp /= cp_nh4) ) THEN
        frac_vol(:,imode,icp) = pvol(:,imode,icp) /                            &
                                     total_vol(:,imode)
        IF (PrintStatus >= PrStatus_Diag) THEN
          WRITE(umMessage,'(2I0,3E15.6)')                                      &
              imode,icp,MAXVAL(frac_vol(:,imode,icp)),                         &
              MINVAL(frac_vol(:,imode,icp)),                                   &
              SUM(frac_vol(:,imode,icp))/SIZE(frac_vol(:,imode,icp))
          CALL umPrint(umMessage,src='ukca_trop_hetchem')
        END IF
      END IF
    END DO  ! end of loop over components
  END IF
END DO  ! end of loop over modes

! Calculate mode_gamma by weighting according to frac_vols,
!   e.g. gammaSO4*frac_volSO4, for each component, then sum up over each mode.

mode_gamma(:,:,:) = 0.0
DO i=1,nhet
  DO imode=1,nmodes
    IF (mode(imode)) THEN
      DO icp=1,ncp
        IF ( (component(imode,icp)) .AND. (icp /= cp_nh4) ) THEN
          mode_gamma(:,i,imode) = mode_gamma(:,i,imode) +                      &
              (frac_vol(:,imode,icp) * sticking_cf(:,i,icp))
        END IF
      END DO ! end loop over ncp
    END IF
  END DO ! end loop over nmodes
END DO  ! end loop over nhet

! So now we have mode_gamma(nbox,nhet,nmodes)

! ===========================================================
! 3. Calculate first order rate coefficients for each mode
! ===========================================================

IF (PrintStatus >= PrStatus_Diag) THEN
  CALL umPrint( 'Mode_rate:')
END IF
mode_rate(:,:,:) = 0.0e0
DO i=1,nhet
  DO imode=1,nmodes
    IF (mode(imode)) THEN
      mode_rate(:,i,imode) = cal_loss1k(nbox, sarea(:,imode),                  &
                             wetdp(:,imode), aird(:),                          &
                             mode_gamma(:,i,imode), temp(:),                   &
                             molec_weight(i))
      IF (PrintStatus >= PrStatus_Diag) THEN
        WRITE(umMessage,'(I0,3E15.6)')                                         &
            imode,MAXVAL(mode_rate(:,i,imode)),                                &
            MINVAL(mode_rate(:,i,imode)),                                      &
            SUM(mode_rate(:,i,imode))/SIZE(mode_rate(:,i,imode))
        CALL umPrint(umMessage,src='ukca_trop_hetchem')
      END IF
    END IF
  END DO
END DO

! =====================================================
! 4. Calc average rate for the box by adding mode rates
! =====================================================


total_sa(:) = 0.0e0
DO imode=1,nmodes
  IF (mode(imode)) THEN
    total_sa(:) = total_sa(:) + sarea(:,imode)
  END IF
END DO

het_rates(:,:)=0.0e0
DO i=1,nhet
  DO imode=1,nmodes
    IF (mode(imode)) THEN
      WHERE (total_sa > 1.0e-30)
        het_rates(:,i) = het_rates(:,i) + mode_rate(:,i,imode)
      END WHERE
    END IF
  END DO
END DO


IF (PrintStatus >= PrStatus_Diag) THEN
  WRITE(umMessage,'(A,3E15.6)')                                                &
                     'rc_n2o5: ',MAXVAL(het_rates(:,1)),                       &
                       MINVAL(het_rates(:,1)),                                 &
                       SUM(het_rates(:,1))/SIZE(het_rates(:,1))
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
  WRITE(umMessage,'(A,3E15.6)')                                                &
                     'rc_ho2:  ',MAXVAL(het_rates(:,2)),                       &
                       MINVAL(het_rates(:,2)),                                 &
                       SUM(het_rates(:,2))/SIZE(het_rates(:,2))
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE ukca_trop_hetchem

! ===============================================================
!                N 2 O 5   F U N C T I O N
! ===============================================================

      !==============================================================
      ! Function N2O5 computes the GAMMA sticking factor
      ! for N2O5 hydrolysis.
      !
      ! NOTES:
      !==============================================================

FUNCTION cal_n2o5(aerotype, nbox, temp, rh)

!==============================================================
! Internal function CAL_N2O5 computes the GAMMA sticking factor
! for N2O5 hydrolysis.
!==============================================================
IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN) :: aerotype     ! # denoting aerosol cpt
INTEGER, INTENT(IN) :: nbox         ! No of points

REAL,    INTENT(IN) :: temp(nbox)   ! Temperature [K]
REAL,    INTENT(IN) :: rh(nbox)     ! Relative Humidity [fraction]

! Function return value
REAL :: cal_n2o5(nbox)

! Local variables
REAL                :: ttemp(nbox)   ! dummy array for temp
REAL                :: rh_p(nbox)    ! RH as a percent
REAL                :: maxtemp       ! maximum
REAL                :: lambda(nbox)  !

CHARACTER(LEN=errormessagelength) :: cmessage ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CAL_N2O5'

!==============================================================
! N2O5 begins here!
!==============================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Convert RH to % (max = 100%)
! Values already limited to <= 1.0 - see rel_humid_frac in UKCA_MAIN1
rh_p(:)  = rh(:) * 100.0

! Special handling for various aerosols
SELECT CASE ( aerotype )

   !----------------
   ! Dust
   !----------------
CASE ( cp_du )

   !======================================================
   ! Based on unpublished Crowley work, was 0.01E0
   !------------------------------------------------------
   ! Now gamma=0.03-0.15 (Mogili et al 2006) so use
   ! a value of 0.1.
   ! Also 0.2 on saharan test dust (Karagulian et al 2006)
   ! though this is now thought to be too high (method)
   ! (hlm 2008)
   !======================================================
  cal_n2o5(:) = 0.1

  !----------------
  ! Sulfate
  !----------------
CASE ( cp_su )

   !========================================================
   ! RH dependence from Kane et al., Heterogeneous uptake of
   ! gaseous N2O5 by (NH4)2SO4, NH4HSO4 and H2SO4 aerosols
   ! J. Phys. Chem. A , 2001, 105, 6465-6470
   !========================================================
  cal_n2o5(:) = 2.79e-4 + rh_p(:)*(  1.30e-4 +                                 &
                       rh_p(:)*( -3.43e-6 +                                    &
                       rh_p(:)*(  7.52e-8 ) ) )

  !========================================================
  ! Temperature dependence factor (Cox et al, Cambridge UK)
  ! is of the form:
  !
  !          10^( LOG10( G294 ) - 0.04 * ( ttemp - 294 ) )
  ! fact = ------------------------------------------------
  !                     10^( LOG10( G294 ) )
  !
  ! Where G294 = 1e-2 and ttemp is MAX( temp, 282 ).
  !
  ! For computational speed, replace LOG10( 1e-2 ) with -2
  ! and replace 10^( LOG10( G294 ) ) with G294
  !========================================================
  ! ttemp(:) = MAX( temp(:), 282E0 ) - 294.0  [MAX isn't array valued]
  maxtemp = 282.0 - 294.0
  WHERE (temp < 282.0e0)
    ttemp = maxtemp
  ELSE WHERE
    ttemp = temp - 294.0
  END WHERE

  ! Apply temperature dependence
  cal_n2o5(:) = cal_n2o5(:)*10.0**(-2.0-4.0e-2*ttemp(:))/1.0e-2

  !----------------
  ! Black Carbon
  !----------------
CASE ( cp_bc )

   !======================================================
   ! From IUPAC
   !======================================================
  cal_n2o5(:) = 0.005e0

  !----------------
  ! Organic Carbon
  !----------------
CASE ( cp_oc , cp_so )

   !========================================================
   ! Based on Thornton, Braban and Abbatt, 2003
   ! N2O5 hydrolysis on sub-micron organic aerosol: the effect
   ! of relative humidity, particle phase and particle size
   !========================================================
  WHERE ( rh_p >= 57.0e0 )
    cal_n2o5(:) = 0.03e0
  ELSE WHERE
    cal_n2o5(:) = rh_p(:) * 5.2e-4
  END WHERE

  !----------------
  ! Sea salt
  ! accum & coarse
  !----------------
CASE ( cp_cl )

   !=======================================================
   ! Based on IUPAC recommendation
   !-------------------------------------------------------
   ! May'07 lechlm Now ELSE is changed from G=0.005E0 to
   ! gamma = rh_p * 0.0005E0 with a study by Thornton and
   ! Abbatt 2005, JOURNAL.
   !=======================================================
  WHERE ( rh_p >= 62.0 )
    cal_n2o5(:) = 0.03e0
  ELSE WHERE
    cal_n2o5(:) = rh_p(:) * 0.0005e0
  END WHERE

CASE ( cp_no3 )

   !=======================================================
   ! Based on Davis et al, 2008
   ! Atmos. Chem. Phys., 8, 5295-5311, 2008
   !=======================================================
  lambda(:) = -8.107744 + 0.04902 * rh_p(:)
  cal_n2o5(:) = 1.0 / ( 1.0 + EXP(-1.0*lambda(:)) )
  WHERE ( cal_n2o5(:) > 0.0154 )
    cal_n2o5(:) = 0.0154
  END WHERE

CASE ( cp_nn )

   !=======================================================
   ! Based on Wahner et al, 1998
   ! JGR, 103, D23, 31,103-31,112, 1998
   ! Linear interpolation between RH points
   ! Could also use Hallquist et al, 2003
   ! Phys. Chem. Chem. Phys., 2003, 5, 3453-3463
   !=======================================================
  WHERE ( rh_p(:) < 50.0 )
    cal_n2o5(:) = 0.0018
  ELSE WHERE ( ( rh_p(:) >= 50.0 ) .AND. ( rh_p(:) < 60.0 ) )
    cal_n2o5(:) = 0.0018 + (0.0032 - 0.0018) * (rh_p(:) - 50.0) / 10.0
  ELSE WHERE ( ( rh_p(:) >= 60.0 ) .AND. ( rh_p(:) < 90.0 ) )
    cal_n2o5(:) = 0.0032 + (0.0230 - 0.0032) * (rh_p(:) - 60.0) / 30.0
  ELSE WHERE ( rh_p(:) >= 90.0)
    cal_n2o5(:) = 0.023
  END WHERE

  !----------------
  ! Default
  !----------------
CASE DEFAULT
  cmessage=' Unknown aerosol surface for N2O5 hydrolysis'
  WRITE(umMessage,'(A)') cmessage
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
  WRITE(umMessage,'(A,I0)') 'AEROSOL TYPE =',aerotype
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
  errcode = aerotype
  CALL ereport('UKCA_TROP_HETCHEM_MOD:CAL_N2O5',errcode,cmessage)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END FUNCTION cal_n2o5

! ======================================================================
!                  H O 2   F U N C T I O N
! ======================================================================

FUNCTION cal_ho2(aerotype, nbox, temp, rh)

!================================================================
! Function HO2 computes the gamma sticking factor
! for HO2 self reaction.
!
! NOTES:
! (1 ) Dust values for RH<50% are needed.
! (2 ) Data for organics is from a conference poster that was not
!        published. Replace with data as it becomes available.
! (3 ) Temperature dependence on salt seems fairly good (N.B. this
!        is only for SOLID NaCl, not aq). Fits available data
!        (though this data was not published to our knowledge).
!        This temperature dependence fits well for sulphate data,
!        so has been used here due to lack of data.
!        Update when available.
! (4 ) Data for various relative humidities available now (04/08)
!        from Taketani et al 2008 on sulfate aerosol.
! (5 ) It should be noted that it may be possible for gamma to
!        reach very high values if TMIs are present. However,
!        it is unclear if there are any observations of this.
!        See notes on dust in the code further down.
! (hlm 3/3/09)
!================================================================

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN) :: aerotype         ! # denoting aerosol cpt
INTEGER, INTENT(IN) :: nbox             ! No of points

REAL,    INTENT(IN) :: temp(nbox)     ! Temperature [K]
REAL,    INTENT(IN) :: rh(nbox)       ! Relative Humidity [fraction]

! Function return value
REAL :: cal_ho2(nbox)

! Local variables
REAL                :: rh_p(nbox)     ! dummy array for RH
REAL                :: factt(nbox)    ! for calc temp dependence
REAL                :: facth(nbox)    ! for calc rh dependence

CHARACTER(LEN=errormessagelength) :: cmessage ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CAL_HO2'

!==============================================================
! HO2 begins here!
!==============================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Convert rh to % (max = 100%)  !! RH is already <= 1.0
rh_p(:)  = rh(:) * 100.0

! Special handling for various aerosols
SELECT CASE ( aerotype )

   !----------------
   ! Dust
   !----------------
CASE ( cp_du )

   !=========================================================
   ! Hanel '76 (reported in Dentener etal '96). If rh>50%,
   ! then dust has enough water to allow gamma=0.1.
   !-------------------------------------------------------
   ! N.B. this should be updated as more information becomes
   ! available. Note that if TMIs (such as Cu and Fe) are
   ! present in high enough concentrations, then large uptake
   ! may be seen (gamma~0.8). See Thornton et al, 2008, JGR.
   !=========================================================
  WHERE ( rh_p >= 50.0 )
    cal_ho2(:) = 0.1e0
  ELSE WHERE
     !update when data available
    cal_ho2(:) = 0.05e0
  END WHERE

  !----------------
  ! Sulfate
  !----------------
CASE ( cp_su )

   !===========================================================
   ! Applying the temperature dependence for NaCl, it seems to
   ! fit the data very well, so use this.
   !
   ! Data from Cooper & Abbatt '96 (Strat conditions) and
   ! Thornton & Abbatt '05.
   !
   ! Now data from Taketani et al 2008.
   ! Measurements on (NH4)2SO4 at 296K, tentatively recommend
   ! gamma = 0.15, though there is dependence with rh.
   !===========================================================

   ! rh_p factor doesn't work below 35%rh
  WHERE ( rh_p <= 35.0 )
    cal_ho2(:) = 0.01e0
  ELSE WHERE
    ! Calculate temperature factor. Fix rel to 0.11
    factt(:) = (5.66e-5 * EXP( 1560.0 / temp(:) ) ) / 0.11
    ! Calc rh_p factor from Taketani et al 2008
    facth(:) = (-2.86 * EXP(-0.078*rh_p(:)) + 0.192 ) / 0.11
    ! Calc gamma
    cal_ho2(:) = 0.11 * factt(:) * facth(:)
  END WHERE

  !Test with recommendations by Thornton and Abbatt 2008
  ! cal_ho2(:) = 0.15

  !----------------
  ! Black Carbon
  !----------------
CASE ( cp_bc )

   !===========================================================
   ! Saathoff et al 2001. gamma=<0.01 on soot.
   !----------------------------------------------------------
   ! Ivanov (ref as dust) has measurements on soot, but can't
   ! find available anywhere.
   !===========================================================
  cal_ho2(:) = 0.01e0

  !----------------
  ! Organic Carbon
  !----------------
CASE ( cp_oc , cp_so )

   !===========================================================
   ! Ivanov et al. Conference poster. Range of organic surfaces
   ! gamma=1d-4 to 5d-2 @ room temp.
   !===========================================================
   ! mid-value of available data
  cal_ho2(:) = 2.5e-2

  !----------------
  ! Sea salt
  ! accum & coarse
  !----------------
CASE ( cp_cl )

   !===========================================================
   ! Taketani etal ?. For NaCl, gamma<0.01 @ 20%rh solid
   ! (at room temp)                   0.05 @ 45%rh aq
   !---------------------------------------------------------
   ! Remorov etal 2002. Study on solid NaCl
   ! gamma=(5.66 pm 3.62)d-5 exp((1560 pm 140)/temp)
   !===========================================================

   ! Temperature dependence (for solid NaCl)
  WHERE ( rh_p >= 62.0 )
    cal_ho2(:) = 0.05e0
  ELSE WHERE
    cal_ho2(:) = 5.66e-5 * EXP( 1560.0 / temp(:) )
  END WHERE

CASE ( cp_no3 )

   !===========================================================
   ! George et al. 2013
   ! Phys. Chem. Chem. Phys.,2013, 15, 12829
   !===========================================================
  cal_ho2(:) = 0.005

CASE ( cp_nn )

   !===========================================================
   ! Nominal value based on NaCl and NH4NO3 composite
   ! George et al. 2013
   ! Phys. Chem. Chem. Phys.,2013, 15, 12829
   !===========================================================
  cal_ho2(:) = 0.015

  !----------------
  ! Default
  !----------------
CASE DEFAULT
  cmessage='Unknown aerosol surface for HO2 self reaction'
  WRITE(umMessage,'(A)') cmessage,' aerotype=',aerotype
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
  WRITE(umMessage,'(A,I0)') 'AEROSOL TYPE =',aerotype
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
  errcode = aerotype
  CALL ereport('UKCA_TROP_HETCHEM_MOD:CAL_HO2',errcode,cmessage)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION cal_ho2

! =========================================================================
!                F I R S T   O R D E R   L O S S
! =========================================================================

FUNCTION cal_loss1k( nbox, sarea, wetdp, aird, mode_gamma,                     &
                     temp, molec_wt_s)

! =================================================================
!  Function CAL_LOSS1K calculates the 1st-order loss rate of
!  species on wet aerosol surface.
! =================================================================
! The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
! is computed as:
!
!      K [1/s] = sarea / [ wetdp/dfkg + 4./(gamma * xmms) ]
!
! where XMMS = Mean molecular speed [cm/s] = sqrt(8R*T/pi*M) for Maxwell
!              [Formulae gives [m/s] if R in J/(mol.K) and M in kg/mol]
!       DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
!
! NOTES:
! (1 ) Return a default value if RADIUS is zero (i.e. is smaller
!      than a very small number), except when the calculation is
!      done on the surface of CLASSIC aerosols.
! *******************************************************************

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN) :: nbox             ! No of points

REAL, INTENT(IN) :: sarea(nbox)         ! surface area of wet
                                        !  aerosol/vol of air [cm2/cm3]
REAL, INTENT(IN) :: wetdp(nbox)         ! radius of wet aerosol [cm],
                                        ! order of 0.01-10 um
REAL, INTENT(IN) :: aird(nbox)          ! density of air [#/cm3]
REAL, INTENT(IN) :: mode_gamma(nbox)    ! Uptake coefficient [unitless] ~ 0.1
REAL, INTENT(IN) :: temp(nbox)          ! temperature [K]
REAL, INTENT(IN) :: molec_wt_s          ! molecular weight of species [g/mole]

! Input from module
!      rmol                ! gas constant [J/(mol.K)]
!      avogadro            ! avogadros number
!      m_air               ! Molec Wt. of air kg/mol

! Function return value
REAL :: cal_loss1k(nbox)       ! Reaction coefficient on aerosol surfaces

! Local variables
REAL     :: cal_loss1k_min     ! Minimum value returned by function
REAL     :: sarea_min          ! Minimum surface area to return cal_loss1k_min

REAL     :: dfkg(nbox)                  ! Gas phase diffusion coefficient
REAL     :: rho_a(nbox)                 ! Air density in [kg/m^3]
REAL     :: mw_s                        ! molec_weight of species [kg/mol]
REAL, PARAMETER :: cf=1.0e-2            ! conversion factor 1/(m/s) to 1/(cm/s)
REAL, PARAMETER :: dq=4.5e-10           ! Diameter of air molecule [m]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CAL_LOSS1K'

!=============================================================
! CAL_LOSS1K begins here!
!=============================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Convert molecular weight from g/mol to kg/mol
mw_s = molec_wt_s/1.0e3

! Convert air density from cm-3 (array aird) to kg/m3 (array rho_a)
rho_a(:) = aird(:)*m_air*1.0e6/avogadro

IF (PrintStatus >= PrStatus_Diag) THEN
  WRITE(umMessage,'(A,2E15.6)') 'RHO_A: ',MAXVAL(rho_a),MINVAL(rho_a)
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
END IF
dfkg(:)=0.0  ! for debug only

! Set minimum threshold of surface area for which the reaction rate is
! expected to be very small. Also set the corresponding return value of
! cal_loss1k. Note that this is slightly different for GLOMAP-mode and
! CLASSIC aerosols.
IF (ukca_config%l_ukca_classic_hetchem) THEN
  ! Reaction on the surface of CLASSIC aerosols. If a given aerosol
  ! type/mode is not being modelled then the surface area is already
  ! filled with zeros and the resulting reaction rate should be set
  ! to zero too.
  sarea_min      = 1.0e-30
  cal_loss1k_min = 0.0
ELSE
  ! Set a minimum value of the reaction rate for GLOMAP-mode aerosols
  ! if the surface area is negative.
  sarea_min      = 0.0
  cal_loss1k_min = 1.0e-3
END IF

WHERE (sarea < sarea_min .OR. wetdp < 1.0e-30 .OR. mode_gamma < 1.0e-30)
  ! Use default value of reaction rate for the grid cells with
  ! surface area, radius or uptake coefficient below given thresholds
  ! (basically around zero or negative)
  cal_loss1k = cal_loss1k_min

ELSE WHERE
  ! dfkg = Gas phase diffusion coeff [cm2/s] (order of 0.1)
  ! see  Bauer et al., 2004, JGR 109, D02304 (eqn 7), with
  ! factor of 1e4 => [cm^2/s]
  dfkg(:) = (3.0e4/(8.0 * avogadro * dq**2 * rho_a(:))) *                      &
       SQRT((rmol*temp(:)*m_air/(2.0*pi))*(mw_s + m_air)/mw_s)
  !
  cal_loss1k(:) = sarea(:)/( wetdp(:)/dfkg(:) +                                &
          (cf/mode_gamma(:))*SQRT(2.0*pi*mw_s/(rmol*temp(:))))
END WHERE

IF (PrintStatus >= PrStatus_Diag) THEN
  WRITE(umMessage,'(A,2E15.6)') 'dfkg: ',MAXVAL(dfkg),MINVAL(dfkg)
  CALL umPrint(umMessage,src='ukca_trop_hetchem')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END FUNCTION cal_loss1k

END MODULE ukca_trop_hetchem_mod
