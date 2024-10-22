! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
! *********************************************************************
! To calculate the vapour pressure of H2SO4 before the
! condensation routine. From Ayers (1980).
!
! Inputs
! NBOX Number of gridboxes
! T    Temperature (K)
! S    Specific humidity (kg/kg)
! SM   grid box mass of air (kg)
! PMID Air pressure at mid point (Pa)
! RP   ?
!
!
! Calculated by subroutine
! TLOG  Natural log of T
! UST   1/T
! BH2O  Partial pressure of water vapour (atm)
! PWLOG Natural log of BH2O
! PATM  Pressure of air (atm)
! WS    Composition of aerosol (wt fraction H2SO4)
! A,B,XSB,MSB,MUH2SO4 all avriables used to calculate PH2SO4
! PH2SO4 Vapour pressure of H2SO4 (atm)
!
! Outputs
! WTS   Catch for WS being less than 40%
! RHOSOL_STRAT Density of particle solution [excl. insoluble cpts] (kg/m3)
! **********************************************************************

MODULE ukca_vapour_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_VAPOUR_MOD'

CONTAINS

SUBROUTINE ukca_vapour(nbox,t,pmid,s,rp,wts,rhosol_strat)

USE ukca_config_constants_mod, ONLY: rmol, rho_so4
USE ukca_um_legacy_mod,  ONLY: log_v, exp_v, powr_v, l_glomap_clim_radaer
USE ukca_constants,      ONLY: mmsul
USE parkind1,            ONLY: jprb, jpim
USE yomhook,             ONLY: lhook, dr_hook

USE ukca_config_specification_mod, ONLY: glomap_config

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN) :: nbox               ! Number of gridboxes

REAL, INTENT(IN)    :: t(nbox)            ! Temperature (K)
REAL, INTENT(IN)    :: pmid(nbox)         ! Air pressure at mid point (Pa)
REAL, INTENT(IN)    :: s(nbox)            ! Specific humidity (kg/kg)
REAL, INTENT(IN)    :: rp(nbox)           ! ?
REAL, INTENT(OUT)   :: wts(nbox)          ! MAX(41.0, WS*100)
                                          ! or MIN(99.0, MAX(41.0,WS*100))
REAL, INTENT(OUT)   :: rhosol_strat(nbox) ! Density of particle solution (kg/m3)

! Local variables
INTEGER :: jl

REAL :: bh2o(nbox), patm
REAL :: ws(nbox)
REAL :: a,b,xsb(nbox),msb(nbox),pwlog(nbox),tlog(nbox),ust(nbox)
REAL :: muh2so4(nbox)
REAL :: ph2so4(nbox)

! select density value from LUT

INTEGER :: k
INTEGER :: round(nbox)
REAL    :: t_diff(nbox)
!  Data for density @253K, per Kelvin difference at weight % of H2SO4
INTEGER, PARAMETER :: ndata = 12 ! no of data points

! percentage of H2SO4 present in aerosol
INTEGER, PARAMETER :: percent(ndata) =                                         &
         [40,45,50,55,60,65,70,75,80,85,90,95]

! per Kelvin difference in density
REAL, PARAMETER :: k_diff(ndata) =                                             &
 [0.80,0.82,0.84,0.86,0.87,0.89,0.92,0.95,0.98,1.02,1.06,1.10]

! density @ 253K kgm-3
REAL, PARAMETER :: data253(ndata) =                                            &
 [1333.8,1381.1,1431.0,1483.3,1537.3,1592.4,1647.6,                            &
   1701.6,1753.0,1800.1,1840.9,1873.2]

REAL :: c,d
REAL, PARAMETER :: xsb_eps=1.0e-6
REAL, PARAMETER :: p0=101325.0

REAL, PARAMETER :: ks1=-21.661
REAL, PARAMETER :: ks2= 2724.2
REAL, PARAMETER :: ks3= 51.81
REAL, PARAMETER :: ks4=-15732.0
REAL, PARAMETER :: ks5=47.004
REAL, PARAMETER :: ks6=-6969.0
REAL, PARAMETER :: ks7=-4.6183
REAL, PARAMETER :: surften=0.0728
REAL, PARAMETER :: bminatm=2.0e-8
REAL, PARAMETER :: bmaxatm=2.0e-6

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_VAPOUR'

REAL :: term1
REAL :: kelvin(nbox)
REAL :: kelvin_out(nbox)
REAL :: tmp2(nbox)
REAL :: tmp2_out(nbox)
REAL :: wts_m40(nbox)
REAL :: tmp4(nbox)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Loop around all the boxes
DO jl=1,nbox

  patm=pmid(jl)/p0

  ! calculate vmr and then atm of water vapour

  bh2o(jl)=1.609*s(jl)*patm
  IF (bh2o(jl) < bminatm) bh2o(jl)=bminatm
  IF (bh2o(jl) > bmaxatm) bh2o(jl)=bmaxatm
  ust(jl)=1/t(jl)

END DO

!  log(T) and  water vapour
CALL log_v(nbox,t,tlog)
CALL log_v(nbox,bh2o,pwlog)

DO jl=1,nbox

  ! Calculate H2SO4 vapour pressure from expression for the H2SO4 vp
  !  over H2SO4/H2O solutions from paper of Ayers + al., GRL, 7, 433-436, 1980
  !  The Giaugue expression for chemical potential has been fitted
  !  to a simple expression (see the paper). It is likely to be
  !  a bit inaccurate for use at low strat conditions, but should be reasonable

  !  log(T)and water vapour - now done outside loop


  !  The H2SO4/H2O pure solution concentration

  !  Mole fraction of H2SO4 in the binary solution
  a=ks1 + ks2*ust(jl)
  b=ks3 + ks4*ust(jl)

  c=ks5 + ks6*ust(jl) + ks7*tlog(jl) - pwlog(jl)
  d=a*a - 4.0*b*c
  IF (d < 0.0) d=0.0

  !!   XSB(JL)=(-A -SQRT(A*A -
  !!             4.0*B*(KS5 + KS6*UST(JL) + KS7*TLOG(JL) - PWLOG(JL))))/(2.0*B)

  xsb(jl)=(-a -SQRT(d))/(2.0*b)
  IF (xsb(jl) < xsb_eps) xsb(jl)=xsb_eps

  msb(jl)=55.51*xsb(jl)/(1.0-xsb(jl))
  ws(jl) = msb(jl)*0.098076/(1.0 + msb(jl)*0.098076)

  ! wt % of H2SO4
  IF (glomap_config%l_fix_neg_pvol_wat .OR. l_glomap_clim_radaer) THEN
    ! Upper limit of 99% prevents negative concentration fields
    wts(jl)=MIN( 99.0, MAX(41.0,ws(jl)*100.0) )
  ELSE
    ! 41.0% threshold originates in the following paper
    ! Carslaw et al., GRL, VOL. 22, NO. 14, PAGES 1877-1880, 1995
    wts(jl)=MAX( 41.0, ws(jl)*100.0 )
  END IF
  wts_m40(jl)=wts(jl)-40.0

END DO  ! end loop around boxes

CALL powr_v(nbox,wts_m40,0.1,tmp4)
tmp2(:)=360.0/t(:)
CALL log_v(nbox,tmp2,tmp2_out)
term1 = 2.0*surften*mmsul
kelvin(:)= term1/(rho_so4*rmol*t(:)*rp(:))
CALL exp_v(nbox,kelvin,kelvin_out)

DO jl=1,nbox

  muh2so4(jl)=4.184*(1.514e4 - 286.0*wts_m40(jl) + 1.080*                      &
                    wts_m40(jl)**2 - 3941.0/tmp4(jl) )

  ! Equilibrium H2SO4 vapour pressure (atm) using Ayers et al 1980.GRL.7.433-436
  !  PH2SO4(JL)=EXP(-10156.0*UST(JL) + 16.2590 - MUH2SO4(JL)/(8.314*T(JL)))

  ! Addendum to Ayers et al (1980)
  !  by Kulmala and Laaksonen.1990.J.Chem.Phys.93.696-701
  !   K&L say Ayers cannot be used outside T range 338-445K.
  !   Correction using Ayers at 360K (ln vp=-10156/360+16.259)=-11.9521

  ph2so4(jl)= -11.9521+10156.0*(-(1.0/t(jl))+(1.0/360.0)+                      &
         0.38/(905.0-360.0)*(1.0+tmp2_out(jl)-(360.0/t(jl))))
  ph2so4(jl)= EXP(ph2so4(jl) - muh2so4(jl)/(8.314*t(jl)))

  !  New bit of code to calculate the density of an H2SO4-H2O mixture, using
  !   Martin et al. 2000. GRL.27. 197-200.
  !   Uses weight % of H2SO4 calculated above.

  !  1300 kg m-3 is about the minimum RHOSOL_STRAT can be.
  rhosol_strat(jl)=1300.0

  !  need to take WTS and work out nearest integer of 5 to correspond to data
  round(jl)=(NINT(wts(jl)/5))*5

  ! find temperature difference - T is allowed to be negative here.
  t_diff(jl)=253.0-t(jl)

  ! read data, find value of n that corresponds to % H2SO4 data number
  DO k=1,ndata
    IF (round(jl)==percent(k)) THEN
      rhosol_strat(jl)=data253(k)+k_diff(k)*t_diff(jl)
    END IF ! round(jl)=percent(k)
  END DO ! loop around ndata

END DO  ! end loop around boxes


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_vapour
END MODULE ukca_vapour_mod
