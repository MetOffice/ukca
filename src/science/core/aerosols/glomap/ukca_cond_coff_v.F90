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
!    Subroutine which calculates the condensation coefficients
!    for condensable cpt vapour condensing on a particle with
!    radius RP. Includes a switch to either use Fuchs (1964)
!    or Fuchs-Sutugin (1971).
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
MODULE ukca_cond_coff_v_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_COND_COFF_V_MOD'

CONTAINS

SUBROUTINE ukca_cond_coff_v(nv,mask,rp,tsqrt,airdm3,rhoa,                      &
             mmcg,se,dmol,ifuchs,cc,sinkarr,pmid,t,difvol,idcmfp)
!----------------------------------------------------------------
!
! Purpose
! -------
! Subroutine which calculates the condensation coefficients
! for condensable cpt vapour condensing on a particle with
! radius RP. Includes a switch to either use Fuchs (1964)
! or Fuchs-Sutugin (1971).
!
! Inputs
! ------
! NV     : total number of values
! MASK   : Mask where to calculate values
! RP     : Radius of aerosol particle
! TSQRT  : Square-root of mid-level air temperature (K)
! AIRDM3 : Number density of air (m3)
! RHOA   : Air density (kg/m3)
! MMCG   : Molar mass of condensing gas (kg per mole)
! SE     : Sticking efficiency [taken as 0.3 from Raes et al 1992]
! DMOL   : Molecular diameter of condensable (m)
! IFUCHS : Switch for Fuchs (1964) or Fuchs-Sutugin (1971)
! PMID   : Centre level pressure (Pa)
! T      : Centre level temperature (K)
! DIFVOL : Diffusion volume for H2SO4 or Sec_Org
! IDCMFP : Switch : diffusion/mfp  (1=as Gbin v1, 2=as Gbin v1_1)
!
! Outputs
! -------
! CC     : Condensation coeff for cpt onto ptcl (m^3s^-1)
!
! Local Variables
! ---------------
! VEL_CP   : Thermal velocity of condensable gas (ms-1)
! MFP_CP   : Mean free path of condensable gas (m)
! DCOFF_CP : Diffusion coefficient of condensable gas in air (m2s-1)
! KN       : Knudsen number
! FKN      : Correction factor for molecular effects
! AKN      : Corr fac for limitations in interfacial mass transport
! DENOM    : Denominator in Fuchs (1964) condensation coeff expression
! ZZ       : Ratio of molar masses of condensable gas and air
! TERM1,TERM2,.. : Terms in evaluation of condensation coefficient
!
! Constants from shared modules
! -----------------------------
! BOLTZMANN : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
! RGAS      : Dry air gas constant = 287.05 Jkg^-1 K^-1
! RMOL      : Universal gas constant = 8.314 K mol^-1 K^-1
! AVOGADRO  : Avogadros constant (mol-1)
!
!--------------------------------------------------------------------
USE chemistry_constants_mod, ONLY: boltzmann, avogadro
USE ukca_um_legacy_mod, ONLY: rmol, pi, rgas => r
USE ukca_types_mod,   ONLY: log_small
USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine interface:
INTEGER, INTENT(IN) :: nv
INTEGER, INTENT(IN) :: ifuchs
INTEGER, INTENT(IN) :: idcmfp
LOGICAL (KIND=log_small), INTENT(IN) :: mask(nv)
REAL, INTENT(IN)    :: rp(nv)
REAL, INTENT(IN)    :: tsqrt(nv)
REAL, INTENT(IN)    :: airdm3(nv)
REAL, INTENT(IN)    :: rhoa(nv)
REAL, INTENT(IN)    :: mmcg
REAL, INTENT(IN)    :: se
REAL, INTENT(IN)    :: dmol
REAL, INTENT(IN)    :: pmid(nv)
REAL, INTENT(IN)    :: t(nv)
REAL, INTENT(IN)    :: difvol
REAL, INTENT(OUT)   :: cc(nv)
REAL, INTENT(OUT)   :: sinkarr(nv)
!
! .. Local variables
REAL    :: vel_cp(nv)
REAL    :: mfp_cp(nv)
REAL    :: dcoff_cp(nv)
REAL    :: kn(nv)
REAL    :: fkn(nv)
REAL    :: akn(nv)
REAL    :: denom(nv)
REAL    :: zz
REAL    :: term1
REAL    :: term2
REAL    :: term3
REAL    :: term4
REAL    :: term5
REAL    :: term6
REAL    :: term7
REAL    :: term8
REAL    :: dair

! Molar mass of dry air (kg/mol)
REAL :: mm_da  ! =avogadro*boltzmann/rgas

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_COND_COFF_V'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

term1=SQRT(8.0*rmol/(pi*mmcg))
! .. used in calcn of thermal velocity of condensable gas

mm_da = avogadro*boltzmann/rgas
zz=mmcg/mm_da
term2=1.0/(pi*SQRT(1.0+zz)*dmol*dmol)
! .. used in calcn of MFP of condensable gas (S&P, pg 457, eq 8.11)

term3=(3.0/(8.0*avogadro*dmol*dmol))
term4=SQRT((rgas*mm_da*mm_da/(2.0*pi))*((mmcg+mm_da)/mmcg))
term5=term3*term4
! .. used in calcn of diffusion coefficient of condensable gas
!
term6=4.0e6*pi
! .. used in calculation of condensation coefficient
!
term7=SQRT((1.0/(mm_da*1000.0))+(1.0/(mmcg*1000.0)))
dair=19.7 ! diffusion volume of air molecule (Fuller et al, Reid et al)
term8=(dair**(1.0/3.0)+difvol**(1.0/3.0))**2
! .. used in new calculation of diffusion coefficient

cc(:)=0.0
sinkarr(:)=0.0

! .. Calc. thermal velocity of condensable gas
WHERE (mask(:)) vel_cp(:)=term1*tsqrt(:)

! .. Calculate diffusion coefficient of condensable gas
IF (idcmfp == 1) WHERE(mask(:)) dcoff_cp(:)=term5*tsqrt(:)/rhoa(:)
IF (idcmfp == 2) WHERE(mask(:)) dcoff_cp(:)=                                   &
          (1.0e-7*(t(:)**1.75)*term7)/(pmid(:)/101325.0*term8)

!! n.b. IDCMFP=2 method is as below as coded in GLOMAP-bin:
!!        DCOFF_CP(:)=(1.0E-7*(T(:)**1.75)*sqrt((1.0/(MM_DA*1000.0))+     &
!!          (1.0/(MMCG*1000.0))))/                                        &
!!          (PMID(:)/101325.0*(DAIR**(1.0/3.0)+DIFVOL**(1.0/3.0))**2)

! .. Calc. mean free path of condensable gas
IF (idcmfp == 1) WHERE(mask(:)) mfp_cp(:)=term2/airdm3(:) ! Gb v1
IF (idcmfp == 2) WHERE(mask(:)) mfp_cp(:)=                                     &
                            3.0*dcoff_cp(:)/vel_cp(:) ! as Gb v1_1

IF (ifuchs == 1) THEN
  !      If IFUCHS=1 use basic Fuchs (1964) expression
  WHERE (mask(:))
    denom(:)=                                                                  &
   (4.0*dcoff_cp(:)/(se*vel_cp(:)*rp(:)))+(rp(:)/(rp(:)+mfp_cp(:)))
    !       Calc. condensation coefficient
    cc(:)=term6*dcoff_cp(:)*rp(:)/denom(:)
    sinkarr(:)=rp(:)*1.0e6/denom(:)
  END WHERE
END IF

IF (ifuchs == 2) THEN
  !      If IFUCHS=2 use basic Fuchs-Sutugin (1971) expression
  WHERE (mask(:))
    !       Calculate Knudsen number of condensable gas w.r.t. particle
    kn(:)=mfp_cp(:)/rp(:)
    !       Calc. corr. factor for molecular effects
    fkn(:)=(1.0+kn(:))/(1.0+1.71*kn(:)+1.33*kn(:)*kn(:))
    !       Calc. corr. factor for limitations in interfacial mass transport
    akn(:)=1.0/(1.0+1.33*kn(:)*fkn(:)*(1.0/se-1.0))
    !       Calc. condensation coefficient
    cc(:)=term6*dcoff_cp(:)*rp(:)*fkn(:)*akn(:)
    sinkarr(:)=rp(:)*1.0e6*fkn(:)*akn(:)
  END WHERE
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_cond_coff_v
END MODULE ukca_cond_coff_v_mod
