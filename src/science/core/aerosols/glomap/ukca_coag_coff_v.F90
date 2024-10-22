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
!    Calculate coagulation coefficient for modes coagulating
!    with radii RI, RJ, volumes VI, VJ.
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
MODULE ukca_coag_coff_v_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_COAG_COFF_V_MOD'

CONTAINS

SUBROUTINE ukca_coag_coff_v(nv,mask,ri,rj,vi,vj,rhoi,rhoj,                     &
     mfpa,dvisc,t,kij,coag_on,icoag)
!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculate coagulation coefficient for modes coagulating
! with radii RI, RJ, volumes VI, VJ.
!
! 4 possible methods of calculation chosen via switch ICOAG.
!
! 1) Method as Jacobson, "Fundamentals of Atmospheric Modelling", p. 446
!    for Brownian coagulation between ptls in the transition regime
!
!      K_ij=4*PI*(r_i+r_j)*(D_i+D_j)/(TERM1+TERM2)
!
!    where r_i,r_j are (wet) particle radii
!          D_i,D_j are particle diffusion coefficients
!
!          D_i=ACunn_i*kT/(6*PI*r_i*DVISC)
!
!          ACunn_i is the Cunningham slip correction
!                = 1.0+Kn_i*(1.257+0.4*EXP{-1.1/Kn_i})
!
!          Kn_i is the Knudsen number = MFPA/r_i
!
!          MFPA=2*DVISC/(RHOA*vel_i)
!
!          TERM1=(r_i+r_j)/(r_i+r_j+sqrt{delta_i^2+delta_j^2}
!          TERM2=4*(D_i+D_j)/( sqrt{vel_i^2+vel_j^2}*(r_i+r_j) )
!
!          delta_i,delta_j are mean distances from sphere centre after
!                          travelling distance of particle mfp (mfpp).
!          delta_i=( (2*r_i+mfpp_i)^3 - (4*r_i^2 + mfpp_i^2)^{3/2} )
!
!          mfpp_i= mean free path of particle = 2*D_i/PI/vel_i
!
!          vel_i,vel_j are thermal speed of particles
!                vel_i = sqrt{8kT/PI/mass_i}
!
!          where k is Boltzmann's constant
!                mass_i is the mass of a particle = rho_i*vol_i
!
! 2) Method as in HAM/M7, see Stier et al (2005), Vignati et al (2004)
!
!      K_ij=16*pi*RMID*DMID/(4*DMID/VMID/RMID + RMID/(RMID+DELTAMID))
!
!    where DMID,VMID,DELTAMID are D_i, vel_i and delta_i
!          evaluated for a particle with mean radius RMID=(r_i+r_j)/2
!
!          n.b. I have been informed since then that although this
!          is what is described in the Stier et al (2005) and
!          Vignati et al (2004), that in the latest HAM module,
!          DMID,VMID and DELTAMID are now evaluated more completely as:
!
!          DMID={D(r_i)+D(r_j)}/2.0
!          VMID=sqrt{vel_i^2+vel_j^2}
!          DELTAMID=sqrt{delta_i^2+delta_j^2}/sqrt{2.0}
!
! 3) Method as in UM sulfur scheme as in UMDP 20
!
!      K_ij=(2*k*T/3*DVISC)*(r_i+r_j)/{1/r_i + 1/r_j +
!                                      MFP*ACunn*((1/r_i)^2+(1/r_j)^2)}
!
!    where k=Boltzmann's constant
!          DVISC=dynamic viscosity of air
!          T=air temperature
!          MFP=MFPA=mean free path of air=6.6e-8*p0*T/(p*T0)
!                             [p0=1.01325e5 Pa,T0=293.15K]
!          ACunn=Cunningham slip correction assumed = 1.59 for all modes
!
! 4) Method as in 3) but with MFPP_i,MFPPj and ACunn_i,ACunn_j and MFPA
!                             calculated explicitly as in 1)
!
!      K_ij=(2*k*T/3*DVISC)*(r_i+r_j)/{1/r_i + 1/r_j +
!                         (MFPP_i*ACunn_i/r_i^2)+(MFPP_j*ACunn_j/r_j^2)}
!
! Vector version !!!!!!!!!!!!!!!!!!!!!
!
! Parameters
! ----------
! None
!
! Inputs
! ------
! NV     : Number of values
! MASK   : Mask where to calculate values
! RI     : Geometric mean radius for mode I
! RJ     : Geometric mean radius for mode J
! VI     : Volume for ptcl with r=RI
! VJ     : Volume for ptcl with r=RJ
! T      : Air temperature (K)
! RHOI   : Density of aerosol particle for mode I (kgm^-3)
! RHOJ   : Density of aerosol particle for mode J (kgm^-3)
! MFPA   : Mean free path of air (m)
! DVISC  : Dynamic viscosity of air (kg m^-1 s^-1)
! COAG_ON: Switch for turning coagulation off
! ICOAG  : KIJ method (1:GLOMAP, 2: M7, 3: UMorig, 4:UMorig MFPP)
!
! Outputs
! -------
! KIJ    : Coagulation coefficient for I-J (cm^3 s^-1)
!
! Local Variables
! ---------------
! RMID   : Arithmetic mean of I and J mode mean radii (m)
! VMID   : Volume corresponding to radius RMID
! VEL    : Thermal velocity of ptcl of size RMID (m/s)
! DCOEF  : Diffusion coefficient of ptcl of size RMID (m^2 s^-1)
! MFPP   : Mean free path of aerosol ptcl of size RMID (m)
! VELI   : Thermal velocity of ptcl in mode I (m/s)
! DCOEFI : Diffusion coefficient of ptcl in mode I (m^2 s^-1)
! MFPPI  : Mean free path of aerosol ptcl in mode I (m)
! DELI   : Coeff. in calc. of coag kernel for mode I
! VELJ   : Thermal velocity of ptcl in mode J (m/s)
! DCOEFJ : Diffusion coefficient of ptcl in mode J (m^2 s^-1)
! MFPPJ  : Mean free path of aerosol ptcl in mode J (m)
! DELJ   : Coeff. in calc. of coag kernel for mode J
! KNI    : Knudsen number of particle in mode I
! CCI    : Cunningham correction factor for ptcl in mode I
! KNJ    : Knudsen number of particle in mode J
! CCJ    : Cunningham correction factor for ptcl in mode J
! TERMV1 : Vector term in calculation of KIJ
! TERMV2 : Vector term in calculation of KIJ
! TERMV3 : Vector term in calculation of KIJ
! TERMV4 : Vector term in calculation of KIJ
!
! Inputted by module CHEMISTRY_CONSTANTS
! ---------------------------------
! BOLTZMANN : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
!
!--------------------------------------------------------------------
USE ukca_constants,     ONLY: pi
USE ukca_types_mod,     ONLY: log_small
USE ukca_config_constants_mod, ONLY: boltzmann
USE parkind1,           ONLY: jprb, jpim
USE yomhook,            ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface:
INTEGER, INTENT(IN) :: nv
INTEGER, INTENT(IN) :: coag_on
INTEGER, INTENT(IN) :: icoag
REAL, INTENT(IN)    :: ri(nv)
REAL, INTENT(IN)    :: rj(nv)
REAL, INTENT(IN)    :: vi(nv)
REAL, INTENT(IN)    :: vj(nv)
REAL, INTENT(IN)    :: rhoi(nv)
REAL, INTENT(IN)    :: rhoj(nv)
REAL, INTENT(IN)    :: mfpa(nv)
REAL, INTENT(IN)    :: dvisc(nv)
REAL, INTENT(IN)    :: t(nv)
LOGICAL(KIND=log_small), INTENT(IN) :: mask(nv)
REAL, INTENT(OUT)   :: kij(nv)
!
! Local variables:
REAL :: dpi(nv)
REAL :: dpj(nv)
REAL :: veli(nv)
REAL :: dcoefi(nv)
REAL :: mfppi(nv)
REAL :: deli(nv)
REAL :: kni(nv)
REAL :: cci(nv)
REAL :: velj(nv)
REAL :: dcoefj(nv)
REAL :: mfppj(nv)
REAL :: delj(nv)
REAL :: knj(nv)
REAL :: ccj(nv)
REAL :: vel(nv)
REAL :: dcoef(nv)
REAL :: mfpp(nv)
REAL :: kn(nv)
REAL :: cc(nv)
REAL :: rtot(nv)
REAL :: dtot(nv)
REAL :: rmid(nv)
REAL :: vmid(nv)
REAL :: rhomid(nv)
REAL :: termv1(nv)
REAL :: termv2(nv)
REAL :: termv3(nv)
REAL :: termv4(nv)
REAL :: term1
REAL :: term2
REAL :: term3
REAL :: term4
REAL :: term5

REAL, PARAMETER :: ukca_mfp_ref=6.6e-8

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_COAG_COFF_V'

!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
kij(:) = 0.0
!
IF (coag_on == 0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
!
term1=8.0*boltzmann/pi
term2=boltzmann/(6.0*pi)
term3=8.0/pi
term4=4.0e6*pi
term5=16.0e6*pi
!
!  .. For ICOAG=1, calculate KIJ following the full method
!  .. as GLOMAP does for each bin (see header above).
IF (icoag == 1) THEN
  WHERE (mask(:))

    dpi(:)=2.0*ri(:)
    veli(:)=SQRT(term1*t(:)/(rhoi(:)*vi(:)))
    kni(:)=mfpa(:)/ri(:)
    cci(:)=1.0+kni(:)*(1.257+0.4*EXP(-1.1/kni(:)))
    ! .. n.b. use 1.257,0.4,1.1 as in Seinfeld & Pandis pg. 465
    ! .. following Allen & Raabe (1982) reanalysis of Millikan (1923)
    dcoefi(:)=term2*cci(:)*t(:)/(ri(:)*dvisc(:))
    mfppi(:)=term3*dcoefi(:)/veli(:)
    deli(:)=((dpi(:)+mfppi(:))**3-                                             &
             SQRT((dpi(:)*dpi(:)+mfppi(:)*mfppi(:))**3))/                      &
             (3.0*dpi(:)*mfppi(:))-dpi(:)
    !
    dpj(:)=2.0*rj(:)
    velj(:)=SQRT(term1*t(:)/(rhoj(:)*vj(:)))
    knj(:)=mfpa(:)/rj(:)
    ccj(:)=1.0+knj(:)*(1.257+0.4*EXP(-1.1/knj(:)))
    ! .. n.b. use 1.257,0.4,1.1 as in Seinfeld & Pandis pg. 465
    ! .. following Allen & Raabe (1982) reanalysis of Millikan (1923)
    dcoefj(:)=term2*ccj(:)*t(:)/(rj(:)*dvisc(:))
    mfppj(:)=term3*dcoefj(:)/velj(:)
    delj(:)=((dpj(:)+mfppj(:))**3-                                             &
             SQRT((dpj(:)*dpj(:)+mfppj(:)*mfppj(:))**3))/                      &
             (3.0*dpj(:)*mfppj(:))-dpj(:)
    !
    rtot(:)=ri(:)+rj(:)
    dtot(:)=dcoefi(:)+dcoefj(:)
    termv1(:)=                                                                 &
     rtot(:)/(rtot(:)+SQRT(deli(:)*deli(:)+delj(:)*delj(:)))
    termv2(:)=                                                                 &
     4.0*dtot(:)/SQRT(veli(:)*veli(:)+velj(:)*velj(:))/rtot(:)
    kij(:)=term4*rtot(:)*dtot(:)/(termv1(:)+termv2(:))
  END WHERE
END IF
!
!  .. For ICOAG=2, calculates KIJ following the M7 method as
!  .. described in Stier et al (2005) in ACP (see header above)
IF (icoag == 2) THEN
  WHERE (mask(:))
    rmid(:)=0.5*(ri(:)+rj(:))
    vmid(:)=(pi/0.75)*rmid(:)**3
    rhomid(:)=0.5*(rhoi(:)+rhoj(:))
    vel(:)=SQRT(term1*t(:)/(rhomid(:)*vmid(:)))
    kn(:)=mfpa(:)/rmid(:)
    cc(:)=1.0+kn(:)*(1.257+0.4*EXP(-1.1/kn(:)))
    ! .. n.b. use 1.257,0.4,1.1 as in Seinfeld & Pandis pg. 465
    ! .. following Allen & Raabe (1982) reanalysis of Millikan (1923)
    dcoef(:)=term2*cc(:)*t(:)/(rmid(:)*dvisc(:))
    termv1(:)=4.0*dcoef(:)/(vel(:)*rmid(:))
    mfpp(:)=2.0*dcoef(:)/(pi*vel(:))
    termv2(:)=rmid(:)/(rmid(:)+mfpp(:))
    kij(:)=term5*rmid(:)*dcoef(:)/(termv1(:)+termv2(:))
  END WHERE
END IF
!
!  .. For ICOAG=3, KIJ is calculated as in the original UM sulphate
!  .. aerosol scheme where the Cunningham slip correction was assumed
!  .. the same for each mode and the mean free path length for
!  .. particles in each mode assumed to that of air
IF (icoag == 3) THEN
  WHERE (mask(:))
    termv3(:)=(2.0e6*boltzmann*t(:)/3.0/dvisc(:))
    termv4(:)=1.591*ukca_mfp_ref*(1.0/ri(:)/ri(:)+1.0/rj(:)/rj(:))
    kij(:)=termv3(:)*(ri(:)+rj(:))*(1.0/ri(:)+1.0/rj(:)+termv4(:))
  END WHERE
END IF
!
!  .. For ICOAG=4, KIJ is calculated as in the original UM sulphate
!  .. aerosol scheme but values of CCI,CCJ,MFPPI,MFPPJ are computed
!  .. rather than just setting MFPPI=MFPPJ=MFPA and CCI=CCJ=1.591
IF (icoag == 4) THEN
  WHERE (mask(:))
    kni(:)=mfpa(:)/ri(:)
    cci(:)=1.0+kni(:)*(1.257+0.4*EXP(-1.1/kni(:)))
    ! .. n.b. use 1.257,0.4,1.1 as in Seinfeld & Pandis pg. 465
    ! .. following Allen & Raabe (1982) reanalysis of Millikan (1923)
    ! .. UM originally used 1.249, 0.42 and -0.87 as in
    ! .. Jacobson page 445 [following Kasten (1968) reanalysis of Millikan]
    knj(:)=mfpa(:)/rj(:)
    ccj(:)=1.0+knj(:)*(1.257+0.4*EXP(-1.1/knj(:)))
    ! .. n.b. use 1.257,0.4,1.1 as in Seinfeld & Pandis pg. 465
    ! .. following Allen & Raabe (1982) reanalysis of Millikan (1923)
    ! .. UM originally used 1.249, 0.42 and -0.87 as in
    ! .. Jacobson page 445 [following Kasten (1968) reanalysis of Millikan]
    termv3(:)=(2.0e6*boltzmann*t(:)/3.0/dvisc(:))
    termv4(:)=                                                                 &
     mfppi(:)*cci(:)/ri(:)/ri(:)+mfppj(:)*ccj(:)/rj(:)/rj(:)
    kij(:)=termv3(:)*(ri(:)+rj(:))*(1.0/ri(:)+1.0/rj(:)+termv4(:))
  END WHERE
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_coag_coff_v
END MODULE ukca_coag_coff_v_mod
