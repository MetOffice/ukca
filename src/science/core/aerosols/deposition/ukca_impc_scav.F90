!*****************************COPYRIGHT*******************************
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
!    Subroutine to calculate impaction scavenging of aerosols
!    by falling raindrops.
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
!    Language:  Fortran
!
! ######################################################################
!
! Subroutine Interface:
MODULE ukca_impc_scav_mod

IMPLICIT NONE

! Variables for impaction scavenging (as in Pringle, 2006 PhD thesis)
INTEGER, PARAMETER, PRIVATE  :: ncoll=20 ! # of columns in LUT (aer. bins)
INTEGER, PARAMETER, PRIVATE  :: nrow =19 ! # of rows in LUT (raindrop bins)

! raindrop bins
REAL, PARAMETER, PRIVATE :: raddrop(nrow)  =                                   &
         [   1.0, 1.587,  2.52,   4.0,  6.35, 10.08,                           &
             16.0,  25.4, 40.32,  64.0, 101.6, 161.3,                          &
            256.0, 406.4, 645.1,1024.0,1625.0,2580.0,4096.0]

REAL, PRIVATE, SAVE :: colleff4(ncoll,nrow) ! collision efficiency

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_IMPC_SCAV_MOD'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE ukca_mode_imscavcoff(verbose)

! Set values of collision efficiencies

USE umPrintMgr,         ONLY: umPrint, umMessage
USE yomhook,            ONLY: lhook, dr_hook
USE parkind1,           ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: verbose ! flag to indicate level of verbosity

INTEGER :: i, j  ! loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_IMSCAVCOFF'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up collision efficiencies for all aerosol bins and raindrop bins
colleff4( 1,1:nrow)=[ 0.522e+05,0.139e+05,0.328e+04,0.775e+03,                 &
                       0.183e+03,0.432e+02,0.102e+02,0.291e+01,                &
                       0.108e+01,0.439e+00,0.201e+00,0.110e+00,                &
                       0.633e-01,0.366e-01,0.815e-02,0.168e-02,                &
                       0.394e-03,0.591e-04,0.132e-04]

colleff4( 2,1:nrow)=[ 0.126e+05,0.373e+04,0.985e+03,0.260e+03,                 &
                       0.687e+02,0.182e+02,0.480e+01,0.150e+01,                &
                       0.536e+00,0.224e+00,0.107e+00,0.608e-01,                &
                       0.364e-01,0.236e-01,0.731e-02,0.167e-02,                &
                       0.395e-03,0.592e-04,0.133e-04]

colleff4( 3,1:nrow)=[ 0.445e+04,0.139e+04,0.390e+03,0.110e+03,                 &
                       0.308e+02,0.864e+01,0.243e+01,0.783e+00,                &
                       0.270e+00,0.117e+00,0.564e-01,0.325e-01,                &
                       0.203e-01,0.136e-01,0.595e-02,0.166e-02,                &
                       0.396e-03,0.594e-04,0.133e-04]

colleff4( 4,1:nrow)=[ 0.259e+04,0.810e+03,0.227e+03,0.639e+02,                 &
                       0.179e+02,0.503e+01,0.141e+01,0.440e+00,                &
                       0.146e+00,0.662e-01,0.309e-01,0.177e-01,                &
                       0.115e-01,0.728e-02,0.446e-02,0.164e-02,                &
                       0.399e-03,0.597e-04,0.134e-04]

colleff4( 5,1:nrow)=[ 0.196e+04,0.602e+03,0.166e+03,0.457e+02,                 &
                       0.126e+02,0.346e+01,0.953e+00,0.280e+00,                &
                       0.927e-01,0.410e-01,0.183e-01,0.990e-02,                &
                       0.655e-02,0.444e-02,0.389e-02,0.164e-02,                &
                       0.402e-03,0.604e-04,0.135e-04]

colleff4( 6,1:nrow)=[ 0.192e+04,0.569e+03,0.151e+03,0.401e+02,                 &
                       0.106e+02,0.282e+01,0.749e+00,0.206e+00,                &
                       0.689e-01,0.280e-01,0.119e-01,0.580e-02,                &
                       0.384e-02,0.304e-02,0.393e-02,0.165e-02,                &
                       0.406e-03,0.612e-04,0.138e-04]

colleff4( 7,1:nrow)=[ 0.208e+04,0.604e+03,0.156e+03,0.405e+02,                 &
                       0.105e+02,0.271e+01,0.703e+00,0.185e+00,                &
                       0.572e-01,0.217e-01,0.903e-02,0.462e-02,                &
                       0.280e-02,0.208e-02,0.398e-02,0.168e-02,                &
                       0.415e-03,0.628e-04,0.142e-04]

colleff4( 8,1:nrow)=[ 0.233e+04,0.636e+03,0.162e+03,0.410e+02,                 &
                       0.104e+02,0.264e+01,0.671e+00,0.174e+00,                &
                       0.513e-01,0.184e-01,0.747e-02,0.384e-02,                &
                       0.222e-02,0.161e-02,0.408e-02,0.173e-02,                &
                       0.430e-03,0.657e-04,0.149e-04]

colleff4( 9,1:nrow)=[ 0.235e+04,0.659e+03,0.165e+03,0.412e+02,                 &
                       0.103e+02,0.257e+01,0.643e+00,0.168e+00,                &
                       0.490e-01,0.168e-01,0.661e-02,0.326e-02,                &
                       0.188e-02,0.140e-02,0.422e-02,0.180e-02,                &
                       0.452e-03,0.698e-04,0.160e-04]

colleff4(10,1:nrow)=[ 0.165e+04,0.457e+03,0.112e+03,0.277e+02,                 &
                       0.680e+01,0.167e+01,0.412e+00,0.106e+00,                &
                       0.304e-01,0.999e-02,0.386e-02,0.186e-02,                &
                       0.124e-02,0.140e-02,0.447e-02,0.193e-02,                &
                       0.491e-03,0.771e-04,0.179e-04]

colleff4(11,1:nrow)=[ 0.899e+03,0.246e+03,0.597e+02,0.145e+02,                 &
                       0.352e+01,0.856e+00,0.208e+00,0.524e-01,                &
                       0.145e-01,0.466e-02,0.179e-02,0.860e-03,                &
                       0.719e-03,0.165e-02,0.486e-02,0.213e-02,                &
                       0.554e-03,0.891e-04,0.211e-04]

colleff4(12,1:nrow)=[ 0.117e+04,0.326e+03,0.807e+02,0.200e+02,                 &
                       0.496e+01,0.123e+01,0.305e+00,0.777e-01,                &
                       0.219e-01,0.720e-02,0.281e-02,0.137e-02,                &
                       0.941e-03,0.330e-02,0.563e-02,0.255e-02,                &
                       0.686e-03,0.116e-03,0.283e-04]

colleff4(13,1:nrow)=[ 0.130e+04,0.371e+03,0.938e+02,0.237e+02,                 &
                       0.601e+01,0.152e+01,0.385e+00,0.979e-01,                &
                       0.276e-01,0.926e-02,0.364e-02,0.173e-02,                &
                       0.101e-02,0.406e-02,0.694e-02,0.327e-02,                &
                       0.930e-03,0.167e-03,0.429e-04]

colleff4(14,1:nrow)=[ 0.118e+04,0.333e+03,0.842e+02,0.213e+02,                 &
                       0.537e+01,0.136e+01,0.342e+00,0.876e-01,                &
                       0.250e-01,0.841e-02,0.330e-02,0.153e-02,                &
                       0.801e-03,0.260e-02,0.973e-02,0.490e-02,                &
                       0.152e-02,0.303e-03,0.842e-04]

colleff4(15,1:nrow)=[ 0.774e+03,0.223e+03,0.572e+02,0.147e+02,                 &
                       0.378e+01,0.970e+00,0.249e+00,0.636e-01,                &
                       0.180e-01,0.606e-02,0.238e-02,0.110e-02,                &
                       0.658e-03,0.107e-02,0.167e-01,0.940e-02,                &
                       0.335e-02,0.791e-03,0.249e-03]

colleff4(16,1:nrow)=[ 0.372e+01,0.177e+01,0.781e+00,0.345e+00,                 &
                       0.153e+00,0.675e-01,0.299e-01,0.130e-01,                &
                       0.624e-02,0.346e-02,0.179e-02,0.142e-02,                &
                       0.164e-02,0.401e-02,0.413e-01,0.277e-01,                &
                       0.124e-01,0.389e-02,0.152e-02]

colleff4(17,1:nrow)=[ 0.234e-18,0.108e-16,0.705e-15,0.462e-13,                 &
                       0.302e-11,0.198e-09,0.129e-07,0.844e-06,                &
                       0.568e-04,0.386e-02,0.286e-01,0.448e-01,                &
                       0.569e-01,0.859e-01,0.183e+00,0.165e+00,                &
                       0.108e+00,0.536e-01,0.295e-01]

colleff4(18,1:nrow)=[ 0.902e-37,0.794e-33,0.160e-28,0.324e-24,                 &
                       0.655e-20,0.132e-15,0.267e-11,0.540e-07,                &
                       0.816e-03,0.482e-01,0.245e+00,0.372e+00,                &
                       0.436e+00,0.473e+00,0.493e+00,0.497e+00,                &
                       0.414e+00,0.299e+00,0.225e+00]

colleff4(19,1:nrow)=[ 0.275e-30,0.136e-26,0.146e-22,0.156e-18,                 &
                       0.167e-14,0.178e-10,0.191e-06,0.202e-02,                &
                       0.203e+00,0.427e+00,0.586e+00,0.669e+00,                &
                       0.708e+00,0.730e+00,0.746e+00,0.738e+00,                &
                       0.679e+00,0.588e+00,0.520e+00]

colleff4(20,1:nrow)=[ 0.136e-33,0.238e-29,0.102e-24,0.436e-20,                 &
                       0.186e-15,0.797e-11,0.341e-06,0.143e-01,                &
                       0.722e+00,0.805e+00,0.869e+00,0.902e+00,                &
                       0.915e+00,0.932e+00,0.104e+01,0.927e+00,                &
                       0.904e+00,0.871e+00,0.842e+00]

IF (verbose >= 2) THEN
  WRITE(umMessage,'(A50,2I5)') 'Set up impaction scavenging coeffs,'//         &
                        'NROW,NCOLL=',nrow,ncoll
  CALL umPrint(umMessage,src=RoutineName)
  DO i = 1,nrow
    WRITE(umMessage,'(A15,I5,E10.2)') 'I,RADDROP(I)=',i,raddrop(i)
    CALL umPrint(umMessage,src=RoutineName)
    DO j = 1,ncoll
      WRITE(umMessage,'(A,2I5,E12.4)')'I,J,COLLEFF4(J,I)=',i,j,colleff4(j,i)
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_mode_imscavcoff
! ----------------------------------------------------------------------
SUBROUTINE ukca_impc_scav(nbox,nbudaer,nd,md,                                  &
 crain,drain,csnow,dsnow,wetdp,dtc,l_dust_slinn_impc_scav,bud_aer_mas)
! ----------------------------------------------------------------------
!
!     Subroutine to calculate impaction scavenging of aerosols
!     by falling raindrops and snow.
!
!     Uses empirical expression for the terminal velocity of raindrops
!     of Easter & Hales (1984)
!
!     Uses look-up table of aerosol-raindrop collision efficiencies
!     provided by Yan Yin, University of Aberystwyth. These were
!     originally generated for a specific NROW=19-bin raindrop radius
!     grid colliding with the NCOLL=20-bin aerosol size grid currently
!     used in GLOMAP.
!
!     Currently the routine is only set for a NRBINS=7-bin raindrop
!     size grid consisting of the 4th,6th,8th,10th,12th,14th,16th bin
!     centres of the original 19-bin raindrop size grid. These
!     correspond to raindrop radii of 4.0, 10.08, 25.4, 64.0, 161.3,
!     406.4 and 1024.0 microns and a geometric scaling factor of 2.52.
!
!     The raindrop radii and corresponding aerosol collision
!     efficiencies have already been read in from file and are
!     passed into IMPC_SCAV explicitly as the arrays RADDROP(NROW)
!     and COLLEFF4(NCOLL,NROW)
!
!     Raindrop size distribution assumed is Marshall-Palmer distribution
!     as modified by Sekhon & Srivastava (1971) to take into account
!     rainfall intensity.
!
!     MP have dN/dDp = n_0 * exp (-Psi*Dp) with Psi=4.1*p0^-0.21 mm^-1
!
!       where p0 is precipitation rate (mm/hour)
!             Dp is particle diameter (mm)
!             n_0 is a constant = 8000 drops/m3
!
!     SS modified this to set n_0 to also vary with precipitation rate:
!
!             n_0=7000*p0^0.37 m^-3 mm^-1
!         and Psi= 3.8*p0^-0.14 mm^-1
!
!     Scheme developed using the BCS relation of Slinn 1983, scavenging
!     coeffs supplied by Yan Yin, from Flossmann et al (1985)
!
!     This scheme was developed by Kirsty Pringle as part of PhD thesis.
!
!     The impaction scavenging by snowfall is based on the simple
!     parameterization of Wang et al. (2011)
!                             k = aP^b ,
!     where P is the total snowfall rate and the coefficients a and b
!     are dependent on aerosol size and are taken from Wang et al. (2011)
!     and Feng (2009).
!
!     References
!     ----------
!     * Flossmann, A. I.; Hall, W. D. & Pruppacher, H. R. (1985)
!       "A theoretical study of the wet removal of atmospheric pollutants:
!       Part 1: The redistribution of aerosol particles captured through
!       nucleation and impaction scavenging by growing cloud drops."
!       J. Atmos. Sci., vol. 42, pp. 582--606.
!
!     * Andronache, C. (2003),
!       "Estimated variability of below-cloud aerosol removal by rainfall
!       for observed aerosol size distribution"
!       Atmos. Chem. Phys., vol. 3, pp. 131--143.
!
!     * Gong, S.-L.; Barry, L. A. & Blanchet J.-P. (1997)
!       "Modelling sea-salt aerosols in the atmosphere,
!        1: Model development!"
!        J. Geophys. Res., vol. 102, no. D3, pp 3,805-3,818
!
!     * Easter, R.C. & Hales, J. M. (1983)
!       "Precipitation scavenging, dry deposition and resuspension,"
!       Chapter Interpretation of the OSCAR data for reactive gas
!       scavenging, pp. 649-662.
!
!     * Sekhon & Srivastava (1971)
!       "Doppler observations of drop size distributions in a thunderstorm"
!       J. Atmos. Sci., vol. 28, pp. 983--984.
!
!     * Feng, J. (2009),
!       "A size-resolved model for below-cloud scavenging of aerosols
!        by snowfall." JGR, 114(D8), doi: 10.1029/ 2008JD011012 .
!
!     * Wang, Q. et al. (2011),
!       "Sources of carbonaceous aerosols and deposited black carbon in
!        the Arctic in winter-spring: implications for radiative forcing."
!        ACP, 11, 12453-12473.
!
!
!     Parameters
!     ----------
!
!     Inputs
!     ------
!     NBOX        : Number of grid boxes
!     nbudaer     : Number of aerosol budget fields
!     ND          : Aerosol ptcl number density for size bin (cm^-3)
!     MD          : Avg cpt mass of aerosol ptcl in size bin (particle^-1)
!     CRAIN       : Convective rain rate array (kgm^-2s^-1)
!     DRAIN       : Dynamic rain rate array (kgm^-2s^-1)
!     csnow       : Convective snowfall rate array (kgm^-2s^-1)
!     dsnow       : Dynamic snowfall rate array (kgm^-2s^-1)
!     WETDP       : Wet diameter corresponding to DRYDP (m)
!     DTC         : Time step of process (s)
!     L_DUST_SLINN_IMPC_SCAV : Flag for new impaction scavenging
!                              scheme for dust
!
!     Outputs
!     -------
!     ND          : new aerosol number conc (cm^-3)
!     BUD_AER_MAS : Updated aerosol budgets
!
!     Local variables
!     ---------------
!     TOTRAIN     : Total combined rain rate array (conv + dyn) (mm/hr)
!     RNDPDIAM_cm : Diameter of raindrop (cm)
!     RNDPDIAM_mm : Diameter of raindrop (mm)
!     RNDPDIAM_m  : Diameter of raindrop (m)
!     VELDR_cms   : Terminal velocity of raindrop (cm/s)
!     VELDR_ms    : Terminal velocity of raindrop (m/s)
!     SCAV        : Holds scavenging coefficient for each of the 7 rain bins
!     SCAVCOEFF_COUNT : Holds sum of calculated scavenging coeffs
!                        over all rain bins for each aerosol bin
!     SCAVCOEFF   : Holds final calculated total scavenging coeff for each
!                   aerosol bin summed over all 7 rain bins
!     COUNTCOLL   : Holds index of column (aerosol size) in aerosol-raindrop
!                   collision l-u table for calculated aerosol bin wet radius
!     NRAINMAX    : dN/dD_p [=n0 in Seinfeld & Pandis pg. 832] (m^-3 mm-1)
!     NDRAIN      : dN/d(log D_p) * delta(log(D_p)) [D_p=rndrp diam] (m^-3)
!                   n.b. delta(log(D_p))=ln(2.52)=0.924, where 2.52 is the
!                   geometric scaling factor for the raindrop size grid
!     INTERZZ     : Holds (pi/4)*RNDPDIAM_m^2*VELDR_ms*NDRAIN
!     INTERC      : Holds Psi*D_p in Marshall-Palmer distribution
!     INTERB      : Holds term in determination of column in LUT.
!     R1          : mid-point of first aerosol particle size bin (microns)
!     FACTOR      : geometric scaling factor for aerosol particle size grid
!     FC          : Fraction of grid box over which convective precip occurs.
!     FD          : Fraction of grid box over which dynamic precip occurs.
!     NRBINS      : No. of raindrop size bins used in GLOMAP raindrop spectrum
!     asnow       : Size dependent parameter used to calc rscav_snow
!     bsnow       : Used to calculate rscav_snow
!     rscav_snow  : Scavenging coefficient (hr-1) for snow droplets
!                   (RSCAV=A*SR^B, where A=asnow, B=bsnow,
!                   SR=snowfall rate (mm hr-1))
!                   uses Wang et al 2011, ACP, 11 12453-12473

!     allrain     : conv & dyn rain rates in mm h-1, in precip area
!     allsnow     : conv & dyn snow rates in mm h-1, in precip area
!     allscavcoeff: single array holding scav coeffs for each
!                   aerosol bin summed over all 7 rain bins, for each rain type
!     allrscav_snow: Scavenging coefficient for snow for each precip type
!     allfrac     : fraction of gridbox where precip occurs
!     DELN1       : Change in aerosol no. conc. by impaction scav by CRAIN
!     DELN2       : Change in aerosol no. conc. by impaction scav by DRAIN
!     deln3       : Change in aerosol no. conc. by impaction scav by csnow
!     deln4       : Change in aerosol no. conc. by impaction scav by dsnow
!     DELN        : Change in aerosol no. conc. by impaction scav (total)
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES      : Number of modes set
!     NCP         : Number of components set
!     MODE        : Logical variable defining which modes are set.
!     COMPONENT   : Logical variable defining which components are set.
!     NUM_EPS     : Value of NEWN below which don't recalculate MD
!                                                   or carry out process
!     CP_SU       : Index of component in which SO4    cpt is stored
!     CP_BC       : Index of component in which BC     cpt is stored
!     CP_OC       : Index of component in which 1st OC cpt is stored
!     CP_CL       : Index of component in which NaCl   cpt is stored
!     CP_DU       : Index of component in which dust   cpt is stored
!     CP_SO       : Index of component in which 2nd OC cpt is stored
!
!     Inputted by module UKCA_SETUP_INDICES
!     -------------------------------------
!     Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
USE ukca_constants,     ONLY: pi

USE ukca_mode_setup,    ONLY:                                                  &
    cp_su, cp_bc, cp_oc, cp_cl, cp_so, cp_du, nmodes,                          &
    mode_nuc_sol, mode_ait_sol, mode_acc_sol, mode_cor_sol,                    &
    mode_ait_insol, mode_acc_insol, mode_cor_insol, mode_sup_insol,            &
    cp_no3, cp_nh4, cp_nn

USE ukca_setup_indices, ONLY:                                                  &
    nmasimscbcaccsol, nmasimscbcaitins, nmasimscbcaitsol,                      &
    nmasimscbccorsol, nmasimscduaccins, nmasimscduaccsol,                      &
    nmasimscducorins, nmasimscducorsol, nmasimscocaccsol,                      &
    nmasimscocaitins, nmasimscocaitsol, nmasimscoccorsol,                      &
    nmasimscocnucsol, nmasimscsoaccsol, nmasimscsoaitsol,                      &
    nmasimscsocorsol, nmasimscsonucsol, nmasimscssaccsol,                      &
    nmasimscsscorsol, nmasimscsuaccsol, nmasimscsuaitsol,                      &
    nmasimscsucorsol, nmasimscsunucsol,                                        &
    nmasimscntaitsol, nmasimscntaccsol, nmasimscntcorsol,                      &
    nmasimscnhaitsol, nmasimscnhaccsol, nmasimscnhcorsol,                      &
    nmasimscnnaccsol, nmasimscnncorsol, nmasimscdusupins

USE yomhook,            ONLY: lhook, dr_hook
USE parkind1,           ONLY: jprb, jpim
USE ukca_config_specification_mod, ONLY: glomap_config, glomap_variables

IMPLICIT NONE

! .. Subroutine interface
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: nbudaer
REAL, INTENT(IN)    :: md(nbox,nmodes,glomap_variables%ncp)
REAL, INTENT(IN)    :: wetdp(nbox,nmodes)
REAL, INTENT(IN)    :: dtc
REAL, INTENT(IN)    :: crain(nbox)
REAL, INTENT(IN)    :: drain(nbox)
REAL, INTENT(IN)    :: csnow(nbox)
REAL, INTENT(IN)    :: dsnow(nbox)
LOGICAL, INTENT(IN)  :: l_dust_slinn_impc_scav
REAL, INTENT(IN OUT) :: nd(nbox,nmodes)
REAL, INTENT(IN OUT) :: bud_aer_mas(nbox,0:nbudaer)

! .. Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: num_eps(:)

INTEGER :: imode
INTEGER :: icp
INTEGER :: jl
INTEGER :: jvr
INTEGER :: irow
INTEGER :: icoll
INTEGER :: iprecip
INTEGER :: countcoll(nbox,nmodes)
INTEGER :: topmode
INTEGER, PARAMETER :: nrbins=7
REAL    :: r1
REAL    :: factor
REAL    :: interzz(nbox,nrbins)
REAL    :: scav(nbox,nmodes,nrbins)
REAL    :: scavcoeff(nbox,nmodes)
REAL    :: scavcoeff_count(nbox,nmodes,nrbins)
REAL    :: interc(nbox,nrbins)
REAL    :: interb
REAL    :: VELDR_ms(nrbins)
REAL    :: VELDR_cms(nrbins)
REAL    :: RNDPDIAM_cm(nrbins)
REAL    :: RNDPDIAM_mm(nrbins)
REAL    :: RNDPDIAM_m(nrbins)
REAL    :: ndrain(nbox,nrbins)
REAL    :: nrainmax(nbox)
REAL    :: rscav_snow(nbox,nmodes)
REAL    :: rscav_snowh(nbox,nmodes)
REAL    :: asnow(nmodes)
REAL, PARAMETER :: bsnow=0.96
REAL    :: totrain(nbox)
REAL    :: totsnow(nbox)
REAL    :: totppn(nbox)
REAL    :: allrain(2,nbox)
REAL    :: allsnow(2,nbox)
REAL    :: allscavcoeff(2,nbox,nmodes)
REAL    :: allrscav_snow(2,nbox,nmodes)
REAL    :: allfrac(2)
REAL    :: deln
REAL    :: deln1
REAL    :: deln2
REAL    :: deln3
REAL    :: deln4
REAL    :: dm(glomap_variables%ncp)
REAL, PARAMETER :: fc=0.3
REAL, PARAMETER :: fd=1.0

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_IMPC_SCAV'

allfrac(1)=fc
allfrac(2)=fd

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
num_eps     => glomap_variables%num_eps

! Prior to UM10.4 this routine contained two bugs.
! A fixed version of the code is enabled by setting l_fix_ukca_impscav true.
!
!The bugs in question are:
!1)  Convective fraction (FC) was used in the calculation of removal rate, but
!not in the calculation of rain properties used for the scavenging coefficients.
!2) Scavenging coefficients were calculated using total precip rates,
!(convective plus dynamic) but then subsequently used to calculate removal by
!convective and dynamic precip separately, effectively "double counting".
!

! Below cloud scavenging for dust by rain for the accumulation/coarse
! insoluble modes can now be dealt with in a separate module. Snow is
! still accounted for in this routine though
IF (l_dust_slinn_impc_scav) THEN
  topmode = mode_ait_insol
ELSE
  topmode = nmodes
END IF

IF (glomap_config%l_fix_ukca_impscav) THEN

  !.........................................................................
  ! Fixed version of code
  !

  ! .. Here convert rain rate from kg/m2/s to mm/hr
  ! .. Combine the two rain rates (con (1) and dyn (2) )
  DO jl=1,nbox
    totrain(jl)=(crain(jl)*3600.0)+(drain(jl)*3600.0)
    totsnow(jl)=(csnow(jl)*3600.0)+(dsnow(jl)*3600.0)

    totppn(jl)=totrain(jl)+totsnow(jl)

    allrain(1,jl)=crain(jl)*3600.0/allfrac(1)
    allrain(2,jl)=drain(jl)*3600.0/allfrac(2)
    allsnow(1,jl)=csnow(jl)*3600.0/allfrac(1)
    allsnow(2,jl)=dsnow(jl)*3600.0/allfrac(2)

  END DO

  ! .. In LUT there are 19 Rows but we only have NRBINS=7 raindrop bins,
  ! .. so use every 2nd row starting from row 4
  ! .. Only consider raindrops radius >4 and <1024 um
  ! .. NB/ HAVE TO CHANGE THIS IF NRBINS IS ALTERED

  jvr=0
  DO irow=4,16,2 ! pick out 7 sizes corresponding to 7 raindrop bins
    jvr=jvr+1
    RNDPDIAM_cm(jvr)=raddrop(irow)*2.0/10000.0
    RNDPDIAM_mm(jvr)=RNDPDIAM_cm(jvr)*10.0
    RNDPDIAM_m(jvr)=RNDPDIAM_mm(jvr)/1000.0
    ! RADDROP contains raindrop bin radii in microns
    ! RNDPDIAM_cm contains raindrop bin diameters in cm
    ! RNDPDIAM_mm & RNDPDIAM_m contain rndrp bin diams in mm & m respect.
  END DO

  ! .. Initialise Arrays
  scavcoeff_count(:,:,:)=0.0
  scavcoeff(:,:)=0.0
  allscavcoeff(:,:,:)=0.0
  allrscav_snow(:,:,:)=0.0
  countcoll(:,:)=0
  interzz(:,:) = 0.0

  DO imode=1,topmode
    IF (mode(imode)) THEN
      DO jl=1,nbox
        IF (totrain(jl) > 0.0) THEN
          !From aerosol wet radius find which column of LU table is needed
          IF (0.5*wetdp(jl,imode) < 0.001e-6) THEN
            countcoll(jl,imode)=1
          ELSE IF (0.5*wetdp(jl,imode) > 9.0e-6) THEN
            countcoll(jl,imode)=10
          ELSE
            factor=0.457
            r1=0.001
            interb=((0.5*wetdp(jl,imode)*1.0e6)/r1)
            ! 1.0E6 in above to convert wet radius to um
            countcoll(jl,imode)=INT((LOG(interb)/factor)+1.0)
          END IF
        END IF
      END DO
    END IF
  END DO

  ! .. Loop to calculate the number of rain drops at each rdrop bin-centre
  ! .. following Marshall Palmer distribution

  DO jvr=1,nrbins
    ! Empirical relationship from Easter and Hales (1984) to calulate
    ! terminal velocity of raindrop in cm/s
    IF (RNDPDIAM_cm(jvr) <= 0.10) THEN
      VELDR_cms(jvr)=4055.0*RNDPDIAM_cm(jvr)
    ELSE IF (RNDPDIAM_cm(jvr) > 0.10) THEN
      VELDR_cms(jvr)=13000.0*(RNDPDIAM_m(jvr)**0.5)
    END IF
    ! .. Convert raindrop velocity from cm/s to m/s
    VELDR_ms(jvr)=VELDR_cms(jvr)/100.0
  END DO

  ! Loop over the 2 rain types
  ! Rain properties are calculated using the precip rates in the part of the
  ! gridbox where is raining

  DO iprecip=1,2
    DO jl=1,nbox
      IF (allrain(iprecip,jl) > 0.0) THEN
        DO jvr=1,nrbins

          ! Marshall Palmer but with sophistication of NRAINMAX calculated
          ! according to Sekhon & Srivastava (1971) (Seinfeld & Pandis pg 832)

          nrainmax(jl)=7000.0*allrain(iprecip,jl)**0.37
          ! NRAINMAX is n_0 in m^-3 per mm
          interc(jl,jvr)=3.8*allrain(iprecip,jl)**(-0.14)*RNDPDIAM_mm(jvr)
          ! INTERC is Psi*D_p in MP distribution
          ndrain(jl,jvr)=nrainmax(jl)*RNDPDIAM_mm(jvr)                         &
                    *0.924*EXP(-interc(jl,jvr))
          ! Includes various conversions as Nd given in terms of mm-1
          ! => *diameter(mm)*width of bin
          ! Geometric scaling factor = 2.52 (change if NRBINS is changed)
          ! (ln(2.52)=0.924)
          !
          ! Leave NDRAIN in m-3 for use in SCAV calculation below
          ! Convert Raindrop size from diameter in mm to diameter in m
          interzz(jl,jvr)=                                                     &
            (pi/4.0)*((RNDPDIAM_m(jvr)*RNDPDIAM_m(jvr))                        &
            *VELDR_ms(jvr)*ndrain(jl,jvr))

        END DO ! over NRBINS
      END IF ! TOTRAIN>0
    END DO ! over NBOX

    ! .. Calculate scavenging coefficients
    DO jl=1,nbox
      IF (allrain(iprecip,jl) > 0.0) THEN
        DO imode=1,topmode
          IF (mode(imode)) THEN
            IF (nd(jl,imode) > num_eps(imode)) THEN
              ! Loop over raindrop bins
              DO jvr=1,nrbins
                irow=(jvr*2)+2 ! hard-wired for 7 raindrop size bins
                !COUNTCOLL contains index of relevant column of the LU table for
                !aerosol-rndrp coll efficiency for each aerosol bin wet radius
                icoll=countcoll(jl,imode)
                IF (countcoll(jl,imode) == 0) icoll=1
                scav(jl,imode,jvr)=colleff4(icoll,irow)*interzz(jl,jvr)

                !For each aerosol bin, sum scav coeff over all rain bins
                IF (jvr == 1) THEN
                  scavcoeff_count(jl,imode,jvr)=scav(jl,imode,jvr)
                END IF
                IF (jvr > 1) THEN
                  scavcoeff_count(jl,imode,jvr)=                               &
                  scavcoeff_count(jl,imode,jvr-1)+scav(jl,imode,jvr)
                END IF
                IF (jvr == nrbins) THEN
                  allscavcoeff(iprecip,jl,imode)=scavcoeff_count(jl,imode,jvr)
                END IF

              END DO ! LOOP OVER RAINDROP BINS
            END IF ! IF ND>ND_EPS
          END IF ! IF MODE PRESENT
        END DO ! LOOP OVER MODES
      END IF ! IF TOTRAIN>0

      IF (totsnow(jl) > 0.0) THEN

        ! Tuned values of ASNOW: lower than in used Feng (2009), but
        ! within theoretical and observational uncertainties
        ! (Fig 8 of Feng 2009)

        asnow(:) = [ 0.014, 0.014, 0.014, 0.10, 0.014, 0.014, 0.10, 0.50 ]

        DO imode=1,nmodes
          IF (mode(imode)) THEN
            IF (nd(jl,imode) > num_eps(imode)) THEN
              ! Calculate scavenging coefficients for snow
              ! RSCAV = aR^b  (UNIT: hr-1)
              allrscav_snow(iprecip,jl,imode) =                                &
               ( asnow(imode) * allsnow(iprecip,jl) ** bsnow )/3600.0
              ! Above converts to unit of s-1
            END IF ! IF ND>ND_EPS
          END IF ! IF MODE PRESENT
        END DO ! LOOP OVER MODES

      END IF ! IF totsnow>0
    END DO ! LOOP OVER BOXES
  END DO ! loop over precip types

  ! .. Below calculates the rate of removal
  !
  ! .. Calculate removal (from each bin) following 1st order rate loss
  ! .. Apply BCS only to the portion of the gridbox where rain is occuring

  DO jl=1,nbox
    IF (totppn(jl) > 0.0) THEN
      DO imode=1,nmodes
        IF (mode(imode)) THEN
          IF (nd(jl,imode) > num_eps(imode)) THEN
            IF (crain(jl) > 0.0) THEN
              deln1 = allfrac(1)*nd(jl,imode) *                                &
                (1.0-EXP(-allscavcoeff(1,jl,imode)*dtc))
              nd(jl,imode)=nd(jl,imode)-deln1
              IF (nd(jl,imode) < 0.0) nd(jl,imode)=0.0
            ELSE
              deln1=0.0
            END IF
            IF (drain(jl) > 0.0) THEN
              deln2 = allfrac(2)*nd(jl,imode) *                                &
                (1.0-EXP(-allscavcoeff(2,jl,imode)*dtc))
              nd(jl,imode)=nd(jl,imode)-deln2
              IF (nd(jl,imode) < 0.0) nd(jl,imode)=0.0
            ELSE
              deln2=0.0
            END IF

            IF (csnow(jl) > 0.0) THEN
              deln3 = allfrac(1)*nd(jl,imode) *                                &
                (1.0-EXP(-allrscav_snow(1,jl,imode)*dtc))
              nd(jl,imode)=nd(jl,imode)-deln3
              IF (nd(jl,imode) < 0.0) nd(jl,imode)=0.0
            ELSE
              deln3=0.0
            END IF
            IF (dsnow(jl) > 0.0) THEN
              deln4 = allfrac(2)*nd(jl,imode) *                                &
                (1.0-EXP(-allrscav_snow(2,jl,imode)*dtc))
              nd(jl,imode)=nd(jl,imode)-deln4
              IF (nd(jl,imode) < 0.0) nd(jl,imode)=0.0
            ELSE
              deln4=0.0
            END IF

            deln=deln1+deln2+deln3+deln4

            ! Calculate cpt mass tendencies
            dm(:)=0.0
            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                dm(icp)=deln*md(jl,imode,icp)
              END IF
            END DO

            ! .. Store cpt imp scav mass fluxes for budget calculations
            DO icp=1,ncp

              IF (component(imode,icp)) THEN

                IF (icp == cp_su) THEN
                  IF ((imode == mode_nuc_sol) .AND. (nmasimscsunucsol > 0))    &
                    bud_aer_mas(jl,nmasimscsunucsol)=                          &
                      bud_aer_mas(jl,nmasimscsunucsol)+dm(icp)
                  IF ((imode == mode_ait_sol) .AND. (nmasimscsuaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscsuaitsol)=                          &
                      bud_aer_mas(jl,nmasimscsuaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscsuaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscsuaccsol)=                          &
                      bud_aer_mas(jl,nmasimscsuaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscsucorsol > 0))    &
                    bud_aer_mas(jl,nmasimscsucorsol)=                          &
                      bud_aer_mas(jl,nmasimscsucorsol)+dm(icp)
                END IF
                IF (icp == cp_bc) THEN
                  IF ((imode == mode_ait_sol) .AND. (nmasimscbcaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscbcaitsol)=                          &
                      bud_aer_mas(jl,nmasimscbcaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscbcaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscbcaccsol)=                          &
                      bud_aer_mas(jl,nmasimscbcaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscbccorsol > 0))    &
                    bud_aer_mas(jl,nmasimscbccorsol)=                          &
                      bud_aer_mas(jl,nmasimscbccorsol)+dm(icp)
                  IF ((imode == mode_ait_insol) .AND. (nmasimscbcaitins > 0))  &
                    bud_aer_mas(jl,nmasimscbcaitins)=                          &
                      bud_aer_mas(jl,nmasimscbcaitins)+dm(icp)
                END IF
                IF (icp == cp_oc) THEN
                  IF ((imode == mode_nuc_sol) .AND. (nmasimscocnucsol > 0))    &
                    bud_aer_mas(jl,nmasimscocnucsol)=                          &
                      bud_aer_mas(jl,nmasimscocnucsol)+dm(icp)
                  IF ((imode == mode_ait_sol) .AND. (nmasimscocaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscocaitsol)=                          &
                      bud_aer_mas(jl,nmasimscocaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscocaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscocaccsol)=                          &
                      bud_aer_mas(jl,nmasimscocaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscoccorsol > 0))    &
                    bud_aer_mas(jl,nmasimscoccorsol)=                          &
                      bud_aer_mas(jl,nmasimscoccorsol)+dm(icp)
                  IF ((imode == mode_ait_insol) .AND. (nmasimscocaitins > 0))  &
                    bud_aer_mas(jl,nmasimscocaitins)=                          &
                      bud_aer_mas(jl,nmasimscocaitins)+dm(icp)
                END IF
                IF (icp == cp_cl) THEN
                  IF ((imode == mode_acc_sol) .AND. (nmasimscssaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscssaccsol)=                          &
                      bud_aer_mas(jl,nmasimscssaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscsscorsol > 0))    &
                    bud_aer_mas(jl,nmasimscsscorsol)=                          &
                      bud_aer_mas(jl,nmasimscsscorsol)+dm(icp)
                END IF
                IF (icp == cp_so) THEN
                  IF ((imode == mode_nuc_sol) .AND. (nmasimscsonucsol > 0))    &
                    bud_aer_mas(jl,nmasimscsonucsol)=                          &
                      bud_aer_mas(jl,nmasimscsonucsol)+dm(icp)
                  IF ((imode == mode_ait_sol) .AND. (nmasimscsoaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscsoaitsol)=                          &
                      bud_aer_mas(jl,nmasimscsoaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscsoaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscsoaccsol)=                          &
                      bud_aer_mas(jl,nmasimscsoaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscsocorsol > 0))    &
                    bud_aer_mas(jl,nmasimscsocorsol)=                          &
                      bud_aer_mas(jl,nmasimscsocorsol)+dm(icp)
                END IF
                IF (icp == cp_du) THEN
                  IF ((imode == mode_acc_sol) .AND. (nmasimscduaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscduaccsol)=                          &
                      bud_aer_mas(jl,nmasimscduaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscducorsol > 0))    &
                    bud_aer_mas(jl,nmasimscducorsol)=                          &
                      bud_aer_mas(jl,nmasimscducorsol)+dm(icp)
                  IF ((imode == mode_acc_insol) .AND. (nmasimscduaccins > 0))  &
                    bud_aer_mas(jl,nmasimscduaccins)=                          &
                      bud_aer_mas(jl,nmasimscduaccins)+dm(icp)
                  IF ((imode == mode_cor_insol) .AND. (nmasimscducorins > 0))  &
                    bud_aer_mas(jl,nmasimscducorins)=                          &
                      bud_aer_mas(jl,nmasimscducorins)+dm(icp)
                  IF ((imode == mode_sup_insol) .AND. (nmasimscdusupins > 0))  &
                    bud_aer_mas(jl,nmasimscdusupins)=                          &
                     bud_aer_mas(jl,nmasimscdusupins)+dm(icp)
                END IF
                IF (icp == cp_no3) THEN
                  IF ((imode == mode_ait_sol) .AND. (nmasimscntaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscntaitsol)=                          &
                      bud_aer_mas(jl,nmasimscntaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscntaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscntaccsol)=                          &
                      bud_aer_mas(jl,nmasimscntaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscntcorsol > 0))    &
                    bud_aer_mas(jl,nmasimscntcorsol)=                          &
                      bud_aer_mas(jl,nmasimscntcorsol)+dm(icp)
                END IF
                IF (icp == cp_nn) THEN
                  IF ((imode == mode_acc_sol) .AND. (nmasimscnnaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscnnaccsol)=                          &
                      bud_aer_mas(jl,nmasimscnnaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscnncorsol > 0))    &
                    bud_aer_mas(jl,nmasimscnncorsol)=                          &
                      bud_aer_mas(jl,nmasimscnncorsol)+dm(icp)
                END IF
                IF (icp == cp_nh4) THEN
                  IF ((imode == mode_ait_sol) .AND. (nmasimscnhaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscnhaitsol)=                          &
                      bud_aer_mas(jl,nmasimscnhaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscnhaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscnhaccsol)=                          &
                      bud_aer_mas(jl,nmasimscnhaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscnhcorsol > 0))    &
                    bud_aer_mas(jl,nmasimscnhcorsol)=                          &
                      bud_aer_mas(jl,nmasimscnhcorsol)+dm(icp)
                END IF

              END IF ! if component present

            END DO ! loop over components
          END IF ! if ND>NUM_EPS
        END IF ! if distribution defined
      END DO ! loop over distributions
    END IF ! if rain is present in box
  END DO ! loop over boxes



ELSE ! l_fix_ukca_impscav

  !...............................................................................

  ! The original buggy version of the code


  ! .. Here convert rain rate from kg/m2/s to mm/hr
  DO jl=1,nbox
    totrain(jl)=(crain(jl)*3600.0)+(drain(jl)*3600.0)
    totsnow(jl)=(csnow(jl)*3600.0)+(dsnow(jl)*3600.0)

    totppn(jl)=totrain(jl)+totsnow(jl)
  END DO

  ! .. In LUT there are 19 Rows but we only have NRBINS=7 raindrop bins,
  ! .. so use every 2nd row starting from row 4
  ! .. Only consider raindrops radius >4 and <1024 um
  ! .. NB/ HAVE TO CHANGE THIS IF NRBINS IS ALTERED

  jvr=0
  DO irow=4,16,2 ! pick out 7 sizes corresponding to 7 raindrop bins
    jvr=jvr+1
    RNDPDIAM_cm(jvr)=raddrop(irow)*2.0/10000.0
    RNDPDIAM_mm(jvr)=RNDPDIAM_cm(jvr)*10.0
    RNDPDIAM_m(jvr)=RNDPDIAM_mm(jvr)/1000.0
    ! RADDROP contains raindrop bin radii in microns
    ! RNDPDIAM_cm contains raindrop bin diameters in cm
    ! RNDPDIAM_mm & RNDPDIAM_m contain rndrp bin diams in mm & m respect.
  END DO

  ! Initialise Arrays
  scavcoeff_count(:,:,:)=0.0
  scavcoeff(:,:)=0.0
  countcoll(:,:)=0

  DO imode=1,topmode
    IF (mode(imode)) THEN
      DO jl=1,nbox
        IF (totrain(jl) > 0.0) THEN

          ! From aerosol wet radius find which column of LU table is needed
          IF (0.5*wetdp(jl,imode) < 0.001e-6) THEN
            countcoll(jl,imode)=1
          ELSE IF (0.5*wetdp(jl,imode) > 9.0e-6) THEN
            countcoll(jl,imode)=10
          ELSE
            factor=0.457
            r1=0.001
            interb=((0.5*wetdp(jl,imode)*1.0e6)/r1)
            ! .. 1.0E6 in above to convert wet radius to um
            countcoll(jl,imode)=INT((LOG(interb)/factor)+1.0)
          END IF
        END IF
      END DO
    END IF
  END DO

  ! .. Loop to calculate the number of rain drops at each rdrop bin-centre
  ! .. following Marshall Palmer distribution
  DO jvr=1,nrbins
    !
    ! Empirical relationship from Easter and Hales (1984) to calulate
    ! terminal velocity of raindrop in cm/s
    IF (RNDPDIAM_cm(jvr) <= 0.10) THEN
      VELDR_cms(jvr)=4055.0*RNDPDIAM_cm(jvr)
    ELSE IF (RNDPDIAM_cm(jvr) > 0.10) THEN
      VELDR_cms(jvr)=13000.0*(RNDPDIAM_m(jvr)**0.5)
    END IF
    ! .. Convert raindrop velocity from cm/s to m/s
    VELDR_ms(jvr)=VELDR_cms(jvr)/100.0
  END DO

  DO jl=1,nbox
    IF (totrain(jl) > 0.0) THEN
      DO jvr=1,nrbins

        ! Marshall Palmer but with sophistication of NRAINMAX calculated
        ! according to Sekhon & Srivastava (1971) (Seinfeld & Pandis pg 832)

        nrainmax(jl)=7000.0*totrain(jl)**0.37
        ! NRAINMAX is n_0 in m^-3 per mm

        interc(jl,jvr)=3.8*totrain(jl)**(-0.14)*RNDPDIAM_mm(jvr)
        ! INTERC is Psi*D_p in MP distribution


        ndrain(jl,jvr)=nrainmax(jl)*RNDPDIAM_mm(jvr)                           &
                    *0.924*EXP(-interc(jl,jvr))
        ! Includes various conversions as Nd given in terms of mm-1
        ! => *diameter(mm)*width of bin
        ! Geometric scaling factor = 2.52 (change if NRBINS is changed)
        ! (ln(2.52)=0.924)
        !
        ! Leave NDRAIN in m-3 for use in SCAV calculation below
        ! Convert Raindrop size from diameter in mm to diameter in m

        interzz(jl,jvr)=                                                       &
           (pi/4.0)*((RNDPDIAM_m(jvr)*RNDPDIAM_m(jvr))                         &
           *VELDR_ms(jvr)*ndrain(jl,jvr))

      END DO ! over NRBINS
    END IF ! TOTRAIN>0
  END DO ! over NBOX

  ! .. Calculate scavenging coefficients
  DO jl=1,nbox
    IF (totrain(jl) > 0.0) THEN
      DO imode=1,topmode
        IF (mode(imode)) THEN
          IF (nd(jl,imode) > num_eps(imode)) THEN
            ! .. Loop over raindrop bins
            DO jvr=1,nrbins
              irow=(jvr*2)+2 ! hard-wired for 7 raindrop size bins
              !COUNTCOLL contains the index of the relevant column of LU table
              !for aerosol-rndrp coll efficiency for each aerosol bin wet radius
              icoll=countcoll(jl,imode)
              IF (countcoll(jl,imode) == 0) icoll=1
              scav(jl,imode,jvr)=colleff4(icoll,irow)*interzz(jl,jvr)

              !For each aerosol bin, sum scav coeff over all rain bins
              IF (jvr == 1) THEN
                scavcoeff_count(jl,imode,jvr)=scav(jl,imode,jvr)
              END IF
              IF (jvr > 1) THEN
                scavcoeff_count(jl,imode,jvr)=                                 &
                 scavcoeff_count(jl,imode,jvr-1)+scav(jl,imode,jvr)
              END IF
              IF (jvr == nrbins) THEN
                scavcoeff(jl,imode)=scavcoeff_count(jl,imode,jvr)
              END IF

            END DO ! LOOP OVER RAINDROP BINS
          END IF ! IF ND>ND_EPS
        END IF ! IF MODE PRESENT
      END DO ! LOOP OVER MODES
    END IF ! IF TOTRAIN>0

    IF (totsnow(jl) > 0.0) THEN

      ! Tuned values of ASNOW: lower than in used Feng (2009), but
      ! within theoretical and observational uncertainties
      ! (Fig 8 of Feng 2009)

      asnow(:) = [ 0.014, 0.014, 0.014, 0.10, 0.014, 0.014, 0.10, 0.50 ]

      DO imode=1,nmodes
        IF (mode(imode)) THEN
          IF (nd(jl,imode) > num_eps(imode)) THEN
            ! Calculate scavenging coefficients for snow
            ! RSCAV = aR^b  (UNIT: hr-1)

            rscav_snowh(jl,imode) = asnow(imode)*totsnow(jl)**bsnow
            rscav_snow(jl,imode) = rscav_snowh(jl,imode)/3600
            ! Above converts to unit of s-1
          END IF ! IF ND>ND_EPS
        END IF ! IF MODE PRESENT
      END DO ! LOOP OVER MODES

    END IF ! IF totsnow>0
  END DO ! LOOP OVER BOXES

  ! .. Below calculates the rate of removal
  !
  ! .. Calculate removal (from each bin) following 1st order rate loss
  ! .. Apply BCS only to the portion of the gridbox where rain is occuring
  ! .. If dynamic rain then cloud cover = 1.0 => apply BCS to all aerosols
  ! .. If convective rain then cloud cover = FC = 0.3

  DO jl=1,nbox
    IF (totppn(jl) > 0.0) THEN
      DO imode=1,nmodes
        IF (mode(imode)) THEN
          IF (nd(jl,imode) > num_eps(imode)) THEN
            IF (crain(jl) > 0.0) THEN
              deln1=fc*nd(jl,imode)*(1.0-EXP(-scavcoeff(jl,imode)*dtc))
              nd(jl,imode)=nd(jl,imode)-deln1
              IF (nd(jl,imode) < 0.0) nd(jl,imode)=0.0
              ! Convective rain occurs only over fraction FC of box
            ELSE
              deln1=0.0
            END IF
            IF (drain(jl) > 0.0) THEN
              deln2=nd(jl,imode)*(1.0-EXP(-scavcoeff(jl,imode)*dtc))
              nd(jl,imode)=nd(jl,imode)-deln2
              IF (nd(jl,imode) < 0.0) nd(jl,imode)=0.0
              ! Dynamic rain occurs over whole of box
            ELSE
              deln2=0.0
            END IF

            IF (csnow(jl) > 0.0) THEN
              deln3=fc*nd(jl,imode)*(1.0-EXP(-rscav_snow(jl,imode)*dtc))
              nd(jl,imode)=nd(jl,imode)-deln3
              IF (nd(jl,imode) < 0.0) nd(jl,imode)=0.0
              ! Convective SNOW occurs only over fraction FC of box
            ELSE
              deln3=0.0
            END IF
            IF (dsnow(jl) > 0.0) THEN
              deln4=nd(jl,imode)*(1.0-EXP(-rscav_snow(jl,imode)*dtc))
              nd(jl,imode)=nd(jl,imode)-deln4
              IF (nd(jl,imode) < 0.0) nd(jl,imode)=0.0
              ! .. Dynamic rain occurs over whole of box
            ELSE
              deln4=0.0
            END IF

            deln=deln1+deln2+deln3+deln4

            ! Calculate cpt mass tendencies
            dm(:)=0.0
            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                dm(icp)=deln*md(jl,imode,icp)
              END IF
            END DO

            ! Store cpt impaction scavenging mass fluxes for budget calculations
            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                IF (icp == cp_su) THEN
                  IF ((imode == mode_nuc_sol) .AND. (nmasimscsunucsol > 0))    &
                    bud_aer_mas(jl,nmasimscsunucsol)=                          &
                     bud_aer_mas(jl,nmasimscsunucsol)+dm(icp)
                  IF ((imode == mode_ait_sol) .AND. (nmasimscsuaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscsuaitsol)=                          &
                     bud_aer_mas(jl,nmasimscsuaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscsuaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscsuaccsol)=                          &
                     bud_aer_mas(jl,nmasimscsuaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscsucorsol > 0))    &
                    bud_aer_mas(jl,nmasimscsucorsol)=                          &
                     bud_aer_mas(jl,nmasimscsucorsol)+dm(icp)
                END IF
                IF (icp == cp_bc) THEN
                  IF ((imode == mode_ait_sol) .AND. (nmasimscbcaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscbcaitsol)=                          &
                     bud_aer_mas(jl,nmasimscbcaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscbcaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscbcaccsol)=                          &
                     bud_aer_mas(jl,nmasimscbcaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscbccorsol > 0))    &
                    bud_aer_mas(jl,nmasimscbccorsol)=                          &
                     bud_aer_mas(jl,nmasimscbccorsol)+dm(icp)
                  IF ((imode == mode_ait_insol) .AND. (nmasimscbcaitins > 0))  &
                    bud_aer_mas(jl,nmasimscbcaitins)=                          &
                     bud_aer_mas(jl,nmasimscbcaitins)+dm(icp)
                END IF
                IF (icp == cp_oc) THEN
                  IF ((imode == mode_nuc_sol) .AND. (nmasimscocnucsol > 0))    &
                    bud_aer_mas(jl,nmasimscocnucsol)=                          &
                     bud_aer_mas(jl,nmasimscocnucsol)+dm(icp)
                  IF ((imode == mode_ait_sol) .AND. (nmasimscocaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscocaitsol)=                          &
                     bud_aer_mas(jl,nmasimscocaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscocaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscocaccsol)=                          &
                     bud_aer_mas(jl,nmasimscocaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscoccorsol > 0))    &
                    bud_aer_mas(jl,nmasimscoccorsol)=                          &
                     bud_aer_mas(jl,nmasimscoccorsol)+dm(icp)
                  IF ((imode == mode_ait_insol) .AND. (nmasimscocaitins > 0))  &
                    bud_aer_mas(jl,nmasimscocaitins)=                          &
                     bud_aer_mas(jl,nmasimscocaitins)+dm(icp)
                END IF
                IF (icp == cp_cl) THEN
                  IF ((imode == mode_acc_sol) .AND. (nmasimscssaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscssaccsol)=                          &
                     bud_aer_mas(jl,nmasimscssaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscsscorsol > 0))    &
                    bud_aer_mas(jl,nmasimscsscorsol)=                          &
                     bud_aer_mas(jl,nmasimscsscorsol)+dm(icp)
                END IF
                IF (icp == cp_so) THEN
                  IF ((imode == mode_nuc_sol) .AND. (nmasimscsonucsol > 0))    &
                    bud_aer_mas(jl,nmasimscsonucsol)=                          &
                     bud_aer_mas(jl,nmasimscsonucsol)+dm(icp)
                  IF ((imode == mode_ait_sol) .AND. (nmasimscsoaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscsoaitsol)=                          &
                     bud_aer_mas(jl,nmasimscsoaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscsoaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscsoaccsol)=                          &
                     bud_aer_mas(jl,nmasimscsoaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscsocorsol > 0))    &
                    bud_aer_mas(jl,nmasimscsocorsol)=                          &
                     bud_aer_mas(jl,nmasimscsocorsol)+dm(icp)
                END IF
                IF (icp == cp_du) THEN
                  IF ((imode == mode_acc_sol) .AND. (nmasimscduaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscduaccsol)=                          &
                     bud_aer_mas(jl,nmasimscduaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscducorsol > 0))    &
                    bud_aer_mas(jl,nmasimscducorsol)=                          &
                     bud_aer_mas(jl,nmasimscducorsol)+dm(icp)
                  IF ((imode == mode_acc_insol) .AND. (nmasimscduaccins > 0))  &
                    bud_aer_mas(jl,nmasimscduaccins)=                          &
                     bud_aer_mas(jl,nmasimscduaccins)+dm(icp)
                  IF ((imode == mode_cor_insol) .AND. (nmasimscducorins > 0))  &
                    bud_aer_mas(jl,nmasimscducorins)=                          &
                     bud_aer_mas(jl,nmasimscducorins)+dm(icp)
                  IF ((imode == mode_sup_insol) .AND. (nmasimscdusupins > 0))  &
                    bud_aer_mas(jl,nmasimscdusupins)=                          &
                     bud_aer_mas(jl,nmasimscdusupins)+dm(icp)
                END IF
                IF (icp == cp_no3) THEN
                  IF ((imode == mode_ait_sol) .AND. (nmasimscntaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscntaitsol)=                          &
                      bud_aer_mas(jl,nmasimscntaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscntaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscntaccsol)=                          &
                      bud_aer_mas(jl,nmasimscntaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscntcorsol > 0))    &
                    bud_aer_mas(jl,nmasimscntcorsol)=                          &
                      bud_aer_mas(jl,nmasimscntcorsol)+dm(icp)
                END IF
                IF (icp == cp_nn) THEN
                  IF ((imode == mode_acc_sol) .AND. (nmasimscnnaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscnnaccsol)=                          &
                      bud_aer_mas(jl,nmasimscnnaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscnncorsol > 0))    &
                    bud_aer_mas(jl,nmasimscnncorsol)=                          &
                      bud_aer_mas(jl,nmasimscnncorsol)+dm(icp)
                END IF
                IF (icp == cp_nh4) THEN
                  IF ((imode == mode_ait_sol) .AND. (nmasimscnhaitsol > 0))    &
                    bud_aer_mas(jl,nmasimscnhaitsol)=                          &
                      bud_aer_mas(jl,nmasimscnhaitsol)+dm(icp)
                  IF ((imode == mode_acc_sol) .AND. (nmasimscnhaccsol > 0))    &
                    bud_aer_mas(jl,nmasimscnhaccsol)=                          &
                      bud_aer_mas(jl,nmasimscnhaccsol)+dm(icp)
                  IF ((imode == mode_cor_sol) .AND. (nmasimscnhcorsol > 0))    &
                    bud_aer_mas(jl,nmasimscnhcorsol)=                          &
                      bud_aer_mas(jl,nmasimscnhcorsol)+dm(icp)
                END IF

              END IF ! if component present
            END DO ! loop over components
          END IF ! if ND>NUM_EPS
        END IF ! if distribution defined
      END DO ! loop over distributions
    END IF ! if rain is present in box
  END DO ! loop over boxes

  !...............................................................................

END IF ! l_fix_ukca_impscav

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_impc_scav
END MODULE ukca_impc_scav_mod
