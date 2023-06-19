! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Subroutine to calculate impaction scavenging of aerosols
!    by falling raindrops. This is a two-moment scavenging scheme
!    for liquid raindroplets only, and currently only used
!    for insoluble dust
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
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ######################################################################
!
! .. Subroutine Interface:
MODULE ukca_impc_scav_dust_mod

USE yomhook,            ONLY: lhook, dr_hook
USE parkind1,           ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

! .. Variables for impaction scavenging
INTEGER, PARAMETER  :: ncol = 22 ! # of aero diameter integration pts
INTEGER, PARAMETER  :: nrow = 22 ! # of rainfall rate integration pts
INTEGER, PARAMETER  :: nlev = 5  ! # of aero stdev integration pts

! .. Parameters for st. dev. integration points
! .. Formulae: (i-1)*dt_stdev + stdev_low for i=1,..,nlev
REAL, PARAMETER :: stdev_low = 1.2
REAL, PARAMETER :: dt_stdev = 0.2

! .. Parameters for median diameter integration points
! .. Formulae: 2 * 10**((i-1)*dt_ldp + ldplow) for i=1,..,ncol
REAL, PARAMETER :: ldplow = -9.0
REAL, PARAMETER :: dt_ldp = 0.2

! .. Parameters to define rainrate integration points
! .. Formulate: 10**((i-1)*dt_lrf + lrflow) for i=1,..,nrow
REAL, PARAMETER :: lrflow = -1.0
REAL, PARAMETER :: dt_lrf = 1.0/7.0

! .. Scavenging  coefficients (s-1) for aerosol number and mass
! .. as a function of median diameter, rainrate and modal stdev
! .. Populated in ukca_impc_scav_dust_init() below
REAL, ALLOCATABLE, SAVE :: scav_coeff_num(:,:,:)
REAL, ALLOCATABLE, SAVE :: scav_coeff_mass(:,:,:)

! Routines available from this module
PUBLIC :: ukca_impc_scav_dust, ukca_impc_scav_dust_init,                       &
          ukca_impc_scav_dust_dealloc

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_IMPC_SCAV_DUST_MOD'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE ukca_impc_scav_dust(nbox,nbudaer,nd,md,mdt,crain,drain,             &
                               wetdp,dtc,bud_aer_mas)
! ----------------------------------------------------------------------
!
! Purpose:
! -------
!     Subroutine to calculate impaction scavenging of aerosols
!     by falling raindrops.
!
!     The method differs to ukca_impc_scav and depends on the values
!     in the scav_coeff_num and scav_coeff_mass arrays. Scavenging
!     coefficients are either derived from Slinn (1983) formulae
!     for collision efficiency, integrated over the rain droplet
!     size distribution (Abel and Boutle, 2012) and using droplet
!     fall speeds from Beard (1976) or provided directly from a
!     parameterisation in Laakso et al (2003)
!
!     The rainfall rate and the geometric median diameters integration
!     points are dispersed logarithmically while the standard deviation
!     integration points are (1.2, 1.4, ..., 2). Interpolation of the
!     scavenging coefficient arrays is log10-log10 with diameter (i.e.
!     that both diameter and scav coefficient are interpolated in log10
!     space) and linear with rain rate, and using a nearest-neighbour
!     approach for standard deviation
!
!     Parameters
!     ----------
!
!     Inputs
!     ------
!     NBOX        : Number of grid boxes
!     NBUDAER     : Number of aerosol budget fields
!     ND          : Aerosol ptcl number density (cm^-3)
!     MD          : Avg cpt mass of aerosol ptcl (particle^-1)
!     MDT         : Total median aerosol mass (molecules per ptcl)
!     CRAIN       : Convective rain rate array (kgm^-2s^-1)
!     DRAIN       : Dynamic rain rate array (kgm^-2s^-1)
!     WETDP       : Wet diameter corresponding to DRYDP (m)
!     DTC         : Time step of process (s)
!
!     Outputs
!     -------
!     ND          : new aerosol number conc (cm^-3)
!     MD          : new cpt mass of aerosol ptcl (particle^-1)
!     MDT         : Total median aerosol mass (molecules per ptcl)
!     BUD_AER_MAS : Updated aerosol budgets
!
!----------------------------------------------------------------------

USE ukca_config_specification_mod, ONLY: glomap_variables

USE ukca_mode_setup,    ONLY: cp_du, nmodes, mode_acc_insol, mode_cor_insol,   &
                              mode_sup_insol

USE ukca_setup_indices, ONLY: nmasimscduaccins, nmasimscducorins,              &
                              nmasimscdusupins

IMPLICIT NONE

! .. Subroutine interface
INTEGER, INTENT(IN)  :: nbox
INTEGER, INTENT(IN)  :: nbudaer
REAL, INTENT(IN)     :: wetdp(nbox,nmodes)
REAL, INTENT(IN)     :: dtc
REAL, INTENT(IN)     :: crain(nbox)
REAL, INTENT(IN)     :: drain(nbox)
REAL, INTENT(IN OUT) :: nd(nbox,nmodes)
REAL, INTENT(IN OUT) :: md(nbox,nmodes,glomap_variables%ncp)
REAL, INTENT(IN OUT) :: mdt(nbox,nmodes)
REAL, INTENT(IN OUT) :: bud_aer_mas(nbox,0:nbudaer)

! .. Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: mfrac_0(:,:)
REAL,    POINTER :: mmid(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: num_eps(:)
REAL,    POINTER :: sigmag(:)

REAL, PARAMETER :: fc=0.3
REAL, PARAMETER :: fd=1.0
REAL, PARAMETER :: secs_per_hr=3600.0
REAL    :: allfrac(2)
REAL    :: allrain(2,nbox)
REAL    :: totrain(nbox)
REAL    :: scavn(2,nbox,nmodes)
REAL    :: scavm(2,nbox,nmodes)
REAL    :: deln, deln1, deln2
REAL    :: dm1(glomap_variables%ncp)
REAL    :: dm2(glomap_variables%ncp)
REAL    :: dm(glomap_variables%ncp)
REAL    :: ndnew
REAL    :: ii, dplow, dpupp, jj, rflow, rfupp
REAL    :: logdp, logdplow, logdpupp, fac1, fac2
REAL    :: sc_num_dp_rf1, sc_num_dp_rf2
REAL    :: sc_mass_dp_rf1, sc_mass_dp_rf2
LOGICAL :: l_interp_dp, l_interp_RF
INTEGER :: ilow,iupp,jlow,jupp,k
INTEGER :: imode, icp, jl, iprecip

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_IMPC_SCAV_DUST'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mfrac_0     => glomap_variables%mfrac_0
mmid        => glomap_variables%mmid
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
num_eps     => glomap_variables%num_eps
sigmag      => glomap_variables%sigmag

! .. Combine the convective and dynamic rain in an array
! .. to loop over and convert rain from kgm-2s-1 to mm/hr
allfrac(1)=fc
allfrac(2)=fd
allrain(:,:)=0.0
totrain(:)=0.0
DO jl=1,nbox
  allrain(1,jl)=crain(jl)*secs_per_hr/allfrac(1)
  allrain(2,jl)=drain(jl)*secs_per_hr/allfrac(2)
  totrain(jl)=(crain(jl)*secs_per_hr)+(drain(jl)*secs_per_hr)
END DO

! .. Initialise 3D scavenging arrays
scavn(:,:,:) = 0.0
scavm(:,:,:) = 0.0

! .. Loop over mode, rain type, and gridcell to populate the
! .. 3D real-time scavenging arrays
DO imode = mode_acc_insol, mode_sup_insol
  IF (mode(imode)) THEN

    ! .. Nearest neighbour interpolation is used for standard
    ! .. deviation in scavenging arrays so only need to do
    ! .. this once for each mode
    k = MAX(NINT((sigmag(imode)-stdev_low+dt_stdev)/dt_stdev),1)
    k = MIN(k,nlev)

    DO jl=1,nbox
      DO iprecip=1,2
        IF (allrain(iprecip,jl) > 0.0) THEN

          ! .. Find the i index (diameter) for interpolation of
          ! .. scavenging coefficient arrays
          l_interp_dp = .FALSE.
          ii = 1 + (LOG10(wetdp(jl,imode)/2.0)-ldplow) / dt_ldp
          IF (ii <= 1) THEN
            ilow = 1
            iupp = 1
            dplow = 2.0*(10.0**ldplow)
            dpupp = dplow
          ELSE IF (ii >= ncol) THEN
            ilow = ncol
            iupp = ncol
            dplow = 2.0*(10.0**(ldplow+(ilow-1)*dt_ldp))
            dpupp = dplow
          ELSE
            ilow = INT(ii)
            iupp = ilow+1
            dplow = 2.0*(10.0**(ldplow+(ilow-1)*dt_ldp))
            dpupp = 2.0*(10.0**(ldplow+(iupp-1)*dt_ldp))
            l_interp_dp = .TRUE.
          END IF

          ! .. Find the j index (rainfall) for interpolation of
          ! .. scavenging coefficient arrays
          l_interp_RF = .FALSE.
          jj = 1 + (LOG10(allrain(iprecip,jl))-lrflow) / dt_lrf
          IF (jj <= 1) THEN
            jlow = 1
            jupp = 1
            rflow = 10.0**lrflow
            rfupp = rflow
          ELSE IF (jj >= nrow) THEN
            jlow = nrow
            jupp = nrow
            rflow = 10.0**(lrflow+(jlow-1)*dt_lrf)
            rfupp = rflow
          ELSE
            jlow = INT(jj)
            jupp = jlow+1
            rflow = 10.0**(lrflow+(jlow-1)*dt_lrf)
            rfupp = 10.0**(lrflow+(jupp-1)*dt_lrf)
            l_interp_RF = .TRUE.
          END IF

          ! .. Now interpolate in i index (LOG10-LOG10)
          IF (l_interp_dp) THEN
            logdp = LOG10(wetdp(jl,imode))
            logdplow = LOG10(dplow)
            logdpupp = LOG10(dpupp)
            fac1 = (logdpupp-logdp)/(logdpupp-logdplow)
            fac2 = (logdp-logdplow)/(logdpupp-logdplow)
            sc_num_dp_rf1 =                                                    &
                     10.0**(fac1*LOG10(scav_coeff_num(ilow,jlow,k)) +          &
                            fac2*LOG10(scav_coeff_num(iupp,jlow,k)))
            sc_num_dp_rf2 =                                                    &
                     10.0**(fac1*LOG10(scav_coeff_num(ilow,jupp,k)) +          &
                            fac2*LOG10(scav_coeff_num(iupp,jupp,k)))
            sc_mass_dp_rf1 =                                                   &
                     10.0**(fac1*LOG10(scav_coeff_mass(ilow,jlow,k)) +         &
                            fac2*LOG10(scav_coeff_mass(iupp,jlow,k)))
            sc_mass_dp_rf2 =                                                   &
                     10.0**(fac1*LOG10(scav_coeff_mass(ilow,jupp,k)) +         &
                            fac2*LOG10(scav_coeff_mass(iupp,jupp,k)))
          ELSE
            sc_num_dp_rf1 = scav_coeff_num(ilow,jlow,k)
            sc_num_dp_rf2 = scav_coeff_num(ilow,jupp,k)
            sc_mass_dp_rf1 = scav_coeff_mass(ilow,jlow,k)
            sc_mass_dp_rf2 = scav_coeff_mass(ilow,jupp,k)
          END IF

          ! .. Lastly interpolate in the j index (LIN-LIN)
          IF (l_interp_RF) THEN
            fac1 = (rfupp-allrain(iprecip,jl))/(rfupp-rflow)
            fac2 = (allrain(iprecip,jl)-rflow)/(rfupp-rflow)
            scavn(iprecip,jl,imode) = fac1 * sc_num_dp_rf1 +                   &
                                      fac2 * sc_num_dp_rf2
            scavm(iprecip,jl,imode) = fac1 * sc_mass_dp_rf1 +                  &
                                      fac2 * sc_mass_dp_rf2
          ELSE
            scavn(iprecip,jl,imode) = sc_num_dp_rf1
            scavm(iprecip,jl,imode) = sc_mass_dp_rf1
          END IF

        END IF ! .. IF (allrain(iprecip,jl) > 0.0)
      END DO ! .. iprecip
    END DO ! .. nbox
  END IF ! .. IF (mode(imode))
END DO ! .. imode

! .. Apply derived scavenging coefficients to input mass/number
DO jl=1,nbox
  IF (totrain(jl) > 0.0) THEN
    DO imode=mode_acc_insol, mode_sup_insol
      IF (mode(imode)) THEN

        ! .. Only do anything if the initial number is greater than
        ! .. a threshold value (num_eps)
        IF (nd(jl,imode) > num_eps(imode)) THEN

          ! .. Nullify mass and number changes
          DO icp=1,ncp
            dm1(icp)=0.0
            dm2(icp)=0.0
          END DO
          deln1=0.0
          deln2=0.0

          ! .. Convective rain
          IF (crain(jl) > 0.0) THEN
            deln1 = allfrac(1)*nd(jl,imode)*(1.0-EXP(-scavn(1,jl,imode)*dtc))
            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                dm1(icp) = allfrac(1)*nd(jl,imode)*md(jl,imode,icp)*           &
                           (1.0-EXP(-scavm(1,jl,imode)*dtc))
              END IF
            END DO
          END IF

          ! .. Dynamical rain
          IF (drain(jl) > 0.0) THEN
            deln2 = allfrac(2)*nd(jl,imode)*(1.0-EXP(-scavn(2,jl,imode)*dtc))
            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                dm2(icp) = allfrac(2)*nd(jl,imode)*md(jl,imode,icp)*           &
                           (1.0-EXP(-scavm(2,jl,imode)*dtc))
              END IF
            END DO
          END IF

          ! .. Sum mass and number changes from dynamic and convective rain
          deln=MIN(nd(jl,imode),deln1+deln2)
          DO icp=1,ncp
            IF (component(imode,icp)) THEN
              dm(icp) = MIN(nd(jl,imode)*md(jl,imode,icp),dm1(icp)+dm2(icp))
            END IF
          END DO

          ! .. Update number and mass concentrations
          mdt(jl,imode)=0.0
          ndnew = nd(jl,imode)-deln
          IF (ndnew > num_eps(imode)) THEN
            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                md(jl,imode,icp)=(nd(jl,imode)*md(jl,imode,icp)-dm(icp))/ndnew
                mdt(jl,imode)=mdt(jl,imode)+md(jl,imode,icp)
              END IF
            END DO
            nd(jl,imode) = ndnew
          ELSE
            deln = nd(jl,imode)
            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                dm(icp)=nd(jl,imode)*md(jl,imode,icp)
                md(jl,imode,icp) = mmid(imode)*mfrac_0(imode,icp)
              END IF
            END DO
            nd(jl,imode) = 0.0
            mdt(jl,imode) = mmid(imode)
          END IF

          ! .. Store cpt imp scav mass fluxes for budget calculations
          DO icp=1,ncp
            IF (component(imode,icp)) THEN
              IF (icp == cp_du) THEN
                IF ((imode == mode_acc_insol) .AND. (nmasimscduaccins > 0))    &
                  bud_aer_mas(jl,nmasimscduaccins)=                            &
                    bud_aer_mas(jl,nmasimscduaccins)+dm(icp)
                IF ((imode == mode_cor_insol) .AND. (nmasimscducorins > 0))    &
                  bud_aer_mas(jl,nmasimscducorins)=                            &
                    bud_aer_mas(jl,nmasimscducorins)+dm(icp)
                IF ((imode == mode_sup_insol) .AND. (nmasimscdusupins > 0))    &
                  bud_aer_mas(jl,nmasimscdusupins)=                            &
                    bud_aer_mas(jl,nmasimscdusupins)+dm(icp)
              END IF
            END IF ! .. if component present
          END DO ! .. icp

        END IF ! .. IF (nd(jl,imode) > num_eps(imode))
      END IF ! .. IF (mode(imode))
    END DO ! .. imode
  END IF !.. IF (totrain(jl) > 0.0)
END DO ! .. jl

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_impc_scav_dust

! ----------------------------------------------------------------------
SUBROUTINE ukca_impc_scav_dust_init(verbose)
! ----------------------------------------------------------------------
! Purpose:
! -------
! Initialise scavenging coefficient lookup tables
! ----------------------------------------------------------------------

USE umPrintMgr,                    ONLY: umPrint, umMessage

IMPLICIT NONE

INTEGER, INTENT(IN) :: verbose ! flag to indicate level of verbosity

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_IMPC_SCAV_DUST_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! New dust impaction scavenging scheme - this module
IF (.NOT. ALLOCATED(scav_coeff_num))                                           &
  ALLOCATE (scav_coeff_num(1:ncol,1:nrow,1:nlev))
IF (.NOT. ALLOCATED(scav_coeff_mass))                                          &
  ALLOCATE (scav_coeff_mass(1:ncol,1:nrow,1:nlev))

IF (verbose >= 2) THEN
  WRITE(umMessage,'(A50,2I5)') 'New dust impaction scavenging scheme on '//    &
                               'in this setup'
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Populate the new scavenging arrays
scav_coeff_num(1:ncol, 1, 1) = [                                               &
      6.1037148e-05,3.4363694e-05,1.9950263e-05,1.1905033e-05,7.2960632e-06,   &
      4.6040514e-06,3.0113628e-06,2.0628355e-06,1.4976336e-06,1.1632750e-06,   &
      9.7163390e-07,8.7728381e-07,8.6524578e-07,9.4244138e-07,1.1398267e-06,   &
      1.6519928e-06,1.2224743e-05,4.3962760e-05,5.5330641e-05,5.4394218e-05,   &
      6.0671539e-05,9.0363565e-05 ]
scav_coeff_num(1:ncol, 1, 2) = [                                               &
      6.5706309e-05,3.6689561e-05,2.1157495e-05,1.2556026e-05,7.6590291e-06,   &
      4.8121111e-06,3.1332367e-06,2.1353433e-06,1.5412965e-06,1.1902741e-06,   &
      9.9004062e-07,8.9285139e-07,8.8212565e-07,9.6559994e-07,1.2355951e-06,   &
      3.4073185e-06,1.6383684e-05,4.0091787e-05,5.3066348e-05,5.5751309e-05,   &
      6.5116983e-05,1.0058589e-04 ]
scav_coeff_num(1:ncol, 1, 3) = [                                               &
      7.2835969e-05,4.0195060e-05,2.2954348e-05,1.3514131e-05,8.1881546e-06,   &
      5.1130763e-06,3.3084738e-06,2.2391734e-06,1.6037626e-06,1.2290530e-06,   &
      1.0165489e-06,9.1507999e-07,9.0711056e-07,1.0318727e-06,1.7965781e-06,   &
      6.0411697e-06,1.8899684e-05,3.7585026e-05,5.0641188e-05,5.7479621e-05,   &
      7.1931053e-05,1.1693697e-04 ]
scav_coeff_num(1:ncol, 1, 4) = [                                               &
      8.2525145e-05,4.4882987e-05,2.5319472e-05,1.4756997e-05,8.8659539e-06,   &
      5.4946394e-06,3.5288539e-06,2.3690317e-06,1.6817423e-06,1.2775862e-06,   &
      1.0498017e-06,9.4370470e-07,9.5525410e-07,1.2631306e-06,2.8580057e-06,   &
      8.4259391e-06,2.0457544e-05,3.6127692e-05,4.9081412e-05,5.9891848e-05,   &
      8.1328636e-05,1.4037703e-04 ]
scav_coeff_num(1:ncol, 1, 5) = [                                               &
      9.5177508e-05,5.0897106e-05,2.8299513e-05,1.6296551e-05,9.6929965e-06,   &
      5.9543943e-06,3.7917755e-06,2.5228860e-06,1.7738382e-06,1.3349742e-06,   &
      1.0898021e-06,9.8577057e-07,1.0692644e-06,1.7122465e-06,4.1157658e-06,   &
      1.0357743e-05,2.1534589e-05,3.5427627e-05,4.8691800e-05,6.3622664e-05,   &
      9.4055582e-05,1.7273542e-04 ]
scav_coeff_num(1:ncol, 2, 1) = [                                               &
      6.3630249e-05,3.5905266e-05,2.0864780e-05,1.2451834e-05,7.6307078e-06,   &
      4.8191336e-06,3.1615601e-06,2.1794014e-06,1.5964712e-06,1.2498167e-06,   &
      1.0454355e-06,9.3824472e-07,9.1931965e-07,1.0024244e-06,1.2237700e-06,   &
      1.8837638e-06,1.7546948e-05,5.7798820e-05,6.9506943e-05,6.7048700e-05,   &
      7.3235012e-05,1.0653566e-04 ]
scav_coeff_num(1:ncol, 2, 2) = [                                               &
      6.8425657e-05,3.8311325e-05,2.2121036e-05,1.3132244e-05,8.0111851e-06,   &
      5.0375335e-06,3.2893441e-06,2.2549185e-06,1.6411685e-06,1.2767445e-06,   &
      1.0638332e-06,9.5493724e-07,9.3861199e-07,1.0292958e-06,1.3619858e-06,   &
      4.4893250e-06,2.2390521e-05,5.2423241e-05,6.6889678e-05,6.8672629e-05,   &
      7.8383078e-05,1.1833013e-04 ]
scav_coeff_num(1:ncol, 2, 3) = [                                               &
      7.5728623e-05,4.1929178e-05,2.3987342e-05,1.4132263e-05,8.5653106e-06,   &
      5.3532468e-06,3.4730013e-06,2.3630820e-06,1.7052896e-06,1.3157713e-06,   &
      1.0906415e-06,9.7873724e-07,9.6731271e-07,1.1188625e-06,2.1728461e-06,   &
      8.0682766e-06,2.5227189e-05,4.8777900e-05,6.3873277e-05,7.0709515e-05,   &
      8.6274862e-05,1.3719369e-04 ]
scav_coeff_num(1:ncol, 2, 4) = [                                               &
      8.5620121e-05,4.6752579e-05,2.6437645e-05,1.5426999e-05,9.2741656e-06,   &
      5.7531346e-06,3.7038482e-06,2.4984124e-06,1.7855925e-06,1.3650462e-06,   &
      1.1246306e-06,1.0098031e-06,1.0283858e-06,1.4431214e-06,3.6365487e-06,   &
      1.1170250e-05,2.6928043e-05,4.6565478e-05,6.1758467e-05,7.3460738e-05,   &
      9.7141839e-05,1.6423156e-04 ]
scav_coeff_num(1:ncol, 2, 5) = [                                               &
      9.8488088e-05,5.2918404e-05,2.9515462e-05,1.7026882e-05,1.0137556e-05,   &
      6.2343828e-06,3.9790803e-06,2.6588234e-06,1.8807314e-06,1.4237662e-06,   &
      1.1661790e-06,1.0587062e-06,1.1815658e-06,2.0598127e-06,5.3174238e-06,   &
      1.3617522e-05,2.8071647e-05,4.5370769e-05,6.0995795e-05,7.7665489e-05,   &
      1.1182145e-04,2.0154434e-04 ]
scav_coeff_num(1:ncol, 3, 1) = [                                               &
      6.6764925e-05,3.7768020e-05,2.1972740e-05,1.3118339e-05,8.0435134e-06,   &
      5.0905131e-06,3.3582677e-06,2.3395106e-06,1.7383607e-06,1.3775219e-06,   &
      1.1545280e-06,1.0251831e-06,9.9051368e-07,1.0748815e-06,1.3205944e-06,   &
      2.1970309e-06,2.4881232e-05,7.5190835e-05,8.7106001e-05,8.2696223e-05,   &
      8.8475609e-05,1.2553365e-04 ]
scav_coeff_num(1:ncol, 3, 2) = [                                               &
      7.1717208e-05,4.0271805e-05,2.3288238e-05,1.3834254e-05,8.4452229e-06,   &
      5.3215353e-06,3.4932499e-06,2.4185050e-06,1.7838321e-06,1.4036331e-06,   &
      1.1722734e-06,1.0430865e-06,1.0130085e-06,1.1066022e-06,1.5214396e-06,   &
      5.9807808e-06,3.0316801e-05,6.8018367e-05,8.4075898e-05,8.4608445e-05,   &
      9.4414322e-05,1.3908811e-04 ]
scav_coeff_num(1:ncol, 3, 3) = [                                               &
      7.9237922e-05,4.4027369e-05,2.5238717e-05,1.4884919e-05,9.0296396e-06,   &
      5.6552044e-06,3.6871153e-06,2.5316617e-06,1.8493319e-06,1.4420764e-06,   &
      1.1986962e-06,1.0685709e-06,1.0466121e-06,1.2291029e-06,2.6768051e-06,   &
      1.0771180e-05,3.3433724e-05,6.2938711e-05,8.0352298e-05,8.6977394e-05,   &
      1.0351285e-04,1.6076000e-04 ]
scav_coeff_num(1:ncol, 3, 4) = [                                               &
      8.9388215e-05,4.9018315e-05,2.7792733e-05,1.6242450e-05,9.7761106e-06,   &
      6.0773440e-06,3.9305950e-06,2.6732936e-06,1.9317611e-06,1.4913466e-06,   &
      1.2328316e-06,1.1024455e-06,1.1254178e-06,1.6804280e-06,4.6703803e-06,   &
      1.4765215e-05,3.5245300e-05,5.9751670e-05,7.7543680e-05,9.0085815e-05,   &
      1.1601922e-04,1.9181368e-04 ]
scav_coeff_num(1:ncol, 3, 5) = [                                               &
      1.0253986e-04,5.5374294e-05,3.0990421e-05,1.7915601e-05,1.0683545e-05,   &
      6.5846386e-06,4.2206325e-06,2.8412692e-06,2.0298892e-06,1.5508223e-06,   &
      1.2756222e-06,1.1601750e-06,1.3315458e-06,2.5194053e-06,6.8952972e-06,   &
      1.7843860e-05,3.6426548e-05,5.7903430e-05,7.6285200e-05,9.4786815e-05,   &
      1.3287082e-04,2.3465061e-04 ]
scav_coeff_num(1:ncol, 4, 1) = [                                               &
      7.0469649e-05,3.9968562e-05,2.3284417e-05,1.3911934e-05,8.5413885e-06,   &
      5.4264945e-06,3.6126266e-06,2.5579084e-06,1.9412913e-06,1.5657118e-06,   &
      1.3163211e-06,1.1505808e-06,1.0857031e-06,1.1628835e-06,1.4318934e-06,   &
      2.6245843e-06,3.4781982e-05,9.6934846e-05,1.0895754e-04,1.0205357e-04,   &
      1.0698456e-04,1.4786945e-04 ]
scav_coeff_num(1:ncol, 4, 2) = [                                               &
      7.5611234e-05,4.2588656e-05,2.4670026e-05,1.4669901e-05,8.9684048e-06,   &
      5.6726922e-06,3.7562388e-06,2.6407795e-06,1.9869166e-06,1.5896531e-06,   &
      1.3322056e-06,1.1696531e-06,1.1125274e-06,1.2011142e-06,1.7234466e-06,   &
      8.0071948e-06,4.0650586e-05,8.7636738e-05,1.0542731e-04,1.0427307e-04,   &
      1.1381060e-04,1.6338662e-04 ]
scav_coeff_num(1:ncol, 4, 3) = [                                               &
      8.3396747e-05,4.6508803e-05,2.6720336e-05,1.5780583e-05,9.5888859e-06,   &
      6.0278768e-06,3.9622694e-06,2.7594745e-06,2.0530482e-06,1.6259471e-06,   &
      1.3569096e-06,1.1967726e-06,1.1526806e-06,1.3696665e-06,3.3459956e-06,   &
      1.4333410e-05,4.3985798e-05,8.0771823e-05,1.0084737e-04,1.0699049e-04,   &
      1.2425512e-04,1.8818449e-04 ]
scav_coeff_num(1:ncol, 4, 4) = [                                               &
      9.3865615e-05,5.1701287e-05,2.9397774e-05,1.7212622e-05,1.0380101e-05,   &
      6.4765864e-06,4.2207127e-06,2.9080839e-06,2.1368933e-06,1.6737241e-06,   &
      1.3899735e-06,1.2337282e-06,1.2556044e-06,1.9922748e-06,6.0297663e-06,   &
      1.9430968e-05,4.5863510e-05,7.6342969e-05,9.7169834e-05,1.1046360e-04,   &
      1.3858109e-04,2.2370008e-04 ]
scav_coeff_num(1:ncol, 4, 5) = [                                               &
      1.0737297e-04,5.8288164e-05,3.2738823e-05,1.8972887e-05,1.1339904e-05,   &
      7.0148886e-06,4.5282008e-06,3.0844549e-06,2.2374451e-06,1.7326860e-06,   &
      1.4331763e-06,1.3026657e-06,1.5324700e-06,3.1229744e-06,8.9489268e-06,   &
      2.3281278e-05,4.7043931e-05,7.3644645e-05,9.5258454e-05,1.1567048e-04,   &
      1.5783441e-04,2.7266614e-04 ]
scav_coeff_num(1:ncol, 5, 1) = [                                               &
      7.4777212e-05,4.2526200e-05,2.4811988e-05,1.4841745e-05,9.1334172e-06,   &
      5.8386359e-06,3.9407403e-06,2.8564586e-06,2.2324790e-06,1.8441613e-06,   &
      1.5581269e-06,1.3343411e-06,1.2159424e-06,1.2712274e-06,1.5598280e-06,   &
      3.2084270e-06,4.7897526e-05,1.2400555e-04,1.3608372e-04,1.2600302e-04,   &
      1.2948627e-04,1.7415961e-04 ]
scav_coeff_num(1:ncol, 5, 2) = [                                               &
      8.0142562e-05,4.5282389e-05,2.6279349e-05,1.5648879e-05,9.5902893e-06,   &
      6.1029521e-06,4.0946420e-06,2.9435275e-06,2.2771123e-06,1.8636369e-06,   &
      1.5700150e-06,1.3542698e-06,1.2487074e-06,1.3184204e-06,1.9800422e-06,   &
      1.0719179e-05,5.3973254e-05,1.1220209e-04,1.3193118e-04,1.2854522e-04,   &
      1.3730564e-04,1.9185783e-04 ]
scav_coeff_num(1:ncol, 5, 3) = [                                               &
      8.8242873e-05,4.9395731e-05,2.8446256e-05,1.6829749e-05,1.0253250e-05,   &
      6.4837137e-06,4.3150590e-06,3.0681599e-06,2.3424317e-06,1.8950741e-06,   &
      1.5906028e-06,1.3826408e-06,1.2976860e-06,1.5505448e-06,4.2260924e-06,   &
      1.8974534e-05,5.7443782e-05,1.0313369e-04,1.2630235e-04,1.3161768e-04,   &
      1.4924604e-04,2.2012185e-04 ]
scav_coeff_num(1:ncol, 5, 4) = [                                               &
      9.9093835e-05,5.4825946e-05,3.1268211e-05,1.8348980e-05,1.1097084e-05,   &
      6.9638621e-06,4.5910498e-06,3.2242133e-06,2.4262125e-06,1.9387085e-06,   &
      1.6203443e-06,1.4226955e-06,1.4330785e-06,2.4009823e-06,7.7996029e-06,   &
      2.5431806e-05,5.9330725e-05,9.7138366e-05,1.2153435e-04,1.3544791e-04,   &
      1.6558163e-04,2.6057389e-04 ]
scav_coeff_num(1:ncol, 5, 5) = [                                               &
      1.1303342e-04,6.1687094e-05,3.4777747e-05,2.0211323e-05,1.2118403e-05,   &
      7.5386978e-06,4.9188563e-06,3.4095557e-06,2.5278463e-06,1.9947970e-06,   &
      1.6622435e-06,1.5055617e-06,1.8033157e-06,3.9106121e-06,1.1598594e-05,   &
      3.0224368e-05,6.0461811e-05,9.3347932e-05,1.1876960e-04,1.4115248e-04,   &
      1.8747535e-04,3.1630875e-04 ]
scav_coeff_num(1:ncol, 6, 1) = [                                               &
      7.9724125e-05,4.5462621e-05,2.6569414e-05,1.5918739e-05,9.8312923e-06,   &
      6.3427052e-06,4.3654193e-06,3.2669310e-06,2.6522982e-06,2.2579003e-06,   &
      1.9220492e-06,1.6077184e-06,1.3988216e-06,1.4074458e-06,1.7074232e-06,   &
      3.9986504e-06,6.4965897e-05,1.5758097e-04,1.6972980e-04,1.5562152e-04,   &
      1.5686130e-04,2.0514049e-04 ]
scav_coeff_num(1:ncol, 6, 2) = [                                               &
      8.5349942e-05,4.8376032e-05,2.8131064e-05,1.6782857e-05,1.0323201e-05,   &
      6.6286409e-06,4.5316290e-06,3.3584494e-06,2.6940494e-06,2.2691918e-06,   &
      1.9263204e-06,1.6276620e-06,1.4398123e-06,1.4672280e-06,2.3059610e-06,   &
      1.4290427e-05,7.0962396e-05,1.4282305e-04,1.6478891e-04,1.5849632e-04,   &
      1.6578976e-04,2.2525382e-04 ]
scav_coeff_num(1:ncol, 6, 3) = [                                               &
      9.3818301e-05,5.2713127e-05,3.0432616e-05,1.8045069e-05,1.1035906e-05,   &
      7.0397556e-06,4.7690696e-06,3.4892541e-06,2.7561133e-06,2.2913297e-06,   &
      1.9386878e-06,1.6562619e-06,1.5007213e-06,1.7858313e-06,5.3708708e-06,   &
      2.4949868e-05,7.4467829e-05,1.3105141e-04,1.5786268e-04,1.6191659e-04,   &
      1.7938623e-04,2.5734624e-04 ]
scav_coeff_num(1:ncol, 6, 4) = [                                               &
      1.0511899e-04,5.8419720e-05,3.3421782e-05,1.9665366e-05,1.1941218e-05,   &
      7.5569871e-06,5.0655865e-06,3.6529533e-06,2.8372300e-06,2.3264001e-06,   &
      1.9611958e-06,1.6988872e-06,1.6792738e-06,2.9351767e-06,1.0079747e-05,   &
      3.3077694e-05,7.6296483e-05,1.2309654e-04,1.5172389e-04,1.6607541e-04,   &
      1.9792970e-04,3.0323612e-04 ]
scav_coeff_num(1:ncol, 6, 5) = [                                               &
      1.1957226e-04,6.5601444e-05,3.7126813e-05,2.1646109e-05,1.3034262e-05,   &
      8.1746664e-06,5.4169187e-06,3.8474987e-06,2.9374625e-06,2.3755825e-06,   &
      1.9986249e-06,1.7985548e-06,2.1711692e-06,4.9318990e-06,1.4986332e-05,   &
      3.9020775e-05,7.7320255e-05,1.1791794e-04,1.4785099e-04,1.7224469e-04,   &
      2.2270834e-04,3.6641079e-04 ]
scav_coeff_num(1:ncol, 7, 1) = [                                               &
      8.5349386e-05,4.8801170e-05,2.8572065e-05,1.7155635e-05,1.0649567e-05,   &
      6.9594621e-06,4.9177858e-06,3.8337627e-06,3.2584061e-06,2.8725172e-06,   &
      2.4706447e-06,2.0180730e-06,1.6614375e-06,1.5831273e-06,1.8789634e-06,   &
      5.0497533e-06,8.6796719e-05,1.9905849e-04,2.1139246e-04,1.9220948e-04,   &
      1.9016876e-04,2.4168060e-04 ]
scav_coeff_num(1:ncol, 7, 2) = [                                               &
      9.1274751e-05,5.1894385e-05,3.0241546e-05,1.8085401e-05,1.1182516e-05,   &
      7.2713001e-06,5.0988747e-06,3.9299606e-06,3.2943883e-06,2.8698488e-06,   &
      2.4613987e-06,2.0362388e-06,1.7138264e-06,1.6608087e-06,2.7184571e-06,   &
      1.8910916e-05,9.2383956e-05,1.8080566e-04,2.0544199e-04,1.9541800e-04,   &
      2.0033145e-04,2.6445765e-04 ]
scav_coeff_num(1:ncol, 7, 3) = [                                               &
      1.0016785e-04,5.6487890e-05,3.2697233e-05,1.9441290e-05,1.1953334e-05,   &
      7.7185407e-06,5.3566189e-06,4.0670274e-06,3.3493907e-06,2.8757215e-06,   &
      2.4588177e-06,2.0629131e-06,1.7908008e-06,2.0950949e-06,6.8409496e-06,   &
      3.2545081e-05,9.5815048e-05,1.6573371e-04,1.9689849e-04,1.9915918e-04,   &
      2.1575299e-04,3.0076022e-04 ]
scav_coeff_num(1:ncol, 7, 4) = [                                               &
      1.1199026e-04,6.2512161e-05,3.5878100e-05,2.1177951e-05,1.2930251e-05,   &
      8.2795754e-06,5.6772501e-06,4.2382785e-06,3.4236941e-06,2.8952142e-06,   &
      2.4676256e-06,2.1065494e-06,2.0258354e-06,3.6308759e-06,1.2983474e-05,   &
      4.2720770e-05,9.7511966e-05,1.5534644e-04,1.8903482e-04,2.0358889e-04,   &
      2.3670731e-04,3.5261136e-04 ]
scav_coeff_num(1:ncol, 7, 5) = [                                               &
      1.2704389e-04,7.0063927e-05,3.9807744e-05,2.3295031e-05,1.4106594e-05,   &
      8.9474846e-06,6.0558363e-06,4.4418341e-06,3.5184402e-06,2.9309969e-06,   &
      2.4951517e-06,2.2256544e-06,2.6739612e-06,6.2467770e-06,1.9274675e-05,   &
      5.0069398e-05,9.8363493e-05,1.4842087e-04,1.8373163e-04,2.1015577e-04,   &
      2.6461554e-04,4.2392017e-04 ]
scav_coeff_num(1:ncol, 8, 1) = [                                               &
      9.1692652e-05,5.2565763e-05,3.0836030e-05,1.8566434e-05,1.1605351e-05,   &
      7.7145577e-06,5.6375551e-06,4.6150183e-06,4.1277044e-06,3.7772409e-06,   &
      3.2907362e-06,2.6324678e-06,2.0428613e-06,1.8150975e-06,2.0803621e-06,   &
      6.4140491e-06,1.1423577e-04,2.5005939e-04,2.6284412e-04,2.3731865e-04,   &
      2.3066778e-04,2.8479056e-04 ]
scav_coeff_num(1:ncol, 8, 2) = [                                               &
      9.7959062e-05,5.5862862e-05,3.2627972e-05,1.9571493e-05,1.2186370e-05,   &
      8.0576191e-06,5.8368973e-06,4.7162114e-06,4.1538019e-06,3.7520708e-06,   &
      3.2588204e-06,2.6455182e-06,2.1108923e-06,1.9182494e-06,3.2365571e-06,   &
      2.4774448e-05,1.1907092e-04,2.2765480e-04,2.5559325e-04,2.4084811e-04,   &
      2.4219796e-04,3.1049170e-04 ]
scav_coeff_num(1:ncol, 8, 3) = [                                               &
      1.0733710e-04,6.0747602e-05,3.5258836e-05,2.1034759e-05,1.3025039e-05,   &
      8.5480825e-06,6.1191638e-06,4.8595729e-06,4.1962247e-06,3.7312988e-06,   &
      3.2308151e-06,2.6662951e-06,2.2092338e-06,2.5043317e-06,8.7008471e-06,   &
      4.2064557e-05,1.2232563e-04,2.0857181e-04,2.4502213e-04,2.4485513e-04,   &
      2.5961898e-04,3.5140014e-04 ]
scav_coeff_num(1:ncol, 8, 4) = [                                               &
      1.1975766e-04,6.7133588e-05,3.8657828e-05,2.2904725e-05,1.4085254e-05,   &
      9.1610648e-06,6.4683941e-06,5.0380189e-06,4.2575838e-06,3.7246743e-06,   &
      3.2158018e-06,2.7076739e-06,2.5167911e-06,4.5317114e-06,1.6633366e-05,   &
      5.4745580e-05,1.2382066e-04,1.9518872e-04,2.3498780e-04,2.4945761e-04,   &
      2.8318570e-04,4.0974864e-04 ]
scav_coeff_num(1:ncol, 8, 5) = [                                               &
      1.3550362e-04,7.5108131e-05,4.2843469e-05,2.5177929e-05,1.5358143e-05,   &
      9.8880078e-06,6.8787322e-06,5.2499369e-06,4.3406813e-06,3.7372008e-06,   &
      3.2247186e-06,2.8480012e-06,3.3622742e-06,7.9250703e-06,2.4642273e-05,   &
      6.3812480e-05,1.2443360e-04,1.8608663e-04,2.2784961e-04,2.5630830e-04,   &
      3.1445892e-04,4.8989476e-04 ]
scav_coeff_num(1:ncol, 9, 1) = [                                               &
      9.8792116e-05,5.6779523e-05,3.3377076e-05,2.0165340e-05,1.2716829e-05,   &
      8.6363299e-06,6.5698626e-06,5.6782247e-06,5.3515253e-06,5.0802428e-06,   &
      4.4897788e-06,3.5356410e-06,2.5934370e-06,2.1253227e-06,2.3191693e-06,   &
      8.1326825e-06,1.4811269e-04,3.1242197e-04,3.2615384e-04,2.9277995e-04,   &
      2.7983981e-04,3.3563328e-04 ]
scav_coeff_num(1:ncol, 9, 2) = [                                               &
      1.0544338e-04,6.0306044e-05,3.5307202e-05,2.1256386e-05,1.3354109e-05,   &
      9.0171855e-06,6.7918806e-06,5.7849874e-06,5.3623067e-06,5.0207332e-06,   &
      4.4219177e-06,3.5380478e-06,2.6824747e-06,2.2642372e-06,3.8794032e-06,   &
      3.2061088e-05,1.5188904e-04,2.8506476e-04,3.1722377e-04,2.9659618e-04,   &
      2.9287626e-04,3.6452592e-04 ]
scav_coeff_num(1:ncol, 9, 3) = [                                               &
      1.1536992e-04,6.5518961e-05,3.8135858e-05,2.2842187e-05,1.4271925e-05,   &
      9.5595652e-06,7.1041255e-06,5.9347961e-06,5.3847496e-06,4.9588412e-06,   &
      4.3532545e-06,3.5464027e-06,2.8086175e-06,3.0451035e-06,1.1013981e-05,   &
      5.3813812e-05,1.5489835e-04,2.6113051e-04,3.0410095e-04,3.0077404e-04,   &
      3.1247105e-04,4.1044123e-04 ]
scav_coeff_num(1:ncol, 9, 4) = [                                               &
      1.2846932e-04,7.2313373e-05,4.1781353e-05,2.4864150e-05,1.5428915e-05,   &
      1.0234360e-05,7.4876356e-06,6.1199531e-06,5.4248078e-06,4.9113047e-06,   &
      4.2996570e-06,3.5798139e-06,3.2072117e-06,5.6867988e-06,2.1154166e-05,   &
      6.9553198e-05,1.5613987e-04,2.4408862e-04,2.9133725e-04,3.0539580e-04,   &
      3.3884091e-04,4.7582175e-04 ]
scav_coeff_num(1:ncol, 9, 5) = [                                               &
      1.4500471e-04,8.0766615e-05,4.6256664e-05,2.7315238e-05,1.6813496e-05,   &
      1.1030852e-05,7.9352630e-06,6.3392103e-06,5.4877231e-06,4.8866386e-06,   &
      4.2770719e-06,3.7415765e-06,4.2974550e-06,1.0043087e-05,3.1275782e-05,   &
      8.0721629e-05,1.5645605e-04,2.3230325e-04,2.8185934e-04,3.1235364e-04,   &
      3.7369175e-04,5.6549571e-04 ]
scav_coeff_num(1:ncol,10, 1) = [                                               &
      1.0668249e-04,6.1463378e-05,3.6209302e-05,2.1964856e-05,1.4000026e-05,   &
      9.7502489e-06,7.7563596e-06,7.0873948e-06,7.0187546e-06,6.8895537e-06,   &
      6.1778736e-06,4.8166411e-06,3.3673735e-06,2.5379705e-06,2.6037438e-06,   &
      1.0225711e-05,1.8917582e-04,3.8818530e-04,4.0370676e-04,3.6073380e-04,   &
      3.3941444e-04,3.9553861e-04 ]
scav_coeff_num(1:ncol,10, 2) = [                                               &
      1.1376453e-04,6.5246178e-05,3.8294328e-05,2.3153548e-05,1.4702876e-05,   &
      1.0176736e-05,8.0066347e-06,7.2007661e-06,7.0077040e-06,6.7806597e-06,   &
      6.0564580e-06,4.8002906e-06,3.4836746e-06,2.7256515e-06,4.6633233e-06,   &
      4.0916431e-05,1.9169289e-04,3.5490103e-04,3.9260766e-04,3.6477089e-04,   &
      3.5409764e-04,4.2789103e-04 ]
scav_coeff_num(1:ncol,10, 3) = [                                               &
      1.2430616e-04,7.0826121e-05,4.1344811e-05,2.4878396e-05,1.5712638e-05,   &
      1.0781354e-05,8.3557008e-06,7.3575511e-06,7.0011478e-06,6.6591524e-06,   &
      5.9268970e-06,4.7867136e-06,3.6448548e-06,3.7500306e-06,1.3835454e-05,   &
      6.8077718e-05,1.9445813e-04,3.2513199e-04,3.7626717e-04,3.6896890e-04,   &
      3.7603260e-04,4.7920819e-04 ]
scav_coeff_num(1:ncol,10, 4) = [                                               &
      1.3816896e-04,7.8078088e-05,4.5266966e-05,2.7072683e-05,1.6981607e-05,   &
      1.1529595e-05,8.7805359e-06,7.5490731e-06,7.0095688e-06,6.5517014e-06,   &
      5.8150636e-06,4.8035718e-06,4.1547858e-06,7.1444065e-06,2.6662471e-05,   &
      8.7540722e-05,1.9543548e-04,3.0366251e-04,3.6007800e-04,3.7338194e-04,   &
      4.0537294e-04,5.5213555e-04 ]
scav_coeff_num(1:ncol,10, 5) = [                                               &
      1.5559547e-04,8.7068839e-05,5.0067714e-05,2.9725266e-05,1.8494879e-05,   &
      1.2407918e-05,9.2722022e-06,7.7745203e-06,7.0416074e-06,6.4718497e-06,   &
      5.7436571e-06,4.9849042e-06,5.5425265e-06,1.2675454e-05,3.9357848e-05,   &
      1.0127903e-04,1.9541929e-04,2.8860595e-04,3.4763604e-04,3.8018721e-04,   &
      4.4397360e-04,6.5198786e-04 ]
scav_coeff_num(1:ncol,11, 1) = [                                               &
      1.1539352e-04,6.6634868e-05,3.9343653e-05,2.3973191e-05,1.5463884e-05,   &
      1.1069993e-05,9.2203326e-06,8.8805685e-06,9.1855733e-06,9.2775361e-06,   &
      8.4328677e-06,6.5416777e-06,4.4069240e-06,3.0727900e-06,2.9413125e-06,   &
      1.2683254e-05,2.3801839e-04,4.7956381e-04,4.9822205e-04,4.4366336e-04,   &
      4.1140041e-04,4.6602559e-04 ]
scav_coeff_num(1:ncol,11, 2) = [                                               &
      1.2295401e-04,7.0701863e-05,4.1601022e-05,2.5271844e-05,1.6242368e-05,   &
      1.1550867e-05,9.5054119e-06,9.0022191e-06,9.1458000e-06,9.1022030e-06,   &
      8.2370802e-06,6.4961469e-06,4.5570119e-06,3.3240749e-06,5.5976117e-06,   &
      5.1430260e-05,2.3927559e-04,4.3917425e-04,4.8432479e-04,4.4780997e-04,   &
      4.2786747e-04,5.0210005e-04 ]
scav_coeff_num(1:ncol,11, 3) = [                                               &
      1.3417991e-04,7.6689226e-05,4.4898431e-05,2.7153186e-05,1.7357909e-05,   &
      1.2229246e-05,9.8993285e-06,9.1671632e-06,9.1005431e-06,8.9000127e-06,   &
      8.0226658e-06,6.4484987e-06,4.7605542e-06,4.6437911e-06,1.7202965e-05,   &
      8.5096620e-05,2.4191831e-04,4.0243329e-04,4.6392566e-04,4.5180174e-04,   &
      4.5229115e-04,5.5919518e-04 ]
scav_coeff_num(1:ncol,11, 4) = [                                               &
      1.4889389e-04,8.4449836e-05,4.9128745e-05,2.9541275e-05,1.8755280e-05,   &
      1.3063894e-05,1.0373743e-05,9.3652155e-06,9.0660440e-06,8.7109080e-06,   &
      7.8293165e-06,6.4376866e-06,5.4027459e-06,8.9406897e-06,3.3253740e-05,   &
      1.0907787e-04,2.4269200e-04,3.7565957e-04,4.4344975e-04,4.5567988e-04,   &
      4.8473041e-04,6.4014377e-04 ]
scav_coeff_num(1:ncol,11, 5) = [                                               &
      1.6731714e-04,9.4039258e-05,5.4292292e-05,3.2420295e-05,2.0415574e-05,   &
      1.4037694e-05,1.0917321e-05,9.5960199e-06,9.0553942e-06,8.5552073e-06,   &
      7.6884774e-06,6.6348011e-06,7.1444391e-06,1.5881557e-05,4.9051677e-05,   &
      1.2595532e-04,2.4235011e-04,3.5666206e-04,4.2727856e-04,4.6196571e-04,   &
      5.2719166e-04,7.5075406e-04 ]
scav_coeff_num(1:ncol,12, 1) = [                                               &
      1.2494917e-04,7.2307336e-05,4.2786618e-05,2.6191647e-05,1.7104993e-05,   &
      1.2587535e-05,1.0949735e-05,1.1043588e-05,1.1839291e-05,1.2237397e-05,   &
      1.1256558e-05,8.7191548e-06,5.7214585e-06,3.7361222e-06,3.3352825e-06,   &
      1.5459683e-05,2.9499600e-04,5.8890539e-04,6.1276249e-04,5.4442619e-04,   &
      4.9811865e-04,5.4883058e-04 ]
scav_coeff_num(1:ncol,12, 2) = [                                               &
      1.3303718e-04,7.6687159e-05,4.5234092e-05,2.7612668e-05,1.7969195e-05,   &
      1.3131666e-05,1.1276526e-05,1.1175861e-05,1.1764734e-05,1.1979072e-05,   &
      1.0965209e-05,8.6329653e-06,5.9111077e-06,4.0657348e-06,6.6797894e-06,   &
      6.3616477e-05,2.9531134e-04,5.3999983e-04,5.9526487e-04,5.4850804e-04,   &
      5.1649700e-04,5.8887651e-04 ]
scav_coeff_num(1:ncol,12, 3) = [                                               &
      1.4501849e-04,8.3123375e-05,4.8804002e-05,2.9668119e-05,1.9204426e-05,   &
      1.3895527e-05,1.1723826e-05,1.1351006e-05,1.1672051e-05,1.1675653e-05,   &
      1.0641543e-05,8.5379665e-06,6.1633677e-06,5.7314359e-06,2.1127042e-05,   &
      1.0504138e-04,2.9813606e-04,4.9499265e-04,5.6975328e-04,5.5196652e-04,   &
      5.4352825e-04,6.5209409e-04 ]
scav_coeff_num(1:ncol,12, 4) = [                                               &
      1.6067398e-04,9.1445044e-05,5.3374608e-05,3.2271756e-05,2.0746774e-05,   &
      1.4829829e-05,1.2256670e-05,1.1556672e-05,1.1584304e-05,1.1383610e-05,   &
      1.0342967e-05,8.4872951e-06,6.9578161e-06,1.1085603e-05,4.0987899e-05,   &
      1.3448111e-04,2.9887620e-04,4.6193416e-04,5.4393416e-04,5.5485718e-04,   &
      5.7913795e-04,7.4147685e-04 ]
scav_coeff_num(1:ncol,12, 5) = [                                               &
      1.8020256e-04,1.0169590e-04,5.8939114e-05,3.5402528e-05,2.2572669e-05,   &
      1.5913117e-05,1.2860678e-05,1.1792902e-05,1.1520021e-05,1.1131773e-05,   &
      1.0111698e-05,8.6956159e-06,9.1114355e-06,1.9688920e-05,6.0483434e-05,   &
      1.5518415e-04,2.9828283e-04,4.3824814e-04,5.2310552e-04,5.6012184e-04,   &
      6.2548539e-04,8.6332247e-04 ]
scav_coeff_num(1:ncol,13, 1) = [                                               &
      1.3536731e-04,7.8489564e-05,4.6539544e-05,2.8613090e-05,1.8904302e-05,   &
      1.4266675e-05,1.2885652e-05,1.3491601e-05,1.4871935e-05,1.5650334e-05,   &
      1.4540484e-05,1.1271453e-05,7.2700341e-06,4.5131877e-06,3.7828833e-06,   &
      1.8472238e-05,3.6012934e-04,7.1861299e-04,7.5071694e-04,6.6626916e-04,   &
      6.0222583e-04,6.4593161e-04 ]
scav_coeff_num(1:ncol,13, 2) = [                                               &
      1.4403284e-04,8.3211160e-05,4.9194688e-05,3.0168240e-05,1.9863322e-05,   &
      1.4881888e-05,1.3260462e-05,1.3637377e-05,1.4758848e-05,1.5296439e-05,   &
      1.4136228e-05,1.1134224e-05,7.5030554e-06,4.9329737e-06,7.8918994e-06,   &
      7.7394984e-05,3.6028370e-04,6.5952598e-04,7.2860482e-04,6.7002698e-04,   &
      6.2262547e-04,6.9017947e-04 ]
scav_coeff_num(1:ncol,13, 3) = [                                               &
      1.5684208e-04,9.0138123e-05,5.3062429e-05,3.2414579e-05,2.1230911e-05,   &
      1.5741695e-05,1.3769075e-05,1.3825672e-05,1.4613087e-05,1.4876071e-05,   &
      1.3683335e-05,1.0979979e-05,7.8082223e-06,6.9886410e-06,2.5582449e-05,   &
      1.2798637e-04,3.6385330e-04,6.0480933e-04,6.9667277e-04,6.7249432e-04,   &
      6.5233931e-04,7.5982034e-04 ]
scav_coeff_num(1:ncol,13, 4) = [                                               &
      1.7353108e-04,9.9073860e-05,5.8005215e-05,3.5254602e-05,2.2933504e-05,   &
      1.6787657e-05,1.4368765e-05,1.4041276e-05,1.4465178e-05,1.4464727e-05,   &
      1.3260146e-05,1.0878654e-05,8.7722536e-06,1.3551026e-05,4.9875366e-05,   &
      1.6398331e-04,3.6488636e-04,5.6439525e-04,6.6422828e-04,6.7378579e-04,   &
      6.9111357e-04,8.5796880e-04 ]
scav_coeff_num(1:ncol,13, 5) = [                                               &
      1.9427565e-04,1.1004966e-04,6.4008637e-05,3.8661549e-05,2.4942314e-05,   &
      1.7993341e-05,1.5041527e-05,1.4284490e-05,1.4339786e-05,1.4101128e-05,   &
      1.2921397e-05,1.1094823e-05,1.1394633e-05,2.4078746e-05,7.3724120e-05,   &
      1.8933088e-04,3.6421464e-04,5.3520643e-04,6.3763001e-04,6.7736206e-04,   &
      7.4126190e-04,9.9139497e-04 ]
scav_coeff_num(1:ncol,14, 1) = [                                               &
      1.4665935e-04,8.5185650e-05,5.0598777e-05,3.1222408e-05,2.0827933e-05,   &
      1.6044155e-05,1.4923451e-05,1.6069650e-05,1.8079433e-05,1.9282590e-05,   &
      1.8060954e-05,1.4029343e-05,8.9570951e-06,5.3658472e-06,4.2744094e-06,   &
      2.1604734e-05,4.3298373e-04,8.7099626e-04,9.1572173e-04,8.1279857e-04,   &
      7.2670573e-04,7.5954936e-04 ]
scav_coeff_num(1:ncol,14, 2) = [                                               &
      1.5595289e-04,9.0277877e-05,5.3478455e-05,3.2922068e-05,2.1888874e-05,   &
      1.6736023e-05,1.5350914e-05,1.6232054e-05,1.7927734e-05,1.8828007e-05,   &
      1.7534766e-05,1.3834301e-05,9.2342669e-06,5.8819119e-06,9.1995094e-06,   &
      9.2575177e-05,4.3438741e-04,7.9980056e-04,8.8772650e-04,8.1586101e-04,   &
      7.4920928e-04,8.0820503e-04 ]
scav_coeff_num(1:ncol,14, 3) = [                                               &
      1.6966326e-04,9.7737356e-05,5.7668398e-05,3.5374241e-05,2.3398910e-05,   &
      1.7699458e-05,1.5926910e-05,1.6437210e-05,1.7728268e-05,1.8284150e-05,   &
      1.6942272e-05,1.3613214e-05,9.5933486e-06,8.3590724e-06,3.0502614e-05,   &
      1.5387786e-04,4.3960973e-04,7.3381076e-04,8.4777032e-04,8.1671283e-04,   &
      7.8161980e-04,8.8451501e-04 ]
scav_coeff_num(1:ncol,14, 4) = [                                               &
      1.8747867e-04,1.0734003e-04,6.3014133e-05,3.8469410e-05,2.5274204e-05,   &
      1.8866220e-05,1.6600238e-05,1.6666449e-05,1.7519031e-05,1.7746613e-05,   &
      1.6384477e-05,1.3454668e-05,1.0739940e-05,1.6266619e-05,5.9864929e-05,   &
      1.9769470e-04,4.4147345e-04,6.8490878e-04,8.0716586e-04,8.1559767e-04,   &
      8.2345186e-04,9.9166039e-04 ]
scav_coeff_num(1:ncol,14, 5) = [                                               &
      2.0955089e-04,1.1910414e-04,6.9493229e-05,4.2174785e-05,2.7480431e-05,   &
      2.0204529e-05,1.7348882e-05,1.6920112e-05,1.7331041e-05,1.7264669e-05,   &
      1.5929752e-05,1.3678807e-05,1.3883738e-05,2.8979219e-05,8.8771825e-05,   &
      2.2865011e-04,4.4103307e-04,6.4935828e-04,7.7348573e-04,8.1662052e-04,   &
      8.7717702e-04,1.1368523e-03 ]
scav_coeff_num(1:ncol,15, 1) = [                                               &
      1.5882892e-04,9.2394684e-05,5.4956433e-05,3.3998998e-05,2.2832537e-05,   &
      1.7839485e-05,1.6928826e-05,1.8576160e-05,2.1192169e-05,2.2819966e-05,   &
      2.1511549e-05,1.6756139e-05,1.0645875e-05,6.2379302e-06,4.7947191e-06,   &
      2.4716191e-05,5.1252354e-04,1.0480100e-03,1.1114708e-03,9.8786049e-04,   &
      8.7479179e-04,8.9209396e-04 ]
scav_coeff_num(1:ncol,15, 2) = [                                               &
      1.6880098e-04,9.7885945e-05,5.8076427e-05,3.5851663e-05,2.3999793e-05,   &
      1.8610448e-05,1.7411059e-05,1.8758192e-05,2.1006204e-05,2.2269408e-05,   &
      2.0866015e-05,1.6502176e-05,1.0964374e-05,6.8485932e-06,1.0554259e-05,   &
      1.0884124e-04,5.1739297e-04,9.6253946e-04,1.0760272e-03,9.8971173e-04,   &
      8.9944117e-04,9.4533247e-04 ]
scav_coeff_num(1:ncol,15, 3) = [                                               &
      1.8348569e-04,1.0591904e-04,6.2611430e-05,3.8522065e-05,2.5658837e-05,   &
      1.9681283e-05,1.8057645e-05,1.8984381e-05,2.0758593e-05,2.1607689e-05,   &
      2.0136958e-05,1.6213322e-05,1.1374678e-05,9.7613540e-06,3.5777963e-05,   &
      1.8249535e-04,5.2561355e-04,8.8365361e-04,1.0261126e-03,9.8811679e-04,   &
      9.3448027e-04,1.0284904e-03 ]
scav_coeff_num(1:ncol,15, 4) = [                                               &
      2.0252034e-04,1.1624069e-04,6.8389078e-05,4.1888218e-05,2.7715403e-05,   &
      2.0973918e-05,1.8808799e-05,1.9232203e-05,2.0494048e-05,2.0949673e-05,   &
      1.9447680e-05,1.5997472e-05,1.2711058e-05,1.9125985e-05,7.0833831e-05,   &
      2.3555018e-04,5.2911775e-04,8.2512184e-04,9.7554533e-04,9.8355336e-04,   &
      9.7913520e-04,1.1447442e-03 ]
scav_coeff_num(1:ncol,15, 5) = [                                               &
      2.2603182e-04,1.2885552e-04,7.5378608e-05,4.5911166e-05,3.0129600e-05,   &
      2.2451216e-05,1.9638398e-05,1.9501771e-05,2.0249545e-05,2.0354880e-05,   &
      1.8880410e-05,1.6234823e-05,1.6421328e-05,3.4268984e-05,1.0553376e-04,   &
      2.7322519e-04,5.2939817e-04,7.8234491e-04,9.3326632e-04,9.8092890e-04,   &
      1.0360442e-03,1.3016995e-03 ]
scav_coeff_num(1:ncol,16, 1) = [                                               &
      1.7186835e-04,1.0010923e-04,5.9600946e-05,3.6920208e-05,2.4873424e-05,   &
      1.9570507e-05,1.8764058e-05,2.0802573e-05,2.3928196e-05,2.5930042e-05,   &
      2.4564003e-05,1.9194901e-05,1.2185828e-05,7.0666815e-06,5.3265629e-06,   &
      2.7653628e-05,5.9695781e-04,1.2508364e-03,1.3413472e-03,1.1952641e-03,   &
      1.0497638e-03,1.0460104e-03 ]
scav_coeff_num(1:ncol,16, 2) = [                                               &
      1.8256876e-04,1.0602710e-04,6.2975722e-05,3.8932305e-05,2.6148484e-05,   &
      2.0419606e-05,1.9300318e-05,2.1006777e-05,2.3716636e-05,2.5298478e-05,   &
      2.3814298e-05,1.8887488e-05,1.2539366e-05,7.7619342e-06,1.1899437e-05,   &
      1.2574140e-04,6.0847266e-04,1.1487602e-03,1.2965602e-03,1.1952081e-03,   &
      1.0765423e-03,1.1039675e-03 ]
scav_coeff_num(1:ncol,16, 3) = [                                               &
      1.9830004e-04,1.1467371e-04,6.7876796e-05,4.1830528e-05,2.7959290e-05,   &
      2.1597345e-05,2.0017400e-05,2.1258246e-05,2.3432576e-05,2.4537268e-05,   &
      2.2966267e-05,1.8537637e-05,1.2994640e-05,1.1103748e-05,4.1257910e-05,   &
      2.1340643e-04,6.2156102e-04,1.0554041e-03,1.2344059e-03,1.1900867e-03,   &
      1.1140333e-03,1.1940699e-03 ]
scav_coeff_num(1:ncol,16, 4) = [                                               &
      2.1864553e-04,1.2576485e-04,7.4113087e-05,4.5480246e-05,3.0201469e-05,   &
      2.3016439e-05,2.0847500e-05,2.1530497e-05,2.3125963e-05,2.3777738e-05,   &
      2.2162856e-05,1.8271798e-05,1.4521280e-05,2.2001113e-05,8.2579368e-05,   &
      2.7724084e-04,6.2784800e-04,9.8617489e-04,1.1718121e-03,1.1807641e-03,   &
      1.1611168e-03,1.3194014e-03 ]
scav_coeff_num(1:ncol,16, 5) = [                                               &
      2.4370633e-04,1.3929101e-04,8.1645293e-05,4.9836410e-05,3.2829856e-05,   &
      2.4634787e-05,2.1760622e-05,2.1823121e-05,2.2838706e-05,2.3088318e-05,   &
      2.1499368e-05,1.8533009e-05,1.8831595e-05,3.9788512e-05,1.2380637e-04,   &
      3.2288479e-04,6.2956374e-04,9.3536230e-04,1.1192288e-03,1.1731456e-03,   &
      1.2206158e-03,1.4879005e-03 ]
scav_coeff_num(1:ncol,17, 1) = [                                               &
      1.8574906e-04,1.0831034e-04,6.4515524e-05,3.9963546e-05,2.6912086e-05,   &
      2.1169005e-05,2.0315226e-05,2.2575256e-05,2.6050481e-05,2.8330224e-05,   &
      2.6935913e-05,2.1121843e-05,1.3444200e-05,7.7960982e-06,5.8543996e-06,   &
      3.0267476e-05,6.8361873e-04,1.4792784e-03,1.6077858e-03,1.4382459e-03,   &
      1.2545275e-03,1.2234476e-03 ]
scav_coeff_num(1:ncol,17, 2) = [                                               &
      1.9722585e-04,1.1468097e-04,6.8158062e-05,4.2139555e-05,2.8293859e-05,   &
      2.2092274e-05,2.0902144e-05,2.2803501e-05,2.5825369e-05,2.7641216e-05,   &
      2.6108162e-05,2.0772531e-05,1.3823830e-05,8.5587295e-06,1.3176970e-05,   &
      1.4268626e-04,7.0600414e-04,1.3582534e-03,1.5514234e-03,1.4353734e-03,   &
      1.2833372e-03,1.2862062e-03 ]
scav_coeff_num(1:ncol,17, 3) = [                                               &
      2.1407326e-04,1.2397895e-04,7.3444102e-05,4.5272521e-05,3.0255772e-05,   &
      2.3372615e-05,2.1686563e-05,2.3084202e-05,2.5521418e-05,2.6809602e-05,   &
      2.5171408e-05,2.0375528e-05,1.4315017e-05,1.2301304e-05,4.6755695e-05,   &
      2.4592078e-04,7.2640836e-04,1.2490684e-03,1.4744255e-03,1.4253645e-03,   &
      1.3229640e-03,1.3832461e-03 ]
scav_coeff_num(1:ncol,17, 4) = [                                               &
      2.3581777e-04,1.3588750e-04,8.0163166e-05,4.9215243e-05,3.2684074e-05,   &
      2.4914713e-05,2.2593812e-05,2.3387101e-05,2.5191770e-05,2.5978704e-05,   &
      2.4283597e-05,2.0073672e-05,1.6025033e-05,2.4759716e-05,9.4811544e-05,   &
      3.2213236e-04,7.3699805e-04,1.1682699e-03,1.3975289e-03,1.4096827e-03,   &
      1.3718897e-03,1.5174536e-03 ]
scav_coeff_num(1:ncol,17, 5) = [                                               &
      2.6253385e-04,1.5038253e-04,8.8267299e-05,5.3916816e-05,3.5528883e-05,   &
      2.6672232e-05,2.3590563e-05,2.3710993e-05,2.4881737e-05,2.5223633e-05,   &
      2.3551357e-05,2.0373495e-05,2.0953136e-05,4.5353892e-05,1.4325229e-04,   &
      3.7709579e-04,7.4112679e-04,1.1087561e-03,1.3327998e-03,1.3954653e-03,   &
      1.4331547e-03,1.6970246e-03 ]
scav_coeff_num(1:ncol,18, 1) = [                                               &
      2.0039841e-04,1.1695429e-04,6.9671293e-05,4.3104732e-05,2.8919381e-05,   &
      2.2590542e-05,2.1511092e-05,2.3785801e-05,2.7409312e-05,2.9839122e-05,   &
      2.8442731e-05,2.2387971e-05,1.4330900e-05,8.3873659e-06,6.3670326e-06,   &
      3.2427266e-05,7.6893481e-04,1.7309666e-03,1.9112638e-03,1.7185348e-03,   &
      1.4908465e-03,1.4256391e-03 ]
scav_coeff_num(1:ncol,18, 2) = [                                               &
      2.1269528e-04,1.2380119e-04,7.3592635e-05,4.5447347e-05,3.0404810e-05,   &
      2.3581850e-05,2.2143358e-05,2.4039194e-05,2.7184516e-05,2.9121588e-05,   &
      2.7570317e-05,2.2012692e-05,1.4726258e-05,9.1953700e-06,1.4333805e-05,   &
      1.5896277e-04,8.0738796e-04,1.5888883e-03,1.8408059e-03,1.7117077e-03,   &
      1.5214820e-03,1.4932076e-03 ]
scav_coeff_num(1:ncol,18, 3) = [                                               &
      2.3072210e-04,1.3378451e-04,7.9279754e-05,4.8819472e-05,3.2514521e-05,   &
      2.4957701e-05,2.2989616e-05,2.4352456e-05,2.6879946e-05,2.8255381e-05,   &
      2.6583482e-05,2.1587431e-05,1.5242636e-05,1.3289646e-05,5.2055464e-05,   &
      2.7905428e-04,8.3811542e-04,1.4629569e-03,1.7461338e-03,1.6951643e-03,   &
      1.5627590e-03,1.5970430e-03 ]
scav_coeff_num(1:ncol,18, 4) = [                                               &
      2.5394583e-04,1.4655349e-04,8.6502242e-05,5.3061605e-05,3.5126394e-05,   &
      2.6616310e-05,2.3969910e-05,2.4691997e-05,2.6549615e-05,2.7390304e-05,   &
      2.5649197e-05,2.1268091e-05,1.7121554e-05,2.7279679e-05,1.0714773e-04,   &
      3.6917994e-04,8.5490872e-04,1.3700704e-03,1.6525605e-03,1.6712508e-03,   &
      1.6127215e-03,1.7397159e-03 ]
scav_coeff_num(1:ncol,18, 5) = [                                               &
      2.8241354e-04,1.6206905e-04,9.5203458e-05,5.8117341e-05,3.8186561e-05,   &
      2.8508065e-05,2.5048332e-05,2.5055574e-05,2.6240405e-05,2.6605346e-05,   &
      2.4882915e-05,2.1622805e-05,2.2664712e-05,5.0767648e-05,1.6337550e-04,   &
      4.3483856e-04,8.6270419e-04,1.3014473e-03,1.5738129e-03,1.6486048e-03,   &
      1.6746791e-03,1.9295927e-03 ]
scav_coeff_num(1:ncol,19, 1) = [                                               &
      2.1565206e-04,1.2594481e-04,7.5011083e-05,4.6308676e-05,3.0871541e-05,   &
      2.3814869e-05,2.2328266e-05,2.4401710e-05,2.7958895e-05,3.0397937e-05,   &
      2.9020409e-05,2.2937800e-05,1.4809764e-05,8.8231511e-06,6.8577445e-06,   &
      3.4034502e-05,8.4856290e-04,2.0004484e-03,2.2488616e-03,2.0349045e-03,   &
      1.7581057e-03,1.6518904e-03 ]
scav_coeff_num(1:ncol,19, 2) = [                                               &
      2.2880313e-04,1.3328599e-04,7.9219048e-05,4.8818556e-05,3.2456094e-05,   &
      2.4866828e-05,2.2999444e-05,2.4680571e-05,2.7748340e-05,2.9682253e-05,   &
      2.8139299e-05,2.2554499e-05,1.5210396e-05,9.6526998e-06,1.5325956e-05,   &
      1.7376963e-04,9.0893417e-04,1.8358031e-03,2.1616463e-03,2.0227819e-03,   &
      1.7902314e-03,1.7241690e-03 ]
scav_coeff_num(1:ncol,19, 3) = [                                               &
      2.4805874e-04,1.4398070e-04,8.5318697e-05,5.2431327e-05,3.4708206e-05,   &
      2.6329380e-05,2.3900596e-05,2.5028817e-05,2.7462663e-05,2.8819143e-05,   &
      2.7144022e-05,2.2122147e-05,1.5740800e-05,1.4031497e-05,5.6921595e-05,   &
      3.1151798e-04,9.5340346e-04,1.6929133e-03,2.0464519e-03,1.9978273e-03,   &
      1.8324833e-03,1.8344783e-03 ]
scav_coeff_num(1:ncol,19, 4) = [                                               &
      2.7282423e-04,1.5764308e-04,9.3059370e-05,5.6975603e-05,3.7498531e-05,   &
      2.8095986e-05,2.4948240e-05,2.5410125e-05,2.7154415e-05,2.7959077e-05,   &
      2.6204085e-05,2.1805588e-05,1.7766656e-05,2.9457124e-05,1.1911238e-04,   &
      4.1685736e-04,9.7860633e-04,1.5879501e-03,1.9339383e-03,1.9636227e-03,   &
      1.8824510e-03,1.9849474e-03 ]
scav_coeff_num(1:ncol,19, 5) = [                                               &
      3.0311868e-04,1.7421877e-04,1.0237578e-04,6.2389900e-05,4.0770022e-05,   &
      3.0114815e-05,2.6104743e-05,2.5821041e-05,2.6870161e-05,2.7181926e-05,   &
      2.5440804e-05,2.2229841e-05,2.3897132e-05,5.5824087e-05,1.8349985e-04,   &
      4.9448159e-04,9.9155929e-04,1.5101947e-03,1.8394457e-03,1.9305921e-03,   &
      1.9437861e-03,2.1840186e-03 ]
scav_coeff_num(1:ncol,20, 1) = [                                               &
      2.3117793e-04,1.3508827e-04,8.0422520e-05,4.9512619e-05,3.2738657e-05,   &
      2.4836794e-05,2.2782615e-05,2.4457380e-05,2.7747804e-05,3.0061125e-05,   &
      2.8717522e-05,2.2803860e-05,1.4895159e-05,9.1051123e-06,7.3215377e-06,   &
      3.5029800e-05,9.1771641e-04,2.2783716e-03,2.6125524e-03,2.3813125e-03,   &
      2.0516697e-03,1.8982161e-03 ]
scav_coeff_num(1:ncol,20, 2) = [                                               &
      2.4519801e-04,1.4293091e-04,8.4918938e-05,5.2187248e-05,3.4416258e-05,   &
      2.5941320e-05,2.3485789e-05,2.4761207e-05,2.7563854e-05,2.9375459e-05,   &
      2.7861779e-05,2.2430083e-05,1.5291630e-05,9.9333383e-06,1.6119397e-05,   &
      1.8627595e-04,1.0058884e-03,2.0906514e-03,2.5060595e-03,2.3624462e-03,   &
      2.0848075e-03,1.9749486e-03 ]
scav_coeff_num(1:ncol,20, 3) = [                                               &
      2.6570381e-04,1.5434738e-04,9.1433893e-05,5.6037522e-05,3.6803062e-05,   &
      2.7480638e-05,2.4434073e-05,2.5145633e-05,2.7314507e-05,2.8550436e-05,   &
      2.6897474e-05,2.2011237e-05,1.5825474e-05,1.4515220e-05,6.1110834e-05,   &
      3.4174951e-04,1.0676047e-03,1.9315513e-03,2.3678290e-03,2.3271255e-03,   &
      2.1271744e-03,2.0911713e-03 ]
scav_coeff_num(1:ncol,20, 4) = [                                               &
      2.9203809e-04,1.6891629e-04,9.9696520e-05,6.0880755e-05,3.9763679e-05,   &
      2.9345145e-05,2.5542200e-05,2.5572283e-05,2.7048785e-05,2.7731786e-05,   &
      2.5990477e-05,2.1716619e-05,1.7969257e-05,3.1207893e-05,1.3014713e-04,   &
      4.6312908e-04,1.1035364e-03,1.8152148e-03,2.2345387e-03,2.2805703e-03,   &
      2.1759195e-03,2.2484487e-03 ]
scav_coeff_num(1:ncol,20, 5) = [                                               &
      3.2419160e-04,1.8656788e-04,1.0963305e-04,6.6650855e-05,4.3238688e-05,   &
      3.1481821e-05,2.6771672e-05,2.6036387e-05,2.6811205e-05,2.6997425e-05,   &
      2.5264825e-05,2.2220460e-05,2.4630624e-05,6.0310713e-05,2.0276195e-04,   &
      5.5369184e-04,1.1232561e-03,1.7288082e-03,2.1229777e-03,2.2352627e-03,   &
      2.2351268e-03,2.4552003e-03 ]
scav_coeff_num(1:ncol,21, 1) = [                                               &
      2.4638974e-04,1.4404260e-04,8.5706955e-05,5.2605837e-05,3.4469203e-05,   &
      2.5651284e-05,2.2911719e-05,2.4031742e-05,2.6891294e-05,2.8965362e-05,   &
      2.7665566e-05,2.2083678e-05,1.4637903e-05,9.2465126e-06,7.7504340e-06,   &
      3.5392958e-05,9.7167576e-04,2.5511446e-03,2.9877629e-03,2.7451491e-03,   &
      2.3612760e-03,2.1559896e-03 ]
scav_coeff_num(1:ncol,21, 2) = [                                               &
      2.6125961e-04,1.5237482e-04,9.0483154e-05,5.5437260e-05,3.6231367e-05,   &
      2.6799548e-05,2.3639719e-05,2.4359071e-05,2.6743502e-05,2.8332887e-05,   &
      2.6863855e-05,2.1734578e-05,1.5022481e-05,1.0053560e-05,1.6688278e-05,   &
      1.9569951e-04,1.0926712e-03,2.3412415e-03,2.8600408e-03,2.7181656e-03,   &
      2.3948154e-03,2.2367055e-03 ]
scav_coeff_num(1:ncol,21, 3) = [                                               &
      2.8298758e-04,1.6449579e-04,9.7400989e-05,5.9513955e-05,3.8741495e-05,   &
      2.8404341e-05,2.4626680e-05,2.4779250e-05,2.6544048e-05,2.7574741e-05,   &
      2.5963648e-05,2.1346922e-05,1.5550408e-05,1.4747501e-05,6.4388139e-05,   &
      3.6800367e-04,1.1747206e-03,2.1678144e-03,2.6970793e-03,2.6707087e-03,   &
      2.4362956e-03,2.3579748e-03 ]
scav_coeff_num(1:ncol,21, 4) = [                                               &
      3.1085491e-04,1.7994935e-04,1.0617008e-04,6.4642810e-05,4.1859375e-05,   &
      3.0354562e-05,2.5787122e-05,2.5252537e-05,2.6336809e-05,2.6827427e-05,   &
      2.5121929e-05,2.1090005e-05,1.7776913e-05,3.2465140e-05,1.3963823e-04,   &
      5.0550610e-04,1.2234962e-03,2.0416021e-03,2.5420090e-03,2.6100407e-03,   &
      2.4824760e-03,2.5206912e-03 ]
scav_coeff_num(1:ncol,21, 5) = [                                               &
      3.4482428e-04,1.9865062e-04,1.1670852e-04,7.0753653e-05,4.5523970e-05,   &
      3.2596925e-05,2.7082457e-05,2.5772823e-05,2.6162841e-05,2.6164510e-05,   &
      2.4462410e-05,2.1676707e-05,2.4882164e-05,6.4010661e-05,2.2013507e-04,   &
      6.0944219e-04,1.2515008e-03,1.9476079e-03,2.4127763e-03,2.5509116e-03,   &
      2.5379672e-03,2.7331538e-03 ]
scav_coeff_num(1:ncol,22, 1) = [                                               &
      2.6039575e-04,1.5228683e-04,9.0560439e-05,5.5416246e-05,3.5978089e-05,   &
      2.6239495e-05,2.2755823e-05,2.3221889e-05,2.5537079e-05,2.7290033e-05,   &
      2.6040320e-05,2.0909476e-05,1.4106944e-05,9.2633958e-06,8.1290462e-06,   &
      3.5137250e-05,1.0064044e-03,2.8015314e-03,3.3530810e-03,3.1065367e-03,   &
      2.6702933e-03,2.4113047e-03 ]
scav_coeff_num(1:ncol,22, 2) = [                                               &
      2.7604435e-04,1.6106764e-04,9.5591630e-05,5.8388012e-05,3.7812396e-05,   &
      2.7421300e-05,2.3501040e-05,2.3570065e-05,2.5431419e-05,2.6727295e-05,   &
      2.5313815e-05,2.0596502e-05,1.4473826e-05,1.0033517e-05,1.7012244e-05,   &
      2.0139560e-04,1.1633780e-03,2.5719962e-03,3.2032573e-03,3.0703877e-03,   &
      2.7035281e-03,2.4952641e-03 ]
scav_coeff_num(1:ncol,22, 3) = [                                               &
      2.9889130e-04,1.7383375e-04,1.0287660e-04,6.2667682e-05,4.0428580e-05,   &
      2.9077944e-05,2.4517128e-05,2.4023487e-05,2.5290427e-05,2.6056538e-05,   &
      2.4502081e-05,2.0253261e-05,1.4987955e-05,1.4744084e-05,6.6546797e-05,   &
      3.8851178e-04,1.2678141e-03,2.3872883e-03,3.0152277e-03,3.0095561e-03,   &
      2.7430632e-03,2.6203491e-03 ]
scav_coeff_num(1:ncol,22, 4) = [                                               &
      3.2816048e-04,1.9009674e-04,1.1210714e-04,6.8053065e-05,4.3682932e-05,   &
      3.1098244e-05,2.5719750e-05,2.4541712e-05,2.5151756e-05,2.5401769e-05,   &
      2.3749277e-05,2.0044944e-05,1.7257137e-05,3.3177060e-05,1.4696676e-04,   &
      5.4122571e-04,1.3309402e-03,2.2534896e-03,2.8386298e-03,2.9336818e-03,   &
      2.7853536e-03,2.7867097e-03 ]
scav_coeff_num(1:ncol,22, 5) = [                                               &
      3.6378751e-04,2.0975741e-04,1.2319388e-04,7.4470560e-05,4.7513526e-05,   &
      3.3429452e-05,2.7070889e-05,2.5117344e-05,2.5052013e-05,2.4830415e-05,   &
      2.3176014e-05,2.0709996e-05,2.4689649e-05,6.6712154e-05,2.3450105e-04,   &
      6.5818696e-04,1.3683810e-03,2.1535597e-03,2.6921600e-03,2.8598736e-03,   &
      2.8356157e-03,3.0024301e-03 ]

scav_coeff_mass(1:ncol, 1, 1) = [                                              &
      5.3745231e-05,3.0471292e-05,1.7799763e-05,1.0682883e-05,6.5863179e-06,   &
      4.1856180e-06,2.7625087e-06,1.9144385e-06,1.4094690e-06,1.1117903e-06,   &
      9.4386013e-07,8.6775625e-07,8.7377425e-07,9.7328374e-07,1.2046318e-06,   &
      2.0995472e-06,1.8703475e-05,4.8626027e-05,5.5368202e-05,5.4622025e-05,   &
      6.4438542e-05,1.0226403e-04 ]
scav_coeff_mass(1:ncol, 1, 2) = [                                              &
      4.2612415e-05,2.4375993e-05,1.4360801e-05,8.6965261e-06,5.4194431e-06,   &
      3.4929053e-06,2.3494812e-06,1.6686298e-06,1.2649820e-06,1.0313155e-06,   &
      9.0977269e-07,8.7656299e-07,9.3306178e-07,1.1239419e-06,2.3087080e-06,   &
      1.1354039e-05,3.4115185e-05,5.1223839e-05,5.5116574e-05,6.1034258e-05,   &
      8.7268862e-05,1.6307426e-04 ]
scav_coeff_mass(1:ncol, 1, 3) = [                                              &
      3.1305905e-05,1.8125496e-05,1.0807173e-05,6.6326575e-06,4.2030985e-06,   &
      2.7701101e-06,1.9191745e-06,1.4142874e-06,1.1196859e-06,9.6005398e-07,   &
      8.9997776e-07,9.3756281e-07,1.2074771e-06,2.9182819e-06,1.0455557e-05,   &
      2.7065133e-05,4.4472562e-05,5.3813904e-05,6.1775700e-05,8.6073684e-05,   &
      1.5637231e-04,3.3835058e-04 ]
scav_coeff_mass(1:ncol, 1, 4) = [                                              &
      2.2049785e-05,1.2951012e-05,7.8401186e-06,4.8995734e-06,3.1787770e-06,   &
      2.1616183e-06,1.5591140e-06,1.2068049e-06,1.0129911e-06,9.3394044e-07,   &
      9.8784714e-07,1.4688329e-06,3.7425976e-06,1.0822082e-05,2.4302812e-05,   &
      3.9804694e-05,5.1701430e-05,6.3501816e-05,9.1095100e-05,1.6668155e-04,   &
      3.6185721e-04,8.5068497e-04 ]
scav_coeff_mass(1:ncol, 1, 5) = [                                              &
      1.5205923e-05,9.0816803e-06,5.6033221e-06,3.5866887e-06,2.4019250e-06,   &
      1.7025239e-06,1.2939047e-06,1.0687080e-06,9.8299404e-07,1.1063105e-06,   &
      1.8864481e-06,4.6635416e-06,1.1540412e-05,2.3269516e-05,3.7223857e-05,   &
      5.0376600e-05,6.6239498e-05,1.0053311e-04,1.8959552e-04,4.1671403e-04,   &
      9.7631825e-04,2.2715296e-03 ]
scav_coeff_mass(1:ncol, 2, 1) = [                                              &
      5.6063579e-05,3.1847880e-05,1.8617014e-05,1.1173214e-05,6.8888188e-06,   &
      4.3829845e-06,2.9033928e-06,2.0262403e-06,1.5054641e-06,1.1957061e-06,   &
      1.0147457e-06,9.2645788e-07,9.2785140e-07,1.0367287e-06,1.2969745e-06,   &
      2.5796045e-06,2.6417692e-05,6.3026124e-05,6.9206335e-05,6.7067651e-05,   &
      7.7368874e-05,1.2006288e-04 ]
scav_coeff_mass(1:ncol, 2, 2) = [                                              &
      4.4471208e-05,2.5480453e-05,1.5019381e-05,9.0959547e-06,5.6712865e-06,   &
      3.6632230e-06,2.4764776e-06,1.7727172e-06,1.3547692e-06,1.1084590e-06,   &
      9.7460454e-07,9.3355668e-07,9.9316647e-07,1.2163740e-06,2.9000762e-06,   &
      1.5607006e-05,4.5154477e-05,6.5078995e-05,6.8267674e-05,7.3942370e-05,   &
      1.0321728e-04,1.8971105e-04 ]
scav_coeff_mass(1:ncol, 2, 3) = [                                              &
      3.2687488e-05,1.8948913e-05,1.1303007e-05,6.9398664e-06,4.4043401e-06,   &
      2.9136490e-06,2.0317622e-06,1.5088221e-06,1.2007138e-06,1.0291133e-06,   &
      9.6079036e-07,1.0019051e-06,1.3532496e-06,3.7446933e-06,1.4059691e-05,   &
      3.5717509e-05,5.6969283e-05,6.7119231e-05,7.5164192e-05,1.0214278e-04,   &
      1.8225454e-04,3.9112904e-04 ]
scav_coeff_mass(1:ncol, 2, 4) = [                                              &
      2.3032570e-05,1.3541404e-05,8.2024565e-06,5.1322591e-06,3.3396168e-06,   &
      2.2832421e-06,1.6582995e-06,1.2908383e-06,1.0851219e-06,9.9924430e-07,   &
      1.0711032e-06,1.7253685e-06,4.8472812e-06,1.4358597e-05,3.1841397e-05,   &
      5.1002928e-05,6.4651509e-05,7.7367913e-05,1.0814713e-04,1.9431944e-04,   &
      4.1837226e-04,9.8078995e-04 ]
scav_coeff_mass(1:ncol, 2, 5) = [                                              &
      1.5889885e-05,9.4993086e-06,5.8680113e-06,3.7654846e-06,2.5331363e-06,   &
      1.8066192e-06,1.3808048e-06,1.1438978e-06,1.0571559e-06,1.2323535e-06,   &
      2.2973392e-06,6.0542218e-06,1.5169707e-05,3.0273798e-05,4.7544866e-05,   &
      6.2907743e-05,8.0573298e-05,1.1916107e-04,2.2085457e-04,4.8168287e-04,   &
      1.1255995e-03,2.6169100e-03 ]
scav_coeff_mass(1:ncol, 3, 1) = [                                              &
      5.8865063e-05,3.3511807e-05,1.9608103e-05,1.1772057e-05,7.2633920e-06,   &
      4.6337085e-06,3.0897236e-06,2.1814105e-06,1.6443404e-06,1.3198936e-06,   &
      1.1191012e-06,1.0090169e-06,9.9761086e-07,1.1123232e-06,1.4030193e-06,   &
      3.2742033e-06,3.6751432e-05,8.1008027e-05,8.6383969e-05,8.2415738e-05,   &
      9.2956540e-05,1.4082319e-04 ]
scav_coeff_mass(1:ncol, 3, 2) = [                                              &
      4.6718202e-05,2.6817213e-05,1.5820184e-05,9.5863739e-06,5.9863461e-06,   &
      3.8832942e-06,2.6480621e-06,1.9198921e-06,1.4858752e-06,1.2221200e-06,   &
      1.0678705e-06,1.0103460e-06,1.0674572e-06,1.3287404e-06,3.7143753e-06,   &
      2.1296290e-05,5.9238788e-05,8.2347562e-05,8.4531440e-05,8.9650033e-05,   &
      1.2203916e-04,2.2020836e-04 ]
scav_coeff_mass(1:ncol, 3, 3) = [                                              &
      3.4359274e-05,1.9948103e-05,1.1909119e-05,7.3210088e-06,4.6607807e-06,   &
      3.1040085e-06,2.1879692e-06,1.6448437e-06,1.3192211e-06,1.1288584e-06,   &
      1.0443371e-06,1.0841229e-06,1.5427714e-06,4.8564869e-06,1.8809883e-05,   &
      4.6806845e-05,7.2653944e-05,8.3604141e-05,9.1491523e-05,1.2118555e-04,   &
      2.1198388e-04,4.5059897e-04 ]
scav_coeff_mass(1:ncol, 3, 4) = [                                              &
      2.4224351e-05,1.4261328e-05,8.6497348e-06,5.4261269e-06,3.5501215e-06,   &
      2.4494350e-06,1.7990594e-06,1.4125482e-06,1.1889538e-06,1.0895776e-06,   &
      1.1806136e-06,2.0650154e-06,6.3110116e-06,1.8969859e-05,4.1486846e-05,   &
      6.5095674e-05,8.0719393e-05,9.4265876e-05,1.2833979e-04,2.2605569e-04,   &
      4.8204648e-04,1.1262084e-03 ]
scav_coeff_mass(1:ncol, 3, 5) = [                                              &
      1.6722850e-05,1.0013141e-05,6.2002053e-06,3.9971731e-06,2.7101916e-06,   &
      1.9524972e-06,1.5053662e-06,1.2514240e-06,1.1597288e-06,1.4001977e-06,   &
      2.8403835e-06,7.8777049e-06,1.9866494e-05,3.9209357e-05,6.0530731e-05,   &
      7.8447370e-05,9.7994258e-05,1.4113923e-04,2.5665878e-04,5.5480118e-04,   &
      1.2923894e-03,3.0018139e-03 ]
scav_coeff_mass(1:ncol, 4, 1) = [                                              &
      6.2174955e-05,3.5477932e-05,2.0782433e-05,1.2486514e-05,7.7171032e-06,   &
      4.9466042e-06,3.3333879e-06,2.3954760e-06,1.8445944e-06,1.5035616e-06,   &
      1.2734719e-06,1.1266545e-06,1.0888987e-06,1.2026400e-06,1.5243501e-06,   &
      4.2706047e-06,5.0354259e-05,1.0339176e-04,1.0771350e-04,1.0135387e-04,   &
      1.1177062e-04,1.6504866e-04 ]
scav_coeff_mass(1:ncol, 4, 2) = [                                              &
      4.9373664e-05,2.8398545e-05,1.6771585e-05,1.0174995e-05,6.3726958e-06,   &
      4.1635909e-06,2.8780101e-06,2.1271299e-06,1.6769869e-06,1.3899341e-06,   &
      1.2032660e-06,1.1154200e-06,1.1602290e-06,1.4659631e-06,4.8225208e-06,   &
      2.8804503e-05,7.7081977e-05,1.0383739e-04,1.0465185e-04,1.0878505e-04,   &
      1.4427413e-04,2.5509006e-04 ]
scav_coeff_mass(1:ncol, 4, 3) = [                                              &
      3.6336622e-05,2.1133059e-05,1.2633391e-05,7.7841450e-06,4.9823298e-06,   &
      3.3539922e-06,2.4036357e-06,1.8401584e-06,1.4927363e-06,1.2739262e-06,   &
      1.1607486e-06,1.1905484e-06,1.7889954e-06,6.3352703e-06,2.5004438e-05,   &
      6.0922175e-05,9.2280461e-05,1.0402064e-04,1.1142176e-04,1.4377745e-04,   &
      2.4610724e-04,5.1740537e-04 ]
scav_coeff_mass(1:ncol, 4, 4) = [                                              &
      2.5636893e-05,1.5119589e-05,9.1903592e-06,5.7909781e-06,3.8226000e-06,   &
      2.6752631e-06,1.9983790e-06,1.5889699e-06,1.3393063e-06,1.2160884e-06,   &
      1.3258378e-06,2.5119789e-06,8.2310650e-06,2.4930036e-05,5.3752433e-05,   &
      8.2771822e-05,1.0063685e-04,1.1487383e-04,1.5227678e-04,2.6246917e-04,   &
      5.5356416e-04,1.2880232e-03 ]
scav_coeff_mass(1:ncol, 4, 5) = [                                              &
      1.7714692e-05,1.0632124e-05,6.6098416e-06,4.2938386e-06,2.9476229e-06,   &
      2.1564451e-06,1.6840448e-06,1.4060173e-06,1.3029455e-06,1.6240283e-06,   &
      3.5527390e-06,1.0247752e-05,2.5899479e-05,5.0546911e-05,7.6816918e-05,   &
      9.7696509e-05,1.1918002e-04,1.6709231e-04,2.9762442e-04,6.3681881e-04,   &
      1.4779011e-03,3.4286285e-03 ]
scav_coeff_mass(1:ncol, 5, 1) = [                                              &
      6.6022410e-05,3.7763585e-05,2.2151225e-05,1.3325469e-05,8.2593716e-06,   &
      5.3340584e-06,3.6516762e-06,2.6915668e-06,2.1343235e-06,1.7764272e-06,   &
      1.5039112e-06,1.2973164e-06,1.2111032e-06,1.3116243e-06,1.6630317e-06,   &
      5.6779651e-06,6.7986717e-05,1.3118650e-04,1.3419538e-04,1.2472807e-04,   &
      1.3450674e-04,1.9334211e-04 ]
scav_coeff_mass(1:ncol, 5, 2) = [                                              &
      5.2461067e-05,3.0238894e-05,1.7883799e-05,1.0871149e-05,6.8414381e-06,   &
      4.5191418e-06,3.1866488e-06,2.4200056e-06,1.9567670e-06,1.6394102e-06,   &
      1.4024814e-06,1.2623565e-06,1.2782647e-06,1.6346185e-06,6.3102336e-06,   &
      3.8586728e-05,9.9546866e-05,1.3054136e-04,1.2954456e-04,1.3211714e-04,   &
      1.7057435e-04,2.9496663e-04 ]
scav_coeff_mass(1:ncol, 5, 3) = [                                              &
      3.8637482e-05,2.2515805e-05,1.3485837e-05,8.3401847e-06,5.3830705e-06,   &
      3.6823789e-06,2.7025122e-06,2.1219429e-06,1.7484391e-06,1.4873349e-06,   &
      1.3261173e-06,1.3310312e-06,2.1085397e-06,8.2791366e-06,3.3000309e-05,   &
      7.8774565e-05,1.1677011e-04,1.2929290e-04,1.3576924e-04,1.7061558e-04,   &
      2.8526600e-04,5.9226060e-04 ]
scav_coeff_mass(1:ncol, 5, 4) = [                                              &
      2.7284249e-05,1.6127239e-05,9.8355802e-06,6.2406433e-06,4.1750323e-06,   &
      2.9832494e-06,2.2821031e-06,1.8464167e-06,1.5593669e-06,1.3963239e-06,   &
      1.5207807e-06,3.0966640e-06,1.0724527e-05,3.2568514e-05,6.9258715e-05,   &
      1.0487119e-04,1.2530086e-04,1.4002222e-04,1.8068889e-04,3.0423832e-04,   &
      6.3367996e-04,1.4673464e-03 ]
scav_coeff_mass(1:ncol, 5, 5) = [                                              &
      1.8877709e-05,1.1368131e-05,7.1108642e-06,4.6730790e-06,3.2671104e-06,   &
      2.4431948e-06,1.9421881e-06,1.6306157e-06,1.5057215e-06,1.9237466e-06,   &
      4.4805973e-06,1.3301649e-05,3.3592505e-05,6.4856458e-05,9.7177514e-05,   &
      1.2151264e-04,1.4495849e-04,1.9777419e-04,3.4446996e-04,7.2855628e-04,   &
      1.6833721e-03,3.8996803e-03 ]
scav_coeff_mass(1:ncol, 6, 1) = [                                              &
      7.0439935e-05,4.0388284e-05,2.3727442e-05,1.4299750e-05,8.9024966e-06,   &
      5.8131340e-06,4.0692331e-06,3.1034519e-06,2.5553863e-06,2.1836199e-06,   &
      1.8507485e-06,1.5492847e-06,1.3791757e-06,1.4453949e-06,1.8218224e-06,   &
      7.6242244e-06,9.0523819e-05,1.6561720e-04,1.6704913e-04,1.5357086e-04,   &
      1.6200857e-04,2.2641857e-04 ]
scav_coeff_mass(1:ncol, 6, 2) = [                                              &
      5.6006685e-05,3.2354708e-05,1.9168931e-05,1.1686650e-05,7.4075464e-06,   &
      4.9710926e-06,3.6033820e-06,2.8363024e-06,2.3683610e-06,2.0127028e-06,   &
      1.6993240e-06,1.4726046e-06,1.4322516e-06,1.8433755e-06,8.2766856e-06,   &
      5.1170271e-05,1.2766158e-04,1.6366547e-04,1.6032605e-04,1.6058233e-04,   &
      2.0172137e-04,3.4054102e-04 ]
scav_coeff_mass(1:ncol, 6, 3) = [                                              &
      4.1282153e-05,2.4110392e-05,1.4478891e-05,9.0036195e-06,5.8826422e-06,   &
      4.1161813e-06,3.1195249e-06,2.5309030e-06,2.1277281e-06,1.8046937e-06,   &
      1.5656560e-06,1.5208025e-06,2.5223494e-06,1.0802177e-05,4.3214389e-05,   &
      1.0120838e-04,1.4723480e-04,1.6054509e-04,1.6552390e-04,2.0253736e-04,   &
      3.3020389e-04,6.7593268e-04 ]
scav_coeff_mass(1:ncol, 6, 4) = [                                              &
      2.9182718e-05,1.7297846e-05,1.0600186e-05,6.7942896e-06,4.6331920e-06,   &
      3.4064296e-06,2.6888317e-06,2.2248140e-06,1.8848193e-06,1.6575214e-06,   &
      1.7861090e-06,3.8566024e-06,1.3928607e-05,4.2272406e-05,8.8743474e-05,   &
      1.3240423e-04,1.5579902e-04,1.7071865e-04,2.1445291e-04,3.5214916e-04,   &
      7.2320599e-04,1.6652607e-03 ]
scav_coeff_mass(1:ncol, 6, 5) = [                                              &
      2.0226878e-05,1.2236652e-05,7.7225089e-06,5.1600602e-06,3.7004028e-06,   &
      2.8496167e-06,2.3181542e-06,1.9603547e-06,1.7969520e-06,2.3272234e-06,   &
      5.6804210e-06,1.7201303e-05,4.3328050e-05,8.2816938e-05,1.2254287e-04,   &
      1.5093285e-04,1.7633129e-04,2.3408617e-04,3.9802276e-04,8.3088771e-04,   &
      1.9099949e-03,4.4170477e-03 ]
scav_coeff_mass(1:ncol, 7, 1) = [                                              &
      7.5462290e-05,4.3373111e-05,2.5525473e-05,1.5422099e-05,9.6620052e-06,   &
      6.4065020e-06,4.6198849e-06,3.6785871e-06,3.1678175e-06,2.7911550e-06,   &
      2.3741637e-06,1.9255898e-06,1.6161996e-06,1.6132992e-06,2.0044454e-06,   &
      1.0247291e-05,1.1894579e-04,2.0814418e-04,2.0774270e-04,1.8912905e-04,   &
      1.9528823e-04,2.6511492e-04 ]
scav_coeff_mass(1:ncol, 7, 2) = [                                              &
      6.0038726e-05,3.4763963e-05,2.0640814e-05,1.2635970e-05,8.0905340e-06,   &
      5.5481249e-06,4.1691664e-06,3.4297652e-06,2.9743206e-06,2.5720687e-06,   &
      2.1446289e-06,1.7789285e-06,1.6385921e-06,2.1034137e-06,1.0830137e-05,   &
      6.7144536e-05,1.6262634e-04,2.0465266e-04,1.9834164e-04,1.9530671e-04,   &
      2.3864055e-04,3.9260744e-04 ]
scav_coeff_mass(1:ncol, 7, 3) = [                                              &
      4.4292677e-05,2.5932650e-05,1.5627507e-05,9.7930809e-06,6.5074794e-06,   &
      4.6928411e-06,3.7041413e-06,3.1257684e-06,2.6913716e-06,2.2791354e-06,   &
      1.9175302e-06,1.7827837e-06,3.0561577e-06,1.4031368e-05,5.6117708e-05,   &
      1.2920370e-04,1.8499432e-04,1.9912670e-04,2.0187539e-04,2.4053675e-04,   &
      3.8176829e-04,7.6921620e-04 ]
scav_coeff_mass(1:ncol, 7, 4) = [                                              &
      3.1350519e-05,1.8647538e-05,1.1503017e-05,7.4775889e-06,5.2327082e-06,   &
      3.9914931e-06,3.2741103e-06,2.7825510e-06,2.3686302e-06,2.0405424e-06,   &
      2.1516604e-06,4.8369580e-06,1.7998390e-05,5.4483254e-05,1.1306403e-04,   &
      1.6656545e-04,1.9343153e-04,2.0817060e-04,2.5460779e-04,4.0709608e-04,   &
      8.2298034e-04,1.8827134e-03 ]
scav_coeff_mass(1:ncol, 7, 5) = [                                              &
      2.1779887e-05,1.3257300e-05,8.4704528e-06,5.7895173e-06,4.2923221e-06,   &
      3.4286882e-06,2.8679014e-06,2.4471583e-06,2.2193904e-06,2.8728252e-06,   &
      7.2195375e-06,2.2131767e-05,5.5545634e-05,1.0521923e-04,1.5401117e-04,   &
      1.8719289e-04,2.1449446e-04,2.7709150e-04,4.5921707e-04,9.4470199e-04,   &
      2.1587925e-03,4.9822431e-03 ]
scav_coeff_mass(1:ncol, 8, 1) = [                                              &
      8.1124862e-05,4.6739721e-05,2.7560506e-05,1.6706750e-05,1.0556382e-05,   &
      7.1424042e-06,5.3470405e-06,4.4792670e-06,4.0520341e-06,3.6892288e-06,   &
      3.1580461e-06,2.4874210e-06,1.9555619e-06,1.8288691e-06,2.2158268e-06,   &
      1.3679779e-05,1.5431208e-04,2.6047316e-04,2.5801885e-04,2.3289077e-04,   &
      2.3554581e-04,3.1039640e-04 ]
scav_coeff_mass(1:ncol, 8, 2) = [                                              &
      6.4586011e-05,3.7485349e-05,2.2314482e-05,1.3735902e-05,8.9143148e-06,   &
      6.2866636e-06,4.9373242e-06,4.2718107e-06,3.8593594e-06,3.4034053e-06,   &
      2.8098296e-06,2.2281511e-06,1.9209552e-06,2.4285659e-06,1.4080018e-05,   &
      8.7138758e-05,2.0580781e-04,2.5520029e-04,2.4519144e-04,2.3762899e-04,   &
      2.8241322e-04,4.5204416e-04 ]
scav_coeff_mass(1:ncol, 8, 3) = [                                              &
      4.7691832e-05,2.7999546e-05,1.6948751e-05,1.0731157e-05,7.2909505e-06,   &
      5.4609028e-06,4.5218383e-06,3.9857103e-06,3.5227172e-06,2.9847551e-06,   &
      2.4357788e-06,2.1494364e-06,3.7402784e-06,1.8100176e-05,7.2220838e-05,   &
      1.6386726e-04,2.3158551e-04,2.4663419e-04,2.4623546e-04,2.8577825e-04,   &
      4.4090589e-04,8.7288901e-04 ]
scav_coeff_mass(1:ncol, 8, 4) = [                                              &
      3.3807059e-05,2.0194561e-05,1.2566763e-05,8.3228389e-06,6.0196915e-06,   &
      4.8001355e-06,4.1126405e-06,3.5994301e-06,3.0842909e-06,2.6027591e-06,   &
      2.6582676e-06,6.0899368e-06,2.3101242e-05,6.9686434e-05,1.4318976e-04,   &
      2.0873848e-04,2.3972812e-04,2.5380571e-04,3.0236797e-04,4.7007658e-04,   &
      9.3382004e-04,2.1203718e-03 ]
scav_coeff_mass(1:ncol, 8, 5) = [                                              &
      2.3556632e-05,1.4453598e-05,9.3869416e-06,6.6063761e-06,5.1020622e-06,   &
      4.2515719e-06,3.6677387e-06,3.1627799e-06,2.8323825e-06,3.6110112e-06,   &
      9.1752366e-06,2.8296294e-05,7.0733733e-05,1.3296059e-04,1.9285204e-04,   &
      2.3174067e-04,2.6085556e-04,3.2802632e-04,5.2908475e-04,1.0708457e-03,   &
      2.4304500e-03,5.5957880e-03 ]
scav_coeff_mass(1:ncol, 9, 1) = [                                              &
      8.7461745e-05,5.0509088e-05,2.9847508e-05,1.8168284e-05,1.1605462e-05,   &
      8.0522544e-06,6.3002945e-06,5.5782757e-06,5.3039705e-06,4.9876748e-06,   &
      4.3066891e-06,3.3124368e-06,2.4404438e-06,2.1097768e-06,2.4620617e-06,   &
      1.8028442e-05,1.9771836e-04,3.2455599e-04,3.1991924e-04,2.8661284e-04,   &
      2.8418951e-04,3.6336325e-04 ]
scav_coeff_mass(1:ncol, 9, 2) = [                                              &
      6.9676370e-05,4.0537116e-05,2.4205063e-05,1.5004147e-05,9.9051666e-06,   &
      7.2279543e-06,5.9696656e-06,5.4469510e-06,5.1257107e-06,4.6123910e-06,   &
      3.7844900e-06,2.8800428e-06,2.3099690e-06,2.8346282e-06,1.8125757e-05,   &
      1.1178826e-04,2.5872087e-04,3.1727190e-04,3.0275522e-04,2.8912364e-04,   &
      3.3428923e-04,5.1980553e-04 ]
scav_coeff_mass(1:ncol, 9, 3) = [                                              &
      5.1501781e-05,3.0328019e-05,1.8460446e-05,1.1842503e-05,8.2706841e-06,   &
      6.4764502e-06,5.6498169e-06,5.2058362e-06,4.7236944e-06,4.0137413e-06,   &
      3.1887309e-06,2.6620306e-06,4.6078234e-06,2.3138887e-05,9.2050973e-05,   &
      2.0641237e-04,2.8876364e-04,3.0492965e-04,3.0026009e-04,3.3961106e-04,   &
      5.0865769e-04,9.8766574e-04 ]
scav_coeff_mass(1:ncol, 9, 4) = [                                              &
      3.6571637e-05,2.1957870e-05,1.3816114e-05,9.3664117e-06,7.0473688e-06,   &
      5.9050092e-06,5.2939667e-06,4.7727084e-06,4.1227798e-06,3.4161418e-06,   &
      3.3564581e-06,7.6718989e-06,2.9407596e-05,8.8393537e-05,1.8018522e-04,   &
      2.6049314e-04,2.9646081e-04,3.0929057e-04,3.5913648e-04,5.4218507e-04,   &
      1.0564707e-03,2.3784686e-03 ]
scav_coeff_mass(1:ncol, 9, 5) = [                                              &
      2.5577743e-05,1.5851111e-05,1.0508292e-05,7.6625165e-06,6.1993265e-06,   &
      5.4034901e-06,4.8104695e-06,4.1957082e-06,3.7096793e-06,4.6023956e-06,   &
      1.1630948e-05,3.5907415e-05,8.9415276e-05,1.6703088e-04,2.4050260e-04,   &
      2.8624490e-04,3.1704878e-04,3.8830970e-04,6.0874462e-04,1.2100621e-03,   &
      2.7251344e-03,6.2567626e-03 ]
scav_coeff_mass(1:ncol,10, 1) = [                                              &
      9.4503917e-05,5.4700162e-05,3.2399813e-05,1.9819515e-05,1.2826775e-05,   &
      9.1644428e-06,7.5256819e-06,7.0450071e-06,7.0175175e-06,6.7967953e-06,   &
      5.9275452e-06,4.4826864e-06,3.1175764e-06,2.4755370e-06,2.7497674e-06,   &
      2.3351376e-05,2.5023929e-04,4.0258401e-04,3.9580811e-04,3.5235096e-04,   &
      3.4285950e-04,4.2526209e-04 ]
scav_coeff_mass(1:ncol,10, 2) = [                                              &
      7.5335028e-05,4.3935656e-05,2.6325950e-05,1.6456342e-05,1.1086690e-05,   &
      8.4099702e-06,7.3245911e-06,7.0370831e-06,6.8749104e-06,6.3065285e-06,   &
      5.1618646e-06,3.7982397e-06,2.8388964e-06,3.3371116e-06,2.3043577e-05,   &
      1.4169263e-04,3.2300052e-04,3.9310504e-04,3.7321833e-04,3.5162717e-04,   &
      3.9570373e-04,5.9692085e-04 ]
scav_coeff_mass(1:ncol,10, 3) = [                                              &
      5.5742505e-05,3.2933162e-05,2.0178400e-05,1.3149277e-05,9.4812690e-06,   &
      7.7923218e-06,7.1625847e-06,6.8800873e-06,6.3972453e-06,5.4611997e-06,   &
      4.2483097e-06,3.3645392e-06,5.6902336e-06,2.9262298e-05,1.1612305e-04,   &
      2.5813040e-04,3.5849826e-04,3.7615875e-04,3.6587389e-04,4.0358739e-04,   &
      5.8616078e-04,1.1141649e-03 ]
scav_coeff_mass(1:ncol,10, 4) = [                                              &
      3.9661530e-05,2.3954366e-05,1.5273334e-05,1.0641805e-05,8.3659125e-06,   &
      7.3761707e-06,6.9063228e-06,6.4002083e-06,5.5773958e-06,4.5557802e-06,   &
      4.2988880e-06,9.6367534e-06,3.7078435e-05,1.1111950e-04,2.2518623e-04,   &
      3.2357698e-04,3.6565486e-04,3.7655123e-04,4.2652202e-04,6.2461384e-04,   &
      1.1915708e-03,2.6566812e-03 ]
scav_coeff_mass(1:ncol,10, 5) = [                                              &
      2.7861705e-05,1.7473000e-05,1.1868066e-05,9.0069450e-06,7.6513762e-06,   &
      6.9683530e-06,6.3892151e-06,5.6363022e-06,4.9275818e-06,5.9089737e-06,   &
      1.4667872e-05,4.5174097e-05,1.1212815e-04,2.0849267e-04,2.9855833e-04,   &
      3.5260152e-04,3.8495088e-04,4.5955770e-04,6.9939975e-04,1.3629470e-03,   &
      3.0423535e-03,6.9624439e-03 ]
scav_coeff_mass(1:ncol,11, 1) = [                                              &
      1.0227785e-04,5.9328663e-05,3.5227476e-05,2.1668507e-05,1.4229920e-05,   &
      1.0494287e-05,9.0492659e-06,8.9212622e-06,9.2527908e-06,9.1913699e-06,   &
      8.0974450e-06,6.0598152e-06,4.0237507e-06,2.9422749e-06,3.0846356e-06,   &
      2.9637124e-05,3.1285925e-04,4.9697419e-04,4.8839585e-04,4.3249350e-04,   &
      4.1345738e-04,4.9750583e-04 ]
scav_coeff_mass(1:ncol,11, 2) = [                                              &
      8.1583280e-05,4.7694013e-05,2.8686394e-05,1.8101602e-05,1.2471760e-05,   &
      9.8539437e-06,9.0365904e-06,9.0934876e-06,9.1741654e-06,8.5607184e-06,   &
      7.0100663e-06,5.0313034e-06,3.5342087e-06,3.9471213e-06,2.8872922e-05,   &
      1.7736938e-04,4.0036510e-04,4.8521558e-04,4.5909934e-04,4.2726902e-04,   &
      4.6830140e-04,6.8450663e-04 ]
scav_coeff_mass(1:ncol,11, 3) = [                                              &
      6.0430254e-05,3.5825931e-05,2.2112343e-05,1.4663894e-05,1.0942159e-05,   &
      9.4396087e-06,9.1064047e-06,9.0699142e-06,8.6140289e-06,7.3953448e-06,   &
      5.6682263e-06,4.2908402e-06,7.0097131e-06,3.6555913e-05,1.4490748e-04,   &
      3.2035582e-04,4.4296420e-04,4.6276854e-04,4.4529828e-04,4.7948772e-04,   &
      6.7466288e-04,1.2529034e-03 ]
scav_coeff_mass(1:ncol,11, 4) = [                                              &
      4.3089621e-05,2.6194871e-05,1.6951267e-05,1.2168127e-05,1.0004993e-05,   &
      9.2570923e-06,9.0071510e-06,8.5485126e-06,7.5142733e-06,6.0768631e-06,   &
      5.5251192e-06,1.2025173e-05,4.6250276e-05,1.3835645e-04,2.7937103e-04,   &
      3.9990233e-04,4.4959849e-04,4.5779641e-04,5.0636283e-04,7.1866754e-04,   &
      1.3396452e-03,2.9540811e-03 ]
scav_coeff_mass(1:ncol,11, 5) = [                                              &
      3.0420691e-05,1.9332986e-05,1.3485736e-05,1.0668902e-05,9.5000863e-06,   &
      9.0007311e-06,8.4670622e-06,7.5480957e-06,6.5415136e-06,7.5771451e-06,   &
      1.8351570e-05,5.6285583e-05,1.3940226e-04,2.5845715e-04,3.6876023e-04,   &
      4.3293874e-04,4.6669959e-04,5.4360431e-04,8.0234988e-04,1.5299382e-03,   &
      3.3808954e-03,7.7081360e-03 ]
scav_coeff_mass(1:ncol,12, 1) = [                                              &
      1.1080475e-04,6.4406221e-05,3.8335775e-05,2.3715529e-05,1.5810478e-05,   &
      1.2032807e-05,1.0858323e-05,1.1192844e-05,1.1998009e-05,1.2166364e-05,   &
      1.0819968e-05,8.0529470e-06,5.1678713e-06,3.5155814e-06,3.4694541e-06,   &
      3.6789490e-05,3.8638738e-04,6.1034052e-04,6.0075583e-04,5.2979458e-04,   &
      4.9817835e-04,5.8170076e-04 ]
scav_coeff_mass(1:ncol,12, 2) = [                                              &
      8.8437639e-05,5.1820652e-05,3.1289114e-05,1.9937765e-05,1.4053614e-05,   &
      1.1549048e-05,1.1092373e-05,1.1603486e-05,1.2015392e-05,1.1374358e-05,   &
      9.3352838e-06,6.5879787e-06,4.4029831e-06,4.6661582e-06,3.5604087e-05,   &
      2.1920397e-04,4.9256540e-04,5.9639140e-04,5.6327388e-04,5.1850350e-04,   &
      5.5396584e-04,7.8379186e-04 ]
scav_coeff_mass(1:ncol,12, 3) = [                                              &
      6.5576317e-05,3.9010936e-05,2.4261634e-05,1.6381040e-05,1.2643950e-05,   &
      1.1406169e-05,1.1468976e-05,1.1766317e-05,1.1371214e-05,9.8198254e-06,   &
      7.4558361e-06,5.4479161e-06,8.5699421e-06,4.5061811e-05,1.7879487e-04,   &
      3.9442179e-04,5.4452192e-04,5.6751933e-04,5.4107906e-04,5.6935040e-04,   &
      7.7554919e-04,1.4043225e-03 ]
scav_coeff_mass(1:ncol,12, 4) = [                                              &
      4.6862145e-05,2.8679889e-05,1.8845614e-05,1.3937024e-05,1.1953511e-05,   &
      1.1536219e-05,1.1587448e-05,1.1213754e-05,9.9353308e-06,7.9852216e-06,   &
      7.0419647e-06,1.4851465e-05,5.7018804e-05,1.7054424e-04,3.4392386e-04,   &
      4.9152411e-04,5.5084626e-04,5.5553920e-04,6.0075466e-04,8.2579070e-04,   &
      1.5011326e-03,3.2691664e-03 ]
scav_coeff_mass(1:ncol,12, 5) = [                                              &
      3.3256178e-05,2.1427563e-05,1.5353780e-05,1.2638205e-05,1.1734722e-05,   &
      1.1491981e-05,1.1039721e-05,9.9318832e-06,8.5562546e-06,9.6160927e-06,   &
      2.2715685e-05,6.9393061e-05,1.7173298e-04,3.1805351e-04,4.5297385e-04,   &
      5.2961593e-04,5.6471015e-04,6.4252685e-04,9.1901892e-04,1.7113447e-03,   &
      3.7388646e-03,8.4872283e-03 ]
scav_coeff_mass(1:ncol,13, 1) = [                                              &
      1.2010022e-04,6.9939963e-05,4.1724398e-05,2.5951233e-05,1.7546167e-05,   &
      1.3739324e-05,1.2888450e-05,1.3769371e-05,1.5141449e-05,1.5603226e-05,   &
      1.3991911e-05,1.0392583e-05,6.5159097e-06,4.1843481e-06,3.9023928e-06,   &
      4.4620635e-05,4.7134503e-04,7.4543041e-04,7.3631662e-04,6.4739236e-04,   &
      5.9953507e-04,6.7967115e-04 ]
scav_coeff_mass(1:ncol,13, 2) = [                                              &
      9.5909453e-05,5.6318797e-05,3.4128908e-05,2.1948447e-05,1.5800047e-05,   &
      1.3442019e-05,1.3414107e-05,1.4466222e-05,1.5284491e-05,1.4638140e-05,   &
      1.2052446e-05,8.4168828e-06,5.4222745e-06,5.4818249e-06,4.3168068e-05,   &
      2.6739198e-04,6.0130423e-04,7.2965705e-04,6.8897893e-04,6.2813027e-04,   &
      6.5484390e-04,8.9614674e-04 ]
scav_coeff_mass(1:ncol,13, 3) = [                                              &
      7.1186390e-05,4.2485158e-05,2.6612595e-05,1.8272477e-05,1.4539098e-05,   &
      1.3621591e-05,1.4157525e-05,1.4861567e-05,1.4560880e-05,1.2644028e-05,   &
      9.5493681e-06,6.8019436e-06,1.0348295e-05,5.4765114e-05,2.1805430e-04,   &
      4.8159460e-04,6.6567016e-04,6.9347433e-04,6.5609986e-04,6.7549454e-04,   &
      8.9037118e-04,1.5688311e-03 ]
scav_coeff_mass(1:ncol,13, 4) = [                                              &
      5.0977372e-05,3.1396961e-05,2.0929887e-05,1.5903813e-05,1.4145397e-05,   &
      1.4126432e-05,1.4545286e-05,1.4291699e-05,1.2749285e-05,1.0213734e-05,   &
      8.8073902e-06,1.8092276e-05,6.9422136e-05,2.0803447e-04,4.1998165e-04,   &
      6.0059365e-04,6.7219993e-04,6.7260307e-04,7.1207086e-04,9.4759711e-04,   &
      1.6764319e-03,3.5999548e-03 ]
scav_coeff_mass(1:ncol,13, 5) = [                                              &
      3.6356184e-05,2.3730881e-05,1.7428893e-05,1.4851430e-05,1.4272260e-05,   &
      1.4345106e-05,1.4007120e-05,1.2697763e-05,1.0902256e-05,1.1980110e-05,   &
      2.7747553e-05,8.4590293e-05,2.0954782e-04,3.8838332e-04,5.5314737e-04,   &
      6.4520133e-04,6.8167579e-04,7.5866311e-04,1.0509852e-03,1.9073984e-03,   &
      4.1137885e-03,9.2914325e-03 ]
scav_coeff_mass(1:ncol,14, 1) = [                                              &
      1.3017396e-04,7.5932455e-05,4.5387641e-05,2.8357187e-05,1.9397740e-05,   &
      1.5542620e-05,1.5024650e-05,1.6484636e-05,1.8470214e-05,1.9266448e-05,   &
      1.7398017e-05,1.2925152e-05,7.9870425e-06,4.9189243e-06,4.3764819e-06,   &
      5.2853756e-05,5.6780803e-04,9.0499195e-04,8.9879648e-04,7.8878561e-04,   &
      7.2035162e-04,7.9346229e-04 ]
scav_coeff_mass(1:ncol,14, 2) = [                                              &
      1.0400472e-04,6.1186517e-05,3.7193032e-05,2.4103750e-05,1.7654429e-05,   &
      1.5438245e-05,1.5860078e-05,1.7492150e-05,1.8758964e-05,1.8129729e-05,   &
      1.4979965e-05,1.0401947e-05,6.5362989e-06,6.3666734e-06,5.1429121e-05,   &
      3.2186617e-04,7.2810009e-04,8.8817611e-04,8.3976775e-04,7.5927888e-04,   &
      7.7334428e-04,1.0230940e-03 ]
scav_coeff_mass(1:ncol,14, 3) = [                                              &
      7.7260577e-05,4.6238237e-05,2.9139122e-05,2.0287939e-05,1.6542923e-05,   &
      1.5957880e-05,1.6998539e-05,1.8147357e-05,1.7966393e-05,1.5678316e-05,   &
      1.1813329e-05,8.2749557e-06,1.2293126e-05,6.5581381e-05,2.6277878e-04,   &
      5.8296778e-04,8.0894174e-04,8.4393679e-04,7.9355541e-04,8.0051376e-04,   &
      1.0208557e-03,1.7468410e-03 ]
scav_coeff_mass(1:ncol,14, 4) = [                                              &
      5.5425851e-05,3.4321212e-05,2.3156223e-05,1.7988380e-05,1.6460232e-05,   &
      1.6864828e-05,1.7684203e-05,1.7574574e-05,1.5767363e-05,1.2617928e-05,   &
      1.0726767e-05,2.1681614e-05,8.3424052e-05,2.5104161e-04,5.0854700e-04,   &
      7.2926267e-04,8.1663807e-04,8.1208537e-04,8.4295096e-04,1.0858781e-03,   &
      1.8659403e-03,3.9440892e-03 ]
scav_coeff_mass(1:ncol,14, 5) = [                                              &
      3.9695768e-05,2.6195476e-05,1.9632774e-05,1.7192427e-05,1.6957125e-05,   &
      1.7373266e-05,1.7170528e-05,1.5661225e-05,1.3431421e-05,1.4564229e-05,   &
      3.3380726e-05,1.0189344e-04,2.5315989e-04,4.7044478e-04,6.7122527e-04,   &
      7.8240138e-04,8.2052540e-04,8.9459557e-04,1.1999898e-03,2.1182990e-03,   &
      4.5027397e-03,1.0111084e-02 ]
scav_coeff_mass(1:ncol,15, 1) = [                                              &
      1.4102872e-04,8.2381580e-05,4.9315491e-05,3.0908875e-05,2.1315157e-05,   &
      1.7351966e-05,1.7118927e-05,1.9121738e-05,2.1702014e-05,2.2838228e-05,   &
      2.0742141e-05,1.5434841e-05,9.4649535e-06,5.6752830e-06,4.8807576e-06,   &
      6.1135890e-05,6.7519159e-04,1.0915255e-03,1.0920305e-03,9.5772555e-04,   &
      8.6369216e-04,9.2529163e-04 ]
scav_coeff_mass(1:ncol,15, 2) = [                                              &
      1.1272353e-04,6.6417260e-05,4.0463339e-05,2.6365015e-05,1.9544514e-05,   &
      1.7416292e-05,1.8246201e-05,2.0431469e-05,2.2140862e-05,2.1546222e-05,   &
      1.7865830e-05,1.2378838e-05,7.6642543e-06,7.2813044e-06,6.0180477e-05,   &
      3.8220265e-04,8.7406430e-04,1.0750446e-03,1.0193695e-03,9.1531856e-04,   &
      9.1207752e-04,1.1662723e-03 ]
scav_coeff_mass(1:ncol,15, 3) = [                                              &
      8.3793658e-05,5.0254278e-05,3.1806860e-05,2.2362968e-05,1.8546628e-05,   &
      1.8248916e-05,1.9763847e-05,2.1345736e-05,2.1294169e-05,1.8661468e-05,   &
      1.4057835e-05,9.7556918e-06,1.4327862e-05,7.7344954e-05,3.1281158e-04,   &
      6.9929280e-04,9.7670233e-04,1.0222897e-03,9.5684356e-04,9.4720578e-04,   &
      1.1688625e-03,1.9387575e-03 ]
scav_coeff_mass(1:ncol,15, 4) = [                                              &
      6.0192095e-05,3.7419340e-05,2.5462840e-05,2.0087495e-05,1.8741543e-05,   &
      1.9537216e-05,2.0742429e-05,2.0781608e-05,1.8730715e-05,1.4996658e-05,   &
      1.2665522e-05,2.5514295e-05,9.8896628e-05,2.9957414e-04,6.1034770e-04,   &
      8.7950301e-04,9.8715241e-04,9.7723734e-04,9.9622103e-04,1.2425555e-03,   &
      2.0700430e-03,4.2988945e-03 ]
scav_coeff_mass(1:ncol,15, 5) = [                                              &
      4.3241039e-05,2.8759670e-05,2.1864152e-05,1.9509950e-05,1.9584579e-05,   &
      2.0327576e-05,2.0261806e-05,1.8569998e-05,1.5938178e-05,1.7217362e-05,   &
      3.9496506e-05,1.2121890e-04,3.0270041e-04,5.6500796e-04,8.0898640e-04,   &
      9.4390267e-04,9.8429839e-04,1.0530652e-03,1.3678871e-03,2.3442070e-03,   &
      4.9024070e-03,1.0935374e-02 ]
scav_coeff_mass(1:ncol,16, 1) = [                                              &
      1.5265715e-04,8.9279445e-05,5.3494687e-05,3.3579958e-05,2.3247062e-05,   &
      1.9074724e-05,1.9019262e-05,2.1455763e-05,2.4540929e-05,2.5981443e-05,   &
      2.3706130e-05,1.7686673e-05,1.0821211e-05,6.4040610e-06,5.4027420e-06,   &
      6.9061124e-05,7.9198513e-04,1.3068683e-03,1.3196213e-03,1.1579575e-03,   &
      1.0326699e-03,1.0774028e-03 ]
scav_coeff_mass(1:ncol,16, 2) = [                                              &
      1.2205799e-04,7.1999927e-05,4.3919114e-05,2.8691960e-05,2.1396345e-05,   &
      1.9251611e-05,2.0382183e-05,2.3023391e-05,2.5115750e-05,2.4564267e-05,   &
      2.0437804e-05,1.4167845e-05,8.7166449e-06,8.1809548e-06,6.9144275e-05,   &
      4.4750869e-04,1.0395654e-03,1.2929166e-03,1.2313876e-03,1.0996335e-03,   &
      1.0736874e-03,1.3273085e-03 ]
scav_coeff_mass(1:ncol,16, 3) = [                                              &
      9.0774592e-05,5.4514054e-05,3.4579392e-05,2.4431235e-05,2.0438478e-05,   &
      2.0323041e-05,2.2215585e-05,2.4164073e-05,2.4231896e-05,2.1312630e-05,   &
      1.6076525e-05,1.1121661e-05,1.6360749e-05,8.9798706e-05,3.6765215e-04,   &
      8.3072783e-04,1.1708061e-03,1.2316764e-03,1.1493153e-03,1.1183847e-03,   &
      1.3362431e-03,2.1448788e-03 ]
scav_coeff_mass(1:ncol,16, 4) = [                                              &
      6.5256531e-05,4.0655505e-05,2.7785784e-05,2.2094866e-05,2.0827481e-05,   &
      2.1920304e-05,2.3444612e-05,2.3614678e-05,2.1361864e-05,1.7132119e-05,   &
      1.4474660e-05,2.9456203e-05,1.1560169e-04,3.5334194e-04,7.2562487e-04,   &
      1.0528028e-03,1.1864354e-03,1.1712012e-03,1.1746913e-03,1.4195304e-03,   &
      2.2890077e-03,4.6613273e-03 ]
scav_coeff_mass(1:ncol,16, 5) = [                                              &
      4.6955055e-05,3.1359233e-05,2.4018432e-05,2.1647297e-05,2.1941032e-05,   &
      2.2946276e-05,2.2996531e-05,2.1154117e-05,1.8199956e-05,1.9769087e-05,   &
      4.5932412e-05,1.4235732e-04,3.5802417e-04,6.7242304e-04,9.6777116e-04,   &
      1.1320760e-03,1.1758800e-03,1.2367601e-03,1.5564856e-03,2.5851350e-03,   &
      5.3090489e-03,1.1752415e-02 ]
scav_coeff_mass(1:ncol,17, 1) = [                                              &
      1.6503346e-04,9.6608170e-05,5.7907931e-05,3.6345463e-05,2.5149792e-05,   &
      2.0634162e-05,2.0599794e-05,2.3299098e-05,2.6737597e-05,2.8408777e-05,   &
      2.6015557e-05,1.9475443e-05,1.1942262e-05,7.0610998e-06,5.9312232e-06,   &
      7.6203696e-05,9.1547866e-04,1.5515636e-03,1.5843175e-03,1.3927115e-03,   &
      1.2300486e-03,1.2517504e-03 ]
scav_coeff_mass(1:ncol,17, 2) = [                                              &
      1.3198590e-04,7.7916436e-05,4.7538433e-05,3.1049054e-05,2.3148062e-05,   &
      2.0840982e-05,2.2109651e-05,2.5048990e-05,2.7416947e-05,2.6906627e-05,   &
      2.2459852e-05,1.5611409e-05,9.6141523e-06,9.0232889e-06,7.7977696e-05,   &
      5.1630837e-04,1.2237707e-03,1.5433919e-03,1.4787393e-03,1.3151704e-03,   &
      1.2604942e-03,1.5075312e-03 ]
scav_coeff_mass(1:ncol,17, 3) = [                                              &
      9.8182867e-05,5.8995383e-05,3.7423441e-05,2.6436577e-05,2.2125527e-05,   &
      2.2037317e-05,2.4154280e-05,2.6354723e-05,2.6512669e-05,2.3389377e-05,   &
      1.7689493e-05,1.2264473e-05,1.8296885e-05,1.0258693e-04,4.2634880e-04,   &
      9.7649713e-04,1.3920585e-03,1.4744385e-03,1.3737909e-03,1.3164922e-03,   &
      1.5245300e-03,2.3651410e-03 ]
scav_coeff_mass(1:ncol,17, 4) = [                                              &
      7.0595256e-05,4.3996019e-05,3.0070260e-05,2.3921647e-05,2.2583017e-05,   &
      2.3826770e-05,2.5557899e-05,2.5819447e-05,2.3421952e-05,1.8834929e-05,   &
      1.6020017e-05,3.3357219e-05,1.3317178e-04,4.1164219e-04,8.5383998e-04,   &
      1.2496985e-03,1.4163463e-03,1.3965176e-03,1.3807469e-03,1.6183531e-03,   &
      2.5227140e-03,5.0277459e-03 ]
scav_coeff_mass(1:ncol,17, 5) = [                                              &
      5.0802345e-05,3.3938559e-05,2.6007732e-05,2.3473396e-05,2.3847107e-05,   &
      2.5007999e-05,2.5132264e-05,2.3181506e-05,2.0022653e-05,2.2060326e-05,   &
      5.2492699e-05,1.6494300e-04,4.1858680e-04,7.9234861e-04,1.1480598e-03,   &
      1.3484781e-03,1.3975163e-03,1.4478964e-03,1.7672037e-03,2.8406656e-03,   &
      5.7182532e-03,1.2549013e-02 ]
scav_coeff_mass(1:ncol,18, 1) = [                                              &
      1.7809283e-04,1.0432826e-04,6.2528146e-05,3.9180876e-05,2.6991825e-05,   &
      2.1981027e-05,2.1782000e-05,2.4534397e-05,2.8134010e-05,2.9935111e-05,   &
      2.7490369e-05,2.0663967e-05,1.2750716e-05,7.6156230e-06,6.4580016e-06,   &
      8.2158532e-05,1.0415600e-03,1.8239835e-03,1.8870028e-03,1.6638010e-03,   &
      1.4575097e-03,1.4494069e-03 ]
scav_coeff_mass(1:ncol,18, 2) = [                                              &
      1.4245458e-04,8.4133038e-05,5.1294933e-05,3.3407479e-05,2.4758195e-05,   &
      2.2119248e-05,2.3329529e-05,2.6370169e-05,2.8873834e-05,2.8393021e-05,   &
      2.3775880e-05,1.6603543e-05,1.0302431e-05,9.7746819e-06,8.6287555e-05,   &
      5.8645576e-04,1.4240812e-03,1.8260933e-03,1.7627120e-03,1.5636257e-03,   &
      1.4738323e-03,1.7074241e-03 ]
scav_coeff_mass(1:ncol,18, 3) = [                                              &
      1.0597718e-04,6.3668168e-05,4.0309408e-05,2.8339685e-05,2.3548151e-05,   &
      2.3301897e-05,2.5454118e-05,2.7759672e-05,2.7963748e-05,2.4732441e-05,   &
      1.8776725e-05,1.3109722e-05,2.0048766e-05,1.1525472e-04,4.8739560e-04,   &
      1.1344722e-03,1.6394420e-03,1.7512058e-03,1.6317148e-03,1.5428874e-03,   &
      1.7343491e-03,2.5986092e-03 ]
scav_coeff_mass(1:ncol,18, 4) = [                                              &
      7.6173633e-05,4.7408949e-05,3.2276505e-05,2.5509902e-05,2.3922694e-05,   &
      2.5138196e-05,2.6933770e-05,2.7231686e-05,2.4754665e-05,1.9979094e-05,   &
      1.7205327e-05,3.7062502e-05,1.5109374e-04,4.7323486e-04,9.9330346e-04,   &
      1.4691059e-03,1.6770700e-03,1.6542876e-03,1.6156121e-03,1.8396035e-03,   &
      2.7701166e-03,5.3933893e-03 ]
scav_coeff_mass(1:ncol,18, 5) = [                                              &
      5.4747810e-05,3.6456097e-05,2.7773864e-05,2.4904672e-05,2.5189024e-05,   &
      2.6371427e-05,2.6512596e-05,2.4500495e-05,2.1275831e-05,2.3967519e-05,   &
      5.8956393e-05,1.8842082e-04,4.8329990e-04,9.2339659e-04,1.3488678e-03,   &
      1.5930743e-03,1.6500033e-03,1.6874760e-03,2.0004287e-03,3.1093862e-03,   &
      6.1243845e-03,1.3310020e-02 ]
scav_coeff_mass(1:ncol,19, 1) = [                                              &
      1.9168898e-04,1.1235384e-04,6.7304105e-05,4.2054330e-05,2.8750759e-05,   &
      2.3094934e-05,2.2540989e-05,2.5126532e-05,2.8681312e-05,3.0499569e-05,   &
      2.8067120e-05,2.1200355e-05,1.3214927e-05,8.0533770e-06,6.9774789e-06,   &
      8.6583777e-05,1.1646868e-03,2.1192455e-03,2.2252210e-03,1.9702127e-03,   &
      1.7144678e-03,1.6695920e-03 ]
scav_coeff_mass(1:ncol,19, 2) = [                                              &
      1.5334703e-04,9.0580621e-05,5.5146679e-05,3.5739753e-05,2.6204950e-05,   &
      2.3063287e-05,2.4011312e-05,2.6944617e-05,2.9431675e-05,2.8961988e-05,   &
      2.4329150e-05,1.7103055e-05,1.0758473e-05,1.0412955e-05,9.3654628e-05,   &
      6.5511498e-04,1.6355263e-03,2.1374149e-03,2.0815400e-03,1.8441574e-03,   &
      1.7129744e-03,1.9257231e-03 ]
scav_coeff_mass(1:ncol,19, 3) = [                                              &
      1.1407056e-04,6.8480203e-05,4.3204060e-05,3.0115962e-05,2.4682692e-05,   &
      2.4087777e-05,2.6076111e-05,2.8328551e-05,2.8527190e-05,2.5285278e-05,   &
      1.9292986e-05,1.3625821e-05,2.1542790e-05,1.2725796e-04,5.4866089e-04,   &
      1.3007199e-03,1.9091089e-03,2.0595752e-03,1.9218465e-03,1.7967050e-03,   &
      1.9644533e-03,2.8426160e-03 ]
scav_coeff_mass(1:ncol,19, 4) = [                                              &
      8.1929155e-05,5.0855074e-05,3.4376542e-05,2.6834398e-05,2.4817482e-05,   &
      2.5817110e-05,2.7524699e-05,2.7796659e-05,2.5305242e-05,2.0517446e-05,   &
      1.7982840e-05,4.0419783e-05,1.6870077e-04,5.3622996e-04,1.1407564e-03,   &
      1.7074572e-03,1.9659209e-03,1.9428990e-03,1.8781831e-03,2.0818789e-03,   &
      3.0283368e-03,5.7514477e-03 ]
scav_coeff_mass(1:ncol,19, 5) = [                                              &
      5.8746038e-05,3.8880166e-05,2.9289491e-05,2.5911305e-05,2.5929870e-05,   &
      2.6990919e-05,2.7085421e-05,2.5058062e-05,2.1908225e-05,2.5413902e-05,   &
      6.5081051e-05,2.1201730e-04,5.5038246e-04,1.0627157e-03,1.5669596e-03,   &
      1.8631447e-03,1.9314717e-03,1.9541224e-03,2.2544742e-03,3.3879365e-03,   &
      6.5195931e-03,1.4017044e-02 ]
scav_coeff_mass(1:ncol,20, 1) = [                                              &
      2.0552655e-04,1.2051263e-04,7.2136228e-05,4.4911228e-05,3.0402518e-05,   &
      2.3975475e-05,2.2896877e-05,2.5113453e-05,2.8430207e-05,3.0156386e-05,   &
      2.7791574e-05,2.1113051e-05,1.3345911e-05,8.3741706e-06,7.4837191e-06,   &
      8.9237650e-05,1.2781326e-03,2.4281236e-03,2.5913814e-03,2.3062741e-03,   &
      1.9964904e-03,1.9083615e-03 ]
scav_coeff_mass(1:ncol,20, 2) = [                                              &
      1.6442758e-04,9.7122332e-05,5.9016209e-05,3.8006511e-05,2.7476266e-05,   &
      2.3683194e-05,2.4184145e-05,2.6816457e-05,2.9142336e-05,2.8662710e-05,   &
      2.4156148e-05,1.7129508e-05,1.0987776e-05,1.0925914e-05,9.9666822e-05,   &
      7.1884859e-04,1.8502739e-03,2.4691078e-03,2.4286202e-03,2.1516940e-03,   &
      1.9736889e-03,2.1581908e-03 ]
scav_coeff_mass(1:ncol,20, 3) = [                                              &
      1.2228976e-04,7.3332540e-05,4.6054712e-05,3.1744245e-05,2.5532097e-05,   &
      2.4417831e-05,2.6058928e-05,2.8109467e-05,2.8251491e-05,2.5087365e-05,   &
      1.9263054e-05,1.3820907e-05,2.2721135e-05,1.3798822e-04,6.0738545e-04,   &
      1.4691289e-03,2.1932944e-03,2.3925075e-03,2.2385740e-03,2.0733403e-03,   &
      2.2104176e-03,3.0915805e-03 ]
scav_coeff_mass(1:ncol,20, 4) = [                                              &
      8.7742108e-05,5.4269397e-05,3.6341435e-05,2.7892516e-05,2.5285579e-05,   &
      2.5897873e-05,2.7375137e-05,2.7560793e-05,2.5113530e-05,2.0476512e-05,   &
      1.8350638e-05,4.3283073e-05,1.8518179e-04,5.9802566e-04,1.2909976e-03,   &
      1.9577753e-03,2.2759211e-03,2.2564148e-03,2.1634953e-03,2.3404343e-03,   &
      3.2914209e-03,6.0917602e-03 ]
scav_coeff_mass(1:ncol,20, 5) = [                                              &
      6.2720054e-05,4.1174646e-05,3.0547220e-05,2.6507692e-05,2.6100459e-05,   &
      2.6907617e-05,2.6894701e-05,2.4892848e-05,2.1942605e-05,2.6367369e-05,   &
      7.0604466e-05,2.3472897e-04,6.1724958e-04,1.2055970e-03,1.7959987e-03,   &
      2.1519940e-03,2.2358672e-03,2.2425669e-03,2.5241881e-03,3.6697062e-03,   &
      6.8924256e-03,1.4646616e-02 ]
scav_coeff_mass(1:ncol,21, 1) = [                                              &
      2.1908395e-04,1.2850022e-04,7.6848478e-05,4.7655422e-05,3.1906320e-05,   &
      2.4626989e-05,2.2896333e-05,2.4582640e-05,2.7502284e-05,2.9043693e-05,   &
      2.6790039e-05,2.0489711e-05,1.3185023e-05,8.5853992e-06,7.9658604e-06,   &
      9.0002160e-05,1.3745584e-03,2.7363928e-03,2.9711936e-03,2.6599098e-03,   &
      2.2937466e-03,2.1573088e-03 ]
scav_coeff_mass(1:ncol,21, 2) = [                                              &
      1.7528008e-04,1.0351635e-04,6.2766829e-05,4.0139373e-05,2.8554465e-05,   &
      2.4005215e-05,2.3915669e-05,2.6090122e-05,2.8134534e-05,2.7625441e-05,   &
      2.3362038e-05,1.6746624e-05,1.1014864e-05,1.1306745e-05,1.0395807e-04,   &
      7.7384755e-04,2.0575174e-03,2.8071875e-03,2.7908872e-03,2.4753034e-03,   &
      2.2468198e-03,2.3964051e-03 ]
scav_coeff_mass(1:ncol,21, 3) = [                                              &
      1.3032880e-04,7.8050572e-05,4.8769346e-05,3.3190354e-05,2.6108970e-05,   &
      2.4346696e-05,2.5494384e-05,2.7220826e-05,2.7262641e-05,2.4248640e-05,   &
      1.8762751e-05,1.3731080e-05,2.3540348e-05,1.4681551e-04,6.6029641e-04,   &
      1.6313132e-03,2.4795467e-03,2.7369353e-03,2.5703267e-03,2.3629731e-03,   &
      2.4633539e-03,3.3358486e-03 ]
scav_coeff_mass(1:ncol,21, 4) = [                                              &
      9.3401328e-05,5.7538256e-05,3.8123289e-05,2.8686795e-05,2.5372597e-05,   &
      2.5463114e-05,2.6594594e-05,2.6643724e-05,2.4288488e-05,1.9936544e-05,   &
      1.8340396e-05,4.5515510e-05,1.9961925e-04,6.5535631e-04,1.4367372e-03,   &
      2.2090285e-03,2.5946042e-03,2.5830911e-03,2.4612501e-03,2.6058505e-03,   &
      3.5491240e-03,6.3995666e-03 ]
scav_coeff_mass(1:ncol,21, 5) = [                                              &
      6.6535252e-05,4.3279183e-05,3.1541287e-05,2.6732520e-05,2.5776293e-05,   &
      2.6223505e-05,2.6053515e-05,2.4110183e-05,2.1455832e-05,2.6829258e-05,   &
      7.5249797e-05,2.5534814e-04,6.8050955e-04,1.3452775e-03,2.0259457e-03,   &
      2.4478821e-03,2.5515784e-03,2.5422126e-03,2.7996007e-03,3.9435675e-03,   &
      7.2265035e-03,1.5168506e-02 ]
scav_coeff_mass(1:ncol,22, 1) = [                                              &
      2.3156777e-04,1.3585274e-04,8.1170971e-05,5.0136448e-05,3.3192574e-05,   &
      2.5043735e-05,2.2592101e-05,2.3642986e-05,2.6054297e-05,2.7343561e-05,   &
      2.5231953e-05,1.9449422e-05,1.2788074e-05,8.6945987e-06,8.4040313e-06,   &
      8.8889600e-05,1.4468730e-03,3.0251956e-03,3.3432472e-03,3.0119083e-03,   &
      2.5902814e-03,2.4029527e-03 ]
scav_coeff_mass(1:ncol,22, 2) = [                                              &
      1.8527095e-04,1.0939330e-04,6.6187340e-05,4.2028279e-05,2.9402523e-05,   &
      2.4053824e-05,2.3287504e-05,2.4898399e-05,2.6576247e-05,2.6023309e-05,   &
      2.2088528e-05,1.6040646e-05,1.0871527e-05,1.1548706e-05,1.0624810e-04,   &
      8.1631053e-04,2.2440446e-03,3.1318983e-03,3.1482536e-03,2.7974770e-03,   &
      2.5176235e-03,2.6272062e-03 ]
scav_coeff_mass(1:ncol,22, 3) = [                                              &
      1.3772091e-04,8.2365801e-05,5.1202591e-05,3.4393205e-05,2.6418529e-05,   &
      2.3938065e-05,2.4497897e-05,2.5816067e-05,2.5727036e-05,2.2916272e-05,   &
      1.7894237e-05,1.3405421e-05,2.3969208e-05,1.5314899e-04,7.0387457e-04,   &
      1.7770311e-03,2.7508573e-03,3.0734041e-03,2.8989422e-03,2.6499019e-03,   &
      2.7093267e-03,3.5612203e-03 ]
scav_coeff_mass(1:ncol,22, 4) = [                                              &
      9.8582809e-05,6.0483578e-05,3.9640602e-05,2.9207916e-05,2.5129638e-05,   &
      2.4615656e-05,2.5324222e-05,2.5202721e-05,2.2975258e-05,1.9005615e-05,   &
      1.8001614e-05,4.6994912e-05,2.1106425e-04,7.0451213e-04,1.5689061e-03,   &
      2.4462674e-03,2.9037811e-03,2.9048456e-03,2.7551834e-03,2.8634449e-03,   &
      3.7864224e-03,6.6551643e-03 ]
scav_coeff_mass(1:ncol,22, 5) = [                                              &
      6.9981710e-05,4.5093472e-05,3.2250139e-05,2.6627101e-05,2.5050415e-05,   &
      2.5069424e-05,2.4709959e-05,2.2849985e-05,2.0553105e-05,2.6820626e-05,   &
      7.8740479e-05,2.7254654e-04,7.3615275e-04,1.4731772e-03,2.2431651e-03,   &
      2.7338469e-03,2.8609875e-03,2.8365528e-03,3.0653482e-03,4.1933884e-03,   &
      7.5001904e-03,1.5545671e-02 ]

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_impc_scav_dust_init

! ----------------------------------------------------------------------
SUBROUTINE ukca_impc_scav_dust_dealloc()
! ----------------------------------------------------------------------
! Purpose:
! -------
! Deallocate scavenging coefficient lookup tables
! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_IMPC_SCAV_DUST_DEALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(scav_coeff_num)) DEALLOCATE (scav_coeff_num)
IF (ALLOCATED(scav_coeff_mass)) DEALLOCATE (scav_coeff_mass)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_impc_scav_dust_dealloc

END MODULE ukca_impc_scav_dust_mod
