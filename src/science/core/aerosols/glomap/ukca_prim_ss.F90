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
!    Calculate emissions of sea salt aerosol as sectional representation
!    of Gong-Monahan and add to soluble accum and coarse modes,
!    changing ND, MDT and MD accordingly (only sea surface boxes).
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
MODULE ukca_prim_ss_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_PRIM_SS_MOD'

CONTAINS

SUBROUTINE ukca_prim_ss(row_length, rows, model_levels,                        &
       verbose, land_fraction, seaice_frac, u_scalar_10m,                      &
       aer_mas_primss, aer_num_primss)

!----------------------------------------------------------------------
!
! Calculate emissions of sea salt aerosol as sectional representation
! of Gong-Monahan and add to soluble accum and coarse modes.
! changing ND, MDT and MD accordingly (only sea surface boxes)
!
! Purpose
! -------
! Use flux parametrisation (Smith,93; Smith,98; Monahan, 86;
! Gong-Monahan, 03) combined with
! wind stress to calculate flux of sea salt aerosol.
!
! Inputs
! ------
! ROW_LENGTH    : No of columns
! ROWS          : No of
! MODEL LEVELS  : No of model levels
! LAND_FRACTION : Fraction of horizontal gridbox area covered by land
! SEAICE_FRAC   : Fraction of horizontal gridbox area containing seaice
! U_SCALAR_10M  : Scalar wind at 10m (ms-1)
! VERBOSE       : Switch for level of verbosity
!
! Outputs
! -------
! AER_MAS_PRIMSS : Updated mass budget terms
! AER_NUM_PRIMSS : Updated number budget terms
!
! Local Variables
! ---------------
! S98FLUX   : Smith sea salt aerosol ems rate [dF/dr] (m-2 s-1 um-1)
! M86FLUX   : Gong/Monahan ssalt aerosol ems rate [dF/dr] (m-2 s-1 um-1)
! COMBFLUX  : Combined Smith-Monahan [dF/dr] (m-2 s-1 um-1)
! DFDR      : Chosen (from above 3) [dF/dr] (m-2 s-1 um-1)
! FLUX      : Sea salt aerosol emissions rate [(dF/dr)*deltar] (m-2 s-1)
! BOXFLUX   : Grid box sea salt aerosol flux [(dF/dr)*deltar] (m-2 s-1)
! DELN      : Change in aerosol bin no. dens. due to ssalt ems (cm^-3)
! DELN_MFLUX : Change in grid box aerosol number density (kg/m2/s)
! BET       : Parameter in Monahan formulation of ssalt aerosol ems
! AGONG     : Parameter in Gong extension of Monahan flux formulation
! MODEMT    : Tmpry varble -- which mode to emit bin-seasalt into
! CUTOFF    : Smallest dry radius that gong-monahan scheme is applicable
!             use 35nm cut-off for r at rh=80 (~17.5 nm cutoff for dryr)
! MBSMALL   : Lower bin edge for smallest sea-salt emission bin
! MBLARGE   : Upper bin edge for largest  sea-salt emission bin
! NBINS     : Number of bins for sea-salt distribution
! DELTAR    : Difference in dry radii edges for each bin (microns)
! MBLO/MBMID/MBHI : Dry mass of lower/mid/upper particle bin edge (m)
! DRLO/DRMID/DRHI : Dry radius of lower/mid/upper particle bin edge (m)
! MODEMT    : Index of mode to emit this bin-resolved ems flux into
! DUM       : Dummy variable in calculation of DRLO,DRHI,DELTAR etc.
!
! References
! ----------
! Smith 1993, QJRMS, 119, 809-824
! Smith 1998, JAS, 29, S189-190
! Monahan 1986, Oceanic Whitecaps, 167-174
! Gong, 2003, Global Biogeochemical Cycles, 17(4), 8-1,6.
!
! Inputted by module UKCA_UM_LEGACY_MOD, CHEMISTRY_CONSTANTS
! ----------------------------------------------------------
! avogadro   : Avogadro's constant (molecules per mole)
! boltzmann  : Boltzmann constant
! rgas       : Specific gas constant for dry air
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! MM        : Molar mass of each of the aerosol components (kg/mol)
! RHOCOMP   : Mass density of each of the aerosol components (kgm^-3)
! DDPLIM0   : Lower limit for dry diameter in mode (m)
!                                              or carry out process
! CP_CL     : index of cpt in which sea-salt mass is stored
!
!--------------------------------------------------------------------
USE ukca_um_legacy_mod,             ONLY: pi, rgas => r
USE chemistry_constants_mod,        ONLY: avogadro, boltzmann
USE ukca_mode_setup,                ONLY: mm, cp_cl, rhocomp, nmodes,          &
                                          ddplim0
USE ukca_config_specification_mod,  ONLY: glomap_config
USE ereport_mod,                    ONLY: ereport
USE yomhook,                        ONLY: lhook, dr_hook
USE parkind1,                       ONLY: jprb, jpim
USE umPrintMgr,                     ONLY: umPrint

USE errormessagelength_mod,         ONLY: errormessagelength

USE ukca_ingridg_mod,               ONLY: ukca_ingridg
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: row_length
! No of rows
INTEGER, INTENT(IN) :: rows
! No of columns
INTEGER, INTENT(IN) :: model_levels
! No of model levels
INTEGER, INTENT(IN) :: verbose
! Verbosity switch

REAL, INTENT(IN)    :: seaice_frac(row_length,rows)
! Sea ice fraction
REAL, INTENT(IN)    :: land_fraction(row_length,rows)
! Land fraction
REAL, INTENT(IN)    :: u_scalar_10m(row_length,rows)
! Scalar 10m wind
REAL, INTENT(IN OUT) :: aer_mas_primss(row_length,rows,                        &
                                          model_levels,nmodes)
! mass ems flux
REAL, INTENT(IN OUT) :: aer_num_primss(row_length,rows,                        &
                                          model_levels,nmodes)
! number ems flux

! Local variables
INTEGER :: jv
INTEGER :: modemt
INTEGER, PARAMETER :: nbins=20
REAL    :: mbmid(nbins)
REAL    :: mblo(nbins)
REAL    :: mbhi(nbins)
REAL    :: drmid(nbins)
REAL    :: deltar(nbins)
REAL    :: drhi(nbins)
REAL    :: drlo(nbins)
REAL    :: flux(row_length,rows)
REAL    :: boxflux(row_length,rows)
REAL    :: deln_mflux(row_length,rows)
REAL    :: bet
REAL    :: agong
REAL    :: s98flux(row_length,rows)
REAL    :: m86flux(row_length,rows)
REAL    :: combflux(row_length,rows)
REAL    :: dfdr(row_length,rows)
REAL    :: dum
REAL    :: mbsmall       ! grid limits (molecules per particle)
REAL    :: mblarge       !
REAL, PARAMETER :: mbsmall_original = 156.0
REAL, PARAMETER :: mblarge_original = 4.6e13
!! changed settings above to match bin -- realised that difference in
!! sea-salt emissions between bin and mode was due to bin upper-limit
!! for sea-salt emissions being ~ 7 microns, whereas mode is ~14 microns
!! Setting MBSMALL,MBLARGE as above results in equivalent emissions
!! grid being used for GLOMAP-mode as in the standard GLOMAP_bin run.
REAL, PARAMETER :: mbsmall_corr_rho = 211.0
REAL, PARAMETER :: mblarge_corr_rho = 6.22e13
! For use with corrected sea-salt density
!
! .. use 35 nm cut-off for r at rh=80 (approx 17.5 nm cutoff for dryr)
REAL, PARAMETER :: cutoff=0.0175e-6

! Molar mass of dry air (kg/mol)
REAL :: mm_da  ! =avogadro*boltzmann/rgas
! Scaling factor for emssion of sea salt
REAL :: scale_emiss

INTEGER, PARAMETER :: i_method_smith = 1
INTEGER, PARAMETER :: i_method_monahan = 2
INTEGER, PARAMETER :: i_method_combined = 3
INTEGER :: i_ems_method             ! Defines emission method
CHARACTER(LEN=errormessagelength) :: cmessage       ! Error message
INTEGER           :: errcode        ! Error code


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PRIM_SS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Molar mass of dry air (kg/mol)
mm_da = avogadro*boltzmann/rgas

! .. Set emission method
i_ems_method = i_method_monahan

IF (glomap_config%l_fix_nacl_density) THEN
  mbsmall = mbsmall_corr_rho
  mblarge = mblarge_corr_rho
ELSE
  mbsmall = mbsmall_original
  mblarge = mblarge_original
END IF
dum=3.0*mm(cp_cl)/(avogadro*rhocomp(cp_cl)*4.0*pi)

CALL ukca_ingridg(nbins,mbsmall,mblarge,mbmid,mblo,mbhi)

DO jv=1,nbins
  deltar(jv)=((dum*mbhi(jv))**(1.0/3.0)-                                       &
              (dum*mblo(jv))**(1.0/3.0))*1.0e6 ! in microns
  drmid(jv)=(dum*mbmid(jv))**(1.0/3.0)         ! in m
  drlo(jv) =(dum*mblo(jv ))**(1.0/3.0)         ! in m
  drhi(jv) =(dum*mbhi(jv ))**(1.0/3.0)         ! in m
END DO

! Determine the scaling factor for emission of sea salt
IF (glomap_config%l_ukca_scale_sea_salt_ems) THEN
  scale_emiss = glomap_config%sea_salt_ems_scaling
ELSE
  ! No scaling.
  scale_emiss = 1.0
END IF

! .. Loop over sea-salt emission bins
DO jv = 1,nbins
  IF (drhi(jv) > cutoff) THEN       ! if bin size is > cutoff for flux
    !                 (cutoff is lowest dry radius for the Gong-Monahan scheme)

    modemt=0
    IF (drmid(jv) <  0.5*ddplim0(4)) modemt=3    ! soluble accum. mode
    IF (drmid(jv) >= 0.5*ddplim0(4)) modemt=4    ! soluble coarse mode

    ! .. now emit according to widths as set in DDPLIM0 rather than
    ! .. as hard-coded at 1 micron dry diameter (equivalent for usual setup)

    IF (modemt > 0) THEN

      IF (i_ems_method == i_method_smith .OR.                                  &
          i_ems_method == i_method_combined) THEN
        s98flux(:,:) = 0.0
        WHERE (land_fraction < 0.999 .AND. seaice_frac < 0.999)
          ! .. Smith et al. (1998)
          s98flux(:,:) = 1.4367*0.2*(u_scalar_10m(:,:)**3.5)*                  &
            EXP(-1.5*(LOG(2.0*drmid(jv)/2.088e-6 ))**2) +                      &
            1.4367 * 6.8e-3*(u_scalar_10m(:,:)**3  )*                          &
            EXP(-1.0*(LOG(2.0*drmid(jv)/20.88e-6))**2)
        END WHERE
      END IF

      IF (i_ems_method == i_method_monahan .OR.                                &
          i_ems_method == i_method_combined) THEN
        m86flux(:,:) = 0.0
        bet = (0.433 - LOG10(drmid(jv)*2.0e6))/0.433
        agong = 4.7*(1.0 + 30.0*(drmid(jv)*2.0e6))**(-0.017*                   &
               (drmid(jv)*2.0e6)**(-1.44))
        WHERE (land_fraction < 0.999 .AND. seaice_frac < 0.999)
          ! ..  Gong (2003) ---- extension from Monahan et al. (1986)
          !          (factor of 2.0 to convert dF/dr(80) to dF/dr(dry) )
          m86flux(:,:) = 2.0*1.373*(u_scalar_10m(:,:)**3.41)*                  &
                         (drmid(jv)*2.0e6)**(-agong)*                          &
                         (1.0+0.057*((2.0e6*drmid(jv))**(3.45)))*              &
                         10.0**(1.607*EXP(-bet*bet))
        END WHERE
      END IF


      ! .. Choose Smith, Monahan, or combined flux
      IF (i_ems_method == i_method_smith) THEN
        dfdr(:,:) = s98flux(:,:)
        IF (verbose >= 2) CALL umPrint(                                        &
            'Smith sea-salt flux scheme selected',                             &
            src='ukca_prim_ss')
      ELSE IF (i_ems_method == i_method_monahan) THEN
        dfdr(:,:) = m86flux(:,:)
        IF (verbose >= 2) CALL umPrint(                                        &
            'Monahan sea-salt flux scheme selected',                           &
            src='ukca_prim_ss')
      ELSE IF (i_ems_method == i_method_combined) THEN
        ! .. If using a combination of Smith and Monahan, then the
        !     recommendation is to use the larger of the two fluxes
        !     at any radius
        IF (verbose >= 2) CALL umPrint(                                        &
            'Combined sea-salt flux scheme selected',                          &
            src='ukca_prim_ss')
        WHERE (s98flux > m86flux)
          combflux = s98flux
        ELSE WHERE
          combflux = m86flux
        END WHERE
        dfdr(:,:) = combflux(:,:)
      ELSE
        cmessage = ' No method for sea-salt flux specified'
        errcode = 1
        CALL ereport('ukca_prim_ss',errcode,cmessage)
      END IF

      ! .. Calculate sea-salt emission flux in particles m-2 s-1
      IF (drlo(jv) > cutoff) THEN
        flux(:,:) = dfdr(:,:)*deltar(jv)
      ELSE
        flux(:,:) = dfdr(:,:)*(drhi(jv) - cutoff)*1.0e6
      END IF

      ! .. Calculate  sea-salt aerosol source (particles m-2 s-1)
      boxflux(:,:) = scale_emiss * flux(:,:) * (1.0 - seaice_frac(:,:)) *      &
                     (1.0 - land_fraction(:,:))

      ! .. Calculate change in grid box aerosol number density (kg/m2/s)
      deln_mflux(:,:) = boxflux(:,:)*mm_da/avogadro

      ! .. sum up each emitted mass into mode in kg/m2/s
      aer_mas_primss(:,:,1,modemt) =                                           &
          aer_mas_primss(:,:,1,modemt) + deln_mflux*mbmid(jv)
      ! .. sum up each emitted number into mode in equiv-kg/m2/s
      aer_num_primss(:,:,1,modemt) =                                           &
          aer_num_primss(:,:,1,modemt) + deln_mflux(:,:)

    END IF    ! if modemt > 0
  END IF     ! if drhi(jv) > cutoff
END DO      ! jv = 1,nbins

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_prim_ss
END MODULE ukca_prim_ss_mod
