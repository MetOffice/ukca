! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Calculates NO3/NH4/NaNO3 emissions (mass in each mode)
!    Uses AER_NO3NH4, based on AER_NO3NH4 from CIFS-AER
!
!  References:
!    Hauglustaine et al ACP 2014
!    Code from Samuel Remy (samuel.remy at lmd.jussieu.fr)
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and The Met Office.
!  See:  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ######################################################################

MODULE ukca_prod_no3_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_PROD_NO3_MOD'

CONTAINS

SUBROUTINE ukca_prod_no3_fine(                                                 &
                        row_length, rows, model_levels, timestep,              &
                        zghno3, t_theta_levels, rel_humid_frac_clr,            &
                        air_density, air_burden, mode_tracers, nh3_mmr,        &
                        hono2_mmr, dust_div1, dust_div2, dust_div3,            &
                        dmas_no3, dnum_no3, dmas_nh4, dnum_nh4)

USE ukca_mode_setup,         ONLY: nmodes, cp_no3, cp_nh4, cp_su,              &
                                   mode_ait_sol, mode_acc_sol, mode_cor_sol,   &
                                   mode_acc_insol

USE ukca_um_legacy_mod,      ONLY: drep, rho_dust => rhop, pi, rgas => r
USE chemistry_constants_mod, ONLY: avc => avogadro, zboltz => boltzmann,       &
                                   rho_so4
USE asad_mod,                ONLY: jpctr
USE ukca_config_defs_mod,          ONLY: n_mode_tracers
USE ukca_config_specification_mod, ONLY: glomap_config, glomap_variables
USE ukca_mode_tracer_maps_mod,     ONLY: mmr_index, nmr_index
USE ereport_mod,              ONLY: ereport
USE errormessagelength_mod,   ONLY: errormessagelength

USE parkind1,   ONLY: jpim, jprb      ! DrHook
USE yomhook,    ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

!
!     Input/Output variables
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
REAL,    INTENT(IN) :: timestep
REAL,    INTENT(IN) :: zghno3 ! HNO3 uptake coefficient (default=0.193)

!     Atmosphere state variables
REAL, INTENT(IN)    :: t_theta_levels(1:row_length,1:rows,1:model_levels)
                              ! Temperature on theta levels (K)
REAL, INTENT(IN)    :: rel_humid_frac_clr(1:row_length,1:rows,1:model_levels)
                              ! Clear-sky fraction relative humidity (%)
REAL, INTENT(IN)    :: air_density(1:row_length,1:rows,1:model_levels)
                              ! Air density (kg m-3)
REAL, INTENT(IN)    :: air_burden(1:row_length,1:rows,1:model_levels)
                              ! Air burden (kg m-3)
REAL, INTENT(IN)    :: mode_tracers(1:row_length,1:rows,1:model_levels,        &
                                    1:n_mode_tracers)
                              ! All mode tracer array in kg/kg

!     Aerosol/gas mixing ratios

REAL, INTENT(IN OUT) :: nh3_mmr(1:row_length,1:rows,1:model_levels)
                              ! Ammonia mass mixing ratio
REAL, INTENT(IN OUT) :: hono2_mmr(1:row_length,1:rows,1:model_levels)
                              ! Nitric acid mass mixing ratio
REAL, INTENT(IN)     :: dust_div1(1:row_length,1:rows,1:model_levels)
                              ! CLASSIC Dust MMR in bin 1
REAL, INTENT(IN)     :: dust_div2(1:row_length,1:rows,1:model_levels)
                              ! CLASSIC Dust MMR in bin 2
REAL, INTENT(IN)     :: dust_div3(1:row_length,1:rows,1:model_levels)
                              ! CLASSIC Dust MMR in bin 3
REAL, INTENT(IN OUT) :: dmas_no3(1:row_length,1:rows,                          &
                                 1:model_levels,1:nmodes)
                              ! NO3 mass tendency (kg m-2 s-1)
REAL, INTENT(IN OUT) :: dnum_no3(1:row_length,1:rows,                          &
                                 1:model_levels,1:nmodes)
                              ! NO3 number tendency (eq-kg m-2 s-1)
REAL, INTENT(IN OUT) :: dmas_nh4(1:row_length,1:rows,                          &
                                 1:model_levels,1:nmodes)
                              ! NH4 mass tendency (kg m-2 s-1)
REAL, INTENT(IN OUT) :: dnum_nh4(1:row_length,1:rows,                          &
                                 1:model_levels,1:nmodes)
                              ! NH4 number tendency (eq-kg m-2 s-1)

!*            LOCAL VARIABLES
!              ---------------

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: ddplim0(:)
REAL,    POINTER :: ddplim1(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: rhocomp(:)
REAL,    POINTER :: sigmag(:)
REAL,    POINTER :: x(:)

! Molar masses
REAL, PARAMETER :: zmwnh3  = 17.0e-3 !Kg/mol
REAL, PARAMETER :: zmwhno3 = 63.0e-3
REAL, PARAMETER :: zmwnh4  = 18.0e-3
REAL, PARAMETER :: zmwno3  = 62.0e-3
REAL, PARAMETER :: zmwso4  = 96.0e-3

! dummy variables
INTEGER :: i, j, k, ifirst, icp, imode, updmode, idx

! Store the tracer index for mass
INTEGER :: midx(1:nmodes,1:glomap_variables%ncp)

INTEGER :: nidx(1:nmodes)         ! Store the tracer index for number

! Variables for main primary production calculation
REAL :: t_inv, t_log, rh_local, rh_inv, aird_local
REAL :: drh, kps, kpl1, kpl2, kpl3, kpl
REAL :: no3_mol_tot, nh4_mol_tot, so4_mol_tot, hno3_mol_tot, nh3_mol_tot
REAL :: ztn, zts, zta, zso4, ztadisp, ztnta
REAL :: nh4plus, nh4inso4tot, nh3inso4
REAL :: zwrk1, zwrk2, nh3_mol_up, hno3_mol_eq, nh3_mol_eq
REAL :: hno3_mol_tdep, nh3_mol_tdep, no3_mol_tdep, nh4_mol_tdep
REAL :: ptno3_tot, pt_no3_mol
REAL :: ptnh4_tot,zhno3_test, zthno3_tmp, znh3_test, ztnh3_tmp
REAL :: dtno3, dtno3lim, dtnh4, dtnh4lim

! Variables for tendency conversion from kg/kg/s to model units
REAL :: mm_da, dnumtot, test_tot_mmr, minnum, tmpnum, vol2num

! Variables for NH4/NO3 tendency allocation between aitken and accum modes
! Parameters for particle growth (C1-C4) from Gerber, H. E. (1985),
! Relative-humidity parameterization of the Navy Aerosol Model (NAM),
! NRL Rep. 8956, Naval Res. Lab., Washington, D. C.,
REAL :: sixovrpix(nmodes)
REAL :: dens2num, meanfac, tmpdiam
REAL :: zghno3_dd, zrh1580, zkn1, zkn2, zkn3
REAL :: mode_diam_dry, mode_diam_wet, mode_nconc
REAL :: zntotdust1, zntotdust2, zntotdust3
REAL :: zdg0, zkn, totrate, tau_tot, tau_fac, nk1, nk2 , nk3
REAL, PARAMETER :: zrgas = 8.314         ! Gas constant [J/mol/K]
REAL, PARAMETER :: zmgas = 29.0e-3        ! Molecular weight dry air [Kg mol-1]
REAL, PARAMETER :: zdmolec = 4.5e-10     ! Molec diameter air [m]
REAL, PARAMETER :: zratiohno3 = (63.0e-3+29e-3) / 63.0e-3
REAL, PARAMETER :: minconc = 1.0e-40      ! Min mass concentration
REAL, PARAMETER :: c1 = 0.6628           ! Parameters for hygroscopic growth
REAL, PARAMETER :: c2 = 3.082            ! using Gerber's scheme
REAL, PARAMETER :: c3 = 1.1658e-13
REAL, PARAMETER :: c4 = -1.428
REAL, PARAMETER :: zghno315_dd = 1.0e-5
REAL, PARAMETER :: zghno330_dd = 1.0e-4
REAL, PARAMETER :: zghno370_dd = 6.0e-4
REAL, PARAMETER :: zghno380_dd = 1.05e-3
REAL, PARAMETER :: athird = 1.0/3.0
REAL, PARAMETER :: drep_dust_2bin(2) = [ 0.299287e-06, 5.806452e-06 ]
REAL, PARAMETER :: drep_dust_6bin(6) = [ 0.881697e-07, 0.278817e-06,           &
                                          0.881697e-06, 0.278817e-05,          &
                                          0.881697e-05, 0.278817e-04 ]
                                         ! Number mean diameters for dust

! Local arrays
REAL :: mode_khno3(row_length,rows,model_levels,2)
REAL :: zdghno3(row_length,rows,model_levels)
REAL :: zmfp(row_length,rows,model_levels)
REAL :: no3_mol_mode(row_length,rows,model_levels,3)
REAL :: nh4_mol_mode(row_length,rows,model_levels,3)
REAL :: so4_mol_mode(row_length,rows,model_levels,3)
REAL :: nh4inso4(row_length,rows,model_levels,3)
REAL :: nh4inno3(row_length,rows,model_levels,3)
REAL :: relrate(row_length,rows,model_levels,3)

! Local tendencies (kg(sub)/kg(air)/s)
! for main primary production calculation
REAL :: tot_mode_mmr(row_length,rows,model_levels,3)
REAL :: tot_mode_vmr(row_length,rows,model_levels,nmodes)
REAL :: ptno3(row_length,rows,model_levels,3)
REAL :: ptnh4(row_length,rows,model_levels,3)
REAL :: pthno3(row_length,rows,model_levels)
REAL :: ptnh3(row_length,rows,model_levels)

INTEGER                            :: errcode    ! Error code for ereport
CHARACTER (LEN=errormessagelength) :: cmessage   ! Error message
INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PROD_NO3_FINE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
ddplim0     => glomap_variables%ddplim0
ddplim1     => glomap_variables%ddplim1
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
rhocomp     => glomap_variables%rhocomp
sigmag      => glomap_variables%sigmag
x           => glomap_variables%x

!---------------------------------------
! 0. Set up arrays for mode/component mass and number indices in mode_tracers
ifirst = jpctr + 1
midx(:,:)=-1
nidx(:)=-1
DO imode = 1, nmodes
  IF (mode(imode)) THEN
    nidx(imode) = nmr_index(imode) - ifirst + 1
    DO icp = 1,ncp
      IF (component(imode,icp)) THEN
        midx(imode,icp) = mmr_index(imode,icp) - ifirst + 1
      END IF
    END DO
  END IF
END DO

!---------------------------------------
! 1. Sum up total mass in Aitken to coarse soluble modes
tot_mode_mmr(:,:,:,:)=0.0
DO imode = mode_ait_sol, mode_cor_sol
  IF (mode(imode)) THEN
    idx=imode-mode_ait_sol+1
    DO icp = 1,ncp
      IF (component(imode,icp)) THEN
        tot_mode_mmr(:,:,:,idx) = tot_mode_mmr(:,:,:,idx) +                    &
                                  mode_tracers(:,:,:,midx(imode,icp))
      END IF ! component
    END DO ! icp
  END IF ! mode
END DO ! imode

!---------------------------------------
! 2. Calculate the time taken to reach equilibrium for the NH4NO3 reaction
!
!    ----- may want to move this to a seperate module -----
!
!    Calculate geometric mean diameter of Aitken and accumulation modes
!    and grow according to RH-dependent ammonium sulfate growth factor
!    Method from Makar (1998) JGR, 103,D 11, 13,095-13,110
!
mode_khno3(:,:,:,:) = 0.0
tot_mode_vmr(:,:,:,:) = 0.0
sixovrpix(:) = 6.0/(pi*x(:))
dens2num=zboltz/rgas

! 2.0 Calculate the diffusion constant and mean free path
DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      zdg0 = 3.0/(8.0*avc*air_density(i,j,k)*(zdmolec**2.0))
      zdghno3(i,j,k) = zdg0 * SQRT(zrgas*zmgas/2.0/pi*t_theta_levels(i,j,k)*   &
                                   zratiohno3)
      zmfp(i,j,k) = 3.0*zdghno3(i,j,k)/SQRT(8.0*zrgas*t_theta_levels(i,j,k)/   &
                                           pi/zmwhno3)
    END DO
  END DO
END DO

! 2.1 Mean wet diameter (m) and number concentration (N)
!     Use wet diameter to calculate uptake rate coefficent k_hno3
!     using Fuchs and Sutugin (1970) uptake formulation
DO imode = mode_ait_sol, mode_acc_sol
  IF (mode(imode)) THEN
    idx=imode-mode_ait_sol+1
    meanfac=EXP(0.5*LOG(sigmag(imode))*LOG(sigmag(imode)))

    ! Calculate the total modal volume
    DO icp = 1,ncp
      IF (component(imode,icp)) THEN
        tot_mode_vmr(:,:,:,imode) = tot_mode_vmr(:,:,:,imode) +                &
                                    mode_tracers(:,:,:,midx(imode,icp)) /      &
                                    rhocomp(icp)
      END IF !component
    END DO !icp

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          IF (mode_tracers(i,j,k,nidx(imode)) > 0.0) THEN
            tmpdiam = (tot_mode_vmr(i,j,k,imode) /                             &
                       mode_tracers(i,j,k,nidx(imode)) *                       &
                       sixovrpix(imode) * dens2num)**athird
            tmpdiam = MAX(ddplim0(imode),tmpdiam)
            tmpdiam = MIN(tmpdiam,ddplim1(imode))
            mode_diam_dry = tmpdiam * meanfac
            mode_nconc = mode_tracers(i,j,k,nidx(imode)) * air_density(i,j,k) /&
                         dens2num
          ELSE
            mode_diam_dry = ddplim0(imode) * meanfac
            mode_nconc = 6.0*minconc/pi/rho_so4/(mode_diam_dry**3.0)
          END IF
          ! Grow the particles with humidity using Gerber's model
          rh_local = rel_humid_frac_clr(i,j,k)
          ! LOG of -ve number is NaN. LOG of zero is -inf.
          ! Both are undesirable.
          IF (rh_local <= 0.0) rh_local = EPSILON(0.0)
          IF (rh_local > 0.98) rh_local = 0.98
          mode_diam_wet =  (c1*(mode_diam_dry**c2) /                           &
                           (c3*(mode_diam_dry**c4)-LOG10(rh_local)) +          &
                            mode_diam_dry**3.0)**athird
          ! Calculate the uptake coefficient
          zkn  = 2.0*zmfp(i,j,k) / mode_diam_wet
          nk1  = (2.0*pi*mode_diam_wet*zdghno3(i,j,k)) / (1.0+(4.0*zkn/        &
                  zghno3/3.0)*(1.0-0.47*zghno3/(1.0+zkn)))
          mode_khno3(i,j,k,idx) = mode_khno3(i,j,k,idx) + mode_nconc * nk1
        END DO
      END DO
    END DO
  END IF ! mode
END DO !imode

! Now include dust
SELECT CASE (glomap_config%i_dust_scheme)
CASE (1)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        !RH dependent HNO3 uptake coefficient based on Fairlie et al. 2010
        !and scale for SS
        zrh1580 = rel_humid_frac_clr(i,j,k)
        IF (zrh1580 < 0.15) zrh1580 = 0.15
        IF (zrh1580 > 0.80) zrh1580 = 0.80
        IF ( zrh1580 < 0.30 ) THEN
          zghno3_dd = zghno315_dd + (zghno330_dd-zghno315_dd)/                 &
                      0.15*(zrh1580-0.15)
        ELSE IF (zrh1580 > 0.30 .AND. zrh1580 < 0.70) THEN
          zghno3_dd = zghno330_dd + (zghno370_dd-zghno330_dd)/                 &
                      0.40*(zrh1580-0.30)
        ELSE IF (zrh1580 > 0.70 ) THEN
          zghno3_dd = zghno370_dd + (zghno380_dd-zghno370_dd)/                 &
                      0.10*(zrh1580-0.70)
        END IF

        ! Calculate number density
        zntotdust1 = dust_div1(i,j,k)*air_density(i,j,k)  /                    &
                     (rho_dust*4.0*pi*(drep(1)*0.5)**3/3.0)
        zntotdust2 = dust_div2(i,j,k)*air_density(i,j,k)  /                    &
                     (rho_dust*4.0*pi*(drep(2)*0.5)**3/3.0)
        zntotdust3 = dust_div3(i,j,k)*air_density(i,j,k)  /                    &
                     (rho_dust*4.0*pi*(drep(3)*0.5)**3/3.0)

        ! Calculate uptake coefficient and apportion to aitken / accum modes
        zkn1 = 2.0*zmfp(i,j,k) / drep_dust_6bin(1)
        zkn2 = 2.0*zmfp(i,j,k) / drep_dust_6bin(2)
        zkn3 = 2.0*zmfp(i,j,k) / drep_dust_6bin(3)
        nk1  = (2.0*pi*drep_dust_6bin(1)*zdghno3(i,j,k)) / (1.0+(4.0*zkn1/     &
                zghno3_dd/3.0)*(1.0-0.47*zghno3_dd/(1.0+zkn1)))
        nk2  = (2.0*pi*drep_dust_6bin(2)*zdghno3(i,j,k)) / (1.0+(4.0*zkn2/     &
                zghno3_dd/3.0)*(1.0-0.47*zghno3_dd/(1.0+zkn2)))
        nk3  = (2.0*pi*drep_dust_6bin(3)*zdghno3(i,j,k)) / (1.0+(4.0*zkn3/     &
                zghno3_dd/3.0)*(1.0-0.47*zghno3_dd/(1.0+zkn3)))
        mode_khno3(i,j,k,mode_ait_sol-1) = mode_khno3(i,j,k,mode_ait_sol-1) +  &
                                           zntotdust1 * nk1
        mode_khno3(i,j,k,mode_acc_sol-1) = mode_khno3(i,j,k,mode_acc_sol-1) +  &
                                           zntotdust2 * nk2 +                  &
                                     0.5 * zntotdust3 * nk3
      END DO
    END DO
  END DO

CASE (2)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        !RH dependent HNO3 uptake coefficient based on Fairlie et al. 2010
        !and scale for SS
        zrh1580 = rel_humid_frac_clr(i,j,k)
        IF (zrh1580 < 0.15) zrh1580 = 0.15
        IF (zrh1580 > 0.80) zrh1580 = 0.80
        IF ( zrh1580 < 0.30 ) THEN
          zghno3_dd = zghno315_dd + (zghno330_dd-zghno315_dd)/                 &
                      0.15*(zrh1580-0.15)
        ELSE IF (zrh1580 > 0.30 .AND. zrh1580 < 0.70) THEN
          zghno3_dd = zghno330_dd + (zghno370_dd-zghno330_dd)/                 &
                      0.40*(zrh1580-0.30)
        ELSE IF (zrh1580 > 0.70 ) THEN
          zghno3_dd = zghno370_dd + (zghno380_dd-zghno370_dd)/                 &
                      0.10*(zrh1580-0.70)
        END IF

        ! Calculate number density
        zntotdust1 = dust_div1(i,j,k)*air_density(i,j,k)  /                    &
                     (rho_dust*4.0*pi*(drep(1)*0.5)**3/3.0)

        ! Calculate uptake coefficient and add to accumulation mode
        zkn1 = 2.0*zmfp(i,j,k) / drep_dust_2bin(1)
        nk1  = (2.0*pi*drep_dust_2bin(1)*zdghno3(i,j,k)) / (1.0+(4.0*zkn1/     &
                zghno3_dd/3.0)*(1.0-0.47*zghno3_dd/(1.0+zkn1)))
        mode_khno3(i,j,k,mode_acc_sol-1) =                                     &
                              mode_khno3(i,j,k,mode_acc_sol-1) +               &
                              0.145 * zntotdust1 * nk1
      END DO
    END DO
  END DO

CASE (3)
  imode = mode_acc_insol
  IF (mode(imode)) THEN
    meanfac=EXP(0.5*LOG(sigmag(imode))*LOG(sigmag(imode)))

    ! Calculate the total modal volume
    DO icp = 1,ncp
      IF (component(imode,icp)) THEN
        tot_mode_vmr(:,:,:,imode) = tot_mode_vmr(:,:,:,imode) +                &
                                    mode_tracers(:,:,:,midx(imode,icp))/       &
                                    rhocomp(icp)
      END IF !component
    END DO !icp

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          !RH dependent HNO3 uptake coefficient based on Fairlie et al. 2010
          !and scale for SS
          zrh1580 = rel_humid_frac_clr(i,j,k)
          IF (zrh1580 < 0.15) zrh1580 = 0.15
          IF (zrh1580 > 0.80) zrh1580 = 0.80
          IF ( zrh1580 < 0.30 ) THEN
            zghno3_dd = zghno315_dd + (zghno330_dd-zghno315_dd)/               &
                        0.15*(zrh1580-0.15)
          ELSE IF (zrh1580 > 0.30 .AND. zrh1580 < 0.70) THEN
            zghno3_dd = zghno330_dd + (zghno370_dd-zghno330_dd)/               &
                        0.40*(zrh1580-0.30)
          ELSE IF (zrh1580 > 0.70 ) THEN
            zghno3_dd = zghno370_dd + (zghno380_dd-zghno370_dd)/               &
                        0.10*(zrh1580-0.70)
          END IF
          IF (mode_tracers(i,j,k,nidx(imode)) > 0.0) THEN
            tmpdiam = (tot_mode_vmr(i,j,k,imode) /                             &
                      mode_tracers(i,j,k,nidx(imode)) *                        &
                      sixovrpix(imode) * dens2num)**athird
            tmpdiam = MAX(ddplim0(imode),tmpdiam)
            tmpdiam = MIN(tmpdiam,ddplim1(imode))
            mode_diam_dry = tmpdiam * meanfac
            mode_nconc = mode_tracers(i,j,k,nidx(imode)) * air_density(i,j,k) /&
                         dens2num
          ELSE
            mode_diam_dry = ddplim0(imode) * meanfac
            mode_nconc = 6.0*minconc/pi/rho_dust/(mode_diam_dry**3.0)
          END IF
          ! Calculate the uptake coefficient and add to accumulation mode
          zkn  = 2.0*zmfp(i,j,k) / mode_diam_dry
          nk1  = (2.0*pi*mode_diam_dry*zdghno3(i,j,k)) / (1.0+(4.0*zkn/        &
                  zghno3_dd/3.0)*(1.0-0.47*zghno3_dd/(1.0+zkn)))
          mode_khno3(i,j,k,mode_acc_sol-1) =                                   &
                                mode_khno3(i,j,k,mode_acc_sol-1) +             &
                                mode_nconc * nk1
        END DO
      END DO
    END DO
  END IF ! mode

CASE DEFAULT
  errcode=-1
  cmessage = 'No dust scheme defined - needed for primary NO3 production'
  CALL ereport (ModuleName//':'//RoutineName, errcode, cmessage)
END SELECT

! Set a minimum value for KHNO3 to avoid dividing by zero
WHERE (mode_khno3(:,:,:,:) < 1.0e-40)
  mode_khno3(:,:,:,:) = 1.0e-40
END WHERE

!---------------------------------------
! 3. Calculate equilibrium concentrations of HNO3 and NH3 and apportion to
!    ammonium, nitrate, and sulfate
!    Scheme follows Hauglustaine (2014) Atmos. Chem. Phys., 14, 11031-11063

! Initialise arrays to hold mass mixing ratio tendencies
ptno3(:,:,:,:) = 0.0
ptnh4(:,:,:,:) = 0.0
pthno3(:,:,:) = 0.0
ptnh3(:,:,:) = 0.0
no3_mol_mode(:,:,:,:)=0.0
nh4_mol_mode(:,:,:,:)=0.0
so4_mol_mode(:,:,:,:)=0.0
nh4inso4(:,:,:,:) = 0.0
nh4inno3(:,:,:,:) = 0.0
relrate(:,:,:,:) = 0.0

updmode = mode_ait_sol - 1 ! any remaining NH4 is emitted into Aitken mode

! Begin loop over gridcells
DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length

      ! 3.1 Parameters related to temperature / humidity
      t_inv = 1.0 / t_theta_levels(i,j,k)
      t_log = LOG(t_theta_levels(i,j,k))
      rh_local = rel_humid_frac_clr(i,j,k)
      aird_local = air_density(i,j,k)
      IF (rh_local < 0.0) rh_local = 0
      IF (rh_local > 0.98) rh_local = 0.98
      rh_inv = 1.0 - rh_local

      ! 3.2 Equilibrium constant for [NH4.NO3] based on Mozurkewich, 1993
      drh  = EXP(723.7 * t_inv + 1.6954) * 1.0e-2
      kps  = EXP( 118.87 - 24084.0 * t_inv - 6.025 * t_log)
      kpl1 = EXP(-135.94 + 8763.0  * t_inv + 19.12 * t_log)
      kpl2 = EXP(-122.65 + 9969.0  * t_inv + 16.22 * t_log)
      kpl3 = EXP(-182.61 + 13875.0 * t_inv + 24.46 * t_log)
      kpl  = kps
      IF (rh_local >= drh ) THEN
        kpl = (kpl1-kpl2*rh_inv+kpl3*rh_inv**2.0) * rh_inv**1.75 * kpl
      END IF

      ! 3.3 Convert aerosol mass mixing ratios to molar concentrations
      no3_mol_tot = 0.0
      nh4_mol_tot = 0.0
      so4_mol_tot = 0.0
      DO imode = mode_ait_sol, mode_cor_sol
        idx = imode - mode_ait_sol + 1
        no3_mol_mode(i,j,k,idx) = mode_tracers(i,j,k,midx(imode,cp_no3))*      &
                                  aird_local*1.0e9/zmwno3
        nh4_mol_mode(i,j,k,idx) = mode_tracers(i,j,k,midx(imode,cp_nh4))*      &
                                  aird_local*1.0e9/zmwnh4
        so4_mol_mode(i,j,k,idx) = mode_tracers(i,j,k,midx(imode,cp_su))*       &
                                  aird_local*1.0e9/zmwso4
        no3_mol_tot = no3_mol_tot + no3_mol_mode(i,j,k,idx)
        nh4_mol_tot = nh4_mol_tot + nh4_mol_mode(i,j,k,idx)
        so4_mol_tot = so4_mol_tot + so4_mol_mode(i,j,k,idx)
      END DO
      hno3_mol_tot = hono2_mmr(i,j,k)*aird_local*1.0e9/zmwhno3
      nh3_mol_tot = nh3_mmr(i,j,k)*aird_local*1.0e9/zmwnh3
      ztn = hno3_mol_tot + no3_mol_tot
      zta = nh4_mol_tot + nh3_mol_tot
      zts = so4_mol_tot

      ! 3.4 SO4 state. Metzger et al. 2002
      zso4 = 2.0
      IF (zts > 0.5  *zta) THEN
        zso4 = 1.5
      END IF
      IF (zts > zta) THEN
        zso4 = 1.0
      END IF

      ! 3.5 Neutralize SO4 using NH4 and if necessary NH3
      nh4plus = zta
      nh4inso4tot = 0.0
      DO imode = mode_ait_sol, mode_cor_sol
        idx = imode - mode_ait_sol + 1
        nh4inso4(i,j,k,idx) = MIN(nh4plus,zso4*so4_mol_mode(i,j,k,idx))
        nh4plus=MAX(0.0,nh4plus-nh4inso4(i,j,k,idx))
        nh4inso4tot = nh4inso4tot + nh4inso4(i,j,k,idx)
      END DO
      ztadisp = MAX(0.0,zta-nh4inso4tot)
      nh3inso4 = MAX(0.0,nh4inso4tot-nh4_mol_tot)

      ! 3.6 Determine whether NO3 is formed or disassociates
      ztnta = ztn * ztadisp
      ptno3_tot = 0.0
      IF (ztnta > kpl ) THEN
        ! 3.6.1 Equilibrium NH4.NO3 concentrations
        zwrk1  = (ztadisp+ztn)**2.0-4.0*(ztnta-kpl)
        zwrk1  = MAX(zwrk1,0.0)
        zwrk2  = 0.5*(ztadisp+ztn-SQRT(zwrk1))
        zwrk2  = MAX(zwrk2,0.0)
        zwrk2  = MIN(zwrk2,ztn)
        zwrk2  = MIN(zwrk2,ztadisp)

        ! 3.6.2 Calculate time dependence of uptake in Aitken
        !       and Accumulation modes = khno3 * number conc
        totrate = 0.0
        DO imode = mode_ait_sol, mode_acc_sol
          idx = imode - mode_ait_sol + 1
          totrate = totrate + mode_khno3(i,j,k,idx)
        END DO
        DO imode = mode_ait_sol, mode_acc_sol
          idx = imode - mode_ait_sol + 1
          relrate(i,j,k,idx) = mode_khno3(i,j,k,idx)/totrate
        END DO
        tau_tot = 1.0/totrate
        tau_fac = 1.0-EXP(-1.0*timestep/tau_tot)

        ! 3.6.3 Time dependence gas phase concentrations
        !       Firstly update NH3 to account for SO4 neutralisation
        !       Then calculate equilibrium gas concentrations
        !       Then add time-dependence
        nh3_mol_up   = MAX(nh3_mol_tot-nh3inso4, 0.0)
        hno3_mol_eq  = MAX(ztn-zwrk2, 0.0 )
        nh3_mol_eq   = MAX(ztadisp-zwrk2, 0.0)
        hno3_mol_tdep= MAX(hno3_mol_tot-tau_fac*(hno3_mol_tot-hno3_mol_eq),0.0)
        nh3_mol_tdep = MAX(nh3_mol_up-tau_fac*(nh3_mol_up-nh3_mol_eq),0.0)

        ! 3.6.4 Update aerosol phase - firstly in moles, then convert to MMRs
        no3_mol_tdep = MAX(ztn-hno3_mol_tdep, 0.0)
        nh4_mol_tdep = MAX(ztadisp-nh3_mol_tdep, 0.0)
        nh4plus = nh4_mol_tdep
        DO imode = mode_ait_sol, mode_cor_sol
          idx = imode - mode_ait_sol + 1
          pt_no3_mol = MAX(0.0,(no3_mol_tdep-no3_mol_tot))*                    &
                                    relrate(i,j,k,idx)/timestep
          nh4inno3(i,j,k,idx) = MIN(nh4plus,no3_mol_mode(i,j,k,idx)+           &
                                    pt_no3_mol*timestep)
          nh4plus = MAX(0.0,nh4plus-nh4inno3(i,j,k,idx))
          ptno3(i,j,k,idx) = pt_no3_mol*zmwno3/aird_local/1.0e9
          ptno3_tot = ptno3_tot+ptno3(i,j,k,idx)
        END DO
        nh4inno3(i,j,k,updmode)=nh4inno3(i,j,k,updmode)+nh4plus

      ELSE
        ! 3.6.5 Remove all nitrate and set NH4(NO3) concentration to zero
        DO imode = mode_ait_sol, mode_cor_sol
          idx = imode - mode_ait_sol + 1
          ptno3(i,j,k,idx) = -1.0*mode_tracers(i,j,k,midx(imode,cp_no3))/      &
                             timestep
          ptno3_tot = ptno3_tot + ptno3(i,j,k,idx)
        END DO
      END IF

      ! 3.7 Absorb remaining NH4 concentration differences int0
      !    total NH4 tendency
      ptnh4_tot = 0.0
      DO imode = mode_ait_sol, mode_cor_sol
        idx = imode - mode_ait_sol + 1
        ptnh4(i,j,k,idx) = (nh4inso4(i,j,k,idx)+nh4inno3(i,j,k,idx)-           &
                            nh4_mol_mode(i,j,k,idx)) *                         &
                            zmwnh4/aird_local/1.0e9/timestep
        ptnh4_tot = ptnh4_tot + ptnh4(i,j,k,idx)
      END DO

      ! 3.8 Update gas phase tendencies and check for negative masses
      pthno3(i,j,k) = -1.0*ptno3_tot*zmwhno3/zmwno3
      ptnh3(i,j,k)  = -1.0*ptnh4_tot*zmwnh3/zmwnh4

      zhno3_test = hono2_mmr(i,j,k) + pthno3(i,j,k)*timestep
      IF (zhno3_test < 0.0) THEN
        zthno3_tmp = pthno3(i,j,k)
        pthno3(i,j,k)  = -1.0*hono2_mmr(i,j,k)/timestep
        zthno3_tmp = pthno3(i,j,k) -zthno3_tmp
        DO imode = mode_ait_sol, mode_acc_sol
          idx = imode - mode_ait_sol + 1
          ptno3(i,j,k,idx) = ptno3(i,j,k,idx)-zthno3_tmp*zmwno3/zmwhno3*       &
                             relrate(i,j,k,idx)
        END DO
      END IF

      znh3_test = nh3_mmr(i,j,k) + ptnh3(i,j,k) * timestep
      IF (znh3_test < 0.0) THEN
        ztnh3_tmp = ptnh3(i,j,k)
        ptnh3(i,j,k) = -1.0*nh3_mmr(i,j,k)/timestep
        ztnh3_tmp = ptnh3(i,j,k)-ztnh3_tmp
        DO imode = mode_ait_sol, mode_acc_sol
          idx = imode - mode_ait_sol + 1
          ptnh4(i,j,k,idx) = ptnh4(i,j,k,idx)-ztnh3_tmp*zmwnh4/zmwnh3*         &
                             relrate(i,j,k,idx)
        END DO
      END IF
    END DO
  END DO
END DO

! 4. Convert NO3 and NH4 mass tendencies to kg/m2/s from kg/kg/s
!    Then fill output mass tendency arrays
!    Currently removing all mass from modes if negative emissions
!    and removing number as a ratio of increment to total mass
!    to conserve diameter
DO imode = mode_ait_sol, mode_cor_sol
  idx = imode - mode_ait_sol + 1
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        ! Local copies of tendencies
        dtno3 = ptno3(i,j,k,idx) * air_burden(i,j,k)
        dtno3lim = -1.0*mode_tracers(i,j,k,midx(imode,cp_no3))*                &
                                    air_burden(i,j,k)/timestep
        dtnh4 = ptnh4(i,j,k,idx) * air_burden(i,j,k)
        dtnh4lim = -1.0*mode_tracers(i,j,k,midx(imode,cp_nh4))*                &
                                    air_burden(i,j,k)/timestep

        dmas_no3(i,j,k,imode) = dmas_no3(i,j,k,imode) + MAX(dtno3lim, dtno3)
        dmas_nh4(i,j,k,imode) = dmas_nh4(i,j,k,imode) + MAX(dtnh4lim, dtnh4)

        ! CHECK NUMBER 1
        ! If all aerosol mass is removed from mode (test_tot_mmr)
        ! then we must also remove number concentration
        test_tot_mmr = tot_mode_mmr(i,j,k,idx) +                               &
                       (dmas_nh4(i,j,k,imode) + dmas_no3(i,j,k,imode)) *       &
                       timestep / air_burden(i,j,k)

        IF (test_tot_mmr < 1.0e-40) THEN
          dnumtot = -1.0*mode_tracers(i,j,k,nidx(imode))*                      &
                        air_burden(i,j,k)/timestep
          ! Arbitrary apportionment of number tendency between NO3 and NH4
          dnum_no3(i,j,k,imode) = dnum_no3(i,j,k,imode) + dnumtot*0.5
          dnum_nh4(i,j,k,imode) = dnum_nh4(i,j,k,imode) + dnumtot*0.5
        END IF

      END DO
    END DO
  END DO
END DO

! 5. CHECK NUMBER 2 - Emissions mode
!    If adding NH4NO3 mass where current number concentrations
!    are negligible we also add number (assuming the growth of
!    nucleation size particles to form Aitken sizes)
mm_da = avc * zboltz / rgas
DO imode = mode_ait_sol, mode_acc_sol
  vol2num = (mm_da/avc)*sixovrpix(imode)/(ddplim1(imode)**3)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length

        dtno3 = dmas_no3(i,j,k,imode)/rhocomp(cp_no3)/air_burden(i,j,k)
        dtnh4 = dmas_nh4(i,j,k,imode)/rhocomp(cp_nh4)/air_burden(i,j,k)
        minnum = (MAX(0.0,dtno3) + MAX(0.0,dtnh4)) * vol2num * timestep
        tmpnum = mode_tracers(i,j,k,nidx(imode))

        IF (tmpnum < minnum) THEN
          dnumtot = (minnum-tmpnum)*air_burden(i,j,k)/timestep
          ! Arbitrary apportionment of number tendency between NO3 and NH4
          dnum_no3(i,j,k,imode) = dnum_no3(i,j,k,imode) + dnumtot*0.5
          dnum_nh4(i,j,k,imode) = dnum_nh4(i,j,k,imode) + dnumtot*0.5
        END IF

      END DO
    END DO
  END DO
END DO

! 6. Update the NH3 and HONO2 MMRs and correct for negative values
nh3_mmr(:,:,:) = nh3_mmr(:,:,:) + ptnh3(:,:,:) * timestep
hono2_mmr(:,:,:) = hono2_mmr(:,:,:) + pthno3(:,:,:) * timestep

WHERE (nh3_mmr(:,:,:) < 0.0)
  nh3_mmr(:,:,:) = 0.0
END WHERE

WHERE (hono2_mmr(:,:,:) < 0.0)
  hono2_mmr(:,:,:) = 0.0
END WHERE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE ukca_prod_no3_fine

!-----------------------------------------------------------------------

SUBROUTINE ukca_prod_no3_coarse(                                               &
                        row_length, rows, model_levels, timestep,              &
                        t_theta_levels, rel_humid_frac_clr, air_density,       &
                        air_burden, tropopause_level, mode_tracers, hono2_mmr, &
                        dust_div1, dust_div2, dust_div3, dust_div4,            &
                        dust_div5, dust_div6, dmas_no3, dmas_dust, dmas_nacl)

USE chemistry_constants_mod,  ONLY: zboltz => boltzmann
USE ukca_um_legacy_mod,       ONLY: pi, rgas => r

USE ukca_mode_setup,          ONLY: nmodes,                                    &
                                    cp_nn, cp_cl, cp_du,                       &
                                    mode_acc_sol, mode_cor_sol,                &
                                    mode_acc_insol, mode_cor_insol

USE asad_mod,                 ONLY: jpctr
USE ukca_aer_no3_mod,         ONLY: aer_no3_2bindu, aer_no3_6bindu,            &
                                    aer_no3_ukcadu
USE ereport_mod,              ONLY: ereport
USE errormessagelength_mod,   ONLY: errormessagelength

USE ukca_mode_tracer_maps_mod,     ONLY: mmr_index, nmr_index
USE ukca_config_defs_mod,          ONLY: n_mode_tracers

USE ukca_config_specification_mod, ONLY: glomap_config, glomap_variables

USE parkind1,   ONLY: jpim, jprb      ! DrHook
USE yomhook,    ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

!     Input/Output variables
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
REAL,    INTENT(IN) :: timestep

!     Atmosphere state variables
REAL, INTENT(IN)    :: t_theta_levels(1:row_length,1:rows,1:model_levels)
                              ! Temperature on theta levels (K)
REAL, INTENT(IN)    :: rel_humid_frac_clr(1:row_length,1:rows,1:model_levels)
                              ! Clear-sky fraction relative humidity (%)
REAL, INTENT(IN)    :: air_density(1:row_length,1:rows,1:model_levels)
                              ! Air density (kg m-3)
REAL, INTENT(IN)    :: air_burden(1:row_length,1:rows,1:model_levels)
                              ! Air burden (kg m-3)
REAL, INTENT(IN)    :: mode_tracers(1:row_length,1:rows,1:model_levels,        &
                                    1:n_mode_tracers)
                              ! All mode tracer array in kg/kg
INTEGER, INTENT(IN) :: tropopause_level(1:row_length, 1:rows)

!     Aerosol/gas mixing ratios
REAL, INTENT(IN) :: dust_div1(1:row_length,1:rows,1:model_levels)
                              ! CLASSIC Dust MMR in bin 1
REAL, INTENT(IN) :: dust_div2(1:row_length,1:rows,1:model_levels)
                              ! CLASSIC Dust MMR in bin 2
REAL, INTENT(IN) :: dust_div3(1:row_length,1:rows,1:model_levels)
                              ! CLASSIC Dust MMR in bin 3
REAL, INTENT(IN) :: dust_div4(1:row_length,1:rows,1:model_levels)
                              ! CLASSIC Dust MMR in bin 4
REAL, INTENT(IN) :: dust_div5(1:row_length,1:rows,1:model_levels)
                              ! CLASSIC Dust MMR in bin 5
REAL, INTENT(IN) :: dust_div6(1:row_length,1:rows,1:model_levels)
                              ! CLASSIC Dust MMR in bin 6
REAL, INTENT(IN OUT) :: hono2_mmr(1:row_length,1:rows,1:model_levels)
                              ! Nitric acid mass mixing ratio
REAL, INTENT(IN OUT) :: dmas_no3(1:row_length,1:rows,                          &
                                 1:model_levels,1:nmodes)
                              ! NO3 mass tendency (kg m-2 s-1)
REAL, INTENT(IN OUT) :: dmas_dust(1:row_length,1:rows,                         &
                                  1:model_levels,1:nmodes)
                              ! UKCA dust mass tendency (kg m-2 s-1)
                              ! Only applicable if dust is on
REAL, INTENT(IN OUT) :: dmas_nacl(1:row_length,1:rows,                         &
                                  1:model_levels,1:nmodes)
                              ! Sea-salt mass tendency (kg m-2 s-1)

! Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: ddplim0(:)
REAL,    POINTER :: ddplim1(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: rhocomp(:)
REAL,    POINTER :: sigmag(:)
REAL,    POINTER :: x (:)

REAL :: dmas_hono2(1:row_length,1:rows,1:model_levels)
REAL :: nano3_mmr(1:row_length,1:rows,1:model_levels)
REAL :: seasalt_diam(1:row_length,1:rows,1:model_levels,1:2)
REAL :: seasalt_mconc(1:row_length,1:rows,1:model_levels,1:2)
REAL :: seasalt_nconc(1:row_length,1:rows,1:model_levels,1:2)
REAL, ALLOCATABLE :: dust_diam(:,:,:,:)
REAL, ALLOCATABLE :: dust_mconc(:,:,:,:)
REAL, ALLOCATABLE :: dust_nconc(:,:,:,:)
REAL :: tot_mode_vmr(1:row_length,1:rows,1:model_levels)
      ! Used to calculate modal diameter - MMR weighted by cpt density

REAL :: sixovrpix(nmodes)
REAL :: dens2num, meanfac, tmpdiam
REAL :: ems_no3, ems_dust, ems_nacl
REAL, PARAMETER :: loc_ems_eps = 1.0e-40
REAL, PARAMETER :: athird = 1.0/3.0

! Local variables
INTEGER :: imode, i, j, k, ifirst, icp, idx

! Store the tracer index for mass
INTEGER :: midx(1:nmodes,1:glomap_variables%ncp)

INTEGER :: nidx(1:nmodes)         ! Store the tracer index for number

INTEGER                            :: errcode    ! Error code for ereport
CHARACTER (LEN=errormessagelength) :: cmessage   ! Error message
INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PROD_NO3_COARSE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
ddplim0     => glomap_variables%ddplim0
ddplim1     => glomap_variables%ddplim1
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp
rhocomp     => glomap_variables%rhocomp
sigmag      => glomap_variables%sigmag
x           => glomap_variables%x

! Set up arrays for mode/component mass and number indices in mode_tracers
ifirst = jpctr + 1
midx(:,:)=-1
nidx(:)=-1
DO imode = 1, nmodes
  IF (mode(imode)) THEN
    nidx(imode) = nmr_index(imode) - ifirst + 1
    DO icp = 1,ncp ! Begin loop over components
      IF (component(imode,icp)) THEN
        midx(imode,icp) = mmr_index(imode,icp) - ifirst + 1
      END IF
    END DO
  END IF
END DO

! Sum up NaNO3 mass in accum / coarse soluble modes
nano3_mmr(:,:,:)=0.0
DO imode = mode_acc_sol, mode_cor_sol
  IF (mode(imode)) THEN
    IF (component(imode,cp_nn)) THEN
      nano3_mmr(:,:,:) = nano3_mmr(:,:,:)+                                     &
                         mode_tracers(:,:,:,midx(imode,cp_nn))
    ELSE  ! Component cp_nn
      WRITE(cmessage,'(A44,I2,A10,I2,A14)')                                    &
                          'l_ukca_coarse_no3_prod = TRUE but component ',      &
                           cp_nn,' and mode ',imode,' not turned on'
      errcode=2
      CALL ereport(RoutineName, errcode, cmessage)
    END IF
  ELSE ! Mode imode
    WRITE(cmessage,'(A39,I2,A14)') 'l_ukca_coarse_no3_prod = TRUE but mode ',  &
                                   imode,' not turned on'
    errcode=1
    CALL ereport(RoutineName, errcode, cmessage)
  END IF
END DO

! Firstly determine geometric mean diameter, number concentration and mass
! concentration for UKCA-mode sea-salt in accum and coarse soluble modes
! Initialise arrays to zero
seasalt_mconc(:,:,:,:)=0.0
seasalt_nconc(:,:,:,:)=0.0
seasalt_diam(:,:,:,:)=0.0

sixovrpix(:) = 6.0/(pi*x(:))
dens2num=zboltz/rgas
DO imode = mode_acc_sol, mode_cor_sol
  tot_mode_vmr(:,:,:) = 0.0
  idx=imode-mode_acc_sol+1
  meanfac=EXP(0.5*LOG(sigmag(imode))*LOG(sigmag(imode)))
  ! Check sea-salt exists for mode
  IF (.NOT. component(imode,cp_cl)) THEN
    WRITE(cmessage,'(A44,I2,A10,I2,A14)')                                      &
                        'l_ukca_coarse_no3_prod = TRUE but component ',        &
                         cp_cl,' and mode ',imode,' not turned on'
    errcode=2
    CALL ereport(RoutineName, errcode, cmessage)
  END IF
  ! 1. Extract modal sea-salt MMR and calculate total volume mixing ratio
  DO icp = 1,ncp
    IF (component(imode,icp)) THEN
      SELECT CASE (icp)
      CASE (cp_cl)
        seasalt_mconc(:,:,:,idx) = mode_tracers(:,:,:,midx(imode,icp))
        tot_mode_vmr(:,:,:) = tot_mode_vmr(:,:,:) +                            &
                              mode_tracers(:,:,:,midx(imode,icp))/rhocomp(icp)
      CASE DEFAULT
        tot_mode_vmr(:,:,:) = tot_mode_vmr(:,:,:) +                            &
                              mode_tracers(:,:,:,midx(imode,icp))/rhocomp(icp)
      END SELECT
    END IF !component
  END DO !icp
  ! 2. Determine geometric mean diameter for mode and sea-salt number
  !    concentration assuming sea-salt is externally mixed
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        IF (mode_tracers(i,j,k,nidx(imode)) > 0.0) THEN
          tmpdiam = (tot_mode_vmr(i,j,k) / mode_tracers(i,j,k,nidx(imode)) *   &
                    sixovrpix(imode) * dens2num)**athird
          tmpdiam = MAX(ddplim0(imode),tmpdiam)
          tmpdiam = MIN(tmpdiam,ddplim1(imode))
          seasalt_diam(i,j,k,idx) = tmpdiam*meanfac
          seasalt_nconc(i,j,k,idx) = seasalt_mconc(i,j,k,idx) /                &
                                     rhocomp(cp_cl) * sixovrpix(imode) *       &
                                     air_density(i,j,k) / (tmpdiam**3.0)
        ELSE
          seasalt_diam(i,j,k,idx) = ddplim0(imode)*meanfac
          seasalt_nconc(i,j,k,idx) = 0.0
        END IF
      END DO
    END DO
  END DO
END DO !imode

SELECT CASE (glomap_config%i_dust_scheme)
CASE (1)
  CALL aer_no3_6bindu (row_length, rows, model_levels, timestep,               &
                       mode_cor_sol-mode_acc_sol+1,mode_cor_sol-mode_acc_sol+1,&
                       dust_div2, dust_div3,dust_div4, dust_div5, dust_div6,   &
                       seasalt_mconc, seasalt_nconc, seasalt_diam,             &
                       hono2_mmr, nano3_mmr, tropopause_level,                 &
                       rel_humid_frac_clr, t_theta_levels, air_density,        &
                       dmas_nacl(:,:,:,mode_acc_sol:mode_cor_sol),             &
                       dmas_hono2, dmas_no3(:,:,:,mode_acc_sol:mode_cor_sol))

  ! Convert tendencies from kg/kg/s to kg/m2/s
  DO imode = mode_acc_sol, mode_cor_sol
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          ! NO3
          ems_no3 = dmas_no3(i,j,k,imode) * air_burden(i,j,k)
          IF (ems_no3 > loc_ems_eps) THEN
            dmas_no3(i,j,k,imode)  =  ems_no3
          ELSE
            dmas_no3(i,j,k,imode)  =  0.0
          END IF
          ! NaCl
          ems_nacl = dmas_nacl(i,j,k,imode) * air_burden(i,j,k)
          IF (ems_nacl < -1.0*loc_ems_eps) THEN
            dmas_nacl(i,j,k,imode) =  ems_nacl
          ELSE
            dmas_nacl(i,j,k,imode) =  0.0
          END IF
        END DO
      END DO
    END DO
  END DO

  ! Update HONO2 array
  hono2_mmr(:,:,:) = hono2_mmr(:,:,:) + dmas_hono2(:,:,:) * timestep

CASE (2)
  CALL aer_no3_2bindu (row_length, rows, model_levels, timestep,               &
                       mode_cor_sol-mode_acc_sol+1,mode_cor_sol-mode_acc_sol+1,&
                       dust_div1, dust_div2,                                   &
                       seasalt_mconc, seasalt_nconc, seasalt_diam,             &
                       hono2_mmr, nano3_mmr, tropopause_level,                 &
                       rel_humid_frac_clr, t_theta_levels, air_density,        &
                       dmas_nacl(:,:,:,mode_acc_sol:mode_cor_sol),             &
                       dmas_hono2, dmas_no3(:,:,:,mode_acc_sol:mode_cor_sol))

  ! Convert tendencies from kg/kg/s to kg/m2/s
  DO imode = mode_acc_sol, mode_cor_sol
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          ! NO3
          ems_no3 = dmas_no3(i,j,k,imode) * air_burden(i,j,k)
          IF (ems_no3 > loc_ems_eps) THEN
            dmas_no3(i,j,k,imode)  =  ems_no3
          ELSE
            dmas_no3(i,j,k,imode)  =  0.0
          END IF
          ! NaCl
          ems_nacl = dmas_nacl(i,j,k,imode) * air_burden(i,j,k)
          IF (ems_nacl < -1.0*loc_ems_eps) THEN
            dmas_nacl(i,j,k,imode) =  ems_nacl
          ELSE
            dmas_nacl(i,j,k,imode) =  0.0
          END IF
        END DO
      END DO
    END DO
  END DO

  ! Update HONO2 array
  hono2_mmr(:,:,:) = hono2_mmr(:,:,:) + dmas_hono2(:,:,:) * timestep

CASE (3)

  ! As for UKCA seasalt, calculate dust diameter, mass and number concentration
  ALLOCATE (dust_diam(1:row_length,1:rows,1:model_levels,1:2) )
  ALLOCATE (dust_mconc(1:row_length,1:rows,1:model_levels,1:2))
  ALLOCATE (dust_nconc(1:row_length,1:rows,1:model_levels,1:2))
  dust_mconc(:,:,:,:)=0.0
  dust_nconc(:,:,:,:)=0.0
  dust_diam(:,:,:,:)=0.0

  DO imode = mode_acc_insol, mode_cor_insol
    tot_mode_vmr(:,:,:) = 0.0
    idx=imode-mode_acc_insol+1
    meanfac=EXP(0.5*LOG(sigmag(imode))*LOG(sigmag(imode)))

    ! Check dust exists for mode
    IF (.NOT. component(imode,cp_du)) THEN
      WRITE(cmessage,'(A44,I2,A10,I2,A14)')                                    &
                          'l_ukca_coarse_no3_prod = TRUE but component ',      &
                           cp_du,' and mode ',imode,' not turned on'
      errcode=2
      CALL ereport(RoutineName, errcode, cmessage)
    END IF

    ! 1. Extract modal dust MMR and calculate total volume mixing ratio
    DO icp = 1,ncp
      IF (component(imode,icp)) THEN
        SELECT CASE (icp)
        CASE (cp_du)
          dust_mconc(:,:,:,idx) = mode_tracers(:,:,:,midx(imode,icp))
          tot_mode_vmr(:,:,:) = tot_mode_vmr(:,:,:) +                          &
                                (mode_tracers(:,:,:,midx(imode,icp)) /         &
                                 rhocomp(icp))
        CASE DEFAULT
          tot_mode_vmr(:,:,:) = tot_mode_vmr(:,:,:) +                          &
                                (mode_tracers(:,:,:,midx(imode,icp)) /         &
                                 rhocomp(icp))
        END SELECT
      END IF !component
    END DO !icp

    ! 2. Determine geometric mean diameter for mode and dust number
    !    concentration assuming dust is externally mixed
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          IF (mode_tracers(i,j,k,nidx(imode)) > 0.0) THEN
            tmpdiam = (tot_mode_vmr(i,j,k) / mode_tracers(i,j,k,nidx(imode)) * &
                      sixovrpix(imode) * dens2num)**athird
            tmpdiam = MAX(ddplim0(imode),tmpdiam)
            tmpdiam = MIN(tmpdiam,ddplim1(imode))
            dust_diam(i,j,k,idx) = tmpdiam*meanfac
            dust_nconc(i,j,k,idx) = dust_mconc(i,j,k,idx) /                    &
                                    rhocomp(cp_du) * sixovrpix(imode) *        &
                                    air_density(i,j,k) / (tmpdiam**3.0)
          ELSE
            dust_diam(i,j,k,idx) = ddplim0(imode)*meanfac
            dust_nconc(i,j,k,idx) = 0.0
          END IF
        END DO
      END DO
    END DO
  END DO !imode

  CALL aer_no3_ukcadu (row_length, rows, model_levels, timestep,               &
                       mode_cor_sol-mode_acc_sol+1,mode_cor_sol-               &
                       mode_acc_sol+1, mode_cor_insol-mode_acc_insol+1,        &
                       dust_mconc, dust_nconc, dust_diam,                      &
                       seasalt_mconc, seasalt_nconc, seasalt_diam,             &
                       hono2_mmr, tropopause_level,                            &
                       rel_humid_frac_clr, t_theta_levels, air_density,        &
                       dmas_dust(:,:,:,mode_acc_insol:mode_cor_insol),         &
                       dmas_nacl(:,:,:,mode_acc_sol:mode_cor_sol),             &
                       dmas_hono2, dmas_no3(:,:,:,mode_acc_sol:mode_cor_sol))

  ! Convert tendencies from kg/kg/s to kg/m2/s
  DO imode = mode_acc_sol, mode_cor_sol
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          ! NO3
          ems_no3 = dmas_no3(i,j,k,imode) * air_burden(i,j,k)
          IF (ems_no3 > loc_ems_eps) THEN
            dmas_no3(i,j,k,imode)  =  ems_no3
          ELSE
            dmas_no3(i,j,k,imode)  =  0.0
          END IF
          ! NaCl
          ems_nacl = dmas_nacl(i,j,k,imode) * air_burden(i,j,k)
          IF (ems_nacl < -1.0*loc_ems_eps) THEN
            dmas_nacl(i,j,k,imode) =  ems_nacl
          ELSE
            dmas_nacl(i,j,k,imode) =  0.0
          END IF
        END DO
      END DO
    END DO
  END DO
  DO imode = mode_acc_insol,mode_cor_insol
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          ! Dust
          ems_dust = dmas_dust(i,j,k,imode) * air_burden(i,j,k)
          IF (ems_dust < -1.0*loc_ems_eps) THEN
            dmas_dust(i,j,k,imode) =  ems_dust
          ELSE
            dmas_dust(i,j,k,imode) =  0.0
          END IF
        END DO
      END DO
    END DO
  END DO

  ! Update HONO2 array
  hono2_mmr(:,:,:) = hono2_mmr(:,:,:) + dmas_hono2(:,:,:) * timestep

CASE DEFAULT
  errcode=-1
  cmessage = 'No dust scheme defined - needed for secondary NO3 production'
  CALL ereport (ModuleName//':'//RoutineName, errcode, cmessage)
END SELECT

! Check for negative HONO2 MMRs
WHERE (hono2_mmr(:,:,:) < 0.0)
  hono2_mmr(:,:,:) = 0.0
END WHERE

! Deallocate dust arrays if allocated
IF (ALLOCATED(dust_diam))  DEALLOCATE (dust_diam)
IF (ALLOCATED(dust_mconc)) DEALLOCATE (dust_mconc)
IF (ALLOCATED(dust_nconc)) DEALLOCATE (dust_nconc)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE ukca_prod_no3_coarse

!-----------------------------------------------------------------------

SUBROUTINE ukca_no3_check_values(row_length, rows, model_levels,               &
                                 mode_tracers)

USE ukca_mode_tracer_maps_mod, ONLY: nmr_index, mmr_index
USE ukca_config_specification_mod, ONLY: glomap_variables
USE ukca_mode_setup,           ONLY: nmodes
USE ukca_config_defs_mod,      ONLY: n_mode_tracers
USE asad_mod,                  ONLY: jpctr

USE parkind1, ONLY: jpim, jprb      ! DrHook
USE yomhook, ONLY: lhook, dr_hook   ! DrHook

IMPLICIT NONE

INTEGER, INTENT(IN)  :: row_length, rows, model_levels
REAL, INTENT(IN OUT) :: mode_tracers(1:row_length,1:rows,1:model_levels,       &
                                     1:n_mode_tracers)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
LOGICAL, POINTER :: mode (:)
INTEGER, POINTER :: ncp

INTEGER :: ifirst, nitra, mitra, imode, icp

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_NO3_CHECK_VALUES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp

ifirst = jpctr + 1
DO imode = 1,nmodes
  IF (mode(imode)) THEN

    nitra = nmr_index(imode) - ifirst + 1
    WHERE (mode_tracers(:,:,:,nitra) < 0.0)
      mode_tracers(:,:,:,nitra) = 0.0
    END WHERE

    DO icp = 1,ncp
      IF (component(imode,icp)) THEN
        mitra = mmr_index(imode,icp) - ifirst + 1
        WHERE (mode_tracers(:,:,:,mitra) < 0.0)
          mode_tracers(:,:,:,mitra) = 0.0
        END WHERE
      END IF
    END DO

  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE ukca_no3_check_values

END MODULE ukca_prod_no3_mod
