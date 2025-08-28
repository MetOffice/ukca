! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Main routine for calculating online photolysis rates
!   using fast-jx. Borrows a lot from the fast-j main routine.
!   Further developments that are still required include:
!    (i) use of MODE aerosol optical depths
!    (ii) use of more complete blocking for potential speedup
!    (iii) update to fast-jx 6.6 which will have some small
!    improvements (eg will remove the requirement that jaceto is last)
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE ukca_fastjx_mod

IMPLICIT NONE

! Default private
PRIVATE

! Subroutines available outside this module
PUBLIC :: ukca_fastjx

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_FASTJX_MOD'

CONTAINS

SUBROUTINE ukca_fastjx(error_code_ptr,                                         &
  row_length, rows, model_levels, jppj,                                        &
  p_layer_boundaries,                                                          &
  t_theta_levels,                                                              &
  r_theta_levels, r_rho_levels,                                                &
  longitude, sin_latitude,                                                     &
  z_top_of_model,                                                              &
  sulph_aitk, sulph_accu,                                                      &
  qcl, qcf, area_cloud_fraction,                                               &
  conv_cloud_lwp, conv_cloud_top, conv_cloud_base,                             &
  conv_cloud_amount, aod_sulph_aitk, aod_sulph_accum,                          &
  surf_albedo,                                                                 &
  ozone,                                                                       &
  land_fraction,                                                               &
  current_time,                                                                &
  photol_rates_fastjx, error_message, error_routine)

USE fastjx_data,    ONLY: fastjx_set_limits,                                   &
                          fastjx_allocate_memory,                              &
                          fastjx_deallocate_memory, nsl,                       &
                          Blocking_Mode,kpcx,                                  &
                          tau, daynumber,                                      &
                          sza, szafac, sza_2d, szafac_2d, u0,                  &
                          sa_block, fl_cyc,                                    &
                          rz_3d, rz_all, zzht,                                 &
                          pz_3d, pz_all, tz_3d, sa_2d,                         &
                          dm_3d, o3_3d,                                        &
                          ods_3d, odw_3d, odi_3d,                              &
                          solcyc_spec, solcyc_av, solcyc_quanta, solcyc_ts

USE photol_config_specification_mod, ONLY: photol_config
USE photol_constants_mod,   ONLY: c_o3 => const_o3_mmr_vmr,                    &
                                  m_s => const_molemass_sulp,                  &
                                  m_nh42so4 => const_molemass_nh42so4,         &
                                  gg => const_g, tm => const_tm,               &
                                  pi_over_180 => const_pi_over_180,            &
                                  avogadro => const_avogadro,                  &
                                  m_air => const_molemass_air
USE ukca_error_mod,          ONLY: maxlen_message, maxlen_procname,            &
                                   error_report, errcode_value_unknown
USE umPrintMgr,              ONLY: umMessage, umPrint, PrintStatus,            &
                                   PrStatus_Diag

USE yomhook,                 ONLY: lhook, dr_hook
USE parkind1,                ONLY: jprb, jpim
USE fastjx_set_aer_mod,      ONLY: fastjx_set_aer
USE fastjx_solar2_mod,       ONLY: fastjx_solar2
USE fastjx_photoj_mod,       ONLY: fastjx_photoj
USE photol_solflux_mod,      ONLY: photol_solflux

IMPLICIT NONE

! Dimensions of UKCA domain
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: jppj ! number of photolytic reactions

! Pressure on rho levels
REAL, INTENT(IN)    :: p_layer_boundaries(row_length,rows,model_levels+1)
! Temperature
REAL, INTENT(IN)    :: t_theta_levels(row_length,rows,model_levels)
! - height of theta levels
REAL, INTENT(IN)  :: r_theta_levels(row_length,rows,0:model_levels)
! - height of rho levels
REAL, INTENT(IN)  :: r_rho_levels(row_length,rows,model_levels)
! longitude (degrees)
REAL, INTENT(IN)    :: longitude(row_length,rows)
! SIN(latitude)
REAL, INTENT(IN)    :: sin_latitude(row_length,rows)
! Sulphate aerosol (aitken mode)
REAL, INTENT(IN) :: sulph_aitk(row_length,rows,model_levels)
! Sulphate aerosol (accumulation mode)
REAL, INTENT(IN) :: sulph_accu(row_length,rows,model_levels)
! Top of model
REAL, INTENT(IN)    :: z_top_of_model
! liquid water cloud
REAL, INTENT(IN)    :: qcl(row_length,rows,model_levels)
! ice water cloud
REAL, INTENT(IN)    :: qcf(row_length,rows,model_levels)
! cloud area fraction
REAL, INTENT(IN)    :: area_cloud_fraction(row_length,rows,model_levels)
! convective cloud amount
REAL, INTENT(IN)    :: conv_cloud_amount(row_length,rows,model_levels)
! ozone mmr
REAL, INTENT(IN)    :: ozone(row_length,rows,model_levels)
! Convective cloud LWP
REAL, INTENT(IN)    :: conv_cloud_lwp(row_length,rows)
! Aerosol Optical Depths - sulphate
REAL, INTENT(IN) :: aod_sulph_aitk(row_length, rows, model_levels)
REAL, INTENT(IN) :: aod_sulph_accum(row_length, rows, model_levels)
! Surface albedo
REAL, INTENT(IN)    :: surf_albedo(row_length,rows)
! Land fraction
REAL, INTENT(IN)    :: land_fraction(row_length,rows)
! Convective cloud top
INTEGER, INTENT(IN) :: conv_cloud_top(row_length,rows)
! Convective cloud bottom
INTEGER, INTENT(IN) :: conv_cloud_base(row_length,rows)
! Current model time
INTEGER, INTENT(IN) :: current_time(7)
! Output photolysis rates
REAL, INTENT(IN OUT) :: photol_rates_fastjx(row_length,rows,model_levels,jppj)

! error handling arguments
INTEGER, POINTER, INTENT(IN) :: error_code_ptr
CHARACTER(LEN=maxlen_message), OPTIONAL, INTENT(OUT) :: error_message
                                                       ! Error return message
CHARACTER(LEN=maxlen_procname), OPTIONAL, INTENT(OUT) :: error_routine
                                         ! Routine in which error was trapped
! Local variables

! used for solar cycle to inform calulation method
LOGICAL, PARAMETER :: l_lookup = .FALSE.

! Loop variables
INTEGER                       :: i,j,l,ix
! Effective diameter
REAL                          :: d_eff
! Timestep in hours
REAL                          :: timej
! Latitude and Longitude of current point
INTEGER                       :: nslat
INTEGER                       :: nslon
! Local time variables
INTEGER                       :: i_day_number
INTEGER                       :: i_hour
INTEGER                       :: i_minute
INTEGER                       :: i_second

CHARACTER(LEN=maxlen_message) :: err_message

! First time?
LOGICAL                       :: first=.TRUE.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_FASTJX'


! *************************************
! End of Header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code_ptr = 0
err_message = ''
IF (PRESENT(error_message)) error_message = ''
IF (PRESENT(error_routine)) error_routine = ''

IF (PrintStatus >= PrStatus_Diag) THEN
  WRITE(umMessage,'(A)') 'UKCA_FASTJX inputs:'
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 'p_layer_boundaries: ',                        &
                   MINVAL(p_layer_boundaries), MAXVAL(p_layer_boundaries),     &
                   SUM(p_layer_boundaries) / SIZE(p_layer_boundaries)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 't_theta_levels: ',                            &
                   MINVAL(t_theta_levels), MAXVAL(t_theta_levels),             &
                   SUM(t_theta_levels) / SIZE(t_theta_levels)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 'SO4_aitken: ',                                &
                   MINVAL(sulph_aitk), MAXVAL(sulph_aitk),                     &
                   SUM(sulph_aitk) / SIZE(sulph_aitk)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 'SO4_accum: ',                                 &
                   MINVAL(sulph_accu), MAXVAL(sulph_accu),                     &
                   SUM(sulph_accu) / SIZE(sulph_accu)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,E12.3)') 'z_top_of_model: ', z_top_of_model
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 'QCL: ',                                       &
                   MINVAL(qcl), MAXVAL(qcl),                                   &
                   SUM(qcl) / SIZE(qcl)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 'QCF: ',                                       &
                   MINVAL(qcf), MAXVAL(qcf),                                   &
                   SUM(qcf) / SIZE(qcf)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 'Area_cloud_frac:',                            &
                   MINVAL(area_cloud_fraction), MAXVAL(area_cloud_fraction),   &
                   SUM(area_cloud_fraction) / SIZE(area_cloud_fraction)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 'Conv_cloud_amount:',                          &
                   MINVAL(conv_cloud_amount), MAXVAL(conv_cloud_amount),       &
                   SUM(conv_cloud_amount) / SIZE(conv_cloud_amount)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 'Ozone MMR:',                                  &
                   MINVAL(ozone), MAXVAL(ozone),                               &
                   SUM(ozone) / SIZE(ozone)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 'Conv_cloud_LWP:',                             &
                   MINVAL(conv_cloud_lwp), MAXVAL(conv_cloud_lwp),             &
                   SUM(conv_cloud_lwp) / SIZE(conv_cloud_lwp)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 'Surf_albedo:',                                &
                   MINVAL(surf_albedo), MAXVAL(surf_albedo),                   &
                   SUM(surf_albedo) / SIZE(surf_albedo)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3E12.3)') 'Land_fraction:',                              &
                   MINVAL(land_fraction), MAXVAL(land_fraction),               &
                   SUM(land_fraction) / SIZE(land_fraction)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3I6)') 'Conv_cloud_top:',                                &
                   MINVAL(conv_cloud_top), MAXVAL(conv_cloud_top),             &
                   SUM(conv_cloud_top) / SIZE(conv_cloud_top)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,'(A,3I6)') 'Conv_cloud_base:',                               &
                   MINVAL(conv_cloud_base), MAXVAL(conv_cloud_base),           &
                   SUM(conv_cloud_base) / SIZE(conv_cloud_base)
  CALL umPrint(umMessage,src='ukca_fastjx')
END IF

! Initialise
photol_rates_fastjx(:,:,:,:) = 0.0

! Set local time variables
i_day_number          = current_time(7)
i_hour                = current_time(4)
i_minute              = current_time(5)
i_second              = current_time(6)

! Set Blocking mode:          0) Column-by-column
        !                     1) blocking 1 row
        !                     2) blocking domain
        !                     3) compressed  (not implemented)
        !                     4) load balancing (not implemented)
Blocking_Mode = 2

! Allocate arrays etc.
CALL fastjx_set_limits(error_code_ptr, row_length, rows, model_levels, jppj,   &
        error_message=error_message, error_routine=error_routine)
IF ( error_code_ptr > 0 ) THEN

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
CALL fastjx_allocate_memory

! Initialise arrays for levels/units appropriate for fast-j
! Need to update to include aerosols
CALL fastjx_set_arrays
! Set aerosol/ cloud indices
CALL fastjx_set_aer(error_code_ptr, error_message=error_message,               &
                    error_routine=error_routine)

IF ( error_code_ptr > 0 ) THEN
  CALL fastjx_deallocate_memory
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Set variables concerning model time
! Convert timestep into hours
timej              = photol_config%timestep/3600.0

  ! Day of the year
daynumber          = i_day_number

  ! Time in hours
tau                = i_hour*1.0+i_minute/60.0+i_second/3600.0                  &
                     - timej*0.5

! Calculate solar zenith angles

CALL fastjx_solar2 (timej, sin_latitude, pi_over_180*longitude,                &
                    photol_config%l_cal360)

! Modify solar flux if using solar cycle
IF (photol_config%i_solcylc_type > 0) THEN
  CALL photol_solflux(current_time, tau, l_lookup, solcyc_spec, fl_cyc)
ELSE
  fl_cyc = 0.0
END IF

! Block the data appropriately and call photolysis routines
! Still need to add modes 3 (compressed) and 4 (load balancing)

! Initialise longitude counter to 1
nslon = 1

SELECT CASE (Blocking_Mode)

  ! if blocking point by point
CASE (0)

  ! Loop over rows
  DO j=1,rows

    ! Loop over row, setting longitude counter
    DO i=1,row_length

      ! Set latitude counter to row number
      nslat  = j
      nslon  = i

      ! Loop over row, setting longitude counter
      nsl(1,1)=nslon
      nsl(2,1)=nslat

      ! block using consistent approach
      DO ix=1, kpcx
        sza(ix) = sza_2d(nsl(1,ix),nsl(2,ix))
        szafac(ix) = szafac_2d(nsl(1,ix),nsl(2,ix))

        u0(ix) = sza(ix)*pi_over_180
        u0(ix) = COS(u0(ix))
      END DO

      CALL fastjx_photoj (photol_rates_fastjx)

    END DO
  END DO

  ! if blocking row by row
CASE (1)

  ! Loop over rows
  DO j=1,rows

    ! Set latitude counter to row number
    nslat  = j

    ! Loop over row, setting longitude counter
    DO i=1,row_length
      nsl(1,i)=nslon+(i-1)
      nsl(2,i)=nslat
    END DO

    ! block using consistent approach
    DO ix=1, kpcx
      sza(ix) = sza_2d(nsl(1,ix),nsl(2,ix))
      szafac(ix) = szafac_2d(nsl(1,ix),nsl(2,ix))

      u0(ix) = sza(ix)*pi_over_180
      u0(ix) = COS(u0(ix))

      sa_block(ix) = sa_2d(nsl(1,ix),nsl(2,ix))
      sa_block(ix) = MIN(1.0, sa_block(ix))
      sa_block(ix) = MAX(0.0, sa_block(ix))

    END DO

    CALL fastjx_photoj (photol_rates_fastjx)

  END DO

  ! *********************************
  ! If blocking whole domain
CASE (2)

  ! initialise latitude counter to 1
  nslat = 1

  ! Loop over rows
  DO j=1,rows

    ! Loop over columns
    DO i=1,row_length

      ! Calculate positions of longitude and latitude in blocked arrays
      l=i+(j-1)*row_length

      nsl(1,l)=nslon+(i-1)
      nsl(2,l)=nslat+(j-1)
    END DO ! columns
  END DO ! rows

  ! block using consistent approach
  DO ix=1, kpcx
    sza(ix) = sza_2d(nsl(1,ix),nsl(2,ix))
    szafac(ix) = szafac_2d(nsl(1,ix),nsl(2,ix))

    u0(ix) = sza(ix)*pi_over_180
    u0(ix) = COS(u0(ix))

    sa_block(ix) = sa_2d(nsl(1,ix),nsl(2,ix))
    sa_block(ix) = MIN(1.0, sa_block(ix))
    sa_block(ix) = MAX(0.0, sa_block(ix))

  END DO

  CALL fastjx_photoj (photol_rates_fastjx)

  ! *********************************
  ! No DEFAULT Case needed since Blocking_Mode is set in this routine itself
  ! and already checked in fastjx_set_limits
END SELECT

  ! Tidy up at the end
CALL fastjx_deallocate_memory

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

CONTAINS

! ######################################################################
SUBROUTINE fastjx_set_arrays

IMPLICIT NONE

! Loop variables
INTEGER :: i,j,k

! Total mass in column
REAL :: total_mass(1:row_length,1:rows)

! Mass column per layer
REAL :: d_mass(1:row_length,1:rows,1:model_levels)

! Local copies of SO4 aitken and accum to enable unit conversion
REAL :: sulph_aitk_local(1:row_length,1:rows,1:model_levels)
REAL :: sulph_accu_local(1:row_length,1:rows,1:model_levels)

! Sulphate total, in accumulation and aitken modes
REAL :: sulphur(1:row_length,1:rows,1:model_levels)

  ! Cloud optical depths
REAL :: odi(1:row_length,1:rows,1:model_levels)
REAL :: odw(1:row_length,1:rows,1:model_levels)
REAL :: ods(1:row_length,1:rows,1:model_levels)

  ! Conversion factor from kg to molecules
REAL :: masfac

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_SET_ARRAYS'


  ! *************************************************
  ! EOH
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  ! rz in cm
rz_3d(:,:,1)              =                                                    &
    r_Theta_levels(1:row_length,1:rows,0)*100.0
rz_3d(:,:,2:model_levels) =                                                    &
    r_rho_levels  (1:row_length,1:rows,2:model_levels)*100.0
rz_3d(:,:,model_levels+1) =                                                    &
    r_Theta_levels(1:row_length,1:rows,model_levels)*100.0

  ! calculate pressure at box edges
pz_3d(:,:,1)              = p_layer_boundaries(:,:,1)
pz_3d(:,:,2:model_levels) = p_layer_boundaries(:,:,2:model_levels)
pz_3d(:,:,model_levels+1) = 0.0

  ! Calculate mass in box from pressure differences
  ! using hydrostatic approximation
d_mass = (pz_3d(:,:,1:model_levels)                                            &
       - pz_3d(:,:,2:model_levels+1))/gg

  ! Calculate total mass within convective clouds
total_mass = 0.0
DO k = 1,model_levels
  WHERE ( k <= conv_cloud_top .AND. k >= conv_cloud_base)
    total_mass = total_mass + d_mass(:,:,k)
  END WHERE
END DO


  ! Initialise water vapour to large scale precipitation
odw(:,:,1:model_levels)   = qcl(:,:,1:model_levels)
odi(:,:,1:model_levels)   = qcf(:,:,1:model_levels)

!-----------------------------------------------------------------
! Only add the convective cloud if we are not using the PC2 cloud scheme
IF (.NOT. photol_config%l_cloud_pc2) THEN

  ! I L-3d_cca (cloud levels?)
  IF (photol_config%l_3d_cca) THEN
    DO k = 1, photol_config%n_cca_lev

      WHERE (conv_cloud_top > 0 .AND. t_theta_levels(:,:,k) > tm)
        ! If above freezing point evaluate liquid water
        odw(:,:,k) = odw(:,:,k) + (conv_cloud_lwp(:,:)*                        &
                   conv_cloud_amount(:,:,k))/total_mass(:,:)
      ELSE WHERE (conv_cloud_top > 0)
        ! Else evaluate ice water
        odi(:,:,k) = odi(:,:,k) + (conv_cloud_lwp(:,:)*                        &
                   conv_cloud_amount(:,:,k))/total_mass(:,:)
      END WHERE

    END DO !  k levels

  ELSE

    ! Else loop over model levels
    DO k=1,model_levels

      WHERE (conv_cloud_top >= k .AND. conv_cloud_base <= k .AND.              &
             t_theta_levels(:,:,k) > tm)
        ! If above freezing point evaluate liquid water
        odw(:,:,k) = odw(:,:,k)+(conv_cloud_lwp(:,:)*                          &
                   conv_cloud_amount(:,:,1))/total_mass(:,:)
      ELSE WHERE (conv_cloud_top >= k .AND. conv_cloud_base <= k)
        ! Else evaluate ice water
        odi(:,:,k) = odi(:,:,k)+(conv_cloud_lwp(:,:)*                          &
                   conv_cloud_amount(:,:,1))/total_mass(:,:)
      END WHERE

    END DO ! k levels
  END IF ! l_3d_cca

END IF ! NOT l_cld_pc2

! Convert mass mixing ratios to column densities.
odw = odw*d_mass
odi = odi*d_mass

!--------------------------------------------------------------
! set effective radii for water drops, different for land&sea
! set effective diameter for ice crystals
! Formulae from John Edwards to convert from column densities to
! optical depths

d_eff=100.0 ! in microns

DO k=1,model_levels

  WHERE (land_fraction > 0.5)
    odw(:,:,k) = odw(:,:,k)*(-8.86964+1.67373e3/6.0)
  ELSE WHERE
    odw(:,:,k) = odw(:,:,k)*(-8.86964+1.67373e3/12.0)
  END WHERE

END DO

odi = odi*(-2.189e-3+3.311e3/d_eff+3.611/d_eff**2)

! Use approach of Briegleb(1992), JGR 97, 7603.
! to account for random overlap of cloud layers
odw(:,:,1:model_levels)    =  odw(:,:,1:model_levels)                          &
             * (area_cloud_fraction(:,:,1:model_levels))**(1.5)
odi(:,:,1:model_levels)    =  odi(:,:,1:model_levels)                          &
             * (area_cloud_fraction(:,:,1:model_levels))**(1.5)

! Set optical depth in top layer to be the same as top+1
! could just set to 0. For top level shouldn't make any difference
odw_3d(:,:,1:model_levels)      = odw
odw_3d(:,:,(model_levels+1))    = odw_3d(:,:,(model_levels))
odi_3d(:,:,1:model_levels)      = odi
odi_3d(:,:,(model_levels+1))    = odi_3d(:,:,(model_levels))

! **********************************************************
! Apply aerosol od columns @600nm to the concentrations

! Firsti, multiply by molecular weight ratio to convert from mass mixing ratio
! of sulphur atoms to mass mixing ratio of ammonium sulphate. Then, use d_mass
! to convert to mass of ammonium sulphate.
sulph_accu_local = sulph_accu * d_mass * m_nh42so4/m_s
sulph_aitk_local = sulph_aitk * d_mass * m_nh42so4/m_s

DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      sulph_accu_local(i,j,k) = sulph_accu_local(i,j,k) * aod_sulph_accum(i,j,k)
      sulph_aitk_local(i,j,k) = sulph_aitk_local(i,j,k) * aod_sulph_aitk(i,j,k)
    END DO
  END DO
END DO

! Sum the aitkin and accumulation type optical depths
sulphur = sulph_aitk_local + sulph_accu_local

! Copy sulphate od to global array
ods_3d(:,:,1:model_levels)  = sulphur(:,:,1:model_levels)
ods_3d(:,:,(model_levels+1))= sulphur(:,:,model_levels)

! ********************************************************
! Set other variables here

! Set surface albedo
sa_2d                           = surf_albedo

! Calculate pressure at box edges (include box to TOA)
! NB convert Pa to hPa for fastj routines
pz_all(:,:,1:(model_levels+1))  = p_layer_boundaries(:,:,:)/100.0
pz_all(:,:,(model_levels+2))    = 0.0e0

! Calculate heights of box edges
rz_all(:,:,1:model_levels+1)    = rz_3d(:,:,1:model_levels+1)
rz_all(:,:,model_levels+2)      = rz_all(:,:,model_levels+1)                   &
                                  + zzht

! Set temperature to be temperature
! Use top temperature for top+1 level
tz_3d(:,:,1:model_levels)       = t_theta_levels(:,:,:)
tz_3d(:,:,(model_levels+1))     = t_theta_levels(:,:,model_levels)

! Copy air mass from local array
dm_3d(:,:,1:model_levels)       = d_mass(:,:,1:model_levels)
! Evaluate top box from difference with TOA p i.e. 0!
dm_3d(:,:,(model_levels+1))     = pz_3d(:,:,model_levels)/gg

! Convert air mass from kg/m2 to molecules/cm2
masfac                          = avogadro/(m_air*1.0e4)
dm_3d(:,:,:)                    = dm_3d(:,:,:)*masfac

! Convert ozone mmr to vmr and write to fj_ozone array
o3_3d(:,:,1:model_levels)       = ozone/c_o3
o3_3d(:,:,(model_levels+1))     = o3_3d(:,:,model_levels)

! Multiply by molecules per m2 to get it in units fastj wants
o3_3d                           = o3_3d*dm_3d

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_set_arrays

!#######################################################################

END SUBROUTINE ukca_fastjx

END MODULE ukca_fastjx_mod
