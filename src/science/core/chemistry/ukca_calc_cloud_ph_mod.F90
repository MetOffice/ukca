! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Module containing subroutine for calculating a value for cloud water pH
!  on all model grid points and levels.
!
! Method:
!
!  An offline relationship has been determined between annual mean cloud pH
!  observations from Pye et al., (2020) and annual mean SO2 observations
!  (in micro grammes per metre cubed) from EMEP, CASTNET and EANET over the
!  period 1980 to 2015. NB that the SO2 concentrations are converted to a
!  logarithmic value (log10[SO2]) first before being used in the equation.
!  The calculated relationship takes the form of a 2nd order polynomial
!  as follows:
!
!  y = a x**2 + b * x + c
!
!  where the default derived coefficients are:
!    a = 0.142
!    b = -0.931
!    c = 4.65
!
!  When using this relationship to calculate cloud ph at all grid points,
!  the so2 concentrations must first be converted from mass mixing
!  ratios to micro grammes per metre cubed and then to a logarithmic value.
!  Then SO2 concentrations at each time step can be used in the above formula
!  to calculate a cloud pH value on each lat/lon/level. This can then be used
!  in other routines to calculate the fraction of gases dissolved in the
!  aqueous phase, aqueous phase reaction rates and wet deposition rates.
!
!  NB pH values calculated based on model SO2 concentrations are limited to be
!  between the minimum and maximum of observed cloud pH values which are 2.0
!  and 8.3
!
!  Reference: Pye et al., (2020) - The acidity of atmospheric particles and
!  clouds. Atmos. Chem. Phys. https://doi.org/10.5194/acp-20-4809-2020
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ######################################################################
!
! ----------------------------------------------------------------------

MODULE ukca_calc_cloud_ph_mod

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_CALC_CLOUD_PH_MOD'

PUBLIC :: ukca_calc_cloud_ph

CONTAINS

SUBROUTINE ukca_calc_cloud_ph(tracer, row_length, rows, model_levels,          &
                              pres, temp, H_plus_3d_arr)

! default global pH value of 5 from [H+] = 1.0e-5
USE ukca_constants, ONLY: H_plus, m_so2
USE ukca_config_constants_mod, ONLY: avogadro, boltzmann
USE ukca_cspecies,  ONLY: c_species,n_so2 ! to convert so2 concs
! import fitting coefficients from input Rose metadata
USE ukca_config_specification_mod, ONLY: ukca_config
USE parkind1,            ONLY: jprb, jpim
USE yomhook,             ONLY: lhook, dr_hook

IMPLICIT NONE

! variables from outside
INTEGER, INTENT(IN) :: row_length     ! No of points per row
INTEGER, INTENT(IN) :: rows           ! No of rows
INTEGER, INTENT(IN) :: model_levels   ! No of levels

! Pressure on theta levels
REAL, INTENT(IN)    :: pres(row_length, rows, model_levels)

! Temperature on theta levels
REAL, INTENT(IN)    :: temp(row_length, rows, model_levels)

! tracer MMR (row_length,rows,model_levels,ntracers)
REAL, INTENT(IN)    :: tracer(:,:,:,:)

! 3D array for calculated pH values at all grid points
REAL, INTENT(IN OUT) :: H_plus_3d_arr(:,:,:)

! local variables
INTEGER    :: i,j,k,l ! loop variables
REAL :: so2_conc ! Single value of current so2 concentrations
REAL :: so2_conc_micro ! so2 converted to ugpm3
REAL :: so2_conc_micro_log ! log10 so2 ugpm3
REAL :: ph_val ! single pH value calculated based on so2 concentrations
REAL :: cur_h_plus ! single value of [H+] convert from init pH
REAL :: new_h_plus ! single value of [H+] convert from pH
REAL, PARAMETER :: kgpcm3_to_ugpm3 = 1.0e15
                            ! Conversion from kg cm^-3 to ug m^-3
REAL :: aird(row_length, rows, model_levels) ! Number density of air (cm^-3)
! coefficients from 2nd order polynomial fit to observations of pH and SO2
! set based on values as input from Rose metadata
! default 1st fit parameter = 0.142
REAL :: fit_1 != 0.142
! default 2nd fit parameter = -0.931
REAL :: fit_2 != -0.931
! default intercept value = 4.65
REAL :: fit_3 != 4.65

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALC_CLOUD_PH'

! End of Header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Replace default fitting parameters with those input from rose metadata
fit_1=ukca_config%ph_fit_coeff_a
fit_2=ukca_config%ph_fit_coeff_b
fit_3=ukca_config%ph_fit_intercept

! calculate molecular concentration of air (/cm3) at all points
aird(:,:,:)=pres(:,:,:)/(temp(:,:,:)*boltzmann*1.0e6)

! now for each lat/lon/level calculate cloud water ph values based on current
! model SO2 concentrations using offline observationally defined relationship
DO k = 1, model_levels
  DO j=1,rows
    DO i=1,row_length
      ! l = (j-1) * row_length + i ! index for 2D array

      ! Extract so2 conc at lat/lon/level from tracer array and convert
      ! from mmr to vmr
      so2_conc = tracer(i,j,k,n_so2) / c_species(n_so2)

      ! Convert So2 from vmr into ug/m3
      so2_conc_micro =                                                         &
         so2_conc * aird(i,j,k) * m_so2 * 1.0e-3 * kgpcm3_to_ugpm3 / avogadro

      ! convert so2 concentrations in log10(so2)
      so2_conc_micro_log = LOG10(so2_conc_micro)

      ! Now calculate new value of pH based on model So2 concentrations
      ! using fit of observed relationship
      ! y = (fit_1 * x^2) + (fit_2 * x) + fit_3
      ph_val = (fit_1 * so2_conc_micro_log * so2_conc_micro_log) +             &
               (fit_2 * so2_conc_micro_log) + fit_3

      ! Convert pH to [H+] for use in rest of model (ph = -log10[H+])
      cur_h_plus = 10.0 ** (-1.0 * ph_val)

      ! Now test if pH value is within observed min (2.0) and max (8.3) range
      IF (ph_val < 2.0) THEN
        ph_val = 2.00
      END IF
      IF (ph_val > 8.3) THEN
        ph_val = 8.30
      END IF

      ! Convert pH to [H+] for use in rest of model (ph = -log10[H+])
      new_h_plus = 10.0 ** (-1.0 * ph_val)

      ! Put current pH values back into 2D array for use in ASAD
      H_plus_3d_arr(i,j,k) = new_h_plus

    END DO   ! i
  END DO    ! j
END DO     ! k

! Finished defined new values of cloud pH to use elsewhere
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_calc_cloud_ph

END MODULE ukca_calc_cloud_ph_mod
