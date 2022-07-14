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
!    Defines the aerosol dry deposition coefficients.
!    Based on the parameterisation of Zhang et al (2001) which
!    uses the method in the model of Slinn (1982).
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
MODULE ukca_ddepaer_coeff_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_DDEPAER_COEFF_MOD'

INTEGER, PARAMETER :: n_zhg = 17
! Number of Zhang et al recognised surface type categories

! Land surface types from Table:2 of Zhang et al
! ==============================================
!
! 1. Evergreen needle leaf trees
! 2. Evergreen broad leaf trees
! 3. Deciduous needle leaf trees
! 4. Deciduous broadleaf trees
! 5. Mixed  broad leaf  and  needle leaf trees
! 6. Grass
! 7. Crops, mixed farming
! 8. Desert
! 9. Tundra
! 10. Shrubs  and  interrupted  woodlands
! 11. Wet land with plants
! 12. Ice cap and glacier
! 13. Inland water
! 14. Ocean
! 15. Urban
!
! Derived categories:
!  If evergreen and deciduous not separated: (use average of respec. categories)
! 16. Needle leaf trees
! 17. Broad leaf trees
!
INTEGER, PARAMETER, PUBLIC :: zhg_eg_nedleaf  = 1
INTEGER, PARAMETER, PUBLIC :: zhg_eg_brdleaf  = 2
INTEGER, PARAMETER, PUBLIC :: zhg_dec_nedleaf = 3
INTEGER, PARAMETER, PUBLIC :: zhg_dec_brdleaf = 4
INTEGER, PARAMETER, PUBLIC :: zhg_mix_brdned_leaf = 5
INTEGER, PARAMETER, PUBLIC :: zhg_grass       = 6
INTEGER, PARAMETER, PUBLIC :: zhg_crop        = 7
INTEGER, PARAMETER, PUBLIC :: zhg_desert      = 8
INTEGER, PARAMETER, PUBLIC :: zhg_tundra      = 9
INTEGER, PARAMETER, PUBLIC :: zhg_shrub       = 10
INTEGER, PARAMETER, PUBLIC :: zhg_wetl_veg    = 11
INTEGER, PARAMETER, PUBLIC :: zhg_ice         = 12
INTEGER, PARAMETER, PUBLIC :: zhg_inl_water   = 13
INTEGER, PARAMETER, PUBLIC :: zhg_ocean       = 14
INTEGER, PARAMETER, PUBLIC :: zhg_urban       = 15
INTEGER, PARAMETER, PUBLIC :: zhg_ned_leaf    = 16
INTEGER, PARAMETER, PUBLIC :: zhg_brd_leaf    = 17

! Indices for 'smooth' surface types - i.e. soil, water, ice
! Values differ for 'old' and 'new' method
INTEGER, PUBLIC :: ls_soil, ls_water, ls_ocean, ls_ice

! Coefficients used in Zhang et al scheme
! CR    : Characteristic radius of collectors (m) -conv from A(mm) in the paper
! Y     : Parameter for calculating Brownian diffusion
! ALPHA : Parameter for calculating EIM

REAL, ALLOCATABLE, PUBLIC :: alpha(:)
REAL, ALLOCATABLE, PUBLIC :: cr(:)
REAL, ALLOCATABLE, PUBLIC :: yr(:)

CONTAINS

! Subroutine Interface:
SUBROUTINE set_ddepaer_coeff(l_improve_aero_drydep)
!
! Sets up the Arrays for Land-type based Coefficients used in the
! parametrisation: Table-3, Zhang et al
! https://doi.org/10.1016/S1352-2310(00)00326-5
!
! Currently contains two methods for specification:
! - 'Old Method' that expands the arrays based on number of model surface types
! - 'New' method that is independent of the number of surface types
!

USE umPrintMgr,              ONLY: umPrint, umMessage
USE ereport_mod,             ONLY: ereport
USE errormessagelength_mod,  ONLY: errormessagelength
USE ukca_config_specification_mod, ONLY: ukca_config

USE parkind1,                ONLY: jpim, jprb
USE yomhook,                 ONLY: lhook, dr_hook

IMPLICIT NONE

! Switch to use Old or New methods
LOGICAL, INTENT(IN) :: l_improve_aero_drydep

!    Local Variables
!
!!!! New method !!!!
!
! Coefficient arrays defined and indexed only as per Zhang et al categories.
! irrespective of number of surface types active in this configuration.
! Assumes that each grid box has a single 'dominant' surface type and
! the coefficients then accessed as e.g. CR(lscat(gridbox)).
!
! NOTE: Values used represent 3 out of the 5 'Seasonal categories' described.
! i.e consistent with:
! 'Midsummer with lush vegetation' or 'Autumn with cropland not harvested'; or
! 'Transitional spring, part green short annuals'
! but neglect:
! 'Late autumn aft frost, no snow' and 'Winter, snow on ground and sub-freezing'
!

REAL, PARAMETER :: CR_fix(n_zhg)=     [0.002, 0.005, 0.002, 0.005, 0.005,      &
                                       0.002, 0.002, 0.000, 0.000,  0.01,      &
                                       0.01,  0.000, 0.000, 0.000,  0.01,      &
                                       0.0035, 0.0035]
!                                      16 = AVG(1,2), 17= AVG(3,4)
REAL, PARAMETER :: ALPHA_fix(n_zhg) = [1.0,  0.6,   1.1,   0.8, 0.8,           &
                                       1.2,  1.2,  50.0,  50.0, 1.3,           &
                                       2.0, 50.0, 100.0, 100.0, 1.5,           &
                                       0.8, 0.85]
REAL, PARAMETER :: YR_fix(n_zhg) =    [0.56, 0.58, 0.56, 0.56, 0.56,           &
                                       0.54, 0.54, 0.54, 0.54, 0.54,           &
                                       0.54, 0.54, 0.50, 0.50, 0.56,           &
                                       0.57, 0.56]

INTEGER, PARAMETER :: ntype_ddepaer = 9
! Number of recognised surface type categories (in old method)
! (see details in comment block below)

INTEGER :: ntype
INTEGER :: errcode                ! Variable passed to ereport

CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_DDEPAER_COEFF'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! New method
IF (l_improve_aero_drydep) THEN

  ALLOCATE(alpha( SIZE(alpha_fix) ))
  ALLOCATE(cr( SIZE(cr_fix) ))
  ALLOCATE(yr( SIZE(yr_fix) ))

  alpha(:) = alpha_fix(:)
  cr(:)    = cr_fix(:)
  yr(:)    = yr_fix(:)

  ! Specify indices for 'smooth' surface types as defined in the parent
  ls_soil  = zhg_desert
  ls_ice   = zhg_ice
  ls_water = zhg_inl_water
  ls_ocean = zhg_ocean

ELSE

  !!! 'Old Method''
  ! Specification for 9 surface tiles (5 veg + 4 nnveg)
  ! This is extended for 13,17,27 surface tile configurations by adding or
  ! repeating values for new types/ sub-types
  !
  !! here changed land-use categories 1/2 (trees) to match values used
  !! for forest in 5-LS-category representation (should match then).
  !
  ! Now set to match 9 UM Land-use types (ILSCAT)      YR  ALPHA  CR
  ! 1=BL tree  (Zhang cat=avg of 2,4) [Evrgrn,Dec BL] 0.56  1.00 0.005
  ! 2=NL tree  (Zhang cat=avg of 1,3) [Evrgrn,Dec NL] 0.56  1.00 0.005
  ! 3=C3 grass (Zhang cat=6)          [grass        ] 0.54  1.20 0.002
  ! 4=C4 grass (Zhang cat=6)          [grass        ] 0.54  1.20 0.002
  ! 5=Shrub    (Zhang cat=10)         [shrub i. wood] 0.54  1.30 0.010
  ! 6=Urban    (Zhang cat=15)         [urban        ] 0.56  1.50 0.010
  ! 7=Water    (Zhang cat=13/14)      [inl wat/ocean] 0.50 100.0 0.000
  ! 8=Soil     (Zhang cat=8)          [desert       ] 0.54 50.00 0.000
  ! 9=Ice      (Zhang cat=12)         [ice cap/glac.] 0.54 50.00 0.000
  !
  ! The category (water,forest,grass,desert) is determined earlier
  ! based on roughness length znot (desert not used at present)
  !--------------------------------------------------------------------

  IF ((ukca_config%ntype == 27) .OR.                                           &
    (ukca_config%ntype == 17) .OR.                                             &
    (ukca_config%ntype == 13)) THEN
    ntype = ukca_config%ntype
  ELSE
    ntype = ntype_ddepaer
  END IF

  ALLOCATE(alpha(ntype))
  ALLOCATE(cr(ntype))
  ALLOCATE(yr(ntype))

  SELECT CASE (ntype)
  CASE (9)
    yr    = [0.56,0.56,0.54,0.54,0.54,0.56,0.50,0.54,0.54]
    cr    = [5.0e-3,5.0e-3,2.0e-3,2.0e-3,1.0e-2,1.0e-2,0.0e0,0.0e0,0.0e0]
    alpha = [1.00,1.00,1.20,1.02,1.30,1.50,100.0,50.0,50.0]
  CASE (13, 17, 27)
    yr(1:6)    = [   0.56,   0.58,   0.58,   0.56,   0.56,   0.54]
    cr(1:6)    = [ 7.0e-3, 5.0e-3, 5.0e-3, 3.2e-3, 2.0e-3, 3.2e-3]
    alpha(1:6) = [   0.80,   0.60,   0.60,   1.10,   1.00,   1.20]
  CASE DEFAULT
    WRITE(umMessage,'(A)') 'NTYPE must equal 9 or 13 or 17 or 27'
    CALL umPrint(umMessage,src=RoutineName)
    cmessage = 'Unexpected value of NTYPE'
    errcode = 1000 + ntype
    CALL ereport(RoutineName,errcode,cmessage)
  END SELECT

  SELECT CASE (ntype)
  CASE (13)
    yr(7:13)    = [   0.54,   0.54,   0.54,   0.56,   0.50,   0.54,  0.54]
    cr(7:13)    = [ 3.2e-3, 1.0e-2, 1.0e-2, 1.0e-2,  0.0e0,  0.0e0, 0.0e0]
    alpha(7:13) = [   1.20,   1.30,   1.30,   1.50,  100.0,   50.0,  50.0]
  CASE (17, 27)
    yr(7:17)    = [   0.54,   0.54,   0.54,   0.54,   0.54,   0.54,            &
                      0.54,   0.56,   0.50,   0.54,   0.54]
    cr(7:17)    = [ 3.2e-3, 3.2e-3, 3.2e-3, 3.2e-3, 3.2e-3, 1.0e-2,            &
                    1.0e-2, 1.0e-2,  0.0e0,  0.0e0,  0.0e0]
    alpha(7:17) = [   1.20,   1.20,   1.20,   1.20,   1.20,   1.30,            &
                      1.30,   1.50,  100.0,   50.0,   50.0]
  END SELECT

  IF (ntype == 27 ) THEN
    yr(18:27)    = [  0.54,  0.54,  0.54,  0.54,  0.54,                        &
                      0.54,  0.54,  0.54,  0.54,  0.54 ]
    cr(18:27)    = [ 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0,                        &
                     0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0 ]
    alpha(18:27) = [  50.0,  50.0,  50.0,  50.0,  50.0,                        &
                      50.0,  50.0,  50.0,  50.0,  50.0 ]
  END IF

  ! Indices for smooth categories (See table above)
  ls_water = 7
  ls_ocean = 7
  ls_soil  = 8
  ls_ice   = 9

END IF    ! l_improve_aero_drydep

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_ddepaer_coeff

END MODULE ukca_ddepaer_coeff_mod
