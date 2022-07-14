! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Calculates aerodynamic and quasi-laminar resistances.
!  Returns Resa, Rb
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_CHEMISTRY_CTL.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
MODULE ukca_aerod_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_AEROD_MOD'

CONTAINS

SUBROUTINE ukca_aerod(row_length, rows, ntype, npft, nwater,                   &
  t0, p0, hflx, u_s, canht, gsf,                                               &
  zbl, z0tile, resa, rb, so4_vd)

USE asad_mod,                ONLY: ndepd, speci, nldepd, jpdd

USE ukca_config_specification_mod, ONLY: ukca_config

USE ukca_constants,          ONLY:                                             &
    m_alkaooh,                                                                 &
    m_aromooh,                                                                 &
    m_asoa,                                                                    &
    m_asvoc1,                                                                  &
    m_asvoc2,                                                                  &
    m_bsoa,                                                                    &
    m_bsvoc1,                                                                  &
    m_bsvoc2,                                                                  &
    m_buooh,                                                                   &
    m_dmso,                                                                    &
    m_etcho,                                                                   &
    m_h2o,                                                                     &
    m_h2so4,                                                                   &
    m_hbr,                                                                     &
    m_hobr,                                                                    &
    m_hcho,                                                                    &
    m_hcl,                                                                     &
    m_hocl,                                                                    &
    m_hono,                                                                    &
    m_isooh,                                                                   &
    m_isosoa,                                                                  &
    m_isosvoc1,                                                                &
    m_isosvoc2,                                                                &
    m_macr,                                                                    &
    m_macrooh,                                                                 &
    m_mecho,                                                                   &
    m_mecoch2ooh,                                                              &
    m_mekooh,                                                                  &
    m_meoh,                                                                    &
    m_gly,                                                                     &
    m_mgly,                                                                    &
    m_monoterp,                                                                &
    m_mpan,                                                                    &
    m_msa,                                                                     &
    m_mvkooh,                                                                  &
    m_nald,                                                                    &
    m_onitu,                                                                   &
    m_pan,                                                                     &
    m_ppan,                                                                    &
    m_sec_org,                                                                 &
    m_etoh,                                                                    &
    m_proh,                                                                    &
    m_hoc2h4ooh,                                                               &
    m_hoch2cho,                                                                &
    m_etco3h,                                                                  &
    m_hoch2co3h,                                                               &
    m_noa,                                                                     &
    m_hoc2h4no3,                                                               &
    m_phan,                                                                    &
    m_rtx24no3,                                                                &
    m_rn9no3,                                                                  &
    m_rn12no3,                                                                 &
    m_rn15no3,                                                                 &
    m_rn18no3,                                                                 &
    m_rtn28no3,                                                                &
    m_rtn25no3,                                                                &
    m_rtn23no3,                                                                &
    m_rtx28no3,                                                                &
    m_rtx22no3,                                                                &
    m_rn16ooh,                                                                 &
    m_rn19ooh,                                                                 &
    m_rn14ooh,                                                                 &
    m_rn17ooh,                                                                 &
    m_nru14ooh,                                                                &
    m_nru12ooh,                                                                &
    m_rn9ooh,                                                                  &
    m_rn12ooh,                                                                 &
    m_rn15ooh,                                                                 &
    m_rn18ooh,                                                                 &
    m_nrn6ooh,                                                                 &
    m_nrn9ooh,                                                                 &
    m_nrn12ooh,                                                                &
    m_ra13ooh,                                                                 &
    m_ra16ooh,                                                                 &
    m_ra19ooh,                                                                 &
    m_rtn28ooh,                                                                &
    m_nrtn28ooh,                                                               &
    m_rtn26ooh,                                                                &
    m_rtn25ooh,                                                                &
    m_rtn24ooh,                                                                &
    m_rtn23ooh,                                                                &
    m_rtn14ooh,                                                                &
    m_rtn10ooh,                                                                &
    m_rtx28ooh,                                                                &
    m_rtx24ooh,                                                                &
    m_rtx22ooh,                                                                &
    m_nrtx28ooh,                                                               &
    m_carb14,                                                                  &
    m_carb17,                                                                  &
    m_carb11a,                                                                 &
    m_carb10,                                                                  &
    m_carb13,                                                                  &
    m_carb16,                                                                  &
    m_carb9,                                                                   &
    m_carb12,                                                                  &
    m_carb15,                                                                  &
    m_ucarb12,                                                                 &
    m_carb15,                                                                  &
    m_ucarb12,                                                                 &
    m_nucarb12,                                                                &
    m_udcarb8,                                                                 &
    m_udcarb11,                                                                &
    m_udcarb14,                                                                &
    m_tncarb26,                                                                &
    m_tncarb10,                                                                &
    m_tncarb12,                                                                &
    m_tncarb11,                                                                &
    m_ccarb12,                                                                 &
    m_tncarb15,                                                                &
    m_rcooh25,                                                                 &
    m_txcarb24,                                                                &
    m_txcarb22,                                                                &
    m_ru12pan,                                                                 &
    m_rtn26pan,                                                                &
    m_rtn26pan,                                                                &
    m_aroh14,                                                                  &
    m_aroh17,                                                                  &
    m_arnoh14,                                                                 &
    m_arnoh17,                                                                 &
    m_anhy,                                                                    &
    m_cri,                                                                     &
    m_air,                                                                     &
    m_iepox,                                                                   &
    m_hmml,                                                                    &
    m_hucarb9,                                                                 &
    m_hpucarb12,                                                               &
    m_dhpcarb9,                                                                &
    m_dhpr12ooh,                                                               &
    m_dhcarb9,                                                                 &
    m_ru12no3,                                                                 &
    m_ru10no3

USE ukca_um_legacy_mod,      ONLY: rmol, vkman, gg => g
USE parkind1,                ONLY: jprb, jpim
USE yomhook,                 ONLY: lhook, dr_hook
USE ereport_mod,             ONLY: ereport
USE errormessagelength_mod,  ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: ntype
INTEGER, INTENT(IN) :: npft
INTEGER, INTENT(IN) :: nwater
REAL, INTENT(IN) :: t0(row_length,rows)
                             ! Surface temperature (K)
REAL, INTENT(IN) :: p0(row_length,rows)
                             ! Surface pressure (Pa)
REAL, INTENT(IN) :: hflx(row_length,rows)
                             ! Surface heat flux (W m-2)
REAL, INTENT(IN) :: u_s(row_length,rows)
                             ! Surface friction velocity (m s-1)
REAL, INTENT(IN) :: canht(row_length,rows,npft)
                             ! Canopy height of vegetation (m)
REAL, INTENT(IN) :: gsf(row_length,rows,ntype)
                             ! Surface tile fractions
REAL, INTENT(IN) :: zbl(row_length,rows)
                             ! Boundary layer depth
REAL, INTENT(IN OUT) :: z0tile(row_length,rows,ntype)
                             ! Roughness length on tiles (m)
REAL, INTENT(OUT) :: resa(row_length,rows,ntype)
                             ! Aerodynamic resistance (s m-1)
REAL, INTENT(OUT) :: rb(row_length,rows,jpdd)
                             ! Quasi-laminar resistance (s m-1)
REAL, INTENT(OUT) :: so4_vd(row_length,rows)
                             ! Sulphate aerosol dep velocity (m s-1)

LOGICAL, SAVE :: first = .TRUE.

INTEGER :: i, j, k, n     ! Loop counts

REAL :: cp = 1004.67         ! Heat capacity of air under const p (J K-1 kg-1)
REAL :: pr = 0.72            ! Prandtl number of air
REAL :: vair296 = 1830.0e-8  ! Dynamic viscosity of air at 296 K (kg m-1 s-1).
REAL :: twoth = 2.0 / 3.0    ! 2/3
REAL :: zref = 50.0          ! Reference height for dry deposition (m).
REAL :: d_h2o = 2.08e-5      ! Diffusion coefficent of water in air (m2 s-1)

REAL :: b0                   ! Temporary store
REAL :: b1                   ! Temporary store
REAL :: sc_pr                ! Ratio of Schmidt and Prandtl numbers

REAL, SAVE :: reus                ! Holds -cp/(vkman*g)
REAL, ALLOCATABLE, SAVE :: d0(:)  ! Diffusion coefficients

REAL :: l(row_length,rows)         ! Monin-Obukhov length (m)
REAL :: rho_air(row_length,rows)   ! Density of air at surface (kg m-3)
REAL :: kva(row_length,rows)       ! Kinematic velocity of air (m2 s-1)
REAL :: ustar(row_length,rows)     ! Friction velocity [checked] (m s-1)

REAL :: d(row_length,rows,npft)    ! Zero-plane displacement (m)

REAL :: z(row_length,rows,ntype)   ! Reference height (m)
REAL :: psi(row_length,rows,ntype) ! Businger function (dimensionless)

REAL, PARAMETER :: r_null = 1.0e50           ! Null resist to depos (1/r_null~0)

REAL :: bl_l                                  ! Boundary layer height divided by
                                              ! Monin-Obukhov length

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_AEROD'

! Ereport variables
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: icode

!     Assign diffusion coefficients, units m2 s-1. Set to -1
!     unless species dry deposits. If no value found in literature,
!     D0 calculated using: D(X) = D(H2O) * SQRT[RMM(H2O)/RMM(X)], where
!     X is the species in question and D(H2O) = 2.08 x 10^-5 m2 s-1
!     (Marrero & Mason, J Phys Chem Ref Dat, 1972).

!     The values of d0 will be used to flag those species
!     that dry deposit

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (first) THEN
  ALLOCATE(d0(jpdd))
  d0(:) = -1.0
  DO j = 1, ndepd
    SELECT CASE (speci(nldepd(j)))
    CASE ('O3        ','NO2       ','O3S       ','NO3       ')
      d0(j) = 1.4e-5
    CASE ('NO        ')
      IF (ukca_config%l_fix_improve_drydep) THEN
        ! Tang et al 2014 https://doi.org/10.5194/acp-14-9233-2014
        d0(j) = 2.3e-5
      ELSE
        ! NO = 6*NO2 (following Gian. 1998)
        d0(j) = 8.4e-5
      END IF
      ! RU14NO3 ~ ISON
    CASE ('HNO3      ','HONO2     ','ISON      ','B2ndry    ',                 &
          'A2ndry    ','N2O5      ','HO2NO2    ','HNO4      ',                 &
          'RU14NO3   ')
      d0(j) = 1.2e-5
    CASE ('H2O2      ','HOOH      ')
      d0(j) = 1.46e-5
    CASE ('CH3OOH    ','MeOOH     ','HCOOH     ')
      d0(j) = 1.27e-5
      ! CARB7 ~ HACET
    CASE ('C2H5OOH   ','EtOOH     ','MeCO3H    ','MeCO2H    ',                 &
          'HACET     ','CARB7     ')
      d0(j) = 1.12e-5
      ! RN10OOH ~ n-PrOOH
    CASE ('n_C3H7OOH ','i_C3H7OOH ','n-PrOOH   ','i-PrOOH   ',                 &
          'PropeOOH  ','RN10OOH   ')
      d0(j) = 1.01e-5
    CASE ('MeCOCH2OOH')
      d0(j) = d_h2o * SQRT(m_h2o / m_mecoch2ooh)
    CASE ('ISOOH     ','RU14OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_isooh)
    CASE ('HONO      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hono)
    CASE ('MACROOH   ','RU10OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_macrooh)
    CASE ('MEKOOH    ','RN11OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_mekooh)
    CASE ('ALKAOOH   ','RN8OOH    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_alkaooh)
    CASE ('AROMOOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_aromooh)
    CASE ('BSVOC1    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_bsvoc1)
    CASE ('BSVOC2    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_bsvoc2)
    CASE ('ASVOC1    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_asvoc1)
    CASE ('ASVOC2    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_asvoc2)
    CASE ('ISOSVOC1  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_isosvoc1)
    CASE ('ISOSVOC2  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_isosvoc2)
    CASE ('PAN       ')
      IF (ukca_config%l_fix_improve_drydep) THEN
        d0(j) = d_h2o * SQRT(m_h2o / m_pan)
      ELSE
        d0(j) = 0.31e-5
      END IF
    CASE ('PPAN      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_ppan)
    CASE ('MPAN      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_mpan)
    CASE ('ONITU     ')
      d0(j) = d_h2o * SQRT(m_h2o / m_onitu)
    CASE ('CO        ')
      d0(j) = 1.86e-5
    CASE ('CH4       ')
      IF (ukca_config%l_fix_improve_drydep) THEN
        ! Tang et al 2015 https://doi.org/10.5194/acp-15-5585-2015
        d0(j) = 2.2e-5
      ELSE
        d0(j) = 5.74e-5
      END IF
    CASE ('NH3       ')
      d0(j) = 2.08e-5
    CASE ('H2        ')
      d0(j) = 6.7e-5
    CASE ('SO2       ')
      d0(j) = 1.2e-5
    CASE ('DMSO      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_dmso)
    CASE ('Sec_Org   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_sec_org)
    CASE ('H2SO4     ')
      d0(j) = d_h2o * SQRT(m_h2o / m_h2so4)
    CASE ('MSA       ')
      d0(j) = d_h2o * SQRT(m_h2o / m_msa)
    CASE ('MVKOOH    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_mvkooh)
    CASE ('s-BuOOH   ','RN13OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_buooh)
    CASE ('ORGNIT    ')
      d0(j) = 1.0e-5
    CASE ('BSOA      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_bsoa)
    CASE ('ASOA      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_asoa)
    CASE ('ISOSOA    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_isosoa)
      ! Species added for sonsitency with 2D scheme
      ! Use simple approximation a la ExTC
    CASE ('HCHO      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hcho)
    CASE ('MeCHO     ')
      d0(j) = d_h2o * SQRT(m_h2o / m_mecho)
    CASE ('EtCHO     ')
      d0(j) = d_h2o * SQRT(m_h2o / m_etcho)
    CASE ('MACR      ','UCARB10   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_macr)
    CASE ('NALD      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_nald)
    CASE ('GLY       ','CARB3     ')
      d0(j) = d_h2o * SQRT(m_h2o / m_gly)
    CASE ('MGLY      ','CARB6     ')
      d0(j) = d_h2o * SQRT(m_h2o / m_mgly)
    CASE ('MeOH      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_meoh)
    CASE ('Monoterp  ', 'APINENE   ', 'BPINENE   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_monoterp)
    CASE ('HBr       ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hbr)
    CASE ('HOBr      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hobr)
    CASE ('HCl       ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hcl)
    CASE ('HOCl      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hocl)
      ! CRImech-only species:
    CASE ('EtOH      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_etoh)
    CASE ('i-PrOH    ', 'n-PrOH    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_proh)
    CASE ('HOC2H4OOH ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hoc2h4ooh)
    CASE ('HOCH2CHO  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hoch2cho)
    CASE ('EtCO3H    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_etco3h)
    CASE ('HOCH2CO3H ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hoch2co3h)
    CASE ('PHAN      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_phan)
    CASE ('NOA       ')
      d0(j) = d_h2o * SQRT(m_h2o / m_noa)
    CASE ('HOC2H4NO3 ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hoc2h4no3)
    CASE ('RTX24NO3  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtx24no3)
    CASE ('RN9NO3    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn9no3)
    CASE ('RN12NO3   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn12no3)
    CASE ('RN15NO3   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn15no3)
    CASE ('RN18NO3   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn18no3)
    CASE ('RTN28NO3  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtn28no3)
    CASE ('RTN25NO3  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtn25no3)
    CASE ('RTN23NO3  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtn23no3)
    CASE ('RTX28NO3  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtx28no3)
    CASE ('RTX22NO3  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtx22no3)
    CASE ('RN16OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn16ooh)
    CASE ('RN19OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn19ooh)
    CASE ('RN14OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn14ooh)
    CASE ('RN17OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn17ooh)
    CASE ('NRU14OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_nru14ooh)
    CASE ('NRU12OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_nru12ooh)
    CASE ('RN9OOH    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn9ooh)
    CASE ('RN12OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn12ooh)
    CASE ('RN15OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn15ooh)
    CASE ('RN18OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rn18ooh)
    CASE ('NRN6OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_nrn6ooh)
    CASE ('NRN9OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_nrn9ooh)
    CASE ('NRN12OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_nrn12ooh)
    CASE ('RA13OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_ra13ooh)
    CASE ('RA16OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_ra16ooh)
    CASE ('RA19OOH   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_ra19ooh)
    CASE ('RTN28OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtn28ooh)
    CASE ('NRTN28OOH ')
      d0(j) = d_h2o * SQRT(m_h2o / m_nrtn28ooh)
    CASE ('RTN26OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtn26ooh)
    CASE ('RTN25OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtn25ooh)
    CASE ('RTN24OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtn24ooh)
    CASE ('RTN23OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtn23ooh)
    CASE ('RTN14OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtn14ooh)
    CASE ('RTN10OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtn10ooh)
    CASE ('RTX28OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtx28ooh)
    CASE ('RTX24OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtx24ooh)
    CASE ('RTX22OOH  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtx22ooh)
    CASE ('NRTX28OOH ')
      d0(j) = d_h2o * SQRT(m_h2o / m_nrtx28ooh)
    CASE ('CARB14    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_carb14)
    CASE ('CARB17    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_carb17)
    CASE ('CARB11A   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_carb11a)
    CASE ('CARB10    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_carb10)
    CASE ('CARB13    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_carb13)
    CASE ('CARB16    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_carb16)
    CASE ('CARB9     ')
      d0(j) = d_h2o * SQRT(m_h2o / m_carb9)
    CASE ('CARB12    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_carb12)
    CASE ('CARB15    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_carb15)
    CASE ('UCARB12   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_ucarb12)
    CASE ('NUCARB12  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_nucarb12)
    CASE ('UDCARB8   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_udcarb8)
    CASE ('UDCARB11  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_udcarb11)
    CASE ('UDCARB14  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_udcarb14)
    CASE ('TNCARB26  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_tncarb26)
    CASE ('TNCARB10  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_tncarb10)
    CASE ('TNCARB12  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_tncarb12)
    CASE ('TNCARB11  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_tncarb11)
    CASE ('CCARB12   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_ccarb12)
    CASE ('TNCARB15  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_tncarb15)
    CASE ('RCOOH25   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rcooh25)
    CASE ('TXCARB24  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_txcarb24)
    CASE ('TXCARB22  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_txcarb22)
    CASE ('RU12PAN   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_ru12pan)
    CASE ('RTN26PAN  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_rtn26pan)
    CASE ('AROH14    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_aroh14)
    CASE ('AROH17    ')
      d0(j) = d_h2o * SQRT(m_h2o / m_aroh17)
    CASE ('ARNOH14   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_arnoh14)
    CASE ('ARNOH17   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_arnoh17)
    CASE ('ANHY      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_anhy)
    CASE ('IEPOX     ')
      d0(j) = d_h2o * SQRT(m_h2o / m_iepox)
    CASE ('HMML      ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hmml)
    CASE ('HUCARB9   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hucarb9)
    CASE ('HPUCARB12 ')
      d0(j) = d_h2o * SQRT(m_h2o / m_hpucarb12)
    CASE ('DHPCARB9  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_dhpcarb9)
    CASE ('DHPR12OOH ')
      d0(j) = d_h2o * SQRT(m_h2o / m_dhpr12ooh)
    CASE ('DHCARB9   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_dhcarb9)
    CASE ('RU12NO3   ')
      d0(j) = d_h2o * SQRT(m_h2o / m_ru12no3)
    CASE ('RU10NO3  ')
      d0(j) = d_h2o * SQRT(m_h2o / m_ru10no3)
      ! RU12OOH happens to have same mass as default (150g/mol)
    CASE ('RU12OOH   ')
         ! m_cri = 150 (same as sec_org) => d0(j) = 0.72e-5
      d0(j) = d_h2o * SQRT(m_h2o / m_cri)
    END SELECT
  END DO
  !
  DO j = 1, ndepd
    IF (d0(j) < 0.0) THEN
      cmessage = 'Warning: No dry deposition for ' //                          &
                  speci(nldepd(j)) // ' will be calculated.'
      icode = -1
      CALL ereport(routinename, icode, cmessage)
    END IF
  END DO

  !       Set up other constants
  reus = -cp / (vkman * gg)

  first = .FALSE.

END IF

!     Calculate air density, rho-air. Set ustar to minimum value if
!     undefined

DO k = 1, rows
  DO i = 1, row_length
    rho_air(i,k) = m_air * p0(i,k) / (rmol * t0(i,k))
    ustar(i,k) = u_s(i,k)
    IF (ustar(i,k) <= 0.0) ustar(i,k) = 1.0e-2
  END DO
END DO

!     Set zero-plane displacement to 0.7 of canopy height
!     (Smith et al., 2000).

DO n = 1, npft
  DO k = 1, rows
    DO i = 1, row_length
      IF (canht(i,k,n) > 0.0) THEN
        d(i,k,n) = canht(i,k,n) * 0.70
      ELSE
        d(i,k,n) = 0.0
      END IF
    END DO
  END DO
END DO

!     Ensure undefined parts of Z0 are set to 0.001 to avoid any
!     floating-point exceptions. Initialise z to reference height.

DO n = 1, ntype
  DO k = 1, rows
    DO i = 1, row_length
      IF (z0tile(i,k,n) <= 0.0) z0tile(i,k,n) = 0.001
      z(i,k,n) = zref
    END DO
  END DO
END DO

!     Set Z to height above surface minus zero plane displacement
!     for vegetated tiles.

DO n = 1, npft
  DO k = 1, rows
    DO i = 1, row_length
      z(i,k,n) = zref - d(i,k,n)
    END DO
  END DO
END DO

!     Calculate kinematic viscosity of air,
!     KVA; = dynamic viscosity / density
!     Formula from Kaye and Laby, 1952, p.43.

DO k = 1, rows
  DO i = 1, row_length
    kva(i,k) = (vair296 - 4.83e-8 * (296.0-t0(i,k))) /                         &
                rho_air(i,k)

    !         Calculate roughness length over oceans using Charnock
    !         formula in form given by Hummelshoj et al., 1992.

    IF (gsf(i,k,nwater) > 0.0) THEN
      z0tile(i,k,nwater) = (kva(i,k) / (9.1 * ustar(i,k))) +                   &
        0.016 * ustar(i,k) * ustar(i,k) / gg
    END IF
  END DO
END DO

!     Calculate Monin-Obukhov length L

DO k = 1, rows
  DO i = 1, row_length
    IF (hflx(i,k) /= 0.0) THEN
      !           Stable or unstable b.l.
      l(i,k) = reus * t0(i,k) * rho_air(i,k) * (ustar(i,k)**3) /               &
        hflx(i,k)
    ELSE
      !           Neutral b.l.
      l(i,k) = 10000.0
    END IF
  END DO
END DO

!     Calculate Businger functions (PSI(ZETA))

DO k = 1, rows
  DO i = 1, row_length
    IF (l(i,k) > 0.0) THEN
      !           Stable b.l.
      !CDIR EXPAND=9
      DO n = 1, ntype
        psi(i,k,n) = -5.0 * (z(i,k,n) - z0tile(i,k,n)) / l(i,k)
      END DO
    ELSE
      !           Unstable b.l.
      !CDIR EXPAND=9
      DO n = 1, ntype
        b1 = (1.0 - 16.0 * (z(i,k,n) / l(i,k)))**0.5
        b0 = (1.0 - 16.0 * (z0tile(i,k,n) / l(i,k)))**0.5
        psi(i,k,n) = 2.0 * LOG((1.0+b1) / (1.0+b0))
      END DO
    END IF
  END DO
END DO

!   1. Resa.  Calculate aerodynamic resistance (Resa) for all species

DO n = 1, ntype
  DO k = 1, rows
    DO i = 1, row_length
      resa(i,k,n) = (LOG(z(i,k,n) / z0tile(i,k,n)) - psi(i,k,n)) /             &
        (vkman * ustar(i,k))
    END DO
  END DO
END DO

!   2. Rb

!    Calculate Schmidt number divided by Prandtl number for each
!    species. Schmidt number = Kinematic viscosity of air divided
!    by molecular diffusivity. Calculate quasi-laminar boundary
!    layer resistance (Hicks et al., 1987)

DO j = 1, ndepd

  ! First set rb for aerosols to 1.0 as done in STOCHEM code
  IF ( (speci(nldepd(j)) == 'ORGNIT    ') .OR.                                 &
       (speci(nldepd(j)) == 'ASOA      ') .OR.                                 &
       (speci(nldepd(j)) == 'BSOA      ') .OR.                                 &
       (speci(nldepd(j)) == 'ISOSOA    ')                                      &
     ) THEN
    DO k = 1, rows
      DO i = 1, row_length
        rb(i,k,j) =  1.0
      END DO
    END DO

  ELSE IF (d0(j) > 0.0) THEN
    DO k = 1, rows
      DO i = 1, row_length
        sc_pr = kva(i,k) / (pr * d0(j))
        rb(i,k,j) = (sc_pr**twoth) / (vkman * ustar(i,k))
      END DO
    END DO
  ELSE
    DO k = 1, rows
      DO i = 1, row_length
        rb(i,k,j) = 0.0
      END DO
    END DO
  END IF
END DO

!    3. SO4_VD: Needed for ORGNIT or other species treated as aerosols
!
!     Calculate surface deposition velocity term for sulphate particles
!     using simple parameterisation of Wesely in the form given by
!     Zhang et al. (2001), Atmos Environ, 35, 549-560.
!     These values can be used for different aerosol types
!     and ORGNIT is considered to be an aerosol.
!
so4_vd = 1.0 / r_null

outer_loop: DO j = 1, ndepd
  IF ( (speci(nldepd(j)) == 'ORGNIT    ') .OR.                                 &
       (speci(nldepd(j)) == 'ASOA      ') .OR.                                 &
       (speci(nldepd(j)) == 'BSOA      ') .OR.                                 &
       (speci(nldepd(j)) == 'ISOSOA    ')                                      &
     ) THEN

    DO k = 1, rows
      DO i = 1, row_length
        so4_vd(i,k) = 0.0
        bl_l = zbl(i,k) / l(i,k)
        IF (l(i,k) >= 0.0) THEN
          so4_vd(i,k) = 0.002 * ustar(i,k)
        END IF
        IF (bl_l < -30.0) THEN
          so4_vd(i,k) = 0.0009 * ustar(i,k) * ((-bl_l)**twoth)
        END IF
        IF (bl_l >= -30.0 .AND. l(i,k) < 0.0) THEN
          so4_vd(i,k) = 0.002*ustar(i,k)*(1.0 +                                &
                           (-300.0/l(i,k))**twoth)
        END IF
      END DO
    END DO

    EXIT outer_loop
  END IF
END DO outer_loop

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_aerod
END MODULE ukca_aerod_mod
