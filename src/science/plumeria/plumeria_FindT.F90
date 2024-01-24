! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module containing the subroutines plumeria_FindT and
!    plumeria_FindT_init.
!
!  Methods:
!    Description of each subroutine contained here:
!    1. plumeria_FindT: subroutine that finds the mixture temperature
!              given a mixture enthalpy and composition.
!    2. plumeria_FindT_init: subroutine that finds the mixture temperature
!                   at the vent given a mixture enthalpy and composition.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE plumeria_FindT_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PLUMERIA_FINDT_MOD'

CONTAINS

!=======================================================================!

SUBROUTINE plumeria_FindT_init(m_m, m_a, m_w, h_mix, pnow, Tmix, m_v, m_l, m_i)
  ! ---------------------------------------------------------------------
  ! Description:
  !   Finds the mixture temperature at the vent given a mixture
  !   enthalpy and composition.
  ! ---------------------------------------------------------------------

USE plumeria_param_mod,     ONLY: kgmole_w, kgmole_air, Cp_m
USE plumeria_functions_mod, ONLY: plumeria_psat, plumeria_h_a, plumeria_cp_a,  &
                                  plumeria_cp_l, plumeria_h_l, plumeria_h_m,   &
                                  plumeria_h_v, plumeria_Tsat
USE umPrintMgr,             ONLY: umPrint, umMessage
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,               ONLY: jpim, jprb      ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Input variables
REAL, INTENT(IN)  :: m_m   ! Mass fraction magma in column
REAL, INTENT(IN)  :: m_a   ! Mass fraction dry air in column
REAL, INTENT(IN)  :: m_w   ! Total water in column
REAL, INTENT(IN)  :: h_mix ! Mixture enthalpy
REAL, INTENT(IN)  :: pnow  ! Partial pressure of the column (Pa)
REAL, INTENT(OUT) :: Tmix  ! Temperature of mixture column (K)
REAL, INTENT(OUT) :: m_v   ! Mass fraction water vapor in column
REAL, INTENT(OUT) :: m_l   ! Mass fraction liquid water in column
REAL, INTENT(OUT) :: m_i   ! Mass fraction ice in water column

! Local real variables
REAL :: Cp_mix       ! Specific heat of mixture column
REAL :: Cpa_avg      ! Average value of cp_air between 270 and 1000 K
REAL :: Cpwi_avg     ! Average cp for ice
REAL :: Cpwl_avg     ! Approximate average water specific heat
REAL :: Cpwv_avg     ! Approximate average specific heat of water vapor
                     ! between 100 and 900 C
REAL  :: H_boiling   ! Enthalpy at boiling, assumine all water is liquid
REAL  :: H_ColdWater ! Enthalpy at freezing, assuming all water is liquid
REAL  :: H_toboil    ! Enthalpy at boiling, assumine all water is vapor
REAL  :: hmixnow     ! Enthalpy of mixture column
REAL  :: m_vfreezing ! Mass fraction of water vapor at freezing point
REAL  :: Tboil       ! Saturation temperature of water
REAL  :: w_freezing  ! Specific humidity at the freezing point

INTEGER                            :: errcode    ! Error code for ereport
CHARACTER (LEN=errormessagelength) :: cmessage   ! Error message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PLUMERIA_FINDT_INIT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Set default values
Cpa_avg  = 1060.0
Cpwv_avg = 2150.0
Cpwl_avg = 4190.0
Cpwi_avg = 1850.0

Tboil = plumeria_Tsat(pnow)

H_toboil = ( m_m * plumeria_h_m(Tboil, pnow) )  +                              &
           ( m_a * plumeria_h_a(Tboil) )        +                              &
           ( m_w * plumeria_h_v(Tboil) )

H_boiling = ( m_m * plumeria_h_m(Tboil, pnow) ) +                              &
            ( m_a * plumeria_h_a(Tboil) )       +                              &
            ( m_w * plumeria_h_l(Tboil) )

w_freezing = ( kgmole_w / kgmole_air ) *                                       &
             ( plumeria_psat(273.15) / ( pnow - plumeria_psat(273.15) ) )

m_vfreezing = m_a * ( w_freezing / ( 1.0 + w_freezing ) )

IF (m_vfreezing > m_w) THEN
  m_vfreezing = m_w
END IF

H_ColdWater = ( m_m * plumeria_h_m(273.15, pnow) )   +                         &
              ( m_a * plumeria_h_a(273.15) )         +                         &
              ( m_vfreezing * plumeria_h_v(273.15) ) +                         &
              ( ( m_w - m_vfreezing ) * plumeria_h_l(273.15) )

!If we're above the boiling regime
IF (h_mix > H_toboil) THEN
  m_v = m_w
  m_l = 0.0
  m_i = 0.0
  Cp_mix = ( m_m * Cp_m )    +                                                 &
           ( m_a * Cpa_avg ) +                                                 &
           ( m_v * Cpwv_avg )
  !Estimate temperature, based on average specific heats.
  Tmix = Tboil + (h_mix - H_toboil) / Cp_mix
  hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                               &
            ( m_v * plumeria_h_v(Tmix) )       +                               &
            ( m_a * plumeria_h_a(Tmix) )
  !Iterate on final solution
  check_hmixnow1: DO
    IF ((ABS(hmixnow - h_mix) / h_mix) >= 0.001) THEN
      Tmix = Tmix + (h_mix - hmixnow) / Cp_mix
      hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                           &
                ( m_v * plumeria_h_v(Tmix) )       +                           &
                ( m_a * plumeria_h_a(Tmix) )
    ELSE
      EXIT check_hmixnow1
    END IF
  END DO check_hmixnow1

  !If we're within the boiling regime
ELSE IF (h_mix > H_boiling) THEN
  Tmix = Tboil
  m_v = m_w * (h_mix - H_boiling) / (H_toboil - H_boiling)
  m_l = m_w - m_v
  m_i = 0.0

  !If we're below the boiling regime but above the freezing regime
ELSE IF (h_mix > H_ColdWater) THEN
  m_l = m_w
  m_i = 0.0
  m_v = 0.0
  Cp_mix = ( m_m * Cp_m )    +                                                 &
           ( m_a * Cpa_avg ) +                                                 &
           ( m_w * Cpwl_avg )
  !Take a first stab at temperature
  Tmix = 273.15 + (h_mix - H_ColdWater) / Cp_mix
  hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                               &
            ( m_l * plumeria_h_l(Tmix) )       +                               &
            ( m_a * plumeria_h_a(Tmix) )
  !Iterate on final solution
  check_hmixnow2: DO
    IF ((ABS(hmixnow - h_mix) / h_mix) >= 0.001) THEN
      Cp_mix = ( m_m * Cp_m )                +                                 &
               ( m_a * plumeria_cp_a(Tmix) ) +                                 &
               ( m_l * plumeria_cp_l(Tmix) )
      Tmix = Tmix + (h_mix - hmixnow) / Cp_mix
      hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                           &
                ( m_v * plumeria_h_v(Tmix) )       +                           &
                ( m_l * plumeria_h_l(Tmix) )       +                           &
                ( m_a * plumeria_h_a(Tmix) )
    ELSE
      EXIT check_hmixnow2
    END IF
  END DO check_hmixnow2
END IF

END SUBROUTINE plumeria_FindT_init

!=======================================================================!

SUBROUTINE plumeria_FindT(m_m, m_a, m_w, h_mix, pnow, Tmix, m_v, m_l, m_i)
  ! ---------------------------------------------------------------------
  ! Description:
  !   Finds the mixture temperature given a mixture
  !   enthalpy and composition.
  ! ---------------------------------------------------------------------

USE plumeria_param_mod,     ONLY: T_ice, T_ColdWater, kgmole_w,                &
                                  kgmole_air, Cp_m
USE plumeria_functions_mod, ONLY: plumeria_psat, plumeria_h_a, plumeria_h_i,   &
                                  plumeria_h_l, plumeria_h_m, plumeria_h_v,    &
                                  plumeria_Tsat
USE umPrintMgr,             ONLY: umPrint, umMessage
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,               ONLY: jpim, jprb      ! DrHook
USE yomhook,                ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Input variables
REAL, INTENT (IN) :: m_m    ! Mass fraction magma in column
REAL, INTENT (IN) :: m_a    ! Mass fraction dry air in column
REAL, INTENT (IN) :: m_w    ! Total water in column
REAL, INTENT (IN) :: h_mix  ! Mixture enthalpy
REAL, INTENT (IN) :: pnow   ! Partial pressure of the column (Pa)
REAL, INTENT (OUT):: Tmix   ! Temperature of column (K)
REAL, INTENT (OUT):: m_v    ! Mass fraction water vapor in column
REAL, INTENT (OUT):: m_l    ! Mass fraction liquid water in column
REAL, INTENT (OUT):: m_i    ! Mass fraction ice in water column

! Local variables
REAL :: Cp_mix      ! Specific heat of mixture column
REAL :: Cpa_avg     ! Average value of cp_air between 270 and 1000 K
REAL :: Cpwi_avg    ! Average cp for ice
REAL :: Cpwl_avg    ! Approximate average water specific heat
REAL :: Cpwv_avg    ! Approximate average specific heat of water vapor
                    ! between 100 and 900 C
REAL :: H_ColdWater ! Enthalpy at freezing, assuming all water is liquid
REAL :: H_freezing  ! Enthalpy at freezing, asssuming all water is ice
REAL :: H_sat       ! Enthalpy at saturation, assumine all water is vapor
REAL :: T_sat       ! Saturation temperature
REAL :: hmixnow     ! Enthalpy of mixture column
REAL :: m_vColdWater! Mass fraction of water vapor at T_ColdWater
REAL :: m_vfreezing ! Mass fraction of water vapor at freezing point
REAL :: w_ColdWater ! Specific humidity at T_ColdWater
REAL :: w_freezing  ! Specific humidity at freezing point
REAL :: w_s         ! Specific humidity at T_mix
REAL :: x_a         ! Mole fraction air at saturation
REAL :: x_w         ! Mole fraction water vapor at saturation

! Local variables used in calculation
REAL :: hdif        ! Difference in enthalpy
REAL :: Hmixhi      ! Enthalpy of mixture, higher limit
REAL :: Hmixlo      ! Enthalpy of mixture column, lower limit
REAL :: TmixHi      ! Temperature at mixture column, higher limit
REAL :: TmixLst     ! Last value of Tmix in the iteration
REAL :: TmixLo      ! Temperature of mixture column, lower limit
INTEGER :: i

INTEGER                            :: errcode    ! Error code for ereport
CHARACTER (LEN=errormessagelength) :: cmessage   ! Error message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PLUMERIA_FINDT'

! end of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Set default values
kgmole_air = 8.314 / 286.98         ! Molar weight of air
kgmole_w = 0.0180152                ! Molar weight of water
Cpa_avg = 1060.0
Cpwv_avg = 2150.0
Cpwl_avg = 4190.0
Cpwi_avg = 1850.0

x_w = (m_w / kgmole_w) /                                                       &
       ( (m_w / kgmole_w) + (m_a / kgmole_air) )

x_a = 1.0 - x_w

T_sat = plumeria_Tsat( x_w * pnow )

!Enthalpy at saturation, assumine all water is vapor
H_sat = ( m_m * plumeria_h_m(T_sat, pnow) ) +                                  &
        ( m_a * plumeria_h_a(T_sat) )       +                                  &
        ( m_w * plumeria_h_v(T_sat) )

!FIND H_freezing
!if the ambient pressure exceeds the boiling pressure at T_ice
IF (pnow > plumeria_psat(T_ice)) THEN
  w_freezing = ( kgmole_w / kgmole_air ) *                                     &
               ( plumeria_psat(T_ice) /                                        &
                 (pnow - plumeria_psat(T_ice) ) )
  m_vfreezing = m_a * w_freezing

  !if there!s enough water vapor to saturate the plume
  IF ( m_w > m_vfreezing ) THEN
    H_freezing = ( m_m * plumeria_h_m(T_ice, pnow) )     +                     &
                 ( m_a * plumeria_h_a(T_ice) )           +                     &
                 ( m_vfreezing * plumeria_h_v(T_ice) )   +                     &
                 ( ( m_w - m_vfreezing ) * plumeria_h_i(T_ice) )
    !if air is non water-saturated at T_ice
  ELSE
    m_vfreezing = m_w
    H_freezing = ( m_m * plumeria_h_m(T_ice, pnow) ) +                         &
                 ( m_a * plumeria_h_a(T_ice) )       +                         &
                 ( m_vfreezing * plumeria_h_v(T_ice) )
  END IF
ELSE
  !if pnow is less than the boiling pressure at T_ice
  m_vfreezing = m_w
  H_freezing = ( m_m * plumeria_h_m(T_ice, pnow) )   +                         &
               ( m_a * plumeria_h_a(T_ice) )         +                         &
               ( m_vfreezing * plumeria_h_v(T_ice) )
END IF

!FIND H_coldwater
IF (pnow > plumeria_psat(T_ColdWater)) THEN
  !if pnow exceeds the boiling pressure at T_ColdWater (i.e. pnow > psat)
  w_ColdWater = ( kgmole_w / kgmole_air ) *                                    &
                ( plumeria_psat(T_ColdWater) /                                 &
                  ( pnow - plumeria_psat(T_ColdWater) ) )
  m_vColdWater = m_a * w_ColdWater
  !if there!s enough water vapor to saturate the plume
  IF (m_w > m_vColdWater) THEN
    !Enthalpy at top of mixed water-ice temperature range
    H_ColdWater = ( m_m * plumeria_h_m(T_ColdWater, pnow) )    +               &
                  ( m_a * plumeria_h_a(T_ColdWater) )          +               &
                  ( m_vColdWater * plumeria_h_v(T_ColdWater) ) +               &
                  ( (m_w - m_vColdWater) * plumeria_h_l(T_ColdWater) )
  ELSE
    !if there!s not enough water vapor to saturate the plume
    m_vColdWater = m_w
    H_ColdWater = ( m_m * plumeria_h_m(T_ColdWater, pnow) ) +                  &
                  ( m_a * plumeria_h_a(T_ColdWater) )       +                  &
                  ( m_vColdWater * plumeria_h_v(T_ColdWater) )
  END IF
ELSE
  !if we're above the boiling point at this pressure
  m_vColdWater = m_w
  H_ColdWater = ( m_m * plumeria_h_m(T_ColdWater, pnow) ) +                    &
                ( m_a * plumeria_h_a(T_ColdWater) )       +                    &
                ( m_vColdWater * plumeria_h_v(T_ColdWater) )
END IF

!The following criteria are met only as a result of inaccuracies
!in calculating T_sat, when it!s close
!to 273.15 K.  if (T_sat < 273.15 K) and (h_mix > H_freezing)
!and (h_mix < h_sat), the program tries to calculate
!the enthalpy of liquid water at T < 273.15, which causes it to blow up.
IF ((h_mix < H_sat) .AND. (h_mix > H_freezing) .AND. (T_sat < 273.15)) THEN
  T_sat = 273.15
END IF

IF ((h_mix > H_sat) .AND. (h_mix > H_ColdWater)) THEN
  !if we're above water saturation, and above freezing
  m_v = m_w
  m_l = 0.0
  m_i = 0.0
  Cp_mix = ( m_m * Cp_m )    +                                                 &
           ( m_a * Cpa_avg ) +                                                 &
           ( m_v * Cpwv_avg )
  !Estimate temperature, based on average specific heats.
  IF (H_sat > H_ColdWater) THEN
    Tmix = T_sat + (h_mix - H_sat) / Cp_mix
  ELSE
    Tmix = T_ColdWater + (h_mix - H_ColdWater) / Cp_mix
  END IF
  hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                               &
            ( m_v * plumeria_h_v(Tmix) )       +                               &
            ( m_a * plumeria_h_a(Tmix) )
  check_hmixnow3: DO
    IF (ABS(hmixnow - h_mix) / h_mix > 0.001) THEN
      !Iterate on final solution
      Tmix = Tmix + (h_mix - hmixnow) / Cp_mix
      hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                           &
                ( m_v * plumeria_h_v(Tmix) )       +                           &
                ( m_a * plumeria_h_a(Tmix) )
    ELSE
      EXIT check_hmixnow3
    END IF
  END DO check_hmixnow3

  !if we're within the saturated regime but above the freezing regime
ELSE IF ((h_mix > H_ColdWater) .AND. (T_sat > T_ColdWater)) THEN
  i = 1
  m_i = 0.0
  !Take a first stab at temperature
  Tmix = T_ColdWater + (T_sat - T_ColdWater) *                                 &
         (h_mix - H_ColdWater) / (H_sat - H_ColdWater)
  w_s = ( kgmole_w / kgmole_air ) *                                            &
        ( plumeria_psat(Tmix)     /                                            &
          ( pnow - plumeria_psat(Tmix) ) )

  !if we're water saturated
  IF (m_w >= m_a * w_s) THEN
    m_v = m_a * w_s
    m_l = m_w - m_v
  ELSE
    m_v = m_w
    m_l = 0.0
  END IF

  Hmixhi = H_sat
  TmixHi = T_sat
  Hmixlo = H_ColdWater
  TmixLo = T_ColdWater

  IF (Hmixhi < Hmixlo) THEN
    cmessage = 'Problem with enthalpy calculations. Hmixhi < Hmixlo.'
    errcode = 1
    CALL ereport(RoutineName,errcode,cmessage)
  END IF

  hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                               &
            ( m_v * plumeria_h_v(Tmix) )       +                               &
            ( m_l * plumeria_h_l(Tmix) )       +                               &
            ( m_a * plumeria_h_a(Tmix) )
  !Iterate on final solution
  check_hmixnow4: DO
    IF ((ABS(hmixnow - h_mix) / h_mix > 0.001) .AND. i < 100) THEN

      IF (hmixnow > h_mix) THEN
        Hmixhi = hmixnow
        TmixHi = Tmix
      ELSE
        Hmixlo = hmixnow
        TmixLo = Tmix
      END IF

      IF (Hmixhi < Hmixlo) THEN
        cmessage = 'Problem with enthalpy calculations. Hmixhi < Hmixlo.'
        errcode = 2
        CALL ereport(RoutineName,errcode,cmessage)
      END IF

      Tmix = TmixLo + (TmixHi - TmixLo) * (h_mix - Hmixlo) / (Hmixhi - Hmixlo)

      IF (Tmix < T_ColdWater) THEN
        cmessage = 'Problem with enthalpy calculations. T_mix < T_ColdWater.'
        errcode = 3
        CALL ereport(RoutineName,errcode,cmessage)
      END IF

      !recalculate m_v
      w_s = ( kgmole_w / kgmole_air ) *                                        &
            ( plumeria_psat(Tmix) / ( pnow - plumeria_psat(Tmix) ) )

      IF (m_w >= m_a * w_s) THEN
        m_v = m_a * w_s
        m_l = m_w - m_v
      ELSE
        m_v = m_w
        m_l = 0.0
      END IF

      hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                           &
                ( m_v * plumeria_h_v(Tmix) )       +                           &
                ( m_l * plumeria_h_l(Tmix) )       +                           &
                ( m_a * plumeria_h_a(Tmix) )
      i = i + 1

    ELSE
      EXIT check_hmixnow4
    END IF
  END DO check_hmixnow4

  !if we're at T_ice < tnow < T_ColdWater
ELSE IF (h_mix > H_freezing) THEN
  i = 1
  !Take a first stab at temperature
  Tmix = T_ice + (T_ColdWater-T_ice) * (h_mix-H_freezing) /                    &
         (H_ColdWater - H_freezing)
  w_s = ( kgmole_w / kgmole_air) *                                             &
        ( plumeria_psat(Tmix) / ( pnow - plumeria_psat(Tmix) ) )

  !w_s < 0 only very high, perhaps at z>40 km
  IF (w_s > 0.0) w_s = m_w / m_a

  !If it is water saturated
  IF (m_w >= m_a * w_s) THEN
    m_v = m_a * w_s
    m_l = (m_w - m_v) * (Tmix - T_ice) / (T_ColdWater - T_ice)
    IF (m_l > 1.0) THEN
      cmessage = 'Problem with enthalpy calculations. m_l > 1.'
      errcode = 4
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
    m_i = m_w - m_v - m_l
  ELSE
    m_v = m_w
    m_l = 0.0
    m_i = 0.0
  END IF

  Hmixhi = H_ColdWater
  TmixHi = T_ColdWater
  Hmixlo = H_freezing
  TmixLo = T_ice

  IF (Hmixhi < Hmixlo) THEN
    cmessage = 'Problem with enthalpy calculations. Hmixhi < Hmixlo.'
    errcode = 5
    CALL ereport(RoutineName,errcode,cmessage)
  END IF

  hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                               &
            ( m_v * plumeria_h_v(Tmix) )       +                               &
            ( m_l * plumeria_h_l(Tmix) )       +                               &
            ( m_i * plumeria_h_i(Tmix) )       +                               &
            ( m_a * plumeria_h_a(Tmix) )
  hdif = ABS(hmixnow - h_mix) / h_mix

  ! initialize TmixLst for the loop
  TmixLst = Tmix + 0.1

  !Iterate on final solution
  check_hmixnow5: DO
    IF ((ABS(hmixnow - h_mix) / h_mix > 0.001) .AND.                           &
       ((TmixHi - TmixLo) > 0.25) .AND. ABS(TmixLst - Tmix) > 0.1) THEN

      IF (hmixnow > h_mix) THEN
        Hmixhi = hmixnow
        TmixHi = Tmix
      ELSE
        Hmixlo = hmixnow
        TmixLo = Tmix
      END IF

      IF (Hmixhi < Hmixlo) THEN
        cmessage = 'Problem with enthalpy calculations. Hmixhi < Hmixlo.'
        errcode = 6
        CALL ereport(RoutineName,errcode,cmessage)
      END IF

      TmixLst = Tmix
      Tmix = TmixLo + (TmixHi-TmixLo)*(h_mix-Hmixlo)/(Hmixhi-Hmixlo)

      IF (Tmix < T_ice) THEN
        cmessage = 'Problem with enthalpy calculations. Tmix < T_ice.'
        errcode = 7
        CALL ereport(RoutineName,errcode,cmessage)
      END IF

      !recalculate m_v
      w_s = ( kgmole_w / kgmole_air ) *                                        &
            ( plumeria_psat(Tmix) / ( pnow - plumeria_psat(Tmix) ) )

      IF (m_w >= m_a * w_s) THEN
        m_v = m_a * w_s
        m_l = (m_w - m_v) * (Tmix - T_ice) / (T_ColdWater - T_ice)
        m_i = m_w - m_v - m_l
      ELSE
        m_v = m_w
        m_l = 0.0
        m_i = 0.0
      END IF

      hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                           &
                ( m_v * plumeria_h_v(Tmix) )       +                           &
                ( m_l * plumeria_h_l(Tmix) )       +                           &
                ( m_i * plumeria_h_i(Tmix) )       +                           &
                ( m_a * plumeria_h_a(Tmix) )
      i = i + 1
      hdif = ABS(hmixnow - h_mix) / h_mix

    ELSE
      EXIT check_hmixnow5
    END IF
  END DO check_hmixnow5

  !if we're below T_ice
ELSE
  i = 1
  m_l = 0.0
  Hmixhi = H_freezing
  Hmixlo = ( m_m * plumeria_h_m(100.0, pnow) ) +                               &
           ( m_w * plumeria_h_i(100.0) )       +                               &
           ( m_a * plumeria_h_a(100.0) )
  TmixHi = T_ice
  TmixLo = 100.0
  !Take a first stab at temperature
  Tmix = 100.0 + (TmixHi - TmixLo) * (h_mix - Hmixlo) / (Hmixhi - Hmixlo)
  w_s = ( kgmole_w / kgmole_air ) *                                            &
        ( plumeria_psat(Tmix) / ( pnow - plumeria_psat(Tmix) ) )
  IF (m_w >= m_a * w_s) THEN
    m_v = m_a * w_s
    m_i = m_w - m_v
  ELSE
    m_v = m_w
    m_i = 0.0
  END IF

  hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                               &
            ( m_v * plumeria_h_v(Tmix) )       +                               &
            ( m_i * plumeria_h_i(Tmix) )       +                               &
            ( m_a * plumeria_h_a(Tmix) )
  !Iterate on final solution
  check_hmixnow6: DO
    IF ((ABS(hmixnow - h_mix) / h_mix > 0.001) .AND.                           &
       ((TmixHi - TmixLo) > 0.25)) THEN

      IF (hmixnow > h_mix) THEN
        TmixHi = Tmix
        Hmixhi = hmixnow
      ELSE
        TmixLo = Tmix
        Hmixlo = hmixnow
      END IF

      Tmix = TmixLo + (TmixHi - TmixLo) * (h_mix - Hmixlo) / (Hmixhi - Hmixlo)

      ! recalculate m_v
      w_s = ( kgmole_w / kgmole_air ) *                                        &
            ( plumeria_psat(Tmix) / ( pnow - plumeria_psat(Tmix) ) )

      IF (m_w >= m_a * w_s) THEN
        m_v = m_a * w_s
        m_i = m_w - m_v
      ELSE
        m_v = m_w
        m_i = 0.0
      END IF

      hmixnow = ( m_m * plumeria_h_m(Tmix, pnow) ) +                           &
                ( m_v * plumeria_h_v(Tmix) )       +                           &
                ( m_i * plumeria_h_i(Tmix) )       +                           &
                ( m_a * plumeria_h_a(Tmix) )

      IF (i > 20) THEN
        cmessage = 'i > 20. Calculations stopped.'
        errcode = 8
        CALL ereport(RoutineName,errcode,cmessage)
      END IF

    ELSE
      EXIT check_hmixnow6
    END IF
  END DO check_hmixnow6
END IF

END SUBROUTINE plumeria_FindT
END MODULE plumeria_FindT_mod
