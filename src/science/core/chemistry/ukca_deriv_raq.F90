! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Perform a chemical integration using the Backward Euler
!  solver for the UKCA model with RAQ chemistry (based on STOCHEM).
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method: Chemistry integration (Backward Euler) for RAQ chemistry
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_deriv_raq_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_DERIV_RAQ_MOD'

CONTAINS

SUBROUTINE ukca_deriv_raq (nr_therm, n_be_calls,                               &
                    n_pnts, dts,                                               &
                    h2o, vol,                                                  &
                    rc, hrc, dd, dw, dj,                                       &
                    ldepd, ldepw,                                              &
                    y, ddflux, wdflux)

USE ukca_chemco_raq_mod, ONLY: in2o5, in2o5_ss

USE asad_mod, ONLY: jpspec, jphk, jppj, jpdd, jpdw
USE ukca_config_specification_mod, ONLY: ukca_config

USE ukca_config_constants_mod, ONLY: avogadro
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ukca_um_legacy_mod, ONLY: oneover_v

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: nr_therm   ! No. thermal reactions
INTEGER, INTENT(IN) :: n_be_calls ! No. chemical steps
INTEGER, INTENT(IN) :: n_pnts     ! Actual no. calculations

REAL, INTENT(IN)    :: dts                 ! BE timestep

REAL, INTENT(IN)    :: h2o (n_pnts)         ! h2o concn
REAL, INTENT(IN)    :: vol (n_pnts)         ! volume of gridboxes (cm3)

REAL, INTENT(IN)    :: rc (n_pnts,nr_therm) ! rxn rate coeffs
REAL, INTENT(IN)    :: hrc(n_pnts,jphk)     ! heterog rxn rate coeffs (s-1)
REAL, INTENT(IN)    :: dd (n_pnts,jpspec)   ! dry dep rates (s-1)
REAL, INTENT(IN)    :: dw (n_pnts,jpspec)   ! wet dep rates (s-1)
REAL, INTENT(IN)    :: dj (n_pnts,jppj)     ! photol rates

LOGICAL, INTENT(IN) :: ldepd (jpspec) ! T for dry dep species
LOGICAL, INTENT(IN) :: ldepw (jpspec) ! T for wet dep species

REAL, INTENT(IN OUT) :: y (n_pnts,jpspec)    ! species concn (molec cm-3)

REAL, INTENT(OUT)   :: ddflux (n_pnts,jpdd)
                                      ! dry deposition flux (mol s-1)
REAL, INTENT(OUT)   :: wdflux (n_pnts,jpdw)
                                      ! wet deposition flux (mol s-1)

!     Local variables

INTEGER :: i, j, n            ! loop variables
INTEGER :: icnt_dd, icnt_wd   ! indices for dry/wet dep fluxes

REAL :: recip_avogad          ! 1 / no.of molecules in 1 mole
REAL :: p(n_pnts)             ! production rate
REAL :: l(n_pnts)             ! loss rate
REAL :: yp(n_pnts,jpspec)     ! previous species concn
REAL :: tmp1(n_pnts)          ! work array for vector reciprocal
REAL :: tmp2(n_pnts)          ! work array for vector reciprocal

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DERIV_RAQ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

recip_avogad = 1.0 / avogadro

!     Initialise flux diagnostics to zero
tmp1 = 0.0
tmp2 = 0.0

ddflux (:, :) = 0.0
wdflux (:, :) = 0.0

DO n = 1, n_be_calls           ! loop over chemical timesteps

  DO j = 1,jpspec
    yp(:,j) = y(:,j)
  END DO

  !       Iteration start

  DO i=1,ukca_config%nit                   ! loop over iterations
    !
    ! This section written automatically by MECH9GEN
    ! from the file newmech9.txt
    ! with 178 equations and 59 species (including h2o).
    !
    ! Note: Loss terms of O3 [y(4)] have been copied manually for O3S [y(28)].
    ! After that (DD(:,4)  ) has to be changed to (DD(:,28)  ) there.
    !
    ! Note: After the automatic generation, both the loss of N2O5 and the
    ! formation of HNO3 following the heterogeneous hydrolysis of N2O5
    ! need to be added manually.
    !
    ! Note: After the automatic generation, the code has been re-written to
    ! split additions for chemical production and loss terms to six per line
    ! maximum.
    ! This significantly increases the speed of the additions on the IBM
    !
    !
    !          O(3P)        y( 1)

    p = 0.0                                                                    &
    +(dj(:,3)  *y(:,7 ))        +(dj(:,15) *y(:,6 ))                           &
    +(rc(:,7)  *y(:,2 ))        +(dj(:,1)  *y(:,4 ))
    l = 0.0                                                                    &
    +(rc(:,1)  )        +(rc(:,5)  *y(:,5 ))+(rc(:,16) *y(:,7 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:, 1) = (yp(:, 1)+dts*p)*tmp1

    !
    !          O(1D)        y( 2)
    p = 0.0                                                                    &
    +(dj(:,2)  *y(:,4 ))
    l = 0.0                                                                    &
    +(rc(:,7)  )        +(rc(:,8)  *h2o)

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:, 2) = (yp(:, 2)+dts*p)*tmp1
    !
    !          OH           y( 3)
    p = 0.0                                                                    &
    +(dj(:,23) *y(:,35))                                                       &
    +(dj(:,21) *y(:,52))        +(dj(:,22) *y(:,31))                           &
    +(dj(:,19) *y(:,21))        +(dj(:,20) *y(:,26))                           &
    +(dj(:,5)  *y(:,10))
    p = p                                                                      &
    +(dj(:,17) *y(:,17))                                                       &
    +(dj(:,4)  *y(:,11))        +(dj(:,4)  *y(:,11))                           &
    +(rc(:,137)*y(:,35)*y(:,3 ))+(rc(:,182)*y(:,15)*y(:,55))                   &
    +(rc(:,129)*y(:,4 )*y(:,34)*0.36)
    p = p                                                                      &
    +(rc(:,135)*y(:,31)*y(:,3 ))                                               &
    +(rc(:,123)*y(:,4 )*y(:,50)*0.28)+(rc(:,128)*y(:,4 )*y(:,29)*0.27)         &
    +(rc(:,102)*y(:,3 )*y(:,26))+(rc(:,104)*y(:,3 )*y(:,52))                   &
    +(rc(:,45) *y(:,3 )*y(:,17))
    p = p                                                                      &
    +(rc(:,100)*y(:,3 )*y(:,21))                                               &
    +(rc(:,17) *y(:,5 )*y(:,15))+(rc(:,34) *y(:,6 )*y(:,15))                   &
    +(rc(:,8)  *y(:,2 )*h2o*2.00)    +(rc(:,14) *y(:,15)*y(:,4 ))
    l = 0.0                                                                    &
    +(rc(:,177)*y(:,39))                                                       &
    +(rc(:,174)*y(:,54))+(rc(:,175)*y(:,54))+(rc(:,176)*y(:,58))               &
    +(rc(:,160)*y(:,32))+(rc(:,170)*y(:,58))
    l = l                                                                      &
    +(rc(:,172)*y(:,56))                                                       &
    +(rc(:,143)*y(:,34))+(rc(:,154)*y(:,45))+(rc(:,156)*y(:,47))               &
    +(rc(:,138)*y(:,33))+(rc(:,139)*y(:,57))
    l = l                                                                      &
    +(rc(:,141)*y(:,29))                                                       &
    +(rc(:,125)*y(:,50))+(rc(:,135)*y(:,31))+(rc(:,137)*y(:,35))               &
    +(rc(:,102)*y(:,26))+(rc(:,104)*y(:,52))
    l = l                                                                      &
    +(rc(:,109)*y(:,49))                                                       &
    +(rc(:,94) *y(:,27))+(rc(:,98) *y(:,23))+(rc(:,100)*y(:,21))               &
    +(rc(:,81) *y(:,51))+(rc(:,86) *y(:,53))
    l = l                                                                      &
    +(rc(:,92) *y(:,25))                                                       &
    +(rc(:,70) *y(:,13))+(rc(:,71) *y(:,19))+(rc(:,75) *y(:,22))               &
    +(rc(:,59) *y(:,12))+(rc(:,63) *y(:,41))
    l = l                                                                      &
    +(rc(:,66) *y(:,14))                                                       &
    +(rc(:,44) *y(:,21))+(rc(:,44) *y(:,26))+(rc(:,44) *y(:,52))               &
    +(rc(:,35) *y(:,10))+(rc(:,44) *y(:,17))
    l = l                                                                      &
    +(rc(:,45) *y(:,17))                                                       &
    +(rc(:,30) *y(:,15))+(rc(:,31) *y(:,11))
    l = l                                                                      &
    +(rc(:,33) *y(:,43))                                                       &
    +(rc(:,13) *y(:,4 ))+(rc(:,21) *y(:,7 ))+(rc(:,24) *y(:,9 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:, 3) = (yp(:, 3)+dts*p)*tmp1
    !
    !          O3           y( 4)
    p = 0.0                                                                    &
    +(rc(:,1)  *y(:,1 ))        +(rc(:,131)*y(:,20)*y(:,15)*0.3 )
    l = 0.0                                                                    &
    +(dj(:,1)  )        +(dj(:,2)  )        +(dd(:,4)  )                       &
    +(rc(:,124)*y(:,50))+(rc(:,128)*y(:,29))+(rc(:,129)*y(:,34))
    l = l                                                                      &
    +(rc(:,14) *y(:,15))+(rc(:,112)*y(:,49))+(rc(:,123)*y(:,50))               &
    +(rc(:,11) *y(:,5 ))+(rc(:,12) *y(:,7 ))+(rc(:,13) *y(:,3 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:, 4) = (yp(:, 4)+dts*p)*tmp1
    !
    !          NO           y( 5)
    p = 0.0                                                                    &
    +(dj(:,3)  *y(:,7 ))        +(dj(:,14) *y(:,6 ))                           &
    +(rc(:,16) *y(:,7 )*y(:,1 ))+(rc(:,19) *y(:,7 )*y(:,6 ))
    l = 0.0                                                                    &
    +(rc(:,173)*y(:,44))                                                       &
    +(rc(:,126)*y(:,40))+(rc(:,142)*y(:,46))+(rc(:,144)*y(:,48))               &
    +(rc(:,95) *y(:,36))+(rc(:,105)*y(:,37))
    l = l                                                                      &
    +(rc(:,110)*y(:,38))                                                       &
    +(rc(:,79) *y(:,20))+(rc(:,83) *y(:,24))+(rc(:,93) *y(:,30))               &
    +(rc(:,17) *y(:,15))
    l = l                                                                      &
    +(rc(:,60) *y(:,16))+(rc(:,72) *y(:,18))                                   &
    +(rc(:,5)  *y(:,1 ))+(rc(:,11) *y(:,4 ))+(rc(:,15) *y(:,6 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:, 5) = (yp(:, 5)+dts*p)*tmp1
    !
    !          NO3          y( 6)
    p = 0.0                                                                    &
    +(dj(:,16) *y(:,8 ))                                                       &
    +(rc(:,35) *y(:,3 )*y(:,10))+(rc(:,98) *y(:,3 )*y(:,23))                   &
    +(rc(:,12) *y(:,7 )*y(:,4 ))+(rc(:,29) *y(:,8 ))
    l = 0.0                                                                    &
    +(dj(:,14) )        +(dj(:,15) )        +(dw(:,6)  )                       &
    +(rc(:,158)*y(:,22))+(rc(:,159)*y(:,29))+(rc(:,178)*y(:,39))
    l = l                                                                      &
    +(rc(:,152)*y(:,51))+(rc(:,153)*y(:,49))+(rc(:,155)*y(:,50))               &
    +(rc(:,34) *y(:,15))+(rc(:,67) *y(:,14))+(rc(:,151)*y(:,19))
    l = l                                                                      &
    +(rc(:,27) *y(:,6 ))+(rc(:,27) *y(:,6 ))+(rc(:,32) *y(:,15))               &
    +(rc(:,15) *y(:,5 ))+(rc(:,19) *y(:,7 ))+(rc(:,20) *y(:,7 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:, 6) = (yp(:, 6)+dts*p)*tmp1
    !
    !          NO2          y( 7)
    p = 0.0                                                                    &
    +(dj(:,18) *y(:,23))                                                       &
    +(dj(:,15) *y(:,6 ))        +(dj(:,16) *y(:,8 ))                           &
    +(dj(:,5)  *y(:,10))        +(dj(:,11) *y(:,9 ))                           &
    +(rc(:,177)*y(:,3 )*y(:,39))
    p = p                                                                      &
    +(rc(:,178)*y(:,6 )*y(:,39)*2.00)                                          &
    +(rc(:,160)*y(:,32)*y(:,3 ))+(rc(:,173)*y(:,44)*y(:,5 ))                   &
    +(rc(:,154)*y(:,45)*y(:,3 ))+(rc(:,156)*y(:,47)*y(:,3 ))                   &
    +(rc(:,142)*y(:,46)*y(:,5 ))
    p = p                                                                      &
    +(rc(:,144)*y(:,48)*y(:,5 ))                                               &
    +(rc(:,110)*y(:,38)*y(:,5 ))+(rc(:,126)*y(:,40)*y(:,5 ))                   &
    +(rc(:,95) *y(:,36)*y(:,5 ))+(rc(:,105)*y(:,5 )*y(:,37))                   &
    +(rc(:,83) *y(:,24)*y(:,5 ))
    p = p                                                                      &
    +(rc(:,93) *y(:,30)*y(:,5 ))                                               &
    +(rc(:,78) *y(:,23))        +(rc(:,79) *y(:,20)*y(:,5 ))                   &
    +(rc(:,60) *y(:,5 )*y(:,16))+(rc(:,72) *y(:,18)*y(:,5 ))                   &
    +(rc(:,29) *y(:,8 ))
    p = p                                                                      &
    +(rc(:,34) *y(:,6 )*y(:,15))                                               &
    +(rc(:,24) *y(:,9 )*y(:,3 ))+(rc(:,27) *y(:,6 )*y(:,6 )*2.00)              &
    +(rc(:,19) *y(:,7 )*y(:,6 ))+(rc(:,23) *y(:,9 ))
    p = p                                                                      &
    +(rc(:,15) *y(:,5 )*y(:,6 )*2.00)+(rc(:,17) *y(:,5 )*y(:,15))              &
    +(rc(:,5)  *y(:,1 )*y(:,5 ))+(rc(:,11) *y(:,5 )*y(:,4 ))
    l = 0.0                                                                    &
    +(dj(:,3)  )        +(dd(:,7)  )                                           &
    +(rc(:,77) *y(:,20))+(rc(:,171)*y(:,42))+(rc(:,183)*y(:,55))
    l = l                                                                      &
    +(rc(:,20) *y(:,6 ))+(rc(:,21) *y(:,3 ))+(rc(:,22) *y(:,15))               &
    +(rc(:,12) *y(:,4 ))+(rc(:,16) *y(:,1 ))+(rc(:,19) *y(:,6 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:, 7) = (yp(:, 7)+dts*p)*tmp1
    !
    !          N2O5         y( 8)
    p = 0.0                                                                    &
    +(rc(:,20) *y(:,7 )*y(:,6 ))
    l = 0.0                                                                    &
    +(rc(:,29) )        +(dj(:,16) )        +(dw(:,8)  )

    ! Optionally add N2O5 loss by heterogeneous hydrolysis if this is on.
    ! N2O5 + H2O (aerosols) --> 2 HNO3.
    ! Note that if sea salt aerosol is considered, the sea-salt
    ! reaction is accounted for separately because the product (HNO3) is
    ! thrown away.
    IF (in2o5 > 0) THEN
      l = l + hrc (:,in2o5)
    END IF

    IF (in2o5_ss > 0) THEN
      l = l + hrc (:,in2o5_ss)
    END IF

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:, 8) = (yp(:, 8)+dts*p)*tmp1
    !
    !          HO2NO2       y( 9)
    p = 0.0                                                                    &
    +(rc(:,22) *y(:,7 )*y(:,15))
    l = 0.0                                                                    &
    +(dw(:,9)  )                                                               &
    +(rc(:,23) )        +(rc(:,24) *y(:,3 ))+(dj(:,11) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:, 9) = (yp(:, 9)+dts*p)*tmp1
    !
    !          HONO2        y(10)
    p = 0.0                                                                    &
    +(rc(:,152)*y(:,6 )*y(:,51))+(rc(:,158)*y(:,6 )*y(:,22))                   &
    +(rc(:,67) *y(:,6 )*y(:,14))+(rc(:,151)*y(:,6 )*y(:,19))                   &
    +(rc(:,21) *y(:,7 )*y(:,3 ))+(rc(:,32) *y(:,6 )*y(:,15))

    ! Optionally add formation of HNO3 by heterogeneous hydrolysis
    ! of N2O5 if this is on:  N2O5 + H2O (aerosols) --> 2 HNO3
    IF (in2o5 > 0) THEN
      p = p + hrc (:,in2o5) * y(:,8) * 2.0
    END IF

    l = 0.0                                                                    &
    +(dd(:,10) )                                                               &
    +(rc(:,35) *y(:,3 ))+(dj(:,5)  )        +(dw(:,10) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,10) = (yp(:,10)+dts*p)*tmp1
    !
    !          H2O2         y(11)
    p = 0.0                                                                    &
    +(rc(:,36) *y(:,15)*y(:,15))
    l = 0.0                                                                    &
    +(dd(:,11) )                                                               &
    +(rc(:,31) *y(:,3 ))+(dj(:,4)  )        +(dw(:,11) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,11) = (yp(:,11)+dts*p)*tmp1
    !
    !          CH4          y(12)
    p = 0.0                                                                    &
    +(rc(:,123)*y(:,4 )*y(:,50)*0.30)
    l = 0.0                                                                    &
    +(rc(:,59) *y(:,3 ))+(dd(:,12) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,12) = (yp(:,12)+dts*p)*tmp1
    !
    !          CO           y(13)
    p = 0.0                                                                    &
    +(dj(:,13) *y(:,57))                                                       &
    +(dj(:,12) *y(:,33))        +(dj(:,13) *y(:,57))                           &
    +(dj(:,7)  *y(:,14))        +(dj(:,8)  *y(:,22))                           &
    +(rc(:,139)*y(:,57)*y(:,3 )*2.0 )
    p = p                                                                      &
    +(dj(:,6)  *y(:,14))                                                       &
    +(rc(:,129)*y(:,4 )*y(:,34)*0.76)+(rc(:,138)*y(:,33)*y(:,3 ))              &
    +(rc(:,124)*y(:,4 )*y(:,50)*0.58)+(rc(:,128)*y(:,4 )*y(:,29)*0.78)
    p = p                                                                      &
    +(rc(:,112)*y(:,4 )*y(:,49)*0.31)+(rc(:,123)*y(:,4 )*y(:,50)*0.40)         &
    +(rc(:,66) *y(:,3 )*y(:,14))+(rc(:,67) *y(:,6 )*y(:,14))
    l = 0.0                                                                    &
    +(rc(:,70) *y(:,3 ))+(dd(:,13) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,13) = (yp(:,13)+dts*p)*tmp1
    !
    !          HCHO         y(14)
    p = 0.0                                                                    &
    +(dj(:,22) *y(:,31))        +(dj(:,23) *y(:,35))                           &
    +(rc(:,181)*y(:,44)*y(:,16))+(dj(:,17) *y(:,17))                           &
    +(rc(:,144)*y(:,48)*y(:,5 ))+(rc(:,154)*y(:,45)*y(:,3 ))
    p = p                                                                      &
    +(rc(:,137)*y(:,35)*y(:,3 ))+(rc(:,142)*y(:,46)*y(:,5 ))                   &
    +(rc(:,133)*y(:,48)*y(:,16)*1.0 )+(rc(:,135)*y(:,31)*y(:,3 ))              &
    +(rc(:,129)*y(:,4 )*y(:,34)*0.24)+(rc(:,132)*y(:,46)*y(:,16)*1.0 )
    p = p                                                                      &
    +(rc(:,127)*y(:,16)*y(:,40)*2.00)+(rc(:,128)*y(:,4 )*y(:,29)*0.22)         &
    +(rc(:,123)*y(:,4 )*y(:,50))+(rc(:,126)*y(:,40)*y(:,5 ))                   &
    +(rc(:,112)*y(:,4 )*y(:,49))+(rc(:,112)*y(:,4 )*y(:,49)*0.47)
    p = p                                                                      &
    +(rc(:,110)*y(:,38)*y(:,5 )*2.00)+(rc(:,111)*y(:,16)*y(:,38)*3.00)         &
    +(rc(:,98) *y(:,3 )*y(:,23))+(rc(:,106)*y(:,16)*y(:,37))                   &
    +(rc(:,96) *y(:,36)*y(:,16)*2.00)+(rc(:,97) *y(:,30)*y(:,16))
    p = p                                                                      &
    +(rc(:,84) *y(:,24)*y(:,16))+(rc(:,95) *y(:,36)*y(:,5 ))                   &
    +(rc(:,74) *y(:,16)*y(:,20)*2.00)+(rc(:,80) *y(:,16)*y(:,20))              &
    +(rc(:,63) *y(:,41)*y(:,3 ))
    p = p                                                                      &
    +(rc(:,73) *y(:,18)*y(:,16))                                               &
    +(rc(:,61) *y(:,16)*y(:,16)*2.00)+(rc(:,62) *y(:,16)*y(:,16))              &
    +(rc(:,45) *y(:,3 )*y(:,17))+(rc(:,60) *y(:,5 )*y(:,16))
    l = 0.0                                                                    &
    +(dj(:,7)  )        +(dw(:,14) )                                           &
    +(rc(:,66) *y(:,3 ))+(rc(:,67) *y(:,6 ))+(dj(:,6)  )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,14) = (yp(:,14)+dts*p)*tmp1
    !
    !          HO2          y(15)
    p = 0.0                                                                    &
    +(dj(:,22) *y(:,31))        +(dj(:,23) *y(:,35))                           &
    +(dj(:,20) *y(:,26))        +(dj(:,21) *y(:,52))                           &
    +(dj(:,17) *y(:,17))        +(dj(:,19) *y(:,21))
    p = p                                                                      &
    +(dj(:,13) *y(:,57))        +(dj(:,13) *y(:,57))                           &
    +(dj(:,11) *y(:,9 ))        +(dj(:,12) *y(:,33))                           &
    +(dj(:,6)  *y(:,14))        +(dj(:,8)  *y(:,22))
    p = p                                                                      &
    +(rc(:,181)*y(:,44)*y(:,16)*2.00)+(dj(:,6)  *y(:,14))                      &
    +(rc(:,173)*y(:,44)*y(:,5 ))+(rc(:,174)*y(:,3 )*y(:,54))                   &
    +(rc(:,159)*y(:,6 )*y(:,29))+(rc(:,170)*y(:,3 )*y(:,58))
    p = p                                                                      &
    +(rc(:,153)*y(:,6 )*y(:,49))+(rc(:,155)*y(:,6 )*y(:,50))                   &
    +(rc(:,142)*y(:,46)*y(:,5 ))+(rc(:,144)*y(:,48)*y(:,5 ))                   &
    +(rc(:,133)*y(:,48)*y(:,16)*2.0 )+(rc(:,139)*y(:,57)*y(:,3 )*1.0 )
    p = p                                                                      &
    +(rc(:,129)*y(:,4 )*y(:,34)*0.36)+(rc(:,132)*y(:,46)*y(:,16)*2.0 )         &
    +(rc(:,127)*y(:,16)*y(:,40)*2.00)+(rc(:,128)*y(:,4 )*y(:,29)*0.27)         &
    +(rc(:,124)*y(:,4 )*y(:,50)*0.18)+(rc(:,126)*y(:,40)*y(:,5 ))
    p = p                                                                      &
    +(rc(:,112)*y(:,4 )*y(:,49)*0.20)+(rc(:,123)*y(:,4 )*y(:,50)*0.30)         &
    +(rc(:,110)*y(:,38)*y(:,5 ))+(rc(:,111)*y(:,16)*y(:,38)*2.00)              &
    +(rc(:,97) *y(:,30)*y(:,16)*2.00)+(rc(:,106)*y(:,16)*y(:,37))
    p = p                                                                      &
    +(rc(:,93) *y(:,30)*y(:,5 ))+(rc(:,96) *y(:,36)*y(:,16))                   &
    +(rc(:,84) *y(:,24)*y(:,16)*2.00)+(rc(:,90) *y(:,18)*y(:,18)*2.00)         &
    +(rc(:,80) *y(:,16)*y(:,20))+(rc(:,83) *y(:,24)*y(:,5 ))
    p = p                                                                      &
    +(rc(:,72) *y(:,18)*y(:,5 ))+(rc(:,73) *y(:,18)*y(:,16)*2.00)              &
    +(rc(:,67) *y(:,6 )*y(:,14))+(rc(:,70) *y(:,3 )*y(:,13))                   &
    +(rc(:,63) *y(:,41)*y(:,3 ))+(rc(:,66) *y(:,3 )*y(:,14))
    p = p                                                                      &
    +(rc(:,60) *y(:,5 )*y(:,16))+(rc(:,61) *y(:,16)*y(:,16)*2.00)              &
    +(rc(:,31) *y(:,3 )*y(:,11))+(rc(:,33) *y(:,3 )*y(:,43))                   &
    +(rc(:,13) *y(:,3 )*y(:,4 ))+(rc(:,23) *y(:,9 ))
    l = 0.0                                                                    &
    +(rc(:,180)*y(:,42))+(rc(:,182)*y(:,55))+(dw(:,15) )                       &
    +(rc(:,131)*y(:,20))+(rc(:,134)*y(:,46))+(rc(:,136)*y(:,48))
    l = l                                                                      &
    +(rc(:,99) *y(:,18))+(rc(:,101)*y(:,30))+(rc(:,103)*y(:,24))               &
    +(rc(:,36) *y(:,15))+(rc(:,36) *y(:,15))+(rc(:,65) *y(:,16))
    l = l                                                                      &
    +(rc(:,30) *y(:,3 ))+(rc(:,32) *y(:,6 ))+(rc(:,34) *y(:,6 ))               &
    +(rc(:,14) *y(:,4 ))+(rc(:,17) *y(:,5 ))+(rc(:,22) *y(:,7 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,15) = (yp(:,15)+dts*p)*tmp1
    !
    !          MeOO         y(16)
    p = 0.0                                                                    &
    +(dj(:,10) *y(:,27))                                                       &
    +(rc(:,131)*y(:,20)*y(:,15)*0.8 )+(dj(:,8)  *y(:,22))                      &
    +(rc(:,91) *y(:,20)*y(:,20)*2.00)+(rc(:,123)*y(:,4 )*y(:,50)*0.58)
    p =p                                                                       &
    +(rc(:,79) *y(:,20)*y(:,5 ))+(rc(:,80) *y(:,16)*y(:,20))                   &
    +(rc(:,44) *y(:,3 )*y(:,17))+(rc(:,59) *y(:,3 )*y(:,12))
    l = 0.0                                                                    &
    +(dw(:,16) )                                                               &
    +(rc(:,132)*y(:,46))+(rc(:,133)*y(:,48))+(rc(:,181)*y(:,44))
    l = l                                                                      &
    +(rc(:,106)*y(:,37))+(rc(:,111)*y(:,38))+(rc(:,127)*y(:,40))               &
    +(rc(:,84) *y(:,24))+(rc(:,96) *y(:,36))+(rc(:,97) *y(:,30))
    l = l                                                                      &
    +(rc(:,73) *y(:,18))+(rc(:,74) *y(:,20))+(rc(:,80) *y(:,20))               &
    +(rc(:,62) *y(:,16))
    l = l                                                                      &
    +(rc(:,62) *y(:,16))+(rc(:,65) *y(:,15))                                   &
    +(rc(:,60) *y(:,5 ))+(rc(:,61) *y(:,16))+(rc(:,61) *y(:,16))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,16) = (yp(:,16)+dts*p)*tmp1
    !
    !          MeOOH        y(17)
    p = 0.0                                                                    &
    +(rc(:,65) *y(:,16)*y(:,15))
    l = 0.0                                                                    &
    +(dw(:,17) )        +(dd(:,17) )                                           &
    +(rc(:,44) *y(:,3 ))+(rc(:,45) *y(:,3 ))+(dj(:,17) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,17) = (yp(:,17)+dts*p)*tmp1
    !
    !          EtOO         y(18)
    p = 0.0                                                                    &
    +(rc(:,151)*y(:,6 )*y(:,19))+(dj(:,9)  *y(:,53))                           &
    +(rc(:,44) *y(:,3 )*y(:,21))+(rc(:,71) *y(:,3 )*y(:,19))
    l = 0.0                                                                    &
    +(rc(:,90) *y(:,18))+(rc(:,99) *y(:,15))                                   &
    +(rc(:,72) *y(:,5 ))+(rc(:,73) *y(:,16))+(rc(:,90) *y(:,18))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,18) = (yp(:,18)+dts*p)*tmp1
    !
    !          C2H6         y(19)
    p = 0.0
    l = 0.0                                                                    &
    +(rc(:,71) *y(:,3 ))+(rc(:,151)*y(:,6 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,19) = (yp(:,19)+dts*p)*tmp1
    !
    !          MeCO3        y(20)
    p = 0.0                                                                    &
    +(dj(:,18) *y(:,23))                                                       &
    +(dj(:,10) *y(:,27))        +(dj(:,12) *y(:,33))                           &
    +(rc(:,158)*y(:,6 )*y(:,22))+(dj(:,9)  *y(:,53))                           &
    +(rc(:,131)*y(:,20)*y(:,15)*0.2 )
    p = p                                                                      &
    +(rc(:,138)*y(:,33)*y(:,3 ))                                               &
    +(rc(:,105)*y(:,5 )*y(:,37))+(rc(:,106)*y(:,16)*y(:,37))
    p = p                                                                      &
    +(rc(:,95) *y(:,36)*y(:,5 ))+(rc(:,96) *y(:,36)*y(:,16))                   &
    +(rc(:,75) *y(:,3 )*y(:,22))+(rc(:,78) *y(:,23))
    l = 0.0                                                                    &
    +(rc(:,131)*y(:,15))                                                       &
    +(rc(:,80) *y(:,16))+(rc(:,91) *y(:,20))+(rc(:,91) *y(:,20))               &
    +(rc(:,74) *y(:,16))+(rc(:,77) *y(:,7 ))+(rc(:,79) *y(:,5 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,20) = (yp(:,20)+dts*p)*tmp1
    !
    !          EtOOH        y(21)
    p = 0.0                                                                    &
    +(rc(:,99) *y(:,15)*y(:,18))
    l = 0.0                                                                    &
    +(dw(:,21) )        +(dd(:,21) )                                           &
    +(rc(:,44) *y(:,3 ))+(rc(:,100)*y(:,3 ))+(dj(:,19) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,21) = (yp(:,21)+dts*p)*tmp1
    !
    !          MeCHO        y(22)
    p = 0.0                                                                    &
    +(dj(:,19) *y(:,21))                                                       &
    +(rc(:,127)*y(:,16)*y(:,40))+(rc(:,156)*y(:,47)*y(:,3 ))                   &
    +(rc(:,124)*y(:,4 )*y(:,50))+(rc(:,126)*y(:,40)*y(:,5 ))
    p = p                                                                      &
    +(rc(:,105)*y(:,5 )*y(:,37))+(rc(:,106)*y(:,16)*y(:,37))                   &
    +(rc(:,90) *y(:,18)*y(:,18)*2.00)+(rc(:,100)*y(:,3 )*y(:,21))              &
    +(rc(:,72) *y(:,18)*y(:,5 ))+(rc(:,73) *y(:,18)*y(:,16))
    l = 0.0                                                                    &
    +(rc(:,75) *y(:,3 ))+(rc(:,158)*y(:,6 ))+(dj(:,8)  )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,22) = (yp(:,22)+dts*p)*tmp1
    !
    !          PAN          y(23)
    p = 0.0                                                                    &
    +(rc(:,77) *y(:,20)*y(:,7 ))
    l = 0.0                                                                    &
    +(dd(:,23) )                                                               &
    +(rc(:,78) )        +(rc(:,98) *y(:,3 ))+(dj(:,18) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,23) = (yp(:,23)+dts*p)*tmp1
    !
    !          s-BuOO       y(24)
    p = 0.0                                                                    &
    +(rc(:,152)*y(:,6 )*y(:,51))                                               &
    +(rc(:,44) *y(:,3 )*y(:,52))+(rc(:,81) *y(:,3 )*y(:,51))
    l = 0.0                                                                    &
    +(rc(:,83) *y(:,5 ))+(rc(:,84) *y(:,16))+(rc(:,103)*y(:,15))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,24) = (yp(:,24)+dts*p)*tmp1
    !
    !          C3H8         y(25)
    p = 0.0
    l = 0.0                                                                    &
    +(rc(:,92) *y(:,3 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,25) = (yp(:,25)+dts*p)*tmp1
    !
    !          i-PrOOH      y(26)
    p = 0.0                                                                    &
    +(rc(:,101)*y(:,15)*y(:,30))
    l = 0.0                                                                    &
    +(dw(:,26) )        +(dd(:,26) )                                           &
    +(rc(:,44) *y(:,3 ))+(rc(:,102)*y(:,3 ))+(dj(:,20) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,26) = (yp(:,26)+dts*p)*tmp1
    !
    !          Me2CO        y(27)
    p = 0.0                                                                    &
    +(rc(:,102)*y(:,3 )*y(:,26))+(dj(:,20) *y(:,26))                           &
    +(rc(:,93) *y(:,30)*y(:,5 ))+(rc(:,97) *y(:,30)*y(:,16))
    l = 0.0                                                                    &
    +(rc(:,94) *y(:,3 ))+(dj(:,10) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,27) = (yp(:,27)+dts*p)*tmp1
    !
    !          O3S          y(28)
    p = 0.0
    l = 0.0                                                                    &
    +(dj(:,1)  )        +(dj(:,2)  )        +(dd(:,28)  )                      &
    +(rc(:,124)*y(:,50))+(rc(:,128)*y(:,29))+(rc(:,129)*y(:,34))
    l = l                                                                      &
    +(rc(:,14) *y(:,15))+(rc(:,112)*y(:,49))+(rc(:,123)*y(:,50))               &
    +(rc(:,11) *y(:,5 ))+(rc(:,12) *y(:,7 ))+(rc(:,13) *y(:,3 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,28) = (yp(:,28)+dts*p)*tmp1
    !
    !          C5H8         y(29)
    p = 0.0
    l = 0.0                                                                    &
    +(rc(:,128)*y(:,4 ))+(rc(:,141)*y(:,3 ))+(rc(:,159)*y(:,6 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,29) = (yp(:,29)+dts*p)*tmp1
    !
    !          i-PrOO       y(30)
    p = 0.0                                                                    &
    +(rc(:,44) *y(:,3 )*y(:,26))+(rc(:,92) *y(:,3 )*y(:,25))
    l = 0.0                                                                    &
    +(rc(:,93) *y(:,5 ))+(rc(:,97) *y(:,16))+(rc(:,101)*y(:,15))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,30) = (yp(:,30)+dts*p)*tmp1
    !
    !          ISOOH        y(31)
    p = 0.0                                                                    &
    +(rc(:,134)*y(:,46)*y(:,15))
    l = 0.0                                                                    &
    +(dd(:,31) )                                                               &
    +(rc(:,135)*y(:,3 ))+(dj(:,22) )        +(dw(:,31) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,31) = (yp(:,31)+dts*p)*tmp1
    !
    !          ISON         y(32)
    p = 0.0                                                                    &
    +(rc(:,159)*y(:,6 )*y(:,29))
    l = 0.0                                                                    &
    +(rc(:,160)*y(:,3 ))+(dw(:,32) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,32) = (yp(:,32)+dts*p)*tmp1
    !
    !          MGLY         y(33)
    p = 0.0                                                                    &
    +(dj(:,23) *y(:,35))                                                       &
    +(rc(:,180)*y(:,42)*y(:,15))+(rc(:,181)*y(:,44)*y(:,16))                   &
    +(rc(:,170)*y(:,3 )*y(:,58)*0.80)+(rc(:,173)*y(:,44)*y(:,5 ))
    p = p                                                                      &
    +(rc(:,137)*y(:,35)*y(:,3 ))+(rc(:,144)*y(:,48)*y(:,5 ))                   &
    +(rc(:,129)*y(:,4 )*y(:,34))+(rc(:,133)*y(:,48)*y(:,16)*1.0 )
    l = 0.0                                                                    &
    +(rc(:,138)*y(:,3 ))+(dj(:,12) )        +(dw(:,33) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,33) = (yp(:,33)+dts*p)*tmp1
    !
    !          MVK          y(34)
    p = 0.0                                                                    &
    +(rc(:,160)*y(:,32)*y(:,3 ))+(dj(:,22) *y(:,31))                           &
    +(rc(:,135)*y(:,31)*y(:,3 ))+(rc(:,142)*y(:,46)*y(:,5 ))                   &
    +(rc(:,128)*y(:,4 )*y(:,29))+(rc(:,132)*y(:,46)*y(:,16)*1.0 )
    l = 0.0                                                                    &
    +(rc(:,129)*y(:,4 ))+(rc(:,143)*y(:,3 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,34) = (yp(:,34)+dts*p)*tmp1
    !
    !          MVKOOH       y(35)
    p = 0.0                                                                    &
    +(rc(:,136)*y(:,48)*y(:,15))
    l = 0.0                                                                    &
    +(dd(:,35) )                                                               &
    +(rc(:,137)*y(:,3 ))+(dj(:,23) )        +(dw(:,35) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,35) = (yp(:,35)+dts*p)*tmp1
    !
    !          ACETO2       y(36)
    p = 0.0                                                                    &
    +(rc(:,94) *y(:,3 )*y(:,27))
    l = 0.0                                                                    &
    +(rc(:,95) *y(:,5 ))+(rc(:,96) *y(:,16))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,36) = (yp(:,36)+dts*p)*tmp1
    !
    !          MEKO2        y(37)
    p = 0.0                                                                    &
    +(rc(:,86) *y(:,3 )*y(:,53))
    l = 0.0                                                                    &
    +(rc(:,105)*y(:,5 ))+(rc(:,106)*y(:,16))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,37) = (yp(:,37)+dts*p)*tmp1
    !
    !          HOC2H4O2     y(38)
    p = 0.0                                                                    &
    +(rc(:,109)*y(:,3 )*y(:,49))
    l = 0.0                                                                    &
    +(rc(:,110)*y(:,5 ))+(rc(:,111)*y(:,16))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,38) = (yp(:,38)+dts*p)*tmp1
    !
    !          ORGNIT       y(39)
    p = 0.0                                                                    &
    +(rc(:,171)*y(:,7 )*y(:,42))+(rc(:,183)*y(:,7 )*y(:,55))
    l = 0.0                                                                    &
    +(dd(:,39) )                                                               &
    +(rc(:,177)*y(:,3 ))+(rc(:,178)*y(:,6 ))+(dw(:,39) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,39) = (yp(:,39)+dts*p)*tmp1
    !
    !          HOC3H6O2     y(40)
    p = 0.0                                                                    &
    +(rc(:,125)*y(:,3 )*y(:,50))
    l = 0.0                                                                    &
    +(rc(:,126)*y(:,5 ))+(rc(:,127)*y(:,16))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,40) = (yp(:,40)+dts*p)*tmp1
    !
    !          CH3OH        y(41)
    p = 0.0                                                                    &
    +(rc(:,62) *y(:,16)*y(:,16))+(rc(:,123)*y(:,4 )*y(:,50)*0.12)
    l = 0.0                                                                    &
    +(rc(:,63) *y(:,3 ))+(dw(:,41) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,41) = (yp(:,41)+dts*p)*tmp1
    !
    !          OXYL1        y(42)
    p = 0.0                                                                    &
    +(rc(:,176)*y(:,3 )*y(:,58))
    l = 0.0                                                                    &
    +(rc(:,171)*y(:,7 ))+(rc(:,180)*y(:,15))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,42) = (yp(:,42)+dts*p)*tmp1
    !
    !          H2           y(43)
    p = 0.0                                                                    &
    +(dj(:,7)  *y(:,14))                                                       &
    +(rc(:,112)*y(:,4 )*y(:,49)*0.13)+(rc(:,124)*y(:,4 )*y(:,50)*0.24)
    l = 0.0                                                                    &
    +(rc(:,33) *y(:,3 ))+(dd(:,43) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,43) = (yp(:,43)+dts*p)*tmp1
    !
    !          MEMALD1      y(44)
    p = 0.0                                                                    &
    +(rc(:,172)*y(:,3 )*y(:,56))
    l = 0.0                                                                    &
    +(rc(:,173)*y(:,5 ))+(rc(:,181)*y(:,16))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,44) = (yp(:,44)+dts*p)*tmp1
    !
    !          RNC2H4       y(45)
    p = 0.0                                                                    &
    +(rc(:,153)*y(:,6 )*y(:,49))
    l = 0.0                                                                    &
    +(rc(:,154)*y(:,3 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,45) = (yp(:,45)+dts*p)*tmp1
    !
    !          HOIPO2       y(46)
    p = 0.0                                                                    &
    +(rc(:,141)*y(:,3 )*y(:,29))
    l = 0.0                                                                    &
    +(rc(:,132)*y(:,16))+(rc(:,134)*y(:,15))+(rc(:,142)*y(:,5 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,46) = (yp(:,46)+dts*p)*tmp1
    !
    !          RNC3H6       y(47)
    p = 0.0                                                                    &
    +(rc(:,155)*y(:,6 )*y(:,50))
    l = 0.0                                                                    &
    +(rc(:,156)*y(:,3 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,47) = (yp(:,47)+dts*p)*tmp1
    !
    !          HOMVKO2      y(48)
    p = 0.0                                                                    &
    +(rc(:,143)*y(:,3 )*y(:,34))
    l = 0.0                                                                    &
    +(rc(:,133)*y(:,16))+(rc(:,136)*y(:,15))+(rc(:,144)*y(:,5 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,48) = (yp(:,48)+dts*p)*tmp1
    !
    !          C2H4         y(49)
    p = 0.0
    l = 0.0                                                                    &
    +(rc(:,109)*y(:,3 ))+(rc(:,112)*y(:,4 ))+(rc(:,153)*y(:,6 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,49) = (yp(:,49)+dts*p)*tmp1
    !
    !          C3H6         y(50)
    p = 0.0
    l = 0.0                                                                    &
    +(rc(:,155)*y(:,6 ))                                                       &
    +(rc(:,123)*y(:,4 ))+(rc(:,124)*y(:,4 ))+(rc(:,125)*y(:,3 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,50) = (yp(:,50)+dts*p)*tmp1
    !
    !          NC4H10       y(51)
    p = 0.0
    l = 0.0                                                                    &
    +(rc(:,81) *y(:,3 ))+(rc(:,152)*y(:,6 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,51) = (yp(:,51)+dts*p)*tmp1
    !
    !          s-BuOOH      y(52)
    p = 0.0                                                                    &
    +(rc(:,103)*y(:,15)*y(:,24))
    l = 0.0                                                                    &
    +(dw(:,52) )        +(dd(:,52) )                                           &
    +(rc(:,44) *y(:,3 ))+(rc(:,104)*y(:,3 ))+(dj(:,21) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,52) = (yp(:,52)+dts*p)*tmp1
    !
    !          MEK          y(53)
    p = 0.0                                                                    &
    +(rc(:,104)*y(:,3 )*y(:,52))+(dj(:,21) *y(:,52))                           &
    +(rc(:,83) *y(:,24)*y(:,5 ))+(rc(:,84) *y(:,24)*y(:,16))
    l = 0.0                                                                    &
    +(rc(:,86) *y(:,3 ))+(dj(:,9)  )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,53) = (yp(:,53)+dts*p)*tmp1
    !
    !          TOLUEN       y(54)
    p = 0.0
    l = 0.0                                                                    &
    +(rc(:,174)*y(:,3 ))+(rc(:,175)*y(:,3 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,54) = (yp(:,54)+dts*p)*tmp1
    !
    !          TOLP1        y(55)
    p = 0.0                                                                    &
    +(rc(:,175)*y(:,3 )*y(:,54))
    l = 0.0                                                                    &
    +(rc(:,182)*y(:,15))+(rc(:,183)*y(:,7 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,55) = (yp(:,55)+dts*p)*tmp1
    !
    !          MEMALD       y(56)
    p = 0.0                                                                    &
    +(rc(:,180)*y(:,42)*y(:,15))+(rc(:,182)*y(:,15)*y(:,55))                   &
    +(rc(:,177)*y(:,3 )*y(:,39))+(rc(:,178)*y(:,6 )*y(:,39))                   &
    +(rc(:,170)*y(:,3 )*y(:,58)*0.80)+(rc(:,174)*y(:,3 )*y(:,54))
    l = 0.0                                                                    &
    +(rc(:,172)*y(:,3 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,56) = (yp(:,56)+dts*p)*tmp1
    !
    !          GLY          y(57)
    p = 0.0                                                                    &
    +(rc(:,181)*y(:,44)*y(:,16))+(rc(:,182)*y(:,15)*y(:,55))                   &
    +(rc(:,177)*y(:,3 )*y(:,39))+(rc(:,178)*y(:,6 )*y(:,39))                   &
    +(rc(:,173)*y(:,44)*y(:,5 ))+(rc(:,174)*y(:,3 )*y(:,54))
    l = 0.0                                                                    &
    +(rc(:,139)*y(:,3 ))+(dj(:,13) )        +(dw(:,57) )

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,57) = (yp(:,57)+dts*p)*tmp1
    !
    !          OXYL         y(58)
    p = 0.0
    l = 0.0                                                                    &
    +(rc(:,170)*y(:,3 ))+(rc(:,176)*y(:,3 ))

    tmp2(:) = 1.0+dts*l(:)

    CALL oneover_v(n_pnts, tmp2, tmp1)

    y(:,58) = (yp(:,58)+dts*p)*tmp1
    !
  END DO ! End of iteration loop

  !       --------------------
  !       DRY DEPOSITION FLUXES
  !       --------------------
  !       Dry deposition fluxes in moles s-1.
  !       Note on units: s-1 * (molec * cm-3) * cm3 * (mol / molec) --> mol s-1
  !
  icnt_dd = 0
  DO j=1,jpspec
    IF (ldepd(j)) THEN
      icnt_dd = icnt_dd + 1
      DO i = 1, n_pnts
        ddflux(i,icnt_dd) = ddflux(i,icnt_dd) +                                &
                         dd(i,j)*y(i,j)*vol(i)                                 &
                           * recip_avogad
      END DO
    END IF
  END DO

  !         --------------------
  !         WET DEPOSITION FLUXES
  !         --------------------
  !         Wet deposition fluxes in moles s-1
  !
  icnt_wd = 0
  DO j=1,jpspec
    IF (ldepw(j)) THEN
      icnt_wd = icnt_wd + 1
      DO i = 1, n_pnts
        wdflux(i,icnt_wd) = wdflux(i,icnt_wd) +                                &
                         dw(i,j)*y(i,j)*vol(i)                                 &
                           * recip_avogad
      END DO
    END IF
  END DO

END DO  ! n_be_calls

! Dry & wet deposition fluxes have been added n_be_calls times (i.e.
! once per each chemical timestep). Need to divide by that number
! to get the average flux (in mol s-1) within the model timestep.
DO j = 1, icnt_dd
  ddflux (:,j) = ddflux (:,j) / n_be_calls
END DO

DO j = 1, icnt_wd
  wdflux (:,j) = wdflux (:,j) / n_be_calls
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_deriv_raq
END MODULE ukca_deriv_raq_mod
