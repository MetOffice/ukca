! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Purpose: Heterogeneous chemistry routine (includes aqueous phase reactions).
!
!     The purpose of this routine is to set and return the heterogeneous
!     reaction rates. If the user has heterogeneous chemistry turned on
!     then this subroutine will be called. The user must supply their
!     own version of this routine to compute the heterogeneous rates.
!
!     Note that this subroutine is called repeatedly. It should not
!     therefore be used to do any I/O unless absolutely necessary. The
!     routine inihet is provided to initialise the heterogeneous chemist
!     by reading in files etc.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_hetero_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_HETERO_MOD'

CONTAINS

SUBROUTINE asad_hetero(n_points, cld_f, cld_l, rc_het, H_plus_1d_arr)

USE asad_findreaction_mod, ONLY: asad_findreaction
USE asad_mod,        ONLY: t, p, tnd, rk, ih_o3, ih_h2o2, ih_so2,              &
                           ih_hno3, ihso3_h2o2, iho2_h, in2o5_h,               &
                           iso3_o3, ihso3_o3, ih2o2_oh,                        &
                           ihno3_oh, spb, sph, nbrkx, nhrkx,                   &
                           jpspb, jpsph, jpeq, ih_o3_const,                    &
                           jpbk, jphk, jpdw
USE ukca_config_specification_mod, ONLY: ukca_config
USE ukca_chem_offline, ONLY: nwet_constant
USE ukca_fdiss_constant_mod, ONLY: ukca_fdiss_constant
USE chemistry_constants_mod,  ONLY: avogadro
USE water_constants_mod,  ONLY: rho_water
USE ukca_constants,  ONLY: m_air, H_plus
USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook
USE ereport_mod,     ONLY: ereport
USE umPrintMgr, ONLY:                                                          &
    umPrint,                                                                   &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

USE ukca_fdiss_mod, ONLY: ukca_fdiss
IMPLICIT NONE


INTEGER, INTENT(IN) :: n_points         ! No of spatial points

REAL, INTENT(IN) :: cld_f(n_points)     ! Cloud fraction
REAL, INTENT(IN) :: cld_l(n_points)     ! Cloud liquid water (kg/kg)
REAL, INTENT(IN) :: rc_het(n_points,2)  ! Heterog. Chem. Rates (tropospheric)
! 1-D pH array to be used in asad_hetero
REAL, INTENT(IN) :: H_plus_1d_arr(n_points)

!     Local variables

REAL, PARAMETER    :: qcl_min = 1.0e-12 ! do calcs when qcl > qcl_min
REAL               :: vr(n_points)      ! volume ratio
REAL               :: fdiss(n_points, jpdw, jpeq+1)
                                        ! fractional dissociation array
                                        ! final index: 1) dissolved
                                        !              2) 1st dissociation
                                        !              3) 2nd dissociation
REAL, ALLOCATABLE  :: fdiss_constant(:,:,:)
                                    ! As fdiss, but for constant species
REAL               :: fdiss_o3(n_points) ! fractional dissociation for O3

INTEGER            :: icode = 0         ! Error code
CHARACTER (LEN=errormessagelength) :: cmessage          ! Error message
CHARACTER(LEN=10)  :: prods(2)          ! Products
LOGICAL, SAVE      :: first = .TRUE.    ! Identifies first call
LOGICAL, SAVE      :: first_pass = .TRUE.    ! Identifies if thread has
                                             ! been through CRITICAL region
LOGICAL            :: todo(n_points)    ! T where cloud frac above threshold

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_HETERO'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!       1. Identify equations and calculate fractional dissociation
!          --------------------------------------------------------

IF (first_pass) THEN
  ! OMP CRITICAL will only allow one thread through this code at a time,
  ! while the other threads are held until completion.
!$OMP CRITICAL (asad_hetero_init)
  IF (first) THEN

    IF (ukca_config%l_ukca_achem) THEN
      ! Check that the indicies of the aqueous arrays are identified
      IF (ih_o3 == 0 .OR. ih_h2o2 == 0 .OR. ih_so2 == 0 .OR.                   &
          ih_hno3 == 0 ) THEN
        cmessage=' Indicies for Aqueous chemistry uninitialised'//             &
                 ' - O3, H2O2, SO2, and HNO3 must be made '//                  &
                 ' soluble species in chch_defs array'
        WRITE(umMessage,'(A9,I5)') 'ih_o3:   ',ih_o3
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A9,I5)') 'ih_h2o2: ',ih_h2o2
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A9,I5)') 'ih_hno3: ',ih_hno3
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A9,I5)') 'ih_so2:  ',ih_so2
        CALL umPrint(umMessage,src='asad_hetero')
        icode=1
        CALL ereport('ASAD_HETERO',icode,cmessage)
      END IF
    END IF

    IF (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) THEN
      ! Check that the indices of the aqueous arrays are identified
      IF (ih_h2o2 == 0 .OR. ih_so2 == 0 ) THEN
        cmessage=' Indices for Aqueous chemistry uninitialised'//              &
                 ' - H2O2, and SO2, must be made '//                           &
                 ' soluble species in chch_defs array'
        icode=1
      END IF

      IF (icode /= 0) THEN
        WRITE(umMessage,'(A10,I5)') 'ih_h2o2: ',ih_h2o2
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A10,I5)') 'ih_so2:  ',ih_so2
        CALL umPrint(umMessage,src='asad_hetero')
        CALL ereport('ASAD_HETERO',icode,cmessage)
      END IF
    END IF

    IF (ukca_config%l_ukca_trophet .AND. .NOT. ukca_config%l_ukca_mode) THEN
      cmessage=' Tropospheric heterogeneous chemistry is flagged'//            &
             ' but MODE aerosol scheme is not in use'
      icode=1
      CALL ereport('ASAD_HETERO',icode,cmessage)
    END IF

    ! Find reaction locations
    ihso3_h2o2 = 0
    iso3_o3    = 0
    ihso3_o3   = 0
    ih2o2_oh   = 0
    ihno3_oh   = 0
    in2o5_h    = 0
    iho2_h     = 0


    IF (ukca_config%l_ukca_nr_aqchem .OR. ukca_config%l_ukca_offline_be) THEN

      prods = ['NULL0     ','          ']
      ihso3_h2o2 = asad_findreaction( 'SO2       ', 'H2O2      ',              &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )
      prods = ['NULL1     ','          ']   ! Identifies HSO3- + O3(aq)
      ihso3_o3 = asad_findreaction( 'SO2       ', 'O3        ',                &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )
      prods = ['NULL2     ','          ']   ! Identifies SO3-- + O3(aq)
      iso3_o3 = asad_findreaction( 'SO2       ', 'O3        ',                 &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )

      IF (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) THEN
        prods = ['H2O       ','          ']
      ELSE
        prods = ['H2O       ','HO2       ']
      END IF
      ih2o2_oh = asad_findreaction( 'H2O2      ', 'OH        ',                &
                               prods, 2, spb, nbrkx, jpbk+1, jpspb )
      prods = ['H2O       ','NO3       ']
      ihno3_oh = asad_findreaction( 'HONO2     ', 'OH        ',                &
                               prods, 2, spb, nbrkx, jpbk+1, jpspb )

      icode = 0
      IF (ihso3_h2o2 == 0 .OR. iso3_o3 == 0 .OR. ih2o2_oh == 0 .OR.            &
          ihso3_o3 == 0 ) THEN
        WRITE(umMessage,'(A12,I5)') 'ihso3_h2o2: ',ihso3_h2o2
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A12,I5)') 'ihso3_o3: ',ihso3_o3
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A12,I5)') 'iso3_o3: ',iso3_o3
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A12,I5)') 'ih2o2_oh: ',ih2o2_oh
        CALL umPrint(umMessage,src='asad_hetero')
        icode = 1
      END IF
      IF (ukca_config%l_ukca_achem .AND. ihno3_oh == 0) THEN
        icode = 1
        WRITE(umMessage,'(A12,I5)') 'ihno3_oh: ',ihno3_oh
        CALL umPrint(umMessage,src='asad_hetero')
      END IF
      IF (icode > 0) THEN
        cmessage=' Heterogeneous chemistry called, but eqns'//                 &
                  ' not found - see output'
        CALL ereport('ASAD_HETERO',icode,cmessage)
      END IF
    END IF   ! l_ukca_achem.....

    ! Search for tropospheric heterogeneous reactions
    IF (ukca_config%l_ukca_trophet) THEN
      prods = ['HONO2     ','          ']
      in2o5_h = asad_findreaction( 'N2O5      ', '          ',                 &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )
      prods = ['H2O2      ','          ']
      iho2_h = asad_findreaction( 'HO2       ', '          ',                  &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )

      IF (iho2_h == 0 .OR. in2o5_h == 0) THEN
        WRITE(umMessage,'(A9,I5)') 'in2o5_h: ',in2o5_h
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A9,I5)') 'iho2_h: ',iho2_h
        CALL umPrint(umMessage,src='asad_hetero')
        cmessage=' Tropospheric heterogeneous chemistry is flagged,'//         &
                 ' but equations not found - see output'
        icode = 1
        CALL ereport('ASAD_HETERO',icode,cmessage)
      END IF   ! iho3_h=0 etc

    END IF      ! l_ukca_trophet

    first = .FALSE.

  END IF      ! first
!$OMP END CRITICAL (asad_hetero_init)
END IF        ! first_pass

! calculate fraction of dissolved species for online NR and offline BE
IF ((ukca_config%l_ukca_nr_aqchem .OR. ukca_config%l_ukca_offline_be) .AND.    &
    ANY(cld_l > qcl_min)) THEN
  ! send new H_plus array to calculate fraction dissolved in solvers
  CALL ukca_fdiss(n_points, qcl_min, t, p, cld_l, fdiss, H_plus_1d_arr)
  todo(:) = cld_l(:) > qcl_min
ELSE
  fdiss(:,:,:) = 0.0
  todo(:) = .FALSE.
END IF

! Assign fraction of O3 dissolved
IF (ANY(cld_l > qcl_min)) THEN
  IF ((ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) .AND.    &
      nwet_constant > 0 ) THEN
    ALLOCATE(fdiss_constant(n_points, nwet_constant, jpeq+1))
    ! send H_plus array to calculate fraction dissolved in offline oxidants
    CALL ukca_fdiss_constant(n_points, qcl_min, t, p, cld_l,                   &
                           fdiss_constant, H_plus_1d_arr)
    fdiss_o3(:) = fdiss_constant(:,ih_o3_const,1)
    DEALLOCATE(fdiss_constant)
  ELSE
    fdiss_o3(:) = fdiss(:,ih_o3,1)
  END IF
END IF

!    2. Calculate heterogeneous rates and reduce rates due to aqueous fraction
!       ----------------------------------------------------------------------

IF (ukca_config%l_ukca_nr_aqchem .OR. ukca_config%l_ukca_offline_be) THEN

  ! use new H_plus array in calculating reaction rates
  WHERE (todo(:))

    ! Convert clw in kg/kg to volume ratio
    vr(:) = cld_l(:)*tnd(:)*m_air*(1e6/avogadro)/rho_water

    ! HSO3- + H2O2(aq) => SO4--  [Kreidenweis et al. (2003), optimised]
    ! optimised means incorporated the EXP[(E/R)/298] part of the rate
    ! expression with K298 (7.45E7) from Kreidenweis (2003)
    rk(:,ihso3_h2o2) = 2.1295e+14*EXP(-4430.0/t(:))*                           &
       (H_plus_1d_arr(:)/(1.0 + 13.0*H_plus_1d_arr(:)))                        &
       *cld_f(:)*fdiss(:,ih_so2,2)*fdiss(:,ih_h2o2,1)*1000.0/(avogadro*vr(:))

    ! HSO3- + O3(aq) => SO4--  [Kreidenweis et al. (2003), optimised]
    ! optimised means incorporated the EXP[(E/R)/298] part of the rate
    ! expression with K298 (3.5E5) from Kreidenweis (2003)
    rk(:,ihso3_o3) = 4.0113e+13*EXP(-5530.0/t(:))*                             &
                     cld_f(:)*fdiss(:,ih_so2,2)*fdiss_o3(:)*                   &
                     1000.0/(avogadro*vr(:))

    ! SO3-- + O3(aq) => SO4-- [Kreidenweis et al. (2003), optimised]
    ! optimised means incorporated the EXP[(E/R)/298] part of the rate
    ! expression with K298 (1.5E9) from Kreidenweis (2003)
    rk(:,iso3_o3) = 7.43e+16*EXP(-5280.0/t(:))*cld_f(:)*                       &
                       fdiss(:,ih_so2,3)*fdiss_o3(:)*                          &
                       1000.0/(avogadro*vr(:))

    ! H2O2 + OH: reduce to take account of dissolved fraction
    rk(:,ih2o2_oh) = rk(:,ih2o2_oh)*                                           &
          (1.0 - (fdiss(:,ih_h2o2,1)+fdiss(:,ih_h2o2,2))*cld_f(:))
  ELSE WHERE
    rk(:,ihso3_h2o2) = 0.0
    rk(:,ihso3_o3) = 0.0
    rk(:,iso3_o3) = 0.0
  END WHERE

END IF      ! l_ukca_achem .....

IF (ukca_config%l_ukca_achem) THEN
  !  HNO3 + OH : reduce to take account of dissolved fraction
  WHERE (todo(:))
    rk(:,ihno3_oh) = rk(:,ihno3_oh)*                                           &
        (1.0 - (fdiss(:,ih_hno3,1)+fdiss(:,ih_hno3,2))*cld_f(:))
  END WHERE
END IF

IF (ukca_config%l_ukca_trophet) THEN
  ! N2O5 => HNO3 (heterogeneous)
  rk(:,in2o5_h) = rc_het(:,1)

  ! HO2 + HO2 => H2O2 (heterogeneous)
  rk(:,iho2_h) = rc_het(:,2)
ELSE
  IF (in2o5_h > 0) rk(:,in2o5_h) = 0.0
  IF (iho2_h > 0)  rk(:,iho2_h) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_hetero
END MODULE asad_hetero_mod
