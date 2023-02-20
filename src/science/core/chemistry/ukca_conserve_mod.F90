! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************

!  Description:
!   Module to contain subroutine ukca_conserve

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################

MODULE ukca_conserve_mod

USE ukca_config_specification_mod, ONLY: ukca_config

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus,                         &
                      PrStatus_Normal, PrStatus_Oper
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CONSERVE_MOD'

CONTAINS

SUBROUTINE ukca_conserve(row_length, rows, model_levels,                       &
  ntracers, tracers,                                                           &
  pres, drain, crain, direction)
! Description:

! This routine calculates and conserves total chlorine, bromine, and
! hydrogen. For these elements closed chemistry should be prescribed.
! Called before chemistry, with direction = .TRUE., it calculates
! total bromine, chlorine, and hydrogen as 3-D fields. Called afer
! chemistry, with direction = .FALSE., it rescales the chlorine, bromine
! and hydrogen containing compounds so that total chlorine, bromine
! and hydrogen are conserved under chemistry. Where a compound contains
! more than one of the 3 elements. e.g, BrCl, it is only scaled for the
! less abundant of the two constituents, Br. It is then subtracted from
! the total chlorine.

! Method: Rescaling of tracer variables.


! Code Description:
! Language: FORTRAN 90 + common extensions.

! Declarations:
! These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable


USE asad_mod,           ONLY: advt, jpctr
USE ukca_constants,     ONLY: c_cf2cl2, c_cfcl3, c_clo, c_cl,                  &
  c_cl2o2, c_hcl, c_hocl, c_oclo, c_brcl,                                      &
  c_clono2, c_chf2cl, c_meccl3, c_mecl,                                        &
  c_cf2clbr, c_ccl4, c_ch2br2, c_cf2clcfcl2,                                   &
  c_bro, c_br, c_hobr, c_brono2, c_mebr, c_hbr,                                &
  c_cf3br, c_h2o, c_ch4, c_h2,                                                 &
  c_ho2no2, c_hono2, c_h2o2, c_hcho, c_meooh,                                  &
  c_hono, c_c2h6, c_etooh, c_mecho, c_pan,                                     &
  c_c3h8, c_prooh, c_etcho, c_mecoch2ooh,                                      &
  c_ppan, c_meono2, c_h, c_oh, c_ho2, c_meoo,                                  &
  c_etoo, c_proo, c_etco3, c_mecoch2oo

USE ereport_mod,        ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: row_length        ! no of points E-W
INTEGER, INTENT(IN) :: rows              ! no of points N-S
INTEGER, INTENT(IN) :: model_levels      ! no of levels
INTEGER, INTENT(IN) :: ntracers          ! no of tracers

LOGICAL, INTENT(IN) :: direction         ! T to calculate total Cl etc
! F to rescale total Cl into compouds

REAL, INTENT(IN) :: pres(row_length, rows, model_levels)  ! pressure (Pa)

REAL, INTENT(IN) :: drain(row_length, rows, model_levels)
REAL, INTENT(IN) :: crain(row_length, rows, model_levels)

REAL, INTENT(IN OUT) :: tracers(row_length, rows, model_levels,                &
  ntracers)    ! tracer array

! Local variables

! Maximum number of Br, Cl, and H compounds permitted
INTEGER, PARAMETER :: max_comp = 50

REAL, PARAMETER :: adjust_pres = 500.0 ! pressure below which allow
! for changes of total hydrogen due to dehydration, and above which only
! advective changes are allowed.

REAL, PARAMETER :: washout_limit = 10000.0 ! pressure limit below which
! hydrogen conservation is not enforced.

INTEGER :: m,i,j,k

! Reservoir tracer for bromine
INTEGER, SAVE :: n_toth=0      ! position of total hydrogen tracer
INTEGER, SAVE :: n_n2o=0         !             N2O

LOGICAL, SAVE :: firstcall = .TRUE. ! flag for first call of subr.

INTEGER, SAVE ::       ncl_tracers     ! number of Cl tracers
INTEGER, SAVE ::       nbr_tracers     ! number of Br tracers
INTEGER, SAVE ::       nh_tracers      ! number of H tracers

REAL, ALLOCATABLE, SAVE :: total_cl(:,:,:) ! total chlorine VMR
REAL, ALLOCATABLE, SAVE :: total_br(:,:,:) ! total bromine VMR
REAL, ALLOCATABLE, SAVE :: total_h (:,:,:) ! total hydrogen VMR


INTEGER, SAVE ::  cl_tracers(max_comp) ! positions of Cl tracers
INTEGER, SAVE ::  br_tracers(max_comp) ! positions of Br tracers
INTEGER, SAVE ::   h_tracers(max_comp) ! positions of hydrogen tracers

REAL   , SAVE ::c_cl_tracers(max_comp) ! conversion factors VMR/MMR
REAL   , SAVE ::c_br_tracers(max_comp) ! conversion factors
REAL   , SAVE ::c_h_tracers(max_comp)  ! conversion factors

REAL   , SAVE :: cl_validity(max_comp) ! number of Cl atoms per mol.
REAL   , SAVE :: br_validity(max_comp) ! number of Br atoms per mol.

REAL   , SAVE :: h_validity(max_comp)  ! number of H atoms per mol.

LOGICAL, SAVE :: contains_bromine(max_comp) ! flag for Cl compounds which
! also contain bromine.
LOGICAL, SAVE :: do_not_change(max_comp)    ! leave unchanged.

! Correction factor to achieve Cl or Br conservation
REAL, ALLOCATABLE :: corrfac(:,:,:)

INTEGER           :: errcode        ! error code
CHARACTER(LEN=errormessagelength) :: cmessage       !   "   message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CONSERVE'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (firstcall) THEN

  IF (ntracers /= jpctr) THEN
    errcode = 1
    cmessage=' Number of tracers does not agree with jpctr'
    WRITE(umMessage,'(1X,A,1X,I4,1X,I4)') cmessage,ntracers,jpctr
    CALL umPrint(umMessage,src='ukca_conserve_mod')
    CALL ereport('UKCA_CONSERVE_MOD:UKCA_CONSERVE',errcode,                    &
      cmessage)
  END IF

  ! Calculate the number of Cl and Br tracers, their positions, mmr/vmr
  ! conversion ratios, and validities (numbers of Cl/Br atoms per molecule).

  ncl_tracers = 0
  nbr_tracers = 0
  nh_tracers = 0
  contains_bromine = .FALSE.
  do_not_change = .FALSE.
  cl_tracers = 0
  br_tracers = 0
  h_tracers = 0
  cl_validity = 1.0
  br_validity = 1.0
  h_validity = 1.0
  c_br_tracers = 0.0
  c_cl_tracers = 0.0
  c_h_tracers = 0.0
  DO m=1,jpctr
    SELECT CASE (advt(m))

      ! chlorine tracers
    CASE ('CF2Cl2    ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cf2cl2
      cl_validity(ncl_tracers) = 2.0
    CASE ('CFCl3     ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cfcl3
      cl_validity(ncl_tracers) = 3.0
    CASE ('Clx       ','ClO       ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_clo
    CASE ('Cl        ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cl
    CASE ('Cl2O2     ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cl2o2
      cl_validity(ncl_tracers) = 2.0
    CASE ('HCl       ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_hcl
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hcl
      do_not_change(nh_tracers) = .TRUE.
    CASE ('HOCl      ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_hocl
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hocl
      do_not_change(nh_tracers) = .TRUE.
    CASE ('OClO      ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_oclo
    CASE ('BrCl      ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_brcl
      contains_bromine(ncl_tracers) = .TRUE.
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_brcl
    CASE ('ClONO2    ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_clono2
    CASE ('CF2ClCFCl2')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cf2clcfcl2
      cl_validity(ncl_tracers) = 3.0
    CASE ('CHF2Cl    ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_chf2cl
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_chf2cl
      do_not_change(nh_tracers) = .TRUE.
    CASE ('MeCCl3    ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_meccl3
      cl_validity(ncl_tracers) = 3.0
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_meccl3
      do_not_change(nh_tracers) = .TRUE.
    CASE ('CCl4      ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_ccl4
      cl_validity(ncl_tracers) = 4.0
    CASE ('MeCl      ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_mecl
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_mecl
      do_not_change(nh_tracers) = .TRUE.
    CASE ('CF2ClBr   ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cf2clbr
      contains_bromine(ncl_tracers) = .TRUE.
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_cf2clbr
      ! Bromine tracers
    CASE ('Brx       ','BrO       ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_bro
    CASE ('Br        ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_br
    CASE ('HOBr      ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_hobr
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hobr
      do_not_change(nh_tracers) = .TRUE.
    CASE ('BrONO2    ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_brono2
    CASE ('MeBr      ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_mebr
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_mebr
      do_not_change(nh_tracers) = .TRUE.
    CASE ('HBr       ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_hbr
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hbr
      do_not_change(nh_tracers) = .TRUE.
    CASE ('CF3Br     ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_cf3br
    CASE ('CH2Br2    ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_ch2br2
      br_validity(nbr_tracers) = 2.0
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_ch2br2
      h_validity(nh_tracers) = 2.0
      do_not_change(nh_tracers) = .TRUE.
      ! hydrogen tracers
    CASE ('TOTH      ')
      n_toth = m
    CASE ('H2O       ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_h2o
      h_validity(nh_tracers) = 2.0
    CASE ('CH4       ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_ch4
      h_validity(nh_tracers) = 4.0
      do_not_change(nh_tracers) = .TRUE.
    CASE ('H2        ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_h2
      h_validity(nh_tracers) = 2.0
    CASE ('HO2NO2    ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_ho2no2
      do_not_change(nh_tracers) = .TRUE.
    CASE ('HONO2     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hono2
      do_not_change(nh_tracers) = .TRUE.
    CASE ('H2O2      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_h2o2
      h_validity(nh_tracers) = 2.0
    CASE ('HCHO      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hcho
      h_validity(nh_tracers) = 2.0
    CASE ('MeOOH     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_meooh
      h_validity(nh_tracers) = 4.0
    CASE ('HONO      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hono
      do_not_change(nh_tracers) = .TRUE.
    CASE ('C2H6      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_c2h6
      h_validity(nh_tracers) = 6.0
    CASE ('EtOOH     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_etooh
      h_validity(nh_tracers) = 6.0
    CASE ('MeCHO     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_mecho
      h_validity(nh_tracers) = 4.0
    CASE ('PAN       ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_pan
      h_validity(nh_tracers) = 3.0
      do_not_change(nh_tracers) = .TRUE.
    CASE ('C3H8      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_c3h8
      h_validity(nh_tracers) = 8.0
    CASE ('n-PrOOH   ','i-PrOOH   ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_prooh
      h_validity(nh_tracers) = 8.0
    CASE ('EtCHO     ','Me2CO     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_etcho
      h_validity(nh_tracers) = 6.0
    CASE ('MeCOCH2OOH')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_mecoch2ooh
      h_validity(nh_tracers) = 6.0
    CASE ('PPAN      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_ppan
      h_validity(nh_tracers) = 5.0
      do_not_change(nh_tracers) = .TRUE.
    CASE ('MeONO2    ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_meono2
      h_validity(nh_tracers) = 3.0
      do_not_change(nh_tracers) = .TRUE.
    CASE ('H         ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_h
    CASE ('OH        ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_oh
    CASE ('HO2       ','HOx       ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_ho2
    CASE ('MeOO      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_meoo
      h_validity(nh_tracers) = 3.0
    CASE ('EtOO      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_etoo
      h_validity(nh_tracers) = 5.0
    CASE ('i-PrOO    ','n-PrOO    ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_proo
      h_validity(nh_tracers) = 7.0
    CASE ('EtCO3     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_etco3
      h_validity(nh_tracers) = 5.0
    CASE ('MeCOCH2OO ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_mecoch2oo
      h_validity(nh_tracers) = 5.0
      ! special tracers
    CASE ('Ox        ','O3        ')
    CASE ('N2O       ')
      n_n2o = m
    CASE DEFAULT
      IF ( PrintStatus >= PrStatus_Oper ) THEN
        cmessage = ' CONSERVE: SPECIES '//advt(m)//                            &
        ' not treated in CASE'
        errcode = -1*m
        CALL ereport('UKCA_CONSERVE_MOD:UKCA_CONSERVE',errcode,cmessage)
      END IF
    END SELECT
  END DO

  firstcall = .FALSE.

END IF  ! if firstcall

IF (direction) THEN
  ! Calculate total chlorine and bromine vmrs.
  IF (.NOT. ALLOCATED(total_br))                                               &
    ALLOCATE(total_br(row_length, rows, model_levels))
  IF (.NOT. ALLOCATED(total_cl))                                               &
    ALLOCATE(total_cl(row_length, rows, model_levels))

  total_br = 0.0
  total_cl = 0.0

  DO m=1,nbr_tracers
    total_br = total_br + tracers(:,:,:,br_tracers(m))*                        &
      br_validity(m) / c_br_tracers(m)
  END DO
  DO m=1,ncl_tracers
    total_cl = total_cl + tracers(:,:,:,cl_tracers(m))*                        &
      cl_validity(m) / c_cl_tracers(m)
  END DO

  ! Do not do hydrogen conservation unless at least two hydrogen
  ! reservoirs (H2O, CH4) are defined, and water vapour feedback is on.
  ! H conservation is currently experimental, so is normally set to be off.

  IF (ukca_config%l_ukca_h2o_feedback .AND. ukca_config%l_ukca_conserve_h) THEN
    ALLOCATE(total_h(row_length, rows, model_levels))
    total_h = 0.0
    DO m=1,nh_tracers
      total_h = total_h + tracers(:,:,:,h_tracers(m))*                         &
        h_validity(m) / c_h_tracers(m)
    END DO
    IF (n_toth > 0) THEN
      ! if a seperate total hydrogen tracer is present
      WHERE (pres > adjust_pres)                                               &
        tracers(:,:,:,n_toth) = total_h
      WHERE (pres <= adjust_pres)                                              &
        total_h = tracers(:,:,:,n_toth)
    END IF
  END IF

ELSE    ! direction = F
  ! Set N2O to 0 at the top, to avoid initialization problem
  ! ensure to not do this in case of box model (model_levels==1)
  IF (n_n2o > 0 .AND. (model_levels /= 1))                                     &
    tracers(:,:,model_levels,n_n2o) = 0.0

  ! Adjust tracers to match
  ALLOCATE(corrfac(row_length, rows, model_levels))
  IF (nbr_tracers > 0) THEN
    corrfac = 0.0

    ! Calculate new total bromine
    DO m=1,nbr_tracers
      corrfac = corrfac + tracers(:,:,:,br_tracers(m))*                        &
        br_validity(m) / c_br_tracers(m)
    END DO

    ! Adjust bromine tracers to match total bromine computed before
    ! chemistry

    corrfac = total_br / corrfac
    DO m=1,nbr_tracers
      DO i=1,row_length
        DO j=1,rows
          DO k=2,model_levels
            IF (drain(i,j,k) + crain(i,j,k) == 0.0) THEN
              tracers(i,j,k,br_tracers(m)) =                                   &
                tracers(i,j,k,br_tracers(m)) * corrfac(i,j,k)
            ELSE
              corrfac(i,j,k) = 1.0
            END IF
          END DO
        END DO
      END DO
    END DO
    IF (((MINVAL(corrfac) < 0.9) .OR. (MAXVAL(corrfac) > 1.1))                 &
      .AND. (printstatus >= prstatus_normal)) THEN
      WRITE(umMessage,'(A,2F8.4)') 'Correct bromine ',MINVAL(corrfac),         &
                                                      MAXVAL(corrfac)
      CALL umPrint(umMessage,src='ukca_conserve_mod')
    END IF
  END IF     ! nbr_tracers > 0

  ! Calculate new total chlorine, excluding BrCl and CF2ClBr

  IF (ncl_tracers > 0) THEN
    corrfac = 0.0

    DO m=1,ncl_tracers
      IF (.NOT. (contains_bromine(m))) THEN
        corrfac = corrfac + tracers(:,:,:,cl_tracers(m))*                      &
          cl_validity(m) / c_cl_tracers(m)
      ELSE
        total_cl = total_cl - tracers(:,:,:,cl_tracers(m))*                    &
          cl_validity(m) / c_cl_tracers(m)
      END IF
    END DO

    ! Adjust chlorine species to match total chlorine computed before.
    ! Leave BrCl and CF2ClBr alone.

    corrfac = total_cl / corrfac
    DO m=1,ncl_tracers
      IF (.NOT. (contains_bromine(m))) THEN
        DO i=1,row_length
          DO j=1,rows
            DO k=2,model_levels
              IF (crain(i,j,k) + drain(i,j,k) == 0.0) THEN
                tracers(i,j,k,cl_tracers(m)) =                                 &
                  tracers(i,j,k,cl_tracers(m)) * corrfac(i,j,k)
              ELSE
                corrfac(i,j,k) = 1.0
              END IF
            END DO
          END DO
        END DO
      END IF
    END DO
    IF (((MINVAL(corrfac) < 0.9) .OR. (MAXVAL(corrfac) > 1.1))                 &
      .AND. (printstatus >= prstatus_normal)) THEN
      WRITE(umMessage,'(A,2F8.4)')'Correct chlorine ',MINVAL(corrfac),         &
                                                      MAXVAL(corrfac)
      CALL umPrint(umMessage,src='ukca_conserve_mod')
    END IF
  END IF   ! ncl_tracers > 0

  IF (ukca_config%l_ukca_h2o_feedback .AND. ukca_config%l_ukca_conserve_h) THEN
    corrfac = 0.0

    ! Calculate new total hydrogen
    DO m=1,nh_tracers
      IF (.NOT. (do_not_change(m))) THEN
        corrfac = corrfac + tracers(:,:,:,h_tracers(m))*                       &
          h_validity(m) / c_h_tracers(m)
      ELSE
        total_h = total_h - tracers(:,:,:,h_tracers(m))*                       &
          h_validity(m) / c_h_tracers(m)
      END IF
    END DO

    ! adjust upper boundary for hydrogen tracers only. This is needed to
    ! prevent model instability if too much water is present at the model
    ! top.
    ! Adjust hydrogen tracers to match total hydrogen computed before
    ! chemistry

    corrfac = total_h / corrfac

    ! Do not enforce hydrogen conservation in the troposphere, due to
    ! washout and dry deposition.

    WHERE (pres > washout_limit) corrfac = 1.0

    DO m=1,nh_tracers
      IF (.NOT. (do_not_change(m)))                                            &
        tracers(:,:,:,h_tracers(m)) = tracers(:,:,:,h_tracers(m))              &
        * corrfac
    END DO
  END IF ! l_ukca_h2o_feedback

  ! Deallocate fields
  IF (ALLOCATED(total_br)) DEALLOCATE(total_br)
  IF (ALLOCATED(total_cl)) DEALLOCATE(total_cl)
  IF (ALLOCATED(total_h)) DEALLOCATE(total_h)
  IF (ALLOCATED(corrfac)) DEALLOCATE(corrfac)
END IF   ! direction

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_conserve

END MODULE ukca_conserve_mod
