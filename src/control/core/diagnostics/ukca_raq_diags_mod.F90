! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Produce diagnostics for RAQ chemistry.
!
!  Method:
!    Go through the list of diagnostics for RAQ chemistry.
!    If selected then copy them to stashwork via calls to COPYDIAG_3D.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------

MODULE ukca_raq_diags_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_RAQ_DIAGS_MOD'

CONTAINS

SUBROUTINE ukca_raq_diags (row_length, rows, model_levels, ntracers,           &
    dry_dep_3d, wet_dep_3d, tracer, pres, temp,                                &
    len_stashwork, stashwork)

USE asad_mod, ONLY: jpdd, jpdw

!     Indices for dry and wet deposition of species.
USE ukca_chem_raq,   ONLY: ndd_o3,   ndd_no2,    ndd_hono2,  ndd_h2o2,         &
    ndd_ch4,  ndd_co,    ndd_meooh,  ndd_etooh,  ndd_pan,    ndd_iprooh,       &
    ndd_o3s,  ndd_isooh, ndd_mvkooh, ndd_orgnit, ndd_h2,     ndd_sbuooh,       &
    ndw_no3,  ndw_n2o5,  ndw_ho2no2, ndw_hono2,  ndw_h2o2,   ndw_hcho,         &
    ndw_ho2,  ndw_meoo,  ndw_meooh,  ndw_etooh,  ndw_iprooh, ndw_isooh,        &
    ndw_ison, ndw_mgly,  ndw_mvkooh, ndw_orgnit, ndw_ch3oh,  ndw_sbuooh,       &
    ndw_gly

USE get_noy_mod,  ONLY: get_noy
USE get_nmvoc_mod,  ONLY: get_nmvoc
USE parkind1, ONLY: jprb,  jpim
USE yomhook,  ONLY: lhook, dr_hook

USE ukca_um_legacy_mod, ONLY: len_stlist, stindex, stlist, num_stash_levels,   &
                              stash_levels, si, sf, si_last, copydiag_3d

IMPLICIT NONE


! Subroutine arguments
INTEGER, INTENT(IN) :: row_length        ! Model dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: ntracers          ! no. of tracers

REAL, INTENT(IN) :: dry_dep_3d(row_length, rows, model_levels, jpdd)
                    ! dry deposition rates (mol s-1)
REAL, INTENT(IN) :: wet_dep_3d(row_length, rows, model_levels, jpdw)
                    ! wet deposition rates (mol s-1)
REAL, INTENT(IN) :: tracer    (row_length, rows, model_levels, ntracers)
                    ! tracers (MMR)
REAL, INTENT(IN) :: pres      (row_length, rows, model_levels)
                    ! Pressure (Pa)
REAL, INTENT(IN) :: temp      (row_length, rows, model_levels)
                    ! Temperature (K)

INTEGER, INTENT(IN) :: len_stashwork     ! Length of diagnostics array

REAL, INTENT(IN OUT) :: stashwork (len_stashwork) ! Diagnostics array

! Local variables
INTEGER    :: item         ! stash item
INTEGER    :: im_index     ! internal model index

INTEGER, PARAMETER :: sect = 50   ! stash section for UKCA diagnostics

REAL :: tmp_diags(row_length, rows, model_levels)

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RAQ_DIAGS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

im_index = 1

! ----------------------------------------------------------------------
! SECT. 50, ITEM 21: O3 DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 21

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_o3)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 79: CO DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 79

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_co)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 173: NO2 DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 173

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_no2)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 174: HONO2 DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 174

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_hono2)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 175: H2O2 DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 175

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_h2o2)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 176: CH4 DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 176

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_ch4)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 177: MeOOH DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 177

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_meooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 178: EtOOH DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 178

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_etooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 179: PAN DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 179

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_pan)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 180: i-PrOOH DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 180

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_iprooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 181: O3S DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 181

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_o3s)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 182: ISOOH DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 182

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_isooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 183: MVKOOH DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 183

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_mvkooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 184: ORGNIT DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 184

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_orgnit)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 185: H2 DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 185

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_h2)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 186: s_BuOOH DRY DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 186

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = dry_dep_3d(:,:,:,ndd_sbuooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 187: NO3 WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 187

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_no3)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 188: N2O5 WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 188

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_n2o5)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 189: HO2NO2 WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 189

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_ho2no2)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 190: HONO2 WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 190

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_hono2)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 191: H2O2 WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 191

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_h2o2)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 192: HCHO WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 192

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_hcho)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 193: HO2 WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 193

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_ho2)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 194: MeOO WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 194

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_meoo)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 195: MeOOH WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 195

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_meooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 196: EtOOH WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 196

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_etooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 197: i-PrOOH WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 197

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_iprooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 198: ISOOH WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 198

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_isooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 199: ISON WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 199

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_ison)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 200: MGLY WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 200

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_mgly)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 201: MVKOOH WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 201

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_mvkooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 202: ORGNIT WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 202

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_orgnit)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 203: CH3OH WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 203

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_ch3oh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 204: s-BuOOH WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 204

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_sbuooh)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 205: GLY WET DEPOSITION (3D)
! ----------------------------------------------------------------------
item = 205

IF (sf(item,sect)) THEN

  tmp_diags(:,:,:) = wet_dep_3d(:,:,:,ndw_gly)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 227: NOy DIAGNOSTIC (3D)
! ----------------------------------------------------------------------
item = 227

IF (sf(item,sect)) THEN

  ! Call routine that returns the array tmp_diags (:,:,:) holding the
  ! volume mixing ratios (mol mol-1) of total reactive nitrogen (NOy).
  CALL get_noy (row_length, rows, model_levels, ntracers, tracer, tmp_diags)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

! ----------------------------------------------------------------------
! SECT. 50, ITEM 230: Non-methane VOC (3D)
! ----------------------------------------------------------------------
item = 230

IF (sf(item,sect)) THEN

  ! Get the concentration (ug C m-3) of total non-methane volatile organic
  ! carbon, excluding PAN.

  CALL get_nmvoc (row_length, rows, model_levels, ntracers, tracer,            &
                  pres, temp, tmp_diags)

  CALL copydiag_3d (stashwork(si(item,sect,im_index):                          &
       si_last(item,sect,im_index)),                                           &
       tmp_diags(:,:,:),                                                       &
       row_length,rows,model_levels,                                           &
       stlist(:,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_raq_diags

END MODULE ukca_raq_diags_mod
