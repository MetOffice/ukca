! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Combine surface resistance rc with aerodynamic resistance ra and
!  quasi-laminar resistance rb to get overall dry deposition velocity.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_DDEPCTL.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_ddcalc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_DDCALC_MOD'

CONTAINS

SUBROUTINE ukca_ddcalc(row_length, rows, bl_levels, ntype,                     &
                       npft, timestep, dzl, zbl, gsf, ra, rb, rc,              &
                       nlev_with_ddep, zdryrt, len_stashwork, stashwork)


USE asad_mod,                ONLY:                                             &
    ndepd,                                                                     &
    nldepd,                                                                    &
    speci,                                                                     &
    jpdd

USE parkind1,                ONLY:                                             &
    jprb,                                                                      &
    jpim

USE ukca_config_specification_mod, ONLY:                                       &
    ukca_config

USE ukca_um_legacy_mod,      ONLY:                                             &
    len_stlist, stash_pseudo_levels, num_stash_pseudo, stindex, stlist,        &
    si, si_last, sf, pdims, copydiag, set_pseudo_list

USE yomhook,                 ONLY:                                             &
    lhook,                                                                     &
    dr_hook


IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: bl_levels
INTEGER, INTENT(IN) :: ntype
INTEGER, INTENT(IN) :: npft

REAL, INTENT(IN) :: timestep
REAL, INTENT(IN) :: zbl(row_length,rows)             ! boundary layer depth
REAL, INTENT(IN) :: dzl(row_length, rows, bl_levels) ! thickness of BL levels
REAL, INTENT(IN) :: gsf(row_length,rows,ntype)       ! surface heat flux
REAL, INTENT(IN) :: ra(row_length,rows,ntype)        ! aerodynamic resistance
REAL, INTENT(IN) :: rb(row_length,rows,jpdd)         ! quasi-laminar resistance
REAL, INTENT(IN) :: rc(row_length,rows,ntype,jpdd)   ! surface resistance

! Diagnostics array
INTEGER, INTENT(IN)  :: len_stashwork
REAL, INTENT(IN OUT) :: stashwork (len_stashwork)

! no of levels over which dry deposition acts
INTEGER, INTENT(OUT) :: nlev_with_ddep(row_length,rows)
REAL, INTENT(OUT)    :: zdryrt(row_length,rows,jpdd) ! dry deposition rate

INTEGER :: i, j, loop, n
INTEGER :: pslevel   !  loop counter for pseudolevels
INTEGER :: pslevel_out   !  index for pseudolevels sent to STASH
INTEGER :: si_start, si_stop   !Stashwork bounds for calling copydiag

LOGICAL :: plltile(ntype)   ! pseudolevel list for surface types

REAL :: dd
REAL :: layer_depth(row_length,rows)
REAL :: vd(row_length,rows,ntype,jpdd)   ! deposition velocity

REAL :: r_nodep = 1.0e40

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DDCALC'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!     Set all arrays to zero
DO loop = 1, jpdd
  DO j = 1, rows
    DO i = 1, row_length
      zdryrt(i,j,loop) = 0.0
    END DO
  END DO
END DO

DO loop = 1, jpdd
  DO n = 1, ntype
    DO j = 1, rows
      DO i = 1, row_length
        vd(i,j,n,loop) = 0.0
      END DO
    END DO
  END DO
END DO


! Two options:
!
! (A) dry deposition confined to lowest model layer:
IF (ukca_config%l_ukca_ddep_lev1) THEN
  nlev_with_ddep(:,:) = 1
  layer_depth(:,:) = dzl(:,:,1)
ELSE

  ! (B) If dry deposition applied everywhere within the boundary
  !     layer then look for the highest model level completely
  !     contained in it.
  DO j = 1, rows
    DO i = 1, row_length
      layer_depth(i,j) = dzl(i,j,1)
      nlev_with_ddep(i,j) = 1
    END DO
  END DO

  DO loop = 2, bl_levels
    DO j = 1, rows
      DO i = 1, row_length
        dd = layer_depth(i,j) + dzl(i,j,loop)
        IF (dd < zbl(i,j)) THEN
          layer_depth(i,j) = dd
          nlev_with_ddep(i,j) = loop
        END IF
      END DO
    END DO
  END DO

END IF

!     Calculate overall dry deposition velocity [vd = 1/(ra + rb + rc)]
!     Do vegetated tiles first. Quasi-laminar resistance pre-multiplied by
!     ln[z0m/z0] = 2.0 for vegetated areas, or 1.0 for smooth surfaces
!     See Ganzeveld & Lelieveld, JGR 1995 Vol.100 No. D10 pp.20999-21012.

DO loop = 1, ndepd
  DO n = 1, npft
    DO j = 1, rows
      DO i = 1, row_length
        IF (rc(i,j,n,loop) < r_nodep .AND. gsf(i,j,n) > 0.0) THEN
          vd(i,j,n,loop) = 1.0 /                                               &
            (ra(i,j,n) + 2.0*rb(i,j,loop) + rc(i,j,n,loop))
        END IF
      END DO
    END DO
  END DO
END DO

    !         Now do calculation for non-vegetated tiles

DO loop = 1, ndepd
  DO n = npft+1, ntype
    DO j = 1, rows
      DO i = 1, row_length
        IF (rc(i,j,n,loop) < r_nodep .AND. gsf(i,j,n) > 0.0) THEN
          vd(i,j,n,loop) = 1.0 /                                               &
            (ra(i,j,n) + rb(i,j,loop) + rc(i,j,n,loop))
        END IF
      END DO
    END DO
  END DO
END DO

!     Now write dry deposition velocities for each tile
!     to diagnostic stream into STASHitem m1s50i434

DO loop = 1, ndepd
  ! process vd diagnostics for CH4 (m1s50i432), O3 (m1s50i433), and
  ! HONO2 (m1s50i434). These diagnostics are mainly for the purpose of
  ! model debugging.
  IF (speci(nldepd(loop)) == 'CH4       ') THEN
    ! vd on tiles for CH4
    IF (sf(432,50)) THEN
      CALL set_pseudo_list(ntype,len_stlist,                                   &
        stlist(:,stindex(1,432,50,1)),                                         &
        plltile,stash_pseudo_levels,num_stash_pseudo)
      pslevel_out=0
      DO pslevel = 1, ntype
        IF (plltile(pslevel)) THEN
          pslevel_out=pslevel_out+1
          si_start = si(432,50,1)+(pslevel_out-1)*pdims%i_end*pdims%j_end
          si_stop  = si(432,50,1)+(pslevel_out)*pdims%i_end*pdims%j_end-1
          CALL copydiag(stashwork(si_start:si_stop),vd(:,:,pslevel,loop),      &
            pdims%i_end,pdims%j_end)
        END IF
      END DO
    END IF
  ELSE IF (speci(nldepd(loop)) == 'O3        ') THEN
    ! vd on tiles for O3
    IF (sf(433,50)) THEN
      CALL set_pseudo_list(ntype,len_stlist,                                   &
        stlist(:,stindex(1,433,50,1)),                                         &
        plltile,stash_pseudo_levels,num_stash_pseudo)
      pslevel_out=0
      DO pslevel = 1, ntype
        IF (plltile(pslevel)) THEN
          pslevel_out=pslevel_out+1
          si_start = si(433,50,1)+(pslevel_out-1)*pdims%i_end*pdims%j_end
          si_stop  = si(433,50,1)+(pslevel_out)*pdims%i_end*pdims%j_end-1
          CALL copydiag(stashwork(si_start:si_stop),vd(:,:,pslevel,loop),      &
            pdims%i_end,pdims%j_end)
        END IF
      END DO
    END IF
  ELSE IF (speci(nldepd(loop)) == 'HONO2     ') THEN
    ! vd on tiles for HONO2
    IF (sf(434,50)) THEN
      CALL set_pseudo_list(ntype,len_stlist,                                   &
        stlist(:,stindex(1,434,50,1)),                                         &
        plltile,stash_pseudo_levels,num_stash_pseudo)
      pslevel_out=0
      DO pslevel = 1, ntype
        IF (plltile(pslevel)) THEN
          pslevel_out=pslevel_out+1
          si_start = si(434,50,1)+(pslevel_out-1)*pdims%i_end*pdims%j_end
          si_stop  = si(434,50,1)+(pslevel_out)*pdims%i_end*pdims%j_end-1
          CALL copydiag(stashwork(si_start:si_stop),vd(:,:,pslevel,loop),      &
            pdims%i_end,pdims%j_end)
        END IF
      END DO
    END IF
  END IF
END DO

!     VD() now contains dry deposition velocities for each tile
!     in each grid sq. Calculate overall first-order loss rate
!     over time "timestep" for each tile and sum over all tiles
!     to obtain overall first-order loss rate zdryrt().

DO loop = 1, ndepd
  DO n = 1, ntype
    DO j = 1, rows
      DO i = 1, row_length
        IF (vd(i,j,n,loop) > 0.0) THEN
          zdryrt(i,j,loop) = zdryrt(i,j,loop) + gsf(i,j,n) *                   &
            (1.0-EXP(-vd(i,j,n,loop) * timestep / layer_depth(i,j)))
        END IF
      END DO
    END DO
  END DO
END DO

!     ZDRYRT() contains loss rate over time "timestep".
!     Divide by timestep to get rate in s-1.

DO loop = 1, ndepd
  DO j = 1, rows
    DO i = 1, row_length
      zdryrt(i,j,loop) = -LOG(1.0 - zdryrt(i,j,loop)) / timestep
    END DO
  END DO
END DO

!     Now write 2D dry deposition rate (1/s) for entire gridbox
!     for a number of selected species to diagnostic stream into
!     STASHitem m1s50i435-37
!     This is a special diagnostic, mainly for debugging

DO loop = 1, ndepd
  ! process ddep diagnostics for CH4
  IF (speci(nldepd(loop)) == 'CH4       ') THEN
    ! zdryrt on gridbox for CH4
    IF (sf(435,50)) THEN
      CALL copydiag (                                                          &
        stashwork(si(435,50,1):si_last(435,50,1)),                             &
        zdryrt(:,:,loop),                                                      &
        row_length,rows)
    END IF
    ! process ddep diagnostics for O3
  ELSE IF (speci(nldepd(loop)) == 'O3        ') THEN
    ! zdryrt on gridbox for O3
    IF (sf(436,50)) THEN
      CALL copydiag (                                                          &
        stashwork(si(436,50,1):si_last(436,50,1)),                             &
        zdryrt(:,:,loop),                                                      &
        row_length,rows)
    END IF
    ! process ddep diagnostics for HONO2
  ELSE IF (speci(nldepd(loop)) == 'HONO2     ') THEN
    ! zdryrt on gridbox for HONO2
    IF (sf(437,50)) THEN
      CALL copydiag (                                                          &
        stashwork(si(437,50,1):si_last(437,50,1)),                             &
        zdryrt(:,:,loop),                                                      &
        row_length,rows)
    END IF
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_ddcalc
END MODULE ukca_ddcalc_mod
