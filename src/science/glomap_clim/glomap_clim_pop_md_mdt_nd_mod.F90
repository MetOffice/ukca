! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Populate md, mdt, nd fields
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE glomap_clim_pop_md_mdt_nd_mod

USE um_types,                                ONLY:                             &
    real_umphys

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='GLOMAP_CLIM_POP_MD_MDT_ND_MOD'

CONTAINS

SUBROUTINE glomap_clim_pop_md_mdt_nd ( i_glomap_clim_setup_in, n_points, aird, &
                                       mmr1d, nmr1d, md, mdt, nd )

USE ereport_mod,                             ONLY:                             &
    ereport

USE errormessagelength_mod,                  ONLY:                             &
    errormessagelength

USE glomap_clim_option_mod,                  ONLY:                             &
    i_gc_sussocbc_5mode,                                                       &
    i_gc_sussocbcdu_7mode

USE parkind1,                                ONLY:                             &
    jpim,                                                                      &
    jprb

USE ukca_constants,                          ONLY:                             &
    m_air

USE ukca_mode_setup,                         ONLY:                             &
    component_list_by_cp_sussbcoc_5mode,                                       &
    component_list_by_cp_sussbcocdu_7mode,                                     &
    component_list_by_mode_sussbcoc_5mode,                                     &
    component_list_by_mode_sussbcocdu_7mode,                                   &
    mfrac_0,                                                                   &
    mlo,                                                                       &
    mm,                                                                        &
    mmid,                                                                      &
    mode_list_sussbcoc_5mode,                                                  &
    mode_list_sussbcocdu_7mode,                                                &
    ncp,                                                                       &
    ncp_list_sussbcoc_5mode,                                                   &
    ncp_list_sussbcocdu_7mode,                                                 &
    nmodes,                                                                    &
    nmodes_list_sussbcoc_5mode,                                                &
    nmodes_list_sussbcocdu_7mode,                                              &
    num_eps

USE umPrintMgr,                              ONLY:                             &
    newline

USE yomhook,                                 ONLY:                             &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments
! aird      : Dry air density
! md        : Component median aerosol mass (molecules per ptcl)
! mdt       : Total median aerosol mass (molecules per ptcl)
! nd        : Aerosol ptcl no. concentration (ptcls per cc)

INTEGER, INTENT(IN) :: i_glomap_clim_setup_in
INTEGER, INTENT(IN) :: n_points
REAL, INTENT(IN)    :: aird(n_points)
REAL, INTENT(IN)    :: mmr1d(n_points,nmodes,ncp)
REAL, INTENT(IN)    :: nmr1d(n_points,nmodes)
REAL, INTENT(OUT)   :: md(n_points,nmodes,ncp)
REAL, INTENT(OUT)   :: mdt(n_points,nmodes)
REAL, INTENT(OUT)   :: nd(n_points,nmodes)

! Local variables

LOGICAL :: mask(n_points,nmodes)

INTEGER :: imode                      ! counter for modes
INTEGER :: icp                        ! counter for components
INTEGER :: m                          ! counter
INTEGER :: i

INTEGER              :: nmodes_list_local
INTEGER              :: ncp_list_local
INTEGER, ALLOCATABLE :: mode_list_local(:)
INTEGER, ALLOCATABLE :: component_list_by_mode_local(:)
INTEGER, ALLOCATABLE :: component_list_by_cp_local(:)

INTEGER                           :: ierrcode
CHARACTER(LEN=errormessagelength) :: cmessage

REAL(KIND=real_umphys) :: mdtmin(nmodes)      ! Minimum value for mdt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName='GLOMAP_CLIM_POP_MD_MDT_ND'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( ANY ( nmr1d(:,:) < 0.0 ) ) THEN
  ierrcode = 1
  WRITE(cmessage,'(A)') 'nmr1d contains negative values.'                      &
        //newline// 'Try setting l_ignore_ancil_grid_check=.false. in suite.'  &
        //newline// 'Check if NetCDF input contains negative values.'
  CALL ereport(Modulename//':'//RoutineName,ierrcode,cmessage)
END IF

IF ( ANY ( mmr1d(:,:,:) < 0.0 ) ) THEN
  ierrcode = 1
  WRITE(cmessage,'(A,I0,A,I0)') 'mmr1d contains negative values.'              &
        //newline// 'Try setting l_ignore_ancil_grid_check=.false. in suite.'  &
        //newline// 'Check if NetCDF input contains negative values.'
  CALL ereport(Modulename//':'//RoutineName,ierrcode,cmessage)
END IF

nd(:,:)=0.0
mdt(:,:)=0.0
md(:,:,:)=0.0
mask(:,:)=.FALSE.

!==============================================================================
! Define local loops depending on GLOMAP-mode setup

SELECT CASE(i_glomap_clim_setup_in)
CASE (i_gc_sussocbc_5mode)
  nmodes_list_local = nmodes_list_sussbcoc_5mode
  ncp_list_local    = ncp_list_sussbcoc_5mode
CASE (i_gc_sussocbcdu_7mode)
  nmodes_list_local = nmodes_list_sussbcocdu_7mode
  ncp_list_local    = ncp_list_sussbcocdu_7mode
CASE DEFAULT
  ierrcode = 1
  WRITE(cmessage,'(A,I0,A)') 'i_glomap_clim_setup_in = ',                      &
                              i_glomap_clim_setup_in,                          &
                              newline // 'This option not available.'
  CALL ereport(RoutineName,ierrcode,cmessage)
END SELECT

ALLOCATE( mode_list_local(nmodes_list_local) )
ALLOCATE( component_list_by_mode_local(ncp_list_local) )
ALLOCATE( component_list_by_cp_local(ncp_list_local) )

SELECT CASE(i_glomap_clim_setup_in)
CASE (i_gc_sussocbc_5mode)
  mode_list_local              = mode_list_sussbcoc_5mode
  component_list_by_mode_local = component_list_by_mode_sussbcoc_5mode
  component_list_by_cp_local   = component_list_by_cp_sussbcoc_5mode
CASE (i_gc_sussocbcdu_7mode)
  mode_list_local              = mode_list_sussbcocdu_7mode
  component_list_by_mode_local = component_list_by_mode_sussbcocdu_7mode
  component_list_by_cp_local   = component_list_by_cp_sussbcocdu_7mode
END SELECT

!==============================================================================
! Calculate nd

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(m, imode, i)          &
!$OMP SHARED(n_points, mode_list_local, nmodes_list_local, nd, nmr1d, aird,    &
!$OMP        mask, num_eps)
DO m = 1, nmodes_list_local
  imode = mode_list_local(m)

  DO i = 1, n_points
    nd(i,imode) = nmr1d(i,imode) * aird(i)
    ! Mask for ND threshold
    mask(i,imode) = ( nd(i,imode) > num_eps(imode) )
  END DO
END DO
!$OMP END PARALLEL DO

!==============================================================================
! Calculate md and mdt

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(m, imode, icp, i)     &
!$OMP SHARED(n_points, ncp_list_local, component_list_by_mode_local,           &
!$OMP        component_list_by_cp_local, mask, md, mmr1d, mm, aird, nd,        &
!$OMP        mmid, mfrac_0, mdt)
DO m = 1, ncp_list_local
  imode = component_list_by_mode_local(m)
  icp = component_list_by_cp_local(m)

  DO i = 1, n_points
    IF (mask(i,imode)) THEN
      md(i,imode,icp) = mmr1d(i,imode,icp) * ( m_air / mm(icp) ) * aird(i) /   &
                        nd(i,imode)
    ELSE
      md(i,imode,icp) = mmid(imode) * mfrac_0(imode,icp)
    END IF

    ! Set total mass array MDT from SUM over individual component MDs
    mdt(i,imode) = mdt(i,imode) + md(i,imode,icp)
  END DO
END DO
!$OMP END PARALLEL DO

!==============================================================================
! Force minimum values of mdt and nd

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(m, imode, i)          &
!$OMP SHARED(n_points, mode_list_local, nmodes_list_local, mlo, mdtmin, mask,  &
!$OMP        mdt, nd, mmid)
DO m = 1, nmodes_list_local
  imode = mode_list_local(m)

  ! set equiv. to DPLIM0*0.1
  mdtmin(imode) = mlo(imode) * 0.001

  DO i = 1, n_points
    ! Set ND -> 0 where MDT too low and set MDT -> MMID
    mask(i,imode) = ( mdt(i,imode) < mdtmin(imode) )

    IF ( mask(i,imode) ) THEN
      nd(i,imode)  = 0.0
      mdt(i,imode) = mmid(imode)
    END IF
  END DO
END DO
!$OMP END PARALLEL DO

!==============================================================================
! Deallocate local arrays

DEALLOCATE( component_list_by_cp_local )
DEALLOCATE( component_list_by_mode_local )
DEALLOCATE( mode_list_local )

!==============================================================================


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE glomap_clim_pop_md_mdt_nd

END MODULE glomap_clim_pop_md_mdt_nd_mod
