! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!  Calculate and return md and mdt and nd
!
! Code Owner:
!  Please refer to the UM file CodeOwners.txt
!
! This file belongs in section:
!  UKCA
!
! Arguments:
!   aird   : Dry air density
!   mmr1d  : Avg cpt mass mixing ratio of aerosol particle in mode (particle^-1)
!   nmr1d  : Aerosol ptcl (number density/ air density) for mode (cm^-3)
!   md     : Component median aerosol mass (molecules per ptcl)
!   mdt    : Total median aerosol mass (molecules per ptcl)
!   nd     : Aerosol ptcl no. concentration (ptcls per cc)
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE ukca_calc_md_mdt_nd_mod

USE um_types,                                ONLY:                             &
    real_umphys

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='UKCA_CALC_MD_MDT_ND_MOD'

CONTAINS

SUBROUTINE ukca_calc_md_mdt_nd ( i_mode_setup_local, n_points,                 &
                                 glomap_variables_local,                       &
                                 aird, mmr1d, nmr1d,                           &
                                 md, mdt, nd )

USE ereport_mod,                             ONLY:                             &
    ereport

USE errormessagelength_mod,                  ONLY:                             &
    errormessagelength

USE parkind1,                                ONLY:                             &
    jpim,                                                                      &
    jprb

USE ukca_constants,                          ONLY:                             &
    m_air

USE ukca_config_specification_mod,           ONLY:                             &
    i_sussbcoc_5mode,                                                          &
    i_sussbcocdu_7mode

USE ukca_mode_setup,                         ONLY:                             &
    component_list_by_cp_sussbcoc_5mode,                                       &
    component_list_by_cp_sussbcocdu_7mode,                                     &
    component_list_by_mode_sussbcoc_5mode,                                     &
    component_list_by_mode_sussbcocdu_7mode,                                   &
    glomap_variables_type,                                                     &
    mode_list_sussbcoc_5mode,                                                  &
    mode_list_sussbcocdu_7mode,                                                &
    ncp_list_sussbcoc_5mode,                                                   &
    ncp_list_sussbcocdu_7mode,                                                 &
    nmodes,                                                                    &
    nmodes_list_sussbcoc_5mode,                                                &
    nmodes_list_sussbcocdu_7mode

USE umPrintMgr,                              ONLY:                             &
    newline

USE yomhook,                                 ONLY:                             &
    lhook,                                                                     &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: i_mode_setup_local
INTEGER, INTENT(IN) :: n_points
TYPE(glomap_variables_type), TARGET, INTENT(IN) :: glomap_variables_local

!   aird   : Dry air density
REAL, INTENT(IN)    :: aird(  n_points )

!   mmr1d  : Avg cpt mass mixing ratio of aerosol particle in mode (particle^-1)
REAL, INTENT(IN)    :: mmr1d( n_points, nmodes, glomap_variables_local%ncp )

!   nmr1d  : Aerosol ptcl (number density/ air density) for mode (cm^-3)
REAL, INTENT(IN)    :: nmr1d( n_points, nmodes )

!   md     : Component median aerosol mass (molecules per ptcl)
REAL, INTENT(OUT)   :: md(    n_points, nmodes, glomap_variables_local%ncp )

!   mdt    : Total median aerosol mass (molecules per ptcl)
REAL, INTENT(OUT)   :: mdt(   n_points, nmodes )

!   nd     : Aerosol ptcl no. concentration (ptcls per cc)
REAL, INTENT(OUT)   :: nd(    n_points, nmodes )

! Local variables

LOGICAL :: mask(n_points,nmodes)

INTEGER :: imode                      ! counter for modes
INTEGER :: icp                        ! counter for components
INTEGER :: loop                       ! counter for n_points
INTEGER :: m                          ! counter for local list

INTEGER              :: nmodes_list_local
INTEGER              :: ncp_list_local
INTEGER, ALLOCATABLE :: mode_list_local(:)
INTEGER, ALLOCATABLE :: component_list_by_mode_local(:)
INTEGER, ALLOCATABLE :: component_list_by_cp_local(:)

REAL :: m_air_div_mm( glomap_variables_local%ncp )
REAL :: mmid_times_mfrac_0( nmodes , glomap_variables_local%ncp )

INTEGER                           :: ierrcode
CHARACTER(LEN=errormessagelength) :: cmessage

REAL(KIND=real_umphys) :: mdtmin(nmodes)      ! Minimum value for mdt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName='UKCA_CALC_MD_MDT_ND'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!==============================================================================
! Error traps for negative concentrations of mmr1d and nmr1d

IF ( ANY ( nmr1d(:,:) < 0.0 ) ) THEN
  ierrcode = 1
  WRITE(cmessage,'(A)') 'nmr1d contains negative values.'                      &
        //newline// 'Try setting l_ignore_ancil_grid_check=.false. in suite.'  &
        //newline// 'Check if NetCDF input contains negative values.'
  CALL ereport(Modulename//':'//RoutineName,ierrcode,cmessage)
END IF

IF ( ANY ( mmr1d(:,:,:) < 0.0 ) ) THEN
  ierrcode = 1
  WRITE(cmessage,'(A)') 'mmr1d contains negative values.'                      &
        //newline// 'Try setting l_ignore_ancil_grid_check=.false. in suite.'  &
        //newline// 'Check if NetCDF input contains negative values.'
  CALL ereport(Modulename//':'//RoutineName,ierrcode,cmessage)
END IF

!==============================================================================
! Define local loops depending on GLOMAP-mode setup

SELECT CASE(i_mode_setup_local)
CASE (i_sussbcoc_5mode)
  nmodes_list_local = nmodes_list_sussbcoc_5mode
  ncp_list_local    = ncp_list_sussbcoc_5mode

CASE (i_sussbcocdu_7mode)
  nmodes_list_local = nmodes_list_sussbcocdu_7mode
  ncp_list_local    = ncp_list_sussbcocdu_7mode

CASE DEFAULT
  ierrcode = 1
  WRITE(cmessage,'(A,I0,A)') 'i_mode_setup_local = ',                          &
                              i_mode_setup_local,                              &
                              newline // 'This option not available.'
  CALL ereport(RoutineName,ierrcode,cmessage)

END SELECT

ALLOCATE( mode_list_local(nmodes_list_local) )
ALLOCATE( component_list_by_mode_local(ncp_list_local) )
ALLOCATE( component_list_by_cp_local(ncp_list_local) )

SELECT CASE(i_mode_setup_local)
CASE (i_sussbcoc_5mode)
  mode_list_local              = mode_list_sussbcoc_5mode
  component_list_by_mode_local = component_list_by_mode_sussbcoc_5mode
  component_list_by_cp_local   = component_list_by_cp_sussbcoc_5mode

CASE (i_sussbcocdu_7mode)
  mode_list_local              = mode_list_sussbcocdu_7mode
  component_list_by_mode_local = component_list_by_mode_sussbcocdu_7mode
  component_list_by_cp_local   = component_list_by_cp_sussbcocdu_7mode

END SELECT

!==============================================================================
! Calculations not made over n_points

DO m = 1, ncp_list_local
  imode = component_list_by_mode_local(m)
  icp = component_list_by_cp_local(m)

  m_air_div_mm(icp) = m_air / glomap_variables_local%mm(icp)

  mmid_times_mfrac_0(imode,icp) = glomap_variables_local%mmid(imode) *         &
                                  glomap_variables_local%mfrac_0(imode,icp)

END DO

DO m = 1, nmodes_list_local
  imode = mode_list_local(m)

  ! set equiv. to DPLIM0*0.1
  mdtmin(imode) = glomap_variables_local%mlo(imode) * 0.001
END DO

!==============================================================================
! Initialise arrays

! Initialise mdt(:,:) and nd(:,:) to zero everywhere
DO imode = 1, nmodes
  DO loop = 1, n_points
    mdt(loop,imode)  = 0.0
    nd(loop,imode)   = 0.0
    mask(loop,imode) = .FALSE.
  END DO
END DO

! Initialise md(:,:,:) to zero everywhere
DO icp = 1, glomap_variables_local%ncp
  DO imode = 1, nmodes
    DO loop = 1, n_points
      md(loop,imode,icp) = 0.0
    END DO
  END DO
END DO

!==============================================================================
! Calculate nd

DO m = 1, nmodes_list_local
  imode = mode_list_local(m)

  DO loop = 1, n_points
    nd(loop,imode) = nmr1d(loop,imode) * aird(loop)

    ! Mask for ND threshold
    mask(loop,imode) = (nd(loop,imode) > glomap_variables_local%num_eps(imode))
  END DO

END DO

!==============================================================================
! Calculate md and mdt

DO m = 1, ncp_list_local
  icp = component_list_by_cp_local(m)
  imode = component_list_by_mode_local(m)

  DO loop = 1, n_points
    IF (mask(loop,imode)) THEN
      md(loop,imode,icp) = mmr1d(loop,imode,icp) * m_air_div_mm(icp) *         &
                           aird(loop) / nd(loop,imode)
    ELSE
      md(loop,imode,icp) = mmid_times_mfrac_0(imode,icp)
    END IF

    ! Set total mass array MDT from SUM over individual component MDs
    mdt(loop,imode) = mdt(loop,imode) + md(loop,imode,icp)
  END DO

END DO

!==============================================================================
! Force minimum values of mdt and nd

DO m = 1, nmodes_list_local
  imode = mode_list_local(m)

  DO loop = 1, n_points
    ! Set ND -> 0 where MDT too low and set MDT -> MMID
    mask(loop,imode) = ( mdt(loop,imode) < mdtmin(imode) )

    IF ( mask(loop,imode) ) THEN
      nd(loop,imode)  = 0.0
      mdt(loop,imode) = glomap_variables_local%mmid(imode)
    END IF
  END DO
END DO

!==============================================================================
! Deallocate local arrays

DEALLOCATE( component_list_by_cp_local )
DEALLOCATE( component_list_by_mode_local )
DEALLOCATE( mode_list_local )

!==============================================================================


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calc_md_mdt_nd

END MODULE ukca_calc_md_mdt_nd_mod
