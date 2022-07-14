! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Module to contain photolysis rate arrays
!   Contains subroutines: strat_photol_init,
!   strat_photol_dealloc
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
MODULE ukca_um_dissoc_mod


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

REAL, ALLOCATABLE, SAVE :: ajhno3 (:)
REAL, ALLOCATABLE, SAVE :: ajpna  (:)
REAL, ALLOCATABLE, SAVE :: ajh2o2 (:)
REAL, ALLOCATABLE, SAVE :: aj2a   (:)
REAL, ALLOCATABLE, SAVE :: aj2b   (:)
REAL, ALLOCATABLE, SAVE :: aj3    (:)
REAL, ALLOCATABLE, SAVE :: aj3a   (:)
REAL, ALLOCATABLE, SAVE :: ajcnita(:)
REAL, ALLOCATABLE, SAVE :: ajcnitb(:)
REAL, ALLOCATABLE, SAVE :: ajbrno3(:)
REAL, ALLOCATABLE, SAVE :: ajbrcl (:)
REAL, ALLOCATABLE, SAVE :: ajoclo (:)
REAL, ALLOCATABLE, SAVE :: ajcl2o2(:)
REAL, ALLOCATABLE, SAVE :: ajhocl (:)
REAL, ALLOCATABLE, SAVE :: ajno   (:)
REAL, ALLOCATABLE, SAVE :: ajno2  (:)
REAL, ALLOCATABLE, SAVE :: ajn2o5 (:)
REAL, ALLOCATABLE, SAVE :: ajno31 (:)
REAL, ALLOCATABLE, SAVE :: ajno32 (:)
REAL, ALLOCATABLE, SAVE :: ajbro  (:)
REAL, ALLOCATABLE, SAVE :: ajhcl  (:)
REAL, ALLOCATABLE, SAVE :: ajn2o  (:)
REAL, ALLOCATABLE, SAVE :: ajhobr (:)
REAL, ALLOCATABLE, SAVE :: ajf11  (:)
REAL, ALLOCATABLE, SAVE :: ajf12  (:)
REAL, ALLOCATABLE, SAVE :: ajh2o  (:)
REAL, ALLOCATABLE, SAVE :: ajccl4 (:)
REAL, ALLOCATABLE, SAVE :: ajf113 (:)
REAL, ALLOCATABLE, SAVE :: ajf22  (:)
REAL, ALLOCATABLE, SAVE :: ajch3cl(:)
REAL, ALLOCATABLE, SAVE :: ajc2oa (:)
REAL, ALLOCATABLE, SAVE :: ajc2ob (:)
REAL, ALLOCATABLE, SAVE :: ajmhp  (:)
REAL, ALLOCATABLE, SAVE :: ajch3br(:)
REAL, ALLOCATABLE, SAVE :: ajmcfm (:)
REAL, ALLOCATABLE, SAVE :: ajch4  (:)
REAL, ALLOCATABLE, SAVE :: ajf12b1(:)
REAL, ALLOCATABLE, SAVE :: ajf13b1(:)
REAL, ALLOCATABLE, SAVE :: ajcof2 (:)
REAL, ALLOCATABLE, SAVE :: ajcofcl(:)
REAL, ALLOCATABLE, SAVE :: ajco2  (:)
REAL, ALLOCATABLE, SAVE :: ajcos  (:)
REAL, ALLOCATABLE, SAVE :: ajhono (:)
REAL, ALLOCATABLE, SAVE :: ajmena (:)
REAL, ALLOCATABLE, SAVE :: ajchbr3(:)
REAL, ALLOCATABLE, SAVE :: ajdbrm (:)
REAL, ALLOCATABLE, SAVE :: ajcs2  (:)
REAL, ALLOCATABLE, SAVE :: ajh2so4(:)
REAL, ALLOCATABLE, SAVE :: ajso3  (:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_UM_DISSOC_MOD'

CONTAINS

!
! Subroutine Interface:
!
!---------------------------------------------------------------------------
! Subroutine STRAT_PHOTOL_INIT
!------------------------------------------------------------------------
!
! This routine computes stratospheric photolysis rates and merges the
! rates, where necessary, with the tropospheric rates. This is done for
! one level at a time. The stratospheric photolysis routines are taken
! from SLIMCAT.

SUBROUTINE strat_photol_init(theta_field_size)

IMPLICIT NONE

INTEGER, INTENT(IN) :: theta_field_size

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STRAT_PHOTOL_INIT'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(ajhno3 )) ALLOCATE(ajhno3 (theta_field_size))
IF (.NOT. ALLOCATED(ajpna  )) ALLOCATE(ajpna  (theta_field_size))
IF (.NOT. ALLOCATED(ajh2o2 )) ALLOCATE(ajh2o2 (theta_field_size))
IF (.NOT. ALLOCATED(aj2a   )) ALLOCATE(aj2a   (theta_field_size))
IF (.NOT. ALLOCATED(aj2b   )) ALLOCATE(aj2b   (theta_field_size))
IF (.NOT. ALLOCATED(aj3    )) ALLOCATE(aj3    (theta_field_size))
IF (.NOT. ALLOCATED(aj3a   )) ALLOCATE(aj3a   (theta_field_size))
IF (.NOT. ALLOCATED(ajcnita)) ALLOCATE(ajcnita(theta_field_size))
IF (.NOT. ALLOCATED(ajcnitb)) ALLOCATE(ajcnitb(theta_field_size))
IF (.NOT. ALLOCATED(ajbrno3)) ALLOCATE(ajbrno3(theta_field_size))
IF (.NOT. ALLOCATED(ajbrcl )) ALLOCATE(ajbrcl (theta_field_size))
IF (.NOT. ALLOCATED(ajoclo )) ALLOCATE(ajoclo (theta_field_size))
IF (.NOT. ALLOCATED(ajcl2o2)) ALLOCATE(ajcl2o2(theta_field_size))
IF (.NOT. ALLOCATED(ajhocl )) ALLOCATE(ajhocl (theta_field_size))
IF (.NOT. ALLOCATED(ajno   )) ALLOCATE(ajno   (theta_field_size))
IF (.NOT. ALLOCATED(ajno2  )) ALLOCATE(ajno2  (theta_field_size))
IF (.NOT. ALLOCATED(ajn2o5 )) ALLOCATE(ajn2o5 (theta_field_size))
IF (.NOT. ALLOCATED(ajno31 )) ALLOCATE(ajno31 (theta_field_size))
IF (.NOT. ALLOCATED(ajno32 )) ALLOCATE(ajno32 (theta_field_size))
IF (.NOT. ALLOCATED(ajbro  )) ALLOCATE(ajbro  (theta_field_size))
IF (.NOT. ALLOCATED(ajhcl  )) ALLOCATE(ajhcl  (theta_field_size))
IF (.NOT. ALLOCATED(ajn2o  )) ALLOCATE(ajn2o  (theta_field_size))
IF (.NOT. ALLOCATED(ajhobr )) ALLOCATE(ajhobr (theta_field_size))
IF (.NOT. ALLOCATED(ajf11  )) ALLOCATE(ajf11  (theta_field_size))
IF (.NOT. ALLOCATED(ajf12  )) ALLOCATE(ajf12  (theta_field_size))
IF (.NOT. ALLOCATED(ajh2o  )) ALLOCATE(ajh2o  (theta_field_size))
IF (.NOT. ALLOCATED(ajccl4 )) ALLOCATE(ajccl4 (theta_field_size))
IF (.NOT. ALLOCATED(ajf113 )) ALLOCATE(ajf113 (theta_field_size))
IF (.NOT. ALLOCATED(ajf22  )) ALLOCATE(ajf22  (theta_field_size))
IF (.NOT. ALLOCATED(ajch3cl)) ALLOCATE(ajch3cl(theta_field_size))
IF (.NOT. ALLOCATED(ajc2oa )) ALLOCATE(ajc2oa (theta_field_size))
IF (.NOT. ALLOCATED(ajc2ob )) ALLOCATE(ajc2ob (theta_field_size))
IF (.NOT. ALLOCATED(ajmhp  )) ALLOCATE(ajmhp  (theta_field_size))
IF (.NOT. ALLOCATED(ajch3br)) ALLOCATE(ajch3br(theta_field_size))
IF (.NOT. ALLOCATED(ajmcfm )) ALLOCATE(ajmcfm (theta_field_size))
IF (.NOT. ALLOCATED(ajch4  )) ALLOCATE(ajch4  (theta_field_size))
IF (.NOT. ALLOCATED(ajf12b1)) ALLOCATE(ajf12b1(theta_field_size))
IF (.NOT. ALLOCATED(ajf13b1)) ALLOCATE(ajf13b1(theta_field_size))
IF (.NOT. ALLOCATED(ajcof2 )) ALLOCATE(ajcof2 (theta_field_size))
IF (.NOT. ALLOCATED(ajcofcl)) ALLOCATE(ajcofcl(theta_field_size))
IF (.NOT. ALLOCATED(ajco2  )) ALLOCATE(ajco2  (theta_field_size))
IF (.NOT. ALLOCATED(ajcos  )) ALLOCATE(ajcos  (theta_field_size))
IF (.NOT. ALLOCATED(ajhono )) ALLOCATE(ajhono (theta_field_size))
IF (.NOT. ALLOCATED(ajmena )) ALLOCATE(ajmena (theta_field_size))
IF (.NOT. ALLOCATED(ajchbr3)) ALLOCATE(ajchbr3(theta_field_size))
IF (.NOT. ALLOCATED(ajdbrm )) ALLOCATE(ajdbrm (theta_field_size))
IF (.NOT. ALLOCATED(ajcs2  )) ALLOCATE(ajcs2  (theta_field_size))
IF (.NOT. ALLOCATED(ajh2so4)) ALLOCATE(ajh2so4(theta_field_size))
IF (.NOT. ALLOCATED(ajso3  )) ALLOCATE(ajso3  (theta_field_size))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE strat_photol_init

!---------------------------------------------------------------------------
! Subroutine STRAT_PHOTOL_DEALLOC
!------------------------------------------------------------------------
!
! This subroutine deallocates the stratospheric photolysis rate arrays

SUBROUTINE strat_photol_dealloc()

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STRAT_PHOTOL_DEALLOC'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Deallocating in reverse order to which the arrays were allocated
IF (ALLOCATED(ajso3  )) DEALLOCATE(ajso3  )
IF (ALLOCATED(ajh2so4)) DEALLOCATE(ajh2so4)
IF (ALLOCATED(ajcs2  )) DEALLOCATE(ajcs2  )
IF (ALLOCATED(ajdbrm )) DEALLOCATE(ajdbrm )
IF (ALLOCATED(ajchbr3)) DEALLOCATE(ajchbr3)
IF (ALLOCATED(ajmena )) DEALLOCATE(ajmena )
IF (ALLOCATED(ajhono )) DEALLOCATE(ajhono )
IF (ALLOCATED(ajcos  )) DEALLOCATE(ajcos  )
IF (ALLOCATED(ajco2  )) DEALLOCATE(ajco2  )
IF (ALLOCATED(ajcofcl)) DEALLOCATE(ajcofcl)
IF (ALLOCATED(ajcof2 )) DEALLOCATE(ajcof2 )
IF (ALLOCATED(ajf13b1)) DEALLOCATE(ajf13b1)
IF (ALLOCATED(ajf12b1)) DEALLOCATE(ajf12b1)
IF (ALLOCATED(ajch4  )) DEALLOCATE(ajch4  )
IF (ALLOCATED(ajmcfm )) DEALLOCATE(ajmcfm )
IF (ALLOCATED(ajch3br)) DEALLOCATE(ajch3br)
IF (ALLOCATED(ajmhp  )) DEALLOCATE(ajmhp  )
IF (ALLOCATED(ajc2ob )) DEALLOCATE(ajc2ob )
IF (ALLOCATED(ajc2oa )) DEALLOCATE(ajc2oa )
IF (ALLOCATED(ajch3cl)) DEALLOCATE(ajch3cl)
IF (ALLOCATED(ajf22  )) DEALLOCATE(ajf22  )
IF (ALLOCATED(ajf113 )) DEALLOCATE(ajf113 )
IF (ALLOCATED(ajccl4 )) DEALLOCATE(ajccl4 )
IF (ALLOCATED(ajh2o  )) DEALLOCATE(ajh2o  )
IF (ALLOCATED(ajf12  )) DEALLOCATE(ajf12  )
IF (ALLOCATED(ajf11  )) DEALLOCATE(ajf11  )
IF (ALLOCATED(ajhobr )) DEALLOCATE(ajhobr )
IF (ALLOCATED(ajn2o  )) DEALLOCATE(ajn2o  )
IF (ALLOCATED(ajhcl  )) DEALLOCATE(ajhcl  )
IF (ALLOCATED(ajbro  )) DEALLOCATE(ajbro  )
IF (ALLOCATED(ajno32 )) DEALLOCATE(ajno32 )
IF (ALLOCATED(ajno31 )) DEALLOCATE(ajno31 )
IF (ALLOCATED(ajn2o5 )) DEALLOCATE(ajn2o5 )
IF (ALLOCATED(ajno2  )) DEALLOCATE(ajno2  )
IF (ALLOCATED(ajno   )) DEALLOCATE(ajno   )
IF (ALLOCATED(ajhocl )) DEALLOCATE(ajhocl )
IF (ALLOCATED(ajcl2o2)) DEALLOCATE(ajcl2o2)
IF (ALLOCATED(ajoclo )) DEALLOCATE(ajoclo )
IF (ALLOCATED(ajbrcl )) DEALLOCATE(ajbrcl )
IF (ALLOCATED(ajbrno3)) DEALLOCATE(ajbrno3)
IF (ALLOCATED(ajcnitb)) DEALLOCATE(ajcnitb)
IF (ALLOCATED(ajcnita)) DEALLOCATE(ajcnita)
IF (ALLOCATED(aj3a   )) DEALLOCATE(aj3a   )
IF (ALLOCATED(aj3    )) DEALLOCATE(aj3    )
IF (ALLOCATED(aj2b   )) DEALLOCATE(aj2b   )
IF (ALLOCATED(aj2a   )) DEALLOCATE(aj2a   )
IF (ALLOCATED(ajh2o2 )) DEALLOCATE(ajh2o2 )
IF (ALLOCATED(ajpna  )) DEALLOCATE(ajpna  )
IF (ALLOCATED(ajhno3 )) DEALLOCATE(ajhno3 )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE strat_photol_dealloc

END MODULE ukca_um_dissoc_mod
