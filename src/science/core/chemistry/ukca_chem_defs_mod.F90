! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module providing UKCA with type definitions describing chemistry
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
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_chem_defs_mod

IMPLICIT NONE

TYPE :: chch_t
  INTEGER           :: item_no   ! dummy integer
  CHARACTER(LEN=10) :: speci     ! species name
  INTEGER           :: nodd      ! no of odd atoms
  CHARACTER(LEN=10) :: ctype     ! Species type
  CHARACTER(LEN=10) :: family    ! Family
  INTEGER           :: switch1   ! 1 if dry deposits
  INTEGER           :: switch2   ! 1 if wet deposits
END TYPE chch_t
PUBLIC chch_t

! Describes the bimolecular reaction rates
TYPE :: ratb_t
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL              :: k0        ! rate coeff param K0
  REAL              :: alpha     ! rate coeff param alpha
  REAL              :: beta      ! rate coeff param beta
  REAL              :: pyield1   ! product1 fractional yield
  REAL              :: pyield2   ! product2 fractional yield
  REAL              :: pyield3   ! product3 fractional yield
  REAL              :: pyield4   ! product4 fractional yield
END TYPE ratb_t
PUBLIC ratb_t

! Describes heterogeneous reactions
TYPE :: rath_t
  CHARACTER(LEN=10) :: react1    ! reactant1 name
  CHARACTER(LEN=10) :: react2    ! reactant2 name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL              :: pyield1   ! product yield
  REAL              :: pyield2   ! product yield
  REAL              :: pyield3   ! product yield
  REAL              :: pyield4   ! product yield
END TYPE rath_t
PUBLIC rath_t

! describes photolytic reactions
INTEGER, PARAMETER, PUBLIC :: photol_varname_len = 10
TYPE :: ratj_t
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL              :: pyield1   ! product yield
  REAL              :: pyield2   ! product yield
  REAL              :: pyield3   ! product yield
  REAL              :: pyield4   ! product yield
  REAL              :: jfacta    ! quantum yield
  ! file name/label
  CHARACTER(LEN=photol_varname_len) :: fname
END TYPE ratj_t
PUBLIC ratj_t

! Describes termolecular reactions
TYPE :: ratt_t
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  REAL              :: f         ! rate coeff param F
  REAL              :: k1        ! rate coeff param K1
  REAL              :: alpha1    ! rate coeff param alpha1
  REAL              :: beta1     ! rate coeff param beta1
  REAL              :: k2        ! rate coeff param K2
  REAL              :: alpha2    ! rate coeff param alpha2
  REAL              :: beta2     ! rate coeff param beta2
  REAL              :: pyield1   ! product1 fractional yield
  REAL              :: pyield2   ! product2 fractional yield
END TYPE ratt_t
PUBLIC ratt_t

! Describes the chemical species used in the model
TYPE :: chch_t1
  INTEGER           :: item_no   ! dummy integer
  CHARACTER(LEN=10) :: speci     ! species name
  INTEGER           :: nodd      ! no of odd atoms
  CHARACTER(LEN=10) :: ctype     ! Species type
  CHARACTER(LEN=10) :: family    ! Family
  INTEGER           :: switch1   ! 1 if dry deposits
  INTEGER           :: switch2   ! 1 if wet deposits
  INTEGER           :: scheme    ! defines the chemistry schemes the species
                                 ! is part of
  INTEGER           :: qual      ! qualifier e.g. to identify aerosol reactions
  INTEGER           :: disqual   ! disqualifier (used to de-select entries)
  INTEGER           :: version
END TYPE chch_t1
PUBLIC chch_t1

! Describes the bimolecular reaction rates
TYPE :: ratb_t1
  INTEGER           :: item_no
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL              :: k0        ! rate coeff param K0
  REAL              :: alpha     ! rate coeff param alpha
  REAL              :: beta      ! rate coeff param beta
  REAL              :: pyield1   ! product1 fractional yield
  REAL              :: pyield2   ! product2 fractional yield
  REAL              :: pyield3   ! product3 fractional yield
  REAL              :: pyield4   ! product4 fractional yield
  INTEGER           :: scheme    ! defines the chemistry schemes
  INTEGER           :: qual
  INTEGER           :: disqual
  INTEGER           :: version
END TYPE ratb_t1
PUBLIC ratb_t1

! Describes heterogeneous reactions
TYPE :: rath_t1
  INTEGER           :: item_no
  CHARACTER(LEN=10) :: react1    ! reactant1 name
  CHARACTER(LEN=10) :: react2    ! reactant2 name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL              :: pyield1   ! product yield
  REAL              :: pyield2   ! product yield
  REAL              :: pyield3   ! product yield
  REAL              :: pyield4   ! product yield
  INTEGER           :: scheme    ! chemistry scheme
  INTEGER           :: qual
  INTEGER           :: disqual
  INTEGER           :: version
END TYPE rath_t1
PUBLIC rath_t1

! describes photolytic reactions
TYPE :: ratj_t1
  INTEGER           :: item_no
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  CHARACTER(LEN=10) :: prod3     ! product3 name
  CHARACTER(LEN=10) :: prod4     ! product4 name
  REAL              :: pyield1   ! product yield
  REAL              :: pyield2   ! product yield
  REAL              :: pyield3   ! product yield
  REAL              :: pyield4   ! product yield
  REAL              :: jfacta    ! quantum yield
  CHARACTER(LEN=10) :: fname     ! file name/label
  INTEGER           :: scheme    ! chemistry scheme
  INTEGER           :: qual
  INTEGER           :: disqual
  INTEGER           :: version
END TYPE ratj_t1
PUBLIC ratj_t1

! Describes termolecular reactions
TYPE :: ratt_t1
  INTEGER           :: item_no
  CHARACTER(LEN=10) :: react1    ! reactant name
  CHARACTER(LEN=10) :: react2    ! reactant name
  CHARACTER(LEN=10) :: prod1     ! product1 name
  CHARACTER(LEN=10) :: prod2     ! product2 name
  REAL              :: f         ! rate coeff param F
  REAL              :: k1        ! rate coeff param K1
  REAL              :: alpha1    ! rate coeff param alpha1
  REAL              :: beta1     ! rate coeff param beta1
  REAL              :: k2        ! rate coeff param K2
  REAL              :: alpha2    ! rate coeff param alpha2
  REAL              :: beta2     ! rate coeff param beta2
  REAL              :: pyield1   ! product1 fractional yield
  REAL              :: pyield2   ! product2 fractional yield
  INTEGER           :: scheme    ! chemistry scheme
  INTEGER           :: qual
  INTEGER           :: disqual
  INTEGER           :: version
END TYPE ratt_t1
PUBLIC ratt_t1

! Describes dry deposition reactions
TYPE :: depvel_t
  INTEGER             :: item_no
  CHARACTER(LEN=10)   :: speci        ! species name to be dry deposited
  REAL                :: velocity(30) ! dry deposition velocity
  INTEGER             :: scheme       ! chemistry scheme
  INTEGER             :: qual
  INTEGER             :: disqual
  INTEGER             :: version
END TYPE depvel_t
PUBLIC depvel_t

! Describes wet deposition parameters.
TYPE :: wetdep
  INTEGER           :: item_no
  CHARACTER(LEN=10) :: speci     ! species name to be wet deposited
  REAL              :: params(6) ! wet deposition parameters
  INTEGER           :: scheme    ! chemistry scheme
  INTEGER           :: qual
  INTEGER           :: disqual
  INTEGER           :: version
END TYPE wetdep
PUBLIC wetdep

! These are set, depending on chemistry selection, in routine
! chem1_init (contained in ukca_chem1_dat)
TYPE(chch_t),   ALLOCATABLE, PUBLIC, SAVE :: chch_defs(:)
TYPE(ratb_t),   ALLOCATABLE, PUBLIC, SAVE :: ratb_defs(:)
TYPE(rath_t),   ALLOCATABLE, PUBLIC, SAVE :: rath_defs(:)
TYPE(ratj_t),   ALLOCATABLE, PUBLIC, SAVE :: ratj_defs(:)
TYPE(ratt_t),   ALLOCATABLE, PUBLIC, SAVE :: ratt_defs(:)
REAL,           ALLOCATABLE, PUBLIC, SAVE :: depvel_defs(:,:,:)
REAL,           ALLOCATABLE, PUBLIC, SAVE :: henry_defs(:,:)
TYPE(chch_t1),  ALLOCATABLE, PUBLIC, SAVE :: chch_defs_new(:)
TYPE(ratb_t1),  ALLOCATABLE, PUBLIC, SAVE :: ratb_defs_new(:)
TYPE(rath_t1),  ALLOCATABLE, PUBLIC, SAVE :: rath_defs_new(:)
TYPE(ratj_t1),  ALLOCATABLE, PUBLIC, SAVE :: ratj_defs_new(:)
TYPE(ratt_t1),  ALLOCATABLE, PUBLIC, SAVE :: ratt_defs_new(:)
REAL,           ALLOCATABLE, PUBLIC, SAVE :: depvel_defs_new(:,:,:)
REAL,           ALLOCATABLE, PUBLIC, SAVE :: henry_defs_new(:,:)

! List of photolysis file/variable names and yield required for the
! UKCA configuration
CHARACTER(LEN=10), ALLOCATABLE, TARGET, SAVE :: ratj_varnames(:)
CHARACTER(LEN=10), ALLOCATABLE, TARGET, SAVE :: ratj_data(:,:)
REAL, ALLOCATABLE, TARGET, SAVE :: ratj_jfacta(:)

CONTAINS

SUBROUTINE ukca_get_photol_reaction_data(data_ptr, varnames_ptr, jfacta_ptr)
!
! Purpose: UKCA API subroutine, used to pass photolysis chemistry reaction
! data from ratj_defs to parent model
! ------------------------------------------------------------------------------
!
USE asad_mod, ONLY: jpspj ! asad_mod PARAMETER, number of reactants/products

IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=10), POINTER, INTENT(OUT) :: data_ptr(:,:)
CHARACTER(LEN=photol_varname_len), POINTER, INTENT(OUT) :: varnames_ptr(:)
! Yield rates- optional for backward compatibility for parent models calling
! this routine.
REAL, POINTER, OPTIONAL, INTENT(OUT) :: jfacta_ptr(:)

! Local variables
INTEGER :: i
INTEGER :: jppj ! number of photolysis reactions

jppj = SIZE(ratj_defs)

IF (.NOT. ALLOCATED(ratj_data)) THEN
  ALLOCATE(ratj_data(jppj, jpspj))
  DO i = 1, jppj
    ratj_data(i, 1) = ratj_defs(i)%react1
    ratj_data(i, 2) = ratj_defs(i)%react2
    ratj_data(i, 3) = ratj_defs(i)%prod1
    ratj_data(i, 4) = ratj_defs(i)%prod2
    ratj_data(i, 5) = ratj_defs(i)%prod3
    ratj_data(i, 6) = ratj_defs(i)%prod4
  END DO
END IF

IF (.NOT. ALLOCATED(ratj_varnames)) THEN
  ALLOCATE(ratj_varnames(jppj))
  DO i = 1, jppj
    ratj_varnames(i) = ratj_defs(i)%fname
  END DO
END IF

! Assign pointers to the data list, variable list and quantum yields (optional)
data_ptr => ratj_data(:,:)
varnames_ptr => ratj_varnames(:)

IF (PRESENT(jfacta_ptr)) THEN
  IF (.NOT. ALLOCATED(ratj_jfacta)) THEN
    ALLOCATE(ratj_jfacta(jppj))
    DO i = 1, jppj
      ratj_jfacta(i) = ratj_defs(i)%jfacta
    END DO
  END IF
  jfacta_ptr => ratj_jfacta(:)
END IF

RETURN
END SUBROUTINE ukca_get_photol_reaction_data
END MODULE ukca_chem_defs_mod
