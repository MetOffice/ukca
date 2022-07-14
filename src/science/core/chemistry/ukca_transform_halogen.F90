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
!  Description:
!    To transform one compound each for Br, Cl and N to become
!    the tracer for the total elemental abundance, and to reverse that
!    transformation.
!    Hydrogen lumping is not included at this version
!
!  Method:
!    If the logical flag forward is true then the code unlumps the tracers.
!    The tracers input are the mass mixing ratios of the individual
!    family members except the 'major' member of each family which is
!    a lumped tracer representing the total elemental mass of the whole
!    family, convered to mass mixing ratio as if it had the molecular
!    mass of the major family member.
!    The code then rescales the concentrations of the individual members
!    to ensure mass conservation of the family and changes the concentration
!    of the major family member to just be the concentration of that member.
!
!    If logical flag forward is false then the code takes in the
!    MMR of the individual tracers, sums them up and modifies the major
!    tracer to a lumped tracer representing the total elemental mass of the
!    whole family, converting to mass mixing ratio as if it had the molecular
!    mass of the major family member.
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
!    Language:  Fortran
!
! ######################################################################
!
MODULE ukca_transform_halogen_mod

USE ukca_constants, ONLY: c_hcl, c_cl, c_clo, c_cl2o2, c_hocl,                 &
                          c_oclo, c_clono2, c_cfcl3,                           &
                          c_cf2cl2, c_mecl, c_mecl, c_ccl4,                    &
                          c_cf2clcfcl2, c_chf2cl,  c_hbr,                      &
                          c_meccl3, c_brcl, c_cf2clbr, c_br,                   &
                          c_hobr, c_brono2, c_mebr, c_cf3br,                   &
                          c_ch2br2, c_n, c_no, c_no3, c_n2o5,                  &
                          c_ho2no2, c_hono2, c_hono, c_meono2,                 &
                          c_pan, c_ppan, c_ison, c_mpan,                       &
                          c_bro, c_no2

IMPLICIT NONE

PRIVATE

PUBLIC :: ukca_transform_halogen

! Type of compounds set in ukca_setup_lumping

! Unrestricted - no special behaviour
INTEGER, PARAMETER :: unrestricted   = 0

! Gases which have a source (emission/lower
! boundary condition)
INTEGER, PARAMETER :: source_gas     = 1

! Mixed compound - in multiple families,
! set this flag for all except one
! family so that its mass is only
! adjusted on unlumping in one family
! If in Br family, then should be unrestricted
! in that family, if in Cl and not Br, then
! Cl is unrestricted. For mixed compounds
! containing N, it should always be a mixed_compound
! in the N family
!
INTEGER, PARAMETER :: mixed_compound = 2

! Positions in tracer array of selected tracers
INTEGER, SAVE :: n_hcl =0
INTEGER, SAVE :: n_bro =0
INTEGER, SAVE :: n_no2 =0

! Number of compounds of conserved elements
INTEGER, SAVE :: n_cl_tracers=0
INTEGER, SAVE :: n_br_tracers=0
INTEGER, SAVE :: n_n_tracers =0

! the maximum number of tracers in any one
! family
INTEGER, PARAMETER :: tr_max = 30

! Position in tracer array of compounds
INTEGER, SAVE :: cl_tracers(tr_max)=0
INTEGER, SAVE :: br_tracers(tr_max)=0
INTEGER, SAVE :: n_tracers (tr_max)=0

! Status of compounds
INTEGER, SAVE :: s_cl_tracers(tr_max)=unrestricted
INTEGER, SAVE :: s_br_tracers(tr_max)=unrestricted
INTEGER, SAVE :: s_n_tracers (tr_max)=unrestricted

! Conversion factors for VMR to MMR
REAL, SAVE :: c_cl_tracers(tr_max)
REAL, SAVE :: c_br_tracers(tr_max)
REAL, SAVE :: c_n_tracers (tr_max)

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName='UKCA_TRANSFORM_HALOGEN_MOD'

CONTAINS

SUBROUTINE ukca_transform_halogen(rows, row_length, model_levels, tr_ukca,     &
                                  forward, tracers)

USE parkind1,       ONLY: jprb, jpim
USE yomhook,        ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine interface

INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: tr_ukca
LOGICAL, INTENT(IN) :: forward
REAL, INTENT(IN OUT) :: tracers(row_length, rows, model_levels, tr_ukca)

! Local variables

LOGICAL, SAVE :: first=.TRUE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_TRANSFORM_HALOGEN'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first) THEN
  CALL ukca_setup_lumping()
  first = .FALSE.
END IF !first

! Call transform for each family in turn. Note that order is important
! and different depending on whether unlumping (forward) or lumping
IF (forward) THEN

  ! ************************************************************
  ! Call Transform for nitrogen family.
  ! ************************************************************
  CALL transform(row_length, rows, model_levels, tr_ukca,                      &
                 n_n_tracers, n_tracers,                                       &
                 c_n_tracers, s_n_tracers, n_no2, c_no2,                       &
                 forward, tracers)

  ! ************************************************************
  ! Call Transform for chlorine family.
  ! ************************************************************
  CALL transform(row_length, rows, model_levels, tr_ukca,                      &
                 n_cl_tracers, cl_tracers,                                     &
                 c_cl_tracers, s_cl_tracers, n_hcl, c_hcl,                     &
                 forward, tracers)

  ! ************************************************************
  ! Call Transform for bromine family.
  ! ************************************************************
  CALL transform(row_length, rows, model_levels, tr_ukca,                      &
                 n_br_tracers, br_tracers,                                     &
                 c_br_tracers, s_br_tracers, n_bro, c_bro,                     &
                 forward, tracers)

ELSE

  ! ************************************************************
  ! Call Transform for bromine family.
  ! ************************************************************
  CALL transform(row_length, rows, model_levels, tr_ukca,                      &
                 n_br_tracers, br_tracers,                                     &
                 c_br_tracers, s_br_tracers, n_bro, c_bro,                     &
                 forward, tracers)

  ! ************************************************************
  ! Call Transform for chlorine family.
  ! ************************************************************
  CALL transform(row_length, rows, model_levels, tr_ukca,                      &
                 n_cl_tracers, cl_tracers,                                     &
                 c_cl_tracers, s_cl_tracers, n_hcl, c_hcl,                     &
                 forward, tracers)

  ! ************************************************************
  ! Call Transform for nitrogen family.
  ! ************************************************************
  CALL transform(row_length, rows, model_levels, tr_ukca,                      &
                 n_n_tracers, n_tracers,                                       &
                 c_n_tracers, s_n_tracers, n_no2, c_no2,                       &
                 forward, tracers)

END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_transform_halogen

! ************************************************************
SUBROUTINE transform(row_length, rows, model_levels, tr_ukca,                  &
                     n_tracers, iposition,                                     &
                     c_tracers, s_tracers, major, c_major,                     &
                     forward, tracers)
!
! Purpose:
! Performs the unlumping (when forward is true)
! or lumping (when forward is false) for a single family
!
! Method:
! Called for each family in turn
! N.B. order is important due to some species being in
! multiple families)
!
! Unlumping
!
! Initially we have MMR of individual tracers except
! for the major tracer which contains the total MMR of all
! family members. We need to rescale all members of
! the family so that their sum is the same as the major
! tracer and then subtract then from major tracer so
! that the major member represents the actual concentration
! of the single species instead of the sum of all of them.
! All calculations are done in VMR.
!
! 1. calculate the sum of the source_gas, mixed_compound and
!    unrestricted tracers
! 2. subtract the mass of mixed tracers from the family
!    leaving the source_gas and unrestricted tracers
! 3. calculate what fraction of the major (transported) tracer
!    the 'source_gas' tracers represent.
! 4. at those points where the source_gas tracers are
!    more than 100% of the transported tracer
! a. set the individual unrestricted tracers at that point
!    to a very small number (as all the mass is considered to be
!    in the source gases)
! b. rescale all the source_gas tracers so their sum is equal to that
!    of the transported tracer
! c. reset the total of the transported tracer to a very small
!    number
! 5. re-calculate the major tracer by subtracting the
!    sum of all the 'source_gas' tracers).
! 6. at those points where the fraction of the 'unrestricted' tracers within
!    the major tracer, (after the 'source_gas' species have been removed)
!    is greater than one rescale the all the unrestricted tracers
!    so that their mass is equal to that of the major tracer
! 7. recalculate the major (transported) tracer by removing the sum of all the
!    'unrestricted' tracers. At this point the major tracer
!    is the MMR of the single species chosen to be
!
! Lumping
!
! Initially we have the MMR of all tracers individually
! and we need to calculate the MMR of the sum of all
! family members and set the MMR of the major tracer to that
! sum
!
! Convert each species in the family from
! MMR to VMR and do a sum to give VMR of that family.
! Then convert the family tracer back to MMR for use in the
! transport scheme. If a species is in multiple families
! it is considered to contribute to both.
!
!
! ************************************************************

USE parkind1,       ONLY: jprb, jpim
USE yomhook,        ONLY: lhook, dr_hook

IMPLICIT NONE

! I/O variables
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: tr_ukca
INTEGER, INTENT(IN) :: n_tracers

INTEGER, INTENT(IN) :: major
REAL, INTENT(IN) :: c_major
LOGICAL, INTENT(IN) :: forward

INTEGER, INTENT(IN) :: iposition(n_tracers)
REAL, INTENT(IN) :: c_tracers(n_tracers)
INTEGER, INTENT(IN) :: s_tracers(n_tracers)
REAL, INTENT(IN OUT) :: tracers(row_length, rows, model_levels, tr_ukca)

! Local variables

REAL, PARAMETER :: rafeps = 1e-150

INTEGER :: i

REAL, ALLOCATABLE :: total_all(:,:,:)
REAL, ALLOCATABLE :: total_source(:,:,:)
REAL, ALLOCATABLE :: total_nochange(:,:,:)
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TRANSFORM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE(total_all(row_length,rows,model_levels))
ALLOCATE(total_source(row_length,rows,model_levels))
ALLOCATE(total_nochange(row_length,rows,model_levels))

! Sum up fractional contributions to total.
! Make sure they do not sum up to more than 1. If they do, rescale
! chlorine tracers.

total_all = 0.0
total_source = 0.0
total_nochange = 0.0


! add up tracers (in mmr (species)) converting to vmr (Br/Cl/N) when
! summing into the total_all ('unrestricted' species), total_source
! ('source_gas' species), and total_nochange ('mixed_compound' species), as
! defined by the particular element being considered (e.g. Br/Cl/N).
DO i=1,n_tracers
  SELECT CASE (s_tracers(i))
  CASE (unrestricted)
    total_all = total_all +                                                    &
      tracers(:,:,:,iposition(i)) / c_tracers(i)
  CASE (source_gas)
    total_source = total_source +                                              &
      tracers(:,:,:,iposition(i)) / c_tracers(i)
  CASE (mixed_compound)
    total_nochange = total_nochange +                                          &
      tracers(:,:,:,iposition(i)) / c_tracers(i)
  END SELECT
END DO

IF (forward) THEN
  !-----------------------------------------------------------------------------
  ! Unlumping (as called at start of UKCA_MAIN1) - subtract all tracers
  ! Subtract non-changing compounds and adjust where necessary.
  ! total_nochange contains 'mixed_compound' tracers, e.g.
  ! ClONO2, MeCl, CHF2Cl, MeCCl3, BrCl, CF2ClBr, etc. It should be noted that
  ! for Br tracers these may be set to source_gas if they have sources.
  ! By default all tracers are set to be 'unrestricted' and are only set to
  ! 'source_gas' or 'mixed_compound' when needed.
  !-----------------------------------------------------------------------------

  ! Remove mixed_compound tracers from transported tracer (i.e. total)
  tracers(:,:,:,major) = tracers(:,:,:,major) -                                &
     c_major * total_nochange

  !-----------------------------------------------------------------------------
  ! if 'major' (i.e. transported) tracer is less than rafeps (currently 1.0e-150,
  ! set at top of routine) then set either 'unrestricted' or 'source_gas'
  ! tracers which are added to major tracer equal to rafeps in the same gridcell
  !
  ! Reset the individual source_gas and unrestricted tracers to {a very
  ! small number} if the transported tracer is now <{a very small number}
  ! since all mass at that point is from mixed_compound tracers
  !-----------------------------------------------------------------------------
  DO i=1,n_tracers
    IF (s_tracers(i) < mixed_compound) THEN
      WHERE (tracers(:,:,:,major) < rafeps)                                    &
        tracers(:,:,:,iposition(i)) = rafeps
    END IF
  END DO

  !-----------------------------------------------------------------------------
  ! NOTE: the major tracer is not held within the iposition(i) array,
  !       and so it will not be re-set in the block above
  !-----------------------------------------------------------------------------
  ! where the major tracer is smaller than rafeps then reset it to rafeps.
  ! Similarly reset the totals for the 'source_gas' and 'unrestricted' tracers
  ! at the same point, since all the tracer in that gridcell is classed as
  ! coming from the 'mixed_compound' species.
  !
  ! Now reset transported tracer to {a very small number} if it is <{a very
  ! small number} after the mixed_compound tracers have been removed
  ! reset the totals for the source_gas and unrestricted tracers at the same
  ! point
  !------------------------------------------------------------------------------
  WHERE (tracers(:,:,:,major) < rafeps)
    tracers(:,:,:,major) = rafeps
    total_source = rafeps
    total_all = rafeps
  END WHERE

  ! Subtract source gas compounds and rescale where necessary
  ! NOTE: total_nochange is now the fraction of the 'source_gas' species
  !       within the major (transported) tracer (calculated from mmr, not
  !       that this makes a difference).
  total_nochange = total_source * c_major /                                    &
    MAX(tracers(:,:,:,major),rafeps)

  ! loop over all tracers to reset unrestricted and source gases
  DO i=1,n_tracers
    SELECT CASE (s_tracers(i))

      ! where the fraction of the 'source_gas' species
      ! within the major (transported) tracer (calculated above)
      ! is >1 then 'unrestricted' tracers
      ! are set to rafeps since all the mass at that point is classed as coming
      ! from the 'source_gas' tracers
      ! Reset the individual unrestricted tracers to {a very small number} since
      ! all mass at that point is from source_gas tracers
    CASE (unrestricted)
      WHERE (total_nochange > 1.0)                                             &
        tracers(:,:,:,iposition(i)) = rafeps

      ! where the fraction of the 'source_gas' species
      ! within the major (transported) tracer (calculated above)
      ! >1 then 'source_gas' tracers
      ! are rescaled by dividing by that fraction (this only affects the
      ! individual gridcells where this occurs)
      ! Rescale all source_gas tracers equally if their total is greater than
      ! that in the transported tracer after the mixed_compounds have been
      ! removed
    CASE (source_gas)
      WHERE (total_nochange > 1.0)                                             &
        tracers(:,:,:,iposition(i)) =                                          &
          tracers(:,:,:,iposition(i)) / total_nochange
    END SELECT
  END DO

  ! reset total_all to be rafeps where the fraction is >1, i.e. all points
  ! in that gridcell are 'source_gas'
  WHERE (total_nochange > 1.0) total_all = rafeps

  ! re-calculate the major tracer by subtracting the total_source (the
  ! sum of all the 'source_gas' tracers). If this value is less than rafeps
  ! then set to rafeps.
  tracers(:,:,:,major) = MAX(tracers(:,:,:,major) -                            &
    c_major * total_source,rafeps)

  ! Rescale remaining gases
  ! NOTE: total_nochange is now changing its definition again
  ! total_nochange is now the fraction of the 'unrestricted' tracers within
  ! the major tracer, after the 'source_gas' species have been removed.
  total_nochange = total_all * c_major /                                       &
     MAX(tracers(:,:,:,major),rafeps)

  DO i=1,n_tracers

    ! now rescale each of the 'unrestricted' tracers equally by the fraction of
    ! all 'unrestricted' tracers within the major tracer. This is done on a
    ! gridpoint basis.
    IF (s_tracers(i) == unrestricted)                                          &

    ! Rescale all unrestricted tracers equally if their total is greater than
    ! that in the transported tracer after the mixed_compound and source_gas
    ! tracers have been removed
          WHERE (total_nochange > 1.0)                                         &
            tracers(:,:,:,iposition(i)) =                                      &
            tracers(:,:,:,iposition(i)) / total_nochange
  END DO

  ! now recalculate the major (transported) tracer by removing the sum of all
  ! the 'unrestricted' tracers. If this makes the major tracer drop below rafeps
  ! then the major tracer is set to rafeps, again done on a gridpoint basis.
  ! Calculate new true value of major tracer (held within transported tracer)
  tracers(:,:,:,major) = MAX(tracers(:,:,:,major) -                            &
    c_major * total_all,rafeps)

ELSE

  !-----------------------------------------------------------------------------
  ! lumping (as called at end of UKCA_MAIN1) - add up all tracers
  ! (in VMR of the compound being added into) and add into transported
  ! tracer to give MMR of the major tracer (as used in transport scheme)
  !-----------------------------------------------------------------------------
  tracers(:,:,:,major) = tracers(:,:,:,major) +                                &
     c_major * (total_all + total_source + total_nochange)
END IF

DEALLOCATE (total_all)
DEALLOCATE (total_source)
DEALLOCATE (total_nochange)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE transform

SUBROUTINE ukca_setup_lumping()

! ************************************************************
! Description: set up arrays for each family used in the
!              transform routine
! Method:
! On the first call, loop over all the tracers and set up information
! to be used for lumping. Note that a tracer may belong in no
! family, one family or more than one family.
!
! For each family to which a tracer begins:
! a. increment the n_fam_tracers integer
! b. put the conversion factor from elemental mole fraction to MMR
!    in the c_fam_tracers array. If there are two atoms of the
!    element in this molecule the conversion factor is the conversion
!    VMR to MMR divided by two.
! c. but the index of this species in the tracer structure into
!    the fam_tracers array
! d. if the type of tracer is not 'unrestricted' (the default)
!    then set the type in the s_fam_tracers array to the appropriate type.
!
! NOTE: the 'major' tracer (i.e. the one which is used to lump the species
!       into) is not given an n_SPECIES_tracers value, instead it is given
!       a singular n_SPECIES integer
!
!       Not all the species below may be active tracer for a given model
!       configuration
!
!       If additional N/Cl/Br species are added as tracers then they also
!       need to be added here, and the initial condition of the lumped
!       species would need to be re-calculated accordingly
!
! ************************************************************

USE asad_mod,       ONLY: advt, jpctr

USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper

USE parkind1,       ONLY: jprb, jpim
USE yomhook,        ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SETUP_LUMPING'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i=1,jpctr
  SELECT CASE (advt(i))

    ! ************************************************************
    ! Chlorine tracers
    ! ************************************************************
  CASE ('HCl       ')
    ! MAJOR TRACER FOR Cl
    n_hcl  = i
  CASE ('Cl        ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_cl
    cl_tracers(n_cl_tracers) = i
  CASE ('ClO       ','Clx       ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_clo
    cl_tracers(n_cl_tracers) = i
  CASE ('Cl2O2     ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_cl2o2/2.0
    cl_tracers(n_cl_tracers) = i
  CASE ('HOCl      ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_hocl
    cl_tracers(n_cl_tracers) = i
  CASE ('OClO      ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_oclo
    cl_tracers(n_cl_tracers) = i
  CASE ('ClONO2    ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_clono2
    cl_tracers(n_cl_tracers) = i
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_clono2
    n_tracers(n_n_tracers) = i
    s_n_tracers(n_n_tracers) = mixed_compound
  CASE ('CFCl3     ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_cfcl3/3.0
    cl_tracers(n_cl_tracers) = i
    s_cl_tracers(n_cl_tracers) = source_gas
  CASE ('CF2Cl2    ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_cf2cl2/2.0
    cl_tracers(n_cl_tracers) = i
    s_cl_tracers(n_cl_tracers) = source_gas
  CASE ('MeCl      ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_mecl
    cl_tracers(n_cl_tracers) = i
    s_cl_tracers(n_cl_tracers) = source_gas
  CASE ('CCl4      ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_ccl4/4.0
    cl_tracers(n_cl_tracers) = i
    s_cl_tracers(n_cl_tracers) = source_gas
  CASE ('CF2ClCFCl2')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_cf2clcfcl2/3.0
    cl_tracers(n_cl_tracers) = i
    s_cl_tracers(n_cl_tracers) = source_gas
  CASE ('CHF2Cl    ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_chf2cl
    cl_tracers(n_cl_tracers) = i
    s_cl_tracers(n_cl_tracers) = source_gas
  CASE ('MeCCl3    ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_meccl3/3.0
    cl_tracers(n_cl_tracers) = i
    s_cl_tracers(n_cl_tracers) = source_gas
  CASE ('BrCl      ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_brcl
    cl_tracers(n_cl_tracers) = i
    s_cl_tracers(n_cl_tracers) = mixed_compound
    n_br_tracers = n_br_tracers + 1
    c_br_tracers(n_br_tracers) = c_brcl
    br_tracers(n_br_tracers) = i
  CASE ('CF2ClBr   ')
    n_cl_tracers = n_cl_tracers + 1
    c_cl_tracers(n_cl_tracers) = c_cf2clbr
    cl_tracers(n_cl_tracers) = i
    s_cl_tracers(n_cl_tracers) = mixed_compound
    n_br_tracers = n_br_tracers + 1
    c_br_tracers(n_br_tracers) = c_cf2clbr
    br_tracers(n_br_tracers) = i
    s_br_tracers(n_br_tracers) = source_gas

    ! ************************************************************
    ! Bromine tracers
    ! ************************************************************
  CASE ('BrO       ','Brx       ')
    ! Major tracer for Bromine
    n_bro = i
  CASE ('Br        ')
    n_br_tracers = n_br_tracers + 1
    c_br_tracers(n_br_tracers) = c_br
    br_tracers(n_br_tracers) = i
  CASE ('HBr       ')
    n_br_tracers = n_br_tracers + 1
    c_br_tracers(n_br_tracers) = c_hbr
    br_tracers(n_br_tracers) = i
  CASE ('HOBr      ')
    n_br_tracers = n_br_tracers + 1
    c_br_tracers(n_br_tracers) = c_hobr
    br_tracers(n_br_tracers) = i
  CASE ('BrONO2    ')
    n_br_tracers = n_br_tracers + 1
    c_br_tracers(n_br_tracers) = c_brono2
    br_tracers(n_br_tracers) = i
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_brono2
    n_tracers(n_n_tracers) = i
    s_n_tracers(n_n_tracers) = mixed_compound
  CASE ('MeBr      ')
    n_br_tracers = n_br_tracers + 1
    c_br_tracers(n_br_tracers) = c_mebr
    br_tracers(n_br_tracers) = i
    s_br_tracers(n_br_tracers) = source_gas
  CASE ('CF3Br     ')
    n_br_tracers = n_br_tracers + 1
    c_br_tracers(n_br_tracers) = c_cf3br
    br_tracers(n_br_tracers) = i
    s_br_tracers(n_br_tracers) = source_gas
  CASE ('CH2Br2    ')
    n_br_tracers = n_br_tracers + 1
    c_br_tracers(n_br_tracers) = c_ch2br2/2.0
    br_tracers(n_br_tracers) = i
    s_br_tracers(n_br_tracers) = 1

    ! ************************************************************
    ! nitrogen tracers
    ! ************************************************************
  CASE ('N         ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_n
    n_tracers(n_n_tracers) = i
  CASE ('NO        ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_no
    n_tracers(n_n_tracers) = i
  CASE ('NO2       ','NOx       ')
    ! Major tracer for Nitrogen
    n_no2 = i
  CASE ('NO3       ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_no3
    n_tracers(n_n_tracers) = i
  CASE ('N2O5      ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_n2o5/2.0
    n_tracers(n_n_tracers) = i
  CASE ('HO2NO2    ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_ho2no2
    n_tracers(n_n_tracers) = i
  CASE ('HONO2     ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_hono2
    n_tracers(n_n_tracers) = i
  CASE ('HONO      ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_hono
    n_tracers(n_n_tracers) = i
  CASE ('MeONO2    ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_meono2
    n_tracers(n_n_tracers) = i
  CASE ('PAN       ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_pan
    n_tracers(n_n_tracers) = i
  CASE ('PPAN      ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_ppan
    n_tracers(n_n_tracers) = i
  CASE ('ISON      ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_ison
    n_tracers(n_n_tracers) = i
  CASE ('MPAN      ')
    n_n_tracers = n_n_tracers + 1
    c_n_tracers(n_n_tracers) = c_mpan
    n_tracers(n_n_tracers) = i

    ! ************************************************************
  CASE DEFAULT
    ! Many species may not be in any family. This will
    ! only print if doing debugging
    umMessage='Species: '//advt(i)//' not included in CASE'
    IF ( PrintStatus > PrStatus_Oper )                                         &
     CALL umPrint(umMessage,src=ModuleName//':'//RoutineName)
  END SELECT
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE ukca_setup_lumping

END MODULE ukca_transform_halogen_mod

