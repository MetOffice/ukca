! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Calculates water content of each mode given component
!    concentrations (in air) using ZSR and binary molalities
!    evaluated using water activity data from Jacobson,
!    "Fundamentals of Atmospheric Modelling", page 610 Table B
!    ("Water Activity Data" Table). Equations are from Chapter
!    18 ("Chemical Equilibrium and Dissolution Processes"
!    Chapter. Be aware that in some editions the chapter numbers
!    may be different.
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
! Subroutine Interface:
MODULE ukca_water_content_v_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_WATER_CONTENT_V_MOD'

CONTAINS

SUBROUTINE ukca_water_content_v(nv,mask,cl,rh,ions,wc)

!-----------------------------------------------------------------------
!
! Purpose:
! -------
! Calculates water content of each mode given component
! concentrations (in air) using ZSR and binary molalities
! evaluated using water activity data from Jacobson,
! "Fundamentals of Atmospheric Modelling", page 610 Table B
! ("Water Activity Data" Table).
!
! Inputs
! ------
! NV     : Total number of gridboxes in domain
! IONS   : Logical indicating presence of each ion
! CL     : Concentration of each ion (moles/cc air)
! RH     : Relative humidity (fraction)
! MASK   : Logical array for where in domain to do calculation.
!
! Outputs
! -------
! WC     : Water content for aerosol (moles/cm3 of air)
!
! Local Variables
! ---------------
! IC     : Loop variable for cations
! IA     : Loop variable for anions
! AW     : Water activity (local copy of RH fraction)
! CLI    : Internal copy of CL
! CLP    : Ion pair concentrations (moles/cc air)
! MB     : Ion pair solution molalities (moles/cc water)
! N      : Ion stoiciometries for each electrolyte
! Z      : Charge for each ion
! Y      : Coefficients in expressions for binary molalities from
!        : Jacobson page 610 (Table B.10) for each electrolyte
!        : ("Water Activity Data" Table)
! RH_MIN : Lowest rh for which expression is valid
! MOLAL_MAX : Highest molality for which expression is valid.
!
!----------------------------------------------------------------------
USE ukca_mode_setup,   ONLY: ncation, nanion
USE parkind1,          ONLY: jprb, jpim
USE yomhook,           ONLY: lhook, dr_hook
USE ukca_types_mod,    ONLY: log_small
USE ukca_config_specification_mod, ONLY: glomap_config

IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: nv                       ! No of points
LOGICAL(KIND=log_small), INTENT(IN) :: mask(nv)                 ! Domain mask
LOGICAL(KIND=log_small), INTENT(IN) :: ions(nv,-nanion:ncation) ! Ion presence
                                                                ! switches
REAL, INTENT(IN)    :: rh(nv)                   ! Relative humidity
REAL, INTENT(IN)    :: cl(nv,-nanion:ncation)   ! Ion conc. (mol/cm3 of air)
REAL, INTENT(OUT)   :: wc(nv)                   ! water content (mol/cm3 of air)

! Local variables
INTEGER :: i                   ! Loop counter
INTEGER :: j                   ! Loop counter
INTEGER :: ic                  ! Cation loop variable
INTEGER :: ia                  ! Anion  loop variable
INTEGER :: m                   ! Counter
INTEGER :: idx(nv)             ! Index

!WATER ACTIVITY (%RH EXPRESSED AS A FRACTION)
REAL    :: aw(nv)
REAL    :: dum(nv)
!INTERNAL COPY OF CL
REAL    :: cli(nv,-nanion:ncation)
!ION PAIR CONCENTRATIONS
REAL    :: clp(nv,ncation,-nanion:-1)
!ION PAIR BINARY SOLUTION MOLALITIES AT R
REAL    :: mb(nv,ncation,-nanion:-1)
!ION STOICIOMETRIES FOR EACH ELECTROLYTE
REAL    :: n(-nanion:ncation)
!ION CHARGES
REAL    :: z(-nanion:ncation)
REAL    :: y(3,-4:-1,0:7)
REAL    :: rh_min(3,-4:-1)
REAL    :: molal_max(3,-4:-1)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_WATER_CONTENT_V'

DATA (z(i),i=-nanion,ncation)/1.0,1.0,2.0,1.0,0.0,1.0,1.0,1.0/
!                                    Cl NO3 SO4 HSO4 H2O  H NH4  Na

! As this subroutine is in an OMP parallel section,
! the use of the DATA statements cause y to be shared between threads,
! which creates issues when modifying the array outside of a DATA statement.
! Setting it as threadprivate avoids that.

!$OMP THREADPRIVATE(y)

! Set y coefficients for the calculations of solute molalities from
! Jacobson page 610 (Table B.10) for each electrolyte.
! ("Water Activity Data" Table)
! Also set the min rh and max molality for which the expressions are valid.
!     H+ HSO4- (1,-1)
DATA (y(1,-1,j),j=0,7)/                                                        &
      3.0391387536e1,-1.8995058929e2, 9.7428231047e2,                          &
     -3.1680155761e3, 6.1400925314e3,-6.9116348199e3,                          &
      4.1631475226e3,-1.0383424491e3/
DATA rh_min(1,-1),molal_max(1,-1)/0.0e0,30.4e0/

!     2H+ SO42- (1,-2)
DATA (y(1,-2,j),j=0,7)/                                                        &
      3.0391387536e1,-1.8995058929e2, 9.7428231047e2,                          &
     -3.1680155761e3, 6.1400925314e3,-6.9116348199e3,                          &
      4.1631475226e3,-1.0383424491e3/
DATA rh_min(1,-2),molal_max(1,-2)/0.0e0,30.4e0/

!     H+ NO3- (1,-3)
! One of the y co-efficients is incorrect here (j=6)
! Fixed with glomap_config%l_fix_ukca_water_content below.
DATA (y(1,-3,j),j=0,7)/                                                        &
      2.306844303e1,-3.563608869e1,-6.210577919e1,                             &
      5.510176187e2,-1.460055286e3, 1.894467542e3,                             &
     -1.220611402e2, 3.098597737e2/
DATA rh_min(1,-3),molal_max(1,-3)/0.0e0,22.6e0/

!     H+ Cl- (1,-4)
DATA (y(1,-4,j),j=0,7)/                                                        &
      1.874637647e1,-2.052465972e1,-9.485082073e1,                             &
      5.362930715e2,-1.223331346e3, 1.427089861e3,                             &
     -8.344219112e2, 1.90992437e2/
DATA rh_min(1,-4),molal_max(1,-4)/0.0e0,18.5e0/

!     Na+ HSO4- (2,-1)
DATA (y(2,-1,j),j=0,7)/                                                        &
      1.8457001681e2,-1.6147765817e3, 8.444076586e3,                           &
     -2.6813441936e4, 5.0821277356e4,-5.5964847603e4,                          &
      3.2945298603e4,-8.002609678e3/
DATA rh_min(2,-1),molal_max(2,-1)/1.9e0,158.0e0/

!     2Na+ SO42- (2,-2)
DATA (y(2,-2,j),j=0,7)/                                                        &
      5.5983158e2,-2.56942664e3,4.47450201e3,                                  &
     -3.45021842e3, 9.8527913e2, 0.0e0,                                        &
      0.0e0,        0.0e0/
DATA rh_min(2,-2),molal_max(2,-2)/58.0e0,13.1e0/

!     Na+ NO3- (2,-3)
DATA (y(2,-3,j),j=0,7)/                                                        &
      3.10221762e2,-1.82975944e3, 5.13445395e3,                                &
     -8.01200018e3, 7.07630664e3,-3.33365806e3,                                &
      6.5442029e2,  0.0e0/
DATA rh_min(2,-3),molal_max(2,-3)/30.0e0,56.8e0/

!     Na+ Cl- (2,-4)
DATA (y(2,-4,j),j=0,7)/                                                        &
      5.875248e1,-1.8781997e2, 2.7211377e2,                                    &
     -1.8458287e2, 4.153689e1, 0.0e0,                                          &
      0.0e0,       0.0e0/
DATA rh_min(2,-4),molal_max(2,-4)/47.0e0,13.5e0/

!     NH4+ HSO4- (3,-1)
DATA (y(3,-1,j),j=0,7)/                                                        &
      2.9997156464e2,-2.8936374637e3, 1.4959985537e4,                          &
     -4.5185935292e4, 8.110895603e4, -8.4994863218e4,                          &
      4.7928255412e4,-1.1223105556e4/
DATA rh_min(3,-1),molal_max(3,-1)/6.5e0,165.0e0/

!     2NH4+ SO42- (3,-2)
DATA (y(3,-2,j),j=0,7)/                                                        &
      1.1065495e2,-3.6759197e2, 5.0462934e2,                                   &
     -3.1543839e2, 6.770824e1,  0.0e0,                                         &
      0.0e0,       0.0e0/
DATA rh_min(3,-2),molal_max(3,-2)/37.0e0,29.0e0/

!     NH4+ NO3- (3,-3)
DATA (y(3,-3,j),j=0,7)/                                                        &
      3.983916445e3, 1.153123266e4,-2.13956707e5,                              &
      7.926990533e5,-1.407853405e6, 1.351250086e6,                             &
     -6.770046795e5, 1.393507324e5/
DATA rh_min(3,-3),molal_max(3,-3)/62.0e0,28.0e0/

!     NH4+ Cl- (3,-4)
DATA (y(3,-4,j),j=0,7)/                                                        &
     -7.110541604e3, 7.217772665e4,-3.071054075e5,                             &
      7.144764216e5,-9.840230371e5, 8.03407288e5,                              &
     -3.603924022e5, 6.856992393e4/
DATA rh_min(3,-4),molal_max(3,-4)/47.0e0,23.2e0/

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Put the correct value into the y array (N.B., can't do with DATA statment).
IF (glomap_config%l_fix_ukca_water_content) y(1,-3,6) = -1.220611402e3

m = 0
DO i=1,nv
  IF (mask(i)) THEN
    m = m+1
    idx(m) = i
  END IF
END DO

! Write cl to internal variable to be adjusted in this subroutine only
DO i=-nanion,ncation
  cli(:m,i)=cl(idx(:m),i)
END DO
! Calculate mole concentrations of hypothetical ion pairs
DO ic=1,ncation
  DO ia=-nanion,-1
    ! ..Calculate stoichiometries for each ion pair
    n(ic)=z(ia)
    n(ia)=z(ic)
    IF (ABS(z(ic) - z(ia)) < EPSILON(0.0) .AND.                                &
        ABS(z(ia) - 1.0) > EPSILON(0.0)) THEN
      n(ic)=n(ic)/z(ic)
      n(ia)=n(ia)/z(ia)
    END IF
    WHERE (ions(idx(:m),ia) .AND. ions(idx(:m),ic))
      ! .. Calculate minimum ion pair concentration and subtract from
      ! .. Ion concentration. Eqn. 18.72 of Jacobson.
      clp(:m,ic,ia)=MIN(cli(:m,ic)/n(ic),cli(:m,ia)/n(ia))
      cli(:m,ic)=cli(:m,ic)-n(ic)*clp(:m,ic,ia)
      cli(:m,ia)=cli(:m,ia)-n(ia)*clp(:m,ic,ia)
    END WHERE
  END DO
END DO
! Calculate binary electrolyte molalities at given aw Eqn. 18.66 of
! Jacobson.
IF (glomap_config%l_fix_ukca_water_content) THEN
  DO ic=1,ncation
    DO ia=-nanion,-1
      !If the water content bug is being fixed then copy the
      !original rh values into aw each loop since they may get erroneously
      !overwritten each time around the loop.
      !We want to apply a different min aw for each anion/cation pair:-
      aw(:m)=rh(idx(:m))
      !Prevent aw from going below rh_min (from Table B10 of Jacobson)
      !for the anion/cation pair
      WHERE (aw(:m) < rh_min(ic,ia)/1.0e2)
        aw(:m)=rh_min(ic,ia)/1.0e2
      END WHERE
      mb(:m,ic,ia)=0.0
      WHERE (ions(idx(:m),ia) .AND. ions(idx(:m),ic))
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,0)*aw(:m)**0
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,1)*aw(:m)**1
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,2)*aw(:m)**2
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,3)*aw(:m)**3
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,4)*aw(:m)**4
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,5)*aw(:m)**5
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,6)*aw(:m)**6
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,7)*aw(:m)**7
        mb(:m,ic,ia)=MIN(mb(:m,ic,ia),molal_max(ic,ia))
      END WHERE
    END DO
  END DO
ELSE
  ! Copy fractional relative humidity to water activity (local)
  aw(:m)=rh(idx(:m))
  DO ic=1,ncation
    DO ia=-nanion,-1
      !Prevent aw from going below rh_min (from Table B10 of Jacobson)
      !for the anion/cation pair
      WHERE (aw(:m) < rh_min(ic,ia)/1.0e2)
        aw(:m)=rh_min(ic,ia)/1.0e2
      END WHERE
      mb(:m,ic,ia)=0.0
      WHERE (ions(idx(:m),ia) .AND. ions(idx(:m),ic))
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,0)*aw(:m)**0
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,1)*aw(:m)**1
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,2)*aw(:m)**2
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,3)*aw(:m)**3
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,4)*aw(:m)**4
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,5)*aw(:m)**5
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,6)*aw(:m)**6
        mb(:m,ic,ia)=mb(:m,ic,ia)+y(ic,ia,7)*aw(:m)**7
        mb(:m,ic,ia)=MIN(mb(:m,ic,ia),molal_max(ic,ia))
      END WHERE
    END DO
  END DO
END IF
! Calculate water content (mol/cm3 air); Eqn. 18.71 of Jacobson.
dum(:m)=0.0
DO ic=1,ncation
  DO ia=-nanion,-1
    WHERE (ions(idx(:m),ia) .AND. ions(idx(:m),ic))
      dum(:m)=dum(:m)+clp(:m,ic,ia)/mb(:m,ic,ia)
    END WHERE
  END DO
END DO
wc(idx(:m))=(1.0/18.0e-3)*dum(:m)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_water_content_v
END MODULE ukca_water_content_v_mod
