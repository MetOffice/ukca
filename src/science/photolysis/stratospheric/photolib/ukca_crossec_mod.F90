! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
MODULE ukca_crossec_mod
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!
! Description:
!  Absorption cross sections (UKCA)
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
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.
!

REAL :: accl4 (jpwav)
REAL :: ach2o (jpwav)
REAL :: acl2o2(jpwav)
REAL :: amcfm (jpwav)
REAL :: acnita(jpwav)
REAL :: acnitb(jpwav)
REAL :: aco2  (jpwav)
REAL :: af11  (jpwav)
REAL :: af113 (jpwav)
REAL :: af12  (jpwav)
REAL :: af22  (jpwav)
REAL :: ah2o  (jpwav)
REAL :: ah2o2 (jpwav)
REAL :: ahcl  (jpwav)
REAL :: ahno2 (jpwav)
REAL :: ahno3 (jpwav)
REAL :: ahocl (jpwav)
REAL :: ach3cl(jpwav)
REAL :: an2o  (jpwav)
REAL :: an2o5 (jpwav)
REAL :: ano   (jpwav)
REAL :: ano2  (jpwav)
REAL :: ano31 (jpwav)
REAL :: ano32 (jpwav)
REAL :: ao2   (jpwav)
REAL :: ao2a  (jpwav)
REAL :: ao2b  (jpwav)
REAL :: abrcl (jpwav)
REAL :: abrno3(jpwav)
REAL :: abro  (jpwav)
REAL :: ahobr (jpwav)
REAL :: aoclo (jpwav)
REAL :: ac2oa (jpwav)
REAL :: ac2ob (jpwav)
REAL :: amhp  (jpwav)
REAL :: ao2sr (jpwav)
REAL :: ao3   (jpwav)
REAL :: apna  (jpwav)
REAL :: acof2 (jpwav)
REAL :: acofcl(jpwav)
REAL :: ach3br(jpwav)
REAL :: af12b1(jpwav)
REAL :: af13b1(jpwav)
REAL :: scs   (jpwav)
REAL :: qeno2 (jpwav)
REAL :: qeo1d (jpwav)
REAL :: quanta(jpwav)
REAL :: ach4  (jpwav)
REAL :: wavecm(jpwav)
REAL :: wavenm(jpwav)
REAL :: aobro (jpwav)
REAL :: ahono (jpwav)
REAL :: aocs  (jpwav)
REAL :: aso2  (jpwav)
REAL :: amena (jpwav)
REAL :: achbr3(jpwav)
REAL :: adbrm (jpwav)
REAL :: acs2  (jpwav)
REAL :: ah2so4(jpwav)
REAL :: aso3  (jpwav)


END MODULE ukca_crossec_mod
