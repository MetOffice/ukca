! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE ukca_tbjs_mod
USE ukca_parpho_mod, ONLY: jplev,jpchi,jptem,jpo3p
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!  Tabulation of photolysis rates

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.

REAL :: tabj2a   (jplev,jpchi,jptem,jpo3p)
REAL :: tabj2b   (jplev,jpchi,jptem,jpo3p)
REAL :: tabj3    (jplev,jpchi,jptem,jpo3p)
REAL :: tabj3a   (jplev,jpchi,jptem,jpo3p)
REAL :: tabjno   (jplev,jpchi,jptem,jpo3p)
REAL :: tabjno31 (jplev,jpchi,jptem,jpo3p)
REAL :: tabjno32 (jplev,jpchi,jptem,jpo3p)
REAL :: tabjno2  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjn2o5 (jplev,jpchi,jptem,jpo3p)
REAL :: tabjhno3 (jplev,jpchi,jptem,jpo3p)
REAL :: tabjh2o  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjh2o2 (jplev,jpchi,jptem,jpo3p)
REAL :: tabjf11  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjf12  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjf22  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjf113 (jplev,jpchi,jptem,jpo3p)
REAL :: tabjch3cl(jplev,jpchi,jptem,jpo3p)
REAL :: tabjccl4 (jplev,jpchi,jptem,jpo3p)
REAL :: tabjcnita(jplev,jpchi,jptem,jpo3p)
REAL :: tabjcnitb(jplev,jpchi,jptem,jpo3p)
REAL :: tabjhcl  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjhocl (jplev,jpchi,jptem,jpo3p)
REAL :: tabjpna  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjcl2o2(jplev,jpchi,jptem,jpo3p)
REAL :: tabjn2o  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjbrcl (jplev,jpchi,jptem,jpo3p)
REAL :: tabjbrno3(jplev,jpchi,jptem,jpo3p)
REAL :: tabjhobr (jplev,jpchi,jptem,jpo3p)
REAL :: tabjbro  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjoclo (jplev,jpchi,jptem,jpo3p)
REAL :: tabjc2oa (jplev,jpchi,jptem,jpo3p)
REAL :: tabjc2ob (jplev,jpchi,jptem,jpo3p)
REAL :: tabjmhp  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjmcfm (jplev,jpchi,jptem,jpo3p)
REAL :: tabjch3br(jplev,jpchi,jptem,jpo3p)
REAL :: tabjf12b1(jplev,jpchi,jptem,jpo3p)
REAL :: tabjf13b1(jplev,jpchi,jptem,jpo3p)
REAL :: tabjcof2 (jplev,jpchi,jptem,jpo3p)
REAL :: tabjcofcl(jplev,jpchi,jptem,jpo3p)
REAL :: tabjch4  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjcos  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjso2  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjco2  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjhono (jplev,jpchi,jptem,jpo3p)
REAL :: tabjmena (jplev,jpchi,jptem,jpo3p)
REAL :: tabjchbr3(jplev,jpchi,jptem,jpo3p)
REAL :: tabjdbrm (jplev,jpchi,jptem,jpo3p)
REAL :: tabjcs2  (jplev,jpchi,jptem,jpo3p)
REAL :: tabjh2so4(jplev,jpchi,jptem,jpo3p)
REAL :: tabjso3  (jplev,jpchi,jptem,jpo3p)
REAL :: sao3c    (jplev)
REAL :: tabpres  (jplev)
REAL :: tabang   (jpchi)
REAL :: tabt     (jptem)
REAL :: tabo3    (jpo3p)

END MODULE ukca_tbjs_mod
