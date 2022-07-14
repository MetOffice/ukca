! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
MODULE ukca_parpho_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!
! Description:
!  Parameters for the photolysis scheme
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

! Number of levels in jtable.
INTEGER, PARAMETER :: jplevp1=64
INTEGER, PARAMETER :: jplev=jplevp1-1
!Number of zenith angles in jtable.
INTEGER, PARAMETER :: jpchi=20
!Number of za greater than 90 deg in jtable.
INTEGER, PARAMETER :: jps90=5
!Index of 90 degrees in jtable.
INTEGER, PARAMETER :: jpchin=jpchi-jps90
!Number of wavelength intervals.
INTEGER, PARAMETER :: jpwav=203
!Range of wavelength intervals to use.
INTEGER, PARAMETER :: jplo=46
INTEGER, PARAMETER :: jphi=203
!Number of temperatures in jtable
INTEGER, PARAMETER :: jptem=3
!Number of O3 profiles in jtable
INTEGER, PARAMETER :: jpo3p=5

!Maximum zenith angle in jtable :: degrees
REAL, PARAMETER :: szamax=98.0
!Min and max temperatures
REAL, PARAMETER :: tmin=200.0
REAL, PARAMETER :: tmax=250.0
!Min and max O3 profile factors
REAL, PARAMETER :: o3min=0.3
REAL, PARAMETER :: o3max=2.0
!

END MODULE ukca_parpho_mod
