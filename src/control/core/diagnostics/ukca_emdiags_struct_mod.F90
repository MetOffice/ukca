! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Structure for UKCA emission diagnostics in section 50.
!   This allows easy addition of new emission diagnostics
!   without additional passing of arguments.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE ukca_emdiags_struct_mod

IMPLICIT NONE
TYPE :: emdiags_struct

  ! Flags indicating whether diagnostics are selected
  LOGICAL :: l_em_no
  LOGICAL :: l_em_ch4
  LOGICAL :: l_em_co
  LOGICAL :: l_em_hcho
  LOGICAL :: l_em_c2h6
  LOGICAL :: l_em_c3h8
  LOGICAL :: l_em_me2co
  LOGICAL :: l_em_mecho
  LOGICAL :: l_em_c5h8
  LOGICAL :: l_em_c4h10
  LOGICAL :: l_em_c2h4
  LOGICAL :: l_em_c3h6
  LOGICAL :: l_em_tol
  LOGICAL :: l_em_oxyl
  LOGICAL :: l_em_ch3oh
  LOGICAL :: l_em_h2
  LOGICAL :: l_em_no_air
  LOGICAL :: l_em_montrp
  LOGICAL :: l_em_meoh
  LOGICAL :: l_em_nh3
  LOGICAL :: l_em_dms
  LOGICAL :: l_em_so2low
  LOGICAL :: l_em_so2hi
  LOGICAL :: l_em_so2nat
  ! Logicals for CRI emission species
  LOGICAL :: l_em_etoh
  LOGICAL :: l_em_c2h2
  LOGICAL :: l_em_tbut2ene
  LOGICAL :: l_em_apinene
  LOGICAL :: l_em_bpinene
  LOGICAL :: l_em_benzene
  LOGICAL :: l_em_hcooh
  LOGICAL :: l_em_meco2h
  LOGICAL :: l_em_etcho
  LOGICAL :: l_em_hoch2cho
  LOGICAL :: l_em_mek


  ! Pointers to hold emission diagnostics
  REAL, POINTER :: em_no     (:,:)
  REAL, POINTER :: em_ch4    (:,:)
  REAL, POINTER :: em_co     (:,:)
  REAL, POINTER :: em_hcho   (:,:)
  REAL, POINTER :: em_c2h6   (:,:)
  REAL, POINTER :: em_c3h8   (:,:)
  REAL, POINTER :: em_me2co  (:,:)
  REAL, POINTER :: em_mecho  (:,:)
  REAL, POINTER :: em_c5h8   (:,:)
  REAL, POINTER :: em_c4h10  (:,:)
  REAL, POINTER :: em_c2h4   (:,:)
  REAL, POINTER :: em_c3h6   (:,:)
  REAL, POINTER :: em_tol    (:,:)
  REAL, POINTER :: em_oxyl   (:,:)
  REAL, POINTER :: em_ch3oh  (:,:)
  REAL, POINTER :: em_h2     (:,:)
  REAL, POINTER :: em_no_air (:,:,:)
  REAL, POINTER :: em_montrp (:,:)
  REAL, POINTER :: em_meoh   (:,:)
  REAL, POINTER :: em_nh3    (:,:)
  REAL, POINTER :: em_dms    (:,:)
  REAL, POINTER :: em_so2low (:,:)
  REAL, POINTER :: em_so2hi  (:,:)
  REAL, POINTER :: em_so2nat (:,:,:)
  ! Pointers for CRI emissions
  REAL, POINTER :: em_etoh    (:,:)
  REAL, POINTER :: em_c2h2    (:,:)
  REAL, POINTER :: em_tbut2ene(:,:)
  REAL, POINTER :: em_apinene (:,:)
  REAL, POINTER :: em_bpinene (:,:)
  REAL, POINTER :: em_benzene (:,:)
  REAL, POINTER :: em_hcooh   (:,:)
  REAL, POINTER :: em_meco2h  (:,:)
  REAL, POINTER :: em_etcho   (:,:)
  REAL, POINTER :: em_hoch2cho(:,:)
  REAL, POINTER :: em_mek     (:,:)

END TYPE emdiags_struct

END MODULE ukca_emdiags_struct_mod
