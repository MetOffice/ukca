! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module containing the subroutine ukca_check_radaer_coupling.
!   This routine is called once to check that all the coupling fields
!   for RADAER required in the all_ntp structure are present.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds, University
! of Oxford, and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_check_radaer_coupling_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName =                           &
                   'UKCA_CHECK_RADAER_COUPLING_MOD'

CONTAINS

SUBROUTINE ukca_check_radaer_coupling(all_ntp)

USE ukca_config_specification_mod, ONLY: glomap_config, glomap_variables

USE ukca_ntp_mod,      ONLY: ntp_type, dim_ntp,                                &
                             name2ntpindex

USE ukca_mode_setup,   ONLY: nmodes,                                           &
                             mode_ait_sol,                                     &
                             mode_acc_sol,                                     &
                             mode_cor_sol,                                     &
                             mode_ait_insol,                                   &
                             mode_acc_insol,                                   &
                             mode_cor_insol,                                   &
                             mode_sup_insol,                                   &
                             cp_su, cp_bc, cp_oc,                              &
                             cp_cl, cp_du, cp_so,                              &
                             cp_no3, cp_nh4, cp_nn

USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,          ONLY: jprb, jpim
USE yomhook,           ONLY: lhook, dr_hook
USE ereport_mod,       ONLY: ereport
USE umPrintMgr, ONLY: umMessage, umPrint

IMPLICIT NONE

! Non transported prognostics
TYPE(ntp_type), INTENT(IN) :: all_ntp(dim_ntp)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
LOGICAL, POINTER :: mode (:)
INTEGER, POINTER :: ncp

INTEGER :: imode       ! loop counter for modes
INTEGER :: icp         ! loop counter for components
INTEGER :: i           ! index of all_ntp array

INTEGER :: errcode     ! error code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName =                                   &
                   'UKCA_CHECK_RADAER_COUPLING'

CHARACTER(LEN=errormessagelength) :: cmessage     ! Error message

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp

IF (glomap_config%l_ukca_radaer) THEN
  errcode = 0
  DO imode = mode_ait_sol,nmodes
    IF (mode(imode)) THEN
      SELECT CASE(imode)
      CASE (mode_ait_sol)
        i = name2ntpindex('drydiam_ait_sol     ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 701
          WRITE(umMessage,'(A40)') 'Error condition for '//'drydiam_ait_sol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
        i = name2ntpindex('wetdiam_ait_sol     ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 702
          WRITE(umMessage,'(A40)') 'Error condition for '//'wetdiam_ait_sol'
          CALL umPrint(umMessage,src=Routinename)
        END IF
        i = name2ntpindex('aerdens_ait_sol     ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 703
          WRITE(umMessage,'(A40)') 'Error condition for '//'aerdens_ait_sol'
          CALL umPrint(umMessage,src=Routinename)
        END IF
        i = name2ntpindex('pvol_h2o_ait_sol    ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 704
          WRITE(umMessage,'(A40)') 'Error condition for '//'pvol_h2o_ait_sol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      CASE (mode_acc_sol)
        i = name2ntpindex('drydiam_acc_sol     ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 705
          WRITE(umMessage,'(A40)') 'Error condition for '//'drydiam_acc_sol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
        i = name2ntpindex('wetdiam_acc_sol     ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 706
          WRITE(umMessage,'(A40)') 'Error condition for '//'wetdiam_acc_sol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
        i = name2ntpindex('aerdens_acc_sol     ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 707
          WRITE(umMessage,'(A40)') 'Error condition for '//'aerdens_acc_sol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
        i = name2ntpindex('pvol_h2o_acc_sol    ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 708
          WRITE(umMessage,'(A40)') 'Error condition for '//'pvol_h2o_acc_sol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      CASE (mode_cor_sol)
        i = name2ntpindex('drydiam_cor_sol     ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 709
          WRITE(umMessage,'(A40)') 'Error condition for '//'drydiam_cor_sol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
        i = name2ntpindex('wetdiam_cor_sol     ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 710
          WRITE(umMessage,'(A40)') 'Error condition for '//'wetdiam_cor_sol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
        i = name2ntpindex('aerdens_cor_sol     ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 711
          WRITE(umMessage,'(A40)') 'Error condition for '//'aerdens_cor_sol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
        i = name2ntpindex('pvol_h2o_cor_sol    ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 712
          WRITE(umMessage,'(A40)') 'Error condition for '//'pvol_h2o_cor_sol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      CASE (mode_ait_insol)
        i = name2ntpindex('drydiam_ait_insol   ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 713
          WRITE(umMessage,'(A40)') 'Error condition for '//                    &
                                                   'drydiam_ait_insol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
        i = name2ntpindex('aerdens_ait_insol   ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 714
          WRITE(umMessage,'(A40)') 'Error condition for '//                    &
                                                   'aerdens_ait_insol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      CASE (mode_acc_insol)
        i = name2ntpindex('drydiam_acc_insol   ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 715
          WRITE(umMessage,'(A40)') 'Error condition for '//                    &
                                                   'drydiam_acc_insol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
        i = name2ntpindex('aerdens_acc_insol   ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 716
          WRITE(umMessage,'(A40)') 'Error condition for '//                    &
                                                   'aerdens_acc_insol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      CASE (mode_cor_insol)
        i = name2ntpindex('drydiam_cor_insol   ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 717
          WRITE(umMessage,'(A40)') 'Error condition for '//                    &
                                                   'drydiam_cor_insol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
        i = name2ntpindex('aerdens_cor_insol   ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 718
          WRITE(umMessage,'(A40)') 'Error condition for '//                    &
                                                   'aerdens_cor_insol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      CASE (mode_sup_insol)
        i = name2ntpindex('drydiam_sup_insol   ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 719
          WRITE(umMessage,'(A40)') 'Error condition for '//                    &
                                                   'drydiam_sup_insol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
        i = name2ntpindex('aerdens_sup_insol   ')
        IF (.NOT. all_ntp(i)%l_required) THEN
          errcode = 720
          WRITE(umMessage,'(A40)') 'Error condition for '//                    &
                                                   'aerdens_sup_insol'
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      CASE DEFAULT
        cmessage=' Mode not found in RADAER coupling CASE statement'
        errcode = ABS(imode)
        CALL ereport(RoutineName,errcode,cmessage)
      END SELECT
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          IF (imode == mode_ait_sol) THEN
            SELECT CASE (icp)
            CASE (cp_su)
              i = name2ntpindex('pvol_su_ait_sol    ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 721
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_su_ait_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_bc)
              i = name2ntpindex('pvol_bc_ait_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 722
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_bc_ait_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_oc)
              i = name2ntpindex('pvol_oc_ait_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 723
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_oc_ait_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_so)
              i = name2ntpindex('pvol_so_ait_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 724
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_so_ait_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_no3)
              i = name2ntpindex('pvol_no3_ait_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 725
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_no3_ait_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_nh4)
              i = name2ntpindex('pvol_nh4_ait_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 726
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_nh4_ait_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE DEFAULT
              cmessage = ' Component not found in RADAER coupling CASE'//      &
                         ' statement'
              errcode = ABS(imode*100) + ABS(icp)
              CALL ereport(RoutineName,errcode,cmessage)
            END SELECT
          ELSE IF (imode == mode_acc_sol) THEN
            SELECT CASE (icp)
            CASE (cp_su)
              i = name2ntpindex('pvol_su_acc_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 727
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_su_acc_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_bc)
              i = name2ntpindex('pvol_bc_acc_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 728
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_bc_acc_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_oc)
              i = name2ntpindex('pvol_oc_acc_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 729
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_oc_acc_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_cl)
              i = name2ntpindex('pvol_ss_acc_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 730
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_ss_acc_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_du)
              i = name2ntpindex('pvol_du_acc_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 731
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_du_acc_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_so)
              i = name2ntpindex('pvol_so_acc_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 732
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_so_acc_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_no3)
              i = name2ntpindex('pvol_no3_acc_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 733
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_no3_acc_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_nh4)
              i = name2ntpindex('pvol_nh4_acc_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 734
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_nh4_acc_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_nn)
              i = name2ntpindex('pvol_nn_acc_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 735
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_nn_acc_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE DEFAULT
              cmessage = ' Component not found in RADAER coupling CASE'//      &
                         ' statement'
              errcode = ABS(imode*100) + ABS(icp)
              CALL ereport(RoutineName,errcode,cmessage)
            END SELECT
          ELSE IF (imode == mode_cor_sol) THEN
            SELECT CASE (icp)
            CASE (cp_su)
              i = name2ntpindex('pvol_su_cor_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 736
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_su_cor_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_bc)
              i = name2ntpindex('pvol_bc_cor_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 737
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_bc_cor_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_oc)
              i = name2ntpindex('pvol_oc_cor_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 738
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_oc_cor_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_cl)
              i = name2ntpindex('pvol_ss_cor_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 739
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_ss_cor_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_du)
              i = name2ntpindex('pvol_du_cor_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 740
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_du_cor_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_so)
              i = name2ntpindex('pvol_so_cor_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 741
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_so_cor_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_no3)
              i = name2ntpindex('pvol_no3_cor_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 742
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_no3_cor_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_nh4)
              i = name2ntpindex('pvol_nh4_cor_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 743
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_nh4_cor_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_nn)
              i = name2ntpindex('pvol_nn_cor_sol     ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 744
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_nn_cor_sol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE DEFAULT
              cmessage = ' Component not found in RADAER coupling CASE'//      &
                         ' statement'
              errcode = ABS(imode*100) + ABS(icp)
              CALL ereport(RoutineName,errcode,cmessage)
            END SELECT
          ELSE IF (imode == mode_ait_insol) THEN
            SELECT CASE (icp)
            CASE (cp_bc)
              i = name2ntpindex('pvol_bc_ait_insol   ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 745
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_bc_ait_insol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE (cp_oc)
              i = name2ntpindex('pvol_oc_ait_insol   ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 746
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_oc_ait_insol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE DEFAULT
              cmessage = ' Component not found in RADAER coupling CASE'//      &
                         ' statement'
              errcode = ABS(imode*100) + ABS(icp)
              CALL ereport(RoutineName,errcode,cmessage)
            END SELECT
          ELSE IF (imode == mode_acc_insol) THEN
            SELECT CASE (icp)
            CASE (cp_du)
              i = name2ntpindex('pvol_du_acc_insol   ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 747
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_du_acc_insol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE DEFAULT
              cmessage = ' Component not found in RADAER coupling CASE'//      &
                         ' statement'
              errcode = ABS(imode*100) + ABS(icp)
              CALL ereport(RoutineName,errcode,cmessage)
            END SELECT
          ELSE IF (imode == mode_cor_insol) THEN
            SELECT CASE (icp)
            CASE (cp_du)
              i = name2ntpindex('pvol_du_cor_insol   ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 748
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_du_cor_insol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE DEFAULT
              cmessage = ' Component not found in RADAER coupling CASE'//      &
                         ' statement'
              errcode = ABS(imode*100) + ABS(icp)
              CALL ereport(RoutineName,errcode,cmessage)
            END SELECT
          ELSE IF (imode == mode_sup_insol) THEN
            SELECT CASE (icp)
            CASE (cp_du)
              i = name2ntpindex('pvol_du_sup_insol   ')
              IF (.NOT. all_ntp(i)%l_required) THEN
                errcode = 749
                WRITE(umMessage,'(A40)') 'Error condition for '//              &
                                                    'pvol_du_sup_insol'
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            CASE DEFAULT
              cmessage = ' Component not found in RADAER coupling CASE'//      &
                         ' statement'
              errcode = ABS(imode*100) + ABS(icp)
              CALL ereport(RoutineName,errcode,cmessage)
            END SELECT
          ELSE
            cmessage = ' mode out of range in RADAER coupling IF clause'
            errcode = ABS(imode)
            CALL ereport(RoutineName,errcode,cmessage)
          END IF        ! imode == ?
        END IF         ! component
      END DO    ! icp

    END IF    ! mode(imode)
  END DO  ! imode

  IF ( errcode /= 0) THEN
    cmessage = ' Element of all_ntp array uninitialised'
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF   ! l_ukca_radaer

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_check_radaer_coupling
END MODULE ukca_check_radaer_coupling_mod
