! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE lymana_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Photodissociation due to Lyman alpha for O2, H2O, CH4.

! Method:
!     Parameterisation based on Nicolet (see Brasseur + Solomon).

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LYMANA_MOD'

CONTAINS
SUBROUTINE lymana(tspo2,tabj2,tabjh2o,tabjch4,quanta,jo)
USE ukca_parpho_mod, ONLY: jplev, jpchi, jpchin,                               &
                           jpwav, jptem, jpo3p

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: jo
REAL, INTENT(IN) :: tspo2  (jplev,jpchi)
REAL, INTENT(IN) :: quanta(jpwav)
REAL, INTENT(IN OUT) :: tabj2  (jplev,jpchi,jptem,jpo3p)
REAL, INTENT(IN OUT) :: tabjh2o(jplev,jpchi,jptem,jpo3p)
REAL, INTENT(IN OUT) :: tabjch4(jplev,jpchi,jptem,jpo3p)

! Local variables
REAL :: expfac
REAL :: quly
REAL :: siglya
REAL :: ply

INTEGER :: jc
INTEGER :: jl
INTEGER :: jt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LYMANA'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Only for zenith angles less than 90 degrees
DO jt=1,jptem
  DO jc=1,jpchin-1
    DO jl=1,jplev

      !       slant column of O2
      ply=tspo2(jl,jc)

      expfac= 4.17e-19*(ply**0.917)
      siglya=expfac/ply
      quly=quanta(1)*EXP(-expfac)

      tabj2  (jl,jc,jt,jo) = siglya*quly
      tabjh2o(jl,jc,jt,jo) = 1.40e-17*0.85*quly
      tabjch4(jl,jc,jt,jo) = 1.37e-17*0.85*quly

    END DO
  END DO
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE lymana
END MODULE lymana_mod
