! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Calculates NaNO3 emissions on dust and sea-salt
!
!  References:
!    Hauglustaine et al ACP 2014
!    Code from Samuel Remy (samuel.remy at lmd.jussieu.fr)
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and The Met Office.
!  See:  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ######################################################################

MODULE ukca_aer_no3_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_AER_NO3_MOD'

CONTAINS

SUBROUTINE aer_no3_2bindu (row_length, rows, nlevs, ptsphy, nssmodes,          &
                           nno3modes, pmdust1, pmdust2, pmss, pnss, pdss,      &
                           phno3, pnano3, ptroplev, prhcl, pt, prho,           &
                           ! INPUT above , OUTPUT below
                           ptss, pthno3, ptno3)

USE ukca_um_legacy_mod,   ONLY: drep, rho_dust => rhop
USE ukca_constants,       ONLY: rpi => pi
USE parkind1,             ONLY: jpim, jprb      ! DrHook
USE yomhook,              ONLY: lhook, dr_hook  ! DrHook
IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: nlevs      ! number of vertical levels
REAL, INTENT(IN)    :: ptsphy     ! Timestep length (s)
INTEGER, INTENT(IN) :: nssmodes   ! No. sea-salt modes
INTEGER, INTENT(IN) :: nno3modes  ! No. nitrate modes

REAL, INTENT(IN)    :: pmss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt modal MMR
REAL, INTENT(IN)    :: pnss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt modal number concentration (m-3)
REAL, INTENT(IN)    :: pdss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt modal diameter (m)
REAL, INTENT(IN)    :: phno3(1:row_length,1:rows,1:nlevs) ! HNO3 MMR
REAL, INTENT(IN)    :: pnano3(1:row_length,1:rows,1:nlevs)
                       ! Total NaNO3 MMR (kg/kg)
REAL, INTENT(IN)    :: prho(1:row_length,1:rows,1:nlevs)
                       ! Air density (kg m-3)
REAL, INTENT(IN)    :: prhcl(1:row_length,1:rows,1:nlevs)
                       ! Clear-sky fraction relative humidity (%)
REAL, INTENT(IN)    :: pt(1:row_length,1:rows,1:nlevs)
                       ! Temperature on theta levels (K)
INTEGER, INTENT(IN) :: ptroplev(1:row_length,1:rows)
                       ! Tropoause pressure (Pa)

REAL, INTENT(IN)    :: pmdust1(1:row_length,1:rows,1:nlevs)  ! Dust MMR 1
REAL, INTENT(IN)    :: pmdust2(1:row_length,1:rows,1:nlevs)  ! Dust MMR 2
REAL, INTENT(IN OUT)   :: ptss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt MMR tendency (kg/kg/s)
REAL, INTENT(IN OUT)   :: pthno3(1:row_length,1:rows,1:nlevs)
                       ! HNO3 MMR tendency (kg/kg/s)
REAL, INTENT(IN OUT)   :: ptno3(1:row_length,1:rows,1:nlevs,1:nno3modes)
                       ! NO3 MMR tendency (kg/kg/s)

! Local variables
INTEGER :: irh(1:row_length,1:rows,1:nlevs)
INTEGER :: iirh
INTEGER :: ji, jj, jk, jtab, jss, jdu
REAL    :: rssgrowth_rhtab(12) = [ 1.0, 1.0, 1.0, 1.0, 1.442, 1.555,           &
                                 1.666, 1.799, 1.988, 2.131, 2.361, 2.876]
REAL    :: rrhtab(12) = [ 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, &
                        85.0, 90.0, 95.0 ]
REAL    :: zdg0, zdghno3, zntotdust1, zntotdust2, zmfp, zrh1580
REAL    :: zghno3_dd, zghno3_ss
REAL    :: zkn_dd, zkhno3_dd(1:2)
REAL    :: zrhloc, zcino3, zcavmr, zhno3_test, zthno3_tmp
REAL    :: zaccno3, zcorno3, dtno3
REAL    :: zkhno3_ss(1:nssmodes), zkn_ss, zfracm1, zfracm2
REAL    :: zdd1, zdd2, zss1, zss2, sstmp, ddtmp, dtno3old, mcf
REAL, PARAMETER :: zmwnacl  = 58.44e-3
REAL, PARAMETER :: zmwhno3 = 63.0e-3
REAL, PARAMETER :: zmwnano3  = 84.0e-3
REAL, PARAMETER :: zmwcano3_2  = 164.0e-3
REAL, PARAMETER :: zmwca = 40.0e-3
REAL, PARAMETER :: zratiohno3  = (63.0e-3+29e-3) / 63.0e-3
REAL, PARAMETER :: zghno315_dd = 1.0e-5
REAL, PARAMETER :: zghno330_dd = 1.0e-4
REAL, PARAMETER :: zghno370_dd = 6.0e-4
REAL, PARAMETER :: zghno380_dd = 1.05e-3
REAL, PARAMETER :: zghno315_ss = 1.0e-3
REAL, PARAMETER :: zghno330_ss = 1.0e-2
REAL, PARAMETER :: zghno370_ss = 6.0e-2
REAL, PARAMETER :: zghno380_ss = 1.05e-1
REAL, PARAMETER :: zrgas = 8.314               !J/mol/K
REAL, PARAMETER :: zmgas = 29.0e-3              !Molecular weight in Kg
REAL, PARAMETER :: zdmolec = 4.5e-10           !molec diameter in m
REAL, PARAMETER :: zavg = 6.02217e+23          ! Avogadro number [mol-1]
REAL, PARAMETER :: pfrac_ca=0.05               ! Fraction of dust that is Ca2+

! DUST PARAMETERS
! Now using the geometric means for the dust diameters
! rather than the volume medians
REAL, PARAMETER :: drep_N_mean(2) = [ 0.299287e-6, 5.806452e-6 ]


INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='AER_NO3_2BINDU'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ptss(:,:,:,:) = 0.0
pthno3(:,:,:) = 0.0
ptno3(:,:,:,:) = 0.0
! Initialise arrays
DO jk = 1,nlevs
  DO jj = 1,rows
    DO ji = 1,row_length
      irh(ji,jj,jk)=1
      DO jtab=1,12
        IF (prhcl(ji,jj,jk)*100.0 > rrhtab(jtab)) THEN
          irh(ji,jj,jk)=jtab
        END IF
      END DO
    END DO
  END DO
END DO

DO jk = 1,nlevs
  DO jj = 1,rows
    DO ji = 1,row_length
      IF (jk <= ptroplev(ji,jj)) THEN
        iirh=irh(ji,jj,jk)

        zrhloc = prhcl(ji,jj,jk)
        IF (zrhloc < 0.0 ) zrhloc = 0.0
        IF (zrhloc > 0.98 ) zrhloc = 0.98

        ! needed for HNO3 update on dust and ss
        zdg0= 3.0/(8.0*zavg*prho(ji,jj,jk)*(zdmolec**2.0))
        zdghno3 = zdg0 * (zrgas*zmgas/2.0/rpi*pt(ji,jj,jk)*zratiohno3)**0.5

        ! Number of dust particles: mass/size
        zntotdust1 = MAX(pmdust1(ji,jj,jk),0.0)*prho(ji,jj,jk)  /              &
                     (rho_dust*4.0*rpi*(drep(1)*0.5)**3/3.0)
        zntotdust2 = MAX(pmdust2(ji,jj,jk),0.0)*prho(ji,jj,jk)  /              &
                     (rho_dust*4.0*rpi*(drep(2)*0.5)**3/3.0)

        !RH dependent HNO3 uptake coefficient based on Fairlie et al. 2010
        !and scale for SS
        zrh1580 = MAX(MIN(zrhloc,0.80),0.15)
        IF ( zrh1580 < 0.30 ) THEN
          zghno3_dd = zghno315_dd + (zghno330_dd-zghno315_dd)/                 &
                      0.15*(zrh1580-0.15)
          zghno3_ss = zghno315_ss + (zghno330_ss-zghno315_ss)/                 &
                      0.15*(zrh1580-0.15)
        ELSE IF (zrh1580 > 0.30 .AND. zrh1580 < 0.70) THEN
          zghno3_dd = zghno330_dd + (zghno370_dd-zghno330_dd)/                 &
                      0.40*(zrh1580-0.30)
          zghno3_ss = zghno330_ss + (zghno370_ss-zghno330_ss)/                 &
                      0.40*(zrh1580-0.30)
        ELSE IF (zrh1580 > 0.70 ) THEN
          zghno3_dd = zghno370_dd + (zghno380_dd-zghno370_dd)/                 &
                      0.10*(zrh1580-0.70)
          zghno3_ss = zghno370_ss + (zghno380_ss-zghno370_ss)/                 &
                      0.10*(zrh1580-0.70)
        END IF

        ! MFP parameter
        zmfp=3.0 * zdghno3 / SQRT(8.0*zrgas*pt(ji,jj,jk)/rpi/zmwhno3)

        ! Uptake of HNO3 on DUST
        DO jdu = 1,2
          zkn_dd=2.0 * zmfp / drep_N_mean(jdu)
          zkhno3_dd(jdu) = (2.0*rpi*drep_N_mean(jdu)*zdghno3) /                &
                           (1.0+(4.0*zkn_dd/zghno3_dd/3.0)*(1.0-0.47*zghno3_dd/&
                           (1.0+zkn_dd)))
        END DO

        ! Uptake of HNO3 on SEA-SALT
        DO jss = 1,nssmodes
          zkn_ss = 2.0 * zmfp / (pdss(ji,jj,jk,jss) * rssgrowth_rhtab(iirh))
          zkhno3_ss(jss) = (2.0*rpi*pdss(ji,jj,jk,jss)*rssgrowth_rhtab(iirh)*  &
                            zdghno3) /(1.0+(4.0*zkn_ss/zghno3_ss/3.0)*         &
                           (1.0-0.47*zghno3_ss/(1.0+zkn_ss)))
        END DO

        ! Calcium limitation for HNO3 condensation on dust
        zcavmr =(pmdust1(ji,jj,jk)+pmdust2(ji,jj,jk))*pfrac_ca/zmwca*zmgas*0.5
                 !Moles of Ca available for Ca(NO3)2 formation
        zcino3 = pnano3(ji,jj,jk)/zmwnano3*zmgas
        IF (zcino3 > zcavmr)  zkhno3_dd(:) = 0.0

        ! Check to ensure that dust and sea-salt tendencies do not
        ! result in negative aerosol mass
        ddtmp = MAX(0.0,(zntotdust1*zmwca/zmwhno3*phno3(ji,jj,jk)*ptsphy))
        IF ((pmdust1(ji,jj,jk)*pfrac_ca-ddtmp*zkhno3_dd(1)< 0.0) .AND.         &
            (ddtmp>0.0))                                                       &
            zkhno3_dd(1) = pmdust1(ji,jj,jk)*pfrac_ca/ddtmp

        ddtmp = MAX(0.0,(zntotdust2*zmwca/zmwhno3*phno3(ji,jj,jk)*ptsphy))
        IF ((pmdust2(ji,jj,jk)*pfrac_ca-ddtmp*zkhno3_dd(2)< 0.0) .AND.         &
            (ddtmp>0.0))                                                       &
            zkhno3_dd(2) = pmdust2(ji,jj,jk)*pfrac_ca/ddtmp

        DO jss = 1,nssmodes
          sstmp = MAX(0.0,(pnss(ji,jj,jk,jss)*zmwnacl/zmwhno3*                 &
                          phno3(ji,jj,jk)*ptsphy))
          IF ((pmss(ji,jj,jk,jss)-sstmp*zkhno3_ss(jss)<0.0) .AND. (sstmp>0.0)) &
            zkhno3_ss(jss) = pmss(ji,jj,jk,jss)/sstmp
        END DO

        ! Divide NO3 tendencies into accum and coarse modes
        zdd1 = MAX(0.0,zkhno3_dd(1)*zntotdust1*0.145)*zmwcano3_2/zmwhno3
        zdd2 = MAX(0.0,zkhno3_dd(1)*zntotdust1*0.855+zkhno3_dd(2)*zntotdust2)* &
                      zmwcano3_2/zmwhno3
        zss1 = MAX(0.0,zkhno3_ss(1)*pnss(ji,jj,jk,1))*zmwnano3/zmwhno3
        zss2 = MAX(0.0,zkhno3_ss(2)*pnss(ji,jj,jk,2))*zmwnano3/zmwhno3
        zaccno3 = zdd1 + zss1
        zcorno3 = zdd2 + zss2

        ! If the NO3 tendency is greater than 0 then store in output array
        IF ((zaccno3+zcorno3) > 0.0) THEN
          zfracm1 = zaccno3 / (zaccno3+zcorno3)
          zfracm2 = zcorno3 / (zaccno3+zcorno3)

          ! Update NO3 and HNO3 tendencies
          dtno3=(zaccno3+zcorno3)*MAX(0.0,phno3(ji,jj,jk))
          pthno3(ji,jj,jk) = -1.0*MAX(0.0,phno3(ji,jj,jk))*                    &
                             ((zdd1+zdd2)*zmwhno3*2.0/zmwcano3_2 +             &
                              (zss1+zss2)*zmwhno3/zmwnano3)

          ! Check for negative HNO3
          zhno3_test = phno3(ji,jj,jk) + pthno3(ji,jj,jk)*ptsphy
          mcf=1.0
          IF (zhno3_test < 0.0) THEN
            dtno3old=dtno3
            zthno3_tmp = pthno3(ji,jj,jk)
            pthno3(ji,jj,jk) = -1.0*phno3(ji,jj,jk)/ptsphy
            zthno3_tmp = pthno3(ji,jj,jk)-zthno3_tmp
            dtno3 = dtno3-zthno3_tmp*zmwnano3/zmwhno3
            mcf=dtno3/dtno3old
          END IF
          ptno3(ji,jj,jk,1) = dtno3 * zfracm1 ! Accum NO3 emissions
          ptno3(ji,jj,jk,2) = dtno3 * zfracm2 ! Coarse NO3 emissions

          ! Update sea-salt tendencies
          DO jss = 1,nssmodes
            ptss(ji,jj,jk,jss)= -1.0 * zkhno3_ss(jss) * pnss(ji,jj,jk,jss) *   &
                                phno3(ji,jj,jk)* zmwnacl / zmwhno3 * mcf
          END DO

        END IF  ! Reaction rates are greater than zero
      END IF ! < tropopause
    END DO !ji
  END DO !jj
END DO !jk


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE aer_no3_2bindu

!-----------------------------------------------------------------------

SUBROUTINE aer_no3_6bindu (row_length, rows, nlevs, ptsphy, nssmodes,          &
                           nno3modes, pmdust2, pmdust3, pmdust4,               &
                           pmdust5, pmdust6, pmss, pnss, pdss,                 &
                           phno3, pnano3, ptroplev, prhcl, pt, prho,           &
                           ! INPUT above , OUTPUT below
                           ptss, pthno3, ptno3)

USE ukca_um_legacy_mod,   ONLY: drep, rho_dust => rhop
USE ukca_constants,       ONLY: rpi => pi
USE parkind1,             ONLY: jpim, jprb      ! DrHook
USE yomhook,              ONLY: lhook, dr_hook  ! DrHook
IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: nlevs      ! number of vertical levels
REAL, INTENT(IN)    :: ptsphy     ! Timestep length (s)
INTEGER, INTENT(IN) :: nssmodes   ! No. sea-salt modes
INTEGER, INTENT(IN) :: nno3modes  ! No. nitrate modes

REAL, INTENT(IN)    :: pmss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt modal MMR
REAL, INTENT(IN)    :: pnss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt modal number concentration (m-3)
REAL, INTENT(IN)    :: pdss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt modal diameter (m)
REAL, INTENT(IN)    :: phno3(1:row_length,1:rows,1:nlevs) ! HNO3 MMR
REAL, INTENT(IN)    :: pnano3(1:row_length,1:rows,1:nlevs)
                       ! Total NaNO3 MMR (kg/kg)
REAL, INTENT(IN)    :: prho(1:row_length,1:rows,1:nlevs)
                       ! Air density (kg m-3)
REAL, INTENT(IN)    :: prhcl(1:row_length,1:rows,1:nlevs)
                       ! Clear-sky fraction relative humidity (%)
REAL, INTENT(IN)    :: pt(1:row_length,1:rows,1:nlevs)
                       ! Temperature on theta levels (K)
INTEGER, INTENT(IN) :: ptroplev(1:row_length,1:rows)
                       ! Tropoause pressure (Pa)

REAL, INTENT(IN)    :: pmdust2(1:row_length,1:rows,1:nlevs)  ! Dust MMR 2
REAL, INTENT(IN)    :: pmdust3(1:row_length,1:rows,1:nlevs)  ! Dust MMR 3
REAL, INTENT(IN)    :: pmdust4(1:row_length,1:rows,1:nlevs)  ! Dust MMR 4
REAL, INTENT(IN)    :: pmdust5(1:row_length,1:rows,1:nlevs)  ! Dust MMR 5
REAL, INTENT(IN)    :: pmdust6(1:row_length,1:rows,1:nlevs)  ! Dust MMR 6
REAL, INTENT(IN OUT)   :: ptss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt MMR tendency (kg/kg/s)
REAL, INTENT(IN OUT)   :: pthno3(1:row_length,1:rows,1:nlevs)
                       ! HNO3 MMR tendency (kg/kg/s)
REAL, INTENT(IN OUT)   :: ptno3(1:row_length,1:rows,1:nlevs,1:nno3modes)
                       ! NO3 MMR tendency (kg/kg/s)

! Local variables
INTEGER :: irh(1:row_length,1:rows,1:nlevs)
INTEGER :: iirh
INTEGER :: ji, jj, jk, jtab, jss, jdu
REAL    :: rssgrowth_rhtab(12) = [ 1.0, 1.0, 1.0, 1.0, 1.442, 1.555,           &
                                 1.666, 1.799, 1.988, 2.131, 2.361, 2.876]
REAL    :: rrhtab(12) = [ 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, &
                        85.0, 90.0, 95.0 ]
REAL    :: zdg0, zdghno3, zmfp, zrh1580
REAL    :: zntotdust2, zntotdust3, zntotdust4, zntotdust5, zntotdust6
REAL    :: zghno3_dd, zghno3_ss
REAL    :: zkn_dd, zkhno3_dd(1:5)
REAL    :: zrhloc, zcino3, zcavmr
REAL    :: zhno3_test, zthno3_tmp
REAL    :: zaccno3, zcorno3, dtno3
REAL    :: zkhno3_ss(1:nssmodes), zkn_ss, zfracm1, zfracm2
REAL    :: zdd1, zdd2, zss1, zss2, sstmp, ddtmp, dtno3old, mcf
REAL, PARAMETER :: zmwnacl  = 58.44e-3
REAL, PARAMETER :: zmwhno3 = 63.0e-3
REAL, PARAMETER :: zmwnano3  = 84.0e-3
REAL, PARAMETER :: zmwcano3_2  = 164.0e-3
REAL, PARAMETER :: zmwca = 40.0e-3
REAL, PARAMETER :: zratiohno3  = (63.0e-3+29e-3) / 63.0e-3
REAL, PARAMETER :: zghno315_dd = 1.0e-5
REAL, PARAMETER :: zghno330_dd = 1.0e-4
REAL, PARAMETER :: zghno370_dd = 6.0e-4
REAL, PARAMETER :: zghno380_dd = 1.05e-3
REAL, PARAMETER :: zghno315_ss = 1.0e-3
REAL, PARAMETER :: zghno330_ss = 1.0e-2
REAL, PARAMETER :: zghno370_ss = 6.0e-2
REAL, PARAMETER :: zghno380_ss = 1.05e-1
REAL, PARAMETER :: zrgas = 8.314               !J/mol/K
REAL, PARAMETER :: zmgas = 29.0e-3              !Molecular weight in Kg
REAL, PARAMETER :: zdmolec = 4.5e-10           !molec diameter in m
REAL, PARAMETER :: zavg = 6.02217e+23          ! Avogadro number [mol-1]
REAL, PARAMETER :: pfrac_ca=0.05               ! Fraction of dust that is Ca2+

! DUST PARAMETERS
! Now using the geometric means for the dust diameters
! rather than the volume medians
REAL, PARAMETER :: drep_N_mean(6) =                                            &
                   [ 0.881697e-07, 0.278817e-06, 0.881697e-06,                 &
                      0.278817e-05, 0.881697e-05, 0.278817e-04 ]

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='AER_NO3_6BINDU'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ptss(:,:,:,:) = 0.0
pthno3(:,:,:) = 0.0
ptno3(:,:,:,:) = 0.0
! Initialise arrays
DO jk = 1,nlevs
  DO jj = 1,rows
    DO ji = 1,row_length
      irh(ji,jj,jk)=1
      DO jtab=1,12
        IF (prhcl(ji,jj,jk)*100.0 > rrhtab(jtab)) THEN
          irh(ji,jj,jk)=jtab
        END IF
      END DO
    END DO
  END DO
END DO

DO jk = 1,nlevs
  DO jj = 1,rows
    DO ji = 1,row_length
      IF (jk <= ptroplev(ji,jj)) THEN
        iirh=irh(ji,jj,jk)

        zrhloc = prhcl(ji,jj,jk)
        IF (zrhloc < 0.0 ) zrhloc = 0.0
        IF (zrhloc > 0.98 ) zrhloc = 0.98

        ! needed for HNO3 update on dust and ss
        zdg0= 3.0/(8.0*zavg*prho(ji,jj,jk)*(zdmolec**2.0))
        zdghno3 = zdg0 * (zrgas*zmgas/2.0/rpi*pt(ji,jj,jk)*zratiohno3)**0.5

        ! Number of dust particles: mass/size
        zntotdust2 = MAX(pmdust2(ji,jj,jk),0.0)*prho(ji,jj,jk)  /              &
                     (rho_dust*4.0*rpi*(drep(2)*0.5)**3/3.0)
        zntotdust3 = MAX(pmdust3(ji,jj,jk),0.0)*prho(ji,jj,jk)  /              &
                     (rho_dust*4.0*rpi*(drep(3)*0.5)**3/3.0)
        zntotdust4 = MAX(pmdust4(ji,jj,jk),0.0)*prho(ji,jj,jk)  /              &
                     (rho_dust*4.0*rpi*(drep(4)*0.5)**3/3.0)
        zntotdust5 = MAX(pmdust5(ji,jj,jk),0.0)*prho(ji,jj,jk)  /              &
                     (rho_dust*4.0*rpi*(drep(5)*0.5)**3/3.0)
        zntotdust6 = MAX(pmdust6(ji,jj,jk),0.0)*prho(ji,jj,jk)  /              &
                     (rho_dust*4.0*rpi*(drep(6)*0.5)**3/3.0)

        !RH dependent HNO3 uptake coefficient based on Fairlie et al. 2010
        !and scale for SS
        zrh1580 = MAX(MIN(zrhloc,0.80),0.15)
        IF ( zrh1580 < 0.30 ) THEN
          zghno3_dd = zghno315_dd + (zghno330_dd-zghno315_dd)/                 &
                      0.15*(zrh1580-0.15)
          zghno3_ss = zghno315_ss + (zghno330_ss-zghno315_ss)/                 &
                      0.15*(zrh1580-0.15)
        ELSE IF (zrh1580 > 0.30 .AND. zrh1580 < 0.70) THEN
          zghno3_dd = zghno330_dd + (zghno370_dd-zghno330_dd)/                 &
                      0.40*(zrh1580-0.30)
          zghno3_ss = zghno330_ss + (zghno370_ss-zghno330_ss)/                 &
                      0.40*(zrh1580-0.30)
        ELSE IF (zrh1580 > 0.70 ) THEN
          zghno3_dd = zghno370_dd + (zghno380_dd-zghno370_dd)/                 &
                      0.10*(zrh1580-0.70)
          zghno3_ss = zghno370_ss + (zghno380_ss-zghno370_ss)/                 &
                      0.10*(zrh1580-0.70)
        END IF

        ! MFP parameter
        zmfp=3.0 * zdghno3 / SQRT(8.0*zrgas*pt(ji,jj,jk)/rpi/zmwhno3)

        ! Uptake of HNO3 on DUST
        DO jdu = 1,5
          zkn_dd=2.0 * zmfp / drep_N_mean(jdu+1)
          zkhno3_dd(jdu) = (2.0*rpi*drep_N_mean(jdu+1)*zdghno3) /              &
                           (1.0+(4.0*zkn_dd/zghno3_dd/3.0)*(1.0-0.47*zghno3_dd/&
                           (1.0+zkn_dd)))
        END DO

        ! Uptake of HNO3 on SEA-SALT
        DO jss = 1,nssmodes
          zkn_ss = 2.0 * zmfp / (pdss(ji,jj,jk,jss) * rssgrowth_rhtab(iirh))
          zkhno3_ss(jss) = (2.0*rpi*pdss(ji,jj,jk,jss)*rssgrowth_rhtab(iirh)*  &
                            zdghno3) /(1.0+(4.0*zkn_ss/zghno3_ss/3.0)*         &
                           (1.0-0.47*zghno3_ss/(1.0+zkn_ss)))
        END DO

        ! Calcium limitation for HNO3 condensation on dust
        zcavmr =(pmdust2(ji,jj,jk)+pmdust3(ji,jj,jk)+pmdust4(ji,jj,jk)+        &
                 pmdust5(ji,jj,jk)+pmdust6(ji,jj,jk))*pfrac_ca/zmwca*zmgas*0.5
                 !Moles of Ca available for Ca(NO3)2 formation
        zcino3 = pnano3(ji,jj,jk)/zmwnano3*zmgas
        IF (zcino3 > zcavmr)  zkhno3_dd(:) = 0.0

        ! Check to ensure that dust and sea-salt tendencies do not
        ! result in negative aerosol mass
        ddtmp=MAX(0.0,(zntotdust2*zmwca/zmwhno3*phno3(ji,jj,jk)*ptsphy))
        IF ((pmdust2(ji,jj,jk)*pfrac_ca-ddtmp*zkhno3_dd(1)< 0.0) .AND.         &
            (ddtmp>0.0))                                                       &
            zkhno3_dd(1) = pmdust2(ji,jj,jk)*pfrac_ca/ddtmp

        ddtmp=MAX(0.0,(zntotdust3*zmwca/zmwhno3*phno3(ji,jj,jk)*ptsphy))
        IF ((pmdust3(ji,jj,jk)*pfrac_ca-ddtmp*zkhno3_dd(2)< 0.0) .AND.         &
            (ddtmp>0.0))                                                       &
            zkhno3_dd(2) = pmdust3(ji,jj,jk)*pfrac_ca/ddtmp

        ddtmp=MAX(0.0,(zntotdust4*zmwca/zmwhno3*phno3(ji,jj,jk)*ptsphy))
        IF ((pmdust4(ji,jj,jk)*pfrac_ca-ddtmp*zkhno3_dd(3)< 0.0) .AND.         &
            (ddtmp>0.0))                                                       &
            zkhno3_dd(3) = pmdust4(ji,jj,jk)*pfrac_ca/ddtmp

        ddtmp=MAX(0.0,(zntotdust5*zmwca/zmwhno3*phno3(ji,jj,jk)*ptsphy))
        IF ((pmdust5(ji,jj,jk)*pfrac_ca-ddtmp*zkhno3_dd(4)< 0.0) .AND.         &
            (ddtmp>0.0))                                                       &
            zkhno3_dd(4) = pmdust5(ji,jj,jk)*pfrac_ca/ddtmp

        ddtmp=MAX(0.0,(zntotdust6*zmwca/zmwhno3*phno3(ji,jj,jk)*ptsphy))
        IF ((pmdust6(ji,jj,jk)*pfrac_ca-ddtmp*zkhno3_dd(5)< 0.0) .AND.         &
            (ddtmp>0.0))                                                       &
          zkhno3_dd(5) = pmdust6(ji,jj,jk)*pfrac_ca/ddtmp

        DO jss = 1,nssmodes
          sstmp = MAX(0.0,(pnss(ji,jj,jk,jss)*zmwnacl/zmwhno3*                 &
                          phno3(ji,jj,jk)*ptsphy))
          IF ((pmss(ji,jj,jk,jss)-sstmp*zkhno3_ss(jss)<0.0) .AND. (sstmp>0.0)) &
            zkhno3_ss(jss) = pmss(ji,jj,jk,jss)/sstmp
        END DO

        ! Divide NO3 tendencies into accum and coarse modes
        zdd1 = MAX(0.0,zkhno3_dd(1)*zntotdust2+zkhno3_dd(2)*zntotdust3*0.5) *  &
                      zmwcano3_2/zmwhno3
        zdd2 = MAX(0.0,zkhno3_dd(2)*zntotdust3*0.5+zkhno3_dd(3)*zntotdust4+    &
                      zkhno3_dd(4)*zntotdust5+zkhno3_dd(5)*zntotdust6) *       &
                      zmwcano3_2/zmwhno3
        zss1 = MAX(0.0,zkhno3_ss(1)*pnss(ji,jj,jk,1))*zmwnano3/zmwhno3
        zss2 = MAX(0.0,zkhno3_ss(2)*pnss(ji,jj,jk,2))*zmwnano3/zmwhno3
        zaccno3 = zdd1 + zss1
        zcorno3 = zdd2 + zss2

        ! If the NO3 tendency is greater than 0 then store in output array
        IF ((zaccno3+zcorno3) > 0.0) THEN
          zfracm1 = zaccno3 / (zaccno3+zcorno3)
          zfracm2 = zcorno3 / (zaccno3+zcorno3)

          ! Update NO3 and HNO3 tendencies
          dtno3=(zaccno3+zcorno3)*MAX(0.0,phno3(ji,jj,jk))
          pthno3(ji,jj,jk) = -1.0*MAX(0.0,phno3(ji,jj,jk))*                    &
                             ((zdd1+zdd2)*zmwhno3*2.0/zmwcano3_2 +             &
                              (zss1+zss2)*zmwhno3/zmwnano3)

          ! Check for negative HNO3
          zhno3_test = phno3(ji,jj,jk) + pthno3(ji,jj,jk)*ptsphy
          mcf=1.0
          IF (zhno3_test < 0.0) THEN
            dtno3old=dtno3
            zthno3_tmp = pthno3(ji,jj,jk)
            pthno3(ji,jj,jk) = -1.0*phno3(ji,jj,jk)/ptsphy
            zthno3_tmp = pthno3(ji,jj,jk)-zthno3_tmp
            dtno3 = dtno3-zthno3_tmp*zmwnano3/zmwhno3
            mcf=dtno3/dtno3old
          END IF
          ptno3(ji,jj,jk,1) = dtno3 * zfracm1 ! Accum NO3 emissions
          ptno3(ji,jj,jk,2) = dtno3 * zfracm2 ! Coarse NO3 emissions

          ! Update sea-salt tendencies
          DO jss = 1,nssmodes
            ptss(ji,jj,jk,jss)= -1.0 * zkhno3_ss(jss) * pnss(ji,jj,jk,jss) *   &
                                phno3(ji,jj,jk)* zmwnacl / zmwhno3 * mcf
          END DO

        END IF  ! Rates are greater than zero
      END IF  ! < tropopause
    END DO !ji
  END DO !jj
END DO !jk

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE aer_no3_6bindu


!-----------------------------------------------------------------------

SUBROUTINE aer_no3_ukcadu (row_length, rows, nlevs, ptsphy, nssmodes,          &
                           nno3modes, ndumodes, pmdu, pndu, pddu, pmss, pnss,  &
                           pdss, phno3, ptroplev, prhcl, pt, prho,             &
                           ! INPUT above , OUTPUT below
                           ptdu, ptss, pthno3, ptno3)

USE ukca_constants,       ONLY: rpi => pi
USE parkind1,             ONLY: jpim, jprb      ! DrHook
USE yomhook,              ONLY: lhook, dr_hook  ! DrHook
IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: nlevs      ! number of vertical levels
REAL, INTENT(IN)    :: ptsphy     ! Timestep length (s)
INTEGER, INTENT(IN) :: nssmodes   ! No. sea-salt modes
INTEGER, INTENT(IN) :: nno3modes  ! No. nitrate modes
INTEGER, INTENT(IN) :: ndumodes   ! No. dust modes

REAL, INTENT(IN)    :: pmdu(1:row_length,1:rows,1:nlevs,1:ndumodes)
                       ! Dust modal MMR
REAL, INTENT(IN)    :: pndu(1:row_length,1:rows,1:nlevs,1:ndumodes)
                       ! Dust modal number concentration (m-3)
REAL, INTENT(IN)    :: pddu(1:row_length,1:rows,1:nlevs,1:ndumodes)
                       ! Dust modal diameter (m)
REAL, INTENT(IN)    :: pmss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt modal MMR
REAL, INTENT(IN)    :: pnss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt modal number concentration (m-3)
REAL, INTENT(IN)    :: pdss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt modal diameter (m)
REAL, INTENT(IN)    :: phno3(1:row_length,1:rows,1:nlevs) ! HNO3 MMR
                       ! Total NaNO3 MMR (kg/kg)
REAL, INTENT(IN)    :: prho(1:row_length,1:rows,1:nlevs)
                       ! Air density (kg m-3)
REAL, INTENT(IN)    :: prhcl(1:row_length,1:rows,1:nlevs)
                       ! Clear-sky fraction relative humidity (%)
REAL, INTENT(IN)    :: pt(1:row_length,1:rows,1:nlevs)
                       ! Temperature on theta levels (K)
INTEGER, INTENT(IN) :: ptroplev(1:row_length,1:rows)
                       ! Tropoause pressure (Pa)

REAL, INTENT(IN OUT)   :: ptdu(1:row_length,1:rows,1:nlevs,1:ndumodes)
                       ! Dust  MMR tendency (kg/kg/s)
REAL, INTENT(IN OUT)   :: ptss(1:row_length,1:rows,1:nlevs,1:nssmodes)
                       ! Sea-salt MMR tendency (kg/kg/s)
REAL, INTENT(IN OUT)   :: pthno3(1:row_length,1:rows,1:nlevs)
                       ! HNO3 MMR tendency (kg/kg/s)
REAL, INTENT(IN OUT)   :: ptno3(1:row_length,1:rows,1:nlevs,1:nno3modes)
                       ! NO3 MMR tendency (kg/kg/s)

! Local variables
INTEGER :: irh(1:row_length,1:rows,1:nlevs)
INTEGER :: iirh
INTEGER :: ji, jj, jk, jss, jdu, jtab
REAL    :: rssgrowth_rhtab(12) = [ 1.0, 1.0, 1.0, 1.0, 1.442, 1.555,           &
                                 1.666, 1.799, 1.988, 2.131, 2.361, 2.876]
REAL    :: rrhtab(12) = [ 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, &
                        85.0, 90.0, 95.0 ]
REAL    :: zdg0, zdghno3, zmfp, zrh1580
REAL    :: zghno3_dd, zghno3_ss, sumdu
REAL    :: zkn_dd, zhno3_test, zthno3_tmp
REAL    :: zkn_ss, zrhloc
REAL    :: zkhno3_ss(1:nssmodes)
REAL    :: zkhno3_dd(1:ndumodes)
REAL    :: zaccno3, zcorno3, dtno3, zfracm1, zfracm2
REAL    :: zdd1, zdd2, zss1, zss2, sstmp, ddtmp, dtno3old, mcf
REAL, PARAMETER :: zmwnacl  = 58.44e-3
REAL, PARAMETER :: zmwhno3 = 63.0e-3
REAL, PARAMETER :: zmwnano3  = 84.0e-3
REAL, PARAMETER :: zmwcano3_2  = 164.0e-3
REAL, PARAMETER :: zmwcaco3  = 100.0e-3
REAL, PARAMETER :: zmwca = 40.0e-3
REAL, PARAMETER :: zratiohno3  = (63.0e-3+29e-3) / 63.0e-3
REAL, PARAMETER :: zghno315_dd = 1.0e-5
REAL, PARAMETER :: zghno330_dd = 1.0e-4
REAL, PARAMETER :: zghno370_dd = 6.0e-4
REAL, PARAMETER :: zghno380_dd = 1.05e-3
REAL, PARAMETER :: zghno315_ss = 1.0e-3
REAL, PARAMETER :: zghno330_ss = 1.0e-2
REAL, PARAMETER :: zghno370_ss = 6.0e-2
REAL, PARAMETER :: zghno380_ss = 1.05e-1
REAL, PARAMETER :: zrgas = 8.314               !J/mol/K
REAL, PARAMETER :: zmgas = 29.0e-3              !Molecular weight in Kg
REAL, PARAMETER :: zdmolec = 4.5e-10           !molec diameter in m
REAL, PARAMETER :: zavg = 6.02217e+23          ! Avogadro number [mol-1]
REAL, PARAMETER :: pfrac_ca=0.05               ! Fraction of dust that is Ca2+

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='AER_NO3_UKCADU'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise arrays
ptss(:,:,:,:) = 0.0
pthno3(:,:,:) = 0.0
ptno3(:,:,:,:) = 0.0
ptdu(:,:,:,:) = 0.0
DO jk = 1,nlevs
  DO jj = 1,rows
    DO ji = 1,row_length
      irh(ji,jj,jk)=1
      DO jtab=1,12
        IF (prhcl(ji,jj,jk)*100.0 > rrhtab(jtab)) THEN
          irh(ji,jj,jk)=jtab
        END IF
      END DO
    END DO
  END DO
END DO

DO jk = 1,nlevs
  DO jj = 1,rows
    DO ji = 1,row_length
      IF (jk <= ptroplev(ji,jj)) THEN
        iirh=irh(ji,jj,jk)

        zrhloc = prhcl(ji,jj,jk)
        IF (zrhloc < 0.0 ) zrhloc = 0.0
        IF (zrhloc > 0.98 ) zrhloc = 0.98

        ! needed for HNO3 update on dust and ss
        zdg0= 3.0/(8.0*zavg*prho(ji,jj,jk)*(zdmolec**2.0))
        zdghno3 = zdg0 * (zrgas*zmgas/2.0/rpi*pt(ji,jj,jk)*zratiohno3)**0.5

        !RH dependent HNO3 uptake coefficient based on Fairlie et al. 2010
        !and scale for SS
        zrh1580 = MAX(MIN(zrhloc,0.80),0.15)
        IF ( zrh1580 < 0.30 ) THEN
          zghno3_dd = zghno315_dd + (zghno330_dd-zghno315_dd)/                 &
                      0.15*(zrh1580-0.15)
          zghno3_ss = zghno315_ss + (zghno330_ss-zghno315_ss)/                 &
                      0.15*(zrh1580-0.15)
        ELSE IF (zrh1580 > 0.30 .AND. zrh1580 < 0.70) THEN
          zghno3_dd = zghno330_dd + (zghno370_dd-zghno330_dd)/                 &
                      0.40*(zrh1580-0.30)
          zghno3_ss = zghno330_ss + (zghno370_ss-zghno330_ss)/                 &
                      0.40*(zrh1580-0.30)
        ELSE IF (zrh1580 > 0.70 ) THEN
          zghno3_dd = zghno370_dd + (zghno380_dd-zghno370_dd)/                 &
                      0.10*(zrh1580-0.70)
          zghno3_ss = zghno370_ss + (zghno380_ss-zghno370_ss)/                 &
                      0.10*(zrh1580-0.70)
        END IF

        ! MFP parameter
        zmfp=3.0 * zdghno3 / SQRT(8.0*zrgas*pt(ji,jj,jk)/rpi/zmwhno3)

        ! Uptake of HNO3 on SEA-SALT
        sumdu = 0.0
        DO jdu = 1,ndumodes
          zkn_dd = 2.0 * zmfp / pddu(ji,jj,jk,jdu)
          zkhno3_dd(jdu) = (2.0*rpi*pddu(ji,jj,jk,jdu)*zdghno3) /              &
                           (1.0+(4.0*zkn_dd/zghno3_dd/3.0)*(1.0-0.47*zghno3_dd/&
                           (1.0+zkn_dd)))
          sumdu = sumdu + pmdu(ji,jj,jk,jdu)
        END DO

        ! Uptake of HNO3 on SEA-SALT
        DO jss = 1,nssmodes
          zkn_ss = 2.0 * zmfp / (pdss(ji,jj,jk,jss) * rssgrowth_rhtab(iirh))
          zkhno3_ss(jss) = (2.0*rpi*pdss(ji,jj,jk,jss)*rssgrowth_rhtab(iirh)*  &
                            zdghno3) /(1.0+(4.0*zkn_ss/zghno3_ss/3.0)*         &
                           (1.0-0.47*zghno3_ss/(1.0+zkn_ss)))
        END DO

        !! NOT APPLICABLE FOR UKCA DUST
        !! Calcium limitation for HNO3 condensation on dust
        !zcavmr = sumdu*pfrac_ca/zmwca*zmgas*0.5
        !         ! Moles of Ca available for Ca(NO3)2 formation
        !zcino3 = pnano3(ji,jj,jk)/zmwnano3*zmgas
        !IF (zcino3 > zcavmr)  zkhno3_dd(:) = 0.

        ! Check to ensure that dust and sea-salt tendencies do not
        ! result in negative aerosol mass
        DO jdu = 1,ndumodes
          ddtmp = MAX(0.0,(pndu(ji,jj,jk,jdu)*zmwca/zmwhno3*                   &
                          phno3(ji,jj,jk)*ptsphy))
          IF ((pmdu(ji,jj,jk,jdu)*pfrac_ca-ddtmp*zkhno3_dd(jdu)<0.0) .AND.     &
              (ddtmp>0.0))                                                     &
            zkhno3_dd(jdu) = pmdu(ji,jj,jk,jdu)*pfrac_ca/ddtmp
        END DO

        DO jss = 1,nssmodes
          sstmp = MAX(0.0,(pnss(ji,jj,jk,jss)*zmwnacl/zmwhno3*                 &
                          phno3(ji,jj,jk)*ptsphy))
          IF ((pmss(ji,jj,jk,jss)-sstmp*zkhno3_ss(jss)<0.0) .AND. (sstmp>0.0)) &
            zkhno3_ss(jss) = pmss(ji,jj,jk,jss)/sstmp
        END DO

        ! Divide NO3 tendencies into accum and coarse modes
        zdd1 = MAX(0.0,zkhno3_dd(1)*pndu(ji,jj,jk,1))*zmwcano3_2/zmwhno3
        zdd2 = MAX(0.0,zkhno3_dd(2)*pndu(ji,jj,jk,2))*zmwcano3_2/zmwhno3
        zss1 = MAX(0.0,zkhno3_ss(1)*pnss(ji,jj,jk,1))*zmwnano3/zmwhno3
        zss2 = MAX(0.0,zkhno3_ss(2)*pnss(ji,jj,jk,2))*zmwnano3/zmwhno3
        zaccno3 = zdd1 + zss1
        zcorno3 = zdd2 + zss2

        ! If the NO3 tendency is greater than 0 then store in output array
        IF ((zaccno3+zcorno3) > 0.0) THEN
          zfracm1 = zaccno3 / (zaccno3+zcorno3)
          zfracm2 = zcorno3 / (zaccno3+zcorno3)

          ! Update NO3 and HNO3 tendencies
          dtno3=(zaccno3+zcorno3)*MAX(0.0,phno3(ji,jj,jk))
          pthno3(ji,jj,jk) = -1.0*MAX(0.0,phno3(ji,jj,jk))*                    &
                             ((zdd1+zdd2)*zmwhno3*2.0/zmwcano3_2 +             &
                              (zss1+zss2)*zmwhno3/zmwnano3)

          ! Check for negative HNO3
          zhno3_test = phno3(ji,jj,jk) + pthno3(ji,jj,jk)*ptsphy
          mcf=1.0
          IF (zhno3_test < 0.0) THEN
            dtno3old=dtno3
            zthno3_tmp = pthno3(ji,jj,jk)
            pthno3(ji,jj,jk) = -1.0*phno3(ji,jj,jk)/ptsphy
            zthno3_tmp = pthno3(ji,jj,jk)-zthno3_tmp
            dtno3 = dtno3-zthno3_tmp*zmwnano3/zmwhno3
            mcf=dtno3/dtno3old
          END IF
          ptno3(ji,jj,jk,1) = dtno3 * zfracm1
          ptno3(ji,jj,jk,2) = dtno3 * zfracm2

          ! Update dust tendencies
          DO jdu = 1,ndumodes
            ptdu(ji,jj,jk,jdu)= -1.0 * zkhno3_dd(jdu) * pndu(ji,jj,jk,jdu) *   &
                                phno3(ji,jj,jk) * zmwcaco3 / zmwhno3 * mcf
          END DO

          ! Update sea-salt tendencies
          DO jss = 1,nssmodes
            ptss(ji,jj,jk,jss)= -1.0 * zkhno3_ss(jss) * pnss(ji,jj,jk,jss) *   &
                                phno3(ji,jj,jk) * zmwnacl / zmwhno3 * mcf
          END DO
        END IF  ! Rates are greater than zero
      END IF ! < tropopause
    END DO !ji
  END DO !jj
END DO !jk

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE aer_no3_ukcadu

END MODULE ukca_aer_no3_mod
