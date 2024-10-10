! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Main routine for reading in fast-jx initialisation files
!   The fastjx_rd_xxx_file subroutine reads the spectral data file that
!   contains cross sections of various species for given wavelength bins.
!   The cross sections used are derived from experimental measurements of
!   the photolysis rates. Please refer to this paper for further details:
!   doi:10.5194/gmd-6-161-2013
!
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!

MODULE fastjx_read_ascii_mod

USE fastjx_data, ONLY: a_, atau, atau0, daa, fl, jtaumx,                       &
                       naa, njval, nw1, nw2,                                   &
                       paa, q1d, qaa, qo2, qo3, qqq, qrayl, raa,               &
                       saa, solcyc_av, solcyc_quanta, solcyc_spec,             &
                       solcyc_ts, tqq,  w_, waa, wl, wx_, x_,                  &
                       n_solcyc_ts, n_solcyc_av

USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper,          &
                      PrStatus_Diag
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Default private
PRIVATE

! Subroutines available outside this module
PUBLIC :: fastjx_rd_mie_file, fastjx_rd_xxx_file, fastjx_rd_sol_file


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FASTJX_READ_ASCII_MOD'

CONTAINS

SUBROUTINE fastjx_rd_mie_file(nj1,namfil,aerosol_cloud_title)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nj1  ! Channel number for reading data file
CHARACTER(LEN=*), INTENT(IN) :: namfil ! Name of scattering data file
                                       ! (e.g., FJX_scat.dat)
! String containing cloud/aerosol scattering
CHARACTER(LEN=7), INTENT(IN OUT)  ::  aerosol_cloud_title(a_)

! String containing description of data set
CHARACTER(LEN=78) :: title0

! USEd variables from fastjx_data
! a_ is the max no. of aerosol/cloud mie data sets

CHARACTER (LEN=errormessagelength) :: cmessage
                                      ! Contains string for error handling
INTEGER :: errcode                    ! error code

INTEGER :: i, j, k                      ! Loop variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_RD_MIE_FILE'

! ***********************************
! End of Header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Open data file containing aerosol/cloud data
OPEN (nj1,FILE=namfil,STATUS='OLD',FORM='FORMATTED',ACTION='READ')

! Read number of data types and title
READ (UNIT=nj1,FMT='(I2,A78)') naa,title0
! If the number of data types exceeds maximum allowed exit with an error
IF (naa > a_) THEN
  cmessage = 'Too many scattering data sets'
  errcode = 100
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF (printstatus >= prstatus_oper) THEN
  ! Output file title
  WRITE(umMessage,'(A)') title0
  CALL umPrint(umMessage,src='fastjx_read_ascii')
END IF

! Read Cloud layering variables:
! Maximum number of cloud sub layers (jtaumx)
! Cloud sub layers factor (atau)
! minimum dtau (atau0)
READ (UNIT=nj1,FMT='(5X,I5,2F10.5)') jtaumx,atau,atau0
IF (printstatus >= prstatus_oper) THEN
  WRITE(umMessage,'(A,2F9.5,I5)')                                              &
      ' atau/atau0/jmx',atau,atau0,jtaumx
  CALL umPrint(umMessage,src='fastjx_read_ascii')
END IF

! Read blank line
READ (UNIT=nj1,FMT=*)

! Loop over aerosol types
DO j = 1,naa
  ! Read title (aerosol_cloud_title), effective radius (raa) and density (daa)
  ! of scattering types
  READ (UNIT=nj1,FMT='(3X,A7,45X,F5.3,15X,F5.3)')  aerosol_cloud_title(j),     &
                                                   raa(j),daa(j)
  ! Loop over 5 wavelength bins
  DO k = 1,5
    ! read wavelength (waa), q (qaa), scattering albedo (saa) and phases (paa)
    ! of scattering types
    READ (UNIT=nj1,FMT='(F4.0,F7.4,F7.4,7F6.3,1X,F7.3,F8.4)')                  &
      waa(k,j),qaa(k,j),saa(k,j),(paa(i,k,j),i=2,8)
    ! set first phase to 0
    paa(1,k,j) = 1.0e0
  END DO ! wavelengths
END DO ! aerosols

! Close file
CLOSE(nj1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_rd_mie_file

SUBROUTINE fastjx_rd_xxx_file(nj1,namfil,titlej)
IMPLICIT NONE

INTEGER, INTENT(IN)       :: nj1 ! Channel number for reading data file
CHARACTER(LEN=*), INTENT(IN)  :: namfil    ! Name of spectral data file
                                           ! (e.g. JX_spec.dat)
! String containing species being photolysed
CHARACTER(LEN=7), INTENT(IN OUT)  ::  titlej(x_)

! Dummy strings containing duplicate species strings
CHARACTER(LEN=7)  ::  titlej2,titlej3

! String containing description of data set
CHARACTER(LEN=78) :: title0

! USEd variables from fastjx_data
! w_ is the number of wavelength bins (can choose smaller for speed)
!      - Set using RUN_UKCA namelist
! wx_ is the max number of wavelength bins
! x_ is the max no. of cross sections to be read

CHARACTER (LEN=errormessagelength)        :: cmessage
                                           ! String for error handling
INTEGER                   :: errcode       ! errror code

INTEGER ::  i, j, jj, k, iw                ! Loop variables

INTEGER ::  n_xsections_read               ! No of cross sections read in (nqqq)
INTEGER ::  n_wl_bins                      ! No of wavelength bins (nwww)
INTEGER ::  n_xsections                    ! number of x-sections (nqrd)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_RD_XXX_FILE'

  ! *****************************

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  !----------spectral data----set for new format data J-ver8.3------------------
!         note that NJVAL = # J-values, but
!              n_xsections_read (>NJVAL) = # Xsects read in
!         for 2005a data, NJVAL = 62 (including a spare XXXX) and
!              n_xsections_read = 64 so that 4 wavelength datasets read in for
!              acetone
!         note n_xsections_read is not used outside this subroutine!
! >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects

! Open file containing cross sections
OPEN (nj1,FILE=namfil,STATUS='old',FORM='formatted',ACTION='READ')

! Read title of file
READ (UNIT=nj1,FMT='(A)') title0

IF (printstatus >= prstatus_oper) THEN
  ! Output file title
  WRITE(umMessage,'(A)') title0
  CALL umPrint(umMessage,src='fastjx_read_ascii')
END IF

! Read number of photolysed species, number of x-sections & number
! of wavelength bins
READ (UNIT=nj1,FMT='(10X,5I5)') njval, n_xsections, n_wl_bins

! set maximum and minimum wavelngth bins
nw1 = 1
nw2 = n_wl_bins

! Check that number of photolysed species and number of cross sections
! doesn't exceed maximum allowed. If either do then exit
IF (njval > x_ .OR. n_xsections > x_) THEN
  cmessage = 'Number of Cross Sections exceeds Maximum Allowed'
  errcode = 100
  CALL ereport(RoutineName,errcode,cmessage)
END IF

!----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
! Read the effective wavelengths
READ (UNIT=nj1,FMT='(10X,    6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      (wl(iw),iw=1,n_wl_bins)
! Read the top of atmosphere solar flux
READ (UNIT=nj1,FMT='(10X,    6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      (fl(iw),iw=1,n_wl_bins)
! Read the Rayleigh parameters (effective cross-section) (cm2)
READ (UNIT=nj1,FMT='(10X,    6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      (qrayl(iw),iw=1,n_wl_bins)

!---READ O2 X-sects, O3 X-sects, O3=>O(1D) quant yields
!     (each at 3 temperatures(tqq))
READ (UNIT=nj1,FMT='(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      titlej(1),tqq(1,1), (qo2(iw,1),iw=1,n_wl_bins)
READ (UNIT=nj1,FMT='(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      titlej2,  tqq(2,1), (qo2(iw,2),iw=1,n_wl_bins)
READ (UNIT=nj1,FMT='(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      titlej3,  tqq(3,1), (qo2(iw,3),iw=1,n_wl_bins)

READ (UNIT=nj1,FMT='(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      titlej(2),tqq(1,2), (qo3(iw,1),iw=1,n_wl_bins)
READ (UNIT=nj1,FMT='(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      titlej2,  tqq(2,2), (qo3(iw,2),iw=1,n_wl_bins)
READ (UNIT=nj1,FMT='(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      titlej3,  tqq(3,2), (qo3(iw,3),iw=1,n_wl_bins)

READ (UNIT=nj1,FMT='(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      titlej(3),tqq(1,3), (q1d(iw,1),iw=1,n_wl_bins)
READ (UNIT=nj1,FMT='(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      titlej2,  tqq(2,3), (q1d(iw,2),iw=1,n_wl_bins)
READ (UNIT=nj1,FMT='(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))')               &
      titlej3,  tqq(3,3), (q1d(iw,3),iw=1,n_wl_bins)

! WRITE information to stdout
IF (printstatus >= prstatus_oper) THEN
  DO j = 1,3
    WRITE(umMessage,'(I6,A7,3E10.3)') j,titlej(j),(tqq(i,j),i=1,3)
    CALL umPrint(umMessage,src='fastjx_read_ascii')
  END DO
END IF

!---READ remaining species:  X-sections (qqq) at 2 temperatures (tqq)
jj = 4
DO j = 4,n_xsections
  READ (UNIT=nj1,FMT='(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))')             &
      titlej(jj),tqq(1,jj),(qqq(iw,1,jj),iw=1,n_wl_bins)
  READ (UNIT=nj1,FMT='(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))')             &
      titlej2,  tqq(2,jj),(qqq(iw,2,jj),iw=1,n_wl_bins)

  !---include stratospheric J's (this also includes Cl and Br compounds!)
  IF (w_ == 18 .OR. titlej2(7:7) /= 'x') THEN
    IF (printstatus >= prstatus_oper) THEN
      WRITE(umMessage,'(I6,A7,2E10.3)') jj,titlej(jj), (tqq(i,jj),i=1,2)
      CALL umPrint(umMessage,src='fastjx_read_ascii')
    END IF
    jj = jj+1
  END IF

END DO
n_xsections_read = jj-1
njval = njval + (n_xsections_read - n_xsections)

! Close file
CLOSE(nj1)

!---truncate number of wavelengths to DO troposphere-only
IF (w_ /= wx_) THEN

  !---TROP-ONLY
  IF (w_ == 12) THEN
    IF (printstatus >= prstatus_oper) THEN
      CALL umPrint(' >>>TROP-ONLY reduce wavelengths to 12,'//                 &
          ' drop strat X-sects',                                               &
          src='fastjx_read_ascii')
    END IF
    nw2 = 12

    ! Remove first four wavelength bins from  total
    DO iw = 1,4
      wl(iw) = wl(iw+4)
      fl(iw) = fl(iw+4)
      qrayl(iw) = qrayl(iw+4)

      DO k = 1,3
        qo2(iw,k) = qo2(iw+4,k)
        qo3(iw,k) = qo3(iw+4,k)
        q1d(iw,k) = q1d(iw+4,k)
      END DO

      DO j = 4,n_xsections_read
        qqq(iw,1,j) = qqq(iw+4,1,j)
        qqq(iw,2,j) = qqq(iw+4,2,j)
      END DO
    END DO

    ! Remove 9/10 wavelength bins from total
    DO iw = 5,12
      wl(iw) = wl(iw+6)
      fl(iw) = fl(iw+6)
      qrayl(iw) = qrayl(iw+6)

      DO k = 1,3
        qo2(iw,k) = qo2(iw+6,k)
        qo3(iw,k) = qo3(iw+6,k)
        q1d(iw,k) = q1d(iw+6,k)
      END DO
      DO j = 4,n_xsections_read
        qqq(iw,1,j) = qqq(iw+6,1,j)
        qqq(iw,2,j) = qqq(iw+6,2,j)
      END DO
    END DO

    !---TROP-QUICK  (must scale solar flux for W=5)
  ELSE IF (w_ == 8) THEN
    IF (printstatus >= prstatus_oper) THEN
      CALL umPrint(' >>>TROP-QUICK reduce wavelengths to 8, '//                &
          'drop strat X-sects',                                                &
          src='fastjx_read_ascii')
    END IF
    nw2 = 8

    DO iw = 1,1
      wl(iw) = wl(iw+4)
      fl(iw) = fl(iw+4)*2.0e0
      qrayl(iw) = qrayl(iw+4)

      DO k = 1,3
        qo2(iw,k) = qo2(iw+4,k)
        qo3(iw,k) = qo3(iw+4,k)
        q1d(iw,k) = q1d(iw+4,k)
      END DO

      DO j = 4,n_xsections_read
        qqq(iw,1,j) = qqq(iw+4,1,j)
        qqq(iw,2,j) = qqq(iw+4,2,j)
      END DO
    END DO

    DO iw = 2,8
      wl(iw) = wl(iw+10)
      fl(iw) = fl(iw+10)
      qrayl(iw) = qrayl(iw+10)

      DO k = 1,3
        qo2(iw,k) = qo2(iw+10,k)
        qo3(iw,k) = qo3(iw+10,k)
        q1d(iw,k) = q1d(iw+10,k)
      END DO

      DO j = 4,n_xsections_read
        qqq(iw,1,j) = qqq(iw+10,1,j)
        qqq(iw,2,j) = qqq(iw+10,2,j)
      END DO
    END DO

  ELSE
    cmessage = 'Incorrect Number of Wavelength Bins, must be 8, 12 or 18'
    errcode = 100
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_rd_xxx_file

SUBROUTINE fastjx_rd_sol_file(nj1,namfil)

USE ukca_parpho_mod, ONLY: jpwav
IMPLICIT NONE

INTEGER, INTENT(IN)       :: nj1 ! Channel number for reading data file
CHARACTER(LEN=*), INTENT(IN)  :: namfil    ! Name of solar cycle data file

CHARACTER (LEN=errormessagelength) :: cmessage
                                      ! Contains string for error handling
INTEGER :: errcode                    ! error code

INTEGER ::  i, j                           ! Loop variables
INTEGER ::  n                              ! Array length

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_RD_SOL_FILE'


OPEN (nj1,FILE=namfil,STATUS='old',FORM='formatted',ACTION='READ')

! Read in the spectral component of the solar cycle.
READ (UNIT=nj1,FMT='(50X,I2)') n
! If the number of data points does not equal correct size, raise error
IF (n /= wx_) THEN
  cmessage = 'Incorrect number of wavelengths'
  errcode = n
  CALL ereport(RoutineName,errcode,cmessage)
END IF

DO i = 1,CEILING(REAL(n)/6.0)
  READ (UNIT=nj1,FMT='(6E10.3)')                                               &
          (solcyc_spec(j+(i-1)*6),j=1,6)
END DO

! Read in the quanta component of the solar cycle.
READ (UNIT=nj1,FMT='(50X,I3)') n
! If the number of data points does not equal correct size, raise error
IF (n /= jpwav) THEN
  cmessage = 'Incorrect number of quanta'
  errcode = n
  CALL ereport(RoutineName,errcode,cmessage)
END IF

DO i = 1,CEILING(REAL(n)/7.0)
  READ (UNIT=nj1,FMT='(7E10.3)')                                               &
        (solcyc_quanta(j+(i-1)*7),j=1,7)
END DO

! Read in the obs. time series of the solar cycle.
READ (UNIT=nj1,FMT='(50X,I4)') n
n_solcyc_ts = n

! Now allocate array to be the correct size
IF (.NOT. ALLOCATED(solcyc_ts)) THEN
  ALLOCATE(solcyc_ts(n_solcyc_ts))
ELSE
   ! already allocated, check size
  IF (n_solcyc_ts /= SIZE(solcyc_ts)) THEN
    cmessage = 'Incorrect number of times for full solar cycle'
    errcode = n
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF

DO i = 1,CEILING(REAL(n)/6.0)
  READ (UNIT=nj1,FMT='(6F9.5)')                                                &
        (solcyc_ts(j+(i-1)*6),j=1,6)
END DO

! Read in the average solar cycle.
READ (UNIT=nj1,FMT='(50X,I3)') n
! If the number of data points does not equal correct size, raise error
IF (n /= n_solcyc_av) THEN
  cmessage = 'Incorrect number of times for average cycle'
  errcode = n
  CALL ereport(RoutineName,errcode,cmessage)
END IF

DO i = 1,CEILING(REAL(n)/8.0)
  READ (UNIT=nj1,FMT='(8F9.5)')                                                &
          (solcyc_av(j+(i-1)*8),j=1,8)
END DO

CLOSE(nj1)

END SUBROUTINE fastjx_rd_sol_file

END MODULE fastjx_read_ascii_mod