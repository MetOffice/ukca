! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description: Calculates number concentration of aerosol particles
!               which become "activated" into cloud droplets from MODE
!               aerosol results for mass mixing ratio, number and dry radius
!               together with updraught velocities.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge,
!  University of Leeds, University of Oxford and The Met Office.
!  See www.ukca.ac.uk
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
! ----------------------------------------------------------------------

MODULE ukca_abdulrazzak_ghan_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_ABDULRAZZAK_GHAN_MOD'

CONTAINS

SUBROUTINE ukca_abdulrazzak_ghan(kbdim,   klev,                                &
                                 glomap_variables_local,                       &
                                 pesw,                                         &
                                 pn,      pxtm1,                               &
                                 ptm1,    papm1,                               &
                                 pqm1,    prdry,                               &
                                 nwbins,  pwarr,                               &
                                 pwpdf,   pwbin,                               &
                                 l_fix_ukca_hygroscopicities_local,            &
                                 psmax,   pwchar,                              &
                                 pcdncactm, pcdncact )


! *aero_activ* calculates the number of activated aerosol
!              particles from the aerosol size-distribution,
!              composition and ambient supersaturation
!
! Author:
! -------
! Philip Stier, MPI-MET                  2002/2003
! Rosalind West, AOPP, Oxford            2008
!
! Method:
! -------
! The calculation of the activation can be reduced to 3 tasks:
!
! I)   Calculate the maximum supersaturation
! II)  Calculate the corresponding radius of activation
!      for each mode
! III) Calculate the number of particles that are larger
!      then the radius of activation for each mode.
!
! III) Calculation of the number of activated particles:
!      See the routine aero_activ_tail below.
!
! References:
! -----------
! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998.
! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000.
! Pruppbacher and Klett, Kluewer Ac. Pub., 1997.
! West et al., ACP, 14, 2014. doi:10.5194/acp-14-6369-2014

USE ukca_mode_setup,       ONLY: nmodes, glomap_variables_type

USE ukca_activ_mod,        ONLY: activclosest
USE ukca_um_legacy_mod,    ONLY: gg => g,  & ! Acceleration due to gravity
                                 cp,       & ! Specific heat at const p J/(kg.K)
                                 umErf

USE ukca_config_constants_mod, ONLY: rmol,      & ! gas constant
                                     rho_water, & ! density H2O kg/m^3
                                     lc,        & ! latent heat of condensation
                                                  ! J/kg
                                     l_ukca_constants_available

USE ukca_constants,        ONLY: zerodegc, &  ! 0 degrees C in K
                                 pi,       &  ! Pi
                                 mmw,      &  ! H2O molecular weight kg/mol
                                 m_air,    &  ! Dry air  ------ " ---------
                                 zsten,    &  ! surface tension of H2O [J m-2]
                                 zosm         ! Osmotic coefficient
USE umPrintMgr, ONLY: umMessage, umPrint, PrintStatus, PrStatus_Diag
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

IMPLICIT NONE

!--- Arguments:

INTEGER, INTENT(IN) :: kbdim                    !
INTEGER, INTENT(IN) :: klev                     !
TYPE(glomap_variables_type), TARGET, INTENT(IN) :: glomap_variables_local
INTEGER, INTENT(IN) :: nwbins                   !

REAL, INTENT(IN)  :: ptm1(kbdim,klev)         ! temperature
REAL, INTENT(IN)  :: papm1(kbdim,klev)        ! pressure
REAL, INTENT(IN)  :: pqm1(kbdim,klev)         ! specific humidity
REAL, INTENT(IN)  :: pesw(kbdim,klev)         ! saturation water vapour pressure
REAL, INTENT(IN)  :: prdry(kbdim,klev,nmodes) ! dry radius for each mode
REAL, INTENT(IN)  :: pn (kbdim,klev,nmodes)   ! aerosol number concentration
!                                                     ! for each mode [m-3]
! aerosol tracers mass
REAL, INTENT(IN)  :: pxtm1(kbdim,klev,nmodes,glomap_variables_local%ncp)
!                                                        ! mixing ratio
! pdf of updraught velocities:
REAL,INTENT(IN)     :: pwarr(kbdim,klev,nwbins)  ! lin array of vert vel [m s-1]
REAL,INTENT(IN)     :: pwpdf(kbdim,klev,nwbins)  ! lin array of pdf of w [m s-1]
REAL,INTENT(IN)     :: pwbin(kbdim,klev,nwbins)  ! w bin width [m s-1]

! local copy of temp logical from glomap_config for hygroscopicity fix
LOGICAL,INTENT(IN) :: l_fix_ukca_hygroscopicities_local

! max supersaturation [fraction]
REAL, INTENT(OUT)   :: psmax(kbdim,klev)

! calculated characteristic updraught vel [m s-1]
REAL, INTENT(OUT)   :: pwchar(kbdim,klev)

! expected number of activated particles in each mode
REAL, INTENT(OUT)   :: pcdncactm(kbdim,klev,nmodes)

! expected number of activated particles over all modes
REAL, INTENT(OUT)   :: pcdncact(kbdim,klev)

!--- Local variables:

! Caution - pointers to TYPE glomap_variables_local%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
LOGICAL, POINTER :: component(:,:)
REAL,    POINTER :: mm(:)
LOGICAL, POINTER :: mode(:)
INTEGER, POINTER :: ncp
REAL,    POINTER :: no_ions(:)
REAL,    POINTER :: rhocomp(:)
REAL,    POINTER :: sigmag(:)
INTEGER, POINTER :: modesol(:)
INTEGER, POINTER :: topmode

REAL :: zsigmaln(nmodes)   ! ln(geometric std dev)

INTEGER :: jmod            ! loop counters
INTEGER :: jcp
INTEGER :: jl
INTEGER :: jk
INTEGER :: jw

REAL :: zmassfrac
REAL :: zalpha
REAL :: zeps
REAL :: zgamma
REAL :: zf
REAL :: zg
REAL :: zxi
REAL :: zeta
REAL :: zsum
REAL :: zgrowth
REAL :: zdif
REAL :: zk
REAL :: zka
REAL :: zerf_ratio
REAL :: zwpwdw                      !
REAL :: zpwdw                       !

REAL :: zsmax(kbdim,klev)           ! maximum supersaturation
REAL :: zsumtop(kbdim,klev)         !
REAL :: zsumbot(kbdim,klev)         !
REAL :: zw                          ! total vertical velocity[m s-1]

REAL :: zrc                         ! critical radius of a dry aerosol particle
!                                   ! that becomes activated at the ambient
!                                   ! radius of activation
REAL :: zmasssum(kbdim,klev,nmodes) !
REAL :: za(kbdim,klev,nmodes)       ! curvature parameter A of the Koehler eqn
REAL :: zb(kbdim,klev,nmodes)       ! hygroscopicity parameter B of Koehler eqn
REAL :: zsm(kbdim,klev,nmodes)      ! critical supersaturation for activating
!                                   ! particles with the mode No. median radius

REAL :: zcdnc(kbdim,klev,nwbins)    ! CDNC calculated at each increment of w
REAL :: zcdnc_tmp(nwbins)           ! Temporary array for CDNC
REAL :: zcdncm                      ! CDNC calculated at each increment of w
!                                         ! by mode
REAL :: psmax2(kbdim,klev)          ! max supersaturation [fraction]
REAL :: zndtopm(kbdim,klev,nmodes,nwbins)  ! top line integral for calc
!                                               ! pcdncactm
REAL :: zndbotm(kbdim,klev,nmodes,nwbins)  ! bottom line integral for calc
!                                               ! pcdncactm
REAL :: zndtop(kbdim,klev,nwbins)   ! top line integral for calc pcdncact
REAL :: zndbot(kbdim,klev,nwbins)   ! bottom line integral for calc pcdncact

! local variables used for finding wchar
INTEGER :: zclsloc
REAL :: zcls

REAL, PARAMETER :: p0=101325.0     ! in Pa

REAL :: cthomi
REAL :: lc_sq

REAL :: no_ions_div_mm( glomap_variables_local%ncp )

INTEGER                           :: errcode  ! error code
CHARACTER(LEN=errormessagelength) :: cmessage ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_ABDULRAZZAK_GHAN'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allow for the possibility that this routine is called from outside UKCA,
! by-passing the UKCA API. This is only allowed if UKCA constants are properly
! set up.
IF (.NOT. l_ukca_constants_available) THEN
  cmessage = 'Configurable UKCA constants have not been set up'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Caution - pointers to TYPE glomap_variables_local%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables_local%component
mm          => glomap_variables_local%mm
mode        => glomap_variables_local%mode
ncp         => glomap_variables_local%ncp
no_ions     => glomap_variables_local%no_ions
rhocomp     => glomap_variables_local%rhocomp
sigmag      => glomap_variables_local%sigmag
modesol     => glomap_variables_local%modesol
topmode     => glomap_variables_local%topmode

cthomi  = zerodegc - 35.0
lc_sq=lc**2

!--- 0) Initializations:

zmasssum(:,:,:) = 0.0
za(:,:,:)       = 0.0
zb(:,:,:)       = 0.0
pcdncact(:,:)   = 0.0
pcdncactm(:,:,:)= 0.0
zcdnc(:,:,:)    = 0.0
zcdncm          = 0.0

zeps=EPSILON(1.0)

!--- calculate ln(sigmag)
DO jmod=1, topmode
  zsigmaln(jmod) = LOG(sigmag(jmod))
END DO

!--- 1) Calculate properties for each aerosol mode:

     !--- 1.1) Calculate the auxiliary parameters A & B of the Koehler equation:

DO jmod=1, topmode
  IF (mode(jmod)) THEN

    !--- 1.1.0 Initializations:

    zsumtop(:,:)    = 0.0
    zsumbot(:,:)    = 0.0

    !--- 1.1.1) Calculate the mean hygroscopicity parameter B:

    !--- 1) Calculate weighted properties:
    !       (Abdul-Razzak & Ghan, 2000; see also Eqn. A2 of Ghan (2011;
    !       doi:10.1029/2011MS000074))
    ! N.B. Molar masses (mm) of components in UKCA already given in kg mol-1

    IF (l_fix_ukca_hygroscopicities_local) THEN
      !Here the hygroscopicity fix is switched on which applies a bug-fix and
      !updates to hygroscopicity values. The 'soluble mass fraction' (epsilon
      !in ARG2000 and Ghan, 2011) is not required since this stems from the case
      !of a two-component insoluble core surrounded by soluble material (see
      !the Pruppacher and Klett Microphysics of Clouds and Precipitation
      !textbook), whereas we have several (internally mixed) components here with
      !varying solubilities that are taken into account through the no_ions array.

      DO jcp=1,ncp
        no_ions_div_mm(jcp)=no_ions(jcp)/mm(jcp)
        IF (component(jmod,jcp)) THEN
          DO jk=1, klev
            DO jl=1, kbdim
              zsumtop(jl,jk)=zsumtop(jl,jk)+                                   &
                             (pxtm1(jl,jk,jmod,jcp)*no_ions_div_mm(jcp))

              zsumbot(jl,jk)=zsumbot(jl,jk)+                                   &
                             (pxtm1(jl,jk,jmod,jcp)/rhocomp(jcp))
            END DO        ! jl
          END DO           ! jk
        END IF              ! component
      END DO                 ! loop over cpts

    ELSE
      !Hygroscopicity fix is off.

      !--- Sum properties over all composition components for the calcuation
      ! of the 'soluble mass fraction' (zmassfrac; epsilon in ARG2000
      ! and Ghan, 2011). N.B., the use of zmassfrac represents a bug that is
      ! fixed by the l_fix_ukca_hygroscopicities temporary logical.
      DO jcp=1,ncp
        IF (component(jmod,jcp)) THEN
          DO jk=1, klev
            DO jl=1, kbdim
              zmasssum(jl,jk,jmod) = zmasssum(jl,jk,jmod) +                    &
                                      pxtm1(jl,jk,jmod,jcp)
            END DO
          END DO
        END IF
      END DO

      DO jcp=1,ncp
        IF (component(jmod,jcp) .AND.                                          &
            NINT(no_ions(jcp)) > 0) THEN
          !Sum properties over only the soluble compounds
          !(nion=0 for insoluble compounds)
          DO jk=1, klev
            DO jl=1, kbdim
              IF (zmasssum(jl,jk,jmod) > zeps) THEN

                zmassfrac=pxtm1(jl,jk,jmod,jcp)/                               &
                           zmasssum(jl,jk,jmod)

                zsumtop(jl,jk)=zsumtop(jl,jk)+                                 &
                               (pxtm1(jl,jk,jmod,jcp)*                         &
                               no_ions(jcp)*zosm*                              &
                               zmassfrac/mm(jcp))
                zsumbot(jl,jk)=zsumbot(jl,jk)+                                 &
                               (pxtm1(jl,jk,jmod,jcp)/                         &
                               rhocomp(jcp))
              END IF
            END DO        ! jl
          END DO           ! jk
        END IF              ! no_ions>0 and component
      END DO                 ! loop over cpts

    END IF
    ! END IF for hygroscopicity fix

    DO jk=1, klev
      DO jl=1, kbdim
        IF (zsumbot(jl,jk)>zeps) THEN
           !--- 1.1.1) Hygroscopicity parameter B (Eq. 4):

          zb(jl,jk,jmod)=(mmw*zsumtop(jl,jk))/                                 &
                         (rho_water*zsumbot(jl,jk))

          !--- 1.1.2) Calculate the curvature parameter A:

          za(jl,jk,jmod)=(2.0*zsten*mmw)/                                      &
                         (rho_water*rmol*ptm1(jl,jk))

        END IF         !zsumbot > 0
      END DO            !jl
    END DO               !jk
  END IF                  !mode
END DO                    !jmod=1, topmode

!2) Calculate maximum supersaturation at each increment of vertical velocity pdf

!$OMP PARALLEL DEFAULT(NONE)                                                   &
!$OMP PRIVATE(jk, jl, jmod, jw, psmax2,                                        &
!$OMP         zalpha, zcdncm, zdif, zerf_ratio, zeta,                          &
!$OMP         zf, zg, zgamma, zgrowth, zk, zka, zpwdw,                         &
!$OMP         zrc, zsm, zsmax, zsum, zw, zwpwdw, zxi)                          &
!$OMP SHARED(cp, cthomi, gg, kbdim, klev, lc,lc_sq, mode, nwbins, topmode,     &
!$OMP        modesol, papm1, pesw, pn, pqm1, prdry, Printstatus,               &
!$OMP        psmax, ptm1, pwarr, pwbin, pwpdf, rho_water, rmol,                &
!$OMP        za, zb, zcdnc, zeps,                                              &
!$OMP        zndbot, zndbotm, zndtop, zndtopm,                                 &
!$OMP        zsigmaln, l_fix_ukca_hygroscopicities_local)

!$OMP DO SCHEDULE(DYNAMIC)
DO jw=1, nwbins

  !    --- 2.2) Abbdul-Razzak and Ghan (2000):
  !         (Equations numbers from this paper unless otherwise quoted)

  zsm(:,:,:)      = 0.0

  DO jk=1, klev
    DO jl=1, kbdim
       ! get updraught velocity from pdf array
      zw=pwarr(jl,jk,jw)

       !--- Water vapour pressure:

      IF ( (zw > zeps) .AND. (pqm1(jl,jk) > zeps) .AND.                        &
           (ptm1(jl,jk) > cthomi) ) THEN

         !--- Abdul-Razzak et al. (1998) (Eq. 11):

        zalpha=(gg*mmw*lc)/(cp*rmol*ptm1(jl,jk)**2) -                          &
               (gg*m_air)/(rmol*ptm1(jl,jk))

        zgamma=(rmol*ptm1(jl,jk))/(pesw(jl,jk)*mmw) +                          &
               (mmw*lc_sq)/(cp* papm1(jl,jk)*m_air*ptm1(jl,jk))

        !--- Diffusivity of water vapour in air (P&K, 13.3) [m2 s-1]:
        !---  RW Code from Steve Ghan

        zdif=0.211*(p0/papm1(jl,jk))*                                          &
             ((ptm1(jl,jk)/zerodegc)**1.94)*1.0e-4 !in m^2/s

        !--- Thermal conductivity zk (P&K, 13.18) [cal cm-1 s-1 K-1]:

        ! Mole fraction of water:
        zka=(5.69+0.017*(ptm1(jl,jk)-zerodegc))*1.0e-5

        ! Moist air, convert to [J m-1 s-1 K-1]:

        zk=zka*4.186*1.0e2

        ! --- Abdul-Razzak et al. (1998) (Eq. 16):

        zgrowth=1.0/                                                           &
                  ( (rho_water*rmol*ptm1(jl,jk))/                              &
                    (pesw(jl,jk)*zdif*mmw) +                                   &
                    (lc*rho_water)/(zk*ptm1(jl,jk)) *                          &
                  ((lc*mmw)/(ptm1(jl,jk)*rmol) -1.0) )

        !--- Summation for equation (6):

        zsum=0.0

        DO jmod=1, topmode
          IF (mode(jmod) .AND. (modesol(jmod) == 1 .OR. .NOT.                  &
             l_fix_ukca_hygroscopicities_local)) THEN
            IF (pn(jl,jk,jmod)   >zeps     .AND.                               &
                 prdry(jl,jk,jmod)>1.0e-9   .AND.                              &
                 za(jl,jk,jmod)   >zeps    .AND.                               &
                 zb(jl,jk,jmod)   >zeps        ) THEN

               ! (7):
              zf=0.5*EXP(2.5*zsigmaln(jmod)**2)

              ! (8):
              zg=1.0+0.25*zsigmaln(jmod)

              ! (10):
              zxi=2.0*za(jl,jk,jmod)/3.0 *                                     &
                   SQRT(zalpha*zw/zgrowth)

              ! (11):
              zeta=((zalpha*zw/zgrowth)**1.5) /                                &
                   (2.0*pi*rho_water*zgamma*pn(jl,jk,jmod))

              ! (9):
              zsm(jl,jk,jmod)=2.0/SQRT(zb(jl,jk,jmod)) *                       &
                   (za(jl,jk,jmod)/                                            &
                   (3.0*prdry(jl,jk,jmod)))**1.5

              ! (6):
              zsum=zsum + (1.0/zsm(jl,jk,jmod)**2 *                            &
                   ( zf*(zxi/zeta)**1.5 +                                      &
                   zg*(zsm(jl,jk,jmod)**2/                                     &
                   (zeta+3.0*zxi))**0.75) )
            END IF
          END IF !mode
        END DO ! jmod

        IF (zsum > zeps) THEN
          zsmax(jl,jk)=1.0/SQRT(zsum)
        ELSE
          zsmax(jl,jk)=0.0
        END IF

      ELSE
        zsmax(jl,jk)=0.0
      END IF
      psmax2(jl,jk)=zsmax(jl,jk)
    END DO ! jl
  END DO ! jk

  !--- 3) Calculate activation:

  !   ---3.1) Calculate the critical radius (12):

  zndtopm(:,:,:,jw)= 0.0
  zndbotm(:,:,:,jw)= 0.0
  zndtop(:,:,jw)   = 0.0
  zndbot(:,:,jw)   = 0.0

  DO jk=1, klev
    DO jl=1, kbdim
      ! Calculate p(w)*dw
      zpwdw=pwpdf(jl,jk,jw)*pwbin(jl,jk,jw)
      ! Calculate w*p(w)*dw
      zwpwdw= pwarr(jl,jk,jw)*zpwdw
      IF (zwpwdw < zeps .AND. Printstatus == prstatus_diag) THEN
        WRITE(umMessage,'(A,I0,A,I0,A,F0.5)')                                  &
                          'RW: zwpwdw(',jl,',',jk,')= ',  zwpwdw
        CALL umPrint(umMessage,src='ukca_abdulrazzak_ghan')
        WRITE(umMessage,'(A,I0,A,I0,A,I0,A,F0.5)')                             &
                          'RW: pwpdf(',jl,',',jk,',',jw,')=  ', pwpdf(jl,jk,jw)
        CALL umPrint(umMessage,src='ukca_abdulrazzak_ghan')
        WRITE(umMessage,'(A,I0,A,I0,A,I0,A,F0.5)')                             &
                          'RW: pwbin(',jl,',',jk,',',jw,')=  ', pwbin(jl,jk,jw)
        CALL umPrint(umMessage,src='ukca_abdulrazzak_ghan')
        WRITE(umMessage,'(A,I0,A,I0,A,I0,A,F0.5)')                             &
                          'RW: pwarr(',jl,',',jk,',',jw,')=  ', pwarr(jl,jk,jw)
        CALL umPrint(umMessage,src='ukca_abdulrazzak_ghan')
      END IF
      DO jmod=1, topmode
        IF (mode(jmod) .AND. (modesol(jmod) == 1 .OR. .NOT.                    &
            l_fix_ukca_hygroscopicities_local)) THEN
          IF (psmax2(jl,jk)          >zeps  .AND.                              &
               zsm(jl,jk,jmod)       >zeps  .AND.                              &
               pn(jl,jk,jmod)        >zeps  .AND.                              &
               prdry(jl,jk,jmod)     >1.0e-9  ) THEN

            zrc = prdry(jl,jk,jmod)*(zsm(jl,jk,jmod)/                          &
                                     psmax2(jl,jk))**(2.0/3.0)
                 !--- 3.2) Calculate the total number of activated droplets
                 !         larger than the critical radii for each mode
            zerf_ratio=LOG(zrc/prdry(jl,jk,jmod))/SQRT(2.0)/zsigmaln(jmod)
            zcdnc(jl,jk,jw)=zcdnc(jl,jk,jw) + 0.5*pn(jl,jk,jmod)*              &
                             (1-umErf(zerf_ratio))
                 ! separate out by mode
            zcdncm = 0.5*pn(jl,jk,jmod)*(1-umErf(zerf_ratio))
                 !! Calculate the expected value of CDNC
                 !! for each mode over the whole pdf of w
                 !!pcdncactm(jl,jk,jmod)=pcdncactm(jl,jk,jmod)+    &
                 !!     zcdncm*pwpdf(jl,jk,jw)*pwbin(jl,jk,jw)
                 ! Calculate top integral for
                 ! expected value of CDNC for each mode
                 ! over the whole pdf of w, weighted by w:
            zndtopm(jl,jk,jmod,jw) = zcdncm*zpwdw
          END IF
              ! Calculate bottom integral for
              ! expected value of CDNC for each mode
              ! over the whole pdf of w, weighted by w:
          zndbotm(jl,jk,jmod,jw) = zpwdw
        END IF !mode
      END DO ! jmod

       !! Calculate the expected value of CDNC
       !! over the pdf of w
       !!pcdncact(jl,jk)=pcdncact(jl,jk) + &
       !!     zcdnc(jl,jk,jw)*pwpdf(jl,jk,jw)*pwbin(jl,jk,jw)
       ! Calculate top and bottom integrals for
       ! expected value of CDNC
       ! over the pdf of w, weighted by w:
      zndtop(jl,jk,jw) = zcdnc(jl,jk,jw)*zpwdw
      zndbot(jl,jk,jw) = zpwdw
    END DO ! jl
  END DO ! jk

  ! The INTENT(OUT) array psmax needs to be from last jw
  IF (jw == nwbins) THEN
    psmax(:,:) = psmax2(:,:)
  END IF

END DO !jw
!$OMP END DO
!$OMP END PARALLEL

! For the output arrays from the DO loop above, sum over all nwbins
! and store summation in (...,1) of each array. The order of summation
! alters results, and so this cannot be done within the OpenMP DO loop
! above.
!$OMP PARALLEL DEFAULT(NONE)                                                   &
!$OMP PRIVATE(jmod, jk, jl, jw)                                                &
!$OMP SHARED(kbdim, klev, nwbins, mode, modesol, topmode,                      &
!$OMP        l_fix_ukca_hygroscopicities_local,                                &
!$OMP        zndbot, zndbotm, zndtop, zndtopm)
DO jw=2,nwbins
  DO jmod=1, topmode
    IF ((mode(jmod) .AND. modesol(jmod) == 1) .OR. .NOT.                       &
        l_fix_ukca_hygroscopicities_local) THEN
!$OMP DO SCHEDULE(STATIC)
      DO jk=1, klev
        DO jl=1, kbdim
          zndtopm(jl,jk,jmod,1) = zndtopm(jl,jk,jmod,1) + zndtopm(jl,jk,jmod,jw)
          zndbotm(jl,jk,jmod,1) = zndbotm(jl,jk,jmod,1) + zndbotm(jl,jk,jmod,jw)
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF
  END DO
!$OMP DO SCHEDULE(STATIC)
  DO jk=1, klev
    DO jl=1, kbdim
      zndtop(jl,jk,1)  = zndtop(jl,jk,1)  + zndtop(jl,jk,jw)
      zndbot(jl,jk,1)  = zndbot(jl,jk,1)  + zndbot(jl,jk,jw)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO
!$OMP END PARALLEL

! Calculate the normalised expected value of CDNC for each mode
! weighted by w
DO jk=1, klev
  DO jl=1, kbdim
    DO jmod=1, topmode
      IF (mode(jmod) .AND. (modesol(jmod) == 1 .OR. .NOT.                      &
          l_fix_ukca_hygroscopicities_local)) THEN
        IF (zndbotm(jl,jk,jmod,1) < zeps) THEN
          ! Initialise to zero and print warning
          IF (PrintStatus == Prstatus_diag) THEN
            WRITE(umMessage,'(A)') 'ARG_Activate: zndbotm < zeps'
            CALL umPrint(umMessage,src='ukca_abdulrazzak_ghan')
            WRITE(umMessage,'(A,3I4,E15.6)') 'jl,jk,jmod,zndbotm= ',           &
                   jl, jk, jmod, zndbotm(jl,jk,jmod,1)
            CALL umPrint(umMessage,src='ukca_abdulrazzak_ghan')
          END IF
          pcdncactm(jl,jk,jmod) = 0.0
        ELSE
          pcdncactm(jl,jk,jmod)= zndtopm(jl,jk,jmod,1)/                        &
                                    zndbotm(jl,jk,jmod,1)
        END IF
      END IF !mode
    END DO ! jmod
    ! Calculate the normalised expected value of total CDNC
    ! weighted by w
    IF (zndbot(jl,jk,1) > zeps) THEN
      pcdncact(jl,jk)= zndtop(jl,jk,1)/zndbot(jl,jk,1)
    ELSE
      IF (PrintStatus == Prstatus_diag) THEN
        WRITE(umMessage,'(A11,I4,A1,I3,A5,E15.6)')                             &
                   'RW: zndbot(',jl,',',jk,',1)= ', zndbot(jl,jk,1)
        CALL umPrint(umMessage,src='ukca_abdulrazzak_ghan')
      END IF
    END IF
  END DO ! jl
END DO ! jk

! Calculate the normalised expected value of total CDNC
! weighted by w
! pcdncact(:,:)= zndtop(:,:,1)/zndbot(:,:,1)

! Now calculate the characteristic updraught velocity:
! i.e. find w* for which E(N_d(w))=N_d(w*)

DO jk=1, klev
  DO jl=1, kbdim
    zcdnc_tmp(:) = zcdnc(jl,jk, :)
    CALL activclosest(zcdnc_tmp, nwbins, pcdncact(jl,jk),                      &
                        zcls, zclsloc)
    pwchar(jl,jk) = pwarr(jl,jk,zclsloc)
  END DO ! jl
END DO ! jk

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_abdulrazzak_ghan
END MODULE ukca_abdulrazzak_ghan_mod
