! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Computes steady-state species for Newton-Raphson integrator.
!    Part of the ASAD chemical solver.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     ASAD: ycn                      Version: steady.f 4.1 20/7/07
!
!     Purpose
!     -------
!
!  Routine to explicitly define steady state expressions - prevents
!  generality, but removes need for sluggish iteration round family
!  stuff.
!
!  Important Notes:
!     1) Ordering of calculations is important - need to avoid feedbacks!
!     2) Needs to be rewritten whenever reaction numbering is changed.
!
!  Additions:
!     To improve generality, reactions involving steady state species
!  are selected in 'setsteady' and loaded into nss* integer arrays,
!  which are then used in this routine. This causes a very slight
!  increase in CPU time, but removes the need to rewrite the routine
!  whenever a new species is added or a reaction changed.
!
!                                            Oliver   (3 Feb 1998)
! We add more general terms for the steady-state species. It is assumed
! that O(1D), O(3P), H and N can be put in steady state.
!
!
!     Method
!     ------
!
!   We sum up production and loss terms for SS species, and divide.
!   Moreover, corresponding terms in the Jacobian are calculated that
!   account for the dependence of steady-state variables on tracer
!   variables.
!
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE asad_steady_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_STEADY_MOD'

CONTAINS

SUBROUTINE asad_steady( kl )

USE asad_mod,               ONLY:  deriv, y, rk, peps,                         &
                                   nspi, nssi, nssrt, nssrx,                   &
                                   nssri, nsspt, nsspi, nsst,                  &
                                   nspo1d, nspo3, nspoh,                       &
                                   nspo3p, nsph, nuni,                         &
                                   nspho2, nspno, nspn, nss_o3p,               &
                                   nss_o1d, nss_n, nss_h,                      &
                                   o3p_in_ss, n_in_ss, h_in_ss
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook


IMPLICIT NONE


! Subroutine interface
INTEGER, INTENT(IN) :: kl            ! No. of points

! Local variables
INTEGER, PARAMETER :: n_o3=1           ! indicies for dssden etc
INTEGER, PARAMETER :: n_oh=2
INTEGER, PARAMETER :: n_ho2=3
INTEGER, PARAMETER :: n_no=4

INTEGER :: jr
INTEGER :: ix
INTEGER :: i
INTEGER :: j

REAL :: ssnum(kl)
REAL :: ssden(kl)

! Add here derivatives w.r.t. O3, OH, and HO2, and NO of numerator and
! denominator of steady state species
REAL :: dssnum(kl,4)
REAL :: dssden(kl,4)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_STEADY'

!
! Set up loops correctly
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise DERIV this timestep, first done in ASAD_INIT
deriv(:,:,:) = 1.0

! Loop through steady state species

DO ix = 1,nsst
  ssnum = 0.0
  ssden = 0.0
  dssden = 0.0
  dssnum = 0.0

  ! Production terms
  DO jr = 1,nsspt(ix)
    i = nsspi(ix,jr)
    IF (i <= nuni) THEN
      ssnum(1:kl) = ssnum(1:kl) + rk(1:kl,i)*y(1:kl,nspi(i,1))

      IF ((ix < 5) .AND. (nspi(i,1) == nspo3 ))                                &
        ! add terms to derivative for d(j[O3])/d[O3] = j_o3
        dssnum(1:kl,n_o3) = dssnum(1:kl,n_o3) + rk(1:kl,i)

      IF ((ix < 5) .AND. (nspi(i,1) == nspno ))                                &
        ! add terms to derivative for d(j[NO])/d[NO] = j_no
        dssnum(1:kl,n_no) = dssnum(1:kl,n_no) + rk(1:kl,i)

    ELSE
      ssnum(1:kl) = ssnum(1:kl) + rk(1:kl,i)*y(1:kl,nspi(i,1))*y(1:kl,nspi(i,2))
      IF (ix < 5) THEN

        ! add terms for derivative w.r.t. ozone.
        IF (nspi(i,1) == nspo1d)                                               &
          dssnum(1:kl,n_o3) = dssnum(1:kl,n_o3) +                              &
              rk(1:kl,i)*y(1:kl,nspi(i,2))*deriv(1:kl,nss_o1d,n_o3)

        IF (nspi(i,2) == nspo1d)                                               &
            dssnum(1:kl,n_o3) = dssnum(1:kl,n_o3) +                            &
              rk(1:kl,i)*y(1:kl,nspi(i,1))*deriv(1:kl,nss_o1d,n_o3)

        IF (nspi(i,1) == nspo3p)                                               &
            dssnum(1:kl,n_o3) = dssnum(1:kl,n_o3) +                            &
               rk(1:kl,i)*y(1:kl,nspi(i,2))*deriv(1:kl,nss_o3p,n_o3)

        IF (nspi(i,2) == nspo3p)                                               &
            dssnum(1:kl,n_o3) = dssnum(1:kl,n_o3) +                            &
              rk(1:kl,i)*y(1:kl,nspi(i,1))*deriv(1:kl,nss_o3p,n_o3)

        IF (nspi(i,1) == nspo3)                                                &
          dssnum(1:kl,n_o3) = dssnum(1:kl,n_o3) + rk(1:kl,i)*y(1:kl,nspi(i,2))

        IF (nspi(i,2) == nspo3)                                                &
          dssnum(1:kl,n_o3) = dssnum(1:kl,n_o3) + rk(1:kl,i)*y(1:kl,nspi(i,1))

        ! add terms for derivative w.r.t OH
        IF (nspi(i,1) == nspo3p)                                               &
          ! add terms to derivates for d(a[A][B])
          dssnum(1:kl,n_oh) = dssnum(1:kl,n_oh) +                              &
            rk(1:kl,i)*y(1:kl,nspi(i,2))*deriv(1:kl,nss_o3p,n_oh)

        IF (nspi(i,2) == nspo3p)                                               &
          ! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
          dssnum(1:kl,n_oh) = dssnum(1:kl,n_oh) +                              &
            rk(1:kl,i)*y(1:kl,nspi(i,1))*deriv(1:kl,nss_o3p,n_oh)

        IF (nspi(i,1) == nspoh)                                                &
          ! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
          dssnum(1:kl,n_oh) = dssnum(1:kl,n_oh) + rk(1:kl,i)*y(1:kl,nspi(i,2))

        IF (nspi(i,2) == nspoh)                                                &
          ! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
          dssnum(1:kl,n_oh) = dssnum(1:kl,n_oh) + rk(1:kl,i)*y(1:kl,nspi(i,1))

        ! add terms for derivative w.r.t HO2
        IF (nspi(i,1) == nspo3p)                                               &
          dssnum(1:kl,n_ho2) = dssnum(1:kl,n_ho2) +                            &
            rk(1:kl,i)*y(1:kl,nspi(i,2))*deriv(1:kl,nss_o3p,n_ho2)

        IF (nspi(i,2) == nspo3p)                                               &
          dssnum(1:kl,n_ho2) = dssnum(1:kl,n_ho2) +                            &
            rk(1:kl,i)*y(1:kl,nspi(i,1))*deriv(1:kl,nss_o3p,n_ho2)

        IF (nspi(i,1) == nspho2)                                               &
          dssnum(1:kl,n_ho2) = dssnum(1:kl,n_ho2) + rk(1:kl,i)*y(1:kl,nspi(i,2))

        IF (nspi(i,2) == nspho2)                                               &
          dssnum(1:kl,n_ho2) = dssnum(1:kl,n_ho2) + rk(1:kl,i)*y(1:kl,nspi(i,1))

        ! add terms for derivative w.r.t NO
        IF (nspi(i,1) == nspno)                                                &
          dssnum(1:kl,n_no) = dssnum(1:kl,n_no) + rk(1:kl,i)*y(1:kl,nspi(i,2))

        IF (nspi(i,2) == nspno)                                                &
          dssnum(1:kl,n_no) = dssnum(1:kl,n_no) + rk(1:kl,i)*y(1:kl,nspi(i,1))

      END IF
    END IF
  END DO ! jr
  !
  ! Destruction terms
  DO jr = 1,nssrt(ix)
    i = nssri(ix,jr)
    j = nssrx(ix,jr)
    IF (i <= nuni) THEN
      ssden(1:kl) = ssden(1:kl) + rk(1:kl,i)
    ELSE
      ssden(1:kl) = ssden(1:kl) + rk(1:kl,i) * y(1:kl,nspi(i,j))
      IF (ix < 5) THEN
        IF (nspi(i,j) == nspo3 )                                               &
          dssden(1:kl,n_o3 ) = dssden(1:kl,n_o3 ) + rk(1:kl,i)
        IF (nspi(i,j) == nspoh )                                               &
          dssden(1:kl,n_oh ) = dssden(1:kl,n_oh ) + rk(1:kl,i)
        IF (nspi(i,j) == nspho2)                                               &
          dssden(1:kl,n_ho2) = dssden(1:kl,n_ho2) + rk(1:kl,i)
        IF (nspi(i,j) == nspno )                                               &
          dssden(1:kl,n_no ) = dssden(1:kl,n_no ) + rk(1:kl,i)
      END IF
    END IF
  END DO ! jr
  !
  ! Steady state and derivatives of steady state
  y(1:kl,nssi(ix)) = ssnum(1:kl)/ssden(1:kl)
  IF (ix < 5) THEN
    DO jr =1,4
      deriv(1:kl,ix,jr) =                                                      &
        (ssden(1:kl)*dssnum(1:kl,jr) - ssnum(1:kl)*dssden(1:kl,jr)) /          &
        (ssden(1:kl)*ssden(1:kl))
    END DO ! jr
  END IF
END DO ! ix

! rescale deriv to mean [O3]/[O] * d[O]/d[O3], where [O] = [O(1D)] or [O(3P)]
! for O(1D), and O(3P), N and H when these are SS species

WHERE (y(1:kl,nspo1d) > peps)
  deriv(1:kl,nss_o1d,n_o3) = deriv(1:kl,nss_o1d,n_o3 )*                        &
                            y(1:kl,nspo3 )/y(1:kl,nspo1d)
  deriv(1:kl,nss_o1d,n_oh) = deriv(1:kl,nss_o1d,n_oh )*                        &
                            y(1:kl,nspoh )/y(1:kl,nspo1d)
  deriv(1:kl,nss_o1d,n_ho2)= deriv(1:kl,nss_o1d,n_ho2)*                        &
                            y(1:kl,nspho2)/y(1:kl,nspo1d)
  deriv(1:kl,nss_o1d,n_no )= deriv(1:kl,nss_o1d,n_no )*                        &
                            y(1:kl,nspno )/y(1:kl,nspo1d)
ELSE WHERE
  deriv(1:kl,nss_o1d,n_o3)  = 1.0
  deriv(1:kl,nss_o1d,n_oh)  = 1.0
  deriv(1:kl,nss_o1d,n_ho2) = 1.0
  deriv(1:kl,nss_o1d,n_no)  = 1.0
END WHERE

IF (o3p_in_ss) THEN
  WHERE (y(1:kl,nspo3p) > peps)
    deriv(1:kl,nss_o3p,n_o3 )= deriv(1:kl,nss_o3p,n_o3 )*                      &
                            y(1:kl,nspo3 )/y(1:kl,nspo3p)
    deriv(1:kl,nss_o3p,n_oh )= deriv(1:kl,nss_o3p,n_oh )*                      &
                            y(1:kl,nspoh )/y(1:kl,nspo3p)
    deriv(1:kl,nss_o3p,n_ho2)= deriv(1:kl,nss_o3p,n_ho2)*                      &
                            y(1:kl,nspho2)/y(1:kl,nspo3p)
    deriv(1:kl,nss_o3p,n_no )= deriv(1:kl,nss_o3p,n_no )*                      &
                            y(1:kl,nspno )/y(1:kl,nspo3p)
  ELSE WHERE
    deriv(1:kl,nss_o3p,n_o3)  = 1.0
    deriv(1:kl,nss_o3p,n_oh)  = 1.0
    deriv(1:kl,nss_o3p,n_ho2) = 1.0
    deriv(1:kl,nss_o3p,n_no)  = 1.0
  END WHERE
END IF

IF (n_in_ss) THEN
  WHERE (y(1:kl,nspn  ) > peps)
    deriv(1:kl,nss_n,n_o3) = deriv(1:kl,nss_n,n_o3 )*                          &
                          y(1:kl,nspo3 )/y(1:kl,nspn)
    deriv(1:kl,nss_n,n_oh) = deriv(1:kl,nss_n,n_oh )*                          &
                          y(1:kl,nspoh )/y(1:kl,nspn)
    deriv(1:kl,nss_n,n_ho2)= deriv(1:kl,nss_n,n_ho2)*                          &
                          y(1:kl,nspho2)/y(1:kl,nspn)
    deriv(1:kl,nss_n,n_no )= deriv(1:kl,nss_n,n_no )*                          &
                          y(1:kl,nspno )/y(1:kl,nspn)
  ELSE WHERE
    deriv(1:kl,nss_n,n_o3)  = 1.0
    deriv(1:kl,nss_n,n_oh)  = 1.0
    deriv(1:kl,nss_n,n_ho2) = 1.0
    deriv(1:kl,nss_n,n_no)  = 1.0
  END WHERE
END IF

IF (h_in_ss) THEN
  WHERE (y(1:kl,nsph  ) > peps)
    deriv(1:kl,nss_h,n_o3) = deriv(1:kl,nss_h,n_o3 )*                          &
                        y(1:kl,nspo3 )/y(1:kl,nsph  )
    deriv(1:kl,nss_h,n_oh) = deriv(1:kl,nss_h,n_oh )*                          &
                        y(1:kl,nspoh )/y(1:kl,nsph  )
    deriv(1:kl,nss_h,n_ho2)= deriv(1:kl,nss_h,n_ho2)*                          &
                        y(1:kl,nspho2)/y(1:kl,nsph  )
    deriv(1:kl,nss_h,n_no )= deriv(1:kl,nss_h,n_no )*                          &
                        y(1:kl,nspno )/y(1:kl,nsph  )
  ELSE WHERE
    deriv(1:kl,nss_h,n_o3)  = 1.0
    deriv(1:kl,nss_h,n_oh)  = 1.0
    deriv(1:kl,nss_h,n_ho2) = 1.0
    deriv(1:kl,nss_h,n_no)  = 1.0
  END WHERE
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_steady
END MODULE asad_steady_mod
