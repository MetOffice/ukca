! *****************************COPYRIGHT*******************************
! (c) [University of California] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution.
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Fast-jx is a routine for calculating online photolysis rates
!   This module stores various parameters and arrays required for the
!   rest of the Fast-JX routines. It is very much based upon the fastj_data
!   module, with modifications due to  differences between fast-j and fast-Jx.
!   This version is a preliminary one to allow scattering phase factors etc
!   to be read in so doesn't include all variables as yet.
!   And definitely needs tidying up!
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE fastjx_data

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim
USE asad_mod,            ONLY: jppj
USE ukca_parpho_mod,     ONLY: jpwav
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
PUBLIC
SAVE

! Setup Blocking of columns for complete fast-jx model
! For time being only tested blocking mode 0

!      Blocking mode          1: rows
!                             2: whole domain
!                             3: compressed   (Not implemented yet)
!                             0: point by point

INTEGER :: Blocking_Mode=0

CHARACTER(LEN=errormessagelength) :: cmessage   ! Error message

! ***********************************************************************
       ! Defaults and Parameters

       ! Height of the atmosphere
REAL, PARAMETER               :: zzht = 5.0e5
! MNaximum solar zenith angle
REAL, PARAMETER               :: szamax = 98.0e0

! Retain these from fast-j implementation
! Inelegant to duplicate rows/ row lengths?
INTEGER                       :: ipcx      ! Number of longitudes
INTEGER                       :: jpcx      ! Number of latitudes
INTEGER                       :: lpar      ! Number of levels
INTEGER                       :: jpcl      ! Number of chem. levels

INTEGER                       :: kpcx      ! No. points (xy) in photoj
INTEGER                       :: jjpnl     ! No. points (xyz)in photoj

REAL                          :: tau       ! Time in hours
INTEGER                       :: daynumber ! Day of the year

! SW band in UM at which we calculate aerosol OD
INTEGER, PARAMETER            :: aer_swband = 4

! No. SW bands for aerosol. (Was 4 in Fastj-j i.e. NK=4 in fastj_specs)
! Increased in Fast-jx to 5 for proper treatment of Pinatubo aerosol
INTEGER, PARAMETER            :: sw_band_aer = 5
! No phases used in mie scattering
INTEGER, PARAMETER            :: sw_phases  = 8

!      Variables taken from parm_CTM (for fast-JX code v5.3+)

INTEGER, PARAMETER :: wx_ = 18   ! Max number of wavelength bins
INTEGER, PARAMETER :: x_  = 100  ! Max no. of cross sections to be read
INTEGER, PARAMETER :: a_ = 40    ! Max no. of aerosol/cloud mie data sets

! Variables taken from parm_MIE (for fast-JX code v5.3+)

! Number of wavelength bins (can choose smaller for speed)
! Set using RUN_UKCA namelist
INTEGER            :: w_

! Number of Mie Scattering arrays
! = 4*LPAR + 5 + 2*sum(JADDLV)
INTEGER, PARAMETER :: n_  = 501

! Number of Gauss points used. HAS to be 4 in fast-JX
INTEGER, PARAMETER :: m_  = 4
INTEGER, PARAMETER :: m2_ = 2*m_

! 4 Gauss pts = 8-stream
REAL, PARAMETER  :: emu(m_) = [.06943184420297e0,                              &
   .33000947820757e0, .66999052179243e0, .93056815579703e0]
REAL, PARAMETER  :: wt(m_)  = [.17392742256873e0,                              &
   .32607257743127e0, .32607257743127e0, .17392742256873e0]

! length of average solar cycle
INTEGER, PARAMETER :: n_solcyc_av = 128
! total number of months in solar cycle data
INTEGER :: n_solcyc_ts

! Number of cloud/aerosol types for scattering
! At present set to 2 (water/ice clouds)
! Will need to increase to include aerosols
INTEGER,PARAMETER             :: mx =  3    !kk look at set_aer!

! Index within the loaded array of scattering types
INTEGER                       :: miedx(mx)

! Arrays that depend on GCM parameters

      ! Array containing quantum yields
REAL,  ALLOCATABLE        :: jfacta(:)
INTEGER, ALLOCATABLE     :: jind(:)

! Array containing labels for photolysed species    ! was len=10!
CHARACTER(LEN=7), ALLOCATABLE     :: jlabel(:)

! Array containing indices in domain (xy) arrays
INTEGER, ALLOCATABLE     :: nsl(:,:)

! air mass fraction
REAL, ALLOCATABLE        :: amf2(:,:,:)

REAL, ALLOCATABLE        :: od_3d(:,:,:,:)
REAL, ALLOCATABLE        :: od_block(:,:,:)
REAL, ALLOCATABLE        :: od600(:,:)

! ozone (mol/m2)
REAL, ALLOCATABLE        :: dm_3d(:,:,:)
! dm(blocked)
REAL, ALLOCATABLE        :: dm_block(:,:)

! ozone (mol/m2)
REAL, ALLOCATABLE        :: o3_3d(:,:,:)
! o3 (blocked)
REAL, ALLOCATABLE        :: o3_block(:,:)


REAL, ALLOCATABLE        :: od(:,:,:)         ! Optical Depth
REAL, ALLOCATABLE        :: odw_3d(:,:,:)     ! OD for water
REAL, ALLOCATABLE        :: odi_3d(:,:,:)     ! OD for ice
REAL, ALLOCATABLE        :: odw_block(:,:)    ! OD for water
REAL, ALLOCATABLE        :: odi_block(:,:)    ! OD for ice

REAL, ALLOCATABLE        :: ods_3d(:,:,:)     ! OD for sulfate aerosol
REAL, ALLOCATABLE        :: ods_block(:,:)    ! OD for sulfate aerosol (blocked)

! Array containing heights of box edges
! Also contains surface pressure and TOA pressure (=0)
REAL, ALLOCATABLE        :: pz_all(:,:,:)
! pz all (blocked)
REAL, ALLOCATABLE        :: pz_block(:,:)

! Pressure at box edges (excludes TOA)
REAL, ALLOCATABLE        :: pz_3d(:,:,:)

! Array containing heights of box edges
REAL, ALLOCATABLE        :: rz_3d(:,:,:)

REAL, ALLOCATABLE        :: rz_all(:,:,:)
REAL, ALLOCATABLE        :: rz_block(:,:)

REAL, ALLOCATABLE        :: re_all(:,:,:)

! Surface albedo
REAL, ALLOCATABLE        :: sa_2d(:,:)
REAL, ALLOCATABLE        :: sa_block(:)

! Solar zenith angle (blocked)
REAL, ALLOCATABLE        :: sza(:)
REAL, ALLOCATABLE        :: szafac(:)

! Solar zenith angle (2D)
REAL, ALLOCATABLE        :: sza_2d(:,:)
REAL, ALLOCATABLE        :: szafac_2d(:,:)

! cosine of solar zenith angle (blocked)
REAL, ALLOCATABLE        :: u0(:)

! Temperature on model levels
REAL, ALLOCATABLE        :: tz_3d(:,:,:)
! Tz (blocked)
REAL, ALLOCATABLE        :: tz_block(:,:)


REAL :: daa(a_)                         ! Density of scattering type
REAL :: fl(wx_)     ! Top of atmosphere solar flux read from spec file
REAL :: fl_cyc(wx_) ! solar cycle component of TOA flux (if i_ukca_solcyc>0)
REAL :: paa(sw_phases, sw_band_aer,a_)  ! phases of scattering types
REAL :: qaa(sw_band_aer,a_)             ! `Q' of scattering types

! Photolysis rates for O2, O3 and O1D given at 3 different temperatures
REAL :: qo2(wx_,3)
REAL :: qo3(wx_,3)
REAL :: q1d(wx_,3)

! Photolysis rates for all other species
REAL :: qqq(wx_,2,x_)

! Solar cycle variables
REAL              :: solcyc_spec(wx_)       ! spectral component of solar cyle
REAL              :: solcyc_quanta(jpwav)   ! quanta component of solar cycle
REAL, ALLOCATABLE :: solcyc_ts(:)           ! obs. time series of solar cyle
REAL              :: solcyc_av(n_solcyc_av) ! average solar cycle

! Rayleigh parameters (effective cross-section) (cm2)
REAL :: qrayl(wx_+1)

! Effective radius of scattering type
REAL :: raa(a_)

! Single scattering albedos
REAL :: saa(sw_band_aer,a_)

! Temperatures read from file to allow interpolation
REAL :: tqq(3, x_)

! Wavelengths at which scattering co-efficents are calculated
REAL :: waa(sw_band_aer,a_)

! Effective wavelengths
REAL :: wl(wx_)

! Cloud sub layers factor and minimum dtau
REAL :: atau
REAL :: atau0

! SIZES

      ! Maximum number of cloud sub layers
INTEGER      :: jtaumx

! Number of aerosol/ cloud data types (must be less than a_)
INTEGER      :: naa

! Number of species to read x-sections for (must be less than x_)
INTEGER      :: njval

! Maximum and minimum wavelngth bins (set from file)
INTEGER      :: nw1
INTEGER      :: nw2



! **********************************************************************
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FASTJX_DATA'

CONTAINS


! ######################################################################
! subroutine that will set limits on lon, lat etc.
! Based on approach of fast-j routines
! Included as a dummy routine at present.
SUBROUTINE fastjx_set_limits (row_length, rows, model_levels,                  &
                              chemlev)

USE ereport_mod, ONLY: ereport
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim

IMPLICIT  NONE

INTEGER,INTENT(IN)            :: row_length        ! Model dimensions
INTEGER,INTENT(IN)            :: rows
INTEGER,INTENT(IN)            :: model_levels
INTEGER,INTENT(IN),OPTIONAL   :: chemlev

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_SET_LIMITS'


! ***********************************************************************
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! load size of GCM patch into module
ipcx   = row_length
jpcx   = rows
lpar   = model_levels

IF (PRESENT(chemlev)) THEN
  jpcl = chemlev
ELSE
  jpcl = model_levels
END IF

! Construct arrays depending on how they are passed to fast-jx
SELECT CASE (Blocking_Mode)
CASE (0)
  kpcx  = 1                    ! Blocking row-wise
  jjpnl  = jpcl*kpcx
CASE (1)
  kpcx  = ipcx                 ! Blocking row-wise
  jjpnl  = jpcl*kpcx
CASE (2)
  kpcx  = ipcx*jpcx            ! Blocking domain
  jjpnl  = jpcl*kpcx
CASE DEFAULT
  cmessage='Blocking Mode does not Exist'
  CALL ereport('FASTJX_DATA.FASTJX_SET_LIMITS', Blocking_Mode, cmessage)
END SELECT


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_set_limits

! ######################################################################
      ! subroutine that allocates arrays that depend on GCM
      ! Based on approach of fast-j routines
SUBROUTINE fastjx_allocate_memory

IMPLICIT  NONE

LOGICAL :: first=.TRUE.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_ALLOCATE_MEMORY'


! ***********************************************************************
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate arrays that depend on number of species
IF (first) THEN
  ALLOCATE(jind(jppj))
  ALLOCATE(jlabel(jppj))
  ALLOCATE(jfacta(jppj))
  first=.FALSE.
END IF

! Allocate arrays that depend on GCM size
IF (.NOT. ALLOCATED(amf2))                                                     &
  ALLOCATE(amf2      ((2*(lpar+1)+1), (2*(lpar+1)+1), kpcx))
IF (.NOT. ALLOCATED(dm_3d))                                                    &
  ALLOCATE(dm_3d     (ipcx, jpcx, (lpar+1)))
IF (.NOT. ALLOCATED(dm_block))                                                 &
  ALLOCATE(dm_block  (kpcx, (lpar+1)))
IF (.NOT. ALLOCATED(o3_3d     ))                                               &
  ALLOCATE(o3_3d     (ipcx, jpcx, (lpar+1)))
IF (.NOT. ALLOCATED(o3_block  ))                                               &
  ALLOCATE(o3_block  (kpcx, (lpar+1)))
IF (.NOT. ALLOCATED(nsl       ))                                               &
  ALLOCATE(nsl       (2,    kpcx))
IF (.NOT. ALLOCATED(odw_3d    ))                                               &
  ALLOCATE(odw_3d    (ipcx, jpcx, (lpar+1)))
IF (.NOT. ALLOCATED(odi_3d    ))                                               &
  ALLOCATE(odi_3d    (ipcx, jpcx, (lpar+1)))
IF (.NOT. ALLOCATED(odw_block ))                                               &
  ALLOCATE(odw_block (kpcx, (lpar+1)))
IF (.NOT. ALLOCATED(odi_block ))                                               &
  ALLOCATE(odi_block (kpcx, (lpar+1)))
IF (.NOT. ALLOCATED(ods_3d ))                                                  &
  ALLOCATE(ods_3d    (ipcx, jpcx, (lpar+1)))
IF (.NOT. ALLOCATED(ods_block  ))                                              &
 ALLOCATE(ods_block (kpcx, (lpar+1)))
IF (.NOT. ALLOCATED(od_block  ))                                               &
  ALLOCATE(od_block  (kpcx, sw_band_aer, (lpar+1)))
IF (.NOT. ALLOCATED(od600     ))                                               &
  ALLOCATE(od600     (kpcx, (lpar+1)))
IF (.NOT. ALLOCATED(pz_all    ))                                               &
  ALLOCATE(pz_all    (ipcx, jpcx, (lpar+2)))
IF (.NOT. ALLOCATED(pz_block  ))                                               &
  ALLOCATE(pz_block  (kpcx, (lpar+2)))
IF (.NOT. ALLOCATED(pz_3d     ))                                               &
  ALLOCATE(pz_3d     (ipcx, jpcx, (lpar+1)))
IF (.NOT. ALLOCATED(re_all    ))                                               &
  ALLOCATE(re_all    (ipcx, jpcx, (lpar+2)))
IF (.NOT. ALLOCATED(rz_all    ))                                               &
  ALLOCATE(rz_all    (ipcx, jpcx, (lpar+2)))
IF (.NOT. ALLOCATED(rz_block  ))                                               &
  ALLOCATE(rz_block  (kpcx, (lpar+2)))
IF (.NOT. ALLOCATED(rz_3d     ))                                               &
  ALLOCATE(rz_3d     (ipcx, jpcx, (lpar+1)))
IF (.NOT. ALLOCATED(sa_2d     ))                                               &
  ALLOCATE(sa_2d     (ipcx, jpcx))
IF (.NOT. ALLOCATED(sa_block  ))                                               &
  ALLOCATE(sa_block  (kpcx))
IF (.NOT. ALLOCATED(sza       ))                                               &
  ALLOCATE(sza       (kpcx ))
IF (.NOT. ALLOCATED(szafac    ))                                               &
  ALLOCATE(szafac    (kpcx ))
IF (.NOT. ALLOCATED(sza_2d    ))                                               &
  ALLOCATE(sza_2d    (ipcx, jpcx))
IF (.NOT. ALLOCATED(szafac_2d ))                                               &
  ALLOCATE(szafac_2d (ipcx, jpcx))
IF (.NOT. ALLOCATED(tz_3d     ))                                               &
  ALLOCATE(tz_3d     (ipcx, jpcx, (lpar+1)))
IF (.NOT. ALLOCATED(tz_block  ))                                               &
  ALLOCATE(tz_block  (kpcx, (lpar+1)))
IF (.NOT. ALLOCATED(u0        ))                                               &
  ALLOCATE(u0        (kpcx ))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_allocate_memory

! ######################################################################
      ! subroutine that deallocates arrays that depend on GCM
      ! Based on approach of fast-j routines
      ! Dummy routine at the moment
SUBROUTINE fastjx_deallocate_memory

IMPLICIT  NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_DEALLOCATE_MEMORY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(amf2))        DEALLOCATE(amf2)
IF (ALLOCATED(dm_3d))       DEALLOCATE(dm_3d)
IF (ALLOCATED(dm_block))    DEALLOCATE(dm_block)
IF (ALLOCATED(o3_3d))       DEALLOCATE(o3_3d)
IF (ALLOCATED(o3_block))    DEALLOCATE(o3_block)
IF (ALLOCATED(nsl))         DEALLOCATE(nsl)
IF (ALLOCATED(od))          DEALLOCATE(od)
IF (ALLOCATED(od600))       DEALLOCATE(od600)
IF (ALLOCATED(od_3d))       DEALLOCATE(od_3d)
IF (ALLOCATED(od_block))    DEALLOCATE(od_block)
IF (ALLOCATED(odi_3d))      DEALLOCATE(odi_3d)
IF (ALLOCATED(ods_block))   DEALLOCATE(ods_block)
IF (ALLOCATED(ods_3d))      DEALLOCATE(ods_3d)
IF (ALLOCATED(odw_3d))      DEALLOCATE(odw_3d)
IF (ALLOCATED(odi_block))   DEALLOCATE(odi_block)
IF (ALLOCATED(odw_block))   DEALLOCATE(odw_block)
IF (ALLOCATED(pz_all))      DEALLOCATE(pz_all)
IF (ALLOCATED(pz_block))    DEALLOCATE(pz_block)
IF (ALLOCATED(pz_3d))       DEALLOCATE(pz_3d)
IF (ALLOCATED(re_all))      DEALLOCATE(re_all)
IF (ALLOCATED(rz_all))      DEALLOCATE(rz_all)
IF (ALLOCATED(rz_block))    DEALLOCATE(rz_block)
IF (ALLOCATED(rz_3d))       DEALLOCATE(rz_3d)
IF (ALLOCATED(sa_2d))       DEALLOCATE(sa_2d)
IF (ALLOCATED(sa_block))    DEALLOCATE(sa_block)
IF (ALLOCATED(sza))         DEALLOCATE(sza)
IF (ALLOCATED(szafac))      DEALLOCATE(szafac)
IF (ALLOCATED(sza_2d))      DEALLOCATE(sza_2d)
IF (ALLOCATED(szafac_2d))   DEALLOCATE(szafac_2d)
IF (ALLOCATED(u0))          DEALLOCATE(u0)
IF (ALLOCATED(tz_3d))       DEALLOCATE(tz_3d)
IF (ALLOCATED(tz_block))    DEALLOCATE(tz_block)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_deallocate_memory

! ######################################################################
! Function to do 3 point linear interpolation
!-----------------------------------------------------------------------
REAL FUNCTION flint (tint, t1, t2, t3, f1, f2, f3)
IMPLICIT NONE

REAL, INTENT(IN) :: tint
REAL, INTENT(IN) :: t1
REAL, INTENT(IN) :: t2
REAL, INTENT(IN) :: t3
REAL, INTENT(IN) :: f1
REAL, INTENT(IN) :: f2
REAL, INTENT(IN) :: f3

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FLINT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (tint <= t2) THEN
  IF (tint <= t1) THEN
    flint = f1
  ELSE
    flint = f1 + (f2 - f1)*(tint - t1)/(t2 -t1)
  END IF
ELSE
  IF (tint >= t3) THEN
    flint = f3
  ELSE
    flint = f2 + (f3 - f2)*(tint - t2)/(t3 -t2)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END FUNCTION flint
!########################################################################

END MODULE fastjx_data
