! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Field name definitions for Photolysis.
!
! Method:
!
!  This module provides field name parameters for referring to fields
!  used in Photolysis. The literals defined here will be used in communications
!  with parent applications via the Photolysis API. They are thus considered
!  part of the API definition and cannot be changed without affecting the API.
!
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA/Photolysis
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
MODULE photol_fieldname_mod

IMPLICIT NONE

PUBLIC

INTEGER, PARAMETER :: fieldname_len = 20       ! Length of env field names and
INTEGER, PARAMETER :: photol_varname_len = 10  ! names of photolysis species.

! Photolysis environment fields (driving fields), arranged by type
! Scalar - real
CHARACTER(LEN=*), PARAMETER :: fldname_equation_of_time = 'equation_of_time'
CHARACTER(LEN=*), PARAMETER :: fldname_sec_since_midnight =                    &
                                                      'sec_since_midnight'
CHARACTER(LEN=*), PARAMETER :: fldname_sin_declination = 'sin_declination'
CHARACTER(LEN=*), PARAMETER :: fldname_z_top_of_model = 'z_top_of_model'

! Flat - integer
CHARACTER(LEN=*), PARAMETER :: fldname_conv_cloud_base = 'conv_cloud_base'
CHARACTER(LEN=*), PARAMETER :: fldname_conv_cloud_top = 'conv_cloud_top'

! Flat - real
CHARACTER(LEN=*), PARAMETER :: fldname_conv_cloud_lwp = 'conv_cloud_lwp'
CHARACTER(LEN=*), PARAMETER :: fldname_cos_latitude = 'cos_latitude'
CHARACTER(LEN=*), PARAMETER :: fldname_land_fraction = 'land_fraction'
CHARACTER(LEN=*), PARAMETER :: fldname_longitude = 'longitude'
CHARACTER(LEN=*), PARAMETER :: fldname_sin_latitude = 'sin_latitude'
CHARACTER(LEN=*), PARAMETER :: fldname_surf_albedo = 'surf_albedo'
CHARACTER(LEN=*), PARAMETER :: fldname_tan_latitude = 'tan_latitude'

! Full height - real
CHARACTER(LEN=*), PARAMETER :: fldname_aod_sulph_aitk = 'aod_sulph_aitk'
CHARACTER(LEN=*), PARAMETER :: fldname_aod_sulph_accum = 'aod_sulph_accum'
CHARACTER(LEN=*), PARAMETER :: fldname_area_cloud_fraction =                   &
                                                   'area_cloud_fraction'
CHARACTER(LEN=*), PARAMETER :: fldname_conv_cloud_amount = 'conv_cloud_amount'
CHARACTER(LEN=*), PARAMETER :: fldname_ozone_mmr = 'ozone_mmr'
CHARACTER(LEN=*), PARAMETER :: fldname_p_theta_levels = 'p_theta_levels'
CHARACTER(LEN=*), PARAMETER :: fldname_qcf = 'qcf'
CHARACTER(LEN=*), PARAMETER :: fldname_qcl = 'qcl'
CHARACTER(LEN=*), PARAMETER :: fldname_rad_ctl_jo2 = 'rad_ctl_jo2'
CHARACTER(LEN=*), PARAMETER :: fldname_rad_ctl_jo2b = 'rad_ctl_jo2b'
CHARACTER(LEN=*), PARAMETER :: fldname_r_rho_levels = 'r_rho_levels'
CHARACTER(LEN=*), PARAMETER :: fldname_so4_aitken = 'so4_aitken'
CHARACTER(LEN=*), PARAMETER :: fldname_so4_accum = 'so4_accum'
CHARACTER(LEN=*), PARAMETER :: fldname_t_theta_levels = 't_theta_levels'

! Fullht plus level=0 -real
CHARACTER(LEN=*), PARAMETER :: fldname_p_layer_boundaries =                    &
                                                    'p_layer_boundaries'
CHARACTER(LEN=*), PARAMETER :: fldname_r_theta_levels = 'r_theta_levels'

! Fullht_phot real
CHARACTER(LEN=*), PARAMETER :: fldname_photol_rates_2d = 'photol_rates_2d'

END MODULE photol_fieldname_mod
