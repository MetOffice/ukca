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
! Description:
!  Routine to impose a top boundary condition for species that have
!  a thermospheric source or sink, i.e. NO, CO, H2O.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_CHEMISTRY_CTL.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_topboundary_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_TOPBOUNDARY_MOD'

CONTAINS

SUBROUTINE ukca_top_boundary(row_length, rows, model_levels, ntracers,         &
                             latitude, tracer)

USE ukca_config_specification_mod, ONLY: ukca_config, i_top_BC_H2O
USE asad_mod,               ONLY: advt, jpctr, peps
USE ukca_constants,         ONLY: pi, c_no, c_co, c_h2o, c_h2, c_h, c_oh, c_o3
USE ukca_time_mod,          ONLY: i_day_number, i_year, days_in_year

USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: ntracers
REAL, INTENT(IN) :: latitude(:)  ! Latitude of each row (degrees)
REAL, INTENT(IN OUT) :: tracer(row_length, rows, model_levels, ntracers)

! local variables
! Positions of NO and CO in UKCA tracers array
INTEGER, SAVE :: i_no  = -1
INTEGER, SAVE :: i_co  = -1
INTEGER, SAVE :: i_h2o = -1
INTEGER, SAVE :: i_h2  = -1
INTEGER, SAVE :: i_h   = -1
INTEGER, SAVE :: i_oh  = -1
INTEGER, SAVE :: i_o3  = -1

INTEGER, PARAMETER :: n_top = 4
INTEGER, PARAMETER :: j_no = 1
INTEGER, PARAMETER :: j_co = 2
INTEGER, PARAMETER :: j_h2o= 3
INTEGER, PARAMETER :: j_o3 = 4
! Number of latitudes in ACE FTS dataset
INTEGER, PARAMETER :: n_ace = 36
! Latitude spacing in ACE-FTS dataset
REAL, PARAMETER :: dellat_ace = 180.0/n_ace

LOGICAL, SAVE :: firstcall = .TRUE.
INTEGER, SAVE :: year_filled = 0

REAL, ALLOCATABLE :: clim(:,:,:)
REAL, ALLOCATABLE, SAVE :: clim_interp(:,:,:)
REAL :: acelat(n_ace)
REAL :: frac
INTEGER :: i
INTEGER :: j
INTEGER :: idx
INTEGER :: yearlen
REAL, ALLOCATABLE :: change_h2o(:)
REAL, ALLOCATABLE :: toth(:)

!ACE-FTS 3-monthly CO dataset given as constant and 1st and 2nd.
!Fourier transform coefficients.
!Olaf Morgenstern, Sept 2018
!Coeffs 1 2 3 4 5
REAL, PARAMETER :: cocoeffs(5,n_ace) = RESHAPE([                               &
    1.0343e-05, -5.7139e-06, -1.0139e-06, -3.0417e-06, -1.5135e-06,            &
    1.0351e-05, -5.7028e-06, -1.0336e-06, -3.0498e-06, -1.5180e-06,            &
    9.7896e-06, -5.6210e-06, -1.3558e-06, -2.4121e-06, -1.1959e-06,            &
    9.3485e-06, -5.5150e-06, -2.1491e-06, -1.7859e-06, -8.8550e-07,            &
    6.4714e-06, -5.7903e-06,  1.6208e-06,  6.8119e-07,  4.2493e-07,            &
    6.5282e-06, -5.2432e-06,  9.1142e-07,  5.6361e-07,  3.5104e-07,            &
    6.7975e-06, -4.6991e-06,  5.3621e-07,  1.8180e-07,  1.4572e-07,            &
    7.3118e-06, -4.4828e-06, -2.9673e-08,  2.4394e-07,  1.6882e-07,            &
    7.6905e-06, -4.0141e-06, -4.2918e-07,  1.7204e-07,  1.2271e-07,            &
    8.3275e-06, -3.7961e-06, -7.4087e-07, -1.6181e-07, -5.4085e-08,            &
    8.7063e-06, -2.9312e-06, -1.3273e-06, -5.7266e-07, -2.7991e-07,            &
    9.6419e-06, -2.2294e-06, -9.7496e-07, -1.6869e-07, -7.5792e-08,            &
    9.8245e-06, -1.7929e-06, -8.1412e-07, -8.7181e-07, -4.3839e-07,            &
    1.0560e-05, -1.1445e-06,  1.0073e-07, -2.1751e-06, -1.1017e-06,            &
    1.1308e-05, -8.0858e-07,  4.9043e-07, -2.2772e-06, -1.1528e-06,            &
    1.1124e-05, -3.1799e-07,  1.4939e-07, -2.3904e-06, -1.2196e-06,            &
    1.1294e-05,  7.5480e-08,  2.2989e-08, -2.1227e-06, -1.0878e-06,            &
    1.1266e-05,  5.2460e-09,  5.1398e-07, -2.2039e-06, -1.1230e-06,            &
    1.1335e-05, -2.2720e-07, -6.5253e-08, -2.6727e-06, -1.3676e-06,            &
    1.0933e-05, -4.5845e-08, -1.8094e-07, -2.6850e-06, -1.3770e-06,            &
    1.0745e-05,  1.0299e-07,  2.8773e-08, -2.7200e-06, -1.3940e-06,            &
    1.0618e-05,  6.9697e-07, -4.1458e-07, -2.3916e-06, -1.2368e-06,            &
    9.9231e-06,  1.2579e-06, -1.9564e-08, -2.1797e-06, -1.1292e-06,            &
    9.1729e-06,  2.0863e-06,  9.3015e-07, -9.2643e-07, -4.8428e-07,            &
    8.3174e-06,  1.9475e-06,  1.5553e-06, -3.5080e-07, -1.8076e-07,            &
    7.9472e-06,  2.9051e-06,  1.6206e-06,  4.3667e-07,  2.1392e-07,            &
    8.0176e-06,  3.8365e-06,  1.9448e-06,  4.5154e-07,  2.1614e-07,            &
    7.4619e-06,  4.5268e-06,  1.2280e-06,  7.5508e-07,  3.5646e-07,            &
    6.9411e-06,  4.8276e-06,  7.1664e-07,  9.7551e-07,  4.6044e-07,            &
    6.6307e-06,  5.5499e-06, -3.9214e-08,  1.3711e-06,  6.4714e-07,            &
    6.1959e-06,  5.5581e-06, -3.2250e-07,  1.4772e-06,  6.9812e-07,            &
    6.1178e-06,  6.0612e-06, -1.0747e-06,  1.4358e-06,  6.6319e-07,            &
    9.2436e-06,  5.1172e-06,  3.9615e-06, -1.8447e-06, -9.4921e-07,            &
    9.4305e-06,  5.1483e-06,  3.3518e-06, -2.2105e-06, -1.1440e-06,            &
    1.0374e-05,  5.4509e-06,  2.1734e-06, -3.2752e-06, -1.7061e-06,            &
    1.0585e-05,  5.5457e-06,  1.8000e-06, -3.5124e-06, -1.8329e-06],           &
    [5,n_ace])

!ACE-FTS 3-monthly NO dataset decomposed into constant, 1st and 2nd
!Fourier transform coefficients.
!Olaf Morgenstern, Sept 2018
!Latitude Coeffs 1 2 3 4
REAL, PARAMETER :: nocoeffs(5,n_ace) = RESHAPE([                               &
    3.7216e-07, -4.3727e-07, -5.5703e-08,  1.6251e-07,  8.6909e-08,            &
    3.5485e-07, -4.2486e-07, -7.9603e-08,  1.7401e-07,  9.2396e-08,            &
    3.3461e-07, -4.0067e-07, -1.2452e-07,  1.8099e-07,  9.5208e-08,            &
    2.8363e-07, -3.4596e-07, -8.8981e-08,  1.6209e-07,  8.5403e-08,            &
    4.7979e-07, -5.4090e-07,  1.7642e-07,  1.0260e-07,  5.9949e-08,            &
    4.4516e-07, -5.0789e-07,  1.1355e-07,  1.2274e-07,  6.9204e-08,            &
    2.5505e-07, -2.9443e-07,  3.2521e-08,  9.6749e-08,  5.2841e-08,            &
    1.2788e-07, -1.5069e-07, -2.8978e-08,  8.1171e-08,  4.2727e-08,            &
    9.2971e-08, -9.9975e-08, -1.7173e-08,  5.5500e-08,  2.9215e-08,            &
    9.8711e-08, -1.0974e-07, -2.0047e-08,  5.9627e-08,  3.1392e-08,            &
    5.7713e-08, -3.1688e-08,  1.2337e-09,  1.6967e-08,  9.0179e-09,            &
    4.7535e-08, -1.4523e-08,  4.5995e-10,  8.7556e-09,  4.6335e-09,            &
    4.6706e-08, -1.3386e-08,  7.6621e-09,  3.5865e-09,  2.0585e-09,            &
    4.1139e-08, -4.4998e-09, -2.0247e-09, -2.0725e-09, -1.0409e-09,            &
    4.1288e-08, -1.5117e-08,  3.0356e-09, -4.9061e-09, -2.3286e-09,            &
    4.4847e-08,  2.8043e-09,  1.0589e-08, -2.8860e-09, -1.3825e-09,            &
    4.3541e-08, -8.4126e-09,  4.8296e-10,  3.8074e-09,  2.0389e-09,            &
    4.7605e-08, -5.2713e-09, -1.3116e-09,  3.4990e-09,  1.8290e-09,            &
    5.1771e-08,  7.5440e-09, -7.9053e-09, -3.5942e-09, -2.0077e-09,            &
    4.1352e-08, -6.0758e-10, -7.5309e-09, -7.3710e-09, -3.8575e-09,            &
    4.3036e-08,  4.1618e-09, -1.1526e-08, -5.0942e-09, -2.7849e-09,            &
    4.1913e-08,  9.4329e-09, -1.9973e-08, -1.0213e-08, -5.5573e-09,            &
    4.0121e-08,  9.3359e-09, -8.6942e-09, -1.9265e-09, -1.1803e-09,            &
    4.5106e-08,  1.5444e-08, -4.6932e-09,  3.6038e-09,  1.6390e-09,            &
    4.6951e-08,  2.7356e-08,  5.5145e-10,  1.2740e-08,  6.2626e-09,            &
    4.7500e-08,  2.9097e-08,  5.8535e-09,  1.8281e-08,  9.1457e-09,            &
    6.2896e-08,  5.0146e-08,  8.3025e-09,  2.9605e-08,  1.4767e-08,            &
    8.2373e-08,  7.1688e-08,  2.9553e-09,  3.7750e-08,  1.8665e-08,            &
    1.0801e-07,  1.1840e-07,  1.2458e-08,  6.7333e-08,  3.3468e-08,            &
    1.8258e-07,  2.0183e-07, -2.3098e-08,  7.8927e-08,  3.8169e-08,            &
    2.5703e-07,  2.9406e-07, -5.0007e-08,  9.5218e-08,  4.5291e-08,            &
    2.2001e-07,  2.2272e-07, -5.9467e-08,  5.5834e-08,  2.5711e-08,            &
    2.1587e-07,  2.1169e-07, -9.0752e-09,  6.1104e-08,  2.9106e-08,            &
    2.3610e-07,  2.2014e-07,  5.0957e-08,  6.5617e-08,  3.2035e-08,            &
    1.7327e-07,  1.1810e-07,  5.4228e-08,  1.1941e-08,  5.5842e-09,            &
    1.9367e-07,  1.2734e-07,  1.8236e-08, -1.0867e-08, -6.6098e-09],           &
    [5,n_ace])

!ACE-FTS 3-monthly H2O dataset decomposed into constant, 1st and 2nd
!Fourier transform coefficients.
!Olaf Morgenstern, Sept 2018
!Latitude Coeffs 1 2 3 4
REAL, PARAMETER :: h2ocoeffs(5,n_ace) = RESHAPE([                              &
    1.8275e-06,  1.9385e-06,  2.7918e-07,  1.1339e-06,  5.6496e-07,            &
    1.7482e-06,  1.9556e-06,  1.5120e-07,  1.2145e-06,  6.0461e-07,            &
    1.9014e-06,  1.8690e-06,  4.0900e-07,  1.0202e-06,  5.0890e-07,            &
    1.9885e-06,  1.7454e-06,  6.5235e-07,  9.1125e-07,  4.5717e-07,            &
    2.2510e-06,  1.9395e-06, -5.1958e-07,  6.0419e-07,  2.8430e-07,            &
    2.1314e-06,  1.6389e-06, -1.8774e-07,  9.6654e-07,  4.7675e-07,            &
    1.9797e-06,  1.5062e-06,  9.0113e-08,  8.8843e-07,  4.4129e-07,            &
    1.9711e-06,  1.2862e-06,  5.1189e-07,  5.0692e-07,  2.5295e-07,            &
    1.8712e-06,  1.1365e-06,  5.9111e-07,  3.5170e-07,  1.7584e-07,            &
    1.7060e-06,  1.0002e-06,  6.0783e-07,  3.2339e-07,  1.6288e-07,            &
    1.6449e-06,  6.3478e-07,  5.4028e-07,  3.3542e-07,  1.7185e-07,            &
    1.4023e-06,  5.5732e-07,  3.3814e-07,  3.1593e-07,  1.6028e-07,            &
    1.1985e-06,  4.3654e-07,  8.5842e-08,  4.5344e-07,  2.2897e-07,            &
    1.1847e-06,  2.5249e-07, -3.7605e-08,  5.8919e-07,  2.9888e-07,            &
    1.1189e-06,  2.5122e-07,  6.7708e-09,  5.6701e-07,  2.8805e-07,            &
    1.0350e-06,  6.5658e-08, -5.3691e-08,  5.7691e-07,  2.9424e-07,            &
    9.6944e-07, -1.7360e-08, -5.1035e-08,  5.0968e-07,  2.6065e-07,            &
    1.0475e-06, -2.0939e-08, -5.5097e-08,  5.8445e-07,  2.9894e-07,            &
    1.0109e-06,  1.1156e-07,  2.8096e-09,  6.2913e-07,  3.2120e-07,            &
    1.0611e-06,  3.2325e-08,  1.6240e-08,  5.7072e-07,  2.9222e-07,            &
    1.1450e-06, -8.9371e-09,  2.1557e-08,  6.0580e-07,  3.1066e-07,            &
    1.1811e-06, -1.5969e-07,  2.0227e-08,  5.8695e-07,  3.0247e-07,            &
    1.3049e-06, -3.4595e-07,  1.0663e-07,  6.0265e-07,  3.1336e-07,            &
    1.4610e-06, -3.9396e-07, -9.9152e-08,  4.6947e-07,  2.4321e-07,            &
    1.6839e-06, -5.2825e-07, -3.5202e-07,  3.3374e-07,  1.7206e-07,            &
    1.7722e-06, -7.4279e-07, -5.9396e-07,  1.2294e-07,  6.3375e-08,            &
    1.7731e-06, -9.1753e-07, -6.8094e-07,  2.3005e-07,  1.1895e-07,            &
    1.9041e-06, -1.2040e-06, -7.6971e-07,  2.6305e-07,  1.3765e-07,            &
    2.0364e-06, -1.3966e-06, -6.2757e-07,  2.4001e-07,  1.2940e-07,            &
    1.9232e-06, -1.6701e-06, -1.2561e-07,  5.7363e-07,  3.0884e-07,            &
    1.8864e-06, -1.8311e-06,  1.2302e-07,  7.5913e-07,  4.0835e-07,            &
    2.0659e-06, -1.8734e-06,  5.7529e-07,  3.9619e-07,  2.2813e-07,            &
    1.7182e-06, -1.5548e-06, -7.4811e-07,  8.4559e-07,  4.3975e-07,            &
    1.7492e-06, -1.6212e-06, -8.4923e-07,  7.3101e-07,  3.8054e-07,            &
    1.5048e-06, -1.7416e-06, -4.1964e-07,  1.0097e-06,  5.2947e-07,            &
    1.4903e-06, -1.7519e-06, -3.9055e-07,  1.0293e-06,  5.3995e-07],           &
    [5,n_ace])

!ACE-FTS 3-monthly O3 dataset.
!Fourier transform coefficients to 2nd order to make climatology.
!Olaf Morgenstern, Sept 2018
REAL, PARAMETER :: o3coeffs(5,n_ace) = RESHAPE([                               &
    5.2093e-07, -3.0501e-07,  2.7602e-07, -1.3353e-07, -6.2171e-08,            &
    4.4081e-07, -2.6863e-07,  1.3607e-07, -4.4413e-08, -1.8514e-08,            &
    4.0440e-07, -2.3802e-07,  3.9693e-08, -6.1140e-09, -3.2167e-10,            &
    3.7188e-07, -1.6954e-07, -1.9669e-08, -2.4127e-08, -1.0916e-08,            &
    4.0728e-07, -2.4002e-07,  1.6186e-07, -1.6336e-08, -4.1134e-09,            &
    3.7676e-07, -1.4354e-07,  1.4899e-07, -6.8827e-08, -3.2103e-08,            &
    3.6341e-07, -1.0716e-07,  6.6159e-08, -2.0992e-08, -8.9248e-09,            &
    3.6908e-07, -1.1438e-07, -3.2458e-08,  3.7764e-08,  2.0094e-08,            &
    3.5972e-07, -8.2912e-08, -4.9043e-08,  4.0317e-08,  2.0898e-08,            &
    4.0258e-07, -1.1774e-07, -7.2229e-08,  4.3080e-08,  2.2386e-08,            &
    3.5882e-07, -5.7688e-08, -7.8390e-08, -1.1117e-08, -6.0398e-09,            &
    4.3687e-07, -1.0382e-07, -7.3913e-08, -2.2041e-08, -1.1128e-08,            &
    4.7485e-07, -6.6653e-08, -1.2075e-08, -1.2933e-07, -6.5730e-08,            &
    4.2807e-07,  1.2086e-09,  1.3166e-08, -1.3935e-07, -7.1241e-08,            &
    4.3500e-07, -8.1355e-08, -8.3265e-09, -8.5678e-08, -4.3182e-08,            &
    4.7448e-07, -7.4928e-08,  5.6054e-09, -1.3205e-07, -6.6838e-08,            &
    4.5172e-07, -4.6306e-08,  3.1445e-08, -1.4693e-08, -6.7030e-09,            &
    4.2373e-07, -4.6250e-08,  4.5800e-08, -1.0232e-08, -4.2510e-09,            &
    4.4011e-07, -3.0652e-08,  3.8499e-08, -7.2965e-08, -3.6624e-08,            &
    3.5624e-07, -6.6562e-08,  1.1966e-07, -4.2304e-08, -1.9618e-08,            &
    4.3206e-07, -4.7455e-08,  4.9168e-09, -1.3392e-07, -6.8075e-08,            &
    3.6728e-07, -3.5805e-08,  8.4126e-08, -7.9432e-08, -3.9354e-08,            &
    3.5121e-07, -4.0872e-08,  5.4060e-08, -1.2650e-07, -6.3764e-08,            &
    3.6722e-07, -4.9058e-08,  4.7771e-08, -8.2817e-08, -4.1381e-08,            &
    3.2403e-07,  4.2926e-10,  6.7643e-08, -4.0692e-08, -2.0060e-08,            &
    3.4871e-07,  5.0233e-08,  1.1321e-07,  5.8957e-08,  3.1025e-08,            &
    3.7537e-07,  6.8456e-08,  7.4805e-08,  6.4409e-08,  3.3190e-08,            &
    3.5921e-07,  7.5955e-08,  6.9903e-08,  8.6319e-08,  4.4282e-08,            &
    3.6728e-07,  9.1783e-08,  4.0672e-08,  1.2052e-07,  6.1304e-08,            &
    4.1543e-07,  1.2018e-07, -3.4157e-08,  8.7355e-08,  4.3163e-08,            &
    4.3857e-07,  1.3323e-07, -1.1552e-07,  3.1365e-08,  1.3406e-08,            &
    4.3702e-07,  1.8558e-07, -1.0531e-07,  6.7390e-08,  3.1462e-08,            &
    4.2658e-07,  1.3384e-07,  5.7944e-08,  5.5368e-08,  2.7717e-08,            &
    4.4652e-07,  1.7079e-07,  6.9623e-08,  7.5859e-08,  3.7985e-08,            &
    5.3677e-07,  2.5131e-07,  7.7937e-08,  7.1042e-08,  3.4821e-08,            &
    5.6287e-07,  2.6383e-07,  3.0677e-08,  4.1262e-08,  1.8892e-08],           &
    [5,n_ace])

CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'UKCA_TOP_BOUNDARY'
REAL(KIND=jprb)               :: zhook_handle
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! make sure we initialize the array either if the years no longer match
! i.e. the job has crossed into a new year, of at the start of each leg
IF ((firstcall) .OR. (i_year /= year_filled)) THEN
  ! determine the positions of the NO and CO tracers; this may need to be
  ! expanded.
  DO i=1,jpctr
    IF (advt(i) == 'CO        ') i_co  = i
    IF (advt(i) == 'NO        ') i_no  = i
    IF (advt(i) == 'H2O       ') i_h2o = i
    IF (advt(i) == 'H2        ') i_h2  = i
    IF (advt(i) == 'H         ') i_h   = i
    IF (advt(i) == 'OH        ') i_oh  = i
    IF (advt(i) == 'O3        ') i_o3  = i
  END DO

  DO i=1,n_ace
    acelat(i) = -90.0 + dellat_ace*(REAL(i) -0.5)
  END DO

  ! Reconstitute field
  yearlen = days_in_year(i_year)
  ! allocate climatology of CO, NO, H2O, O3
  ALLOCATE(clim(n_ace,yearlen,n_top))

  ! expand climatology from Fourier coefficients; makes annually periodic field
  DO i=1,n_ace
    DO j=1,yearlen
      clim(i,j,j_co) =cocoeffs(1,i) +                                          &
                      cocoeffs(2,i) * COS((REAL(j)-1.0)/REAL(yearlen)*2.0*pi) +&
                      cocoeffs(3,i) * SIN((REAL(j)-1.0)/REAL(yearlen)*2.0*pi) +&
                      cocoeffs(4,i) * COS((REAL(j)-1.0)/REAL(yearlen)*4.0*pi) +&
                      cocoeffs(5,i) * SIN((REAL(j)-1.0)/REAL(yearlen)*4.0*pi)
      clim(i,j,j_no) =nocoeffs(1,i) +                                          &
                      nocoeffs(2,i) * COS((REAL(j)-1.0)/REAL(yearlen)*2.0*pi) +&
                      nocoeffs(3,i) * SIN((REAL(j)-1.0)/REAL(yearlen)*2.0*pi) +&
                      nocoeffs(4,i) * COS((REAL(j)-1.0)/REAL(yearlen)*4.0*pi) +&
                      nocoeffs(5,i) * SIN((REAL(j)-1.0)/REAL(yearlen)*4.0*pi)
      clim(i,j,j_h2o)=h2ocoeffs(1,i) +                                         &
                      h2ocoeffs(2,i) * COS((REAL(j)-1.0)/REAL(yearlen)*2.0*pi)+&
                      h2ocoeffs(3,i) * SIN((REAL(j)-1.0)/REAL(yearlen)*2.0*pi)+&
                      h2ocoeffs(4,i) * COS((REAL(j)-1.0)/REAL(yearlen)*4.0*pi)+&
                      h2ocoeffs(5,i) * SIN((REAL(j)-1.0)/REAL(yearlen)*4.0*pi)
      clim(i,j,j_o3) =o3coeffs(1,i) +                                          &
                      o3coeffs(2,i) * COS((REAL(j)-1.0)/REAL(yearlen)*2.0*pi) +&
                      o3coeffs(3,i) * SIN((REAL(j)-1.0)/REAL(yearlen)*2.0*pi) +&
                      o3coeffs(4,i) * COS((REAL(j)-1.0)/REAL(yearlen)*4.0*pi) +&
                      o3coeffs(5,i) * SIN((REAL(j)-1.0)/REAL(yearlen)*4.0*pi)
    END DO
  END DO

  ! prevent negatives
  WHERE (clim < 10.0 * peps) clim = 10.0 * peps

  ! interpolate to UM grid
  IF (ALLOCATED(clim_interp)) DEALLOCATE(clim_interp)
  ALLOCATE(clim_interp(rows,yearlen,n_top))
  DO i=1,rows
    frac = 0.0
    idx = FLOOR((latitude(i) - acelat(1))/dellat_ace) + 1
    IF ((idx >= 1) .AND. (idx < n_ace)) THEN
      frac = (latitude(i) - acelat(idx))/dellat_ace
      clim_interp(i,:,:) = (1.0-frac)*clim(idx,:,:)+frac*clim(idx+1,:,:)
    ELSE IF (idx < 1) THEN
      clim_interp(i,:,:) = clim(1,:,:)
    ELSE
      clim_interp(i,:,:) = clim(n_ace,:,:)
    END IF
  END DO

  ! prevent negatives
  WHERE (clim_interp < 10.0 * peps) clim_interp = 10.0 * peps

  ! rescale to make MMR
  clim_interp(:,:,j_no) = clim_interp(:,:,j_no )*c_no
  clim_interp(:,:,j_co) = clim_interp(:,:,j_co )*c_co
  clim_interp(:,:,j_h2o)= clim_interp(:,:,j_h2o)*c_h2o
  clim_interp(:,:,j_o3) = clim_interp(:,:,j_o3 )*c_o3
  firstcall = .FALSE.
  year_filled = i_year
  DEALLOCATE(clim)
END IF

! impose upper-boundary mixing ratio for NO, CO, O3, and H2O.
DO i=1,row_length
  tracer(i,:,model_levels,i_no ) = clim_interp(:,i_day_number,j_no )
  tracer(i,:,model_levels,i_co ) = clim_interp(:,i_day_number,j_co )
  tracer(i,:,model_levels,i_o3 ) = clim_interp(:,i_day_number,j_o3 )
END DO

! Take note of change of water vapour; reduce H2, H and OH accordingly to
! preserve hydrogen
IF (ukca_config%i_ukca_topboundary == i_top_BC_H2O) THEN
  ALLOCATE(change_h2o(rows))
  ALLOCATE(toth(rows))
  DO i=1,row_length
    change_h2o = (clim_interp(:,i_day_number,j_h2o) -                          &
      tracer(i,:,model_levels,i_h2o)) / c_h2o
    tracer(i,:,model_levels,i_h2o) = clim_interp(:,i_day_number,j_h2o)
    ! calculate total hydrogen for OH, H, H2 (the other significant hydrogen
    ! compounds at 85 km).
    toth = tracer(i,:,model_levels,i_oh) *0.5 / c_oh                           &
         + tracer(i,:,model_levels,i_h ) *0.5 / c_h                            &
         + tracer(i,:,model_levels,i_h2)      / c_h2
    ! calculate rescaling factor for these compounds
    toth = MAX(1.0 - change_h2o/MAX(toth,1.0e-20),0.0)
    tracer(i,:,model_levels,i_h2) = tracer(i,:,model_levels,i_h2) * toth
    tracer(i,:,model_levels,i_h ) = tracer(i,:,model_levels,i_h ) * toth
    tracer(i,:,model_levels,i_oh) = tracer(i,:,model_levels,i_oh) * toth
  END DO
  DEALLOCATE(change_h2o)
  DEALLOCATE(toth)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_top_boundary

END MODULE ukca_topboundary_mod
