! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
MODULE ukca_stdto3_mod
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE
!
! Description:
!  Holds standard atmosphere data used in routine initjtab
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.
!
! Number of altitudes in standard atmosphere (file or
! this module)
INTEGER, PARAMETER :: nin=44  !

CHARACTER(LEN=80) :: stdto3(nin) = [                                           &
'  0.1013E+04   295.60  0.3082E-07 ',                                          &
'  0.9161E+03   281.60  0.3082E-07 ',                                          &
'  0.7153E+03   268.80  0.3144E-07 ',                                          &
'  0.5517E+03   256.00  0.3641E-07 ',                                          &
'  0.4199E+03   242.80  0.5097E-07 ',                                          &
'  0.3147E+03   229.60  0.1044E-06 ',                                          &
'  0.2327E+03   221.80  0.2388E-06 ',                                          &
'  0.1709E+03   219.40  0.4084E-06 ',                                          &
'  0.1251E+03   217.00  0.6852E-06 ',                                          &
'  0.9139E+02   217.00  0.1268E-05 ',                                          &
'  0.6677E+02   217.00  0.2100E-05 ',                                          &
'  0.4882E+02   218.00  0.3021E-05 ',                                          &
'  0.3577E+02   220.00  0.3985E-05 ',                                          &
'  0.2628E+02   222.00  0.4897E-05 ',                                          &
'  0.1937E+02   224.00  0.5464E-05 ',                                          &
'  0.1431E+02   226.00  0.5866E-05 ',                                          &
'  0.1061E+02   229.00  0.6358E-05 ',                                          &
'  0.7899E+01   233.00  0.6831E-05 ',                                          &
'  0.5912E+01   237.00  0.7234E-05 ',                                          &
'  0.4449E+01   242.20  0.7184E-05 ',                                          &
'  0.3368E+01   247.40  0.6774E-05 ',                                          &
'  0.2566E+01   252.80  0.6049E-05 ',                                          &
'  0.1965E+01   258.40  0.5465E-05 ',                                          &
'  0.1514E+01   264.00  0.4609E-05 ',                                          &
'  0.1172E+01   266.80  0.3685E-05 ',                                          &
'  0.9088E+00   269.60  0.3011E-05 ',                                          &
'  0.7058E+00   269.00  0.2314E-05 ',                                          &
'  0.5469E+00   265.00  0.1861E-05 ',                                          &
'  0.4221E+00   261.00  0.1527E-05 ',                                          &
'  0.3242E+00   255.40  0.1302E-05 ',                                          &
'  0.2476E+00   249.80  0.1096E-05 ',                                          &
'  0.1879E+00   244.20  0.9237E-06 ',                                          &
'  0.1417E+00   238.60  0.7866E-06 ',                                          &
'  0.1062E+00   233.00  0.5897E-06 ',                                          &
'  0.7901E-01   227.80  0.3627E-06 ',                                          &
'  0.5839E-01   222.60  0.2879E-06 ',                                          &
'  0.4292E-01   220.00  0.2950E-06 ',                                          &
'  0.2400E-01   200.00  0.2444E-06 ',                                          &
'  0.1040E-01   181.00  0.16E-06   ',                                          &
'  0.4100E-02   181.00  0.11E-06   ',                                          &
'  0.1800E-02   176.70  0.90E-07   ',                                          &
'  0.7900E-03   193.00  0.80E-07   ',                                          &
'  0.1900E-03   209.20  0.70E-07   ',                                          &
'  1.0389E-05   230.90  0.60E-07   '                                           &
   ]
END MODULE ukca_stdto3_mod
