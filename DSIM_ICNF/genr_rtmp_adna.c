
/*********************************************************************
**********************************************************************
****        <<<<< MSIM :: Mesoscale Simulation Code >>>>>         ****
****                                                              ****
****        Edward J. Sambriski                                   ****
****        Postdoctoral Associate                                ****
****        Department of Chemical and Biological Engineering     ****
****        University of Wisconsin-Madison                       ****
****        Wisconsin, Madison 53706                              ****
****                                                              ****
****        Copyright (C) 2008 UW-Madison [2008.08.14]            ****
**********************************************************************
*********************************************************************/

#include "decl_ilib.h"
#include "decl_ivar.h"
#include "decl_ifun.h"

void genr_rtmp_adna(void)
{
    twist = 32.7;
    rise = 2.56;

    P.xval = 3.315; 
    P.yval = 8.372;
    P.zval = -3.935;
    P.rval = sqrt(P.xval * P.xval + P.yval * P.yval);
    //P.thta = 94.0380 * _PI_ / 180.0;
    P.thta = 68.4010 * _PI_ / 180.0;
    P.mass = 94.9696;
    P.stid = 1;
    P.chrg = -0.6;
    P.chrg2 = -1.0;

    S.xval = 6.713;
    S.yval = 6.218;
    S.zval = -2.889;
    S.rval = sqrt(S.xval * S.xval + S.yval * S.yval);
    //S.thta = 70.1971 * _PI_ / 180;
    S.thta = 42.807 * _PI_ / 180;
    S.mass = 83.1104;
    S.stid = 2;
    S.chrg = 0.0;

    A.xval = 4.995;
    A.yval = 2.338;
    A.zval = -0.875;
    A.rval = sqrt(A.xval * A.xval + A.yval * A.yval);
    //A.thta = 82.3742 * _PI_ / 180.0;
    A.thta = 25.080 * _PI_ / 180.0;
    A.mass = 134.1220;
    A.stid = 3;
    A.stid2 = 3; // Site ID 2 is used in LAMMPS notation
    A.chrg = 0.0;

    T.xval = 4.499; 
    T.yval = 3.225;
    T.zval = -0.902;
    T.rval = sqrt(T.xval * T.xval + T.yval * T.yval);
    //T.thta = 92.9990 * _PI_ / 180.0;
    T.thta = 35.634 * _PI_ / 180.0;
    T.mass = 125.1078;
    T.stid = 5;
    T.stid2 = 4;
    T.chrg = 0.0;

    G.xval = 5.240;
    G.yval = 2.087;
    G.zval = -0.812;
    G.rval = sqrt(G.xval * G.xval + G.yval * G.yval);
    //G.thta = 75.551 * _PI_ / 180.0;
    G.thta = 21.720 * _PI_ / 180.0;
    G.mass = 150.1214;
    G.stid = 4;
    G.stid2 = 5;
    G.chrg = 0.0;

    C.xval = 4.833;
    C.yval = 3.076;
    C.zval = -1.120;
    C.rval = sqrt(C.xval * C.xval + C.yval * C.yval);
    //C.thta = 88.992 * _PI_ / 180.0;
    C.thta = 32.475 * _PI_ / 180.0;
    C.mass = 110.0964;
    C.stid = 6;
    C.stid2 = 6;
    C.chrg = 0.0;

    Na.mass = 22.989769;
    Na.stid = 7;
    Na.stid2 = 15;
    Na.chrg = 1.0;

    Mg.mass = 24.305;
    Mg.stid = 8;
    Mg.stid2 = 16;
    Mg.chrg = 2.0;

    Cl.mass = 35.453;
    Cl.stid = 9;
    Cl.stid2 = 17;
    Cl.chrg = -1.0;
}
