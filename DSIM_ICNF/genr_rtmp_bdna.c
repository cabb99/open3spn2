
/*********************************************************************
**********************************************************************
****        <<<<< 3SPN.2 Configuration Generator    >>>>>         ****
****                                                              ****  
****        Dan Hinckley    `                                     ****
****        Instutite for Molecular Engineering                   ****
****        University of Chicago                                 ****
****        Chicago, IL 60637                                     ****
****                                                              ****
**********************************************************************
*********************************************************************/

#include "decl_ilib.h"
#include "decl_ivar.h"
#include "decl_ifun.h"

void genr_rtmp_bdna(void)
{
    twist = 36.0;
    rise = 3.38;
    P.xval = -0.628;
    P.yval = 8.896;
    P.zval = 2.186;
    P.rval = sqrt(P.xval * P.xval + P.yval * P.yval);
    P.thta = 94.0350 * _PI_ / 180.0;
    P.mass = 94.9696;
    P.stid = 1;
    P.chrg = -0.6;
    P.chrg2 = -1.0;

    S.xval = 2.365;
    S.yval = 6.568;
    S.zval = 1.280;
    S.rval = sqrt(S.xval * S.xval + S.yval * S.yval);
    S.thta = 70.196 * _PI_ / 180;
    S.mass = 83.1104;
    S.stid = 2;
    S.chrg = 0.0;

    A.xval = 0.296;
    A.yval = 2.489;
    A.zval = 0.204;
    A.rval = sqrt(A.xval * A.xval + A.yval * A.yval);
    A.thta = 83.207 * _PI_ / 180.0;
    A.mass = 134.1220;
    A.stid = 3;
    A.stid2 = 3; // Site ID 2 is used in LAMMPS notation
    A.chrg = 0.0;

    T.xval = -0.198;
    T.yval = 3.412;
    T.zval = 0.272;
    T.rval = sqrt(T.xval * T.xval + T.yval * T.yval);
    T.thta = 93.327 * _PI_ / 180.0;
    T.mass = 125.1078;
    T.stid = 5;
    T.stid2 = 4;
    T.chrg = 0.0;

    G.xval = 0.542;
    G.yval = 2.232;
    G.zval = 0.186;
    G.rval = sqrt(G.xval * G.xval + G.yval * G.yval);
    G.thta = 76.349 * _PI_ / 180.0;
    G.mass = 150.1214;
    G.stid = 4;
    G.stid2 = 5;
    G.chrg = 0.0;

    C.xval = 0.137;
    C.yval = 3.265;
    C.zval = 0.264;
    C.rval = sqrt(C.xval * C.xval + C.yval * C.yval);
    C.thta = 87.599 * _PI_ / 180.0;
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
