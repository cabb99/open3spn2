
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
#define ALLW 1.15

void genr_boxd(long nrmx)
{
    double dlng, bdim;

    dlng = rise * nrmx;
    bdim = 2.0 * ALLW * dlng;
    boxs.lenx = bdim;
    boxs.leny = bdim;
    boxs.lenz = bdim;
}
