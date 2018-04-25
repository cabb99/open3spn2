
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

void genr_trns(void)
{
    long i;
    double xcom, ycom, zcom, mcom;

    xcom = 0.0;
    ycom = 0.0;
    zcom = 0.0;
    mcom = 0.0;
    for (i = 0; i < site.dna_totl; i++)
    {
        xcom += atom[i].xval * atom[i].mass;
        ycom += atom[i].yval * atom[i].mass;
        zcom += atom[i].zval * atom[i].mass;
        mcom += atom[i].mass;
    }
    xcom /= mcom;
    ycom /= mcom;
    zcom /= mcom;

    for (i = 0; i < site.dna_totl; i++)
    {
        atom[i].xval -= xcom;
        atom[i].yval -= ycom;
        atom[i].zval -= zcom;
    }
}
