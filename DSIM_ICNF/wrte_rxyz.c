
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

void wrte_rxyz(char *dnme)
{
    long i, ncyc;
    char fnme[_TW_];
    FILE *fptr;

    sprintf(fnme, "%s/in00_conf.xyz", dnme);
    fptr = fopen(fnme, "w");

    ncyc = site.totl;
    fprintf(fptr, "%ld\n", ncyc);
    fprintf(fptr, "\n");

    for (i = 0; i < site.totl; i++)
    {
        fprintf(fptr, "%4s  ", atom[i].snme);
        fprintf(fptr, "%12.6lf  ", atom[i].xval);
        fprintf(fptr, "%12.6lf  ", atom[i].yval);
        fprintf(fptr, "%12.6lf\n", atom[i].zval);
    }
    fclose(fptr);
}
