
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

void wrte_rpsf(char *dnme)
{
    long i, nsit, ndna;
    char fnme[_TW_], moln[_TW_];
    FILE *fptr;

    nsit = site.totl;
    ndna = site.dna_totl;
    sprintf(fnme, "%s/in00_cvmd.psf", dnme);
    fptr = fopen(fnme, "w");

    fprintf(fptr, "PSF STANDARD\n\n");
    fprintf(fptr,"%7ld !NTITLE\n",(long)1);
    fprintf(fptr,"REMARKS .psf file generated using USER-3SPN2-JJDP\n\n");
    fprintf(fptr, "%7ld !NATOM\n", site.totl);

    sprintf(moln,"%s","DNA1 ");
    for (i = 0; i < nsit; i++)
    {
        if (comp)
        {
            if (i == ndna/2)
            {
                sprintf(moln,"%s","DNA2 ");
            }
        }
        fprintf(fptr, "%8ld ", i + 1);
        fprintf(fptr, "%s",moln);
        fprintf(fptr, "%-4ld ", atom[i].rnum + 1);
        fprintf(fptr, "%s  ", atom[i].rnme);
        fprintf(fptr, "%-5s ", atom[i].snme);
        fprintf(fptr, "%3ld  ", atom[i].stid);
        fprintf(fptr, "%13.6e   ", atom[i].chrg);
        fprintf(fptr, "%7.3lf           0\n", atom[i].mass);
    }

    fprintf(fptr, "\n%8ld !NBOND\n", cbon);
    for (i = 0; i < cbon; i++)
    {
        fprintf(fptr, "%8ld%8ld", bond[i].aste + 1, bond[i].bste + 1);
        if (i < cbon - 1)
        {
            if ((i + 1) % 4 == 0)
            {
                fprintf(fptr, "\n");
            }
        }
        else
        {
            fprintf(fptr, "\n\n");
        }
    }
    if (!sequ.nste)
    {
        fprintf(fptr, "\n");
    }

    fprintf(fptr, "%8ld !NTHETA: angles\n", cben);
    for (i = 0; i < cben; i++)
    {
        fprintf(fptr, "%8ld%8ld%8ld", bend[i].aste + 1,
            bend[i].bste + 1, bend[i].cste + 1);
        if (i < cben - 1)
        {
            if ((i + 1) % 3 == 0)
            {
                fprintf(fptr, "\n");
            }
        }
        else
        {
            fprintf(fptr, "\n\n");
        }
    }
    if (!sequ.nste)
    {
        fprintf(fptr, "\n");
    }

    fprintf(fptr, "%8ld !NPHI: dihedrals\n", ctor);
    for (i = 0; i < ctor; i++)
    {
        fprintf(fptr, "%8ld%8ld%8ld%8ld", tors[i].aste + 1,
            tors[i].bste + 1, tors[i].cste + 1, tors[i].dste + 1);
        if (i < ctor - 1)
        {
            if ((i + 1) % 2 == 0)
            {
                fprintf(fptr, "\n");
            }
        }
        else
        {
            fprintf(fptr, "\n\n");
        }
    }
    if (!sequ.nste)
    {
        fprintf(fptr, "\n");
    }

    fprintf(fptr, "%8ld !NIMPHI\n\n", (long)0);
    fprintf(fptr, "%8ld !NDON\n\n", (long)0);
    fprintf(fptr, "%8ld !NACC\n\n", (long)0);
    fprintf(fptr, "%8ld !NNB\n\n", (long)0);
    fclose(fptr);
}
