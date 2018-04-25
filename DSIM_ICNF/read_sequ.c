
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

void read_sequ(char *finp)
{
    char tlin[_TW_];
    long siz1, siz2, sqts, nrmx;
    FILE *fptr;

    nrmx = 0;
    if ((fptr = fopen(finp, "r")) != NULL)
    {
        fscanf(fptr, "%ld", &sequ.nste);
        fgets(tlin, _TW_, fptr);

        if (sequ.nste)
        {
            sequ.base = (char *)calloc(sequ.nste + 1, sizeof(char));
            if (sequ.base == NULL)
            {
                fprintf(stdout,
                    "ERROR: Cannot allocate memory for sequence #1 array.\n");
                exit(1);
            }

            seqv.base = (char *)calloc(sequ.nste + 1, sizeof(char));
            if (seqv.base == NULL)
            {
                fprintf(stdout,
                    "ERROR: Cannot allocate memory for sequence #2 array.\n");
                exit(1);
            }

            fscanf(fptr, "%s", sequ.base);
            fgets(tlin, _TW_, fptr);
            sqts = fscanf(fptr, "%s", seqv.base);
            fgets(tlin, _TW_, fptr);
            fclose(fptr);

            if (sqts == 1)
            {
                siz2 = strlen(seqv.base);
            }
            else
            {
                siz2 = 0;
            }

            seqv.nste = siz2;
            siz1 = strlen(sequ.base);

            if (siz1 != sequ.nste)
            {
                fprintf(stdout,
                    "ERROR: Sequence size fails to match sequence number.\n");
                exit(1);
            }
            if (sqts == 1)
            {
                if (siz2 != siz1)
                {
                    fprintf(stdout,
                        "ERROR: Base numbers fail to match sequences.\n");
                    exit(1);
                }
            }

            if (nrmx < sequ.nste)
            {
                nrmx = sequ.nste;
            }
            if (nrmx < seqv.nste)
            {
                nrmx = seqv.nste;
            }
        }
    }
    else
    {
        fprintf(stdout,
            "ERROR: Sequence input file cannot be opened.\n");
        exit(1);
    }

    genr_boxd(nrmx);

    site.dna1 = 0;
    if (comp)
    {
        site.dna2 = 3 * sequ.nste - 1;
        site.totl = 3 * sequ.nste * 2 - 2;
    }
    else
    {  
        site.dna2 = 3 * sequ.nste - 1;
        site.totl = 3 * sequ.nste * 1 - 1;
    }
    site.dna_totl = site.totl;

    atom =
        (struct atom_data *)calloc(2*site.totl,
        sizeof(struct atom_data));
    if (atom == NULL)
    {
        fprintf(stdout,
            "ERROR: Cannot allocate memory for atom array.\n");
        exit(1);
    }
}
