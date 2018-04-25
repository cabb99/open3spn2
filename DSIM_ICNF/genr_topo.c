
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

int is_base(int stea) 
{
    if (!strcmp(atom[stea].snme,"A"))
    {
        return 1;
    }
    else if (!strcmp(atom[stea].snme,"T"))
    {
        return 1;
    }
    else if (!strcmp(atom[stea].snme,"G"))
    {
        return 1;
    }
    else if (!strcmp(atom[stea].snme,"C"))
    {
        return 1;
    }
    else 
    {
        return 0;
    }
}

int is_sugar(int stea) 
{
    if (!strcmp(atom[stea].snme,"S"))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int is_phosphate(int stea) 
{
    if (!strcmp(atom[stea].snme,"P"))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


void genr_topo(void)
{
    long i, ii, jj, j, j1, k1, j2, k2;

    j1 = 0;
    k1 = 0;
    k2 = 3 * sequ.nste - 1;
    cbon = sequ.nste * 3 - 2;
    if (comp)
    {
        j2 = cbon;
        cbon *= 2;
    }

    bond = (struct bond_data *)calloc(cbon, sizeof(struct bond_data));
    if (bond == NULL && ions_flag != 2)
    {
        fprintf(stdout, "ERROR: Cannot allocate memory for bonds.\n");
        exit(1);
    }

    cben = 6 * cbon;
    ctor = 3 * cbon;
    bend = (struct bend_data *)calloc(cben, sizeof(struct bend_data));
    if (bend == NULL && ions_flag != 2)
    {
        fprintf(stdout, "ERROR: Cannot allocate memory for bends.\n");
        exit(1);
    }
    tors = (struct tors_data *)calloc(ctor, sizeof(struct tors_data));
    if (tors == NULL && ions_flag != 2)
    {
        fprintf(stdout, "ERROR: Cannot allocate memory for torss.\n");
        exit(1);
    }

    for (i = 0; i < sequ.nste; i++)
    {
        if (i == 0)
        {
            bond[j1].aste = k1;
            bond[j1].bste = k1 + 1;
            j1 = j1 + 1;

            bond[j1].aste = k1;
            bond[j1].bste = k1 + 2;
            j1 = j1 + 1;

            k1 = k1 + 2;

            if (comp)
            {
                bond[j2].aste = k2;
                bond[j2].bste = k2 + 1;
                j2 = j2 + 1;

                bond[j2].aste = k2;
                bond[j2].bste = k2 + 2;
                j2 = j2 + 1;

                k2 = k2 + 2;
            }

        }
        else
        {
            bond[j1].aste = k1;
            bond[j1].bste = k1 + 1;
            j1 = j1 + 1;

            bond[j1].aste = k1 + 1;
            bond[j1].bste = k1 + 2;
            j1 = j1 + 1;

            if (comp)
            {
                bond[j2].aste = k2;
                bond[j2].bste = k2 + 1;
                j2 = j2 + 1;

                bond[j2].aste = k2 + 1;
                bond[j2].bste = k2 + 2;
                j2 = j2 + 1;
            }


            if (i != sequ.nste - 1)
            {
                bond[j1].aste = k1 + 1;
                bond[j1].bste = k1 + 3;
                j1 = j1 + 1;

                if (comp)
                {
                    bond[j2].aste = k2 + 1;
                    bond[j2].bste = k2 + 3;
                    j2 = j2 + 1;
                }
            }
            k1 = k1 + 3;
            if (comp) 
            k2 = k2 + 3;
        }
    }
    //printf("%ld\t%ld\n",j2,cbon); 
    if (comp)
    {
        if (j2 != cbon)
        {
            fprintf(stdout, "ERROR: Inconsistent number bonds.");
            exit(1);
        }
    }
    else
    {
        if (j1 != cbon)
        {
            fprintf(stdout, "ERROR: Inconsistent number bonds.");
            exit(1);
        }
    }


    j = 0;
    for (ii = 0; ii < cbon - 1; ii++)
    {
        for (jj = ii + 1; jj < cbon; jj++)
        {
            if (bond[ii].aste == bond[jj].aste)
            {
                bend[j].aste = bond[ii].bste;
                bend[j].bste = bond[ii].aste;
                bend[j].cste = bond[jj].bste;
                j++;
            }
            else if (bond[ii].aste == bond[jj].bste)
            {
                bend[j].aste = bond[ii].bste;
                bend[j].bste = bond[ii].aste;
                bend[j].cste = bond[jj].aste;
                j++;
            }
            else if (bond[ii].bste == bond[jj].aste)
            {
                bend[j].aste = bond[ii].aste;
                bend[j].bste = bond[ii].bste;
                bend[j].cste = bond[jj].bste;
                j++;
            }
            else if (bond[ii].bste == bond[jj].bste)
            {
                bend[j].aste = bond[ii].aste;
                bend[j].bste = bond[ii].bste;
                bend[j].cste = bond[jj].aste;
                j++;
            }
        }
    }

    // Loop through all of the bends and create the new stacking angle
    long id_match;
    for (i = 0; i < cben; i++)
    {
        id_match = is_base(bend[i].aste) + is_sugar(bend[i].bste) + is_phosphate(bend[i].cste);
        if (id_match == 3)
        {
            // Create angle
            bend[j].aste = bend[i].bste;
            bend[j].bste = bend[i].aste;
            bend[j].cste = bend[i].aste + 3;
            j++;
        }
    }
    cben = j;


    j = 0;
    if (dna_type == 1)  // Curved DNA
    { 
        // Removing dihedrals formed by stacking "angles" but keeping those involving only one base
        for (ii = 0; ii < cben - 1; ii++)
        {
            for (jj = ii + 1; jj < cben; jj++)
            {
                if (bend[ii].cste == bend[jj].bste)
                {
                    if (bend[ii].bste == bend[jj].aste)
                    {
                        if (!(is_base(bend[ii].bste) + is_base(bend[ii].cste)))
                        {
                            tors[j].aste = bend[ii].aste;
                            tors[j].bste = bend[ii].bste;
                            tors[j].cste = bend[ii].cste;
                            tors[j].dste = bend[jj].cste;
                            j++;
                        }
                    }
                    else if (bend[ii].bste == bend[jj].cste)
                    {
                        if (!(is_base(bend[ii].bste) + is_base(bend[ii].cste)))
                        {
                            tors[j].aste = bend[ii].aste;
                            tors[j].bste = bend[ii].bste;
                            tors[j].cste = bend[ii].cste;
                            tors[j].dste = bend[jj].aste;
                            j++;
                        }
                    }
                }
                else if (bend[ii].aste == bend[jj].bste)
                {
                    if (bend[ii].bste == bend[jj].cste)
                    {
                        if (!(is_base(bend[jj].bste) + is_base(bend[jj].cste)))
                        {
                            tors[j].aste = bend[jj].aste;
                            tors[j].bste = bend[jj].bste;
                            tors[j].cste = bend[jj].cste;
                            tors[j].dste = bend[ii].cste;
                            j++;
                        }
                    }
                    else if (bend[ii].bste == bend[jj].aste)
                    {
                        if (!(is_base(bend[jj].bste) + is_base(bend[jj].aste)))
                        {
                            tors[j].aste = bend[jj].cste;
                            tors[j].bste = bend[jj].bste;
                            tors[j].cste = bend[jj].aste;
                            tors[j].dste = bend[ii].cste;
                            j++;
                        }
                    }
                }
            }
        }
    }
    else  // Discarding all dihedrals but those that form included only sugar and phosphate sites
    {
        for (ii = 0; ii < cben - 1; ii++)
        {
            for (jj = ii + 1; jj < cben; jj++)
            {
                if (bend[ii].cste == bend[jj].bste)
                {
                    if (bend[ii].bste == bend[jj].aste)
                    {
                        if (!(is_base(bend[ii].aste) + is_base(bend[jj].cste)))
                        {
                            tors[j].aste = bend[ii].aste;
                            tors[j].bste = bend[ii].bste;
                            tors[j].cste = bend[ii].cste;
                            tors[j].dste = bend[jj].cste;
                            j++;
                        }
                    }
                    else if (bend[ii].bste == bend[jj].cste)
                    {
                        if (!(is_base(bend[ii].aste) + is_base(bend[jj].aste)))
                        {
                            tors[j].aste = bend[ii].aste;
                            tors[j].bste = bend[ii].bste;
                            tors[j].cste = bend[ii].cste;
                            tors[j].dste = bend[jj].aste;
                            j++;
                        }
                    }
                }
                else if (bend[ii].aste == bend[jj].bste)
                {
                    if (bend[ii].bste == bend[jj].cste)
                    {
                        if (!(is_base(bend[jj].aste) + is_base(bend[ii].cste)))
                        {
                            tors[j].aste = bend[jj].aste;
                            tors[j].bste = bend[jj].bste;
                            tors[j].cste = bend[jj].cste;
                            tors[j].dste = bend[ii].cste;
                            j++;
                        }
                    }
                    else if (bend[ii].bste == bend[jj].aste)
                    {
                        if (!(is_base(bend[jj].cste) + is_base(bend[ii].aste)))
                        {
                            tors[j].aste = bend[jj].cste;
                            tors[j].bste = bend[jj].bste;
                            tors[j].cste = bend[jj].aste;
                            tors[j].dste = bend[ii].cste;
                            j++;
                        }
                    }
                }
            }
        }

    }
    ctor = j;
}
