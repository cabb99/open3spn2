
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
#define mM2Ions 6.0221415E-7

void read_ions(char *finp)
{
    char tlin[_TW_];
    long i, compFactor;
    FILE *fptr;
    double boxScaling = 1.0;

    if ((fptr = fopen(finp, "r")) != NULL)
    {
        fscanf(fptr, "%lf", &conc_NaCl);    // Concentration in
        fgets(tlin, _TW_, fptr);
        fscanf(fptr, "%lf", &conc_MgCl2);   // Concentration
        fgets(tlin, _TW_, fptr);
        fscanf(fptr, "%lf", &boxScaling);
        fgets(tlin, _TW_, fptr);
        fscanf(fptr, "%lf", &ionBoxSize);   // Only used when simulating an empty box of ions
        fgets(tlin, _TW_, fptr);
        fclose(fptr);
    }
    else
    {
        fprintf(stdout,
            "ERROR: Sequence input file cannot be opened.\n");
        exit(1);
    }
    if (ions_flag == 2) {
        boxs.lenx = ionBoxSize;
        boxs.leny = ionBoxSize;
        boxs.lenz = ionBoxSize;
    } else {
        // Adjust the size of the box to ensure that a vast, unneeded amount of ions are simulated.
        boxs.lenx *= boxScaling;
        boxs.leny *= boxScaling;
        boxs.lenz *= boxScaling;
    }


    // Calculate the number of ions needed
    site.nions = 0;

    if (!comp) compFactor = 1;
    else compFactor = 2;

    // Get volume of the box
    double vBox = boxs.lenx * boxs.leny * boxs.lenz;
    long nNa, nCl, nMg;



    nNa = compFactor * (sequ.nste - 1) + (long)(mM2Ions * conc_NaCl * vBox);
    nCl = (long)(mM2Ions * conc_NaCl * vBox) + (long)(2.0 * mM2Ions * conc_MgCl2 * vBox);
    nMg = (long)(mM2Ions * conc_MgCl2 * vBox);

    // Check to see if this will make the box charge neutral
    site.extra_Na = 0;
    site.extra_Cl = 0;
    
    long total_charge = nNa - nCl + 2 * nMg - compFactor * (sequ.nste - 1);
    //printf("total_charge = %ld\n",total_charge);
    if (total_charge)
    {
        if (total_charge < 0)
        {
            // Add sodium
            for (i = 0; i < abs(total_charge); i++)
            {
                site.extra_Na++;
            }
        }
        else
        {
            // Add Chlorine
            for (i = 0; i < abs(total_charge); i++)
            {
                site.extra_Cl++;
            }

        }
    }
    // Add 2 extra sodiums if the box is only salt
    if (sequ.nste == 0) site.extra_Na += 2;

    site.nions = nNa + nCl + nMg + site.extra_Na + site.extra_Cl;

    site.dna_totl = site.totl;
    site.totl += site.nions;

    atom =
        (struct atom_data *)calloc(site.totl,
        sizeof(struct atom_data));
    if (atom == NULL)
    {
        fprintf(stdout,
            "ERROR: Cannot allocate memory for ions array.\n");
        exit(1);
    }
}
