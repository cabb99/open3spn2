
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
#define mM2Ions 6.0221415E-7

// Apply periodic boundary conditions
void func_pbcs(double *delr, double boxd)
{
    if (*delr > 0.5 * boxd)
    {
        *delr = *delr - boxd;
    }
    else if (*delr < - 0.5 * boxd)
    {
        *delr = *delr - boxd;
    }
}

void genr_ions(void)
{
    long i, flag, nIonsAdded, j;
    double rx, ry, rz, dx, dy, dz,dr;
    srand(time(NULL)); // Seed the random number generator

    // Determining the amount of each ion
    double vBox = boxs.lenx * boxs.leny * boxs.lenz;
    long nNa, nCl, nMg;
    long ndna = site.dna_totl;
    if (!comp) nNa = 1 * (sequ.nste - 1) + (long)(mM2Ions * conc_NaCl * vBox) + site.extra_Na;
    else nNa = 2 * (sequ.nste - 1) + (long)(mM2Ions * conc_NaCl * vBox) + site.extra_Na;
    nCl = (long)(mM2Ions * conc_NaCl * vBox) + (long)(2.0 * mM2Ions * conc_MgCl2 * vBox) + site.extra_Cl;
    nMg = (long)(mM2Ions * conc_MgCl2 * vBox);

    double dcut = 5.0;

    // Adding the sodium
    nIonsAdded = 0;
    long offset = ndna;
    for (i = 0; i < nNa; i++)
    {   
        flag = 1;
        while (flag)
        {
            flag = 0;  // Assume this is a success
            rx = ((double) rand()/RAND_MAX-0.5) * boxs.lenx;
            ry = ((double) rand()/RAND_MAX-0.5) * boxs.leny;
            rz = ((double) rand()/RAND_MAX-0.5) * boxs.lenz;

            // Loop through all DNa sites
            for (j = 0; j < site.dna_totl + nIonsAdded; j++)
            {   
                dx = rx - atom[j].xval;
                func_pbcs(&dx,boxs.lenx);
                dy = ry - atom[j].yval;
                func_pbcs(&dy,boxs.leny);
                dz = rz - atom[j].zval;
                func_pbcs(&dz,boxs.lenz);

                dr = sqrt(dx * dx + dy * dy + dz * dz);
                if (dr < dcut)
                {
                    flag = 1;
                    break;
                }
            }
        }

        // If no overlap, the write
        atom[i + offset].mass = Na.mass;
        atom[i + offset].xval = rx;
        atom[i + offset].yval = ry;
        atom[i + offset].zval = rz;
        sprintf(atom[i + offset].snme, "Na");
        sprintf(atom[i + offset].rnme, "SOD");
        atom[i + offset].rnum = nIonsAdded + sequ.nste * (comp + 1);
        atom[i + offset].stid = Na.stid;
        atom[i + offset].stid2 = Na.stid2;
        atom[i + offset].chrg = Na.chrg;
        nIonsAdded++;
    }


    // Adding the Magsesium
    offset += nNa;
    for (i = 0; i < nMg; i++)
    {   
        flag = 1;
        while (flag)
        {
            flag = 0;  // Assume this is a success
            rx = ((double) rand()/RAND_MAX-0.5) * boxs.lenx;
            ry = ((double) rand()/RAND_MAX-0.5) * boxs.leny;
            rz = ((double) rand()/RAND_MAX-0.5) * boxs.lenz;

            // Loop through all DNa sites
            for (j = 0; j < site.dna_totl + nIonsAdded; j++)
            {   
                dx = rx - atom[j].xval;
                func_pbcs(&dx,boxs.lenx);
                dy = ry - atom[j].yval;
                func_pbcs(&dy,boxs.leny);
                dz = rz - atom[j].zval;
                func_pbcs(&dz,boxs.lenz);

                dr = sqrt(dx * dx + dy * dy + dz * dz);
                if (dr < dcut)
                {
                    flag = 1;
                    break;
                }
            }
        }

        // If no overlap, then write
        atom[i + offset].mass = Mg.mass;
        atom[i + offset].xval = rx;
        atom[i + offset].yval = ry;
        atom[i + offset].zval = rz;
        sprintf(atom[i + offset].snme, "Mg");
        sprintf(atom[i + offset].rnme, "MAG");
        atom[i + offset].rnum = nIonsAdded + sequ.nste * (comp + 1);
        atom[i + offset].stid = Mg.stid;
        atom[i + offset].stid2 =Mg.stid2;
        atom[i + offset].chrg = Mg.chrg;
        nIonsAdded++;
    }

    // Adding the Chlorine
    offset += nMg;
    for (i = 0; i < nCl; i++)
    {   
        flag = 1;
        while (flag)
        {
            flag = 0;  // Assume this is a success
            rx = ((double) rand()/RAND_MAX-0.5) * boxs.lenx;
            ry = ((double) rand()/RAND_MAX-0.5) * boxs.leny;
            rz = ((double) rand()/RAND_MAX-0.5) * boxs.lenz;

            // Loop through all DNa sites
            for (j = 0; j < site.dna_totl + nIonsAdded; j++)
            {   
                dx = rx - atom[j].xval;
                func_pbcs(&dx,boxs.lenx);
                dy = ry - atom[j].yval;
                func_pbcs(&dy,boxs.leny);
                dz = rz - atom[j].zval;
                func_pbcs(&dz,boxs.lenz);

                dr = sqrt(dx * dx + dy * dy + dz * dz);
                if (dr < dcut)
                {
                    flag = 1;
                    break;
                }
            }
        }

        // If no overlap, the write
        atom[i + offset].mass = Cl.mass;
        atom[i + offset].xval = rx;
        atom[i + offset].yval = ry;
        atom[i + offset].zval = rz;
        sprintf(atom[i + offset].snme, "Cl");
        sprintf(atom[i + offset].rnme, "CHL");
        atom[i + offset].rnum = nIonsAdded + sequ.nste * (comp + 1);
        atom[i + offset].stid = Cl.stid;
        atom[i + offset].stid2 = Cl.stid2;
        atom[i + offset].chrg = Cl.chrg;
        nIonsAdded++;
    }
}
