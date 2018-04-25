/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Daniel Hinckley (Wisconsin/UChicago) dhinckley@wisc.edu
   ------------------------------------------------------------------------- */

#include "math.h"
#include "atom.h"
#include "stdio.h"
#include "stdlib.h"
#include "math_const.h"
#include "base_pair.h"
#include "memory.h"
#include "error.h"
#include "string.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SIX 6
#define THREE 3
#define SMALL 10E-6
#define EMPTY -12345

/* ---------------------------------------------------------------------- */

// My constructor
BasePair::BasePair(LAMMPS *lmp) : Pointers(lmp){}

/* ---------------------------------------------------------------------- */

// My destructor
BasePair::~BasePair()
{
}

/* ---------------------------------------------------------------------- */

void BasePair::assign(int *indices, int *site_type, double ***param, double ***angle,double **x)
{
    double athe[6], cnum;
    int debug = 0;
    int i,j;
    int typa,typb,typc,typd,type, typf;

    // Schematic of layout for all sites

    /* A drawing of the system
    a====b ---- d====c   --> The Hydrogen Bond
          \    /
       Q1  \  /  Q2
            \/  Cross stacking interactions
            /\ 
           /  \
          /    \
         f      e
    */

    // Storing the indices in a local array
    stea = indices[0];
    steb = indices[1];
    stec = indices[2];
    sted = indices[3];
    stee = indices[4];
    stef = indices[5];

    typa = site_type[0];
    typb = site_type[1];
    typc = site_type[2];
    typd = site_type[3];
    type = site_type[4];
    typf = site_type[5];

    // Get relevant angles, sigm, and epsilon
    if ((typb != EMPTY) && (typd != EMPTY))
      ateq[0] = angle[0][typb][typd];  // Angle between ba and dc
    if ((typb != EMPTY) && (type != EMPTY))
      ateq[1] = angle[1][typb][type]; // Angle between ba and be
    if ((typb != EMPTY) && (typf != EMPTY))
      ateq[2] = angle[2][typd][typf]; // Angle between dc and df
    if ((typb != EMPTY) && (typd != EMPTY))
      ateq[3] = angle[3][typb][typd]; // Torsion angle created by a-b-d-c
    if ((typb != EMPTY) && (typd != EMPTY))
      ateq[4] = angle[4][typb][typd]; // Angle between  a-b-d
    if ((typb != EMPTY) && (typd != EMPTY))
      ateq[5] = angle[5][typb][typd]; // Angle between  b-d-c
    for (i = 0; i < 5; i++)
    {
        ////printf("angle %d = %lf ",i+1,ateq[i]);
    }

    // Checking to make sure there aren't any issues
    for (i = 0; i < SIX; i++)
    {
       // if (ateq[i] == EMPTY) error->all(FLERR,"EMPTY value for equilibrium angle %d",i); // Throw error
    }

    // Storing the sigmas that we need
    j = 0;  // First cross-stacking interaction
    if ((typb != EMPTY) && (type != EMPTY)){
      sigm[j]  = param[j][typb][type];
      alpha[j] = param[j+3][typb][type];
      range[j] = param[j+6][typb][type];
      epsi[j] = param[j+9][typb][type];
    }

    j++;    // Second cross-stacking interaction
    if ((typd != EMPTY) && (typf != EMPTY)){
      sigm[j]  = param[j][typd][typf];
      alpha[j] = param[j+3][typd][typf];
      range[j] = param[j+6][typd][typf];
      epsi[j] = param[j+9][typd][typf];
    }

    j++;    // Base pairing interaction
    if ((typb != EMPTY) && (typd != EMPTY)){
      sigm[j]  = param[j][typb][typd];
      alpha[j] = param[j+3][typb][typd];
      range[j] = param[j+6][typb][typd];
      epsi[j] = param[j+9][typb][typd];
    }

    // Checking to make sure there aren't any issues with the values that I have stored
    for (i = 0; i < THREE; i++)
    {
        //if (sigm[i] == EMPTY) error->all(FLERR,"EMPTY value for equilibrium sigma %d",i);
        //if (alpha[i] == EMPTY) error->all(FLERR,"EMPTY value for alpha value %d",i);
        //if (range[i] == EMPTY) error->all(FLERR,"EMPTY value for range value %d",i);
        //if (epsi[i] == EMPTY) error->all(FLERR,"EMPTY value for epsilon %d",i);
    }
    
    // These ar

    // Compute distances for all interactions and store
    // if steX = -1 then arbitratily set to values 1. Value shouldn't be used...
    if ((typa != EMPTY) && (typb != EMPTY)){
      dbax = x[steb][0] - x[stea][0];
      dbay = x[steb][1] - x[stea][1];
      dbaz = x[steb][2] - x[stea][2];
    }
    else dbax = dbay = dbaz = 1;
    domain->minimum_image(dbax,dbay,dbaz);
    
    if ((typd != EMPTY) && (typc != EMPTY)){
      ddcx = x[sted][0] - x[stec][0];
      ddcy = x[sted][1] - x[stec][1];
      ddcz = x[sted][2] - x[stec][2];
    }
    else ddcx = ddcy = ddcz = 1;
    domain->minimum_image(ddcx,ddcy,ddcz);

    if ((typd != EMPTY) && (typb != EMPTY)){
      ddbx = x[sted][0] - x[steb][0];
      ddby = x[sted][1] - x[steb][1];
      ddbz = x[sted][2] - x[steb][2];
    }
    else ddbx = ddby = ddbz = 1;
    domain->minimum_image(ddbx,ddby,ddbz);

    if ((typb != EMPTY) && (type != EMPTY)){
      dbex = x[steb][0] - x[stee][0];
      dbey = x[steb][1] - x[stee][1];
      dbez = x[steb][2] - x[stee][2];
    }
    else dbex = dbey = dbez = 1;
    domain->minimum_image(dbex,dbey,dbez);

    if ((typd != EMPTY) && (typf != EMPTY)){
      ddfx = x[sted][0] - x[stef][0];
      ddfy = x[sted][1] - x[stef][1];
      ddfz = x[sted][2] - x[stef][2];
    }
    else ddfx = ddfy = ddfz = 1;
    domain->minimum_image(ddfx,ddfy,ddfz);

    // Constructing the relevant unit vectors
    dbai = 1.0 / sqrt(dbax * dbax + dbay * dbay + dbaz * dbaz);
    ebax = dbax * dbai;
    ebay = dbay * dbai;
    ebaz = dbaz * dbai;
    ddci = 1.0 / sqrt(ddcx * ddcx + ddcy * ddcy + ddcz * ddcz);
    edcx = ddcx * ddci;
    edcy = ddcy * ddci;
    edcz = ddcz * ddci;
    ddbi = 1.0 / sqrt(ddbx * ddbx + ddby * ddby + ddbz * ddbz);
    edbx = ddbx * ddbi;
    edby = ddby * ddbi;
    edbz = ddbz * ddbi;
    dbei = 1.0 / sqrt(dbex * dbex + dbey * dbey + dbez * dbez);
    ebex = dbex * dbei;
    ebey = dbey * dbei;
    ebez = dbez * dbei;
    ddfi = 1.0 / sqrt(ddfx * ddfx + ddfy * ddfy + ddfz * ddfz);
    edfx = ddfx * ddfi;
    edfy = ddfy * ddfi;
    edfz = ddfz * ddfi;

    // Store the length of each of these vectors in the class
    dba = 1.0/dbai;
    ddc = 1.0/ddci;
    ddb = 1.0/ddbi;
    dbe = 1.0/dbei;
    ddf = 1.0/ddfi;
    ddbsi = ddbi * ddbi;
    dbesi = dbei * dbei;
    ddfsi = ddfi * ddfi;

    // Calculating the relevant angles
    cphi[0] = ebax * edcx + ebay * edcy + ebaz * edcz;
    cphi[1] = ebax * ebex + ebay * ebey + ebaz * ebez;
    cphi[2] = edcx * edfx + edcy * edfy + edcz * edfz;

    padx = ebay * edbz - ebaz * edby;
    pady = ebaz * edbx - ebax * edbz;
    padz = ebax * edby - ebay * edbx;
    cosb = -(ebax * edbx + ebay * edby + ebaz * edbz);  
    if (cosb > (1. - SMALL)) cosb = (1. - SMALL);
    if (cosb < (-1. + SMALL)) cosb = (-1. + SMALL);
    isb2 = 1.0 / (1.0 - cosb * cosb);
    isnb = sqrt(isb2);

    pbcx = edcz * edby - edcy * edbz;
    pbcy = edcx * edbz - edcz * edbx;
    pbcz = edcy * edbx - edcx * edby;
    cosd = (edcx * edbx + edcy * edby + edcz * edbz); 
    if (cosd > (1. - SMALL)) cosd = (1. - SMALL);
    if (cosd < (-1. + SMALL)) cosd = (-1. + SMALL);
    isd2 = 1.0 / (1.0 - cosd * cosd);
    isnd = sqrt(isd2);

    pacx = (padz * pbcy - pady * pbcz);
    pacy = (padx * pbcz - padz * pbcx);
    pacz = (pady * pbcx - padx * pbcy);
    cnum = (pacx * edbx + pacy * edby + pacz * edbz);

    cphi[3] =
        -(padx * pbcx + pady * pbcy + padz * pbcz) * isnb * isnd;
    cphi[4] = ebax * -edbx + ebay * -edby + ebaz * -edbz;
    cphi[5] = edbx * edcx + edby * edcy + edbz * edcz;

    // Making sure that we avoid anamolous behavior
    for (i = 0; i < SIX; i++)
    {
        if (cphi[i] > (1. - SMALL)) cphi[i] = (1. - SMALL);
        if (cphi[i] < (-1. + SMALL)) cphi[i] = (-1. + SMALL);
        athe[i] = acos(cphi[i]);
    }

    // Residual code... not sure if necessary
    if (cnum < 0.0)
    {
        athe[3] = -athe[3];
    }

    for (i = 0; i < SIX; i++)
    {
        dtha[i] = athe[i] - ateq[i];
    }

    // Finished initializing all values need for interactions
}

/* ---------------------------------------------------------------------- */

double BasePair::cross_stacking(int c,  double**f)
{
    // Calculate the possible cross stacking interactions

    /*
    c == 0  cstk1
    c == 1  cstk2
    */

    int ibp;
    double cos2, sin2, ctrm, bptrm, pref, pref2, ecstk;
    double frce, engy;
    double fr1x, fr1y, fr1z, fr2x, fr2y, fr2z, fr3x, fr3y, fr3z,
        fr4x, fr4y, fr4z, fr5x, fr5y, fr5z;

    fr1x = fr1y = fr1z = 0.0;
    fr2x = fr2y = fr2z = 0.0;
    fr3x = fr3y = fr3z = 0.0;
    fr4x = fr4y = fr4z = 0.0;
    fr5x = fr5y = fr5z = 0.0;
    ibp = 2; 

    ecstk = 0.0;

    // Calculate the a - b - e cross-stacking interaction
    if (!c)
    {
        if ( (dtha[0] >= -MY_PI/(range[ibp]*2.0)) && (dtha[0] <= MY_PI/(range[ibp]*2.0)))
        {
            if ( (dtha[c+1] >= -MY_PI/(range[c]*2.0)) && (dtha[c+1] <= MY_PI/(range[c]*2.0)))
            {
                mors_norp(dbesi, alpha[c], epsi[c], sigm[c], &frce, &engy);
                ecstk += engy;

                f[stee][0] -= frce * dbex;
                f[stee][1] -= frce * dbey;
                f[stee][2] -= frce * dbez;

                f[steb][0] += frce * dbex;
                f[steb][1] += frce * dbey;
                f[steb][2] += frce * dbez;

            }
            else if (((dtha[c+1] >= MY_PI/(range[c]*2.0)) && (dtha[c+1] <= MY_PI/range[c])) || ((dtha[c+1] <= -MY_PI/(range[c]*2.0)) && (dtha[c+1] >= -MY_PI/range[c])))
            {
                cos2 = cos(range[c]*dtha[c+1]);
                sin2 = sin(range[c]*dtha[c+1]);
                ctrm = 1.0 - cos2 * cos2;

                mors_norp(dbesi, alpha[c], epsi[c], sigm[c], &frce, &engy);
                ecstk += ctrm * engy;
                pref2 = 2.0 * range[c] * cos2 * sin2 * 1.0 /sqrt(1.0-cphi[c+1]*cphi[c+1]);
                // 1 site A
                fr1x = pref2 * (dbai * (cphi[c+1] * ebax - ebex)) * engy;
                fr1y = pref2 * (dbai * (cphi[c+1] * ebay - ebey)) * engy;
                fr1z = pref2 * (dbai * (cphi[c+1] * ebaz - ebez)) * engy;
                // 2 site B
                fr2x = pref2 * (dbai * (ebex - ebax * cphi[c+1] ) + dbei * (ebax - cphi[c+1] * ebex)) * engy + ctrm * dbex * frce;
                fr2y = pref2 * (dbai * (ebey - ebay * cphi[c+1] ) + dbei * (ebay - cphi[c+1] * ebey)) * engy + ctrm * dbey * frce;
                fr2z = pref2 * (dbai * (ebez - ebaz * cphi[c+1] ) + dbei * (ebaz - cphi[c+1] * ebez)) * engy + ctrm * dbez * frce;
                // 3 site E      
                fr3x = pref2 * (dbei * (ebex * cphi[c+1] - ebax)) * engy - ctrm * dbex * frce;
                fr3y = pref2 * (dbei * (ebey * cphi[c+1] - ebay)) * engy - ctrm * dbey * frce;
                fr3z = pref2 * (dbei * (ebez * cphi[c+1] - ebaz)) * engy - ctrm * dbez * frce;

                f[stea][0] += fr1x;
                f[stea][1] += fr1y;
                f[stea][2] += fr1z;

                f[steb][0] += fr2x;
                f[steb][1] += fr2y;
                f[steb][2] += fr2z;

                f[stee][0] += fr3x;
                f[stee][1] += fr3y;
                f[stee][2] += fr3z;
            }
            else if ((dtha[c+1] > MY_PI/range[c]) ||  (dtha[c+1] < -MY_PI/range[c])) 
            {
            }
        }
        else if (((dtha[0] >= MY_PI/(range[ibp]*2.0)) && (dtha[0] <= MY_PI/range[ibp])) || ((dtha[0] <= -MY_PI/(range[ibp]*2.0)) && (dtha[0] >= -MY_PI/range[ibp])))
        {
            cos2 = cos(range[ibp]*dtha[0]);
            sin2 = sin(range[ibp]*dtha[0]);
            bptrm = 1.0 - cos2 * cos2;
            pref = 2.0 * range[ibp] * cos2 * sin2 * 1.0 /sqrt(1.0-cphi[0]*cphi[0]);
            
            if ( (dtha[c+1] >= -MY_PI/(range[c]*2.0)) && (dtha[c+1] <= MY_PI/(range[c]*2.0)))
            {
                mors_norp(dbesi, alpha[c], epsi[c], sigm[c], &frce, &engy);
                ecstk += bptrm * engy;

                // Site A
                fr1x = pref * dbai * (ebax * cphi[0] - edcx) * engy;
                fr1y = pref * dbai * (ebay * cphi[0] - edcy) * engy;
                fr1z = pref * dbai * (ebaz * cphi[0] - edcz) * engy;
                // Site B
                fr2x = pref * dbai * (edcx - ebax * cphi[0]) * engy + bptrm * dbex * frce;
                fr2y = pref * dbai * (edcy - ebay * cphi[0]) * engy + bptrm * dbey * frce;
                fr2z = pref * dbai * (edcz - ebaz * cphi[0]) * engy + bptrm * dbez * frce;
                // Site C
                fr3x = pref * ddci * (edcx * cphi[0] - ebax) * engy;
                fr3y = pref * ddci * (edcy * cphi[0] - ebay) * engy;
                fr3z = pref * ddci * (edcz * cphi[0] - ebaz) * engy;
                // Site D
                fr4x = pref * ddci * (ebax - edcx * cphi[0]) * engy;
                fr4y = pref * ddci * (ebay - edcy * cphi[0]) * engy; 
                fr4z = pref * ddci * (ebaz - edcz * cphi[0]) * engy; 
                // Site E
                fr5x = -bptrm * dbex * frce;
                fr5y = -bptrm * dbey * frce;
                fr5z = -bptrm * dbez * frce;

                f[stea][0] += fr1x;
                f[stea][1] += fr1y;
                f[stea][2] += fr1z;

                f[steb][0] += fr2x;
                f[steb][1] += fr2y;
                f[steb][2] += fr2z;

                f[stec][0] += fr3x;
                f[stec][1] += fr3y;
                f[stec][2] += fr3z;

                f[sted][0] += fr4x;
                f[sted][1] += fr4y;
                f[sted][2] += fr4z;

                f[stee][0] += fr5x;
                f[stee][1] += fr5y;
                f[stee][2] += fr5z;
            }
            else if (((dtha[c+1] >= MY_PI/(range[c]*2.0)) && (dtha[c+1] <= MY_PI/range[c])) || ((dtha[c+1] <= -MY_PI/(range[c]*2.0)) && (dtha[c+1] >= -MY_PI/range[c])))
            {
                cos2 = cos(range[c]*dtha[c+1]);
                sin2 = sin(range[c]*dtha[c+1]);
                ctrm = 1.0 - cos2 * cos2;
                pref2 = 2.0 * range[c] * cos2 * sin2 * 1.0 /sqrt(1.0-cphi[c+1]*cphi[c+1]);

                mors_norp(dbesi, alpha[c], epsi[c], sigm[c], &frce, &engy);
                ecstk += bptrm * ctrm * engy;

                // The potential is modulated by the angle for the hydrogen bonding

                // Site A
                fr1x = pref * dbai * (ebax * cphi[0] - edcx) * ctrm * engy + pref2 * (dbai * (ebax * cphi[c+1] - ebex)) * bptrm * engy;
                fr1y = pref * dbai * (ebay * cphi[0] - edcy) * ctrm * engy + pref2 * (dbai * (ebay * cphi[c+1] - ebey)) * bptrm * engy;
                fr1z = pref * dbai * (ebaz * cphi[0] - edcz) * ctrm * engy + pref2 * (dbai * (ebaz * cphi[c+1] - ebez)) * bptrm * engy;
                // Site B
                fr2x = pref * dbai * (edcx - ebax * cphi[0]) * ctrm * engy + pref2 * (dbai * (ebex - cphi[c+1] * ebax) + dbei * (ebax - cphi[c+1] * ebex)) * bptrm * engy + bptrm * ctrm * dbex * frce;
                fr2y = pref * dbai * (edcy - ebay * cphi[0]) * ctrm * engy + pref2 * (dbai * (ebey - cphi[c+1] * ebay) + dbei * (ebay - cphi[c+1] * ebey)) * bptrm * engy + bptrm * ctrm * dbey * frce;
                fr2z = pref * dbai * (edcz - ebaz * cphi[0]) * ctrm * engy + pref2 * (dbai * (ebez - cphi[c+1] * ebaz) + dbei * (ebaz - cphi[c+1] * ebez)) * bptrm * engy + bptrm * ctrm * dbez * frce;
                // Site C
                fr3x = pref * ddci * (edcx * cphi[0] - ebax) * ctrm * engy;
                fr3y = pref * ddci * (edcy * cphi[0] - ebay) * ctrm * engy;
                fr3z = pref * ddci * (edcz * cphi[0] - ebaz) * ctrm * engy;
                // Site D
                fr4x = pref * ddci * (ebax - edcx * cphi[0]) * ctrm * engy;
                fr4y = pref * ddci * (ebay - edcy * cphi[0]) * ctrm * engy; 
                fr4z = pref * ddci * (ebaz - edcz * cphi[0]) * ctrm * engy; 
                // Site E
                fr5x = pref2 * (dbei * (ebex * cphi[c+1] - ebax)) * bptrm * engy - bptrm * ctrm * dbex * frce ;
                fr5y = pref2 * (dbei * (ebey * cphi[c+1] - ebay)) * bptrm * engy - bptrm * ctrm * dbey * frce ;
                fr5z = pref2 * (dbei * (ebez * cphi[c+1] - ebaz)) * bptrm * engy - bptrm * ctrm * dbez * frce ;

                f[stea][0] += fr1x;
                f[stea][1] += fr1y;
                f[stea][2] += fr1z;

                f[steb][0] += fr2x;
                f[steb][1] += fr2y;
                f[steb][2] += fr2z;

                f[stec][0] += fr3x;
                f[stec][1] += fr3y;
                f[stec][2] += fr3z;

                f[sted][0] += fr4x;
                f[sted][1] += fr4y;
                f[sted][2] += fr4z;

                f[stee][0] += fr5x;
                f[stee][1] += fr5y;
                f[stee][2] += fr5z;
            }
            else if ((dtha[c+1] > MY_PI/range[c])  ||  (dtha[c+1] < -MY_PI/range[c]))
            {
            }
        }
    }

    // Calculate the other cross stacking interction
    if (c == 1)
    {
        if ( (dtha[0] >= -MY_PI/(range[ibp]*2.0)) && (dtha[0] <= MY_PI/(range[ibp]*2.0)))
        {
            // Now examining the other cross stacking interaction
            if ( (dtha[c+1] >= -MY_PI/(range[c]*2.0)) && (dtha[c+1] <= MY_PI/(range[c]*2.0)))
            {
                mors_norp(ddfsi, alpha[c], epsi[c], sigm[c], &frce, &engy);
                ecstk += engy;

                // Apply the force
                f[stef][0] -= frce * ddfx;
                f[stef][1] -= frce * ddfy;
                f[stef][2] -= frce * ddfz;

                f[sted][0] += frce * ddfx;
                f[sted][1] += frce * ddfy;
                f[sted][2] += frce * ddfz;
            }
            else if (((dtha[c+1] >= MY_PI/(range[c]*2.0)) && (dtha[c+1] <= MY_PI/range[c])) || ((dtha[c+1] <= -MY_PI/(range[c]*2.0)) && (dtha[c+1] >= -MY_PI/range[c])))
            {
                cos2 = cos(range[c]*dtha[c+1]);
                sin2 = sin(range[c]*dtha[c+1]);
                ctrm = 1.0 - cos2 * cos2;

                mors_norp(ddfsi, alpha[c], epsi[c], sigm[c], &frce, &engy);
                ecstk += ctrm * engy;
                pref2 = 2.0 * range[c] * cos2 * sin2 * 1.0 /sqrt(1.0-cphi[c+1]*cphi[c+1]);
                // 1 site C
                fr1x = pref2 * (ddci * (cphi[c+1] * edcx - edfx)) * engy;
                fr1y = pref2 * (ddci * (cphi[c+1] * edcy - edfy)) * engy;
                fr1z = pref2 * (ddci * (cphi[c+1] * edcz - edfz)) * engy;
                // 2 site D
                fr2x = pref2 * (ddci * (edfx - cphi[c+1] * edcx) + ddfi * (edcx - cphi[c+1] * edfx)) * engy + ctrm * ddfx * frce;
                fr2y = pref2 * (ddci * (edfy - cphi[c+1] * edcy) + ddfi * (edcy - cphi[c+1] * edfy)) * engy + ctrm * ddfy * frce;
                fr2z = pref2 * (ddci * (edfz - cphi[c+1] * edcz) + ddfi * (edcz - cphi[c+1] * edfz)) * engy + ctrm * ddfz * frce;
                // 3 site F      
                fr3x = pref2 * (ddfi * (cphi[c+1] * edfx - edcx)) * engy - ctrm * ddfx * frce;
                fr3y = pref2 * (ddfi * (cphi[c+1] * edfy - edcy)) * engy - ctrm * ddfy * frce;
                fr3z = pref2 * (ddfi * (cphi[c+1] * edfz - edcz)) * engy - ctrm * ddfz * frce;

                f[stec][0] += fr1x;
                f[stec][1] += fr1y;
                f[stec][2] += fr1z;

                f[sted][0] += fr2x;
                f[sted][1] += fr2y;
                f[sted][2] += fr2z;

                f[stef][0] += fr3x;
                f[stef][1] += fr3y;
                f[stef][2] += fr3z;
            }
            else if ((dtha[c+1] > MY_PI/range[c]) || (dtha[c+1] < -MY_PI/range[c])) 
            {
            }
        }
        else if (((dtha[0] >= MY_PI/(range[ibp]*2.0)) && (dtha[0] <= MY_PI/range[ibp])) || ((dtha[0] <= -MY_PI/(range[ibp]*2.0)) && (dtha[0] >= -MY_PI/range[ibp])))
        {
            cos2 = cos(range[ibp]*dtha[0]);
            sin2 = sin(range[ibp]*dtha[0]);
            bptrm = 1.0 - cos2 * cos2;
            pref = 2.0 * range[ibp] * cos2 * sin2 * 1.0 /sqrt(1.0-cphi[0]*cphi[0]);

            if ( (dtha[c+1] >= -MY_PI/(range[c]*2.0)) && (dtha[c+1] <= MY_PI/(range[c]*2.0)))
            {
                mors_norp(ddfsi, alpha[c], epsi[c], sigm[c], &frce, &engy);
                ecstk += bptrm * engy;

                // The potential is modulated by the angle for the hydrogen bonding
                // Site A
                fr1x = pref * dbai * (ebax * cphi[0] - edcx) * engy;
                fr1y = pref * dbai * (ebay * cphi[0] - edcy) * engy;
                fr1z = pref * dbai * (ebaz * cphi[0] - edcz) * engy;
                // Site B
                fr2x = pref * dbai * (edcx - ebax * cphi[0]) * engy;
                fr2y = pref * dbai * (edcy - ebay * cphi[0]) * engy;
                fr2z = pref * dbai * (edcz - ebaz * cphi[0]) * engy;
                // Site C
                fr3x = pref * ddci * (edcx * cphi[0] - ebax) * engy;
                fr3y = pref * ddci * (edcy * cphi[0] - ebay) * engy;
                fr3z = pref * ddci * (edcz * cphi[0] - ebaz) * engy;
                // Site D
                fr4x = pref * ddci * (ebax - edcx * cphi[0]) * engy + bptrm * ddfx * frce;
                fr4y = pref * ddci * (ebay - edcy * cphi[0]) * engy + bptrm * ddfy * frce; 
                fr4z = pref * ddci * (ebaz - edcz * cphi[0]) * engy + bptrm * ddfz * frce; 
                // Site F
                fr5x = -bptrm * ddfx * frce;
                fr5y = -bptrm * ddfy * frce;
                fr5z = -bptrm * ddfz * frce;

                f[stea][0] += fr1x;
                f[stea][1] += fr1y;
                f[stea][2] += fr1z;

                f[steb][0] += fr2x;
                f[steb][1] += fr2y;
                f[steb][2] += fr2z;

                f[stec][0] += fr3x;
                f[stec][1] += fr3y;
                f[stec][2] += fr3z;

                f[sted][0] += fr4x;
                f[sted][1] += fr4y;
                f[sted][2] += fr4z;

                f[stef][0] += fr5x;
                f[stef][1] += fr5y;
                f[stef][2] += fr5z;
            }
            else if (((dtha[c+1] >= MY_PI/(range[c]*2.0)) && (dtha[c+1] <= MY_PI/range[c])) || ((dtha[c+1] <= -MY_PI/(range[c]*2.0)) && (dtha[c+1] >= -MY_PI/range[c])))
            {
                cos2 = cos(range[c]*dtha[c+1]);
                sin2 = sin(range[c]*dtha[c+1]);
                ctrm = 1.0 - cos2 * cos2;
                pref2 = 2.0 * range[c] * cos2 * sin2 * 1.0 /sqrt(1.0-cphi[c+1]*cphi[c+1]);

                mors_norp(ddfsi, alpha[c], epsi[c], sigm[c], &frce, &engy);
                ecstk += bptrm * ctrm * engy;

                // The potential is modulated by the angle for the hydrogen bonding
                // Site A
                fr1x = pref * dbai * (ebax * cphi[0] - edcx) * ctrm * engy;
                fr1y = pref * dbai * (ebay * cphi[0] - edcy) * ctrm * engy;
                fr1z = pref * dbai * (ebaz * cphi[0] - edcz) * ctrm * engy;
                // Site B
                fr2x = pref * dbai * (edcx - ebax * cphi[0]) * ctrm * engy;
                fr2y = pref * dbai * (edcy - ebay * cphi[0]) * ctrm * engy;
                fr2z = pref * dbai * (edcz - ebaz * cphi[0]) * ctrm * engy;
                // Site C
                fr3x = pref * ddci * (edcx * cphi[0] - ebax) * ctrm * engy + pref2 * (ddci * (edcx * cphi[c+1] - edfx)) * bptrm * engy;
                fr3y = pref * ddci * (edcy * cphi[0] - ebay) * ctrm * engy + pref2 * (ddci * (edcy * cphi[c+1] - edfy)) * bptrm * engy;
                fr3z = pref * ddci * (edcz * cphi[0] - ebaz) * ctrm * engy + pref2 * (ddci * (edcz * cphi[c+1] - edfz)) * bptrm * engy;
                // Site D 
                fr4x = pref * ddci * (ebax - edcx * cphi[0]) * ctrm * engy + pref2 * (ddci * (edfx - cphi[c+1] * edcx) + ddfi * (edcx - cphi[c+1] * edfx)) * bptrm * engy + bptrm * ctrm * ddfx * frce;
                fr4y = pref * ddci * (ebay - edcy * cphi[0]) * ctrm * engy + pref2 * (ddci * (edfy - cphi[c+1] * edcy) + ddfi * (edcy - cphi[c+1] * edfy)) * bptrm * engy + bptrm * ctrm * ddfy * frce; 
                fr4z = pref * ddci * (ebaz - edcz * cphi[0]) * ctrm * engy + pref2 * (ddci * (edfz - cphi[c+1] * edcz) + ddfi * (edcz - cphi[c+1] * edfz)) * bptrm * engy + bptrm * ctrm * ddfz * frce; 
                // Site F
                fr5x = pref2 * (ddfi * (edfx * cphi[c+1] - edcx)) * bptrm * engy - bptrm * ctrm * ddfx * frce ;
                fr5y = pref2 * (ddfi * (edfy * cphi[c+1] - edcy)) * bptrm * engy - bptrm * ctrm * ddfy * frce ;
                fr5z = pref2 * (ddfi * (edfz * cphi[c+1] - edcz)) * bptrm * engy - bptrm * ctrm * ddfz * frce ;

                f[stea][0] += fr1x;
                f[stea][1] += fr1y;
                f[stea][2] += fr1z;

                f[steb][0] += fr2x;
                f[steb][1] += fr2y;
                f[steb][2] += fr2z;

                f[stec][0] += fr3x;
                f[stec][1] += fr3y;
                f[stec][2] += fr3z;

                f[sted][0] += fr4x;
                f[sted][1] += fr4y;
                f[sted][2] += fr4z;

                f[stef][0] += fr5x;
                f[stef][1] += fr5y;
                f[stef][2] += fr5z;
            }
            else if ((dtha[c+1] > MY_PI/range[c]) || (dtha[c+1] <= -MY_PI/range[c]))
            {
            }
        }
    }
    return ecstk;
}

/* ---------------------------------------------------------------------- */

double BasePair::base_pairing(double **f)
{
    double phi_factor, ftor, ebasepair, first_term, second_term,
    cosine, sine, cosine2, sine2, hbon_cosine_term, hbon_cosine_term2;
    double prefactor, prefactor2;

    double fra1, frb1, frb2, frc1, frc2, frd1;
    double fr1x, fr1y, fr1z, fr2x, fr2y, fr2z, fr3x, fr3y, fr3z,
        fr4x, fr4y, fr4z, fr5x, fr5y, fr5z, fr6x, fr6y, fr6z;
        double frce, engy;
    int iphi = 3, itha1 = 4, itha2 = 5, ibp = 2;
    ebasepair = 0.0;
    phi_factor = 0.5 * (1+cos(dtha[iphi]));
    ftor = 0.5 * sin(dtha[iphi]);

    // The repulsion is always present
    if (ddb < sigm[ibp])
    {
        mors_rp(ddbsi, alpha[ibp], epsi[ibp], sigm[ibp], &frce, &engy);
        ebasepair += engy;

        f[steb][0] -= frce * ddbx;
        f[steb][1] -= frce * ddby;
        f[steb][2] -= frce * ddbz;

        f[sted][0] += frce * ddbx;
        f[sted][1] += frce * ddby;
        f[sted][2] += frce * ddbz;
    }



    if ( (dtha[itha1] >= -MY_PI/(range[ibp]*2.0)) && (dtha[itha1] <= MY_PI/(range[ibp]*2.0)))
    {
        if ( (dtha[itha2] >= -MY_PI/(range[ibp]*2.0)) && (dtha[itha2] <= MY_PI/(range[ibp]*2.0)))
        {
            // Modulated attractive term
            mors_norp(ddbsi, alpha[ibp], epsi[ibp], sigm[ibp], &frce, &engy);
            ebasepair += phi_factor * engy;

            first_term = engy;
            second_term = phi_factor * frce;
            fra1 = -ftor * dbai * isb2 * first_term;
            f[stea][0] += fra1 * padx;
            f[stea][1] += fra1 * pady;
            f[stea][2] += fra1 * padz;

            frb1 = ftor * (ddb - dba * cosb) * dbai * ddbi * isb2 * first_term;
            frb2 = ftor * cosd * ddbi * isd2 * first_term;

            f[steb][0] += (frb1 * padx + frb2 * pbcx) - second_term * ddbx;
            f[steb][1] += (frb1 * pady + frb2 * pbcy) - second_term * ddby;
            f[steb][2] += (frb1 * padz + frb2 * pbcz) - second_term * ddbz;

            frc1 = ftor * (ddb - ddc * cosd) * ddci * ddbi * isd2 * first_term;
            frc2 = ftor * cosb * ddbi * isb2 * first_term;
            f[sted][0] += (frc1 * pbcx + frc2 * padx) + second_term * ddbx;
            f[sted][1] += (frc1 * pbcy + frc2 * pady) + second_term * ddby;
            f[sted][2] += (frc1 * pbcz + frc2 * padz) + second_term * ddbz;

            frd1 = -ftor * ddci * isd2 * first_term;
            f[stec][0] += frd1 * pbcx;
            f[stec][1] += frd1 * pbcy;
            f[stec][2] += frd1 * pbcz;
        }
        else if (((dtha[itha2] >= MY_PI/(range[ibp]*2.0)) && (dtha[itha2] <= MY_PI/range[ibp])) || ((dtha[itha2] <= -MY_PI/(range[ibp]*2.0)) && (dtha[itha2] >= -MY_PI/range[ibp])))
        {
            // Now I modulated using thata 2 
            cosine2 = cos(range[ibp]*dtha[itha2]);
            sine2 = sin(range[ibp]*dtha[itha2]);
            hbon_cosine_term2 = 1.0 - cosine2 * cosine2;

            mors_norp(ddbsi, alpha[ibp], epsi[ibp], sigm[ibp], &frce, &engy);
            ebasepair += phi_factor * hbon_cosine_term2 * engy;

            prefactor2 = 2.0 * range[ibp] * cosine2 * sine2 * 1.0 /sqrt(1.0-cphi[5]*cphi[5]);
            // 1 site B
            fr1x = prefactor2 * (ddbi * (cphi[5] * edbx - edcx)) * engy - hbon_cosine_term2 * ddbx * frce;
            fr1y = prefactor2 * (ddbi * (cphi[5] * edby - edcy)) * engy - hbon_cosine_term2 * ddby * frce;
            fr1z = prefactor2 * (ddbi * (cphi[5] * edbz - edcz)) * engy - hbon_cosine_term2 * ddbz * frce;
            // 2 site D
            fr2x = prefactor2 * (ddci * (edbx - edcx * cphi[5] ) + ddbi * (edcx - cphi[5] * edbx)) * engy + hbon_cosine_term2 * ddbx * frce;
            fr2y = prefactor2 * (ddci * (edby - edcy * cphi[5] ) + ddbi * (edcy - cphi[5] * edby)) * engy + hbon_cosine_term2 * ddby * frce;
            fr2z = prefactor2 * (ddci * (edbz - edcz * cphi[5] ) + ddbi * (edcz - cphi[5] * edbz)) * engy + hbon_cosine_term2 * ddbz * frce;
            // 3 site C     
            fr3x = prefactor2 * (ddci * (edcx * cphi[5] - edbx)) * engy;
            fr3y = prefactor2 * (ddci * (edcy * cphi[5] - edby)) * engy;
            fr3z = prefactor2 * (ddci * (edcz * cphi[5] - edbz)) * engy;

            f[steb][0] += fr1x * phi_factor;
            f[steb][1] += fr1y * phi_factor;
            f[steb][2] += fr1z * phi_factor;
            f[sted][0] += fr2x * phi_factor;
            f[sted][1] += fr2y * phi_factor;
            f[sted][2] += fr2z * phi_factor;
            f[stec][0] += fr3x * phi_factor;
            f[stec][1] += fr3y * phi_factor;
            f[stec][2] += fr3z * phi_factor;

            first_term = engy * hbon_cosine_term2;
            fra1 = -ftor * dbai * isb2 * first_term;
            f[stea][0] += fra1 * padx;
            f[stea][1] += fra1 * pady;
            f[stea][2] += fra1 * padz;

            frb1 = ftor * (ddb - dba * cosb) * dbai * ddbi * isb2 * first_term;
            frb2 = ftor * cosd * ddbi * isd2 * first_term;

            f[steb][0] += (frb1 * padx + frb2 * pbcx);
            f[steb][1] += (frb1 * pady + frb2 * pbcy);
            f[steb][2] += (frb1 * padz + frb2 * pbcz);

            frc1 = ftor * (ddb - ddc * cosd) * ddci * ddbi * isd2 * first_term;
            frc2 = ftor * cosb * ddbi * isb2 * first_term;
            f[sted][0] += (frc1 * pbcx + frc2 * padx);
            f[sted][1] += (frc1 * pbcy + frc2 * pady);
            f[sted][2] += (frc1 * pbcz + frc2 * padz);

            frd1 = -ftor * ddci * isd2 * first_term;
            f[stec][0] += frd1 * pbcx;
            f[stec][1] += frd1 * pbcy;
            f[stec][2] += frd1 * pbcz;
        }
        else
        {
                //f//printf(stdout,"PANIC!: angle difference %lf is outside of permitted range for CS1, %ld, %ld\n",dtha[itha2],stee, steb);
        }
    }
    else if (((dtha[itha1] >= MY_PI/(range[ibp]*2.0)) && (dtha[itha1] <= MY_PI/range[ibp])) || ((dtha[itha1] <= -MY_PI/(range[ibp]*2.0)) && (dtha[itha1] >= -MY_PI/range[ibp])))
    {
        cosine = cos(range[ibp]*dtha[itha1]);
        sine = sin(range[ibp]*dtha[itha1]);
        hbon_cosine_term = 1.0 - cosine * cosine;

        // Full for dtha[itha2] and modulated dtha1
        if ( (dtha[itha2] >= -MY_PI/(range[ibp]*2.0)) && (dtha[itha2] <= MY_PI/(range[ibp]*2.0)))
        {
            mors_norp(ddbsi, alpha[ibp], epsi[ibp], sigm[ibp], &frce, &engy);

            ebasepair += phi_factor * hbon_cosine_term * engy;

            prefactor = 2.0 * range[ibp] * cosine * sine * 1.0 /sqrt(1.0-cphi[4]*cphi[4]);
            // 1 site A
            fr1x = prefactor * (dbai * (cphi[4] * ebax + edbx)) * engy;
            fr1y = prefactor * (dbai * (cphi[4] * ebay + edby)) * engy;
            fr1z = prefactor * (dbai * (cphi[4] * ebaz + edbz)) * engy;
            // 2 site B
            fr2x = prefactor * (dbai * (-edbx - ebax * cphi[4] ) + ddbi * (ebax + cphi[4] * edbx)) * engy - hbon_cosine_term * ddbx * frce;
            fr2y = prefactor * (dbai * (-edby - ebay * cphi[4] ) + ddbi * (ebay + cphi[4] * edby)) * engy - hbon_cosine_term * ddby * frce;
            fr2z = prefactor * (dbai * (-edbz - ebaz * cphi[4] ) + ddbi * (ebaz + cphi[4] * edbz)) * engy - hbon_cosine_term * ddbz * frce;
            // 3 site D     
            fr3x = prefactor * (ddbi * (-edbx * cphi[4] - ebax)) * engy + hbon_cosine_term * ddbx * frce;
            fr3y = prefactor * (ddbi * (-edby * cphi[4] - ebay)) * engy + hbon_cosine_term * ddby * frce;
            fr3z = prefactor * (ddbi * (-edbz * cphi[4] - ebaz)) * engy + hbon_cosine_term * ddbz * frce;

            f[stea][0] += fr1x * phi_factor;
            f[stea][1] += fr1y * phi_factor;
            f[stea][2] += fr1z * phi_factor;
            f[steb][0] += fr2x * phi_factor;
            f[steb][1] += fr2y * phi_factor;
            f[steb][2] += fr2z * phi_factor;
            f[sted][0] += fr3x * phi_factor;
            f[sted][1] += fr3y * phi_factor;
            f[sted][2] += fr3z * phi_factor;

            first_term = engy * hbon_cosine_term;
            fra1 = -ftor * dbai * isb2 * first_term;
            f[stea][0] += fra1 * padx;
            f[stea][1] += fra1 * pady;
            f[stea][2] += fra1 * padz;

            frb1 = ftor * (ddb - dba * cosb) * dbai * ddbi * isb2 * first_term;
            frb2 = ftor * cosd * ddbi * isd2 * first_term;

            f[steb][0] += (frb1 * padx + frb2 * pbcx);
            f[steb][1] += (frb1 * pady + frb2 * pbcy);
            f[steb][2] += (frb1 * padz + frb2 * pbcz);

            frc1 = ftor * (ddb - ddc * cosd) * ddci * ddbi * isd2 * first_term;
            frc2 = ftor * cosb * ddbi * isb2 * first_term;
            f[sted][0] += (frc1 * pbcx + frc2 * padx);
            f[sted][1] += (frc1 * pbcy + frc2 * pady);
            f[sted][2] += (frc1 * pbcz + frc2 * padz);

            frd1 = -ftor * ddci * isd2 * first_term;
            f[stec][0] += frd1 * pbcx;
            f[stec][1] += frd1 * pbcy;
            f[stec][2] += frd1 * pbcz;
        }

        // Both potentials are modulated
        else if (((dtha[itha2] >= MY_PI/(range[ibp]*2.0)) && (dtha[itha2] <= MY_PI/range[ibp])) || ((dtha[itha2] <= -MY_PI/(range[ibp]*2.0)) && (dtha[itha2] >= -MY_PI/range[ibp])))
        {
            cosine2 = cos(range[ibp]*dtha[itha2]);
            sine2 = sin(range[ibp]*dtha[itha2]);
            hbon_cosine_term2 = 1.0 - cosine2 * cosine2;

            mors_norp(ddbsi, alpha[ibp], epsi[ibp], sigm[ibp], &frce, &engy);

            ebasepair += phi_factor * hbon_cosine_term * hbon_cosine_term2 * engy;
            prefactor2 = 2.0 * range[ibp] * cosine2 * sine2 * 1.0 /sqrt(1.0-cphi[5]*cphi[5]);
            prefactor = 2.0 * range[ibp] * cosine * sine * 1.0 /sqrt(1.0-cphi[4]*cphi[4]);
            // Site A 
            fr1x = prefactor * (dbai * (cphi[4] * ebax + edbx)) * hbon_cosine_term2 * engy; 
            fr1y = prefactor * (dbai * (cphi[4] * ebay + edby)) * hbon_cosine_term2 * engy; 
            fr1z = prefactor * (dbai * (cphi[4] * ebaz + edbz)) * hbon_cosine_term2 * engy; 
            // Site B
            fr2x = prefactor * (dbai * (-edbx - ebax * cphi[4] ) + ddbi * (ebax + cphi[4] * edbx)) * hbon_cosine_term2 * engy + prefactor2 * (ddbi * (cphi[5] * edbx - edcx)) * hbon_cosine_term * engy - hbon_cosine_term * hbon_cosine_term2 * ddbx * frce; 
            fr2y = prefactor * (dbai * (-edby - ebay * cphi[4] ) + ddbi * (ebay + cphi[4] * edby)) * hbon_cosine_term2 * engy + prefactor2 * (ddbi * (cphi[5] * edby - edcy)) * hbon_cosine_term * engy - hbon_cosine_term * hbon_cosine_term2 * ddby * frce; 
            fr2z = prefactor * (dbai * (-edbz - ebaz * cphi[4] ) + ddbi * (ebaz + cphi[4] * edbz)) * hbon_cosine_term2 * engy + prefactor2 * (ddbi * (cphi[5] * edbz - edcz)) * hbon_cosine_term * engy - hbon_cosine_term * hbon_cosine_term2 * ddbz * frce; 
            // Site C
            fr3x = prefactor2 * (ddci * (edcx * cphi[5] - edbx)) * hbon_cosine_term * engy; 
            fr3y = prefactor2 * (ddci * (edcy * cphi[5] - edby)) * hbon_cosine_term * engy; 
            fr3z = prefactor2 * (ddci * (edcz * cphi[5] - edbz)) * hbon_cosine_term * engy; 
            // Site D
            fr4x = prefactor * (ddbi * (-edbx * cphi[4] - ebax)) * hbon_cosine_term2 * engy + prefactor2 * (ddci * (edbx - edcx * cphi[5] ) + ddbi * (edcx - cphi[5] * edbx)) * hbon_cosine_term * engy + hbon_cosine_term * hbon_cosine_term2 * ddbx * frce; 
            fr4y = prefactor * (ddbi * (-edby * cphi[4] - ebay)) * hbon_cosine_term2 * engy + prefactor2 * (ddci * (edby - edcy * cphi[5] ) + ddbi * (edcy - cphi[5] * edby)) * hbon_cosine_term * engy + hbon_cosine_term * hbon_cosine_term2 * ddby * frce; 
            fr4z = prefactor * (ddbi * (-edbz * cphi[4] - ebaz)) * hbon_cosine_term2 * engy + prefactor2 * (ddci * (edbz - edcz * cphi[5] ) + ddbi * (edcz - cphi[5] * edbz)) * hbon_cosine_term * engy + hbon_cosine_term * hbon_cosine_term2 * ddbz * frce; 

            f[stea][0] += fr1x * phi_factor;
            f[stea][1] += fr1y * phi_factor;
            f[stea][2] += fr1z * phi_factor;
            f[steb][0] += fr2x * phi_factor;
            f[steb][1] += fr2y * phi_factor;
            f[steb][2] += fr2z * phi_factor;
            f[stec][0] += fr3x * phi_factor;
            f[stec][1] += fr3y * phi_factor;
            f[stec][2] += fr3z * phi_factor;
            f[sted][0] += fr4x * phi_factor;
            f[sted][1] += fr4y * phi_factor;
            f[sted][2] += fr4z * phi_factor;

            // Code block
            first_term = engy * hbon_cosine_term * hbon_cosine_term2;
            fra1 = -ftor * dbai * isb2 * first_term;
            f[stea][0] += fra1 * padx;
            f[stea][1] += fra1 * pady;
            f[stea][2] += fra1 * padz;

            frb1 = ftor * (ddb - dba * cosb) * dbai * ddbi * isb2 * first_term;
            frb2 = ftor * cosd * ddbi * isd2 * first_term;

            f[steb][0] += (frb1 * padx + frb2 * pbcx);
            f[steb][1] += (frb1 * pady + frb2 * pbcy);
            f[steb][2] += (frb1 * padz + frb2 * pbcz);

            frc1 = ftor * (ddb - ddc * cosd) * ddci * ddbi * isd2 * first_term;
            frc2 = ftor * cosb * ddbi * isb2 * first_term;
            f[sted][0] += (frc1 * pbcx + frc2 * padx);
            f[sted][1] += (frc1 * pbcy + frc2 * pady);
            f[sted][2] += (frc1 * pbcz + frc2 * padz);

            frd1 = -ftor * ddci * isd2 * first_term;
            f[stec][0] += frd1 * pbcx;
            f[stec][1] += frd1 * pbcy;
            f[stec][2] += frd1 * pbcz;
        }
        else
        {
        }
    }
    else if ((dtha[itha1] >= MY_PI/range[ibp]) || (dtha[itha1] <= -MY_PI/range[ibp])) 
    {
    }
    else
    {
    }

    // Now return total energy of interaction
    return ebasepair;
}

/* ---------------------------------------------------------------------- */

void BasePair::mors_norp( double drsi, double alfa, double epsi, double sigm,
    double *forc, double *ener)
{
    double disi, dist, argu;
    
    disi = sqrt(drsi);
    dist = 1.0 / disi;
    if (dist > sigm)
    {
    argu = alfa * (dist - sigm);
    *ener = epsi * (1.0 - exp(-argu)) * (1.0 - exp(-argu)) - epsi;
    *forc =
        -2.0 * alfa * epsi * disi * exp(-argu) * (1.0 - exp(-argu));
    }
    else
    {
        *ener = -epsi;
        *forc = 0.0;
    }
}

/* ---------------------------------------------------------------------- */

void BasePair::mors_rp(double drsi, double alfa, double epsi, double sigm, 
    double *forc, double *ener)
{
    double disi, dist, argu;
    
    disi = sqrt(drsi);
    dist = 1.0 / disi;
    if (dist < sigm)
    {
    argu = alfa * (dist - sigm);
    *ener = epsi * (1.0 - exp(-argu)) * (1.0 - exp(-argu));
    *forc =
        -2.0 * alfa * epsi * disi * exp(-argu) * (1.0 - exp(-argu));
    }
    else
    {
        *ener = 0.0;
        *forc = 0.0;
    }
}

