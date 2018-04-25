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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_3spn2.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "fix.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
//#include "base_pair.h"

#define _SF_ 8.908987181403393E-1
#define _EEF1_ 0.17958712212516656 // 1/pi^(3/2);used in solvation calculation
#define _KB_ 1.3806505E-23
#define _NA_ 6.0221415E23
#define _EC_ 1.60217653E-19
#define _PV_ 8.8541878176E-22
#define EMPTY -12345
#define FOUR 4
#define SIX 6

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

Pair3spn2::Pair3spn2(LAMMPS *lmp) : Pair(lmp) {
  // An array to store additional quantities for energy calculations  
  nextra = 6;
  pvector = new double[nextra];
  basepair  = new BasePair(lmp);
    }

/* ---------------------------------------------------------------------- */

Pair3spn2::~Pair3spn2()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(ljeps);
    memory->destroy(param);
    memory->destroy(angle);
    memory->destroy(bp_array);
    memory->destroy(R_solv);
    memory->destroy(lambda);
    memory->destroy(DGsolv);
    memory->destroy(volu);
  
    delete [] pvector;
  }
}

/* ---------------------------------------------------------------------- */

void Pair3spn2::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul, ebp, ecstk, ebasepair,
  ecrossstack, fpair, eexcl,engy;
  double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  int stea, steb, stec, sted, stee, stef, myitype, myjtype;
  int *ilist,*jlist,*numneigh,**firstneigh;

  eexcl = ecoul = ebp = ecstk = ebasepair = ecrossstack = engy = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int *tag = atom->tag;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  //double *special_coul = force->special_coul;
  //double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqr2e = force->qqr2e;
  double dr, drsi, r, rinv, screening;
  int itag, jtag, iflag,bpflag;
  double bp_cutoff = 0.8;
  nbps = 0;
  nbps_inter = 0; // Interstrand base pairs

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    itag = tag[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      //factor_lj = special_lj[sbmask(j)];
      //factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      dr = sqrt(rsq);
      drsi = 1.0 / rsq;
      jtype = type[j];
      jtag = tag[j];
      iflag = 0;
      bpflag = 0;

      // Set flags for determining pair interactions
      if (molecule[i] != molecule[j]) // If different molecules
      {
          iflag = 1;
          if (bp_array[itype-1][jtype-1]) bpflag = 1;
      }
      else
      {
        if (abs(jtag - itag) > 4) iflag = 1;
        if (bp_array[itype-1][jtype-1]) 
        {
            if (abs(jtag-itag) > 10) bpflag = 1;
        }
      }
      ebasepair = 0;
      ecrossstack = 0.0;
      forcecoul = 0.0;
      forcelj = 0.0;
      // If not within one nucleotide of each other
      if (iflag) {
        if (rsq < cutsq[itype][jtype]) {
            r2inv = 1.0/rsq;
            if (rsq < cut_ljsq[itype][jtype])
            {
                // Calculate the base pairing interaction
                if (bpflag)
                {
                    // Set the global IDs of the sites we are working with
                    steb = itag;
                    sted = jtag;

                    if (steb < sted)
                    {
                        stee = sted - 3;
                        stef = steb + 3;
                        stea = steb - 1;
                        stec = sted - 1;
                        myitype = itype;
                        myjtype = jtype;
                    }
                    else
                    {
                        int tmp = sted;
                        sted = steb;
                        steb = tmp;
                        stea = steb - 1;
                        stec = sted - 1;
                        stee = sted - 3;
                        stef = steb + 3;
                        myitype = jtype;
                        myjtype = itype;
                    }
                    // Now I map them to local coordinates
                    int index[6] = {atom->map(stea), atom->map(steb), atom->map(stec),atom->map(sted),atom->map(stee), atom->map(stef)};
                    int site_type[6], mytype;

                    for (int q = 0; q < SIX; q++)
                    {
                        if (index[q] != -1){
                          mytype = atom->type[index[q]];  
                          if (mytype <= 14){
                            //make sure its a dna site
                            site_type[q] = (mytype - 3) % 4; 
                          }
                          else{
                            site_type[q] = EMPTY;
                          }
                        }
                        else{
                          site_type[q] = EMPTY;
                        }
                    }
                   
                    // Initializing the distances used for cross stacking
                    basepair->assign(index, site_type ,param, angle, x);

                    // Calculating each cross stacking interaction

                    // Calculate first cross stacking interaction if site d is
                    // not a terminal base
                    if (!(bp_array[myitype-1][myjtype-1] % 2))
                    {
                        ecrossstack = basepair->cross_stacking(0,f);  // Interaction 1
                    }

                    // Calcuate second cross stacking interaction if site b is
                    // not a terminal base
                    if (!(bp_array[myitype-1][myjtype-1] % 3))
                    {
                        ecrossstack += basepair->cross_stacking(1,f);  // Interaction 2
                    }
                    
                    // Calculating the base pair interaction
                    ebasepair = basepair->base_pairing(f);  // Interaction 3

                    // Determine whether or not this is indeed a base pair
                    if (ebasepair <  -bp_cutoff * basepair->epsi[2]) nbps++;
                    if ((ebasepair <  -bp_cutoff * basepair->epsi[2]) && (molecule[i] != molecule[j])) nbps_inter++;
                }
                else
                {

                    // Excluded volume interaction
                    if (dr * _SF_ < sigma[itype][jtype])
                    {
                        r6inv = r2inv*r2inv*r2inv;
                        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
                        if (eflag) eexcl += r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) + ljeps[itype][jtype];
                    } 
                    else forcelj = 0.0;

                }
                if (eflag)
                {
                    ebp += ebasepair;
                    ecstk += ecrossstack;
                }
            }

            // Calculate the screened Coulombic interactions
            if (qtmp)
            {
                if (q[j])
                {
                    if (rsq < cut_coulsq[itype][jtype])
                    {
                          rinv = 1.0/dr;
                          double ldbi = 1.0 / ldby;
                          screening = exp(-dr * ldbi)/ dielectric;
                          engy = qqr2e * qtmp*q[j]*screening * rinv;
                          forcecoul = engy * dr * (ldbi + rinv);
                          if (eflag) ecoul += engy;
                    }
                }
            }
        }
        //if ((itype == 1) && (jtype == 1)) printf("forcecoul=%lf\n", forcecoul);

        fpair = (forcecoul + forcelj) * r2inv;
        /*
        if (fpair > 100.0)
        {
            printf("Large force i=%d, j=%d,itag=%d, jtag=%d,r=%lf\n",i,j,itag,jtag,dr);
            printf("forecoul=%lf, forcelj=%lf\n",forcecoul,forcelj);
        }
        */
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }
      }
     }
  }

  pvector[0] = ebp;
  pvector[1] = ecstk;
  pvector[2] = eexcl;
  pvector[3] = ecoul;
  pvector[4] = nbps;
  pvector[5] = nbps_inter;
  evdwl = ebp + ecstk + eexcl;

  // This correctly adds the potential energy.  However, the virial is not calculated properly
  if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,ecoul,fpair,delx,dely,delz);
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void Pair3spn2::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  n = 40; // Large enought to accomodate extra site types, including proteins

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(cut_coul,n+1,n+1,"pair:cut_coul");
  memory->create(cut_coulsq,n+1,n+1,"pair:cut_coulsq");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(ljeps,n+1,n+1,"pair:ljeps");
  memory->create(param,13,5,5,"pair:param");
  memory->create(angle,7,5,5,"pair:angle");
  memory->create(bp_array,n+1,n+1,"pair:bp_array");
  memory->create(R_solv,n+1,"pair:R_solv");
  memory->create(lambda,n+1,"pair:lambda");
  memory->create(DGsolv,n+1,"pair:DGsolv");
  memory->create(volu,n+1,"pair:volu");
  //
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void Pair3spn2::settings(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Illegal pair_style command");
  if (!allocated) allocate();
  int i,j,k;
  int n = 40; // Hardcoded number of sites

    /* Args: 

    Temperature - 
    Salt Concentration (mM)
    Short-range cutoff
    Long-range cutoff
    */
  sprintf(dna_type,"%s",arg[0]); // The y
  temp = force->numeric(FLERR,arg[1]);
  salt_conc = force->numeric(FLERR,arg[2]) / 1000.0; // Convert to [M]

  // Calculate solution dielectric (salt concentration is in M)
  dielectric = 249.4 - 7.88E-01 * temp + 7.20E-04 * temp * temp;
  dielectric *= (1.000 - (0.2551 * salt_conc) + (0.05151 * salt_conc * salt_conc) - (0.006889 * salt_conc * salt_conc * salt_conc));

  // Calculate Deybe length
  ldby = sqrt(dielectric * _PV_ * _KB_ * temp * 1.0E27 / (2.0 * _NA_ * _EC_ * _EC_ * salt_conc));

  cut_lj_global = force->numeric(FLERR,arg[3]);
  cut_coul_global = force->numeric(FLERR,arg[4]); // This cutoff must be the same as that for coul/long

    // Hard-coded values (used to be pair_style arguments
  scale_Santalucia = 0.70; //scale factor applied to Santalucia parameters
  hbond_scaling = 0.694; // scale factor applied to epsilon to get AT base pairing strength
  bp_ratio = 1.266; // ratio between G:C and A:T base pairing


  /* 
  With this information, we can proceeding to initialize a number 
  of important variable that must be hard-coded.
  */

    // Santalucia base step enthalpies as taken from J. Santalucia, PNAS, 1996.
    // A - 0, T - 1, G - 2, C - 3
  double Santalucia[4][4];
  Santalucia[0][0] = -7.6; // AA base step
  Santalucia[0][1] = -7.2;
  Santalucia[0][2] = -7.8;
  Santalucia[0][3] = -8.4;
  Santalucia[1][0] = -7.2;
  Santalucia[1][1] = -7.6;
  Santalucia[1][2] = -8.5;
  Santalucia[1][3] = -8.2;
  Santalucia[2][0] = -8.2;
  Santalucia[2][1] = -8.4;
  Santalucia[2][2] = -8.0;
  Santalucia[2][3] = -9.8;
  Santalucia[3][0] = -8.5;
  Santalucia[3][1] = -7.8;
  Santalucia[3][2] = -10.5;
  Santalucia[3][3] = -8.0;

  eps = 0.0;
  for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
          Santalucia[i][j] = -Santalucia[i][j] * scale_Santalucia;
          eps += Santalucia[i][j];
      }
  }
  eps /= 16.0;

 // Assign angles
 assign_angles(dna_type,angle);

 // Setting all of the equilibrium distances
 assign_distances(dna_type,param);


 // Now setting the values of alphas for the Morse potential
 for (i = 0; i < FOUR; i++) {
     for (j = 0; j < FOUR; j++) {
         param[3][i][j] = 4.0; // Cross stacking 1
         param[4][i][j] = 4.0; // Cross stacking 2
         param[5][i][j] = 2.0; // Base pairing
     }
 }

 // Now setting the ranges for angle-modulation
 for (i = 0; i < FOUR; i++) {
     for (j = 0; j < FOUR; j++) {
         param[6][i][j] = 8.0; // Cross stacking 1
         param[7][i][j] = 8.0; // Cross stacking 2
         param[8][i][j] = 12.0; // Base pairing
     }
 }

 // An array that we use to set the cross stacking values
 int WC[4] = {1,0,3,2};

 double hbond[4];
 double ref_hbond_energy, avg_bp_energy, cstk_scale_factor;
 
 double curv_scale_factor = 1.0;
 if (strcmp("bdna/curv",dna_type) == 0) curv_scale_factor = 0.861;

 ref_hbond_energy = eps * hbond_scaling;

 hbond[0] = ref_hbond_energy;
 hbond[1] = ref_hbond_energy;
 hbond[2] = ref_hbond_energy * bp_ratio;
 hbond[3] = ref_hbond_energy * bp_ratio;

 // Calculating the average base pairing strength
 avg_bp_energy = ref_hbond_energy * (1.0 + bp_ratio) / 2.0;

    // Specifying epsilon_BP
    param[11][0][0] = EMPTY;
    param[11][0][1] = ref_hbond_energy * curv_scale_factor;
    param[11][0][2] = EMPTY;
    param[11][0][3] = EMPTY;
    param[11][1][0] = ref_hbond_energy * curv_scale_factor;
    param[11][1][1] = EMPTY;
    param[11][1][2] = EMPTY;
    param[11][1][3] = EMPTY;
    param[11][2][0] = EMPTY;
    param[11][2][1] = EMPTY;
    param[11][2][2] = EMPTY;
    param[11][2][3] = bp_ratio * ref_hbond_energy * curv_scale_factor;
    param[11][3][0] = EMPTY;
    param[11][3][1] = EMPTY;
    param[11][3][2] = bp_ratio * ref_hbond_energy * curv_scale_factor;
    param[11][3][3] = EMPTY;
    //printf("AT=%lf\n", param[11][0][1]*4.184);

    double cstk[FOUR][FOUR];
   for (i = 0; i < FOUR; i++) {
        for (j = 0; j < FOUR; j++) {
              cstk[i][j] = Santalucia[i][j] - 0.5 * (hbond[i] + hbond[j]);
        }
   }

   double avg_cstk = 0.0;
   for (i = 0; i < FOUR; i++) {
        for (j = 0; j < FOUR; j++) {
            avg_cstk += 0.5 * cstk[WC[i]][j];
            avg_cstk += 0.5 * cstk[i][WC[j]];
        }
   }
   avg_cstk /= 32.0;

   cstk_scale_factor = avg_bp_energy * 100/88.0 * 0.12/avg_cstk * curv_scale_factor;
    
   for (i = 0; i < FOUR; i++) {
        for (j = 0; j < FOUR; j++) {
            // cstk_1
            param[9][i][j] = cstk_scale_factor * 0.5 * cstk[i][WC[j]];
            // cstk_2
            param[10][i][j] = cstk_scale_factor * 0.5 * cstk[WC[i]][j];
            //printf("cstk1[%d][%d]=%lf\n", i,j,param[9][i][j]*4.184);
        }
   }

  /* Initializing the base pair array (contains flags that determine whether 
  or not cross stacking interactions are calculated) */

    // Initialize all to zero
  for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            bp_array[i][j] = 0;
        }
  }
    
    /* 
    flag non-zero - calculate base pairing
       flag divisible by 2 calculate the base pair interaction b--e
       flag divisible by 3 calculate the base pair interaction d--f
    */

    bp_array[2][3] = bp_array[2][11] = 6; bp_array[2][7] = 3; // A
    bp_array[3][2] = bp_array[3][10] = 6; bp_array[3][6] = 3; // T
    bp_array[4][5] = bp_array[4][13] = 6; bp_array[4][9] = 3; // G
    bp_array[5][4] = bp_array[5][12] = 6; bp_array[5][8] = 3; // C
    bp_array[6][3] = bp_array[6][11] = 6; bp_array[6][7] = 3; // 5A
    bp_array[7][2] = bp_array[7][10] = 6; bp_array[7][6] = 3; // 5T
    bp_array[8][5] = bp_array[8][13] = 6; bp_array[8][9] = 3; // 5G
    bp_array[9][4] = bp_array[9][12] = 6; bp_array[9][8] = 3; // 5C
    bp_array[10][3] = bp_array[10][11] = 2; bp_array[10][7] = 7; // 3A
    bp_array[11][2] = bp_array[11][10] = 2; bp_array[11][6] = 7; // 3T
    bp_array[12][5] = bp_array[12][13] = 2; bp_array[12][9] = 7; // 3G
    bp_array[13][4] = bp_array[13][12] = 2; bp_array[13][8] = 7; // 3C
/*
  bp_array[14][14] =  {{0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // P
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // S
  {0,0,0,6,0,0,0,3,0,0,0,6,0,0}, // A 
  {0,0,6,0,0,0,3,0,0,0,6,0,0,0}, // T 
  {0,0,0,0,0,6,0,0,0,3,0,0,0,6}, // G 
  {0,0,0,0,6,0,0,0,3,0,0,0,6,0}, // C 
  {0,0,0,6,0,0,0,3,0,0,0,6,0,0}, // 5A
  {0,0,6,0,0,0,3,0,0,0,6,0,0,0}, // 5T
  {0,0,0,0,0,6,0,0,0,3,0,0,0,6}, // 5G
  {0,0,0,0,6,0,0,0,3,0,0,0,6,0}, // 5C
  {0,0,0,2,0,0,0,7,0,0,0,2,0,0}, //3A
  {0,0,2,0,0,0,7,0,0,0,2,0,0,0}, //3T
  {0,0,0,0,0,2,0,0,0,7,0,0,0,2}, //3G
  {0,0,0,0,2,0,0,0,7,0,0,0,2,0}}; //3C
  */

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_lj[i][j] = cut_lj_global;
          cut_coul[i][j] = cut_coul_global;

        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void Pair3spn2::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int n = atom->ntypes;
  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);

  double cut_lj_one = cut_lj_global;
  double cut_coul_one = cut_coul_global;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one * _SF_;
      cut_lj[i][j] = cut_lj_one;
      cut_coul[i][j] = cut_coul_one;
      setflag[i][j] = 1; 
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void Pair3spn2::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style 3spn2 requires atom attribute q");

  neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double Pair3spn2::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {

    // Apply mixing rules to sigma and epsilon
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut_lj[i][j] = mix_distance(cut_lj[i][i],cut_lj[j][j]);
    cut_coul[i][j] = mix_distance(cut_coul[i][i],cut_coul[j][j]);
  }

  double cut = MAX(cut_lj[i][j],cut_coul[i][j]);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];
  cut_coulsq[i][j] = cut_coul[i][j] * cut_coul[i][j];

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  ljeps[i][j] = epsilon[i][j];

  cut_ljsq[j][i] = cut_ljsq[i][j];
  cut_coulsq[j][i] = cut_coulsq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  ljeps[j][i] = ljeps[i][j];
  sigma[j][i] = sigma[i][j];

  // Set cross terms

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut_lj[i][j]*cut_lj[i][j]*cut_lj[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void Pair3spn2::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
        fwrite(&cut_coul[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void Pair3spn2::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&cut_lj[i][j],sizeof(double),1,fp);
          fread(&cut_coul[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_coul[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void Pair3spn2::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void Pair3spn2::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul_global,sizeof(double),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */


// This function is
double Pair3spn2::single(int i, int j, int itype, int jtype,
                                double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r2inv,r6inv,forcecoul,forcelj,phicoul,philj;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq[itype][jtype]) {
    forcecoul = force->qqrd2e * atom->q[i]*atom->q[j]*sqrt(r2inv);
  } else forcecoul = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  } else forcelj = 0.0;
  fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq[itype][jtype]) {
    phicoul = force->qqrd2e * atom->q[i]*atom->q[j]*sqrt(r2inv);
    eng += factor_coul*phicoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
    eng += factor_lj*philj;
  }

  return eng;
}

/* ---------------------------------------------------------------------- */

void *Pair3spn2::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul[1][1];
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  return NULL;
}

void Pair3spn2::assign_angles(char *dna_type, double***angles)
{
    if (strcmp(dna_type,"bdna") == 0) {
        // Angle between ba and dc
        angle[0][0][0] = EMPTY;
        angle[0][0][1] = 116.09;
        angle[0][0][2] = EMPTY;
        angle[0][0][3] = EMPTY;
        angle[0][1][0] = 116.09;
        angle[0][1][1] = EMPTY;
        angle[0][1][2] = EMPTY;
        angle[0][1][3] = EMPTY;
        angle[0][2][0] = EMPTY;
        angle[0][2][1] = EMPTY;
        angle[0][2][2] = EMPTY;
        angle[0][2][3] = 124.93;
        angle[0][3][0] = EMPTY;
        angle[0][3][1] = EMPTY;
        angle[0][3][2] = 124.93;
        angle[0][3][3] = EMPTY;

        // Defining the angles for cross stacking 1
        angle[1][0][0] = 154.38;
        angle[1][0][1] = 159.10;
        angle[1][0][2] = 152.46;
        angle[1][0][3] = 158.38;
        angle[1][1][0] = 147.10;
        angle[1][1][1] = 153.79;
        angle[1][1][2] = 144.44;
        angle[1][1][3] = 151.48;
        angle[1][2][0] = 154.69;
        angle[1][2][1] = 157.83;
        angle[1][2][2] = 153.43;
        angle[1][2][3] = 158.04;
        angle[1][3][0] = 152.99;
        angle[1][3][1] = 159.08;
        angle[1][3][2] = 150.53;
        angle[1][3][3] = 157.17;

        // Defining the angles for cross stacking 2
        angle[2][0][0] = 116.88;
        angle[2][0][1] = 121.74;
        angle[2][0][2] = 114.23;
        angle[2][0][3] = 119.06;
        angle[2][1][0] = 109.42;
        angle[2][1][1] = 112.95;
        angle[2][1][2] = 107.32;
        angle[2][1][3] = 110.56;
        angle[2][2][0] = 119.34;
        angle[2][2][1] = 124.72;
        angle[2][2][2] = 116.51;
        angle[2][2][3] = 121.98;
        angle[2][3][0] = 114.60;
        angle[2][3][1] = 118.26;
        angle[2][3][2] = 112.45;
        angle[2][3][3] = 115.88;

        // Torsion angles used in base pairing
        angle[3][0][0] = EMPTY;
        angle[3][0][1] = -38.35;
        angle[3][0][2] = EMPTY;
        angle[3][0][3] = EMPTY;
        angle[3][1][0] = -38.35;
        angle[3][1][1] = EMPTY;
        angle[3][1][2] = EMPTY;
        angle[3][1][3] = EMPTY;
        angle[3][2][0] = EMPTY;
        angle[3][2][1] = EMPTY;
        angle[3][2][2] = EMPTY;
        angle[3][2][3] = -42.98;
        angle[3][3][0] = EMPTY;
        angle[3][3][1] = EMPTY;
        angle[3][3][2] = -42.98;
        angle[3][3][3] = EMPTY;

        // Theta 1 angle - (base pairing)
        angle[4][0][0] = EMPTY;
        angle[4][0][1] = 156.54;
        angle[4][0][2] = EMPTY;
        angle[4][0][3] = EMPTY;
        angle[4][1][0] = 135.78;
        angle[4][1][1] = EMPTY;
        angle[4][1][2] = EMPTY;
        angle[4][1][3] = EMPTY;
        angle[4][2][0] = EMPTY;
        angle[4][2][1] = EMPTY;
        angle[4][2][2] = EMPTY;
        angle[4][2][3] = 159.81;
        angle[4][3][0] = EMPTY;
        angle[4][3][1] = EMPTY;
        angle[4][3][2] = 141.16;
        angle[4][3][3] = EMPTY;

        // Theta 2 angle (base pairing)
        angle[5][0][0] = EMPTY;
        angle[5][0][1] = 135.78;
        angle[5][0][2] = EMPTY;
        angle[5][0][3] = EMPTY;
        angle[5][1][0] = 156.54;
        angle[5][1][1] = EMPTY;
        angle[5][1][2] = EMPTY;
        angle[5][1][3] = EMPTY;
        angle[5][2][0] = EMPTY;
        angle[5][2][1] = EMPTY;
        angle[5][2][2] = EMPTY;
        angle[5][2][3] = 141.16;
        angle[5][3][0] = EMPTY;
        angle[5][3][1] = EMPTY;
        angle[5][3][2] = 159.81;
        angle[5][3][3] = EMPTY;
    } else if (strcmp("bdna/curv",dna_type) == 0) {
        angle[0][0][0] = EMPTY;
        angle[0][0][1] = 110.92; //was bug
        angle[0][0][2] = EMPTY;
        angle[0][0][3] = EMPTY;
        angle[0][1][0] = 110.92; //was bug
        angle[0][1][1] = EMPTY;
        angle[0][1][2] = EMPTY;
        angle[0][1][3] = EMPTY;
        angle[0][2][0] = EMPTY;
        angle[0][2][1] = EMPTY;
        angle[0][2][2] = EMPTY;
        angle[0][2][3] = 120.45; //was bug
        angle[0][3][0] = EMPTY;
        angle[0][3][1] = EMPTY;
        angle[0][3][2] = 120.45; //was bug
        angle[0][3][3] = EMPTY;


        angle[1][0][0] = 154.04;
        angle[1][0][1] = 158.77;
        angle[1][0][2] = 153.88;
        angle[1][0][3] = 157.69;
        angle[1][1][0] = 148.62;
        angle[1][1][1] = 155.05;
        angle[1][1][2] = 147.54;
        angle[1][1][3] = 153.61;
        angle[1][2][0] = 153.91;
        angle[1][2][1] = 155.72;
        angle[1][2][2] = 151.84;
        angle[1][2][3] = 157.80;
        angle[1][3][0] = 152.04;
        angle[1][3][1] = 157.72;
        angle[1][3][2] = 151.65;
        angle[1][3][3] = 154.49;

        angle[2][0][0] = 116.34;
        angle[2][0][1] = 119.61;
        angle[2][0][2] = 115.19;
        angle[2][0][3] = 120.92;
        angle[2][1][0] = 107.40;
        angle[2][1][1] = 110.76;
        angle[2][1][2] = 106.33;
        angle[2][1][3] = 111.57;
        angle[2][2][0] = 121.61;
        angle[2][2][1] = 124.92;
        angle[2][2][2] = 120.52;
        angle[2][2][3] = 124.88;
        angle[2][3][0] = 112.45;
        angle[2][3][1] = 115.43;
        angle[2][3][2] = 110.51;
        angle[2][3][3] = 115.80;

        // Torsion angles used in base pairing
        angle[3][0][0] = EMPTY;
        //angle[3][0][1] = -38.35;
        angle[3][0][1] = -38.18;
        angle[3][0][2] = EMPTY;
        angle[3][0][3] = EMPTY;
        //angle[3][1][0] = -38.35;
        angle[3][1][0] = -38.18;
        angle[3][1][1] = EMPTY;
        angle[3][1][2] = EMPTY;
        angle[3][1][3] = EMPTY;
        angle[3][2][0] = EMPTY;
        angle[3][2][1] = EMPTY;
        angle[3][2][2] = EMPTY;
        //angle[3][2][3] = -45.81;
        angle[3][2][3] = -35.75;
        angle[3][3][0] = EMPTY;
        angle[3][3][1] = EMPTY;
        //angle[3][3][2] = -45.81;
        angle[3][3][2] = -35.75;
        angle[3][3][3] = EMPTY;

        // Theta 1 angle - (base pairing)
        angle[4][0][0] = EMPTY;
        //angle[4][0][1] = 156.54;
        angle[4][0][1] = 153.17;
        angle[4][0][2] = EMPTY;
        angle[4][0][3] = EMPTY;
        //angle[4][1][0] = 135.78;
        angle[4][1][0] = 133.51;
        angle[4][1][1] = EMPTY;
        angle[4][1][2] = EMPTY;
        angle[4][1][3] = EMPTY;
        angle[4][2][0] = EMPTY;
        angle[4][2][1] = EMPTY;
        angle[4][2][2] = EMPTY;
        //angle[4][2][3] = 154.62;
        angle[4][2][3] = 159.50;
        angle[4][3][0] = EMPTY;
        angle[4][3][1] = EMPTY;
        //angle[4][3][2] = 152.74;
        angle[4][3][2] = 138.08;
        angle[4][3][3] = EMPTY;

        // Theta 2 angle (base pairing)
        angle[5][0][0] = EMPTY;
        //angle[5][0][1] = 135.78;
        angle[5][0][1] = 133.51;
        angle[5][0][2] = EMPTY;
        angle[5][0][3] = EMPTY;
        //angle[5][1][0] = 156.54;
        angle[5][1][0] = 153.17;
        angle[5][1][1] = EMPTY;
        angle[5][1][2] = EMPTY;
        angle[5][1][3] = EMPTY;
        angle[5][2][0] = EMPTY;
        angle[5][2][1] = EMPTY;
        angle[5][2][2] = EMPTY;
        //angle[5][2][3] = 152.74;
        angle[5][2][3] = 138.08;
        angle[5][3][0] = EMPTY;
        angle[5][3][1] = EMPTY;
        //angle[5][3][2] = 154.62;
        angle[5][3][2] = 159.50;
        angle[5][3][3] = EMPTY;

    
    } else if (strcmp("adna",dna_type) == 0) {
        // Angle between ba and dc
        angle[0][0][0] = EMPTY;
        angle[0][0][1] = 126.57;
        angle[0][0][2] = EMPTY;
        angle[0][0][3] = EMPTY;
        angle[0][1][0] = 126.57;
        angle[0][1][1] = EMPTY;
        angle[0][1][2] = EMPTY;
        angle[0][1][3] = EMPTY;
        angle[0][2][0] = EMPTY;
        angle[0][2][1] = EMPTY;
        angle[0][2][2] = EMPTY;
        angle[0][2][3] = 134.71;
        angle[0][3][0] = EMPTY;
        angle[0][3][1] = EMPTY;
        angle[0][3][2] = 134.71;
        angle[0][3][3] = EMPTY;

        // Defining the angles for cross stacking 1
        angle[1][0][0] = 147.44;
        angle[1][0][1] = 148.97;
        angle[1][0][2] = 146.21;
        angle[1][0][3] = 150.17;
        angle[1][1][0] = 138.42;
        angle[1][1][1] = 141.67;
        angle[1][1][2] = 136.64;
        angle[1][1][3] = 141.64;
        angle[1][2][0] = 147.67;
        angle[1][2][1] = 148.28;
        angle[1][2][2] = 146.84;
        angle[1][2][3] = 150.02;
        angle[1][3][0] = 145.83;
        angle[1][3][1] = 148.39;
        angle[1][3][2] = 144.24;
        angle[1][3][3] = 148.74;

        // Defining the angles for cross stacking 2
        angle[2][0][0] = 130.50;
        angle[2][0][1] = 138.73;
        angle[2][0][2] = 126.68;
        angle[2][0][3] = 134.18;
        angle[2][1][0] = 130.41;
        angle[2][1][1] = 134.68;
        angle[2][1][2] = 127.69;
        angle[2][1][3] = 131.38;
        angle[2][2][0] = 130.57;
        angle[2][2][1] = 140.17;
        angle[2][2][2] = 126.44;
        angle[2][2][3] = 135.31;
        angle[2][3][0] = 132.69;
        angle[2][3][1] = 138.21;
        angle[2][3][2] = 129.73;
        angle[2][3][3] = 134.45;

        // Torsion angles used in base pairing
        angle[3][0][0] = EMPTY;
        angle[3][0][1] = 50.17;
        angle[3][0][2] = EMPTY;
        angle[3][0][3] = EMPTY;
        angle[3][1][0] = 50.17;
        angle[3][1][1] = EMPTY;
        angle[3][1][2] = EMPTY;
        angle[3][1][3] = EMPTY;
        angle[3][2][0] = EMPTY;
        angle[3][2][1] = EMPTY;
        angle[3][2][2] = EMPTY;
        angle[3][2][3] = 38.33;
        angle[3][3][0] = EMPTY;
        angle[3][3][1] = EMPTY;
        angle[3][3][2] = 38.33;
        angle[3][3][3] = EMPTY;

        // Theta 1 angle - (base pairing)
        angle[4][0][0] = EMPTY;
        angle[4][0][1] = 160.91;
        angle[4][0][2] = EMPTY;
        angle[4][0][3] = EMPTY;
        angle[4][1][0] = 140.49;
        angle[4][1][1] = EMPTY;
        angle[4][1][2] = EMPTY;
        angle[4][1][3] = EMPTY;
        angle[4][2][0] = EMPTY;
        angle[4][2][1] = EMPTY;
        angle[4][2][2] = EMPTY;
        angle[4][2][3] = 165.25;
        angle[4][3][0] = EMPTY;
        angle[4][3][1] = EMPTY;
        angle[4][3][2] = 147.11;
        angle[4][3][3] = EMPTY;

        // Theta 2 angle (base pairing)
        angle[5][0][0] = EMPTY;
        angle[5][0][1] = 140.49;
        angle[5][0][2] = EMPTY;
        angle[5][0][3] = EMPTY;
        angle[5][1][0] = 160.91;
        angle[5][1][1] = EMPTY;
        angle[5][1][2] = EMPTY;
        angle[5][1][3] = EMPTY;
        angle[5][2][0] = EMPTY;
        angle[5][2][1] = EMPTY;
        angle[5][2][2] = EMPTY;
        angle[5][2][3] = 147.11;
        angle[5][3][0] = EMPTY;
        angle[5][3][1] = EMPTY;
        angle[5][3][2] = 165.25;
        angle[5][3][3] = EMPTY;

    } else error->all(FLERR,"Incorrect dna_type in pair_style 3spn2");
        // Convert all angles from degrees to radians
        int i, j, k;
    for (i = 0; i < SIX; i++) {
        for (j = 0; j < FOUR; j++) {
            for (k = 0; k < FOUR; k++) {
                angle[i][j][k] *= MY_PI/180.0;
            }
        }
    }
}

void Pair3spn2::assign_distances(char *dna_type, double***param)
{
    if (strcmp(dna_type,"bdna") == 0 ) {
        // Equilibrium distance for cross stacking 1
        param[0][0][0] = 6.208;
        param[0][0][1] = 6.876;
        param[0][0][2] = 6.072;
        param[0][0][3] = 6.811;
        param[0][1][0] = 6.876;
        param[0][1][1] = 7.480;
        param[0][1][2] = 6.771;
        param[0][1][3] = 7.453;
        param[0][2][0] = 6.072;
        param[0][2][1] = 6.771;
        param[0][2][2] = 5.921;
        param[0][2][3] = 6.688;
        param[0][3][0] = 6.811;
        param[0][3][1] = 7.453;
        param[0][3][2] = 6.688;
        param[0][3][3] = 7.409;


        // Equilibrium distance for cross stacking 2
        param[1][0][0] = 5.435;
        param[1][0][1] = 6.295;
        param[1][0][2] = 5.183;
        param[1][0][3] = 6.082;
        param[1][1][0] = 6.295;
        param[1][1][1] = 7.195;
        param[1][1][2] = 6.028;
        param[1][1][3] = 6.981;
        param[1][2][0] = 5.183;
        param[1][2][1] = 6.028;
        param[1][2][2] = 4.934;
        param[1][2][3] = 5.811;
        param[1][3][0] = 6.082;
        param[1][3][1] = 6.981;
        param[1][3][2] = 5.811;
        param[1][3][3] = 6.757;

        // Equilibrium distances for base pairing
        param[2][0][0] = EMPTY;
        param[2][0][1] = 5.941;
        param[2][0][2] = EMPTY;
        param[2][0][3] = EMPTY;
        param[2][1][0] = 5.941;
        param[2][1][1] = EMPTY;
        param[2][1][2] = EMPTY;
        param[2][1][3] = EMPTY;
        param[2][2][0] = EMPTY;
        param[2][2][1] = EMPTY;
        param[2][2][2] = EMPTY;
        param[2][2][3] = 5.530;
        param[2][3][0] = EMPTY;
        param[2][3][1] = EMPTY;
        param[2][3][2] = 5.530;
        param[2][3][3] = EMPTY;
    } else if (strcmp("bdna/curv",dna_type) == 0) {
        printf("getting bdna/curv dists!\n");
         //Cross stacking 1
        param[0][0][0] = 6.420;
        param[0][0][1] = 6.770;
        param[0][0][2] = 6.270;
        param[0][0][3] = 6.840;
        param[0][1][0] = 6.770;
        param[0][1][1] = 7.210;
        param[0][1][2] = 6.530;
        param[0][1][3] = 7.080;
        param[0][2][0] = 6.270;
        param[0][2][1] = 6.530;
        param[0][2][2] = 5.740;
        param[0][2][3] = 6.860;
        param[0][3][0] = 6.840;
        param[0][3][1] = 7.080;
        param[0][3][2] = 6.860;
        param[0][3][3] = 6.790;

        //Cross stacking 2
        param[1][0][0] = 5.580;
        param[1][0][1] = 6.140;
        param[1][0][2] = 5.630;
        param[1][0][3] = 6.180;
        param[1][1][0] = 6.140;
        param[1][1][1] = 6.800;
        param[1][1][2] = 6.070;
        param[1][1][3] = 6.640;
        param[1][2][0] = 5.630;
        param[1][2][1] = 6.070;
        param[1][2][2] = 5.870;
        param[1][2][3] = 5.660;
        param[1][3][0] = 6.180;
        param[1][3][1] = 6.640;
        param[1][3][2] = 5.660;
        param[1][3][3] = 6.800;


        param[2][0][0] = EMPTY;
        param[2][0][1] = 5.82;
        param[2][0][2] = EMPTY;
        param[2][0][3] = EMPTY;
        param[2][1][0] = 5.82;
        param[2][1][1] = EMPTY;
        param[2][1][2] = EMPTY;
        param[2][1][3] = EMPTY;
        param[2][2][0] = EMPTY;
        param[2][2][1] = EMPTY;
        param[2][2][2] = EMPTY;
        param[2][2][3] = 5.52;
        param[2][3][0] = EMPTY;
        param[2][3][1] = EMPTY;
        param[2][3][2] = 5.52;
        param[2][3][3] = EMPTY;
    } else if (strcmp("adna",dna_type) == 0) {
        // Equilibrium distance for cross stacking 1
        param[0][0][0] = 7.344;
        param[0][0][1] = 8.081;
        param[0][0][2] = 7.187;
        param[0][0][3] = 7.990;
        param[0][1][0] = 8.081;
        param[0][1][1] = 8.755;
        param[0][1][2] = 7.952;
        param[0][1][3] = 8.697;
        param[0][2][0] = 7.187;
        param[0][2][1] = 7.952;
        param[0][2][2] = 7.019;
        param[0][2][3] = 7.844;
        param[0][3][0] = 7.990;
        param[0][3][1] = 8.697;
        param[0][3][2] = 7.844;
        param[0][3][3] = 8.630;

        // Equilibrium distance for cross stacking 2
        param[1][0][0] = 4.624;
        param[1][0][1] = 5.095;
        param[1][0][2] = 4.464;
        param[1][0][3] = 5.162;
        param[1][1][0] = 5.095;
        param[1][1][1] = 5.693;
        param[1][1][2] = 4.896;
        param[1][1][3] = 5.724;
        param[1][2][0] = 4.464;
        param[1][2][1] = 4.896;
        param[1][2][2] = 4.315;
        param[1][2][3] = 4.968;
        param[1][3][0] = 5.162;
        param[1][3][1] = 5.724;
        param[1][3][2] = 4.968;
        param[1][3][3] = 5.759;

        // Equilibrium distances for base pairing
        param[2][0][0] = EMPTY;
        param[2][0][1] = 5.861;
        param[2][0][2] = EMPTY;
        param[2][0][3] = EMPTY;
        param[2][1][0] = 5.861;
        param[2][1][1] = EMPTY;
        param[2][1][2] = EMPTY;
        param[2][1][3] = EMPTY;
        param[2][2][0] = EMPTY;
        param[2][2][1] = EMPTY;
        param[2][2][2] = EMPTY;
        param[2][2][3] = 5.528;
        param[2][3][0] = EMPTY;
        param[2][3][1] = EMPTY;
        param[2][3][2] = 5.528;
        param[2][3][3] = EMPTY;
    } else error->all(FLERR,"Incorrect dna_type in pair_style 3spn2");
}

