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
#include "stdlib.h"
#include "angle_3spn2_stacking.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001
#define _PI_ 3.14159265358979

/* ---------------------------------------------------------------------- */

Angle3spn2Stacking::Angle3spn2Stacking(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

Angle3spn2Stacking::~Angle3spn2Stacking()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(epsi);
    memory->destroy(sigm);
    memory->destroy(alpha);
    memory->destroy(theta0);
    memory->destroy(range);
  }
}

/* ---------------------------------------------------------------------- */

void Angle3spn2Stacking::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double estack,f1[3],f3[3];
  double dtheta,tk;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;
  double argu, fmorse, emorse, cosine, sine, cosine_term, prefactor, estck, rng, 
    dtha, eangle;

  estack = 0.0;
  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nanglelist; n++) {
    estck = 0.0;
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond
    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond
    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // angle (cos and sin)

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    // force & energy
    rng = range[type];

    dtha = acos(c) - theta0[type];
    
    if (r2 < sigm[type])
    {
        // A purely repulsive Morse potential...
        argu = alpha[type] * (r2 - sigm[type]);
        fmorse = -2.0 * alpha[type] * epsi[type] * exp(-argu) * (1.0 - exp(-argu)) / r2;
        // Apply the repulsive force to each of the sites
        if (newton_bond || i2 < nlocal) {
            f[i2][0] -= fmorse * delx2;
            f[i2][1] -= fmorse * dely2;
            f[i2][2] -= fmorse * delz2;
        }

        if (newton_bond || i3 < nlocal) {
            f[i3][0] += fmorse * delx2;
            f[i3][1] += fmorse * dely2;
            f[i3][2] += fmorse * delz2;
        }

        if (eflag) estck += epsi[type] * (1.0 - exp(-argu)) * (1.0 - exp(-argu));
    }

    // Now we calculate a modulated attractive interaction between these sites

    // Calculate the unmodulate Morse potential
    if ( (dtha >= -_PI_/(rng*2.0)) && (dtha <= _PI_/(rng*2.0)))
    {
        // Calculate pair-wise, attractive-only Morse potential
        if (r2 >= sigm[type])
        {
            argu = alpha[type] * (r2 - sigm[type]);
            fmorse = -2.0 * alpha[type] * epsi[type] * exp(-argu) * (1.0 - exp(-argu)) / r2;
            emorse = epsi[type] * (1.0 - exp(-argu)) * (1.0 - exp(-argu)) - epsi[type];
        }
        else
        {
            fmorse = 0.0;
            emorse = -epsi[type];
        }
        // Apply the attractive force to each of the sites
        if (newton_bond || i2 < nlocal) {
            f[i2][0] -= fmorse * delx2;
            f[i2][1] -= fmorse * dely2;
            f[i2][2] -= fmorse * delz2;
        }

        if (newton_bond || i3 < nlocal) {
            f[i3][0] += fmorse * delx2;
            f[i3][1] += fmorse * dely2;
            f[i3][2] += fmorse * delz2;
        }

        if (eflag) estck += emorse;
    }
    // If the angle falls within the "cone"
    else if (((dtha >= _PI_/(rng*2.0)) && (dtha <= _PI_/rng)) 
        || ((dtha <= -_PI_/(rng*2.0)) && (dtha >= -_PI_/rng)))
    {
        // Calculate pair-wise, attractive-only Morse potential
        if (r2 >= sigm[type])
        {
            argu = alpha[type] * (r2 - sigm[type]);
            fmorse = -2.0 * alpha[type] * epsi[type] * exp(-argu) * (1.0 - exp(-argu)) / r2;
            emorse = epsi[type] * (1.0 - exp(-argu)) * (1.0 - exp(-argu)) - epsi[type];
        }
        else
        {
            fmorse = 0.0;
            emorse = -epsi[type];
        }
        cosine = cos(rng * dtha);
        sine = sin(rng*dtha);
        cosine_term = 1.0 - cosine * cosine;
        
        if (eflag) estck += cosine_term * emorse; 
        
        prefactor = 2.0 * rng * cosine * sine * 1.0 /sqrt(1.0-c*c);

        a = -prefactor * emorse;
        a11 = a*c / rsq1;
        a12 = -a / (r1*r2);
        a22 = a*c / rsq2;

        f1[0] = a11*delx1 + a12*delx2;
        f1[1] = a11*dely1 + a12*dely2;
        f1[2] = a11*delz1 + a12*delz2;
        f3[0] = a22*delx2 + a12*delx1;
        f3[1] = a22*dely2 + a12*dely1;
        f3[2] = a22*delz2 + a12*delz1;

        if (newton_bond || i1 < nlocal) {
            f[i1][0] += f1[0];
            f[i1][1] += f1[1];
            f[i1][2] += f1[2];
        }

        if (newton_bond || i2 < nlocal) {
            f[i2][0] -= (f1[0] + f3[0] + cosine_term * delx2 * fmorse);
            f[i2][1] -= (f1[1] + f3[1] + cosine_term * dely2 * fmorse);
            f[i2][2] -= (f1[2] + f3[2] + cosine_term * delz2 * fmorse);
        }

        if (newton_bond || i3 < nlocal) {
            f[i3][0] += f3[0] + cosine_term * delx2 * fmorse;
            f[i3][1] += f3[1] + cosine_term * dely2 * fmorse;
            f[i3][2] += f3[2] + cosine_term * delz2 * fmorse;
        }
    }
    // If the angle falls outside the "cone"
    else
    {
        // Do nothing...
    }
    if (eflag) eangle = estck;

    // Need to determine how to calculate this virial
    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle, f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
}

/* ---------------------------------------------------------------------- */

void Angle3spn2Stacking::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(epsi,n+1,"angle:epsi");
  memory->create(sigm,n+1,"angle:sigm");
  memory->create(alpha,n+1,"angle:alpha");
  memory->create(theta0,n+1,"angle:theta0");
  memory->create(range,n+1,"angle:range");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void Angle3spn2Stacking::coeff(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nangletypes,ilo,ihi);

  double epsi_one = force->numeric(FLERR,arg[1]);
  double sigm_one = force->numeric(FLERR,arg[2]);
  double theta0_one = force->numeric(FLERR,arg[3]);
  double alpha_one = force->numeric(FLERR,arg[4]);
  double range_one = force->numeric(FLERR,arg[5]);

  // convert theta0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    epsi[i] = epsi_one;
    sigm[i] = sigm_one;
    alpha[i] = alpha_one;
    theta0[i] = theta0_one/180.0 * MY_PI;
    range[i] = range_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double Angle3spn2Stacking::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void Angle3spn2Stacking::write_restart(FILE *fp)
{
  fwrite(&epsi[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&sigm[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&alpha[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&theta0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&range[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void Angle3spn2Stacking::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&epsi[1],sizeof(double),atom->nangletypes,fp);
    fread(&sigm[1],sizeof(double),atom->nangletypes,fp);
    fread(&alpha[1],sizeof(double),atom->nangletypes,fp);
    fread(&theta0[1],sizeof(double),atom->nangletypes,fp);
    fread(&range[1],sizeof(double),atom->nangletypes,fp);
  }
  MPI_Bcast(&epsi[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigm[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&alpha[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&range[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */
double Angle3spn2Stacking::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
  double emorse, rng, argu, cosine, cosine_term;
  double estck =0.0;

  double dtha = acos(c) - theta0[type];
  rng = range[type];

    if (r2 < sigm[type])
    {
        // A purely repulsive Morse potential...
        argu = alpha[type] * (r2 - sigm[type]);
        estck = epsi[type] * (1.0 - exp(-argu)) * (1.0 - exp(-argu));
    }

    // Now we calculate a modulated attractive interaction between these sites

    // Calculate the unmodulate Morse potential
    if ( (dtha >= -_PI_/(rng*2.0)) && (dtha <= _PI_/(rng*2.0)))
    {
        // Calculate pair-wise, attractive-only Morse potential
        if (r2 < sigm[type])
        {
            emorse = epsi[type];
        }
        else
        {
            argu = alpha[type] * (r2 - sigm[type]);
            emorse = epsi[type] * (1.0 - exp(-argu)) * (1.0 - exp(-argu)) - epsi[type];
        }
        estck += emorse;

    }
    // If the angle falls within the "cone"
    else if (((dtha >= _PI_/(rng*2.0)) && (dtha <= _PI_/rng)) 
        || ((dtha <= -_PI_/(rng*2.0)) && (dtha >= -_PI_/rng)))
    {
        // Calculate pair-wise, attractive-only Morse potential
        if (r2 < sigm[type])
        {
            emorse = -epsi[type];
        }
        else
        {
            argu = alpha[type] * (r2 - sigm[type]);
            emorse = epsi[type] * (1.0 - exp(-argu)) * (1.0 - exp(-argu)) - epsi[type];
        }
        cosine = cos(rng * dtha);
        cosine_term = 1.0 - cosine * cosine;
        emorse *= cosine_term;
        
        estck += emorse; // I think that estack needs to be added to possible energies to write
    }
    else
    {
        // Do nothing...
    }

  return estck;
}
