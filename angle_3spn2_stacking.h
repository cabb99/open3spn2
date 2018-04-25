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

#ifdef ANGLE_CLASS
AngleStyle(stacking/3spn2,Angle3spn2Stacking)
#else

#ifndef LMP_ANGLE_3SPN2_STACKING_H
#define LMP_ANGLE_3SPN2_STACKING_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class Angle3spn2Stacking : public Angle {
 public:
  Angle3spn2Stacking(class LAMMPS *);
  virtual ~Angle3spn2Stacking();
  virtual void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);

 protected:
  double *epsi, *sigm, *alpha, *theta0, *range;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
