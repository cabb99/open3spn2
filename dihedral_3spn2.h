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

#ifdef DIHEDRAL_CLASS

DihedralStyle(3spn2,Dihedral3spn2)

#else

#ifndef LMP_DIHEDRAL_3SPN2_H
#define LMP_DIHEDRAL_3SPN2_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class Dihedral3spn2 : public Dihedral {
 public:
  Dihedral3spn2(class LAMMPS *);
  virtual ~Dihedral3spn2();
  virtual void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 protected:
  double *k, *phi,  *sigm;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Dihedral problem: %d %ld %d %d %d %d

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

E: Incorrect args for dihedral coefficients

Self-explanatory.  Check the input script or data file.

E: Incorrect sign arg for dihedral coefficients

Self-explanatory.  Check the input script or data file.

E: Incorrect multiplicity arg for dihedral coefficients

Self-explanatory.  Check the input script or data file.

*/
