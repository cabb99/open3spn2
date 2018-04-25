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

#ifdef PAIR_CLASS

PairStyle(3spn2,Pair3spn2)

#else

#ifndef LMP_PAIR_3SPN2_H
#define LMP_PAIR_3SPN2_H

#include "pair.h"
#include "base_pair.h"

namespace LAMMPS_NS {

class Pair3spn2 : public Pair {
 public:
  Pair3spn2(class LAMMPS *);
  virtual ~Pair3spn2();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);
  int nbps, nbps_inter;
 
 private:
  void assign_angles(char*, double***);
  void assign_distances(char*, double***);

 protected:
  double kappa_global, cut_lj_global,cut_coul_global, ebp_global,
  bp_ratio_global;
  double scale_Santalucia, hbond_scaling, bp_ratio,eps;
  double dielectric, ldby, salt_conc, temp;
  char dna_type[1024];
  double **cut_lj,**cut_ljsq;
  double **cut_vdW,**cut_vdWsq;
  double **cut_coul,**cut_coulsq;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4, **ljeps;
  double ***param;
  double ***angle;
  int **bp_array;
  double *R_solv, *lambda, *DGsolv, *volu;
  class BasePair *basepair;
  


  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style 3spn2 requires atom attribute q

The atom style defined does not have this attribute.

*/
