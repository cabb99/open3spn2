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

#ifndef LMP_EXTRA_BASEPAIR_H
#define LMP_EXTRA_BASEPAIR_H
#endif

#include "pointers.h"
namespace LAMMPS_NS {

class BasePair : protected Pointers {
    public:
    BasePair(class LAMMPS *);
        virtual ~BasePair();
        void assign(int*, int *, double ***, double***, double**);
        double base_pairing(double**);
        double cross_stacking(int,double**);
        double epsi[3];
            
    protected:
        int stea, steb, stec, sted, stee, stef;
        double dbax, dbay, dbaz, ddcx, ddcy, ddcz, ddbx, ddby, ddbz,
        dbex, dbey, dbez, ddfx, ddfy, ddfz, ebax, ebay, ebaz;
        double dbai, eabx, eaby, eabz, ddci, edcx, edcy, edcz,
        ddbi, edbx, edby, edbz, dbei, ebex, ebey, ebez, 
        ddfi, edfx, edfy, edfz;
        double dba, ddc, ddb, dbe, ddf, ddbsi,  dbesi, ddfsi, dbi;
        double padx, pady, padz, cosb, isb2, isnb, pbcx, pbcy, 
        pbcz, cosd, isd2, isnd, pacx, pacy, pacz, cnum;
        double cphi[6], sigm[3], dtha[6], ateq[6],range[6],alpha[6];
        void mors_norp(double , double, double, double, double*, double*);
        void mors_rp(double , double, double, double, double*, double*);
    private:
};

}

