#!/bin/sh

# Contains all of the code that is used to optimize the parameters
TOOLS=/media/d27443fd-635b-479e-b466-f65859915fae/Ions/relative_entropy/tools
NPAIRS=15
DIELECTRIC=78.0
PAIRLIST=(P-P P-Na P-Mg P-Cl Na-Na Na-Mg Na-Cl Mg-Mg Mg-Cl Cl-Cl P-N Na-N Mg-N Cl-N  N-N)
QIQJ=(1.0 -1.0 -2.0 1.0 1.0 2.0 -1.0 4.0 -2.0 1.0 -1.0 1.0 2.0 -1.0 1.0)

for j in `seq 0 $(($NPAIRS-1))`; do
    PAIR=${PAIRLIST[${j}]}
    Q=${QIQJ[${j}]}
    ${TOOLS}/make_lammps_table_full_pote.exe ../shifted_param/${PAIR}.param 550 ${DIELECTRIC} ${Q}
done
