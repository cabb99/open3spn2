#!/bin/sh

# Calls to


if [ $# -ne 1 ]; then
    echo "Usage: $0 <sequence file>"
    exit 1
fi

ICNF=../DSIM_ICNF/icnf.exe
UTILS=../utils/

echo "Making parameter file"
python ${UTILS}/make_bp_params.py $1

echo "Running X3DNA"
x3dna_utils cp_std BDNA
rebuild -atomic bp_step.par atomistic.pdb
mkdir basis
mv Atomic* basis/

echo "Mapping to CG coordinates"
python ${UTILS}/pdb2cg_dna.py atomistic.pdb
mv in00_conf.xyz bdna_curv.xyz

echo "Building LAMMPS input file"
${ICNF} $1 1 1 . 0

echo "Replacing atoms in configuration file"
${UTILS}/replace_atoms.sh conf_lammps.in dna_conf.in bdna_curv_conf.in

echo "Making list files"
python ${UTILS}/make_list_files.py bdna_curv_conf.in

