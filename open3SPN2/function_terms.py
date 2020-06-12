
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np


def dna_constraint_term(dna, k=10*kilocalorie_per_mole, forceGroup=20):
    print("dna_constraint_term is on")
    k = k.value_in_unit(kilojoule_per_mole)
    data = dna.atoms
    bonds_list = []
    DNA_resNames = ["DA", "DC", "DT", "DG"]
    interacting_atom = ["A", "C", "G", "T"]
    a = data.query("resname in @DNA_resNames").query("name in @interacting_atom").reset_index(drop=True)
    n_atom = len(a)
    for i in range(n_atom):
        dna_res_1 = a.iloc[i]
        chainID_1 = dna_res_1["chainID"]
        for j in range(i+1, n_atom):
            dna_res_2 = a.iloc[j]
            chainID_2 = dna_res_2["chainID"]
            if chainID_1 == chainID_2:
                continue
            atom1_pos = dna_res_1[["x", "y", "z"]].values
            atom2_pos = dna_res_2[["x", "y", "z"]].values
            dis = atom1_pos - atom2_pos
            r = (dis[0]**2 + dis[1]**2 + dis[2]**2)**0.5
            r /= 10   # convert to nm
            if r > 0.65:
                continue
            # print(i, j, r)
            bonds_list.append([dna_res_1["serial"], dna_res_2["serial"], r])
    bondForce = CustomBondForce(f"{k}*(r-r0)^2")
    bondForce.addPerBondParameter('r0')
    for (atom1, atom2, r) in bonds_list:
        # print(atom1, atom2, r)
        bondForce.addBond(int(atom1), int(atom2), [r])

    bondForce.setForceGroup(forceGroup)
    return bondForce
