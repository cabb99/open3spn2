#!/usr/bin/env python3
import os
import argparse
import mdtraj as md
import pandas as pd
import openmm
import open3SPN2
import openawsem
from functools import partial
import importlib.util

import openmm.app
import openmm.unit
import sys
import numpy as np

import argparse

def run(args):
    trajectoryPath = os.path.abspath(args.trajectory)
    if args.output is None:
        outFile = os.path.join(os.path.dirname(trajectoryPath), "info.dat")
    else:
        outFile = os.path.join(os.path.dirname(trajectoryPath), args.output)

    simulation_platform = args.platform
    platform = openmm.Platform.getPlatformByName(simulation_platform)

    #aries specific block
    if simulation_platform == "OpenCL":
        platform.setPropertyDefaultValue('OpenCLPlatformIndex', '0')
        platform.setPropertyDefaultValue('DeviceIndex', args.device)

    # fix=open3SPN2.fixPDB(args.protein)
    fix=open3SPN2.fixPDB(args.proteinDNA)

    #Create a table containing both the proteins and the DNA
    complex_table=open3SPN2.pdb2table(fix)

    #Generate a coarse-grained model of the Protein molecules
    protein_atoms=openawsem.Protein.CoarseGrain(complex_table)

    #Create the merged system
    pdb=openmm.app.PDBFile(args.proteinDNA)
    top=pdb.topology
    coord=pdb.positions
    forcefield=openmm.app.ForceField(openawsem.xml,open3SPN2.xml)
    s=forcefield.createSystem(top)

    #Create the DNA and Protein Objects
    dna=open3SPN2.DNA.fromCoarsePDB(args.proteinDNA)
    with open('protein.seq') as ps:
        protein_seq=ps.readlines()[0]
    protein=openawsem.Protein.fromCoarsePDB(args.proteinDNA,
                                        sequence=protein_seq)
    dna.periodic=False
    protein.periodic=False

    #Initialize the force dictionary
    forceSetupFile = args.forces
    #forces={}

    print(f"using force setup file from {forceSetupFile}")
    spec = importlib.util.spec_from_file_location("forces", forceSetupFile)
    # print(spec)
    forces_file = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(forces_file)
    forces = forces_file.set_up_forces(s,protein, dna, computeQ = True)

    #Initialize the simulation
    temperature=300 * openmm.unit.kelvin
    platform_name=args.platform #'Reference','CPU','CUDA', 'OpenCL'
    integrator = openmm.LangevinIntegrator(temperature, 
                            1 / openmm.unit.picosecond,
                            2 * openmm.unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName(platform_name)
    simulation = openmm.app.Simulation(top,s, integrator, platform)
    simulation.context.setPositions(coord)
    energy_unit=openmm.unit.kilojoule_per_mole

    trajectory = md.load(args.trajectory, top=args.proteinDNA)

    forceGroupTable = {
                "Q": 1,
                "Rg": 2,
                # "Rg_bias": 5,
                "Bond": 6,
                "Angle": 7,
                "Stacking": 8,
                "Dihedral": 9,
                "BasePair": 10,
                "CrossStacking": 11,
                "ExclusionDNA": 12,
                "ElectrostaticsDNA": 13,
                "ExclusionProteinDNA": 14,
                "ElectrostaticsProteinDNA": 15,
                # "Reserved_for_direct_protein_DNA_readout_interactions": 16,
                # "Burial": 17,
                # "Helix_orientation": 18,
                # "Pulling": 19,
                "Backbone": 20,
                "Rama": 21,
                "Contact": 22,
                "Fragment": 23,
                # "Membrane": 24,
                # "ER": 25,
                # "TBM_Q": 26,
                "Beta": 27,
                "Pap": 28,
                "Helical": 29,
                "Debye_Huckel": 30,
                "AMH-Go": 31,
                "Total_Energy": list(range(5, 32))
                }
    print("Verify forceGroupTable in protein_DNA_analysis is set up correctly.")
    print("Total Energy includes forceGroups from 5 to 31, both sides inclusive.")
    showValue = ["Q", "Rg"]
    showEnergy = [
                # "Rg_bias",
                "Bond", 
                "Angle",
                "Stacking",
                "Dihedral",
                "BasePair",
                "CrossStacking",
                "ExclusionDNA",
                "ElectrostaticsDNA",
                "ExclusionProteinDNA",
                "ElectrostaticsProteinDNA",
                # "Reserved_for_direct_protein_DNA_readout_interactions", 
                # "Burial",
                # "Helix_orientation",
                # "Pulling", 
                "Backbone",
                "Rama",
                "Contact",
                "Fragment",
                # "Membrane",
                # "ER",
                # "TBM_Q",
                "Beta",
                "Pap",
                "Helical",
                "Debye_Huckel",
                "AMH-Go",
                "Total_Energy"
    ]
    showAll = showValue + showEnergy

    print("Printing energies")

    with open(outFile, "w") as out:
        line = " ".join(["{0:<8s}".format(i) for i in ["Steps"] + showAll])
        print(line)
        out.write(line+"\n")
        # for step, pdb in enumerate(pdb_trajectory):
        #     simulation.context.setPositions(pdb.positions)
        for step in range(len(trajectory)):
            simulation.context.setPositions(trajectory.openmm_positions(step))
            e = []
            for term in showAll:
                if type(forceGroupTable[term]) == list:
                    g = set(forceGroupTable[term])
                elif forceGroupTable[term] == -1:
                    g = -1
                else:
                    g = {forceGroupTable[term]}
                state = simulation.context.getState(getEnergy=True, groups=g)
                # if term == "Q" or term == "Rg" or term == "Qc" or term == "Q_wat" or term == "Q_mem":
                if term in showValue:
                    termEnergy = state.getPotentialEnergy().value_in_unit(openmm.unit.kilojoule_per_mole)
                else:
                    termEnergy = state.getPotentialEnergy().value_in_unit(openmm.unit.kilocalories_per_mole)
                e.append(termEnergy)
        #     print(*e)
            line = " ".join([f"{step:<8}"] + ["{0:<8.2f}".format(i) for i in e])
            print(line)
            out.write(line+"\n")
        #         print(forceGroupTable[term], state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("proteinDNA", help="The name of the protein", default="./clean.pdb")
    parser.add_argument("-t", "--trajectory", type=str, default="./output.dcd")
    parser.add_argument("-o", "--output", type=str, default=None, help="The Name of file that show your energy.")
    parser.add_argument("-p", "--platform", type=str, default="OpenCL", help="Could be OpenCL, CUDA and CPU")
    parser.add_argument('--device',default='0')
    parser.add_argument("-f", "--forces", default="forces_setup.py", type=str, help="forces setup file")
    #parser.add_argument("-l", "--fragment", type=str, default="./frags.mem", help="Fragment memory")  #temporary placeholder
    #parser.add_argument("-a", "--AWSEM", type=str, default="./", help="protein-only AWSEM folder, should have fragment library") #not temporary
    args = parser.parse_args()

    with open('analysis_commandline_args.txt', 'a') as f:
        f.write(' '.join(sys.argv))
        f.write('\n')
    print(' '.join(sys.argv))

    run(args)

if __name__=="__main__":
    main()