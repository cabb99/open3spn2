#!/usr/bin/env python
# coding: utf-8

# If you want to specify the package address
# you can add them to the PYTHONPATH environment variable.
# Also you can add them on the run time uncommenting the lines below
import sys
import os
# open3SPN2_HOME = '/Users/weilu/open3spn2/'
openAWSEM_HOME = '/home/sl206/Programs/openawsem'
# sys.path.insert(0,open3SPN2_HOME)
sys.path.insert(0,openAWSEM_HOME)

#sys.path.append('/home/sl206/miniconda3/envs/openmm/lib/python3.6')
#sys.path.append('/home/sl206/miniconda3/envs/openmm/lib/python3.6/site-packages')
#sys.path.append('/home/sl206/miniconda3/pkgs')

import argparse

#Import openAWSEM, open3SPN2 and other libraries
import pandas as pd
import numpy as np
import openmm
#import openmm

from functools import partial
import sys

import open3SPN2
import openawsem

import openmm.app
import openmm.unit


def printEnergy(simulation, forces):
    # #Total energy
    energy_unit=openmm.unit.kilocalorie_per_mole
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(energy_unit)
    print('TotalEnergy',round(energy,6),energy_unit.get_symbol())

    # #Detailed energy
    energies = {}
    for force_name, force in forces.items():
        group=force.getForceGroup()
        state = simulation.context.getState(getEnergy=True, groups=2**group)
        energies[force_name] =state.getPotentialEnergy().value_in_unit(energy_unit)

    for force_name in forces.keys():
        print(force_name, round(energies[force_name],6),energy_unit.get_symbol())

def write(message, output):
    with open(output, 'a') as f:
        f.write(message)
        f.write('\n')


def writeEnergy(simulation, forces, output):
    # #Total energy
    energy_unit=openmm.unit.kilocalorie_per_mole
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(energy_unit)
    print('TotalEnergy',round(energy,6),energy_unit.get_symbol())
    with open(output, 'a') as f:
        f.write(f'TotalEnergy {round(energy,6)} {energy_unit.get_symbol()}')
        f.write('\n')

    # #Detailed energy
    energies = {}

    for force_name, force in forces.items():
        group=force.getForceGroup()
        state = simulation.context.getState(getEnergy=True, groups=2**group)
        energies[force_name] =state.getPotentialEnergy().value_in_unit(energy_unit)

    for force_name in forces.keys():
        with open(output, 'a') as f:
            f.write(f'{force_name} {round(energies[force_name],6)} {energy_unit.get_symbol()}')
            f.write('\n')

def savePDB(toPath, simulation, PDBfile_name):
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    with open(os.path.join(toPath, PDBfile_name), "w") as pdb_file:
        openmm.app.PDBFile.writeFile(simulation.topology, positions, file=pdb_file)

def run(args):
    proteinDNA = args.proteinDNA


    #Create the merged system
    pdb=openmm.app.PDBFile(f'{proteinDNA}.pdb')
    top=pdb.topology
    coord=pdb.positions
    forcefield=openmm.app.ForceField(openawsem.xml,open3SPN2.xml)
    s=forcefield.createSystem(top)

    #Create the DNA and Protein Objects
    dna=open3SPN2.DNA.fromCoarsePDB(f'{proteinDNA}.pdb')
    #dna.computeTopology(template_from_X3DNA=True)
    with open('protein.seq') as ps:
        protein_sequence_one=ps.readlines()[0]
    protein=openawsem.Protein.fromCoarsePDB(f'{proteinDNA}.pdb',sequence=protein_sequence_one)
    dna.periodic=False
    protein.periodic=False
    #Don't activate this below. Appears not to apply if you have Protein.
    #s=open3SPN2.System(dna, periodicBox=None) 

    print(s.getForces())

    forces={}
    for i in range(s.getNumForces()):
        force = s.getForce(i)
        force_name="CMMotionRemover"

    #Add 3SPN2 forces
    for force_name in open3SPN2.forces:
        print(force_name)
        force = open3SPN2.forces[force_name](dna)
        if force_name in ['BasePair','CrossStacking']:
            force.addForce(s)
        else:
            s.addForce(force)
        forces.update({force_name:force})

    #Add AWSEM forces. Fragment memories are in the protein residue only AWSEM-created folder
    frags_dir = args.AWSEM
    openAWSEMforces = dict(Connectivity=openawsem.functionTerms.basicTerms.con_term,
                        Chain=openawsem.functionTerms.basicTerms.chain_term,
                        Chi=openawsem.functionTerms.basicTerms.chi_term,
                        Excl=openawsem.functionTerms.basicTerms.excl_term,
                        rama=openawsem.functionTerms.basicTerms.rama_term,
                        rama_pro=openawsem.functionTerms.basicTerms.rama_proline_term,
                        contact=openawsem.functionTerms.contactTerms.contact_term,
                        frag  = partial(openawsem.functionTerms.templateTerms.fragment_memory_term, 
                                        frag_file_list_file=f"{frags_dir}/{args.fragment}", 
                                        UseSavedFragTable=False, 
                                        k_fm=0.04184/3),
                        beta1 = openawsem.functionTerms.hydrogenBondTerms.beta_term_1,
                        beta2 = openawsem.functionTerms.hydrogenBondTerms.beta_term_2,
                        beta3 = openawsem.functionTerms.hydrogenBondTerms.beta_term_3,
                        pap1 = partial(openawsem.functionTerms.hydrogenBondTerms.pap_term_1,
                                        ssweightFileName=f"{frags_dir}/ssweight"),
                        pap2 = partial(openawsem.functionTerms.hydrogenBondTerms.pap_term_2,
                                        ssweightFileName=f"{frags_dir}/ssweight"),
                        DH = partial(openawsem.functionTerms.debyeHuckelTerms.debye_huckel_term, 
                                        chargeFile=f"{frags_dir}/charge.txt")
                        )

    protein.setup_virtual_sites(s)

    #Add DNA-protein interaction forces
    for force_name in open3SPN2.protein_dna_forces:
        print(force_name)
        force = open3SPN2.protein_dna_forces[force_name](dna,protein)
        s.addForce(force)
        forces.update({force_name: force})  
        
    #OpenAWSEM forces with exclusions
    for force_name in openAWSEMforces:
        print(force_name)
        if force_name in ['contact']:
            force = openAWSEMforces[force_name](protein, withExclusion=False,periodic=False)
            print(force_name, "pre-add #Exclusions", force.getNumExclusions())
            open3SPN2.addNonBondedExclusions(dna,force)
            print(force_name, "post-add #Exclusions", force.getNumExclusions())
        elif force_name in ['Excl']:
            force = openAWSEMforces[force_name](protein)
            print(force_name, "pre-add #Exclusions", force.getNumExclusions())
            open3SPN2.addNonBondedExclusions(dna,force)
            print(force_name, "post-add #Exclusions", force.getNumExclusions())
        #continue
        else:
            force = openAWSEMforces[force_name](protein)
        s.addForce(force)
        forces.update({force_name: force})

    # #Initialize Molecular Dynamics simulations

    temperature=args.tempStart * openmm.unit.kelvin
    Tstart = args.tempStart
    output = args.output
    platform_name=args.Platform #'Reference','CPU','CUDA', 'OpenCL'

    integrator = openmm.LangevinIntegrator(temperature, 1 / openmm.unit.picosecond, 2 * openmm.unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName(platform_name)
    #platform.setPropertyDefaultValue('OpenCLPlatformIndex', '0')
    #platform.setPropertyDefaultValue('DeviceIndex', args.device)

    simulation = openmm.app.Simulation(top,s, integrator, platform)
    simulation.context.setPositions(coord)
    printEnergy(simulation, forces)
    write("Starting structure", output)
    writeEnergy(simulation, forces, output)

    #reporter_frequency = 1000
    #append = False
    toPath = "."
    #print("reporter_frequency", reporter_frequency)
    #pdb_reporter=openmm.app.PDBReporter(os.path.join(toPath, "movie.pdb"), reporter_frequency)
    #dcd_reporter=openmm.app.DCDReporter(os.path.join(toPath, "output.dcd"), reporter_frequency, append=False)
    #energy_reporter=openmm.app.StateDataReporter(sys.stdout, reporter_frequency, step=True,time=True, potentialEnergy=True, temperature=True)
    #output_reporter=openmm.app.StateDataReporter(os.path.join(toPath, "output.log"), reporter_frequency, step=True,time=True, potentialEnergy=True, temperature=True)
    #checkpoint_reporter=openmm.app.CheckpointReporter(os.path.join(toPath, "checkpoint.chk"), reporter_frequency)
    #simulation.reporters.append(pdb_reporter)
    #simulation.reporters.append(dcd_reporter)
    #simulation.reporters.append(energy_reporter)
    #simulation.reporters.append(output_reporter)
    #simulation.reporters.append(checkpoint_reporter)

    # initial minimization block
    print("minization start")
    integrator = openmm.CustomIntegrator(0.001)
    simulation = openmm.app.Simulation(top,s, integrator, platform)
    simulation.context.setPositions(coord)
    print("Initial energies")
    printEnergy(simulation, forces)
    write("", output)
    write("Initial energies", output)
    writeEnergy(simulation, forces, output)
    savePDB(toPath, simulation, PDBfile_name = "init.pdb")
    simulation.minimizeEnergy() 
    print("Initial min energies")
    printEnergy(simulation, forces)
    write("", output)
    write("Initial min energies", output)
    writeEnergy(simulation, forces, output)
    savePDB(toPath, simulation, PDBfile_name = "init_min.pdb")
    dcd_reporter=openmm.app.DCDReporter(os.path.join(toPath, "output.dcd"), 1)
    simulation.reporters.append(dcd_reporter)
    simulation.step(1)
    # MD minimization block
    integrator = openmm.LangevinIntegrator(Tstart*openmm.unit.kelvin, 1/openmm.unit.picosecond, args.timeStep*openmm.unit.femtoseconds)
    simulation = openmm.app.Simulation(top,s, integrator, platform)
    #simulation.context.setPositions(coord)  # set the initial positions of the atoms
    print("Now T = 300 K energies")
    init_min_pdb = openmm.app.PDBFile(os.path.join(toPath, "init_min.pdb"))
    simulation.context.setPositions(init_min_pdb.positions)
    simulation.context.setVelocitiesToTemperature(Tstart*openmm.unit.kelvin)
    printEnergy(simulation, forces)
    write("", output)
    write("Now T = 300 K energies", output)
    writeEnergy(simulation, forces, output)
    simulation.minimizeEnergy()  # first, minimize the energy to a local minimum to reduce any large forces that might be present
    savePDB(toPath, simulation, PDBfile_name = "MD_min.pdb")
    print("Now T = 300 K min energies")
    printEnergy(simulation, forces)
    write("", output)
    write("Now T = 300 K min energies", output)
    writeEnergy(simulation, forces, output)
    simulation.step(1)
    print("minization end")
    write("minimization end", output)

def main():
    # from run_parameter import *
    parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatic copy the template file, \
        run simulations")

    parser.add_argument("proteinDNA", help="The name of the proteinDNA system")
    parser.add_argument("--tempStart", type=float, default=300, help="Starting temperature")
    parser.add_argument("-l", "--fragment", type=str, default="./frags.mem", help="Fragment memory (single or std)")  #temporary placeholder
    parser.add_argument("-a", "--AWSEM", type=str, default="./", help="protein-only AWSEM folder, should have fragment library") #not temporary
    parser.add_argument("-o", "--output", type=str, default="energy_output.log", help="Output file.")
    parser.add_argument("-p", "--Platform", type=str, default="OpenCL", help="platform to use")
    parser.add_argument("--timeStep", type=int, default=2)
    args = parser.parse_args()


    with open('energy_args.txt', 'a') as f:
        f.write(' '.join(sys.argv))
        f.write('\n')
    print(' '.join(sys.argv))

    run(args)

if __name__=="__main__":
    main()