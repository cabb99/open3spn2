import sys
import os

#Import openAWSEM, open3SPN2 and other libraries
import pandas as pd
import numpy as np
import openmm

from functools import partial

import open3SPN2
import openawsem

import openmm.app
import openmm.unit

import time

import importlib.util

import argparse

def printEnergy(simulation, forces):
    # #Total energy
    energy_unit=openmm.unit.kilocalorie_per_mole
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(energy_unit)
    print('Caution! The energy terms with identical energy values are in the same forceGroup!')
    print('TotalEnergy',round(energy,6),energy_unit.get_symbol())

    # #Detailed energy
    energies = {}
    for force_name, force in forces.items():
        group=force.getForceGroup()
        state = simulation.context.getState(getEnergy=True, groups=2**group)
        energies[force_name] =state.getPotentialEnergy().value_in_unit(energy_unit)

    for force_name in forces.keys():
        print(force_name, round(energies[force_name],6),energy_unit.get_symbol())

def savePDB(toPath, simulation, PDBfile_name):
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    with open(os.path.join(toPath, PDBfile_name), "w") as pdb_file:
        openmm.app.PDBFile.writeFile(simulation.topology, positions, file=pdb_file)

def run(args):
    simulation_platform = args.platform
    platform = openmm.Platform.getPlatformByName(simulation_platform)
    if simulation_platform == "CPU":
        if args.thread != -1:
            platform.setPropertyDefaultValue("Threads", str(args.thread))
        print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")

    #aries specific block
    elif simulation_platform == "OpenCL":
        platform.setPropertyDefaultValue('OpenCLPlatformIndex', '0')
        platform.setPropertyDefaultValue('DeviceIndex', args.device)
    
    pwd = os.getcwd()
    toPath = os.path.abspath(args.to)
    forceSetupFile = args.forces

    if args.to != "./":
        # os.system(f"mkdir -p {args.to}")
        os.makedirs(toPath, exist_ok=True)
        os.system(f"cp {forceSetupFile} {toPath}/{forceSetupFile}")

    checkPointPath = None if args.fromCheckPoint is None else os.path.abspath(args.fromCheckPoint)
    
    proteinDNA = args.proteinDNA

    pdb=openmm.app.PDBFile(f'{proteinDNA}.pdb')
    top=pdb.topology
    coord=pdb.positions
    forcefield=openmm.app.ForceField(openawsem.xml,open3SPN2.xml)
    s=forcefield.createSystem(top, removeCMMotion=not args.removeCMMotionRemover)

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

    #forces={}

    print(f"using force setup file from {forceSetupFile}")
    spec = importlib.util.spec_from_file_location("forces", forceSetupFile)
    # print(spec)
    forces_file = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(forces_file)
    forces = forces_file.set_up_forces(s,protein, dna, computeQ = False)

    # #Initialize Molecular Dynamics simulations
    snapShotCount = args.Frames
    stepsPerT = int(args.steps/snapShotCount)
    Tstart = args.tempStart
    Tend = args.tempEnd
    if args.reportFrequency == -1:
        if stepsPerT == 0:
            reporter_frequency = 4000
        else:
            reporter_frequency = stepsPerT
    else:
        reporter_frequency = args.reportFrequency
    # reporter_frequency = 4000
    #temperature=300 * openmm.unit.kelvin

    if args.fromCheckPoint:
        integrator = openmm.LangevinIntegrator(Tstart*openmm.unit.kelvin, 1/openmm.unit.picosecond, args.timeStep*openmm.unit.femtoseconds)
        simulation = openmm.app.Simulation(top,s, integrator, platform)
        simulation.loadCheckpoint(checkPointPath)
    else:
        # initial minimization block
        print("minization start")
        integrator = openmm.CustomIntegrator(0.001)
        simulation = openmm.app.Simulation(top,s, integrator, platform)
        simulation.context.setPositions(coord)
        print("Initial energies")
        printEnergy(simulation, forces)
        savePDB(toPath, simulation, PDBfile_name = "init.pdb")
        simulation.minimizeEnergy() 
        print("Initial min energies")
        printEnergy(simulation, forces)
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
        simulation.minimizeEnergy()  # first, minimize the energy to a local minimum to reduce any large forces that might be present
        savePDB(toPath, simulation, PDBfile_name = "MD_min.pdb")
        print("Now T = 300 K min energies")
        printEnergy(simulation, forces)
        simulation.step(1)
        print("minization end")

    #simulation = openmm.app.Simulation(top,s, integrator, platform)
    
    energy_unit=openmm.unit.kilocalorie_per_mole
    #state = simulation.context.getState(getEnergy=True)
    #energy = state.getPotentialEnergy().value_in_unit(energy_unit)
    #print(energy)
    # #Set initial positions
    #simulation.context.setPositions(s.coord.getPositions())

    print("reporter_frequency", reporter_frequency)
    pdb_reporter=openmm.app.PDBReporter(os.path.join(toPath, "movie.pdb"), reporter_frequency)
    dcd_reporter=openmm.app.DCDReporter(os.path.join(toPath, "output.dcd"), reporter_frequency, append=True)
    energy_reporter=openmm.app.StateDataReporter(sys.stdout, reporter_frequency, step=True,time=True, potentialEnergy=True, temperature=True)
    output_reporter=openmm.app.StateDataReporter(os.path.join(toPath, "output.log"), reporter_frequency, step=True,time=True, potentialEnergy=True, temperature=True, append=True)
    checkpoint_reporter=openmm.app.CheckpointReporter(os.path.join(toPath, "checkpoint.chk"), reporter_frequency)
    simulation.reporters.append(pdb_reporter)
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(energy_reporter)
    simulation.reporters.append(output_reporter)
    simulation.reporters.append(checkpoint_reporter)

    print("Simulation Start")
    if args.simulation_mode == 0:
        simulation.step(int(args.steps))
    elif args.simulation_mode == 1:
        deltaT = (Tend - Tstart) / snapShotCount
        for i in range(snapShotCount):
            integrator.setTemperature((Tstart + deltaT*i)*openmm.unit.kelvin)
            simulation.step(stepsPerT)
    print("Simulation finished")

    print("Analysis start")
    #os.system(f"{sys.executable} protein_DNA_analysis.py {args.proteinDNA}.pdb -t {os.path.join(toPath, "output.dcd")} -a {args.AWSEM} -l {args.fragment} -o {os.path.join(toPath, "info.dat")}")
    os.system(f"{sys.executable} protein_DNA_analysis.py {args.proteinDNA}.pdb -t {os.path.join(toPath, 'output.dcd')} -f {args.forces} -o {os.path.join(toPath, 'info.dat')}")
    print("Analysis finished")

def main():
    # from run_parameter import *
    parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatic copy the template file, \
        run simulations")

    parser.add_argument("proteinDNA", help="The name of the proteinDNA system")
    parser.add_argument("--to", default="./", help="location of movie file")
    #parser.add_argument("-c", "--chain", type=str, default="-1")
    parser.add_argument("-t", "--thread", type=int, default=-1, help="default is using all that is available")
    parser.add_argument("-p", "--platform", type=str, default="OpenCL")
    parser.add_argument("-s", "--steps", type=int, default=10000000, help="step size, default 10,000,000")
    parser.add_argument("--tempStart", type=float, default=300, help="Starting temperature")
    parser.add_argument("--tempEnd", type=float, default=300, help="Ending temperature")
    parser.add_argument("--fromCheckPoint", type=str, default=None, help="The checkpoint file you want to start from")
    parser.add_argument("-m", "--simulation_mode", type=int, default=0,
                    help="default 1,\
                            0: constant temperature,\
                            1: temperature annealing")
    parser.add_argument("--subMode", type=int, default=-1)
    #parser.add_argument("-f", "--forces", default="forces_setup.py")
    #parser.add_argument("--parameters", default=None)
    parser.add_argument("-r", "--reportFrequency", type=int, default=-1, help="default value step/400")
    #parser.add_argument("--fromOpenMMPDB", action="store_true", default=False)
    #parser.add_argument("--fasta", type=str, default="crystal_structure.fasta")
    parser.add_argument("--timeStep", type=int, default=2)
    #parser.add_argument("--includeLigands", action="store_true", default=False)
    parser.add_argument('--Frames', type=int, default=400, help="Number of frames")
    parser.add_argument('--device',default='0')
    parser.add_argument("-f", "--forces", default="forces_setup.py", type=str, help="forces setup file")
    parser.add_argument('--removeCMMotionRemover', action="store_true", default=False, help='Removes CMMotionRemover. Recommended for periodic boundary conditions and membrane simulations')
    #parser.add_argument("-l", "--fragment", type=str, default="./frags.mem", help="Fragment memory (single or std)")  #temporary placeholder
    #parser.add_argument("-a", "--AWSEM", type=str, default="./", help="protein-only AWSEM folder, should have fragment library") #not temporary
    args = parser.parse_args()


    with open('commandline_args.txt', 'a') as f:
        f.write(' '.join(sys.argv))
        f.write('\n')
    print(' '.join(sys.argv))

    run(args)

if __name__=="__main__":

    main()
