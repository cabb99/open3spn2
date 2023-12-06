#!/usr/bin/env python3
import os
import argparse
import mdtraj as md
import pandas as pd
import simtk.openmm
import open3SPN2
import openawsem
from functools import partial

parser = argparse.ArgumentParser()
parser.add_argument("--protein", help="The name of the protein", default="./clean.pdb")
parser.add_argument("-t", "--trajectory", type=str, default="./output.dcd")
parser.add_argument("-o", "--output", type=str, default=None, help="The Name of file that show your energy.")
args = parser.parse_args()

trajectoryPath = os.path.abspath(args.trajectory)
if args.output is None:
    outFile = os.path.join(os.path.dirname(trajectoryPath), "info.dat")
else:
    outFile = os.path.join(os.path.dirname(trajectoryPath), args.output)

# fix=open3SPN2.fixPDB(args.protein)
fix=open3SPN2.fixPDB(args.protein)

#Create a table containing both the proteins and the DNA
complex_table=open3SPN2.pdb2table(fix)

#Generate a coarse-grained model of the Protein molecules
protein_atoms=openawsem.Protein.CoarseGrain(complex_table)

#Create the merged system
pdb=simtk.openmm.app.PDBFile(args.protein)
top=pdb.topology
coord=pdb.positions
forcefield=simtk.openmm.app.ForceField(openawsem.xml,open3SPN2.xml)
s=forcefield.createSystem(top)

#Create the DNA and Protein Objects
dna=open3SPN2.DNA.fromCoarsePDB(args.protein)
with open('protein.seq') as ps:
    protein_seq=ps.readlines()[0]
protein=openawsem.Protein.fromCoarsePDB(args.protein,
                                      sequence=protein_seq)
dna.periodic=False
protein.periodic=False

#Initialize the force dictionary
forces={}
for i in range(s.getNumForces()):
    force = s.getForce(i)
    force_name="CMMotionRemover"

#Add 3SPN2 forces
for force_name in open3SPN2.forces:
    # print(force_name)
    force = open3SPN2.forces[force_name](dna)
    if force_name in ['BasePair','CrossStacking']:
        force.addForce(s)
    else:
        s.addForce(force)
    forces[force_name] = force

#Add AWSEM forces
ft=openawsem.functionTerms
openAWSEMforces = dict(Connectivity=ft.basicTerms.con_term,
                       Chain=ft.basicTerms.chain_term,
                       Chi=ft.basicTerms.chi_term,
                       Excl=ft.basicTerms.excl_term,
                       rama=ft.basicTerms.rama_term,
                       rama_pro=ft.basicTerms.rama_proline_term,
                       contact=ft.contactTerms.contact_term,
                       frag  = partial(ft.templateTerms.fragment_memory_term,
                                       frag_file_list_file = "./single_frags.mem",
                                       npy_frag_table = "./single_frags.npy",
                                       UseSavedFragTable = False,
                                       k_fm = 0.04184/3),
                       beta1 = ft.hydrogenBondTerms.beta_term_1,
                       beta2 = ft.hydrogenBondTerms.beta_term_2,
                       beta3 = ft.hydrogenBondTerms.beta_term_3,
                       pap1 = ft.hydrogenBondTerms.pap_term_1,
                       pap2 = ft.hydrogenBondTerms.pap_term_2,
                       qval = partial(ft.biasTerms.q_value,
                                      reference_pdb_file = "crystal-structure.pdb")
                      )
protein.setup_virtual_sites(s)

#Add DNA-protein interaction forces
for force_name in open3SPN2.protein_dna_forces:
    # print(force_name)
    force = open3SPN2.protein_dna_forces[force_name](dna,protein)
    s.addForce(force)
    forces[force_name] = force

#Fix exclussions
for force_name in openAWSEMforces:
    # print(force_name)
    if force_name in ['contact']:
        force = openAWSEMforces[force_name](protein, 
                                            withExclusion=False,
                                            periodic=False)
        # print(force.getNumExclusions())
        open3SPN2.addNonBondedExclusions(dna,force)
        # print(force.getNumExclusions())
    elif force_name in ['Excl']:
        force = openAWSEMforces[force_name](protein)
        # print(force.getNumExclusions())
        open3SPN2.addNonBondedExclusions(dna,force)
        # print(force.getNumExclusions())
    else:
        force = openAWSEMforces[force_name](protein)
    s.addForce(force)
    forces[force_name] = force

#Initialize the simulation
temperature=300 * simtk.openmm.unit.kelvin
platform_name='OpenCL' #'Reference','CPU','CUDA', 'OpenCL'
# platform_name='Reference'
integrator = simtk.openmm.LangevinIntegrator(temperature, 
					     1 / simtk.openmm.unit.picosecond,
					     2 * simtk.openmm.unit.femtoseconds)
platform = simtk.openmm.Platform.getPlatformByName(platform_name)
simulation = simtk.openmm.app.Simulation(top,s, integrator, platform)
simulation.context.setPositions(coord)
energy_unit=simtk.openmm.unit.kilojoule_per_mole

trajectory = md.load(args.trajectory, top=args.protein)

energy_data = []
for frame in trajectory:
    simulation.context.setPositions(frame.xyz[0])
    #Obtain total energy
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(energy_unit)

    # Collect energies
    energies = {}

    for force_name, force in forces.items():
        group = force.getForceGroup()
        state = simulation.context.getState(getEnergy=True, groups=2**group)
        energies[force_name] = state.getPotentialEnergy().value_in_unit(energy_unit)

    energy_data.append({"TotalEnergy": energy, **energies})

energy_df = pd.DataFrame(energy_data)
showAll = {"TotalEnergy": energy, **energies}

# write energy_df into info.dat
with open(outFile, "w") as out:
    line = " ".join(["{0:<8s}".format(i) for i in ["Steps"] + list(showAll.keys())])
    print(line)
    out.write(line+"\n")

    for step, e in enumerate(energy_data):
        line = " ".join([f"{step:<8}"] + ["{0:<8.2f}".format(i) for i in e.values()])
        print(line)
        out.write(line+"\n")