# If you want to specify the package address
# you can add them to the PYTHONPATH environment variable.
# Also you can add them on the run time uncommenting the lines below
# import sys
# open3SPN2_HOME = '/Users/weilu/open3spn2/'
# openAWSEM_HOME = '/Users/weilu/openmmawsem/'
# sys.path.insert(0,open3SPN2_HOME)
# sys.path.insert(0,openAWSEM_HOME)

#Import openAWSEM, open3SPN2 and other libraries
import open3SPN2
import ffAWSEM
import pandas
import numpy as np
import simtk.openmm
from functools import partial
import sys

#Fix the system (adds missing atoms)
fix=open3SPN2.fixPDB("1lmb.pdb")

#Create a table containing both the proteins and the DNA
complex_table=open3SPN2.pdb2table(fix)

# Create a single memory file
ffAWSEM.create_single_memory(fix)

#Generate a coarse-grained model of the DNA molecules
dna_atoms=open3SPN2.DNA.CoarseGrain(complex_table)

#Generate a coarse-grained model of the Protein molecules
protein_atoms=ffAWSEM.Protein.CoarseGrain(complex_table)

#Merge the models
Coarse=pandas.concat([protein_atoms,dna_atoms],sort=False)
Coarse.index=range(len(Coarse))
Coarse['serial']=list(Coarse.index)

#Save the protein_sequence
ffAWSEM.save_protein_sequence(Coarse,sequence_file='protein.seq')

# Create a merged PDB
ffAWSEM.writePDB(Coarse,'clean.pdb')

#Create the merged system
pdb=simtk.openmm.app.PDBFile('clean.pdb')
top=pdb.topology
coord=pdb.positions
forcefield=simtk.openmm.app.ForceField(ffAWSEM.xml,open3SPN2.xml)
s=forcefield.createSystem(top)

#Create the DNA and Protein Objects
dna=open3SPN2.DNA.fromCoarsePDB('clean.pdb')
with open('protein.seq') as ps:
    protein_seq=ps.readlines()[0]
protein=ffAWSEM.Protein.fromCoarsePDB('clean.pdb',
                                      sequence=protein_seq)
dna.periodic=False
protein.periodic=False

#Copy the AWSEM parameter files
ffAWSEM.copy_parameter_files()

#Clear Forces from the system (optional)
keepCMMotionRemover=True
j=0
for i, f in enumerate(s.getForces()):
    if keepCMMotionRemover and i == 0 and f.__class__ == simtk.openmm.CMMotionRemover:
        # print('Kept ', f.__class__)
        j += 1
        continue
    else:
        # print('Removed ', f.__class__)
        s.removeForce(j)
if keepCMMotionRemover == False:
    assert len(s.getForces()) == 0, 'Not all the forces were removed'
else:
    assert len(s.getForces()) <= 1, 'Not all the forces were removed'

#Initialize the force dictionary
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

#Add AWSEM forces
ft=ffAWSEM.functionTerms
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
                      )
protein.setup_virtual_sites(s)

#Add DNA-protein interaction forces
for force_name in open3SPN2.protein_dna_forces:
    print(force_name)
    force = open3SPN2.protein_dna_forces[force_name](dna,protein)
    s.addForce(force)
    forces.update({force_name: force})

#Fix exclussions
for force_name in openAWSEMforces:
    print(force_name)
    if force_name in ['contact']:
        force = openAWSEMforces[force_name](protein, 
                                            withExclusion=False,
                                            periodic=False)
        print(force.getNumExclusions())
        open3SPN2.addNonBondedExclusions(dna,force)
        print(force.getNumExclusions())
    elif force_name in ['Excl']:
        force = openAWSEMforces[force_name](protein)
        print(force.getNumExclusions())
        open3SPN2.addNonBondedExclusions(dna,force)
        print(force.getNumExclusions())
    else:
        force = openAWSEMforces[force_name](protein)
    s.addForce(force)
    forces.update({force_name: force})

#Initialize the simulation
temperature=300 * simtk.openmm.unit.kelvin
platform_name='OpenCL' #'Reference','CPU','CUDA', 'OpenCL'
integrator = simtk.openmm.LangevinIntegrator(temperature, 
					     1 / simtk.openmm.unit.picosecond,
					     2 * simtk.openmm.unit.femtoseconds)
platform = simtk.openmm.Platform.getPlatformByName(platform_name)
simulation = simtk.openmm.app.Simulation(top,s, integrator, platform)
simulation.context.setPositions(coord)
energy_unit=simtk.openmm.unit.kilojoule_per_mole
state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(energy_unit)
print(energy)

#Obtain total energy
energy_unit=simtk.openmm.unit.kilojoule_per_mole
state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(energy_unit)
print('TotalEnergy',round(energy,6),energy_unit.get_symbol())

#Obtain detailed energy
energies = {}
for force_name, force in forces.items():
    group=force.getForceGroup()
    state = simulation.context.getState(getEnergy=True, 
                                        groups=2**group)
    energies[force_name] =state.getPotentialEnergy().value_in_unit(energy_unit)

for force_name in forces.keys():
    print(force_name, round(energies[force_name],6),
          energy_unit.get_symbol())

#Add simulation reporters
dcd_reporter=simtk.openmm.app.DCDReporter(f'output.dcd', 10000)
energy_reporter=simtk.openmm.app.StateDataReporter(sys.stdout, 10000, step=True,time=True, potentialEnergy=True, temperature=True)
simulation.reporters.append(dcd_reporter)
simulation.reporters.append(energy_reporter)

#Run simulation
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(100000)
