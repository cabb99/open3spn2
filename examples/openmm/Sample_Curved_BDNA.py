import sys
#print(sys.path)
sys.path.insert(0, '../..')
#print(sys.path)
import ff3SPN2
from simtk.unit import *
import simtk.openmm
import numpy as np
import sys

dna=ff3SPN2.DNA.fromXYZ('in00_conf.xyz','B_curved')
system=ff3SPN2.System(dna)
size=2*200*nanometer
centered=np.array(system.coord.getPositions())-np.array(system.coord.getPositions()).mean(axis=0)+np.array(size/2)
system.setDefaultPeriodicBoxVectors(*np.diag([size]*3))
system.add3SPN2forces()

temperature=300*kelvin
integrator = simtk.openmm.LangevinIntegrator(temperature, 1/500 / picosecond, 0.02 * picoseconds)
platform = simtk.openmm.Platform.getPlatformByName('Reference')
simulation = simtk.openmm.app.Simulation(system.coord.topology, system._wrapped_system, integrator, platform)
simulation.context.setPositions(system.coord.positions)
simulation.context.setPositions(centered)

simulation.reporters.append(simtk.openmm.app.PDBReporter(f'output.pdb', 1000),)
simulation.reporters.append(simtk.openmm.app.DCDReporter(f'output.dcd', 1000),)
#sim_out=open('sim_out.txt','w+')
simulation.reporters.append(simtk.openmm.app.StateDataReporter(sys.stdout, 1000, step=True,time=True,potentialEnergy=True, temperature=True,separator='\t',))
simulation.reporters.append(simtk.openmm.app.StateDataReporter('sim.log', 1000, step=True,time=True,totalEnergy=True,
                                              kineticEnergy=True,potentialEnergy=True, temperature=True))

simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(100000)
