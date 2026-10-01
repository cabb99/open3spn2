"""Simulation of a supercoiled plasmid with open3SPN2.

Builds a 500 bp circular DNA from its sequence, removing 3 turns from its relaxed linking number
(negative supercoiling), and runs a short Langevin dynamics while checking that the linking number
of the two strands, a topological invariant, does not change."""
import numpy as np
import openmm.app
import openmm.unit as unit
import open3SPN2

rng = np.random.default_rng(0)
sequence = ''.join(rng.choice(list('ATGC'), 500))

dna = open3SPN2.DNA.fromSequence(sequence, circular=True, delta_linking_number=-3,
                                 output_pdb='plasmid.pdb', atomistic_template='plasmid_atomistic.pdb')
dna.periodic = False
system = open3SPN2.System(dna, periodicBox=None)
system.add3SPN2forces()

platform = 'CPU'  # 'OpenCL' or 'CUDA' are much faster
simulation = system.initializeMD(temperature=300 * unit.kelvin, platform_name=platform)
simulation.minimizeEnergy()
simulation.reporters.append(openmm.app.DCDReporter('plasmid.dcd', 1000))
simulation.reporters.append(openmm.app.StateDataReporter('plasmid.log', 1000, step=True, potentialEnergy=True,
                                                         temperature=True))
for block in range(10):
    simulation.step(1000)
    positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    lk = dna.computeLinkingNumber(positions.value_in_unit(unit.angstrom))
    print(f'step {simulation.currentStep}: linking number {lk:.2f}')
