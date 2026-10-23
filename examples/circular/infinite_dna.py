"""Simulation of an infinite DNA with open3SPN2.

Builds a straight 105 bp duplex (10 helical turns) along z whose ends are bonded to their own periodic images
through the box, so the DNA has no ends, and runs a short Langevin dynamics in that periodic box. Each strand winds
an integer number of turns per period; delta_linking_number adds or removes turns, which the periodic box keeps
fixed during the simulation (torsionally constrained DNA)."""
import numpy as np
import openmm.app
import openmm.unit as unit
import open3SPN2

rng = np.random.default_rng(0)
sequence = ''.join(rng.choice(list('ATGC'), 105))

dna = open3SPN2.DNA.fromSequence(sequence, infinite=True, delta_linking_number=0,
                                 output_pdb='infinite.pdb', atomistic_template='infinite_atomistic.pdb')
print(f'Periodic box: {dna.periodic_box[0]:.1f} x {dna.periodic_box[1]:.1f} x {dna.periodic_box[2]:.2f} nm')
system = open3SPN2.System(dna, periodicBox=dna.periodic_box)
system.add3SPN2forces()

platform = 'CPU'  # 'OpenCL' or 'CUDA' are much faster
simulation = system.initializeMD(temperature=300 * unit.kelvin, platform_name=platform)
simulation.minimizeEnergy()
# Unwrapped coordinates keep the helix continuous when viewed; in VMD, `pbc box` draws the box
simulation.reporters.append(openmm.app.DCDReporter('infinite.dcd', 1000, enforcePeriodicBox=False))
simulation.reporters.append(openmm.app.StateDataReporter('infinite.log', 1000, step=True, potentialEnergy=True,
                                                         temperature=True))
simulation.step(10000)
