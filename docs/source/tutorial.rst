Tutorial
================================================

    .. index:: single: Feynman
    
    *“What I cannot create, I do not understand.”* — Feynman
    
DNA from sequence
------------------

This section allows you to quickly simulate a piece of DNA. 

First import the ff3SPN2 module.

.. testcode::
    
    import ff3SPN2

If you obtain a ``ModuleNotFoundError`` error then you may need to check that ``ff3SPN2`` is in the installation path as detailed in :doc:`installation`.

After importing the module, follow the next steps to create a DNA system and add the 3SPN2 forces.
Make sure you have installed X3DNA before this step.

.. testsetup:: *

    import ff3SPN2

.. testcode::
    
    # Initialize the DNA from a sequence.
    # DNA type can be changed to 'A' or 'B'
    
    seq='ATACAAAGGTGCGAGGTTTCTATGCTCCCACG'
    dna=ff3SPN2.DNA.fromSequence(seq,dna_type='B_curved')
    
    # Compute the topology for the DNA structure.
    # Since the dna was generated from the sequence using X3DNA,
    # it is not necesary to recompute the geometry.
    
    dna.computeTopology(template_from_X3DNA=False)
    
    # Create the system.
    # To set periodic boundary conditions (periodicBox=[50,50,50]).
    # The periodic box size is in nanometers.
    dna.periodic=False
    s=ff3SPN2.System(dna, periodicBox=None)
    
    #Add 3SPN2 forces
    s.add3SPN2forces(verbose=True)

.. testoutput::
    
    Bond
    Angle
    Stacking
    Dihedral
    BasePair
    CrossStacking
    Exclusion
    Electrostatics

To set up the simulation you will need the openmm package as detailed in :doc:`installation`.

.. testcode::

    import simtk.openmm
    import simtk.openmm.app
    import simtk.unit
    import sys
    import numpy as np
    
    #Initialize Molecular Dynamics simulations
    s.initializeMD(temperature=300 * simtk.unit.kelvin,platform_name='OpenCL')
    simulation=s.simulation
    
    #Set initial positions
    simulation.context.setPositions(s.coord.getPositions())
    
    energy_unit=simtk.openmm.unit.kilojoule_per_mole
    #Total energy
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(energy_unit)
    print('TotalEnergy',round(energy,6),energy_unit.get_symbol())
    
    #Detailed energy
    energies = {}
    for force_name, force in s.forces.items():
        group=force.getForceGroup()
        state = simulation.context.getState(getEnergy=True, groups=2**group)
        energies[force_name] =state.getPotentialEnergy().value_in_unit(energy_unit)

    for force_name in s.forces.keys():
        print(force_name, round(energies[force_name],6),energy_unit.get_symbol())

.. testoutput::
    
    TotalEnergy -2046.579102 kJ/mol
    Bond 2.1e-05 kJ/mol
    Angle 0.001745 kJ/mol
    Stacking -596.408569 kJ/mol
    Dihedral -839.999756 kJ/mol
    BasePair -518.252075 kJ/mol
    CrossStacking -135.335464 kJ/mol
    Exclusion 0.11901 kJ/mol
    Electrostatics 43.296059 kJ/mol

Please make sure that the energies obtained coincide with the energies shown here. Also you can check the energy obtained using other platforms_ by changing ``OpenCL`` to ``Reference``, ``CUDA`` or ``CPU``. 

.. testcode::
    
    #Add simulation reporters
    dcd_reporter=simtk.openmm.app.DCDReporter(f'output.dcd', 1000)
    energy_reporter=simtk.openmm.app.StateDataReporter(sys.stdout, 1000, step=True,time=True,
                                                       potentialEnergy=True, temperature=True)
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(energy_reporter)

    #Run simulation
    simulation.step(10000)

.. testoutput::
    
    #"Step","Time (ps)","Potential Energy (kJ/mole)","Temperature (K)"
    1000,2.0000000000000013,-1651.121826171875,304.6066812070446
    2000,3.999999999999781,-1646.61328125,309.5945230237376
    3000,5.999999999999561,-1661.6788330078125,318.46432160647703
    4000,7.999999999999341,-1676.956298828125,268.04874840144447
    5000,10.000000000000009,-1629.8892822265625,271.8654648104738
    6000,12.000000000000677,-1622.474853515625,312.1083112301662
    7000,14.000000000001345,-1704.033203125,283.5259033832464
    8000,16.00000000000201,-1608.751708984375,281.82371603990293
    9000,18.000000000000902,-1623.486572265625,283.86225823944585
    10000,19.999999999999794,-1671.9105224609375,296.18167366285144

The forces and system in open3SPN can be treated as forces and systems from openmm. Please refer to openmm_ documentation to learn more about running the simulations or adding forces.

.. warning::

   This module only handles three-sided polygons;
   five-sided figures are right out.
   
   
DNA from atomistic PDB
-----------------------


Protein DNA system
------------------


DNA from sequence
------------------

.. _openmm: http://docs.openmm.org/latest/api-python/index.html
.. _platforms: http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.simulation.Simulation.html#simtk.openmm.app.simulation.Simulation.platform
