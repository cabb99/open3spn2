Tutorial
================================================

    .. index:: single: Feynman
    
    *“What I cannot create, I do not understand.”* — Feynman
    
DNA from sequence
------------------

This section allows you to quickly simulate a piece of DNA. 

First import the open3SPN2 module.

.. testcode::
    
    import open3SPN2

If you obtain a ``ModuleNotFoundError`` error then you may need to check that ``open3SPN2`` is in the installation path as detailed in :doc:`installation`.

After importing the module, follow the next steps to create a DNA system and add the 3SPN2 forces.
Make sure you have installed X3DNA before this step.

.. testsetup:: *

    import open3SPN2

.. testcode::
    
    # Initialize the DNA from a sequence.
    # DNA type can be changed to 'A' or 'B'
    
    seq='ATACAAAGGTGCGAGGTTTCTATGCTCCCACG'
    dna=open3SPN2.DNA.fromSequence(seq,dna_type='B_curved')
    
    # Compute the topology for the DNA structure.
    # Since the dna was generated from the sequence using X3DNA,
    # it is not necesary to recompute the geometry.
    
    dna.computeTopology(template_from_X3DNA=False)
    
    # Create the system.
    # To set periodic boundary conditions (periodicBox=[50,50,50]).
    # The periodic box size is in nanometers.
    dna.periodic=False
    s=open3SPN2.System(dna, periodicBox=None)
    
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

.. DNA from atomistic PDB
.. ----------------------


Protein DNA system
------------------
You can find this example on the `examples/Protein_DNA <https://github.com/cabb99/open3spn2/tree/master/examples/Protein_DNA>`_ folder. For this example you need to download the structure of the Lambda repressor-operator complex `1mlb.pdb <https://www.rcsb.org/structure/1LMB>`_. You will also need to have installed the `openAWSEM <https://github.com/npschafer/openawsem>`_ library.

.. code:: ipython3

    # If you want to specify the package address
    # you can add them to the PYTHONPATH environment variable.
    # Also you can add them on the run time uncommenting the lines below
    # import sys
    # open3SPN2_HOME = '/Users/weilu/open3spn2/'
    # openAWSEM_HOME = '/Users/weilu/openmmawsem/'
    # sys.path.insert(0,open3SPN2_HOME)
    # sys.path.insert(0,openAWSEM_HOME)

.. code:: ipython3

    #Import openAWSEM, open3SPN2 and other libraries
    import open3SPN2
    import ffAWSEM

    import pandas
    import numpy as np
    import simtk.openmm

    from functools import partial
    import sys

.. code:: ipython3

    #Fix the system (adds missing atoms)
    fix=open3SPN2.fixPDB("1lmb.pdb")

.. code:: ipython3

    #Create a table containing both the proteins and the DNA
    complex_table=open3SPN2.pdb2table(fix)

    # Create a single memory file
    ffAWSEM.create_single_memory(fix)

.. code:: ipython3

    #Generate a coarse-grained model of the DNA molecules
    dna_atoms=open3SPN2.DNA.CoarseGrain(complex_table)

    #Generate a coarse-grained model of the Protein molecules
    protein_atoms=ffAWSEM.Protein.CoarseGrain(complex_table)

.. code:: ipython3

    #Merge the models
    Coarse=pandas.concat([protein_atoms,dna_atoms],sort=False)
    Coarse.index=range(len(Coarse))
    Coarse['serial']=list(Coarse.index)

.. code:: ipython3

    #Save the protein_sequence
    ffAWSEM.save_protein_sequence(Coarse,sequence_file='protein.seq')

.. code:: ipython3

    # Create a merged PDB
    ffAWSEM.writePDB(Coarse,'clean.pdb')

.. code:: ipython3

    #Create the merged system

    pdb=simtk.openmm.app.PDBFile('clean.pdb')
    top=pdb.topology
    coord=pdb.positions
    forcefield=simtk.openmm.app.ForceField(ffAWSEM.xml,open3SPN2.xml)
    s=forcefield.createSystem(top)

.. code:: ipython3

    dna=open3SPN2.DNA.fromCoarsePDB('clean.pdb')
    with open('protein.seq') as ps:
        protein_sequence_one=ps.readlines()[0]
    protein=ffAWSEM.Protein.fromCoarsePDB('clean.pdb',sequence=protein_sequence_one)
    dna.periodic=False
    protein.periodic=False

.. code:: ipython3

    #Create the DNA and Protein Objects
    dna=open3SPN2.DNA.fromCoarsePDB('clean.pdb')
    with open('protein.seq') as ps:
        protein_sequence_one=ps.readlines()[0]
    protein=ffAWSEM.Protein.fromCoarsePDB('clean.pdb',sequence=protein_sequence_one)
    dna.periodic=False
    protein.periodic=False
    
    #Copy the AWSEM parameter files
    ffAWSEM.copy_parameter_files()

.. code:: ipython3

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
    
.. code:: ipython3

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
    openAWSEMforces = dict(Connectivity=ffAWSEM.functionTerms.basicTerms.con_term,
                           Chain=ffAWSEM.functionTerms.basicTerms.chain_term,
                           Chi=ffAWSEM.functionTerms.basicTerms.chi_term,
                           Excl=ffAWSEM.functionTerms.basicTerms.excl_term,
                           rama=ffAWSEM.functionTerms.basicTerms.rama_term,
                           rama_pro=ffAWSEM.functionTerms.basicTerms.rama_proline_term,
                           contact=ffAWSEM.functionTerms.contactTerms.contact_term,
                           frag  = partial(ffAWSEM.functionTerms.templateTerms.fragment_memory_term, 
                                           frag_file_list_file="./single_frags.mem", 
                                           npy_frag_table="./single_frags.npy", 
                                           UseSavedFragTable=False, 
                                           k_fm=0.04184/3),
                           beta1 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_1,
                           beta2 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_2,
                           beta3 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_3,
                           pap1 = ffAWSEM.functionTerms.hydrogenBondTerms.pap_term_1,
                           pap2 = ffAWSEM.functionTerms.hydrogenBondTerms.pap_term_2,
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
            force = openAWSEMforces[force_name](protein, withExclusion=False,periodic=False)
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


.. parsed-literal::

    Bond
    Angle
    Stacking
    Dihedral
    BasePair
    CrossStacking
    Exclusion
    Electrostatics
    ExclusionProteinDNA
    ElectrostaticsProteinDNA
    Connectivity
    Chain
    Chi
    Excl
    1205
    1844
    rama
    rama_pro
    contact
    Number of atom:  1171 Number of residue:  179
    Contact cutoff  1.0 nm
    NonbondedMethod:  1
    0
    639
    frag
    Loading Fragment files(Gro files)
    Saving fragment table as npy file to speed up future calculation.
    All gro files information have been stored in the ./single_frags.npy.             
    You might want to set the 'UseSavedFragTable'=True to speed up the loading next time.             
    But be sure to remove the .npy file if you modify the .mem file. otherwise it will keep using the old frag memeory.
    beta1
    beta_1 term ON
    beta2
    beta_2 term ON
    beta3
    beta_3 term ON
    pap1
    pap_1 term ON
    No ssweight given, assume all zero
    pap2
    pap_2 term ON
    No ssweight given, assume all zero


.. code:: ipython3

    # Set up the simulation
    temperature=300 * simtk.openmm.unit.kelvin
    platform_name='OpenCL' #'Reference','CPU','CUDA', 'OpenCL'

    integrator = simtk.openmm.LangevinIntegrator(temperature, 1 / simtk.openmm.unit.picosecond, 2 * simtk.openmm.unit.femtoseconds)
    platform = simtk.openmm.Platform.getPlatformByName(platform_name)
    simulation = simtk.openmm.app.Simulation(top,s, integrator, platform)
    simulation.context.setPositions(coord)
    energy_unit=simtk.openmm.unit.kilojoule_per_mole
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(energy_unit)
    print(energy)


.. parsed-literal::

    -899.144763339983

.. code:: ipython3

    #Obtain total energy

    energy_unit=simtk.openmm.unit.kilojoule_per_mole
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(energy_unit)
    print('TotalEnergy',round(energy,6),energy_unit.get_symbol())

    #Obtain detailed energy

    energies = {}
    for force_name, force in forces.items():
        group=force.getForceGroup()
        state = simulation.context.getState(getEnergy=True, groups=2**group)
        energies[force_name] =state.getPotentialEnergy().value_in_unit(energy_unit)

    for force_name in forces.keys():
        print(force_name, round(energies[force_name],6),energy_unit.get_symbol())


.. parsed-literal::

    TotalEnergy -899.147583 kJ/mol
    Bond 327.558105 kJ/mol
    Angle 973.859009 kJ/mol
    Stacking 203.565979 kJ/mol
    Dihedral -385.277161 kJ/mol
    BasePair -284.232208 kJ/mol
    CrossStacking -47.586143 kJ/mol
    Exclusion 23.991552 kJ/mol
    Electrostatics 23.268274 kJ/mol
    ExclusionProteinDNA 296.033478 kJ/mol
    ElectrostaticsProteinDNA -10.459805 kJ/mol
    Connectivity 1899.296875 kJ/mol
    Chain 1899.296875 kJ/mol
    Chi 1899.296875 kJ/mol
    Excl 1899.296875 kJ/mol
    rama -1363.522705 kJ/mol
    rama_pro -1363.522705 kJ/mol
    contact -1041.547607 kJ/mol
    frag -1213.29834 kJ/mol
    beta1 -300.796692 kJ/mol
    beta2 -300.796692 kJ/mol
    beta3 -300.796692 kJ/mol
    pap1 0.0 kJ/mol
    pap2 0.0 kJ/mol


.. code:: ipython3

    #Add simulation reporters
    dcd_reporter=simtk.openmm.app.DCDReporter(f'output.dcd', 10000)
    energy_reporter=simtk.openmm.app.StateDataReporter(sys.stdout, 10000, step=True,time=True,
                                                       potentialEnergy=True, temperature=True)
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(energy_reporter)

.. code:: ipython3

    #Run simulation
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(100000)


.. parsed-literal::

    #"Step","Time (ps)","Potential Energy (kJ/mole)","Temperature (K)"
    1000,2.0000000000000013,-3507.7216796875,291.6758662031427
    2000,3.999999999999781,-3233.395751953125,299.73657773689456
    3000,5.999999999999561,-3472.61083984375,314.04882500167656
    4000,7.999999999999341,-3199.33251953125,309.8776608226663
    5000,10.000000000000009,-3275.126220703125,321.27214200690065
    6000,12.000000000000677,-3167.643798828125,295.1803751740272
    7000,14.000000000001345,-3187.2998046875,310.9197062404284
    8000,16.00000000000201,-3296.966064453125,305.38198987465074
    9000,18.000000000000902,-3182.54443359375,316.07187422798313
    10000,19.999999999999794,-3229.941650390625,309.6002450725328
        
.. code:: ipython3

    #Get the detailed energy after the simulation
    energies = {}
    for force_name, force in forces.items():
        group=force.getForceGroup()
        state = simulation.context.getState(getEnergy=True, groups=2**group)
        energies[force_name] =state.getPotentialEnergy().value_in_unit(energy_unit)

    for force_name in forces.keys():
        print(force_name, round(energies[force_name],6),energy_unit.get_symbol())


.. parsed-literal::

    Bond 102.787193 kJ/mol
    Angle 169.343231 kJ/mol
    Stacking -411.438232 kJ/mol
    Dihedral -452.744202 kJ/mol
    BasePair -242.723328 kJ/mol
    CrossStacking -48.275833 kJ/mol
    Exclusion 1.823019 kJ/mol
    Electrostatics 23.37781 kJ/mol
    ExclusionProteinDNA -6.992864 kJ/mol
    ElectrostaticsProteinDNA -10.225787 kJ/mol
    Connectivity 1487.126465 kJ/mol
    Chain 1487.126587 kJ/mol
    Chi 1487.126465 kJ/mol
    Excl 1487.126587 kJ/mol
    rama -1456.633789 kJ/mol
    rama_pro -1456.633789 kJ/mol
    contact -1326.107178 kJ/mol
    frag -922.991699 kJ/mol
    beta1 -136.266449 kJ/mol
    beta2 -136.266449 kJ/mol
    beta3 -136.266449 kJ/mol
    pap1 0.0 kJ/mol
    pap2 0.0 kJ/mol

.. _openmm: http://docs.openmm.org/latest/api-python/index.html
.. _platforms: http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.simulation.Simulation.html#simtk.openmm.app.simulation.Simulation.platform
