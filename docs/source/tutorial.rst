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

First Generate a Coarse Grained model

.. code:: ipython3

    import open3SPN2
    import ffAWSEM

.. code:: ipython3

    #Fix the system (adds missing atoms)
    fix=open3SPN2.fixPDB("1lmb.pdb")

.. code:: ipython3

    #Create a table containing both the proteins and the DNA
    complex_table=open3SPN2.pdb2table(fix)

.. code:: ipython3

    #Coarse Grain the system
    dna_atoms=open3SPN2.DNA.CoarseGrain(complex_table)
    protein_atoms=ffAWSEM.Protein.CoarseGrain(complex_table)

.. code:: ipython3

    #Merge the models
    import pandas
    Coarse=pandas.concat([protein_atoms,dna_atoms],sort=False)
    Coarse.index=range(len(Coarse))
    Coarse.serial=list(Coarse.index)

.. code:: ipython3

    #Save the protein_sequence
    from Bio.PDB.Polypeptide import three_to_one
    _AWSEMresidues=['IPR','IGL','NGP']
    protein_data=Coarse[Coarse.resname.isin(_AWSEMresidues)].copy()
    resix = (protein_data.chainID + '_' + protein_data.resSeq.astype(str))
    res_unique = resix.unique()
    protein_data['resID'] = resix.replace(dict(zip(res_unique, range(len(res_unique)))))
    protein_sequence=[r.iloc[0]['real_resname'] for i, r in protein_data.groupby('resID')]      
    protein_sequence_one = [three_to_one(a) for a in protein_sequence]
    
    with open('protein.seq','w+') as ps:
        ps.write(''.join(protein_sequence_one))

.. code:: ipython3

    # Create a merged PDB
    def writePDB(atoms,pdb_file):
        with open(pdb_file, 'w+') as pdb:
            for i, atom in atoms.iterrows():
                pdb_line = f'{atom.recname:<6}{atom.serial:>5} {atom["name"]:^4}{atom.altLoc:1}'+\
                           f'{atom.resname:<3} {atom.chainID:1}{atom.resSeq:>4}{atom.iCode:1}   '+\
                           f'{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}' +\
                           f'{atom.occupancy:>6.2f}{atom.occupancy:>6.2f}'+' ' * 10 +\
                           f'{atom.element:>2}{atom.charge:>2}'
                assert len(pdb_line) == 80, f'An item in the atom table is longer than expected ({len(pdb_line)})\n{pdb_line}'
                pdb.write(pdb_line + '\n')
    writePDB(Coarse,'clean.pdb')

Then generate the system.

.. code:: ipython3

    #Create the merged system
    import simtk.openmm
    pdb=simtk.openmm.app.PDBFile('clean.pdb')
    top=pdb.topology
    coord=pdb.positions
    forcefield=simtk.openmm.app.ForceField(ffAWSEM.xml,open3SPN2.xml)
    s=forcefield.createSystem(top)

Then add the forces

.. code:: ipython3

    dna=open3SPN2.DNA.fromCoarsePDB('clean.pdb')
    with open('protein.seq') as ps:
        protein_sequence_one=ps.readlines()[0]
    protein=ffAWSEM.Protein.fromCoarsePDB('clean.pdb',sequence=protein_sequence_one)
    dna.periodic=False
    protein.periodic=False

.. code:: ipython3

    ffAWSEM.copy_parameter_files()

.. code:: ipython3

    #Clear Forces from the system
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
        forces.update({force_name: force})
        
    #Add AWSEM forces
    openAWSEMforces = dict(Connectivity=ffAWSEM.functionTerms.basicTerms.con_term,
                           Chain=ffAWSEM.functionTerms.basicTerms.chain_term,
                           Chi=ffAWSEM.functionTerms.basicTerms.chi_term,
                           Excl=ffAWSEM.functionTerms.basicTerms.excl_term_v2,
                           rama=ffAWSEM.functionTerms.basicTerms.rama_term,
                           rama_pro=ffAWSEM.functionTerms.basicTerms.rama_proline_term,
                           #rama_ss=ffAWSEM.functionTerms.basicTerms.rama_ssweight_term,
                           contact=ffAWSEM.functionTerms.contactTerms.contact_term,
                           beta1 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_1,
                           beta2 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_2,
                           beta3 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_3,
                           pap1 = ffAWSEM.functionTerms.hydrogenBondTerms.pap_term_1,
                           pap2 = ffAWSEM.functionTerms.hydrogenBondTerms.pap_term_2,
                          )
    protein.setup_virtual_sites(s)
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
    
    #Add DNA-protein interaction forces
    for force_name in open3SPN2.protein_dna_forces:
        print(force_name)
        force = open3SPN2.protein_dna_forces[force_name](dna,protein)
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
    Connectivity
    Chain
    Chi
    Excl
    rama
    rama_pro
    contact
    Number of atom:  1171 Number of residue:  179
    Contact cutoff  1.0 nm
    NonbondedMethod:  1
    0
    639
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
    ExclusionProteinDNA
    ElectrostaticsProteinDNA


Then set-up the simulation

.. code:: ipython3

    import numpy as np
    temperature=300 * simtk.openmm.unit.kelvin
    #platform_name='CUDA'
    
    platform_name='OpenCL'
    
    integrator = simtk.openmm.LangevinIntegrator(temperature, 1 / simtk.openmm.unit.picosecond, 2 * simtk.openmm.unit.femtoseconds)
    platform = simtk.openmm.Platform.getPlatformByName(platform_name)
    simulation = simtk.openmm.app.Simulation(top,s, integrator, platform)
    simulation.context.setPositions(coord)
    energy_unit=simtk.openmm.unit.kilojoule_per_mole
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(energy_unit)
    print(energy)


.. parsed-literal::

    -1406.7860107421875


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

    TotalEnergy -1406.786011 kJ/mol
    Bond 0.0 kJ/mol
    Angle 0.0 kJ/mol
    Stacking 203.56601 kJ/mol
    Dihedral -503.999969 kJ/mol
    BasePair -284.232208 kJ/mol
    CrossStacking -47.58614 kJ/mol
    Exclusion 23.991552 kJ/mol
    Electrostatics 23.268293 kJ/mol
    Connectivity 1899.296875 kJ/mol
    Chain 1899.296875 kJ/mol
    Chi 1899.296875 kJ/mol
    Excl 1899.296875 kJ/mol
    rama -1363.522705 kJ/mol
    rama_pro -1363.522705 kJ/mol
    contact -1041.547729 kJ/mol
    beta1 -601.593384 kJ/mol
    beta2 -601.593384 kJ/mol
    beta3 -601.593384 kJ/mol
    pap1 0.0 kJ/mol
    pap2 0.0 kJ/mol
    ExclusionProteinDNA 296.033478 kJ/mol
    ElectrostaticsProteinDNA -10.459808 kJ/mol


.. code:: ipython3

    #Add simulation reporters
    import sys
    dcd_reporter=simtk.openmm.app.DCDReporter(f'output.dcd', 1000)
    energy_reporter=simtk.openmm.app.StateDataReporter(sys.stdout, 1000, step=True,time=True,
                                                       potentialEnergy=True, temperature=True)
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(energy_reporter)

.. code:: ipython3

    #Run simulation
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(10000)


.. parsed-literal::

    #"Step","Time (ps)","Potential Energy (kJ/mole)","Temperature (K)"
    1000,2.0000000000000013,-2821.62255859375,285.56926619310553
    2000,3.999999999999781,-2839.069580078125,315.4594998778808
    3000,5.999999999999561,-2805.5634765625,323.20843224380087
    4000,7.999999999999341,-2801.427734375,314.2296292320282
    5000,10.000000000000009,-2491.56884765625,305.8918233488531
    6000,12.000000000000677,-2708.4228515625,293.2187270882386
    7000,14.000000000001345,-2837.1474609375,296.9420474031121
    8000,16.00000000000201,-2738.265380859375,306.50966191307515
    9000,18.000000000000902,-2657.7529296875,312.0026886071036
    10000,19.999999999999794,-2780.637451171875,318.9033562086061


.. _openmm: http://docs.openmm.org/latest/api-python/index.html
.. _platforms: http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.simulation.Simulation.html#simtk.openmm.app.simulation.Simulation.platform
