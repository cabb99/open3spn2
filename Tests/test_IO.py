from open3SPN2 import *

# Unit testing
def test_DNA_from_pdb():
    """ Test correct DNA initialization from PDB"""
    mol = DNA.fromPDB("Tests/1svc/1svc.pdb", template_from_X3DNA=False)


def test_DNA_from_gro():
    """ Test correct DNA initialization from gromacs files"""
    pass


def test_DNA_from_seq():
    """ Test correct DNA initialization from sequence files"""
    #return True #Needs X3DNA
    seq = 'ATACAAAGGTGCGAGGTTTCTATGCTCCCACG'
    dna = DNA.fromSequence(seq, dna_type='B_curved')

    # Compute the topology for the DNA structure.
    # Since the dna was generated from the sequence using X3DNA,
    # it is not necesary to recompute the geometry.

    dna.computeTopology(template_from_X3DNA=True)

    # Create the system.
    # To set periodic boundary conditions (periodicBox=[50,50,50]).
    # The periodic box size is in nanometers.
    dna.periodic = False
    s = System(dna, periodicBox=None)

    # Add 3SPN2 forces
    s.add3SPN2forces(verbose=True)

    # Initialize Molecular Dynamics simulations
    s.initializeMD(temperature=300 * openmm.unit.kelvin, platform_name='CPU')
    simulation = s.simulation

    # Set initial positions
    simulation.context.setPositions(s.coord.getPositions())

    energy_unit = openmm.unit.kilojoule_per_mole
    # Total energy
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(energy_unit)
    print('TotalEnergy', round(energy, 6), energy_unit.get_symbol())

    # Detailed energy
    energies = {}
    for force_name, force in s.forces.items():
        group = force.getForceGroup()
        state = simulation.context.getState(getEnergy=True, groups=2 ** group)
        energies[force_name] = state.getPotentialEnergy().value_in_unit(energy_unit)

    for force_name in s.forces.keys():
        print(force_name, round(energies[force_name], 6), energy_unit.get_symbol())


def test_DNA_from_xyz():
    """Tests the correct parsing from an xyz file"""
    mol = DNA.fromXYZ('Tests/adna/in00_conf.xyz', template_from_X3DNA=False)
    assert mol.atoms.at[8, 'name'] == 'P'
    assert round(mol.atoms.at[188, 'y'], 6) == -8.779343


# Functional testing

def test_parse_xyz():
    """Tests the example trajectory parsing"""
    xyz_data = parse_xyz('Tests/adna/traj.xyz')
    assert xyz_data.at[1, 'name'] == 7
    assert xyz_data.at[1, 'x'] == 4.34621


def test_parse_log():
    """Tests the example log parsing"""
    log_data = parse_log('Tests/adna/sim.log')
    assert log_data.at[1, 'Step'] == 2000
    assert log_data.at[1, 'eexcl'] == 0.45734636


def test_transform_zero():
    transform = Transform()
    
    np.testing.assert_almost_equal(transform.full_rotation, np.eye(3), decimal=7)
    np.testing.assert_almost_equal(transform.half_rotation, np.eye(3), decimal=7)
    np.testing.assert_almost_equal(transform.full_displacement, [0, 0, 0], decimal=7)
    np.testing.assert_almost_equal(transform.half_displacement, [0, 0, 0], decimal=7)

def test_transform_known():
    transform = Transform(0.07, -0.19, 0.07, 
                          1.8/180*np.pi, -15.0/180*np.pi, 1.5/180*np.pi)
    
    expected_full_rotation = np.array([[0.965592, -0.029813, -0.258348], 
                                       [0.021636,  0.999173, -0.034438], 
                                       [0.259161,  0.027663,  0.965438]])
    expected_half_rotation = np.array([[0.991374, -0.014114, -0.130305],
                                       [0.011951,  0.999778, -0.017370],
                                       [0.130521,  0.015662,  0.991322]])
    expected_full_displacement = [0.062957, -0.190337, 0.075553]
    expected_half_displacement = [0.031478, -0.095169, 0.037777]
    
    np.testing.assert_almost_equal(transform.full_rotation, expected_full_rotation, decimal=6)
    np.testing.assert_almost_equal(transform.half_rotation, expected_half_rotation, decimal=6)
    np.testing.assert_almost_equal(transform.full_displacement, expected_full_displacement, decimal=6)
    np.testing.assert_almost_equal(transform.half_displacement, expected_half_displacement, decimal=6)