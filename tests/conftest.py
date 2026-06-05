import numpy as np
import pytest
import openmm
import openmm.unit as unit


def pytest_addoption(parser):
    parser.addoption("--skip-platform", action="store", default=None, help="Skip tests for a specific platform")

def pytest_configure(config):
    config.addinivalue_line("markers", "skip_platform: mark test to be skipped for a specific platform")


# --------------------------------------------------------------------------- #
# Real-structure fixtures for the opt-in bias / collective-variable tests.
#
# Each fixture exposes a ForceHarness that evaluates one force at a time on a
# fixed real structure, reusing a single openmm.System (only a lightweight
# Context is rebuilt per call) so the suite does not create many systems.
# --------------------------------------------------------------------------- #

_PLATFORM = openmm.Platform.getPlatformByName('Reference')
_DNA_SEQUENCE = 'ATACAAAGGTGCGAGGTTTCTATGCTCCCACG'
_CLEAN_PDB = 'examples/Protein_DNA/clean.pdb'
_PROTEIN_SEQ = 'examples/Protein_DNA/protein.seq'


class ForceHarness:
    """Evaluate a single opt-in force on a fixed structure.

    ``dna.atoms`` row order matches the system particle order, so atom indices in
    ``dna.atoms`` are valid particle indices for the force and for the centroid
    helpers below.
    """

    def __init__(self, dna, raw_system, positions, protein=None):
        self.dna = dna
        self.protein = protein
        self._raw = raw_system
        self.positions = positions
        self.posnm = np.array(positions.value_in_unit(unit.nanometer))
        self.masses = np.array([raw_system.getParticleMass(i).value_in_unit(unit.dalton)
                                for i in range(raw_system.getNumParticles())])

    def _context(self, force_obj):
        raw = self._raw
        while raw.getNumForces() > 0:
            raw.removeForce(raw.getNumForces() - 1)
        raw.addForce(force_obj.force)
        integrator = openmm.VerletIntegrator(1 * unit.femtosecond)
        context = openmm.Context(raw, integrator, _PLATFORM)
        context.setPositions(self.positions)
        return context

    def energy(self, force_obj, group):
        """Potential energy (kJ/mol) of ``force_obj`` read from its force group."""
        state = self._context(force_obj).getState(getEnergy=True, groups={int(group)})
        return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    def max_force_error(self, force_obj, group, h=1e-5):
        """Largest |analytic force - central finite difference| over active atoms."""
        context = self._context(force_obj)
        analytic = context.getState(getForces=True, groups={int(group)}).getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer)

        def energy_at(coords):
            context.setPositions(coords * unit.nanometer)
            return context.getState(getEnergy=True, groups={int(group)}).getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole)

        worst = 0.0
        active = np.where(np.abs(analytic).sum(axis=1) > 1e-9)[0]
        for i in active:
            for dim in range(3):
                up = self.posnm.copy(); up[i, dim] += h
                dn = self.posnm.copy(); dn[i, dim] -= h
                numeric = -(energy_at(up) - energy_at(dn)) / (2 * h)
                worst = max(worst, abs(numeric - analytic[i, dim]))
        return worst

    def mass_centroid(self, indices):
        weights = self.masses[indices]
        return (weights[:, None] * self.posnm[indices]).sum(axis=0) / weights.sum()

    def geometric_centroid(self, indices):
        return self.posnm[indices].mean(axis=0)

    def atom_indices(self, resSeq=None, chainID=None, name=None, nonzero_mass=False):
        atoms = self.dna.atoms
        mask = np.ones(len(atoms), dtype=bool)
        if resSeq is not None:
            mask &= (atoms['resSeq'].values == resSeq)
        if chainID is not None:
            mask &= (atoms['chainID'].values == chainID)
        if name is not None:
            mask &= (atoms['name'].values == name)
        indices = [int(i) for i in atoms.index[mask]]
        if nonzero_mass:
            indices = [i for i in indices if self.masses[i] > 0]
        return indices


@pytest.fixture(scope='module')
def dna_harness():
    """DNA-only real structure built from a sequence (no AWSEM needed)."""
    from open3SPN2 import DNA, System
    dna = DNA.fromSequence(_DNA_SEQUENCE)
    dna.periodic = False
    system = System(dna)
    return ForceHarness(dna, system._wrapped_system, system.coord.getPositions())


@pytest.fixture(scope='module')
def protein_dna_harness():
    """Real merged protein-DNA structure from examples/Protein_DNA/clean.pdb.

    Skipped when openawsem is not installed (it is needed to build the protein
    object and the merged force field). Chains 1,2 are DNA; 3,4 are protein.
    """
    openawsem = pytest.importorskip('openawsem')
    import openmm.app
    import open3SPN2

    dna = open3SPN2.DNA.fromCoarsePDB(_CLEAN_PDB)
    dna.periodic = False
    with open(_PROTEIN_SEQ) as seq_file:
        sequence = seq_file.readlines()[0].strip()
    protein = openawsem.Protein.fromCoarsePDB(_CLEAN_PDB, sequence=sequence)
    protein.periodic = False

    pdb = openmm.app.PDBFile(_CLEAN_PDB)
    forcefield = openmm.app.ForceField(openawsem.xml, open3SPN2.xml)
    system = forcefield.createSystem(pdb.topology)
    return ForceHarness(dna, system, pdb.getPositions(), protein=protein)
