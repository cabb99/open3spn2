"""Tests for the sequence separation masks of the nonbonded and hydrogen bond forces.

OpenMM requires every CustomNonbondedForce of a system to carry the same exclusions, so the forces
mask pairs inside their energy expressions instead of adding exclusions. The behavioural
requirements tested here:

* No force adds pair exclusions, so the forces compose with each other and with other force fields
  on the platforms that enforce identical exclusion lists.
* Pairs closer than ``min_seq_sep`` residues along a chain are masked, and pairs on different
  chains are never masked, for both the nonbonded and the hydrogen bond forces.
* The mask follows the residue numbering rather than the order of the particles, so residues that
  are not consecutive in a chain are not masked as neighbors.
"""
import itertools
import pytest
from pathlib import Path

from open3SPN2 import DNA, forces
from open3SPN2.force.dna import chained_sequence_index

test_path = Path(__file__).parent

NONBONDED_FORCES = ['Exclusion', 'Electrostatics']
HBOND_FORCES = ['BasePair', 'CrossStacking']


@pytest.fixture(scope='module')
def dna():
    """A two chain B-DNA duplex."""
    return DNA.fromXYZ(str(test_path / 'bdna' / 'in00_conf.xyz'), 'B', template_from_X3DNA=False)


def _openmm_forces(force):
    """The OpenMM forces of a wrapper, which may hold several of them."""
    if hasattr(force, 'forces'):
        return list(force.forces.values())
    if hasattr(force, 'crossStackingForces'):
        return [f for pair in force.crossStackingForces.values() for f in pair]
    return [force.force]


def _seqids(force, kind='particle'):
    """Sequence index of every particle, donor or acceptor of an OpenMM force."""
    if kind == 'particle':
        names = [force.getPerParticleParameterName(i)
                 for i in range(force.getNumPerParticleParameters())]
        return [force.getParticleParameters(i)[names.index('seqid')]
                for i in range(force.getNumParticles())]
    if kind == 'donor':
        names = [force.getPerDonorParameterName(i) for i in range(force.getNumPerDonorParameters())]
        return [force.getDonorParameters(i)[3][names.index('seqid_d')]
                for i in range(force.getNumDonors())]
    names = [force.getPerAcceptorParameterName(i)
             for i in range(force.getNumPerAcceptorParameters())]
    return [force.getAcceptorParameters(i)[3][names.index('seqid_a')]
            for i in range(force.getNumAcceptors())]


@pytest.mark.parametrize('force_name', NONBONDED_FORCES)
def test_nonbonded_forces_have_no_exclusions(dna, force_name):
    """Every CustomNonbondedForce must share the same (empty) exclusion list."""
    for force in _openmm_forces(forces[force_name](dna)):
        assert force.getNumExclusions() == 0


@pytest.mark.parametrize('force_name', NONBONDED_FORCES)
def test_nonbonded_mask_follows_sequence_separation(dna, force_name):
    """Particles closer than min_seq_sep residues in a chain are masked, other pairs are not."""
    wrapper = forces[force_name](dna)
    min_seq_sep = wrapper.min_seq_sep
    atoms = dna.atoms
    for force in _openmm_forces(wrapper):
        seqid = _seqids(force, 'particle')
        assert len(seqid) == len(atoms)
        for i, j in itertools.combinations(range(len(atoms)), 2):
            same_chain = atoms['chainID'][i] == atoms['chainID'][j]
            close = abs(atoms['resSeq'][i] - atoms['resSeq'][j]) < min_seq_sep
            assert (abs(seqid[i] - seqid[j]) < min_seq_sep) == (same_chain and close)


@pytest.mark.parametrize('force_name', HBOND_FORCES)
def test_hbond_mask_follows_sequence_separation(dna, force_name):
    """Donors and acceptors closer than min_seq_sep residues in a chain are masked.

    Their sequence indices come from different parameter lists, so a donor and an acceptor of the
    same residue must still get the same index."""
    wrapper = forces[force_name](dna)
    min_seq_sep = wrapper.min_seq_sep
    seqid = chained_sequence_index(dna.atoms)
    for force in _openmm_forces(wrapper):
        donors = _seqids(force, 'donor')
        acceptors = _seqids(force, 'acceptor')
        assert set(donors) <= set(seqid) and set(acceptors) <= set(seqid)
        masked = [(d, a) for d in donors for a in acceptors if abs(d - a) < min_seq_sep]
        assert masked, 'the mask never applies, it can not be tested'
        for d, a in masked:
            residues = dna.atoms[seqid.isin([d, a])]
            assert residues['chainID'].nunique() == 1
            assert residues['resSeq'].max() - residues['resSeq'].min() < min_seq_sep


def test_sequence_index_separates_chains_and_keeps_numbering_gaps(dna):
    """Chains stay further apart than any pair inside a chain and numbering gaps are preserved."""
    atoms = dna.atoms.copy()
    # Renumber the second half of the first chain to leave a gap in the residue numbering
    chain = atoms['chainID'] == atoms['chainID'].iloc[0]
    gapped = chain & (atoms['resSeq'] > atoms.loc[chain, 'resSeq'].median())
    atoms.loc[gapped, 'resSeq'] += 10
    seqid = chained_sequence_index(atoms)

    separation = abs(seqid.values[:, None] - seqid.values[None, :])
    different_chain = atoms['chainID'].values[:, None] != atoms['chainID'].values[None, :]
    assert separation[different_chain].min() > separation[~different_chain].max()

    assert seqid[gapped].min() - seqid[chain & ~gapped].max() == 11
