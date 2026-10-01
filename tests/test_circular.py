"""Tests for circular DNA (plasmids).

* The bonded interactions of a circular chain wrap around its end: every strand of a closed duplex
  has as many bonds as beads, and the closing bonds are near their equilibrium length.
* The sequence separation masks of the nonbonded and hydrogen bond forces count residues around
  the circle.
* The bent circle has the requested linking number, and minimizing it keeps it closed.
* Circular structures can be loaded from atomistic and coarse grained pdb files, and a straight
  duplex can be bent into a circle.
"""
import itertools
import numpy as np
import openmm
import openmm.unit as unit
import pytest

from open3SPN2 import DNA, System, forces
from open3SPN2.circular import linking_number, residue_wrapper
from open3SPN2.force.dna import chained_sequence_index, circular_sequence_length

SEQUENCE = 'ATACAAAGGTGCGAGGTTTCTATGCTCCCACGGTCAATGCATTGAGCCATTAGGCTACGTTAGCATCGATT'  # 72 bp


@pytest.fixture(scope='module')
def workdir(tmp_path_factory):
    return tmp_path_factory.mktemp('circular')


@pytest.fixture(scope='module')
def plasmid(workdir):
    return DNA.fromSequence(SEQUENCE, circular=True, output_pdb=str(workdir / 'plasmid.pdb'),
                            atomistic_template=str(workdir / 'plasmid_atomistic.pdb'))


def _bond_deviation(dna, xyz=None):
    xyz = dna.atoms[['x', 'y', 'z']].values / 10 if xyz is None else xyz
    r = np.linalg.norm(xyz[dna.bonds['aai']] - xyz[dna.bonds['aaj']], axis=1)
    return np.abs(r - dna.bonds['r0'].values)


def _minimize(dna, iterations=500):
    system = System(dna)
    system.add3SPN2forces()
    system.initializeMD(platform_name='Reference')
    context = system.simulation.context
    initial = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    system.simulation.minimizeEnergy(maxIterations=iterations)
    state = context.getState(getEnergy=True, getPositions=True)
    final = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return initial, final, state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)


def test_circular_topology(plasmid):
    n = len(SEQUENCE)
    assert plasmid.circular_chains == ['A', 'B']
    assert len(plasmid.atoms) == 6 * n  # every residue has a phosphate
    assert len(plasmid.bonds) == 6 * n
    assert len(plasmid.angles) == 8 * n
    assert len(plasmid.dihedrals) == 8 * n
    assert len(plasmid.stackings) == 2 * n
    assert not plasmid.bonds['r0'].isna().any()
    assert not plasmid.angles['t0'].isna().any()
    assert not plasmid.dihedrals['t0'].isna().any()
    # The closing bonds exist and are close to equilibrium
    atoms = plasmid.atoms
    for chain in ['A', 'B']:
        last_sugar = atoms[(atoms.chainID == chain) & (atoms.resSeq == n) & (atoms.name == 'S')].index[0]
        first_phosphate = atoms[(atoms.chainID == chain) & (atoms.resSeq == 1) & (atoms.name == 'P')].index[0]
        closing = plasmid.bonds[(plasmid.bonds.aai == last_sugar) & (plasmid.bonds.aaj == first_phosphate)]
        assert len(closing) == 1
    assert _bond_deviation(plasmid).max() < 0.1  # strain of a 72 bp minicircle (radius ~4 nm)


def test_linear_topology_is_unchanged(workdir):
    dna = DNA.fromSequence(SEQUENCE, output_pdb=str(workdir / 'linear.pdb'),
                           atomistic_template=str(workdir / 'linear_atomistic.pdb'))
    n = len(SEQUENCE)
    assert len(dna.atoms) == 6 * n - 2
    assert len(dna.bonds) == 6 * n - 4


def test_closing_parameters_match_the_linear_template(plasmid, workdir):
    """The equilibrium values of the interactions around the closure come from the same template
    as those of a linear DNA with the circular sequence rotated."""
    rotated = SEQUENCE[-10:] + SEQUENCE[:-10]
    linear = DNA.fromSequence(rotated, output_pdb=str(workdir / 'rotated.pdb'),
                              atomistic_template=str(workdir / 'rotated_atomistic.pdb'))

    def table(dna, shift):
        n = len(SEQUENCE)
        atoms = dna.atoms
        label = [(c, (r - 1 + shift) % n if c == 'A' else (r - 1 - shift) % n, name)
                 for c, r, name in zip(atoms.chainID, atoms.resSeq, atoms.name)]
        return {(label[i], label[j]): r0 for i, j, r0 in zip(dna.bonds.aai, dna.bonds.aaj, dna.bonds.r0)}

    circular = table(plasmid, 0)
    rotated_bonds = table(linear, -10)
    common = set(circular) & set(rotated_bonds)
    assert len(common) == len(rotated_bonds)
    # The closure of the plasmid (residue n -> 1) is the middle of the rotated linear DNA
    n = len(SEQUENCE)
    assert (('A', n - 1, 'S'), ('A', 0, 'P')) in common
    for k in common:
        assert circular[k] == pytest.approx(rotated_bonds[k], abs=1e-6)


@pytest.mark.parametrize('delta', [0, -2, 1])
def test_linking_number(workdir, delta):
    dna = DNA.fromSequence(SEQUENCE, circular=True, delta_linking_number=delta,
                           output_pdb=str(workdir / f'lk{delta}.pdb'),
                           atomistic_template=str(workdir / f'lk{delta}_atomistic.pdb'))
    relaxed = round(len(SEQUENCE) / 10.5)
    assert dna.bend.relaxed_linking_number == relaxed
    assert dna.computeLinkingNumber() == pytest.approx(relaxed + delta, abs=1e-6)
    initial, final, xyz = _minimize(dna)
    assert final < initial
    assert _bond_deviation(dna, xyz).max() < 0.05
    assert dna.computeLinkingNumber(xyz) == pytest.approx(relaxed + delta, abs=1e-6)


def test_linking_number_of_simple_curves():
    t = np.linspace(0, 2 * np.pi, 50, endpoint=False)
    ring = np.stack([np.cos(t), np.sin(t), 0 * t], axis=1)
    linked = np.stack([1 + np.cos(t), 0 * t, np.sin(t)], axis=1)
    apart = linked + [5, 0, 0]
    assert abs(linking_number(ring, linked)) == pytest.approx(1)
    assert linking_number(ring, apart) == pytest.approx(0, abs=1e-9)


def test_wrapped_sequence_separation(plasmid):
    """Residues near the closure of a circular chain are close in sequence, residues on different
    chains are never close."""
    atoms = plasmid.atoms
    seqid = chained_sequence_index(atoms, circular_chains=plasmid.circular_chains)
    seqlen = circular_sequence_length(atoms, plasmid.circular_chains)
    n = len(SEQUENCE)
    assert set(seqlen) == {n}

    def separation(i, j):
        d = abs(seqid[i] - seqid[j])
        length = max(seqlen[i], seqlen[j])
        return min(d, length - d) if d < length else d

    wrap = residue_wrapper(atoms, plasmid.circular_chains)
    for i, j in itertools.combinations(range(0, len(atoms), 7), 2):
        same_chain = atoms.chainID[i] == atoms.chainID[j]
        d = abs(atoms.resSeq[i] - atoms.resSeq[j])
        expected = min(d, n - d) if same_chain else None
        if same_chain:
            assert separation(i, j) == expected
        else:
            assert separation(i, j) >= 3
    assert wrap('A', n + 1) == 1 and wrap('A', 0) == n


@pytest.mark.parametrize('force_name', ['Exclusion', 'Electrostatics'])
def test_nonbonded_mask_wraps_around(plasmid, force_name):
    """The nonbonded energy between the closing residues is masked, like for any neighbors."""
    wrapper = forces[force_name](plasmid)
    force = wrapper.force
    names = [force.getPerParticleParameterName(i) for i in range(force.getNumPerParticleParameters())]
    assert 'seqlen' in names
    assert force.getNumExclusions() == 0
    # Evaluate the mask expression with the parameters of the closing residues
    atoms = plasmid.atoms
    n = len(SEQUENCE)
    first = atoms[(atoms.chainID == 'A') & (atoms.resSeq == 1) & (atoms.name == 'S')].index[0]
    last = atoms[(atoms.chainID == 'A') & (atoms.resSeq == n) & (atoms.name == 'S')].index[0]
    middle = atoms[(atoms.chainID == 'A') & (atoms.resSeq == n // 2) & (atoms.name == 'S')].index[0]
    expression = force.getEnergyFunction()
    mask = expression.split('seqsep=')[1].split(';')[0]

    def evaluate(i, j):
        p1 = force.getParticleParameters(int(i))
        p2 = force.getParticleParameters(int(j))
        values = {f'{name}1': v for name, v in zip(names, p1)}
        values.update({f'{name}2': v for name, v in zip(names, p2)})
        return eval(mask, {'abs': abs, 'max': max, 'step': lambda x: float(x >= 0)}, values)

    assert evaluate(first, last) == 0
    assert evaluate(first, middle) == 1


def test_from_atomistic_and_coarse_pdb(plasmid, workdir):
    """A circular structure is recognized from atomistic and coarse grained pdb files, also when the
    phosphate of the first residue is missing."""
    n = len(SEQUENCE)
    atomistic = DNA.fromPDB(str(workdir / 'plasmid_atomistic.pdb'), circular=True,
                            output_pdb=str(workdir / 'from_atomistic.pdb'))
    assert len(atomistic.atoms) == 6 * n
    assert len(atomistic.bonds) == 6 * n
    assert _bond_deviation(atomistic).max() < 0.1
    assert not hasattr(atomistic, 'bend')  # already closed, not bent

    coarse = plasmid.atoms
    coarse = coarse[~((coarse.chainID == 'A') & (coarse.resSeq == 1) & (coarse.name == 'P'))]
    no_phosphate = DNA()
    no_phosphate.atoms = coarse.reset_index(drop=True)
    no_phosphate.writePDB(str(workdir / 'no_phosphate.pdb'))
    dna = DNA.fromCoarsePDB(str(workdir / 'no_phosphate.pdb'), circular=True,
                            output_pdb=str(workdir / 'from_coarse.pdb'))
    assert len(dna.atoms) == 6 * n
    assert len(dna.bonds) == 6 * n
    rebuilt = dna.atoms[(dna.atoms.chainID == 'A') & (dna.atoms.resSeq == 1) & (dna.atoms.name == 'P')]
    original = plasmid.atoms[(plasmid.atoms.chainID == 'A') & (plasmid.atoms.resSeq == 1) & (plasmid.atoms.name == 'P')]
    assert np.linalg.norm(rebuilt[['x', 'y', 'z']].values - original[['x', 'y', 'z']].values) < 1.5  # Angstrom


def test_bend_linear_structure(workdir):
    """A straight duplex is bent into a closed circle when loaded as circular."""
    DNA.fromSequence(SEQUENCE, output_pdb=str(workdir / 'straight.pdb'),
                     atomistic_template=str(workdir / 'straight_atomistic.pdb'))
    dna = DNA.fromPDB(str(workdir / 'straight_atomistic.pdb'), circular=True, delta_linking_number=-1,
                      output_pdb=str(workdir / 'bent.pdb'))
    n = len(SEQUENCE)
    assert len(dna.bonds) == 6 * n
    lk = dna.bend.linking_number
    assert lk == round(n / 10.5) - 1
    assert dna.computeLinkingNumber() == pytest.approx(lk, abs=1e-6)
    initial, final, xyz = _minimize(dna)
    assert _bond_deviation(dna, xyz).max() < 0.05
    assert dna.computeLinkingNumber(xyz) == pytest.approx(lk, abs=1e-6)


def test_joined_strands(workdir):
    """Joining the strands gives a single circular chain running along both sides of the duplex."""
    dna = DNA.fromSequence(SEQUENCE, circular=True, join=[['A', 'B']], output_pdb=str(workdir / 'joined.pdb'),
                           atomistic_template=str(workdir / 'joined_atomistic.pdb'))
    n = len(SEQUENCE)
    assert dna.circular_chains == ['A']
    assert list(dna.atoms.chainID.unique()) == ['A']
    assert sorted(dna.atoms.resSeq.unique()) == list(range(1, 2 * n + 1))
    assert len(dna.bonds) == 6 * n
    initial, final, xyz = _minimize(dna, iterations=2000)
    assert final < initial
    assert _bond_deviation(dna, xyz).max() < 0.1
