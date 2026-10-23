"""Tests for circular DNA (plasmids).

* The bonded interactions of a circular chain wrap around its end: every strand of a closed duplex
  has as many bonds as beads, and the closing bonds are near their equilibrium length.
* The equilibrium values of the closing interactions come from the same template as a linear DNA
  with the circular sequence rotated.
* The sequence separation masks of the nonbonded and hydrogen bond forces follow the topology: they
  count residues around a chain that the topology closes, and only then.
* CrossStacking takes its neighboring bases from the topology, so a closed chain has cross stacking
  across its closure.
* The bent circle has the requested linking number, and minimizing it keeps it closed.
* A circular structure loaded as linear (atomistic or coarse grained, e.g. the PDB minicircle 13GQ)
  is closed by DNA.circularize, which rebuilds the phosphate of the first residue at the position of
  the phosphate group bridging the last and first residues.
* The core of the force field (readers, coarse graining, topology and forces) has no code specific to
  circular DNA: circular DNA is an add-on that patches the topology.
"""
import inspect
import itertools
import re
from pathlib import Path

import numpy as np
import openmm.unit as unit
import pandas
import pytest
from openmm.app import PDBFile

import open3SPN2.force.dna
from open3SPN2 import DNA, System, forces, pdb2table
from open3SPN2.circular import linking_number
from open3SPN2.force.dna import sequence_parameters

SEQUENCE = 'ATACAAAGGTGCGAGGTTTCTATGCTCCCACGGTCAATGCATTGAGCCATTAGGCTACGTTAGCATCGATT'  # 72 bp
MINICIRCLE = Path(__file__).parent / '13gq' / '13GQ.pdb'


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


def _bead(dna, chain, resseq, name):
    atoms = dna.atoms
    return atoms[(atoms.chainID == chain) & (atoms.resSeq == resseq) & (atoms.name == name)].index[0]


def _center_of_mass(atoms):
    masses = atoms['name'].str[0].map({'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'P': 30.973762}).to_numpy()
    return (atoms[['x', 'y', 'z']].to_numpy() * masses[:, None]).sum(axis=0) / masses.sum()


def test_circular_topology(plasmid):
    """Feature: a closed duplex has a bond, two angles... per bead, and closing bonds near equilibrium."""
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
    for chain in ['A', 'B']:
        closing = plasmid.bonds[(plasmid.bonds.aai == _bead(plasmid, chain, n, 'S')) &
                                (plasmid.bonds.aaj == _bead(plasmid, chain, 1, 'P'))]
        assert len(closing) == 1
    assert _bond_deviation(plasmid).max() < 0.1  # strain of a 72 bp minicircle (radius ~4 nm)


def test_linear_topology_is_unchanged(workdir):
    """Feature: a linear DNA keeps its open ends."""
    dna = DNA.fromSequence(SEQUENCE, output_pdb=str(workdir / 'linear.pdb'),
                           atomistic_template=str(workdir / 'linear_atomistic.pdb'))
    n = len(SEQUENCE)
    assert len(dna.atoms) == 6 * n - 2
    assert len(dna.bonds) == 6 * n - 4


def test_closing_parameters_match_the_linear_template(plasmid, workdir):
    """Feature: the equilibrium values of the interactions around the closure come from the same template
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
    n = len(SEQUENCE)
    assert (('A', n - 1, 'S'), ('A', 0, 'P')) in common
    for k in common:
        assert circular[k] == pytest.approx(rotated_bonds[k], abs=1e-6)


@pytest.mark.parametrize('delta', [0, -2, 1])
def test_linking_number(workdir, delta):
    """Feature: a circle built from sequence has the requested linking number, which minimization keeps."""
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
    """Feature: the linking number of two closed curves is their Gauss linking integral."""
    t = np.linspace(0, 2 * np.pi, 50, endpoint=False)
    ring = np.stack([np.cos(t), np.sin(t), 0 * t], axis=1)
    linked = np.stack([1 + np.cos(t), 0 * t, np.sin(t)], axis=1)
    apart = linked + [5, 0, 0]
    assert abs(linking_number(ring, linked)) == pytest.approx(1)
    assert linking_number(ring, apart) == pytest.approx(0, abs=1e-9)


def test_sequence_separation_follows_topology(plasmid):
    """Feature: residues near the closure of a chain that the topology closes are close in sequence,
    residues on different chains are never close."""
    parameters = sequence_parameters(plasmid)
    seqid, period = parameters['seqid'], parameters['period']
    n = len(SEQUENCE)
    assert set(period) == {n}
    atoms = plasmid.atoms

    def separation(i, j):
        d = abs(seqid[i] - seqid[j])
        length = max(period[i], period[j])
        return min(d, length - d) if d < length else d

    for i, j in itertools.combinations(range(0, len(atoms), 7), 2):
        if atoms.chainID[i] == atoms.chainID[j]:
            d = abs(atoms.resSeq[i] - atoms.resSeq[j])
            assert separation(i, j) == min(d, n - d)
        else:
            assert separation(i, j) >= 3


def _mask_terms(force_name, dna):
    """The sequence mask expression of a force and the per-particle names it reads."""
    wrapper = forces[force_name](dna)
    openmm_forces = [wrapper.force] if force_name in ['Exclusion', 'Electrostatics'] else (
        list(wrapper.forces.values()) if force_name == 'BasePair' else [f for pair in wrapper.crossStackingForces.values() for f in pair])
    return re.search(r'seqsep\s*=\s*([^;]*);', openmm_forces[0].getEnergyFunction()).group(1)


@pytest.mark.parametrize('force_name', ['Exclusion', 'Electrostatics', 'BasePair', 'CrossStacking'])
def test_masks_follow_topology(plasmid, force_name):
    """Feature: the sequence mask of the forces counts residues around a chain closed by the topology, and
    stops doing so if the topology is linear."""
    n = len(SEQUENCE)
    first, last, middle = _bead(plasmid, 'A', 1, 'S'), _bead(plasmid, 'A', n, 'S'), _bead(plasmid, 'A', n // 2, 'S')
    mask = _mask_terms(force_name, plasmid)
    names = ('seqid_d', 'seqid_a', 'period_d', 'period_a') if 'seqid_d' in mask else (
        'seqid1', 'seqid2', 'period1', 'period2')

    def evaluate(parameters, i, j):
        values = dict(zip(names, [parameters['seqid'][i], parameters['seqid'][j],
                                  parameters['period'][i], parameters['period'][j]]))
        return eval(mask, {'abs': abs, 'max': max, 'step': lambda x: float(x >= 0)}, values)

    closed = sequence_parameters(plasmid)
    assert evaluate(closed, first, last) == 0
    assert evaluate(closed, first, middle) == 1

    linear = DNA.fromAtoms(plasmid.atoms, plasmid.DNAtype, output_pdb=None)
    assert evaluate(sequence_parameters(linear), first, last) == 1


def _crossstacking_donors(dna):
    wrapper = forces['CrossStacking'](dna)
    triplets = []
    for donor_force, _ in wrapper.crossStackingForces.values():
        for i in range(donor_force.getNumDonors()):
            triplets.append(tuple(donor_force.getDonorParameters(i)[:3]))
    return triplets


def test_crossstacking_sites_follow_topology(plasmid, workdir):
    """Feature: CrossStacking has one donor per stacked base step of the topology, so a closed chain also
    cross stacks across its closure."""
    linear = DNA.fromSequence(SEQUENCE, output_pdb=str(workdir / 'linear_cs.pdb'), atomistic_template=None)
    assert len(_crossstacking_donors(linear)) == len(linear.stackings)
    donors = _crossstacking_donors(plasmid)
    assert len(donors) == len(plasmid.stackings)
    n = len(SEQUENCE)
    base = lambda resseq: plasmid.atoms[(plasmid.atoms.chainID == 'A') & (plasmid.atoms.resSeq == resseq)
                                        & plasmid.atoms.name.isin(list('ATGC'))].index[0]
    assert (base(n), _bead(plasmid, 'A', n, 'S'), base(1)) in donors


@pytest.mark.parametrize('with_atomistic, tolerance', [(True, 1e-3), (False, 1.5)])
def test_circularize_patches_a_loaded_minicircle(workdir, with_atomistic, tolerance):
    """Feature: a circular structure loaded as linear (the PDB minicircle 13GQ) is closed by circularize. Given
    the atomistic structure, the phosphate of the first residue is the phosphate group bridging the last and
    first residues (as every phosphate bead); without it, it is rebuilt from a phosphate step of the chain."""
    dna = DNA.fromPDB(str(MINICIRCLE), output_pdb=str(workdir / '13gq_linear.pdb'))
    linear_bonds = len(dna.bonds)
    dna.circularize(atomistic=str(MINICIRCLE) if with_atomistic else None,
                    output_pdb=str(workdir / '13gq_circular.pdb'))
    atomistic = pdb2table(PDBFile(str(MINICIRCLE)))

    assert dna.circular_chains == ['i', 'j']
    assert len(dna.atoms) == 6 * 95
    assert len(dna.bonds) == 6 * 95 == linear_bonds + 4
    deviation = _bond_deviation(dna)
    for chain in ['i', 'j']:
        closing = (dna.bonds.aai == _bead(dna, chain, 95, 'S')) & (dna.bonds.aaj == _bead(dna, chain, 1, 'P'))
        assert closing.sum() == 1
        assert deviation[closing.to_numpy()].max() <= deviation.max()
        residues = atomistic[atomistic.chainID == chain]
        group = pandas.concat([residues[(residues.resSeq == 95) & (residues.name == "O3'")],
                               residues[(residues.resSeq == 1) & residues.name.isin(['P', 'OP1', 'OP2', "O5'"])]])
        rebuilt = dna.atoms.loc[_bead(dna, chain, 1, 'P'), ['x', 'y', 'z']].to_numpy(dtype=float)
        assert np.linalg.norm(rebuilt - _center_of_mass(group)) < tolerance  # Angstrom
    lk = dna.computeLinkingNumber()
    assert lk == pytest.approx(round(lk), abs=1e-6)


def test_circularize_a_coarse_grained_circle_without_first_phosphate(plasmid, workdir):
    """Feature: a coarse grained circle that lost the phosphate of its first residue is closed again, with
    the phosphate rebuilt near its original position."""
    n = len(SEQUENCE)
    coarse = plasmid.atoms
    coarse = coarse[~((coarse.chainID == 'A') & (coarse.resSeq == 1) & (coarse.name == 'P'))]
    no_phosphate = DNA()
    no_phosphate.atoms = coarse.reset_index(drop=True)
    no_phosphate.writePDB(str(workdir / 'no_phosphate.pdb'))

    dna = DNA.fromCoarsePDB(str(workdir / 'no_phosphate.pdb'))
    dna.circularize(output_pdb=str(workdir / 'from_coarse.pdb'))

    assert len(dna.atoms) == 6 * n
    assert len(dna.bonds) == 6 * n
    rebuilt = dna.atoms.loc[_bead(dna, 'A', 1, 'P'), ['x', 'y', 'z']].to_numpy(dtype=float)
    original = plasmid.atoms.loc[_bead(plasmid, 'A', 1, 'P'), ['x', 'y', 'z']].to_numpy(dtype=float)
    assert np.linalg.norm(rebuilt - original) < 1.5  # Angstrom


def test_bend_and_close_a_linear_structure(workdir):
    """Feature: a straight duplex can be bent into a circle (DNA.bendToCircle) and then closed."""
    DNA.fromSequence(SEQUENCE, output_pdb=str(workdir / 'straight.pdb'),
                     atomistic_template=str(workdir / 'straight_atomistic.pdb'))
    dna = DNA.fromPDB(str(workdir / 'straight_atomistic.pdb'), output_pdb=None)
    dna.bendToCircle(delta_linking_number=-1)
    dna.circularize(output_pdb=str(workdir / 'bent.pdb'))
    n = len(SEQUENCE)
    assert len(dna.bonds) == 6 * n
    lk = dna.bend.linking_number
    assert lk == round(n / 10.5) - 1
    assert dna.computeLinkingNumber() == pytest.approx(lk, abs=1e-6)
    initial, final, xyz = _minimize(dna)
    assert _bond_deviation(dna, xyz).max() < 0.05
    assert dna.computeLinkingNumber(xyz) == pytest.approx(lk, abs=1e-6)


def test_joined_strands(workdir):
    """Feature: joining the strands gives a single circular chain running along both sides of the duplex."""
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


def test_core_is_unaware_of_circular_dna():
    """Feature: circular DNA is an add-on. The readers, coarse graining, topology and forces of the core have no
    code specific to circular DNA; they read the connectivity from the topology."""
    core = [open3SPN2.force.dna, DNA.computeTopology, DNA.computeGeometry, DNA.mapTemplate, DNA.CoarseGrain,
            DNA.checkAtoms, DNA.fromAtoms, DNA.fromPDB, DNA.fromCoarsePDB, DNA.fromXYZ, DNA.create_mutant]
    for item in core:
        assert 'circular' not in inspect.getsource(item).lower(), item
