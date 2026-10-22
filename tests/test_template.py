"""Tests for the equilibrium geometry of curved B-DNA.

3SPN2.C takes the equilibrium bonds, angles and dihedrals of B_curved DNA from an ideal atomistic
template that is built from the sequence and coarse-grained in memory. The behavioural requirements
tested here:

* The built-in atomistic template uses the standard (PDB v3) atom names that pdbfixer and OpenMM use.
* checkAtoms reports heavy atoms that have no coarse-grained bead (unmapped) and expected heavy atoms
  that are absent (missing), ignoring hydrogens and the atoms that 3SPN2 drops at the chain ends.
* Coarse-graining warns, naming them, about unmapped and missing heavy atoms.
* fromPDB always sends external structures through pdbfixer, which standardizes legacy atom names and
  rebuilds missing heavy atoms and missing internal residues.
* computeGeometry repairs the atomistic template with pdbfixer only when checkAtoms finds problems, so
  a custom template with legacy names gives the same equilibrium geometry and the ideal one is never
  sent through pdbfixer.
* fromAtoms builds a DNA from a coarse-grained table; fromSequence gives the same beads as the
  atomistic route through pdbfixer.
* A B_curved DNA built from sequence sits at the minimum of its bonds, angles and dihedrals.
* The equilibrium values depend only on the sequence, not on the coordinates of the input structure.
"""
import warnings
from pathlib import Path

import numpy as np
import pandas
import pdbfixer
import pytest
from openmm.app import PDBFile

import open3SPN2
from open3SPN2 import DNA, fixPDB, pdb2table, writePDB

SEQUENCE = 'AACAGATCCGCTGGTTA'  # contains all 16 dinucleotide steps


def _ideal_template():
    """Returns the coarse-grained beads and the atomistic table of the ideal template for SEQUENCE."""
    dna = DNA()
    dna.parseConfigurationFile()
    beads = dna.computeGeometry(pandas.Series(list(SEQUENCE)))
    return beads, dna.template


def _from_sequence(folder):
    return DNA.fromSequence(SEQUENCE, 'B_curved', output_pdb=str(folder / 'clean.pdb'),
                            atomistic_template=str(folder / 'atomistic.pdb'))


def _wrapped(degrees):
    return (degrees + 180) % 360 - 180


def _terms_off(table, deviation, atol):
    """Names the interaction types (e.g. 'P-S-T') whose deviation exceeds atol."""
    atoms = [c for c in ['i', 'j', 'k', 'l'] if c in table]
    off = np.abs(np.asarray(deviation)) > atol
    return sorted({'-'.join(map(str, row)) for row in table.loc[off, atoms].itertuples(index=False)})


@pytest.mark.parametrize('base', ['A', 'C', 'G', 'T'])
def test_atomic_template_uses_standard_atom_names(base):
    """Feature: the built-in atomistic template names its heavy atoms like pdbfixer's standard residues."""
    atomic = pandas.read_csv(Path(open3SPN2.__file__).parent / 'DNA_atomic.csv')
    names = atomic.loc[atomic['base'] == base, 'name']
    standard = PDBFile(str(Path(pdbfixer.__file__).parent / 'templates' / f'D{base}.pdb'))
    standard_names = [atom.name for atom in standard.topology.atoms() if atom.element.symbol != 'H']
    assert sorted(names) == sorted(standard_names)


def test_check_atoms_accepts_ideal_template():
    """Feature: the ideal atomistic template has every heavy atom of every bead and nothing else."""
    _, atomistic = _ideal_template()
    assert DNA.checkAtoms(atomistic).empty


def test_check_atoms_reports_unmapped_and_missing():
    """Feature: checkAtoms reports heavy atoms without a bead as unmapped and absent expected heavy
    atoms as missing, and ignores hydrogens."""
    _, atomistic = _ideal_template()
    table = atomistic.reset_index(drop=True)
    table['name'] = table['name'].replace({'C7': 'C5M', 'OP1': 'O1P'})
    table = table.drop(np.flatnonzero(table['name'] == 'N3')[0])
    hydrogen = table.iloc[[0]].assign(name='HX1', element='H')
    table = pandas.concat([table, hydrogen], ignore_index=True)

    problems = DNA.checkAtoms(table)

    assert set(zip(problems['name'], problems['problem'])) == {('C5M', 'unmapped'), ('C7', 'missing'),
                                                               ('O1P', 'unmapped'), ('OP1', 'missing'),
                                                               ('N3', 'missing')}


def test_check_atoms_allows_chain_termini():
    """Feature: the 5'-terminal phosphate bead and the 3'-terminal O3', which 3SPN2 drops, may be absent."""
    _, atomistic = _ideal_template()
    table = atomistic.reset_index(drop=True)
    first = table.groupby('chainID')['resSeq'].transform('min') == table['resSeq']
    last = table.groupby('chainID')['resSeq'].transform('max') == table['resSeq']
    trimmed = table[~(first & table['name'].isin(['P', 'OP1', 'OP2', "O5'"])) & ~(last & (table['name'] == "O3'"))]
    assert DNA.checkAtoms(trimmed).empty


def test_coarse_grain_warns_about_missing_atoms():
    """Feature: coarse-graining warns, naming them, about expected heavy atoms that are absent."""
    _, atomistic = _ideal_template()
    table = atomistic.reset_index(drop=True)
    table = table.drop(np.flatnonzero(table['name'] == 'C7')[0])
    with pytest.warns(UserWarning, match='C7'):
        DNA.CoarseGrain(table)


def test_coarse_grain_warns_about_atoms_without_bead():
    """Feature: coarse-graining warns, naming them, about heavy atoms that have no coarse-grained bead."""
    _, atomistic = _ideal_template()
    unknown = atomistic.copy()
    unknown['name'] = unknown['name'].replace({'N3': 'C9X'})
    with pytest.warns(UserWarning, match='C9X'):
        DNA.CoarseGrain(unknown)


def test_from_pdb_standardizes_legacy_atom_names(tmp_path):
    """Feature: fromPDB standardizes atom names with pdbfixer, so a structure with legacy names
    (C5M, O1P, O2P) gives the same beads as the ideal template, without coarse-graining warnings."""
    beads, atomistic = _ideal_template()
    legacy = atomistic.reset_index(drop=True)
    legacy['name'] = legacy['name'].replace({'C7': 'C5M', 'OP1': 'O1P', 'OP2': 'O2P'})
    writePDB(legacy, str(tmp_path / 'legacy.pdb'))

    with warnings.catch_warnings():
        warnings.filterwarnings('error', message='Atoms do not match the coarse-grained beads')
        dna = DNA.fromPDB(str(tmp_path / 'legacy.pdb'), output_pdb=str(tmp_path / 'clean.pdb'),
                          compute_topology=False)

    assert list(dna.atoms['name']) == list(beads['name'])
    np.testing.assert_allclose(dna.atoms[['x', 'y', 'z']], beads[['x', 'y', 'z']], atol=1e-2)


def test_from_pdb_rebuilds_missing_atoms(tmp_path):
    """Feature: fromPDB rebuilds missing heavy atoms with pdbfixer, so a structure missing a base atom
    gives the beads of the complete structure, without coarse-graining warnings."""
    beads, atomistic = _ideal_template()
    incomplete = atomistic.reset_index(drop=True)
    incomplete = incomplete.drop(np.flatnonzero(incomplete['name'] == 'N3')[:3])
    writePDB(incomplete, str(tmp_path / 'incomplete.pdb'))

    with warnings.catch_warnings():
        warnings.filterwarnings('error', message='Atoms do not match the coarse-grained beads')
        dna = DNA.fromPDB(str(tmp_path / 'incomplete.pdb'), output_pdb=str(tmp_path / 'clean.pdb'),
                          compute_topology=False)

    assert list(dna.atoms['name']) == list(beads['name'])
    np.testing.assert_allclose(dna.atoms[['x', 'y', 'z']], beads[['x', 'y', 'z']], atol=0.1)


def test_from_pdb_rebuilds_missing_residues(tmp_path):
    """Feature: fromPDB always sends external structures through pdbfixer, so an internal residue listed in
    SEQRES but absent from the coordinates is rebuilt."""
    beads, atomistic = _ideal_template()
    table = atomistic.reset_index(drop=True)
    writePDB(table[~((table['chainID'] == 'A') & (table['resSeq'] == 6))].copy(), str(tmp_path / 'gapped.pdb'))
    seqres = []
    for chain, residues in table.groupby('chainID'):
        names = residues.drop_duplicates('resSeq').sort_values('resSeq')['resname'].tolist()
        for k in range(0, len(names), 13):
            seqres.append(f'SEQRES {k // 13 + 1:>3} {chain} {len(names):>4}  ' + ' '.join(f'{n:>3}' for n in names[k:k + 13]))
    (tmp_path / 'gapped_seqres.pdb').write_text('\n'.join(seqres) + '\n' + (tmp_path / 'gapped.pdb').read_text())

    dna = DNA.fromPDB(str(tmp_path / 'gapped_seqres.pdb'), output_pdb=str(tmp_path / 'clean.pdb'),
                      compute_topology=False)

    assert list(dna.atoms['name']) == list(beads['name'])


def test_fix_pdb_can_skip_hydrogens(tmp_path):
    """Feature: fixPDB can standardize a structure without adding hydrogens."""
    _, atomistic = _ideal_template()
    pdb_file = writePDB(atomistic.reset_index(drop=True), str(tmp_path / 'template.pdb'))
    fixer = fixPDB(pdb_file, add_hydrogens=False)
    assert 'H' not in {atom.element.symbol for atom in fixer.topology.atoms()}


def test_compute_geometry_standardizes_custom_template(tmp_path, monkeypatch):
    """Feature: computeGeometry standardizes the atomistic template with pdbfixer before coarse-graining,
    so a custom template with legacy atom names (C5M, O1P, O2P) gives the same beads as the built-in one."""
    beads, _ = _ideal_template()
    atomic = pandas.read_csv(Path(open3SPN2.__file__).parent / 'DNA_atomic.csv')
    atomic['name'] = atomic['name'].replace({'C7': 'C5M', 'OP1': 'O1P', 'OP2': 'O2P'})
    atomic.to_csv(tmp_path / 'DNA_atomic.csv', index=False)
    monkeypatch.setattr(open3SPN2.ff3SPN2, '__location__', tmp_path)

    with warnings.catch_warnings():
        warnings.filterwarnings('error', message='Atoms do not match the coarse-grained beads')
        custom_beads, _ = _ideal_template()

    assert list(custom_beads['name']) == list(beads['name'])
    np.testing.assert_allclose(custom_beads[['x', 'y', 'z']], beads[['x', 'y', 'z']], atol=1e-2)


def test_compute_geometry_rebuilds_missing_atoms_of_custom_template(tmp_path, monkeypatch):
    """Feature: computeGeometry rebuilds heavy atoms missing from a custom template with pdbfixer, so the
    beads match the built-in template, without coarse-graining warnings."""
    beads, _ = _ideal_template()
    atomic = pandas.read_csv(Path(open3SPN2.__file__).parent / 'DNA_atomic.csv')
    atomic = atomic[~((atomic['base'] == 'T') & (atomic['name'] == 'N3'))]
    atomic.to_csv(tmp_path / 'DNA_atomic.csv', index=False)
    monkeypatch.setattr(open3SPN2.ff3SPN2, '__location__', tmp_path)

    with warnings.catch_warnings():
        warnings.filterwarnings('error', message='Atoms do not match the coarse-grained beads')
        custom_beads, _ = _ideal_template()

    assert list(custom_beads['name']) == list(beads['name'])
    np.testing.assert_allclose(custom_beads[['x', 'y', 'z']], beads[['x', 'y', 'z']], atol=0.1)


def test_ideal_template_is_not_repaired(tmp_path, monkeypatch):
    """Feature: the ideal template and fromSequence are coarse-grained without pdbfixer, which only runs
    when checkAtoms finds problems."""
    def no_pdbfixer(*args, **kwargs):
        raise AssertionError('pdbfixer should not run on the ideal template')

    monkeypatch.setattr(open3SPN2.ff3SPN2, 'fixPDB', no_pdbfixer)
    monkeypatch.setattr(open3SPN2.ff3SPN2, 'fixTable', no_pdbfixer, raising=False)

    _ideal_template()
    _from_sequence(tmp_path)


def test_from_atoms_matches_from_coarse_pdb(tmp_path):
    """Feature: fromAtoms builds from a coarse-grained table the same DNA that fromCoarsePDB builds from
    the same beads in a file."""
    beads, _ = _ideal_template()
    writePDB(beads.copy(), str(tmp_path / 'beads.pdb'))

    from_file = DNA.fromCoarsePDB(str(tmp_path / 'beads.pdb'))
    from_table = DNA.fromAtoms(beads.copy(), output_pdb=str(tmp_path / 'clean.pdb'))

    assert list(from_table.atoms['name']) == list(from_file.atoms['name'])
    np.testing.assert_allclose(from_table.atoms[['x', 'y', 'z']], from_file.atoms[['x', 'y', 'z']], atol=1e-3)
    np.testing.assert_allclose(from_table.bonds['r0'], from_file.bonds['r0'], atol=1e-6)
    np.testing.assert_allclose(from_table.angles['t0'], from_file.angles['t0'], atol=1e-6)
    np.testing.assert_allclose(_wrapped(from_table.dihedrals['t0'] - from_file.dihedrals['t0']), 0, atol=1e-6)


def test_from_sequence_matches_atomistic_route(tmp_path):
    """Feature: fromSequence gives the same beads as coarse-graining its atomistic file after pdbfixer."""
    dna = _from_sequence(tmp_path)
    reference = DNA.CoarseGrain(pdb2table(fixPDB(str(tmp_path / 'atomistic.pdb'))))
    assert list(dna.atoms['name']) == list(reference['name'])
    np.testing.assert_allclose(dna.atoms[['x', 'y', 'z']], reference[['x', 'y', 'z']], atol=1e-3)


def test_b_curved_structure_is_at_bonded_minimum(tmp_path):
    """Feature: a B_curved DNA built from sequence sits at the minimum of its bonds, angles and
    dihedrals, because the structure and the equilibrium template share the same ideal geometry."""
    dna = _from_sequence(tmp_path)
    r0, angle_t0, dihedral_t0 = dna.bonds['r0'].copy(), dna.angles['t0'].copy(), dna.dihedrals['t0'].copy()

    dna.computeTopology(template_from_X3DNA=False)  # equilibrium values measured on the structure

    off = {'bonds': _terms_off(dna.bonds, dna.bonds['r0'] - r0, atol=1e-3),
           'angles': _terms_off(dna.angles, dna.angles['t0'] - angle_t0, atol=0.1),
           'dihedrals': _terms_off(dna.dihedrals, _wrapped(dna.dihedrals['t0'] - dihedral_t0), atol=0.1)}
    assert off == {'bonds': [], 'angles': [], 'dihedrals': []}


def test_b_curved_equilibrium_ignores_input_coordinates(tmp_path):
    """Feature: B_curved equilibrium values depend only on the sequence, so a distorted structure
    loaded with fromCoarsePDB gets the same r0 and t0 as the ideal one."""
    ideal = _from_sequence(tmp_path)
    distorted = ideal.atoms.copy()
    distorted[['x', 'y', 'z']] += np.random.default_rng(0).normal(scale=1.0, size=(len(distorted), 3))
    writePDB(distorted, str(tmp_path / 'distorted.pdb'))

    loaded = DNA.fromCoarsePDB(str(tmp_path / 'distorted.pdb'))

    np.testing.assert_allclose(loaded.bonds['r0'], ideal.bonds['r0'], atol=1e-6)
    np.testing.assert_allclose(loaded.angles['t0'], ideal.angles['t0'], atol=1e-6)
    np.testing.assert_allclose(_wrapped(loaded.dihedrals['t0'] - ideal.dihedrals['t0']), 0, atol=1e-6)
