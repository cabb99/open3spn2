"""Tests for infinite DNA: a straight duplex whose ends are bonded to each other through the periodic box.

* A System built with a periodic box uses that box (default box vectors and topology unit cell) and turns on
  periodic boundary conditions for the DNA forces.
* DNA.fromSequence(infinite=True) builds a straight duplex along z whose length is the box length along z, with
  both strands closed through the box: as many bonds as beads, and every bond near its equilibrium length when
  measured to the nearest periodic image.
* Over one period each strand winds an integer number of turns around the axis, the relaxed number of helical
  turns plus `delta_linking_number`, so the strands join their periodic images.
* The infinite DNA can be minimized and keeps its bonds.
"""
import numpy as np
import openmm
import openmm.unit as unit
import pytest

from open3SPN2 import DNA, System

SEQUENCE = 'ATACAAAGGTGCGAGGTTTCTATGCTCCCACGGTCAATGCATTG'  # 44 bp, about 4 turns


@pytest.fixture(scope='module')
def workdir(tmp_path_factory):
    return tmp_path_factory.mktemp('infinite')


@pytest.fixture(scope='module')
def infinite(workdir):
    return DNA.fromSequence(SEQUENCE, infinite=True, output_pdb=str(workdir / 'infinite.pdb'),
                            atomistic_template=str(workdir / 'infinite_atomistic.pdb'))


def _minimum_image_bond_deviation(dna, box, xyz=None):
    xyz = dna.atoms[['x', 'y', 'z']].to_numpy(dtype=float) / 10 if xyz is None else xyz
    d = xyz[dna.bonds['aai']] - xyz[dna.bonds['aaj']]
    d -= np.round(d / box) * box
    return np.abs(np.linalg.norm(d, axis=1) - dna.bonds['r0'].to_numpy())


def test_system_uses_the_periodic_box(workdir):
    """Feature: a System built with a periodic box (nm) uses it as default box and topology unit cell, and the
    DNA forces use periodic boundary conditions."""
    dna = DNA.fromSequence(SEQUENCE, output_pdb=str(workdir / 'linear.pdb'), atomistic_template=None)
    system = System(dna, periodicBox=[15, 16, 17])
    vectors = [v.value_in_unit(unit.nanometer) for v in system.getDefaultPeriodicBoxVectors()]
    np.testing.assert_allclose(np.array([[v[0], v[1], v[2]] for v in vectors]), np.diag([15, 16, 17]))
    np.testing.assert_allclose(system.top.getUnitCellDimensions().value_in_unit(unit.nanometer), [15, 16, 17])
    assert dna.periodic
    system.add3SPN2forces()
    assert system.forces['Bond'].force.usesPeriodicBoundaryConditions()


def test_infinite_dna_closes_through_the_box(infinite):
    """Feature: the infinite duplex has both strands closed through the box, with every bond near equilibrium
    when measured to the nearest periodic image."""
    n = len(SEQUENCE)
    box = np.array(infinite.periodic_box)
    assert infinite.circular_chains == ['A', 'B']
    assert len(infinite.atoms) == 6 * n
    assert len(infinite.bonds) == 6 * n
    assert len(infinite.stackings) == 2 * n
    assert box[2] == pytest.approx(n * 0.338, rel=0.05)
    assert min(box[:2]) > 10  # the electrostatics cutoff is 5 nm
    xyz = infinite.atoms[['x', 'y', 'z']].to_numpy(dtype=float) / 10
    np.testing.assert_allclose(xyz[:, :2].mean(axis=0), box[:2] / 2, atol=0.1)
    # Straightening the intrinsically curved template strains the steps a little; the closing bonds through the
    # box are no more strained than the other steps
    deviation = _minimum_image_bond_deviation(infinite, box)
    assert deviation.max() < 0.06
    atoms = infinite.atoms
    for chain, last_sugar, first_phosphate in [('A', n, 1), ('B', n, 1)]:
        closing = (infinite.bonds.aai == atoms[(atoms.chainID == chain) & (atoms.resSeq == last_sugar) & (atoms.name == 'S')].index[0]) & \
                  (infinite.bonds.aaj == atoms[(atoms.chainID == chain) & (atoms.resSeq == first_phosphate) & (atoms.name == 'P')].index[0])
        assert closing.sum() == 1
        assert deviation[closing.to_numpy()].max() <= deviation.max()


@pytest.mark.parametrize('delta', [0, 1, -1])
def test_strands_wind_an_integer_number_of_turns(workdir, delta):
    """Feature: over one period each strand winds around the axis the relaxed number of helical turns plus
    delta_linking_number."""
    n = len(SEQUENCE)
    dna = DNA.fromSequence(SEQUENCE, infinite=True, delta_linking_number=delta,
                           output_pdb=str(workdir / f'turns{delta}.pdb'), atomistic_template=None)
    box = np.array(dna.periodic_box) * 10
    atoms = dna.atoms
    sugars = atoms[(atoms.chainID == 'A') & (atoms.name == 'S')].sort_values('resSeq')
    angles = np.arctan2(sugars.y - box[1] / 2, sugars.x - box[0] / 2).to_numpy()
    steps = np.diff(np.concatenate([angles, angles[:1]]))
    turns = np.sum((steps + np.pi) % (2 * np.pi) - np.pi) / (2 * np.pi)
    assert dna.helix.relaxed_turns == round(n / 10.5)
    assert turns == pytest.approx(dna.helix.relaxed_turns + delta, abs=1e-6)


def test_infinite_dna_minimizes_and_keeps_its_bonds(infinite):
    """Feature: the infinite DNA can be minimized in its periodic box and its bonds stay near equilibrium."""
    system = System(infinite, periodicBox=infinite.periodic_box)
    system.add3SPN2forces()
    simulation = system.initializeMD(platform_name='Reference')
    initial = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    simulation.minimizeEnergy(maxIterations=500)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    final = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    assert np.isfinite(initial) and final < initial
    xyz = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    assert _minimum_image_bond_deviation(infinite, np.array(infinite.periodic_box), xyz).max() < 0.05
