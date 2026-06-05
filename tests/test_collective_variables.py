"""Executable specification for the opt-in collective-variable readouts.

Each readout reports a geometric value as its "energy". Each test states one
contract in its docstring and pins it with an analytical oracle on a real
structure. ``_OFFSET`` has three distinct nonzero components so a wrong axis or
sign fails a single witness.
"""

import numpy as np
import pytest
import openmm.unit as unit

from open3SPN2.force.collective_variables import (DistanceFromPoint, StringLength,
                                                  BasePairProteinPosition)
from open3SPN2.force.protein_dna import select_string_groups
from open3SPN2.force import force_groups as fg

_OFFSET = np.array([0.13, -0.21, 0.37])


def test_distance_from_point_reports_centroid_distance(dna_harness):
    """Reports the Euclidean distance from the selected-atom centroid to r0."""
    centroid = dna_harness.geometric_centroid(dna_harness.atom_indices())
    cv = DistanceFromPoint(dna_harness.dna,
                           x0=(centroid[0] + _OFFSET[0]) * unit.nanometer,
                           y0=(centroid[1] + _OFFSET[1]) * unit.nanometer,
                           z0=(centroid[2] + _OFFSET[2]) * unit.nanometer)
    assert dna_harness.energy(cv, fg.MEASUREMENT) == pytest.approx(np.linalg.norm(_OFFSET))


def test_string_length_reports_centroid_distance(protein_dna_harness):
    """Reports the mass-weighted centroid-centroid distance of the protein and DNA groups."""
    harness = protein_dna_harness
    protein_idx, dna_idx = select_string_groups(harness.dna, '34', '12')
    expected = np.linalg.norm(harness.mass_centroid(protein_idx) - harness.mass_centroid(dna_idx))
    cv = StringLength(harness.dna, None, chain_protein='34', chain_DNA='12')
    assert harness.energy(cv, fg.MEASUREMENT) == pytest.approx(expected, rel=1e-4)


def test_string_length_pools_the_requested_chains(protein_dna_harness):
    """Chain spec accepts one or many chains: 'AB' equals ['A','B'] and differs from a single chain."""
    harness = protein_dna_harness
    pooled = harness.energy(StringLength(harness.dna, None, chain_protein='34', chain_DNA='12'), fg.MEASUREMENT)
    listed = harness.energy(StringLength(harness.dna, None, chain_protein=['3', '4'], chain_DNA=['1', '2']), fg.MEASUREMENT)
    single = harness.energy(StringLength(harness.dna, None, chain_protein='3', chain_DNA='1'), fg.MEASUREMENT)
    assert pooled == pytest.approx(listed, rel=1e-9)
    assert pooled != pytest.approx(single, rel=1e-9)


def test_base_pair_position_reports_projection_in_base_pairs(protein_dna_harness):
    """Reports the protein point's projection onto the base-pair axis as a base-pair offset."""
    harness = protein_dna_harness
    sugars = harness.atom_indices(chainID='1', name='S', nonzero_mass=True)
    left, right = [sugars[0]], [sugars[5]]
    protein = [harness.atom_indices(chainID='3', name='CA', nonzero_mass=True)[0]]
    bpl, bpr, point = harness.posnm[left[0]], harness.posnm[right[0]], harness.posnm[protein[0]]
    sep = 4
    expected = sep * ((point - bpl).dot(bpr - bpl) / (bpr - bpl).dot(bpr - bpl) - 0.5)
    cv = BasePairProteinPosition(harness.dna, None, left, right, protein, base_pair_sep=sep)
    assert harness.energy(cv, fg.MEASUREMENT) == pytest.approx(expected, rel=1e-4)


def test_readout_stays_out_of_the_energy_total(dna_harness):
    """A readout's force group is excluded from the groups summed into the total energy."""
    cv = DistanceFromPoint(dna_harness.dna)
    assert cv.force.getForceGroup() not in fg.TOTAL_ENERGY
