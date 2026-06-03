"""Executable specification for the opt-in bias forces.

Each test states one contract in its docstring and pins it with an analytical
oracle on a real structure. Target offsets use ``_OFFSET``, three distinct
nonzero components, so a wrong axis, component, or sign fails a single witness.
"""

import numpy as np
import pytest
import openmm.unit as unit

from open3SPN2.force.bias import (PositionRestraint, StringProteinDNA,
                                  BasePairProteinHarmonicBias, ElectrostaticBias)
from open3SPN2.force.protein_dna import ElectrostaticsProteinDNA, select_string_groups
from open3SPN2.force import force_groups as fg

_kJ_nm2 = unit.kilojoule_per_mole / unit.nanometer ** 2
_OFFSET = np.array([0.13, -0.21, 0.37])


def _target(point):
    return dict(x0=point[0] * unit.nanometer, y0=point[1] * unit.nanometer, z0=point[2] * unit.nanometer)


def test_position_restraint_is_half_k_distance_squared(dna_harness):
    """Restraint energy is 0.5*k*|centroid - r0|^2 about the centroid of all atoms."""
    centroid = dna_harness.geometric_centroid(dna_harness.atom_indices())
    k = 1000.0
    force = PositionRestraint(dna_harness.dna, k=k * _kJ_nm2, **_target(centroid + _OFFSET))
    assert dna_harness.energy(force, fg.MEMBRANE) == pytest.approx(0.5 * k * _OFFSET.dot(_OFFSET))


def test_position_restraint_acts_on_selected_residues_only(dna_harness):
    """The restrained point is the centroid of the selected residues, not of all atoms."""
    subset_centroid = dna_harness.geometric_centroid(dna_harness.atom_indices(resSeq=1))
    assert not np.allclose(subset_centroid,
                           dna_harness.geometric_centroid(dna_harness.atom_indices()))
    force = PositionRestraint(dna_harness.dna, k=1000 * _kJ_nm2,
                              appliedToResidues=[1], **_target(subset_centroid))
    assert dna_harness.energy(force, fg.MEMBRANE) == pytest.approx(0.0, abs=1e-6)


def test_position_restraint_force_is_negative_energy_gradient(dna_harness):
    """The restraint force equals -dE/dx (finite-difference check)."""
    centroid = dna_harness.geometric_centroid(dna_harness.atom_indices(resSeq=2))
    force = PositionRestraint(dna_harness.dna, k=500 * _kJ_nm2,
                              appliedToResidues=[2], **_target(centroid + _OFFSET))
    assert dna_harness.max_force_error(force, fg.MEMBRANE) < 1e-4


def test_position_restraint_rejects_empty_selection(dna_harness):
    """Selecting no atoms is rejected with a clear error rather than a silent empty force."""
    with pytest.raises(ValueError):
        PositionRestraint(dna_harness.dna, appliedToResidues=[9999])


def test_string_potential_is_half_k_distance_from_r0_squared(protein_dna_harness):
    """String energy is 0.5*k*(d - r0)^2 in the centroid-centroid distance d."""
    harness = protein_dna_harness
    protein_idx, dna_idx = select_string_groups(harness.dna, '34', '12')
    distance = np.linalg.norm(harness.mass_centroid(protein_idx) - harness.mass_centroid(dna_idx))
    k, r0 = 41.84, distance - 0.7
    force = StringProteinDNA(harness.dna, None, r0=r0, chain_protein='34', chain_DNA='12', k_string_PD=k)
    assert harness.energy(force, fg.PULLING) == pytest.approx(0.5 * k * (distance - r0) ** 2, rel=1e-4)


def test_base_pair_bias_is_half_k_projection_squared(protein_dna_harness):
    """Base-pair bias energy is 0.5*k*i^2 in the protein's base-pair projection offset i."""
    harness = protein_dna_harness
    sugars = harness.atom_indices(chainID='1', name='S', nonzero_mass=True)
    left, right = [sugars[0]], [sugars[5]]
    protein = [harness.atom_indices(chainID='3', name='CA', nonzero_mass=True)[0]]
    bpl, bpr, point = harness.posnm[left[0]], harness.posnm[right[0]], harness.posnm[protein[0]]
    sep, k = 4, 4.184
    offset = sep * ((point - bpl).dot(bpr - bpl) / (bpr - bpl).dot(bpr - bpl) - 0.5)
    force = BasePairProteinHarmonicBias(harness.dna, None, left, right, protein, base_pair_sep=sep, k=k)
    assert harness.energy(force, fg.PROTEIN_DNA_BIAS) == pytest.approx(0.5 * k * offset ** 2, rel=1e-4)


def test_electrostatic_bias_squares_the_scaled_energy_deviation(protein_dna_harness):
    """Electrostatic bias is 0.5*k_ebias*((E_elec - center)/4.184)^2 over the protein-DNA electrostatic energy."""
    harness = protein_dna_harness
    k_elec, ldby = 1.0, 1.2 * unit.nanometer
    e_elec = harness.energy(ElectrostaticsProteinDNA(harness.dna, harness.protein, k=k_elec, ldby=ldby),
                            fg.ELECTROSTATICS_PROTEIN_DNA)
    k_ebias, center = 2.0, e_elec - 5.0
    bias = ElectrostaticBias(harness.dna, harness.protein, k_ebias=k_ebias, center=center,
                             k_elec=k_elec, ldby=ldby)
    assert harness.energy(bias, fg.PROTEIN_DNA_BIAS) == pytest.approx(
        0.5 * k_ebias * ((e_elec - center) / 4.184) ** 2, rel=1e-4)
