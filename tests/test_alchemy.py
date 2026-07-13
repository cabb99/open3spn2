"""Tests for alchemical base transformations.

The contract for an :class:`AlchemicalTransformation` is a single-topology hybrid Hamiltonian
``U(lambda) = (1 - lambda) * U_initial + lambda * U_target``. The behavioural requirements tested here:

* Building a plain DNA (no transformation) is unchanged (the LAMMPS energy tests still hold).
* ``DNA.create_mutant`` returns a mutant copy without touching the original.
* At ``lambda = 0`` every force-group energy equals a natively built initial structure; at
  ``lambda = 1`` it equals a natively built target structure (the core correctness contract).
* The energy interpolates linearly in ``lambda`` and ``energy_derivative`` equals ``U(1) - U(0)``.
* The hybrid builds a Context on the CPU platform (shared exclusion list) and runs without NaN forces.
"""
import tempfile
import numpy as np
import openmm.unit as unit
import pytest
from pathlib import Path

from open3SPN2 import DNA, System, AlchemicalTransformation

test_path = Path(__file__).parent
# Scratch directory for the intermediate PDBs the builders write, kept out of the repo tree.
_outdir = Path(tempfile.mkdtemp(prefix='open3spn2_alchemy_'))

FORCE_GROUPS = {'Bond': 6, 'Angle': 7, 'Stacking': 8, 'Dihedral': 9,
                'BasePair': 10, 'CrossStacking': 11, 'Exclusion': 12, 'Electrostatics': 13}

# (reference structure, DNA type)
STRUCTURES = [
    (test_path / 'adna' / 'in00_conf.xyz', 'A'),
    (test_path / 'bdna' / 'in00_conf.xyz', 'B'),
    (test_path / 'bdna_curv' / 'in00_conf.xyz', 'B_curved'),
]

# transition (purine<->purine / pyrimidine<->pyrimidine) and transversions
_TARGETS = {'A': ['G', 'C', 'T'], 'T': ['C', 'A', 'G'],
            'G': ['A', 'T', 'C'], 'C': ['T', 'G', 'A']}
_COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def _load(xyz, dna_type, tag):
    return DNA.fromXYZ(str(xyz), dna_type, template_from_X3DNA=False,
                       output_pdb=str(_outdir / f'{tag}.pdb'))


def _coords(dna):
    return (np.array(dna.atoms[['x', 'y', 'z']]) / 10.0) * unit.nanometer


def _native(xyz, dna_type, mutations, tag):
    """A natively built structure with the given base identities (via DNA.create_mutant)."""
    dna = _load(xyz, dna_type, tag)
    if mutations:
        dna = dna.create_mutant(mutations)
        dna.writePDB(str(_outdir / f'{tag}.pdb'))
    return dna


def _group_energies(dna, coords):
    """Per-force-group energy of a natively built (non-alchemical) DNA at the given coordinates."""
    system = System(dna, periodicBox=None)
    system.add3SPN2forces()
    system.initializeMD(platform_name='Reference')
    system.simulation.context.setPositions(coords)
    return _by_group(system.simulation.context)


def _alchemical(xyz, dna_type, mutations, tag):
    """Return (transformation, system) for an alchemical hybrid on the Reference platform."""
    dna = _load(xyz, dna_type, tag)
    transformation = AlchemicalTransformation(dna, mutations)
    system = System(dna, periodicBox=None)
    transformation.add_forces(system)
    system.initializeMD(platform_name='Reference')
    return transformation, system


def _by_group(context):
    return {name: context.getState(getEnergy=True, groups=2 ** g)
            .getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            for name, g in FORCE_GROUPS.items()}


def _group_energies_at(transformation, system, coords, lam):
    system.simulation.context.setPositions(coords)
    transformation.set_lambda(system.simulation.context, lam)
    return _by_group(system.simulation.context)


def _total(system, energy_unit=unit.kilojoule_per_mole):
    return system.simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(energy_unit)


def _pick_base(dna, position=4):
    """Returns (chainID, resSeq, base) of a mutable base bead on the first chain."""
    bases = dna.atoms[dna.atoms['name'].isin(['A', 'T', 'G', 'C'])]
    first_chain = bases[bases['chainID'] == bases['chainID'].min()]
    row = first_chain.iloc[position]
    return row['chainID'], int(row['resSeq']), row['name']


_CASES = [(xyz, dna_type) for xyz, dna_type in STRUCTURES]
_IDS = [s[1] for s in STRUCTURES]


@pytest.mark.parametrize("xyz, dna_type", _CASES, ids=_IDS)
def test_plain_dna_is_unchanged(xyz, dna_type):
    """A plain DNA (no transformation) builds the eight standard forces and has no lambda parameter."""
    dna = _load(xyz, dna_type, 'plain')
    system = System(dna, periodicBox=None)
    system.add3SPN2forces()
    assert set(system.forces) == set(FORCE_GROUPS)
    system.initializeMD(platform_name='Reference')
    system.simulation.context.setPositions(_coords(dna))
    with pytest.raises(Exception):
        system.simulation.context.getParameter('lambda_target')


def test_create_mutant_returns_a_new_dna():
    """DNA.create_mutant returns a new DNA with the target identity, leaving the original untouched."""
    xyz, dna_type = STRUCTURES[1]
    dna = _load(xyz, dna_type, 'mut')
    chainID, resSeq, xbase = _pick_base(dna)
    target = _TARGETS[xbase][0]
    mutant = dna.create_mutant([(chainID, resSeq, target)])
    assert mutant is not dna
    picked = (mutant.atoms['chainID'] == chainID) & (mutant.atoms['resSeq'] == resSeq) \
        & mutant.atoms['name'].isin(['A', 'T', 'G', 'C'])
    assert mutant.atoms.loc[picked, 'name'].iloc[0] == target
    original = (dna.atoms['chainID'] == chainID) & (dna.atoms['resSeq'] == resSeq) \
        & dna.atoms['name'].isin(['A', 'T', 'G', 'C'])
    assert dna.atoms.loc[original, 'name'].iloc[0] == xbase


def test_invalid_target_raises():
    xyz, dna_type = STRUCTURES[1]
    dna = _load(xyz, dna_type, 'bad')
    chainID, resSeq, _ = _pick_base(dna)
    with pytest.raises(ValueError):
        dna.create_mutant([(chainID, resSeq, 'X')])


@pytest.mark.parametrize("xyz, dna_type", _CASES, ids=_IDS)
def test_single_mutation_endpoints(xyz, dna_type):
    """lambda=0 reproduces the initial structure and lambda=1 the target, per force group, for a
    transition and both transversions."""
    ref = _load(xyz, dna_type, 'ref')
    coords = _coords(ref)
    chainID, resSeq, xbase = _pick_base(ref)

    for target in _TARGETS[xbase]:
        mutations = [(chainID, resSeq, target)]
        native_x = _group_energies(_native(xyz, dna_type, None, 'nx'), coords)
        native_y = _group_energies(_native(xyz, dna_type, mutations, 'ny'), coords)

        transformation, system = _alchemical(xyz, dna_type, mutations, 'alch')
        e0 = _group_energies_at(transformation, system, coords, 0.0)
        e1 = _group_energies_at(transformation, system, coords, 1.0)

        for name in FORCE_GROUPS:
            assert e0[name] == pytest.approx(native_x[name], abs=1e-4), \
                f"{dna_type} {xbase}->{target} lambda=0 mismatch in {name}"
            assert e1[name] == pytest.approx(native_y[name], abs=1e-4), \
                f"{dna_type} {xbase}->{target} lambda=1 mismatch in {name}"


@pytest.mark.parametrize("xyz, dna_type", _CASES, ids=_IDS)
def test_base_pair_comutation_endpoints(xyz, dna_type):
    """Mutating both strands of a Watson-Crick pair keeps pairing intact and matches native endpoints."""
    ref = _load(xyz, dna_type, 'ref')
    coords = _coords(ref)
    bases = ref.atoms[ref.atoms['name'].isin(['A', 'T', 'G', 'C'])]
    first, last = bases['chainID'].min(), bases['chainID'].max()
    if first == last:
        pytest.skip("single-strand structure")
    chainID, resSeq, xbase = _pick_base(ref)
    partner_res = int(bases[bases['chainID'] == last]['resSeq'].iloc[4])
    target = _TARGETS[xbase][0]
    mutations = [(chainID, resSeq, target), (last, partner_res, _COMPLEMENT[target])]

    native_x = _group_energies(_native(xyz, dna_type, None, 'nx'), coords)
    native_y = _group_energies(_native(xyz, dna_type, mutations, 'ny'), coords)
    transformation, system = _alchemical(xyz, dna_type, mutations, 'alch')
    e0 = _group_energies_at(transformation, system, coords, 0.0)
    e1 = _group_energies_at(transformation, system, coords, 1.0)
    for name in FORCE_GROUPS:
        assert e0[name] == pytest.approx(native_x[name], abs=1e-4)
        assert e1[name] == pytest.approx(native_y[name], abs=1e-4)


@pytest.mark.parametrize("xyz, dna_type", _CASES, ids=_IDS)
def test_energy_is_linear_in_lambda(xyz, dna_type):
    ref = _load(xyz, dna_type, 'ref')
    coords = _coords(ref)
    chainID, resSeq, xbase = _pick_base(ref)
    transformation, system = _alchemical(xyz, dna_type, [(chainID, resSeq, _TARGETS[xbase][1])], 'alch')
    system.simulation.context.setPositions(coords)
    ctx = system.simulation.context

    transformation.set_lambda(ctx, 0.0); u0 = _total(system)
    transformation.set_lambda(ctx, 1.0); u1 = _total(system)
    for lam in [0.1, 0.25, 0.5, 0.75, 0.9]:
        transformation.set_lambda(ctx, lam)
        assert _total(system) == pytest.approx((1 - lam) * u0 + lam * u1, abs=1e-4)


@pytest.mark.parametrize("xyz, dna_type", _CASES, ids=_IDS)
def test_ti_derivative(xyz, dna_type):
    """energy_derivative equals U(1)-U(0) and matches a central finite difference at several lambda."""
    ref = _load(xyz, dna_type, 'ref')
    coords = _coords(ref)
    chainID, resSeq, xbase = _pick_base(ref)
    transformation, system = _alchemical(xyz, dna_type, [(chainID, resSeq, _TARGETS[xbase][0])], 'alch')
    system.simulation.context.setPositions(coords)
    ctx = system.simulation.context

    def total(lam):
        transformation.set_lambda(ctx, lam)
        return _total(system)

    analytic = transformation.energy_derivative(ctx)
    assert analytic == pytest.approx(total(1.0) - total(0.0), abs=1e-6)
    d = 1e-4
    for lam in [0.2, 0.5, 0.8]:
        fd = (total(lam + d) - total(lam - d)) / (2 * d)
        assert fd == pytest.approx(analytic, abs=1e-2)


def test_from_sequence_bcurved_endpoints():
    """The primary build path (fromSequence, B_curved with template geometry) reproduces native
    endpoints, including the sequence-dependent equilibrium geometry rebuilt for the target."""
    seq = 'ATACAAAGGTGCGAGGTTTCTATGCTCCCACG'
    dna_x = DNA.fromSequence(seq, dna_type='B_curved', output_pdb=str(_outdir / 'seqx.pdb'),
                             atomistic_template=str(_outdir / 'atomistic_dna'))
    coords = _coords(dna_x)
    chainID, resSeq, xbase = _pick_base(dna_x, position=6)
    mutations = [(chainID, resSeq, _TARGETS[xbase][0])]

    dna_y = dna_x.create_mutant(mutations)
    dna_y.writePDB(str(_outdir / 'seqy.pdb'))
    native_x = _group_energies(dna_x, coords)
    native_y = _group_energies(dna_y, coords)

    dna_a = DNA.fromSequence(seq, dna_type='B_curved', output_pdb=str(_outdir / 'seqa.pdb'),
                             atomistic_template=str(_outdir / 'atomistic_dna'))
    transformation = AlchemicalTransformation(dna_a, mutations)
    system = System(dna_a, periodicBox=None)
    transformation.add_forces(system)
    system.initializeMD(platform_name='Reference')
    e0 = _group_energies_at(transformation, system, coords, 0.0)
    e1 = _group_energies_at(transformation, system, coords, 1.0)
    for name in FORCE_GROUPS:
        assert e0[name] == pytest.approx(native_x[name], abs=1e-4)
        assert e1[name] == pytest.approx(native_y[name], abs=1e-4)


def test_builds_on_cpu_platform():
    """All CustomNonbondedForce objects must share identical exclusions; the CPU platform enforces
    this (the Reference platform does not). The in-expression base-pair mask lets the two Exclusion
    copies share one exclusion list, so the hybrid must build a Context on CPU."""
    import openmm
    try:
        openmm.Platform.getPlatformByName('CPU')
    except Exception:
        pytest.skip("CPU platform unavailable")
    xyz, dna_type = STRUCTURES[1]
    ref = _load(xyz, dna_type, 'ref')
    coords = _coords(ref)
    chainID, resSeq, xbase = _pick_base(ref)
    dna = _load(xyz, dna_type, 'cpu')
    transformation = AlchemicalTransformation(dna, [(chainID, resSeq, _TARGETS[xbase][0])])
    system = System(dna, periodicBox=None)
    transformation.add_forces(system)
    system.initializeMD(platform_name='CPU')  # raises if exclusions are inconsistent
    system.simulation.context.setPositions(coords)
    for lam in (0.0, 0.5, 1.0):
        transformation.set_lambda(system.simulation.context, lam)
        assert np.isfinite(_total(system))


def test_md_runs_without_nan_forces():
    """A short simulation at fixed lambda produces no undefined forces."""
    xyz, dna_type = STRUCTURES[1]
    ref = _load(xyz, dna_type, 'ref')
    coords = _coords(ref)
    chainID, resSeq, xbase = _pick_base(ref)
    transformation, system = _alchemical(xyz, dna_type, [(chainID, resSeq, _TARGETS[xbase][0])], 'alch')
    system.simulation.context.setPositions(coords)
    transformation.set_lambda(system.simulation.context, 0.5)
    force_unit = unit.kilojoule_per_mole / unit.nanometer
    nan_particles = 0
    for _ in range(10):
        forces_state = system.simulation.context.getState(getForces=True).getForces(asNumpy=True)
        magnitudes = (np.asarray(forces_state.value_in_unit(force_unit)) ** 2).sum(axis=1) ** 0.5
        nan_particles += int(np.isnan(magnitudes).sum())
        system.simulation.step(1)
    assert nan_particles == 0
