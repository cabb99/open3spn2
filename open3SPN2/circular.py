"""Circular DNA (plasmids, minicircles) as an add-on to the linear force field.

A circular chain is a chain whose last residue is bonded to its first one. The core of open3SPN2 builds the
topology of linear chains and reads the connectivity of every force from that topology. Circular DNA only patches
the topology: :func:`circularize` registers :class:`ClosedChains` in ``DNA.topology_patches``, which adds the bonds,
angles, stackings and dihedrals that cross the closure each time ``DNA.computeTopology`` runs. The sequence
separation masks and the cross stacking of the forces then follow the closed topology.

* :func:`circularize` closes chains of any DNA, also a circular structure that was loaded as linear.
* :func:`bend_to_circle` and :class:`CircularBend` map a duplex onto a planar circle, adding or removing turns of the
  helix so the two ends join with a chosen linking number (relaxed, under- or overwound).
* :func:`make_periodic_helix` and :class:`PeriodicHelix` map a duplex onto a straight helix whose strands join their
  images through a periodic box, an infinite DNA.
* :func:`from_sequence` builds a closed duplex from a sequence, a circle or an infinite helix (used by
  ``DNA.fromSequence(circular=True)`` and ``DNA.fromSequence(infinite=True)``).
* :func:`linking_number` measures the linking number of two closed curves, e.g. the two strands of a plasmid along
  a trajectory.
"""
import numpy as np
import pandas

from .parser import fixPDB, pdb2table, writePDB

_dna_residues = ['DA', 'DC', 'DT', 'DG']


def circular_chain_info(atoms, circular_chains):
    """Returns ``{chainID: (first resSeq, number of residues)}`` for every circular chain.

    The residues of a circular chain must be numbered consecutively."""
    info = {}
    for chain in circular_chains:
        residues = np.unique(atoms.loc[atoms['chainID'] == chain, 'resSeq'].astype(int))
        if len(residues) == 0:
            raise ValueError(f'Circular chain {chain} is not in the structure')
        start, length = int(residues.min()), int(residues.max() - residues.min() + 1)
        if len(residues) != length:
            raise ValueError(f'The residues of the circular chain {chain} must be numbered consecutively')
        info[chain] = (start, length)
    return info


def dna_chains(dna):
    """DNA chain identifiers in the order they appear in the structure"""
    return list(dna.atoms.loc[dna.atoms['resname'].isin(_dna_residues), 'chainID'].unique())


def circularize(dna, chains=True, join=None, atomistic=None, output_pdb='clean.pdb', compute_topology=None):
    """Bonds the last residue of each circular chain to its first residue (e.g. a plasmid).

    The bonds, angles, stackings and dihedrals that cross the end of a circular chain are added to the topology
    by a topology patch (:class:`ClosedChains`), so they are also present after computeTopology runs again.

    Parameters
    ----------
    chains : True, str or list
        Chains to close. True closes every DNA chain (both strands of a plasmid).
    join : list of lists, optional
        Groups of chains joined end to end into a single chain before closing them. For example
        ``[['A', 'B']]`` bonds the 3' end of strand A to the 5' end of strand B, so the closed molecule is a
        single strand that runs along both sides of the duplex. The joined chain keeps the identifier of the
        first chain of the group.
    atomistic : str or DataFrame, optional
        The atomistic structure the DNA was coarse grained from (a pdb file or an atom table). Coarse graining
        a linear chain drops the phosphate of its first residue; with the atomistic structure it is rebuilt as
        the phosphate group that bridges the last and first residues, as every other phosphate bead.
        Otherwise it is placed by superposing a phosphate step of the chain.
    output_pdb : str or None
        Coarse grained pdb file to write the structure to (it is the topology read by System).
    compute_topology : bool, optional
        Recompute the bonded interactions. By default they are recomputed only if they had already been
        computed."""
    if join and atomistic is not None:
        raise ValueError('The phosphates of joined chains can not be taken from the atomistic structure')
    renamed = {}
    for group in join or []:
        join_chains(dna, group)
        renamed.update({c: group[0] for c in group})
    if chains is True:
        chains = dna_chains(dna)
    elif not chains:
        chains = []
    elif isinstance(chains, str):
        chains = [chains]
    chains = list(dict.fromkeys(renamed.get(c, c) for c in chains))

    if atomistic is not None and not isinstance(atomistic, pandas.DataFrame):
        atomistic = pdb2table(fixPDB(str(atomistic), add_hydrogens=False))
    minimum = 2 * (_interaction_extent(dna) + 2)
    for chain, (start, length) in circular_chain_info(dna.atoms, chains).items():
        if length < minimum:
            raise ValueError(f'Circular chain {chain} has {length} residues, at least {minimum} are needed')
        position = None if atomistic is None else _closing_phosphate(type(dna), atomistic, chain, start, length)
        add_closing_phosphate(dna, chain, start, length, position)
    dna.atoms.index = range(len(dna.atoms))
    dna.circular_chains = chains

    for chain, (start, length) in circular_chain_info(dna.atoms, chains).items():
        in_chain = dna.atoms['chainID'] == chain
        last_sugar = dna.atoms[in_chain & (dna.atoms['resSeq'] == start + length - 1) & (dna.atoms['name'] == 'S')]
        first_phosphate = dna.atoms[in_chain & (dna.atoms['resSeq'] == start) & (dna.atoms['name'] == 'P')]
        separation = last_sugar[['x', 'y', 'z']].to_numpy(dtype=float)[0] - first_phosphate[['x', 'y', 'z']].to_numpy(dtype=float)[0]
        if getattr(dna, 'periodic_box', None) is not None:
            box = np.asarray(dna.periodic_box, dtype=float) * 10
            separation -= np.round(separation / box) * box
        gap = np.linalg.norm(separation) / 10
        if gap > 1:
            print(f'Warning: the ends of the circular chain {chain} are {gap:.1f} nm apart. '
                  f'Bend the DNA (DNA.bendToCircle) or minimize the energy before running dynamics.')

    patches = [patch for patch in dna.topology_patches if not isinstance(patch, ClosedChains)]
    dna.topology_patches = patches + ([ClosedChains(chains)] if chains else [])
    if compute_topology is None:
        compute_topology = hasattr(dna, 'bonds')
    if compute_topology:
        dna.computeTopology(template_from_X3DNA=getattr(dna, 'template_from_X3DNA', True))
    if output_pdb:
        dna.writePDB(output_pdb)
    return dna


class ClosedChains:
    """Topology patch (see ``DNA.topology_patches``) that adds the interactions across the closure of `chains`."""

    def __init__(self, chains):
        self.chains = list(chains)

    def __call__(self, dna):
        close_topology(dna, self.chains)

    def __repr__(self):
        return f'ClosedChains({self.chains})'


def close_topology(dna, chains):
    """Adds to the bonds, angles, stackings and dihedrals of `dna` the interactions that close each of the
    `chains`: those that cross the closure and those of the phosphate of the first residue, which belongs to the
    closing step.

    The interactions come from the linear topology of a junction window: the last and first residues of the
    chain, numbered consecutively across the closure, built with the same DNA type, configuration and template
    mode. Their equilibrium values therefore come from the same local sequence (or coordinates) as any other
    interaction of the chain."""
    margin = _interaction_extent(dna) + 2
    tables = {'bonds': ['aai', 'aaj'], 'angles': ['aai', 'aaj', 'aak', 'aax'], 'stackings': ['aai', 'aaj', 'aak'],
              'dihedrals': ['aai', 'aaj', 'aak', 'aal']}
    for chain, (start, length) in circular_chain_info(dna.atoms, chains).items():
        last = start + length - 1
        in_chain = dna.atoms['chainID'] == chain
        tail = dna.atoms[in_chain & (dna.atoms['resSeq'] > last - margin)]
        head = dna.atoms[in_chain & (dna.atoms['resSeq'] < start + margin)]
        window_atoms = pandas.concat([tail.assign(resSeq=tail['resSeq'] - (last - margin)),
                                      head.assign(resSeq=head['resSeq'] - start + margin + 1)])
        origin = window_atoms.index.to_numpy()
        side = np.array([0] * len(tail) + [1] * len(head))
        phosphate = dna.atoms.index[in_chain & (dna.atoms['resSeq'] == start) & (dna.atoms['name'] == 'P')]
        window_phosphate = np.flatnonzero(np.isin(origin, phosphate))

        window = type(dna).fromAtoms(window_atoms.reset_index(drop=True), dna.DNAtype, output_pdb=None,
                                     compute_topology=False)
        for attr in ['config', 'particle_definition', 'bond_definition', 'angle_definition', 'dihedral_definition',
                     'stacking_definition', 'pair_definition', 'cross_definition']:
            if hasattr(dna, attr):
                setattr(window, attr, getattr(dna, attr))
        window.computeTopology(template_from_X3DNA=getattr(dna, 'template_from_X3DNA', True))

        for name, columns in tables.items():
            linear = getattr(dna, name)
            linear = linear[~linear[columns].isin(phosphate).any(axis=1)]
            table = getattr(window, name)
            atoms = table[columns].to_numpy(dtype=int)
            closes = ((side[atoms].min(axis=1) == 0) & (side[atoms].max(axis=1) == 1)) | \
                np.isin(atoms, window_phosphate).any(axis=1)
            closing = table[closes].copy()
            closing[columns] = origin[atoms[closes]]
            setattr(dna, name, pandas.concat([linear, closing], ignore_index=True))


def _interaction_extent(dna):
    """Largest number of residues between two atoms of a bonded interaction of the configuration"""
    if not hasattr(dna, 'bond_definition'):
        dna.parseConfigurationFile()
    extent = 0
    for definition, offsets in [(dna.bond_definition, ['s1']), (dna.angle_definition, ['s1', 's2', 'sB']),
                                (dna.stacking_definition, ['s1', 's2']), (dna.dihedral_definition, ['s1', 's2', 's3'])]:
        values = definition[offsets].astype(int).to_numpy()
        extent = max(extent, int((np.maximum(values.max(axis=1), 0) - np.minimum(values.min(axis=1), 0)).max()))
    return extent


def join_chains(dna, group):
    """Joins the chains in `group` end to end into a single chain named as the first one."""
    atoms = dna.atoms
    chain_order = list(atoms['chainID'].unique())
    if len(set(group)) != len(group) or any(c not in chain_order for c in group):
        raise ValueError(f'Can not join the chains {group}, the structure has the chains {chain_order}')
    rank = {c: float(i) for i, c in enumerate(chain_order)}
    for i, c in enumerate(group[1:], 1):
        rank[c] = rank[group[0]] + i / len(group)
    atoms['_rank'] = atoms['chainID'].map(rank)
    end = atoms.loc[atoms['chainID'] == group[0], 'resSeq'].max()
    for c in group[1:]:
        in_chain = atoms['chainID'] == c
        atoms.loc[in_chain, 'resSeq'] += end + 1 - atoms.loc[in_chain, 'resSeq'].min()
        atoms.loc[in_chain, 'chainID'] = group[0]
        end = atoms.loc[in_chain, 'resSeq'].max()
    dna.atoms = atoms.sort_values('_rank', kind='stable').drop(columns='_rank').reset_index(drop=True)


def _closing_phosphate(dna_class, atomistic, chain, start, length):
    """Position of the phosphate bead of the first residue of a closed chain, coarse grained from the atomistic
    structure: the chain is reordered to start at its last residue, so the O3' of the last residue joins the
    phosphate group of the first one."""
    last = start + length - 1
    table = atomistic.copy()
    in_chain = table['chainID'] == chain
    table.loc[in_chain & (table['resSeq'] == last), 'resSeq'] = start - 1
    table['_chain_order'] = table['chainID'].map({c: i for i, c in enumerate(table['chainID'].unique())})
    table = table.sort_values(['_chain_order', 'resSeq'], kind='stable').drop(columns='_chain_order')
    beads = dna_class.CoarseGrain(table.reset_index(drop=True))
    bead = beads[(beads['chainID'] == chain) & (beads['resSeq'] == start) & (beads['name'] == 'P')]
    return bead[['x', 'y', 'z']].to_numpy(dtype=float)[0]


def add_closing_phosphate(dna, chain, start, length, position=None):
    """Adds the phosphate bead of the first residue of a circular chain if it is missing.

    It is placed at `position` if given; otherwise the phosphate of a residue near the middle of the chain is
    copied, superposing the two sugars and bases it connects onto the last and first residues of the chain."""
    atoms = dna.atoms
    xyz = atoms[['x', 'y', 'z']].values
    in_chain = atoms['chainID'] == chain
    is_base = atoms['name'].isin(['A', 'T', 'G', 'C', 'U'])

    def bead(resSeq, name):
        sel = in_chain & (atoms['resSeq'] == resSeq) & (is_base if name == 'B' else atoms['name'] == name)
        return np.flatnonzero(sel.values)[0] if sel.any() else None

    first, last = start, start + length - 1
    if bead(first, 'P') is not None:
        return
    target_ids = [bead(last, 'S'), bead(last, 'B'), bead(first, 'S'), bead(first, 'B')]
    if None in target_ids:
        raise ValueError(f'The first and last residues of the circular chain {chain} need a sugar and a base')
    target = xyz[target_ids]

    if position is None:
        position = target[[0, 2]].mean(axis=0)
        for r in sorted(range(first + 1, last + 1), key=lambda r: abs(r - first - length / 2)):
            ids = [bead(r - 1, 'S'), bead(r - 1, 'B'), bead(r, 'S'), bead(r, 'B'), bead(r, 'P')]
            if None not in ids:
                rotation, translation = kabsch(xyz[ids[:4]], target)
                rmsd = np.sqrt(np.mean(np.sum((xyz[ids[:4]] @ rotation.T + translation - target) ** 2, axis=1)))
                if rmsd < 2.0:  # Angstrom, otherwise the ends do not look like a dinucleotide step
                    position = xyz[ids[4]] @ rotation.T + translation
                break

    phosphate = atoms.iloc[[target_ids[2]]].copy()
    phosphate['name'] = 'P'
    if 'type' in phosphate:
        phosphate['type'] = 'P'
    phosphate['element'] = 'P'
    phosphate[['x', 'y', 'z']] = position
    position_in_table = np.flatnonzero((in_chain & (atoms['resSeq'] == first)).values)[0]
    dna.atoms = pandas.concat([atoms.iloc[:position_in_table], phosphate,
                               atoms.iloc[position_in_table:]]).reset_index(drop=True)


def duplex_base_pairs(dna):
    """Base pairs of a two stranded DNA, ordered along the first strand.

    Returns the two chain identifiers, the index of the base pair of every atom of the duplex (-1 for other
    atoms) and four (n, 3) arrays with the sugar and base of the first strand and the sugar and base of the
    second strand of each base pair."""
    chains = dna_chains(dna)
    if len(chains) != 2:
        raise ValueError(f'A double stranded DNA with two chains is needed, the DNA has the chains {chains}')
    a, b = chains
    residues_a = sorted(dna.atoms.loc[dna.atoms['chainID'] == a, 'resSeq'].unique())
    residues_b = sorted(dna.atoms.loc[dna.atoms['chainID'] == b, 'resSeq'].unique())[::-1]
    if len(residues_a) != len(residues_b):
        raise ValueError(f'The chains {a} and {b} have a different number of residues')
    pair_index = {(a, r): i for i, r in enumerate(residues_a)}
    pair_index.update({(b, r): i for i, r in enumerate(residues_b)})
    base_pair = np.array([pair_index.get((c, r), -1) for c, r in zip(dna.atoms['chainID'], dna.atoms['resSeq'])])

    lookup = {}
    for c, r, n, x, y, z in zip(dna.atoms['chainID'], dna.atoms['resSeq'], dna.atoms['name'],
                                dna.atoms['x'], dna.atoms['y'], dna.atoms['z']):
        lookup[(c, r, n if n in ['P', 'S'] else 'B')] = (x, y, z)
    beads = [np.array([lookup[(c, r, n)] for r in residues]) for c, residues in [(a, residues_a), (b, residues_b)]
             for n in ['S', 'B']]
    return (a, b), base_pair, beads


def bend_to_circle(dna, delta_linking_number=0, superhelical_density=None, closed=False, first=0, last=None):
    """Bends a double stranded DNA into a planar circle so its ends can be joined.

    The circle lies on the xy plane, centered at the origin. The base pairs are placed rigidly along the circle
    keeping their rise and twist, and the twist is adjusted so the circle has an integer linking number
    Lk = Lk0 + delta_linking_number, where Lk0 is the number of helical turns of the relaxed duplex. A positive
    delta_linking_number overwinds the DNA (positive supercoiling) and a negative one underwinds it (negative
    supercoiling, like most plasmids in cells). In the planar circle all the excess goes to twist, and it converts
    partially to writhe (plectonemes) during the simulation.

    Parameters
    ----------
    delta_linking_number : int
        Turns added to the relaxed linking number.
    superhelical_density : float, optional
        Alternative to delta_linking_number, sigma = delta_Lk / Lk0 (e.g. -0.06).
    closed, first, last :
        See CircularBend. Used by from_sequence to close the circle exactly.

    Returns the CircularBend object, with the radius (Angstrom) and linking numbers; it is also stored in
    ``dna.bend``."""
    (a, b), base_pair, (sugar_a, base_a, sugar_b, base_b) = duplex_base_pairs(dna)
    centers = (sugar_a + base_a + sugar_b + base_b) / 4
    bend = CircularBend(centers, sugar_b - sugar_a, closed=closed, delta_linking_number=delta_linking_number,
                        superhelical_density=superhelical_density, first=first, last=last)
    duplex = base_pair >= 0
    dna.atoms.loc[duplex, ['x', 'y', 'z']] = bend.apply(dna.atoms.loc[duplex, ['x', 'y', 'z']].values,
                                                       base_pair[duplex])
    dna.bend = bend
    print(f'DNA bent into a circle of radius {bend.radius / 10:.2f} nm with linking number {bend.linking_number} '
          f'(relaxed {bend.relaxed_linking_number}, delta {bend.delta_linking_number:+d})')
    return bend


def make_periodic_helix(dna, delta_linking_number=0, superhelical_density=None, closed=False, first=0, last=None):
    """Maps a double stranded DNA onto a straight helix along z that repeats every period, so its ends can be joined
    through a periodic box (an infinite DNA).

    The base pairs keep their rise and twist, and the twist is adjusted so each strand winds an integer number of
    turns per period, turns = relaxed turns + delta_linking_number (see PeriodicHelix). Returns the PeriodicHelix,
    with the period `length` (Angstrom) and the numbers of turns; it is also stored in ``dna.helix``."""
    (a, b), base_pair, (sugar_a, base_a, sugar_b, base_b) = duplex_base_pairs(dna)
    centers = (sugar_a + base_a + sugar_b + base_b) / 4
    helix = PeriodicHelix(centers, sugar_b - sugar_a, closed=closed, delta_linking_number=delta_linking_number,
                          superhelical_density=superhelical_density, first=first, last=last)
    duplex = base_pair >= 0
    dna.atoms.loc[duplex, ['x', 'y', 'z']] = helix.apply(dna.atoms.loc[duplex, ['x', 'y', 'z']].values,
                                                        base_pair[duplex])
    dna.helix = helix
    print(f'DNA straightened into a periodic helix of {helix.length / 10:.2f} nm with {helix.turns} turns '
          f'(relaxed {helix.relaxed_turns}, delta {helix.delta_turns:+d})')
    return helix


def trim_circle(atoms, first, length):
    """Keeps the base pairs first ... first + length - 1 (0 based, along the first strand) of a duplex, with the
    phosphates that join them to the previous base pairs, and renumbers both strands from 1."""
    a, b = atoms['chainID'].unique()[:2]
    residues_a = sorted(atoms.loc[atoms['chainID'] == a, 'resSeq'].unique())
    residues_b = sorted(atoms.loc[atoms['chainID'] == b, 'resSeq'].unique())[::-1]
    new_a = {residues_a[first + i]: i + 1 for i in range(length)}
    new_b = {residues_b[first + i]: length - i for i in range(length)}
    keep = (((atoms['chainID'] == a) & atoms['resSeq'].isin(new_a)) |
            ((atoms['chainID'] == b) & atoms['resSeq'].isin(new_b)))
    atoms = atoms[keep].copy()
    atoms['resSeq'] = [new_a[r] if c == a else new_b[r] for c, r in zip(atoms['chainID'], atoms['resSeq'])]
    atoms['_chain'] = (atoms['chainID'] == b).astype(int)
    atoms = atoms.sort_values(['_chain', 'resSeq'], kind='stable')
    return atoms.drop(columns=['_chain']).reset_index(drop=True)


def compute_linking_number(dna, positions=None):
    """Linking number of a circular duplex, from the sugar beads of its two strands.

    `positions` is an optional (n_atoms, 3) array (e.g. a frame of a trajectory) in the same units as the
    structure; by default the current coordinates are used. Lk is a topological invariant: it only changes if a
    strand breaks, so it is a good check that a plasmid keeps its topology during a simulation."""
    atoms = dna.atoms
    xyz = atoms[['x', 'y', 'z']].values if positions is None else np.asarray(positions)
    curves = []
    for chain in dna_chains(dna):
        sugars = atoms[(atoms['chainID'] == chain) & (atoms['name'] == 'S')].sort_values('resSeq')
        curves += [xyz[sugars.index]]
    if len(curves) != 2:
        raise ValueError('The linking number is defined for two strands')
    return linking_number(*curves)


def from_sequence(dna_class, sequence, dna_type='B_curved', output_pdb='clean.pdb', compute_topology=True,
                  atomistic_template=None, infinite=False, box_size=20, join=None, delta_linking_number=0,
                  superhelical_density=None, padding=6):
    """Builds a closed duplex of `dna_class` from a sequence: a planar circle, or with `infinite` a straight helix
    along z whose ends are bonded through a periodic box.

    A straight duplex is built with the sequence, one extra base pair that copies the first one, and `padding`
    base pairs at each end that continue the circular sequence, so the first and last base pairs keep the same
    neighbors they have in the closed molecule. It is bent into a circle (or straightened into a periodic helix) so
    the extra base pair falls on the first one (or on its periodic image), which closes it exactly, and the extra
    base pairs are removed. The infinite helix is centered in a box of `box_size` x `box_size` x period nm, stored in
    ``dna.periodic_box``. The atomistic structure is written to `atomistic_template` if given."""
    n = len(sequence)
    extended = [sequence[i % n] for i in range(-padding, n + padding + 1)]
    extended = pandas.Series(extended, index=[('A', i) for i in range(len(extended))])
    builder = dna_class()
    builder.parseConfigurationFile()
    atoms = builder.computeGeometry(extended)

    dna = dna_class.fromAtoms(atoms, dna_type, output_pdb=None, compute_topology=False)
    shape = make_periodic_helix if infinite else bend_to_circle
    mapping = shape(dna, delta_linking_number, superhelical_density, closed=True, first=padding, last=padding + n)
    dna.atoms = trim_circle(dna.atoms, padding, n)
    center = np.zeros(3)
    if infinite:
        dna.periodic_box = [box_size, box_size, mapping.length / 10]
        center = np.array([box_size, box_size, 0]) * 10 / 2
        dna.atoms[['x', 'y', 'z']] += center

    if atomistic_template:
        atomistic = builder.template.reset_index(drop=True)
        total = len(extended)
        base_pair = np.where(atomistic['chainID'] == 'A', atomistic['resSeq'] - 1, total - atomistic['resSeq'])
        atomistic[['x', 'y', 'z']] = mapping.apply(atomistic[['x', 'y', 'z']].values, base_pair) + center
        atomistic = trim_circle(atomistic, padding, n)
        atomistic['serial'] = atomistic.index + 1
        writePDB(atomistic, pdb_file=atomistic_template)

    circularize(dna, chains=True, join=join, output_pdb=None, compute_topology=False)
    if compute_topology:
        dna.computeTopology()
    if output_pdb:
        dna.writePDB(output_pdb)
    return dna


def kabsch(mobile, reference):
    """Optimal rotation ``R`` and translation ``t`` so that ``mobile @ R.T + t`` fits `reference`."""
    mobile = np.asarray(mobile, dtype=float)
    reference = np.asarray(reference, dtype=float)
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    h = (mobile - mobile_center).T @ (reference - reference_center)
    u, s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1, 1, d]) @ u.T
    translation = reference_center - mobile_center @ rotation.T
    return rotation, translation


class CircularBend:
    """Maps a DNA duplex onto a planar circle, one rigid base pair at a time.

    Each base pair gets a local frame: its origin and tangent follow the helical axis (a line fitted
    to the base pair centers over about one helical turn) and its first axis points from one sugar to
    the other. The base pairs are then placed rigidly along a circle lying on the xy plane, centered at the
    origin, keeping the rise and the twist between consecutive base pairs. Any intrinsic curvature of
    the duplex is replaced by the uniform curvature of the circle. An extra twist, distributed
    uniformly along the circle, makes the end of the duplex join its start with an integer linking
    number.

    Parameters
    ----------
    centers : (n, 3) array
        Center of each base pair, ordered along the duplex.
    references : (n, 3) array
        A vector fixed to each base pair, from the sugar of the first strand to the sugar of the
        second one, used to measure the twist.
    closed : bool
        If True, base pair `last` is a copy of base pair `first` translated one step along the helix,
        so the length and twist of the step that closes the circle are known exactly. If False, the
        closing step is the average step of the duplex.
    delta_linking_number : int
        Turns added to (positive, overwound) or removed from (negative, underwound) the relaxed
        linking number Lk0, the number of helical turns of the duplex rounded to an integer.
    superhelical_density : float, optional
        Alternative to `delta_linking_number`; the change is ``round(sigma * Lk0)``.
    first, last : int, optional
        Base pairs that form the circle (by default all of them). Base pairs outside this range are
        still mapped, continuing the circle, which is useful to pad the ends of the duplex so all the
        base pairs of the circle have a complete neighborhood to compute their frames.
    window : int
        Number of base pairs used to fit the helical axis.
    """

    def __init__(self, centers, references, closed=False, delta_linking_number=0, superhelical_density=None,
                 first=0, last=None, window=11):
        last = len(centers) - 1 if last is None else last
        self.origins, self.frames, self.contour, self.turns = _base_pair_frames(centers, references, first, last,
                                                                                window)
        self.circumference, self.twist = _period(self.contour, self.turns, first, last, closed)
        self.radius = self.circumference / (2 * np.pi)
        n = len(self.contour)

        self.relaxed_linking_number = int(np.round(self.twist / (2 * np.pi)))
        if superhelical_density is not None:
            delta_linking_number = int(np.round(superhelical_density * self.relaxed_linking_number))
        self.delta_linking_number = int(delta_linking_number)
        self.linking_number = self.relaxed_linking_number + self.delta_linking_number
        self.extra_twist = 2 * np.pi * self.linking_number - self.twist

        # Target frames on the circle; e1 -> z, e2 -> radial direction, t -> tangent of the circle when untwisted
        theta = 2 * np.pi * self.contour / self.circumference
        alpha = self.turns + self.extra_twist * self.contour / self.circumference
        z = np.array([0, 0, 1.])
        radial = np.stack([np.cos(theta), np.sin(theta), np.zeros(n)], axis=1)
        tangent = np.stack([-np.sin(theta), np.cos(theta), np.zeros(n)], axis=1)
        target_e1 = np.cos(alpha)[:, None] * z + np.sin(alpha)[:, None] * radial
        target_e2 = -np.sin(alpha)[:, None] * z + np.cos(alpha)[:, None] * radial
        self.target_origins = self.radius * radial
        self.target_frames = np.stack([target_e1, target_e2, tangent], axis=2)

    def apply(self, xyz, base_pair):
        """Returns the bent coordinates of `xyz` (same units as the base pair centers). `base_pair` is the
        index of the base pair each coordinate belongs to."""
        return _move_base_pairs(self, xyz, base_pair)


class PeriodicHelix:
    """Maps a DNA duplex onto a straight helix along z that repeats with period `length`, one rigid base pair at a
    time.

    The base pair frames, rise and twist are measured as in CircularBend. The base pairs are placed along the z axis,
    keeping their rise and twist; an extra twist, distributed uniformly, makes each strand wind an integer number of
    `turns` per period, so the strands join their images in a periodic box of that length (an infinite DNA). Any
    intrinsic curvature of the duplex is removed. The parameters are those of CircularBend, with the linking number
    of the circle replaced by the number of turns per period."""

    def __init__(self, centers, references, closed=False, delta_linking_number=0, superhelical_density=None,
                 first=0, last=None, window=11):
        last = len(centers) - 1 if last is None else last
        self.origins, self.frames, self.contour, self.turns_along = _base_pair_frames(centers, references, first,
                                                                                      last, window)
        self.length, self.twist = _period(self.contour, self.turns_along, first, last, closed)
        self.relaxed_turns = int(np.round(self.twist / (2 * np.pi)))
        if superhelical_density is not None:
            delta_linking_number = int(np.round(superhelical_density * self.relaxed_turns))
        self.delta_turns = int(delta_linking_number)
        self.turns = self.relaxed_turns + self.delta_turns
        self.extra_twist = 2 * np.pi * self.turns - self.twist

        alpha = self.turns_along + self.extra_twist * self.contour / self.length
        n = len(self.contour)
        axis = np.array([0, 0, 1.])
        target_e1 = np.stack([np.cos(alpha), np.sin(alpha), np.zeros(n)], axis=1)
        target_e2 = np.stack([-np.sin(alpha), np.cos(alpha), np.zeros(n)], axis=1)
        self.target_origins = self.contour[:, None] * axis
        self.target_frames = np.stack([target_e1, target_e2, np.tile(axis, (n, 1))], axis=2)

    def apply(self, xyz, base_pair):
        """Returns the straightened coordinates of `xyz` (same units as the base pair centers). `base_pair` is the
        index of the base pair each coordinate belongs to."""
        return _move_base_pairs(self, xyz, base_pair)


def _base_pair_frames(centers, references, first, last, window):
    """Local frames of the base pairs of a duplex and the contour length and accumulated twist along it, measured
    from base pair `first` (see CircularBend)."""
    centers = np.asarray(centers, dtype=float)
    references = np.asarray(references, dtype=float)
    n = len(centers)
    if last - first < 2:
        raise ValueError('At least three base pairs are needed to bend a duplex')

    # The helical axis is the line that best fits the centers of the `window` base pairs around each base pair
    # (shifted inwards at the ends of the duplex)
    window = min(window, n)
    origins = np.zeros((n, 3))
    tangents = np.zeros((n, 3))
    for i in range(n):
        start = min(max(i - window // 2, 0), n - window)
        points = centers[start:start + window]
        mean = points.mean(axis=0)
        direction = np.linalg.svd(points - mean)[2][0]
        if (points[-1] - points[0]) @ direction < 0:
            direction = -direction
        origins[i] = mean + ((centers[i] - mean) @ direction) * direction
        tangents[i] = direction
    e1 = references - np.sum(references * tangents, axis=1, keepdims=True) * tangents
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(tangents, e1)
    frames = np.stack([e1, e2, tangents], axis=2)  # columns e1, e2, t

    # Rise and twist of each step; the twist is the rotation of e1 around the axis after moving it along the axis
    # without twisting it (parallel transport)
    rise = np.linalg.norm(np.diff(origins, axis=0), axis=1)
    twist = np.zeros(n - 1)
    for i in range(n - 1):
        moved = _minimal_rotation(tangents[i], tangents[i + 1]) @ e1[i]
        twist[i] = np.arctan2(np.cross(moved, e1[i + 1]) @ tangents[i + 1], moved @ e1[i + 1])
    contour = np.concatenate([[0], np.cumsum(rise)])
    turns = np.concatenate([[0], np.cumsum(twist)])
    return origins, frames, contour - contour[first], turns - turns[first]


def _period(contour, turns, first, last, closed):
    """Length and total twist of one period of the base pairs first ... last. If `closed`, base pair `last` repeats
    base pair `first`; otherwise the closing step is the average step."""
    if closed:
        return contour[last], turns[last]
    scale = (last - first + 1) / (last - first)
    return contour[last] * scale, turns[last] * scale


def _move_base_pairs(mapping, xyz, base_pair):
    """Moves each coordinate rigidly with its base pair, from the mapping's frames to its target frames"""
    xyz = np.asarray(xyz, dtype=float)
    base_pair = np.asarray(base_pair, dtype=int)
    local = np.einsum('nji,nj->ni', mapping.frames[base_pair], xyz - mapping.origins[base_pair])
    return np.einsum('nij,nj->ni', mapping.target_frames[base_pair], local) + mapping.target_origins[base_pair]


def _minimal_rotation(a, b):
    """Rotation matrix that takes the unit vector `a` to the unit vector `b` around their common normal"""
    v = np.cross(a, b)
    c = a @ b
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx / (1 + c)


def linking_number(curve1, curve2):
    """Gauss linking number of two closed polygonal curves (``(n, 3)`` arrays, the last point is
    connected back to the first one).

    For the two strands of a circular duplex, e.g. the sugar beads of each strand, this is the
    linking number Lk = Tw + Wr of the plasmid, which is constant along an unbroken simulation. It is
    positive for a right handed double helix. Uses the exact solid angle of every pair of segments
    (Klenin & Langowski, Biopolymers 54, 307 (2000))."""
    p1 = np.asarray(curve1, dtype=float)
    p3 = np.asarray(curve2, dtype=float)
    p2 = np.roll(p1, -1, axis=0)
    p4 = np.roll(p3, -1, axis=0)
    p1, p2 = p1[:, None], p2[:, None]
    p3, p4 = p3[None, :], p4[None, :]

    r13, r14, r23, r24 = p3 - p1, p4 - p1, p3 - p2, p4 - p2

    def normal(a, b):
        n = np.cross(a, b)
        norm = np.linalg.norm(n, axis=-1, keepdims=True)
        return np.divide(n, norm, out=np.zeros_like(n), where=norm > 0)

    n1, n2, n3, n4 = normal(r13, r14), normal(r14, r24), normal(r24, r23), normal(r23, r13)

    def asin(x):
        return np.arcsin(np.clip(x, -1, 1))

    omega = (asin(np.sum(n1 * n2, -1)) + asin(np.sum(n2 * n3, -1)) +
             asin(np.sum(n3 * n4, -1)) + asin(np.sum(n4 * n1, -1)))
    sign = np.sign(np.sum(np.cross(p4 - p3, p2 - p1) * r13, -1))
    return float(-np.sum(omega * sign) / (4 * np.pi))
