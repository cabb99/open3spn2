"""Collective-variable readouts for protein-DNA systems.

Each class reports a geometric quantity (a distance, a position, a base-pair
offset) as an OpenMM force whose "energy" is that quantity -- there is no biasing
of the dynamics. They are meant to be read per force group (or embedded in a
``CustomCVForce``) for analysis or umbrella sampling, and like the biases in
:mod:`open3SPN2.force.bias` they are opt-in, never auto-added.

Because the reported value is carried as an identity energy, these forces use
``force_groups.MEASUREMENT`` (outside ``force_groups.TOTAL_ENERGY``) so they do
not contaminate the total potential energy.
"""

import logging

import openmm
import openmm.unit as unit

from .template import DNAForce, ProteinDNAForce
from .dna import select_residue_atom_indices
from .protein_dna import select_string_groups
from . import force_groups

logger = logging.getLogger(__name__)


class DistanceFromPoint(DNAForce):
    """Report the distance from the centroid of selected atoms to a point.

    Reports ``|centroid - r0|`` (Euclidean distance, nm), where the centroid is
    the unweighted mean of the selected atoms and ``r0 = (x0, y0, z0)``. Readout
    companion of :class:`open3SPN2.force.bias.PositionRestraint`.

    Parameters
    ----------
    appliedToResidues : sequence of int or None
        ``resSeq`` values of the residues to include; ``None`` uses every atom.
    """

    def __init__(self, dna, x0=1.0 * unit.nanometer, y0=1.0 * unit.nanometer,
                 z0=1.0 * unit.nanometer, appliedToResidues=None,
                 force_group=force_groups.MEASUREMENT):
        self.force_group = force_group
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.appliedToResidues = appliedToResidues
        super().__init__(dna)

    def reset(self):
        x0 = self.x0.value_in_unit(unit.nanometer)
        y0 = self.y0.value_in_unit(unit.nanometer)
        z0 = self.z0.value_in_unit(unit.nanometer)

        cx = openmm.CustomExternalForce("weight * x")
        cy = openmm.CustomExternalForce("weight * y")
        cz = openmm.CustomExternalForce("weight * z")
        for cv in (cx, cy, cz):
            cv.addPerParticleParameter("weight")

        distance = openmm.CustomCVForce(
            f"sqrt((cx-({x0}))^2 + (cy-({y0}))^2 + (cz-({z0}))^2)"
        )
        distance.addCollectiveVariable("cx", cx)
        distance.addCollectiveVariable("cy", cy)
        distance.addCollectiveVariable("cz", cz)
        distance.setForceGroup(self.force_group)
        self.force = distance

    def defineInteraction(self):
        indices = select_residue_atom_indices(self.dna, self.appliedToResidues)
        if len(indices) == 0:
            raise ValueError("No atoms selected for DistanceFromPoint; check appliedToResidues.")
        weight = 1.0 / len(indices)
        for cv in range(3):
            collective_variable = self.force.getCollectiveVariable(cv)
            for i in indices:
                collective_variable.addParticle(i, [weight])


class StringLength(ProteinDNAForce):
    """Report the centroid-centroid distance of a protein-DNA string.

    Reports ``distance(g1, g2)`` between the protein ``CA`` centroid and the DNA
    ``S`` centroid. Readout companion of
    :class:`open3SPN2.force.bias.StringProteinDNA`; the chain arguments accept one
    or many chains.
    """

    def __init__(self, dna, protein, chain_protein='A', chain_DNA='B',
                 protein_seg=False, group=None, force_group=force_groups.MEASUREMENT):
        self.chain_protein = chain_protein
        self.chain_DNA = chain_DNA
        self.protein_seg = protein_seg
        self.group = group or []
        self.force_group = force_group
        super().__init__(dna, protein)

    def reset(self):
        force = openmm.CustomCentroidBondForce(2, "distance(g1,g2)")
        force.setForceGroup(self.force_group)
        self.force = force

    def defineInteraction(self):
        protein_index, dna_index = select_string_groups(
            self.dna, self.chain_protein, self.chain_DNA, self.protein_seg, self.group)
        self.force.addGroup(protein_index)
        self.force.addGroup(dna_index)
        self.force.addBond([0, 1])


class BasePairProteinPosition(ProteinDNAForce):
    """Report a protein point's position along a DNA stretch, in base pairs.

    Reports the base-pair offset ``i`` of the protein point's projection onto the
    left->right base-pair vector (``i = 0`` at the midpoint). Readout companion of
    :class:`open3SPN2.force.bias.BasePairProteinHarmonicBias`.
    """

    def __init__(self, dna, protein, base_pair_left_indicies, base_pair_right_indicies,
                 protein_indices, base_pair_sep=4, force_group=force_groups.MEASUREMENT):
        self.base_pair_left_indicies = base_pair_left_indicies
        self.base_pair_right_indicies = base_pair_right_indicies
        self.protein_indices = protein_indices
        self.base_pair_sep = base_pair_sep
        self.force_group = force_group
        super().__init__(dna, protein)

    def reset(self):
        # Reported value: i = base-pair offset of the protein projection.
        energy = "i;"
        ratio_bp = f"i={self.base_pair_sep}*(ratio-0.5);"
        ratio_dots = ("ratio=VAVDdots/VDVDdots;"
                      "VAVDdots=VAx*VDx+VAy*VDy+VAz*VDz;"
                      "VDVDdots=VDx*VDx+VDy*VDy+VDz*VDz;")
        vectors = ("VAx=Px-BPLx;VAy=Py-BPLy;VAz=Pz-BPLz;"
                   "VDx=BPRx-BPLx;VDy=BPRy-BPLy;VDz=BPRz-BPLz;")
        particles = "BPLx=x1;BPLy=y1;BPLz=z1;BPRx=x2;BPRy=y2;BPRz=z2;Px=x3;Py=y3;Pz=z3"
        expression = f"{energy}{ratio_bp}{ratio_dots}{vectors}{particles}"

        force = openmm.CustomCentroidBondForce(3, expression)
        force.addGroup(self.base_pair_left_indicies)
        force.addGroup(self.base_pair_right_indicies)
        force.addGroup(self.protein_indices)
        force.addBond([0, 1, 2])
        force.setForceGroup(self.force_group)
        self.force = force

    def defineInteraction(self):
        pass
