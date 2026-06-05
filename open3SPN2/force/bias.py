"""Bias forces for protein-DNA systems.

These forces add energy to steer or restrain a simulation (position restraints,
string/umbrella pulling, harmonic biases on a collective variable). They are not
part of the physical force field and are never auto-added; a user opts in by
adding them explicitly.

The identity-energy *readout* versions (report a value without biasing the
system) live in :mod:`open3SPN2.force.collective_variables`.
"""

import logging

import numpy as np
import openmm
import openmm.unit as unit

from .template import DNAForce, ProteinDNAForce
from .dna import select_residue_atom_indices
from .protein_dna import (ElectrostaticsProteinDNA, select_string_groups,
                          _proteinResidues, _dnaResidues)
from . import force_groups

logger = logging.getLogger(__name__)


class PositionRestraint(DNAForce):
    """Harmonic restraint on the centroid of selected atoms.

    Energy ``0.5 * k * |centroid - r0|**2`` pulls the geometric centroid of the
    selected atoms toward ``r0 = (x0, y0, z0)``. The centroid is the unweighted
    mean ``sum(x_i)/N``, implemented as three ``CustomExternalForce("weight*x")``
    collective variables with ``weight = 1/N``.

    Parameters
    ----------
    k : openmm.unit.Quantity
        Force constant, in energy / length**2 (e.g. kJ/mol/nm**2).
    x0, y0, z0 : openmm.unit.Quantity
        Target position of the centroid.
    appliedToResidues : sequence of int or None
        ``resSeq`` values of the residues whose atoms are restrained; ``None``
        (default) uses the centroid of every atom.
    """

    def __init__(self, dna, k=1000 * unit.kilojoule_per_mole / unit.nanometer ** 2,
                 x0=1.0 * unit.nanometer, y0=1.0 * unit.nanometer, z0=1.0 * unit.nanometer,
                 appliedToResidues=None, force_group=force_groups.MEMBRANE, OpenCLPatch=True):
        self.force_group = force_group
        self.k = k
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.appliedToResidues = appliedToResidues
        super().__init__(dna, OpenCLPatch=OpenCLPatch)

    def reset(self):
        k = self.k.value_in_unit(unit.kilojoule_per_mole / unit.nanometer ** 2)
        x0 = self.x0.value_in_unit(unit.nanometer)
        y0 = self.y0.value_in_unit(unit.nanometer)
        z0 = self.z0.value_in_unit(unit.nanometer)

        cx = openmm.CustomExternalForce("weight * x")
        cy = openmm.CustomExternalForce("weight * y")
        cz = openmm.CustomExternalForce("weight * z")
        for cv in (cx, cy, cz):
            cv.addPerParticleParameter("weight")

        harmonic = openmm.CustomCVForce(
            f"0.5*{k}*((cx-({x0}))^2 + (cy-({y0}))^2 + (cz-({z0}))^2)"
        )
        harmonic.addCollectiveVariable("cx", cx)
        harmonic.addCollectiveVariable("cy", cy)
        harmonic.addCollectiveVariable("cz", cz)
        harmonic.setForceGroup(self.force_group)
        self.force = harmonic

    def defineInteraction(self):
        indices = select_residue_atom_indices(self.dna, self.appliedToResidues)
        if len(indices) == 0:
            raise ValueError("No atoms selected for PositionRestraint; check appliedToResidues.")
        weight = 1.0 / len(indices)
        for cv in range(3):
            collective_variable = self.force.getCollectiveVariable(cv)
            for i in indices:
                collective_variable.addParticle(i, [weight])


class StringProteinDNA(ProteinDNAForce):
    """Harmonic string/umbrella potential between a protein and a DNA segment.

    Energy ``0.5 * k_string_PD * (distance(g1, g2) - r0)**2`` where ``g1`` is the
    centroid of the protein ``CA`` atoms and ``g2`` the centroid of the DNA ``S``
    atoms. The chain arguments accept one or many chains (``'A'``, ``'AB'`` or
    ``['A', 'B']``).

    Parameters
    ----------
    r0 : float
        Equilibrium centroid-centroid distance, in nanometers.
    k_string_PD : float
        Force constant, in kJ/mol/nm**2.
    """

    def __init__(self, dna, protein, r0, chain_protein='A', chain_DNA='B',
                 k_string_PD=10 * 4.184, protein_seg=False, group=None,
                 force_group=force_groups.PULLING):
        self.k_string_PD = k_string_PD
        self.chain_protein = chain_protein
        self.chain_DNA = chain_DNA
        self.r0 = r0
        self.protein_seg = protein_seg
        self.group = group or []
        self.force_group = force_group
        super().__init__(dna, protein)

    def reset(self):
        stringForce = openmm.CustomCentroidBondForce(
            2, f"0.5*{self.k_string_PD}*(distance(g1,g2)-{self.r0})^2")
        stringForce.setForceGroup(self.force_group)
        self.force = stringForce
        logger.debug("StringProteinDNA: r0=%s, k_string_PD=%s", self.r0, self.k_string_PD)

    def defineInteraction(self):
        protein_index, dna_index = select_string_groups(
            self.dna, self.chain_protein, self.chain_DNA, self.protein_seg, self.group)
        self.force.addGroup(protein_index)
        self.force.addGroup(dna_index)
        self.force.addBond([0, 1])


class ElectrostaticBias(ProteinDNAForce):
    """Harmonic bias on the protein-DNA electrostatic energy.

    Adds ``0.5 * k_ebias * ((E_elec - center)/4.184)**2`` where ``E_elec`` is the
    protein-DNA Debye-Huckel electrostatic energy (from
    :class:`open3SPN2.force.protein_dna.ElectrostaticsProteinDNA`) used as a
    collective variable.

    Units: ``E_elec`` is reported by OpenMM in kJ/mol, so ``center`` is in kJ/mol.
    The ``/4.184`` converts the deviation ``(E_elec - center)`` into kcal/mol, so
    ``k_ebias`` has units of kJ/mol per (kcal/mol)**2.
    """

    def __init__(self, dna, protein, k_ebias, center, k_elec, ldby,
                 cutoff_distance=None, force_group=force_groups.PROTEIN_DNA_BIAS):
        self.k_ebias = k_ebias
        self.center = center
        self.k_elec = k_elec
        self.ldby = ldby
        self.force_group = force_group
        self.cutoff_distance = cutoff_distance
        super().__init__(dna, protein)

    def reset(self):
        ebiasForce = openmm.CustomCVForce("0.5*k_ebias*((E_elec-center)/4.184)^2")
        E_elec = ElectrostaticsProteinDNA(self.dna, self.protein, k=self.k_elec,
                                          ldby=self.ldby, cutoff_distance=self.cutoff_distance)
        ebiasForce.addCollectiveVariable("E_elec", E_elec.force)  # kJ/mol
        ebiasForce.addGlobalParameter("k_ebias", self.k_ebias)
        ebiasForce.addGlobalParameter("center", self.center)     # kJ/mol
        ebiasForce.setForceGroup(self.force_group)
        self.force = ebiasForce

    def defineInteraction(self):
        logger.debug("ElectrostaticBias: center=%s, k_ebias=%s, k_elec=%s, ldby=%s",
                     self.center, self.k_ebias, self.k_elec, self.ldby)


class BasePairProteinHarmonicBias(ProteinDNAForce):
    """Harmonic bias that holds a protein point near a target base-pair index.

    Two centroid groups span a stretch of DNA (``base_pair_left`` /
    ``base_pair_right``, ``base_pair_sep`` base pairs apart) and a third is a
    protein point. The protein's projection onto the left->right vector,
    expressed as a base-pair offset ``i``, is harmonically restrained with energy
    ``0.5 * k * i**2`` (``i = 0`` at the midpoint).

    The readout version is
    :class:`open3SPN2.force.collective_variables.BasePairProteinPosition`.
    """

    def __init__(self, dna, protein, base_pair_left_indicies, base_pair_right_indicies,
                 protein_indices, base_pair_sep=4, k=4.184,
                 force_group=force_groups.PROTEIN_DNA_BIAS):
        self.k = k
        self.base_pair_left_indicies = base_pair_left_indicies
        self.base_pair_right_indicies = base_pair_right_indicies
        self.protein_indices = protein_indices
        self.base_pair_sep = base_pair_sep
        self.force_group = force_group
        super().__init__(dna, protein)

    def reset(self):
        # k is in kJ/mol per base-pair**2; i is the deviation in base pairs.
        energy = f"0.5*{self.k}*(i)^2;"
        # ratio (0..1 along left->right) -> base-pair offset; ratio 0.5 -> i = 0.
        ratio_bp = f"i={self.base_pair_sep}*(ratio-0.5);"
        # ratio = (VA . VD) / (VD . VD): projection of the protein vector VA onto
        # the base-pair-separation vector VD, normalized by |VD|^2.
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


class AMHgoProteinDNA(ProteinDNAForce):
    """Structure-based (AMH-Go) contact potential between a protein and DNA. (Xinyu)

    Adds a Gaussian well at each native contact listed in ``contact_file``
    (columns: protein resSeq, DNA resSeq, native distance in angstrom, optional
    weight) between a protein CB (or CA) atom and a DNA base atom. The well depth
    at the native distance is ``k_3spn2 * k_amhgo_PD * gamma_ij``.
    """
    def __init__(self, dna, protein, chain_protein='A', chain_DNA='B', k_amhgo_PD=1*unit.kilocalorie_per_mole,
                 k_3spn2=1.0, sigma_sq=0.05*unit.nanometers**2, aaweight=False, globalct=True, cutoff=1.8,
                 contact_file="contact_protein_DNA.dat", force_group=force_groups.AMH_GO):
        self.force_group = force_group
        self.k_amhgo_PD = k_amhgo_PD
        self.k_3spn2 = k_3spn2
        self.sigma_sq = sigma_sq
        self.chain_protein = chain_protein
        self.chain_DNA = chain_DNA
        self.aaweight = aaweight
        self.cutoff = cutoff
        self.globalct = globalct
        self.contact_file = contact_file
        super().__init__(dna, protein)

    def reset(self):
        cutoff = self.cutoff
        k_3spn2 = self.k_3spn2
        if self.globalct:
            amhgoForce = openmm.CustomBondForce(f"-k_amhgo_PD*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_sq))*step({cutoff}-r)")
        else:
            amhgoForce = openmm.CustomBondForce(f"-k_amhgo_PD*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_sq))*step(r_ijN+{cutoff}-r)")
        amhgoForce.addGlobalParameter("k_amhgo_PD", k_3spn2 * self.k_amhgo_PD)
        amhgoForce.addGlobalParameter("sigma_sq", self.sigma_sq)
        amhgoForce.addPerBondParameter("gamma_ij")
        amhgoForce.addPerBondParameter("r_ijN")
        amhgoForce.setUsesPeriodicBoundaryConditions(self.periodic)
        amhgoForce.setForceGroup(self.force_group)  # one cutoff per force group
        self.force = amhgoForce

    def defineInteraction(self):
        atoms = self.dna.atoms.copy()
        atoms['index'] = atoms.index
        atoms.index = zip(atoms['chainID'], atoms['resSeq'], atoms['name'])

        contact_list = np.loadtxt(self.contact_file)
        for i in range(len(contact_list)):
            if self.aaweight:
                gamma_ij = contact_list[i][3]
            else:
                gamma_ij = 1.0
            if (self.chain_protein, int(contact_list[i][0]), 'CB') in atoms.index:
                CB_protein = atoms[(atoms['chainID'] == self.chain_protein) & (atoms['resSeq'] == int(contact_list[i][0])) & (atoms['name'] == 'CB') & atoms['resname'].isin(_proteinResidues)].copy()
            else:
                CB_protein = atoms[(atoms['chainID'] == self.chain_protein) & (atoms['resSeq'] == int(contact_list[i][0])) & (atoms['name'] == 'CA') & atoms['resname'].isin(_proteinResidues)].copy()
            base_DNA = atoms[(atoms['chainID'] == self.chain_DNA) & (atoms['resSeq'] == int(contact_list[i][1])) & (atoms['name'].isin(['A', 'T', 'G', 'C'])) & atoms['resname'].isin(_dnaResidues)].copy()
            r_ijN = contact_list[i][2] / 10.0 * unit.nanometers
            self.force.addBond(int(CB_protein['index'].values[0]), int(base_DNA['index'].values[0]), [gamma_ij, r_ijN])
            logger.debug('AMHgo contact: %s %s %s', int(CB_protein['index'].values[0]),
                         int(base_DNA['index'].values[0]), [gamma_ij, r_ijN])
