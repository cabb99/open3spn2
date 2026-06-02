import logging

import numpy as np

import openmm
import openmm.unit as unit
import pandas
from .template import ProteinDNAForce
from .dna import addNonBondedExclusions
from . import force_groups

logger = logging.getLogger(__name__)

_af = 1 * unit.degree / unit.radian  # angle scaling factor
_dnaResidues = ['DA', 'DC', 'DT', 'DG']
_proteinResidues = ['IPR', 'IGL', 'NGP']


def select_string_groups(dna, chain_protein='A', chain_DNA='B', protein_seg=False, group=None):
    """Pick the two centroid groups for a protein-DNA string/pulling potential.

    Returns ``(protein_CA_indices, dna_S_indices)``: the protein alpha-carbon
    (``CA``) atoms and the DNA sugar (``S``) atoms that define the ends of the
    pulling vector. ``chain_protein`` / ``chain_DNA`` accept a single chain id
    (``'A'``) or several (``'AB'`` or ``['A', 'B']``); every matching chain is
    pooled into one group. When ``protein_seg`` is true only the ``CA`` atoms at
    the positions listed in ``group`` (indices into the selected ``CA`` list) are
    used.
    """
    atoms = dna.atoms.copy()
    atoms['index'] = atoms.index
    CA_atoms = atoms[
        atoms['chainID'].isin(list(chain_protein)) &
        (atoms['name'] == 'CA') &
        atoms['resname'].isin(_proteinResidues)
    ]
    S_atoms = atoms[
        atoms['chainID'].isin(list(chain_DNA)) &
        (atoms['name'] == 'S') &
        atoms['resname'].isin(_dnaResidues)
    ]
    CA_index = [int(i) for i in CA_atoms['index']]
    S_index = [int(i) for i in S_atoms['index']]
    if protein_seg:
        CA_index = [CA_index[x] for x in (group or [])]
    return CA_index, S_index


class ExclusionProteinDNA(ProteinDNAForce):
    """ Protein-DNA exclusion potential"""
    def __init__(self, dna, protein, k=1, radius_override = None, cutoff = 1.55, force_group=force_groups.EXCLUSION_PROTEIN_DNA):
        self.k = k
        self.force_group = force_group
        self.radius_override = radius_override
        self.cutoff = cutoff    #cutoff is in nm
        super().__init__(dna, protein)

    def reset(self):
        exclusionForce = openmm.CustomNonbondedForce(f"""k_exclusion_protein_DNA*energy;
                         energy=(4*epsilon*((sigma/r)^12-(sigma/r)^6)-offset)*step(cutoff-r);
                         offset=4*epsilon*((sigma/cutoff)^12-(sigma/cutoff)^6);
                         sigma=0.5*(sigma1+sigma2); 
                         epsilon=sqrt(epsilon1*epsilon2);
                         cutoff=sqrt(cutoff1*cutoff2)""")
        exclusionForce.addGlobalParameter('k_exclusion_protein_DNA', self.k)
        exclusionForce.addPerParticleParameter('epsilon')
        exclusionForce.addPerParticleParameter('sigma')
        exclusionForce.addPerParticleParameter('cutoff')
        exclusionForce.setCutoffDistance(self.cutoff)
        # exclusionForce.setUseLongRangeCorrection(True)
        exclusionForce.setForceGroup(self.force_group)  # There can not be multiple cutoff distance on the same force group
        if self.periodic:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffPeriodic)
        else:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffNonPeriodic)
        logger.debug('protein-DNA exclusion cutoff %s', exclusionForce.getCutoffDistance())
        self.force = exclusionForce

    def defineInteraction(self):

        particle_definition = self.dna.config['Protein-DNA particles']
        dna_particle_definition=particle_definition[(particle_definition['molecule'] == 'DNA') &
                                                    (particle_definition['DNA'] == self.dna.DNAtype)]
        protein_particle_definition = particle_definition[(particle_definition['molecule'] == 'Protein')]

        # Merge DNA and protein particle definitions
        particle_definition = pandas.concat([dna_particle_definition, protein_particle_definition], sort=False)
        particle_definition.index = particle_definition.molecule + particle_definition.name
        self.particle_definition = particle_definition

        is_dna = self.dna.atoms['resname'].isin(_dnaResidues)
        is_protein = self.dna.atoms['resname'].isin(_proteinResidues)
        atoms = self.dna.atoms.copy()
        atoms['is_dna'] = is_dna
        atoms['is_protein'] = is_protein
        atoms['epsilon']=np.nan
        atoms['radius']=np.nan
        atoms['cutoff'] = np.nan
        DNA_list = []
        protein_list = []
        for i, atom in atoms.iterrows():
            if atom.is_dna:
                param = particle_definition.loc['DNA' + atom['name']]
                if self.radius_override == None:
                    parameters = [param.epsilon,
                                param.radius,
                                param.cutoff]
                else:
                    parameters = [param.epsilon,
                                self.radius_override,
                                param.cutoff]
                DNA_list += [i]
            elif atom.is_protein:
                param = particle_definition.loc['Protein' + atom['name']]
                if self.radius_override == None:
                    parameters = [param.epsilon,
                                param.radius,
                                param.cutoff]
                else:
                    parameters = [param.epsilon,
                                self.radius_override,
                                param.cutoff]
                protein_list += [i]
            else:
                logger.warning('Residue %s not included in protein-DNA interactions', i)
                parameters = [0, .1,.1]
            atoms.loc[i, ['epsilon', 'radius', 'cutoff']] = parameters
            self.atoms = atoms
            self.force.addParticle(parameters)
        self.force.addInteractionGroup(DNA_list, protein_list)

        # addExclusions
        addNonBondedExclusions(self.dna, self.force)


class ElectrostaticsProteinDNA(ProteinDNAForce):
    """DNA-protein and protein-protein electrostatics."""
    def __init__(self, dna, protein, k=1, ldby = 1.2 * unit.nanometer, cutoff_distance = None, force_group=force_groups.ELECTROSTATICS_PROTEIN_DNA):
        self.k = k
        self.force_group = force_group
        self.ldby = ldby
        self.cutoff_distance = cutoff_distance
        super().__init__(dna, protein)

    def reset(self):
        dielectric = 78 # e * a
        #print(dielectric)
        # Debye length
        Na = unit.AVOGADRO_CONSTANT_NA  # Avogadro number
        ec = 1.60217653E-19 * unit.coulomb  # proton charge
        pv = 8.8541878176E-12 * unit.farad / unit.meter  # dielectric permittivity of vacuum

        #ldby = 1.2 * unit.nanometer # np.sqrt(dielectric * pv * kb * T / (2.0 * Na * ec ** 2 * C))
        denominator = 4 * np.pi * pv * dielectric / (Na * ec ** 2)
        denominator = denominator.in_units_of(unit.kilocalorie_per_mole**-1 * unit.nanometer**-1)
        #print(ldby, denominator)
        k = self.k
        ldby = self.ldby
        electrostaticForce = openmm.CustomNonbondedForce(f"""k_electro_protein_DNA*energy;
                             energy=q1*q2*exp(-r/inter_dh_length)/inter_denominator/r;""")
        electrostaticForce.addPerParticleParameter('q')
        electrostaticForce.addGlobalParameter('k_electro_protein_DNA', k)
        electrostaticForce.addGlobalParameter('inter_dh_length', ldby)
        electrostaticForce.addGlobalParameter('inter_denominator', denominator)

        if self.cutoff_distance == None:
            cutoff_distance = 4 * unit.nanometer # for backward compatibility
        else:
            cutoff_distance = self.cutoff_distance

        cutoff_nm = cutoff_distance.value_in_unit(unit.nanometer)

        electrostaticForce.setCutoffDistance(cutoff_nm)
        logger.debug('protein-DNA screening length %s nm', ldby)
        logger.debug('protein-DNA electrostatic cutoff %s nm', electrostaticForce.getCutoffDistance())
        if self.periodic:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffPeriodic)
        else:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffNonPeriodic)
        electrostaticForce.setForceGroup(self.force_group)
        self.force = electrostaticForce

    def defineInteraction(self):
        # Merge DNA and protein particle definitions
        particle_definition = self.dna.config['Protein-DNA particles']
        dna_particle_definition=particle_definition[(particle_definition['molecule'] == 'DNA') &
                                                    (particle_definition['DNA'] == self.dna.DNAtype)]
        protein_particle_definition = particle_definition[(particle_definition['molecule'] == 'Protein')]

        # Merge DNA and protein particle definitions
        particle_definition = pandas.concat([dna_particle_definition, protein_particle_definition], sort=False)
        particle_definition.index = particle_definition.molecule + particle_definition.name
        self.particle_definition = particle_definition

        # Open Sequence dependent electrostatics
        sequence_electrostatics = self.dna.config['Sequence dependent electrostatics']
        sequence_electrostatics.index = sequence_electrostatics.resname

        # Select only dna and protein atoms
        is_dna = self.protein.atoms['resname'].isin(_dnaResidues)
        is_protein = self.protein.atoms['resname'].isin(_proteinResidues)
        atoms = self.protein.atoms.copy()
        atoms['is_dna'] = is_dna
        atoms['is_protein'] = is_protein
        DNA_list = []
        protein_list = []

        for i, atom in atoms.iterrows():
            if atom.is_dna:
                param = particle_definition.loc['DNA' + atom['name']]
                charge = param.charge
                parameters = [charge]
                if charge != 0:
                    DNA_list += [i]
                    #print(atom.chainID, atom.resSeq, atom.resname, atom['name'], charge)
            elif atom.is_protein:
                atom_param = particle_definition.loc['Protein' + atom['name']]
                seq_param = sequence_electrostatics.loc[atom.real_resname]
                charge = atom_param.charge * seq_param.charge
                parameters = [charge]
                if charge != 0:
                    protein_list += [i]
                    #print(atom.chainID, atom.resSeq, atom.resname, atom['name'], charge)
            else:
                logger.warning('Residue %s not included in protein-DNA electrostatics', i)
                parameters = [0]  # No charge if it is not DNA
            # print (i,parameters)
            self.force.addParticle(parameters)
        self.force.addInteractionGroup(DNA_list, protein_list)
        # self.force.addInteractionGroup(protein_list, protein_list) #protein-protein electrostatics should be included using debye Huckel Terms

        # addExclusions
        addNonBondedExclusions(self.dna, self.force)


class AMHgoProteinDNA(ProteinDNAForce):
    """Structure-based (AMH-Go) contact potential between a protein and DNA. (Xinyu)

    Adds a Gaussian well at each native contact listed in ``contact_file``
    (columns: protein resSeq, DNA resSeq, native distance in angstrom, optional
    weight) between a protein CB (or CA) atom and a DNA base atom.
    """
    def __init__(self, dna, protein, chain_protein='A', chain_DNA='B', k_amhgo_PD=1*unit.kilocalorie_per_mole,
                 k_3spn2=1.0, sigma_sq=0.05*unit.nanometers**2, aaweight=False, globalct=True, cutoff=1.8,
                 contact_file="contact_protein_DNA.dat", force_group=force_groups.AMH_GO):
        self.force_group = force_group
        self.k_amhgo_PD = k_amhgo_PD
        self.k_3spn2 = k_3spn2
        self.sigma_sq= sigma_sq
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
        amhgoForce.addGlobalParameter("k_amhgo_PD", k_3spn2*self.k_amhgo_PD)
        amhgoForce.addGlobalParameter("sigma_sq", self.sigma_sq)
        amhgoForce.addPerBondParameter("gamma_ij")
        amhgoForce.addPerBondParameter("r_ijN")
        amhgoForce.setUsesPeriodicBoundaryConditions(self.periodic)
        amhgoForce.setForceGroup(self.force_group)  # There can not be multiple cutoff distance on the same force group
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
            r_ijN = contact_list[i][2]/10.0*unit.nanometers
            self.force.addBond(int(CB_protein['index'].values[0]), int(base_DNA['index'].values[0]), [gamma_ij, r_ijN])
            logger.debug('AMHgo contact: %s %s %s', int(CB_protein['index'].values[0]),
                         int(base_DNA['index'].values[0]), [gamma_ij, r_ijN])

