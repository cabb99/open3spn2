import numpy as np

import openmm
import openmm.unit as unit
from .template import ProteinDNAForce
_af = 1 * unit.degree / unit.radian  # angle scaling factor
_dnaResidues = ['DA', 'DC', 'DT', 'DG']
_proteinResidues = ['IPR', 'IGL', 'NGP']

class ExclusionProteinDNA(ProteinDNAForce):
    """ Protein-DNA exclusion potential"""
    def __init__(self, dna, protein, k=1, force_group=14):
        self.k = k
        self.force_group = force_group
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
        exclusionForce.setCutoffDistance(1.55)
        # exclusionForce.setUseLongRangeCorrection(True)
        exclusionForce.setForceGroup(self.force_group)  # There can not be multiple cutoff distance on the same force group
        if self.periodic:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffPeriodic)
        else:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffNonPeriodic)
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
                parameters = [param.epsilon,
                              param.radius,
                              param.cutoff]
                DNA_list += [i]
            elif atom.is_protein:
                param = particle_definition.loc['Protein' + atom['name']]
                parameters = [param.epsilon,
                              param.radius,
                              param.cutoff]
                protein_list += [i]
            else:
                print(f'Residue {i} not included in protein-DNA interactions')
                parameters = [0, .1,.1]
            atoms.loc[i, ['epsilon', 'radius', 'cutoff']] = parameters
            self.atoms = atoms
            self.force.addParticle(parameters)
        self.force.addInteractionGroup(DNA_list, protein_list)

        # addExclusions
        addNonBondedExclusions(self.dna, self.force)


class ElectrostaticsProteinDNA(ProteinDNAForce):
    """DNA-protein and protein-protein electrostatics."""
    def __init__(self, dna, protein, k=1, force_group=15):
        self.k = k
        self.force_group = force_group
        super().__init__(dna, protein)

    def reset(self):
        dielectric = 78 # e * a
        #print(dielectric)
        # Debye length
        Na = unit.AVOGADRO_CONSTANT_NA  # Avogadro number
        ec = 1.60217653E-19 * unit.coulomb  # proton charge
        pv = 8.8541878176E-12 * unit.farad / unit.meter  # dielectric permittivity of vacuum

        ldby = 1.2 * unit.nanometer # np.sqrt(dielectric * pv * kb * T / (2.0 * Na * ec ** 2 * C))
        denominator = 4 * np.pi * pv * dielectric / (Na * ec ** 2)
        denominator = denominator.in_units_of(unit.kilocalorie_per_mole**-1 * unit.nanometer**-1)
        #print(ldby, denominator)
        k = self.k
        electrostaticForce = openmm.CustomNonbondedForce(f"""k_electro_protein_DNA*energy;
                             energy=q1*q2*exp(-r/inter_dh_length)/inter_denominator/r;""")
        electrostaticForce.addPerParticleParameter('q')
        electrostaticForce.addGlobalParameter('k_electro_protein_DNA', k)
        electrostaticForce.addGlobalParameter('inter_dh_length', ldby)
        electrostaticForce.addGlobalParameter('inter_denominator', denominator)

        electrostaticForce.setCutoffDistance(4)
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
                print(f'Residue {i} not included in protein-DNA electrostatics')
                parameters = [0]  # No charge if it is not DNA
            # print (i,parameters)
            self.force.addParticle(parameters)
        self.force.addInteractionGroup(DNA_list, protein_list)
        # self.force.addInteractionGroup(protein_list, protein_list) #protein-protein electrostatics should be included using debye Huckel Terms

        # addExclusions
        addNonBondedExclusions(self.dna, self.force)


class AMHgoProteinDNA(ProteinDNAForce):
    """ Protein-DNA amhgo potential"""
    def __init__(self, dna, protein, chain_protein='A', chain_DNA='B', k_amhgo_PD=1*unit.kilocalorie_per_mole, sigma_sq=0.05*unit.nanometers**2, aaweight=False, globalct=True, cutoff=1.8, force_group=16):
        self.force_group = force_group
        self.k_amhgo_PD = k_amhgo_PD
        self.sigma_sq= sigma_sq
        self.chain_protein = chain_protein
        self.chain_DNA = chain_DNA
        self.aaweight = aaweight
        self.cutoff = cutoff
        self.globalct = globalct
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

        contact_list = np.loadtxt("contact_protein_DNA.dat")
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
            print(int(CB_protein['index'].values[0]), int(base_DNA['index'].values[0]), [gamma_ij, r_ijN])


#class AMHgoProteinDNA(ProteinDNAForce):
#    """ Protein-DNA amhgo potential (Xinyu)"""
#    def __init__(self, dna, protein, chain_protein='A', chain_DNA='B', k_amhgo_PD=1*unit.kilocalorie_per_mole,
#                 sigma_sq=0.05*unit.nanometers**2, aaweight=False, cutoff=1.8, force_group=16):
#        self.force_group = force_group
#        self.k_amhgo_PD = k_amhgo_PD
#        self.sigma_sq= sigma_sq
#        self.chain_protein = chain_protein
#        self.chain_DNA = chain_DNA
#        self.aaweight = aaweight
#        self.cutoff = cutoff
#        super().__init__(dna, protein)
#
#    def reset(self):
#        cutoff = self.cutoff
#        amhgoForce = openmm.CustomBondForce(f"-k_amhgo_PD*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_sq))*step({cutoff}-r)")
#        amhgoForce.addGlobalParameter("k_amhgo_PD", self.k_amhgo_PD)
#        amhgoForce.addGlobalParameter("sigma_sq", self.sigma_sq)
#        amhgoForce.addPerBondParameter("gamma_ij")
#        amhgoForce.addPerBondParameter("r_ijN")
#        amhgoForce.setUsesPeriodicBoundaryConditions(self.periodic)
#        amhgoForce.setForceGroup(self.force_group)  # There can not be multiple cutoff distance on the same force group
#        self.force = amhgoForce
#
#    def defineInteraction(self):
#        atoms = self.dna.atoms.copy()
#        atoms['index'] = atoms.index
#        atoms.index = zip(atoms['chainID'], atoms['resSeq'], atoms['name'])
#
#        contact_list = np.loadtxt("contact_protein_DNA.dat")
#        for i in range(len(contact_list)):
#            if self.aaweight: gamma_ij = contact_list[i][3]
#            else:   gamma_ij = 1.0
#            if (self.chain_protein, int(contact_list[i][0]), 'CB') in atoms.index:
#                 CB_protein = atoms[(atoms['chainID'] == self.chain_protein) & (atoms['resSeq'] == int(contact_list[i][0])) &
#                                    (atoms['name'] == 'CB') & atoms['resname'].isin(_proteinResidues)].copy()
#            else:
#                 CB_protein = atoms[(atoms['chainID'] == self.chain_protein) & (atoms['resSeq'] == int(contact_list[i][0])) &
#                                    (atoms['name'] == 'CA') & atoms['resname'].isin(_proteinResidues)].copy()
#            base_DNA = atoms[(atoms['chainID'] == self.chain_DNA) & (atoms['resSeq'] == int(contact_list[i][1])) & (atoms['name'].isin(['A', 'T', 'G', 'C'])) & atoms['resname'].isin(_dnaResidues)].copy()
#            r_ijN = contact_list[i][2]/10.0*unit.nanometers
#            self.force.addBond(int(CB_protein['index'].values[0]), int(base_DNA['index'].values[0]), [gamma_ij, r_ijN])
#            print(int(CB_protein['index'].values[0]), int(base_DNA['index'].values[0]), [gamma_ij, r_ijN])

class StringProteinDNA(ProteinDNAForce):
    """ Protein-DNA string potential (Xinyu)"""
    def __init__(self, dna, protein, r0, chain_protein='A', chain_DNA='B', k_string_PD=10*4.184, protein_seg=False, group=[]):
        self.k_string_PD = k_string_PD
        self.chain_protein = chain_protein
        self.chain_DNA = chain_DNA
        self.r0 = r0
        self.protein_seg = protein_seg
        self.group = group
        super().__init__(dna, protein)

    def reset(self):
        r0=self.r0
        k_string_PD=self.k_string_PD
        stringForce = openmm.CustomCentroidBondForce(2, f"0.5*{k_string_PD}*(distance(g1,g2)-{r0})^2")
        self.force = stringForce
        print("String_PD bias on: r0, k_string = ", r0, k_string_PD)

    def defineInteraction(self):
        atoms = self.dna.atoms.copy()
        atoms['index'] = atoms.index
        CA_atoms = atoms[(atoms['chainID'] == self.chain_protein) & (atoms['name'] == 'CA') & atoms['resname'].isin(_proteinResidues)].copy()
        S_atoms = atoms[(atoms['chainID'] == self.chain_DNA) & (atoms['name'] == 'S') & atoms['resname'].isin(_dnaResidues)].copy()
        CA_index = [int(atom.index) for atom in CA_atoms.itertuples()]
        if self.protein_seg: self.force.addGroup([CA_index[x] for x in self.group])
        else:   self.force.addGroup(CA_index)
        self.force.addGroup([int(atom.index) for atom in S_atoms.itertuples()])
        bondGroups = [0, 1]
        print(self.force.getGroupParameters(0))
        print(self.force.getGroupParameters(1))

        self.force.addBond(bondGroups)


class String_length_ProteinDNA(ProteinDNAForce):
    """ Protein-DNA string potential (Xinyu)"""
    def __init__(self, dna, protein, chain_protein='A', chain_DNA='B', protein_seg=False, group=[], force_group=17):
        self.force_group = force_group
        self.chain_protein = chain_protein
        self.chain_DNA = chain_DNA
        self.protein_seg = protein_seg
        self.group = group
        super().__init__(dna, protein)

    def reset(self):
        length = openmm.CustomCentroidBondForce(2, "distance(g1,g2)")
        length.setForceGroup(self.force_group)
        self.force = length

    def defineInteraction(self):
        atoms = self.dna.atoms.copy()
        atoms['index'] = atoms.index
        CA_atoms = atoms[(atoms['chainID'] == self.chain_protein) & (atoms['name'] == 'CA') & atoms['resname'].isin(_proteinResidues)].copy()
        S_atoms = atoms[(atoms['chainID'] == self.chain_DNA) & (atoms['name'] == 'S') & atoms['resname'].isin(_dnaResidues)].copy()
        CA_index = [int(atom.index) for atom in CA_atoms.itertuples()]
        if self.protein_seg: self.force.addGroup([CA_index[x] for x in self.group])
        else:   self.force.addGroup(CA_index)
        self.force.addGroup([int(atom.index) for atom in S_atoms.itertuples()])
        bondGroups = [0, 1]
        print(self.force.getGroupParameters(0))
        print(self.force.getGroupParameters(1))

        self.force.addBond(bondGroups)