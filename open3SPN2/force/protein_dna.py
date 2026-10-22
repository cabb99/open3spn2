import numpy as np

import openmm
import openmm.unit as unit
import pandas
from .template import ProteinDNAForce
_af = 1 * unit.degree / unit.radian  # angle scaling factor
_dnaResidues = ['DA', 'DC', 'DT', 'DG']
_proteinResidues = ['IPR', 'IGL', 'NGP']

class ExclusionProteinDNA(ProteinDNAForce):
    """ Protein-DNA exclusion potential"""
    def __init__(self, dna, protein, k=1, radius_override = None, cutoff = 1.55, force_group=14):
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
        print(f"protein dna cutoff {exclusionForce.getCutoffDistance()}")
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
                print(f'Residue {i} not included in protein-DNA interactions')
                parameters = [0, .1, .1]
            atoms.loc[i, ['epsilon', 'radius', 'cutoff']] = parameters
            self.atoms = atoms
            self.force.addParticle(parameters)
        self.force.addInteractionGroup(DNA_list, protein_list)


class ElectrostaticsProteinDNA(ProteinDNAForce):
    """DNA-protein and protein-protein electrostatics."""
    def __init__(self, dna, protein, k=1, ldby = 1.2 * unit.nanometer, cutoff_distance = None, force_group=15):
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
            cutoff_distance = 4 * unit.nanometer
        else:
            cutoff_distance = self.cutoff_distance

        cutoff_nm = cutoff_distance.value_in_unit(unit.nanometer)

        electrostaticForce.setCutoffDistance(cutoff_nm)
        print(f"protein dna screening length {ldby} nm")
        print(f"protein dna cutoff {electrostaticForce.getCutoffDistance()} nm")
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


class BiasElectrostaticsProteinDNA(ProteinDNAForce):
    """ Protein-DNA string potential"""
    #k_ebias and center should be inputted
    def __init__(self, dna, protein, k_ebias,center, k_elec, ldby, cutoff_distance = None, forceGroup=16):
        self.k_ebias = k_ebias
        self.center = center
        self.k_elec = k_elec
        self.ldby = ldby
        self.forceGroup = forceGroup
        self.cutoff_distance = cutoff_distance
        super().__init__(dna, protein)

    def reset(self):
        #k_ebias=self.k_ebias.value_in_unit(unit.kilojoule_per_mole)
        #center=self.center.value_in_unit(unit.kilojoule_per_mole)
        k_ebias = self.k_ebias
        center = self.center
        ebiasForce = openmm.CustomCVForce(f"0.5*k_ebias*((E_elec-center)/4.184)^2")
        E_elec = ElectrostaticsProteinDNA(self.dna, self.protein, k = self.k_elec, ldby = self.ldby, cutoff_distance = self.cutoff_distance)
        elec = E_elec.force
        #ebiasForce.addCollectiveVariable("E_elec", E_elec)
        ebiasForce.addCollectiveVariable("E_elec", elec)    #Is in kJ/mol
        ebiasForce.addGlobalParameter("k_ebias", k_ebias)   #
        ebiasForce.addGlobalParameter("center", center)
        ebiasForce.setForceGroup(self.forceGroup)
        print(f"k_ebias = {k_ebias}")
        print(f"center = {center}")
        #print (E_elec)
        self.force = ebiasForce

    def defineInteraction(self):
        print(f"ElectrostaticsProteinDNA bias on: center at {self.center}, k_ebias = {self.k_ebias}, with electrostatic parameters k_elec = {self.k_elec} and screening length {self.ldby}")

class BasePairProteinHarmonicBias(ProteinDNAForce):
    def __init__(self, dna, protein, base_pair_left_indicies, base_pair_right_indicies, protein_indices, base_pair_sep = 4, k=4.184, forceGroup=16):
        #what indicies are members of protein_indicies, base_pair_left_indicies, and base_pair_right_indicies
        self.k = k
        self.base_pair_left_indicies = base_pair_left_indicies
        self.base_pair_right_indicies = base_pair_right_indicies
        self.protein_indices = protein_indices
        self.base_pair_sep = base_pair_sep
        self.forceGroup = forceGroup
        super().__init__(dna,protein)
    def reset(self):
        #Group 1 centered around base pairs i-base_pair_sep/2
        #Group 2 centered around base pairs i+base_pair_sep/2
        #Group 3 centered around protein point

        # Harmonic parameters
        k = self.k      # kJ/mol (adjust as needed)
        #k = self.k * 0.34 ** 2      # kJ/mol (adjust as needed)  (if need to adjust to kJ/(bp^2*mol))

        #Define the harmonic bias
        energy = f"0.5*{k}*(i)^2;"    # i is deviation base pair from target; be careful when processing for WHAM

        #Converting ratio to base pair
        ratio_bp = f"i={self.base_pair_sep}*(ratio-0.5);"    #base_pair_sep coefficient corresponds to base pairs i-base_pair_sep/2; i+base_pair_sep/2 where i is target base pair; midpoint type approximation 

        #Calculate ratio of dot products
        ratio_dots = "ratio=VAVDdots/VDVDdots;VAVDdots=VAx*VDx+VAy*VDy+VAz*VDz;VDVDdots=VDx*VDx+VDy*VDy+VDz*VDz;"

        '''
        # ChatGPT suggested simplification with cancellation out and removal of sqrt function; left as commented out as placeholder code

        #Trigonometry
        trig = "ratio=lengthVA*theta/lengthVD;"    #theta is a dot product ratio; arccos(theta) is the angle between VA and VD (not angle)

        #Angle Definitions
        angles = "theta=dot/(lengthVA*lengthVD);dot=VAx*VDx+VAy*VDy+VAz*VDz;lengthVA=sqrt(VAx^2+VAy^2+VAz^2);lengthVD=sqrt(VDx^2+VDy^2+VDz^2);"
        
        #If restoring this block of code; update expression initialization
        '''

        #Vector definitions
        vectors = "VAx=Px-BPLx;VAy=Py-BPLy;VAz=Pz-BPLz;VDx=BPRx-BPLx;VDy=BPRy-BPLy;VDz=BPRz-BPLz;"
        
        #Particle definitions
        particles= "BPLx=x1;BPLy=y1;BPLz=z1;BPRx=x2;BPRy=y2;BPRz=z2;Px=x3;Py=y3;Pz=z3"

        expression = f"{energy}{ratio_bp}{ratio_dots}{vectors}{particles}"

        #print(expression)

        force = openmm.CustomCentroidBondForce(3, expression)
        force.addGroup(self.base_pair_left_indicies)
        force.addGroup(self.base_pair_right_indicies)
        force.addGroup(self.protein_indices)
        force.addBond([0,1,2])
        force.setForceGroup(self.forceGroup)

        self.force = force
    def defineInteraction(self):
        pass

class proteinBasePairPosition(ProteinDNAForce):
    def __init__(self, dna, protein, base_pair_left_indicies, base_pair_right_indicies, protein_indices, base_pair_sep = 4, k=4.184, forceGroup=16):
        #what indicies are members of protein_indicies, base_pair_left_indicies, and base_pair_right_indicies
        self.k = k
        self.base_pair_left_indicies = base_pair_left_indicies
        self.base_pair_right_indicies = base_pair_right_indicies
        self.protein_indices = protein_indices
        self.base_pair_sep = base_pair_sep
        self.forceGroup = forceGroup
        super().__init__(dna,protein)
    def reset(self):
        #Group 1 centered around base pairs i-base_pair_sep/2
        #Group 2 centered around base pairs i+base_pair_sep/2
        #Group 3 centered around protein point

        # Harmonic parameters
        k = self.k      # kJ/mol (adjust as needed)
        #k = self.k * 0.34 ** 2      # kJ/mol (adjust as needed)  (if need to adjust to kJ/(bp^2*mol))

        #Define the harmonic bias
        energy = f"i;"    # i is deviation base pair from target; be careful when processing for WHAM

        #Converting ratio to base pair
        ratio_bp = f"i={self.base_pair_sep}*(ratio-0.5);"    #base_pair_sep coefficient corresponds to base pairs i-base_pair_sep/2; i+base_pair_sep/2 where i is target base pair; midpoint type approximation 

        #Calculate ratio of dot products
        ratio_dots = "ratio=VAVDdots/VDVDdots;VAVDdots=VAx*VDx+VAy*VDy+VAz*VDz;VDVDdots=VDx*VDx+VDy*VDy+VDz*VDz;"

        '''
        # ChatGPT suggested simplification with cancellation out and removal of sqrt function; left as commented out as placeholder code

        #Trigonometry
        trig = "ratio=lengthVA*theta/lengthVD;"    #theta is a dot product ratio; arccos(theta) is the angle between VA and VD (not angle)

        #Angle Definitions
        angles = "theta=dot/(lengthVA*lengthVD);dot=VAx*VDx+VAy*VDy+VAz*VDz;lengthVA=sqrt(VAx^2+VAy^2+VAz^2);lengthVD=sqrt(VDx^2+VDy^2+VDz^2);"
        
        #If restoring this block of code; update expression initialization
        '''

        #Vector definitions
        vectors = "VAx=Px-BPLx;VAy=Py-BPLy;VAz=Pz-BPLz;VDx=BPRx-BPLx;VDy=BPRy-BPLy;VDz=BPRz-BPLz;"
        
        #Particle definitions
        particles= "BPLx=x1;BPLy=y1;BPLz=z1;BPRx=x2;BPRy=y2;BPRz=z2;Px=x3;Py=y3;Pz=z3"

        expression = f"{energy}{ratio_bp}{ratio_dots}{vectors}{particles}"

        #print(expression)

        force = openmm.CustomCentroidBondForce(3, expression)
        force.addGroup(self.base_pair_left_indicies)
        force.addGroup(self.base_pair_right_indicies)
        force.addGroup(self.protein_indices)
        force.addBond([0,1,2])
        force.setForceGroup(self.forceGroup)

        self.force = force
    def defineInteraction(self):
        pass

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
#            base_DNA = atoms[(atoms['chainID'] == self.chain_DNA) & (atoms['resSeq'] == int(contact_list[i][1])) & (atoms['name'].isin(['A', 'T', 'G', 'C'])) & atoms['resname'].isin(_dnaResidues].copy()
#            r_ijN = contact_list[i][2]/10.0*unit.nanometers
#            self.force.addBond(int(CB_protein['index'].values[0]), int(base_DNA['index'].values[0]), [gamma_ij, r_ijN])
#            print(int(CB_protein['index'].values[0]), int(base_DNA['index'].values[0]), [gamma_ij, r_ijN])

class StringProteinDNA(ProteinDNAForce):
    """ Protein-DNA string potential (Xinyu)"""
    def __init__(self, dna, protein, r0, chain_protein='A', chain_DNA='B', k_string_PD=10*4.184, protein_seg=False, group=None):
        self.k_string_PD = k_string_PD
        self.chain_protein = chain_protein
        self.chain_DNA = chain_DNA
        self.r0 = r0
        self.protein_seg = protein_seg
        self.group = list(group) if group is not None else []
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

class MultiChainProteinDNA(ProteinDNAForce):
    """ Protein-DNA string potential for multiple chains (added and amended by Steven and ChatGPT from Xinyu)"""
    def __init__(self, dna, protein, r0, chain_protein='AB', chain_DNA='CD', k_string_PD=10*4.184, protein_seg=False, group=[], force_group=19): #forceGroup placeholder is 19.
        self.k_string_PD = k_string_PD
        self.chain_protein = chain_protein
        self.chain_DNA = chain_DNA
        self.r0 = r0
        self.protein_seg = protein_seg
        self.group = group
        self.force_group = force_group
        super().__init__(dna, protein)

    def reset(self):
        r0=self.r0
        k_string_PD=self.k_string_PD
        stringForce = openmm.CustomCentroidBondForce(2, f"0.5*{k_string_PD}*(distance(g1,g2)-{r0})^2")
        stringForce.setForceGroup(self.force_group)
        self.force = stringForce
        print("String_PD bias on: r0, k_string = ", r0, k_string_PD)

    def defineInteraction(self):
        atoms = self.dna.atoms.copy()
        atoms['index'] = atoms.index

        CA_atoms = atoms[
            atoms['chainID'].isin(list(self.chain_protein)) &
            (atoms['name'] == 'CA') &
            atoms['resname'].isin(_proteinResidues)
        ].copy()

        S_atoms = atoms[
            atoms['chainID'].isin(list(self.chain_DNA)) &
            (atoms['name'] == 'S') &
            atoms['resname'].isin(_dnaResidues)
        ].copy()

        if CA_atoms.empty or S_atoms.empty:
            raise ValueError("No CA or S atoms found in the specified chains.")

        CA_index = [int(atom.index) for atom in CA_atoms.itertuples()]
        if self.protein_seg:
            self.force.addGroup([CA_index[x] for x in self.group])
        else:
            self.force.addGroup(CA_index)

        self.force.addGroup([int(atom.index) for atom in S_atoms.itertuples()])
        self.force.addBond([0, 1])

        print("Protein group:", self.force.getGroupParameters(0))
        print("DNA group:", self.force.getGroupParameters(1))

class String_length_ProteinDNA(ProteinDNAForce):
    """ Protein-DNA string potential (Xinyu)"""
    def __init__(self, dna, protein, chain_protein='A', chain_DNA='B', protein_seg=False, group=None, force_group=17):
        self.force_group = force_group
        self.chain_protein = chain_protein
        self.chain_DNA = chain_DNA
        self.protein_seg = protein_seg
        self.group = list(group) if group is not None else []
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

class Multi_length_ProteinDNA(ProteinDNAForce):
    """ Protein-DNA string length for multiple chains (added and amended by Steven and ChatGPT from Xinyu)"""
    def __init__(self, dna, protein, chain_protein='AB', chain_DNA='CD', protein_seg=False, group=[], force_group=4): #force group placeholder is 4
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

        CA_atoms = atoms[
            atoms['chainID'].isin(list(self.chain_protein)) &
            (atoms['name'] == 'CA') &
            atoms['resname'].isin(_proteinResidues)
        ].copy()

        S_atoms = atoms[
            atoms['chainID'].isin(list(self.chain_DNA)) &
            (atoms['name'] == 'S') &
            atoms['resname'].isin(_dnaResidues)
        ].copy()

        if CA_atoms.empty or S_atoms.empty:
            raise ValueError("No CA or S atoms found in the specified chains.")

        CA_index = [int(atom.index) for atom in CA_atoms.itertuples()]
        if self.protein_seg:
            self.force.addGroup([CA_index[x] for x in self.group])
        else:
            self.force.addGroup(CA_index)

        self.force.addGroup([int(atom.index) for atom in S_atoms.itertuples()])
        self.force.addBond([0, 1])

        print("Protein group:", self.force.getGroupParameters(0))
        print("DNA group:", self.force.getGroupParameters(1))
