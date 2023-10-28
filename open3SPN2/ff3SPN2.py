#!/usr/bin/python3
"""
This module implements the 3SPN2 forcefield and 3SPN2.C forcefield in openmm.
It also contains Protein-DNA interaction potentials to be used with openAWSEM.
"""
# TODO Curved BDNA is currently using the original pdb template, it should use X3DNA if possible,
#  so I need to make it able to take parameters from another pdb if necessary


__author__ = 'Carlos Bueno'
__version__ = '0.3.2'

import configparser
import numpy as np

import openmm.app
import openmm
import openmm.unit as unit
import scipy.spatial.distance as sdist

import pandas
from pathlib import Path

from .geometry import Transform
from .force import *
from .parser import *

__location__ = Path(__file__).parent.resolve()
_ef = 1 * unit.kilocalorie / unit.kilojoule  # energy scaling factor
_df = 1 * unit.angstrom / unit.nanometer  # distance scaling factor
_af = 1 * unit.degree / unit.radian  # angle scaling factor
_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
_dnaResidues = ['DA', 'DC', 'DT', 'DG']
_proteinResidues = ['IPR', 'IGL', 'NGP']
xml = __location__/'3SPN2.xml'

class DNA(object):
    """ A Coarse Grained DNA object."""
    def __init__(self, periodic=True):
        """Initializes an DNA object"""
        self.periodic = periodic

    def __repr__(self):
        return f'3SPN2 DNA object ({len(self.atoms)} atoms)'
        # print the sequence and the identity of the DNA object

    def parseConfigurationFile(self, configuration_file=__location__ /'3SPN2.conf'):
        """Reads the configuration file for the forcefield. The default configuration file is 3SPN2.conf
        and it contains most of the parameters used in the simulation."""
        self.configuration_file = configuration_file
        config = configparser.ConfigParser()
        config.read(configuration_file)

        # Parse all sections of the configuration file
        self.config = {}
        for c in config.sections():
            self.config.update({c: parseConfigTable(config[c])})

        # Assign main sections to variables
        self.particle_definition = self.config['Particles']
        self.bond_definition = self.config['Bonds']
        self.angle_definition = self.config['Harmonic Angles']
        self.dihedral_definition = self.config['Dihedrals']
        self.stacking_definition = self.config['Base Stackings']
        self.pair_definition = self.config['Base Pairs']
        self.cross_definition = self.config['Cross Stackings']

    def getSequences(self):
        """ Returns the DNA sequence as a Pandas Series. The index of the Series is (Chain, resid)"""
        dna_data = self.atoms[self.atoms.resname.isin(_dnaResidues)].copy()
        sequences = {}
        for c, chain in dna_data.groupby('chainID'):
            chain = chain.copy()
            resix = chain.resSeq
            res_unique = resix.unique()
            # chain['resID'] = resix.replace(dict(zip(res_unique, range(len(res_unique)))))
            sequences.update({(c, i): r.iloc[0]['resname'][1] for i, r in chain.groupby('resSeq')})
        self.sequence = pandas.Series(sequences)
        return self.sequence

    def computeGeometry(self, sequence=None, temp_name='temp'):
        """ This function requires X3DNA. It returns a pdb table containing the expected DNA structure"""
        #print("Computing geometry")
        pair = self.config['Base Pair Geometry']
        step = self.config['Base Step Geometry']
        pair.index = pair['stea']
        step.index = step['stea'] + step['steb']
        data = []
        _s = None
        if sequence is None:
            sequence = self.getSequences()
        seq = ''.join(sequence.values)
        for s in seq:
            pair_s = pair.loc[s, ['shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening']]
            if _s:
                step_s = step.loc[_s + s, ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']]
            else:
                step_s = (step.loc['AA', ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']] + 100) * 0
            data += [pandas.concat([pandas.Series([f'{s}-{_complement[s]}'], index=['Sequence']), pair_s, step_s])]
            _s = s
        data = pandas.concat(data, axis=1).T
        # try:
        #     location_x3dna = os.environ["X3DNA"]
        # except KeyError as ex:
        #     raise X3DNAnotFound from ex

        # with open(f'{temp_name}_parameters.par', 'w+') as par:
        #     par.write(f' {len(data)} # Number of base pairs\n')
        #     par.write(f' 0 # local base-pair & step parameters\n')
        #     par.write('#')
        #     par.write(data.to_csv(sep=' ', index=False))

        # # Attempting to call rebuild multiple times
        # # This function fails sometimes when called by multiple because a file with the same name is created
        # attempt = 0
        # max_attempts = 10
        # while attempt < max_attempts:
        #     try:
        #         subprocess.check_output([f'{location_x3dna}/bin/x3dna_utils',
        #                                  'cp_std', 'BDNA'])
        #         subprocess.check_output([f'{location_x3dna}/bin/rebuild',
        #                                  '-atomic', f'{temp_name}_parameters.par',
        #                                  f'{temp_name}_template.pdb'])
        #         break
        #     except subprocess.CalledProcessError as e:
        #         attempt += 1
        #         if attempt == max_attempts:
        #             print(f"subprocess.CalledProcessError failed {max_attempts} times {e.args[0]}: {e.args[1]}")
        # template_dna = self.fromPDB(f'{temp_name}_template.pdb', output_pdb=f'{temp_name}_temp.pdb', compute_topology=False)
        # template = template_dna.atoms.copy()

        new_orien = np.eye(3)
        new_pos = np.zeros(3)

        forward_strand=[]
        reverse_strand=[]
        xyz = pandas.read_csv(__location__/'DNA_atomic.csv')

        for row in data.itertuples():
            base1 = row.Sequence[0]
            base2 = row.Sequence[-1]
            base_parameters = np.array(row[2:8])
            step_parameters = np.array(row[8:16])
            
            base_parameters[3:]*=np.pi/180
            step_parameters[3:]*=np.pi/180
            
            t=Transform(*base_parameters)

            b1 = xyz[xyz['base']==base1].copy()
            b2 = xyz[xyz['base']==base2].copy()
            
            xyz1=b1[['x','y','z']].copy()
            xyz1 = np.dot(b1[['x','y','z']],t.full_rotation.T)+t.full_displacement
            xyz1 = np.dot((xyz1-t.half_displacement),t.half_rotation)
            
            xyz2 = b2[['x','y','z']].copy()
            xyz2[['y','z']]=-xyz2[['y','z']]
            xyz2 = np.dot((xyz2-t.half_displacement),t.half_rotation)
            
            # Step parameters
            t=Transform(*step_parameters)
            new_pos+=np.dot(t.full_displacement,new_orien)
            new_orien=np.dot(t.full_rotation.T,new_orien)
            
            xyz1 = np.dot(xyz1,new_orien)+new_pos
            #xyz2 = np.dot(xyz2,new_orien)+new_pos
            
            b1[['x','y','z']] = xyz1
            b2[['x','y','z']] = xyz2
            
            for b in [b1,b2]:
                b['recname']='ATOM'
                b['resname']='D'+b['base']
                b['resSeq']=row[0]+1
                b['occupancy']=1.00
                b['tempFactor']=1.00
                b['element']=b['name'].str[0]
            
            b1['chainID']='A'
            b2['chainID']='B'

            forward_strand=forward_strand + [b1]
            reverse_strand=[b2] + reverse_strand
            

        template=pandas.concat(forward_strand+reverse_strand)
        template.reindex()
        template['serial']=template.index+1
        writePDB(template,pdb_file='temp_template.pdb')

        try:
            self.atoms
        except AttributeError:
            return template
        
        template = template[template['chainID'] == 'A']
        original = self.atoms.copy()
        original.index = original['chainID'].astype(str) + '_' + original['resSeq'].astype(str) + '_' + original['name']
        ii = []
        for i, r in template.iterrows():
            name = r['name']
            seq_template = r['resname'][1:]
            j = r['resSeq'] - 1
            seq_seq = self.sequence.iloc[j]
            assert seq_seq == seq_template
            chain, resseq = self.sequence.index[j]
            ii += [f'{chain}_{resseq}_{name}']
        template.index = ii
        merge = pandas.merge(original, template, left_index=True, right_index=True, how='left', suffixes=['_old', ''])
        original[['x', 'y', 'z']] = merge[['x', 'y', 'z']]
        original.index = self.atoms.index
        return original

    def computeTopology(self, template_from_X3DNA=True, temp_name='temp'):
        """ Creates tables of bonds, angles and dihedrals with their respective parameters (bonded interactions).
        3SPN2.C requires a template structure to calculate the equilibrium bonds, angles and dihedrals.
        If template_from_structure is True, it will try to compute the equilibrium geometry using X3DNA.
        If template_from_structure is False, then the initial structure is expected to be the equilibrium geometry"""
        # Parse configuration file if not already done
        try:

            self.bond_definition
        except AttributeError:
            self.parseConfigurationFile()

        DNAtype = self.DNAtype
        if DNAtype not in self.angle_definition['DNA'].unique():
            raise DNATypeError(self)

        # Rewrite index in case it is not ordered
        self.atoms.index = range(len(self.atoms))

        # Compute B_curved geometry if needed
        if DNAtype == 'B_curved' and template_from_X3DNA:
            self.template_atoms = self.computeGeometry(temp_name=temp_name)
        else:
            self.template_atoms = self.atoms
        # Make an index to build the topology
        index = {}
        cr_list = set()  # Chain residue list
        for i, atom in self.atoms.iterrows():
            index.update({(atom['chainID'], atom['resSeq'], atom['name']): i})
            cr_list.update([(atom['chainID'], atom['resSeq'])])
        cr_list = list(cr_list)
        cr_list.sort()
        # max_chain = self.atoms['chain'].max()
        # max_residue = self.atoms['resSeq'].max()
        assert len(index) == len(self.atoms), 'Atom index was repeated'

        # Select ADNA bond definitions

        bond_types = self.bond_definition[self.bond_definition['DNA'] == DNAtype]
        angle_types = self.angle_definition[self.angle_definition['DNA'] == DNAtype]
        stacking_types = self.stacking_definition[self.stacking_definition['DNA'] == DNAtype]
        dihedral_types = self.dihedral_definition[self.dihedral_definition['DNA'] == DNAtype]
        # print(bond_types)
        # print(index)
        # Make a table with bonds
        data = []
        for i, ftype in bond_types.iterrows():
            # print(bond_type)
            ai = ftype['i']
            aj = ftype['j']
            s1 = ftype['s1']
            for c, r in cr_list:
                k1 = (c, r, ai)
                k2 = (c, r + s1, aj)
                if k1 in index and k2 in index:
                    data += [[i, index[k1], index[k2]]]
        data = pandas.DataFrame(data, columns=['name', 'aai', 'aaj'])
        self.bonds = data.merge(bond_types, left_on='name', right_index=True)

        if DNAtype == 'B_curved':
            # Make default distances the same as the initial distance
            x1 = self.template_atoms.loc[self.bonds['aai']][['x', 'y', 'z']]
            x2 = self.template_atoms.loc[self.bonds['aaj']][['x', 'y', 'z']]
            self.bonds['r0'] = np.diag(sdist.cdist(x1, x2))/10

        # Make a table with angles
        data = []
        base = self.atoms['resname'].str[1:2]
        for i, ftype in angle_types.iterrows():
            # if ftype.name != 37:
            #    continue
            # print(bond_type)
            ai = ftype['i']
            aj = ftype['j']
            ak = ftype['k']
            s1 = ftype['s1']
            s2 = ftype['s2']
            b1 = ftype['Base1']
            b2 = ftype['Base2']
            sb = ftype['sB']
            for c, r in cr_list:
                # print(ftype)
                k1 = (c, r, ai)
                k2 = (c, r + s1, aj)
                k3 = (c, r + s2, ak)
                k4 = (c, r + sb, 'S')
                if (
                    k1 in index
                    and k2 in index
                    and k3 in index
                    and k4 in index
                    and (b1 == '*' or base[index[k1]] == b1)
                    and (b2 == '*' or base[index[k4]] == b2)
                ):
                    data += [[i, index[k1], index[k2], index[k3], index[k4], sb]]
        data = pandas.DataFrame(data, columns=['name', 'aai', 'aaj', 'aak', 'aax', 'sB'])
        self.angles = data.merge(angle_types, left_on='name', right_index=True)

        if DNAtype == 'B_curved':
            # Make initial angles the default angles
            v1 = self.template_atoms.loc[self.angles['aai']].reset_index()[['x', 'y', 'z']]
            v2 = self.template_atoms.loc[self.angles['aaj']].reset_index()[['x', 'y', 'z']]
            v3 = self.template_atoms.loc[self.angles['aak']].reset_index()[['x', 'y', 'z']]
            a = v1 - v2
            a = np.array(a) / np.linalg.norm(a, keepdims=True, axis=1, )
            b = v3 - v2
            b = np.array(b) / np.linalg.norm(b, keepdims=True, axis=1, )
            self.angles['t0'] = np.arccos(np.einsum('ij,ij->i', a, b)) / np.pi * 180

        # Make a table with stackings
        data = []
        for i, ftype in stacking_types.iterrows():
            # print(bond_type)
            ai = ftype['i']
            aj = ftype['j']
            ak = ftype['k']
            s1 = ftype['s1']
            s2 = ftype['s2']
            for c, r in cr_list:
                k1 = (c, r, ai)
                k2 = (c, r + s1, aj)
                k3 = (c, r + s2, ak)
                if k1 in index and k2 in index and k3 in index:
                    data += [[i, index[k1], index[k2], index[k3]]]
        data = pandas.DataFrame(data, columns=['name', 'aai', 'aaj', 'aak'])
        self.stackings = data.merge(stacking_types, left_on='name', right_index=True)

        # Make a table with dihedrals
        data = []
        for i, ftype in dihedral_types.iterrows():
            # print(bond_type)
            ai = ftype['i']
            aj = ftype['j']
            ak = ftype['k']
            al = ftype['l']
            s1 = ftype['s1']
            s2 = ftype['s2']
            s3 = ftype['s3']
            for c, r in cr_list:
                k1 = (c, r, ai)
                k2 = (c, r + s1, aj)
                k3 = (c, r + s2, ak)
                k4 = (c, r + s3, al)
                if k1 in index and k2 in index and k3 in index and k4 in index:
                    data += [[i, index[k1], index[k2], index[k3], index[k4]]]
        data = pandas.DataFrame(data, columns=['name', 'aai', 'aaj', 'aak', 'aal'])
        self.dihedrals = data.merge(dihedral_types, left_on='name', right_index=True)

        if DNAtype == 'B_curved':
            # Make initial dihedrals the default dihedrals

            a1 = self.template_atoms.loc[self.dihedrals['aai']].reset_index()[['x', 'y', 'z']]
            a2 = self.template_atoms.loc[self.dihedrals['aaj']].reset_index()[['x', 'y', 'z']]
            a3 = self.template_atoms.loc[self.dihedrals['aak']].reset_index()[['x', 'y', 'z']]
            a4 = self.template_atoms.loc[self.dihedrals['aal']].reset_index()[['x', 'y', 'z']]

            b1 = np.array(a2 - a1)
            b2 = np.array(a3 - a2)
            b3 = np.array(a4 - a3)

            n1 = np.cross(b1, b2)
            n2 = np.cross(b2, b3)

            n1 /= np.linalg.norm(n1, axis=1, keepdims=True)
            n2 /= np.linalg.norm(n2, axis=1, keepdims=True)
            b2 /= np.linalg.norm(b2, axis=1, keepdims=True)

            m1 = np.cross(n1, b2)
            x = np.einsum('ij,ij->i', n1, n2)
            y = np.einsum('ij,ij->i', m1, n2)

            d = np.arctan2(y, x) / np.pi * 180
            self.dihedrals['t0'] = -d - 180

    def writePDB(self, pdb_file='clean.pdb'):
        # Compute element fields
        element_ix = {'P': 'P', 'S': 'H', 'A': 'N', 'T': 'S', 'C': 'O', 'G': 'C'}  # Elements choosen to keep VMD colors
        self.atoms.loc[:, 'element'] = [element_ix[atomType] for atomType in self.atoms['name']]
        
        self.pdb_file = writePDB(self.atoms,pdb_file)
        return self.pdb_file

    @classmethod
    def fromCoarsePDB(cls, pdb_file, dna_type='B_curved', template_from_X3DNA=True, temp_name='temp',
                      compute_topology=True):
        """Initializes a DNA object from a pdb file containing the Coarse Grained atoms"""
        self = cls()

        self.atoms = parsePDB(pdb_file)
        self.atoms.loc[:, 'type'] = self.atoms['name']
        # Initialize the system from the pdb
        self.DNAtype = dna_type
        if compute_topology:
            self.parseConfigurationFile()
            self.computeTopology(temp_name=temp_name, template_from_X3DNA=template_from_X3DNA)
        self.pdb_file = pdb_file
        return self

    @classmethod
    def fromPDB(cls, pdb_file, dna_type='B_curved', template_from_X3DNA=True, output_pdb='clean.pdb', temp_name='temp',
                compute_topology=True):
        """Creates a DNA object from a complete(atomistic) pdb file"""
        self = cls()
        pdb = fixPDB(pdb_file)
        pdb_table = pdb2table(pdb)

        self.atoms = self.CoarseGrain(pdb_table)
        self.DNAtype = dna_type
        if compute_topology:
            self.parseConfigurationFile()
            self.computeTopology(temp_name=temp_name, template_from_X3DNA=template_from_X3DNA)
        self.writePDB(output_pdb)
        #self.atomistic_model=temp
        return self

    @staticmethod
    def CoarseGrain(pdb_table, rna=False):
        """ Selects DNA atoms from a pdb table and returns a table containing only the coarse-grained atoms for 3SPN2"""
        masses = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "P": 30.973762, }
        CG = {"O5\'": 'P', "C5\'": 'S', "C4\'": 'S', "O4\'": 'S', "C3\'": 'S', "O3\'": 'P',
              "C2\'": 'S', "C1\'": 'S', "O5*": 'P', "C5*": 'S', "C4*": 'S', "O4*": 'S',
              "C3*": 'S', "O3*": 'P', "C2*": 'S', "C1*": 'S', "N1": 'B', "C2": 'B', "O2": 'B',
              "N2": 'B', "N3": 'B', "C4": 'B', "N4": 'B', "C5": 'B', "C6": 'B', "N9": 'B',
              "C8": 'B', "O6": 'B', "N7": 'B', "N6": 'B', "O4": 'B', "C7": 'B', "P": 'P',
              "OP1": 'P', "OP2": 'P', "O1P": 'P', "O2P": 'P', "OP3": 'P', "HO5'": 'P',
              "H5'": 'S', "H5''": 'S', "H4'": 'S', "H3'": 'S', "H2'": 'S', "H2''": 'S',
              "H1'": 'S', "H8": 'B', "H61": 'B', "H62": 'B', 'H2': 'B', 'H1': 'B', 'H21': 'B',
              'H22': 'B', 'H3': 'B', 'H71': 'B', 'H72': 'B', 'H73': 'B', 'H6': 'B', 'H41': 'B',
              'H42': 'B', 'H5': 'B', "HO3'": 'P'}
        if rna:
            CG.update({"HO2\'": 'S', "O2\'": 'S'})
        cols = ['recname', 'serial', 'name', 'altLoc',
                'resname', 'chainID', 'resSeq', 'iCode',
                'x', 'y', 'z', 'occupancy', 'tempFactor',
                'element', 'charge', 'type']
        temp = pdb_table.copy()

        # Select DNA residues
        if rna:
            temp = temp[temp['resname'].isin(['A', 'U', 'G', 'C'])]
        else:
            temp = temp[temp['resname'].isin(['DA', 'DT', 'DG', 'DC'])]

        # Group the atoms by sugar, phosphate or base
        temp['group'] = temp['name'].replace(CG)
        temp = temp[temp['group'].isin(['P', 'S', 'B'])]

        # Move the O3' to the next residue
        for c in temp['chainID'].unique():
            sel = temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resSeq"]
            temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resSeq"] = list(sel)[1:] + [-1]
            sel = temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resname"]
            temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resname"] = list(sel)[1:] + ["remove"]
        #temp = temp[temp['resSeq'] > 0]
        temp = temp[temp['resname'] != 'remove']

        # Calculate center of mass
        temp['element'] = temp['element'].str.strip()
        temp['mass'] = temp.element.replace(masses).astype(float)
        temp[['x', 'y', 'z']] = (temp[['x', 'y', 'z']].T * temp['mass']).T[['x', 'y', 'z']]
        temp = temp[temp['element'] != 'H']  # Exclude hydrogens
        Coarse = temp.groupby(['chainID', 'resSeq', 'resname', 'group']).sum().reset_index()
        Coarse[['x', 'y', 'z']] = (Coarse[['x', 'y', 'z']].T / Coarse['mass']).T[['x', 'y', 'z']]

        # Set pdb columns
        Coarse['recname'] = 'ATOM'
        Coarse['name'] = Coarse['group']
        Coarse['altLoc'] = ''
        Coarse['iCode'] = ''
        Coarse['charge'] = ''
        # Change name of base to real base
        mask = (Coarse.name == 'B')
        Coarse.loc[mask, 'name'] = Coarse[mask].resname.str[-1]  # takes last letter from the residue name
        Coarse['type'] = Coarse['name']
        # Set element (depends on base)
        if rna:
            Coarse['name'] = Coarse['name'].replace({'S', 'Sr'})
        Coarse['element'] = Coarse['name'].replace({'P': 'P', 'S': 'H', 'Sr': 'D',
                                                    'A': 'N', 'T': 'S', 'U': 'S', 'G': 'C', 'C': 'O'})
        # Remove P from the beginning
        drop_list = []
        for chain in Coarse.chainID.unique():
            sel = Coarse[Coarse.chainID == chain]
            drop_list += list(sel[(sel.resSeq == sel.resSeq.min()) & sel['name'].isin(['P'])].index)
        Coarse = Coarse.drop(drop_list)

        # Renumber
        Coarse.index = range(len(Coarse))
        Coarse['serial'] = Coarse.index
        return Coarse[cols]

#    @classmethod
#    def fromGRO(cls, gro_file):
#        """Initializes a DNA object from a gromacs input file"""
#        # Parse the gromacs file
#
#        # Make a clean pdb file
#
#        # Initialize the system from the pdb
#        pass

    @classmethod
    def fromSequence(cls, sequence, dna_type='B_curved', output_pdb='clean.pdb', temp_name='temp',
                     compute_topology=True):
        """ Initializes a DNA object from a DNA sequence """
        self = cls()
        self.parseConfigurationFile()
        sequence = pandas.Series([a for a in sequence], index=[('A', i) for i in range(len(sequence))])
        # Make a possible structure
        self.computeGeometry(sequence, temp_name=temp_name)
        self.DNAtype=dna_type

        # Make a clean pdb file
        self = self.fromPDB(f'{temp_name}_template.pdb', dna_type=dna_type, output_pdb=output_pdb,
                            compute_topology=compute_topology)
        return self

    @classmethod
    def fromXYZ(cls, xyz_file, dnatype='B_curved', template_from_X3DNA=True, output_pdb='clean.pdb', temp_name='temp',
                compute_topology=True):
        """ Initializes DNA object from xyz file (as seen on the examples) """
        # Parse the file
        self = cls()
        self.atoms = pandas.read_csv(xyz_file, delim_whitespace=True, skiprows=2,
                                     names=['name', 'x', 'y', 'z'])

        # Compute residues and chains
        residue = 0
        residues = []
        chain = 0
        chains = []
        sugar = True
        for t in self.atoms['name']:
            if t == 'P':
                residue += 1
                sugar = False
            elif t == 'S':
                if sugar:
                    residue += 1
                    chain += 1
                sugar = True
            residues += [residue]
            chains += [chain]
        self.atoms['resSeq'] = residues
        self.atoms['chainID'] = chains

        # compute resid and resname fields
        res_ix = {}
        # min_res = self.atoms.groupby('chain')['resSeq'].min()
        # max_res = self.atoms.groupby('chain')['resSeq'].max()
        for i, res in self.atoms[~self.atoms['name'].isin(['S', 'P'])].iterrows():
            resname = 'D' + res['name']
            # if res['resSeq'] == min_res[res['chain']]:
            #    resname += 'i'
            # if res['resSeq'] == max_res[res['chain']]:
            #    resname += 'f'
            res_ix.update({(res['chainID'], res['resSeq']): resname})
        self.atoms['resname'] = [res_ix[(r['chainID'], r['resSeq'])] for i, r in self.atoms.iterrows()]
        self.DNAtype = dnatype
        if compute_topology:
            self.parseConfigurationFile()
            self.computeTopology(temp_name=temp_name, template_from_X3DNA=template_from_X3DNA)
        self.writePDB(output_pdb)
        return self


class System(openmm.System):
    """ Wrapper of openmm system class, adds some openmm simulation attributes"""

    def __init__(self, dna, forcefieldFiles=[__location__/'3SPN2.xml'], periodicBox=None):
        self.dna = dna
        self.top = openmm.app.PDBFile(dna.pdb_file).getTopology()
        if self.top.getUnitCellDimensions() is None:
            x = dna.atoms[['x', 'y', 'z']]
            d = np.round((x.max() - x.min()) * 2 + 5, -1)
            self.top.setUnitCellDimensions(d)

        self.coord = openmm.app.PDBFile(dna.pdb_file)
        self.forcefield = openmm.app.ForceField(*forcefieldFiles)
        self._wrapped_system = self.forcefield.createSystem(self.top)
        self.periodicBox = periodicBox
        if periodicBox is not None:
            self._wrapped_system.setDefaultPeriodicBoxVectors(*np.diag(self.periodic_box))
            self.dna.periodic = True
        elif self.dna.periodic == True:
            self.dna.periodic = False
            print('Periodic boundary conditions not defined, system will be non periodic')
        self.forces = {}

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self._wrapped_system, attr)

    def clearForces(self, keepCMMotionRemover=True):
        """ Removes all the forces from the system.
        openMM usually adds a "CMMotionRemover" force to keep the center of mass of the system from drifting."""
        j = 0
        for i, f in enumerate(self.getForces()):
            if keepCMMotionRemover and i == 0 and f.__class__ == openmm.CMMotionRemover:
                # print('Kept ', f.__class__)
                j += 1
                continue
            else:
                # print('Removed ', f.__class__)
                self.removeForce(j)
        if keepCMMotionRemover == False:
            assert len(self.getForces()) == 0, 'Not all the forces were removed'
        else:
            assert len(self.getForces()) <= 1, 'Not all the forces were removed'

    def add3SPN2forces(self, verbose=False):
        """ Adds all DNA forces"""
        for force_name in forces:
            if verbose:
                print(force_name)
            force = forces[force_name](self.dna)
            if force_name in ['BasePair', 'CrossStacking']:
                force.addForce(self)
            else:
                self.addForce(force)
            self.forces.update({force_name: force})

    def addProteinDNAforces(self, verbose=False):
        """ Adds protein - DNA interaction forces"""
        for force_name in protein_dna_forces:
            if verbose:
                print(force_name)
            force = forces[force_name](self.dna)
            self.addForce(force)
            self.forces.update({force_name: force})

    def initializeMD(self, temperature=300 * unit.kelvin, platform_name='Reference', damping=2/unit.picosecond, timestep=2*unit.femtoseconds):
        """Starts a simple simulation using the selected system"""
        self.integrator = openmm.LangevinIntegrator(temperature, damping, timestep)
        self.platform = openmm.Platform.getPlatformByName(platform_name)
        self.simulation = openmm.app.Simulation(self.top, self._wrapped_system, self.integrator, self.platform)
        self.simulation.context.setPositions(self.coord.positions)
        return self.simulation

    def setPositions(self, coords=None):
        """Sets the particle positions in the simulation"""
        # Initialize trial MD if not setup
        try:
            self.simulation
        except AttributeError:
            self.initializeMD()

        # Set up coords for MD
        if coords is None:
            self.simulation.context.setPositions(self.coord.positions)
        else:
            self.simulation.context.setPositions(coords)

    def getPotentialEnergy(self, coords=None, energy_unit=unit.kilojoule_per_mole):
        """Returns the potential energy of the current state of the system (default unit KJ/mol)"""
        # Initialize trial MD if not setup
        try:
            self.simulation
        except AttributeError:
            self.initializeMD()
        self.setPositions(coords)
        state = self.simulation.context.getState(getEnergy=True)
        return state.getPotentialEnergy().value_in_unit(energy_unit)

    def recomputeEnergy(self, trajectory, platform_name='Reference'):
        """Returns the potential energy of each snapshot in a xyz trajectory"""
        traj = parse_xyz(trajectory)
        self.initializeMD(platform_name=platform_name)
        energies = []
        for time, snapshot in traj.groupby('timestep'):
            energy = self.getPotentialEnergy(np.array(snapshot[['x', 'y', 'z']]) * _df)
            energies += [energy]
        return np.array(energies)

# List forces
forces = dict(Bond=Bond,
              Angle=Angle,
              Stacking=Stacking,
              Dihedral=Dihedral,
              BasePair=BasePair,
              CrossStacking=CrossStacking,
              Exclusion=Exclusion,
              Electrostatics=Electrostatics)

protein_dna_forces=dict(ExclusionProteinDNA=ExclusionProteinDNA,
                        ElectrostaticsProteinDNA=ElectrostaticsProteinDNA)





# Some typical errors
class BaseError(Exception):
    pass


class DNATypeError(BaseError):
    """Only some DNA types are defined (A, B, B_curved)"""

    def __init__(self, dna):
        self.dna = dna
        self.message = f'DNA type {dna.DNAtype} not defined in the configuration file\n'
        defined_types = dna.angle_definition['DNA'].unique()
        self.message += f'Only the types {str(defined_types)} were defined'
        print(self.message)

class X3DNAnotFound(BaseError):
    """Only some DNA types are defined (A, B, B_curved)"""

    def __init__(self):
        self.message = f'The $X3DNA variable not found in the environment.\n Make sure X3DNA is installed and the environment ' \
                        f'variable $X3DNA is defined.'
        print(self.message)
