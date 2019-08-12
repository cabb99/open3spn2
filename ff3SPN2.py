#!/usr/bin/python3
"""
DNA simulation using the 3SPN2 and 3SPN2.C forcefields in openmm.
"""
# TODO Curved BDNA is currently using the original pdb template, it should use X3DNA if possible,
#  so I need to make it able to take parameters from another pdb if necessary

# TODO put tests in their own folder

__author__ = 'Carlos Bueno'
__version__ = '0.3'

import simtk.openmm.app
import simtk.openmm
import simtk.unit as unit
import configparser
import numpy as np
import itertools
import scipy.spatial.distance as sdist
import os
import simtk.openmm.app
import pdbfixer
import pandas
import nose

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
_ef = 1 * unit.kilocalorie / unit.kilojoule  # energy scaling factor
_df = 1 * unit.angstrom / unit.nanometer  # distance scaling factor
_af = 1 * unit.degree / unit.radian  # angle scaling factor
_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
_dnaResidues = ['DA', 'DC', 'DT', 'DG']
_proteinResidues = ['IPR', 'IGL', 'NGP']
xml = f'{__location__}/3SPN2.xml'


def parseConfigTable(config_section):
    """ Parses a section of the configuration file as a table.
        This function is used to parse the 3SPN2.conf file"""

    def readData(config_section, a):
        """Filters comments and returns values as a list"""
        temp = config_section.get(a).split('#')[0].split()
        l = []
        for val in temp:
            val = val.strip()
            try:
                x = int(val)
                l += [x]
            except ValueError:
                try:
                    y = float(val)
                    l += [y]
                except ValueError:
                    l += [val]
        return l

    data = []
    for a in config_section:
        if a == 'name':
            columns = readData(config_section, a)
        elif len(a) > 3 and a[:3] == 'row':
            data += [readData(config_section, a)]
        else:
            print(f'Unexpected row {readData(config_section, a)}')
    return pandas.DataFrame(data, columns=columns)


def parsePDB(pdb_file):
    """ Transforms the pdb file into a pandas table for easy access and data editing"""

    def pdb_line(line):
        return dict(recname=str(line[0:6]).strip(),
                    serial=int(line[6:11]),
                    name=str(line[12:16]).strip(),
                    altLoc=str(line[16:17]),
                    resname=str(line[17:20]).strip(),
                    chainID=str(line[21:22]),
                    resSeq=int(line[22:26]),
                    iCode=str(line[26:27]),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    occupancy=float(line[54:60]),
                    tempFactor=float(line[60:66]),
                    element=str(line[76:78]),
                    charge=str(line[78:80]))

    with open(pdb_file, 'r') as pdb:
        lines = []
        for line in pdb:
            if len(line) > 6 and line[:6] in ['ATOM  ', 'HETATM']:
                lines += [pdb_line(line)]
    pdb_atoms = pandas.DataFrame(lines)
    pdb_atoms = pdb_atoms[['recname', 'serial', 'name', 'altLoc',
                           'resname', 'chainID', 'resSeq', 'iCode',
                           'x', 'y', 'z', 'occupancy', 'tempFactor',
                           'element', 'charge']]
    return pdb_atoms


def fixPDB(pdb_file):
    """Uses the pdbfixer library to fix a pdb file, replacing non standard residues, removing
    hetero-atoms and adding missing hydrogens. The input is a pdb file location,
    the output is a fixer object, which is a pdb in the openawsem format."""
    fixer = pdbfixer.PDBFixer(filename=pdb_file, )
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    for key in list(keys):
        chain_tmp = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain_tmp.residues())):
            del fixer.missingResidues[key]

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer


def pdb2table(pdb):
    """ Parses a pdb in the openmm format and outputs a table that contains all the information
    on a pdb file """
    cols = ['recname', 'serial', 'name', 'altLoc',
            'resname', 'chainID', 'resSeq', 'iCode',
            'x', 'y', 'z', 'occupancy', 'tempFactor',
            'element', 'charge']
    data = []
    for atom, pos in zip(pdb.topology.atoms(), pdb.positions):
        residue = atom.residue
        chain = residue.chain
        pos = pos.value_in_unit(simtk.unit.angstrom)
        data += [dict(zip(cols, ['ATOM', int(atom.id), atom.name, '',
                                 residue.name, chain.id, int(residue.id), '',
                                 pos[0], pos[1], pos[2], 0, 0,
                                 atom.element.symbol, '']))]
    atom_list = pandas.DataFrame(data)
    atom_list = atom_list[cols]
    atom_list.index = atom_list['serial']
    return atom_list


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


class DNA(object):
    def __init__(self, periodic=True):
        """Initializes an DNA object"""
        self.periodic = periodic

    def __repr__(self):
        return f'3SPN2 DNA object ({len(self.atoms)} atoms)'
        # print the sequence and the identity of the DNA object

    def parseConfigurationFile(self, configuration_file=f'{__location__}/3SPN2.conf'):
        """Reads the configuration file for the forcefield"""
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


    def computeTopology(self):
        """Creates tables of bonds, angles and dihedrals with their respective parameters (bonded interactions)"""
        # Parse configuration file if not already done
        try:
            
            self.bond_definition
        except AttributeError:
            self.parseConfigurationFile()

        DNAtype = self.DNAtype
        if DNAtype not in self.angle_definition['DNA'].unique():
            raise DNATypeError(self)

        #Rewrite index in case it is not ordered
        self.atoms.index = range(len(self.atoms))

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
            x1 = self.atoms.loc[self.bonds['aai']][['x', 'y', 'z']]
            x2 = self.atoms.loc[self.bonds['aaj']][['x', 'y', 'z']]
            self.bonds['r0'] = np.diag(sdist.cdist(x1, x2))

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
                if k1 in index and k2 in index and k3 in index and k4 in index:
                    if b1 == '*' or base[index[k1]] == b1:
                        if b2 == '*' or base[index[k4]] == b2:
                            data += [[i, index[k1], index[k2], index[k3], index[k4], sb]]
        data = pandas.DataFrame(data, columns=['name', 'aai', 'aaj', 'aak', 'aax', 'sB'])
        self.angles = data.merge(angle_types, left_on='name', right_index=True)

        if DNAtype == 'B_curved':
            # Make initial angles the default angles
            v1 = self.atoms.loc[self.angles['aai']].reset_index()[['x', 'y', 'z']]
            v2 = self.atoms.loc[self.angles['aaj']].reset_index()[['x', 'y', 'z']]
            v3 = self.atoms.loc[self.angles['aak']].reset_index()[['x', 'y', 'z']]
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

            a1 = self.atoms.loc[self.dihedrals['aai']].reset_index()[['x', 'y', 'z']]
            a2 = self.atoms.loc[self.dihedrals['aaj']].reset_index()[['x', 'y', 'z']]
            a3 = self.atoms.loc[self.dihedrals['aak']].reset_index()[['x', 'y', 'z']]
            a4 = self.atoms.loc[self.dihedrals['aal']].reset_index()[['x', 'y', 'z']]

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
        """ Writes a minimal version of the pdb file needed for openmm """
        # Compute chain field
        if type(self.atoms['chainID'].iloc[0]) is not str:
            chain_ix = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
            self.atoms['chainID'] = [chain_ix[i - 1] for i in self.atoms['chainID']]

        # Compute element fields
        element_ix = {'P': 'P', 'S': 'H', 'A': 'N', 'T': 'S', 'C': 'O', 'G': 'C'}  # Elements choosen to keep VMD colors
        self.atoms['element'] = [element_ix[atomType] for atomType in self.atoms['name']]

        # Write pdb file
        with open(pdb_file, 'w+') as pdb:
            for i, atom in self.atoms.iterrows():
                pdb_line = f'ATOM  {i + 1:>5} {atom["name"]:^4} {atom.resname:<3}  {atom.chainID}{atom.resSeq:>3}    {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}' + ' ' * 22 + f'{atom.element:2}' + ' ' * 2
                assert len(pdb_line) == 80, 'An item in the atom table is longer than expected'
                pdb.write(pdb_line + '\n')
        self.pdb_file = pdb_file

    @classmethod
    def fromCoarsePDB(cls, pdb_file, dna_type='B_curved'):
        """Initializes a DNA object from a pdb file"""
        self = cls()

        def pdb_line(line):
            return dict(recname=str(line[0:6]).strip(),
                        serial=int(line[6:11]),
                        name=str(line[12:16]).strip(),
                        altLoc=str(line[16:17]),
                        resname=str(line[17:20]).strip(),
                        chainID=str(line[21:22]),
                        resSeq=int(line[22:26]),
                        iCode=str(line[26:27]),
                        x=float(line[30:38]),
                        y=float(line[38:46]),
                        z=float(line[46:54]),
                        occupancy=float(line[54:60]),
                        tempFactor=float(line[60:66]),
                        element=str(line[76:78]).strip(),
                        charge=str(line[78:80]).strip())

        with open(pdb_file, 'r') as pdb:
            lines = []
            for line in pdb:
                if len(line) > 6 and line[:6] in ['ATOM  ', 'HETATM']:
                    lines += [pdb_line(line)]
        pdb_atoms = pandas.DataFrame(lines)
        self.atoms = pdb_atoms[['recname', 'serial', 'name', 'altLoc',
                                'resname', 'chainID', 'resSeq', 'iCode',
                                'x', 'y', 'z', 'occupancy', 'tempFactor',
                                'element', 'charge']]
        # print(self.atoms.columns)
        # self.atoms.loc[:, 'chain'] = self.atoms['chainID']
        # self.atoms.loc[:, 'residue'] = self.atoms['resSeq']
        self.atoms.loc[:, 'type'] = self.atoms['name']
        # print(self.atoms.columns)
        # Initialize the system from the pdb
        self.DNAtype = dna_type
        self.parseConfigurationFile()
        self.computeTopology()
        self.pdb_file = pdb_file
        return self

    @classmethod
    def fromPDB(cls, pdb_file, dna_type='B_curved'):
        """Creates a DNA object from a complete(atomistic) pdb file"""
        self = cls()
        pdb = fixPDB(pdb_file)
        pdb_table = pdb2table(pdb)
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
        cols = ['recname', 'serial', 'name', 'altLoc',
                'resname', 'chainID', 'resSeq', 'iCode',
                'x', 'y', 'z', 'occupancy', 'tempFactor',
                'element', 'charge', 'type']
        temp = pdb_table.copy()

        # Select DNA residues
        temp = temp[temp['resname'].isin(['DA', 'DT', 'DG', 'DC'])]

        # Groupd the atoms by sugar, phosphate or base
        temp['group'] = temp.name.replace(CG)
        temp = temp[temp['group'].isin(['P', 'S', 'B'])]

        # Calculate center of mass
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
        Coarse['element'] = Coarse['name'].replace({'P': 'P', 'S': 'H', 'A': 'N', 'T': 'S', 'G': 'C', 'C': 'O'})
        # Remove P from the beggining
        drop_list = []
        for chain in Coarse.chainID.unique():
            sel = Coarse[Coarse.chainID == chain]
            drop_list += list(sel[(sel.resSeq == sel.resSeq.min()) & sel['name'].isin(['P'])].index)
        Coarse = Coarse.drop(drop_list)
        # Renumber
        Coarse.index = range(len(Coarse))
        Coarse['serial'] = Coarse.index
        self.atoms = Coarse[cols]
        self.DNAtype = dna_type
        self.parseConfigurationFile()
        self.computeTopology()
        self.pdb_file = pdb_file
        return self

    @classmethod
    def fromGRO(cls, gro_file):
        """Initializes a DNA object from a gromacs input file"""
        # Parse the gromacs file

        # Make a clean pdb file

        # Initialize the system from the pdb

        pass

    @classmethod
    def fromSequence(cls, sequence, center=(0, 0, 0)):
        """Initializes a DNA object from a DNA sequence"""
        # Make a possible structure

        # Make a clean pdb file

        # Initialize the system from the pdb

        pass

    @classmethod
    def fromXYZ(cls, xyz_file, dnatype='B_curved'):
        """Initializes DNA object from xyz file (as seen on the examples)"""
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
            if t == 'S':
                if sugar:
                    residue += 1
                    chain += 1
                sugar = True
            if t == 'P':
                residue += 1
                sugar = False
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
        self.parseConfigurationFile()
        self.computeTopology()
        self.writePDB()
        return self


class System(simtk.openmm.System):
    """ Wrapper of openmm system class, adds some openmm simulation attributes"""

    def __init__(self, dna, forcefieldFiles=[f'{__location__}/3SPN2.xml'], periodicBox=None):
        self.dna = dna
        self.top = simtk.openmm.app.PDBFile(dna.pdb_file)
        self.coord = simtk.openmm.app.PDBFile(dna.pdb_file)
        self.forcefield = simtk.openmm.app.ForceField(*forcefieldFiles)
        self._wrapped_system = self.forcefield.createSystem(self.top.topology)
        self.periodicBox = periodicBox
        if periodicBox is not None:
            self._wrapped_system.setDefaultPeriodicBoxVectors(*np.diag(self.periodic_box))

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self._wrapped_system, attr)

    def clearForces(self, keepCMMotionRemover=True):
        """ Removes all the forces from the system """
        j = 0
        for i, f in enumerate(self.getForces()):
            if keepCMMotionRemover and i == 0 and f.__class__ == simtk.openmm.CMMotionRemover:
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

    def add3SPN2forces(self):
        forces = dict(Bond=Bond3SPN2,
                      Angle=Angle3SPN2,
                      Stacking=Stacking3SPN2,
                      Dihedral=Dihedral3SPN2,
                      BasePair=BasePair3SPN2,
                      CrossStacking=CrossStacking3SPN2,
                      Exclusion=Exclusion3SPN2,
                      Electrostatics=Electrostatics3SPN2)
        for force_name in forces:

            force = forces[force_name](self.dna)
            if force_name == 'CrossStacking':
                force.addForce(self)
            else:
                self.addForce(force)

    def initializeMD(self, temperature=300 * unit.kelvin, platform_name='Reference'):
        """Starts a sample simulation using the selected system"""
        self.integrator = simtk.openmm.LangevinIntegrator(temperature, 1E-4 / unit.picosecond, 2 * unit.femtoseconds)
        platform = simtk.openmm.Platform.getPlatformByName(platform_name)
        self.simulation = simtk.openmm.app.Simulation(self.top.topology, self._wrapped_system, self.integrator,
                                                      platform)
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
        energy = state.getPotentialEnergy().value_in_unit(energy_unit)
        return energy

    def recomputeEnergy(self, trajectory, platform_name='Reference'):
        """Returns the potential energy of each snapshot in a trajectory"""
        traj = parse_xyz(trajectory)
        self.initializeMD(platform_name=platform_name)
        energies = []
        for time, snapshot in traj.groupby('timestep'):
            energy = self.getPotentialEnergy(np.array(snapshot[['x', 'y', 'z']]) * _df)
            energies += [energy]
        return np.array(energies)


class Force(object):

    def __init__(self, dna):
        self.periodic = dna.periodic
        self.force = None
        self.dna = dna
        # The patch allows the crosstacking force to run in OpenCL
        # introducing a small difference in the crosstacking energy
        self.OpenCLPatch = True

        # Define the dna force
        self.reset()

        # Define the interaction pairs
        self.defineInteraction()

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.force, attr)

    def computeEnergy(self, system, trajectory):
        # Parse trajectory
        traj = parse_xyz('Tests/adna/traj.xyz')

        # clear all forces on the system
        system.clearForces()
        # setup the force
        self.setUpInteraction()
        # for each item of the table:
        # add the force item

        # compute the energy for every frame

        # return a Table with the energy

    def computeSingleEnergy(self, system, trajectory):
        # Parse trajectory
        traj = parse_xyz('Tests/adna/traj.xyz')
        # for each item of the table:
        # clear all forces on the system
        # system.clearForces()
        # setup the force
        # self.setUpInteraction()
        # add the force item

        # compute the energy for every frame

        # return a table with the energy
    def reset(self):
        pass
    def defineInteraction(self):
        pass

class Bond3SPN2(Force, simtk.openmm.CustomBondForce):
    def getParameterNames(self):
        self.perInteractionParameters = []
        self.GlobalParameters = []
        for i in range(self.force.getNumPerBondParameters()):
            self.perInteractionParameters += [self.force.getPerBondParameterName(i)]
        for i in range(self.force.getNumGlobalParameters()):
            self.GlobalParameters += [self.force.getGlobalParameterName(i)]
        return [self.perInteractionParameters, self.GlobalParameters]

    def reset(self):
        bondForce = simtk.openmm.CustomBondForce("Kb2*(r-r0)^2+Kb3*(r-r0)^3+Kb4*(r-r0)^4")
        bondForce.addPerBondParameter('r0')
        bondForce.addPerBondParameter('Kb2')
        bondForce.addPerBondParameter('Kb3')
        bondForce.addPerBondParameter('Kb4')
        bondForce.setUsesPeriodicBoundaryConditions(self.periodic)
        bondForce.setForceGroup(3)
        self.force = bondForce

    def defineInteraction(self):
        for i, b in self.dna.bonds.iterrows():
            # Units converted from
            parameters = [b['r0'] * _df,
                          b['Kb2'] / _df ** 2 * _ef,
                          b['Kb3'] / _df ** 3 * _ef,
                          b['Kb4'] / _df ** 4 * _ef]
            self.force.addBond(int(b['aai']), int(b['aaj']), parameters)


class Angle3SPN2(Force, simtk.openmm.HarmonicAngleForce):

    def reset(self):
        angleForce = simtk.openmm.HarmonicAngleForce()
        angleForce.setUsesPeriodicBoundaryConditions(self.periodic)
        angleForce.setForceGroup(4)
        self.force = angleForce

    def defineInteraction(self):
        for i, a in self.dna.angles.iterrows():
            parameters = [a['t0'] * _af,
                          a['epsilon'] * 2]
            self.force.addAngle(int(a['aai']), int(a['aaj']), int(a['aak']), *parameters)


class Stacking3SPN2(Force, simtk.openmm.CustomCompoundBondForce):
    def reset(self):
        stackingForce = simtk.openmm.CustomCompoundBondForce(3, """rep+f2*attr;
                                                                rep=epsilon*(1-exp(-alpha*(dr)))^2*step(-dr);
                                                                attr=epsilon*(1-exp(-alpha*(dr)))^2*step(dr)-epsilon;
                                                                dr=distance(p2,p3)-sigma;
                                                                f2=max(f*pair2,pair1);
                                                                pair1=step(dt+pi/2)*step(pi/2-dt);
                                                                pair2=step(dt+pi)*step(pi-dt);
                                                                f=1-cos(dt)^2;
                                                                dt=rng*(angle(p1,p2,p3)-t0);""")
        stackingForce.setUsesPeriodicBoundaryConditions(self.periodic)
        stackingForce.addPerBondParameter('epsilon')
        stackingForce.addPerBondParameter('sigma')
        stackingForce.addPerBondParameter('t0')
        stackingForce.addPerBondParameter('alpha')
        stackingForce.addPerBondParameter('rng')
        stackingForce.addGlobalParameter('pi', np.pi)
        stackingForce.setForceGroup(5)
        self.force = stackingForce

    def defineInteraction(self):
        for i, a in self.dna.stackings.iterrows():
            parameters = [a['epsilon'] * _ef,
                          a['sigma'] * _df,
                          a['t0'] * _af,
                          a['alpha'] / _df,
                          a['rng']]
            self.force.addBond([a['aai'], a['aaj'], a['aak']], parameters)


class Dihedral3SPN2(Force, simtk.openmm.CustomTorsionForce):
    def reset(self):
        dihedralForce = simtk.openmm.CustomTorsionForce("""K_periodic*(1-cs)-K_gaussian*exp(-atan(tan(dt))^2/2/sigma^2);
                                                      cs=cos(dt);
                                                      dt=theta-t0""")
        # dihedralForce=simtk.openmm.CustomTorsionForce("theta/60.")
        dihedralForce.setUsesPeriodicBoundaryConditions(self.periodic)
        dihedralForce.addPerTorsionParameter('K_periodic')
        dihedralForce.addPerTorsionParameter('K_gaussian')
        dihedralForce.addPerTorsionParameter('sigma')
        dihedralForce.addPerTorsionParameter('t0')
        dihedralForce.addGlobalParameter('pi', np.pi)
        dihedralForce.setForceGroup(6)
        self.force = dihedralForce

    def defineInteraction(self):
        for i, a in self.dna.dihedrals.iterrows():
            parameters = [a['K_dihedral'] * _ef,
                          a['K_gaussian'] * _ef,
                          a['sigma'],
                          (180 + a['t0']) * _af]
            particles = [a['aai'], a['aaj'], a['aak'], a['aal']]
            self.force.addTorsion(*particles, parameters)


class BasePair3SPN2(Force, simtk.openmm.CustomHbondForce):
    def reset(self):
        def basePairForce():
            pairForce = simtk.openmm.CustomHbondForce('''temp;temp=rep+1/2*(1+cos(dphi))*fdt1*fdt2*attr;
                                                         rep  = epsilon*(1-exp(-alpha*dr))^2*(1-step(dr));
                                                         attr = epsilon*(1-exp(-alpha*dr))^2*step(dr)-epsilon;
                                                         fdt1 = max(f1*pair0t1,pair1t1);
                                                         fdt2 = max(f2*pair0t2,pair1t2);
                                                         pair1t1 = step(pi/2+dt1)*step(pi/2-dt1);
                                                         pair1t2 = step(pi/2+dt2)*step(pi/2-dt2);
                                                         pair0t1 = step(pi+dt1)*step(pi-dt1);
                                                         pair0t2 = step(pi+dt2)*step(pi-dt2);
                                                         f1 = 1-cos(dt1)^2;
                                                         f2 = 1-cos(dt2)^2;
                                                         dphi = dihedral(d2,d1,a1,a2)-phi0;
                                                         dr    = distance(d1,a1)-sigma;
                                                         dt1   = rng*(angle(d2,d1,a1)-t01);
                                                         dt2   = rng*(angle(a2,a1,d1)-t02);''')
            if self.periodic:
                pairForce.setNonbondedMethod(pairForce.CutoffPeriodic)
            else:
                pairForce.setNonbondedMethod(pairForce.CutoffNonPeriodic)
            pairForce.setCutoffDistance(1.8)  # Paper
            pairForce.addPerDonorParameter('phi0')
            pairForce.addPerDonorParameter('sigma')
            pairForce.addPerDonorParameter('t01')
            pairForce.addPerDonorParameter('t02')
            pairForce.addPerDonorParameter('rng')
            pairForce.addPerDonorParameter('epsilon')
            pairForce.addPerDonorParameter('alpha')
            pairForce.addGlobalParameter('pi', np.pi)
            self.force = pairForce
            pairForce.setForceGroup(7)
            return pairForce

        basePairForces = {}
        pair_definition = self.dna.pair_definition[self.dna.pair_definition['DNA'] == self.dna.DNAtype]
        for i, pair in pair_definition.iterrows():
            basePairForces.update({i: basePairForce()})
        self.forces = basePairForces

    def defineInteraction(self):
        pair_definition = self.dna.pair_definition[self.dna.pair_definition['DNA'] == self.dna.DNAtype]
        is_dna = self.dna.atoms['resname'].isin(_dnaResidues)

        for i, pair in pair_definition.iterrows():
            D1 = list(self.dna.atoms[(self.dna.atoms['name'] == pair['Base1']) & is_dna].index)
            A1 = list(self.dna.atoms[(self.dna.atoms['name'] == pair['Base2']) & is_dna].index)

            D2 = [j - 1 for j in D1]
            A2 = [j - 1 for j in A1]

            # Define parameters
            parameters = [pair.torsion * _af,
                          pair.sigma * _df,
                          pair.t1 * _af,
                          pair.t2 * _af,
                          pair.rang,
                          pair.epsilon * _ef,
                          pair.alpha / _df]

            # Add donors and acceptors
            # Here I am including the same atom twice,
            # it doesn't seem to break things
            for d1, d2 in zip(D1, D2):
                self.forces[i].addDonor(d1, d2, d2, parameters)
            for a1, a2 in zip(A1, A2):
                self.forces[i].addAcceptor(a1, a2, a2)

            # Exclude interactions
            atoms = self.dna.atoms.copy()
            atoms.loc[D1, 'donor_id'] = [i for i in range(len(D1))]
            atoms.loc[A1, 'aceptor_id'] = [i for i in range(len(A1))]

            interaction = np.zeros([len(self.dna.atoms.resSeq.unique()) // 2, len(self.dna.atoms.resSeq.unique()) // 2])
            interaction[:, :] = 1
            for k, chain in atoms.loc[D1 + A1].groupby('chainID'):
                dds = chain.donor_id.dropna()
                aas = chain.aceptor_id.dropna()
                for d, a in itertools.product(dds.iteritems(), aas.iteritems()):
                    di = d[0]
                    d = int(d[1])
                    ai = a[0]
                    a = int(a[1])
                    # The sequence exclusion was reduced to two residues
                    # since the maximum number of exclusions in OpenCL is 4.
                    # In the original 3SPN2 it was 3 residues (6 to 9)
                    # This change has no noticeable effect
                    if abs(ai - di) <= 6:
                        self.forces[i].addExclusion(d, a)
                        # print(di, ai, d, a)
                        interaction[d, a] = 0

    def addForce(self, system):
        for f in self.forces:
            system.addForce(self.forces[f])


class CrossStacking3SPN2(Force):
    def reset(self):
        def crossStackingForce(parametersOnDonor=False):
            if self.OpenCLPatch:
                # Divide the energy over two to make the interaction symetric
                patch_string = '/2'
            else:
                patch_string = ''
            crossForce = simtk.openmm.CustomHbondForce(f'''energy{patch_string};temp=dtCS/rng_CS;
                                                        energy   = fdt3*fdtCS*attr;
                                                        attr     = epsilon*(1-exp(-alpha*dr))^2*step(dr)-epsilon;
                                                        fdt3     = max(f1*pair0t3,pair1t3);
                                                        fdtCS    = max(f2*pair0tCS,pair1tCS);
                                                        pair0t3  = step(pi+dt3)*step(pi-dt3);
                                                        pair0tCS = step(pi+dtCS)*step(pi-dtCS);
                                                        pair1t3  = step(pi/2+dt3)*step(pi/2-dt3);
                                                        pair1tCS = step(pi/2+dtCS)*step(pi/2-dtCS);
                                                        f1       = 1-cos(dt3)^2;
                                                        f2       = 1-cos(dtCS)^2;
                                                        dr       = distance(d1,a3)-sigma;
                                                        dt3      = rng_BP*(t3-t03);
                                                        dtCS     = rng_CS*(tCS-t0CS);
                                                        tCS      = angle(d2,d1,a3);
                                                        t3       = acos(sin(t1)*sin(t2)*cos(phi)-cos(t1)*cos(t2));
                                                        t1       = angle(d2,d1,a1);
                                                        t2       = angle(d1,a1,a2);
                                                        phi      = dihedral(d2,d1,a1,a2);''')
            if self.periodic:
                crossForce.setNonbondedMethod(crossForce.CutoffPeriodic)
            else:
                crossForce.setNonbondedMethod(crossForce.CutoffNonPeriodic)
            crossForce.setCutoffDistance(1.8)  # Paper
            parameters = ['t03', 't0CS', 'rng_CS', 'rng_BP', 'epsilon', 'alpha', 'sigma']
            for p in parameters:
                if parametersOnDonor:
                    crossForce.addPerDonorParameter(p)
                else:
                    crossForce.addPerAcceptorParameter(p)
            crossForce.addGlobalParameter('pi', np.pi)
            crossForce.setForceGroup(8)
            return crossForce

        crossStackingForces = {}
        for base in ['A', 'T', 'G', 'C']:
            crossStackingForces.update({base: (crossStackingForce(), crossStackingForce())})
        self.crossStackingForces = crossStackingForces

    def defineInteraction(self):
        # Donators
        bases = np.array(self.dna.atoms[self.dna.atoms['name'].isin(['A', 'T', 'G', 'C'])].index)
        D1 = self.dna.atoms.reindex(bases)
        D2 = self.dna.atoms.reindex(bases - 1)
        D3 = self.dna.atoms.reindex(bases + 3)
        A1 = self.dna.atoms.reindex(bases)
        A2 = self.dna.atoms.reindex(bases - 1)
        A3 = self.dna.atoms.reindex(bases - 3)

        # Select only bases where the other atoms exist
        D2.index = D1.index
        D3.index = D1.index
        temp = pandas.concat([D1, D2, D3], axis=1, keys=['D1', 'D2', 'D3'])
        sel = np.array(temp[temp['D3', 'name'].isin(['A', 'T', 'G', 'C']) & temp['D2', 'name'].isin(['S'])].index)
        D1 = self.dna.atoms.reindex(sel)
        D2 = self.dna.atoms.reindex(sel - 1)
        D3 = self.dna.atoms.reindex(sel + 3)

        # Aceptors
        A2.index = A1.index
        A3.index = A1.index
        temp = pandas.concat([A1, A2, A3], axis=1, keys=['A1', 'A2', 'A3'])
        sel = np.array(temp[temp['A3', 'name'].isin(['A', 'T', 'G', 'C']) & temp['A2', 'name'].isin(['S'])].index)
        A1 = self.dna.atoms.reindex(sel)
        A2 = self.dna.atoms.reindex(sel - 1)
        A3 = self.dna.atoms.reindex(sel - 3)

        # Parameters
        cross_definition = self.dna.cross_definition[self.dna.cross_definition['DNA'] == self.dna.DNAtype].copy()
        i = [a for a in zip(cross_definition['Base_d1'], cross_definition['Base_a1'], cross_definition['Base_a3'])]
        cross_definition.index = i

        donors = {i: [] for i in ['A', 'T', 'G', 'C']}
        for donator, donator2, d1, d2, d3 in zip(D1.itertuples(), D3.itertuples(), D1.index, D2.index, D3.index):
            # if d1!=4:
            #    continue
            d1t = donator.name
            d3t = donator2.name
            c1, c2 = self.crossStackingForces[d1t]
            a1t = _complement[d1t]
            # print(d1, d2, d3)
            param = cross_definition.loc[[(a1t, d1t, d3t)]].squeeze()
            # parameters=[param1['t03']*af,param1['T0CS_1']*af,param1['rng_cs1'],param1['rng_bp'],param1['eps_cs1']*ef,param1['alpha_cs1']/df,param1['Sigma_1']*df]
            parameters = [param['t03'] * _af,
                          param['T0CS_2'] * _af,
                          param['rng_cs2'],
                          param['rng_bp'],
                          param['eps_cs2'] * _ef,
                          param['alpha_cs2'] / _df,
                          param['Sigma_2'] * _df]
            # print(param)
            c1.addDonor(d1, d2, d3)
            c2.addAcceptor(d1, d2, d3, parameters)
            donors[d1t] += [d1]

        aceptors = {i: [] for i in ['A', 'T', 'G', 'C']}
        for aceptor, aceptor2, a1, a2, a3 in zip(A1.itertuples(), A3.itertuples(), A1.index, A2.index, A3.index):
            # if a1!=186:
            #    continue
            a1t = aceptor.name
            a3t = aceptor2.name
            c1, c2 = self.crossStackingForces[_complement[a1t]]
            d1t = _complement[a1t]
            param = cross_definition.loc[[(d1t, a1t, a3t)]].squeeze()
            # print(param)
            # print(a1, a2, a3)
            parameters = [param['t03'] * _af,
                          param['T0CS_1'] * _af,
                          param['rng_cs1'],
                          param['rng_bp'],
                          param['eps_cs1'] * _ef,
                          param['alpha_cs1'] / _df,
                          param['Sigma_1'] * _df]
            # parameters=[param1['t03']*af,param1['T0CS_2']*af,param1['rng_cs2'],param1['rng_bp'],param1['eps_cs2']*ef,param1['alpha_cs2']/df,param1['Sigma_2']*df]
            c1.addAcceptor(a1, a2, a3, parameters)
            c2.addDonor(a1, a2, a3)
            aceptors[_complement[a1t]] += [a1]

        # Exclusions
        for base in ['A', 'T', 'G', 'C']:
            c1, c2 = self.crossStackingForces[base]
            for ii, i in enumerate(donors[base]):
                for jj, j in enumerate(aceptors[base]):
                    # The sequence exclusion was reduced to two residues
                    # since the maximum number of exclusions in OpenCL is 4.
                    # In the original 3SPN2 it was 3 residues (6 to 9)
                    # This change has a small effect in B-DNA and curved B-DNA
                    # The second change is to make the interaction symetric and dividing the energy over 2
                    # This also reduces the number of exclusions in the force
                    if self.OpenCLPatch:
                        maxn = 6
                    else:
                        maxn = 9
                    if (self.dna.atoms.at[i, 'chainID'] == self.dna.atoms.at[j, 'chainID'] and abs(i - j) <= maxn) or \
                            (not self.OpenCLPatch and i > j):
                        c1.addExclusion(ii, jj)
                        c2.addExclusion(jj, ii)

    def addForce(self, system):
        for c1, c2 in self.crossStackingForces.values():
            system.addForce(c1)
            system.addForce(c2)

    def getForceGroup(self):
        fg = 0
        for c1, c2 in self.crossStackingForces.values():
            fg = c1.getForceGroup()
            break
        for c1, c2 in self.crossStackingForces.values():
            assert fg == c1.getForceGroup()
            assert fg == c2.getForceGroup()
        return fg


def addNonBondedExclusions(dna, force):
    is_dna = dna.atoms['resname'].isin(_dnaResidues)
    atoms = dna.atoms.copy()
    selection = atoms[is_dna]
    for (i, atom_a), (j, atom_b) in itertools.product(selection.iterrows(), repeat=2):
        if j >= i:
            continue
        # Neighboring residues
        if atom_a['chainID'] == atom_b['chainID'] and (atom_a.resSeq - atom_b.resSeq <= 1):
            force.addExclusion(i, j)
            #print(i, j)
        # Base-pair residues
        elif (atom_a['name'] in _complement.keys()) and (atom_b['name'] in _complement.keys()) and (
                atom_a['name'] == _complement[atom_b['name']]):
            force.addExclusion(i, j)
            #print(i, j)


class Exclusion3SPN2(Force, simtk.openmm.CustomNonbondedForce):
    def reset(self):
        exclusionForce = simtk.openmm.CustomNonbondedForce("""energy;
                                                            energy=(epsilon*((sigma/r)^12-2*(sigma/r)^6)+epsilon)*step(sigma-r);
                                                            sigma=0.5*(sigma1+sigma2); 
                                                            epsilon=sqrt(epsilon1*epsilon2)""")
        exclusionForce.addPerParticleParameter('epsilon')
        exclusionForce.addPerParticleParameter('sigma')
        exclusionForce.setCutoffDistance(1.8)
        exclusionForce.setForceGroup(9)  # There can not be multiple cutoff distance on the same force group
        if self.periodic:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffPeriodic)
        else:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffNonPeriodic)
        self.force = exclusionForce

    def defineInteraction(self):
        # addParticles
        particle_definition = self.dna.particle_definition[self.dna.particle_definition['DNA'] == self.dna.DNAtype]
        particle_definition.index = particle_definition.name

        # Reduces or increases the cutoff to the maximum particle radius
        self.force.setCutoffDistance(particle_definition.radius.max() * _df)

        # Select only dna atoms
        is_dna = self.dna.atoms['resname'].isin(_dnaResidues)
        atoms = self.dna.atoms.copy()
        atoms['is_dna'] = is_dna
        for i, atom in atoms.iterrows():
            if atom.is_dna:
                param = particle_definition.loc[atom['name']]
                parameters = [param.epsilon * _ef,
                              param.radius * _df]
            else:
                parameters = [0, .1]  # Null energy and some radius)
            # print(i, parameters)
            self.force.addParticle(parameters)

        # addExclusions
        selection = atoms[is_dna]
        for (i, atom_a), (j, atom_b) in itertools.product(selection.iterrows(), repeat=2):
            if j >= i:
                continue
            # Neighboring residues
            if atom_a['chainID'] == atom_b['chainID'] and (atom_a.resSeq - atom_b.resSeq <= 1):
                self.force.addExclusion(i, j)
            # Base-pair residues
            elif (atom_a['name'] in _complement.keys()) and (atom_b['name'] in _complement.keys()) and (
                    atom_a['name'] == _complement[atom_b['name']]):
                self.force.addExclusion(i, j)


class Electrostatics3SPN2(Force, simtk.openmm.CustomNonbondedForce):
    def reset(self):
        T = 300 * unit.kelvin
        C = 100 * unit.millimolar
        e = 249.4 - 0.788 * (T / unit.kelvin) + 7.2E-4 * (T / unit.kelvin) ** 2
        a = 1 - 0.2551 * (C / unit.molar) + 5.151E-2 * (C / unit.molar) ** 2 - 6.889E-3 * (C / unit.molar) ** 3
        print(e, a)
        dielectric = e * a
        # Debye length
        kb = simtk.unit.BOLTZMANN_CONSTANT_kB  # Bolztmann constant
        Na = simtk.unit.AVOGADRO_CONSTANT_NA  # Avogadro number
        ec = 1.60217653E-19 * unit.coulomb  # proton charge
        pv = 8.8541878176E-12 * unit.farad / unit.meter  # dielectric permittivity of vacuum

        ldby = np.sqrt(dielectric * pv * kb * T / (2.0 * Na * ec ** 2 * C))
        ldby = ldby.in_units_of(unit.nanometer)
        denominator = 4 * np.pi * pv * dielectric / Na * ec ** 2
        print(ldby, denominator)

        electrostaticForce = simtk.openmm.CustomNonbondedForce("""energy;
                                                                energy=q1*q2*exp(-r/dh_length)/denominator/r;""")
        electrostaticForce.addPerParticleParameter('q')
        electrostaticForce.addGlobalParameter('dh_length', ldby)
        electrostaticForce.addGlobalParameter('denominator', denominator)

        electrostaticForce.setCutoffDistance(5)
        if self.periodic:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffPeriodic)
        else:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffNonPeriodic)
        electrostaticForce.setForceGroup(10)
        self.force = electrostaticForce

    def defineInteraction(self):
        # addParticles
        particle_definition = self.dna.particle_definition[self.dna.particle_definition['DNA'] == self.dna.DNAtype]
        particle_definition.index = particle_definition.name

        # Select only dna atoms
        is_dna = self.dna.atoms['resname'].isin(_dnaResidues)
        atoms = self.dna.atoms.copy()
        atoms['is_dna'] = is_dna

        for i, atom in atoms.iterrows():
            if atom.is_dna:
                param = particle_definition.loc[atom['name']]
                parameters = [param.charge]
            else:
                parameters = [0]  # No charge if it is not DNA
            # print (i,parameters)
            self.force.addParticle(parameters)

        # add neighbor exclusion
        selection = atoms[is_dna]
        for (i, a1), (j, a2) in itertools.product(selection.iterrows(), repeat=2):
            if j >= i:
                continue
            if (a1.chainID == a2.chainID) and (abs(a1.resSeq - a2.resSeq) <= 1):
                self.force.addExclusion(i, j)
                # print (i,j,a1.chain,a2.chain,a1.resSeq,a2.resSeq,(a1.resSeq-a2.resSeq))
            # All Non-bonded forces must have the same exclusions
            elif self.OpenCLPatch and (a1['name'] in _complement.keys()) and (a2['name'] in _complement.keys()) and (
                    a1['name'] == _complement[a2['name']]):
                self.force.addExclusion(i, j)


class ProteinDNAForce(Force):
    def __init__(self, dna, protein):
        self.protein = protein
        super().__init__(dna)


class ExclusionProteinDNA(ProteinDNAForce,Exclusion3SPN2):
    def reset(self):
        super().reset()
        self.force.setForceGroup(2)

    def defineInteraction(self):
        # Merge DNA and protein particle definitions
        dna_particle_definition = self.dna.particle_definition[self.dna.particle_definition['DNA'] == self.dna.DNAtype].copy()
        dna_particle_definition.loc[:, 'molecule'] = 'DNA'
        protein_particle_definition = self.dna.config['Protein Particles']
        protein_particle_definition.loc[:, 'molecule'] = 'protein'
        particle_definition = pandas.concat([dna_particle_definition, protein_particle_definition], sort=False)
        particle_definition.index = particle_definition.molecule + particle_definition.name

        self.force.setCutoffDistance(particle_definition.radius.max() * _df)
        # Select only dna atoms
        is_dna = self.dna.atoms['resname'].isin(_dnaResidues)
        is_protein = self.dna.atoms['resname'].isin(_proteinResidues)
        atoms = self.dna.atoms.copy()
        atoms['is_dna'] = is_dna
        atoms['is_protein'] = is_protein
        DNA_list = []
        protein_list = []
        for i, atom in atoms.iterrows():
            if atom.is_dna:
                param = particle_definition.loc['DNA' + atom['name']]
                parameters = [param.epsilon * _ef,
                              param.radius * _df]
                DNA_list += [i]
            elif atom.is_protein:
                param = particle_definition.loc['protein' + atom['name']]
                parameters = [param.epsilon * _ef,
                              param.radius * _df]
                protein_list += [i]
            else:
                print(f'Residue {i} not included in protein-DNA interactions')
                parameters = [0, .1]
            self.force.addParticle(parameters)
        self.force.addInteractionGroup(DNA_list, protein_list)


class ElectrostaticsProteinDNA(Electrostatics3SPN2, ProteinDNAForce):
    def reset(self):
        super().reset()
        self.force.setForceGroup(1)

    def defineInteraction(self):
        # Merge DNA and protein particle definitions
        dna_particle_definition = self.dna.particle_definition[self.dna.particle_definition['DNA'] == self.dna.DNAtype].copy()
        dna_particle_definition.loc[:, 'molecule'] = 'DNA'
        protein_particle_definition = self.dna.config['Protein Particles']
        protein_particle_definition.loc[:, 'molecule'] = 'protein'
        particle_definition = pandas.concat([dna_particle_definition, protein_particle_definition], sort=False)
        particle_definition.index = particle_definition.molecule + particle_definition.name

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
            elif atom.is_protein:
                atom_param = particle_definition.loc['protein' + atom['name']]
                seq_param = sequence_electrostatics.loc[atom.real_resname]
                charge = atom_param.charge * seq_param.charge
                parameters = [charge]
                if charge != 0:
                    protein_list += [i]
            else:
                print(f'Residue {i} not included in protein-DNA electrostatics')
                parameters = [0]  # No charge if it is not DNA
            # print (i,parameters)
            self.force.addParticle(parameters)
        self.force.addInteractionGroup(DNA_list, protein_list)


# Unit testing
def test_DNA_from_pdb():
    mol = DNA.fromPDB("Tests/1svc/1svc.pdb")


def test_DNA_from_gro():
    pass


def test_DNA_from_seq():
    pass


def test_DNA_from_xyz():
    """Tests the correct parsing from an xyz file"""
    mol = DNA.fromXYZ('Tests/adna/in00_conf.xyz')
    assert mol.atoms.at[8, 'name'] == 'P'
    assert round(mol.atoms.at[188, 'y'], 6) == -8.779343


# Functional testing

def parse_xyz(filename=''):
    columns = ['N', 'timestep', 'id', 'name', 'x', 'y', 'z']
    data = []
    with open(filename, 'r') as traj_file:
        atom = pandas.Series(index=columns)
        atom['id'] = None
        for line in traj_file:
            s = line.split()
            if len(s) == 1:
                atom['N'] = int(s[0])
                if atom['id'] > -1:
                    assert atom['id'] == atoms
                atoms = int(s[0])
            elif len(s) == 3:
                atom['timestep'] = int(s[2])
                atom['id'] = 0
            elif len(s) > 3:
                atom['name'] = int(s[0])
                atom['x'], atom['y'], atom['z'] = [float(a) for a in s[1:4]]
                data += [atom.copy()]
                atom['id'] += 1
    xyz_data = pandas.concat(data, axis=1).T
    for i in ['N', 'timestep', 'id', 'name']:
        xyz_data[i] = xyz_data[i].astype(int)
    return xyz_data


def parse_log(filename=''):
    columns = ''
    log_data = []
    with open(filename, 'r') as log_file:
        start = False
        for line in log_file:
            if line[:4] == 'Step':
                columns = line.split()
                start = True
                continue
            if start:
                try:
                    log_data += [[float(a) for a in line.split()]]
                except ValueError:
                    break
    log_data = pandas.DataFrame(log_data, columns=columns)

    for i in ['Step', 'nbp']:
        log_data[i] = log_data[i].astype(int)
    return log_data


def test_parse_xyz():
    """Tests the example trajectory parsing"""
    xyz_data = parse_xyz('Tests/adna/traj.xyz')
    assert xyz_data.at[1, 'name'] == 7
    assert xyz_data.at[1, 'x'] == 4.34621


def test_parse_log():
    """Tests the example log parsing"""
    log_data = parse_log('Tests/adna/sim.log')
    assert log_data.at[1, 'Step'] == 2000
    assert log_data.at[1, 'eexcl'] == 0.45734636


class TestEnergies:
    """Tests that the energies are the same as the example outputs from lammps"""
    force3SPN2 = dict(Bond=Bond3SPN2,
                      Angle=Angle3SPN2,
                      Stacking=Stacking3SPN2,
                      Dihedral=Dihedral3SPN2,
                      BasePair=BasePair3SPN2,
                      CrossStacking=CrossStacking3SPN2,
                      Exclusion=Exclusion3SPN2,
                      Electrostatics=Electrostatics3SPN2)

    def _test_energy(self,
                     log_energy='E_bond',
                     log_file='Tests/adna/sim.log',
                     traj_file='Tests/adna/traj.xyz',
                     force='Bond', periodic_size=94.2,
                     platform_name='Reference',dna=None, system=None):
        self.dna=dna
        self.system=system
        self.system.clearForces()
        self.system.setDefaultPeriodicBoxVectors(*np.diag([2 * periodic_size / 10] * 3))

        log = parse_log(log_file)

        self.system.clearForces()
        f = self.force3SPN2[force]
        tempforce = f(self.dna)
        try:
            tempforce.addForce(self.system)
        except AttributeError:
            self.system.addForce(tempforce)
        energies = self.system.recomputeEnergy(traj_file, platform_name=platform_name)
        d = (energies / _ef - log[log_energy])
        diff = np.sqrt((d ** 2).sum() / len(energies))
        print(f'The difference in the energy term {log_energy} is {diff} Kcal/mol')
        print(f'The DNA type of the system is {self.dna.DNAtype}')
        data = np.array([energies / _ef, np.array(log[log_energy]), np.array(d)])
        results = pandas.DataFrame(data.T, columns=['Openmm energy', 'Lammps energy', 'Difference'])
        # print(data)
        print(pandas.DataFrame(data.T, columns=['Openmm energy', 'Lammps energy', 'Difference']))
        # The maximum error seems to be bonds in curved BDNA (0.002)
        assert diff < 2E-3, diff
        assert len(results.dropna()) == len(results), results

    def test_energies(self):
        test_sets = pandas.read_csv('Tests/test_cases.csv', comment='#')
        for i, tests in test_sets.groupby(['Folder', 'DNA type']):
            folder = i[0]
            dna_type = i[1]
            self.dna = DNA.fromXYZ(f'{folder}/in00_conf.xyz', dna_type)
            self.system = System(self.dna)
            for j, test in tests.iterrows():
                print(j)
                yield self._test_energy, test['Energy term'], f'{folder}/{test.Log}', f'{folder}/{test.Trajectory}',\
                      test['Name'], test['periodic size'], test['Platform'], self.dna, self.system

    def _test_energies_slow(self):
        test_sets = pandas.read_csv('Tests/test_cases.csv', comment='#')
        for i, test in test_sets.iterrows():
            dna_type = test['DNA type']
            folder = test['Folder']
            yield self._test_energy, test['Energy term'], f'{folder}/{test.Log}', f'{folder}/{test.Trajectory}', \
                  test['Name'], folder, dna_type, test['periodic size'], test['Platform']
