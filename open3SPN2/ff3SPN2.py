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
from .circular import circular_chain_info, residue_wrapper, CircularBend, kabsch, linking_number
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
    # Chains whose last residue is bonded to the first one (see DNA.circularize)
    circular_chains = ()

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

    def to_pdb(self, pdb_file):
        writePDB(self.atoms,pdb_file)

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

    def computeGeometry(self, sequence=None, atomistic_template=None, Coarse_template=None):
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
            xyz2 = np.dot(xyz2,new_orien)+new_pos
            
            b1[['x','y','z']] = xyz1
            b2[['x','y','z']] = xyz2
            
            for b in [b1,b2]:
                b['recname']='ATOM'
                b['resname']='D'+b['base']
                b['occupancy']=1.00
                b['tempFactor']=1.00
                b['element']=b['name'].str[0]

            b1['resSeq']=row[0]+1
            b2['resSeq']=len(data)-row[0]

            b1['chainID']='A'
            b2['chainID']='B'

            forward_strand=forward_strand + [b1]
            reverse_strand=[b2] + reverse_strand
            

        template=pandas.concat(forward_strand+reverse_strand)
        template.reindex()
        template['serial']=template.index+1
        self.template=template
        if atomistic_template:
            writePDB(template,pdb_file=atomistic_template)
        
        template = self.CoarseGrain(template)
        if Coarse_template:
            writePDB(template,pdb_file=Coarse_template)
        
        try:
            self.atoms
        except AttributeError:
            return template
        
        template = template[template['chainID'] == 'A']
        original = self.atoms.copy()
        original.index = original['chainID'].astype(str) + '_' + original['resSeq'].astype(str) + '_' + original['name']
        ii = []
        keys = []
        for i, r in template.iterrows():
            name = r['name']
            seq_template = r['resname'][1:]
            j = r['resSeq'] - 1
            seq_seq = sequence.iloc[j]
            assert seq_seq == seq_template
            chain, resseq = sequence.index[j]
            ii += [f'{chain}_{resseq}_{name}']
            keys += [(chain, resseq, name)]
        # Template coordinates by (chain, resSeq, name), including the residues that pad circular chains
        self.template_positions = dict(zip(keys, template[['x', 'y', 'z']].values))
        template.index = ii
        merge = pandas.merge(original, template, left_index=True, right_index=True, how='left', suffixes=['_old', ''])
        original[['x', 'y', 'z']] = merge[['x', 'y', 'z']]
        original.index = self.atoms.index
        return original

    def computeTopology(self, template_from_X3DNA=True):
        """ Creates tables of bonds, angles and dihedrals with their respective parameters (bonded interactions).
        3SPN2.C requires a template structure to calculate the equilibrium bonds, angles and dihedrals.
        If template_from_structure is True, it will try to compute the equilibrium geometry using X3DNA.
        If template_from_structure is False, then the initial structure is expected to be the equilibrium geometry.

        Residue offsets of the interactions are wrapped around the circular chains (see DNA.circularize),
        so the last residue of a circular chain is bonded, stacked, etc. to its first residue."""
        # Parse configuration file if not already done
        try:

            self.bond_definition
        except AttributeError:
            self.parseConfigurationFile()

        self.template_from_X3DNA = template_from_X3DNA
        DNAtype = self.DNAtype
        if DNAtype not in self.angle_definition['DNA'].unique():
            raise DNATypeError(self)

        # Rewrite index in case it is not ordered
        self.atoms.index = range(len(self.atoms))
        wrap = residue_wrapper(self.atoms, self.circular_chains)

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

        def key(c, r, name):
            """Index key of an atom, wrapping the residue number around circular chains"""
            return c, wrap(c, r), name

        # Compute B_curved geometry if needed
        if DNAtype == 'B_curved' and template_from_X3DNA:
            self.template_atoms = self.computeGeometry(sequence=self._templateSequence())

            def template_xyz(keys):
                nan = np.full(3, np.nan)
                return pandas.DataFrame([self.template_positions.get(k, nan) for k in keys], columns=['x', 'y', 'z'])
        else:
            self.template_atoms = self.atoms

            def template_xyz(keys):
                return self.atoms.loc[[index[key(*k)] for k in keys], ['x', 'y', 'z']].reset_index(drop=True)

        # Select ADNA bond definitions

        bond_types = self.bond_definition[self.bond_definition['DNA'] == DNAtype]
        angle_types = self.angle_definition[self.angle_definition['DNA'] == DNAtype]
        stacking_types = self.stacking_definition[self.stacking_definition['DNA'] == DNAtype]
        dihedral_types = self.dihedral_definition[self.dihedral_definition['DNA'] == DNAtype]
        # print(bond_types)
        # print(index)
        # Make a table with bonds
        # The template keys do not wrap the residue number, so the equilibrium geometry of the interactions that
        # close a circular chain is taken from the residues that pad the template (see DNA._templateSequence)
        data = []
        template_keys = []
        for i, ftype in bond_types.iterrows():
            # print(bond_type)
            ai = ftype['i']
            aj = ftype['j']
            s1 = ftype['s1']
            for c, r in cr_list:
                k1 = key(c, r, ai)
                k2 = key(c, r + s1, aj)
                if k1 in index and k2 in index:
                    data += [[i, index[k1], index[k2]]]
                    template_keys += [[(c, r, ai), (c, r + s1, aj)]]
        data = pandas.DataFrame(data, columns=['name', 'aai', 'aaj'])
        self.bonds = data.merge(bond_types, left_on='name', right_index=True)

        if DNAtype == 'B_curved':
            # Make default distances the same as the initial distance
            x1 = template_xyz([k[0] for k in template_keys])
            x2 = template_xyz([k[1] for k in template_keys])
            self.bonds['r0'] = np.linalg.norm(x1.values - x2.values, axis=1) / 10

        # Make a table with angles
        data = []
        template_keys = []
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
                k1 = key(c, r, ai)
                k2 = key(c, r + s1, aj)
                k3 = key(c, r + s2, ak)
                k4 = key(c, r + sb, 'S')
                if (
                    k1 in index
                    and k2 in index
                    and k3 in index
                    and k4 in index
                    and (b1 == '*' or base[index[k1]] == b1)
                    and (b2 == '*' or base[index[k4]] == b2)
                ):
                    data += [[i, index[k1], index[k2], index[k3], index[k4], sb]]
                    template_keys += [[(c, r, ai), (c, r + s1, aj), (c, r + s2, ak)]]
        data = pandas.DataFrame(data, columns=['name', 'aai', 'aaj', 'aak', 'aax', 'sB'])
        self.angles = data.merge(angle_types, left_on='name', right_index=True)

        if DNAtype == 'B_curved':
            # Make initial angles the default angles
            v1 = template_xyz([k[0] for k in template_keys])
            v2 = template_xyz([k[1] for k in template_keys])
            v3 = template_xyz([k[2] for k in template_keys])
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
                k1 = key(c, r, ai)
                k2 = key(c, r + s1, aj)
                k3 = key(c, r + s2, ak)
                if k1 in index and k2 in index and k3 in index:
                    data += [[i, index[k1], index[k2], index[k3]]]
        data = pandas.DataFrame(data, columns=['name', 'aai', 'aaj', 'aak'])
        self.stackings = data.merge(stacking_types, left_on='name', right_index=True)

        # Make a table with dihedrals
        data = []
        template_keys = []
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
                k1 = key(c, r, ai)
                k2 = key(c, r + s1, aj)
                k3 = key(c, r + s2, ak)
                k4 = key(c, r + s3, al)
                if k1 in index and k2 in index and k3 in index and k4 in index:
                    data += [[i, index[k1], index[k2], index[k3], index[k4]]]
                    template_keys += [[(c, r, ai), (c, r + s1, aj), (c, r + s2, ak), (c, r + s3, al)]]
        data = pandas.DataFrame(data, columns=['name', 'aai', 'aaj', 'aak', 'aal'])
        self.dihedrals = data.merge(dihedral_types, left_on='name', right_index=True)

        if DNAtype == 'B_curved':
            # Make initial dihedrals the default dihedrals

            a1 = template_xyz([k[0] for k in template_keys])
            a2 = template_xyz([k[1] for k in template_keys])
            a3 = template_xyz([k[2] for k in template_keys])
            a4 = template_xyz([k[3] for k in template_keys])

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

    def _templateSequence(self, padding=3):
        """Sequence used to build the B_curved template. Each circular chain is padded with the last
        `padding` residues before its start and the first `padding` residues after its end, numbered
        outside the chain, so the template also contains the geometry of the steps that close the chain."""
        sequence = self.getSequences()
        info = circular_chain_info(self.atoms, self.circular_chains)
        if not info:
            return sequence
        padded = {}
        for chain in sequence.index.get_level_values(0).unique():
            chain_sequence = sequence[chain]
            if chain in info:
                start, length = info[chain]
                bases = [chain_sequence[start + i] for i in range(length)]
                pad = min(padding, length)
                for i in range(pad):
                    padded[(chain, start - pad + i)] = bases[length - pad + i]
                for i in range(length):
                    padded[(chain, start + i)] = bases[i]
                for i in range(pad):
                    padded[(chain, start + length + i)] = bases[i]
            else:
                for r, b in chain_sequence.items():
                    padded[(chain, r)] = b
        return pandas.Series(padded)

    # Circular DNA
    def _dnaChains(self):
        """DNA chain identifiers in the order they appear in the structure"""
        return list(self.atoms.loc[self.atoms['resname'].isin(_dnaResidues), 'chainID'].unique())

    def circularize(self, chains=True, join=None, output_pdb='clean.pdb', compute_topology=None):
        """Bonds the last residue of each circular chain to its first residue (e.g. a plasmid).

        The bonds, angles, stackings, dihedrals and cross stackings that cross the end of a circular
        chain are added by computeTopology and the forces, and the sequence separation masks of the
        nonbonded forces count residues around the circle.

        Parameters
        ----------
        chains : True, str or list
            Chains to close. True closes every DNA chain (both strands of a plasmid).
        join : list of lists, optional
            Groups of chains joined end to end into a single chain before closing them. For example
            ``[['A', 'B']]`` bonds the 3' end of strand A to the 5' end of strand B, so the closed
            molecule is a single strand that runs along both sides of the duplex. The joined chain
            keeps the identifier of the first chain of the group.
        output_pdb : str or None
            Coarse grained pdb file to write the structure to (it is the topology read by System).
        compute_topology : bool, optional
            Recompute the bonded interactions. By default they are recomputed only if they had already
            been computed.

        The first residue of a linear chain has no phosphate. If the first residue of a circular chain
        has none, a phosphate is added between the two sugars that it connects."""
        renamed = {}
        for group in join or []:
            self._joinChains(group)
            renamed.update({c: group[0] for c in group})
        if chains is True:
            chains = self._dnaChains()
        elif not chains:
            chains = []
        elif isinstance(chains, str):
            chains = [chains]
        chains = list(dict.fromkeys(renamed.get(c, c) for c in chains))

        for chain, (start, length) in circular_chain_info(self.atoms, chains).items():
            if length < 4:
                raise ValueError(f'Circular chain {chain} has {length} residues, at least 4 are needed')
            self._addClosingPhosphate(chain, start, length)
        self.atoms.index = range(len(self.atoms))
        self.circular_chains = chains

        # Warn if the ends are not close enough to be bonded
        for chain, (start, length) in circular_chain_info(self.atoms, chains).items():
            in_chain = self.atoms['chainID'] == chain
            last_sugar = self.atoms[in_chain & (self.atoms['resSeq'] == start + length - 1) & (self.atoms['name'] == 'S')]
            first_phosphate = self.atoms[in_chain & (self.atoms['resSeq'] == start) & (self.atoms['name'] == 'P')]
            gap = np.linalg.norm(last_sugar[['x', 'y', 'z']].values[0] - first_phosphate[['x', 'y', 'z']].values[0]) / 10
            if gap > 1:
                print(f'Warning: the ends of the circular chain {chain} are {gap:.1f} nm apart. '
                      f'Bend the DNA (DNA.bendToCircle) or minimize the energy before running dynamics.')

        if compute_topology is None:
            compute_topology = hasattr(self, 'bonds')
        if compute_topology:
            self.computeTopology(template_from_X3DNA=getattr(self, 'template_from_X3DNA', True))
        if output_pdb:
            self.writePDB(output_pdb)
        return self

    def _joinChains(self, group):
        """Joins the chains in `group` end to end into a single chain named as the first one."""
        atoms = self.atoms
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
        self.atoms = atoms.sort_values('_rank', kind='stable').drop(columns='_rank').reset_index(drop=True)

    def _addClosingPhosphate(self, chain, start, length):
        """Adds the phosphate of the first residue of a circular chain if it is missing.

        The phosphate is copied from a residue in the middle of the chain, superposing the two sugars
        and bases it connects onto the last and first residues of the chain."""
        atoms = self.atoms
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

        # Copy the phosphate of a residue near the middle of the chain
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
        self.atoms = pandas.concat([atoms.iloc[:position_in_table], phosphate,
                                    atoms.iloc[position_in_table:]]).reset_index(drop=True)

    def _duplexBasePairs(self):
        """Base pairs of a two stranded DNA, ordered along the first strand.

        Returns the two chain identifiers, the index of the base pair of every atom of the duplex (-1
        for other atoms) and four (n, 3) arrays with the sugar and base of the first strand and the
        sugar and base of the second strand of each base pair."""
        chains = self._dnaChains()
        if len(chains) != 2:
            raise ValueError(f'A double stranded DNA with two chains is needed, the DNA has the chains {chains}')
        a, b = chains
        residues_a = sorted(self.atoms.loc[self.atoms['chainID'] == a, 'resSeq'].unique())
        residues_b = sorted(self.atoms.loc[self.atoms['chainID'] == b, 'resSeq'].unique())[::-1]
        if len(residues_a) != len(residues_b):
            raise ValueError(f'The chains {a} and {b} have a different number of residues')
        pair_index = {(a, r): i for i, r in enumerate(residues_a)}
        pair_index.update({(b, r): i for i, r in enumerate(residues_b)})
        base_pair = np.array([pair_index.get((c, r), -1) for c, r in zip(self.atoms['chainID'], self.atoms['resSeq'])])

        lookup = {}
        for c, r, n, x, y, z in zip(self.atoms['chainID'], self.atoms['resSeq'], self.atoms['name'],
                                    self.atoms['x'], self.atoms['y'], self.atoms['z']):
            lookup[(c, r, n if n in ['P', 'S'] else 'B')] = (x, y, z)
        beads = [np.array([lookup[(c, r, n)] for r in residues]) for c, residues in [(a, residues_a), (b, residues_b)]
                 for n in ['S', 'B']]
        return (a, b), base_pair, beads

    def _duplexEndsApart(self):
        """True if the two ends of a duplex are further apart than three base pair steps"""
        _, _, (sugar_a, base_a, sugar_b, base_b) = self._duplexBasePairs()
        centers = (sugar_a + base_a + sugar_b + base_b) / 4
        steps = np.linalg.norm(np.diff(centers, axis=0), axis=1)
        return np.linalg.norm(centers[-1] - centers[0]) > 3 * np.median(steps)

    def bendToCircle(self, delta_linking_number=0, superhelical_density=None, closed=False, first=0, last=None):
        """Bends a double stranded DNA into a planar circle so its ends can be joined.

        The circle lies on the xy plane, centered at the origin. The base pairs are placed rigidly along
        the circle keeping their rise and twist, and the twist is adjusted so the circle has an integer
        linking number Lk = Lk0 + delta_linking_number, where Lk0 is the number of helical turns of the
        relaxed duplex. A positive delta_linking_number overwinds the DNA (positive supercoiling) and a
        negative one underwinds it (negative supercoiling, like most plasmids in cells). In the planar
        circle all the excess goes to twist, and it converts partially to writhe (plectonemes) during
        the simulation.

        Parameters
        ----------
        delta_linking_number : int
            Turns added to the relaxed linking number.
        superhelical_density : float, optional
            Alternative to delta_linking_number, sigma = delta_Lk / Lk0 (e.g. -0.06).
        closed, first, last :
            See CircularBend. Used by fromSequence to close the circle exactly.

        Returns the CircularBend object, with the radius (Angstrom) and linking numbers."""
        (a, b), base_pair, (sugar_a, base_a, sugar_b, base_b) = self._duplexBasePairs()
        centers = (sugar_a + base_a + sugar_b + base_b) / 4
        bend = CircularBend(centers, sugar_b - sugar_a, closed=closed, delta_linking_number=delta_linking_number,
                            superhelical_density=superhelical_density, first=first, last=last)
        duplex = base_pair >= 0
        self.atoms.loc[duplex, ['x', 'y', 'z']] = bend.apply(self.atoms.loc[duplex, ['x', 'y', 'z']].values,
                                                            base_pair[duplex])
        self.bend = bend
        print(f'DNA bent into a circle of radius {bend.radius / 10:.2f} nm with linking number {bend.linking_number} '
              f'(relaxed {bend.relaxed_linking_number}, delta {bend.delta_linking_number:+d})')
        return bend

    @staticmethod
    def _trimCircle(atoms, first, length):
        """Keeps the base pairs first ... first + length - 1 (0 based, along the first strand) of a duplex,
        with the phosphates that join them to the previous base pairs, and renumbers both strands from 1."""
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

    def computeLinkingNumber(self, positions=None):
        """Linking number of a circular duplex, from the sugar beads of its two strands.

        `positions` is an optional (n_atoms, 3) array (e.g. a frame of a trajectory) in the same
        units as the structure; by default the current coordinates are used. Lk is a topological
        invariant: it only changes if a strand breaks, so it is a good check that a plasmid keeps its
        topology during a simulation."""
        atoms = self.atoms
        xyz = atoms[['x', 'y', 'z']].values if positions is None else np.asarray(positions)
        curves = []
        for chain in self._dnaChains():
            sugars = atoms[(atoms['chainID'] == chain) & (atoms['name'] == 'S')].sort_values('resSeq')
            curves += [xyz[sugars.index]]
        if len(curves) != 2:
            raise ValueError('The linking number is defined for two strands')
        return linking_number(*curves)

    def _makeCircular(self, circular, join, bend, delta_linking_number, superhelical_density):
        """Bends (if requested or if the ends are apart) and closes the DNA, used by the constructors"""
        if circular and bend == 'auto':
            bend = len(self._dnaChains()) == 2 and self._duplexEndsApart()
        if circular and bend:
            self.bendToCircle(delta_linking_number, superhelical_density)
        elif delta_linking_number or superhelical_density:
            print('Warning: the linking number is only set when the DNA is bent (bend=True)')
        self.circularize(chains=circular, join=join, output_pdb=None, compute_topology=False)

    def create_mutant(self, mutations):
        """Returns a new DNA identical to this one but with the given base identities.

        Single-topology: the mutant keeps the same beads and coordinates, only the base labels and
        the potential parameters change. `mutations` is an iterable of ``(chainID, resSeq, target)``
        (dicts with those keys or a DataFrame are also accepted) with ``target`` in {A, T, G, C}.
        Building the mutant is the general operation the alchemical transformation uses to obtain its
        target state, but it is also useful on its own (e.g. to score a sequence variant)."""
        element = {'A': 'N', 'T': 'S', 'G': 'C', 'C': 'O'}
        base_names = ['A', 'T', 'G', 'C']

        if isinstance(mutations, pandas.DataFrame):
            mutations = mutations[['chainID', 'resSeq', 'target']]
        else:
            mutations = pandas.DataFrame([(m['chainID'], m['resSeq'], m['target']) if isinstance(m, dict)
                                          else (m[0], m[1], m[2]) for m in mutations],
                                         columns=['chainID', 'resSeq', 'target'])
        for target in mutations['target']:
            if target not in base_names:
                raise ValueError(f"Mutation target '{target}' is not one of {base_names}")

        mutant = DNA(periodic=self.periodic)
        mutant.atoms = self.atoms.copy()
        mutant.DNAtype = self.DNAtype
        mutant.circular_chains = list(self.circular_chains)
        mutant.config = self.config
        for attr in ['particle_definition', 'bond_definition', 'angle_definition',
                     'dihedral_definition', 'stacking_definition', 'pair_definition', 'cross_definition']:
            setattr(mutant, attr, getattr(self, attr))

        for _, mut in mutations.iterrows():
            residue = (mutant.atoms['chainID'] == mut['chainID']) & (mutant.atoms['resSeq'] == mut['resSeq'])
            base_bead = residue & mutant.atoms['name'].isin(base_names)
            if int(base_bead.sum()) != 1:
                raise ValueError(f"Mutation {mut['chainID']}:{mut['resSeq']} matched {int(base_bead.sum())} "
                                 f"base beads (expected exactly 1)")
            # The residue name is shared by all three beads (read by the angle base-matching).
            mutant.atoms.loc[residue, 'resname'] = 'D' + mut['target']
            mutant.atoms.loc[base_bead, 'name'] = mut['target']
            mutant.atoms.loc[base_bead, 'type'] = mut['target']
            mutant.atoms.loc[base_bead, 'element'] = element[mut['target']]

        mutant.computeTopology(template_from_X3DNA=getattr(self, 'template_from_X3DNA', False))
        return mutant

    def writePDB(self, pdb_file='clean.pdb'):
        # Compute element fields
        element_ix = {'P': 'P', 'S': 'H', 'A': 'N', 'T': 'S', 'C': 'O', 'G': 'C'}  # Elements choosen to keep VMD colors
        self.atoms.loc[:, 'element'] = [element_ix[atomType] for atomType in self.atoms['name']]
        
        self.pdb_file = writePDB(self.atoms,pdb_file)
        return self.pdb_file

    @classmethod
    def fromCoarsePDB(cls, pdb_file, dna_type='B_curved', template_from_X3DNA=True, compute_topology=True,
                      circular=False, join=None, bend='auto', delta_linking_number=0, superhelical_density=None,
                      output_pdb='clean.pdb'):
        """Initializes a DNA object from a pdb file containing the Coarse Grained atoms

        Circular DNA (e.g. plasmids):
            circular: True closes every DNA chain, or a list of the chains to close.
            join: groups of chains joined end to end into a single chain, e.g. [['A', 'B']].
            bend: 'auto' bends a straight duplex into a circle only if its ends are apart, True always
                bends it and False never does.
            delta_linking_number, superhelical_density: supercoiling of the bent circle (see bendToCircle).
            output_pdb: the circular structure is written to this file, as it may gain phosphates or new
                coordinates. A linear structure keeps using the input file."""
        self = cls()

        self.atoms = parsePDB(pdb_file)
        self.atoms.loc[:, 'type'] = self.atoms['name']
        # Initialize the system from the pdb
        self.DNAtype = dna_type
        self.pdb_file = pdb_file
        if circular or join:
            self._makeCircular(circular, join, bend, delta_linking_number, superhelical_density)
        if compute_topology:
            self.parseConfigurationFile()
            self.computeTopology(template_from_X3DNA=template_from_X3DNA)
        if circular or join:
            self.writePDB(output_pdb)
        return self

    @classmethod
    def fromPDB(cls, pdb_file, dna_type='B_curved', template_from_X3DNA=True, output_pdb='clean.pdb', compute_topology=True,
                circular=False, join=None, bend='auto', delta_linking_number=0, superhelical_density=None):
        """Creates a DNA object from a complete(atomistic) pdb file

        Circular DNA (e.g. plasmids):
            circular: True closes every DNA chain, or a list of the chains to close.
            join: groups of chains joined end to end into a single chain, e.g. [['A', 'B']].
            bend: 'auto' bends a straight duplex into a circle only if its ends are apart, True always
                bends it and False never does.
            delta_linking_number, superhelical_density: supercoiling of the bent circle (see bendToCircle)."""
        self = cls()
        pdb = fixPDB(pdb_file)
        pdb_table = pdb2table(pdb)

        if circular is True or join:
            keep_phosphate = True
        else:
            keep_phosphate = [circular] if isinstance(circular, str) else (circular or [])
        self.atoms = self.CoarseGrain(pdb_table, circular_chains=keep_phosphate)
        self.DNAtype = dna_type
        if circular or join:
            self._makeCircular(circular, join, bend, delta_linking_number, superhelical_density)
        if compute_topology:
            self.parseConfigurationFile()
            self.computeTopology(template_from_X3DNA=template_from_X3DNA)
        if output_pdb:
            self.writePDB(output_pdb)
        #self.atomistic_model=temp
        return self

    @staticmethod
    def CoarseGrain(pdb_table, rna=False, circular_chains=()):
        """ Selects DNA atoms from a pdb table and returns a table containing only the coarse-grained atoms for 3SPN2

        The phosphate of the first residue of each chain is removed, except for the `circular_chains`
        (True for all of them) whose first residue has a phosphorus atom."""
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
        has_phosphorus = set(zip(temp.loc[temp['name'] == 'P', 'chainID'], temp.loc[temp['name'] == 'P', 'resSeq']))
        drop_list = []
        for chain in Coarse.chainID.unique():
            sel = Coarse[Coarse.chainID == chain]
            first = sel.resSeq.min()
            if (circular_chains is True or chain in circular_chains) and (chain, first) in has_phosphorus:
                continue
            drop_list += list(sel[(sel.resSeq == first) & sel['name'].isin(['P'])].index)
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
    def fromSequence(cls, sequence, dna_type='B_curved', output_pdb='clean.pdb', compute_topology=True,atomistic_template='atomistic_dna',
                     circular=False, join=None, delta_linking_number=0, superhelical_density=None):
        """ Initializes a DNA object from a DNA sequence

        Circular DNA (e.g. plasmids):
            circular: if True the duplex is bent into a planar circle and both strands are closed.
            join: groups of chains joined end to end into a single chain, e.g. [['A', 'B']].
            delta_linking_number: turns added to (positive) or removed from (negative) the relaxed
                linking number of the circle, to build supercoiled DNA.
            superhelical_density: alternative to delta_linking_number, sigma = delta_Lk / Lk0."""
        if circular:
            return cls._circularFromSequence(sequence, dna_type=dna_type, output_pdb=output_pdb,
                                             compute_topology=compute_topology, atomistic_template=atomistic_template,
                                             circular=circular, join=join, delta_linking_number=delta_linking_number,
                                             superhelical_density=superhelical_density)
        if join:
            self = cls.fromSequence(sequence, dna_type=dna_type, output_pdb=None, compute_topology=False,
                                    atomistic_template=atomistic_template)
            self.circularize(chains=False, join=join, output_pdb=None, compute_topology=False)
            if compute_topology:
                self.parseConfigurationFile()
                self.computeTopology()
            self.writePDB(output_pdb)
            return self
        self = cls()
        self.parseConfigurationFile()
        sequence = pandas.Series([a for a in sequence], index=[('A', i) for i in range(len(sequence))])
        # Make a possible structure
        self.computeGeometry(sequence, atomistic_template=atomistic_template)
        self.DNAtype=dna_type

        # Make a clean pdb file
        self = self.fromPDB(atomistic_template, dna_type=dna_type, output_pdb=output_pdb,
                            compute_topology=compute_topology)
        return self

    @classmethod
    def _circularFromSequence(cls, sequence, dna_type, output_pdb, compute_topology, atomistic_template,
                              circular, join, delta_linking_number, superhelical_density, padding=6):
        """Builds a circular duplex from a sequence.

        A straight duplex is built with the sequence, one extra base pair that copies the first one, and
        `padding` base pairs at each end that continue the circular sequence, so the first and last base
        pairs keep the same neighbors they have in the circle. It is bent so the extra base pair falls on
        the first one, which closes the circle exactly, and the extra base pairs are removed."""
        n = len(sequence)
        if n < 4:
            raise ValueError('A circular DNA needs at least 4 base pairs')
        extended = [sequence[i % n] for i in range(-padding, n + padding + 1)]
        extended = pandas.Series(extended, index=[('A', i) for i in range(len(extended))])
        builder = cls()
        builder.parseConfigurationFile()
        atomistic_file = atomistic_template or 'atomistic_dna'
        builder.computeGeometry(extended, atomistic_template=atomistic_file)

        self = cls.fromPDB(atomistic_file, dna_type=dna_type, output_pdb=None, compute_topology=False)
        bend = self.bendToCircle(delta_linking_number, superhelical_density, closed=True,
                                 first=padding, last=padding + n)
        self.atoms = self._trimCircle(self.atoms, padding, n)

        if atomistic_template:
            # Bent atomistic structure
            atomistic = builder.template.reset_index(drop=True)
            total = len(extended)
            base_pair = np.where(atomistic['chainID'] == 'A', atomistic['resSeq'] - 1, total - atomistic['resSeq'])
            atomistic[['x', 'y', 'z']] = bend.apply(atomistic[['x', 'y', 'z']].values, base_pair)
            atomistic = self._trimCircle(atomistic, padding, n)
            atomistic['serial'] = atomistic.index + 1
            writePDB(atomistic, pdb_file=atomistic_template)

        self.circularize(chains=circular, join=join, output_pdb=None, compute_topology=False)
        if compute_topology:
            self.parseConfigurationFile()
            self.computeTopology()
        self.writePDB(output_pdb)
        return self

    @classmethod
    def fromXYZ(cls, xyz_file, dnatype='B_curved', template_from_X3DNA=True, output_pdb='clean.pdb', compute_topology=True):
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
            self.computeTopology(template_from_X3DNA=template_from_X3DNA)
        self.writePDB(output_pdb)
        return self


class System(openmm.System):
    """ Wrapper of openmm system class, adds some openmm simulation attributes"""

    def __init__(self, dna, forcefieldFiles=[__location__/'3SPN2.xml'], periodicBox=None):
        self.dna = dna
        self.top = openmm.app.PDBFile(dna.pdb_file).getTopology()
        if self.top.getUnitCellDimensions() is None:
            x = dna.atoms[['x', 'y', 'z']]
            d = np.round((x.max() - x.min()).values * 2 + 5, -1)
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
            force.addForce(self)
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
