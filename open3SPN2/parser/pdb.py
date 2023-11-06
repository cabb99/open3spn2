import pandas
import pdbfixer
import openmm.unit as unit

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
        return dict(
            recname=str(line[:6]).strip(),
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
            occupancy=1.0 if line[54:60].strip() == '' else float(line[54:60]),
            tempFactor=1.0
            if line[60:66].strip() == ''
            else float(line[60:66]),
            element=str(line[76:78]),
            charge=str(line[78:80]),
        )

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
        if key[1] in [0, len(list(chain_tmp.residues()))]:
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
        pos = pos.value_in_unit(unit.angstrom)
        data += [dict(zip(cols, ['ATOM', int(atom.id), atom.name, '',
                                 residue.name, chain.id, int(residue.id), '',
                                 pos[0], pos[1], pos[2], 0, 0,
                                 atom.element.symbol, '']))]
    atom_list = pandas.DataFrame(data)
    atom_list = atom_list[cols]
    atom_list.index = atom_list['serial']
    return atom_list

def writePDB(atoms, pdb_file='clean.pdb'):
    """ Writes a minimal version of the pdb file needed for openmm """
    # Compute chain field
    if type(atoms['chainID'].iloc[0]) is not str:
        chain_ix = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
        atoms['chainID'] = [chain_ix[i - 1] for i in atoms['chainID']]

    #Complete missing fields
    for field, value in {'occupancy':0.0, 'beta':0.0, 'segment':'', 'charge':'', 'insertion':''}.items():
        if field not in atoms:
            atoms[field]=value

    # Write pdb file
    with open(pdb_file, 'w+') as pdb:
        for i, atom in atoms.iterrows():
            pdb_line = f'ATOM  {i + 1:>5} {atom["name"]:^4} {atom.resname:<3} {atom.chainID:1}{atom.resSeq:>4}{atom.insertion:1}   {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}{atom.occupancy:>6.2f}{atom.beta:>6.2f}      {atom.segment:4}{atom.element:2}{atom.charge:2}'
            assert len(pdb_line) == 80, 'An item in the atom table is longer than expected'
            pdb.write(pdb_line + '\n')
    return pdb_file