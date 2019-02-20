import simtk.openmm.app
import simtk.openmm
from simtk.unit import *
import configparser
import numpy as np
import itertools

"""
This tool has been constructed to simulate DNA using the 3SPN2 forcefield 
in openmm

Author: Carlos Bueno
"""

_ef = 1 * kilocalorie / kilojoule  # energy scaling factor
_df = 1 * angstrom / nanometer  # distance scaling factor
_af = 1 * degree / radian  # angle scaling factor


def parseConfigTable(config_section):
    '''Parses a section of the configuration file as a table'''

    def readData(config_section, a):
        '''Filters comments and returns values as a list'''
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


## Errors
class BaseError(Exception):
    pass


class DNATypeError(BaseError):
    def __init__(self, dna):
        self.dna = dna
        self.message = f'DNA type {dna.DNAType} not defined in the configuration file\n'
        defined_types = dna.angle_definition['DNA'].unique()
        self.message += f'Only the types {str(defined_types)} were defined'
        print(self.message)


class DNA(object):
    def __init__(self,periodic=True):
        '''Initializes an empty DNA object'''
        self.periodic = periodic
        pass

    def __repr__(self):
        return f'DNA object'
        # print the sequence and the identity of the DNA object

    def parseConfigurationFile(self, configuration_file='3SPN2.conf'):
        self.configuration_file = configuration_file
        config = configparser.ConfigParser()
        config.read(configuration_file)
        self.particle_definition = parseConfigTable(config['Particles'])
        self.bond_definition = parseConfigTable(config['Bonds'])
        self.angle_definition = parseConfigTable(config['Angles'])
        self.dihedral_definition = parseConfigTable(config['Dihedrals'])
        self.pair_definition = parseConfigTable(config['Base Pairs'])
        self.cross_definition = parseConfigTable(config['Cross stackings'])

    def computeTopology(self, DNAType='A'):
        # Parse configuration file if not already done
        try:
            self.bond_definition
        except AttributeError:
            self.parseConfigurationFile()

        self.DNAType = DNAType
        if DNAType not in self.angle_definition['DNA'].unique():
            raise DNATypeError(self)

        # Make an index to build the topology
        index = {}
        cr_list = set()  # Chain residue list
        for i, atom in self.atoms.iterrows():
            index.update({(atom['chain'], atom['residue'], atom['type']): i})
            cr_list.update([(atom['chain'], atom['residue'])])
        cr_list = list(cr_list)
        cr_list.sort()
        max_chain = self.atoms['chain'].max()
        max_residue = self.atoms['residue'].max()
        assert len(index) == len(self.atoms), 'Atom index was repeated'

        # Select ADNA bond definitions

        angle_types = self.angle_definition[self.angle_definition['DNA'] == DNAType]
        if DNAType == 'B_curved':
            DNAType = 'B'
        bond_types = self.bond_definition[self.bond_definition['DNA'] == DNAType]
        dihedral_types = self.dihedral_definition[self.dihedral_definition['DNA'] == DNAType]

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
        data = pandas.DataFrame(data, columns=['type', 'aai', 'aaj'])
        self.bonds = data.merge(bond_types, left_on='type', right_index=True)

        # Make a table with angles
        data = []
        for i, ftype in angle_types.iterrows():
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
        data = pandas.DataFrame(data, columns=['type', 'aai', 'aaj', 'aak'])
        self.angles = data.merge(angle_types, left_on='type', right_index=True)

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
        data = pandas.DataFrame(data, columns=['type', 'aai', 'aaj', 'aak', 'aal'])
        self.dihedrals = data.merge(dihedral_types, left_on='type', right_index=True)

    def writePDB(self, pdb_file='clean.pdb'):
        ''' Writes a minimal version of the pdb file needed for openmm '''
        # Compute chain field
        chain_ix = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
        self.atoms['chainid'] = [chain_ix[i - 1] for i in self.atoms['chain']]
        # compute resid and resname fields
        res_ix = {}
        min_res = self.atoms.groupby('chain')['residue'].min()
        max_res = self.atoms.groupby('chain')['residue'].max()
        for i, res in self.atoms[~self.atoms['type'].isin(['S', 'P'])].iterrows():
            resname = 'D' + res['type']
            if res['residue'] == min_res[res['chain']]:
                resname += 'i'
            if res['residue'] == max_res[res['chain']]:
                resname += 'f'
            res_ix.update({(res['chain'], res['residue']): resname})
        self.atoms['resname'] = [res_ix[(r['chain'], r['residue'])] for i, r in self.atoms.iterrows()]

        # Compute element fields
        element_ix = {'P': 'P', 'S': 'H', 'A': 'N', 'T': 'S', 'C': 'O', 'G': 'C'}  # Elements choosen to keep VMD colors
        self.atoms['element'] = [element_ix[atomType] for atomType in self.atoms['type']]

        # Write pdb file
        with open(pdb_file, 'w+') as pdb:
            for i, atom in self.atoms.iterrows():
                pdb_line = f'ATOM  {i + 1:>5} {atom.type:^4} {atom.resname:<3}  {atom.chainid}{atom.residue:>3}    {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}' + ' ' * 22 + f'{atom.element:2}' + ' ' * 2
                assert len(pdb_line) == 80, 'An item in the atom table is longer than expected'
                pdb.write(pdb_line + '\n')
        self.pdb_file = pdb_file

    @classmethod
    def fromPDB(cls, pdb_file):
        '''Initializes a DNA object from a pdb file'''
        # Parse the pdb

        # Make a clean pdb file

        # Initialize the system from the pdb

        pass

    @classmethod
    def fromGRO(cls, gro_file):
        '''Initializes a DNA object from a gromacs input file'''
        # Parse the gromacs file

        # Make a clean pdb file

        # Initialize the system from the pdb

        pass

    @classmethod
    def fromSequence(cls, sequence, center=[0, 0, 0]):
        '''Initializes a DNA object from a DNA sequence'''
        # Make a possible structure

        # Make a clean pdb file

        # Initialize the system from the pdb

        pass

    @classmethod
    def fromXYZ(cls, xyz_file):
        '''Initializes DNA object from xyz file (as seen on the examples)'''
        # Parse the file
        self = cls()
        self.atoms = pandas.read_csv('examples/adna/in00_conf.xyz', delim_whitespace=True, skiprows=2,
                                     names=['type', 'x', 'y', 'z'])

        # Compute residues and chains
        residue = 0
        residues = []
        chain = 0
        chains = []
        sugar = True
        for t in self.atoms['type']:
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
        self.atoms['residue'] = residues
        self.atoms['chain'] = chains
        return self


class System(simtk.openmm.System):
    ''' Wrapper of openmm system class, adds some openmm simulation attributes'''

    def __init__(self, mol, forcefieldFiles=['3SPN2.xml'], periodicBox=None):
        self.top = simtk.openmm.app.PDBFile(mol.pdb_file)
        self.coord = simtk.openmm.app.PDBFile(mol.pdb_file)
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
        ''' Removes all the forces from the system '''
        for i, f in enumerate(self.getForces()):
            if keepCMMotionRemover and i == 0 and f.__class__ == simtk.openmm.CMMotionRemover:
                continue
            else:
                self.removeForce(0)
        if keepCMMotionRemover == False:
            assert len(self.getForces()) == 0, 'Not all the forces were removed'
        else:
            assert len(self.getForces()) <= 1, 'Not all the forces were removed'

    def initializeMD(self, temperature=300 * kelvin):
        self.integrator = simtk.openmm.LangevinIntegrator(temperature, 1E-4/picosecond, 2*femtoseconds)
        platform = simtk.openmm.Platform.getPlatformByName('Reference')
        self.simulation = simtk.openmm.app.Simulation(self.top.topology, self._wrapped_system, self.integrator,
                                                      platform)
        self.simulation.context.setPositions(self.coord.positions)
        return self.simulation

    def setPositions(self, coords=None):
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

    def getPotentialEnergy(self, coords=None, unit=kilojoule_per_mole):
        # Initialize trial MD if not setup
        try:
            self.simulation
        except AttributeError:
            self.initializeMD()
        self.setPositions(coords)
        state = self.simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit)
        return energy

    def recomputeEnergy(self, trajectory):
        traj = parse_xyz(trajectory)
        self.initializeMD()
        energies = []
        for time, snapshot in traj.groupby('timestep'):
            energy = self.getPotentialEnergy(np.array(snapshot[['x', 'y', 'z']]) * _df)
            energies += [energy]
        return np.array(energies)


class Force(object):

    def __init__(self, dna):
        periodic = dna.periodic
        self.reset(periodic)
        self.defineInteraction(dna)

    def computeEnergy(self, system, trajectory):
        # Parse trajectory
        traj = ff3SPN2.parse_xyz('examples/adna/traj.xyz')

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
        traj = ff3SPN2.parse_xyz('examples/adna/traj.xyz')
        # for each item of the table:
        # clear all forces on the system
        # system.clearForces()
        # setup the force
        # self.setUpInteraction()
        # add the force item

        # compute the energy for every frame

        # return a table with the energy

    def setUpInteraction(self,periodic):
        ''' Creates a new force instance '''
        pass

    def defineInteraction(self, dna):
        ''' Adds the parameter for the force '''
        pass


class Bond3SPN2(Force, simtk.openmm.CustomBondForce):

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.force, attr)

    def getParameterNames():
        self.perInteractionParameters = []
        self.GlobalParameters = []
        for i in range(self.force.getNumPerBondParameters()):
            self.perInteractionParameters += [self.force.getPerBondParameterName(i)]
        for i in range(self.force.getNumGlobalParameters()):
            self.GlobalParameters += [self.force.getGlobalParameterName(i)]
        self.perInteractionParameters, self.GlobalParameters

    def reset(self, periodic=False):
        bondForce = simtk.openmm.CustomBondForce("Kb2*(r-r0)^2+Kb3*(r-r0)^3+Kb4*(r-r0)^4")
        bondForce.addPerBondParameter('r0')
        bondForce.addPerBondParameter('Kb2')
        bondForce.addPerBondParameter('Kb3')
        bondForce.addPerBondParameter('Kb4')
        bondForce.setUsesPeriodicBoundaryConditions(periodic)
        self.force = bondForce

    def defineInteraction(self, dna):
        for i, b in dna.bonds.iterrows():
            # Units converted from
            parameters = [b['r0'] * _df,
                          b['Kb2'] / _df ** 2 * _ef,
                          b['Kb3'] / _df ** 3 * _ef,
                          b['Kb4'] / _df ** 4 * _ef]
            self.force.addBond(int(b['aai']), int(b['aaj']), parameters)


class Angle3SPN2(Force, simtk.openmm.HarmonicAngleForce):

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.force, attr)

    def reset(self, periodic=False):
        angleForce = simtk.openmm.HarmonicAngleForce()
        angleForce.setUsesPeriodicBoundaryConditions(periodic)
        self.force = angleForce

    def defineInteraction(self, dna):
        for i, a in dna.angles.iterrows():
            if a['type_y'] == 'harmonic':
                parameters = [a['t0'] * _af,
                              a['epsilon'] * _ef * 2]
                self.force.addAngle(int(a['aai']), int(a['aaj']), int(a['aak']), *parameters)


class Stacking3SPN2(Force, simtk.openmm.CustomCompoundBondForce):
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.force, attr)

    def reset(self,periodic=False):
        stackingForce = simtk.openmm.CustomCompoundBondForce(3, """rep+f2*attr;
                                                                rep=epsilon*(1-exp(-alpha*(dr)))^2*step(-dr);
                                                                attr=epsilon*(1-exp(-alpha*(dr)))^2*step(dr)-epsilon;
                                                                dr=distance(p2,p3)-sigma;
                                                                f2=max(f*pair2,pair1);
                                                                pair1=step(dt+pi/2)*step(pi/2-dt);
                                                                pair2=step(dt+pi)*step(pi-dt);
                                                                f=1-cos(dt)^2;
                                                                dt=rng*(angle(p1,p2,p3)-t0);""")
        stackingForce.setUsesPeriodicBoundaryConditions(periodic)
        stackingForce.addPerBondParameter('epsilon')
        stackingForce.addPerBondParameter('sigma')
        stackingForce.addPerBondParameter('t0')
        stackingForce.addPerBondParameter('alpha')
        stackingForce.addPerBondParameter('rng')
        stackingForce.addGlobalParameter('pi', np.pi)
        self.force=stackingForce

    def defineInteraction(self, dna):
        for i, a in dna.angles.iterrows():
            if a['type_y'] == 'stacking/3spn2':
                parameters = [a['epsilon'] * _ef,
                              a['sigma'] * _df,
                              a['t0'] * _af,
                              a['alpha'] / _df,
                              a['rng']]
                self.force.addBond([a['aai'], a['aaj'], a['aak']], parameters)


class Dihedral3SPN2(Force, simtk.openmm.CustomTorsionForce):
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.force, attr)

    def reset(self, periodic=False):
        dihedralForce = simtk.openmm.CustomTorsionForce("-K_dihedral*exp(-acos(cos(theta-t0))^2/2/sigma^2)")
        # dihedralForce=simtk.openmm.CustomTorsionForce("theta/60.")
        dihedralForce.setUsesPeriodicBoundaryConditions(periodic)
        dihedralForce.addPerTorsionParameter('K_dihedral')
        dihedralForce.addPerTorsionParameter('t0')
        dihedralForce.addPerTorsionParameter('sigma')
        dihedralForce.addGlobalParameter('pi', np.pi)
        self.force=dihedralForce

    def defineInteraction(self, dna):
        for i, a in dna.dihedrals.iterrows():
            parameters = [a['K_dihedral'] * _ef,
                          (180 + a['t0']) * _af,
                          a['sigma']]
            particles = [a['aai'], a['aaj'], a['aak'], a['aal']]
            self.force.addTorsion(*particles, parameters)

class BasePair3SPN2(Force,simtk.openmm.CustomHbondForce):
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.force, attr)

    def reset(self, periodic=False):
        pairForce = simtk.openmm.CustomHbondForce(''' temp;temp=rep+1/2*(1+cos(dphi))*fdt1*fdt2*attr;
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

        if periodic:
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
        self.force=pairForce

    def defineInteraction(self, dna):
        DD, AA = [], []
        for i, pair in dna.pair_definition.iterrows():
            D1 = list(dna.atoms[dna.atoms['type'] == pair['Base1']].index)
            A1 = list(dna.atoms[dna.atoms['type'] == pair['Base2']].index)
            # D1=[1]
            # A1=[189]
            D2 = [j - 1 for j in D1]
            A2 = [j - 1 for j in A1]
            # parameters=[0.1,0.2,0.3,0.4,0.5,0.6,0.7]
            parameters = [pair.torsion * _af,
                          pair.sigma * _df,
                          pair.t1 * _af,
                          pair.t2 * _af,
                          pair.rang,
                          pair.epsilon * _ef,
                          pair.alpha / _df]
            # parameters=[50.17*af,pair.sigma*df,pair.t1*af,pair.t2*af,pair.rang,pair.epsilon*ef,pair.alpha/df]
            # not sure if including the same particle twice will break things
            d_num = self.force.getNumDonors()
            for d1, d2 in zip(D1, D2):
                self.force.addDonor(d1, d2, d2, parameters)
                # break
            a_num = self.force.getNumAcceptors()
            for a1, a2 in zip(A1, A2):
                self.force.addAcceptor(a1, a2, a2)
                # break

            DD += [D1]
            AA += [A1]
            # break

        # Map donors and aceptors
        d_num = 0
        a_num = 0
        atoms = dna.atoms.copy()
        for D1, A1 in zip(DD, AA):
            donors = atoms.iloc[D1]
            aceptors = atoms.iloc[A1]
            atoms.loc[D1, 'donor_id'] = [i + d_num for i in range(len(donors))]
            atoms.loc[A1, 'aceptor_id'] = [i + a_num for i in range(len(aceptors))]
            d_num += len(D1)
            a_num += len(A1)

        # Exclude interactions
        interaction = np.zeros([32, 32]) + 1
        for i, D1 in enumerate(DD):
            for j, A1 in enumerate(AA):
                print('Hi', ['A', 'G'][i], ['T', 'C'][j])
                if i == j:  # Exclude same chain interactions
                    for k, chain in atoms.loc[D1 + A1].groupby('chain'):
                        di = list(chain.donor_id.dropna().astype(int))
                        ai = list(chain.aceptor_id.dropna().astype(int))
                        for d, a in itertools.product(di, ai):
                            if abs(a - d) <= 3:
                                self.force.addExclusion(d, a)
                                # print(d,a)
                                interaction[d, a] = 0
                else:  # Exclude non-Watson-Crick interactions
                    di = list(atoms.loc[D1, 'donor_id'].astype(int))
                    ai = list(atoms.loc[A1, 'aceptor_id'].astype(int))
                    for d, a in itertools.product(di, ai):
                        self.force.addExclusion(d, a)
                        # print (d,a)
                        interaction[d, a] = 0


class CrossStacking3SPN2(Force):
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.force, attr)

    def reset(self, periodic=False):
        def crossStackingForce(parametersOnDonor=False):
            crossForce = simtk.openmm.CustomHbondForce('''energy;temp=dtCS/rng_CS;
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
            if periodic:
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
            return crossForce

        crossStackingForces = {}
        for base in ['A', 'T', 'G', 'C']:
            crossStackingForces.update({base: (crossStackingForce(), crossStackingForce())})
        self.crossStackingForces=crossStackingForces

    def defineInteraction(self, dna):
        # Donators
        bases = np.array(dna.atoms[dna.atoms['type'].isin(['A', 'T', 'G', 'C'])].index)
        D1 = dna.atoms.reindex(bases)
        D2 = dna.atoms.reindex(bases - 1)
        D3 = dna.atoms.reindex(bases + 3)
        A1 = dna.atoms.reindex(bases)
        A2 = dna.atoms.reindex(bases - 1)
        A3 = dna.atoms.reindex(bases - 3)

        # Select only bases where the other atoms exist
        D2.index = D1.index
        D3.index = D1.index
        temp = pandas.concat([D1, D2, D3], axis=1, keys=['D1', 'D2', 'D3'])
        sel = np.array(temp[temp['D3', 'type'].isin(['A', 'T', 'G', 'C']) & temp['D2', 'type'].isin(['S'])].index)
        D1 = dna.atoms.reindex(sel)
        D2 = dna.atoms.reindex(sel - 1)
        D3 = dna.atoms.reindex(sel + 3)

        # Aceptors
        A2.index = A1.index
        A3.index = A1.index
        temp = pandas.concat([A1, A2, A3], axis=1, keys=['A1', 'A2', 'A3'])
        sel = np.array(temp[temp['A3', 'type'].isin(['A', 'T', 'G', 'C']) & temp['A2', 'type'].isin(['S'])].index)
        A1 = dna.atoms.reindex(sel)
        A2 = dna.atoms.reindex(sel - 1)
        A3 = dna.atoms.reindex(sel - 3)

        # Parameters
        cross_definition = dna.cross_definition.copy()
        i = [a for a in zip(cross_definition['Base_d1'], cross_definition['Base_a1'], cross_definition['Base_a3'])]
        cross_definition.index = i

        donors = {i: [] for i in ['A', 'T', 'G', 'C']}
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        for donator, donator2, d1, d2, d3 in zip(D1.itertuples(), D3.itertuples(), D1.index, D2.index, D3.index):
            # if d1!=4:
            #    continue
            d1t = donator.type
            d3t = donator2.type
            c1, c2 = self.crossStackingForces[d1t]
            a1t = complement[d1t]
            #print(d1, d2, d3)
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
            a1t = aceptor.type
            a3t = aceptor2.type
            c1, c2 = self.crossStackingForces[complement[a1t]]
            d1t = complement[a1t]
            param = cross_definition.loc[[(d1t, a1t, a3t)]].squeeze()
            # print(param)
            #print(a1, a2, a3)
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
            aceptors[complement[a1t]] += [a1]

        # Exclusions
        for base in ['A', 'T', 'G', 'C']:
            c1, c2 = self.crossStackingForces[base]
            for ii, i in enumerate(donors[base]):
                for jj, j in enumerate(aceptors[base]):
                    if (dna.atoms.at[i, 'chain'] == dna.atoms.at[j, 'chain'] and abs(i - j) <= 9) or i > j:
                        c1.addExclusion(ii, jj)
                        c2.addExclusion(jj, ii)

    def addForce(self,system):
        for c1, c2 in self.crossStackingForces.values():
            system.addForce(c1)
            system.addForce(c2)


class Exclusion3SPN2(Force, simtk.openmm.CustomNonbondedForce):
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.force, attr)

    def reset(self, periodic=False):
        exclusionForce = simtk.openmm.CustomNonbondedForce("""energy;
                                                            energy=(epsilon*((sigma/r)^12-2*(sigma/r)^6)+epsilon)*step(sigma-r);
                                                            sigma=0.5*(sigma1+sigma2); 
                                                            epsilon=sqrt(epsilon1*epsilon2)""")
        exclusionForce.addPerParticleParameter('epsilon')
        exclusionForce.addPerParticleParameter('sigma')
        exclusionForce.setCutoffDistance(1.8)
        if periodic:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffPeriodic)
        else:
            exclusionForce.setNonbondedMethod(exclusionForce.CutoffNonPeriodic)
        self.force=exclusionForce

    def defineInteraction(self, dna):
        # addParticles
        particle_definition = dna.particle_definition[dna.particle_definition['DNA'] == 'A']
        particle_definition.index = particle_definition.name
        self.force.setCutoffDistance(particle_definition.radius.max() * _df)
        for i, atom in dna.atoms.iterrows():
            param = particle_definition.loc[atom.type]
            parameters = [param.epsilon * _ef,
                          param.radius * _df]
            #print(i, parameters)
            self.force.addParticle(parameters)

        # addExclusions
        for i, atom_a in dna.atoms.iterrows():
            for j, atom_b in dna.atoms.iterrows():
                if j > i:
                    continue
                # Neighboring residues
                if atom_a['chain'] == atom_b['chain'] and (atom_a.residue - atom_b.residue <= 1):
                    self.force.addExclusion(i, j)


class Electrostatics3SPN2(Force, simtk.openmm.CustomNonbondedForce):
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.force, attr)

    def reset(self, periodic=False):
        T = 300 * kelvin
        C = 100 * millimolar
        e = 249.4 - 0.788 * (T / kelvin) + 7.2E-4 * (T / kelvin) ** 2
        a = 1 - 0.2551 * (C / molar) + 5.151E-2 * (C / molar) ** 2 - 6.889E-3 * (C / molar) ** 3
        print(e, a)
        dielectric = e * a
        # Debye length
        kb = simtk.unit.BOLTZMANN_CONSTANT_kB  # Bolztmann constant
        Na = simtk.unit.AVOGADRO_CONSTANT_NA  # Avogadro number
        ec = 1.60217653E-19 * coulomb  # proton charge
        pv = 8.8541878176E-12 * farad / meter  # dielectric permittivity of vacuum

        ldby = sqrt(dielectric * pv * kb * T / (2.0 * Na * ec ** 2 * C))
        ldby = ldby.in_units_of(nanometer)
        denominator = 4 * np.pi * pv * dielectric / Na * ec ** 2
        print(ldby, denominator)

        electrostaticForce = simtk.openmm.CustomNonbondedForce("""energy;
                                                                energy=q1*q2*exp(-r/dh_length)/denominator/r;""")
        electrostaticForce.addPerParticleParameter('q')
        electrostaticForce.addGlobalParameter('dh_length', ldby)
        electrostaticForce.addGlobalParameter('denominator', denominator)

        electrostaticForce.setCutoffDistance(5)
        if periodic:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffPeriodic)
        else:
            electrostaticForce.setNonbondedMethod(electrostaticForce.CutoffNonPeriodic)
        self.force=electrostaticForce

    def defineInteraction(self, dna):
        # addParticles
        particle_definition = dna.particle_definition[dna.particle_definition['DNA'] == 'A']
        particle_definition.index = particle_definition.name
        for i, atom in dna.atoms.iterrows():
            param = particle_definition.loc[atom.type]
            parameters = [param.charge]
            # print (i,parameters)
            self.force.addParticle(parameters)

        # add neighbor exclusion
        for i, a1 in dna.atoms.iterrows():
            for j, a2 in dna.atoms.iterrows():
                if j < i:
                    continue
                if (a1.chain == a2.chain) and (abs(a1.residue - a2.residue) <= 1):
                    self.force.addExclusion(i, j)
                    # print (i,j,a1.chain,a2.chain,a1.residue,a2.residue,(a1.residue-a2.residue))


# Unit testing
def test_DNA_from_pdb():
    pass


def test_DNA_from_gro():
    pass


def test_DNA_from_seq():
    pass


def test_DNA_from_xyz():
    mol = DNA.fromXYZ('examples/adna/in00_conf.xyz')
    assert mol.atoms.at[8, 'type'] == 'P'
    assert round(mol.atoms.at[188, 'y'], 6) == -8.779343


# Functional testing
import pandas


def parse_xyz(filename=''):
    import pandas
    columns = ['N', 'timestep', 'id', 'type', 'x', 'y', 'z']
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
                atom['type'] = int(s[0])
                atom['x'], atom['y'], atom['z'] = [float(a) for a in s[1:4]]
                data += [atom.copy()]
                atom['id'] += 1
    xyz_data = pandas.concat(data, axis=1).T
    for i in ['N', 'timestep', 'id', 'type']:
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
    xyz_data = parse_xyz('examples/adna/traj.xyz')
    assert xyz_data.at[1, 'type'] == 7
    assert xyz_data.at[1, 'x'] == 4.34621


def test_parse_log():
    log_data = parse_log('examples/adna/sim.log')
    assert log_data.at[1, 'Step'] == 2000
    assert log_data.at[1, 'eexcl'] == 0.45734636


class TestForces:

    @classmethod
    def setup(cls, dna_type='A'):
        cls.dna = DNA.fromXYZ('examples/adna/in00_conf.xyz')

        cls.dna.parseConfigurationFile()
        cls.dna.computeTopology(dna_type)
        cls.dna.writePDB()
        cls.system = System(cls.dna)
        cls.system.clearForces()
        cls.system.setDefaultPeriodicBoxVectors(*np.diag([2 * 9.4208000] * 3))

    def _test_energy(self,
                     log_energy='E_bond',
                     log_file='examples/adna/sim.log',
                     traj='examples/adna/traj.xyz'):
        log = parse_log(log_file)

        self.system.clearForces()
        tempforce = Bond3SPN2(self.dna)
        self.system.addForce(tempforce)
        energies = self.system.recomputeEnergy(traj)
        diff = ((energies / _ef - log[log_energy]) ** 2).sum()
        print(f'Testing {log_energy} on DNA ', diff)
        assert diff < 1E-3, diff

    def test_energies(self):
        pandas.read_csv
        for i in range(10):
            yield self._test_energy
        for i in range(10):
            yield self._test_energy


def test_Eangle_harmonic():
    pass


def test_Eangle_stacking():
    pass


def test_Edihedral():
    pass


def test_Ebp():
    pass


def test_Ecstk():
    pass


def test_Eelectro():
    pass
