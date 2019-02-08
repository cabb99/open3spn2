import simtk.openmm.app
import simtk.openmm
from simtk.unit import *
import configparser

"""
This tool has been constructed to simulate DNA using the 3SPN2 forcefield 
in openmm

Author: Carlos Bueno
"""

def parseConfigTable(config_section):
    '''Parses a section of the configuration file as a table'''
    def readData(config_section,a):
        '''Filters comments and returns values as a list'''
        temp=config_section.get(a).split('#')[0].split()
        l=[]
        for val in temp:
            val=val.strip()
            try:
                x=int(val)
                l+=[x]
            except ValueError:
                try:
                    y=float(val)
                    l+=[y]
                except ValueError:
                    l+=[val]
        return l
    data=[]
    for a in config_section:
        if a=='name':
            columns=readData(config_section,a)
        elif len(a)>3 and a[:3]=='row':
            data+=[readData(config_section,a)]
        else:
            print (f'Unexpected row {readData(config_section,a)}')
    return pandas.DataFrame(data,columns=columns)


## Errors
class BaseError(Exception):
    pass

class DNATypeError(BaseError):
    def __init__(self, dna):
        self.dna=dna
        self.message=f'DNA type {dna.DNAType} not defined in the configuration file\n'
        defined_types=dna.angle_definition['DNA'].unique()
        self.message+=f'Only the types {str(defined_types)} were defined'
        print(self.message)


class DNA(object):
    def __init__(self):
        '''Initializes an empty DNA object'''
        pass
    
    def __repr__(self):
        return f'DNA object'
        #print the sequence and the identity of the DNA object

    
    def parseConfigurationFile(self,configuration_file='3SPN2.conf'):
        self.configuration_file=configuration_file
        config = configparser.ConfigParser()
        config.read(configuration_file)
        self.bond_definition=parseConfigTable(config['Bonds'])
        self.angle_definition=parseConfigTable(config['Angles'])
        self.dihedral_definition=parseConfigTable(config['Dihedrals'])

    def computeTopology(self,DNAType='A'):
        #Parse configuration file if not already done
        try:
            self.bond_definition
        except AttributeError:
            self.parseConfigurationFile()
        
        self.DNAType=DNAType
        if DNAType not in self.angle_definition['DNA'].unique():
            raise DNATypeError(self)
        
        
        #Make an index to build the topology
        index={}
        cr_list=set()#Chain residue list
        for i,atom in self.atoms.iterrows():
            index.update({(atom['chain'],atom['residue'],atom['type']):i})
            cr_list.update([(atom['chain'],atom['residue'])])
        cr_list=list(cr_list)
        cr_list.sort()
        max_chain=self.atoms['chain'].max()
        max_residue=self.atoms['residue'].max()
        assert len(index)==len(self.atoms),'Atom index was repeated'

        #Select ADNA bond definitions
        
        angle_types=self.angle_definition[self.angle_definition['DNA']==DNAType]
        if DNAType=='B_curved':
            DNAType='B'
        bond_types=self.bond_definition[self.bond_definition['DNA']==DNAType]
        dihedral_types=self.dihedral_definition[self.dihedral_definition['DNA']==DNAType]

        #Make a table with bonds
        data=[]
        for i,ftype in bond_types.iterrows():
            #print(bond_type)
            ai=ftype['i']
            aj=ftype['j']
            s1=ftype['s1']
            for c,r in cr_list:
                k1=(c,r,ai)
                k2=(c,r+s1,aj)
                if k1 in index and k2 in index:
                    data+=[[i,index[k1],index[k2]]]
        data=pandas.DataFrame(data,columns=['type','aai','aaj'])
        self.bonds=data.merge(bond_types,left_on='type',right_index=True)

        #Make a table with angles
        data=[]
        for i,ftype in angle_types.iterrows():
            #print(bond_type)
            ai=ftype['i']
            aj=ftype['j']
            ak=ftype['k']
            s1=ftype['s1']
            s2=ftype['s2']
            for c,r in cr_list:
                k1=(c,r,ai)
                k2=(c,r+s1,aj)
                k3=(c,r+s2,ak)
                if k1 in index and k2 in index and k3 in index:
                    data+=[[i,index[k1],index[k2],index[k3]]]
        data=pandas.DataFrame(data,columns=['type','aai','aaj','aak'])
        self.angles=data.merge(angle_types,left_on='type',right_index=True)

        #Make a table with dihedrals
        data=[]
        for i,ftype in dihedral_types.iterrows():
            #print(bond_type)
            ai=ftype['i']
            aj=ftype['j']
            ak=ftype['k']
            al=ftype['l']
            s1=ftype['s1']
            s2=ftype['s2']
            s3=ftype['s3']
            for c,r in cr_list:
                k1=(c,r,ai)
                k2=(c,r+s1,aj)
                k3=(c,r+s2,ak)
                k4=(c,r+s3,al)
                if k1 in index and k2 in index and k3 in index and k4 in index:
                    data+=[[i,index[k1],index[k2],index[k3],index[k4]]]
        data=pandas.DataFrame(data,columns=['type','aai','aaj','aak','aal'])
        self.dihedrals=data.merge(dihedral_types,left_on='type',right_index=True)
    
    def writePDB(self,pdb_file='clean.pdb'):
        ''' Writes a minimal version of the pdb file needed for openmm '''
        #Compute chain field
        chain_ix='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
        self.atoms['chainid']=[chain_ix[i-1] for i in self.atoms['chain']]
        #compute resid and resname fields
        res_ix={}
        min_res=self.atoms.groupby('chain')['residue'].min()
        max_res=self.atoms.groupby('chain')['residue'].max()
        for i,res in self.atoms[~self.atoms['type'].isin(['S','P'])].iterrows():
            resname='D'+res['type']
            if res['residue']==min_res[res['chain']]:
                resname+='i'
            if res['residue']==max_res[res['chain']]:
                resname+='f'
            res_ix.update({(res['chain'],res['residue']):resname})
        self.atoms['resname']=[res_ix[(r['chain'],r['residue'])] for i,r in self.atoms.iterrows()]
        
        #Compute element fields
        element_ix={'P':'P','S':'H','A':'N','T':'S','C':'O','G':'C'} #Elements choosen to keep VMD colors
        self.atoms['element'] = [element_ix[atomType] for atomType in self.atoms['type']]

        #Write pdb file
        with open(pdb_file,'w+') as pdb:
            for i,atom in self.atoms.iterrows():
                pdb_line=f'ATOM  {i+1:>5} {atom.type:^4} {atom.resname:<3}  {atom.chainid}{atom.residue:>3}    {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}'+' '*22+f'{atom.element:2}'+' '*2
                assert len(pdb_line)==80,'An item in the atom table is longer than expected'
                pdb.write(pdb_line+'\n')
        self.pdb_file=pdb_file
        
    @classmethod
    def fromPDB(cls,pdb_file):
        '''Initializes a DNA object from a pdb file'''
        #Parse the pdb
        
        #Make a clean pdb file
        
        #Initialize the system from the pdb
        
        pass        
    
    @classmethod
    def fromGRO(cls,gro_file):
        '''Initializes a DNA object from a gromacs input file'''
        #Parse the gromacs file
        
        #Make a clean pdb file
        
        #Initialize the system from the pdb
        
        pass        
    
    @classmethod
    def fromSequence(cls,sequence, center=[0,0,0]):
        '''Initializes a DNA object from a DNA sequence'''
        #Make a possible structure
        
        #Make a clean pdb file
        
        #Initialize the system from the pdb
        
        pass
        
    @classmethod
    def fromXYZ(cls, xyz_file):
        '''Initializes DNA object from xyz file (as seen on the examples)'''
        #Parse the file
        self=cls()
        self.atoms=pandas.read_csv('examples/adna/in00_conf.xyz',delim_whitespace=True,skiprows=2,names=['type','x','y','z'])
        
        #Compute residues and chains
        residue=0
        residues=[]
        chain=0
        chains=[]
        sugar=True
        for t in self.atoms['type']:
            if t=='S':
                if sugar:
                    residue+=1
                    chain+=1
                sugar=True
            if t=='P':
                residue+=1
                sugar=False
            residues+=[residue]
            chains+=[chain]
        self.atoms['residue']=residues
        self.atoms['chain']=chains
        return self

class System(simtk.openmm.System):
    ''' Wrapper of openmm system class, adds some openmm simulation attributes'''    
    def __init__(self,mol,forcefieldFiles=['3SPN2.xml'],periodicBox=None):
        self.top=simtk.openmm.app.PDBFile(mol.pdb_file)
        self.coord=simtk.openmm.app.PDBFile(mol.pdb_file)
        self.forcefield=simtk.openmm.app.ForceField(*forcefieldFiles)
        self._wrapped_system = self.forcefield.createSystem(self.top.topology)
        self.periodicBox=periodicBox
        if periodicBox is not None:
            self._wrapped_system.setDefaultPeriodicBoxVectors(*np.diag(self.periodic_box))
    
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self._wrapped_system, attr)
    
    def clearForces(self,keepCMMotionRemover=True):
        ''' Removes all the forces from the system '''
        for i,f in enumerate(self.getForces()):
            if keepCMMotionRemover and i==0 and f.__class__==simtk.openmm.CMMotionRemover:
                continue
            else:
                self.removeForce(0)
        if keepCMMotionRemover==False:
            assert len(self.getForces())==0, 'Not all the forces were removed'
        else:
            assert len(self.getForces())<=1, 'Not all the forces were removed'
            
    def initializeMD(self,temperature=300*kelvin):
        self.integrator=simtk.openmm.LangevinIntegrator(temperature, 1E-4/picosecond, 0.1*picoseconds)
        platform = simtk.openmm.Platform.getPlatformByName('Reference')
        self.simulation = simtk.openmm.app.Simulation(self.top.topology, self._wrapped_system, self.integrator,platform)
        self.simulation.context.setPositions(self.coord.positions)
        return self.simulation
    
    def setPositions(self, coords=None):
        #Initialize trial MD if not setup
        try:
            self.simulation
        except AttributeError:
            self.initializeMD()
        
        #Set up coords for MD
        if coords is None:
            self.simulation.context.setPositions(self.coord.positions)
        else:
            self.simulation.context.setPositions(coords)
            
    def getPotentialEnergy(self, coords=None, unit=kilojoule_per_mole):
        #Initialize trial MD if not setup
        try:
            self.simulation
        except AttributeError:
            self.initializeMD()
        self.setPositions(coords)
        state = self.simulation.context.getState(getEnergy=True)
        energy= state.getPotentialEnergy().value_in_unit(unit)
        return energy
        
class Force(object):
    def computeEnergy(self, system, trajectory):
        #Parse trajectory
        traj=ff3SPN2.parse_xyz('examples/adna/traj.xyz')
        
        #clear all forces on the system
        system.clearForces()
        #setup the force
        self.setUpInteraction()
        #for each item of the table:
            #add the force item
            
        #compute the energy for every frame
        
        #return a Table with the energy
        
    def computeSingleEnergy(self, system, trajectory):
        #Parse trajectory
        traj=ff3SPN2.parse_xyz('examples/adna/traj.xyz')
        #for each item of the table:
            #clear all forces on the system
            #system.clearForces()
            #setup the force
            #self.setUpInteraction()
            #add the force item
            
        #compute the energy for every frame
        
        #return a table with the energy
        
   
    def setUpInteraction(self):
        ''' Creates a new force instance '''
        pass
    
    def defineInteractions(self, parameters):
        ''' Adds the parameter for the force '''
        pass
        
    def addForce(system):
        try:        
            system.addForce(self.force)
        except AttributeError:
            print (' The interaction has not been set up ')
            raise AttributeError
    
        
class Bond3SPN2(Force,simtk.openmm.CustomBondForce):
    def __init__(self,parameters,periodic=False):
        self.reset()
        self.defineInteraction(parameters,periodic)
    
    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.force, attr)
        
    def defineInteraction(self,parameters,periodic):
        cf=4.184 #Conversion factor from Kcal (lammps) to KJ (openmm)
        for i,b in parameters.iterrows():
            #Units converted from 
            parameters=[b['x0']/10., b['K2']*100*cf, b['K3']*1000*cf, b['K4']*10000*cf]
            self.force.addBond(int(b['aai']), int(b['aaj']), parameters)
        bondForce.setUsesPeriodicBoundaryConditions(periodic)
        
    def reset():
        bondForce=simtk.openmm.CustomBondForce("Kb2*(r-r0)^2+Kb3*(r-r0)^3+Kb4*(r-r0)^4")
        bondForce.addPerBondParameter('r0')
        bondForce.addPerBondParameter('Kb2')
        bondForce.addPerBondParameter('Kb3')
        bondForce.addPerBondParameter('Kb4')
        self.force=bondForce
        
        self.perInteractionParameters=[]
        self.GlobalParameters=[]
        for i in range(bondForce.getNumPerBondParameters()):
            self.perInteractionParameter+=[bondForce.getPerBondParameterName(i)]
        for i in range(bondForce.getNumGlobalParameters()):
            self.GlobalParameter+=[bondForce.getGlobalParameterName(i)]
        self.perInteractionParameters, self.GlobalParameters
        

class Eangle_harmonic(Energy):
    pass

class EAngle_3SPN2_Stacking(Energy):
    pass
    
class Edihedral(Energy):
    pass

class Ebp(Energy):
    pass

class Ecstk(Energy):
    pass

class Eelectro(Energy):
    pass
    

#Unit testing
def test_DNA_from_pdb():
    pass
    
def test_DNA_from_gro():
    pass
    
def test_DNA_from_seq():
    pass
    
def test_DNA_from_xyz():
    mol=DNA.from_xyz('examples/adna/in00_conf.xyz')
    assert mol.atoms.at[8,'type']=='P'
    assert round(mol.atoms.at[188,'y'],6)==-8.779343
 
#Functional testing
import pandas
def parse_xyz(filename=''):
    import pandas
    columns=['N','timestep','id','type','x','y','z']
    data=[]
    with open(filename,'r') as traj_file:
        atom=pandas.Series(index=columns)
        atom['id']=None
        for line in traj_file:
            s=line.split()
            if len(s)==1:
                atom['N']=int(s[0])
                if atom['id']>-1:
                    assert atom['id']==atoms
                atoms=int(s[0])
            elif len(s)==3:
                atom['timestep']=int(s[2])
                atom['id']=0
            elif len(s)>3:
                atom['type']=int(s[0])
                atom['x'],atom['y'],atom['z']=[float(a) for a in s[1:4]]
                data+=[atom.copy()]
                atom['id']+=1
    xyz_data=pandas.concat(data,axis=1).T
    for i in ['N','timestep','id','type']:
        xyz_data[i]=xyz_data[i].astype(int)
    return xyz_data

def parse_log(filename=''):
    columns=''
    log_data=[]
    with open(filename,'r') as log_file:
        start=False
        for line in log_file:
            if line[:4]=='Step':
                columns=line.split()
                start=True
                continue
            if start:
                try:
                    log_data+=[[float(a) for a in line.split()]]
                except ValueError:
                    break
    log_data=pandas.DataFrame(log_data,columns=columns)
    
    for i in ['Step','nbp']:
        log_data[i]=log_data[i].astype(int)
    return log_data

def test_parse_xyz():
    xyz_data=parse_xyz('examples/adna/traj.xyz')
    assert xyz_data.at[1,'type']==7
    assert xyz_data.at[1,'x']==4.34621

def test_parse_log():
    log_data=parse_log('examples/adna/sim.log')
    assert log_data.at[1,'Step']==2000
    assert log_data.at[1,'eexcl']==0.45734636

def test_Ebond():
    traj=ff3SPN2.parse_xyz('examples/adna/traj.xyz')
    log=ff3SPN2.parse_log('examples/adna/sim.log')
    simulation_platform='Reference'
    temperature=300*kelvin
    integrator = simtk.openmm.LangevinIntegrator(temperature, 1E-4/picosecond, 0.1*picoseconds)

    platform = simtk.openmm.Platform.getPlatformByName(simulation_platform)
    simulation = simtk.openmm.app.Simulation(DNA.top.topology, DNA.system, integrator,platform)
    energies=[]
    for time,snapshot in traj.groupby('timestep'):
        #print (time,snapshot[['x','y','z']])
        simulation.context.setPositions(np.array(snapshot[['x','y','z']])*nanometers/10)
        state = simulation.context.getState(getEnergy=True)
        energy=state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        #print (time,energy)
        energies+=[energy/4.184]
    log['openmm_bond']=energies
    log.plot(x='Step',y=['E_bond','openmm_bond'])
    plt.savefig('bond_energy.png')
    
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




        


