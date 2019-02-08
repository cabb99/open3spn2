import simtk.openmm.app
import simtk.openmm
import simtk.unit

class DNA(object):
    def __init__(self):
        '''Initializes an empty DNA object'''
        pass
    
    def __repr__(self):
        print(f'DNA object')
        #print the sequence and the identity of the DNA object
        pass
    
    @classmethod
    def from_pdb(cls,pdb_file):
        '''Initializes a DNA object from a pdb file'''
        #Parse the pdb
        
        #Make a clean pdb file
        
        #Initialize the system from the pdb
        
        pass        
    
    @classmethod
    def from_gro(cls,gro_file):
        '''Initializes a DNA object from a gromacs input file'''
        #Parse the gromacs file
        
        #Make a clean pdb file
        
        #Initialize the system from the pdb
        
        pass        
    
    @classmethod
    def from_seq(cls,sequence, center=[0,0,0]):
        '''Initializes a DNA object from a DNA sequence'''
        #Make a possible structure
        
        #Make a clean pdb file
        
        #Initialize the system from the pdb
        
        pass
        
    @classmethod
    def from_xyz(cls, xyz_file):
        '''Initializes DNA object from xyz file (as seen on the examples)'''
        #Parse the file
        self.atoms=pandas.read_csv('examples/adna/in00_conf.xyz',delim_whitespace=True,skiprows=2,names=['type','x','y','z'])
        
    
            

class Energy:
    def compute_energy(self,system):
        pass
        
    def compute_single_Energy(self,system,i1,i2,i3):
        pass

class Ebond(Energy):
    ''' '''
    def set_forces():
        force = simtk.openmm.CustomBondForce('K2*(r-r0)^2+K3*(r-r0)^2+K4*(r-r0)^2')
        force.addPerBondParameter('K2')
        force.addPerBondParameter('K3')
        force.addPerBondParameter('K4')
        harmonic_bond.setUsesPeriodicBoundaryConditions(True)
        for i,b in self.sbonds.iterrows():
            harmonic_bond.addBond(int(b['i']),int(b['j']),b['x0']/10.,b['k']*4.184*100)
        self.system.addForce(harmonic_bond)

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
    pass
    
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




        


