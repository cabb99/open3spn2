#!/usr/bin/env python

"""
This tool has been constructed to map a PDB file to a CG representation.  It
prints .xyz and .psf files for visualization, as well as the Atoms section for
inclusion in a LAMMPS data file.  It also can print the bonds, angles, and
dihedrals sections (uncomment the lines in main().

Author: Dan Hinckley
Date: 12/16/11
"""

import sys, re

# Defining a molecule and atom object for easy manipulation

# Dictionaries for masses and maping to CG beads

residues = ["DC","DG","DA","DT","A","T","G","C","DC5","DC3","DA5","DA3","DG5","DG3","DT5", "DT3","ADE","THY","GUA","CYT"]

masses = {"H":1.00794,"C":12.0107,"N":14.0067,"O":15.9994,"P":30.973762}
bead_mass = {
    "P":94.9696,
    "S":83.1104,
    "A":134.1220,
    "G":150.1214,
    "T":125.1078,
    "C":110.0964
}
# 0 - sugar
# 1 - base
# 2 - phosphate
CG = {
    "O5\'":2,
    "C5\'":0,
    "C4\'":0,
    "O4\'":0,
    "C3\'":0,
    "O3\'":2,
    "C2\'":0,
    "C1\'":0,
    "O5*":2,
    "C5*":0,
    "C4*":0,
    "O4*":0,
    "C3*":0,
    "O3*":2,
    "C2*":0,
    "C1*":0,

    "N1":1,
    "C2":1,
    "O2":1,
    "N2":1,
    "N3":1,
    "C4":1,
    "N4":1,
    "C5":1,
    "C6":1,
    "N9":1,
    "C8":1,
    "O6":1,
    "N7":1,
    "N6":1,
    "O4":1,
    "C7":1,
    "P":2,
    "OP1":2,
    "OP2":2,
    "O1P":2,
    "O2P":2,
    "OP3":2
}

restype = {
    "A":"ADE",
    "G":"GUA",
    "T":"THY",
    "C":"CYT"
}
siteNumber = {
    "P":1,
    "S":2,
    "A":3,
    "T":4,
    "G":5,
    "C":6
}

bt = {
    "P":{"S":1},
    "S":{
        "P":2,
        "A":3,
        "T":4,
        "G":5,
        "C":6
        }
}
at = {
    "P":{
        "S":{
            "P":2,
            "A":7,
            "T":8,
            "G":9,
            "C":10
            }
        },  
    "S":{"P":{"S":1},
        "A":{
            "A":11,
            "T":12,
            "G":13,
            "C":14},
        "T":{
            "A":15,
            "T":16,
            "G":17,
            "C":18},
        "G":{
            "A":19,
            "T":20,
            "G":21,
            "C":22},
        "C":{
            "A":23,
            "T":24,
            "G":25,
            "C":26}
    },
    "A":{"S":{"P":3}},
    "T":{"S":{"P":4}},
    "G":{"S":{"P":5}},
    "C":{"S":{"P":6}}
}

dt = {
    "S":{"P":{"S":{"P":1,
                    "A":-1,
                    "T":-1,
                    "G":-1,
                    "C":-1
                    }}},
    "P":{"S":{"P":{"S":2}}},
    "A":{"S":{"P":{"S":-1}}},
    "T":{"S":{"P":{"S":-1}}},
    "G":{"S":{"P":{"S":-1}}},
    "C":{"S":{"P":{"S":-1}}}
}

pCharge = -0.6
    
class System:
    def __init__(self):
        self.molecules = []
        self.cgmodel = []
        self.nbonds = None
        self.nbends = None
        self.ntorsions = None
        self.nimprs = 0
        self.hdon = 0
        self.hacc = 0
        self.nnb = 0
        self.nbeads = 0
    def number_molecules(self):
        return len(self.molecules)
   
    def generate_statistics(self):
        for i in range(self.number_molecules()):
            self.molecules[i].sequence()

    def do_coarse_graining(self):
        for i in range(self.number_molecules()):
            self.molecules[i].coarse_grain()

    def generate_topology(self):
        nbeads = 0
        nbonds = 0
        nbends = 0
        ntorsions = 0
        shift = 0
        for i  in range(self.number_molecules()):
            self.molecules[i].shift = shift
            self.molecules[i].topology()

            # Calcualting the total number of each topogical interaction
            nbonds = nbonds + self.molecules[i].number_bonds()
            nbends = nbends + self.molecules[i].number_bends()
            ntorsions = ntorsions + self.molecules[i].number_torsions()
            nbeads = nbeads + self.molecules[i].number_beads()
            shift = shift + self.molecules[i].number_beads()
        self.nbonds =  nbonds  
        self.nbends =  nbends  
        self.nbeads =  nbeads  
        self.ntorsions = ntorsions   
            
    def write_psf(self):
        fptr = open("in00_cvmd.psf","w")
        fptr.write("*\n*\n*\n*\n*\n\n")
        fptr.write("%7ld !NATOMS\n" % self.nbeads)
        bead_shift = 0
        residue_shift = 0
        for i in range(self.number_molecules()):
            for j in range(self.molecules[i].number_beads()):
                bead = self.molecules[i].beads[j]
                charge = 0.0
                if siteNumber[bead.beadtype] == 1:
                    charge = pCharge
                fptr.write("%8ld  DNA%d %-4ld %s  %-5s %3ld  %13.6e   %7.3lf          0\n" % (
                    bead.beadnumber + bead_shift, i, bead.residuenumber + residue_shift, bead.residuetype, bead.beadtype, siteNumber[bead.beadtype], charge, bead_mass[bead.beadtype]))
            if self.molecules[i].number_beads() == 0:
                lastbead = 0
                lastresidue = 0
            else:
                lastbead = bead.beadnumber
                lastresidue = bead.residuenumber
            bead_shift = bead_shift + lastbead
            residue_shift = residue_shift + lastresidue
        
        fptr.write("\n%8ld !NBONDS\n" % self.nbonds)
        count = 0 
        for i in range(self.number_molecules()):
            for j in range(self.molecules[i].number_bonds()):
                count = count + 1
                bond = self.molecules[i].bonds[j]
                fptr.write("%8ld%8ld" %( bond.stea, bond.steb))
                if count < self.nbonds:
                    if (count % 4) == 0:
                        fptr.write("\n")
                else:       
                    fptr.write("\n\n")

        count = 0
        fptr.write("%8ld !NTHETA: angles\n" % self.nbends)
        for i in range(self.number_molecules()):
            for j in range(self.molecules[i].number_bends()):
                count = count + 1
                bend = self.molecules[i].bends[j]
                fptr.write("%8ld%8ld%8ld" % ( bend.stea, bend.steb, bend.stec))
                if count < self.nbends:
                    if (count % 3) == 0:
                        fptr.write("\n")
                else:         
                    fptr.write("\n\n")
    
        count = 0
        fptr.write("%8ld !NPHI: dihedrals\n" % self.ntorsions)
        for i in range(self.number_molecules()):
            for j in range(self.molecules[i].number_torsions()):
                count = count + 1
                torsion = self.molecules[i].torsions[j]
                fptr.write("%8ld%8ld%8ld%8ld" % ( torsion.stea, torsion.steb, torsion.stec, torsion.sted))
                if count < self.ntorsions:
                    if (count % 2) == 0:
                        fptr.write("\n")
                else:        
                    fptr.write("\n\n")

        fptr.write("%8ld !NIMPR\n\n" % 0)
        fptr.write("%8ld !HDON\n\n" % 0)
        fptr.write("%8ld !HACC\n\n" % 0)
        fptr.write("%8ld !NNB\n\n" % 0)

        fptr.close()

    def write_xyz(self):
        fptr = open("in00_conf.xyz",'w')
        nbeads = 0
        for i in range(self.number_molecules()):
            nbeads = nbeads + self.molecules[i].number_beads()
        
        fptr.write("%ld\n\n" % nbeads)
        for i in range(self.number_molecules()):
            for j in range(self.molecules[i].number_beads()):
                bead = self.molecules[i].beads[j]
                fptr.write("%s\t%lf\t%lf\t%lf\n" % (bead.beadtype, bead.x, bead.y, bead.z))
        fptr.close()            

    def write_lammps_cord(self):
        fptr = open("dna_conf.in","w")
        fptr.write("Atoms\n\n")
        nbeads = 0 
        for i in range(self.number_molecules()):
            for j in range(self.molecules[i].number_beads()):
                bead = self.molecules[i].beads[j]
                shift = 0
                if j == 1:
                    shift = 4
                if j == self.molecules[i].number_beads()-1:
                    shift = 8
                fptr.write("\t%ld\t%ld\t%ld\t%lf\t%lf\t%lf\t%lf\n" % (j+nbeads+1, i+1,siteNumber[bead.beadtype]+shift, bead.charge,bead.x, bead.y, bead.z))
            nbeads = nbeads + self.molecules[i].number_beads()
        fptr.close()

    def write_lammps_bonds(self):
        count = 0
        fptr = open("dna_bonds.in","w")
        fptr.write("Bonds\n\n")
        for i in range(self.number_molecules()):
            for j in range(self.molecules[i].number_bonds()):
                count = count + 1
                bond = self.molecules[i].bonds[j]
                fptr.write("%ld\t%ld\t%ld\t%ld\n" % (count,bond.type, bond.stea, bond.steb))
        fptr.close()

    def write_lammps_angles(self):
        count = 0
        fptr = open("dna_angles.in","w")
        fptr.write("Angles\n\n")
        for i in range(self.number_molecules()):
            for j in range(self.molecules[i].number_bends()):
                count = count + 1
                bend = self.molecules[i].bends[j]
                fptr.write("%ld\t%ld\t%ld\t%ld\t%ld\n" % (count,bend.type, bend.stea, bend.steb, bend.stec))
        fptr.close()

    def write_lammps_dihedrals(self):
        count = 0
        fptr = open("dna_dihedrals.in","w")
        fptr.write("Dihedrals\n\n")
        for i in range(self.number_molecules()):
            for j in range(self.molecules[i].number_torsions()):
                torsion = self.molecules[i].torsions[j]
                if (torsion.type > 0):
                    count = count + 1
                    torsion = self.molecules[i].torsions[j]
                    fptr.write("%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n" % (count,torsion.type, torsion.stea, torsion.steb, torsion.stec, torsion.sted))
        fptr.close()
    
    def write_crd(self):
        fptr = open("in00_conf.crd",'w')
        nbeads = 0
        for i in range(self.number_molecules()):
            nbeads = nbeads + self.molecules[i].number_beads()
        
        fptr.write("%ld\n\n" % nbeads)
        bead_shift = 0
        residue_shift = 0
        for i in range(self.number_molecules()):
            for j in range(self.molecules[i].number_beads()):
                bead = self.molecules[i].beads[j]
                fptr.write("%ld\t%ld\t%s\t%s\t%lf\t%lf\t%lf\n" % (bead.beadnumber + bead_shift, bead.residuenumber + residue_shift, bead.residuetype,bead.beadtype, bead.x, bead.y, bead.z))
            if self.molecules[i].number_beads() == 0:
                lastbead = 0
                lastresidue = 0
            else:
                lastbead = bead.beadnumber
                lastresidue = bead.residuenumber
            bead_shift = bead_shift + lastbead
            residue_shift = residue_shift + lastresidue
        fptr.close()            
      

    def print_statistics(self):
        print "There are %d molecules" % self.number_molecules()
        print "Molecule statistics:"
        for i in range(self.number_molecules()):
            print "Molecule - %d :" % (i+1),
            print "Sequence: %s, NumAtoms: %d" % (self.molecules[i].sequence,self.molecules[i].number_atoms())
        filename = "seq"
        fptr = open(filename,"w")
        fptr.write("%d\n" % len(self.molecules[0].sequence))
        for i in range(len(self.molecules[0].sequence)):
            fptr.write("%s" % self.molecules[0].sequence[i])
        fptr.close()

class Molecule:
    def __init__(self):
        self.atoms = []
        self.beads = []
        self.bonds = []
        self.bends = []
        self.torsions = []
        self.shift = 0
    def add_atom(self,atom):
        assert isinstance(atom,Atom)
        if atom.atomtype != "H":
            self.atoms.append(atom)
    def number_atoms(self):
        return len(self.atoms)
    def number_beads(self):
        return len(self.beads)
    def number_bonds(self):
        return len(self.bonds)
    def number_bends(self):
        return len(self.bends)
    def number_torsions(self):
        return len(self.torsions)
    def number_nucleotides(self):
        return len(self.sequence)

    def coarse_grain(self):
        coarsen = 0
        sugar = []
        base = []
        phosphate = []
        beadnumber = 0
        residuenumber = 1
        for i in range(self.number_atoms()):
            # Skipping the capping oxygen
            if coarsen == 0:
                if self.atoms[i].name =="O5\'" or self.atoms[i].name =="O5*":
                    coarsen = 1
                    pass
            else:
                name = self.atoms[i].name
                if CG[name] == 0:
                    sugar.append(self.atoms[i])
                if CG[name] == 1:
                    base.append(self.atoms[i])
                if CG[name] == 2:
                    phosphate.append(self.atoms[i])
                if residuenumber < len(self.sequence):
                    if len(phosphate) == 5:
                        # Write coordinates to file

                        # Sugar
                        residuetype = restype[self.sequence[residuenumber -1]]
                        beadnumber = beadnumber + 1
                        beadtype = 'S'
                        c = get_com(sugar)
                        self.beads.append(Bead(beadnumber,residuenumber,residuetype,beadtype, c[0],c[1],c[2],0))
            
                        # Base
                        residuetype = restype[self.sequence[residuenumber-1]]
                        beadnumber = beadnumber + 1
                        beadtype = self.sequence[residuenumber-1]
                        c = get_com(base)
                        self.beads.append(Bead(beadnumber,residuenumber,residuetype,beadtype, c[0],c[1],c[2],0))

                        residuenumber = residuenumber + 1
                        # Phosphate
                        residuetype = restype[self.sequence[residuenumber-1]]
                        beadnumber = beadnumber + 1
                        beadtype = 'P'
                        c = get_com(phosphate)
                        self.beads.append(Bead(beadnumber,residuenumber,residuetype,beadtype, c[0],c[1],c[2],pCharge))
                        
                        # Resetting the atom lists
                        sugar = []
                        base = []
                        phosphate = []
                else:
                    if (i == self.number_atoms() - 1):
                        # Sugar
                        residuetype = restype[self.sequence[residuenumber -1]]
                        beadnumber = beadnumber + 1
                        beadtype = 'S'
                        c = get_com(sugar)
                        self.beads.append(Bead(beadnumber,residuenumber,residuetype,beadtype, c[0],c[1],c[2],0))
            
                        # Base
                        residuetype = restype[self.sequence[residuenumber-1]]
                        beadnumber = beadnumber + 1
                        beadtype = self.sequence[residuenumber-1]
                        c = get_com(base)
                        self.beads.append(Bead(beadnumber,residuenumber,residuetype,beadtype, c[0],c[1],c[2],0))

    def topology(self):
        k1 = 1 # Indexing in PSF file starts at 1
        shift = self.shift 
        for i in range(self.number_nucleotides()):
            if i == 0:
                # S - Base
                stea = k1
                steb = k1 + 1
                type = bt[self.beads[stea-1].beadtype][self.beads[steb-1].beadtype]
                self.bonds.append(Bond(stea,steb, shift,type))
                # S - P
                stea = k1
                steb = k1 + 2
                type = bt[self.beads[stea-1].beadtype][self.beads[steb-1].beadtype]
                self.bonds.append(Bond(stea ,steb, shift,type))
                k1 = k1 + 2
            else:
                # P - S
                stea = k1
                steb = k1 + 1
                type = bt[self.beads[stea-1].beadtype][self.beads[steb-1].beadtype]
                self.bonds.append(Bond(stea,steb, shift,type))
                # S - Base
                stea = k1 +1
                steb = k1 + 2
                type = bt[self.beads[stea-1].beadtype][self.beads[steb-1].beadtype]
                self.bonds.append(Bond(stea ,steb, shift,type))
                if i != self.number_nucleotides() - 1:
                    # S - P
                    stea = k1 +1
                    steb = k1 + 3
                    type = bt[self.beads[stea-1].beadtype][self.beads[steb-1].beadtype]
                    self.bonds.append(Bond(stea ,steb, shift,type))
                k1 = k1 + 3     
        cben = 0 # Bend counter
        for i in range(self.number_bonds()-1):
            for j in range(i+1,self.number_bonds()):
                istea = self.beads[self.bonds[i].stea- 1 - shift].beadtype
                isteb = self.beads[self.bonds[i].steb- 1 - shift].beadtype
                jstea = self.beads[self.bonds[j].stea- 1 - shift].beadtype
                jsteb = self.beads[self.bonds[j].steb- 1 - shift].beadtype

                if self.bonds[i].stea == self.bonds[j].stea:
                    type = at[isteb][istea][jsteb]
                    self.bends.append(Bend(self.bonds[i].steb,self.bonds[i].stea,self.bonds[j].steb,type))
                    cben = cben + 1
                elif self.bonds[i].stea == self.bonds[j].steb:
                    type = at[isteb][istea][jstea]
                    self.bends.append(Bend(self.bonds[i].steb,self.bonds[i].stea,self.bonds[j].stea,type))
                    cben = cben + 1
                elif self.bonds[i].steb == self.bonds[j].stea:
                    type = at[istea][isteb][jsteb]
                    self.bends.append(Bend(self.bonds[i].stea,self.bonds[i].steb,self.bonds[j].steb,type))
                    cben = cben + 1
                elif self.bonds[i].steb == self.bonds[j].steb:
                    type = at[istea][isteb][jstea]
                    self.bends.append(Bend(self.bonds[i].stea,self.bonds[i].steb,self.bonds[j].stea,type))
                    cben = cben + 1
        
        ctor = 0 
        for i in range(cben-1):
            for j in range(i+1,cben):
                istea = self.beads[self.bends[i].stea-1 - shift].beadtype
                isteb = self.beads[self.bends[i].steb-1 - shift].beadtype
                istec = self.beads[self.bends[i].stec-1 - shift].beadtype
                jstea = self.beads[self.bends[j].stea-1 - shift].beadtype
                jsteb = self.beads[self.bends[j].steb-1 - shift].beadtype
                jstec = self.beads[self.bends[j].stec-1 - shift].beadtype

                if self.bends[i].stec == self.bends[j].steb:
                    if self.bends[i].steb == self.bends[j].stea:
                        type = dt[istea][isteb][istec][jstec]
                        self.torsions.append(Torsion(self.bends[i].stea, self.bends[i].steb,self.bends[i].stec,self.bends[j].stec,type))
                        ctor  = ctor + 1
                    elif self.bends[i].steb == self.bends[j].stec:
                        type = dt[istea][isteb][istec][jstea]
                        self.torsions.append(Torsion(self.bends[i].stea, self.bends[i].steb,self.bends[i].stec,self.bends[j].stea,type))
                        ctor  = ctor + 1
                elif self.bends[i].stea == self.bends[j].steb:
                    if self.bends[i].steb == self.bends[j].stec:
                        type = dt[jstea][jsteb][jstec][istec]
                        self.torsions.append(Torsion(self.bends[j].stea, self.bends[j].steb,self.bends[j].stec,self.bends[i].stec,type))
                        ctor  = ctor + 1
                    elif self.bends[i].steb == self.bends[j].stea:
                        type = dt[jstec][jsteb][jstea][istec]
                        self.torsions.append(Torsion(self.bends[j].stec, self.bends[j].steb,self.bends[j].stea,self.bends[i].stec,type))
                        ctor  = ctor + 1

        
        # Building stacking interactions
        for i in range(self.number_nucleotides()-1):
            stea = 3 * i + 1 + shift
            steb = 3 * i + 2+ shift
            stec = 3 * i + 5 + shift
            istea = self.beads[stea - 1 - shift].beadtype
            isteb = self.beads[steb - 1 - shift].beadtype
            istec = self.beads[stec - 1 - shift].beadtype
            type = at[istea][isteb][istec]
            self.bends.append(Bend(stea,steb,stec,type))


    def sequence(self):
        sequence = []
        residue = None
        reSeq = re.compile(r"^\w(\w)")
        for i in range(len(self.atoms)):
            tmp = self.atoms[i].residuenumber
            if residue != tmp:
                residue = tmp
                matchSeq = reSeq.search(self.atoms[i].residuetype)
                sequence.append(matchSeq.group(1))
        self.nNucleotides = len(sequence)        
        string = ''.join(sequence)
        self.sequence = string 
    
class Atom(object):
    def __init__(self,number,name,residuetype,moleculeletter,residuenumber,x,y,z, atomtype):
        self.number = number
        self.name = name
        self.residuetype = residuetype
        self.moleculeletter = moleculeletter
        self.residuenumber = residuenumber
        self.x = x
        self.y = y
        self.z = z
        self.atomtype = atomtype
        self.mass = masses[atomtype]
   
class Bead(object):
    def __init__(self,beadnumber,residuenumber,residuetype,beadtype, x,y,z,charge):
        self.beadnumber = beadnumber
        self.residuenumber = residuenumber
        self.residuetype = residuetype
        self.beadtype = beadtype
        self.x = x
        self.y = y
        self.z = z
        self.charge = charge

class Bond(object):
    def __init__(self, stea, steb, shift,type):
        self.stea = stea + shift
        self.steb = steb + shift
        self.type = type


class Bend(object):
    def __init__(self, stea, steb,stec,type):
        self.stea = stea
        self.steb = steb
        self.stec = stec
        self.type = type

class Torsion(object):
    def __init__(self, stea, steb, stec, sted,type):
        self.stea = stea
        self.steb = steb
        self.stec = stec
        self.sted = sted
        self.type = type

def get_com(atomlist):
    xcom = 0.0
    ycom = 0.0
    zcom = 0.0
    mtot = 0.0
    for i in range(len(atomlist)):
        xcom = xcom + atomlist[i].x * atomlist[i].mass
        ycom = ycom + atomlist[i].y * atomlist[i].mass
        zcom = zcom +atomlist[i].z * atomlist[i].mass
        mtot = mtot + atomlist[i].mass
    xcom = xcom/mtot
    ycom = ycom/mtot
    zcom = zcom/mtot
    return (xcom,ycom,zcom)

def check_args():
    if len(sys.argv) != 2:
        print "Incorrect command line arguments supplied!"
        print "Correct usage is:"
        print "\t %s <PDB file>" % sys.argv[0]
        sys.exit(1)
    else:
        return sys.argv

def read_pdb(filename,system):
    fptr = open(filename,'r')
    mol_id = None
    count = 0
    for l in fptr.readlines():
        line = l.strip().split()
       # Here I hard code in the relevant space as explained here: http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM
        if line[0] == "ATOM":
            count = count + 1
            number = int(l[6:11])
            name = l[12:16].strip()
            altlocindicator = l[16]
            residuetype = l[17:20].strip()
            moleculeletter = l[21]
            if moleculeletter != mol_id:
                system.molecules.append(Molecule())
                mol_id = moleculeletter
            residuenumber = l[22:26].strip()
            x = float(l[30:38].strip())
            y = float(l[38:46].strip())
            z = float(l[46:54].strip())
            #atomtype = l[76:78].strip()
            atomtype = name[0]
            if name in CG and residuetype in residues and (altlocindicator == "A" or altlocindicator == " "):
                system.molecules[-1].add_atom(Atom(count,name,residuetype,moleculeletter,residuenumber,x,y,z, atomtype))
            else:
                count = count - 1
        elif line[0] == "ENDMDL":
            break
        else:    
            None
    fptr.close()

def main():
    args = check_args()
    filename = args[1]
    # Initializing class object
    system = System()
    read_pdb(filename,system)
    system.generate_statistics()
    #system.print_statistics()
    system.do_coarse_graining()
    system.generate_topology()
    system.write_xyz()
    system.write_crd()
    system.write_psf()
    system.write_lammps_cord()
    system.write_lammps_bonds()
    system.write_lammps_angles()
    system.write_lammps_dihedrals()

if __name__ == "__main__":
    main()
