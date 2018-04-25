#!/usr/bin/env python

"""
This scripts reads in an existing LAMMPS configuration file contains a 3SPN.2
configuration and adds ions according to the commands specified on the 
command line.

Usage:

Please use at your own risk.  No guarantee or warranty is implied. 
Usage: <LAMMPS configuation file> <NaCl [mM]> <MgCl2 mM> <SpdCl3 [mM] > <output file>

By Dan Hinckley (hinckley.dan@gmail.com)
"""

import sys, math
import random

# Coordianates of Spd
spdCoords = [[5.79,0.0,0.0],[0.0,0.0,0.0],[-4.72,0.0,0.0]]

# Defining a few constant values
siteDict = {
    1:"P",
    2:"S",
    3:"A",
    4:"T",
    5:"G",
    6:"C",
    7:"A",
    8:"T",
    9:"G",
    10:"C",
    11:"A",
    12:"T",
    13:"G",
    14:"C",
    15:"Na",
    16:"Mg",
    17:"Cl",
    18:"Spd",
}

baseDict = {"A":0,"T":1,"C":2,"G":3}
FOUR = 4
THREE = 3
kJ2kCal = 1.0/4.184
mM2Ions = 6.0221415E-7
printFreq = 20


# Defining classes to facilitate analysis
class Atom:
    def __init__(self,mol,type,q,x,y,z):
        self.mol = mol
        self.type = type
        self.q = q
        self.x = x
        self.y = y
        self.z = z

class Bond:
    def __init__(self,type,a,b):
        self.type = type
        self.a = a
        self.b = b

class Angle:
    def __init__(self,type,a,b,c):
        self.type = type
        self.a = a
        self.b = b
        self.c = c

class Dihedral:
    def __init__(self,type,a,b,c,d):
        self.type = type
        self.a = a
        self.b = b
        self.c = c
        self.d = d


class LmpConf:
    def __init__(self,args):
        self.xlo = None
        self.xhi = None
        self.ylo = None
        self.yhi = None
        self.zlo = None
        self.yhi = None
        self.nAtoms = None
        self.nBonds = None
        self.nAngles = None
        self.nDihedrals = None
        self.nAtomTypes = None
        self.nBondTypes = None
        self.nAngleTypes = None
        self.nDihedralTypes = None
        self.confFile = args[1]
        self.NaClConc = float(args[2])
        self.MgCl2Conc = float(args[3])
        self.SpdConc = float(args[4])
        self.outFile = args[5]
        self.mass = []
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.bonddata = []
        self.angledata = []
        self.dihedraldata = []
        self.blen = []
        self.bhlf = []
        self.blni = []
        self.nNa_orig = 0
        self.nMg_orig = 0
        self.nCl_orig = 0
        self.nSpd_orig = 0

    def initialize_box(self):
        self.blen.append(abs(self.xhi-self.xlo))
        self.blen.append(abs(self.yhi-self.ylo))
        self.blen.append(abs(self.zhi-self.zlo))
        for i in range(THREE):
            self.bhlf.append(self.blen[i]/2.0)
            self.blni.append(1.0/self.blen[i])
    def apply_pbcs(self,dr,dim):
        if dr > self.bhlf[dim]:
            dr = dr - self.blen[dim] * float(int(dr * self.blni[dim] + 0.5))
        elif dr < -self.bhlf[dim]:
            dr = dr - self.blen[dim] * float(int(dr * self.blni[dim] - 0.5))
        return dr

    def read_configuration_file(self):
        try:
            file = open(self.confFile,"r")
            line = file.readline()
            while line:
                if "atoms" in line:
                    self.nAtoms = int(line.strip().split()[0])
                if "bonds" in line:
                    self.nBonds = int(line.strip().split()[0])
                if "angles" in line:
                    self.nAngles = int(line.strip().split()[0])
                if "dihedrals" in line:
                    self.nDihedrals = int(line.strip().split()[0])
                if "atom types" in line:
                    self.nAtomTypes = int(line.strip().split()[0])
                if "bond types" in line:
                    self.nBondTypes = int(line.strip().split()[0])
                if "angle types" in line:
                    self.nAngleTypes = int(line.strip().split()[0])
                if "dihedral types" in line:
                    self.nDihedralTypes = int(line.strip().split()[0])
                if "xlo xhi" in line:
                    self.xlo = float(line.strip().split()[0])
                    self.xhi = float(line.strip().split()[1])
                if "ylo yhi" in line:
                    self.ylo = float(line.strip().split()[0])
                    self.yhi = float(line.strip().split()[1])
                if "zlo zhi" in line:
                    self.zlo = float(line.strip().split()[0])
                    self.zhi = float(line.strip().split()[1])
                if "Masses" in line :
                    file.readline()
                    for i in range(self.nAtomTypes):
                        l = file.readline().strip().split()
                        self.mass.append(float(l[1]))
                if "Pair Coeffs" in line:
                    print "This script is not designed to read Pair Coeffs data!"
                    sys.exit(1)
                if "Bond Coeffs" in line:
                    file.readline()
                    for i in range(self.nBondTypes):
                        l = file.readline()
                        self.bonddata.append(l)
                if "Angle Coeffs" in line:
                    file.readline()
                    for i in range(self.nAngleTypes):
                        l = file.readline()
                        self.angledata.append(l)
                if "Dihedral Coeffs" in line:
                    file.readline()
                    for i in range(self.nDihedralTypes):
                        l = file.readline()
                        self.dihedraldata.append(l)

                if "Atoms" in line:
                    file.readline()
                    for i in range(self.nAtoms):
                        l = file.readline().strip().split()
                        self.atoms.append(Atom(int(l[1]),int(l[2]),float(l[3]),float(l[4]), float(l[5]),float(l[6])))
                if "Bonds" in line:
                    file.readline()
                    for i in range(self.nBonds):
                        l = file.readline().strip().split()
                        self.bonds.append(Bond(int(l[1]),int(l[2]),int(l[3])))

                if "Angles" in line:
                    file.readline()
                    for i in range(self.nAngles):
                        l = file.readline().strip().split()
                        self.angles.append(Angle(int(l[1]),int(l[2]),int(l[3]),int(l[4])))

                if "Dihedrals" in line:
                    file.readline()
                    for i in range(self.nDihedrals):
                        l = file.readline().strip().split()
                        self.dihedrals.append(Dihedral(int(l[1]),int(l[2]),int(l[3]),int(l[4]),int(l[5])))
                line = file.readline()

        except IOError:
            print "Count not open %s" % self.confFile
            sys.exit(2)
        self.initialize_box()
        file.close()


    def change_phosphate_charge(self):
        print "Changing phosphate charge to -1.0 (if necessary)"
        nChanged = 0
        for i in range(self.nAtoms):
            if siteDict[self.atoms[i].type] == "P":
                if self.atoms[i].q != -1.0:
                    self.atoms[i].q = -1.0
                    nChanged += 1
                    
        if nChanged:
            "Adjusted charge of {} phosphate".format(nChanged)

    def count_existing_ions(self):
        for i in range(len(self.atoms)):
            if siteDict[self.atoms[i].type] == "Na":
                self.nNa_orig += 1
            if siteDict[self.atoms[i].type] == "Mg":
                self.nMg_orig += 1
            if siteDict[self.atoms[i].type] == "Cl":
                self.nCl_orig += 1
            if siteDict[self.atoms[i].type] == "Spd":
                self.nSpd_orig += 1
    
    def report_concentration(self):
        nNa = 0
        nCl = 0
        nMg = 0
        nSpd = 0
        for i in range(len(self.atoms)):
            if siteDict[self.atoms[i].type] == "Na":
                nNa += 1
            if siteDict[self.atoms[i].type] == "Mg":
                nMg += 1
            if siteDict[self.atoms[i].type] == "Cl":
                nCl += 1
            if siteDict[self.atoms[i].type] == "Spd":
                nSpd += 1

        vBox = self.blen[0] * self.blen[1] * self.blen[2]

        print "Final Box concentrations [mM] are as follows (includes counterions):"
        print "Na: {}".format(nNa/(mM2Ions * vBox))
        print "Mg: {}".format(nMg/(mM2Ions * vBox))
        print "Cl: {}".format(nCl/(mM2Ions * vBox))
        print "Spd: {}".format(nSpd/3/(mM2Ions * vBox))

    def add_ions(self):
        dcut = 5

        # Get the box volume
        vBox = self.blen[0] * self.blen[1] * self.blen[2]

        # Get the total charge of the box
        qtot = 0.0
        for i in range(len(self.atoms)):
            qtot += self.atoms[i].q

        print "Total Charge = %lf" % qtot
        if qtot < 0:
            nNa = int(abs(qtot) + int(mM2Ions * self.NaClConc *vBox))
            nCl = int(mM2Ions * self.NaClConc * vBox) + 2*int(mM2Ions * self.MgCl2Conc * vBox)+3*int(mM2Ions * self.SpdConc *vBox)
        else:
            nNa = int(mM2Ions * self.NaClConc *vBox)
            nCl = int(abs(qtot) + int(mM2Ions * self.NaClConc * vBox) + 2*int(mM2Ions * self.MgCl2Conc * vBox) + 3*int(mM2Ions * self.SpdConc *vBox))

        nMg = int(mM2Ions * self.MgCl2Conc * vBox)

        nSpd = int(mM2Ions * self.SpdConc *vBox)


        totCharge = nNa + -1.0 * nCl + 2.0 * nMg + 3.0 * nSpd + qtot 

        nNa2Add = nNa - self.nNa_orig
        if not nNa2Add:
            print "not adding any Na"

        nMg2Add = nMg - self.nMg_orig
        if not nMg2Add:
            print "not adding any Mg"

        nCl2Add = nCl - self.nCl_orig
        if not nCl2Add:
            print "not adding any Cl"

        nSpd2Add = nSpd - self.nSpd_orig
        if not nSpd2Add:
            print "not adding any Spd"


        if totCharge != 0.0:
            print "final charge will not be zero! (%f)" % totCharge
            if totCharge < 0.0:
                for i in range(int(abs(totCharge))):
                    nNa += 1
            if totCharge > 0.0:
                for i in range(int(totCharge)):
                    nCl += 1

        x = []
        dx = []
        for i in range(THREE):
            x.append(0.0)
            dx.append(0.0)

        for i in range(nNa2Add):
            flag = 1;
            if not (i % printFreq):
                print "added %d Na" % i
            while (flag):
                flag = 0
                x[0] = float(random.random()) * self.blen[0] + self.xlo
                x[1] = float(random.random()) * self.blen[1] + self.ylo
                x[2] = float(random.random()) * self.blen[2] + self.zlo
                r = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
                if not flag:
                    for j in range(len(self.atoms)):
                        dx[0] = x[0] - self.atoms[j].x
                        dx[1] = x[1] - self.atoms[j].y
                        dx[2] = x[2] - self.atoms[j].z
                        for k in range(THREE):
                            dx[k] = self.apply_pbcs(dx[k],k)

                        dr = math.sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2])
                        if dr < dcut:
                            flag = 1
            self.atoms.append(Atom(self.atoms[-1].mol+1,15,1.0,x[0],x[1],x[2]))

        for i in range(nMg2Add):
            flag = 1;
            if not (i % printFreq):
                print "added %d Mg" % i
            while (flag):
                flag = 0
                x[0] = float(random.random()) * self.blen[0] + self.xlo
                x[1] = float(random.random()) * self.blen[1] + self.ylo
                x[2] = float(random.random()) * self.blen[2] + self.zlo
                r = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
                if not flag:
                    for j in range(len(self.atoms)):
                        dx[0] = x[0] - self.atoms[j].x
                        dx[1] = x[1] - self.atoms[j].y
                        dx[2] = x[2] - self.atoms[j].z
                        for k in range(THREE):
                            dx[k] = self.apply_pbcs(dx[k],k)

                        dr = math.sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2])
                        if dr < dcut:
                            flag = 1
            self.atoms.append(Atom(self.atoms[-1].mol+1,16,2.0,x[0],x[1],x[2]))

        for i in range(nSpd2Add):
            flag = 1;
            if not (i % printFreq):
                print "added %d Spd" % i
            while (flag):
                flag = 0

                # Copy coords locally
                tmpCoords = [[0,0,0],[0,0,0],[0,0,0]]
                for j in range(THREE):
                    for k in range(THREE):
                        tmpCoords[j][k] = spdCoords[j][k]
                
                # Rotate around each of the axes
                for j in range(THREE):
                    angle = float(random.random())*2.0 * math.pi - math.pi
                    for k in range(THREE):
                        tmpCoords[k] = rotate(j,tmpCoords[k],angle)

                # Apply a random displacement
                x[0] = float(random.random()) * self.blen[0] + self.xlo
                x[1] = float(random.random()) * self.blen[1] + self.ylo
                x[2] = float(random.random()) * self.blen[2] + self.zlo
                for j in range(THREE):
                    for k in range(THREE):
                        tmpCoords[j][k] += x[k]
    
                # Calculate distance between adjacent ions
                for j in range(THREE):
                    for k in range(len(self.atoms)):
                        dx[0] = tmpCoords[j][0] - self.atoms[k].x
                        dx[1] = tmpCoords[j][1] - self.atoms[k].y
                        dx[2] = tmpCoords[j][2] - self.atoms[k].z
                        for l in range(THREE):
                            dx[l] = self.apply_pbcs(dx[l],l)
                        dr = math.sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2])
                        if dr < dcut:
                            flag = 1
            
            # Add the three atoms
            molNumber = self.atoms[-1].mol + 1
            for j in range(THREE):
                self.atoms.append(Atom(molNumber,18,1.0,tmpCoords[j][0],tmpCoords[j][1],tmpCoords[j][2]))
            # Add bonds and angles, hmmmm..
            self.nBonds += 2
            self.bonds.append(Bond(7,len(self.atoms)-2, len(self.atoms)-1))
            self.bonds.append(Bond(8,len(self.atoms)-1, len(self.atoms)-0))

            self.nAngles += 1
            self.angles.append(Angle(27,len(self.atoms)-2, len(self.atoms)-1, len(self.atoms)-0))

        for i in range(nCl2Add):
            flag = 1;
            if not (i % printFreq):
                print "added %d Cl" % i
            while (flag):
                flag = 0
                x[0] = float(random.random()) * self.blen[0] + self.xlo
                x[1] = float(random.random()) * self.blen[1] + self.ylo
                x[2] = float(random.random()) * self.blen[2] + self.zlo
                r = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
                if not flag:
                    for j in range(len(self.atoms)):
                        dx[0] = x[0] - self.atoms[j].x
                        dx[1] = x[1] - self.atoms[j].y
                        dx[2] = x[2] - self.atoms[j].z
                        for k in range(THREE):
                            dx[k] = self.apply_pbcs(dx[k],k)

                        dr = math.sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2])
                        if dr < dcut:
                            flag = 1
            self.atoms.append(Atom(self.atoms[-1].mol+1,17,-1.0,x[0],x[1],x[2]))





        # Get the total charge of the box
        qtot = 0.0
        for i in range(len(self.atoms)):
            qtot += self.atoms[i].q
        if qtot != 0.0:
            print "Box is not neutral after adding ions! (q={}".format(qtot)
            sys.exit(1)

    def write_configuration_file(self):
        file = open(self.outFile,"w")
        file.write("LAMMPS configuration file with ions added\n\n")
        file.write("\t%ld atoms\n" % len(self.atoms))
        file.write("\t%ld bonds\n" % self.nBonds)
        file.write("\t%ld angles\n" % self.nAngles)
        file.write("\t%ld dihedrals\n\n" % self.nDihedrals)

        file.write("\t%ld atom types\n" % self.nAtomTypes)
        file.write("\t%ld bond types\n" % self.nBondTypes)
        file.write("\t%ld angle types\n" % self.nAngleTypes)
        file.write("\t%ld dihedral types\n\n" % self.nDihedralTypes)
        file.write("\t%lf\t%lf xlo xhi\n" % (self.xlo,self.xhi))
        file.write("\t%lf\t%lf ylo yhi\n" % (self.ylo,self.yhi))
        file.write("\t%lf\t%lf zlo zhi\n\n" % (self.zlo,self.zhi))

        file.write("Masses\n\n")
        for i in range(len(self.mass)):
            file.write("\t%d\t%lf\n" %(i+1,self.mass[i]))
        file.write("\n")
        file.write("Bond Coeffs\n\n")
        for i in range(len(self.bonddata)):
            file.write("%s" % self.bonddata[i])
        file.write("\n")
        file.write("Angle Coeffs\n\n")
        for i in range(len(self.angledata)):
            file.write("%s" % self.angledata[i])
        file.write("\n")
        file.write("Dihedral Coeffs\n\n")
        for i in range(len(self.dihedraldata)):
            file.write("%s" % self.dihedraldata[i])
        file.write("\n")
        file.write("Atoms\n\n")
        for i in range(len(self.atoms)):
            a = self.atoms[i]
            file.write("\t%ld\t%ld\t%ld\t%lf\t%lf\t%lf\t%lf\n" % (i+1,a.mol,a.type,a.q,a.x,a.y,a.z))
        file.write("\n")
        file.write("Bonds\n\n")
        for i in range(len(self.bonds)):
            b = self.bonds[i]
            file.write("\t%ld\t%ld\t%ld\t%ld\n" % (i+1,b.type,b.a,b.b))
        file.write("\n")
        file.write("Angles\n\n")
        for i in range(len(self.angles)):
            a = self.angles[i]
            file.write("\t%ld\t%ld\t%ld\t%ld\t%ld\n" % (i+1,a.type,a.a,a.b,a.c))
        file.write("\n")
        file.write("Dihedrals\n\n")
        for i in range(len(self.dihedrals)):
            d = self.dihedrals[i]
            file.write("\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n" % (i+1,d.type,d.a,d.b,d.c,d.d))
        file.close()
    
    def write_xyz(self):
        file = open("in00_conf.xyz","w")
        file.write("%ld\n\n" % len(self.atoms))
        for i in range(len(self.atoms)):
            atom = self.atoms[i]
            file.write("\t%s\t%lf\t%lf\t%lf\n" % (siteDict[atom.type],atom.x,atom.y,atom.z))
        file.close()
        
def rotate(axis,coords,angle):
    cosT = math.cos(angle)
    sinT = math.sin(angle)
    rotated = [0.0,0.0,0.0]
    if (axis == 0): # Rotate x
        matrix = [[1.0,0.0,0.0],[0.0,cosT,-sinT],[0,sinT,cosT]]
    elif (axis == 1): # Rotate y
        matrix = [[cosT,0.0,sinT],[0.0,1.0,0.0],[-sinT,0.0,cosT]]
    elif (axis == 2): # Rotate z
        matrix = [[cosT,-sinT,0.0],[sinT,cosT,0],[0,0,1.0]]
    else:
        "Did not recognize axis value in rotate()!"
        sys.exit(1)
    for i in range(THREE):
        for j in range(THREE):
            rotated[i] += matrix[i][j] *coords[j]
    return rotated

def get_command_line_args(args):
    if len(args) != 6:
        print "Usage: %s <LAMMPS configuation file> <NaCl [mM]> <MgCl2 mM> <SpdCl3 [mM] > <output file>" % args[0]
        sys.exit(1)
    return args

def main():
    args = get_command_line_args(sys.argv)
    lmpconf = LmpConf(args)
    lmpconf.read_configuration_file()
    lmpconf.change_phosphate_charge()
    lmpconf.count_existing_ions()
    lmpconf.add_ions()
    lmpconf.report_concentration()
    lmpconf.write_configuration_file()
    lmpconf.write_xyz()

if __name__ == "__main__":
    main()
