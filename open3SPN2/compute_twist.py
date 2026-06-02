#define these functions

import numpy as np
import math

def vector(p1, p2):
    return [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]

def vadd(p1,p2):
    return [p2[0]/2+p1[0]/2, p2[1]/2+p1[1]/2, p2[2]/2+p1[2]/2]

def vmulti(p1,k):
    return [p1[0]*k,p1[1]*k,p1[2]*k]

def vabs(a):
    return math.sqrt(pow(float(a[0]),2)+pow(float(a[1]),2)+pow(float(a[2]),2))

def v_product(p1,p2):
    return p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2]

def vcross_product(a, b):
    cx = a[1]*b[2]-a[2]*b[1]
    cy = a[2]*b[0]-a[0]*b[2]
    cz = a[0]*b[1]-a[1]*b[0]
    return [cx, cy, cz]

#Various methods of calculations of twist. These are similar, using slightly different methods. But they should yield the same results.
def Xun_twist_OG(b_atoms):  #original Xun's code
    dna_length = int(len(b_atoms)/2)
    b1_atoms  =  b_atoms[0:dna_length]
    b2_atoms  =  b_atoms[dna_length:dna_length*2]
    Twist = ""
    for i in range(dna_length-1):
        b_atom_1 = b_atoms[i]
        b_atom_2 = b_atoms[-1-i]
        b_atomn_1 = b_atoms[i+1]
        b_atomn_2 = b_atoms[-1-i-1]
        b_center = vadd(b_atom_1,b_atom_2)
        b_centern = vadd(b_atomn_1,b_atomn_2)
        zmst=vector(b_center,b_centern)
        ys = vector(b_atom_1,b_atom_2)
        ysn = vector(b_atomn_1,b_atomn_2)
        yr = vector(ys,vmulti(zmst,v_product(ys,zmst)))
        yrn = vector(ysn,vmulti(zmst,v_product(ysn,zmst)))
        zmsto = vcross_product(yr,yrn)
        ozmst = v_product(zmst,zmsto)
        #print v_product(yr,yrn)/(vabs(yr)*vabs(yrn))
        twist = math.acos(v_product(yr,yrn)/(vabs(yr)*vabs(yrn)))/3.15159*180
        if twist > 180:
           twist = twist - 180
        elif twist < -180:
           twist = twist + 180
        if ozmst > 0:
           twist = twist
        else:
           twist = 0 - twist
        Twist += str(twist) + " "
    return Twist

def Xun_twist(Vec_1, Vec_2, Vec_3, Vec_4):  #Steven's implementation of Xun's code!! 
    #WARNING: Sign is not taken into account but under in vivo conditions, should be but an edge case!
    b_atom_1 = [float(x)/vabs(Vec_1[0]) for x in Vec_1[0]]
    b_atom_2 = [float(x)/vabs(Vec_2[0]) for x in Vec_2[0]]
    b_atomn_1 = [float(x)/vabs(Vec_3[0]) for x in Vec_3[0]]
    b_atomn_2 = [float(x)/vabs(Vec_4[0]) for x in Vec_4[0]]
    b_center = vadd(b_atom_1,b_atom_2)
    b_centern = vadd(b_atomn_1,b_atomn_2)
    zmst=vector(b_center,b_centern)
    ys = vector(b_atom_1,b_atom_2)
    ysn = vector(b_atomn_1,b_atomn_2)
    yr = vector(ys,vmulti(zmst,v_product(ys,zmst)))
    yrn = vector(ysn,vmulti(zmst,v_product(ysn,zmst)))
    #zmsto = vcross_product(yr,yrn)
    #ozmst = v_product(zmst,zmsto)
    #print v_product(yr,yrn)/(vabs(yr)*vabs(yrn))
    twist = math.acos(v_product(yr,yrn)/(vabs(yr)*vabs(yrn))) #/3.15159*180 #Conversion from radians to degrees?
    return twist

def Steven_twist(Vec_1, Vec_2, Vec_3, Vec_4):    #original Steven's code
    x1, y1, z1 = Vec_1[0].astype(float)
    x2, y2, z2 = Vec_2[0].astype(float)
    x3, y3, z3 = Vec_3[0].astype(float)
    x4, y4, z4 = Vec_4[0].astype(float)

    # Calculate vector A, B, and M
    Ax, Ay, Az = x2 - x1, y2 - y1, z2 - z1
    Bx, By, Bz = x4 - x3, y4 - y3, z4 - z3
    Mx, My, Mz = (x4 + x3)/2 - (x1 + x2)/2, (y4 + y3)/2 - (y1 + y2)/2, (z4 + z3)/2 - (z1 + z2)/2

    # Calculate vector C and D
    Cx, Cy, Cz = Ay*Mz - Az*My, Az*Mx - Ax*Mz, Ax*My - Ay*Mx
    Dx, Dy, Dz = By*Mz - Bz*My, Bz*Mx - Bx*Mz, Bx*My - By*Mx

    # Calculate magnitudes of vectors C and D
    Cm = np.sqrt(Cx**2 + Cy**2 + Cz**2)
    Dm = np.sqrt(Dx**2 + Dy**2 + Dz**2)

    # Calculate dot product and theta
    dot = (Cx*Dx + Cy*Dy + Cz*Dz) / (Cm * Dm)
    theta = np.arccos(np.clip(dot, -1, 1))

    # Define Ex, Ey, and Ez
    Ex, Ey, Ez = Cy*Dz-Cz*Dy, Cz*Dx-Cx*Dz, Cx*Dy-Cy*Dx
    dot_sign = Ex*Mx+Ey*My+Ez*Mz

    # Calculate the angle according to the provided formula
    angle = theta - 2*(np.sign(dot_sign) < 0)*(theta - np.pi)
    
    return angle
    
def diff_twist(Vec_1, Vec_2, Vec_3, Vec_4):
    return Steven_twist(Vec_1, Vec_2, Vec_3, Vec_4) - Xun_twist(Vec_1, Vec_2, Vec_3, Vec_4)

def test_twist(Vec_1, Vec_2, Vec_3, Vec_4):     
    return Steven_twist(Vec_1, Vec_2, Vec_3, Vec_4)
    
if __name__=='__main__':
    test_twist()
