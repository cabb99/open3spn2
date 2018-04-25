#!/usr/bin/env python

"""
A script to determine the cross-terms for a LAMMPS simulation 
of 3SPN.2 with ions
"""
outfile = "cross_interactions.dat"

# Lists of charged and uncharged sites
charged = [15,16,17,18]
uncharged = range(2,15)
#uncharged.append(18)
#uncharged.append(19)
#uncharged.append(20)

# Corresponding sigmas

# These are taken from the all-atom force-field

# Original AMBER99
#charged_sigma = [3.3284,1.41225,4.40104]

# Joung--Cheatham AMBER99 Mg2+
charged_sigma = [2.493928,4.0,4.47766,3.3]
# From 3SPN.2
#uncharged_sigma = [4.5,6.2,5.4,7.1,4.9,6.4,5.4,7.1,4.9,6.4,5.4,7.1,4.9,6.4,8.0,2.0,8.0]
uncharged_sigma = [6.2,5.4,7.1,4.9,6.4,5.4,7.1,4.9,6.4,5.4,7.1,4.9,6.4]

# Specifying the epsilon
epsilon = 0.239000 # kcal/mol; 1 kJ/mol

# Cutoff factor for all interfactions (2^1/6)
cutoff = 1.12246204831

# pair style string
pairstyle = "lj/cut"

# Printing out all of the cross terms to file
#pair_coeff i j epsilon simga cutoff 1
file = open(outfile,"w")
for i in range(len(charged)):
    for j in range(len(uncharged)):
        if charged[i] >  uncharged[j]:
            myI = uncharged[j]
            myJ = charged[i]
        else:
            myI = charged[i]
            myJ = uncharged[j]
        #file.write("%s %d %d %s %3f %3f %f\n" % ("pair_coeff",charged[i],uncharged[j],pairstyle,epsilon,0.5 * (charged_sigma[i] + uncharged_sigma[j]), epsilon*cutoff))
        mySigma = 0.5 * (charged_sigma[i] + uncharged_sigma[j])
        file.write("%s %d %d %s %3f %3f %f\n" % ("pair_coeff",myI,myJ,pairstyle,epsilon,mySigma, cutoff*mySigma))
        #file.write("%s %d %d %s %3f %3f %f\n" % ("pair_coeff",myJ,myI,pairstyle,epsilon,mySigma, cutoff*mySigma))
        



