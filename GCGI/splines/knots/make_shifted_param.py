#!/usr/bin/env python
import os
types = ["P","Na","Mg","Cl","N"]
nPoints = 3
for i in range(len(types)):
    for j in range(i,len(types)):
        command = "python ../shift_pair_params.py ../orig_param/%s-%s.param %s-%s.param %d" % (types[i],types[j],types[i],types[j],nPoints)
        print command
        os.system(command)

bondType = ["B1","B2","A1"]
for i in range(len(bondType)):
    command = "python ../shift_bond_params.py ../orig_param/%s.param %s.param" % (bondType[i], bondType[i])
    print command
    os.system(command)
