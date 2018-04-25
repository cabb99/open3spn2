#/usr/bin/env python

"""
This file contains the base pair and base step parameters used in 3SPN.2C.
The parameters come from the following references:

Morozov, Alexandre V., et al. "Using DNA mechanics to predict in vitro 
nucleosome positions and formation energies." Nucleic acids research
37.14 (2009): 4707-4722.

Olson, Wilma K., et al. "DNA simulation benchmarks as revealed by X-ray
structures." Computational studies of RNA and DNA. Springer Netherlands, 2006.
235-257.
"""

class BasePair:
    def __init__(self,stea,steb,buckle, propeller, opening, shear, stretch, stagger):
        self.stea = stea
        self.steb = steb
        self.buckle = buckle
        self.propeller = propeller
        self.opening = opening
        self.shear = shear
        self.stretch = stretch
        self.stagger = stagger

class BaseStep:
    def __init__(self, stea, steb,twist, roll, tilt, shift, slide, rise):
        self.stea = stea
        self.steb = steb
        self.twist = twist
        self.roll = roll
        self.tilt = tilt
        self.shift = shift
        self.slide = slide
        self.rise = rise


# Defining the base pair parameters
bp = []
# A-T
bp.append(BasePair("A","T", 1.8, -15.0, 1.5, 0.07, -0.19, 0.07))
# T-A
bp.append(BasePair("T","A",-1.8, -15.0, 1.5, -0.07, -0.19, 0.07))
# C-G
bp.append(BasePair("C","G",-4.9, -8.7,-0.6, 0.16, -0.17, 0.15))
# G-C
bp.append(BasePair("G","C", 4.9, -8.7,-0.6, -0.16, -0.17, 0.15))

# Defining the base step parameters
step = []
# ApA
step.append(BaseStep("A","A",35.31,0.76,-1.84,-0.05,-0.21,3.27))
# ApT
step.append(BaseStep("A","T",31.21,-1.39, 0.0, 0.0, -0.56, 3.39))
# ApC
step.append(BaseStep("A","C",31.52,0.91,-0.64,0.21, -0.54,3.39))
# ApG
step.append(BaseStep("A","G",33.05, 3.15,-1.48,0.12,-0.27,3.38))
# TpA
step.append(BaseStep("T","A",36.2,5.25, 0.0, 0.0, 0.03,3.34))
# TpT
step.append(BaseStep("T","T",35.31,0.91,1.84,0.05,-0.21,3.27))
# TpC
step.append(BaseStep("T","C",34.80, 3.87,1.52, 0.27,-0.03,3.35))
# TpG
step.append(BaseStep("T","G",35.02, 5.95, 0.05,0.16,0.18, 3.38))
# CpA
step.append(BaseStep("C","A",35.02,5.95,-0.05,-0.16,0.18,3.38))
# CpT
step.append(BaseStep("C","T",33.05, 3.15,1.48,-0.12, -0.27,3.38))
# CpC
step.append(BaseStep("C","C",33.17, 3.86, 0.4, 0.02, -0.47,3.28))
# CpG
step.append(BaseStep("C","G",34.38,4.29,0.0,0.0,0.57,3.49))
# GpA
step.append(BaseStep("G","A",34.8,3.87,-1.52, -0.27,-0.03,3.35))
# GpT
step.append(BaseStep("G","T",31.52, 0.91, 0.64, -0.21, -0.54,3.39))
# GpC
step.append(BaseStep("G","C",34.38,0.67,0.0,0.0,-0.07,3.38))
# GpG
step.append(BaseStep("G","G",33.17, 3.86, -0.4, -0.02, -0.47,3.28))
