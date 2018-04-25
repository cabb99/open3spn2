#!/usr/bin/env  python

""" 
Force constants for bonded interactions, as described by Freeman et al. in [arXiv reference]
"""

# Angle force constants (kJ/mol/rad^2)

# Base-Sugar-Phosphate
BSP = {
"A": {
    "A":460,
    "T":370,
    "C":442,
    "G":358},
"T": {
    "A":120,
    "T":460,
    "C":383,
    "G":206},
"C": {
    "A":206,
    "T":358,
    "C":278,
    "G":278},
"G": {
    "A":383,
    "T":442,
    "C":336,
    "G":278}
}

# Phosphate-Sugar-Base
PSB = {
"A": {
    "A":460,
    "T":120,
    "C":206,
    "G":383},
"T": {
    "A":370,
    "T":460,
    "C":358,
    "G":442},
"C": {
    "A":442,
    "T":383,
    "C":278,
    "G":336},
"G": {
    "A":358,
    "T":206,
    "C":278,
    "G":278}
}

# Phosphate-Sugar-Phosphate
PSP = {
"A": {
    "A":300,
    "T":300,
    "C":300,
    "G":300},
"T": {
    "A":300,
    "T":300,
    "C":300,
    "G":300},
"C": {
    "A":300,
    "T":300,
    "C":300,
    "G":300},
"G": {
    "A":300,
    "T":300,
    "C":300,
    "G":300}
}

# Sugar-Phosphate-Sugar
SPS = {
"A": {
    "A":355,
    "T":147,
    "C":464,
    "G":368},
"T": {
    "A":230,
    "T":355,
    "C":442,
    "G":273},
"C": {
    "A":273,
    "T":368,
    "C":165,
    "G":478},
"G": {
    "A":442,
    "T":464,
    "C":228,
    "G":165}
}

# Bond force constant
kbond = 0.6

# Dihedral parameters
kgauss = 7.0 # kJ/mol/rad^2
kperiodic = 2.0 # kJ/mol/rad^2
sigma = 0.3

