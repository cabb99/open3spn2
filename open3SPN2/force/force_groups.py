"""Default OpenMM force-group numbers used across open3SPN2.

Every force accepts a ``force_group`` argument; the values here are only the
defaults. Centralizing them keeps the analysis energy decomposition in sync with
the force classes and documents which of the 32 slots are in use.

``TOTAL_ENERGY`` (groups 5-31) are summed to report the total potential energy.
Measurement-only readouts use ``MEASUREMENT`` (below 5) so they do not contribute
to that total. A ``CustomNonbondedForce`` carries one cutoff distance per force
group, so two non-bonded forces with different cutoffs must use different groups.
"""

# Measurement / analysis readouts (excluded from the energy total)
Q = 1
Rg = 2
MEASUREMENT = 4

# DNA (3SPN2) terms
BOND = 6
ANGLE = 7
STACKING = 8
DIHEDRAL = 9
BASE_PAIR = 10
CROSS_STACKING = 11
EXCLUSION_DNA = 12
ELECTROSTATICS_DNA = 13

# Protein-DNA terms
EXCLUSION_PROTEIN_DNA = 14
ELECTROSTATICS_PROTEIN_DNA = 15
PROTEIN_DNA_BIAS = 16

# AWSEM protein terms
BURIAL = 17
HELIX = 18
PULLING = 19
BACKBONE = 20
RAMA = 21
CONTACT = 22
FRAGMENT = 23
MEMBRANE = 24
ER = 25
TBM_Q = 26
BETA = 27
PAP = 28
HELICAL = 29
DEBYE_HUCKEL = 30
AMH_GO = 31

# Force groups summed to give the total potential energy.
TOTAL_ENERGY = range(5, 32)
