# USER-3SPN2
Coarse-grained molecular model of DNA (for LAMMPS)

The USER-3SPN2 directory contains the source code needed to run simulations of
the 3SPN.2 coarse-grained DNA model within the LAMMPS MD package. 3SPN.2 is a
highly-robust coarse-grained DNA model that is capable of studying processes
such as hybridization and DNA protein interactions.  For details about the
model itself and it's development, please see Hinckley et al., JCP, 139, 144903
,2013 (DOI: 10.1063/1.4822042)

The directory also contains sample input files for a normal simulation of a DNA
duplex, a restart, and a parallel tempering calculation. 

Also include is a directory containing a configuration generator that generates
B-DNA according to the crystal structure of Arnott et al.  This configuration
generator will generate data files that can be read in as data files using the
sample .in files.

Also include is a .pdf containing additional documentation and instructions for
compiling and visualization.

### Oct. 12, 2015 ###
The GCGI/ directory was added, which contains necessary files for using 
the general coarse- grained ions model presented by Hinckley and de Pablo, 
JCTC (2015) (DOI:10.1021/acs.jctc.5b00341).

If you have any questions and concerns associated with bugs (not with general 
(non-3SPN.2) compilation or LAMMPS issuyes!) please direct your inquiries to 
Dan Hinckley (hinckley.dan@gmail.com).

The following website has historically held information regarding on-going improvements to
the 3SPN.2 model, as well as recent publications by the de Pablo group and
others.

http://ime.uchicago.edu/de_pablo_lab/research/dna_folding_and_hybridization/3spn.2/

Moving forward, we hope to coordinate all updates and communication through github.
