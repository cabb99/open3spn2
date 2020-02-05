[![Build Status](https://travis-ci.org/cabb99/open3spn2.svg?branch=master)](https://travis-ci.org/cabb99/open3spn2?branch=master)
[![Documentation Status](https://readthedocs.org/projects/open3spn2/badge/?version=latest)](https://open3spn2.readthedocs.io/en/latest/?badge=latest)
[![Anaconda-Server Badge](https://anaconda.org/wolynes-lab/open3spn2/badges/installer/conda.svg)](https://conda.anaconda.org/wolynes-lab)

# Open-3SPN2
A Implementation of the 3SPN.2 and 3SPN.2C coarse-grained molecular model of DNA in OpenMM.

3SPN.2 and 3SPN.2C are DNA coarse-grained forcefields developed by the [de Pablo Group](https://pme.uchicago.edu/de_pablo_lab/research/dna_folding_and_hybridization/3spn.2/). Each DNA nucleotide is modelled by 3 beads: one bead for the phosphate, a second one for the sugar and a third one nucleobase. This forcefield was adapted by the Wolynes group to model protein-DNA interactions as a complement for the [AWSEM](https://github.com/npschafer/openawsem) coarse-grained protein forcefield.

## Installation

Installation of the open3SPN2 repository is available through anaconda. Some of the dependencies (openmm, pdbfixer) are contained in the omnia channel.

```conda config --append channels omnia```

```conda install -c wolynes-lab open3spn2```

It is also necessary to install [X3DNA](http://x3dna.org/) and set the environment variable `$X3DNA` to the correct location. 
For protein-DNA simulations you will also need to install [openAWSEM](http://openawsem.org/) and add the openAWSEM path to the `$PYTHONPATH` environment variable.

## Documentation

Further documentation and tutorials are hosted in [readthedocs](https://open3spn2.readthedocs.io/en/latest/).

## Citations

If you publish any work using the open3SPN2 package, please include the following references:

3SPN.2
Hinckley, D. M., Freeman, G. S., Whitmer, J. K., & De Pablo, J. J. (2013). An experimentally-informed coarse-grained 3-site-per-nucleotide model of DNA: Structure, thermodynamics, and dynamics of hybridization. Journal of Chemical Physics, 139(14). https://doi.org/10.1063/1.4822042

3SPN.2C
Freeman, G. S., Hinckley, D. M., Lequieu, J. P., Whitmer, J. K., & De Pablo, J. J. (2014). Coarse-grained modeling of DNA curvature. Journal of Chemical Physics, 141(16). https://doi.org/10.1063/1.4897649
