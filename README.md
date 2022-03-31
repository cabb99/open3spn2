[![Build Status](https://travis-ci.org/cabb99/open3spn2.svg?branch=master)](https://travis-ci.org/cabb99/open3spn2?branch=master)
[![Documentation Status](https://readthedocs.org/projects/open3spn2/badge/?version=latest)](https://open3spn2.readthedocs.io/en/latest/?badge=latest)
[![Anaconda-Server Badge](https://anaconda.org/wolynes-lab/open3spn2/badges/installer/conda.svg)](https://conda.anaconda.org/wolynes-lab)

# Open-3SPN2
A Implementation of the 3SPN.2 and 3SPN.2C coarse-grained molecular model of DNA in OpenMM.

3SPN.2 and 3SPN.2C are DNA coarse-grained forcefields developed by the [de Pablo group](https://pme.uchicago.edu/group/de-pablo-group). Each DNA nucleotide is modelled by 3 beads: one bead for the phosphate, a second one for the sugar and a third one nucleobase. These forcefields were adapted by the [Wolynes group](https://wolynes.rice.edu/) to model protein-DNA interactions as a complement for the [AWSEM](https://github.com/npschafer/openawsem) coarse-grained protein forcefield.

## Installation

Installation of the open3SPN2 repository is available through anaconda. Some of the dependencies  are contained in the omnia channel (openmm, pdbfixer) or in the conda-forge channel (mdtraj).

```conda config --append channels omnia```

```conda config --append channels conda-forge```

```conda install -c wolynes-lab open3spn2```

It is also necessary to install [X3DNA](http://x3dna.org/) >= 2.4 and set the environment variable `$X3DNA` to the location of the installation folder.For protein-DNA simulations you will also need to install [openAWSEM](https://github.com/npschafer/openawsem) and add the openAWSEM path to the `$PYTHONPATH` environment variable. In linux you can set the path variables on `~/.bashrc`.

```export X3DNA=/path/to/x3dna-v2.4```

```export PYTHONPATH=/path/to/openAWSEM:$PYTHONPATH```

Note: open3SPN2 requires the installation of openMM. As of Feb 10 2020, openMM requires 3.5 <= python <= 3.7 and 7.5 <= CUDA <= 10.1.

## Documentation

Further documentation and tutorials are hosted in [readthedocs](https://open3spn2.readthedocs.io/en/latest/).

## Citations

If you publish any work using the open3SPN2 package, please include the following references:

Open3SPN2
Lu, W., Bueno, C., Schafer, N. P., Moller, J., Jin, S., Chen, X., ... & Wolynes, P. G. (2021). OpenAWSEM with Open3SPN2: A fast, flexible, and accessible framework for large-scale coarse-grained biomolecular simulations. PLoS computational biology, 17(2), e1008308. https://doi.org/10.1371/journal.pcbi.1008308

3SPN.2
Hinckley, D. M., Freeman, G. S., Whitmer, J. K., & De Pablo, J. J. (2013). An experimentally-informed coarse-grained 3-site-per-nucleotide model of DNA: Structure, thermodynamics, and dynamics of hybridization. Journal of Chemical Physics, 139(14). https://doi.org/10.1063/1.4822042

3SPN.2C
Freeman, G. S., Hinckley, D. M., Lequieu, J. P., Whitmer, J. K., & De Pablo, J. J. (2014). Coarse-grained modeling of DNA curvature. Journal of Chemical Physics, 141(16). https://doi.org/10.1063/1.4897649
