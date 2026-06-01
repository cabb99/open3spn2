[![Actions Status](https://github.com/cabb99/open3spn2/workflows/CI/badge.svg)](https://github.com/cabb99/open3spn2/actions)
[![Documentation Status](https://readthedocs.org/projects/open3spn2/badge/?version=latest)](https://open3spn2.readthedocs.io/en/latest/?badge=latest)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/open3spn2/badges/version.svg)](https://anaconda.org/conda-forge/open3spn2)

# Open-3SPN2
A Implementation of the 3SPN.2 and 3SPN.2C coarse-grained molecular model of DNA in OpenMM.

3SPN.2 and 3SPN.2C are DNA coarse-grained forcefields developed by the [de Pablo group](https://pme.uchicago.edu/group/de-pablo-group). Each DNA nucleotide is modelled by 3 beads: one bead for the phosphate, a second one for the sugar and a third one nucleobase. These forcefields were adapted by the [Wolynes group](https://wolynes.rice.edu/) to model protein-DNA interactions as a complement for the [AWSEM](https://github.com/npschafer/openawsem) coarse-grained protein forcefield.

## Installation

Installation of the open3SPN2 repository is available through anaconda.

```conda install -c conda-forge open3spn2```

For protein-DNA simulations you will also need to install [openAWSEM](https://github.com/npschafer/openawsem) and add the openAWSEM path to the `$PYTHONPATH` environment variable. In linux you can set the path variables on `~/.bashrc`.

```export PYTHONPATH=/path/to/openAWSEM:$PYTHONPATH```

## Documentation

Further documentation and tutorials are hosted in [readthedocs](https://open3spn2.readthedocs.io/en/latest/).

## Acknowledgment
Carlos Bueno was supported by the MolSSI Software Fellowship, under the mentorship of Jessica Nash. We thank AMD (Advanced Micro Devices, Inc.) for the donation of high-performance computing hardware and HPC resources. This project is also supported by the Center for Theoretical Biological Physics (NSF Grants PHY-2019745 and PHY-1522550), with additional support from the D.R. Bullard Welch Chair at Rice University (Grant No. C-0016 to PGW).

## Citations

If you publish any work using the open3SPN2 package, please include the following references:

Open3SPN2
Lu, W., Bueno, C., Schafer, N. P., Moller, J., Jin, S., Chen, X., ... & Wolynes, P. G. (2021). OpenAWSEM with Open3SPN2: A fast, flexible, and accessible framework for large-scale coarse-grained biomolecular simulations. PLoS computational biology, 17(2), e1008308. https://doi.org/10.1371/journal.pcbi.1008308

3SPN.2C
Freeman, G. S., Hinckley, D. M., Lequieu, J. P., Whitmer, J. K., & De Pablo, J. J. (2014). Coarse-grained modeling of DNA curvature. Journal of Chemical Physics, 141(16). https://doi.org/10.1063/1.4897649

3SPN.2
Hinckley, D. M., Freeman, G. S., Whitmer, J. K., & De Pablo, J. J. (2013). An experimentally-informed coarse-grained 3-site-per-nucleotide model of DNA: Structure, thermodynamics, and dynamics of hybridization. Journal of Chemical Physics, 139(14). https://doi.org/10.1063/1.4822042
