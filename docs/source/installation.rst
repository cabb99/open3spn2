Installation
=======================

Requirements
------------
X3DNA_ is needed on the 3SPN2.C forcefield in order to calculate the equilibrium bonds, angles and dihedrals. It is also needed to create a structure from sequence. Please follow their instructions for installation and add the installation directory to the X3DNA environment variable.

open3SPN also requires the following python libraries:

* **biopython**
* **pandas**
* **numpy**
* **scipy**
* **mdtraj**
* **openmm**
* **parmed**
* **pdbfixer**
* **nose**

From source code
----------------

The source code is available at https://github.com/cabb99/open3spn2/


From pip
--------

.. code-block:: bash

    $ pip install open3spn2

.. _X3DNA: https://x3dna.org/
