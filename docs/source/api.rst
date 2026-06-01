open3SPN2 Module Documentation
==============================

.. automodule:: open3SPN2
    :noindex:

Classes and Functions
---------------------

DNA Class
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: open3SPN2.DNA
    :member-order: bysource
    :members:

System Class
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: open3SPN2.System
    
    .. automethod:: __init__
    
    .. rubric:: Methods

DNA Forces
^^^^^^^^^^

.. autoclass:: open3SPN2.Bond
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.Angle
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.Dihedral
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.Stacking
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.BasePair
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.CrossStacking
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.Exclusion
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.Electrostatics
    :member-order: bysource
    :members:

DNA-Protein Forces
^^^^^^^^^^^^^^^^^^

.. autoclass:: open3SPN2.ExclusionProteinDNA
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.ElectrostaticsProteinDNA
    :member-order: bysource
    :members:

Utility Functions
^^^^^^^^^^^^^^^^^

.. autofunction:: open3SPN2.parseConfigTable

.. autofunction:: open3SPN2.parsePDB

.. autofunction:: open3SPN2.fixPDB

.. autofunction:: open3SPN2.pdb2table

Tests Functions
^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: tests.TestEnergies
    :member-order: bysource
    :members:

.. autofunction:: tests.test_DNA_from_xyz

.. autofunction:: tests.test_parse_xyz

.. autofunction:: tests.test_parse_log

Exceptions
^^^^^^^^^^

.. autoclass:: open3SPN2.DNATypeError