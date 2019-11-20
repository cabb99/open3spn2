Classes and functions
===============================

.. automodule:: ff3SPN2
    :noindex:
    
DNA
------------------
.. autoclass:: ff3SPN2.DNA
    :member-order: bysource
    :members:

System
------------------
.. autoclass:: ff3SPN2.System
    
    .. automethod:: __init__
    
    .. rubric:: Methods
    
DNA Forces
------------------
.. autoclass:: ff3SPN2.Bond
    :member-order: bysource
    :members:

.. autoclass:: ff3SPN2.Angle
    :member-order: bysource
    :members:

.. autoclass:: ff3SPN2.Dihedral
    :member-order: bysource
    :members:

.. autoclass:: ff3SPN2.Stacking
    :member-order: bysource
    :members:
    
.. autoclass:: ff3SPN2.BasePair
    :member-order: bysource
    :members:

.. autoclass:: ff3SPN2.CrossStacking
    :member-order: bysource
    :members:

.. autoclass:: ff3SPN2.Exclusion
    :member-order: bysource
    :members:

.. autoclass:: ff3SPN2.Electrostatics
    :member-order: bysource
    :members:

DNA - Protein Forces
--------------------
.. autoclass:: ff3SPN2.ExclusionProteinDNA
    :member-order: bysource
    :members:

.. autoclass:: ff3SPN2.ElectrostaticsProteinDNA
    :member-order: bysource
    :members:

Utils
-----
.. autofunction:: ff3SPN2.parseConfigTable
    
.. autofunction:: ff3SPN2.parsePDB
    
.. autofunction:: ff3SPN2.fixPDB
    
.. autofunction:: ff3SPN2.pdb2table

Tests
-----
.. autoclass:: ff3SPN2.TestEnergies
    :member-order: bysource
    :members:
    
.. autofunction:: ff3SPN2.test_DNA_from_xyz
    
.. autofunction:: ff3SPN2.test_parse_xyz
    
.. autofunction:: ff3SPN2.test_parse_log
    
Exceptions
----------
.. autoclass:: ff3SPN2.DNATypeError
