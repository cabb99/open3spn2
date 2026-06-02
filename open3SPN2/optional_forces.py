"""Opt-in bias and collective-variable forces, resolved lazily.

These forces are not part of the core force field and are never auto-added. They
are resolved on demand through the :data:`optional_forces` mapping so that an
import error in one of the optional modules cannot break ``import open3SPN2``:
the error surfaces only when that specific force is requested.

Examples
--------
>>> import open3SPN2
>>> open3SPN2.optional_forces.available()
['BasePairProteinHarmonicBias', 'BasePairProteinPosition', 'DistanceFromPoint', ...]
>>> PositionRestraint = open3SPN2.optional_forces['PositionRestraint']

Equivalently, the force can be imported straight from its submodule:

>>> from open3SPN2.force.bias import PositionRestraint
"""

import importlib

try:
    from collections.abc import Mapping
except ImportError:  # pragma: no cover - Python 2 fallback
    from collections import Mapping

# force name -> (submodule, attribute)
_REGISTRY = {
    # biases (apply energy to steer/restrain the system)
    "PositionRestraint": ("open3SPN2.force.bias", "PositionRestraint"),
    "StringProteinDNA": ("open3SPN2.force.bias", "StringProteinDNA"),
    "ElectrostaticBias": ("open3SPN2.force.bias", "ElectrostaticBias"),
    "BasePairProteinHarmonicBias": ("open3SPN2.force.bias", "BasePairProteinHarmonicBias"),
    # collective-variable readouts (identity energy)
    "DistanceFromPoint": ("open3SPN2.force.collective_variables", "DistanceFromPoint"),
    "StringLength": ("open3SPN2.force.collective_variables", "StringLength"),
    "BasePairProteinPosition": ("open3SPN2.force.collective_variables", "BasePairProteinPosition"),
}


class _LazyForceRegistry(Mapping):
    """Mapping of force name -> class, importing the submodule on first access."""

    def __init__(self, registry):
        self._registry = dict(registry)

    def __getitem__(self, name):
        if name not in self._registry:
            raise KeyError(name)
        module_name, attribute = self._registry[name]
        return getattr(importlib.import_module(module_name), attribute)

    def __getattr__(self, name):
        # Attribute-style access (optional_forces.PositionRestraint); never
        # intercept private/dunder lookups, which are not force names.
        if name.startswith('_'):
            raise AttributeError(name)
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    def __iter__(self):
        return iter(self._registry)

    def __len__(self):
        return len(self._registry)

    def available(self):
        """Names of the available optional forces (does not import anything)."""
        return sorted(self._registry)


optional_forces = _LazyForceRegistry(_REGISTRY)
