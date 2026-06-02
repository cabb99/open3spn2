"""Tests for the lazy optional-forces registry."""

import subprocess
import sys

import pytest

import open3SPN2
from open3SPN2.optional_forces import optional_forces
import open3SPN2.force.bias as bias
import open3SPN2.force.collective_variables as cv

_EXPECTED = {
    "PositionRestraint", "StringProteinDNA", "ElectrostaticBias",
    "BasePairProteinHarmonicBias", "DistanceFromPoint", "StringLength",
    "BasePairProteinPosition",
}


def test_available_lists_all_without_importing():
    assert set(optional_forces.available()) == _EXPECTED
    assert optional_forces.available() == sorted(_EXPECTED)
    assert len(optional_forces) == len(_EXPECTED)
    assert set(optional_forces) == _EXPECTED


def test_getitem_resolves_to_the_real_class():
    assert optional_forces["PositionRestraint"] is bias.PositionRestraint
    assert optional_forces["StringProteinDNA"] is bias.StringProteinDNA
    assert optional_forces["DistanceFromPoint"] is cv.DistanceFromPoint
    assert optional_forces["BasePairProteinPosition"] is cv.BasePairProteinPosition


def test_attribute_access():
    assert optional_forces.ElectrostaticBias is bias.ElectrostaticBias
    assert optional_forces.StringLength is cv.StringLength


def test_exposed_on_package():
    assert open3SPN2.optional_forces is optional_forces


def test_unknown_name_raises():
    with pytest.raises(KeyError):
        optional_forces["NotAForce"]
    with pytest.raises(AttributeError):
        optional_forces.NotAForce


def test_core_import_does_not_pull_in_optional_modules():
    # Importing open3SPN2 must not import the bias / CV submodules; their import
    # cost (and any error in them) is deferred until a force is requested.
    code = (
        "import sys, open3SPN2; "
        "assert open3SPN2.optional_forces.available(); "
        "assert 'open3SPN2.force.bias' not in sys.modules, 'bias imported eagerly'; "
        "assert 'open3SPN2.force.collective_variables' not in sys.modules, 'cv imported eagerly'"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
