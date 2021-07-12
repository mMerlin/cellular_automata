#!/usr/bin/env python
# coding=utf-8

"""
verify typehints match expected objects

no method found to directly check object instances against non-trivial typehints.
isinstance(«object», «typehint») fails with
  `TypeError: isinstance() argument 2 cannot be a parameterized generic`
Some tools exist to parse typehint details, but that then relies on the parsing
code interpretation.
IE: is tuple() a valid match for tuple[int, ...]?
"""

# import unittest
# import pytest
# import os
# import sys
# import inspect
import automata_typehints as AHint

# remove need to add parent to path?
# `pipenv run python -m pytest «»`
# ? -q


# class ConstructorArgumentsTestCase(unittest.TestCase):
#     def setUp(self):
#         pass # run before each individual test
#     def tearDown(self):
#         pass
#     def test_neighbourhood_data_type(self):
#         self.assertRaises(TypeError)

def test_not_cell_address():
    assert not isinstance(None, AHint.CellAddressType)
    assert not isinstance(tuple(), AHint.CellAddressType)

def test_is_cell_address():
    assert isinstance((0,), AHint.CellAddressType)

# if __name__ == '__main__':
#     unittest.main()
