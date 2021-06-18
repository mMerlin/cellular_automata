#!/usr/bin/env python
# coding=utf-8

"""
unit testing description
"""

import unittest
import os
import sys
import inspect

# remove need to add parent to path?
# `python -m unittest discover --start-directory tests --top-level-directory .`
# ???

# def add_parent_to_path():
#     """Add the parent of the testing folder to the import search path"""
#     test_path = os.path.split(inspect.getfile(inspect.currentframe()))[0]
#     parent_folder = os.path.realpath(os.path.abspath(test_path + '/..'))
#     # print('sys.path', sys.path)  # DEBUG
#     if parent_folder not in sys.path:
#         sys.path.insert(0, parent_folder)
#     # print('sys.path', sys.path)  # DEBUG


# add_parent_to_path()
# # pylint: disable=wrong-import-position
# import «module_from_parent_folder»
# # pylint: enable=wrong-import-position
# ../../unittesting/shell_test.py

class ConstructorArgumentsTestCase(unittest.TestCase):
    # @classmethod
    # def setUpClass(cls):
    # …
    # @classmethod
    # def tearDownClass(cls):
    # …
    def setUp(self):
        pass # run before each individual test
    def tearDown(self):
        pass
    # def test«descriptive_unit_test_name»(self):
    #     self.assert«True¦False¦Is«Not»None»(x,msg=None)
    #     self.assert««Not»«Equal¦In¦IsInstance»¦Is«Not»»(a, b, msg=None)
    #     self.assert«Raises¦Warns¦Logs»«Regex»()
    #     self.assert««Not»AlmostEqual»«Greater¦Less¦Count»«Equal»(a,b)
    #     self.assert«Not»Regex(s,r)
    #     self.assert«MultLine¦Sequence¦List¦Tuple¦Set¦Dict»Equal(a,b)
    # def test«descriptive_unit_test_name2»(self): …
    # @unittest.skip«If¦Unless»(«condition, » "message")
    # @unittest.expectedFailure
    # def test«descriptive_unit_test_name3»(self): …
    def test_neighbourhood_data_type(self):
        self.assertRaises(TypeError)
# @unittest.skip«If¦Unless»(«condition, » "message")
# class «TestingClassName2»(unittest.TestCase):


if __name__ == '__main__':
    unittest.main()
