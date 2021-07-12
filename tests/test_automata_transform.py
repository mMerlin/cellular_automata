#!/usr/bin/env python
# coding=utf-8
# pylint: disable=W0212,R0914

"""
regression tests for cellular automata universe transform management
"""

# import sys
import re
from typing import Hashable
# from collections import namedtuple
# from _pytest._code.code import ExceptionInfo
import pytest
# reduce duplication
from test_create_universe import (base_universe_instance_1d,
    base_universe_instance_2d,
    base_universe_instance_3d,
)
from common_test_data import (
    IDENTITY_MATRIX_1D,
    IDENTITY_MATRIX_2D,
    IDENTITY_MATRIX_3D,
    REFLECT_MATRIX_1D,
    ZERO_MATRIX_1D,
    ZERO_MATRIX_2D,
    ZERO_MATRIX_3D,
    REFLECTION_MATRIX_2D_HORIZONTAL,
    REFLECTION_MATRIX_2D_VERTICAL,
    ROTATE_MATRIX_2D_90,
    ROTATE_MATRIX_2D_180,
    ROTATE_MATRIX_2D_270,
    TRANPOSE_MATRIX_2D,
    REFLECT_AND_ROTATE_2D,
    BAD_TYPE_NOT_ITERABLE,
    BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE,
    NOT_INTEGER_TYPE_SAMPLES,
)
from automata_transforms import (AutomataTransforms, TransformSequence, XformSequence,
    _base_index, _check_transform_key
)
# avoid need to add parent directory to path
# `pipenv run python -m pytest «»`
# ? -q -v

BASE_INDEX_GOOD_KEYS = (
    None,
    "test",
    ("try", 2),
    frozenset(["this"]),
)
BASE_INDEX_BAD_KEYS = (
    set(["this"]),
    [1, [2], 3],
    {'1': True},
)
XFORM_GOOD_SEQUENCE = (
    0,
    1,
    99,
    12345,
)
# 'normal' sequence will be positive integer, but zero and negative are «currently»
# technically valid, and might be used for marker cases
XFORM_BAD_SEQUENCE = (
    None,
    "test",
    tuple(),
    list(),
    1.5,
)

# The only good 1D transform matrix is reflection. There is no distinction between
# rotation and reflection. Other than 0 and identity, all other 1D transformation
# matrices are scaling, with possible reflection
GOOD_MATRICES_1D = (
    REFLECT_MATRIX_1D,
)
GOOD_MATRICES_2D = (
    REFLECTION_MATRIX_2D_HORIZONTAL,
    ROTATE_MATRIX_2D_90,
    ROTATE_MATRIX_2D_180,
    ROTATE_MATRIX_2D_270,
    REFLECTION_MATRIX_2D_VERTICAL,
)
NOT_MATRICES_TYPE_1D = (
    *NOT_INTEGER_TYPE_SAMPLES,
    (1,2),
    (('a',)),
)
NOT_MATRICES_TYPE_2D = NOT_MATRICES_TYPE_1D
NOT_MATRICES_1D = (
    IDENTITY_MATRIX_2D,
    IDENTITY_MATRIX_3D,
    ZERO_MATRIX_2D,
    ZERO_MATRIX_3D,
    ((1,2),(2,3)),
    ((1,2,3),(2,3,4),(4,5,6)),
)
NOT_MATRICES_2D = (
    IDENTITY_MATRIX_1D,
    IDENTITY_MATRIX_3D,
    ZERO_MATRIX_1D,
    ZERO_MATRIX_3D,
    ((7,),),
    ((7,),(3,)),
    ((1,2,3),),
    ((1,2,3),(2,3,4)),
    ((1,2,3),(2,3,4),(4,5,6)),
)
BAD_TRANSFORM_1D = (
    ((3,),),
    # IDENTITY_MATRIX_1D, # identity matrix passed initial neighbourhood -> neighbourhood check
    ZERO_MATRIX_1D,
    ((-7,),),
)
BAD_TRANSFORM_2D = (
    ((3,0),(0,1)),
    ((-3,0),(0,1)),
    ZERO_MATRIX_2D,
    ((0,-7,),(1,0),),
    ((0,7,),(1,0),),
    ((0,7,),(1,0),),
    ((0,7,),(0,0),),
)
BAD_MATRIX_VECTOR_COUNT_1D = (
    ((1,),(2,)),
)
BAD_MATRIX_VECTOR_COUNT_2D = (
    ((1,2),),
    ((1,2),(2,5),(3,7)),
)

# utility functions to use in tests

def not_matrix_typeerror_1d() -> tuple[object]:
    """combine test cases for invalid 1d matrices that raise TypeError"""
    assert len(BAD_TYPE_NOT_ITERABLE) > 0
    assert len(BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE) > 0
    assert len(BAD_MATRIX_VECTOR_COUNT_1D) > 0
    assert len(NOT_MATRICES_TYPE_1D) > 0
    return (
        *BAD_TYPE_NOT_ITERABLE,
        *(case for (case, _err_details) in BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE),
        *(case for case in BAD_MATRIX_VECTOR_COUNT_1D),
        *(case for case in NOT_MATRICES_TYPE_1D),
    )
def not_matrix_typeerror_2d() -> tuple[object]:
    """combine test cases for invalid 2d matrices that raise TypeError"""
    assert len(BAD_TYPE_NOT_ITERABLE) > 0
    assert len(BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE) > 0
    assert len(BAD_MATRIX_VECTOR_COUNT_2D) > 0
    assert len(NOT_MATRICES_TYPE_2D) > 0
    return (
        *BAD_TYPE_NOT_ITERABLE,
        *(case for (case, _err_details) in BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE),
        *(case for case in BAD_MATRIX_VECTOR_COUNT_2D),
        *(case for case in NOT_MATRICES_TYPE_2D),
    )
def bad_matrix_valueerror_1d() -> tuple[tuple[int]]:
    """combine test cases for invalid 1d matrices that raise ValueError"""
    assert len(NOT_MATRICES_1D) > 0
    assert len(BAD_TRANSFORM_1D) > 0
    return (
        *NOT_MATRICES_1D,
        *BAD_TRANSFORM_1D,
        IDENTITY_MATRIX_1D,
    )
def bad_matrix_valueerror_2d() -> tuple[tuple[int, int]]:
    """combine test cases for invalid 2d matrices that raise ValueError"""
    assert len(NOT_MATRICES_2D) > 0
    assert len(BAD_TRANSFORM_2D) > 0
    return (
        *NOT_MATRICES_2D,
        *BAD_TRANSFORM_2D,
        IDENTITY_MATRIX_2D,
    )
def match_added_transform_reflect_1d(instance: AutomataTransforms, key: Hashable) -> None:
    """verify instance properties with (only) 1d reflection transform added"""
    test_matrix = REFLECT_MATRIX_1D
    assert instance._transform_cycles == {None: 0, key: 2}
    assert _base_index(key) == (key, 1)
    assert instance._index_to_transform == {_base_index(key): test_matrix}
    assert instance._transform_to_index == {test_matrix: _base_index(key)}
    assert instance._scratch.primes == {
        instance._prime_cells: TransformSequence(None, 0),
        frozenset((-cell[0],) for cell in instance._prime_cells): _base_index(key)}
def match_added_transform_reflect_horizontal_2d(instance: AutomataTransforms,
        key: Hashable) -> None:
    """verify instance properties with (only) 2d reflection transform added"""
    test_matrix = REFLECTION_MATRIX_2D_HORIZONTAL
    assert instance._transform_cycles == {None: 0, key: 2}
    assert _base_index(key) == (key, 1)
    assert instance._index_to_transform == {_base_index(key): test_matrix}
    assert instance._transform_to_index == {test_matrix: _base_index(key)}
    assert instance._scratch.primes == {
        instance._prime_cells: TransformSequence(None, 0),
        frozenset((-cell[0],cell[1]) for cell in instance._prime_cells): _base_index(key)}
def match_added_transform_rotate_2d(instance: AutomataTransforms, key: Hashable) -> None:
    """verify instance properties with (only) 2d rotation transform added"""
    # assert len(instance._transform_cycles) == 2
    # assert key in instance._transform_cycles
    # assert instance._transform_cycles[key] == 4
    assert instance._transform_cycles == {None: 0, key: 4} # 4 90° rotation steps in cycle
    assert _base_index(key) == (key, 1)
    assert frozenset(instance._index_to_transform.keys()) == frozenset(
        ((key, 1), (key, 2), (key, 3)))
    assert len(instance._index_to_transform) == 3 # 3 transforms in cycle
    assert instance._index_to_transform[(key, 1)] == ROTATE_MATRIX_2D_90
    assert instance._index_to_transform[(key, 2)] == ROTATE_MATRIX_2D_180
    assert instance._index_to_transform[(key, 3)] == ROTATE_MATRIX_2D_270
    assert frozenset(instance._index_to_transform.keys()) == frozenset(
        instance._transform_to_index.values())
    assert frozenset(instance._index_to_transform.values()) == frozenset(
        instance._transform_to_index.keys())
    assert frozenset(instance._scratch.primes.values()) == frozenset(
        ((None, 0), *instance._index_to_transform.keys()))
    assert instance._transform_to_index[ROTATE_MATRIX_2D_90] == (key, 1)
    assert instance._transform_to_index[ROTATE_MATRIX_2D_180] == (key, 2)
    assert instance._transform_to_index[ROTATE_MATRIX_2D_270] == (key, 3)

# checks for the standalone (non-class) functions in the module

def test_base_index() -> None:
    """verify result for base index creation"""
    assert len(BASE_INDEX_GOOD_KEYS) > 0
    for test_key in BASE_INDEX_GOOD_KEYS:
        idx = _base_index(test_key)
        assert idx == TransformSequence(test_key, 1)
        assert isinstance(idx, TransformSequence)
    assert len(BASE_INDEX_BAD_KEYS) > 0
    for test_key in BASE_INDEX_BAD_KEYS:
        with pytest.raises(TypeError, match=re.compile(
                "unhashable type: '{}'".format(type(test_key).__name__))):
            idx = _base_index(test_key)
def test_check_transform_key() -> None:
    """verify exception for (only) non hashable key values"""
    assert len(BASE_INDEX_GOOD_KEYS) > 0
    for test_key in BASE_INDEX_GOOD_KEYS:
        _check_transform_key(test_key)
    assert len(BASE_INDEX_BAD_KEYS) > 0
    for test_key in BASE_INDEX_BAD_KEYS:
        with pytest.raises(TypeError, match=re.compile(
                "unhashable type: '{}'".format(type(test_key).__name__))):
            _check_transform_key(test_key)

# checks for the mini TransformSequence class, which is really a namedtuple that
# verifies the datatype of its elements on (actually before) instantiation.

def test_transform_sequence() -> None:
    """check good and bad values for creating TransformSequence keys"""
    assert len(BASE_INDEX_GOOD_KEYS) > 0
    assert len(XFORM_GOOD_SEQUENCE) > 0
    for good_key in BASE_INDEX_GOOD_KEYS:
        for good_seq in XFORM_GOOD_SEQUENCE:
            instance = TransformSequence(good_key, good_seq)
            assert isinstance(instance, tuple)
            # assert isinstance(instance, namedtuple)
            assert isinstance(instance, XformSequence)
            assert instance.key == good_key
            assert instance.seq == good_seq
            assert hash(instance) == hash((good_key, good_seq))
    assert len(BASE_INDEX_BAD_KEYS) > 0
    good_seq = XFORM_GOOD_SEQUENCE[0]
    for bad_key in BASE_INDEX_BAD_KEYS:
        with pytest.raises(TypeError, match=re.compile(
                "unhashable type: '{}'".format(type(bad_key).__name__))):
            instance = TransformSequence(bad_key, good_seq)
    assert len(XFORM_BAD_SEQUENCE) > 0
    good_key = BASE_INDEX_GOOD_KEYS[0]
    for bad_seq in XFORM_BAD_SEQUENCE:
        with pytest.raises(TypeError, match=re.compile(
                "not integer type: '{}'".format(type(bad_seq).__name__))):
            instance = TransformSequence(good_key, bad_seq)

# checks for the AutomataTransforms class. The main class for the module

def test_create_transforms_1d() -> None:
    """verify creation and basic properties of a 1d transform management instance"""
    xform = AutomataTransforms(base_universe_instance_1d())
    assert isinstance(xform, AutomataTransforms)
    assert xform._universe.dimensions == 1
    assert isinstance(hash(xform), int)
def test_create_transforms_2d() -> None:
    """verify creation and basic properties of a 2d transform management instance"""
    xform = AutomataTransforms(base_universe_instance_2d())
    assert isinstance(xform, AutomataTransforms)
    assert xform._universe.dimensions == 2
    assert isinstance(hash(xform), int)
def test_create_transforms_3d() -> None:
    """verify creation and basic properties of a 3d transform management instance"""
    xform = AutomataTransforms(base_universe_instance_3d())
    assert isinstance(xform, AutomataTransforms)
    assert xform._universe.dimensions == 3
    assert isinstance(hash(xform), int)

def test_check_transform_not_matrix_1d() -> None:
    """verify raised exception for 1d non-matrix arguments"""
    xform = AutomataTransforms(base_universe_instance_1d())
    not_matrix = not_matrix_typeerror_1d()
    for case in not_matrix:
        # minimal check: details verified in tests for validate_matrix, valid_address
        with pytest.raises(TypeError):
            xform._check_transform_matrix(case)
    assert len(NOT_MATRICES_1D) > 0
    for case in NOT_MATRICES_1D:
        # minimal check: details verified in tests for validate_matrix, valid_address
        with pytest.raises(ValueError):
            xform._check_transform_matrix(case)
def test_check_transform_not_matrix_2d() -> None:
    """verify raised exception for 2d non-matrix arguments"""
    xform = AutomataTransforms(base_universe_instance_2d())
    not_matrix = not_matrix_typeerror_2d()
    for case in not_matrix:
        # minimal check: details verified in tests for validate_matrix, valid_address
        with pytest.raises(TypeError):
            xform._check_transform_matrix(case)
    assert len(NOT_MATRICES_2D) > 0
    for case in NOT_MATRICES_2D:
        # minimal check: details verified in tests for validate_matrix, valid_address
        with pytest.raises(ValueError):
            xform._check_transform_matrix(case)
def test_check_transform_matrix_not_transform_1d() -> None:
    """verify raised exception for 1d matrix that is not a transform"""
    xform = AutomataTransforms(base_universe_instance_1d())
    assert len(BAD_TRANSFORM_1D) > 0
    for matrix in BAD_TRANSFORM_1D:
        assert isinstance(matrix[0][0], int)
        expected_message = "transformation of neighbourhood with {} is not the " \
            "neighbourhood".format(matrix)
        with pytest.raises(ValueError) as excinfo:
            xform._check_transform_matrix(matrix)
        assert excinfo.value.args == (
            frozenset({(matrix[0][0],), (-matrix[0][0],)}), expected_message)
def test_check_transform_matrix_not_transform_2d() -> None:
    """verify raised exception for 2d matrix that is not a transform"""
    xform = AutomataTransforms(base_universe_instance_2d())
    assert len(BAD_TRANSFORM_2D) > 0
    for matrix in BAD_TRANSFORM_2D:
        assert isinstance(matrix[0][0], int)
        expected_message = "transformation of neighbourhood with {} is not the " \
            "neighbourhood".format(matrix)
        with pytest.raises(ValueError) as excinfo:
            xform._check_transform_matrix(matrix)
        # can not reproduce the transformed coordinate set without (effectively) using
        # the same code as the method being tested. Use more indirect sanity checks
        assert isinstance(excinfo.value.args, tuple)
        assert len(excinfo.value.args) == 2
        assert isinstance(excinfo.value.args[0], frozenset)
        assert isinstance(excinfo.value.args[1], str)
        # assert len(excinfo.value.args[0]) in (1, xform._universe.neighbourhood_population)
        assert 1 <= len(excinfo.value.args[0]) <= xform._universe.neighbourhood_population
        # length matches neighbourhood_population ONLY if the matrix does not produce
        # duplicates from the neighbourhood (IE multiplication by zero)
        assert excinfo.value.args[1] == expected_message
def test_check_transform_identity_matrix_1d() -> None:
    """verify raised exception for 1d identity matrix"""
    xform = AutomataTransforms(base_universe_instance_1d())
    expected_message = "transform {} matches the identity matrix".format(IDENTITY_MATRIX_1D)
    # expanded message contains regex pattern, causing match to fail
    # with pytest.raises(ValueError, match=re.compile(expected_message)):
    with pytest.raises(ValueError) as excinfo:
        xform._check_transform_matrix(IDENTITY_MATRIX_1D)
    assert excinfo.value.args[0] == expected_message
def test_check_transform_identity_matrix_2d() -> None:
    """verify raised exception for 2d identity matrix"""
    xform = AutomataTransforms(base_universe_instance_2d())
    expected_message = "transform {} matches the identity matrix".format(IDENTITY_MATRIX_2D)
    # expanded message contains regex pattern, causing match to fail
    # with pytest.raises(ValueError, match=re.compile(expected_message)):
    with pytest.raises(ValueError) as excinfo:
        xform._check_transform_matrix(IDENTITY_MATRIX_2D)
    assert excinfo.value.args[0] == expected_message

def test_make_transform_hashable_bad_transform_1d() -> None:
    """verify exception raised for bad 1d transformation matrix"""
    xform = AutomataTransforms(base_universe_instance_1d())
    not_matrix = not_matrix_typeerror_1d()
    for case in not_matrix:
        # minimal check: details verified in tests for check_transform_matrix (and others)
        with pytest.raises(TypeError):
            xform._make_transform_hashable(case)
    bad_matrix = bad_matrix_valueerror_1d()
    for case in bad_matrix:
        # minimal check: details verified in tests for check_transform_matrix (and others)
        with pytest.raises(ValueError):
            xform._make_transform_hashable(case)
def test_make_transform_hashable_bad_transform_2d() -> None:
    """verify exception raised for bad 2d transformation matrix"""
    xform = AutomataTransforms(base_universe_instance_2d())
    not_matrix = not_matrix_typeerror_2d()
    for case in not_matrix:
        # minimal check: details verified in tests for check_transform_matrix (and others)
        with pytest.raises(TypeError):
            xform._make_transform_hashable(case)
    bad_matrix = bad_matrix_valueerror_2d()
    for case in bad_matrix:
        # minimal check: details verified in tests for check_transform_matrix (and others)
        with pytest.raises(ValueError):
            xform._make_transform_hashable(case)
# def test_make_transform_hashable_not_hashable_1d() -> None:
    # """verify exception raised for bad transformation matrix"""
    # # Unable to create a test for this case. Any matrix that matches this condition also
    # # matches earlier test cases
    # xform = AutomataTransforms(base_universe_instance_1d())
def test_make_transform_hashable_results_1d() -> None:
    """verify correct 1d hashable transform results"""
    xform = AutomataTransforms(base_universe_instance_1d())
    assert len(GOOD_MATRICES_1D) > 0
    for case in GOOD_MATRICES_1D:
        hashable = xform._make_transform_hashable(case)
        assert len(hashable) == xform._universe.dimensions
        assert isinstance(hashable, tuple)
        assert isinstance(hash(hashable), int)
def test_make_transform_hashable_results_2d() -> None:
    """verify correct 2d hashable transform results"""
    xform = AutomataTransforms(base_universe_instance_2d())
    assert len(GOOD_MATRICES_2D) > 0
    for case in GOOD_MATRICES_2D:
        hashable = xform._make_transform_hashable(case)
        assert isinstance(hashable, tuple)
        assert len(hashable) == xform._universe.dimensions
        assert isinstance(hash(hashable), int)

def test_add_transform_cycle_bad_arguments_1d() -> None:
    """verify exception raised for bad 1d key or transform arguments"""
    xform = AutomataTransforms(base_universe_instance_1d())
    assert len(BASE_INDEX_BAD_KEYS) > 0
    test_matrix = REFLECT_MATRIX_1D
    for test_key in BASE_INDEX_BAD_KEYS:
        # minimal check: details verified in tests for check_transform_key
        with pytest.raises(TypeError):
            xform.add_transform_cycle(test_key, test_matrix)
    not_matrix = not_matrix_typeerror_1d()
    for case in not_matrix:
        # minimal check: details verified in tests for make_transform_hashable (and others)
        with pytest.raises(TypeError):
            xform.add_transform_cycle("test", case)
    bad_matrix = bad_matrix_valueerror_1d()
    for case in bad_matrix:
        # minimal check: details verified in tests for make_transform_hashable (and others)
        with pytest.raises(ValueError):
            xform.add_transform_cycle("test", case)
def test_add_transform_cycle_bad_arguments_2d() -> None:
    """verify exception raised for bad 2d key or transform arguments"""
    xform = AutomataTransforms(base_universe_instance_2d())
    assert len(BASE_INDEX_BAD_KEYS) > 0
    test_matrix = REFLECTION_MATRIX_2D_HORIZONTAL
    for test_key in BASE_INDEX_BAD_KEYS:
        # minimal check: details verified in tests for check_transform_key
        with pytest.raises(TypeError):
            xform.add_transform_cycle(test_key, test_matrix)
    not_matrix = not_matrix_typeerror_2d()
    for case in not_matrix:
        # minimal check: details verified in tests for make_transform_hashable (and others)
        with pytest.raises(TypeError):
            xform.add_transform_cycle("test", case)
    bad_matrix = bad_matrix_valueerror_2d()
    for case in bad_matrix:
        # minimal check: details verified in tests for make_transform_hashable (and others)
        with pytest.raises(ValueError):
            xform.add_transform_cycle("test", case)
def test_add_transform_cycle_duplicate_key_1d() -> None:
    """verify raised exceptions when matching 1d data already stored"""
    xform = AutomataTransforms(base_universe_instance_1d())
    test_key = 'test01'
    test_matrix = REFLECT_MATRIX_1D
    xform._transform_cycles[test_key] = None # setup to cause duplicate error
    assert test_key in xform._transform_cycles
    expected_message = "Cycle identifier '{}' already in use".format(test_key)
    with pytest.raises(ValueError, match=re.compile(expected_message)):
        xform.add_transform_cycle(test_key, test_matrix)
    xform._transform_cycles.pop(test_key)
    assert test_key not in xform._transform_cycles

    test_index = _base_index(test_key)
    xform._index_to_transform[test_index] = None # setup to cause duplicate error
    assert test_index in xform._index_to_transform
    expected_message = "Base cycle index {} already in use".format(test_index)
    with pytest.raises(ValueError) as excinfo:
        xform.add_transform_cycle(test_key, test_matrix)
    assert excinfo.value.args == (expected_message,)
    xform._index_to_transform.pop(test_index)
    assert test_index not in xform._index_to_transform

    xform._transform_to_index[test_matrix] = test_index
    assert test_matrix in xform._transform_to_index
    expected_message = "'{}' transform matches {}: {}".format(
        test_key, test_index, test_matrix)
    # with pytest.raises(ValueError, match=re.compile(expected_message)):
    with pytest.raises(ValueError) as excinfo:
        xform.add_transform_cycle(test_key, test_matrix)
    assert excinfo.value.args == (expected_message,)
    xform._transform_to_index.pop(test_matrix)
    assert test_matrix not in xform._transform_to_index
def test_add_transform_cycle_duplicate_key_2d() -> None:
    """verify raised exceptions when matching 2d data already stored"""
    xform = AutomataTransforms(base_universe_instance_2d())
    test_key = 'test01'
    test_matrix = REFLECTION_MATRIX_2D_HORIZONTAL
    xform._transform_cycles[test_key] = None # setup to cause duplicate error
    assert test_key in xform._transform_cycles
    expected_message = "Cycle identifier '{}' already in use".format(test_key)
    with pytest.raises(ValueError, match=re.compile(expected_message)):
        xform.add_transform_cycle(test_key, test_matrix)
    xform._transform_cycles.pop(test_key)
    assert test_key not in xform._transform_cycles

    test_index = _base_index(test_key)
    xform._index_to_transform[test_index] = None # setup to cause duplicate error
    assert test_index in xform._index_to_transform
    expected_message = "Base cycle index {} already in use".format(test_index)
    with pytest.raises(ValueError) as excinfo:
        xform.add_transform_cycle(test_key, test_matrix)
    assert excinfo.value.args == (expected_message,)
    xform._index_to_transform.pop(test_index)
    assert test_index not in xform._index_to_transform

    xform._transform_to_index[test_matrix] = test_index
    assert test_matrix in xform._transform_to_index
    expected_message = "'{}' transform matches {}: {}".format(
        test_key, test_index, test_matrix)
    # with pytest.raises(ValueError, match=re.compile(expected_message)):
    with pytest.raises(ValueError) as excinfo:
        xform.add_transform_cycle(test_key, test_matrix)
    assert excinfo.value.args == (expected_message,)
    xform._transform_to_index.pop(test_matrix)
    assert test_matrix not in xform._transform_to_index
def test_add_transfrom_cycle_reflection_result_1d() -> None:
    """verify stored results after adding 1d reflection transform"""
    xform = AutomataTransforms(base_universe_instance_1d())
    test_key = 'reflection'
    test_matrix = REFLECT_MATRIX_1D
    xform.add_transform_cycle(test_key, test_matrix)
    match_added_transform_reflect_1d(xform, test_key)

    expected_message = "Cycle identifier '{}' already in use".format(test_key)
    with pytest.raises(ValueError, match=re.compile(expected_message)):
        xform.add_transform_cycle(test_key, test_matrix)
    match_added_transform_reflect_1d(xform, test_key)

    test_key2 = 'other'
    test_index = _base_index(test_key)
    expected_message = "'{}' transform matches {}: {}".format(
        test_key2, test_index, test_matrix)
    with pytest.raises(ValueError) as excinfo:
        xform.add_transform_cycle(test_key2, test_matrix)
    assert excinfo.value.args == (expected_message,)
    match_added_transform_reflect_1d(xform, test_key)
def test_add_transfrom_cycle_reflection_result_2d() -> None:
    """verify stored results after adding 2d reflection transform"""
    xform = AutomataTransforms(base_universe_instance_2d())
    test_key = 'reflection'
    test_matrix = REFLECTION_MATRIX_2D_HORIZONTAL
    xform.add_transform_cycle(test_key, test_matrix)
    match_added_transform_reflect_horizontal_2d(xform, test_key)

    expected_message = "Cycle identifier '{}' already in use".format(test_key)
    with pytest.raises(ValueError, match=re.compile(expected_message)):
        xform.add_transform_cycle(test_key, test_matrix)
    match_added_transform_reflect_horizontal_2d(xform, test_key)

    test_key2 = 'other'
    test_index = _base_index(test_key)
    expected_message = "'{}' transform matches {}: {}".format(
        test_key2, test_index, test_matrix)
    with pytest.raises(ValueError) as excinfo:
        xform.add_transform_cycle(test_key2, test_matrix)
    assert excinfo.value.args == (expected_message,)
    match_added_transform_reflect_horizontal_2d(xform, test_key)

def test_generate_combinations_not_possible_1d() -> None:
    """verify exceptions attempting to search combinations for 1d transforms"""
    xform = AutomataTransforms(base_universe_instance_1d())
    with pytest.raises(ValueError, match=re.compile(
            "No transforms have been added yet. Nothing to base generated transforms on")):
        xform.generate_combination_transforms()
    test_key = 'reflection'
    test_matrix = REFLECT_MATRIX_1D
    xform.add_transform_cycle(test_key, test_matrix)
    match_added_transform_reflect_1d(xform, test_key)
    expected_message = "Only a single transform cycle '{}' exists. Nothing to base" \
        " combinations on".format(test_key)
    with pytest.raises(ValueError, match=re.compile(expected_message)):
        xform.generate_combination_transforms()
    match_added_transform_reflect_1d(xform, test_key)
def test_generate_combinations_not_possible_2d() -> None:
    """verify exceptions attempting to search combinations for 2d transforms"""
    xform = AutomataTransforms(base_universe_instance_2d())
    with pytest.raises(ValueError, match=re.compile(
            "No transforms have been added yet. Nothing to base generated transforms on")):
        xform.generate_combination_transforms()
    test_key = 'reflection'
    test_matrix = REFLECTION_MATRIX_2D_HORIZONTAL
    xform.add_transform_cycle(test_key, test_matrix)
    match_added_transform_reflect_horizontal_2d(xform, test_key)
    expected_message = "Only a single transform cycle '{}' exists. Nothing to base" \
        " combinations on".format(test_key)
    with pytest.raises(ValueError, match=re.compile(expected_message)):
        xform.generate_combination_transforms()
    match_added_transform_reflect_horizontal_2d(xform, test_key)

    xform = AutomataTransforms(base_universe_instance_2d())
    test_key = 'rotate90'
    test_matrix = ROTATE_MATRIX_2D_90
    xform.add_transform_cycle(test_key, test_matrix)
    match_added_transform_rotate_2d(xform, test_key)
    expected_message = "Only a single transform cycle '{}' exists. Nothing to base" \
        " combinations on".format(test_key)
    with pytest.raises(ValueError, match=re.compile(expected_message)):
        xform.generate_combination_transforms()
    match_added_transform_rotate_2d(xform, test_key)
# def test_generate_combinations_results_1d() -> None:
    # """verify generated results for 1d universe"""
    # xform = AutomataTransforms(base_universe_instance_1d())
    # # not possible with currently available neighbourhoods and 1d transforms.
    # # only a single transform is possible, to no way to combine multiple transforms.
def test_generate_combinations_reflection_results_2d() -> None:
    """verify combination results for 2d universe with reflections"""
    xform = AutomataTransforms(base_universe_instance_2d())
    # add 2 transforms (cycles), to be able to combine to generate additional patterns
    test_key1 = 'horizontal'
    test_matrix1 = REFLECTION_MATRIX_2D_HORIZONTAL
    test_key2 = 'vertical'
    test_matrix2 = REFLECTION_MATRIX_2D_VERTICAL
    xform.add_transform_cycle(test_key1, test_matrix1)
    match_added_transform_reflect_horizontal_2d(xform, test_key1)
    xform.add_transform_cycle(test_key2, test_matrix2)
    # before generating combinations
    # assert len(xform._transform_cycles) == 3 # base `None` plus 2 cycles added here
    # assert tuple(xform._transform_cycles.keys()) == (None, test_key1, test_key2)
    assert xform._transform_cycles == {None: 0, test_key1: 2, test_key2: 2}
    assert len(xform._index_to_transform) == 2
    assert xform._index_to_transform == {
        (test_key1, 1): REFLECTION_MATRIX_2D_HORIZONTAL,
        (test_key2, 1): REFLECTION_MATRIX_2D_VERTICAL,
    }
    assert xform._transform_to_index == {
        REFLECTION_MATRIX_2D_HORIZONTAL: (test_key1, 1),
        REFLECTION_MATRIX_2D_VERTICAL: (test_key2, 1),
    }
    assert frozenset(xform._index_to_transform.keys()) == frozenset(
        xform._transform_to_index.values())
    assert frozenset(xform._index_to_transform.values()) == frozenset(
        xform._transform_to_index.keys())
    assert len(xform._scratch.primes) == 3 # original plus steps in cycles
    xform.generate_combination_transforms()
    # after generating combinations
    generated_key = frozenset((((test_key1, 1), (test_key2, 1)),
                               ((test_key2, 1), (test_key1, 1))))
    assert xform._transform_cycles == {None: 0, test_key1: 2, test_key2: 2, generated_key : 0}
    # assert len(xform._scratch.primes.values()) == 4 # one extra generated
    assert len(xform._scratch.primes) == 4 # one extra generated
    assert frozenset(xform._index_to_transform.keys()) == frozenset((
        (test_key1, 1), (test_key2, 1), generated_key))
    assert xform._index_to_transform[generated_key] == ((-1, 0), (0, -1)) # horiz+vert reflection
def test_generate_combinations_rotation_horizontal_reflect_results_2d() -> None:
    """verify combination results for 2d universe with rotations and horizontal reflection"""
    xform = AutomataTransforms(base_universe_instance_2d())
    # add 2 transforms (cycles), to be able to combine to generate additional patterns
    r_key = 'rotate 90°'
    r_matrix= ROTATE_MATRIX_2D_90
    h_key = 'horizontal'
    h_matrix = REFLECTION_MATRIX_2D_HORIZONTAL
    start_cycles = {None: 0, r_key: 4, h_key: 2}
    start_i2t = {
        (r_key, 1): r_matrix,
        (r_key, 2): ROTATE_MATRIX_2D_180,
        (r_key, 3): ROTATE_MATRIX_2D_270,
        (h_key, 1): h_matrix,
    }
    generated_key1 = frozenset((
        ((r_key, 1), (h_key, 1)),
        ((h_key, 1), (r_key, 3)),
    ))
    generated_key2 = frozenset((
        ((r_key, 2), (h_key, 1)),
        ((h_key, 1), (r_key, 2))
    ))
    generated_key3 = frozenset((
        ((r_key, 3), (h_key, 1)),
        ((h_key, 1), (r_key, 1))
    ))
    after_cycles = start_cycles.copy()
    after_cycles.update({generated_key1: 0, generated_key2: 0, generated_key3: 0})

    xform.add_transform_cycle(r_key, r_matrix)
    match_added_transform_rotate_2d(xform, r_key)
    xform.add_transform_cycle(h_key, h_matrix)
    # before generating combinations
    # assert len(xform._transform_cycles) == 3 # base `None` plus 2 cycles added here
    # assert tuple(xform._transform_cycles.keys()) == (None, test_key1, test_key2)
    assert xform._transform_cycles == start_cycles
    # assert len(xform._index_to_transform) == 4
    assert xform._index_to_transform == start_i2t
    start_t2i = dict()
    for (key, value) in start_i2t.items():
        start_t2i[value] = key
    assert xform._transform_to_index == start_t2i
    assert frozenset(xform._index_to_transform.keys()) == frozenset(
        xform._transform_to_index.values())
    assert frozenset(xform._index_to_transform.values()) == frozenset(
        xform._transform_to_index.keys())
    assert len(xform._scratch.primes) == 5 # original plus steps in cycles
    xform.generate_combination_transforms()
    # after generating combinations
    assert len(xform._transform_cycles) == 6
    assert xform._transform_cycles == after_cycles
    assert len(xform._scratch.primes) == 8 # three extra generated
    assert len(xform._index_to_transform) == len(start_i2t) + 3
    assert frozenset(xform._index_to_transform.keys()) == frozenset((*start_i2t,
        generated_key1, generated_key2, generated_key3))
    assert xform._index_to_transform[generated_key1] == TRANPOSE_MATRIX_2D
    assert xform._index_to_transform[generated_key2] == REFLECTION_MATRIX_2D_VERTICAL
    assert xform._index_to_transform[generated_key3] == REFLECT_AND_ROTATE_2D
    assert frozenset(xform._scratch.primes.values()) == frozenset((
        (None, 0), *xform._index_to_transform.keys()))
def test_generate_combinations_rotation_2reflect_results_2d() -> None:
    """verify combination results for 2d universe with rotations and 2 reflections"""
    xform = AutomataTransforms(base_universe_instance_2d())
    # add 2 transforms (cycles), to be able to combine to generate additional patterns
    r_key = 'rotate 90°'
    r_matrix= ROTATE_MATRIX_2D_90
    h_key = 'horizontal'
    h_matrix = REFLECTION_MATRIX_2D_HORIZONTAL
    v_key = 'vertical'
    v_matrix = REFLECTION_MATRIX_2D_VERTICAL
    start_cycles = {None: 0, r_key: 4, h_key: 2, v_key: 2}
    start_i2t = {
        (r_key, 1): r_matrix,
        (r_key, 2): ROTATE_MATRIX_2D_180,
        (r_key, 3): ROTATE_MATRIX_2D_270,
        (h_key, 1): h_matrix,
        (v_key, 1): v_matrix,
    }
    generated_key1 = frozenset((
        ((h_key, 1), (r_key, 3)),
        ((r_key, 1), (h_key, 1)),
        ((r_key, 3), (v_key, 1)),
        ((v_key, 1), (r_key, 1)),
    ))
    generated_key2 = frozenset((
        ((h_key, 1), (r_key, 1)),
        ((r_key, 1), (v_key, 1)),
        ((r_key, 3), (h_key, 1)),
        ((v_key, 1), (r_key, 3)),
    ))
    after_cycles = start_cycles.copy()
    after_cycles.update({generated_key1: 0, generated_key2: 0})

    xform.add_transform_cycle(r_key, r_matrix)
    match_added_transform_rotate_2d(xform, r_key)
    xform.add_transform_cycle(h_key, h_matrix)
    xform.add_transform_cycle(v_key, v_matrix)
    # before generating combinations
    # assert len(xform._transform_cycles) == 4 # base `None` plus 3 cycles added here
    # assert tuple(xform._transform_cycles.keys()) == (None, test_key1, test_key2, test_key3)
    assert xform._transform_cycles == start_cycles
    # assert len(xform._index_to_transform) == 5
    assert xform._index_to_transform == start_i2t
    start_t2i = dict()
    for (key, value) in start_i2t.items():
        start_t2i[value] = key
    assert xform._transform_to_index == start_t2i
    assert frozenset(xform._index_to_transform.keys()) == frozenset(
        xform._transform_to_index.values())
    assert frozenset(xform._index_to_transform.values()) == frozenset(
        xform._transform_to_index.keys())
    assert len(xform._scratch.primes) == 6 # original plus steps in cycles
    xform.generate_combination_transforms()
    # after generating combinations
    assert len(xform._transform_cycles) == 6
    assert xform._transform_cycles == after_cycles
    assert len(xform._scratch.primes) == 8 # two extra generated
    assert len(xform._index_to_transform) == len(start_i2t) + 2
    assert frozenset(xform._index_to_transform.keys()) == frozenset(
        (*start_i2t, generated_key1, generated_key2))
    assert xform._index_to_transform[generated_key1] == TRANPOSE_MATRIX_2D
    assert xform._index_to_transform[generated_key2] == REFLECT_AND_ROTATE_2D
    assert frozenset(xform._scratch.primes.values()) == frozenset(
        ((None, 0), *xform._index_to_transform.keys()))
def test_generate_combinations_rot_ref_transpose_results_2d() -> None:
    """verify combination results for 2d universe with rotations and 2 reflections"""
    xform = AutomataTransforms(base_universe_instance_2d())
    # add 2 transforms (cycles), to be able to combine to generate additional patterns
    r_key = 'rotate 90°'
    r_matrix= ROTATE_MATRIX_2D_90
    h_key = 'horizontal'
    h_matrix = REFLECTION_MATRIX_2D_HORIZONTAL
    v_key = 'vertical'
    v_matrix = REFLECTION_MATRIX_2D_VERTICAL
    t_key = 'transpose'
    t_matrix = TRANPOSE_MATRIX_2D
    start_cycles = {None: 0, r_key: 4, h_key: 2, v_key: 2, t_key: 2}
    start_i2t = {
        (r_key, 1): r_matrix,
        (r_key, 2): ROTATE_MATRIX_2D_180,
        (r_key, 3): ROTATE_MATRIX_2D_270,
        (h_key, 1): h_matrix,
        (v_key, 1): v_matrix,
        (t_key, 1): t_matrix,
    }
    generated_key = frozenset((
        ((h_key, 1), (r_key, 1)),
        ((r_key, 1), (v_key, 1)),
        ((r_key, 2), (t_key, 1)),
        ((r_key, 3), (h_key, 1)),
        ((v_key, 1), (r_key, 3)),
        ((t_key, 1), (r_key, 2)),
    ))
    after_cycles = start_cycles.copy()
    after_cycles[generated_key] = 0

    xform.add_transform_cycle(r_key, r_matrix)
    match_added_transform_rotate_2d(xform, r_key)
    xform.add_transform_cycle(h_key, h_matrix)
    xform.add_transform_cycle(v_key, v_matrix)
    xform.add_transform_cycle(t_key, t_matrix)
    # before generating combinations
    # assert len(xform._transform_cycles) == 4 # base `None` plus 3 cycles added here
    # assert tuple(xform._transform_cycles.keys()) == (None, test_key1, test_key2, test_key3)
    assert xform._transform_cycles == start_cycles
    # assert len(xform._index_to_transform) == 5
    assert xform._index_to_transform == start_i2t
    start_t2i = dict()
    for (key, value) in start_i2t.items():
        start_t2i[value] = key
    assert xform._transform_to_index == start_t2i
    assert frozenset(xform._index_to_transform.keys()) == frozenset(
        xform._transform_to_index.values())
    assert frozenset(xform._index_to_transform.values()) == frozenset(
        xform._transform_to_index.keys())
    assert len(xform._scratch.primes) == 7 # original plus steps in cycles
    xform.generate_combination_transforms()
    # after generating combinations
    assert len(xform._transform_cycles) == 6
    assert xform._transform_cycles == after_cycles
    assert len(xform._scratch.primes) == 8 # two extra generated
    assert len(xform._index_to_transform) == len(start_i2t) + 1
    assert frozenset(xform._index_to_transform.keys()) == frozenset(
        (*start_i2t, generated_key))
    assert xform._index_to_transform[generated_key] == REFLECT_AND_ROTATE_2D
    assert frozenset(xform._scratch.primes.values()) == frozenset((
        (None, 0), *xform._index_to_transform.keys()))

def no_test_manual() -> None:
    """verify testing logic and «exception» results"""
    xform = AutomataTransforms(base_universe_instance_1d())
    test_key = 'test02'
    test_matrix = REFLECT_MATRIX_1D
    xform.add_transform_cycle(test_key, test_matrix)
    # xform = AutomataTransforms(base_universe_instance_2d())
    # assert len(BAD_TRANSFORM_1D) > 0
    # for matrix in BAD_TRANSFORM_1D:
    #     # minimal check: details verified in tests for validate_matrix, valid_address
    #     with pytest.raises(ValueError):
    #         xform._check_transform_matrix(matrix)

# cSpell:ignore horiz
