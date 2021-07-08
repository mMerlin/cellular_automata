#!/usr/bin/env python
# coding=utf-8
# pylint: disable=W0212

'''
regression tests for cellular automata universe transform management
'''

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
NOT_MATRICES = (
    *NOT_INTEGER_TYPE_SAMPLES,
    (1,2),
    (('a',)),
)
NOT_MATRICES_1D = (
    IDENTITY_MATRIX_2D,
    IDENTITY_MATRIX_3D,
    ZERO_MATRIX_2D,
    ZERO_MATRIX_3D,
    ((1,2),(2,3)),
    ((1,2,3),(2,3,4),(4,5,6)),
)
BAD_TRANSFORM_1D = (
    ((3,),),
    # IDENTITY_MATRIX_1D, # identity matrix passed initial neighbourhood -> neighbourhood check
    ZERO_MATRIX_1D,
    ((-7,),),
)
BAD_MATRIX_VECTOR_COUNT_1D = (
    ((1,),(2,)),
)

# utility functions to use in tests

def not_matrix_typeerror_1d() -> tuple[object]:
    '''combine test cases for invalid 1d matrices that raise TypeError'''
    assert len(BAD_TYPE_NOT_ITERABLE) > 0
    assert len(BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE) > 0
    assert len(BAD_MATRIX_VECTOR_COUNT_1D) > 0
    assert len(NOT_MATRICES) > 0
    return (
        *BAD_TYPE_NOT_ITERABLE,
        *(case for (case, _err_details) in BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE),
        *(case for case in BAD_MATRIX_VECTOR_COUNT_1D),
        *(case for case in NOT_MATRICES),
    )
def bad_matrix_valueerror_1d() -> tuple[tuple[int]]:
    '''combine test cases for invalid 1d matrices that raise ValueError'''
    assert len(NOT_MATRICES_1D) > 0
    assert len(BAD_TRANSFORM_1D) > 0
    return (
        *NOT_MATRICES_1D,
        *BAD_TRANSFORM_1D,
        IDENTITY_MATRIX_1D,
    )
def match_added_transform_reflect_1d(instance: AutomataTransforms, key: Hashable) -> None:
    '''verify instance properties with (only) reflection transform added'''
    # assert len(xform._transform_cycles) == 2
    # assert test_key in xform._transform_cycles
    # assert xform._transform_cycles[test_key] == 2
    assert instance._transform_cycles == {None: 0, key: 2}
    # assert xform._transform_cycles is None
    # assert len(xform._transform_to_index) == 1
    assert instance._transform_to_index[REFLECT_MATRIX_1D] == _base_index(key)
    # assert xform._transform_to_index is None
    assert instance._primes_transformed == {
        instance._prime_cells: TransformSequence(None, 0),
        frozenset((-cell[0],) for cell in instance._prime_cells): _base_index(key)}
    # assert xform._primes_transformed is None
    # assert len(xform._index_to_transform) == 1
    # assert _base_index(test_key) in xform._index_to_transform
    # assert xform._index_to_transform[_base_index(test_key)] == REFLECT_MATRIX_1D
    assert instance._index_to_transform == {_base_index(key): REFLECT_MATRIX_1D}
    # assert instance._index_to_transform is None

# checks for the standalone (non-class) functions in the module

def test_base_index() -> None:
    '''verify result for base index creation'''
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
    '''verify exception for (only) non hashable key values'''
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
    '''check good and bad values for creating TransformSequence keys'''
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
    '''verify creation and basic properties of a 1d transform management instance'''
    xform = AutomataTransforms(base_universe_instance_1d())
    assert isinstance(xform, AutomataTransforms)
    assert xform._universe.dimensions == 1
    assert isinstance(hash(xform), int)
def test_create_transforms_2d() -> None:
    '''verify creation and basic properties of a 2d transform management instance'''
    xform = AutomataTransforms(base_universe_instance_2d())
    assert isinstance(xform, AutomataTransforms)
    assert xform._universe.dimensions == 2
    assert isinstance(hash(xform), int)
def test_create_transforms_3d() -> None:
    '''verify creation and basic properties of a 3d transform management instance'''
    xform = AutomataTransforms(base_universe_instance_3d())
    assert isinstance(xform, AutomataTransforms)
    assert xform._universe.dimensions == 3
    assert isinstance(hash(xform), int)

def test_check_transform_not_matrix_1d() -> None:
    '''verify raised exception for non-matrix arguments'''
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
def test_check_transform_matrix_not_transform_1d() -> None:
    '''verify raised exception for matrix that is not a transform'''
    xform = AutomataTransforms(base_universe_instance_1d())
    assert len(BAD_TRANSFORM_1D) > 0
    # manual_test = ( ((1,),), ) # identity matrix passes neighbourhood -> neighbourhood check
    for matrix in BAD_TRANSFORM_1D:
        assert isinstance(matrix[0][0], int)
        expected_message = "transformation of neighbourhood with {} is not the " \
            "neighbourhood".format(matrix)
        with pytest.raises(ValueError) as excinfo:
            xform._check_transform_matrix(matrix)
        assert excinfo.value.args == (
            frozenset({(matrix[0][0],), (-matrix[0][0],)}), expected_message)
def test_check_transform_identity_matrix_1d() -> None:
    '''verify raised exception for identity matrix'''
    xform = AutomataTransforms(base_universe_instance_1d())
    expected_message = "transform {} matches the identity matrix".format(IDENTITY_MATRIX_1D)
    # expanded message contains regex pattern, causing match to fail
    # with pytest.raises(ValueError, match=re.compile(expected_message)):
    with pytest.raises(ValueError) as excinfo:
        xform._check_transform_matrix(IDENTITY_MATRIX_1D)
    assert excinfo.value.args[0] == expected_message

def test_make_transform_hashable_bad_transform_1d() -> None:
    '''verify exception raised for bad transformation matrix'''
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
# def test_make_transform_hashable_not_hashable_1d() -> None:
    # '''verify exception raised for bad transformation matrix'''
    # # Unable to create a test for this case. Any matrix that matches this condition also
    # # matches earlier test cases
    # xform = AutomataTransforms(base_universe_instance_1d())
def test_make_transform_hashable_results_1d() -> None:
    '''verify correct hashable transform results'''
    xform = AutomataTransforms(base_universe_instance_1d())
    assert len(GOOD_MATRICES_1D) > 0
    for case in GOOD_MATRICES_1D:
        hashable = xform._make_transform_hashable(case)
        assert isinstance(hashable, tuple)
        assert isinstance(hash(hashable), int)

def test_add_transform_cycle_bad_arguments_1d() -> None:
    '''verify exception raised for bad key or transform arguments'''
    xform = AutomataTransforms(base_universe_instance_1d())
    assert len(BASE_INDEX_BAD_KEYS) > 0
    for test_key in BASE_INDEX_BAD_KEYS:
        # minimal check: details verified in tests for check_transform_key
        with pytest.raises(TypeError):
            xform.add_transform_cycle(test_key, REFLECT_MATRIX_1D)
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
def test_add_transform_cycle_duplicate_key_1d() -> None:
    '''verify raised exceptions when matching data already stored'''
    xform = AutomataTransforms(base_universe_instance_1d())
    test_key = 'test01'
    xform._transform_cycles[test_key] = None # setup to cause duplicate error
    assert test_key in xform._transform_cycles
    expected_message = "Cycle identifier '{}' already in use".format(test_key)
    with pytest.raises(ValueError, match=re.compile(expected_message)):
        xform.add_transform_cycle(test_key, REFLECT_MATRIX_1D)
    xform._transform_cycles.pop(test_key)
    assert test_key not in xform._transform_cycles

    test_index = _base_index(test_key)
    xform._index_to_transform[test_index] = None # setup to cause duplicate error
    assert test_index in xform._index_to_transform
    expected_message = "Base cycle index {} already in use".format(test_index)
    with pytest.raises(ValueError) as excinfo:
        xform.add_transform_cycle(test_key, REFLECT_MATRIX_1D)
    assert excinfo.value.args == (expected_message,)
    xform._index_to_transform.pop(test_index)
    assert test_index not in xform._index_to_transform

    xform._transform_to_index[REFLECT_MATRIX_1D] = test_index
    assert REFLECT_MATRIX_1D in xform._transform_to_index
    expected_message = "'{}' transform matches {}: {}".format(
        test_key, test_index, REFLECT_MATRIX_1D)
    # with pytest.raises(ValueError, match=re.compile(expected_message)):
    with pytest.raises(ValueError) as excinfo:
        xform.add_transform_cycle(test_key, REFLECT_MATRIX_1D)
    assert excinfo.value.args == (expected_message,)
    xform._transform_to_index.pop(REFLECT_MATRIX_1D)
    assert REFLECT_MATRIX_1D not in xform._transform_to_index
def test_add_transfrom_cycle_reflection_result_1d() -> None:
    '''verify stored results after adding reflection transform'''
    xform = AutomataTransforms(base_universe_instance_1d())
    test_key = 'reflection'
    xform.add_transform_cycle(test_key, REFLECT_MATRIX_1D)
    match_added_transform_reflect_1d(xform, test_key)

    expected_message = "Cycle identifier '{}' already in use".format(test_key)
    with pytest.raises(ValueError, match=re.compile(expected_message)):
        xform.add_transform_cycle(test_key, REFLECT_MATRIX_1D)
    match_added_transform_reflect_1d(xform, test_key)

    test_key2 = 'other'
    test_index = _base_index(test_key)
    expected_message = "'{}' transform matches {}: {}".format(
        test_key2, test_index, REFLECT_MATRIX_1D)
    with pytest.raises(ValueError) as excinfo:
        xform.add_transform_cycle(test_key2, REFLECT_MATRIX_1D)
    assert excinfo.value.args == (expected_message,)
    match_added_transform_reflect_1d(xform, test_key)

def test_generate_combinations_not_possible_1d() -> None:
    '''verify exceptions attempting to search combinations for 1d transforms'''
    xform = AutomataTransforms(base_universe_instance_1d())
    with pytest.raises(ValueError, match=re.compile(
            "No transforms have been added yet. Nothing to base generated transforms on")):
        xform.generate_combination_transforms()
    test_key = 'reflection'
    xform.add_transform_cycle(test_key, REFLECT_MATRIX_1D)
    match_added_transform_reflect_1d(xform, test_key)
    expected_message = "Only a single transform cycle '{}' exists. Nothing to base" \
        " combinations on".format(test_key)
    with pytest.raises(ValueError, match=re.compile(expected_message)):
        xform.generate_combination_transforms()
    match_added_transform_reflect_1d(xform, test_key)

def no_test_manual() -> None:
    '''verify testing logic and «exception» results'''
    xform = AutomataTransforms(base_universe_instance_1d())
    test_key = 'test02'
    xform.add_transform_cycle(test_key, REFLECT_MATRIX_1D)
    # xform = AutomataTransforms(base_universe_instance_2d())
    # assert len(BAD_TRANSFORM_1D) > 0
    # for matrix in BAD_TRANSFORM_1D:
    #     # minimal check: details verified in tests for validate_matrix, valid_address
    #     with pytest.raises(ValueError):
    #         xform._check_transform_matrix(matrix)
