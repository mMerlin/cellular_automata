#!/usr/bin/env python
# coding=utf-8
# pylint: disable=E1120,E1121

'''
regression tests for creation of cellular automata universe instances
'''

# import unittest
from typing import Iterable
import re
from _pytest._code.code import ExceptionInfo
import pytest
from automata_universe import AutomataUniverse

# remove need to add parent to path?
# `pipenv run python -m pytest «»`
# ? -q

NEIGHBOURHOOD_1D = ((-1,),(1,))
NEIGHBOURHOOD_2D = (
    (-1,-1), (-1,0), (-1,1),
    ( 0,-1),         ( 0,1),
    ( 1,-1), ( 1,0), ( 1,1),
)
# NEIGHBOURHOOD_3D = ()
GOOD_NEIGHBOURHOOD = (
    NEIGHBOURHOOD_1D,
    NEIGHBOURHOOD_2D,
    list(NEIGHBOURHOOD_1D),
    list(NEIGHBOURHOOD_2D),
    set(NEIGHBOURHOOD_1D),
    set(NEIGHBOURHOOD_2D),
    frozenset(NEIGHBOURHOOD_1D),
    frozenset(NEIGHBOURHOOD_2D),
)
BAD_NEIGHBOURHOOD_TYPE_NOT_ITER = (
    (None, 'NoneType'),
    (1, 'int'),
    (-1, 'int'),
    (2.5, 'float'),
    (True, 'bool'),
)
BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC = (
    tuple(),
    (0,),
    (tuple()),
    ((1,0),),
)
BAD_NEIGHBOURHOOD_TOO_FEW = (
    "",
    {1: 'test'},
    {(1,0): 'test'},
    BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[0],
    BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[1],
    BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[2],
    BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[3],
    list(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[0]),
    list(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[1]),
    list(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[2]),
    list(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[3]),
    set(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[0]),
    set(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[1]),
    set(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[2]),
    set(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[3]),
    frozenset(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[0]),
    frozenset(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[1]),
    frozenset(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[2]),
    frozenset(BAD_NEIGHBOURHOOD_TOO_FEW_ITERABLE_SRC[3]),
)
BAD_NEIGHBOURHOOD_NOT_UNIQUE_ITERABLE_SRC = (
    (0, 0),
    (tuple(), tuple()),
)
BAD_NEIGHBOURHOOD_NOT_UNIQUE = (
    "test",
    BAD_NEIGHBOURHOOD_NOT_UNIQUE_ITERABLE_SRC[0],
    BAD_NEIGHBOURHOOD_NOT_UNIQUE_ITERABLE_SRC[1],
    list(BAD_NEIGHBOURHOOD_NOT_UNIQUE_ITERABLE_SRC[0]),
    list(BAD_NEIGHBOURHOOD_NOT_UNIQUE_ITERABLE_SRC[1]),
)
BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE = (
    ("step", type("")),
    ((1, 0), type(0)),
    (frozenset([0, 1]), type(0)),
    (set([0, 1]), type(0)),
    ({1: 'test', 2: 'test'}, type(0)),
    ((1, (1,2)), type(0)),
)
BAD_NEIGHBOURHOOD_ADDR_NOT_TUPLE = (
    ((99999, (1,2)), type(0)),
    # (("zz", (1,2)), type("")), # sometimes reports as "element not tuple"
    # (("zzz", (1,2), (-8, 0)), type("")), # sometimes reports as element not tuple
    ((frozenset([99999]), (1,2)), type(frozenset())),
)
BAD_NEIGHBOURHOOD_DIMENSIONS = (
    ((tuple(), (1, 2)), 2, 0),
    (((1, 2), tuple()), 2, 0),
)
BAD_NEIGHBOURHOOD_COORDINATE = (
    ((("A", 1), (1,2)), type("")),
    (((0, "A"), (1,2)), type("")),
    (((0, tuple()), (1,2)), type(tuple())),
    (((0, frozenset((9,8))), (1,2)), type(frozenset())),
    (((1.5, 1), (1,2)), type(0.1)),
)
BAD_NEIGHBOURHOOD_WITH_ORIGIN = (
    (*NEIGHBOURHOOD_1D, tuple([0]*len(NEIGHBOURHOOD_1D[0]))),
    (tuple([0]*len(NEIGHBOURHOOD_2D[0])), *NEIGHBOURHOOD_2D),
)
BAD_NEIGHBOURHOOD_WITHOUT_SYMMETRIC = (
    ((*NEIGHBOURHOOD_1D, (-2,)), (-2,)),
    ((*NEIGHBOURHOOD_2D, (-2,0)), (-2,0)),
)
BAD_NEIGHBOURHOOD_NO_COORDINATE = (
    (tuple(), (-1,)),
    (tuple(), (-1,-1)),
)

GOOD_COUNTS_SRC = (
    tuple(),
    (1,),
    (2, 3),
    (1, 2, 3, 4, 5, 6, 7, 8),
)
GOOD_BIRTH_COUNTS = (
    GOOD_COUNTS_SRC[0],
    GOOD_COUNTS_SRC[1],
    GOOD_COUNTS_SRC[2],
    GOOD_COUNTS_SRC[3],
    list(GOOD_COUNTS_SRC[0]),
    list(GOOD_COUNTS_SRC[1]),
    list(GOOD_COUNTS_SRC[2]),
    list(GOOD_COUNTS_SRC[3]),
    set(GOOD_COUNTS_SRC[0]),
    set(GOOD_COUNTS_SRC[1]),
    set(GOOD_COUNTS_SRC[2]),
    set(GOOD_COUNTS_SRC[3]),
    frozenset(GOOD_COUNTS_SRC[0]),
    frozenset(GOOD_COUNTS_SRC[1]),
    frozenset(GOOD_COUNTS_SRC[2]),
    frozenset(GOOD_COUNTS_SRC[3]),
)
GOOD_SURVIVAL_COUNTS_SRC = (
    (0, 1, 2, 3, 4, 5, 6, 7, 8),
)
GOOD_SURVIVAL_COUNTS = (
    *GOOD_BIRTH_COUNTS,
    GOOD_SURVIVAL_COUNTS_SRC[0],
    list(GOOD_SURVIVAL_COUNTS_SRC[0]),
    set(GOOD_SURVIVAL_COUNTS_SRC[0]),
    frozenset(GOOD_SURVIVAL_COUNTS_SRC[0]),
)
BAD_COUNT_RULES_SRC = (
    (-1,),
    (9,),
    (0.5,),
    (1, 1),
)
BAD_COUNT_RULES_TYPE_NOT_ITER = (
    (None, 'NoneType'),
    (1, 'int'),
    (-1, 'int'),
    (2.5, 'float'),
    (True, 'bool'),
)
BAD_COUNT_RULES_NOT_UNIQUE = (
    "test",
    (*GOOD_COUNTS_SRC[1], GOOD_COUNTS_SRC[1][0]),
    (*GOOD_COUNTS_SRC[2], GOOD_COUNTS_SRC[2][1]),
    (GOOD_COUNTS_SRC[3][3], *GOOD_COUNTS_SRC[3]),
    BAD_COUNT_RULES_SRC[3],
    list(BAD_COUNT_RULES_SRC[3]),
)
BAD_COUNT_OUTSIDE_RANGE = (
    ((*GOOD_COUNTS_SRC[1], -1), -1),
    ((9, *GOOD_COUNTS_SRC[2]), 9),
    (BAD_COUNT_RULES_SRC[0], -1),
    (BAD_COUNT_RULES_SRC[1], 9),
    (list(BAD_COUNT_RULES_SRC[0]), -1),
    (list(BAD_COUNT_RULES_SRC[1]), 9),
    (set(BAD_COUNT_RULES_SRC[0]), -1),
    (set(BAD_COUNT_RULES_SRC[1]), 9),
    (frozenset(BAD_COUNT_RULES_SRC[0]), -1),
    (frozenset(BAD_COUNT_RULES_SRC[1]), 9),
)
BAD_COUNT_NOT_INTEGER = (
    ("step", type("")),
    ((None,), type(None)),
    ((2, None), type(None)),
    ((2, "test"), type("")),
    ((2, 1.5), type(0.1)),
    (BAD_COUNT_RULES_SRC[2], type(0.1)),
    (list(BAD_COUNT_RULES_SRC[2]), type(0.1)),
    (set(BAD_COUNT_RULES_SRC[2]), type(0.1)),
    (frozenset(BAD_COUNT_RULES_SRC[2]), type(0.1)),
)
BAD_BIRTH_COUNTS_SRC = (
    (0,),
    (0, 3),
    (*GOOD_COUNTS_SRC[0], 0),
    (*GOOD_COUNTS_SRC[1], 0),
    (*GOOD_COUNTS_SRC[2], 0),
    (*GOOD_COUNTS_SRC[3], 0),
)
BAD_BIRTH_COUNTS = (
    BAD_BIRTH_COUNTS_SRC[0],
    BAD_BIRTH_COUNTS_SRC[1],
    BAD_BIRTH_COUNTS_SRC[2],
    BAD_BIRTH_COUNTS_SRC[3],
    BAD_BIRTH_COUNTS_SRC[4],
    BAD_BIRTH_COUNTS_SRC[5],
    list(BAD_BIRTH_COUNTS_SRC[0]),
    list(BAD_BIRTH_COUNTS_SRC[1]),
    list(BAD_BIRTH_COUNTS_SRC[2]),
    list(BAD_BIRTH_COUNTS_SRC[3]),
    list(BAD_BIRTH_COUNTS_SRC[4]),
    list(BAD_BIRTH_COUNTS_SRC[5]),
    set(BAD_BIRTH_COUNTS_SRC[0]),
    set(BAD_BIRTH_COUNTS_SRC[1]),
    set(BAD_BIRTH_COUNTS_SRC[2]),
    set(BAD_BIRTH_COUNTS_SRC[3]),
    set(BAD_BIRTH_COUNTS_SRC[4]),
    set(BAD_BIRTH_COUNTS_SRC[5]),
    frozenset(BAD_BIRTH_COUNTS_SRC[0]),
    frozenset(BAD_BIRTH_COUNTS_SRC[1]),
    frozenset(BAD_BIRTH_COUNTS_SRC[2]),
    frozenset(BAD_BIRTH_COUNTS_SRC[3]),
    frozenset(BAD_BIRTH_COUNTS_SRC[4]),
    frozenset(BAD_BIRTH_COUNTS_SRC[5]),
)

# Data for testing automata cell «address» group operations
GOOD_POPULATION_SRC = (
)
GOOD_POPULATION = (
)
BAD_POPULATION_SRC = (
)
BAD_POPULATION = (
)

BAD_EXTRA_ARGUMENTS = (
    (None,),
    (0,),
    (1.2,),
    (tuple(),),
    (None, None),
)

def verify_general_exception_tuple(ex: ExceptionInfo, length: int) -> None:
    '''check the structure of a custom exception'''
    assert isinstance(ex, ExceptionInfo)
    assert isinstance(ex.type, type)
    # assert isinstance(excinfo.type, ValueError)
    # assert isinstance(ex.value, ValueError)
    # assert excinfo.type == "<class 'ValueError'>"
    # assert str(ex.type) == "<class 'ValueError'>"
    assert isinstance(ex.value.args, tuple)
    assert len(ex.value.args) == 1
    assert isinstance(ex.value.args[0], tuple)
    assert len(ex.value.args[0]) == length
def plural_suffix(count: int) -> str:
    '''"s" when count is not one'''
    suffix = ''
    if count != 1:
        suffix = 's'
    return suffix
def required_positional_message(missing: int) -> str:
    return re.compile("__init__[(][)] missing {} required positional argument{}: .*".format(
        missing, plural_suffix(missing)))

def test_no_arguments() -> None:
    '''no arguments supplied when 3 expected'''
    with pytest.raises(TypeError, match=required_positional_message(3)):
        AutomataUniverse()
def test_1_argument() -> None:
    '''only one argument supplied when 3 expected'''
    bad_neighbourhoods = [case for (case, _err) in BAD_NEIGHBOURHOOD_TYPE_NOT_ITER]
    bad_neighbourhoods.extend(BAD_NEIGHBOURHOOD_TOO_FEW)
    for (case, _err) in BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE:
        bad_neighbourhoods.append(case)
    for (case, _dim1, _dim2) in BAD_NEIGHBOURHOOD_DIMENSIONS:
        bad_neighbourhoods.append(case)
    bad_neighbourhoods.extend(GOOD_NEIGHBOURHOOD)
    assert len(bad_neighbourhoods) > 0
    for neighbourhood in bad_neighbourhoods:
        with pytest.raises(TypeError, match=required_positional_message(2)):
            AutomataUniverse(neighbourhood)
def test_2_arguments() -> None:
    '''two arguments supplied when 3 expected'''
    bad_survivals = [case for (case, _err) in BAD_COUNT_RULES_TYPE_NOT_ITER]
    bad_survivals.extend(BAD_COUNT_RULES_NOT_UNIQUE)
    for (case, _coord) in BAD_COUNT_OUTSIDE_RANGE:
        bad_survivals.append(case)
    for (case, _err) in BAD_COUNT_NOT_INTEGER:
        bad_survivals.append(case)
    bad_survivals.extend(BAD_BIRTH_COUNTS)
    bad_survivals.extend(GOOD_SURVIVAL_COUNTS)
    bad_survivals.extend(GOOD_BIRTH_COUNTS)
    bad_survivals.extend(GOOD_SURVIVAL_COUNTS)
    assert len(bad_survivals) > 0
    for survival in bad_survivals:
        with pytest.raises(TypeError, match=required_positional_message(1)):
            AutomataUniverse(NEIGHBOURHOOD_2D, survival)
def test_4_arguments() -> None:
    '''four arguments supplied when 3 expected'''
    for extra in BAD_EXTRA_ARGUMENTS:
        with pytest.raises(TypeError, match=re.compile(
            "__init__[(][)] takes 4 positional arguments but {} were given".format(
                len(extra) + 4))):
            AutomataUniverse(NEIGHBOURHOOD_2D, [2,3], [3], *extra)
def test_non_iterable_neighbourhood() -> None:
    '''neighbourhood parameter is not iterable'''
    assert len(BAD_NEIGHBOURHOOD_TYPE_NOT_ITER) > 0
    for (neighbourhood, err_obj) in BAD_NEIGHBOURHOOD_TYPE_NOT_ITER:
        expected = "'{}' object is not iterable".format(err_obj)
        with pytest.raises(TypeError, match=re.compile(expected)):
            AutomataUniverse(neighbourhood, GOOD_SURVIVAL_COUNTS[0], GOOD_BIRTH_COUNTS[0])
def test_too_few_neighbours() -> None:
    '''neighbourhood parameter is iterable but too short'''
    assert len(BAD_NEIGHBOURHOOD_TOO_FEW) > 0
    for neighbourhood in BAD_NEIGHBOURHOOD_TOO_FEW:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_SURVIVAL_COUNTS[0], GOOD_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (len(neighbourhood),
            "neighbourhood does not contain enough neighbours: at least 2 required")
def test_neighbours_not_unique() -> None:
    '''iterable neighbourhood contains duplicate elements'''
    assert len(BAD_NEIGHBOURHOOD_NOT_UNIQUE) > 0
    for neighbourhood in BAD_NEIGHBOURHOOD_NOT_UNIQUE:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_SURVIVAL_COUNTS[0], GOOD_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0][0] == len(neighbourhood)
        assert excinfo.value.args[0][1] < len(neighbourhood)
        assert excinfo.value.args[0][2] == "Neighbourhood addresses are not unique"
def test_neighbour_element_not_tuple() -> None:
    '''iterable neighbourhood contains non-tuple element'''
    assert len(BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE) > 0
    for (neighbourhood, err_type) in BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(TypeError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_SURVIVAL_COUNTS[0], GOOD_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type, "element of neighbourhood is not a tuple")
def do_addr_test(nei, err_type) -> None:
    '''see if can get report showing neighbourhood set'''
    assert isinstance(nei, Iterable)
    with pytest.raises(TypeError) as excinfo:
        AutomataUniverse(nei, GOOD_SURVIVAL_COUNTS[0], GOOD_BIRTH_COUNTS[0])
    verify_general_exception_tuple(excinfo, 2)
    assert excinfo.value.args[0] == (err_type, "automata universe address is not a tuple")
def test_neighbour_address_not_tuple() -> None:
    '''iterable neighbourhood contains non-tuple address'''
    assert len(BAD_NEIGHBOURHOOD_ADDR_NOT_TUPLE) > 0
    for (neighbourhood, err_type) in BAD_NEIGHBOURHOOD_ADDR_NOT_TUPLE:
        do_addr_test(neighbourhood, err_type)
def test_bad_neighbour_dimensions() -> None:
    '''iterable neighbourhood address element different dimensions'''
    assert len(BAD_NEIGHBOURHOOD_DIMENSIONS) > 0
    for (neighbourhood, dimensions_uni, dimensions_addr) in BAD_NEIGHBOURHOOD_DIMENSIONS:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_SURVIVAL_COUNTS[0], GOOD_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0] == (dimensions_uni, dimensions_addr, "automata universe"
            " address does not have the same number of dimensions as the universe")
def test_coordinate_not_integer() -> None:
    '''iterable neighbourhood contains non-integer coordinate'''
    assert len(BAD_NEIGHBOURHOOD_COORDINATE) > 0
    for (neighbourhood, err_type) in BAD_NEIGHBOURHOOD_COORDINATE:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(TypeError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_SURVIVAL_COUNTS[0], GOOD_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type,
            "automata universe address coordinate is not an integer")
def test_origin_in_neighbourhood() -> None:
    '''iterable neighbourhood contains universe origin'''
    assert len(BAD_NEIGHBOURHOOD_WITH_ORIGIN) > 0
    for neighbourhood in BAD_NEIGHBOURHOOD_WITH_ORIGIN:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_SURVIVAL_COUNTS[0], GOOD_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (tuple([0]*len(neighbourhood[0])),
            "the universe origin is not a valid neighbourhood address")
def test_no_symmetric_address() -> None:
    '''iterable neighbourhood missing symmetric address'''
    assert len(BAD_NEIGHBOURHOOD_WITHOUT_SYMMETRIC) > 0
    for (neighbourhood, address) in BAD_NEIGHBOURHOOD_WITHOUT_SYMMETRIC:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_SURVIVAL_COUNTS[0], GOOD_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (address, "no symmetric address in neighbourhood")
def test_no_coordinate() -> None:
    '''iterable neighbourhood first address without coordinate'''
    assert len(BAD_NEIGHBOURHOOD_NO_COORDINATE) > 0
    for neighbourhood in BAD_NEIGHBOURHOOD_NO_COORDINATE:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(TypeError, match=re.compile(
                "neighbourhood address does not contain at least 1 coordinate")):
            AutomataUniverse(neighbourhood, GOOD_SURVIVAL_COUNTS[0], GOOD_BIRTH_COUNTS[0])

def test_not_iterable_survial() -> None:
    '''survival rule parameter is not iterable'''
    assert len(BAD_COUNT_RULES_TYPE_NOT_ITER) > 0
    for (rule, err_obj) in BAD_COUNT_RULES_TYPE_NOT_ITER:
        expected = "'{}' object is not iterable".format(err_obj)
        with pytest.raises(TypeError, match=re.compile(expected)):
            AutomataUniverse(NEIGHBOURHOOD_2D, rule, GOOD_BIRTH_COUNTS[0])
def test_not_iterable_birth() -> None:
    '''birth rule parameter is not iterable'''
    assert len(BAD_COUNT_RULES_TYPE_NOT_ITER) > 0
    for (rule, err_obj) in BAD_COUNT_RULES_TYPE_NOT_ITER:
        expected = "'{}' object is not iterable".format(err_obj)
        with pytest.raises(TypeError, match=re.compile(expected)):
            AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_SURVIVAL_COUNTS[0], rule)
def test_not_unique_survial() -> None:
    '''survival rule parameter contains duplicate entry'''
    assert len(BAD_COUNT_RULES_NOT_UNIQUE) > 0
    for rule in BAD_COUNT_RULES_NOT_UNIQUE:
        assert isinstance(rule, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, rule, GOOD_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0][0] == len(rule)
        assert excinfo.value.args[0][1] < len(rule)
        assert excinfo.value.args[0][2] == "cell propagation rule counts are not unique"
def test_not_unique_birth() -> None:
    '''birth rule parameter contains duplicate entry'''
    assert len(BAD_COUNT_RULES_NOT_UNIQUE) > 0
    for rule in BAD_COUNT_RULES_NOT_UNIQUE:
        assert isinstance(rule, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_SURVIVAL_COUNTS[0], rule)
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0][0] == len(rule)
        assert excinfo.value.args[0][1] < len(rule)
        assert excinfo.value.args[0][2] == "cell propagation rule counts are not unique"
def test_survival_outside_neighbourhood() -> None:
    '''survival rule parameter entry outside neighbourhood size'''
    assert len(BAD_COUNT_OUTSIDE_RANGE) > 0
    for (rule, outside) in BAD_COUNT_OUTSIDE_RANGE:
        assert isinstance(rule, Iterable)
        assert isinstance(outside, int)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, rule, GOOD_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0][0] == outside
        assert excinfo.value.args[0][1] == len(NEIGHBOURHOOD_2D)
        assert excinfo.value.args[0][2] == \
            "cell propagation count case is not between 0 and neighbourhood size"
def test_birth_outside_neighbourhood() -> None:
    '''birth rule parameter entry outside neighbourhood size'''
    assert len(BAD_COUNT_OUTSIDE_RANGE) > 0
    for (rule, outside) in BAD_COUNT_OUTSIDE_RANGE:
        assert isinstance(rule, Iterable)
        assert isinstance(outside, int)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_SURVIVAL_COUNTS[0], rule)
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0][0] == outside
        assert excinfo.value.args[0][1] == len(NEIGHBOURHOOD_2D)
        assert excinfo.value.args[0][2] == \
            "cell propagation count case is not between 0 and neighbourhood size"
def test_survival_not_integer() -> None:
    '''survival rule parameter entry not an integer'''
    assert len(BAD_COUNT_NOT_INTEGER) > 0
    for (rule, err_type) in BAD_COUNT_NOT_INTEGER:
        assert isinstance(rule, Iterable)
        assert isinstance(err_type, type)
        with pytest.raises(TypeError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, rule, GOOD_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type,
            "cell propagation count case is not an integer")
def test_birth_not_integer() -> None:
    '''birth rule parameter entry not an integer'''
    assert len(BAD_COUNT_NOT_INTEGER) > 0
    for (rule, err_type) in BAD_COUNT_NOT_INTEGER:
        assert isinstance(rule, Iterable)
        assert isinstance(err_type, type)
        with pytest.raises(TypeError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_SURVIVAL_COUNTS[0], rule)
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type,
            "cell propagation count case is not an integer")
def test_zero_birth_rule() -> None:
    '''birth rule contains entry for zero neighbours'''
    assert len(BAD_BIRTH_COUNTS) > 0
    for rule in BAD_BIRTH_COUNTS:
        assert isinstance(rule, Iterable)
        with pytest.raises(ValueError, match=re.compile(
                "zero is not a valid birth propagation rule value")):
            AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_SURVIVAL_COUNTS[0], rule)

def test_properties() -> None:
    '''verify the properties of a successful instance creation'''
    test_neighbourhood = NEIGHBOURHOOD_2D
    test_survival = GOOD_SURVIVAL_COUNTS[0]
    test_birth = GOOD_BIRTH_COUNTS[0]
    uni = AutomataUniverse(test_neighbourhood, test_survival, test_birth)
    assert isinstance(uni, AutomataUniverse)
    assert isinstance(uni.dimensions, int)
    assert uni.dimensions == 2 # previously caused pylint crash/trace
    # assert isinstance(uni.neighbourhood, frozenset)
    # assert uni.neighbourhood == frozenset(test_neighbourhood)
    # assert isinstance(uni.neighbourhood_population, int)
    # assert uni.neighbourhood_population == len(test_neighbourhood)
    # assert isinstance(uni.survival_rules, frozenset)
    # assert uni.survival_rules == frozenset(test_survival)
    # assert isinstance(uni.birth_rules, frozenset)
    # assert uni.birth_rules == frozenset(test_birth)
    # assert isinstance(hash(uni), int)
    # assert hasattr(uni, "__hash__")


# def test_exploration() -> None:
#     '''view exact failure for test case'''
#     # AutomataUniverse(NEIGHBOURHOOD_2D, ((-1,)), GOOD_BIRTH_COUNTS[0])
#     AutomataUniverse(NEIGHBOURHOOD_2D, BAD_[0], GOOD_BIRTH_COUNTS[0])
#     # AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_SURVIVAL_COUNTS[0], BAD_COUNT_RULES_NOT_UNIQUE[0])

    # isinstance(AutomataUniverse(), AutomataUniverse)

# if __name__ == '__main__':
#     unittest.main()
