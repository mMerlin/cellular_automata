#!/usr/bin/env python
# coding=utf-8
# pylint: disable=E1120,E1121,W0212,C0302

"""
regression tests for creation of cellular automata universe instances
"""

from typing import Iterable
import re
from _pytest._code.code import ExceptionInfo
import pytest
from common_test_data import (
    NEIGHBOURHOOD_1D,
    NEIGHBOURHOOD_2D,
    NEIGHBOURHOOD_3D,
    IDENTITY_MATRIX_1D,
    IDENTITY_MATRIX_2D,
    IDENTITY_MATRIX_3D,
    NOT_INTEGER_TYPE_SAMPLES,
    ORIGIN_2D,
    TRANPOSE_MATRIX_2D,
    UNIT_VECTOR_2D,

    ZERO_MATRIX_2D,
    ROTATE_MATRIX_2D_90,
    ROTATE_MATRIX_2D_180,
    ROTATE_MATRIX_2D_270,
    REFLECTION_MATRIX_2D_HORIZONTAL,
    REFLECTION_MATRIX_2D_VERTICAL,
    SCALE_3X_5Y_MATRIX,

    GOOD_NEIGHBOURHOOD,
    BAD_TYPE_NOT_ITERABLE,
    BAD_NEIGHBOURHOOD_TOO_FEW,
    BAD_NEIGHBOURHOOD_NOT_UNIQUE,
    BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE,
    BAD_NEIGHBOURHOOD_ADDR_NOT_TUPLE,
    BAD_NEIGHBOURHOOD_DIMENSIONS,
    BAD_NEIGHBOURHOOD_COORDINATE,
    BAD_NEIGHBOURHOOD_WITH_ORIGIN,
    BAD_NEIGHBOURHOOD_WITHOUT_SYMMETRIC,
    BAD_NEIGHBOURHOOD_NO_COORDINATE,

    SQUARE_MATRIX_TRANSFORM_OPERATIONS_2D,
)
from automata_universe import AutomataUniverse

# remove need to add parent to path?
# `pipenv run python -m pytest «»`
# ? -q -v

GOOD_COUNTS_SRC = (
    tuple(),
    (1,),
    (2, 3),
    (1, 2, 3, 4, 5, 6, 7, 8),
)
GOOD_1D_BIRTH_COUNTS = (
    GOOD_COUNTS_SRC[0],
    GOOD_COUNTS_SRC[1],
)
GOOD_1D_SURVIVAL_COUNTS = (
    GOOD_COUNTS_SRC[0],
    GOOD_COUNTS_SRC[1],
)
GOOD_2D_BIRTH_COUNTS = (
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
GOOD_2D_SURVIVAL_COUNTS = (
    *GOOD_2D_BIRTH_COUNTS,
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
    *(((sample,), type(sample)) for sample in NOT_INTEGER_TYPE_SAMPLES
        if not isinstance(sample, (set, dict, list))),
        # if not isinstance(sample, Iterable)),
    ("step", type("s")),
    ((None,), type(None)),
    ((2, None), type(None)),
    ((2, "test"), type("t")),
    ((2, 1.5), type(1.5)),
    (BAD_COUNT_RULES_SRC[2], type(BAD_COUNT_RULES_SRC[2][0])),
    (list(BAD_COUNT_RULES_SRC[2]), type(BAD_COUNT_RULES_SRC[2][0])),
    (set(BAD_COUNT_RULES_SRC[2]), type(BAD_COUNT_RULES_SRC[2][0])),
    (frozenset(BAD_COUNT_RULES_SRC[2]), type(BAD_COUNT_RULES_SRC[2][0])),
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

GOOD_2D_ADDRESS = (
    (1, 2),
    (-1555, 987777),
)
BAD_ADDRESS_NOT_TUPLE = (
    ("step", type("")),
    (-1, type(0)),
    (frozenset((0,)), type(frozenset())),
    (set([0,]), type(set())),
    ({1: 'test', 2: 'test'}, type(dict())),
    (99999, type(0)),
    (frozenset([99999]), type(frozenset())),
)
BAD_ADDRESS_NOT_2_DIMENSIONS = (
    (tuple(), 0),
    ((8,), 1),
    ((2.1,), 1),
    ((1, -2, 4), 3),
)
BAD_ADDRESS_COORD_NOT_INTEGER = (
    ((1, None), type(None)),
    ((None, -1), type(None)),
    ((2, "z"), type("")),
    (("a", -2), type("")),
    ((1, 1.1), type(0.1)),
    ((1.2, -1), type(0.1)),
    ((1, set()), type(set())),
    ((set(), -1), type(set())),
)

GOOD_2D_MATRIX = (
    ((1, 0), (0, 1)),
    ((0, -1), (1, 0)),
    ((99, -9), (11, 0)),
)
BAD_2D_MATRIX_VECTOR_COUNT = (
    ((1,2),),
    ((1,2),(3,4),(5,6)),
)

Q1_7_11 = (7,11)
Q2_7_11 = (-7,11)
Q3_7_11 = (-7,-11)
Q4_7_11 = (7,-11)
Q1_11_7 = (11,7)
Q2_11_7 = (-11,7)
Q3_11_7 = (-11,-7)
Q4_11_7 = (11,-7)
DOT_PRODUCT_2D_VECTOR_TESTS = (
    # («vector», «matrix», «expected_result»),
    (Q1_7_11, ZERO_MATRIX_2D, ORIGIN_2D),
    (ORIGIN_2D, IDENTITY_MATRIX_2D, ORIGIN_2D),
    (Q1_7_11, IDENTITY_MATRIX_2D, Q1_7_11),
    (Q2_7_11, IDENTITY_MATRIX_2D, Q2_7_11),
    (Q3_7_11, IDENTITY_MATRIX_2D, Q3_7_11),
    (Q4_7_11, IDENTITY_MATRIX_2D, Q4_7_11),
    (Q1_7_11, ROTATE_MATRIX_2D_90, Q2_11_7),
    (Q2_7_11, ROTATE_MATRIX_2D_90, Q3_11_7),
    (Q3_7_11, ROTATE_MATRIX_2D_90, Q4_11_7),
    (Q4_7_11, ROTATE_MATRIX_2D_90, Q1_11_7),
    (Q1_11_7, ROTATE_MATRIX_2D_90, Q2_7_11),
    (Q2_11_7, ROTATE_MATRIX_2D_90, Q3_7_11),
    (Q3_11_7, ROTATE_MATRIX_2D_90, Q4_7_11),
    (Q4_11_7, ROTATE_MATRIX_2D_90, Q1_7_11),
    (Q1_7_11, REFLECTION_MATRIX_2D_HORIZONTAL, Q2_7_11),
    (Q2_7_11, REFLECTION_MATRIX_2D_HORIZONTAL, Q1_7_11),
    (Q3_7_11, REFLECTION_MATRIX_2D_HORIZONTAL, Q4_7_11),
    (Q4_7_11, REFLECTION_MATRIX_2D_HORIZONTAL, Q3_7_11),
    (Q1_7_11, REFLECTION_MATRIX_2D_VERTICAL, Q4_7_11),
    (Q2_7_11, REFLECTION_MATRIX_2D_VERTICAL, Q3_7_11),
    (Q3_7_11, REFLECTION_MATRIX_2D_VERTICAL, Q2_7_11),
    (Q4_7_11, REFLECTION_MATRIX_2D_VERTICAL, Q1_7_11),
    (UNIT_VECTOR_2D, SCALE_3X_5Y_MATRIX, (3,5)),
    (UNIT_VECTOR_2D, SCALE_3X_5Y_MATRIX, (3,5)),
    (ORIGIN_2D, SCALE_3X_5Y_MATRIX, ORIGIN_2D),
    (Q1_11_7, TRANPOSE_MATRIX_2D, (7,11)),
    (Q1_11_7, TRANPOSE_MATRIX_2D, (Q1_11_7[1], Q1_11_7[0])),
    (Q2_11_7, TRANPOSE_MATRIX_2D, (Q2_11_7[1], Q2_11_7[0])),
    (Q3_11_7, TRANPOSE_MATRIX_2D, (Q3_11_7[1], Q3_11_7[0])),
    (Q4_11_7, TRANPOSE_MATRIX_2D, (Q4_11_7[1], Q4_11_7[0])),
)

DOT_PRODUCT_2D_SQUARE_EXPECTED = (
    # («matrix», «transform», «expected_result»),
    *SQUARE_MATRIX_TRANSFORM_OPERATIONS_2D,
)

# Data for testing automata cell «address» group operations
GOOD_2D_POPULATION_SRC = (
    set(((1,2), (2,3))),
)
GOOD_2D_POPULATION = (
    GOOD_2D_POPULATION_SRC[0],
    frozenset(GOOD_2D_POPULATION_SRC[0]),
)
BAD_POPULATION_NOT_SET = (
    (None, type(None)),
    (1, type(0)),
    (-1, type(0)),
    (2.5, type(0.1)),
    (True, type(False)),
    (((1,2)), type((1,))),
    ([(1,2)], type([])),
    ({1: (1,2)}, type({1:0})),
)

GOOD_CELL_ADDRESS_NEIGHBOURS = (
    (ORIGIN_2D, frozenset(NEIGHBOURHOOD_2D)),
    (UNIT_VECTOR_2D, frozenset((
        (0,0), (0,1), (0,2),
        (1,0),        (1,2),
        (2,0), (2,1), (2,2),
    ))),
    ((49,99), frozenset((
        (48,98), (48,99), (48,100),
        (49,98),          (49,100),
        (50,98), (50,99), (50,100),
    ))),
)
GOOD_CELL_GROUP_TRANSLATE = (
    *((set(population), ORIGIN_2D, set(population)) for population in GOOD_2D_POPULATION_SRC),
    *((frozenset(population), ORIGIN_2D, set(population)) for population in GOOD_2D_POPULATION_SRC),
    *((set(population), ORIGIN_2D, frozenset(population)) for population in GOOD_2D_POPULATION_SRC),
    *((frozenset(population), ORIGIN_2D, frozenset(population))
        for population in GOOD_2D_POPULATION_SRC),
    (GOOD_2D_POPULATION_SRC[0], UNIT_VECTOR_2D, frozenset(((2, 3),(3, 4)))),
    (GOOD_2D_POPULATION_SRC[0], (-1,-5), frozenset(((0,-3),(1,-2)))),
)
GOOD_CELL_GROUP_TRANSFORM = (
    *((set(population), IDENTITY_MATRIX_2D, set(population))
        for population in GOOD_2D_POPULATION_SRC),
    *((frozenset(population), IDENTITY_MATRIX_2D, set(population))
        for population in GOOD_2D_POPULATION_SRC),
    *((set(population), IDENTITY_MATRIX_2D, frozenset(population))
        for population in GOOD_2D_POPULATION_SRC),
    *((frozenset(population), IDENTITY_MATRIX_2D, frozenset(population))
        for population in GOOD_2D_POPULATION_SRC),
    *((frozenset((vector,)), matrix, frozenset((expected,)))
        for (vector, matrix, expected) in DOT_PRODUCT_2D_VECTOR_TESTS)
)

GOOD_ROTATION_MATRIX = (
    IDENTITY_MATRIX_2D,
    ROTATE_MATRIX_2D_90,
    ROTATE_MATRIX_2D_180,
    ROTATE_MATRIX_2D_270,
)
NOT_ROTATION_MATRIX = (
    ZERO_MATRIX_2D,
    REFLECTION_MATRIX_2D_HORIZONTAL,
    REFLECTION_MATRIX_2D_VERTICAL,
    SCALE_3X_5Y_MATRIX,
)

GOOD_STEP_POPULATIONS = (
    (set(), set()),
    (set([(0,0)]), set()),
    (set([(0,0),(0,1)]), set()),
    (set([(0,0),(0,1),(0,2)]), set([(-1,1),(0,1),(1,1)])),
)

BAD_EXTRA_ARGUMENTS = (
    (None,),
    (0,),
    (1.2,),
    (tuple(),),
    (None, None),
)

def verify_general_exception_tuple(ex: ExceptionInfo, length: int) -> None:
    """check the structure of a custom exception"""
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
    """"s" when count is not one"""
    suffix = ''
    if count != 1:
        suffix = 's'
    return suffix
def required_positional_message(missing: int) -> str:
    return re.compile("__init__[(][)] missing {} required positional argument{}: .*".format(
        missing, plural_suffix(missing)))
def base_universe_instance_1d() -> AutomataUniverse:
    """create and return 1 dimensional universe instance

    Use when further steps do not need the parameters used to create the universe instance

    :returns: standard 1 dimensional automata universe
    :rtype: AutomataUniverse
    """
    test_neighbourhood = NEIGHBOURHOOD_1D
    test_survival = GOOD_1D_SURVIVAL_COUNTS[0]
    test_birth = GOOD_1D_BIRTH_COUNTS[0]
    universe = AutomataUniverse(test_neighbourhood, test_survival, test_birth)
    assert universe.dimensions == 1
    return universe
def base_universe_instance_2d() -> AutomataUniverse:
    """create and return 2 dimensional universe instance

    Use when further steps do not need the parameters used to create the universe instance

    :returns: standard 2 dimensional rectangular grid automata universe
    :rtype: AutomataUniverse
    """
    test_neighbourhood = NEIGHBOURHOOD_2D
    test_survival = frozenset((2,3))
    test_birth = frozenset((3,))
    universe = AutomataUniverse(test_neighbourhood, test_survival, test_birth)
    assert universe.dimensions == 2
    return universe
def base_universe_instance_3d() -> AutomataUniverse:
    """create and return 3 dimensional universe instance

    Use when further steps do not need the parameters used to create the universe instance

    :returns: standard 3 dimensional cube grid automata universe
    :rtype: AutomataUniverse
    """
    test_neighbourhood = NEIGHBOURHOOD_3D
    test_survival = GOOD_2D_SURVIVAL_COUNTS[2] # good enough for base universe
    test_birth = GOOD_2D_BIRTH_COUNTS[2]
    universe = AutomataUniverse(test_neighbourhood, test_survival, test_birth)
    assert universe.dimensions == 3
    return universe

def test_no_arguments() -> None:
    """no arguments supplied when 3 expected"""
    with pytest.raises(TypeError, match=required_positional_message(3)):
        AutomataUniverse()
def test_1_argument() -> None:
    """only one argument supplied when 3 expected"""
    bad_neighbourhoods = [
        *BAD_TYPE_NOT_ITERABLE,
        *BAD_NEIGHBOURHOOD_TOO_FEW]
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
    """two arguments supplied when 3 expected"""
    assert len(BAD_TYPE_NOT_ITERABLE) > 0
    assert len(BAD_COUNT_RULES_NOT_UNIQUE) > 0
    assert len(BAD_COUNT_OUTSIDE_RANGE) > 0
    assert len(BAD_COUNT_NOT_INTEGER) > 0
    assert len(BAD_BIRTH_COUNTS) > 0
    assert len(GOOD_2D_SURVIVAL_COUNTS) > 0
    assert len(GOOD_2D_BIRTH_COUNTS) > 0
    bad_survivals = [
        *BAD_TYPE_NOT_ITERABLE,
        *BAD_COUNT_RULES_NOT_UNIQUE,
        (case for (case, _coord) in BAD_COUNT_OUTSIDE_RANGE),
        (case for (case, _err_type) in BAD_COUNT_NOT_INTEGER),
        *BAD_BIRTH_COUNTS,
        *GOOD_2D_SURVIVAL_COUNTS,
        *GOOD_2D_BIRTH_COUNTS,
    ]
    for survival in bad_survivals:
        with pytest.raises(TypeError, match=required_positional_message(1)):
            AutomataUniverse(NEIGHBOURHOOD_2D, survival)
def test_4_arguments() -> None:
    """four arguments supplied when 3 expected"""
    for extra in BAD_EXTRA_ARGUMENTS:
        with pytest.raises(TypeError, match=re.compile(
            "__init__[(][)] takes 4 positional arguments but {} were given".format(
                len(extra) + 4))):
            AutomataUniverse(NEIGHBOURHOOD_2D, [2,3], [3], *extra)
def test_non_iterable_neighbourhood() -> None:
    """neighbourhood parameter is not iterable"""
    assert len(BAD_TYPE_NOT_ITERABLE) > 0
    for neighbourhood in BAD_TYPE_NOT_ITERABLE:
        expected = "'{}' object is not iterable".format(type(neighbourhood).__name__)
        with pytest.raises(TypeError, match=re.compile(expected)):
            AutomataUniverse(neighbourhood, GOOD_2D_SURVIVAL_COUNTS[0], GOOD_2D_BIRTH_COUNTS[0])
def test_too_few_neighbours() -> None:
    """neighbourhood parameter is iterable but too short"""
    assert len(BAD_NEIGHBOURHOOD_TOO_FEW) > 0
    for neighbourhood in BAD_NEIGHBOURHOOD_TOO_FEW:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_2D_SURVIVAL_COUNTS[0], GOOD_2D_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (len(neighbourhood),
            "neighbourhood does not contain enough neighbours: at least 2 required")
def test_neighbours_not_unique() -> None:
    """iterable neighbourhood contains duplicate elements"""
    assert len(BAD_NEIGHBOURHOOD_NOT_UNIQUE) > 0
    for neighbourhood in BAD_NEIGHBOURHOOD_NOT_UNIQUE:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_2D_SURVIVAL_COUNTS[0], GOOD_2D_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0][0] == len(neighbourhood)
        assert excinfo.value.args[0][1] < len(neighbourhood)
        assert excinfo.value.args[0][2] == "Neighbourhood addresses are not unique"
def test_neighbour_element_not_tuple() -> None:
    """iterable neighbourhood contains non-tuple element"""
    assert len(BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE) > 0
    for (neighbourhood, err_type) in BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(TypeError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_2D_SURVIVAL_COUNTS[0], GOOD_2D_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type, "element of neighbourhood is not a tuple")
def do_addr_test(nei, err_type) -> None:
    """see if can get report showing neighbourhood set"""
    assert isinstance(nei, Iterable)
    with pytest.raises(TypeError) as excinfo:
        AutomataUniverse(nei, GOOD_2D_SURVIVAL_COUNTS[0], GOOD_2D_BIRTH_COUNTS[0])
    verify_general_exception_tuple(excinfo, 2)
    assert excinfo.value.args[0] == (err_type, "automata universe address is not a tuple")
def test_neighbour_address_not_tuple() -> None:
    """iterable neighbourhood contains non-tuple address"""
    assert len(BAD_NEIGHBOURHOOD_ADDR_NOT_TUPLE) > 0
    for (neighbourhood, err_type) in BAD_NEIGHBOURHOOD_ADDR_NOT_TUPLE:
        do_addr_test(neighbourhood, err_type)
def test_bad_neighbour_dimensions() -> None:
    """iterable neighbourhood address element different dimensions"""
    assert len(BAD_NEIGHBOURHOOD_DIMENSIONS) > 0
    for (neighbourhood, dimensions_uni, dimensions_addr) in BAD_NEIGHBOURHOOD_DIMENSIONS:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_2D_SURVIVAL_COUNTS[0], GOOD_2D_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0] == (dimensions_uni, dimensions_addr, "automata universe"
            " address does not have the same number of dimensions as the universe")
def test_coordinate_not_integer() -> None:
    """iterable neighbourhood contains non-integer coordinate"""
    assert len(BAD_NEIGHBOURHOOD_COORDINATE) > 0
    for (neighbourhood, err_type) in BAD_NEIGHBOURHOOD_COORDINATE:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(TypeError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_2D_SURVIVAL_COUNTS[0], GOOD_2D_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type,
            "automata universe address coordinate is not an integer")
def test_origin_in_neighbourhood() -> None:
    """iterable neighbourhood contains universe origin"""
    assert len(BAD_NEIGHBOURHOOD_WITH_ORIGIN) > 0
    for neighbourhood in BAD_NEIGHBOURHOOD_WITH_ORIGIN:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_2D_SURVIVAL_COUNTS[0], GOOD_2D_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (tuple([0]*len(neighbourhood[0])),
            "the universe origin is not a valid neighbourhood address")
def test_no_symmetric_address() -> None:
    """iterable neighbourhood missing symmetric address"""
    assert len(BAD_NEIGHBOURHOOD_WITHOUT_SYMMETRIC) > 0
    for (neighbourhood, address) in BAD_NEIGHBOURHOOD_WITHOUT_SYMMETRIC:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(neighbourhood, GOOD_2D_SURVIVAL_COUNTS[0], GOOD_2D_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (address, "no symmetric address in neighbourhood")
def test_no_coordinate() -> None:
    """iterable neighbourhood first address without coordinate"""
    assert len(BAD_NEIGHBOURHOOD_NO_COORDINATE) > 0
    for neighbourhood in BAD_NEIGHBOURHOOD_NO_COORDINATE:
        assert isinstance(neighbourhood, Iterable)
        with pytest.raises(TypeError, match=re.compile(
                "neighbourhood address does not contain at least 1 coordinate")):
            AutomataUniverse(neighbourhood, GOOD_2D_SURVIVAL_COUNTS[0], GOOD_2D_BIRTH_COUNTS[0])

def test_not_iterable_survial() -> None:
    """survival rule parameter is not iterable"""
    assert len(BAD_TYPE_NOT_ITERABLE) > 0
    for rule in BAD_TYPE_NOT_ITERABLE:
        expected = "'{}' object is not iterable".format(type(rule).__name__)
        with pytest.raises(TypeError, match=re.compile(expected)):
            AutomataUniverse(NEIGHBOURHOOD_2D, rule, GOOD_2D_BIRTH_COUNTS[0])
def test_not_iterable_birth() -> None:
    """birth rule parameter is not iterable"""
    assert len(BAD_TYPE_NOT_ITERABLE) > 0
    for rule in BAD_TYPE_NOT_ITERABLE:
        expected = "'{}' object is not iterable".format(type(rule).__name__)
        with pytest.raises(TypeError, match=re.compile(expected)):
            AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_2D_SURVIVAL_COUNTS[0], rule)
def test_not_unique_survial() -> None:
    """survival rule parameter contains duplicate entry"""
    assert len(BAD_COUNT_RULES_NOT_UNIQUE) > 0
    for rule in BAD_COUNT_RULES_NOT_UNIQUE:
        assert isinstance(rule, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, rule, GOOD_2D_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0][0] == len(rule)
        assert excinfo.value.args[0][1] < len(rule)
        assert excinfo.value.args[0][2] == "cell propagation rule counts are not unique"
def test_not_unique_birth() -> None:
    """birth rule parameter contains duplicate entry"""
    assert len(BAD_COUNT_RULES_NOT_UNIQUE) > 0
    for rule in BAD_COUNT_RULES_NOT_UNIQUE:
        assert isinstance(rule, Iterable)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_2D_SURVIVAL_COUNTS[0], rule)
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0][0] == len(rule)
        assert excinfo.value.args[0][1] < len(rule)
        assert excinfo.value.args[0][2] == "cell propagation rule counts are not unique"
def test_survival_outside_neighbourhood() -> None:
    """survival rule parameter entry outside neighbourhood size"""
    assert len(BAD_COUNT_OUTSIDE_RANGE) > 0
    for (rule, outside) in BAD_COUNT_OUTSIDE_RANGE:
        assert isinstance(rule, Iterable)
        assert isinstance(outside, int)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, rule, GOOD_2D_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0][0] == outside
        assert excinfo.value.args[0][1] == len(NEIGHBOURHOOD_2D)
        assert excinfo.value.args[0][2] == \
            "cell propagation count case is not between 0 and neighbourhood size"
def test_birth_outside_neighbourhood() -> None:
    """birth rule parameter entry outside neighbourhood size"""
    assert len(BAD_COUNT_OUTSIDE_RANGE) > 0
    for (rule, outside) in BAD_COUNT_OUTSIDE_RANGE:
        assert isinstance(rule, Iterable)
        assert isinstance(outside, int)
        with pytest.raises(ValueError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_2D_SURVIVAL_COUNTS[0], rule)
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0][0] == outside
        assert excinfo.value.args[0][1] == len(NEIGHBOURHOOD_2D)
        assert excinfo.value.args[0][2] == \
            "cell propagation count case is not between 0 and neighbourhood size"
def test_survival_not_integer() -> None:
    """survival rule parameter entry not an integer"""
    assert len(BAD_COUNT_NOT_INTEGER) > 0
    for (rule, err_type) in BAD_COUNT_NOT_INTEGER:
        assert isinstance(rule, Iterable)
        with pytest.raises(TypeError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, rule, GOOD_2D_BIRTH_COUNTS[0])
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type,
            "cell propagation count case is not an integer")
def test_birth_not_integer() -> None:
    """birth rule parameter entry not an integer"""
    assert len(BAD_COUNT_NOT_INTEGER) > 0
    for (rule, err_type) in BAD_COUNT_NOT_INTEGER:
        assert isinstance(rule, Iterable)
        with pytest.raises(TypeError) as excinfo:
            AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_2D_SURVIVAL_COUNTS[0], rule)
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type,
            "cell propagation count case is not an integer")
def test_zero_birth_rule() -> None:
    """birth rule contains entry for zero neighbours"""
    assert len(BAD_BIRTH_COUNTS) > 0
    for rule in BAD_BIRTH_COUNTS:
        assert isinstance(rule, Iterable)
        with pytest.raises(ValueError, match=re.compile(
                "zero is not a valid birth propagation rule value")):
            AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_2D_SURVIVAL_COUNTS[0], rule)

def test_properties() -> None:
    """verify the properties of a successful instance creation"""
    test_neighbourhood = NEIGHBOURHOOD_2D
    test_survival = GOOD_2D_SURVIVAL_COUNTS[0]
    test_birth = GOOD_2D_BIRTH_COUNTS[0]
    uni = AutomataUniverse(test_neighbourhood, test_survival, test_birth)
    assert isinstance(uni, AutomataUniverse)
    assert isinstance(uni.dimensions, int)
    assert uni.dimensions == 2 # previously caused pylint crash/trace
    assert isinstance(uni.neighbourhood, frozenset)
    assert uni.neighbourhood == frozenset(test_neighbourhood)
    assert isinstance(uni.neighbourhood_population, int)
    assert uni.neighbourhood_population == len(test_neighbourhood)
    assert isinstance(uni.survival_rules, frozenset)
    assert uni.survival_rules == frozenset(test_survival)
    assert isinstance(uni.birth_rules, frozenset)
    assert uni.birth_rules == frozenset(test_birth)
    assert uni.identity_matrix == IDENTITY_MATRIX_2D
    assert isinstance(hash(uni), int)
    assert hasattr(uni, "__hash__")
def test_identity() -> None:
    """verify identity matrix for different universe dimensions"""
    uni = base_universe_instance_1d()
    assert uni.identity_matrix == IDENTITY_MATRIX_1D
    uni = base_universe_instance_2d()
    assert uni.identity_matrix == IDENTITY_MATRIX_2D
    uni = base_universe_instance_3d()
    assert uni.identity_matrix == IDENTITY_MATRIX_3D

def test_validate_address_not_tuple() -> None:
    """address parameter is not a tuple"""
    uni = base_universe_instance_2d()
    assert len(BAD_ADDRESS_NOT_TUPLE) > 0
    for (address, err_type) in BAD_ADDRESS_NOT_TUPLE:
        with pytest.raises(TypeError) as excinfo:
            uni.validate_address(address)
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type, "automata universe address is not a tuple")
def test_validate_address_bad_dimensions() -> None:
    """address parameter coordinates do not match universe dimensions"""
    uni = base_universe_instance_2d()
    assert len(BAD_ADDRESS_NOT_2_DIMENSIONS) > 0
    for (address, dim) in BAD_ADDRESS_NOT_2_DIMENSIONS:
        with pytest.raises(ValueError) as excinfo:
            uni.validate_address(address)
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0] == (uni.dimensions, dim, "automata universe address "
            "does not have the same number of dimensions as the universe")
def test_validate_address_not_integer() -> None:
    """coordinate in address parameter is not an integer"""
    uni = base_universe_instance_2d()
    assert len(BAD_ADDRESS_COORD_NOT_INTEGER) > 0
    for (address, err_type) in BAD_ADDRESS_COORD_NOT_INTEGER:
        with pytest.raises(TypeError) as excinfo:
            uni.validate_address(address)
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type,
            "automata universe address coordinate is not an integer")
def test_validate_good_address() -> None:
    """good 2D address"""
    uni = base_universe_instance_2d()
    assert len(GOOD_2D_ADDRESS) > 0
    for address in GOOD_2D_ADDRESS:
        uni.validate_address(address)

def test_universe_address_false() -> None:
    """parameter is not a valid (2D) universe address"""
    uni = base_universe_instance_2d()
    bad_address = [case for (case, _err) in BAD_ADDRESS_NOT_TUPLE]
    for (case, _err) in BAD_ADDRESS_NOT_2_DIMENSIONS:
        bad_address.append(case)
    for (case, _err) in BAD_ADDRESS_COORD_NOT_INTEGER:
        bad_address.append(case)
    assert len(bad_address) > 0
    for address in bad_address:
        assert uni.is_universe_address(address) is False
def test_universe_address_true() -> None:
    """parameter is a valid (2D) universe address"""
    uni = base_universe_instance_2d()
    assert len(GOOD_2D_ADDRESS) > 0
    for address in GOOD_2D_ADDRESS:
        assert uni.is_universe_address(address) is True

def test_neighbours_bad_address() -> None:
    """cell address is not valid for universe"""
    uni = base_universe_instance_2d()
    bad_addresses = [address for (address, _error_details) in BAD_ADDRESS_NOT_TUPLE]
    bad_addresses.extend([address for (address, _error_details) in BAD_ADDRESS_COORD_NOT_INTEGER])
    assert len(bad_addresses) > 0
    for address in bad_addresses:
        # simple check: details verified with specific tests for .validate_address
        with pytest.raises(TypeError):
            uni.neighbours(address)
    for (address, _err_details) in BAD_ADDRESS_NOT_2_DIMENSIONS:
        # simple check: details verified with specific tests for .validate_address
        with pytest.raises(ValueError):
            uni.neighbours(address)
def test_neighbours_result() -> None:
    """verify result with good cell address tuple"""
    uni = base_universe_instance_2d()
    assert len(GOOD_CELL_ADDRESS_NEIGHBOURS) > 0
    for (address, neighbours) in GOOD_CELL_ADDRESS_NEIGHBOURS:
        assert uni.neighbours(address) == neighbours

def test_validate_matrix_not_iterable() -> None:
    """matrix parameter is not iterable"""
    uni = base_universe_instance_2d()
    assert len(BAD_TYPE_NOT_ITERABLE) > 0
    for matrix in BAD_TYPE_NOT_ITERABLE:
        expected = "'{}' object is not iterable".format(type(matrix).__name__)
        with pytest.raises(TypeError, match=re.compile(expected)):
            uni.validate_matrix(matrix)
def test_validate_matrix_not_tuple() -> None:
    """matrix parameter contains address that is not a tuple"""
    uni = base_universe_instance_2d()
    assert len(BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE) > 0
    for (matrix, err_type) in BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE:
        assert isinstance(matrix, Iterable)
        with pytest.raises(TypeError) as excinfo:
            uni.validate_matrix(matrix)
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type, "automata universe address is not a tuple")
def test_validate_matrix_bad_vector_dimensions() -> None:
    """matrix contains address with dimensions different from universe"""
    uni = base_universe_instance_2d()
    assert len(BAD_NEIGHBOURHOOD_DIMENSIONS) > 0
    for (matrix, u_dim, a_dim) in BAD_NEIGHBOURHOOD_DIMENSIONS:
        assert u_dim == uni.dimensions
        with pytest.raises(ValueError) as excinfo:
            assert isinstance(matrix, Iterable)
            uni.validate_matrix(matrix)
        verify_general_exception_tuple(excinfo, 3)
        assert excinfo.value.args[0] == (uni.dimensions, a_dim, "automata universe address "
            "does not have the same number of dimensions as the universe")
def test_validate_matrix_not_integer() -> None:
    """matrix contains address with non-integer coordinate"""
    uni = base_universe_instance_2d()
    assert len(BAD_NEIGHBOURHOOD_COORDINATE) > 0
    for (matrix, err_type) in BAD_NEIGHBOURHOOD_COORDINATE:
        assert isinstance(matrix, Iterable)
        with pytest.raises(TypeError) as excinfo:
            uni.validate_matrix(matrix)
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type,
            "automata universe address coordinate is not an integer")
def test_validate_matrix_bad_matrix_dimensions() -> None:
    """matrix contains number of vectors different from universe dimensions"""
    uni = base_universe_instance_2d()
    assert len(BAD_2D_MATRIX_VECTOR_COUNT) > 0
    for matrix in BAD_2D_MATRIX_VECTOR_COUNT:
        expected = "matrix contains {} vectors, which does not match {} dimensions " \
            "for the configured universe".format(len(matrix), uni.dimensions)
        with pytest.raises(TypeError, match=re.compile(expected)):
            uni.validate_matrix(matrix)
def test_validate_matrix_good_data() -> None:
    """matrix is valid for the current (2D) universe"""
    uni = base_universe_instance_2d()
    assert len(GOOD_2D_MATRIX) > 0
    for matrix in GOOD_2D_MATRIX:
        uni.validate_matrix(matrix)

def test_vector_dot_product() -> None:
    """verify 2d dot product vector cases"""
    uni = base_universe_instance_2d()
    assert len(DOT_PRODUCT_2D_VECTOR_TESTS) > 0
    for (vector, matrix, expected) in DOT_PRODUCT_2D_VECTOR_TESTS:
        assert uni.vector_dot_product(vector, matrix) == expected
def test_square_dot_product() -> None:
    """verify square 2d dot product cases"""
    uni = base_universe_instance_2d()
    assert len(DOT_PRODUCT_2D_SQUARE_EXPECTED) > 0
    for (transform, base, expected) in DOT_PRODUCT_2D_SQUARE_EXPECTED:
        assert uni.matrix_transform(transform, base) == expected

def test_check_cell_group_not_set() -> None:
    """verify exception raised if cell group is not the expected datatype"""
    uni = base_universe_instance_2d()
    assert len(BAD_POPULATION_NOT_SET) > 0
    for (cells, err_type) in BAD_POPULATION_NOT_SET:
        with pytest.raises(TypeError) as excinfo:
            uni._check_cell_group(cells)
        verify_general_exception_tuple(excinfo, 2)
        assert excinfo.value.args[0] == (err_type, "automata cell group is not a set")
    with pytest.raises(TypeError) as excinfo:
        uni._check_cell_group(uni)
    verify_general_exception_tuple(excinfo, 2)
    assert excinfo.value.args[0] == (type(uni), "automata cell group is not a set")
def test_check_cell_group_not_tuple() -> None:
    """verify exception raised if cell in group is not a tuple"""
    uni = base_universe_instance_2d()
    bad_populations = (
        *(((case,)) for (case, _err_details) in BAD_ADDRESS_NOT_TUPLE),
        *(((case,)) for (case, _err_details) in BAD_ADDRESS_NOT_2_DIMENSIONS),
        *(((case,)) for (case, _err_details) in BAD_ADDRESS_COORD_NOT_INTEGER),
    )
    assert len(bad_populations) > 0
    # assert len(bad_populations) == 7 + 4 + 8
    for cells in bad_populations:
        # simple check: details verified with specific tests for validate_address
        with pytest.raises(TypeError):
            uni._check_cell_group(cells)
def test_check_cell_group_good() -> None:
    """verify success with good cell group"""
    uni = base_universe_instance_2d()
    assert len(GOOD_2D_POPULATION) > 0
    for cells in GOOD_2D_POPULATION:
        assert uni._check_cell_group(cells) is None

def test_group_translate_bad_arguments() -> None:
    """verify exception raised if cell group or vector is not valid for universe"""
    uni = base_universe_instance_2d()
    bad_populations = (
        *(((case,)) for (case, _err_details) in BAD_ADDRESS_NOT_TUPLE),
        *(((case,)) for (case, _err_details) in BAD_ADDRESS_NOT_2_DIMENSIONS),
        *(((case,)) for (case, _err_details) in BAD_ADDRESS_COORD_NOT_INTEGER),
        *(((case,)) for (case, _err_details) in BAD_POPULATION_NOT_SET),
    )
    assert len(bad_populations) > 0
    # assert len(bad_populations) == 7 + 4 + 8 + 8
    for cells in bad_populations:
        # simple check: details verified with specific tests for _check_cell_group
        with pytest.raises(TypeError):
            uni.cell_group_translate(cells, ORIGIN_2D)
    bad_offsets = [address for (address, _error_details) in BAD_ADDRESS_NOT_TUPLE]
    bad_offsets.extend([address for (address, _error_details) in BAD_ADDRESS_COORD_NOT_INTEGER])
    assert len(bad_offsets) > 0
    good_matrix = frozenset(IDENTITY_MATRIX_2D)
    for vector in bad_offsets:
        # simple check: details verified with specific tests for .validate_address
        with pytest.raises(TypeError):
            uni.cell_group_translate(good_matrix, vector)
    for (vector, _err_details) in BAD_ADDRESS_NOT_2_DIMENSIONS:
        # simple check: details verified with specific tests for .validate_address
        with pytest.raises(ValueError):
            uni.cell_group_translate(good_matrix, vector)
def test_group_translate_result() -> None:
    """verify result with good matrix and offset vector"""
    uni = base_universe_instance_2d()
    assert len(GOOD_CELL_GROUP_TRANSLATE) > 0
    # assert len(GOOD_CELL_GROUP_TRANSLATE) == 4 * 1 + 2
    for (cells, offset, result) in GOOD_CELL_GROUP_TRANSLATE:
        assert uni.cell_group_translate(cells, offset) == result

def test_group_transform_bad_arguments() -> None:
    """verify exception raised if cell group or matrix is not valid for universe"""
    uni = base_universe_instance_2d()
    bad_transform = (
        *(case for case in BAD_TYPE_NOT_ITERABLE),
        *(case for (case, _err_details) in BAD_NEIGHBOURHOOD_ELE_NOT_TUPLE),
        *(case for (case, _err_details) in BAD_NEIGHBOURHOOD_COORDINATE),
        *(case for case in BAD_2D_MATRIX_VECTOR_COUNT),
    )
    assert len(bad_transform) > 0
    # for matrix in BAD_TYPE_NOT_ITERABLE:
    for matrix in bad_transform:
        # simple check: details verified with specific tests for .validate_matrix
        with pytest.raises(TypeError):
            uni.cell_group_transform(GOOD_2D_POPULATION[0], matrix)
    assert len(BAD_NEIGHBOURHOOD_DIMENSIONS) > 0
    for (matrix, _dim, _err_details) in BAD_NEIGHBOURHOOD_DIMENSIONS:
        with pytest.raises(ValueError):
            uni.cell_group_transform(GOOD_2D_POPULATION[0], matrix)
    bad_populations = (
        *(((case,)) for (case, _err_details) in BAD_ADDRESS_NOT_TUPLE),
        *(((case,)) for (case, _err_details) in BAD_ADDRESS_NOT_2_DIMENSIONS),
        *(((case,)) for (case, _err_details) in BAD_ADDRESS_COORD_NOT_INTEGER),
        *(((case,)) for (case, _err_details) in BAD_POPULATION_NOT_SET),
    )
    assert len(bad_populations) > 0
    # assert len(bad_populations) == 7 + 4 + 8 + 8
    for cells in bad_populations:
        # simple check: details verified with specific tests for _check_cell_group
        with pytest.raises(TypeError):
            uni.cell_group_transform(cells, IDENTITY_MATRIX_2D)
def test_group_transform_result() -> None:
    """verify result with good cell data and transformation matrix"""
    uni = base_universe_instance_2d()
    assert len(GOOD_CELL_GROUP_TRANSFORM) > 0
    # assert len(GOOD_CELL_GROUP_TRANSFORM) == 4 * 1 + ?
    for (cells, matrix, result) in GOOD_CELL_GROUP_TRANSFORM:
        assert uni.cell_group_transform(cells, matrix) == result

def test_rotate_matrix_non_iterable_matrix() -> None:
    """matrix test parameter is not iterable"""
    uni = base_universe_instance_2d()
    assert len(BAD_TYPE_NOT_ITERABLE) > 0
    for matrix in BAD_TYPE_NOT_ITERABLE:
        expected = "'{}' object is not iterable".format(type(matrix).__name__)
        with pytest.raises(TypeError, match=re.compile(expected)):
            uni.is_rotation_matrix(matrix)
# wrong dimensions for universe «validate_address»
# not square ?? «validate_matrix» ((1,2),)
def test_is_rotation_matrix_truth() -> None:
    """verify result with good matrix values"""
    uni = base_universe_instance_2d()
    for case in NOT_ROTATION_MATRIX:
        assert uni.is_rotation_matrix(case) is False
    for case in GOOD_ROTATION_MATRIX:
        assert uni.is_rotation_matrix(case) is True

def test_step_result() -> None:
    """verify result with good cells"""
    uni = base_universe_instance_2d()
    for (cells, expected) in GOOD_STEP_POPULATIONS:
        assert uni.step(cells) == expected

# def test_exploration() -> None:
    # """view exact failure for test case"""
    # # AutomataUniverse(NEIGHBOURHOOD_2D, ((-1,)), GOOD_BIRTH_COUNTS[0])
    # # AutomataUniverse(NEIGHBOURHOOD_2D, BAD_[0], GOOD_BIRTH_COUNTS[0])
    # # AutomataUniverse(NEIGHBOURHOOD_2D, GOOD_SURVIVAL_COUNTS[0], BAD_COUNT_RULES_NOT_UNIQUE[0])
    # uni = base_2d_universe_instance()
    # test_mat = (
    #     (1,2,0),(2,2.-1),(3,1,-2)
    # )
    # assert uni.is_rotation_matrix(test_mat) is True

    # isinstance(AutomataUniverse(), AutomataUniverse)
