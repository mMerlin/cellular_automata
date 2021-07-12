#!/usr/bin/env python
# coding=utf-8
# pylint: disable=R0801

"""
regression tests for creation of cellular automata universe instances
"""

import re
# from _pytest._code.code import ExceptionInfo
import pytest
from common_test_data import (
    ORIGIN_1D,
    ORIGIN_2D,
    ORIGIN_3D,
    IDENTITY_MATRIX_1D,
    IDENTITY_MATRIX_2D,
    IDENTITY_MATRIX_3D,
    UNIT_VECTOR_1D,
    UNIT_VECTOR_2D,
    UNIT_VECTOR_3D,
    ZERO_MATRIX_1D,
    ZERO_MATRIX_2D,
    ZERO_MATRIX_3D,
    REFLECT_MATRIX_1D,
    ROTATE_MATRIX_2D_90,
    ROTATE_MATRIX_2D_180,
    # ROTATE_MATRIX_2D_270,
    REFLECTION_MATRIX_2D_HORIZONTAL,
    REFLECTION_MATRIX_2D_VERTICAL,
    TRANPOSE_MATRIX_2D,
    SCALE_3X_5Y_MATRIX,
    NOT_INTEGER_TYPE_SAMPLES,
    NON_POSITIVE_INTEGERS,
    SQUARE_MATRIX_TRANSFORM_OPERATIONS_2D,
)
from math_tools import (identity_matrix, vector_dot_product, matrix_transpose, matrix_transform,
    matrix_determinant, matrix_determinant2)

# remove need to add parent to path?
# `pipenv run python -m pytest «»`
# ? -q -v

GOOD_IDENTITY_MATRIX = (
    IDENTITY_MATRIX_1D,
    IDENTITY_MATRIX_2D,
    IDENTITY_MATRIX_3D,
    ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)),
)

Q1_5_11 = (5,11)
Q2_5_11 = (-5,11)
Q3_5_11 = (-5,-11)
Q4_5_11 = (5,-11)
Q1_11_5 = (11,5)
Q2_11_5 = (-11,5)
Q3_11_5 = (-11,-5)
Q4_11_5 = (11,-5)
VECTOR_DOT_PRODUCT_OPERATIONS = (
    # («vector», «matrix», «expected_result»),
    (ORIGIN_1D, IDENTITY_MATRIX_1D, ORIGIN_1D),
    (UNIT_VECTOR_1D, IDENTITY_MATRIX_1D, UNIT_VECTOR_1D),
    ((9,), ZERO_MATRIX_1D, ORIGIN_1D),
    ((-19,), ZERO_MATRIX_1D, ORIGIN_1D),
    (UNIT_VECTOR_1D, ZERO_MATRIX_1D, ORIGIN_1D),
    (UNIT_VECTOR_1D, REFLECT_MATRIX_1D, ((-1),)),
    ((-19,), REFLECT_MATRIX_1D, ((19),)),

    (ORIGIN_2D, IDENTITY_MATRIX_2D, ORIGIN_2D),
    (UNIT_VECTOR_2D, IDENTITY_MATRIX_2D, UNIT_VECTOR_2D),
    (UNIT_VECTOR_2D, IDENTITY_MATRIX_2D, UNIT_VECTOR_2D),
    (UNIT_VECTOR_2D, ZERO_MATRIX_2D, ORIGIN_2D),
    (Q1_11_5, ZERO_MATRIX_2D, ORIGIN_2D),
    (Q1_5_11, ZERO_MATRIX_2D, ORIGIN_2D),
    (Q2_5_11, ZERO_MATRIX_2D, ORIGIN_2D),
    (Q1_5_11, IDENTITY_MATRIX_2D, Q1_5_11),
    (Q2_5_11, IDENTITY_MATRIX_2D, Q2_5_11),
    (Q3_5_11, IDENTITY_MATRIX_2D, Q3_5_11),
    (Q4_5_11, IDENTITY_MATRIX_2D, Q4_5_11),
    ((1.2, 2.5), IDENTITY_MATRIX_2D, (1.2, 2.5)),
    (Q1_5_11, ROTATE_MATRIX_2D_90, Q2_11_5),
    (Q2_5_11, ROTATE_MATRIX_2D_90, Q3_11_5),
    (Q3_5_11, ROTATE_MATRIX_2D_90, Q4_11_5),
    (Q4_5_11, ROTATE_MATRIX_2D_90, Q1_11_5),
    (Q1_11_5, ROTATE_MATRIX_2D_90, Q2_5_11),
    (Q2_11_5, ROTATE_MATRIX_2D_90, Q3_5_11),
    (Q3_11_5, ROTATE_MATRIX_2D_90, Q4_5_11),
    (Q4_11_5, ROTATE_MATRIX_2D_90, Q1_5_11),
    (Q1_5_11, REFLECTION_MATRIX_2D_HORIZONTAL, Q2_5_11),
    (Q2_5_11, REFLECTION_MATRIX_2D_HORIZONTAL, Q1_5_11),
    (Q3_5_11, REFLECTION_MATRIX_2D_HORIZONTAL, Q4_5_11),
    (Q4_5_11, REFLECTION_MATRIX_2D_HORIZONTAL, Q3_5_11),
    (Q1_5_11, REFLECTION_MATRIX_2D_VERTICAL, Q4_5_11),
    (Q2_5_11, REFLECTION_MATRIX_2D_VERTICAL, Q3_5_11),
    (Q3_5_11, REFLECTION_MATRIX_2D_VERTICAL, Q2_5_11),
    (Q4_5_11, REFLECTION_MATRIX_2D_VERTICAL, Q1_5_11),
    (UNIT_VECTOR_2D, SCALE_3X_5Y_MATRIX, (3,5)),
    (UNIT_VECTOR_2D, SCALE_3X_5Y_MATRIX, (3,5)),
    (ORIGIN_2D, SCALE_3X_5Y_MATRIX, ORIGIN_2D),

    (ORIGIN_3D, IDENTITY_MATRIX_3D, ORIGIN_3D),
    (UNIT_VECTOR_3D, IDENTITY_MATRIX_3D, UNIT_VECTOR_3D),
    ((1,2,3), IDENTITY_MATRIX_3D, (1,2,3)),
    (UNIT_VECTOR_3D, ZERO_MATRIX_3D, ORIGIN_3D),
    ((1,2,3), ZERO_MATRIX_3D, ORIGIN_3D),
)

GOOD_MATRIX_TRANSPOSE = (
    # («original», «transposed»),
    (((3,5),(7,11)), ((3,7),(5,11))),
    (((3,5),), ((3,),(5,))),
    (((3,),(7,)), ((3,7),)),
    (((3,),(7,),(11,)), ((3,7,11),)),
    (((3,5,13),(7,11,17)), ((3,7),(5,11),(13,17))),
    (((3,5),(7,11),(13,17)), ((3,7,13),(5,11,17))),
    ((("test", 5.5), ("step", 11.7)), (("test", "step"), (5.5, 11.7))),
    (((3,[5,6]),(7,11),([13,17],2)), ((3,7,[13,17]),([5,6],11,2))),
    (TRANPOSE_MATRIX_2D, TRANPOSE_MATRIX_2D),
)

MATRIX_TRANSFORM_AND_RESULTS = (
    # («matrix», «transform», «expected»),
    *SQUARE_MATRIX_TRANSFORM_OPERATIONS_2D,
    (IDENTITY_MATRIX_3D, IDENTITY_MATRIX_3D, IDENTITY_MATRIX_3D),
    (
        ((3,2,1,5),(9,1,3,0)),
        ((2,9,0),(1,3,5),(2,4,7),(8,1,5)),
        ((50,42,42),(25,96,26))
    ),
    (
        ((1,2,3),(4,5,6)),
        ((10,11),(20,21),(30,31)),
        ((140,146),(320,335))
    ),
)

MATRIX_AND_DETERMINANT = (
    # («matrix», «determinant»),
    (((3,4,-2),(3,5,0),(-1,4,1)), -31),
    (((-5,-4),(-2,-3)), 7),
    (((2,-3,1),(2,0,-1),(1,4,5)), 49),
    (((1,3,1,4),(3,9,5,15),(0,2,1,1),(0,4,2,3)), -4),
    (((-5,0,-1),(1,2,-1),(-3,4,1)), -40),
    (IDENTITY_MATRIX_1D, 1),
    (IDENTITY_MATRIX_2D, 1),
    (IDENTITY_MATRIX_3D, 1),
    (ROTATE_MATRIX_2D_90, 1),
    (ROTATE_MATRIX_2D_180, 1),
    (ZERO_MATRIX_1D, 0),
    (ZERO_MATRIX_2D, 0),
    (ZERO_MATRIX_3D, 0),
    (((8,-6,2),(-6,7,-4),(2,-4,3)), 0),
    (((1,2,3),(2,4,6),(0,2,7)), 0),
)

def test_identity_matrix_good() -> None:
    """verify results creating different identity matrices"""
    assert len(GOOD_IDENTITY_MATRIX) > 0
    for dim in range(1, len(GOOD_IDENTITY_MATRIX) + 1):
        assert identity_matrix(dim) == GOOD_IDENTITY_MATRIX[dim - 1]
def test_identity_matrix_not_integer() -> None:
    """verify exception raised for non-integer dimension"""
    assert len(NOT_INTEGER_TYPE_SAMPLES) > 0
    for dim in NOT_INTEGER_TYPE_SAMPLES:
        expected_error = "Identity Matrix dimensions must be an intger, not {}".format(type(dim))
        with pytest.raises(TypeError, match=re.compile(expected_error)):
            identity_matrix(dim)
def test_identity_matrix_not_positive() -> None:
    """verify exception raised for invalid dimension value"""
    assert len(NON_POSITIVE_INTEGERS) > 0
    for dim in NON_POSITIVE_INTEGERS:
        expected_error = "{} is not a valid Identity Matrix dimension " \
            "must be greater than zero".format(dim)
        with pytest.raises(ValueError, match=re.compile(expected_error)):
            identity_matrix(dim)

def test_vector_dot_product_results() -> None:
    """verify results of vector dot product with good inputs"""
    assert len(VECTOR_DOT_PRODUCT_OPERATIONS) > 0
    for (vector, matrix, expected) in VECTOR_DOT_PRODUCT_OPERATIONS:
        assert vector_dot_product(vector, matrix) == expected

def test_transpose_matrix_results() -> None:
    """verify mxn matrix transpose cases"""
    assert len(GOOD_MATRIX_TRANSPOSE) > 0
    for (matrix, expected) in GOOD_MATRIX_TRANSPOSE:
        assert matrix_transpose(matrix) == expected

def test_matrix_transform_results() -> None:
    """verify matrix tranform results"""
    assert len(MATRIX_TRANSFORM_AND_RESULTS) > 0
    for (matrix, transform, expected) in MATRIX_TRANSFORM_AND_RESULTS:
        assert matrix_transform(matrix, transform) == expected

def test_matrix_determinant_results() -> None:
    """verify matrix tranform results"""
    assert len(MATRIX_AND_DETERMINANT) > 0
    for (matrix, determinant) in MATRIX_AND_DETERMINANT:
        assert matrix_determinant(matrix) == determinant

def test_matrix_determinant2_results() -> None:
    """verify matrix tranform results"""
    assert len(MATRIX_AND_DETERMINANT) > 0
    for (matrix, determinant) in MATRIX_AND_DETERMINANT:
        assert matrix_determinant2(matrix) == determinant
