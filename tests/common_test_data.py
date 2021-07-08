#!/usr/bin/env python
# coding=utf-8
# py lint: disable=

'''
data that can be used for multiple test conditions
'''

# from typing import Iterable
# import re
# from _pytest._code.code import ExceptionInfo
# import pytest
# from automata_universe import AutomataUniverse

NEIGHBOURHOOD_1D = ((-1,),(1,))
NEIGHBOURHOOD_2D = (
    (-1,-1), (-1,0), (-1,1),
    ( 0,-1),         ( 0,1),
    ( 1,-1), ( 1,0), ( 1,1),
)
NEIGHBOURHOOD_3D = (
    (-1,-1,-1), (-1,-1,0), (-1,-1, 1),
    (-1, 0,-1), (-1, 0,0), (-1, 0, 1),
    (-1, 1,-1), (-1, 1,0), (-1, 1, 1),

    (0,-1,-1), (0,-1,0), (0,-1, 1),
    (0, 0,-1),           (0, 0, 1),
    (0, 1,-1), (0, 1,0), (0, 1, 1),

    (1,-1,-1), (1,-1,0), (1,-1, 1),
    (1, 0,-1), (1, 0,0), (1, 0, 1),
    (1, 1,-1), (1, 1,0), (1, 1, 1),
)
ORIGIN_1D = (0,)
ORIGIN_2D = (0,0)
ORIGIN_3D = (0,0,0)
UNIT_VECTOR_1D = (1,)
UNIT_VECTOR_2D = (1,1)
UNIT_VECTOR_3D = (1,1,1)
IDENTITY_MATRIX_1D = ((1,),)
IDENTITY_MATRIX_2D = ((1, 0), (0, 1))
IDENTITY_MATRIX_3D = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
ZERO_MATRIX_1D = ((0,),)
ZERO_MATRIX_2D = ((0, 0), (0, 0))
ZERO_MATRIX_3D = ((0, 0, 0), (0, 0, 0), (0, 0, 0))

REFLECT_MATRIX_1D = ((-1,),)
ROTATE_MATRIX_2D_90 = ((0, -1), (1, 0))
ROTATE_MATRIX_2D_180 = ((-1, 0), (0, -1))
ROTATE_MATRIX_2D_270 = ((0, 1), (-1, 0))
REFLECTION_MATRIX_2D_HORIZONTAL = ((-1, 0), (0, 1))
REFLECTION_MATRIX_2D_VERTICAL = ((1, 0), (0, -1))

# SCALE_7X = ((7,),)
SCALE_3X_5Y_MATRIX = ((3, 0), (0, 5))

NOT_INTEGER_TYPE_SAMPLES = (
    None,
    1.2,
    1/3,
    [],
    set(),
    frozenset(),
    dict(),
    (1,),
    "test",
)
NON_POSITIVE_INTEGERS = (
    0,
    -1,
    -99,
)

SQUARE_MATRIX_TRANSFORM_OPERATIONS_2D = (
    (IDENTITY_MATRIX_2D, IDENTITY_MATRIX_2D, IDENTITY_MATRIX_2D),
    (ROTATE_MATRIX_2D_90, ROTATE_MATRIX_2D_270, IDENTITY_MATRIX_2D),
    (ROTATE_MATRIX_2D_270, ROTATE_MATRIX_2D_90, IDENTITY_MATRIX_2D),
    (ROTATE_MATRIX_2D_180, ROTATE_MATRIX_2D_180, IDENTITY_MATRIX_2D),
    (IDENTITY_MATRIX_2D, ROTATE_MATRIX_2D_90, ROTATE_MATRIX_2D_90),
    (ROTATE_MATRIX_2D_90, IDENTITY_MATRIX_2D, ROTATE_MATRIX_2D_90),
    (ROTATE_MATRIX_2D_90, ROTATE_MATRIX_2D_90, ROTATE_MATRIX_2D_180),
    (ROTATE_MATRIX_2D_180, ROTATE_MATRIX_2D_90, ROTATE_MATRIX_2D_270),
    (
        ((3,  5),  # a=3 b=5
         (7, 11)), # c=7 d=11
        (          (13,   17),  # e=13 f=17
                   (29,   23)), # g=19 h=23
        ((3*13+ 5*29, 3*17+ 5*23), # ae + bg, af + bh
         (7*13+11*29, 7*17+11*23)) # ce + dg, cf + dh
    ),
)

GOOD_NEIGHBOURHOOD = (
    NEIGHBOURHOOD_1D,
    NEIGHBOURHOOD_2D,
    NEIGHBOURHOOD_3D,
    list(NEIGHBOURHOOD_1D),
    list(NEIGHBOURHOOD_2D),
    list(NEIGHBOURHOOD_3D),
    set(NEIGHBOURHOOD_1D),
    set(NEIGHBOURHOOD_2D),
    set(NEIGHBOURHOOD_3D),
    frozenset(NEIGHBOURHOOD_1D),
    frozenset(NEIGHBOURHOOD_2D),
    frozenset(NEIGHBOURHOOD_3D),
)
BAD_TYPE_NOT_ITERABLE = (
    None,
    1,
    -1,
    2.5,
    True,
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
    *BAD_NEIGHBOURHOOD_NOT_UNIQUE_ITERABLE_SRC,
    *(list(case) for case in BAD_NEIGHBOURHOOD_NOT_UNIQUE_ITERABLE_SRC),
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
