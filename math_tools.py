#!/usr/bin/env python
# coding=utf-8
# pylint: disable=C0415

'''
some purely mathematical operations
'''

# standard library imports
# start needed to use `kill -SIGUSR1 «pid»` to attach pdb to running session
import os  # DEBUG
import signal
# end needed to use `kill -SIGUSR1 «pid»` to attach pdb to running session
from typing import Union
from collections.abc import Iterable
# import os
# import sys


# def local_function_definition(«arguments»…):
# def local_function_definition(«[argument«: «type»»[, …]»):
def identity_matrix(dimensions: int) -> tuple[tuple[int, ...]]:
    '''identity matrix for the specified number of dimensions

    :param dimensions:
    :type dimensions: int
    :returns: identity matrix for the specified number of dimensions
    :rtype: tuple containing dimension tuples of dimension integers
    '''
    # matrix of all zeros that does NOT have linked rows
    assert isinstance(dimensions, int)
    assert dimensions > 0
    identity = [list(row) for row in ((0,) * dimensions,) * dimensions]
    for i in range(dimensions):
        identity[i][i] = 1
    return tuple(tuple(row) for row in identity)

def vector_dot_product(vector: tuple[int, ...], matrix: Iterable[Iterable[int]]) \
        -> tuple[int, ...]:
    '''transform a cell address vector through a rotation or reflection matrix

    matrix is not validated here. Expected to be validated once by caller, then
    used multiple times with different vectors

    :param vector: cell address
    :type vector: tuple of integer coordinates for universe dimensions
    :param matrix: rotation matrix
    :type matrix: tuple of «dimension» cell address coordinates
    '''
    # Does this still work with a one dimensional vector and matrix?
    return tuple(sum(x * y for x, y in zip(row, vector)) for row in matrix)
    # return tuple([sum(x * y for x, y in zip(row, vector)) for row in matrix])
    # product = []
    # for ele in matrix:
    #     product.append(x*y for x,y in zip(ele, vector))
    # result = [sum(x * y for x, y in zip(row, vector)) for row in matrix]
    # return tuple(result)

def matrix_transpose(matrix: Iterable[Iterable[int]]) -> tuple[tuple[int, ...]]:
    '''transpose matrix (rotate 90°)

    :param matrix: matrix to transpose
    :type matrix: iterable of «count1» tuples of «count2» elements
    :returns: transposed (reversed x and y) matrix
    :rtype: tuple of «count2» tuples of «count1» elements
    '''
    return tuple(tuple(matrix[j][i] for j in range(len(matrix))) for i in range(len(matrix[0])))

def matrix_transform(matrix: Iterable[Iterable[int]], transform: Iterable[Iterable[int]]) \
        -> tuple[tuple[int, ...]]:
    '''transform a matrix through a rotation or reflection matrix

    C[i,j] = A[i,*] · B[*,j]

    :param matrix: «square?» matrix
    :type matrix: tuple of «dimension» cell address coordinates
    :param transform: square matrix
    :type transform: tuple of «dimension» cell address coordinates
    :returns: matrix dot product
    :rtype: tuple of «dimension» cell address tuples
    '''
    return tuple(vector_dot_product(row, matrix_transpose(transform)) for row in matrix)

def matrix_determinant(matrix: Iterable[Iterable[int]]) -> Union[int, float]:
    '''calculated the determinant of a square matrix

    The reference code used does not look at all pythonic. Direct translation from C.

    :param matrix: the matrix
    :type matrix: any iterable containing iterables of numbers.
                The inner iterables are all the same size
    :returns: determinant of the matrix
    :rtype: number
    '''
    dim = len(matrix)
    mat = [list(row) for row in matrix] # make a fully indexable copy
    temp = [0] * dim # temporary array for storing row
    total = 1
    det = 1 # initialize result

    for i in range(dim): # traverse diagonal elements
        index = i
        while(mat[index][i] == 0 and index < dim):
            index += 1
        if index == dim:
            continue # the determinant is zero
        if index != i: # swap diagonal element and index rows
            for j in range(dim):
                mat[index][j], mat[i][j] = mat[i][j], mat[index][j]
            # determinant sign changes when shifting rows
            det = det * int(pow(-1, index - i))
        for j in range(dim): # diagonal row elements
            temp[j] = mat[i][j]
        for j in range(i + 1, dim): # traverse rows below diagonal element
            num1 = temp[i]   # diagonal element
            num2 = mat[j][i] # next row element
            for k in range(dim): # traverse every column of row
                # multiply to make diagonal element and next row element equal
                mat[j][k] = (num1 * mat[j][k]) - (num2 * temp[k])
            total = total * num1 # Det(kA) = kDet(A)
    for i in range(dim): # product of diagonal elements is determinant
        det = det * mat[i][i]
    return int(det/total) # Det(kA)/k = Det(A)

def my_main(): # pragma: no cover
    '''wrapper for test/start code so that variables do not look like constants'''
    dummy_start = identity_matrix(4)


def start_pdb(_sig, frame):  # DEBUG # pragma: no cover
    '''handler to start pdb session on interrupt'''
    import pdb
    pdb.Pdb().set_trace(frame)


# Standalone module execution
if __name__ == "__main__": # pragma: no cover
    signal.signal(signal.SIGUSR1, start_pdb)  # DEBUG
    print(  # DEBUG
        'Session process ID: {pid}\n\nuse: `kill -SIGUSR1 {pid}` to start pdb\n'  # DEBUG
        'NOTE: .pdbrc gets executed when signal received'.format(pid=os.getpid()))  # DEBUG
    my_main()
