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
from typing import Iterable, Hashable, Union


# def local_function_definition(«arguments»…):
# def local_function_definition(«[argument«: «type»»[, …]»):
def identity_matrix(dimensions: int) -> tuple[tuple[int, ...]]:
    '''identity matrix for the specified number of dimensions

    :param dimensions:
    :type dimensions: int
    :returns: identity matrix for the specified number of dimensions
    :rtype: tuple containing dimension tuples of dimension integers
    :raises: TypeError, ValueError
    '''
    if not isinstance(dimensions, int):
        raise TypeError("Identity Matrix dimensions must be an intger, not {}".format(
            type(dimensions)))
    if dimensions < 1:
        raise ValueError("{} is not a invalid; Identity Matrix dimension "
            "must be greater than zero".format(dimensions))
    # matrix of all zeros that does NOT have linked rows
    identity = [list(row) for row in ((0,) * dimensions,) * dimensions]
    for i in range(dimensions):
        identity[i][i] = 1
    return tuple(tuple(row) for row in identity)

def vector_dot_product(vector: tuple[Union[int, float], ...],
        matrix: Iterable[Iterable[Union[int, float]]]) -> tuple[Union[int, float], ...]:
    '''transform a vector through a transformation matrix

    :param vector: cell address
    :type vector: tuple of integer coordinates for universe dimensions
    :param matrix: rotation matrix
    :type matrix: tuple of «dimension» cell address coordinates
    :raises: TypeError
    '''
    # This code does not work with one dimensional data
    return tuple(sum(x * y for x, y in zip(row, vector)) for row in matrix)
    # return tuple([sum(x * y for x, y in zip(row, vector)) for row in matrix])
    # product = []
    # for ele in matrix:
    #     product.append(x*y for x,y in zip(ele, vector))
    # result = [sum(x * y for x, y in zip(row, vector)) for row in matrix]
    # return tuple(result)

def matrix_transpose(matrix: Iterable[Iterable[Hashable]]) \
        -> tuple[tuple[Hashable, ...]]:
    '''transpose matrix (rotate 90°)

    :param matrix: matrix to transpose
    :type matrix: iterable of «count1» iterables of «count2» Immutable elements
    :returns: transposed (reversed x and y) matrix
    :rtype: tuple of «count2» tuples of «count1» Immutable elements
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
        while(index < dim and mat[index][i] == 0):
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

def matrix_determinant2(matrix: Iterable[Iterable[int]]) -> Union[int, float]:
    '''calculated the determinant of a square matrix

    :param matrix: the matrix
    :type matrix: any iterable containing iterables of numbers.
                The inner iterables are all the same size
    :returns: determinant of the matrix
    :rtype: number
    '''
    mat = [list(row) for row in matrix] # make a fully indexable copy
    return _recursive_matrix_determinant(mat, 1)

def _recursive_matrix_determinant(matrix: list[list[Union[int, float]]],
        mul: Union[int, float] = 1) -> Union[int, float]:
    '''recursive calculation of matrix determinant using minors and cofactors

    :param matrix: the matrix
    :type matrix: any iterable containing iterables of numbers.
                The inner iterables are all the same size
    :param mul: accumulated multiplier
    :type mul: number
    :returns: determinant of the matrix
    :rtype: number
    '''
    width = len(matrix)
    if width == 1:
        return mul * matrix[0][0]

    sign = -1
    total = 0
    for i in range(width):
        mat = []
        for j in range(1, width):
            buff = []
            for k in range(width):
                if k != i:
                    buff.append(matrix[j][k])
            mat.append(buff)
        sign *= -1
        total += mul * _recursive_matrix_determinant(mat, sign * matrix[0][i])
    return total


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
