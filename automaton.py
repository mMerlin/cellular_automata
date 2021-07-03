#!/usr/bin/env python
# coding=utf-8

'''
manipulation of cellular automata generations
'''

# pipenv shell

# standard library imports
# import os
# import sys
from collections.abc import Iterable
from typing import Hashable
# from threading import Lock

# related third party imports

# local application/library specific imports
import automata_typehints as AHint
from automata_universe import AutomataUniverse
from automata_transforms import AutomataTransforms

# class AutomataCells:
#     '''Storage and operations for a group of cells in an AutomataUniverse
#     '''

class Automaton:
    '''Storage for and operations on a cellular automata

    :property universe: the automata universe configuration
    :type universe: AutomataUniverse
    :property dimensions: the number of dimension for the universe
    :type: int
    :property neighbourhood: all cell addresses that are neighbors of the universe origin
    :type: frozenset of coordinate tuples for the origin cell address
    :property generation: living cells
    :type generation: set of automata universe cell address tuples
    :property generation_extent:
    :type generation_extent: tuple of 2 universe cell address tuples
    :property iteration: the current generation sequence number (starts at 0)
    :type iteration: int
    :property population: the number of living cells in the current generation
    :type population: int
    '''

    def __init__(self, universe: AutomataUniverse) -> None:
        '''constructor

        :param universe: parent cellular automata universe configuration
        :type universe: AutomataUniverse
        '''
        self._universe = universe
        self._generation = set()
        self._iteration = 0
        self._transforms = AutomataTransforms(universe)

    # properties : getter, setter, deleter methods

    @property
    def dimensions(self) -> int:
        return self._universe.dimensions

    @property
    def iteration(self) -> int:
        return self._iteration

    @property
    def population(self) -> int:
        '''the number of living cells in the current universe generation'''
        return len(self._generation)

    @iteration.setter
    def iteration(self, new_iteration) -> None:
        '''set generation sequence number for the current set

        :param new_iteration: generation sequence number to use
        :rtype new_iteration: int
        '''
        if not (isinstance(new_iteration, int) and new_iteration >= 0):
            raise TypeError(new_iteration,
                "generation sequence number must be an integer equal to or greater than zero")
        self._iteration = new_iteration
    # end iteration() property setter

    @property
    def neighbourhood(self) -> AHint.NeighbourhoodType:
        return self._universe.neighbourhood

    @property
    def neighbourhood_population(self) -> int:
        '''get size of neighbourhood

        :returns population: number of neighbours for a single cell
        :rtype: int
        '''
        return len(self._universe.neighbourhood_population)

    @property
    def generation(self) -> AHint.CellGroupSnapshotType:
        '''get the universe content for the current generation

        :returns current: living cells in the current generation
        :rtype: frozenset of cell address tuples
        '''
        return frozenset(self._generation)
    # end generation property getter

    @property
    def generation_extent(self) -> AHint.BoundingBoxType:
        return self._get_extent(self.generation)

    # end of property methods

    # general methods

    def clear(self) -> None:
        self._generation.clear()

    def __hash__(self) -> int:
        # hash of the configuration and dynamic data of the automaton
        return hash((self._universe.survival_rules, self._universe.birth_rules,
            self._universe.neighbourhood, self._transforms, self.iteration, self.generation))

    def _get_extent(self, cells: AHint.CellGroupWorkingType) -> AHint.BoundingBoxType:
        '''determine the n-dimensional bounding box for a set of cells

        :param cells:
        :type cells: «forzen»set of cell coordinate tuples
        :returns: minimum and maximum corner coordinates of bounding box
        :rtype: tuple[tuple[int, ...],tuple[int, ...]]
        '''
        box_min = [float('inf')] * self.dimensions
        box_max = [-box_min[0]] * self.dimensions
        for coord in cells:
            box_min = [min(new, cur) for new, cur in zip(coord, box_min)]
            box_max = [max(new, cur) for new, cur in zip(coord, box_max)]
        return tuple((tuple(box_min), tuple(box_max)))
    # end def _get_extent()

    def _normalize_cells(self, cells: AHint.CellGroupWorkingType) \
            -> AHint.BoundingBoxType:
        # '''shift group of cells to fit against all axis in the first quadrant

        # The lowest coordinate value in every dimension will be zero.

        # The bounding box for the normalized set is the origin (0, 0, …) and the returned
        # extent address.

        # :param cells: (in|out)
        # :type cells: set of cell address tuples
        # :return extent: bounding box of the original block of cells
        # :rtype: tuple containing minimum and maximum coordinate tuples
        # NOTE  return_value[0] is the offset vector need to move the normalized
        #       block back to its original location
        # '''
        bounding_box = self._get_extent(cells)
        offset_vector = [-1 * delta for delta in bounding_box[0]]
        normalized = self._universe.cell_group_translate(cells, offset_vector)
        # n_max = self._universe.cell_group_translate(bounding_box, offset_vector)
        # print("bb", bounding_box, "offset bb", n_max) # DEBUG
        cells.clear()
        cells.update(normalized)
        # return (tuple(offset_vector), n_max[1])
        return bounding_box
    # end def _normalize_cells()

    def merge_cells(self, cells: AHint.CellorCellsType) -> None:
        '''add living cells to the current generation

        Cells that already exist in the current generation are ignored

        :param cells:
        :type cells: single cell address tuple or iterable of cell address tuples
        '''
        # validate that input cells are «all» address tuples
        # a tuple is iterable, so need to be careful with the single case test
        if self._universe.is_universe_address(cells):
            self._generation.add(cells)
            return
        if not isinstance(cells, Iterable):
            raise TypeError((type(cells), "cells object must be iterable"))
        for addr in cells:
            self._universe.validate_address(addr)
        self._generation.update(cells)
    # end merge_cells()

    def step(self) -> None:
        '''iterate from the current generation to the next'''
        next_generation = self._universe.step(self.generation)
        # self.generation = next_generation
        self._generation.clear()
        self._generation.update(next_generation)
    # end def step()

    def add_transform(self, key: Hashable, transform: AHint.TransformInputType) -> None:
        self._transforms.add_transform_cycle(key, transform)
    # end def add_transform()

    # def erase_cells(self, cells: AutomataTypeHints.cell_or_cells_type):
        # '''delete living cells from the current generation

        # Cells that do not exist in the current generation are ignored

        # :param cells:
        # :type cells: single cell address tuple or iterable of cell address tuples
        # '''
    # def _get_expanded_neighborhood(self) -> AutomataTypeHints.neighbourhood_type:
        # '''neighbourhood that covers as far as it is possible of a cell to interact

        # Is the union of the standard neighbourhood of every neighbourhood address an accurate
        # representation of the extended neighbourhood? It should be for a `normal` cellular
        # automata neighbourhood. Will it be for all possible special cases? With the
        # neighbour symmetry rule enforced in this code, it SHOULD work.
    #     :returns increased_neighbourhood: address of all cells that could interact with the
        #         origin cell in the next generation
        # :rtype: frozenset of cell address tuples
        # '''
    # def get_connected_cells(self, address: AutomataTypeHints.address_type) -> AutomataTypeHints.neighbourhood_type:
        # '''set of cells that interact with the start cell and each other

        # Collect all cells that are in the extended neighbourhood of the starting cell, and the
        # extended neighborhoods of those neighbours recursively.

        # The extended neighbourhood needs to be used, because empty standard neighbourhood cells
        # at the edge of the group can be affected by existing cell one step (unit) further away.

        # :param address:
        # :type: tuple of integer dimension coordinates
        # :returns connected_set:
        # :rtype: frozenset of cell address tuples
        # '''
# end class Automaton


SQUARE_GRID_NEIGHBORS = [
    (-1,-1), (-1,0), (-1,1),
    (0,1),           (0,-1),
    (1,-1),  (1,0),  (1,1)
]
ROTATE90 = [(0, -1), (1, 0)]

def _validate_and_expand_matrices(atm: Automaton, uni: AutomataUniverse,
        r_input: AHint.RotateReflectInputType) -> None:
    '''check supplied matrices and expand to all permutations

    Verify that the supplied rotation and reflection matrices are consistent with the universe.
    Each matrix must be cyclic. That is, repeated operations using the matrix must return
    the data to the starting pattern in a few steps.

    The raw input needs to be an iterable object containing one or more square matrices.
    Each matrix must have a number or elements matching the dimensions of the universe, and
    the elements are cell addresses.

    After validation, build nested tuples from the input, to be compatible with hashing.

    :param r_input: one or more rotation and reflection matrices
    :rtype r_input: iterable of iterable of cell address tuples
    :output self._rotate_reflect: all matrices that generate equivalent rotated and
        reflected (flipped) cell patterns
    :otype self._rotate_reflect: tuple of tuples of cell address tuples
    '''
    print("start _validate_and_expand_matrices") # DEBUG
    # RotateReflectInputType = Iterable[Iterable[CellAddressType]]
    r_matrices = set()
    r_mat_count = len(r_matrices)
    if not isinstance(r_input, Iterable):
        # TypeError: 'NoneType' object is not iterable
        # raise TypeError("{} object is not iterable".format(type(r_input)))
        raise TypeError(type(r_input), "object is not iterable")
    if len(r_input) < 1:
        raise ValueError("universe must have at least one rotation or reflection symmetry")
    for mat in r_input:
        if not isinstance(mat, Iterable):
            raise TypeError(type(mat), "object is not iterable")
        if not len(mat) == atm.dimensions:
            raise TypeError("{} rows is not valid for matrix operation in universe "
                "with {} dimensions".format(len(mat), atm.dimensions))
        working_matrix = []
        for adr in mat:
            uni.validate_address(adr)
            working_matrix.append(adr)
        # r_matrices.update(tuple(adr) for adr in mat)
        r_matrices.add(tuple(working_matrix))
        if len(r_matrices) != 1 + r_mat_count:
            raise ValueError(tuple(working_matrix), "duplicate symmetry matrix encountered")
        r_mat_count += 1
    print("working r_and_m", r_matrices, hash(frozenset(r_matrices))) # DEBUG
    for ele in r_matrices:
        print(hash(ele)) # DEBUG
    # all of the rotation or reflection matrices are at least structurally valid for the
    # configured universe. They have the correct shape.
    # check for duplicates and cycle sizes
    # first 1000 prime numbers
    #   https://en.wikipedia.org/wiki/List_of_prime_numbers#The_first_1000_prime_numbers
    # not all counting numbers might be valid as cell address coordinates. Currently there is
    # no process to validate that case.
    mat = list(r_matrices)[0]
    test = (5, 37) # DEBUG
    print(type(test), test, hash(test)) # DEBUG
    test = uni.vector_dot_product(test, mat) # DEBUG
    print(type(test), test, hash(test)) # DEBUG
    test = uni.vector_dot_product(test, mat) # DEBUG
    print(type(test), test, hash(test)) # DEBUG
    test = uni.vector_dot_product(test, mat) # DEBUG
    print(type(test), test, hash(test)) # DEBUG
    test = uni.vector_dot_product(test, mat) # DEBUG
    print(type(test), test, hash(test)) # DEBUG
    print("") # DEBUG

    test = None # _test_prime_set(uni.dimensions)
    print("prime cells", test, hash(test)) # DEBUG
    test = uni.cell_group_transform(test, mat) # DEBUG
    print(test, hash(test)) # DEBUG
    test = uni.cell_group_transform(test, mat) # DEBUG
    print(test, hash(test)) # DEBUG
    test = uni.cell_group_transform(test, mat) # DEBUG
    print(test, hash(test)) # DEBUG
    test = uni.cell_group_transform(test, mat) # DEBUG
    print(test, hash(test)) # DEBUG
    print("") # DEBUG

    # self._rotate_reflect = frozenset(r_matrices)
# end def _validate_and_expand_matrices()

def matrix_test(instance: Automaton) -> None:
    '''exercise transform validation and setup'''
    tr_mat = ROTATE90 # GOOD
    tr_key = "Rotate 90°"
    instance.add_transform(tr_key, tr_mat)

def my_main() -> None:
    '''wrapper for test/start code so that variables do not look like constants'''
    # code to minimally exercise the class methods
    # universe = AutomataUniverse(SQUARE_GRID_NEIGHBORS, [2,3], [3], [rotate90, horizontal_mirror])
    universe = AutomataUniverse(SQUARE_GRID_NEIGHBORS, [2,3], [3])
    instance = Automaton(universe)
    print("dimensions", instance.dimensions)
    print("neighbourhood", instance.neighbourhood)
    print("generation", instance.iteration)
    instance.iteration = 2
    print("generation", instance.iteration)
    instance.merge_cells((0, 0))
    # print(instance.generation)
    # for cell in instance.generation:
    #     print(cell)
    print("current generation", list(instance.generation))
    print()
    matrix_test(instance)


# Standalone module execution
if __name__ == "__main__":
    my_main()
