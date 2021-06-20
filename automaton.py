#!/usr/bin/env python
# coding=utf-8

'''
manipulation of cellular automata generations
'''

# pipenv shell

# standard library imports
# import os
# import sys
# from typing import Tuple
from collections.abc import Iterable
from typing import Union

# related third party imports
# import tkinter.filedialog as tkf

# local application/library specific imports
# from commontools import constant

# validate data against type hints
# https://stackoverflow.com/questions/50563546/validating-detailed-types-in-python-dataclasses

class AutomataUniverse:
    '''properties for a cellular automata universe

    :property dimensions: the number of dimension for the universe
    :type: int
    :property neighbourhood: all cell addresses that are neighbors of the universe origin
    :type: frozenset of coordinate tuples for the origin cell address
    :property survival: neighbour counts required for living cell to survive to the next generation
    :type: frozen set of integers
    :property birth: neighbour counts required for an empty cell to spawn into the next generation
    :type: frozen set of integers
    '''

    automaton_address_type = tuple[int, ...]
    neighbourhood_input_type = Iterable[automaton_address_type]
    propagation_rule_input_type = Iterable[int]
    automaton_cell_group_type = frozenset[automaton_address_type]
    automaton_rule_count_type = frozenset[int]

    def __init__(self,
                neighbourhood: neighbourhood_input_type,
                survival_counts: propagation_rule_input_type,
                birth_counts: propagation_rule_input_type):
        self._origin_neighbourhood = frozenset(neighbourhood)
        self._survive = frozenset(survival_counts)
        self._birth = frozenset(birth_counts)
        self._dimensions = None
        self._validate_universe_data(len(neighbourhood))
        self._check_propagation_type(self._survive, len(survival_counts))
        self._check_propagation_type(self._birth, len(birth_counts))
        if 0 in self._birth:
            # Given the way the processing is done, "zero" neighbours is not a valid 'birth'
            # propagation value. It would fill the whole universe that was not a neighbor of
            # the starting generation at the first iteration
            raise ValueError("zero is not a valid birth propagation rule value")

    # properties : getter, setter, deleter methods

    @property
    def dimensions(self) -> int:
        return self._dimensions

    @property
    def neighbourhood(self) -> automaton_cell_group_type:
        return self._origin_neighbourhood

    @property
    def neighbourhood_size(self) -> int:
        '''get the size (number of cells) in the neighbourhood

        :returns neighbourhood_size: number of possible neighbours for a single cell
        :rtype: int
        '''
        return len(self._origin_neighbourhood)

    @property
    def survival_rules(self) -> automaton_rule_count_type:
        return self._survive

    @property
    def birth_rules(self) -> automaton_rule_count_type:
        return self._birth

    def __hash__(self):
        # hash of the automata universe configuration
        return hash((self.survival_rules, self.birth_rules, self.neighbourhood))

    # end of property methods

    def is_universe_address(self, address: automaton_address_type) -> bool:
        '''check whether argument is a valid universe address

        :param address: coordinates for cell in automata universe
        :type address: tuple of integers of size self.dimensions
        :return: True when address is a valid universe cell location
        :rtype: bool
        '''
        if not isinstance(address, tuple):
            return False
        dimensions = len(address)
        if not dimensions == self.dimensions:
            return False
        for coord in address:
            if not isinstance(coord, int):
                return False
        return True

    def validate_address(self, address: automaton_address_type):
        '''verify address is valid for the universe

        :param address: coordinates for cell in automata universe
        :type address: tuple of integers of size self.dimensions
        :raises: TypeError, ValueError
        '''
        if not isinstance(address, tuple):
            raise TypeError((type(address), "automata universe address is not a tuple"))
        dimensions = len(address)
        if not dimensions == self.dimensions:
            raise ValueError((self.dimensions, dimensions, "automata universe address "
                "does not have the same number of dimensions as the universe"))
        for coord in address:
            if not isinstance(coord, int):
                raise TypeError((type(coord),
                    "automata universe address coordinate is not an integer"))
    # end def validate_address()

    def get_neighbours(self, address: automaton_address_type) -> automaton_cell_group_type:
        '''get the set of neighbours for a cell address

        :param address: coordinates for cell in n-space
        :type address: tuple of integers
        :returns neighbours:
        :rtype neighbours: frozenset of coordinate tuples with same dimensionality as address
        '''
        self.validate_address(address)
        base_address = list(address)
        neighbourhood = []
        for neighbour in self._origin_neighbourhood:
            relative_address = list(neighbour)
            neighbourhood.append(tuple(sum(dim) for dim in zip(base_address, relative_address)))
        return frozenset(neighbourhood)
    # end def get_neighbours()

    def _validate_universe_data(self, input_size: int):
        '''check that the cell neighbourhood is structurally correct

        Also extract the number of dimensions for the automata universe

        :param input_size: the length of the neighbourhood Iterable used to create the universe
        :type input_size: int
        :output self._dimensions: the number of dimensions for the automata universe
        :otype: int
        :raises: TypeError, ValueError
        '''
        self._check_neighbourhood_type(input_size)

        # the origin point can not a neighbour (of the origin)
        origin = tuple([0] * self.dimensions)
        if origin in self._origin_neighbourhood:
            raise ValueError((origin, "the universe origin is not a valid neighbourhood address"))

        # every address in the neighbourhood must have the origin as one of its neighbours
        for addr in self._origin_neighbourhood:
            neighbours = self.get_neighbours(addr)
            if not origin in neighbours:
                raise ValueError((addr, "no symmetric address in neighbourhood"))
    # end def _validate_universe_data()

    def _check_neighbourhood_type(self, input_size: int):
        '''check that the cell neighbourhood is structurally correct

        Also extract the number of dimensions for the automata universe

        :param input_size: the length of the neighbourhood Iterable used to create the universe
        :type input_size: int
        :output self._dimensions: the number of dimensions for the automata universe
        :otype: int
        :raises: TypeError, ValueError
        '''
        # match self._origin_neighbourhood to the constructor neighourhood parameter typehint
        neighbourhood_size = len(self._origin_neighbourhood)
        if not neighbourhood_size == input_size:
            raise ValueError((input_size, neighbourhood_size,
                "Neighbourhood addresses are not unique"))
        if neighbourhood_size < 2:
            raise ValueError((neighbourhood_size,
                "neighbourhood does not contain at least 2 neighbours"))
            # Two neighbours is the absolute minimum to maintain symmetry. It would be valid
            # for a simple one dimensional universe
        # use first element to get the number of dimension to match for the rest
        neighbourhood_iter = iter(self._origin_neighbourhood)
        address = next(neighbourhood_iter)
        if not isinstance(address, tuple):
            raise TypeError((type(address), "element of neighbourhood is not a tuple"))
        self._dimensions = len(address)
        if self.dimensions < 1: # support for 1 dimensional automata?? or only 2+?
            raise TypeError((self.dimensions,
                "neighbourhood address does not contain at least 1 coordinate"))
        while True:
            self.validate_address(address)
            try:
                address = next(neighbourhood_iter)
            except StopIteration:
                break
        # end while True:
        # self._origin_neighbourhood confirmed to be a frozenset with at least 2 members, with all
        # members being tuples with only integer elements. All of the tuples have at least
        # one element, and all tuples have the same number of elements.
    # end def _check_neighbourhood_type(input_size):

    def _check_propagation_type(self, counts_rule: automaton_rule_count_type, input_size: int):
        '''check that a cell propagation rule is structurally correct with valid values

        :param counts_rule:
        :type counts_rule: frozen set of integers
        :param input_size: the length of the Iterable used to create the propagation rule
        :type input_size: int
        :raises: TypeError, ValueError
        '''
        rule_size = len(counts_rule)
        if not rule_size == input_size:
            raise ValueError((input_size, rule_size,
                "cell propagation rule counts are not unique"))
        # it would be "odd", but technically valid to have no counts for the rule
        # if rule_size < 1:
        for count in counts_rule:
            if not isinstance(count, int):
                raise TypeError((type(count), "cell propagation count case is not an integer"))
            if not 0 <= count <= self.neighbourhood_size:
                raise ValueError((count, self.neighbourhood_size,
                    "cell propagation count case is not between 0 and neighbourhood size"))
    # end def _check_propagation_type()
# end class AutomataUniverse


class Automaton:
    '''Storage for and operations on a cellular automata

    :property dimensions: the number of dimension for the current automaton
    :type: int
    :property iteration: the current generation sequence number (starts at 0)
    :type: int
    :property neighbourhood: all cells that are neighbors or the universe origin
    :type: frozenset of coordinate tuples for the origin cell address
    '''

    # aliases for typehint patterns
    automaton_address_type = tuple[int, ...]
    automaton_bounding_box_type = tuple[tuple[int, ...], tuple[int, ...]]
    automaton_neighbourhood_type = frozenset[automaton_address_type]
    automaton_cell_group_type = set[automaton_address_type]
    automaton_rules_type = tuple[frozenset[int], frozenset[int]]
    automaton_cells_type = Union[Iterable[automaton_address_type], automaton_address_type]

    def __init__(self, neighbourhood: automaton_neighbourhood_type,
            rules: tuple[frozenset[int], frozenset[int]]):
        '''constructor

        :param neighourhood: the neighbors of the origin point (0, 0, …)
        :type neighourhood: frozen set of tuples of integers; all tuples with same number elements
        :param rules:
        :type rules: tuple with 2 elements, each of which is a frozenset of integers
        '''
        self._universe = AutomataUniverse(neighbourhood, rules[0], rules[1])
        self._generation = set()
        self._iteration = 0

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
    def iteration(self, new_iteration):
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
    def neighbourhood(self) -> automaton_neighbourhood_type:
        return self._universe.neighbourhood

    @property
    def neighbourhood_size(self) -> int:
        '''get size of neighbourhood

        :returns neighbourhood_size: number of neighbours for a single cell
        :rtype: int
        '''
        return len(self._universe.neighbourhood_size)

    @property
    def generation(self) -> automaton_neighbourhood_type:
        '''get the universe for the current generation

        :returns current: living cells in the current generation
        :rtype: frozenset of cell address tuples
        '''
        return frozenset(self._generation)
    # end generation property getter

    @property
    def generation_extent(self) -> automaton_bounding_box_type:
        return self._get_extent(self.generation)

    # end of property methods

    # general methods

    def __hash__(self):
        # hash of the configuration and dynamic data of the automaton
        return hash((self._universe.survival_rules, self._universe.birth_rules,
            self._universe.neighbourhood, self.iteration, self.generation))

    def _get_extent(self, cells: automaton_cell_group_type) \
            -> automaton_bounding_box_type:
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

    def _offset_cell_block(self, cells: automaton_cell_group_type,
            offset_vector: automaton_address_type) -> automaton_cell_group_type:
        '''add offset «vector» to each cell coordinate'''
        assert len(offset_vector) == self.dimensions
        # all offset integers
        return [tuple(base + offset for base, offset in zip(c, offset_vector)) for c in cells]

    def _normalize_cells(self, cells: automaton_cell_group_type) -> automaton_bounding_box_type:
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
        normalized = self._offset_cell_block(cells, offset_vector)
        n_max = self._offset_cell_block(bounding_box, offset_vector)
        # print("bb", bounding_box, "offset bb", n_max) # DEBUG
        cells.clear()
        cells.update(normalized)
        # return (tuple(offset_vector), n_max[1])
        return bounding_box
    # end def _normalize_cells()

    def merge_cells(self, cells: automaton_cells_type):
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

    # def erase_cells(self, cells: automaton_cells_type):
        # '''delete living cells from the current generation

        # Cells that do not exist in the current generation are ignored

        # :param cells:
        # :type cells: single cell address tuple or iterable of cell address tuples
        # '''
    # def _get_expanded_neighborhood(self) -> automaton_neighbourhood_type:
        # '''neighbourhood that covers as far as it is possible of a cell to interact

        # Is the union of the standard neighbourhood of every neighbourhood address an accurate
        # representation of the extended neighbourhood? It should be for a `normal` cellular
        # automata neighbourhood. Will it be for all possible special cases? With the
        # neighbour symmetry rule enforced in this code, it SHOULD work.
    #     :returns increased_neighbourhood: address of all cells that could interact with the
        #         origin cell in the next generation
        # :rtype: frozenset of cell address tuples
        # '''
    # def get_connected_cells(self, address: automaton_address_type) -> automaton_neighbourhood_type:
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
    # def flip_rotate(self)
        # rules and code to flip and rotate «sets of» cells in ways that to not change the
        # properties of the set
        # -- in general «not always» any of these operation on the cell neighourhood will return
        #   the cell neighbourhood unchanged.

    def step(self):
        '''iterate from the current generation to the next

        Any cell in the current generation with a number of neighbors that is in the _survial rule
        continues to the next generation

        Any neighbor location of a cell in the current generation that does hold a cell is a
        candidate to spawn a cell in«to» the next generation

        Any candidate location with a number of neighbor cells in the current generation that is
        in the _birth rule becomes a new living cell in the next generation
        '''
        next_generation = set() # start with an empty next generation cell set
        womb_candidates = set() # no initial candidates for new cells either
        for living_cell in self.generation:
            cell_neighbourhood = self._universe.get_neighbours(living_cell)
            # print(living_cell, "neighbourhood", cell_neighbourhood) # DEBUG
            cell_neighbors = cell_neighbourhood.intersection(self._generation)
            neighbour_count = len(cell_neighbors)
            # print(living_cell, neighbour_count, "neighbors", cell_neighbors) # DEBUG
            print(living_cell, neighbour_count, "neighbors") # DEBUG
            if neighbour_count in self._universe.neighbourhood:
                next_generation.add(living_cell)
            empty_cell = cell_neighbourhood.difference(cell_neighbors)
            # print(living_cell, len(empty_cell), "no neighbours", empty_cell) # DEBUG
            assert neighbour_count + len(empty_cell) == self.neighbourhood_size, \
                "living and dead neighbours should add up to neighbourhood size"
            # print(living_cell, neighbour_count, "neighbors,",
            #     len(empty_cell), "no neighbours") # DEBUG
            womb_candidates.update(empty_cell)
        # print(len(womb_candidates), "womb candidates", womb_candidates) # DEBUG
        print(len(womb_candidates), "womb candidates") # DEBUG
        for womb_cell in womb_candidates:
            womb_neighbourhood = self._universe.get_neighbours(womb_cell)
            # print(womb_cell, "womb neighbourhood", womb_neighbourhood) # DEBUG
            womb_neighbors = womb_neighbourhood.intersection(self._generation)
            parent_count = len(womb_neighbors)
            # print(womb_cell, parent_count, "parents", womb_neighbors) # DEBUG
            if parent_count in self._universe.birth_rules:
                next_generation.add(womb_cell)
        self._generation = next_generation
        self._iteration += 1
# end class Automaton

def my_main():
    '''wrapper for test/start code so that variables do not look like constants'''
    # code to minimally exercise the class methods
    neighbours = [
        (-1,-1), (-1,0), (-1,1),
        (0,1),           (0,-1),
        (1,-1),  (1,0),  (1,1)
    ]
    instance = Automaton(neighbours, [[2,3], [3]])
    print("dimensions", instance.dimensions)
    print("neighbourhood", instance.neighbourhood)
    print("generation", instance.iteration)
    instance.iteration = 2
    print("generation", instance.iteration)
    instance.merge_cells((0, 0))
    # print(instance.current)
    # for cell in instance.current:
    #     print(cell)
    print("current generation", list(instance.generation))


# Standalone module execution
if __name__ == "__main__":
    my_main()
