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

class Automaton:
    '''Storage for and operations on a cellular automata

    :property dimensions: the number of dimension for the current automaton
    :type: int
    :property iteration: the current generation sequence number (starts at 0)
    :type: int
    :property neighbourhood: all cells that are neighbors or the universe origin
    :type: frozenset of coordinate tuples for the origin cell address
    '''

    # aliases for hinttype patterns
    automaton_address_type = tuple[int, ...]
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
        self._cell_neighbourhood = neighbourhood
        # self._dimensions = None
        self._validate_neighbourhood_data()
        self._validate_rule_data(rules)
        self._survive = rules[0]
        self._birth = rules[1]
        self._generation = set()
        self._iteration = 0

    # properties : getter, setter, deleter methods

    @property
    def dimensions(self) -> int:
        return self._dimensions

    @property
    def iteration(self) -> int:
        return self._iteration

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
        return self._cell_neighbourhood

    @property
    def neighbourhood_size(self) -> int:
        '''get size of neighbourhood

        :returns neighbourhood_size: number of neighbours for a single cell
        :rtype: int
        '''
        return len(self._cell_neighbourhood)
    # end neighbourhood_size() property getter

    @property
    def generation(self) -> automaton_neighbourhood_type:
        '''get the universe for the current generation

        :returns current: living cells in the current generation
        :rtype: frozenset of cell address tuples
        '''
        return frozenset(self._generation)
    # end current() property getter

    # end of property methodsNone

    # general methods

    def merge_cells(self, cells: automaton_cells_type):
        '''add living cells to the current generation

        Cells that already exist in the current generation are ignored

        :param cells:
        :type cells: single cell address tuple or iterable of cell address tuples
        '''
        # validate that input cells are «all» address tuples
        # a tuple is iterable, so need to be careful with the single case test
        addr_state = self.address_state(cells)
        if addr_state == 0: # single cell address
            self._generation.add(cells)
            return
        if not isinstance(cells, Iterable):
            raise TypeError((type(cells), "cells object must be iterable"))
        for addr in cells:
            addr_state = self.address_state(addr)
            if not addr_state == 0:
                self.address_state_fail(addr_state)
        self._generation.update(cells)
    # end merge_cells()

    # def erase_cells(self, cells: automaton_cells_type):
        # '''delete living cells from the current generation

        # Cells that do not exist in the current generation are ignored

        # :param cells:
        # :type cells: single cell address tuple or iterable of cell address tuples
        # '''
    # def _get_expanded_neighborhood(self) -> automaton_neighbourhood_type:
    #     '''neighbourhood that covers as far as it is possible of a cell to interact

    #     Is the union of the standard neighbourhood of every neighbourhood address an accurate
    #     representation of the extended neighbourhood? It should be for a `normal` cellular
    #     automata neighbourhood. Will it be for all possible special cases? With the
    #     neighbour symmetry rule enforced in this code, it SHOULD work.

    #     :returns increased_neighbourhood: address of all cells that could interact with the
    #             origin cell in the next generation
    #     :rtype: frozenset of cell address tuples
    #     '''
    # def get_connected_cells(self, address: automaton_address_type) -> automaton_neighbourhood_type:
    #     '''set of cells that interact with the start cell and each other

    #     Collect all cells that are in the extended neighbourhood of the starting cell, and the
    #     extended neighborhoods of those neighbours recursively.

    #     The extended neighbourhood needs to be used, because empty standard neighbourhood cells
    #     at the edge of the group can be affected by existing cell one step (unit) further away.

    #     :param address:
    #     :type: tuple of integer dimension coordinates
    #     :returns connected_set:
    #     :rtype: frozenset of cell address tuples
    #     '''
    # def normalize(self, cells: automaton_cell_group_type) -> automaton_address_type:
    #     '''shift group of cells to fit against all axis in the first quadrant

    #     The lowest coordinate value in every dimension will be zero.

    #     The bounding box for the normalized set is the origin (0, 0, …) and the returned
    #     extent address.

    #     :param cells: (in|out)
    #     :type: set of tuples where every tuple contains integer values, and every tuple contains
    #            the number of elements that matches the dimensions for the automaton configuration
    #     :return extent: the largest coordinate in each dimension for the normalized set
    #     :rtype: tuple of integers
    #     IDEA  instead of new single reference point extent, return original lower and upper
    #           extent address, so that the offset can be collected
    #       Alternatively, return the delta was well as the new extent maximum
    #       -- effectively the original minimum plus normalized maximum
    #     '''
    # def ?
    #     rules and code to flip and rotate «sets of» cells in ways that to not change the
    #     properties of the set
    #     -- in general «not always» any of these operation on the cell neighourhood will return
    #       the cell neighbourhood unchanged.
    # def ?
    #     hash of «normalized» set

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
            cell_neighbourhood = self.get_neighbours(living_cell)
            print(living_cell, "neighbourhood", cell_neighbourhood) # DEBUG
            cell_neighbors = cell_neighbourhood.intersection(self._generation)
            neighbour_count = len(cell_neighbors)
            print(living_cell, neighbour_count, "neighbors", cell_neighbors) # DEBUG
            if neighbour_count in self._survive:
                next_generation.add(living_cell)
            empty_cell = cell_neighbourhood.difference(cell_neighbors)
            print(living_cell, len(empty_cell), "no neighbours", empty_cell) # DEBUG
            womb_candidates.update(empty_cell)
        print(len(womb_candidates), "womb candidates", womb_candidates)
        for womb_cell in womb_candidates:
            womb_neighbourhood = self.get_neighbours(womb_cell)
            # print(womb_cell, "womb neighbourhood", womb_neighbourhood) # DEBUG
            womb_neighbors = womb_neighbourhood.intersection(self._generation)
            parent_count = len(womb_neighbors)
            print(womb_cell, parent_count, "parents", womb_neighbors) # DEBUG
            if parent_count in self._birth:
                next_generation.add(womb_cell)
        self._generation = next_generation
        self._iteration += 1

    def _validate_rule_data(self, rules: automaton_rules_type):
        '''extract and validate rules for cell survival and creation in the next generation

        :param rules:
        :type rules: tuple with 2 elements, each of which is a frozenset of integers
        '''
        if not isinstance(rules, tuple):
            raise TypeError((type(rules), "automaton rules are not a tuple"))
        if not len(rules) == 2:
            raise TypeError((len(rules), "automaton rules do not contain 2 conditions"))
        # print("rule sets", rules) # DEBUG
        for rule_counts in rules:
            if not isinstance(rule_counts, frozenset):
                raise TypeError((type(rule_counts), "automaton rule is not a frozenset"))
            live_cases = len(rule_counts)
            if not 0 < live_cases < self.neighbourhood_size:
                print("rule counts", rule_counts) # DEBUG
                raise TypeError((live_cases, self.neighbourhood_size,
                    "living cell count not between 1 and neighbourhood size"))
            for count in rule_counts:
                if not isinstance(count, int):
                    raise TypeError((type(count), "rule case is not an integer"))
                if not 0 < count <= self.neighbourhood_size:
                    raise TypeError((count, self.neighbourhood_size,
                        "rule case is not between 1 and neighbourhood size"))
    # end _validate_rule_data()

    def _validate_neighbourhood_data(self):
        '''check that the cell neighbourhood is structurally correct

        Also extract the number of dimensions for the automaton
        :output self._dimensions: the number of dimensions for the current automaton
        :type: int
        :raises: TypeError
        '''
        # match self._cell_neighbourhood to the constructor neighourhood parameter typehint
        if not isinstance(self._cell_neighbourhood, frozenset):
            raise TypeError((type(frozenset),
                "neighbourhood does not contain at least 2 neighbours"))
        neighbourhood_iter = iter(self._cell_neighbourhood)
        address = next(neighbourhood_iter)
        if not isinstance(address, tuple):
            raise TypeError((type(address), "element of neighbourhood is not a tuple"))
        self._dimensions = len(address)
        if self.dimensions < 1: # support for 1 dimensional automata?? or only 2+?
            # print(address) # DEBUG
            raise TypeError((self.dimensions,
                "neighbourhood address does not contain at least 1 coordinate"))
        while True:
            addr_state = self.address_state(address)
            if isinstance(addr_state, tuple):
                if addr_state[0] == 1:
                    raise TypeError((addr_state[1], "element of neighbourhood is not a tuple"))
                if addr_state[0] == 2:
                    raise TypeError((self.dimensions, addr_state[1],
                        "neighbourhood address coordinate is not an integer"))
                raise ValueError((addr_state[1], "unknown address state error code"))
            if not addr_state == 0:
                raise ValueError(addr_state, "unknown address state encountered")

            try:
                address = next(neighbourhood_iter)
            except StopIteration:
                break
        # end while True:
        # self.cell_neighbour confirmed to be a frozenset with at least 2 members, with all
        # members being tuples with only integer elements. All of the tuples have at least
        # one element, and all tuples have the same number of elements.

        # the origin point can not a neighbour (of the origin)
        origin = tuple([0] * self.dimensions)
        # print("neighbourhood", self._cell_neighbourhood) # DEBUG
        # print("origin", origin) # DEBUG
        if origin in self._cell_neighbourhood:
            raise TypeError((origin, "the origin is not a valid neighbourhood address"))

        # every address in the neighbourhood must have the origin as one of its neighbours
        for addr in self._cell_neighbourhood:
            neighbours = self.get_neighbours(addr)
            # print("neighbours", addr, neighbours) # DEBUG
            if not origin in neighbours:
                raise TypeError((addr, "no symmetric address in neighbourhood"))
    # end _validate_neighbourhood_data()

    def address_state(self, address: automaton_address_type) -> Union[int, tuple[int, any]]:
        '''validate an address for use with the current automaton configuration

        :param address: coordinates for cell in n-space
        :type address: tuple of integers of size self.dimensions
        :returns state: 0 = valid address;
            (1, type) = not a tuple;
            (2, dim) = wrong dimensionality;
            (3, type) = not integer coordinate
        :rtype: int | tuple[int, int | type]
        '''
        if not isinstance(address, tuple):
            return (1, type(address))
        dimensions = len(address)
        if not dimensions == self.dimensions:
            return (2, dimensions)
        for coord in address:
            if not isinstance(coord, int):
                return (3, type(coord))
        return 0
    # end address_state()

    def get_neighbours(self, address: automaton_address_type) -> automaton_neighbourhood_type:
        '''get the set of neighbours for a cell address

        :param address: coordinates for cell in n-space
        :type address: tuple of integers
        :returns neighbours:
        :rtype neighbours: frozenset of coordinate tuples with same dimensionality as address
        '''
        addr_state = self.address_state(address)
        if not addr_state == 0:
            self.address_state_fail(addr_state)
        base_address = list(address)
        neighbourhood = []
        # print("base", base_address) # DEBUG
        for neighbour in self._cell_neighbourhood:
            relative_address = list(neighbour)
            # print("relative", neighbour, relative_address) # DEBUG
            # print("zip", list(zip(base_address, relative_address))) # DEBUG
            neighbourhood.append(tuple(sum(dim) for dim in zip(base_address, relative_address)))
        return frozenset(neighbourhood)
    # end def get_neighbours()

    def address_state_fail(self, detected_state: Union[int, tuple[int, any]]):
        '''throw exceptions based on the detected address state

        This is only «to be» called when it is already known that there is a problem with the state

        :param detected_state:
        :type detected_state: integer or a tuple 2 elements; an integer plus an integer or a type
        '''
        if not isinstance(detected_state, tuple):
            raise ValueError(detected_state, "unknown address state encountered")
        if detected_state[0] == 1:
            raise TypeError((detected_state[1], "address is not a tuple"))
        if detected_state[0] == 2:
            raise TypeError((self.dimensions, detected_state[1],
                "address does not have the same number of dimension coordinates as automaton"))
        if detected_state[0] == 3:
            raise TypeError((detected_state[1], "address coordinate is not an integer"))
        raise ValueError((detected_state[1], "unknown address state error code"))
    # end address_state_fail()
# end class Automaton

def my_main():
    '''wrapper for test/start code so that variables do not look like constants'''
    # code to minimally exercise the class methods
    neighbours = frozenset((
        (-1,-1), (-1,0), (-1,1),
        (0,1),           (0,-1),
        (1,-1),  (1,0),  (1,1)
    ))
    instance = Automaton(neighbours, (frozenset((2,3)), frozenset([3])))
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
