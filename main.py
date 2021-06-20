#!/usr/bin/env python
# coding=utf-8

'''
shell to work with automaton class
'''

# pipenv shell


# standard library imports
# import os
# import sys
# from typing import Tuple

# related third party imports
# import tkinter.filedialog as tkf

# local application/library specific imports
# from commontools import constant
from automaton import AutomataUniverse, Automaton

# validate data against type hints
# https://stackoverflow.com/questions/50563546/validating-detailed-types-in-python-dataclasses

# setup functions to add patterns of cell to and automaton universe generation

def report_generation(instance: Automaton):
    '''show information about an Automaton generation'''
    # print(amn.generation)
    # for cell in amn.generation:
    #     print(cell)
    print("generation", instance.iteration, "population", instance.population,
        "extent", instance.generation_extent,
        "hashes", hash(instance), hash(instance.generation),
        "contains", list(instance.generation))

def offset_cells(cells: list[tuple[int,int]], offset_x: int = 0, offset_y: int = 0) -> \
        list[tuple[int,int]]:
    '''add x and/or y offsets to cell coordinate tuples in a list'''
    return [tuple(x + y for x, y in zip(c, (offset_x, offset_y))) for c in cells]

def _one_cell(instance: Automaton, offset_x: int = 0, offset_y: int = 0):
    '''add single living cell to existing universe'''
    instance.merge_cells((0 + offset_x, 0 + offset_y))

def _two_cell_h_pair(instance: Automaton, offset_x: int = 0, offset_y: int = 0):
    '''add 2 horizontal adjacent cells to existing universe'''
    template = [(0, 0), (1, 0)]
    # instance.merge_cells((0 + offset_x, 0 + offset_y))
    # instance.merge_cells((1 + offset_x, 0 + offset_y))
    # instance.merge_cells([(0 + offset_x, 0 + offset_y), (1 + offset_x, 0 + offset_y)])
    instance.merge_cells(offset_cells(template, offset_x, offset_y))

def _two_cell_v_pair(instance: Automaton, offset_x: int = 0, offset_y: int = 0):
    '''add 2 vertical adjacent cells to existing universe'''
    template = [(0, 0), (0, 1)]
    instance.merge_cells(offset_cells(template, offset_x, offset_y))

def _three_cell_row(instance: Automaton, offset_x: int = 0, offset_y: int = 0):
    '''add 3 horizontal adjacent cells to existing universe'''
    # spinner
    template = [(0, 0), (1, 0), (2, 0)]
    instance.merge_cells(offset_cells(template, offset_x, offset_y))

def _three_cell_diagonal(instance: Automaton, offset_x: int = 0, offset_y: int = 0):
    '''add 3 horizontal adjacent cells to existing universe'''
    # spinner
    template = [(0, 0), (1, 1), (2, 2)]
    instance.merge_cells(offset_cells(template, offset_x, offset_y))

NEIGHBOURS_LIST = [
    (-1,-1),
    (-1,0),
    (-1,1),
    (0,1),
    (0,-1),
    (1,-1),
    (1,0),
    (1,1)
]
NEIGHBOURS = frozenset((
    (-1,-1),
    (-1,0),
    (-1,1),
    (0,1),
    (0,-1),
    (1,-1),
    (1,0),
    (1,1)
))
assert NEIGHBOURS == frozenset(NEIGHBOURS_LIST)

def my_main():
    '''wrapper for test/start code so that variables do not look like constants'''
    # old instantiation test lines
        # uni = AutomataUniverse(None)
        # uni = AutomataUniverse([])
        # uni = AutomataUniverse("tesr")
        # uni = AutomataUniverse("test")
        # uni = AutomataUniverse([tuple(), tuple()])
        # uni = AutomataUniverse([tuple(), (1,2)])
        # uni = AutomataUniverse([(-1,),(1,)]) # good «single argument»
        # uni = AutomataUniverse([(-1,-1),(1,1)]) # good «single argument»
        # uni = AutomataUniverse(NEIGHBOURS_LIST) # good «single argument»
        # uni = AutomataUniverse(frozenset(NEIGHBOURS_LIST)) # good «single argument»
        # uni = AutomataUniverse(NEIGHBOURS) # good «single argument»
        # uni = AutomataUniverse(NEIGHBOURS, None)
        # uni = AutomataUniverse(NEIGHBOURS, "test")
        # uni = AutomataUniverse(NEIGHBOURS, "step")
        # uni = AutomataUniverse(NEIGHBOURS, []) # good «two arguments»
        # uni = AutomataUniverse(NEIGHBOURS, [-1])
        # uni = AutomataUniverse(NEIGHBOURS, [9])
        # uni = AutomataUniverse(NEIGHBOURS, [0.5])
        # uni = AutomataUniverse(NEIGHBOURS, [0, 1, 2, 3, 4, 5, 6, 7, 8]) # good «two arguments»
        # uni = AutomataUniverse(NEIGHBOURS, [2,3]) # good «two arguments»
        # uni = AutomataUniverse(NEIGHBOURS, [2,3], [3]) # good
        # uni = AutomataUniverse(NEIGHBOURS, (2,3), (3,)) # good
        # uni = AutomataUniverse(NEIGHBOURS, set((2,3)), frozenset((3,))) # good
        # uni = AutomataUniverse(NEIGHBOURS, [2,3], [3]) # good
    uni = AutomataUniverse(NEIGHBOURS_LIST, [2,3], [3]) # good
    print("universe dimensions", uni.dimensions)

def my_main2():
    '''wrapper for test/start code so that variables do not look like constants'''
    # old instantiation test lines
        # amn = Automaton(None)
        # amn = Automaton(1)
        # amn = Automaton(frozenset())
        # amn = Automaton(frozenset((1,)))
        # amn = Automaton(frozenset((1,2)))
        # amn = Automaton(frozenset(((),)))
        # amn = Automaton(frozenset((tuple(),tuple('a'))))
        # amn = Automaton(frozenset((tuple(),(1,))))
        # amn = Automaton(frozenset(((0,),(1,))))
        # amn = Automaton(frozenset(((-2,),(1,))))
        # amn = Automaton(frozenset(((-1,),(1,)))) # good «single argument»
        # amn = Automaton(frozenset(((-1,-1),(1,1)))) # good «single argument»
        # amn = Automaton(NEIGHBOURS) # good «single argument»
        # amn = Automaton(NEIGHBOURS, None)
        # amn = Automaton(NEIGHBOURS, ('a',))
        # amn = Automaton(NEIGHBOURS, ('a', None))
        # amn = Automaton(NEIGHBOURS, (frozenset(), None))
        # amn = Automaton(NEIGHBOURS, (frozenset(), frozenset()))
        # amn = Automaton(NEIGHBOURS, (frozenset('a'), frozenset((None,))))
        # amn = Automaton(NEIGHBOURS, (frozenset((0,)), frozenset([99])))
        # amn = Automaton(NEIGHBOURS, (frozenset((2,)), frozenset([99])))

    amn = Automaton(NEIGHBOURS, (frozenset((2,3)), frozenset([3]))) # good «2 arguments»
    print("dimensions", amn.dimensions)
    print("neighbourhood", amn.neighbourhood)
    # print("extent", amn._get_extent(amn.neighbourhood))
    c_block = set(amn.neighbourhood)
    print("neighbourhood block", c_block)
    print("normalized result", amn._normalize_cells(c_block))
    print("normalized neighbourhood", c_block)

    # amn.iteration = 1
    # _one_cell(amn)
    # _one_cell(amn, 1, 2)
    # _two_cell_h_pair(amn)
    # _two_cell_h_pair(amn)
    # _two_cell_v_pair(amn)
    _three_cell_row(amn)
    # _three_cell_diagonal(amn)

    report_generation(amn)
    # amn.step()
    # report_generation(amn)
    # amn.step()
    # report_generation(amn)
    # amn.step()
    # report_generation(amn)


# Standalone module execution
if __name__ == "__main__":
    # my_main()
    my_main2()
