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
from automaton import Automaton

# validate data against type hints
# https://stackoverflow.com/questions/50563546/validating-detailed-types-in-python-dataclasses

# setup functions to add patterns of cell to and automaton universe generation

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

def my_main():
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
        # amn = Automaton(neighbours) # good «single argument»
        # amn = Automaton(neighbours, None)
        # amn = Automaton(neighbours, ('a',))
        # amn = Automaton(neighbours, ('a', None))
        # amn = Automaton(neighbours, (frozenset(), None))
        # amn = Automaton(neighbours, (frozenset(), frozenset()))
        # amn = Automaton(neighbours, (frozenset('a'), frozenset((None,))))
        # amn = Automaton(neighbours, (frozenset((0,)), frozenset([99])))
        # amn = Automaton(neighbours, (frozenset((2,)), frozenset([99])))

    neighbours = frozenset((
        (-1,-1),
        (-1,0),
        (-1,1),
        (0,1),
        (0,-1),
        (1,-1),
        (1,0),
        (1,1)
    ))
    amn = Automaton(neighbours, (frozenset((2,3)), frozenset([3]))) # good «2 arguments»
    print("dimensions", amn.dimensions)
    # amn.iteration = 1
    # _one_cell(amn)
    # _one_cell(amn, 1, 2)
    # _two_cell_h_pair(amn)
    # _two_cell_h_pair(amn)
    # _two_cell_v_pair(amn)
    _three_cell_row(amn)
    # _three_cell_diagonal(amn)

    print("generation", amn.iteration, "is", list(amn.generation))
        # print(amn.generation)
        # for cell in amn.generation:
        #     print(cell)
    amn.step()
    print("generation", amn.iteration, "is", list(amn.generation))
    amn.step()
    print("generation", amn.iteration, "is", list(amn.generation))


# Standalone module execution
if __name__ == "__main__":
    my_main()
