#!/usr/bin/env python
# coding=utf-8

"""typehint aliases for the Automata¦Automaton families of classes"""

# pipenv shell

# standard library imports
# from typing import Tuple
from collections.abc import Iterable
from typing import Union

# validate data against type hints
# https://stackoverflow.com/questions/50563546/validating-detailed-types-in-python-dataclasses

CellAddressType = tuple[int, ...]
NeighbourhoodInputType = Iterable[CellAddressType]
NeighbourhoodType = frozenset[CellAddressType]
PropagationRuleInputType = Iterable[int]
PropagationRuleType = frozenset[int]
CellGroupWorkingType = set[CellAddressType]
CellGroupSnapshotType = frozenset[CellAddressType]
CellGroupType = Union[CellGroupWorkingType, CellGroupSnapshotType]
BoundingBoxType = tuple[CellAddressType, CellAddressType]
CellorCellsType = Union[Iterable[CellAddressType], CellAddressType]
TransformInputType = Iterable[CellAddressType]
TransformType = tuple[CellAddressType, ...]
RotateReflectInputType = Iterable[TransformInputType]

# def my_main() -> None:
#     """wrapper for test/start code so that variables do not look like constants"""

# # Standalone module execution
# if __name__ == "__main__":
#     my_main()
