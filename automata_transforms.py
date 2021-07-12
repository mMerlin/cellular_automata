#!/usr/bin/env python
# coding=utf-8

"""
rotation and reflection transformation matrix management for cellular automata
"""

# pipenv shell

# standard library imports
from collections import namedtuple
from typing import Hashable, Iterable

# local application/library specific imports
import automata_typehints as AHint
from automata_universe import AutomataUniverse

XformSequence = namedtuple('TransformSequence', 'key seq')

class TransformSequence(XformSequence):
    """smarter namedtuple with parameter datatype validation

    tuples (namedtuple) are (with isinstance) Hashable, but will still fail to hash
    if a content element is not hashable. This mini class rejects elements that are
    not truly hashable.
    """

    def __new__(cls, key: Hashable, seq: int):
        """«pre» constructor

        :raises: TypeError
        """
        if not isinstance(seq, int):
            raise TypeError("not integer type: '{}'".format(type(seq).__name__))
        _test_hashable = hash(key) # only way to find out if 'really' hashable
        # self = super(TransformSequence, cls).__new__(cls, key, seq)
        # return self
        return super(TransformSequence, cls).__new__(cls, key, seq)
# end class TransformSequence()

def _base_index(key: Hashable) -> TransformSequence:
    return TransformSequence(key, 1)

def _check_transform_key(key: Hashable) -> None:
    """check that the key is valid to use for transform lookup

    :param key: unique identifier for stored transformation matrix
    :type key: any Hashable
    :raises: TypeError
    """
    _test_hashable = hash(key) # The only 'real' way to make sure is hashable
    # if not isinstance(key, Hashable):
    #     raise TypeError((type(key), "transformation lookup key is not hashable"))
# end def _check_transform_key()

class AutomataTransforms:
    """Storage for cellular automata transforms"""
    WorkArea = namedtuple('WorkArea', 'primes cycles combinations')

    def __init__(self, universe: AutomataUniverse) -> None:
        self._universe = universe
        self._prime_cells = self._build_transform_test_set()
        self._transform_cycles = dict()
        self._index_to_transform = dict()
        self._transform_to_index = dict()
        self._transform_cycles[None] = 0
        # self._neighbourhood_hash = hash(self._universe.neighbourhood)
        # instance shared work area
        self._scratch = self.WorkArea(dict(), dict(), dict())
        self._scratch.primes[self._prime_cells] = TransformSequence(None, 0)
        self._scratch.cycles[None] = 0
        self._scratch.combinations.clear()

    def __hash__(self) -> int:
        # hash of the set of all stored transforms, without key information
        # The transformation matrices themselves matter, not the lookup or count details
        return hash(frozenset(self._index_to_transform.values()))

    def _build_transform_test_set(self) -> AHint.CellGroupSnapshotType:
        some_prime_numbers = frozenset((2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53))
        prime_iter = iter(some_prime_numbers)
        cells = set()
        while True:
            cell_address = []
            try:
                for _number in range(self._universe.dimensions):
                    cell_address.append(next(prime_iter))
            except StopIteration:
                break
            cells.add(tuple(cell_address))
        return frozenset(cells)
    # end def _test_prime_set()

    def _check_transform_matrix(self, transform: AHint.TransformInputType) -> None:
        """check that the transformation matrix is structurally correct

        For a 'standard' universe, the extent should be inside ((-1,…),(1,…))
        -- each dimension extent should be either -1,0 or 0,1 ??
        -- a transform of any cell in the origin neighbourhood should be in the neighbourhood
        -- a transform of the origin neighbourhood should equal the neighbourhood

        :param transform: square rotation or reflection matrix
        :type transform: subscriptable sequence of «dimension» universe address tuples
        :raises: TypeError, ValueError
        """
        transform_test = self._universe.cell_group_transform(
            self._universe.neighbourhood, transform)
        if self._universe.neighbourhood != transform_test:
            raise ValueError(transform_test, "transformation of neighbourhood with "
                "{} is not the neighbourhood".format(transform))
        transform_test = self._universe.cell_group_transform(
            self._prime_cells, transform)
        if transform_test == self._prime_cells:
            raise ValueError("transform {} matches the identity matrix".format(transform))
    # end def _check_transform_matrix()

    def _make_transform_hashable(self, transform: AHint.TransformInputType) -> AHint.TransformType:
        """create a hashable version of a transformation matrix

        :param transform: square rotation or reflection matrix
        :type transform: sequence of «dimension» universe address tuples
        :raises: TypeError, ValueError
        """
        self._check_transform_matrix(transform)
        try_it = tuple(row for row in transform)
        _actual_hash = hash(try_it) # sanity check: should never actually raise exception
        return try_it
    # end def _make_transform_hashable()

    def _add_transform_setup(self, key: Hashable,
            requested_transform: AHint.TransformInputType) -> tuple[
                TransformSequence, AHint.TransformType]:
        """'setup before adding transformation matrices

        :param key: lookup key for cycle of transformation matrices
        :type key: any Hashable
        :param requested_transform: square rotation or reflection matrix
        :type requested_transform: sequence of «dimension» universe address tuples
        :raises: TypeError, ValueError
        """
        _check_transform_key(key)
        transform = self._make_transform_hashable(requested_transform)
        # print("hashable transform", transform, hash(transform),
        #     "extent", self._dbg_extent(transform)) # DEBUG
        if key in self._transform_cycles:
            raise ValueError("Cycle identifier '{}' already in use".format(key))
        cycle_index = _base_index(key)
        if cycle_index in self._index_to_transform:
            # Code logic error?: should never get here in normal usage
            raise ValueError("Base cycle index {} already in use".format(cycle_index))
        if transform in self._transform_to_index:
            raise ValueError("'{}' transform matches {}: {}".format(
                key, self._transform_to_index[transform], transform))
        return (cycle_index, transform)
    def add_transform_cycle(self, key: Hashable,
            requested_transform: AHint.TransformInputType) -> None:
        """add transform cycle to automaton

        :param key: lookup key for cycle of transformation matrices
        :type key: any Hashable
        :param requested_transform: square rotation or reflection matrix operation
        :type requested_transform: sequence of «dimension» universe address tuples
        :raises: TypeError, ValueError, ReferenceError
        """
        (cycle_index, transform) = self._add_transform_setup(key, requested_transform)
        new_transforms = 0
        new_patterns = 0
        transform_cycle = self._universe.identity_matrix
        transformed_pattern = self._universe.cell_group_transform(self._prime_cells, transform)
        while True: # tail recursion equivalent
            if transformed_pattern == self._prime_cells:
                self._end_cycle_reached((cycle_index, new_transforms, new_patterns),
                    transform_cycle, transform)
                break
            if transformed_pattern not in self._scratch.primes:
                self._scratch.primes[transformed_pattern] = cycle_index
                new_patterns += 1
            else: # DEBUG normal should just continue
                raise ValueError("{} pattern for {} matches {}".format(transformed_pattern,
                    cycle_index, self._scratch.primes[transformed_pattern]))

            transform_cycle = self._universe.matrix_transform(transform, transform_cycle)
            if transform_cycle == self._universe.identity_matrix:
                raise ValueError("transformation for {} matches the identity matrix".format(
                    cycle_index))
            test_pattern = self._universe.cell_group_transform(self._prime_cells, transform_cycle)
            if transformed_pattern != test_pattern:
                # transform_cycle is supposed to be a one step operation creating
                # transformed_pattern from self._prime_cells
                cycle_start = _base_index(cycle_index.key)
                # this is really a code logic error. Code is broken to get here
                raise ReferenceError((transformed_pattern, test_pattern, "{} transform {} does not "
                    "produce the same result as multiple operations of the {} transform {}".format(
                    cycle_index, transform_cycle, cycle_start,
                    self._index_to_transform[cycle_start])))

            if transform_cycle not in self._transform_to_index:
                self._transform_to_index[transform_cycle] = cycle_index
                new_transforms += 1
            else: # DEBUG normally should just continue
                raise ValueError("index {} already matches transform {}".format(
                    self._transform_to_index[transform_cycle], transform_cycle))

            self._index_to_transform[cycle_index] = transform_cycle
            transformed_pattern = self._universe.cell_group_transform(
                transformed_pattern, transform)
            cycle_index = TransformSequence(key, cycle_index.seq + 1)
            if cycle_index in self._index_to_transform:
                raise ValueError("cycle index {} already in use".format(cycle_index))
            if cycle_index.seq > 10:
                # Should be impossible to get here. The initial validation of a transformation
                # matrix should guarantee that it cycles in a few iterations. At least for a
                # small number of dimensions
                raise RecursionError("cycle {} did not wrap to the beginning before {} steps "
                    "starting from {}".format(cycle_index.key, cycle_index.seq, transform))
        # end while True
    # end def add_transform_cycle()
    def _end_cycle_reached(self, context: tuple[TransformSequence, int, int],
            cycle_result: AHint.TransformType, operation: AHint.TransformType) -> None:
        """save information for the completed cycle

        :param context: data accumulated for cycle
        :type context: tuple containing
            identifier index for final transformation matrix in cycle
            number of new transformation matrices found for cycle
            number of new patterns found for cycle
        :param cycle_result: matrix transformation for final step in cycle
        :type cycle_result: tuple of «dimension» cellular automata cell address tuples
        :param operation: matrix transformation that advances between cycle steps
        :type operation: tuple of «dimension» cellular automata cell address tuples
        """
        (end_cycle, added_transforms, added_patterns) = context
        wrap_transform = self._universe.matrix_transform(operation, cycle_result)
        if wrap_transform != self._universe.identity_matrix:
            raise ValueError("pattern ended at {} before transform cycled: {}".format(
                end_cycle, cycle_result))
        if added_transforms == 0:
            # safety net: logic should never be able to reach here
            raise TypeError("no new transformations for '{}' starting from {}".format(
                end_cycle.key, operation))
        if added_patterns == 0:
            raise TypeError("no new patterns for '{}' starting from {}".format(
                end_cycle.key, operation))
        self._transform_cycles[end_cycle.key] = end_cycle.seq
    # end def _end_cycle_reached()

    def _verify_new_transform(self, transform: AHint.TransformType,
            target_pattern: AHint.CellGroupWorkingType,
            sequence: tuple[TransformSequence, TransformSequence]) -> None:
        """sanity check a new transformation matrix and expected pattern

        :param transform: single step transformation matrix
        :type transform: tuple of «dimension» cellular automata cell address tuples
        :param target_pattern: expected pattern after transforming base cells
        :type target_pattern: set of cellular automata cell address tuples
        :param sequence: identifers for equivalent sequential tranformation matrices
        :type sequence: tuple of transformation cycle identifiers
        :raises: ValueError, ReferenceError
        """
        if transform in self._index_to_transform:
            raise ValueError((sequence, transform,
                "new pattern with old transform"))
        test_pattern = self._universe.cell_group_transform(self._prime_cells, transform)
        if test_pattern != target_pattern:
            # this is really a code logic error. Code is broken to get here
            raise ReferenceError((target_pattern, test_pattern, "{} transform does not "
                "produce the same result as sequential operation of {}".format(
                transform, sequence)))
    def _combinations_from_cycle(self, source_index: TransformSequence,
            source_pattern: AHint.CellGroupSnapshotType,
            cycle: Hashable, count: int) -> None:
        """locate new unique patterns using transfrom from cycle

        :param source_index:
        :type source_index: tuple containing cycle key and sequence number
        :param source_pattern: existing pattern after previous transforms
        :type source_pattern: frozenset of cellular automata universe coordinates
        :param cycle: identifier for transformation matrix cycle
        :type cycle: Hashable
        :param count: number of steps in cycle
        :type count: int
        :raises: ValueError, ReferenceError «via _verify_new_transform»
        """
        work_pattern = source_pattern
        work_transform = self._index_to_transform[source_index]
        cycle_transform = self._index_to_transform[_base_index(cycle)]
        for sequence in range(1, count):
            combined_index = (source_index, TransformSequence(cycle, sequence))
            work_pattern = self._universe.cell_group_transform(work_pattern, cycle_transform)
            work_transform = self._universe.matrix_transform(cycle_transform, work_transform)
            if work_pattern in self._scratch.primes:
                continue # not a new pattern
            self._verify_new_transform(work_transform, work_pattern, combined_index)
            if work_pattern not in self._scratch.combinations:
                self._scratch.combinations[work_pattern] = (work_transform, [])
            if self._scratch.combinations[work_pattern][0] != work_transform:
                raise ValueError((work_transform, source_index, (cycle, sequence),
                    "existing pattern with new transform"))
            self._scratch.combinations[work_pattern][1].append(combined_index)
    IDENTITY_KEY = TransformSequence(None, 0)
    def _combinations_after_pattern(self, starting_pattern: AHint.CellGroupSnapshotType,
            starting_index: TransformSequence) -> None:
        """locate new unique patterns transformed from specified pattern

        :param starting_pattern: existing pattern after previous transforms
        :type starting_pattern: frozenset of cellular automata universe coordinates
        :param starting_index:
        :type starting_index: tuple containing cycle key and sequence number
        :raises: ValueError, ReferenceError
        """
        # only transforms from different cycles are really valid to produce new combinations
        # this ignores that, by just filtering them out based on duplication of an existing
        # pattern.
        if starting_index == self.IDENTITY_KEY:
            return # skip the identity transformation case
        # limit to transforms for base_index, but repeat for cycle length
        for cycle, count in self._transform_cycles.items():
            # print("cycle {}, length {}".format(cycle, count)) # DEBUG
            if count < 2:
                continue
            self._combinations_from_cycle(starting_index, starting_pattern, cycle, count)
    def _find_combination_targets(self) -> None:
        """locate new patterns for loaded transforms

        :raises: ValueError, ReferenceError
        """
        for cell_pattern, pattern_index in self._scratch.primes.items():
            self._combinations_after_pattern(cell_pattern, pattern_index)
    def _unique_combination_index(self, pattern: AHint.CellGroupSnapshotType,
            sequence: Iterable[tuple[TransformSequence, ...]],
            transform: AHint.TransformType) -> frozenset[Hashable]:
        """create index for unique transformation sequence

        :param pattern: new unique transformation result pattern
        :type pattern: frozenset of cellular automata universe coordinates
        :param sequence: sequences of cycle index keys that (each) produce pattern
        :type sequence: iterable of tuples of TransformationSequence instances
        :param transform: transformation matrix that directly produces pattern from base
        :type transform: tuple of cellular automata universe coordinates
        :raises: ValueError, TypeError «unhashable type»
        """
        sequence_key = frozenset(sequence)
        # print("routes", len(path_key)) # DEBUG
        # for path in path_key:
        #     print((tuple(path[0]), tuple(path[1]))) # DEBUG
        if sequence_key in self._transform_cycles:
            raise ValueError(("generated combined key already has cycle", sequence_key))
        if sequence_key in self._index_to_transform:
            raise ValueError(("generated combined key already in use", sequence_key))
        if transform in self._transform_to_index:
            raise ValueError(("generated transform already in use",
                transform, sequence))
        if pattern in self._scratch.primes:
            raise ValueError(("generated pattern already in use", pattern))
        assert isinstance(pattern, frozenset)
        return sequence_key
    def _add_unique_combinations(self) -> None:
        """add paths for combination patters that are unique

        :raises: ValueError, TypeError
        """
        for single_target, (combined_transform, paths) in self._scratch.combinations.items():
            path_key = self._unique_combination_index(single_target, paths, combined_transform)
            self._transform_cycles[path_key] = 0
            self._scratch.primes[single_target] = path_key
            self._transform_to_index[combined_transform] = path_key
            self._index_to_transform[path_key] = combined_transform
    def generate_combination_transforms(self) -> None:
        """locate additional transform possibilities by combining existing transforms

        :raises: ValueError, TypeError, ReferenceError
        """
        if len(self._index_to_transform) == 0:
            raise ValueError(
                "No transforms have been added yet. Nothing to base generated transforms on")
        if len(self._transform_cycles) == 2:
            cycles = [cycle for cycle in self._transform_cycles if cycle is not None]
            assert len(cycles) == 1
            raise ValueError("Only a single transform cycle '{}' exists. Nothing to base "
                "combinations on".format(cycles[0]))
        self._scratch.combinations.clear()
        # self._combination.transforms.clear()
        self._find_combination_targets()
        self._add_unique_combinations()
    # end def generate_combination_transforms()
# end class AutomataTransforms


SQUARE_GRID_NEIGHBORS = [
    (-1,-1), (-1,0), (-1, 1),
    ( 0, 1),         ( 0,-1),
    ( 1,-1), ( 1,0), ( 1, 1)
]

def my_main() -> None: # pragma: no cover
    """wrapper for test/start code so that variables do not look like constants"""
    universe = AutomataUniverse(SQUARE_GRID_NEIGHBORS, [2,3], [3])
    instance = AutomataTransforms(universe)
    assert isinstance(instance, AutomataTransforms)
    # # _is_rot_mat_test(instance)
    # # _rotations_check(instance)
    # # _prime_cells_check(instance)
    # _check_transform_test(instance)
    # # _hashable_transform_test(instance)
    # _duplicate_test(instance)
    # _collision_test(instance)
    # _end_cycle_test(instance)
    # _add_transform_test(instance)
    # instance.generate_combination_transforms()

    # # _matrix_rotate_test(instance)
    # # _duplicate_test(instance) # test again after transform(s) added
    # # _collision_test(instance) # test again after transform(s) added «also refactoring»
    # instance.dbg_report_instance() # DEBUG

# Standalone module execution
if __name__ == "__main__": # pragma: no cover
    my_main()
