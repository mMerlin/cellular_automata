#!/usr/bin/env python
# coding=utf-8

'''
rotation and reflection transformation matrix managate for cellular automata
'''

# pipenv shell

# standard library imports
# from collections.abc import Iterable
from collections import namedtuple
from typing import Hashable
# , Sequence
import math
# from main import report_generation

# local application/library specific imports
import automata_typehints as AHint
from automata_universe import AutomataUniverse

TransformSequence = namedtuple('TransformSequence', 'key seq')

class AutomataTransforms:
    '''Storage for cellular automata transforms'''

    def __init__(self, universe: AutomataUniverse) -> None:
        self._universe = universe
        self._prime_cells = self._build_transform_test_set()
        self._transform_cycles = dict()
        self._index_to_transform = dict()
        self._transform_to_index = dict()
        # self._neighbourhood_hash = hash(self._universe.neighbourhood)
        # instance shared work area
        self._cycle_cells = None
        self._primes_transformed = dict()
        self._primes_transformed[self._prime_cells] = TransformSequence(None, 0)
        self._transform_cycles[None] = 0

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
        '''check that the transformation matrix is structurally correct

        For a 'standard' universe, the extent should be inside ((-1,…),(1,…))
        -- each dimension extent should be either -1,0 or 0,1 ??
        -- a transform of any cell in the origin neighbourhood should be in the neighbourhood
        -- a transform of the origin neighbourhood should equal the neighbourhood

        :param transform: square rotation or reflection matrix
        :type transform: subscriptable sequence of «dimension» universe address tuples
        :raises: TypeError, ValueError
        '''
        self._universe.validate_matrix(transform)
        transform_test = self._universe.cell_group_transform(self._universe.neighbourhood,
            transform)
        if self._universe.neighbourhood != transform_test:
            print("neighbourhood", self._universe.neighbourhood) # DEBUG
            print("transformed", transform_test) # DEBUG
            raise ValueError(transform_test, "transformation of neighbourhood with "
                "{} is not the neighbourhood".format(transform))
        transform_test = self._universe.cell_group_transform(self._prime_cells,
            transform)
        if transform_test == self._prime_cells:
            print("primes cells", type(self._prime_cells),
                self._prime_cells) # DEBUG
            print("transformed cell", type(transform_test), transform_test) # DEBUG
            raise ValueError("transform {} matches the identity matrix".format(transform))
    # end def _check_transform_matrix()

    def _make_transform_hashable(self, transform: AHint.TransformInputType) -> AHint.TransformType:
        '''create a hashable version of a transformation matrix

        :param transform: square rotation or reflection matrix
        :type transform: sequence of «dimension» universe address tuples
        :raises: TypeError
        '''
        self._check_transform_matrix(transform)
        return tuple(row for row in transform)
    # end def _make_transform_hashable()

    @staticmethod
    def _base_index(key: Hashable) -> TransformSequence:
        return TransformSequence(key, 1)

    @staticmethod
    def _check_transform_key(key: Hashable) -> None:
        '''check that the key is valid to use for transform lookup

        :param key: unique identifier for stored transformation matrix
        :type key: any Hashable
        :raises: TypeError
        '''
        if not isinstance(key, Hashable):
            raise TypeError((type(key), "transformation lookup key is not hashable"))
    # end def _check_transform_key()

    def _check_duplication(self, key: Hashable, transform: AHint.TransformType) -> None:
        '''check that neither key or transform are already being used

        :param key: unique identifier for stored transformation matrix
        :type key: any Hashable
        :param transform: square rotation or reflection matrix
        :type transform: sequence of «dimension» universe address tuples
        :raises: ValueError
        '''
        if key in self._transform_cycles:
            raise ValueError("transform key {} already exists".format(key))
        if transform in self._transform_to_index:
            index = self._transform_to_index[transform]
            raise ValueError("matrix {} already exists in transformations as sequence "
                "{} of '{}'".format(transform, index.seq, index.key))
    # end def _check_duplication()

    def _check_transform_collision(self, sequence: int, index: TransformSequence,
            transform: AHint.TransformType) -> None:
        '''check if the cycle transform collided with a different cycle

        :param sequence: sequence number of the cycle being checked
        :type sequence: int
        :param index: lookup key for the first (base) entry of the cycle being checked
        :type index: TransformSequence namedtuple
        :param transform: cycle sequence transformation matrix
        :type transform: tuple of «dimension» cell address tuples
        :raises: ValueError
        '''
        if transform in self._transform_to_index:
            existing_index = self._transform_to_index[transform]
            if existing_index != index:
                print("\ncollision", index.key, sequence, existing_index, transform)
                raise ValueError("Sequence {} of {}: {} intersected existing sequence {} of {} "
                    "cycle".format(sequence, index.key, transform,
                    existing_index.seq, existing_index.key))
    # end def _check_transform_collision()

    def _check_end_cycle(self, key: Hashable, sequence: int,
            transform: AHint.TransformType) -> bool:
        '''check if the end of a transform cycle (loop) has been reached

        :param key: unique identifier for stored transformation matrix
        :type key: any Hashable
        :param sequence: sequence number of the cycle being checked
        :type sequence: int
        :param transform: cycle sequence transformation matrix
        :type transform: tuple of cell address tuples
        :raises: ValueError «in _check_transform_collision»
        '''
        # print("check end cycle '{}' {}: {}".format(key, sequence, transform)) # DEBUG
        base_idx = self._base_index(key)
        self._cycle_cells = self._universe.cell_group_transform(
            self._cycle_cells, self._index_to_transform[base_idx])
        print("cycle cells", self._cycle_cells, "\n -hash", hash(self._cycle_cells)) # DEBUG
        if self._cycle_cells == self._prime_cells:
            self._check_transform_collision(sequence, base_idx, transform)
            print("is end cycle") # TRACE
            return True
        if not transform in self._transform_to_index:
            return False
        existing_index = self._transform_to_index[transform]
        raise ValueError("Sequence {} of {}: {} intersected existing sequence {} of {} "
            "cycle".format(sequence, key, self._index_to_transform[base_idx],
            existing_index.seq, existing_index.key))
    # end def _check_end_cycle()

    def _add_cycle_references(self, key: Hashable, sequence: int,
            transform: AHint.TransformType) -> None:
        '''recursively create keys for remaining transformations in a cycle

        :param key: unique identifier for stored transformation matrix cycle
        :type key: any Hashable
        :param sequence: sequence number in cycle
        :type sequence: int
        :param transform: square rotation or reflection matrix
        :type transform: sequence of «dimension» universe address tuples
        '''
        self._cycle_cells = self._universe.cell_group_transform(self._cycle_cells, transform)
        if self._cycle_cells == self._prime_cells:
            # end of cycle, just save the cycle length
            self._transform_cycles[key] = sequence
            return
        if sequence > 10:
            raise ValueError("transform did not cycle in {} steps for {} = {}".format(
                sequence, key, self._index_to_transform[self._base_index(key)]))
        if self._cycle_cells in self._primes_transformed:
            raise NotImplementedError("not ready for duplicate cycle transform")
        cycle_key = TransformSequence(key, sequence)
        self._index_to_transform[cycle_key] = None
        self._primes_transformed[self._cycle_cells] = cycle_key
        self._add_cycle_references(key, sequence + 1, transform)
    # end def _add_cycle_references()

    def add_transform_cycle(self, key: Hashable,
            requested_transform: AHint.TransformInputType) -> None:
        '''add transform cycle to automaton

        :param key: lookup key for cycle of transformation matrices
        :type key: any Hashable
        :param requested_transform: square rotation or reflection matrix
        :type requested_transform: sequence of «dimension» universe address tuples
        :raises: TypeError, ValueError
        '''
        self._check_transform_key(key)
        transform = self._make_transform_hashable(requested_transform)
        # print("hashable transform", transform, hash(transform),
        #     "extent", self._dbg_extent(transform)) # DEBUG
        if key in self._transform_cycles:
            raise ValueError("Cycle identifier '{}' already in use".format(key))
        base_index = self._base_index(key)
        if base_index in self._index_to_transform:
            # Code logic error: should never get here
            raise ValueError("Base cycle index {} already in use".format(base_index))
        if transform in self._transform_to_index:
            raise ValueError("Transform matches {}".format(self._transform_to_index[transform]))

        self._cycle_cells = self._universe.cell_group_transform(self._prime_cells, transform)
        if self._cycle_cells in self._primes_transformed:
            prev_source = self._primes_transformed[self._cycle_cells]
            # print("starting transform for '{}' matched previous cycle instance {}:".format(
            #     key, prev_source), transform) # DEBUG
            prev_transform = self._index_to_transform[prev_source]
            if prev_transform is not None:
                raise ValueError("{} cycle transform already populated with {}".format(
                    prev_source, prev_transform))
            self._index_to_transform[prev_source] = transform
            self._transform_to_index[transform] = prev_source
            return
        self._index_to_transform[base_index] = transform
        self._transform_to_index[transform] = base_index
        self._primes_transformed[self._cycle_cells] = base_index
        # print("start of new cycle", base_index) # DEBUG
        self._add_cycle_references(key, base_index.seq + 1, transform)
    # end def add_transform()

    IDENTITY_KEY = TransformSequence(None, 0)
    def generate_combination_transforms(self) -> None:
        '''locate additional transform possibilities by combining existing transforms'''
        new_targets = dict()
        for cell_pattern, first in self._primes_transformed.items():
            # print("cell pattern", cell_pattern, first) # DEBUG
            # print("first", first) # DEBUG
            if first == self.IDENTITY_KEY:
                continue # skip the identity transformation case
            # limit to base_index, but repeat for cycle length
            for cycle, count in self._transform_cycles.items():
                # print("cycle {}, length {}".format(cycle, count)) # DEBUG
                if count < 2:
                    continue
                transform = self._index_to_transform[self._base_index(cycle)]
                self._cycle_cells = cell_pattern
                for sequence in range(1, count):
                    combined = (first, TransformSequence(cycle, sequence))
                    # print("sequence", sequence, combined) # DEBUG
                    self._cycle_cells = self._universe.cell_group_transform(
                        self._cycle_cells, transform)
                    if self._cycle_cells in self._primes_transformed:
                        continue # not a new result
                    # print("NEW", combined) # DEBUG
                    if self._cycle_cells not in new_targets:
                        new_targets[self._cycle_cells] = list()
                    new_targets[self._cycle_cells].append(combined)
        generated_instance = 0
        for new_target, paths in new_targets.items():
            # paths = new_targets[new_target]
            path_key = frozenset(paths)
            print("routes", len(path_key)) # DEBUG
            for path in path_key:
                print((tuple(path[0]), tuple(path[1]))) # DEBUG
            # transform_key = self._base_index(path_key)
            generated_instance += 1
            transform_key = "generated {}".format(generated_instance)
            while transform_key in self._index_to_transform:
                generated_instance += 1
                transform_key = "generated {}".format(generated_instance)
            if path_key in self._transform_cycles:
                raise ValueError(("generated combined key already in use", path_key))
            if new_target in self._index_to_transform:
                raise ValueError(("generated transform already in use", new_target))
            # merge into self._primes_transformed
            self._primes_transformed[new_target] = transform_key
            self._index_to_transform[transform_key] = None
    # end def generate_combination_transforms()

    # DEBUG methods

    @property
    def dbg_prime_cells(self) -> AHint.CellGroupSnapshotType:
        return self._prime_cells
    def dbg_check_transform_matrix(self, transform: AHint.TransformInputType) -> None:
        return self._check_transform_matrix(transform)
    def dbg_make_transform_hashable(self,
            transform: AHint.TransformInputType) -> AHint.CellGroupSnapshotType:
        return self._make_transform_hashable(transform)
    def dbg_base_index(self, key: Hashable) -> TransformSequence:
        return self._base_index(key)
    # def dbg_check_transform_key(self, key: Hashable) -> None:
        # return self._check_transform_key(key)
    def dbg_check_duplication(self, key: Hashable, transform: AHint.TransformType) -> None:
        return self._check_duplication(key, transform)
    def dbg_check_transform_collision(self, sequence: int, index: TransformSequence,
            transform: AHint.TransformType) -> None:
        return self._check_transform_collision(sequence, index, transform)
    def dbg_check_end_cycle(self, key: Hashable, sequence: int,
            transform: AHint.TransformType) -> bool:
        return self._check_end_cycle(key, sequence, transform)
    def dbg_add_cycle_references(self, key: Hashable, sequence: int,
            transform: AHint.TransformType) -> None:
        return self._add_cycle_references(key, sequence, transform)

    def dbg_is_rot_mat(self, matrix: AHint.TransformInputType) -> bool:
        return self._universe.is_rotation_matrix(matrix)
    def dbg_data_rotate(self, cells: AHint.CellGroupType,
            transform: AHint.TransformInputType) -> AHint.CellGroupSnapshotType:
        return self._universe.cell_group_transform(cells, transform)
    def dbg_mat_rotate(self, matrix: AHint.TransformInputType,
            transform: AHint.TransformInputType) -> AHint.TransformType:
        return self._universe.matrix_transform(matrix, transform)
    def dbg_report_instance(self) -> None:
        '''dump instance data'''
        print("\ninstance dump")
        for transform, index in self._transform_to_index.items():
            print("transform to index entry", index, transform) # DEBUG
        for index, transform in self._index_to_transform.items():
            print("index to transform entry", index, transform) # DEBUG
        print("cycles", self._transform_cycles) # DEBUG
        print("transformed primes")
        for prime_cells, index in self._primes_transformed.items():
            print(index, prime_cells)
        print("cycled cells", self._cycle_cells)
    def _dbg_extent(self, cells: AHint.CellGroupWorkingType) -> AHint.BoundingBoxType: # DEBUG
        '''determine the n-dimensional bounding box for a set of cells

        :param cells:
        :type cells: «forzen»set of cell coordinate tuples
        :returns: minimum and maximum corner coordinates of bounding box
        :rtype: tuple[tuple[int, ...],tuple[int, ...]]
        '''
        box_min = [float('inf')] * self._universe.dimensions
        box_max = [-box_min[0]] * self._universe.dimensions
        for coord in cells:
            box_min = [min(new, cur) for new, cur in zip(coord, box_min)]
            box_max = [max(new, cur) for new, cur in zip(coord, box_max)]
        return tuple((tuple(box_min), tuple(box_max)))
# end class Transforms


SQUARE_GRID_NEIGHBORS = [
    (-1,-1), (-1,0), (-1, 1),
    ( 0, 1),         ( 0,-1),
    ( 1,-1), ( 1,0), ( 1, 1)
]
HORIZONTAL_REFLECTION = [(-1, 0), (0, 1)]
VERTICAL_REFLECTION = [(1, 0), (0, -1)]
OPPOSITE_REFLECTION = [(-1, 0), (0, -1)]

def _angle_to_matrix(degrees_angle: int) -> AHint.TransformType:
    radians_angle = math.radians(degrees_angle)
    return (
        (int(math.cos(radians_angle)), -int(math.sin(radians_angle))),
        (int(math.sin(radians_angle)), int(math.cos(radians_angle)))
    )
ROTATIONS = tuple(_angle_to_matrix(rot_angle) for rot_angle in [0, 90, 180, 270, 360])
# ROTATE90 = [(0, -1), (1, 0)]
ROTATE90 = ROTATIONS[1]

def transforms_create_test() -> AutomataTransforms:
    universe = AutomataUniverse(SQUARE_GRID_NEIGHBORS, [2,3], [3])
    return AutomataTransforms(universe)
def _is_rot_mat_test(instance: AutomataTransforms) -> None: # DEBUG
    '''exercise is_rot_mat'''
    # folder marker for tst_mat cases that raise exceptions
        # tst_mat = None
        # tst_mat = 1
        # tst_mat = []
        # tst_mat = [3,7]
    # instance.dbg_is_rot_mat(tst_mat)
    false_tests = [
        [(1, 5),(3, 4)],
        [(1, 1),(1, 1)],
        [(0, 1),(0, 1)],
        [(0, 0),(0, 1)]
    ]
    for test_case in false_tests:
        assert instance.dbg_is_rot_mat(test_case) is False
    for rot in ROTATIONS:
        assert instance.dbg_is_rot_mat(rot) is True
def _rotations_check(instance: AutomataTransforms) -> None: # DEBUG
    print("ROTATIONS entries")
    for rot in ROTATIONS:
        print(rot)
    print("is_rot_matrix ?", tuple(instance.dbg_is_rot_mat(rot) for rot in ROTATIONS))
def _prime_cells_check(instance: AutomataTransforms) -> None: # DEBUG
    print("prime test", instance.dbg_prime_cells, "\n -hash",
        hash(instance.dbg_prime_cells))
def _check_transform_test(instance: AutomataTransforms) -> None:
    '''exercise _check_transform_matrix'''
    # folder marker for tst_mat cases
        # tst_mat = None
        # tst_mat = []
        # tst_mat = [1, 2]
        # tst_mat = [(1, 2),(3, 4)]
        # tst_mat = ROTATIONS[0] # Identity, correctly fails
    tst_mat = ROTATIONS[1]
    instance.dbg_check_transform_matrix(tst_mat)
    tst_mat = HORIZONTAL_REFLECTION
    instance.dbg_check_transform_matrix(tst_mat)
def _hashable_transform_test(instance: AutomataTransforms) -> None:
    '''exercise validate_matrix, _make_transform_hashable'''
    # folder marker for tst_mat cases
        # tst_mat = set([(0, -1), (1, 0)]) # correctly fails not subscriptable
        # tst_mat = frozenset([(0, -1), (1, 0)]) # correctly fails not subscriptable
        # tst_mat = [(0, -1), (1, 0)] # good
    tst_mat = ((0, -1), (1, 0)) # good
    tst_result = instance.dbg_make_transform_hashable(tst_mat)
    assert tst_mat == tst_result
    print(tst_mat, "as hashable =", tst_result)
def _duplicate_test(instance: AutomataTransforms) -> None:
    '''exercise _check_duplication'''
    tst_transform = ROTATIONS[1] # expects already validated (hashable) matrix
    # tst_key fold marker
        # tst_key = None
        # tst_key = []
        # tst_key = 1 # Good
    tst_key = "Rotate 90°"
    instance.dbg_check_duplication(tst_key, tst_transform)
    # verify and add more checks once transformation matrices added
def _collision_test(instance: AutomataTransforms) -> None:
    '''exercise _check_transform_collision'''
    tst_mat = ROTATIONS[0] # Identity test
    tst_seq = 2
    tst_idx = TransformSequence(None, 0)
    instance.dbg_check_transform_collision(tst_seq, tst_idx, tst_mat)
    tst_idx = instance.dbg_base_index(None)
    instance.dbg_check_transform_collision(tst_seq, tst_idx, tst_mat)
    tst_mat = ROTATIONS[1] # Rotate +90°
    tst_idx = instance.dbg_base_index(999)
    instance.dbg_check_transform_collision(tst_seq, tst_idx, tst_mat)
    tst_idx = instance.dbg_base_index("Rotate 90°")
    instance.dbg_check_transform_collision(tst_seq, tst_idx, tst_mat)
    # verify and add more checks once transformation matrices added
def _end_cycle_test(instance: AutomataTransforms) -> None:
    '''exercise _check_end_cycle'''
    tst_key = 998
    tst_seq = 2
    tst_mat = ROTATIONS[2] # validated hashable matrix expected
    # print("check end cycle '{}' {}: {}".format(tst_key, tst_seq, tst_mat)) # DEBUG
    try:
        print(instance.dbg_check_end_cycle(tst_key, tst_seq, tst_mat))
        print("ERROR: should have raised KeyError")
    except KeyError as ex:
        # print("got expected KeyError:", ex)
        assert ex.args == (instance.dbg_base_index(tst_key),)
    # add more checks once transformation matrices added
def _add_transform_test(instance: AutomataTransforms) -> None:
    '''exercise add_transfrom'''
    # tst_key = None
    tst_key = "Rotate 90°"
    # tst_mat = ROTATIONS[0] # identity matrix correctly fails
    tst_mat = ROTATIONS[1]
    instance.add_transform_cycle(tst_key, tst_mat)
    # instance.add_transform_cycle(tst_key, tst_mat)
    tst_mat = HORIZONTAL_REFLECTION
    tst_key = "Horizontal Reflection"
    instance.add_transform_cycle(tst_key, tst_mat)
    # tst_key = None
    # tst_key = "filler"
    # tst_mat = ROTATIONS[2]
    # instance.add_transform_cycle(tst_key, tst_mat)
    # # instance.add_transform_cycle(tst_key, tst_mat)
    # tst_mat = ROTATIONS[3]
    # instance.add_transform_cycle(tst_key, tst_mat)
    tst_mat = VERTICAL_REFLECTION
    tst_key = "Vertical Reflection"
    instance.add_transform_cycle(tst_key, tst_mat)
    # tst_mat = OPPOSITE_REFLECTION
    # tst_key = "Opposite Reflection"
    # instance.add_transform_cycle(tst_key, tst_mat)

def _matrix_rotate_test(instance: AutomataTransforms) -> None:
    '''exercise transform validation and setup'''
    print("\nmatrix test")
    # fold marker
        # # matrix test data cases
        #     # tr_key = []
        #     # tr_mat = None
        #     # tr_mat = tuple()
        #     # tr_mat = (5,7)
        #     # tr_mat = (tuple(),(7,5,3))
        #     # tr_mat = ((None, None),(7,5,3))
        #     # tr_mat = ((0, -1),(7,5,3))
        #     # tr_mat = [(1, 0), (0, 1)] # identity matrix
        # # tr_mat = [(1, -1), (1, 0)]
        # tr_mat = ROTATE90 # GOOD
        # # matrix unique identifier cases
        #     # tr_key = []
        #     # tr_key = None # GOOD
        #     # tr_key = 1 # GOOD
        # tr_key = "Rotate 90°"
        # instance.add_transform(tr_key, tr_mat)
        # # tr_key = 1
        # # # tr_mat = ROTATE90
        # # tr_mat = [(0, 1), (-1, 0)]
        # tr_mat = HORIZONTAL_REFLECTION # GOOD
        # tr_key = "Horizontal Reflection"
        # instance.add_transform(tr_key, tr_mat)
    tst_mat1 = ((0, -1), (1, 0))
    tst_dat1 = instance.dbg_data_rotate(instance.dbg_prime_cells, tst_mat1)
    print(instance.dbg_prime_cells)
    print(tst_mat1)
    print(tst_dat1)
    tst_dat2 = instance.dbg_data_rotate(tst_dat1, tst_mat1)
    print(tst_dat1)
    print(tst_mat1)
    print(tst_dat2)
    tst_mat2 = ((-1, 0), (0, -1))
    tst_dat3 = instance.dbg_data_rotate(instance.dbg_prime_cells, tst_mat2)
    print(instance.dbg_prime_cells)
    print(tst_mat2)
    print(tst_dat3)
    print(tst_dat2 == tst_dat3)
    # ·
    print()
    for rot_angle in [0, 90, 180, 270]:
        rot = _angle_to_matrix(rot_angle)
        print(rot_angle, rot, instance.dbg_is_rot_mat(rot))
    rot = _angle_to_matrix(45)
    print(45, rot, instance.dbg_is_rot_mat(rot))
    print("reflection", HORIZONTAL_REFLECTION, instance.dbg_is_rot_mat(HORIZONTAL_REFLECTION))
    print()
    tst_mat3 = instance.dbg_mat_rotate(tst_mat1, tst_mat1)
    print(tst_mat3)
def _hold_comments() -> None:
    '''comment block

    rotate (x,y) by 90° = ((0, -1), (1, 0)) · (x, y)
    rotate (x,y) by 180° = ((0, -1), (1, 0)) · [((0, -1), (1, 0)) · (x, y)]
    rotate (x,y) by 180° = ((-1, 0), (0, -1)) · (x, y)

    ((-1, 0), (0, -1)) = (0, -1), (1, 0)) ? (0, -1), (1, 0))

    (0, -1), (1, 0)) ? (0, -1), (1, 0))
    a1       a2        b1      b2
    (0, -1)*(0, -1)  (0, -1)*(1, 0)
    (1, 0)*(0, -1)   (1, 0)*(1, 0)

    rotation matrix
    [ cos θ   - sin θ]
    [ sin θ     cos θ]

    φ = 2 * arccos(cos θ)

    [ cos φ   - sin φ]
    [ sin φ     cos φ]
    '''
    raise NotImplementedError("dummy method")

def my_main() -> None:
    '''wrapper for test/start code so that variables do not look like constants
    '''
    instance = transforms_create_test()
    assert isinstance(instance, AutomataTransforms)
    # _is_rot_mat_test(instance)
    # _rotations_check(instance)
    # _prime_cells_check(instance)
    _check_transform_test(instance)
    # _hashable_transform_test(instance)
    _duplicate_test(instance)
    _collision_test(instance)
    _end_cycle_test(instance)
    _add_transform_test(instance)
    instance.generate_combination_transforms()

    # _matrix_rotate_test(instance)
    # _duplicate_test(instance) # test again after transform(s) added
    # _collision_test(instance) # test again after transform(s) added «also refactoring»
    instance.dbg_report_instance() # DEBUG

# Standalone module execution
if __name__ == "__main__":
    my_main()
