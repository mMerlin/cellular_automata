#!/usr/bin/env python
# coding=utf-8

'''
manipulation of cellular automata generations
'''

# pipenv shell

# standard library imports
import math
# local application/library specific imports
import automata_typehints as AHint

class AutomataUniverse:
    '''properties for a cellular automata universe

    Instances are Immutable

    :property dimensions: the number of dimension for the universe
    :type dimensions: int
    :property neighbourhood: all cell addresses that are neighbors of the universe origin
    :type neighbourhood: frozenset of coordinate tuples for the origin cell address
    :property neighbourhood_population: the number of cell address in a neighbourhood
    :type neighbourhood_population: int
    :property survival_rules: neighbour counts required for living cell to survive to the next
        generation
    :type survival_rules: frozen set of integers
    :property birth_rules: neighbour counts required for an empty cell to spawn into the next
        generation
    :type birth_rules: frozen set of integers

    :property rotate_reflect: matrices to generate equivalent cell patterns
    :type: tuple of tuples «of tuples»
    '''

    def __init__(self,
            neighbourhood: AHint.NeighbourhoodInputType,
            survival_counts: AHint.PropagationRuleInputType,
            birth_counts: AHint.PropagationRuleInputType):
        '''constructor'''
        self._origin_neighbourhood = frozenset(neighbourhood)
        self._survive = frozenset(survival_counts)
        self._birth = frozenset(birth_counts)
        self._dimensions = None
        # self._rotate_reflect = None
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
    def neighbourhood(self) -> AHint.CellGroupSnapshotType:
        return self._origin_neighbourhood

    @property
    def neighbourhood_population(self) -> int:
        '''get the size (number of cells) in the neighbourhood

        :returns population: number of possible neighbours for a single cell
        :rtype: int
        '''
        return len(self._origin_neighbourhood)

    @property
    def survival_rules(self) -> AHint.PropagationRuleType:
        return self._survive

    @property
    def birth_rules(self) -> AHint.PropagationRuleType:
        return self._birth

    # end of property methods

    def __hash__(self) -> int:
        # hash of the automata universe configuration
        # return hash((self.survival_rules, self.birth_rules, self.neighbourhood,
        #     self._rotate_reflect))
        return hash((self.survival_rules, self.birth_rules, self.neighbourhood))

    def is_rotation_matrix(self, matrix: AHint.TransformInputType):
        '''check if the matrix is a valid rotation (about the origin) matrix for the universe

        https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions

        :param matrix: possible rotation matrix
        :type matrix: Iterable of «dimension» cell address tuples
        '''
        self.validate_matrix(matrix)
        if self.dimensions != 2:
            raise NotImplementedError("rotation matrix check only implemented for 2 dimensions")
        # print("is_rotation_matrix?", len(matrix), len(matrix[0]), self.dimensions) # DEBUG
        try:
            rotation_radians = [math.acos(matrix[0][0]), math.asin(-matrix[0][1]),
                math.asin(matrix[1][0]), math.acos(matrix[1][1])]
        except ValueError as ex:
            if ex.args != ("math domain error",):
                raise
            return False
        rotation_degrees = tuple(math.degrees(rad_angle) for rad_angle in rotation_radians)
        # print(rotation_radians) # DEBUG
        # print(rotation_degrees,
        #     rotation_degrees[0] == rotation_degrees[3],
        #     rotation_degrees[1] == rotation_degrees[2],
        #     rotation_degrees[0] == rotation_degrees[1],
        #     max(rotation_degrees[0], rotation_degrees[1]) - 180 ==
        #     min(rotation_degrees[0], rotation_degrees[1])
        # ) # DEBUG
        return rotation_degrees[0] == rotation_degrees[3] and \
            rotation_degrees[1] == rotation_degrees[2] and \
            (rotation_degrees[0] == rotation_degrees[1] or
            max(rotation_degrees[0], rotation_degrees[1]) - 180 ==
            min(rotation_degrees[0], rotation_degrees[1]))
    # end def is_rotation_matrix()

    def step(self, cells: set[AHint.CellAddressType]) -> AHint.CellGroupWorkingType:
        '''iterate from the current generation to the next

        Any cell with a number of neighbors that is in the survial rule continues to the next
        generation

        Any neighbor location of a cell that does not hold a cell is a candidate to spawn a cell
        in«to» the next generation

        Any candidate location with a number of neighbor cells that is in the birth rule becomes
        a new living cell in the next generation

        :param cells: living cells
        :type cells: set of universe cell address tuples
        :returns: next generation of cells for universe configuration
        :rtype: set of universe cell address tuples
        '''
        self._check_cell_group(cells)
        new_generation = set() # start with an empty next generation cell set
        womb_candidates = set() # no initial candidates for new cells either
        for living_cell in cells:
            cell_neighbourhood = self.neighbours(living_cell)
            # print(living_cell, "neighbourhood", cell_neighbourhood) # DEBUG
            cell_neighbors = cell_neighbourhood.intersection(cells)
            neighbour_count = len(cell_neighbors)
            # print(living_cell, neighbour_count, "neighbors", cell_neighbors) # DEBUG
            print(living_cell, neighbour_count, "neighbors") # DEBUG
            if neighbour_count in self.neighbourhood:
                new_generation.add(living_cell)
            empty_cells = cell_neighbourhood.difference(cell_neighbors)
            # print(living_cell, len(empty_cells), "no neighbours", empty_cells) # DEBUG
            assert neighbour_count + len(empty_cells) == self.neighbourhood_population, \
                "living and dead neighbours should add up to neighbourhood size"
            # print(living_cell, neighbour_count, "neighbors,",
            #     len(empty_cells), "no neighbours") # DEBUG
            womb_candidates.update(empty_cells)
        # print(len(womb_candidates), "womb candidates", womb_candidates) # DEBUG
        print(len(womb_candidates), "womb candidates") # DEBUG
        for womb_cell in womb_candidates:
            womb_neighbourhood = self.neighbours(womb_cell)
            # print(womb_cell, "womb neighbourhood", womb_neighbourhood) # DEBUG
            womb_neighbors = womb_neighbourhood.intersection(cells)
            parent_count = len(womb_neighbors)
            # print(womb_cell, parent_count, "parents", womb_neighbors) # DEBUG
            if parent_count in self.birth_rules:
                new_generation.add(womb_cell)
        return new_generation
    # end def step(self)

    def neighbours(self, address: AHint.CellAddressType) -> \
            AHint.CellGroupSnapshotType:
        '''the set of neighbours for a cell address

        :param address: coordinates for cell in n-space
        :type address: tuple of integers
        :returns neighbours:
        :rtype neighbours: frozenset of coordinate tuples with same dimensionality as address
        '''
        self.validate_address(address)
        return frozenset(tuple(sum(ele) for ele in zip(address, neighbour))
            for neighbour in self._origin_neighbourhood)
        # neighbourhood = []
        # for neighbour in self._origin_neighbourhood:
        #     neighbourhood.append(tuple(sum(ele) for ele in \
        #         zip(address, neighbour)))
        # return frozenset(neighbourhood)
    # end def get_neighbours()

    def cell_group_translate(self, cells: AHint.CellGroupType,
            offset_vector: AHint.CellAddressType) -> AHint.CellGroupWorkingType:
        '''add offset «vector» to each cell coordinate

        :param cells: group of cells to translate (move)
        :type cells: «frozen»set of cell address tuples
        :param offset_vector: coordinate delta values
        :type offset_vector: tuple of integers
        :returns translated (moved) cell coordinates
        :rtype: set of tuples of «dimension» integers
        '''
        self._check_cell_group(cells)
        self.validate_address(offset_vector)
        return [tuple(base + delta for base, delta in zip(cell, offset_vector)) for cell in cells]
    # end def cell_group_translate()

    def cell_group_transform(self, cells: AHint.CellGroupSnapshotType,
            matrix: AHint.TransformInputType) -> AHint.CellGroupSnapshotType:
        '''perform rotation or reflection on a set of cells around the origin

        :param cells: group of cells to rotate or reflect around¦across the origin
        :type cells: «frozen»set of cell address tuples
        :param matrix: rotation or reflection matrix
        :type matrix: iterable of cell address vector tuples
        :returns rotated or reflected set of cell addresses
        :rtype frozenset of cell address tuples
        '''
        self.validate_matrix(matrix)
        self._check_cell_group(cells)
        return frozenset(self._vector_dot_product(cell, matrix) for cell in cells)
        # return frozenset(self.matrix_transform(cells, matrix))
    # end def cell_group_rotate_reflect()

    def matrix_transform(self, matrix: AHint.TransformInputType,
            transform: AHint.TransformInputType) -> AHint.TransformType:
        '''transform a cell address matrix through a rotation or reflection matrix

        :param matrix: rotation matrix
        :type matrix: tuple of «dimension» cell address coordinates
        :param transform: rotation matrix
        :type transform: tuple of «dimension» cell address coordinates
        :returns: rotated or reflected matrix
        :rtype: tuple of cell address tuples
        '''
        self.validate_matrix(transform)
        return tuple(self._vector_dot_product(row, transform) for row in matrix)
    # end def matrix_transform()

    def vector_dot_product(self, vector: AHint.CellAddressType, matrix: AHint.TransformInputType) \
            -> AHint.CellAddressType:
        '''transform a cell address vector through a rotation or reflection matrix

        :param vector: cell address
        :type vector: tuple of integer coordinates for universe dimensions
        :param matrix: rotation matrix
        :type matrix: tuple of «dimension» cell address coordinates
        :returns: rotated or reflected cell address
        :rtype: universe cell address tuple
        '''
        self.validate_matrix(matrix)
        return self._vector_dot_product(vector, matrix)
    # end def vector_dot_product()

    def _vector_dot_product(self, vector: AHint.CellAddressType, matrix: AHint.TransformInputType) \
            -> AHint.CellAddressType:
        '''transform a cell address vector through a rotation or reflection matrix

        matrix is not validated here. Expected to be validated once by caller, then
        used multiple times with different vectors

        :param vector: cell address
        :type vector: tuple of integer coordinates for universe dimensions
        :param matrix: rotation matrix
        :type matrix: tuple of «dimension» cell address coordinates
        '''
        self.validate_address(vector)
        # Does this still work with a one dimensional universe?
        return tuple(sum(x * y for x, y in zip(row, vector)) for row in matrix)
        # product = []
        # for ele in matrix:
        #     product.append(x*y for x,y in zip(ele, vector))
        # result = [sum(x * y for x, y in zip(row, vector)) for row in matrix]
        # return tuple(result)
        # return tuple([sum(x * y for x, y in zip(row, vector)) for row in matrix])
    # end def _vector_dot_product()

    def is_universe_address(self, address: AHint.CellAddressType) -> bool:
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

    def validate_address(self, address: AHint.CellAddressType):
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

    def validate_matrix(self, matrix: AHint.TransformInputType):
        '''verify matrix is valid to perform dot product with universe cell address

        :param matrix: (ordered) series of «dimension» cell address vectors
        :type matrix: any iterable of «dimension» cell address tuples
        :raises: TypeError, ValueError
        '''
        for ele in matrix:
            self.validate_address(ele)
        if len(matrix) != self.dimensions:
            raise TypeError("matrix contains {} vectors, which does not match {} dimensions "
                "for the configured universe".format(len(matrix), self.dimensions))
        _test_subscriptable = matrix[0]
    # end def validate_matrix()

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
            neighbours = self.neighbours(addr)
            if not origin in neighbours:
                raise ValueError((addr, "no symmetric address in neighbourhood"))
    # end def _validate_universe_data()

    def _check_cell_group(self, cells: AHint.CellGroupType):
        '''check that cells are valid universe address tuples

        :param cells: universe cell addresses
        :type cells: set of universe cell address tuples
        :raises: TypeError
        '''
        if not isinstance(cells, (set, frozenset)):
            raise TypeError((type(cells), "automata cell group is not a set"))
        for cell in cells:
            self.validate_address(cell)
    # end def _check_cell_group()

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
        if not self.neighbourhood_population == input_size:
            raise ValueError((input_size, self.neighbourhood_population,
                "Neighbourhood addresses are not unique"))
        if self.neighbourhood_population < 2:
            raise ValueError((self.neighbourhood_population, "neighbourhood does not contain"
                " enough neighbours: at least 2 required"))
            # Two neighbours is the absolute minimum to maintain symmetry. It would be valid
            # for a simple one dimensional universe
        # use first element to get the number of dimension to match for the rest
        neighbourhood_iter = iter(self._origin_neighbourhood)
        address = next(neighbourhood_iter)
        if not isinstance(address, tuple):
            raise TypeError((type(address), "element of neighbourhood is not a tuple"))
        self._dimensions = len(address)
        if self.dimensions < 1: # support for 1 dimensional automata?? or only 2+?
            raise TypeError("neighbourhood address does not contain at least 1 coordinate")
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

    def _check_propagation_type(self, counts_rule: AHint.PropagationRuleType,
            input_size: int):
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
            if not 0 <= count <= self.neighbourhood_population:
                raise ValueError((count, self.neighbourhood_population,
                    "cell propagation count case is not between 0 and neighbourhood size"))
    # end def _check_propagation_type()
# end class AutomataUniverse


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
ROTATE90 = [(0, -1), (1, 0)]
TEST_CELLS = frozenset(((7, 23),(41, 5)))
ROTATED_CELLS = frozenset(((-23, 7),(-5, 41)))

def vector_test(universe: AutomataUniverse):
    '''exercise vector dot product code'''
    # mat = None
    # mat = []
    # mat = [1, 2]
    # mat = [(1,), 2]
    mat = ROTATE90
    # in_xy = None
    # in_xy = (1,)
    in_xy = (7, 23)
    out_xy = universe.vector_dot_product(in_xy, mat)
    print("{} translated {}".format(in_xy, out_xy)) # DEBUG

def group_test(universe: AutomataUniverse):
    '''exercise cell group rotate¦reflect code'''
    # mat = None
    mat = ROTATE90
    # in_cells = None
    # in_cells = []
    # in_cells = frozenset() # good
    # in_cells = set() # good
    in_cells = set()
    # in_cells.add(None)
    # in_cells.add(tuple())
    # in_cells.add((1.2, 2))
    # in_cells.add((7, 23)) # good
    # in_cells = set([(7, 23)]) # good
    # in_cells = set(((7, 23),)) # good
    # in_cells = frozenset(((7, 23),)) # good
    # in_cells = set(((7, 23),(41, 5))) # good
    # in_cells = set(TEST_CELLS) # good
    in_cells = TEST_CELLS # good
    out_cells = universe.cell_group_transform(in_cells, mat)
    print("{} transforms to {}".format(in_cells, out_cells)) # DEBUG
    assert out_cells == ROTATED_CELLS

def my_main():
    '''wrapper for test/start code so that variables do not look like constants'''
    # code to minimally exercise the class methods
    uni = AutomataUniverse(NEIGHBOURS_LIST, [2,3], [3])
    print("universe dimensions", uni.dimensions)
    print("neighbourhood", uni.neighbourhood)
    print("neighbourhood size", uni.neighbourhood_population)
    print("hash", hash(uni))
    print("(0,0) neighbours", uni.neighbours((0,0)))
    print()
    vector_test(uni)
    group_test(uni)

# Standalone module execution
if __name__ == "__main__":
    my_main()
