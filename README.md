<!-- cSpell:enable -->
# Cellular Automata

<link href="/home/phil/Documents/data_files/markdown.css" rel="stylesheet"/>

Exploration of concept code to handle multiple¦variable cellular automata dimensions and rules with a single set of functions

* [here](./)
* [top](/home/phil/index.md)
</br></br>

* [Concepts](#link_concepts)
* [References](#link_reference)

<!--
* [Link](#link_link)
## <a name="link_link">⚓</a> Link
-->

## <a name="link_reference">⚓</a> References

* [dot product](https://stackoverflow.com/questions/25592838/generics-of-raw-types-int-float-double-create-weird-error-messages)
* [Multiplying Matrices and Vectors](https://hadrienj.github.io/posts/Deep-Learning-Book-Series-2.2-Multiplying-Matrices-and-Vectors/)
* [multiply two matrices](https://www.geeksforgeeks.org/python-program-multiply-two-matrices/)
* [n dimensional vector dot product](https://stackoverflow.com/questions/28307008/how-to-write-a-function-to-do-the-dot-product-on-n-dimensional-vector)

* [pytest docs](https://docs.pytest.org/en/latest/)
* [pytest getting started](https://docs.pytest.org/en/latest/getting-started.html#get-started)
* [pytest how to assert](https://docs.pytest.org/en/latest/how-to/assert.html)
* [pytest good practices](https://docs.pytest.org/en/latest/explanation/goodpractices.html#test-discovery)

* [pytest-cov plugin reference](https://pytest-cov.readthedocs.io/en/latest/config.html#reference)

* [rotation matrix in 3d](https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions)
* [matrix orthogonal rotation](http://www.euclideanspace.com/maths/algebra/matrix/orthogonal/rotation/)
* [determinant index](https://www.euclideanspace.com/maths/algebra/matrix/functions/determinant/index.htm)
* [determinant of a matrix](https://www.geeksforgeeks.org/determinant-of-a-matrix/) **

* [inheriting from a namedtuple base class](https://stackoverflow.com/questions/42385916/inheriting-from-a-namedtuple-base-class)

## <a name="link_concepts">⚓</a> Concepts

An automaton `cell` is in one of 2 states. Alive or dead. Its only property is its location. Its address in n-dimensional space.

A living cell will either propagate to the next generation or not based only on its neighbours. Specifically the number of living neighbours in the current generation.

A currently dead cell will come to life «be created» in the next generation based only on its neighbours. Specifically the number of living neighbours in the current generation.

'birthing' a cell that has zero living neighbours is not valid. Allowing that would immediately (first generation step) fill the whole universe that is not included in the neighbourhood of the initial populations of cells.

Neighbours are symmetrical. A cell that has a neighbour at «address» is a neighbour of the cell at «address». «address == n-dimensional coordinates»

It is not necessary to search the whole dimensional space to check for dead cells that need to be created in the next generation. It is sufficient to check the dead cells that are neighbours of living cells in the current generation. A dead cell in the current generation that is not a neighbour of a living cell in the current generation can not possibly be born into the next generation.

The extent of the cell neighbourhood determines how much the extent of all living cells can expand in a single generations. That is the maximum distance to a neighbor cell in each dimension is the maximum possible increase in each dimension (in both plus and minus directions). The generation can shrink much faster, potentially going from a very large living population to zero in a single generation.

Using «multiples of» an integer cell to cell distance in each dimension, a tuple of integers can be used to specify any cell address. A tuple is immutable, can be added to a python set. A genertion can then be represented as a single set of the address tuples of the living cells.

The neighbourhood of a cell can be represented by the relative addresses of each of the neighbors of the cell. That is the same as the address of all of the neighbours of the origin cell at address (0, 0, …).

The size «living cell count» of the intersection of the «living» generation set and the neighbors of a cell determine if that cell survives in the next generation.

All of the neighbour cells that are dead (do not exist) in the current generation are candidates to be born in the next generation. The union of the dead neighbours of all living cells in one generation are the only cells that need be checked for creation in the next generation «based on the count of living neighbors of the dead cells».

The union of the surviving and «to be» created cells becomes the next generation «living set».

## functional comment block

* The above header prevents the comments here from being hidden if the previous block is folded in the editor
  * and this prevents the editor from switching to 4 space tabs as well

<!-- cSpell:disable -->
<!-- cSpell:enable -->
<!--
# cSpell:disable
# cSpell:enable
cSpell:words
cSpell:ignore
cSpell:enableCompoundWords
-->
