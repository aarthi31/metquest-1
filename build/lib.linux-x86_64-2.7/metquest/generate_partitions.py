# -*- coding: utf-8 -*-

from __future__ import absolute_import
from itertools import product


def generate_partitions(maximumvalue, lbnumlist, columnvalue):
    """
    This code takes as input the columnvalue (j), values of the shortest path
    of each of the metabolites (given as a list) and the sum that has to be
    obtained using these combination of numbers.

    Parameters
    ----------
    maximumvalue : int
        Maximum values which the numbers can take
    lbnumlist : list
        a list of values pertaining to the length of shortest paths of
        every metabolite
    columnvalue : int
        Desired sum to be obtained
        All the partitions of numbers which will generate the desired sum
        whose values are between the values for shortest paths and the
        maximum values.
    Returns
    -------
    all_partitions : List of tuples

    Notes
    -----
    For instance, if the column value is 7, the number of imputs is 2,
    and the shortest path of the metabolites is 4,3 respectively, and the
    maximum sum that has to be obtained is 8, then

    >>> generate_partitions(7,[4,3],8)
    [(4, 4), (5, 3)]

    >>> generate_partitions(4, [2,1,1], 5)
    [(2, 1, 2), (2, 2, 1), (3, 1, 1)]

    """
    all_combinations = []
    for entry in lbnumlist:
        all_combinations.append(list(range(entry, columnvalue + 1)))
    all_partitions = []
    for partitions in product(*all_combinations):
        if sum(partitions) == maximumvalue:
            all_partitions.append(partitions)
    return all_partitions
