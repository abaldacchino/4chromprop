import json
import numpy as np

tiles = ['AAL', 'ABL', 'AVL', 'ADL', 'BAL', 'BBL', 'BVL', 'BDL', 'VAL', 'VBL', 'VVL',
         'VDL', 'DAL', 'DBL', 'DVL', 'DDL', 'AIBL', 'AIVL', 'BIAL', 'VIAL', 'HL', 'AAdL',
         'ABdL', 'AVdL', 'ADdL', 'BAdL', 'BBdL', 'BVdL', 'BDdL', 'VAdL', 'VBdL', 'VVdL',
         'VDdL', 'DAdL', 'DBdL', 'DVdL', 'DDdL', 'AIBdL', 'AIVdL', 'BIAdL', 'VIAdL', 'HdL']


def getAdjProp(concretePropagations):
    """
    Method finds and returns the adjacency matrix that represents the propagation graph of
    a tile given its concrete propagations, ie a list of permissible transitions from a colouring
    to a colouring using said tile.
    :param concretePropagations: list of propagations in concrete form ie (1,1) ~> (2,3)
    :return: returns an adjacency matrix that represents the graph with vertices being all
    3 colourings of 2 vertices, with edges if a propagation can take you from one colouring
    to the next.
    """
    m = np.zeros((9, 9), dtype=int)
    for edge in concretePropagations:
        # x1, x2 are input colours which lead to output colours y1, y2
        x1 = int(edge[0])
        x2 = int(edge[1])
        y1 = int(edge[3])  # recall we flip end vertices before joining!
        y2 = int(edge[2])

        # positions are as following:
        # (1, 1), (1, 2), (1, 3), (2, 1), (2, 2) ... (3, 3)
        # equations below give the above ordering
        row_position = (3 * (x1 - 1)) + (x2 - 1)
        col_position = (3 * (y1 - 1)) + (y2 - 1)

        # since there is an edge, set value to 1
        m[row_position][col_position] = 1
    return m


def getStringRep(prop):
    """
    Method converts a propagation adjacency matrix to a string which contains all relevant information
    :param prop: aggregated propagation in question, in the form of an adjacency matrix
    :return: returns string representation of propagation to be used in dictionaries, summaries etc.
    """
    # due to symmetry of our array, displaying the first and last columns suffices
    return np.array2string(prop[0]) + " " + np.array2string(prop[1])


def getTilePropagations(filename="info.json"):
    """

    :param filename: file name of tile_info file
    :return: returns a list of all propagations extracted from each tile
             and a dictionary that maps each tile to its propagation matrix
    """
    tile_info = json.load(open(filename))
    tile_props = []  # stores all unique propagations that come out of the 42 tiles
    tile_prop_mapping = {}  # dictionary that maps a tile to its propagation matrix

    for tile in tiles:
        # getting adjacency matrix for the tile
        new_prop = getAdjProp(tile_info[tile]['concrete propagations'])
        # storing in dictionary
        tile_prop_mapping[tile] = new_prop

        # checking if new_prop is already stored in tile_props
        stored = False
        for prop in tile_props:
            if np.equal(prop, new_prop).all():
                stored = True
        # if not, add it to tile_props
        if not stored:
            tile_props += [new_prop]

    return tile_props, tile_prop_mapping


def get_combination(prop1, prop2, all_props):
    # multiplying matrices prop1 and prop2
    combination = np.matmul(prop1, prop2)
    # normalising matrix so all positive integers become 1s
    combination[(combination > 0)] = 1

    # checking if propagation is newly encountered
    stored = False
    for prop in all_props:
        if np.equal(prop, combination).all():
            stored = True

    return combination, stored


if __name__ == "__main__":
    # getting tile propagations
    tile_props, tile_prop_mapping = getTilePropagations()
    # dict of dicts where propagation_combinations[prop1][prop2] = prop1 * prop2
    # ie stores the outcome of multiplying the two adjacency matrices together
    propagation_combinations = {}

    all_props = tile_props.copy()
    new_props = tile_props.copy()
    old_props = []
    # do this until no new propagations are found
    i = 0
    while new_props:
        i += 1
        print("Iteration", i)
        print("Number of old tiles:", len(old_props))
        print("Number of new tiles:", len(new_props))
        next_props = []

        # initialising all new propagations
        for prop in new_props:
            propagation_combinations[getStringRep(prop)] = {}

        # finding pairs of old + new propagations
        for prop1 in old_props:
            for prop2 in new_props:
                # this shouldn't have been calculated already
                if getStringRep(prop2) in propagation_combinations[getStringRep(prop1)]:
                    print("oh no!! 1")

                # order prop1 * prop2
                combination1, stored1 = get_combination(prop1, prop2, all_props)
                propagation_combinations[getStringRep(prop1)][getStringRep(prop2)] = combination1
                # if propagation newly encountered, add it to all_props
                if not stored1:
                    all_props += [combination1]
                    next_props += [combination1]

                # order prop2 *prop1
                combination2, stored2 = get_combination(prop2, prop1, all_props)
                propagation_combinations[getStringRep(prop2)][getStringRep(prop1)] = combination2
                # if propagation newly encountered, add it to all_props
                if not stored2:
                    all_props += [combination2]
                    next_props += [combination2]

        # finding combinations of new * new
        for prop1 in new_props:
            for prop2 in new_props:
                # this shouldn't be already calculated
                if getStringRep(prop2) in propagation_combinations[getStringRep(prop1)]:
                    print("oh no!! 2")

                combination, stored = get_combination(prop1, prop2, all_props)
                # storing the propagation in propagation_combinations
                propagation_combinations[getStringRep(prop1)][getStringRep(prop2)] = combination

                # checking if propagation is newly encountered
                if not stored:
                    all_props += [combination]  # adding new propagation to all_props
                    next_props += [combination]

        # updating old_props to include the now no longer new ones
        old_props = old_props + new_props
        # changing the new_props to be the most recently found new ones
        new_props = next_props
    # It might be interesting to note that this usually ends after 4 loops
    # meaning that every possible propagation is 4 iterations away from the initial tiles
    print("This was done in", i, "iterations")
    # Check for duplicate propagations (just in case)
    count = 0
    for prop1 in all_props:
        for prop2 in all_props:
            if np.equal(prop1, prop2).all():
                count += 1
    print("Number of unique large propagations = ", len(all_props))
    print("Number of duplicates = ", (len(all_props) - count))

    # all_props = tile_props.copy()
    # prev_props = 0
    # new_props = len(all_props)
    # while not prev_props == new_props:
    #     prev_props = new_props
    #     for prop1 in all_props:
    #         for prop2 in all_props:
    #             if getStringRep(prop2) not in propagation_combinations[getStringRep(prop1)]:
    #                 combination = np.matmul(prop1, prop2)
    #                 combination[(combination > 0)] = 1
    #                 propagation_combinations[getStringRep(prop1)][getStringRep(prop2)] = combination
    #                 stored = False
    #                 for prop in all_props:
    #                     if np.equal(prop, combination).all():
    #                         stored = True
    #                 if not stored:
    #                     all_props += [combination]
    #                     propagation_combinations[getStringRep(combination)] = {}
    #     new_props = len(all_props)

    # finding propagations which do not allow a 4-colouring were you to loop them
    naughty_tiles = []
    # finding overlap with the initial tile propagations
    naughty_initials = []
    for prop in all_props:
        if prop[0][0] == 0 and prop[1][1] == 0:
            naughty_tiles += [prop]
            for tile_prop in tile_props:
                if np.equal(tile_prop, prop).all():
                    naughty_initials += [prop]

    print(len(naughty_tiles), "propagations are 4-chromatic when looping back")
    print(len(naughty_initials), "tiles are 4-chromatic when looping back")

    # Note: the following datasets contain redundant info, in different formats
    # some refactoring is in order...

    # finding combinations of 2 propagations which turn out 4-chromatic ie end with a
    # 4-chromatic tile
    naughty_combinations = []
    # dictionary that stores count of how many such combinations end with a given propagation
    end_tile_count = {}
    # dictionary that given a 4-chromatic propagation, points to a list of all the pairs of
    # propagations whose multiplication results in said propagation
    four_chromatic_pairs = {}
    # dictionary that given a propagation that may lead to a 4-chromatic propagation, stores the
    # number of times it will lead to a 4-chromatic propagation
    lead_to_4 = {}

    # initialising dictionaries
    for prop in all_props:
        four_chromatic_pairs[getStringRep(prop)] = []
    for four_colour_tile in naughty_tiles:
        end_tile_count[getStringRep(four_colour_tile)] = 0

    # traversing all propagations
    for prop1 in all_props:
        for prop2 in all_props:
            combination = propagation_combinations[getStringRep(prop1)][getStringRep(prop2)]
            # checking if the combination is 4-chromatic
            naughty_combination = False
            for naughty in naughty_tiles:
                if np.equal(naughty, combination).all():
                    naughty_combination = True
                    end_tile_count[getStringRep(naughty)] += 1
            # if combination if 4-chromatic store data
            if naughty_combination:
                naughty_combinations += [(prop1, prop2)]
                four_chromatic_pairs[getStringRep(combination)] += [(getStringRep(prop1), getStringRep(prop2))]
                # initialising lead_to_4 to 0
                if getStringRep(prop1) not in lead_to_4:
                    lead_to_4[getStringRep(prop1)] = 0
                if getStringRep(prop2) not in lead_to_4:
                    lead_to_4[getStringRep(prop2)] = 0
                lead_to_4[getStringRep(prop1)] += 1
                if not np.equal(prop1, prop2).all(): #checking if propagations are equal
                    lead_to_4[getStringRep(prop2)] += 1

    print("Number of combinations that lead to a 4-chromatic sequence:", len(naughty_combinations))
    print("How are these distributed against 4-chromatic propagations")
    print(end_tile_count)

    i = 1
    for naughty in naughty_tiles:
        print("Propagation", i)
        print("Represented by: ", getStringRep(naughty))
        print("This one is led to by", len(four_chromatic_pairs[getStringRep(naughty)]), "tile combinations")
        print(four_chromatic_pairs[getStringRep(naughty)])
        i += 1

    print("Which propagations can result in 4-chromatic combinations?")
    print("How often do they do so?")
    print(lead_to_4)
