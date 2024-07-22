import json
import numpy as np
import dfa

tiles = ['AAL', 'ABL', 'AVL', 'ADL', 'BAL', 'BBL', 'BVL', 'BDL', 'VAL', 'VBL', 'VVL',
         'VDL', 'DAL', 'DBL', 'DVL', 'DDL', 'AIBL', 'AIVL', 'BIAL', 'VIAL', 'HL', 'AAdL',
         'ABdL', 'AVdL', 'ADdL', 'BAdL', 'BBdL', 'BVdL', 'BDdL', 'VAdL', 'VBdL', 'VVdL',
         'VDdL', 'DAdL', 'DBdL', 'DVdL', 'DDdL', 'AIBdL', 'AIVdL', 'BIAdL', 'VIAdL', 'HdL']


def get_trans(triples, state, input_symbol):
    if state == -1:
        return input_symbol

    match = -1
    # trying to find if (state, input_symbol, y) is a triple in triples
    for t in triples:
        if triples[0] == state and triples[1] == input_symbol:
            match = triples[2]

    if match == -1:  # case no match found
        return -2
    else:
        return match


def getAdjProp(concretePropagations):
    """
    Method finds and returns the adjacency matrix that represents the propagation graph of
    a tile given its concrete propagations, ie a list of permissible transitions from a colouring
    to a colouring using said tile. We join (a, b) to (c,d) by identifying a with d and b with c.

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
            if props_equal(prop, new_prop):
                stored = True
        # if not, add it to tile_props
        if not stored:
            tile_props += [new_prop]

    return tile_props, tile_prop_mapping


def get_combination(prop1, prop2, all_props):
    """
    Method multiplies prop1 * prop2 to get another propagation and checks if
    the propagation is newly encountered by checking with all_props
    :param prop1: Adjacency matrix to be multiplied first
    :param prop2: Adjacency matrix to be multiplied second
    :param all_props: List of all propagations
    :return: Returns the resulting propagation and whether it is stored in the all_props list or not
    """
    # multiplying matrices prop1 and prop2
    combination = np.matmul(prop1, prop2)
    # normalising matrix so all positive integers become 1s
    combination[(combination > 0)] = 1

    # checking if propagation is newly encountered
    stored = prop_in_list(combination, all_props)

    return combination, stored


def props_equal(prop1, prop2):
    """
    Checks if prop1 is equal to prop2
    :param prop1: Adjacency matrix for a propagation
    :param prop2: Adjacency matrix for a propagation
    :return: True/False depending on whether propagations are equal or not
    """
    return np.equal(prop1, prop2).all()


def prop_in_list(prop, proplist):
    """
    Method checks if a propagation is in a list of propagations
    :param prop: propagation to be checked
    :param proplist: list of propagations
    :return: True if found, False if not found in list
    """
    stored = False
    for prop2 in proplist:
        if props_equal(prop, prop2):
            stored = True
    return stored


def remove_duplicates(prop_list):
    """
    Removes duplicates from a propagation list
    :param prop_list: list of propagations
    :return: filtered list of propagations such that no two are equal
    """
    unique_list = []
    for prop in prop_list:
        if not prop_in_list(prop, unique_list):
            unique_list += [prop]

    return unique_list


def all_successors(prop, all_props, propagation_combinations):
    """
    Finds all propagations that arise from propagation prop
    :param prop:
    :param all_props:
    :param propagation_combinations:
    :return:
    """
    successors = []
    for prop2 in all_props:
        # finding propagations that can arise from prop
        successor = propagation_combinations[getStringRep(prop)][getStringRep(prop2)]
        if not prop_in_list(successor, successors):
            successors += [successor]

    new_successors = successors.copy()
    while new_successors:
        next_successors = []
        for succ in new_successors:
            for prop2 in all_props:
                successor = propagation_combinations[getStringRep(succ)][getStringRep(prop2)]
                if not prop_in_list(successor, successors):
                    successors += [successor]
                    next_successors += [successor]
        new_successors = next_successors
    return successors


def find_and_print_succ_data(all_props, tile_props, tile_prop_mapping, naughty_propagations, propagation_combinations):
    """
    Method finds and prints data regarding successors and 4-chromatic tiles
    :param all_props: All propagations
    :param tile_props: Propagations for the tiles
    :param tile_prop_mapping: Mapping of a tile label to its propagation
    :param naughty_propagations: 4-chromatic propagations
    :param propagation_combinations: dictionary which specifies all combinations of propagations
    :return:
    """
    tiles_can_get_to_4 = []
    for prop in tile_props:
        successors = all_successors(prop, tile_props, propagation_combinations)
        get_to_4 = False
        # number of propagations don't change even if this is ommitted
        if prop_in_list(prop, naughty_propagations):
            get_to_4 = True
        # this means every 4-chromatic tile also has 4-chromatic predecessors
        for successor in successors:
            if prop_in_list(successor, naughty_propagations):
                get_to_4 = True
        if get_to_4:
            tiles_can_get_to_4 += [prop]
    print("No of tile propagations that feature in a 4 chromatic sequence = ", len(tiles_can_get_to_4))

    # finding the corresponding tiles to above propagations
    could_be_4chrom_tiles = []
    for tile in tiles:
        for tile_prop in tiles_can_get_to_4:
            if props_equal(tile_prop_mapping[tile], tile_prop):
                could_be_4chrom_tiles += [tile]
                print(tile, getStringRep(tile_prop))
    print("Tiles that can be included in a 4 chromatic sequence =", could_be_4chrom_tiles)
    print("Number of tiles of interest =", len(could_be_4chrom_tiles))

    # finding all propagations that can be included in a 4-chromatic sequence
    props_can_get_to_4 = []
    for prop in all_props:
        successors = all_successors(prop, all_props, propagation_combinations)
        get_to_4 = False
        if prop_in_list(prop, naughty_propagations):
            get_to_4 = True
        for successor in successors:
            if prop_in_list(successor, naughty_propagations):
                get_to_4 = True
        if get_to_4:
            props_can_get_to_4 += [prop]
    print("No of propagations that feature in a 4 chromatic sequence = ", len(props_can_get_to_4))

    naughty_can_get_to_4 = []
    stopping4 = []
    for prop in naughty_propagations:
        successors = all_successors(prop, all_props, propagation_combinations)
        get_to_4 = False
        for successor in successors:
            if prop_in_list(successor, naughty_propagations):
                get_to_4 = True
        if get_to_4:
            naughty_can_get_to_4 += [prop]
        if not get_to_4:
            stopping4 += [prop]
    print("From the ", len(naughty_propagations), " propagations that are 4-chromatic")
    print(len(naughty_can_get_to_4), "have 4-chromatic successors")

    # finding which propagations can end in a 4-colourable propagation that can still be 4
    # colourable after we add more tiles!
    tiles_can_get_to_4 = []
    for prop in tile_props:
        successors = all_successors(prop, tile_props, propagation_combinations)
        get_to_4 = False
        for successor in successors:
            if prop_in_list(successor, naughty_can_get_to_4):
                get_to_4 = True
        if get_to_4:
            tiles_can_get_to_4 += [prop]
    print("No of tile propagations that feature in a 4 chromatic sequence  that do not end= ", len(tiles_can_get_to_4))


def find_and_print_4chrom_data(all_props, naughty_propagations, propagation_combinations):
    """
    Messy method that will display nicely some statistics on 4-chromatic propagations
    :param all_props:
    :param naughty_propagations:
    :param propagation_combinations:
    :return:
    """

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
    for four_colour_tile in naughty_propagations:
        end_tile_count[getStringRep(four_colour_tile)] = 0

    # traversing all propagations to populate data
    for prop1 in all_props:
        for prop2 in all_props:
            combination = propagation_combinations[getStringRep(prop1)][getStringRep(prop2)]
            # checking if the combination is 4-chromatic
            naughty_combination = False

            for naughty in naughty_propagations:
                if props_equal(naughty, combination):
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
                if not props_equal(prop1, prop2):  # checking if propagations are equal
                    lead_to_4[getStringRep(prop2)] += 1

    print("Number of combinations that lead to a 4-chromatic sequence:", len(naughty_combinations))
    print("How are these distributed against 4-chromatic propagations")
    print(end_tile_count)

    i = 1
    countpropcombs = 0
    for naughty in naughty_propagations:
        print("Propagation", i)
        print("Represented by: ", getStringRep(naughty))
        print("This one is led to by", len(four_chromatic_pairs[getStringRep(naughty)]), "tile combinations")
        countpropcombs += len(four_chromatic_pairs[getStringRep(naughty)])
        print(four_chromatic_pairs[getStringRep(naughty)])
        i += 1
    print("countpropcombs", countpropcombs)

    print("Which propagations can result in 4-chromatic combinations?")
    print("How often do they do so?")
    print("No of propagations that can lead to a 4 chromatic sequence =", len(lead_to_4.keys()))
    print(lead_to_4)


def get_index(prop, prop_list):
    for i, p in enumerate(prop_list):
        if props_equal(prop, p):
            return i
    return -1


def get_combination_subset(all_props, naughty_props, propagation_combinations, edge_labels):
    prop_list = naughty_props.copy()
    edge_list = []

    for prop1 in all_props:
        for prop2 in all_props:
            combination = propagation_combinations[getStringRep(prop1)][(getStringRep(prop2))]
            if prop_in_list(combination, naughty_props):
                # adding components if they are not already in the prop list
                for prop in [prop1, prop2]:
                    if not prop_in_list(prop, prop_list):
                        prop_list += [prop]

    filt_edge_label = []
    for edge_label in edge_labels:
        if prop_in_list(edge_label, prop_list):
            filt_edge_label += [edge_label]

            # vertices will just be numbers 0, 1, 2 ... m for simplicity's sake
    vertices = list(range(0, len(prop_list)))
    vertex_prop_map = {}
    edge_indices = []
    final_states = []
    for i, prop in enumerate(prop_list):
        vertex_prop_map[i] = prop
        if prop_in_list(prop, filt_edge_label):
            edge_indices += [i]
        if prop_in_list(prop, naughty_props):
            final_states += [i]

    for prop1 in prop_list:
        for prop2 in filt_edge_label:
            combination = propagation_combinations[getStringRep(prop1)][(getStringRep(prop2))]
            if prop_in_list(combination, prop_list):
                edge_list += [(get_index(prop1, prop_list), get_index(prop2, prop_list),
                               get_index(combination, prop_list))]

    return vertices, edge_list, edge_indices, vertex_prop_map, final_states


if __name__ == "__main__":
    # getting tile propagations
    tile_props, tile_prop_mapping = getTilePropagations()
    # dict of dicts where propagation_combinations[prop1][prop2] = prop1 * prop2
    # ie stores the outcome of multiplying the two adjacency matrices together

    count =0
    for tile_prop in tile_props:

        if tile_prop[1][1] == 0 and tile_prop[0][3] == 0:
            count += 1
            print(next(key for key, value in tile_prop_mapping.items() if props_equal(value, tile_prop)))
    print(count)
    count = 0
    for tile_prop in tile_props:
        if tile_prop[1][1] == 0 and tile_prop[0][0] ==0:
            count += 1
            print(next(key for key, value in tile_prop_mapping.items() if props_equal(value, tile_prop)))
    print(count)

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
    # print("This was done in", i, "iterations")
    # Check for duplicate propagations (just in case)
    count = 0
    for prop1 in all_props:
        for prop2 in all_props:
            if props_equal(prop1, prop2):
                count += 1
    print("Number of unique large propagations = ", len(all_props))
    print("Number of duplicates = ", (len(all_props) - count))

    # finding propagations which do not allow a 4-colouring were you to loop them
    naughty_propagations = []
    # finding overlap with the initial tile propagations
    naughty_tile_props = []
    for prop in all_props:
        # if prop[0][0] == 1 and prop[1][1] == 0:
            # print(prop)
            # if prop_in_list(prop, tile_props):
            #     print("is a tile propagation")
        if prop[0][0] == 0 and prop[1][1] == 0:
            naughty_propagations += [prop]
            if prop_in_list(prop, tile_props):
                naughty_tile_props += [prop]

    print(len(naughty_propagations), "propagations are 4-chromatic when looping back")
    print(len(naughty_tile_props), "tiles are 4-chromatic when looping back")

    vertices, edge_list, edge_indices, vertex_prop_map, final_states = (get_combination_subset(all_props,
                                                                                               naughty_propagations,
                                                                                               propagation_combinations,
                                                                                               tile_props))
    tile_vertex_map = {}
    for edge_index in edge_indices:
        prop = vertex_prop_map[edge_index]
        for tile in tile_prop_mapping.keys():
            if props_equal(tile_prop_mapping[tile], prop):
                tile_vertex_map[tile] = edge_index
    print(vertices, edge_list, edge_indices, vertex_prop_map, final_states, tile_vertex_map)

    automata = dfa.DFA(
        start=-1,
        inputs=tile_vertex_map.keys(),
        label=lambda s: s in final_states,
        transition=lambda s, c: get_trans(edge_list, s, tile_vertex_map[c]),
    )

    print(dfa.dfa2dict(automata))
    print("********************************")

    automata = automata.minimize()
    print(dfa.dfa2dict(automata))
    #
    # for i in range(50):
    #     print(i, vertex_prop_map[i][0])
    #     if i == 22:
    #         print(vertex_prop_map[22])
