import json
import copy
import numpy as np
import utils
import dfa

tiles = ['AAL', 'ABL', 'AVL', 'ADL', 'BAL', 'BBL', 'BVL', 'BDL', 'VAL', 'VBL', 'VVL',
         'VDL', 'DAL', 'DBL', 'DVL', 'DDL', 'AIBL', 'AIVL', 'BIAL', 'VIAL', 'HL', 'AAdL',
         'ABdL', 'AVdL', 'ADdL', 'BAdL', 'BBdL', 'BVdL', 'BDdL', 'VAdL', 'VBdL', 'VVdL',
         'VDdL', 'DAdL', 'DBdL', 'DVdL', 'DDdL', 'AIBdL', 'AIVdL', 'BIAdL', 'VIAdL', 'HdL']


def findPropagation(queue, curr_colouring, colours, adjacency_list):
    """
    Method findPropagation attempts to find a colouring of the graph as specified by the
    adjacency list, starting from the colouring specified in currColouring. If a colouring
    is found, the propagation of the graph considering vertices 1, 2 on the left wall and
    3, 4 on the right wall is returned.

    :param queue: vertices to be coloured next
    :param curr_colouring: current (partial) colouring of graph
    :param colours: colours available
    :param adjacency_list: adjacency list of variables
    :return: A list of all propagations with left wall (1, 2) right wall (3, 4)
    """

    if len(queue) == 0:   # If queue empty, no more vertices left to colour (graph is connected)
        return [((curr_colouring[1], curr_colouring[2]), (curr_colouring[3], curr_colouring[4]))]
    else:
        # Colouring vertex v at top of queue
        v = queue[0]
        queue.remove(v)

        free_colours = copy.deepcopy(colours)   # stores free colours left to colour v

        for u in adjacency_list[v]:
            if curr_colouring[u] in free_colours:
                free_colours.remove(curr_colouring[u])

            if curr_colouring[u] == '' and u not in queue:
                queue.append(u)

        if len(free_colours) == 0:   # if no free colours left, there is no colouring with this configuration
            return []

        result = []

        # Attempting to colour graph with all possible colourings for v
        for c in free_colours:
            new_colouring = copy.deepcopy(curr_colouring)
            new_colouring[v] = c

            new_queue = copy.deepcopy(queue)
            result += findPropagation(new_queue, new_colouring, colours, adjacency_list)

        return sorted(list(set(result)), key=lambda t: t[1][0])


def getConcretePropagations(tile):
    """
    Method getConcretePropagations returns all the 'concrete' 3-propagations of the tile specified in
    tileID, that is all possible colourings of the input and output vertices.

    :param tile:
    :return:
    """
    propagations = []

    for _input in [('1', '1'), ('1', '2'), ('1', '3'),
                   ('2', '1'), ('2', '2'), ('2', '3'),
                   ('3', '1'), ('3', '2'), ('3', '3')]:
        if 2 in tile[1] and _input[0] == _input[1]:
            continue

        colouring = dict.fromkeys(tile.keys(), '')
        colouring[1] = _input[0]
        colouring[2] = _input[1]

        queue = []
        for v in tile[1]:
            if v != 2: queue.append(v)

        for p in findPropagation(queue, colouring, ['1', '2', '3'], tile):
            propagations.append((p[0][0], p[0][1], p[1][0], p[1][1]))
    return propagations


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


def getTilePropagations(tile_adjacencies):
    """
    Method getTilePropagations takes a file containing tiles and their structure, and returns
    all unique propagations that come out of the 42 tiles, as well as a mapping from tile to
    propagation.
    :param tile_adjacencies: json file containing adjacency information of each tile
    :return: returns a list of all propagations extracted from each tile
             and a dictionary that maps each tile to its propagation matrix
    """

    tile_props = []  # stores all unique propagations that come out of the 42 tiles
    tile_prop_mapping = {}  # dictionary that maps a tile to its propagation matrix

    for tile in tiles:
        # getting adjacency matrix for the tile
        tile_set = {int(k): v for k, v in tile_adjacencies[tile].items()}

        new_prop = getAdjProp(getConcretePropagations(tile_set))
        # storing in dictionary
        tile_prop_mapping[tile] = new_prop

        # checking if new_prop is already stored in tile_props
        stored = False
        for prop in tile_props:
            if utils.props_equal(prop, new_prop):
                stored = True
        # if not, add it to tile_props
        if not stored:
            tile_props += [new_prop]

    return tile_props, tile_prop_mapping


def getPropagationClosure(starting_props):
    """
    Method getPropagationClosure calculates the closure of matrices in starting_props.

    :param starting_props: Starting set of propagations
    :return: closure_props: a list containing all propagations that arise from normalised matrix multiplication
                            starting from starting_props
             propagation_combinations: a dictionary which encodes the result of the multiplication of any
                            two propagations in all_props
    """
    propagation_combinations = {}

    closure_props = starting_props.copy()
    new_props = starting_props.copy()
    old_props = []
    # Iterate until no new propagations are found
    while new_props:
        next_props = []

        # initialising all new propagations
        for prop in new_props:
            propagation_combinations[utils.getStringRep(prop)] = {}

        # Calculating old * new and new * old
        for prop1 in old_props:
            for prop2 in new_props:

                # Sanity check that new propagation is indeed new
                if utils.getStringRep(prop2) in propagation_combinations[utils.getStringRep(prop1)]:
                    raise Exception(f'Unexpected behaviour! Propagation {prop2} is new, but {prop1}*{prop2} was '
                                    f'already calculated')

                # Order prop1 * prop2
                combination1, stored1 = utils.get_combination(prop1, prop2, closure_props)
                propagation_combinations[utils.getStringRep(prop1)][utils.getStringRep(prop2)] = combination1
                # If propagation newly encountered, add it to all_props
                if not stored1:
                    closure_props += [combination1]
                    next_props += [combination1]

                # order prop2 *prop1
                combination2, stored2 = utils.get_combination(prop2, prop1, closure_props)
                propagation_combinations[utils.getStringRep(prop2)][utils.getStringRep(prop1)] = combination2
                # if propagation newly encountered, add it to all_props
                if not stored2:
                    closure_props += [combination2]
                    next_props += [combination2]

        # Finding combinations of new * new
        for prop1 in new_props:
            for prop2 in new_props:
                # Sanity check for new propagation
                if utils.getStringRep(prop2) in propagation_combinations[utils.getStringRep(prop1)]:
                    raise Exception(f'Unexpected behaviour! Propagations {prop1} and {prop2} are new, but '
                                    f'{prop1}*{prop2} has already been calculated.')

                combination, stored = utils.get_combination(prop1, prop2, closure_props)
                # Storing the propagation in propagation_combinations
                propagation_combinations[utils.getStringRep(prop1)][utils.getStringRep(prop2)] = combination

                # Checking if propagation is newly encountered
                if not stored:
                    closure_props += [combination]  # adding new propagation to all_props
                    next_props += [combination]

        # Updating old_props to include the now calculated propagations
        old_props = old_props + new_props
        # Changing the new_props to be the most recently found new ones
        new_props = next_props

    return closure_props, propagation_combinations


def get_predecessors(propagations, resulting_props, propagation_combinations, edge_labels):
    predecessors = resulting_props.copy()
    edge_list = []

    for prop1 in propagations:
        for prop2 in propagations:
            # Checking if prop1 * prop2 results in resulting_props
            combination = propagation_combinations[utils.getStringRep(prop1)][(utils.getStringRep(prop2))]
            if utils.prop_in_list(combination, resulting_props):
                # adding components if they are not already in the prop list
                for prop in [prop1, prop2]:
                    if not utils.prop_in_list(prop, predecessors):
                        predecessors += [prop]

    # filtered edge label keeps as edges only successor propagations in edge_label
    filt_edge_label = []
    for edge_label in edge_labels:
        if utils.prop_in_list(edge_label, predecessors):
            filt_edge_label += [edge_label]

    # vertices will just be numbers 0, 1, 2 ... m for simplicity's sake
    vertices = list(range(0, len(predecessors)))
    vertex_prop_map = {}    # mapping propagations with vertices
    edge_indices = []
    final_states = []
    for i, prop in enumerate(predecessors):
        vertex_prop_map[i] = prop
        # vertex is an edge if it is in filt_edge_label
        if utils.prop_in_list(prop, filt_edge_label):
            edge_indices += [i]
        # vertex is a final state if it is in resulting_props
        if utils.prop_in_list(prop, resulting_props):
            final_states += [i]

    # defining edge list (v1, l, v2) where v1 -l-> v2 if v1*l = v2
    for prop1 in predecessors:
        for prop2 in filt_edge_label:
            combination = propagation_combinations[utils.getStringRep(prop1)][(utils.getStringRep(prop2))]
            if utils.prop_in_list(combination, predecessors):
                edge_list += [(utils.get_index(prop1, predecessors), utils.get_index(prop2, predecessors),
                               utils.get_index(combination, predecessors))]

    return vertices, edge_list, edge_indices, vertex_prop_map, final_states


def find_corresponding_labels(indices, map_index_to_prop, map_prop_to_label):
    """
    Method find_corresponding_labels takes integer indices, a mapping from indices to propagations and
    a matching for propagations to labels, and returns a mapping from labels to indices, where each label
    is mapped to the index corresponding to its propagation.

    :param indices: list of integer indices
    :param map_index_to_prop: map from integers to propagations
    :param map_prop_to_label: map from propagations to labels
    :return: label_vertex_map: a map from labels to indices (integers)
    """
    label_vertex_map = {}
    for index in indices:
        # finding corresponding propagation for index
        prop = map_index_to_prop[index]
        # storing any label that has propagation prop
        for label in map_prop_to_label.keys():
            if utils.props_equal(map_prop_to_label[label], prop):
                label_vertex_map[label] = index
    return label_vertex_map


def get_trans(triples, state, input_symbol):
    if state == -1:
        return input_symbol

    match = -1
    # trying to find if (state, input_symbol, y) is a triple in triples
    for t in triples:
        if t[0] == state and t[1] == input_symbol:
            match = t[2]

    if match == -1:  # case no match found
        return -2
    else:
        return match


if __name__ == "__main__":
    tile_data = json.load(open("tilings_2.json"))
    # getting tile propagations
    tile_props, tile_prop_mapping = getTilePropagations(tile_data)

    # finding closure of propagations under multiplication
    all_props, mult_result = getPropagationClosure(tile_props)

    # finding all 4 chromatic propagations
    four_chrom_props = []
    for prop in all_props:
        if prop[0][0] == 0 and prop[1][1] == 0:
            four_chrom_props += [prop]

    # calculating vertices, edges and final states of automata, representing propagations as integers
    vertices, edge_list, edge_indices, vertex_prop_map, final_states = \
        get_predecessors(all_props, four_chrom_props, mult_result, tile_props)

    # finding map from tiles to integers representing propagations
    tile_vertex_map = find_corresponding_labels(edge_indices, vertex_prop_map, tile_prop_mapping)

    print("States ", vertices)
    print("Final states", final_states)
    print("Alphabet ", edge_indices)
    print("Transition function", edge_list)
    print("Tiles to numbers", tile_vertex_map)

    automata = dfa.DFA(
        start=-1,
        inputs=tile_vertex_map.keys(),
        label=lambda s: s in final_states,
        transition=lambda s, c: get_trans(edge_list, s, tile_vertex_map[c]),
    )

    automata = automata.minimize()
    dfa_dict = dfa.dfa2dict(automata)
    start_state = dfa_dict[1]
    transition_info = dfa_dict[0]
    # removing sink state
    del transition_info[2]
    for state in transition_info.keys():
        state_transitions = transition_info[state][1]

        for tile in tile_vertex_map.keys():
            if state_transitions[tile] == 2:
                # removing edges to the sink state
                del state_transitions[tile]

    print(start_state)
    print(transition_info)







