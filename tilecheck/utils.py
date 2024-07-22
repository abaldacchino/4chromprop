
import numpy as np


def getStringRep(prop):
    """
    Method converts a propagation adjacency matrix to a string which contains all relevant information
    :param prop: aggregated propagation in question, in the form of an adjacency matrix
    :return: returns string representation of propagation to be used in dictionaries, summaries etc.
    """
    # due to symmetry of our array, displaying the first and second columns suffices
    return np.array2string(prop[0]) + " " + np.array2string(prop[1])


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


def get_index(prop, prop_list):
    """
    Method takes a propagation and returns its index in the specified list.
    :param prop: propagation matrix
    :param prop_list: list of propagation matrices
    :return: An integer corresponding to the position of prop in prop_list
    """
    for i, p in enumerate(prop_list):
        if props_equal(prop, p):
            return i
    return -1
