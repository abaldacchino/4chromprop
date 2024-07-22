
import json
import copy
import numpy as np

import findpropagation

insides = ["AA", "AB", "AV", "AD",
           "BA", "BB", "BV", "BD",
           "VA", "VB", "VV", "VD",
           "DA", "DB", "DV", "DD",
           "AIB", "AIV", "BIA", "VIA",
           "H"]

frames = ["L", "dL"]

def containsAll(tilings):
    for f in frames:
        for i in insides:
            tile = i + f

            if tile not in tilings:
                return tile
    return "ok"


def isUndirected(tilings):
    for f in frames:
        for i in insides:
            tile = i + f
            tileDict = tilings[tile]
            for node in tileDict.keys():
                adjacencies = tileDict[node]
                for adjacent in adjacencies:
                    adj_adjacencies = tileDict[str(adjacent)]
                    if int(node) not in adj_adjacencies:
                        return tile + " " + node + " " + str(adjacent)

    return "ok"


def createMatrices(tileDict):
    # print(tileDict)
    vertices = list(tileDict.keys())
    dim = len(vertices)
    m = np.zeros((dim, dim), dtype=int)
    for i in range(dim):
        vertex = vertices[i]
        adjacencies = tileDict[vertex]
        for adjacent in adjacencies:
            j = vertices.index(str(adjacent))
            m[i][j] = 1
    return m


def getMathematicaAdj(m):
    s = "{"
    for j in range(len(m)):
        l = m[j]
        if j !=0:
            s += ", "
        s += "{"
        s += str(l[0])
        for i in range(1, len(l)):
            s += ", " + str(l[i])
        s += "}"

    s += "}"
    return s

def getMathematicaAdjRemoving4(m):
    s = "{"
    for j in range(len(m)):
        if j == 3:
            continue
        l = m[j]
        if j !=0:
            s += ", "
        s += "{"
        s += str(l[0])
        for i in range(1, len(l)):
            if i == 3:
                continue
            s += ", " + str(l[i])
        s += "}"

    s += "}"
    return s


def getMathematicaAdjRemoving4and3(m):
    s = "{"
    for j in range(len(m)):
        if j == 3 or j==2:
            continue

        l = m[j]
        if j !=0:
            s += ", "
        s += "{"
        s += str(l[0])
        for i in range(1, len(l)):
            if i == 3 or i==2:
                continue
            s += ", " + str(l[i])
        s += "}"

    s += "}"
    return s

if __name__ == "__main__":
    with open('tilings_2.json') as json_file:
        tilings = json.load(json_file)
        print(containsAll(tilings))
        print(isUndirected(tilings))

        tile_data = {}
        for f in frames:
            for i in insides:
                tile = i + f

                tile_data[tile] = {}
                tile_data[tile]["adjacency matrix"] = getMathematicaAdj(createMatrices(tilings[tile]))
                tile_data[tile]["dict adjacencies"] = tilings[tile]
                tile_data[tile]["aggregated propagations"] = findpropagation.getPropagations(tile)
                tile_data[tile]["simple propagations"] = findpropagation.getSimplePropagations(tile)
                tile_data[tile]["concrete propagations"] = findpropagation.getConcretePropagations(tile)

        with open("info.json", "w") as outfile:
            json.dump(tile_data, outfile, indent=4)
