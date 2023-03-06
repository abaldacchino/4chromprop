import json
import copy


def findPropagation(queue, currColouring, colours, adjacencyList):
    if len(queue) == 0:
        return [((currColouring[1], currColouring[2]), (currColouring[3], currColouring[4]))]
    else:
        v = queue[0]
        queue.remove(v)

        freeColours = copy.deepcopy(colours)

        for u in adjacencyList[v]:
            if currColouring[u] in freeColours:
                freeColours.remove(currColouring[u])

            if currColouring[u] == '' and u not in queue:
                queue.append(u)

        if len(freeColours) == 0:
            return []

        result = []

        for c in freeColours:
            newColouring = copy.deepcopy(currColouring)
            newColouring[v] = c

            newQueue = copy.deepcopy(queue)
            result += findPropagation(newQueue, newColouring, colours, adjacencyList)

        return sorted(list(set(result)), key=lambda t: t[1][0])


def simplifyPropagations(propagations, colours):
    for _input in set([p[0] for p in propagations]):
        _output = set([p[1] for p in propagations if p[0] == _input])

        cY1 = set(o[0] for o in _output)

        for c in cY1:
            _output_cY1 = set([o[1] for o in _output if o[0] == c])

            if len(_output_cY1) == len(colours):
                for c2 in colours:
                    propagations.remove((_input, (c, c2)))

                propagations.append((_input, (c, '*')))

        _output = set([p[1] for p in propagations if p[0] == _input])
        cY1 = set(o[0] for o in _output)

        for c in cY1:
            _output_cY1 = set([o[1] for o in _output if o[0] == c])

            if len(_output_cY1) == len(colours) - 1:
                for c2 in colours:
                    if c2 not in _output_cY1:
                        propagations.append((_input, (c, '!' + c2)))
                    else:
                        try:
                            propagations.remove((_input, (c, c2)))
                        except Exception:
                            continue

        _output = set([p[1] for p in propagations if p[0] == _input])
        cY2 = set(o[1] for o in _output)

        for c in cY2:
            _output_cY2 = set([o[0] for o in _output if o[1] == c])

            if len(_output_cY2) == len(colours):
                for c2 in colours:
                    propagations.remove((_input, (c2, c)))

                propagations.append((_input, ('*', c)))

        _output = set([p[1] for p in propagations if p[0] == _input])
        cY2 = set(o[1] for o in _output)

        for c in cY2:
            _output_cY2 = set([o[0] for o in _output if o[1] == c])

            if len(_output_cY2) == len(colours) - 1:
                for c2 in colours:
                    if c2 not in _output_cY2:
                        propagations.append((_input, ('!' + c2, c)))
                    else:
                        try:
                            propagations.remove((_input, (c2, c)))
                        except Exception:
                            continue

    return propagations


def getPropagations(tileID):
    tiles = json.load(open("tilings_2.json"))
    propagations = []

    tile = {int(k): v for k, v in tiles[tileID].items()}

    for _input in [('a', 'b'), ('a', 'a')]:
        if 2 in tile[1] and _input[0] == _input[1]:
            continue

        colouring = dict.fromkeys(tile.keys(), '')
        colouring[1] = _input[0]
        colouring[2] = _input[1]

        queue = []
        for v in tile[1]:
            if v != 2: queue.append(v)

        for p in simplifyPropagations(findPropagation(queue, colouring, ['a', 'b', 'c'], tile), ['a', 'b', 'c']):
            propagations.append((p[0][0], p[0][1], p[1][0], p[1][1]))
    return propagations


def getSimplePropagations(tileID):
    tiles = json.load(open("tilings_2.json"))
    propagations = []

    tile = {int(k): v for k, v in tiles[tileID].items()}

    for _input in [('a', 'b'), ('a', 'a')]:
        if 2 in tile[1] and _input[0] == _input[1]:
            continue

        colouring = dict.fromkeys(tile.keys(), '')
        colouring[1] = _input[0]
        colouring[2] = _input[1]

        queue = []
        for v in tile[1]:
            if v != 2: queue.append(v)

        for p in findPropagation(queue, colouring, ['a', 'b', 'c'], tile):
            # here I am removing the duplicate propagations that mean the exact same thing
            # ie (a,a) ~> (b, b) is the same as (a, a) ~> (c, c)
            # and (a,a) ~> (b, c) is the same as (a, a) ~> (c, b)
            # and (a,a) ~> (a, b) is the same as (a, a) ~> (a, c)
            # and (a,a) ~> (b, a) is the same as (a, a) ~> (c, a)
            # hence when one is present the other always is
            if p[0][0] == 'a' and p[0][1] == 'a':
                if p[1][0] == 'c' or (p[1][0] == 'a' and p[1][1] == 'c'):
                    continue
            propagations.append((p[0][0], p[0][1], p[1][0], p[1][1]))
    return propagations


def getConcretePropagations(tileID):
    tiles = json.load(open("tilings_2.json"))
    propagations = []

    tile = {int(k): v for k, v in tiles[tileID].items()}

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
            print(p)
            propagations.append((p[0][0], p[0][1], p[1][0], p[1][1]))
    return propagations


if __name__ == "__main__":
    tiles = json.load(open("tilings_2.json"))

    for tileID in tiles:
        tile = {int(k): v for k, v in tiles[tileID].items()}

        print("\n---------------------------------\n3-colour Propagations for " + tileID)

        for _input in [('a', 'b'), ('a', 'a')]:
            if 2 in tile[1] and _input[0] == _input[1]:
                continue

            colouring = dict.fromkeys(tile.keys(), '')
            colouring[1] = _input[0]
            colouring[2] = _input[1]

            queue = []
            for v in tile[1]:
                if v != 2: queue.append(v)

            for p in simplifyPropagations(findPropagation(queue, colouring, ['a', 'b', 'c'], tile), ['a', 'b', 'c']):
                print("(" + p[0][0] + ", " + p[0][1] + ") ~> (" + p[1][0] + ", " + p[1][1] + ")")
