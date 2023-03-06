import json

insides = ["AA", "AB", "AV", "AD",
           "BA", "BB", "BV", "BD",
           "VA", "VB", "VV", "VD",
           "DA", "DB", "DV", "DD",
           "AIB", "AIV", "BIA", "VIA",
           "H"]

frames = ["L", "dL"]





def findAllUniquePropagations(tilings, whichType):
    propagations = []
    for f in frames:
        for i in insides:
            tile = i + f
            for prop in tilings[tile][whichType]:
                propagations += [tuple(prop)]
    return set(propagations)


def findTilesSatisfying(tilings, whichType):
    satisfyingTiles = []
    for prop in findAllUniquePropagations(tilings, whichType):
        propTiles = []
        for f in frames:
            for i in insides:
                tile = i + f
                if list(prop) in tilings[tile][whichType]:
                    propTiles += [tile]
        satisfyingTiles += [(prop, propTiles)]

    return satisfyingTiles


def partitionTilesOnPropagations(tilings):
    # type shouldn't matter here!!!!!!!
    covered_tiles = []
    partitions = []
    for tile in tiles:
        tilings[tile]["simple propagations"].sort()

    for t1 in tiles:
        if t1 in covered_tiles:
            continue
        eqv_t1 = [t1]
        covered_tiles += [t1]
        for t2 in tiles:
            if t1 == t2:
                continue
            if tilings[t1]["simple propagations"] == tilings[t2]["simple propagations"]:
                eqv_t1 += [t2]
                covered_tiles += [t2]
        partitions += [eqv_t1]
    return partitions


if __name__ == "__main__":
    print("Weeee data!!")
    with open('info.json') as json_file:
        tilings = json.load(json_file)

        uniqueProps = (findAllUniquePropagations(tilings, "simple propagations"))
        for pair in findTilesSatisfying(tilings, "aggregated propagations"):
            print("Tiles for", pair[0])
            print(pair[1])
        partitions = partitionTilesOnPropagations(tilings)
        total = 0
        for p in partitions:
            print(p)
            total += len(p)
        print(total)    # should be 42

        print("Tiles with only 1 propagation")
        for tile in tiles:
            if len(tilings[tile]["aggregated propagations"]) == 1:
                print(tile, tilings[tile]["aggregated propagations"])

        print("Tiles with 2 propagations")
        for tile in tiles:
            if len(tilings[tile]["aggregated propagations"]) == 2:
                print(tile, tilings[tile]["aggregated propagations"])

        print("Tiles with 3 propagations")
        for tile in tiles:
            if len(tilings[tile]["aggregated propagations"]) == 3:
                print(tile, tilings[tile]["aggregated propagations"])

        print("Tiles with 4 or more?!?!?!?!")
        for tile in tiles:
            if len(tilings[tile]["aggregated propagations"]) > 4:
                print(tile, tilings[tile]["aggregated propagations"])

        print("")
        print(len(findAllUniquePropagations(tilings, "aggregated propagations")))
