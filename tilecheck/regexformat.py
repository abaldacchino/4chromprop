from sympy import *


def regexToMult(regex, tile_alphabet_map):
    """ Method converts a regular expression into a parsable python expres
    :param regex: String regular expression (not containing * Kleene operator)
    :param tile_alphabet_map: Map from alphabet to corresponding tile
    :return: Regular expression modified to be in terms of * for concatenation and + for choice
    """
    result = ""
    for i, letter in enumerate(regex):
        result += letter
        if (letter in tile_alphabet_map.keys() or letter in ['x', 'y', ')']) and (
                i != len(regex) - 1 and regex[i + 1] != ')' and regex[i + 1] != '+'):
            result += "*"

    return result


def listLettersToTiles(l, tile_alphabet_map):
    """
    Converts a list of letters to their corresponding tiles
    :param l: list of letters
    :param tile_alphabet_map: map from letters to tiles
    :return: list of strings (tiles)
    """
    result = []
    for i in l:
        result += [tile_alphabet_map[i]]
    return result


def listPairsToTiles(l, tile_alphabet_map):
    """
    Converts a list of pairs of letters to corresponding pairs of tiles
    :param l: list of pairs of letters
    :param tile_alphabet_map: map from letters to tiles
    :return: list of pairs of strings (tiles)
    """
    result = []
    for (a, b) in l:
        result += [(tile_alphabet_map[a], tile_alphabet_map[b])]
    return result


if __name__ == "__main__":
    # regex containing the expression for 4-chrom cyclizations, using letters for tiles
    regex = "(lx+pxg)(j+i)+qxgr+((q+k+gx(e+gk))x+s+t+m+o+kxn+rq+exd+gx(a+h+c+b+rg+pxg(" \
            "j+i)+ne+(d+ng)x(k+ge)+exn+kxd+g(s+q+t+m+o+kxn+rq+exd)+(j+i)l)+(j+i+gxg(j+i))(p+(" \
            "g)xgp)+r)g+gx(r+n+(d+ng)xg+(j+i)(gx(g+g(" \
            "p+g)+y)+g+p)+j+i+s+q+t+m+o+rq+exd+kx(n+y+gg)+g(a+h+c+b+kxd+rg+pxg(j+i)+exn+n(" \
            "e+g)+(d+ng)x(k+ge)+j+i+g+d+(d+ng+e)x(y+gg)+e)+y+k)+e+ex(y+gg+n)+d+n(" \
            "g+e)+(d+ng)x(y+g(g+e)+k)+g+(j+i+gxg(j+i))(gx(g+y+g(" \
            "g))+g+l)+j+i+a+h+c+b+kxd+(a+h+c+b+kxd+(r+s+t+m+o+kxn+rq+exd+gx(a+h+c+b+rg+pxg(j+i)+ne+(" \
            "d+ng)x(k+ge)+exn+kxd+g(s+q+t+m+o+kxn+rq+exd)+(j+i)l)+(j+i+gxg(j+i))(p+gx(" \
            "g)p))g+pxg(j+i)+exn+ne+(d+ng)x(k+ge)+gx(s+q+t+m+o+rq+exd+kxn+g(a+h+c+b+kxd+r(" \
            "g)+pxg(j+i)+exn+ne+(d+ng)x(k+ge))+(j+i)(p+gxgp))+(j+i+gxg(j+i))l)x(y+(" \
            "g)g)"
    # correspondence between letters ant tiles.
    # Note g actually stands for both f and g since they have the same propagation: this makes simplification easier
    tile_alphabet_map = {'a': "ABL", 'b': "BAL", 'c': "VVL", 'd': "VDL", 'e': "DVL", 'f': "AIVL", 'g': "VIAL",
                         'h': "DDdL", 'i': "AIVdL", 'j': "VIAdL", 'k': "AVL", 'l': "ADL", 'm': "BDL", 'n': "VAL",
                         'o': "DBL",
                         'p': "DDL", 'q': "AIBL", 'r': "BIAL", 's': "HL", 't': "DAdL"}

    # x is an even amount of Vis, y is $ symbol
    a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, x, y = symbols("a, b, c, d, e, f, g, h, i, j, k, l, m, n, "
                                                                             "o, p, q, r, s, t, x, y", commutative=False)
    expr = (l*x+p*x*g)*(j+i)+q*x*g*r+((q+k+g*x*(e+g*k))*x+s+t+m+o+k*x*n+r*q+e*x*d+g*x*(a+h+c+b+r*g+p*x*g*(j+i)+n*e+(d+n*g)*x*(k+g*e)+e*x*n+k*x*d+g*(s+q+t+m+o+k*x*n+r*q+e*x*d)+(j+i)*l)+(j+i+g*x*g*(j+i))*(p+g*x*g*p)+r)*g+g*x*(r+n+(d+n*g)*x*g+(j+i)*(g*x*(g+g*(p+g)+y)+g+p)+j+i+s+q+t+m+o+r*q+e*x*d+k*x*(n+y+g*g)+g*(a+h+c+b+k*x*d+r*g+p*x*g*(j+i)+e*x*n+n*(e+g)+(d+n*g)*x*(k+g*e)+j+i+g+d+(d+n*g+e)*x*(y+g*g)+e)+y+k)+e+e*x*(y+g*g+n)+d+n*(g+e)+(d+n*g)*x*(y+g*(g+e)+k)+g+(j+i+g*x*g*(j+i))*(g*x*(g+y+g*g)+g+l)+j+i+a+h+c+b+k*x*d+(a+h+c+b+k*x*d+(r+s+t+m+o+k*x*n+r*q+e*x*d+g*x*(a+h+c+b+r*g+p*x*g*(j+i)+n*e+(d+n*g)*x*(k+g*e)+e*x*n+k*x*d+g*(s+q+t+m+o+k*x*n+r*q+e*x*d)+(j+i)*l)+(j+i+g*x*g*(j+i))*(p+g*x*g*p))*g+p*x*g*(j+i)+e*x*n+n*e+(d+n*g)*x*(k+g*e)+g*x*(s+q+t+m+o+r*q+e*x*d+k*x*n+g*(a+h+c+b+k*x*d+r*g+p*x*g*(j+i)+e*x*n+n*e+(d+n*g)*x*(k+g*e))+(j+i)*(p+g*x*g*p))+(j+i+g*x*g*(j+i))*l)*x*(y+g*g)
    expr = expr.expand()

    strExpr = str(expr)
    # storing each multiplicative term in expressions
    expressions = []
    buffer = ""
    for i, letter in enumerate(strExpr):
        if letter == '+':
            expressions += [buffer]
            buffer = ""
        elif letter != ' ' and (not (letter == '*' and strExpr[i+1] == 'y')) and not letter == 'y':
            buffer += letter
    expressions += [buffer]
    expressions = list(set(expressions))  # removing any duplicates
    print("initial: ")
    print(expressions)
    print(len(expressions))

    # Language generated by e is contained in language generated by e*x
    toRemove = []
    for e in expressions:
        match = e+"*x"
        for f in expressions:
            if f == match:
                toRemove += [e]

    for e in toRemove:
        expressions.remove(e)
    print(len(expressions))

    # language s1*x*g**2*s2 is contained in language s1*x*s2
    toRemove = []
    for e in expressions:
        index = e.find("x*g**2")
        if index != -1:
            match = e[:index] + "x" + e[index+6:]
            for f in expressions:
                if f == match:
                    toRemove += [e]

    for e in toRemove:
        expressions.remove(e)
    print(len(expressions))
    print(expressions)

    # language s1*x*g**2*s2 is contained in language s1*x*s2
    toRemove = []
    for e in expressions:
        found = false
        for letter in tile_alphabet_map.keys():
            if letter in e and letter != 'g':
                found = true

        if found and e[0] == 'g':
            toRemove += [e]
    print(toRemove)
    for e in toRemove:
        expressions.remove(e)
    print(len(expressions))

    # s1 g*x*g s2 and s1 s2 can be described using s1 x s2
    toRemove = []
    toAdd = []
    for e in expressions:
        index = e.find("g*x*g")
        if index != -1:
            match = e[:index-1] + e[index + 5:]
            print(e)
            print(match)
            for f in expressions:
                if f == match:
                    toRemove += [e]
                    toRemove += [f]
                    toAdd += [e[:index] + "x" + e[index + 5:]]
    for e in toRemove:
        expressions.remove(e)
    expressions += toAdd
    print(len(expressions))

    # switching x*g to g*x
    toRemove = []
    toAdd = []
    for i, e in enumerate(expressions):
        index = e.find("x*g")
        if index != -1:
            expressions[i] = e[:index] + "g*x" + e[index + 3:]

    expressions = list(set(expressions))
    print(len(expressions))

    expressions.sort()
    print(expressions)

    # We match expression based on the cases defined in the characterization
    remaining = []
    matchCaseI = []
    matchCaseII = []
    matchCaseIII = []
    matchCaseIV = []
    matchCaseV = []
    matchCaseVI = []
    for e in expressions:
        if e[2:] == "x":        # t * x
            matchCaseI += [e[0]]
        # t1 * g * x * t2 * x
        elif len(e) == 9 and e[2:5] == "g*x" and e[8] == 'x':
            matchCaseII += [(e[0], e[6])]
        # t1 * t2 * g * x
        elif len(e) == 7 and e[4:] == "g*x":
            matchCaseIII += [(e[0], e[2])]
        # t * g * x
        elif len(e) == 5 and e[2:] == "g*x":
            matchCaseIV += [e[0]]
        # t1 * x * t2 * x
        elif len(e) == 7 and e[2] == 'x' and e[6:] == "x":
            matchCaseV += [(e[0], e[4])]
        # t1 * t2 * x
        elif len(e) == 5 and e[4:] == "x":
            matchCaseVI += [(e[0], e[2])]
        else:
            remaining += [e]

    print("Case i")
    print(listLettersToTiles(matchCaseI, tile_alphabet_map))
    print("Case ii")
    print(listPairsToTiles(matchCaseII, tile_alphabet_map))
    print("Case iii")
    print(listPairsToTiles(matchCaseIII, tile_alphabet_map))
    print("Case iv")
    print(listLettersToTiles(matchCaseIV, tile_alphabet_map))
    print("Case v")
    print(listPairsToTiles(matchCaseV, tile_alphabet_map))
    print("Case vi")
    print(listPairsToTiles(matchCaseVI, tile_alphabet_map))
    print("Remaining")
    print(remaining)

