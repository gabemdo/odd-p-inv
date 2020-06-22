import braid as b
import sys

def prepare_list(file):
    braid_list = []
    with open(file) as f:
        for line in f:
            sl = line.strip().split()
            if sl:
                braid_list.append((sl[0],[int(sl[i]) for i in range(1,len(sl))]))
    return braid_list

def flypable(braid_word):
    m = comp_m(braid_word)
    i = 0
    minus_count = 0
    for gen in braid_word:
        if gen == -m:
            minus_count += 1
    if minus_count != 1:
        return False, [], []
    minus_index = braid_word.index(-m)
    ordered = braid_word[minus_index:] + braid_word[:minus_index]
    A = []
    B = []
    i = 1
    while i < len(ordered) and ordered[i] != m:
        A.append(ordered[i])
        i += 1
    k = 0
    while i < len(ordered) and ordered[i] == m:
        i += 1
        k += 1
    while i < len(ordered) and ordered[i] != m:
        B.append(ordered[i])
        i += 1
    if i != len(ordered):
        return False, [], []
    if A == B:
        return False, [], []
    return True, ordered, [-m] + B + [m]*k + A

def mirror(braid_pair):
    name = braid_pair[0]
    if name[0] == 'm':
        mname = name[1:]
    else:
        mname = 'm' + name
    mir = [-gen for gen in braid_pair[1]]
    return mname, mir


def comp_m(braid):
        #return the width - 1 of the braid
        m = 0
        for s in braid:
            c = abs(s)
            if c > m:
                m = c
        return m


def driver(braid_list):
    for braid in braid_list:
        flyable, br, flype = flypable(braid[1])
        if flyable:
            print(braid[0])
            print(br)
            brd = b.Braid(br)
            mult,qhom,ihom = brd.comp_inv(False)
            if mult == 1:
                print("    Odd  Inv Zero\n")
            elif mult == 0:
                print("    Odd  Inv NonZero\n")
            else:
                print("    Odd  Inv Torsion: {}\n".format(mult))
            print(flype)
            flp = b.Braid(flype)
            mult,qhom,ihom = flp.comp_inv(False)
            if mult == 1:
                print("    Odd  Inv Zero\n")
            elif mult == 0:
                print("    Odd  Inv NonZero\n")
            else:
                print("    Odd  Inv Torsion: {}\n".format(mult))
            print()
        mname, mir = mirror(braid)
        flyable, br, flype = flypable(mir)
        if flyable:
            print(mname)
            print(br)
            brd = b.Braid(br)
            mult,qhom,ihom = brd.comp_inv(False)
            if mult == 1:
                print("    Odd  Inv Zero\n")
            elif mult == 0:
                print("    Odd  Inv NonZero\n")
            else:
                print("    Odd  Inv Torsion: {}\n".format(mult))
            print(flype)
            flp = b.Braid(flype)
            mult,qhom,ihom = flp.comp_inv(False)
            if mult == 1:
                print("    Odd  Inv Zero\n")
            elif mult == 0:
                print("    Odd  Inv NonZero\n")
            else:
                print("    Odd  Inv Torsion: {}\n".format(mult))
            print()

if __name__ == "__main__":
    file = sys.argv[1]
    braid_list = prepare_list(file)
    driver(braid_list)
