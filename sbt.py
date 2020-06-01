import braid as b
import sys

def driver():
    braid_word = []
    n = int(sys.argv[1])
    for i in range(2,n+2):
        braid_word.append(int(sys.argv[i]))
    braid = b.Braid(braid_word)
    mult,qhom,ihom = braid.comp_inv(False)
    if mult == 1:
        print("Inv Zero")
    elif mult == 0:
        print("Inv NonZero")
    else:
        print("Inv Torsion")


if __name__ == "__main__":
    driver()