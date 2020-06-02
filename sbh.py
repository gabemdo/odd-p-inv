import braid as b
import sys

def driver():
    braid_word = []
    n = int(sys.argv[1])
    for i in range(2,n+2):
        braid_word.append(int(sys.argv[i]))
    braid = b.Braid(braid_word)
    braid.comp_homology(1,-1,False)

if __name__ == "__main__":
    driver()