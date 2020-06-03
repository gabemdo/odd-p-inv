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

def driver(braid_list):
    for braid in braid_list:
        n = len(braid[1])
        nmin = 0
        for s in braid[1]:
            if s < 0:
                nmin += 1
        nplus = n - nmin
        ratio = min(nmin,nplus)/n
        print("Braid: {}, {}, negative ratio = {:.5}".format(braid[0],"" if n == 1 else braid[1], ratio))
        br = b.Braid(braid[1])
        br.comp_full_graded_homology()

if __name__ == "__main__":
    file = sys.argv[1]
    braid_list = prepare_list(file)
    driver(braid_list)