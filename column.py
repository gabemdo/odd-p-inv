import algebra as alg

class Row:
    #This is a container for the labeling of matrices in the chain complex.
    #It depends on the number of circles in the resolutions at each vertex.
    #These sizes are provided by the list res_len.
    #The iterator goes through all values ordered by vertex number, then group element

    def __init__(self,d,res_len):
        self.r = 0
        self.d = d
        self.res_len = res_len

    def __call__(self,r):
        self.r = r
        
    def v(self,v,i):
        return (v>>i)%2

    def __iter__(self):
        v = (1<<self.r)-1
        end = (1<<self.d)
        while v < end:
            for g in range(1<<self.res_len[v]):
                yield (v,g)
            if not v:
                break
            z = 0   #number of terminal zeros in binary representation of i
            u = -1   #number of ones ina row at end of binary representation of i
            while (v>>z) % 2 == 0:
                z += 1
            while (v>>(z+u+1)) % 2 == 1:
                u += 1
            v += (1<<z)+(1<<u)-1



def bin_pr(i,n):
    s = ""

    while i:
        b = "1" if i%2 else "0"
        s = b + s
        i >>= 1
    return "0"*(n-len(s))+s

def iterate_choose(n,k):
    i = (1<<k)-1
    end = (1<<n)
    count = 0
    if i == 0:
        print(bin_pr(i,n))
        print(1)
    else:
        while i < end:
            print(bin_pr(i,n))
            z = 0
            u = -1
            j = i
            while j % 2 == 0:
                z += 1
                j >>= 1
            while j % 2 == 1:
                u += 1
                j >>= 1
            i += (1<<z)+(1<<u)-1
            if z:
                print()
            count += 1
        print(count)
