class EdgeStruct:
    #This provides a container and iterator for the edges in a cube complex.

    def __init__(self,d,value = None):
        self.d = d
        self.default = value
        self.e = [{i:value for i in range(self.d) if not v&(1<<i)} for v in range(1<<d)]
        self.max_size = (1<<(d-1)) * d

    def str_v(self,v):
        s = ""
        for i in range(self.d):
            s += str((v>>i)%2)
        return s

    def __str__(self):
        s = ""
        for ed in self:
            s += "\nVertex: "+self.str_v(ed[0])+" d: "+self.str_v(1<<ed[1])+" val: "+str(self[ed])+"\n"
        return s

    def it_print(self):
        for edge in self:
            print("Vertex:", edge[0], "del:", edge[1],"value:",self.e[edge[0]][edge[1]])

    def __getitem__(self,edge):
        assert (len(edge) == 2), "An edge should be 2-tuple of integers: (vertex,change)."
        assert (edge[1] < self.d), "The changed is out of range."
        assert ((~edge[0])&(1<<edge[1])), "The change should be a place where the vertex is zero (in binary)."
        return self.e[edge[0]][edge[1]]

    def __contains__(self,edge):
        assert (len(edge) == 2), "An edge should be 2-tuple of integers: (vertex,change)."
        assert (edge[1] < self.d), "The changed is out of range."
        assert ((~edge[0])&(1<<edge[1])), "The change should be a place where the vertex is zero (in binary)."
        return self[edge] != None

    def __iter__(self):
        for v in range(1<<self.d):
            for i in self.e[v]:
                yield (v,i)

    def __len__(self):
        #update later to an attribute that is updated with different methods
        i = 0
        for edge in self:
            if self[edge] != None:
                i += 1
        return i

    def __delitem__(self, edge):
        assert (len(edge) == 2), "An edge should be 2-tuple of integers: (vertex,change)."
        assert (edge[1] < self.d), "The changed is out of range."
        assert ((~edge[0])&(1<<edge[1])), "The change should be a place where the vertex is zero (in binary)."
        self.e[edge[0]][edge[1]] = None

    def __setitem__(self, edge, val):
        assert (len(edge) == 2), "An edge should be 2-tuple of integers: (vertex,change)."
        assert (edge[1] < self.d), "The changed is out of range."
        assert ((~edge[0])&(1<<edge[1])), "The change should be a place where the vertex is zero (in binary)."
        self.e[edge[0]][edge[1]] = val

    def get(self, edge, value):
        if edge in self:
            return self[edge]
        return value

    def is_full(self):
        #This version until controlled __len__() is implemented
        for edge in self:
            if self[edge] == None:
                return False
        return True

    def is_fully_changed(self):
        for edge in self:
            if self[edge] == self.default:
                return False
        return True

    def nonzero(self):
        assert (self.is_full()), "Not all edges have values."
        for edge in self:
            if self[edge] != 0:
                return True
        return False


class SquareStruct:
    #This is a container for the faces of the cube for the Khovanov homology.
    #It provides an iterator that is used in the inductive definition of the sign assignment for edges in the odd Khovanov homology.
    #The iterator does not iterator through all faces.

    def __init__(self,d):
        self.d = d
        self.e = [{(j,k):None for j in range(d-1) for k in range(j+1,d) if not (i)&((1<<j)+(1<<k)) } for i in range(1<<d)]
        self.max_size = (1<<(d-2)) * d * (d-1) // 2 if d >= 2 else 0

    def str_v(self,v):
        s = ""
        for i in range(self.d):
            s += str((v>>i)%2)
        return s

    def __getitem__(self,square):
        assert (len(square) == 3), "A square should be a 3-tuple of integers: (vertex, delta1, delta2)."
        assert (square[1] < self.d and square[2] < self.d), "The change is out of range."
        assert (square[1] != square[2]), "The two edges must be different."
        assert ((~square[0])&((1<<square[1])+(1<<square[2]))), "The changes should be at places where the vertex is zero."
        if square[1] < square[2]:
            return self.e[square[0]][(square[1],square[2])]
        return self.e[square[0]][(square[2],square[1])]

    def __contains__(self,square):
        assert (len(square) == 3), "A square should be a 3-tuple of integers: (vertex, delta1, delta2)."
        assert (square[1] < self.d and square[2] < self.d), "The change is out of range."
        assert (square[1] != square[2]), "The two edges must be different."
        assert ((~square[0])&((1<<square[1])+(1<<square[2]))), "The changes should be at places where the vertex is zero."
        if square[1] < square[2]:
            return self.e[square[0]][(square[1],square[2])] != None
        return self.e[square[0]][(square[2],square[1])] != None

    def __iter__(self):
        #change to inductive order
        for j in range(1,self.d):
            for i in range(j):
                for b in range(1<<i):
                    for k in range(1<<(j-i-1)):
                        v = b + (k<<(i+1))
                        yield v, i, j
        #v = A0B0C where |C| = i, |B| + 1 + |C| = j
        #A = 0
        #B = k
        #C = b


    def __len__(self):
        #update later to attribute that is self updating
        i = 0
        for sq in self:
            if self[sq] != None:
                i += 1
        return i

    def __str__(self):
        s = ""
        for sq in self:
            s += "Vertex: "+self.str_v(sq[0])+" delta: "+self.str_v((1<<sq[1])+(1<<sq[2]))+" val: "+str(self[sq])+"\n"
        return s

    def __delitem__(self,square):
        if square in self:
            self[square] = None

    def __setitem__(self,square,value):
        assert (len(square) == 3), "A square should be a 3-tuple of integers: (vertex, delta1, delta2)."
        assert (square[1] < self.d and square[2] < self.d), "The change is out of range."
        assert (square[1] != square[2]), "The two edges must be different."
        assert ((~square[0])&((1<<square[1])+(1<<square[2]))), "The changes should be at places where the vertex is zero."
        if square[1] < square[2]:
            self.e[square[0]][(square[1],square[2])] = value
        else: 
            self.e[square[0]][(square[2],square[1])] = value

    def get(self,square,value):
        if square in self:
            return self[square]
        return value

    def is_full(self):
        return self.max_size == len(self)

    def nonzero(self):
        assert (self.is_full()), "Not all squares have values."
        for sq in self:
            if self[sq] != 0:
                return True
        return False


class Vertical:
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


"""
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
"""
