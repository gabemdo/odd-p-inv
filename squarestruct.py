class SquareStruct:

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
        for i in range(1<<self.d):
            for j in self.e[i]:
                yield (i,*j)

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


