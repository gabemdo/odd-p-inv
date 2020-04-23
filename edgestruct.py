class EdgeStruct:

    def __init__(self,d):
        self.d = d
        self.e = [{j:None for j in range(self.d) if not i&(1<<j)} for i in range(1<<d)]
        self.max_size = (1<<(d-1)) * d

    def str_v(self,v):
        s = ""
        for i in range(self.d):
            s += str((v>>i)%2)
        return s

    def __str__(self):
        s = ""
        for ed in self:
            s += "Vertex: "+self.str_v(ed[0])+" d: "+self.str_v(1<<ed[1])+" val: "+str(self[ed])+"\n"
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
        for i in range(1<<self.d):
            for j in self.e[i]:
                yield (i,j)

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
        return self.__len__() == self.max_size

    def nonzero(self):
        assert (self.is_full()), "Not all edges have values."
        for edge in self:
            if self[edge] != 0:
                return True
        return False





