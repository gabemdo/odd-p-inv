class Braid:
    def __init__(self,pres):
        assert (len(pres) > 0)
        self.presentation = pres
        self.braid_index = self.braid_to_braid_index(self.presentation)
        self.crossing_number = len(self.presentation)
        self.planar_diagram = []
        self.inv_r = -1
        self.inv_vertex = []
        self.inv_resolution = []
        self.vertices_before_inv = []
        self.vertices_after_inv = []
    
    def crossing_locus(self,x):
        if x < 0:
            return -1*x
        return x
    
    def braid_to_braid_index(self,braid):
        i = 0
        for x in braid:
            c = self.crossing_locus(x)
            if c > i:
                i = c
        return i+1
    
    def get_index(self):
        return self.braid_index
    
    def get_planar_diagram(self):
        if not self.planar_diagram:
            n = self.braid_index
            a = [i+1 for i in range(n)]
            b = [0 for i in range(n)]
            last_b = [0 for i in range(n)]
            i = n + 1
            d = []
            l = len(self.presentation)
            for j in range(l):
                x = self.presentation[j]
                b[:] = a[:]
                c = self.crossing_locus(x)
                b[c-1] = i
                i += 1
                b[c] = i
                i += 1
                if j == l-1:
                    last_b[:] = b[:]
                    b = [i+1 for i in range(n)]
                if x > 0:
                    d.append([a[c-1],b[c-1],b[c],a[c]])
                else:
                    d.append([a[c],a[c-1],b[c-1],b[c]])
                a[:] = b[:]
            for x in range(l):
                for j in range(1,4):
                    if d[x][j] in last_b:
                        if self.presentation[x] > 0:
                            d[x][j] = j - 1 + self.presentation[x] 
                        else:
                            d[x][j] = j - 2 - self.presentation[x]
            self.planar_diagram = d
        return self.planar_diagram
    
    def vertices(self,n=0):
        assert (n >= 0)
        if n == 0 :
            n = self.crossing_number
        vertices = []
        for i in range(2**n):
            x = i
            vertex = []
            for j in range(n):
                vertex.append(x%2)
                x //= 2
            vertices.append(vertex)
        return vertices

    def zero_smoothing(self,crossing):
        assert (len(crossing) == 4)
        return (crossing[0],crossing[1]),(crossing[2],crossing[3])
    
    def one_smoothing(self,crossing):
        assert (len(crossing) == 4)
        return (crossing[1],crossing[2]),(crossing[3],crossing[0])
    
    def sym_to_cycles(self,sym):
        cycles = []
        while sym:
            a,b = sym.pop()
            cycle = [a]
            while cycle[0] != b:
                cycle.append(b)
                for i in range(len(sym)):
                    if sym[i][0] == b:
                        a,b = sym.pop(i)
                        break
                    elif sym[i][1] == b:
                        b,a = sym.pop(i)
                        break
            cycles.append(cycle)
        return cycles
    
    def resolution(self,vertex):
        assert (len(vertex) == self.crossing_number)
        if not self.planar_diagram:
            self.get_planar_diagram()
        sym = []
        for i in range(self.crossing_number):
            if vertex[i] == 0:
                a,b = self.zero_smoothing(self.planar_diagram[i])
            else:
                a,b = self.one_smoothing(self.planar_diagram[i])
            sym.append(a)
            sym.append(b)
        return self.sym_to_cycles(sym)
    
    def get_inv_r(self):
        if self.inv_r < 0:
            i = 0
            for x in self.presentation:
                if x < 0:
                    i += 1
            self.inv_r = i
        return self.inv_r
    
    def get_inv_vertex(self):
        if not self.inv_vertex:
            vertex = []
            for x in self.presentation:
                if x > 0:
                    vertex.append(0)
                else:
                    vertex.append(1)
            self.inv_vertex = vertex
        return self.inv_vertex
    
    def get_inv_resolution(self):
        if not self.inv_resolution:
            vertex = self.get_inv_vertex()
            self.inv_resolution = self.resolution(vertex)
        return self.inv_resolution
    
    def to_vertices(self,vertices):
        edges = []
        for vertex in vertices:
            for i in range(len(vertex)):
                edge = [0 for i in range(len(vertex))]
                if vertex[i] == 0:
                    edge[:] = vertex[:]
                    edge[i] = 1
                    edges.append(edge)
        return edges
    
    
    def from_vertices(self,vertices):
        edges = []
        for vertex in vertices:
            for i in range(len(vertex)):
                edge = []
                if vertex[i] == 1:
                    edge[:] = vertex[:]
                    edge[i] = 0
                    edges.append(edge)
        return edges
    
    def get_vertices_after_inv(self):
        return self.to_vertices([self.get_inv_vertex()])
    
    def get_vertices_before_inv(self):
        return self.from_vertices([self.get_inv_vertex()])
    
    def vertex_explanation(self):
        print("Braid:",self.presentation)
        print("Vertices that map to cycle's vertex:")
        from_v = self.get_vertices_before_inv()
        for vertex in from_v:
            print("    Vertex:", vertex,
                  "Resolution:", 
                  self.resolution(vertex))
        if not from_v:
            print("    NONE")
        print("The distinguished vertex:")
        print("    Vertex:", self.get_inv_vertex(), "Resolution:", self.get_inv_resolution())

        print("Vertices that the vertex maps to:")
        to_v = self.get_vertices_after_inv()
        for vertex in to_v:
            print("    Vertex:", vertex,
                  "Resolution:", 
                  self.resolution(vertex))
        if not to_v:
            print("    NONE")

b = Braid([2,-3,2])
#print(b.get_planar_diagram())
for vertex in b.vertices():
    print(b.resolution(vertex))
    
print(123,b.get_inv_resolution())
b.vertex_explanation()