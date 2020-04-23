class Braid:
    def __init__(self,pres):
        #The presentation of a braid is given as a list of non-zero integers
        assert (len(pres) > 0)
        self.presentation = pres
        #The braid index is inferred. Perhaps it is better to take it as an argument.
        self.braid_index = self.braid_to_braid_index(self.presentation)
        self.crossing_number = len(self.presentation)
        self.planar_diagram = []
        self.inv_r = -1
        self.inv_vertex = []
        self.inv_resolution = []
        self.vertices_before_inv = []
        self.vertices_after_inv = []
        self.resolutions = [[] for _ in range(1 << self.crossing_number)]
    
    def crossing_locus(self,x):
        #This is the absolute value function
        #It converts a braid group generator to the location of the crossing
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
        #deprecate ?
        return self.braid_index
    
    def get_planar_diagram(self):
        if not self.planar_diagram:
            n = self.braid_index
            #The segments of braid are labeled.
            # a: segment labels, before the first crossing we start [1,2,...,n]
            a = [i+1 for i in range(n)]
            #After a crossing the new labels for each strand are given in b
            b = []
            #last_b: used to relabel the final segments to "close" the braid
            last_b = []
            # i: the index of the new segments, initialized for the first new segment, n+1
            i = n + 1
            # d: the planar diagram
            d = []
            
            l = len(self.presentation)
            for j in range(l):
                # x: current crossing
                x = self.presentation[j]
                # compute current strand labels
                b[:] = a[:]
                c = self.crossing_locus(x)
                b[c-1] = i
                i += 1
                b[c] = i
                i += 1
                
                #If last strand, save values of b before
                #make b initial strand ("close" the braid)
                if j == l-1:
                    last_b[:] = b[:]
                    b = [i+1 for i in range(n)]
                    
                #if positive crossing (`/,) label from upper left counterclockwise 
                if x > 0:
                    d.append([a[c-1],b[c-1],b[c],a[c]])
                #if negative crossing (,\') label from upper right counterclockwise
                else:
                    d.append([a[c],a[c-1],b[c-1],b[c]])
                    
                #shift labeling
                a[:] = b[:]
                
            #relabel middle segments at bottom to match relabelling at bottom to close the braid
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
        #produce the vertices as lists of ones and zeros 
        #index matches power, i.e., i vertex = sum (vertices[i][j]<<j)
        #vertices labelled 0,...,(2^n)-1 
        #n: dimension of hypercube
        #if no n, it becomes the crossing_number
        assert (n >= 0)
        if n == 0 :
            n = self.crossing_number
        vertices = []
        for i in range(1<<n):
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
            a,b = sym.pop(0)
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
        s = 0
        for i in range(self.crossing_number):
            s += vertex[i] << i 
        if self.resolutions[s]:
            return self.resolutions[s]
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
        self.resolutions[s] = self.sym_to_cycles(sym)
        return self.resolutions[s]
    
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
        if not self.vertices_after_inv:
            self.vertices_after_inv = self.to_vertices([self.get_inv_vertex()])
        return self.vertices_after_inv
    
    def get_vertices_before_inv(self):
        if not self.vertices_before_inv:
            self.vertices_before_inv = self.from_vertices([self.get_inv_vertex()])
        return self.vertices_before_inv
    
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
            
    def exterior_algebra(self,resolution):
        algebra = []
        for selection in self.vertices(len(resolution)):
            gen = []
            for i in range(len(selection)):
                if selection[i] == 1:
                    gen.append(resolution[i])
            algebra.append(gen)
        return algebra
    
    def induced_map(self,vmap,algebra=[]):
        if not algebra:
            m = 0
            for nm in vmap:
                if m < nm[0]:
                    m = nm[0]
            m += 1
            algebra = self.exterior_algebra(range(m))
        algebra_map = {(0,0):1}
        for i in range(1,len(algebra)):
            image = []
            p = 1
            for vec in algebra[i]:
                for arrow in vmap:
                    if vec == arrow[0]:
                        image.append(arrow[1])
                        p *= vmap[arrow]
            s = 1
            #signed sort image (s = 0 if degen:)    
            for j in range(len(image)):
                for k in range(j+1,len(image)):
                    if image[j] > image[k]:
                        s *= -1
                        image[j],image[k]=image[k],image[j]
                    elif image[j] == image[k]:
                        s = 0
                        image = []
                        break
                if not image:
                    break
            s *= p
            target = 0
            for vec in image:
                target += 2**vec
            algebra_map[(i,target)] = s
        return algebra_map
    
    def add_algebra_maps(self,map1,map2):
        for key in map2:
            map1[key] = map1.get(key,0) + map2[key]
        return map1
    
    def edge_to_split_map(self,edge):
        assert (len(edge) == self.crossing_number)
        delta = edge.index(-1)
        source = []
        target = []
        source[:] = edge[:]
        target[:] = edge[:]
        source[delta] = 0
        target[delta] = 1
        source_res = self.resolution(source)
        target_res = self.resolution(target)
        a = self.get_planar_diagram()[delta][0:2]
        bmap = {}
        d = [-1,-1]
        for j in range(len(target_res)):
            if a[1] in target_res[j]:
                assert (d[1] == -1)
                d[1] = j
            for i in range(len(source_res)):
                    if a[0] not in source_res[i]:
                        if source_res[i][0] in target_res[j]:
                            bmap[(i,j)] = 1
                    else:
                        if a[0] in target_res[j]:
                            bmap[(i,j)] = 1
                            assert (d[0] == -1)
                            d[0] = j
        imap = self.induced_map(bmap)
        map0 = { (key[0],key[1]+(1<<d[0]) if (key[1]>>d[0])%2 == 0 else 0) :
                -1*imap[key] if (key[1]>>d[0])%2 == 0 else 0 
                for key in imap}
        map1 = { (key[0],key[1]+(1<<d[1]) if (key[1]>>d[1])%2 == 0 else 0) :
                imap[key] if (key[1]>>d[1])%2 == 0 else 0 
                for key in imap}
        total_map = self.add_algebra_maps(map0,map1)
        return total_map
    
    def get_cycle(self,resolution,segment):
        for i in range(len(resolution)):
            if segment in resolution[i]:
                return i
        #Raise exception if here
        
    def square_sign(self,v00,a1,a2):
        #The input should be a vertex and the indices of two arrow changing from 0 to 1
        assert (v00[a1] == 0 and v00[a2] == 0 and a1 != a2)
        r00 = self.resolution(v00)
        c00 = len(r00)
        v11 = []
        v11[:] = v00[:]
        v11[a1] = v11[a2] = 1
        r11 = self.resolution(v11)
        c11 = len(r11)
        if c00 < c11:
            return 1 #Type A
        elif c11 < c11:
            return -1 #Type C
        pd = self.get_planar_diagram()
        t1 = pd[a1][2]
        h1 = pd[a1][0]
        t2 = pd[a2][2]
        h2 = pd[a2][0]
        ct1 = self.get_cycle(r00,t1)
        if t2 not in r00[ct1]:
            return -1 #Type C
        if h1 not in r00[ct1]:
            return 1 #Type A
        t1i = r00[ct1].index(t1)
        h1i = r00[ct1].index(h1)
        if t1i < h1i:
            if h2 in r00[ct1][t1i:h1i]:
                print("Look at resolution {}".format(r00))
                print("Look at cylce {}".format(ct1))
                print("That is: {}".format(r00[ct1]))
                print("Look at slice [{}:{}]".format(t1i,h1i))
                print("That is: {}".format(r00[ct1][t1i:h1i]))
                print("Looking for {}".format(h2))
                print(pd[2],pd[3])
                print((t1,h1),(t2,h2))
                return "Y" #-1Type Y
            return "X" #1Type X
        else:
            print(self.crossing_number)
            if h2 not in r00[ct1][t1i:h1i]:
                return "Y" #-1Type Y
            return "X" #1Type X
        
    
    def edge_map(self,source_res,target_res):
        #get resolutions of vertices
        if len(source_res) < len(target_res):
            org = []
            for i in range(len(target_res)):
                for j in range(len(source_res)):
                    if target_res[i][0] in source_res[j]:
                        org.append(j)
                        break
            doubled = []
            for i in range(len(org)):
                if org[i] in (org[:i] + org[i+1:]):
                    doubled.append(i)
            map1 = {(i,i+2**doubled[0]):1 for i in range(2**len(source_res))}
            map2 = {(i,i+2**doubled[1]):-1 for i in range(2**len(source_res))}
            return self.add_algebra_maps(map1,map2)
            pass #do split map
        else: 
            vmap = {}
            for i in range(len(source_res)):
                for j in range(len(target_res)):
                    if source_res[i][0] in target_res[j]:
                        vmap[(i,j)] = vmap.get((i,j),0)+1
                        break
            print(vmap)
            return self.induced_map(vmap)#,self.exterior_algebra(source_res))
