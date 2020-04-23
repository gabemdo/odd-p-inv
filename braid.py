class Braid:
    # Some convetions:
    # Implement x used for crossing index number in presentation
    # x ranges in self.d, the number of crossing or the dimention of the hypercube
    # Implement s used for crossing strand number in presentation
    # s ranges in self.b, the braid index or the number strands in the braid
    #           | self.presentation[x] | = s
    # height of vertex is number of 1s in binary, ranges 0 .. self.d
    # Vextex convention:
    ### vertices are numbered 0 .. (1<<self.d)-1
    # Edge convention:
    ### edges are pairs(vertex, i)
    ### it is assumed that ((vertex>>i) % 2 == 0) (that is, in binary the ith place is 0)
    ### this is the edge between (vertex) and vertex + (1<<i)
    # ABBREVIATIONS:
    ### inv : the invariant or associated vertices, heights, etc.
    ### res : resolution
    
    def __init__(self,pres):
        #The presentation of a braid is given as a list of non-zero integers
        assert (len(pres) > 0)
        self.presentation = pres
        #The braid index is inferred. Perhaps it is better to take it as an argument.
        self.b = self.comp_b(self.presentation)
        self.d = len(self.presentation)
        self.x_diagram = []
        self.closure = []
        self.inv_r = -1
        self.inv_v = -1
        self.edges_out_of_inv = []
        self.edges_into_inv = []
        self.res = [[] for _ in range(1 << self.d)]
        self.res_len = [0 for _ in range(1 << self.d)]
        self.edge_signs = []
        #print("Initialized")
    
    def abss(self,s):
        #This is the absolute value function
        #It converts a braid group generator to the location of the crossing
        if s < 0:
            return -1*s
        return s
    
    def comp_b(self,braid):
        i = 0
        for s in braid:
            c = self.abss(s)
            if c > i:
                i = c
        return i+1

    def height(self,v):
        r = 0
        while v:
            r += v%2
            v == v >> 1
        return r
        
    def get_x_diagram(self):
        if not self.x_diagram:
            n = self.b
            #The segments of braid are labeled.
            # a: segment labels, before the first crossing we start [1,2,...,n]
            a = [i+1 for i in range(n)]
            #After a crossing the new labels for each strand are given in b
            b = []
            # i: segment index, starting at the next segment, n+1
            i = n + 1
            # d: the planar diagram
            d = []
            
            l = len(self.presentation)
            for j in range(l):
                ####print(a)
                # x: current crossing
                x = self.presentation[j]
                # compute current strand labels
                b[:] = a[:]
                c = self.abss(x)
                b[c-1] = i
                i += 1
                b[c] = i
                i += 1
                
                #if positive crossing (`/,) label from upper left counterclockwise 
                if x > 0:
                    d.append([a[c-1],b[c-1],b[c],a[c]])
                #if negative crossing (,\') label from upper right counterclockwise
                else:
                    d.append([a[c],a[c-1],b[c-1],b[c]])
                #shift labeling
                a[:] = b[:]
            ####print(b)
            #close braid (initializes which segments at top connect to those at bottom)
            #side-effect: resolutions are listed in clockwise order
            ##under the convention that braids are read top to bottom and oriened upwardly
            for j in range(n):
                self.closure.append((j+1,b[j]))
            self.x_diagram = d
            #print(d)
        return self.x_diagram
        
    def v(self,v,i):
        return (v >> i) % 2
    
    def edges(self):
        #Check if edges_signs has already be initilized or computed
        if self.edge_signs:
            return self.edge_signs
        n = self.d
        vertices = []
        for i in range(1<<n):
            e = {}
            for j in range(n):
                if self.v(i,j) == 0:
                    e[j] = (1 if i < (1<<j) else 0)
            vertices.append(e)
        self.edge_signs = vertices
        return vertices
    
    def comp_edge_signs(self):
        ## Algorithm due to Shumakovitch
        #Make sure edge_signs is initialized
        if not self.edge_signs:
            self.edges()
        n = self.d
        #Check if sign assignment is trival (no squares)
        if n < 2:
            return self.edge_signs
        #Check if sign assignment has already been computed
        #print(self.edge_signs)
        assert (len(self.edge_signs) > 3)
        if self.edge_signs[2][0] != 0:
            return self.edge_signs
        for i in range(1,n):
            for j in range(i):
                for k in range(1<<j):
                    for m in range(1<<(i-j-1)):
                        v = k + (m << (j+1))
                        #print("v,j,i, [], k,m")
                        #print(v,j,i,[(v>>s)%2 for s in range(n)],k,m)
                        self.edge_signs[v + (1<<i)][j] = -1 * self.edge_signs[v][j] * self.square_sign(v,j,i)
        return self.edge_signs
    
    def get_edge_sign(self,edge):
        return self.comp_edge_signs()[edge[0]][edge[1]]
        
    def zero_smoothing(self,crossing):
        assert (len(crossing) == 4)
        return (crossing[0],crossing[1]),(crossing[2],crossing[3])
    
    def one_smoothing(self,crossing):
        assert (len(crossing) == 4)
        return (crossing[1],crossing[2]),(crossing[3],crossing[0])
    
    def sym_to_circles(self,sym):
        #Convert a list of segment attachments into full circles
        circles = []
        while sym:
            a,b = sym.pop(0)
            circle = [a]
            while circle[0] != b:
                circle.append(b)
                for i in range(len(sym)):
                    if sym[i][0] == b:
                        a,b = sym.pop(i)
                        break
                    elif sym[i][1] == b:
                        b,a = sym.pop(i)
                        break
            circles.append(circle)
        return circles
    
    def get_res(self,v):
        assert (v < (1<<self.d))
        if self.res[v]:
            return self.res[v]
        if not self.x_diagram:
            self.get_x_diagram()
        sym = []
        sym[:] = self.closure[:]
        for i in range(self.d):
            if self.v(v,i) == 0:
                a,b = self.zero_smoothing(self.x_diagram[i])
            else:
                a,b = self.one_smoothing(self.x_diagram[i])
            sym.append(a)
            sym.append(b)
        self.res[v] = self.sym_to_circles(sym)
        return self.res[v]
    
    def get_res_len(self,v):
        assert (v < (1<<self.d))
        if self.res_len[v]:
            return self.res_len[v]
        if self.res[v]:
            self.res_len[v] = len(self.res[v])
            return self.res_len[v]
        self.res_len[v] = len(self.get_res(v))
        return self.res_len[v]
        
    def get_inv_v(self):
        #return the vertex number for the invariant
        #if not self.inv_vertex:
        if self.inv_v < 0:
            #vertex = []
            v = 0
            #for x in self.presentation:
            #    if x > 0:
            #        vertex.append(0)
            #    else:
            #        vertex.append(1)
            for i in range(self.d):
                if self.presentation[i] < 0:
                    v += (1<<i)
            self.inv_v = v
        return self.inv_v

    def edges_out(self,v):
        #Returns the list of edges out of vertex 
        #format: list of tuples (vertex,i)
        #this is the edge between vertex and vertex + (1<<i)
        i = 0
        edges = []
        for i in range(self.d):
            if (v >> i) % 2 == 0:
                edges.append((v,i))
        return edges

    def edges_in(self,v):
        #Returns the lsit of edges in of vertex
        #format: list of tuples (v,i)
        #this is the edge between v and vertex 
        #such that v + (1<<i) = vertex
        i = 0
        edges = []
        for i in range(self.d):
            if (v>>i) % 2 == 1:
                edges.append((v - (1<<i),i))
        return edges

    def get_edges_out_of_inv(self):
        if not self.edges_out_of_inv:
            self.edges_out_of_inv = self.edges_out(self.get_inv_v())
        return self.edges_out_of_inv

    def get_edges_into_inv(self):
        if not self.edges_into_inv:
            self.edges_into_inv = self.edges_in(self.get_inv_v())
        return self.edges_into_inv
    
    def str_v(self,v):
        s = ""
        for i in range (self.d):
            s += str((v>>i)%2)
        return s
    
    def explanation(self):
        print("Braid:",self.presentation)
        print("Vertices that map to cycle's vertex:")
        edges_in = self.get_edges_into_inv()
        for edge in edges_in:
            print("  Vertex:", self.str_v(edge[0]),"del:",edge[1],
                  "Resolution:",
                  self.get_res(edge[0]))
            print("  Map:")
            self.print_map(self.edge_map(edge),edge)
        if not edges_in:
            print("    NONE")
        print("The distinguished vertex:")
        print("    Vertex:",self.str_v(self.get_inv_v()),
              "Resolution:",self.get_res(self.get_inv_v()))
        print("Vertices that the vertex maps to:")
        edges_out = self.get_edges_out_of_inv()
        for edge in edges_out:
            print("  Vertex:", self.str_v(edge[0]+(1<<edge[1])),"del:",edge[1],
                  "Resolution:",
                  self.get_res(edge[0]+(1<<edge[1])))
            print("  Map:")
            self.print_map(self.edge_map(edge),edge)

    def induced_map(self,v,vmap):
        xd = self.get_res_len(v)
        assert (xd == len(vmap))
        #print(vmap)
        algebra_map = [{} for _ in range(1<<xd)]
        for i in range(1<<xd):
            ls = []
            p = 1
            for j in range(xd):
                if (i>>j)%2 == 1:
                    ls.append(vmap[j][0])
                    p *= vmap[j][1]
            sign = 1
            #signed sort ls
            for k in range(len(ls)):
                for m in range(k+1,len(ls)):
                    if ls[k] > ls[m]:
                        sign *= -1
                        ls[k],ls[m] = ls[m],ls[k]
                    elif ls[k] == ls[m]:
                        sign = 0
                        ls = []
                        break
            target = 0
            #print(ls)
            for t in ls:
                target += (1<<t)
            algebra_map[i][target] = algebra_map[i].get(target,0) + p*sign
        return algebra_map
    
    def scalar_mult_algebra_map(self,scalar,mp):
        for i in range(len(mp)):
            for key in mp[i]:
                mp[i][key] *= scalar
        #not returned, stored in place in mp
    
    def add_algebra_maps(self,map1,map2):
        assert (len(map1) == len(map2))
        for i in range(len(map2)):
            for key in map2[i]:
                map1[i][key] = map1[i].get(key,0) + map2[i][key]
        #not returned, stored in place in map1
    
    def split_map(self,edge):
        assert (self.v(edge[0],edge[1]) == 0)
        source = edge[0]
        target = edge[0] + (1<<edge[1])
        source_res = self.get_res(source)
        target_res = self.get_res(target)
        a = self.get_x_diagram()[edge[1]]
        bmap = [None for _ in range(self.get_res_len(source))]
        #head,tail
        d = [-1,-1]
        for j in range(self.get_res_len(target)):
            #arrow tail in target ?
            if a[2] in target_res[j]:
                
                assert(d[1] == -1)
                d[1] = j
            for i in range(self.get_res_len(source)):
                #if not involved in the split
                if a[0] not in source_res[i]:
                    if source_res[i][0] in target_res[j]:
                        bmap[i] = (j,1)
                #if involved in the split
                else: 
                    #arrow head in target ?
                    if a[0] in target_res[j]:
                        bmap[i] = (j,1)
                        assert(d[0] == -1)
                        d[0] = j
        #print(bmap)
        imap = self.induced_map(edge[0],bmap)
        map0 = [{} for _ in range(len(imap))]
        map1 = [{} for _ in range(len(imap))]
        for i in range(len(imap)):
            for key in imap[i]:
                if (key>>d[0])%2 == 0:
                    sign = 1
                    for j in range(d[0]):
                        if (key>>j)%2 == 1:
                            sign *= -1
                    map0[i][key+(1<<d[0])] = -1*sign*imap[i][key]
                #else:
                #    map0[i][0] = 0
                if (key>>d[1])%2 == 0:
                    sign = 1
                    for j in range(d[1]):
                        if (key>>j)%2 == 1:
                            sign *= -1
                    #print("i,key,imap[i],d[1]")
                    #print(i,key,imap[i],d[1])
                    #print("map1[i][j+(1<<d[1])] = -1*imap[i][key]")
                    map1[i][key+(1<<d[1])] = sign*imap[i][key]
                #else:
                #    map1[i][0] = 0
        #print(imap)
        #print(d)
        #print(map0,map1)
        self.add_algebra_maps(map0,map1)
        return map0
        
    def old_edge_to_split_map(self,edge):
        assert (len(edge) == self.d)
        delta = edge.index(-1)
        source = []
        target = []
        source[:] = edge[:]
        target[:] = edge[:]
        source[delta] = 0
        target[delta] = 1
        source_res = self.get_res(source)
        target_res = self.get_res(target)
        a = self.get_x_diagram()[delta][0:2]
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
    
    #old imap {(i,j):k} where i --k--> j
    #old map0 = { (i,    j+(1<<d[0])) if (j>>d[0])%2 == 0) else 0    ):
    #              -1*k if (j>>d[0])%2 == 0 else 0  }
    #new imap [(j,k)] = [(t[0],t[1])]
    #new map0 = [ {t[0]+1<<d[0] if (t[0]>>d[0])%2 == 0 else 0 :
    #              -1*t[1] if (t[0]>>d[0])%2 == 0 else 0 }  for t in imap]
    # NEW VIE  = [(t[0]+1<<d[0],-1*t[1]) if (t[0]>>d[0])%2 == 0, else (0,0)]

    
    def get_circle(self,re,segment):
        resolution = self.get_res(re)
        for i in range(len(resolution)):
            if segment in resolution[i]:
                return i
        #Raise exception if here

    def square_sign(self,v00,a1,a2):
        #Algorithm due to Shumakovitch
        #Adapted to current set up
        #The input should be a vertex and the indices of two arrow changing from 0 to 1
        assert (self.v(v00,a1) == 0 and self.v(v00,a2) == 0 and a1 != a2)
        c00 = self.get_res_len(v00)
        c11 = self.get_res_len(v00 + (1<<a1) + (1<<a2))
        if c00 < c11:
            return -1 #Type A
        elif c11 < c00:
            return 1  #Type C
        pd = self.get_x_diagram()
        t1 = pd[a1][2]
        h1 = pd[a1][0]
        t2 = pd[a2][2]
        h2 = pd[a2][0]
        r00 = self.get_res(v00)
        lct1 = self.get_circle(v00,t1)
        lct2 = self.get_circle(v00,t2)

        if h1 not in r00[lct2] and h2 not in r00[lct1]:
            return -1 #Type A
        if h1 not in r00[lct1] or h2 not in r00[lct2]:
            return 1  #Type C
        assert(lct1 == lct2)
        assert(h1 in r00[lct1] and h2 in r00[lct2])
        #r01 = self.get_res(v00 + (1<<a2))
        rct1 = self.get_circle(v00 + (1<<a2),t1)
        rct2 = self.get_circle(v00 + (1<<a2),t2)
        if rct1 == rct2:
            return 1  #Type Y
        return -1     #Type X

    def edge_map(self,edge):
        assert (self.v(edge[0],edge[1]) == 0)
        scir = self.get_res_len(edge[0])
        next = edge[0] + (1<<edge[1])
        if scir < self.get_res_len(next):
            mp = self.split_map(edge)
            self.scalar_mult_algebra_map(self.get_edge_sign(edge),mp)
            return mp
        else:
            vmap = [(self.get_circle(next,self.get_res(edge[0])[i][0]),1) for i in range(scir)]
            #print(vmap)
            mp = self.induced_map(edge[0],vmap)
            self.scalar_mult_algebra_map(self.get_edge_sign(edge),mp)
            return mp
    
    def tex_braid(self):
        print("\\begin{center}\\begin{tikzpicture}")
        i = 0
        for x in self.presentation:
            if x > 0:
                print("\t\\braidline{{{}}}{{{}}}{{0}}{{{}}}{{1}}".format(self.b,x,i))
            else: 
                print("\t\\braidline{{{}}}{{{}}}{{0}}{{{}}}{{0}}".format(self.b,-1*x,i))
            i += 1
        print("\end{tikzpicture}\end{center}")
            
    def print_map(self,algebra_map,edge):
        print("  The source resolution, index 0..")
        print("      ",self.get_res(edge[0]))
        print("  The target resolution, index 0..")
        print("      ",self.get_res(edge[0] + (1<<edge[1])),"\n")
        for i in range(len(algebra_map)):
            print("       ", self.str_v(i),"maps to")
            for key in algebra_map[i]:
                if algebra_map[i][key] == 0:
                    print("             ZERO")                    
                else:
                    print("            ",self.str_v(key),"with factor",algebra_map[i][key])
        print()
        
    def tex_map(self,algebra_map,edge):
        print("The source resolution:")
        print("\\[",self.get_res(edge[0]),"\\]\n")
        print("The target resolution:")
        print("\\[",self.get_res(edge[0] + (1<<edge[1])),"\\]\n")
        print("\\begin{center}\\begin{tabular}{|r|r|r|}\\hline")
        print("Source&Target&Factor\\\\\\hline")
        for i in range(len(algebra_map)):
            print(self.str_v(i),"&&\\\\")
            for key in algebra_map[i]:
                print("&",self.str_v(key),"&",algebra_map[i][key],"\\\\")
            print("\\hline")
        print("\\end{tabular}\\end{center}")
        
    def texplanation(self):
        print("\\section{Braid: $",self.presentation ,"$}\n")
        self.tex_braid()
        print("\\vspace{.5cm}\n")
        print("\\subsection*{The distinguished vertex}\n")
        print("Vertex:",self.str_v(self.get_inv_v()))
        print("\nThe distinguished vertex's resolution:")
        print("\\[",self.get_res(self.get_inv_v()),"\\]")
        print("\\subsection*{Vertices that map to the vertex}\n")
        edges_in = self.get_edges_into_inv()
        for edge in edges_in:
            print("Vertex:", self.str_v(edge[0]),"resolving index:",edge[1])
            print("\nMap:")
            self.tex_map(self.edge_map(edge),edge)
        if not edges_in:
            print("NONE")
        
        print("\\subsection*{Vertices that the vertex maps to}\n")
        edges_out = self.get_edges_out_of_inv()
        for edge in edges_out:
            print("Vertex:", self.str_v(edge[0]+(1<<edge[1])),"resolving index:",edge[1])
            print("\nMap:")
            self.tex_map(self.edge_map(edge),edge)
        if not edges_out:
            print("NONE")
        print("\\newpage\n")
        



