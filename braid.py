import sys
import fields as fe
import edgestruct as es
import squarestruct as ss
import algebra as alg
import column as col

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
        self.nmin = self.neg()
        self.d = len(self.presentation)
        self.x_diagram = []
        self.closure = []
        self.inv_r = -1
        self.inv_v = -1
        self.edges_out_of_inv = []
        self.edges_into_inv = []
        self.res = [[] for _ in range(1 << self.d)]
        self.res_len = [0 for _ in range(1 << self.d)]
        self.edge_signs = es.EdgeStruct(self.d,0)
        self.squares = ss.SquareStruct(self.d)
        self.maps = es.EdgeStruct(self.d)
        self.char = 0
        #print("Initialized")



    #Basic helper functions

    def comp_b(self,braid):
        i = 0
        for s in braid:
            c = alg.abs(s)
            if c > i:
                i = c
        return i+1

    def neg(self):
        n = 0
        for x in self.presentation:
            if x < 0:
                n += 1
        return n

    def v(self,v,i):
        return (v >> i) % 2

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
            self.inv_r = alg.height(v)
        return self.inv_v


    #Basic Math Checks

    def raw_group_dims(self):
        dims = [0 for i in range(self.d + 1)]
        for v in range(1<<self.d):
            r = alg.height(v)
            dims[r] += (1<<self.get_res_len(v))
        return dims

    def raw_euler_characteristic(self):
        dims = self.raw_group_dims()
        chi = 0
        for i in range(len(dims)):
            chi += dims[i] * (-1) ** i
        return chi * (-1) ** self.nmin

    def euler_characteristic(self):
        dims = self.quantum_group_dims()
        shift = 0
        if self.nmin % 2 == 1:
            shift = 1
            for key in dims[0]:
                dims[0][key] = dims[0][key] -1
        for i in range(1,len(dims)):
            for key in dims[i]:
                dims[0][key] = dims[0].get(key,0) + ((-1) ** (i+shift))*dims[i][key]
        return dims[0]

    def quantum_group_dims(self):
        dims = [{} for i in range(self.d + 1)]
        for v in range(1<<self.d):
            r = alg.height(v)
            dv = self.get_res_len(v)
            for k in range(dv+1):
                Q = dv - 2*k + self.d - 3*self.nmin + r
                dims[r][Q] = dims[r].get(Q,0) + alg.choose(dv,k)
        return dims 

    def dd(self,r):
        assert (r >= 0 and r < self.d - 1)
        if not self.maps.is_fully_changed():
            self.comp_maps()
        p = [[] for _ in range(1<<self.d)]
        #Possibly to be deprecated so field stuff is handled later.
        zero = 0 #fe.FE(0,self.char)
        for v,i in self.maps:
            if alg.height(v) == r:
                u = v + (1<<i)
                alg_size = 1<<self.get_res_len(v)
                if not p[v]:
                    p[v] = [{} for _ in range(alg_size)]
                for s in range(alg_size):
                    for t1 in self.maps[v,i][s]:
                        for j in range(self.d):
                            if not u & (1<<j):#(u>>j)%2 == 0:
                                if i<j:
                                    a = (i,j)
                                else:
                                    a = (j,i)
                                p[v][s][a] = p[v][s].get(a,{})
                                for t2 in self.maps[u,j][t1]:
                                    prod = self.maps[v,i][s][t1] * self.maps[u,j][t1][t2]
                                    p[v][s][a][t2] = p[v][s][a].get(t2,zero) + prod
        return p

    def check_dd(self):
        print("\n\nCHECKING $dd = 0$\n\n")
        #error checking
        flag = True
        #error_list = []
        error_dict = {}
        for r in range(self.d - 1):
            prod = self.dd(r)
            #print(prod)
            #print(prod,"\n")
            for v in range((1<<r)-1,(((1<<r)-1)<<(self.d-r))+1):
                if alg.height(v) == r:
                    for g in range(len(prod[v])):
                        for de in prod[v][g]:
                            for t in prod[v][g][de]:
                                #print("$"+str(prod[v][g][de][t])+"$, ",end="")
                                #print(v,g,de,t,"$",str(prod[v][g][de][t]),"$",end="")
                                if prod[v][g][de][t].n:
                                    flag = False
                                    #error_list.append((self.str_v(v),self.str_v(g),de,self.str_v(t)))
                                    error_dict[(v,de[0],de[1])] = error_dict.get((v,de[0],de[1]),[]) + [(g,t)]
                                    #assert False
                else:
                    assert (not prod[v])
        if not flag:
            print("ERROR!!!!! (vector,i,j) : [$g\\to t$],\n\n")
            for v,i,j in error_dict:
                s = ""
                for map in error_dict[(v,i,j)]:
                    s += self.str_v(map[0]) + "\\to" + self.str_v(map[1]) + ", "
                print("({},{},{}) : $[{}]$\n\n".format(self.str_v(v),i,j,s))


    #Check if Invariant is Zero

    def make_column_dictionary(self,r,column):
        #
        rev_dict = {}
        i = 0
        column(r)
        for col in column:
            rev_dict[col] = i
            i += 1
        return i, rev_dict

    def make_row_dictionary(self,r,row):
        rev_dict = {}
        i = 0
        if r <= self.d and r >= 0:
            row(r)
            for v,g in row:
                rev_dict[(v,g)] = i
                i += 1
        return i, rev_dict

    def make_entry_dictionary(self,r,verticals):
        rev_dict = {}
        i = 0
        if r <= self.d and r >= 0:
            verticals(r)
            for coord in verticals:
                rev_dict[coord] = i
                i += 1
        return i, rev_dict

    def make_matrix(self,r,column,row,augmented = False):
        n, row_i = self.make_row_dictionary(r+1,row)
        #print(n,row_i)
        m, col_i = self.make_column_dictionary(r,column)
        col_labels = [None for _ in range(m)]
        for key in col_i:
            col_labels[col_i[key]] = key
        #print(m,col_i)
        #Possibly to be deprecated 
        #M = [[fe.FE(0,self.char) for _ in range(m)] for _ in range(n)]
        M = [[0 for _ in range(m)] for _ in range(n)]
        column(r)
        for v,g in column:
            for i in range(self.d):
                if self.v(v,i) == 0:
                    u = v + (1<<i)
                    for h in self.maps[v,i][g]:
                        M[row_i[u,h]][col_i[v,g]] = self.maps[v,i][g][h]
        if augmented:
            col_labels.append((0,0))
            inv_v = self.get_inv_v()
            inv_g = (1<<self.get_res_len(inv_v))-1
            for row in row_i:
                if row == (inv_v, inv_g):
                    M[row_i[row]].append(1) #fe.FE(1,self.char)) #possibly to be deprecated
                else:
                    M[row_i[row]].append(0) #fe.FE(0,self.char)) #possibly to be deprecated
        return M, col_labels

    def make_matrix_pair(self,r,verticals):
        k, col_i = self.make_entry_dictionary(r-1,verticals)
        m, mid_i = self.make_entry_dictionary(r,verticals)
        n, row_i = self.make_entry_dictionary(r+1,verticals)
        # In this convention C_(r-1) --B--> C_(r) --A--> C_(r+1) so that we have AB = 0
        # A is an (n x m) matrix and B is an (m x k) matrix.
        A = [[0 for _ in range(m)] for _ in range(n)]
        B = [[0 for _ in range(k)] for _ in range(m)]
        if r > 0:
            verticals(r-1)
            for v,g in verticals:
                for i in range(self.d):
                    if self.v(v,i) == 0:
                        u = v + (1<<i)
                        for h in self.maps[v,i][g]:
                            B[mid_i[u,h]][col_i[v,g]] = self.maps[v,i][g][h]
        inv_v = self.get_inv_v()
        inv_g= (1<<self.get_res_len(inv_v))-1
        for mid in mid_i:
            if mid == (inv_v,inv_g):
                inv_index = mid_i[mid]
                B[mid_i[mid]].append(1)
            else:
                B[mid_i[mid]].append(0)
        if r < self.d:
            verticals(r)
            for v,g in verticals:
                for i in range(self.d):
                    if self.v(v,i) == 0:
                        u = v + (1<<i)
                        for h in self.maps[v,i][g]:
                            A[row_i[u,h]][mid_i[v,g]] = self.maps[v,i][g][h]
        return A,B,k,m,n,inv_index

    def comp_inv(self):
        #prepare invariant vertex and height
        v = self.get_inv_v()
        r = self.inv_r

        if r == 0:
            print("\n\nThe invariant is non-zero and in the lowest level homology group.\n")
            return 0
        #prepare maps
        self.comp_maps()
        for i in range(1<<self.d):
            if not self.res_len[v]:
                assert False, "We should not get here because res_len should be set by comp_maps."
                self.get_res_len(v)
        #prepare matrices
        verticals = col.Row(self.d,self.res_len)
        A,B,k,m,n,inv_index = self.make_matrix_pair(r,verticals)
        TempB = [[B[i][j] for j in range(k+1)] for i in range(m)]
        in_image,factor = alg.i_row_reduce(TempB)
        if in_image:
            print("\n\nThe invariant is in the image of $d$ over $\\mathbb Z$, and is thus is $0$ over $\\mathbb Z,\\mathbb Z/p,$ and $\\mathbb Q$.\n")
            return 0
        for char in [0,2,3,5]:
            TempA = [[fe.FE(A[i][j],char) for j in range(m)] for i in range(n)]
            TempB = [[fe.FE(B[i][j],char) for j in range(k)] for i in range(m)]
            S,psi = alg.simultaneous_reduce(TempA,TempB,alg.print_mat,inv_index)
            assert len(S) == len(psi)
            #self.print_eq(S,psi)
            #for entry in psi:
            #    print("${}$ ".format(entry),end="")
            #if statements once it's returning something
        S,D,_ = alg.smith_normal_form(B)
        d = []
        i = 0
        while i < min(m,k) and D[i][i]:
            d.append(D[i][i])
            i += 1
        y = [1 if i == inv_index else 0 for i in range(m)]
        mult = alg.solve_mat_mult(S,d,y) 
        if mult == 0:
            print("\n\n There is no $n\\in\\mathbb Z$ and $x$ such that $dx=n\\psi(L)$. Thus $\\psi(L)$ is non-torsion.")
        elif mult == 1:
            print("\n\n There is an $x$ such that $dx=\\psi(L)$. Thus $\\psi(L)=0$.")
        else:
            print("\n\n The smallest positive integer $n$ such that there is an $x$ such that $dx=n\\psi(L)$ is {}. Thus $\\psi(L)$ is torsion.".format(mult))
        #compute simultaneous reduction of integer stuff




    def inv_nonzero(self):
        v = self.get_inv_v()
        r = self.inv_r
        self.comp_maps()
        if r == 0:
            print("\n\nThe Invariant is Non-ZERO: Zero Map\n")
            return 0 
        for i in range(1<<self.d):
            if not self.res_len[v]:
                self.get_res_len(v)
        column = col.Row(self.d,self.res_len)
        row = col.Row(self.d,self.res_len)
        M, _ = self.make_matrix(r-1,column,row,True)
        #alg.print_mat(M)
        if alg.dumb_row_reduce(M):
            print("\n\nThe Invariant is Non-ZERO: Pivot in last column\n")
        else:
            print("\n\nThe Invariant is ZERO: All pivots before last column\n")
        #alg.print_mat(M)
        print("\n\n")
        return 1

    def inv_factors(self):
        v = self.get_inv_v()
        r = self.inv_r
        self.comp_maps()
        if r == 0:
            print("\n\nThe invariant Non-Zero, and in the lowest level homology group.\n")
            return 0
        for i in range(1<<self.d):
            if not self.res_len[v]:
                self.get_res_len(v)
        column = col.Row(self.d,self.res_len)
        row = col.Row(self.d,self.res_len)
        M, col_labels = self.make_matrix(r-1,column,row,True)
        for char in [0,2,3,5]:
            #copy matrix
            TempM = [[fe.FE(M[i][j],char) for j in range(len(M[0]))] for i in range(len(M))]
            assert len(col_labels) == len(TempM[0]), "{},{},{}".format(self.presentation,len(col_labels),len(TempM[0]))
            in_image, factor = alg.dumber_row_reduce(TempM, col_labels,self.d)
            assert (isinstance(factor,int))
            if not in_image:
                if char == 0:
                    print("\n\nThe invariant is non-zero over $\\mathbb Q$, {}\n".format(factor))
                else:
                    print("\n\nThe invariant is non-zero over $\\mathbb Z/{}$, {}\n".format(char,factor))
            """
            else:
                if char == 0:
                    print("\n\nThere is an element $\\alpha$ such that $d\\alpha=\\alpha=\\psi$ over $\\mathbb Q$.\n")
                else:
                    print("\n\nThere is an element $\\alpha$ such that $d\\alpha={}\\alpha=\\psi$ over $\\mathbb Z/{}$.\n".format(factor,char))
            """
        in_image, factor = alg.int_row_reduce(M,col_labels,self.d)
        if not in_image:
            print("\n\nThe invariant is non-zero over $\\mathbb Z$, {}\n".format(factor))
        """
        else:
            print("\n\nThere is an element $\\alpha$ such that $d\\alpha={}\\alpha=\\psi$ over $\\mathbb Z$.\n".format(factor))
        """


    #TeX output functions

    def str_q(self,a,p):
        if a == 0:
            return ""
        if p == 0:
            if a > 0:
                return " + " + str(a)
            return " - " + str(-a)
        if a == 1: 
            if p == 1:
                return " + q"
            return " + q^{" + str(p) + "}"
        if a == -1:
            if p == 1:
                return " - q"
            return " - q^{" + str(p) + "}"
        if a > 1:
            if p == 1:
                return " + " + str(a) + "q"
            return " + " + str(a) + "q^{" + str(p) + "}"
        if p == 1:
            return " - " + str(a) + "q"
        return " - " + str(-a) + "q^{" + str(p) + "}"
        #Raise exception if here

    def print_Jones(self):
        poly = self.euler_characteristic()
        s = ""
        for key in sorted(poly.keys()):
            s += self.str_q(poly[key],key)
        if s[1] == '-':
            return '-'+s[3:]
        return s[3:]

    def str_v(self,v):
        s = ""
        for i in range(self.d):
            s += str((v>>i)%2)
        return s

    def str_e(self,g,l):
        s = ""
        for i in range(l):
            s += str((g>>i)%2)
        return s

    def print_eq(self,S,psi):
        print("\\[\\left[\\begin{{array}}{{{}}}".format("r"*len(S)))
        for row in S:
            print("&".join(map(str,row)),"\\\\")
        print("\\end{array}\\right]x=")
        print("\\left[\\begin{array}{r}")
        for e in psi:
            print(e,"\\\\",end="")
        print("\\end{array}\\right]\\]")

    def tex_map(self,algebra_map,edge):
        sv = edge[0]
        sl = self.get_res_len(sv)
        tv = edge[0] + (1<<edge[1])
        tl = self.get_res_len(tv)
        print("The source resolution:")
        print("\\[",self.get_res(sv),"\\]\n")
        print("The target resolution:")
        print("\\[",self.get_res(tv),"\\]\n")
        print("\\begin{center}\\begin{tabular}{|r|r|r|}\\hline")
        print("Source&Target&Factor\\\\\\hline")
        for i in range(len(algebra_map)):
            print(self.str_e(i,sl),"&&\\\\")
            for key in algebra_map[i]:
                print("&",self.str_e(key,tl),"&",algebra_map[i][key],"\\\\")
            print("\\hline")
        print("\\end{tabular}\\end{center}")

    def texplanation(self,s = "", maps = True, ddCheck = False):
        print("\\section{Braid: ",s,"$",self.presentation ,"$}\n")
        self.tex_braid()
        print("\nUngraded Euler Characteristic: $",self.raw_euler_characteristic(),"$\n")
        print("\n(Unreduced) Jones Polynomial: $\\hat{J}(L) = $\n\\[",self.print_Jones(),"\\]\n")
        print("\\subsection*{The distinguished vertex}\n")
        print("Vertex:",self.str_v(self.get_inv_v()))
        print("\nThe distinguished vertex's resolution:")
        print("\\[",self.get_res(self.get_inv_v()),"\\]")
        print("The algebra element of the invariant is", "1"*self.b)
        if maps:
            #### Fix maps so correct number of 0s
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
        if ddCheck:
            self.check_dd()
        print("\n\\noindent ",end="")
        self.comp_inv()
        print("\\newpage\n")

    def tex_braid(self):
        print("\\begin{center}\\begin{tikzpicture}[scale = 0.4]")
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



    #Resolution Computation

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
                c = alg.abs(x)
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
        self.res_len[v] = len(self.res[v])
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

    def get_circle(self,re,segment):
        resolution = self.get_res(re)
        for i in range(len(resolution)):
            if segment in resolution[i]:
                return i
        #Raise exception if here



    #Map Computation

    def comp_maps(self):
        if self.maps.is_fully_changed():
            return None
        self.comp_edge_signs()
        for edge in self.maps:
            self.maps[edge] = self.edge_map(edge)

    def comp_edge_signs(self):
        ## Algorithm due to Shumakovitch
        #Make sure edge_signs is initialized
        if self.edge_signs.is_fully_changed():
            return self.edge_signs
        n = self.d
        for delta in range(self.d):
            for base in range(1<<delta):
                self.edge_signs[base,delta] = 1
        for v,i,j in self.squares:
            self.edge_signs[v + (1<<j),i] = -self.edge_signs[v,i] * self.square_sign(v,i,j)
        return self.edge_signs

    def get_edge_sign(self,edge):
        return self.comp_edge_signs()[edge[0],edge[1]]

    def square_sign(self,v,i,j):
        #Algorithm due to Shumakovitch
        #Adapted to current set up
        #The input should be a vertex and the indices of two arrow changing from 0 to 1
        assert (self.v(v,i) == 0 and self.v(v,j) == 0 and i != j)
        c00 = self.get_res_len(v)
        c11 = self.get_res_len(v + (1<<i) + (1<<j))
        if c00 < c11:
            self.squares[v,i,j] = "Type A, increasing cycles"
            return -1 #Type A
        elif c11 < c00:
            self.squares[v,i,j] = "Type C, decreasing cycles"
            return 1  #Type C
        pd = self.get_x_diagram()
        t1 = pd[i][2]
        h1 = pd[i][0]
        t2 = pd[j][2]
        h2 = pd[j][0]
        r00 = self.get_res(v)
        lct1 = self.get_circle(v,t1) #index of circle t1 is in
        lch2 = self.get_circle(v,h2) #index of circle h2 is in

        if (h1 in r00[lch2] and t2 in r00[lct1]) and (h1 not in r00[lct1]):
            self.squares[v,i,j] = "Type A, same no., arrows point likewise"
            return -1 #Type A
        if h1 not in r00[lct1] or t2 not in r00[lch2]:
            self.squares[v,i,j] = "Type C, same no., arrows point differently"
            return 1  #Type C
        assert(lct1 == lch2)
        assert(h1 in r00[lct1] and t2 in r00[lct1])
        #r01 = self.get_res(v00 + (1<<a2))
        rct1 = self.get_circle(v + (1<<j),t1)
        rct2 = self.get_circle(v + (1<<j),t2)
        if rct1 == rct2:
            self.squares[v,i,j] = "Type Y"
            return 1  #Type Y
        self.squares[v,i,j] = "Type X"
        return -1     #Type X

    def edge_map(self,edge):
        assert (self.v(edge[0],edge[1]) == 0)
        scir = self.get_res_len(edge[0])
        next = edge[0] + (1<<edge[1])
        if scir < self.get_res_len(next):
            mp = self.split_map(edge) #TAG FIELD UPDATE
            self.scalar_mult_algebra_map(self.get_edge_sign(edge),mp)
            return mp
        else:
            #Possibly to be deprecated
            #vmap = [(self.get_circle(next,self.get_res(edge[0])[i][0]),fe.FE(1,self.char)) for i in range(scir)]
            vmap = [(self.get_circle(next,self.get_res(edge[0])[i][0]),1) for i in range(scir)]
            #print(vmap)
            mp = self.induced_map(edge[0],vmap) #TAG FIELD UPDATE
            self.scalar_mult_algebra_map(self.get_edge_sign(edge),mp)
            return mp

    def scalar_mult_algebra_map(self,scalar,mp):
        for i in range(len(mp)):
            for key in mp[i]:
                mp[i][key] = mp[i][key] * scalar
        #not returned, stored in place in mp

    def add_algebra_maps(self,map1,map2):
        assert (len(map1) == len(map2))
        #Possibly to be deprecated so field stuff is handled later
        zero = 0 #fe.FE(0,self.char)
        for i in range(len(map2)):
            for key in map2[i]:
                map1[i][key] = map1[i].get(key,0) + map2[i][key]
        #not returned, stored in place in map1

    def induced_map(self,v,vmap):
        xd = self.get_res_len(v)
        assert (xd == len(vmap))
        #Possibly to be deprecated so field stuff is handled later
        zero = 0 #fe.FE(0,self.char)
        #print(vmap)
        algebra_map = [{} for _ in range(1<<xd)]
        for i in range(1<<xd):
            ls = []
            #Possibly to be deprecated so field stuff is handled later
            p = 1 #fe.FE(1,self.char)
            for j in range(xd):
                if (i>>j)%2 == 1:
                    ls.append(vmap[j][0])
                    p = p * vmap[j][1]
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
            if sign:
                algebra_map[i][target] = algebra_map[i].get(target,zero) + p*sign
        return algebra_map

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
                        bmap[i] = (j,1) #fe.FE(1,self.char)) #possibly to be deprecated
                #if involved in the split
                else: 
                    #arrow head in target ?
                    if a[0] in target_res[j]:
                        bmap[i] = (j,1) #fe.FE(1,self.char)) #poassibly to be deprecated
                        assert(d[0] == -1)
                        d[0] = j
        #print(bmap)
        imap = self.induced_map(edge[0],bmap) #TAG for UPDATE
        map0 = [{} for _ in range(len(imap))]
        map1 = [{} for _ in range(len(imap))]
        for i in range(len(imap)):
            for key in imap[i]:
                if (key>>d[0])%2 == 0:
                    sign = 1
                    for j in range(d[0]):
                        if (key>>j)%2 == 1:
                            sign *= -1
                    map0[i][key+(1<<d[0])] = -imap[i][key]*sign
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
                    map1[i][key+(1<<d[1])] = imap[i][key]*sign
                #else:
                #    map1[i][0] = 0
        #print(imap)
        #print(d)
        #print(map0,map1)
        self.add_algebra_maps(map0,map1)
        return map0



    #Possibly Defunct Functions 
    """ 
    def get_edges_out_of_inv(self):
        if not self.edges_out_of_inv:
            self.edges_out_of_inv = self.edges_out(self.get_inv_v())
        return self.edges_out_of_inv

    def get_edges_into_inv(self):
        if not self.edges_into_inv:
            self.edges_into_inv = self.edges_in(self.get_inv_v())
        return self.edges_into_inv

    def edges(self):
        #Check if edges_signs has already be initilized or computed
        if self.edge_signs.is_fully_changed():
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

    def height(self,v):
        r = 0 
        while v:
            r += v%2
            v = v >> 1
        return r

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
    """




def make_braid_list(file):
    ls = []
    with open(file) as f:
        for line in f:
            sl = line.strip().split()
            ls.append([int(i) for i in sl])
    return ls

def make_fancy_braid_list(file):
    ls = []
    with open(file) as f:
        for line in f:
            sl = line.strip().split()
            ls.append((sl[0],[int(sl[i]) for i in range(1,len(sl))]))
    return ls

def driver(fancy = False):
    if fancy:
        braid_list = make_fancy_braid_list(sys.argv[1])
    else:
        braid_list = make_braid_list(sys.argv[1])
    print("""\\documentclass{article}
\\usepackage{thesispkg}
\\usepackage{fancyhdr}

\\pagestyle{fancy}
\\fancyhf{}
\\renewcommand{\\sectionmark}[1]{\\markright{#1}}
\\lhead{\\fancyplain{}{\\rightmark }} 

\\begin{document}""")
    for braid in braid_list:
        if fancy:
            b = Braid(braid[1])
            b.texplanation(braid[0],False)
            #print(b.edge_signs)
            #for map in b.maps:
            #    print(map,b.maps[map],"\n\n")
            #for v,i,j in b.squares:
            #    print((b.str_v(v),i,j),b.squares[v,i,j],"\n\n")
            #print("\\newpage\n\n")
        else:
            b = Braid(braid)
            b.texplanation()
    print("\\end{document}")

def test_driver():
    braid_list = make_fancy_braid_list(sys.argv[1])
    for braid in braid_list:
        b = Braid(braid[1])
        #print(braid[0], b.raw_euler_characteristic())
        print(braid[0], b.print_Jones())

if __name__ == "__main__":
    #test_driver()
    driver(True)