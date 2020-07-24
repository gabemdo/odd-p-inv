import sys
import braid as br
import fields as fe
import edgestruct as es
import squarestruct as ss
import algebra as alg
import column as col

class KhovanovHomology:
    # Some convetions:
    # Implement x used for crossing index number in the braid word
    # x ranges in self.d, the number of crossing or the dimention of the hypercube
    # Implement s used for crossing strand number in the braid word
    # s ranges in self.b, the braid index or the number strands in the braid
    #           | self.word[x] | = s
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

    def __init__(self,braid_word):
        #The braid word of a braid is given as a list of non-zero integers
        if True: #TODO check if braid_word
            assert (len(braid_word) > 0)
            self.braid = br.Braid(braid_word)
            self.d = self.braid.get_crossing_number()
            self.res = [[] for _ in range(1 << self.d)]
            self.res_len = [0 for _ in range(1 << self.d)]
            for i in range(1 << self.d):
                self.res[i] = self.braid.get_res(i)
                self.res_len[i] = len(self.res[i])
                print(self.res_len[i])
            self.inv_v = self.braid.get_inv_v()
            self.inv_r = self.braid.get_inv_r()
            self.verticals = col.Row(self.d,self.res_len)
        if False:
            self.transverse = False
            pass #TODO check if other type of input, e.g. grid or just x_diagram
        self.edge_signs = es.EdgeStruct(self.d,0)
        self.squares = ss.SquareStruct(self.d)
        #TODO: split by grading
        self.maps = es.EdgeStruct(self.d)
        self.computed = False
        #self.maps = {}
        #self.char = 0
        self.transverse = True #TODO, implement doing all of the non invariant stuff from planar diagrams more generally
        #TODO, also implement doing all of this from initializing with a grid or a braid directly. 
        self.odd = True
        #TODO: add in storage for grading information
        #Note: Number of components (hence dim V_alpha) is at most self.d, and thus between 1 and self.d
        #Note: Exterior algebra alg-grading, k, is limited by dim V_alpha, and thus between 0 and dim V_alpha
        #Note: Thus Q_0 = (dim V_alpha) - 2k is     -self.d <= Q_0 <= self.d
        #Note: Shift S = n+ - 2* n- = self.d - 3 * self.nmin
        #TODO: After testing compare to M_0 below
        #Note: M_0 = "r" and thus between 0 and self.d 
        #Note: And Q = Q_0 + S + M_0
        #Hence: -self.d <= Q <= 2*self.d
        #Storage convention: store in dictionary by key Q_0+M_0
        #print("Initialized")

    #TODO
    def set_odd(self):
        if not self.odd:
            self.reset()
        self.odd = True

    def set_even(self):
        if self.odd:
            self.reset()
        self.odd = False

    def reset(self):
        self.maps = es.EdgeStruct(self.d)
        self.computed = False

    #Basic helper functions
    def target_vertex(self,edge):
        assert self.is_edge(edge)
        return edge[0] + (1<<edge[1])

    def is_edge(self,edge):
        if alg.v(*edge):
            return False
        return True

    def get_res(self,v):
        return self.res[v]

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
        return chi * (-1) ** self.braid.get_nmin()

    def euler_characteristic(self):
        dims = self.quantum_group_dims()
        shift = 0
        if self.braid.get_nmin() % 2 == 1:
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
                Q = dv - 2*k + self.d - 3*self.braid.get_nmin() + r
                dims[r][Q] = dims[r].get(Q,0) + alg.choose(dv,k)
        return dims 

    def dd(self,r):
        assert (r >= 0 and r < self.d - 1)
        if not self.computed:#self.maps.is_fully_changed():
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
                                if prod[v][g][de][t]:#.n:
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

    def make_entry_dictionary(self,r):
        rev_dict = {}
        i = 0
        if r <= self.d and r >= 0:
            self.verticals(r)
            for coord in self.verticals:
                rev_dict[coord] = i
                i += 1
        return i, rev_dict

    def grade(self,v,g):
        dimV = self.get_res_len(v)
        k = alg.height(g)
        r = alg.height(v)
        #print("grading of {},{}, dimV= {}, k={},r={}. Grading = {}".format(v,g,dimV,k,r,dimV+r-2*k))
        if self.odd:
            return dimV + r - 2*k
        return 2*k - dimV + r

    def make_graded_entry_dictionary(self,r,grade):
        rev_dict = {}
        i = 0
        if r <= self.d and r >= 0:
            self.verticals(r)
            for v,g in self.verticals:
                #print("v,g={},{}, r={}, grade={} =?= {}".format(v,g,r,grade,self.grade(v,g)))
                if grade == self.grade(v,g):
                    rev_dict[v,g] = i
                    i += 1
        return i, rev_dict

    def make_graded_entry_dictionary_set(self):
        rev_dict = [{} for _ in range(self.d+1)]
        i = [{} for _ in range(self.d+2)] #possibly + 2 for terminal 0s
        for r in range(self.d+1):
            self.verticals(r)
            for v,g in self.verticals:
                grade = self.grade(v,g)
                rev_dict[r][grade] = rev_dict[r].get(grade,{})
                i[r][grade] = i[r].get(grade,0)
                rev_dict[r][grade][v,g] = i[r][grade]
                i[r][grade] += 1
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
                if alg.v(v,i) == 0:
                    u = v + (1<<i)
                    for h in self.maps[v,i][g]:
                        M[row_i[u,h]][col_i[v,g]] = self.maps[v,i][g][h]
        if augmented:
            col_labels.append((0,0))
            inv_v = self.inv_v
            inv_g = (1<<self.get_res_len(inv_v))-1
            for row in row_i:
                if row == (inv_v, inv_g):
                    M[row_i[row]].append(1) #fe.FE(1,self.char)) #possibly to be deprecated
                else:
                    M[row_i[row]].append(0) #fe.FE(0,self.char)) #possibly to be deprecated
        return M, col_labels

    def make_graded_matrix_pair(self,r,grade):
        k, col_i = self.make_graded_entry_dictionary(r-1,grade)
        m, mid_i = self.make_graded_entry_dictionary(r,grade)
        n, row_i = self.make_graded_entry_dictionary(r+1,grade)
        #Later: turn A to non-field option?
        A = []
        B = []
        if r>0:
            B = [[0 for _ in range(k)] for _ in range(m)]
            self.verticals(r-1)
            for v,g in self.verticals:
                if self.grade(v,g) == grade:
                    for i in range(self.d):
                        if alg.v(v,i) == 0:
                            u = v + (1<<i)
                            for h in self.maps[v,i][g]:
                                B[mid_i[u,h]][col_i[v,g]] = self.maps[v,i][g][h]
        inv_v = self.inv_v
        inv_g = 0
        if self.odd:
            inv_g = (1<<self.get_res_len(inv_v))-1
        #print(k,m,n,inv_v,inv_g,mid_i)
        #print(mid_i)
        inv_index = mid_i[inv_v,inv_g]
        if r < self.d:
            A = [[fe.FE(0,0) for _ in range(m)] for _ in range(n)]
            self.verticals(r)
            for v,g in self.verticals:
                if self.grade(v,g) == grade:
                    for i in range(self.d):
                        if alg.v(v,i) == 0:
                            u = v + (1<<i)
                            for h in self.maps[v,i][g]:
                                value = self.maps[v,i][g][h]
                                A[row_i[u,h]][mid_i[v,g]] = fe.FE(value,0)
        return A,B,k,m,n,inv_index

    def make_matrix_pair(self,r,augmented = True,field_A = False):
        k, col_i = self.make_entry_dictionary(r-1)
        m, mid_i = self.make_entry_dictionary(r)
        n, row_i = self.make_entry_dictionary(r+1)
        # In this convention C_(r-1) --B--> C_(r) --A--> C_(r+1) so that we have AB = 0
        # A is an (n x m) matrix and B is an (m x k) matrix.
        if field_A:
            A = [[fe.FE(0,0) for _ in range(m)] for _ in range(n)]
        else:
            A = [[0 for _ in range(m)] for _ in range(n)]
        B = [[0 for _ in range(k)] for _ in range(m)]
        if r > 0:
            self.verticals(r-1)
            for v,g in self.verticals:
                for i in range(self.d):
                    if alg.v(v,i) == 0:
                        u = v + (1<<i)
                        for h in self.maps[v,i][g]:
                            B[mid_i[u,h]][col_i[v,g]] = self.maps[v,i][g][h]
        inv_v = self.inv_v
        inv_g= (1<<self.get_res_len(inv_v))-1
        inv_index = mid_i[inv_v,inv_g]
        if augmented:
            for mid in mid_i:
                if mid == (inv_v,inv_g):
                    B[mid_i[mid]].append(1)
                else:
                    B[mid_i[mid]].append(0)
        if r < self.d:
            self.verticals(r)
            for v,g in self.verticals:
                for i in range(self.d):
                    if alg.v(v,i) == 0:
                        u = v + (1<<i)
                        for h in self.maps[v,i][g]:
                            value = self.maps[v,i][g][h]
                            if field_A:
                                A[row_i[u,h]][mid_i[v,g]] = fe.FE(value,0)
                            else:
                                A[row_i[u,h]][mid_i[v,g]] = value
        return A,B,k,m,n,inv_index

    def make_graded_matrix(self,r,grade):
        m, col_i = self.make_graded_entry_dictionary(r,grade)
        n, row_i = self.make_graded_entry_dictionary(r+1,grade)
        A = [[0 for _ in range(m)] for _ in range(n)]
        self.verticals(r)
        for v,g in self.verticals:
            if self.grade(v,g) == grade:
                for i in range(self.d):
                    if alg.v(v,i) == 0:
                        u = v + (1<<i)
                        for h in self.maps[v,i][g]:
                            A[row_i[u,h]][col_i[v,g]] = self.maps[v,i][g][h]
        return A

    def make_graded_matrix_set(self):
        matrix_set = [{} for _ in range(self.d+2)]
        dim, index = self.make_graded_entry_dictionary_set()
        #for grade in indices[0]:
        #    matrix_set[-1][grade] = []
        #    indices[-1][grade] = 0
        for r in range(self.d):
            #NOTE spot of r+1 is rows
            #NOTE spot of  r  is cols
            for grade in dim[r]:
                matrix_set[r][grade] = [[0 for _ in range(dim[r][grade])] for _ in range(dim[r+1].get(grade,0))]
                #indices[-1][grade] = 0 
            self.verticals(r)
            for v,g in self.verticals:
                grade = self.grade(v,g)
                for i in range(self.d):
                    if alg.v(v,i) == 0:
                        u = v + (1<<i)
                        for h in self.maps[v,i][g]:
                            value = self.maps[v,i][g][h]
                            matrix_set[r][grade][index[r+1][grade][u,h]][index[r][grade][v,g]] = value
        #for grade in indices[-1]:
        #    print("Grade = {:>2},  {:>2}".format(grade, " ".join([" {:>2}".format(indices[r].get(grade,0)) for r in range(self.d+1)])))
        return dim, matrix_set

    def make_matrix(self,r):
        m, col_i = self.make_entry_dictionary(r)
        n, row_i = self.make_entry_dictionary(r+1)
        A = [[0 for _ in range(m)] for _ in range(n)]
        self.verticals(r)
        for v,g in self.verticals:
            for i in range(self.d):
                if alg.v(v,i) == 0:
                    u = v + (1<<i)
                    for h in self.maps[v,i][g]:
                        A[row_i[u,h]][col_i[v,g]] = self.maps[v,i][g][h]
        return A

    def comp_full_graded_homology(self):
        self.comp_maps()
        #self.check_dd()
        #Grading and Thin-ness stuff
        if self.odd:
            group_str = "KH'"
        else:
            group_str = "KH"

        shift = self.d - 3*self.braid.get_nmin()
        print("Q-Grade is shifting up by {}. H-Grade is shifted down by {}.".format(shift,self.braid.get_nmin()))
        diagonal = {}
        sl = self.d - (2*self.braid.get_nmin() + self.braid.get_b())   #the self-linking number
        #Set up
        self.verticals = col.Row(self.d,self.res_len)
        #print("Odd?",self.odd)
        #print(self.maps)
        indices, matrix_set = self.make_graded_matrix_set()
        ##print(indices)
        started = False
        ended = False
        zero_start = -1
        for r in range(self.d+1):
            ##print("RRRR:", r)
            s = []
            for grade in indices[r]:
                if r > 0 and grade in matrix_set[r-1]:
                    B = [row[:] for row in matrix_set[r-1][grade]]
                else:
                    B = []
                if r < self.d and grade in matrix_set[r]:
                    A = [row[:] for row in matrix_set[r].get(grade,[])]
                else:
                    A = []
                #TEST to deal with m8_19 anomoly 
                #if grade == 5 and r == 5:
                #    print(B)
                #    print("\n",A)
                ##print("\nHeight: {}; Grading: {}".format(r,grade))
                ##print(B,"\n")
                ##print(A)
                k = 0 if r == 0 else indices[r-1].get(grade,0)
                m = indices[r].get(grade,0)
                n = 0 if r == self.d else indices[r+1].get(grade,0) 
                ##print("k= {}, m={}, n= {}".format(k,m,n))
                #print("r = {}, grade = {}, (k,m,n) = ({},{},{})".format(r,grade,k,m,n))
                betti,tor = alg.integer_homology(A,B,k,m,n)
                #print("\nBetti: {}; Torsion: {}".format(betti,tor))
                #if tor and tor[0] == 3:
                #    print(grade,r)
                #    print(B,"\n")
                #    print(A,"\n")
                #    print(k,m,n)
                groups = (["Z^{}".format(betti)] if betti else []) + ["Z/{}".format(p) for p in tor]
                #groups += ["Z/{}".format(p) for p in tor]
                if groups:
                    d = (grade+shift) - 2*(r- self.braid.get_nmin())
                    diagonal[d] = diagonal.get(d,0) + 1
                    if len(groups) == 1:
                        s.append("{}[{:>2}]".format(groups[0],grade+shift))
                    else:
                        s.append("({})[{:>2}]".format(" + ".join(groups),grade+shift))
            if not s:
                if started and not ended:
                    ended = True
                    zero_start = r
                continue
            elif ended:
                ended = False
                for i in range(zero_start,r):
                    print("{}_({:>2})(L) = 0".format(group_str,i-self.braid.get_nmin()))
                zero_start = -1
            else:
                started = True
            print("{}_({:>2})(L) = {}".format(group_str,r-self.braid.get_nmin()," + ".join(s)))
        if len(diagonal) not in [0,2]:
            max_width = 0
            mode_width = -1
            for width in diagonal:
                if diagonal[width] > max_width:
                    max_width = diagonal[width]
                    mode_width = width
            if diagonal.get(mode_width+2,0) > diagonal.get(mode_width-2,0):
                sigma = mode_width+1
            elif diagonal.get(mode_width+2,0) < diagonal.get(mode_width-2,0):
                sigma = mode_width-1
            else:
                sigma = "No guess."
            print("Wide knot, sigma = {}, sl = {}.\n".format(sigma, sl))
        else:
            sigma = 0
            for width in diagonal:
                sigma += width
            sigma //= 2
            print("Thin knot, sigma = {}, sl = {}.\n".format(sigma, sl))

    def comp_homology(self,char = -1,r = -1,tex = True):
        #prepare maps
        self.comp_maps()
        for i in range(1<<self.d):
            if not self.res_len[i]:
                assert False, "We should not get here because res_len should be set by comp_maps."
                self.get_res_len(i)
        self.verticals = col.Row(self.d,self.res_len)
        if r == -1:
            chain_maps = [None for _ in range(self.d)]
            for r in range(self.d):
                chain_maps[r] = self.make_matrix(r,self.verticals)
                #alg.print_mat_(chain_maps[r])
                #print()
                ##### STUFF
            #this adds the map from r=-1 and r=self.d the map []
            chain_maps.append([])
            for r in range(self.d+1):
                if char != 1:
                    B = []
                    A = []
                    if chain_maps[r-1]:
                        B = [[fe.FE(col,0) for col in row] for row in chain_maps[r-1]]
                    if chain_maps[r]:
                        A = [[fe.FE(col,0) for col in row] for row in chain_maps[r]]
                    H = alg.field_homology(A,B)
                    if tex:
                        if H:
                            if char == 0:
                                print("\n$KH_{{{}}}'(L;\\mathbb Q)=\\mathbb Q^{{{}}}$".format(r,H))
                            else:
                                print("\n$KH_{{{0}}}'(L;\\mathbb Z/{1})=(\\mathbb Z/{1})^{{{2}}}$".format(r,char,H))
                        else:
                            if char == 0:
                                print("\n$KH_{{{}}}'(L;\\mathbb Q)=0$".format(r))
                            else:
                                print("\n$KH_{{{}}}'(L;\\mathbb Z/{}=0".format(r,char))
                    else:
                        if H:
                            if char == 0:
                                print("KH'_{}(L,Q) = Q^{}".format(r,H))
                            else:
                                print("KH'_{0}(L,Z/{1}) = (Z/{1})^{2}".format(r,char,H))
                        else:
                            if char == 0:
                                print("KH'_{}(L,Q) = 0".format(r))
                            else:
                                print("KH'_{}(L,Z/{}) = 0".format(r,char))
                if abs(char) == 1:
                    B = []
                    A = []
                    if chain_maps[r-1]:
                        B = [row[:] for row in chain_maps[r-1]]
                    if chain_maps[r]:
                        A = [row[:] for row in chain_maps[r]]
                    betti,tor = alg.integer_homology(A,B)
                    if tex:
                        torsion = ["\\mathbb Z^{{{}}}".format(betti)] if betti else [] + ["\\mathbb Z/{}".format(p) for p in tor]
                        print("\n$KH_{{{}}}'(L;\\mathbb Z)={}$".format(r," \\oplus ".join(torsion)))
                    else:
                        torsion = ["Z^{}".format(betti)] if betti else [] + ["Z/{}".format(p) for p in tor]
                        if not torsion:
                            torsion = ["0"]
                        print("KH'_{}(L) = {}".format(r," + ".join(torsion)))
            return 0
        elif r == 0:
            B_ = []
            A_ = self.make_matrix(r)
        elif r == self.d:
            B_ = self.make_matrix(r-1)
            A_ = []
        elif 0 < r < self.d:
            B_ = self.make_matrix(r-1)
            A_ = self.make_matrix(r)
        else:
            assert False, "Your value of r is out of range. The raw grading must be between 0 and the number of crossings."
        if char != 1:
            B = []
            if B_:
                B = [[fe.FE(col,0) for col in row] for row in B_]
            A = []
            if A_:
                A = [[fe.FE(col,0) for col in row] for row in A_]
            H = alg.field_homology(A,B)
            if tex:
                if char == 0:
                    print("\n$KH_{{{}}}'(L;\\mathbb Q)=\\mathbb Q^{{{}}}$".format(r,H))
                else:
                    print("\n$KH_{{{}}}'(L;\\mathbb Z/{})=(\\mathbb Z/{})^{{{}}}$".format(r,H))
            else:
                if char == 0:
                    print("KH'_{}(L;Q) = {}".format(r,"Q^"+str(H) if H else "0"))
                else:
                    print("KH'_{}(L;Z/{}) = {}".format(r,char,"(Z/{})^{}".format(char,H) if H else "0"))
        if abs(char) == 1:
            B = []
            if B_:
                B = [row[:] for row in B_]
            A = []
            if A_:
                A = [row[:] for row in A_]
            betti,tor = alg.integer_homology(A,B)
            if tex:
                torsion = ["\\mathbb Z^{}".format(betti)] if betti else [] + ["\\mathbb Z/{}".format(p) for p in tor]
                print("\n$KH_{{{}}}'(L;\\mathbb Z)={}$".format(r," \\oplus ".join(torsion)))
            else:
                torsion = ["Z^p".format(betti)] if betti else [] + ["Z/{}".format(p) for p in tor]
                if not torsion:
                    torsion = ["0"]
                print("KH'_{}(L) = {}".format(r," + ".join(torsion)))

    def comp_inv(self,print_output = True):
        v = self.inv_v
        r = self.get_res_len(v)
        g = 0
        if self.odd:
            g = (1<<r)-1
        r = self.inv_r

        print(v,r,g)
        if self.odd:
            group_str = "KH'"
        else:
            group_str = "KH"

        #print(v,g)
        grade = self.grade(v,g)
        #make a comp_maps that ignores extra levels
        self.comp_maps()
        for i in range(1<<self.d):
            if not self.res_len[i]:
                assert False, "We should not get here because res_len should be set by comp_maps."
                self.get_res_len(i)
        self.verticals = col.Row(self.d,self.res_len)
        A_Q,B,k,m,n,inv_index = self.make_graded_matrix_pair(r,grade)
        #print(k,m,n,B)
        if m and k:
            assert len(B) == m
            assert len(B[0]) == k
        assert len(A_Q) == n
        if n:
            assert len(A_Q[0]) == m
        if not (k and m):
            B = []
        if not (m and n):
            A = []
        d = []
        tor = []
        im = torsion = rank_A = 0
        S = []
        if k and m:
            S,D,_ = alg.smith_normal_form(B)
            while im < min(m,k) and D[im][im]:
                d.append(D[im][im])
                if abs(D[im][im]) != 1:
                    tor.append(D[im][im])
                    torsion += 1
                im += 1
        nontor = im - torsion
        if m and n:
            rank_A = alg.field_row_echelon(A_Q)
        ker = m - rank_A
        betti = ker - im
        if print_output:
            print("\n${{{}}}_0'(L;\\mathbb Q)={}$".format(group_str,"\\mathbb Q^{{{}}}".format(ker-im) if (ker-im) else "0"))
        betti_str = ["\\mathbb Z^{{{}}}".format(betti)] if betti else []
        homology_groups = betti_str + ["\\mathbb Z/{}".format(p) for p in tor]
        if not homology_groups:
            homology_groups = ["0"]
        if print_output:
            print("\n${{{}}}_0(L;\\mathbb Z)={}$".format(group_str," \\oplus ".join(homology_groups)))
            if self.inv_r == 0:
                print("\nThe invariant is non-zero, non-torsion, indivisible, and in the lowest level homology group.\n")
        mult = 0
        if k and m:
            y = [1 if i == inv_index else 0 for i in range(m)]
            mult = alg.solve_mat_mult(S,d,y) 
        if print_output:
            if mult == 0:
                print("\nThere is no $x$ and no $n\\in\\mathbb Z$ such that $dx=n\\psi(L)$. ")
                print("Thus $\\psi(L)$ is non-zero and non-torsion in $\\mathbb Z$ and $\\mathbb Q$. ")
            elif mult == 1:
                print("\nIn the chain complex $\\psi(L)$ is a boundary. Thus, in homology $\\psi(L)=0$ in all coefficients.")
            else:
                print("\nThe smallest positive integer $n$ such that there is an $x$ such that $dx=n\\psi(L)$ is {}. Thus $\\psi(L)$ is torsion.".format(mult))
        rational_homology = ker - im 
        integer_homology = (betti,tor)
        return mult, rational_homology, integer_homology


    def inv_nonzero(self):
        v = self.inv_v
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

    #TODO: deprecated??
    def inv_factors(self):
        v = self.inv_v
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
            assert len(col_labels) == len(TempM[0]), "{},{},{}".format(self.word,len(col_labels),len(TempM[0]))
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

    def texplanation(self,s = "", ddCheck = False):
        print("\\section{{Braid: {}, ${}$}}\n".format(s,self.braid.word))
        self.tex_braid()
        print("\nUngraded Euler Characteristic: ${}$\n".format(self.raw_euler_characteristic()))
        print("\n(Unreduced) Jones Polynomial: $\\hat J(L) = $\n\\[{}\\]\n".format(self.print_Jones()))
        print("\n\n\n\\subsection*{The distinguished vertex}\n")
        print("Vertex:",self.str_v(self.inv_v))
        print("\nThe distinguished vertex's resolution:")
        print("\\[{}\\]".format(self.get_res(self.inv_v)))
        #print("The algebra element of the invariant is", "1"*self.b)
        if ddCheck:
            self.check_dd()
        print("\n\n\n\\subsection*{The Invariant and its Homology Group}\n")
        self.comp_inv()
        print("\\newpage\n")

    def tex_braid(self):
        print("\\begin{center}\\begin{tikzpicture}[scale = 0.4]")
        i = 0
        for x in self.braid.word:
            if x > 0:
                print("\t\\braidline{{{}}}{{{}}}{{0}}{{{}}}{{1}}".format(self.braid.get_b(),x,i))
            else: 
                print("\t\\braidline{{{}}}{{{}}}{{0}}{{{}}}{{0}}".format(self.braid.get_b(),-1*x,i))
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
    def get_res_len(self,v):
        if isinstance(v,int):
            assert (v < (1<<self.d))
            return self.res_len[v]
        assert self.is_edge(v)
        v_len = self.get_res_len(v[0])
        u_len = self.get_res_len(self.target_vertex(v))
        return (v_len,u_len)

    def get_circle(self,re,segment):
        resolution = self.get_res(re)
        for i in range(len(resolution)):
            if segment in resolution[i]:
                return i
        #Raise exception if here



    #Map Computation

    def comp_maps(self):
        if self.computed: #self.maps.is_fully_changed():
            return None
        if self.odd:
            self.comp_edge_signs()
        for edge in self.maps:
            self.maps[edge] = self.edge_map(edge)
        self.computed = True

    def comp_edge_signs(self):
        ## Algorithm due to Shumakovitch
        #Make sure edge_signs is initialized
        if self.computed: #self.edge_signs.is_fully_changed():
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
        assert (alg.v(v,i) == 0 and alg.v(v,j) == 0 and i != j)
        c00 = self.get_res_len(v)
        c11 = self.get_res_len(v + (1<<i) + (1<<j))
        if c00 < c11:
            self.squares[v,i,j] = "Type A, increasing cycles"
            return -1 #Type A
        elif c11 < c00:
            self.squares[v,i,j] = "Type C, decreasing cycles"
            return 1  #Type C
        pd = self.braid.get_x_diagram()
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
        # v (edge[0]) and u (self.target_vertex(*edge) are connected by edge
        v = edge[0]
        u = self.target_vertex(edge)
        v_res_len, u_res_len = self.get_res_len(edge)
        #'edge[0] + (1<<edge[1])
        if v_res_len < u_res_len:
            if self.odd:
                edge_map = self.split_map(edge) #TAG FIELD UPDATE
                self.scalar_mult_algebra_map(self.get_edge_sign(edge),edge_map)
            else:
                edge_map = self.even_split_map(edge)
                self.scalar_mult_algebra_map(self.even_edge_sign(edge),edge_map)
            return edge_map
        else:
            if self.odd:
                base_map = [(self.get_circle(u,self.get_res(v)[i][0]),1) for i in range(v_res_len)]
                edge_map = self.induced_map(edge[0],base_map)
                self.scalar_mult_algebra_map(self.get_edge_sign(edge),edge_map)
            else:
                edge_map = self.even_merge_map(edge)
                self.scalar_mult_algebra_map(self.even_edge_sign(edge),edge_map)
            return edge_map

    def even_merge_map(self,edge):
        #print("Not assigned, Merge")
        v = edge[0]
        v_res = self.get_res(v)
        v_res_len = len(v_res)
        u = self.target_vertex(edge)
        target_circles = []
        i = j = -1
        for k in range(v_res_len):
            circle = self.get_circle(u,v_res[k][0])
            if i < 0 and circle in target_circles:
                #print("Assigned, Merge")
                i = target_circles.index(circle)
                j = len(target_circles)
            target_circles.append(circle)
        edge_map = [{} for _ in range(1<<v_res_len)]
        for g in range(1<<v_res_len):
            targets = self.even_merge(g,i,j)
            for target in targets:
                edge_map[g][target] = 1
        return edge_map

    def even_split_map(self,edge):
        #print("Not assigned, Split",edge)
        v = edge[0]
        v_res_len = self.get_res_len(v)
        u = self.target_vertex(edge)
        u_res = self.get_res(u)
        u_res_len = len(u_res)
        source_circles = []
        i = -1
        j = -1
        for k in range(u_res_len):
            circle = self.get_circle(v,u_res[k][0])
            if i < 0 and circle in source_circles:
                #print("Assigned, Split")
                i = source_circles.index(circle)
                j = len(source_circles)
                #print(i,j)
                break
            source_circles.append(circle)
        #print(source_circles)
        edge_map = [{} for _ in range(1<<v_res_len)]
        for g in range(1<<v_res_len):
            targets = self.even_split(g,i,j)
            for target in targets:
                edge_map[g][target] = 1
        return edge_map



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
                assert algebra_map[i][target] != 0
        return algebra_map

    def split_map(self,edge):
        assert (alg.v(edge[0],edge[1]) == 0)
        source = edge[0]
        target = edge[0] + (1<<edge[1])
        source_res = self.get_res(source)
        target_res = self.get_res(target)
        a = self.braid.get_x_diagram()[edge[1]]
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
                    #print_mapnt("map1[i][j+(1<<d[1])] = -1*imap[i][key]")
                    map1[i][key+(1<<d[1])] = imap[i][key]*sign
                #else:
                #    map1[i][0] = 0
        #print(imap)
        #print(d)
        #print(map0,map1)
        self.add_algebra_maps(map0,map1)
        return map0


    #Even Khovanov Homology

    def even_edge_sign(self,edge):
        i = 1
        for k in range(edge[1]):
            if (edge[0] >> k) % 2 == 1:
                i *= -1
        return i 

    def even_merge(self,g,i,j):
        #g is binary number < (1<<res_len)
        #g>>i % 2 = 0 if i tensor factor is v-, = 1 if i tensor factor is v+
        #i is part merge goes to
        #j is part merged away
        if g&((1<<i)|(1<<j)):
            if not(g & (1<<i) and  g & (1<<j)):
                g = g & ~(1<<i)
            return [((g&~((2<<j)-1))>>1) | (g&((1<<j)-1))]
        return []

    def even_split(self,g,i,j):
        assert i < j, "i = {}, j = {}".format(i,j)
        #g is binary number < (1<<res_len)
        #g>>i % 2 = 0 if i tensor factor is v-, = 1 if i tensor factor is v+
        #i is part to be split
        #j is tensor factor added in
        #by ordering on resolutions it is assumed i < j
        g_ = g & ~(1<<i)
        split = ((g_ & ~((1<<j)-1))<<1) | g_ & ((1<<j)-1)
        if (g >> i) % 2:
            return[split+(1<<i),split+(1<<j)]
        return [split]


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
        print("Braid:",self.word)
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