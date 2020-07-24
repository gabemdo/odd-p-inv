import sys
import fields as fe
import edgestruct as es
import squarestruct as ss
import algebra as alg
import column as col

class Braid:
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
        assert (len(braid_word) > 0)
        self.word = braid_word
        #The braid index is inferred. Perhaps it is better to take it as an argument.
        self.b = self.comp_b(self.word)
        self.nmin = self.count_negative_crossings()
        self.d = len(self.word)
        self.x_diagram = []
        self.res = [[] for _ in range(1 << self.d)]
        self.res_len = [0 for _ in range(1 << self.d)]
        self.closure = []
        self.inv_v = -1



    def self_linking_number(self):
        return self.d - self.b - 2*self.nmin

    def mirror(self):
        word = [-el for el in self.word]
        return Braid(word)

    def reverse(self):
        word = self.word[::-1]
        return Braid(word)

    def __add__(self,other):
        #Warning! This is intended only for sums of knots not for links
        ##TODO: add component check and return add only if both component numbers are 1
        word = self.word
        b = self.b
        new = other.word
        for el in new:
            if el > 0:
                word.append(el + b - 1)
            else:
                word.append(el - b + 1)
        return Braid(word)

    def __str__(self):
        sl = self.self_linking_number()
        return "Braid word: B = {}, sl(B) = {}".format(self.word,sl)

    def get_word(self):
        return self.word

    def get_crossing_number(self):
        return self.d

    def get_b(self):
        return self.b

    def get_nmin(self):
        return self.nmin

    #Basic helper functions

    def comp_b(self,braid):
        #return the width of the braid
        i = 0
        for s in braid:
            c = alg.abs(s)
            if c > i:
                i = c
        return i+1

    def count_negative_crossings(self):
        #return n-
        n = 0
        for x in self.word:
            if x < 0:
                n += 1
        return n

    def get_inv_v(self):
        #return the vertex number for the invariant
        #if not self.inv_vertex:
        if self.inv_v < 0:
            #vertex = []
            v = 0
            #for x in self.word:
            #    if x > 0:
            #        vertex.append(0)
            #    else:
            #        vertex.append(1)
            for i in range(self.d):
                if self.word[i] < 0:
                    v += (1<<i)
            self.inv_v = v
            self.inv_r = alg.height(v)
        return self.inv_v

    def get_inv_r(self):
        if self.inv_v < 0:
            self.get_inv_v()
        return self.inv_r

    def tex_braid(self):
        print("\\begin{center}\\begin{tikzpicture}[scale = 0.4]")
        i = 0
        for x in self.word:
            if x > 0:
                print("\t\\braidline{{{}}}{{{}}}{{0}}{{{}}}{{1}}".format(self.b,x,i))
            else: 
                print("\t\\braidline{{{}}}{{{}}}{{0}}{{{}}}{{0}}".format(self.b,-1*x,i))
            i += 1
        print("\end{tikzpicture}\end{center}")


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
            
            l = len(self.word)
            for j in range(l):
                ####print(a)
                # x: current crossing
                x = self.word[j]
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
            if alg.v(v,i) == 0:
                a,b = self.zero_smoothing(self.x_diagram[i])
            else:
                a,b = self.one_smoothing(self.x_diagram[i])
            sym.append(a)
            sym.append(b)
        self.res[v] = self.sym_to_circles(sym)
        self.res_len[v] = len(self.res[v])
        #print(v,self.res[v])
        return self.res[v]

    def get_res_len(self,v):
        if isinstance(v,int):
            assert (v < (1<<self.d))
            if self.res_len[v]:
                return self.res_len[v]
            if self.res[v]:
                self.res_len[v] = len(self.res[v])
                return self.res_len[v]
            self.res_len[v] = len(self.get_res(v))
            return self.res_len[v]
        else:
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


