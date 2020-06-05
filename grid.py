class Grid:
    
    def __init__(self, orow, xrow):
        self.orow = orow
        self.xrow = xrow
        assert (len(self.orow) == len(self.xrow))
        self.n = len(self.orow)
        self.ocol = [None for _ in range(self.n)]
        self.xcol = [None for _ in range(self.n)]
        for i in range(self.n):
            assert self.orow[i] != self.xrow[i]
            self.ocol[self.orow[i]] = i
            self.xcol[self.xrow[i]] = i
        
    def braid_R(self):
        braid = []
        strands = []
        #get starting strands
        for i in range(self.n):
            if self.ocol[i] > self.xcol[i]:
                strands.append(i)
        for i in range(self.n):
            assert self.xrow[i] in strands
            o = strands.index(self.xrow[i])
            strands[o] = self.orow[i]
            s = 1 if self.xrow[i] < self.orow[i] else -1
            for j in range(1,len(strands)):
                if strands[j-1] > strands [j]:
                    strands[j-1],strands[j] = strands[j],strands[j-1]
                    braid.append(s*j)
        return braid
    
    def __str__(self):
        s = ""
        for i in range(self.n):
            r = ['⬜️' for _ in range(self.n)]
            r[self.xcol[i]] = '❌'
            r[self.ocol[i]] = '⭕️'
            s += ''.join(r) + '\n'
        return s
    
    def tikz(self,centered=False):
        s = "\\begin{tikzpicture}\n"
        x = "\\draw ({},{}) ++(-.3,-.3) -- ++ (.6,.6) ++(-.6,0) --++(.6,-.6);\n"
        o = "\\draw ({},{}) circle (.3cm);\n"
        ha = "\\draw[->] ({},{}) ++ ({}*0.5,0) --++({},0);\n"
        vb = "\\draw[white, line width = 5pt] ({},{}) ++ (0,{}*0.5) --++ (0,{});\n"
        va = "\\draw[->] ({},{}) ++ (0,{}*0.5) --++(0,{});\n"
        for i in range(self.n):
            oc = self.ocol[i]
            xc = self.xcol[i]
            s += o.format(oc,self.n-i)
            s += x.format(xc,self.n-i)
            sgn = 1 if xc > oc else -1
            s += ha.format(oc,self.n-i,sgn,xc - oc - sgn)
        for i in range(self.n):
            orow = self.orow[i]
            xrow = self.xrow[i]
            sgn = 1 if xrow > orow else -1
            s += vb.format(i,self.n-xrow,sgn,xrow - orow - sgn)
            s += va.format(i,self.n-xrow,sgn,xrow - orow - sgn)
        s += "\\end{tikzpicture}"
        if centered:
            s = "\\begin{center}" + s + "\\end{center}"
        return s
    
    def tikz_legendrian(self,centered=False):
        s = "\\begin{tikzpicture}\n"
        brc = "\\draw[white, shorten <= -5pt, shorten >= -5pt, line width = 5pt] ({},{}) to[out=0, in=180] ({},{});\n"
        arc = "\\draw ({},{}) to[out=0, in=180] ({},{});\n"
        for i in range(self.n):
            l = min(self.ocol[i],self.xcol[i])
            r = max(self.ocol[i],self.xcol[i])
            s += arc.format(l+i,l-i,r+i,r-i)
        for i in range(self.n):
            l = min(self.orow[i],self.xrow[i])
            r = max(self.orow[i],self.xrow[i])
            s += brc.format(i+l,i-l,i+r,i-r)
            s += arc.format(i+l,i-l,i+r,i-r)
        s += "\\end{tikzpicture}"
        if centered:
            s = "\\begin{center}" + s + "\\end{center}"
        return s
    
L1xrow = [9,2,7,3,0,6,8,4,5,1]#[10,3,8,4,1,7,9,5,6,2]
L1orow = [4,8,0,1,2,9,5,7,3,6]#[5,9,1,2,3,10,6,8,4,7]

L2xrow = [9,4,7,5,2,6,1,3,8,0]#[10,5,8,6,3,7,2,4,9,1]
L2orow = [6,8,2,3,4,0,5,9,1,7]#[7,9,3,4,5,1,6,10,2,8]

L1 = Grid(L1orow,L1xrow)
print("L1")
print(L1)
bL1 = L1.braid_R()
print(bL1)

L2 = Grid(L2orow,L2xrow)
print("L2")
print(L2)
bL2 = L2.braid_R()
print(bL2)


