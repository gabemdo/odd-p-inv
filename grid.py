import braid as br

class Grid:
    #NOTE: The convention here is that the xcol and ocol lists are the indices (indexed from 0) if the x's and o's as given row by row.
    def __init__(self, xcol, ocol):
        self.xcol = xcol
        self.ocol = ocol
        assert (len(self.ocol) == len(self.xcol))
        self.n = len(self.xcol)
        self.xrow = [None for _ in range(self.n)]
        self.orow = [None for _ in range(self.n)]
        for i in range(self.n):
            assert self.xcol[i] != self.ocol[i]
            self.xrow[self.xcol[i]] = i
            self.orow[self.ocol[i]] = i

    def rotate_grid(self,n=1):
        #rotates the grid clockwise by n*90degrees.
        if n%4 == 0:
            return self
        elif n%4 == 1:
            new_xcol = [self.n - self.orow[i] - 1 for i in range(self.n)]
            new_ocol = [self.n - self.xrow[i] - 1 for i in range(self.n)]
        elif n%4 == 2:
            new_xcol = [self.n - self.xcol[self.n-i-1] - 1 for i in range(self.n)]
            new_ocol = [self.n - self.ocol[self.n-i-1] - 1 for i in range(self.n)]
        else:
            new_xcol = self.orow[::-1]
            new_ocol = self.xrow[::-1]
        return Grid(new_xcol,new_ocol)

    def braid(self,direction='r'):
        if isinstance(direction,int):
            braid = self.rotate_grid(direction).braid_R()
            #return br.Braid(self.rotate_grid(direction).braid_R())
        else:
            d = direction.strip().lower()[0]
            if d == 'r':
                braid = self.braid_R()
                #return br.Braid(self.braid_R())
            if d == 'u':
                braid = self.rotate_grid().braid_R()
                #return br.Braid(self.rotate_grid().braid_R())
            if d == 'l':
                braid = self.rotate_grid(2).braid_R()
                #return br.Braid(self.rotate_grid(2).braid_R())
            if d == 'd':
                braid = self.rotate_grid(3).braid_R()
                #return br.Braid(self.rotate_grid(3).braid_R())
        #return braid 
        return br.Braid(braid)
        #assert False, "The argument should be an integer representing the number of +90 degree rotations of the direction, or a string such as 'r', 'u', 'l', 'd'. 'Rightward', 'right', 'RIGHT', etc., also work."
        
    def braid_R(self):
        braid = []
        strands = []
        #get starting strands
        #look for strands as leftward arrows that have been split
        for i in range(self.n):
            if self.ocol[i] > self.xcol[i]:
                strands.append(i)
        s = len(strands)
        #print(strands)
        #after each column in the grid the strands are in new positions
        #first we see which strand is moving: the strand in the new x row
        #then we see which row that strand is going to: the new o row
        #if orow < xrow, the crossings are negative and happen higher strand index to lower
        #if orow > xrow, the crossings are positive and happen lower strand index to higher
        for i in range(self.n):
            assert self.xrow[i] in strands
            strand_row = self.xrow[i]
            s_index = strands.index(self.xrow[i])
            new_strand_row = self.orow[i]
            strands[s_index] = new_strand_row
            while s_index and strands[s_index] < strands[s_index - 1]:
                strands[s_index-1], strands[s_index] = strands[s_index], strands[s_index-1]
                braid.insert(0,-s_index)
                s_index -= 1
            while s_index < s - 1 and strands[s_index] > strands[s_index + 1]:
                strands[s_index], strands[s_index+1] = strands[s_index+1], strands[s_index]
                braid.insert(0,s_index+1)
                s_index += 1
        #Reduce braid
        reduced = True
        while reduced:
            reduced = False
            i = 0
            stop = len(braid) - 1
            while i < stop:
                if braid[i] == -braid[i+1]:
                    braid.pop(i)
                    braid.pop(i)
                    stop -= 2
                    reduced = True
                else:
                    i += 1
        return braid
    
    def __str__(self):
        s = ""
        for i in range(self.n):
            r = ['⬜️' for _ in range(self.n)]
            r[self.xcol[i]] = '❌'
            r[self.ocol[i]] = '⭕️'
            s += ''.join(r) + '\n'
        return s
    
    def tex_grid(self,centered=False):
        s = "\\begin{tikzpicture}\n\\begin{scope}[shift = {(-0.5,0.5)}]\n"
        s += "\\draw (0,0) grid ({0},{0});\n\\end{{scope}}\n".format(self.n)
        x = "\\draw ({},{}) ++(-.3,-.3) -- ++ (.6,.6) ++(-.6,0) --++(.6,-.6);\n"
        o = "\\draw ({},{}) circle (.3cm);\n"
        for i in range(self.n):
            oc = self.ocol[i]
            xc = self.xcol[i]
            s += o.format(oc,self.n-i)
            s += x.format(xc,self.n-i)
        s += "\\end{tikzpicture}"
        if centered:
            s = "\\begin{center}" + s + "\\end{center}"
        return s

    def tex_knot(self,centered=False):
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

    def tex_braid(self,centered=False):
        s = "\\begin{tikzpicture}\n"
        ha = "\\draw ({},{}) --++({},0);\n"
        hx = "\\draw[->] (-1,{1}) -- ({3},{1}) ({0},{1}) -- ({2},{1});"
        vb = "\\draw[shorten >=0.25cm,shorten <=0.25cm, white, line width = 5pt] ({},{}) --++ (0,{});\n"
        va = "\\draw ({},{}) --++(0,{});\n"
        for i in range(self.n):
            oc = self.ocol[i]
            xc = self.xcol[i] 
            if xc > oc: 
                s += ha.format(oc,self.n-i,xc - oc)
            else: 
                s += hx.format(oc,self.n-i,self.n+1,xc)
        for i in range(self.n):
            orow = self.orow[i]
            xrow = self.xrow[i]
            s += vb.format(i,self.n-xrow,xrow - orow)
            s += va.format(i,self.n-xrow,xrow - orow)
        s += "\\end{tikzpicture}"
        if centered:
            s = "\\begin{center}" + s + "\\end{center}"
        return s
    
    def tex_Legendrian_front(self,centered=False):
        s = "\\begin{tikzpicture}\n"
        brc = "\\draw[white, line width = 5pt] ({},{}) to[out=0, in=180] ({},{});\n"
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
   

def main(): 
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

if __name__ == '__main__':
    main()


