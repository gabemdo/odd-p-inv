

def choose(n,k):
        a = b = 1
        if k > n//2:
            return choose(n,n-k)
        for i in range(k):
            a *= (n-i)
            b *= i+1
        return a//b

def abs(x):
    if x >= 0:
        return x
    return -x

def gcd(a,b):
    if b == 0:
        return a
    return gcd(b, a%b)

def v(v,i):
    return (v >> i) % 2

def height(b):
    #       0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
    hght = [0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4]
    sum = 0
    while b:
        sum += hght[b%16]
        b //= 16
    return sum

def difference(v,u):
    assert(u > v)
    return diff(u-v)

def diff(delta,mod = 0):
    assert (delta >= 0)
    edge = {(1<<i): i for i in range(16)}
    if delta % (1<<16) == 0 & delta >> 16:
        return diff(delta>>16,mod + 16)
    if delta in edge:
        return edge[delta] + mod
    return -1

def swap_rows(M,i,j):
    if i != j:
        m = len(M[0])
        for col in range(m):
            M[i][col],M[j][col] = M[j][col],M[i][col]

def div_row(M,row,val):
    m = len(M[0])
    for col in range(m):
        M[row][col] /= val

def reduce_row(M,row):
    m = len(M[0])
    d = 0
    #Find smallest magnitude nonzero entry in row.
    for i in range(m):
        if M[row][i] != 0:
            if d == 0:
                d = abs(M[row][i])
            else: 
                d = min(d,abs(M[row][i]))
    #If all zero, done.
    if d == 0:
        return d
    #Find gcd of row.
    for i in range(m):
        d = gcd(M[row][i],d)
    #Do nothing if gcd = 1
    if d == 1:
        return d
    #Divide by gcd
    for i in range(m):
        M[row][i] //= d
        return d

def int_row_reduce(M,col_labels,dim):
    #M is an nxm matrix
    n = len(M)
    m = len(M[0])
    pivot_col = 0
    last_pivot_col = -1
    factor = 1
    for pivot_row in range(n):
        #check if we are beyond final column
        if pivot_col >= m:
            break
        #search for next non-zero entry in matrix below and right of last pivot, read downward by columns left to right.
        k = pivot_row
        while M[k][pivot_col] == 0:
            k += 1
            #if we went beyond last row, go to next column
            if k >= n:
                k = pivot_row
                pivot_col += 1
                #if we went beyond last col, break
                if pivot_col >= m:
                    pivot_col = -1
                    break
        #if we went beyond last col, second break
        if pivot_col < 0:
            break
        #if here we eventually found a pivot. last_pivot_col is here now.
        last_pivot_col = pivot_col
        #if new pivot row is not in next row position swap
        if k != pivot_row:
            swap_rows(M,pivot_row,k)
        #make pivot positive
        if M[pivot_row][pivot_col] < 0:
            for j in range(pivot_col,m):
                M[pivot_row][j] = -M[pivot_row][j]
        #make sure gcd of pivot row is 1
        if pivot_row == m-1:
            factor = M[pivot_row][pivot_col]
        reduce_row(M,pivot_row)
        #This is for reduced REF so we cancel all rows
        for k in range(n):
            d = M[k][pivot_col]
            pivot = M[pivot_row][pivot_col]
            if k != pivot_row and d != 0:
                for j in range(m):
                    M[k][j] = pivot*M[k][j] - d*M[pivot_row][j]
                #Make sure changed row has gcd 1
                reduce_row(M,k)
        pivot_col += 1
    if last_pivot_col == m-1:
        return False, factor
    assert last_pivot_col >= 0
    for i in range(n):
        if M[i][m-1] != 0:
            print("The invariant over $\\mathbb Z$ is zero as the image:")
            alpha = []
            for j in range(m-1):
                if M[i][j]:
                    alpha.append((M[i][j],j))
            print("\\[",' + '.join('{} ({},{}) '.format(a,str_v(col_labels[b][0],dim),str_g(col_labels[b][1])) for a,b in alpha),"=\\psi\\]")
            for j in range(m-1):
                if M[i][j] != 0:
                    return True, int(M[i][j])
            return -3
    return -5

def sub_row_mult(M,row,p_row,val):
    m = len(M[0])
    for col in range(m):
        M[row][col] = M[row][col] - (M[p_row][col] * val)



def dumber_row_reduce(M,col_labels,d):
    n = len(M)
    m = len(M[0])
    i = 0
    factor = 1
    last_pivot_col = -1
    #check each column for a pivot
    for col in range(m):
        #look for pivot in altered rows only
        for row in range(i,n):
            val = M[row][col]
            #first non-zero value in current column in a row that is not already a pivot
            if val.nonzero():
                last_pivot_col = col
                #move row into next pivot position
                if row != i:
                    swap_rows(M,row,i)
                #change new pivot row value so pivot is 1
                if val != val/val:
                    if col == m-1:
                        factor = val.n
                    div_row(M,i,val)
                #kill off values after new pivot in the same column
                for row in range(n):
                    if row != i:
                        val_ = M[row][col]
                        if val_.nonzero():
                            sub_row_mult(M,row,i,val_)
                i += 1
                break
    if last_pivot_col == m-1:
        return False, factor # (False,0)
    for i in range(n):
        if M[i][m-1]:
            if M[0][0].c == 0:
                print("The invariant over $\\mathbb Q$ is zero as the image:")
            else:
                print("The invariant over $\\mathbb Z/{}$ is zero as the inage:".format(M[0][0].c))
            alpha = []
            for j in range(m-1):
                if M[i][j]:
                    alpha.append((M[i][j],j))
            print("\\[",' + '.join('{} ({},{}) '.format(a,str_v(col_labels[b][0],d),str_g(col_labels[b][1])) for a,b in alpha),"=\\psi\\]")
            for j in range(m-1):
                if M[i][j]:
                    return True, int(M[i][j]) #(True, int(M[i][j]))
            return -3 #(True, "ERROR no invariant in col")
    return -5 #(True, "ERROR no factor in row")

def str_v(v,d):
    s = ""
    for i in range(d):
        s += '1' if v & (1<<i) else '0'
    return s

def str_g(g):
    if g == 0:
        return '1'
    l = []
    i = 1
    while g:
        if g%2 == 1:
            l.append('v_{{{}}}'.format(i))
        g = g>>1
        i += 1
    return '\\wedge '.join(l)

def print_mat(M):
    n = len(M)
    m = len(M[0])
    if m > 17:
        print("LARGE")
    else:
        print("\\[ \\begin{array}{" + "r"*m + "}")
        for row in range(n):
            for col in range(m-1):
                print(int(M[row][col]),end="&")
            print(int(M[row][m-1]),"\\\\")
        print("\\end{array}\\]")

def dumb_row_reduce(M):
    n = len(M)
    m = len(M[0])
    i = 0
    for col in range(m):
        for row in range(i,n):
            val = M[row][col]
            if val.nonzero():
                last_pivot_col = col
                #print("\nPivot {}: ({},{}) = {}".format(i,row,col,val))
                #print_mat(M)
                if row != i:
                    #print("\nSwap row {} and row {}".format(i,row))
                    swap_rows(M,row,i)
                    #print_mat(M)
                if val != val/val:
                    #print("\nDivide row {} by {}".format(i,val))
                    div_row(M,i,val)
                    #print_mat(M)
                for row in range(n):
                    if row != i:
                        val = M[row][col]
                        if val.nonzero():
                            #print("\nSubtract {} times row {} from row {}".format(val,i,row))                            
                            sub_row_mult(M,row,i,val)
                            #print_mat(M)
                i += 1
                break
    return last_pivot_col == m - 1
    print()
    print_mat(M)

"""
def print_mat(M):
    for row in M:
        for col in row:
            print(int(col), end="  ")
        print()
"""