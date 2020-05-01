

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

def int_row_reduce(M):
    n = len(M)
    m = len(M[0])
    pivot = 0
    det = 1
    for i in range(n):
        if pivot >= m:
            break
        k = i
        while M[i,pivot] == 0:
            k += 1
            if k >= n:
                k = r
                pivot += 1
                if pivot >= m:
                    break
        if pivot < 0:
            break
        if k != i:
            swap_rows(M,i,k)
        if M[i][pivot] < 0:
            div_row(M,i,-1)
            det = -det
        det *= M[i][pivot]
        d = reduce_row(M,[i])
        for k in range(n):
            if i != k and M[k][pivot] != 0:
                for j in range(n):
                    M[k][j] = M[i][pivot]*M[k][j] - M[k][pivot]*M[i][j]
                d = reduce_row(M,k)
        pivot += 1
    return det

def sub_row_mult(M,row,p_row,val):
    m = len(M[0])
    for col in range(m):
        M[row][col] = M[row][col] - (M[p_row][col] * val)



def dumber_row_reduce(M):
    n = len(M)
    m = len(M[0])
    i = 0
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
                    div_row(M,i,val)
                #kill off values after new pivot in the same column
                for row in range(n):
                    if row != i:
                        val_ = M[row][col]
                        if val_.nonzero():
                            sub_row_mult(M,row,i,val_)
                i += 1
                break
    print("IN NEW ALG")
    if last_pivot_col == m-1:
        return 0# (False,0)
    for i in range(n):
        if M[i][m-1]:
            for j in range(m-1):
                if M[i][j]:
                    return int(M[i][j]) #(True, int(M[i][j]))
            return -3 #(True, "ERROR no invariant in col")
    return -5 #(True, "ERROR no factor in row")


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