import fields as fe

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
            M[i][col], M[j][col] = M[j][col], M[i][col]

def swap_cols(M,i,j):
    if i != j:
        n = len(M)
        for row in range(n):
            M[row][i], M[row][j] = M[row][j], M[row][i]

def div_row(M,row,val):
    m = len(M[0])
    for col in range(m):
        M[row][col] = M[row][col]/val

def div_col(M,col,val):
    n = len(M)
    for row in range(n):
        M[row][col] = M[row][col]/val

def mul_row(M,row,val):
    n = len(M[0])
    for mid in range(n):
        M[row][mid] = M[row][mid]*val

def sub_row_mult(M,row,p_row,val):
    m = len(M[0])
    for col in range(m):
        M[row][col] = M[row][col] - (M[p_row][col] * val)

def sub_col_mult(M,col,p_col,val):
    n = len(M)
    for row in range(n):
        M[row][col] = M[row][col] - (M[row][p_col] * val)

def add_col_mult(M,col,p_col,val):
    n = len(M)
    for row in range(n):
        M[row][col] = M[row][col] + val * M[row][p_col]

def add_row_mult(M,row,p_row,val):
    m = len(M[0])
    for col in range(m):
        M[row][col] =  M[row][col] + val * M[p_row][col]

def anti_sub_row_mult(M,row,p_row,val):
    n = len(M[0])
    for col in range(n):
        M[row][col] = M[row][col] + (M[p_row][col] * val)

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

def field_row_echelon(A):
    n = len(A)
    m = len(A[0])
    pivot_row = 0
    for col in range(m):
        pivot_not_found = True
        #find non-zero entry:
        for row in range(pivot_row,n):
            if A[row][col].nonzero():
                if pivot_not_found:
                    #print_mat_(A)
                    #print()
                    #print("Swap {} and {}".format(pivot_row,row))
                    swap_rows(A,pivot_row,row)
                    
                    #print_mat_(A)
                    #print()
                    #print("Divide row {} by {}".format(pivot_row,A[pivot_row][col]))
                    div_row(A,pivot_row,A[pivot_row][col])

                    pivot_not_found = False
                    pivot_row += 1
                else:
                    #print_mat_(A)
                    #print()
                    #print("Sub divide row {} with pivot {} and val {}".format(row,pivot_row,A[row][col]))
                    sub_row_mult(A,row,pivot_row-1,A[row][col])
    #REMOVE
    #print_mat_(A)
    #######
    return pivot_row
        #put into pivot position

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

def i_row_reduce(M):
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
        last_pivot_col = pivot_col
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
        if pivot_col > last_pivot_col:
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
            for j in range(m-1):
                if M[i][j] != 0:
                    return True, int(M[i][j])
            return -3
    return -5

def mat_str(A):
    s = ""
    for row in A:
        s += " ".join(["{:>3}".format(col) for col in row]) + "\n"
    return s

def print_mat(A):
    print(mat_str(A))

def sim_int_red(A,B):
    return A, B 

def field_homology(A,B,not_field = False):
    if not_field:
        A = [[fe.FE(el,0) for el in row] for row in A]
        B = [[fe.FE(el,0) for el in row] for row in B]
    n = len(A)
    m = len(B)
    k = rank_A = rank_B = 0
    if B:
        k = len(B[0])
    elif A:
        m = len(A[0])
    if n and m:
        rank_A = field_row_echelon(A)
    if m and k:
        rank_B = field_row_echelon(B)
    ker = m - rank_A
    im = rank_B
    H = ker - im
    return H

def integer_homology(A,B,k,m,n):
    assert n == len(A) or len(A) == 0, "n = {}, len(A) = {}".format(n,len(A))
    assert m == len(B) or len(B) == 0, "m = {}, len(B) = {}".format(m,len(B))
    assert not A or m == len(A[0])
    assert not B or k == len(B[0])
    rank_A = rank_B = 0
    tor = []
    if A:
        A_Q = [[fe.FE(col,0) for col in row] for row in A]
        rank_A = field_row_echelon(A_Q)
    if B:
        S,D,_ = smith_normal_form(B)
        d = [D[i][i] for i in range(min(m,k)) if D[i][i] != 0]
        rank_B = len(d)
        tor = [i for i in d if (i != 1)]
    ker = m - rank_A
    betti = ker - rank_B
    #print("Rank A: {}, Rank_B: {}, ker: {}".format(rank_A,rank_B,ker))
    return betti, tor

def homology(A,B):
    print ("\n\n\nB:")
    print_mat(B)
    print("A:")
    print_mat(A)
    if not B:
        if not A:
            return "0"
        #A = int_reduce(A)
        last_pivot_col = -1
        for i in range(min(len(A),len(A[0]))):
            if A[i][i] != 0:
                last_pivot_col = i
        ker = len(A[0]) - last_pivot_col - 1
        if ker:
            return " + ".join(["Z" for _ in range(ker)])
        return "0"
        #homology of one, ie. dim of ker of A raised over Z
    A,B = sim_int_red(A,B)
    n = len(A)
    m = len(B)
    k = len(B[0])
    assert len(B) == len(A[0])
    last_pivot_col = -1
    for i in range(min(n,m)):
        if A[i][i] != 0:
            last_pivot_col = i
    ker = m - last_pivot_col - 1
    im  = [B[i][i] for i in range(min(m,k)) if B[i][i] != 0]
    s = ""
    for i in range(len(im)):
        if im[i] == 1:
            ker -= 1
            im[i] = 0
        else:
            s += "Z_{} + ".format(im[i])
            ker -= 1
    s += "Z + "*ker
    if s:
        s = s[:-3]
    else:
        s = "0"
    return s

def functi(A,B):
    n = len(A)
    m = len(B)
    assert m == len(A[0])
    k = len(B[0])
    pivot_row = 0
    for pivot_mid in range(m):
        for col in range(pivot_mid, m):
            if A[pivot_mid][col] != 0:
                if col != pivot_mid:
                    swap_cols(A,pivot_mid,col)
                    swap_rows(B,pivot_mid,col)
                    for col in range(pivot_mid+1,m):
                        pass
                break

def neg_row(A,row):
    m = len(A[0])
    for col in range(m):
        A[row][col] = -A[row][col]

def neg_col(A,col):
    n = len(A)
    for row in range(n):
        A[row][col] = -A[row][col]

def pivoted(A,pivot):
    n = len(A)
    m = len(A[0])
    for row in range(pivot+1,n):
        if A[row][pivot] != 0:
            return False
    for col in range(pivot+1,m):
        if A[pivot][col] != 0:
            return False
    return True

def next_loc(A,pivot):
    n = len(A)
    m = len(A[0])
    for row in range(pivot+1,n):
        for col in range(pivot+1,m):
            if A[row][col] % A[pivot][pivot] != 0:
                return (row,col)
    return None

def min_loc(A,pivot):
    n = len(A)
    m = len(A[0])
    loc = [pivot,pivot]
    min = float('inf')
    for row in range(pivot,n):
        for col in range(pivot,m):
            if A[row][col] != 0 and (abs(A[row][col]) <= min):
                loc = [row,col]
                min = abs(A[row][col])
    return loc

def I(n):
    return [[1 if i == j else 0 for i in range(n)] for j in range(n)]

def divisible(A,pivot):
    n = len(A)
    m = len(A[0])
    for i in range(pivot+1,n):
        for j in range(pivot+1,m):
            if A[i][j] % A[pivot][pivot] != 0:
                return (False,i,j,A[i][j]//A[pivot][pivot])
    return (True,0,0,0)

def smallest(v,ind):
    alpha = min(abs(v[i]) for i in range(ind,len(v)) if v[i] != 0)
    i = min(i for i in range(len(v)) if abs(v[i]) == alpha)
    return (alpha,i)

def min_nonzero_(A,pivot):
    n = len(A)
    m = len(A[0])
    v = []
    q = []
    for row in range(n):
        if row < pivot:
            v.append(0)
            q.append(0)
        else:
            a,b = smallest(A[row][:],pivot)
            v.append(a)
            q.append(b)
    alpha, i = smallest(v,pivot)
    return alpha, i, q[i]

def min_nonzero(A,pivot):
    m_row = m_col = pivot
    c_min = float('inf')
    for row in range(pivot):
        for col in range(pivot):
            if abs(A[row][col]) < c_min:
                m_row = row
                m_col = col
    return m_row, m_col

def completely_pivoted(A,pivot):
    n = len(A)
    m = len(A[0])
    for row in range(pivot+1,n):
        for col in range(pivot+1,m):
            if A[row][col] != 0:
                return False
    return True

def smith_normal_form(A):
    n = len(A)
    m = len(A[0])
    S = I(n)
    T = I(m)
    t = -1
    pivot = -1
    while not completely_pivoted(A,pivot):
        pivot += 1
        while True:
            #print_mat(A)
            #print()
            row, col = min_loc(A,pivot)
            swap_rows(A,pivot,row)
            swap_rows(S,pivot,row)
            swap_cols(A,pivot,col)
            swap_cols(T,pivot,col)
            for row in range(pivot+1,n):
                if A[row][pivot] != 0:
                    q = A[row][pivot]//A[pivot][pivot]
                    add_row_mult(A,row,pivot,-q)
                    add_row_mult(S,row,pivot,-q)
            if [A[row][pivot] for row in range(pivot+1,n) if A[row][pivot] != 0]:
                continue
                assert False
            for col in range(pivot+1,m):
                if A[pivot][col] != 0:
                    q = A[pivot][col]//A[pivot][pivot]
                    add_col_mult(A,col,pivot,-q)
                    add_col_mult(T,col,pivot,-q)
            if [A[pivot][col] for col in range(pivot+1,m) if A[pivot][col] !=0]:
                continue
                assert False
            div,row,col,q = divisible(A,pivot)
            if div:
                break
                assert False
            add_row_mult(A,row,pivot,1)
            add_row_mult(S,row,pivot,1)
            add_col_mult(A,col,pivot,-q)
            add_col_mult(T,col,pivot,-q)
        if A[pivot][pivot] < 0:
            neg_row(A,pivot)
            neg_row(S,pivot)
        if A[pivot][pivot] == 1:
            t += 1
    #print("S:")
    #tex_mat(S)
    #print("A:")
    #tex_mat(A)
    #print("T:")
    #tex_mat(T)
    return S,A,T

def _smith_normal_form(A):
    n = len(A)
    m = len(A[0])
    S = I(n)
    T = I(m)
    #print_mat(A)
    for pivot in range(min(n,m)):
        while not pivoted(A,pivot):
            row,col = min_loc(A,pivot)
            swap_rows(A,pivot,row)
            swap_rows(S,pivot,row)
            swap_cols(A,pivot,col)
            swap_cols(T,pivot,col)
            for row in range(pivot+1,n):
                if A[row][pivot] != 0:
                    val = A[row][pivot] // A[pivot][pivot]
                    add_row_mult(A,row,pivot,-val)
                    add_row_mult(S,row,pivot,-val)
            for col in range(pivot+1,m):
                if A[pivot][col] != 0:
                    val = A[pivot][col] // A[pivot][pivot]
                    add_col_mult(A,col,pivot,-val)
                    add_col_mult(T,col,pivot,-val)
            if pivoted(A,pivot):
                loc = next_loc(A,pivot)
                if loc:
                    row, _ = loc
                    add_row_mult(A,row,pivot,1)
                    add_row_mult(S,row,pivot,1)
                else:
                    if A[pivot][pivot] < 0:
                        neg_row(A,pivot)
                        neg_row(S,pivot)
    print("S:")
    tex_mat(S)
    print("A:")
    tex_mat(A)
    print("T:")
    tex_mat(T)
    return S,A,T        
    #print_mat(A)
    #print_mat(S)
    #print_mat(T)

def solve_mat_eq(S,D,T,r,y):
    #Solves Ax = y where SAT = D is the Smith normal form of A. and r is number of nonzero pivots of A
    #A is m x k
    m = len(D)
    k = len(D[0])
    #S is m x m
    assert m == len(S)
    #T is k x k
    assert k == len(T)
    #y is m x 1
    assert m == len(y)
    #x is k x 1
    #Check Sy (m x 1) is 0 for i >= r (indexed from 0)
    for i in range(r,m):
        sum = 0
        for j in range(m):
            sum += S[i][j]*y[j]
        if sum != 0:
            return False
    P = []
    for i in range(r):
        sum = 0
        for j in range(m):
            sum += S[i][j]*y[j]
        if sum % D[i][i] != 0:
            return False
        P.append(sum//D[i][i])
    #print("P=D^-1 ^Sy=",P)
    #print("r=",r)
    T_ = []
    for i in range(r):
        sum = 0
        for j in range(k):
            sum += T[i][j]*P[j]
        T_.append(sum)
    #print("^T=",T_)
    #print("Ax=y when x is a linear combination of:")
    for i in range(r):
        l = [0 if i >= r else T_[i]] 
        l += [T[i][j] for j in range(r,k)]
        s = "  ".join(["[{:>2}]".format(e) for e in l])
        #

        print(s)
    return True
    #Compute ^Sy
    #Check that d[i] = A[i][i] divides all in row [i] of Sy

def solve_mat_eq_gcd(S,T,d,y):
    pass

def solve_mat_mult(S,d,y):
    #Finds smallest positive integer s such that Ax = sy has a solution
    #If there is no such solution, returns 0
    #Assumes SAT = D is the Smith normal form of A and r is the number of nonzero pivots of A
    #D = diag(d_1,..,d_r,0,..,0)
    #d = [d_1,..,d_r]
    r = len(d)
    #A is m x k, S is m x m, T is k x k
    m = len(S)
    #y is m x 1
    assert m == len(y)
    #Check Sy (m x 1) is 0 for i >= r (indexed from 0)
    for i in range(r,m):
        sum = 0
        for j in range(m):
            sum += S[i][j]*y[j]
        if sum != 0:
            #print(i,sum)
            return 0
    s = 1
    #we need for each i in range(r) that d[i]| s(Sy)[i]
    for i in range(r):
        assert d[i] != 0, "i = {}".format(i)
        sum = 0
        for j in range(m):
            sum += S[i][j]*y[j]
        #print("i={}, s={}, sum={}, d[i]={},(s*sum)%d[i]={}".format(i,s,sum,d[i],(s*sum)%d[i]))
        if (s*sum) % d[i] != 0:
            #print(s,sum,d[i])
            s *= (d[i]//gcd(s*sum,d[i]))
    return abs(s)


def simultaneous_smith(A,B):
    n = len(A)
    m = len(B)
    k = len(B[0])
    S = I(n)
    S_ = I(n)
    T = I(m)
    print("A")
    print_mat(A)
    print("B")
    print_mat(B)
    for pivot in range(min(m,k)):
        while not pivoted(B,pivot):
            row,col = min_loc(B,pivot)
            swap_rows(B,pivot,row)
            swap_rows(S,pivot,row)
            swap_cols(B,pivot,col)
            swap_cols(T,pivot,col)
            swap_rows(A,pivot,col)
            swap_cols(S_,pivot,row)
            print("Swap, pivot = {}, other = {}".format(pivot,row))
            print_pair(S,S_)
            for row in range(pivot+1,m):
                if B[row][pivot] != 0:
                    val = B[row][pivot] // B[pivot][pivot]
                    add_row_mult(B,row,pivot,-val)
                    add_row_mult(S,row,pivot,-val)
                    add_col_mult(S_,pivot,row,val)
                    print("Add mult, pivot = {}, other= {}, val = {}".format(pivot,row,-val))
                    print_pair(S,S_)
            for col in range(pivot+1,k):
                if B[pivot][col] != 0:
                    val = B[pivot][col] // B[pivot][pivot]
                    add_col_mult(B,col,pivot,-val)
                    add_row_mult(A,pivot,col,val)
                    add_col_mult(T,col,pivot,-val)
            if pivoted(B,pivot):
                loc = next_loc(B,pivot)
                if loc:
                    row, _ = loc
                    add_row_mult(B,row,pivot,1)
                    add_row_mult(S,row,pivot,1)
                    add_col_mult(S_,pivot,row,-1)
                    print("All mult pivoted, pivot = {}, other = {}, val = {}".format(pivot,row,1))
                    print_pair(S,S_)
                else:
                    if B[pivot][pivot] < 0:
                        neg_row(B,pivot)
                        neg_row(S,pivot)
                        neg_col(S_,pivot)
                        print("Neg, pivot = {}".format(pivot))
                        print_pair(S,S_)
    print("B")                  
    print_mat(B)
    print("A")
    print_mat(A)
    print("S")
    print_mat(S)
    print("S^-1")
    print_mat(S_)
    print("SS^-1")
    print_mat(mat_mult(S,S_))
    print("T")
    print_mat(T)

def print_pair(A,B):
    n = len(A)
    assert n == len(B)
    print("S " + " "*(len(A[0]))*3 + "S^-1")
    for row in range(n):
        print("".join(["{:>3}".format(col) for col in A[row]]) + "     " + "".join(["{:3}".format(col) for col in B[row]]))
    print("")

def mat_mult(A,B):
    n = len(A)
    m = len(B)
    k = len(B[0])
    AB = [[0 for j in range(k)] for i in range(n)]
    for row in range(n):
        for col in range(k):
            for j in range(m):
                AB[row][col] += A[row][j]*B[j][col]
    return AB

def simultaneous_int_reduce(A,B,f,index):
    n = len(A)
    m = len(B)
    if B:
        k = len(B[0])
    else:
        k = 0
    assert len(B) == len(A[0])
    S = [[0 for _ in range(m)] for _ in range(m)]
    S_ = [[0 for _ in range(m)] for _ in range(m)]
    for i in range(m):
        S[i][i] = 1
        S_[i][i] = 1
    pivot_row = 0
    for pivot_mid in range(m):
        if pivot_row >= m:
            break
        j = pivot_mid
        while A[pivot_row][j] == 0:
            j += 1
            if j >= m:
                j = pivot_mid
                pivot_row += 1
                if pivot_row >= k:
                    pivot_row = -1 
                    break
        if pivot_row < 0:
            break
        if j != pivot_mid:
            swap_cols(A,pivot_mid,j)
            swap_rows(B,pivot_mid,j)
        if A[pivot_row][pivot_mid] < 0:
            for l in range(pivot_row,j):
                A[j][pivot_mid] = -A[j][pivot_mid]
        #d = reduce_col(A,pivot_mid)
        #if d != 1:
        #    mul_row(B,pivot_mid, d)
        for j in range(m):
            d = A[pivot_mid][j] 
            pivot = A[pivot_row][pivot_mid]
            if j != pivot_mid and d != 0:
                for l in range(): 
                    A[l][j] = pivot*A[l][j] - d*A[l][pivot_mid]
        pivot_row += 1
    return A,B     
        
def simultaneous_reduce(A,B):#    ,f,index):
    # Convention:  C_(i-1) -- B --> C_i -- A --> C_(i+1), so AB = 0
    # Dimensions:     k     (mxk)    m   (nxm)      n        (nxk)
    n = len(A)
    m = len(B)
    #assert len(A[0]) == m + 1
    k = len(B[0])
    tex_mat(A)
    tex_mat(B)
    ### TESTING STUFF
    ##char = A[0][0].c
    ##S = [[fe.FE(0,char) for _ in range(m)] for _ in range(m)]
    ##S_= [[fe.FE(0,char) for _ in range(m)] for _ in range(m)]
    ##for i in range(m):
    ##    S[i][i].n = 1
    ##    S_[i][i].n = 1
    #assert (m == len(A[0])), "These matrices do not compose together."
    #first writing out the column reduction by transposition
    i = 0
    #check each column for a pivot
    """
    f(A)
    print()
    f(B)
    print()
    """
    for row in range(n):
        #look for pivot in altered rows only:
        for mid in range(i,m):
            val = A[row][mid]
            #first non-zero value in current row in a col not already a pivot
            if val.nonzero():
                #move col into next pivot position
                if mid != i:
                    swap_cols(A,mid,i)
                    ##swap_cols(S,mid,i) #TEST STUFF
                    swap_rows(B,mid,i)
                    ##swap_rows(S_,mid,i) #TEST STUFF
                #change new pivot row value so pivot is 1
                if val != val/val:
                    div_col(A,i,val)
                    ##div_col(S,i,val) #TEST STUFF
                    #This is mul and not div right?
                    mul_row(B,i,val)
                    ##mul_row(S_,i,val) #TEST STUFF
                #kill off values after new pivot in same column
                for _mid in range(m):
                    if _mid != i:
                        _val = A[row][_mid]
                        if _val.nonzero():
                            sub_col_mult(A,_mid,i,_val) 
                            ##sub_col_mult(S,_mid,i,_val) #TEST STUFF
                            anti_sub_row_mult(B,i,_mid,_val)
                            ##anti_sub_row_mult(S_,i,_mid,_val) #TEST STUFF
                i += 1
                break
    ker_A = 0
    for j in range(m):
        zero_col = True
        for i in range(n):
            if A[i][j].nonzero():
                zero_col = False
                break
        if zero_col:
            ker_A += 1
    im_B = 0
    for row in B:
        if [el for el in row if el.nonzero()]:
            im_B += 1
    #print("After reduction:")
    """
    f(A)
    print()
    f(B)
    print()
    f(S)
    print()
    f(S_)
    print()
    """
    """
    for i in range(m):
        for j in range(m):
            s = fe.FE(0,char)
            for l in range(m):
                s += S[i][l] * S_[l][j]
            if i == j:
                assert s.n == 1, "Inverse problem in product {},{}. Yields {}, not 1".format(i,j,s.n)
            else:
                assert s.n == 0, "Inverse problem in product {},{}. Yields {}, not 0".format(i,j,s.n)
    """
    #transformed_inv = [S_[i][index] for i in range(m)]
    try:
        print("m = {}, n = {}, k = {}, ker = {}, im = {}".format(m,n,k,ker_A,im_B))
    except:
        pass
    tex_mat(A)
    tex_mat(B)
    return ker_A, im_B
        #print_mat(A)

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


def tex_mat(M):
    n = len(M)
    m = 0
    if n:
        m = len(M[0])
    else:
        return 0
    if m > 17:
        print("LARGE")
    else:
        print("\\[ \\begin{array}{" + "r"*m + "}")
        for row in range(n):
            for col in range(m-1):
                print(M[row][col],end="&")
            print(M[row][m-1],"\\\\")
        print("\\end{array}\\]")

"""
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


def print_mat_(M):
    for row in M:
        for col in row:
            print(int(col), end="  ")
        print()
