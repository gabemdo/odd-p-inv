import fields as fe
_0 = fe.FE(0,3)
_1 = fe.FE(1,3)
_2 = fe.FE(2,3)

#A = [[_1,_2,_0,_2],[_2,_1,_1,_0],[_0,_0,_1,_0]]
A = [[fe.FE(1,3),fe.FE(0,3),fe.FE(1,3)],[fe.FE(0,3),-fe.FE(1,3),fe.FE(0,3)],[-fe.FE(1,3),fe.FE(0,3),fe.FE(0,3)]]


def change(M,row,col,d):
    assert (row < len(M))
    assert (col < len(M[0]))
    if d[1] == -1:
        M[row][col] *= d[0]
    else:
        M[row][col] += M[row][d[1]]*d[0]

def swap_rows(M,i,j):
    if i != j:
        m = len(M[0])
        print("Swap",i,j)
        for col in range(m):
            M[i][col],M[j][col] = M[j][col],M[i][col]

def div_row(M,row,val):
    m = len(M[0])
    for col in range(m):
        M[row][col] /= val

def sub_row_mult(M,row,mrow,val):
    m = len(M[0])
    for col in range(m):
        print(col,M[row][col],M[mrow][col],val)
        M[row][col] = M[row][col] - (M[mrow][col] * val)

def dumb_row_reduce(M):
    n = len(M)
    m = len(M[0])
    i = 0
    for col in range(m):
        for row in range(i,n):
            val = M[row][col]
            if val != _0:
                print("\nPivot:({},{}) = {}".format(i,col,val))
                print_mat(M)
                if row != i:
                    print("\nSwap row {} and row {}".format(i,row))
                    swap_rows(M,row,i)
                    print_mat(M)
                if val != val/val:
                    print("\nDivide row {} by {}".format(i,val))
                    div_row(M,i,val)
                    print_mat(M)
                for row in range(n):
                    if row != i:
                        val = M[row][col]
                        if val != _0:
                            print("\nSubtract {} times row {} from row {}".format(val,i,row))                            
                            sub_row_mult(M,row,i,val)
                            print_mat(M)
                i += 1
                break
    print()
    print_mat(M)




def col_reduce(M):
    n = len(M)
    m = len(M[0])
    pivots = {}
    delta = [[] for _ in range(m)]
    for row in range(n):
        assert (len(M[row]) == m)
        current_pivot = None
        found_pivot = False 

        for col in range(m):

            for d in delta[col]:
                change(M,row,col,d)
                print("Delta",row,col,d)
                print_mat(M)

            if (not found_pivot) and (col not in pivots) and not(-0.001 < M[row][col] < 0.001):

                found_pivot = True
                current_pivot = col
                pivots[col] = M[row][col]
                M[row][col] = 1
                assert (int(M[row][col] == 1))

                print("Pivot change",row,col,d)
                print_mat(M)

            for col in range(m):
                if col not in pivots and not(-0.001 < M[row][col] < 0.001):

                    d = (-M[row][col],current_pivot)
                    change(M,row,col,d)
                    delta[col].append(d)
                    assert (-0.001 < M[row][col] <0.001)

                    print("Row change",row,col,d)
                    print_mat(M)
                    
                    
    print(pivots)
    print(delta)
    print_mat(M)


def print_mat(M):
    for row in M:
        for col in row:
            print(int(col), end="  ")
        print()

print_mat(A)
dumb_row_reduce(A)
#print_mat(A)

