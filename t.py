import squarestruct as sq

dim = 4

s = sq.SquareStruct(dim)

def bs(v):
    s = ""
    for i in range(dim):
        s += str((v>>i)%2)
    return s

for v,i,j in s:
    print(v,i,j)