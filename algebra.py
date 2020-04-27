

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
    edge = {(1<<i): i for i in range(16)} #1: 0, 2: 1, 4: 2, 8: 3, 16: 4, 32: 5, 64: 6, 128: 7, 256: 8, 512: 9, 1024: 10, 
            #2048: 11, 4096: 12, 8192: 13, 16384: 14, 32768: 15}
    if delta % (1<<16) == 0 & delta >> 16:
        return diff(delta>>16,mod + 16)
    if delta in edge:
        return edge[delta] + mod
    return -1

