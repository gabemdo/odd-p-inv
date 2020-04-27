import algebra as alg

class Column:

    def __init__(self,d,res_lens):
        self.r = 0
        self.d = d
        self.res_lens = res_lens

    def __call__(self,r):
        self.r = r

    def __iter__(self):
        for v in range((1<<r)-1,(((1<<r)-1)<<(self.d-r))+1):
            if alg.height(v) == self.r:
                for g in range(1<<self.res_lens[v]):
                    yield (v,g)
