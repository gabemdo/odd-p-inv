class FE:
    #This is a class that contains field elements of variable characteristic.
    #For finite fields, characteristic = p, we initialize FE(value,p)
    #For rationals, we initialize FE(n,0) for integer n or FE((n,d),0) for rational n/d.
    def __init__(self, value, characteristic):
        self.c = characteristic
        if characteristic == 0:
            if isinstance(value, int):
                self.n = value
                self.d = 1
            else:
                assert(len(value) == 2)
                n = value[0]
                d = value[1]
                if n == 0:
                    self.n, self.d = 0, 1
                if d == 0:
                    raise ZeroDivisionError
                else:
                    g = self.gcd(n,d)
                    self.n, self.d = n//g, d//g
                    assert (self.d > 0), "Error with my gcd"
        elif characteristic == 2:
            assert (isinstance(value, int) or isinstance(value, bool))
            self.n = bool(value)
        else:
            assert (isinstance(value, int))
            assert (isinstance(characteristic, int))
            assert (self.is_prime(characteristic))
            self.n = value % characteristic
            self.i = 0

    def is_prime(self, p):
        if p in [2,3,5,7,11]:
            return True
        if p <= 1:
            return False
        for i in range(2, 1 + (p-1)//2):
            if p%i == 0: 
                return False
        return True

    def gcd(self, a, b):
        if b == 0:
            return a
        return self.gcd(b, a%b)

    def nonzero(self):
        return self.n != 0

    def __eq__(self,other):
        if self.c != other.c:
            return False
        if self.c != 0:
            return self.n == other.n
        return self.n == other.n and self.d == other.d

    def __repr__(self):
        if self.c == 0 and self.d != 1:
            return "{}(({},{}),{})".format(self.__class__.__name__, self.n, self.d, self.c)
        return "{}({},{})".format(self.__class__.__name__, self.n,self.c)

    def __str__(self):
        if self.c == 0:
            if self.d == 1:
                return str(self.n)
            return "{}/{}".format(self.n, self.d)
        elif self.c == 2:
            return "I" if self.n else "O"
        return "{}_{}".format(self.n, self.c)

    def __int__(self):
        if self.c != 0 or self.d == 1:
            return self.n
        if self.d == -1:
            return -self.d
        return self.d*1000

    def __bool__(self):
        return bool(self.n)

    def __add__(self, other):
        if self.c == 0:
            return self.addr(other)
        if self.c == 2:
            return self.addb(other)
        return self.addp(other)

    def addr(self, other):
        assert (self.c == 0 and other.c == 0)
        n, d = self.n * other.d + other.n * self.d, self.d * other.d
        return FE((n, d), 0)

    def addb(self, other):
        assert (self.c == 2 and other.c == 2)
        return FE(self.n ^ bool(other), 2)

    def addp(self, other):
        assert (self.c == other.c)
        return FE(self.n + other.n, self.c)

    def __iadd__(self, other):
        if self.c == 0:
            return self.iaddr(other)
        if self.c == 2:
            return self.iaddb(other)
        return self.iaddp(other)

    def iaddr(self, other):
        assert (self.c == 0 and other.c == 0)
        n, d = self.n * other.d + other.n * self.d, self.d * other.d
        g = self.gcd(n, d)
        self.n, self.d = n // g, d // g
        assert (self.d > 0)
        return self

    def iaddb(self, other):
        #works if other is int, bool or FE
        self.n = self.n ^ bool(other)
        return self

    def iaddp(self, other):
        self.n = (self.n + other.n) % self.c
        return self

    def __neg__(self):
        if self.c == 0:
            return FE((-self.n, self.d), 0)
        if self.c == 2:
            return self
        return FE(self.c-self.n, self.c)

    def __pos__(self):
        return self

    def __sub__(self, other):
        if self.c == 2:
            return self.addb(other)
        return self.__add__(other.__neg__())

    def __isub__(self, other):
        if self.c == 2:
            return self.iaddb(other)
        return self.__iadd__(other.__neg__())

    def __mul__(self, other):
        if self.c == 2:
            return FE(self.n and bool(other), 2)
        if isinstance(other, int):
            if self.c == 0:
                return FE((self.n*other, self.d), 0)
            return FE(self.n*other, self.c)
        if self.c == 0:
            assert (other.c == 0)
            n, d = self.n * other.n, self.d * other.d
            return FE((n, d), 0)
        assert (self.c == other.c)
        return FE(self.n * other.n, self.c)

    def __imul__(self, other):
        if self.c == 2:
            self.n &= bool(other)
            return self
        if isinstance(other, int):
            if self.c == 0:
                g = self.gcd(self.n, other)
                self.n, self.d = self.n * other // g, self.d // g
                return self
            self.n = (self.n * other) % self.c
            return self
        if self.c == 0:
            assert (other.c == 0)
            n, d = self.n * other.n, self.d * other.d
            g = self.gcd(n,d)
            self.n, self.d = n // g, d // g
            assert (self.d > 0)
            return self
        assert (self.c == other.c)
        self.n = (self.n * other.n) % self.c
        return self

    def __invert__(self):
        if self.c == 0:
            return FE((self.d, self.n), 0)
        if self.n == 0:
            raise ZeroDivisionError
        if self.c == 2:
            return self
        if self.i:
            return FE(self.i, self.c)
        for i in range(1, self.c):
            if (self.n * i) % self.c == 1:
                self.i = i
                return FE(self.i, self.c)

    def __truediv__(self,other):
        if self.c == 2:
            return self.__mul__(other)
        return self.__mul__(other.__invert__())

    def __idiv__(self,other):
        if self.c == 2:
            return self.__imul__(other)
        return self.__imul__(other.__invert__())

    def __pow__(self,power):
        assert (isinstance(power,int))
        if self.c == 2:
            return FE(self.n or not bool(power),2)
        if self.c == 0:
            n = self.n ** power
            d = self.d ** power
            return FE((self.n ** power, self.d ** power), 0)
        return FE(self.n ** power, self.c)

    def __ipow__(self,power):
        assert (isinstance(power,int))
        if self.c == 2:
            self.n |= not bool(power)
        elif self.c == 0:
            self.n **= power
            self.d **= power
        else:
            self.n = (self.n ** power) % self.c
        return self

