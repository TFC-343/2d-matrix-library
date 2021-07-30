# 2d matrix library
# by Ariel (tfc-343)

from copy import copy


class matrix:
    """the storage and calculating of 2d matrices"""

    class BaseError(Exception):
        """for errors where the dimensions of the matrix will not work with the matrix"""
        pass

    class TypeError(Exception):
        """for errors where the data type is incorrect"""
        pass

    class InvertError(Exception):
        pass

    def __init__(self, m, n, *data):
        """
        :param m: numbers down
        :param n: number left
        :param data: all the values in the matrix
        """
        if not(isinstance(m, int) and isinstance(n, int)) or (m < 1 or n < 1):
            raise self.TypeError("cannot create matrix with those dimensions")
        self.m = m
        self.n = n
        self.data = data

    def __str__(self):
        """will print the data when called"""
        return str(tuple([self.m, self.n, *self.data]))

    def __add__(self, other):
        """will add 2 matrices together"""
        m = self.m
        n = self.n

        if not(isinstance(self, matrix)) or not(isinstance(other, matrix)):
            # ^^^ if the user tried to add a matrix and a type that isn't a matrix an error will be called
            raise self.TypeError("Must add matrix to another matrix")

        if (n != other.n) or (m != other.m):  # check that both bases are equal, if not, raise an error
            raise self.BaseError("Base must be the same when adding")

        t = ()
        for i in range(m*n):
            t += tuple([self.get(i // m + 1, i % n + 1) + other.get(i // m + 1, i % n + 1)])
        return matrix(m, n, *t)

    def __radd__(self, other):
        """refer to __add__"""
        self.__add__(other)

    def __sub__(self, other):
        """will subtract matrices from each other"""
        m = self.m
        n = self.n

        if not(isinstance(self, matrix)) or not(isinstance(other, matrix)):
            raise self.TypeError("Must subtract matrix from another matrix")

        if (n != other.n) or (m != other.m):  # check that both bases are equal, if not, raise an error
            raise self.BaseError("Base must be the same when subtracting")

        return self+other*-1

    def __rsub__(self, other):
        """refer to __sub__"""
        self.__sub__(other)

    def __mul__(self, other):
        """will multiply matrices"""
        if isinstance(other, int) or isinstance(other, float):  # a matrix times a scalar
            t = tuple([i*other for i in self.data])
            return matrix(self.m, self.n, *t)

        if isinstance(other, matrix):  # a matrix times a matrix
            if not(self.n == other.m):
                raise self.BaseError("Bases do not work for multiplication")

            t = ()
            for i in range(self.m):
                for j in range(other.n):
                    o = 0
                    for k in range(self.n):
                        p = self.get(i+1, k+1) * other.get(k+1, j+1)
                        o += p
                    t += tuple([o])

            return matrix(self.m, other.n, *t)

        raise self.TypeError("type error")

    def __rmul__(self, other):
        """refer to __mul__"""
        return self.__mul__(other)

    def __truediv__(self, other):
        """multiplies 'self' by the inverse of 'other'"""
        return self * (other**-1)

    def __rtruediv__(self, other):
        """inverse of 'self' times 'other'"""
        return (self**-1) * other

    def __pow__(self, power):
        """'self' raised to the power of 'power'"""
        m = self.m
        n = self.n

        if n != m:
            raise self.BaseError("Must be a square matrix")

        if not(isinstance(power, int)):
            raise self.TypeError("power must be an int")

        def exp_by_sqr(x, p):
            if p == -1:
                return x.get_inverse()
            elif p == 0:
                return matrix.identity_matrix(m)
            elif p == 1:
                return x
            elif p % 2 == 0:
                return exp_by_sqr(x * x, p // 2)
            elif p % 2 != 0:
                return x * exp_by_sqr(x * x, (p - 1) // 2)
        return exp_by_sqr(self, power)

    def __abs__(self):
        return self.get_determinant()

    def __len__(self):
        """returns the amount of points in a matrix"""
        return len(self.data)

    def __round__(self, r=None):
        """returns the rounded points in matrix"""
        m = self.m
        n = self.n
        new_data = list(copy(self.data))
        for i in range(m * n):
            new_data[i] = round(new_data[i], r)
        return matrix(m, n, *tuple(new_data))

    def __eq__(self, other):
        """returns True if all points in a matrix are equal"""
        return self.m == other.m and self.n == other.n and self.data == other.data

    def __invert__(self):
        """refer to get_inverse()"""
        return self.get_inverse()

    def __pos__(self):
        """positive matrix is just a matrix, so returns original matrix"""
        return self

    def __neg__(self):
        """returns matrix * -1"""
        return self*-1

    def set(self, m, n, new):
        """set data in matrix to new value based on coords"""
        n_base = self.n
        data = self.data

        c = data[(m - 1) * n_base + (n - 1)]

        self.data = data[:(m - 1) * n_base + (n - 1)] + tuple([new]) + data[(m - 1) * n_base + (n - 1) + 1:]

        return c

    def fill(self):
        """fills a matrix with zeroes, or deletes overflow"""
        m = self.m
        n = self.n
        data = self.data

        if len(data) == n * m:  # does nothing
            c = 0
        elif len(data) > n * m:  # removes extra data
            self.data = data[0: n * m]
            c = -len(data) + m * n
        else:  # adds zeroes if not filled
            self.data = data + tuple(0 for _ in range(n * m - len(data)))
            c = n * m - len(data)

        return c

    def flip(self, o='n'):
        """flips a matrix along an axis"""
        m = self.m
        n = self.n
        data = self.data

        if o == 'n':  # flips a matrix horizontally (left becomes right)
            t = ()
            for i in range(m):
                t += data[(m-i-1)*n: (m-i-1)*n+n]
            self.data = t

        if o == 'm':  # flips a matrix vertically (top becomes bottom)
            t = ()
            for i in range(m):
                t += data[i*n: i*n+n][::-1]
            self.data = t

    def transpose(self):
        """transposes a matrix"""
        m = self.m
        n = self.n
        get = self.get

        new = matrix(n, m)
        new.fill()
        for i in range(n):
            for j in range(m):
                new.set(i+1, j+1, get(j+1, i+1))
        self.m = n
        self.n = m
        self.data = new.data

    def output(self):
        """outputs the matrix"""
        m = self.m
        n = self.n
        data = self.data

        for i in range(m):
            for j in data[i*n:i*n+n]:
                print(j, end=" ")
            print()

    def get_list(self, p='tuple'):
        """returns matrix in list form"""
        m_base = self.m
        n_base = self.n
        data = self.data

        if p == 'list':  # returns the matrix as an n dimensional list
            o = []
            for i in range(m_base):
                p = []
                for j in data[i * n_base:i * n_base + n_base]:
                    p.append(j)
                o.append(p)
            return o

        if p == 'tuple':  # return the matrix as a tuple
            return data

    def get(self, m, n):
        """returns data in matrix based on coords"""
        n_base = self.n
        data = self.data

        return data[(m - 1) * n_base + (n - 1)]  # gets an int from the tuple based on coords

    def get_determinant(self):
        """will find the determinant of a matrix"""
        m = self.m
        n = self.n
        get = self.get
        data = self.data

        if m != n:
            raise self.BaseError("must be square matrix")
        elif m == 1:
            return data[0]
        elif m == 2:
            return get(1, 1) * get(2, 2) - get(2, 1) * get(1, 2)
        else:
            rt = 0
            for i in range(n):
                t = ()
                for j in range(n, m*n):
                    if (j+1 + n-(i+1)) % n != 0:
                        t += tuple([data[j]])
                rt += (-1)**i * data[i] * abs(matrix(m - 1, n - 1, *t))
            return rt

    def get_inverse(self):
        """returns the inverse"""
        n = self.n
        get = self.get
        get_cofactors = self.get_cofactors
        get_determinant = self.get_determinant

        print(get_determinant())

        if get_determinant() == 0:
            text = "The determinant is 0, the matrix is not invertible"
            raise self.InvertError(text)

        if n == 1:
            return matrix(1, 1, 1/get(1, 1))

        else:
            new_matrix = get_cofactors()
            transposed = copy(new_matrix)
            transposed.transpose()
            return get_determinant()**-1 * transposed

    def get_minors(self):
        """will return a matrix of minors"""
        m = self.m
        n = self.n
        get = self.get

        if m != n:
            raise self.BaseError("must be square matrix")

        if m == 1:
            return self

        if m == 2:
            return matrix(2, 2, get(2, 2), get(2, 1), get(1, 2), get(1, 1))

        new_matrix = copy(self)
        for i in range(m ** 2):  # cycle though the matrix to find the matrix of minors of each coord
            t = ()
            for j in range(m ** 2):
                if j // m + 1 == i // m + 1 or j % n + 1 == i % n + 1:  # creating the minor matrix
                    pass
                else:
                    t += tuple([get(j // m + 1, j % n + 1)])
            minor_matrix = matrix(m - 1, n - 1, *t)
            new_matrix.set(i // m + 1, i % n + 1, minor_matrix.get_determinant())
        return new_matrix

    def get_cofactors(self):
        """will return the cofactors of the minors"""
        m = self.m
        n = self.n

        if m != n:
            raise self.BaseError("must be square matrix")

        new_matrix = copy(self)
        new_matrix = new_matrix.get_minors()
        for i in range(m**2):
            if ((i // m + 1) + (i % n + 1)) % 2 != 0:
                new_matrix.set(i // m + 1, i % n + 1,
                               -new_matrix.get(i // m + 1, i % n + 1))

        return new_matrix

    @staticmethod
    def identity_matrix(m):
        """returns an identity matrix of m size"""
        t = tuple([1])
        for _ in range(m):
            t = t + tuple(0 for _ in range(m)) + tuple([1])
        return matrix(m, m, *t)

    def test(self):
        print(self.data)
