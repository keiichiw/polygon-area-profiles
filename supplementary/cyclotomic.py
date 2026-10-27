"""Exact Q(zeta_N) arithmetic; integer coefficient vectors/common denominator.

No floating-point decisions occur here. Numerical embeddings are used only by
explore_cells.py to order distinct intersections and incident edges.
"""
from fractions import Fraction
from functools import reduce
from math import gcd, lcm
from verify_myerson import cyclotomic


def trim(a):
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def divide(a, b):
    a = list(map(Fraction, a))
    b = trim(list(map(Fraction, b)))
    q = [Fraction(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        k, c = len(a) - len(b), a[-1] / b[-1]
        q[k] = c
        for j, x in enumerate(b):
            a[k+j] -= c*x
        trim(a)
    return trim(q), trim(a)


def poly_sub_product(a, q, b):
    out = list(a) + [Fraction(0)] * max(0, len(q)+len(b)-1-len(a))
    for i, x in enumerate(q):
        for j, y in enumerate(b):
            out[i+j] -= x*y
    return trim(out)


class Field:
    def __init__(self, order):
        self.order = order
        self.modulus = cyclotomic(order)
        self.degree = len(self.modulus)-1
        self.reduction = [(i, -c) for i, c in enumerate(self.modulus[:-1]) if c]
        self.zero = Element(self, [0])
        self.one = Element(self, [1])
        self.powers = []
        for j in range(order):
            self.powers.append(Element(self, [0]*j+[1]))
        self.conjugate_powers = [tuple((k,c) for k,c in enumerate(self.root(-j).coeff) if c)
                                 for j in range(self.degree)]
        self.binomial_inverses = {}

    def root(self, exponent):
        return self.powers[exponent % self.order]

    def binomial_inverse(self, exponent):
        r = exponent % self.order
        if not r:
            raise ZeroDivisionError('parallel/coincident chord lines')
        if r not in self.binomial_inverses:
            self.binomial_inverses[r] = (self.one-self.root(r)).inverse()
        return self.binomial_inverses[r]

    def intersection(self, a, b, c, d):
        # Endpoints are powers of this field's primitive root.
        # z + zeta^(a+b) conjugate(z) = zeta^a + zeta^b.
        r = a+b-c-d
        numerator = self.root(a)+self.root(b)-self.root(r+c)-self.root(r+d)
        return numerator*self.binomial_inverse(r)


class Element:
    __slots__ = ('field', 'coeff', 'den')

    def __init__(self, field, coeff, den=1):
        if not den:
            raise ZeroDivisionError
        a = list(coeff)
        # Phi_n divides x^n-1. Fold first; for prime n this leaves just
        # one coefficient above degree Phi_n, instead of quadratic division.
        if len(a)>field.order:
            for k in range(field.order,len(a)):
                a[k%field.order] += a[k]
            a = a[:field.order]
        for k in range(len(a)-1, field.degree-1, -1):
            c = a[k]
            if c:
                for j, v in field.reduction:
                    a[k-field.degree+j] += c*v
        a = (a[:field.degree]+[0]*field.degree)[:field.degree]
        g = reduce(gcd, a, den)
        if den < 0:
            g = -g
        self.field, self.coeff, self.den = field, tuple(c//g for c in a), den//g

    def scalar(self, value):
        if isinstance(value, Element):
            if value.field is not self.field:
                raise ValueError('different fields')
            return value
        q = Fraction(value)
        return Element(self.field, [q.numerator], q.denominator)

    def __add__(self, other):
        other = self.scalar(other)
        g = gcd(self.den, other.den)
        a, b = other.den//g, self.den//g
        return Element(self.field, [x*a+y*b for x,y in zip(self.coeff, other.coeff)], self.den*a)

    __radd__ = __add__

    def __neg__(self):
        return Element(self.field, [-x for x in self.coeff], self.den)

    def __sub__(self, other):
        return self + (-self.scalar(other))

    def __rsub__(self, other):
        return self.scalar(other)-self

    def __mul__(self, other):
        other = self.scalar(other)
        a = [0]*(2*self.field.degree-1)
        left = [(i,x) for i,x in enumerate(self.coeff) if x]
        right = [(i,x) for i,x in enumerate(other.coeff) if x]
        if len(left)*len(right)>256:
            # Exact Kronecker substitution. Balanced base digits recover the
            # signed convolution; the bound prevents any ambiguous carry.
            bound = min(len(left),len(right))*max(abs(x) for i,x in left)*max(abs(y) for j,y in right)
            bits = bound.bit_length()+2
            base, mask = 1<<bits, (1<<bits)-1
            first = second = 0
            for x in reversed(self.coeff):
                first = (first<<bits)+x
            for y in reversed(other.coeff):
                second = (second<<bits)+y
            product = first*second
            for k in range(len(a)):
                digit = product&mask
                if digit>=base//2:
                    digit -= base
                a[k] = digit
                product = (product-digit)>>bits
            if product:
                raise ArithmeticError('polynomial product carry bound failed')
        else:
            for i,x in left:
                for j,y in right:
                    a[i+j] += x*y
        return Element(self.field, a, self.den*other.den)

    __rmul__ = __mul__

    def inverse(self):
        if not self:
            raise ZeroDivisionError
        r0, r1 = list(self.field.modulus), trim(list(self.coeff))
        t0, t1 = [Fraction(0)], [Fraction(1)]
        while any(r1):
            q, r = divide(r0, r1)
            r0, r1 = r1, r
            t0, t1 = t1, poly_sub_product(t0, q, t1)
        if len(r0) != 1:
            raise ArithmeticError('noninvertible field element')
        values = [x*self.den/r0[0] for x in t0]
        den = lcm(*(x.denominator for x in values))
        return Element(self.field, [int(x*den) for x in values], den)

    def __truediv__(self, other):
        return self*self.scalar(other).inverse()

    def conjugate(self):
        coeff = [0]*self.field.degree
        for j, c in enumerate(self.coeff):
            if c:
                for k, x in self.field.conjugate_powers[j]:
                    coeff[k] += c*x
        return Element(self.field, coeff, self.den)

    def rotate(self, exponent):
        """Multiply by zeta^exponent using a cyclic coefficient shift."""
        n = self.field.order
        coeff = [0]*n
        for j,c in enumerate(self.coeff):
            coeff[(j+exponent)%n] = c
        return Element(self.field, coeff, self.den)

    def rational(self):
        return Fraction(self.coeff[0], self.den) if not any(self.coeff[1:]) else None

    def rational_ratio(self, other):
        """Return q in Q with self=q*other, or None (other must be nonzero)."""
        other = self.scalar(other)
        j = next((j for j,x in enumerate(other.coeff) if x), None)
        if j is None:
            raise ZeroDivisionError
        a, b = self.coeff[j], other.coeff[j]
        if all(x*b == y*a for x,y in zip(self.coeff, other.coeff)):
            return Fraction(a*other.den, b*self.den)
        return None

    def __bool__(self):
        return any(self.coeff)

    def __eq__(self, other):
        other = self.scalar(other)
        return self.coeff == other.coeff and self.den == other.den

    def __hash__(self):
        return hash((id(self.field), self.coeff, self.den))

    def key(self):
        return self.coeff, self.den


def wedge(a, b):
    """Twice the oriented area times 2i (the factor cancels in area ratios)."""
    return a.conjugate()*b-a*b.conjugate()


def profile(n, k, m, field=None):
    """Independent geometry: four boundary arcs followed by the intersection."""
    f = field or Field(2*n)
    if f.order != 2*n or not (1 <= k < m < n/2):
        raise ValueError('inadmissible configuration')
    end = n-m+k
    x = f.intersection(0,2*m,2*k,2*end)
    step = f.root(2)-f.root(-2)
    # Use the known binomial denominator instead of repeated Euclidean inversion.
    normalizer = -f.root(2)*f.binomial_inverse(4)*Fraction(1,n)
    areas = []
    for start, stop in ((0,k),(k,m),(m,end),(end,n)):
        doubled = (stop-start)*step+wedge(f.root(2*stop),x)+wedge(x,f.root(2*start))
        areas.append(doubled*normalizer)
    if sum(areas) != 1 or areas[1] != areas[3]:
        raise ArithmeticError('area sum or reflection failed')
    return areas
