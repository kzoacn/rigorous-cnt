"""Exact number-field embeddings and local ray conditions."""

from dataclasses import dataclass
from sage.all import AA, QQ, QQbar, ZZ, matrix, vector
from .runtime import EmptySupport
from .numerics import rational, nearest


class NumberFieldContext:
    """Absolute number field with deterministic real/upper-half-plane places.

    The integral basis is preparatory data. Its computation is not charged to the
    paper's algorithms, whose complexity statements take an integral basis as input.
    """

    def __init__(self, field, integral_basis=None):
        if not field.is_absolute() or field.degree() < 2:
            raise ValueError("an absolute number field of degree at least two is required")
        self.K = field
        self.n = int(field.degree())
        self.discriminant = abs(ZZ(field.discriminant()))
        self.integral_basis = tuple(field.integral_basis() if integral_basis is None else integral_basis)
        given = matrix(QQ,[field(a).vector() for a in self.integral_basis]).row_module(ZZ)
        if given != field.ideal(1).free_module() or len(self.integral_basis) != self.n:
            raise ValueError("the supplied elements must generate the maximal order")
        maps = list(field.embeddings(QQbar))
        real = [e for e in maps if e(field.gen()).imag() == 0]
        complex_ = [e for e in maps if e(field.gen()).imag() > 0]
        real.sort(key=lambda e: e(field.gen()).real())
        complex_.sort(key=lambda e: (e(field.gen()).real(), e(field.gen()).imag()))
        self.places = tuple(real+complex_)
        self.r1, self.r2 = len(real), len(complex_)
        self.weights = (1,)*self.r1+(2,)*self.r2
        self.unit_ideal = field.ideal(1)

    def scales(self, values=None):
        values = [1]*len(self.places) if values is None else values
        if len(values) != len(self.places):
            raise ValueError("one scale per infinite place is required")
        out = []
        for i, v in enumerate(values):
            re, im = v if isinstance(v, (tuple, list)) else (v, 0)
            re, im = rational(re), rational(im)
            if (i < self.r1 and im) or (not re and not im):
                raise ValueError("real scales must be real; every scale must be nonzero")
            out.append((re, im))
        return tuple(out)

    def scale_norm(self, scales):
        out = QQ(1)
        for i, (re, im) in enumerate(scales):
            out *= abs(re) if i < self.r1 else re*re+im*im
        return out

    def embed(self, element, scales=None):
        scales = self.scales() if scales is None else scales
        out = []
        for i, (e, (re, im)) in enumerate(zip(self.places, scales)):
            z = e(self.K(element))*(QQbar(re)+QQbar.gen()*im)
            out.append(AA(z.real()))
            if i >= self.r1:
                out.append(AA(z.imag()))
        return vector(AA, out)

    def embedding_matrix(self, basis, scales=None):
        return matrix(AA, [self.embed(a, scales) for a in basis]).transpose()

    def weighted_norm2(self, coordinates):
        return sum(coordinates[i]**2 for i in range(self.r1)) + 2*sum(
            coordinates[i]**2 for i in range(self.r1, self.n))

    def coordinates_in_basis(self, element, basis):
        rows = matrix(QQ, [self.K(a).vector() for a in basis])
        return rows.solve_left(self.K(element).vector())


class RayConditions:
    def __init__(self, context, modulus=None, real_places=(), tau=1):
        self.context = context
        self.modulus = context.K.ideal(1 if modulus is None else modulus)
        if not self.modulus or not self.modulus.is_integral():
            raise ValueError("the finite modulus must be a nonzero integral ideal")
        self.real_places = tuple(sorted(set(int(i) for i in real_places)))
        if any(not 0 <= i < context.r1 for i in self.real_places):
            raise ValueError("ray sign constraints index real places only")
        self.tau = context.K(tau)
        if not self.tau:
            raise ValueError("tau must be nonzero")
        self.factorization = tuple(self.modulus.factor())
        if any(self.tau.valuation(p) != 0 for p, _ in self.factorization):
            raise ValueError("tau must be a local unit at every prime dividing the modulus")

    def coprime(self, ideal):
        return all(ideal.valuation(p) == 0 for p, _ in self.factorization)

    def contains(self, element):
        a = self.context.K(element)
        if not a:
            return False
        delta = a-self.tau
        if delta and any(delta.valuation(p) < e for p, e in self.factorization):
            return False
        return all((self.context.places[i](a)/self.context.places[i](self.tau)).real() > 0
                   for i in self.real_places)

    def affine_shift(self, ideal, gamma=0):
        """CRT for fractional ideals, imposing congruence locally at the modulus."""
        K = self.context.K
        gamma = K(gamma)
        if not self.coprime(ideal):
            raise ValueError("the input ideal must be coprime to the modulus")
        if gamma and any(gamma.valuation(p) < 0 for p, _ in self.factorization):
            raise EmptySupport("the affine shift is not integral at the modulus")
        if self.modulus == K.ideal(1):
            return gamma
        # J is locally O_K at the modulus, but permits the other denominators.
        J = K.ideal(1)+ideal+K.ideal(gamma)+K.ideal(self.tau)
        left, right = list(ideal.basis()), list((self.modulus*J).basis())
        A = matrix(QQ, [a.vector() for a in left+right])
        target = (self.tau-gamma).vector()
        den = ZZ(1)
        for x in list(A.list())+list(target):
            den = den.lcm(x.denominator())
        H, U = matrix(ZZ, den*A).hermite_form(transformation=True, include_zero_rows=False)
        w = H.change_ring(QQ).solve_left(den*target)
        if any(x.denominator() != 1 for x in w):
            raise EmptySupport("incompatible affine and ray conditions")
        z = vector(ZZ, w)*U
        return gamma+sum((z[i]*left[i] for i in range(len(left))), K(0))
