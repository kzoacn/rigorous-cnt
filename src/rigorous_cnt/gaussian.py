"""Lemma 2.22: a GPV sampler on rational lattices with explicit tail/error data."""

from dataclasses import dataclass
from sage.all import QQ, ZZ, matrix, vector, RealBallField
from .runtime import Limits, RandomBits, ResourceLimit
from .numerics import rational, precisions, endpoints, bernoulli_exp_pi, gaussian_integer_bounds


@dataclass(frozen=True)
class GaussianSample:
    value: object
    coefficients: tuple
    statistical_distance_bound: object
    squared_radius_log_argument: object
    attempts: int


def rational_gso(B):
    orthogonal = []
    for row in B.rows():
        u = vector(QQ, row)
        for b in orthogonal:
            u -= (row.dot_product(b)/b.dot_product(b))*b
        if not u:
            raise ValueError("linearly dependent lattice basis")
        orthogonal.append(u)
    return orthogonal


def discrete_gaussian(B, s, epsilon, center=None, rng=None, limits=Limits()):
    """Sample close to density exp(-pi*||x-c||^2/s^2).

    Rows of B are a rational basis. The width hypothesis is checked, not assumed.
    One-dimensional rejection uses exact Bernoulli decisions, so it adds no
    floating-point probability error. Truncation follows Lemma 2.22.
    """
    B = matrix(QQ, B)
    n, d = B.nrows(), B.ncols()
    s, epsilon = rational(s), rational(epsilon)
    if s <= 0 or not 0 < epsilon <= 1:
        raise ValueError("s>0 and 0<epsilon<=1 required")
    rng = rng or RandomBits()
    c = vector(QQ, [0]*d if center is None else center)
    if len(c) != d:
        raise ValueError("center dimension mismatch")
    if not n:
        if c:
            raise ValueError("center is outside the zero lattice span")
        return GaussianSample(c, (), QQ(0), QQ(1), 0)
    orthogonal = rational_gso(B)
    if c not in B.row_space():
        raise ValueError("center must lie in the rational span of B")
    max_norm2 = max(v.dot_product(v) for v in B.rows())
    for p in precisions(limits):
        R = RealBallField(p)
        required = (R(1/epsilon).log()+2*R(n).log()+3)/R.pi()*R(max_norm2)
        lo, hi = endpoints(required)
        if s*s >= hi:
            break
        if s*s < lo:
            raise ValueError("Gaussian width violates the hypothesis of Lemma 2.22")
    log_arg = QQ(2*n*n)/epsilon
    residual = vector(QQ, c)
    coeff = [ZZ(0)]*n
    attempts = 0
    for i in reversed(range(n)):
        norm2 = orthogonal[i].dot_product(orthogonal[i])
        ci = residual.dot_product(orthogonal[i])/norm2
        width2 = s*s/norm2
        lo, hi = gaussian_integer_bounds(ci, width2, log_arg, limits)
        while True:
            attempts += 1
            if attempts > limits.attempts:
                raise ResourceLimit("discrete Gaussian rejection budget exhausted")
            z = ZZ(rng.integer(lo, hi))
            if bernoulli_exp_pi((z-ci)**2/width2, rng, limits):
                break
        coeff[i] = z
        residual -= z*B.row(i)
    out = vector(ZZ, coeff)*B
    return GaussianSample(out, tuple(coeff), epsilon, log_arg, attempts)
