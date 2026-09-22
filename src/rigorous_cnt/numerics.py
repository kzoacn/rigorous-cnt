"""Certified transcendental enclosures, with exact rational endpoints."""

from sage.all import AA, QQ, ZZ, RealBallField
from .runtime import Limits, ResourceLimit


def rational(x):
    # QQ(float) is a rational reconstruction, not the exact requested input.
    if isinstance(x, float) or type(x).__name__ in ("RealNumber", "RealDoubleElement"):
        raise TypeError("use an exact integer, rational, or rational string")
    return QQ(x)


def endpoints(ball):
    # QQ(ball.lower()) can reconstruct a rational OUTSIDE the interval.
    return ball.lower().exact_rational(), ball.upper().exact_rational()


def precisions(limits=Limits()):
    p = 64
    while p <= limits.precision:
        yield p
        p *= 2
    raise ResourceLimit("certified numerical comparison exceeded the precision limit")


def nearest(x):
    """Translation-equivariant rounding, including half-integer ties."""
    return ZZ((x + QQ(1)/2).floor())


def dyadic_lower(evaluate, limits=Limits()):
    """Return 2**(-k) below a positive real given by a ball evaluator."""
    for p in precisions(limits):
        lo, hi = endpoints(evaluate(RealBallField(p)))
        if lo > 0:
            k = max(0, int(lo.denominator().nbits() - lo.numerator().nbits() + 1))
            q = QQ(1) / ZZ(2)**k
            while q > lo:
                q /= 2
            return q


def bernoulli_exp_pi(q, rng, limits=Limits()):
    """Exact Bernoulli(exp(-pi*q)), almost surely, from a lazy uniform real."""
    q = rational(q)
    if q < 0:
        raise ValueError("negative exponent parameter")
    if q == 0:
        return True
    prefix, width = 0, 0
    for p in precisions(limits):
        prefix = (prefix << 32) + rng.bits(32)
        width += 32
        ulo, uhi = QQ(prefix)/ZZ(2)**width, QQ(prefix+1)/ZZ(2)**width
        R = RealBallField(p)
        plo, phi = endpoints((-R.pi()*R(q)).exp())
        if uhi <= plo:
            return True
        if ulo >= phi:
            return False


def gaussian_integer_bounds(center, variance_parameter, log_argument, limits=Limits()):
    """Integers in c +/- sqrt(variance_parameter * log(log_argument))."""
    center, variance_parameter, log_argument = map(rational, (center, variance_parameter, log_argument))
    if variance_parameter <= 0 or log_argument <= 1:
        raise ValueError("invalid Gaussian truncation")
    for p in precisions(limits):
        R = RealBallField(p)
        radius = (R(variance_parameter)*R(log_argument).log()).sqrt()
        alo, ahi = endpoints(R(center)-radius)
        blo, bhi = endpoints(R(center)+radius)
        if alo.ceil() == ahi.ceil() and blo.floor() == bhi.floor():
            return ZZ(alo.ceil()), ZZ(blo.floor())


def norm_one_exp_approx(log_coordinates, weights, relative_error, limits=Limits()):
    """Positive rational A with product A_i**w_i=1, close to exp(a_i/w_i)."""
    a = tuple(map(rational, log_coordinates))
    w = tuple(map(int, weights))
    eta = rational(relative_error)
    if len(a) != len(w) or not a or sum(a) != 0 or any(x not in (1, 2) for x in w):
        raise ValueError("expected zero-sum log coordinates and real/complex weights")
    if not 0 < eta < 1:
        raise ValueError("relative error must lie in (0,1)")
    last = w.index(1) if 1 in w else len(w)-1
    for p in precisions(limits):
        R = RealBallField(p)
        target = [(R(ai)/wi).exp() for ai, wi in zip(a, w)]
        values = [v.mid().exact_rational() for v in target]
        if min(values) <= 0:
            continue
        product = QQ(1)
        for i in range(len(a)):
            if i != last:
                product *= values[i] ** (w[i] // w[last])
        values[last] = 1/product
        if all(endpoints(abs(R(v)/t-1))[1] <= eta for v, t in zip(values, target)):
            return tuple(values)
