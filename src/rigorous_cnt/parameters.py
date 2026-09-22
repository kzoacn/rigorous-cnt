"""Algorithm 2 parameters and the explicit boundary of its mixing guarantee."""

from dataclasses import dataclass
from functools import lru_cache
from sage.all import QQ, ZZ, RealBallField
from .numerics import rational, dyadic_lower, endpoints
from .runtime import Limits, UnresolvedBound


@dataclass(frozen=True)
class WalkParameters:
    prime_bound: int
    steps: int
    subgroup_index: int = 1
    mixing_l1_bound: object = None
    justification: str = ""
    certificate: object = None

    def __post_init__(self):
        for name in ("prime_bound","steps","subgroup_index"):
            value = getattr(self,name)
            if ZZ(value) != value:
                raise ValueError(f"{name} must be an integer")
        if self.prime_bound < 2 or self.steps < 0 or self.subgroup_index < 1:
            raise ValueError("invalid random-walk parameters")
        if self.mixing_l1_bound is not None:
            error = rational(self.mixing_l1_bound)
            if not 0 <= error <= 2 or not self.justification.strip():
                raise ValueError("a mixing bound needs an error in [0,2] and a justification")

    def check(self,epsilon,require_bound):
        if self.mixing_l1_bound is None:
            if require_bound:
                raise UnresolvedBound(
                    "Corollary 6.5's effective prime bound is not supplied. Provide "
                    "an externally justified mixing_l1_bound, or explicitly request "
                    "require_mixing_bound=False for a run without the mixing guarantee.")
            return False
        if rational(self.mixing_l1_bound) > epsilon/2:
            raise ValueError("the walk must use at most epsilon/2 of the error budget")
        return True

    @classmethod
    def from_prime_bound(cls,context,ray,prime_bound,epsilon=QQ(1)/100,
                         subgroup_index=1,mixing_l1_bound=None,justification=""):
        """Compute N from Lemma 5.1 and Corollary 6.5; B still needs justification."""
        if ray.context is not context:
            raise ValueError("ray and parameters must share a number-field context")
        steps=corollary65_steps(context.n,context.discriminant,
            ray.modulus.norm()*ZZ(2)**len(ray.real_places),epsilon,subgroup_index)
        return cls(prime_bound,steps,subgroup_index,mixing_l1_bound,justification)


def corollary65_steps(degree,discriminant,modulus_norm,epsilon,subgroup_index=1):
    """Conservative N, conditional on a sufficient prime bound B.

    Lemma 5.1 gives Vol(Pic^0_m) <= |Delta|*N(m). Use epsilon/2 for
    the walk's L1 budget. This computation requires no class or unit group.
    """
    n=ZZ(degree)
    D,M,index=map(ZZ,(discriminant,modulus_norm,subgroup_index))
    eps=rational(epsilon)
    if n < 2 or D < 1 or M < 1 or index < 1 or not 0 < eps < 1:
        raise ValueError("invalid Corollary 6.5 parameters")
    argument=QQ(D*M)/index*(2/eps)**2
    # Taking the ceiling of the upper endpoint is always conservative, even
    # when the exact logarithmic expression is very close to an integer.
    R=RealBallField(128)
    upper=endpoints(R(7*n+2)+R(argument).log())[1]
    return int(max(1,upper.ceil()))


@lru_cache(maxsize=256)
def distortion_grid(degree,discriminant,epsilon,omega=1,limits=Limits()):
    """A dyadic delta below Algorithm 2 line 4's bound (s=1/n^2).

    The proof uses an upper bound on delta, permitting downward rounding.
    """
    n = int(degree)
    eps,w = rational(epsilon),rational(omega)
    if n < 2 or not 0 < eps < min(QQ(1),QQ(n)/20) or w < 1:
        raise ValueError("Algorithm 2 requires n>=2, 0<epsilon<min(1,n/20), omega>=1")
    s = QQ(1)/(n*n)
    return dyadic_lower(lambda R: (R(eps/40)**5*R(s))/(
        R(w)**n*R(10*n*n).exp()*R(abs(discriminant))*(R(n)*R(80*n/eps).log()).sqrt()),limits)
