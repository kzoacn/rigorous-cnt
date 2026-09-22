"""Lemma 5.4 / Appendix A.2: exactly uniform bounded prime ideals."""

from sage.all import ZZ
from .runtime import Limits, RandomBits, ResourceLimit


class PrimeIdealSampler:
    def __init__(self,context,ray,bound,subgroup_oracle=None,limits=Limits()):
        self.context,self.ray,self.bound,self.limits = context,ray,ZZ(bound),limits
        self.subgroup_oracle=subgroup_oracle
        self.oracle = subgroup_oracle or (lambda ideal: True)
        if self.bound < 2:
            raise ValueError("prime bound must be at least 2")
        self._cache = {}

    def eligible_above(self,p):
        p = ZZ(p)
        if p not in self._cache:
            self._cache[p] = tuple(q for q in self.context.K.primes_above(p)
                if q.norm() <= self.bound and self.ray.coprime(q) and self.oracle(q))
        return self._cache[p]

    def sample(self,rng=None):
        rng = rng or RandomBits()
        for attempt in range(1,self.limits.attempts+1):
            p = ZZ(rng.integer(1,self.bound))
            if not p.is_prime(proof=True):
                continue
            eligible = self.eligible_above(p)
            k = len(eligible)
            if not k:
                continue
            q = eligible[rng.below(k)]
            if rng.rational_event(k,self.context.n):
                return q,attempt
        raise ResourceLimit("prime-ideal rejection budget exhausted; the eligible set may be empty")
