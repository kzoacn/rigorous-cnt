"""Explicit resource failures and random-bit access (never floating-point RNG)."""

from dataclasses import dataclass
import random


class ResourceLimit(RuntimeError):
    """Computation was interrupted; no mathematical output is being certified."""


class EmptySupport(ValueError):
    """The requested sampling set is empty."""


class UnresolvedBound(ValueError):
    """A theorem hypothesis has not been supplied or verified."""


@dataclass(frozen=True)
class Limits:
    attempts: int = 100000
    precision: int = 32768
    lattice_steps: int = 10000
    enumeration_nodes: int = 1000000

    def __post_init__(self):
        if min(self.attempts, self.lattice_steps, self.enumeration_nodes) < 1:
            raise ValueError("resource limits must be positive")
        if self.precision < 64:
            raise ValueError("precision limit must be at least 64 bits")


class RandomBits:
    """A uniform-bit oracle interface; a seed enables deterministic replay.

    The mathematical analysis assumes independent unbiased bits. Python's seeded
    PRNG is an experimental realization of that interface, not such a proof.
    """

    def __init__(self, seed=None, source=None):
        if source is not None and seed is not None:
            raise ValueError("choose a source or a seed")
        self.source = source or (random.SystemRandom() if seed is None else random.Random(seed))
        self.seed = seed
        self.bits_used = 0

    def bits(self, n):
        n = int(n)
        self.bits_used += n
        return self.source.getrandbits(n)

    def below(self, n):
        n = int(n)
        if n < 1:
            raise ValueError("empty integer range")
        if n == 1:
            return 0
        k = (n - 1).bit_length()
        while True:
            x = self.bits(k)
            if x < n:
                return x

    def integer(self, lo, hi):
        lo, hi = int(lo), int(hi)
        if hi < lo:
            raise EmptySupport("empty integer interval")
        return lo + self.below(hi - lo + 1)

    def rational_event(self, numerator, denominator):
        numerator, denominator = int(numerator), int(denominator)
        if not 0 <= numerator <= denominator or denominator < 1:
            raise ValueError("invalid probability")
        return self.below(denominator) < numerator
