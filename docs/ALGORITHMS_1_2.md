# Algorithms 1–2: sampling guide

Reference implementations of Algorithm 1 (`Sample in a Box`) and Algorithm 2 (`Sample in an Ideal Class`) from [arXiv:2512.01588v2](https://arxiv.org/abs/2512.01588v2). The code runs as ordinary Python inside SageMath and has been tested with SageMath 10.7 and Python 3.12.

Supported inputs include absolute number fields with real, complex, or mixed signatures; integral and fractional ideals; rational shifts; finite moduli; sign conditions at real embeddings; rational real or complex scaling at each infinite place; and a caller-supplied membership oracle for a subgroup of finite index.

## Supported behavior and guarantees

- Algorithm 1 implements exact uniform sampling and retains the output symbolically as `x * alpha`. Ideal membership, local congruences, strict sign conditions, and geometric boundaries are checked exactly.
- Algorithm 2 implements the random prime-ideal product, the discrete Gaussian of Lemma 2.22, rational distortion, and the call to Algorithm 1. It returns an exact field element `beta`, intermediate data, and a report.
- For small fields, `calibrate_walk(...)` can choose and certify `B,N`. It explicitly computes class and unit groups, bounds the eigenvalues of finitely many characters, and combines them with the Gaussian tail bound of Appendix A.1 to certify the mixing error. This separate preprocessing supports the full Arakelov ray class group; general subgroups still need an external mixing bound.
- **The efficient general parameter bound of Corollary 6.5 has not been fully instantiated.** `calibrate_walk` is an expensive direct calibration, not an arbitrary assignment to hidden constants. With uncalibrated `B,N`, the default interface refuses to claim a mixing guarantee; explicitly disabling that requirement permits experiments. `WalkParameters.from_prime_bound(...)` only chooses `N` from a volume upper bound.
- **The paper's full bit-complexity bounds are not claimed.** The implementation first proposes an LLL transformation using an integer approximation and certifies the lengths of the resulting exact basis vectors. It can fall back to exact LLL and block HKZ enumeration. The enumeration backend is not a complete implementation of the Kannan subroutine used in the complexity proof.

Probability statements use the model of independent uniform random bits. Fixed seeds support reproducible experiments. Resource exhaustion raises an exception instead of returning a fabricated sample; an uninterrupted distribution must be distinguished from a truncated execution. See [MATHEMATICS.md](MATHEMATICS.md) for mathematical details, implementation differences, and outstanding theoretical guarantees.

## Run

From the repository root, without modifying the existing Sage environment:

```bash
PYTHONPATH=src sage -python examples/demo.py --calibrate
PYTHONPATH=src sage -python examples/demo.py --field imaginary --calibrate
PYTHONPATH=src sage -python examples/demo.py --field mixed --calibrate
PYTHONPATH=src sage -python -m unittest discover -s tests -v
```

Alternatively, install the package in the Sage environment:

```bash
sage -pip install -e . --no-build-isolation
```

## Algorithm 1

```python
from sage.all import QQ, PolynomialRing, NumberField
from rigorous_cnt import NumberFieldContext, RayConditions, RandomBits, Algorithm1

R = PolynomialRing(QQ, 't')
t = R.gen()
K = NumberField(t*t - 2, 'a')
ctx = NumberFieldContext(K)

# Real embeddings are ordered by the image of the generator.
# Require a negative value at real embedding 0.
ray = RayConditions(ctx, modulus=K.ideal(2), real_places=[0], tau=-1)
sampler = Algorithm1(ctx, K.ideal(3), ray=ray, block_size=2)
sample = sampler.sample(RandomBits(seed=7))
assert sampler.verify(sample.element)
print(sample.element)       # Exact alpha
print(sample.point(ctx))    # Algorithm 1's x*alpha
print(sample.report)
```

Reuse an `Algorithm1` object to sample repeatedly with the same prepared basis and grid. Supply general scaling through `x`, with one entry per infinite place: an exact rational at a real place, or `(real_part, imaginary_part)` at a complex place. From each complex conjugate pair, the embedding with positive imaginary part is selected. Full real coordinates are ordered as real places followed by the real and imaginary parts at complex places.

`gamma` may be fractional. If it is nonintegral at a prime dividing the finite modulus, the support is empty and `EmptySupport` is raised. Congruence is a local valuation condition; membership of the difference in an integral modulus ideal is not an equivalent replacement for all fractional inputs.

## Algorithm 2 with automatic calibration

```python
from rigorous_cnt import algorithm2, calibrate_walk

# Run the explicit, more expensive preprocessing; its result is reusable.
walk = calibrate_walk(ctx, ray, epsilon=QQ(1)/100)
result = algorithm2(ctx, K.ideal(3), walk, ray=ray, rng=RandomBits(seed=7))
assert result.report['mixing_guarantee'] == 'verified_spectral_bound'
print(walk.prime_bound, walk.steps, walk.mixing_l1_bound)
```

The calibrator constructs ray units and characters from PARI data certified with `proof=True`. Directed ball arithmetic covers all relevant low-frequency characters and bounds their eigenvalues. Explicit budgets limit group size, frequency count, and the prime-ideal bound; exceeding a budget raises `ResourceLimit`. A certificate is bound to its number-field context, modulus, step count, and norm bound, and cannot be reused for different parameters.

## Algorithm 2 with uncalibrated parameters

```python
from rigorous_cnt import algorithm2, WalkParameters

result = algorithm2(
    ctx, K.ideal(3),
    walk=WalkParameters.from_prime_bound(ctx, ray, prime_bound=19),
    ray=ray, epsilon=QQ(1)/100,
    rng=RandomBits(seed=7),
    require_mixing_bound=False,  # Explicitly allow B without a mixing guarantee
)
assert result.element in K.ideal(3)
assert ray.contains(result.element)
assert K.ideal(result.element) == result.output_ideal() * result.input_ideal
print(result.report['mixing_guarantee'])  # not_established
```

`result.output_ideal()` is `(beta) * input_ideal^(-1)` and lies in the inverse of the input ideal class. It is a candidate ideal: an individual output need not be smooth or near-prime. The application checks membership in its target family of ideals.

For parameters derived externally, supply `WalkParameters(..., mixing_l1_bound=..., justification=...)`, with error at most `epsilon/2`. The report uses `conditional_on_supplied_bound`: the caller's mathematical justification is a premise, and descriptive text alone is not treated as a verified proof.

A general subgroup also requires `subgroup_index` and `subgroup_oracle(ideal)`. Its index, the oracle's homomorphic interpretation, and the corresponding closed subgroup of finite index are caller-supplied assumptions. The program checks that sampled prime ideals satisfy the oracle and checks the condition associated with `tau`.

## Validation and implementation map

Tests enumerate proposal masses on small grids; check uniform prime-ideal masses for splitting, inert, and ramified primes; and cover lattice transformations, HKZ conditions, Gaussian support, exact ball endpoints, rational distortion normalization, fractional CRT, complex embeddings, replay, resource failures, and rejection of unjustified mixing claims. Tests also forbid `pari_bnf` in the core sampling path. Class and unit computations inside `calibrate_walk` are explicit, separate preprocessing.

| Files | Purpose |
| --- | --- |
| `arithmetic.py` | Number-field context, embeddings, local ray conditions, fractional CRT |
| `geometry.py` | Proposition 8.8 and grid sampling in real intervals and complex disks |
| `gaussian.py`, `numerics.py` | Lemma 2.22, exact random decisions, bounded transcendental evaluation |
| `lattice.py` | Exact LLL, HKZ, and construction and checking of the required short basis |
| `primes.py` | Lemma 5.4 and Appendix A.2 |
| `parameters.py` | Rational grids, walk parameters, and unverified assumptions |
| `calibration.py` | Explicit character and Gaussian-tail calibration for small fields |
| `sampling.py` | Algorithms 1–2 and their result objects |

This guide covers Algorithms 1–2; see the [main documentation](../README.md) for Algorithms 3–6. Nonmaximal orders are not supported. Convert relative fields to an absolute representation first. The underlying field, maximal-order, and exact algebraic-number operations use SageMath; their preprocessing cost should be distinguished from the paper's core algorithmic cost.
