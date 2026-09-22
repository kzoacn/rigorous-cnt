# rigorous-cnt

A SageMath reference implementation of Algorithms 1–6 in [*Rigorous methods for computational number theory*, arXiv:2512.01588v2](https://arxiv.org/abs/2512.01588v2), by Koen de Boer, Alice Pellet-Mary, and Benjamin Wesolowski.

The accompanying implementation and validation paper by Kaiyi Zhang is published as [Markdownxiv 2609.00011v1](https://markdownxiv.github.io/abs/2609.00011v1/). [Markdown source](paper/manuscript.md) · [Publication record](docs/PUBLICATION.md)

The implementation connects ideal sampling, smooth relations, S-unit observations, exceptional units, lattice postprocessing, and restriction to an arbitrary set of finite primes. It retains units as compact power products and exports exact data that can be checked again without repeating the random search.

## Guarantees and scope

Ideal identities, valuations, and integer transformations use exact arithmetic. Real quantities use directed ball bounds where certification is required. The complete-generation test combines an independently bounded Dedekind zeta residue with a lattice-index test and is conditional on GRH for the zeta functions in the Belabas–Friedman bound.

Default sampling parameters are intended for reference experiments: their general mixing and running-time guarantees have not been established. A result is returned only after the separate completeness test succeeds. Exhausting a resource budget raises an error. This implementation does not establish the original paper's full bit-complexity bounds or a performance advantage over SageMath.

The default solving path does not call an existing class-group, unit-group, or S-unit-group solver. Optional spectral calibration explicitly uses existing class and unit data; its preprocessing should be accounted for separately. SageMath supplies the underlying number-field, ideal, lattice, and ball-arithmetic operations.

## Run

Tested with SageMath 10.7 and Python 3.12. Ordinary Python without SageMath is insufficient. From the repository root:

```bash
PYTHONPATH=src sage -python -m unittest discover -s tests -v
PYTHONPATH=src sage -python examples/demo_sunits.py --field real
PYTHONPATH=src sage -python examples/demo_sunits.py --field imaginary --primes 3 --seed 7
PYTHONPATH=src sage -python examples/demo.py --calibrate
```

The S-unit examples report relation counts, current rank, and possible lattice indices. More relations may be needed after full rank is first attained. Retry counts depend on the field, factor base, and seed; examples commonly take tens of seconds to several minutes.

Alternatively install in the Sage environment with `sage -pip install -e . --no-build-isolation`.

The saved results are quicker to check than to regenerate:

```bash
PYTHONPATH=src sage -python examples/verify_sunits.py examples/results/real_units.json
PYTHONPATH=src sage -python examples/verify_sunits.py examples/results/class2_sunits.json
```

The verifier recomputes the residue interval, the complete working lattice, and the restriction to the requested primes. It does not trust success flags or volume claims stored in the JSON. The test suite also squares a correct generator and verifies that the resulting full-rank subgroup is rejected with index two.

## API example

```python
from sage.all import QQ, PolynomialRing, NumberField
from rigorous_cnt import NumberFieldContext, RandomBits, SUnitSettings, algorithm6

R = PolynomialRing(QQ, "t")
t = R.gen()
K = NumberField(t*t - 2, "a")
result = algorithm6(
    NumberFieldContext(K), S=[],
    settings=SUnitSettings(smooth_bound=59, generation_radius=2),
    rng=RandomBits(seed=2026),
)
print(result.materialize_generators())
print(result.torsion_order, result.torsion_generator)
print(result.report)
```

An empty `S` requests ordinary units. Use, for example, `S=K.primes_above(3)` for S-units.

The solver enlarges the working set, certifies its complete S-unit group, and then restricts to the requested set using the saturated integer kernel of extra valuations. `result.generators` contains the compact free generators; `torsion_generator` generates the roots of unity. `materialize_generators(max_bits=...)` recovers explicit elements within a certified height budget. `result.to_dict()` exports exact base elements, exponents, valuations, and verification intervals. The demonstration script supports `--output result.json`.

See the [sampling guide](docs/ALGORITHMS_1_2.md) for detailed APIs for Algorithms 1–2.

## Run Algorithms 3–5 separately

The following continues the API example above:

```python
from rigorous_cnt import (
    RelationEngine, FactorBase, prime_ideals_up_to,
    algorithm3, algorithm4, algorithm5,
)

ctx = NumberFieldContext(K)
settings = SUnitSettings(smooth_bound=59, generation_radius=2)
engine = RelationEngine(ctx, settings)
T = FactorBase(ctx, [
    P for P in prime_ideals_up_to(ctx, 59)
    if P not in engine.exceptional_primes
])

relation = algorithm3(engine, K.ideal(1), T, rng=RandomBits(2))
assert relation.verify()  # (alpha) = input_ideal * product(P**v)

observation = algorithm4(engine, T, rng=RandomBits(3))
print(observation.valuations)

# For nontrivial engine.modulus, apply this to each prime factor q:
# exceptional = algorithm5(engine, q, T, rng=RandomBits(4))
```

Algorithm 3 selects its modulus from the residue approximation according to the paper's cases. `modulus_override` supplies an explicit squarefree modulus and records that choice in the report. To exercise the exceptional-unit branch, run:

```bash
PYTHONPATH=src sage -python examples/demo_sunits.py --field imaginary --primes 2,3 --modulus-prime 2
```

## Sampling parameters and output correctness

`SUnitSettings` records the factor-base bound, walk parameters, candidate generation radius, and resource budgets. Defaults are intended for reference experiments, so `sampling_mixing_guarantee` is usually `not_established`. A result must still pass the independent completeness test before it is returned.

After completion, the generated group can be used to check the candidate generation radius and the ray class generation condition. The report records these checks as separate booleans. Passing them does not establish the asymptotic complexity of the entire program.

To request certified sampling parameters, explicitly supply a parameter factory:

```python
from rigorous_cnt import calibrate_walk

settings = SUnitSettings(
    require_mixing_bound=True,
    walk_factory=lambda context, ray, epsilon:
        calibrate_walk(context, ray, epsilon),
)
```

This factory performs the more expensive spectral preprocessing and computes existing class and unit data. Default Algorithm 6 does not use it. The factory supplies sampling parameters; this package still computes the S-unit relations and the evidence for complete generation.

## Implementation map

| Paper component | Main modules |
| --- | --- |
| Algorithms 1–2: lattice and ideal sampling | `sampling.py`, `geometry.py`, `gaussian.py`, `primes.py` |
| Specialized sampling and Algorithms 3–5 | `relations.py`, `factor_base.py` |
| Independent residue interval | `residue.py` |
| Two-pass BKP lattice basis extraction | `postprocess.py` |
| Algorithm 6 and completeness checking | `sunits.py` |
| Compact arithmetic and modular recovery | `compact.py` |
| Exported-record verification | `verification.py` |
| Optional spectral preprocessing | `calibration.py` |

The mathematical notes are in [Algorithms 1–2](docs/MATHEMATICS.md) and [Algorithms 3–6](docs/ALGORITHMS_3_6.md). [Validation records](docs/VALIDATION.md) describe the 48-test suite and the completed examples. The implementation uses maximal orders and absolute number-field representations. End-to-end S-unit validation is concentrated on quadratic fields; sampler tests also include a mixed cubic field and a totally complex quartic field. Nonmaximal orders and large-degree performance are outside the validated scope.

This is an independent implementation. Credit for the underlying algorithms belongs to the original authors. Development and documentation used OpenAI's GPT-6 through Codex; the exact model version is unknown.
