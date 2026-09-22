# rigorous-cnt

A SageMath reference implementation of Algorithms 1–6 in [*Rigorous methods for computational number theory*, arXiv:2512.01588v2](https://arxiv.org/abs/2512.01588v2), by Koen de Boer, Alice Pellet-Mary, and Benjamin Wesolowski.

[中文说明](README.md)

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

An empty `S` requests ordinary units. Use, for example, `S=K.primes_above(3)` for S-units. The free generators are returned separately from the root-of-unity generator. Large outputs can remain compact; explicit materialization has a certified height bound and a resource budget.

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

The mathematical notes are in [Algorithms 1–2](docs/MATHEMATICS.md) and [Algorithms 3–6](docs/ALGORITHMS_3_6.md). [Validation records](docs/VALIDATION.md) describe the 48-test suite and the completed examples. End-to-end S-unit validation is concentrated on quadratic fields; sampler tests also include a mixed cubic field and a totally complex quartic field. Nonmaximal orders and large-degree performance are outside the validated scope.

This is an independent implementation. Credit for the underlying algorithms belongs to the original authors. Development and documentation used OpenAI's GPT-6 through Codex; the exact model version is unknown.
