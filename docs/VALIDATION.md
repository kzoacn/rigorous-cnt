# Validation record: Algorithms 1–6

Environment: SageMath 10.7 and Python 3.12. Date: September 22, 2026. The [publication environment record](validation-environment.json) includes exact version numbers and fixture hashes; the [test log](validation-tests.txt) records the publication check.

All 48 tests passed. Coverage includes sampling probability masses, lattice transformations, Gaussian support bounds, local congruences, residue error intervals, full-rank proper subgroups, exceptional units, compact representations, and record verification. End-to-end solving tests disable `pari_bnf` during the solve and call Sage's S-unit interface for comparison only afterwards.

A ray class test also covers Q(i) modulo (3): roots of unity alone do not span the residue-unit group, but adjoining 1+i makes the image surjective. This distinguishes ordinary class generation from ray class generation. A version 0.2.0 wheel was built and passed an import check.

## Completed end-to-end examples

| Field and target | Settings | Observed result |
| --- | --- | --- |
| Q(sqrt(2)), ordinary unit group | smooth_bound=59, generation_radius=2, seed=2026 | 20 relations and 5,452 sampling calls; free generator sqrt(2)-1; torsion generator -1 |
| Q(sqrt(-5)), S the prime above 2 | smooth_bound=97, generation_radius=2, seed=11 | 25 relations and 5,244 sampling calls; free generator -2; torsion generator -1 |
| Q(i), S the prime above 3 | smooth_bound=59, generation_radius=2, seed=7 | Automated test checks that the free generator has exponent ±1 in Sage's basis |
| Q(i), S the primes above 2 and 3 | smooth_bound=59, generation_radius=2, seed=31; modulus explicitly set to the prime above 2 | Exceptional-unit branch completed; generators -3i and -2187(i+1); determinant -1 in Sage's free coordinates |

Exact outputs for the first two examples are saved as [real_units.json](../examples/results/real_units.json) and [class2_sunits.json](../examples/results/class2_sunits.json). Their completeness conclusions are conditional on GRH. Their default reference walk parameters have no certified mixing-time guarantee. The class-number-two example's candidate generation-radius check remains false in its saved report.

In the real quadratic example, after full rank was attained, the possible indices still included a value near 10^13, then 52–56, and then 6. Further relations reduced the interval to `{1}`. This checks the distinction between full rank and complete generation. In the class-number-two example, the generator -2 has valuation 2 at the nonprincipal prime ideal, as required by the class-group constraint.

## Recheck saved outputs

```bash
PYTHONPATH=src sage -python examples/verify_sunits.py examples/results/real_units.json
PYTHONPATH=src sage -python examples/verify_sunits.py examples/results/class2_sunits.json
```

The verifier recomputes the residue interval, the basis of the working lattice, the integer kernel of extra valuations, and the final lattice index. It ignores completion flags and volume claims stored in the JSON. A test squares the final generator in a record while retaining full rank and the original success flag; the verifier must obtain index 2 and reject the record.

## Local integration check with a certified mixing bound

An Algorithm 3 run over Q(sqrt(2)) used explicit spectral calibration with B=38 and N=114. Its specialized sampler used epsilon=1/2199023255552. It found an exact relation after 22 candidates, reported `verified_spectral_bound`, and satisfied the combinatorial success-probability conditions. This calibration uses additional class/unit preprocessing and is not a dependency of default Algorithm 6.

## Scope not established

The general asymptotic bit-complexity bounds have not been certified. Candidate generation-radius bounds remain unverified in some examples. End-to-end S-unit validation is concentrated on quadratic fields; sampler tests also cover a mixed cubic field and a totally complex quartic field. These individual relation counts and observed running times do not establish performance guarantees for higher degrees or larger discriminants.
