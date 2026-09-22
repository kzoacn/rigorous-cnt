# rigorous-cnt: A SageMath Reference Implementation of Rigorous Ideal Sampling and S-Unit Computation

Kaiyi Zhang

## Abstract

We present rigorous-cnt, a SageMath reference implementation of Algorithms 1–6 from de Boer, Pellet-Mary, and Wesolowski's *Rigorous methods for computational number theory*. The software connects ideal-lattice sampling, smooth relations, S-unit observations, exceptional units, and lattice postprocessing. It separates exact algebraic checks, certified numerical bounds, conditional completeness, and unverified sampling-time assumptions. Complete generation is checked using an independent Belabas–Friedman interval for the Dedekind zeta residue and an integer lattice-index test, conditional on the Generalized Riemann Hypothesis. Compact power products avoid premature expansion, while a verifier reconstructs the relevant evidence from exported records. Validation comprises 48 passing tests and small quadratic-field computations, including an ordinary unit group and an S-unit group constrained by a nontrivial ideal class. Tests reject full-rank proper subgroups and exercise the exceptional-unit branch. The default solving path does not invoke an existing class-group, unit-group, or S-unit-group solver. The implementation does not establish the source paper's full bit-complexity bounds, and its default sampling parameters have no general certified mixing-time guarantee. The contribution is an executable, inspectable implementation and validation artifact, rather than a new number-theoretic algorithm or a performance claim against SageMath.

## 1. Scope and provenance

Turning a theoretical number-field algorithm into executable software requires explicit choices about numerical error, lattice coordinates, stopping conditions, and output representation. An implementation can preserve an algebraic identity while failing to establish a distributional premise or an asymptotic running-time bound. These distinctions are especially relevant when a randomized search is followed by an independent completeness check.

This work implements Algorithms 1–6 of the February 19, 2026 version of [de Boer, Pellet-Mary, and Wesolowski](https://arxiv.org/abs/2512.01588v2). Algorithm numbering throughout refers to that version. The original authors supply the underlying mathematical algorithms and analysis. Our contribution is the implementation, explicit handling of several representation and certification details, and reproducible validation. We do not implement every algorithm in their 154-page paper or claim its full class-group algorithm and complexity theorem as a software result.

The public artifact is [rigorous-cnt](https://github.com/kzoacn/rigorous-cnt), version 0.2.0. The code and validation data described here are pinned to commit [`d75af574dd9de4dcd7dbfef5f4ff87dec5264a13`](https://github.com/kzoacn/rigorous-cnt/tree/d75af574dd9de4dcd7dbfef5f4ff87dec5264a13). SageMath provides exact number-field arithmetic, ideal operations, integer lattices, and ball arithmetic. Using these facilities does not make the implementation an independent arithmetic library. The narrower independence claim concerns the default solving path: it does not obtain its answer from a pre-existing class-group, unit-group, or S-unit-group solver.

## 2. Mathematical interface and evidence

Let $K$ be an absolute number field of degree $n$, signature $(r_1,r_2)$, and absolute discriminant $D$. Set $s=r_1+r_2$. For a finite set $S$ of finite prime ideals, the target group is

$$
\mathcal O_{K,S}^{\times}=\{\alpha\in K^{\times}:v_{\mathfrak p}(\alpha)=0\text{ for every }\mathfrak p\notin S\}.
$$

Its free rank is $|S|+s-1$. The software returns compact generators for this free part and a separate generator for the roots of unity. It works in the maximal order.

The implementation reports different kinds of evidence separately.

| Property | Evidence and scope |
| --- | --- |
| Ideal identities and finite support | Exact ideal arithmetic and integer valuations |
| Lattice transformations | Exact integer matrices and reconstruction checks |
| Required numerical inequalities | Directed real-ball bounds with precision escalation |
| Complete generation at successful termination | Independent residue interval, class-generation witnesses, and index-one test; conditional on GRH |
| Default sampling mixing bound | Not established in general |
| Optional spectral mixing calibration | Certified for supported instances, using explicit class/unit preprocessing |
| Full asymptotic bit complexity | Not established for this implementation |

Here GRH refers, for the residue certificate, to the hypotheses on the Dedekind zeta function of $K$ and the Riemann zeta function in [Belabas–Friedman, Theorem 1](https://arxiv.org/abs/1305.0035). The computations are not formal proofs in a proof assistant: they also rely on the correctness of the implementation and its arithmetic dependencies. Resource exhaustion raises an exception rather than returning an uncertified partial generating set as complete.

## 3. From ideal sampling to S-unit observations

### 3.1 Algorithms 1 and 2

Algorithm 1 samples in a shifted ideal lattice intersected with the required bounded Minkowski region and ray conditions. The implementation uses the rational-grid proposal and rejection construction underlying Proposition 8.8. Exact algebraic arithmetic decides membership, including complex-disk boundaries and real signs. Coordinate rounding uses $\lfloor x+1/2\rfloor$, preserving equivariance under integer translation at half-integer ties.

An approximate integer LLL computation can propose a useful basis transformation. The transformation is applied to exact field elements and accepted only after the required geometric bound is certified. Exact lattice routines provide a fallback. This engineering choice does not establish the theoretical cost of the paper's short-basis subroutine.

Algorithm 2 combines bounded prime-ideal sampling, a random walk, a rational implementation of the Gaussian construction in Lemma 2.22, and norm-one distortion before invoking Algorithm 1. Prime selection includes the acceptance adjustment needed to compensate for different splitting patterns. Tests enumerate proposal masses on small cases; reproducible random seeds also check replay and exact ideal identities. These tests complement the mathematical construction but do not replace its proof.

The generic efficient prime bound contains constants that are not fully instantiated here. Default reference parameters therefore carry no general mixing certificate. The optional `calibrate_walk` routine instead bounds low-frequency characters and a Gaussian tail for supported small instances. It uses existing class and unit data during preprocessing. That dependency is explicit and is not part of the default S-unit solving path.

### 3.2 Specialized sampling and Algorithms 3–5

The specialized sampler for Theorem 13.2 draws a uniform element of $(\mathcal O_K/\mathfrak m)^\times$ using prime-power quotient lattices, rejection of nonunits, and the Chinese remainder theorem. It selects a dyadic error parameter satisfying the prescribed upper bound. Exact ideal-content and common-scale normalizations reduce coefficient growth while transforming the ideal, residue condition, and sampling region together.

Algorithm 3 repeatedly samples until it finds an exact smooth relation

$$
(\alpha)=\mathfrak a\prod_{\mathfrak p\in T}\mathfrak p^{e_{\mathfrak p}},\qquad e_{\mathfrak p}\geq 0,
$$

for the working factor base $T$. The input ideal may be fractional. Successful candidates are checked by exact reconstruction. The report records applicable prime-count conditions separately from the mixing status, so satisfying a combinatorial condition alone is not presented as a proved success probability.

Algorithm 4 samples a discretized divisor, with finite integer coordinates and infinite coordinates on the paper's rational grid. It uses certified rational approximations to the required exponentials and calls Algorithm 3 to produce an S-unit observation. The finite valuations and weighted logarithmic coordinates are retained together. The supplied `generation_radius` is a candidate bound; a later check can validate it, but an unsuccessful check remains visible in the report.

Algorithm 5 supplies the exceptional elements associated with primes dividing the modulus. It samples using the modulus with the target prime removed and requires valuation exactly one at that prime. A modulus override permits direct testing of this branch and is recorded as an override, rather than being attributed to the paper's automatic modulus-selection rule.

## 4. Basis extraction and complete generation

### 4.1 Numerical lattice postprocessing

The implementation follows the two-pass BKP construction in Lemma 23.4 and Theorem 23.5 of the source paper. It constructs augmented integer matrices from certified approximations and retains the unimodular transformations produced by LLL. A lower bound on nonzero logarithmic-lattice lengths permits certification of zero relations. A positive interval Gram determinant certifies independence of the retained rows.

The second pass requests sufficiently precise source data again. Compact products may involve severe cancellation, so working precision accounts for exponent sizes and the required final interval width. Approximate real columns are never treated as exact integer data for Hermite normal form. A numerical failure at the selected fpLLL reduction boundary triggers an exact NTL LLL fallback.

### 4.2 Independent target covolume

The stopping test needs a bound for $h_KR_K$, where $h_K$ is the class number and $R_K$ the regulator. The implementation evaluates the explicit finite prime-power formula and analytic error bound of Belabas–Friedman with ball arithmetic. The analytic class number formula then gives an interval for $h_KR_K$. This calculation uses prime splitting, discriminant and signature data, and roots of unity; it does not compute fundamental units or a class group. Agreement between successive numerical truncations is not used as an error certificate.

The working factor base must generate the ordinary ideal class group. The implementation establishes this by accounting for every prime ideal up to the Minkowski bound, either directly in the base or by an exact principal relation. This explicit enumeration is useful for the validated examples but can be expensive as the degree and discriminant grow.

A coordinate convention matters when comparing volumes. Write

$$
\operatorname{Log}(\alpha)=(w_\nu\log|\sigma_\nu(\alpha)|)_\nu,
$$

where $w_\nu=1$ at real places and $w_\nu=2$ at complex places. The implemented rows are

$$
\ell_T(\alpha)=\bigl((-v_{\mathfrak p}(\alpha))_{\mathfrak p\in T},\operatorname{Log}(\alpha)\bigr).
$$

They lie in a hyperplane whose normal vector is $((\log N\mathfrak p)_{\mathfrak p\in T},1,\ldots,1)$. When $T$ generates the class group, the full lattice's covolume in this ambient Euclidean metric is

$$
V_T=h_KR_K\sqrt{s+\sum_{\mathfrak p\in T}(\log N\mathfrak p)^2}.
$$

The slope factor converts the normalized volume used in the paper's discussion to the Euclidean Gram volume actually computed by the software. For example, in $K=\mathbb Q(i)$ with $T=\{(3)\}$, the generator $3$ gives the row $(-1,\log 9)$, of length $\sqrt{1+(\log 9)^2}$. A regression test checks this normalization.

For a full-rank candidate lattice $L'\subseteq L_T$, the ratio $\operatorname{covol}(L')/V_T$ is the positive integer $[L_T:L']$. Interval arithmetic bounds that ratio. Complete generation is accepted only when its possible positive integers form the singleton $\{1\}$. Full rank alone does not suffice.

### 4.3 Algorithm 6 and arbitrary requested prime sets

Algorithm 6 collects Algorithm 4 observations until basis extraction and the index test establish the complete working group. It adds the exceptional elements and then restricts to the user's requested $S$. The restriction is the saturated integer kernel of valuations at working primes outside $S$; its transformation acts on compact exponent vectors. Thus the restriction is exact even when logarithmic coordinates require increasing precision.

Ordinary class generation and ray class generation are distinct conditions. After completion, the implementation can additionally test the images of working units in the finite residue-unit group. One regression example uses $\mathbb Q(i)$ modulo $(3)$: roots of unity alone do not account for the residue-unit group, while adjoining $1+i$ does. The reports retain this distinction.

## 5. Compact outputs and record verification

A compact element is a power product $\beta=\prod_j\alpha_j^{e_j}$ with exact field elements and integer exponents. Expanding its factors independently can create enormous intermediates even when the final result is small.

For explicit recovery, exact valuations first provide a denominator $d$ making $d\beta$ integral. Certified logarithmic bounds and the inverse integral-basis embedding matrix bound its integral-basis coefficients. Modular evaluation at a sufficiently large prime-power modulus then recovers the unique coefficients within those bounds. Negative exponents use exact inverse field factors. Output-size limits preserve the compact representation when recovery would be too large. Tests include cancellation with exponents of size $10^{50}$, fractional results, roots-of-unity phases, and a non-power integral basis.

The exported JSON records contain exact source elements, exponents, ideals, valuations, and diagnostic information. The `verify_sunit_record` entry point recomputes the residue interval, the complete working lattice, the restriction kernel, and the final index. Stored success flags and stored volume claims are deliberately untrusted. This is a second checking path within the same software and arithmetic stack, not an independently implemented formal verifier.

## 6. Validation

The publication check on September 22, 2026 passed all 48 tests under SageMath 10.7, Python 3.12.13, and fpylll 0.6.4. The repository includes the complete test log, environment record, and SHA-256 hashes of the saved outputs. Test coverage includes exact small-case sampling masses, boundary rounding, fractional ideals, mixed and complex signatures, Gaussian bounds, residue intervals, unknown lattice rank, proper subgroups, compact recovery, and exported-record verification.

End-to-end solving tests disable `pari_bnf` during the solve. Sage's S-unit interface is called afterwards for comparison. This checks the intended solver boundary while retaining Sage as the underlying arithmetic system.

| Field and requested group | Settings | Observed result |
| --- | --- | --- |
| $\mathbb Q(\sqrt 2)$, $S=\varnothing$ | Smoothness bound 59, radius candidate 2, seed 2026 | 20 relations, 5,452 sampling calls; free generator $\sqrt 2-1$; torsion order 2 |
| $\mathbb Q(\sqrt{-5})$, $S$ the prime above 2 | Smoothness bound 97, radius candidate 2, seed 11 | 25 relations, 5,244 sampling calls; free generator $-2$; torsion order 2 |
| $\mathbb Q(i)$, $S$ the prime above 3 | Smoothness bound 59, radius candidate 2, seed 7 | Automated comparison gives exponent $\pm1$ in Sage's free basis; torsion order 4 |
| $\mathbb Q(i)$, $S$ the primes above 2 and 3 | Smoothness bound 59, radius candidate 2, seed 31; modulus forced to the prime above 2 | Exceptional-unit branch completed; generators $-3i$ and $-2187(i+1)$; determinant $-1$ in Sage's free coordinates |

The first two rows have saved exact records. They are individual reproducible runs, not averages or a performance benchmark. Both used uncalibrated default mixing parameters. In the real quadratic example, attaining full rank still left an index near $10^{13}$; later reported integer intervals included $[52,56]$ and $[6,6]$ before reaching $[1,1]$. Continuing beyond full rank was therefore essential.

In the second example, the prime above $2$ is nonprincipal and its square is principal. The returned generator has valuation $2$ at that prime, as the class-group obstruction requires. Its generation-radius check remains false in the saved report, although the separate GRH-conditional completeness test succeeds. This is evidence about output correctness, not verification of the corresponding sampling-time premise.

A corruption test squares the final generator in the real quadratic record while preserving its stored success flag. The modified output remains full rank, but the verifier obtains index $2$ and rejects it. A separate calibrated Algorithm 3 experiment over $\mathbb Q(\sqrt 2)$ used walk parameters $B=38$, $N=114$, and error budget $2^{-41}$; it returned an exact relation after 22 attempts with a verified spectral bound. That experiment used the optional class/unit preprocessing and should not be conflated with the default solving path.

## 7. Reproduction and limitations

With SageMath available, the pinned artifact can be checked using:

```bash
git clone https://github.com/kzoacn/rigorous-cnt.git
cd rigorous-cnt
git checkout d75af574dd9de4dcd7dbfef5f4ff87dec5264a13
PYTHONPATH=src sage -python -m unittest discover -s tests -v
PYTHONPATH=src sage -python examples/verify_sunits.py examples/results/real_units.json
PYTHONPATH=src sage -python examples/verify_sunits.py examples/results/class2_sunits.json
```

The two saved searches can be regenerated with:

```bash
PYTHONPATH=src sage -python examples/demo_sunits.py \
  --field real --smooth-bound 59 --generation-radius 2 --seed 2026
PYTHONPATH=src sage -python examples/demo_sunits.py \
  --field class2 --primes 2 --smooth-bound 97 --generation-radius 2 --seed 11
```

The exceptional branch can be exercised with `--field imaginary --primes 2,3 --modulus-prime 2 --seed 31`. Randomized searches are materially more expensive than rechecking the saved records. Fixed seeds support reproducibility within the tested software environment; identical behavior across different dependency versions is not promised.

The principal limitations are explicit. Completeness relies on GRH through the residue bound. Default reference parameters lack a general mixing certificate. A candidate generation-radius premise may remain unverified. Exact fallback lattice routines, explicit small-prime enumeration, and optional spectral preprocessing do not establish the source paper's asymptotic cost bounds. End-to-end S-unit validation is concentrated on quadratic fields; higher-degree sampler tests do not constitute higher-degree S-unit benchmarks. There is no measured speed advantage over SageMath and no claim of formal verification.

Useful next steps are to instantiate further explicit parameter bounds, expand validation across degrees and signatures, and benchmark complete computations with preprocessing and certification costs reported separately. The present artifact provides the modular interfaces and saved evidence needed to investigate those questions.

## AI assistance and availability

Development, test construction, and manuscript preparation used OpenAI's GPT-6 through the Codex client. The exact model version is unknown. This disclosure does not authenticate a model or imply independent human or peer review. Authorship is Kaiyi Zhang; no institutional affiliation is asserted.

The manuscript is distributed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/). This manuscript license does not assign a license to the linked software or its dependencies. Code, examples, and validation records are available in the [public repository](https://github.com/kzoacn/rigorous-cnt).

## References

1. Koen de Boer, Alice Pellet-Mary, and Benjamin Wesolowski. *Rigorous methods for computational number theory*. arXiv:2512.01588v2, February 19, 2026. [Versioned preprint](https://arxiv.org/abs/2512.01588v2).
2. Karim Belabas and Eduardo Friedman. *Computing the residue of the Dedekind zeta function*. arXiv:1305.0035v1, 2013. [Preprint](https://arxiv.org/abs/1305.0035v1).
3. The Sage Developers. *SageMath*, version 10.7, used for the reported computations. [Project](https://www.sagemath.org/); [number-field reference documentation](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html).
4. The PARI Group. *PARI/GP: General number fields*. [Official reference documentation](https://pari.math.u-bordeaux.fr/dochtml/html-stable/General_number_fields.html).
