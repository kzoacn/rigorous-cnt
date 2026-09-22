# Algorithms 1–2: mathematical specification and limitations

Reference version: [arXiv:2512.01588v2](https://arxiv.org/pdf/2512.01588v2), revised February 19, 2026. These notes record the implementation arguments and checks; they are not a formal proof.

## Exact objects and coordinates

Number-field elements belong to Sage's exact number fields. Embeddings are retained as exact algebraic numbers through `QQbar` and `AA`. Infinite places comprise all real embeddings and the embedding with positive imaginary part from each complex conjugate pair. Complex scaling factors are stored as pairs of rationals, and their physical coordinates are real and imaginary parts.

The paper's Minkowski squared norm weights complex places by 2. Short-basis checks use this weighted norm. Geometric sampling takes place in unweighted real coordinates, where each complex coordinate forms a disk. The unweighted length is at most the paper's weighted length, so the certified short basis satisfies the geometric condition used by the sampler.

A supplied integral basis is checked by comparing its **additive Z-span** with the maximal order. Checking only that its elements generate the unit ideal is insufficient. Automatic integral-basis computation is preprocessing: the paper takes a reduced integral basis as input, and the implementation does not include this preprocessing cost in a complexity claim.

## Local congruences and fractional CRT

For a finite modulus m, the nonzero element tau must be a unit at every p dividing m. The ray conditions are `ord_p(alpha-tau) >= ord_p(m)` and `alpha/tau > 0` at the selected real places. With the unit ideal as modulus, only nonvanishing and the selected sign conditions remain.

For a fractional ideal b coprime to m and a shift gamma integral at m, construct the fractional ideal `J = O_K + b + (gamma) + (tau)`. It has valuation zero at every p dividing m and accommodates the denominators needed elsewhere.

Integer HNF solves `tau-gamma = u+v`, with u in b and v in mJ. Set `gamma_m = gamma+u`; it satisfies the local congruences. Since `b ∩ mJ = bm`, all elements satisfying them form `gamma_m+bm`. Finally, subtract a lattice vector by rounding along an exact basis of bm to obtain a short shift.

If gamma has negative valuation at a prime dividing m, every element of gamma+b still has negative valuation there and cannot satisfy the unit congruence prescribed by tau. The program reports empty support.

## Short bases and geometric guarantees

The radius uses the algebraic expression from Algorithm 1:

```text
r = 48 * omega * block_size^(2*n/block_size) * n^(7/2)
    * abs(discriminant)^(3/(2*n)) * abs(N(x*b*m))^(1/n).
```

First, LLL on an integer matrix obtained from a rational approximation supplies an integer transformation. Apply that transformation to the original exact field elements and accept it only if ball arithmetic certifies that every true basis-vector length is at most `r/(24n²)`. If this fails, use LLL on an exact algebraic Gram matrix and, when needed, the block HKZ schedule of Algorithm 7. Every path retains the unimodular transformation and checks the geometric bound. The output distribution depends only on the target set, so another basis satisfying the same condition preserves uniformity. A failure of fpLLL at the `eta=1/2` boundary triggers a fallback to exact NTL LLL.

The SVP backend is a checkable enumeration routine for small instances. It is not the Kannan subroutine used to obtain the paper's `poly(input)*b^b` bound, and termination is controlled by the actual geometric bound. Consequently, the full bit-complexity conclusion of Lemma 8.5 has not been established. If the backend cannot certify the geometry within its budget, it raises an error instead of sampling without the required bound.

## Uniformity in Proposition 8.8

Use a column basis B for the exact lattice. The target is S+t, where S is a product of real intervals and complex disks. A real interval with a sign restriction has radius r/2 and center ±r/2; the sign of its center accounts for both the image of tau and the real scaling factor. A final strict ray check excludes boundary cases such as zero.

Let `D = 2 sum ||b_i||` and `epsilon_geom = 1/(6n)`. Here D is a geometric quantity, not the field discriminant. The code verifies that `D/epsilon_geom` is at most the inradius of S and chooses U to strictly bound the lengths of points in S and the shift t. The Frobenius norm in real coordinates bounds the inverse matrix.

Choose a dyadic grid `1/N`, round entries to obtain C and t_tilde, and verify:

```text
n/N <= D
||C-B||_F * ||B^(-1)||_F <= D/U
sum ||c_i|| <= D
||t-t_tilde||_infinity <= 1/(2N)
||t_tilde|| < 2U
det(C) != 0
```

The sampler then performs the three steps of Proposition 8.8:

1. Set `w0 = round(C^(-1)*t_tilde)`.
2. Sample u uniformly from `(1+4 epsilon_geom)S ∩ (1/N)Z^n`, set `v = round(C^(-1)u)`, and reject if Cv is outside `(1+3 epsilon_geom)S`.
3. Return only if the exact point `B(v+w0)` lies in S+t and satisfies the strict ray conditions.

`round(q) = floor(q+1/2)` commutes with integer translation, including half-integer ties. Each half-open fundamental cell of C contains the same number of grid points because the columns of C belong to `(1/N)Z^n`. The first acceptance stage is therefore uniform. The geometric conditions cover every target lattice point, and final rejection using exact membership preserves uniformity on the target set. Small tests enumerate all proposal grid points and compare the mass assigned to each target point.

Complex disks are sampled separately by two-dimensional rejection and then combined as a product. This avoids the deteriorating acceptance rate of joint rejection across all complex coordinates.

## Discrete Gaussians and finite precision

The implementation of Lemma 2.22 uses rational Gram–Schmidt data. In each direction, the center and squared width are rational, and the one-dimensional support is:

```text
Z ∩ [c_i - s_i*sqrt(log(2*r^2/epsilon_G)),
     c_i + s_i*sqrt(log(2*r^2/epsilon_G))]
```

Here r is the lattice rank; rectangular bases use the same metric on their span. The program checks the width hypothesis of Lemma 2.22. It proposes uniformly on this finite interval and accepts with probability `exp(-pi*(z-c_i)^2/s_i²)`. This follows the lemma's truncation and GPV error allocation and retains its output-length bound.

Acceptance compares an Arb ball enclosure of the probability with a lazy uniform binary real. If the intervals overlap, the program reads more random bits and increases precision. With unlimited resources, these decisions introduce no additional bias from ordinary floating-point rounding. Resource interruptions can depend on the random path, so an experiment that retains only successful runs is not automatically evidence for the original distribution.

Sage may interpret `QQ(real_float)` through rational reconstruction, which can destroy enclosure. Ball endpoints and midpoints therefore use `exact_rational()`, rather than default rational reconstruction. A dedicated test checks the direction of square-root interval endpoints.

## Discretization and distortion in Algorithm 2

Set `s = 1/n²` and `epsilon_G = epsilon/4`. The exponent `4*n²*s+1` in Algorithm 2's grid scale is 5. Certified ball arithmetic computes the paper's upper bound, and delta is rounded down to a dyadic rational. The error argument in Section 9 uses an upper bound on delta, so rounding down is permitted.

Sampling on the row basis `(delta/n)*(e_i-e_(i+1))` gives a rational vector a whose coordinates sum to zero. Ball arithmetic constructs positive rationals A_i approximating `exp(a_i/w_i)`, where w_i is 1 or 2.

To impose exact multiplicative norm normalization, reserve one real place, if any, and determine its factor from the product of the others. If all places are complex, determine one complex factor from the product of the other positive rational factors. After construction, certify every relative error as at most `delta/(2n)` and verify the exact equality `product(A_i^w_i) = 1`. Independently rounding every factor would not preserve this equality.

Pass the scale A*y and the ideal bbar obtained from the random walk to Algorithm 1. Its exact field element alpha is also the beta obtained by cancelling A*y at the end of Algorithm 2. The intermediate scaled object remains available in `BoxSample`.

## Uniform prime-ideal sampling

The implementation of Lemma 5.4 samples p uniformly from the integers `{1,...,B}` and retries if p is not prime. If k prime ideals above p satisfy the norm, modulus, and subgroup conditions, choose one uniformly and accept it with probability k/n.

Each eligible prime ideal has proposal mass `(1/B)*(1/k)*(k/n) = 1/(nB)`, so the accepted distribution is uniform. `K.primes_above(p)` includes inert and ramified primes and primes dividing the index of the defining power order, rather than only residue-degree-one primes. Exhaustive probability-mass tests include these cases.

## Limits of the mixing guarantees

`WalkParameters(B,N)` fully defines a random walk but does not establish Corollary 6.5 from two integers alone. A sufficient B also depends on hidden constants in the analytic estimates of Appendix A.1 and on subgroup information.

The step count can be made explicit separately. Lemma 5.1 gives `Vol(Pic^0_m) <= |Delta| * N(m)`. Substituting this in Corollary 6.5 and assigning an L1 budget of epsilon/2 to the walk yields:

```text
N >= ceil(7*n + 2*log(2/epsilon)
          + log(|Delta|*N(m)) - log(subgroup_index) + 2).
```

`WalkParameters.from_prime_bound` computes this conservative N using directed ball endpoints and does not call a class-group or unit-group solver. It still requires a sufficiently large B; selecting N alone is not a mixing proof.

Three states are distinguished: a bound certified by direct spectral calibration for a supported small field; an externally justified bound supplied by the caller; and an explicitly permitted run without an established mixing bound. External justification remains a mathematical premise. Descriptive text is not a verified proof, and the program does not automatically establish the algebraic properties of a general subgroup oracle.

Proposition 4.6 of [Koen de Boer's doctoral thesis](https://ir.cwi.nl/pub/32136/32136D.pdf) also states its Hecke eigenvalue bound using big-O notation. An explicit prime-count threshold does not replace the missing eigenvalue constant.

Membership, geometry, and rational distortion can therefore be checked for each run, while closeness to stationarity and the density-based success bound of Theorem 9.5 remain `not_established` without calibration or an external mixing premise. Full bit complexity remains uncertified in every run. Passing tests, observed sampling frequencies, and descriptive configuration text do not change these statuses.

## Direct spectral calibration for small fields

`calibrate_walk` uses the Fourier decomposition and tail bound of Appendix A.1, replacing the analytic estimate with hidden constants by directly computed eigenvalue bounds. This is additional preprocessing supplied by this project. It calls PARI class-group, unit-group, and ray class-group facilities and requires the underlying bnf data to be certified with `proof=True`. Its cost is not covered by a claim to the paper's complexity bound. The calibrator currently handles the full group `G = Pic^0_m`; the core sampler continues to support general subgroups.

For cyclic ray class generators g_j with orders d_j, exact ideal arithmetic verifies:

```text
g_j^d_j = (alpha_j)
P = (alpha_P) * product(g_j^e_j)
```

Both alpha_j and alpha_P must be principal ray elements, with congruence and positivity checked. PARI and Sage ideal coordinates are converted through explicit field elements instead of assuming identical internal integral bases. Conditions on selected real places are translated by matching the bnf ordering of real roots with exact real embeddings. Failure to distinguish the embeddings within the budget is reported as a resource failure.

An integer exponent matrix represents the map from the full unit group, including roots of unity, to `(O_K/m)^*` and the real sign groups. Take its integer kernel, project to the free-unit coordinates, and apply HNF to obtain a basis of the complete ray-unit lattice `Lambda = Log(O_K^* ∩ K^{m,1})`. Retaining torsion coordinates in the kernel prevents a congruence that roots of unity can cancel from becoming an unnecessary constraint on free units.

Let L be a row basis of Lambda and r its rank. For every k in Z^r,

```text
u = L^T * (L*L^T)^(-1) * k
```

lies in the dual lattice. Let Log_0 be the logarithmic embedding with the coordinate mean subtracted. Every character is described by u and finite parameters `a_j in {0,...,d_j-1}`. Its phase at a prime ideal P is:

```text
<u, Log_0(alpha_P)>
  + sum_j e_j * (<u, Log_0(alpha_j)> + a_j) / d_j.
```

Averaging `exp(2*pi*i*phase)` over all prime ideals of norm at most B coprime to m gives the corresponding Hecke eigenvalue. Ball arithmetic encloses every phase, trigonometric evaluation, and absolute value. Exclude the constant character with zero frequency and all finite parameters zero.

To cover all characters with `||u|| <= R`, use `k_i = <L_i,u>` and hence `|k_i| <= ||L_i|| R` to construct a finite integer box with directed upper bounds. Discard a frequency only when ball arithmetic proves its norm exceeds R; retain uncertain boundary cases. Explicit budgets limit group size, box size, and the total character count.

Set `s = 1/n²` and `s0 = 1/(2800n²)`. The smoothing-parameter bound in the proof of Corollary 6.5 allows s0 as a lower bound for s'. Let V be the upper bound on the volume of Pic^0_m from Lemma 5.1. Appendix A.1 gives:

```text
L1_error^2 <= 2 * V * s0^(-r) * (c^(2*N) + exp(-2*R^2*s^2)),
```

where c bounds the eigenvalues of all nonconstant characters with `||u|| <= R`, and R satisfies `sqrt(2) R s >= sqrt(r)`. With `e = epsilon/2`, first choose rational R so the tail contribution is at most e²/2. Increase B until c is certified to be at most 3/4, then choose N and directly check that the entire right-hand side is at most e². The Gaussian tail bounds the frequencies outside the finite enumeration; apparent mixing after frequency truncation is not used as a proof. When r is zero there is no high-frequency tail. A trivial group has zero mixing error.

The returned `SpectralMixingCertificate` records the spectral bound, cutoff radius, volume bound, character count, precision, B,N, and L1 bound, and is bound to the original field context and modulus. It checks consistency between the run and its calibration; it is an evidence record, not an independent formal proof file.

For an imaginary quadratic field of class number two, tests compare with the exact error of a two-state random walk. For a real quadratic field, the known unit length checks frequency coverage and the total error inequality. Further tests cover a mixed cubic field, a totally complex quartic field, and a nontrivial modulus with sign conditions at only some real places.

## Implementation conventions

- Geometric maps use column bases; lattice reduction uses row transformations. Exact matrices are retained at conversions.
- Proposition 8.8's N is chosen by the actual matrix inequalities, avoiding reliance on an easily mistranscribed inverse-norm inequality.
- Negative ray signs, fractional local congruences, and strict nonvanishing are handled explicitly.
- Exact algebraic comparisons replace a separate implementation of root isolation. This does not assert that every Sage backend operation automatically has the bit complexity required by the paper.
