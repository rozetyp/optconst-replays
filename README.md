# Independent replays of two certificates on the Optimization Constants list

These are replays of two bounds recorded with an asterisk on [teorth/optimizationproblems](https://github.com/teorth/optimizationproblems), both by D. Koshta with Claude (23 Oct 2026). Each replay starts from the certificate bytes the constant page cites. It shares no code with the certificate's own checkers, which were not run.

| constant | claim | certificate | result |
|---|---|---|---|
| [1c](https://teorth.github.io/optimizationproblems/constants/1c.html) Barnard–Steinerberger | $C_{1c} \ge 0.4103333$ | `K24576_0.4103333.json`, commit `a4ee94ef29ff`, SHA-256 `963cb727…c9e` | exact value 134648548624804778248790016/328144294688458468928334649 = 0.41033335274849…, attained at k = 21264; the claimed fraction is matched exactly by two implementations |
| [86a](https://teorth.github.io/optimizationproblems/constants/86a.html) Schur–Siegel–Smyth trace | $C_{86} \le 1.8084990$ | `c86_A19d7_robust_1.8084990.json`, commit `90c57242`, SHA-256 `1a9e108e…45c` | all six conditions certified; f ≥ 10⁻⁷ on the whole support (rigorous), every ∫log\|Q\|dμ > 0 (smallest 2.0000·10⁻⁶), bound = T/m exactly |

## Run

    ./fetch.sh                      # clones both certificate repositories at the cited commits and prints SHA-256

**C_1c.** Python 3 only.

    python3 c1c/bigint_1c.py   certs/c1c/K24576_0.4103333.json      # about 2 s
    python3 c1c/termwise_1c.py certs/c1c/K24576_0.4103333.json      # slower: plain term-by-term sums
    python3 c1c/red_controls_1c.py certs/c1c/K24576_0.4103333.json

**C_86.** Needs [python-flint](https://pypi.org/project/python-flint/) (tested with 0.9.0).

    python3 c86/c86_replay.py certs/c86/certificate/c86_A19d7_robust_1.8084990.json 1e-7 4   # about 25 s on 4 cores
    python3 c86/red_controls_86a.py certs/c86/certificate/c86_A19d7_robust_1.8084990.json red86
    python3 c86/c86_replay.py red86/<file>.json 0 4                                           # each must fail

The outputs are in `logs/`.

## Method

### C_1c

For a step function with equal steps 1/K, the autocorrelation is linear between multiples of 1/K. So the bound is K·min_{0≤k≤K} A_k / S², with A_k = Σᵢ wᵢwᵢ₊ₖ and S = Σᵢ wᵢ, computed in two ways.
- **`bigint_1c.py`** gets every A_k from one big-integer product of the two digit-packed weight sequences (Kronecker substitution). It is written for this replay.
- **`termwise_1c.py`** sums each A_k term by term. It is the verifier from [rozetyp/barnard-steinerberger-lower-bound](https://github.com/rozetyp/barnard-steinerberger-lower-bound), written before this certificate existed. The certificate was obtained by refining that repository's construction.

**Red controls:** a claim 2·10⁻¹⁵ above the exact value is refuted; a negative weight and a non-integer weight are refused.

### C_86

`c86_replay.py` decides the conditions of Smith's criterion [Smi24, Prop. 5.15 and Thm 1.5] in the form used by [SO24, Cor. 1.6]. It reads only `edges`, `densities`, `A`, `b` and `bound`, as described in the certificate's `FORMAT.md`.
- **E1, E3, E6** in exact rationals.
- **E2** (primitive, irreducible, pairwise non-associate) by FLINT factorisation.
- **E4** (∫log|Q|dμ > 0) in Arb ball arithmetic on certified root enclosures. Rational roots are handled exactly, including x = 2, which is also a cell edge.
- **E5** (f = U_μ − Σ b log|Q| ≥ 0 on the support) by adaptive bisection of each support cell. The checker uses the closed form U_μ(x) = Σ_k c_k G(x − e_k), with c_k the density jumps and G(u) = u log|u| − u.
  - On each piece, the two edge terms of the cell are replaced by a linear minorant: the tangent where c·G is convex, the chord where it is concave.
  - The rest is bounded by a second-order Taylor bound with a ball enclosure of |R''| over the piece.
  - Where a polynomial has a root in the piece, a termwise bound is used instead: the exact range of G on each argument interval, and −b log max|Q|.
  - 76,743 pieces; the smallest certified lower bound on a piece is 1.014·10⁻⁷.

**Red controls**, each rejected:
- one density lowered by 1% (bound recomputed): E4 fails;
- the largest weight b raised by 0.1%: E5 fails on 432 cells;
- the bound lowered by 10⁻¹⁵: E6 fails.

## References

- [Smi24] A. Smith, *Algebraic integers with conjugates in a prescribed distribution*, Ann. of Math. 200 (2024); arXiv:2111.12660.
- [SO24] N. Talebizadeh Sardari, B. J. Orloski, *A quantitative converse of Fekete's theorem*, arXiv:2304.10021.

## AI disclosure

The replay code and this text were prepared with the help of an AI assistant (Claude, Anthropic). The checks were run and their outputs recorded in `logs/`.
