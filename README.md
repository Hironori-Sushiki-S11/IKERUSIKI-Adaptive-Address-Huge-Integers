# IKERUSIKI Adaptive Address — Finite Identification of Huge Integers

English | [日本語](README.ja.md)

Experiment date: 2026-09-09. This is a new application extending the existing work in its publication and external validation phase.
This repository publishes the code, execution results, validation records, and reproduction instructions for finite identification of huge integers.

## Position within Adaptive Address

Adaptive Address is the principal research object; the Mersenne and ultra-huge-integer corpora here are application cases, not external benchmarks that confer validity on AA. The operative chain is **Current Held-Information State → Collision → Structural Saved Information → Blind Spot → Pick Up → Select → Address Update → Complete Identification → EXPAND when required**. The reported success is complete distinction *within the stated finite corpora*, with the observed axes and costs.

This documentation keeps independent claims separate: collision-free distinction, the computational effect of modular exponentiation, adaptive selection cost, formal convergence conditions, primality/compositeness, and cryptographic consequences. Neither RSA factorization nor an unrelated primality test is an upstream certificate for the observed AA identification; nor does corpus distinction alone establish factoring or primality. The evidence tables below retain those boundaries.

The upstream repository's historical RSA-oriented name records an earlier application, not the current AA evaluation hierarchy.

## Results

The experiment completely distinguished **52 known Mersenne primes** using axes **467 → 11**, without expanding the integer bodies.
Collision pairs decreased **1326 → 4 → 0**. The axes were selected from candidates 2–512; EXPAND was not needed.

It also completely distinguished 52 integers `2**p-1`, defined by `p = 10**100 + i (i=0..51)`, using axis **53**.
Each integer has exactly `p` bits, with a maximum bit length of `10**100+51`.
These exponents are explicitly chosen consecutive integers, not a random sample. No primality claim is made for this set.
Within this generating formula, the exponent itself is also an identifier. This experiment verifies that Adaptive Addresses can be formed by observing exact residues of the represented integers, rather than merely comparing exponent labels.

Extending to the combined set of 104 integers while retaining the two axes from the known-prime set left two collision pairs.
Adding axis **5** completely distinguished the combined set with **[467,11,5]**.
Re-evaluating the same 52 known primes with their existing Address left zero collisions and required no additional axis.

| Experiment | Count | Final axes | Collision pairs from the starting state | Median of 5 runs (seconds) |
|---|---:|---|---|---:|
| Known Mersenne primes | 52 | [467, 11] | 1326 → 4 → 0 | 0.029489 |
| Exponents at the 10^100 scale | 52 | [53] | 1326 → 0 | 0.176675 |
| Retain the existing Address | 52 | [467, 11] | 0 | 0.000074 |
| Extend to 104 while retaining the Address | 104 | [467, 11, 5] | 2 → 0 | 0.026203 |
| Insufficient-candidate control | 2 | [47] | 1 → 0 | 0.000228 |

The 104-integer row starts after applying the two existing axes.
These are observations in the recorded execution environment, not guarantees of depth for arbitrary sets or of runtime on a smartphone.

## Primary source and compatibility with the existing code

- Exponents: [GIMPS official list](https://www.mersenne.org/primes/), checked 2026-09-09.
  The list contains 52 known Mersenne primes, with maximum exponent 136,279,841 and maximum decimal length 41,024,320 digits. All exponents are fixed in the script.
  GIMPS marks the 52nd position as a provisional ranking.
- Existing work: [RSA-2048-Collapse-S11](https://github.com/Hironori-Sushiki-S11/RSA-2048-Collapse-S11).
- DOI for the existing work, supplied in the author's project handoff: [10.5281/zenodo.22378466](https://doi.org/10.5281/zenodo.22378466).
- Reference code: `ADAPTIVE_ADDRESS_SCALING_REPRODUCER.py`.
- Retrieved Git blob SHA: `f00a6160787750aae76e18791f8f31600ff32834`.
- Latest change commit for that reference file: `13b77310a1d0a8747bf5a60d75317a036ad03e2d`.
- The retrieved original is included as `upstream_reproducer.py`; its blob SHA was checked locally.

The selection score is the same as in the original: `(before-after, entropy, -b)`.
Entropy uses floating-point addition in the same order as the original. Rounding differences between implementations or accumulation orders can change very close tie decisions, even when entropies are mathematically equal. The code and environment are therefore recorded.
The corpus order in this experiment is fixed by the order of the listed exponents.

For LCM state L, the original evaluates `q=b/gcd(L,b)` and partitions each collision group by `(n-ref)//L mod q`, relative to a reference value ref in that group.
This implementation uses the following exactly equivalent expression without expanding `n = 2**p-1`:

```text
R_i = (pow(2, p_i, L*q) - 1) % (L*q)
relative_class = ((R_i - R_ref) % (L*q)) // L
```

Within a current collision group, `L | (n-ref)`, so the numerator is a multiple of L.
The identity `((n-ref) mod (L*q))/L = ((n-ref)/L) mod q` holds.
Evaluations with the same q are reused within each depth, and already distinguished singleton rows are excluded from partition evaluation.
Final Addresses are independently calculated as `[(pow(2,p,b)-1)%b for b in selected_axes]`.
These raw `n mod b` coordinates are distinct from IKERUSIKI's folded boundary distances.
The shifted coordinate `(n-1) mod b` and the raw coordinate `n mod b` have different values but induce identical equality partitions.

The original stops when no progress is made. As part of this extension, EXPAND doubles the candidate upper bound only when no progress is possible with the current candidates.
If collisions remain at the default cap of 8192, they are reported as unresolved.
An axis that adds no distinction is not appended to the adaptively selected Address.
The control uses p=1000 and p=1000+lcm(1,...,16). It expands the candidate upper bound 16→32→64, then selects 47.
This Blind Spot in the initial candidate range does not mean that the two represented integers are equal.

## Separating two effects

### 1. Resource savings from the generating formula and modular exponentiation

The exponent p is stored instead of the huge integer `2**p-1`. Residues are calculated through modular exponentiation with small moduli.
The exponent still requires storage and computation, but expansion of the p-bit integer body is unnecessary.
This is a property of modular exponentiation and is not counted as a speedup unique to Adaptive Address.

Only a separate, bounded comparison explicitly expanded integers, up to 1,048,576 bits.
It calculated residues for the same 511 bases at each size and verified that every residue matched.
Times are medians of five runs. The expanded-integer method includes integer construction time.

| p (bit length) | Modular exponentiation (seconds) | Expansion + residues (seconds) | Expanded integer object (bytes) |
|---:|---:|---:|---:|
| 32,768 | 0.000128 | 0.002352 | 4,396 |
| 131,072 | 0.000146 | 0.009442 | 17,504 |
| 1,048,576 | 0.000179 | 0.075535 | 139,836 |

Object sizes are measured with Python's sys.getsizeof, not as whole-process peak memory.
Runtime and peak memory for explicitly expanding the huge-integer corpus were not measured.

### 2. Adaptive selection and short Addresses

A comparison applied bases in ascending order from 2, using the same corpus, residue observations, and LCM partitioning.
This simple ascending baseline retains even bases that add no distinction to its Address.
It is not a comparison against the best existing algorithms or an assessment of novelty.

| Corpus and method | Stored axes | pow calls during selection | Median selection time (seconds) | Recompute all Addresses with selected axes (seconds) |
|---|---:|---:|---:|---:|
| Known 52 — adaptive selection | 2 | 30,652 | 0.029489 | 0.000045 |
| Known 52 — ascending | 28 | 600 | 0.000985 | 0.000422 |
| Huge 52 — adaptive selection | 1 | 26,572 | 0.176675 | 0.000285 |
| Huge 52 — ascending | 10 | 416 | 0.002565 | 0.002118 |

**In this experiment, adaptive selection produced shorter Addresses, but selection itself was slower than the ascending baseline.**
Candidate search and evaluation costs are not replaced by the cost of applying only the final two axes or single axis.
Recomputation with already selected axes is measured separately and excludes axis selection costs.
Selection time includes initial Address application, candidate evaluation, and partition refinement; it excludes the final independent verification and JSON output.
Evaluating the candidates requires tens of thousands of pow calls, even with small moduli; that search cost must be accounted for.
Selection skips singleton rows, so its pow count differs from the count for reconstructing Addresses for all rows.

## Validation scope

- Independently compared all final groups with groups formed from raw residue tuples.
- Checked at runtime that every update refines existing groups without merging them.
- Compared all candidate partitions and entropies between the original and extended code for one seed at each of the three existing bit sizes, plus a small Mersenne set.
- Matched the published EXPECTED axes and collision trajectories for the three existing cases.
- Verified 39,858 residue identities on small integers.
- Checked that unresolved collisions at the candidate cap are not reported as successful identification.
- Matched axes and collision trajectories across five repeated runs.

These are execution records from this ChatGPT session. Independent third-party reproduction of this experiment is a subsequent step.
Observed depth, implementation compatibility, mathematical equivalence of expressions, and the author's broader structural interpretation are distinct claims.
No new largest prime discovery, primality test, computation of unknown factors, or depth-two guarantee for arbitrary sets is claimed.
The completed status of the existing work is preserved.

## Reproduction on Android

The standard experiments run from the single file `AA_HUGE_INTEGER.py`, with no additional libraries or network connection required.
Open and run the file in a Python 3.9 or later environment on Android.
The output is `results.json` in the execution environment's working directory; use a writable directory.

If command execution is available:

```sh
python AA_HUGE_INTEGER.py
```

Use `--repeats 1` for a short run. Add `--benchmark` to include the bounded expanded-integer comparison.
The reported run used:

```sh
python AA_HUGE_INTEGER.py --repeats 5 --benchmark
python validate.py
```

The default run does not expand the huge integer bodies. For compatibility checks, validate.py explicitly generates the existing cases up to 131,072 bits and small Mersenne integers.
Runtime on a physical smartphone has not been measured. The recorded environment is Linux x86_64 / Python 3.12.14.

## Repository files and next steps

- `README.md`: English edition.
- `README.ja.md`: Japanese edition.
- `AA_HUGE_INTEGER.py`: standalone script including the exponent list.
- `results.json`: full Addresses for the 52/104-integer cases, steps, computation counts, timing, and environment.
- `validate.py` / `validation.json`: compatibility checks and their execution record.
- `upstream_reproducer.py`: the verified existing original.
- `SHA256SUMS.txt`: file checksums.

Dedicated repository: [IKERUSIKI-Adaptive-Address-Huge-Integers](https://github.com/Hironori-Sushiki-S11/IKERUSIKI-Adaptive-Address-Huge-Integers).
The code and experiment results are published. Next steps are reproduction on smartphones and other environments, and documentation of independent external validation.
The links above connect to the existing repository and its DOI. A new DOI dedicated to this experiment has not yet been obtained.
