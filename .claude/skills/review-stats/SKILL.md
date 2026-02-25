---
name: review-stats
description: Review statistical methods for correctness — degrees of freedom, null distributions, FDR, edge cases
---

Perform a thorough review of statistical method implementations in this tool. This is a code review focused on mathematical/statistical correctness.

## Checklist

For each statistical method in `src/bactassocnet/`:

1. **Degrees of freedom**: Are df correctly computed? (chi2 df=1 for 2x2, CMH df=1)
2. **Null distribution**: Is the null model correct? (chi2 for phi, permutation for surrogate MI)
3. **Multiple testing correction**: Is FDR (Benjamini-Hochberg) applied after pairwise tests?
4. **Edge cases**:
   - Zero cells in contingency table → Jeffreys prior or Fisher exact?
   - All-NA column → graceful skip (not crash)?
   - Single unique value → phi=0, MI=0?
   - n < min_n → skip pair?
5. **PID constraint**: Does redundancy + unique_A + unique_B + synergy == I_Total?
6. **Transfer entropy guard**: Is cross-sectional TE blocked/warned?
7. **CI correctness**: Fisher-z for phi CI, Wald on log(OR) for OR CI?

## Key files to review:
- `src/bactassocnet/information_theory/core.py` — chi2_phi, mutual_information
- `src/bactassocnet/information_theory/advanced.py` — surrogate_mi, PID, TE
- `src/bactassocnet/pipeline.py` — _or_jeffreys_ci, _cmh_2x2xk, FDR application
- `src/bactassocnet/innovation_v3.py` — topology, Ising

Report findings as a structured checklist with PASS/WARN/FAIL per item.
