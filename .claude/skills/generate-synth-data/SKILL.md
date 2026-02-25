---
name: generate-synth-data
description: Generate synthetic data with known ground truth for signal recovery testing
---

Generate synthetic binary matrices with known statistical properties for validating bact-assoc-net analyses.

## Scenarios to generate

Each scenario produces: `synth_<scenario>.csv` + `ground_truth_<scenario>.json` + `synth_manifest.json`

1. **MI signal recovery** (seed=42, n=200, p=20):
   - 2 feature pairs with known MI > 0 (correlated via XOR with noise)
   - Remaining features independent (Bernoulli p=0.5)
   - Ground truth: known MI values, expected p < 0.05 for signal pairs

2. **Phi coefficient edge cases** (seed=43, n=200):
   - Feature pairs with known phi in {-1, -0.5, 0, 0.5, 1.0}
   - Include sparse pairs (low prevalence) and balanced pairs
   - Ground truth: exact phi values and expected CI coverage

3. **OR with zero cells** (seed=44, n=100):
   - Pairs where one cell in the 2x2 table is 0 → test Jeffreys prior (+0.5)
   - Ground truth: expected OR direction and CI

4. **CMH Simpson's paradox** (seed=45, n=300, K=3 strata):
   - Marginal association positive, but within each stratum negative (or vice versa)
   - Ground truth: marginal vs stratified direction

5. **PID decomposition** (seed=46, n=500):
   - 3 binary variables with known redundancy, unique, synergy
   - Ground truth: exact PID components summing to I_Total

6. **NA patterns** (seed=47, n=200, p=20):
   - MCAR at 5%, 20%, 50%; MNAR on high-prevalence features
   - Ground truth: complete-data statistics for comparison

7. **Scale test** (seed=48, p in {50, 200, 1000, 5000}, n=100):
   - Fixed signal (2 correlated pairs) + noise features
   - Ground truth: signal recovery should be scale-invariant

## Implementation approach

Write a Python script using numpy/pandas with deterministic seeds. Save all outputs to a `synth_data/` directory. The `ground_truth_*.json` files encode the known parameters for each scenario.
