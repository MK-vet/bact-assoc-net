---
name: validate-on-synth
description: Validate signal recovery on synthetic data with known ground truth
---

Run the pipeline on synthetic data (from `/generate-synth-data`) and verify signal recovery.

## Prerequisites
- Synthetic data must exist in `synth_data/` (run `/generate-synth-data` first)
- Each scenario has `synth_<name>.csv` + `ground_truth_<name>.json`

## Validation per scenario

1. **MI signal recovery**: Known MI pairs detected with p < 0.05; noise pairs have p > 0.05
2. **Phi edge cases**: Estimated phi within 0.1 of ground truth; CI coverage >= 90%
3. **OR zero cells**: Jeffreys-corrected OR is finite; direction matches ground truth
4. **CMH Simpson's**: Stratified result reverses marginal direction (as designed)
5. **PID decomposition**: |sum(components) - I_Total| < 1e-6
6. **NA patterns**: Results stable across MCAR levels; MNAR flagged appropriately
7. **Scale test**: Signal pair detected at all scales; runtime scales reasonably

## Multi-seed validation
For scenarios 1, 2, 3: re-run with seeds {42, 43, 44, 45, 46} and compute:
- Bias: mean(estimated - true)
- RMSE: sqrt(mean((estimated - true)^2))
- CI coverage: fraction of runs where true value falls in CI

## Report
Structured JSON-like output with scenario name, metric, expected, observed, PASS/FAIL.
