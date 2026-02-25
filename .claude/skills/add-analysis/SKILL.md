---
name: add-analysis
description: Guided workflow to add a new analysis step to the pipeline
---

Walk through the checklist for adding a new analysis step to bact-assoc-net.

## Pre-flight

1. **Scope check**: Does this analysis belong to bact-assoc-net? Check against `docs/anti_salami_checklist.json`.
   - If it involves phylogenetics → belongs to bact-phylo-trait
   - If it involves MDR classification or causal discovery → belongs to bact-mdr-profiler
   - If it involves clustering or consensus → belongs to bact-trait-cluster

2. **Method validation**: Is the statistical method appropriate for binary data?

## Implementation checklist

3. **Config**: Add any new parameters to `src/bactassocnet/config.py` Config dataclass
4. **Pipeline**: Add the step to `src/bactassocnet/pipeline.py` in the correct position
5. **Tests**: Add tests in `tests/` covering:
   - Signal recovery on synthetic data with known ground truth
   - Edge cases (NA, zero cells, degenerate input)
   - Reproducibility (same seed → same result)
6. **Validate config**: Run `/validate-config` to verify schema still parses
7. **Run tests**: Run `/run-tests` to verify no regressions
8. **Selfcheck**: Run `/selfcheck` to verify internal consistency
9. **Reliability**: Run `/reliability-check` on example data
10. **Cross-tool check**: Run `/cross-tool-check` to verify no scope violation

## Post-implementation

11. **Update README.md** if the analysis is user-facing
12. **Verify provenance**: Check that `run_manifest.json` includes the new step
