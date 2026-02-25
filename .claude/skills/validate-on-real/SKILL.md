---
name: validate-on-real
description: Validate pipeline on real data from examples/
---

Run the pipeline on real example data and verify basic invariants.

## Steps

1. **Run pipeline**:
   ```bash
   python -m bactassocnet.cli --config examples/config.yaml
   ```

2. **Check invariants**:
   - No crash (exit code 0)
   - Output files exist in the configured output directory
   - `run_manifest.json` generated with input SHA-256 hashes
   - `config_used.yaml` snapshot exists
   - No NA values were silently converted to 0 (check log for warnings)
   - All p-values are in [0, 1]
   - All phi values are in [-1, 1]
   - All MI values are >= 0
   - OR values are > 0
   - FDR-corrected q-values are >= raw p-values
   - Sample counts in output match input

3. **Report**: Structured checklist of invariants PASS/FAIL.
