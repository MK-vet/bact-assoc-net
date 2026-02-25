---
name: reliability-check
description: Run reliability preflight QC and quality gate check
---

Run the reliability preflight to verify data quality before analysis.

```bash
python -m bactassocnet.cli --config examples/config.yaml --reliability-only
```

## What to check:
1. ID overlap across layers (if multi-layer)
2. NA proportion sanity (flag columns with >50% missing)
3. Degenerate features (zero variance, all-NA)
4. Quality gate status (PASS/WARN/FAIL)

Report the quality gate result and any warnings. If FAIL, list the failing checks with remediation suggestions.
