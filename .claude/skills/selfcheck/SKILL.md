---
name: selfcheck
description: Run internal validation self-check
---

Run the built-in self-check that validates internal consistency.

```bash
python -m bactassocnet.cli --self-check
```

Report:
1. Self-check status (PASS/FAIL)
2. Individual check results from the JSON output
3. If any checks fail, list them with details
