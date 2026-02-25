---
name: run-tests
description: Run the test suite for bact-assoc-net
---

Run the full test suite with verbose output and short tracebacks.

```bash
python -m pytest tests/ -v --tb=short
```

Report:
1. Total tests passed/failed/skipped
2. Any failures with file:line and assertion message
3. If all tests pass, confirm with "All tests PASS"
