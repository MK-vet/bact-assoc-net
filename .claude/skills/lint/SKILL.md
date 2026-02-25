---
name: lint
description: Run ruff linter and formatter check for bact-assoc-net
---

Run ruff linting and format check on the source code.

```bash
ruff check src/
ruff format --check src/
```

Report:
1. Number of lint errors/warnings (if any) with file:line
2. Format check result (all formatted / N files would be reformatted)
3. If clean, confirm "Lint PASS — no issues found"
