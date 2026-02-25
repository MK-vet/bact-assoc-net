---
name: full-ci
description: Run full CI pipeline — lint, tests, selfcheck, placeholder scan
---

Run the full CI validation pipeline in order. Stop on first failure.

## Steps

1. **Lint**: `ruff check src/ && ruff format --check src/`
2. **Tests**: `python -m pytest tests/ -v --tb=short`
3. **Selfcheck**: `python -m bactassocnet.cli --self-check`
4. **Placeholder scan**: Search for leftover placeholders that should not exist in production code:
   - `grep -rn '<PKG>\|<TOOL\|{PKG}\|TODO\|FIXME\|HACK\|XXX' src/ tests/ CLAUDE.md .claude/ .github/ --include='*.py' --include='*.md' --include='*.json' --include='*.yaml'`

Report each step as PASS/FAIL. If all pass, confirm "Full CI PASS".
