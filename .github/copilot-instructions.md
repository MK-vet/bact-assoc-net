# bact-assoc-net — AI Coding Assistant Instructions

## Context
This is **bact-assoc-net**, a Python tool for information-theoretic multi-layer network analysis of binary antimicrobial resistance (AMR) features. It computes MI, phi, OR, CMH, PID, and simplicial topology. It is published as an individual SoftwareX paper.

## CRITICAL: Scientific Correctness Rules (NEVER violate)
1. NEVER fabricate p-values, test statistics, CI, or effect sizes
2. NEVER generate fake citations (valid: Williams & Beer 2010, Schreiber 2000, Benjamini & Hochberg 1995)
3. NEVER convert NA to 0 (NA = "not tested", 0 = "tested and susceptible")
4. NEVER apply continuous-variable methods to binary (0/1) data
5. NEVER skip multiple testing correction (FDR/Benjamini-Hochberg)
6. NEVER assume sample independence without phylogenetic justification
7. NEVER add functionality from sister tools (see scope below)

## Scope Boundary (Anti-Salami)
- THIS tool: MI, phi, OR, CMH, PID, transfer entropy, simplicial topology
- NOT this tool: phylogenetic methods, MDR classification, clustering/consensus

## Code Conventions
- Python 3.10+, type hints with `X | Y` syntax
- Linter: ruff (check + format)
- Tests: pytest in `tests/`, parametrized
- Logging: `logging.getLogger(__name__)`, never `print()`
- Data: binary 0/1/NA, nullable Int8, CSV wide/long

## Dev Commands
```bash
pip install -e ".[dev]"
ruff check src/ && ruff format --check src/
python -m pytest tests/ -v --tb=short
python -m bactassocnet.cli --self-check
```

## Read-Only Files (DO NOT MODIFY)
- CITATION.cff, LICENSE, docs/anti_salami_checklist.json, examples/*
