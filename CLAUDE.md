# bact-assoc-net — Claude Code Context

## Project Identity

- **Tool**: bact-assoc-net v1.0.0
- **Package**: `bactassocnet` (in `src/bactassocnet/`)
- **Author**: Maciej Kochanowski, National Veterinary Research Institute (PIWet-PIB), Pulawy, Poland
- **Domain**: Antimicrobial resistance (AMR) epidemiology, *Streptococcus suis* and related bacteria
- **Purpose**: Information-theoretic multi-layer network analysis for binary features with transfer entropy, partial information decomposition (PID), and simplicial topology
- **Publication**: Individual SoftwareX paper — do NOT mix scope with sister tools

## Sister Tools (Anti-Salami Boundaries)

| Tool | Scope | THIS tool does NOT do this |
|------|-------|---------------------------|
| **bact-assoc-net** (THIS) | MI, phi, OR, CMH, PID, TE, simplicial topology | — |
| bact-phylo-trait | Pagel lambda, Fritz-Purvis D, Fitch ASR, phylo-logistic | NO phylogenetic comparative methods |
| bact-mdr-profiler | PC-algorithm causal, hypergraph MDR, Shapley, EVPI | NO causal discovery, NO MDR classification |
| bact-trait-cluster | K-Modes consensus, CKA fusion, NVI, TDA, SHAP | NO clustering, NO consensus fusion |

Reference: `docs/anti_salami_checklist.json`

## Architecture

```
src/bactassocnet/
  config.py          — YAML config v1.1, dataclass schema, strict mode
  pipeline.py        — Main orchestrator: load → QC → associations → topology → PID → manifest
  cli.py             — CLI entry: --config, --self-check, --benchmark, --reliability-*
  io/loader.py       — CSV loading, wide/long auto-detect, NA-preserving Int8
  information_theory/
    core.py          — chi2_phi, mutual_information, adaptive_phi_threshold
    advanced.py      — surrogate_mi, PID (Williams & Beer), transfer entropy
  innovation_v3.py   — Ising PLM-L1, topology persistence, multilayer dispersion
  reliability/core.py — Preflight QC, consistency checks, marimo warnings
  selfcheck.py       — Internal validation
  benchmark.py       — Synthetic benchmark
  gui/               — Marimo interactive app
  dashboards/        — Plotly dashboard
```

**Key patterns**: config-driven pipeline, NA-preserving (never NA→0), provenance via `run_manifest.json` + `config_used.yaml`, deterministic seeds.

## Development Quickstart

```bash
pip install -e ".[dev]"
ruff check src/
ruff format --check src/
python -m pytest tests/ -v --tb=short
python -m bactassocnet.cli --self-check
python -m bactassocnet.cli --benchmark
python -m bactassocnet.cli --config examples/config.yaml
```

## CRITICAL: Scientific Correctness Rules

**These rules are ABSOLUTE. Violating any of them produces scientifically invalid results.**

1. **NEVER fabricate** p-values, test statistics, confidence intervals, or effect sizes. All must come from actual computation on actual data.
2. **NEVER generate fake citations**. Valid references for this tool:
   - Williams & Beer (2010) — PID I_min
   - Schreiber (2000) — Transfer entropy
   - Benjamini & Hochberg (1995) — FDR correction
   - Cochran (1954) / Mantel & Haenszel (1959) — CMH test
3. **NEVER convert NA to 0**. NA means "not tested"; 0 means "tested and susceptible". These are biologically different.
4. **NEVER apply continuous-variable methods** to binary (0/1) data. All methods must handle categorical/binary inputs.
5. **NEVER skip multiple testing correction** when computing pairwise statistics. Always apply FDR (Benjamini-Hochberg) or equivalent.
6. **NEVER assume independence** of samples with shared phylogenetic origin without explicit justification.
7. **NEVER add functionality belonging to a sister tool** (see anti-salami table above).

## Statistical Constraints — THIS TOOL ONLY

- **MI**: Computed from joint probability table; must be non-negative; significance via surrogate permutation null distribution
- **Phi coefficient**: Complete-case only; Fisher-z transform for CI; Fisher exact test fallback for small cells
- **Odds ratio**: Jeffreys prior (add 0.5 to all cells); Wald CI on log(OR)
- **CMH**: Requires K >= 1 strata with n >= 5 per stratum; report K and per-stratum N
- **PID**: Williams & Beer I_min decomposition; redundancy + unique_A + unique_B + synergy MUST equal I_Total
- **Transfer entropy**: Cross-sectional data does NOT support temporal TE — use surrogate-tested MI instead
- **Simplicial topology**: Euler characteristic chi = V - E + T; persistence sweep across phi thresholds

## Data Format Constraints

- Input: CSV with binary features (0 = absent, 1 = present, NA = not tested)
- Internal dtype: nullable Int8 (`pd.Int8Dtype()`)
- Formats: wide (samples x features) or long (ID, feature, value) — auto-detected
- NA handling: preserve throughout; complete-case per pair for statistics
- Output: CSV with rounded floats (6 decimal places), provenance JSON

## Testing Conventions

- **Real data**: Run on `examples/` CSV files with example configs
- **Synthetic data**: Generated with known ground truth (`ground_truth.json`), deterministic seed
- **Signal recovery**: Verify tool recovers known MI, phi, OR from synthetic data
- **Multi-seed**: Run with multiple seeds to check bias/RMSE
- **Edge cases**: Zero cells, all-NA columns, single-value columns, n < min_n
- **Reproducibility**: Same seed + config → identical output (bit-for-bit)

## Publication Context (SoftwareX)

- CITATION.cff, LICENSE, README.md exist — **DO NOT MODIFY** these files
- Scope boundary defined in `docs/anti_salami_checklist.json` — **DO NOT MODIFY**
- Reproducibility artifacts: `run_manifest.json` (input SHA-256), `config_used.yaml` (config snapshot)
- All results must be fully reproducible given seed + config + input data

## Code Style

- Python 3.10+ (type hints with `X | Y` union syntax)
- Linter: ruff (check + format)
- Logging: `logging.getLogger(__name__)`, never print()
- Tests: pytest in `tests/`, parametrized where possible
