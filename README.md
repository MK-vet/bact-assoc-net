# bact-assoc-net

[![CI](https://github.com/MK-vet/bact-assoc-net/actions/workflows/ci.yaml/badge.svg)](https://github.com/MK-vet/bact-assoc-net/actions/workflows/ci.yaml)
[![PyPI](https://img.shields.io/pypi/v/bact-assoc-net.svg)](https://pypi.org/project/bact-assoc-net/)
[![Python](https://img.shields.io/pypi/pyversions/bact-assoc-net.svg)](https://pypi.org/project/bact-assoc-net/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)


Information-theoretic multi-layer network analysis with transfer entropy, partial information decomposition, and simplicial topology.

## Novel Contributions

1. **Transfer entropy** with surrogate permutation testing — directional information flow: TE(gene→phenotype) ≠ TE(pheno→gene) (Schreiber 2000)
2. **Partial Information Decomposition** — decomposes I(S₁,S₂;T) into unique, redundant, synergistic components (Williams & Beer 2010)
3. **Simplicial complex construction** from significant associations — Betti numbers and Euler characteristic capture higher-order topology beyond pairwise graphs
4. **Layer-resolved information flow** — directed MI/TE between entire feature categories

## Installation

```bash
pip install bact-assoc-net
pip install bact-assoc-net[gui]         # interactive dashboard (marimo)
```

## Interactive dashboard (marimo)

```bash
bact-assoc-net-dashboard
```

Edit mode:

```bash
bact-assoc-net-dashboard --edit
```

## Quick Start

```bash
bact-assoc-net config.yaml -v
```

### Python API

```python
from bactassocnet.config import Config
from bactassocnet.pipeline import Pipeline

cfg = Config.from_yaml("config.yaml")
results = Pipeline(cfg).run()
```

## Configuration

Layers can be specified explicitly or auto-discovered from a directory:

```yaml
# Option 1: explicit
layers:
  - name: MIC
    path: data/MIC.csv
  - name: AMR_genes
    path: data/AMR_genes.csv

# Option 2: auto-discover all CSVs in directory
input_dir: data/
```

See `examples/config.yaml` for a complete template.

## Outputs

| File | Description |
|------|-------------|
| `pairwise_associations.csv` | Chi², phi, NMI for all feature pairs |
| `association_network.graphml` | Network above adaptive threshold |
| `transfer_entropy.csv` | TE forward/reverse/net for all pairs |
| `pid.csv` | PID decomposition per source pair × target |
| `simplicial_complex.csv` | Betti numbers and Euler characteristic |
| `layer_info_flow.csv` | Directed TE between feature categories |
| `mutual_exclusivity.csv` | Feature pairs that never co-occur |

## Testing

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

## License

MIT
