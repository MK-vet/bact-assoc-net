---
name: validate-config
description: Validate a YAML pipeline configuration file against schema v1.1
---

Validate the given YAML config file against the bact-assoc-net config schema.

If `$ARGUMENTS` is provided, use it as the config path. Otherwise use `examples/config.yaml`.

```bash
python -c "
from bactassocnet.config import Config
import sys
path = sys.argv[1] if len(sys.argv) > 1 else 'examples/config.yaml'
cfg = Config.from_yaml(path)
print(f'Config loaded OK: schema_version={cfg.schema_version}')
print(f'Layers: {len(cfg.layers)}')
print(f'Output dir: {cfg.output_dir}')
" $ARGUMENTS
```

Report:
1. Whether config parsed successfully
2. Schema version detected
3. Number of layers configured
4. Any unknown keys or warnings
