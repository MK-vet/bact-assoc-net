"""YAML configuration — auto-detect CSV layers, configurable thresholds.

v1.1 changes
----------
- Adds schema_version and config validation with explicit handling of unknown keys.
- Unknown keys are reported (WARN by default) and ignored; enable strict mode to fail.

Strict mode triggers
--------------------
- YAML key: config_strict: true
- Env var: SSUIS_CONFIG_STRICT=1
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List
import dataclasses as _dc
import logging
import os
from typing import get_args, get_origin

import yaml

logger = logging.getLogger(__name__)

SCHEMA_VERSION = "1.1"
SUPPORTED_SCHEMA_VERSIONS = {"1.0", "1.1"}


def _is_dataclass_type(tp: Any) -> bool:
    try:
        return hasattr(tp, "__dataclass_fields__")
    except Exception:
        return False


def _dataclass_allowed_fields(dc_cls: Any) -> set[str]:
    return {f.name for f in _dc.fields(dc_cls)}


def _validate_and_filter(
    dc_cls: Any, raw: Any, prefix: str, unknown_paths: List[str]
) -> Any:
    """Return filtered raw structure matching dc_cls; accumulate unknown key paths."""
    if raw is None:
        raw = {}

    # dataclass instance
    if isinstance(raw, dict):
        allowed = _dataclass_allowed_fields(dc_cls)
        filtered: dict = {}
        for k, v in raw.items():
            if k not in allowed:
                unknown_paths.append(f"{prefix}.{k}" if prefix else str(k))
                continue
            # recurse when field type is dataclass or list of dataclasses
            f = next((ff for ff in _dc.fields(dc_cls) if ff.name == k), None)
            if f is None:
                filtered[k] = v
                continue
            tp = f.type
            if _is_dataclass_type(tp) and isinstance(v, dict):
                filtered[k] = _validate_and_filter(
                    tp, v, f"{prefix}.{k}" if prefix else k, unknown_paths
                )
            else:
                origin = get_origin(tp)
                args = get_args(tp)
                if (
                    origin in (list, List)
                    and args
                    and _is_dataclass_type(args[0])
                    and isinstance(v, list)
                ):
                    filtered_list = []
                    for i, item in enumerate(v):
                        if isinstance(item, dict):
                            filtered_list.append(
                                _validate_and_filter(
                                    args[0],
                                    item,
                                    f"{prefix}.{k}[{i}]" if prefix else f"{k}[{i}]",
                                    unknown_paths,
                                )
                            )
                        else:
                            filtered_list.append(item)
                    filtered[k] = filtered_list
                else:
                    filtered[k] = v
        return filtered

    # unexpected raw type
    return raw


@dataclass
class LayerSpec:
    name: str
    path: str
    id_column: str = "Strain_ID"
    format: str = "auto"  # auto | wide | long
    feature_column: str | None = None
    value_column: str | None = None


@dataclass
class InputSpec:
    align_mode: str = "union"  # union | intersection
    drop_samples_with_missing: bool = False
    max_missing_sample: float = 0.0
    max_missing_feature: float = 0.0


@dataclass
class FeatureQCSpec:
    min_prev: float = 0.01
    max_prev: float = 0.99
    max_missing_frac: float = 0.0


@dataclass
class ConditionalSpec:
    enabled: bool = True
    confounders: List[str] = field(
        default_factory=list
    )  # column names from confounder CSV(s)
    confounder_files: List[str] = field(
        default_factory=list
    )  # CSV paths with id_column + confounder cols


@dataclass
class InfoTheorySpec:
    surrogate_mi_enabled: bool = True  # cross-sectional surrogate-tested MI
    mi_n_surrogate: int = 200
    mi_alpha: float = 0.05
    pid_enabled: bool = True
    simplicial_enabled: bool = True
    simplicial_max_dim: int = 3


@dataclass
class NetworkSpec:
    phi_method: str = "adaptive"  # adaptive | fixed
    phi_percentile: int = 90
    phi_fixed: float = 0.3
    fdr_method: str = "fdr_bh"
    alpha: float = 0.05
    min_pairwise_n: int = 20  # minimum complete-case N for a pair


@dataclass
class Config:
    # schema + validation
    schema_version: str = SCHEMA_VERSION
    config_strict: bool = False

    layers: List[LayerSpec] = field(default_factory=list)
    input_dir: str = ""  # if set, auto-discover CSVs in dir
    output_dir: str = "network_results"
    id_column: str = "Strain_ID"
    input: InputSpec = field(default_factory=InputSpec)
    feature_qc: FeatureQCSpec = field(default_factory=FeatureQCSpec)
    conditional: ConditionalSpec = field(default_factory=ConditionalSpec)
    info_theory: InfoTheorySpec = field(default_factory=InfoTheorySpec)
    network: NetworkSpec = field(default_factory=NetworkSpec)
    n_jobs: int = -1
    seed: int = 42

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Config":
        raw = yaml.safe_load(Path(path).read_text()) or {}

        # schema version check
        schema_in = str(raw.get("schema_version", "1.0"))
        strict = bool(raw.get("config_strict", False)) or (
            os.environ.get("SSUIS_CONFIG_STRICT", "0") == "1"
        )
        if schema_in not in SUPPORTED_SCHEMA_VERSIONS:
            msg = f"Unsupported schema_version={schema_in!r}. Supported: {sorted(SUPPORTED_SCHEMA_VERSIONS)}"
            if strict:
                raise ValueError(msg)
            logger.warning(msg)

        unknown_paths: List[str] = []
        filtered = _validate_and_filter(
            cls, raw, prefix="", unknown_paths=unknown_paths
        )

        # instantiate
        layers = [
            LayerSpec(**(lyr or {})) for lyr in (filtered.pop("layers", []) or [])
        ]
        inp = InputSpec(**(filtered.pop("input", {}) or {}))
        fqc = FeatureQCSpec(**(filtered.pop("feature_qc", {}) or {}))
        cond = ConditionalSpec(**(filtered.pop("conditional", {}) or {}))
        it = InfoTheorySpec(**(filtered.pop("info_theory", {}) or {}))
        net = NetworkSpec(**(filtered.pop("network", {}) or {}))

        cfg = cls(
            layers=layers,
            input=inp,
            feature_qc=fqc,
            conditional=cond,
            info_theory=it,
            network=net,
            **filtered,
        )

        # attach validation payload for CLI/pipeline
        cfg._config_validation = {
            "schema_version_in": schema_in,
            "schema_version_effective": cfg.schema_version,
            "supported_schema_versions": sorted(SUPPORTED_SCHEMA_VERSIONS),
            "unknown_keys": sorted(set(unknown_paths)),
            "strict": strict,
            "status": "PASS"
            if (schema_in in SUPPORTED_SCHEMA_VERSIONS and not unknown_paths)
            else ("WARN" if not strict else "FAIL"),
        }

        if unknown_paths:
            msg = (
                f"Unknown config keys ignored ({len(set(unknown_paths))}): {sorted(set(unknown_paths))[:20]}"
                + (" ..." if len(set(unknown_paths)) > 20 else "")
            )
            if strict:
                raise ValueError(msg)
            logger.warning(msg)

        return cfg

    def resolve_layers(self) -> List[LayerSpec]:
        """If input_dir set and layers empty, auto-discover CSVs."""
        if self.layers:
            return self.layers
        if self.input_dir:
            p = Path(self.input_dir)
            return [
                LayerSpec(name=f.stem, path=str(f), id_column=self.id_column)
                for f in sorted(p.glob("*.csv"))
            ]
        return []

    def to_yaml(self, path: str | Path) -> None:
        import dataclasses as dc

        def ser(o):
            if dc.is_dataclass(o):
                return {k: ser(v) for k, v in dc.asdict(o).items()}
            if isinstance(o, list):
                return [ser(v) for v in o]
            return o

        Path(path).write_text(
            yaml.dump(ser(self), default_flow_style=False, sort_keys=False)
        )
