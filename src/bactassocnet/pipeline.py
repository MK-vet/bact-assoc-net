"""Pipeline: binary layer association network with defensible information-theory statistics.

Key design decisions in v1.0.0
- Missingness: complete-case for each pair; never impute NA→0.
- Long-format layers supported (ID,feature[,value]) with explicit observed/unobserved.
- Methodological fix: cross-sectional "transfer entropy" removed; replaced by surrogate-tested MI.
- Optional confounder control: CMH (2×2×K) and conditional MI I(X;Y|Z).
- Signed effect sizes: OR (Jeffreys CI) + sign (co-occur vs exclusive) alongside phi/MI.
- Topological sweep: Euler characteristic and triangle counts across association thresholds.
- run_manifest.json with input hashes for reproducibility.
"""
from __future__ import annotations
from . import __version__

import json
import logging
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import chi2

from .config import Config, LayerSpec
from .io import load_layer_csv, layer_coverage, qc_binary_features, sha256_file
from .information_theory.core import chi2_phi, mutual_information, adaptive_phi_threshold
from .information_theory.advanced import surrogate_mi
from .innovation_v3 import glasso_ising_proxy, topology_persistence_from_sweep, multilayer_dispersion

logger = logging.getLogger(__name__)


def _save(df: pd.DataFrame, p: str) -> None:
    num = df.select_dtypes(include=[np.number]).columns
    out = df.copy()
    if len(num):
        out[num] = out[num].round(6)
    out.to_csv(p, index=False)
    logger.info("  → %s (%d rows)", p, len(out))


def _or_jeffreys_ci(a: int, b: int, c: int, d: int, alpha: float = 0.05) -> Tuple[float, float, float]:
    """Odds ratio with Jeffreys prior (0.5) and Wald CI on log(OR)."""
    z = 1.959963984540054  # norm.ppf(1-alpha/2)
    a2, b2, c2, d2 = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orr = (d2 * a2) / (b2 * c2)
    se = np.sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
    lo = np.exp(np.log(orr) - z * se)
    hi = np.exp(np.log(orr) + z * se)
    return float(orr), float(lo), float(hi)


def _contingency_2x2(x: pd.Series, y: pd.Series) -> Tuple[np.ndarray, int]:
    m = x.notna() & y.notna()
    xv = x[m].astype(int).to_numpy()
    yv = y[m].astype(int).to_numpy()
    n = int(len(xv))
    a = int(np.sum((xv == 0) & (yv == 0)))
    b = int(np.sum((xv == 0) & (yv == 1)))
    c = int(np.sum((xv == 1) & (yv == 0)))
    d = int(np.sum((xv == 1) & (yv == 1)))
    return np.array([[a, b], [c, d]], dtype=int), n


def _cmh_2x2xk(x: pd.Series, y: pd.Series, z: pd.Series) -> Dict[str, float]:
    """Cochran–Mantel–Haenszel test for binary x,y stratified by categorical z."""
    m = x.notna() & y.notna() & z.notna()
    if int(m.sum()) < 20:
        return {"CMH": np.nan, "P": 1.0, "OR_MH": np.nan, "N": int(m.sum()), "K": 0}
    x = x[m].astype(int).to_numpy()
    y = y[m].astype(int).to_numpy()
    z = z[m].astype(str).to_numpy()

    strata = {}
    for xi, yi, zi in zip(x, y, z):
        strata.setdefault(zi, []).append((xi, yi))

    num = 0.0
    den = 0.0
    s = 0.0
    v = 0.0
    K = 0
    for key, pairs in strata.items():
        pairs = np.array(pairs, dtype=int)
        if len(pairs) < 5:
            continue
        K += 1
        xv = pairs[:, 0]; yv = pairs[:, 1]
        a = np.sum((xv == 0) & (yv == 0))
        b = np.sum((xv == 0) & (yv == 1))
        c = np.sum((xv == 1) & (yv == 0))
        d = np.sum((xv == 1) & (yv == 1))
        n = a + b + c + d
        if n <= 1:
            continue
        row1 = a + b
        row2 = c + d
        col1 = a + c
        col2 = b + d
        exp_a = row1 * col1 / n
        var_a = (row1 * row2 * col1 * col2) / (n**2 * (n - 1) + 1e-12)

        s += (a - exp_a)
        v += var_a

        num += (a * d) / max(n, 1)
        den += (b * c) / max(n, 1)

    if K == 0 or v <= 0:
        return {"CMH": np.nan, "P": 1.0, "OR_MH": np.nan, "N": int(m.sum()), "K": K}

    cmh = (s**2) / v
    p = float(1.0 - chi2.cdf(cmh, df=1))
    or_mh = float(num / (den + 1e-12))
    return {"CMH": float(cmh), "P": p, "OR_MH": or_mh, "N": int(m.sum()), "K": K}


def _cmi_binary(x: pd.Series, y: pd.Series, z: pd.Series) -> float:
    """Conditional MI I(X;Y|Z) for binary X,Y and discrete Z."""
    m = x.notna() & y.notna() & z.notna()
    x = x[m].astype(int).to_numpy()
    y = y[m].astype(int).to_numpy()
    z = z[m].astype(str).to_numpy()
    n = len(x)
    if n < 20:
        return float("nan")
    cmi = 0.0
    for level in np.unique(z):
        idx = (z == level)
        if idx.sum() < 10:
            continue
        cmi += (idx.mean()) * mutual_information(x[idx], y[idx])
    return float(cmi)


def _simplicial_sweep(G: nx.Graph, thresholds: List[float]) -> pd.DataFrame:
    """Count V/E/T and Euler characteristic across phi thresholds."""
    rows = []
    for thr in thresholds:
        H = nx.Graph(((u, v, d) for u, v, d in G.edges(data=True) if abs(d.get("phi", 0.0)) >= thr))
        V = H.number_of_nodes()
        E = H.number_of_edges()
        # triangles
        T = sum(nx.triangles(H).values()) // 3
        chi = V - E + T
        rows.append({"phi_threshold": thr, "V": V, "E": E, "Triangles": int(T), "Euler": int(chi)})
    return pd.DataFrame(rows)


class Pipeline:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.out = Path(cfg.output_dir)
        self.out.mkdir(parents=True, exist_ok=True)

    def _write_manifest(self, used_layers: List[LayerSpec], extra_inputs: List[str]) -> None:
        import platform
        import sys
        inp = [{"name": sp.name, "path": sp.path, "sha256": sha256_file(sp.path)} for sp in used_layers]
        for p in extra_inputs:
            inp.append({"name": "confounder", "path": p, "sha256": sha256_file(p)})
        import dataclasses as dc
        try:
            config_snapshot = dc.asdict(self.cfg)
        except Exception:
            config_snapshot = None
        manifest = {
            "tool": "bact-assoc-net",
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "seed": self.cfg.seed,
            "tool_version": __version__,
            "inputs": inp,
        }
        manifest["config_snapshot"] = config_snapshot
        try:
            import numpy
            import pandas
            manifest["package_versions"]={"numpy":numpy.__version__,"pandas":pandas.__version__}
            try:
                import scipy; manifest["package_versions"]["scipy"]=scipy.__version__
            except Exception: pass
            try:
                import networkx; manifest["package_versions"]["networkx"]=networkx.__version__
            except Exception: pass
        except Exception:
            pass
        (self.out / "run_manifest.json").write_text(json.dumps(manifest, indent=2))

    def run(self) -> Dict[str, pd.DataFrame]:
        cfg = self.cfg
        out = self.out

        # 1) global ID universe
        id_sets = []
        for sp in cfg.layers:
            tmp = pd.read_csv(sp.path, usecols=[sp.id_column])
            tmp[sp.id_column] = tmp[sp.id_column].astype(str)
            id_sets.append(set(tmp[sp.id_column].unique()))
        all_ids = sorted(set.intersection(*id_sets)) if cfg.input.align_mode.lower() == "intersection" else sorted(set.union(*id_sets))
        if not all_ids:
            raise ValueError("No IDs across inputs.")
        id_col = cfg.layers[0].id_column

        # 2) load layers
        layers = {}
        cov_rows = []
        qc_rows = []
        for sp in cfg.layers:
            lr = load_layer_csv(sp.name, sp.path, id_col=sp.id_column, all_ids=all_ids,
                                feature_col=sp.feature_column if sp.format == "long" else sp.feature_column,
                                value_col=sp.value_column)
            df = lr.data
            # Keep only numeric/binary columns for network
            df = df.select_dtypes(include=[np.number])
            fqc = cfg.feature_qc
            df_qc, rep = qc_binary_features(df, observed=lr.observed,
                                            min_prev=fqc.min_prev, max_prev=fqc.max_prev, max_missing_frac=fqc.max_missing_frac)
            rep.insert(0, "Layer", sp.name)
            qc_rows.append(rep)
            # sample/feature missingness trimming (no imputation)
            if cfg.input.drop_samples_with_missing:
                obs = lr.observed & df_qc.notna().any(axis=1)
                obs = obs & (df_qc.isna().mean(axis=1) <= cfg.input.max_missing_sample)
                df_qc = df_qc.loc[obs]
            df_qc = df_qc.loc[:, (df_qc.isna().mean(axis=0) <= cfg.input.max_missing_feature)]
            layers[sp.name] = df_qc
            cov_rows.append(layer_coverage(lr))

        cov_df = pd.DataFrame(cov_rows)
        _save(cov_df, str(out / "layer_coverage.csv"))
        qc_df = pd.concat(qc_rows, ignore_index=True) if qc_rows else pd.DataFrame()
        if not qc_df.empty:
            _save(qc_df, str(out / "feature_qc.csv"))

        # 3) merge layers (outer join on union IDs)
        # Use NA to represent missing features in missing strains; association is computed pairwise with complete-case.
        merged = pd.DataFrame(index=all_ids)
        layer_cols = {}
        for lname, df in layers.items():
            pref = lname + "::"
            cols = [pref + c for c in df.columns]
            merged = merged.join(df.rename(columns=dict(zip(df.columns, cols))), how="left")
            layer_cols[lname] = cols

        merged.index.name = id_col

        # 4) load confounders (optional)
        conf = None
        extra_inp = []
        if cfg.conditional.enabled and cfg.conditional.confounder_files and cfg.conditional.confounders:
            parts = []
            for p in cfg.conditional.confounder_files:
                extra_inp.append(p)
                cdf = pd.read_csv(p)
                if id_col not in cdf.columns:
                    raise ValueError(f"Confounder file {p} missing '{id_col}'.")
                cdf[id_col] = cdf[id_col].astype(str)
                avail = [c for c in cfg.conditional.confounders if c in cdf.columns]
                if not avail:
                    continue
                parts.append(cdf.set_index(id_col)[avail])
            if parts:
                conf = pd.concat(parts, axis=1)
                conf = conf.reindex(all_ids)
                conf = conf.loc[:, ~conf.columns.duplicated()]
            else:
                conf = None

        # 5) pairwise associations
        net = cfg.network
        it = cfg.info_theory

        cols = merged.columns.tolist()
        edges = []
        for i, j in combinations(range(len(cols)), 2):
            c1, c2 = cols[i], cols[j]
            x = merged[c1]
            y = merged[c2]
            tab, n = _contingency_2x2(x, y)
            if n < net.min_pairwise_n:
                continue

            a, b = int(tab[0, 0]), int(tab[0, 1])
            c_, d = int(tab[1, 0]), int(tab[1, 1])

            st = chi2_phi(x, y)
            phi = st["phi"]
            p = st["p"]
            mxy = x.notna() & y.notna()
            mi = float(mutual_information(x[mxy].astype(int).to_numpy(), y[mxy].astype(int).to_numpy()))
            # Surrogate-tested MI is optional because it can be expensive for wide feature sets.
            if it.surrogate_mi_enabled:
                smi = surrogate_mi(x, y, n_surrogates=it.mi_n_surrogate, seed=cfg.seed)
                mi_z = smi["Z"]
                mi_p = smi["P"]
            else:
                mi_z = float("nan")
                mi_p = 1.0
            orr, or_lo, or_hi = _or_jeffreys_ci(a, b, c_, d)

            row = {
                "Feature_A": c1,
                "Feature_B": c2,
                "N": n,
                "phi": phi,
                "p": p,
                "MI": mi,
                "MI_Z": mi_z,
                "MI_P": mi_p,
                "OR": orr,
                "OR_CI_Lo": or_lo,
                "OR_CI_Hi": or_hi,
                "Effect": "co-occur" if orr > 1 else "exclusive",
                "a00": a, "a01": b, "a10": c_, "a11": d,
            }

            # confounder control (CMH + CMI) for each named confounder column
            if conf is not None:
                for zcol in conf.columns:
                    cmh = _cmh_2x2xk(x, y, conf[zcol])
                    row[f"CMH_{zcol}"] = cmh["CMH"]
                    row[f"CMH_P_{zcol}"] = cmh["P"]
                    row[f"OR_MH_{zcol}"] = cmh["OR_MH"]
                    row[f"CMI_{zcol}"] = _cmi_binary(x, y, conf[zcol])

            edges.append(row)

        edges_df = pd.DataFrame(edges)
        if edges_df.empty:
            raise ValueError("No edges computed (check min_pairwise_n / QC filters).")
        _save(edges_df.sort_values(["MI_P", "p"]), str(out / "associations_all.csv"))

        # 6) thresholded network
        # Determine phi threshold either adaptively (percentile) or fixed.
        if net.phi_method.lower() == "adaptive":
            thr_phi = adaptive_phi_threshold(edges_df["phi"].to_numpy(), method="percentile", percentile=int(net.phi_percentile))
        else:
            thr_phi = float(net.phi_fixed)
        alpha = float(net.alpha)
        keep = (edges_df["N"] >= net.min_pairwise_n) & (edges_df["p"] <= alpha) & (edges_df["phi"].abs() >= thr_phi)
        E = edges_df[keep].copy()
        _save(E.sort_values(["p", "MI_P"]), str(out / "associations_filtered.csv"))

        G = nx.Graph()
        for _, r in E.iterrows():
            G.add_edge(r["Feature_A"], r["Feature_B"], phi=float(r["phi"]), p=float(r["p"]), MI=float(r["MI"]), OR=float(r["OR"]))

        # Graph summaries
        summ = [{
            "Nodes": G.number_of_nodes(),
            "Edges": G.number_of_edges(),
            "Components": nx.number_connected_components(G) if G.number_of_nodes() else 0,
            "AvgDegree": float(np.mean([d for _, d in G.degree()])) if G.number_of_nodes() else 0.0,
        }]
        _save(pd.DataFrame(summ), str(out / "network_summary.csv"))

        # 7) topological sweep (Euler characteristic)
        sweep = _simplicial_sweep(G, thresholds=[0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
        _save(sweep, str(out / "simplicial_sweep.csv"))

        # 8) layer-level "flow": aggregate MI across layer pairs
        rows = []
        for la, lb in combinations(layer_cols.keys(), 2):
            m = edges_df[(edges_df["Feature_A"].str.startswith(la + "::") & edges_df["Feature_B"].str.startswith(lb + "::")) |
                         (edges_df["Feature_A"].str.startswith(lb + "::") & edges_df["Feature_B"].str.startswith(la + "::"))]
            if m.empty:
                continue
            rows.append({"Layer_A": la, "Layer_B": lb, "N_pairs": int(len(m)),
                         "MI_mean": float(m["MI"].mean()), "MI_median": float(m["MI"].median()),
                         "phi_mean_abs": float(m["phi"].abs().mean())})
        if rows:
            _save(pd.DataFrame(rows), str(out / "layer_flow.csv"))


        # 9) v3 innovations: Ising-like deconvolution and multilayer dispersion
        try:
            direct = glasso_ising_proxy(E if not E.empty else edges_df, feature_matrix=merged, thr=0.0, alpha=0.02)
            if not direct.empty:
                _save(direct, str(out / 'ising_direct_edges.csv'))
            sweep_p = topology_persistence_from_sweep(sweep)
            if not sweep_p.empty:
                _save(sweep_p, str(out / 'topology_persistence.csv'))
            mld = multilayer_dispersion(edges_df)
            if not mld.empty:
                _save(mld, str(out / 'multilayer_dispersion.csv'))
        except Exception as e:
            logger.warning('  v3 innovation outputs skipped: %s', e)

        # Persist config + manifest
        cfg.to_yaml(str(out / "config_used.yaml"))
        self._write_manifest(cfg.layers, extra_inp)
        logger.info("Done — results in %s", out)

        return {
            "coverage": cov_df,
            "feature_qc": qc_df,
            "associations_all": edges_df,
            "associations_filtered": E,
            "simplicial_sweep": sweep,
        }
