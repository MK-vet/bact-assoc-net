import marimo
from pathlib import Path


app = marimo.App(width="full")


@app.cell
def _():
    import os
    import sys
    import json
    import time
    import glob
    import subprocess
    from pathlib import Path

    import marimo as mo
    import pandas as pd
    import yaml

    return Path, glob, json, mo, os, pd, subprocess, sys, time, yaml


@app.cell
def _(mo):
    mo.md(
        """
        # bact-assoc-net — interactive dashboard (marimo)

        Use this dashboard to generate and explore association networks.
        It supports config editing, non-blocking runs, and scalable table previews.
        """
    )
    return


@app.cell
def _(mo):
    data_root = mo.ui.text(value=str(Path.cwd()), label="Data root")
    out_root = mo.ui.text(
        value=str(Path.cwd() / "bactassocnet_out"), label="Output directory"
    )

    mic = mo.ui.text(value="MIC.csv", label="MIC CSV (optional)")
    amr = mo.ui.text(value="AMR_genes.csv", label="AMR genes CSV (optional)")
    vir = mo.ui.text(value="Virulence.csv", label="Virulence CSV (optional)")
    mlst = mo.ui.text(value="MLST.csv", label="MLST CSV (optional)")
    sero = mo.ui.text(value="Serotype.csv", label="Serotype CSV (optional)")

    phi_method = mo.ui.dropdown(
        options=["adaptive", "percentile", "fixed"],
        value="adaptive",
        label="Edge threshold method",
    )
    surrogate_mi = mo.ui.switch(value=False, label="Enable surrogate MI (slow)")
    mi_n = mo.ui.number(value=200, label="MI surrogate n", step=50)
    fdr = mo.ui.number(value=0.05, label="FDR alpha", step=0.01)

    config_strict = mo.ui.switch(value=False, label="Config strict mode")
    schema_version = mo.ui.text(value="1.1", label="schema_version")

    mo.md("## Config builder")
    mo.vstack(
        [
            data_root,
            out_root,
            mo.hstack([mic, amr, vir]),
            mo.hstack([mlst, sero]),
            mo.hstack([phi_method, surrogate_mi, mi_n, fdr]),
            mo.hstack([config_strict, schema_version]),
        ]
    )
    return (
        amr,
        config_strict,
        data_root,
        fdr,
        mic,
        mi_n,
        mlst,
        out_root,
        phi_method,
        schema_version,
        sero,
        surrogate_mi,
        vir,
    )


@app.cell
def _(
    Path,
    mo,
    yaml,
    amr,
    config_strict,
    data_root,
    fdr,
    mic,
    mi_n,
    mlst,
    out_root,
    phi_method,
    schema_version,
    sero,
    surrogate_mi,
    vir,
):
    def _abs(p: str) -> str:
        p = p.strip()
        if not p:
            return ""
        pp = Path(p)
        if pp.is_absolute():
            return str(pp)
        return str(Path(data_root.value) / pp)

    layers = []
    for name, widget, fmt in [
        ("MIC", mic, "wide"),
        ("AMR_genes", amr, "wide"),
        ("Virulence", vir, "wide"),
        ("MLST", mlst, "wide"),
        ("Serotype", sero, "wide"),
    ]:
        p = _abs(widget.value)
        if p:
            layers.append({"name": name, "path": p, "format": fmt})

    cfg = {
        "schema_version": schema_version.value,
        "config_strict": bool(config_strict.value),
        "output_dir": str(Path(out_root.value)),
        "align_mode": "union",
        "layers": layers,
        "association": {
            "phi_method": phi_method.value,
            "fdr_alpha": float(fdr.value),
            "surrogate_mi_enabled": bool(surrogate_mi.value),
            "mi_n_surrogate": int(mi_n.value),
        },
        "topology": {"enabled": True, "simplicial_dim": 2},
        "ising": {"enabled": True},
        "reliability": {"enabled": True, "fail_fast": False},
    }

    editor = mo.ui.text_area(
        value=yaml.safe_dump(cfg, sort_keys=False),
        label="Config YAML (editable)",
        rows=20,
    )
    mo.md("## Config editor")
    editor
    return editor


@app.cell
def _(Path, json, mo, subprocess, sys, time, editor, out_root):
    mo.md("## Run controls")
    start = mo.ui.button(label="Start analysis", kind="success")
    refresh = mo.ui.button(label="Refresh status")
    mo.hstack([start, refresh])

    out_dir = Path(out_root.value)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = out_dir / "config_used.yaml"
    log_path = out_dir / "run.log"
    state_path = out_dir / "run_state.json"

    if start.value:
        cfg_path.write_text(editor.value, encoding="utf-8")
        cmd = [sys.executable, "-m", "bactassocnet.cli", "--config", str(cfg_path)]
        with open(log_path, "ab") as f:
            proc = subprocess.Popen(
                cmd, stdout=f, stderr=subprocess.STDOUT, cwd=str(out_dir)
            )
        state_path.write_text(
            json.dumps(
                {
                    "pid": proc.pid,
                    "cmd": cmd,
                    "started_at": time.time(),
                    "cwd": str(out_dir),
                },
                indent=2,
            ),
            encoding="utf-8",
        )

    return cfg_path, log_path, refresh, state_path


@app.cell
def _(mo, log_path, refresh):
    if not log_path.exists():
        mo.md("No log yet.")
    else:
        with open(log_path, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 20000), 0)
            tail = f.read().decode("utf-8", errors="replace")
        mo.md("### Log (tail)")
        mo.code(tail)
    return


@app.cell
def _(Path, mo, pd, out_root):
    mo.md("## Results explorer")
    out_dir = Path(out_root.value)
    pattern = mo.ui.text(value="*.csv", label="File glob")
    limit = mo.ui.number(value=200, label="Preview rows", step=50)
    label_map = mo.ui.text_area(
        value="{}", label="Label map (JSON) — applied to displayed tables", rows=4
    )
    mo.vstack([mo.hstack([pattern, limit]), label_map])

    files = (
        sorted([p for p in out_dir.rglob(pattern.value) if p.is_file()])
        if out_dir.exists()
        else []
    )
    if not files:
        mo.md("No files yet.")
        return

    display = files[:200]
    selector = mo.ui.dropdown(
        options=[str(p.relative_to(out_dir)) for p in display],
        value=str(display[0].relative_to(out_dir)),
        label=f"Select a file ({len(display)}/{len(files)})",
    )
    file_path = out_dir / selector.value
    selector

    if file_path.suffix.lower() == ".csv":
        df = pd.read_csv(file_path, nrows=int(limit.value))
        try:
            mapping = __import__("json").loads(label_map.value or "{}")
            if isinstance(mapping, dict) and mapping:
                df = df.rename(
                    columns={k: v for k, v in mapping.items() if k in df.columns}
                )
                for c in df.columns:
                    if df[c].dtype == object:
                        df[c] = df[c].replace(mapping)
        except Exception:
            pass
        mo.ui.table(df, pagination=True, label=f"Preview: {selector.value}")
    else:
        mo.md(f"Selected: `{selector.value}`")
    return


@app.cell
def _(Path, mo, out_root, pd):
    mo.md("## Quick plots")
    out_dir = Path(out_root.value)
    edges_path = out_dir / "associations_filtered.csv"
    if not edges_path.exists():
        mo.md("`associations_filtered.csv` not found yet.")
        return
    try:
        edges = pd.read_csv(edges_path)
    except Exception as e:
        mo.md(f"Failed to read associations_filtered.csv: `{e}`")
        return

    # Try to infer source/target columns
    src = (
        "source"
        if "source" in edges.columns
        else ("Feature_A" if "Feature_A" in edges.columns else None)
    )
    dst = (
        "target"
        if "target" in edges.columns
        else ("Feature_B" if "Feature_B" in edges.columns else None)
    )
    if src is None or dst is None:
        mo.md(
            "Could not infer edge columns (expected source/target or Feature_A/Feature_B)."
        )
        return

    # Degree from edges (sample if huge)
    if len(edges) > 200000:
        edges = edges.sample(200000, random_state=0)
    deg = pd.concat([edges[src], edges[dst]]).value_counts()

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.hist(deg.values, bins=30)
    ax.set_xlabel("Degree")
    ax.set_ylabel("Count of nodes")
    ax.set_title("Degree distribution (from filtered edges)")
    mo.mpl.interactive(fig)
    return


@app.cell
def _(mo):
    mo.md("## Output reference")
    artifacts = [
        {
            "file": "associations_all.csv",
            "meaning": "All tested feature pairs + association statistics.",
            "interpretation": "Use to audit false positives; filter by prevalence and effect size.",
        },
        {
            "file": "associations_filtered.csv",
            "meaning": "Edges after thresholding and (optionally) FDR.",
            "interpretation": "Main edge list; verify robustness with sensitivity settings.",
        },
        {
            "file": "network_summary.csv",
            "meaning": "Network-level metrics.",
            "interpretation": "Check components/hubs; a single giant component may indicate a low threshold.",
        },
        {
            "file": "ising_direct_edges.csv",
            "meaning": "Approximate direct edges (reduces indirect correlations).",
            "interpretation": "Prefer for mechanistic hypotheses; still not causal.",
        },
        {
            "file": "simplicial_sweep.csv",
            "meaning": "Higher-order topology across thresholds.",
            "interpretation": "Betti/Euler changes indicate modularity vs dense cores.",
        },
        {
            "file": "topology_persistence.csv",
            "meaning": "Topological persistence vs threshold.",
            "interpretation": "Features stable across thresholds are more robust.",
        },
        {
            "file": "multilayer_dispersion.csv",
            "meaning": "Layer dispersion/informativeness metrics.",
            "interpretation": "Helps identify layers dominating the graph.",
        },
        {
            "file": "layer_flow.csv",
            "meaning": "Cross-layer information flow summary.",
            "interpretation": "Shows which layer drives associations.",
        },
        {
            "file": "config_validation.json",
            "meaning": "Config schema + unknown-key report.",
            "interpretation": "Fix config if strict mode fails.",
        },
        {
            "file": "run_manifest.json",
            "meaning": "Input hashes + environment.",
            "interpretation": "Reproducibility.",
        },
    ]
    mo.ui.table(artifacts, pagination=True, label="Artifacts")

    mo.md(
        """
        ## Statistical assumptions (practical)
        - Association ≠ causation; edges capture dependence and may be driven by clonal structure.
        - Sparse features inflate MI/phi variance; use prevalence filtering.
        - Surrogate/permutation MI controls type-I error but is computationally expensive.
        """
    )
    return


if __name__ == "__main__":
    app.run()
