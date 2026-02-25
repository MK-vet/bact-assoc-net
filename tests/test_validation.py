"""Signal-recovery and integration validation for bact-assoc-net."""

import subprocess
import sys
import json
import numpy as np
import pandas as pd
import pytest


# ── fixtures ──────────────────────────────────────────────────────────────


@pytest.fixture
def correlated_pair():
    """Strongly correlated binary pair (P(Y=1|X=1)=0.90)."""
    rng = np.random.RandomState(5)
    x = rng.binomial(1, 0.5, 500)
    y = np.where(x == 1, rng.binomial(1, 0.90, 500), rng.binomial(1, 0.05, 500))
    return pd.Series(x), pd.Series(y)


@pytest.fixture
def independent_pair():
    rng = np.random.RandomState(6)
    return pd.Series(rng.binomial(1, 0.5, 500)), pd.Series(rng.binomial(1, 0.5, 500))


@pytest.fixture
def planted_df():
    """300 samples. gene_0 → pheno_0 (P=0.90). All others independent."""
    rng = np.random.RandomState(3)
    n = 300
    g0 = rng.binomial(1, 0.5, n)
    p0 = np.where(g0 == 1, rng.binomial(1, 0.90, n), rng.binomial(1, 0.05, n))
    return pd.DataFrame(
        {
            "gene_0": g0,
            "gene_1": rng.binomial(1, 0.4, n),
            "gene_2": rng.binomial(1, 0.3, n),
            "pheno_0": p0,
            "pheno_1": rng.binomial(1, 0.35, n),
        }
    )


# ── signal recovery ───────────────────────────────────────────────────────


class TestSignalRecovery:
    def test_mi_correlated_exceeds_independent(self, correlated_pair, independent_pair):
        """MI(correlated) substantially > MI(independent)."""
        from bactassocnet.information_theory.core import mi

        assert mi(*correlated_pair) > mi(*independent_pair) + 0.1

    def test_nmi_self_equals_one(self):
        """NMI(X, X) = 1.0."""
        from bactassocnet.information_theory.core import nmi

        x = np.array([0, 0, 1, 1, 1, 0])
        assert abs(nmi(x, x) - 1.0) < 1e-9

    def test_nmi_independent_near_zero(self, independent_pair):
        """NMI of independent binary variables ≈ 0 (< 0.05)."""
        from bactassocnet.information_theory.core import nmi

        assert nmi(independent_pair[0].values, independent_pair[1].values) < 0.05

    def test_surrogate_mi_significant_for_correlated(self, correlated_pair):
        """Surrogate MI z-score > 3 for strongly correlated pair."""
        from bactassocnet.information_theory.advanced import surrogate_mi

        # Result keys: MI, Z, P, N, Mean_Surr, SD_Surr
        r = surrogate_mi(*correlated_pair, n_surrogates=199, seed=0)
        assert r["Z"] > 3.0, (
            f"Expected Z > 3 for strongly correlated pair, got Z={r['Z']:.2f}"
        )

    def test_surrogate_mi_nonsignificant_for_independent(self, independent_pair):
        """Surrogate MI z-score < 2 for independent pair (false-positive control)."""
        from bactassocnet.information_theory.advanced import surrogate_mi

        r = surrogate_mi(*independent_pair, n_surrogates=199, seed=0)
        assert r["Z"] < 2.0, f"False positive: Z={r['Z']:.2f} for independent pair"

    def test_pid_gene0_unique_exceeds_gene1(self, planted_df):
        """gene_0 has higher Unique information for pheno_0 than gene_1."""
        from bactassocnet.information_theory.advanced import pid

        r = pid(
            planted_df["pheno_0"].values,
            planted_df["gene_0"].values,
            planted_df["gene_1"].values,
        )
        assert r["Unique_S1"] > r["Unique_S2"] + 0.02, (
            f"Unique(gene_0)={r['Unique_S1']:.4f} not > Unique(gene_1)={r['Unique_S2']:.4f}"
        )

    def test_pid_decomposition_sums_to_total(self, planted_df):
        """PID components sum to I_Total within floating-point tolerance."""
        from bactassocnet.information_theory.advanced import pid

        r = pid(
            planted_df["pheno_0"].values,
            planted_df["gene_0"].values,
            planted_df["gene_1"].values,
        )
        total = r["Redundancy"] + r["Unique_S1"] + r["Unique_S2"] + r["Synergy"]
        assert abs(total - r["I_Total"]) < 0.01

    def test_euler_characteristic(self):
        """Simplicial complex: χ = V - E + F for triangle graph."""
        from bactassocnet.information_theory.advanced import build_simplicial_complex

        edf = pd.DataFrame(
            {
                "Feature1": ["A", "A", "B"],
                "Feature2": ["B", "C", "C"],
                "phi": [0.8, 0.8, 0.8],
            }
        )
        sc = build_simplicial_complex(edf, "phi", threshold=0.5, max_dim=2)
        v = sc["counts"].get(0, 0)
        e = sc["counts"].get(1, 0)
        f = sc["counts"].get(2, 0)
        assert sc["euler"] == v - e + f


# ── reproducibility ───────────────────────────────────────────────────────


class TestReproducibility:
    def test_surrogate_mi_reproducible(self, correlated_pair):
        from bactassocnet.information_theory.advanced import surrogate_mi

        r1 = surrogate_mi(*correlated_pair, n_surrogates=50, seed=42)
        r2 = surrogate_mi(*correlated_pair, n_surrogates=50, seed=42)
        assert r1["MI"] == r2["MI"] and r1["Z"] == r2["Z"]


# ── CLI integration ───────────────────────────────────────────────────────


class TestCLI:
    def test_help_exits_zero(self):
        r = subprocess.run(
            [sys.executable, "-m", "bactassocnet.cli", "--help"], capture_output=True
        )
        assert r.returncode == 0

    def test_version_output(self):
        from bactassocnet import __version__

        r = subprocess.run(
            [sys.executable, "-m", "bactassocnet.cli", "--version"],
            capture_output=True,
            text=True,
        )
        assert __version__ in r.stdout

    def test_self_check_passes(self):
        r = subprocess.run(
            [sys.executable, "-m", "bactassocnet.cli", "--self-check"],
            capture_output=True,
            text=True,
        )
        assert r.returncode == 0
        assert json.loads(r.stdout)["status"] == "PASS"

    def test_run_produces_mandatory_outputs(self, planted_df, tmp_path):
        import yaml

        data_dir = tmp_path / "data"
        data_dir.mkdir()
        out_dir = tmp_path / "results"
        planted_df["Strain_ID"] = [f"S{i:04d}" for i in range(len(planted_df))]
        planted_df.to_csv(data_dir / "Layer.csv", index=False)

        config = {
            "schema_version": "1.1",
            "layers": [
                {
                    "name": "Layer",
                    "path": str(data_dir / "Layer.csv"),
                    "id_column": "Strain_ID",
                }
            ],
            "output_dir": str(out_dir),
            "association": {
                "method": "phi",
                "n_surrogates": 50,
                "threshold_method": "percentile",
                "threshold_value": 75,
            },
            "n_jobs": 1,
            "seed": 0,
        }
        with open(tmp_path / "config.yaml", "w") as f:
            yaml.dump(config, f)

        r = subprocess.run(
            [
                sys.executable,
                "-m",
                "bactassocnet.cli",
                "--config",
                str(tmp_path / "config.yaml"),
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert r.returncode == 0, f"CLI failed:\n{r.stderr}"
        assert (out_dir / "associations_all.csv").exists()
        assert (out_dir / "run_manifest.json").exists()
