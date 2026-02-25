"""Tests for bact-assoc-net information-theoretic and topological modules."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def binary_pair():
    """Two correlated binary features."""
    rng = np.random.RandomState(42)
    x = rng.binomial(1, 0.5, 200)
    y = x.copy()
    flip = rng.choice(200, 30, replace=False)
    y[flip] = 1 - y[flip]
    return x, y


@pytest.fixture
def binary_df():
    """100 samples × 8 features across 2 layers."""
    rng = np.random.RandomState(42)
    n = 100
    df = pd.DataFrame({f"gene_{i}": rng.binomial(1, 0.4, n) for i in range(4)})
    for i in range(4):
        df[f"pheno_{i}"] = rng.binomial(1, 0.3, n)
    # Plant signal: gene_0 → pheno_0
    df.loc[df["gene_0"] == 1, "pheno_0"] = rng.binomial(
        1, 0.9, (df["gene_0"] == 1).sum()
    )
    return df


class TestCoreEntropy:
    def test_shannon_nonneg(self, binary_pair):
        from bactassocnet.information_theory.core import shannon_H

        assert shannon_H(binary_pair[0]) >= 0

    def test_mi_nonneg(self, binary_pair):
        from bactassocnet.information_theory.core import mi

        assert mi(*binary_pair) >= 0

    def test_nmi_range(self, binary_pair):
        from bactassocnet.information_theory.core import nmi

        v = nmi(*binary_pair)
        assert 0 <= v <= 1.01  # allow float precision

    def test_chi2_phi_shape(self, binary_pair):
        from bactassocnet.information_theory.core import chi2_phi

        r = chi2_phi(pd.Series(binary_pair[0]), pd.Series(binary_pair[1]))
        assert "phi" in r and "p" in r

    def test_adaptive_threshold(self):
        from bactassocnet.information_theory.core import adaptive_phi_threshold

        phis = np.array([0.1, 0.2, 0.3, 0.5, 0.8])
        t = adaptive_phi_threshold(phis, "percentile", 90)
        assert 0 < t < 1


class TestSurrogateMI:
    def test_smi_positive_for_correlated(self, binary_pair):
        from bactassocnet.information_theory.advanced import surrogate_mi

        r = surrogate_mi(
            pd.Series(binary_pair[0]),
            pd.Series(binary_pair[1]),
            n_surrogates=50,
            seed=42,
        )
        assert r["MI"] >= 0

    def test_smi_symmetry(self, binary_df):
        from bactassocnet.information_theory.advanced import surrogate_mi

        x = binary_df["gene_0"]
        y = binary_df["pheno_0"]
        a = surrogate_mi(x, y, n_surrogates=20, seed=1)["MI"]
        b = surrogate_mi(y, x, n_surrogates=20, seed=1)["MI"]
        assert abs(a - b) < 1e-12


class TestPID:
    def test_pid_decomposition(self, binary_df):
        from bactassocnet.information_theory.advanced import pid

        t = binary_df["pheno_0"].values
        s1 = binary_df["gene_0"].values
        s2 = binary_df["gene_1"].values
        r = pid(t, s1, s2)
        # Decomposition should sum to I_Total (approximately)
        total = r["Redundancy"] + r["Unique_S1"] + r["Unique_S2"] + r["Synergy"]
        assert abs(total - r["I_Total"]) < 0.01

    def test_pid_dominant_field(self, binary_df):
        from bactassocnet.information_theory.advanced import pid

        r = pid(
            binary_df["pheno_0"].values,
            binary_df["gene_0"].values,
            binary_df["gene_1"].values,
        )
        assert r["Dominant"] in ("redundant", "synergistic")


class TestSimplicial:
    def test_build_simplicial(self):
        from bactassocnet.information_theory.advanced import build_simplicial_complex

        edf = pd.DataFrame(
            {
                "Feature1": ["A", "A", "B"],
                "Feature2": ["B", "C", "C"],
                "phi": [0.5, 0.4, 0.6],
            }
        )
        sc = build_simplicial_complex(edf, "phi", threshold=0.3, max_dim=2)
        assert 0 in sc["counts"]  # nodes
        assert 1 in sc["counts"]  # edges
        assert "euler" in sc


class TestConfig:
    def test_resolve_layers_from_dir(self, tmp_path):
        from bactassocnet.config import Config

        # Create dummy CSVs
        for name in ["AMR", "VIR"]:
            pd.DataFrame({"Strain_ID": ["S1"], "f": [1]}).to_csv(
                tmp_path / f"{name}.csv", index=False
            )
        cfg = Config(input_dir=str(tmp_path))
        layers = cfg.resolve_layers()
        assert len(layers) == 2

    def test_roundtrip_yaml(self, tmp_path):
        from bactassocnet.config import Config

        cfg = Config(output_dir="/tmp/out")
        p = tmp_path / "c.yaml"
        cfg.to_yaml(p)
        cfg2 = Config.from_yaml(p)
        assert cfg2.output_dir == "/tmp/out"
