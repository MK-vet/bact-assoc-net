"""Advanced information-theoretic and topological utilities.

Methodological note:
- For cross-sectional data, "transfer entropy" is not identifiable without a temporal ordering + lag.
  The legacy implementation effectively reduced to MI. We therefore expose surrogate-tested MI
  (permutation surrogates) as the defensible default for cross-sectional layers.

- A true TE is exposed only for genuine time series with an explicit lag and ordered observations.
"""
from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd

from .core import mutual_information


def surrogate_mi(
    x: pd.Series,
    y: pd.Series,
    n_surrogates: int = 999,
    seed: int = 42,
) -> Dict[str, float]:
    """Mutual information with a permutation-based surrogate test.

    Missing values are excluded (complete-case).
    P-value is one-sided: P(MI_surrogate >= MI_observed).

    n_surrogates
    ------------
    The default of 999 gives a minimum achievable p-value of 1/1000 = 0.001,
    which is sufficient for FDR correction at alpha=0.05 with up to ~100 pairs.
    When the caller passes n_surrogates explicitly, that value is used as-is.
    For large panels (many pairs tested simultaneously), consider setting
    n_surrogates=9999 to achieve p_min = 0.0001.
    """
    x = pd.Series(x)
    y = pd.Series(y)
    m = x.notna() & y.notna()
    x = x[m].astype(int).to_numpy()
    y = y[m].astype(int).to_numpy()
    n = int(len(x))
    if n < 10:
        return {"MI": np.nan, "Z": np.nan, "P": 1.0, "N": n, "Mean_Surr": np.nan, "SD_Surr": np.nan}

    mi0 = float(mutual_information(x, y))

    rng = np.random.RandomState(seed)
    surr = np.empty(n_surrogates, dtype=float)
    for i in range(n_surrogates):
        yp = y.copy()
        rng.shuffle(yp)
        surr[i] = mutual_information(x, yp)

    mean = float(np.mean(surr))
    sd = float(np.std(surr, ddof=1)) if n_surrogates > 1 else np.nan
    z = (mi0 - mean) / (sd + 1e-12) if np.isfinite(sd) else np.nan
    p = float((np.sum(surr >= mi0) + 1) / (n_surrogates + 1))

    return {"MI": mi0, "Z": float(z), "P": p, "N": n, "Mean_Surr": mean, "SD_Surr": sd}


def transfer_entropy(
    source: np.ndarray,
    target: np.ndarray,
    lag: int = 1,
    time_series: bool = False,
) -> Dict[str, float]:
    """Transfer entropy wrapper.

    For cross-sectional layers, TE is not defined: set time_series=True to enable.
    """
    if not time_series:
        raise ValueError("Transfer entropy requires time_series=True and an explicit lag on ordered observations.")
    te = transfer_entropy_timeseries(source, target, lag=lag)
    return {"TE": float(te), "Lag": int(lag)}


def transfer_entropy_timeseries(
    source: np.ndarray,
    target: np.ndarray,
    lag: int = 1,
) -> float:
    """True transfer entropy for time series (binary), TE = I(S_{t-lag}; T_t | T_{t-lag})."""
    source = np.asarray(source).astype(int)
    target = np.asarray(target).astype(int)
    if lag < 1 or len(source) != len(target) or len(source) <= lag + 5:
        raise ValueError("Invalid lag or series length.")

    s_lag = source[:-lag]
    t_now = target[lag:]
    t_lag = target[:-lag]

    te = 0.0
    for u in (0, 1):
        m = (t_lag == u)
        if m.sum() < 5:
            continue
        te += (m.mean()) * mutual_information(s_lag[m], t_now[m])
    return float(te)


def pid(t: np.ndarray, s1: np.ndarray, s2: np.ndarray) -> Dict[str, float]:
    """Partial Information Decomposition using I_min redundancy (Williams & Beer 2010).

    Decomposes the total mutual information I(T; S1, S2) into four non-negative
    components: Redundancy (shared by both sources), Unique_S1 (unique to S1),
    Unique_S2 (unique to S2), and Synergy (only available jointly).

    Redundancy measure
    ------------------
    We use I_min = min(I(T;S1), I(T;S2)) as the redundancy measure.  This is
    the original Williams & Beer (2010) proposal and satisfies the consistency
    axiom: R ≥ 0 and Unique ≥ 0 are guaranteed.

    Alternative redundancy measures exist — notably Bertschinger et al. (2014)
    BROJA-PID which is operationally motivated and Makkeh et al. (2021) PROJ
    which is information-geometric — and may yield different decompositions,
    particularly when S1 and S2 are correlated.  I_min is appropriate for
    exploratory analysis of binary AMR gene-phenotype pairs where computational
    efficiency and non-negativity are priorities.

    References
    ----------
    Williams & Beer (2010) arXiv:1004.2515.
    Bertschinger et al. (2014) Entropy 16:2161.
    Makkeh et al. (2021) Entropy 23:141.
    """
    t = np.asarray(t).astype(int)
    s1 = np.asarray(s1).astype(int)
    s2 = np.asarray(s2).astype(int)

    I1 = float(mutual_information(t, s1))
    I2 = float(mutual_information(t, s2))
    # joint source as 4-state variable
    joint = (s1.astype(int) << 1) | s2.astype(int)
    Itot = float(mutual_information(t, joint))

    red = min(I1, I2)
    u1 = max(0.0, I1 - red)
    u2 = max(0.0, I2 - red)
    syn = Itot - red - u1 - u2

    dominant = "redundant" if red >= syn else "synergistic"
    return {
        "I_S1": I1, "I_S2": I2, "I_Total": Itot,
        "Redundancy": float(red),
        "Unique_S1": float(u1),
        "Unique_S2": float(u2),
        "Synergy": float(syn),
        "Dominant": dominant,
    }


def build_simplicial_complex(
    edges: pd.DataFrame,
    weight_col: str = "phi",
    threshold: float = 0.2,
    max_dim: int = 2,
) -> Dict[str, object]:
    """Build a simplicial summary from an edge list by thresholding |weight|.

    Returns counts by simplex dimension and Euler characteristic (using triangles for dim=2).
    """
    import networkx as nx
    edf = edges.copy()
    if "Feature1" in edf.columns and "Feature2" in edf.columns:
        ucol, vcol = "Feature1", "Feature2"
    else:
        ucol, vcol = "Feature_A", "Feature_B"
    edf = edf[edf[weight_col].abs() >= threshold]
    G = nx.Graph()
    for _, r in edf.iterrows():
        G.add_edge(str(r[ucol]), str(r[vcol]), weight=float(r[weight_col]))

    counts = {0: G.number_of_nodes(), 1: G.number_of_edges()}
    euler = counts[0] - counts[1]
    if max_dim >= 2:
        T = sum(nx.triangles(G).values()) // 3
        counts[2] = int(T)
        euler = euler + int(T)

    return {"counts": counts, "euler": int(euler)}
