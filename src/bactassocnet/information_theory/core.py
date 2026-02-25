"""
Core information-theoretic and association metrics.

Retained from v1: entropy, conditional entropy, NMI, chi2_phi, adaptive threshold.
These are the building blocks; advanced.py adds the novel contributions.
"""
from __future__ import annotations
from itertools import combinations
from typing import Dict, List
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, entropy as sp_entropy

# ── entropy functions ────────────────────────────────────────────────────

def shannon_H(x: np.ndarray) -> float:
    _, c = np.unique(x, return_counts=True)
    return float(sp_entropy(c, base=np.e))

def joint_H(x: np.ndarray, y: np.ndarray) -> float:
    combo = np.column_stack([x, y])
    _, c = np.unique(combo, axis=0, return_counts=True)
    return float(sp_entropy(c, base=np.e))

def cond_H(x: np.ndarray, y: np.ndarray) -> float:
    """H(X|Y) = H(X,Y) − H(Y)."""
    return joint_H(x, y) - shannon_H(y)

def mi(x: np.ndarray, y: np.ndarray) -> float:
    """I(X;Y) = H(X) + H(Y) − H(X,Y)."""
    return shannon_H(x) + shannon_H(y) - joint_H(x, y)

# Backwards-compatible alias used by the pipeline.
def mutual_information(x: np.ndarray, y: np.ndarray) -> float:
    """Alias for :func:`mi` (kept for API stability)."""
    return mi(x, y)

def nmi(x: np.ndarray, y: np.ndarray) -> float:
    """Normalised MI = I(X;Y) / sqrt(H(X)·H(Y))."""
    hx, hy = shannon_H(x), shannon_H(y)
    return mi(x, y) / np.sqrt(hx * hy) if hx > 0 and hy > 0 else 0.0

# ── chi-square / phi ─────────────────────────────────────────────────────

def chi2_phi(x: pd.Series, y: pd.Series) -> Dict[str, float]:
    """Chi-square/Fisher association for binary data with explicit missingness handling.

    Important: missing values are excluded (complete-case) and are NOT imputed to 0.
    """
    x = pd.Series(x)
    y = pd.Series(y)
    m = x.notna() & y.notna()
    if int(m.sum()) < 4:
        return {"chi2": np.nan, "p": 1.0, "phi": np.nan, "phi_ci_lo": np.nan, "phi_ci_hi": np.nan, "n": int(m.sum())}

    xv = x[m].astype(int).to_numpy()
    yv = y[m].astype(int).to_numpy()

    a = int(np.sum((xv == 0) & (yv == 0)))
    b = int(np.sum((xv == 0) & (yv == 1)))
    c = int(np.sum((xv == 1) & (yv == 0)))
    d = int(np.sum((xv == 1) & (yv == 1)))
    n = a + b + c + d

    denom = max((a + b) * (c + d) * (a + c) * (b + d), 1)
    phi = (a * d - b * c) / np.sqrt(denom)

    tab = np.array([[a, b], [c, d]], dtype=int)
    exp = np.outer(tab.sum(1), tab.sum(0)) / max(n, 1)
    if (exp < 5).any():
        _, p = fisher_exact(tab)
        chi2 = np.nan
    else:
        chi2, p, _, _ = chi2_contingency(tab, correction=False)

    # Fisher-z CI (approx.)
    z = np.arctanh(np.clip(phi, -0.999, 0.999))
    se = 1 / np.sqrt(max(n - 3, 1))
    lo = np.tanh(z - 1.96 * se)
    hi = np.tanh(z + 1.96 * se)

    return {"chi2": float(chi2) if np.isfinite(chi2) else np.nan,
            "p": float(p),
            "phi": float(phi),
            "phi_ci_lo": float(lo),
            "phi_ci_hi": float(hi),
            "n": int(n)}

# ── adaptive threshold ───────────────────────────────────────────────────

def adaptive_phi_threshold(
    phis: np.ndarray, method: str = "percentile", percentile: int = 90,
) -> float:
    """Adaptive edge threshold: percentile, IQR, or statistical."""
    phis = np.abs(phis[np.isfinite(phis)])
    if len(phis) == 0: return 0.3
    if method == "percentile":
        return float(np.percentile(phis, percentile))
    elif method == "iqr":
        q1, q3 = np.percentile(phis, [25, 75])
        return float(q3 + 1.5 * (q3 - q1))
    elif method == "statistical":
        return float(np.mean(phis) + 2 * np.std(phis))
    return 0.3

# ── mutual exclusivity ──────────────────────────────────────────────────

def find_mutually_exclusive(
    df: pd.DataFrame, features: List[str], k: int = 2,
) -> pd.DataFrame:
    """Detect feature pairs/triplets that never co-occur."""
    rows = []
    for sub in combinations(features, k):
        cols = df[list(sub)]
        co = (cols.sum(axis=1) == k).sum()
        if co == 0:
            rows.append({"Features": ", ".join(sub), "Order": k, "Co_occurrences": 0})
    return pd.DataFrame(rows)
