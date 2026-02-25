from __future__ import annotations
import numpy as np
import pandas as pd

def _wide_from_edges(edges_df: pd.DataFrame):
    feats=sorted(set(edges_df['Feature_A']).union(edges_df['Feature_B']))
    S=np.eye(len(feats), dtype=float)
    idx={f:i for i,f in enumerate(feats)}
    wcol = 'phi' if 'phi' in edges_df.columns else ('MI' if 'MI' in edges_df.columns else None)
    if wcol is None:
        return feats, S
    for _,r in edges_df.iterrows():
        i,j=idx[r['Feature_A']],idx[r['Feature_B']]
        try: w=float(r.get(wcol,0.0))
        except Exception: w=0.0
        if np.isfinite(w): S[i,j]=S[j,i]=w
    return feats, (S+S.T)/2

def _plm_ising(binary_df: pd.DataFrame, alpha: float = 0.01):
    """Pseudolikelihood Method with L1 regularisation for Ising model estimation.

    Fits one L1-penalised logistic regression per feature (y = feature_j,
    predictors = all other features) and symmetrises the coupling matrix:
    W = 0.5 * (B + B^T).
    """
    try:
        from sklearn.linear_model import LogisticRegression
        import sklearn
        _sklearn_version = tuple(int(x) for x in sklearn.__version__.split('.')[:2])
    except Exception:
        return pd.DataFrame(columns=['Feature_A', 'Feature_B', 'IsingDirect', 'Method'])

    X = binary_df.copy()
    if X.empty or X.shape[1] < 2:
        return pd.DataFrame(columns=['Feature_A', 'Feature_B', 'IsingDirect', 'Method'])
    X = X.dropna(axis=0, how='any')
    if X.shape[0] < 20:
        X = binary_df.copy()
        for c in X.columns:
            s = X[c].dropna()
            X[c] = X[c].fillna(int(round(float(s.mean()))) if not s.empty else 0)
    X = X.astype(int)
    feats = X.columns.tolist()
    n, p = X.shape
    if p < 2:
        return pd.DataFrame(columns=['Feature_A', 'Feature_B', 'IsingDirect', 'Method'])

    B = np.zeros((p, p), float)
    Cinv = max(alpha, 1e-4)

    for j in range(p):
        y = X.iloc[:, j].to_numpy()
        if np.unique(y).size < 2:
            continue
        Xj = X.drop(columns=[feats[j]]).to_numpy()
        try:
            # Use l1_ratio=1.0 with elasticnet for sklearn >= 1.2 (no penalty= deprecation)
            if _sklearn_version >= (1, 2):
                lr = LogisticRegression(
                    penalty='elasticnet',
                    l1_ratio=1.0,
                    solver='saga',
                    C=1.0 / Cinv,
                    max_iter=1000,
                )
            else:
                lr = LogisticRegression(
                    penalty='l1',
                    solver='liblinear',
                    C=1.0 / Cinv,
                    max_iter=1000,
                )
            lr.fit(Xj, y)
            coef = lr.coef_.ravel()
        except Exception:
            # Final fallback: liblinear L2 (weaker but stable)
            try:
                lr = LogisticRegression(solver='liblinear', C=1.0 / Cinv, max_iter=1000)
                lr.fit(Xj, y)
                coef = lr.coef_.ravel()
            except Exception:
                continue
        k = 0
        for i in range(p):
            if i == j:
                continue
            B[i, j] = coef[k]
            k += 1

    W = 0.5 * (B + B.T)
    rows = []
    for i in range(p):
        for j in range(i + 1, p):
            w = float(W[i, j])
            if np.isfinite(w) and w != 0:
                rows.append({'Feature_A': feats[i], 'Feature_B': feats[j],
                             'IsingDirect': w, 'Method': 'PLM-L1'})
    return pd.DataFrame(rows)

def glasso_ising_proxy(edges_df: pd.DataFrame, feature_matrix: pd.DataFrame | None = None, thr: float = 0.0, alpha: float = 0.01) -> pd.DataFrame:
    """Estimate direct (conditional) interactions between binary AMR features.

    Method: Pseudolikelihood Method with L1 regularisation (PLM-L1).
    Per feature j, fit L1-penalised logistic regression with all other features
    as predictors; symmetrise: W = 0.5*(B + B^T).

    Note on Graphical LASSO
    -----------------------
    Graphical LASSO assumes multivariate Gaussian data and is therefore
    methodologically inappropriate for binary (Bernoulli) resistance data —
    normalising binary indicators to N(0,1) does not satisfy the Gaussian
    copula requirement. PLM-L1 is the statistically valid alternative for
    binary data and is the sole method used here.

    Parameters
    ----------
    edges_df : pd.DataFrame
        Edge list with columns Feature_A, Feature_B (used only when
        feature_matrix is None and a fallback empty result is needed).
    feature_matrix : pd.DataFrame or None
        Binary feature matrix (rows = isolates, columns = AMR features).
        When provided, PLM-L1 is fitted directly on this matrix.
    thr : float
        Minimum absolute IsingDirect weight to retain an edge (default 0 = keep all).
    alpha : float
        L1 regularisation strength for logistic regression (C = 1/alpha).

    Returns
    -------
    pd.DataFrame with columns Feature_A, Feature_B, IsingDirect, Method.
    """
    if feature_matrix is not None:
        X = feature_matrix.select_dtypes(include=[np.number]).copy()
        keep = [c for c in X.columns if X[c].notna().sum() >= 10 and X[c].dropna().nunique() > 1]
        X = X[keep]
        if X.shape[1] >= 2:
            plm = _plm_ising(X, alpha=alpha)
            if not plm.empty:
                if thr > 0:
                    plm = plm[plm['IsingDirect'].abs() >= thr].copy()
                return plm.sort_values('IsingDirect', key=lambda s: s.abs(), ascending=False).reset_index(drop=True)
    # Fallback: no feature_matrix or PLM returned empty — return empty frame
    return pd.DataFrame(columns=['Feature_A', 'Feature_B', 'IsingDirect', 'Method'])

def topology_persistence_from_sweep(sweep_df: pd.DataFrame) -> pd.DataFrame:
    if sweep_df.empty: return pd.DataFrame()
    s=sweep_df.sort_values('phi_threshold').copy()
    for col in ['V','E','Triangles','Euler']:
        if col in s.columns: s[f'd{col}']=s[col].diff()
    if {'E','Triangles','V'}.issubset(set(s.columns)):
        s['TriangleDensity']=s['Triangles']/s['E'].clip(lower=1)
        s['EdgePerNode']=s['E']/s['V'].clip(lower=1)
    return s

def multilayer_dispersion(edges_df: pd.DataFrame) -> pd.DataFrame:
    if edges_df.empty: return pd.DataFrame()
    def split_layer(f): return str(f).split('::',1)[0]
    rows=[]
    all_layers=sorted(set(split_layer(x) for x in pd.concat([edges_df['Feature_A'], edges_df['Feature_B']],axis=0)))
    for f in sorted(set(edges_df['Feature_A']).union(edges_df['Feature_B'])):
        sub = edges_df[(edges_df['Feature_A']==f)|(edges_df['Feature_B']==f)]
        if sub.empty: continue
        cnt={L:0.0 for L in all_layers}
        for _,r in sub.iterrows():
            other = r['Feature_B'] if r['Feature_A']==f else r['Feature_A']
            w = abs(float(r.get('MI',0.0))) + abs(float(r.get('phi',0.0)))
            if not np.isfinite(w): w=0.0
            cnt[split_layer(other)] += w
        p=np.array([cnt[L] for L in all_layers],float)
        if p.sum()==0: continue
        p/=p.sum(); m=np.ones_like(p)/len(p)
        jsd=0.5*np.sum(p*np.log2((p+1e-12)/(0.5*(p+m)+1e-12)))+0.5*np.sum(m*np.log2((m+1e-12)/(0.5*(p+m)+1e-12)))
        part=1.0-np.sum(p**2)
        rows.append({'Feature':f,'JSD_to_uniform_layers':float(jsd),'ParticipationCoeff':float(part),'Degree':int(len(sub))})
    return pd.DataFrame(rows).sort_values(['ParticipationCoeff','Degree'], ascending=[False,False]) if rows else pd.DataFrame()
