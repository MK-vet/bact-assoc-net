
from __future__ import annotations
import json
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
def _assoc(df, max_pairs=2000):
    a=df.to_numpy(dtype=np.uint8); n=a.shape[1]; cnt=0; checksum=0.0
    for i in range(min(n,2000)):
        xi=a[:,i]
        for j in range(i+1,min(n,2000)):
            xj=a[:,j]
            q=((xi==1)&(xj==1)).sum(); w=((xi==1)&(xj==0)).sum(); e=((xi==0)&(xj==1)).sum(); r=((xi==0)&(xj==0)).sum()
            den=((q+w)*(e+r)*(q+e)*(w+r))**0.5
            if den: checksum += (q*r-w*e)/den
            cnt += 1
            if cnt>=max_pairs: return cnt, float(checksum)
    return cnt, float(checksum)
def run_benchmark(out_path=None, n_rows=300, p_values=(1000,5000,20000), seed=123):
    rng=np.random.default_rng(seed); rows=[]
    for p in p_values:
        X=(rng.random((n_rows,p))<0.08).astype(np.uint8)
        df=pd.DataFrame(X)
        t=time.perf_counter(); prev=df.mean(0); informative=int(((prev>0.01)&(prev<0.99)).sum()); qc=time.perf_counter()-t
        t=time.perf_counter(); pairs, chk=_assoc(df); assoc=time.perf_counter()-t
        rows.append({'n_rows':n_rows,'p':int(p),'informative_features':informative,'qc_sec':round(qc,4),
                      'assoc2000pairs_sec':round(assoc,4),'pairs_eval':pairs,'checksum':round(chk,6)})
    out={'tool':'bact-assoc-net','results':rows}
    if out_path: Path(out_path).write_text(json.dumps(out,indent=2))
    return out


def run_benchmark_paper_ready(out_dir, n_rows=300, p_values=(1000,5000,20000), seed=123):
    from pathlib import Path
    import json
    import pandas as pd
    import matplotlib.pyplot as plt
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    js = run_benchmark(None, n_rows=n_rows, p_values=p_values, seed=seed)
    (out_dir / "benchmark_synthetic.json").write_text(json.dumps(js, indent=2))
    df = pd.DataFrame(js.get("results", []))
    if not df.empty:
        df.insert(0, "tool", js.get("tool", "unknown"))
        df.to_csv(out_dir / "benchmark_synthetic.csv", index=False)
        fig = plt.figure(figsize=(6,4))
        plt.plot(df["p"], df["qc_sec"], marker="o", label="QC")
        plt.plot(df["p"], df["assoc2000pairs_sec"], marker="o", label="Assoc (2000 pairs)")
        plt.xlabel("Number of features (p)")
        plt.ylabel("Time (s)")
        plt.xscale("log")
        plt.legend()
        plt.tight_layout()
        fig.savefig(out_dir / "benchmark_runtime_vs_p.png", dpi=160)
        plt.close(fig)
        fig = plt.figure(figsize=(6,4))
        plt.plot(df["p"], df["informative_features"], marker="o")
        plt.xlabel("Number of features (p)")
        plt.ylabel("Informative features")
        plt.xscale("log")
        plt.tight_layout()
        fig.savefig(out_dir / "benchmark_informative_vs_p.png", dpi=160)
        plt.close(fig)
    try:
        run_benchmark_compare_models(out_dir, n_rows=n_rows, p_values=p_values, seed=seed)
    except Exception:
        pass
    return js

def _edge_set_from_df(df, top_k=100):
    if df is None or getattr(df, "empty", True):
        return set()
    d = df.copy()
    if "IsingDirect" in d.columns:
        col = "IsingDirect"
    elif "weight" in d.columns:
        col = "weight"
    else:
        return set()
    try:
        d = d.assign(_abs=d[col].astype(float).abs()).sort_values("_abs", ascending=False)
    except Exception:
        d = d.sort_values(col, ascending=False)
    d = d.head(int(top_k))
    return {tuple(sorted((str(a), str(b)))) for a, b in zip(d["Feature_A"], d["Feature_B"])}

def _glasso_edges_from_binary(Xsub, alpha=0.02):
    import numpy as np
    import pandas as pd
    try:
        from sklearn.covariance import GraphicalLasso
    except Exception:
        return pd.DataFrame(columns=["Feature_A","Feature_B","IsingDirect"])
    X2 = Xsub.astype(float).copy()
    X2 = (X2 - X2.mean(0)) / (X2.std(0, ddof=1) + 1e-9)
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    gl = GraphicalLasso(alpha=float(alpha), max_iter=100).fit(X2.to_numpy())
    Prec = gl.precision_
    feats = Xsub.columns.tolist()
    rows=[]
    for i in range(len(feats)):
        for j in range(i+1, len(feats)):
            w = -Prec[i,j]/max((Prec[i,i]*Prec[j,j])**0.5, 1e-12)
            if np.isfinite(w) and abs(w) > 0:
                rows.append({"Feature_A":feats[i], "Feature_B":feats[j], "IsingDirect":float(w)})
    return pd.DataFrame(rows)

def run_benchmark_compare_models(out_dir, n_rows=180, p_values=(1000,5000,20000), seed=123, boot_reps=1):
    from pathlib import Path
    import json
    import time
    import numpy as np
    import pandas as pd
    from .innovation_v3 import _plm_ising, glasso_ising_proxy
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    rows=[]
    for p in p_values:
        X = (rng.random((n_rows, p)) < 0.08).astype(np.uint8)
        # weak dependency blocks
        for b in range(0, min(p, 120), 12):
            base = (rng.random(n_rows) < 0.1).astype(np.uint8)
            for j in range(b, min(b+6, p)):
                flip = (rng.random(n_rows) < 0.05).astype(np.uint8)
                X[:,j] = np.bitwise_xor(base, flip)
        feat_cap = min(p, 40)  # keep benchmark stable in constrained environments
        Xsub = pd.DataFrame(X[:, :feat_cap], columns=[f"f{i}" for i in range(feat_cap)])
        methods=[]
        t=time.perf_counter(); plm=_plm_ising(Xsub, alpha=0.02); methods.append(("PLM-L1", plm, time.perf_counter()-t))
        try:
            t=time.perf_counter(); gdf=_glasso_edges_from_binary(Xsub, alpha=0.02); tg=time.perf_counter()-t
        except Exception:
            gdf=pd.DataFrame(columns=["Feature_A","Feature_B","IsingDirect"]); tg=float("nan")
        methods.append(("GraphicalLasso", gdf, tg))
        corr = np.corrcoef(Xsub.to_numpy(float), rowvar=False)
        ers=[]; feats=Xsub.columns.tolist()
        for i in range(len(feats)):
            for j in range(i+1, len(feats)):
                ers.append({"Feature_A":feats[i], "Feature_B":feats[j], "phi": float(corr[i,j]) if np.isfinite(corr[i,j]) else 0.0})
        t=time.perf_counter(); ldf=glasso_ising_proxy(pd.DataFrame(ers), feature_matrix=None, thr=0.0, alpha=0.02); methods.append(("LegacyPrecisionOnPhi", ldf, time.perf_counter()-t))
        refs={m:_edge_set_from_df(d) for m,d,_ in methods}
        for mname, d0, rt in methods:
            jaccs=[]
            for _ in range(int(max(1, boot_reps))):
                idx = rng.integers(0, n_rows, size=n_rows)
                Xb = Xsub.iloc[idx].reset_index(drop=True)
                try:
                    if mname=="PLM-L1":
                        db = _plm_ising(Xb, alpha=0.02)
                    elif mname=="GraphicalLasso":
                        db = _glasso_edges_from_binary(Xb, alpha=0.02)
                    else:
                        c2 = np.corrcoef(Xb.to_numpy(float), rowvar=False)
                        er=[]; fb=Xb.columns.tolist()
                        for i in range(len(fb)):
                            for j in range(i+1, len(fb)):
                                er.append({"Feature_A":fb[i], "Feature_B":fb[j], "phi": float(c2[i,j]) if np.isfinite(c2[i,j]) else 0.0})
                        db = glasso_ising_proxy(pd.DataFrame(er), feature_matrix=None, thr=0.0, alpha=0.02)
                    a = refs.get(mname,set()); b = _edge_set_from_df(db)
                    jaccs.append(len(a & b)/max(1, len(a | b)))
                except Exception:
                    jaccs.append(float("nan"))
            rows.append({
                "tool":"bact-assoc-net","benchmark":"ising_compare","p":int(p),"p_used":int(feat_cap),"method":mname,
                "runtime_sec": None if rt!=rt else round(float(rt),4),
                "n_edges": int(len(d0)) if d0 is not None else 0,
                "sparsity_edges_per_feature": round((len(d0)/max(1,feat_cap)) if d0 is not None else 0.0, 4),
                "stability_jaccard_mean": round(float(np.nanmean(jaccs)),4) if np.isfinite(np.nanmean(jaccs)) else None,
                "stability_jaccard_sd": round(float(np.nanstd(jaccs)),4) if np.isfinite(np.nanstd(jaccs)) else None,
            })
    out={"tool":"bact-assoc-net","report":"ising_model_comparison","results":rows}
    (out_dir/"benchmark_ising_model_comparison.json").write_text(json.dumps(out, indent=2))
    import pandas as pd
    df = pd.DataFrame(rows)
    if not df.empty:
        df.to_csv(out_dir/"benchmark_ising_model_comparison.csv", index=False)
        try:
            import matplotlib.pyplot as plt
            for y,fn in [("runtime_sec","benchmark_ising_runtime_comparison.png"),("sparsity_edges_per_feature","benchmark_ising_sparsity_comparison.png"),("stability_jaccard_mean","benchmark_ising_stability_comparison.png")]:
                fig = plt.figure(figsize=(6,4))
                for m, sub in df.groupby("method"):
                    plt.plot(sub["p"], sub[y], marker="o", label=m)
                plt.xscale("log"); plt.xlabel("p (synthetic features)"); plt.ylabel(y); plt.legend(); plt.tight_layout()
                fig.savefig(out_dir/fn, dpi=160); plt.close(fig)
        except Exception:
            pass
    return out
