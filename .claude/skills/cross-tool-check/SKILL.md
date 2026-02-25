---
name: cross-tool-check
description: Verify anti-salami scope boundary and consistency with sister tools
---

Check that this tool (bact-assoc-net) stays within its defined scope and does not duplicate functionality from sister tools.

## Steps

1. **Read scope boundary**: Parse `docs/anti_salami_checklist.json` and verify:
   - `scope_boundary` lists what this tool does NOT do
   - `novelty_axis` describes the unique contribution

2. **Grep for scope violations** — search `src/bactassocnet/` for terms that belong to sister tools:
   - Clustering terms: `k_modes`, `kmodes`, `consensus_matrix`, `NVI`, `CKA`, `kernel_fusion`
   - Phylo terms: `pagel`, `fritz_purvis`, `ancestral_state`, `newick`, `phylo_logistic`, `vcv_matrix`
   - MDR terms: `mdr_class`, `xdr`, `pdr`, `pc_algorithm`, `shapley`, `evpi`, `hypergraph`

3. **Check shared interfaces** (if sister repos exist at `../bact-*`):
   - Config schema version consistency (all should be v1.1)
   - CSV format compatibility (all use 0/1/NA binary)
   - CLI flag naming conventions

Report violations as WARN (potential overlap) or FAIL (clear duplication). Report PASS if scope is clean.
