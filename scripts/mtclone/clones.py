"""clones.py — clone inference from a binarized cell x variant matrix.

Two callers, both writing obs['clone_id'] (int; -1 = unassigned/singleton):

  method='graph'  (default): build a cell-cell graph weighted by the number of shared
                   informative variants, drop weak edges, take communities (Leiden if
                   available, else connected components) as clones. This is the
                   mgatk/Liu-style approach.
  method='variant_group': cells sharing a rare high-confidence variant form one clone;
                   used as an independent cross-check of the graph caller.

Clones are donor-private. Always call per donor (call_clones_per_donor) unless the object is
already single-donor — sharing an mtDNA variant across donors is coincidental, not lineage.
"""
from __future__ import annotations

import warnings

import numpy as np
import scipy.sparse as sp
import anndata as ad


def _binary(adata: ad.AnnData) -> sp.csr_matrix:
    if "binary" not in adata.layers:
        raise KeyError("no layers['binary']; run qc.binarize_heteroplasmy / "
                       "select_informative_variants first.")
    return adata.layers["binary"].tocsr().astype(np.float32)


def call_clones(
    adata: ad.AnnData,
    *,
    method: str = "graph",
    min_shared_variants: int = 1,
    edge_weight_cutoff: float = 0.5,
    min_clone_size: int = 2,
    resolution: float = 1.0,
    inplace: bool = True,
) -> ad.AnnData:
    """Assign obs['clone_id']. See module docstring for methods.

    graph caller:
      - S = B @ B.T gives # shared alt variants between each cell pair;
      - keep edges with shared >= min_shared_variants AND jaccard >= edge_weight_cutoff;
      - communities via Leiden (python-igraph+leidenalg) or connected components fallback.
    """
    a = adata if inplace else adata.copy()
    B = _binary(a)

    if method == "graph":
        labels = _graph_clones(B, min_shared_variants, edge_weight_cutoff, resolution)
    elif method == "variant_group":
        labels = _variant_group_clones(B)
    else:
        raise ValueError(f"unknown method {method!r}")

    labels = _apply_min_size(labels, min_clone_size)
    a.obs["clone_id"] = labels.astype(int)
    a.uns["clone_params"] = dict(method=method, min_shared_variants=min_shared_variants,
                                 edge_weight_cutoff=edge_weight_cutoff,
                                 min_clone_size=min_clone_size, resolution=resolution)
    return a


def call_clones_per_donor(adata: ad.AnnData, *, donor_key: str = "donor", **kwargs) -> ad.AnnData:
    """Run call_clones separately within each donor; clone_ids are made globally unique."""
    a = adata.copy()
    a.obs["clone_id"] = -1
    offset = 0
    for donor, idx in a.obs.groupby(donor_key, observed=True).groups.items():
        sub = a[idx].copy()
        sub = call_clones(sub, inplace=True, **kwargs)
        lab = sub.obs["clone_id"].values.copy()
        assigned = lab >= 0
        lab[assigned] += offset
        a.obs.loc[idx, "clone_id"] = lab
        if assigned.any():
            offset = int(a.obs["clone_id"].max()) + 1
    a.obs["clone_id"] = a.obs["clone_id"].astype(int)
    a.uns.setdefault("clone_params", {})["per_donor"] = True
    return a


# ----------------------------------------------------------------------------- graph caller
def _graph_clones(B: sp.csr_matrix, min_shared: int, jac_cut: float, resolution: float):
    n = B.shape[0]
    shared = (B @ B.T).tocoo()  # # shared alt variants
    per_cell = np.asarray(B.sum(axis=1)).ravel()  # # alt variants per cell

    rows, cols, weights = [], [], []
    for i, j, s in zip(shared.row, shared.col, shared.data):
        if i >= j:
            continue
        if s < min_shared:
            continue
        union = per_cell[i] + per_cell[j] - s
        jac = s / union if union > 0 else 0.0
        if jac >= jac_cut:
            rows.append(i); cols.append(j); weights.append(jac)

    if not rows:
        return -np.ones(n, dtype=int)

    try:
        import igraph as ig
        import leidenalg
        g = ig.Graph(n=n, edges=list(zip(rows, cols)))
        g.es["weight"] = weights
        part = leidenalg.find_partition(
            g, leidenalg.RBConfigurationVertexPartition,
            weights="weight", resolution_parameter=resolution, seed=0,
        )
        labels = np.array(part.membership, dtype=int)
        # singletons (degree-0 vertices) -> -1
        deg = np.zeros(n, dtype=int)
        for i, j in zip(rows, cols):
            deg[i] += 1; deg[j] += 1
        labels[deg == 0] = -1
        return _compact(labels)
    except ImportError:
        warnings.warn("leidenalg/igraph unavailable; using connected components.")
        adj = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
        adj = adj + adj.T
        ncomp, comp = sp.csgraph.connected_components(adj, directed=False)
        # mark isolated vertices (no edges) as -1
        deg = np.asarray((adj > 0).sum(axis=1)).ravel()
        comp = comp.astype(int)
        comp[deg == 0] = -1
        return _compact(comp)


def _variant_group_clones(B: sp.csr_matrix):
    """Cells sharing their rarest common variant form a clone (independent cross-check)."""
    n, v = B.shape
    per_var = np.asarray(B.sum(axis=0)).ravel()  # cells per variant
    labels = -np.ones(n, dtype=int)
    order = np.argsort(per_var)  # rarest variants first
    Bc = B.tocsc()
    clone = 0
    for vi in order:
        if per_var[vi] < 2:
            continue
        cells = Bc.getcol(vi).nonzero()[0]
        unassigned = cells[labels[cells] == -1]
        if len(unassigned) >= 2:
            labels[unassigned] = clone
            clone += 1
    return _compact(labels)


# ----------------------------------------------------------------------------- utils
def _apply_min_size(labels: np.ndarray, min_size: int):
    labels = labels.copy()
    uniq, counts = np.unique(labels[labels >= 0], return_counts=True)
    too_small = set(uniq[counts < min_size].tolist())
    if too_small:
        labels[np.isin(labels, list(too_small))] = -1
    return _compact(labels)


def _compact(labels: np.ndarray):
    """Relabel assigned clones to 0..K-1 contiguous; keep -1."""
    labels = labels.astype(int).copy()
    uniq = sorted(set(labels[labels >= 0].tolist()))
    remap = {old: new for new, old in enumerate(uniq)}
    out = np.array([remap.get(x, -1) for x in labels], dtype=int)
    return out


def clone_qc(adata: ad.AnnData, clone_key: str = "clone_id") -> dict:
    """Summary stats for a clone assignment."""
    lab = adata.obs[clone_key].values
    assigned = lab[lab >= 0]
    uniq, counts = np.unique(assigned, return_counts=True)
    return dict(
        n_cells=int(adata.n_obs),
        n_clones=int(len(uniq)),
        n_assigned=int(assigned.size),
        frac_assigned=float(assigned.size / adata.n_obs) if adata.n_obs else 0.0,
        n_singletons_or_unassigned=int((lab < 0).sum()),
        max_clone_size=int(counts.max()) if counts.size else 0,
        median_clone_size=float(np.median(counts)) if counts.size else 0.0,
    )
