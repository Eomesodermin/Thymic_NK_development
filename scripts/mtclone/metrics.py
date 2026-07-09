"""metrics.py — clonal-structure statistics shared by both scientific questions.

- clone_size_distribution / gini / shannon  : how clonal is a group of cells?
- clone_sharing_matrix                       : fraction of clones shared across cell groups
                                               (the Liu cross-group statistic; group = cell
                                               type for the tumor Q, or NK subset for the
                                               development Q).
- between_vs_within_sharing                  : the development-question core test.
- permutation_null                           : donor-stratified label shuffle -> p-value for
                                               any statistic.

Donor stratification is mandatory in the null: clones are donor-private, so shuffling labels
across donors would manufacture structure. Pass stratify_by='donor'.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import pandas as pd
import anndata as ad


# ----------------------------------------------------------------------------- diversity
def clone_size_distribution(adata: ad.AnnData, *, group=None, clone_key="clone_id") -> pd.Series:
    """Clone sizes (assigned clones only). If `group` (obs col) given, restrict to those cells."""
    obs = adata.obs
    if group is not None:
        col, val = group
        obs = obs[obs[col] == val]
    lab = obs[clone_key].values
    lab = lab[lab >= 0]
    if lab.size == 0:
        return pd.Series(dtype=int)
    uniq, counts = np.unique(lab, return_counts=True)
    return pd.Series(counts, index=uniq, name="clone_size").sort_values(ascending=False)


def gini(sizes) -> float:
    """Gini coefficient of a clone-size distribution. 0 = even, ->1 = one dominant clone."""
    x = np.sort(np.asarray(sizes, dtype=float))
    n = x.size
    if n == 0 or x.sum() == 0:
        return np.nan
    cum = np.cumsum(x)
    return float((n + 1 - 2 * (cum / cum[-1]).sum()) / n)


def shannon(sizes) -> float:
    """Shannon entropy (nats) of the clone-size distribution."""
    x = np.asarray(sizes, dtype=float)
    x = x[x > 0]
    if x.size == 0:
        return np.nan
    p = x / x.sum()
    return float(-(p * np.log(p)).sum())


def normalized_shannon(sizes) -> float:
    """Shannon / log(n_clones): 1 = maximally even, 0 = single clone. (Pielou evenness)."""
    x = np.asarray(sizes, dtype=float)
    x = x[x > 0]
    if x.size <= 1:
        return np.nan if x.size == 0 else 0.0
    return float(shannon(x) / np.log(x.size))


def clonality_summary(adata: ad.AnnData, *, group=None, clone_key="clone_id") -> dict:
    sizes = clone_size_distribution(adata, group=group, clone_key=clone_key).values
    return dict(n_clones=int(sizes.size), n_cells=int(sizes.sum()),
                gini=gini(sizes), shannon=shannon(sizes),
                normalized_shannon=normalized_shannon(sizes),
                max_clone_frac=float(sizes.max() / sizes.sum()) if sizes.sum() else np.nan)


# ----------------------------------------------------------------------------- sharing
def clone_sharing_matrix(adata: ad.AnnData, *, group_key: str, clone_key="clone_id",
                         min_clone_size=2) -> pd.DataFrame:
    """Fraction of clones shared between each pair of groups.

    For each clone, the set of groups it contains is computed. Entry (A,B) = fraction of
    clones present in group A that are also present in group B (directional; diagonal = 1).
    This is the cross-group statistic used to gauge lineage relatedness.
    """
    obs = adata.obs[[group_key, clone_key]].copy()
    obs = obs[obs[clone_key] >= 0]
    groups = sorted(obs[group_key].dropna().unique().tolist())
    # clone -> set(groups)
    clone_groups = obs.groupby(clone_key)[group_key].agg(lambda s: set(s.dropna()))
    # size filter
    sizes = obs.groupby(clone_key).size()
    clone_groups = clone_groups[sizes >= min_clone_size]

    M = pd.DataFrame(0.0, index=groups, columns=groups)
    for A in groups:
        clones_in_A = [c for c, gs in clone_groups.items() if A in gs]
        nA = len(clones_in_A)
        for B in groups:
            if nA == 0:
                M.loc[A, B] = np.nan
            else:
                shared = sum(1 for c in clones_in_A if B in clone_groups[c])
                M.loc[A, B] = shared / nA
    return M


def between_vs_within_sharing(adata: ad.AnnData, group_a, group_b, *, group_key: str,
                              clone_key="clone_id", min_clone_size=2) -> dict:
    """Core development-question test: are two groups (e.g. bright/dim NK) clonally mixed?

    Returns:
      frac_mixed_clones  : of clones containing >=1 cell from A or B, fraction containing BOTH.
      n_clones_a, n_clones_b, n_clones_mixed
      A high frac_mixed => shared lineages (continuum). ~0 => disjoint (two lineages).
    """
    obs = adata.obs[[group_key, clone_key]].copy()
    obs = obs[obs[clone_key] >= 0]
    obs = obs[obs[group_key].isin([group_a, group_b])]
    clone_groups = obs.groupby(clone_key)[group_key].agg(lambda s: set(s.dropna()))
    sizes = obs.groupby(clone_key).size()
    clone_groups = clone_groups[sizes >= min_clone_size]

    has_a = [group_a in gs for gs in clone_groups]
    has_b = [group_b in gs for gs in clone_groups]
    mixed = [(a and b) for a, b in zip(has_a, has_b)]
    n_total = len(clone_groups)
    return dict(
        frac_mixed_clones=float(np.mean(mixed)) if n_total else np.nan,
        n_clones_total=int(n_total),
        n_clones_a_only=int(sum(a and not b for a, b in zip(has_a, has_b))),
        n_clones_b_only=int(sum(b and not a for a, b in zip(has_a, has_b))),
        n_clones_mixed=int(sum(mixed)),
    )


# ----------------------------------------------------------------------------- null model
def permutation_null(
    adata: ad.AnnData,
    statistic: Callable[[ad.AnnData], float],
    *,
    group_key: str,
    stratify_by: str = "donor",
    n: int = 1000,
    seed: int = 0,
) -> dict:
    """Permute `group_key` labels WITHIN each `stratify_by` level and recompute `statistic`.

    `statistic` takes an AnnData and returns a float (e.g. a lambda wrapping
    between_vs_within_sharing). Returns observed, null mean/sd, and a two-sided-ish p (fraction
    of null >= observed, and <= observed; the relevant tail depends on the question).
    """
    rng = np.random.default_rng(seed)
    obs_val = statistic(adata)

    a = adata.copy()
    labels = a.obs[group_key].values.copy()
    strat = a.obs[stratify_by].values
    idx_by_stratum = {s: np.where(strat == s)[0] for s in np.unique(strat)}

    null = np.empty(n, dtype=float)
    for k in range(n):
        perm = labels.copy()
        for s, idx in idx_by_stratum.items():
            perm[idx] = rng.permutation(labels[idx])
        a.obs[group_key] = perm
        null[k] = statistic(a)

    null = null[~np.isnan(null)]
    if null.size == 0:
        return dict(observed=obs_val, null_mean=np.nan, null_sd=np.nan,
                    p_greater=np.nan, p_less=np.nan, n_null=0)
    p_greater = float((np.sum(null >= obs_val) + 1) / (null.size + 1))
    p_less = float((np.sum(null <= obs_val) + 1) / (null.size + 1))
    return dict(observed=float(obs_val), null_mean=float(null.mean()),
                null_sd=float(null.std()), p_greater=p_greater, p_less=p_less,
                n_null=int(null.size), null=null)
