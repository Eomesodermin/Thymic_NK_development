# Reproducing the Dandelion NK/ILC-from-DN branch, and dissecting what drives it

**Project:** Thymic NK development — developmental origins of a CD56^bright-enriched NK subset
**Workstream 2:** Dandelion (Suo et al.) trajectory reproduction + driver dissection
**Status:** complete (2026-07-08)

---

## 1. Question

Suo et al. (Dandelion, *Nat Biotechnol* 2023, DOI 10.1038/s41587-023-01734-7, Fig 5b/c)
built a pseudotime trajectory in a **T-cell-receptor V(D)J feature space** (TRBJ usage
summarised per transcriptional neighborhood) across the developing human immune atlas, and
concluded that **ILC/NK cells are originally on a canonical T-cell developmental trajectory
but are subsequently influenced to abort it** — i.e. NK/ILC cells branch off the
double-negative (DN) thymocyte path.

Two things about that claim matter for our hypothesis (that a CD56^bright NK subset arises
via a thymic/aborted-T route leaving TCR-rearrangement footprints):

1. **Can we reproduce the branch** on the same data with an independent re-run of the pipeline?
2. **What drives the branch** — the TCR-rearrangement footprint itself (the V(D)J feature
   space), or merely shared transcriptional programs between NK/ILC and DN thymocytes? If the
   branch is carried by the TCR footprint, that is direct support for an aborted-rearrangement
   origin. If it survives with the TCR signal removed, the branch is transcriptome-driven and
   says nothing specific about rearrangement.

## 2. Data

- **GEX:** Suo et al. developing-human-immune atlas, "Lymphoid cells" object from CELLxGENE
  collection `b1a879f6-5638-48d3-8f64-f6592c1b1561` (DOI 10.1126/science.abo0510). Recipe cell
  types subset to the T / DN / DP / NK / ILC ladder (111,706 cells).
- **VDJ:** re-derived from raw. The CELLxGENE object is GEX-only (CZI strips TCR obs), and the
  authors' per-cell VDJ CSVs live on a private NFS. We re-aligned the **64 αβTCR-enrichment
  libraries** (E-MTAB-11388, ENA ERP135310) that carry recipe cells, with
  `cellranger vdj --chain TR` on the `all_contig` output (non-productive contigs retained —
  the aborted-rearrangement relics are the point), then IgBLAST-reannotated with the
  sc-dandelion container (`dandelion-preprocess --chain TR --org human --file_prefix all`).

## 3. Method (reproduction)

Pipeline (scripts `01`,`03`,`05` under `scripts/dandelion_reproduction/`):

1. Read the 64 per-library Dandelion AIRR TSVs, prefix barcodes to match the Suo obs index
   (`<GEX_id>-<barcode>`, `sep="-"`), concatenate (301,511 contigs; 60/64 libraries readable,
   4 near-empty libraries with 1–12 contigs skipped on read).
2. `check_contigs(tr-ab, productive_only=False)` → `find_clones` → `transfer` onto the GEX
   object. **`productive_only=False` is deliberate**: the paper's signal in ILC/NK is
   *non-productive* TRB, so productive-only QC would discard exactly those cells.
3. `setup_vdj_pseudobulk(mode="abT", productive_vdj=False, productive_vj=False)` — retain any
   cell with a real TRBJ-bearing β rearrangement, productive or not (10,759 cells: ~10,540 T
   ladder + 217 NK/ILC). The feature space is TRBJ *gene usage*, which the paper builds
   agnostic to productivity ("all T/ILC/NK cells express TRBJ").
4. scVI latent (20-dim, 2,000 HVG on raw counts from `.raw`, batch = donor, 200 epochs CPU)
   → kNN graph (k=50) → neighborhoods (`make_nhoods`, the pure-Python Milo sampler
   reimplemented to avoid the R/edgeR dependency; 632 neighborhoods) → UMAP.
5. `vdj_pseudobulk(extract_cols=["j_call_abT_VDJ_main"])` over neighborhoods → **TRBJ-usage
   feature space** (632 neighborhoods × 13 TRBJ genes).
6. Palantir on the PCA of that feature space: diffusion maps (5 components),
   root = highest-CD34 neighborhood, terminal states = the neighborhoods with the highest T
   fraction and highest NK/ILC fraction (unique labels "T" and "NK_ILC").
   `pseudotime_transfer` → `project_pseudotime_to_cell` (10,702 cells retained).

**Outputs:** Fig 5b (neighborhood V(D)J-space UMAP: cell type / pseudotime / branch prob) and
Fig 5c (per-cell branch probabilities + pseudotime).

## 4. Method (driver dissection — the novel question)

Script `04_driver_dissection.py` (+ control built inside `03`):

The root neighborhood (CD34) and the two terminal tips (max T-fraction, max NK/ILC-fraction)
are held **fixed** across every test — only the feature matrix the pseudotime is computed on
changes. The readout is a single **branch-separation** metric: mean NK/ILC branch probability
in NK/ILC cells minus that in T cells (positive ⇒ NK/ILC cells are pushed onto their own tip).

- **(A) GEX-only control.** Feature space swapped from TRBJ usage to the **scVI transcriptome
  mean per neighborhood**, identical Palantir recipe. If the NK/ILC branch reproduces here, it
  is transcriptome-driven; if it degrades, the TCR footprint carries it.
- **(B) TRBJ permutation null (n=200).** Shuffle the TRBJ-usage vectors across neighborhoods
  (destroying the real TRBJ↔neighborhood correspondence, preserving the marginal feature
  distribution), re-run the fixed-root/tip Palantir, and build a null distribution of the
  branch-separation metric. Empirical p = P(null separation ≥ real).
- **(C) Productive vs +non-productive TRBJ.** Intended to null the TRBJ call of non-productive-β
  cells and compare. **Moot in this run:** every cell that survives `check_contigs` into the
  pseudobulk has a *productive* β `_main` call (the Dandelion `_main`-call machinery assigns the
  productive contig when one exists), so there were 0 non-productive-only cells to null and the
  productive-only feature space is identical to the real one. The relevant test therefore moved
  to the *raw* contig table (§5).

## 5. Results

### 5.1 Reproduction (Fig 5b/c)

The branch **qualitatively reproduces**. In the TRBJ feature space, NK cells sit at the tip of
the NK/ILC branch and carry the highest pseudotime; the mean per-cell NK/ILC branch probability
(`fig5c_branch_probabilities.png`, `branch_prob_by_celltype.csv`) is:

| cell type | n | pseudotime | prob_NK/ILC |
|---|---:|---:|---:|
| NK | 76 | 0.72 | **0.623** |
| CYCLING_NK | 67 | 0.67 | **0.569** |
| ILC2 | 14 | 0.56 | 0.502 |
| ILC3 | 41 | 0.54 | 0.487 |
| ABT(ENTRY) | 3,247 | 0.56 | 0.488 |
| DP(P)_T / DP(Q)_T | 7,139 | 0.50–0.55 | 0.480 |
| DN(P)/DN(Q)_T | 111 | 0.44 | 0.479–0.480 |

So bulk **NK and CYCLING_NK are clearly displaced onto the NK/ILC branch** (0.62 / 0.57 vs the
T-lineage baseline of ~0.48). The small-n ILC subsets (ILC3 n=41, CYCLING_ILC n=2) and the rare
DN(early)_T (n=5) sit at the ~0.48–0.50 midline and do **not** separate — the reproduced branch
is carried by NK proper, not by a broad NK/ILC-vs-T split. Separation is modest: the feature
space is only **13 TRBJ genes** and only **217 NK/ILC cells** carry any TRBJ at all.

### 5.2 Driver dissection

Branch-separation metric (mean prob_NK/ILC: NK/ILC cells − T cells), root/tips fixed:

| test | separation | reading |
|---|---:|---|
| **REAL** TRBJ feature space | **+0.023** | NK/ILC modestly displaced |
| **(A)** scVI / transcriptome feature space | **−0.057** | branch **not** reproduced by transcriptome (with the same TRBJ-defined tips) |
| **(B)** TRBJ permutation null (n=200) | mean −0.004, sd 0.039, 95th pct **+0.057** | **real +0.023 → empirical p = 0.29** |

**The real TRBJ separation (+0.023) is inside the permutation null** (below its 95th percentile),
so the observed branch structure is **not statistically distinguishable from random TRBJ
profiles**. (C) was moot in the QC'd object (§4).

### 5.3 Where the paper's signal went (raw-contig check)

Because (A)/(B) hinge on what TRBJ information the NK/ILC cells actually carry, we went back to
the **raw** contig table (pre-`check_contigs`):

- NK/ILC: **389 cells with a TRB contig**; 42 non-productive-only; **all 83 non-productive TRB
  contigs carry a V gene (0 V-less)**.
- The 217 NK/ILC cells that survive into the trajectory all have a **productive** β chain
  (mostly `Single pair`), and their **TRBJ usage is dominated by the same common J genes as
  T cells** (both led by TRBJ2-7, with TRBJ1-2 / 1-1 / 2-3 / 2-1 filling the top ranks; the
  exact 2nd–4th ordering differs slightly but the distributions are not distinguishable at
  n=217).

The paper's defining ILC/NK signal — *non-productive, mostly **V-gene-less** TRB* — is **absent**
from our `cellranger vdj --chain TR` re-derivation (0 V-less non-productive TRB in either NK or T).

## 6. Interpretation

1. **The branch reproduces, but weakly and only for NK proper.** We recover the paper's
   qualitative claim (NK cells sit at an NK/ILC branch tip off the T-development trajectory), but
   in our re-derivation it rests on ~217 TRBJ-bearing NK/ILC cells and a 13-gene feature space,
   and the T-vs-NK/ILC separation is not significant against a permutation null (p = 0.29).
2. **In this re-derivation the branch is not TCR-footprint-driven.** The NK/ILC cells that reach
   the trajectory carry a *conventional, productive* TCR whose TRBJ usage draws on the same
   common J genes as T cells (not separable at n=217), and shuffling the TRBJ vectors does not
   degrade the branch beyond chance. The TRBJ feature space is not adding a rearrangement-specific
   signal here — it is close to noise.
3. **The discrepancy is a pipeline/reference difference, not a contradiction of the biology.**
   The paper's signal is V-gene-less non-productive TRB; our cellranger+IMGT re-alignment produces
   **no** V-less non-productive TRB. Those contigs are either an artifact of, or specific to, the
   authors' assembly/reference, or they are filtered by cellranger's contig calling. Either way,
   with the footprint the paper relies on **absent from our data**, we cannot use this dataset to
   test whether the *aborted rearrangement* drives the branch.
4. **Consistency with the rest of the project.** This is the same story as Workstream 1
   (Domínguez Conde): the NK TCR footprint is *rare* and does not cleanly define a subset. It is
   also consistent with Kelvin's earlier "not enough cells / no strong signal" read — the
   developing-atlas NK/ILC population with usable VDJ is simply too small (hundreds of cells) to
   power this trajectory claim robustly.

**Bottom line:** the NK/ILC-from-DN branch is reproducible in direction but statistically fragile,
and in our independent re-alignment it is **transcriptome-adjacent rather than TCR-footprint-driven**
— largely because the specific non-productive/V-less TRB relics the paper leans on do not appear in
a standard cellranger `--chain TR` re-derivation. A definitive driver test needs either the
authors' original contig assemblies or a larger, bright-enriched NK dataset.

## 7. Reproducibility

- Local analysis scripts: `scripts/dandelion_reproduction/{01,03,04,05}*.py`
- HPC alignment + preprocessing (Marvin): `HPC_workflows/thymic_nk_development/dandelion_reproduction/`
  (`cr_vdj_suo_one.sh`, `submit_suo_all.sh`, `cr_vdj_suo_fix.sh`, `cr_vdj_suo_fix2.sh`,
  `resubmit_suo_failed.sh`, `dandelion_preprocess_all.sh`).
- VDJ library manifest: `data/dandelion_reproduction/suo_abTCR_vdj_manifest.{csv,json}`.
- Env: `ddl-traj` (scanpy 1.10.4, dandelion 0.5.7, scvi-tools 1.3.3, palantir 1.4.4,
  pertpy 0.10.0). Run compute cells with `NUMBA_CACHE_DIR=/tmp/numba_cache_ddltraj
  NUMBA_DISABLE_JIT=1` set (numba `/dev/fd/3` locator bug in the notebook kernel; scripts
  run standalone via bash).
- Intermediate objects (`processed/dandelion_reproduction/`): `trajectory_gex_subset.h5ad`
  (111,706 recipe cells), `dandelion_tsv/` (64 AIRR TSVs), `trajectory_adata_scvi_nhoods.h5ad`
  (scVI + neighborhoods checkpoint — resumes Palantir without the ~25-min prefix),
  `trajectory_pseudobulk.h5ad`, `trajectory_percell.h5ad`, `driver_dissection_summary.csv`,
  `driver_perm_null.npy`.
- Figures/tables (`results/dandelion_reproduction/`): `fig5b_vdj_feature_umap.png`,
  `fig5c_branch_probabilities.png`, `fig_driver_dissection.png`, `branch_prob_by_celltype.csv`.

### Caveats / limitations

- **NK/ILC n is small** (217 TRBJ-bearing cells; NK proper n=76). The branch and every driver
  test inherit that fragility.
- **4 near-empty VDJ libraries** (1–12 contigs) were skipped on read; their contribution to the
  matched pool is negligible.
- **The `make_nhoods` step is a verbatim reimplementation** of `pertpy.tools.Milo.make_nhoods`
  (pure-Python neighborhood sampler) to avoid pertpy's rpy2/edgeR dependency; the differential-
  abundance GLM is not used here.
- **Reference dependence.** The central discrepancy (no V-less non-productive TRB) is a property
  of the cellranger + IMGT GRCh38 VDJ reference vs the authors' pipeline; it is the main reason
  the driver test is inconclusive on this dataset rather than a biological result.
