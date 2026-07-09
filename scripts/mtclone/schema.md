# `mtclone` canonical data schema

Every ingest adapter returns an `AnnData` with the **identical** structure below, so all
downstream code (QC, clone inference, metrics, classification) is assay-agnostic. This is the
single contract the whole package depends on.

## The object

`adata` : `anndata.AnnData`

- **`adata.X`** — cell × variant **heteroplasmy** matrix, `scipy.sparse` (CSR), float32,
  values in **[0, 1]** (fraction of mtDNA molecules at that cell carrying the alt allele).
  0 = not detected / reference. This is the *continuous* matrix; binarization happens in `qc`.
- **`adata.layers["binary"]`** — (added by `qc.binarize_heteroplasmy`) uint8 {0,1}, alt call
  above the binarization threshold. Absent until QC is run.

## `adata.obs` (one row per cell) — required columns

| column        | dtype    | units / meaning                                                        |
|---------------|----------|------------------------------------------------------------------------|
| `cell_id`     | str      | unique barcode, globally unique across samples (prefix with sample)    |
| `donor`       | str      | **patient/individual id** — clones are donor-private; never pool donors |
| `sample`      | str      | library/sample id (a donor may have several)                           |
| `tissue`      | str      | e.g. `tumor`, `normal`, `PBMC`, `BMMC`, `HSPC`                          |
| `site`        | str      | coarse site for overlap tests: `tumor` / `normal` / `blood` / `marrow` |
| `assay`       | str      | `mtscATAC` \| `redeem` \| `maester` \| `mgatk`                          |
| `dataset`     | str      | source accession, e.g. `GSE302113`, `GSE219014`                        |
| `coverage`    | float32  | per-cell **mean mtDNA coverage** (×). From depth file if shipped, else |
|               |          | derived (mean per-cell nonzero-variant read support / fragments proxy) |
| `coverage_source` | str  | `depth_file` \| `derived_from_matrix` \| `fragments` — provenance flag  |

### Added later by pipeline stages (not required at ingest)
| column        | added by      | meaning                                                    |
|---------------|---------------|------------------------------------------------------------|
| `cell_type`   | `classify`    | NK / T / B / myeloid / HSPC / tumor / …                     |
| `nk_subset`   | `classify`    | `bright` \| `dim` \| `unassigned` (NK only)                 |
| `clone_id`    | `clones`      | clone label; `-1` / NaN = unassigned/singleton             |

## `adata.var` (one row per variant) — required columns

| column                   | dtype   | units / meaning                                             |
|--------------------------|---------|-------------------------------------------------------------|
| `variant_id`             | str     | canonical `chrM:POS:REF>ALT` (e.g. `chrM:3243:A>G`)         |
| `pos`                    | int     | 1-based rCRS position                                       |
| `ref`, `alt`             | str     | reference / alternate base                                  |
| `pseudobulk_heteroplasmy`| float32 | mean heteroplasmy across all cells (0–1) — for the ceiling filter |
| `n_cells_detected`       | int     | # cells with heteroplasmy > 0                               |

### Optional `adata.var` columns (present when the source provides them; not enforced)
| column                   | dtype   | units / meaning                                             |
|--------------------------|---------|-------------------------------------------------------------|
| `strand_correlation`     | float32 | strand-balance stat where available (mtscATAC/mgatk), else NaN. Adapters always create it (NaN when absent), but `validate_schema` does **not** require it. |

## `adata.uns` — metadata bag
- `uns['mtclone_schema_version']` = `"1.0"`
- `uns['thresholds']` — dict of QC parameters actually applied (filled by `qc`)
- `uns['clone_params']` — dict of clone-caller parameters (filled by `clones`)

## Variant id convention
Always `chrM:POS:REF>ALT`, 1-based, rCRS (NC_012920.1) coordinates, bases uppercase.
Adapters must normalize whatever the source uses (mgatk `POS_REF_ALT`, ReDeeM `Variants`
strings, Liu `variant` column) to this exact form so variants are comparable **across
datasets** — essential if the optional re-calling module ever puts both on one coordinate frame.

## Invariants (checked by `io.validate_schema`)
1. `X` is sparse, float32, `0 ≤ X ≤ 1`.
2. All required `obs` and `var` columns present, no all-NaN required column.
3. `cell_id` unique; `donor` non-null for every cell.
4. `var_names == var['variant_id']` and are unique.
5. `uns['mtclone_schema_version']` set.
