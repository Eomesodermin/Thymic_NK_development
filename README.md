# Thymic NK development — developmental origins of NK cells

**Do NK cells carry a TCR‑rearrangement footprint of a thymic / aborted‑T origin, and are CD56^bright and CD56^dim NK developmentally distinct?**

## 📄 Full project report

**→ Read the comprehensive report: [`report.md`](report.md)** — research questions, hypotheses, datasets, every analysis stream, results, interpretation, and honest limits.

A standalone rendered version with all figures embedded is at **[`report.html`](report.html)** (self‑contained — open directly in any browser).

### The project in one paragraph

The idea: a subset of NK cells (enriched in CD56^bright) might arise via a **thymic / aborted‑T pathway** and carry a relic of it as **non‑productive TCR rearrangements**. The first approach — reading that footprint from single‑cell TCR (mRNA/VDJ) data — turned out to be **structurally the wrong assay**: non‑productive rearrangements are NMD‑silenced/degraded, so mRNA‑based VDJ measures a deflated floor, not the true rate; the footprint lives in genomic DNA where this assay can't see it. So the effort **pivoted to somatic‑mtDNA lineage tracing (mtSNP)** — but those datasets are **sparse** (no thymic dataset exists at all), and from the one decent tumor dataset there *appears* to be a clonal‑restriction signal in **one tumor type only (high‑grade serous ovarian), in ~2 patients** — too few to be confident. The one place the data spoke with any strength: bright and dim NK are **not cleanly two lineages** — though even that rests on very few informative clones.

### Headline findings

| Question | Verdict | Confidence |
|----------|---------|:----------:|
| Can mRNA‑based VDJ measure the NK TCR footprint? | **No** — NMD silences the transcripts; assay sees a deflated floor | — |
| Do NK carry a non‑productive TCR footprint (as far as we can see)? | Detectable but weak (~0.4–0.5%), **not** CD56^bright‑enriched, under‑measured | moderate |
| Do NK trace a thymic/DN‑branch trajectory? | Reproducible in direction but not significant (perm p = 0.29) | low |
| Are CD56^bright and CD56^dim NK two lineages? | **No clear evidence of two** — clean segregation excluded, but built on ~38 informative clones | moderate (strongest, power‑limited) |
| Are tumor‑infiltrating NK clonally restricted? | Apparent in **HGSC ovarian only** (not lung), but ~2 patients | low–moderate |
| Are NK clonally closer to T or myeloid? | Myeloid‑leaning (5/6 patients, p ≈ 0.06); human evidence mixed | low–moderate |

An honest negative on the headline hypothesis: the direct assay is defeated by measurement biology and the alternative by data scarcity. See the report for the full arc across both the TCR‑footprint and mtSNP lineage‑tracing streams, and §10 for what would be needed to test it properly.

---

## Repository

Project repository — code and configs only. Large/raw/processed data and HPC run outputs
live outside this repo at `~/Documents/HPC_data/Thymic_NK_development/` (not version-controlled).

## Structure

- `report.md` / `report.html` — the full project report (start here)
- `docs/` — supporting memos: literature review, NK↔T origin probe, Six 2020 re‑analysis
- `scripts/mtclone/` — frozen mtDNA clone‑calling package (shared by both mtDNA plans) + tests
- `scripts/analysis/development_brightdim/` — Plan 2: bright‑vs‑dim one‑lineage test (ingestion + steps 3–7)
- `scripts/analysis/tumor_clonal_restriction/` — Plan 3: tumor NK clonal restriction (01_download → 08_robustness)
- `scripts/dandelion_reproduction/` — Stream A: Dandelion NK‑from‑DN trajectory reproduction
- `scripts/analysis/footprint_pipeline.py`, `harvest_footprint.py` — Q1/Q2 bright‑vs‑dim TCR footprint
- `results/` — figures, per‑stream result memos, per‑plan rendered HTML reports
- `report_assets/` — figures embedded in the report
- `configs/`, `notebooks/` — small parameter files and exploratory notebooks

HPC job scripts (cellranger VDJ, sc‑dandelion) live in the companion
[`HPC_workflows`](https://github.com/Eomesodermin/HPC_workflows) repo under `thymic_nk_development/`.

## HPC

See [`../HPC_workflows/marvin_hpc_reference.md`](../HPC_workflows/marvin_hpc_reference.md) for
cluster access, partitions, and submission conventions. This project's SLURM account:
`ag_iei_abdullah`.

Harvested job outputs land in `~/Documents/HPC_data/Thymic_NK_development/hpc_runs/`.

## Remote

`git remote add origin https://github.com/Eomesodermin/Thymic_NK_development.git`
