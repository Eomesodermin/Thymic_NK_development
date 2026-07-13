# Thymic NK development — developmental origins of NK cells

**Do NK cells carry a TCR‑rearrangement footprint of a thymic / aborted‑T origin, and are CD56^bright and CD56^dim NK developmentally distinct?**

## 📄 Full project report

**→ Read the comprehensive report: [`report.md`](report.md)** — research questions, hypotheses, datasets, every analysis stream, results, interpretation, and honest limits.

A standalone rendered version with all figures embedded is at **[`report.html`](report.html)** (self‑contained — open directly in any browser).

### Headline findings

| Question | Verdict | Confidence |
|----------|---------|:----------:|
| Are CD56^bright and CD56^dim NK one lineage or two? | **One clonal lineage / continuum** (mtDNA sharing) | high |
| Do NK carry a non‑productive TCR footprint? | Yes, but weak (~0.4–0.5%) and **not** CD56^bright‑enriched | moderate |
| Do NK trace a thymic/DN‑branch trajectory? | Reproducible in direction but statistically fragile (perm p = 0.29) | low |
| Are tumor‑infiltrating NK clonally restricted? | **Tumor‑type‑dependent** — restricted in HGSC ovarian, passive in NSCLC | moderate–high |
| Are NK clonally closer to T or myeloid? | Myeloid‑leaning (5/6 patients, p ≈ 0.06); human evidence mixed | low–moderate |

The strongest result — that bright and dim NK are one clonal lineage — argues against the two‑origin premise as originally framed. See the report for the full picture across both the TCR‑footprint and mtDNA lineage‑tracing streams.

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
