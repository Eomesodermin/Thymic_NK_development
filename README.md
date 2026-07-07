# Thymic_NK_development

Project repository — code and configs only. Large/raw/processed data and HPC run outputs
live outside this repo at `~/Documents/HPC_data/Thymic_NK_development/` (not version-controlled).

## Structure

- `scripts/preprocessing/` — data preprocessing scripts
- `scripts/analysis/` — downstream analysis scripts
- `scripts/hpc_jobs/` — Marvin SLURM job scripts for this project (built from templates in
  [`../HPC_workflows/templates/`](../HPC_workflows/templates/))
- `notebooks/` — exploratory notebooks
- `configs/` — small parameter files (YAML/JSON), not data

## HPC

See [`../HPC_workflows/marvin_hpc_reference.md`](../HPC_workflows/marvin_hpc_reference.md) for
cluster access, partitions, and submission conventions. This project's SLURM account:
`ag_iei_abdullah`.

Harvested job outputs land in `~/Documents/HPC_data/Thymic_NK_development/hpc_runs/`.

## Remote

`git remote add origin https://github.com/Eomesodermin/Thymic_NK_development.git`
