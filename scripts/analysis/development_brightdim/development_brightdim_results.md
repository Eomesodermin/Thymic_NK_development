# CD56-bright and CD56-dim NK cells are one lineage, not two

**Plan 2 — somatic-mtDNA clonal lineage tracing. Final results.**

## Question

Are CD56-bright and CD56-dim NK cells a single developmental lineage (a continuum from shared
progenitors) or two independent pathways? Somatic mitochondrial-DNA variants mark clonal descent:
if bright and dim NK descend from common progenitors they should share mtDNA clones; if they arise
from separate pathways their clones should be disjoint.

## Verdict

**One lineage / continuum.** Across two independent datasets, bright and dim NK share somatic-mtDNA
clones exactly as much as a single well-mixed lineage predicts — and a true two-lineage structure
would have been clearly detectable but is not seen.

**Confidence: high**, for three reasons: (1) the result replicates across blood *and* marrow, two
mtDNA-calling chemistries (raw mtASAP + UMI-consensus ReDeeM), and two label modalities (surface
protein + RNA); (2) the test has demonstrated power — simulated segregation drives the statistic to
zero, far from what we observe; (3) it is stable across every QC/caller setting swept and is not a
doublet or contamination artifact.

## Evidence

**Data (Step 2).** Two complementary arms, all objects schema-valid:
- **Rückert (blood, GSE197008):** 29,536 NK across 3 ASAP-seq samples / 6 unique donors; bright/dim
  gated directly on CD56/CD16 surface protein.
- **ReDeeM (marrow, GSE219014):** 4,189 annotated NK across 4 donors; paired RNA labels bright/dim;
  progenitor-proximal complement.

**The core test (Step 5).** Per donor, `frac_mixed_clones` = fraction of testable clones containing
*both* a bright and a dim NK, compared to a within-donor label-shuffle null (1000×) that fixes
clone structure and the bright base-rate. Observed sits on the null in **all 12 donor-units and
both pooled arms**:

| Arm | observed | one-lineage null (mean [95% CI]) | two-lineage floor | p |
|-----|:--------:|:--------------------------------:|:-----------------:|:--:|
| ReDeeM (marrow, RNA) | 0.789 | 0.822 [0.711, 0.921] | 0.000 | ns |
| Rückert (blood, protein) | 0.180 | 0.179 [0.177, 0.182] | 0.000 | ns |

The two arms differ in absolute value by clone-definition design (ReDeeM's few large lineages
nearly all span both subsets; Rückert's many size-2 groups in ~90%-dim blood mix less), but **each
observed value lands on its own null.** A two-lineage world drives `frac_mixed → 0.000` — below
both observed values and both null CIs — so the test could detect segregation and does not.

**The converse (Step 6).** Taking the clonal partition label-free and asking whether the major
lineages carry any bright/dim transcriptional signal: clonal lineage explains a **negligible**
share of the bright/dim axis (max η² = 0.022; the two largest ReDeeM lineages differ by Cohen's
d = 0.16). Where Rückert has large-enough clones (mtASAP2) η² equals its null. There are not "two
clonal lineages mapping to bright and dim."

**Robustness (Step 7).** The continuum holds across binarization thresholds (0.03–0.10), coverage
floors (3–10×), and both clone callers. Mixing is not doublet-driven (mixed:pure variants/cell
ratio 0.95–1.04 in Rückert; dropping the top-5% variant-load cells leaves observed unchanged,
0.857 → 0.857).

## Interpretation and limits

Somatic mtDNA clones are **shared** between CD56-bright and CD56-dim NK — positive evidence for a
common progenitor and a bright↔dim continuum, consistent with the classic bright→dim maturation
model rather than two independent origins. Per the evidence asymmetry, shared clones are the strong
direction of inference, and the power check places this result on that strong side.

Honest limits: ReDeeM marrow under-samples CD56-dim and covers only 4 donors (Old1 has 0 annotated
NK); Rückert is blood (~90% dim), so bright is the limiting compartment. The **ENKP/ILCP
transcriptional origin** (2024 Nat Immunol) was not tested — the exported marker panels lack those
progenitor signatures; connecting clonal lineage to ENKP vs ILCP would need a full-transcriptome
re-export from the Seurat objects (a clean follow-up, not required for the one-vs-two verdict).

## Artifacts

| Step | Figure / table | Memo |
|------|----------------|------|
| 3 | `plan2_step3_power.png`, `plan2_step3_classification.png`, `plan2_step3_counts.csv` | `plan2_step3_power.md` |
| 4 | `plan2_step4_clone_qc.png`, `plan2_step4_clone_qc.csv` | `plan2_step4_clones.md` |
| 5 | `plan2_step5_sharing.png`, `plan2_step5_{perdonor,pooled}.csv` | `plan2_step5_sharing.md` |
| 6 | `plan2_step6_clonal_partition.png`, `plan2_step6_eta2.csv` | `plan2_step6_clonal_partition.md` |
| 7 | `plan2_step7_robustness.png`, `plan2_step7_{sweep,doublet_check}.csv` | `plan2_step7_robustness.md` |

Ingestion summary: `plan2_ingestion_summary.csv`. Dataset selection + gap resolution:
`plan2_dataset_selection.md`. Analysis scripts (ordered, documented):
`scripts/analysis/development_brightdim/` (README there); ingestion adapters in `scripts/mtclone/`.
All scripts committed to `Eomesodermin/Thymic_NK_development` (through commit `f287563`).
