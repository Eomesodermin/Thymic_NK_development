# NK↔T shared-clone origin test on ReDeeM marrow (GSE219014) — feasibility probe

**Question.** Do NK cells and T cells in the same donor carry the same somatic-mtDNA
clone? A shared clone (a rare mtDNA variant / ClonalGroup carried by both an NK and a T
cell) implies a common ancestor upstream of the NK-vs-T fate decision — a direct test of
shared developmental origin, and the analysis mtDNA can do that the abandoned TCR-footprint
approach cannot.

**Data.** ReDeeM deep-mtDNA scATAC BMMC, 4 donors with NK present (Old1 skipped: 0 NK).
Per-cell binarized heteroplasmy (Liu threshold 0.07), informative variants = pseudobulk
heteroplasmy in (0, 0.01] and detected in ≥5 cells. Coverage is real depth-file data
(median 12–25×). Clones NEVER pooled across donors.

## 1. Power gate — GO
NK cells carrying ≥1 informative mtDNA variant, per donor:
Old2 219/235 (93%), Young1.T1 1192/1346 (89%), Young1.T2 1403/1591 (88%), Young2 902/1017 (89%).
T cells: 63–92%. Median informative variants/cell = 3–4 for NK. **Coverage is not the binding
constraint** — the data are rich enough per cell to attempt the test.

## 2–3. The money numbers — NK↔T shared clones (author ClonalGroup)
| donor | clones w/ NK | NK–T shared | frac of NK-clones | vs null (mean) | p(obs≥null) |
|-------|-------------:|------------:|------------------:|---------------:|------------:|
| Old2      | 43 | 43 | 1.00 | 52.6 | 0.997 |
| Young1.T1 | 76 | 76 | 1.00 | 77.7 | 1.000 |
| Young1.T2 | 78 | 78 | 1.00 | 77.9 | 0.922 |
| Young2    | 38 | 38 | 1.00 | 37.9 | 0.923 |

**Every** NK-containing author clone also contains T cells (frac = 1.00), and the observed
count sits **at or below** the within-donor label-shuffle null (p for obs≥null ≈ 0.92–1.00).
This is the signature of clones that are too coarse to be lineage-informative: author
ClonalGroups have median size 46–150 cells (per-donor largest 261–1562), so nearly every clone
trivially spans multiple mature lineages.

## 4. Null / chance & rare-variant re-call
Re-calling fine clones from rare private variants (mtclone `variant_group`, median clone
size 2) does **not** rescue signal: NK–T shared = 70/467/571/418, all sitting **at the
permutation null** (p = 0.39–0.98). No donor shows above-chance NK–T co-membership under
either clone definition.

## 5. Variant-level test and the decisive control
Counting individual **private** variants (≤1% of cells, ≤20 detected) co-carried by an NK
and a T cell gives an apparently above-chance signal vs an all-cell-type shuffle
(159/977/1378/1197, p = 0.002 each donor). **This is a coverage artifact, not lineage.**
The all-cell-type null mixes high-coverage lymphocytes with low-coverage myeloid/erythroid
cells, so any two high-quality subsets co-detect rare variants above that null.

The lineage-specific control settles it. Private-variant sharing with NK, normalized per
1000 partner cells:
| donor | NK–T | NK–B | NK–Mono |
|-------|-----:|-----:|--------:|
| Old2      | 78.6  | 87.4  | 80.8  |
| Young1.T1 | 289.5 | 277.3 | 401.3 |
| Young1.T2 | 327.4 | 322.0 | 487.2 |
| Young2    | 265.6 | 251.3 | 312.0 |

NK shares private variants with B cells (a different lymphoid lineage) at essentially the
same rate as with T cells, and with **monocytes (a myeloid lineage) at a HIGHER rate**. If
shared private variants reflected NK↔T common ancestry, NK–T would exceed NK–B and NK–Mono.
It does not. The apparent sharing tracks partner-cell coverage/quality, not developmental
relatedness.

## 6. Progenitor rooting (Old2, the only object with HSC)
Of 43 author NK–T shared clones, 39 also contain HSC and 42 contain GMP. That a
granulocyte–monocyte progenitor sits in almost every "shared clone" confirms the author
clones are too coarse to root anything — they are not resolving a lymphoid common ancestor.
No CLP/MPP/LMPP are annotated in any object, so a clean lymphoid-progenitor rooting is not
available in this data regardless.

## Verdict: (c) NOT viable with these mtscATAC snapshots — with one caveat

The NK↔T shared-clone origin test is **not feasible** on the ReDeeM marrow data as it
stands, and the reason is structural rather than a matter of more cells:

- **Binding constraint is clone resolution, not NK cell number or coverage.** Power is fine
  (88–93% of NK informative). But author clones are too large to be lineage-specific
  (everything shares), and rare-variant clones are too sparse to exceed chance.
- **The one signal that looks positive (private-variant co-detection) fails its lineage
  control**: NK shares with B and myeloid cells as much as or more than with T, so it is
  coverage-driven co-detection, not common ancestry.
- **Mechanistic reason it may be unrecoverable even with more cells:** the NK-vs-T (indeed
  lymphoid-vs-myeloid) common ancestor is an HSC/early-progenitor that divided early in
  ontogeny. A clone-defining somatic mtDNA mutation must have arisen *in that specific
  ancestor and been inherited by both an NK and a T descendant that both survive to the
  marrow snapshot*. In an adult marrow cross-section that lineage bottleneck is deep and
  most such lineages are unsampled — the shared mutation is either too old (near-public,
  uninformative) or private to one branch. mtscATAC snapshots see the leaves, not the split.

### Recommendation
1. **Do not build a large marrow atlas to chase NK↔T marrow sharing** — the test is
   confounded by coverage and limited by ancestral depth; more of the same data will not
   move it off the null.
2. **If the origin question is pursued via mtDNA, change the substrate, not the scale.**
   The clean version needs (a) a captured lymphoid progenitor compartment (CLP/ETP/thymic
   DN) sequenced *in the same individual* as the mature NK and T, so the shared mutation
   can be caught at the split rather than inferred between two leaves; and (b) a
   coverage-matched cross-lineage control (NK–B / NK–myeloid) reported alongside every
   NK–T number — which this probe shows is mandatory to avoid a co-detection false positive.
3. **The permutation-null-plus-lineage-control design here is reusable** and should gate any
   future mtDNA origin claim: an NK–T sharing count is only evidence if it exceeds both the
   label-shuffle null AND the NK–B / NK–myeloid controls.

*Clean negative: the honest read is that adult-marrow mtscATAC snapshots cannot certify an
NK–T common ancestor, because the ancestor is too deep to have left a shared, private,
coverage-robust mark on surviving leaves of both lineages.*
