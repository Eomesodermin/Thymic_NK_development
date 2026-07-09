"""recall/ — OPTIONAL consistency re-calling module (DEFERRED, stub).

Purpose: when we later want both datasets called by identical software (to remove the
"different authors, different thresholds" confound), re-call mtDNA variants from raw reads:
  - mgatk   for mtscATAC (needs cellranger-atac possorted BAMs);
  - redeemR  for ReDeeM (needs the raw-genotype tables).

This is NOT on the critical path for a first result — the deposited per-cell calls are the
primary input (see the top-level README). Do not build this until a question actually needs
identical-software calling.

FIRST TASK when this is built (do this before writing any wrapper):
    check_raw_availability(dataset) -> report whether the required raw inputs exist and are
    accessible. mgatk BAMs are usually in SRA/EGA and may be controlled-access; the ReDeeM
    raw-genotype tables' presence for the *human* atlas (GSE219014) is UNCONFIRMED — GSE219014
    ships only the called matrix, RNA h5, and ATAC fragments.
"""

def check_raw_availability(dataset: str):
    raise NotImplementedError(
        "recall/ is a deferred stub. Its first task is to verify raw-input availability "
        "(mgatk BAMs via SRA/EGA; redeemR raw-genotype tables) before any re-calling wrapper "
        "is written. See recall/__init__.py docstring and the top-level README."
    )
