# Code for "Near-tight values of the Lonely Runner spectrum for five and six speeds"

Everything is exact (integer or rational arithmetic; no floating point in any
decision).  Python 3 (standard library; `numpy` optional, only to speed up the
witness filter in `classify6.py`) and a C99 compiler with OpenMP.

    make              # compiles the C programs
    ./run_all.sh      # reproduces every number in the paper (about 20 minutes)

| paper | program | what it does | time |
|---|---|---|---|
| Theorem 3.7 | `classify6.py` (uses `lr_ml.py`, `lr_subtorus.py`, `lr_classify.py`, `lr_lemma_s5.py`, `lr_lemma_s6.py`) | the three critical subtori of (R/Z)^6 | 1 min |
| Theorem 3.8 | `lr5.py` | the two critical subtori of (R/Z)^5 | 10 s |
| Theorem 4.2 | `fastrunner.c`, `fastrunner5.c` | exact ML of V u {c} for c up to 30000 / 20000 | 1 min |
| Proposition 4.1 | `prejump.c` | numerical check of the pre-jump proposition | 1 s |
| Section 5 / Theorem 5.6 | `uc_theorem.py` (uses `jk_uc.py`, `jk_uc2.py`, `jk_uc3.py`, `jk_uc4.py`, `ucpoints.c`) | the complete analysis of U_C: components, thresholds, period, cones, strips, direct checks, ratio | 2.5 min |
| Remarks after Theorem 5.6 | `ucdump.c` + `python3 jk_uc2.py 300`; `ucscan.c` | consistency of (5.2) with exact ML up to 300; direct scan of U_C | 2 min |
| Section 6 | `survey.c` | exhaustive searches for n = 4..8 with membership in the critical subtori | 15 min |
| everywhere | `lrk.h`, `mltest.c` | the exact ML kernel in C (Lemma 2.3) | |

`uc_theorem_output.txt` is the output of `uc_theorem.py` on the author's machine.

Archived at Zenodo: https://doi.org/10.5281/zenodo.XXXXXXX (replace with the DOI of the record).
License: MIT (see LICENSE). Cite via CITATION.cff.
