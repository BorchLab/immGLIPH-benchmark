# Huang 2020 — supplementary data provenance

**Citation:** Huang H, Wang C, Rubelt F, Scriba TJ, Davis MM. Analyzing the *Mycobacterium
tuberculosis* immune response by T-cell receptor clustering with GLIPH2 and genome-wide
antigen screening. *Nat Biotechnol* **38**, 1194–1202 (2020).
**DOI:** [10.1038/s41587-020-0505-4](https://doi.org/10.1038/s41587-020-0505-4)
**Retrieved:** 2026-05-02
**Source:** Springer Nature supplementary information for the paper above.

## Files

| Filename | SHA256 | Used by benchmark? |
|---|---|---|
| `41587_2020_505_MOESM3_ESM.xlsx` | `e5802812c0e3dc1344acc85ea2b02ced92a1870f81a92c4e7cd0bb3b8392a005` | yes — sheet `bulk TCR` is the input |
| `41587_2020_505_MOESM4_ESM.xlsx` | `8e3d2afe7ae53a227e3587c3ab9d124d60d06cac986e5213c9bf57bd3d3ffb4e` | yes — sheet `known_CDR3` enriches antigen labels |
| `41587_2020_505_MOESM5_ESM.xlsx` | `19f0baae842e2c1c9df8c4ec677e48f821a01fd9a4fd249dc3c7476788d050e7` | yes — published GLIPH2 cluster output for direct comparison |
| `41587_2020_505_MOESM6_ESM.xlsx` | `4768ecf12277e62e258dbef09b0a839a87a2535fa1caf6900458681d8a49b1b4` | no (Mtb antigen library annotations) |

## Sheet → schema map

`MOESM3.xlsx`, sheet **`bulk TCR`** (input TCR table):
- Columns: `CDR3b`, `Vb`, `Jb`, `CDR3a`, `Va`, `Ja`, `Individual`, `Counts`.
- Antigen specificity is implicit: all rows are Mtb tetramer-sorted from the TB cohort.

`MOESM4.xlsx`, sheet **`known_CDR3`** (epitope-labeled subset):
- Columns: `CDR3b`, `TRBV`, `TRBJ`, `peptide` (no header in source; assigned by the script).
- Used to upgrade `antigen = "Mtb"` to `antigen = "Mtb (known epitope)"` for matching CDR3s.

`MOESM5.xlsx`, sheet **`GLIPH_group_member`** (published GLIPH2 output):
- Columns: `index`, `pattern`, `Fisher_score`, `number_individual`, `number_unique_cdr3`,
  `final_score`, `hla_score`, `vb_score`, ... plus per-cluster member CDR3s.
- This is the cluster assignment from the original GLIPH2 binary run. Used as the
  reference cluster vector for the concordance comparison.

## Filters applied

- CDR3b starts with `C`, ends with `F`.
- Length ∈ [8, 30].
- Drop rows missing CDR3b.

## Licensing

Springer Nature supplementary materials are distributed alongside the article
under the journal's standard reuse terms. Files are kept locally only; this
directory is `.gitignore`d.
