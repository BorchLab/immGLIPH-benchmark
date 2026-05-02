# Glanville 2017 — supplementary data provenance

**Citation:** Glanville J, Huang H, Nau A, *et al.* Identifying specificity groups in
the T cell receptor repertoire. *Nature* **547**, 94–98 (2017).
**DOI:** [10.1038/nature22976](https://doi.org/10.1038/nature22976)
**Retrieved:** 2026-05-02
**Source:** Springer Nature supplementary information for the paper above.

## Files

| Filename | SHA256 | Used by benchmark? |
|---|---|---|
| `41586_2017_BFnature22976_MOESM2_ESM.xlsx` | `343075bb63ed878a5ea1952c5365ca587a1a63fe696d22c8168e7b985b9d0689` | no (overlap with published GLIPH groups was 1/172) |
| `41586_2017_BFnature22976_MOESM3_ESM.xlsx` | `5e4d2b4591eb5d8750247701c8ee37abc3b87118b1e36a557ce78a882325d4b2` | no (structural alignment) |
| `41586_2017_BFnature22976_MOESM4_ESM.xlsx` | `3ad6f619a5ed18eb97ad87424a9f55f5ffc976b4005d7768cd084fefb7e8e27a` | yes — sheet `Sheet1` is the input single-cell paired-chain TCR table; 170/172 reference CDR3s present |
| `41586_2017_BFnature22976_MOESM5_ESM.xlsx` | `237a00f94492a2d1bf4f61fd37bab51097950968f6de28f89988ee02b378cdf2` | no (annotated convergence groups) |
| `41586_2017_BFnature22976_MOESM6_ESM.xlsx` | `85593487890d211e9d303ca204b2d89c92cb43f13a7c8cf2adcfc69eb708b3dc` | yes — published GLIPH cluster output for direct comparison |
| `41586_2017_BFnature22976_MOESM7_ESM.xlsx` | `e4ba23bbe8b97e4fe69862a6c4e212ac9c58e468a80f697c2ec2d14841f73f0b` | no (HLA typing) |
| `41586_2017_BFnature22976_MOESM8_ESM.xlsx` | `f4b2550a9bddebaa6f5ae21018d962ce79d8304a6ebb08acf9c0b65ab1c97134` | no (replicates) |

## Sheet → schema map

`MOESM4.xlsx`, sheet **`Sheet1`** (input TCR table):
- Single-cell paired-chain TCR table from tetramer-sorted T cells.
- Columns of interest: `Donor`, `Stim`, `Vbeta`, `Jbeta`, `CDR3beta`, `BetaReads`.
- 170 of 172 published-reference CDR3s are present (98.8%). Selected over MOESM2
  `Raw` (which had only 1/172 overlap) because the published GLIPH groups in
  MOESM6 were derived from this single-cell dataset.

`MOESM6.xlsx`, sheet **`all GLIPH Group Scoring`** (published GLIPH output):
- Columns of interest: `Group` (Roman numeral cluster ID), `Name`,
  `CDR3bs, all clones` (space-delimited member CDR3 list).
- 51 rows total; 8 are footer/criteria rows, leaving 43 actual GLIPH groups
  spanning 172 unique CDR3s.
- This is the cluster assignment from the original GLIPH Perl run. Used as the
  reference cluster vector for the concordance comparison.

## Filters applied (matching paper Methods)

- CDR3b starts with `C`, ends with `F`.
- Length ∈ [8, 30].
- Drop rows missing CDR3b.

## Licensing

Springer Nature supplementary materials are distributed alongside the article
under the journal's standard reuse terms. Files are kept locally only; this
directory is `.gitignore`d.
