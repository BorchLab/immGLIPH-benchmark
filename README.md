# immGLIPH-benchmark

Concordance benchmark of [immGLIPH](https://github.com/BorchLab/immGLIPH) against the published
cluster vectors from the original [GLIPH](https://doi.org/10.1038/nature22976)
(Glanville et al., 2017) and [GLIPH2](https://doi.org/10.1038/s41587-020-0505-4)
(Huang et al., 2020) papers.

## Why this exists

immGLIPH is an R reimplementation of GLIPH and GLIPH2. This repo answers, with numbers, the
question: *does immGLIPH reproduce the original implementations' published cluster output?*

## Headline result

When `runGLIPH()` is invoked with each paper's documented parameters and (for Huang) the
paper's post-hoc filter, immGLIPH reproduces the published cluster vectors at high concordance
on the intersection of CDR3s present in both runs:

| Dataset | immGLIPH configuration | n | ARI | NMI | Pairwise F1 | Precision | Recall |
|:--------|:-----------------------|--:|----:|----:|------------:|----------:|-------:|
| Glanville 2017 | `gliph1` + paper params | 144 | **0.985** | 0.994 | **0.985** | 1.000 | 0.971 |
| Huang 2020 | `gliph2` + paper params + filter | 171 | **0.863** | 0.968 | **0.867** | 0.931 | 0.812 |

![](results/figures/concordance_summary.png)

Full numbers (including default-parameter and full-input-universe variants) are in
[`results/metrics.tsv`](results/metrics.tsv); the rendered narrative including methodology,
caveats, and reproducibility steps is in [`report.qmd`](report.qmd).

## Approach

Rather than re-running the original GLIPH Perl scripts and the GLIPH2 binary inside Docker
containers, this benchmark compares immGLIPH against the **published reference cluster vectors**
that each paper shipped as supplementary data:

- **Glanville:** MOESM6 sheet `all GLIPH Group Scoring` — 43 published GLIPH groups, 172 unique CDR3s.
- **Huang:** MOESM5 sheet `Filtered_list(354)` joined to `GLIPH_group_member` — 354 published
  GLIPH2 groups, 1,263 unique CDR3s. (The unfiltered 4,185-cluster set is also extracted for
  completeness.)

This is the canonical "what GLIPH/GLIPH2 said on this data" reference and avoids the
Bioconductor-incompatible burden of shipping reference binaries.

## Datasets

| Dataset | Input source | Reference cluster source |
|---|---|---|
| Glanville 2017 | MOESM4 sheet `Sheet1` (single-cell paired chains, 5,661 CDR3b) | MOESM6 sheet `all GLIPH Group Scoring` |
| Huang 2020 | MOESM3 sheet `bulk TCR` (10,501 CDR3b) | MOESM5 sheet `Filtered_list(354)` + `GLIPH_group_member` |

Provenance, SHA256 hashes, and per-sheet schema notes are in
[`data/glanville2017/SOURCE.md`](data/glanville2017/SOURCE.md) and
[`data/huang2020/SOURCE.md`](data/huang2020/SOURCE.md).

## Pipeline

```
01_prep_data.R           Normalize each paper's TCR table to a common schema.
02_extract_reference.R   Extract published cluster vectors (raw + curated).
03_run_immgliph.R        Run runGLIPH() with --method and --params {default,paper}.
03b_filter_immgliph.R    Apply each paper's documented post-hoc filter.
04_compute_metrics.R     Compute ARI, NMI, pairwise precision/recall/F1.
05_make_figures.R        Render headline figure plus auxiliary distributions.
report.qmd               Narrative report with metrics tables and figures.
```

## Metrics

For each (immGLIPH, reference) pair, two universes are reported:

- **`intersection`** — CDR3s in both the immGLIPH input and the published reference cluster
  output. Cleanest answer to "do the tools cluster the same objects the same way".
- **`input_full`** — every CDR3b in either input or reference (singletons fill the gaps).
  Penalises over- and under-clustering relative to the reference.

Cluster-pair metrics:

- **Adjusted Rand Index (ARI)** — chance-corrected agreement between cluster assignments.
- **Normalized Mutual Information (NMI)** — information-theoretic agreement.
- **Pairwise precision / recall / F1** — over all CDR3 pairs, treating "co-clustered by reference"
  as truth.

## Running the benchmark

```bash
# one-time
Rscript -e 'renv::restore()'

# place the supplementary XLSX files into:
#   data/glanville2017/raw/   data/huang2020/raw/
# (see SOURCE.md in each for filenames and SHA256s)

# end-to-end
make all
```

Artifacts produced: `results/metrics.tsv`, `results/figures/*.png`, `report.html`.

## Citation

If you use these results, please cite immGLIPH and the original GLIPH and GLIPH2 papers.

## License

MIT — see [LICENSE](LICENSE).
