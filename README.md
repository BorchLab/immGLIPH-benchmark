# immGLIPH-benchmark

Concordance benchmark of [immGLIPH](https://github.com/BorchLab/immGLIPH) against the original
[GLIPH](https://doi.org/10.1038/nature22976) (Glanville et al., 2017) and
[GLIPH2](https://doi.org/10.1038/s41587-020-0505-4) (Huang et al., 2020) implementations.

## Why this exists

immGLIPH is an R reimplementation of GLIPH and GLIPH2. This repo answers, with numbers, the
question: *does immGLIPH produce results consistent with the reference implementations?*

## Datasets

| Dataset | Source | Used for |
|---|---|---|
| Glanville 2017 | Nature Supp. Tables 1 and 3 | immGLIPH (`method = "gliph1"`) vs. GLIPH Perl |
| Huang 2020 | Nat. Biotechnol. Supp. Tables | immGLIPH (`method = "gliph2"`) vs. GLIPH2 standalone |

Provenance for each dataset is recorded in `data/<dataset>/SOURCE.md`.

## Reference implementations

Both reference tools run inside Docker containers (`linux/amd64` platform; required for the
GLIPH2 binary on Apple Silicon).

- **GLIPH** — pinned commit of the upstream Perl scripts; see `reference/gliph/install.sh`.
- **GLIPH2** — pinned binary download with verified SHA256; see `reference/gliph2/install.sh`.

## Metrics

For each (immGLIPH, reference) pair on each dataset:

1. **Adjusted Rand Index** and **Normalized Mutual Information** between cluster assignments.
2. **Pairwise concordance** — precision, recall, F1 over all CDR3 pairs.
3. **Antigen-label recovery** — fraction of antigens whose labeled members are recovered as a
   purity ≥ 0.7, Fisher *p* < 0.05 cluster.

Definitions and exact formulas are in `report.qmd`.

## Running the benchmark

```bash
# one-time
Rscript -e 'renv::restore()'
bash reference/gliph/install.sh
bash reference/gliph2/install.sh

# end-to-end
make all
```

Artifacts produced: `results/metrics.tsv`, `results/figures/*.png`, `report.html`.

## Results

*Populated after the first run. See `results/metrics.tsv` and `report.html`.*

## Citation

If you use these results, please cite immGLIPH and the original GLIPH/GLIPH2 papers.

## License

MIT — see [LICENSE](LICENSE).
