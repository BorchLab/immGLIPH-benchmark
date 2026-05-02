.PHONY: all data immgliph reference metrics figures report clean

all: report

# ---- 1. Data prep --------------------------------------------------------

data: data/glanville2017/cdr3_input.tsv data/huang2020/cdr3_input.tsv

data/glanville2017/cdr3_input.tsv data/glanville2017/antigen_labels.tsv: scripts/01_prep_data.R data/glanville2017/raw
	Rscript scripts/01_prep_data.R --dataset glanville

data/huang2020/cdr3_input.tsv data/huang2020/antigen_labels.tsv: scripts/01_prep_data.R data/huang2020/raw
	Rscript scripts/01_prep_data.R --dataset huang

# ---- 2. Reference cluster extraction (Approach B) -----------------------

reference: results/raw/reference_gliph_glanville.tsv results/raw/reference_gliph2_huang.tsv

results/raw/reference_gliph_glanville.tsv results/raw/reference_gliph2_huang.tsv: scripts/02_extract_reference.R
	Rscript scripts/02_extract_reference.R

# ---- 3. immGLIPH runs ----------------------------------------------------

immgliph: results/raw/immgliph_glanville_assignments.tsv results/raw/immgliph_huang_assignments.tsv

results/raw/immgliph_glanville_assignments.tsv: scripts/03_run_immgliph.R data/glanville2017/cdr3_input.tsv
	Rscript scripts/03_run_immgliph.R --dataset glanville

results/raw/immgliph_huang_assignments.tsv: scripts/03_run_immgliph.R data/huang2020/cdr3_input.tsv
	Rscript scripts/03_run_immgliph.R --dataset huang

# ---- 4. Metrics ----------------------------------------------------------

metrics: results/metrics.tsv

results/metrics.tsv: scripts/04_compute_metrics.R immgliph reference
	Rscript scripts/04_compute_metrics.R

# ---- 5. Figures ----------------------------------------------------------

figures: results/figures/concordance_summary.png

results/figures/concordance_summary.png: scripts/05_make_figures.R results/metrics.tsv
	Rscript scripts/05_make_figures.R

# ---- 6. Report -----------------------------------------------------------

report: report.html

report.html: report.qmd figures
	quarto render report.qmd

# ---- Cleanup -------------------------------------------------------------

clean:
	rm -rf results/raw/immgliph_*.tsv results/raw/immgliph_*.rds \
	       results/figures/*.png results/metrics.tsv report.html

clean-all: clean
	rm -rf results/raw/reference_*.tsv \
	       data/glanville2017/cdr3_input.tsv data/glanville2017/antigen_labels.tsv \
	       data/huang2020/cdr3_input.tsv data/huang2020/antigen_labels.tsv
