.PHONY: all data immgliph reference metrics figures report clean

all: report

data: data/glanville2017/cdr3_input.tsv data/huang2020/cdr3_input.tsv

data/glanville2017/cdr3_input.tsv data/glanville2017/antigen_labels.tsv: scripts/01_prep_data.R
	Rscript scripts/01_prep_data.R --dataset glanville

data/huang2020/cdr3_input.tsv data/huang2020/antigen_labels.tsv: scripts/01_prep_data.R
	Rscript scripts/01_prep_data.R --dataset huang

immgliph: data
	Rscript scripts/02_run_immgliph.R --dataset glanville
	Rscript scripts/02_run_immgliph.R --dataset huang

reference: data
	bash scripts/03_run_reference.sh

results/metrics.tsv: immgliph reference scripts/04_compute_metrics.R
	Rscript scripts/04_compute_metrics.R

metrics: results/metrics.tsv

results/figures/concordance_summary.png: results/metrics.tsv scripts/05_make_figures.R
	Rscript scripts/05_make_figures.R

figures: results/figures/concordance_summary.png

report.html: figures report.qmd
	quarto render report.qmd

report: report.html

clean:
	rm -rf results/raw/* results/figures/*.png results/metrics.tsv report.html
