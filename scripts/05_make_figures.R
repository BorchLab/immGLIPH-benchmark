#!/usr/bin/env Rscript
# Build static PNGs from results/metrics.tsv and the assignment TSVs.
#
# Outputs:
#   results/figures/concordance_summary.png  — bar chart, ARI/NMI/F1 by tool pair
#   results/figures/cluster_size_distribution.png  — overlaid density per dataset
#   results/figures/cluster_count_comparison.png   — # clusters and # clustered CDR3s

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(fs)
})

theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank()))

# ---- 1. Concordance summary ------------------------------------------------

m <- read_tsv("results/metrics.tsv", show_col_types = FALSE) %>%
  filter(metric %in% c("ARI", "NMI", "pairwise_F1"))

p1 <- ggplot(m,
             aes(x = metric, y = value, fill = universe)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", value)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3) +
  facet_wrap(~ tool_pair, ncol = 2, scales = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  labs(
    title = "immGLIPH vs. published reference cluster vectors",
    subtitle = "Higher = stronger agreement with original GLIPH/GLIPH2 output",
    x = NULL, y = NULL, fill = "Comparison universe"
  )

dir_create("results/figures")
ggsave("results/figures/concordance_summary.png", p1,
       width = 8.5, height = 4.5, dpi = 300)

# ---- 2. Cluster-size distribution -----------------------------------------

read_assign <- function(path, tool, dataset) {
  read_tsv(path, show_col_types = FALSE) %>%
    mutate(tool = tool, dataset = dataset)
}

assigns <- bind_rows(
  read_assign("results/raw/immgliph_glanville_assignments.tsv",  "immGLIPH",   "glanville"),
  read_assign("results/raw/reference_gliph_glanville.tsv",        "GLIPH",      "glanville"),
  read_assign("results/raw/immgliph_huang_assignments.tsv",       "immGLIPH",   "huang"),
  read_assign("results/raw/reference_gliph2_huang.tsv",           "GLIPH2",     "huang")
)

sizes <- assigns %>%
  count(dataset, tool, cluster_id, name = "size") %>%
  filter(size >= 2)   # singletons elide the visual signal

p2 <- ggplot(sizes, aes(x = size, fill = tool)) +
  geom_histogram(position = "identity", alpha = 0.55,
                 binwidth = 1, boundary = 1) +
  facet_wrap(~ dataset, scales = "free", ncol = 2) +
  scale_x_continuous(limits = c(1, 30)) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Cluster size distribution (multi-member clusters only)",
    x = "Members per cluster", y = "Number of clusters", fill = NULL
  )

ggsave("results/figures/cluster_size_distribution.png", p2,
       width = 8.5, height = 4.5, dpi = 300)

# ---- 3. Cluster-count comparison ------------------------------------------

counts <- assigns %>%
  group_by(dataset, tool) %>%
  summarise(
    n_clusters = n_distinct(cluster_id),
    n_cdr3_clustered = n_distinct(CDR3b),
    .groups = "drop"
  ) %>%
  pivot_longer(c(n_clusters, n_cdr3_clustered),
               names_to = "stat", values_to = "value") %>%
  mutate(stat = recode(stat,
                       n_clusters = "# clusters",
                       n_cdr3_clustered = "# CDR3s clustered"))

p3 <- ggplot(counts, aes(x = tool, y = value, fill = tool)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = scales::comma(value)), vjust = -0.35, size = 3) +
  facet_grid(stat ~ dataset, scales = "free_y") +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Cluster and member counts", x = NULL, y = NULL)

ggsave("results/figures/cluster_count_comparison.png", p3,
       width = 8.5, height = 5, dpi = 300)

cat("Figures written to results/figures/\n")
