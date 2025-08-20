# =============================================================================
# D-statistics (ABBA–BABA) with Block Jackknife
# -----------------------------------------------------------------------------
# This script:
#   1) Loads allele-frequency data (one row per site).
#   2) Defines a D-statistic function.
#   3) Builds genomic blocks for jackknife (via helper functions you source).
#   4) Computes D and a jackknife Z-score across multiple P1 populations.
#   5) Saves a CSV and a PDF plot.
# =============================================================================

# ==========================
# 1) Load required packages
# ==========================
library(ggplot2)
library(viridis)

# ===========================================
# 2) Output paths and file name conventions
# ===========================================
run_date <- format(Sys.Date(), "%d%b%Y")
out_dir  <- "results_dstats"
dir.create(out_dir, showWarnings = FALSE)

csv_A <- file.path(out_dir, paste0("D_stat_results_setA_", run_date, ".csv"))
pdf_A <- file.path(out_dir, paste0("D_stat_plot_setA_", run_date, ".pdf"))

csv_file_name <- csv_A
pdf_file_name <- pdf_A

# ======================
# 3) Load input dataset
# ======================
# Expecting columns:
#   - scaffold (chr/scaffold ID)
#   - position (genomic coordinate)
#   - one numeric column per population (allele freq in [0,1])
freq_table <- read.table(
  "Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.150ind_amo_populations.tsv.gz",
  header = TRUE,
  as.is  = TRUE
)

# ===========================================
# 4) Define D-statistic (ABBA–BABA) function
# ===========================================
# Inputs: three numeric vectors p1, p2, p3 (allele frequencies)
# Output: scalar D value
D.stat <- function(p1, p2, p3) {
  ABBA <- (1 - p1) * p2 * p3
  BABA <- p1 * (1 - p2) * p3
  (sum(ABBA, na.rm = TRUE) - sum(BABA, na.rm = TRUE)) /
    (sum(ABBA, na.rm = TRUE) + sum(BABA, na.rm = TRUE))
}

# ===============================================
# 5) Load jackknife helpers
# ===============================================
# Expects functions:
#   - get.block.indices(block_size, positions, chromosomes)
#   - block.jackknife(block_indices, FUN, p1, p2, p3)
source("/home/nh253642e/sandbox/genomicsGeneral/genomics_general-master/jackknife.R")

# ==================================================
# 6) Build genomic blocks (for block-jackknife SEs)
# ==================================================
block_indices <- get.block.indices(
  block_size  = 1e6,               # 1 Mb blocks
  positions   = freq_table$position,
  chromosomes = freq_table$scaffold
)
cat("Genome divided into", length(block_indices), "blocks.\n")


# ====================================================
# 7) Set population labels (used in the loop below)
# ====================================================
P1 <- "cra_A"
P2 <- "ruficaudus"
P3 <- "amoenus"

# ===========================
# 8) Compute D across P1 set
# ===========================
# Loop over P1 = "cra_1"..."cra_20" while keeping (current) P2, P3 fixed.
results <- data.frame(
  P1                 = character(),
  D                  = numeric(),
  D_jackknife_mean   = numeric(),
  D_Z                = numeric(),
  stringsAsFactors   = FALSE
)

for (i in 1:20) {
  P1 <- paste0("cra_", i)
  cat("\n=== Analyzing", P1, "vs", P2, "and", P3, "===\n")

  # -- Compute D on full data
  D <- D.stat(freq_table[[P1]], freq_table[[P2]], freq_table[[P3]])

  # -- Jackknife estimate (leave-one-block-out)
  D_jackknife <- block.jackknife(
    block_indices = block_indices,
    FUN           = D.stat,
    freq_table[[P1]], freq_table[[P2]], freq_table[[P3]]
  )

  D_Z <- D_jackknife$mean / D_jackknife$standard_error

  # -- Append a row to results
  results <- rbind(
    results,
    data.frame(
      P1               = P1,
      D                = round(D, 4),
      D_jackknife_mean = round(D_jackknife$mean, 4),
      D_Z              = round(D_Z, 3)
    )
  )
}


# ====================
# 9) Save results CSV
# ====================
write.csv(results, csv_file_name, row.names = FALSE)
cat("Results saved to", csv_file_name, "\n")


# ====================
# 10) Plotting section
# ====================
# Read the CSV we just wrote, relabel columns for plotting, and render a PDF.
data <- read.csv(csv_file_name)
colnames(data)[colnames(data) == "P1"] <- "SAMPLE_ID"
colnames(data)[colnames(data) == "D"]  <- "D_stat"

pdf(pdf_file_name, width = 10, height = 6)

ggplot(data, aes(SAMPLE_ID, D_stat, color = D_stat)) +
  geom_point(size = 5) +
  scale_color_viridis(
    option    = "plasma",   # palette; other options include "D", "C", "A", "B", "E"
    direction = -1,         # reverse scale (darker for smaller values)
    limits    = c(-0.32, -0.20),
    breaks    = seq(-0.32, -0.20, by = 0.04),
    name      = "D-stat"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(y = "D-statistic")

dev.off()
cat("Plot saved to", pdf_file_name, "\n")
