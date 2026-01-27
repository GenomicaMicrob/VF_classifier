#!/usr/bin/env Rscript

# plot_vf_heatmap.R
# Usage: Rscript plot_vf_heatmap.R [input.csv] [output.pdf]

# --- Dependency Management ---
message("Checking dependencies...")

# 1. Check if we are in a Conda environment
in_conda <- Sys.getenv("CONDA_PREFIX") != ""

if (in_conda) {
  # FORCE R to prioritize the conda library and potentially ignore system libraries
  # that lead to "undefined symbol" or "installed before R 4.0.0" errors.
  conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")
  if (dir.exists(conda_lib)) {
    # Keep only the conda library and the base R library to avoid system pollution
    current_paths <- .libPaths()
    # Find the base R library (usually ends in /R/library)
    base_lib <- current_paths[grep("R/library$", current_paths)]
    # Set the path strictly: [Conda Lib, Base R Lib]
    .libPaths(unique(c(conda_lib, base_lib)))
  }
}

# 2. Define required packages
packages <- c("RColorBrewer", "dplyr", "tidyr", "GetoptLong", "BiocManager", "ComplexHeatmap")
missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]

if (length(missing_packages) > 0) {
  if (in_conda) {
    message("---------------------------------------------------------")
    message("MISSING PACKAGES DETECTED IN CONDA ENVIRONMENT:")
    message(paste(missing_packages, collapse = ", "))
    message("\nPlease update your conda environment using:")
    message("  conda env update -f environment.yml")
    message("---------------------------------------------------------")
    stop("Dependencies missing from Conda environment. Execution halted.")
  } else {
    # Non-conda logic: Try to find a writable library
    get_writable_lib <- function() {
      for (path in .libPaths()) {
        if (file.access(path, 2) == 0) {
          return(path)
        }
      }
      user_lib <- Sys.getenv("R_LIBS_USER")
      if (user_lib == "" || user_lib == "NA") user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
      if (!dir.exists(user_lib)) try(dir.create(user_lib, recursive = TRUE), silent = TRUE)
      if (dir.exists(user_lib) && file.access(user_lib, 2) == 0) {
        if (!(user_lib %in% .libPaths())) .libPaths(c(user_lib, .libPaths()))
        return(user_lib)
      }
      return(NULL)
    }

    target_lib <- get_writable_lib()
    message(paste("Installing missing packages to:", ifelse(is.null(target_lib), "Default", target_lib)))

    # Install CRAN packages
    cran_packages <- intersect(missing_packages, c("RColorBrewer", "dplyr", "tidyr", "GetoptLong", "BiocManager"))
    if (length(cran_packages) > 0) {
      install.packages(cran_packages, repos = "http://cran.us.r-project.org", lib = target_lib)
    }

    # Install ComplexHeatmap if still missing
    if ("ComplexHeatmap" %in% missing_packages) {
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "http://cran.us.r-project.org", lib = target_lib)
      }
      BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE, lib = target_lib)
    }
  }
}

# --- Dependency Loading & Error Handling ---
load_packages <- function() {
  # We try to load the packages and catch the specific "R 4.0.0" version error
  tryCatch(
    {
      suppressPackageStartupMessages({
        library(ComplexHeatmap)
        library(RColorBrewer)
        library(dplyr)
        library(tidyr)
      })
    },
    error = function(e) {
      message("\n[ERROR] Failed to load R dependencies.")
      message(paste("Original Error:", e$message))

      if (grepl("installed before R 4.0.0", e$message)) {
        message("\n---------------------------------------------------------")
        message("DIAGNOSIS: Version Mismatch")
        message("R is finding a very old version of a package (pre-2020) in")
        message("your system or user folders that is shadowing the one in")
        message("your Conda environment.")
        message("\nACTION REQUIRED:")
        message("Please run the following command to find and remove it:")
        message("  Rscript -e '.libPaths(); find.package(\"ComplexHeatmap\")'")
        message("\nOnce you see the path, you can remove that specific folder.")
        message("---------------------------------------------------------")
      }

      message("\nCurrent .libPaths() being searched:")
      for (p in .libPaths()) message(paste(" -", p))

      stop("Execution halted due to library conflict.")
    }
  )
}

load_packages()

# --- Input Handling ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  input_file <- args[1]
} else {
  # Priority logic for file finding
  paths_to_check <- c("VF_results/multiple_samples_best_hits.csv", "multiple_samples_best_hits.csv")
  input_file <- paths_to_check[file.exists(paths_to_check)][1]
  if (is.na(input_file)) input_file <- "VF_results/multiple_samples_best_hits.csv"
}
output_file <- ifelse(length(args) >= 2, args[2], "vf_heatmap.pdf")

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file, "\nRun analysis first to generate the matrix."))
}

# --- Data Loading ---
message("Loading and sorting data...")
df <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

# Sort strictly for biological grouping: Broad -> Specific -> Gene (All Alphabetical)
df <- df %>% arrange(VF_Broad_Category, VF_Specific_Name, VF_Gene_Symbol)

# CRITICAL: Convert Specific_Name to a factor preserving the sorted order.
# This ensures that when ComplexHeatmap splits by "Specific_Name", it maintains
# the heirarchy of the Broad_Category and doesn't just sort blocks A-Z.
df$VF_Specific_Name <- factor(df$VF_Specific_Name, levels = unique(df$VF_Specific_Name))
df$VF_Broad_Category <- factor(df$VF_Broad_Category, levels = unique(df$VF_Broad_Category))

# Identify samples
metadata_cols <- c("VF_Broad_Category", "VF_Specific_Name", "VF_Gene_Symbol")
sample_cols <- setdiff(names(df), metadata_cols)
mat <- as.matrix(df[, sample_cols])
rownames(mat) <- make.unique(df$VF_Gene_Symbol)

# --- Define Colors ---
heat_colors <- colorRampPalette(c("#f7f7f7", "#ebf5ff", "#084594"))(100)

cat_unique <- unique(df$VF_Broad_Category)
cat_colors <- setNames(
  colorRampPalette(brewer.pal(min(8, length(cat_unique)), "Set3"))(length(cat_unique)),
  cat_unique
)

spec_unique <- unique(df$VF_Specific_Name)
spec_colors <- setNames(
  colorRampPalette(brewer.pal(min(8, length(spec_unique)), "Pastel1"))(length(spec_unique)),
  spec_unique
)

# --- Layout Configuration ---

# Left-side annotations: Broad (leftmost) then Specific (to the right of broad)
left_ann <- rowAnnotation(
  Broad_Category = df$VF_Broad_Category,
  Specific_Mechanism = df$VF_Specific_Name,
  col = list(Broad_Category = cat_colors, Specific_Mechanism = spec_colors),
  show_annotation_name = TRUE,
  annotation_label = c("Broad Category", "Specific Mechanism")
)

# --- Plotting ---
message("Generating heatmap: ", output_file, "...")

# Define fixed cell size to ensure "square" cells
cell_size <- 5 # mm

# Calculate heatmap body dimensions
ht_width <- ncol(mat) * unit(cell_size, "mm")
ht_height <- nrow(mat) * unit(cell_size, "mm")

# Calculate PDF dimensions (adding margins for annotations, labels and legends)
# Convert mm to inches (1mm = 0.03937 inches)
# Padding: 10 inches for side annotations/labels/legends, 7 inches for top/bottom
pdf_width <- (ncol(mat) * cell_size * 0.03937) + 10
pdf_height <- (nrow(mat) * cell_size * 0.03937) + 7

pdf(output_file, width = pdf_width, height = pdf_height)

# We use 'row_split' to group by Specific_Name and put labels on the right
ht <- Heatmap(
  mat,
  name = "Identity %",
  col = heat_colors,
  width = ht_width,
  height = ht_height,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  column_dend_height = unit(30, "mm"), # Give the dendrogram explicit space
  row_split = df$VF_Specific_Name,
  row_title_side = "right",
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  left_annotation = left_ann,
  column_title = paste("Virulence Gene Identity across Samples\nSource:", basename(input_file)),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_names_side = "right",
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10)
)

# draw the heatmap with legends on the right and extra padding at the top
draw(ht,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(30, 20, 20, 20), "mm") # Top, Right, Bottom, Left padding
)

invisible(dev.off())
message("Success! Heatmap saved to: ", output_file)
