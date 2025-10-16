# =============================================================================
# DATA PREPARATION SCRIPT (SERVER & FULL DATASET VERSION)
#
# Purpose: This script performs all heavy, one-time computations for the app.
# It loads the full Seurat object, pre-screens all viable genes with geneCOCOA,
# and saves the final plot-ready results to a compact file.
# This script is designed to be run once on the server.
# =============================================================================

# --- 1. Load Libraries ---
print("--- Loading Libraries ---")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(msigdbr)
  library(geneCOCOA)
  library(pbapply)
  library(parallel)
})

# --- 2. Configuration ---
print("--- Setting Parameters ---")
# Input/Output Files
CACHED_DATA_PATH <- "data/AMI.data.sets_V5_updated.RDS"
OUTPUT_RESULTS_FILE <- "precomputed_results.rds"

# Analysis Parameters
MIN_CELLS_PER_GENE <- 20 # Minimum number of cells a gene must be expressed in to be analyzed
N_SIMULATIONS <- 250      # Number of simulations for geneCOCOA null distribution

# Parallel Processing
NUM_CORES <- detectCores() - 1

# --- 3. Check for Existing Results ---
if (file.exists(OUTPUT_RESULTS_FILE)) {
  print("Precomputed results file already exists. Skipping.")
  quit(save = "no")
}

# --- 4. Load Full Dataset ---
print(paste("--- Loading full dataset from:", CACHED_DATA_PATH, "---"))
if (!file.exists(CACHED_DATA_PATH)) {
  stop("FATAL: Cached data file not found at the specified path!")
}
full_ami_data <- readRDS(CACHED_DATA_PATH)

# --- 5. Prepare Data Splits and Gene Sets ---
print("--- Preparing Data Splits and Gene Sets ---")
control_expr_matrix <- as.data.frame(as.matrix(GetAssayData(full_ami_data, assay = "RNA", layer = "data")[, full_ami_data$Condition == "Control"]))
ami_expr_matrix <- as.data.frame(as.matrix(GetAssayData(full_ami_data, assay = "RNA", layer = "data")[, full_ami_data$Condition == "AMI"]))

# Fetch Hallmark gene sets
print("--- Fetching MSigDB Hallmark Gene Sets ---")
m_df <- msigdbr(species = "Mus musculus", category = "H")
hallmark_sets_raw <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# Filter gene sets to only include genes present in each expression matrix
filter_hallmark_sets <- function(expr_matrix, hallmark_sets) {
  available_genes <- rownames(expr_matrix)
  filtered_sets <- lapply(hallmark_sets, function(genes) intersect(genes, available_genes))
  filtered_sets[sapply(filtered_sets, length) > 1] # Ensure gene sets are not empty
}

ami_hallmark_sets <- filter_hallmark_sets(ami_expr_matrix, hallmark_sets_raw)
control_hallmark_sets <- filter_hallmark_sets(control_expr_matrix, hallmark_sets_raw)

# --- 6. Identify Viable Genes for Analysis ---
print("--- Identifying viable genes for pre-analysis ---")
all_genes <- intersect(rownames(ami_expr_matrix), rownames(control_expr_matrix))
viable_genes_mask <- pbsapply(all_genes, function(g) {
  (sum(ami_expr_matrix[g, ] > 0) > MIN_CELLS_PER_GENE) &&
  (sum(control_expr_matrix[g, ] > 0) > MIN_CELLS_PER_GENE) &&
  (sd(ami_expr_matrix[g, ]) > 0) &&
  (sd(control_expr_matrix[g, ]) > 0)
})
genes_to_run <- names(viable_genes_mask[viable_genes_mask])
print(paste("Found", length(genes_to_run), "viable genes to analyze out of", length(all_genes), "total genes."))

# --- 7. Define Core Analysis Function ---
run_cocoa_analysis <- function(goi, expr_matrix, hallmark_sets) {
  tryCatch({
    expr_info <- get_expr_info(expr = expr_matrix, GOI = goi)
    get_stats(
      geneset_collection = hallmark_sets,
      GOI = goi,
      GOI_expr = expr_info$GOI_expr,
      expr_df = expr_info$expr_df,
      nsims = N_SIMULATIONS
    )
  }, error = function(e) {
    # message("Error processing gene '", goi, "': ", e$message)
    NULL # Return NULL on error
  })
}

# --- 8. Run Pre-analysis for all Viable Genes ---
print(paste("--- Starting pre-analysis for", length(genes_to_run), "viable genes. This will take a long time. ---"))
print(paste("--- Using", NUM_CORES, "cores for parallel processing. ---"))

# Setup for parallel processing
cl <- makeCluster(max(1, NUM_CORES))
clusterEvalQ(cl, { library(geneCOCOA) })
clusterExport(cl, c("get_expr_info", "get_stats", "N_SIMULATIONS", "run_cocoa_analysis"))

# Main analysis loop
all_results <- pblapply(genes_to_run, function(goi) {
  # Run for AMI condition
  res_ami <- run_cocoa_analysis(goi, ami_expr_matrix, ami_hallmark_sets)
  
  # Run for Control condition
  res_control <- run_cocoa_analysis(goi, control_expr_matrix, control_hallmark_sets)
  
  # Only return if both analyses were successful
  if (!is.null(res_ami) && !is.null(res_control)) {
    return(list(ami = res_ami$stats_summary, control = res_control$stats_summary))
  } else {
    return(NULL)
  }
}, cl = cl)

stopCluster(cl)

# --- 9. Format and Save Results ---
print("--- Formatting and saving results ---")
names(all_results) <- genes_to_run
successful_results <- all_results[!sapply(all_results, is.null)]

saveRDS(successful_results, file = OUTPUT_RESULTS_FILE)

print("-------------------------------------------------")
print("Data preparation and pre-analysis complete!")
print(paste("Successfully analyzed", length(successful_results), "out of", length(genes_to_run), "genes."))
print(paste("Results saved to:", OUTPUT_RESULTS_FILE))
print("You can now launch the Shiny app.")
print("-------------------------------------------------")