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
  library(pbapply) # For the progress bar
})

# --- 2. Define File Paths ---
print("--- Setting Parameters ---")
cached_data_path <- "data/AMI.data.sets_V5_updated.RDS"
output_results_file <- "precomputed_results.rds"

# --- 3. Check for Existing Results ---
if (file.exists(output_results_file)) {
  print("Precomputed results file already exists. Skipping.")
  quit(save = "no")
}

# --- 4. Load Full Dataset ---
print(paste("--- Loading full dataset from:", cached_data_path, "---"))
if (!file.exists(cached_data_path)) stop("Cached data file not found!")
full_ami_data <- readRDS(cached_data_path)

# --- 5. Prepare Data Splits and Gene Sets ---
print("--- Preparing Data Splits and Gene Sets ---")
control_expr_matrix <- as.data.frame(as.matrix(GetAssayData(full_ami_data, assay = "RNA", layer = "data")[, full_ami_data$Condition == "Control"]))
ami_expr_matrix <- as.data.frame(as.matrix(GetAssayData(full_ami_data, assay = "RNA", layer = "data")[, full_ami_data$Condition == "AMI"]))

m_df <- msigdbr(species = "Mus musculus", category = "H")
hallmark_sets_raw <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

ami_available_genes <- rownames(ami_expr_matrix)
ami_hallmark_sets <- lapply(hallmark_sets_raw, function(genes) intersect(genes, ami_available_genes))
ami_hallmark_sets <- ami_hallmark_sets[sapply(ami_hallmark_sets, length) > 1]

control_available_genes <- rownames(control_expr_matrix)
control_hallmark_sets <- lapply(hallmark_sets_raw, function(genes) intersect(genes, control_available_genes))
control_hallmark_sets <- control_hallmark_sets[sapply(control_hallmark_sets, length) > 1]

# --- 6. Identify Viable Genes for Analysis ---
print("--- Identifying viable genes for pre-analysis ---")
all_genes <- intersect(rownames(ami_expr_matrix), rownames(control_expr_matrix))
viable_genes <- pbsapply(all_genes, function(g) {
  sum(ami_expr_matrix[g, ] > 0) > 20 && sum(control_expr_matrix[g, ] > 0) > 20 &&
  sd(ami_expr_matrix[g, ]) > 0 && sd(control_expr_matrix[g, ]) > 0
})
genes_to_run <- names(viable_genes[viable_genes])
print(paste("Found", length(genes_to_run), "viable genes to analyze."))

# --- 7. Run Pre-analysis for all Viable Genes ---
print("--- Starting pre-analysis for all viable genes. This will take a long time. ---")

# Setup for parallel processing if possible
num_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(max(1, num_cores))
parallel::clusterEvalQ(cl, { library(geneCOCOA) })

all_results <- pblapply(genes_to_run, function(goi) {
  
  res_ami <- tryCatch({
    expr_info <- get_expr_info(expr = ami_expr_matrix, GOI = goi)
    get_stats(geneset_collection = ami_hallmark_sets, GOI = goi, 
              GOI_expr = expr_info$GOI_expr, expr_df = expr_info$expr_df, nsims = 250)
  }, error = function(e) { NULL })

  res_control <- tryCatch({
    expr_info <- get_expr_info(expr = control_expr_matrix, GOI = goi)
    get_stats(geneset_collection = control_hallmark_sets, GOI = goi, 
              GOI_expr = expr_info$GOI_expr, expr_df = expr_info$expr_df, nsims = 250)
  }, error = function(e) { NULL })
  
  # Only return if both analyses were successful
  if (!is.null(res_ami) && !is.null(res_control)) {
    return(list(ami = res_ami$stats_summary, control = res_control$stats_summary))
  } else {
    return(NULL)
  }
}, cl = cl)

parallel::stopCluster(cl)

# --- 8. Format and Save Results ---
print("--- Formatting and saving results ---")
names(all_results) <- genes_to_run
# Filter out any genes that failed during the run
successful_results <- all_results[!sapply(all_results, is.null)]

saveRDS(successful_results, file = output_results_file)

print("-------------------------------------------------")
print("Data preparation and pre-analysis complete!")
print(paste("Successfully analyzed", length(successful_results), "genes."))
print(paste("Results saved to:", output_results_file))
print("You can now launch the Shiny app.")
print("-------------------------------------------------")

