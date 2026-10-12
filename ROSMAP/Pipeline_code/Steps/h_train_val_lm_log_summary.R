#Load libraries
library(readr)
library(pheatmap)  # Only keep if pheatmap is still used (is it??? checl later)


# Capture arguments from the main wrapper script
args <- commandArgs(trailingOnly = TRUE)
main_directory <- args[1]
cell_type <- args[2] # Major cell type being processed i.e. Astrocytes

# Working directory from wrapper
setwd(main_directory)

# Construct patterns for clustered genes per subtype file imports
train_pattern <- paste0(cell_type, "_.*_train_gene_clust_metacell_zscorefilt.csv")
validation_pattern <- paste0(cell_type, "_.*_validation_gene_clust_metacell_zscorefilt.csv")

# Define paths where to look for those gene cluster files
train_folder <- file.path(main_directory, "Gene_clusters")
validation_folder <- file.path(main_directory, "Gene_clusters")

# List all files matching the pattern for training and validation sets
training_file_list <- list.files(path=train_folder, pattern=train_pattern, full.names=TRUE)
validation_file_list <- list.files(path=validation_folder, pattern=validation_pattern, full.names=TRUE)

# CHECK: Print the files found for training and validation
print("Training files:")
print(training_file_list)
print("Validation files:")
print(validation_file_list)

# Define phenotypes to be used in lm
pheno <- c("E2", "E4", "new_cogdx", "age_death", "braaksc", "ceradsc", "gpath", 
           "niareagansc", "amyloid", "tangles", "arteriol_scler", "caa_4gp", "cvda_4gp2", "ci_num2_mct")

# Phenotypes to log1p-transform BEFORE modeling
log_pheno <- c("tangles", "amyloid", "gpath")

safe_log1p <- function(x, colname = "") {
  x <- suppressWarnings(as.numeric(x))
  if (any(!is.na(x) & x < 0)) {
    warning(sprintf("Column '%s' has negative values; skipping log1p transform.", colname))
    return(x)
  }
  return(log1p(x))
}


# Ensure the output directory for LM results exists, create if not
output_dir <- file.path(main_directory, "LMER_res_nosex", cell_type)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create empty summary dataframe once at the very start
# to store summaries of each subtype dataset
summary_df <- data.frame(
  file_name = character(),
  n_metacells = numeric(),
  n_donors = numeric(),
  n_clusters = numeric(),
  braak_categories = character(),
  cerad_categories = character(),
  new_cogdx_categories = character(),
  mean_sparsity = numeric(),
  sd_sparsity = numeric(),
  mean_age_death = numeric(),
  sd_age_death = numeric(),
  sex_ratio = numeric(),
  missing_data_cols = character(), # Pheno that has ANY missing data
  data_set = character(),          # column to identify training vs validation
  stringsAsFactors = FALSE
)

# Define a function to log messages with a structured format
log_entry <- function(type, base_name, cluster, phenotype, message) {
  file_path <- file.path(output_dir, paste0(cell_type, "_model_fitting_log_nosex.csv"))
  full_message <- sprintf('%s, %s, %s, %s, %s, "%s"', 
                          Sys.time(), type, base_name, cluster, phenotype, message)
  # Check if the log file exists; if not, write the header
  if (!file.exists(file_path)) {
    writeLines("Date and Time, Type, Model Base Name, Cluster, Phenotype, Message", file_path)
  }
  # Append the full log message
  writeLines(full_message, file_path, append = TRUE)
}

# Function to run simple linear models on the given file list
run_lm_analysis <- function(file_list, data_set_label) {
  for (file_path in file_list) {
    metacells <- read_csv(file_path, show_col_types = FALSE)
    # Apply log1p transform to selected phenotypes (train + validation)
    for (col in log_pheno) {
      if (col %in% colnames(metacells)) {
        metacells[[col]] <- safe_log1p(metacells[[col]], colname = col)
      } else {
        log_entry("Warning", base_name, NA, col, "Log-transform column not found in dataset; skipping.")
      }
    }
    
    clusters <- grep("Cluster_", colnames(metacells), value = TRUE)
    base_name <- gsub("_gene_clust_metacell_zscorefilt.csv", "", basename(file_path))
    
    cat("Processing", base_name, "\n")
    
    # Create summary for this file (subtype) and add to the summary_df
    summary_row <- data.frame(
      file_name = base_name,
      n_metacells = nrow(metacells),
      n_donors = length(unique(metacells$projid)),
      n_clusters = length(clusters),
      braak_categories = paste(sort(unique(na.omit(metacells$braaksc))), collapse = "|"),
      cerad_categories = paste(sort(unique(na.omit(metacells$ceradsc))), collapse = "|"),
      new_cogdx_categories = paste(sort(unique(na.omit(metacells$new_cogdx))), collapse = "|"),
      mean_sparsity = mean(metacells$sparsity, na.rm = TRUE),
      sd_sparsity = sd(metacells$sparsity, na.rm = TRUE),
      mean_age_death = mean(metacells$age_death, na.rm = TRUE),
      sd_age_death = sd(metacells$age_death, na.rm = TRUE),
      sex_ratio = sum(metacells$msex == 1, na.rm = TRUE) / nrow(metacells),
      missing_data_cols = paste(names(which(colSums(is.na(metacells[, pheno])) > 0)), collapse = "|"),
      data_set = data_set_label,
      stringsAsFactors = FALSE
    )
    
    summary_df <<- rbind(summary_df, summary_row)
    
    # Initialize matrix to store LM p-values and coefficients
    p_value_matrix <- matrix(NA, nrow = length(clusters), ncol = length(pheno), 
                             dimnames = list(clusters, pheno))
    coef_matrix <- matrix(NA, nrow = length(clusters), ncol = length(pheno), 
                          dimnames = list(clusters, pheno))
    
    # For each cluster and each phenotype:
    for (cluster in clusters) {
      for (phenotype in pheno) {
        # Create an environment to store warnings that may be generated during model fitting
        warning_collector <- new.env()
        warning_collector$warnings <- character()
        
        # Fit LM using tryCatch to handle errors and warnings
        model <- tryCatch({
          withCallingHandlers(
            {
              # Define formula dynamically: cluster ~ phenotype + msex + sparsity
              formula <- as.formula(paste(cluster, "~", phenotype, "+ msex + sparsity"))
              fitted_model <- lm(formula, data = metacells)
              
              # Return fitted model
              fitted_model
            },
            warning = function(w) {
              # Collect warnings without interrupting execution
              warning_collector$warnings <- c(warning_collector$warnings, w$message)
              # Suppress the printing of warnings to console
              invokeRestart("muffleWarning")
            }
          )
        },
        error = function(e) {
          # Log errors encountered during model fitting
          log_entry("Error", base_name, cluster, phenotype, e$message)
          return(NULL)
        })
        
        # Process model if successful
        if (!is.null(model)) {
          # Log any collected warnings
          if (length(warning_collector$warnings) > 0) {
            for (warn in warning_collector$warnings) {
              log_entry("Warning", base_name, cluster, phenotype, warn)
            }
          }
          # Extract and store results
          coef_summary <- summary(model)$coefficients
          if (phenotype %in% rownames(coef_summary)) {
            p_value_matrix[cluster, phenotype] <- coef_summary[phenotype, "Pr(>|t|)"]
            coef_matrix[cluster, phenotype] <- coef_summary[phenotype, "Estimate"]
          }
        }
      }
    }
    
    # Save results
    write.csv(p_value_matrix, file = paste0(output_dir, "/", base_name, "_lmer_pvals_nosex.csv"))
    write.csv(coef_matrix, file = paste0(output_dir, "/", base_name, "_lmer_coef_nosex.csv"))
  }
}

# Run LM analysis for the training set and validation set
run_lm_analysis(training_file_list, "train")
run_lm_analysis(validation_file_list, "validation")

write.csv(summary_df, 
          file = paste0(output_dir, "/", cell_type, "_combined_data_summary_nosex.csv"), 
          row.names = FALSE)
