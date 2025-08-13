# Load required libraries
library(SnowballC)
library(lsa)

# Calculate weighted cosine similarity
weighted_cosine <- function(mass_q, intensity_q, mass_p, intensity_p, a, b) {
  Q_prime <- (mass_q ^ a) * (intensity_q ^ b)
  P_prime <- (mass_p ^ a) * (intensity_p ^ b)
  
  dot_product <- sum(Q_prime * P_prime)
  
  norm_Q <- sqrt(sum(Q_prime ^ 2))
  norm_P <- sqrt(sum(P_prime ^ 2))
  
  if (norm_Q > 0 & norm_P > 0) {
    similarity <- dot_product / (norm_Q * norm_P)
  } else {
    similarity <- 0
  }
  
  return(similarity)
}

# Function to merge two spectra based on m/z tolerance
merge_spectra <- function(mz_A, intensities_A, mz_B, intensities_B, tolerance = 0.02) {
  
  combined_mz <- mz_A  
  combined_intensities_A <- intensities_A 
  combined_intensities_B <- rep(0, length(mz_A))  
  
  for (i in seq_along(mz_B)) {
    close_to_A <- abs(mz_A - mz_B[i]) <= tolerance
    if (any(close_to_A, na.rm = TRUE))  {
      combined_intensities_B[which(close_to_A)] <- intensities_B[i]
    } else {
      combined_mz <- c(combined_mz, mz_B[i])
      combined_intensities_A <- c(combined_intensities_A, 0)  
      combined_intensities_B <- c(combined_intensities_B, intensities_B[i]) 
    }
  }
  
  order_index <- order(combined_mz)
  combined_mz <- combined_mz[order_index]
  combined_intensities_A <- combined_intensities_A[order_index]
  combined_intensities_B <- combined_intensities_B[order_index]
  
  return(list(mz = combined_mz, combined_intensities_A = combined_intensities_A, combined_intensities_B = combined_intensities_B))
}

# Calculate cosine similarity between two spectra
calculate_cosine_similarity <- function(mz_A, intensities_A, mz_B, intensities_B, tolerance = 0.01) {
  # Merge spectra
  merged <- merge_spectra(mz_A, intensities_A, mz_B, intensities_B, tolerance)
  
  mz_A_merged <- merged$mz
  intensities_A_merged <- merged$combined_intensities_A
  mz_B_merged <- merged$mz
  intensities_B_merged <- merged$combined_intensities_B
  
  # Apply intensity weighting
  Q_prime <- intensities_A_merged ^ 0.5
  P_prime <- intensities_B_merged ^ 0.5
  
  dot_product <- sum(Q_prime * P_prime)
  norm_Q <- sqrt(sum(Q_prime ^ 2))
  norm_P <- sqrt(sum(P_prime ^ 2))
  similarity <- ifelse(norm_Q > 0 & norm_P > 0, dot_product / (norm_Q * norm_P), 0)
  
  return(similarity)
}

# Calculate similarity between two intensity matrices
calculate_matrix_similarity <- function(matrix_A, matrix_B, tolerance = 0.03) {
  mz_A <- matrix_A[, 1]
  mz_B <- matrix_B[, 1]
  
  intensity_cols <- setdiff(colnames(matrix_A), colnames(matrix_A)[1])
  similarities <- c()
  
  for (col_name in intensity_cols) {
    intensities_A <- matrix_A[[col_name]]
    intensities_B <- matrix_B[[col_name]]
    
    similarity <- calculate_cosine_similarity(
      mz_A, intensities_A,
      mz_B, intensities_B,
      tolerance
    )
    
    similarities <- c(similarities, similarity)
  }
  print(similarities)
  
  # Compute geometric mean of similarities
  similarity_product <- prod(similarities)
  final_similarity <- similarity_product^(1 / length(similarities))
  
  return(final_similarity)
}

# Modified matrix similarity calculation for comparison
calculate_matrix_similarity <- function(matrix_A, matrix_B, tolerance = 0.02, target_cols = NULL) {
  mz_A <- as.numeric(matrix_A[[1]])
  mz_B <- as.numeric(matrix_B[[1]])
  
  all_intensity_cols <- colnames(matrix_A)[-1]
  
  if (is.null(target_cols)) {
    target_cols <- all_intensity_cols
  } else {
    target_cols <- intersect(target_cols, all_intensity_cols)
    if (length(target_cols) == 0) {
      stop
    }
  }
  
  similarities <- c()
  for (col_name in target_cols) {
    intensities_A <- as.numeric(matrix_A[[col_name]])
    intensities_B <- as.numeric(matrix_B[[col_name]])
    
    similarity <- calculate_cosine_similarity(
      mz_A, intensities_A,
      mz_B, intensities_B,
      tolerance
    )
    similarities <- c(similarities, similarity)
  }
  print(similarities)
  
  similarity_product <- prod(similarities)
  final_similarity_all <- similarity_product^(1 / length(similarities))
  if ("15EV" %in% all_intensity_cols) {
    intensities_A_15EV <- as.numeric(matrix_A[["15EV"]])
    intensities_B_15EV <- as.numeric(matrix_B[["15EV"]])
    
    similarity_15EV <- calculate_cosine_similarity(
      mz_A, intensities_A_15EV,
      mz_B, intensities_B_15EV,
      tolerance
    )
  } else {
    stop
  }
  
  return(list(
    similarity_15EV = similarity_15EV,
    final_similarity_all = final_similarity_all
  ))
}

# Calculate pairwise similarity matrix for multiple matrix
calculate_muti_similarity <- function(matrices) {
  num_matrices <- length(matrices)
  result_matrix <- matrix(0, nrow = num_matrices, ncol = num_matrices)
  
  for (i in 1:num_matrices) {
    for (j in 1:num_matrices) {
      if (i != j) {
        result_matrix[i, j] <- calculate_matrix_similarity(matrices[[i]],matrices[[j]],tolerance = 0.025)
      }else{
        result_matrix[i, j] <- 1
      }
    }
  }
  
  return(result_matrix)
}

# Call function and save results
matrix_X <- read.csv("path/to/your/ER-MS matrix.csv")
matrices <- list(matrix_A, ……,matrix_X)
similarity_results <- sapply(matrices, function(mat) calculate_matrix_similarity(matrix_X, mat, tolerance = 0.025))
print(similarity_results)

# Bulk similarity result calculation
process_similarity <- function(file_name) {
  file_path_A <- file.path("path/to/your/folder-1", file_name)
  file_path_B <- file.path("path/to/your/folder-2", file_name)
  
  if (!file.exists(file_path_A) || !file.exists(file_path_B)) {
    message("Skipping ", file_name, " (missing in one directory)")
    return(NULL)
  }
  
  matrix_A <- read.csv(file_path_A)
  matrix_B <- read.csv(file_path_B)
  similarity <- calculate_matrix_similarity(matrix_A, matrix_B, tolerance = 0.025)
  
  return(data.frame(File = file_name, Similarity = similarity))
}

csv_files <- list.files("path/to/your/folder-1", pattern = "\\.csv$", full.names = FALSE)
similarity_results <- map_dfr(csv_files, process_similarity)

output_file <- "path/to/your/similarity result.csv"
write.csv(similarity_results, file = output_file, row.names = FALSE)
message("All similarity calculations completed! Results saved to: ", output_file)