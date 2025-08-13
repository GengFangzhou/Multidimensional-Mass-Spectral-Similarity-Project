# Load required libraries
library(MSnbase)
library(dplyr)

# Read MS data file and extract MS2 spectra
file_path <- "path/to/your/file.mzML"  # Generic path placeholder
ms_data <- readMSData(file_path, mode = "onDisk")
ms2_data <- filterMsLevel(ms_data, msLevel = 2)
ms2_metadata <- fData(ms2_data)

# Visualize single spectrum (example)
single_spectrum <- ms2_data[[255]]
plot(single_spectrum, centroided = TRUE, full = TRUE)

# Preprocess metadata
ms2_metadata$cycle <- as.numeric(sub(".*cycle=(\\d+).*", "\\1", ms2_metadata$spectrumId))

# Filter metadata based on RT and cycles
filtered_metadata <- ms2_metadata[ms2_metadata$retentionTime >= a & ms2_metadata$retentionTime <= b, ]
filtered_metadata <- ms2_metadata[ms2_metadata$cycle >= c & ms2_metadata$cycle <= d, ]

# Extract filtered spectra
filtered_spectrum_ids <- filtered_metadata$spectrumId
filtered_metadata$experiment <- as.numeric(
  sub(".*experiment=(\\d+).*", "\\1", filtered_metadata$spectrumId)
)
filtered_ms2_data <- ms2_data[fData(ms2_data)$spectrumId %in% filtered_spectrum_ids, ]

# Process energy scan using defined parameters
ranges <- list(2,3,4,5,6)  # Example energy ranges
final_results <- process_energy_scan(
  filtered_ms2_data, 
  filtered_metadata, 
  ranges
)

# Extract and save results
aligned_avg_intensity_results <- final_results$avg_intensity_results
aligned_rsd_intensity_results <- final_results$rsd_intensity_results

write.csv(aligned_avg_intensity_results, "output_avg_intensity.csv")
write.csv(aligned_rsd_intensity_results, "output_rsd_intensity.csv")

# Core processing functions ---------------------------------------------------
process_cycle <- function(selected_spectra_indices, 
                          filtered_ms2_data, 
                          filtered_metadata, 
                          mz_range_start = 100, 
                          intensity_threshold = 32, 
                          bin_size = 0.025, 
                          delta_mz = 0.035, 
                          tolerance = 0.0001) {
  
  # Extract and trim spectra
  selected_spectra <- filtered_ms2_data[selected_spectra_indices]
  parent_mz <- filtered_metadata$precursorMZ[selected_spectra_indices]
  trimmed_spectrum <- filterMz(selected_spectra, mz = c(mz_range_start, max(parent_mz) + 0.5))
  
  # Noise filtering
  mz_values_list <- lapply(spectra(trimmed_spectrum), function(spectrum) {
    mz_values <- MSnbase::mz(spectrum)
    intensity_values <- MSnbase::intensity(spectrum)
    keep <- intensity_values > intensity_threshold
    mz_values[keep]
  })
  
  intensity_values_list <- lapply(spectra(trimmed_spectrum), function(spectrum) {
    mz_values <- MSnbase::mz(spectrum)
    intensity_values <- MSnbase::intensity(spectrum)
    keep <- intensity_values > intensity_threshold
    intensity_values[keep]
  })
  
  # Merge spectra
  combined_mz <- numeric(0)
  combined_intensity <- numeric(0)
  
  for (i in 1:length(mz_values_list)) {
    current_mz <- mz_values_list[[i]]
    current_intensity <- intensity_values_list[[i]]
    
    for (j in 1:length(current_mz)) {
      mz_tolerance <- abs(combined_mz - current_mz[j]) < tolerance
      if (any(mz_tolerance)) {
        matched_index <- which(mz_tolerance)
        combined_intensity[matched_index] <- combined_intensity[matched_index] + current_intensity[j]
      } else {
        combined_mz <- c(combined_mz, current_mz[j])
        combined_intensity <- c(combined_intensity, current_intensity[j])
      }
    }
  }
  
  # Binning and peak merging
  get_bin_representative_mz <- function(mz_values, intensity_values, binSize) {
    bins <- floor(mz_values / binSize) * binSize
    binned_data <- data.frame(mz = mz_values, intensity = intensity_values, bin = bins)
    
    representative_data <- binned_data %>%
      group_by(bin) %>%
      reframe(
        max_intensity_mz = mz[which.max(intensity)],
        total_intensity = sum(intensity),
        .groups = "drop"
      ) %>%
      filter(total_intensity > 0)
    
    updated_mz_values <- representative_data$max_intensity_mz
    updated_intensity_values <- representative_data$total_intensity
    
    return(list(updated_mz_values = updated_mz_values, updated_intensity_values = updated_intensity_values))
  }
  
  result <- get_bin_representative_mz(combined_mz, combined_intensity, bin_size)
  updated_mz_values <- result$updated_mz_values
  updated_intensity_values <- result$updated_intensity_values
  
  # Peak merging
  merge_adjacent_peaks <- function(mz_values, intensity_values, delta_mz) {
    merged_mz <- c()
    merged_intensity <- c()
    i <- 1
    n <- length(mz_values)
    
    while (i <= n) {
      current_bin_mz <- mz_values[i]
      current_bin_intensity <- intensity_values[i]
      j <- i + 1
      
      while (j <= n && mz_values[j] - current_bin_mz <= delta_mz) {
        if (intensity_values[j] > current_bin_intensity) {
          current_bin_mz <- mz_values[j]
        }
        current_bin_intensity <- current_bin_intensity + intensity_values[j]
        j <- j + 1
      }
      
      merged_mz <- c(merged_mz, current_bin_mz)
      merged_intensity <- c(merged_intensity, current_bin_intensity)
      i <- j
    }
    
    return(data.frame(mz = merged_mz, intensity = merged_intensity))
  }
  
  merged_data <- merge_adjacent_peaks(updated_mz_values, updated_intensity_values, delta_mz)
  
  
  return(merged_data)
  
}

align_peaks <- function(results, tolerance = 0.025) {
  # Initialize alignment with first result
  aligned_peaks <- data.frame(mz = as.numeric(results[[1]]$mz))
  intensity_col <- names(results[[1]])[2]
  aligned_peaks[[paste0("Intensity_1")]] <- as.numeric(results[[1]][[intensity_col]])
  
  # Iterate through remaining results
  for (i in 2:length(results)) {
    single_result <- results[[i]]
    single_result$mz <- as.numeric(single_result$mz)
    
    intensity_col <- names(single_result)[2]
    single_result[[intensity_col]] <- as.numeric(single_result[[intensity_col]])
    
    for (j in seq_along(single_result$mz)) {
      mz_val <- single_result$mz[j]
      intensity_val <- single_result[[intensity_col]][j]
    
      # Find matching mz
      match_idx <- which(abs(aligned_peaks$mz - mz_val) <= tolerance)
      if (length(match_idx) > 0) {
        col_name <- paste0("Intensity_", i)
        if (!col_name %in% colnames(aligned_peaks)) {
          aligned_peaks[[col_name]] <- 0
        }
        aligned_peaks[match_idx, col_name] <- aligned_peaks[match_idx, col_name] + intensity_val
      } else {
        new_row <- setNames(as.list(rep(0, ncol(aligned_peaks))), colnames(aligned_peaks))
        new_row$mz <- mz_val
        new_row[[paste0("Intensity_", i)]] <- intensity_val
        aligned_peaks <- rbind(aligned_peaks, new_row)
      }
    }
  }
  
  # Clean and sort
  aligned_peaks[is.na(aligned_peaks)] <- 0
  aligned_peaks <- aligned_peaks[order(aligned_peaks$mz), ]
  
  return(aligned_peaks)
}

process_peaks <- function(filtered_peaks, threshold = 0.001) {
  # Filter and normalize
  range_peaks <- filter(filtered_peaks, mz >= a)
  for (col in colnames(filtered_peaks)[-1]) {  # 跳过 "mz" 列
    max_val <- max(range_peaks[[col]], na.rm = TRUE)
    if (max_val > 0) {
      filtered_peaks[[col]] <- filtered_peaks[[col]] / max_val
    }
  }
  
  # Filter valid peaks
  valid_counts <- apply(normalized_peaks[, -1], 1, function(x) sum(x > 0, na.rm = TRUE))
  filtered_peaks <- normalized_peaks[valid_counts >= 8, ]
  
  # Calculate statistics
  avg_intensity <- apply(filtered_peaks[, -1], 1, mean, na.rm = TRUE)
  rsd_intensity <- apply(filtered_peaks[, -1], 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100)
  result <- data.frame(mz = filtered_peaks$mz, 
                       avg_intensity = avg_intensity, 
                       rsd_intensity = rsd_intensity)
  result <- result[avg_intensity >= threshold, ]
  
  return(result)
}

process_energy_scan <- function(filtered_ms2_data, filtered_metadata, ranges, tolerance = 0.025, threshold = 0.001) {
  aligned_avg_intensity_results <- list()
  aligned_rsd_intensity_results <- list()
  
  avg_intensity_lists <- list()
  rsd_intensity_lists <- list()
  
  for (range_idx in seq_along(ranges)) {
    current_range <- ranges[[range_idx]]
    range_indices <- which(filtered_metadata$experiment %in% current_range)
    
    if (length(range_indices) > 0) {
      unique_cycles <- unique(filtered_metadata$cycle[range_indices])
      cycle_results <- list()
      
      for (cycle_id in unique_cycles) {
        cycle_indices <- which(
          filtered_metadata$cycle == cycle_id & 
            filtered_metadata$experiment %in% current_range
        )
        
        if (length(cycle_indices) > 0) {
          processed_data <- process_cycle(
            selected_spectra_indices = cycle_indices,
            filtered_ms2_data = filtered_ms2_data,
            filtered_metadata = filtered_metadata,
            mz_range_start = 100,
            intensity_threshold = 32,
            bin_size = 0.025,
            delta_mz = 0.03,
            tolerance = 0.0001
          )
          cycle_results[[as.character(cycle_id)]] <- processed_data
        }
      }
      
      aligned_data <- align_peaks(cycle_results, tolerance = tolerance)
      cat("Aligned Peaks for Range Index", range_idx, ":\n")
      print(tail(aligned_data))
      
      processed_peaks <- process_peaks(aligned_data, threshold = threshold)
      cat("Processed Peaks for Range Index", range_idx, ":\n")
      print(processed_peaks)
      
      avg_list[[range_idx]] <- data.frame(mz = processed$mz, avg_intensity = processed$avg_intensity)
      rsd_list[[range_idx]] <- data.frame(mz = processed$mz, rsd_intensity = processed$rsd_intensity)
    }
  }
  
  # Final alignment and labeling
  if (length(avg_intensity_lists) > 0) {
    aligned_avg_intensity_results <- align_peaks(avg_intensity_lists, tolerance = tolerance)
  }
  if (length(rsd_intensity_lists) > 0) {
    aligned_rsd_intensity_results <- align_peaks(rsd_intensity_lists, tolerance = tolerance)
  }
  
  colnames(aligned_avg_intensity_results) <- c("mz", "8EV", "10EV", "15EV", "20EV", "25EV")
  colnames(aligned_rsd_intensity_results) <- c("mz", "8EV", "10EV", "15EV", "20EV", "25EV")
  
  return(list(avg_intensity_results = aligned_avg_intensity_results, rsd_intensity_results = aligned_rsd_intensity_results))
}