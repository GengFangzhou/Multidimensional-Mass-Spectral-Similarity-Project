# Load required packages
library(dplyr)

# Function to merge and align peak files
merge_and_align_peaks <- function(file_list, mz_tolerance = 0.02) {
  if (!all(file.exists(file_list))) {
    stop("Missing files:", paste(file_list[!file.exists(file_list)], collapse = ", "))
  }
  
  # Read and preprocess files
  data_list <- lapply(file_list, function(file) {
    df <- read.csv(file)
    colnames(df) <- tolower(gsub("\\s+", "_", colnames(df)))
    if (!"mz" %in% colnames(df)) {
      stop(paste("Missing 'mz' column in file：", file))
    }
    df <- df %>% mutate(across(everything(), as.numeric))
    return(df)
  })
  
  # Extract all unique m/z values
  all_mz <- unique(do.call(c, lapply(data_list, function(df) df$mz)))
  aligned_peaks <- data.frame(mz = sort(all_mz))
  
  # Align peaks across files
  for (i in seq_along(data_list)) {
    single_result <- data_list[[i]]
    energy_columns <- setdiff(colnames(single_result), "mz")
    
    for (energy in energy_columns) {
      col_name <- paste0(energy, "_", i)
      aligned_peaks[[col_name]] <- NA 
      
      for (j in seq_along(single_result$mz)) {
        mz_val <- single_result$mz[j]
        intensity_val <- single_result[[energy]][j]
        
        match_idx <- which(abs(aligned_peaks$mz - mz_val) <= mz_tolerance)
        if (length(match_idx) > 0) {
          aligned_peaks[match_idx, col_name] <- intensity_val
        }
      }
    }
  }
  
  aligned_peaks[is.na(aligned_peaks)] <- 0
  
  # Aggregate nearby peaks
  aggregated_peaks <- aligned_peaks %>%
    mutate(group_id = cumsum(c(1, diff(mz) > mz_tolerance))) %>% 
    group_by(group_id) %>%
    summarise(
      mz = median(mz),# mz 取中位数
      across(starts_with("x8ev"), \(x) mean(x, na.rm = TRUE)),
      across(starts_with("x10ev"), \(x) mean(x, na.rm = TRUE)),
      across(starts_with("x15ev"), \(x) mean(x, na.rm = TRUE)),
      across(starts_with("x20ev"), \(x) mean(x, na.rm = TRUE)),
      across(starts_with("x25ev"), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    select(-group_id)  
  
  # Calculate average intensities
  avg_peaks <- aggregated_peaks %>%
    mutate(
      `8EV` = rowMeans(select(., starts_with("x8ev")), na.rm = TRUE),
      `10EV` = rowMeans(select(., starts_with("x10ev")), na.rm = TRUE),
      `15EV` = rowMeans(select(., starts_with("x15ev")), na.rm = TRUE),
      `20EV` = rowMeans(select(., starts_with("x20ev")), na.rm = TRUE),
      `25EV` = rowMeans(select(., starts_with("x25ev")), na.rm = TRUE)
    ) %>%
    select(mz, `8EV`, `10EV`, `15EV`, `20EV`, `25EV`) 
  
  # Filter weak peaks
  avg_peaks <- avg_peaks %>%
    filter(
      `8EV` >= 0.001 | `10EV` >= 0.001 | `15EV` >= 0.001 | `20EV` >= 0.001 | `25EV` >= 0.001
    )
  
  return(avg_peaks)
}

# Ecample usage
{
  file_list <- c("path/to/file1.csv",
                 "path/to/file2.csv",
                 "path/to/file3.csv")
  
  result <- merge_and_align_peaks(file_list, mz_tolerance = 0.02)
  write.csv(result, "merge_peaks.csv", row.names = FALSE)
}
