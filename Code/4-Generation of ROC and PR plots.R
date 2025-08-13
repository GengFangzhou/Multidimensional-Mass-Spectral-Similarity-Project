# Load required libraries
library(pROC)      
library(PRROC)     
library(ggplot2)   
library(reshape2)  

# Read data - REPLACE WITH YOUR ACTUAL FILE PATH
data <- read.csv("path/to/your/data.csv")
data$Labels <- as.factor(data$Labels)

# ---- ROC Analysis ----
# Calculate ROC curves
roc_single <- roc(data$Labels, data$S.single)
roc_multi <- roc(data$Labels, data$S.multi)

roc_single_smooth <- smooth(roc_single)
roc_multi_smooth <- smooth(roc_multi)

roc_data <- rbind(
  data.frame(
    FPR = 1 - roc_single_smooth$specificities,
    TPR = roc_single_smooth$sensitivities,
    Method = "S-single"
  ),
  data.frame(
    FPR = 1 - roc_multi_smooth$specificities,
    TPR = roc_multi_smooth$sensitivities,
    Method = "S-multi"
  )
)

write.csv(roc_data, "roc_data.csv", row.names = FALSE)

cat("S-single AUC:", auc(roc_single), "\n")
cat("S-multi AUC:", auc(roc_multi), "\n")

# ---- Precision-Recall Analysis ----
# Calculate PR curves
pr_single <- pr.curve(
  scores.class0 = data$S.single[data$Labels == 1],
  scores.class1 = data$S.single[data$Labels == 0],
  curve = TRUE
)

pr_multi <- pr.curve(
  scores.class0 = data$S.multi[data$Labels == 1],
  scores.class1 = data$S.multi[data$Labels == 0],
  curve = TRUE
)

pr_data <- rbind(
  data.frame(
    Recall = pr_single$curve[, 1],
    Precision = pr_single$curve[, 2],
    Threshold = pr_single$curve[, 3],
    Method = "Single"
  ),
  data.frame(
    Recall = pr_multi$curve[, 1],
    Precision = pr_multi$curve[, 2],
    Threshold = pr_multi$curve[, 3],
    Method = "Multi"
  )
)

write.csv(pr_data, "pr_data.csv", row.names = FALSE)

# ---- Optimal F1 Calculation ----
calculate_f1_metrics <- function(pr_curve) {
  recall <- pr_curve$curve[, 1]
  precision <- pr_curve$curve[, 2]
  threshold <- pr_curve$curve[, 3]
  
  f1_scores <- 2 * precision * recall / (precision + recall)
  best_idx <- which.max(f1_scores)
  
  list(
    f1 = f1_scores[best_idx],
    precision = precision[best_idx],
    recall = recall[best_idx],
    threshold = threshold[best_idx]
  )
}

metrics_single <- calculate_f1_metrics(pr_single)
metrics_multi <- calculate_f1_metrics(pr_multi)

cat("Single-energy model:\n")
cat("  Best F1-score:", round(metrics_single$f1, 4), "\n")
cat("  Precision:", round(metrics_single$precision, 4), "\n")
cat("  Recall:", round(metrics_single$recall, 4), "\n")
cat("  Threshold:", round(metrics_single$threshold, 4), "\n\n")

cat("Multidimensional model:\n")
cat("  Best F1-score:", round(metrics_multi$f1, 4), "\n")
cat("  Precision:", round(metrics_multi$precision, 4), "\n")
cat("  Recall:", round(metrics_multi$recall, 4), "\n")
cat("  Threshold:", round(metrics_multi$threshold, 4), "\n")