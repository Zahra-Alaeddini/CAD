setwd("C:/Users/ZAHRA/Desktop/Microbiome")

library(randomForest)
library(caret)
library(dplyr)
library(ggplot2)
library(pROC)
library(e1071) 
library(xgboost)


otu_data <- read.csv("main_OTU.csv", row.names = 1) 
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1) 

rownames(metadata) <- tolower(rownames(metadata))
rownames(otu_data) <- tolower(rownames(otu_data))

common_samples <- intersect(rownames(metadata), rownames(otu_data))
otu_data <- otu_data[common_samples, ]
metadata <- metadata[common_samples, ]

otu_data[is.na(otu_data)] <- 0 

print(paste("Total samples:", nrow(metadata)))

otu_data <- otu_data / rowSums(otu_data)

otu_data[is.nan(as.matrix(otu_data))] <- 0

print("Class distribution in y:")
print(table(metadata$healthstatus))
metadata$healthstatus <- factor(metadata$healthstatus, levels = c("control", "CAD"))
X <- otu_data 
y <- metadata$healthstatus 
data_train_full <- cbind(X, healthstatus = y)

set.seed(42) 
train_index <- createDataPartition(y, p = 0.7, list = FALSE)
data_train <- data_train_full[train_index, ]
data_test <- data_train_full[-train_index, ]

ctrl <- trainControl(method = "repeatedcv",
                     number = 5, 
                     repeats = 3, 
                     sampling = "smote", 
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary, 
                     savePredictions = "final",
                     verboseIter = TRUE)

rf_grid <- expand.grid(mtry = seq(2, sqrt(ncol(X)), length.out = 5)) 
rf_model <- train(healthstatus ~ ., data = data_train,
                  method = "rf",
                  trControl = ctrl,
                  tuneGrid = rf_grid,
                  ntree = 500,
                  importance = TRUE)

svm_grid <- expand.grid(C = c(0.1, 1, 10), sigma = c(0.01, 0.1, 1))
svm_model <- train(healthstatus ~ ., data = data_train,
                   method = "svmRadial",
                   trControl = ctrl,
                   tuneGrid = svm_grid,
                   prob.model = TRUE)

xgb_grid <- expand.grid(nrounds = c(50, 100, 200),
                        max_depth = c(3, 6),
                        eta = c(0.01, 0.1),
                        gamma = 0,
                        colsample_bytree = 0.8,
                        min_child_weight = 1,
                        subsample = 0.8)
xgb_model <- train(healthstatus ~ ., data = data_train,
                   method = "xgbTree",
                   trControl = ctrl,
                   tuneGrid = xgb_grid,
                   verbose = FALSE)

rf_prob <- predict(rf_model, data_test, type = "prob")[, "CAD"]
rf_class <- predict(rf_model, data_test)
svm_prob <- predict(svm_model, data_test, type = "prob")[, "CAD"]
svm_class <- predict(svm_model, data_test)
xgb_prob <- predict(xgb_model, data_test, type = "prob")[, "CAD"]
xgb_class <- predict(xgb_model, data_test)

rf_roc <- roc(data_test$healthstatus, rf_prob, levels = c("control", "CAD"))
svm_roc <- roc(data_test$healthstatus, svm_prob, levels = c("control", "CAD"))
xgb_roc <- roc(data_test$healthstatus, xgb_prob, levels = c("control", "CAD"))

get_metrics <- function(class_pred, true_labels, roc_obj) {
  cm <- confusionMatrix(class_pred, true_labels, positive = "CAD")
  data.frame(
    AUC = as.numeric(roc_obj$auc),
    Accuracy = cm$overall["Accuracy"],
    Precision = cm$byClass["Precision"],
    Recall = cm$byClass["Recall"],
    F1_Score = cm$byClass["F1"],
    Kappa = cm$overall["Kappa"]
  )
}
rf_metrics <- get_metrics(rf_class, data_test$healthstatus, rf_roc)
print("RF Confusion Matrix (Test Set):")
print(confusionMatrix(rf_class, data_test$healthstatus, positive = "CAD"))
svm_metrics <- get_metrics(svm_class, data_test$healthstatus, svm_roc)
print("SVM Confusion Matrix (Test Set):")
print(confusionMatrix(svm_class, data_test$healthstatus, positive = "CAD"))
xgb_metrics <- get_metrics(xgb_class, data_test$healthstatus, xgb_roc)
print("XGBoost Confusion Matrix (Test Set):")
print(confusionMatrix(xgb_class, data_test$healthstatus, positive = "CAD"))

rf_prob_full <- predict(rf_model, data_train_full, type = "prob")[, "CAD"]
rf_class_full <- predict(rf_model, data_train_full)
svm_prob_full <- predict(svm_model, data_train_full, type = "prob")[, "CAD"]
svm_class_full <- predict(svm_model, data_train_full)
xgb_prob_full <- predict(xgb_model, data_train_full, type = "prob")[, "CAD"]
xgb_class_full <- predict(xgb_model, data_train_full)

rf_roc_full <- roc(data_train_full$healthstatus, rf_prob_full, levels = c("control", "CAD"))
svm_roc_full <- roc(data_train_full$healthstatus, svm_prob_full, levels = c("control", "CAD"))
xgb_roc_full <- roc(data_train_full$healthstatus, xgb_prob_full, levels = c("control", "CAD"))

rf_metrics_full <- get_metrics(rf_class_full, data_train_full$healthstatus, rf_roc_full)
print("RF Confusion Matrix (Full Dataset):")
print(confusionMatrix(rf_class_full, data_train_full$healthstatus, positive = "CAD"))
svm_metrics_full <- get_metrics(svm_class_full, data_train_full$healthstatus, svm_roc_full)
print("SVM Confusion Matrix (Full Dataset):")
print(confusionMatrix(svm_class_full, data_train_full$healthstatus, positive = "CAD"))
xgb_metrics_full <- get_metrics(xgb_class_full, data_train_full$healthstatus, xgb_roc_full)
print("XGBoost Confusion Matrix (Full Dataset):")
print(confusionMatrix(xgb_class_full, data_train_full$healthstatus, positive = "CAD"))

metrics_table <- rbind(rf_metrics, svm_metrics, xgb_metrics)
metrics_table$Model <- c("RF", "SVM", "XGBoost")
metrics_table <- metrics_table[, c("Model", "AUC", "Accuracy", "Precision", "Recall", "F1_Score", "Kappa")]
print("Table: Performance Metrics (Test Set)")
print(metrics_table)
metrics_table_full <- rbind(rf_metrics_full, svm_metrics_full, xgb_metrics_full)
metrics_table_full$Model <- c("RF", "SVM", "XGBoost")
metrics_table_full <- metrics_table_full[, c("Model", "AUC", "Accuracy", "Precision", "Recall", "F1_Score", "Kappa")]
print("Table: Performance Metrics (Full Dataset)")
print(metrics_table_full)

tiff("ROC_Curve_Comparison.tiff", width = 600, height = 600)
plot(rf_roc, main = "ROC Curve Comparison", col = "green", lwd = 3, legacy.axes = FALSE)
lines(svm_roc, col = "purple", lwd = 2)
lines(xgb_roc, col = "maroon1", lwd = 2)
legend("bottomright", legend = c("Random Forest", "SVM", "XGBoost"),
       col = c("green", "purple", "maroon1"), lwd = 2)
dev.off()

tiff("ROC_Curve_Comparison_Full.tiff", width = 600, height = 600)
plot(rf_roc_full, main = "ROC Curve Comparison (Full Dataset)", col = "green", lwd = 3, legacy.axes = FALSE)
lines(svm_roc_full, col = "purple", lwd = 2)
lines(xgb_roc_full, col = "maroon1", lwd = 2)
legend("bottomright", legend = c("Random Forest", "SVM", "XGBoost"),
       col = c("green", "purple", "maroon1"), lwd = 2)
dev.off()

rf_final <- rf_model$finalModel

gini_importance <- importance(rf_final)

importance_scores <- gini_importance[, "MeanDecreaseGini"]

min_score <- min(importance_scores, na.rm = TRUE)
max_score <- max(importance_scores, na.rm = TRUE)
scaled_importance <- 4 * (importance_scores - min_score) / (max_score - min_score)

important_taxa <- data.frame(
  Feature = names(importance_scores),
  Importance = scaled_importance,
  Genus = sub(".*G__", "", names(importance_scores))
) %>%
  arrange(desc(Importance)) %>%
  head(10) 

print("Top 10 Important Taxa for Diagnosis (Scaled Mean Decrease Gini [0-4]):")
print(important_taxa)

ggplot(important_taxa, aes(x = reorder(Genus, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1.0)) + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, color = "black")
  ) +
  labs(x = "Genus")

ggsave("Feature_Importance_Scaled.tiff", width = 6, height = 6)

